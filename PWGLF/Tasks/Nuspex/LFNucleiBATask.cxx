// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

///
/// \file LFNucleiBATask.cxx
///
/// \brief  Analysis task for the measurement of the coalescence parameter B2/B3 in pp collisions for (anti)deuteron/(anti)helium-3
///
/// \author Giovanni Malfattore <giovanni.malfattore@cern.ch> and Rutuparna Rath <rutuparna.rath@cern.ch>
///

#include "PWGLF/DataModel/LFNucleiTables.h"
#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "PWGLF/DataModel/mcCentrality.h"
#include "PWGLF/Utils/inelGt.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "EventFiltering/Zorro.h"
#include "EventFiltering/ZorroSummary.h"

#include "CCDB/BasicCCDBManager.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/Track.h"

#include "TMCProcess.h"
#include <TF1.h>

#include <gsl/span>
#include <memory>
#include <string>
#include <unordered_set>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct LFNucleiBATask {
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Service<o2::framework::O2DatabasePDG> pdgDB;

  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  // Efficiency configurator
  std::unordered_set<int> effEvtSet;
  bool effEvtSetReady = false;
  Configurable<bool> enableEffEvtSet{"enableEffEvtSet", true, "If true, MCGen uses the event-set built by MCReco; if false, MCGen runs stand-alone."};

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry spectraGen{"spectraGen", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry debugHistos{"debugHistos", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  HistogramRegistry evtimeHistos{"evtimeHistos", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  HistogramRegistry evLossHistos{"evLossHistos", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  HistogramRegistry histoGen{"histoGen", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};

  // Enable particle for analysis
  Configurable<bool> enablePr{"enablePr", true, "Flag to enable proton analysis."};
  Configurable<bool> enableDe{"enableDe", true, "Flag to enable deuteron analysis."};
  Configurable<bool> enableTr{"enableTr", true, "Flag to enable triton analysis."};
  Configurable<bool> enableHe{"enableHe", true, "Flag to enable helium-3 analysis."};
  Configurable<bool> enableAl{"enableAl", true, "Flag to enable alpha analysis."};

  Configurable<bool> enableTrackingEff{"enableTrackingEff", 0, "Flag to enable tracking efficiency histos."};
  Configurable<std::string> ccdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  // Set the triggered events skimming scheme
  struct : ConfigurableGroup {
    Configurable<bool> applySkimming{"applySkimming", false, "Skimmed dataset processing"};
    Configurable<std::string> cfgSkimming{"cfgSkimming", "fHe", "Configurable for skimming"};
  } skimmingOptions;

  // Set the event selection cuts
  struct : ConfigurableGroup {
    Configurable<bool> useSel8{"useSel8", true, "Use Sel8 for run3 Event Selection"};
    Configurable<bool> useTVXtrigger{"useTVXtrigger", false, "Use TVX for Event Selection (default w/ Sel8)"};
    Configurable<bool> removeTFBorder{"removeTFBorder", false, "Remove TimeFrame border (default w/ Sel8)"};
    Configurable<bool> removeITSROFBorder{"removeITSROFBorder", false, "Remove ITS Read-Out Frame border (default w/ Sel8)"};
  } evselOptions;

  // Set the multiplity event limits
  Configurable<float> cfgMultCutLow{"cfgMultCutLow", 0.0f, "Accepted multiplicity percentage lower limit"};
  Configurable<float> cfgMultCutHigh{"cfgMultCutHigh", 100.0f, "Accepted multiplicity percentage higher limit"};

  // Set the z-vertex event cut limits
  Configurable<float> cfgVzCutLow{"cfgVzCutLow", -10.0f, "Accepted z-vertex lower limit"};
  Configurable<float> cfgVzCutHigh{"cfgVzCutHigh", 10.0f, "Accepted z-vertex upper limit"};

  // Set the quality cuts for tracks
  struct : ConfigurableGroup {
    Configurable<bool> rejectFakeTracks{"rejectFakeTracks", false, "Flag to reject ITS-TPC fake tracks (for MC)"};
    Configurable<float> cfgCutITSClusters{"cfgCutITSClusters", -1.f, "Minimum number of ITS clusters"};
    Configurable<float> cfgCutTPCXRows{"cfgCutTPCXRows", -1.f, "Minimum number of crossed TPC rows"};
    Configurable<float> cfgCutTPCClusters{"cfgCutTPCClusters", -1.f, "Minimum number of found TPC clusters"};
    Configurable<float> cfgCutTPCCROFnd{"cfgCutTPCCROFnd", 0.8, "Minimum ratio of crossed TPC clusters over findable"};
    Configurable<std::vector<float>> tpcChi2NclCuts{"tpcChi2NclCuts", {0.5, 4}, "Range of accepted of Chi2/TPC clusters"};
    Configurable<std::vector<float>> itsChi2NclCuts{"itsChi2NclCuts", {0.f, 36}, "Range of accepted of Chi2/ITS clusters"};
    Configurable<int> nITSLayer{"nITSLayer", 0, "ITS Layer (0-6)"};
  } trkqcOptions;

  // Set the kinematic and PID cuts for tracks
  struct : ConfigurableGroup {
    Configurable<float> cfgMomentumCut{"cfgMomentumCut", 0.3f, "Value of the p selection for spectra (default 0.3)"};
    Configurable<float> cfgEtaCut{"cfgEtaCut", 0.8f, "Value of the eta selection for spectra (default 0.8)"};
    Configurable<float> cfgRapidityCutLow{"cfgRapidityCutLow", -1.0f, "Value of the low rapidity selection for spectra (default -1.0)"};
    Configurable<float> cfgRapidityCutHigh{"cfgRapidityCutHigh", 1.0f, "Value of the high rapidity selection for spectra (default 1.0)"};
  } kinemOptions;

  Configurable<bool> isPVContributorCut{"isPVContributorCut", false, "Flag to enable isPVContributor cut."};
  Configurable<bool> initITSPID{"initITSPID", false, "Flag to init the ITS PID response"};

  struct : ConfigurableGroup {
    Configurable<float> nsigmaTPCPr{"nsigmaTPCPr", 3.f, "Value of the Nsigma TPC cut for protons"};
    Configurable<float> nsigmaTPCDe{"nsigmaTPCDe", 3.f, "Value of the Nsigma TPC cut for deuterons"};
    Configurable<float> nsigmaTPCTr{"nsigmaTPCTr", 3.f, "Value of the Nsigma TPC cut for tritons"};
    Configurable<float> nsigmaTPCHe{"nsigmaTPCHe", 3.f, "Value of the Nsigma TPC cut for helium-3"};
    Configurable<float> nsigmaTPCAl{"nsigmaTPCAl", 3.f, "Value of the Nsigma TPC cut for alpha"};
  } nsigmaTPCvar;

  struct : ConfigurableGroup {
    Configurable<bool> useITSDeCut{"useITSDeCut", false, "Select Deuteron if compatible with deuteron hypothesis (via SigmaITS)"};
    Configurable<bool> useITSHeCut{"useITSHeCut", false, "Select Helium if compatible with helium hypothesis (via SigmaITS)"};
    Configurable<float> nsigmaITSDe{"nsigmaITSDe", -1.f, "Value of the Nsigma ITS cut for deuteron ( > nSigmaITSHe)"};
    Configurable<float> nsigmaITSHe{"nsigmaITSHe", -1.f, "Value of the Nsigma ITS cut for helium-3 ( > nSigmaITSHe)"};
    Configurable<bool> showAverageClusterSize{"showAverageClusterSize", false, "Show average cluster size"};
  } nsigmaITSvar;

  // Set additional cuts (used for debug)
  Configurable<float> cfgBetaCut{"cfgBetaCut", 0.4f, "Value of the beta selection for TOF cut (default 0.4)"};

  // Set the axis used in this task
  ConfigurableAxis binsPercentile{"binsPercentile", {100, 0, 100}, "Centrality FT0M"};
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.425, 0.45, 0.475, 0.5, 0.5125, 0.525, 0.5375, 0.55, 0.5625, 0.575, 0.5875, 0.6, 0.6125, 0.625, 0.6375, 0.65, 0.6625, 0.675, 0.6875, 0.7, 0.7125, 0.725, 0.7375, 0.75, 0.7625, 0.775, 0.7875, 0.8, 0.8125, 0.825, 0.8375, 0.85, 0.8625, 0.875, 0.8875, 0.9, 0.9125, 0.925, 0.9375, 0.95, 0.9625, 0.975, 0.9875, 1.0, 1.0125, 1.025, 1.0375, 1.05, 1.0625, 1.075, 1.0875, 1.1, 1.1125, 1.125, 1.1375, 1.15, 1.1625, 1.175, 1.1875, 1.2, 1.2125, 1.225, 1.2375, 1.25, 1.2625, 1.275, 1.2875, 1.3, 1.3125, 1.325, 1.3375, 1.35, 1.3625, 1.375, 1.3875, 1.4, 1.4125, 1.425, 1.4375, 1.45, 1.4625, 1.475, 1.4875, 1.5, 1.5125, 1.525, 1.5375, 1.55, 1.5625, 1.575, 1.5875, 1.6, 1.6125, 1.625, 1.6375, 1.65, 1.6625, 1.675, 1.6875, 1.7, 1.7125, 1.725, 1.7375, 1.75, 1.7625, 1.775, 1.7875, 1.8, 1.8125, 1.825, 1.8375, 1.85, 1.8625, 1.875, 1.8875, 1.9, 1.9125, 1.925, 1.9375, 1.95, 1.9625, 1.975, 1.9875, 2.0, 2.0625, 2.125, 2.1875, 2.25, 2.3125, 2.375, 2.4375, 2.5, 2.625, 2.75, 2.875, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0}, ""};
  ConfigurableAxis binsPtHe{"binsPtHe", {VARIABLE_WIDTH, 1.0, 1.25, 1.50, 1.75, 2.0, 2.25, 2.50, 2.75, 3.0, 3.25, 3.50, 3.75, 4.0, 4.50, 5.0, 6.0, 7.0, 8.0}, ""};
  ConfigurableAxis binsPtZHe{"binsPtZHe", {VARIABLE_WIDTH, 0.5, 0.625, 0.75, 0.875, 1.0, 1.125, 1.25, 1.375, 1.5, 1.625, 1.75, 1.875, 2.0, 2.25, 2.5, 3.0, 3.5, 4.0}, ""};

  ConfigurableAxis binsdEdx{"binsdEdx", {600, 0.f, 3000.f}, ""};
  ConfigurableAxis binsBeta{"binsBeta", {120, 0.0, 1.2}, ""};
  ConfigurableAxis binsDCA{"binsDCA", {400, -1.f, 1.f}, ""};
  ConfigurableAxis binsSigmaITS{"binsSigmaITS", {200, -20, 20}, ""};
  ConfigurableAxis binsSigmaTPC{"binsSigmaTPC", {1000, -100, 100}, ""};
  ConfigurableAxis binsSigmaTOF{"binsSigmaTOF", {1000, -100, 100}, ""};
  ConfigurableAxis binsMassPr{"binsMassPr", {100, -1., 1.f}, ""};
  ConfigurableAxis binsMassDe{"binsMassDe", {180, -1.8, 1.8f}, ""};
  ConfigurableAxis binsMassTr{"binsMassTr", {250, -2.5, 2.5f}, ""};
  ConfigurableAxis binsMassHe{"binsMassHe", {300, -3., 3.f}, ""};
  ConfigurableAxis avClsBins{"avClsBins", {200, 0, 20}, "Binning in average cluster size"};

  // Enable custom cuts/debug functions
  struct : ConfigurableGroup {
    Configurable<bool> enableFiltering{"enableFiltering", false, "Flag to enable filtering for p,d,t,He only -- disable if launch on skimmed dataset!"};
    Configurable<bool> enableIsGlobalTrack{"enableIsGlobalTrack", true, "Flag to enable IsGlobalTrackWoDCA"};
    Configurable<bool> enableEvTimeSplitting{"enableEvTimeSplitting", false, "Flag to enable histograms splitting depending on the Event Time used"};
  } filterOptions;

  Configurable<bool> enableCustomDCACut{"enableCustomDCACut", false, "Flag to enable DCA custom cuts - unflag to use standard isGlobalCut DCA cut"};
  struct : ConfigurableGroup {
    Configurable<int> cfgCustomDCA{"cfgCustomDCA", 0, "Select to use: pT independent DCAxy and DCAz CustomCut (0), pT dependent DCAxy and DCAz cut (1),  pt dependent DCAxy, DCAz CustomCut (2) DCAxy CustomCut, pT dependent DCAz (3) or a circular DCAxy,z cut  (4) for tracks. Need 'enableCustomDCACut' to be enabled."};
    Configurable<float> cfgCustomDCAxy{"cfgCustomDCAxy", 0.05f, "Value of the DCAxy selection for spectra (default 0.05 cm)"};
    Configurable<float> cfgCustomDCAz{"cfgCustomDCAz", 0.5f, "Value of the DCAz selection for spectra (default 0.5 cm)"};
  } dcaConfOptions;

  Configurable<std::vector<float>> parDCAxycuts{"parDCAxycuts", {0.004f, 0.013f, 1, 1}, "Parameters for Pt dependent DCAxy cut (if enabled): |DCAxy| < [3] * ([O] + [1]/Pt^[2])."};
  Configurable<std::vector<float>> parDCAzcuts{"parDCAzcuts", {0.004f, 0.013f, 1, 1}, "Parameters for Pt dependent DCAz cut (if enabled): |DCAz| < [3] * ([O] + [1]/Pt^[2])."};

  // Enable output histograms
  struct : ConfigurableGroup {
    Configurable<bool> makeDCABeforeCutPlots{"makeDCABeforeCutPlots", false, "Flag to enable plots of DCA before cuts"};
    Configurable<bool> makeDCAAfterCutPlots{"makeDCAAfterCutPlots", false, "Flag to enable plots of DCA after cuts"};
    Configurable<bool> makeFakeTracksPlots{"makeFakeTracksPlots", false, "Flag to enable plots of misidentified particles"};
    Configurable<bool> makeWrongEventPlots{"makeWrongEventPlots", false, "Flag to enable plots of particles from wrong event"};
    Configurable<bool> doTOFplots{"doTOFplots", true, "Flag to export plots of tracks with 1 hit on TOF."};
    Configurable<bool> enableExpSignalTPC{"enableExpSignalTPC", true, "Flag to export dEdX - dEdX(exp) plots."};
    Configurable<bool> enableExpSignalTOF{"enableExpSignalTOF", false, "Flag to export T - T(exp) plots."};
    Configurable<bool> enableBetaCut{"enableBetaCut", false, "Flag to enable TOF histograms with beta cut for debug"};
    Configurable<bool> enablePIDplot{"enablePIDplot", false, "Flag to enable PID histograms for debug"};
    Configurable<bool> enableEffPlots{"enableEffPlots", false, "Flag to enable histograms for efficiency debug."};
    Configurable<bool> enableNoTOFPlots{"enableNoTOFPlots", false, "Flag to enable histograms for TOF debug."};

  } outFlagOptions;

  Configurable<bool> enableDebug{"enableDebug", false, "Flag to enable histograms for debug"};

  Configurable<bool> usenITSLayer{"usenITSLayer", false, "Flag to enable ITS layer hit"};
  Configurable<int> useHasTRDConfig{"useHasTRDConfig", 0, "No selections on TRD (0); With TRD (1); Without TRD (2)"};
  Configurable<int> massTOFConfig{"massTOFConfig", 0, "Estimate massTOF using beta with (0) TPC momentum (1) TOF expected momentum (2) p momentum."};
  Configurable<int> helium3Pt{"helium3Pt", 0, "Select use default pT (0) or use instead 2*pT (1) for helium-3"};
  Configurable<int> unableDPtShift{"unableDPtShift", 0, "Select (0) to apply deuteron pT shift or (1) to use default pT."};
  Configurable<int> unableAntiDPtShift{"unableAntiDPtShift", 0, "Select (0) to apply antideuteron pT shift or (1) to use default pT."};

  // Additional function used for pT-shift calibration
  TF1* fShiftPtHe = 0;
  TF1* fShiftPtantiHe = 0;
  TF1* fShiftAntiD = 0;
  TF1* fShiftD = 0;

  Configurable<bool> enablePtShiftD{"enablePtShiftD", true, "Flag to enable Pt shift (for Deuteron only)"};
  Configurable<bool> enablePtShiftAntiD{"enablePtShiftAntiD", true, "Flag to enable Pt shift (for antiDeuteron only)"};
  Configurable<std::vector<float>> parShiftPtD{"parShiftPtD", {-0.0955412, 0.798164, -0.536111, 0.0887876, -1.11022e-13}, "Parameters for Pt shift (if enabled)."};
  Configurable<std::vector<float>> parShiftPtAntiD{"parShiftPtAntiD", {-0.0955412, 0.798164, -0.536111, 0.0887876, -1.11022e-13}, "Parameters for Pt shift (if enabled)."};

  Configurable<bool> enablePtShiftHe{"enablePtShiftHe", false, "Flag to enable Pt shift (for He only)"};
  Configurable<std::vector<float>> parShiftPtHe{"parShiftPtHe", {0.0f, 0.1f, 0.1f, 0.1f, 0.1f}, "Parameters for helium3-Pt shift (if enabled)."};
  Configurable<std::vector<float>> parShiftPtAntiHe{"parShiftPtAntiHe", {0.0f, 0.1f, 0.1f, 0.1f, 0.1f}, "Parameters for anti-helium3-Pt shift (if enabled)."};

  Configurable<bool> enableCentrality{"enableCentrality", true, "Flag to enable centrality 3D histos)"};

  // ITS to TPC - Fake hit loop
  static constexpr int kFakeLoop = 10; // Fixed O2Linter error
  // TPC low/high momentum range
  static constexpr float kCfgTpcClasses[] = {0.5f, 0.1f};
  static constexpr float kCfgKaonCut = 5.f;

  // PDG codes and masses used in this analysis
  static constexpr int PDGPion = PDG_t::kPiPlus;
  static constexpr int PDGKaon = PDG_t::kKPlus;
  static constexpr int PDGProton = PDG_t::kProton;
  static constexpr int PDGDeuteron = o2::constants::physics::Pdg::kDeuteron;
  static constexpr int PDGTriton = o2::constants::physics::Pdg::kTriton;
  static constexpr int PDGHelium = o2::constants::physics::Pdg::kHelium3;
  static constexpr int PDGAlpha = o2::constants::physics::Pdg::kAlpha;
  static constexpr int PDGHyperTriton = o2::constants::physics::Pdg::kHyperTriton;
  static constexpr float MassProtonVal = o2::constants::physics::MassProton;
  static constexpr float MassDeuteronVal = o2::constants::physics::MassDeuteron;
  static constexpr float MassTritonVal = o2::constants::physics::MassTriton;
  static constexpr float MassHeliumVal = o2::constants::physics::MassHelium3;
  static constexpr float MassAlphaVal = o2::constants::physics::MassAlpha;

  // PDG of Mothers
  static constexpr int kPdgMotherList[] = {
    PDG_t::kPiPlus,
    PDG_t::kKPlus,
    PDG_t::kK0Short,
    PDG_t::kNeutron,
    PDG_t::kProton,
    PDG_t::kLambda0,
    o2::constants::physics::Pdg::kDeuteron,
    o2::constants::physics::Pdg::kHelium3,
    o2::constants::physics::Pdg::kTriton,
    o2::constants::physics::Pdg::kHyperTriton,
    o2::constants::physics::Pdg::kAlpha};

  static constexpr int kNumMotherList = sizeof(kPdgMotherList) / sizeof(kPdgMotherList[0]);

  static constexpr const char* kMotherNames[kNumMotherList] = {
    "#pi^{+}",
    "K^{+}",
    "K^{0}_{S}",
    "n",
    "p",
    "#Lambda",
    "d",
    "He3",
    "t",
    "^{3}_{#Lambda}H",
    "He4"};

  static constexpr int kMaxNumMom = 2; // X: 0..4, overflow=5

  template <typename TrackType>
  float averageClusterSizeTrk(const TrackType& track)
  {
    return o2::aod::ITSResponse::averageClusterSize(track.itsClusterSizes());
  }

  float averageClusterSizePerCoslInv(uint32_t itsClusterSizes, float eta) { return o2::aod::ITSResponse::averageClusterSize(itsClusterSizes) * std::cosh(eta); }

  template <typename TrackType>
  float averageClusterSizePerCoslInv(const TrackType& track)
  {
    return averageClusterSizePerCoslInv(track.itsClusterSizes(), track.eta());
  }

  void initCCDB(o2::aod::BCsWithTimestamps::iterator const& bc)
  {
    if (skimmingOptions.applySkimming) {
      zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), skimmingOptions.cfgSkimming.value);
      zorro.populateHistRegistry(histos, bc.runNumber());
    }
  }

  void init(o2::framework::InitContext& context)
  {
    if (initITSPID) {
      o2::aod::ITSResponse::setParameters(context);
    }
    if (skimmingOptions.applySkimming) {
      zorroSummary.setObject(zorro.getZorroSummary());
    }

    effEvtSet.clear();
    effEvtSetReady = false;

    const AxisSpec pAxis{binsPt, "#it{p} (GeV/#it{c})"};
    const AxisSpec ptAxis{binsPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec ptHeAxis{binsPtHe, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec pZAxis{binsPt, "#it{p}/z (GeV/#it{c})"};
    const AxisSpec ptZHeAxis{binsPtZHe, "#it{p}_{T}/z (GeV/#it{c})"};
    const AxisSpec dedxAxis{binsdEdx, "d#it{E}/d#it{x} A.U."};
    const AxisSpec betaAxis{binsBeta, "TOF #beta"};
    const AxisSpec dcaxyAxis{binsDCA, "DCAxy (cm)"};
    const AxisSpec dcazAxis{binsDCA, "DCAz (cm)"};
    const AxisSpec massPrAxis{binsMassPr, ""};
    const AxisSpec massDeAxis{binsMassDe, ""};
    const AxisSpec massTrAxis{binsMassTr, ""};
    const AxisSpec massHeAxis{binsMassHe, ""};
    const AxisSpec sigmaITSAxis{binsSigmaITS, ""};
    const AxisSpec sigmaTPCAxis{binsSigmaTPC, ""};
    const AxisSpec sigmaTOFAxis{binsSigmaTOF, ""};
    const AxisSpec avClsAxis{avClsBins, "<ITS Cls. Size>"};
    const AxisSpec avClsEffAxis{avClsBins, "<ITS Cls. Size> / cosh(#eta)"};

    if (doprocessData == true && doprocessMCReco == true) {
      LOG(fatal) << "Can't enable processData and processMCReco in the same time, pick one!";
    }
    if (doprocessEvSgLossMC) {
      evLossHistos.add<TH1>("evLoss/hEvent", "Event loss histograms; ; counts", HistType::kTH1F, {{4, 0., 4.}});
      evLossHistos.get<TH1>(HIST("evLoss/hEvent"))->GetXaxis()->SetBinLabel(1, "All Gen.");
      evLossHistos.get<TH1>(HIST("evLoss/hEvent"))->GetXaxis()->SetBinLabel(2, "TVX (reco.)");
      evLossHistos.get<TH1>(HIST("evLoss/hEvent"))->GetXaxis()->SetBinLabel(3, "MC Sel8 (TVX + NoTFB) (reco.)");
      evLossHistos.get<TH1>(HIST("evLoss/hEvent"))->GetXaxis()->SetBinLabel(4, "Sel8 (reco.)");

      evLossHistos.add<TH1>("evLoss/pt/hDeuteronTriggeredTVX", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{100, 0., 5.}});
      evLossHistos.add<TH1>("evLoss/pt/hDeuteronTriggeredSel8", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{100, 0., 5.}});
      evLossHistos.add<TH1>("evLoss/pt/hDeuteronTriggeredMCSel8", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{100, 0., 5.}});
      evLossHistos.add<TH1>("evLoss/pt/hDeuteronGen", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{100, 0., 5.}});
      evLossHistos.add<TH1>("evLoss/pt/hAntiDeuteronTriggeredTVX", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{100, 0., 5.}});
      evLossHistos.add<TH1>("evLoss/pt/hAntiDeuteronTriggeredMCSel8", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{100, 0., 5.}});
      evLossHistos.add<TH1>("evLoss/pt/hAntiDeuteronTriggeredSel8", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{100, 0., 5.}});
      evLossHistos.add<TH1>("evLoss/pt/hAntiDeuteronGen", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{100, 0., 5.}});

      evLossHistos.add<TH1>("evLoss/pt/hHeliumTriggeredTVX", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{100, 0., 5.}});
      evLossHistos.add<TH1>("evLoss/pt/hHeliumTriggeredSel8", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{100, 0., 5.}});
      evLossHistos.add<TH1>("evLoss/pt/hHeliumTriggeredMCSel8", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{100, 0., 5.}});
      evLossHistos.add<TH1>("evLoss/pt/hHeliumGen", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{100, 0., 5.}});
      evLossHistos.add<TH1>("evLoss/pt/hAntiHeliumTriggeredTVX", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{100, 0., 5.}});
      evLossHistos.add<TH1>("evLoss/pt/hAntiHeliumTriggeredSel8", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{100, 0., 5.}});
      evLossHistos.add<TH1>("evLoss/pt/hAntiHeliumTriggeredMCSel8", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{100, 0., 5.}});
      evLossHistos.add<TH1>("evLoss/pt/hAntiHeliumGen", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{100, 0., 5.}});
    }
    if (doprocessMCRecoLfPidEv) {
      spectraGen.add<TH1>("LfEv/pT_nocut", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{ptHeAxis}});
      spectraGen.add<TH1>("LfEv/pT_TVXtrigger", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{ptHeAxis}});
      spectraGen.add<TH1>("LfEv/pT_TFrameBorder", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{ptHeAxis}});
      spectraGen.add<TH1>("LfEv/pT_ITSROFBorder", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{ptHeAxis}});
      spectraGen.add<TH1>("LfEv/pT_sel8", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{ptHeAxis}});
      spectraGen.add<TH1>("LfEv/pT_MCsel8", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{ptHeAxis}});

      spectraGen.add<TH1>("LfEv/helium/pT_nocut_He", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{ptHeAxis}});
      spectraGen.add<TH1>("LfEv/helium/pT_TVXtrigger_He", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{ptHeAxis}});
      spectraGen.add<TH1>("LfEv/helium/pT_TFrameBorder_He", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{ptHeAxis}});
      spectraGen.add<TH1>("LfEv/helium/pT_ITSROFBorder_He", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{ptHeAxis}});
      spectraGen.add<TH1>("LfEv/helium/pT_sel8_He", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{ptHeAxis}});
      spectraGen.add<TH1>("LfEv/helium/pT_MCsel8_He", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{ptHeAxis}});

      spectraGen.add<TH1>("LfEv/helium/prim/pT_nocut_He", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{ptHeAxis}});
      spectraGen.add<TH1>("LfEv/helium/prim/pT_TVXtrigger_He", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{ptHeAxis}});
      spectraGen.add<TH1>("LfEv/helium/prim/pT_TFrameBorder_He", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{ptHeAxis}});
      spectraGen.add<TH1>("LfEv/helium/prim/pT_ITSROFBorder_He", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{ptHeAxis}});
      spectraGen.add<TH1>("LfEv/helium/prim/pT_sel8_He", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{ptHeAxis}});
      spectraGen.add<TH1>("LfEv/helium/prim/pT_MCsel8_He", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{ptHeAxis}});

      spectraGen.add<TH1>("LfEv/helium/pT_nocut_antiHe", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{ptHeAxis}});
      spectraGen.add<TH1>("LfEv/helium/pT_TVXtrigger_antiHe", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{ptHeAxis}});
      spectraGen.add<TH1>("LfEv/helium/pT_TFrameBorder_antiHe", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{ptHeAxis}});
      spectraGen.add<TH1>("LfEv/helium/pT_ITSROFBorder_antiHe", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{ptHeAxis}});
      spectraGen.add<TH1>("LfEv/helium/pT_sel8_antiHe", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{ptHeAxis}});
      spectraGen.add<TH1>("LfEv/helium/pT_MCsel8_antiHe", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{ptHeAxis}});

      spectraGen.add<TH1>("LfEv/helium/prim/pT_nocut_antiHe", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{ptHeAxis}});
      spectraGen.add<TH1>("LfEv/helium/prim/pT_TVXtrigger_antiHe", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{ptHeAxis}});
      spectraGen.add<TH1>("LfEv/helium/prim/pT_TFrameBorder_antiHe", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{ptHeAxis}});
      spectraGen.add<TH1>("LfEv/helium/prim/pT_ITSROFBorder_antiHe", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{ptHeAxis}});
      spectraGen.add<TH1>("LfEv/helium/prim/pT_sel8_antiHe", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{ptHeAxis}});
      spectraGen.add<TH1>("LfEv/helium/prim/pT_MCsel8_antiHe", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{ptHeAxis}});
    }

    if (doprocessMCGenLosses) {
      histoGen.add("events/hMCGen", "hMCGen", {HistType::kTH1D, {{3, 0.f, 3.f}}});
      histoGen.get<TH1>(HIST("events/hMCGen"))->GetXaxis()->SetBinLabel(1, "All");
      histoGen.get<TH1>(HIST("events/hMCGen"))->GetXaxis()->SetBinLabel(2, "Vtz");
      histoGen.get<TH1>(HIST("events/hMCGen"))->GetXaxis()->SetBinLabel(3, "INELgt0");

      histoGen.add("events/hMCGenReco", "hMCGenReco", {HistType::kTH1D, {{2, 0.f, 2.f}}});
      histoGen.get<TH1>(HIST("events/hMCGenReco"))->GetXaxis()->SetBinLabel(1, "INEL");
      histoGen.get<TH1>(HIST("events/hMCGenReco"))->GetXaxis()->SetBinLabel(2, "INELgt0");

      histoGen.add("events/hMCReco", "hMCReco", {HistType::kTH1D, {{3, 0.f, 3.f}}});
      histoGen.get<TH1>(HIST("events/hMCReco"))->GetXaxis()->SetBinLabel(1, "All");
      histoGen.get<TH1>(HIST("events/hMCReco"))->GetXaxis()->SetBinLabel(2, "Ev sel passed");
      histoGen.get<TH1>(HIST("events/hMCReco"))->GetXaxis()->SetBinLabel(3, "INELgt0");

      histoGen.add("helium/MCGen/ptGen_INEL_Prim_He", "generated particles", HistType::kTH1F, {ptHeAxis});
      histoGen.add("helium/MCGen/ptGen_INEL_Prim_antiHe", "generated particles", HistType::kTH1F, {ptHeAxis});

      histoGen.add("helium/MCGen/ptGen_INELgt0_Prim_He", "generated particles", HistType::kTH1F, {ptHeAxis});
      histoGen.add("helium/MCGen/ptGen_INELgt0_Prim_antiHe", "generated particles", HistType::kTH1F, {ptHeAxis});

      histoGen.add("helium/MCGenReco/ptGen_INEL_Prim_He", "generated particles", HistType::kTH1F, {ptHeAxis});
      histoGen.add("helium/MCGenReco/ptGen_INEL_Prim_antiHe", "generated particles", HistType::kTH1F, {ptHeAxis});

      histoGen.add("helium/MCGenReco/ptGen_INELgt0_Prim_He", "generated particles", HistType::kTH1F, {ptHeAxis});
      histoGen.add("helium/MCGenReco/ptGen_INELgt0_Prim_antiHe", "generated particles", HistType::kTH1F, {ptHeAxis});

      histoGen.add("helium/MCReco/ptGen_INEL_Prim_He", "generated particles", HistType::kTH1F, {ptHeAxis});
      histoGen.add("helium/MCReco/ptGen_INEL_Prim_antiHe", "generated particles", HistType::kTH1F, {ptHeAxis});

      histoGen.add("helium/MCReco/ptGen_INELgt0_Prim_He", "generated particles", HistType::kTH1F, {ptHeAxis});
      histoGen.add("helium/MCReco/ptGen_INELgt0_Prim_antiHe", "generated particles", HistType::kTH1F, {ptHeAxis});

      if (enableCentrality) {
        histoGen.add("helium/MCGen/ptGenVsMult_INEL_Prim_He", "generated particles", HistType::kTH2F, {{ptHeAxis}, {binsPercentile}});
        histoGen.add("helium/MCGen/ptGenVsMult_INEL_Prim_antiHe", "generated particles", HistType::kTH2F, {{ptHeAxis}, {binsPercentile}});
        histoGen.add("helium/MCGen/ptGenVsMult_INELgt0_Prim_He", "generated particles", HistType::kTH2F, {{ptHeAxis}, {binsPercentile}});
        histoGen.add("helium/MCGen/ptGenVsMult_INELgt0_Prim_antiHe", "generated particles", HistType::kTH2F, {{ptHeAxis}, {binsPercentile}});

        histoGen.add("helium/MCGenReco/ptGenVsMult_INEL_Prim_He", "generated particles", HistType::kTH2F, {{ptHeAxis}, {binsPercentile}});
        histoGen.add("helium/MCGenReco/ptGenVsMult_INEL_Prim_antiHe", "generated particles", HistType::kTH2F, {{ptHeAxis}, {binsPercentile}});
        histoGen.add("helium/MCGenReco/ptGenVsMult_INELgt0_Prim_He", "generated particles", HistType::kTH2F, {{ptHeAxis}, {binsPercentile}});
        histoGen.add("helium/MCGenReco/ptGenVsMult_INELgt0_Prim_antiHe", "generated particles", HistType::kTH2F, {{ptHeAxis}, {binsPercentile}});

        histoGen.add("helium/MCReco/ptGenVsMult_INEL_Prim_He", "generated particles", HistType::kTH2F, {{ptHeAxis}, {binsPercentile}});
        histoGen.add("helium/MCReco/ptGenVsMult_INEL_Prim_antiHe", "generated particles", HistType::kTH2F, {{ptHeAxis}, {binsPercentile}});
        histoGen.add("helium/MCReco/ptGenVsMult_INELgt0_Prim_He", "generated particles", HistType::kTH2F, {{ptHeAxis}, {binsPercentile}});
        histoGen.add("helium/MCReco/ptGenVsMult_INELgt0_Prim_antiHe", "generated particles", HistType::kTH2F, {{ptHeAxis}, {binsPercentile}});
      }
    }

    if (enableDebug) {
      debugHistos.add<TH1>("qa/h1VtxZ_nocut", "V_{z};V_{z} (in cm); counts", HistType::kTH1F, {{1500, -15, 15}});
      debugHistos.add<TH1>("qa/h1VtxZ_TVXtrigger", "V_{z};V_{z} (in cm); counts", HistType::kTH1F, {{1500, -15, 15}});
      debugHistos.add<TH1>("qa/h1VtxZ_TFrameBorder", "V_{z};V_{z} (in cm); counts", HistType::kTH1F, {{1500, -15, 15}});
      debugHistos.add<TH1>("qa/h1VtxZ_ITSROFBorder", "V_{z};V_{z} (in cm); counts", HistType::kTH1F, {{1500, -15, 15}});
      debugHistos.add<TH1>("qa/h1VtxZ_sel8", "V_{z};V_{z} (in cm); counts", HistType::kTH1F, {{1500, -15, 15}});
      debugHistos.add<TH1>("qa/h1VtxZ_Centrality", "V_{z};V_{z} (in cm); counts", HistType::kTH1F, {{1500, -15, 15}});

      if (enableCentrality) {
        debugHistos.add<TH1>("event/hFT0M", "hFT0M", HistType::kTH1F, {{binsPercentile, "Centrality FT0M"}});
        debugHistos.add<TH1>("event/hFV0M", "hFV0M", HistType::kTH1F, {{binsPercentile, "Centrality FV0M"}});
      }
    }

    histos.add<TH1>("event/eventSkimming", "eventSkimming", HistType::kTH1D, {{2, 0.0, 2.0}});
    auto hSkim = histos.get<TH1>(HIST("event/eventSkimming"));
    hSkim->GetXaxis()->SetBinLabel(1, "Total");
    hSkim->GetXaxis()->SetBinLabel(2, "Skimmed events");

    if (enableCentrality) {
      histos.add<TH2>("event/eventSelection", "eventSelection", HistType::kTH2D, {{8, -0.5, 7.5}, {binsPercentile, "Centrality FT0M"}});
      auto h2d = histos.get<TH2>(HIST("event/eventSelection"));
      if (skimmingOptions.applySkimming)
        h2d->GetXaxis()->SetBinLabel(1, "Skimmed events");
      else
        h2d->GetXaxis()->SetBinLabel(1, "Total");

      h2d->GetXaxis()->SetBinLabel(2, "TVX trigger cut");
      h2d->GetXaxis()->SetBinLabel(3, "TF border cut");
      h2d->GetXaxis()->SetBinLabel(4, "ITS ROF cut");
      h2d->GetXaxis()->SetBinLabel(5, "TVX + TF + ITS ROF");
      h2d->GetXaxis()->SetBinLabel(6, "Sel8 cut");
      h2d->GetXaxis()->SetBinLabel(7, "Z-vert Cut");
      h2d->GetXaxis()->SetBinLabel(8, "Multiplicity cut");
    } else {
      histos.add<TH1>("event/eventSelection", "eventSelection", HistType::kTH1D, {{8, -0.5, 7.5}});
      auto h1d = histos.get<TH1>(HIST("event/eventSelection"));
      if (skimmingOptions.applySkimming)
        h1d->GetXaxis()->SetBinLabel(1, "Skimmed events");
      else
        h1d->GetXaxis()->SetBinLabel(1, "Total");

      h1d->GetXaxis()->SetBinLabel(2, "TVX trigger cut");
      h1d->GetXaxis()->SetBinLabel(3, "TF border cut");
      h1d->GetXaxis()->SetBinLabel(4, "ITS ROF cut");
      h1d->GetXaxis()->SetBinLabel(5, "TVX + TF + ITS ROF");
      h1d->GetXaxis()->SetBinLabel(6, "Sel8 cut");
      h1d->GetXaxis()->SetBinLabel(7, "Z-vert Cut");
      h1d->GetXaxis()->SetBinLabel(8, "Multiplicity cut");
    }

    if (enableCentrality)
      histos.add<TH2>("event/h1VtxZ", "V_{z};V_{z} (in cm); counts", HistType::kTH2F, {{1500, -15, 15}, {binsPercentile, "Centrality FT0M"}});
    else
      histos.add<TH1>("event/h1VtxZ", "V_{z};V_{z} (in cm); counts", HistType::kTH1F, {{1500, -15, 15}});

    if (outFlagOptions.enablePIDplot) {
      histos.add<TH1>("tracks/h1pT", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{500, 0., 10.}});
      histos.add<TH1>("tracks/h1p", "Track momentum; p (GeV/#it{c}); counts", HistType::kTH1F, {{500, 0., 10.}});
    }

    if (enableDebug) {
      histos.add<TH1>("qa/h1ITSncr", "number of crossed rows in ITS; ITSncr; counts", HistType::kTH1F, {{12, 0, 12}});
      histos.add<TH1>("qa/h1TPCncr", "number of crossed rows in TPC; TPCncr; counts", HistType::kTH1F, {{150, 60, 170}});
      histos.add<TH1>("qa/h1rTPC", "ratio of ncr over findable in TPC; rTPC; counts", HistType::kTH1F, {{200, 0.9, 1.8}});
      histos.add<TH1>("qa/h1TPCnfound", "ratio of found cluster in TPC; TPCnfound; counts", HistType::kTH1F, {{150, 60, 170}});
      histos.add<TH1>("qa/h1chi2ITS", "#chi^{2}_{ITS}/n_{ITS}; #chi^{2}_{ITS}/n_{ITS};counts", HistType::kTH1F, {{51, -0.5, 50.5}});
      histos.add<TH1>("qa/h1chi2TPC", "#chi^{2}_{TPC}/n_{TPC}; #chi^{2}_{TPC}/n_{TPC}; counts", HistType::kTH1F, {{11, -0.5, 10.5}});
    }

    if (outFlagOptions.enableEffPlots) {
      histos.add<TH2>("tracks/eff/h2pVsTPCmomentum", "#it{p}_{TPC} vs #it{p}; #it{p}_{TPC}; #it{p}", HistType::kTH2F, {{200, 0.f, 8.f}, {200, 0.f, 8.f}});
      if (outFlagOptions.doTOFplots)
        histos.add<TH2>("tracks/eff/h2TPCmomentumVsTOFExpMomentum", "#it{p}_{TOF} vs #it{p}_{TPC}; #it{p}_{TOF}; #it{p}_{TPC}", HistType::kTH2F, {{200, 0.f, 8.f}, {200, 0.f, 8.f}});

      if (enablePr) {
        histos.add<TH2>("tracks/eff/proton/h2pVsTPCmomentumPr", "#it{p}_{TPC} vs #it{p}; #it{p}_{TPC}; #it{p}", HistType::kTH2F, {{200, 0.f, 8.f}, {200, 0.f, 8.f}});
        histos.add<TH2>("tracks/eff/proton/h2pVsTPCmomentumantiPr", "#it{p}_{TPC} vs #it{p}; #it{p}_{TPC}; #it{p}", HistType::kTH2F, {{200, 0.f, 8.f}, {200, 0.f, 8.f}});
        if (outFlagOptions.doTOFplots) {
          histos.add<TH2>("tracks/eff/proton/h2pVsTOFExpMomentumPr", "#it{p}_{TOF} vs #it{p}; #it{p}_{TOF}; #it{p}", HistType::kTH2F, {{200, 0.f, 8.f}, {200, 0.f, 8.f}});
          histos.add<TH2>("tracks/eff/proton/h2pVsTOFExpMomentumantiPr", "#it{p}_{TOF} vs #it{p}; #it{p}_{TOF}; #it{p}", HistType::kTH2F, {{200, 0.f, 8.f}, {200, 0.f, 8.f}});
          histos.add<TH2>("tracks/eff/proton/h2TPCmomentumVsTOFExpMomentumPr", "#it{p}_{TOF} vs #it{p}_{TPC}; #it{p}_{TOF}; #it{p}_{TPC}", HistType::kTH2F, {{200, 0.f, 8.f}, {200, 0.f, 8.f}});
          histos.add<TH2>("tracks/eff/proton/h2TPCmomentumVsTOFExpMomentumantiPr", "#it{p}_{TOF} vs #it{p}_{TPC}; #it{p}_{TOF}; #it{p}_{TPC}", HistType::kTH2F, {{200, 0.f, 8.f}, {200, 0.f, 8.f}});
        }
      }
      if (enableDe) {
        histos.add<TH2>("tracks/eff/deuteron/h2pVsTPCmomentumDe", "#it{p}_{TPC} vs #it{p}; #it{p}_{TPC}; #it{p}", HistType::kTH2F, {{200, 0.f, 8.f}, {200, 0.f, 8.f}});
        histos.add<TH2>("tracks/eff/deuteron/h2pVsTPCmomentumantiDe", "#it{p}_{TPC} vs #it{p}; #it{p}_{TPC}; #it{p}", HistType::kTH2F, {{200, 0.f, 8.f}, {200, 0.f, 8.f}});
        if (outFlagOptions.doTOFplots) {
          histos.add<TH2>("tracks/eff/deuteron/h2pVsTOFExpMomentumDe", "#it{p}_{TOF} vs #it{p}; #it{p}_{TOF}; #it{p}", HistType::kTH2F, {{200, 0.f, 8.f}, {200, 0.f, 8.f}});
          histos.add<TH2>("tracks/eff/deuteron/h2pVsTOFExpMomentumantiDe", "#it{p}_{TOF} vs #it{p}; #it{p}_{TOF}; #it{p}", HistType::kTH2F, {{200, 0.f, 8.f}, {200, 0.f, 8.f}});
          histos.add<TH2>("tracks/eff/deuteron/h2TPCmomentumVsTOFExpMomentumDe", "#it{p}_{TOF} vs #it{p}_{TPC}; #it{p}_{TOF}; #it{p}_{TPC}", HistType::kTH2F, {{200, 0.f, 8.f}, {200, 0.f, 8.f}});
          histos.add<TH2>("tracks/eff/deuteron/h2TPCmomentumVsTOFExpMomentumantiDe", "#it{p}_{TOF} vs #it{p}_{TPC}; #it{p}_{TOF}; #it{p}_{TPC}", HistType::kTH2F, {{200, 0.f, 8.f}, {200, 0.f, 8.f}});
        }
      }
      if (enableTr) {
        histos.add<TH2>("tracks/eff/triton/h2pVsTPCmomentumTr", "#it{p}_{TPC} vs #it{p}; #it{p}_{TPC}; #it{p}", HistType::kTH2F, {{200, 0.f, 8.f}, {200, 0.f, 8.f}});
        histos.add<TH2>("tracks/eff/triton/h2pVsTPCmomentumantiTr", "#it{p}_{TPC} vs #it{p}; #it{p}_{TPC}; #it{p}", HistType::kTH2F, {{200, 0.f, 8.f}, {200, 0.f, 8.f}});
        if (outFlagOptions.doTOFplots) {
          histos.add<TH2>("tracks/eff/triton/h2pVsTOFExpMomentumTr", "#it{p}_{TOF} vs #it{p}; #it{p}_{TOF}; #it{p}", HistType::kTH2F, {{200, 0.f, 8.f}, {200, 0.f, 8.f}});
          histos.add<TH2>("tracks/eff/triton/h2pVsTOFExpMomentumantiTr", "#it{p}_{TOF} vs #it{p}; #it{p}_{TOF}; #it{p}", HistType::kTH2F, {{200, 0.f, 8.f}, {200, 0.f, 8.f}});
          histos.add<TH2>("tracks/eff/triton/h2TPCmomentumVsTOFExpMomentumTr", "#it{p}_{TOF} vs #it{p}_{TPC}; #it{p}_{TOF}; #it{p}_{TPC}", HistType::kTH2F, {{200, 0.f, 8.f}, {200, 0.f, 8.f}});
          histos.add<TH2>("tracks/eff/triton/h2TPCmomentumVsTOFExpMomentumantiTr", "#it{p}_{TOF} vs #it{p}_{TPC}; #it{p}_{TOF}; #it{p}_{TPC}", HistType::kTH2F, {{200, 0.f, 8.f}, {200, 0.f, 8.f}});
        }
      }
      if (enableHe) {
        histos.add<TH2>("tracks/eff/helium/h2pVsTPCmomentumHe", "#it{p}_{TPC} vs #it{p}; #it{p}_{TPC}; #it{p}", HistType::kTH2F, {{200, 0.f, 8.f}, {200, 0.f, 8.f}});
        histos.add<TH2>("tracks/eff/helium/h2pVsTPCmomentumantiHe", "#it{p}_{TPC} vs #it{p}; #it{p}_{TPC}; #it{p}", HistType::kTH2F, {{200, 0.f, 8.f}, {200, 0.f, 8.f}});
        if (outFlagOptions.doTOFplots) {
          histos.add<TH2>("tracks/eff/helium/h2pVsTOFExpMomentumHe", "#it{p}_{TOF} vs #it{p}; #it{p}_{TOF}; #it{p}", HistType::kTH2F, {{200, 0.f, 8.f}, {200, 0.f, 8.f}});
          histos.add<TH2>("tracks/eff/helium/h2pVsTOFExpMomentumantiHe", "#it{p}_{TOF} vs #it{p}; #it{p}_{TOF}; #it{p}", HistType::kTH2F, {{200, 0.f, 8.f}, {200, 0.f, 8.f}});
          histos.add<TH2>("tracks/eff/helium/h2TPCmomentumVsTOFExpMomentumHe", "#it{p}_{TOF} vs #it{p}_{TPC}; #it{p}_{TOF}; #it{p}_{TPC}", HistType::kTH2F, {{200, 0.f, 8.f}, {200, 0.f, 8.f}});
          histos.add<TH2>("tracks/eff/helium/h2TPCmomentumVsTOFExpMomentumantiHe", "#it{p}_{TOF} vs #it{p}_{TPC}; #it{p}_{TOF}; #it{p}_{TPC}", HistType::kTH2F, {{200, 0.f, 8.f}, {200, 0.f, 8.f}});
        }
      }
    }

    if (enableDebug) {
      debugHistos.add<TH1>("debug/event/h1CentV0M", "V0M; Multiplicity; counts", HistType::kTH1F, {{27000, 0, 27000}});
      // trackQA
      debugHistos.add<TH1>("debug/tracks/h1Eta", "pseudoRapidity; #eta; counts", HistType::kTH1F, {{200, -2.0, 2.0}});
      debugHistos.add<TH1>("debug/tracks/h1VarPhi", "#phi; #phi; counts", HistType::kTH1F, {{63, 0.0, 6.3}});
      debugHistos.add<TH2>("debug/tracks/h2EtaVsPhi", "#eta vs #phi; #eta; #phi", HistType::kTH2F, {{200, -2.0, 2.0}, {63, 0.0, 6.3}});
      debugHistos.add<TH2>("debug/tracks/h2PionYvsPt", "#it{y} vs #it{p}_{T} (#pi)", HistType::kTH2F, {{200, -2.0, 2.0}, {ptAxis}});
    }

    if (outFlagOptions.enableEffPlots) {
      if (enableDebug) {
        debugHistos.add<TH1>("tracks/eff/hPtP", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{400, 0., 8.}});
        debugHistos.add<TH1>("tracks/eff/hPtantiP", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{400, 0., 8.}});
        if (outFlagOptions.doTOFplots) {
          debugHistos.add<TH1>("tracks/eff/hPtPTOF", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{400, 0., 8.}});
          debugHistos.add<TH1>("tracks/eff/hPtantiPTOF", "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{400, 0., 8.}});
        }
      }

      if (enablePr) {
        histos.add<TH1>("tracks/eff/proton/hPtPr", "Track #it{p}_{T} (p); #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{400, 0., 8.}});
        histos.add<TH1>("tracks/eff/proton/hPtantiPr", "Track #it{p}_{T} (#bar{p}); #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{400, 0., 8.}});
        histos.add<TH1>("tracks/eff/proton/hPtPrTOF", "Track #it{p}_{T} (p); #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{400, 0., 8.}});
        histos.add<TH1>("tracks/eff/proton/hPtantiPrTOF", "Track #it{p}_{T} (#bar{p}); #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{400, 0., 8.}});
      }
      if (enableDe) {
        histos.add<TH1>("tracks/eff/deuteron/hPtDe", "Track #it{p}_{T} (d); #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{400, 0., 8.}});
        histos.add<TH1>("tracks/eff/deuteron/hPtantiDe", "Track #it{p}_{T} (#bar{d}); #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{400, 0., 8.}});
        histos.add<TH1>("tracks/eff/deuteron/hPtDeTOF", "Track #it{p}_{T} (d); #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{400, 0., 8.}});
        histos.add<TH1>("tracks/eff/deuteron/hPtantiDeTOF", "Track #it{p}_{T} (#bar{d}); #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{400, 0., 8.}});
      }
      if (enableTr) {
        histos.add<TH1>("tracks/eff/triton/hPtTr", "Track #it{p}_{T} (t); #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{400, 0., 8.}});
        histos.add<TH1>("tracks/eff/triton/hPtantiTr", "Track #it{p}_{T} (#bar{t}); #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{400, 0., 8.}});
        histos.add<TH1>("tracks/eff/triton/hPtTrTOF", "Track #it{p}_{T} (t); #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{400, 0., 8.}});
        histos.add<TH1>("tracks/eff/triton/hPtantiTrTOF", "Track #it{p}_{T} (#bar{t}); #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{400, 0., 8.}});
      }
      if (enableHe) {
        histos.add<TH1>("tracks/eff/helium/hPtHe", "Track #it{p}_{T} (He); #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{400, 0., 8.}});
        histos.add<TH1>("tracks/eff/helium/hPtantiHe", "Track #it{p}_{T} (#bar{He}); #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{400, 0., 8.}});
        histos.add<TH1>("tracks/eff/helium/hPtHeTOF", "Track #it{p}_{T} (He); #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{400, 0., 8.}});
        histos.add<TH1>("tracks/eff/helium/hPtantiHeTOF", "Track #it{p}_{T} (#bar{He}); #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{400, 0., 8.}});
      }
    }
    // tracks
    // DCAxy,z
    if (outFlagOptions.makeDCABeforeCutPlots) {
      histos.add<TH3>("tracks/dca/before/hDCAxyVsDCAzVsPt", "DCAxy vs DCAz vs Pt/z; DCAxy; DCAz", HistType::kTH3F, {{140, -0.7f, 0.7f}, {160, -0.8f, 0.8f}, {binsPtHe}});

      histos.add<TH2>("tracks/dca/before/hDCAxyVsDCAz", "DCAxy vs DCAz (before cuts)", HistType::kTH2F, {{550, -1.1, 1.1}, {550, -1.1, 1.1}});
      histos.add<TH2>("tracks/dca/before/hDCAxyVsPt", "DCAxy vs Pt", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
      histos.add<TH2>("tracks/dca/before/hDCAzVsPt", "DCAz vs Pt", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

      if (enablePr) {
        histos.add<TH2>("tracks/proton/dca/before/hDCAxyVsPtProton", "DCAxy vs Pt (p)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
        histos.add<TH2>("tracks/proton/dca/before/hDCAxyVsPtantiProton", "DCAxy vs Pt (#bar{p})", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
        histos.add<TH2>("tracks/proton/dca/before/hDCAzVsPtProton", "DCAz vs Pt (p)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
        histos.add<TH2>("tracks/proton/dca/before/hDCAzVsPtantiProton", "DCAz vs Pt (#bar{p})", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
      }
      if (enableDe) {
        if (enableCentrality) {
          histos.add<TH3>("tracks/deuteron/dca/before/hDCAxyVsPtDeuteronVsMult", "DCAxy vs Pt (d)", HistType::kTH3F, {{ptAxis}, {dcaxyAxis}, {binsPercentile}});
          histos.add<TH3>("tracks/deuteron/dca/before/hDCAxyVsPtantiDeuteronVsMult", "DCAxy vs Pt (#bar{d})", HistType::kTH3F, {{ptAxis}, {dcaxyAxis}, {binsPercentile}});
        } else {
          histos.add<TH2>("tracks/deuteron/dca/before/hDCAxyVsPtDeuteron", "DCAxy vs Pt (d)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/deuteron/dca/before/hDCAxyVsPtantiDeuteron", "DCAxy vs Pt (#bar{d})", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
        }
        histos.add<TH2>("tracks/deuteron/dca/before/hDCAzVsPtDeuteron", "DCAz vs Pt (d)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
        histos.add<TH2>("tracks/deuteron/dca/before/hDCAzVsPtantiDeuteron", "DCAz vs Pt (#bar{d})", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

        if (outFlagOptions.enableNoTOFPlots) {
          histos.add<TH2>("tracks/deuteron/dca/before/hDCAxyVsPtDeuteronNoTOF", "DCAxy vs Pt (d)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/deuteron/dca/before/hDCAxyVsPtantiDeuteronNoTOF", "DCAxy vs Pt (#bar{d})", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/deuteron/dca/before/hDCAzVsPtDeuteronNoTOF", "DCAz vs Pt (d)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/deuteron/dca/before/hDCAzVsPtantiDeuteronNoTOF", "DCAz vs Pt (#bar{d})", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
        }
      }
      if (enableTr) {
        histos.add<TH2>("tracks/triton/dca/before/hDCAxyVsPtTriton", "DCAxy vs Pt (t)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
        histos.add<TH2>("tracks/triton/dca/before/hDCAxyVsPtantiTriton", "DCAxy vs Pt (#bar{t})", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
        histos.add<TH2>("tracks/triton/dca/before/hDCAzVsPtTriton", "DCAz vs Pt (t)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
        histos.add<TH2>("tracks/triton/dca/before/hDCAzVsPtantiTriton", "DCAz vs Pt (#bar{t})", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
      }
      if (enableHe) {
        histos.add("tracks/helium/dca/before/h3DCAvsPtHelium", "", HistType::kTHnSparseD, {dcaxyAxis, dcazAxis, ptZHeAxis});
        histos.add("tracks/helium/dca/before/h3DCAvsPtantiHelium", "", HistType::kTHnSparseD, {dcaxyAxis, dcazAxis, ptZHeAxis});

        histos.add<TH2>("tracks/helium/dca/before/hDCAxyVsPtHelium", "DCAxy vs Pt (He)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});
        histos.add<TH2>("tracks/helium/dca/before/hDCAxyVsPtantiHelium", "DCAxy vs Pt (#bar{He})", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});
        histos.add<TH2>("tracks/helium/dca/before/hDCAzVsPtHelium", "DCAz vs Pt (He)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
        histos.add<TH2>("tracks/helium/dca/before/hDCAzVsPtantiHelium", "DCAz vs Pt (#bar{He})", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});

        if (outFlagOptions.enableNoTOFPlots) {
          histos.add<TH2>("tracks/helium/dca/before/hDCAxyVsPtHeliumNoTOF", "DCAxy vs Pt (He)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/helium/dca/before/hDCAxyVsPtantiHeliumNoTOF", "DCAxy vs Pt (#bar{He})", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/helium/dca/before/hDCAzVsPtHeliumNoTOF", "DCAz vs Pt (He)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/helium/dca/before/hDCAzVsPtantiHeliumNoTOF", "DCAz vs Pt (#bar{He})", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
        }

        if (outFlagOptions.doTOFplots) {
          histos.add<TH3>("tracks/helium/dca/before/TOF/hDCAxyVsDCAzVsPtHelium", "DCAxy vs DCAz vs Pt/z (He) (w/TOF); DCAxy; DCAz", HistType::kTH3F, {{140, -0.7f, 0.7f}, {160, -0.8f, 0.8f}, {binsPtZHe}});
          histos.add<TH3>("tracks/helium/dca/before/TOF/hDCAxyVsDCAzVsPtantiHelium", "DCAxy vs DCAz vs Pt/z (#bar{He}) (w/TOF); DCAxy; DCAz", HistType::kTH3F, {{140, -0.7f, 0.7f}, {160, -0.8f, 0.8f}, {binsPtZHe}});
          histos.add<TH2>("tracks/helium/dca/before/TOF/hDCAxyVsPtHelium", "DCAxy vs Pt (He)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/helium/dca/before/TOF/hDCAxyVsPtantiHelium", "DCAxy vs Pt (#bar{He})", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/helium/dca/before/TOF/hDCAzVsPtHelium", "DCAz vs Pt (He)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/helium/dca/before/TOF/hDCAzVsPtantiHelium", "DCAz vs Pt (#bar{He})", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
        }
      }
      if (enableAl) {
        histos.add<TH2>("tracks/alpha/dca/before/hDCAxyVsPtAlpha", "DCAxy vs Pt (#alpha)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
        histos.add<TH2>("tracks/alpha/dca/before/hDCAxyVsPtantiAlpha", "DCAxy vs Pt (#bar{#alpha})", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
        histos.add<TH2>("tracks/alpha/dca/before/hDCAzVsPtAlpha", "DCAz vs Pt (#alpha)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
        histos.add<TH2>("tracks/alpha/dca/before/hDCAzVsPtantiAlpha", "DCAz vs Pt (#bar{#alpha})", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
      }
    }

    if (outFlagOptions.makeDCAAfterCutPlots) {
      histos.add<TH1>("tracks/dca/after/hDCAxy", "DCAxy; DCAxy; counts", HistType::kTH1F, {{dcaxyAxis}});
      histos.add<TH1>("tracks/dca/after/hDCAz", "DCAz; DCAz; counts", HistType::kTH1F, {{dcazAxis}});
      histos.add<TH2>("tracks/dca/after/hDCAxyVsPt", "DCAxy vs Pt", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
      histos.add<TH2>("tracks/dca/after/hDCAzVsPt", "DCAz vs Pt", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
    }
    if (enablePr && outFlagOptions.makeDCAAfterCutPlots) {
      histos.add<TH2>("tracks/proton/dca/after/hDCAxyVsPtProton", "DCAxy vs Pt (p)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
      histos.add<TH2>("tracks/proton/dca/after/hDCAxyVsPtantiProton", "DCAxy vs Pt (#bar{p})", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
      histos.add<TH2>("tracks/proton/dca/after/hDCAzVsPtProton", "DCAz vs Pt (p)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
      histos.add<TH2>("tracks/proton/dca/after/hDCAzVsPtantiProton", "DCAz vs Pt (#bar{p})", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
    }
    if (enableDe && outFlagOptions.makeDCAAfterCutPlots) {
      histos.add<TH2>("tracks/deuteron/dca/after/hDCAxyVsPtDeuteron", "DCAxy vs Pt (d)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
      histos.add<TH2>("tracks/deuteron/dca/after/hDCAxyVsPtantiDeuteron", "DCAxy vs Pt (#bar{d})", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
      histos.add<TH2>("tracks/deuteron/dca/after/hDCAzVsPtDeuteron", "DCAz vs Pt (d)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
      histos.add<TH2>("tracks/deuteron/dca/after/hDCAzVsPtantiDeuteron", "DCAz vs Pt (#bar{d})", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
    }
    if (enableTr && outFlagOptions.makeDCAAfterCutPlots) {
      histos.add<TH2>("tracks/triton/dca/after/hDCAxyVsPtTriton", "DCAxy vs Pt (t)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
      histos.add<TH2>("tracks/triton/dca/after/hDCAxyVsPtantiTriton", "DCAxy vs Pt (#bar{t})", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
      histos.add<TH2>("tracks/triton/dca/after/hDCAzVsPtTriton", "DCAz vs Pt (t)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
      histos.add<TH2>("tracks/triton/dca/after/hDCAzVsPtantiTriton", "DCAz vs Pt (#bar{t})", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
    }
    if (enableHe && outFlagOptions.makeDCAAfterCutPlots) {
      histos.add<TH3>("tracks/helium/dca/after/hDCAxyVsDCAzVsPtHelium", "DCAxy vs DCAz vs Pt/z (He); DCAxy; DCAz", HistType::kTH3F, {{140, -0.7f, 0.7f}, {160, -0.8f, 0.8f}, {binsPtZHe}});
      histos.add<TH3>("tracks/helium/dca/after/hDCAxyVsDCAzVsPtantiHelium", "DCAxy vs DCAz vs Pt/z (#bar{He}); DCAxy; DCAz", HistType::kTH3F, {{140, -0.7f, 0.7f}, {160, -0.8f, 0.8f}, {binsPtZHe}});
      histos.add("tracks/helium/dca/after/h3DCAvsPtHelium", "", HistType::kTHnSparseD, {dcaxyAxis, dcazAxis, ptZHeAxis});
      histos.add("tracks/helium/dca/after/h3DCAvsPtantiHelium", "", HistType::kTHnSparseD, {dcaxyAxis, dcazAxis, ptZHeAxis});
      histos.add<TH2>("tracks/helium/dca/after/hDCAxyVsPtHelium", "DCAxy vs Pt (He)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});
      histos.add<TH2>("tracks/helium/dca/after/hDCAxyVsPtantiHelium", "DCAxy vs Pt (#bar{He})", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});
      histos.add<TH2>("tracks/helium/dca/after/hDCAzVsPtHelium", "DCAz vs Pt (He)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
      histos.add<TH2>("tracks/helium/dca/after/hDCAzVsPtantiHelium", "DCAz vs Pt (#bar{He})", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});

      if (outFlagOptions.doTOFplots) {
        histos.add<TH3>("tracks/helium/dca/after/TOF/hDCAxyVsDCAzVsPtHelium", "DCAxy vs DCAz vs Pt/z (He); DCAxy; DCAz", HistType::kTH3F, {{140, -0.7f, 0.7f}, {160, -0.8f, 0.8f}, {binsPtZHe}});
        histos.add<TH3>("tracks/helium/dca/after/TOF/hDCAxyVsDCAzVsPtantiHelium", "DCAxy vs DCAz vs Pt/z (#bar{He}); DCAxy; DCAz", HistType::kTH3F, {{140, -0.7f, 0.7f}, {160, -0.8f, 0.8f}, {binsPtZHe}});
        histos.add<TH2>("tracks/helium/dca/after/TOF/hDCAxyVsPtHelium", "DCAxy vs Pt (He)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});
        histos.add<TH2>("tracks/helium/dca/after/TOF/hDCAxyVsPtantiHelium", "DCAxy vs Pt (#bar{He})", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});
        histos.add<TH2>("tracks/helium/dca/after/TOF/hDCAzVsPtHelium", "DCAz vs Pt (He)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
        histos.add<TH2>("tracks/helium/dca/after/TOF/hDCAzVsPtantiHelium", "DCAz vs Pt (#bar{He})", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
      }
    }
    if (enableAl && outFlagOptions.makeDCAAfterCutPlots) {
      histos.add<TH2>("tracks/alpha/dca/after/hDCAxyVsPtAlpha", "DCAxy vs Pt (#alpha)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
      histos.add<TH2>("tracks/alpha/dca/after/hDCAxyVsPtantiAlpha", "DCAxy vs Pt (#bar{#alpha})", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
      histos.add<TH2>("tracks/alpha/dca/after/hDCAzVsPtAlpha", "DCAz vs Pt (#alpha)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
      histos.add<TH2>("tracks/alpha/dca/after/hDCAzVsPtantiAlpha", "DCAz vs Pt (#bar{#alpha})", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
    }

    // Spectra (After cuts)
    if (enablePr) {
      histos.add<TH1>("tracks/proton/h1ProtonSpectra", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
      histos.add<TH1>("tracks/proton/h1antiProtonSpectra", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});

      histos.add<TH2>("tracks/proton/h2ProtonYvsPt", "#it{y} vs #it{p}_{T} (p)", HistType::kTH2F, {{96, -1.2, 1.2}, {ptAxis}});
      histos.add<TH2>("tracks/proton/h2antiProtonYvsPt", "#it{y} vs #it{p}_{T} (#bar{p})", HistType::kTH2F, {{96, -1.2, 1.2}, {ptAxis}});
      histos.add<TH2>("tracks/proton/h2ProtonEtavsPt", "#it{#eta} vs #it{p}_{T} (p)", HistType::kTH2F, {{96, -1.2, 1.2}, {ptAxis}});
      histos.add<TH2>("tracks/proton/h2antiProtonEtavsPt", "#it{#eta} vs #it{p}_{T} (#bar{p})", HistType::kTH2F, {{96, -1.2, 1.2}, {ptAxis}});
    }
    if (enableDe) {
      histos.add<TH1>("tracks/deuteron/h1DeuteronSpectra", "#it{p}_{T} (d)", HistType::kTH1F, {ptAxis});
      histos.add<TH1>("tracks/deuteron/h1antiDeuteronSpectra", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});

      histos.add<TH2>("tracks/deuteron/h2DeuteronYvsPt", "#it{y} vs #it{p}_{T} (d)", HistType::kTH2F, {{96, -1.2, 1.2}, {ptAxis}});
      histos.add<TH2>("tracks/deuteron/h2antiDeuteronYvsPt", "#it{y} vs #it{p}_{T} (#bar{d})", HistType::kTH2F, {{96, -1.2, 1.2}, {ptAxis}});
      histos.add<TH2>("tracks/deuteron/h2DeuteronEtavsPt", "#it{#eta} vs #it{p}_{T} (d)", HistType::kTH2F, {{96, -1.2, 1.2}, {ptAxis}});
      histos.add<TH2>("tracks/deuteron/h2antiDeuteronEtavsPt", "#it{#eta} vs #it{p}_{T} (#bar{d})", HistType::kTH2F, {{96, -1.2, 1.2}, {ptAxis}});
    }
    if (enableTr) {
      histos.add<TH1>("tracks/triton/h1TritonSpectra", "#it{p}_{T} (t)", HistType::kTH1F, {ptAxis});
      histos.add<TH1>("tracks/triton/h1antiTritonSpectra", "#it{p}_{T} (#bar{t})", HistType::kTH1F, {ptAxis});
    }
    if (enableHe) {
      histos.add<TH2>("tracks/helium/h2HeliumYvsPt_Z2", "#it{y} vs #it{p}_{T} (He)", HistType::kTH2F, {{96, -1.2, 1.2}, {ptHeAxis}});
      histos.add<TH2>("tracks/helium/h2HeliumEtavsPt_Z2", "#it{#eta} vs #it{p}_{T} (He)", HistType::kTH2F, {{96, -1.2, 1.2}, {ptHeAxis}});

      histos.add<TH1>("tracks/helium/h1HeliumSpectra_Z2", "#it{p}_{T} (He)", HistType::kTH1F, {ptHeAxis});
      histos.add<TH1>("tracks/helium/h1antiHeliumSpectra_Z2", "#it{p}_{T} (#bar{He})", HistType::kTH1F, {ptHeAxis});

      histos.add<TH2>("tracks/helium/h2antiHeliumYvsPt_Z2", "#it{y} vs #it{p}_{T} (#bar{He})", HistType::kTH2F, {{96, -1.2, 1.2}, {ptHeAxis}});
      histos.add<TH2>("tracks/helium/h2antiHeliumEtavsPt_Z2", "#it{#eta} vs #it{p}_{T} (#bar{He})", HistType::kTH2F, {{96, -1.2, 1.2}, {ptHeAxis}});
    }
    if (enableAl) {
      histos.add<TH1>("tracks/alpha/h1AlphaSpectra", "#it{p}_{T} (#alpha)", HistType::kTH1F, {ptAxis});
      histos.add<TH1>("tracks/alpha/h1antiAlphaSpectra", "#it{p}_{T} (#bar{#alpha})", HistType::kTH1F, {ptAxis});
    }
    if (doprocessMCReco || doprocessMCRecoLfPid || doprocessMCRecoFiltered || doprocessMCRecoFilteredLight) {

      histos.add<TH1>("tracks/hItsDeHeChecker", "d and {}^{3}He counters", HistType::kTH1F, {{4, -0.5, 3.5}});
      histos.get<TH1>(HIST("tracks/hItsDeHeChecker"))->GetXaxis()->SetBinLabel(1, "totDe");
      histos.get<TH1>(HIST("tracks/hItsDeHeChecker"))->GetXaxis()->SetBinLabel(2, "totHe");
      histos.get<TH1>(HIST("tracks/hItsDeHeChecker"))->GetXaxis()->SetBinLabel(3, "keptDe");
      histos.get<TH1>(HIST("tracks/hItsDeHeChecker"))->GetXaxis()->SetBinLabel(4, "keptHe");

      // inclusive production
      if (enableTrackingEff) {
        debugHistos.add<TH1>("tracks/trackingEff/h1_its", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
        debugHistos.add<TH1>("tracks/trackingEff/h1_tpc", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
        debugHistos.add<TH1>("tracks/trackingEff/h1_its_tpc", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
        debugHistos.add<TH1>("tracks/trackingEff/h1_its_tpc_tof", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
        debugHistos.add<TH1>("tracks/trackingEff/h1_its_tpc_trd", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
        debugHistos.add<TH1>("tracks/trackingEff/h1_its_tpc_trd_tof", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});

        debugHistos.add<TH1>("tracks/trackingEff/h1GoodHit_its", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
        debugHistos.add<TH1>("tracks/trackingEff/h1GoodHit_tpc", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
        debugHistos.add<TH1>("tracks/trackingEff/h1GoodHit_its_tpc", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
        debugHistos.add<TH1>("tracks/trackingEff/h1GoodHit_its_tpc_tof", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
        debugHistos.add<TH1>("tracks/trackingEff/h1GoodHit_its_tpc_trd", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
        debugHistos.add<TH1>("tracks/trackingEff/h1GoodHit_its_tpc_trd_tof", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});

        debugHistos.add<TH1>("tracks/trackingEff/h1anti_its", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});
        debugHistos.add<TH1>("tracks/trackingEff/h1anti_tpc", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});
        debugHistos.add<TH1>("tracks/trackingEff/h1anti_its_tpc", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});
        debugHistos.add<TH1>("tracks/trackingEff/h1anti_its_tpc_tof", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});
        debugHistos.add<TH1>("tracks/trackingEff/h1anti_its_tpc_trd", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});
        debugHistos.add<TH1>("tracks/trackingEff/h1anti_its_tpc_trd_tof", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});

        debugHistos.add<TH1>("tracks/trackingEff/h1antiGoodHit_its", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});
        debugHistos.add<TH1>("tracks/trackingEff/h1antiGoodHit_tpc", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});
        debugHistos.add<TH1>("tracks/trackingEff/h1antiGoodHit_its_tpc", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});
        debugHistos.add<TH1>("tracks/trackingEff/h1antiGoodHit_its_tpc_tof", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});
        debugHistos.add<TH1>("tracks/trackingEff/h1antiGoodHit_its_tpc_trd", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});
        debugHistos.add<TH1>("tracks/trackingEff/h1antiGoodHit_its_tpc_trd_tof", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});
      }

      // 1D pT
      if (enablePr) {
        histos.add<TH1>("tracks/proton/h1ProtonSpectraTrue", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
        histos.add<TH1>("tracks/proton/h1ProtonSpectraTrueWPID", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
        histos.add<TH1>("tracks/proton/h1ProtonSpectraTruePrim", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
        histos.add<TH1>("tracks/proton/h1ProtonSpectraTrueSec", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
        histos.add<TH1>("tracks/proton/h1ProtonSpectraTrueTransport", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});

        histos.add<TH1>("tracks/proton/h1antiProtonSpectraTrue", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});
        histos.add<TH1>("tracks/proton/h1antiProtonSpectraTrueWPID", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});
        histos.add<TH1>("tracks/proton/h1antiProtonSpectraTruePrim", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});
        histos.add<TH1>("tracks/proton/h1antiProtonSpectraTrueSec", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});
        histos.add<TH1>("tracks/proton/h1antiProtonSpectraTrueTransport", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});

        if (enableTrackingEff) {
          debugHistos.add<TH1>("tracks/proton/trackingEffPID/h1ProtonSpectraPIDTruePrim_its", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/proton/trackingEffPID/h1ProtonSpectraPIDTruePrim_tpc", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/proton/trackingEffPID/h1ProtonSpectraPIDTruePrim_its_tpc", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/proton/trackingEffPID/h1ProtonSpectraPIDTruePrim_its_tpc_tof", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/proton/trackingEffPID/h1ProtonSpectraPIDTruePrim_its_tpc_trd", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/proton/trackingEffPID/h1ProtonSpectraPIDTruePrim_its_tpc_trd_tof", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});

          debugHistos.add<TH1>("tracks/proton/trackingEffPID/h1ProtonSpectraPIDTruePrimGoodHit_its", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/proton/trackingEffPID/h1ProtonSpectraPIDTruePrimGoodHit_tpc", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/proton/trackingEffPID/h1ProtonSpectraPIDTruePrimGoodHit_its_tpc", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/proton/trackingEffPID/h1ProtonSpectraPIDTruePrimGoodHit_its_tpc_tof", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/proton/trackingEffPID/h1ProtonSpectraPIDTruePrimGoodHit_its_tpc_trd", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/proton/trackingEffPID/h1ProtonSpectraPIDTruePrimGoodHit_its_tpc_trd_tof", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});

          debugHistos.add<TH1>("tracks/proton/trackingEffPID/h1antiProtonSpectraPIDTruePrim_its", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/proton/trackingEffPID/h1antiProtonSpectraPIDTruePrim_tpc", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/proton/trackingEffPID/h1antiProtonSpectraPIDTruePrim_its_tpc", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/proton/trackingEffPID/h1antiProtonSpectraPIDTruePrim_its_tpc_tof", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/proton/trackingEffPID/h1antiProtonSpectraPIDTruePrim_its_tpc_trd", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/proton/trackingEffPID/h1antiProtonSpectraPIDTruePrim_its_tpc_trd_tof", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});

          debugHistos.add<TH1>("tracks/proton/trackingEffPID/h1antiProtonSpectraPIDTruePrimGoodHit_its", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/proton/trackingEffPID/h1antiProtonSpectraPIDTruePrimGoodHit_tpc", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/proton/trackingEffPID/h1antiProtonSpectraPIDTruePrimGoodHit_its_tpc", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/proton/trackingEffPID/h1antiProtonSpectraPIDTruePrimGoodHit_its_tpc_tof", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/proton/trackingEffPID/h1antiProtonSpectraPIDTruePrimGoodHit_its_tpc_trd", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/proton/trackingEffPID/h1antiProtonSpectraPIDTruePrimGoodHit_its_tpc_trd_tof", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});

          debugHistos.add<TH1>("tracks/proton/trackingEff/h1ProtonSpectra_its", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/proton/trackingEff/h1ProtonSpectra_tpc", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/proton/trackingEff/h1ProtonSpectra_its_tpc", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/proton/trackingEff/h1ProtonSpectra_its_tpc_tof", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/proton/trackingEff/h1ProtonSpectra_its_tpc_trd", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/proton/trackingEff/h1ProtonSpectra_its_tpc_trd_tof", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});

          debugHistos.add<TH1>("tracks/proton/trackingEff/h1ProtonSpectraGoodHit_its", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/proton/trackingEff/h1ProtonSpectraGoodHit_tpc", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/proton/trackingEff/h1ProtonSpectraGoodHit_its_tpc", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/proton/trackingEff/h1ProtonSpectraGoodHit_its_tpc_tof", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/proton/trackingEff/h1ProtonSpectraGoodHit_its_tpc_trd", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/proton/trackingEff/h1ProtonSpectraGoodHit_its_tpc_trd_tof", "#it{p}_{T} (p)", HistType::kTH1F, {ptAxis});

          debugHistos.add<TH1>("tracks/proton/trackingEff/h1antiProtonSpectra_its", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/proton/trackingEff/h1antiProtonSpectra_tpc", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/proton/trackingEff/h1antiProtonSpectra_its_tpc", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/proton/trackingEff/h1antiProtonSpectra_its_tpc_tof", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/proton/trackingEff/h1antiProtonSpectra_its_tpc_trd", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/proton/trackingEff/h1antiProtonSpectra_its_tpc_trd_tof", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});

          debugHistos.add<TH1>("tracks/proton/trackingEff/h1antiProtonSpectraGoodHit_its", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/proton/trackingEff/h1antiProtonSpectraGoodHit_tpc", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/proton/trackingEff/h1antiProtonSpectraGoodHit_its_tpc", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/proton/trackingEff/h1antiProtonSpectraGoodHit_its_tpc_tof", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/proton/trackingEff/h1antiProtonSpectraGoodHit_its_tpc_trd", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/proton/trackingEff/h1antiProtonSpectraGoodHit_its_tpc_trd_tof", "#it{p}_{T} (#bar{p})", HistType::kTH1F, {ptAxis});
        }
      }
      if (enableDe) {
        histos.add<TH1>("tracks/deuteron/h1DeuteronSpectraTrue", "#it{p}_{T} (d)", HistType::kTH1F, {ptAxis});
        histos.add<TH1>("tracks/deuteron/h1DeuteronSpectraTrueWPID", "#it{p}_{T} (d)", HistType::kTH1F, {ptAxis});
        histos.add<TH1>("tracks/deuteron/h1DeuteronSpectraTrueWPIDPrim", "#it{p}_{T} (d)", HistType::kTH1F, {ptAxis});
        histos.add<TH1>("tracks/deuteron/hPtDeuteronTOFTrue", "#it{p}_{T} (d) with TOF", HistType::kTH1F, {ptAxis});
        histos.add<TH1>("tracks/deuteron/hPtDeuteronTOFTrueWPIDPrim", "#it{p}_{T} (d) with TOF", HistType::kTH1F, {ptAxis});
        histos.add<TH1>("tracks/deuteron/h1DeuteronSpectraTruePrim", "#it{p}_{T} (d)", HistType::kTH1F, {ptAxis});
        histos.add<TH1>("tracks/deuteron/h1DeuteronSpectraTrueSec", "#it{p}_{T} (d)", HistType::kTH1F, {ptAxis});
        histos.add<TH1>("tracks/deuteron/h1DeuteronSpectraTrueTransport", "#it{p}_{T} (d)", HistType::kTH1F, {ptAxis});

        histos.add<TH1>("tracks/deuteron/h1antiDeuteronSpectraTrue", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
        histos.add<TH1>("tracks/deuteron/h1antiDeuteronSpectraTrueWPID", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
        histos.add<TH1>("tracks/deuteron/h1antiDeuteronSpectraTrueWPIDPrim", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
        histos.add<TH1>("tracks/deuteron/hPtantiDeuteronTOFTrue", "#it{p}_{T} (#bar{d}) with TOF", HistType::kTH1F, {ptAxis});
        histos.add<TH1>("tracks/deuteron/hPtantiDeuteronTOFTrueWPIDPrim", "#it{p}_{T} (#bar{d}) with TOF", HistType::kTH1F, {ptAxis});
        histos.add<TH1>("tracks/deuteron/h1antiDeuteronSpectraTruePrim", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
        histos.add<TH1>("tracks/deuteron/h1antiDeuteronSpectraTrueSec", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
        histos.add<TH1>("tracks/deuteron/h1antiDeuteronSpectraTrueTransport", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});

        if (enableTrackingEff) {
          debugHistos.add<TH1>("tracks/deuteron/trackingEffPID/h1DeuteronSpectraPIDTruePrim_its", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/deuteron/trackingEffPID/h1DeuteronSpectraPIDTruePrim_tpc", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/deuteron/trackingEffPID/h1DeuteronSpectraPIDTruePrim_its_tpc", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/deuteron/trackingEffPID/h1DeuteronSpectraPIDTruePrim_its_tpc_tof", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/deuteron/trackingEffPID/h1DeuteronSpectraPIDTruePrim_its_tpc_trd", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/deuteron/trackingEffPID/h1DeuteronSpectraPIDTruePrim_its_tpc_trd_tof", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});

          debugHistos.add<TH1>("tracks/deuteron/trackingEffPID/h1DeuteronSpectraPIDTruePrimGoodHit_its", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/deuteron/trackingEffPID/h1DeuteronSpectraPIDTruePrimGoodHit_tpc", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/deuteron/trackingEffPID/h1DeuteronSpectraPIDTruePrimGoodHit_its_tpc", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/deuteron/trackingEffPID/h1DeuteronSpectraPIDTruePrimGoodHit_its_tpc_tof", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/deuteron/trackingEffPID/h1DeuteronSpectraPIDTruePrimGoodHit_its_tpc_trd", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/deuteron/trackingEffPID/h1DeuteronSpectraPIDTruePrimGoodHit_its_tpc_trd_tof", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});

          debugHistos.add<TH1>("tracks/deuteron/trackingEffPID/h1antiDeuteronSpectraPIDTruePrim_its", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/deuteron/trackingEffPID/h1antiDeuteronSpectraPIDTruePrim_tpc", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/deuteron/trackingEffPID/h1antiDeuteronSpectraPIDTruePrim_its_tpc", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/deuteron/trackingEffPID/h1antiDeuteronSpectraPIDTruePrim_its_tpc_tof", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/deuteron/trackingEffPID/h1antiDeuteronSpectraPIDTruePrim_its_tpc_trd", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/deuteron/trackingEffPID/h1antiDeuteronSpectraPIDTruePrim_its_tpc_trd_tof", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});

          debugHistos.add<TH1>("tracks/deuteron/trackingEffPID/h1antiDeuteronSpectraPIDTruePrimGoodHit_its", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/deuteron/trackingEffPID/h1antiDeuteronSpectraPIDTruePrimGoodHit_tpc", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/deuteron/trackingEffPID/h1antiDeuteronSpectraPIDTruePrimGoodHit_its_tpc", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/deuteron/trackingEffPID/h1antiDeuteronSpectraPIDTruePrimGoodHit_its_tpc_tof", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/deuteron/trackingEffPID/h1antiDeuteronSpectraPIDTruePrimGoodHit_its_tpc_trd", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/deuteron/trackingEffPID/h1antiDeuteronSpectraPIDTruePrimGoodHit_its_tpc_trd_tof", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});

          debugHistos.add<TH1>("tracks/deuteron/trackingEff/h1DeuteronSpectra_its", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/deuteron/trackingEff/h1DeuteronSpectra_tpc", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/deuteron/trackingEff/h1DeuteronSpectra_its_tpc", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/deuteron/trackingEff/h1DeuteronSpectra_its_tpc_tof", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/deuteron/trackingEff/h1DeuteronSpectra_its_tpc_trd", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/deuteron/trackingEff/h1DeuteronSpectra_its_tpc_trd_tof", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});

          debugHistos.add<TH1>("tracks/deuteron/trackingEff/h1DeuteronSpectraGoodHit_its", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/deuteron/trackingEff/h1DeuteronSpectraGoodHit_tpc", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/deuteron/trackingEff/h1DeuteronSpectraGoodHit_its_tpc", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/deuteron/trackingEff/h1DeuteronSpectraGoodHit_its_tpc_tof", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/deuteron/trackingEff/h1DeuteronSpectraGoodHit_its_tpc_trd", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/deuteron/trackingEff/h1DeuteronSpectraGoodHit_its_tpc_trd_tof", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});

          debugHistos.add<TH1>("tracks/deuteron/trackingEff/h1antiDeuteronSpectra_its", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/deuteron/trackingEff/h1antiDeuteronSpectra_tpc", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/deuteron/trackingEff/h1antiDeuteronSpectra_its_tpc", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/deuteron/trackingEff/h1antiDeuteronSpectra_its_tpc_tof", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/deuteron/trackingEff/h1antiDeuteronSpectra_its_tpc_trd", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/deuteron/trackingEff/h1antiDeuteronSpectra_its_tpc_trd_tof", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});

          debugHistos.add<TH1>("tracks/deuteron/trackingEff/h1antiDeuteronSpectraGoodHit_its", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/deuteron/trackingEff/h1antiDeuteronSpectraGoodHit_tpc", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/deuteron/trackingEff/h1antiDeuteronSpectraGoodHit_its_tpc", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/deuteron/trackingEff/h1antiDeuteronSpectraGoodHit_its_tpc_tof", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/deuteron/trackingEff/h1antiDeuteronSpectraGoodHit_its_tpc_trd", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
          debugHistos.add<TH1>("tracks/deuteron/trackingEff/h1antiDeuteronSpectraGoodHit_its_tpc_trd_tof", "#it{p}_{T} (#bar{d})", HistType::kTH1F, {ptAxis});
        }
      }
      if (enableTr) {
        histos.add<TH1>("tracks/triton/h1TritonSpectraTrue", "#it{p}_{T} (t)", HistType::kTH1F, {ptAxis});
        histos.add<TH1>("tracks/triton/h1TritonSpectraTruePrim", "#it{p}_{T} (t)", HistType::kTH1F, {ptAxis});
        histos.add<TH1>("tracks/triton/h1TritonSpectraTrueSec", "#it{p}_{T} (t)", HistType::kTH1F, {ptAxis});
        histos.add<TH1>("tracks/triton/h1TritonSpectraTrueTransport", "#it{p}_{T} (t)", HistType::kTH1F, {ptAxis});

        histos.add<TH1>("tracks/triton/h1antiTritonSpectraTrue", "#it{p}_{T} (#bar{t})", HistType::kTH1F, {ptAxis});
        histos.add<TH1>("tracks/triton/h1antiTritonSpectraTruePrim", "#it{p}_{T} (#bar{t})", HistType::kTH1F, {ptAxis});
        histos.add<TH1>("tracks/triton/h1antiTritonSpectraTrueSec", "#it{p}_{T} (#bar{t})", HistType::kTH1F, {ptAxis});
        histos.add<TH1>("tracks/triton/h1antiTritonSpectraTrueTransport", "#it{p}_{T} (#bar{t})", HistType::kTH1F, {ptAxis});
      }
      if (enableHe) {
        histos.add<TH1>("tracks/helium/h1HeliumSpectraTrue_Z2", "#it{p}_{T} (He)", HistType::kTH1F, {ptHeAxis});
        histos.add<TH1>("tracks/helium/h1HeliumSpectraTrueWPID_Z2", "#it{p}_{T} (He)", HistType::kTH1F, {ptHeAxis});
        histos.add<TH1>("tracks/helium/h1HeliumSpectraTruePrim_Z2", "#it{p}_{T} (He)", HistType::kTH1F, {ptHeAxis});
        histos.add<TH1>("tracks/helium/h1HeliumSpectraTrueSec_Z2", "#it{p}_{T} (He)", HistType::kTH1F, {ptHeAxis});
        histos.add<TH1>("tracks/helium/h1HeliumSpectraTrueTransport_Z2", "#it{p}_{T} (He)", HistType::kTH1F, {ptHeAxis});

        if (enableCentrality) {
          histos.add<TH2>("tracks/helium/h2HeliumSpectraTrueVsMult_Z2", "#it{p}_{T} (He)", HistType::kTH2F, {{ptHeAxis}, {binsPercentile}});
          histos.add<TH2>("tracks/helium/h2HeliumSpectraTrueWPIDVsMult_Z2", "#it{p}_{T} (He)", HistType::kTH2F, {{ptHeAxis}, {binsPercentile}});
          histos.add<TH2>("tracks/helium/h2HeliumSpectraTruePrimVsMult_Z2", "#it{p}_{T} (He)", HistType::kTH2F, {{ptHeAxis}, {binsPercentile}});
          histos.add<TH2>("tracks/helium/h2HeliumSpectraTrueSecVsMult_Z2", "#it{p}_{T} (He)", HistType::kTH2F, {{ptHeAxis}, {binsPercentile}});
        }

        histos.add<TH1>("tracks/helium/h1antiHeliumSpectraTrue_Z2", "#it{p}_{T} (He)", HistType::kTH1F, {ptHeAxis});
        histos.add<TH1>("tracks/helium/h1antiHeliumSpectraTrueWPID_Z2", "#it{p}_{T} (He)", HistType::kTH1F, {ptHeAxis});
        histos.add<TH1>("tracks/helium/h1antiHeliumSpectraTruePrim_Z2", "#it{p}_{T} (He)", HistType::kTH1F, {ptHeAxis});
        histos.add<TH1>("tracks/helium/h1antiHeliumSpectraTrueSec_Z2", "#it{p}_{T} (He)", HistType::kTH1F, {ptHeAxis});
        histos.add<TH1>("tracks/helium/h1antiHeliumSpectraTrueTransport_Z2", "#it{p}_{T} (He)", HistType::kTH1F, {ptHeAxis});

        if (enableCentrality) {
          histos.add<TH2>("tracks/helium/h2antiHeliumSpectraTrueVsMult_Z2", "#it{p}_{T} (He)", HistType::kTH2F, {{ptHeAxis}, {binsPercentile}});
          histos.add<TH2>("tracks/helium/h2antiHeliumSpectraTrueWPIDVsMult_Z2", "#it{p}_{T} (He)", HistType::kTH2F, {{ptHeAxis}, {binsPercentile}});
          histos.add<TH2>("tracks/helium/h2antiHeliumSpectraTruePrimVsMult_Z2", "#it{p}_{T} (He)", HistType::kTH2F, {{ptHeAxis}, {binsPercentile}});
          histos.add<TH2>("tracks/helium/h2antiHeliumSpectraTrueSecVsMult_Z2", "#it{p}_{T} (He)", HistType::kTH2F, {{ptHeAxis}, {binsPercentile}});
        }

        if (outFlagOptions.doTOFplots) {
          histos.add<TH1>("tracks/helium/TOF/h1HeliumSpectraTruePrim_Z2", "#it{p}_{T} (He)", HistType::kTH1F, {ptHeAxis});
          histos.add<TH1>("tracks/helium/TOF/h1antiHeliumSpectraTruePrim_Z2", "#it{p}_{T} (He)", HistType::kTH1F, {ptHeAxis});
          if (enableCentrality) {
            histos.add<TH2>("tracks/helium/TOF/h2HeliumSpectraTruePrimVsMult_Z2", "#it{p}_{T} (He)", HistType::kTH2F, {{ptHeAxis}, {binsPercentile}});
            histos.add<TH2>("tracks/helium/TOF/h2antiHeliumSpectraTruePrimVsMult_Z2", "#it{p}_{T} (He)", HistType::kTH2F, {{ptHeAxis}, {binsPercentile}});
          }
        }
        if (outFlagOptions.enableEffPlots) {
          histos.add<TH1>("tracks/eff/helium/hPtHeTrue_Z2", "Track #it{p}_{T} (He); #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{400, 0., 8.}});
          histos.add<TH1>("tracks/eff/helium/hPtantiHeTrue_Z2", "Track #it{p}_{T} (#bar{He}); #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{400, 0., 8.}});
          if (outFlagOptions.doTOFplots) {
            histos.add<TH1>("tracks/eff/helium/hPtHeTOFTrue_Z2", "Track #it{p}_{T} (He); #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{400, 0., 8.}});
            histos.add<TH1>("tracks/eff/helium/hPtantiHeTOFTrue_Z2", "Track #it{p}_{T} (#bar{He}); #it{p}_{T} (GeV/#it{c}); counts", HistType::kTH1F, {{400, 0., 8.}});
          }
        }
      }
      if (enableAl) {
        histos.add<TH1>("tracks/alpha/h1AlphaSpectraTrue", "#it{p}_{T} (#alpha)", HistType::kTH1F, {ptAxis});
        histos.add<TH1>("tracks/alpha/h1AlphaSpectraTruePrim", "#it{p}_{T} (#alpha)", HistType::kTH1F, {ptAxis});
        histos.add<TH1>("tracks/alpha/h1AlphaSpectraTrueSec", "#it{p}_{T} (#alpha)", HistType::kTH1F, {ptAxis});
        histos.add<TH1>("tracks/alpha/h1AlphaSpectraTrueTransport", "#it{p}_{T} (#alpha)", HistType::kTH1F, {ptAxis});

        histos.add<TH1>("tracks/alpha/h1antiAlphaSpectraTrue", "#it{p}_{T} (#bar{#alpha})", HistType::kTH1F, {ptAxis});
        histos.add<TH1>("tracks/alpha/h1antiAlphaSpectraTruePrim", "#it{p}_{T} (#bar{#alpha})", HistType::kTH1F, {ptAxis});
        histos.add<TH1>("tracks/alpha/h1antiAlphaSpectraTrueSec", "#it{p}_{T} (#bar{#alpha})", HistType::kTH1F, {ptAxis});
        histos.add<TH1>("tracks/alpha/h1antiAlphaSpectraTrueTransport", "#it{p}_{T} (#bar{#alpha})", HistType::kTH1F, {ptAxis});
      }
      // 2D-DCAxy
      if (enablePr) {
        if (outFlagOptions.makeDCABeforeCutPlots) {
          histos.add<TH2>("tracks/proton/dca/before/hDCAxyVsPtProtonTrue", "DCAxy vs Pt (p); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

          histos.add<TH2>("tracks/proton/dca/before/hDCAxyVsPtProtonTruePrim", "DCAxy vs Pt (p); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/proton/dca/before/hDCAxyVsPtProtonTrueSec", "DCAxy vs Pt (p); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/proton/dca/before/hDCAxyVsPtProtonTrueMaterial", "DCAxy vs Pt (p); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

          histos.add<TH1>("tracks/proton/dca/before/hNumMothers", "N mothers per particle; N mothers;counts", HistType::kTH1I, {{7, 1.0, 8.0}});
          histos.add<TH3>("tracks/proton/dca/before/hMomTrueMaterial", "MC mothers;mother index;mother type; mother #it{p}_{T}", HistType::kTH3F, {{2, -2.0, 2.0}, {kNumMotherList + 2, -1.5, static_cast<double>(kNumMotherList) + 0.5}, {150, 0.0, 15.0}});

          std::shared_ptr<TH3> hTempPr = histos.get<TH3>(HIST("tracks/proton/dca/before/hMomTrueMaterial"));
          TH3* hPdgPr = hTempPr.get();

          TAxis* axPdgPr = hPdgPr->GetXaxis();
          axPdgPr->SetBinLabel(1, "antiparticles");
          axPdgPr->SetBinLabel(2, "particles");

          TAxis* ayPdgPr = hPdgPr->GetYaxis();
          ayPdgPr->SetBinLabel(1, "undef.");
          ayPdgPr->SetBinLabel(2, "other");
          for (int i = 0; i < kNumMotherList; i++) {
            ayPdgPr->SetBinLabel(i + 3, kMotherNames[i]);
          }

          histos.add<TH2>("tracks/proton/dca/before/hDCAxyVsPtantiProtonTrue", "DCAxy vs Pt (#bar{p}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

          histos.add<TH2>("tracks/proton/dca/before/hDCAxyVsPtantiProtonTruePrim", "DCAxy vs Pt (#bar{p}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/proton/dca/before/hDCAxyVsPtantiProtonTrueSec", "DCAxy vs Pt (#bar{p}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/proton/dca/before/hDCAxyVsPtantiProtonTrueMaterial", "DCAxy vs Pt (#bar{p}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

          if (outFlagOptions.doTOFplots) {
            histos.add<TH2>("tracks/proton/dca/before/TOF/hDCAxyVsPtProtonTrue", "DCAxy vs Pt (p); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

            histos.add<TH2>("tracks/proton/dca/before/TOF/hDCAxyVsPtProtonTruePrim", "DCAxy vs Pt (p); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
            histos.add<TH2>("tracks/proton/dca/before/TOF/hDCAxyVsPtProtonTrueSec", "DCAxy vs Pt (p); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
            histos.add<TH2>("tracks/proton/dca/before/TOF/hDCAxyVsPtProtonTrueMaterial", "DCAxy vs Pt (p); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

            histos.add<TH2>("tracks/proton/dca/before/TOF/hDCAxyVsPtantiProtonTrue", "DCAxy vs Pt (#bar{p}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

            histos.add<TH2>("tracks/proton/dca/before/TOF/hDCAxyVsPtantiProtonTruePrim", "DCAxy vs Pt (#bar{p}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
            histos.add<TH2>("tracks/proton/dca/before/TOF/hDCAxyVsPtantiProtonTrueSec", "DCAxy vs Pt (#bar{p}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
            histos.add<TH2>("tracks/proton/dca/before/TOF/hDCAxyVsPtantiProtonTrueMaterial", "DCAxy vs Pt (#bar{p}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

            histos.add<TH2>("tracks/proton/dca/before/TOF/hDCAzVsPtProtonTrue", "DCAz vs Pt (p); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

            histos.add<TH2>("tracks/proton/dca/before/TOF/hDCAzVsPtProtonTruePrim", "DCAz vs Pt (p); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
            histos.add<TH2>("tracks/proton/dca/before/TOF/hDCAzVsPtProtonTrueSec", "DCAz vs Pt (p); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
            histos.add<TH2>("tracks/proton/dca/before/TOF/hDCAzVsPtProtonTrueMaterial", "DCAz vs Pt (p); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

            histos.add<TH2>("tracks/proton/dca/before/TOF/hDCAzVsPtantiProtonTrue", "DCAz vs Pt (#bar{p}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

            histos.add<TH2>("tracks/proton/dca/before/TOF/hDCAzVsPtantiProtonTruePrim", "DCAz vs Pt (#bar{p}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
            histos.add<TH2>("tracks/proton/dca/before/TOF/hDCAzVsPtantiProtonTrueSec", "DCAz vs Pt (#bar{p}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
            histos.add<TH2>("tracks/proton/dca/before/TOF/hDCAzVsPtantiProtonTrueMaterial", "DCAz vs Pt (#bar{p}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          }
        }

        if (outFlagOptions.makeDCAAfterCutPlots) {
          histos.add<TH2>("tracks/proton/dca/after/hDCAxyVsPtProtonTrue", "DCAxy vs Pt (p); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

          histos.add<TH2>("tracks/proton/dca/after/hDCAxyVsPtProtonTruePrim", "DCAxy vs Pt (p); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/proton/dca/after/hDCAxyVsPtProtonTrueSec", "DCAxy vs Pt (p); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/proton/dca/after/hDCAxyVsPtProtonTrueTransport", "DCAxy vs Pt (p); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

          histos.add<TH2>("tracks/proton/dca/after/hDCAxyVsPtantiProtonTrue", "DCAxy vs Pt (#bar{p}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

          histos.add<TH2>("tracks/proton/dca/after/hDCAxyVsPtantiProtonTruePrim", "DCAxy vs Pt (#bar{p}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/proton/dca/after/hDCAxyVsPtantiProtonTrueSec", "DCAxy vs Pt (#bar{p}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/proton/dca/after/hDCAxyVsPtantiProtonTrueTransport", "DCAxy vs Pt (#bar{p}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
        }
      }
      if (enableDe) {
        if (outFlagOptions.makeDCABeforeCutPlots) {
          histos.add<TH2>("tracks/deuteron/dca/before/hDCAxyVsPtDeuteronTrue", "DCAxy vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

          histos.add<TH2>("tracks/deuteron/dca/before/hDCAxyVsPtDeuteronTruePrim", "DCAxy vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/deuteron/dca/before/hDCAxyVsPtDeuteronTrueSec", "DCAxy vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/deuteron/dca/before/hDCAxyVsPtDeuteronTrueMaterial", "DCAxy vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

          histos.add<TH1>("tracks/deuteron/dca/before/hNumMothers", "N mothers per particle; N mothers;counts", HistType::kTH1I, {{7, 1.0, 8.0}});
          histos.add<TH3>("tracks/deuteron/dca/before/hMomTrueMaterial", "MC mothers;mother index;mother type; mother #it{p}_{T}", HistType::kTH3F, {{2, -2.0, 2.0}, {kNumMotherList + 2, -1.5, static_cast<double>(kNumMotherList) + 0.5}, {150, 0.0, 15.0}});

          std::shared_ptr<TH3> hTempDe = histos.get<TH3>(HIST("tracks/deuteron/dca/before/hMomTrueMaterial"));
          TH3* hPdgDe = hTempDe.get();

          TAxis* axPdgDe = hPdgDe->GetXaxis();
          axPdgDe->SetBinLabel(1, "antiparticles");
          axPdgDe->SetBinLabel(2, "particles");

          TAxis* ayPdgDe = hPdgDe->GetYaxis();
          ayPdgDe->SetBinLabel(1, "undef.");
          ayPdgDe->SetBinLabel(2, "other");
          for (int i = 0; i < kNumMotherList; i++) {
            ayPdgDe->SetBinLabel(i + 3, kMotherNames[i]);
          }

          histos.add<TH2>("tracks/deuteron/dca/before/hDCAxyVsPtantiDeuteronTrue", "DCAxy vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

          histos.add<TH2>("tracks/deuteron/dca/before/hDCAxyVsPtantiDeuteronTruePrim", "DCAxy vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/deuteron/dca/before/hDCAxyVsPtantiDeuteronTrueSec", "DCAxy vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/deuteron/dca/before/hDCAxyVsPtantiDeuteronTrueMaterial", "DCAxy vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

          // Fake & wrong histos
          if (outFlagOptions.makeFakeTracksPlots) {
            histos.add<TH2>("tracks/deuteron/dca/before/fake/hDCAxyVsPtDeuteronTrue", "DCAxy vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

            histos.add<TH2>("tracks/deuteron/dca/before/fake/hDCAxyVsPtDeuteronTruePrim", "DCAxy vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
            histos.add<TH2>("tracks/deuteron/dca/before/fake/hDCAxyVsPtDeuteronTrueSec", "DCAxy vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
            histos.add<TH2>("tracks/deuteron/dca/before/fake/hDCAxyVsPtDeuteronTrueMaterial", "DCAxy vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

            histos.add<TH2>("tracks/deuteron/dca/before/fake/hDCAxyVsPtantiDeuteronTrue", "DCAxy vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

            histos.add<TH2>("tracks/deuteron/dca/before/fake/hDCAxyVsPtantiDeuteronTruePrim", "DCAxy vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
            histos.add<TH2>("tracks/deuteron/dca/before/fake/hDCAxyVsPtantiDeuteronTrueSec", "DCAxy vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
            histos.add<TH2>("tracks/deuteron/dca/before/fake/hDCAxyVsPtantiDeuteronTrueMaterial", "DCAxy vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          }

          if (outFlagOptions.doTOFplots) {
            histos.add<TH2>("tracks/deuteron/dca/before/TOF/hDCAxyVsPtDeuteronTrue", "DCAxy vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

            histos.add<TH2>("tracks/deuteron/dca/before/TOF/hDCAxyVsPtDeuteronTruePrim", "DCAxy vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
            histos.add<TH2>("tracks/deuteron/dca/before/TOF/hDCAxyVsPtDeuteronTrueSec", "DCAxy vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
            histos.add<TH2>("tracks/deuteron/dca/before/TOF/hDCAxyVsPtDeuteronTrueMaterial", "DCAxy vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

            histos.add<TH2>("tracks/deuteron/dca/before/TOF/hDCAxyVsPtantiDeuteronTrue", "DCAxy vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

            histos.add<TH2>("tracks/deuteron/dca/before/TOF/hDCAxyVsPtantiDeuteronTruePrim", "DCAxy vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
            histos.add<TH2>("tracks/deuteron/dca/before/TOF/hDCAxyVsPtantiDeuteronTrueSec", "DCAxy vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
            histos.add<TH2>("tracks/deuteron/dca/before/TOF/hDCAxyVsPtantiDeuteronTrueMaterial", "DCAxy vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

            histos.add<TH2>("tracks/deuteron/dca/before/TOF/hDCAzVsPtDeuteronTrue", "DCAz vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

            histos.add<TH2>("tracks/deuteron/dca/before/TOF/hDCAzVsPtDeuteronTruePrim", "DCAz vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
            histos.add<TH2>("tracks/deuteron/dca/before/TOF/hDCAzVsPtDeuteronTrueSec", "DCAz vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
            histos.add<TH2>("tracks/deuteron/dca/before/TOF/hDCAzVsPtDeuteronTrueMaterial", "DCAz vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

            histos.add<TH2>("tracks/deuteron/dca/before/TOF/hDCAzVsPtantiDeuteronTrue", "DCAz vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

            histos.add<TH2>("tracks/deuteron/dca/before/TOF/hDCAzVsPtantiDeuteronTruePrim", "DCAz vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
            histos.add<TH2>("tracks/deuteron/dca/before/TOF/hDCAzVsPtantiDeuteronTrueSec", "DCAz vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
            histos.add<TH2>("tracks/deuteron/dca/before/TOF/hDCAzVsPtantiDeuteronTrueMaterial", "DCAz vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

            if (outFlagOptions.makeFakeTracksPlots) {
              histos.add<TH2>("tracks/deuteron/dca/before/fake/TOF/hDCAxyVsPtDeuteronTrue", "DCAxy vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

              histos.add<TH2>("tracks/deuteron/dca/before/fake/TOF/hDCAxyVsPtDeuteronTruePrim", "DCAxy vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
              histos.add<TH2>("tracks/deuteron/dca/before/fake/TOF/hDCAxyVsPtDeuteronTrueSec", "DCAxy vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
              histos.add<TH2>("tracks/deuteron/dca/before/fake/TOF/hDCAxyVsPtDeuteronTrueMaterial", "DCAxy vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

              histos.add<TH2>("tracks/deuteron/dca/before/fake/TOF/hDCAxyVsPtantiDeuteronTrue", "DCAxy vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

              histos.add<TH2>("tracks/deuteron/dca/before/fake/TOF/hDCAxyVsPtantiDeuteronTruePrim", "DCAxy vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
              histos.add<TH2>("tracks/deuteron/dca/before/fake/TOF/hDCAxyVsPtantiDeuteronTrueSec", "DCAxy vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
              histos.add<TH2>("tracks/deuteron/dca/before/fake/TOF/hDCAxyVsPtantiDeuteronTrueMaterial", "DCAxy vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

              histos.add<TH2>("tracks/deuteron/dca/before/fake/TOF/hDCAzVsPtDeuteronTrue", "DCAz vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

              histos.add<TH2>("tracks/deuteron/dca/before/fake/TOF/hDCAzVsPtDeuteronTruePrim", "DCAz vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
              histos.add<TH2>("tracks/deuteron/dca/before/fake/TOF/hDCAzVsPtDeuteronTrueSec", "DCAz vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
              histos.add<TH2>("tracks/deuteron/dca/before/fake/TOF/hDCAzVsPtDeuteronTrueTransport", "DCAz vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

              histos.add<TH2>("tracks/deuteron/dca/before/fake/TOF/hDCAzVsPtantiDeuteronTrue", "DCAz vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

              histos.add<TH2>("tracks/deuteron/dca/before/fake/TOF/hDCAzVsPtantiDeuteronTruePrim", "DCAz vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
              histos.add<TH2>("tracks/deuteron/dca/before/fake/TOF/hDCAzVsPtantiDeuteronTrueSec", "DCAz vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
              histos.add<TH2>("tracks/deuteron/dca/before/fake/TOF/hDCAzVsPtantiDeuteronTrueTransport", "DCAz vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
            }
          }
        }

        if (outFlagOptions.makeDCAAfterCutPlots) {
          histos.add<TH2>("tracks/deuteron/dca/after/hDCAxyVsPtDeuteronTrue", "DCAxy vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

          histos.add<TH2>("tracks/deuteron/dca/after/hDCAxyVsPtDeuteronTruePrim", "DCAxy vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/deuteron/dca/after/hDCAxyVsPtDeuteronTrueWPIDPrim", "DCAxy vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/deuteron/dca/after/hDCAxyVsPtDeuteronTrueSec", "DCAxy vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/deuteron/dca/after/hDCAxyVsPtDeuteronTrueTransport", "DCAxy vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

          histos.add<TH2>("tracks/deuteron/dca/after/hDCAxyVsPtantiDeuteronTrue", "DCAxy vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

          histos.add<TH2>("tracks/deuteron/dca/after/hDCAxyVsPtantiDeuteronTruePrim", "DCAxy vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/deuteron/dca/after/hDCAxyVsPtantiDeuteronTrueWPIDPrim", "DCAxy vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/deuteron/dca/after/hDCAxyVsPtantiDeuteronTrueSec", "DCAxy vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/deuteron/dca/after/hDCAxyVsPtantiDeuteronTrueTransport", "DCAxy vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
        }
      }
      if (enableTr) {
        if (outFlagOptions.makeDCABeforeCutPlots) {
          histos.add<TH2>("tracks/triton/dca/before/hDCAxyVsPtTritonTrue", "DCAxy vs Pt (t); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

          histos.add<TH2>("tracks/triton/dca/before/hDCAxyVsPtTritonTruePrim", "DCAxy vs Pt (t); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/triton/dca/before/hDCAxyVsPtTritonTrueSec", "DCAxy vs Pt (t); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/triton/dca/before/hDCAxyVsPtTritonTrueMaterial", "DCAxy vs Pt (t); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

          histos.add<TH2>("tracks/triton/dca/before/hDCAxyVsPtantiTritonTrue", "DCAxy vs Pt (#bar{t}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

          histos.add<TH2>("tracks/triton/dca/before/hDCAxyVsPtantiTritonTruePrim", "DCAxy vs Pt (#bar{t}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/triton/dca/before/hDCAxyVsPtantiTritonTrueSec", "DCAxy vs Pt (#bar{t}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/triton/dca/before/hDCAxyVsPtantiTritonTrueMaterial", "DCAxy vs Pt (t); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

          if (outFlagOptions.doTOFplots) {
            histos.add<TH2>("tracks/triton/dca/before/TOF/hDCAxyVsPtTritonTrue", "DCAxy vs Pt (t); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

            histos.add<TH2>("tracks/triton/dca/before/TOF/hDCAxyVsPtTritonTruePrim", "DCAxy vs Pt (t); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
            histos.add<TH2>("tracks/triton/dca/before/TOF/hDCAxyVsPtTritonTrueSec", "DCAxy vs Pt (t); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
            histos.add<TH2>("tracks/triton/dca/before/TOF/hDCAxyVsPtTritonTrueMaterial", "DCAxy vs Pt (t); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

            histos.add<TH2>("tracks/triton/dca/before/TOF/hDCAxyVsPtantiTritonTrue", "DCAxy vs Pt (#bar{t}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

            histos.add<TH2>("tracks/triton/dca/before/TOF/hDCAxyVsPtantiTritonTruePrim", "DCAxy vs Pt (#bar{t}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
            histos.add<TH2>("tracks/triton/dca/before/TOF/hDCAxyVsPtantiTritonTrueSec", "DCAxy vs Pt (#bar{t}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
            histos.add<TH2>("tracks/triton/dca/before/TOF/hDCAxyVsPtantiTritonTrueMaterial", "DCAxy vs Pt (#bar{t}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

            // Unused histograms
            // histos.add<TH2>("tracks/triton/dca/before/TOF/hDCAzVsPtTritonTrue", "DCAz vs Pt (t); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

            // histos.add<TH2>("tracks/triton/dca/before/TOF/hDCAzVsPtTritonTruePrim", "DCAz vs Pt (t); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
            // histos.add<TH2>("tracks/triton/dca/before/TOF/hDCAzVsPtTritonTrueSec", "DCAz vs Pt (t); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
            // histos.add<TH2>("tracks/triton/dca/before/TOF/hDCAzVsPtTritonTrueMaterial", "DCAz vs Pt (t); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
            // histos.add<TH2>("tracks/triton/dca/before/TOF/hDCAzVsPtTritonTrueTransport", "DCAz vs Pt (t); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

            // histos.add<TH2>("tracks/triton/dca/before/TOF/hDCAzVsPtantiTritonTrue", "DCAz vs Pt (#bar{t}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

            // histos.add<TH2>("tracks/triton/dca/before/TOF/hDCAzVsPtantiTritonTruePrim", "DCAz vs Pt (#bar{t}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
            // histos.add<TH2>("tracks/triton/dca/before/TOF/hDCAzVsPtantiTritonTrueSec", "DCAz vs Pt (#bar{t}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
            // histos.add<TH2>("tracks/triton/dca/before/TOF/hDCAzVsPtantiTritonTrueMaterial", "DCAz vs Pt (#bar{t}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
            // histos.add<TH2>("tracks/triton/dca/before/TOF/hDCAzVsPtantiTritonTrueTransport", "DCAz vs Pt (#bar{t}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          }
        }

        if (outFlagOptions.makeDCAAfterCutPlots) {
          histos.add<TH2>("tracks/triton/dca/after/hDCAxyVsPtTritonTrue", "DCAxy vs Pt (t); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

          histos.add<TH2>("tracks/triton/dca/after/hDCAxyVsPtTritonTruePrim", "DCAxy vs Pt (t); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/triton/dca/after/hDCAxyVsPtTritonTrueSec", "DCAxy vs Pt (t); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/triton/dca/after/hDCAxyVsPtTritonTrueTransport", "DCAxy vs Pt (t); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

          histos.add<TH2>("tracks/triton/dca/after/hDCAxyVsPtantiTritonTrue", "DCAxy vs Pt (#bar{t}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

          histos.add<TH2>("tracks/triton/dca/after/hDCAxyVsPtantiTritonTruePrim", "DCAxy vs Pt (#bar{t}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/triton/dca/after/hDCAxyVsPtantiTritonTrueSec", "DCAxy vs Pt (#bar{t}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/triton/dca/after/hDCAxyVsPtantiTritonTrueTransport", "DCAxy vs Pt (#bar{t}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

          // Unused histograms
          // if (outFlagOptions.doTOFplots) {
          //   histos.add<TH2>("tracks/triton/dca/after/TOF/hDCAxyVsPtTritonTrue", "DCAxy vs Pt (t); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

          //   histos.add<TH2>("tracks/triton/dca/after/TOF/hDCAxyVsPtTritonTruePrim", "DCAxy vs Pt (t); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          //   histos.add<TH2>("tracks/triton/dca/after/TOF/hDCAxyVsPtTritonTrueSec", "DCAxy vs Pt (t); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          //   histos.add<TH2>("tracks/triton/dca/after/TOF/hDCAxyVsPtTritonTrueTransport", "DCAxy vs Pt (t); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

          //   histos.add<TH2>("tracks/triton/dca/after/TOF/hDCAxyVsPtantiTritonTrue", "DCAxy vs Pt (#bar{t}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

          //   histos.add<TH2>("tracks/triton/dca/after/TOF/hDCAxyVsPtantiTritonTruePrim", "DCAxy vs Pt (#bar{t}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          //   histos.add<TH2>("tracks/triton/dca/after/TOF/hDCAxyVsPtantiTritonTrueSec", "DCAxy vs Pt (#bar{t}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          //   histos.add<TH2>("tracks/triton/dca/after/TOF/hDCAxyVsPtantiTritonTrueTransport", "DCAxy vs Pt (#bar{t}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

          //   histos.add<TH2>("tracks/triton/dca/after/TOF/hDCAzVsPtTritonTrue", "DCAz vs Pt (t); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

          //   histos.add<TH2>("tracks/triton/dca/after/TOF/hDCAzVsPtTritonTruePrim", "DCAz vs Pt (t); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          //   histos.add<TH2>("tracks/triton/dca/after/TOF/hDCAzVsPtTritonTrueSec", "DCAz vs Pt (t); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          //   histos.add<TH2>("tracks/triton/dca/after/TOF/hDCAzVsPtTritonTrueTransport", "DCAz vs Pt (t); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

          //   histos.add<TH2>("tracks/triton/dca/after/TOF/hDCAzVsPtantiTritonTrue", "DCAz vs Pt (#bar{t}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

          //   histos.add<TH2>("tracks/triton/dca/after/TOF/hDCAzVsPtantiTritonTruePrim", "DCAz vs Pt (#bar{t}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          //   histos.add<TH2>("tracks/triton/dca/after/TOF/hDCAzVsPtantiTritonTrueSec", "DCAz vs Pt (#bar{t}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          //   histos.add<TH2>("tracks/triton/dca/after/TOF/hDCAzVsPtantiTritonTrueTransport", "DCAz vs Pt (#bar{t}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          // }
        }
      }
      if (enableHe) {
        // all tracks
        if (outFlagOptions.makeDCABeforeCutPlots) {
          histos.add<TH2>("tracks/helium/dca/before/hDCAxyVsPtHeliumTrue", "DCAxy vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});

          histos.add<TH2>("tracks/helium/dca/before/hDCAxyVsPtHeliumTruePrim", "DCAxy vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/helium/dca/before/hDCAxyVsPtHeliumTrueSec", "DCAxy vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/helium/dca/before/hDCAxyVsPtHeliumTrueMaterial", "DCAxy vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});

          histos.add<TH1>("tracks/helium/dca/before/hNumMothers", "N mothers per particle; N mothers;counts", HistType::kTH1I, {{7, 1.0, 8.0}});
          histos.add<TH3>("tracks/helium/dca/before/hMomTrueMaterial", "MC mothers;mother index;mother type; mother #it{p}_{T}", HistType::kTH3F, {{2, -2.0, 2.0}, {kNumMotherList + 2, -1.5, static_cast<double>(kNumMotherList) + 0.5}, {150, 0.0, 15.0}});

          // Fix for getting TH3 pointer
          std::shared_ptr<TH3> hTempHe = histos.get<TH3>(HIST("tracks/helium/dca/before/hMomTrueMaterial"));
          TH3* hPdgHe = hTempHe.get();

          TAxis* axPdgHe = hPdgHe->GetXaxis();
          axPdgHe->SetBinLabel(1, "antiparticles");
          axPdgHe->SetBinLabel(2, "particles");

          TAxis* ayPdgHe = hPdgHe->GetYaxis();
          ayPdgHe->SetBinLabel(1, "undef.");
          ayPdgHe->SetBinLabel(2, "other");
          for (int i = 0; i < kNumMotherList; i++) {
            ayPdgHe->SetBinLabel(i + 3, kMotherNames[i]);
          }

          histos.add<TH2>("tracks/helium/dca/before/hDCAxyVsPtantiHeliumTrue", "DCAxy vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});

          histos.add<TH2>("tracks/helium/dca/before/hDCAxyVsPtantiHeliumTruePrim", "DCAxy vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/helium/dca/before/hDCAxyVsPtantiHeliumTrueSec", "DCAxy vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/helium/dca/before/hDCAxyVsPtantiHeliumTrueMaterial", "DCAxy vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});

          histos.add<TH2>("tracks/helium/dca/before/hDCAzVsPtHeliumTrue", "DCAz vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});

          histos.add<TH2>("tracks/helium/dca/before/hDCAzVsPtHeliumTruePrim", "DCAz vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/helium/dca/before/hDCAzVsPtHeliumTrueSec", "DCAz vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/helium/dca/before/hDCAzVsPtHeliumTrueMaterial", "DCAz vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});

          histos.add<TH2>("tracks/helium/dca/before/hDCAzVsPtantiHeliumTrue", "DCAz vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});

          histos.add<TH2>("tracks/helium/dca/before/hDCAzVsPtantiHeliumTruePrim", "DCAz vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/helium/dca/before/hDCAzVsPtantiHeliumTrueSec", "DCAz vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/helium/dca/before/hDCAzVsPtantiHeliumTrueMaterial", "DCAz vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});

          if (outFlagOptions.doTOFplots) {
            histos.add<TH2>("tracks/helium/dca/before/TOF/hDCAxyVsPtHeliumTrue", "DCAxy vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});

            histos.add<TH2>("tracks/helium/dca/before/TOF/hDCAxyVsPtHeliumTruePrim", "DCAxy vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});
            histos.add<TH2>("tracks/helium/dca/before/TOF/hDCAxyVsPtHeliumTrueSec", "DCAxy vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});
            histos.add<TH2>("tracks/helium/dca/before/TOF/hDCAxyVsPtHeliumTrueMaterial", "DCAxy vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});

            histos.add<TH2>("tracks/helium/dca/before/TOF/hDCAxyVsPtantiHeliumTrue", "DCAxy vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});

            histos.add<TH2>("tracks/helium/dca/before/TOF/hDCAxyVsPtantiHeliumTruePrim", "DCAxy vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});
            histos.add<TH2>("tracks/helium/dca/before/TOF/hDCAxyVsPtantiHeliumTrueSec", "DCAxy vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});
            histos.add<TH2>("tracks/helium/dca/before/TOF/hDCAxyVsPtantiHeliumTrueMaterial", "DCAxy vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});

            histos.add<TH2>("tracks/helium/dca/before/TOF/hDCAzVsPtHeliumTrue", "DCAz vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});

            histos.add<TH2>("tracks/helium/dca/before/TOF/hDCAzVsPtHeliumTruePrim", "DCAz vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
            histos.add<TH2>("tracks/helium/dca/before/TOF/hDCAzVsPtHeliumTrueSec", "DCAz vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
            histos.add<TH2>("tracks/helium/dca/before/TOF/hDCAzVsPtHeliumTrueMaterial", "DCAz vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});

            histos.add<TH2>("tracks/helium/dca/before/TOF/hDCAzVsPtantiHeliumTrue", "DCAz vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});

            histos.add<TH2>("tracks/helium/dca/before/TOF/hDCAzVsPtantiHeliumTruePrim", "DCAz vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
            histos.add<TH2>("tracks/helium/dca/before/TOF/hDCAzVsPtantiHeliumTrueSec", "DCAz vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
            histos.add<TH2>("tracks/helium/dca/before/TOF/hDCAzVsPtantiHeliumTrueMaterial", "DCAz vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
          }

          // Fake & wrong histos
          if (outFlagOptions.makeFakeTracksPlots) {
            histos.add<TH2>("tracks/helium/dca/before/fake/hDCAxyVsPtHeliumTrue", "DCAxy vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});

            histos.add<TH2>("tracks/helium/dca/before/fake/hDCAxyVsPtHeliumTruePrim", "DCAxy vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});
            histos.add<TH2>("tracks/helium/dca/before/fake/hDCAxyVsPtHeliumTrueSec", "DCAxy vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});
            histos.add<TH2>("tracks/helium/dca/before/fake/hDCAxyVsPtHeliumTrueMaterial", "DCAxy vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});

            histos.add<TH2>("tracks/helium/dca/before/fake/hDCAxyVsPtantiHeliumTrue", "DCAxy vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});

            histos.add<TH2>("tracks/helium/dca/before/fake/hDCAxyVsPtantiHeliumTruePrim", "DCAxy vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});
            histos.add<TH2>("tracks/helium/dca/before/fake/hDCAxyVsPtantiHeliumTrueSec", "DCAxy vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});
            histos.add<TH2>("tracks/helium/dca/before/fake/hDCAxyVsPtantiHeliumTrueMaterial", "DCAxy vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});

            histos.add<TH2>("tracks/helium/dca/before/fake/hDCAzVsPtHeliumTrue", "DCAz vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});

            histos.add<TH2>("tracks/helium/dca/before/fake/hDCAzVsPtHeliumTruePrim", "DCAz vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
            histos.add<TH2>("tracks/helium/dca/before/fake/hDCAzVsPtHeliumTrueSec", "DCAz vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
            histos.add<TH2>("tracks/helium/dca/before/fake/hDCAzVsPtHeliumTrueTransport", "DCAz vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});

            histos.add<TH2>("tracks/helium/dca/before/fake/hDCAzVsPtantiHeliumTrue", "DCAz vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});

            histos.add<TH2>("tracks/helium/dca/before/fake/hDCAzVsPtantiHeliumTruePrim", "DCAz vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
            histos.add<TH2>("tracks/helium/dca/before/fake/hDCAzVsPtantiHeliumTrueSec", "DCAz vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
            histos.add<TH2>("tracks/helium/dca/before/fake/hDCAzVsPtantiHeliumTrueTransport", "DCAz vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
          }

          if (outFlagOptions.makeWrongEventPlots) {
            histos.add<TH2>("tracks/helium/dca/before/wrong/hDCAxyVsPtHeliumTrue", "DCAxy vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});

            histos.add<TH2>("tracks/helium/dca/before/wrong/hDCAxyVsPtHeliumTruePrim", "DCAxy vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});
            histos.add<TH2>("tracks/helium/dca/before/wrong/hDCAxyVsPtHeliumTrueSec", "DCAxy vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});
            histos.add<TH2>("tracks/helium/dca/before/wrong/hDCAxyVsPtHeliumTrueTransport", "DCAxy vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});

            histos.add<TH2>("tracks/helium/dca/before/wrong/hDCAxyVsPtantiHeliumTrue", "DCAxy vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});

            histos.add<TH2>("tracks/helium/dca/before/wrong/hDCAxyVsPtantiHeliumTruePrim", "DCAxy vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});
            histos.add<TH2>("tracks/helium/dca/before/wrong/hDCAxyVsPtantiHeliumTrueSec", "DCAxy vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});
            histos.add<TH2>("tracks/helium/dca/before/wrong/hDCAxyVsPtantiHeliumTrueTransport", "DCAxy vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});

            histos.add<TH2>("tracks/helium/dca/before/wrong/hDCAzVsPtHeliumTrue", "DCAz vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});

            histos.add<TH2>("tracks/helium/dca/before/wrong/hDCAzVsPtHeliumTruePrim", "DCAz vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
            histos.add<TH2>("tracks/helium/dca/before/wrong/hDCAzVsPtHeliumTrueSec", "DCAz vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
            histos.add<TH2>("tracks/helium/dca/before/wrong/hDCAzVsPtHeliumTrueTransport", "DCAz vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});

            histos.add<TH2>("tracks/helium/dca/before/wrong/hDCAzVsPtantiHeliumTrue", "DCAz vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});

            histos.add<TH2>("tracks/helium/dca/before/wrong/hDCAzVsPtantiHeliumTruePrim", "DCAz vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
            histos.add<TH2>("tracks/helium/dca/before/wrong/hDCAzVsPtantiHeliumTrueSec", "DCAz vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
            histos.add<TH2>("tracks/helium/dca/before/wrong/hDCAzVsPtantiHeliumTrueTransport", "DCAz vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
          }
          if (outFlagOptions.doTOFplots) {
            if (outFlagOptions.makeFakeTracksPlots) {
              histos.add<TH2>("tracks/helium/dca/before/fake/TOF/hDCAxyVsPtHeliumTrue", "DCAxy vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});

              histos.add<TH2>("tracks/helium/dca/before/fake/TOF/hDCAxyVsPtHeliumTruePrim", "DCAxy vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});
              histos.add<TH2>("tracks/helium/dca/before/fake/TOF/hDCAxyVsPtHeliumTrueSec", "DCAxy vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});
              histos.add<TH2>("tracks/helium/dca/before/fake/TOF/hDCAxyVsPtHeliumTrueMaterial", "DCAxy vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});

              histos.add<TH2>("tracks/helium/dca/before/fake/TOF/hDCAxyVsPtantiHeliumTrue", "DCAxy vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});

              histos.add<TH2>("tracks/helium/dca/before/fake/TOF/hDCAxyVsPtantiHeliumTruePrim", "DCAxy vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});
              histos.add<TH2>("tracks/helium/dca/before/fake/TOF/hDCAxyVsPtantiHeliumTrueSec", "DCAxy vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});
              histos.add<TH2>("tracks/helium/dca/before/fake/TOF/hDCAxyVsPtantiHeliumTrueMaterial", "DCAxy vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});

              histos.add<TH2>("tracks/helium/dca/before/fake/TOF/hDCAzVsPtHeliumTrue", "DCAz vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});

              histos.add<TH2>("tracks/helium/dca/before/fake/TOF/hDCAzVsPtHeliumTruePrim", "DCAz vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
              histos.add<TH2>("tracks/helium/dca/before/fake/TOF/hDCAzVsPtHeliumTrueSec", "DCAz vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
              histos.add<TH2>("tracks/helium/dca/before/fake/TOF/hDCAzVsPtHeliumTrueTransport", "DCAz vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});

              histos.add<TH2>("tracks/helium/dca/before/fake/TOF/hDCAzVsPtantiHeliumTrue", "DCAz vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});

              histos.add<TH2>("tracks/helium/dca/before/fake/TOF/hDCAzVsPtantiHeliumTruePrim", "DCAz vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
              histos.add<TH2>("tracks/helium/dca/before/fake/TOF/hDCAzVsPtantiHeliumTrueSec", "DCAz vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
              histos.add<TH2>("tracks/helium/dca/before/fake/TOF/hDCAzVsPtantiHeliumTrueTransport", "DCAz vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
            }

            if (outFlagOptions.makeWrongEventPlots) {
              histos.add<TH2>("tracks/helium/dca/before/wrong/TOF/hDCAxyVsPtHeliumTrue", "DCAxy vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});

              histos.add<TH2>("tracks/helium/dca/before/wrong/TOF/hDCAxyVsPtHeliumTruePrim", "DCAxy vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});
              histos.add<TH2>("tracks/helium/dca/before/wrong/TOF/hDCAxyVsPtHeliumTrueSec", "DCAxy vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});
              histos.add<TH2>("tracks/helium/dca/before/wrong/TOF/hDCAxyVsPtHeliumTrueTransport", "DCAxy vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});

              histos.add<TH2>("tracks/helium/dca/before/wrong/TOF/hDCAxyVsPtantiHeliumTrue", "DCAxy vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});

              histos.add<TH2>("tracks/helium/dca/before/wrong/TOF/hDCAxyVsPtantiHeliumTruePrim", "DCAxy vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});
              histos.add<TH2>("tracks/helium/dca/before/wrong/TOF/hDCAxyVsPtantiHeliumTrueSec", "DCAxy vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});
              histos.add<TH2>("tracks/helium/dca/before/wrong/TOF/hDCAxyVsPtantiHeliumTrueTransport", "DCAxy vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});

              histos.add<TH2>("tracks/helium/dca/before/wrong/TOF/hDCAzVsPtHeliumTrue", "DCAz vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});

              histos.add<TH2>("tracks/helium/dca/before/wrong/TOF/hDCAzVsPtHeliumTruePrim", "DCAz vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
              histos.add<TH2>("tracks/helium/dca/before/wrong/TOF/hDCAzVsPtHeliumTrueSec", "DCAz vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
              histos.add<TH2>("tracks/helium/dca/before/wrong/TOF/hDCAzVsPtHeliumTrueTransport", "DCAz vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});

              histos.add<TH2>("tracks/helium/dca/before/wrong/TOF/hDCAzVsPtantiHeliumTrue", "DCAz vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});

              histos.add<TH2>("tracks/helium/dca/before/wrong/TOF/hDCAzVsPtantiHeliumTruePrim", "DCAz vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
              histos.add<TH2>("tracks/helium/dca/before/wrong/TOF/hDCAzVsPtantiHeliumTrueSec", "DCAz vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
              histos.add<TH2>("tracks/helium/dca/before/wrong/TOF/hDCAzVsPtantiHeliumTrueTransport", "DCAz vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
            }
          }
        }

        if (outFlagOptions.makeDCAAfterCutPlots) {
          histos.add<TH2>("tracks/helium/dca/after/hDCAxyVsPtHeliumTrue", "DCAxy vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});

          histos.add<TH2>("tracks/helium/dca/after/hDCAxyVsPtHeliumTruePrim", "DCAxy vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/helium/dca/after/hDCAxyVsPtHeliumTrueSec", "DCAxy vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/helium/dca/after/hDCAxyVsPtHeliumTrueTransport", "DCAxy vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});

          histos.add<TH2>("tracks/helium/dca/after/hDCAxyVsPtantiHeliumTrue", "DCAxy vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});

          histos.add<TH2>("tracks/helium/dca/after/hDCAxyVsPtantiHeliumTruePrim", "DCAxy vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/helium/dca/after/hDCAxyVsPtantiHeliumTrueSec", "DCAxy vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/helium/dca/after/hDCAxyVsPtantiHeliumTrueTransport", "DCAxy vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});

          histos.add<TH2>("tracks/helium/dca/after/hDCAzVsPtHeliumTrue", "DCAz vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});

          histos.add<TH2>("tracks/helium/dca/after/hDCAzVsPtHeliumTruePrim", "DCAz vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/helium/dca/after/hDCAzVsPtHeliumTrueSec", "DCAz vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/helium/dca/after/hDCAzVsPtHeliumTrueTransport", "DCAz vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});

          histos.add<TH2>("tracks/helium/dca/after/hDCAzVsPtantiHeliumTrue", "DCAz vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});

          histos.add<TH2>("tracks/helium/dca/after/hDCAzVsPtantiHeliumTruePrim", "DCAz vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/helium/dca/after/hDCAzVsPtantiHeliumTrueSec", "DCAz vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/helium/dca/after/hDCAzVsPtantiHeliumTrueTransport", "DCAz vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});

          if (outFlagOptions.doTOFplots) {
            histos.add<TH2>("tracks/helium/dca/after/TOF/hDCAxyVsPtHeliumTrue", "DCAxy vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});

            histos.add<TH2>("tracks/helium/dca/after/TOF/hDCAxyVsPtHeliumTruePrim", "DCAxy vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});
            histos.add<TH2>("tracks/helium/dca/after/TOF/hDCAxyVsPtHeliumTrueSec", "DCAxy vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});
            histos.add<TH2>("tracks/helium/dca/after/TOF/hDCAxyVsPtHeliumTrueTransport", "DCAxy vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});

            histos.add<TH2>("tracks/helium/dca/after/TOF/hDCAxyVsPtantiHeliumTrue", "DCAxy vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});

            histos.add<TH2>("tracks/helium/dca/after/TOF/hDCAxyVsPtantiHeliumTruePrim", "DCAxy vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});
            histos.add<TH2>("tracks/helium/dca/after/TOF/hDCAxyVsPtantiHeliumTrueSec", "DCAxy vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});
            histos.add<TH2>("tracks/helium/dca/after/TOF/hDCAxyVsPtantiHeliumTrueTransport", "DCAxy vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcaxyAxis}});

            histos.add<TH2>("tracks/helium/dca/after/TOF/hDCAzVsPtHeliumTrue", "DCAz vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});

            histos.add<TH2>("tracks/helium/dca/after/TOF/hDCAzVsPtHeliumTruePrim", "DCAz vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
            histos.add<TH2>("tracks/helium/dca/after/TOF/hDCAzVsPtHeliumTrueSec", "DCAz vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
            histos.add<TH2>("tracks/helium/dca/after/TOF/hDCAzVsPtHeliumTrueTransport", "DCAz vs Pt (He); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});

            histos.add<TH2>("tracks/helium/dca/after/TOF/hDCAzVsPtantiHeliumTrue", "DCAz vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});

            histos.add<TH2>("tracks/helium/dca/after/TOF/hDCAzVsPtantiHeliumTruePrim", "DCAz vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
            histos.add<TH2>("tracks/helium/dca/after/TOF/hDCAzVsPtantiHeliumTrueSec", "DCAz vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
            histos.add<TH2>("tracks/helium/dca/after/TOF/hDCAzVsPtantiHeliumTrueTransport", "DCAz vs Pt (#bar{He}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptZHeAxis}, {dcazAxis}});
          }
        }
      }
      if (enableAl) {
        if (outFlagOptions.makeDCABeforeCutPlots) {
          histos.add<TH2>("tracks/alpha/dca/before/hDCAxyVsPtAlphaTrue", "DCAxy vs Pt (#alpha); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

          histos.add<TH2>("tracks/alpha/dca/before/hDCAxyVsPtAlphaTruePrim", "DCAxy vs Pt (#alpha); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/alpha/dca/before/hDCAxyVsPtAlphaTrueSec", "DCAxy vs Pt (#alpha); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/alpha/dca/before/hDCAxyVsPtAlphaTrueMaterial", "DCAxy vs Pt (#alpha); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

          histos.add<TH2>("tracks/alpha/dca/before/hDCAxyVsPtantiAlphaTrue", "DCAxy vs Pt (#bar{#alpha}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

          histos.add<TH2>("tracks/alpha/dca/before/hDCAxyVsPtantiAlphaTruePrim", "DCAxy vs Pt (#bar{#alpha}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/alpha/dca/before/hDCAxyVsPtantiAlphaTrueSec", "DCAxy vs Pt (#bar{#alpha}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/alpha/dca/before/hDCAxyVsPtantiAlphaTrueMaterial", "DCAxy vs Pt (#bar{#alpha}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
        }

        if (outFlagOptions.makeDCAAfterCutPlots) {
          histos.add<TH2>("tracks/alpha/dca/after/hDCAxyVsPtAlphaTrue", "DCAxy vs Pt (#alpha); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

          histos.add<TH2>("tracks/alpha/dca/after/hDCAxyVsPtAlphaTruePrim", "DCAxy vs Pt (#alpha); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/alpha/dca/after/hDCAxyVsPtAlphaTrueSec", "DCAxy vs Pt (#alpha); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          // histos.add<TH2>("tracks/alpha/dca/after/hDCAxyVsPtAlphaTrueMaterial", "DCAxy vs Pt (#alpha); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/alpha/dca/after/hDCAxyVsPtAlphaTrueTransport", "DCAxy vs Pt (#alpha); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

          histos.add<TH2>("tracks/alpha/dca/after/hDCAxyVsPtantiAlphaTrue", "DCAxy vs Pt (#bar{#alpha}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});

          histos.add<TH2>("tracks/alpha/dca/after/hDCAxyVsPtantiAlphaTruePrim", "DCAxy vs Pt (#bar{#alpha}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/alpha/dca/after/hDCAxyVsPtantiAlphaTrueSec", "DCAxy vs Pt (#bar{#alpha}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          // histos.add<TH2>("tracks/alpha/dca/after/hDCAxyVsPtantiAlphaTrueMaterial", "DCAxy vs Pt (#bar{#alpha}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
          histos.add<TH2>("tracks/alpha/dca/after/hDCAxyVsPtantiAlphaTrueTransport", "DCAxy vs Pt (#bar{#alpha}); #it{p}_{T} (GeV/#it{c}); DCAxy (cm)", HistType::kTH2F, {{ptAxis}, {dcaxyAxis}});
        }
      }
      // 2D-DCAz
      if (enablePr) {
        if (outFlagOptions.makeDCABeforeCutPlots) {
          histos.add<TH2>("tracks/proton/dca/before/hDCAzVsPtProtonTrue", "DCAz vs Pt (p); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

          histos.add<TH2>("tracks/proton/dca/before/hDCAzVsPtProtonTruePrim", "DCAz vs Pt (p); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/proton/dca/before/hDCAzVsPtProtonTrueSec", "DCAz vs Pt (p); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/proton/dca/before/hDCAzVsPtProtonTrueMaterial", "DCAz vs Pt (p); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

          histos.add<TH2>("tracks/proton/dca/before/hDCAzVsPtantiProtonTrue", "DCAz vs Pt (#bar{p}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

          histos.add<TH2>("tracks/proton/dca/before/hDCAzVsPtantiProtonTruePrim", "DCAz vs Pt (#bar{p}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/proton/dca/before/hDCAzVsPtantiProtonTrueSec", "DCAz vs Pt (#bar{p}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/proton/dca/before/hDCAzVsPtantiProtonTrueMaterial", "DCAz vs Pt (#bar{p}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
        }
        if (outFlagOptions.makeDCAAfterCutPlots) {
          histos.add<TH2>("tracks/proton/dca/after/hDCAzVsPtProtonTrue", "DCAz vs Pt (p); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

          histos.add<TH2>("tracks/proton/dca/after/hDCAzVsPtProtonTruePrim", "DCAz vs Pt (p); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/proton/dca/after/hDCAzVsPtProtonTrueSec", "DCAz vs Pt (p); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/proton/dca/after/hDCAzVsPtProtonTrueTransport", "DCAz vs Pt (p); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

          histos.add<TH2>("tracks/proton/dca/after/hDCAzVsPtantiProtonTrue", "DCAz vs Pt (#bar{p}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

          histos.add<TH2>("tracks/proton/dca/after/hDCAzVsPtantiProtonTruePrim", "DCAz vs Pt (#bar{p}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/proton/dca/after/hDCAzVsPtantiProtonTrueSec", "DCAz vs Pt (#bar{p}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/proton/dca/after/hDCAzVsPtantiProtonTrueTransport", "DCAz vs Pt (#bar{p}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
        }
      }
      if (enableDe) {
        if (outFlagOptions.makeDCABeforeCutPlots) {
          histos.add<TH2>("tracks/deuteron/dca/before/hDCAzVsPtDeuteronTrue", "DCAz vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

          histos.add<TH2>("tracks/deuteron/dca/before/hDCAzVsPtDeuteronTruePrim", "DCAz vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/deuteron/dca/before/hDCAzVsPtDeuteronTrueSec", "DCAz vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/deuteron/dca/before/hDCAzVsPtDeuteronTrueMaterial", "DCAz vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

          histos.add<TH2>("tracks/deuteron/dca/before/hDCAzVsPtantiDeuteronTrue", "DCAz vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

          histos.add<TH2>("tracks/deuteron/dca/before/hDCAzVsPtantiDeuteronTruePrim", "DCAz vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/deuteron/dca/before/hDCAzVsPtantiDeuteronTrueSec", "DCAz vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/deuteron/dca/before/hDCAzVsPtantiDeuteronTrueMaterial", "DCAz vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

          if (outFlagOptions.makeFakeTracksPlots) {
            histos.add<TH2>("tracks/deuteron/dca/before/fake/hDCAzVsPtDeuteronTrue", "DCAz vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

            histos.add<TH2>("tracks/deuteron/dca/before/fake/hDCAzVsPtDeuteronTruePrim", "DCAz vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
            histos.add<TH2>("tracks/deuteron/dca/before/fake/hDCAzVsPtDeuteronTrueSec", "DCAz vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
            histos.add<TH2>("tracks/deuteron/dca/before/fake/hDCAzVsPtDeuteronTrueTransport", "DCAz vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

            histos.add<TH2>("tracks/deuteron/dca/before/fake/hDCAzVsPtantiDeuteronTrue", "DCAz vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

            histos.add<TH2>("tracks/deuteron/dca/before/fake/hDCAzVsPtantiDeuteronTruePrim", "DCAz vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
            histos.add<TH2>("tracks/deuteron/dca/before/fake/hDCAzVsPtantiDeuteronTrueSec", "DCAz vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
            histos.add<TH2>("tracks/deuteron/dca/before/fake/hDCAzVsPtantiDeuteronTrueTransport", "DCAz vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          }
        }

        if (outFlagOptions.makeDCAAfterCutPlots) {
          histos.add<TH2>("tracks/deuteron/dca/after/hDCAzVsPtDeuteronTrue", "DCAz vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

          histos.add<TH2>("tracks/deuteron/dca/after/hDCAzVsPtDeuteronTruePrim", "DCAz vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/deuteron/dca/after/hDCAzVsPtDeuteronTrueWPIDPrim", "DCAz vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/deuteron/dca/after/hDCAzVsPtDeuteronTrueSec", "DCAz vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/deuteron/dca/after/hDCAzVsPtDeuteronTrueTransport", "DCAz vs Pt (d); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

          histos.add<TH2>("tracks/deuteron/dca/after/hDCAzVsPtantiDeuteronTrue", "DCAz vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

          histos.add<TH2>("tracks/deuteron/dca/after/hDCAzVsPtantiDeuteronTruePrim", "DCAz vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/deuteron/dca/after/hDCAzVsPtantiDeuteronTrueWPIDPrim", "DCAz vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/deuteron/dca/after/hDCAzVsPtantiDeuteronTrueSec", "DCAz vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/deuteron/dca/after/hDCAzVsPtantiDeuteronTrueTransport", "DCAz vs Pt (#bar{d}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
        }
      }
      if (enableTr) {
        if (outFlagOptions.makeDCABeforeCutPlots) {
          histos.add<TH2>("tracks/triton/dca/before/hDCAzVsPtTritonTrue", "DCAz vs Pt (t); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

          histos.add<TH2>("tracks/triton/dca/before/hDCAzVsPtTritonTruePrim", "DCAz vs Pt (t); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/triton/dca/before/hDCAzVsPtTritonTrueSec", "DCAz vs Pt (t); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/triton/dca/before/hDCAzVsPtTritonTrueMaterial", "DCAz vs Pt (t); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

          histos.add<TH2>("tracks/triton/dca/before/hDCAzVsPtantiTritonTrue", "DCAz vs Pt (#bar{t}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

          histos.add<TH2>("tracks/triton/dca/before/hDCAzVsPtantiTritonTruePrim", "DCAz vs Pt (#bar{t}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/triton/dca/before/hDCAzVsPtantiTritonTrueSec", "DCAz vs Pt (#bar{t}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/triton/dca/before/hDCAzVsPtantiTritonTrueMaterial", "DCAz vs Pt (#bar{t}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
        }

        if (outFlagOptions.makeDCAAfterCutPlots) {
          histos.add<TH2>("tracks/triton/dca/after/hDCAzVsPtTritonTrue", "DCAz vs Pt (#bar{t}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

          histos.add<TH2>("tracks/triton/dca/after/hDCAzVsPtTritonTruePrim", "DCAz vs Pt (t); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/triton/dca/after/hDCAzVsPtTritonTrueSec", "DCAz vs Pt (t); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/triton/dca/after/hDCAzVsPtTritonTrueTransport", "DCAz vs Pt (t); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

          histos.add<TH2>("tracks/triton/dca/after/hDCAzVsPtantiTritonTrue", "DCAz vs Pt (#bar{t}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

          histos.add<TH2>("tracks/triton/dca/after/hDCAzVsPtantiTritonTruePrim", "DCAz vs Pt (#bar{t}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/triton/dca/after/hDCAzVsPtantiTritonTrueSec", "DCAz vs Pt (#bar{t}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/triton/dca/after/hDCAzVsPtantiTritonTrueTransport", "DCAz vs Pt (#bar{t}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
        }
      }

      if (enableAl) {
        if (outFlagOptions.makeDCABeforeCutPlots) {
          histos.add<TH2>("tracks/alpha/dca/before/hDCAzVsPtAlphaTrue", "DCAz vs Pt (#alpha); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

          histos.add<TH2>("tracks/alpha/dca/before/hDCAzVsPtAlphaTruePrim", "DCAz vs Pt (#alpha); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/alpha/dca/before/hDCAzVsPtAlphaTrueSec", "DCAz vs Pt (#alpha); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/alpha/dca/before/hDCAzVsPtAlphaTrueMaterial", "DCAz vs Pt (#alpha); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

          histos.add<TH2>("tracks/alpha/dca/before/hDCAzVsPtantiAlphaTrue", "DCAz vs Pt (#bar{#alpha}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

          histos.add<TH2>("tracks/alpha/dca/before/hDCAzVsPtantiAlphaTruePrim", "DCAz vs Pt (#bar{#alpha}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/alpha/dca/before/hDCAzVsPtantiAlphaTrueSec", "DCAz vs Pt (#bar{#alpha}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/alpha/dca/before/hDCAzVsPtantiAlphaTrueMaterial", "DCAz vs Pt (#bar{#alpha}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
        }

        if (outFlagOptions.makeDCAAfterCutPlots) {
          histos.add<TH2>("tracks/alpha/dca/after/hDCAzVsPtAlphaTrue", "DCAz vs Pt (#alpha); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

          histos.add<TH2>("tracks/alpha/dca/after/hDCAzVsPtAlphaTruePrim", "DCAz vs Pt (#alpha); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/alpha/dca/after/hDCAzVsPtAlphaTrueSec", "DCAz vs Pt (#alpha); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/alpha/dca/after/hDCAzVsPtAlphaTrueTransport", "DCAz vs Pt (#alpha); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

          histos.add<TH2>("tracks/alpha/dca/after/hDCAzVsPtantiAlphaTrue", "DCAz vs Pt (#bar{#alpha}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});

          histos.add<TH2>("tracks/alpha/dca/after/hDCAzVsPtantiAlphaTruePrim", "DCAz vs Pt (#bar{#alpha}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/alpha/dca/after/hDCAzVsPtantiAlphaTrueSec", "DCAz vs Pt (#bar{#alpha}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
          histos.add<TH2>("tracks/alpha/dca/after/hDCAzVsPtantiAlphaTrueTransport", "DCAz vs Pt (#bar{#alpha}); #it{p}_{T} (GeV/#it{c}); DCAz (cm)", HistType::kTH2F, {{ptAxis}, {dcazAxis}});
        }
      }
    }

    //  Bethe-Bloch TPC distribution and Beta vs pT TOF distribution
    if (nsigmaITSvar.showAverageClusterSize && outFlagOptions.enablePIDplot) {
      histos.add<TH3>("tracks/avgClusterSizePerCoslInvVsITSlayers", "", HistType::kTH3F, {{pZAxis}, {avClsEffAxis}, {8, -0.5, 7.5}});
      histos.add<TH2>("tracks/averageClusterSize", "", HistType::kTH2F, {{pZAxis}, {avClsAxis}});
      histos.add<TH2>("tracks/averageClusterSizePerCoslInv", "", HistType::kTH2F, {{pZAxis}, {avClsEffAxis}});
    }
    if (enableDebug) {
      debugHistos.add<TH2>("debug/h2TPCsignVsTPCmomentum_AllTracks", "TPC <-dE/dX> vs #it{p}/Z (w/o rejection); Signed #it{p} (GeV/#it{c}); TPC <-dE/dx> (a.u.)", HistType::kTH2F, {{400, -8.f, 8.f}, {dedxAxis}});
      debugHistos.add<TH2>("debug/h2TPCsignVsTPCmomentum_FakeHits", "TPC <-dE/dX> vs #it{p}/Z (Fake hits); Signed #it{p} (GeV/#it{c}); TPC <-dE/dx> (a.u.)", HistType::kTH2F, {{400, -8.f, 8.f}, {dedxAxis}});
    }
    if (outFlagOptions.enablePIDplot) {
      histos.add<TH2>("tracks/h2TPCsignVsTPCmomentum", "TPC <-dE/dX> vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TPC <-dE/dx> (a.u.)", HistType::kTH2F, {{400, -8.f, 8.f}, {dedxAxis}});
      if (enablePr) {
        histos.add<TH2>("tracks/proton/h2TPCsignVsTPCmomentumProton", "TPC <-dE/dX> vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TPC <-dE/dx> (a.u.)", HistType::kTH2F, {{200, 0.f, 8.f}, {dedxAxis}});
        histos.add<TH2>("tracks/proton/h2TPCsignVsTPCmomentumantiProton", "TPC <-dE/dX> vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TPC <-dE/dx> (a.u.)", HistType::kTH2F, {{200, 0.f, 8.f}, {dedxAxis}});
      }
      if (enableDe) {
        histos.add<TH2>("tracks/deuteron/h2TPCsignVsTPCmomentumDeuteron", "TPC <-dE/dX> vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TPC <-dE/dx> (a.u.)", HistType::kTH2F, {{200, 0.f, 8.f}, {dedxAxis}});
        histos.add<TH2>("tracks/deuteron/h2TPCsignVsTPCmomentumantiDeuteron", "TPC <-dE/dX> vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TPC <-dE/dx> (a.u.)", HistType::kTH2F, {{200, 0.f, 8.f}, {dedxAxis}});
      }
      if (enableTr) {
        histos.add<TH2>("tracks/triton/h2TPCsignVsTPCmomentumTriton", "TPC <-dE/dX> vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TPC <-dE/dx> (a.u.)", HistType::kTH2F, {{200, 0.f, 8.f}, {dedxAxis}});
        histos.add<TH2>("tracks/triton/h2TPCsignVsTPCmomentumantiTriton", "TPC <-dE/dX> vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TPC <-dE/dx> (a.u.)", HistType::kTH2F, {{200, 0.f, 8.f}, {dedxAxis}});
      }
      if (enableHe) {
        histos.add<TH2>("tracks/helium/h2TPCsignVsTPCmomentumHelium", "TPC <-dE/dX> vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TPC <-dE/dx> (a.u.)", HistType::kTH2F, {{200, 0.f, 8.f}, {dedxAxis}});
        histos.add<TH2>("tracks/helium/h2TPCsignVsTPCmomentumantiHelium", "TPC <-dE/dX> vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TPC <-dE/dx> (a.u.)", HistType::kTH2F, {{200, 0.f, 8.f}, {dedxAxis}});
      }
      if (enableAl) {
        histos.add<TH2>("tracks/alpha/h2TPCsignVsTPCmomentumAlpha", "TPC <-dE/dX> vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TPC <-dE/dx> (a.u.)", HistType::kTH2F, {{200, 0.f, 8.f}, {dedxAxis}});
        histos.add<TH2>("tracks/alpha/h2TPCsignVsTPCmomentumantiAlpha", "TPC <-dE/dX> vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TPC <-dE/dx> (a.u.)", HistType::kTH2F, {{200, 0.f, 8.f}, {dedxAxis}});
      }
    }

    if (enableHe) {
      if (nsigmaITSvar.showAverageClusterSize) {
        histos.add<TH2>("tracks/helium/averageClusterSize", "", HistType::kTH2F, {{pZAxis}, {avClsAxis}});
        histos.add<TH2>("tracks/helium/averageClusterSizePerCoslInv", "", HistType::kTH2F, {{pZAxis}, {avClsEffAxis}});
      }
    }

    if (outFlagOptions.doTOFplots && outFlagOptions.enablePIDplot) {
      histos.add<TH2>("tracks/h2TPCsignVsBetaGamma", "TPC <-dE/dX> vs #beta#gamma/Z; Signed #beta#gamma; TPC <-dE/dx> (a.u.)", HistType::kTH2F, {{250, -5.f, 5.f}, {dedxAxis}});
      histos.add<TH2>("tracks/h2TOFbetaVsP", "TOF #beta vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TOF #beta", HistType::kTH2F, {{250, -5.f, 5.f}, {betaAxis}});
      if (outFlagOptions.enableBetaCut)
        histos.add<TH2>("tracks/h2TOFbetaVsP_BetaCut", "TOF #beta vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TOF #beta", HistType::kTH2F, {{250, -5.f, 5.f}, {betaAxis}});
    }

    if (outFlagOptions.enableExpSignalTPC) {
      //  TPCExpSignal histograms
      if (enablePr) {
        histos.add<TH2>("tracks/proton/h2ProtonTPCExpSignalDiffVsPt", "TPC <-dE/dX> - Exp <-dE/dX> (p) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TPC <-dE/dX> ExpDiff (p)", HistType::kTH2F, {{ptAxis}, {16000, -800, 800.}});
        histos.add<TH2>("tracks/proton/h2antiProtonTPCExpSignalDiffVsPt", "TPC <-dE/dX> - Exp <-dE/dX> (#bar{p}) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TPC <-dE/dX> ExpDiff (#bar{p})", HistType::kTH2F, {{ptAxis}, {16000, -800, 800.}});
      }
      if (enableDe) {
        histos.add<TH2>("tracks/deuteron/h2DeuteronTPCExpSignalDiffVsPt", "TPC <-dE/dX> - Exp <-dE/dX> (d) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TPC <-dE/dX> ExpDiff (d)", HistType::kTH2F, {{ptAxis}, {16000, -800, 800.}});
        histos.add<TH2>("tracks/deuteron/h2antiDeuteronTPCExpSignalDiffVsPt", "TPC <-dE/dX> - Exp <-dE/dX> (#bar{d}) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TPC <-dE/dX> ExpDiff (#bar{d})", HistType::kTH2F, {{ptAxis}, {16000, -800, 800.}});
      }
      if (enableHe) {
        histos.add<TH2>("tracks/helium/h2HeliumTPCExpSignalDiffVsPt", "TPC <-dE/dX> - Exp <-dE/dX> (He) vs #it{p}_{T}/z; #it{p}_{T}/z (GeV/#it{c}); TPC <-dE/dX> ExpDiff (He)", HistType::kTH2F, {{ptZHeAxis}, {16000, -800, 800.}});
        histos.add<TH2>("tracks/helium/h2antiHeliumTPCExpSignalDiffVsPt", "TPC <-dE/dX> - Exp <-dE/dX> (#bar{He}) vs #it{p}_{T}/z; #it{p}_{T}/z (GeV/#it{c}); TPC <-dE/dX> ExpDiff (#bar{He})", HistType::kTH2F, {{ptZHeAxis}, {16000, -800, 800.}});
      }
    }

    //  NSigmasTPC histograms
    if (enableDebug) {
      debugHistos.add<TH2>("debug/tracks/pion/h2PionVspTNSigmaTPC", "NSigmaTPC(pi) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptAxis}, {sigmaTPCAxis}});
      debugHistos.add<TH2>("debug/tracks/kaon/h2KaonVspTNSigmaTPC", "NSigmaTPC(Ka) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptAxis}, {sigmaTPCAxis}});
    }

    if (enablePr) {
      histos.add<TH2>("tracks/proton/h2ProtonVspTNSigmaTPC", "NSigmaTPC(p) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptAxis}, {sigmaTPCAxis}});
      histos.add<TH2>("tracks/proton/h2antiProtonVspTNSigmaTPC", "NSigmaTPC(#bar{p}) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptAxis}, {sigmaTPCAxis}});
    }
    if (enableDe) {
      histos.add<TH2>("tracks/deuteron/h2DeuteronVspNSigmaITSDe", "NSigmaITS(d) vs p/z; #it{p}/z (GeV/#it{c}); NSigmaITS(d)", HistType::kTH2F, {{pZAxis}, {sigmaITSAxis}});
      histos.add<TH2>("tracks/deuteron/h2antiDeuteronVspNSigmaITSDe", "NSigmaITS(#bar{d}) vs p/z; #it{p}/z (GeV/#it{c}); NSigmaITS(#bar{d})", HistType::kTH2F, {{pZAxis}, {sigmaITSAxis}});
      histos.add<TH2>("tracks/deuteron/h2DeuteronVspNSigmaITSDe_wTPCpid", "NSigmaITS(d) vs p/z; #it{p}/z (GeV/#it{c}); NSigmaITS(d)", HistType::kTH2F, {{pZAxis}, {sigmaITSAxis}});
      histos.add<TH2>("tracks/deuteron/h2antiDeuteronVspNSigmaITSDe_wTPCpid", "NSigmaITS(#bar{d}) vs p/z; #it{p}/z (GeV/#it{c}); NSigmaITS(#bar{d})", HistType::kTH2F, {{pZAxis}, {sigmaITSAxis}});

      if (enableCentrality) {
        histos.add<TH3>("tracks/deuteron/h3DeuteronVspTNSigmaTPCVsMult", "NSigmaTPC(d) vs pT; #it{p}_{T} (GeV/#it{c}) vs mult; NSigmaTPC", HistType::kTH3F, {{ptAxis}, {sigmaTPCAxis}, {binsPercentile}});
        histos.add<TH3>("tracks/deuteron/h3antiDeuteronVspTNSigmaTPCVsMult", "NSigmaTPC(#bar{d}) vs pT; #it{p}_{T} (GeV/#it{c}) vs mult; NSigmaTPC", HistType::kTH3F, {{ptAxis}, {sigmaTPCAxis}, {binsPercentile}});
      } else {
        histos.add<TH2>("tracks/deuteron/h2DeuteronVspTNSigmaTPC", "NSigmaTPC(d) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptAxis}, {sigmaTPCAxis}});
        histos.add<TH2>("tracks/deuteron/h2antiDeuteronVspTNSigmaTPC", "NSigmaTPC(#bar{d}) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptAxis}, {sigmaTPCAxis}});
        histos.add<TH2>("tracks/deuteron/h2DeuteronVspTNSigmaTPCTruePrim", "NSigmaTPC(d) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptAxis}, {sigmaTPCAxis}});
        histos.add<TH2>("tracks/deuteron/h2antiDeuteronVspTNSigmaTPCTruePrim", "NSigmaTPC(#bar{d}) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptAxis}, {sigmaTPCAxis}});
      }
    }
    if (enableTr) {
      histos.add<TH2>("tracks/triton/h2TritonVspTNSigmaTPC", "NSigmaTPC(t) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptAxis}, {sigmaTPCAxis}});
      histos.add<TH2>("tracks/triton/h2antiTritonVspTNSigmaTPC", "NSigmaTPC(#bar{t}) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptAxis}, {sigmaTPCAxis}});
    }
    if (enableHe) {
      histos.add<TH2>("tracks/helium/h2HeliumVspTNSigmaITSHe", "NSigmaITS(He) vs p/z; #it{p}/z (GeV/#it{c}); NSigmaITS(He)", HistType::kTH2F, {{pZAxis}, {sigmaITSAxis}});
      histos.add<TH2>("tracks/helium/h2antiHeliumVspTNSigmaITSHe", "NSigmaITS(#bar{He}) vs p/z; #it{p}/z (GeV/#it{c}); NSigmaITS(#bar{He})", HistType::kTH2F, {{pZAxis}, {sigmaITSAxis}});

      histos.add<TH2>("tracks/helium/h2HeliumVspTNSigmaITSHe_wTPCpid", "NSigmaITS(He) vs p/z; #it{p}/z (GeV/#it{c}); NSigmaITS(He)", HistType::kTH2F, {{pZAxis}, {sigmaITSAxis}});
      histos.add<TH2>("tracks/helium/h2antiHeliumVspTNSigmaITSHe_wTPCpid", "NSigmaITS(#bar{He}) vs p/z; #it{p}/z (GeV/#it{c}); NSigmaITS(#bar{He})", HistType::kTH2F, {{pZAxis}, {sigmaITSAxis}});

      if (enableCentrality) {
        histos.add<TH3>("tracks/helium/h3HeliumVspTNSigmaTPCVsMult", "NSigmaTPC(He) vs pT; #it{p}_{T} (GeV/#it{c}) vs mult; NSigmaTPC", HistType::kTH3F, {{ptZHeAxis}, {sigmaTPCAxis}, {binsPercentile}});
        histos.add<TH3>("tracks/helium/h3antiHeliumVspTNSigmaTPCVsMult", "NSigmaTPC(#bar{He}) vs pT; #it{p}_{T} (GeV/#it{c}) vs mult; NSigmaTPC", HistType::kTH3F, {{ptZHeAxis}, {sigmaTPCAxis}, {binsPercentile}});
      }

      histos.add<TH2>("tracks/helium/h2HeliumVspTNSigmaTPC", "NSigmaTPC(He) vs pT/z; #it{p}_{T}/z (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptZHeAxis}, {sigmaTPCAxis}});
      histos.add<TH2>("tracks/helium/h2antiHeliumVspTNSigmaTPC", "NSigmaTPC(#bar{He}) vs pT/z; #it{p}_{T}/z (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptZHeAxis}, {sigmaTPCAxis}});
    }
    if (enableAl) {
      histos.add<TH2>("tracks/alpha/h2AlphaVspTNSigmaTPC", "NSigmaTPC(#alpha) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptAxis}, {sigmaTPCAxis}});
      histos.add<TH2>("tracks/alpha/h2antiAlphaVspTNSigmaTPC", "NSigmaTPC(#bar{#alpha}) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptAxis}, {sigmaTPCAxis}});
    }

    // TOF plots
    if (outFlagOptions.doTOFplots) {
      // TOF beta histograms
      if (enablePr) {
        histos.add<TH2>("tracks/proton/h2ProtonTOFbetaVsP", "TOF #beta (p) vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TOF #beta (p)", HistType::kTH2F, {{250, 0.f, 5.f}, {betaAxis}});
        histos.add<TH2>("tracks/proton/h2antiProtonTOFbetaVsP", "TOF #beta (#bar{p}) vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TOF #beta (#bar{p})", HistType::kTH2F, {{250, 0.f, 5.f}, {betaAxis}});
      }
      if (enableDe) {
        histos.add<TH2>("tracks/deuteron/h2DeuteronTOFbetaVsP", "TOF #beta (d) vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TOF #beta (d)", HistType::kTH2F, {{250, 0.f, 5.f}, {betaAxis}});
        histos.add<TH2>("tracks/deuteron/h2antiDeuteronTOFbetaVsP", "TOF #beta (#bar{d}) vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TOF #beta (#bar{d})", HistType::kTH2F, {{250, 0.f, 5.f}, {betaAxis}});
      }
      if (enableTr) {
        histos.add<TH2>("tracks/triton/h2TritonTOFbetaVsP", "TOF #beta (t) vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TOF #beta (t)", HistType::kTH2F, {{250, 0.f, 5.f}, {betaAxis}});
        histos.add<TH2>("tracks/triton/h2antiTritonTOFbetaVsP", "TOF #beta (#bar{t}) vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TOF #beta (#bar{t})", HistType::kTH2F, {{250, 0.f, 5.f}, {betaAxis}});
      }
      if (enableHe) {
        histos.add<TH2>("tracks/helium/h2HeliumTOFbetaVsP", "TOF #beta (He) vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TOF #beta (He)", HistType::kTH2F, {{250, 0.f, 5.f}, {betaAxis}});
        histos.add<TH2>("tracks/helium/h2antiHeliumTOFbetaVsP", "TOF #beta (#bar{He}) vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TOF #beta (#bar{He})", HistType::kTH2F, {{250, 0.f, 5.f}, {betaAxis}});
      }
      if (outFlagOptions.enableExpSignalTOF) {
        //  TOFExpSignal histograms
        if (enablePr) {
          histos.add<TH2>("tracks/proton/h2ProtonTOFExpSignalDiffVsPt", "TOF t - t_{exp}(p) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (p)", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
          histos.add<TH2>("tracks/proton/h2antiProtonTOFExpSignalDiffVsPt", "TOF t - t_{exp}(#bar{p}) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (#bar{p})", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
          histos.add<TH2>("tracks/proton/h2ProtonTOFExpSignalDiffVsPtCut", "TOF t - t_{exp}(p) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (p)", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
          histos.add<TH2>("tracks/proton/h2antiProtonTOFExpSignalDiffVsPtCut", "TOF t - t_{exp}(#bar{p}) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (#bar{p})", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
        }
        if (enableDe) {
          histos.add<TH2>("tracks/deuteron/h2DeuteronTOFExpSignalDiffVsPt", "TOF t - t_{exp}(d) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (d)", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
          histos.add<TH2>("tracks/deuteron/h2antiDeuteronTOFExpSignalDiffVsPt", "TOF t - t_{exp}(#bar{d}) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (#bar{d})", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
          histos.add<TH2>("tracks/deuteron/h2DeuteronTOFExpSignalDiffVsPtCut", "TOF t - t_{exp}(d) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (d)", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
          histos.add<TH2>("tracks/deuteron/h2antiDeuteronTOFExpSignalDiffVsPtCut", "TOF t - t_{exp}(#bar{d}) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (#bar{d})", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
        }
        if (enableHe) {
          histos.add<TH2>("tracks/helium/h2HeliumTOFExpSignalDiffVsPt", "TOF t - t_{exp}(He) vs #it{p}_{T}/z; #it{p}_{T}/z (GeV/#it{c}); TOF t - t_{exp} (He)", HistType::kTH2F, {{ptZHeAxis}, {2000, -25000, 25000}});
          histos.add<TH2>("tracks/helium/h2antiHeliumTOFExpSignalDiffVsPt", "TOF t - t_{exp}(#bar{He}) vs #it{p}_{T}/z; #it{p}_{T}/z (GeV/#it{c}); TOF t - t_{exp} (#bar{He})", HistType::kTH2F, {{ptZHeAxis}, {2000, -25000, 25000}});
          histos.add<TH2>("tracks/helium/h2HeliumTOFExpSignalDiffVsPtCut", "TOF t - t_{exp}(He) vs #it{p}_{T}/z; #it{p}_{T}/z (GeV/#it{c}); TOF t - t_{exp} (He)", HistType::kTH2F, {{ptZHeAxis}, {2000, -25000, 25000}});
          histos.add<TH2>("tracks/helium/h2antiHeliumTOFExpSignalDiffVsPtCut", "TOF t - t_{exp}(#bar{He}}) vs #it{p}_{T}/z; #it{p}_{T}/z (GeV/#it{c}); TOF t - t_{exp} (#bar{He})", HistType::kTH2F, {{ptZHeAxis}, {2000, -25000, 25000}});
        }
      }

      // NSigmaTOF histograms
      if (enableDebug) {
        histos.add<TH2>("tracks/pion/h2PionVspTNSigmaTOF", "NSigmaTOF(pi) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {sigmaTOFAxis}});
        histos.add<TH2>("tracks/kaon/h2KaonVspTNSigmaTOF", "NSigmaTOF(Ka) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {sigmaTOFAxis}});
      }
      if (enablePr) {
        histos.add<TH2>("tracks/proton/h2ProtonVspTNSigmaTOF", "NSigmaTOF(p) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {sigmaTOFAxis}});
        histos.add<TH2>("tracks/proton/h2antiProtonVspTNSigmaTOF", "NSigmaTOF(#bar{p}) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {sigmaTOFAxis}});
      }
      if (enableDe) {
        histos.add<TH2>("tracks/deuteron/h2DeuteronVspTNSigmaTOF", "NSigmaTOF(d) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {sigmaTOFAxis}});
        histos.add<TH2>("tracks/deuteron/h2antiDeuteronVspTNSigmaTOF", "NSigmaTOF(#bar{d}) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {sigmaTOFAxis}});
      }
      if (enableHe) {
        histos.add<TH2>("tracks/helium/h2HeliumVspTNSigmaTOF", "NSigmaTOF(He) vs #it{p}_{T}/z; #it{p}_{T}/z (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptZHeAxis}, {sigmaTOFAxis}});
        histos.add<TH2>("tracks/helium/h2antiHeliumVspTNSigmaTOF", "NSigmaTOF(#bar{He}) vs #it{p}_{T}/z; #it{p}_{T}/z (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptZHeAxis}, {sigmaTOFAxis}});
      }
      // TOF mass histograms
      if (outFlagOptions.enablePIDplot)
        histos.add<TH2>("tracks/h2TOFmassVsPt", "h2TOFmassVsPt; TOFmass; #it{p}_{T} (GeV)", HistType::kTH2F, {{180, 0.4, 4.}, {250, 0., 5.}});
      if (enablePr) {
        histos.add<TH2>("tracks/proton/h2TOFmassProtonVsPt", "h2TOFmassProtonVsPt; TOFmass; #it{p}_{T} (GeV)", HistType::kTH2F, {{180, 0.4, 4.}, {250, 0., 5.}});
        histos.add<TH2>("tracks/proton/h2TOFmassantiProtonVsPt", "h2TOFmassantiProtonVsPt; TOFmass; #it{p}_{T} (GeV)", HistType::kTH2F, {{180, 0.4, 4.}, {250, 0., 5.}});
        if (outFlagOptions.enableBetaCut) {
          histos.add<TH2>("tracks/proton/h2TOFmassProtonVsPt_BetaCut", "h2TOFmassProtonVsPt_BetaCut; TOFmass; #it{p}_{T} (GeV)", HistType::kTH2F, {{180, 0.4, 4.}, {250, 0., 5.}});
          histos.add<TH2>("tracks/proton/h2TOFmassantiProtonVsPt_BetaCut", "h2TOFmassantiProtonVsPt_BetaCut; TOFmass; #it{p}_{T} (GeV)", HistType::kTH2F, {{180, 0.4, 4.}, {250, 0., 5.}});
        }
      }
      if (enableDe) {
        histos.add<TH2>("tracks/deuteron/h2TOFmassDeuteronVsPt", "h2TOFmassDeuteronVsPt; TOFmass; #it{p}_{T} (GeV)", HistType::kTH2F, {{180, 0.4, 4.}, {250, 0., 5.}});
        histos.add<TH2>("tracks/deuteron/h2TOFmassantiDeuteronVsPt", "h2TOFmassantiDeuteronVsPt; TOFmass; #it{p}_{T} (GeV)", HistType::kTH2F, {{180, 0.4, 4.}, {250, 0., 5.}});
        if (outFlagOptions.enableBetaCut) {
          histos.add<TH2>("tracks/deuteron/h2TOFmassDeuteronVsPt_BetaCut", "h2TOFmassDeuteronVsPt_BetaCut; TOFmass; #it{p}_{T} (GeV)", HistType::kTH2F, {{180, 0.4, 4.}, {250, 0., 5.}});
          histos.add<TH2>("tracks/deuteron/h2TOFmassantiDeuteronVsPt_BetaCut", "h2TOFmassantiDeuteronVsPt_BetaCut; TOFmass; #it{p}_{T} (GeV)", HistType::kTH2F, {{180, 0.4, 4.}, {250, 0., 5.}});
        }
      }
      if (enableTr) {
        histos.add<TH2>("tracks/triton/h2TOFmassTritonVsPt", "h2TOFmassTritonVsPt; TOFmass; #it{p}_{T} (GeV)", HistType::kTH2F, {{180, 0.4, 4.}, {250, 0., 5.}});
        histos.add<TH2>("tracks/triton/h2TOFmassantiTritonVsPt", "h2TOFmassantiTritonVsPt; TOFmass; #it{p}_{T} (GeV)", HistType::kTH2F, {{180, 0.4, 4.}, {250, 0., 5.}});
        if (outFlagOptions.enableBetaCut) {
          histos.add<TH2>("tracks/triton/h2TOFmassTritonVsPt_BetaCut", "h2TOFmassTritonVsPt_BetaCut; TOFmass; #it{p}_{T} (GeV)", HistType::kTH2F, {{180, 0.4, 4.}, {250, 0., 5.}});
          histos.add<TH2>("tracks/triton/h2TOFmassantiTritonVsPt_BetaCut", "h2TOFmassantiTritonVsPt_BetaCut; TOFmass; #it{p}_{T} (GeV)", HistType::kTH2F, {{180, 0.4, 4.}, {250, 0., 5.}});
        }
      }
      if (enableHe) {
        histos.add<TH2>("tracks/helium/h2TOFmassHeliumVsPt", "h2TOFmassHeliumVsPt; TOFmass; #it{p}_{T}/z (GeV)", HistType::kTH2F, {{180, 0.4, 4.}, {ptZHeAxis}});
        histos.add<TH2>("tracks/helium/h2TOFmassantiHeliumVsPt", "h2TOFmassantiHeliumVsPt; TOFmass; #it{p}_{T}/z (GeV)", HistType::kTH2F, {{180, 0.4, 4.}, {ptZHeAxis}});
        if (outFlagOptions.enableBetaCut) {
          histos.add<TH2>("tracks/helium/h2TOFmassHeliumVsPt_BetaCut", "h2TOFmassHeliumVsPt_BetaCut; TOFmass; #it{p}_{T}/z (GeV)", HistType::kTH2F, {{180, 0.4, 4.}, {ptZHeAxis}});
          histos.add<TH2>("tracks/helium/h2TOFmassantiHeliumVsPt_BetaCut", "h2TOFmassantiHeliumVsPt_BetaCut; TOFmass; #it{p}_{T}/z (GeV)", HistType::kTH2F, {{180, 0.4, 4.}, {ptZHeAxis}});
        }
      }
      // TOF mass squared histograms
      if (enablePr) {
        histos.add<TH2>("tracks/proton/h2TOFmass2ProtonVsPt", "#Delta M^{2} (p) vs #it{p}_{T}; #Delta M^{2} (p); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{massPrAxis}, {250, 0., 5.}});
        histos.add<TH2>("tracks/proton/h2TOFmass2antiProtonVsPt", "#Delta M^{2} (#bar{p}) vs #it{p}_{T}; #Delta M^{2} (#bar{p}); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{massPrAxis}, {250, 0., 5.}});
        if (outFlagOptions.enableBetaCut) {
          histos.add<TH2>("tracks/proton/h2TOFmass2ProtonVsPt_BetaCut", "#Delta M^{2} (p) vs #it{p}_{T}; #Delta M^{2} (p); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{massPrAxis}, {250, 0., 5.}});
          histos.add<TH2>("tracks/proton/h2TOFmass2antiProtonVsPt_BetaCut", "#Delta M^{2} (#bar{p}) vs #it{p}_{T}; #Delta M^{2} (#bar{p}); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{massPrAxis}, {250, 0., 5.}});
        }
      }
      if (enableDe) {
        if (enableCentrality) {
          histos.add<TH3>("tracks/deuteron/h3TOFmass2DeuteronVsPtVsMult", "#Delta M^{2} (d) vs #it{p}_{T} vs multiplicity; #Delta M^{2} (d); #it{p}_{T} (GeV/#it{c})", HistType::kTH3F, {{massDeAxis}, {250, 0., 5.}, {binsPercentile}});
          histos.add<TH3>("tracks/deuteron/h3TOFmass2antiDeuteronVsPtVsMult", "#Delta M^{2} (#bar{d}) vs #it{p}_{T} vs multiplicity; #Delta M^{2} (#bar{d}); #it{p}_{T} (GeV/#it{c})", HistType::kTH3F, {{massDeAxis}, {250, 0., 5.}, {binsPercentile}});
        } else {
          histos.add<TH2>("tracks/deuteron/h2TOFmass2DeuteronVsPt", "#Delta M^{2} (d) vs #it{p}_{T}; #Delta M^{2} (d); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{massDeAxis}, {250, 0., 5.}});
          histos.add<TH2>("tracks/deuteron/h2TOFmass2antiDeuteronVsPt", "#Delta M^{2} (#bar{d}) vs #it{p}_{T}; #Delta M^{2} (#bar{d}); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{massDeAxis}, {250, 0., 5.}});
        }
        if (outFlagOptions.enableBetaCut) {
          histos.add<TH2>("tracks/deuteron/h2TOFmass2DeuteronVsPt_BetaCut", "#Delta M^{2} (d) vs #it{p}_{T}; #Delta M^{2} (d); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{massDeAxis}, {250, 0., 5.}});
          histos.add<TH2>("tracks/deuteron/h2TOFmass2antiDeuteronVsPt_BetaCut", "#Delta M^{2} (#bar{d}) vs #it{p}_{T}; #Delta M^{2} (#bar{d}); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{massDeAxis}, {250, 0., 5.}});
        }
      }
      if (enableTr) {
        histos.add<TH2>("tracks/triton/h2TOFmass2TritonVsPt", "#Delta M^{2} (t) vs #it{p}_{T}; #Delta M^{2} (t); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{massTrAxis}, {250, 0., 5.}});
        histos.add<TH2>("tracks/triton/h2TOFmass2antiTritonVsPt", "#Delta M^{2} (#bar{t}) vs #it{p}_{T}; #Delta M^{2} (#bar{t}; #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{massTrAxis}, {250, 0., 5.}});
        if (outFlagOptions.enableBetaCut) {
          histos.add<TH2>("tracks/triton/h2TOFmass2TritonVsPt_BetaCut", "#Delta M^{2} (t) vs #it{p}_{T}; #Delta M^{2} (t); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{massTrAxis}, {250, 0., 5.}});
          histos.add<TH2>("tracks/triton/h2TOFmass2antiTritonVsPt_BetaCut", "#Delta M^{2} (#bar{t}) vs #it{p}_{T}; #Delta M^{2} (#bar{t}; #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{massTrAxis}, {250, 0., 5.}});
        }
      }
      if (enableHe) {
        if (enableCentrality) {
          histos.add<TH3>("tracks/helium/h3TOFmass2HeliumVsPtVsMult", "#Delta M^{2} (He) vs #it{p}_{T}/z; #Delta M^{2} (He); #it{p}_{T}/z (GeV/#it{c})", HistType::kTH3F, {{massHeAxis}, {ptZHeAxis}, {binsPercentile}});
          histos.add<TH3>("tracks/helium/h3TOFmass2antiHeliumVsPtVsMult", "#Delta M^{2} (#bar{He}) vs #it{p}_{T}/z; #Delta M^{2} (#bar{He}); #it{p}_{T}/z (GeV/#it{c})", HistType::kTH3F, {{massHeAxis}, {ptZHeAxis}, {binsPercentile}});
        }
        histos.add<TH2>("tracks/helium/h2TOFmass2antiHeliumVsPt", "#Delta M^{2} (#bar{He}) vs #it{p}_{T}/z; #Delta M^{2} (#bar{He}); #it{p}_{T}/z (GeV/#it{c})", HistType::kTH2F, {{massHeAxis}, {ptZHeAxis}});
        histos.add<TH2>("tracks/helium/h2TOFmass2HeliumVsPt", "#Delta M^{2} (He) vs #it{p}_{T}/z; #Delta M^{2} (He); #it{p}_{T}/z (GeV/#it{c})", HistType::kTH2F, {{massHeAxis}, {ptZHeAxis}});
        histos.add<TH2>("tracks/helium/h2TOFmassDeltaHeliumVsPt", "#Delta M (He) vs #it{p}_{T}/z; #Delta M (He); #it{p}_{T}/z (GeV/#it{c})", HistType::kTH2F, {{massHeAxis}, {ptZHeAxis}});
        histos.add<TH2>("tracks/helium/h2TOFmassDeltaantiHeliumVsPt", "#Delta M (#bar{He}) vs #it{p}_{T}/z; #Delta M (#bar{He}); #it{p}_{T}/z (GeV/#it{c})", HistType::kTH2F, {{massHeAxis}, {ptZHeAxis}});
        if (outFlagOptions.enableBetaCut) {
          histos.add<TH2>("tracks/helium/h2TOFmass2antiHeliumVsPt_BetaCut", "#Delta M^{2} (#bar{He}) vs #it{p}_{T}/z; #Delta M^{2} (#bar{He}); #it{p}_{T}/z (GeV/#it{c})", HistType::kTH2F, {{massHeAxis}, {ptZHeAxis}});
          histos.add<TH2>("tracks/helium/h2TOFmass2HeliumVsPt_BetaCut", "#Delta M^{2} (He) vs #it{p}_{T}/z; #Delta M^{2} (He); #it{p}_{T}/z (GeV/#it{c})", HistType::kTH2F, {{massHeAxis}, {ptZHeAxis}});
        }
      }
      // TOF EvTime Splitting plots
      if (filterOptions.enableEvTimeSplitting) {
        //  Bethe-Bloch TPC distribution - TOF EvTime Splitted
        evtimeHistos.add<TH2>("tracks/evtime/fill/h2TPCsignVsTPCmomentum", "TPC <-dE/dX> vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TPC <-dE/dx> (a.u.)", HistType::kTH2F, {{500, -5.f, 5.f}, {dedxAxis}});
        evtimeHistos.add<TH2>("tracks/evtime/tof/h2TPCsignVsTPCmomentum", "TPC <-dE/dX> vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TPC <-dE/dx> (a.u.)", HistType::kTH2F, {{500, -5.f, 5.f}, {dedxAxis}});
        evtimeHistos.add<TH2>("tracks/evtime/ft0/h2TPCsignVsTPCmomentum", "TPC <-dE/dX> vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TPC <-dE/dx> (a.u.)", HistType::kTH2F, {{500, -5.f, 5.f}, {dedxAxis}});
        evtimeHistos.add<TH2>("tracks/evtime/ft0tof/h2TPCsignVsTPCmomentum", "TPC <-dE/dX> vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TPC <-dE/dx> (a.u.)", HistType::kTH2F, {{500, -5.f, 5.f}, {dedxAxis}});

        // Beta vs pT TOF distribution - TOF EvTime Splitted
        evtimeHistos.add<TH2>("tracks/evtime/fill/h2TOFbetaVsP", "TOF #beta vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TOF #beta", HistType::kTH2F, {{500, -5.f, 5.f}, {betaAxis}});
        evtimeHistos.add<TH2>("tracks/evtime/tof/h2TOFbetaVsP", "TOF #beta vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TOF #beta", HistType::kTH2F, {{500, -5.f, 5.f}, {betaAxis}});
        evtimeHistos.add<TH2>("tracks/evtime/ft0/h2TOFbetaVsP", "TOF #beta vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TOF #beta", HistType::kTH2F, {{500, -5.f, 5.f}, {betaAxis}});
        evtimeHistos.add<TH2>("tracks/evtime/ft0tof/h2TOFbetaVsP", "TOF #beta vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TOF #beta", HistType::kTH2F, {{500, -5.f, 5.f}, {betaAxis}});

        if (outFlagOptions.enableExpSignalTOF) {
          //  TOFExpSignal histograms - TOF EvTime Splitted
          if (enablePr) {
            evtimeHistos.add<TH2>("tracks/evtime/fill/proton/h2ProtonTOFExpSignalDiffVsPt", "TOF t - t_{exp}(p) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (p)", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
            evtimeHistos.add<TH2>("tracks/evtime/fill/proton/h2antiProtonTOFExpSignalDiffVsPt", "TOF t - t_{exp}(#bar{p}) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (#bar{p})", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
            evtimeHistos.add<TH2>("tracks/evtime/tof/proton/h2ProtonTOFExpSignalDiffVsPt", "TOF t - t_{exp}(p) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (p)", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
            evtimeHistos.add<TH2>("tracks/evtime/tof/proton/h2antiProtonTOFExpSignalDiffVsPt", "TOF t - t_{exp}(#bar{p}) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (#bar{p})", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});

            evtimeHistos.add<TH2>("tracks/evtime/ft0/proton/h2ProtonTOFExpSignalDiffVsPt", "TOF t - t_{exp}(p) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (p)", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
            evtimeHistos.add<TH2>("tracks/evtime/ft0/proton/h2antiProtonTOFExpSignalDiffVsPt", "TOF t - t_{exp}(#bar{p}) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (#bar{p})", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
            evtimeHistos.add<TH2>("tracks/evtime/ft0tof/proton/h2ProtonTOFExpSignalDiffVsPt", "TOF t - t_{exp}(p) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (p)", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
            evtimeHistos.add<TH2>("tracks/evtime/ft0tof/proton/h2antiProtonTOFExpSignalDiffVsPt", "TOF t - t_{exp}(#bar{p}) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (#bar{p})", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
          }
          if (enableDe) {
            evtimeHistos.add<TH2>("tracks/evtime/fill/deuteron/h2DeuteronTOFExpSignalDiffVsPt", "TOF t - t_{exp}(d) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (d)", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
            evtimeHistos.add<TH2>("tracks/evtime/fill/deuteron/h2antiDeuteronTOFExpSignalDiffVsPt", "TOF t - t_{exp}(#bar{d}) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (#bar{d})", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});

            evtimeHistos.add<TH2>("tracks/evtime/tof/deuteron/h2DeuteronTOFExpSignalDiffVsPt", "TOF t - t_{exp}(d) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (d)", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
            evtimeHistos.add<TH2>("tracks/evtime/tof/deuteron/h2antiDeuteronTOFExpSignalDiffVsPt", "TOF t - t_{exp}(#bar{d}) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (#bar{d})", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});

            evtimeHistos.add<TH2>("tracks/evtime/ft0/deuteron/h2DeuteronTOFExpSignalDiffVsPt", "TOF t - t_{exp}(d) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (d)", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
            evtimeHistos.add<TH2>("tracks/evtime/ft0/deuteron/h2antiDeuteronTOFExpSignalDiffVsPt", "TOF t - t_{exp}(#bar{d}) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (#bar{d})", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});

            evtimeHistos.add<TH2>("tracks/evtime/ft0tof/deuteron/h2DeuteronTOFExpSignalDiffVsPt", "TOF t - t_{exp}(d) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (d)", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
            evtimeHistos.add<TH2>("tracks/evtime/ft0tof/deuteron/h2antiDeuteronTOFExpSignalDiffVsPt", "TOF t - t_{exp}(#bar{d}) vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (#bar{d})", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
          }
        }

        // NSigmaTOF histograms - TOF EvTime Splitted
        if (enablePr) {
          evtimeHistos.add<TH2>("tracks/evtime/fill/proton/h2ProtonVspTNSigmaTOF", "NSigmaTOF(p) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {sigmaTOFAxis}});
          evtimeHistos.add<TH2>("tracks/evtime/fill/proton/h2antiProtonVspTNSigmaTOF", "NSigmaTOF(#bar{p}) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {sigmaTOFAxis}});
          evtimeHistos.add<TH2>("tracks/evtime/tof/proton/h2ProtonVspTNSigmaTOF", "NSigmaTOF(#bar{p}) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {sigmaTOFAxis}});
          evtimeHistos.add<TH2>("tracks/evtime/tof/proton/h2antiProtonVspTNSigmaTOF", "NSigmaTOF(#bar{p}) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {sigmaTOFAxis}});
          evtimeHistos.add<TH2>("tracks/evtime/ft0/proton/h2ProtonVspTNSigmaTOF", "NSigmaTOF(p) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {sigmaTOFAxis}});
          evtimeHistos.add<TH2>("tracks/evtime/ft0/proton/h2antiProtonVspTNSigmaTOF", "NSigmaTOF(#bar{p}) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {sigmaTOFAxis}});
          evtimeHistos.add<TH2>("tracks/evtime/ft0tof/proton/h2ProtonVspTNSigmaTOF", "NSigmaTOF(p) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {sigmaTOFAxis}});
          evtimeHistos.add<TH2>("tracks/evtime/ft0tof/proton/h2antiProtonVspTNSigmaTOF", "NSigmaTOF(#bar{p}) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {sigmaTOFAxis}});
        }
        if (enableDe) {
          evtimeHistos.add<TH2>("tracks/evtime/fill/deuteron/h2DeuteronVspTNSigmaTOF", "NSigmaTOF(d) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {sigmaTOFAxis}});
          evtimeHistos.add<TH2>("tracks/evtime/fill/deuteron/h2antiDeuteronVspTNSigmaTOF", "NSigmaTOF(#bar{d}) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {sigmaTOFAxis}});
          evtimeHistos.add<TH2>("tracks/evtime/tof/deuteron/h2DeuteronVspTNSigmaTOF", "NSigmaTOF(d) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {sigmaTOFAxis}});
          evtimeHistos.add<TH2>("tracks/evtime/tof/deuteron/h2antiDeuteronVspTNSigmaTOF", "NSigmaTOF(#bar{d}) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {sigmaTOFAxis}});
          evtimeHistos.add<TH2>("tracks/evtime/ft0/deuteron/h2DeuteronVspTNSigmaTOF", "NSigmaTOF(d) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {sigmaTOFAxis}});
          evtimeHistos.add<TH2>("tracks/evtime/ft0/deuteron/h2antiDeuteronVspTNSigmaTOF", "NSigmaTOF(#bar{d}) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {sigmaTOFAxis}});
          evtimeHistos.add<TH2>("tracks/evtime/ft0tof/deuteron/h2DeuteronVspTNSigmaTOF", "NSigmaTOF(d) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {sigmaTOFAxis}});
          evtimeHistos.add<TH2>("tracks/evtime/ft0tof/deuteron/h2antiDeuteronVspTNSigmaTOF", "NSigmaTOF(#bar{d}) vs pT; #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {sigmaTOFAxis}});
        }

        // NSigmaTPC vs NSigmaTOF histograms - TOF EvTime Splitted
        if (enablePr) {
          evtimeHistos.add<TH3>("tracks/evtime/fill/proton/h3ProtonNSigmaTPCvsNSigmaTOFvsPt", "NSigmaTPC (p) vs NSigmaTOF(p); NSigmaTPC; NSigmaTOF; #it{p}_{T} (GeV/#it{c})", HistType::kTH3F, {{200, -20, 20}, {200, -20, 20}, {ptAxis}});
          evtimeHistos.add<TH3>("tracks/evtime/fill/proton/h3antiProtonNSigmaTPCvsNSigmaTOFvsPt", "NSigmaTPC(#bar{p}) vs NSigmaTOF(#bar{p}); NSigmaTPC; NSigmaTOF; #it{p}_{T} (GeV/#it{c})", HistType::kTH3F, {{200, -20, 20}, {200, -20, 20}, {ptAxis}});
          evtimeHistos.add<TH3>("tracks/evtime/tof/proton/h3ProtonNSigmaTPCvsNSigmaTOFvsPt", "NSigmaTPC (p) vs NSigmaTOF(p); NSigmaTPC; NSigmaTOF; #it{p}_{T} (GeV/#it{c})", HistType::kTH3F, {{200, -20, 20}, {200, -20, 20}, {ptAxis}});
          evtimeHistos.add<TH3>("tracks/evtime/tof/proton/h3antiProtonNSigmaTPCvsNSigmaTOFvsPt", "NSigmaTPC(#bar{p}) vs NSigmaTOF(#bar{p}); NSigmaTPC; NSigmaTOF; #it{p}_{T} (GeV/#it{c})", HistType::kTH3F, {{200, -20, 20}, {200, -20, 20}, {ptAxis}});
          evtimeHistos.add<TH3>("tracks/evtime/ft0/proton/h3ProtonNSigmaTPCvsNSigmaTOFvsPt", "NSigmaTPC (p) vs NSigmaTOF(p); NSigmaTPC; NSigmaTOF; #it{p}_{T} (GeV/#it{c})", HistType::kTH3F, {{200, -20, 20}, {200, -20, 20}, {ptAxis}});
          evtimeHistos.add<TH3>("tracks/evtime/ft0/proton/h3antiProtonNSigmaTPCvsNSigmaTOFvsPt", "NSigmaTPC(#bar{p}) vs NSigmaTOF(#bar{p}); NSigmaTPC; NSigmaTOF; #it{p}_{T} (GeV/#it{c})", HistType::kTH3F, {{200, -20, 20}, {200, -20, 20}, {ptAxis}});
          evtimeHistos.add<TH3>("tracks/evtime/ft0tof/proton/h3ProtonNSigmaTPCvsNSigmaTOFvsPt", "NSigmaTPC (p) vs NSigmaTOF(p); NSigmaTPC; NSigmaTOF; #it{p}_{T} (GeV/#it{c})", HistType::kTH3F, {{200, -20, 20}, {200, -20, 20}, {ptAxis}});
          evtimeHistos.add<TH3>("tracks/evtime/ft0tof/proton/h3antiProtonNSigmaTPCvsNSigmaTOFvsPt", "NSigmaTPC(#bar{p}) vs NSigmaTOF(#bar{p}); NSigmaTPC; NSigmaTOF; #it{p}_{T} (GeV/#it{c})", HistType::kTH3F, {{200, -20, 20}, {200, -20, 20}, {ptAxis}});
        }
        if (enableDe) {
          evtimeHistos.add<TH3>("tracks/evtime/fill/deuteron/h3DeuteronNSigmaTPCvsNSigmaTOFvsPt", "NSigmaTPC(d) vs NSigmaTOF(d); NSigmaTPC; NSigmaTOF; #it{p}_{T} (GeV/#it{c})", HistType::kTH3F, {{200, -20, 20}, {200, -20, 20}, {ptAxis}});
          evtimeHistos.add<TH3>("tracks/evtime/fill/deuteron/h3antiDeuteronNSigmaTPCvsNSigmaTOFvsPt", "NSigmaTPC(#bar{d}) vs NSigmaTOF(#bar{p}); NSigmaTPC; NSigmaTOF; #it{p}_{T} (GeV/#it{c})", HistType::kTH3F, {{200, -20, 20}, {200, -20, 20}, {ptAxis}});
          evtimeHistos.add<TH3>("tracks/evtime/tof/deuteron/h3DeuteronNSigmaTPCvsNSigmaTOFvsPt", "NSigmaTPC(d) vs NSigmaTOF(d); NSigmaTPC; NSigmaTOF; #it{p}_{T} (GeV/#it{c})", HistType::kTH3F, {{200, -20, 20}, {200, -20, 20}, {ptAxis}});
          evtimeHistos.add<TH3>("tracks/evtime/tof/deuteron/h3antiDeuteronNSigmaTPCvsNSigmaTOFvsPt", "NSigmaTPC(#bar{d}) vs NSigmaTOF(#bar{p}); NSigmaTPC; NSigmaTOF; #it{p}_{T} (GeV/#it{c})", HistType::kTH3F, {{200, -20, 20}, {200, -20, 20}, {ptAxis}});
          evtimeHistos.add<TH3>("tracks/evtime/ft0/deuteron/h3DeuteronNSigmaTPCvsNSigmaTOFvsPt", "NSigmaTPC(d) vs NSigmaTOF(d); NSigmaTPC; NSigmaTOF; #it{p}_{T} (GeV/#it{c})", HistType::kTH3F, {{200, -20, 20}, {200, -20, 20}, {ptAxis}});
          evtimeHistos.add<TH3>("tracks/evtime/ft0/deuteron/h3antiDeuteronNSigmaTPCvsNSigmaTOFvsPt", "NSigmaTPC(#bar{d}) vs NSigmaTOF(#bar{p}); NSigmaTPC; NSigmaTOF; #it{p}_{T} (GeV/#it{c})", HistType::kTH3F, {{200, -20, 20}, {200, -20, 20}, {ptAxis}});
          evtimeHistos.add<TH3>("tracks/evtime/ft0tof/deuteron/h3DeuteronNSigmaTPCvsNSigmaTOFvsPt", "NSigmaTPC(d) vs NSigmaTOF(d); NSigmaTPC; NSigmaTOF; #it{p}_{T} (GeV/#it{c})", HistType::kTH3F, {{200, -20, 20}, {200, -20, 20}, {ptAxis}});
          evtimeHistos.add<TH3>("tracks/evtime/ft0tof/deuteron/h3antiDeuteronNSigmaTPCvsNSigmaTOFvsPt", "NSigmaTPC(#bar{d}) vs NSigmaTOF(#bar{p}); NSigmaTPC; NSigmaTOF; #it{p}_{T} (GeV/#it{c})", HistType::kTH3F, {{200, -20, 20}, {200, -20, 20}, {ptAxis}});
        }

        // TOF mass histograms - TOF EvTime Splitted
        evtimeHistos.add<TH2>("tracks/evtime/fill/h2TOFmassVsPt", "h2TOFmassVsPt; TOFmass; #it{p}_{T} (GeV)", HistType::kTH2F, {{80, 0.4, 2.}, {250, 0., 5.}});
        evtimeHistos.add<TH2>("tracks/evtime/tof/h2TOFmassVsPt", "h2TOFmassVsPt; TOFmass; #it{p}_{T} (GeV)", HistType::kTH2F, {{80, 0.4, 2.}, {250, 0., 5.}});
        evtimeHistos.add<TH2>("tracks/evtime/ft0/h2TOFmassVsPt", "h2TOFmassVsPt; TOFmass; #it{p}_{T} (GeV)", HistType::kTH2F, {{80, 0.4, 2.}, {250, 0., 5.}});
        evtimeHistos.add<TH2>("tracks/evtime/ft0tof/h2TOFmassVsPt", "h2TOFmassVsPt; TOFmass; #it{p}_{T} (GeV)", HistType::kTH2F, {{80, 0.4, 2.}, {250, 0., 5.}});

        // TOF mass squared histograms - TOF EvTime Splitted
        if (enablePr) {
          evtimeHistos.add<TH2>("tracks/evtime/fill/proton/h2TOFmass2ProtonVsPt", "#Delta M^{2} (p) vs #it{p}_{T}; #Delta M^{2} (p); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{massPrAxis}, {250, 0., 5.}});
          evtimeHistos.add<TH2>("tracks/evtime/fill/proton/h2TOFmass2antiProtonVsPt", "#Delta M^{2} (#bar{p}) vs #it{p}_{T}; #Delta M^{2} (#bar{p}); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{massPrAxis}, {250, 0., 5.}});
          evtimeHistos.add<TH2>("tracks/evtime/tof/proton/h2TOFmass2ProtonVsPt", "#Delta M^{2} (p) vs #it{p}_{T}; #Delta M^{2} (p); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{massPrAxis}, {250, 0., 5.}});
          evtimeHistos.add<TH2>("tracks/evtime/tof/proton/h2TOFmass2antiProtonVsPt", "#Delta M^{2} (#bar{p}) vs #it{p}_{T}; #Delta M^{2} (#bar{p}); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{massPrAxis}, {250, 0., 5.}});
          evtimeHistos.add<TH2>("tracks/evtime/ft0/proton/h2TOFmass2ProtonVsPt", "#Delta M^{2} (p) vs #it{p}_{T}; #Delta M^{2} (p); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{massPrAxis}, {250, 0., 5.}});
          evtimeHistos.add<TH2>("tracks/evtime/ft0/proton/h2TOFmass2antiProtonVsPt", "#Delta M^{2} (#bar{p}) vs #it{p}_{T}; #Delta M^{2} (#bar{p}); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{massPrAxis}, {250, 0., 5.}});
          evtimeHistos.add<TH2>("tracks/evtime/ft0tof/proton/h2TOFmass2ProtonVsPt", "#Delta M^{2} (p) vs #it{p}_{T}; #Delta M^{2} (p); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{massPrAxis}, {250, 0., 5.}});
          evtimeHistos.add<TH2>("tracks/evtime/ft0tof/proton/h2TOFmass2antiProtonVsPt", "#Delta M^{2} (#bar{p}) vs #it{p}_{T}; #Delta M^{2} (#bar{p}); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{massPrAxis}, {250, 0., 5.}});
        }
        if (enableDe) {
          evtimeHistos.add<TH2>("tracks/evtime/fill/deuteron/h2TOFmass2DeuteronVsPt", "#Delta M^{2} (d) vs #it{p}_{T}; #Delta M^{2} (d); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{massDeAxis}, {250, 0., 5.}});
          evtimeHistos.add<TH2>("tracks/evtime/tof/deuteron/h2TOFmass2DeuteronVsPt", "#Delta M^{2} (d) vs #it{p}_{T}; #Delta M^{2} (d); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{massDeAxis}, {250, 0., 5.}});
          evtimeHistos.add<TH2>("tracks/evtime/ft0/deuteron/h2TOFmass2DeuteronVsPt", "#Delta M^{2} (d) vs #it{p}_{T}; #Delta M^{2} (d); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{massDeAxis}, {250, 0., 5.}});
          evtimeHistos.add<TH2>("tracks/evtime/ft0tof/deuteron/h2TOFmass2DeuteronVsPt", "#Delta M^{2} (d) vs #it{p}_{T}; #Delta M^{2} (d); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{massDeAxis}, {250, 0., 5.}});
          evtimeHistos.add<TH2>("tracks/evtime/fill/deuteron/h2TOFmass2antiDeuteronVsPt", "#Delta M^{2} (#bar{d}) vs #it{p}_{T}; #Delta M^{2} (#bar{d}); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{massDeAxis}, {250, 0., 5.}});
          evtimeHistos.add<TH2>("tracks/evtime/tof/deuteron/h2TOFmass2antiDeuteronVsPt", "#Delta M^{2} (#bar{d}) vs #it{p}_{T}; #Delta M^{2} (#bar{d}); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{massDeAxis}, {250, 0., 5.}});
          evtimeHistos.add<TH2>("tracks/evtime/ft0/deuteron/h2TOFmass2antiDeuteronVsPt", "#Delta M^{2} (#bar{d}) vs #it{p}_{T}; #Delta M^{2} (#bar{d}); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{massDeAxis}, {250, 0., 5.}});
          evtimeHistos.add<TH2>("tracks/evtime/ft0tof/deuteron/h2TOFmass2antiDeuteronVsPt", "#Delta M^{2} (#bar{d}) vs #it{p}_{T}; #Delta M^{2} (#bar{d}); #it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{massDeAxis}, {250, 0., 5.}});
        }
      }

      if (enableDebug) {
        debugHistos.add<TH2>("debug/qa/h2TPCncrVsPtPos", "number of crossed rows in TPC vs Pt;  #it{p}_{T} (GeV/#it{c}); TPCncr", HistType::kTH2F, {{ptAxis}, {150, 60, 170}});
        debugHistos.add<TH2>("debug/qa/h2TPCncrVsTPCsignalPos", "number of crossed rows in TPC vs Pt;  TPC <-dE/dX>; TPCncr", HistType::kTH2F, {dedxAxis, {150, 60, 170}});
        debugHistos.add<TH1>("debug/qa/h1TPCncrLowPPos", "number of crossed rows in TPC (p<0.5 GeV/c); TPCncr; counts", HistType::kTH1F, {{150, 60, 170}});
        debugHistos.add<TH1>("debug/qa/h1TPCncrMidPPos", "number of crossed rows in TPC (0.5<p<1.0 GeV/c); TPCncr; counts", HistType::kTH1F, {{150, 60, 170}});
        debugHistos.add<TH1>("debug/qa/h1TPCncrHighPPos", "number of crossed rows in TPC (p>1.0 GeV/c); TPCncr; counts", HistType::kTH1F, {{150, 60, 170}});

        debugHistos.add<TH2>("debug/qa/h2TPCncrVsPtNeg", "number of crossed rows in TPC vs Pt;  #it{p}_{T} (GeV/#it{c}); TPCncr", HistType::kTH2F, {{ptAxis}, {150, 60, 170}});
        debugHistos.add<TH2>("debug/qa/h2TPCncrVsTPCsignalNeg", "number of crossed rows in TPC vs Pt;  TPC <-dE/dX>; TPCncr", HistType::kTH2F, {dedxAxis, {150, 60, 170}});
        debugHistos.add<TH1>("debug/qa/h1TPCncrLowPNeg", "number of crossed rows in TPC (p<0.5 GeV/c); TPCncr; counts", HistType::kTH1F, {{150, 60, 170}});
        debugHistos.add<TH1>("debug/qa/h1TPCncrMidPNeg", "number of crossed rows in TPC (0.5<p<1.0 GeV/c); TPCncr; counts", HistType::kTH1F, {{150, 60, 170}});
        debugHistos.add<TH1>("debug/qa/h1TPCncrHighPNeg", "number of crossed rows in TPC (p>1.0 GeV/c); TPCncr; counts", HistType::kTH1F, {{150, 60, 170}});

        // Beta < 0.5
        //  NSigmasTPC histograms
        if (enablePr) {
          debugHistos.add<TH2>("debug/evtime/fill/proton/h2ProtonVspTNSigmaTPC_BetaCut", "NSigmaTPC(p) vs pT (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptAxis}, {sigmaTPCAxis}});
          debugHistos.add<TH2>("debug/evtime/fill/proton/h2antiProtonVspTNSigmaTPC_BetaCut", "NSigmaTPC(#bar{p}) vs pT (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptAxis}, {sigmaTPCAxis}});
          debugHistos.add<TH2>("debug/evtime/tof/proton/h2ProtonVspTNSigmaTPC_BetaCut", "NSigmaTPC(p) vs pT (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptAxis}, {sigmaTPCAxis}});
          debugHistos.add<TH2>("debug/evtime/tof/proton/h2antiProtonVspTNSigmaTPC_BetaCut", "NSigmaTPC(#bar{p}) vs pT (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptAxis}, {sigmaTPCAxis}});
          debugHistos.add<TH2>("debug/evtime/ft0/proton/h2ProtonVspTNSigmaTPC_BetaCut", "NSigmaTPC(p) vs pT (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptAxis}, {sigmaTPCAxis}});
          debugHistos.add<TH2>("debug/evtime/ft0/proton/h2antiProtonVspTNSigmaTPC_BetaCut", "NSigmaTPC(#bar{p}) vs pT (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptAxis}, {sigmaTPCAxis}});
          debugHistos.add<TH2>("debug/evtime/ft0tof/proton/h2ProtonVspTNSigmaTPC_BetaCut", "NSigmaTPC(p) vs pT (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptAxis}, {sigmaTPCAxis}});
          debugHistos.add<TH2>("debug/evtime/ft0tof/proton/h2antiProtonVspTNSigmaTPC_BetaCut", "NSigmaTPC(#bar{p}) vs pT (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptAxis}, {sigmaTPCAxis}});
        }
        if (enableDe) {
          debugHistos.add<TH2>("debug/evtime/fill/deuteron/h2DeuteronVspTNSigmaTPC_BetaCut", "NSigmaTPC(d) vs pT (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptAxis}, {sigmaTPCAxis}});
          debugHistos.add<TH2>("debug/evtime/fill/deuteron/h2antiDeuteronVspTNSigmaTPC_BetaCut", "NSigmaTPC(#bar{d}) vs pT (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptAxis}, {sigmaTPCAxis}});
          debugHistos.add<TH2>("debug/evtime/tof/deuteron/h2DeuteronVspTNSigmaTPC_BetaCut", "NSigmaTPC(d) vs pT (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptAxis}, {sigmaTPCAxis}});
          debugHistos.add<TH2>("debug/evtime/tof/deuteron/h2antiDeuteronVspTNSigmaTPC_BetaCut", "NSigmaTPC(#bar{d}) vs pT (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptAxis}, {sigmaTPCAxis}});
          debugHistos.add<TH2>("debug/evtime/ft0/deuteron/h2DeuteronVspTNSigmaTPC_BetaCut", "NSigmaTPC(d) vs pT (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptAxis}, {sigmaTPCAxis}});
          debugHistos.add<TH2>("debug/evtime/ft0/deuteron/h2antiDeuteronVspTNSigmaTPC_BetaCut", "NSigmaTPC(#bar{d}) vs pT (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptAxis}, {sigmaTPCAxis}});
          debugHistos.add<TH2>("debug/evtime/ft0tof/deuteron/h2DeuteronVspTNSigmaTPC_BetaCut", "NSigmaTPC(d) vs pT (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptAxis}, {sigmaTPCAxis}});
          debugHistos.add<TH2>("debug/evtime/ft0tof/deuteron/h2antiDeuteronVspTNSigmaTPC_BetaCut", "NSigmaTPC(#bar{d}) vs pT (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); NSigmaTPC", HistType::kTH2F, {{ptAxis}, {sigmaTPCAxis}});
        }

        // NSigmaTOF histograms - TOF EvTime Splitted
        if (enablePr) {
          debugHistos.add<TH2>("debug/evtime/fill/proton/h2ProtonVspTNSigmaTOF_BetaCut", "NSigmaTOF(p) vs pT (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {sigmaTOFAxis}});
          debugHistos.add<TH2>("debug/evtime/fill/proton/h2antiProtonVspTNSigmaTOF_BetaCut", "NSigmaTOF(#bar{p}) vs pT (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {sigmaTOFAxis}});
          debugHistos.add<TH2>("debug/evtime/tof/proton/h2ProtonVspTNSigmaTOF_BetaCut", "NSigmaTOF(p) vs pT (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {sigmaTOFAxis}});
          debugHistos.add<TH2>("debug/evtime/tof/proton/h2antiProtonVspTNSigmaTOF_BetaCut", "NSigmaTOF(#bar{p}) vs pT (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {sigmaTOFAxis}});
          debugHistos.add<TH2>("debug/evtime/ft0/proton/h2ProtonVspTNSigmaTOF_BetaCut", "NSigmaTOF(p) vs pT (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {sigmaTOFAxis}});
          debugHistos.add<TH2>("debug/evtime/ft0/proton/h2antiProtonVspTNSigmaTOF_BetaCut", "NSigmaTOF(#bar{p}) vs pT (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {sigmaTOFAxis}});
          debugHistos.add<TH2>("debug/evtime/ft0tof/proton/h2ProtonVspTNSigmaTOF_BetaCut", "NSigmaTOF(p) vs pT (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {sigmaTOFAxis}});
          debugHistos.add<TH2>("debug/evtime/ft0tof/proton/h2antiProtonVspTNSigmaTOF_BetaCut", "NSigmaTOF(#bar{p}) vs pT (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {sigmaTOFAxis}});
        }
        if (enableDe) {
          debugHistos.add<TH2>("debug/evtime/fill/deuteron/h2DeuteronVspTNSigmaTOF_BetaCut", "NSigmaTOF(d) vs pT (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {sigmaTOFAxis}});
          debugHistos.add<TH2>("debug/evtime/fill/deuteron/h2antiDeuteronVspTNSigmaTOF_BetaCut", "NSigmaTOF(#bar{d}) vs pT (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {sigmaTOFAxis}});
          debugHistos.add<TH2>("debug/evtime/tof/deuteron/h2DeuteronVspTNSigmaTOF_BetaCut", "NSigmaTOF(d) vs pT (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {sigmaTOFAxis}});
          debugHistos.add<TH2>("debug/evtime/tof/deuteron/h2antiDeuteronVspTNSigmaTOF_BetaCut", "NSigmaTOF(#bar{d}) vs pT (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {sigmaTOFAxis}});
          debugHistos.add<TH2>("debug/evtime/ft0/deuteron/h2DeuteronVspTNSigmaTOF_BetaCut", "NSigmaTOF(d) vs pT (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {sigmaTOFAxis}});
          debugHistos.add<TH2>("debug/evtime/ft0/deuteron/h2antiDeuteronVspTNSigmaTOF_BetaCut", "NSigmaTOF(#bar{d}) vs pT (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {sigmaTOFAxis}});
          debugHistos.add<TH2>("debug/evtime/ft0tof/deuteron/h2DeuteronVspTNSigmaTOF_BetaCut", "NSigmaTOF(d) vs pT (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {sigmaTOFAxis}});
          debugHistos.add<TH2>("debug/evtime/ft0tof/deuteron/h2antiDeuteronVspTNSigmaTOF_BetaCut", "NSigmaTOF(#bar{d}) vs pT (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); NSigmaTOF", HistType::kTH2F, {{ptAxis}, {sigmaTOFAxis}});
        }

        if (outFlagOptions.enableExpSignalTOF) {
          //  TOFExpSignal histograms - TOF EvTime Splitted
          if (enablePr) {
            debugHistos.add<TH2>("debug/evtime/fill/proton/h2ProtonTOFExpSignalDiffVsPt_BetaCut", "TOF t - t_{exp}(p) vs #it{p}_{T} (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (p)", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
            debugHistos.add<TH2>("debug/evtime/fill/proton/h2antiProtonTOFExpSignalDiffVsPt_BetaCut", "TOF t - t_{exp}(#bar{p}) vs #it{p}_{T} (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (#bar{p})", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
            debugHistos.add<TH2>("debug/evtime/tof/proton/h2ProtonTOFExpSignalDiffVsPt_BetaCut", "TOF t - t_{exp}(p) vs #it{p}_{T} (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (p)", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
            debugHistos.add<TH2>("debug/evtime/tof/proton/h2antiProtonTOFExpSignalDiffVsPt_BetaCut", "TOF t - t_{exp}(#bar{p}) vs #it{p}_{T} (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (#bar{p})", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
            debugHistos.add<TH2>("debug/evtime/ft0/proton/h2ProtonTOFExpSignalDiffVsPt_BetaCut", "TOF t - t_{exp}(p) vs #it{p}_{T} (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (p)", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
            debugHistos.add<TH2>("debug/evtime/ft0/proton/h2antiProtonTOFExpSignalDiffVsPt_BetaCut", "TOF t - t_{exp}(#bar{p}) vs #it{p}_{T} (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (#bar{p})", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
            debugHistos.add<TH2>("debug/evtime/ft0tof/proton/h2ProtonTOFExpSignalDiffVsPt_BetaCut", "TOF t - t_{exp}(p) vs #it{p}_{T} (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (p)", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
            debugHistos.add<TH2>("debug/evtime/ft0tof/proton/h2antiProtonTOFExpSignalDiffVsPt_BetaCut", "TOF t - t_{exp}(#bar{p}) vs #it{p}_{T} (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (#bar{p})", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
          }
          if (enableDe) {
            debugHistos.add<TH2>("debug/evtime/fill/deuteron/h2DeuteronTOFExpSignalDiffVsPt_BetaCut", "TOF t - t_{exp}(d) vs #it{p}_{T} (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (d)", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
            debugHistos.add<TH2>("debug/evtime/fill/deuteron/h2antiDeuteronTOFExpSignalDiffVsPt_BetaCut", "TOF t - t_{exp}(#bar{d}) vs #it{p}_{T} (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (#bar{d})", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
            debugHistos.add<TH2>("debug/evtime/tof/deuteron/h2DeuteronTOFExpSignalDiffVsPt_BetaCut", "TOF t - t_{exp}(d) vs #it{p}_{T} (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (d)", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
            debugHistos.add<TH2>("debug/evtime/tof/deuteron/h2antiDeuteronTOFExpSignalDiffVsPt_BetaCut", "TOF t - t_{exp}(#bar{d}) vs #it{p}_{T} (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (#bar{d})", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
            debugHistos.add<TH2>("debug/evtime/ft0/deuteron/h2DeuteronTOFExpSignalDiffVsPt_BetaCut", "TOF t - t_{exp}(d) vs #it{p}_{T} (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (d)", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
            debugHistos.add<TH2>("debug/evtime/ft0/deuteron/h2antiDeuteronTOFExpSignalDiffVsPt_BetaCut", "TOF t - t_{exp}(#bar{d}) vs #it{p}_{T} (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (#bar{d})", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
            debugHistos.add<TH2>("debug/evtime/ft0tof/deuteron/h2DeuteronTOFExpSignalDiffVsPt_BetaCut", "TOF t - t_{exp}(d) vs #it{p}_{T} (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (d)", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
            debugHistos.add<TH2>("debug/evtime/ft0tof/deuteron/h2antiDeuteronTOFExpSignalDiffVsPt_BetaCut", "TOF t - t_{exp}(#bar{d}) vs #it{p}_{T} (#beta < 0.5); #it{p}_{T} (GeV/#it{c}); TOF t - t_{exp} (#bar{d})", HistType::kTH2F, {{ptAxis}, {2000, -25000, 25000}});
          }
        }
      }
    }
    // To be optimised
    if (!doprocessMCGen && !doprocessMCReco && !doprocessMCRecoLfPid && !doprocessMCRecoFiltered && !doprocessMCRecoFilteredLight) {
      LOG(info) << "Histograms of LFNucleiBATask:";
      histos.print();
      return;
    }
    // MC histograms  -   all, primary, sec. from weak decay, sec. from material
    if (enableCentrality)
      spectraGen.add("histGenVetxZ", "PosZ generated events", HistType::kTH2F, {{1500, -15.f, 15.f, "Vertex Z (cm)"}, {binsPercentile, "Centrality FT0M"}});
    else
      spectraGen.add("histGenVetxZ", "PosZ generated events", HistType::kTH1F, {{1500, -15.f, 15.f, "Vertex Z (cm)"}});

    spectraGen.add("helium/histPtGenHe", "PtGenHe", HistType::kTH1F, {{800, 0.f, 8.f}});
    spectraGen.add("helium/histPtRecHe", "PtRecHe", HistType::kTH1F, {{800, 0.f, 8.f}});
    spectraGen.add("helium/histPtShiftHe", "PtReco-PtGen vs PtReco", HistType::kTH2F, {{800, 0.f, 8.f}, {400, -4.f, 4.f}});
    spectraGen.add("helium/histPtShiftVsEtaHe", "PtReco-PtGen vs #eta", HistType::kTH2F, {{200, -2.f, 2.f}, {400, -4.f, 4.f}});

    spectraGen.add("helium/histPGenHe", "PGenHe", HistType::kTH1F, {{800, 0.f, 8.f}});
    spectraGen.add("helium/histPRecHe", "PRecHe", HistType::kTH1F, {{800, 0.f, 8.f}});
    spectraGen.add("helium/histPShiftHe", "PReco-PGen vs PReco", HistType::kTH2F, {{800, 0.f, 8.f}, {400, -4.f, 4.f}});
    spectraGen.add("helium/histPShiftVsEtaHe", "PReco-PGen vs #eta", HistType::kTH2F, {{200, -2.f, 2.f}, {400, -4.f, 4.f}});

    spectraGen.add("helium/histPtGenantiHe", "PtGenantiHe", HistType::kTH1F, {{800, 0.f, 8.f}});
    spectraGen.add("helium/histPtRecantiHe", "PtRecantiHe", HistType::kTH1F, {{800, 0.f, 8.f}});
    spectraGen.add("helium/histPtShiftantiHe", "PtReco-PtGen vs PtReco", HistType::kTH2F, {{800, 0.f, 8.f}, {400, -4.f, 4.f}});
    spectraGen.add("helium/histPtShiftVsEtaantiHe", "PtReco-PtGen vs #eta", HistType::kTH2F, {{200, -2.f, 2.f}, {400, -4.f, 4.f}});

    spectraGen.add("helium/histPGenantiHe", "PGenantiHe", HistType::kTH1F, {{800, 0.f, 8.f}});
    spectraGen.add("helium/histPRecantiHe", "PRecantiHe", HistType::kTH1F, {{800, 0.f, 8.f}});
    spectraGen.add("helium/histPShiftantiHe", "PReco-PGen vs PReco", HistType::kTH2F, {{800, 0.f, 8.f}, {400, -4.f, 4.f}});
    spectraGen.add("helium/histPShiftVsEtaantiHe", "PReco-PGen vs #eta", HistType::kTH2F, {{200, -2.f, 2.f}, {400, -4.f, 4.f}});

    spectraGen.add("pion/histGenPtPion", "generated particles", HistType::kTH1F, {ptAxis});
    spectraGen.add("pion/histGenPtPionPrim", "generated particles", HistType::kTH1F, {ptAxis});
    spectraGen.add("pion/histGenPtPionSec", "generated particles", HistType::kTH1F, {ptAxis});
    spectraGen.add("pion/histSecTransportPtPion", "generated particles", HistType::kTH1F, {ptAxis});

    spectraGen.add("pion/histGenPtantiPion", "generated particles", HistType::kTH1F, {ptAxis});
    spectraGen.add("pion/histGenPtantiPionPrim", "generated particles", HistType::kTH1F, {ptAxis});
    spectraGen.add("pion/histGenPtantiPionSec", "generated particles", HistType::kTH1F, {ptAxis});
    spectraGen.add("pion/histSecTransportPtantiPion", "generated particles", HistType::kTH1F, {ptAxis});

    spectraGen.add("kaon/histGenPtKaon", "generated particles", HistType::kTH1F, {ptAxis});
    spectraGen.add("kaon/histGenPtKaonPrim", "generated particles", HistType::kTH1F, {ptAxis});
    spectraGen.add("kaon/histGenPtKaonSec", "generated particles", HistType::kTH1F, {ptAxis});
    spectraGen.add("kaon/histSecTransportPtKaon", "generated particles", HistType::kTH1F, {ptAxis});

    spectraGen.add("kaon/histGenPtantiKaon", "generated particles", HistType::kTH1F, {ptAxis});
    spectraGen.add("kaon/histGenPtantiKaonPrim", "generated particles", HistType::kTH1F, {ptAxis});
    spectraGen.add("kaon/histGenPtantiKaonSec", "generated particles", HistType::kTH1F, {ptAxis});
    spectraGen.add("kaon/histSecTransportPtantiKaon", "generated particles", HistType::kTH1F, {ptAxis});

    if (enablePr) {
      spectraGen.add("proton/histGenPtProton", "generated particles", HistType::kTH1F, {ptAxis});
      spectraGen.add("proton/histGenPtProtonPrim", "generated particles", HistType::kTH1F, {ptAxis});
      spectraGen.add("proton/histGenPtProtonSec", "generated particles", HistType::kTH1F, {ptAxis});
      spectraGen.add("proton/histSecTransportPtProton", "generated particles", HistType::kTH1F, {ptAxis});

      spectraGen.add("proton/histGenPtantiProton", "generated particles", HistType::kTH1F, {ptAxis});
      spectraGen.add("proton/histGenPtantiProtonPrim", "generated particles", HistType::kTH1F, {ptAxis});
      spectraGen.add("proton/histGenPtantiProtonSec", "generated particles", HistType::kTH1F, {ptAxis});
      spectraGen.add("proton/histSecTransportPtantiProton", "generated particles", HistType::kTH1F, {ptAxis});

      spectraGen.add("proton/histGenPtProtonPrim_Y", "generated particles", HistType::kTH2F, {{150, -1.5f, 1.5f}, {ptAxis}});
      spectraGen.add("proton/histGenPtantiProtonPrim_Y", "generated particles", HistType::kTH2F, {{150, -1.5f, 1.5f}, {ptAxis}});
    }
    if (enableDe) {
      spectraGen.add("deuteron/histGenPtD", "generated particles", HistType::kTH1F, {ptAxis});
      spectraGen.add("deuteron/histGenPtDPrim", "generated particles", HistType::kTH1F, {ptAxis});
      spectraGen.add("deuteron/histGenPtDSec", "generated particles", HistType::kTH1F, {ptAxis});
      spectraGen.add("deuteron/histSecTransportPtD", "generated particles", HistType::kTH1F, {ptAxis});
      histos.add("tracks/deuteron/histAntiDeuteronPtShiftRec", "histAntiDeuteronPtShiftRec", HistType::kTH1F, {ptAxis});
      histos.add("tracks/deuteron/histAntiDeuteronPtRec", "histAntiDeuteronPtRec", HistType::kTH1F, {ptAxis});
      spectraGen.add("deuteron/histAntiDeuteronPtShiftCorrection", "histAntiDeuteronPtShiftCorrection", HistType::kTH2F, {{800, 0.f, 8.f}, {400, -4.f, 4.f}});
      spectraGen.add("deuteron/histAntiDeuteronPtShift", "PtReco-PtGen vs PtReco", HistType::kTH2F, {{800, 0.f, 8.f}, {400, -4.f, 4.f}});
      spectraGen.add("deuteron/histAntiDeuteronPtShiftVsEta", "PtReco-PtGen vs #eta", HistType::kTH2F, {{200, -2.f, 2.f}, {400, -4.f, 4.f}});

      histos.add("tracks/deuteron/histDeuteronPtShiftRec", "histDeuteronPtShiftRec", HistType::kTH1F, {ptAxis});
      histos.add("tracks/deuteron/histDeuteronPtRec", "histDeuteronPtRec", HistType::kTH1F, {ptAxis});
      spectraGen.add("deuteron/histDeuteronPtShiftCorrection", "histDeuteronPtShiftCorrection", HistType::kTH2F, {{800, 0.f, 8.f}, {400, -4.f, 4.f}});
      spectraGen.add("deuteron/histDeuteronPtShift", "PtReco-PtGen vs PtReco", HistType::kTH2F, {{800, 0.f, 8.f}, {400, -4.f, 4.f}});
      spectraGen.add("deuteron/histDeuteronPtShiftVsEta", "PtReco-PtGen vs #eta", HistType::kTH2F, {{200, -2.f, 2.f}, {400, -4.f, 4.f}});

      spectraGen.add("deuteron/histGenPtantiD", "generated particles", HistType::kTH1F, {ptAxis});
      spectraGen.add("deuteron/histGenPtantiDPrim", "generated particles", HistType::kTH1F, {ptAxis});
      spectraGen.add("deuteron/histGenPtantiDSec", "generated particles", HistType::kTH1F, {ptAxis});
      spectraGen.add("deuteron/histSecTransportPtantiD", "generated particles", HistType::kTH1F, {ptAxis});
    }
    if (enableTr) {
      spectraGen.add("triton/histGenPtT", "generated particles", HistType::kTH1F, {ptAxis});
      spectraGen.add("triton/histGenPtTPrim", "generated particles", HistType::kTH1F, {ptAxis});
      spectraGen.add("triton/histGenPtTSec", "generated particles", HistType::kTH1F, {ptAxis});
      spectraGen.add("triton/histSecTransportPtT", "generated particles", HistType::kTH1F, {ptAxis});

      spectraGen.add("triton/histGenPtantiT", "generated particles", HistType::kTH1F, {ptAxis});
      spectraGen.add("triton/histGenPtantiTPrim", "generated particles", HistType::kTH1F, {ptAxis});
      spectraGen.add("triton/histGenPtantiTSec", "generated particles", HistType::kTH1F, {ptAxis});
      spectraGen.add("triton/histSecTransportPtantiT", "generated particles", HistType::kTH1F, {ptAxis});
    }
    if (enableHe) {
      spectraGen.add("helium/histGenPtHe", "generated particles", HistType::kTH1F, {ptHeAxis});
      spectraGen.add("helium/histGenPtHePrim", "generated particles", HistType::kTH1F, {ptHeAxis});
      if (enableCentrality)
        spectraGen.add("helium/histGenPtHePrimVsMult", "generated particles", HistType::kTH2F, {{ptHeAxis}, {binsPercentile}});
      spectraGen.add("helium/histGenPtHeSec", "generated particles", HistType::kTH1F, {ptHeAxis});
      spectraGen.add("helium/histSecTransportPtHe", "generated particles", HistType::kTH1F, {ptHeAxis});

      spectraGen.add("helium/histGenPtantiHe", "generated particles", HistType::kTH1F, {ptHeAxis});
      spectraGen.add("helium/histGenPtantiHePrim", "generated particles", HistType::kTH1F, {ptHeAxis});
      if (enableCentrality)
        spectraGen.add("helium/histGenPtantiHePrimVsMult", "generated particles", HistType::kTH2F, {{ptHeAxis}, {binsPercentile}});
      spectraGen.add("helium/histGenPtantiHeSec", "generated particles", HistType::kTH1F, {ptHeAxis});
      spectraGen.add("helium/histSecTransportPtantiHe", "generated particles", HistType::kTH1F, {ptHeAxis});
    }
    if (enableAl) {
      spectraGen.add("alpha/histGenPtAl", "generated particles", HistType::kTH1F, {ptAxis});
      spectraGen.add("alpha/histGenPtAlPrim", "generated particles", HistType::kTH1F, {ptAxis});
      spectraGen.add("alpha/histGenPtAlSec", "generated particles", HistType::kTH1F, {ptAxis});
      spectraGen.add("alpha/histSecTransportPtAl", "generated particles", HistType::kTH1F, {ptAxis});

      spectraGen.add("alpha/histGenPtantiAl", "generated particles", HistType::kTH1F, {ptAxis});
      spectraGen.add("alpha/histGenPtantiAlPrim", "generated particles", HistType::kTH1F, {ptAxis});
      spectraGen.add("alpha/histGenPtantiAlSec", "generated particles", HistType::kTH1F, {ptAxis});
      spectraGen.add("alpha/histSecTransportPtantiAl", "generated particles", HistType::kTH1F, {ptAxis});
    }
    LOG(info) << "Histograms of LFNucleiBATask:";
    histos.print();
    if (doprocessMCGen)
      spectraGen.print();
  }

  template <bool IsMC, bool IsFilteredData, typename CollisionType, typename TracksType, typename ParticleType>
  void fillHistograms(const CollisionType& event,
                      const TracksType& tracks,
                      const ParticleType& particles)
  {
    histos.fill(HIST("event/eventSkimming"), 0.5);
    // Apply skimming
    if constexpr (!IsFilteredData) {
      const auto& bc = event.template bc_as<o2::aod::BCsWithTimestamps>();
      initCCDB(bc);
      if (skimmingOptions.applySkimming) {
        if (!zorro.isSelected(bc.globalBC())) {
          return;
        }
      }
      histos.fill(HIST("event/eventSkimming"), 1.5);
    }

    // Event histos fill
    if (enableCentrality)
      histos.fill(HIST("event/eventSelection"), 0, event.centFT0M());
    else
      histos.fill(HIST("event/eventSelection"), 0);
    if (enableDebug)
      debugHistos.fill(HIST("qa/h1VtxZ_nocut"), event.posZ());

    if constexpr (!IsFilteredData) {
      if (!event.selection_bit(aod::evsel::kIsTriggerTVX)) {
        if (evselOptions.useTVXtrigger)
          return;
      } else {
        if (enableCentrality)
          histos.fill(HIST("event/eventSelection"), 1, event.centFT0M());
        else
          histos.fill(HIST("event/eventSelection"), 1);
        if (enableDebug)
          debugHistos.fill(HIST("qa/h1VtxZ_TVXtrigger"), event.posZ());
      }

      if (!event.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
        if (evselOptions.removeTFBorder)
          return;
      } else {
        if (enableCentrality)
          histos.fill(HIST("event/eventSelection"), 2, event.centFT0M());
        else
          histos.fill(HIST("event/eventSelection"), 2);
        if (enableDebug)
          debugHistos.fill(HIST("qa/h1VtxZ_TFrameBorder"), event.posZ());
      }

      if (!event.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
        if (evselOptions.removeITSROFBorder)
          return;
      } else {
        if (enableCentrality)
          histos.fill(HIST("event/eventSelection"), 3, event.centFT0M());
        else
          histos.fill(HIST("event/eventSelection"), 3);
        if (enableDebug)
          debugHistos.fill(HIST("qa/h1VtxZ_ITSROFBorder"), event.posZ());
      }

      if ((event.selection_bit(aod::evsel::kNoITSROFrameBorder)) &&
          (event.selection_bit(aod::evsel::kNoTimeFrameBorder)) &&
          (event.selection_bit(aod::evsel::kIsTriggerTVX))) {
        if (enableCentrality)
          histos.fill(HIST("event/eventSelection"), 4, event.centFT0M());
        else
          histos.fill(HIST("event/eventSelection"), 4);
      }

      if (evselOptions.useSel8 && !event.sel8())
        return;
      if (enableCentrality)
        histos.fill(HIST("event/eventSelection"), 5, event.centFT0M());
      else
        histos.fill(HIST("event/eventSelection"), 5);
      if (enableDebug)
        debugHistos.fill(HIST("qa/h1VtxZ_sel8"), event.posZ());

      if (event.posZ() < cfgVzCutLow || event.posZ() > cfgVzCutHigh)
        return;
      if (enableCentrality)
        histos.fill(HIST("event/eventSelection"), 6, event.centFT0M());
      else
        histos.fill(HIST("event/eventSelection"), 6);

    } else {
      if (event.posZ() < cfgVzCutLow || event.posZ() > cfgVzCutHigh)
        return;
      if (evselOptions.removeTFBorder && !event.selection_bit(aod::evsel::kNoTimeFrameBorder))
        return;
    }
    if (event.centFT0M() <= cfgMultCutLow || event.centFT0M() > cfgMultCutHigh) {
      return;
    }
    if (enableCentrality)
      histos.fill(HIST("event/eventSelection"), 7, event.centFT0M());
    else
      histos.fill(HIST("event/eventSelection"), 7);

    if (enableCentrality && enableDebug) {
      debugHistos.fill(HIST("event/h1VtxZ_Centrality"), event.posZ());
    }

    float gamma = 0., massTOF = 0., massTOFhe = 0., massTOFantihe = 0., heTPCmomentum = 0.f, antiheTPCmomentum = 0.f, heP = 0.f, antiheP = 0.f, hePt = 0.f, antihePt = 0.f, antiDPt = 0.f, DPt = 0.f;
    bool isTritonTPCpid = false;
    bool prRapCut = false;
    bool deRapCut = false;
    bool trRapCut = false;
    bool heRapCut = false;
    bool alRapCut = false;

    // Event histos fill
    if (enableCentrality)
      histos.fill(HIST("event/h1VtxZ"), event.posZ(), event.centFT0M());
    else
      histos.fill(HIST("event/h1VtxZ"), event.posZ());
    if (enableDebug && enableCentrality)
      debugHistos.fill(HIST("event/hFT0M"), event.centFT0M());

    if constexpr (IsFilteredData) {
      if (enableCentrality)
        debugHistos.fill(HIST("event/hFV0M"), event.centFV0M());
    }

    auto tracksWithITS = soa::Attach<TracksType,
                                     aod::pidits::ITSNSigmaDe,
                                     aod::pidits::ITSNSigmaTr,
                                     aod::pidits::ITSNSigmaHe>(tracks);
    if (tracksWithITS.size() != tracks.size()) {
      LOG(fatal) << "Problem with track size";
    }

    tracks.copyIndexBindings(tracksWithITS);

    for (auto const& track : tracksWithITS) {
      if constexpr (!IsFilteredData) {
        if (!track.isGlobalTrackWoDCA() && filterOptions.enableIsGlobalTrack) {
          continue;
        }
      }
      std::bitset<8> itsClusterMap = track.itsClusterMap();

      if constexpr (!IsFilteredData) {
        if (nsigmaITSvar.showAverageClusterSize && outFlagOptions.enablePIDplot)
          histos.fill(HIST("tracks/avgClusterSizePerCoslInvVsITSlayers"), track.p(), averageClusterSizePerCoslInv(track), track.itsNCls());
      }

      if (track.itsNCls() < trkqcOptions.cfgCutITSClusters ||
          track.tpcNClsCrossedRows() < trkqcOptions.cfgCutTPCXRows ||
          track.tpcNClsFound() < trkqcOptions.cfgCutTPCClusters ||
          track.tpcCrossedRowsOverFindableCls() < trkqcOptions.cfgCutTPCCROFnd) {
        continue;
      }

      auto tpcChi2NclRange = (std::vector<float>)trkqcOptions.tpcChi2NclCuts;
      if ((track.tpcChi2NCl() < tpcChi2NclRange[0]) || (track.tpcChi2NCl() > tpcChi2NclRange[1]))
        continue;
      auto itsChi2NclRange = (std::vector<float>)trkqcOptions.itsChi2NclCuts;
      if ((track.itsChi2NCl() < itsChi2NclRange[0]) || (track.itsChi2NCl() > itsChi2NclRange[1]))
        continue;

      // p & eta cut
      if (std::abs(track.tpcInnerParam()) < kinemOptions.cfgMomentumCut ||
          std::abs(track.eta()) > kinemOptions.cfgEtaCut)
        continue;

      if (outFlagOptions.enablePIDplot) {
        histos.fill(HIST("tracks/h1pT"), track.pt());
        histos.fill(HIST("tracks/h1p"), track.p());
      }

      isTritonTPCpid = std::abs(track.tpcNSigmaTr()) < nsigmaTPCvar.nsigmaTPCTr;

      float shiftPtPos = 0.f;
      float shiftPtNeg = 0.f;

      if (enablePtShiftHe && !fShiftPtHe) {
        fShiftPtHe = new TF1("fShiftPtHe", "[0] * exp([1] + [2] * x) + [3] + [4] * x", 0.f, 8.f);
        auto parHe = (std::vector<float>)parShiftPtHe; // NOLINT
        fShiftPtHe->SetParameters(parHe[0], parHe[1], parHe[2], parHe[3], parHe[4]);
      }

      if (enablePtShiftHe && !fShiftPtantiHe) {
        fShiftPtantiHe = new TF1("fShiftPtantiHe", "[0] * exp([1] + [2] * x) + [3] + [4] * x", 0.f, 8.f);
        auto parAntiHe = (std::vector<float>)parShiftPtAntiHe; // NOLINT
        fShiftPtantiHe->SetParameters(parAntiHe[0], parAntiHe[1], parAntiHe[2], parAntiHe[3], parAntiHe[4]);
      }

      if (enablePtShiftAntiD && !fShiftAntiD) {
        fShiftAntiD = new TF1("fShiftAntiD", "[0] * exp([1] + [2] * x) + [3] + [4] * x", 0.f, 8.f);
        auto parAntiD = (std::vector<float>)parShiftPtAntiD; // NOLINT
        fShiftAntiD->SetParameters(parAntiD[0], parAntiD[1], parAntiD[2], parAntiD[3], parAntiD[4]);
      }

      switch (unableAntiDPtShift) {
        case 0:
          if (enablePtShiftAntiD && fShiftAntiD) {
            auto shiftAntiD = fShiftAntiD->Eval(track.pt());
            antiDPt = track.pt() - shiftAntiD;
          }
          break;
        case 1:
          antiDPt = track.pt();
          break;
      }

      if (enablePtShiftD && !fShiftD) {
        fShiftD = new TF1("fShiftD", "[0] * exp([1] + [2] * x) + [3] + [4] * x", 0.f, 8.f);
        auto parD = (std::vector<float>)parShiftPtD; // NOLINT
        fShiftD->SetParameters(parD[0], parD[1], parD[2], parD[3], parD[4]);
      }

      switch (unableDPtShift) {
        case 0:
          if (enablePtShiftD && fShiftD) {
            auto shiftD = fShiftD->Eval(track.pt());
            DPt = track.pt() - shiftD;
          }
          break;
        case 1:
          DPt = track.pt();
          break;
      }

      switch (helium3Pt) {
        case 0:
          hePt = track.pt();
          if (enablePtShiftHe && fShiftPtHe) {
            shiftPtPos = fShiftPtHe->Eval(2 * track.pt());
            hePt = track.pt() - shiftPtPos / 2.f;
          }
          antihePt = track.pt();
          if (enablePtShiftHe && fShiftPtantiHe) {
            shiftPtNeg = fShiftPtantiHe->Eval(2 * track.pt());
            antihePt = track.pt() - shiftPtNeg / 2.f;
          }
          break;
        case 1:
          hePt = 2.f * track.pt();
          antihePt = 2.f * track.pt();
          break;
      }

      float nITSDe = 99.f;
      float nITSHe = 99.f;

      if (!IsFilteredData) {
        nITSDe = track.itsNSigmaDe();
        nITSHe = track.itsNSigmaHe();
      }
      heP = track.p();
      antiheP = track.p();
      heTPCmomentum = track.tpcInnerParam();
      antiheTPCmomentum = track.tpcInnerParam();

      // auto parDCAxy = (std::vector<float>)parDCAxycuts;
      // auto parDCAz = (std::vector<float>)parDCAzcuts;
      const auto& parDCAxy = parDCAxycuts.value;
      const auto& parDCAz = parDCAzcuts.value;

      bool passDCAxyCut = false;
      bool passDCAzCut = false;
      bool passDCAxyCutDe = false;
      bool passDCAzCutDe = false;
      bool passDCAxyCutAntiDe = false;
      bool passDCAzCutAntiDe = false;
      bool passDCAxyCutHe = false;
      bool passDCAzCutHe = false;
      bool passDCAxyCutAntiHe = false;
      bool passDCAzCutAntiHe = false;

      bool isDeuteron = false;
      bool isHelium = false;
      bool isDe = false;
      bool isAntiDe = false;
      bool isHe = false;
      bool isAntiHe = false;

      bool isDeWoDCAxy = false;
      bool isAntiDeWoDCAxy = false;
      bool isHeWoDCAxy = false;
      bool isAntiHeWoDCAxy = false;

      bool isDeWoDCAz = false;
      bool isAntiDeWoDCAz = false;
      bool isHeWoDCAz = false;
      bool isAntiHeWoDCAz = false;

      bool isDeWoDCAxyWTPCpid = false;
      bool isAntiDeWoDCAxyWTPCpid = false;
      bool isHeWoDCAxyWTPCpid = false;
      bool isAntiHeWoDCAxyWTPCpid = false;

      bool isDeWoDCAzWTPCpid = false;
      bool isAntiDeWoDCAzWTPCpid = false;
      bool isHeWoDCAzWTPCpid = false;
      bool isAntiHeWoDCAzWTPCpid = false;

      bool isDeWoTPCpid = false;
      bool isAntiDeWoTPCpid = false;
      bool isHeWoTPCpid = false;
      bool isAntiHeWoTPCpid = false;

      bool isDeWTPCpid = false;
      bool isAntiDeWTPCpid = false;
      bool isHeWTPCpid = false;
      bool isAntiHeWTPCpid = false;

      bool passDCAxyzCut = false;

      const float dcaXY2 = track.dcaXY() * track.dcaXY();
      const float dcaZ2 = track.dcaZ() * track.dcaZ();

      switch (dcaConfOptions.cfgCustomDCA) {
        case 0:
          passDCAxyCut = passDCAxyCutDe = passDCAxyCutAntiDe = passDCAxyCutHe = passDCAxyCutAntiHe = (std::abs(track.dcaXY()) <= dcaConfOptions.cfgCustomDCAxy);
          passDCAzCut = passDCAzCutDe = passDCAzCutAntiDe = passDCAzCutHe = passDCAzCutAntiHe = (std::abs(track.dcaZ()) <= dcaConfOptions.cfgCustomDCAz);
          break;
        case 1:
          passDCAxyCut = (std::abs(track.dcaXY()) <= parDCAxy[3] * (parDCAxy[0] + parDCAxy[1] / std::pow(track.pt(), parDCAxy[2])));
          passDCAzCut = (std::abs(track.dcaZ()) <= parDCAz[3] * (parDCAz[0] + parDCAz[1] / std::pow(track.pt(), parDCAz[2])));

          passDCAxyCutDe = (std::abs(track.dcaXY()) <= parDCAxy[3] * (parDCAxy[0] + parDCAxy[1] / std::pow(DPt, parDCAxy[2])));
          passDCAzCutDe = (std::abs(track.dcaZ()) <= parDCAz[3] * (parDCAz[0] + parDCAz[1] / std::pow(DPt, parDCAz[2])));
          passDCAxyCutAntiDe = (std::abs(track.dcaXY()) <= parDCAxy[3] * (parDCAxy[0] + parDCAxy[1] / std::pow(antiDPt, parDCAxy[2])));
          passDCAzCutAntiDe = (std::abs(track.dcaZ()) <= parDCAz[3] * (parDCAz[0] + parDCAz[1] / std::pow(antiDPt, parDCAz[2])));

          passDCAxyCutHe = (std::abs(track.dcaXY()) <= parDCAxy[3] * (parDCAxy[0] + parDCAxy[1] / std::pow(hePt, parDCAxy[2])));
          passDCAzCutHe = (std::abs(track.dcaZ()) <= parDCAz[3] * (parDCAz[0] + parDCAz[1] / std::pow(hePt, parDCAz[2])));
          passDCAxyCutAntiHe = (std::abs(track.dcaXY()) <= parDCAxy[3] * (parDCAxy[0] + parDCAxy[1] / std::pow(antihePt, parDCAxy[2])));
          passDCAzCutAntiHe = (std::abs(track.dcaZ()) <= parDCAz[3] * (parDCAz[0] + parDCAz[1] / std::pow(antihePt, parDCAz[2])));
          break;
        case 2:
          passDCAxyCutDe = (std::abs(track.dcaXY()) <= parDCAxy[3] * (parDCAxy[0] + parDCAxy[1] / std::pow(DPt, parDCAxy[2])));
          passDCAxyCutAntiDe = (std::abs(track.dcaXY()) <= parDCAxy[3] * (parDCAxy[0] + parDCAxy[1] / std::pow(antiDPt, parDCAxy[2])));

          passDCAxyCutHe = (std::abs(track.dcaXY()) <= parDCAxy[3] * (parDCAxy[0] + parDCAxy[1] / std::pow(hePt, parDCAxy[2])));
          passDCAxyCutAntiHe = (std::abs(track.dcaXY()) <= parDCAxy[3] * (parDCAxy[0] + parDCAxy[1] / std::pow(antihePt, parDCAxy[2])));

          passDCAxyCut = (std::abs(track.dcaXY()) <= parDCAxy[3] * (parDCAxy[0] + parDCAxy[1] / std::pow(track.pt(), parDCAxy[2])));
          passDCAzCut = passDCAzCutDe = passDCAzCutAntiDe = passDCAzCutHe = passDCAzCutAntiHe = (std::abs(track.dcaZ()) <= dcaConfOptions.cfgCustomDCAz);
          break;
        case 3:
          passDCAxyCut = passDCAxyCutDe = passDCAxyCutAntiDe = passDCAxyCutHe = passDCAxyCutAntiHe = (std::abs(track.dcaXY()) <= dcaConfOptions.cfgCustomDCAxy);
          passDCAzCut = (std::abs(track.dcaZ()) <= parDCAz[3] * (parDCAz[0] + parDCAz[1] / std::pow(track.pt(), parDCAz[2])));

          passDCAzCutDe = (std::abs(track.dcaZ()) <= parDCAz[3] * (parDCAz[0] + parDCAz[1] / std::pow(DPt, parDCAz[2])));
          passDCAzCutAntiDe = (std::abs(track.dcaZ()) <= parDCAz[3] * (parDCAz[0] + parDCAz[1] / std::pow(antiDPt, parDCAz[2])));

          passDCAzCutHe = (std::abs(track.dcaZ()) <= parDCAz[3] * (parDCAz[0] + parDCAz[1] / std::pow(hePt, parDCAz[2])));
          passDCAzCutAntiHe = (std::abs(track.dcaZ()) <= parDCAz[3] * (parDCAz[0] + parDCAz[1] / std::pow(antihePt, parDCAz[2])));
          break;
        case 4:
          passDCAxyCut = passDCAzCut = passDCAxyCutDe = passDCAzCutDe = passDCAxyCutAntiDe = passDCAzCutAntiDe = passDCAxyCutHe = passDCAzCutHe = passDCAxyCutAntiHe = passDCAzCutAntiHe = dcaXY2 / std::pow(dcaConfOptions.cfgCustomDCAxy, 2) + dcaZ2 / std::pow(dcaConfOptions.cfgCustomDCAz, 2) <= 1;
          break;
        case 5:
          passDCAxyCut = dcaXY2 / std::pow(parDCAxy[3] * (parDCAxy[0] + parDCAxy[1] / std::pow(track.pt(), parDCAxy[2])), 2) + dcaZ2 / std::pow(parDCAz[3] * (parDCAz[0] + parDCAz[1] / std::pow(track.pt(), parDCAz[2])), 2) <= 1;
          passDCAzCut = dcaXY2 / std::pow(parDCAxy[3] * (parDCAxy[0] + parDCAxy[1] / std::pow(track.pt(), parDCAxy[2])), 2) + dcaZ2 / std::pow(parDCAz[3] * (parDCAz[0] + parDCAz[1] / std::pow(track.pt(), parDCAz[2])), 2) <= 1;

          passDCAxyCutDe = dcaXY2 / std::pow(parDCAxy[3] * (parDCAxy[0] + parDCAxy[1] / std::pow(DPt, parDCAxy[2])), 2) + dcaZ2 / std::pow(parDCAz[3] * (parDCAz[0] + parDCAz[1] / std::pow(DPt, parDCAz[2])), 2) <= 1;
          passDCAzCutDe = dcaXY2 / std::pow(parDCAxy[3] * (parDCAxy[0] + parDCAxy[1] / std::pow(DPt, parDCAxy[2])), 2) + dcaZ2 / std::pow(parDCAz[3] * (parDCAz[0] + parDCAz[1] / std::pow(DPt, parDCAz[2])), 2) <= 1;
          passDCAxyCutAntiDe = dcaXY2 / std::pow(parDCAxy[3] * (parDCAxy[0] + parDCAxy[1] / std::pow(antiDPt, parDCAxy[2])), 2) + dcaZ2 / std::pow(parDCAz[3] * (parDCAz[0] + parDCAz[1] / std::pow(antiDPt, parDCAz[2])), 2) <= 1;
          passDCAzCutAntiDe = dcaXY2 / std::pow(parDCAxy[3] * (parDCAxy[0] + parDCAxy[1] / std::pow(antiDPt, parDCAxy[2])), 2) + dcaZ2 / std::pow(parDCAz[3] * (parDCAz[0] + parDCAz[1] / std::pow(antiDPt, parDCAz[2])), 2) <= 1;

          passDCAxyCutHe = dcaXY2 / std::pow(parDCAxy[3] * (parDCAxy[0] + parDCAxy[1] / std::pow(hePt, parDCAxy[2])), 2) + dcaZ2 / std::pow(parDCAz[3] * (parDCAz[0] + parDCAz[1] / std::pow(hePt, parDCAz[2])), 2) <= 1;
          passDCAzCutHe = dcaXY2 / std::pow(parDCAxy[3] * (parDCAxy[0] + parDCAxy[1] / std::pow(hePt, parDCAxy[2])), 2) + dcaZ2 / std::pow(parDCAz[3] * (parDCAz[0] + parDCAz[1] / std::pow(hePt, parDCAz[2])), 2) <= 1;
          passDCAxyCutAntiHe = dcaXY2 / std::pow(parDCAxy[3] * (parDCAxy[0] + parDCAxy[1] / std::pow(antihePt, parDCAxy[2])), 2) + dcaZ2 / std::pow(parDCAz[3] * (parDCAz[0] + parDCAz[1] / std::pow(antihePt, parDCAz[2])), 2) <= 1;
          passDCAzCutAntiHe = dcaXY2 / std::pow(parDCAxy[3] * (parDCAxy[0] + parDCAxy[1] / std::pow(antihePt, parDCAxy[2])), 2) + dcaZ2 / std::pow(parDCAz[3] * (parDCAz[0] + parDCAz[1] / std::pow(antihePt, parDCAz[2])), 2) <= 1;
          break;
      }

      // Rapidity cuts
      auto rapCheck = [&](float m2z) {
        const float rap = track.rapidity(m2z);
        return (rap > kinemOptions.cfgRapidityCutLow) && (rap < kinemOptions.cfgRapidityCutHigh);
      };

      prRapCut = rapCheck(MassProtonVal);
      deRapCut = rapCheck(MassDeuteronVal);
      trRapCut = rapCheck(MassTritonVal);
      heRapCut = rapCheck(MassHeliumVal / 2.0);
      alRapCut = rapCheck(MassAlphaVal / 2.0);

      isDeuteron = enableDe && deRapCut;
      isHelium = enableHe && heRapCut;
      isDe = isDeuteron && track.sign() > 0;
      isAntiDe = isDeuteron && track.sign() < 0;

      if constexpr (IsMC && !IsFilteredData) {
        int pdgCheck = track.mcParticle().pdgCode();
        if (std::abs(pdgCheck) == PDGDeuteron)
          histos.fill(HIST("tracks/hItsDeHeChecker"), 0);
        if (std::abs(pdgCheck) == PDGHelium)
          histos.fill(HIST("tracks/hItsDeHeChecker"), 1);
      }

      // nSigmaITSHe cut
      if (nsigmaITSvar.useITSDeCut && (nITSDe <= nsigmaITSvar.nsigmaITSDe)) {
        continue;
      }

      if (nsigmaITSvar.useITSHeCut && (nITSHe <= nsigmaITSvar.nsigmaITSHe)) {
        continue;
      }

      if constexpr (IsMC && !IsFilteredData) {
        int pdgCheck = track.mcParticle().pdgCode();
        if (std::abs(pdgCheck) == PDGDeuteron)
          histos.fill(HIST("tracks/hItsDeHeChecker"), 2);
        if (std::abs(pdgCheck) == PDGHelium)
          histos.fill(HIST("tracks/hItsDeHeChecker"), 3);
      }

      isHe = isHelium && track.sign() > 0;
      isAntiHe = isHelium && track.sign() < 0;

      isDeWoDCAxy = isDe && passDCAzCutDe;
      isAntiDeWoDCAxy = isAntiDe && passDCAzCutAntiDe;
      isHeWoDCAxy = isHe && passDCAzCutHe;
      isAntiHeWoDCAxy = isAntiHe && passDCAzCutAntiHe;

      isDeWoDCAz = isDe && passDCAxyCutDe;
      isAntiDeWoDCAz = isAntiDe && passDCAxyCutAntiDe;
      isHeWoDCAz = isHe && passDCAxyCutHe;
      isAntiHeWoDCAz = isAntiHe && passDCAxyCutAntiHe;

      isDeWoDCAxyWTPCpid = isDeWoDCAxy && std::abs(track.tpcNSigmaDe()) < nsigmaTPCvar.nsigmaTPCDe;
      isAntiDeWoDCAxyWTPCpid = isAntiDeWoDCAxy && std::abs(track.tpcNSigmaDe()) < nsigmaTPCvar.nsigmaTPCDe;
      isHeWoDCAxyWTPCpid = isHeWoDCAxy && std::abs(track.tpcNSigmaHe()) < nsigmaTPCvar.nsigmaTPCHe;
      isAntiHeWoDCAxyWTPCpid = isAntiHeWoDCAxy && std::abs(track.tpcNSigmaHe()) < nsigmaTPCvar.nsigmaTPCHe;

      isDeWoDCAzWTPCpid = isDeWoDCAz && std::abs(track.tpcNSigmaDe()) < nsigmaTPCvar.nsigmaTPCDe;
      isAntiDeWoDCAzWTPCpid = isAntiDeWoDCAz && std::abs(track.tpcNSigmaDe()) < nsigmaTPCvar.nsigmaTPCDe;
      isHeWoDCAzWTPCpid = isHeWoDCAz && std::abs(track.tpcNSigmaHe()) < nsigmaTPCvar.nsigmaTPCHe;
      isAntiHeWoDCAzWTPCpid = isAntiHeWoDCAz && std::abs(track.tpcNSigmaHe()) < nsigmaTPCvar.nsigmaTPCHe;

      isDeWoTPCpid = isDe && passDCAzCutDe && passDCAxyCutDe;
      isAntiDeWoTPCpid = isAntiDe && passDCAzCutAntiDe && passDCAxyCutAntiDe;
      isHeWoTPCpid = isHe && passDCAzCutHe && passDCAxyCutHe;
      isAntiHeWoTPCpid = isAntiHe && passDCAzCutAntiHe && passDCAxyCutAntiHe;

      isDeWTPCpid = isDeWoTPCpid && std::abs(track.tpcNSigmaDe()) < nsigmaTPCvar.nsigmaTPCDe;
      isAntiDeWTPCpid = isAntiDeWoTPCpid && std::abs(track.tpcNSigmaDe()) < nsigmaTPCvar.nsigmaTPCDe;
      isHeWTPCpid = isHeWoTPCpid && std::abs(track.tpcNSigmaHe()) < nsigmaTPCvar.nsigmaTPCHe;
      isAntiHeWTPCpid = isAntiHeWoTPCpid && std::abs(track.tpcNSigmaHe()) < nsigmaTPCvar.nsigmaTPCHe;

      passDCAxyzCut = passDCAxyCut && passDCAzCut;

      // DCAxy vs DCAz plots BEFORE cut
      if (outFlagOptions.makeDCABeforeCutPlots) {
        histos.fill(HIST("tracks/dca/before/hDCAxyVsDCAzVsPt"), track.dcaXY(), track.dcaZ(), track.pt());
        histos.fill(HIST("tracks/dca/before/hDCAxyVsDCAz"), track.dcaZ(), track.dcaXY());

        if (isHe && std::abs(track.tpcNSigmaHe()) < nsigmaTPCvar.nsigmaTPCHe) {
          histos.fill(HIST("tracks/helium/dca/before/h3DCAvsPtHelium"), track.dcaXY(), track.dcaZ(), hePt);
          if (track.hasTOF() && outFlagOptions.doTOFplots) {
            histos.fill(HIST("tracks/helium/dca/before/TOF/hDCAxyVsDCAzVsPtHelium"), track.dcaXY(), track.dcaZ(), hePt);
          }
        }
        if (isAntiHe && std::abs(track.tpcNSigmaHe()) < nsigmaTPCvar.nsigmaTPCHe) {
          histos.fill(HIST("tracks/helium/dca/before/h3DCAvsPtantiHelium"), track.dcaXY(), track.dcaZ(), antihePt);
          if (track.hasTOF() && outFlagOptions.doTOFplots) {
            histos.fill(HIST("tracks/helium/dca/before/TOF/hDCAxyVsDCAzVsPtantiHelium"), track.dcaXY(), track.dcaZ(), antihePt);
          }
        }
        if (passDCAxyCut) {
          histos.fill(HIST("tracks/dca/before/hDCAzVsPt"), track.pt(), track.dcaZ());

          if (enablePr && prRapCut && (std::abs(track.tpcNSigmaPr()) < nsigmaTPCvar.nsigmaTPCPr)) {
            if (track.sign() > 0) {
              histos.fill(HIST("tracks/proton/dca/before/hDCAzVsPtProton"), track.pt(), track.dcaZ());
            } else {
              histos.fill(HIST("tracks/proton/dca/before/hDCAzVsPtantiProton"), track.pt(), track.dcaZ());
            }
          }
          if (enableTr && trRapCut && isTritonTPCpid) {
            if (track.sign() > 0) {
              histos.fill(HIST("tracks/triton/dca/before/hDCAzVsPtTriton"), track.pt(), track.dcaZ());
            } else {
              histos.fill(HIST("tracks/triton/dca/before/hDCAzVsPtantiTriton"), track.pt(), track.dcaZ());
            }
          }
          if (enableAl && alRapCut && (std::abs(track.tpcNSigmaAl()) < nsigmaTPCvar.nsigmaTPCAl)) {
            if (track.sign() > 0) {
              histos.fill(HIST("tracks/alpha/dca/before/hDCAzVsPtAlpha"), track.pt(), track.dcaZ());
            } else {
              histos.fill(HIST("tracks/alpha/dca/before/hDCAzVsPtantiAlpha"), track.pt(), track.dcaZ());
            }
          }
        }

        if (isDeWoDCAzWTPCpid) {
          histos.fill(HIST("tracks/deuteron/dca/before/hDCAzVsPtDeuteron"), DPt, track.dcaZ());
          if (!track.hasTOF() && (outFlagOptions.enableNoTOFPlots))
            histos.fill(HIST("tracks/deuteron/dca/before/hDCAzVsPtDeuteronNoTOF"), DPt, track.dcaZ());
        }

        if (isAntiDeWoDCAzWTPCpid) {
          histos.fill(HIST("tracks/deuteron/dca/before/hDCAzVsPtantiDeuteron"), antiDPt, track.dcaZ());
          if (!track.hasTOF() && (outFlagOptions.enableNoTOFPlots))
            histos.fill(HIST("tracks/deuteron/dca/before/hDCAzVsPtantiDeuteronNoTOF"), antiDPt, track.dcaZ());
        }

        if (isHeWoDCAzWTPCpid) {
          histos.fill(HIST("tracks/helium/dca/before/hDCAzVsPtHelium"), hePt, track.dcaZ());
          if (!track.hasTOF() && (outFlagOptions.enableNoTOFPlots))
            histos.fill(HIST("tracks/helium/dca/before/hDCAzVsPtHeliumNoTOF"), hePt, track.dcaZ());
          if (track.hasTOF() && outFlagOptions.doTOFplots) {
            histos.fill(HIST("tracks/helium/dca/before/TOF/hDCAzVsPtHelium"), hePt, track.dcaZ());
          }
        }

        if (isAntiHeWoDCAzWTPCpid) {
          histos.fill(HIST("tracks/helium/dca/before/hDCAzVsPtantiHelium"), antihePt, track.dcaZ());
          if (!track.hasTOF() && (outFlagOptions.enableNoTOFPlots))
            histos.fill(HIST("tracks/helium/dca/before/hDCAzVsPtantiHeliumNoTOF"), hePt, track.dcaZ());
          if (track.hasTOF() && outFlagOptions.doTOFplots) {
            histos.fill(HIST("tracks/helium/dca/before/TOF/hDCAzVsPtantiHelium"), antihePt, track.dcaZ());
          }
        }
      }

      if constexpr (IsMC) {
        bool isPhysPrim = false;
        bool isProdByGen = false;
        bool isWeakDecay = false;

        // PID
        int pdgCode = 0;
        if constexpr (IsFilteredData) {
          isPhysPrim = track.isPhysicalPrimary();
          isProdByGen = track.producedByGenerator();
          isWeakDecay = (track.getProcess() == TMCProcess::kPDecay);
          pdgCode = track.pdgCode();
        } else {
          if (!track.has_mcParticle()) {
            continue;
          }
          isPhysPrim = track.mcParticle().isPhysicalPrimary();
          isProdByGen = track.mcParticle().producedByGenerator();
          isWeakDecay = (track.mcParticle().getProcess() == TMCProcess::kPDecay);
          pdgCode = track.mcParticle().pdgCode();
        }

        if (outFlagOptions.makeDCABeforeCutPlots) {
          switch (pdgCode) {
            case PDGProton:
              if (enablePr && prRapCut && passDCAxyCut) {
                histos.fill(HIST("tracks/proton/dca/before/hDCAzVsPtProtonTrue"), track.pt(), track.dcaZ());
                if (track.hasTOF() && outFlagOptions.doTOFplots) {
                  histos.fill(HIST("tracks/proton/dca/before/TOF/hDCAzVsPtProtonTrue"), track.pt(), track.dcaZ());
                }
                if (isPhysPrim) {
                  histos.fill(HIST("tracks/proton/dca/before/hDCAzVsPtProtonTruePrim"), track.pt(), track.dcaZ());
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    histos.fill(HIST("tracks/proton/dca/before/TOF/hDCAzVsPtProtonTruePrim"), track.pt(), track.dcaZ());
                  }
                }
                if (!isPhysPrim && !isProdByGen) {
                  if (isWeakDecay) {
                    histos.fill(HIST("tracks/proton/dca/before/hDCAzVsPtProtonTrueSec"), track.pt(), track.dcaZ());
                  } else {
                    histos.fill(HIST("tracks/proton/dca/before/hDCAzVsPtProtonTrueMaterial"), track.pt(), track.dcaZ());
                  }
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    if (isWeakDecay) {
                      histos.fill(HIST("tracks/proton/dca/before/TOF/hDCAzVsPtProtonTrueSec"), track.pt(), track.dcaZ());
                    } else {
                      histos.fill(HIST("tracks/proton/dca/before/TOF/hDCAzVsPtProtonTrueMaterial"), track.pt(), track.dcaZ());
                    }
                  }
                }
              }
              break;
            case -PDGProton:
              if (enablePr && prRapCut && passDCAxyCut) {
                histos.fill(HIST("tracks/proton/dca/before/hDCAzVsPtantiProtonTrue"), track.pt(), track.dcaZ());
                if (track.hasTOF() && outFlagOptions.doTOFplots) {
                  histos.fill(HIST("tracks/proton/dca/before/TOF/hDCAzVsPtantiProtonTrue"), track.pt(), track.dcaZ());
                }
                if (isPhysPrim) {
                  histos.fill(HIST("tracks/proton/dca/before/hDCAzVsPtantiProtonTruePrim"), track.pt(), track.dcaZ());
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    histos.fill(HIST("tracks/proton/dca/before/TOF/hDCAzVsPtantiProtonTruePrim"), track.pt(), track.dcaZ());
                  }
                }
                if (!isPhysPrim && !isProdByGen) {
                  if (isWeakDecay) {
                    histos.fill(HIST("tracks/proton/dca/before/hDCAzVsPtantiProtonTrueSec"), track.pt(), track.dcaZ());
                  } else {
                    histos.fill(HIST("tracks/proton/dca/before/hDCAzVsPtantiProtonTrueMaterial"), track.pt(), track.dcaZ());
                  }
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    if (isWeakDecay) {
                      histos.fill(HIST("tracks/proton/dca/before/TOF/hDCAzVsPtantiProtonTrueSec"), hePt, track.dcaZ());
                    } else {
                      histos.fill(HIST("tracks/proton/dca/before/TOF/hDCAzVsPtantiProtonTrueMaterial"), hePt, track.dcaZ());
                    }
                  }
                }
              }
              break;
            case PDGDeuteron:
              if (isDeWoDCAz) {
                histos.fill(HIST("tracks/deuteron/dca/before/hDCAzVsPtDeuteronTrue"), DPt, track.dcaZ());
                if (track.hasTOF() && outFlagOptions.doTOFplots) {
                  histos.fill(HIST("tracks/deuteron/dca/before/TOF/hDCAzVsPtDeuteronTrue"), DPt, track.dcaZ());
                }
                if (isPhysPrim) {
                  histos.fill(HIST("tracks/deuteron/dca/before/hDCAzVsPtDeuteronTruePrim"), DPt, track.dcaZ());
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    histos.fill(HIST("tracks/deuteron/dca/before/TOF/hDCAzVsPtDeuteronTruePrim"), DPt, track.dcaZ());
                  }
                }
                if (!isPhysPrim && !isProdByGen) {
                  if (isWeakDecay) {
                    histos.fill(HIST("tracks/deuteron/dca/before/hDCAzVsPtDeuteronTrueSec"), DPt, track.dcaZ());
                  } else {
                    histos.fill(HIST("tracks/deuteron/dca/before/hDCAzVsPtDeuteronTrueMaterial"), DPt, track.dcaZ());
                  }
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    if (isWeakDecay) {
                      histos.fill(HIST("tracks/deuteron/dca/before/TOF/hDCAzVsPtDeuteronTrueSec"), DPt, track.dcaZ());
                    } else {
                      histos.fill(HIST("tracks/deuteron/dca/before/TOF/hDCAzVsPtDeuteronTrueMaterial"), DPt, track.dcaZ());
                    }
                  }
                }
              }
              break;
            case -PDGDeuteron:
              if (isAntiDeWoDCAz) {
                histos.fill(HIST("tracks/deuteron/dca/before/hDCAzVsPtantiDeuteronTrue"), antiDPt, track.dcaZ());
                if (track.hasTOF() && outFlagOptions.doTOFplots) {
                  histos.fill(HIST("tracks/deuteron/dca/before/TOF/hDCAzVsPtantiDeuteronTrue"), antiDPt, track.dcaZ());
                }
                if (isPhysPrim) {
                  histos.fill(HIST("tracks/deuteron/dca/before/hDCAzVsPtantiDeuteronTruePrim"), antiDPt, track.dcaZ());
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    histos.fill(HIST("tracks/deuteron/dca/before/TOF/hDCAzVsPtantiDeuteronTruePrim"), antiDPt, track.dcaZ());
                  }
                }
                if (!isPhysPrim && !isProdByGen) {
                  if (isWeakDecay) {
                    histos.fill(HIST("tracks/deuteron/dca/before/hDCAzVsPtantiDeuteronTrueSec"), antiDPt, track.dcaZ());
                  } else {
                    histos.fill(HIST("tracks/deuteron/dca/before/hDCAzVsPtantiDeuteronTrueMaterial"), antiDPt, track.dcaZ());
                  }
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    if (isWeakDecay) {
                      histos.fill(HIST("tracks/deuteron/dca/before/TOF/hDCAzVsPtantiDeuteronTrueSec"), antiDPt, track.dcaZ());
                    } else {
                      histos.fill(HIST("tracks/deuteron/dca/before/TOF/hDCAzVsPtantiDeuteronTrueMaterial"), antiDPt, track.dcaZ());
                    }
                  }
                }
              }
              break;
            case PDGTriton:
              if (enableTr && trRapCut && passDCAxyCut) {
                histos.fill(HIST("tracks/triton/dca/before/hDCAzVsPtTritonTrue"), track.pt(), track.dcaZ());
                if (isPhysPrim) {
                  histos.fill(HIST("tracks/triton/dca/before/hDCAzVsPtTritonTruePrim"), track.pt(), track.dcaZ());
                }
                if (!isPhysPrim && !isProdByGen) {
                  if (isWeakDecay) {
                    histos.fill(HIST("tracks/triton/dca/before/hDCAzVsPtTritonTrueSec"), track.pt(), track.dcaZ());
                  } else {
                    histos.fill(HIST("tracks/triton/dca/before/hDCAzVsPtTritonTrueMaterial"), track.pt(), track.dcaZ());
                  }
                }
              }
              break;
            case -PDGTriton:
              if (enableTr && trRapCut && passDCAxyCut) {
                histos.fill(HIST("tracks/triton/dca/before/hDCAzVsPtantiTritonTrue"), track.pt(), track.dcaZ());
                if (isPhysPrim) {
                  histos.fill(HIST("tracks/triton/dca/before/hDCAzVsPtantiTritonTruePrim"), track.pt(), track.dcaZ());
                }
                if (!isPhysPrim && !isProdByGen) {
                  if (isWeakDecay) {
                    histos.fill(HIST("tracks/triton/dca/before/hDCAzVsPtantiTritonTrueSec"), track.pt(), track.dcaZ());
                  } else {
                    histos.fill(HIST("tracks/triton/dca/before/hDCAzVsPtantiTritonTrueMaterial"), track.pt(), track.dcaZ());
                  }
                }
              }
              break;
            case PDGHelium:
              if (isHeWoDCAz) {
                histos.fill(HIST("tracks/helium/dca/before/hDCAzVsPtHeliumTrue"), hePt, track.dcaZ());
                if (track.hasTOF() && outFlagOptions.doTOFplots) {
                  histos.fill(HIST("tracks/helium/dca/before/TOF/hDCAzVsPtHeliumTrue"), hePt, track.dcaZ());
                }
                if (isPhysPrim) {
                  histos.fill(HIST("tracks/helium/dca/before/hDCAzVsPtHeliumTruePrim"), hePt, track.dcaZ());
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    histos.fill(HIST("tracks/helium/dca/before/TOF/hDCAzVsPtHeliumTruePrim"), hePt, track.dcaZ());
                  }
                }
                if (!isPhysPrim && !isProdByGen) {
                  if (isWeakDecay) {
                    histos.fill(HIST("tracks/helium/dca/before/hDCAzVsPtHeliumTrueSec"), hePt, track.dcaZ());
                  } else {
                    histos.fill(HIST("tracks/helium/dca/before/hDCAzVsPtHeliumTrueMaterial"), hePt, track.dcaZ());
                  }
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    if (isWeakDecay) {
                      histos.fill(HIST("tracks/helium/dca/before/TOF/hDCAzVsPtHeliumTrueSec"), hePt, track.dcaZ());
                    } else {
                      histos.fill(HIST("tracks/helium/dca/before/TOF/hDCAzVsPtHeliumTrueMaterial"), hePt, track.dcaZ());
                    }
                  }
                }
              }
              if constexpr (!IsFilteredData) {
                if ((event.has_mcCollision() && (track.mcParticle().mcCollisionId() != event.mcCollisionId())) || !event.has_mcCollision()) {
                  if (isHeWoDCAz && outFlagOptions.makeWrongEventPlots) {
                    histos.fill(HIST("tracks/helium/dca/before/wrong/hDCAzVsPtHeliumTrue"), hePt, track.dcaZ());
                    if (track.hasTOF() && outFlagOptions.doTOFplots) {
                      histos.fill(HIST("tracks/helium/dca/before/wrong/TOF/hDCAzVsPtHeliumTrue"), hePt, track.dcaZ());
                    }
                    if (isPhysPrim) {
                      histos.fill(HIST("tracks/helium/dca/before/wrong/hDCAzVsPtHeliumTruePrim"), hePt, track.dcaZ());
                      if (track.hasTOF() && outFlagOptions.doTOFplots) {
                        histos.fill(HIST("tracks/helium/dca/before/wrong/TOF/hDCAzVsPtHeliumTruePrim"), hePt, track.dcaZ());
                      }
                    }
                    if (!isPhysPrim && !isProdByGen) {
                      histos.fill(HIST("tracks/helium/dca/before/wrong/hDCAzVsPtHeliumTrueTransport"), hePt, track.dcaZ());
                      if (isWeakDecay) {
                        histos.fill(HIST("tracks/helium/dca/before/wrong/hDCAzVsPtHeliumTrueSec"), hePt, track.dcaZ());
                      }
                      if (track.hasTOF() && outFlagOptions.doTOFplots) {
                        histos.fill(HIST("tracks/helium/dca/before/wrong/TOF/hDCAzVsPtHeliumTrueTransport"), hePt, track.dcaZ());
                        if (isWeakDecay) {
                          histos.fill(HIST("tracks/helium/dca/before/wrong/TOF/hDCAzVsPtHeliumTrueSec"), hePt, track.dcaZ());
                        }
                      }
                    }
                  }
                }
              }
              break;
            case -PDGHelium:
              if (isAntiHeWoDCAz) {
                histos.fill(HIST("tracks/helium/dca/before/hDCAzVsPtantiHeliumTrue"), antihePt, track.dcaZ());
                if (track.hasTOF() && outFlagOptions.doTOFplots) {
                  histos.fill(HIST("tracks/helium/dca/before/TOF/hDCAzVsPtantiHeliumTrue"), antihePt, track.dcaZ());
                }
                if (isPhysPrim) {
                  histos.fill(HIST("tracks/helium/dca/before/hDCAzVsPtantiHeliumTruePrim"), antihePt, track.dcaZ());
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    histos.fill(HIST("tracks/helium/dca/before/TOF/hDCAzVsPtantiHeliumTruePrim"), antihePt, track.dcaZ());
                  }
                }
                if (!isPhysPrim && !isProdByGen) {
                  if (isWeakDecay) {
                    histos.fill(HIST("tracks/helium/dca/before/hDCAzVsPtantiHeliumTrueSec"), antihePt, track.dcaZ());
                  } else {
                    histos.fill(HIST("tracks/helium/dca/before/hDCAzVsPtantiHeliumTrueMaterial"), antihePt, track.dcaZ());
                  }
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    if (isWeakDecay) {
                      histos.fill(HIST("tracks/helium/dca/before/TOF/hDCAzVsPtantiHeliumTrueSec"), antihePt, track.dcaZ());
                    } else {
                      histos.fill(HIST("tracks/helium/dca/before/TOF/hDCAzVsPtantiHeliumTrueMaterial"), antihePt, track.dcaZ());
                    }
                  }
                }
              }
              if constexpr (!IsFilteredData) {
                if ((event.has_mcCollision() && (track.mcParticle().mcCollisionId() != event.mcCollisionId())) || !event.has_mcCollision()) {
                  if (isAntiHeWoDCAz && outFlagOptions.makeWrongEventPlots) {
                    histos.fill(HIST("tracks/helium/dca/before/wrong/hDCAzVsPtantiHeliumTrue"), antihePt, track.dcaZ());
                    if (track.hasTOF() && outFlagOptions.doTOFplots) {
                      histos.fill(HIST("tracks/helium/dca/before/wrong/TOF/hDCAzVsPtantiHeliumTrue"), antihePt, track.dcaZ());
                    }
                    if (isPhysPrim) {
                      histos.fill(HIST("tracks/helium/dca/before/wrong/hDCAzVsPtantiHeliumTruePrim"), antihePt, track.dcaZ());
                      if (track.hasTOF() && outFlagOptions.doTOFplots) {
                        histos.fill(HIST("tracks/helium/dca/before/wrong/TOF/hDCAzVsPtantiHeliumTruePrim"), antihePt, track.dcaZ());
                      }
                    }
                    if (!isPhysPrim && !isProdByGen) {
                      histos.fill(HIST("tracks/helium/dca/before/wrong/hDCAzVsPtantiHeliumTrueTransport"), antihePt, track.dcaZ());
                      if (isWeakDecay) {
                        histos.fill(HIST("tracks/helium/dca/before/wrong/hDCAzVsPtantiHeliumTrueSec"), antihePt, track.dcaZ());
                      }
                      if (track.hasTOF() && outFlagOptions.doTOFplots) {
                        histos.fill(HIST("tracks/helium/dca/before/wrong/TOF/hDCAzVsPtantiHeliumTrueTransport"), antihePt, track.dcaZ());
                        if (isWeakDecay) {
                          histos.fill(HIST("tracks/helium/dca/before/wrong/TOF/hDCAzVsPtantiHeliumTrueSec"), antihePt, track.dcaZ());
                        }
                      }
                    }
                  }
                }
              }
              break;
            case PDGAlpha:
              if (enableAl && alRapCut && passDCAxyCut) {
                histos.fill(HIST("tracks/alpha/dca/before/hDCAzVsPtAlphaTrue"), track.pt(), track.dcaZ());
                if (isPhysPrim) {
                  histos.fill(HIST("tracks/alpha/dca/before/hDCAzVsPtAlphaTruePrim"), track.pt(), track.dcaZ());
                }
                if (!isPhysPrim && !isProdByGen) {
                  if (isWeakDecay) {
                    histos.fill(HIST("tracks/alpha/dca/before/hDCAzVsPtAlphaTrueSec"), track.pt(), track.dcaZ());
                  } else {
                    histos.fill(HIST("tracks/alpha/dca/before/hDCAzVsPtAlphaTrueMaterial"), track.pt(), track.dcaZ());
                  }
                }
              }
              break;
            case -PDGAlpha:
              if (enableAl && alRapCut && passDCAxyCut) {
                histos.fill(HIST("tracks/alpha/dca/before/hDCAzVsPtantiAlphaTrue"), track.pt(), track.dcaZ());
                if (isPhysPrim) {
                  histos.fill(HIST("tracks/alpha/dca/before/hDCAzVsPtantiAlphaTruePrim"), track.pt(), track.dcaZ());
                }
                if (!isPhysPrim && !isProdByGen) {
                  if (isWeakDecay) {
                    histos.fill(HIST("tracks/alpha/dca/before/hDCAzVsPtantiAlphaTrueSec"), track.pt(), track.dcaZ());
                  } else {
                    histos.fill(HIST("tracks/alpha/dca/before/hDCAzVsPtantiAlphaTrueMaterial"), track.pt(), track.dcaZ());
                  }
                }
              }
              break;
            default:
              break;
          }
          switch (std::abs(pdgCode)) {
            case PDGDeuteron:
              //
              break;
            default:
              if (isDeWoDCAzWTPCpid && outFlagOptions.makeFakeTracksPlots) {
                histos.fill(HIST("tracks/deuteron/dca/before/fake/hDCAzVsPtDeuteronTrue"), DPt, track.dcaZ());
                if (track.hasTOF() && outFlagOptions.doTOFplots) {
                  histos.fill(HIST("tracks/deuteron/dca/before/fake/TOF/hDCAzVsPtDeuteronTrue"), DPt, track.dcaZ());
                }
                if (isPhysPrim) {
                  histos.fill(HIST("tracks/deuteron/dca/before/fake/hDCAzVsPtDeuteronTruePrim"), DPt, track.dcaZ());
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    histos.fill(HIST("tracks/deuteron/dca/before/fake/TOF/hDCAzVsPtDeuteronTruePrim"), DPt, track.dcaZ());
                  }
                }
                if (!isPhysPrim && !isProdByGen) {
                  histos.fill(HIST("tracks/deuteron/dca/before/fake/hDCAzVsPtDeuteronTrueTransport"), DPt, track.dcaZ());
                  if (isWeakDecay) {
                    histos.fill(HIST("tracks/deuteron/dca/before/fake/hDCAzVsPtDeuteronTrueSec"), DPt, track.dcaZ());
                  }
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    histos.fill(HIST("tracks/deuteron/dca/before/fake/TOF/hDCAzVsPtDeuteronTrueTransport"), DPt, track.dcaZ());
                    if (isWeakDecay) {
                      histos.fill(HIST("tracks/deuteron/dca/before/fake/TOF/hDCAzVsPtDeuteronTrueSec"), DPt, track.dcaZ());
                    }
                  }
                }
              } else if (isAntiDeWoDCAzWTPCpid && outFlagOptions.makeFakeTracksPlots) {
                histos.fill(HIST("tracks/deuteron/dca/before/fake/hDCAzVsPtantiDeuteronTrue"), antiDPt, track.dcaZ());
                if (track.hasTOF() && outFlagOptions.doTOFplots) {
                  histos.fill(HIST("tracks/deuteron/dca/before/fake/TOF/hDCAzVsPtantiDeuteronTrue"), antiDPt, track.dcaZ());
                }
                if (isPhysPrim) {
                  histos.fill(HIST("tracks/deuteron/dca/before/fake/hDCAzVsPtantiDeuteronTruePrim"), antiDPt, track.dcaZ());
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    histos.fill(HIST("tracks/deuteron/dca/before/fake/TOF/hDCAzVsPtantiDeuteronTruePrim"), antiDPt, track.dcaZ());
                  }
                }
                if (!isPhysPrim && !isProdByGen) {
                  histos.fill(HIST("tracks/deuteron/dca/before/fake/hDCAzVsPtantiDeuteronTrueTransport"), antiDPt, track.dcaZ());
                  if (isWeakDecay) {
                    histos.fill(HIST("tracks/deuteron/dca/before/fake/hDCAzVsPtantiDeuteronTrueSec"), antiDPt, track.dcaZ());
                  }
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    histos.fill(HIST("tracks/deuteron/dca/before/fake/TOF/hDCAzVsPtantiDeuteronTrueTransport"), antiDPt, track.dcaZ());
                    if (isWeakDecay) {
                      histos.fill(HIST("tracks/deuteron/dca/before/fake/TOF/hDCAzVsPtantiDeuteronTrueSec"), antiDPt, track.dcaZ());
                    }
                  }
                }
              }
              break;
          }

          switch (std::abs(pdgCode)) {
            case PDGHelium:
              //
              break;
            default:
              if (isHeWoDCAzWTPCpid && outFlagOptions.makeFakeTracksPlots) {
                histos.fill(HIST("tracks/helium/dca/before/fake/hDCAzVsPtHeliumTrue"), hePt, track.dcaZ());
                if (track.hasTOF() && outFlagOptions.doTOFplots) {
                  histos.fill(HIST("tracks/helium/dca/before/fake/TOF/hDCAzVsPtHeliumTrue"), hePt, track.dcaZ());
                }
                if (isPhysPrim) {
                  histos.fill(HIST("tracks/helium/dca/before/fake/hDCAzVsPtHeliumTruePrim"), hePt, track.dcaZ());
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    histos.fill(HIST("tracks/helium/dca/before/fake/TOF/hDCAzVsPtHeliumTruePrim"), hePt, track.dcaZ());
                  }
                }
                if (!isPhysPrim && !isProdByGen) {
                  histos.fill(HIST("tracks/helium/dca/before/fake/hDCAzVsPtHeliumTrueTransport"), hePt, track.dcaZ());
                  if (isWeakDecay) {
                    histos.fill(HIST("tracks/helium/dca/before/fake/hDCAzVsPtHeliumTrueSec"), hePt, track.dcaZ());
                  }
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    histos.fill(HIST("tracks/helium/dca/before/fake/TOF/hDCAzVsPtHeliumTrueTransport"), hePt, track.dcaZ());
                    if (isWeakDecay) {
                      histos.fill(HIST("tracks/helium/dca/before/fake/TOF/hDCAzVsPtHeliumTrueSec"), hePt, track.dcaZ());
                    }
                  }
                }
              }
              if (isAntiHeWoDCAzWTPCpid && outFlagOptions.makeFakeTracksPlots) {
                histos.fill(HIST("tracks/helium/dca/before/fake/hDCAzVsPtantiHeliumTrue"), antihePt, track.dcaZ());
                if (track.hasTOF() && outFlagOptions.doTOFplots) {
                  histos.fill(HIST("tracks/helium/dca/before/fake/TOF/hDCAzVsPtantiHeliumTrue"), antihePt, track.dcaZ());
                }
                if (isPhysPrim) {
                  histos.fill(HIST("tracks/helium/dca/before/fake/hDCAzVsPtantiHeliumTruePrim"), antihePt, track.dcaZ());
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    histos.fill(HIST("tracks/helium/dca/before/fake/TOF/hDCAzVsPtantiHeliumTruePrim"), antihePt, track.dcaZ());
                  }
                }
                if (!isPhysPrim && !isProdByGen) {
                  histos.fill(HIST("tracks/helium/dca/before/fake/hDCAzVsPtantiHeliumTrueTransport"), antihePt, track.dcaZ());
                  if (isWeakDecay) {
                    histos.fill(HIST("tracks/helium/dca/before/fake/hDCAzVsPtantiHeliumTrueSec"), antihePt, track.dcaZ());
                  }
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    histos.fill(HIST("tracks/helium/dca/before/fake/TOF/hDCAzVsPtantiHeliumTrueTransport"), antihePt, track.dcaZ());
                    if (isWeakDecay) {
                      histos.fill(HIST("tracks/helium/dca/before/fake/TOF/hDCAzVsPtantiHeliumTrueSec"), antihePt, track.dcaZ());
                    }
                  }
                }
              }
              break;
          }
        }
      } else {
        (void)particles;
      }
      // Tracks DCA histos fill
      if (outFlagOptions.makeDCABeforeCutPlots) {
        if (passDCAzCut) {
          histos.fill(HIST("tracks/dca/before/hDCAxyVsPt"), track.pt(), track.dcaXY());

          if (enablePr && prRapCut) {
            if (std::abs(track.tpcNSigmaPr()) < nsigmaTPCvar.nsigmaTPCPr) {
              if (track.sign() > 0) {
                histos.fill(HIST("tracks/proton/dca/before/hDCAxyVsPtProton"), track.pt(), track.dcaXY());
              } else {
                histos.fill(HIST("tracks/proton/dca/before/hDCAxyVsPtantiProton"), track.pt(), track.dcaXY());
              }
            }
          }

          if (enableTr && trRapCut) {
            if (isTritonTPCpid) {
              if (track.sign() > 0) {
                histos.fill(HIST("tracks/triton/dca/before/hDCAxyVsPtTriton"), track.pt(), track.dcaXY());
              } else {
                histos.fill(HIST("tracks/triton/dca/before/hDCAxyVsPtantiTriton"), track.pt(), track.dcaXY());
              }
            }
          }
          if (enableAl && alRapCut) {
            if (std::abs(track.tpcNSigmaAl()) < nsigmaTPCvar.nsigmaTPCAl) {
              if (track.sign() > 0) {
                histos.fill(HIST("tracks/alpha/dca/before/hDCAxyVsPtAlpha"), track.pt(), track.dcaXY());
              } else {
                histos.fill(HIST("tracks/alpha/dca/before/hDCAxyVsPtantiAlpha"), track.pt(), track.dcaXY());
              }
            }
          }
        }

        if (isDeWoDCAxyWTPCpid) {
          if (usenITSLayer && !itsClusterMap.test(trkqcOptions.nITSLayer))
            continue;
          if (enableCentrality)
            histos.fill(HIST("tracks/deuteron/dca/before/hDCAxyVsPtDeuteronVsMult"), DPt, track.dcaXY(), event.centFT0M());
          else
            histos.fill(HIST("tracks/deuteron/dca/before/hDCAxyVsPtDeuteron"), DPt, track.dcaXY());
          if (!track.hasTOF() && (outFlagOptions.enableNoTOFPlots))
            histos.fill(HIST("tracks/deuteron/dca/before/hDCAxyVsPtDeuteronNoTOF"), DPt, track.dcaXY());
        }
        if (isAntiDeWoDCAxyWTPCpid) {
          if (usenITSLayer && !itsClusterMap.test(trkqcOptions.nITSLayer))
            continue;
          if (enableCentrality)
            histos.fill(HIST("tracks/deuteron/dca/before/hDCAxyVsPtantiDeuteronVsMult"), antiDPt, track.dcaXY(), event.centFT0M());
          else
            histos.fill(HIST("tracks/deuteron/dca/before/hDCAxyVsPtantiDeuteron"), antiDPt, track.dcaXY());
          if (!track.hasTOF() && (outFlagOptions.enableNoTOFPlots))
            histos.fill(HIST("tracks/deuteron/dca/before/hDCAxyVsPtantiDeuteronNoTOF"), antiDPt, track.dcaXY());
        }

        if (isHeWoDCAxyWTPCpid) {
          histos.fill(HIST("tracks/helium/dca/before/hDCAxyVsPtHelium"), hePt, track.dcaXY());
          if (!track.hasTOF() && (outFlagOptions.enableNoTOFPlots))
            histos.fill(HIST("tracks/helium/dca/before/hDCAxyVsPtHeliumNoTOF"), hePt, track.dcaXY());
          if (track.hasTOF() && outFlagOptions.doTOFplots) {
            histos.fill(HIST("tracks/helium/dca/before/TOF/hDCAxyVsPtHelium"), hePt, track.dcaXY());
          }
        }
        if (isAntiHeWoDCAxyWTPCpid) {
          histos.fill(HIST("tracks/helium/dca/before/hDCAxyVsPtantiHelium"), antihePt, track.dcaXY());
          if (!track.hasTOF() && (outFlagOptions.enableNoTOFPlots))
            histos.fill(HIST("tracks/helium/dca/before/hDCAxyVsPtantiHeliumNoTOF"), antihePt, track.dcaXY());
          if (track.hasTOF() && outFlagOptions.doTOFplots) {
            histos.fill(HIST("tracks/helium/dca/before/TOF/hDCAxyVsPtantiHelium"), antihePt, track.dcaXY());
          }
        }
      }

      if constexpr (IsMC) {
        // auto const& mcParticles = particles;
        bool isPhysPrim = false;
        bool isProdByGen = false;
        bool isWeakDecay = false;
        bool hasFakeHit = false;

        // PID
        int pdgCode = 0;
        int pdgMom = 0;
        // gen Pt
        float genPt = 0;
        float ptMom = 0;
        // Mothers variables
        [[maybe_unused]] int firstMotherId = -1;
        [[maybe_unused]] int firstMotherPdg = -1;
        [[maybe_unused]] float firstMotherPt = -1.f;
        [[maybe_unused]] int pdgMomList[kMaxNumMom];
        [[maybe_unused]] float ptMomList[kMaxNumMom];
        [[maybe_unused]] int nSaved = 0;

        if constexpr (IsFilteredData) {
          isPhysPrim = track.isPhysicalPrimary();
          isProdByGen = track.producedByGenerator();
          isWeakDecay = (track.getProcess() == TMCProcess::kPDecay);
          pdgCode = track.pdgCode();
          genPt = std::sqrt(track.px() * track.px() + track.py() * track.py());

        } else {
          if (!track.has_mcParticle()) {
            continue;
          }
          isPhysPrim = track.mcParticle().isPhysicalPrimary();
          isProdByGen = track.mcParticle().producedByGenerator();
          isWeakDecay = (track.mcParticle().getProcess() == TMCProcess::kPDecay);
          pdgCode = track.mcParticle().pdgCode();

          // Access to MC particles mother
          o2::aod::McParticles::iterator mc = particles.iteratorAt(track.mcParticleId());
          gsl::span<const int> motherIds = mc.mothersIds();
          const int nMothers = static_cast<int>(motherIds.size());
          firstMotherId = -1;
          firstMotherPdg = -1;
          firstMotherPt = -1.f;
          nSaved = 0;

          for (int iMom = 0; iMom < nMothers; iMom++) {
            int motherId = motherIds[iMom];
            if (motherId < 0 || motherId >= particles.size()) {
              continue; // added check on mother
            }
            o2::aod::McParticles::iterator mother = particles.iteratorAt(motherId);
            pdgMom = mother.pdgCode();
            ptMom = mother.pt();

            if (iMom == 0) {
              firstMotherId = motherId;
              firstMotherPdg = pdgMom;
              firstMotherPt = ptMom;
            }
            if (nSaved < kMaxNumMom) {
              pdgMomList[nSaved] = pdgMom;
              ptMomList[nSaved] = ptMom;
              nSaved++;
            }
          }

          genPt = track.mcParticle().pt();
          for (int i = 0; i < kFakeLoop; i++) { // From ITS to TPC
            if (track.mcMask() & 1 << i) {
              hasFakeHit = true;
              break;
            }
          }
          if (hasFakeHit && passDCAzCut) {
            debugHistos.fill(HIST("debug/h2TPCsignVsTPCmomentum_FakeHits"), track.tpcInnerParam() / (1.f * track.sign()), track.tpcSignal());
          }
        }
        switch (pdgCode) {
          case PDGProton:
            if (enablePr && prRapCut && outFlagOptions.makeDCABeforeCutPlots && passDCAzCut) {
              histos.fill(HIST("tracks/proton/dca/before/hDCAxyVsPtProtonTrue"), track.pt(), track.dcaXY());
              if (track.hasTOF() && outFlagOptions.doTOFplots) {
                histos.fill(HIST("tracks/proton/dca/before/TOF/hDCAxyVsPtProtonTrue"), track.pt(), track.dcaXY());
              }
              if (isPhysPrim) {
                histos.fill(HIST("tracks/proton/dca/before/hDCAxyVsPtProtonTruePrim"), track.pt(), track.dcaXY());
                if (track.hasTOF() && outFlagOptions.doTOFplots) {
                  histos.fill(HIST("tracks/proton/dca/before/TOF/hDCAxyVsPtProtonTruePrim"), track.pt(), track.dcaXY());
                }
              }
              if (!isPhysPrim && !isProdByGen) {
                if (isWeakDecay) {
                  histos.fill(HIST("tracks/proton/dca/before/hDCAxyVsPtProtonTrueSec"), track.pt(), track.dcaXY());
                } else {
                  histos.fill(HIST("tracks/proton/dca/before/hDCAxyVsPtProtonTrueMaterial"), track.pt(), track.dcaXY());
                  if constexpr (!IsFilteredData) {
                    histos.fill(HIST("tracks/proton/dca/before/hNumMothers"), nSaved);
                    if (nSaved > 0) {
                      for (int iMom = 0; iMom < nSaved; iMom++) {
                        int pdgMom = pdgMomList[iMom];
                        float pdgSign = (pdgMom > 0) ? 1.0 : -1.0;
                        float ptMom = ptMomList[iMom];
                        int motherSpeciesBin = -1;
                        if (pdgMom != -1) {
                          motherSpeciesBin = 0;
                          for (int j = 0; j < kNumMotherList; j++) {
                            if (std::abs(kPdgMotherList[j]) == std::abs(pdgMom)) {
                              motherSpeciesBin = j + 1;
                              break;
                            }
                          }
                        }
                        histos.fill(HIST("tracks/proton/dca/before/hMomTrueMaterial"), pdgSign, motherSpeciesBin, ptMom);
                      }
                    }
                  }
                }
                if (track.hasTOF() && outFlagOptions.doTOFplots) {
                  if (isWeakDecay) {
                    histos.fill(HIST("tracks/proton/dca/before/TOF/hDCAxyVsPtProtonTrueSec"), track.pt(), track.dcaXY());
                  } else {
                    histos.fill(HIST("tracks/proton/dca/before/TOF/hDCAxyVsPtProtonTrueMaterial"), track.pt(), track.dcaXY());
                  }
                }
              }
            }
            break;
          case -PDGProton:
            if (enablePr && prRapCut && outFlagOptions.makeDCABeforeCutPlots && passDCAzCut) {
              histos.fill(HIST("tracks/proton/dca/before/hDCAxyVsPtantiProtonTrue"), track.pt(), track.dcaXY());
              if (track.hasTOF() && outFlagOptions.doTOFplots) {
                histos.fill(HIST("tracks/proton/dca/before/TOF/hDCAxyVsPtantiProtonTrue"), track.pt(), track.dcaXY());
              }
              if (isPhysPrim) {
                histos.fill(HIST("tracks/proton/dca/before/hDCAxyVsPtantiProtonTruePrim"), track.pt(), track.dcaXY());
                if (track.hasTOF() && outFlagOptions.doTOFplots) {
                  histos.fill(HIST("tracks/proton/dca/before/TOF/hDCAxyVsPtantiProtonTruePrim"), track.pt(), track.dcaXY());
                }
              }
              if (!isPhysPrim && !isProdByGen) {
                if (isWeakDecay) {
                  histos.fill(HIST("tracks/proton/dca/before/hDCAxyVsPtantiProtonTrueSec"), track.pt(), track.dcaXY());
                } else {
                  histos.fill(HIST("tracks/proton/dca/before/hDCAxyVsPtantiProtonTrueMaterial"), track.pt(), track.dcaXY());
                }
                if (track.hasTOF() && outFlagOptions.doTOFplots) {
                  if (isWeakDecay) {
                    histos.fill(HIST("tracks/proton/dca/before/TOF/hDCAxyVsPtantiProtonTrueSec"), track.pt(), track.dcaXY());
                  } else {
                    histos.fill(HIST("tracks/proton/dca/before/TOF/hDCAxyVsPtantiProtonTrueMaterial"), track.pt(), track.dcaXY());
                  }
                }
              }
            }
            break;
          case PDGDeuteron:
            if (isDeWoDCAxy) {
              if (outFlagOptions.makeDCABeforeCutPlots) {
                histos.fill(HIST("tracks/deuteron/dca/before/hDCAxyVsPtDeuteronTrue"), DPt, track.dcaXY());
                if (track.hasTOF() && outFlagOptions.doTOFplots) {
                  histos.fill(HIST("tracks/deuteron/dca/before/TOF/hDCAxyVsPtDeuteronTrue"), DPt, track.dcaXY());
                }
              }
              if (isPhysPrim) {
                if (outFlagOptions.makeDCABeforeCutPlots) {
                  histos.fill(HIST("tracks/deuteron/dca/before/hDCAxyVsPtDeuteronTruePrim"), DPt, track.dcaXY());
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    histos.fill(HIST("tracks/deuteron/dca/before/TOF/hDCAxyVsPtDeuteronTruePrim"), DPt, track.dcaXY());
                  }
                }
                if constexpr (IsFilteredData) {
                  spectraGen.fill(HIST("deuteron/histDeuteronPtShift"), track.pt(), track.pt() - genPt);
                  spectraGen.fill(HIST("deuteron/histDeuteronPtShiftVsEta"), track.eta(), track.pt() - genPt);

                  histos.fill(HIST("tracks/deuteron/histDeuteronPtShiftRec"), DPt);
                  histos.fill(HIST("tracks/deuteron/histDeuteronPtRec"), track.pt());
                  spectraGen.fill(HIST("deuteron/histDeuteronPtShiftCorrection"), DPt, DPt - genPt);
                } else {
                  spectraGen.fill(HIST("deuteron/histDeuteronPtShift"), track.pt(), track.pt() - genPt);
                  spectraGen.fill(HIST("deuteron/histDeuteronPtShiftVsEta"), track.eta(), track.pt() - genPt);

                  histos.fill(HIST("tracks/deuteron/histDeuteronPtShiftRec"), DPt);
                  histos.fill(HIST("tracks/deuteron/histDeuteronPtRec"), track.pt());
                  spectraGen.fill(HIST("deuteron/histDeuteronPtShiftCorrection"), DPt, DPt - genPt);
                }
              }
              if (!isPhysPrim && !isProdByGen) {
                if (outFlagOptions.makeDCABeforeCutPlots) {
                  if (isWeakDecay) {
                    histos.fill(HIST("tracks/deuteron/dca/before/hDCAxyVsPtDeuteronTrueSec"), DPt, track.dcaXY());
                  } else {
                    histos.fill(HIST("tracks/deuteron/dca/before/hDCAxyVsPtDeuteronTrueMaterial"), DPt, track.dcaXY());
                    if constexpr (!IsFilteredData) {
                      histos.fill(HIST("tracks/deuteron/dca/before/hNumMothers"), nSaved);
                      if (nSaved > 0) {
                        for (int iMom = 0; iMom < nSaved; iMom++) {
                          int pdgMom = pdgMomList[iMom];
                          float pdgSign = (pdgMom > 0) ? 1.0 : -1.0;
                          float ptMom = ptMomList[iMom];
                          int motherSpeciesBin = -1;
                          if (pdgMom != -1) {
                            motherSpeciesBin = 0;
                            for (int j = 0; j < kNumMotherList; j++) {
                              if (std::abs(kPdgMotherList[j]) == std::abs(pdgMom)) {
                                motherSpeciesBin = j + 1;
                                break;
                              }
                            }
                          }
                          histos.fill(HIST("tracks/deuteron/dca/before/hMomTrueMaterial"), pdgSign, motherSpeciesBin, ptMom);
                        }
                      }
                    }
                  }
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    if (isWeakDecay) {
                      histos.fill(HIST("tracks/deuteron/dca/before/TOF/hDCAxyVsPtDeuteronTrueSec"), DPt, track.dcaXY());
                    } else {
                      histos.fill(HIST("tracks/deuteron/dca/before/TOF/hDCAxyVsPtDeuteronTrueMaterial"), DPt, track.dcaXY());
                    }
                  }
                }
              }
            }
            break;
          case -PDGDeuteron:
            if (isAntiDeWoDCAxy) {
              if (outFlagOptions.makeDCABeforeCutPlots) {
                histos.fill(HIST("tracks/deuteron/dca/before/hDCAxyVsPtantiDeuteronTrue"), antiDPt, track.dcaXY());
                if (track.hasTOF() && outFlagOptions.doTOFplots) {
                  histos.fill(HIST("tracks/deuteron/dca/before/TOF/hDCAxyVsPtantiDeuteronTrue"), antiDPt, track.dcaXY());
                }
              }
              if (isPhysPrim) {
                if (outFlagOptions.makeDCABeforeCutPlots) {
                  histos.fill(HIST("tracks/deuteron/dca/before/hDCAxyVsPtantiDeuteronTruePrim"), antiDPt, track.dcaXY());
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    histos.fill(HIST("tracks/deuteron/dca/before/TOF/hDCAxyVsPtantiDeuteronTruePrim"), antiDPt, track.dcaXY());
                  }
                }
                if constexpr (IsFilteredData) {
                  spectraGen.fill(HIST("deuteron/histAntiDeuteronPtShift"), track.pt(), track.pt() - genPt);
                  spectraGen.fill(HIST("deuteron/histAntiDeuteronPtShiftVsEta"), track.eta(), track.pt() - genPt);

                  histos.fill(HIST("tracks/deuteron/histAntiDeuteronPtShiftRec"), antiDPt);
                  histos.fill(HIST("tracks/deuteron/histAntiDeuteronPtRec"), track.pt());
                  spectraGen.fill(HIST("deuteron/histAntiDeuteronPtShiftCorrection"), antiDPt, antiDPt - genPt);
                } else {
                  spectraGen.fill(HIST("deuteron/histAntiDeuteronPtShift"), track.pt(), track.pt() - genPt);
                  spectraGen.fill(HIST("deuteron/histAntiDeuteronPtShiftVsEta"), track.eta(), track.pt() - genPt);

                  histos.fill(HIST("tracks/deuteron/histAntiDeuteronPtShiftRec"), antiDPt);
                  histos.fill(HIST("tracks/deuteron/histAntiDeuteronPtRec"), track.pt());
                  spectraGen.fill(HIST("deuteron/histAntiDeuteronPtShiftCorrection"), antiDPt, antiDPt - genPt);
                }
              }
              if (!isPhysPrim && !isProdByGen) {
                if (outFlagOptions.makeDCABeforeCutPlots) {
                  if (isWeakDecay) {
                    histos.fill(HIST("tracks/deuteron/dca/before/hDCAxyVsPtantiDeuteronTrueSec"), antiDPt, track.dcaXY());
                  } else {
                    histos.fill(HIST("tracks/deuteron/dca/before/hDCAxyVsPtantiDeuteronTrueMaterial"), antiDPt, track.dcaXY());
                  }
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    if (isWeakDecay) {
                      histos.fill(HIST("tracks/deuteron/dca/before/TOF/hDCAxyVsPtantiDeuteronTrueSec"), antiDPt, track.dcaXY());
                    } else {
                      histos.fill(HIST("tracks/deuteron/dca/before/TOF/hDCAxyVsPtantiDeuteronTrueMaterial"), antiDPt, track.dcaXY());
                    }
                  }
                }
              }
            }
            break;
          case PDGTriton:
            if (enableTr && trRapCut && outFlagOptions.makeDCABeforeCutPlots && passDCAzCut) {
              histos.fill(HIST("tracks/triton/dca/before/hDCAxyVsPtTritonTrue"), track.pt(), track.dcaXY());
              if (track.hasTOF() && outFlagOptions.doTOFplots) {
                histos.fill(HIST("tracks/triton/dca/before/TOF/hDCAxyVsPtTritonTrue"), track.pt(), track.dcaXY());
              }
              if (isPhysPrim) {
                histos.fill(HIST("tracks/triton/dca/before/hDCAxyVsPtTritonTruePrim"), track.pt(), track.dcaXY());
                if (track.hasTOF() && outFlagOptions.doTOFplots) {
                  histos.fill(HIST("tracks/triton/dca/before/TOF/hDCAxyVsPtTritonTruePrim"), track.pt(), track.dcaXY());
                }
              }
              if (!isPhysPrim && !isProdByGen) {
                if (isWeakDecay) {
                  histos.fill(HIST("tracks/triton/dca/before/hDCAxyVsPtTritonTrueSec"), track.pt(), track.dcaXY());
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    histos.fill(HIST("tracks/triton/dca/before/TOF/hDCAxyVsPtTritonTrueSec"), track.pt(), track.dcaXY());
                  }
                } else {
                  histos.fill(HIST("tracks/triton/dca/before/hDCAxyVsPtTritonTrueMaterial"), track.pt(), track.dcaXY());
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    histos.fill(HIST("tracks/triton/dca/before/TOF/hDCAxyVsPtTritonTrueMaterial"), track.pt(), track.dcaXY());
                  }
                }
              }
            }
            break;
          case -PDGTriton:
            if (enableTr && trRapCut && outFlagOptions.makeDCABeforeCutPlots && passDCAzCut) {
              histos.fill(HIST("tracks/triton/dca/before/hDCAxyVsPtantiTritonTrue"), track.pt(), track.dcaXY());
              if (track.hasTOF() && outFlagOptions.doTOFplots) {
                histos.fill(HIST("tracks/triton/dca/before/TOF/hDCAxyVsPtantiTritonTrue"), track.pt(), track.dcaXY());
              }

              if (isPhysPrim) {
                histos.fill(HIST("tracks/triton/dca/before/hDCAxyVsPtantiTritonTruePrim"), track.pt(), track.dcaXY());
                if (track.hasTOF() && outFlagOptions.doTOFplots) {
                  histos.fill(HIST("tracks/triton/dca/before/TOF/hDCAxyVsPtantiTritonTruePrim"), track.pt(), track.dcaXY());
                }
              }
              if (!isPhysPrim && isProdByGen) {
                //
              }
              if (!isPhysPrim && !isProdByGen) {
                if (isWeakDecay) {
                  histos.fill(HIST("tracks/triton/dca/before/hDCAxyVsPtantiTritonTrueSec"), track.pt(), track.dcaXY());
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    histos.fill(HIST("tracks/triton/dca/before/TOF/hDCAxyVsPtantiTritonTrueSec"), track.pt(), track.dcaXY());
                  }
                } else {
                  histos.fill(HIST("tracks/triton/dca/before/hDCAxyVsPtantiTritonTrueMaterial"), track.pt(), track.dcaXY());
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    histos.fill(HIST("tracks/triton/dca/before/TOF/hDCAxyVsPtantiTritonTrueMaterial"), track.pt(), track.dcaXY());
                  }
                }
              }
            }
            break;
          case PDGHelium:
            if (isHeWoDCAxy) {
              if (outFlagOptions.makeDCABeforeCutPlots) {
                histos.fill(HIST("tracks/helium/dca/before/hDCAxyVsPtHeliumTrue"), hePt, track.dcaXY());
                if (track.hasTOF() && outFlagOptions.doTOFplots) {
                  histos.fill(HIST("tracks/helium/dca/before/TOF/hDCAxyVsPtHeliumTrue"), hePt, track.dcaXY());
                }
              }
              if (isPhysPrim) {
                if constexpr (!IsFilteredData) {
                  spectraGen.fill(HIST("helium/histPtGenHe"), std::abs(track.mcParticle().pt()));
                  spectraGen.fill(HIST("helium/histPtRecHe"), 2.f * hePt);
                  spectraGen.fill(HIST("helium/histPtShiftHe"), 2.f * hePt, 2.f * hePt - track.mcParticle().pt());
                  spectraGen.fill(HIST("helium/histPtShiftVsEtaHe"), track.eta(), 2.f * hePt - track.mcParticle().pt());

                  spectraGen.fill(HIST("helium/histPGenHe"), std::abs(track.mcParticle().p()));
                  spectraGen.fill(HIST("helium/histPRecHe"), 2.f * heP);
                  spectraGen.fill(HIST("helium/histPShiftHe"), 2.f * heP, 2.f * heP - track.mcParticle().p());
                  spectraGen.fill(HIST("helium/histPShiftVsEtaHe"), track.eta(), 2.f * heP - track.mcParticle().p());
                }
                if (outFlagOptions.makeDCABeforeCutPlots) {
                  histos.fill(HIST("tracks/helium/dca/before/hDCAxyVsPtHeliumTruePrim"), hePt, track.dcaXY());
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    histos.fill(HIST("tracks/helium/dca/before/TOF/hDCAxyVsPtHeliumTruePrim"), hePt, track.dcaXY());
                  }
                }
              }
              if (!isPhysPrim && !isProdByGen && outFlagOptions.makeDCABeforeCutPlots) {
                if (isWeakDecay) {
                  histos.fill(HIST("tracks/helium/dca/before/hDCAxyVsPtHeliumTrueSec"), hePt, track.dcaXY());
                } else {
                  histos.fill(HIST("tracks/helium/dca/before/hDCAxyVsPtHeliumTrueMaterial"), hePt, track.dcaXY());
                  if (!IsFilteredData) {
                    histos.fill(HIST("tracks/helium/dca/before/hNumMothers"), nSaved);
                    if (nSaved > 0) {
                      for (int iMom = 0; iMom < nSaved; iMom++) {
                        int pdgMom = pdgMomList[iMom];
                        float pdgSign = (pdgMom > 0) ? 1.0 : -1.0;
                        float ptMom = ptMomList[iMom];
                        int motherSpeciesBin = -1;
                        if (pdgMom != -1) {
                          motherSpeciesBin = 0;
                          for (int j = 0; j < kNumMotherList; j++) {
                            if (std::abs(kPdgMotherList[j]) == std::abs(pdgMom)) {
                              motherSpeciesBin = j + 1;
                              break;
                            }
                          }
                        }
                        histos.fill(HIST("tracks/helium/dca/before/hMomTrueMaterial"), pdgSign, motherSpeciesBin, ptMom);
                      }
                    }
                  }
                }
                if (track.hasTOF() && outFlagOptions.doTOFplots) {
                  if (isWeakDecay) {
                    histos.fill(HIST("tracks/helium/dca/before/TOF/hDCAxyVsPtHeliumTrueSec"), hePt, track.dcaXY());
                  } else {
                    histos.fill(HIST("tracks/helium/dca/before/TOF/hDCAxyVsPtHeliumTrueMaterial"), hePt, track.dcaXY());
                  }
                }
              }
            }
            if constexpr (!IsFilteredData) {
              if ((event.has_mcCollision() && (track.mcParticle().mcCollisionId() != event.mcCollisionId())) || !event.has_mcCollision()) {
                if (isHeWoDCAxy && outFlagOptions.makeDCABeforeCutPlots && outFlagOptions.makeWrongEventPlots) {
                  histos.fill(HIST("tracks/helium/dca/before/wrong/hDCAxyVsPtHeliumTrue"), hePt, track.dcaXY());
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    histos.fill(HIST("tracks/helium/dca/before/wrong/TOF/hDCAxyVsPtHeliumTrue"), hePt, track.dcaXY());
                  }
                  if (isPhysPrim) {
                    histos.fill(HIST("tracks/helium/dca/before/wrong/hDCAxyVsPtHeliumTruePrim"), hePt, track.dcaXY());
                    if (track.hasTOF() && outFlagOptions.doTOFplots) {
                      histos.fill(HIST("tracks/helium/dca/before/wrong/TOF/hDCAxyVsPtHeliumTruePrim"), hePt, track.dcaXY());
                    }
                  }
                  if (!isPhysPrim && !isProdByGen) {
                    histos.fill(HIST("tracks/helium/dca/before/wrong/hDCAxyVsPtHeliumTrueTransport"), hePt, track.dcaXY());
                    if (isWeakDecay) {
                      histos.fill(HIST("tracks/helium/dca/before/wrong/hDCAxyVsPtHeliumTrueSec"), hePt, track.dcaXY());
                    }
                    if (track.hasTOF() && outFlagOptions.doTOFplots) {
                      histos.fill(HIST("tracks/helium/dca/before/wrong/TOF/hDCAxyVsPtHeliumTrueTransport"), hePt, track.dcaXY());
                      if (isWeakDecay) {
                        histos.fill(HIST("tracks/helium/dca/before/wrong/TOF/hDCAxyVsPtHeliumTrueSec"), hePt, track.dcaXY());
                      }
                    }
                  }
                }
              }
            }
            break;
          case -PDGHelium:
            if (isAntiHeWoDCAxy) {
              if (outFlagOptions.makeDCABeforeCutPlots) {
                histos.fill(HIST("tracks/helium/dca/before/hDCAxyVsPtantiHeliumTrue"), antihePt, track.dcaXY());
                if (track.hasTOF() && outFlagOptions.doTOFplots) {
                  histos.fill(HIST("tracks/helium/dca/before/TOF/hDCAxyVsPtantiHeliumTrue"), antihePt, track.dcaXY());
                }
              }
              if (isPhysPrim) {
                if constexpr (!IsFilteredData) {
                  spectraGen.fill(HIST("helium/histPtGenantiHe"), std::abs(track.mcParticle().pt()));
                  spectraGen.fill(HIST("helium/histPtRecantiHe"), 2.f * antihePt);
                  spectraGen.fill(HIST("helium/histPtShiftantiHe"), 2.f * antihePt, 2.f * antihePt - track.mcParticle().pt());
                  spectraGen.fill(HIST("helium/histPtShiftVsEtaantiHe"), track.eta(), 2.f * antihePt - track.mcParticle().pt());

                  spectraGen.fill(HIST("helium/histPGenantiHe"), std::abs(track.mcParticle().p()));
                  spectraGen.fill(HIST("helium/histPRecantiHe"), 2.f * antiheP);
                  spectraGen.fill(HIST("helium/histPShiftantiHe"), 2.f * antiheP, 2.f * antiheP - track.mcParticle().p());
                  spectraGen.fill(HIST("helium/histPShiftVsEtaantiHe"), track.eta(), 2.f * antiheP - track.mcParticle().p());
                }
                if (outFlagOptions.makeDCABeforeCutPlots) {
                  histos.fill(HIST("tracks/helium/dca/before/hDCAxyVsPtantiHeliumTruePrim"), antihePt, track.dcaXY());
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    histos.fill(HIST("tracks/helium/dca/before/TOF/hDCAxyVsPtantiHeliumTruePrim"), antihePt, track.dcaXY());
                  }
                }
              }
              if (!isPhysPrim && !isProdByGen && outFlagOptions.makeDCABeforeCutPlots) {
                if (isWeakDecay) {
                  histos.fill(HIST("tracks/helium/dca/before/hDCAxyVsPtantiHeliumTrueSec"), antihePt, track.dcaXY());
                } else {
                  histos.fill(HIST("tracks/helium/dca/before/hDCAxyVsPtantiHeliumTrueMaterial"), antihePt, track.dcaXY());
                }
                if (track.hasTOF() && outFlagOptions.doTOFplots) {
                  if (isWeakDecay) {
                    histos.fill(HIST("tracks/helium/dca/before/TOF/hDCAxyVsPtantiHeliumTrueSec"), antihePt, track.dcaXY());
                  } else {
                    histos.fill(HIST("tracks/helium/dca/before/TOF/hDCAxyVsPtantiHeliumTrueMaterial"), antihePt, track.dcaXY());
                  }
                }
              }
            }
            if constexpr (!IsFilteredData) {
              if ((event.has_mcCollision() && (track.mcParticle().mcCollisionId() != event.mcCollisionId())) || !event.has_mcCollision()) {
                if (isAntiHeWoDCAxy && outFlagOptions.makeDCABeforeCutPlots && outFlagOptions.makeWrongEventPlots) {
                  histos.fill(HIST("tracks/helium/dca/before/wrong/hDCAxyVsPtantiHeliumTrue"), antihePt, track.dcaXY());
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    histos.fill(HIST("tracks/helium/dca/before/wrong/TOF/hDCAxyVsPtantiHeliumTrue"), antihePt, track.dcaXY());
                  }
                  if (isPhysPrim) {
                    histos.fill(HIST("tracks/helium/dca/before/wrong/hDCAxyVsPtantiHeliumTruePrim"), antihePt, track.dcaXY());
                    if (track.hasTOF() && outFlagOptions.doTOFplots) {
                      histos.fill(HIST("tracks/helium/dca/before/wrong/TOF/hDCAxyVsPtantiHeliumTruePrim"), antihePt, track.dcaXY());
                    }
                  }
                  if (!isPhysPrim && !isProdByGen) {
                    histos.fill(HIST("tracks/helium/dca/before/wrong/hDCAxyVsPtantiHeliumTrueTransport"), antihePt, track.dcaXY());
                    if (isWeakDecay) {
                      histos.fill(HIST("tracks/helium/dca/before/wrong/hDCAxyVsPtantiHeliumTrueSec"), antihePt, track.dcaXY());
                    }
                    if (track.hasTOF() && outFlagOptions.doTOFplots) {
                      histos.fill(HIST("tracks/helium/dca/before/wrong/TOF/hDCAxyVsPtantiHeliumTrueTransport"), antihePt, track.dcaXY());
                      if (isWeakDecay) {
                        histos.fill(HIST("tracks/helium/dca/before/wrong/TOF/hDCAxyVsPtantiHeliumTrueSec"), antihePt, track.dcaXY());
                      }
                    }
                  }
                }
              }
            }
            break;
          case PDGAlpha:
            if (enableAl && alRapCut && outFlagOptions.makeDCABeforeCutPlots && passDCAzCut) {
              histos.fill(HIST("tracks/alpha/dca/before/hDCAxyVsPtAlphaTrue"), track.pt(), track.dcaXY());
              if (isPhysPrim) {
                histos.fill(HIST("tracks/alpha/dca/before/hDCAxyVsPtAlphaTruePrim"), track.pt(), track.dcaXY());
              }
              if (!isPhysPrim && !isProdByGen) {
                if (isWeakDecay) {
                  histos.fill(HIST("tracks/alpha/dca/before/hDCAxyVsPtAlphaTrueSec"), track.pt(), track.dcaXY());
                } else {
                  histos.fill(HIST("tracks/alpha/dca/before/hDCAxyVsPtAlphaTrueMaterial"), track.pt(), track.dcaXY());
                }
              }
            }
            break;
          case -PDGAlpha:
            if (enableAl && alRapCut && outFlagOptions.makeDCABeforeCutPlots && passDCAzCut) {
              histos.fill(HIST("tracks/alpha/dca/before/hDCAxyVsPtantiAlphaTrue"), track.pt(), track.dcaXY());
              if (isPhysPrim) {
                histos.fill(HIST("tracks/alpha/dca/before/hDCAxyVsPtantiAlphaTruePrim"), track.pt(), track.dcaXY());
              }
              if (!isPhysPrim && !isProdByGen) {
                if (isWeakDecay) {
                  histos.fill(HIST("tracks/alpha/dca/before/hDCAxyVsPtantiAlphaTrueSec"), track.pt(), track.dcaXY());
                } else {
                  histos.fill(HIST("tracks/alpha/dca/before/hDCAxyVsPtantiAlphaTrueMaterial"), track.pt(), track.dcaXY());
                }
              }
            }
            break;
          default:
            break;
        }
        switch (std::abs(pdgCode)) {
          case PDGDeuteron:
            //
            break;
          default:
            if (isDeWoDCAxyWTPCpid && outFlagOptions.makeFakeTracksPlots) {
              if (outFlagOptions.makeDCABeforeCutPlots) {
                histos.fill(HIST("tracks/deuteron/dca/before/fake/hDCAxyVsPtDeuteronTrue"), DPt, track.dcaXY());
                if (track.hasTOF() && outFlagOptions.doTOFplots) {
                  histos.fill(HIST("tracks/deuteron/dca/before/fake/TOF/hDCAxyVsPtDeuteronTrue"), DPt, track.dcaXY());
                }
              }
              if (isPhysPrim) {
                if (outFlagOptions.makeDCABeforeCutPlots) {
                  histos.fill(HIST("tracks/deuteron/dca/before/fake/hDCAxyVsPtDeuteronTruePrim"), DPt, track.dcaXY());
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    histos.fill(HIST("tracks/deuteron/dca/before/fake/TOF/hDCAxyVsPtDeuteronTruePrim"), DPt, track.dcaXY());
                  }
                }
              }
              if (!isPhysPrim && !isProdByGen) {
                if (outFlagOptions.makeDCABeforeCutPlots) {
                  if (isWeakDecay) {
                    histos.fill(HIST("tracks/deuteron/dca/before/fake/hDCAxyVsPtDeuteronTrueSec"), DPt, track.dcaXY());
                  } else {
                    histos.fill(HIST("tracks/deuteron/dca/before/fake/hDCAxyVsPtDeuteronTrueMaterial"), DPt, track.dcaXY());
                  }
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    if (isWeakDecay) {
                      histos.fill(HIST("tracks/deuteron/dca/before/fake/TOF/hDCAxyVsPtDeuteronTrueSec"), DPt, track.dcaXY());
                    } else {
                      histos.fill(HIST("tracks/deuteron/dca/before/fake/TOF/hDCAxyVsPtDeuteronTrueMaterial"), DPt, track.dcaXY());
                    }
                  }
                }
              }
            }
            if (isAntiDeWoDCAxyWTPCpid && outFlagOptions.makeFakeTracksPlots) {
              if (outFlagOptions.makeDCABeforeCutPlots) {
                histos.fill(HIST("tracks/deuteron/dca/before/fake/hDCAxyVsPtantiDeuteronTrue"), antiDPt, track.dcaXY());
                if (track.hasTOF() && outFlagOptions.doTOFplots) {
                  histos.fill(HIST("tracks/deuteron/dca/before/fake/TOF/hDCAxyVsPtantiDeuteronTrue"), antiDPt, track.dcaXY());
                }
              }
              if (isPhysPrim) {
                if (outFlagOptions.makeDCABeforeCutPlots) {
                  histos.fill(HIST("tracks/deuteron/dca/before/fake/hDCAxyVsPtantiDeuteronTruePrim"), antiDPt, track.dcaXY());
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    histos.fill(HIST("tracks/deuteron/dca/before/fake/TOF/hDCAxyVsPtantiDeuteronTruePrim"), antiDPt, track.dcaXY());
                  }
                }
              }
              if (!isPhysPrim && !isProdByGen) {
                if (outFlagOptions.makeDCABeforeCutPlots) {
                  if (isWeakDecay) {
                    histos.fill(HIST("tracks/deuteron/dca/before/fake/hDCAxyVsPtantiDeuteronTrueSec"), antiDPt, track.dcaXY());
                  } else {
                    histos.fill(HIST("tracks/deuteron/dca/before/fake/hDCAxyVsPtantiDeuteronTrueMaterial"), antiDPt, track.dcaXY());
                  }
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    if (isWeakDecay) {
                      histos.fill(HIST("tracks/deuteron/dca/before/fake/TOF/hDCAxyVsPtantiDeuteronTrueSec"), antiDPt, track.dcaXY());
                    } else {
                      histos.fill(HIST("tracks/deuteron/dca/before/fake/TOF/hDCAxyVsPtantiDeuteronTrueMaterial"), antiDPt, track.dcaXY());
                    }
                  }
                }
              }
            }
            break;
        }

        switch (std::abs(pdgCode)) {
          case PDGHelium:
            //
            break;
          default:
            if (isHeWoDCAxyWTPCpid && outFlagOptions.makeFakeTracksPlots) {
              if (outFlagOptions.makeDCABeforeCutPlots) {
                histos.fill(HIST("tracks/helium/dca/before/fake/hDCAxyVsPtHeliumTrue"), hePt, track.dcaXY());
                if (track.hasTOF() && outFlagOptions.doTOFplots) {
                  histos.fill(HIST("tracks/helium/dca/before/fake/TOF/hDCAxyVsPtHeliumTrue"), hePt, track.dcaXY());
                }
                if (isPhysPrim) {
                  histos.fill(HIST("tracks/helium/dca/before/fake/hDCAxyVsPtHeliumTruePrim"), hePt, track.dcaXY());
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    histos.fill(HIST("tracks/helium/dca/before/fake/TOF/hDCAxyVsPtHeliumTruePrim"), hePt, track.dcaXY());
                  }
                }
                if (!isPhysPrim && !isProdByGen) {
                  if (isWeakDecay) {
                    histos.fill(HIST("tracks/helium/dca/before/fake/hDCAxyVsPtHeliumTrueSec"), hePt, track.dcaXY());
                  } else {
                    histos.fill(HIST("tracks/helium/dca/before/fake/hDCAxyVsPtHeliumTrueMaterial"), hePt, track.dcaXY());
                  }
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    if (isWeakDecay) {
                      histos.fill(HIST("tracks/helium/dca/before/fake/TOF/hDCAxyVsPtHeliumTrueSec"), hePt, track.dcaXY());
                    } else {
                      histos.fill(HIST("tracks/helium/dca/before/fake/TOF/hDCAxyVsPtHeliumTrueMaterial"), hePt, track.dcaXY());
                    }
                  }
                }
              }
            }
            if (isAntiHeWoDCAxyWTPCpid && outFlagOptions.makeFakeTracksPlots) {
              if (outFlagOptions.makeDCABeforeCutPlots) {
                histos.fill(HIST("tracks/helium/dca/before/fake/hDCAxyVsPtantiHeliumTrue"), antihePt, track.dcaXY());
                if (track.hasTOF() && outFlagOptions.doTOFplots) {
                  histos.fill(HIST("tracks/helium/dca/before/fake/TOF/hDCAxyVsPtantiHeliumTrue"), antihePt, track.dcaXY());
                }
                if (isPhysPrim) {
                  histos.fill(HIST("tracks/helium/dca/before/fake/hDCAxyVsPtantiHeliumTruePrim"), antihePt, track.dcaXY());
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    histos.fill(HIST("tracks/helium/dca/before/fake/TOF/hDCAxyVsPtantiHeliumTruePrim"), antihePt, track.dcaXY());
                  }
                }
                if (!isPhysPrim && !isProdByGen) {
                  if (isWeakDecay) {
                    histos.fill(HIST("tracks/helium/dca/before/fake/hDCAxyVsPtantiHeliumTrueSec"), antihePt, track.dcaXY());
                  } else {
                    histos.fill(HIST("tracks/helium/dca/before/fake/hDCAxyVsPtantiHeliumTrueMaterial"), antihePt, track.dcaXY());
                  }
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    if (isWeakDecay) {
                      histos.fill(HIST("tracks/helium/dca/before/fake/TOF/hDCAxyVsPtantiHeliumTrueSec"), antihePt, track.dcaXY());
                    } else {
                      histos.fill(HIST("tracks/helium/dca/before/fake/TOF/hDCAxyVsPtantiHeliumTrueMaterial"), antihePt, track.dcaXY());
                    }
                  }
                }
              }
            }
            break;
        }
      } else {
        (void)particles;
      }

      // DCA Cut
      if constexpr (!IsFilteredData) {
        if (filterOptions.enableIsGlobalTrack) {
          if (!enableCustomDCACut) {
            if (!track.isGlobalTrack())
              continue;
          } else {
            if (!track.isGlobalTrackWoDCA())
              continue;
          }
        }
      }

      if (isPVContributorCut && !track.isPVContributor()) {
        continue;
      }

      if (outFlagOptions.makeDCAAfterCutPlots) {
        if (isHeWTPCpid) {
          histos.fill(HIST("tracks/helium/dca/after/hDCAxyVsDCAzVsPtHelium"), track.dcaXY(), track.dcaZ(), hePt);
          histos.fill(HIST("tracks/helium/dca/after/h3DCAvsPtHelium"), track.dcaXY(), track.dcaZ(), hePt);
          if (track.hasTOF() && outFlagOptions.doTOFplots) {
            histos.fill(HIST("tracks/helium/dca/after/TOF/hDCAxyVsDCAzVsPtHelium"), track.dcaXY(), track.dcaZ(), hePt);
          }
        }
        if (isAntiHeWTPCpid) {
          histos.fill(HIST("tracks/helium/dca/after/hDCAxyVsDCAzVsPtantiHelium"), track.dcaXY(), track.dcaZ(), antihePt);
          histos.fill(HIST("tracks/helium/dca/after/h3DCAvsPtantiHelium"), track.dcaXY(), track.dcaZ(), antihePt);
          if (track.hasTOF() && outFlagOptions.doTOFplots) {
            histos.fill(HIST("tracks/helium/dca/after/TOF/hDCAxyVsDCAzVsPtantiHelium"), track.dcaXY(), track.dcaZ(), antihePt);
          }
        }

        if (passDCAxyzCut) {
          histos.fill(HIST("tracks/dca/after/hDCAxy"), track.dcaXY());
          histos.fill(HIST("tracks/dca/after/hDCAz"), track.dcaZ());
          histos.fill(HIST("tracks/dca/after/hDCAxyVsPt"), track.pt(), track.dcaXY());
          histos.fill(HIST("tracks/dca/after/hDCAzVsPt"), track.pt(), track.dcaZ());

          if (enablePr && prRapCut && (std::abs(track.tpcNSigmaPr()) < nsigmaTPCvar.nsigmaTPCPr)) {
            if (track.sign() > 0) {
              histos.fill(HIST("tracks/proton/dca/after/hDCAxyVsPtProton"), track.pt(), track.dcaXY());
              histos.fill(HIST("tracks/proton/dca/after/hDCAzVsPtProton"), track.pt(), track.dcaZ());
            } else {
              histos.fill(HIST("tracks/proton/dca/after/hDCAxyVsPtantiProton"), track.pt(), track.dcaXY());
              histos.fill(HIST("tracks/proton/dca/after/hDCAzVsPtantiProton"), track.pt(), track.dcaZ());
            }
          }
          if (enableTr && trRapCut && isTritonTPCpid) {
            if (track.sign() > 0) {
              histos.fill(HIST("tracks/triton/dca/after/hDCAxyVsPtTriton"), track.pt(), track.dcaXY());
              histos.fill(HIST("tracks/triton/dca/after/hDCAzVsPtTriton"), track.pt(), track.dcaZ());
            } else {
              histos.fill(HIST("tracks/triton/dca/after/hDCAxyVsPtantiTriton"), track.pt(), track.dcaXY());
              histos.fill(HIST("tracks/triton/dca/after/hDCAzVsPtantiTriton"), track.pt(), track.dcaZ());
            }
          }
          if (enableAl && alRapCut && (std::abs(track.tpcNSigmaAl()) < nsigmaTPCvar.nsigmaTPCAl)) {
            if (track.sign() > 0) {
              histos.fill(HIST("tracks/alpha/dca/after/hDCAxyVsPtAlpha"), track.pt(), track.dcaXY());
              histos.fill(HIST("tracks/alpha/dca/after/hDCAzVsPtAlpha"), track.pt(), track.dcaZ());
            } else {
              histos.fill(HIST("tracks/alpha/dca/after/hDCAxyVsPtantiAlpha"), track.pt(), track.dcaXY());
              histos.fill(HIST("tracks/alpha/dca/after/hDCAzVsPtantiAlpha"), track.pt(), track.dcaZ());
            }
          }
        }
        if (isDeWTPCpid) {
          if (usenITSLayer && !itsClusterMap.test(trkqcOptions.nITSLayer))
            continue;
          histos.fill(HIST("tracks/deuteron/dca/after/hDCAxyVsPtDeuteron"), DPt, track.dcaXY());
          histos.fill(HIST("tracks/deuteron/dca/after/hDCAzVsPtDeuteron"), DPt, track.dcaZ());
        }
        if (isAntiDeWTPCpid) {
          if (usenITSLayer && !itsClusterMap.test(trkqcOptions.nITSLayer))
            continue;
          histos.fill(HIST("tracks/deuteron/dca/after/hDCAxyVsPtantiDeuteron"), antiDPt, track.dcaXY());
          histos.fill(HIST("tracks/deuteron/dca/after/hDCAzVsPtantiDeuteron"), antiDPt, track.dcaZ());
        }
        if (isHeWTPCpid) {
          histos.fill(HIST("tracks/helium/dca/after/hDCAxyVsPtHelium"), hePt, track.dcaXY());
          histos.fill(HIST("tracks/helium/dca/after/hDCAzVsPtHelium"), hePt, track.dcaZ());
          if (track.hasTOF() && outFlagOptions.doTOFplots) {
            histos.fill(HIST("tracks/helium/dca/after/TOF/hDCAxyVsPtHelium"), hePt, track.dcaXY());
            histos.fill(HIST("tracks/helium/dca/after/TOF/hDCAzVsPtHelium"), hePt, track.dcaZ());
          }
        }
        if (isAntiHeWTPCpid) {
          histos.fill(HIST("tracks/helium/dca/after/hDCAxyVsPtantiHelium"), antihePt, track.dcaXY());
          histos.fill(HIST("tracks/helium/dca/after/hDCAzVsPtantiHelium"), antihePt, track.dcaZ());
          if (track.hasTOF() && outFlagOptions.doTOFplots) {
            histos.fill(HIST("tracks/helium/dca/after/TOF/hDCAxyVsPtantiHelium"), antihePt, track.dcaXY());
            histos.fill(HIST("tracks/helium/dca/after/TOF/hDCAzVsPtantiHelium"), antihePt, track.dcaZ());
          }
        }
      }

      if (passDCAxyzCut) {
        // QA histos fill
        if (enableDebug) {
          histos.fill(HIST("qa/h1ITSncr"), track.itsNCls());
          histos.fill(HIST("qa/h1TPCnfound"), track.tpcNClsFound());
          histos.fill(HIST("qa/h1TPCncr"), track.tpcNClsCrossedRows());
          histos.fill(HIST("qa/h1rTPC"), track.tpcCrossedRowsOverFindableCls());
          histos.fill(HIST("qa/h1chi2ITS"), track.tpcChi2NCl());
          histos.fill(HIST("qa/h1chi2TPC"), track.itsChi2NCl());
          debugHistos.fill(HIST("debug/h2TPCsignVsTPCmomentum_AllTracks"), track.tpcInnerParam() / (1.f * track.sign()), track.tpcSignal());
        }

        if (enableDebug) {
          debugHistos.fill(HIST("debug/tracks/h1Eta"), track.eta());
          debugHistos.fill(HIST("debug/tracks/h1VarPhi"), track.phi());
          debugHistos.fill(HIST("debug/tracks/h2EtaVsPhi"), track.eta(), track.phi());
          debugHistos.fill(HIST("debug/tracks/h2PionYvsPt"), track.rapidity(o2::track::PID::getMass2Z(o2::track::PID::Pion)), track.pt());

          if (track.sign() > 0) {
            debugHistos.fill(HIST("debug/qa/h2TPCncrVsPtPos"), track.tpcInnerParam(), track.tpcNClsCrossedRows());
            debugHistos.fill(HIST("debug/qa/h2TPCncrVsTPCsignalPos"), track.tpcSignal(), track.tpcNClsCrossedRows());

            if (track.tpcInnerParam() < kCfgTpcClasses[0]) {
              debugHistos.fill(HIST("debug/qa/h1TPCncrLowPPos"), track.tpcNClsCrossedRows());
            }
            if ((track.tpcInnerParam() >= kCfgTpcClasses[0]) && (track.tpcInnerParam() < kCfgTpcClasses[1])) {
              debugHistos.fill(HIST("debug/qa/h1TPCncrMidPPos"), track.tpcNClsCrossedRows());
            }
            if (track.tpcInnerParam() >= kCfgTpcClasses[1]) {
              debugHistos.fill(HIST("debug/qa/h1TPCncrHighPPos"), track.tpcNClsCrossedRows());
            }
          } else {
            debugHistos.fill(HIST("debug/qa/h2TPCncrVsPtNeg"), track.tpcInnerParam(), track.tpcNClsCrossedRows());
            debugHistos.fill(HIST("debug/qa/h2TPCncrVsTPCsignalNeg"), track.tpcSignal(), track.tpcNClsCrossedRows());

            if (track.tpcInnerParam() < kCfgTpcClasses[0]) {
              debugHistos.fill(HIST("debug/qa/h1TPCncrLowPNeg"), track.tpcNClsCrossedRows());
            }
            if ((track.tpcInnerParam() >= kCfgTpcClasses[0]) && (track.tpcInnerParam() < kCfgTpcClasses[1])) {
              debugHistos.fill(HIST("debug/qa/h1TPCncrMidPNeg"), track.tpcNClsCrossedRows());
            }
            if (track.tpcInnerParam() >= kCfgTpcClasses[1]) {
              debugHistos.fill(HIST("debug/qa/h1TPCncrHighPNeg"), track.tpcNClsCrossedRows());
            }
          }

          debugHistos.fill(HIST("debug/tracks/pion/h2PionVspTNSigmaTPC"), track.pt(), track.tpcNSigmaPi());
          debugHistos.fill(HIST("debug/tracks/kaon/h2KaonVspTNSigmaTPC"), track.pt(), track.tpcNSigmaKa());
        }

        if (outFlagOptions.enableEffPlots)
          histos.fill(HIST("tracks/eff/h2pVsTPCmomentum"), track.tpcInnerParam(), track.p());

        if (filterOptions.enableFiltering) {
          if (track.tpcNSigmaKa() < kCfgKaonCut)
            continue;
        }

        if (outFlagOptions.enablePIDplot)
          histos.fill(HIST("tracks/h2TPCsignVsTPCmomentum"), track.tpcInnerParam() / (1.f * track.sign()), track.tpcSignal());

        if constexpr (!IsFilteredData) {
          if (nsigmaITSvar.showAverageClusterSize && outFlagOptions.enablePIDplot) {
            histos.fill(HIST("tracks/averageClusterSize"), track.p(), averageClusterSizeTrk(track));
            histos.fill(HIST("tracks/averageClusterSizePerCoslInv"), track.p(), averageClusterSizePerCoslInv(track));
          }
        }

        if (track.sign() > 0) {
          if (enablePr && prRapCut) {
            if (outFlagOptions.enableExpSignalTPC)
              histos.fill(HIST("tracks/proton/h2ProtonTPCExpSignalDiffVsPt"), track.pt(), track.tpcExpSignalDiffPr());

            switch (useHasTRDConfig) {
              case 0:
                histos.fill(HIST("tracks/proton/h2ProtonVspTNSigmaTPC"), track.pt(), track.tpcNSigmaPr());
                break;
              case 1:
                if (track.hasTRD()) {
                  histos.fill(HIST("tracks/proton/h2ProtonVspTNSigmaTPC"), track.pt(), track.tpcNSigmaPr());
                }
                break;
              case 2:
                if (!track.hasTRD()) {
                  histos.fill(HIST("tracks/proton/h2ProtonVspTNSigmaTPC"), track.pt(), track.tpcNSigmaPr());
                }
                break;
            }
          }
          if (enableTr && trRapCut) {
            histos.fill(HIST("tracks/triton/h2TritonVspTNSigmaTPC"), track.pt(), track.tpcNSigmaTr());
          }
          if (enableAl && alRapCut) {
            histos.fill(HIST("tracks/alpha/h2AlphaVspTNSigmaTPC"), track.pt(), track.tpcNSigmaAl());
          }
        } else {
          if (enablePr && prRapCut) {
            if (outFlagOptions.enableExpSignalTPC)
              histos.fill(HIST("tracks/proton/h2antiProtonTPCExpSignalDiffVsPt"), track.pt(), track.tpcExpSignalDiffPr());
            switch (useHasTRDConfig) {
              case 0:
                histos.fill(HIST("tracks/proton/h2antiProtonVspTNSigmaTPC"), track.pt(), track.tpcNSigmaPr());
                break;
              case 1:
                if (track.hasTRD()) {
                  histos.fill(HIST("tracks/proton/h2antiProtonVspTNSigmaTPC"), track.pt(), track.tpcNSigmaPr());
                }
                break;
              case 2:
                if (!track.hasTRD()) {
                  histos.fill(HIST("tracks/proton/h2antiProtonVspTNSigmaTPC"), track.pt(), track.tpcNSigmaPr());
                }
                break;
            }
          }
          if (enableTr && trRapCut) {
            histos.fill(HIST("tracks/triton/h2antiTritonVspTNSigmaTPC"), track.pt(), track.tpcNSigmaTr());
          }
          if (enableAl && alRapCut) {
            histos.fill(HIST("tracks/alpha/h2antiAlphaVspTNSigmaTPC"), track.pt(), track.tpcNSigmaAl());
          }
        }

        //  TOF
        if (outFlagOptions.doTOFplots) {
          if (enableDebug) {
            histos.fill(HIST("tracks/pion/h2PionVspTNSigmaTOF"), track.pt(), track.tofNSigmaPi());
            histos.fill(HIST("tracks/kaon/h2KaonVspTNSigmaTOF"), track.pt(), track.tofNSigmaKa());
          }
          if (track.sign() > 0) {
            if (enablePr && prRapCut) {

              switch (useHasTRDConfig) {
                case 0:
                  histos.fill(HIST("tracks/proton/h2ProtonVspTNSigmaTOF"), track.pt(), track.tofNSigmaPr());
                  break;
                case 1:
                  if (track.hasTRD()) {
                    histos.fill(HIST("tracks/proton/h2ProtonVspTNSigmaTOF"), track.pt(), track.tofNSigmaPr());
                  }
                  break;
                case 2:
                  if (!track.hasTRD()) {
                    histos.fill(HIST("tracks/proton/h2ProtonVspTNSigmaTOF"), track.pt(), track.tofNSigmaPr());
                  }
                  break;
              }
              if (outFlagOptions.enableExpSignalTOF)
                histos.fill(HIST("tracks/proton/h2ProtonTOFExpSignalDiffVsPt"), track.pt(), track.tofExpSignalDiffPr());
            }

            if (filterOptions.enableEvTimeSplitting && track.hasTOF()) {
              if (track.isEvTimeTOF() && track.isEvTimeT0AC()) {
                if (enablePr)
                  evtimeHistos.fill(HIST("tracks/evtime/ft0tof/proton/h2ProtonVspTNSigmaTOF"), track.pt(), track.tofNSigmaPr());
                if (enableDe)
                  evtimeHistos.fill(HIST("tracks/evtime/ft0tof/deuteron/h2DeuteronVspTNSigmaTOF"), DPt, track.tofNSigmaDe());
                if (outFlagOptions.enableExpSignalTOF) {
                  if (enablePr)
                    evtimeHistos.fill(HIST("tracks/evtime/ft0tof/proton/h2ProtonTOFExpSignalDiffVsPt"), track.pt(), track.tofExpSignalDiffPr());
                  if (enableDe)
                    evtimeHistos.fill(HIST("tracks/evtime/ft0tof/deuteron/h2DeuteronTOFExpSignalDiffVsPt"), DPt, track.tofExpSignalDiffDe());
                }
                if (enablePr)
                  evtimeHistos.fill(HIST("tracks/evtime/ft0tof/proton/h3ProtonNSigmaTPCvsNSigmaTOFvsPt"), track.tpcNSigmaPr(), track.tofNSigmaPr(), track.pt());
                if (enableDe)
                  evtimeHistos.fill(HIST("tracks/evtime/ft0tof/deuteron/h3DeuteronNSigmaTPCvsNSigmaTOFvsPt"), track.tpcNSigmaDe(), track.tofNSigmaDe(), DPt);
                if (enableDebug && (track.beta() > cfgBetaCut)) {
                  if (enablePr)
                    debugHistos.fill(HIST("debug/evtime/ft0tof/proton/h2ProtonVspTNSigmaTPC_BetaCut"), track.pt(), track.tpcNSigmaPr());
                  if (enableDe)
                    debugHistos.fill(HIST("debug/evtime/ft0tof/deuteron/h2DeuteronVspTNSigmaTPC_BetaCut"), DPt, track.tpcNSigmaDe());

                  if (enablePr)
                    debugHistos.fill(HIST("debug/evtime/ft0tof/proton/h2ProtonVspTNSigmaTOF_BetaCut"), track.pt(), track.tofNSigmaPr());
                  if (enableDe)
                    debugHistos.fill(HIST("debug/evtime/ft0tof/deuteron/h2DeuteronVspTNSigmaTOF_BetaCut"), DPt, track.tofNSigmaDe());
                  if (outFlagOptions.enableExpSignalTOF) {
                    debugHistos.fill(HIST("debug/evtime/ft0tof/proton/h2ProtonTOFExpSignalDiffVsPt_BetaCut"), track.pt(), track.tofExpSignalDiffPr());
                    if (enableDe)
                      debugHistos.fill(HIST("debug/evtime/ft0tof/deuteron/h2DeuteronTOFExpSignalDiffVsPt_BetaCut"), DPt, track.tofExpSignalDiffDe());
                  }
                }
              } else if (track.isEvTimeT0AC()) {
                if (enablePr)
                  evtimeHistos.fill(HIST("tracks/evtime/ft0/proton/h2ProtonVspTNSigmaTOF"), track.pt(), track.tofNSigmaPr());
                if (enableDe)
                  evtimeHistos.fill(HIST("tracks/evtime/ft0/deuteron/h2DeuteronVspTNSigmaTOF"), DPt, track.tofNSigmaDe());
                if (outFlagOptions.enableExpSignalTOF) {
                  if (enablePr)
                    evtimeHistos.fill(HIST("tracks/evtime/ft0/proton/h2ProtonTOFExpSignalDiffVsPt"), track.pt(), track.tofExpSignalDiffPr());
                  if (enableDe)
                    evtimeHistos.fill(HIST("tracks/evtime/ft0/deuteron/h2DeuteronTOFExpSignalDiffVsPt"), DPt, track.tofExpSignalDiffDe());
                }
                if (enablePr)
                  evtimeHistos.fill(HIST("tracks/evtime/ft0/proton/h3ProtonNSigmaTPCvsNSigmaTOFvsPt"), track.tpcNSigmaPr(), track.tofNSigmaPr(), track.pt());
                if (enableDe)
                  evtimeHistos.fill(HIST("tracks/evtime/ft0/deuteron/h3DeuteronNSigmaTPCvsNSigmaTOFvsPt"), track.tpcNSigmaDe(), track.tofNSigmaDe(), DPt);
                if (enableDebug && (track.beta() > cfgBetaCut)) {
                  if (enablePr)
                    debugHistos.fill(HIST("debug/evtime/ft0/proton/h2ProtonVspTNSigmaTPC_BetaCut"), track.pt(), track.tpcNSigmaPr());
                  if (enableDe)
                    debugHistos.fill(HIST("debug/evtime/ft0/deuteron/h2DeuteronVspTNSigmaTPC_BetaCut"), DPt, track.tpcNSigmaDe());
                  if (enablePr)
                    debugHistos.fill(HIST("debug/evtime/ft0/proton/h2ProtonVspTNSigmaTOF_BetaCut"), track.pt(), track.tofNSigmaPr());
                  if (enableDe)
                    debugHistos.fill(HIST("debug/evtime/ft0/deuteron/h2DeuteronVspTNSigmaTOF_BetaCut"), DPt, track.tofNSigmaDe());
                  if (outFlagOptions.enableExpSignalTOF) {
                    if (enablePr)
                      debugHistos.fill(HIST("debug/evtime/ft0/proton/h2ProtonTOFExpSignalDiffVsPt_BetaCut"), track.pt(), track.tofExpSignalDiffPr());
                    if (enableDe)
                      debugHistos.fill(HIST("debug/evtime/ft0/deuteron/h2DeuteronTOFExpSignalDiffVsPt_BetaCut"), DPt, track.tofExpSignalDiffDe());
                  }
                }
              } else if (track.isEvTimeTOF()) {
                if (enablePr)
                  evtimeHistos.fill(HIST("tracks/evtime/tof/proton/h2ProtonVspTNSigmaTOF"), track.pt(), track.tofNSigmaPr());
                if (enableDe)
                  evtimeHistos.fill(HIST("tracks/evtime/tof/deuteron/h2DeuteronVspTNSigmaTOF"), DPt, track.tofNSigmaDe());
                if (outFlagOptions.enableExpSignalTOF) {
                  if (enablePr)
                    evtimeHistos.fill(HIST("tracks/evtime/tof/proton/h2ProtonTOFExpSignalDiffVsPt"), track.pt(), track.tofExpSignalDiffPr());
                  if (enableDe)
                    evtimeHistos.fill(HIST("tracks/evtime/tof/deuteron/h2DeuteronTOFExpSignalDiffVsPt"), DPt, track.tofExpSignalDiffDe());
                }
                if (enablePr)
                  evtimeHistos.fill(HIST("tracks/evtime/tof/proton/h3ProtonNSigmaTPCvsNSigmaTOFvsPt"), track.tpcNSigmaPr(), track.tofNSigmaPr(), track.pt());
                if (enableDe)
                  evtimeHistos.fill(HIST("tracks/evtime/tof/deuteron/h3DeuteronNSigmaTPCvsNSigmaTOFvsPt"), track.tpcNSigmaDe(), track.tofNSigmaDe(), DPt);
                if (enableDebug && (track.beta() > cfgBetaCut)) {
                  if (enablePr)
                    debugHistos.fill(HIST("debug/evtime/tof/proton/h2ProtonVspTNSigmaTPC_BetaCut"), track.pt(), track.tpcNSigmaPr());
                  if (enableDe)
                    debugHistos.fill(HIST("debug/evtime/tof/deuteron/h2DeuteronVspTNSigmaTPC_BetaCut"), DPt, track.tpcNSigmaDe());

                  if (enablePr)
                    debugHistos.fill(HIST("debug/evtime/tof/proton/h2ProtonVspTNSigmaTOF_BetaCut"), track.pt(), track.tofNSigmaPr());
                  if (enableDe)
                    debugHistos.fill(HIST("debug/evtime/tof/deuteron/h2DeuteronVspTNSigmaTOF_BetaCut"), DPt, track.tofNSigmaDe());

                  if (outFlagOptions.enableExpSignalTOF) {
                    if (enablePr)
                      debugHistos.fill(HIST("debug/evtime/tof/proton/h2ProtonTOFExpSignalDiffVsPt_BetaCut"), track.pt(), track.tofExpSignalDiffPr());
                    if (enableDe)
                      debugHistos.fill(HIST("debug/evtime/tof/deuteron/h2DeuteronTOFExpSignalDiffVsPt_BetaCut"), DPt, track.tofExpSignalDiffDe());
                  }
                }
              } else {
                if (enablePr)
                  evtimeHistos.fill(HIST("tracks/evtime/fill/proton/h2ProtonVspTNSigmaTOF"), track.pt(), track.tofNSigmaPr());
                if (enableDe)
                  evtimeHistos.fill(HIST("tracks/evtime/fill/deuteron/h2DeuteronVspTNSigmaTOF"), DPt, track.tofNSigmaDe());
                if (outFlagOptions.enableExpSignalTOF) {
                  if (enablePr)
                    evtimeHistos.fill(HIST("tracks/evtime/fill/proton/h2ProtonTOFExpSignalDiffVsPt"), track.pt(), track.tofExpSignalDiffPr());
                  if (enableDe)
                    evtimeHistos.fill(HIST("tracks/evtime/fill/deuteron/h2DeuteronTOFExpSignalDiffVsPt"), DPt, track.tofExpSignalDiffDe());
                }
                if (enablePr)
                  evtimeHistos.fill(HIST("tracks/evtime/fill/proton/h3ProtonNSigmaTPCvsNSigmaTOFvsPt"), track.tpcNSigmaPr(), track.tofNSigmaPr(), track.pt());
                if (enableDe)
                  evtimeHistos.fill(HIST("tracks/evtime/fill/deuteron/h3DeuteronNSigmaTPCvsNSigmaTOFvsPt"), track.tpcNSigmaDe(), track.tofNSigmaDe(), DPt);
                if (enableDebug && (track.beta() > cfgBetaCut)) {
                  if (enablePr)
                    debugHistos.fill(HIST("debug/evtime/fill/proton/h2ProtonVspTNSigmaTPC_BetaCut"), track.pt(), track.tpcNSigmaPr());
                  if (enableDe)
                    debugHistos.fill(HIST("debug/evtime/fill/deuteron/h2DeuteronVspTNSigmaTPC_BetaCut"), DPt, track.tpcNSigmaDe());

                  if (enablePr)
                    debugHistos.fill(HIST("debug/evtime/fill/proton/h2ProtonVspTNSigmaTOF_BetaCut"), track.pt(), track.tofNSigmaPr());
                  if (enableDe)
                    debugHistos.fill(HIST("debug/evtime/fill/deuteron/h2DeuteronVspTNSigmaTOF_BetaCut"), DPt, track.tofNSigmaDe());

                  if (outFlagOptions.enableExpSignalTOF) {
                    if (enablePr)
                      debugHistos.fill(HIST("debug/evtime/fill/proton/h2ProtonTOFExpSignalDiffVsPt_BetaCut"), track.pt(), track.tofExpSignalDiffPr());
                    if (enableDe)
                      debugHistos.fill(HIST("debug/evtime/fill/deuteron/h2DeuteronTOFExpSignalDiffVsPt_BetaCut"), DPt, track.tofExpSignalDiffDe());
                  }
                }
              }
            }
          } else {
            if (enablePr && prRapCut) {
              switch (useHasTRDConfig) {
                case 0:
                  histos.fill(HIST("tracks/proton/h2antiProtonVspTNSigmaTOF"), track.pt(), track.tofNSigmaPr());
                  break;
                case 1:
                  if (track.hasTRD()) {
                    histos.fill(HIST("tracks/proton/h2antiProtonVspTNSigmaTOF"), track.pt(), track.tofNSigmaPr());
                  }
                  break;
                case 2:
                  if (!track.hasTRD()) {
                    histos.fill(HIST("tracks/proton/h2antiProtonVspTNSigmaTOF"), track.pt(), track.tofNSigmaPr());
                  }
                  break;
              }
              if (outFlagOptions.enableExpSignalTOF)
                histos.fill(HIST("tracks/proton/h2antiProtonTOFExpSignalDiffVsPt"), track.pt(), track.tofExpSignalDiffPr());
            }
            if (filterOptions.enableEvTimeSplitting && track.hasTOF()) {
              if (track.isEvTimeTOF() && track.isEvTimeT0AC()) {
                if (enablePr)
                  evtimeHistos.fill(HIST("tracks/evtime/ft0tof/proton/h2antiProtonVspTNSigmaTOF"), track.pt(), track.tofNSigmaPr());
                if (enableDe)
                  evtimeHistos.fill(HIST("tracks/evtime/ft0tof/deuteron/h2antiDeuteronVspTNSigmaTOF"), antiDPt, track.tofNSigmaDe());
                if (outFlagOptions.enableExpSignalTOF) {
                  if (enablePr)
                    evtimeHistos.fill(HIST("tracks/evtime/ft0tof/proton/h2antiProtonTOFExpSignalDiffVsPt"), track.pt(), track.tofExpSignalDiffPr());
                  if (enableDe)
                    evtimeHistos.fill(HIST("tracks/evtime/ft0tof/deuteron/h2antiDeuteronTOFExpSignalDiffVsPt"), antiDPt, track.tofExpSignalDiffDe());
                }
                if (enablePr)
                  evtimeHistos.fill(HIST("tracks/evtime/ft0tof/proton/h3antiProtonNSigmaTPCvsNSigmaTOFvsPt"), track.tpcNSigmaPr(), track.tofNSigmaPr(), track.pt());
                if (enableDe)
                  evtimeHistos.fill(HIST("tracks/evtime/ft0tof/deuteron/h3antiDeuteronNSigmaTPCvsNSigmaTOFvsPt"), track.tpcNSigmaDe(), track.tofNSigmaDe(), antiDPt);
                if (enableDebug && (track.beta() > cfgBetaCut)) {
                  if (enablePr)
                    debugHistos.fill(HIST("debug/evtime/ft0tof/proton/h2antiProtonVspTNSigmaTPC_BetaCut"), track.pt(), track.tpcNSigmaPr());
                  if (enableDe)
                    debugHistos.fill(HIST("debug/evtime/ft0tof/deuteron/h2antiDeuteronVspTNSigmaTPC_BetaCut"), antiDPt, track.tpcNSigmaDe());
                  if (enablePr)
                    debugHistos.fill(HIST("debug/evtime/ft0tof/proton/h2antiProtonVspTNSigmaTOF_BetaCut"), track.pt(), track.tofNSigmaPr());
                  if (enableDe)
                    debugHistos.fill(HIST("debug/evtime/ft0tof/deuteron/h2antiDeuteronVspTNSigmaTOF_BetaCut"), antiDPt, track.tofNSigmaDe());

                  if (outFlagOptions.enableExpSignalTOF) {
                    if (enablePr)
                      debugHistos.fill(HIST("debug/evtime/ft0tof/proton/h2antiProtonTOFExpSignalDiffVsPt_BetaCut"), track.pt(), track.tofExpSignalDiffPr());
                    if (enableDe)
                      debugHistos.fill(HIST("debug/evtime/ft0tof/deuteron/h2antiDeuteronTOFExpSignalDiffVsPt_BetaCut"), antiDPt, track.tofExpSignalDiffDe());
                  }
                }
              } else if (track.isEvTimeT0AC()) {
                if (enablePr)
                  evtimeHistos.fill(HIST("tracks/evtime/ft0/proton/h2antiProtonVspTNSigmaTOF"), track.pt(), track.tofNSigmaPr());
                if (enableDe)
                  evtimeHistos.fill(HIST("tracks/evtime/ft0/deuteron/h2antiDeuteronVspTNSigmaTOF"), antiDPt, track.tofNSigmaDe());
                if (outFlagOptions.enableExpSignalTOF) {
                  if (enablePr)
                    evtimeHistos.fill(HIST("tracks/evtime/ft0/proton/h2antiProtonTOFExpSignalDiffVsPt"), track.pt(), track.tofExpSignalDiffPr());
                  if (enableDe)
                    evtimeHistos.fill(HIST("tracks/evtime/ft0/deuteron/h2antiDeuteronTOFExpSignalDiffVsPt"), antiDPt, track.tofExpSignalDiffDe());
                  if (enablePr)
                    evtimeHistos.fill(HIST("tracks/evtime/ft0/proton/h3antiProtonNSigmaTPCvsNSigmaTOFvsPt"), track.tpcNSigmaPr(), track.tofNSigmaPr(), track.pt());
                  if (enablePr)
                    debugHistos.fill(HIST("debug/evtime/ft0/proton/h2antiProtonVspTNSigmaTPC_BetaCut"), track.pt(), track.tpcNSigmaPr());
                  if (enableDe)
                    debugHistos.fill(HIST("debug/evtime/ft0/deuteron/h2antiDeuteronVspTNSigmaTPC_BetaCut"), antiDPt, track.tpcNSigmaDe());
                  if (enablePr)
                    debugHistos.fill(HIST("debug/evtime/ft0/proton/h2antiProtonVspTNSigmaTOF_BetaCut"), track.pt(), track.tofNSigmaPr());
                  if (enableDe)
                    debugHistos.fill(HIST("debug/evtime/ft0/deuteron/h2antiDeuteronVspTNSigmaTOF_BetaCut"), antiDPt, track.tofNSigmaDe());

                  if (outFlagOptions.enableExpSignalTOF) {
                    if (enablePr)
                      debugHistos.fill(HIST("debug/evtime/ft0/proton/h2antiProtonTOFExpSignalDiffVsPt_BetaCut"), track.pt(), track.tofExpSignalDiffPr());
                    if (enableDe)
                      debugHistos.fill(HIST("debug/evtime/ft0/deuteron/h2antiDeuteronTOFExpSignalDiffVsPt_BetaCut"), antiDPt, track.tofExpSignalDiffDe());
                  }
                }
              } else if (track.isEvTimeTOF()) {
                if (enablePr)
                  evtimeHistos.fill(HIST("tracks/evtime/tof/proton/h2antiProtonVspTNSigmaTOF"), track.pt(), track.tofNSigmaPr());
                if (enableDe)
                  evtimeHistos.fill(HIST("tracks/evtime/tof/deuteron/h2antiDeuteronVspTNSigmaTOF"), antiDPt, track.tofNSigmaDe());
                if (outFlagOptions.enableExpSignalTOF) {
                  if (enablePr)
                    evtimeHistos.fill(HIST("tracks/evtime/tof/proton/h2antiProtonTOFExpSignalDiffVsPt"), track.pt(), track.tofExpSignalDiffPr());
                  if (enableDe)
                    evtimeHistos.fill(HIST("tracks/evtime/tof/deuteron/h2antiDeuteronTOFExpSignalDiffVsPt"), antiDPt, track.tofExpSignalDiffDe());
                  if (enablePr)
                    evtimeHistos.fill(HIST("tracks/evtime/tof/proton/h3antiProtonNSigmaTPCvsNSigmaTOFvsPt"), track.tpcNSigmaPr(), track.tofNSigmaPr(), track.pt());
                  if (enablePr)
                    debugHistos.fill(HIST("debug/evtime/tof/proton/h2antiProtonVspTNSigmaTPC_BetaCut"), track.pt(), track.tpcNSigmaPr());
                  if (enableDe)
                    debugHistos.fill(HIST("debug/evtime/tof/deuteron/h2antiDeuteronVspTNSigmaTPC_BetaCut"), antiDPt, track.tpcNSigmaDe());

                  if (enablePr)
                    debugHistos.fill(HIST("debug/evtime/tof/proton/h2antiProtonVspTNSigmaTOF_BetaCut"), track.pt(), track.tofNSigmaPr());
                  if (enableDe)
                    debugHistos.fill(HIST("debug/evtime/tof/deuteron/h2antiDeuteronVspTNSigmaTOF_BetaCut"), antiDPt, track.tofNSigmaDe());

                  if (outFlagOptions.enableExpSignalTOF) {
                    if (enablePr)
                      debugHistos.fill(HIST("debug/evtime/tof/proton/h2antiProtonTOFExpSignalDiffVsPt_BetaCut"), track.pt(), track.tofExpSignalDiffPr());
                    if (enableDe)
                      debugHistos.fill(HIST("debug/evtime/tof/deuteron/h2antiDeuteronTOFExpSignalDiffVsPt_BetaCut"), antiDPt, track.tofExpSignalDiffDe());
                  }
                }
              } else {
                if (enablePr)
                  evtimeHistos.fill(HIST("tracks/evtime/fill/proton/h2antiProtonVspTNSigmaTOF"), track.pt(), track.tofNSigmaPr());
                if (enableDe)
                  evtimeHistos.fill(HIST("tracks/evtime/fill/deuteron/h2antiDeuteronVspTNSigmaTOF"), antiDPt, track.tofNSigmaDe());
                if (outFlagOptions.enableExpSignalTOF) {
                  if (enablePr)
                    evtimeHistos.fill(HIST("tracks/evtime/fill/proton/h2antiProtonTOFExpSignalDiffVsPt"), track.pt(), track.tofExpSignalDiffPr());
                  if (enableDe)
                    evtimeHistos.fill(HIST("tracks/evtime/fill/deuteron/h2antiDeuteronTOFExpSignalDiffVsPt"), antiDPt, track.tofExpSignalDiffDe());
                  if (enablePr)
                    evtimeHistos.fill(HIST("tracks/evtime/fill/proton/h3antiProtonNSigmaTPCvsNSigmaTOFvsPt"), track.tpcNSigmaPr(), track.tofNSigmaPr(), track.pt());
                  if (enablePr)
                    debugHistos.fill(HIST("debug/evtime/fill/proton/h2antiProtonVspTNSigmaTPC_BetaCut"), track.pt(), track.tpcNSigmaPr());
                  if (enableDe)
                    debugHistos.fill(HIST("debug/evtime/fill/deuteron/h2antiDeuteronVspTNSigmaTPC_BetaCut"), antiDPt, track.tpcNSigmaDe());

                  if (enablePr)
                    debugHistos.fill(HIST("debug/evtime/fill/proton/h2antiProtonVspTNSigmaTOF_BetaCut"), track.pt(), track.tofNSigmaPr());
                  if (enableDe)
                    debugHistos.fill(HIST("debug/evtime/fill/deuteron/h2antiDeuteronVspTNSigmaTOF_BetaCut"), antiDPt, track.tofNSigmaDe());

                  if (outFlagOptions.enableExpSignalTOF) {
                    if (enablePr)
                      debugHistos.fill(HIST("debug/evtime/fill/proton/h2antiProtonTOFExpSignalDiffVsPt_BetaCut"), track.pt(), track.tofExpSignalDiffPr());
                    if (enableDe)
                      debugHistos.fill(HIST("debug/evtime/fill/deuteron/h2antiDeuteronTOFExpSignalDiffVsPt_BetaCut"), antiDPt, track.tofExpSignalDiffDe());
                  }
                }
              }
            }
          }
        }
      }

      if (isDeWoTPCpid) {
        if (outFlagOptions.enableExpSignalTPC)
          histos.fill(HIST("tracks/deuteron/h2DeuteronTPCExpSignalDiffVsPt"), DPt, track.tpcExpSignalDiffDe());

        histos.fill(HIST("tracks/deuteron/h2DeuteronVspNSigmaITSDe"), track.p(), nITSDe);

        switch (useHasTRDConfig) {
          case 0:
            if (enableCentrality)
              histos.fill(HIST("tracks/deuteron/h3DeuteronVspTNSigmaTPCVsMult"), DPt, track.tpcNSigmaDe(), event.centFT0M());
            else
              histos.fill(HIST("tracks/deuteron/h2DeuteronVspTNSigmaTPC"), DPt, track.tpcNSigmaDe());
            break;
          case 1:
            if (track.hasTRD()) {
              histos.fill(HIST("tracks/deuteron/h2DeuteronVspTNSigmaTPC"), DPt, track.tpcNSigmaDe());
            }
            break;
          case 2:
            if (!track.hasTRD()) {
              histos.fill(HIST("tracks/deuteron/h2DeuteronVspTNSigmaTPC"), DPt, track.tpcNSigmaDe());
            }
            break;
        }
      }
      if (isAntiDeWoTPCpid) {
        if (outFlagOptions.enableExpSignalTPC)
          histos.fill(HIST("tracks/deuteron/h2antiDeuteronTPCExpSignalDiffVsPt"), antiDPt, track.tpcExpSignalDiffDe());

        histos.fill(HIST("tracks/deuteron/h2antiDeuteronVspNSigmaITSDe"), track.p(), nITSDe);

        switch (useHasTRDConfig) {
          case 0:
            if (enableCentrality)
              histos.fill(HIST("tracks/deuteron/h3antiDeuteronVspTNSigmaTPCVsMult"), antiDPt, track.tpcNSigmaDe(), event.centFT0M());
            else
              histos.fill(HIST("tracks/deuteron/h2antiDeuteronVspTNSigmaTPC"), antiDPt, track.tpcNSigmaDe());
            break;
          case 1:
            if (track.hasTRD()) {
              histos.fill(HIST("tracks/deuteron/h2antiDeuteronVspTNSigmaTPC"), antiDPt, track.tpcNSigmaDe());
            }
            break;
          case 2:
            if (!track.hasTRD()) {
              histos.fill(HIST("tracks/deuteron/h2antiDeuteronVspTNSigmaTPC"), antiDPt, track.tpcNSigmaDe());
            }
            break;
        }
      }

      if (isHeWoTPCpid) {
        if (outFlagOptions.enableExpSignalTPC)
          histos.fill(HIST("tracks/helium/h2HeliumTPCExpSignalDiffVsPt"), hePt, track.tpcExpSignalDiffHe());
        histos.fill(HIST("tracks/helium/h2HeliumVspTNSigmaITSHe"), track.p(), nITSHe);
        histos.fill(HIST("tracks/helium/h2HeliumVspTNSigmaTPC"), hePt, track.tpcNSigmaHe());
        if (enableCentrality)
          histos.fill(HIST("tracks/helium/h3HeliumVspTNSigmaTPCVsMult"), hePt, track.tpcNSigmaHe(), event.centFT0M());
      }
      if (isAntiHeWoTPCpid) {
        if (outFlagOptions.enableExpSignalTPC)
          histos.fill(HIST("tracks/helium/h2antiHeliumTPCExpSignalDiffVsPt"), antihePt, track.tpcExpSignalDiffHe());
        histos.fill(HIST("tracks/helium/h2antiHeliumVspTNSigmaITSHe"), track.p(), nITSHe);
        histos.fill(HIST("tracks/helium/h2antiHeliumVspTNSigmaTPC"), antihePt, track.tpcNSigmaHe());
        if (enableCentrality)
          histos.fill(HIST("tracks/helium/h3antiHeliumVspTNSigmaTPCVsMult"), antihePt, track.tpcNSigmaHe(), event.centFT0M());
      }
      if (isHeWTPCpid) {
        histos.fill(HIST("tracks/helium/h2HeliumVspTNSigmaITSHe_wTPCpid"), track.p(), nITSHe);
      }
      if (isAntiHeWTPCpid) {
        histos.fill(HIST("tracks/helium/h2antiHeliumVspTNSigmaITSHe_wTPCpid"), track.p(), nITSHe);
      }
      if constexpr (!IsFilteredData) {
        if (isHeWTPCpid || isAntiHeWTPCpid) {
          if (nsigmaITSvar.showAverageClusterSize) {
            histos.fill(HIST("tracks/helium/averageClusterSize"), track.p(), averageClusterSizeTrk(track));
            histos.fill(HIST("tracks/helium/averageClusterSizePerCoslInv"), track.p(), averageClusterSizePerCoslInv(track));
          }
        }
      }

      //  TOF
      if (outFlagOptions.doTOFplots) {

        if (isDeWoTPCpid) {
          switch (useHasTRDConfig) {
            case 0:
              histos.fill(HIST("tracks/deuteron/h2DeuteronVspTNSigmaTOF"), DPt, track.tofNSigmaDe());
              break;
            case 1:
              if (track.hasTRD()) {
                histos.fill(HIST("tracks/deuteron/h2DeuteronVspTNSigmaTOF"), DPt, track.tofNSigmaDe());
              }
              break;
            case 2:
              if (!track.hasTRD()) {
                histos.fill(HIST("tracks/deuteron/h2DeuteronVspTNSigmaTOF"), DPt, track.tofNSigmaDe());
              }
              break;
          }
          if (outFlagOptions.enableExpSignalTOF)
            histos.fill(HIST("tracks/deuteron/h2DeuteronTOFExpSignalDiffVsPt"), DPt, track.tofExpSignalDiffDe());
        }

        if (isAntiDeWoTPCpid) {
          switch (useHasTRDConfig) {
            case 0:
              histos.fill(HIST("tracks/deuteron/h2antiDeuteronVspTNSigmaTOF"), antiDPt, track.tofNSigmaDe());
              break;
            case 1:
              if (track.hasTRD()) {
                histos.fill(HIST("tracks/deuteron/h2antiDeuteronVspTNSigmaTOF"), antiDPt, track.tofNSigmaDe());
              }
              break;
            case 2:
              if (!track.hasTRD()) {
                histos.fill(HIST("tracks/deuteron/h2antiDeuteronVspTNSigmaTOF"), antiDPt, track.tofNSigmaDe());
              }
              break;
          }
          if (outFlagOptions.enableExpSignalTOF)
            histos.fill(HIST("tracks/deuteron/h2antiDeuteronTOFExpSignalDiffVsPt"), antiDPt, track.tofExpSignalDiffDe());
        }

        if (isHeWoTPCpid) {
          histos.fill(HIST("tracks/helium/h2HeliumVspTNSigmaTOF"), hePt, track.tofNSigmaHe());
          if (outFlagOptions.enableExpSignalTOF)
            histos.fill(HIST("tracks/helium/h2HeliumTOFExpSignalDiffVsPt"), hePt, track.tofExpSignalDiffHe());
        }

        if (isAntiHeWoTPCpid) {
          histos.fill(HIST("tracks/helium/h2antiHeliumVspTNSigmaTOF"), antihePt, track.tofNSigmaHe());
          if (outFlagOptions.enableExpSignalTOF)
            histos.fill(HIST("tracks/helium/h2antiHeliumTOFExpSignalDiffVsPt"), antihePt, track.tofExpSignalDiffHe());
        }
      }

      if (passDCAxyzCut) {
        // PID
        if (outFlagOptions.enableEffPlots && enableDebug) {
          if (track.sign() > 0)
            debugHistos.fill(HIST("tracks/eff/hPtP"), track.pt());
          else
            debugHistos.fill(HIST("tracks/eff/hPtantiP"), track.pt());
        }

        if (enablePr) {
          if (std::abs(track.tpcNSigmaPr()) < nsigmaTPCvar.nsigmaTPCPr && prRapCut) {
            if (track.sign() > 0) {
              if (outFlagOptions.enableEffPlots) {
                histos.fill(HIST("tracks/eff/proton/hPtPr"), track.pt());
                histos.fill(HIST("tracks/eff/proton/h2pVsTPCmomentumPr"), track.tpcInnerParam(), track.p());
              }
              histos.fill(HIST("tracks/proton/h1ProtonSpectra"), track.pt());
              histos.fill(HIST("tracks/proton/h2ProtonYvsPt"), track.rapidity(o2::track::PID::getMass2Z(o2::track::PID::Proton)), track.pt());
              histos.fill(HIST("tracks/proton/h2ProtonEtavsPt"), track.eta(), track.pt());

              if (outFlagOptions.enablePIDplot)
                histos.fill(HIST("tracks/proton/h2TPCsignVsTPCmomentumProton"), track.tpcInnerParam(), track.tpcSignal());
            } else {
              if (outFlagOptions.enableEffPlots) {
                histos.fill(HIST("tracks/eff/proton/hPtantiPr"), track.pt());
                histos.fill(HIST("tracks/eff/proton/h2pVsTPCmomentumantiPr"), track.tpcInnerParam(), track.p());
              }
              histos.fill(HIST("tracks/proton/h1antiProtonSpectra"), track.pt());
              histos.fill(HIST("tracks/proton/h2antiProtonYvsPt"), track.rapidity(o2::track::PID::getMass2Z(o2::track::PID::Proton)), track.pt());
              histos.fill(HIST("tracks/proton/h2antiProtonEtavsPt"), track.eta(), track.pt());

              if (outFlagOptions.enablePIDplot)
                histos.fill(HIST("tracks/proton/h2TPCsignVsTPCmomentumantiProton"), track.tpcInnerParam(), track.tpcSignal());
            }
          }
        }
        if (enableTr) {
          if ((isTritonTPCpid) && trRapCut) {
            if (track.sign() > 0) {
              if (outFlagOptions.enableEffPlots) {
                histos.fill(HIST("tracks/eff/triton/hPtTr"), track.pt());
                histos.fill(HIST("tracks/eff/triton/h2pVsTPCmomentumTr"), track.tpcInnerParam(), track.p());
              }
              histos.fill(HIST("tracks/triton/h1TritonSpectra"), track.pt());
              if (outFlagOptions.enablePIDplot)
                histos.fill(HIST("tracks/triton/h2TPCsignVsTPCmomentumTriton"), track.tpcInnerParam(), track.tpcSignal());
            } else {
              if (outFlagOptions.enableEffPlots) {
                histos.fill(HIST("tracks/eff/triton/hPtantiTr"), track.pt());
                histos.fill(HIST("tracks/eff/triton/h2pVsTPCmomentumantiTr"), track.tpcInnerParam(), track.p());
              }
              histos.fill(HIST("tracks/triton/h1antiTritonSpectra"), track.pt());
              if (outFlagOptions.enablePIDplot)
                histos.fill(HIST("tracks/triton/h2TPCsignVsTPCmomentumantiTriton"), track.tpcInnerParam(), track.tpcSignal());
            }
          }
        }
        if (enableAl) {
          if ((std::abs(track.tpcNSigmaAl()) < nsigmaTPCvar.nsigmaTPCAl) && alRapCut) {
            if (track.sign() > 0) {
              histos.fill(HIST("tracks/alpha/h1AlphaSpectra"), track.pt());
              if (outFlagOptions.enablePIDplot)
                histos.fill(HIST("tracks/alpha/h2TPCsignVsTPCmomentumAlpha"), track.tpcInnerParam(), track.tpcSignal());
            } else {
              histos.fill(HIST("tracks/alpha/h1antiAlphaSpectra"), track.pt());
              if (outFlagOptions.enablePIDplot)
                histos.fill(HIST("tracks/alpha/h2TPCsignVsTPCmomentumantiAlpha"), track.tpcInnerParam(), track.tpcSignal());
            }
          }
        }

        if (outFlagOptions.doTOFplots && track.hasTOF()) {
          if (outFlagOptions.enableEffPlots && enableDebug) {
            if (track.sign() > 0)
              debugHistos.fill(HIST("tracks/eff/hPtPTOF"), track.pt());
            else
              debugHistos.fill(HIST("tracks/eff/hPtantiPTOF"), track.pt());
          }

          if (outFlagOptions.enableBetaCut && (track.beta() > cfgBetaCut) && outFlagOptions.enablePIDplot)
            histos.fill(HIST("tracks/h2TOFbetaVsP_BetaCut"), track.p() / (1.f * track.sign()), track.beta());
          if (outFlagOptions.enablePIDplot) {
            switch (useHasTRDConfig) {
              case 0:
                histos.fill(HIST("tracks/h2TOFbetaVsP"), track.p() / (1.f * track.sign()), track.beta());
                break;
              case 1:
                if (track.hasTRD()) {
                  histos.fill(HIST("tracks/h2TOFbetaVsP"), track.p() / (1.f * track.sign()), track.beta());
                }
                break;
              case 2:
                if (!track.hasTRD()) {
                  histos.fill(HIST("tracks/h2TOFbetaVsP"), track.p() / (1.f * track.sign()), track.beta());
                }
                break;
            }
          }

          if (outFlagOptions.enableEffPlots)
            histos.fill(HIST("tracks/eff/h2TPCmomentumVsTOFExpMomentum"), track.tofExpMom(), track.tpcInnerParam());

          if (enablePr && prRapCut) {
            if (std::abs(track.tpcNSigmaPr()) < nsigmaTPCvar.nsigmaTPCPr && track.sign() > 0) {
              histos.fill(HIST("tracks/proton/h2ProtonTOFbetaVsP"), track.p(), track.beta());
              if (outFlagOptions.enableEffPlots) {
                histos.fill(HIST("tracks/eff/proton/h2pVsTOFExpMomentumPr"), track.tofExpMom(), track.p());
                histos.fill(HIST("tracks/eff/proton/h2TPCmomentumVsTOFExpMomentumPr"), track.tofExpMom(), track.tpcInnerParam());
              }
            }
            if (std::abs(track.tpcNSigmaPr()) < nsigmaTPCvar.nsigmaTPCPr && track.sign() < 0) {
              histos.fill(HIST("tracks/proton/h2antiProtonTOFbetaVsP"), track.p(), track.beta());
              if (outFlagOptions.enableEffPlots) {
                histos.fill(HIST("tracks/eff/proton/h2pVsTOFExpMomentumantiPr"), track.tofExpMom(), track.p());
                histos.fill(HIST("tracks/eff/proton/h2TPCmomentumVsTOFExpMomentumantiPr"), track.tofExpMom(), track.tpcInnerParam());
              }
            }
          }
          if (enableTr && trRapCut) {
            if (isTritonTPCpid && track.sign() > 0) {
              histos.fill(HIST("tracks/triton/h2TritonTOFbetaVsP"), track.p(), track.beta());
              if (outFlagOptions.enableEffPlots) {
                histos.fill(HIST("tracks/eff/triton/h2pVsTOFExpMomentumTr"), track.tofExpMom(), track.p());
                histos.fill(HIST("tracks/eff/triton/h2TPCmomentumVsTOFExpMomentumTr"), track.tofExpMom(), track.tpcInnerParam());
              }
            }
            if (isTritonTPCpid && track.sign() < 0) {
              histos.fill(HIST("tracks/triton/h2antiTritonTOFbetaVsP"), track.p(), track.beta());
              if (outFlagOptions.enableEffPlots) {
                histos.fill(HIST("tracks/eff/triton/h2pVsTOFExpMomentumantiTr"), track.tofExpMom(), track.p());
                histos.fill(HIST("tracks/eff/triton/h2TPCmomentumVsTOFExpMomentumantiTr"), track.tofExpMom(), track.tpcInnerParam());
              }
            }
          }
          if (filterOptions.enableEvTimeSplitting) {
            if (track.isEvTimeTOF() && track.isEvTimeT0AC()) {
              evtimeHistos.fill(HIST("tracks/evtime/ft0tof/h2TOFbetaVsP"), track.p() / (1.f * track.sign()), track.beta());
              evtimeHistos.fill(HIST("tracks/evtime/ft0tof/h2TPCsignVsTPCmomentum"), track.tpcInnerParam() / (1.f * track.sign()), track.tpcSignal());
            } else if (track.isEvTimeT0AC()) {
              evtimeHistos.fill(HIST("tracks/evtime/ft0/h2TOFbetaVsP"), track.p() / (1.f * track.sign()), track.beta());
              evtimeHistos.fill(HIST("tracks/evtime/ft0/h2TPCsignVsTPCmomentum"), track.tpcInnerParam() / (1.f * track.sign()), track.tpcSignal());
            } else if (track.isEvTimeTOF()) {
              evtimeHistos.fill(HIST("tracks/evtime/tof/h2TOFbetaVsP"), track.p() / (1.f * track.sign()), track.beta());
              evtimeHistos.fill(HIST("tracks/evtime/tof/h2TPCsignVsTPCmomentum"), track.tpcInnerParam() / (1.f * track.sign()), track.tpcSignal());
            } else {
              evtimeHistos.fill(HIST("tracks/evtime/fill/h2TOFbetaVsP"), track.p() / (1.f * track.sign()), track.beta());
              evtimeHistos.fill(HIST("tracks/evtime/fill/h2TPCsignVsTPCmomentum"), track.tpcInnerParam() / (1.f * track.sign()), track.tpcSignal());
            }
          }
        }
      }

      if (isDeWTPCpid) {
        if (outFlagOptions.enableEffPlots) {
          histos.fill(HIST("tracks/eff/deuteron/hPtDe"), DPt);
          histos.fill(HIST("tracks/eff/deuteron/h2pVsTPCmomentumDe"), track.tpcInnerParam(), track.p());
        }
        histos.fill(HIST("tracks/deuteron/h1DeuteronSpectra"), DPt);
        histos.fill(HIST("tracks/deuteron/h2DeuteronYvsPt"), track.rapidity(o2::track::PID::getMass2Z(o2::track::PID::Deuteron)), DPt);
        histos.fill(HIST("tracks/deuteron/h2DeuteronEtavsPt"), track.eta(), DPt);
        histos.fill(HIST("tracks/deuteron/h2DeuteronVspNSigmaITSDe_wTPCpid"), track.p(), nITSDe);

        if (outFlagOptions.enablePIDplot)
          histos.fill(HIST("tracks/deuteron/h2TPCsignVsTPCmomentumDeuteron"), track.tpcInnerParam(), track.tpcSignal());
      }
      if (isAntiDeWTPCpid) {
        if (outFlagOptions.enableEffPlots) {
          histos.fill(HIST("tracks/eff/deuteron/hPtantiDe"), antiDPt);
          histos.fill(HIST("tracks/eff/deuteron/h2pVsTPCmomentumantiDe"), track.tpcInnerParam(), track.p());
        }
        histos.fill(HIST("tracks/deuteron/h1antiDeuteronSpectra"), antiDPt);
        histos.fill(HIST("tracks/deuteron/h2antiDeuteronYvsPt"), track.rapidity(o2::track::PID::getMass2Z(o2::track::PID::Deuteron)), antiDPt);
        histos.fill(HIST("tracks/deuteron/h2antiDeuteronEtavsPt"), track.eta(), antiDPt);
        histos.fill(HIST("tracks/deuteron/h2antiDeuteronVspNSigmaITSDe_wTPCpid"), track.p(), nITSDe);

        if (outFlagOptions.enablePIDplot)
          histos.fill(HIST("tracks/deuteron/h2TPCsignVsTPCmomentumantiDeuteron"), track.tpcInnerParam(), track.tpcSignal());
      }
      if (isHeWTPCpid) {
        if (outFlagOptions.enableEffPlots) {
          histos.fill(HIST("tracks/eff/helium/hPtHe"), 2 * hePt);
          histos.fill(HIST("tracks/eff/helium/h2pVsTPCmomentumHe"), heTPCmomentum, heP);
        }
        histos.fill(HIST("tracks/helium/h1HeliumSpectra_Z2"), 2 * hePt);
        histos.fill(HIST("tracks/helium/h2HeliumYvsPt_Z2"), track.rapidity(o2::track::PID::getMass2Z(o2::track::PID::Helium3)), 2 * hePt);
        histos.fill(HIST("tracks/helium/h2HeliumEtavsPt_Z2"), track.eta(), 2 * hePt);
        if (outFlagOptions.enablePIDplot)
          histos.fill(HIST("tracks/helium/h2TPCsignVsTPCmomentumHelium"), heTPCmomentum, track.tpcSignal());
      }
      if (isAntiHeWTPCpid) {
        if (outFlagOptions.enableEffPlots) {
          histos.fill(HIST("tracks/eff/helium/hPtantiHe"), 2 * antihePt);
          histos.fill(HIST("tracks/eff/helium/h2pVsTPCmomentumantiHe"), antiheTPCmomentum, antiheP);
        }
        histos.fill(HIST("tracks/helium/h1antiHeliumSpectra_Z2"), 2 * antihePt);
        histos.fill(HIST("tracks/helium/h2antiHeliumYvsPt_Z2"), track.rapidity(o2::track::PID::getMass2Z(o2::track::PID::Helium3)), 2 * antihePt);
        histos.fill(HIST("tracks/helium/h2antiHeliumEtavsPt_Z2"), track.eta(), 2 * antihePt);
        if (outFlagOptions.enablePIDplot)
          histos.fill(HIST("tracks/helium/h2TPCsignVsTPCmomentumantiHelium"), antiheTPCmomentum, track.tpcSignal());
      }

      if (outFlagOptions.doTOFplots && track.hasTOF()) {
        if (isDeWTPCpid) {
          histos.fill(HIST("tracks/deuteron/h2DeuteronTOFbetaVsP"), track.p(), track.beta());
          if (outFlagOptions.enableEffPlots) {
            histos.fill(HIST("tracks/eff/deuteron/h2pVsTOFExpMomentumDe"), track.tofExpMom(), track.p());
            histos.fill(HIST("tracks/eff/deuteron/h2TPCmomentumVsTOFExpMomentumDe"), track.tofExpMom(), track.tpcInnerParam());
          }
        }
        if (isAntiDeWTPCpid) {
          histos.fill(HIST("tracks/deuteron/h2antiDeuteronTOFbetaVsP"), track.p(), track.beta());
          if (outFlagOptions.enableEffPlots) {
            histos.fill(HIST("tracks/eff/deuteron/h2pVsTOFExpMomentumantiDe"), track.tofExpMom(), track.p());
            histos.fill(HIST("tracks/eff/deuteron/h2TPCmomentumVsTOFExpMomentumantiDe"), track.tofExpMom(), track.tpcInnerParam());
          }
        }

        if (isHeWTPCpid) {
          histos.fill(HIST("tracks/helium/h2HeliumTOFbetaVsP"), heP, track.beta());
          if (outFlagOptions.enableEffPlots) {
            histos.fill(HIST("tracks/eff/helium/h2pVsTOFExpMomentumHe"), track.tofExpMom(), heP);
            histos.fill(HIST("tracks/eff/helium/h2TPCmomentumVsTOFExpMomentumHe"), track.tofExpMom(), heTPCmomentum);
          }
        }

        if (isAntiHeWTPCpid) {
          histos.fill(HIST("tracks/helium/h2antiHeliumTOFbetaVsP"), antiheP, track.beta());
          if (outFlagOptions.enableEffPlots) {
            histos.fill(HIST("tracks/eff/helium/h2pVsTOFExpMomentumantiHe"), track.tofExpMom(), antiheP);
            histos.fill(HIST("tracks/eff/helium/h2TPCmomentumVsTOFExpMomentumantiHe"), track.tofExpMom(), antiheTPCmomentum);
          }
        }

        if ((track.beta() * track.beta()) < 1.) {
          gamma = 1.f / std::sqrt(1.f - (track.beta() * track.beta()));

          switch (massTOFConfig) {
            case 0:
              massTOF = track.tpcInnerParam() * std::sqrt(1.f / (track.beta() * track.beta()) - 1.f);
              massTOFhe = heTPCmomentum * std::sqrt(1.f / (track.beta() * track.beta()) - 1.f);
              massTOFantihe = antiheTPCmomentum * std::sqrt(1.f / (track.beta() * track.beta()) - 1.f);
              break;
            case 1:
              massTOF = track.tofExpMom() * std::sqrt(1.f / (track.beta() * track.beta()) - 1.f);
              break;
            case 2:
              massTOF = track.p() * std::sqrt(1.f / (track.beta() * track.beta()) - 1.f);
              massTOFhe = heP * std::sqrt(1.f / (track.beta() * track.beta()) - 1.f);
              massTOFantihe = antiheP * std::sqrt(1.f / (track.beta() * track.beta()) - 1.f);
              break;
          }
          if (passDCAxyzCut && outFlagOptions.doTOFplots && outFlagOptions.enablePIDplot)
            histos.fill(HIST("tracks/h2TPCsignVsBetaGamma"), (track.beta() * gamma) / (1.f * track.sign()), track.tpcSignal());
        } else {
          massTOF = -99.f;
          massTOFhe = -99.f;
          massTOFantihe = -99.f;
        }

        if (passDCAxyzCut) {
          if (outFlagOptions.enablePIDplot)
            histos.fill(HIST("tracks/h2TOFmassVsPt"), massTOF, track.pt());
          if (filterOptions.enableEvTimeSplitting) {
            if (track.isEvTimeTOF() && track.isEvTimeT0AC()) {
              evtimeHistos.fill(HIST("tracks/evtime/ft0tof/h2TOFmassVsPt"), massTOF, track.pt());
            } else if (track.isEvTimeT0AC()) {
              evtimeHistos.fill(HIST("tracks/evtime/ft0/h2TOFmassVsPt"), massTOF, track.pt());
            } else if (track.isEvTimeTOF()) {
              evtimeHistos.fill(HIST("tracks/evtime/tof/h2TOFmassVsPt"), massTOF, track.pt());
            } else {
              evtimeHistos.fill(HIST("tracks/evtime/fill/h2TOFmassVsPt"), massTOF, track.pt());
            }
          }

          if (enablePr) {
            if ((std::abs(track.tpcNSigmaPr()) < nsigmaTPCvar.nsigmaTPCPr) && prRapCut) {
              if (track.sign() > 0) {
                if (outFlagOptions.enableEffPlots)
                  histos.fill(HIST("tracks/eff/proton/hPtPrTOF"), track.pt());
                histos.fill(HIST("tracks/proton/h2TOFmassProtonVsPt"), massTOF, track.pt());
                histos.fill(HIST("tracks/proton/h2TOFmass2ProtonVsPt"), massTOF * massTOF - MassProtonVal * MassProtonVal, track.pt());
                if (outFlagOptions.enableBetaCut && (track.beta() > cfgBetaCut)) {
                  histos.fill(HIST("tracks/proton/h2TOFmassProtonVsPt_BetaCut"), massTOF, track.pt());
                  histos.fill(HIST("tracks/proton/h2TOFmass2ProtonVsPt_BetaCut"), massTOF * massTOF - MassProtonVal * MassProtonVal, track.pt());
                }
                if (outFlagOptions.enableExpSignalTOF)
                  histos.fill(HIST("tracks/proton/h2ProtonTOFExpSignalDiffVsPtCut"), track.pt(), track.tofExpSignalDiffPr());
                if (filterOptions.enableEvTimeSplitting) {
                  if (track.isEvTimeTOF() && track.isEvTimeT0AC()) {
                    evtimeHistos.fill(HIST("tracks/evtime/ft0tof/proton/h2TOFmass2ProtonVsPt"), massTOF * massTOF - MassProtonVal * MassProtonVal, track.pt());
                  } else if (track.isEvTimeT0AC()) {
                    evtimeHistos.fill(HIST("tracks/evtime/ft0/proton/h2TOFmass2ProtonVsPt"), massTOF * massTOF - MassProtonVal * MassProtonVal, track.pt());
                  } else if (track.isEvTimeTOF()) {
                    evtimeHistos.fill(HIST("tracks/evtime/tof/proton/h2TOFmass2ProtonVsPt"), massTOF * massTOF - MassProtonVal * MassProtonVal, track.pt());
                  } else {
                    evtimeHistos.fill(HIST("tracks/evtime/fill/proton/h2TOFmass2ProtonVsPt"), massTOF * massTOF - MassProtonVal * MassProtonVal, track.pt());
                  }
                }
              } else {
                if (outFlagOptions.enableEffPlots)
                  histos.fill(HIST("tracks/eff/proton/hPtantiPrTOF"), track.pt());
                histos.fill(HIST("tracks/proton/h2TOFmassantiProtonVsPt"), massTOF, track.pt());
                histos.fill(HIST("tracks/proton/h2TOFmass2antiProtonVsPt"), massTOF * massTOF - MassProtonVal * MassProtonVal, track.pt());
                if (outFlagOptions.enableBetaCut && (track.beta() > cfgBetaCut)) {
                  histos.fill(HIST("tracks/proton/h2TOFmassantiProtonVsPt_BetaCut"), massTOF, track.pt());
                  histos.fill(HIST("tracks/proton/h2TOFmass2antiProtonVsPt_BetaCut"), massTOF * massTOF - MassProtonVal * MassProtonVal, track.pt());
                }
                if (outFlagOptions.enableExpSignalTOF)
                  histos.fill(HIST("tracks/proton/h2antiProtonTOFExpSignalDiffVsPtCut"), track.pt(), track.tofExpSignalDiffPr());
                if (filterOptions.enableEvTimeSplitting) {
                  if (track.isEvTimeTOF() && track.isEvTimeT0AC()) {
                    evtimeHistos.fill(HIST("tracks/evtime/ft0tof/proton/h2TOFmass2antiProtonVsPt"), massTOF * massTOF - MassProtonVal * MassProtonVal, track.pt());
                  } else if (track.isEvTimeT0AC()) {
                    evtimeHistos.fill(HIST("tracks/evtime/ft0/proton/h2TOFmass2antiProtonVsPt"), massTOF * massTOF - MassProtonVal * MassProtonVal, track.pt());
                  } else if (track.isEvTimeTOF()) {
                    evtimeHistos.fill(HIST("tracks/evtime/tof/proton/h2TOFmass2antiProtonVsPt"), massTOF * massTOF - MassProtonVal * MassProtonVal, track.pt());
                  } else {
                    evtimeHistos.fill(HIST("tracks/evtime/fill/proton/h2TOFmass2antiProtonVsPt"), massTOF * massTOF - MassProtonVal * MassProtonVal, track.pt());
                  }
                }
              }
            }
          }

          if (enableTr && isTritonTPCpid && trRapCut) {
            const float m2diff = massTOF * massTOF - MassTritonVal * MassTritonVal;
            if (track.sign() > 0) {
              if (outFlagOptions.enableEffPlots)
                histos.fill(HIST("tracks/eff/triton/hPtTrTOF"), track.pt());
              histos.fill(HIST("tracks/triton/h2TOFmassTritonVsPt"), massTOF, track.pt());
              histos.fill(HIST("tracks/triton/h2TOFmass2TritonVsPt"), m2diff, track.pt());
              if (outFlagOptions.enableBetaCut && (track.beta() > cfgBetaCut)) {
                histos.fill(HIST("tracks/triton/h2TOFmassTritonVsPt_BetaCut"), massTOF, track.pt());
                histos.fill(HIST("tracks/triton/h2TOFmass2TritonVsPt_BetaCut"), m2diff, track.pt());
              }
            } else {
              if (outFlagOptions.enableEffPlots)
                histos.fill(HIST("tracks/eff/triton/hPtantiTrTOF"), track.pt());
              histos.fill(HIST("tracks/triton/h2TOFmassantiTritonVsPt"), massTOF, track.pt());
              histos.fill(HIST("tracks/triton/h2TOFmass2antiTritonVsPt"), m2diff, track.pt());
              if (outFlagOptions.enableBetaCut && (track.beta() > cfgBetaCut)) {
                histos.fill(HIST("tracks/triton/h2TOFmassantiTritonVsPt_BetaCut"), massTOF, track.pt());
                histos.fill(HIST("tracks/triton/h2TOFmass2antiTritonVsPt_BetaCut"), m2diff, track.pt());
              }
            }
          }
        }

        if (isDeWTPCpid) {
          if (outFlagOptions.enableEffPlots)
            histos.fill(HIST("tracks/eff/deuteron/hPtDeTOF"), DPt);
          histos.fill(HIST("tracks/deuteron/h2TOFmassDeuteronVsPt"), massTOF, DPt);
          if (enableCentrality)
            histos.fill(HIST("tracks/deuteron/h3TOFmass2DeuteronVsPtVsMult"), massTOF * massTOF - MassDeuteronVal * MassDeuteronVal, DPt, event.centFT0M());
          else
            histos.fill(HIST("tracks/deuteron/h2TOFmass2DeuteronVsPt"), massTOF * massTOF - MassDeuteronVal * MassDeuteronVal, DPt);
          if (outFlagOptions.enableBetaCut && (track.beta() > cfgBetaCut)) {
            histos.fill(HIST("tracks/deuteron/h2TOFmassDeuteronVsPt_BetaCut"), massTOF, DPt);
            histos.fill(HIST("tracks/deuteron/h2TOFmass2DeuteronVsPt_BetaCut"), massTOF * massTOF - MassDeuteronVal * MassDeuteronVal, DPt);
          }
          if (outFlagOptions.enableExpSignalTOF)
            histos.fill(HIST("tracks/deuteron/h2DeuteronTOFExpSignalDiffVsPtCut"), DPt, track.tofExpSignalDiffDe());
          if (filterOptions.enableEvTimeSplitting) {
            if (track.isEvTimeTOF() && track.isEvTimeT0AC()) {
              evtimeHistos.fill(HIST("tracks/evtime/ft0tof/deuteron/h2TOFmass2DeuteronVsPt"), massTOF * massTOF - MassDeuteronVal * MassDeuteronVal, DPt);
            } else if (track.isEvTimeT0AC()) {
              evtimeHistos.fill(HIST("tracks/evtime/ft0/deuteron/h2TOFmass2DeuteronVsPt"), massTOF * massTOF - MassDeuteronVal * MassDeuteronVal, DPt);
            } else if (track.isEvTimeTOF()) {
              evtimeHistos.fill(HIST("tracks/evtime/tof/deuteron/h2TOFmass2DeuteronVsPt"), massTOF * massTOF - MassDeuteronVal * MassDeuteronVal, DPt);
            } else {
              evtimeHistos.fill(HIST("tracks/evtime/fill/deuteron/h2TOFmass2DeuteronVsPt"), massTOF * massTOF - MassDeuteronVal * MassDeuteronVal, DPt);
            }
          }
        }
        if (isAntiDeWTPCpid) {
          if (outFlagOptions.enableEffPlots)
            histos.fill(HIST("tracks/eff/deuteron/hPtantiDeTOF"), antiDPt);
          histos.fill(HIST("tracks/deuteron/h2TOFmassantiDeuteronVsPt"), massTOF, antiDPt);
          if (enableCentrality)
            histos.fill(HIST("tracks/deuteron/h3TOFmass2antiDeuteronVsPtVsMult"), massTOF * massTOF - MassDeuteronVal * MassDeuteronVal, antiDPt, event.centFT0M());
          else
            histos.fill(HIST("tracks/deuteron/h2TOFmass2antiDeuteronVsPt"), massTOF * massTOF - MassDeuteronVal * MassDeuteronVal, antiDPt);
          if (outFlagOptions.enableBetaCut && (track.beta() > cfgBetaCut)) {
            histos.fill(HIST("tracks/deuteron/h2TOFmassantiDeuteronVsPt_BetaCut"), massTOF, antiDPt);
            histos.fill(HIST("tracks/deuteron/h2TOFmass2antiDeuteronVsPt_BetaCut"), massTOF * massTOF - MassDeuteronVal * MassDeuteronVal, antiDPt);
          }
          if (outFlagOptions.enableExpSignalTOF)
            histos.fill(HIST("tracks/deuteron/h2antiDeuteronTOFExpSignalDiffVsPtCut"), antiDPt, track.tofExpSignalDiffDe());
          if (filterOptions.enableEvTimeSplitting) {
            if (track.isEvTimeTOF() && track.isEvTimeT0AC()) {
              evtimeHistos.fill(HIST("tracks/evtime/ft0tof/deuteron/h2TOFmass2antiDeuteronVsPt"), massTOF * massTOF - MassDeuteronVal * MassDeuteronVal, antiDPt);
            } else if (track.isEvTimeT0AC()) {
              evtimeHistos.fill(HIST("tracks/evtime/ft0/deuteron/h2TOFmass2antiDeuteronVsPt"), massTOF * massTOF - MassDeuteronVal * MassDeuteronVal, antiDPt);
            } else if (track.isEvTimeTOF()) {
              evtimeHistos.fill(HIST("tracks/evtime/tof/deuteron/h2TOFmass2antiDeuteronVsPt"), massTOF * massTOF - MassDeuteronVal * MassDeuteronVal, antiDPt);
            } else {
              evtimeHistos.fill(HIST("tracks/evtime/fill/deuteron/h2TOFmass2antiDeuteronVsPt"), massTOF * massTOF - MassDeuteronVal * MassDeuteronVal, antiDPt);
            }
          }
        }

        if (isHeWTPCpid) {
          if (outFlagOptions.enableEffPlots)
            histos.fill(HIST("tracks/eff/helium/hPtHeTOF"), 2 * hePt);
          histos.fill(HIST("tracks/helium/h2TOFmassHeliumVsPt"), 2.f * massTOFhe, hePt);
          histos.fill(HIST("tracks/helium/h2TOFmassDeltaHeliumVsPt"), 2.f * massTOFhe - MassHeliumVal, hePt);
          histos.fill(HIST("tracks/helium/h2TOFmass2HeliumVsPt"), 2.f * massTOFhe * 2.f * massTOFhe - MassHeliumVal * MassHeliumVal, hePt);
          if (enableCentrality)
            histos.fill(HIST("tracks/helium/h3TOFmass2HeliumVsPtVsMult"), 2.f * massTOFantihe * 2.f * massTOFantihe - MassHeliumVal * MassHeliumVal, hePt, event.centFT0M());
          if (outFlagOptions.enableBetaCut && (track.beta() > cfgBetaCut)) {
            histos.fill(HIST("tracks/helium/h2TOFmassHeliumVsPt_BetaCut"), 2.f * massTOFhe, hePt);
            histos.fill(HIST("tracks/helium/h2TOFmass2HeliumVsPt_BetaCut"), 2.f * massTOFhe * 2.f * massTOFhe - MassHeliumVal * MassHeliumVal, hePt);
          }
          if (outFlagOptions.enableExpSignalTOF)
            histos.fill(HIST("tracks/helium/h2HeliumTOFExpSignalDiffVsPtCut"), hePt, track.tofExpSignalDiffHe());
        }
        if (isAntiHeWTPCpid) {
          if (outFlagOptions.enableEffPlots)
            histos.fill(HIST("tracks/eff/helium/hPtantiHeTOF"), 2 * antihePt);
          histos.fill(HIST("tracks/helium/h2TOFmassantiHeliumVsPt"), 2.f * massTOFantihe, antihePt);
          histos.fill(HIST("tracks/helium/h2TOFmassDeltaantiHeliumVsPt"), 2.f * massTOFantihe - MassHeliumVal, antihePt);
          histos.fill(HIST("tracks/helium/h2TOFmass2antiHeliumVsPt"), 2.f * massTOFantihe * 2.f * massTOFantihe - MassHeliumVal * MassHeliumVal, antihePt);
          if (enableCentrality)
            histos.fill(HIST("tracks/helium/h3TOFmass2antiHeliumVsPtVsMult"), 2.f * massTOFantihe * 2.f * massTOFantihe - MassHeliumVal * MassHeliumVal, antihePt, event.centFT0M());
          if (outFlagOptions.enableBetaCut && (track.beta() > cfgBetaCut)) {
            histos.fill(HIST("tracks/helium/h2TOFmassantiHeliumVsPt_BetaCut"), 2.f * massTOFantihe, antihePt);
            histos.fill(HIST("tracks/helium/h2TOFmass2antiHeliumVsPt_BetaCut"), 2.f * massTOFantihe * 2.f * massTOFantihe - MassHeliumVal * MassHeliumVal, antihePt);
          }
          if (outFlagOptions.enableExpSignalTOF)
            histos.fill(HIST("tracks/helium/h2antiHeliumTOFExpSignalDiffVsPtCut"), antihePt, track.tofExpSignalDiffHe());
        }
      }

      if constexpr (IsMC) {
        // auto const& mcParticles = particles;
        bool isPhysPrim = false;
        bool isProdByGen = false;
        bool isWeakDecay = false;
        bool isItsPassed = false;
        bool isTpcPassed = false;
        bool isGoodHit = false;
        bool hasFakeHit = false;

        // PID
        int pdgCode = 0;
        if constexpr (IsFilteredData) {
          isPhysPrim = track.isPhysicalPrimary();
          isProdByGen = track.producedByGenerator();
          isWeakDecay = (track.getProcess() == TMCProcess::kPDecay);
          pdgCode = track.pdgCode();
          isItsPassed = track.itsPassed();
          isTpcPassed = track.tpcPassed();
          isGoodHit = !track.fakeHitsFlag();

        } else {
          if (!track.has_mcParticle()) {
            continue;
          }
          isPhysPrim = track.mcParticle().isPhysicalPrimary();
          isProdByGen = track.mcParticle().producedByGenerator();
          isWeakDecay = (track.mcParticle().getProcess() == TMCProcess::kPDecay);
          pdgCode = track.mcParticle().pdgCode();
          isItsPassed = track.passedITSNCls() &&
                        track.passedITSChi2NDF() &&
                        track.passedITSRefit() &&
                        track.passedITSHits() &&
                        track.hasITS();
          isTpcPassed = track.passedTPCNCls() &&
                        track.passedTPCCrossedRows() &&
                        track.passedTPCCrossedRowsOverNCls() &&
                        track.passedTPCChi2NDF() &&
                        track.passedTPCRefit() &&
                        track.hasTPC();

          for (int i = 0; i < kFakeLoop; i++) { // From ITS to TPC
            if (track.mcMask() & 1 << i) {
              hasFakeHit = true;
              break;
            }
          }
          isGoodHit = !hasFakeHit;
        }

        if (passDCAxyzCut) {
          if ((track.sign() > 0) && enableTrackingEff) {
            if (isItsPassed) {
              debugHistos.fill(HIST("tracks/trackingEff/h1_its"), track.pt());
              if (isGoodHit) {
                debugHistos.fill(HIST("tracks/trackingEff/h1GoodHit_its"), track.pt());
              }
            }

            if (isTpcPassed) {
              debugHistos.fill(HIST("tracks/trackingEff/h1_tpc"), track.pt());
              if (isGoodHit) {
                debugHistos.fill(HIST("tracks/trackingEff/h1GoodHit_tpc"), track.pt());
              }
            }

            if (isItsPassed && isTpcPassed) {
              debugHistos.fill(HIST("tracks/trackingEff/h1_its_tpc"), track.pt());
              if (isGoodHit) {
                debugHistos.fill(HIST("tracks/trackingEff/h1GoodHit_its_tpc"), track.pt());
              }
            }

            if (isItsPassed && isTpcPassed && track.hasTOF()) {
              debugHistos.fill(HIST("tracks/trackingEff/h1_its_tpc_tof"), track.pt());
              if (isGoodHit) {
                debugHistos.fill(HIST("tracks/trackingEff/h1GoodHit_its_tpc_tof"), track.pt());
              }
            }

            if (isItsPassed && isTpcPassed && track.hasTRD()) {
              debugHistos.fill(HIST("tracks/trackingEff/h1_its_tpc_trd"), track.pt());
              if (isGoodHit) {
                debugHistos.fill(HIST("tracks/trackingEff/h1GoodHit_its_tpc_trd"), track.pt());
              }
            }

            if (isItsPassed && isTpcPassed && track.hasTOF() && track.hasTRD()) {
              debugHistos.fill(HIST("tracks/trackingEff/h1_its_tpc_trd_tof"), track.pt());
              if (isGoodHit) {
                debugHistos.fill(HIST("tracks/trackingEff/h1GoodHit_its_tpc_trd_tof"), track.pt());
              }
            }
          } else if ((track.sign() < 0) && enableTrackingEff) {
            if (isItsPassed) {
              debugHistos.fill(HIST("tracks/trackingEff/h1anti_its"), track.pt());
              if (isGoodHit) {
                debugHistos.fill(HIST("tracks/trackingEff/h1antiGoodHit_its"), track.pt());
              }
            }

            if (isTpcPassed) {
              debugHistos.fill(HIST("tracks/trackingEff/h1anti_tpc"), track.pt());
              if (isGoodHit) {
                debugHistos.fill(HIST("tracks/trackingEff/h1antiGoodHit_tpc"), track.pt());
              }
            }

            if (isItsPassed && isTpcPassed) {
              debugHistos.fill(HIST("tracks/trackingEff/h1anti_its_tpc"), track.pt());
              if (isGoodHit) {
                debugHistos.fill(HIST("tracks/trackingEff/h1antiGoodHit_its_tpc"), track.pt());
              }
            }

            if (isItsPassed && isTpcPassed && track.hasTOF()) {
              debugHistos.fill(HIST("tracks/trackingEff/h1anti_its_tpc_tof"), track.pt());
              if (isGoodHit) {
                debugHistos.fill(HIST("tracks/trackingEff/h1antiGoodHit_its_tpc_tof"), track.pt());
              }
            }

            if (isItsPassed && isTpcPassed && track.hasTRD()) {
              debugHistos.fill(HIST("tracks/trackingEff/h1anti_its_tpc_trd"), track.pt());
              if (isGoodHit) {
                debugHistos.fill(HIST("tracks/trackingEff/h1antiGoodHit_its_tpc_trd"), track.pt());
              }
            }

            if (isItsPassed && isTpcPassed && track.hasTOF() && track.hasTRD()) {
              debugHistos.fill(HIST("tracks/trackingEff/h1anti_its_tpc_trd_tof"), track.pt());
              if (isGoodHit) {
                debugHistos.fill(HIST("tracks/trackingEff/h1antiGoodHit_its_tpc_trd_tof"), track.pt());
              }
            }
          }
        }
        switch (pdgCode) {
          case PDGProton:
            if (enablePr && prRapCut && passDCAxyzCut) {
              histos.fill(HIST("tracks/proton/h1ProtonSpectraTrue"), track.pt());
              if (outFlagOptions.makeDCAAfterCutPlots) {
                histos.fill(HIST("tracks/proton/dca/after/hDCAxyVsPtProtonTrue"), track.pt(), track.dcaXY());
                histos.fill(HIST("tracks/proton/dca/after/hDCAzVsPtProtonTrue"), track.pt(), track.dcaZ());
              }
              if (std::abs(track.tpcNSigmaPr()) < nsigmaTPCvar.nsigmaTPCPr) {
                histos.fill(HIST("tracks/proton/h1ProtonSpectraTrueWPID"), track.pt());
              }
              if (isPhysPrim) {
                histos.fill(HIST("tracks/proton/h1ProtonSpectraTruePrim"), track.pt());
                if (outFlagOptions.makeDCAAfterCutPlots) {
                  histos.fill(HIST("tracks/proton/dca/after/hDCAxyVsPtProtonTruePrim"), track.pt(), track.dcaXY());
                  histos.fill(HIST("tracks/proton/dca/after/hDCAzVsPtProtonTruePrim"), track.pt(), track.dcaZ());
                }

                if (enableTrackingEff) {
                  if (isItsPassed) {
                    debugHistos.fill(HIST("tracks/proton/trackingEff/h1ProtonSpectra_its"), track.pt());
                    if (isGoodHit) {
                      debugHistos.fill(HIST("tracks/proton/trackingEff/h1ProtonSpectraGoodHit_its"), track.pt());
                    }
                  }

                  if (isTpcPassed) {
                    debugHistos.fill(HIST("tracks/proton/trackingEff/h1ProtonSpectra_tpc"), track.pt());
                    if (isGoodHit) {
                      debugHistos.fill(HIST("tracks/proton/trackingEff/h1ProtonSpectraGoodHit_tpc"), track.pt());
                    }
                  }

                  if (isItsPassed && isTpcPassed) {
                    debugHistos.fill(HIST("tracks/proton/trackingEff/h1ProtonSpectra_its_tpc"), track.pt());
                    if (isGoodHit) {
                      debugHistos.fill(HIST("tracks/proton/trackingEff/h1ProtonSpectraGoodHit_its_tpc"), track.pt());
                    }
                  }

                  if (isItsPassed && isTpcPassed && track.hasTOF()) {
                    debugHistos.fill(HIST("tracks/proton/trackingEff/h1ProtonSpectra_its_tpc_tof"), track.pt());
                    if (isGoodHit) {
                      debugHistos.fill(HIST("tracks/proton/trackingEff/h1ProtonSpectraGoodHit_its_tpc_tof"), track.pt());
                    }
                  }

                  if (isItsPassed && isTpcPassed && track.hasTRD()) {
                    debugHistos.fill(HIST("tracks/proton/trackingEff/h1ProtonSpectra_its_tpc_trd"), track.pt());
                    if (isGoodHit) {
                      debugHistos.fill(HIST("tracks/proton/trackingEff/h1ProtonSpectraGoodHit_its_tpc_trd"), track.pt());
                    }
                  }

                  if (isItsPassed && isTpcPassed && track.hasTOF() && track.hasTRD()) {
                    debugHistos.fill(HIST("tracks/proton/trackingEff/h1ProtonSpectra_its_tpc_trd_tof"), track.pt());
                    if (isGoodHit) {
                      debugHistos.fill(HIST("tracks/proton/trackingEff/h1ProtonSpectraGoodHit_its_tpc_trd_tof"), track.pt());
                    }
                  }
                }

                if (std::abs(track.tpcNSigmaPr()) < nsigmaTPCvar.nsigmaTPCPr) {
                  if (enableTrackingEff) {
                    if (isItsPassed) {
                      debugHistos.fill(HIST("tracks/proton/trackingEffPID/h1ProtonSpectraPIDTruePrim_its"), track.pt());
                      if (isGoodHit) {
                        debugHistos.fill(HIST("tracks/proton/trackingEffPID/h1ProtonSpectraPIDTruePrimGoodHit_its"), track.pt());
                      }
                    }

                    if (isTpcPassed) {
                      debugHistos.fill(HIST("tracks/proton/trackingEffPID/h1ProtonSpectraPIDTruePrim_tpc"), track.pt());
                      if (isGoodHit) {
                        debugHistos.fill(HIST("tracks/proton/trackingEffPID/h1ProtonSpectraPIDTruePrimGoodHit_tpc"), track.pt());
                      }
                    }

                    if (isItsPassed && isTpcPassed) {
                      debugHistos.fill(HIST("tracks/proton/trackingEffPID/h1ProtonSpectraPIDTruePrim_its_tpc"), track.pt());
                      if (isGoodHit) {
                        debugHistos.fill(HIST("tracks/proton/trackingEffPID/h1ProtonSpectraPIDTruePrimGoodHit_its_tpc"), track.pt());
                      }
                    }

                    if (isItsPassed && isTpcPassed && track.hasTOF()) {
                      debugHistos.fill(HIST("tracks/proton/trackingEffPID/h1ProtonSpectraPIDTruePrim_its_tpc_tof"), track.pt());
                      if (isGoodHit) {
                        debugHistos.fill(HIST("tracks/proton/trackingEffPID/h1ProtonSpectraPIDTruePrimGoodHit_its_tpc_tof"), track.pt());
                      }
                    }

                    if (isItsPassed && isTpcPassed && track.hasTRD()) {
                      debugHistos.fill(HIST("tracks/proton/trackingEffPID/h1ProtonSpectraPIDTruePrim_its_tpc_trd"), track.pt());
                      if (isGoodHit) {
                        debugHistos.fill(HIST("tracks/proton/trackingEffPID/h1ProtonSpectraPIDTruePrimGoodHit_its_tpc_trd"), track.pt());
                      }
                    }

                    if (isItsPassed && isTpcPassed && track.hasTOF() && track.hasTRD()) {
                      debugHistos.fill(HIST("tracks/proton/trackingEffPID/h1ProtonSpectraPIDTruePrim_its_tpc_trd_tof"), track.pt());
                      if (isGoodHit) {
                        debugHistos.fill(HIST("tracks/proton/trackingEffPID/h1ProtonSpectraPIDTruePrimGoodHit_its_tpc_trd_tof"), track.pt());
                      }
                    }
                  }
                }
              }
              if (!isPhysPrim && !isProdByGen) {
                histos.fill(HIST("tracks/proton/h1ProtonSpectraTrueTransport"), track.pt());
                if (outFlagOptions.makeDCAAfterCutPlots) {
                  histos.fill(HIST("tracks/proton/dca/after/hDCAxyVsPtProtonTrueTransport"), track.pt(), track.dcaXY());
                  histos.fill(HIST("tracks/proton/dca/after/hDCAzVsPtProtonTrueTransport"), track.pt(), track.dcaZ());
                }
                if (isWeakDecay) {
                  histos.fill(HIST("tracks/proton/h1ProtonSpectraTrueSec"), track.pt());
                  if (outFlagOptions.makeDCAAfterCutPlots) {
                    histos.fill(HIST("tracks/proton/dca/after/hDCAxyVsPtProtonTrueSec"), track.pt(), track.dcaXY());
                    histos.fill(HIST("tracks/proton/dca/after/hDCAzVsPtProtonTrueSec"), track.pt(), track.dcaZ());
                  }
                }
              }
            }
            break;
          case -PDGProton:
            if (enablePr && prRapCut && passDCAxyzCut) {
              histos.fill(HIST("tracks/proton/h1antiProtonSpectraTrue"), track.pt());
              if (outFlagOptions.makeDCAAfterCutPlots) {
                histos.fill(HIST("tracks/proton/dca/after/hDCAxyVsPtantiProtonTrue"), track.pt(), track.dcaXY());
                histos.fill(HIST("tracks/proton/dca/after/hDCAzVsPtantiProtonTrue"), track.pt(), track.dcaZ());
              }
              if (std::abs(track.tpcNSigmaPr()) < nsigmaTPCvar.nsigmaTPCPr) {
                histos.fill(HIST("tracks/proton/h1antiProtonSpectraTrueWPID"), track.pt());
              }
              if (isPhysPrim) {
                histos.fill(HIST("tracks/proton/h1antiProtonSpectraTruePrim"), track.pt());
                if (outFlagOptions.makeDCAAfterCutPlots) {
                  histos.fill(HIST("tracks/proton/dca/after/hDCAxyVsPtantiProtonTruePrim"), track.pt(), track.dcaXY());
                  histos.fill(HIST("tracks/proton/dca/after/hDCAzVsPtantiProtonTruePrim"), track.pt(), track.dcaZ());
                }

                if (enableTrackingEff) {
                  if (isItsPassed) {
                    debugHistos.fill(HIST("tracks/proton/trackingEff/h1antiProtonSpectra_its"), track.pt());
                    if (isGoodHit) {
                      debugHistos.fill(HIST("tracks/proton/trackingEff/h1antiProtonSpectraGoodHit_its"), track.pt());
                    }
                  }

                  if (isTpcPassed) {
                    debugHistos.fill(HIST("tracks/proton/trackingEff/h1antiProtonSpectra_tpc"), track.pt());
                    if (isGoodHit) {
                      debugHistos.fill(HIST("tracks/proton/trackingEff/h1antiProtonSpectraGoodHit_tpc"), track.pt());
                    }
                  }

                  if (isItsPassed && isTpcPassed) {
                    debugHistos.fill(HIST("tracks/proton/trackingEff/h1antiProtonSpectra_its_tpc"), track.pt());
                    if (isGoodHit) {
                      debugHistos.fill(HIST("tracks/proton/trackingEff/h1antiProtonSpectraGoodHit_its_tpc"), track.pt());
                    }
                  }

                  if (isItsPassed && isTpcPassed && track.hasTOF()) {
                    debugHistos.fill(HIST("tracks/proton/trackingEff/h1antiProtonSpectra_its_tpc_tof"), track.pt());
                    if (isGoodHit) {
                      debugHistos.fill(HIST("tracks/proton/trackingEff/h1antiProtonSpectraGoodHit_its_tpc_tof"), track.pt());
                    }
                  }

                  if (isItsPassed && isTpcPassed && track.hasTRD()) {
                    debugHistos.fill(HIST("tracks/proton/trackingEff/h1antiProtonSpectra_its_tpc_trd"), track.pt());
                    if (isGoodHit) {
                      debugHistos.fill(HIST("tracks/proton/trackingEff/h1antiProtonSpectraGoodHit_its_tpc_trd"), track.pt());
                    }
                  }

                  if (isItsPassed && isTpcPassed && track.hasTOF() && track.hasTRD()) {
                    debugHistos.fill(HIST("tracks/proton/trackingEff/h1antiProtonSpectra_its_tpc_trd_tof"), track.pt());
                    if (isGoodHit) {
                      debugHistos.fill(HIST("tracks/proton/trackingEff/h1antiProtonSpectraGoodHit_its_tpc_trd_tof"), track.pt());
                    }
                  }
                }

                if (std::abs(track.tpcNSigmaPr()) < nsigmaTPCvar.nsigmaTPCPr) {
                  if (enableTrackingEff) {
                    if (isItsPassed) {
                      debugHistos.fill(HIST("tracks/proton/trackingEffPID/h1antiProtonSpectraPIDTruePrim_its"), track.pt());
                      if (isGoodHit) {
                        debugHistos.fill(HIST("tracks/proton/trackingEffPID/h1antiProtonSpectraPIDTruePrimGoodHit_its"), track.pt());
                      }
                    }

                    if (isTpcPassed) {
                      debugHistos.fill(HIST("tracks/proton/trackingEffPID/h1antiProtonSpectraPIDTruePrim_tpc"), track.pt());
                      if (isGoodHit) {
                        debugHistos.fill(HIST("tracks/proton/trackingEffPID/h1antiProtonSpectraPIDTruePrimGoodHit_tpc"), track.pt());
                      }
                    }

                    if (isItsPassed && isTpcPassed) {
                      debugHistos.fill(HIST("tracks/proton/trackingEffPID/h1antiProtonSpectraPIDTruePrim_its_tpc"), track.pt());
                      if (isGoodHit) {
                        debugHistos.fill(HIST("tracks/proton/trackingEffPID/h1antiProtonSpectraPIDTruePrimGoodHit_its_tpc"), track.pt());
                      }
                    }

                    if (isItsPassed && isTpcPassed && track.hasTOF()) {
                      debugHistos.fill(HIST("tracks/proton/trackingEffPID/h1antiProtonSpectraPIDTruePrim_its_tpc_tof"), track.pt());
                      if (isGoodHit) {
                        debugHistos.fill(HIST("tracks/proton/trackingEffPID/h1antiProtonSpectraPIDTruePrimGoodHit_its_tpc_tof"), track.pt());
                      }
                    }

                    if (isItsPassed && isTpcPassed && track.hasTRD()) {
                      debugHistos.fill(HIST("tracks/proton/trackingEffPID/h1antiProtonSpectraPIDTruePrim_its_tpc_trd"), track.pt());
                      if (isGoodHit) {
                        debugHistos.fill(HIST("tracks/proton/trackingEffPID/h1antiProtonSpectraPIDTruePrimGoodHit_its_tpc_trd"), track.pt());
                      }
                    }

                    if (isItsPassed && isTpcPassed && track.hasTOF() && track.hasTRD()) {
                      debugHistos.fill(HIST("tracks/proton/trackingEffPID/h1antiProtonSpectraPIDTruePrim_its_tpc_trd_tof"), track.pt());
                      if (isGoodHit) {
                        debugHistos.fill(HIST("tracks/proton/trackingEffPID/h1antiProtonSpectraPIDTruePrimGoodHit_its_tpc_trd_tof"), track.pt());
                      }
                    }
                  }
                }
              }

              if (!isPhysPrim && !isProdByGen) {
                histos.fill(HIST("tracks/proton/h1antiProtonSpectraTrueTransport"), track.pt());
                if (outFlagOptions.makeDCAAfterCutPlots) {
                  histos.fill(HIST("tracks/proton/dca/after/hDCAxyVsPtantiProtonTrueTransport"), track.pt(), track.dcaXY());
                  histos.fill(HIST("tracks/proton/dca/after/hDCAzVsPtantiProtonTrueTransport"), track.pt(), track.dcaZ());
                }
                if (isWeakDecay) {
                  histos.fill(HIST("tracks/proton/h1antiProtonSpectraTrueSec"), track.pt());
                  if (outFlagOptions.makeDCAAfterCutPlots) {
                    histos.fill(HIST("tracks/proton/dca/after/hDCAxyVsPtantiProtonTrueSec"), track.pt(), track.dcaXY());
                    histos.fill(HIST("tracks/proton/dca/after/hDCAzVsPtantiProtonTrueSec"), track.pt(), track.dcaZ());
                  }
                }
              }
            }
            break;
          case PDGDeuteron:
            if (isDeuteron && passDCAzCutDe && passDCAxyCutDe) {
              histos.fill(HIST("tracks/deuteron/h1DeuteronSpectraTrue"), DPt);
              if (outFlagOptions.makeDCAAfterCutPlots) {
                histos.fill(HIST("tracks/deuteron/dca/after/hDCAxyVsPtDeuteronTrue"), DPt, track.dcaXY());
                histos.fill(HIST("tracks/deuteron/dca/after/hDCAzVsPtDeuteronTrue"), DPt, track.dcaZ());
              }
              if (std::abs(track.tpcNSigmaDe()) < nsigmaTPCvar.nsigmaTPCDe) {
                histos.fill(HIST("tracks/deuteron/h1DeuteronSpectraTrueWPID"), DPt);
                if (track.hasTOF() && outFlagOptions.doTOFplots) {
                  histos.fill(HIST("tracks/deuteron/hPtDeuteronTOFTrue"), DPt);
                }
              }
              if (isPhysPrim) {
                histos.fill(HIST("tracks/deuteron/h2DeuteronVspTNSigmaTPCTruePrim"), DPt, track.tpcNSigmaDe());
                histos.fill(HIST("tracks/deuteron/h1DeuteronSpectraTruePrim"), DPt);
                if (outFlagOptions.makeDCAAfterCutPlots) {
                  histos.fill(HIST("tracks/deuteron/dca/after/hDCAxyVsPtDeuteronTruePrim"), DPt, track.dcaXY());
                  histos.fill(HIST("tracks/deuteron/dca/after/hDCAzVsPtDeuteronTruePrim"), DPt, track.dcaZ());
                }

                if (enableTrackingEff) {
                  if (isItsPassed) {
                    debugHistos.fill(HIST("tracks/deuteron/trackingEff/h1DeuteronSpectra_its"), DPt);
                    if (isGoodHit) {
                      debugHistos.fill(HIST("tracks/deuteron/trackingEff/h1DeuteronSpectraGoodHit_its"), DPt);
                    }
                  }

                  if (isTpcPassed) {
                    debugHistos.fill(HIST("tracks/deuteron/trackingEff/h1DeuteronSpectra_tpc"), DPt);
                    if (isGoodHit) {
                      debugHistos.fill(HIST("tracks/deuteron/trackingEff/h1DeuteronSpectraGoodHit_tpc"), DPt);
                    }
                  }

                  if (isItsPassed && isTpcPassed) {
                    debugHistos.fill(HIST("tracks/deuteron/trackingEff/h1DeuteronSpectra_its_tpc"), DPt);
                    if (isGoodHit) {
                      debugHistos.fill(HIST("tracks/deuteron/trackingEff/h1DeuteronSpectraGoodHit_its_tpc"), DPt);
                    }
                  }

                  if (isItsPassed && isTpcPassed && track.hasTOF()) {
                    debugHistos.fill(HIST("tracks/deuteron/trackingEff/h1DeuteronSpectra_its_tpc_tof"), DPt);
                    if (isGoodHit) {
                      debugHistos.fill(HIST("tracks/deuteron/trackingEff/h1DeuteronSpectraGoodHit_its_tpc_tof"), DPt);
                    }
                  }

                  if (isItsPassed && isTpcPassed && track.hasTRD()) {
                    debugHistos.fill(HIST("tracks/deuteron/trackingEff/h1DeuteronSpectra_its_tpc_trd"), DPt);
                    if (isGoodHit) {
                      debugHistos.fill(HIST("tracks/deuteron/trackingEff/h1DeuteronSpectraGoodHit_its_tpc_trd"), DPt);
                    }
                  }

                  if (isItsPassed && isTpcPassed && track.hasTOF() && track.hasTRD()) {
                    debugHistos.fill(HIST("tracks/deuteron/trackingEff/h1DeuteronSpectra_its_tpc_trd_tof"), DPt);
                    if (isGoodHit) {
                      debugHistos.fill(HIST("tracks/deuteron/trackingEff/h1DeuteronSpectraGoodHit_its_tpc_trd_tof"), DPt);
                    }
                  }
                }

                if (std::abs(track.tpcNSigmaDe()) < nsigmaTPCvar.nsigmaTPCDe) {
                  histos.fill(HIST("tracks/deuteron/h1DeuteronSpectraTrueWPIDPrim"), DPt);
                  if (outFlagOptions.makeDCAAfterCutPlots) {
                    histos.fill(HIST("tracks/deuteron/dca/after/hDCAxyVsPtDeuteronTrueWPIDPrim"), DPt, track.dcaXY());
                    histos.fill(HIST("tracks/deuteron/dca/after/hDCAzVsPtDeuteronTrueWPIDPrim"), DPt, track.dcaZ());
                  }
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    histos.fill(HIST("tracks/deuteron/hPtDeuteronTOFTrueWPIDPrim"), DPt);
                  }

                  if (enableTrackingEff) {
                    if (isItsPassed) {
                      debugHistos.fill(HIST("tracks/deuteron/trackingEffPID/h1DeuteronSpectraPIDTruePrim_its"), DPt);
                      if (isGoodHit) {
                        debugHistos.fill(HIST("tracks/deuteron/trackingEffPID/h1DeuteronSpectraPIDTruePrimGoodHit_its"), DPt);
                      }
                    }

                    if (isTpcPassed) {
                      debugHistos.fill(HIST("tracks/deuteron/trackingEffPID/h1DeuteronSpectraPIDTruePrim_tpc"), DPt);
                      if (isGoodHit) {
                        debugHistos.fill(HIST("tracks/deuteron/trackingEffPID/h1DeuteronSpectraPIDTruePrimGoodHit_tpc"), DPt);
                      }
                    }

                    if (isItsPassed && isTpcPassed) {
                      debugHistos.fill(HIST("tracks/deuteron/trackingEffPID/h1DeuteronSpectraPIDTruePrim_its_tpc"), DPt);
                      if (isGoodHit) {
                        debugHistos.fill(HIST("tracks/deuteron/trackingEffPID/h1DeuteronSpectraPIDTruePrimGoodHit_its_tpc"), DPt);
                      }
                    }

                    if (isItsPassed && isTpcPassed && track.hasTOF()) {
                      debugHistos.fill(HIST("tracks/deuteron/trackingEffPID/h1DeuteronSpectraPIDTruePrim_its_tpc_tof"), DPt);
                      if (isGoodHit) {
                        debugHistos.fill(HIST("tracks/deuteron/trackingEffPID/h1DeuteronSpectraPIDTruePrimGoodHit_its_tpc_tof"), DPt);
                      }
                    }

                    if (isItsPassed && isTpcPassed && track.hasTRD()) {
                      debugHistos.fill(HIST("tracks/deuteron/trackingEffPID/h1DeuteronSpectraPIDTruePrim_its_tpc_trd"), DPt);
                      if (isGoodHit) {
                        debugHistos.fill(HIST("tracks/deuteron/trackingEffPID/h1DeuteronSpectraPIDTruePrimGoodHit_its_tpc_trd"), DPt);
                      }
                    }

                    if (isItsPassed && isTpcPassed && track.hasTOF() && track.hasTRD()) {
                      debugHistos.fill(HIST("tracks/deuteron/trackingEffPID/h1DeuteronSpectraPIDTruePrim_its_tpc_trd_tof"), DPt);
                      if (isGoodHit) {
                        debugHistos.fill(HIST("tracks/deuteron/trackingEffPID/h1DeuteronSpectraPIDTruePrimGoodHit_its_tpc_trd_tof"), DPt);
                      }
                    }
                  }
                }
              }
              if (!isPhysPrim && !isProdByGen) {
                histos.fill(HIST("tracks/deuteron/h1DeuteronSpectraTrueTransport"), DPt);
                if (outFlagOptions.makeDCAAfterCutPlots) {
                  histos.fill(HIST("tracks/deuteron/dca/after/hDCAxyVsPtDeuteronTrueTransport"), DPt, track.dcaXY());
                  histos.fill(HIST("tracks/deuteron/dca/after/hDCAzVsPtDeuteronTrueTransport"), DPt, track.dcaZ());
                }
                if (isWeakDecay) {
                  histos.fill(HIST("tracks/deuteron/h1DeuteronSpectraTrueSec"), DPt);
                  if (outFlagOptions.makeDCAAfterCutPlots) {
                    histos.fill(HIST("tracks/deuteron/dca/after/hDCAxyVsPtDeuteronTrueSec"), DPt, track.dcaXY());
                    histos.fill(HIST("tracks/deuteron/dca/after/hDCAzVsPtDeuteronTrueSec"), DPt, track.dcaZ());
                  }
                }
              }
            }
            break;
          case -PDGDeuteron:
            if (isDeuteron && passDCAzCutAntiDe && passDCAxyCutAntiDe) {
              histos.fill(HIST("tracks/deuteron/h1antiDeuteronSpectraTrue"), antiDPt);
              if (outFlagOptions.makeDCAAfterCutPlots) {
                histos.fill(HIST("tracks/deuteron/dca/after/hDCAxyVsPtantiDeuteronTrue"), antiDPt, track.dcaXY());
                histos.fill(HIST("tracks/deuteron/dca/after/hDCAzVsPtantiDeuteronTrue"), antiDPt, track.dcaZ());
              }
              if (std::abs(track.tpcNSigmaDe()) < nsigmaTPCvar.nsigmaTPCDe) {
                histos.fill(HIST("tracks/deuteron/h1antiDeuteronSpectraTrueWPID"), antiDPt);
                if (track.hasTOF() && outFlagOptions.doTOFplots) {
                  histos.fill(HIST("tracks/deuteron/hPtantiDeuteronTOFTrue"), antiDPt);
                }
              }
              if (isPhysPrim) {
                histos.fill(HIST("tracks/deuteron/h2antiDeuteronVspTNSigmaTPCTruePrim"), antiDPt, track.tpcNSigmaDe());
                histos.fill(HIST("tracks/deuteron/h1antiDeuteronSpectraTruePrim"), antiDPt);
                if (outFlagOptions.makeDCAAfterCutPlots) {
                  histos.fill(HIST("tracks/deuteron/dca/after/hDCAxyVsPtantiDeuteronTruePrim"), antiDPt, track.dcaXY());
                  histos.fill(HIST("tracks/deuteron/dca/after/hDCAzVsPtantiDeuteronTruePrim"), antiDPt, track.dcaZ());
                }

                if (enableTrackingEff) {
                  if (isItsPassed) {
                    debugHistos.fill(HIST("tracks/deuteron/trackingEff/h1antiDeuteronSpectra_its"), antiDPt);
                    if (isGoodHit) {
                      debugHistos.fill(HIST("tracks/deuteron/trackingEff/h1antiDeuteronSpectraGoodHit_its"), antiDPt);
                    }
                  }

                  if (isTpcPassed) {
                    debugHistos.fill(HIST("tracks/deuteron/trackingEff/h1antiDeuteronSpectra_tpc"), antiDPt);
                    if (isGoodHit) {
                      debugHistos.fill(HIST("tracks/deuteron/trackingEff/h1antiDeuteronSpectraGoodHit_tpc"), antiDPt);
                    }
                  }

                  if (isItsPassed && isTpcPassed) {
                    debugHistos.fill(HIST("tracks/deuteron/trackingEff/h1antiDeuteronSpectra_its_tpc"), antiDPt);
                    if (isGoodHit) {
                      debugHistos.fill(HIST("tracks/deuteron/trackingEff/h1antiDeuteronSpectraGoodHit_its_tpc"), antiDPt);
                    }
                  }

                  if (isItsPassed && isTpcPassed && track.hasTOF()) {
                    debugHistos.fill(HIST("tracks/deuteron/trackingEff/h1antiDeuteronSpectra_its_tpc_tof"), antiDPt);
                    if (isGoodHit) {
                      debugHistos.fill(HIST("tracks/deuteron/trackingEff/h1antiDeuteronSpectraGoodHit_its_tpc_tof"), antiDPt);
                    }
                  }

                  if (isItsPassed && isTpcPassed && track.hasTRD()) {
                    debugHistos.fill(HIST("tracks/deuteron/trackingEff/h1antiDeuteronSpectra_its_tpc_trd"), antiDPt);
                    if (isGoodHit) {
                      debugHistos.fill(HIST("tracks/deuteron/trackingEff/h1antiDeuteronSpectraGoodHit_its_tpc_trd"), antiDPt);
                    }
                  }

                  if (isItsPassed && isTpcPassed && track.hasTOF() && track.hasTRD()) {
                    debugHistos.fill(HIST("tracks/deuteron/trackingEff/h1antiDeuteronSpectra_its_tpc_trd_tof"), antiDPt);
                    if (isGoodHit) {
                      debugHistos.fill(HIST("tracks/deuteron/trackingEff/h1antiDeuteronSpectraGoodHit_its_tpc_trd_tof"), antiDPt);
                    }
                  }
                }

                if (std::abs(track.tpcNSigmaDe()) < nsigmaTPCvar.nsigmaTPCDe) {
                  histos.fill(HIST("tracks/deuteron/h1antiDeuteronSpectraTrueWPIDPrim"), antiDPt);
                  if (outFlagOptions.makeDCAAfterCutPlots) {
                    histos.fill(HIST("tracks/deuteron/dca/after/hDCAxyVsPtantiDeuteronTrueWPIDPrim"), antiDPt, track.dcaXY());
                    histos.fill(HIST("tracks/deuteron/dca/after/hDCAzVsPtantiDeuteronTrueWPIDPrim"), antiDPt, track.dcaZ());
                  }
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    histos.fill(HIST("tracks/deuteron/hPtantiDeuteronTOFTrueWPIDPrim"), antiDPt);
                  }

                  if (enableTrackingEff) {
                    if (isItsPassed) {
                      debugHistos.fill(HIST("tracks/deuteron/trackingEffPID/h1antiDeuteronSpectraPIDTruePrim_its"), antiDPt);
                      if (isGoodHit) {
                        debugHistos.fill(HIST("tracks/deuteron/trackingEffPID/h1antiDeuteronSpectraPIDTruePrimGoodHit_its"), antiDPt);
                      }
                    }

                    if (isTpcPassed) {
                      debugHistos.fill(HIST("tracks/deuteron/trackingEffPID/h1antiDeuteronSpectraPIDTruePrim_tpc"), antiDPt);
                      if (isGoodHit) {
                        debugHistos.fill(HIST("tracks/deuteron/trackingEffPID/h1antiDeuteronSpectraPIDTruePrimGoodHit_tpc"), antiDPt);
                      }
                    }

                    if (isItsPassed && isTpcPassed) {
                      debugHistos.fill(HIST("tracks/deuteron/trackingEffPID/h1antiDeuteronSpectraPIDTruePrim_its_tpc"), antiDPt);
                      if (isGoodHit) {
                        debugHistos.fill(HIST("tracks/deuteron/trackingEffPID/h1antiDeuteronSpectraPIDTruePrimGoodHit_its_tpc"), antiDPt);
                      }
                    }

                    if (isItsPassed && isTpcPassed && track.hasTOF()) {
                      debugHistos.fill(HIST("tracks/deuteron/trackingEffPID/h1antiDeuteronSpectraPIDTruePrim_its_tpc_tof"), antiDPt);
                      if (isGoodHit) {
                        debugHistos.fill(HIST("tracks/deuteron/trackingEffPID/h1antiDeuteronSpectraPIDTruePrimGoodHit_its_tpc_tof"), antiDPt);
                      }
                    }

                    if (isItsPassed && isTpcPassed && track.hasTRD()) {
                      debugHistos.fill(HIST("tracks/deuteron/trackingEffPID/h1antiDeuteronSpectraPIDTruePrim_its_tpc_trd"), antiDPt);
                      if (isGoodHit) {
                        debugHistos.fill(HIST("tracks/deuteron/trackingEffPID/h1antiDeuteronSpectraPIDTruePrimGoodHit_its_tpc_trd"), antiDPt);
                      }
                    }

                    if (isItsPassed && isTpcPassed && track.hasTOF() && track.hasTRD()) {
                      debugHistos.fill(HIST("tracks/deuteron/trackingEffPID/h1antiDeuteronSpectraPIDTruePrim_its_tpc_trd_tof"), antiDPt);
                      if (isGoodHit) {
                        debugHistos.fill(HIST("tracks/deuteron/trackingEffPID/h1antiDeuteronSpectraPIDTruePrimGoodHit_its_tpc_trd_tof"), antiDPt);
                      }
                    }
                  }
                }
              }
              if (!isPhysPrim && !isProdByGen) {
                histos.fill(HIST("tracks/deuteron/h1antiDeuteronSpectraTrueTransport"), antiDPt);
                if (outFlagOptions.makeDCAAfterCutPlots) {
                  histos.fill(HIST("tracks/deuteron/dca/after/hDCAxyVsPtantiDeuteronTrueTransport"), antiDPt, track.dcaXY());
                  histos.fill(HIST("tracks/deuteron/dca/after/hDCAzVsPtantiDeuteronTrueTransport"), antiDPt, track.dcaZ());
                }
                if (isWeakDecay) {
                  histos.fill(HIST("tracks/deuteron/h1antiDeuteronSpectraTrueSec"), antiDPt);
                  if (outFlagOptions.makeDCAAfterCutPlots) {
                    histos.fill(HIST("tracks/deuteron/dca/after/hDCAxyVsPtantiDeuteronTrueSec"), antiDPt, track.dcaXY());
                    histos.fill(HIST("tracks/deuteron/dca/after/hDCAzVsPtantiDeuteronTrueSec"), antiDPt, track.dcaZ());
                  }
                }
              }
            }
            break;
          case PDGTriton:
            if (enableTr && trRapCut && passDCAxyzCut) {
              histos.fill(HIST("tracks/triton/h1TritonSpectraTrue"), track.pt());
              if (outFlagOptions.makeDCAAfterCutPlots) {
                histos.fill(HIST("tracks/triton/dca/after/hDCAxyVsPtTritonTrue"), track.pt(), track.dcaXY());
                histos.fill(HIST("tracks/triton/dca/after/hDCAzVsPtTritonTrue"), track.pt(), track.dcaZ());
              }
              if (isPhysPrim) {
                histos.fill(HIST("tracks/triton/h1TritonSpectraTruePrim"), track.pt());
                if (outFlagOptions.makeDCAAfterCutPlots) {
                  histos.fill(HIST("tracks/triton/dca/after/hDCAxyVsPtTritonTruePrim"), track.pt(), track.dcaXY());
                  histos.fill(HIST("tracks/triton/dca/after/hDCAzVsPtTritonTruePrim"), track.pt(), track.dcaZ());
                }
              }
              if (!isPhysPrim && !isProdByGen) {
                histos.fill(HIST("tracks/triton/h1TritonSpectraTrueTransport"), track.pt());
                if (outFlagOptions.makeDCAAfterCutPlots) {
                  histos.fill(HIST("tracks/triton/dca/after/hDCAxyVsPtTritonTrueTransport"), track.pt(), track.dcaXY());
                  histos.fill(HIST("tracks/triton/dca/after/hDCAzVsPtTritonTrueTransport"), track.pt(), track.dcaZ());
                }
                if (isWeakDecay) {
                  histos.fill(HIST("tracks/triton/h1TritonSpectraTrueSec"), track.pt());
                  if (outFlagOptions.makeDCAAfterCutPlots) {
                    histos.fill(HIST("tracks/triton/dca/after/hDCAxyVsPtTritonTrueSec"), track.pt(), track.dcaXY());
                    histos.fill(HIST("tracks/triton/dca/after/hDCAzVsPtTritonTrueSec"), track.pt(), track.dcaZ());
                  }
                }
              }
            }
            break;
          case -PDGTriton:
            if (enableTr && trRapCut && passDCAxyzCut) {
              histos.fill(HIST("tracks/triton/h1antiTritonSpectraTrue"), track.pt());
              if (outFlagOptions.makeDCAAfterCutPlots) {
                histos.fill(HIST("tracks/triton/dca/after/hDCAxyVsPtantiTritonTrue"), track.pt(), track.dcaXY());
                histos.fill(HIST("tracks/triton/dca/after/hDCAzVsPtantiTritonTrue"), track.pt(), track.dcaZ());
              }
              if (isPhysPrim) {
                histos.fill(HIST("tracks/triton/h1antiTritonSpectraTruePrim"), track.pt());
                if (outFlagOptions.makeDCAAfterCutPlots) {
                  histos.fill(HIST("tracks/triton/dca/after/hDCAxyVsPtantiTritonTruePrim"), track.pt(), track.dcaXY());
                  histos.fill(HIST("tracks/triton/dca/after/hDCAzVsPtantiTritonTruePrim"), track.pt(), track.dcaZ());
                }
              }
              if (!isPhysPrim && !isProdByGen) {
                histos.fill(HIST("tracks/triton/h1antiTritonSpectraTrueTransport"), track.pt());
                if (outFlagOptions.makeDCAAfterCutPlots) {
                  histos.fill(HIST("tracks/triton/dca/after/hDCAxyVsPtantiTritonTrueTransport"), track.pt(), track.dcaXY());
                  histos.fill(HIST("tracks/triton/dca/after/hDCAzVsPtantiTritonTrueTransport"), track.pt(), track.dcaZ());
                }
                if (isWeakDecay) {
                  histos.fill(HIST("tracks/triton/h1antiTritonSpectraTrueSec"), track.pt());
                  if (outFlagOptions.makeDCAAfterCutPlots) {
                    histos.fill(HIST("tracks/triton/dca/after/hDCAxyVsPtantiTritonTrueSec"), track.pt(), track.dcaXY());
                    histos.fill(HIST("tracks/triton/dca/after/hDCAzVsPtantiTritonTrueSec"), track.pt(), track.dcaZ());
                  }
                }
              }
            }
            break;
          case PDGHelium:
            if (isHelium && passDCAzCutHe && passDCAxyCutHe) {
              histos.fill(HIST("tracks/helium/h1HeliumSpectraTrue_Z2"), 2 * hePt);
              if (enableCentrality)
                histos.fill(HIST("tracks/helium/h2HeliumSpectraTrueVsMult_Z2"), 2 * hePt, event.centFT0M());

              if (outFlagOptions.makeDCAAfterCutPlots) {
                histos.fill(HIST("tracks/helium/dca/after/hDCAxyVsPtHeliumTrue"), hePt, track.dcaXY());
                histos.fill(HIST("tracks/helium/dca/after/hDCAzVsPtHeliumTrue"), hePt, track.dcaZ());

                if (track.hasTOF() && outFlagOptions.doTOFplots) {
                  histos.fill(HIST("tracks/helium/dca/after/TOF/hDCAxyVsPtHeliumTrue"), hePt, track.dcaXY());
                  histos.fill(HIST("tracks/helium/dca/after/TOF/hDCAzVsPtHeliumTrue"), hePt, track.dcaZ());
                }
              }
              if (std::abs(track.tpcNSigmaHe()) < nsigmaTPCvar.nsigmaTPCHe) {
                histos.fill(HIST("tracks/helium/h1HeliumSpectraTrueWPID_Z2"), 2 * hePt);
                if (enableCentrality)
                  histos.fill(HIST("tracks/helium/h2HeliumSpectraTrueWPIDVsMult_Z2"), 2 * hePt, event.centFT0M());
                if (outFlagOptions.enableEffPlots) {
                  histos.fill(HIST("tracks/eff/helium/hPtHeTrue_Z2"), 2 * hePt);
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    histos.fill(HIST("tracks/eff/helium/hPtHeTOFTrue_Z2"), 2 * hePt);
                  }
                }
              }
              if (isPhysPrim) {
                histos.fill(HIST("tracks/helium/h1HeliumSpectraTruePrim_Z2"), 2 * hePt);
                if (enableCentrality)
                  histos.fill(HIST("tracks/helium/h2HeliumSpectraTruePrimVsMult_Z2"), 2 * hePt, event.centFT0M());

                if (std::abs(track.tpcNSigmaHe()) < nsigmaTPCvar.nsigmaTPCHe) {
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    histos.fill(HIST("tracks/helium/TOF/h1HeliumSpectraTruePrim_Z2"), 2 * hePt);
                    if (enableCentrality)
                      histos.fill(HIST("tracks/helium/TOF/h2HeliumSpectraTruePrimVsMult_Z2"), 2 * hePt, event.centFT0M());
                  }
                }

                if (outFlagOptions.makeDCAAfterCutPlots) {
                  histos.fill(HIST("tracks/helium/dca/after/hDCAxyVsPtHeliumTruePrim"), hePt, track.dcaXY());
                  histos.fill(HIST("tracks/helium/dca/after/hDCAzVsPtHeliumTruePrim"), hePt, track.dcaZ());
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    histos.fill(HIST("tracks/helium/dca/after/TOF/hDCAxyVsPtHeliumTruePrim"), hePt, track.dcaXY());
                    histos.fill(HIST("tracks/helium/dca/after/TOF/hDCAzVsPtHeliumTruePrim"), hePt, track.dcaZ());
                  }
                }
              }

              if (!isPhysPrim && !isProdByGen) {
                histos.fill(HIST("tracks/helium/h1HeliumSpectraTrueTransport_Z2"), 2 * hePt);
                if (outFlagOptions.makeDCAAfterCutPlots) {
                  histos.fill(HIST("tracks/helium/dca/after/hDCAxyVsPtHeliumTrueTransport"), hePt, track.dcaXY());
                  histos.fill(HIST("tracks/helium/dca/after/hDCAzVsPtHeliumTrueTransport"), hePt, track.dcaZ());
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    histos.fill(HIST("tracks/helium/dca/after/TOF/hDCAxyVsPtHeliumTrueTransport"), hePt, track.dcaXY());
                    histos.fill(HIST("tracks/helium/dca/after/TOF/hDCAzVsPtHeliumTrueTransport"), hePt, track.dcaZ());
                  }
                }
                if (isWeakDecay) {
                  histos.fill(HIST("tracks/helium/h1HeliumSpectraTrueSec_Z2"), 2 * hePt);
                  if (enableCentrality)
                    histos.fill(HIST("tracks/helium/h2HeliumSpectraTrueSecVsMult_Z2"), 2 * hePt, event.centFT0M());
                  if (outFlagOptions.makeDCAAfterCutPlots) {
                    histos.fill(HIST("tracks/helium/dca/after/hDCAxyVsPtHeliumTrueSec"), hePt, track.dcaXY());
                    histos.fill(HIST("tracks/helium/dca/after/hDCAzVsPtHeliumTrueSec"), hePt, track.dcaZ());
                    if (track.hasTOF() && outFlagOptions.doTOFplots) {
                      histos.fill(HIST("tracks/helium/dca/after/TOF/hDCAxyVsPtHeliumTrueSec"), hePt, track.dcaXY());
                      histos.fill(HIST("tracks/helium/dca/after/TOF/hDCAzVsPtHeliumTrueSec"), hePt, track.dcaZ());
                    }
                  }
                }
              }
            }
            break;
          case -PDGHelium:
            if (isHelium && passDCAzCutAntiHe && passDCAxyCutAntiHe) {
              histos.fill(HIST("tracks/helium/h1antiHeliumSpectraTrue_Z2"), 2 * antihePt);
              if (enableCentrality)
                histos.fill(HIST("tracks/helium/h2antiHeliumSpectraTrueVsMult_Z2"), 2 * antihePt, event.centFT0M());

              if (outFlagOptions.makeDCAAfterCutPlots) {
                histos.fill(HIST("tracks/helium/dca/after/hDCAxyVsPtantiHeliumTrue"), antihePt, track.dcaXY());
                histos.fill(HIST("tracks/helium/dca/after/hDCAzVsPtantiHeliumTrue"), antihePt, track.dcaZ());
                if (track.hasTOF() && outFlagOptions.doTOFplots) {
                  histos.fill(HIST("tracks/helium/dca/after/TOF/hDCAxyVsPtantiHeliumTrue"), antihePt, track.dcaXY());
                  histos.fill(HIST("tracks/helium/dca/after/TOF/hDCAzVsPtantiHeliumTrue"), antihePt, track.dcaZ());
                }
              }
              if (std::abs(track.tpcNSigmaHe()) < nsigmaTPCvar.nsigmaTPCHe) {
                histos.fill(HIST("tracks/helium/h1antiHeliumSpectraTrueWPID_Z2"), 2 * antihePt);
                if (enableCentrality)
                  histos.fill(HIST("tracks/helium/h2antiHeliumSpectraTrueWPIDVsMult_Z2"), 2 * antihePt, event.centFT0M());
                if (outFlagOptions.enableEffPlots) {
                  histos.fill(HIST("tracks/eff/helium/hPtantiHeTrue_Z2"), 2 * antihePt);
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    histos.fill(HIST("tracks/eff/helium/hPtantiHeTOFTrue_Z2"), 2 * antihePt);
                  }
                }
              }
              if (isPhysPrim) {
                histos.fill(HIST("tracks/helium/h1antiHeliumSpectraTruePrim_Z2"), 2 * antihePt);
                if (enableCentrality)
                  histos.fill(HIST("tracks/helium/h2antiHeliumSpectraTruePrimVsMult_Z2"), 2 * antihePt, event.centFT0M());

                if (std::abs(track.tpcNSigmaHe()) < nsigmaTPCvar.nsigmaTPCHe) {
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    histos.fill(HIST("tracks/helium/TOF/h1antiHeliumSpectraTruePrim_Z2"), 2 * antihePt);
                    if (enableCentrality)
                      histos.fill(HIST("tracks/helium/TOF/h2antiHeliumSpectraTruePrimVsMult_Z2"), 2 * antihePt, event.centFT0M());
                  }
                }

                if (outFlagOptions.makeDCAAfterCutPlots) {
                  histos.fill(HIST("tracks/helium/dca/after/hDCAxyVsPtantiHeliumTruePrim"), antihePt, track.dcaXY());
                  histos.fill(HIST("tracks/helium/dca/after/hDCAzVsPtantiHeliumTruePrim"), antihePt, track.dcaZ());
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    histos.fill(HIST("tracks/helium/dca/after/TOF/hDCAxyVsPtantiHeliumTruePrim"), antihePt, track.dcaXY());
                    histos.fill(HIST("tracks/helium/dca/after/TOF/hDCAzVsPtantiHeliumTruePrim"), antihePt, track.dcaZ());
                  }
                }
              }
              if (!isPhysPrim && !isProdByGen) {
                histos.fill(HIST("tracks/helium/h1antiHeliumSpectraTrueTransport_Z2"), 2 * antihePt);
                if (outFlagOptions.makeDCAAfterCutPlots) {
                  histos.fill(HIST("tracks/helium/dca/after/hDCAxyVsPtantiHeliumTrueTransport"), antihePt, track.dcaXY());
                  histos.fill(HIST("tracks/helium/dca/after/hDCAzVsPtantiHeliumTrueTransport"), antihePt, track.dcaZ());
                  if (track.hasTOF() && outFlagOptions.doTOFplots) {
                    histos.fill(HIST("tracks/helium/dca/after/TOF/hDCAxyVsPtantiHeliumTrueTransport"), antihePt, track.dcaXY());
                    histos.fill(HIST("tracks/helium/dca/after/TOF/hDCAzVsPtantiHeliumTrueTransport"), antihePt, track.dcaZ());
                  }
                }

                if (isWeakDecay) {
                  histos.fill(HIST("tracks/helium/h1antiHeliumSpectraTrueSec_Z2"), 2 * antihePt);
                  if (enableCentrality)
                    histos.fill(HIST("tracks/helium/h2antiHeliumSpectraTrueSecVsMult_Z2"), 2 * antihePt, event.centFT0M());
                  if (outFlagOptions.makeDCAAfterCutPlots) {
                    histos.fill(HIST("tracks/helium/dca/after/hDCAxyVsPtantiHeliumTrueSec"), antihePt, track.dcaXY());
                    histos.fill(HIST("tracks/helium/dca/after/hDCAzVsPtantiHeliumTrueSec"), antihePt, track.dcaZ());
                    if (track.hasTOF() && outFlagOptions.doTOFplots) {
                      histos.fill(HIST("tracks/helium/dca/after/TOF/hDCAxyVsPtantiHeliumTrueSec"), antihePt, track.dcaXY());
                      histos.fill(HIST("tracks/helium/dca/after/TOF/hDCAzVsPtantiHeliumTrueSec"), antihePt, track.dcaZ());
                    }
                  }
                }
              }
            }
            break;
          case PDGAlpha:
            if (enableAl && alRapCut && passDCAxyzCut) {
              histos.fill(HIST("tracks/alpha/h1AlphaSpectraTrue"), track.pt());
              if (outFlagOptions.makeDCAAfterCutPlots) {
                histos.fill(HIST("tracks/alpha/dca/after/hDCAxyVsPtAlphaTrue"), track.pt(), track.dcaXY());
                histos.fill(HIST("tracks/alpha/dca/after/hDCAzVsPtAlphaTrue"), track.pt(), track.dcaZ());
              }
              if (isPhysPrim) {
                histos.fill(HIST("tracks/alpha/h1AlphaSpectraTruePrim"), track.pt());
                if (outFlagOptions.makeDCAAfterCutPlots) {
                  histos.fill(HIST("tracks/alpha/dca/after/hDCAxyVsPtAlphaTruePrim"), track.pt(), track.dcaXY());
                  histos.fill(HIST("tracks/alpha/dca/after/hDCAzVsPtAlphaTruePrim"), track.pt(), track.dcaZ());
                }
              }
              if (!isPhysPrim && !isProdByGen) {
                histos.fill(HIST("tracks/alpha/h1AlphaSpectraTrueTransport"), track.pt());
                if (outFlagOptions.makeDCAAfterCutPlots) {
                  histos.fill(HIST("tracks/alpha/dca/after/hDCAxyVsPtAlphaTrueTransport"), track.pt(), track.dcaXY());
                  histos.fill(HIST("tracks/alpha/dca/after/hDCAzVsPtAlphaTrueTransport"), track.pt(), track.dcaZ());
                }

                if (isWeakDecay) {
                  histos.fill(HIST("tracks/alpha/h1AlphaSpectraTrueSec"), track.pt());
                  if (outFlagOptions.makeDCAAfterCutPlots) {
                    histos.fill(HIST("tracks/alpha/dca/after/hDCAxyVsPtAlphaTrueSec"), track.pt(), track.dcaXY());
                    histos.fill(HIST("tracks/alpha/dca/after/hDCAzVsPtAlphaTrueSec"), track.pt(), track.dcaZ());
                  }
                }
              }
            }
            break;
          case -PDGAlpha:
            if (enableAl && alRapCut && passDCAxyzCut) {
              histos.fill(HIST("tracks/alpha/h1antiAlphaSpectraTrue"), track.pt());
              if (outFlagOptions.makeDCAAfterCutPlots) {
                histos.fill(HIST("tracks/alpha/dca/after/hDCAxyVsPtantiAlphaTrue"), track.pt(), track.dcaXY());
                histos.fill(HIST("tracks/alpha/dca/after/hDCAzVsPtantiAlphaTrue"), track.pt(), track.dcaZ());
              }
              if (isPhysPrim) {
                histos.fill(HIST("tracks/alpha/h1antiAlphaSpectraTruePrim"), track.pt());
                if (outFlagOptions.makeDCAAfterCutPlots) {
                  histos.fill(HIST("tracks/alpha/dca/after/hDCAxyVsPtantiAlphaTruePrim"), track.pt(), track.dcaXY());
                  histos.fill(HIST("tracks/alpha/dca/after/hDCAzVsPtantiAlphaTruePrim"), track.pt(), track.dcaZ());
                }
              }
              if (!isPhysPrim && !isProdByGen) {
                histos.fill(HIST("tracks/alpha/h1antiAlphaSpectraTrueTransport"), track.pt());
                if (outFlagOptions.makeDCAAfterCutPlots) {
                  histos.fill(HIST("tracks/alpha/dca/after/hDCAxyVsPtantiAlphaTrueTransport"), track.pt(), track.dcaXY());
                  histos.fill(HIST("tracks/alpha/dca/after/hDCAzVsPtantiAlphaTrueTransport"), track.pt(), track.dcaZ());
                }
                if (isWeakDecay) {
                  histos.fill(HIST("tracks/alpha/h1antiAlphaSpectraTrueSec"), track.pt());
                  if (outFlagOptions.makeDCAAfterCutPlots) {
                    histos.fill(HIST("tracks/alpha/dca/after/hDCAxyVsPtantiAlphaTrueSec"), track.pt(), track.dcaXY());
                    histos.fill(HIST("tracks/alpha/dca/after/hDCAzVsPtantiAlphaTrueSec"), track.pt(), track.dcaZ());
                  }
                }
              }
            }
            break;
          default:
            break;
        }
      } else {
        (void)particles;
      }
    }
  }

  using EventCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFV0As, aod::CentFT0Cs>;
  using EventCandidatesMC = soa::Join<EventCandidates, aod::McCollisionLabels>;

  using TrackCandidates0 = soa::Join<aod::Tracks, aod::TracksIU, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::TrackSelectionExtension,
                                     aod::pidTOFbeta, aod::TOFSignal, aod::pidEvTimeFlags,
                                     aod::pidTPCFullPi, aod::pidTOFFullPi,
                                     aod::pidTPCFullKa, aod::pidTOFFullKa,
                                     aod::pidTOFFullPr,
                                     aod::pidTOFFullDe,
                                     aod::pidTOFFullTr,
                                     aod::pidTOFFullHe,
                                     aod::pidTOFFullAl>;
  using TrackCandidates = soa::Join<TrackCandidates0,
                                    aod::pidTPCFullPr,
                                    aod::pidTPCFullDe,
                                    aod::pidTPCFullTr,
                                    aod::pidTPCFullHe,
                                    aod::pidTPCFullAl>;

  using TrackCandidatesLfPid = soa::Join<TrackCandidates0,
                                         aod::pidTPCLfFullPr,
                                         aod::pidTPCLfFullDe,
                                         aod::pidTPCLfFullTr,
                                         aod::pidTPCLfFullHe,
                                         aod::pidTPCLfFullAl>;

  //////////
  // Data //
  //////////

  // Process function that runs on the original AO2D
  void processData(EventCandidates::iterator const& event,
                   TrackCandidates const& tracks,
                   o2::aod::BCsWithTimestamps const&)
  {
    fillHistograms<false /*MC*/, false /*Filtered*/>(event, tracks, true /*dummy*/);
  }
  PROCESS_SWITCH(LFNucleiBATask, processData, "process data", true);

  // Process function that runs on the original AO2D
  void processDataLfPid(EventCandidates::iterator const& event,
                        TrackCandidatesLfPid const& tracks,
                        o2::aod::BCsWithTimestamps const&)
  {
    fillHistograms<false /*MC*/, false /*Filtered*/>(event, tracks, true /*dummy*/);
  }
  PROCESS_SWITCH(LFNucleiBATask, processDataLfPid, "process data with LF PID", false);

  // Process function that runs on the filtered data
  void processDataFiltered(o2::aod::LfNuclEvents::iterator const& event,
                           o2::aod::LfCandNucleusFull const& tracks,
                           o2::aod::BCsWithTimestamps const&)
  {
    // Runs on data filtered on the fly with LF Tree creator nuclei task
    // Takes as input full AO2Ds
    fillHistograms<false /*MC*/, true /*Filtered*/>(event, tracks, true /*dummy*/);
  }
  PROCESS_SWITCH(LFNucleiBATask, processDataFiltered, "process data on the filtered data", false);

  void processDataLight(o2::aod::LfNuclEvents::iterator const& event,
                        o2::aod::LfCandNucleusDummy const& tracks,
                        o2::aod::BCsWithTimestamps const&)
  {
    // Runs on derived tables produced with LF Tree creator nuclei task
    // Takes as input derived trees
    fillHistograms<false /*MC*/, true /*Filtered*/>(event, tracks, true /*dummy*/);
  }
  PROCESS_SWITCH(LFNucleiBATask, processDataLight, "process data on the derived trees", false);

  /////////////
  // MC Reco //
  /////////////

  Preslice<aod::McParticles> perMCCol = aod::mcparticle::mcCollisionId;
  SliceCache cache;

  // Process function that runs on the original AO2D (for the MC)
  void processMCReco(EventCandidatesMC::iterator const& event,
                     soa::Join<aod::McCollisions, aod::McCentFT0Ms> const&,
                     soa::Join<TrackCandidates, aod::McTrackLabels> const& tracks,
                     aod::McParticles const& mcParticles,
                     o2::aod::BCsWithTimestamps const&)
  {
    bool doRecoSep = true;

    const bool hasTVX = event.selection_bit(aod::evsel::kIsTriggerTVX);
    const bool hasNoTFB = event.selection_bit(aod::evsel::kNoTimeFrameBorder);
    const bool hasNoItsRofFB = event.selection_bit(aod::evsel::kNoITSROFrameBorder);

    if (evselOptions.useSel8) {
      doRecoSep = (hasTVX && hasNoTFB);
    } else {
      if (evselOptions.useTVXtrigger && !hasTVX)
        doRecoSep = false;
      if (evselOptions.removeTFBorder && !hasNoTFB)
        doRecoSep = false;
      if (evselOptions.removeITSROFBorder && !hasNoItsRofFB)
        doRecoSep = false;
    }

    if (doRecoSep && event.has_mcCollision()) {
      const int mcIdx = event.mcCollisionId();
      if (mcIdx >= 0) {
        effEvtSet.insert(mcIdx);
        effEvtSetReady = true;
      }
    }

    fillHistograms<true /*MC*/, false /*Filtered*/>(event, tracks, mcParticles);
  } // CLOSING PROCESS MC RECO
  PROCESS_SWITCH(LFNucleiBATask, processMCReco, "process mc reco", false);

  // Process function that runs on the original AO2D (for the MC) with the LfPIDcalibration
  void processMCRecoLfPid(EventCandidatesMC::iterator const& event,
                          soa::Join<aod::McCollisions, aod::McCentFT0Ms> const&,
                          soa::Join<TrackCandidatesLfPid, aod::McTrackLabels> const& tracks,
                          aod::McParticles const& mcParticles,
                          o2::aod::BCsWithTimestamps const&)
  {
    bool doRecoSep = true;

    const bool hasTVX = event.selection_bit(aod::evsel::kIsTriggerTVX);
    const bool hasNoTFB = event.selection_bit(aod::evsel::kNoTimeFrameBorder);
    const bool hasNoItsRofFB = event.selection_bit(aod::evsel::kNoITSROFrameBorder);

    if (evselOptions.useSel8) {
      doRecoSep = (hasTVX && hasNoTFB);
    } else {
      if (evselOptions.useTVXtrigger && !hasTVX)
        doRecoSep = false;
      if (evselOptions.removeTFBorder && !hasNoTFB)
        doRecoSep = false;
      if (evselOptions.removeITSROFBorder && !hasNoItsRofFB)
        doRecoSep = false;
    }

    if (doRecoSep && event.has_mcCollision()) {
      const int mcIdx = event.mcCollisionId();
      if (mcIdx >= 0) {
        effEvtSet.insert(mcIdx);
        effEvtSetReady = true;
      }
    }

    fillHistograms<true /*MC*/, false /*Filtered*/>(event, tracks, mcParticles);
  } // CLOSING PROCESS MC RECO
  PROCESS_SWITCH(LFNucleiBATask, processMCRecoLfPid, "process mc reco with LfPid", false);

  // Process function that runs on the original AO2D (for the MC) with the LfPIDcalibration
  void processMCRecoLfPidEv(EventCandidatesMC const& collisions,
                            soa::Join<TrackCandidatesLfPid, aod::McTrackLabels> const&,
                            aod::McParticles const& mcParticles,
                            aod::McCollisions const&)
  {
    for (const auto& collision : collisions) {
      if (!collision.has_mcCollision())
        continue;

      const auto& particlesInCollision = mcParticles.sliceByCached(aod::mcparticle::mcCollisionId, collision.mcCollision().globalIndex(), cache);

      const bool hasTVX = collision.selection_bit(aod::evsel::kIsTriggerTVX);
      const bool hasNoTFB = collision.selection_bit(aod::evsel::kNoTimeFrameBorder);
      const bool hasNoItsRofFB = collision.selection_bit(aod::evsel::kNoITSROFrameBorder);
      const bool hasSel8 = collision.sel8();

      for (const auto& mcParticle : particlesInCollision) {
        const float rapidity = mcParticle.y();
        if (rapidity > kinemOptions.cfgRapidityCutHigh || rapidity < kinemOptions.cfgRapidityCutLow)
          continue;

        const int pdg = mcParticle.pdgCode();
        const float pt = mcParticle.pt();
        bool isPhysPrim = mcParticle.isPhysicalPrimary();

        // No cut
        spectraGen.fill(HIST("LfEv/pT_nocut"), pt);
        if (pdg == PDGHelium) {
          spectraGen.fill(HIST("LfEv/helium/pT_nocut_He"), pt);
          if (isPhysPrim)
            spectraGen.fill(HIST("LfEv/helium/prim/pT_nocut_He"), pt);
        }
        if (pdg == -PDGHelium) {
          spectraGen.fill(HIST("LfEv/helium/pT_nocut_antiHe"), pt);
          if (isPhysPrim)
            spectraGen.fill(HIST("LfEv/helium/prim/pT_nocut_antiHe"), pt);
        }
        // Trigger TVX
        if (hasTVX) {
          spectraGen.fill(HIST("LfEv/pT_TVXtrigger"), pt);
          if (pdg == PDGHelium) {
            spectraGen.fill(HIST("LfEv/helium/pT_TVXtrigger_He"), pt);
            if (isPhysPrim)
              spectraGen.fill(HIST("LfEv/helium/prim/pT_TVXtrigger_He"), pt);
          }
          if (pdg == -PDGHelium) {
            spectraGen.fill(HIST("LfEv/helium/pT_TVXtrigger_antiHe"), pt);
            if (isPhysPrim)
              spectraGen.fill(HIST("LfEv/helium/prim/pT_TVXtrigger_antiHe"), pt);
          }
        }
        // No Time Frame Border
        if (hasNoTFB) {
          spectraGen.fill(HIST("LfEv/pT_TFrameBorder"), pt);
          if (pdg == PDGHelium) {
            spectraGen.fill(HIST("LfEv/helium/pT_TFrameBorder_He"), pt);
            if (isPhysPrim)
              spectraGen.fill(HIST("LfEv/helium/prim/pT_TFrameBorder_He"), pt);
          }
          if (pdg == -PDGHelium) {
            spectraGen.fill(HIST("LfEv/helium/pT_TFrameBorder_antiHe"), pt);
            if (isPhysPrim)
              spectraGen.fill(HIST("LfEv/helium/prim/pT_TFrameBorder_antiHe"), pt);
          }
        }
        // No ITS ROF Frame Border
        if (hasNoItsRofFB) {
          spectraGen.fill(HIST("LfEv/pT_ITSROFBorder"), pt);
          if (pdg == PDGHelium) {
            spectraGen.fill(HIST("LfEv/helium/pT_ITSROFBorder_He"), pt);
            if (isPhysPrim)
              spectraGen.fill(HIST("LfEv/helium/prim/pT_ITSROFBorder_He"), pt);
          }
          if (pdg == -PDGHelium) {
            spectraGen.fill(HIST("LfEv/helium/pT_ITSROFBorder_antiHe"), pt);
            if (isPhysPrim)
              spectraGen.fill(HIST("LfEv/helium/prim/pT_ITSROFBorder_antiHe"), pt);
          }
        }
        // Sel8 MC
        if (hasTVX && hasNoTFB) {
          spectraGen.fill(HIST("LfEv/pT_MCsel8"), pt);
          if (pdg == PDGHelium) {
            spectraGen.fill(HIST("LfEv/helium/pT_MCsel8_He"), pt);
            if (isPhysPrim)
              spectraGen.fill(HIST("LfEv/helium/prim/pT_MCsel8_He"), pt);
          }
          if (pdg == -PDGHelium) {
            spectraGen.fill(HIST("LfEv/helium/pT_MCsel8_antiHe"), pt);
            if (isPhysPrim)
              spectraGen.fill(HIST("LfEv/helium/prim/pT_MCsel8_antiHe"), pt);
          }
        }
        // Sel8 tag
        if (hasSel8) {
          spectraGen.fill(HIST("LfEv/pT_sel8"), pt);
          if (pdg == PDGHelium)
            spectraGen.fill(HIST("LfEv/helium/pT_sel8_He"), pt);
          if (isPhysPrim)
            spectraGen.fill(HIST("LfEv/helium/prim/pT_sel8_He"), pt);
          if (pdg == -PDGHelium)
            spectraGen.fill(HIST("LfEv/helium/pT_sel8_antiHe"), pt);
          if (isPhysPrim)
            spectraGen.fill(HIST("LfEv/helium/prim/pT_sel8_antiHe"), pt);
        }
      }
    }
  }
  // CLOSING PROCESS MC RECO
  PROCESS_SWITCH(LFNucleiBATask, processMCRecoLfPidEv, "process mc reco with LfPid w/ Event", false);

  // Process function that runs on the filtered AO2D (for the MC)
  void processMCRecoFiltered(o2::aod::LfNuclEvents::iterator const& event,
                             soa::Join<o2::aod::LfCandNucleusFull, o2::aod::LfCandNucleusMC> const& tracks,
                             o2::aod::BCsWithTimestamps const&)
  {
    fillHistograms<true /*MC*/, true /*Filtered*/>(event, tracks, true /*dummy*/);
  } // CLOSING PROCESS MC RECO ON FILTERED DATA
  PROCESS_SWITCH(LFNucleiBATask, processMCRecoFiltered, "process mc reco on the filtered data", false);

  void processMCRecoFilteredLight(o2::aod::LfNuclEvents::iterator const& event,
                                  soa::Join<o2::aod::LfCandNucleusDummy, o2::aod::LfCandNucleusMC> const& tracks,
                                  o2::aod::BCsWithTimestamps const&)
  {
    fillHistograms<true /*MC*/, true /*Filtered*/>(event, tracks, true /*dummy*/);
  } // CLOSING PROCESS MC RECO ON FILTERED DATA
  PROCESS_SWITCH(LFNucleiBATask, processMCRecoFilteredLight, "process mc reco on the derived trees", false);

  ////////////
  // MC Gen //
  ////////////

  // LOOP OVER GENERATED MC PARTICLES
  void processMCGen(soa::Join<aod::McCollisions, aod::McCentFT0Ms>::iterator const& mcCollision,
                    aod::McParticles const& mcParticles)
  {
    // Only events that are reconstructed
    const int mcIdx = mcCollision.globalIndex();
    if (enableEffEvtSet) {
      if (!effEvtSetReady)
        return;
      if (!effEvtSet.count(mcIdx))
        return;
    }

    if (mcCollision.centFT0M() < cfgMultCutLow || mcCollision.centFT0M() > cfgMultCutHigh)
      return;

    if (enableCentrality)
      spectraGen.fill(HIST("histGenVetxZ"), mcCollision.posZ(), mcCollision.centFT0M());
    else
      spectraGen.fill(HIST("histGenVetxZ"), mcCollision.posZ());

    // const auto& particlesInCollision = mcParticles.sliceByCached(aod::mcparticle::mcCollisionId, mcIdx, cache);
    // const auto particlesInCollision = mcParticles.sliceBy(perMCCol, mcIdx);

    for (const auto& mcParticleGen : mcParticles) {
      if (mcParticleGen.mcCollisionId() != mcIdx)
        continue;

      if (mcParticleGen.y() > kinemOptions.cfgRapidityCutHigh || mcParticleGen.y() < kinemOptions.cfgRapidityCutLow)
        continue;

      const int pdgCode = mcParticleGen.pdgCode();
      const float ptMC = mcParticleGen.pt();
      bool isPhysPrim = mcParticleGen.isPhysicalPrimary();
      bool isProdByGen = mcParticleGen.producedByGenerator();
      bool isWeakDecay = (mcParticleGen.getProcess() == TMCProcess::kPDecay);

      if (pdgCode == PDGPion) {
        spectraGen.fill(HIST("pion/histGenPtPion"), ptMC);
        if (isPhysPrim)
          spectraGen.fill(HIST("pion/histGenPtPionPrim"), ptMC);
        if (!isPhysPrim && isProdByGen) {
          //
        }
        if (!isPhysPrim && !isProdByGen) {
          spectraGen.fill(HIST("pion/histSecTransportPtPion"), ptMC);
          if (isWeakDecay) {
            spectraGen.fill(HIST("pion/histGenPtPionSec"), ptMC);
          }
        }
      }
      if (pdgCode == -PDGPion) {
        spectraGen.fill(HIST("pion/histGenPtantiPion"), ptMC);
        if (isPhysPrim)
          spectraGen.fill(HIST("pion/histGenPtantiPionPrim"), ptMC);
        if (!isPhysPrim && isProdByGen) {
          //
        }
        if (!isPhysPrim && !isProdByGen) {
          spectraGen.fill(HIST("pion/histSecTransportPtantiPion"), ptMC);
          if (isWeakDecay) {
            spectraGen.fill(HIST("pion/histGenPtantiPionSec"), ptMC);
          }
        }
      }
      if (pdgCode == PDGKaon) {
        spectraGen.fill(HIST("kaon/histGenPtKaon"), ptMC);
        if (isPhysPrim)
          spectraGen.fill(HIST("kaon/histGenPtKaonPrim"), ptMC);
        if (!isPhysPrim && isProdByGen) {
          //
        }
        if (!isPhysPrim && !isProdByGen) {
          spectraGen.fill(HIST("kaon/histSecTransportPtKaon"), ptMC);
          if (isWeakDecay) {
            spectraGen.fill(HIST("kaon/histGenPtKaonSec"), ptMC);
          }
        }
      }
      if (pdgCode == -PDGKaon) {
        spectraGen.fill(HIST("kaon/histGenPtantiKaon"), ptMC);
        if (isPhysPrim)
          spectraGen.fill(HIST("kaon/histGenPtantiKaonPrim"), ptMC);
        if (!isPhysPrim && isProdByGen) {
          //
        }
        if (!isPhysPrim && !isProdByGen) {
          spectraGen.fill(HIST("kaon/histSecTransportPtantiKaon"), ptMC);
          if (isWeakDecay) {
            spectraGen.fill(HIST("kaon/histGenPtantiKaonSec"), ptMC);
          }
        }
      }
      if (enablePr) {
        if (pdgCode == PDGProton) {
          spectraGen.fill(HIST("proton/histGenPtProton"), ptMC);
          if (isPhysPrim) {
            spectraGen.fill(HIST("proton/histGenPtProtonPrim"), ptMC);
            spectraGen.fill(HIST("proton/histGenPtProtonPrim_Y"), mcParticleGen.y(), ptMC);
          }
          if (!isPhysPrim && isProdByGen) {
            //
          }
          if (!isPhysPrim && !isProdByGen) {
            spectraGen.fill(HIST("proton/histSecTransportPtProton"), ptMC);
            if (isWeakDecay) {
              spectraGen.fill(HIST("proton/histGenPtProtonSec"), ptMC);
            }
          }
        }
        if (pdgCode == -PDGProton) {
          spectraGen.fill(HIST("proton/histGenPtantiProton"), ptMC);
          if (isPhysPrim) {
            spectraGen.fill(HIST("proton/histGenPtantiProtonPrim"), ptMC);
            spectraGen.fill(HIST("proton/histGenPtantiProtonPrim_Y"), mcParticleGen.y(), ptMC);
          }
          if (!isPhysPrim && isProdByGen) {
            //
          }
          if (!isPhysPrim && !isProdByGen) {
            spectraGen.fill(HIST("proton/histSecTransportPtantiProton"), ptMC);
            if (isWeakDecay) {
              spectraGen.fill(HIST("proton/histGenPtantiProtonSec"), ptMC);
            }
          }
        }
      }
      if (enableDe) {
        if (pdgCode == PDGDeuteron) {
          spectraGen.fill(HIST("deuteron/histGenPtD"), ptMC);
          if (isPhysPrim)
            spectraGen.fill(HIST("deuteron/histGenPtDPrim"), ptMC);
          if (!isPhysPrim && isProdByGen) {
            //
          }
          if (!isPhysPrim && !isProdByGen) {
            spectraGen.fill(HIST("deuteron/histSecTransportPtD"), ptMC);
            if (isWeakDecay) {
              spectraGen.fill(HIST("deuteron/histGenPtDSec"), ptMC);
            }
          }
        }
        if (pdgCode == -PDGDeuteron) {
          spectraGen.fill(HIST("deuteron/histGenPtantiD"), ptMC);
          if (isPhysPrim)
            spectraGen.fill(HIST("deuteron/histGenPtantiDPrim"), ptMC);
          if (!isPhysPrim && isProdByGen) {
            //
          }
          if (!isPhysPrim && !isProdByGen) {
            spectraGen.fill(HIST("deuteron/histSecTransportPtantiD"), ptMC);
            if (isWeakDecay) {
              spectraGen.fill(HIST("deuteron/histGenPtantiDSec"), ptMC);
            }
          }
        }
      }
      if (enableTr) {
        if (pdgCode == PDGTriton) {
          spectraGen.fill(HIST("triton/histGenPtT"), ptMC);
          if (isPhysPrim)
            spectraGen.fill(HIST("triton/histGenPtTPrim"), ptMC);
          if (!isPhysPrim && isProdByGen) {
            //
          }
          if (!isPhysPrim && !isProdByGen) {
            spectraGen.fill(HIST("triton/histSecTransportPtT"), ptMC);
            if (isWeakDecay) {
              spectraGen.fill(HIST("triton/histGenPtTSec"), ptMC);
            }
          }
        }
        if (pdgCode == -PDGTriton) {
          spectraGen.fill(HIST("triton/histGenPtantiT"), ptMC);
          if (isPhysPrim)
            spectraGen.fill(HIST("triton/histGenPtantiTPrim"), ptMC);
          if (!isPhysPrim && isProdByGen) {
            //
          }
          if (!isPhysPrim && !isProdByGen) {
            spectraGen.fill(HIST("triton/histSecTransportPtantiT"), ptMC);
            if (isWeakDecay) {
              spectraGen.fill(HIST("triton/histGenPtantiTSec"), ptMC);
            }
          }
        }
      }
      if (enableHe) {
        if (pdgCode == PDGHelium) {
          spectraGen.fill(HIST("helium/histGenPtHe"), ptMC);
          if (isPhysPrim) {
            spectraGen.fill(HIST("helium/histGenPtHePrim"), ptMC);
            if (enableCentrality)
              spectraGen.fill(HIST("helium/histGenPtHePrimVsMult"), ptMC, mcCollision.centFT0M());
          }
          if (!isPhysPrim && isProdByGen) {
            {
              //
            }
          }
          if (!isPhysPrim && !isProdByGen) {
            spectraGen.fill(HIST("helium/histSecTransportPtHe"), ptMC);
            if (isWeakDecay) {
              spectraGen.fill(HIST("helium/histGenPtHeSec"), ptMC);
            }
          }
        }
        if (pdgCode == -PDGHelium) {
          spectraGen.fill(HIST("helium/histGenPtantiHe"), ptMC);
          if (isPhysPrim) {
            spectraGen.fill(HIST("helium/histGenPtantiHePrim"), ptMC);
            if (enableCentrality)
              spectraGen.fill(HIST("helium/histGenPtantiHePrimVsMult"), ptMC, mcCollision.centFT0M());
          }
          if (!isPhysPrim && isProdByGen) {
            {
              //
            }
          }
          if (!isPhysPrim && !isProdByGen) {
            spectraGen.fill(HIST("helium/histSecTransportPtantiHe"), ptMC);
            if (isWeakDecay) {
              spectraGen.fill(HIST("helium/histGenPtantiHeSec"), ptMC);
            }
          }
        }
      }
      if (enableAl) {
        if (pdgCode == PDGAlpha) {
          spectraGen.fill(HIST("alpha/histGenPtAl"), ptMC);
          if (isPhysPrim)
            spectraGen.fill(HIST("alpha/histGenPtAlPrim"), ptMC);
          if (!isPhysPrim && isProdByGen) {
            //
          }
          if (!isPhysPrim && !isProdByGen) {
            spectraGen.fill(HIST("alpha/histSecTransportPtAl"), ptMC);
            if (isWeakDecay) {
              spectraGen.fill(HIST("alpha/histGenPtAlSec"), ptMC);
            }
          }
        }
        if (pdgCode == -PDGAlpha) {
          spectraGen.fill(HIST("alpha/histGenPtantiAl"), ptMC);
          if (isPhysPrim)
            spectraGen.fill(HIST("alpha/histGenPtantiAlPrim"), ptMC);
          if (!isPhysPrim && isProdByGen) {
            //
          }
          if (!isPhysPrim && !isProdByGen) {
            spectraGen.fill(HIST("alpha/histSecTransportPtantiAl"), ptMC);
            if (isWeakDecay) {
              spectraGen.fill(HIST("alpha/histGenPtantiAlSec"), ptMC);
            }
          }
        }
      }
    }
  } // Close processMCGen
  PROCESS_SWITCH(LFNucleiBATask, processMCGen, "process MC Generated", true);

  void processEvSgLossMC(soa::Join<aod::McCollisions, aod::McCentFT0Ms>::iterator const& mcCollision,
                         aod::McParticles const& mcParticles,
                         const soa::SmallGroups<EventCandidatesMC>& recoColls)
  {
    if (std::abs(mcCollision.posZ()) < cfgVzCutHigh)
      evLossHistos.fill(HIST("evLoss/hEvent"), 0.5);

    bool isSel8Event = false;
    bool isTvxEvent = false;
    bool isNoTFBEvent = false;
    bool isMCSel8Event = false;

    // Check if there is an event reconstructed for a generated event
    for (const auto& recoColl : recoColls) {
      if (std::abs(recoColl.posZ()) > cfgVzCutHigh)
        continue;
      if (recoColl.selection_bit(o2::aod::evsel::kIsTriggerTVX))
        isTvxEvent = true;
      if (recoColl.sel8())
        isSel8Event = true;
      if (recoColl.selection_bit(aod::evsel::kNoTimeFrameBorder))
        isNoTFBEvent = true;
      if (isTvxEvent && isSel8Event && isNoTFBEvent)
        break; // Optimize loop
    }

    if (isTvxEvent && isNoTFBEvent)
      isMCSel8Event = true;

    if (isTvxEvent)
      evLossHistos.fill(HIST("evLoss/hEvent"), 1.5);
    if (isSel8Event)
      evLossHistos.fill(HIST("evLoss/hEvent"), 2.5);
    if (isMCSel8Event)
      evLossHistos.fill(HIST("evLoss/hEvent"), 3.5);

    // Loop over all the Generated level particles
    for (const auto& mcPart : mcParticles) {
      if (!mcPart.isPhysicalPrimary())
        continue;
      if (std::abs(mcPart.y()) >= kCfgTpcClasses[0])
        continue;

      const float pt = mcPart.pt();
      const int pdg = mcPart.pdgCode();

      if (pdg == PDGDeuteron) {
        evLossHistos.fill(HIST("evLoss/pt/hDeuteronGen"), pt);
        if (isTvxEvent)
          evLossHistos.fill(HIST("evLoss/pt/hDeuteronTriggeredTVX"), pt);
        if (isSel8Event)
          evLossHistos.fill(HIST("evLoss/pt/hDeuteronTriggeredSel8"), pt);
        if (isMCSel8Event)
          evLossHistos.fill(HIST("evLoss/pt/hDeuteronTriggeredMCSel8"), pt);
      }
      if (pdg == -PDGDeuteron) {
        evLossHistos.fill(HIST("evLoss/pt/hAntiDeuteronGen"), pt);
        if (isTvxEvent)
          evLossHistos.fill(HIST("evLoss/pt/hAntiDeuteronTriggeredTVX"), pt);
        if (isSel8Event)
          evLossHistos.fill(HIST("evLoss/pt/hAntiDeuteronTriggeredSel8"), pt);
        if (isMCSel8Event)
          evLossHistos.fill(HIST("evLoss/pt/hAntiDeuteronTriggeredMCSel8"), pt);
      }
      if (pdg == PDGHelium) {
        evLossHistos.fill(HIST("evLoss/pt/hHeliumGen"), pt);
        if (isTvxEvent)
          evLossHistos.fill(HIST("evLoss/pt/hHeliumTriggeredTVX"), pt);
        if (isSel8Event)
          evLossHistos.fill(HIST("evLoss/pt/hHeliumTriggeredSel8"), pt);
        if (isMCSel8Event)
          evLossHistos.fill(HIST("evLoss/pt/hHeliumTriggeredMCSel8"), pt);
      }
      if (pdg == -PDGHelium) {
        evLossHistos.fill(HIST("evLoss/pt/hAntiHeliumGen"), pt);
        if (isTvxEvent)
          evLossHistos.fill(HIST("evLoss/pt/hAntiHeliumTriggeredTVX"), pt);
        if (isSel8Event)
          evLossHistos.fill(HIST("evLoss/pt/hAntiHeliumTriggeredSel8"), pt);
        if (isMCSel8Event)
          evLossHistos.fill(HIST("evLoss/pt/hAntiHeliumTriggeredMCSel8"), pt);
      }
    }
  }
  PROCESS_SWITCH(LFNucleiBATask, processEvSgLossMC, "process MC SignLoss", false);

  // EVENT LOSS, SIGNAL LOSS and EFFICIENCY CHECKER process function
  void processMCGenLosses(
    soa::Join<aod::McCollisions, aod::McCentFT0Ms>::iterator const& mcCollision,
    const soa::SmallGroups<soa::Join<EventCandidatesMC, o2::aod::PVMults>>& collisions,
    o2::aod::McParticles const& mcParticles)
  {
    bool isINELgt0true = pwglf::isINELgtNmc(mcParticles, 0, pdgDB);

    // EVENT LOSS DENOMINATOR
    // No cuts
    histoGen.fill(HIST("events/hMCGen"), 0.5);
    // Vtz cut
    if (mcCollision.posZ() < cfgVzCutLow || mcCollision.posZ() > cfgVzCutHigh)
      return;
    histoGen.fill(HIST("events/hMCGen"), 1.5);
    // INEL > 0
    if (isINELgt0true)
      histoGen.fill(HIST("events/hMCGen"), 2.5);

    // SIGNAL LOSS DENOMINATOR
    for (const auto& mcParticle : mcParticles) {
      if (mcParticle.y() > kinemOptions.cfgRapidityCutHigh || mcParticle.y() < kinemOptions.cfgRapidityCutLow)
        continue;

      int pdg = mcParticle.pdgCode();
      float pt = mcParticle.pt();
      bool isPhysPrim = mcParticle.isPhysicalPrimary();

      if (enableHe && isPhysPrim && (std::abs(pdg) == PDGHelium)) {
        if (pdg > 0) {
          histoGen.fill(HIST("helium/MCGen/ptGen_INEL_Prim_He"), pt);
          if (enableCentrality)
            histoGen.fill(HIST("helium/MCGen/ptGenVsMult_INEL_Prim_He"), pt, mcCollision.centFT0M());
        } else {
          histoGen.fill(HIST("helium/MCGen/ptGen_INEL_Prim_antiHe"), pt);
          if (enableCentrality)
            histoGen.fill(HIST("helium/MCGen/ptGenVsMult_INEL_Prim_antiHe"), pt, mcCollision.centFT0M());
        }
        if (isINELgt0true) {
          if (pdg > 0) {
            histoGen.fill(HIST("helium/MCGen/ptGen_INELgt0_Prim_He"), pt);
            if (enableCentrality)
              histoGen.fill(HIST("helium/MCGen/ptGenVsMult_INELgt0_Prim_He"), pt, mcCollision.centFT0M());
          } else {
            histoGen.fill(HIST("helium/MCGen/ptGen_INELgt0_Prim_antiHe"), pt);
            if (enableCentrality)
              histoGen.fill(HIST("helium/MCGen/ptGenVsMult_INELgt0_Prim_antiHe"), pt, mcCollision.centFT0M());
          }
        }
      }
    }

    int recoIdxINEL = 0;
    int recoIdxINELgt0 = 0;
    for (const auto& collision : collisions) {
      bool hasTVX = collision.selection_bit(aod::evsel::kIsTriggerTVX);
      bool hasNoTFB = collision.selection_bit(aod::evsel::kNoTimeFrameBorder);
      bool hasNoItsRofFB = collision.selection_bit(aod::evsel::kNoITSROFrameBorder);

      // Check event selection
      histoGen.fill(HIST("events/hMCReco"), 0.5);
      if (evselOptions.useTVXtrigger && !hasTVX)
        continue;
      if (evselOptions.removeTFBorder && !hasNoTFB)
        continue;
      if (evselOptions.removeITSROFBorder && !hasNoItsRofFB)
        continue;
      histoGen.fill(HIST("events/hMCReco"), 1.5);

      recoIdxINEL++;

      if (collision.isInelGt0() && isINELgt0true) {
        histoGen.fill(HIST("events/hMCReco"), 2.5);
        recoIdxINELgt0++;
      }

      for (const auto& mcParticle : mcParticles) {
        if (mcParticle.y() > kinemOptions.cfgRapidityCutHigh || mcParticle.y() < kinemOptions.cfgRapidityCutLow)
          continue;

        int pdg = mcParticle.pdgCode();
        float pt = mcParticle.pt();
        bool isPhysPrim = mcParticle.isPhysicalPrimary();

        if (enableHe && isPhysPrim && (std::abs(pdg) == PDGHelium)) {
          if (pdg > 0) {
            histoGen.fill(HIST("helium/MCReco/ptGen_INEL_Prim_He"), pt);
            if (enableCentrality)
              histoGen.fill(HIST("helium/MCReco/ptGenVsMult_INEL_Prim_He"), pt, mcCollision.centFT0M());
          } else {
            histoGen.fill(HIST("helium/MCReco/ptGen_INEL_Prim_antiHe"), pt);
            if (enableCentrality)
              histoGen.fill(HIST("helium/MCReco/ptGenVsMult_INEL_Prim_antiHe"), pt, mcCollision.centFT0M());
          }
          if (recoIdxINELgt0 > 0) {
            if (pdg > 0) {
              histoGen.fill(HIST("helium/MCReco/ptGen_INELgt0_Prim_He"), pt);
              if (enableCentrality)
                histoGen.fill(HIST("helium/MCReco/ptGenVsMult_INELgt0_Prim_He"), pt, mcCollision.centFT0M());
            } else {
              histoGen.fill(HIST("helium/MCReco/ptGen_INELgt0_Prim_antiHe"), pt);
              if (enableCentrality)
                histoGen.fill(HIST("helium/MCReco/ptGenVsMult_INELgt0_Prim_antiHe"), pt, mcCollision.centFT0M());
            }
          }
        }
      }
    }

    if (recoIdxINEL < 1) {
      return;
    }

    // EVENT LOSS NUMERATOR
    histoGen.fill(HIST("events/hMCGenReco"), 0.5);
    if (recoIdxINELgt0 > 0)
      histoGen.fill(HIST("events/hMCGenReco"), 1.5);

    // SIGNAL LOSS NUMERATOR
    for (const auto& mcParticle : mcParticles) {
      if (mcParticle.y() > kinemOptions.cfgRapidityCutHigh || mcParticle.y() < kinemOptions.cfgRapidityCutLow)
        continue;

      int pdg = mcParticle.pdgCode();
      float pt = mcParticle.pt();
      bool isPhysPrim = mcParticle.isPhysicalPrimary();

      if (enableHe && isPhysPrim && (std::abs(pdg) == PDGHelium)) {
        if (pdg > 0) {
          histoGen.fill(HIST("helium/MCGenReco/ptGen_INEL_Prim_He"), pt);
          if (enableCentrality)
            histoGen.fill(HIST("helium/MCGenReco/ptGenVsMult_INEL_Prim_He"), pt, mcCollision.centFT0M());
        } else {
          histoGen.fill(HIST("helium/MCGenReco/ptGen_INEL_Prim_antiHe"), pt);
          if (enableCentrality)
            histoGen.fill(HIST("helium/MCGenReco/ptGenVsMult_INEL_Prim_antiHe"), pt, mcCollision.centFT0M());
        }
        if (recoIdxINELgt0 > 0) {
          if (pdg > 0) {
            histoGen.fill(HIST("helium/MCGenReco/ptGen_INELgt0_Prim_He"), pt);
            if (enableCentrality)
              histoGen.fill(HIST("helium/MCGenReco/ptGenVsMult_INELgt0_Prim_He"), pt, mcCollision.centFT0M());
          } else {
            histoGen.fill(HIST("helium/MCGenReco/ptGen_INELgt0_Prim_antiHe"), pt);
            if (enableCentrality)
              histoGen.fill(HIST("helium/MCGenReco/ptGenVsMult_INELgt0_Prim_antiHe"), pt, mcCollision.centFT0M());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(LFNucleiBATask, processMCGenLosses, "process MCGen losses", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<LFNucleiBATask>(cfgc)};
}
