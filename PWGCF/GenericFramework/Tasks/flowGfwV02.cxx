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
//
/// \file flowGfwV02.cxx
/// \brief Skeleton copy of flowGfwLightIons with empty function bodies
/// \author Emil Gorm Nielsen, NBI, emil.gorm.nielsen@cern.ch

#include "FlowContainer.h"
#include "FlowPtContainer.h"
#include "GFW.h"
#include "GFWConfig.h"
#include "GFWCumulant.h"
#include "GFWPowerArray.h"
#include "GFWWeights.h"
#include "GFWWeightsList.h"

#include "PWGCF/DataModel/CorrelationsDerived.h"
#include "PWGCF/JCorran/DataModel/JCatalyst.h"

#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"
#include <CCDB/BasicCCDBManager.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DataFormatsParameters/GRPObject.h>

// PID ADD
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"

#include "ReconstructionDataFormats/PID.h"
// PID ADD

#include <TF1.h>
#include <TPDGCode.h>
#include <TProfile.h>
#include <TRandom3.h>

#include <boost/algorithm/string.hpp>

#include <algorithm>
#include <chrono>
#include <complex>
#include <ctime>
#include <experimental/type_traits>
#include <map>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};
static constexpr float LongArrayFloat[3][20] = {{1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2}, {2.1, 2.2, 2.3, -2.1, -2.2, -2.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2}, {3.1, 3.2, 3.3, -3.1, -3.2, -3.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2}};

namespace o2::analysis::gfw
{
std::vector<double> ptbinning = {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.5, 4, 5, 6, 8, 10};
float ptpoilow = 0.2, ptpoiup = 10.0;
float ptreflow = 0.2, ptrefup = 3.0;
float ptlow = 0.2, ptup = 10.0;
int etabins = 16;
float etalow = -0.8, etaup = 0.8;
int vtxZbins = 40;
float vtxZlow = -10.0, vtxZup = 10.0;
int phibins = 72;
float philow = 0.0;
float phiup = o2::constants::math::TwoPI;
int nchbins = 300;
float nchlow = 0;
float nchup = 3000;
std::vector<double> centbinning(90);
int nBootstrap = 10;
std::vector<std::pair<double, double>> etagapsPtPt;
GFWRegions regions;
GFWCorrConfigs configs;
} // namespace o2::analysis::gfw

struct FlowGfwV02 {
  O2_DEFINE_CONFIGURABLE(cfgNbootstrap, int, 10, "Number of subsamples")
  O2_DEFINE_CONFIGURABLE(cfgMpar, int, 4, "Highest order of pt-pt correlations")
  O2_DEFINE_CONFIGURABLE(cfgCentEstimator, int, 0, "0:FT0C; 1:FT0CVariant1; 2:FT0M; 3:FT0A")
  O2_DEFINE_CONFIGURABLE(cfgEfficiency, std::string, "", "CCDB path to efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgAcceptance, std::string, "", "CCDB path to acceptance object")
  O2_DEFINE_CONFIGURABLE(cfgFixedMultMin, int, 1, "Minimum for fixed nch range");
  O2_DEFINE_CONFIGURABLE(cfgFixedMultMax, int, 3000, "Maximum for fixed nch range");
  // PID ADD
  O2_DEFINE_CONFIGURABLE(cfgTofPtCut, float, 0.5f, "Minimum pt to use TOF N-sigma")
  O2_DEFINE_CONFIGURABLE(cfgUseItsPID, bool, true, "Use ITS PID for particle identification")
  O2_DEFINE_CONFIGURABLE(cfgGetNsigmaQA, bool, true, "Get QA histograms for selection of pions, kaons, and protons")

  // PID ADD
  O2_DEFINE_CONFIGURABLE(cfgUseMultiplicityFlowWeights, bool, true, "Enable or disable the use of multiplicity-based event weighting");
  O2_DEFINE_CONFIGURABLE(cfgConsistentEventFlag, int, 15, "Flag for consistent event selection");

  Configurable<GFWBinningCuts> cfgGFWBinning{"cfgGFWBinning", {40, 16, 72, 300, 0, 3000, 0.2, 10.0, 0.2, 5.0, {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.5, 5}, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90}}, "Configuration for binning"};
  Configurable<GFWRegions> cfgRegions{"cfgRegions", {{"refN", "refP", "refFull", "refMid", "piP", "kaP", "prP"}, {-0.8, 0.5, -0.8, -0.4, 0.5, 0.5, 0.5}, {-0.5, 0.8, 0.8, 0.4, 0.8, 0.8, 0.8}, {0, 0, 0, 0, 1, 1, 1}, {1, 1, 1, 1, 1, 1, 1}}, "Configurations for GFW regions"};

  Configurable<GFWCorrConfigs> cfgCorrConfig{"cfgCorrConfig", {{"refP {2} refN {-2}", "piP {2} refN {-2}", "kaP {2} refN {-2}", "prP {2} refN {-2}", "refFull {2 -2}", "refFull {2 -2}", "refFull {2 -2}", "refFull {2 -2}", "refFull {2 -2}", "refFull {2 -2}", "refFull {2 -2}", "refFull {2 -2}", "refFull {2 -2}"}, {"ChGap22", "PiGap22", "KaGap22", "PrGap22", "ChFull22", "nchCh", "nchPi", "nchKa", "nchPr", "v02ptCh", "v02ptPi", "v02ptKa", "v02ptPr"}, {0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1}, {15, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}, "Configurations for each correlation to calculate"};
  Configurable<LabeledArray<float>> nSigmas{"nSigmas", {LongArrayFloat[0], 6, 3, {"UpCut_pi", "UpCut_ka", "UpCut_pr", "LowCut_pi", "LowCut_ka", "LowCut_pr"}, {"TPC", "TOF", "ITS"}}, "Labeled array for n-sigma values for TPC, TOF, ITS for pions, kaons, protons (positive and negative)"};

  struct : ConfigurableGroup {
    Configurable<float> cfgPtMin{"cfgPtMin", 0.2f, "Minimum pT used for track selection."};
    Configurable<float> cfgPtMax{"cfgPtMax", 5.0f, "Maximum pT used for track selection."};
    Configurable<float> cfgEtaMax{"cfgEtaMax", 0.8f, "Maximum eta used for track selection."};
  } cfgTrackCuts;

  struct : ConfigurableGroup {
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", 10.0f, "Maximum primary vertex cut applied for the events."};
    Configurable<int> cfgMultMin{"cfgMultMin", 10, "Minimum number of particles required for the event to have."};
  } cfgEventCuts;

  // // Filters to be applied to the received data.
  // // The analysis assumes the data has been subjected to a QA of its selection,
  // // and thus only the final distributions of the data for analysis are saved.
  o2::framework::expressions::Filter collFilter = (nabs(aod::collision::posZ) < cfgEventCuts.cfgZvtxMax);
  o2::framework::expressions::Filter trackFilter = (aod::track::pt > cfgTrackCuts.cfgPtMin) && (aod::track::pt < cfgTrackCuts.cfgPtMax) && (nabs(aod::track::eta) < cfgTrackCuts.cfgEtaMax);
  o2::framework::expressions::Filter cftrackFilter = (aod::cftrack::pt > cfgTrackCuts.cfgPtMin) && (aod::cftrack::pt < cfgTrackCuts.cfgPtMax); // eta cuts done by jfluc

  //  Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;

  struct Config {
    TH1D* mEfficiency = nullptr;
    GFWWeights* mAcceptance;
    bool correctionsLoaded = false;
  } cfg;

  // Define output
  OutputObj<FlowContainer> fFC{FlowContainer("FlowContainer")};
  HistogramRegistry registry{"registry"};

  GFW* fGFW = new GFW();
  std::vector<GFW::CorrConfig> corrconfigs;

  TRandom3* fRndm = new TRandom3(0);
  TAxis* fSecondAxis;
  int lastRun = -1;

  // region indices for consistency flag
  int posRegionIndex = -1;
  int negRegionIndex = -1;
  int fullRegionIndex = -1;
  int midRegionIndex = -1;
  // PID

  struct PIDState {
    o2::aod::ITSResponse itsResponse;
    std::array<float, 6> tofNsigmaCut;
    std::array<float, 6> itsNsigmaCut;
    std::array<float, 6> tpcNsigmaCut;
    TH1D* hPtMid[4] = {nullptr, nullptr, nullptr, nullptr};
    ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.2, 0.5, 1, 1.5, 2, 3, 4, 6, 10}, "pt axis for histograms"};
    ConfigurableAxis axisNsigmaTPC{"axisNsigmaTPC", {80, -5, 5}, "nsigmaTPC axis"};
    ConfigurableAxis axisNsigmaTOF{"axisNsigmaTOF", {80, -5, 5}, "nsigmaTOF axis"};
    ConfigurableAxis axisNsigmaITS{"axisNsigmaITS", {80, -5, 5}, "nsigmaITS axis"};
    ConfigurableAxis axisTpcSignal{"axisTpcSignal", {250, 0, 250}, "dEdx axis for TPC"};
  };
  PIDState pidStates;

  using GFWTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFbeta, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>>;

  ~FlowGfwV02()
  {
    delete fGFW;
    fGFW = nullptr;
    delete fRndm;
    fRndm = nullptr;
    delete fSecondAxis;
    fSecondAxis = nullptr;
  }

  // PID ADD
  enum PIDIndex {
    kCharged = 0,
    kPions,
    kKaons,
    kProtons
  };
  enum PiKpArrayIndex {
    iPionUp = 0,
    iKaonUp,
    iProtonUp,
    iPionLow,
    iKaonLow,
    iProtonLow
  };
  enum DetectorType {
    kTPC = 0,
    kTOF,
    kITS
  };

  // PID ADD

  void init(InitContext const&)
  {

    // PID ADD
    pidStates.tpcNsigmaCut[iPionUp] = nSigmas->getData()[iPionUp][kTPC];
    pidStates.tpcNsigmaCut[iKaonUp] = nSigmas->getData()[iKaonUp][kTPC];
    pidStates.tpcNsigmaCut[iProtonUp] = nSigmas->getData()[iProtonUp][kTPC];
    pidStates.tpcNsigmaCut[iPionLow] = nSigmas->getData()[iPionLow][kTPC];
    pidStates.tpcNsigmaCut[iKaonLow] = nSigmas->getData()[iKaonLow][kTPC];
    pidStates.tpcNsigmaCut[iProtonLow] = nSigmas->getData()[iProtonLow][kTPC];

    pidStates.tofNsigmaCut[iPionUp] = nSigmas->getData()[iPionUp][kTOF];
    pidStates.tofNsigmaCut[iKaonUp] = nSigmas->getData()[iKaonUp][kTOF];
    pidStates.tofNsigmaCut[iProtonUp] = nSigmas->getData()[iProtonUp][kTOF];
    pidStates.tofNsigmaCut[iPionLow] = nSigmas->getData()[iPionLow][kTOF];
    pidStates.tofNsigmaCut[iKaonLow] = nSigmas->getData()[iKaonLow][kTOF];
    pidStates.tofNsigmaCut[iProtonLow] = nSigmas->getData()[iProtonLow][kTOF];

    pidStates.itsNsigmaCut[iPionUp] = nSigmas->getData()[iPionUp][kITS];
    pidStates.itsNsigmaCut[iKaonUp] = nSigmas->getData()[iKaonUp][kITS];
    pidStates.itsNsigmaCut[iProtonUp] = nSigmas->getData()[iProtonUp][kITS];
    pidStates.itsNsigmaCut[iPionLow] = nSigmas->getData()[iPionLow][kITS];
    pidStates.itsNsigmaCut[iKaonLow] = nSigmas->getData()[iKaonLow][kITS];
    pidStates.itsNsigmaCut[iProtonLow] = nSigmas->getData()[iProtonLow][kITS];

    if (cfgGetNsigmaQA) {
      if (!cfgUseItsPID) {
        registry.add("TofTpcNsigma_before", "", {HistType::kTHnSparseD, {{pidStates.axisNsigmaTPC, pidStates.axisNsigmaTOF, pidStates.axisPt}}});
        registry.add("TofTpcNsigma_after", "", {HistType::kTHnSparseD, {{pidStates.axisNsigmaTPC, pidStates.axisNsigmaTOF, pidStates.axisPt}}});
      }
      if (cfgUseItsPID) {
        registry.add("TofItsNsigma_before", "", {HistType::kTHnSparseD, {{pidStates.axisNsigmaITS, pidStates.axisNsigmaTOF, pidStates.axisPt}}});
        registry.add("TofItsNsigma_after", "", {HistType::kTHnSparseD, {{pidStates.axisNsigmaITS, pidStates.axisNsigmaTOF, pidStates.axisPt}}});
      }

      registry.add("TpcdEdx_ptwise", "", {HistType::kTH2D, {{pidStates.axisTpcSignal, pidStates.axisPt}}});
      registry.add("TpcdEdx_ptwise_afterCut", "", {HistType::kTH2D, {{pidStates.axisTpcSignal, pidStates.axisPt}}});
    }
    // PID ADD

    o2::analysis::gfw::regions.SetNames(cfgRegions->GetNames());
    o2::analysis::gfw::regions.SetEtaMin(cfgRegions->GetEtaMin());
    o2::analysis::gfw::regions.SetEtaMax(cfgRegions->GetEtaMax());
    o2::analysis::gfw::regions.SetpTDifs(cfgRegions->GetpTDifs());
    o2::analysis::gfw::regions.SetBitmasks(cfgRegions->GetBitmasks());
    o2::analysis::gfw::configs.SetCorrs(cfgCorrConfig->GetCorrs());
    o2::analysis::gfw::configs.SetHeads(cfgCorrConfig->GetHeads());
    o2::analysis::gfw::configs.SetpTDifs(cfgCorrConfig->GetpTDifs());
    o2::analysis::gfw::configs.SetpTCorrMasks(cfgCorrConfig->GetpTCorrMasks());
    o2::analysis::gfw::regions.Print();
    o2::analysis::gfw::configs.Print();
    o2::analysis::gfw::ptbinning = cfgGFWBinning->GetPtBinning();
    o2::analysis::gfw::ptpoilow = cfgGFWBinning->GetPtPOImin();
    o2::analysis::gfw::ptpoiup = cfgGFWBinning->GetPtPOImax();
    o2::analysis::gfw::ptreflow = cfgGFWBinning->GetPtRefMin();
    o2::analysis::gfw::ptrefup = cfgGFWBinning->GetPtRefMax();
    o2::analysis::gfw::ptlow = cfgTrackCuts.cfgPtMin;
    o2::analysis::gfw::ptup = cfgTrackCuts.cfgPtMax;
    o2::analysis::gfw::etabins = cfgGFWBinning->GetEtaBins();
    o2::analysis::gfw::vtxZbins = cfgGFWBinning->GetVtxZbins();
    o2::analysis::gfw::phibins = cfgGFWBinning->GetPhiBins();
    o2::analysis::gfw::philow = 0.0f;
    o2::analysis::gfw::phiup = o2::constants::math::TwoPI;
    o2::analysis::gfw::nchbins = cfgGFWBinning->GetNchBins();
    o2::analysis::gfw::nchlow = cfgGFWBinning->GetNchMin();
    o2::analysis::gfw::nchup = cfgGFWBinning->GetNchMax();
    o2::analysis::gfw::centbinning = cfgGFWBinning->GetCentBinning();
    cfgGFWBinning->Print();

    // Initialise pt spectra histograms for different particles
    pidStates.hPtMid[kCharged] = new TH1D("hPtMid_charged", "hPtMid_charged", o2::analysis::gfw::ptbinning.size() - 1, &o2::analysis::gfw::ptbinning[0]);
    pidStates.hPtMid[kPions] = new TH1D("hPtMid_pions", "hPtMid_pions", o2::analysis::gfw::ptbinning.size() - 1, &o2::analysis::gfw::ptbinning[0]);
    pidStates.hPtMid[kKaons] = new TH1D("hPtMid_kaons", "hPtMid_kaons", o2::analysis::gfw::ptbinning.size() - 1, &o2::analysis::gfw::ptbinning[0]);
    pidStates.hPtMid[kProtons] = new TH1D("hPtMid_protons", "hPtMid_protons", o2::analysis::gfw::ptbinning.size() - 1, &o2::analysis::gfw::ptbinning[0]);
    pidStates.hPtMid[kCharged]->SetDirectory(nullptr);
    pidStates.hPtMid[kPions]->SetDirectory(nullptr);
    pidStates.hPtMid[kKaons]->SetDirectory(nullptr);
    pidStates.hPtMid[kProtons]->SetDirectory(nullptr);

    AxisSpec phiAxis = {o2::analysis::gfw::phibins, o2::analysis::gfw::philow, o2::analysis::gfw::phiup, "#phi"};
    AxisSpec etaAxis = {o2::analysis::gfw::etabins, -cfgTrackCuts.cfgEtaMax, cfgTrackCuts.cfgEtaMax, "#eta"};
    AxisSpec vtxAxis = {o2::analysis::gfw::vtxZbins, -cfgEventCuts.cfgZvtxMax, cfgEventCuts.cfgZvtxMax, "Vtx_{z} (cm)"};
    AxisSpec ptAxis = {o2::analysis::gfw::ptbinning, "#it{p}_{T} GeV/#it{c}"};

    std::string sCentralityEstimator = "FT0C centrality (%)";
    AxisSpec centAxis = {o2::analysis::gfw::centbinning, sCentralityEstimator.c_str()};

    std::vector<double> nchbinning;
    int nchskip = (o2::analysis::gfw::nchup - o2::analysis::gfw::nchlow) / o2::analysis::gfw::nchbins;
    for (int i = 0; i <= o2::analysis::gfw::nchbins; ++i) {
      nchbinning.push_back(nchskip * i + o2::analysis::gfw::nchlow + 0.5);
    }
    AxisSpec nchAxis = {nchbinning, "N_{ch}"};
    registry.add("v02pt", "", {HistType::kTProfile2D, {ptAxis, centAxis}});
    registry.add("nchMid", "", {HistType::kTProfile2D, {ptAxis, centAxis}});

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    int ptbins = o2::analysis::gfw::ptbinning.size() - 1;
    fSecondAxis = new TAxis(ptbins, &o2::analysis::gfw::ptbinning[0]);

    // QA histograms
    registry.add("trackQA/before/phi_eta_vtxZ", "", {HistType::kTH3D, {phiAxis, etaAxis, vtxAxis}});
    registry.add("trackQA/before/nch_pt", "#it{p}_{T} vs multiplicity; N_{ch}; #it{p}_{T}", {HistType::kTH2D, {nchAxis, ptAxis}});
    registry.addClone("trackQA/before/", "trackQA/after/");
    registry.add("trackQA/after/pt_ref", "", {HistType::kTH1D, {{100, o2::analysis::gfw::ptreflow, o2::analysis::gfw::ptrefup}}});
    registry.add("trackQA/after/pt_poi", "", {HistType::kTH1D, {{100, o2::analysis::gfw::ptpoilow, o2::analysis::gfw::ptpoiup}}});

    registry.add("eventQA/before/multiplicity", "", {HistType::kTH1D, {nchAxis}});
    registry.add("eventQA/before/centrality", "", {HistType::kTH1D, {centAxis}});
    registry.addClone("eventQA/before/", "eventQA/after/");

    if (o2::analysis::gfw::regions.GetSize() < 0)
      LOGF(error, "Configuration contains vectors of different size - check the GFWRegions configurable");
    for (auto i(0); i < o2::analysis::gfw::regions.GetSize(); ++i) {
      fGFW->AddRegion(o2::analysis::gfw::regions.GetNames()[i], o2::analysis::gfw::regions.GetEtaMin()[i], o2::analysis::gfw::regions.GetEtaMax()[i], (o2::analysis::gfw::regions.GetpTDifs()[i]) ? ptbins + 1 : 1, o2::analysis::gfw::regions.GetBitmasks()[i]);
    }
    for (auto i = 0; i < o2::analysis::gfw::configs.GetSize(); ++i) {
      corrconfigs.push_back(fGFW->GetCorrelatorConfig(o2::analysis::gfw::configs.GetCorrs()[i], o2::analysis::gfw::configs.GetHeads()[i], o2::analysis::gfw::configs.GetpTDifs()[i]));
    }
    if (corrconfigs.empty())
      LOGF(error, "Configuration contains vectors of different size - check the GFWCorrConfig configurable");
    fGFW->CreateRegions();
    TObjArray* oba = new TObjArray();
    addConfigObjectsToObjArray(oba, corrconfigs);
    LOGF(info, "Number of correlators: %d", oba->GetEntries());
    fFC->SetName("FlowContainer");
    fFC->SetXAxis(fSecondAxis);
    fFC->Initialize(oba, centAxis, cfgNbootstrap);
    delete oba;

    if (cfgConsistentEventFlag) {
      posRegionIndex = [&]() {
        auto begin = cfgRegions->GetNames().begin();
        auto end = cfgRegions->GetNames().end();
        auto it = std::find(begin, end, "refP");
        return (it != end) ? std::distance(begin, it) : -1;
      }();
      negRegionIndex = [&]() {
        auto begin = cfgRegions->GetNames().begin();
        auto end = cfgRegions->GetNames().end();
        auto it = std::find(begin, end, "refN");
        return (it != end) ? std::distance(begin, it) : -1;
      }();
      fullRegionIndex = [&]() {
        auto begin = cfgRegions->GetNames().begin();
        auto end = cfgRegions->GetNames().end();
        auto it = std::find(begin, end, "refFull");
        return (it != end) ? std::distance(begin, it) : -1;
      }();
      midRegionIndex = [&]() {
        auto begin = cfgRegions->GetNames().begin();
        auto end = cfgRegions->GetNames().end();
        auto it = std::find(begin, end, "refMid");
        return (it != end) ? std::distance(begin, it) : -1;
      }();
    }
  }

  static constexpr std::string_view FillTimeName[] = {"before/", "after/"};
  enum QAFillTime {
    kBefore,
    kAfter
  };

  void addConfigObjectsToObjArray(TObjArray* oba, const std::vector<GFW::CorrConfig>& configs)
  {
    for (auto it = configs.begin(); it != configs.end(); ++it) {
      if (it->pTDif) {
        std::string suffix = "_ptDiff";
        for (auto i = 0; i < fSecondAxis->GetNbins(); ++i) {
          std::string index = Form("_pt_%i", i + 1);
          oba->Add(new TNamed(it->Head.c_str() + index, it->Head.c_str() + suffix));
        }
      } else {
        oba->Add(new TNamed(it->Head.c_str(), it->Head.c_str()));
      }
    }
  }

  // PID ADD
  template <typename TTrack>
  int getNsigmaPID(TTrack track)
  {
    // Computing Nsigma arrays for pion, kaon, and protons
    std::array<float, 3> nSigmaTPC = {track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr()};
    std::array<float, 3> nSigmaTOF = {track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr()};
    std::array<float, 3> nSigmaITS = {pidStates.itsResponse.nSigmaITS<o2::track::PID::Pion>(track), pidStates.itsResponse.nSigmaITS<o2::track::PID::Kaon>(track), pidStates.itsResponse.nSigmaITS<o2::track::PID::Proton>(track)};
    int pid = -1; // -1 = not identified, 1 = pion, 2 = kaon, 3 = proton

    std::array<float, 3> nSigmaToUse = cfgUseItsPID ? nSigmaITS : nSigmaTPC;                                 // Choose which nSigma to use: TPC or ITS
    std::array<float, 6> detectorNsigmaCut = cfgUseItsPID ? pidStates.itsNsigmaCut : pidStates.tpcNsigmaCut; // Choose which nSigma to use: TPC or ITS

    bool isPion, isKaon, isProton;
    bool isDetectedPion = nSigmaToUse[iPionUp] < detectorNsigmaCut[iPionUp] && nSigmaToUse[iPionUp] > detectorNsigmaCut[iPionLow];
    bool isDetectedKaon = nSigmaToUse[iKaonUp] < detectorNsigmaCut[iKaonUp] && nSigmaToUse[iKaonUp] > detectorNsigmaCut[iKaonLow];
    bool isDetectedProton = nSigmaToUse[iProtonUp] < detectorNsigmaCut[iProtonUp] && nSigmaToUse[iProtonUp] > detectorNsigmaCut[iProtonLow];

    bool isTofPion = nSigmaTOF[iPionUp] < pidStates.tofNsigmaCut[iPionUp] && nSigmaTOF[iPionUp] > pidStates.tofNsigmaCut[iPionLow];
    bool isTofKaon = nSigmaTOF[iKaonUp] < pidStates.tofNsigmaCut[iKaonUp] && nSigmaTOF[iKaonUp] > pidStates.tofNsigmaCut[iKaonLow];
    bool isTofProton = nSigmaTOF[iProtonUp] < pidStates.tofNsigmaCut[iProtonUp] && nSigmaTOF[iProtonUp] > pidStates.tofNsigmaCut[iProtonLow];

    if (track.pt() > cfgTofPtCut && !track.hasTOF()) {
      return -1;
    } else if (track.pt() > cfgTofPtCut && track.hasTOF()) {
      isPion = isTofPion && isDetectedPion;
      isKaon = isTofKaon && isDetectedKaon;
      isProton = isTofProton && isDetectedProton;
    } else {
      isPion = isDetectedPion;
      isKaon = isDetectedKaon;
      isProton = isDetectedProton;
    }

    if ((isPion && isKaon) || (isPion && isProton) || (isKaon && isProton)) {
      return -1; // more than one particle satisfy the criteria
    }

    if (isPion) {
      pid = kPions;
    } else if (isKaon) {
      pid = kKaons;
    } else if (isProton) {
      pid = kProtons;
    } else {
      return -1; // no particle satisfies the criteria
    }

    return pid; // -1 = not identified, 1 = pion, 2 = kaon, 3 = proton
  }
  // PID ADD
  void loadCorrections(aod::BCsWithTimestamps::iterator const& bc)
  {
    uint64_t timestamp = bc.timestamp();
    if (cfg.correctionsLoaded)
      return;
    if (!cfgAcceptance.value.empty()) {
      cfg.mAcceptance = ccdb->getForRun<GFWWeights>(cfgAcceptance.value, timestamp);
    }
    if (!cfgEfficiency.value.empty()) {
      cfg.mEfficiency = ccdb->getForTimeStamp<TH1D>(cfgEfficiency, timestamp);
      if (cfg.mEfficiency == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram from %s", cfgEfficiency.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram from %s (%p)", cfgEfficiency.value.c_str(), (void*)cfg.mEfficiency);
    }
    cfg.correctionsLoaded = true;
  }

  void loadCorrections(int runnumber)
  {
    if (cfg.correctionsLoaded)
      return;
    if (!cfgAcceptance.value.empty()) {
      cfg.mAcceptance = ccdb->getForRun<GFWWeights>(cfgAcceptance.value, runnumber);
    }
    if (!cfgEfficiency.value.empty()) {
      cfg.mEfficiency = ccdb->getForRun<TH1D>(cfgEfficiency.value, runnumber);
    }
    cfg.correctionsLoaded = true;
  }

  template <typename TTrack>
  double getAcceptance(TTrack track, const double& /*vtxz*/, const int& /*pidInd*/ = 0)
  {
    double wacc = 1;
    if constexpr (requires { track.weightNUA(); })
      wacc = 1. / track.weightNUA();
    return wacc;
  }

  template <typename TTrack>
  double getEfficiency(TTrack track, const int& /*pidInd*/ = 0)
  {
    double eff = 1.;
    if constexpr (requires { track.weightEff(); })
      eff = track.weightEff();
    return eff;
  }

  // Define the data type
  enum DataType {
    kReco,
    kGen
  };

  int getPIDIndex(const std::string& corrconfig)
  {
    if (boost::ifind_first(corrconfig, "pi"))
      return kPions;
    if (boost::ifind_first(corrconfig, "ka"))
      return kKaons;
    if (boost::ifind_first(corrconfig, "pr"))
      return kProtons;
    return kCharged;
  }

  GFW::CorrConfig getRelevantCorrName(const int& pidInd)
  {
    if (pidInd == kPions)
      return fGFW->GetCorrelatorConfig("piP {2} refN {-2}", "PiGap22", kFALSE);
    if (pidInd == kKaons)
      return fGFW->GetCorrelatorConfig("kaP {2} refN {-2}", "KaGap22", kFALSE);
    if (pidInd == kProtons)
      return fGFW->GetCorrelatorConfig("prP {2} refN {-2}", "PrGap22", kFALSE);
    return fGFW->GetCorrelatorConfig("refP {2} refN {-2}", "ChGap22", kFALSE);
  }

  template <DataType dt>
  void fillOutputContainers(const float& centmult, const double& rndm, const int& /*run*/ = 0)
  {
    for (uint l_ind = 0; l_ind < corrconfigs.size(); ++l_ind) {
      if (!corrconfigs.at(l_ind).pTDif) {
        auto dnx = fGFW->Calculate(corrconfigs.at(l_ind), 0, kTRUE).real();
        if (dnx == 0)
          continue;
        auto val = fGFW->Calculate(corrconfigs.at(l_ind), 0, kFALSE).real() / dnx;

        if (std::abs(val) < 1) {
          fFC->FillProfile(corrconfigs.at(l_ind).Head.c_str(), centmult, val, (cfgUseMultiplicityFlowWeights) ? dnx : 1.0, rndm);
        }
        continue;
      }

      // Fill pt profiles for different particles
      int pidInd = getPIDIndex(corrconfigs.at(l_ind).Head.c_str());

      // Find the corresponding non-pT-differential correlation configuration
      GFW::CorrConfig corrName = getRelevantCorrName(pidInd); // May be used later for QA

      auto dnx = fGFW->Calculate(corrconfigs.at(0), 0, kTRUE).real();
      if (dnx == 0)
        continue;
      auto val = fGFW->Calculate(corrconfigs.at(0), 0, kFALSE).real() / dnx;
      for (int i = 1; i <= fSecondAxis->GetNbins(); i++) {
        if (corrconfigs.at(l_ind).Head.find("nch") != std::string::npos)
          val = 1.0;
        double ptFraction = 0;
        if (pidStates.hPtMid[pidInd]->Integral() > 0) {
          ptFraction = pidStates.hPtMid[pidInd]->GetBinContent(i) / pidStates.hPtMid[pidInd]->Integral();
          if (std::abs(val) < 1.01)
            fFC->FillProfile(Form("%s_pt_%i", corrconfigs.at(l_ind).Head.c_str(), i), centmult, val * ptFraction, (cfgUseMultiplicityFlowWeights) ? dnx : 1.0, rndm);
        }
      }
    }
    // Fill the profiles for each pT bin
    // printf("Config name: %s\n", corrconfigs.at(0).Head.c_str());
    auto dnx = fGFW->Calculate(corrconfigs.at(0), 0, kTRUE).real();
    if (dnx == 0)
      return;
    auto val = fGFW->Calculate(corrconfigs.at(0), 0, kFALSE).real() / dnx;
    for (int i = 1; i <= fSecondAxis->GetNbins(); i++) {
      double ptFraction = 0;
      if (pidStates.hPtMid[kCharged]->Integral() > 0) {
        ptFraction = pidStates.hPtMid[kCharged]->GetBinContent(i) / pidStates.hPtMid[kCharged]->Integral();
        if (std::abs(val) < 1)
          registry.fill(HIST("v02pt"), fSecondAxis->GetBinCenter(i), centmult, val * ptFraction, (cfgUseMultiplicityFlowWeights) ? dnx : 1.0);
        // printf("bincenter hPtMid: %f, fsecondaxis: %f\n", hPtMid->GetBinCenter(i), fSecondAxis->GetBinCenter(i));
        registry.fill(HIST("nchMid"), fSecondAxis->GetBinCenter(i), centmult, ptFraction);
      }
    }
    return;
  }

  struct XAxis {
    float centrality;
    int64_t multiplicity;
    double time;
  };

  struct AcceptedTracks {
    int nPos;
    int nNeg;
    int nFull;
    int nMid;
  };

  template <DataType dt, typename TCollision, typename TTracks>
  void processCollision(TCollision collision, TTracks tracks, const XAxis& xaxis, const int& run)
  {
    float vtxz = collision.posZ();
    if (tracks.size() < 1)
      return;
    if (xaxis.centrality >= 0 && (xaxis.centrality < o2::analysis::gfw::centbinning.front() || xaxis.centrality > o2::analysis::gfw::centbinning.back()))
      return;
    if (xaxis.multiplicity < cfgFixedMultMin || xaxis.multiplicity > cfgFixedMultMax)
      return;
    fGFW->Clear();
    pidStates.hPtMid[kCharged]->Reset();
    pidStates.hPtMid[kPions]->Reset();
    pidStates.hPtMid[kKaons]->Reset();
    pidStates.hPtMid[kProtons]->Reset();

    float lRandom = fRndm->Rndm();

    // Loop over tracks and check if they are accepted
    AcceptedTracks acceptedTracks{0, 0, 0, 0};
    for (const auto& track : tracks) {
      processTrack(track, vtxz, xaxis.multiplicity, run, acceptedTracks);
      pidStates.hPtMid[kCharged]->Fill(track.pt(), getEfficiency(track));
      // If PID is identified, fill pt spectrum for the corresponding particle
      int pidInd = getNsigmaPID(track);
      if (pidInd != -1 && track.eta() > -0.4 && track.eta() < 0.4) {
        pidStates.hPtMid[pidInd]->Fill(track.pt(), getEfficiency(track, pidInd));
      }
    }
    if (cfgConsistentEventFlag & 1)
      if (!acceptedTracks.nPos || !acceptedTracks.nNeg)
        return;
    if (cfgConsistentEventFlag & 2)
      if (acceptedTracks.nFull < 4) // o2-linter: disable=magic-number (at least four tracks in full acceptance)
        return;
    if (cfgConsistentEventFlag & 4)
      if (acceptedTracks.nPos < 2 || acceptedTracks.nNeg < 2) // o2-linter: disable=magic-number (at least two tracks in each subevent)
        return;
    if (cfgConsistentEventFlag & 8)
      if (acceptedTracks.nPos < 2 || acceptedTracks.nMid < 2 || acceptedTracks.nNeg < 2) // o2-linter: disable=magic-number (at least two tracks in all three subevents)
        return;
    // Fill output containers
    fillOutputContainers<dt>(xaxis.centrality, lRandom, run);
  }

  template <typename TTrack>
  void fillAcceptedTracks(TTrack track, AcceptedTracks& acceptedTracks)
  {
    if (posRegionIndex >= 0 && track.eta() > o2::analysis::gfw::regions.GetEtaMin()[posRegionIndex] && track.eta() < o2::analysis::gfw::regions.GetEtaMax()[posRegionIndex])
      ++acceptedTracks.nPos;
    if (negRegionIndex >= 0 && track.eta() > o2::analysis::gfw::regions.GetEtaMin()[negRegionIndex] && track.eta() < o2::analysis::gfw::regions.GetEtaMax()[negRegionIndex])
      ++acceptedTracks.nNeg;
    if (fullRegionIndex >= 0 && track.eta() > o2::analysis::gfw::regions.GetEtaMin()[fullRegionIndex] && track.eta() < o2::analysis::gfw::regions.GetEtaMax()[fullRegionIndex])
      ++acceptedTracks.nFull;
    if (midRegionIndex >= 0 && track.eta() > o2::analysis::gfw::regions.GetEtaMin()[midRegionIndex] && track.eta() < o2::analysis::gfw::regions.GetEtaMax()[midRegionIndex])
      ++acceptedTracks.nMid;
  }

  template <typename TTrack>
  inline void processTrack(TTrack const& track, const float& vtxz, const int& multiplicity, const int& /*run*/, AcceptedTracks& acceptedTracks)
  {
    // fillPtSums<kReco>(track); // Fill pT sums
    fillTrackQA<kBefore>(track, vtxz);
    registry.fill(HIST("trackQA/before/nch_pt"), multiplicity, track.pt());

    fillGFW<kReco>(track, vtxz);               // Fill GFW
    fillAcceptedTracks(track, acceptedTracks); // Fill accepted tracks
    fillTrackQA<kAfter>(track, vtxz);
    registry.fill(HIST("trackQA/after/nch_pt"), multiplicity, track.pt());
  }

  template <DataType dt, typename TTrack>
  inline void fillGFW(TTrack track, const double& vtxz)
  {
    int pidInd = getNsigmaPID(track);

    // PID ADD

    bool withinPtRef = (track.pt() > o2::analysis::gfw::ptreflow && track.pt() < o2::analysis::gfw::ptrefup);
    bool withinPtPOI = (track.pt() > o2::analysis::gfw::ptpoilow && track.pt() < o2::analysis::gfw::ptpoiup);

    if (!withinPtPOI && !withinPtRef)
      return;
    double weff = getEfficiency(track, pidInd);
    if (weff < 0)
      return;

    double wacc = getAcceptance(track, vtxz, pidInd);

    // Fill cumulants for different particles
    // ***Need to add proper weights for each particle!***
    if (withinPtRef)
      fGFW->Fill(track.eta(), fSecondAxis->FindBin(track.pt()) - 1, track.phi(), weff * wacc, 0);
    if (withinPtPOI && pidInd == kPions)
      fGFW->Fill(track.eta(), fSecondAxis->FindBin(track.pt()) - 1, track.phi(), weff * wacc, kPions);
    if (withinPtPOI && pidInd == kKaons)
      fGFW->Fill(track.eta(), fSecondAxis->FindBin(track.pt()) - 1, track.phi(), weff * wacc, kKaons);
    if (withinPtPOI && pidInd == kProtons)
      fGFW->Fill(track.eta(), fSecondAxis->FindBin(track.pt()) - 1, track.phi(), weff * wacc, kProtons);
    return;
  }

  template <QAFillTime ft, typename TTrack>
  inline void fillTrackQA(TTrack track, const float vtxz)
  {
    double wacc = getAcceptance(track, vtxz);
    registry.fill(HIST("trackQA/") + HIST(FillTimeName[ft]) + HIST("phi_eta_vtxZ"), track.phi(), track.eta(), vtxz, (ft == kAfter) ? wacc : 1.0);
    if (ft == kAfter) {
      registry.fill(HIST("trackQA/") + HIST(FillTimeName[ft]) + HIST("pt_ref"), track.pt());
      registry.fill(HIST("trackQA/") + HIST(FillTimeName[ft]) + HIST("pt_poi"), track.pt());
    }
    return;
  }

  double getTimeSinceStartOfFill(uint64_t, int) { return 0.0; }

  void processData(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::CentFT0Ms, aod::CentFV0As, aod::CentNTPVs, aod::CentNGlobals, aod::CentMFTs>>::iterator const& collision, aod::BCsWithTimestamps const&, GFWTracks const& tracks)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int run = bc.runNumber();
    if (run != lastRun) {
      lastRun = run;
      LOGF(info, "run = %d", run);
    }

    loadCorrections(bc);
    const XAxis xaxis{collision.centFT0C(), tracks.size(), -1.0};
    processCollision<kReco>(collision, tracks, xaxis, run);
  }
  PROCESS_SWITCH(FlowGfwV02, processData, "Process analysis for non-derived data", true);

  void processCFDerived(aod::CFCollision const& collision, soa::Filtered<aod::CFTracks> const& tracks)
  {
    int run = collision.runNumber();
    if (run != lastRun) {
      lastRun = run;
      LOGF(info, "run = %d", run);
    }
    loadCorrections(run);
    const XAxis xaxis{collision.multiplicity(), tracks.size(), -1.0};

    registry.fill(HIST("eventQA/after/centrality"), xaxis.centrality);
    registry.fill(HIST("eventQA/after/multiplicity"), xaxis.multiplicity);

    // processCollision<kReco>(collision, tracks, xaxis, run);
  }
  PROCESS_SWITCH(FlowGfwV02, processCFDerived, "Process analysis for CF derived data", false);
  void processCFDerivedCorrected(aod::CFCollision const& collision, soa::Filtered<soa::Join<aod::CFTracks, aod::JWeights>> const& tracks)
  {
    int run = collision.runNumber();
    if (run != lastRun) {
      lastRun = run;
      LOGF(info, "run = %d", run);
    }
    const XAxis xaxis{collision.multiplicity(), tracks.size(), -1.0};
    registry.fill(HIST("eventQA/after/centrality"), xaxis.centrality);
    registry.fill(HIST("eventQA/after/multiplicity"), xaxis.multiplicity);
    // processCollision<kReco>(collision, tracks, xaxis, run);
  }
  PROCESS_SWITCH(FlowGfwV02, processCFDerivedCorrected, "Process analysis for CF derived data with corrections", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowGfwV02>(cfgc),
  };
}
