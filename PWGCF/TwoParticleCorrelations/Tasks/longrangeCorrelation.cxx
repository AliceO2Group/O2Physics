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
/// \file longrangeCorrelation.cxx
///
/// \brief task for long range correlation analysis
/// \author Abhi Modak (abhi.modak@cern.ch) and Debojit Sarkar (debojit.sarkar@cern.ch)
/// \since April 22, 2025

#include "PWGCF/Core/CorrelationContainer.h"
#include "PWGCF/Core/PairCuts.h"
#include "PWGCF/DataModel/CorrelationsDerived.h"
#include "PWGMM/Mult/DataModel/bestCollisionTable.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "CommonConstants/MathConstants.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DetectorsCommonDataFormats/AlignParam.h"
#include "FT0Base/Geometry.h"
#include "FV0Base/Geometry.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/Track.h"

#include <TComplex.h>
#include <TH1F.h>
#include <TMath.h>
#include <TPDGCode.h>

#include <chrono>
#include <cstdio>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::track;
using namespace o2::aod::fwdtrack;
using namespace o2::aod::evsel;
using namespace o2::constants::math;

static constexpr TrackSelectionFlags::flagtype TrackSelectionIts =
  TrackSelectionFlags::kITSNCls | TrackSelectionFlags::kITSChi2NDF |
  TrackSelectionFlags::kITSHits;
static constexpr TrackSelectionFlags::flagtype TrackSelectionTpc =
  TrackSelectionFlags::kTPCNCls |
  TrackSelectionFlags::kTPCCrossedRowsOverNCls |
  TrackSelectionFlags::kTPCChi2NDF;
static constexpr TrackSelectionFlags::flagtype TrackSelectionDca =
  TrackSelectionFlags::kDCAz | TrackSelectionFlags::kDCAxy;
static constexpr TrackSelectionFlags::flagtype TrackSelectionDcaxyOnly =
  TrackSelectionFlags::kDCAxy;

enum KindOfEvntType {
  kSE,
  kME
};

enum KindOfCorrType {
  kFT0AGLOBAL,
  kFT0CGLOBAL,
  kMFTGLOBAL,
  kFT0AMFT,
  kFT0AFT0C
};

enum KindOfParticles {
  PIONS,
  KAONS,
  PROTONS
};

static constexpr std::string_view kCorrType[] = {
  "Ft0aGlobal/",
  "Ft0cGlobal/",
  "MftGlobal/",
  "Ft0aMft/",
  "Ft0aFt0c/"};
static constexpr std::string_view kEvntType[] = {"SE/", "ME/"};
auto static constexpr kMinFt0cCell = 96;
auto static constexpr kMinCharge = 3.f;
AxisSpec axisEvent{10, 0.5, 9.5, "#Event", "EventAxis"};

struct LongrangeCorrelation {

  struct : ConfigurableGroup {
    Configurable<std::string> cfgURL{"cfgURL", "http://alice-ccdb.cern.ch", "Address of the CCDB to browse"};
    Configurable<int64_t> noLaterThan{"noLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "Latest acceptable timestamp of creation for the object"};
  } cfgCcdbParam;

  SliceCache cache;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Service<o2::framework::O2DatabasePDG> pdg;
  o2::ccdb::CcdbApi ccdbApi;
  o2::ft0::Geometry ft0Det;
  std::vector<o2::detectors::AlignParam>* offsetFT0;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<float> cfgVtxCut{"cfgVtxCut", 10.0f, "Vertex Z range to consider"};
  Configurable<float> cfgEtaCut{"cfgEtaCut", 1.0f, "Eta range to consider"};
  Configurable<float> dcaZ{"dcaZ", 0.2f, "Custom DCA Z cut (ignored if negative)"};
  Configurable<float> cfgPtCutMin{"cfgPtCutMin", 0.2f, "minimum accepted track pT"};
  Configurable<float> cfgPtCutMax{"cfgPtCutMax", 10.0f, "maximum accepted track pT"};
  Configurable<int> mixingParameter{"mixingParameter", 5, "how many events are mixed"};
  Configurable<float> cfigMftEtaMax{"cfigMftEtaMax", -2.5f, "Maximum MFT eta cut"};
  Configurable<float> cfigMftEtaMin{"cfigMftEtaMin", -3.6f, "Minimum MFT eta cut"};
  Configurable<float> cfigMftDcaxy{"cfigMftDcaxy", 2.0f, "cut on DCA xy for MFT tracks"};
  Configurable<int> cfigMftCluster{"cfigMftCluster", 5, "cut on MFT Cluster"};
  Configurable<double> cfgSampleSize{"cfgSampleSize", 10, "Sample size for mixed event"};
  Configurable<bool> isApplySameBunchPileup{"isApplySameBunchPileup", false, "Enable SameBunchPileup cut"};
  Configurable<bool> isApplyGoodZvtxFT0vsPV{"isApplyGoodZvtxFT0vsPV", false, "Enable GoodZvtxFT0vsPV cut"};
  Configurable<bool> isReadoutCenter{"isReadoutCenter", false, "Enable Readout Center"};
  Configurable<bool> isUseEffCorr{"isUseEffCorr", false, "Enable efficiency correction"};
  Configurable<float> cfgLowEffCut{"cfgLowEffCut", 0.001f, "Low efficiency cut"};
  Configurable<bool> isUseItsPid{"isUseItsPid", false, "Use ITS PID for particle identification"};
  Configurable<float> cfgTofPidPtCut{"cfgTofPidPtCut", 0.3f, "Minimum pt to use TOF N-sigma"};
  Configurable<int> cfgTrackPid{"cfgTrackPid", 0, "1 = pion, 2 = kaon, 3 = proton, 0 for no PID"};
  Configurable<std::vector<double>> itsNsigmaPidCut{"itsNsigmaPidCut", std::vector<double>{3, 2.5, 2, -3, -2.5, -2}, "ITS n-sigma cut for pions_posNsigma, kaons_posNsigma, protons_posNsigma, pions_negNsigma, kaons_negNsigma, protons_negNsigma"};
  Configurable<std::vector<double>> tpcNsigmaPidCut{"tpcNsigmaPidCut", std::vector<double>{1.5, 1.5, 1.5, -1.5, -1.5, -1.5}, "TPC n-sigma cut for pions_posNsigma, kaons_posNsigma, protons_posNsigma, pions_negNsigma, kaons_negNsigma, protons_negNsigma"};
  Configurable<std::vector<double>> tofNsigmaPidCut{"tofNsigmaPidCut", std::vector<double>{1.5, 1.5, 1.5, -1.5, -1.5, -1.5}, "TOF n-sigma cut for pions_posNsigma, kaons_posNsigma, protons_posNsigma, pions_negNsigma, kaons_negNsigma, protons_negNsigma"};
  Configurable<std::string> cfgEffccdbPath{"cfgEffccdbPath", "/alice/data/CCDB/Users/a/abmodak/OO/Efficiency", "Browse track eff object from CCDB"};
  Configurable<std::string> cfgMultccdbPath{"cfgMultccdbPath", "/alice/data/CCDB/Users/a/abmodak/OO/Multiplicity", "Browse mult efficiency object from CCDB"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {72, -PIHalf, PIHalf * 3}, "delta phi axis for histograms"};
  ConfigurableAxis axisDeltaEta{"axisDeltaEta", {40, -6, -2}, "delta eta axis for histograms"};
  ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt trigger axis for histograms"};
  ConfigurableAxis axisPtAssoc{"axisPtAssoc", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt associated axis for histograms"};
  ConfigurableAxis axisMultME{"axisMultME", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 1000}, "Mixing bins - multiplicity"};
  ConfigurableAxis axisVtxZME{"axisVtxZME", {VARIABLE_WIDTH, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10}, "Mixing bins - z-vertex"};
  ConfigurableAxis axisVertexEfficiency{"axisVertexEfficiency", {10, -10, 10}, "vertex axis for efficiency histograms"};
  ConfigurableAxis axisEtaEfficiency{"axisEtaEfficiency", {20, -1.0, 1.0}, "eta axis for efficiency histograms"};
  ConfigurableAxis axisPtEfficiency{"axisPtEfficiency", {VARIABLE_WIDTH, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0}, "pt axis for efficiency histograms"};
  ConfigurableAxis axisVtxZ{"axisVtxZ", {40, -20, 20}, "vertex axis"};
  ConfigurableAxis axisPhi{"axisPhi", {96, 0, TwoPI}, "#phi axis"};
  ConfigurableAxis axisEtaTrig{"axisEtaTrig", {40, -1., 1.}, "#eta trig axis"};
  ConfigurableAxis axisEtaAssoc{"axisEtaAssoc", {96, 3.5, 4.9}, "#eta assoc axis"};
  ConfigurableAxis axisSample{"axisSample", {cfgSampleSize, 0, cfgSampleSize}, "sample axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 10, 15, 25, 50, 60, 1000}, "multiplicity / centrality axis for histograms"};
  ConfigurableAxis amplitudeFt0a{"amplitudeFt0a", {5000, 0, 10000}, "FT0A amplitude"};
  ConfigurableAxis channelFt0aAxis{"channelFt0aAxis", {96, 0.0, 96.0}, "FT0A channel"};

  Configurable<bool> isApplyCentFT0C{"isApplyCentFT0C", true, "Centrality based on FT0C"};
  Configurable<bool> isApplyCentFV0A{"isApplyCentFV0A", false, "Centrality based on FV0A"};
  Configurable<bool> isApplyCentFT0M{"isApplyCentFT0M", false, "Centrality based on FT0A + FT0C"};
  Configurable<bool> isUseCentEst{"isUseCentEst", false, "Centrality based classification"};
  Configurable<float> cfgPtCutMult{"cfgPtCutMult", 3.0f, "maximum track pT for multiplicity classification"};

  using CollTable = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::CentFV0As, aod::CentFT0Ms>;
  using TrksTable = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFbeta, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>>;
  using MftTrkTable = soa::Filtered<aod::MFTTracks>;
  using CollTableMC = soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::CentFV0As, aod::CentFT0Ms>>;
  using TrksTableMC = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels, aod::TrackSelection>>;
  Preslice<TrksTable> perColGlobal = aod::track::collisionId;
  Preslice<TrksTableMC> perColMC = aod::track::collisionId;
  Preslice<MftTrkTable> perColMft = aod::fwdtrack::collisionId;

  OutputObj<CorrelationContainer> sameFt0aGlobal{"sameEventFt0aGlobal"};
  OutputObj<CorrelationContainer> mixedFt0aGlobal{"mixedEventFt0aGlobal"};
  OutputObj<CorrelationContainer> sameFt0cGlobal{"sameEventFt0cGlobal"};
  OutputObj<CorrelationContainer> mixedFt0cGlobal{"mixedEventFt0cGlobal"};
  OutputObj<CorrelationContainer> sameMftGlobal{"sameEventMftGlobal"};
  OutputObj<CorrelationContainer> mixedMftGlobal{"mixedEventMftGlobal"};
  OutputObj<CorrelationContainer> sameFt0aMft{"sameEventFt0aMft"};
  OutputObj<CorrelationContainer> mixedFt0aMft{"mixedEventFt0aMft"};
  OutputObj<CorrelationContainer> sameFt0aFt0c{"sameEventFt0aFt0c"};
  OutputObj<CorrelationContainer> mixedFt0aFt0c{"mixedEventFt0aFt0c"};

  // corrections
  TH3D* hTrkEff = nullptr;
  TH1D* hMultEff = nullptr;
  bool fLoadTrkEffCorr = false;
  bool fLoadMultEffCorr = false;

  std::vector<double> tofNsigmaCut;
  std::vector<double> itsNsigmaCut;
  std::vector<double> tpcNsigmaCut;
  o2::aod::ITSResponse itsResponse;

  template <KindOfCorrType corrType, KindOfEvntType evntType>
  void addhistosQA()
  {
    histos.add(Form("%s%shMult", kCorrType[corrType].data(), kEvntType[evntType].data()), "", kTH1D, {axisMultiplicity});
    histos.add(Form("%s%sTrig_etavsphi", kCorrType[corrType].data(), kEvntType[evntType].data()), "", kTH2D, {axisPhi, axisEtaTrig});
    histos.add(Form("%s%sTrig_eta", kCorrType[corrType].data(), kEvntType[evntType].data()), "", kTH1D, {axisEtaTrig});
    histos.add(Form("%s%sTrig_phi", kCorrType[corrType].data(), kEvntType[evntType].data()), "", kTH1D, {axisPhi});
    histos.add(Form("%s%sTrig_pt", kCorrType[corrType].data(), kEvntType[evntType].data()), "", kTH1D, {axisPtTrigger});
    histos.add(Form("%s%sTrig_hist", kCorrType[corrType].data(), kEvntType[evntType].data()), "", kTHnSparseF, {axisSample, axisVtxZ, axisPtTrigger, axisMultiplicity});
    histos.add(Form("%s%sAssoc_amp", kCorrType[corrType].data(), kEvntType[evntType].data()), "", kTH2D, {channelFt0aAxis, amplitudeFt0a});
    histos.add(Form("%s%sAssoc_eta", kCorrType[corrType].data(), kEvntType[evntType].data()), "", kTH1D, {axisEtaAssoc});
    histos.add(Form("%s%sAssoc_phi", kCorrType[corrType].data(), kEvntType[evntType].data()), "", kTH1D, {axisPhi});
    histos.add(Form("%s%sAssoc_etavsphi", kCorrType[corrType].data(), kEvntType[evntType].data()), "", kTH2D, {axisPhi, axisEtaAssoc});
    histos.add(Form("%s%sdeltaEta_deltaPhi", kCorrType[corrType].data(), kEvntType[evntType].data()), "", kTH2D, {axisDeltaPhi, axisDeltaEta});
  }

  void init(InitContext const&)
  {
    ccdb->setURL(cfgCcdbParam.cfgURL);
    ccdbApi.init("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    LOGF(info, "Getting alignment offsets from the CCDB...");
    offsetFT0 = ccdb->getForTimeStamp<std::vector<o2::detectors::AlignParam>>("FT0/Calib/Align", cfgCcdbParam.noLaterThan.value);
    LOGF(info, "Offset for FT0A: x = %.3f y = %.3f z = %.3f\n", (*offsetFT0)[0].getX(), (*offsetFT0)[0].getY(), (*offsetFT0)[0].getZ());
    LOGF(info, "Offset for FT0C: x = %.3f y = %.3f z = %.3f\n", (*offsetFT0)[1].getX(), (*offsetFT0)[1].getY(), (*offsetFT0)[1].getZ());

    std::vector<AxisSpec> corrAxis = {{axisSample, "Sample"},
                                      {axisVtxZ, "z-vtx (cm)"},
                                      {axisPtTrigger, "p_{T} (GeV/c)"},
                                      {axisPtAssoc, "p_{T} (GeV/c)"},
                                      {axisDeltaPhi, "#Delta#varphi (rad)"},
                                      {axisDeltaEta, "#Delta#eta"}};
    std::vector<AxisSpec> effAxis = {{axisVertexEfficiency, "z-vtx (cm)"},
                                     {axisPtEfficiency, "p_{T} (GeV/c)"},
                                     {axisEtaEfficiency, "#eta"}};

    std::vector<AxisSpec> userAxis;
    userAxis.emplace_back(axisMultiplicity, "multiplicity");

    tofNsigmaCut = tofNsigmaPidCut;
    itsNsigmaCut = itsNsigmaPidCut;
    tpcNsigmaCut = tpcNsigmaPidCut;

    if (doprocessEventStat) {
      histos.add("QA/EventHist", "events", kTH1F, {axisEvent}, false);
      histos.add("QA/VtxZHist", "v_{z} (cm)", kTH1F, {axisVtxZ}, false);

      auto hstat = histos.get<TH1>(HIST("QA/EventHist"));
      auto* x = hstat->GetXaxis();
      x->SetBinLabel(1, "All events");
      x->SetBinLabel(2, "sel8");
      x->SetBinLabel(3, "kNoSameBunchPileup"); // reject collisions in case of pileup with another collision in the same foundBC
      x->SetBinLabel(4, "kIsGoodZvtxFT0vsPV"); // small difference between z-vertex from PV and from FT0
      x->SetBinLabel(5, "|vz|<10");
    }

    if (doprocessFt0aGlobalSE || doprocessFt0aGlobalME) {
      addhistosQA<kFT0AGLOBAL, kSE>();
      histos.add("Ft0aGlobal/ME/deltaEta_deltaPhi", "", kTH2D, {axisDeltaPhi, axisDeltaEta});
      sameFt0aGlobal.setObject(new CorrelationContainer("sameEventFt0aGlobal", "sameEventFt0aGlobal", corrAxis, effAxis, userAxis));
      mixedFt0aGlobal.setObject(new CorrelationContainer("mixedEventFt0aGlobal", "mixedEventFt0aGlobal", corrAxis, effAxis, userAxis));
    }

    if (doprocessFt0cGlobalSE || doprocessFt0cGlobalME) {
      addhistosQA<kFT0CGLOBAL, kSE>();
      histos.add("Ft0cGlobal/ME/deltaEta_deltaPhi", "", kTH2D, {axisDeltaPhi, axisDeltaEta});
      sameFt0cGlobal.setObject(new CorrelationContainer("sameEventFt0cGlobal", "sameEventFt0cGlobal", corrAxis, effAxis, userAxis));
      mixedFt0cGlobal.setObject(new CorrelationContainer("mixedEventFt0cGlobal", "mixedEventFt0cGlobal", corrAxis, effAxis, userAxis));
    }

    if (doprocessMftGlobalSE || doprocessMftGlobalME) {
      addhistosQA<kMFTGLOBAL, kSE>();
      histos.add("MftGlobal/ME/deltaEta_deltaPhi", "", kTH2D, {axisDeltaPhi, axisDeltaEta});
      sameMftGlobal.setObject(new CorrelationContainer("sameEventMftGlobal", "sameEventMftGlobal", corrAxis, effAxis, userAxis));
      mixedMftGlobal.setObject(new CorrelationContainer("mixedEventMftGlobal", "mixedEventMftGlobal", corrAxis, effAxis, userAxis));
    }

    if (doprocessFt0aMftSE || doprocessFt0aMftME) {
      addhistosQA<kFT0AMFT, kSE>();
      histos.add("Ft0aMft/ME/deltaEta_deltaPhi", "", kTH2D, {axisDeltaPhi, axisDeltaEta});
      sameFt0aMft.setObject(new CorrelationContainer("sameEventFt0aMft", "sameEventFt0aMft", corrAxis, effAxis, userAxis));
      mixedFt0aMft.setObject(new CorrelationContainer("mixedEventFt0aMft", "mixedEventFt0aMft", corrAxis, effAxis, userAxis));
    }

    if (doprocessFt0aFt0cSE || doprocessFt0aFt0cME) {
      addhistosQA<kFT0AFT0C, kSE>();
      histos.add("Ft0aFt0c/ME/deltaEta_deltaPhi", "", kTH2D, {axisDeltaPhi, axisDeltaEta});
      sameFt0aFt0c.setObject(new CorrelationContainer("sameEventFt0aFt0c", "sameEventFt0aFt0c", corrAxis, effAxis, userAxis));
      mixedFt0aFt0c.setObject(new CorrelationContainer("mixedEventFt0aFt0c", "mixedEventFt0aFt0c", corrAxis, effAxis, userAxis));
    }

    if (doprocessEff) {
      histos.add("hmcgendndptPrimary", "hmcgendndptPrimary", kTHnSparseD, {axisEtaTrig, axisPtTrigger, axisMultiplicity, axisVtxZ}, false);
      histos.add("hmcrecdndptRecoPrimary", "hmcrecdndptRecoPrimary", kTHnSparseD, {axisEtaTrig, axisPtTrigger, axisMultiplicity, axisVtxZ}, false);
      histos.add("hmcrecdndptRecoAll", "hmcrecdndptRecoAll", kTHnSparseD, {axisEtaTrig, axisPtTrigger, axisMultiplicity, axisVtxZ}, false);
      histos.add("hmcrecdndptFake", "hmcrecdndptFake", kTHnSparseD, {axisEtaTrig, axisPtTrigger, axisMultiplicity, axisVtxZ}, false);
    }
  }

  Filter fTrackSelectionITS = ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::ITS) &&
                              ncheckbit(aod::track::trackCutFlag, TrackSelectionIts);
  Filter fTrackSelectionTPC = ifnode(ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::TPC),
                                     ncheckbit(aod::track::trackCutFlag, TrackSelectionTpc), true);
  Filter fTrackSelectionDCA = ifnode(dcaZ.node() > 0.f, nabs(aod::track::dcaZ) <= dcaZ && ncheckbit(aod::track::trackCutFlag, TrackSelectionDcaxyOnly),
                                     ncheckbit(aod::track::trackCutFlag, TrackSelectionDca));
  Filter fTracksEta = nabs(aod::track::eta) < cfgEtaCut;
  Filter fTracksPt = (aod::track::pt > cfgPtCutMin) && (aod::track::pt < cfgPtCutMax);

  Filter fMftTrackEta = (aod::fwdtrack::eta < cfigMftEtaMax) && (aod::fwdtrack::eta > cfigMftEtaMin);
  Filter fMftTrackColID = (aod::fwdtrack::bestCollisionId >= 0);
  Filter fMftTrackDca = (nabs(aod::fwdtrack::bestDCAXY) < cfigMftDcaxy);

  double getPhiFT0(uint chno, int i)
  {
    ft0Det.calculateChannelCenter();
    auto chPos = ft0Det.getChannelCenter(chno);
    return RecoDecay::phi(chPos.X() + (*offsetFT0)[i].getX(), chPos.Y() + (*offsetFT0)[i].getY());
  }

  double getEtaFT0(uint chno, int i)
  {
    ft0Det.calculateChannelCenter();
    auto chPos = ft0Det.getChannelCenter(chno);
    auto x = chPos.X() + (*offsetFT0)[i].getX();
    auto y = chPos.Y() + (*offsetFT0)[i].getY();
    auto z = chPos.Z() + (*offsetFT0)[i].getZ();
    if (chno >= kMinFt0cCell)
      z = -z;
    auto r = std::sqrt(x * x + y * y);
    auto theta = std::atan2(r, z);
    return -std::log(std::tan(0.5 * theta));
  }

  template <typename CheckCol>
  bool isEventSelected(CheckCol const& col)
  {
    if (!col.sel8()) {
      return false;
    }
    if (isApplySameBunchPileup && !col.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if (isApplyGoodZvtxFT0vsPV && !col.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if (std::abs(col.posZ()) >= cfgVtxCut) {
      return false;
    }
    return true;
  }

  template <typename CheckMftTrack>
  bool isMftTrackSelected(CheckMftTrack const& track)
  {
    if (track.nClusters() < cfigMftCluster) {
      return false;
    }
    return true;
  }

  template <KindOfCorrType corrType, KindOfEvntType evntType, typename TTrack>
  void fillTrackQA(TTrack const& track)
  {
    histos.fill(HIST(kCorrType[corrType]) + HIST(kEvntType[evntType]) + HIST("Trig_etavsphi"), track.phi(), track.eta());
    histos.fill(HIST(kCorrType[corrType]) + HIST(kEvntType[evntType]) + HIST("Trig_eta"), track.eta());
    histos.fill(HIST(kCorrType[corrType]) + HIST(kEvntType[evntType]) + HIST("Trig_phi"), track.phi());
    histos.fill(HIST(kCorrType[corrType]) + HIST(kEvntType[evntType]) + HIST("Trig_pt"), track.pt());
  }

  template <typename countTrk>
  int countNTracks(countTrk const& tracks)
  {
    auto nTrk = 0;
    for (const auto& track : tracks) {
      if (track.pt() < cfgPtCutMin || track.pt() > cfgPtCutMult) {
        continue;
      }
      nTrk++;
    }
    return nTrk;
  }

  template <typename CheckColCent>
  float selColCent(CheckColCent const& col)
  {
    auto cent = -1;
    if (isApplyCentFT0C) {
      cent = col.centFT0C();
    }
    if (isApplyCentFV0A) {
      cent = col.centFV0A();
    }
    if (isApplyCentFT0M) {
      cent = col.centFT0M();
    }
    return cent;
  }

  void loadEffCorrection(uint64_t timestamp)
  {
    if (fLoadTrkEffCorr) {
      return;
    }
    if (cfgEffccdbPath.value.empty() == false) {
      hTrkEff = ccdb->getForTimeStamp<TH3D>(cfgEffccdbPath, timestamp);
      if (hTrkEff == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", cfgEffccdbPath.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram from %s (%p)", cfgEffccdbPath.value.c_str(), (void*)hTrkEff);
    }
    fLoadTrkEffCorr = true;
  }

  void loadMultCorrection(uint64_t timestamp)
  {
    if (fLoadMultEffCorr) {
      return;
    }
    if (cfgMultccdbPath.value.empty() == false) {
      hMultEff = ccdb->getForTimeStamp<TH1D>(cfgMultccdbPath, timestamp);
      if (hMultEff == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for multiplicity from %s", cfgMultccdbPath.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram from %s (%p)", cfgMultccdbPath.value.c_str(), (void*)hMultEff);
    }
    fLoadMultEffCorr = true;
  }

  float getTrkEffCorr(float eta, float pt, float posZ)
  {
    float eff = 1.;
    if (hTrkEff) {
      int etaBin = hTrkEff->GetXaxis()->FindBin(eta);
      int ptBin = hTrkEff->GetYaxis()->FindBin(pt);
      int zBin = hTrkEff->GetZaxis()->FindBin(posZ);
      eff = hTrkEff->GetBinContent(etaBin, ptBin, zBin);
    } else {
      eff = 1.0;
    }
    if (eff < cfgLowEffCut)
      eff = 1.0;

    return eff;
  }

  float getMultEffCorr(float mult)
  {
    float eff = 1.;
    if (hMultEff) {
      int multBin = hMultEff->GetXaxis()->FindBin(mult);
      eff = hMultEff->GetBinContent(multBin);
    } else {
      eff = 1.0;
    }
    if (eff < cfgLowEffCut)
      eff = 1.0;

    return eff;
  }

  template <typename TTrack>
  int getTrackPID(TTrack track)
  {
    // Computing Nsigma arrays for pion, kaon, and protons
    std::array<float, 3> nSigmaTPC = {track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr()};
    std::array<float, 3> nSigmaTOF = {track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr()};
    std::array<float, 3> nSigmaITS = {itsResponse.nSigmaITS<o2::track::PID::Pion>(track), itsResponse.nSigmaITS<o2::track::PID::Kaon>(track), itsResponse.nSigmaITS<o2::track::PID::Proton>(track)};
    int pid = -1;

    std::array<float, 3> nSigmaToUse = isUseItsPid ? nSigmaITS : nSigmaTPC;            // Choose which nSigma to use: TPC or ITS
    std::vector<double> detectorNsigmaCut = isUseItsPid ? itsNsigmaCut : tpcNsigmaCut; // Choose which nSigma to use: TPC or ITS

    bool isPion, isKaon, isProton;
    bool isDetectedPion = nSigmaToUse[0] < detectorNsigmaCut[0] && nSigmaToUse[0] > detectorNsigmaCut[0 + 3];
    bool isDetectedKaon = nSigmaToUse[1] < detectorNsigmaCut[1] && nSigmaToUse[1] > detectorNsigmaCut[1 + 3];
    bool isDetectedProton = nSigmaToUse[2] < detectorNsigmaCut[2] && nSigmaToUse[2] > detectorNsigmaCut[2 + 3];

    bool isTofPion = nSigmaTOF[0] < tofNsigmaCut[0] && nSigmaTOF[0] > tofNsigmaCut[0 + 3];
    bool isTofKaon = nSigmaTOF[1] < tofNsigmaCut[1] && nSigmaTOF[1] > tofNsigmaCut[1 + 3];
    bool isTofProton = nSigmaTOF[2] < tofNsigmaCut[2] && nSigmaTOF[2] > tofNsigmaCut[2 + 3];

    if (track.pt() > cfgTofPidPtCut && !track.hasTOF()) {
      return 0;
    } else if (track.pt() > cfgTofPidPtCut && track.hasTOF()) {
      isPion = isTofPion && isDetectedPion;
      isKaon = isTofKaon && isDetectedKaon;
      isProton = isTofProton && isDetectedProton;
    } else {
      isPion = isDetectedPion;
      isKaon = isDetectedKaon;
      isProton = isDetectedProton;
    }

    if ((isPion && isKaon) || (isPion && isProton) || (isKaon && isProton)) {
      return 0; // more than one particle satisfy the criteria
    }

    if (isPion) {
      pid = PIONS;
    } else if (isKaon) {
      pid = KAONS;
    } else if (isProton) {
      pid = PROTONS;
    } else {
      return 0; // no particle satisfies the criteria
    }

    return pid + 1; // shift the pid by 1, 1 = pion, 2 = kaon, 3 = proton
  }

  template <CorrelationContainer::CFStep step, typename TTarget, typename TTriggers, typename TFT0s>
  void fillCorrFt0aGlobal(TTarget target, TTriggers const& triggers, TFT0s const& ft0, bool mixing, float vz, float multiplicity)
  {
    histos.fill(HIST("Ft0aGlobal/SE/hMult"), multiplicity);
    int fSampleIndex = gRandom->Uniform(0, cfgSampleSize);
    for (auto const& triggerTrack : triggers) {
      if (cfgTrackPid && getTrackPID(triggerTrack) != cfgTrackPid)
        continue; // if PID is selected, check if the track has the right PID

      float trkeffw = 1.0f;
      if (isUseEffCorr)
        trkeffw = getTrkEffCorr(triggerTrack.eta(), triggerTrack.pt(), vz);

      if (!mixing) {
        fillTrackQA<kFT0AGLOBAL, kSE>(triggerTrack);
        histos.fill(HIST("Ft0aGlobal/SE/Trig_hist"), fSampleIndex, vz, triggerTrack.pt(), multiplicity, trkeffw);
      }

      for (std::size_t iCh = 0; iCh < ft0.channelA().size(); iCh++) {
        auto chanelid = ft0.channelA()[iCh];
        float ampl = ft0.amplitudeA()[iCh];

        auto phi = getPhiFT0(chanelid, 0);
        auto eta = getEtaFT0(chanelid, 0);

        if (!mixing) {
          histos.fill(HIST("Ft0aGlobal/SE/Assoc_amp"), chanelid, ampl);
          histos.fill(HIST("Ft0aGlobal/SE/Assoc_eta"), eta);
          histos.fill(HIST("Ft0aGlobal/SE/Assoc_phi"), phi);
          histos.fill(HIST("Ft0aGlobal/SE/Assoc_etavsphi"), phi, eta);
        }
        float deltaPhi = RecoDecay::constrainAngle(triggerTrack.phi() - phi, -PIHalf);
        float deltaEta = triggerTrack.eta() - eta;
        if (mixing)
          histos.fill(HIST("Ft0aGlobal/ME/deltaEta_deltaPhi"), deltaPhi, deltaEta, trkeffw);
        else
          histos.fill(HIST("Ft0aGlobal/SE/deltaEta_deltaPhi"), deltaPhi, deltaEta, trkeffw);
        target->getPairHist()->Fill(step, fSampleIndex, vz, triggerTrack.pt(), triggerTrack.pt(), deltaPhi, deltaEta, multiplicity, trkeffw);
      } // associated ft0 tracks
    } // trigger tracks
  } // fillCorrFt0aGlobal

  template <CorrelationContainer::CFStep step, typename TTarget, typename TTriggers, typename TFT0s>
  void fillCorrFt0cGlobal(TTarget target, TTriggers const& triggers, TFT0s const& ft0, bool mixing, float vz, float multiplicity)
  {
    histos.fill(HIST("Ft0cGlobal/SE/hMult"), multiplicity);
    int fSampleIndex = gRandom->Uniform(0, cfgSampleSize);
    for (auto const& triggerTrack : triggers) {
      if (cfgTrackPid && getTrackPID(triggerTrack) != cfgTrackPid)
        continue; // if PID is selected, check if the track has the right PID

      float trkeffw = 1.0f;
      if (isUseEffCorr)
        trkeffw = getTrkEffCorr(triggerTrack.eta(), triggerTrack.pt(), vz);

      if (!mixing) {
        fillTrackQA<kFT0CGLOBAL, kSE>(triggerTrack);
        histos.fill(HIST("Ft0cGlobal/SE/Trig_hist"), fSampleIndex, vz, triggerTrack.pt(), multiplicity, trkeffw);
      }

      for (std::size_t iCh = 0; iCh < ft0.channelC().size(); iCh++) {
        auto chanelid = ft0.channelC()[iCh] + 96;
        float ampl = ft0.amplitudeC()[iCh];

        auto phi = getPhiFT0(chanelid, 1);
        auto eta = getEtaFT0(chanelid, 1);

        if (!mixing) {
          histos.fill(HIST("Ft0cGlobal/SE/Assoc_amp"), chanelid, ampl);
          histos.fill(HIST("Ft0cGlobal/SE/Assoc_eta"), eta);
          histos.fill(HIST("Ft0cGlobal/SE/Assoc_phi"), phi);
          histos.fill(HIST("Ft0cGlobal/SE/Assoc_etavsphi"), phi, eta);
        }

        float deltaPhi = RecoDecay::constrainAngle(triggerTrack.phi() - phi, -PIHalf);
        float deltaEta = triggerTrack.eta() - eta;
        if (mixing)
          histos.fill(HIST("Ft0cGlobal/ME/deltaEta_deltaPhi"), deltaPhi, deltaEta, trkeffw);
        else
          histos.fill(HIST("Ft0cGlobal/SE/deltaEta_deltaPhi"), deltaPhi, deltaEta, trkeffw);
        target->getPairHist()->Fill(step, fSampleIndex, vz, triggerTrack.pt(), triggerTrack.pt(), deltaPhi, deltaEta, multiplicity, trkeffw);
      } // associated ft0 tracks
    } // trigger tracks
  } // fillCorrFt0cGlobal

  template <CorrelationContainer::CFStep step, typename TTarget, typename TTriggers, typename TMFTs>
  void fillCorrMftGlobal(TTarget target, TTriggers const& triggers, TMFTs const& mft, bool mixing, float vz, float multiplicity)
  {
    histos.fill(HIST("MftGlobal/SE/hMult"), multiplicity);
    int fSampleIndex = gRandom->Uniform(0, cfgSampleSize);
    for (auto const& triggerTrack : triggers) {
      if (cfgTrackPid && getTrackPID(triggerTrack) != cfgTrackPid)
        continue; // if PID is selected, check if the track has the right PID

      float trkeffw = 1.0f;
      if (isUseEffCorr)
        trkeffw = getTrkEffCorr(triggerTrack.eta(), triggerTrack.pt(), vz);

      if (!mixing) {
        fillTrackQA<kMFTGLOBAL, kSE>(triggerTrack);
        histos.fill(HIST("MftGlobal/SE/Trig_hist"), fSampleIndex, vz, triggerTrack.pt(), multiplicity, trkeffw);
      }

      for (auto const& assoTrack : mft) {
        if (!isMftTrackSelected(assoTrack)) {
          continue;
        }
        auto phi = assoTrack.phi();
        o2::math_utils::bringTo02Pi(phi);

        if (!mixing) {
          histos.fill(HIST("MftGlobal/SE/Assoc_eta"), assoTrack.eta());
          histos.fill(HIST("MftGlobal/SE/Assoc_phi"), phi);
          histos.fill(HIST("MftGlobal/SE/Assoc_etavsphi"), phi, assoTrack.eta());
        }

        float deltaPhi = RecoDecay::constrainAngle(triggerTrack.phi() - phi, -PIHalf);
        float deltaEta = triggerTrack.eta() - assoTrack.eta();
        if (mixing)
          histos.fill(HIST("MftGlobal/ME/deltaEta_deltaPhi"), deltaPhi, deltaEta, trkeffw);
        else
          histos.fill(HIST("MftGlobal/SE/deltaEta_deltaPhi"), deltaPhi, deltaEta, trkeffw);
        target->getPairHist()->Fill(step, fSampleIndex, vz, triggerTrack.pt(), assoTrack.pt(), deltaPhi, deltaEta, multiplicity, trkeffw);
      } // associated mft tracks
    } // trigger tracks
  } // fillCorrMftGlobal

  template <CorrelationContainer::CFStep step, typename TTarget, typename TTriggers, typename TFT0s>
  void fillCorrFt0aMft(TTarget target, TTriggers const& triggers, TFT0s const& ft0, bool mixing, float vz, float multiplicity)
  {
    histos.fill(HIST("Ft0aMft/SE/hMult"), multiplicity);
    int fSampleIndex = gRandom->Uniform(0, cfgSampleSize);
    for (auto const& triggerTrack : triggers) {
      if (!isMftTrackSelected(triggerTrack)) {
        continue;
      }

      auto trigphi = triggerTrack.phi();
      o2::math_utils::bringTo02Pi(trigphi);

      if (!mixing) {
        fillTrackQA<kFT0AMFT, kSE>(triggerTrack);
        histos.fill(HIST("Ft0aMft/SE/Trig_hist"), fSampleIndex, vz, triggerTrack.pt(), multiplicity);
      }

      for (std::size_t iCh = 0; iCh < ft0.channelA().size(); iCh++) {
        auto chanelid = ft0.channelA()[iCh];
        float ampl = ft0.amplitudeA()[iCh];

        auto phi = getPhiFT0(chanelid, 0);
        auto eta = getEtaFT0(chanelid, 0);

        if (!mixing) {
          histos.fill(HIST("Ft0aMft/SE/Assoc_amp"), chanelid, ampl);
          histos.fill(HIST("Ft0aMft/SE/Assoc_eta"), eta);
          histos.fill(HIST("Ft0aMft/SE/Assoc_phi"), phi);
          histos.fill(HIST("Ft0aMft/SE/Assoc_etavsphi"), phi, eta);
        }

        float deltaPhi = RecoDecay::constrainAngle(trigphi - phi, -PIHalf);
        float deltaEta = triggerTrack.eta() - eta;
        if (mixing)
          histos.fill(HIST("Ft0aMft/ME/deltaEta_deltaPhi"), deltaPhi, deltaEta);
        else
          histos.fill(HIST("Ft0aMft/SE/deltaEta_deltaPhi"), deltaPhi, deltaEta);
        target->getPairHist()->Fill(step, fSampleIndex, vz, triggerTrack.pt(), triggerTrack.pt(), deltaPhi, deltaEta, multiplicity);
      } // associated ft0 tracks
    } // trigger tracks
  } // fillCorrFt0aMft

  template <CorrelationContainer::CFStep step, typename TTarget, typename TFT0As, typename TFT0Cs>
  void fillCorrFt0aFt0c(TTarget target, TFT0As const& ft0a, TFT0Cs const& ft0c, bool mixing, float vz, float multiplicity)
  {
    histos.fill(HIST("Ft0aFt0c/SE/hMult"), multiplicity);
    int fSampleIndex = gRandom->Uniform(0, cfgSampleSize);
    for (std::size_t iChA = 0; iChA < ft0a.channelA().size(); iChA++) {

      auto chanelidA = ft0a.channelA()[iChA];
      auto phiA = getPhiFT0(chanelidA, 0);
      auto etaA = getEtaFT0(chanelidA, 0);

      if (!mixing) {
        histos.fill(HIST("Ft0aFt0c/SE/Trig_eta"), etaA);
        histos.fill(HIST("Ft0aFt0c/SE/Trig_phi"), phiA);
        histos.fill(HIST("Ft0aFt0c/SE/Trig_etavsphi"), phiA, etaA);
        histos.fill(HIST("Ft0aFt0c/SE/Trig_hist"), fSampleIndex, vz, 1.0, multiplicity);
      }

      for (std::size_t iChC = 0; iChC < ft0c.channelC().size(); iChC++) {

        auto chanelidC = ft0c.channelC()[iChC] + 96;
        float ampl = ft0c.amplitudeC()[iChC];
        auto phiC = getPhiFT0(chanelidC, 1);
        auto etaC = getEtaFT0(chanelidC, 1);

        if (mixing) {
          histos.fill(HIST("Ft0aFt0c/SE/Assoc_amp"), chanelidC, ampl);
          histos.fill(HIST("Ft0aFt0c/SE/Assoc_eta"), etaC);
          histos.fill(HIST("Ft0aFt0c/SE/Assoc_phi"), phiC);
          histos.fill(HIST("Ft0aFt0c/SE/Assoc_etavsphi"), phiC, etaC);
        }

        float deltaPhi = RecoDecay::constrainAngle(phiA - phiC, -PIHalf);
        float deltaEta = etaA - etaC;
        if (mixing)
          histos.fill(HIST("Ft0aFt0c/ME/deltaEta_deltaPhi"), deltaPhi, deltaEta);
        else
          histos.fill(HIST("Ft0aFt0c/SE/deltaEta_deltaPhi"), deltaPhi, deltaEta);
        target->getPairHist()->Fill(step, fSampleIndex, vz, 1.0, 1.0, deltaPhi, deltaEta, multiplicity);
      } // associated ft0 tracks
    } // trigger tracks
  } // fillCorrFt0aFt0c

  void processEventStat(CollTable::iterator const& col)
  {
    histos.fill(HIST("QA/EventHist"), 1);
    if (!col.sel8()) {
      return;
    }
    histos.fill(HIST("QA/EventHist"), 2);
    if (isApplySameBunchPileup && !col.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return;
    }
    histos.fill(HIST("QA/EventHist"), 3);
    if (isApplyGoodZvtxFT0vsPV && !col.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return;
    }
    histos.fill(HIST("QA/EventHist"), 4);
    if (std::abs(col.posZ()) >= cfgVtxCut) {
      return;
    }
    histos.fill(HIST("QA/EventHist"), 5);
    histos.fill(HIST("QA/VtxZHist"), col.posZ());
  }

  void processFt0aGlobalSE(CollTable::iterator const& col, aod::FT0s const&, TrksTable const& tracks, aod::BCsWithTimestamps const&)
  {
    if (!isEventSelected(col)) {
      return;
    }
    if (col.has_foundFT0()) {
      auto bc = col.bc_as<aod::BCsWithTimestamps>();
      loadEffCorrection(bc.timestamp());
      loadMultCorrection(bc.timestamp());
      const auto& ft0 = col.foundFT0();
      auto multiplicity = 1.0f;
      if (isUseCentEst)
        multiplicity = selColCent(col);
      else
        multiplicity = countNTracks(tracks);
      float multw = getMultEffCorr(multiplicity);
      if (isUseEffCorr)
        multiplicity = multiplicity * multw;
      fillCorrFt0aGlobal<CorrelationContainer::kCFStepReconstructed>(sameFt0aGlobal, tracks, ft0, false, col.posZ(), multiplicity);
    }
  } // same event

  void processFt0cGlobalSE(CollTable::iterator const& col, aod::FT0s const&, TrksTable const& tracks, aod::BCsWithTimestamps const&)
  {
    if (!isEventSelected(col)) {
      return;
    }
    if (col.has_foundFT0()) {
      auto bc = col.bc_as<aod::BCsWithTimestamps>();
      loadEffCorrection(bc.timestamp());
      loadMultCorrection(bc.timestamp());
      const auto& ft0 = col.foundFT0();
      auto multiplicity = 1.0f;
      if (isUseCentEst)
        multiplicity = selColCent(col);
      else
        multiplicity = countNTracks(tracks);
      float multw = getMultEffCorr(multiplicity);
      if (isUseEffCorr)
        multiplicity = multiplicity * multw;
      fillCorrFt0cGlobal<CorrelationContainer::kCFStepReconstructed>(sameFt0cGlobal, tracks, ft0, false, col.posZ(), multiplicity);
    }
  } // same event

  void processMftGlobalSE(CollTable::iterator const& col, MftTrkTable const& mfttracks, TrksTable const& tracks, aod::BCsWithTimestamps const&)
  {
    if (!isEventSelected(col)) {
      return;
    }
    auto bc = col.bc_as<aod::BCsWithTimestamps>();
    loadEffCorrection(bc.timestamp());
    loadMultCorrection(bc.timestamp());
    auto multiplicity = 1.0f;
    if (isUseCentEst)
      multiplicity = selColCent(col);
    else
      multiplicity = countNTracks(tracks);
    float multw = getMultEffCorr(multiplicity);
    if (isUseEffCorr)
      multiplicity = multiplicity * multw;
    fillCorrMftGlobal<CorrelationContainer::kCFStepReconstructed>(sameMftGlobal, tracks, mfttracks, false, col.posZ(), multiplicity);
  } // same event

  void processFt0aMftSE(CollTable::iterator const& col, aod::FT0s const&, MftTrkTable const& mfttracks, TrksTable const& tracks, aod::BCsWithTimestamps const&)
  {
    if (!isEventSelected(col)) {
      return;
    }
    if (col.has_foundFT0()) {
      auto bc = col.bc_as<aod::BCsWithTimestamps>();
      loadMultCorrection(bc.timestamp());
      const auto& ft0 = col.foundFT0();
      auto multiplicity = 1.0f;
      if (isUseCentEst)
        multiplicity = selColCent(col);
      else
        multiplicity = countNTracks(tracks);
      float multw = getMultEffCorr(multiplicity);
      if (isUseEffCorr)
        multiplicity = multiplicity * multw;
      fillCorrFt0aMft<CorrelationContainer::kCFStepReconstructed>(sameFt0aMft, mfttracks, ft0, false, col.posZ(), multiplicity);
    }
  } // same event

  void processFt0aFt0cSE(CollTable::iterator const& col, aod::FT0s const&, TrksTable const& tracks, aod::BCsWithTimestamps const&)
  {
    if (!isEventSelected(col)) {
      return;
    }
    if (col.has_foundFT0()) {
      auto bc = col.bc_as<aod::BCsWithTimestamps>();
      const auto& ft0 = col.foundFT0();
      loadMultCorrection(bc.timestamp());
      auto multiplicity = 1.0f;
      if (isUseCentEst)
        multiplicity = selColCent(col);
      else
        multiplicity = countNTracks(tracks);
      float multw = getMultEffCorr(multiplicity);
      if (isUseEffCorr)
        multiplicity = multiplicity * multw;
      fillCorrFt0aFt0c<CorrelationContainer::kCFStepReconstructed>(sameFt0aFt0c, ft0, ft0, false, col.posZ(), multiplicity);
    }
  } // same event

  void processFt0aGlobalME(CollTable const& col, aod::FT0s const&, TrksTable const& tracks, aod::BCsWithTimestamps const&)
  {
    auto getTracksSize = [&tracks, this](CollTable::iterator const& collision) {
      auto associatedTracks = tracks.sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), this->cache);
      return associatedTracks.size();
    };

    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, aod::collision::PosZ, decltype(getTracksSize)>;
    MixedBinning binningOnVtxAndMult{{getTracksSize}, {axisVtxZME, axisMultME}, true};
    for (auto const& [col1, col2] : soa::selfCombinations(binningOnVtxAndMult, mixingParameter, -1, col, col)) {
      if (!isEventSelected(col1) || !isEventSelected(col2)) {
        continue;
      }
      if (col1.globalIndex() == col2.globalIndex()) {
        continue;
      }
      if (col1.has_foundFT0() && col2.has_foundFT0()) {
        auto bc = col1.bc_as<aod::BCsWithTimestamps>();
        loadEffCorrection(bc.timestamp());
        loadMultCorrection(bc.timestamp());
        auto slicedTriggerTracks = tracks.sliceBy(perColGlobal, col1.globalIndex());
        const auto& ft0 = col2.foundFT0();
        auto multiplicity = 1.0f;
        if (isUseCentEst)
          multiplicity = selColCent(col1);
        else
          multiplicity = countNTracks(slicedTriggerTracks);
        float multw = getMultEffCorr(multiplicity);
        if (isUseEffCorr)
          multiplicity = multiplicity * multw;
        fillCorrFt0aGlobal<CorrelationContainer::kCFStepReconstructed>(mixedFt0aGlobal, slicedTriggerTracks, ft0, true, col1.posZ(), multiplicity);
      }
    }
  } // mixed event

  void processFt0cGlobalME(CollTable const& col, aod::FT0s const&, TrksTable const& tracks, aod::BCsWithTimestamps const&)
  {
    auto getTracksSize = [&tracks, this](CollTable::iterator const& collision) {
      auto associatedTracks = tracks.sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), this->cache);
      return associatedTracks.size();
    };

    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, aod::collision::PosZ, decltype(getTracksSize)>;
    MixedBinning binningOnVtxAndMult{{getTracksSize}, {axisVtxZME, axisMultME}, true};
    for (auto const& [col1, col2] : soa::selfCombinations(binningOnVtxAndMult, mixingParameter, -1, col, col)) {
      if (!isEventSelected(col1) || !isEventSelected(col2)) {
        continue;
      }
      if (col1.globalIndex() == col2.globalIndex()) {
        continue;
      }
      if (col1.has_foundFT0() && col2.has_foundFT0()) {
        auto bc = col1.bc_as<aod::BCsWithTimestamps>();
        loadEffCorrection(bc.timestamp());
        loadMultCorrection(bc.timestamp());
        auto slicedTriggerTracks = tracks.sliceBy(perColGlobal, col1.globalIndex());
        const auto& ft0 = col2.foundFT0();
        auto multiplicity = 1.0f;
        if (isUseCentEst)
          multiplicity = selColCent(col1);
        else
          multiplicity = countNTracks(slicedTriggerTracks);
        float multw = getMultEffCorr(multiplicity);
        if (isUseEffCorr)
          multiplicity = multiplicity * multw;
        fillCorrFt0cGlobal<CorrelationContainer::kCFStepReconstructed>(mixedFt0cGlobal, slicedTriggerTracks, ft0, true, col1.posZ(), multiplicity);
      }
    }
  } // mixed event

  void processMftGlobalME(CollTable const& col, MftTrkTable const& mfttracks, TrksTable const& tracks, aod::BCsWithTimestamps const&)
  {
    auto getTracksSize = [&tracks, this](CollTable::iterator const& collision) {
      auto associatedTracks = tracks.sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), this->cache);
      return associatedTracks.size();
    };

    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, aod::collision::PosZ, decltype(getTracksSize)>;
    MixedBinning binningOnVtxAndMult{{getTracksSize}, {axisVtxZME, axisMultME}, true};
    auto tracksTuple = std::make_tuple(tracks, mfttracks);
    Pair<CollTable, TrksTable, MftTrkTable, MixedBinning> pairs{binningOnVtxAndMult, mixingParameter, -1, col, tracksTuple, &cache};
    for (auto const& [col1, tracks1, col2, tracks2] : pairs) {
      if (!isEventSelected(col1) || !isEventSelected(col2)) {
        continue;
      }
      auto bc = col1.bc_as<aod::BCsWithTimestamps>();
      loadEffCorrection(bc.timestamp());
      loadMultCorrection(bc.timestamp());
      auto multiplicity = 1.0f;
      if (isUseCentEst)
        multiplicity = selColCent(col1);
      else
        multiplicity = countNTracks(tracks1);
      float multw = getMultEffCorr(multiplicity);
      if (isUseEffCorr)
        multiplicity = multiplicity * multw;
      fillCorrMftGlobal<CorrelationContainer::kCFStepReconstructed>(mixedMftGlobal, tracks1, tracks2, true, col1.posZ(), multiplicity);
    }
  } // mixed event

  void processFt0aMftME(CollTable const& col, aod::FT0s const&, TrksTable const& tracks, MftTrkTable const& mfttracks, aod::BCsWithTimestamps const&)
  {
    auto getTracksSize = [&tracks, this](CollTable::iterator const& collision) {
      auto associatedTracks = tracks.sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), this->cache);
      return associatedTracks.size();
    };

    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, aod::collision::PosZ, decltype(getTracksSize)>;
    MixedBinning binningOnVtxAndMult{{getTracksSize}, {axisVtxZME, axisMultME}, true};
    for (auto const& [col1, col2] : soa::selfCombinations(binningOnVtxAndMult, mixingParameter, -1, col, col)) {
      if (!isEventSelected(col1) || !isEventSelected(col2)) {
        continue;
      }
      if (col1.globalIndex() == col2.globalIndex()) {
        continue;
      }
      if (col1.has_foundFT0() && col2.has_foundFT0()) {
        auto bc = col1.bc_as<aod::BCsWithTimestamps>();
        loadMultCorrection(bc.timestamp());
        auto slicedTriggerMftTracks = mfttracks.sliceBy(perColMft, col1.globalIndex());
        auto slicedTriggerTracks = tracks.sliceBy(perColGlobal, col1.globalIndex());
        const auto& ft0 = col2.foundFT0();
        auto multiplicity = 1.0f;
        if (isUseCentEst)
          multiplicity = selColCent(col1);
        else
          multiplicity = countNTracks(slicedTriggerTracks);
        float multw = getMultEffCorr(multiplicity);
        if (isUseEffCorr)
          multiplicity = multiplicity * multw;
        fillCorrFt0aMft<CorrelationContainer::kCFStepReconstructed>(mixedFt0aMft, slicedTriggerMftTracks, ft0, true, col1.posZ(), multiplicity);
      }
    }
  } // mixed event

  void processFt0aFt0cME(CollTable const& col, aod::FT0s const&, TrksTable const& tracks, aod::BCsWithTimestamps const&)
  {
    auto getTracksSize = [&tracks, this](CollTable::iterator const& collision) {
      auto associatedTracks = tracks.sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), this->cache);
      return associatedTracks.size();
    };

    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, aod::collision::PosZ, decltype(getTracksSize)>;
    MixedBinning binningOnVtxAndMult{{getTracksSize}, {axisVtxZME, axisMultME}, true};
    for (auto const& [col1, col2] : soa::selfCombinations(binningOnVtxAndMult, mixingParameter, -1, col, col)) {
      if (!isEventSelected(col1) || !isEventSelected(col2)) {
        continue;
      }
      if (col1.globalIndex() == col2.globalIndex()) {
        continue;
      }
      if (col1.has_foundFT0() && col2.has_foundFT0()) {
        auto bc = col1.bc_as<aod::BCsWithTimestamps>();
        loadMultCorrection(bc.timestamp());
        auto slicedTriggerTracks = tracks.sliceBy(perColGlobal, col1.globalIndex());
        auto multiplicity = 1.0f;
        if (isUseCentEst)
          multiplicity = selColCent(col1);
        else
          multiplicity = countNTracks(slicedTriggerTracks);
        const auto& ft0a = col1.foundFT0();
        const auto& ft0c = col2.foundFT0();
        float multw = getMultEffCorr(multiplicity);
        if (isUseEffCorr)
          multiplicity = multiplicity * multw;
        fillCorrFt0aFt0c<CorrelationContainer::kCFStepReconstructed>(mixedFt0aFt0c, ft0a, ft0c, true, col1.posZ(), multiplicity);
      }
    }
  } // mixed event

  template <typename CheckGenTrack>
  bool isGenTrackSelected(CheckGenTrack const& track)
  {
    if (!track.isPhysicalPrimary()) {
      return false;
    }
    if (!track.producedByGenerator()) {
      return false;
    }
    auto pdgTrack = pdg->GetParticle(track.pdgCode());
    if (pdgTrack == nullptr) {
      return false;
    }
    if (std::abs(pdgTrack->Charge()) < kMinCharge) {
      return false;
    }
    if (std::abs(track.eta()) >= cfgEtaCut) {
      return false;
    }
    return true;
  }

  void processEff(aod::McCollisions::iterator const& mcCollision, CollTableMC const& RecCols, aod::McParticles const& GenParticles, TrksTableMC const& RecTracks)
  {
    if (std::abs(mcCollision.posZ()) >= cfgVtxCut) {
      return;
    }

    auto multiplicity = -999.;
    auto numcontributors = -999;
    for (const auto& RecCol : RecCols) {
      if (!isEventSelected(RecCol)) {
        continue;
      }
      if (RecCol.numContrib() <= numcontributors) {
        continue;
      } else {
        numcontributors = RecCol.numContrib();
      }
      if (isUseCentEst) {
        multiplicity = selColCent(RecCol);
      } else {
        auto recTracksPart = RecTracks.sliceBy(perColMC, RecCol.globalIndex());
        multiplicity = countNTracks(recTracksPart);
      }
    }

    for (const auto& particle : GenParticles) {
      if (!isGenTrackSelected(particle)) {
        continue;
      }
      histos.fill(HIST("hmcgendndptPrimary"), particle.eta(), particle.pt(), multiplicity, mcCollision.posZ());
    } // track (mcgen) loop

    for (const auto& RecCol : RecCols) {
      if (!isEventSelected(RecCol)) {
        continue;
      }
      auto recTracksPart = RecTracks.sliceBy(perColMC, RecCol.globalIndex());
      for (const auto& Rectrack : recTracksPart) {
        if (Rectrack.has_mcParticle()) {
          auto mcpart = Rectrack.mcParticle();
          histos.fill(HIST("hmcrecdndptRecoAll"), mcpart.eta(), mcpart.pt(), multiplicity, mcCollision.posZ());
          if (mcpart.isPhysicalPrimary()) {
            histos.fill(HIST("hmcrecdndptRecoPrimary"), mcpart.eta(), mcpart.pt(), multiplicity, mcCollision.posZ());
          }
        } else {
          histos.fill(HIST("hmcrecdndptFake"), Rectrack.eta(), Rectrack.pt(), multiplicity, mcCollision.posZ());
        }
      } // track (mcrec) loop
    } // rec collision
  }

  PROCESS_SWITCH(LongrangeCorrelation, processEventStat, "event stat", false);
  PROCESS_SWITCH(LongrangeCorrelation, processFt0aGlobalSE, "same event FT0a vs global", false);
  PROCESS_SWITCH(LongrangeCorrelation, processFt0aGlobalME, "mixed event FT0a vs global", false);
  PROCESS_SWITCH(LongrangeCorrelation, processFt0cGlobalSE, "same event FT0c vs global", false);
  PROCESS_SWITCH(LongrangeCorrelation, processFt0cGlobalME, "mixed event FT0c vs global", false);
  PROCESS_SWITCH(LongrangeCorrelation, processMftGlobalSE, "same event MFT vs global", false);
  PROCESS_SWITCH(LongrangeCorrelation, processMftGlobalME, "mixed event MFT vs global", false);
  PROCESS_SWITCH(LongrangeCorrelation, processFt0aMftSE, "same event FT0a vs MFT", false);
  PROCESS_SWITCH(LongrangeCorrelation, processFt0aMftME, "mixed event FT0a vs MFT", false);
  PROCESS_SWITCH(LongrangeCorrelation, processFt0aFt0cSE, "same event FT0a vs FT0c", false);
  PROCESS_SWITCH(LongrangeCorrelation, processFt0aFt0cME, "mixed event FT0a vs FT0c", false);
  PROCESS_SWITCH(LongrangeCorrelation, processEff, "Estimate efficiency", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<LongrangeCorrelation>(cfgc)};
}
