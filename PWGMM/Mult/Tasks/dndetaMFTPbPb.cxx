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
///
/// \file dndetaMFTPbPb.cxx
/// \brief  Task for calculating dNdeta in Pb-Pb collisions using MFT detector
/// \author Gyula Bencedi, gyula.bencedi@cern.ch
/// \since  Nov 2024

#include "PWGMM/Mult/Core/include/Functions.h"
#include "PWGMM/Mult/DataModel/Index.h"
#include "PWGMM/Mult/DataModel/bestCollisionTable.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/McCollisionExtra.h"
#include "Common/DataModel/Multiplicity.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/Propagator.h>
#include <Field/MagneticField.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/DataTypes.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TGeoGlobalMagField.h>
#include <TH1.h>
#include <TH2.h>
#include <THnSparse.h>
#include <TString.h>

#include <sys/types.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::track;
using namespace o2::aod::fwdtrack;
using namespace o2::constants::physics;
using namespace o2::constants::math;
using namespace pwgmm::mult;
using namespace o2::aod::rctsel;

auto static constexpr CminCharge = 3.f;
auto static constexpr CintZero = 0;
auto static constexpr CfloatFive = 5.f;
auto static constexpr CminAccFT0A = 3.5f;
auto static constexpr CmaxAccFT0A = 4.9f;
auto static constexpr CminAccFT0C = -3.3f;
auto static constexpr CmaxAccFT0C = -2.1f;

enum class EvtSel {
  evtAll = 0,
  evtSel,
  evtIsGoodZvtx,
  evtNoSameBunchPileup,
  evtZvtxCut,
  evtNoCollInTimeRangeStd,
  evtNoCollInTimeRangeNarrow,
  evtNoCollInTimeRangeStrict,
  evtNoCollInRofStrict,
  evtNoCollInRofStandard,
  evtNoHighMultCollInPrevRof,
  evtGoodITSLayersAll,
  evtBelowMinOccup,
  evtAboveMaxOccup,
  evtRCTFlagChecker,
  evtRCTFlagCheckerExtra,
  nEvtSel
};

enum class TrkSel {
  trkSelAll = 0,
  trkSelNCls,
  trkSelChi2Ncl,
  trkSelEta,
  trkSelPhiCut,
  trkSelPt,
  trkSelCA,
  nTrkSel
};

enum class TrkBestSel {
  trkBestSelAll = 0,
  trkBestSelCollID,
  trkBestSelOrphan,
  trkBestSelDCAxyCut,
  trkBestSelDCAzCut,
  trkBestSelNumReassoc,
  nTrkBestSel
};

enum class AmbTrkType {
  kAll = 0,
  kNonAmb,
  kOrphan,
  kOrphanNull,
  kNonAmbSame,
  kAmb,
  kAmbGt1,
  nAmbTrkType
};

enum class GenTrkType {
  kGenAll = 0,
  kGenRecEvt,
  nGenTrkType
};

enum class GenIdxTrkType {
  kGenIdxAll = 0,
  kGenIdxDupl,
  kGenIdxFake,
  nGenIdxTrkType
};

enum class RecTrkType {
  kRecAll = 0,
  kRecFake,
  nRecTrkType
};

enum class RecIdxTrkType {
  kRecIdxPrim = 0,
  kRecIdxSec,
  kRecIdxDupl,
  nRecIdxTrkType
};

enum class EvtLossType {
  kGenAll = 0,
  kGenSel,
  kGenSplit,
  kGenRecEvt,
  nEvtLossType
};

enum class McEffStatus {
  kMcEffAll = 0,
  kMcEffSel,
  kMcEffHasMcColl,
  kMcEffNoSplitVtx,
  nMcEfftStatus
};

enum class McStatus {
  kMcRecAll = 0,
  kMcRecSel,
  kMcRecHasMcColl,
  kMcRecNoSplitVtx,
  kMcGenAll,
  nMcStatus
};

enum class McTrackStatus {
  kBestTrkAll = 0,
  kBesTrktSel,
  kTrkAsBestSel,
  kTrkHasColl,
  kTrkReassignedRemoved,
  nMcTrackStatus
};

struct DndetaMFTPbPb {
  SliceCache cache;

  enum OccupancyEst { TrkITS = 1,
                      Ft0C };

  HistogramRegistry registryData{"registryData", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryMC{"registryMC", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryQC{"registryQC", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  struct : ConfigurableGroup {
    Configurable<bool> cfgDoIR{"cfgDoIR", false, "Flag to retrieve Interaction rate from CCDB"};
    Configurable<bool> cfgUseIRCut{"cfgUseIRCut", false, "Flag to cut on IR rate"};
    Configurable<bool> cfgIRCrashOnNull{"cfgIRCrashOnNull", false, "Flag to avoid CTP RateFetcher crash"};
    Configurable<std::string> cfgIRSource{"cfgIRSource", "ZNC hadronic", "Estimator of the interaction rate (Pb-Pb: ZNC hadronic)"};
    Configurable<bool> cfgUseTrackSel{"cfgUseTrackSel", false, "Flag to apply track selection"};
    Configurable<bool> cfgUseParticleSel{"cfgUseParticleSel", true, "Flag to apply particle selection"};
    Configurable<bool> cfgUsePrimaries{"cfgUsePrimaries", true, "Select primary particles"};
    Configurable<bool> cfgUseSecondaries{"cfgUseSecondaries", false, "Select secondary particles"};
    Configurable<bool> cfgRemoveReassigned{"cfgRemoveReassigned", false, "Remove reassgined tracks"};
    Configurable<bool> cfgRemoveSplitVertex{"cfgRemoveSplitVertex", true, "Remove split vertices"};
    Configurable<bool> cfgUseTrackParExtra{"cfgUseTrackParExtra", false, "Use table with refitted track parameters"};
    Configurable<bool> cfgUseInelgt0wMFT{"cfgUseInelgt0wMFT", false, "Use INEL > 0 condition with MFT acceptance"};
    Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  } gConf;

  struct : ConfigurableGroup {
    ConfigurableAxis interactionRateBins{"interactionRateBins", {500, 0, 50}, "Binning for the interaction rate (kHz)"};
    ConfigurableAxis occupancyBins{"occupancyBins", {VARIABLE_WIDTH, 0.0f, 250.0f, 500.0f, 750.0f, 1000.0f, 1500.0f, 2000.0f, 3000.0f, 4500.0f, 6000.0f, 8000.0f, 10000.0f, 50000.0f}, "Occupancy"};
    ConfigurableAxis centralityBins{"centralityBins", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100}, "Centrality"};
    Configurable<std::vector<float>> genMultBins{"genMultBins", std::vector<float>{500, 300, 200, 120, 80, 50, 30, 10, 0}, "Generated multiplicity low values matching reco cent bins"};
    Configurable<std::vector<float>> recCentBinCenters{"recCentBinCenters", std::vector<float>{2.5, 7.5, 15.0, 25.0, 35.0, 45.0, 55.0, 75.0, 95.0}, "Reco cent bin centers"};
    ConfigurableAxis irBins{"irBins", {500, 0, 50}, "Interaction rate (kHz)"};
    ConfigurableAxis pvBins{"pvBins", {501, -0.5, 500.5}, ""};
    ConfigurableAxis fv0aMultBins{"fv0aMultBins", {501, -0.5, 500.5}, ""};
    ConfigurableAxis ft0aMultBins{"ft0aMultBins", {501, -0.5, 500.5}, ""};
    ConfigurableAxis ft0cMultBins{"ft0cMultBins", {501, -0.5, 500.5}, ""};
    ConfigurableAxis ptBins{"ptBins", {101, -0.5, 10.5}, "pT binning (GeV/c)"};
    ConfigurableAxis multBins{"multBins", {1001, -0.5, 1000.5}, "Multiplicity binning"};
    ConfigurableAxis zvtxBins{"zvtxBins", {60, -30., 30.}, "Z-vtx binning (cm)"};
    ConfigurableAxis deltaZBins{"deltaZBins", {800, -10., 10.}, "Delta Z-vtx binning (cm)"};
    ConfigurableAxis dcaXYBins{"dcaXYBins", {800, -1., 1.}, "DCAxy binning (cm)"};
    ConfigurableAxis dcaZBins{"dcaZBins", {800, -1., 1.}, "DCAz binning (cm)"};
    ConfigurableAxis phiBins{"phiBins", {629, 0., TwoPI}, "#varphi binning (rad)"};
    ConfigurableAxis etaBins{"etaBins", {20, -4., -2.}, "#eta binning"};
    ConfigurableAxis chiSqPerNclBins{"chiSqPerNclBins", {100, 0, 100}, "#chi^{2} binning"};
    ConfigurableAxis nClBins{"nClBins", {10, 0.5, 10.5}, "number of clusters binning"};
    ConfigurableAxis tanLambdaBins{"tanLambdaBins", {100, -25, 0}, "binning for tan(lambda)"};
    ConfigurableAxis invQPtBins{"invQPtBins", {200, -10, 10}, "binning for q/p_{T}"};
  } binOpt;

  struct : ConfigurableGroup {
    Configurable<bool> requireRCTFlagChecker{"requireRCTFlagChecker", false, "Check event quality in run condition table"};
    Configurable<bool> requireRCTFlagCheckerExtra{"requireRCTFlagCheckerExtra", true, "Check RCT flag extra"};
    Configurable<std::string> cfgEvtRCTFlagCheckerLabel{"cfgEvtRCTFlagCheckerLabel", "CBT_fw", "Evt sel: RCT flag checker label"};
    Configurable<bool> cfgEvtRCTFlagCheckerZDCCheck{"cfgEvtRCTFlagCheckerZDCCheck", false, "Evt sel: RCT flag checker ZDC check"};
    Configurable<bool> cfgEvtRCTFlagCheckerLimitAcceptAsBad{"cfgEvtRCTFlagCheckerLimitAcceptAsBad", true, "Evt sel: RCT flag checker treat Limited Acceptance As Bad"};
  } rctCuts;

  struct : ConfigurableGroup {
    Configurable<bool> usephiCut{"usephiCut", false, "use azimuthal angle cut"};
    Configurable<float> phiCut{"phiCut", 0.1f, "Cut on azimuthal angle of MFT tracks"};
    Configurable<float> minPhi{"minPhi", 0.f, ""};
    Configurable<float> maxPhi{"maxPhi", 6.2832, ""};
    Configurable<float> minEta{"minEta", -3.6f, ""};
    Configurable<float> maxEta{"maxEta", -2.5f, ""};
    Configurable<int> minNclusterMft{"minNclusterMft", 5, "minimum number of MFT clusters"};
    Configurable<bool> useChi2Cut{"useChi2Cut", true, "use track chi2 cut"};
    Configurable<float> maxChi2NCl{"maxChi2NCl", 40.0f, "maximum chi2 per MFT clusters"};
    Configurable<bool> usePtCut{"usePtCut", false, "use track pT cut"};
    Configurable<float> minPt{"minPt", 0., "minimum pT of the MFT tracks"};
    Configurable<bool> requireCA{"requireCA", false, "Use Cellular Automaton track-finding algorithm"};
    Configurable<float> maxDCAxy{"maxDCAxy", 1.0f, "Cut on dca XY"};
    Configurable<bool> useDCAzCut{"useDCAzCut", true, "use dca Z cut"};
    Configurable<float> maxDCAz{"maxDCAz", 1.0f, "Cut on dca Z"};
    Configurable<int> selMcMask{"selMcMask", 0, "McMask for correct match"};
  } trackCuts;

  struct : ConfigurableGroup {
    Configurable<float> minCentrality{"minCentrality", 0.0f, "minimum centrality selection"};
    Configurable<float> maxCentrality{"maxCentrality", 100.0f, "maximum centrality selection"};
    Configurable<float> maxZvtx{"maxZvtx", 10.0f, "maximum cut on z-vtx (cm)"};
    Configurable<float> minZvtx{"minZvtx", -10.0f, "minimum cut on z-vtx (cm)"};
    Configurable<bool> useZDiffCut{"useZDiffCut", false, "use Zvtx reco-mc diff. cut"};
    Configurable<float> maxZvtxDiff{"maxZvtxDiff", 1.0f, "max allowed Z vtx difference for reconstruced collisions (cm)"};
    Configurable<bool> useZVtxCutMC{"useZVtxCutMC", false, "use Zvtx cut in MC"};
    Configurable<bool> requireIsGoodZvtxFT0VsPV{"requireIsGoodZvtxFT0VsPV", true, "require events with PV position along z consistent (within 1 cm) between PV reconstructed using tracks and PV using FT0 A-C time difference"};
    Configurable<bool> requireRejectSameBunchPileup{"requireRejectSameBunchPileup", true, "reject collisions in case of pileup with another collision in the same foundBC"};
    Configurable<bool> requireNoCollInTimeRangeStrict{"requireNoCollInTimeRangeStrict", false, " requireNoCollInTimeRangeStrict"};
    Configurable<bool> requireNoCollInRofStrict{"requireNoCollInRofStrict", true, "requireNoCollInRofStrict"};
    Configurable<bool> requireNoCollInRofStandard{"requireNoCollInRofStandard", true, "requireNoCollInRofStandard"};
    Configurable<bool> requireNoHighMultCollInPrevRof{"requireNoHighMultCollInPrevRof", true, "requireNoHighMultCollInPrevRof"};
    Configurable<bool> requireGoodITSLayersAll{"requireGoodITSLayersAll", true, "requireGoodITSLayersAll"};
    Configurable<bool> requireNoCollInTimeRangeStd{"requireNoCollInTimeRangeStd", true, "reject collisions corrupted by the cannibalism, with other collisions within +/- 10 microseconds"};
    Configurable<bool> requireNoCollInTimeRangeNarrow{"requireNoCollInTimeRangeNarrow", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 10 microseconds"};
    Configurable<uint> occupancyEstimator{"occupancyEstimator", 1, "Occupancy estimator: 1 = trackOccupancyInTimeRange, 2 = ft0cOccupancyInTimeRange"};
    Configurable<float> minOccupancy{"minOccupancy", -1, "minimum occupancy from neighbouring collisions"};
    Configurable<float> maxOccupancy{"maxOccupancy", -1, "maximum occupancy from neighbouring collisions"};
    Configurable<float> minIR{"minIR", -1, "minimum IR (kHz) collisions"};
    Configurable<float> maxIR{"maxIR", -1, "maximum IR (kHz) collisions"};
    Configurable<bool> useInelgt0wTVX{"useInelgt0wTVX", false, "Use INEL > 0 condition with TVX trigger, i.e. FT0A and FT0C acceptance"};
    Configurable<bool> useGenMult{"useGenMult", false, "use MultMC for centrality"};
  } eventCuts;

  Service<o2::framework::O2DatabasePDG> pdg{};
  Service<ccdb::BasicCCDBManager> ccdb{};
  Configurable<int64_t> ccdbNoLaterThan{"ccdbNoLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  int mRunNumber{-1};
  uint64_t mSOR{0};
  float mMinSeconds{-1.};
  std::unordered_map<int, TH2*> gHadronicRate;
  ctpRateFetcher rateFetcher;
  TH2* gCurrentHadronicRate = nullptr;
  RCTFlagsChecker rctChecker;
  RCTFlagsChecker rctCheckerExtra{kFT0Bad, kITSBad, kTPCBadTracking, kMFTBad};

  float bZ = 0;                                 // Magnetic field for MFT
  std::array<double, 3> centerMFT{0, 0, -61.4}; // Field at center of MFT
  o2::parameters::GRPMagField* grpmag = nullptr;

  std::vector<int> ambiguousTrkIds;
  std::vector<int> reassignedTrkIds;

  /// @brief init function, definition of histograms
  void init(InitContext&)
  {
    const AxisSpec pvAxis = {binOpt.pvBins, "PV", "PV axis"};
    const AxisSpec multFV0aAxis = {binOpt.fv0aMultBins, "fv0a", "FV0AMult axis"};
    const AxisSpec multFT0aAxis = {binOpt.ft0aMultBins, "ft0a", "FT0AMult axis"};
    const AxisSpec multFT0cAxis = {binOpt.ft0cMultBins, "ft0c", "FT0CMult axis"};
    const AxisSpec centralityAxis = {binOpt.centralityBins, "Centrality", "centrality axis"};
    const AxisSpec occupancyAxis = {binOpt.occupancyBins, "Occupancy", "occupancy axis"};
    const AxisSpec irAxis = {binOpt.interactionRateBins, "Interaction Rate", "IR axis"};
    const AxisSpec ptAxis = {binOpt.ptBins, "Pt axis (GeV/c)"};
    const AxisSpec multAxis = {binOpt.multBins, "N_{trk} axis"};
    const AxisSpec zAxis = {binOpt.zvtxBins, "Z-vtx axis"};
    const AxisSpec deltaZAxis = {binOpt.deltaZBins, "Delta Z-vtx axis"};
    const AxisSpec dcaxyAxis = {binOpt.dcaXYBins, "DCA-xy axis"};
    const AxisSpec dcazAxis = {binOpt.dcaZBins, "DCA-z axis"};
    const AxisSpec phiAxis = {binOpt.phiBins, "#phi axis"};
    const AxisSpec etaAxis = {binOpt.etaBins, "#eta axis"};
    const AxisSpec chiSqAxis = {binOpt.chiSqPerNclBins, "Chi2 axis"};
    const AxisSpec nclsAxis = {binOpt.nClBins, "Number of clusters axis"};
    const AxisSpec tanLambdaAxis{binOpt.tanLambdaBins, "tan(#lambda)"};
    const AxisSpec invQPtAxis{binOpt.invQPtBins, "q/p_{T} (1/GeV)"};
    const AxisSpec genEvtTypeAxis = {static_cast<int>(GenTrkType::nGenTrkType) - 1, +0.5, +static_cast<int>(GenTrkType::nGenTrkType) - 0.5, "", "Gen-Trk Type axis"};
    const AxisSpec recEvtTypeAxis = {static_cast<int>(RecTrkType::nRecTrkType) - 1, +0.5, +static_cast<int>(RecTrkType::nRecTrkType) - 0.5, "", "Rec-Trk Type axis"};
    const AxisSpec genEvtIdxTypeAxis = {static_cast<int>(GenTrkType::nGenTrkType) - 1, +0.5, +static_cast<int>(GenTrkType::nGenTrkType) - 0.5, "", "Gen-Trk Type axis"};
    const AxisSpec recEvtIdxTypeAxis = {static_cast<int>(RecTrkType::nRecTrkType) - 1, +0.5, +static_cast<int>(RecTrkType::nRecTrkType) - 0.5, "", "Rec-Trk Type axis"};
    const AxisSpec evtLossTypeAxis = {static_cast<int>(EvtLossType::nEvtLossType) - 1, +0.5, +static_cast<int>(EvtLossType::nEvtLossType) - 0.5, "", "Evt-loss Type axis"};

    rctChecker.init(rctCuts.cfgEvtRCTFlagCheckerLabel, rctCuts.cfgEvtRCTFlagCheckerZDCCheck, rctCuts.cfgEvtRCTFlagCheckerLimitAcceptAsBad);

    ccdb->setURL(ccdbUrl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(ccdbNoLaterThan.value);
    ccdb->setFatalWhenNull(false);

    if (static_cast<int>(doprocessDataInclusive) + static_cast<int>(doprocessDatawBestTracksInclusive) > 1) {
      LOGP(fatal, "Either processDataInclusive OR processDatawBestTracksInclusive should be enabled!");
    }
    if (static_cast<int>(doprocessDataCentFT0C) + static_cast<int>(doprocessDatawBestTracksCentFT0C) > 1) {
      LOGP(fatal, "Either processDataCent[ESTIMATOR] OR processDatawBestTracksCent[ESTIMATOR] should be enabled!");
    }
    if (static_cast<int>(doprocessMCInclusive) + static_cast<int>(doprocessMCbestInclusive) > 1) {
      LOGP(fatal, "Either processMCInclusive OR processMCbestInclusive should be enabled!");
    }
    if (static_cast<int>(doprocessMCCentFT0C) + static_cast<int>(doprocessMCbestCentFT0C) > 1) {
      LOGP(fatal, "Either processMCCent[ESTIMATOR] OR processMCbestCent[ESTIMATOR] should be enabled!");
    }
    if (static_cast<int>(doprocessEfficiencyInclusive) + static_cast<int>(doprocessEfficiencyBestInclusive) > 1) {
      LOGP(fatal, "Either doprocessEfficiencyInclusive OR doprocessEfficiencyBestInclusive should be enabled!");
    }
    if (static_cast<int>(doprocessEfficiencyCentFT0C) + static_cast<int>(doprocessEfficiencyBestCentFT0C) > 1) {
      LOGP(fatal, "Either doprocessEfficiencyCentFT0C OR doprocessEfficiencyBestCentFT0C should be enabled!");
    }
    if (static_cast<int>(doprocessEfficiencyIdxInlusive) + static_cast<int>(doprocessEfficiencyIdxBestInlusive) > 1) {
      LOGP(fatal, "Either doprocessEfficiencyIdxInlusive OR doprocessEfficiencyIdxBestInlusive should be enabled!");
    }
    if (static_cast<int>(doprocessEfficiencyIdxCentFT0C) + static_cast<int>(doprocessEfficiencyIdxBestCentFT0C) > 1) {
      LOGP(fatal, "Either doprocessEfficiencyIdxCentFT0C OR doprocessEfficiencyIdxBestCentFT0C should be enabled!");
    }

    // registry.add("Events/hRCTSel", "Event accepted if RCT not selected;RCT Status;", {HistType::kTH1F, {{2, 0.5, 2.5}}});
    // auto hrctSel = registry.get<TH1>(HIST("Events/hRCTSel"));
    // auto* x = hrctSel->GetXaxis();
    // x->SetBinLabel(1, "All");
    // x->SetBinLabel(2, "kFT0Bad");

    // General counters - QC
    auto hBcSel = registryQC.add<TH1>("hBcSel", "hBcSel", HistType::kTH1F, {{3, -0.5f, +2.5f}});
    hBcSel->GetXaxis()->SetBinLabel(1, "Good BCs");
    hBcSel->GetXaxis()->SetBinLabel(2, "BCs with collisions");
    hBcSel->GetXaxis()->SetBinLabel(3, "BCs with pile-up/splitting");

    registryQC.add("Events/hEvtSel", "Number of events; Cut; #Evt Passed Cut", {HistType::kTH1F, {{static_cast<int>(EvtSel::nEvtSel), -0.5, +static_cast<int>(EvtSel::nEvtSel) - 0.5}}});
    std::array<std::string_view, static_cast<int>(EvtSel::nEvtSel)> labelEvtSel{
      "All coll.",
      "Sel 8",
      "kIsGoodZvtxFT0vsPV",
      "NoSameBunchPileup",
      "Z-vtx cut",
      "kNoCollInTimeRangeStd",
      "kNoCollInTimeRangeNarrow",
      "kNoCollInTimeRangeStrict",
      "kNoCollInRofStrict",
      "kNoCollInRofStandard",
      "kNoHighMultCollInPrevRof",
      "kIsGoodITSLayersAll",
      "Below min occup.",
      "Above max occup.",
      "RCT Flag Checker",
      "RCT Flag Checker Extra"};
    registryQC.get<TH1>(HIST("Events/hEvtSel"))->SetMinimum(0.1);
    for (int iBin = 0; iBin < static_cast<int>(EvtSel::nEvtSel); iBin++) {
      registryQC.get<TH1>(HIST("Events/hEvtSel"))->GetXaxis()->SetBinLabel(iBin + 1, labelEvtSel[iBin].data());
    }

    registryQC.add("Tracks/hTrkSel", "Number of tracks; Cut; #Tracks Passed Cut", {HistType::kTH1F, {{static_cast<int>(TrkSel::nTrkSel), -0.5, +static_cast<int>(TrkSel::nTrkSel) - 0.5}}});
    std::array<std::string_view, static_cast<int>(TrkSel::nTrkSel)> labelTrkSel{
      "All",
      "Ncls",
      "Chi2",
      "Eta",
      "Phi cut",
      "Pt",
      "CA"};
    registryQC.get<TH1>(HIST("Tracks/hTrkSel"))->SetMinimum(0.1);
    for (int iBin = 0; iBin < static_cast<int>(TrkSel::nTrkSel); iBin++) {
      registryQC.get<TH1>(HIST("Tracks/hTrkSel"))->GetXaxis()->SetBinLabel(iBin + 1, labelTrkSel[iBin].data());
    }

    if (doprocessDatawBestTracksInclusive || doprocessDatawBestTracksCentFT0C ||
        doprocessMCbestInclusive || doprocessMCbestCentFT0C ||
        doprocessEfficiencyBestInclusive || doprocessEfficiencyBestCentFT0C ||
        doprocessEfficiencyIdxBestInlusive || doprocessEfficiencyIdxBestCentFT0C) {
      registryQC.add("Tracks/hBestTrkSel", "Number of best tracks; Cut; #Tracks Passed Cut", {HistType::kTH1F, {{static_cast<int>(TrkBestSel::nTrkBestSel), -0.5, +static_cast<int>(TrkBestSel::nTrkBestSel) - 0.5}}});
      std::array<std::string_view, static_cast<int>(TrkBestSel::nTrkBestSel)> labelTrkTrkBestSel{
        "All",
        "Assigned (ID>=0)",
        "No orphans",
        "DCA xy cut",
        "DCA z cut",
        "#Reassoc"};
      registryQC.get<TH1>(HIST("Tracks/hBestTrkSel"))->SetMinimum(0.1);
      for (int iBin = 0; iBin < static_cast<int>(TrkBestSel::nTrkBestSel); iBin++) {
        registryQC.get<TH1>(HIST("Tracks/hBestTrkSel"))->GetXaxis()->SetBinLabel(iBin + 1, labelTrkTrkBestSel[iBin].data());
      }
    }
    if (doprocessDataInclusive || doprocessDatawBestTracksInclusive) {
      registryData.add({"Events/hInteractionRate", "; IR (kHz); occupancy", {HistType::kTH2F, {irAxis, occupancyAxis}}});
      registryData.add({"Events/Selection", ";status; occupancy", {HistType::kTH2F, {{2, 0.5, 2.5}, occupancyAxis}}});
      auto hstat = registryData.get<TH2>(HIST("Events/Selection"));
      auto* x = hstat->GetXaxis();
      x->SetBinLabel(1, "All");
      x->SetBinLabel(2, "Selected");

      registryData.add({"Tracks/EtaZvtx", "; #eta; #it{z}_{vtx} (cm); occupancy", {HistType::kTHnSparseF, {etaAxis, zAxis, occupancyAxis}}});
      registryData.add({"Tracks/PhiEta", "; #varphi; #eta; occupancy", {HistType::kTHnSparseF, {phiAxis, etaAxis, occupancyAxis}}});
      registryData.add({"Events/NtrkZvtx", "; N_{trk}; Z_{vtx} (cm); occupancy", {HistType::kTHnSparseF, {multAxis, zAxis, occupancyAxis}}});
      registryData.add({"Tracks/Chi2Eta", "; #chi^{2}; #eta; occupancy", {HistType::kTHnSparseF, {chiSqAxis, etaAxis, occupancyAxis}}});
      registryData.add({"Tracks/NclustersEta", "; nClusters; #eta; occupancy", {HistType::kTHnSparseF, {nclsAxis, etaAxis, occupancyAxis}}});
      registryData.add({"Tracks/TanLambda", "; TanLambda; occupancy", {HistType::kTH2F, {tanLambdaAxis, occupancyAxis}}});
      registryData.add({"Tracks/InvQPt", "; InvQPt; occupancy", {HistType::kTH2F, {invQPtAxis, occupancyAxis}}});

      if (doprocessDatawBestTracksInclusive) {
        registryData.add({"Tracks/DCA3d", "; p_{T} (GeV/c); #eta; DCA_{XY} (cm); DCA_{Z} (cm); occupancy", {HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, occupancyAxis}}});
        registryData.add({"Tracks/ReTracksEtaZvtx", "; #eta; #it{z}_{vtx} (cm); occupancy", {HistType::kTHnSparseF, {etaAxis, zAxis, occupancyAxis}}});
        registryData.add({"Tracks/ReTracksPhiEta", "; #varphi; #eta; occupancy", {HistType::kTHnSparseF, {phiAxis, etaAxis, occupancyAxis}}});
        registryData.add({"Tracks/OrigTracksEtaZvtx", "; #eta; #it{z}_{vtx} (cm); occupancy", {HistType::kTHnSparseF, {etaAxis, zAxis, occupancyAxis}}});
        registryData.add({"Tracks/OrigTracksPhiEta", "; #varphi; #eta; occupancy", {HistType::kTHnSparseF, {phiAxis, etaAxis, occupancyAxis}}});
        registryData.add({"Tracks/RestTracksEtaZvtx", "; #eta; #it{z}_{vtx} (cm); occupancy", {HistType::kTHnSparseF, {etaAxis, zAxis, occupancyAxis}}});
        registryData.add({"Tracks/RestTracksPhiEta", "; #varphi; #eta; occupancy", {HistType::kTHnSparseF, {phiAxis, etaAxis, occupancyAxis}}});
        registryData.add({"Tracks/TrackAmbDegree", "; N_{coll}^{comp}; occupancy", {HistType::kTH2F, {{51, -0.5, 50.5}, occupancyAxis}}});
        registryData.add({"Tracks/TanLambdaExtra", "; TanLambda; occupancy", {HistType::kTH2F, {tanLambdaAxis, occupancyAxis}}});
        registryData.add({"Tracks/InvQPtExtra", "; InvQPt; occupancy", {HistType::kTH2F, {invQPtAxis, occupancyAxis}}});
        registryData.add({"Tracks/EtaExtra", "; #eta; occupancy", {HistType::kTH2F, {etaAxis, occupancyAxis}}});
        registryData.add({"Tracks/PhiExtra", "; #varphi; occupancy", {HistType::kTH2F, {phiAxis, occupancyAxis}}});
      }
    }
    if (doprocessDataCentFT0C || doprocessDatawBestTracksCentFT0C) {
      registryData.add({"Events/Centrality/hInteractionRate", "; IR (kHz); centrality; occupancy", {HistType::kTHnSparseF, {irAxis, centralityAxis, occupancyAxis}}});
      registryData.add({"Events/Centrality/Selection", ";status; centrality; occupancy", {HistType::kTHnSparseF, {{2, 0.5, 2.5}, centralityAxis, occupancyAxis}}});
      auto hstat = registryData.get<TH2>(HIST("Events/Centrality/Selection"));
      hstat->GetXaxis()->SetBinLabel(1, "All");
      hstat->GetXaxis()->SetBinLabel(2, "Selected");

      registryData.add({"Tracks/Centrality/EtaZvtx", "; #eta; #it{z}_{vtx} (cm); centrality; occupancy", {HistType::kTHnSparseF, {etaAxis, zAxis, centralityAxis, occupancyAxis}}});
      registryData.add({"Tracks/Centrality/PhiEta", "; #varphi; #eta; centrality; occupancy", {HistType::kTHnSparseF, {phiAxis, etaAxis, centralityAxis, occupancyAxis}}});
      registryData.add({"Events/Centrality/NtrkZvtx", "; N_{trk}; Z_{vtx} (cm); centrality; occupancy", {HistType::kTHnSparseF, {multAxis, zAxis, centralityAxis, occupancyAxis}}});
      registryData.add({"Tracks/Centrality/Chi2Eta", "; #chi^{2}; #eta; centrality; occupancy", {HistType::kTHnSparseF, {chiSqAxis, etaAxis, centralityAxis, occupancyAxis}}});
      registryData.add({"Tracks/Centrality/NclustersEta", "; nClusters; #eta; centrality; occupancy", {HistType::kTHnSparseF, {nclsAxis, etaAxis, centralityAxis, occupancyAxis}}});
      registryData.add({"Tracks/Centrality/TanLambda", "; TanLambda; centrality; occupancy", {HistType::kTHnSparseF, {tanLambdaAxis, centralityAxis, occupancyAxis}}});
      registryData.add({"Tracks/Centrality/InvQPt", "; InvQPt; centrality; occupancy", {HistType::kTHnSparseF, {invQPtAxis, centralityAxis, occupancyAxis}}});
      registryData.add({"Events/Centrality/hZvtxCent", "; #it{z}_{vtx} (cm); centrality; occupancy", {HistType::kTHnSparseF, {zAxis, centralityAxis, occupancyAxis}}});

      if (doprocessDatawBestTracksCentFT0C) {
        registryData.add({"Tracks/Centrality/DCA3d", "; p_{T} (GeV/c); #eta; DCA_{XY} (cm); DCA_{Z} (cm); centrality; occupancy", {HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, centralityAxis, occupancyAxis}}});
        registryData.add({"Tracks/Centrality/ReTracksEtaZvtx", "; #eta; #it{z}_{vtx} (cm); centrality; occupancy", {HistType::kTHnSparseF, {etaAxis, zAxis, centralityAxis, occupancyAxis}}});
        registryData.add({"Tracks/Centrality/ReTracksPhiEta", "; #varphi; #eta; centrality; occupancy", {HistType::kTHnSparseF, {phiAxis, etaAxis, centralityAxis, occupancyAxis}}});
        registryData.add({"Tracks/Centrality/OrigTracksEtaZvtx", "; #eta; #it{z}_{vtx} (cm); centrality; occupancy", {HistType::kTHnSparseF, {etaAxis, zAxis, centralityAxis, occupancyAxis}}});
        registryData.add({"Tracks/Centrality/OrigTracksPhiEta", "; #varphi; #eta; centrality; occupancy", {HistType::kTHnSparseF, {phiAxis, etaAxis, centralityAxis, occupancyAxis}}});
        registryData.add({"Tracks/Centrality/RestTracksEtaZvtx", "; #eta; #it{z}_{vtx} (cm); centrality; occupancy", {HistType::kTHnSparseF, {etaAxis, zAxis, centralityAxis, occupancyAxis}}});
        registryData.add({"Tracks/Centrality/RestTracksPhiEta", "; #varphi; #eta; centrality; occupancy", {HistType::kTHnSparseF, {phiAxis, etaAxis, centralityAxis, occupancyAxis}}});
        registryData.add({"Tracks/Centrality/TrackAmbDegree", "; N_{coll}^{comp}; centrality; occupancy", {HistType::kTHnSparseF, {{51, -0.5, 50.5}, centralityAxis, occupancyAxis}}});
        registryData.add({"Tracks/Centrality/TanLambdaExtra", "; TanLambda; centrality; occupancy", {HistType::kTHnSparseF, {tanLambdaAxis, centralityAxis, occupancyAxis}}});
        registryData.add({"Tracks/Centrality/InvQPtExtra", "; InvQPt; centrality; occupancy", {HistType::kTHnSparseF, {invQPtAxis, centralityAxis, occupancyAxis}}});
        registryData.add({"Tracks/Centrality/EtaExtra", "; #eta; centrality; occupancy", {HistType::kTHnSparseF, {etaAxis, centralityAxis, occupancyAxis}}});
        registryData.add({"Tracks/Centrality/PhiExtra", "; #varphi; centrality; occupancy", {HistType::kTHnSparseF, {phiAxis, centralityAxis, occupancyAxis}}});
      }
    }
    if (doprocessMCInclusive || doprocessMCbestInclusive) {
      registryMC.add("Events/McStatus", "Number of events; Cut; occupancy", {HistType::kTH2F, {{static_cast<int>(McStatus::nMcStatus), -0.5, +static_cast<int>(McStatus::nMcStatus) - 0.5}, occupancyAxis}});
      std::array<std::string_view, static_cast<int>(McStatus::nMcStatus)> labelMcStatus{
        "Rec all",
        "Rec sel",
        "Rec w/ mc coll",
        "Rec w/o split vtx",
        "Gen all"};
      for (int iBin = 0; iBin < static_cast<int>(McStatus::nMcStatus); iBin++) {
        registryMC.get<TH2>(HIST("Events/McStatus"))->GetXaxis()->SetBinLabel(iBin + 1, labelMcStatus[iBin].data());
      }
      if (doprocessMCbestInclusive) {
        registryMC.add({"Tracks/hMcTrackStatus", "Number of tracks; Cut; occupancy", {HistType::kTH2F, {{static_cast<int>(McTrackStatus::nMcTrackStatus), -0.5, +static_cast<int>(McTrackStatus::nMcTrackStatus) - 0.5}, occupancyAxis}}});
        std::array<std::string_view, static_cast<int>(McTrackStatus::nMcTrackStatus)> labelMcTrackStatus{
          "Best all",
          "Best sel",
          "Trk sel",
          "Has coll",
          "Reas rm"};
        for (int iBin = 0; iBin < static_cast<int>(McTrackStatus::nMcTrackStatus); iBin++) {
          registryMC.get<TH2>(HIST("Tracks/hMcTrackStatus"))->GetXaxis()->SetBinLabel(iBin + 1, labelMcTrackStatus[iBin].data());
        }
        registryMC.add({"Tracks/DCA3d", "; p_{T} (GeV/c); #eta; DCA_{XY} (cm); DCA_{Z} (cm); occupancy", {HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, occupancyAxis}}});
        registryMC.add({"Tracks/TrackAmbDegree", "; N_{coll}^{comp}; occupancy", {HistType::kTH2F, {{51, -0.5, 50.5}, occupancyAxis}}});
      }
      registryMC.add({"Tracks/EtaZvtx", "; #eta; #it{z}_{vtx} (cm); occupancy", {HistType::kTHnSparseF, {etaAxis, zAxis, occupancyAxis}}});
      registryMC.add({"Tracks/EtaZvtxGen", "; #eta; #it{z}_{vtx} (cm); occupancy", {HistType::kTHnSparseF, {etaAxis, zAxis, occupancyAxis}}});
      registryMC.add({"Tracks/PhiEta", "; #varphi; #eta; occupancy", {HistType::kTHnSparseF, {phiAxis, etaAxis, occupancyAxis}}});
      registryMC.add({"Tracks/PhiEtaGen", "; #varphi; #eta; occupancy", {HistType::kTHnSparseF, {phiAxis, etaAxis, occupancyAxis}}});
      registryMC.add({"Events/NtrkZvtxGen_t", "; N_{trk}; #it{z}_{vtx} (cm); occupancy", {HistType::kTHnSparseF, {multAxis, zAxis, occupancyAxis}}});
      registryMC.add({"Events/NtrkZvtxGen", "; N_{trk}; #it{z}_{vtx} (cm); occupancy", {HistType::kTHnSparseF, {multAxis, zAxis, occupancyAxis}}});
      registryMC.add({"Tracks/NclustersEta", "; nClusters; #eta; occupancy", {HistType::kTHnSparseF, {nclsAxis, etaAxis, occupancyAxis}}});
      registryMC.add({"Events/NotFoundEventZvtx", "; #it{z}_{vtx} (cm); occupancy", {HistType::kTH2F, {zAxis, occupancyAxis}}});
      registryMC.add({"Events/ZvtxDiff", "; Z_{rec} - Z_{gen} (cm); occupancy", {HistType::kTH2F, {deltaZAxis, occupancyAxis}}});
    }
    if (doprocessMCCentFT0C || doprocessMCbestCentFT0C) {
      registryMC.add("Events/Centrality/McStatus", "Number of events; Cut; centrality; occupancy", {HistType::kTHnSparseF, {{static_cast<int>(McStatus::nMcStatus), -0.5, +static_cast<int>(McStatus::nMcStatus) - 0.5}, centralityAxis, occupancyAxis}});
      std::array<std::string_view, static_cast<int>(McStatus::nMcStatus)> labelMcStatusCent{
        "Rec all",
        "Rec sel",
        "Rec w/ mc coll",
        "Rec w/o split vtx",
        "Gen all"};
      for (int iBin = 0; iBin < static_cast<int>(McStatus::nMcStatus); iBin++) {
        registryMC.get<THnSparse>(HIST("Events/Centrality/McStatus"))->GetAxis(0)->SetBinLabel(iBin + 1, labelMcStatusCent[iBin].data());
      }
      if (doprocessMCbestCentFT0C) {
        registryMC.add({"Tracks/Centrality/hMcTrackStatus", "Number of tracks; Cut; centrality; occupancy", {HistType::kTHnSparseF, {{static_cast<int>(McTrackStatus::nMcTrackStatus), -0.5, +static_cast<int>(McTrackStatus::nMcTrackStatus) - 0.5}, centralityAxis, occupancyAxis}}});
        std::array<std::string_view, static_cast<int>(McTrackStatus::nMcTrackStatus)> labelMcStatusCentBest{
          "Best all",
          "Best sel",
          "Trk sel",
          "Has coll",
          "Reas rm"};
        for (int iBin = 0; iBin < static_cast<int>(McTrackStatus::nMcTrackStatus); iBin++) {
          registryMC.get<THnSparse>(HIST("Tracks/Centrality/hMcTrackStatus"))->GetAxis(0)->SetBinLabel(iBin + 1, labelMcStatusCentBest[iBin].data());
        }
        registryMC.add({"Tracks/Centrality/DCA3d", "; p_{T} (GeV/c); #eta; DCA_{XY} (cm); DCA_{Z} (cm); centrality; occupancy", {HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, centralityAxis, occupancyAxis}}});
        registryMC.add({"Tracks/Centrality/TrackAmbDegree", "; N_{coll}^{comp}; centrality; occupancy", {HistType::kTHnSparseF, {{51, -0.5, 50.5}, centralityAxis, occupancyAxis}}});
      }
      registryMC.add({"Tracks/Centrality/EtaZvtx", "; #eta; #it{z}_{vtx} (cm); centrality; occupancy", {HistType::kTHnSparseF, {etaAxis, zAxis, centralityAxis, occupancyAxis}}});
      registryMC.add({"Tracks/Centrality/EtaZvtxGen", "; #eta; #it{z}_{vtx} (cm); centrality; occupancy", {HistType::kTHnSparseF, {etaAxis, zAxis, centralityAxis, occupancyAxis}}});
      registryMC.add({"Tracks/Centrality/EtaZvtxGen_t", "; #eta; #it{z}_{vtx} (cm); centrality; occupancy", {HistType::kTHnSparseF, {etaAxis, zAxis, centralityAxis, occupancyAxis}}});
      registryMC.add({"Tracks/Centrality/PhiEta", "; #varphi; #eta; centrality; occupancy", {HistType::kTHnSparseF, {phiAxis, etaAxis, centralityAxis, occupancyAxis}}});
      registryMC.add({"Tracks/Centrality/PhiEtaGen", "; #varphi; #eta; centrality; occupancy", {HistType::kTHnSparseF, {phiAxis, etaAxis, centralityAxis, occupancyAxis}}});
      registryMC.add({"Tracks/Centrality/PhiEtaGen_t", "; #varphi; #eta; centrality; occupancy", {HistType::kTHnSparseF, {phiAxis, etaAxis, centralityAxis, occupancyAxis}}});
      registryMC.add({"Events/Centrality/NtrkZvtxGen_t", "; N_{trk}; #it{z}_{vtx} (cm); centrality; occupancy", {HistType::kTHnSparseF, {multAxis, zAxis, centralityAxis, occupancyAxis}}});
      registryMC.add({"Events/Centrality/NtrkZvtxGen", "; N_{trk}; #it{z}_{vtx} (cm); centrality; occupancy", {HistType::kTHnSparseF, {multAxis, zAxis, centralityAxis, occupancyAxis}}});
      registryMC.add({"Tracks/Centrality/NclustersEta", "; nClusters; #eta; centrality; occupancy", {HistType::kTHnSparseF, {nclsAxis, etaAxis, centralityAxis, occupancyAxis}}});
      registryMC.add({"Events/Centrality/NotFoundEventZvtx", "; #it{z}_{vtx} (cm); centrality; occupancy", {HistType::kTHnSparseF, {zAxis, centralityAxis, occupancyAxis}}});
      registryMC.add({"Events/Centrality/ZvtxDiff", "; Z_{rec} - Z_{gen} (cm); centrality; occupancy", {HistType::kTHnSparseF, {deltaZAxis, centralityAxis, occupancyAxis}}});
      registryMC.add({"Events/Centrality/hRecZvtxCent", "; #it{z}_{vtx} (cm); centrality; occupancy", {HistType::kTHnSparseF, {zAxis, centralityAxis, occupancyAxis}}});
    }
    if (doprocessEfficiencyInclusive || doprocessEfficiencyBestInclusive) {
      registryMC.add("Events/hMcEffStatus", "Number of events; Cut; occupancy", {HistType::kTH2F, {{static_cast<int>(McEffStatus::nMcEfftStatus), -0.5, +static_cast<int>(McEffStatus::nMcEfftStatus) - 0.5}, occupancyAxis}});
      std::array<std::string_view, static_cast<int>(McEffStatus::nMcEfftStatus)> labelMcEffStatus{
        "All",
        "Selected",
        "Has Mc Coll",
        "Split Vtx Removed"};
      for (int iBin = 0; iBin < static_cast<int>(McEffStatus::nMcEfftStatus); iBin++) {
        registryMC.get<TH2>(HIST("Events/hMcEffStatus"))->GetXaxis()->SetBinLabel(iBin + 1, labelMcEffStatus[iBin].data());
      }
      registryMC.add({"Events/hVtxZGen", "; #it{z}_{vtx} (cm); genEvtType; occupancy", {HistType::kTHnSparseF, {zAxis, genEvtTypeAxis, occupancyAxis}}});
      registryMC.add({"Tracks/hEffGen", "; p_{T} (GeV/c); #varphi; #eta; #it{z}_{vtx} (cm); genEvtType; occupancy", {HistType::kTHnSparseF, {ptAxis, phiAxis, etaAxis, zAxis, genEvtTypeAxis, occupancyAxis}}});
      registryMC.add({"Events/hVtxZRec", "#it{z}_{vtx} (cm); occupancy", {HistType::kTH2F, {zAxis, occupancyAxis}}});
      registryMC.add({"Tracks/hEffRec", "; p_{T} (GeV/c); #varphi; #eta; #it{z}_{vtx} (cm); recEvtType; occupancy", {HistType::kTHnSparseF, {ptAxis, phiAxis, etaAxis, zAxis, recEvtTypeAxis, occupancyAxis}}});
      registryMC.add({"Tracks/hEtaRes", "#eta resolution;;(#eta_{rec} - #eta_{gen})/#eta_{gen}; occupancy", {HistType::kTHnSparseF, {etaAxis, {100, -1.0, 1.0}, occupancyAxis}}});
    }
    if (doprocessEfficiencyCentFT0C || doprocessEfficiencyBestCentFT0C) {
      registryMC.add("Events/Centrality/hMcEffStatus", "Number of events; Cut; centrality; occupancy", {HistType::kTHnSparseF, {{static_cast<int>(McEffStatus::nMcEfftStatus), -0.5, +static_cast<int>(McEffStatus::nMcEfftStatus) - 0.5}, centralityAxis, occupancyAxis}});
      std::array<std::string_view, static_cast<int>(McEffStatus::nMcEfftStatus)> labelMcEffStatusCent{
        "All",
        "Selected",
        "Has Mc Coll",
        "Split Vtx Removed"};
      for (int iBin = 0; iBin < static_cast<int>(McEffStatus::nMcEfftStatus); iBin++) {
        registryMC.get<THnSparse>(HIST("Events/Centrality/hMcEffStatus"))->GetAxis(0)->SetBinLabel(iBin + 1, labelMcEffStatusCent[iBin].data());
      }
      registryMC.add({"Events/Centrality/hVtxZGen", "; #it{z}_{vtx} (cm); genEvtType; centrality; occupancy", {HistType::kTHnSparseF, {zAxis, genEvtTypeAxis, centralityAxis, occupancyAxis}}});
      registryMC.add({"Tracks/Centrality/hEffGen", "; p_{T} (GeV/c); #varphi; #eta; #it{z}_{vtx} (cm); genEvtType; centrality; occupancy", {HistType::kTHnSparseF, {ptAxis, phiAxis, etaAxis, zAxis, genEvtTypeAxis, centralityAxis, occupancyAxis}}});
      registryMC.add({"Events/Centrality/hVtxZRec", "#it{z}_{vtx} (cm); centrality; occupancy", {HistType::kTHnSparseF, {zAxis, centralityAxis, occupancyAxis}}});
      registryMC.add({"Tracks/Centrality/hEffRec", "; p_{T} (GeV/c); #varphi; #eta; #it{z}_{vtx} (cm); recEvtType; centrality; occupancy", {HistType::kTHnSparseF, {ptAxis, phiAxis, etaAxis, zAxis, recEvtTypeAxis, centralityAxis, occupancyAxis}}});
    }
    if (doprocessEfficiencyIdxInlusive || doprocessEfficiencyIdxBestInlusive) {
      registryMC.add({"Tracks/hEffIdxGen", "; p_{T} (GeV/c); #eta; recEvtIdxType; occupancy", {HistType::kTHnSparseF, {ptAxis, etaAxis, recEvtIdxTypeAxis, occupancyAxis}}});
      registryMC.add({"Tracks/hEffIdxRec", "; p_{T} (GeV/c); #eta; genEvtIdxType; occupancy", {HistType::kTHnSparseF, {ptAxis, etaAxis, genEvtIdxTypeAxis, occupancyAxis}}});
      registryMC.add({"Tracks/NmftTrkPerPart", "; #it{N}_{mft tracks per particle}; occupancy", {HistType::kTH2F, {multAxis, occupancyAxis}}});
    }
    if (doprocessEfficiencyIdxCentFT0C || doprocessEfficiencyIdxBestCentFT0C) {
      registryMC.add({"Tracks/Centrality/hEffIdxGen", "; p_{T} (GeV/c); #eta; recEvtIdxType; centrality; occupancy", {HistType::kTHnSparseF, {ptAxis, etaAxis, recEvtIdxTypeAxis, centralityAxis, occupancyAxis}}});
      registryMC.add({"Tracks/Centrality/hEffIdxRec", "; p_{T} (GeV/c); #eta; genEvtIdxType; centrality; occupancy", {HistType::kTHnSparseF, {ptAxis, etaAxis, genEvtIdxTypeAxis, centralityAxis, occupancyAxis}}});
      registryMC.add({"Tracks/Centrality/NmftTrkPerPart", "; #it{N}_{mft tracks per particle}; centrality; occupancy", {HistType::kTH2F, {multAxis, centralityAxis, occupancyAxis}}});
    }
    if (doprocessMcSgnEvtLossCentFT0C) {
      registryMC.add("Events/hEvtMcGen", "Events/hEvtMcGen", {HistType::kTH1F, {{4, 0.f, 4.f}}});
      registryMC.get<TH1>(HIST("Events/hEvtMcGen"))->GetXaxis()->SetBinLabel(1, "all");
      registryMC.get<TH1>(HIST("Events/hEvtMcGen"))->GetXaxis()->SetBinLabel(2, "z-vtx");
      registryMC.get<TH1>(HIST("Events/hEvtMcGen"))->GetXaxis()->SetBinLabel(3, "isInelGt0wMft");
      registryMC.get<TH1>(HIST("Events/hEvtMcGen"))->GetXaxis()->SetBinLabel(4, "TVX");
      //
      registryMC.add("Events/EvtSigLossStatus", ";status;centrality", {HistType::kTH1F, {{3, 0.5, 3.5}}});
      auto hstat = registryMC.get<TH2>(HIST("Events/EvtSigLossStatus"));
      hstat->GetXaxis()->SetBinLabel(1, "All MC gen events");
      hstat->GetXaxis()->SetBinLabel(2, "MC gen events with rec event with event selection");
      hstat->GetXaxis()->SetBinLabel(3, "MC gen events with no rec events");
      //
      registryMC.add({"Events/hNchGen", "Evt loss; Gen Nch FT0C; evtLossType", {HistType::kTH2F, {multFT0cAxis, evtLossTypeAxis}}});
      registryMC.add({"Events/hMultGenVsCentSplit", "Split MC events: Gen Nch vs Rec Cent; rec cent; Gen Nch ", {HistType::kTH2F, {centralityAxis, multFT0cAxis}}});
      registryMC.add({"Events/hNchTVX", "; Nch; status", {HistType::kTH2F, {{2, 0, 2}, multAxis}}});
      registryMC.add({"Events/hMultGenVsCent", "event mult MC gen", {HistType::kTH2F, {centralityAxis, multFT0cAxis}}});
      registryMC.add({"Events/hMultGenVsCentNParticlesEta05", "event mult MC gen", {HistType::kTH2F, {centralityAxis, multAxis}}});
      registryMC.add({"Events/hMultGenVsCentNParticlesEtaMFT", "event mult MC gen", {HistType::kTH2F, {centralityAxis, multAxis}}});
      registryMC.add({"Events/hMultGenVsCentRec", "event mult MC gen vs centrality", {HistType::kTH2F, {centralityAxis, multFT0cAxis}}});
      registryMC.add({"Events/hMultGenVsCentRecNParticlesEta05", "event mult MC gen vs centrality", {HistType::kTH2F, {centralityAxis, multAxis}}});
      registryMC.add({"Events/hMultGenVsCentRecNParticlesEtaMFT", "event mult MC gen vs centrality", {HistType::kTH2F, {centralityAxis, multAxis}}});
      registryMC.add({"Tracks/hEtaVsNchGen", "; #eta; mult gen", {HistType::kTH2F, {etaAxis, multFT0cAxis}}});
      registryMC.add({"Tracks/hEtaVsNchGenRecEvt", "; #eta; mult gen w/ Rec evt", {HistType::kTH2F, {etaAxis, multFT0cAxis}}});
    }
    if (doprocessCorrelationwBestTracksInclusive) {
      registryData.add("Events/hMultMFTvsFT0A", "MultMFT_vs_FT0A", {HistType::kTH2F, {multAxis, multFT0aAxis}});
      registryData.add("Events/hMultMFTvsFT0C", "MultMFT_vs_FT0C", {HistType::kTH2F, {multAxis, multFT0cAxis}});
      registryData.add("Events/hNPVtracksVsFT0C", "NPVtracks_vs_FT0C", {HistType::kTH2F, {pvAxis, multFT0cAxis}});
      registryData.add("Events/hMultMFTvsFV0A", "MultMFT_vs_FV0A", {HistType::kTH2F, {multAxis, multFV0aAxis}});
      registryData.add("Events/hNPVtracksVsMultMFT", "NPVtracks_vs_MultMFT", {HistType::kTH2F, {pvAxis, multAxis}});
    }
  }

  /// Filters - particles
  Filter primaries = ifnode(nsqrt(aod::mcparticle::vx * aod::mcparticle::vx + aod::mcparticle::vy * aod::mcparticle::vy) > CfloatFive, false, ncheckbit(aod::mcparticle_v2::storedFlags, (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary));
  /// Joined tables
  using FullBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;
  using CollBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
  using ExtBCs = soa::Join<aod::BCs, aod::Timestamps, aod::MatchedBCCollisionsSparseMulti>;
  /// Collisions
  using Colls = soa::Join<aod::Collisions, aod::EvSels>;
  using CollsCentFT0C = soa::Join<aod::Collisions, aod::CentFT0Cs, aod::EvSels>;
  using CollsGenCentFT0C = soa::Join<aod::McCollisionLabels, aod::Collisions, aod::CentFT0Cs, aod::EvSels>;
  using CollsCorr = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::PVMults, aod::CentFT0Cs, aod::CentFV0As, aod::CentFT0CVariant1s, aod::CentFT0Ms, aod::CentNGlobals, aod::CentMFTs>;
  using CollsMCExtraMult = soa::Join<aod::McCollisions, aod::MultMCExtras, aod::McCollsExtra>;
  /// Tracks
  using MftTracksLabeled = soa::Join<aod::MFTTracks, aod::McMFTTrackLabels>;
  using MftBestTracksLabeled = soa::Join<aod::MFTTracks, aod::BestCollisionsFwd3d, aod::McMFTTrackLabels>;
  /// Particles
  using ParticlesIdx = soa::Join<aod::McParticles, aod::ParticlesToMftTracks>;
  using FiltParticlesIdx = soa::Filtered<ParticlesIdx>;

  template <typename T>
  float getDCAz(const T& track)
  {
    if constexpr (requires { track.bestDCAZ(); }) {
      return track.bestDCAZ();
    } else {
      return 999.;
    }
  }

  std::unordered_map<int64_t, int64_t> mapMcCollIdPerRecColl;
  template <typename C>
  void buildLookupTable(C const& collisions)
  {
    const auto& nRecoColls = collisions.size();
    mapMcCollIdPerRecColl.clear();
    mapMcCollIdPerRecColl.reserve(nRecoColls);
    int maxNcontributors = -1;
    for (auto const& collision : collisions) {
      // const int nContrib = collision.multPVTotalContributors();
      const int nContrib = collision.numContrib();
      if (maxNcontributors < nContrib) {
        maxNcontributors = nContrib;
      }
      mapMcCollIdPerRecColl.emplace(collision.globalIndex(), collision.mcCollisionId());
    }
  }

  void initCCDB(ExtBCs::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(gConf.grpmagPath, bc.timestamp());
    LOG(info) << "Setting magnetic field to current " << grpmag->getL3Current()
              << " A for run " << bc.runNumber()
              << " from its GRPMagField CCDB object";
    o2::base::Propagator::initFieldFromGRP(grpmag);
    mRunNumber = bc.runNumber();

    auto field = dynamic_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
    bZ = field->getBz(centerMFT.data());
    LOG(info) << "The field at the center of the MFT is bZ = " << bZ;
  }

  template <bool fillHis = true, typename B>
  bool isBestTrackSelected(const B& besttrack)
  {
    if constexpr (fillHis) {
      registryQC.fill(HIST("Tracks/hBestTrkSel"), static_cast<int>(TrkBestSel::trkBestSelAll));
    }
    if (besttrack.bestCollisionId() < CintZero) {
      return false;
    }
    if constexpr (fillHis) {
      registryQC.fill(HIST("Tracks/hBestTrkSel"), static_cast<int>(TrkBestSel::trkBestSelCollID));
    }
    if (besttrack.ambDegree() == CintZero) {
      return false;
    }
    if constexpr (fillHis) {
      registryQC.fill(HIST("Tracks/hBestTrkSel"), static_cast<int>(TrkBestSel::trkBestSelOrphan));
    }
    if (std::abs(besttrack.bestDCAXY()) >= trackCuts.maxDCAxy) {
      return false;
    }
    if constexpr (fillHis) {
      registryQC.fill(HIST("Tracks/hBestTrkSel"), static_cast<int>(TrkBestSel::trkBestSelDCAxyCut));
    }
    if (trackCuts.useDCAzCut) {
      const float bestDcaZ = getDCAz(besttrack);
      if (std::abs(bestDcaZ) >= trackCuts.maxDCAz) {
        return false;
      }
    }
    if constexpr (fillHis) {
      registryQC.fill(HIST("Tracks/hBestTrkSel"), static_cast<int>(TrkBestSel::trkBestSelDCAzCut));
    }
    return true;
  }

  template <bool fillHis = true, typename T>
  bool isTrackSelected(const T& track)
  {
    if constexpr (fillHis) {
      registryQC.fill(HIST("Tracks/hTrkSel"), static_cast<int>(TrkSel::trkSelAll));
    }
    if (track.nClusters() < trackCuts.minNclusterMft) {
      return false;
    }
    if constexpr (fillHis) {
      registryQC.fill(HIST("Tracks/hTrkSel"), static_cast<int>(TrkSel::trkSelNCls));
    }
    if (trackCuts.useChi2Cut) {
      float nclMft = std::max(2.0f * track.nClusters() - 5.0f, 1.0f);
      float mftChi2NCl = track.chi2() / nclMft;
      if (mftChi2NCl > trackCuts.maxChi2NCl) {
        return false;
      }
    }
    if constexpr (fillHis) {
      registryQC.fill(HIST("Tracks/hTrkSel"), static_cast<int>(TrkSel::trkSelChi2Ncl));
    }
    if (track.eta() < trackCuts.minEta || track.eta() > trackCuts.maxEta) {
      return false;
    }
    if constexpr (fillHis) {
      registryQC.fill(HIST("Tracks/hTrkSel"), static_cast<int>(TrkSel::trkSelEta));
    }
    if (trackCuts.usephiCut) {
      float phi = track.phi();
      if ((phi < trackCuts.phiCut) ||
          ((phi > PI - trackCuts.phiCut) && (phi < PI + trackCuts.phiCut)) ||
          (phi > TwoPI - trackCuts.phiCut) ||
          ((phi > ((PIHalf - 0.1) * PI) - trackCuts.phiCut) &&
           (phi < ((PIHalf - 0.1) * PI) + trackCuts.phiCut))) {
        return false;
      }
    }
    if constexpr (fillHis) {
      registryQC.fill(HIST("Tracks/hTrkSel"), static_cast<int>(TrkSel::trkSelPhiCut));
    }
    if (trackCuts.usePtCut && track.pt() < trackCuts.minPt) {
      return false;
    }
    if constexpr (fillHis) {
      registryQC.fill(HIST("Tracks/hTrkSel"), static_cast<int>(TrkSel::trkSelPt));
    }
    if (trackCuts.requireCA && !track.isCA()) {
      return false;
    }
    if constexpr (fillHis) {
      registryQC.fill(HIST("Tracks/hTrkSel"), static_cast<int>(TrkSel::trkSelCA));
    }
    return true;
  }

  template <typename C, bool fillHis = false, typename T>
  int countTracks(T const& tracks, float z, float c, float occ)
  {
    auto nTrk = 0;
    for (auto const& track : tracks) {
      if (fillHis) {
        if constexpr (has_reco_cent<C>) {
          registryData.fill(HIST("Tracks/Centrality/Chi2Eta"), track.chi2(), track.eta(), c, occ);
          registryData.fill(HIST("Tracks/Centrality/NclustersEta"), track.nClusters(), track.eta(), c, occ);
        } else {
          registryData.fill(HIST("Tracks/Chi2Eta"), track.chi2(), track.eta(), occ);
          registryData.fill(HIST("Tracks/NclustersEta"), track.nClusters(), track.eta(), occ);
        }
      }
      if (!isTrackSelected(track)) {
        continue;
      }
      if (fillHis) {
        if constexpr (has_reco_cent<C>) {
          registryData.fill(HIST("Tracks/Centrality/EtaZvtx"), track.eta(), z, c, occ);
          registryData.fill(HIST("Tracks/Centrality/PhiEta"), track.phi(), track.eta(), c, occ);
          registryData.fill(HIST("Tracks/Centrality/TanLambda"), track.tgl(), c, occ);
          registryData.fill(HIST("Tracks/Centrality/InvQPt"), track.signed1Pt(), c, occ);
        } else {
          registryData.fill(HIST("Tracks/EtaZvtx"), track.eta(), z, occ);
          registryData.fill(HIST("Tracks/PhiEta"), track.phi(), track.eta(), occ);
          registryData.fill(HIST("Tracks/TanLambda"), track.tgl(), occ);
          registryData.fill(HIST("Tracks/InvQPt"), track.signed1Pt(), occ);
        }
      }
      ++nTrk;
    }
    return nTrk;
  }

  template <typename C, bool fillHis = false, typename B>
  void countBestTracksExtra(B const& besttracksExtra, float c, float occ)
  {
    for (auto const& etrack : besttracksExtra) {
      if (fillHis) {
        if constexpr (has_reco_cent<C>) {
          if (gConf.cfgUseTrackParExtra) {
            registryData.fill(HIST("Tracks/Centrality/TanLambdaExtra"), etrack.tgl(), c, occ);
            registryData.fill(HIST("Tracks/Centrality/InvQPtExtra"), etrack.signed1Pt(), c, occ);
            registryData.fill(HIST("Tracks/Centrality/EtaExtra"), etrack.etas(), c, occ);
            registryData.fill(HIST("Tracks/Centrality/PhiExtra"), etrack.phis(), c, occ);
          }
        } else {
          if (gConf.cfgUseTrackParExtra) {
            registryData.fill(HIST("Tracks/TanLambdaExtra"), etrack.tgl(), occ);
            registryData.fill(HIST("Tracks/InvQPtExtra"), etrack.signed1Pt(), occ);
            registryData.fill(HIST("Tracks/EtaExtra"), etrack.etas(), occ);
            registryData.fill(HIST("Tracks/PhiExtra"), etrack.phis(), occ);
          }
        }
      }
    }
  }

  template <typename C, bool fillHis = false, typename T, typename B>
  int countBestTracks(T const& tracks, B const& besttracks, float z, float c, float occ)
  {
    auto nATrk = 0;
    ambiguousTrkIds.reserve(besttracks.size());
    reassignedTrkIds.reserve(besttracks.size());
    for (auto const& atrack : besttracks) {
      if (!isBestTrackSelected(atrack)) {
        continue;
      }
      const float bestDcaZ = getDCAz(atrack);
      auto itrack = atrack.template mfttrack_as<T>();
      if constexpr (has_reco_cent<C>) {
        registryData.fill(HIST("Tracks/Centrality/Chi2Eta"), itrack.chi2(), itrack.eta(), c, occ);
        registryData.fill(HIST("Tracks/Centrality/NclustersEta"), itrack.nClusters(), itrack.eta(), c, occ);
      } else {
        registryData.fill(HIST("Tracks/Chi2Eta"), itrack.chi2(), itrack.eta(), occ);
        registryData.fill(HIST("Tracks/NclustersEta"), itrack.nClusters(), itrack.eta(), occ);
      }
      if (!isTrackSelected(itrack)) {
        continue;
      }
      ambiguousTrkIds.emplace_back(atrack.mfttrackId());
      ++nATrk;
      if (fillHis) {
        if constexpr (has_reco_cent<C>) {
          registryData.fill(HIST("Tracks/Centrality/EtaZvtx"), itrack.eta(), z, c, occ);
          registryData.fill(HIST("Tracks/Centrality/PhiEta"), itrack.phi(), itrack.eta(), c, occ);
          registryData.fill(HIST("Tracks/Centrality/TanLambda"), itrack.tgl(), c, occ);
          registryData.fill(HIST("Tracks/Centrality/InvQPt"), itrack.signed1Pt(), c, occ);
          registryData.fill(HIST("Tracks/Centrality/DCA3d"), itrack.pt(), itrack.eta(), atrack.bestDCAXY(), bestDcaZ, c, occ);
          registryData.fill(HIST("Tracks/Centrality/TrackAmbDegree"), atrack.ambDegree(), c, occ);
        } else {
          registryData.fill(HIST("Tracks/EtaZvtx"), itrack.eta(), z, occ);
          registryData.fill(HIST("Tracks/PhiEta"), itrack.phi(), itrack.eta(), occ);
          registryData.fill(HIST("Tracks/TanLambda"), itrack.tgl(), occ);
          registryData.fill(HIST("Tracks/InvQPt"), itrack.signed1Pt(), occ);
          registryData.fill(HIST("Tracks/DCA3d"), itrack.pt(), itrack.eta(), atrack.bestDCAXY(), bestDcaZ, occ);
          registryData.fill(HIST("Tracks/TrackAmbDegree"), atrack.ambDegree(), occ);
        }
      }
      if (itrack.has_collision() && itrack.collisionId() != atrack.bestCollisionId()) {
        reassignedTrkIds.emplace_back(atrack.mfttrackId());
        if (fillHis) {
          registryQC.fill(HIST("Tracks/hBestTrkSel"), static_cast<int>(TrkBestSel::trkBestSelNumReassoc));
          if constexpr (has_reco_cent<C>) {
            registryData.fill(HIST("Tracks/Centrality/ReTracksEtaZvtx"), itrack.eta(), itrack.template collision_as<C>().posZ(), c, occ);
            registryData.fill(HIST("Tracks/Centrality/ReTracksPhiEta"), itrack.phi(), itrack.eta(), c, occ);
          } else {
            registryData.fill(HIST("Tracks/ReTracksEtaZvtx"), itrack.eta(), itrack.template collision_as<C>().posZ(), occ);
            registryData.fill(HIST("Tracks/ReTracksPhiEta"), itrack.phi(), itrack.eta(), occ);
          }
        }
      }
    }

    for (auto const& track : tracks) {
      if (!isTrackSelected(track)) {
        continue;
      }
      if (fillHis) {
        if constexpr (has_reco_cent<C>) {
          registryData.fill(HIST("Tracks/Centrality/OrigTracksEtaZvtx"), track.eta(), z, c, occ);
          registryData.fill(HIST("Tracks/Centrality/OrigTracksPhiEta"), track.phi(), track.eta(), c, occ);
        } else {
          registryData.fill(HIST("Tracks/OrigTracksEtaZvtx"), track.eta(), z, occ);
          registryData.fill(HIST("Tracks/OrigTracksPhiEta"), track.phi(), track.eta(), occ);
        }
      }
      if (std::find(ambiguousTrkIds.begin(), ambiguousTrkIds.end(), track.globalIndex()) != ambiguousTrkIds.end()) {
        continue;
      }
      if (std::find(reassignedTrkIds.begin(), reassignedTrkIds.end(), track.globalIndex()) != reassignedTrkIds.end()) {
        continue;
      }
      // ++nATrk; // use for testing purposes only!
      if (fillHis) {
        if constexpr (has_reco_cent<C>) {
          registryData.fill(HIST("Tracks/Centrality/RestTracksEtaZvtx"), track.eta(), z, c, occ);
          registryData.fill(HIST("Tracks/Centrality/RestTracksPhiEta"), track.phi(), track.eta(), c, occ);
          // registryData.fill(HIST("Tracks/Centrality/EtaZvtx"), track.eta(), z, c, occ);
          // registryData.fill(HIST("Tracks/Centrality/PhiEta"), phi, track.eta(), c, occ);
          // registryData.fill(HIST("Tracks/Centrality/NclustersEta"), track.nClusters(), track.eta(), c, occ);
        } else {
          registryData.fill(HIST("Tracks/RestTracksEtaZvtx"), track.eta(), z, occ);
          registryData.fill(HIST("Tracks/RestTracksPhiEta"), track.phi(), track.eta(), occ);
          // registryData.fill(HIST("Tracks/EtaZvtx"), track.eta(), z, occ);
          // registryData.fill(HIST("Tracks/PhiEta"), phi, track.eta(), occ);
          // registryData.fill(HIST("Tracks/NclustersEta"), track.nClusters(), track.eta(), occ);
        }
      }
    }
    ambiguousTrkIds.clear();
    ambiguousTrkIds.shrink_to_fit();
    reassignedTrkIds.clear();
    reassignedTrkIds.shrink_to_fit();
    return nATrk;
  }

  template <typename P>
  bool isInelGt0wMft(P const& particles)
  {
    int nChrgMc = 0;
    int nChrgFT0A = 0;
    int nChrgFT0C = 0;
    for (auto const& particle : particles) {
      if (!isChrgParticle(particle.pdgCode())) {
        continue;
      }
      if (!particle.isPhysicalPrimary()) {
        continue;
      }
      // trigger TVX
      if (particle.eta() > CminAccFT0A && particle.eta() < CmaxAccFT0A) {
        nChrgFT0A++;
      }
      if (particle.eta() > CminAccFT0C && particle.eta() < CmaxAccFT0C) {
        nChrgFT0C++;
      }
      // acceptance MFT
      if (particle.eta() < trackCuts.minEta || particle.eta() > trackCuts.maxEta) {
        continue;
      }
      nChrgMc++;
    }
    if (nChrgFT0A == CintZero || nChrgFT0C == CintZero) {
      registryMC.fill(HIST("Events/hNchTVX"), nChrgMc, 0.5);
      return false;
    }
    registryMC.fill(HIST("Events/hNchTVX"), nChrgMc, 1.5);

    return nChrgMc != CintZero;
  }

  template <typename P>
  int countPart(P const& particles)
  {
    auto nCharged = 0;
    for (auto const& particle : particles) {
      if (!isChrgParticle(particle.pdgCode())) {
        continue;
      }
      if (gConf.cfgUseParticleSel && !isParticleSelected(particle)) {
        continue;
      }
      if (particle.eta() < trackCuts.minEta || particle.eta() > trackCuts.maxEta) {
        continue;
      }
      nCharged++;
    }
    return nCharged;
  }

  template <typename P>
  bool isParticleSelected(P const& particle)
  {
    if (gConf.cfgUsePrimaries && !particle.isPhysicalPrimary()) {
      return false;
    }
    if (!gConf.cfgUsePrimaries && (gConf.cfgUseSecondaries && particle.isPhysicalPrimary())) {
      return false;
    }
    if (particle.eta() < trackCuts.minEta || particle.eta() > trackCuts.maxEta) {
      return false;
    }
    if (trackCuts.usephiCut) {
      float phi = particle.phi();
      if ((phi < trackCuts.phiCut) ||
          ((phi > PI - trackCuts.phiCut) && (phi < PI + trackCuts.phiCut)) ||
          (phi > TwoPI - trackCuts.phiCut) ||
          ((phi > ((PIHalf - 0.1) * PI) - trackCuts.phiCut) &&
           (phi < ((PIHalf - 0.1) * PI) + trackCuts.phiCut))) {
        return false;
      }
    }
    return true;
  }

  template <typename P>
  bool cutInPhi(P const& particle)
  {
    if (!trackCuts.usephiCut) {
      return true;
    }
    auto phi = particle.phi();
    return (phi < trackCuts.phiCut) || ((phi > PI - trackCuts.phiCut) && (phi < PI + trackCuts.phiCut)) || (phi > TwoPI - trackCuts.phiCut) || ((phi > ((PIHalf - 0.1) * PI) - trackCuts.phiCut) && (phi < ((PIHalf - 0.1) * PI) + trackCuts.phiCut));
  }

  template <typename C>
  float getOccupancy(C const& collision, uint occEstimator)
  {
    switch (occEstimator) {
      case OccupancyEst::TrkITS:
        return collision.trackOccupancyInTimeRange();
      case OccupancyEst::Ft0C:
        return collision.ft0cOccupancyInTimeRange();
      default:
        LOG(fatal) << "No valid occupancy estimator ";
        break;
    }
    return -1.f;
  }

  void initHadronicRate(CollBCs::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    mRunNumber = bc.runNumber();
    if (!gHadronicRate.contains(mRunNumber)) {
      auto runDuration = ccdb->getRunDuration(mRunNumber);
      mSOR = runDuration.first;
      mMinSeconds = std::floor(mSOR * 1.e-3);               /// round tsSOR to the highest integer lower than tsSOR
      float maxSec = std::ceil(runDuration.second * 1.e-3); /// round tsEOR to the lowest integer higher than tsEOR
      const AxisSpec axisSeconds{static_cast<int>((maxSec - mMinSeconds) / 20.f), 0, maxSec - mMinSeconds, "Seconds since SOR"};
      int hadronicRateBins = static_cast<int>(eventCuts.maxIR - eventCuts.minIR);
      gHadronicRate[mRunNumber] = registryMC.add<TH2>(Form("HadronicRate/%i", mRunNumber), ";Time since SOR (s);Hadronic rate (kHz)", kTH2D, {axisSeconds, {hadronicRateBins, eventCuts.minIR, eventCuts.maxIR}}).get();
    }
    gCurrentHadronicRate = gHadronicRate[mRunNumber];
  }

  template <bool fillHis = false, typename C>
  bool isGoodEvent(C const& collision)
  {
    if constexpr (fillHis) {
      registryQC.fill(HIST("Events/hEvtSel"), static_cast<int>(EvtSel::evtAll));
    }
    if (!collision.sel8()) {
      return false;
    }
    if constexpr (fillHis) {
      registryQC.fill(HIST("Events/hEvtSel"), static_cast<int>(EvtSel::evtSel));
    }
    if (eventCuts.requireIsGoodZvtxFT0VsPV && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if constexpr (fillHis) {
      registryQC.fill(HIST("Events/hEvtSel"), static_cast<int>(EvtSel::evtIsGoodZvtx));
    }
    if (eventCuts.requireRejectSameBunchPileup && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if constexpr (fillHis) {
      registryQC.fill(HIST("Events/hEvtSel"), static_cast<int>(EvtSel::evtNoSameBunchPileup));
    }
    if (collision.posZ() <= eventCuts.minZvtx || collision.posZ() >= eventCuts.maxZvtx) {
      return false;
    }
    if constexpr (fillHis) {
      registryQC.fill(HIST("Events/hEvtSel"), static_cast<int>(EvtSel::evtZvtxCut));
    }
    if (eventCuts.requireNoCollInTimeRangeStd &&
        !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return false;
    }
    if constexpr (fillHis) {
      registryQC.fill(HIST("Events/hEvtSel"), static_cast<int>(EvtSel::evtNoCollInTimeRangeStd));
    }
    if (eventCuts.requireNoCollInTimeRangeNarrow &&
        !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow)) {
      return false;
    }
    if constexpr (fillHis) {
      registryQC.fill(HIST("Events/hEvtSel"), static_cast<int>(EvtSel::evtNoCollInTimeRangeNarrow));
    }
    if (eventCuts.requireNoCollInTimeRangeStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict)) {
      return false;
    }
    if constexpr (fillHis) {
      registryQC.fill(HIST("Events/hEvtSel"), static_cast<int>(EvtSel::evtNoCollInTimeRangeStrict));
    }
    if (eventCuts.requireNoCollInRofStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStrict)) {
      return false;
    }
    if constexpr (fillHis) {
      registryQC.fill(HIST("Events/hEvtSel"), static_cast<int>(EvtSel::evtNoCollInRofStrict));
    }
    if (eventCuts.requireNoCollInRofStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      return false;
    }
    if constexpr (fillHis) {
      registryQC.fill(HIST("Events/hEvtSel"), static_cast<int>(EvtSel::evtNoCollInRofStandard));
    }
    if (eventCuts.requireNoHighMultCollInPrevRof && !collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
      return false;
    }
    if constexpr (fillHis) {
      registryQC.fill(HIST("Events/hEvtSel"), static_cast<int>(EvtSel::evtNoHighMultCollInPrevRof));
    }
    if (eventCuts.requireGoodITSLayersAll && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      return false;
    }
    if constexpr (fillHis) {
      registryQC.fill(HIST("Events/hEvtSel"), static_cast<int>(EvtSel::evtGoodITSLayersAll));
    }
    if (eventCuts.minOccupancy >= 0 &&
        getOccupancy(collision, eventCuts.occupancyEstimator) <
          eventCuts.minOccupancy) {
      return false;
    }
    if constexpr (fillHis) {
      registryQC.fill(HIST("Events/hEvtSel"), static_cast<int>(EvtSel::evtBelowMinOccup));
    }
    if (eventCuts.maxOccupancy >= 0 &&
        getOccupancy(collision, eventCuts.occupancyEstimator) >
          eventCuts.maxOccupancy) {
      return false;
    }
    if constexpr (fillHis) {
      registryQC.fill(HIST("Events/hEvtSel"), static_cast<int>(EvtSel::evtAboveMaxOccup));
    }
    if (rctCuts.requireRCTFlagChecker && !rctChecker(collision)) {
      return false;
    }
    if constexpr (fillHis) {
      registryQC.fill(HIST("Events/hEvtSel"), static_cast<int>(EvtSel::evtRCTFlagChecker));
    }
    if (rctCuts.requireRCTFlagCheckerExtra && !rctCheckerExtra(collision)) {
      return false;
    }
    if constexpr (fillHis) {
      registryQC.fill(HIST("Events/hEvtSel"), static_cast<int>(EvtSel::evtRCTFlagCheckerExtra));
    }
    return true;
  }

  /// @brief Selection of charged particles
  /// @return true: charged; false: not charged
  bool isChrgParticle(int code)
  {
    auto p = pdg->GetParticle(code);
    auto charge = 0.;
    if (p != nullptr) {
      charge = p->Charge();
    }
    return std::abs(charge) >= CminCharge;
  }

  template <bool isCent, typename P>
  void fillHistMC(P const& particles, float zvtx, float c, float occ)
  {
    for (auto const& particle : particles) {
      if (!isChrgParticle(particle.pdgCode())) {
        continue;
      }
      if (gConf.cfgUseParticleSel && !isParticleSelected(particle)) {
        continue;
      }
      if constexpr (isCent) {
        registryMC.fill(HIST("Tracks/Centrality/EtaZvtxGen"), particle.eta(), zvtx, c, occ);
        registryMC.fill(HIST("Tracks/Centrality/PhiEtaGen"), particle.phi(), particle.eta(), c, occ);
      } else {
        registryMC.fill(HIST("Tracks/EtaZvtxGen"), particle.eta(), zvtx, occ);
        registryMC.fill(HIST("Tracks/PhiEtaGen"), particle.phi(), particle.eta(), occ);
      }
    }
  }

  /// @brief process function for general event statistics
  void processTagging(FullBCs const& bcs, CollsCentFT0C const& collisions)
  {
    std::vector<typename std::decay_t<decltype(collisions)>::iterator> cols;
    for (auto const& bc : bcs) {
      if (static_cast<int>(bc.selection_bit(aod::evsel::kIsBBT0A) &&
                           bc.selection_bit(aod::evsel::kIsBBT0C)) != 0) {
        registryQC.fill(HIST("hBcSel"), 0);
        cols.clear();
        for (auto const& collision : collisions) {
          if (collision.has_foundBC()) {
            if (collision.foundBCId() == bc.globalIndex()) {
              cols.emplace_back(collision);
            }
          } else if (collision.bcId() == bc.globalIndex()) {
            cols.emplace_back(collision);
          }
        }
        LOGP(debug, "BC {} has {} collisions", bc.globalBC(), cols.size());
        if (!cols.empty()) {
          registryQC.fill(HIST("hBcSel"), 1);
          if (cols.size() > 1) {
            registryQC.fill(HIST("hBcSel"), 2);
          }
        }
      }
    }
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processTagging, "Collect event sample stats", true);

  /// @brief process function for counting tracks
  template <typename C>
  void processData(typename C::iterator const& collision,
                   aod::MFTTracks const& tracks,
                   CollBCs const& /*bcs*/)
  {
    float c = getRecoCent(collision);
    float occ = getOccupancy(collision, eventCuts.occupancyEstimator);
    auto bc = collision.template foundBC_as<CollBCs>();
    if constexpr (has_reco_cent<C>) {
      registryData.fill(HIST("Events/Centrality/Selection"), 1., c, occ);
    } else {
      registryData.fill(HIST("Events/Selection"), 1., occ);
    }
    if (gConf.cfgDoIR) {
      initHadronicRate(bc);
      float ir = !gConf.cfgIRSource.value.empty() ? rateFetcher.fetch(ccdb.service, bc.timestamp(), bc.runNumber(), gConf.cfgIRSource, gConf.cfgIRCrashOnNull) * 1.e-3 : -1;
      if constexpr (has_reco_cent<C>) {
        registryData.fill(HIST("Events/Centrality/hInteractionRate"), ir, c, occ);
      } else {
        registryData.fill(HIST("Events/hInteractionRate"), ir, occ);
      }
      float seconds = bc.timestamp() * 1.e-3 - mMinSeconds;
      if (gConf.cfgUseIRCut && (ir < eventCuts.minIR || ir > eventCuts.maxIR)) { // cut on hadronic rate
        return;
      }
      gCurrentHadronicRate->Fill(seconds, ir);
    }
    if (!isGoodEvent<true>(collision)) {
      return;
    }
    auto z = collision.posZ();
    if constexpr (has_reco_cent<C>) {
      registryData.fill(HIST("Events/Centrality/Selection"), 2., c, occ);
      registryData.fill(HIST("Events/Centrality/hZvtxCent"), z, c, occ);
    } else {
      registryData.fill(HIST("Events/Selection"), 2., occ);
    }

    auto nTrk = countTracks<C, true>(tracks, z, c, occ);

    if constexpr (has_reco_cent<C>) {
      registryData.fill(HIST("Events/Centrality/NtrkZvtx"), nTrk, z, c, occ);
    } else {
      registryData.fill(HIST("Events/NtrkZvtx"), nTrk, z, occ);
    }
  }

  /// @brief process function for counting tracks (based on BestCollisionsFwd3d table)
  template <typename C>
  void processDatawBestTracks(typename C::iterator const& collision,
                              aod::MFTTracks const& tracks,
                              soa::SmallGroups<aod::BestCollisionsFwd3d> const& besttracks,
                              aod::BestCollisionsFwd3dExtra const& besttracksExtra,
                              CollBCs const& /*bcs*/)
  {
    float c = getRecoCent(collision);
    float occ = getOccupancy(collision, eventCuts.occupancyEstimator);
    auto bc = collision.template foundBC_as<CollBCs>();
    if constexpr (has_reco_cent<C>) {
      registryData.fill(HIST("Events/Centrality/Selection"), 1., c, occ);
    } else {
      registryData.fill(HIST("Events/Selection"), 1., occ);
    }
    if (gConf.cfgDoIR) {
      initHadronicRate(bc);
      float ir = !gConf.cfgIRSource.value.empty() ? rateFetcher.fetch(ccdb.service, bc.timestamp(), bc.runNumber(), gConf.cfgIRSource, gConf.cfgIRCrashOnNull) * 1.e-3 : -1;
      if constexpr (has_reco_cent<C>) {
        registryData.fill(HIST("Events/Centrality/hInteractionRate"), ir, c, occ);
      } else {
        registryData.fill(HIST("Events/hInteractionRate"), ir, occ);
      }
      float seconds = bc.timestamp() * 1.e-3 - mMinSeconds;
      if (gConf.cfgUseIRCut && (ir < eventCuts.minIR || ir > eventCuts.maxIR)) { // cut on hadronic rate
        return;
      }
      gCurrentHadronicRate->Fill(seconds, ir);
    }
    if (!isGoodEvent<true>(collision)) {
      return;
    }
    auto z = collision.posZ();
    if constexpr (has_reco_cent<C>) {
      registryData.fill(HIST("Events/Centrality/Selection"), 2., c, occ);
      registryData.fill(HIST("Events/Centrality/hZvtxCent"), z, c, occ);
    } else {
      registryData.fill(HIST("Events/Selection"), 2., occ);
      registryData.fill(HIST("Events/hZvtxCent"), z, c, occ);
    }

    auto nBestTrks = countBestTracks<C, true>(tracks, besttracks, z, c, occ);
    countBestTracksExtra<C, true>(besttracksExtra, c, occ);

    if constexpr (has_reco_cent<C>) {
      registryData.fill(HIST("Events/Centrality/NtrkZvtx"), nBestTrks, z, c, occ);
    } else {
      registryData.fill(HIST("Events/NtrkZvtx"), nBestTrks, z, occ);
    }
  }

  void processDataInclusive(Colls::iterator const& collision, aod::MFTTracks const& tracks, CollBCs const& bcs)
  {
    processData<Colls>(collision, tracks, bcs);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processDataInclusive, "Count tracks (inclusive)", false);

  void processDataCentFT0C(CollsCentFT0C::iterator const& collision, aod::MFTTracks const& tracks, CollBCs const& bcs)
  {
    processData<CollsCentFT0C>(collision, tracks, bcs);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processDataCentFT0C, "Count tracks in FT0C centrality bins", false);

  void processDatawBestTracksInclusive(Colls::iterator const& collision, aod::MFTTracks const& tracks, soa::SmallGroups<aod::BestCollisionsFwd3d> const& besttracks, aod::BestCollisionsFwd3dExtra const& besttracksExtra, CollBCs const& bcs)
  {
    processDatawBestTracks<Colls>(collision, tracks, besttracks, besttracksExtra, bcs);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processDatawBestTracksInclusive, "Count tracks based on BestCollisionsFwd3d table (inclusive)", false);

  void processDatawBestTracksCentFT0C(CollsCentFT0C::iterator const& collision, aod::MFTTracks const& tracks, soa::SmallGroups<aod::BestCollisionsFwd3d> const& besttracks, aod::BestCollisionsFwd3dExtra const& besttracksExtra, CollBCs const& bcs)
  {
    processDatawBestTracks<CollsCentFT0C>(collision, tracks, besttracks, besttracksExtra, bcs);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processDatawBestTracksCentFT0C, "Count tracks in FT0C centrality bins based on BestCollisionsFwd3d table", false);

  Preslice<MftTracksLabeled> perCol = o2::aod::fwdtrack::collisionId;
  Partition<aod::McParticles> mcSample = (aod::mcparticle::eta < trackCuts.maxEta) && (aod::mcparticle::eta > trackCuts.minEta);

  template <typename MC, typename C>
  void processMC(typename MC::iterator const& mcCollision,
                 soa::SmallGroups<soa::Join<C, aod::McCollisionLabels>> const& collisions,
                 aod::McParticles const& particles,
                 MftTracksLabeled const& tracks)
  {
    LOGP(debug, "MC col {} has {} reco cols", mcCollision.globalIndex(), collisions.size());
    float occGen = -1.;
    for (const auto& collision : collisions) {
      if (isGoodEvent<false>(collision)) {
        float o = getOccupancy(collision, eventCuts.occupancyEstimator);
        if (o > occGen) {
          occGen = o;
        }
      }
    }
    float cGen = -1;
    if (eventCuts.useGenMult) {
      cGen = mcCollision.multMCFT0C();
    }
    if constexpr (has_reco_cent<C>) {
      float crecMin = 999.;
      for (const auto& collision : collisions) {
        if (isGoodEvent<false>(collision)) {
          float c = getRecoCent(collision);
          if (c < crecMin) {
            crecMin = c;
          }
        }
      }
      if (cGen < 0) {
        cGen = crecMin;
      }
    }

    for (auto const& collision : collisions) {
      float occRec = getOccupancy(collision, eventCuts.occupancyEstimator);
      float cRec = getRecoCent(collision);
      if constexpr (has_reco_cent<C>) {
        registryMC.fill(HIST("Events/Centrality/McStatus"), static_cast<float>(McStatus::kMcRecAll), cRec, occRec);
      } else {
        registryMC.fill(HIST("Events/McStatus"), static_cast<float>(McStatus::kMcRecAll), occRec);
      }
      if (!isGoodEvent<true>(collision)) {
        continue;
      }
      if constexpr (has_reco_cent<C>) {
        registryMC.fill(HIST("Events/Centrality/McStatus"), static_cast<float>(McStatus::kMcRecSel), cRec, occRec);
      } else {
        registryMC.fill(HIST("Events/McStatus"), static_cast<float>(McStatus::kMcRecSel), occRec);
      }
      if (!collision.has_mcCollision()) {
        continue;
      }
      if constexpr (has_reco_cent<C>) {
        registryMC.fill(HIST("Events/Centrality/McStatus"), static_cast<float>(McStatus::kMcRecHasMcColl), cRec, occRec);
      } else {
        registryMC.fill(HIST("Events/McStatus"), static_cast<float>(McStatus::kMcRecHasMcColl), occRec);
      }
      if (gConf.cfgRemoveSplitVertex && collision.globalIndex() != mcCollision.bestCollisionIndex()) {
        continue;
      }
      if constexpr (has_reco_cent<C>) {
        registryMC.fill(HIST("Events/Centrality/McStatus"), static_cast<float>(McStatus::kMcRecNoSplitVtx), cRec, occRec);
        registryMC.fill(HIST("Events/Centrality/hRecZvtxCent"), collision.posZ(), cRec, occRec);
      } else {
        registryMC.fill(HIST("Events/McStatus"), static_cast<float>(McStatus::kMcRecNoSplitVtx), occRec);
      }

      auto nTrkRec = 0;
      auto perColSample = tracks.sliceBy(perCol, collision.globalIndex());
      for (auto const& track : perColSample) {
        if (!isTrackSelected<true>(track)) {
          continue;
        }
        if (track.has_mcParticle() && track.mcMask() == trackCuts.selMcMask) {
          const auto& particle = track.template mcParticle_as<aod::McParticles>();
          if (!isChrgParticle(particle.pdgCode())) {
            continue;
          }
          if (gConf.cfgUseParticleSel && !isParticleSelected(particle)) {
            continue;
          }
          if (collision.mcCollisionId() != particle.mcCollisionId()) {
            continue;
          }
          if constexpr (has_reco_cent<C>) {
            registryMC.fill(HIST("Tracks/Centrality/EtaZvtx"), track.eta(), collision.posZ(), cRec, occRec);
            registryMC.fill(HIST("Tracks/Centrality/PhiEta"), track.phi(), track.eta(), cRec, occRec);
            registryMC.fill(HIST("Tracks/Centrality/NclustersEta"), track.nClusters(), track.eta(), cRec, occRec);
          } else {
            registryMC.fill(HIST("Tracks/EtaZvtx"), track.eta(), collision.posZ(), occRec);
            registryMC.fill(HIST("Tracks/PhiEta"), track.phi(), track.eta(), occRec);
            registryMC.fill(HIST("Tracks/NclustersEta"), track.nClusters(), track.eta(), occRec);
          }
          ++nTrkRec;
        }
      }
      if (eventCuts.useZDiffCut) {
        if (std::abs(collision.posZ() - mcCollision.posZ()) > eventCuts.maxZvtxDiff) {
          continue;
        }
      }
      if constexpr (has_reco_cent<C>) {
        registryMC.fill(HIST("Events/Centrality/NtrkZvtxGen"), nTrkRec, collision.posZ(), cRec, occRec);
        registryMC.fill(HIST("Events/Centrality/ZvtxDiff"), collision.posZ() - mcCollision.posZ(), cRec, occRec);
      } else {
        registryMC.fill(HIST("Events/NtrkZvtxGen"), nTrkRec, collision.posZ(), occRec);
        registryMC.fill(HIST("Events/ZvtxDiff"), collision.posZ() - mcCollision.posZ(), occRec);
      }
    }

    if constexpr (has_reco_cent<C>) {
      registryMC.fill(HIST("Events/Centrality/McStatus"), static_cast<float>(McStatus::kMcGenAll), cGen, occGen);
    } else {
      registryMC.fill(HIST("Events/McStatus"), static_cast<float>(McStatus::kMcGenAll), occGen);
    }

    auto nCharged = countPart(particles);
    if constexpr (has_reco_cent<C>) {
      registryMC.fill(HIST("Events/Centrality/NtrkZvtxGen_t"), nCharged, mcCollision.posZ(), cGen, occGen);
    } else {
      registryMC.fill(HIST("Events/NtrkZvtxGen_t"), nCharged, mcCollision.posZ(), occGen);
    }
    if (collisions.size() == 0) {
      if constexpr (has_reco_cent<C>) {
        registryMC.fill(HIST("Events/Centrality/NotFoundEventZvtx"), mcCollision.posZ(), cGen, occGen);
      } else {
        registryMC.fill(HIST("Events/NotFoundEventZvtx"), mcCollision.posZ(), occGen);
      }
    }
    fillHistMC<has_reco_cent<C>>(particles, mcCollision.posZ(), cGen, occGen);
  }

  void processMCInclusive(CollsMCExtraMult::iterator const& mccollision,
                          soa::SmallGroups<soa::Join<Colls, aod::McCollisionLabels>> const& collisions,
                          aod::McParticles const& particles,
                          MftTracksLabeled const& tracks)
  {
    processMC<CollsMCExtraMult, Colls>(mccollision, collisions, particles, tracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processMCInclusive, "Count MC particles (inclusive)", false);

  void processMCCentFT0C(CollsMCExtraMult::iterator const& mccollision,
                         soa::SmallGroups<soa::Join<CollsCentFT0C, aod::McCollisionLabels>> const& collisions,
                         aod::McParticles const& particles,
                         MftTracksLabeled const& tracks)
  {
    processMC<CollsMCExtraMult, CollsCentFT0C>(mccollision, collisions, particles, tracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processMCCentFT0C, "Count MC particles in FT0C centrality bins", false);

  PresliceUnsorted<aod::BestCollisionsFwd3d> perColU = aod::fwdtrack::bestCollisionId;

  template <typename MC, typename C>
  void processMCbest(typename MC::iterator const& mcCollision,
                     soa::SmallGroups<soa::Join<C, aod::McCollisionLabels>> const& collisions,
                     aod::McParticles const& particles,
                     MftTracksLabeled const& /*tracks*/,
                     MftBestTracksLabeled const& besttracks)
  {
    LOGP(debug, "MC col {} has {} reco cols", mcCollision.globalIndex(), collisions.size());
    float occGen = -1.;
    for (const auto& collision : collisions) {
      if (isGoodEvent<false>(collision)) {
        float o = getOccupancy(collision, eventCuts.occupancyEstimator);
        if (o > occGen) {
          occGen = o;
        }
      }
    }
    float cGen = -1;
    if (eventCuts.useGenMult) {
      cGen = mcCollision.multMCFT0C();
    }
    if constexpr (has_reco_cent<C>) {
      float crecMin = 999.;
      for (const auto& collision : collisions) {
        if (isGoodEvent<false>(collision)) {
          float c = getRecoCent(collision);
          if (c < crecMin) {
            crecMin = c;
          }
        }
      }
      if (cGen < 0) {
        cGen = crecMin;
      }
    }

    buildLookupTable(collisions);

    for (auto const& collision : collisions) {
      float occRec = getOccupancy(collision, eventCuts.occupancyEstimator);
      float cRec = getRecoCent(collision);
      if constexpr (has_reco_cent<C>) {
        registryMC.fill(HIST("Events/Centrality/McStatus"), static_cast<float>(McStatus::kMcRecAll), cRec, occRec);
      } else {
        registryMC.fill(HIST("Events/McStatus"), static_cast<float>(McStatus::kMcRecAll), occRec);
      }
      if (!isGoodEvent<true>(collision)) {
        continue;
      }
      if constexpr (has_reco_cent<C>) {
        registryMC.fill(HIST("Events/Centrality/McStatus"), static_cast<float>(McStatus::kMcRecSel), cRec, occRec);
      } else {
        registryMC.fill(HIST("Events/McStatus"), static_cast<float>(McStatus::kMcRecSel), occRec);
      }
      if (!collision.has_mcCollision()) {
        continue;
      }
      if constexpr (has_reco_cent<C>) {
        registryMC.fill(HIST("Events/Centrality/McStatus"), static_cast<float>(McStatus::kMcRecHasMcColl), cRec, occRec);
      } else {
        registryMC.fill(HIST("Events/McStatus"), static_cast<float>(McStatus::kMcRecHasMcColl), occRec);
      }
      if (gConf.cfgRemoveSplitVertex && collision.globalIndex() != mcCollision.bestCollisionIndex()) {
        continue;
      }
      if constexpr (has_reco_cent<C>) {
        registryMC.fill(HIST("Events/Centrality/McStatus"), static_cast<float>(McStatus::kMcRecNoSplitVtx), cRec, occRec);
        registryMC.fill(HIST("Events/Centrality/hRecZvtxCent"), collision.posZ(), cRec, occRec);
      } else {
        registryMC.fill(HIST("Events/McStatus"), static_cast<float>(McStatus::kMcRecNoSplitVtx), occRec);
      }

      auto nATrk = 0;
      auto perCollisionASample = besttracks.sliceBy(perColU, collision.globalIndex());
      for (auto const& atrack : perCollisionASample) {
        if constexpr (has_reco_cent<C>) {
          registryMC.fill(HIST("Tracks/Centrality/hMcTrackStatus"), static_cast<int>(McTrackStatus::kBestTrkAll), cRec, occRec);
        } else {
          registryMC.fill(HIST("Tracks/hMcTrackStatus"), static_cast<int>(McTrackStatus::kBestTrkAll), occRec);
        }
        if (!isBestTrackSelected(atrack)) {
          continue;
        }
        if constexpr (has_reco_cent<C>) {
          registryMC.fill(HIST("Tracks/Centrality/hMcTrackStatus"), static_cast<int>(McTrackStatus::kBesTrktSel), cRec, occRec);
        } else {
          registryMC.fill(HIST("Tracks/hMcTrackStatus"), static_cast<int>(McTrackStatus::kBesTrktSel), occRec);
        }
        const float bestDcaZ = getDCAz(atrack);
        auto itrack = atrack.template mfttrack_as<MftBestTracksLabeled>();
        if (!isTrackSelected(itrack)) {
          continue;
        }
        if constexpr (has_reco_cent<C>) {
          registryMC.fill(HIST("Tracks/Centrality/hMcTrackStatus"), static_cast<int>(McTrackStatus::kTrkAsBestSel), cRec, occRec);
        } else {
          registryMC.fill(HIST("Tracks/hMcTrackStatus"), static_cast<int>(McTrackStatus::kTrkAsBestSel), occRec);
        }
        if (!itrack.has_collision()) {
          continue;
        }
        if constexpr (has_reco_cent<C>) {
          registryMC.fill(HIST("Tracks/Centrality/hMcTrackStatus"), static_cast<int>(McTrackStatus::kTrkHasColl), cRec, occRec);
        } else {
          registryMC.fill(HIST("Tracks/hMcTrackStatus"), static_cast<int>(McTrackStatus::kTrkHasColl), occRec);
        }
        if (gConf.cfgRemoveReassigned) {
          if (itrack.collisionId() != atrack.bestCollisionId()) {
            continue;
          }
        }
        if constexpr (has_reco_cent<C>) {
          registryMC.fill(HIST("Tracks/Centrality/hMcTrackStatus"), static_cast<int>(McTrackStatus::kTrkReassignedRemoved), cRec, occRec);
        } else {
          registryMC.fill(HIST("Tracks/hMcTrackStatus"), static_cast<int>(McTrackStatus::kTrkReassignedRemoved), occRec);
        }
        if (itrack.collisionId() >= 0 && itrack.has_mcParticle() && itrack.mcMask() == trackCuts.selMcMask) {
          const auto& particle = itrack.template mcParticle_as<aod::McParticles>();
          if (!isChrgParticle(particle.pdgCode())) {
            continue;
          }
          if (gConf.cfgUseParticleSel && !isParticleSelected(particle)) {
            continue;
          }
          if (eventCuts.useZDiffCut) {
            if (std::abs(collision.posZ() - atrack.mcParticle().mcCollision().posZ()) > eventCuts.maxZvtxDiff) {
              continue;
            }
          }
          // if (collision.mcCollisionId() != particle.mcCollisionId()) {
          //   continue;
          // }
          const int bestRecColl = atrack.bestCollisionId();
          if (!mapMcCollIdPerRecColl.contains(bestRecColl)) {
            continue;
          }
          int64_t mcCollIdRec = mapMcCollIdPerRecColl.find(bestRecColl)->second;
          if (mcCollIdRec != particle.mcCollisionId()) {
            continue;
          }
          if constexpr (has_reco_cent<C>) {
            registryMC.fill(HIST("Tracks/Centrality/EtaZvtx"), itrack.eta(), collision.posZ(), cRec, occRec);
            registryMC.fill(HIST("Tracks/Centrality/PhiEta"), itrack.phi(), itrack.eta(), cRec, occRec);
            registryMC.fill(HIST("Tracks/Centrality/NclustersEta"), itrack.nClusters(), itrack.eta(), cRec, occRec);
            registryMC.fill(HIST("Tracks/Centrality/DCA3d"), itrack.pt(), itrack.eta(), atrack.bestDCAXY(), bestDcaZ, cRec, occRec);
            registryMC.fill(HIST("Tracks/Centrality/TrackAmbDegree"), atrack.ambDegree(), cRec, occRec);
          } else {
            registryMC.fill(HIST("Tracks/EtaZvtx"), itrack.eta(), collision.posZ(), occRec);
            registryMC.fill(HIST("Tracks/PhiEta"), itrack.phi(), itrack.eta(), occRec);
            registryMC.fill(HIST("Tracks/NclustersEta"), itrack.nClusters(), itrack.eta(), occRec);
            registryMC.fill(HIST("Tracks/DCA3d"), itrack.pt(), itrack.eta(), atrack.bestDCAXY(), bestDcaZ, occRec);
            registryMC.fill(HIST("Tracks/TrackAmbDegree"), atrack.ambDegree(), occRec);
          }
          ++nATrk;
        }
      }
      if (eventCuts.useZDiffCut) {
        if (std::abs(collision.posZ() - mcCollision.posZ()) > eventCuts.maxZvtxDiff) {
          continue;
        }
      }
      if constexpr (has_reco_cent<C>) {
        registryMC.fill(HIST("Events/Centrality/NtrkZvtxGen"), nATrk, collision.posZ(), cRec, occRec);
        registryMC.fill(HIST("Events/Centrality/ZvtxDiff"), collision.posZ() - mcCollision.posZ(), cRec, occRec);
      } else {
        registryMC.fill(HIST("Events/NtrkZvtxGen"), nATrk, collision.posZ(), occRec);
        registryMC.fill(HIST("Events/ZvtxDiff"), collision.posZ() - mcCollision.posZ(), occRec);
      }
    }

    if constexpr (has_reco_cent<C>) {
      registryMC.fill(HIST("Events/Centrality/McStatus"), static_cast<float>(McStatus::kMcGenAll), cGen, occGen);
    } else {
      registryMC.fill(HIST("Events/McStatus"), static_cast<float>(McStatus::kMcGenAll), occGen);
    }

    auto nCharged = countPart(particles);
    if constexpr (has_reco_cent<C>) {
      registryMC.fill(HIST("Events/Centrality/NtrkZvtxGen_t"), nCharged, mcCollision.posZ(), cGen, occGen);
    } else {
      registryMC.fill(HIST("Events/NtrkZvtxGen_t"), nCharged, mcCollision.posZ(), occGen);
    }
    if (collisions.size() == 0) {
      if constexpr (has_reco_cent<C>) {
        registryMC.fill(HIST("Events/Centrality/NotFoundEventZvtx"), mcCollision.posZ(), cGen, occGen);
      } else {
        registryMC.fill(HIST("Events/NotFoundEventZvtx"), mcCollision.posZ(), occGen);
      }
    }
    fillHistMC<has_reco_cent<C>>(particles, mcCollision.posZ(), cGen, occGen);
  }

  void processMCbestInclusive(CollsMCExtraMult::iterator const& mccollision,
                              soa::SmallGroups<soa::Join<Colls, aod::McCollisionLabels>> const& collisions,
                              aod::McParticles const& particles,
                              MftTracksLabeled const& tracks,
                              MftBestTracksLabeled const& besttracks)
  {
    processMCbest<CollsMCExtraMult, Colls>(mccollision, collisions, particles, tracks, besttracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processMCbestInclusive, "Count MC particles using aod::BestCollisionsFwd3d (inclusive)", false);

  void processMCbestCentFT0C(CollsMCExtraMult::iterator const& mccollision,
                             soa::SmallGroups<soa::Join<CollsCentFT0C, aod::McCollisionLabels>> const& collisions,
                             aod::McParticles const& particles,
                             MftTracksLabeled const& tracks,
                             MftBestTracksLabeled const& besttracks)
  {
    processMCbest<CollsMCExtraMult, CollsCentFT0C>(mccollision, collisions, particles, tracks, besttracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processMCbestCentFT0C, "Count MC particles in FT0C centrality bins using aod::BestCollisionsFwd3d", false);

  Preslice<MftTracksLabeled> filtMcTrkperCol = o2::aod::fwdtrack::collisionId;

  /// @brief process function to calculate tracking efficiency
  template <typename MC, typename C>
  void processEfficiency(typename MC::iterator const& mcCollision,
                         soa::SmallGroups<soa::Join<C, aod::McCollisionLabels>> const& collisions,
                         aod::McParticles const& particles,
                         MftTracksLabeled const& tracks)
  {
    LOGP(debug, "MC col {} has {} reco cols", mcCollision.globalIndex(), collisions.size());
    auto cRec = -999.;
    auto occRec = -999.;
    bool gtOneRec = false;
    for (const auto& collision : collisions) {
      if (!isGoodEvent<false>(collision)) {
        continue;
      }
      if (gConf.cfgRemoveSplitVertex && collision.globalIndex() != mcCollision.bestCollisionIndex()) {
        continue;
      }
      gtOneRec = true;
      cRec = getRecoCent(collision);
      occRec = getOccupancy(collision, eventCuts.occupancyEstimator);
    }

    if constexpr (has_reco_cent<C>) {
      registryMC.fill(HIST("Events/Centrality/hVtxZGen"), mcCollision.posZ(), static_cast<float>(GenTrkType::kGenAll), cRec, occRec);
      if (gtOneRec) {
        registryMC.fill(HIST("Events/Centrality/hVtxZGen"), mcCollision.posZ(), static_cast<float>(GenTrkType::kGenRecEvt), cRec, occRec);
      }
    } else {
      registryMC.fill(HIST("Events/hVtxZGen"), mcCollision.posZ(), static_cast<float>(GenTrkType::kGenAll), occRec);
      if (gtOneRec) {
        registryMC.fill(HIST("Events/hVtxZGen"), mcCollision.posZ(), static_cast<float>(GenTrkType::kGenAll), occRec);
      }
    }
    for (const auto& particle : particles) {
      if (!isChrgParticle(particle.pdgCode())) {
        continue;
      }
      if (gConf.cfgUseParticleSel && !isParticleSelected(particle)) {
        continue;
      }
      if constexpr (has_reco_cent<C>) {
        if (particle.eta() > trackCuts.minEta && particle.eta() < trackCuts.maxEta) {
          registryMC.fill(HIST("Tracks/Centrality/hEffGen"), particle.pt(), particle.phi(), particle.eta(), mcCollision.posZ(), static_cast<float>(GenTrkType::kGenAll), cRec, occRec);
          if (gtOneRec) {
            registryMC.fill(HIST("Tracks/Centrality/hEffGen"), particle.pt(), particle.phi(), particle.eta(), mcCollision.posZ(), static_cast<float>(GenTrkType::kGenRecEvt), cRec, occRec);
          }
        }
      } else {
        if (particle.eta() > trackCuts.minEta && particle.eta() < trackCuts.maxEta) {
          registryMC.fill(HIST("Tracks/hEffGen"), particle.pt(), particle.phi(), particle.eta(), mcCollision.posZ(), static_cast<float>(GenTrkType::kGenAll), occRec);
          if (gtOneRec) {
            registryMC.fill(HIST("Tracks/hEffGen"), particle.pt(), particle.phi(), particle.eta(), mcCollision.posZ(), static_cast<float>(GenTrkType::kGenRecEvt), occRec);
          }
        }
      }
    }
    for (const auto& collision : collisions) {
      if constexpr (has_reco_cent<C>) {
        registryMC.fill(HIST("Events/Centrality/hMcEffStatus"), static_cast<int>(McEffStatus::kMcEffAll), getRecoCent(collision), cRec, occRec);
      } else {
        registryMC.fill(HIST("Events/hMcEffStatus"), static_cast<int>(McEffStatus::kMcEffAll), occRec);
      }
      if (!isGoodEvent<false>(collision)) {
        continue;
      }
      if constexpr (has_reco_cent<C>) {
        registryMC.fill(HIST("Events/Centrality/hMcEffStatus"), static_cast<int>(McEffStatus::kMcEffSel), getRecoCent(collision), cRec, occRec);
      } else {
        registryMC.fill(HIST("Events/hMcEffStatus"), static_cast<int>(McEffStatus::kMcEffSel), occRec);
      }
      if (!collision.has_mcCollision()) {
        continue;
      }
      if constexpr (has_reco_cent<C>) {
        registryMC.fill(HIST("Events/Centrality/hMcEffStatus"), static_cast<int>(McEffStatus::kMcEffHasMcColl), getRecoCent(collision), cRec, occRec);
      } else {
        registryMC.fill(HIST("Events/hMcEffStatus"), static_cast<int>(McEffStatus::kMcEffHasMcColl), occRec);
      }
      if (gConf.cfgRemoveSplitVertex && collision.globalIndex() != mcCollision.bestCollisionIndex()) {
        continue;
      }
      if constexpr (has_reco_cent<C>) {
        registryMC.fill(HIST("Events/Centrality/hMcEffStatus"), static_cast<int>(McEffStatus::kMcEffNoSplitVtx), getRecoCent(collision), cRec, occRec);
      } else {
        registryMC.fill(HIST("Events/hMcEffStatus"), static_cast<int>(McEffStatus::kMcEffNoSplitVtx), occRec);
      }
      if constexpr (has_reco_cent<C>) {
        registryMC.fill(HIST("Events/Centrality/hVtxZRec"), collision.posZ(), getRecoCent(collision), cRec, occRec);
      } else {
        registryMC.fill(HIST("Events/hVtxZRec"), collision.posZ(), occRec);
      }
      auto perColTrks = tracks.sliceBy(filtMcTrkperCol, collision.globalIndex());
      for (auto const& track : perColTrks) {
        if (!isTrackSelected<false>(track)) {
          continue;
        }
        if (track.has_mcParticle() && track.mcMask() == trackCuts.selMcMask) {
          const auto& particle = track.template mcParticle_as<aod::McParticles>();
          if (!isChrgParticle(particle.pdgCode())) {
            continue;
          }
          if (gConf.cfgUseParticleSel && !isParticleSelected(particle)) {
            continue;
          }
          if (collision.mcCollisionId() != particle.mcCollisionId()) {
            continue;
          }
          if constexpr (has_reco_cent<C>) {
            registryMC.fill(HIST("Tracks/Centrality/hEffRec"), particle.pt(), particle.phi(), particle.eta(), collision.posZ(), static_cast<float>(RecTrkType::kRecAll), getRecoCent(collision), getOccupancy(collision, eventCuts.occupancyEstimator));
          } else {
            registryMC.fill(HIST("Tracks/hEffRec"), particle.pt(), particle.phi(), particle.eta(), collision.posZ(), static_cast<float>(RecTrkType::kRecAll), getOccupancy(collision, eventCuts.occupancyEstimator));
            registryMC.fill(HIST("Tracks/hEtaRes"), particle.eta(), (track.eta() - particle.eta()) / particle.eta(), getOccupancy(collision, eventCuts.occupancyEstimator));
          }
        } else {
          if constexpr (has_reco_cent<C>) {
            registryMC.fill(HIST("Tracks/Centrality/hEffRec"), track.pt(), track.phi(), track.eta(), collision.posZ(), static_cast<float>(RecTrkType::kRecFake), getRecoCent(collision), getOccupancy(collision, eventCuts.occupancyEstimator));
          } else {
            registryMC.fill(HIST("Tracks/hEffRec"), track.pt(), track.phi(), track.eta(), collision.posZ(), static_cast<float>(RecTrkType::kRecFake), getOccupancy(collision, eventCuts.occupancyEstimator));
          }
        }
      }
    }
  }

  void processEfficiencyInclusive(CollsMCExtraMult::iterator const& mccollision,
                                  soa::SmallGroups<soa::Join<Colls, aod::McCollisionLabels>> const& collisions,
                                  aod::McParticles const& particles,
                                  MftTracksLabeled const& tracks)
  {
    processEfficiency<CollsMCExtraMult, Colls>(mccollision, collisions, particles, tracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processEfficiencyInclusive, "Process efficiencies (inclusive)", false);

  void processEfficiencyCentFT0C(CollsMCExtraMult::iterator const& mccollision,
                                 soa::SmallGroups<soa::Join<CollsCentFT0C, aod::McCollisionLabels>> const& collisions,
                                 aod::McParticles const& particles,
                                 MftTracksLabeled const& tracks)
  {
    processEfficiency<CollsMCExtraMult, CollsCentFT0C>(mccollision, collisions, particles, tracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processEfficiencyCentFT0C, "Process efficiencies in FT0C centrality bins", false);

  /// @brief process function to calculate tracking efficiency based on BestCollisionsFwd3d in FT0C bins
  template <typename MC, typename C>
  void processEfficiencyBest(typename MC::iterator const& mcCollision,
                             soa::SmallGroups<soa::Join<C, aod::McCollisionLabels>> const& collisions,
                             aod::McParticles const& particles,
                             MftBestTracksLabeled const& besttracks)
  {
    LOGP(debug, "MC col {} has {} reco cols", mcCollision.globalIndex(), collisions.size());

    buildLookupTable(collisions);

    auto cRec = -999.;
    auto occRec = -999.;
    bool gtOneRec = false;
    for (const auto& collision : collisions) {
      if (!isGoodEvent<false>(collision)) {
        continue;
      }
      if (gConf.cfgRemoveSplitVertex && collision.globalIndex() != mcCollision.bestCollisionIndex()) {
        continue;
      }
      gtOneRec = true;
      cRec = getRecoCent(collision);
      occRec = getOccupancy(collision, eventCuts.occupancyEstimator);
    }

    if constexpr (has_reco_cent<C>) {
      registryMC.fill(HIST("Events/hVtxZGen"), mcCollision.posZ(), static_cast<float>(GenTrkType::kGenAll), cRec, occRec);
      if (gtOneRec) {
        registryMC.fill(HIST("Events/hVtxZGen"), mcCollision.posZ(), static_cast<float>(GenTrkType::kGenRecEvt), cRec, occRec);
      }
    } else {
      registryMC.fill(HIST("Events/hVtxZGen"), mcCollision.posZ(), static_cast<float>(GenTrkType::kGenAll), occRec);
      if (gtOneRec) {
        registryMC.fill(HIST("Events/hVtxZGen"), mcCollision.posZ(), static_cast<float>(GenTrkType::kGenAll), occRec);
      }
    }

    for (const auto& particle : particles) {
      if (!isChrgParticle(particle.pdgCode())) {
        continue;
      }
      if (gConf.cfgUseParticleSel && !isParticleSelected(particle)) {
        continue;
      }
      if constexpr (has_reco_cent<C>) {
        if (particle.eta() > trackCuts.minEta && particle.eta() < trackCuts.maxEta) {
          registryMC.fill(HIST("Tracks/Centrality/hEffGen"), particle.pt(), particle.phi(), particle.eta(), mcCollision.posZ(), static_cast<float>(GenTrkType::kGenAll), cRec, occRec);
          if (gtOneRec) {
            registryMC.fill(HIST("Tracks/Centrality/hEffGen"), particle.pt(), particle.phi(), particle.eta(), mcCollision.posZ(), static_cast<float>(GenTrkType::kGenRecEvt), cRec, occRec);
          }
        }
      } else {
        if (particle.eta() > trackCuts.minEta && particle.eta() < trackCuts.maxEta) {
          registryMC.fill(HIST("Tracks/hEffGen"), particle.pt(), particle.phi(), particle.eta(), mcCollision.posZ(), static_cast<float>(GenTrkType::kGenAll), occRec);
          if (gtOneRec) {
            registryMC.fill(HIST("Tracks/hEffGen"), particle.pt(), particle.phi(), particle.eta(), mcCollision.posZ(), static_cast<float>(GenTrkType::kGenRecEvt), occRec);
          }
        }
      }
    }

    for (const auto& collision : collisions) {
      if constexpr (has_reco_cent<C>) {
        registryMC.fill(HIST("Events/Centrality/hMcEffStatus"), static_cast<int>(McEffStatus::kMcEffAll), getRecoCent(collision), cRec, occRec);
      } else {
        registryMC.fill(HIST("Events/hMcEffStatus"), static_cast<int>(McEffStatus::kMcEffAll), occRec);
      }
      if (!isGoodEvent<false>(collision)) {
        continue;
      }
      if constexpr (has_reco_cent<C>) {
        registryMC.fill(HIST("Events/Centrality/hMcEffStatus"), static_cast<int>(McEffStatus::kMcEffSel), getRecoCent(collision), cRec, occRec);
      } else {
        registryMC.fill(HIST("Events/hMcEffStatus"), static_cast<int>(McEffStatus::kMcEffSel), occRec);
      }
      if (!collision.has_mcCollision()) {
        continue;
      }
      if constexpr (has_reco_cent<C>) {
        registryMC.fill(HIST("Events/Centrality/hMcEffStatus"), static_cast<int>(McEffStatus::kMcEffHasMcColl), getRecoCent(collision), cRec, occRec);
      } else {
        registryMC.fill(HIST("Events/hMcEffStatus"), static_cast<int>(McEffStatus::kMcEffHasMcColl), occRec);
      }
      if (gConf.cfgRemoveSplitVertex && collision.globalIndex() != mcCollision.bestCollisionIndex()) {
        continue;
      }
      if constexpr (has_reco_cent<C>) {
        registryMC.fill(HIST("Events/Centrality/hMcEffStatus"), static_cast<int>(McEffStatus::kMcEffNoSplitVtx), getRecoCent(collision), cRec, occRec);
      } else {
        registryMC.fill(HIST("Events/hMcEffStatus"), static_cast<int>(McEffStatus::kMcEffNoSplitVtx), occRec);
      }
      if constexpr (has_reco_cent<C>) {
        registryMC.fill(HIST("Events/Centrality/hVtxZRec"), collision.posZ(), getRecoCent(collision), cRec, occRec);
      } else {
        registryMC.fill(HIST("Events/hVtxZRec"), collision.posZ(), occRec);
      }

      auto perCollisionASample = besttracks.sliceBy(perColU, collision.globalIndex());
      for (auto const& atrack : perCollisionASample) {
        if (!isBestTrackSelected<false>(atrack)) {
          continue;
        }
        auto itrack = atrack.template mfttrack_as<MftBestTracksLabeled>();
        if (!isTrackSelected<false>(itrack)) {
          continue;
        }
        if (!itrack.has_collision()) {
          continue;
        }
        if (gConf.cfgRemoveReassigned) {
          if (itrack.collisionId() != atrack.bestCollisionId()) {
            continue;
          }
        }
        if (itrack.collisionId() >= 0 && itrack.has_mcParticle() && itrack.mcMask() == trackCuts.selMcMask) {
          const auto& particle = itrack.template mcParticle_as<aod::McParticles>();
          if (!isChrgParticle(particle.pdgCode())) {
            continue;
          }
          if (gConf.cfgUseParticleSel && !isParticleSelected(particle)) {
            continue;
          }
          if (eventCuts.useZDiffCut) {
            if (std::abs(collision.posZ() - atrack.mcParticle().mcCollision().posZ()) > eventCuts.maxZvtxDiff) {
              continue;
            }
          }
          const int bestRecColl = atrack.bestCollisionId();
          if (!mapMcCollIdPerRecColl.contains(bestRecColl)) {
            continue;
          }
          int64_t mcCollIdRec = mapMcCollIdPerRecColl.find(bestRecColl)->second;
          if (mcCollIdRec != particle.mcCollisionId()) {
            continue;
          }
          if constexpr (has_reco_cent<C>) {
            registryMC.fill(HIST("Tracks/Centrality/hEffRec"), particle.pt(), particle.phi(), particle.eta(), collision.posZ(), static_cast<float>(RecTrkType::kRecAll), getRecoCent(collision), getOccupancy(collision, eventCuts.occupancyEstimator));
          } else {
            registryMC.fill(HIST("Tracks/hEffRec"), particle.pt(), particle.phi(), particle.eta(), collision.posZ(), static_cast<float>(RecTrkType::kRecAll), getOccupancy(collision, eventCuts.occupancyEstimator));
            registryMC.fill(HIST("Tracks/hEtaRes"), particle.eta(), (itrack.eta() - particle.eta()) / particle.eta(), getOccupancy(collision, eventCuts.occupancyEstimator));
          }
        } else {
          if constexpr (has_reco_cent<C>) {
            registryMC.fill(HIST("Tracks/Centrality/hEffRec"), itrack.pt(), itrack.phi(), itrack.eta(), collision.posZ(), static_cast<float>(RecTrkType::kRecFake), getRecoCent(collision), getOccupancy(collision, eventCuts.occupancyEstimator));
          } else {
            registryMC.fill(HIST("Tracks/hEffRec"), itrack.pt(), itrack.phi(), itrack.eta(), collision.posZ(), static_cast<float>(RecTrkType::kRecFake), getOccupancy(collision, eventCuts.occupancyEstimator));
          }
        }
      }
    }
  }

  void processEfficiencyBestInclusive(CollsMCExtraMult::iterator const& mccollision,
                                      soa::SmallGroups<soa::Join<Colls, aod::McCollisionLabels>> const& collisions,
                                      aod::McParticles const& particles,
                                      MftBestTracksLabeled const& besttracks)
  {
    processEfficiencyBest<CollsMCExtraMult, Colls>(mccollision, collisions, particles, besttracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processEfficiencyBestInclusive, "Process tracking efficiency (inclusive, based on BestCollisionsFwd3d)", false);

  void processEfficiencyBestCentFT0C(CollsMCExtraMult::iterator const& mccollision,
                                     soa::SmallGroups<soa::Join<CollsCentFT0C, aod::McCollisionLabels>> const& collisions,
                                     aod::McParticles const& particles,
                                     MftBestTracksLabeled const& besttracks)
  {
    processEfficiencyBest<CollsMCExtraMult, CollsCentFT0C>(mccollision, collisions, particles, besttracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processEfficiencyBestCentFT0C, "Process tracking efficiency (in FT0 centrality bins, based on BestCollisionsFwd3d)", false);

  // Partition<ParticlesIdx> primariesIdx = ifnode(nsqrt(aod::mcparticle::vx * aod::mcparticle::vx + aod::mcparticle::vy * aod::mcparticle::vy) > 5.f, false, ncheckbit(aod::mcparticle_v2::storedFlags, (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary)) && (aod::mcparticle::eta < trackCuts.maxEta) && (aod::mcparticle::eta > trackCuts.minEta);

  /// @brief process template function to calculate tracking efficiency (indexed as particle-to-MFT-tracks)
  template <typename MC, typename C>
  void processEfficiencyIdx(typename MC::iterator const& mcCollision,
                            soa::SmallGroups<soa::Join<C, aod::McCollisionLabels>> const& collisions,
                            FiltParticlesIdx const& particles,
                            MftTracksLabeled const& /*tracks*/)
  {
    LOGP(debug, "MC col {} has {} reco cols", mcCollision.globalIndex(), collisions.size());
    for (const auto& collision : collisions) {
      if (!isGoodEvent<false>(collision)) {
        continue;
      }
      if (gConf.cfgRemoveSplitVertex && collision.globalIndex() != mcCollision.bestCollisionIndex()) {
        continue;
      }
      float cRec = getRecoCent(collision);
      float occRec = getOccupancy(collision, eventCuts.occupancyEstimator);

      auto partsPerCol = particles.sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);
      for (auto const& particle : partsPerCol) {
        if (!isChrgParticle(particle.pdgCode())) {
          continue;
        }
        if (collision.mcCollisionId() != particle.mcCollisionId()) {
          continue;
        }
        // MC gen
        if constexpr (has_reco_cent<C>) {
          if (particle.eta() > trackCuts.minEta && particle.eta() < trackCuts.maxEta && cutInPhi(particle)) {
            if (std::abs(mcCollision.posZ()) < eventCuts.maxZvtx) {
              registryMC.fill(HIST("Tracks/Centrality/hEffIdxGen"), particle.pt(), particle.eta(), static_cast<float>(GenIdxTrkType::kGenIdxAll), cRec, occRec);
            }
          }
        } else {
          if (particle.eta() > trackCuts.minEta && particle.eta() < trackCuts.maxEta && cutInPhi(particle)) {
            if (std::abs(mcCollision.posZ()) < eventCuts.maxZvtx) {
              registryMC.fill(HIST("Tracks/hEffIdxGen"), particle.pt(), particle.eta(), static_cast<float>(GenIdxTrkType::kGenIdxAll), occRec);
            }
          }
        }
        // MC rec
        if (particle.has_mfttracks()) {
          auto iscounted = false;
          auto ncnt = 0;
          auto relatedTracks = particle.template mfttracks_as<MftTracksLabeled>();
          for (auto const& track : relatedTracks) {
            if (!isTrackSelected<false>(track)) {
              continue;
            }
            if (track.mcMask() != trackCuts.selMcMask) {
              continue;
            }
            ++ncnt;
            if constexpr (has_reco_cent<C>) {
              if (track.eta() > trackCuts.minEta && track.eta() < trackCuts.maxEta) {
                if (!iscounted) { // primaries
                  if (std::abs(mcCollision.posZ()) < eventCuts.maxZvtx) {
                    registryMC.fill(HIST("Tracks/Centrality/hEffIdxRec"), particle.pt(), particle.eta(), static_cast<float>(RecIdxTrkType::kRecIdxPrim), cRec, occRec);
                  }
                  iscounted = true;
                }
              }
              registryMC.fill(HIST("Tracks/Centrality/NmftTrkPerPart"), ncnt, cRec, occRec);
              if (ncnt > 1) { // secondaries
                if (track.eta() > trackCuts.minEta && track.eta() < trackCuts.maxEta) {
                  registryMC.fill(HIST("Tracks/Centrality/hEffIdxRec"), particle.pt(), particle.eta(), static_cast<float>(RecIdxTrkType::kRecIdxSec), cRec, occRec);
                }
              }
            } else {
              if (track.eta() > trackCuts.minEta && track.eta() < trackCuts.maxEta) {
                if (!iscounted) { // primaries
                  if (std::abs(mcCollision.posZ()) < eventCuts.maxZvtx) {
                    registryMC.fill(HIST("Tracks/hEffIdxRec"), particle.pt(), particle.eta(), static_cast<float>(RecIdxTrkType::kRecIdxPrim), occRec);
                  }
                  iscounted = true;
                }
              }
              registryMC.fill(HIST("Tracks/NmftTrkPerPart"), ncnt, occRec);
              if (ncnt > 1) { // secondaries
                if (track.eta() > trackCuts.minEta && track.eta() < trackCuts.maxEta) {
                  registryMC.fill(HIST("Tracks/hEffIdxRec"), particle.pt(), particle.eta(), static_cast<float>(RecIdxTrkType::kRecIdxSec), occRec);
                }
              }
            }
          }
          if (relatedTracks.size() > 1) { // duplicates
            if constexpr (has_reco_cent<C>) {
              registryMC.fill(HIST("Tracks/Centrality/hEffIdxGen"), particle.pt(), particle.eta(), static_cast<float>(GenIdxTrkType::kGenIdxDupl), cRec, occRec);
              for (auto const& track : relatedTracks) {
                registryMC.fill(HIST("Tracks/Centrality/hEffIdxRec"), track.pt(), track.eta(), static_cast<float>(RecIdxTrkType::kRecIdxDupl), cRec, occRec);
              }
            } else {
              registryMC.fill(HIST("Tracks/hEffIdxGen"), particle.pt(), particle.eta(), static_cast<float>(GenIdxTrkType::kGenIdxDupl), occRec);
              for (auto const& track : relatedTracks) {
                registryMC.fill(HIST("Tracks/hEffIdxRec"), track.pt(), track.eta(), static_cast<float>(RecIdxTrkType::kRecIdxDupl), occRec);
              }
            }
          }
        } else {
          // MC FAKES
          if constexpr (has_reco_cent<C>) {
            if (particle.eta() > trackCuts.minEta && particle.eta() < trackCuts.maxEta) {
              if (std::abs(mcCollision.posZ()) < eventCuts.maxZvtx) {
                registryMC.fill(HIST("Tracks/Centrality/hEffIdxGen"), particle.pt(), particle.eta(), static_cast<float>(GenIdxTrkType::kGenIdxFake), cRec, occRec);
              }
            }
          } else {
            if (particle.eta() > trackCuts.minEta && particle.eta() < trackCuts.maxEta) {
              if (std::abs(mcCollision.posZ()) < eventCuts.maxZvtx) {
                registryMC.fill(HIST("Tracks/hEffIdxGen"), particle.pt(), particle.eta(), static_cast<float>(GenIdxTrkType::kGenIdxFake), occRec);
              }
            }
          }
        }
      }
    }
  }

  void processEfficiencyIdxInlusive(CollsMCExtraMult::iterator const& mccollision,
                                    soa::SmallGroups<soa::Join<Colls, aod::McCollisionLabels>> const& collisions,
                                    FiltParticlesIdx const& particles,
                                    MftTracksLabeled const& tracks)
  {
    processEfficiencyIdx<CollsMCExtraMult, Colls>(mccollision, collisions, particles, tracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processEfficiencyIdxInlusive, "Process tracking efficiency (inclusive, indexed)", false);

  void processEfficiencyIdxCentFT0C(CollsMCExtraMult::iterator const& mccollision,
                                    soa::SmallGroups<soa::Join<CollsCentFT0C, aod::McCollisionLabels>> const& collisions,
                                    FiltParticlesIdx const& particles,
                                    MftTracksLabeled const& tracks)
  {
    processEfficiencyIdx<CollsMCExtraMult, CollsCentFT0C>(mccollision, collisions, particles, tracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processEfficiencyIdxCentFT0C, "Process tracking efficiency (in FT0C centrality bins, indexed)", false);

  template <typename MC, typename C>
  void processEfficiencyIdxBest(typename MC::iterator const& mcCollision,
                                soa::SmallGroups<soa::Join<C, aod::McCollisionLabels>> const& collisions,
                                FiltParticlesIdx const& particles,
                                MftBestTracksLabeled const& /*atracks*/)
  {
    LOGP(debug, "MC col {} has {} reco cols", mcCollision.globalIndex(), collisions.size());

    buildLookupTable(collisions);

    for (const auto& collision : collisions) {
      if (!isGoodEvent<false>(collision)) {
        continue;
      }
      if (gConf.cfgRemoveSplitVertex && collision.globalIndex() != mcCollision.bestCollisionIndex()) {
        continue;
      }
      float cRec = getRecoCent(collision);
      float occRec = getOccupancy(collision, eventCuts.occupancyEstimator);

      auto partsPerCol = particles.sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);
      for (auto const& particle : partsPerCol) {
        if (!isChrgParticle(particle.pdgCode())) {
          continue;
        }
        if (collision.mcCollisionId() != particle.mcCollisionId()) {
          continue;
        }
        // MC gen
        if constexpr (has_reco_cent<C>) {
          if (particle.eta() > trackCuts.minEta && particle.eta() < trackCuts.maxEta && cutInPhi(particle)) {
            if (std::abs(mcCollision.posZ()) < eventCuts.maxZvtx) {
              registryMC.fill(HIST("Tracks/Centrality/hEffIdxGen"), particle.pt(), particle.eta(), static_cast<float>(GenIdxTrkType::kGenIdxAll), cRec, occRec);
            }
          }
        } else {
          if (particle.eta() > trackCuts.minEta && particle.eta() < trackCuts.maxEta && cutInPhi(particle)) {
            if (std::abs(mcCollision.posZ()) < eventCuts.maxZvtx) {
              registryMC.fill(HIST("Tracks/hEffIdxGen"), particle.pt(), particle.eta(), static_cast<float>(GenIdxTrkType::kGenIdxAll), occRec);
            }
          }
        }
        // MC rec
        if (particle.has_mfttracks()) {
          auto iscounted = false;
          auto ncnt = 0;
          auto relatedTracks = particle.template mfttracks_as<MftBestTracksLabeled>();
          for (auto const& atrack : relatedTracks) {
            if (!isBestTrackSelected<false>(atrack)) {
              continue;
            }
            const int bestRecColl = atrack.bestCollisionId();
            if (!mapMcCollIdPerRecColl.contains(bestRecColl)) {
              continue;
            }
            int64_t mcCollIdRec = mapMcCollIdPerRecColl.find(bestRecColl)->second;
            if (mcCollIdRec != particle.mcCollisionId()) {
              continue;
            }
            ++ncnt;
            if constexpr (has_reco_cent<C>) {
              if (atrack.eta() > trackCuts.minEta && atrack.eta() < trackCuts.maxEta) {
                if (!iscounted) { // primaries
                  if (std::abs(mcCollision.posZ()) < eventCuts.maxZvtx) {
                    registryMC.fill(HIST("Tracks/Centrality/hEffIdxRec"), particle.pt(), particle.eta(), static_cast<float>(RecIdxTrkType::kRecIdxPrim), cRec, occRec);
                  }
                  iscounted = true;
                }
              }
              registryMC.fill(HIST("Tracks/Centrality/NmftTrkPerPart"), ncnt, cRec, occRec);
              if (ncnt > 1) { // secondaries
                if (atrack.eta() > trackCuts.minEta && atrack.eta() < trackCuts.maxEta) {
                  registryMC.fill(HIST("Tracks/Centrality/hEffIdxRec"), particle.pt(), particle.eta(), static_cast<float>(RecIdxTrkType::kRecIdxSec), cRec, occRec);
                }
              }
            } else {
              if (atrack.eta() > trackCuts.minEta && atrack.eta() < trackCuts.maxEta) {
                if (!iscounted) { // primaries
                  if (std::abs(mcCollision.posZ()) < eventCuts.maxZvtx) {
                    registryMC.fill(HIST("Tracks/hEffIdxRec"), particle.pt(), particle.eta(), static_cast<float>(RecIdxTrkType::kRecIdxPrim), occRec);
                  }
                  iscounted = true;
                }
              }
              registryMC.fill(HIST("Tracks/NmftTrkPerPart"), ncnt, occRec);
              if (ncnt > 1) { // secondaries
                if (atrack.eta() > trackCuts.minEta && atrack.eta() < trackCuts.maxEta) {
                  registryMC.fill(HIST("Tracks/hEffIdxRec"), particle.pt(), particle.eta(), static_cast<float>(RecIdxTrkType::kRecIdxSec), occRec);
                }
              }
            }
          }
          if (relatedTracks.size() > 1) { // duplicates
            if constexpr (has_reco_cent<C>) {
              registryMC.fill(HIST("Tracks/Centrality/hEffIdxGen"), particle.pt(), particle.eta(), static_cast<float>(GenIdxTrkType::kGenIdxDupl), cRec, occRec);
              for (auto const& track : relatedTracks) {
                registryMC.fill(HIST("Tracks/Centrality/hEffIdxRec"), track.pt(), track.eta(), static_cast<float>(RecIdxTrkType::kRecIdxDupl), cRec, occRec);
              }
            } else {
              registryMC.fill(HIST("Tracks/hEffIdxGen"), particle.pt(), particle.eta(), static_cast<float>(GenIdxTrkType::kGenIdxDupl), occRec);
              for (auto const& track : relatedTracks) {
                registryMC.fill(HIST("Tracks/hEffIdxRec"), track.pt(), track.eta(), static_cast<float>(RecIdxTrkType::kRecIdxDupl), occRec);
              }
            }
          }
        } else {
          // MC FAKES
          if constexpr (has_reco_cent<C>) {
            if (particle.eta() > trackCuts.minEta && particle.eta() < trackCuts.maxEta) {
              if (std::abs(mcCollision.posZ()) < eventCuts.maxZvtx) {
                registryMC.fill(HIST("Tracks/Centrality/hEffIdxGen"), particle.pt(), particle.eta(), static_cast<float>(GenIdxTrkType::kGenIdxFake), cRec, occRec);
              }
            }
          } else {
            if (particle.eta() > trackCuts.minEta && particle.eta() < trackCuts.maxEta) {
              if (std::abs(mcCollision.posZ()) < eventCuts.maxZvtx) {
                registryMC.fill(HIST("Tracks/hEffIdxGen"), particle.pt(), particle.eta(), static_cast<float>(GenIdxTrkType::kGenIdxFake), occRec);
              }
            }
          }
        }
      }
    }
  }

  void processEfficiencyIdxBestInlusive(CollsMCExtraMult::iterator const& mccollision,
                                        soa::SmallGroups<soa::Join<Colls, aod::McCollisionLabels>> const& collisions,
                                        FiltParticlesIdx const& particles,
                                        MftBestTracksLabeled const& atracks)
  {
    processEfficiencyIdxBest<CollsMCExtraMult, Colls>(mccollision, collisions, particles, atracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processEfficiencyIdxBestInlusive, "Process tracking efficiency best (inclusive, indexed)", false);

  void processEfficiencyIdxBestCentFT0C(CollsMCExtraMult::iterator const& mccollision,
                                        soa::SmallGroups<soa::Join<CollsCentFT0C, aod::McCollisionLabels>> const& collisions,
                                        FiltParticlesIdx const& particles,
                                        MftBestTracksLabeled const& atracks)
  {
    processEfficiencyIdxBest<CollsMCExtraMult, CollsCentFT0C>(mccollision, collisions, particles, atracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processEfficiencyIdxBestCentFT0C, "Process tracking efficiency best (in FT0C centrality bins, indexed)", false);

  /// @brief process function to calculate signal loss based on MC
  void processMcSgnEvtLossCentFT0C(CollsMCExtraMult::iterator const& mcCollision,
                                   soa::SmallGroups<soa::Join<CollsCentFT0C, aod::McCollisionLabels>> const& collisions,
                                   aod::McParticles const& particles)
  {
    LOGP(debug, "MC col {} has {} reco cols", mcCollision.globalIndex(), collisions.size());
    registryMC.fill(HIST("Events/hNchGen"), mcCollision.multMCFT0C(), static_cast<float>(EvtLossType::kGenAll));
    registryMC.fill(HIST("Events/hEvtMcGen"), 0.5);
    if (eventCuts.useZVtxCutMC && (std::abs(mcCollision.posZ()) >= eventCuts.maxZvtx)) {
      return;
    }
    registryMC.fill(HIST("Events/hEvtMcGen"), 1.5);
    // At least one generated primary in MFT acceptance + TVX triggered collisions
    if (gConf.cfgUseInelgt0wMFT && !isInelGt0wMft(particles)) {
      return;
    }
    registryMC.fill(HIST("Events/hEvtMcGen"), 2.5);
    if (eventCuts.useInelgt0wTVX && !(mcCollision.multMCFT0C() > 0 && mcCollision.multMCFT0A() > 0)) {
      return;
    }
    registryMC.fill(HIST("Events/hEvtMcGen"), 3.5);
    registryMC.fill(HIST("Events/hNchGen"), mcCollision.multMCFT0C(), static_cast<float>(EvtLossType::kGenSel)); // Evt loss den

    bool gtZeroColl = false;
    auto maxNcontributors = -1;
    float cRec = -999.f;
    for (auto const& collision : collisions) {
      if (!isGoodEvent<false>(collision)) {
        continue;
      }
      if (std::abs(collision.posZ()) >= eventCuts.maxZvtx) {
        continue;
      }
      registryMC.fill(HIST("Events/hNchGen"), mcCollision.multMCFT0C(), static_cast<float>(EvtLossType::kGenSplit));
      registryMC.fill(HIST("Events/hMultGenVsCentSplit"), getRecoCent(collision), mcCollision.multMCFT0C());
      if (maxNcontributors < collision.numContrib()) {
        maxNcontributors = collision.numContrib();
        cRec = getRecoCent(collision);
      }
      gtZeroColl = true;
    }

    auto perCollMCsample = mcSample->sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);
    auto multMCNParticlesEtaMFT = countPart(perCollMCsample);

    registryMC.fill(HIST("Events/EvtSigLossStatus"), 1., cRec); // Evt split den
    registryMC.fill(HIST("Events/hMultGenVsCent"), cRec, mcCollision.multMCFT0C());
    registryMC.fill(HIST("Events/hMultGenVsCentNParticlesEta05"), cRec, mcCollision.multMCNParticlesEta05());
    registryMC.fill(HIST("Events/hMultGenVsCentNParticlesEtaMFT"), cRec, multMCNParticlesEtaMFT);

    if (gtZeroColl) {
      registryMC.fill(HIST("Events/hNchGen"), mcCollision.multMCFT0C(), static_cast<float>(EvtLossType::kGenRecEvt)); // Evt loss num
      registryMC.fill(HIST("Events/EvtSigLossStatus"), 2., cRec);                                                     // Evt split num
      registryMC.fill(HIST("Events/hMultGenVsCentRec"), cRec, mcCollision.multMCFT0C());
      registryMC.fill(HIST("Events/hMultGenVsCentRecNParticlesEta05"), cRec, mcCollision.multMCNParticlesEta05());
      registryMC.fill(HIST("Events/hMultGenVsCentRecNParticlesEtaMFT"), cRec, multMCNParticlesEtaMFT);
    }
    if (collisions.size() == 0) {
      registryMC.fill(HIST("Events/EvtSigLossStatus"), 3., cRec);
    }
    for (auto const& particle : particles) {
      if (!isChrgParticle(particle.pdgCode())) {
        continue;
      }
      if (gConf.cfgUseParticleSel && !isParticleSelected(particle)) {
        continue;
      }
      registryMC.fill(HIST("Tracks/hEtaVsNchGen"), particle.eta(), mcCollision.multMCFT0C()); // Sgn loss den
      if (gtZeroColl) {
        registryMC.fill(HIST("Tracks/hEtaVsNchGenRecEvt"), particle.eta(), mcCollision.multMCFT0C()); // Sgn loss num
      }
    }
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processMcSgnEvtLossCentFT0C, "Process Signal/event loss based on MC (in FT0C centrality bins)", false);

  template <typename C>
  void processCorrelationwBestTracks(typename C::iterator const& collision,
                                     aod::MFTTracks const& /*tracks*/,
                                     soa::SmallGroups<aod::BestCollisionsFwd3d> const& besttracks)
  {
    if (!isGoodEvent<false>(collision)) {
      return;
    }

    auto nBestTrks = 0;
    for (auto const& atrack : besttracks) {
      if (gConf.cfgUseTrackSel && !isBestTrackSelected<false>(atrack)) {
        continue;
      }
      auto itrack = atrack.template mfttrack_as<aod::MFTTracks>();
      if (itrack.eta() < trackCuts.minEta || itrack.eta() > trackCuts.maxEta) {
        continue;
      }
      if (gConf.cfgUseTrackSel && !isTrackSelected<false>(itrack)) {
        continue;
      }
      nBestTrks++;
    }
    registryData.fill(HIST("Events/hMultMFTvsFT0A"), nBestTrks, collision.multFT0A());
    registryData.fill(HIST("Events/hMultMFTvsFT0C"), nBestTrks, collision.multFT0C());
    registryData.fill(HIST("Events/hNPVtracksVsFT0C"), collision.multNTracksPV(), collision.multFT0C());
    registryData.fill(HIST("Events/hMultMFTvsFV0A"), nBestTrks, collision.multFV0A());
    registryData.fill(HIST("Events/hNPVtracksVsMultMFT"), collision.multNTracksPV(), nBestTrks);
  }

  void processCorrelationwBestTracksInclusive(CollsCorr::iterator const& collision,
                                              aod::MFTTracks const& tracks,
                                              soa::SmallGroups<aod::BestCollisionsFwd3d> const& besttracks)
  {
    processCorrelationwBestTracks<CollsCorr>(collision, tracks, besttracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processCorrelationwBestTracksInclusive, "Process correlation QA based on BestCollisionsFwd3d table", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<DndetaMFTPbPb>(cfgc)};
}
