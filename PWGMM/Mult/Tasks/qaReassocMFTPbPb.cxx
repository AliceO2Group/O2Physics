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
/// \file qaReassocMFTPbPb.cxx
/// \brief  Task for performing reassociation checks in Pb-Pb collisions using MFT detector
/// \author Gyula Bencedi, gyula.bencedi@cern.ch
/// \since  Aug 2026

#include "PWGMM/Mult/Core/include/Functions.h"
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
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <MathUtils/Utils.h>

#include <TGeoGlobalMagField.h>
#include <TH1.h>
#include <TH2.h>
#include <THnSparse.h>

#include <sys/types.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <memory>
#include <numeric>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
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

enum class TrkTrkBestSel {
  trkTrkBestSelAll = 0,
  trkTrkBestSelCollID,
  trkTrkBestSelOrphan,
  trkTrkBestSelDCAxyCut,
  trkTrkBestSelDCAzCut,
  trkTrkBestSelNumReassoc,
  nTrkTrkBestSel
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

enum AmbTrkTypeAssocFlag {
  kSelAll = 0,
  kSelGoodVtx,
  kSelGoodVtxTrue,
  kSelBadVtx,
  kSelBadVtxTrue,
  kSelNonAmb,
  kSelNonAmbGoodVtx,
  kSelNonAmbGoodVtxTrue,
  kSelNonAmbBadVtx,
  kSelNonAmbBadVtxTrue,
  kSelNonAmbID,
  kSelNonAmbIDGoodVtxTrue,
  kSelNonAmbIDBadVtxTrue,
  kSelAmbID,
  kSelAmbIDGoodVtxTrue,
  kSelAmbIDBadVtxTrue,
  kSelNonAmbIDExtra,
  kSelNonAmbIDExtraGoodVtxTrue,
  kSelNonAmbIDExtraBadVtxTrue,
  kSelAmb,
  kSelAmbGoodVtxTrue,
  kSelAmbBadVtxTrue,
  kSelOrphanNull,
  nSelAmbTrkTypeAssocFlag
};

enum class VertexStatusMC {
  kNull = 0,
  kGood,
  kBad
};

enum class AssocCheckVtxType {
  kAllVtxTrue = 0,
  kAllVtxFalse,
  kAllGoodVtx,
  kAllGoodVtxTrue,
  kAllGoodVtxFalse,
  kAllBadVtx,
  kAllBadVtxTrue,
  kAllBadVtxFalse,
  nAssocVtxType
};

enum class ReassocCheckVtxType {
  kIsTrueVtxAllTrue = 0,
  kIsTrueVtxAllFalse,
  kIsTrueVtxVsGoodVtxTrue,
  kIsTrueVtxVsGoodVtxFalse,
  kIsTrueVtxVsBadVtxTrue,
  kIsTrueVtxVsBadVtxFalse,
  kIsRecGoodAllTrue,
  kIsRecGoodAllFalse,
  kIsRecGoodMatchAllTrue,
  kIsRecGoodMatchAllFalse,
  kIsRecGoodVsGoodVtxTrue,
  kIsRecGoodVsGoodVtxFalse,
  kIsRecGoodMatchVsGoodVtxTrue,
  kIsRecGoodMatchVsGoodVtxFalse,
  kIsRecGoodVsBadVtxTrue,
  kIsRecGoodVsBadVtxFalse,
  kIsRecGoodMatchVsBadVtxTrue,
  kIsRecGoodMatchVsBadVtxFalse,
  nReassocVtxType
};

enum class ReAssocMCEventStatus {
  kEvtReAsAll = 0,
  kEvtReAsSelected,
  kEvtReAsHasMcColl,
  kEvtReAsSplitVtxRemoved,
  kEvtReAsZVtxCutMC,
  nEvtReAsReAssocMCEventStatus
};

enum class ReAssocMCTrackStatus {
  kTrkReAssocAll = 0,
  kTrkBestSel,
  kTrkSel,
  kTrkHasColl,
  kTrkReassignedRemoved,
  kTrkHasMcPart,
  kTrkIdGt0,
  kTrkNonAmb,
  kTrkNonAmbGood,
  kTrkNonAmbBad,
  kTrkNonAmbID,
  kTrkNonAmbIDGood,
  kTrkNonAmbIDBad,
  kTrkAmbID,
  kTrkAmbIDGood,
  kTrkAmbIDBad,
  kTrkNonAmbIDExtra,
  kTrkNonAmbIDExtraGood,
  kTrkNonAmbIDExtraBad,
  kReAssoc,
  kReAssocGood,
  kReAssocGoodIsCompTrue,
  kReAssocGoodIsCompFalse,
  kReAssocBad,
  kReAssocBadIsCompTrue,
  kReAssocBadIsCompFalse,
  kAssoc,
  kAssocGood,
  kAssocGoodIsCompTrue,
  kAssocGoodIsCompFalse,
  kAssocBad,
  kAssocBadIsCompTrue,
  kAssocBadIsCompFalse,
  nReAssocMCTrackStatusCheck
};

enum class HistStatusReAssocVtx {
  kTrkRec = 0,
  kTrkIdGt0,
  kTrkNonAmb,
  kTrkNonAmbGood,
  kTrkNonAmbBad,
  kTrkNonAmbID,
  kTrkNonAmbIDGood,
  kTrkNonAmbIDBad,
  kTrkAmbID,
  kTrkAmbIDGood,
  kTrkAmbIDBad,
  kTrkNonAmbIDExtra,
  kTrkNonAmbIDExtraGood,
  kTrkNonAmbIDExtraBad,
  kReAssoc,
  kReAssocGood,
  kReAssocGoodIsCompTrue,
  kReAssocGoodIsCompFalse,
  kReAssocBad,
  kReAssocBadIsCompTrue,
  kReAssocBadIsCompFalse,
  kAssoc,
  kAssocGood,
  kAssocGoodIsCompTrue,
  kAssocGoodIsCompFalse,
  kAssocBad,
  kAssocBadIsCompTrue,
  kAssocBadIsCompFalse,
  nHistStatusReAssocVtx
};

struct QaReassocMFTPbPb {

  std::array<std::shared_ptr<THnSparse>, 4> hCollAssoc;
  std::array<std::shared_ptr<THnSparse>, 28> hReAssocVtxRes;
  std::array<std::shared_ptr<THnSparse>, 28> hReAssocDCA;
  std::array<std::shared_ptr<THnSparse>, 21> hTimeAssocWithReassocMC;

  enum OccupancyEst { TrkITS = 1,
                      Ft0C };

  HistogramRegistry registry{
    "registry",
    {},
    OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry qaregistry{
    "qaregistry",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    false,
    true};

  struct : ConfigurableGroup {
    Configurable<bool> cfgUseTrackSel{"cfgUseTrackSel", false, "Flag to apply track selection"};
    Configurable<bool> cfgUseParticleSel{"cfgUseParticleSel", true, "Flag to apply particle selection"};
    Configurable<bool> cfgUsePrimaries{"cfgUsePrimaries", true, "Select primary particles"};
    Configurable<bool> cfgUseSecondaries{"cfgUseSecondaries", false, "Select secondary particles"};
    Configurable<bool> cfgRemoveTrivialAssoc{"cfgRemoveTrivialAssoc", false, "Skip trivial associations"};
    Configurable<bool> cfgRemoveAmbiguousTracks{"cfgRemoveAmbiguousTracks", false, "Remove ambiguous tracks"};
    Configurable<bool> cfgRemoveOrphanTracks{"cfgRemoveOrphanTracks", true, "Remove orphan tracks"};
    Configurable<bool> cfgRemoveReassigned{"cfgRemoveReassigned", false, "Remove reassgined tracks"};
    Configurable<bool> cfgRemoveSplitVertex{"cfgRemoveSplitVertex", true, "Remove split vertices"};
    Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  } gConf;

  struct : ConfigurableGroup {
    ConfigurableAxis occupancyBins{"occupancyBins", {VARIABLE_WIDTH, 0.0f, 250.0f, 500.0f, 750.0f, 1000.0f, 1500.0f, 2000.0f, 3000.0f, 4500.0f, 6000.0f, 8000.0f, 10000.0f, 50000.0f}, "Occupancy"};
    ConfigurableAxis centralityBins{"centralityBins", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100}, "Centrality"};
    ConfigurableAxis ptBins{"ptBins", {101, -0.5, 10.5}, "pT binning (GeV/c)"};
    ConfigurableAxis multBins{"multBins", {1001, -0.5, 1000.5}, "Multiplicity binning"};
    ConfigurableAxis deltaZBins{"deltaZBins", {800, -10., 10.}, "Delta Z-vtx binning (cm)"};
    ConfigurableAxis dcaXYBins{"dcaXYBins", {800, -1., 1.}, "DCAxy binning (cm)"};
    ConfigurableAxis dcaZBins{"dcaZBins", {800, -1., 1.}, "DCAz binning (cm)"};
    ConfigurableAxis phiBins{"phiBins", {629, 0., TwoPI}, "#varphi binning (rad)"};
    ConfigurableAxis etaBins{"etaBins", {20, -4., -2.}, "#eta binning"};
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
  } eventCuts;

  Service<o2::framework::O2DatabasePDG> pdg{};
  Service<ccdb::BasicCCDBManager> ccdb{};
  Configurable<int64_t> ccdbNoLaterThan{"ccdbNoLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<bool> cfgApplyZShiftFromCCDB{"cfgApplyZShiftFromCCDB", false, "flag to apply z shift from CCDB"};
  Configurable<std::string> cfgZShiftPath{"cfgZShiftPath", "Users/m/mcoquet/ZShift", "CCDB path for z shift to apply to forward tracks"};
  Configurable<float> cfgManualZShift{"cfgManualZShift", 0.0f, "manual z-shift for propagation of global muon to PV"};

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
  float mZShift = 0;                            // z-vertex shift

  o2::parameters::GRPMagField* grpmag = nullptr;

  /// @brief init function, definition of histograms
  void init(InitContext&)
  {
    const AxisSpec centralityAxis = {binOpt.centralityBins, "Centrality", "centrality axis"};
    const AxisSpec occupancyAxis = {binOpt.occupancyBins, "Occupancy", "occupancy axis"};
    const AxisSpec ptAxis = {binOpt.ptBins, "Pt axis (GeV/c)"};
    const AxisSpec multAxis = {binOpt.multBins, "N_{trk} axis"};
    const AxisSpec deltaZAxis = {binOpt.deltaZBins, "Delta Z-vtx axis"};
    const AxisSpec dcaxyAxis = {binOpt.dcaXYBins, "DCA-xy axis"};
    const AxisSpec dcazAxis = {binOpt.dcaZBins, "DCA-z axis"};
    const AxisSpec etaAxis = {binOpt.etaBins, "#eta axis"};

    rctChecker.init(rctCuts.cfgEvtRCTFlagCheckerLabel, rctCuts.cfgEvtRCTFlagCheckerZDCCheck, rctCuts.cfgEvtRCTFlagCheckerLimitAcceptAsBad);

    ccdb->setURL(ccdbUrl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(ccdbNoLaterThan.value);
    ccdb->setFatalWhenNull(false);

    qaregistry.add("Events/hNchTVX", "; status;", HistType::kTH2F, {multAxis, {2, 0, 2}});

    registry.add("Events/hEvtSel", "Number of events; Cut; #Evt Passed Cut", {HistType::kTH1F, {{static_cast<int>(EvtSel::nEvtSel), -0.5, +static_cast<int>(EvtSel::nEvtSel) - 0.5}}});
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
    registry.get<TH1>(HIST("Events/hEvtSel"))->SetMinimum(0.1);
    for (int iBin = 0; iBin < static_cast<int>(EvtSel::nEvtSel); iBin++) {
      registry.get<TH1>(HIST("Events/hEvtSel"))->GetXaxis()->SetBinLabel(iBin + 1, labelEvtSel[iBin].data());
    }

    registry.add("Tracks/hBestTrkSel", "Number of best tracks; Cut; #Tracks Passed Cut", {HistType::kTH1F, {{static_cast<int>(TrkTrkBestSel::nTrkTrkBestSel), -0.5, +static_cast<int>(TrkTrkBestSel::nTrkTrkBestSel) - 0.5}}});
    std::array<std::string_view, static_cast<int>(TrkTrkBestSel::nTrkTrkBestSel)> labelTrkTrkBestSel{
      "All",
      "Assigned (ID>=0)",
      "No orphans",
      "DCA xy cut",
      "DCA z cut",
      "#Reassoc"};
    registry.get<TH1>(HIST("Tracks/hBestTrkSel"))->SetMinimum(0.1);
    for (int iBin = 0; iBin < static_cast<int>(TrkTrkBestSel::nTrkTrkBestSel); iBin++) {
      registry.get<TH1>(HIST("Tracks/hBestTrkSel"))->GetXaxis()->SetBinLabel(iBin + 1, labelTrkTrkBestSel[iBin].data());
    }

    registry.add("Tracks/hTrkSel", "Number of tracks; Cut; #Tracks Passed Cut", {HistType::kTH1F, {{static_cast<int>(TrkSel::nTrkSel), -0.5, +static_cast<int>(TrkSel::nTrkSel) - 0.5}}});
    std::array<std::string_view, static_cast<int>(TrkSel::nTrkSel)> labelTrkSel{
      "All",
      "Ncls",
      "Chi2",
      "Eta",
      "Phi cut",
      "Pt",
      "CA"};
    registry.get<TH1>(HIST("Tracks/hTrkSel"))->SetMinimum(0.1);
    for (int iBin = 0; iBin < static_cast<int>(TrkSel::nTrkSel); iBin++) {
      registry.get<TH1>(HIST("Tracks/hTrkSel"))->GetXaxis()->SetBinLabel(iBin + 1, labelTrkSel[iBin].data());
    }

    if (doprocessTimeAssocMC) {
      registry.add("TimeAssocMC/hTimeAssocMCEventStatus", ";status", {HistType::kTH1F, {{static_cast<int>(ReAssocMCEventStatus::nEvtReAsReAssocMCEventStatus), -0.5, +static_cast<int>(ReAssocMCEventStatus::nEvtReAsReAssocMCEventStatus) - 0.5}}});
      std::array<std::string_view, static_cast<int>(ReAssocMCEventStatus::nEvtReAsReAssocMCEventStatus)> labelTimeAssocMCEventStatus{
        "All",
        "Selected",
        "Has Mc Coll",
        "Split Vtx Removed",
        "Vtx-z cut MC"};
      registry.get<TH1>(HIST("TimeAssocMC/hTimeAssocMCEventStatus"))->SetMinimum(0.1);
      for (int iBin = 0; iBin < static_cast<int>(ReAssocMCEventStatus::nEvtReAsReAssocMCEventStatus); iBin++) {
        registry.get<TH1>(HIST("TimeAssocMC/hTimeAssocMCEventStatus"))->GetXaxis()->SetBinLabel(iBin + 1, labelTimeAssocMCEventStatus[iBin].data());
      }

      registry.add("TimeAssocMC/VtxStatus", ";status", {HistType::kTH1F, {{2, 0.5, 2.5}}});
      auto hstat = registry.get<TH1>(HIST("TimeAssocMC/VtxStatus"));
      hstat->GetXaxis()->SetBinLabel(1, "Good vtx");
      hstat->GetXaxis()->SetBinLabel(2, "Wrong vtx");

      registry.add({"TimeAssocMC/hVertexResV1", "; #Delta X (cm); #Delta Y (cm); #Delta Z (cm)", {HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis}}});
      registry.add({"TimeAssocMC/hVertexResV2", "; #Delta X (cm); #Delta Y (cm); #Delta Z (cm)", {HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis}}});
      registry.add({"TimeAssocMC/hVTXkSelGoodVtxTrue", "; #Delta X (cm); #Delta Y (cm); #Delta Z (cm)", {HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis}}});
      registry.add({"TimeAssocMC/hVTXkSelBadVtxTrue", "; #Delta X (cm); #Delta Y (cm); #Delta Z (cm)", {HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis}}});
      registry.add({"TimeAssocMC/hVTXkSelNonAmb", "; #Delta X (cm); #Delta Y (cm); #Delta Z (cm)", {HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis}}});
      registry.add({"TimeAssocMC/hVTXkSelNonAmbGoodVtx", "; #Delta X (cm); #Delta Y (cm); #Delta Z (cm)", {HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis}}});
      registry.add({"TimeAssocMC/hVTXkSelNonAmbGoodVtxTrue", "; #Delta X (cm); #Delta Y (cm); #Delta Z (cm)", {HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis}}});
      registry.add({"TimeAssocMC/hVTXkSelNonAmbBadVtxTrue", "; #Delta X (cm); #Delta Y (cm); #Delta Z (cm)", {HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis}}});
      registry.add({"TimeAssocMC/hVTXkSelNonAmbID", "; #Delta X (cm); #Delta Y (cm); #Delta Z (cm)", {HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis}}});
      registry.add({"TimeAssocMC/hVTXkSelNonAmbIDGoodVtxTrue", "; #Delta X (cm); #Delta Y (cm); #Delta Z (cm)", {HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis}}});
      registry.add({"TimeAssocMC/hVTXkSelNonAmbIDBadVtxTrue", "; #Delta X (cm); #Delta Y (cm); #Delta Z (cm)", {HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis}}});
      registry.add({"TimeAssocMC/hVTXkSelAmbID", "; #Delta X (cm); #Delta Y (cm); #Delta Z (cm)", {HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis}}});
      registry.add({"TimeAssocMC/hVTXkSelAmbIDGoodVtxTrue", "; #Delta X (cm); #Delta Y (cm); #Delta Z (cm)", {HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis}}});
      registry.add({"TimeAssocMC/hVTXkSelAmbIDBadVtxTrue", "; #Delta X (cm); #Delta Y (cm); #Delta Z (cm)", {HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis}}});
      registry.add({"TimeAssocMC/hVTXkSelNonAmbIDExtra", "; #Delta X (cm); #Delta Y (cm); #Delta Z (cm)", {HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis}}});
      registry.add({"TimeAssocMC/hVTXkSelNonAmbIDExtraGoodVtxTrue", "; #Delta X (cm); #Delta Y (cm); #Delta Z (cm)", {HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis}}});
      registry.add({"TimeAssocMC/hVTXkSelNonAmbIDExtraBadVtxTrue", "; #Delta X (cm); #Delta Y (cm); #Delta Z (cm)", {HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis}}});
      registry.add({"TimeAssocMC/hVTXkSelAmb", "; #Delta X (cm); #Delta Y (cm); #Delta Z (cm)", {HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis}}});
      registry.add({"TimeAssocMC/hVTXkSelAmbGoodVtxTrue", "; #Delta X (cm); #Delta Y (cm); #Delta Z (cm)", {HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis}}});
      registry.add({"TimeAssocMC/hVTXkSelAmbBadVtxTrue", "; #Delta X (cm); #Delta Y (cm); #Delta Z (cm)", {HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis}}});

      registry.add("TimeAssocMC/hTimeAssocCheckVtxType", ";status", {HistType::kTH1F, {{static_cast<int>(AssocCheckVtxType::nAssocVtxType), -0.5, +static_cast<int>(AssocCheckVtxType::nAssocVtxType) - 0.5}}});
      std::array<std::string_view, static_cast<int>(AssocCheckVtxType::nAssocVtxType)> labelAssocVtxType{
        "kAllVtx=True",
        "kAllVtx=False",
        "kAllGoodVtx",
        "kAllGoodVtx=True",
        "kAllGoodVtx=False",
        "kAllBadVtx",
        "kAllBadVtx=True",
        "kAllBadVtx=False"};
      registry.get<TH1>(HIST("TimeAssocMC/hTimeAssocCheckVtxType"))->SetMinimum(0.1);
      for (int iBin = 0; iBin < static_cast<int>(AssocCheckVtxType::nAssocVtxType); iBin++) {
        registry.get<TH1>(HIST("TimeAssocMC/hTimeAssocCheckVtxType"))->GetXaxis()->SetBinLabel(iBin + 1, labelAssocVtxType[iBin].data());
      }

      registry.add("TimeAssocMC/hAmbTrackType", ";status", {HistType::kTH1F, {{static_cast<int>(AmbTrkType::nAmbTrkType), -0.5, +static_cast<int>(AmbTrkType::nAmbTrkType) - 0.5}}});
      std::array<std::string_view, static_cast<int>(AmbTrkType::nAmbTrkType)> labelAmbiguity{
        "all",
        "orphan",
        "orphanNull",
        "nonAmb",
        "nonAmbSame",
        "Amb",
        "AmbGt1"};
      registry.get<TH1>(HIST("TimeAssocMC/hAmbTrackType"))->SetMinimum(0.1);
      for (int iBin = 0; iBin < static_cast<int>(AmbTrkType::nAmbTrkType); iBin++) {
        registry.get<TH1>(HIST("TimeAssocMC/hAmbTrackType"))->GetXaxis()->SetBinLabel(iBin + 1, labelAmbiguity[iBin].data());
      }

      registry.add("TimeAssocMC/hAmbTrkTypeAssocFlag", ";status", {HistType::kTH1F, {{static_cast<int>(AmbTrkTypeAssocFlag::nSelAmbTrkTypeAssocFlag), -0.5, +static_cast<int>(AmbTrkTypeAssocFlag::nSelAmbTrkTypeAssocFlag) - 0.5}}});
      std::array<std::string_view, static_cast<int>(AmbTrkTypeAssocFlag::nSelAmbTrkTypeAssocFlag)> lAmbTrackType{
        "all sel",
        "all sel good vtx",
        "all sel good vtx true",
        "all sel bad vtx",
        "all sel bad vtx true",
        "all non-amb",
        "non-amb good vtx",
        "non-amb good vtx true",
        "non-amb bad vtx",
        "non-amb bad vtx true",
        "non-amb id",
        "non-amb id good vtx true",
        "non-amb bad vtx true",
        "amb id",
        "amb id good vtx true",
        "amb id bad vtx true",
        "non-amb id ext",
        "non-amb id ext good vtx true",
        "non-amb id ext bad vtx true",
        "amb all",
        "amb good vtx true",
        "amb bad vtx true",
        "orhpan null"};
      registry.get<TH1>(HIST("TimeAssocMC/hAmbTrkTypeAssocFlag"))->SetMinimum(0.1);
      for (int iBin = 0; iBin < static_cast<int>(AmbTrkTypeAssocFlag::nSelAmbTrkTypeAssocFlag); iBin++) {
        registry.get<TH1>(HIST("TimeAssocMC/hAmbTrkTypeAssocFlag"))->GetXaxis()->SetBinLabel(iBin + 1, lAmbTrackType[iBin].data());
      }
    }

    if (doprocessTimeAssocWithReassocMC) {
      registry.add("TimeAssocWithReassocMC/VtxStatus", ";status", {HistType::kTH1F, {{2, 0.5, 2.5}}});
      auto hstat = registry.get<TH1>(HIST("TimeAssocWithReassocMC/VtxStatus"));
      hstat->GetXaxis()->SetBinLabel(1, "Good vtx");
      hstat->GetXaxis()->SetBinLabel(2, "Wrong vtx");

      registry.add("TimeAssocWithReassocMC/hReassocCheckVtxType", ";status", {HistType::kTH1F, {{static_cast<int>(ReassocCheckVtxType::nReassocVtxType), -0.5, +static_cast<int>(ReassocCheckVtxType::nReassocVtxType) - 0.5}}});
      std::array<std::string_view, static_cast<int>(ReassocCheckVtxType::nReassocVtxType)> labelReAssocVtxType{
        "kIsTrueVtxAll=True",
        "kIsTrueVtxAll=False",
        "IsRecGoodAll=True",
        "kIsRecGoodAll=False",
        "kIsRecGoodMatchAll=True",
        "kIsRecGoodMatchAll=False",
        "kIsTrueVtxVsGoodVtx=True",
        "kIsTrueVtxVsGoodVtx=False",
        "kIsRecGoodVsGoodVtx=True",
        "kIsRecGoodVsGoodVtx=False",
        "kIsRecGoodMatchVsGoodVtx=True",
        "kIsRecGoodMatchVsGoodVtx=False",
        "kIsTrueVtxVsBadVtx=True",
        "kIsTrueVtxVsBadVtx=False",
        "kIsRecGoodVsBadVtx=True",
        "kIsRecGoodVsBadVtx=False",
        "kIsRecGoodMatchVsBadVtx=True",
        "kIsRecGoodMatchVsBadVtx=False"};
      registry.get<TH1>(HIST("TimeAssocWithReassocMC/hReassocCheckVtxType"))->SetMinimum(0.1);
      for (int iBin = 0; iBin < static_cast<int>(ReassocCheckVtxType::nReassocVtxType); iBin++) {
        registry.get<TH1>(HIST("TimeAssocWithReassocMC/hReassocCheckVtxType"))->GetXaxis()->SetBinLabel(iBin + 1, labelReAssocVtxType[iBin].data());
      }

      hTimeAssocWithReassocMC[0] = registry.add<THnSparse>("TimeAssocWithReassocMC/hDCAReassocOrig", ";#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis});
      hTimeAssocWithReassocMC[1] = registry.add<THnSparse>("TimeAssocWithReassocMC/hDCAReassocBest", ";#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis});
      hTimeAssocWithReassocMC[2] = registry.add<THnSparse>("TimeAssocWithReassocMC/hDCAReassocTruth", ";#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis});

      hTimeAssocWithReassocMC[3] = registry.add<THnSparse>("TimeAssocWithReassocMC/hDCAReassocOrigTrueVtx", ";#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis});
      hTimeAssocWithReassocMC[4] = registry.add<THnSparse>("TimeAssocWithReassocMC/hDCAReassocOrigRecGood", ";#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis});
      hTimeAssocWithReassocMC[5] = registry.add<THnSparse>("TimeAssocWithReassocMC/hDCAReassocOrigRecGoodMatch", ";#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis});

      hTimeAssocWithReassocMC[6] = registry.add<THnSparse>("TimeAssocWithReassocMC/hDCAReassocOrigTrueVtxVtxFlagGood", ";#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis});
      hTimeAssocWithReassocMC[7] = registry.add<THnSparse>("TimeAssocWithReassocMC/hDCAReassocOrigRecGoodVtxFlagGood", ";#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis});
      hTimeAssocWithReassocMC[8] = registry.add<THnSparse>("TimeAssocWithReassocMC/hDCAReassocOrigRecGoodMatchVtxFlagGood", ";#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis});

      hTimeAssocWithReassocMC[9] = registry.add<THnSparse>("TimeAssocWithReassocMC/hDCAReassocOrigTrueVtxVtxFlagBad", ";#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis});
      hTimeAssocWithReassocMC[10] = registry.add<THnSparse>("TimeAssocWithReassocMC/hDCAReassocOrigRecGoodVtxFlagBad", ";#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis});
      hTimeAssocWithReassocMC[11] = registry.add<THnSparse>("TimeAssocWithReassocMC/hDCAReassocOrigRecGoodMatchVtxFlagBad", ";#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis});

      hTimeAssocWithReassocMC[12] = registry.add<THnSparse>("TimeAssocWithReassocMC/hDCAReassocBestTrueVtx", ";#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis});
      hTimeAssocWithReassocMC[13] = registry.add<THnSparse>("TimeAssocWithReassocMC/hDCAReassocBestRecGood", ";#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis});
      hTimeAssocWithReassocMC[14] = registry.add<THnSparse>("TimeAssocWithReassocMC/hDCAReassocBestRecGoodMatch", ";#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis});

      hTimeAssocWithReassocMC[15] = registry.add<THnSparse>("TimeAssocWithReassocMC/hDCAReassocBestTrueVtxVtxFlagGood", ";#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis});
      hTimeAssocWithReassocMC[16] = registry.add<THnSparse>("TimeAssocWithReassocMC/hDCAReassocBestRecGoodVtxFlagGood", ";#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis});
      hTimeAssocWithReassocMC[17] = registry.add<THnSparse>("TimeAssocWithReassocMC/hDCAReassocBestRecGoodMatchVtxFlagGood", ";#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis});

      hTimeAssocWithReassocMC[18] = registry.add<THnSparse>("TimeAssocWithReassocMC/hDCAReassocBestTrueVtxVtxFlagBad", ";#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis});
      hTimeAssocWithReassocMC[19] = registry.add<THnSparse>("TimeAssocWithReassocMC/hDCAReassocBestRecGoodVtxFlagBad", ";#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis});
      hTimeAssocWithReassocMC[20] = registry.add<THnSparse>("TimeAssocWithReassocMC/hDCAReassocBestRecGoodMatchVtxFlagBad", ";#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis});
    }

    if (doprocessAssocMC) {
      registry.add("TrackToColl/hAmbTrackType", "hAmbTrackType", {HistType::kTH1F, {{static_cast<int>(AmbTrkType::nAmbTrkType), -0.5, +static_cast<int>(AmbTrkType::nAmbTrkType) - 0.5}}});
      std::array<std::string_view, static_cast<int>(AmbTrkType::nAmbTrkType)> labelAmbiguity{
        "all",
        "orphan",
        "orphanNull",
        "nonAmb",
        "nonAmbSame",
        "Amb",
        "AmbGt1"};
      registry.get<TH1>(HIST("TrackToColl/hAmbTrackType"))->SetMinimum(0.1);
      for (int iBin = 0; iBin < static_cast<int>(AmbTrkType::nAmbTrkType); iBin++) {
        registry.get<TH1>(HIST("TrackToColl/hAmbTrackType"))->GetXaxis()->SetBinLabel(iBin + 1, labelAmbiguity[iBin].data());
      }

      // tracks not associated to any collision
      hCollAssoc[0] = registry.add<THnSparse>("TrackToColl/hNonAssocTracks", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis});
      // tracks associasted to a collision
      hCollAssoc[1] = registry.add<THnSparse>("TrackToColl/hAssocTracks", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis});
      // tracks associated to the correct collision considering only first reco collision (based on the MC collision index)
      hCollAssoc[2] = registry.add<THnSparse>("TrackToColl/hGoodAssocTracks", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis});
      // tracks associated to the correct collision considering all ambiguous reco collisions (based on the MC collision index)
      hCollAssoc[3] = registry.add<THnSparse>("TrackToColl/hGoodAssocTracksAmb", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis});

      registry.add("TrackToColl/histFracTracksFakeMcColl", "Fraction of tracks originating from fake collision; fraction; entries", {HistType::kTH1F, {{101, 0., 1.01}}});
      registry.add("TrackToColl/histFracGoodTracks", "Fraction of tracks originating from the correct collision; fraction; entries", {HistType::kTH1F, {{101, 0., 1.01}}});
      registry.add("TrackToColl/histAmbTrackNumColls", "Number of collisions associated to an ambiguous track; no. collisions; entries", {HistType::kTH1F, {{30, -0.5, 29.5}}});
      registry.add("TrackToColl/histTrackNumColls", "Number of collisions associated to track; no. collisions; entries", {HistType::kTH1F, {{30, -0.5, 29.5}}});
      registry.add("TrackToColl/histNonAmbTrackNumColls", "Number of collisions associated to non-ambiguous track; no. collisions; entries", {HistType::kTH1F, {{30, -0.5, 29.5}}});
      registry.add("TrackToColl/histAmbTrackZvtxRMS", "RMS of #it{Z}^{reco} of collisions associated to a track; RMS(#it{Z}^{reco}) (cm); entries", {HistType::kTH1F, {{100, 0., 10.}}});
    }

    if (doprocessReAssoc3dMC) {
      registry.add({"ReAssocMC/EvtGenRecReassoc", ";status;centrality;occupancy", {HistType::kTHnSparseF, {{3, 0.5, 3.5}, centralityAxis, occupancyAxis}}});
      auto heff = registry.get<THnSparse>(HIST("ReAssocMC/EvtGenRecReassoc"));
      heff->GetAxis(0)->SetBinLabel(1, "All generated");
      heff->GetAxis(0)->SetBinLabel(2, "All reconstructed");
      heff->GetAxis(0)->SetBinLabel(3, "Selected reconstructed");

      registry.add({"ReAssocMC/hReAssocMCEventStatus", ";status;centrality;occupancy", {HistType::kTHnSparseF, {{static_cast<int>(ReAssocMCEventStatus::nEvtReAsReAssocMCEventStatus), -0.5, +static_cast<int>(ReAssocMCEventStatus::nEvtReAsReAssocMCEventStatus) - 0.5}, centralityAxis, occupancyAxis}}});
      std::array<std::string_view, static_cast<int>(ReAssocMCEventStatus::nEvtReAsReAssocMCEventStatus)> labelReAssocMCEventStatus{
        "All",
        "Selected",
        "Has Mc Coll",
        "Split Vtx Removed",
        "Vtx-z cut MC"};
      // registry.get<THnSparse>(HIST("ReAssocMC/hReAssocMCEventStatus"))->SetMinimum(0.1);
      for (int iBin = 0; iBin < static_cast<int>(ReAssocMCEventStatus::nEvtReAsReAssocMCEventStatus); iBin++) {
        registry.get<THnSparse>(HIST("ReAssocMC/hReAssocMCEventStatus"))->GetAxis(0)->SetBinLabel(iBin + 1, labelReAssocMCEventStatus[iBin].data());
      }

      registry.add({"ReAssocMC/hReAssocMCTrackStatus", ";status;centrality;occupancy", {HistType::kTHnSparseF, {{static_cast<int>(ReAssocMCTrackStatus::nReAssocMCTrackStatusCheck), -0.5, +static_cast<int>(ReAssocMCTrackStatus::nReAssocMCTrackStatusCheck) - 0.5}, centralityAxis, occupancyAxis}}});
      std::array<std::string_view, static_cast<int>(ReAssocMCTrackStatus::nReAssocMCTrackStatusCheck)> labelReAssocMCTrackStatus{
        "All",
        "Best sel",
        "Trk sel",
        "Has coll",
        "Reas rm",
        "Has part",
        "Trk idGt0",
        "Non-amb",
        "Non-amb good coll.",
        "Non-amb bad coll.",
        "Non-amb id",
        "Non-amb id good",
        "Non-amb id bad",
        "Amb id",
        "Amb id good",
        "Amb id bad",
        "Non-amb id ex",
        "Non-amb id ex good",
        "Non-amb id ex bad",
        "ReAssoc",
        "ReAssoc good",
        "ReAssoc good Comp True",
        "ReAssoc good Comp False",
        "ReAssoc bad",
        "ReAssoc bad Comp True",
        "ReAssoc bad Comp False",
        "Assoc (gt1 amb)",
        "Assoc good",
        "Assoc good Comp True",
        "Assoc good Comp False",
        "Assoc bad",
        "Assoc bad Comp True",
        "Assoc bad Comp False"};
      // registry.get<THnSparse>(HIST("ReAssocMC/hReAssocMCTrackStatus"))->SetMinimum(0.1);
      for (int iBin = 0; iBin < static_cast<int>(ReAssocMCTrackStatus::nReAssocMCTrackStatusCheck); iBin++) {
        registry.get<THnSparse>(HIST("ReAssocMC/hReAssocMCTrackStatus"))->GetAxis(0)->SetBinLabel(iBin + 1, labelReAssocMCTrackStatus[iBin].data());
      }

      // Vertex resolution
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkRec)] = registry.add<THnSparse>("ReAssocMC/hVtxResRec", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis, centralityAxis, occupancyAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkIdGt0)] = registry.add<THnSparse>("ReAssocMC/hVtxResIdGt0", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis, centralityAxis, occupancyAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmb)] = registry.add<THnSparse>("ReAssocMC/hVtxResNonAmb", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis, centralityAxis, occupancyAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbGood)] = registry.add<THnSparse>("ReAssocMC/hVtxResNonAmbGood", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis, centralityAxis, occupancyAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbBad)] = registry.add<THnSparse>("ReAssocMC/hVtxResNonAmbBad", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis, centralityAxis, occupancyAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbID)] = registry.add<THnSparse>("ReAssocMC/hVtxResNonAmbID", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis, centralityAxis, occupancyAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbIDGood)] = registry.add<THnSparse>("ReAssocMC/hVtxResNonAmbIDGood", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis, centralityAxis, occupancyAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbIDBad)] = registry.add<THnSparse>("ReAssocMC/hVtxResNonAmbIDBad", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis, centralityAxis, occupancyAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkAmbID)] = registry.add<THnSparse>("ReAssocMC/hVtxResAmbID", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis, centralityAxis, occupancyAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkAmbIDGood)] = registry.add<THnSparse>("ReAssocMC/hVtxResAmbIDGood", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis, centralityAxis, occupancyAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkAmbIDBad)] = registry.add<THnSparse>("ReAssocMC/hVtxResAmbIDBad", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis, centralityAxis, occupancyAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbIDExtra)] = registry.add<THnSparse>("ReAssocMC/hVtxResNonAmbIDExtra", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis, centralityAxis, occupancyAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbIDExtraGood)] = registry.add<THnSparse>("ReAssocMC/hVtxResNonAmbIDExtraGood", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis, centralityAxis, occupancyAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbIDExtraBad)] = registry.add<THnSparse>("ReAssocMC/hVtxResNonAmbIDExtraBad", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis, centralityAxis, occupancyAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kReAssoc)] = registry.add<THnSparse>("ReAssocMC/hVtxResReAssoc", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis, centralityAxis, occupancyAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kReAssocGood)] = registry.add<THnSparse>("ReAssocMC/hVtxResReAssocGood", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis, centralityAxis, occupancyAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kReAssocGoodIsCompTrue)] = registry.add<THnSparse>("ReAssocMC/hVtxResReAssocGoodIsCompTrue", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis, centralityAxis, occupancyAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kReAssocGoodIsCompFalse)] = registry.add<THnSparse>("ReAssocMC/hVtxResReAssocGoodIsCompFalse", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis, centralityAxis, occupancyAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kReAssocBad)] = registry.add<THnSparse>("ReAssocMC/hVtxResReAssocBad", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis, centralityAxis, occupancyAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kReAssocBadIsCompTrue)] = registry.add<THnSparse>("ReAssocMC/hVtxResReAssocBadIsCompTrue", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis, centralityAxis, occupancyAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kReAssocBadIsCompFalse)] = registry.add<THnSparse>("ReAssocMC/hVtxResReAssocBadIsCompFalse", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis, centralityAxis, occupancyAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kAssoc)] = registry.add<THnSparse>("ReAssocMC/hVtxResAssoc", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis, centralityAxis, occupancyAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kAssocGood)] = registry.add<THnSparse>("ReAssocMC/hVtxResAssocGood", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis, centralityAxis, occupancyAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kAssocGoodIsCompTrue)] = registry.add<THnSparse>("ReAssocMC/hVtxResAssocGoodIsCompTrue", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis, centralityAxis, occupancyAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kAssocGoodIsCompFalse)] = registry.add<THnSparse>("ReAssocMC/hVtxResAssocGoodIsCompFalse", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis, centralityAxis, occupancyAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kAssocBad)] = registry.add<THnSparse>("ReAssocMC/hVtxResAssocBad", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis, centralityAxis, occupancyAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kAssocBadIsCompTrue)] = registry.add<THnSparse>("ReAssocMC/hVtxResAssocBadIsCompTrue", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis, centralityAxis, occupancyAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kAssocBadIsCompFalse)] = registry.add<THnSparse>("ReAssocMC/hVtxResAssocBadIsCompFalse", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis, centralityAxis, occupancyAxis});

      // DCA
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkRec)] = registry.add<THnSparse>("ReAssocMC/hDCARec", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis, centralityAxis, occupancyAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkIdGt0)] = registry.add<THnSparse>("ReAssocMC/hDCAIdGt0", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis, centralityAxis, occupancyAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmb)] = registry.add<THnSparse>("ReAssocMC/hDCANonAmbAll", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis, centralityAxis, occupancyAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbGood)] = registry.add<THnSparse>("ReAssocMC/hDCANonAmbGood", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis, centralityAxis, occupancyAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbBad)] = registry.add<THnSparse>("ReAssocMC/hDCANonAmbBad", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis, centralityAxis, occupancyAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbID)] = registry.add<THnSparse>("ReAssocMC/hDCAAmbAll", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis, centralityAxis, occupancyAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbIDGood)] = registry.add<THnSparse>("ReAssocMC/hDCAAmbGood", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis, centralityAxis, occupancyAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbIDBad)] = registry.add<THnSparse>("ReAssocMC/hDCAAmbBad", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis, centralityAxis, occupancyAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkAmbID)] = registry.add<THnSparse>("ReAssocMC/hDCANonAmbAllE", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis, centralityAxis, occupancyAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkAmbIDGood)] = registry.add<THnSparse>("ReAssocMC/hDCANonAmbGoodE", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis, centralityAxis, occupancyAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkAmbIDBad)] = registry.add<THnSparse>("ReAssocMC/hDCANonAmbBadE", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis, centralityAxis, occupancyAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbIDExtra)] = registry.add<THnSparse>("ReAssocMC/hDCANonAmbIDExtra", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis, centralityAxis, occupancyAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbIDExtraGood)] = registry.add<THnSparse>("ReAssocMC/hDCANonAmbIDExtraGood", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis, centralityAxis, occupancyAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbIDExtraBad)] = registry.add<THnSparse>("ReAssocMC/hDCANonAmbIDExtraBad", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis, centralityAxis, occupancyAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kReAssoc)] = registry.add<THnSparse>("ReAssocMC/hDCAReAssoc", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis, centralityAxis, occupancyAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kReAssocGood)] = registry.add<THnSparse>("ReAssocMC/hDCAReAssocGood", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis, centralityAxis, occupancyAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kReAssocGoodIsCompTrue)] = registry.add<THnSparse>("ReAssocMC/hDCAReAssocGoodIsCompTrue", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis, centralityAxis, occupancyAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kReAssocGoodIsCompFalse)] = registry.add<THnSparse>("ReAssocMC/hDCAReAssocGoodIsCompFalse", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis, centralityAxis, occupancyAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kReAssocBad)] = registry.add<THnSparse>("ReAssocMC/hDCAReAssocBad", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis, centralityAxis, occupancyAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kReAssocBadIsCompTrue)] = registry.add<THnSparse>("ReAssocMC/hDCAReAssocBadIsCompTrue", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis, centralityAxis, occupancyAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kReAssocBadIsCompFalse)] = registry.add<THnSparse>("ReAssocMC/hDCAReAssocBadIsCompFalse", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis, centralityAxis, occupancyAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kAssoc)] = registry.add<THnSparse>("ReAssocMC/hDCAAssoc", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis, centralityAxis, occupancyAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kAssocGood)] = registry.add<THnSparse>("ReAssocMC/hDCAAssocGood", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis, centralityAxis, occupancyAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kAssocGoodIsCompTrue)] = registry.add<THnSparse>("ReAssocMC/hDCAAssocGoodIsCompTrue", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis, centralityAxis, occupancyAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kAssocGoodIsCompFalse)] = registry.add<THnSparse>("ReAssocMC/hDCAAssocGoodIsCompFalse", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis, centralityAxis, occupancyAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kAssocBad)] = registry.add<THnSparse>("ReAssocMC/hDCAAssocBad", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis, centralityAxis, occupancyAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kAssocBadIsCompTrue)] = registry.add<THnSparse>("ReAssocMC/hDCAAssocBadIsCompTrue", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis, centralityAxis, occupancyAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kAssocBadIsCompFalse)] = registry.add<THnSparse>("ReAssocMC/hDCAAssocBadIsCompFalse", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm); centrality; occupancy", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis, centralityAxis, occupancyAxis});
    }
  }

  using ExtBCs = soa::Join<aod::BCs, aod::Timestamps, aod::MatchedBCCollisionsSparseMulti>;
  using CollsGenCentFT0CExtra = soa::Join<aod::McCollisionLabels, aod::Collisions, aod::CentFT0Cs, aod::EvSels, aod::MultMCExtras, aod::McCollsExtra>;
  using MftTracksWCollsMC = soa::Join<aod::MFTTracks, aod::MFTTrkCompColls, aod::McMFTTrackLabels>;
  using BestTracks3dWCollsMC = soa::Join<aod::MFTTrkCompColls, aod::BestCollisionsFwd3d, aod::McMFTTrackLabels>;

  template <typename T>
  float getDCAz(const T& track)
  {
    if constexpr (requires { track.bestDCAZ(); }) {
      return track.bestDCAZ();
    } else {
      return 999.;
    }
  }

  std::unordered_set<int64_t> setRecCollSel;
  std::unordered_map<int64_t, float> mapVtxXrec;
  std::unordered_map<int64_t, float> mapVtxYrec;
  std::unordered_map<int64_t, float> mapVtxZrec;
  std::unordered_map<int64_t, float> mapVtxXgen;
  std::unordered_map<int64_t, float> mapVtxYgen;
  std::unordered_map<int64_t, float> mapVtxZgen;
  std::unordered_map<int64_t, int64_t> mapMcCollIdPerRecColl;
  std::unordered_map<int64_t, int64_t> mapMcToRec;

  template <typename C, typename MC>
  void buildLookupTable(C const& collisions, MC const& mcCollisions)
  {
    const auto& nMcColls = mcCollisions.size();
    const auto& nRecoColls = collisions.size();

    setRecCollSel.clear();
    setRecCollSel.reserve(nRecoColls);
    mapVtxXgen.clear();
    mapVtxXgen.reserve(nMcColls);
    mapVtxYgen.clear();
    mapVtxYgen.reserve(nMcColls);
    mapVtxZgen.clear();
    mapVtxZgen.reserve(nMcColls);
    mapVtxXrec.clear();
    mapVtxXrec.reserve(nRecoColls);
    mapVtxYrec.clear();
    mapVtxYrec.reserve(nRecoColls);
    mapVtxZrec.clear();
    mapVtxZrec.reserve(nRecoColls);
    mapMcCollIdPerRecColl.clear();
    mapMcCollIdPerRecColl.reserve(nRecoColls);
    mapMcToRec.clear();
    mapMcToRec.reserve(nRecoColls);

    int maxNcontributors = -1;
    for (auto const& collision : collisions) {
      // const int nContrib = collision.multPVTotalContributors();
      const int nContrib = collision.numContrib();
      if (maxNcontributors < nContrib) {
        maxNcontributors = nContrib;
        setRecCollSel.insert(collision.globalIndex());
      }
      mapVtxXrec.emplace(collision.globalIndex(), collision.posX());
      mapVtxYrec.emplace(collision.globalIndex(), collision.posY());
      mapVtxZrec.emplace(collision.globalIndex(), collision.posZ());
      mapMcCollIdPerRecColl.emplace(collision.globalIndex(), collision.mcCollisionId());
      mapMcToRec.emplace(collision.mcCollisionId(), collision.globalIndex());
    }
    for (const auto& mcCollision : mcCollisions) {
      const auto mcId = mcCollision.globalIndex();
      mapVtxXgen.emplace(mcId, mcCollision.posX());
      mapVtxYgen.emplace(mcId, mcCollision.posY());
      mapVtxZgen.emplace(mcId, mcCollision.posZ());
    }
  }

  /// \brief check good rec vertex of corresponding mc coll is available in compatible rec coll
  template <typename AT>
  bool isRecInCompColl(AT const& atrack)
  {
    auto const& ids = atrack.compatibleCollIds();
    bool isInCoColl = false;
    if (atrack.ambDegree() != 0) {
      const int mcCollId = atrack.mcParticle().mcCollisionId();
      if (!ids.empty()) {
        for (auto const& id : ids) {
          auto itMcCollId = mapMcCollIdPerRecColl.find(id);
          if (itMcCollId != mapMcCollIdPerRecColl.end()) {
            if (itMcCollId->second == mcCollId) {
              isInCoColl = true;
              return isInCoColl;
            }
          } else {
            return isInCoColl;
          }
        }
      }
    }
    return isInCoColl;
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

    if (cfgApplyZShiftFromCCDB) {
      auto* zShift = ccdb->getForTimeStamp<std::vector<float>>(cfgZShiftPath, bc.timestamp());
      if (zShift != nullptr && !zShift->empty()) {
        LOGF(info, "reading z shift %f from %s", (*zShift)[0], cfgZShiftPath.value);
        mZShift = (*zShift)[0];
      } else {
        LOGF(info, "z shift is not found in ccdb path %s. set to 0 cm", cfgZShiftPath.value);
        mZShift = 0;
      }
    } else {
      LOGF(info, "z shift is manually set to %f cm", cfgManualZShift.value);
      mZShift = cfgManualZShift;
    }
  }

  /// \brief RMS calculation
  /// \param vec  vector of values to compute RMS
  template <typename T>
  T computeRMS(std::vector<T>& vec)
  {
    T sum = std::accumulate(vec.begin(), vec.end(), 0.0);
    T mean = sum / vec.size();

    std::vector<T> diff(vec.size());
    std::transform(vec.begin(), vec.end(), diff.begin(), [mean](T x) { return x - mean; });
    T sqSum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    T stdev = std::sqrt(sqSum / vec.size());

    return stdev;
  }

  template <bool fillHis = true, typename B>
  bool isBestTrackSelected(const B& besttrack)
  {
    if constexpr (fillHis) {
      registry.fill(HIST("Tracks/hBestTrkSel"), static_cast<int>(TrkTrkBestSel::trkTrkBestSelAll));
    }
    if (besttrack.bestCollisionId() < CintZero) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Tracks/hBestTrkSel"), static_cast<int>(TrkTrkBestSel::trkTrkBestSelCollID));
    }
    if (besttrack.ambDegree() == CintZero) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Tracks/hBestTrkSel"), static_cast<int>(TrkTrkBestSel::trkTrkBestSelOrphan));
    }
    if (std::abs(besttrack.bestDCAXY()) >= trackCuts.maxDCAxy) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Tracks/hBestTrkSel"), static_cast<int>(TrkTrkBestSel::trkTrkBestSelDCAxyCut));
    }
    if (trackCuts.useDCAzCut) {
      const float bestDcaZ = getDCAz(besttrack);
      if (std::abs(bestDcaZ) >= trackCuts.maxDCAz) {
        return false;
      }
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Tracks/hBestTrkSel"), static_cast<int>(TrkTrkBestSel::trkTrkBestSelDCAzCut));
    }
    return true;
  }

  template <bool fillHis = true, typename T>
  bool isTrackSelected(const T& track)
  {
    if constexpr (fillHis) {
      registry.fill(HIST("Tracks/hTrkSel"), static_cast<int>(TrkSel::trkSelAll));
    }
    if (track.nClusters() < trackCuts.minNclusterMft) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Tracks/hTrkSel"), static_cast<int>(TrkSel::trkSelNCls));
    }
    if (trackCuts.useChi2Cut) {
      float nclMft = std::max(2.0f * track.nClusters() - 5.0f, 1.0f);
      float mftChi2NCl = track.chi2() / nclMft;
      if (mftChi2NCl > trackCuts.maxChi2NCl) {
        return false;
      }
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Tracks/hTrkSel"), static_cast<int>(TrkSel::trkSelChi2Ncl));
    }
    if (track.eta() < trackCuts.minEta || track.eta() > trackCuts.maxEta) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Tracks/hTrkSel"), static_cast<int>(TrkSel::trkSelEta));
    }
    if (trackCuts.usephiCut) {
      float phi = track.phi();
      o2::math_utils::bringTo02Pi(phi);
      if (phi < trackCuts.minPhi || trackCuts.maxPhi < phi) {
        return false;
      }
      if ((phi < trackCuts.phiCut) ||
          ((phi > PI - trackCuts.phiCut) && (phi < PI + trackCuts.phiCut)) ||
          (phi > TwoPI - trackCuts.phiCut) ||
          ((phi > ((PIHalf - 0.1) * PI) - trackCuts.phiCut) &&
           (phi < ((PIHalf - 0.1) * PI) + trackCuts.phiCut))) {
        return false;
      }
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Tracks/hTrkSel"), static_cast<int>(TrkSel::trkSelPhiCut));
    }
    if (trackCuts.usePtCut && track.pt() < trackCuts.minPt) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Tracks/hTrkSel"), static_cast<int>(TrkSel::trkSelPt));
    }
    if (trackCuts.requireCA && !track.isCA()) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Tracks/hTrkSel"), static_cast<int>(TrkSel::trkSelCA));
    }
    return true;
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

  template <bool fillHis = false, typename C>
  bool isGoodEvent(C const& collision)
  {
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), static_cast<int>(EvtSel::evtAll));
    }
    if (!collision.sel8()) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), static_cast<int>(EvtSel::evtSel));
    }
    if (eventCuts.requireIsGoodZvtxFT0VsPV && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), static_cast<int>(EvtSel::evtIsGoodZvtx));
    }
    if (eventCuts.requireRejectSameBunchPileup && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), static_cast<int>(EvtSel::evtNoSameBunchPileup));
    }
    if (collision.posZ() <= eventCuts.minZvtx || collision.posZ() >= eventCuts.maxZvtx) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), static_cast<int>(EvtSel::evtZvtxCut));
    }
    if (eventCuts.requireNoCollInTimeRangeStd &&
        !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), static_cast<int>(EvtSel::evtNoCollInTimeRangeStd));
    }
    if (eventCuts.requireNoCollInTimeRangeNarrow &&
        !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow)) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), static_cast<int>(EvtSel::evtNoCollInTimeRangeNarrow));
    }
    if (eventCuts.requireNoCollInTimeRangeStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict)) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), static_cast<int>(EvtSel::evtNoCollInTimeRangeStrict));
    }
    if (eventCuts.requireNoCollInRofStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStrict)) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), static_cast<int>(EvtSel::evtNoCollInRofStrict));
    }
    if (eventCuts.requireNoCollInRofStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), static_cast<int>(EvtSel::evtNoCollInRofStandard));
    }
    if (eventCuts.requireNoHighMultCollInPrevRof && !collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), static_cast<int>(EvtSel::evtNoHighMultCollInPrevRof));
    }
    if (eventCuts.requireGoodITSLayersAll && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), static_cast<int>(EvtSel::evtGoodITSLayersAll));
    }
    if (eventCuts.minOccupancy >= 0 &&
        getOccupancy(collision, eventCuts.occupancyEstimator) <
          eventCuts.minOccupancy) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), static_cast<int>(EvtSel::evtBelowMinOccup));
    }
    if (eventCuts.maxOccupancy >= 0 &&
        getOccupancy(collision, eventCuts.occupancyEstimator) >
          eventCuts.maxOccupancy) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), static_cast<int>(EvtSel::evtAboveMaxOccup));
    }
    if (rctCuts.requireRCTFlagChecker && !rctChecker(collision)) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), static_cast<int>(EvtSel::evtRCTFlagChecker));
    }
    if (rctCuts.requireRCTFlagCheckerExtra && !rctCheckerExtra(collision)) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), static_cast<int>(EvtSel::evtRCTFlagCheckerExtra));
    }
    return true;
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
      o2::math_utils::bringTo02Pi(phi);
      if (phi < trackCuts.minPhi || trackCuts.maxPhi < phi) {
        return false;
      }
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

  void processTimeAssocMC(MftTracksWCollsMC const& tracks,
                          aod::McParticles const& /*particles*/,
                          CollsGenCentFT0CExtra const& collisions,
                          aod::McCollisions const& mcCollisions)
  {
    buildLookupTable(collisions, mcCollisions);

    int nNoMC{0};
    for (const auto& collision : collisions) {
      registry.fill(HIST("TimeAssocMC/hTimeAssocMCEventStatus"), static_cast<int>(ReAssocMCEventStatus::kEvtReAsAll));
      if (!isGoodEvent<true>(collision)) {
        continue;
      }
      registry.fill(HIST("TimeAssocMC/hTimeAssocMCEventStatus"), static_cast<int>(ReAssocMCEventStatus::kEvtReAsSelected));
      if (!collision.has_mcCollision()) {
        continue;
      }
      int64_t recCollId = collision.globalIndex();
      auto itMC = mapMcCollIdPerRecColl.find(recCollId);
      if (itMC == mapMcCollIdPerRecColl.end()) {
        nNoMC++;
        LOGP(debug, "collison {} has no MC coll", recCollId);
        continue;
      }
      registry.fill(HIST("TimeAssocMC/hTimeAssocMCEventStatus"), static_cast<int>(ReAssocMCEventStatus::kEvtReAsHasMcColl));
      auto mcCollision = collision.mcCollision_as<aod::McCollisions>();
      if (gConf.cfgRemoveSplitVertex && (!setRecCollSel.contains(collision.globalIndex()))) {
        continue;
      }
      registry.fill(HIST("TimeAssocMC/hTimeAssocMCEventStatus"), static_cast<int>(ReAssocMCEventStatus::kEvtReAsSplitVtxRemoved));
      if (eventCuts.useZVtxCutMC && (std::abs(mcCollision.posZ()) >= eventCuts.maxZvtx)) {
        continue;
      }
      registry.fill(HIST("TimeAssocMC/hTimeAssocMCEventStatus"), static_cast<int>(ReAssocMCEventStatus::kEvtReAsZVtxCutMC));

      // Loop on collision compatible MFT tracks
      // Check: (1) good/bad vertices (2) if bad true vertex is available among rec vertices
      for (const auto& track : tracks) {
        if (!track.has_collision()) {
          continue;
        }
        auto trkCollId = track.has_collision() ? track.collisionId() : -1;
        auto ids = track.compatibleCollIds();
        // check if track is associated to rec coll
        // if (trkCollId != recCollId) {
        //   continue;
        // }
        if (!setRecCollSel.contains(trkCollId)) {
          continue;
        }
        registry.fill(HIST("TimeAssocMC/hAmbTrackType"), static_cast<int>(AmbTrkType::kAll));
        if (ids.empty()) {
          registry.fill(HIST("TimeAssocMC/hAmbTrackType"), static_cast<int>(AmbTrkType::kOrphan));
        }
        if (gConf.cfgRemoveOrphanTracks && ids.empty()) {
          continue;
        }
        if (gConf.cfgRemoveTrivialAssoc) {
          if (ids.empty() || (ids.size() == 1 && trkCollId == ids[0])) {
            continue;
          }
        }
        if (gConf.cfgRemoveAmbiguousTracks && (track.compatibleCollIds().size() != 1)) {
          continue;
        }

        if (!ids.empty()) {
          if (ids.size() == 1) {
            if (trkCollId == ids[0]) {
              registry.fill(HIST("TimeAssocMC/hAmbTrackType"), static_cast<int>(AmbTrkType::kNonAmb));
            } else if (trkCollId != ids[0]) {
              registry.fill(HIST("TimeAssocMC/hAmbTrackType"), static_cast<int>(AmbTrkType::kAmb));
            } else {
              registry.fill(HIST("TimeAssocMC/hAmbTrackType"), static_cast<int>(AmbTrkType::kNonAmbSame));
            }
          } else {
            registry.fill(HIST("TimeAssocMC/hAmbTrackType"), static_cast<int>(AmbTrkType::kAmbGt1));
          }
        } else {
          registry.fill(HIST("TimeAssocMC/hAmbTrackType"), static_cast<int>(AmbTrkType::kOrphanNull));
        }
        if (gConf.cfgUseTrackSel && !isTrackSelected<true>(track)) {
          continue;
        }

        bool isTrueVtx = false;
        int vtxFlag = static_cast<int>(VertexStatusMC::kNull);

        float vtxX = -1.;
        float vtxY = -1.;
        float vtxZ = -1.;
        float deltaXv1 = -1.;
        float deltaYv1 = -1.;
        float deltaZv1 = -1.;
        float deltaXv2 = -1.;
        float deltaYv2 = -1.;
        float deltaZv2 = -1.;

        if (track.collisionId() >= 0 && track.has_mcParticle() && track.mcMask() == trackCuts.selMcMask) {
          const auto& mcPart = track.mcParticle();
          if (!isChrgParticle(mcPart.pdgCode())) {
            continue;
          }
          if (gConf.cfgUseParticleSel && !isParticleSelected(mcPart)) {
            continue;
          }
          int64_t mcPartId = mcPart.mcCollisionId();
          // check if rec vertex is available in MC collisions
          for (const auto& mcTrkId : mapMcToRec) {
            if (mcTrkId.second == mcPartId) {
              isTrueVtx = true;
              break;
            }
          }
          // check if there is good or bad collision
          auto itMCTrk = mapMcCollIdPerRecColl.find(trkCollId);
          if (itMCTrk != mapMcCollIdPerRecColl.end()) {
            int mcTrkCollId = itMCTrk->second;
            if (mcPartId == mcTrkCollId) { // particle.mcCollisionId == collision.mcCollisionId -> good vtx
              vtxFlag = static_cast<int>(VertexStatusMC::kGood);
            } else { // wrong vtx
              vtxFlag = static_cast<int>(VertexStatusMC::kBad);
            }
          }
          if (!mapVtxXrec.contains(trkCollId)) {
            continue;
          }
          if (!mapVtxYrec.contains(trkCollId)) {
            continue;
          }
          if (!mapVtxZrec.contains(trkCollId)) {
            continue;
          }
          if (!mapMcCollIdPerRecColl.contains(trkCollId)) {
            continue;
          }
          vtxX = mapVtxXrec.find(trkCollId)->second;
          vtxY = mapVtxYrec.find(trkCollId)->second;
          vtxZ = mapVtxZrec.find(trkCollId)->second;
          // LOGP(info, "\t ---> \t .... \t vtxZrec: {} - collision.posZ(): {}", vtxZrec, collision.posZ());
          int64_t mcCollIdRec = mapMcCollIdPerRecColl.find(trkCollId)->second;
          // int64_t mcCollId = itMC->second;
          // LOGP(info, "\t ---> \t .... \t mcCollIdRec: {} - mcCollId: {} - bestMCCol: {}", mcCollIdRec, mcCollId, bestMCCol);
          if (!mapVtxXgen.contains(mcCollIdRec)) {
            continue;
          }
          if (!mapVtxYgen.contains(mcCollIdRec)) {
            continue;
          }
          if (!mapVtxZgen.contains(mcCollIdRec)) {
            continue;
          }
          // vertex resolution - ver 1
          // rec coll vtx - mc associated to orig rec coll (first in time)
          deltaXv1 = vtxX - mcCollision.posX();
          deltaYv1 = vtxY - mcCollision.posY();
          deltaZv1 = vtxZ - mcCollision.posZ();
          // vertex resolution - ver 2
          // rec coll vtx - mc associated to orig rec coll (first in time)
          deltaXv2 = vtxX - mapVtxXgen.find(mcCollIdRec)->second;
          deltaYv2 = vtxY - mapVtxYgen.find(mcCollIdRec)->second;
          deltaZv2 = vtxZ - mapVtxZgen.find(mcCollIdRec)->second;

          registry.fill(HIST("TimeAssocMC/VtxStatus"), vtxFlag);
          registry.fill(HIST("TimeAssocMC/hVertexResV1"), deltaXv1, deltaYv1, deltaZv1);
          registry.fill(HIST("TimeAssocMC/hVertexResV2"), deltaXv2, deltaYv2, deltaZv2);

          registry.fill(HIST("TimeAssocMC/hAmbTrkTypeAssocFlag"), static_cast<int>(AmbTrkTypeAssocFlag::kSelAll));
          if (isTrueVtx) {
            registry.fill(HIST("TimeAssocMC/hTimeAssocCheckVtxType"), static_cast<int>(AssocCheckVtxType::kAllVtxTrue));
          } else {
            registry.fill(HIST("TimeAssocMC/hTimeAssocCheckVtxType"), static_cast<int>(AssocCheckVtxType::kAllVtxFalse));
          }
          if (vtxFlag == static_cast<int>(VertexStatusMC::kGood)) {
            registry.fill(HIST("TimeAssocMC/hAmbTrkTypeAssocFlag"), static_cast<int>(AmbTrkTypeAssocFlag::kSelGoodVtx));
            registry.fill(HIST("TimeAssocMC/hTimeAssocCheckVtxType"), static_cast<int>(AssocCheckVtxType::kAllGoodVtx));
            if (isTrueVtx) {
              registry.fill(HIST("TimeAssocMC/hVTXkSelGoodVtxTrue"), deltaXv2, deltaYv2, deltaZv2);
              registry.fill(HIST("TimeAssocMC/hAmbTrkTypeAssocFlag"), static_cast<int>(AmbTrkTypeAssocFlag::kSelGoodVtxTrue));
              registry.fill(HIST("TimeAssocMC/hTimeAssocCheckVtxType"), static_cast<int>(AssocCheckVtxType::kAllGoodVtxTrue));
            } else {
              registry.fill(HIST("TimeAssocMC/hTimeAssocCheckVtxType"), static_cast<int>(AssocCheckVtxType::kAllGoodVtxFalse));
            }
          } else if (vtxFlag == static_cast<int>(VertexStatusMC::kBad)) {
            registry.fill(HIST("TimeAssocMC/hAmbTrkTypeAssocFlag"), static_cast<int>(AmbTrkTypeAssocFlag::kSelBadVtx));
            registry.fill(HIST("TimeAssocMC/hTimeAssocCheckVtxType"), static_cast<int>(AssocCheckVtxType::kAllBadVtx));
            if (isTrueVtx) {
              registry.fill(HIST("TimeAssocMC/hVTXkSelBadVtxTrue"), deltaXv2, deltaYv2, deltaZv2);
              registry.fill(HIST("TimeAssocMC/hAmbTrkTypeAssocFlag"), static_cast<int>(AmbTrkTypeAssocFlag::kSelBadVtxTrue));
              registry.fill(HIST("TimeAssocMC/hTimeAssocCheckVtxType"), static_cast<int>(AssocCheckVtxType::kAllBadVtxTrue));
            } else {
              registry.fill(HIST("TimeAssocMC/hTimeAssocCheckVtxType"), static_cast<int>(AssocCheckVtxType::kAllBadVtxFalse));
            }
          }

          if (!ids.empty()) {
            if (ids.size() == 1) {
              registry.fill(HIST("TimeAssocMC/hVTXkSelNonAmb"), deltaXv2, deltaYv2, deltaZv2);
              registry.fill(HIST("TimeAssocMC/hAmbTrkTypeAssocFlag"), static_cast<int>(AmbTrkTypeAssocFlag::kSelNonAmb));
              if (vtxFlag == static_cast<int>(VertexStatusMC::kGood)) {
                registry.fill(HIST("TimeAssocMC/hVTXkSelNonAmbGoodVtx"), deltaXv2, deltaYv2, deltaZv2);
                registry.fill(HIST("TimeAssocMC/hAmbTrkTypeAssocFlag"), static_cast<int>(AmbTrkTypeAssocFlag::kSelNonAmbGoodVtx));
                if (isTrueVtx) {
                  registry.fill(HIST("TimeAssocMC/hVTXkSelNonAmbGoodVtxTrue"), deltaXv2, deltaYv2, deltaZv2);
                  registry.fill(HIST("TimeAssocMC/hAmbTrkTypeAssocFlag"), static_cast<int>(AmbTrkTypeAssocFlag::kSelNonAmbGoodVtxTrue));
                }
              } else if (vtxFlag == static_cast<int>(VertexStatusMC::kBad)) {
                registry.fill(HIST("TimeAssocMC/hAmbTrkTypeAssocFlag"), static_cast<int>(AmbTrkTypeAssocFlag::kSelNonAmbBadVtx));
                if (isTrueVtx) {
                  registry.fill(HIST("TimeAssocMC/hVTXkSelNonAmbBadVtxTrue"), deltaXv2, deltaYv2, deltaZv2);
                  registry.fill(HIST("TimeAssocMC/hAmbTrkTypeAssocFlag"), static_cast<int>(AmbTrkTypeAssocFlag::kSelNonAmbBadVtxTrue));
                }
              }
              if (trkCollId == ids[0]) { // non ambiguous
                registry.fill(HIST("TimeAssocMC/hVTXkSelNonAmbID"), deltaXv2, deltaYv2, deltaZv2);
                registry.fill(HIST("TimeAssocMC/hAmbTrkTypeAssocFlag"), static_cast<int>(AmbTrkTypeAssocFlag::kSelNonAmbID));
                if (vtxFlag == static_cast<int>(VertexStatusMC::kGood)) {
                  if (isTrueVtx) {
                    registry.fill(HIST("TimeAssocMC/hVTXkSelNonAmbIDGoodVtxTrue"), deltaXv2, deltaYv2, deltaZv2);
                    registry.fill(HIST("TimeAssocMC/hAmbTrkTypeAssocFlag"), static_cast<int>(AmbTrkTypeAssocFlag::kSelNonAmbIDGoodVtxTrue));
                  }
                } else if (vtxFlag == static_cast<int>(VertexStatusMC::kBad)) {
                  if (isTrueVtx) {
                    registry.fill(HIST("TimeAssocMC/hVTXkSelNonAmbIDBadVtxTrue"), deltaXv2, deltaYv2, deltaZv2);
                    registry.fill(HIST("TimeAssocMC/hAmbTrkTypeAssocFlag"), static_cast<int>(AmbTrkTypeAssocFlag::kSelNonAmbIDBadVtxTrue));
                  }
                }
              } else if (trkCollId != ids[0]) {
                registry.fill(HIST("TimeAssocMC/hVTXkSelAmbID"), deltaXv2, deltaYv2, deltaZv2);
                registry.fill(HIST("TimeAssocMC/hAmbTrkTypeAssocFlag"), static_cast<int>(AmbTrkTypeAssocFlag::kSelAmbID));
                if (vtxFlag == static_cast<int>(VertexStatusMC::kGood)) {
                  if (isTrueVtx) {
                    registry.fill(HIST("TimeAssocMC/hVTXkSelAmbIDGoodVtxTrue"), deltaXv2, deltaYv2, deltaZv2);
                    registry.fill(HIST("TimeAssocMC/hAmbTrkTypeAssocFlag"), static_cast<int>(AmbTrkTypeAssocFlag::kSelAmbIDGoodVtxTrue));
                  }
                } else if (vtxFlag == static_cast<int>(VertexStatusMC::kBad)) {
                  if (isTrueVtx) {
                    registry.fill(HIST("TimeAssocMC/hVTXkSelAmbIDBadVtxTrue"), deltaXv2, deltaYv2, deltaZv2);
                    registry.fill(HIST("TimeAssocMC/hAmbTrkTypeAssocFlag"), static_cast<int>(AmbTrkTypeAssocFlag::kSelAmbIDBadVtxTrue));
                  }
                }
              } else { // non ambiguous (extra)
                registry.fill(HIST("TimeAssocMC/hVTXkSelNonAmbIDExtra"), deltaXv2, deltaYv2, deltaZv2);
                registry.fill(HIST("TimeAssocMC/hAmbTrkTypeAssocFlag"), static_cast<int>(AmbTrkTypeAssocFlag::kSelNonAmbIDExtra));
                if (vtxFlag == static_cast<int>(VertexStatusMC::kGood)) {
                  if (isTrueVtx) {
                    registry.fill(HIST("TimeAssocMC/hVTXkSelNonAmbIDExtraGoodVtxTrue"), deltaXv2, deltaYv2, deltaZv2);
                    registry.fill(HIST("TimeAssocMC/hAmbTrkTypeAssocFlag"), static_cast<int>(AmbTrkTypeAssocFlag::kSelNonAmbIDExtraGoodVtxTrue));
                  }
                } else if (vtxFlag == static_cast<int>(VertexStatusMC::kBad)) {
                  if (isTrueVtx) {
                    registry.fill(HIST("TimeAssocMC/hVTXkSelNonAmbIDExtraBadVtxTrue"), deltaXv2, deltaYv2, deltaZv2);
                    registry.fill(HIST("TimeAssocMC/hAmbTrkTypeAssocFlag"), static_cast<int>(AmbTrkTypeAssocFlag::kSelNonAmbIDExtraBadVtxTrue));
                  }
                }
              }
            } else { // ambiguous
              registry.fill(HIST("TimeAssocMC/hVTXkSelAmb"), deltaXv2, deltaYv2, deltaZv2);
              registry.fill(HIST("TimeAssocMC/hAmbTrkTypeAssocFlag"), static_cast<int>(AmbTrkTypeAssocFlag::kSelAmb));
              if (vtxFlag == static_cast<int>(VertexStatusMC::kGood)) {
                if (isTrueVtx) {
                  registry.fill(HIST("TimeAssocMC/hVTXkSelAmbGoodVtxTrue"), deltaXv2, deltaYv2, deltaZv2);
                  registry.fill(HIST("TimeAssocMC/hAmbTrkTypeAssocFlag"), static_cast<int>(AmbTrkTypeAssocFlag::kSelAmbGoodVtxTrue));
                }
              } else if (vtxFlag == static_cast<int>(VertexStatusMC::kBad)) {
                if (isTrueVtx) {
                  registry.fill(HIST("TimeAssocMC/hVTXkSelAmbBadVtxTrue"), deltaXv2, deltaYv2, deltaZv2);
                  registry.fill(HIST("TimeAssocMC/hAmbTrkTypeAssocFlag"), static_cast<int>(AmbTrkTypeAssocFlag::kSelAmbBadVtxTrue));
                }
              }
            }
          } else {
            registry.fill(HIST("TimeAssocMC/hAmbTrkTypeAssocFlag"), static_cast<int>(AmbTrkTypeAssocFlag::kSelOrphanNull));
          }
        }
      }
    }
    LOG(info) << "No MC: " << nNoMC;
  }
  PROCESS_SWITCH(QaReassocMFTPbPb, processTimeAssocMC, "process MC: check MFT tracks in compatible collisions (time-associasted)", false);

  void processTimeAssocWithReassocMC(MftTracksWCollsMC const& /*tracks*/,
                                     BestTracks3dWCollsMC const& besttracks,
                                     aod::McParticles const& /*particles*/,
                                     ExtBCs const& bcs,
                                     CollsGenCentFT0CExtra const& collisions,
                                     aod::McCollisions const& mcCollisions)
  {
    if (bcs.size() == 0) {
      return;
    }
    auto bc = bcs.begin();
    initCCDB(bc);

    buildLookupTable(collisions, mcCollisions);

    int nNoMC{0};
    for (const auto& collision : collisions) {
      int64_t recCollId = collision.globalIndex();
      auto itMC = mapMcCollIdPerRecColl.find(recCollId);
      if (itMC == mapMcCollIdPerRecColl.end()) {
        nNoMC++;
        LOGP(debug, "collison {} has no MC coll", recCollId);
        continue;
      }
      if (!collision.has_mcCollision()) {
        continue;
      }
      auto mcCollision = collision.mcCollision_as<aod::McCollisions>();
      if (gConf.cfgRemoveSplitVertex && (!setRecCollSel.contains(collision.globalIndex()))) {
        continue;
      }
      if (!isGoodEvent<true>(collision)) {
        continue;
      }
      if (eventCuts.useZVtxCutMC && (std::abs(mcCollision.posZ()) >= eventCuts.maxZvtx)) {
        continue;
      }

      for (auto const& atrack : besttracks) {
        if (!isBestTrackSelected<true>(atrack)) {
          continue;
        }
        auto itrack = atrack.mfttrack_as<MftTracksWCollsMC>();
        if (!isTrackSelected(itrack)) {
          continue;
        }
        if (gConf.cfgUseTrackSel && !isTrackSelected<true>(itrack)) {
          continue;
        }
        auto trkBestCollId = itrack.has_collision() ? atrack.bestCollisionId() : -1;
        // if (trkBestCollId != recCollId) { // check if best track is associated to original rec coll
        //   continue;
        // }

        bool isTrueVtx = false;
        bool isRecGoodMatch = false;
        bool isRecGood = false;
        int vtxFlag = static_cast<int>(VertexStatusMC::kNull);

        float vtxXbest = -1.;
        float vtxYbest = -1.;
        float vtxZbest = -1.;
        float deltaXOrigMc = -1.;
        float deltaYOrigMc = -1.;
        float deltaZOrigMc = -1.;
        // dca best rec - dca mc associated to best rec coll
        float deltaX = -1.;
        float deltaY = -1.;
        float deltaZ = -1.;

        if (trkBestCollId >= 0 && atrack.has_mcParticle() && atrack.mcMask() == trackCuts.selMcMask) {
          auto itRecToMc = mapMcCollIdPerRecColl.find(trkBestCollId); // try mfttrackId ???
          int64_t mcPartId = atrack.mcParticle().mcCollisionId();
          auto const& idCompColl = atrack.compatibleCollIds();

          // check if rec vertex is available in MC collisions
          for (const auto& mcTrkId : mapMcCollIdPerRecColl) {
            if (mcTrkId.second == mcPartId) {
              isTrueVtx = true;
              break;
            }
          }
          // check good rec vertex of corresponding mc coll is available in compatible rec coll
          if (!idCompColl.empty()) {
            for (auto const& id : idCompColl) {
              auto itMcCollId = mapMcCollIdPerRecColl.find(id);
              if (itMcCollId != mapMcCollIdPerRecColl.end()) {
                if (itMcCollId->second == mcPartId) {
                  isRecGoodMatch = true;
                  break;
                }
              }
            }
          }
          // check if there is good or bad collision
          if (itRecToMc != mapMcCollIdPerRecColl.end()) {
            int mcTrkCollId = itRecToMc->second;
            if (mcPartId == mcTrkCollId) { // particle.mcCollisionId == collision.mcCollisionId -> good vtx
              vtxFlag = static_cast<int>(VertexStatusMC::kGood);
            } else { // wrong vtx
              vtxFlag = static_cast<int>(VertexStatusMC::kBad);
            }
          }

          //
          // check: vertex resolution of time-to-coll reassoc: pos(rec coll best) - gen (mc coll)
          //
          // const auto& particle = atrack.mcParticle_as<aod::McParticles>();
          // // // auto collision = atrack.collision_as<CollisionsWithMCLabels>();  // not in use
          // const auto& mcColl = particle.mcCollision_as<McCollsEx>();
          if (!mapVtxXrec.contains(trkBestCollId)) {
            continue;
          }
          if (!mapVtxYrec.contains(trkBestCollId)) {
            continue;
          }
          if (!mapVtxZrec.contains(trkBestCollId)) {
            continue;
          }
          if (!mapMcCollIdPerRecColl.contains(trkBestCollId)) {
            continue;
          }
          vtxXbest = mapVtxXrec.find(trkBestCollId)->second;
          vtxYbest = mapVtxYrec.find(trkBestCollId)->second;
          vtxZbest = mapVtxZrec.find(trkBestCollId)->second;
          // LOGP(info, "\t ---> \t .... \t vtxZrec: {} - collision.posZ(): {}", vtxZrec, collision.posZ());
          int64_t mcCollIdRec = mapMcCollIdPerRecColl.find(trkBestCollId)->second;
          // int64_t mcCollId = itMC->second;
          // LOGP(info, "\t ---> \t .... \t mcCollIdRec: {} - mcCollId: {} - bestMCCol: {}", mcCollIdRec, mcCollId, bestMCCol);
          if (!mapVtxXgen.contains(mcCollIdRec)) {
            continue;
          }
          if (!mapVtxYgen.contains(mcCollIdRec)) {
            continue;
          }
          if (!mapVtxZgen.contains(mcCollIdRec)) {
            continue;
          }
          // vertex resolution: best rec - mc associated to orig rec coll (first in time)
          deltaXOrigMc = vtxXbest - mcCollision.posX();
          deltaYOrigMc = vtxYbest - mcCollision.posY();
          deltaZOrigMc = vtxZbest - mcCollision.posZ();
          // vertex resolution: best rec - mc associated to best rec coll
          deltaX = vtxXbest - mapVtxXgen.find(mcCollIdRec)->second;
          deltaY = vtxYbest - mapVtxYgen.find(mcCollIdRec)->second;
          deltaZ = vtxZbest - mapVtxZgen.find(mcCollIdRec)->second;
        } // has_mcParticle

        registry.fill(HIST("TimeAssocWithReassocMC/VtxStatus"), vtxFlag);

        hTimeAssocWithReassocMC[0]->Fill(deltaXOrigMc, deltaYOrigMc, deltaZOrigMc);
        hTimeAssocWithReassocMC[1]->Fill(deltaX, deltaY, deltaZ);

        if (isTrueVtx) {
          hTimeAssocWithReassocMC[3]->Fill(deltaXOrigMc, deltaYOrigMc, deltaZOrigMc);
          hTimeAssocWithReassocMC[12]->Fill(deltaX, deltaY, deltaZ);
          registry.fill(HIST("TimeAssocWithReassocMC/hReassocCheckVtxType"), static_cast<int>(ReassocCheckVtxType::kIsTrueVtxAllTrue));
        } else {
          registry.fill(HIST("TimeAssocWithReassocMC/hReassocCheckVtxType"), static_cast<int>(ReassocCheckVtxType::kIsTrueVtxAllFalse));
        }
        if (isRecGood) {
          hTimeAssocWithReassocMC[4]->Fill(deltaXOrigMc, deltaYOrigMc, deltaZOrigMc);
          hTimeAssocWithReassocMC[13]->Fill(deltaX, deltaY, deltaZ);
          registry.fill(HIST("TimeAssocWithReassocMC/hReassocCheckVtxType"), static_cast<int>(ReassocCheckVtxType::kIsRecGoodAllTrue));
        } else {
          registry.fill(HIST("TimeAssocWithReassocMC/hReassocCheckVtxType"), static_cast<int>(ReassocCheckVtxType::kIsRecGoodAllFalse));
        }
        if (isRecGoodMatch) {
          hTimeAssocWithReassocMC[5]->Fill(deltaXOrigMc, deltaYOrigMc, deltaZOrigMc);
          hTimeAssocWithReassocMC[14]->Fill(deltaX, deltaY, deltaZ);
          registry.fill(HIST("TimeAssocWithReassocMC/hReassocCheckVtxType"), static_cast<int>(ReassocCheckVtxType::kIsRecGoodMatchAllTrue));
        } else {
          registry.fill(HIST("TimeAssocWithReassocMC/hReassocCheckVtxType"), static_cast<int>(ReassocCheckVtxType::kIsRecGoodMatchAllFalse));
        }

        if (vtxFlag == static_cast<int>(VertexStatusMC::kGood)) {
          if (isTrueVtx) {
            hTimeAssocWithReassocMC[6]->Fill(deltaXOrigMc, deltaYOrigMc, deltaZOrigMc);
            hTimeAssocWithReassocMC[15]->Fill(deltaX, deltaY, deltaZ);
            registry.fill(HIST("TimeAssocWithReassocMC/hReassocCheckVtxType"), static_cast<int>(ReassocCheckVtxType::kIsTrueVtxVsGoodVtxTrue));
          } else {
            registry.fill(HIST("TimeAssocWithReassocMC/hReassocCheckVtxType"), static_cast<int>(ReassocCheckVtxType::kIsTrueVtxVsGoodVtxFalse));
          }
          if (isRecGood) {
            hTimeAssocWithReassocMC[7]->Fill(deltaXOrigMc, deltaYOrigMc, deltaZOrigMc);
            hTimeAssocWithReassocMC[16]->Fill(deltaX, deltaY, deltaZ);
            registry.fill(HIST("TimeAssocWithReassocMC/hReassocCheckVtxType"), static_cast<int>(ReassocCheckVtxType::kIsRecGoodVsGoodVtxTrue));
          } else {
            registry.fill(HIST("TimeAssocWithReassocMC/hReassocCheckVtxType"), static_cast<int>(ReassocCheckVtxType::kIsRecGoodVsGoodVtxFalse));
          }
          if (isRecGoodMatch) {
            hTimeAssocWithReassocMC[8]->Fill(deltaXOrigMc, deltaYOrigMc, deltaZOrigMc);
            hTimeAssocWithReassocMC[17]->Fill(deltaX, deltaY, deltaZ);
            registry.fill(HIST("TimeAssocWithReassocMC/hReassocCheckVtxType"), static_cast<int>(ReassocCheckVtxType::kIsRecGoodMatchVsGoodVtxTrue));
          } else {
            registry.fill(HIST("TimeAssocWithReassocMC/hReassocCheckVtxType"), static_cast<int>(ReassocCheckVtxType::kIsRecGoodMatchVsGoodVtxFalse));
          }
        } else if (vtxFlag == static_cast<int>(VertexStatusMC::kBad)) {
          if (isTrueVtx) {
            hTimeAssocWithReassocMC[9]->Fill(deltaXOrigMc, deltaYOrigMc, deltaZOrigMc);
            hTimeAssocWithReassocMC[18]->Fill(deltaX, deltaY, deltaZ);
            registry.fill(HIST("TimeAssocWithReassocMC/hReassocCheckVtxType"), static_cast<int>(ReassocCheckVtxType::kIsTrueVtxVsBadVtxTrue));
          } else {
            registry.fill(HIST("TimeAssocWithReassocMC/hReassocCheckVtxType"), static_cast<int>(ReassocCheckVtxType::kIsTrueVtxVsBadVtxFalse));
          }
          if (isRecGood) {
            hTimeAssocWithReassocMC[10]->Fill(deltaXOrigMc, deltaYOrigMc, deltaZOrigMc);
            hTimeAssocWithReassocMC[19]->Fill(deltaX, deltaY, deltaZ);
            registry.fill(HIST("TimeAssocWithReassocMC/hReassocCheckVtxType"), static_cast<int>(ReassocCheckVtxType::kIsRecGoodVsBadVtxTrue));
          } else {
            registry.fill(HIST("TimeAssocWithReassocMC/hReassocCheckVtxType"), static_cast<int>(ReassocCheckVtxType::kIsRecGoodVsBadVtxFalse));
          }
          if (isRecGoodMatch) {
            hTimeAssocWithReassocMC[11]->Fill(deltaXOrigMc, deltaYOrigMc, deltaZOrigMc);
            hTimeAssocWithReassocMC[20]->Fill(deltaX, deltaY, deltaZ);
            registry.fill(HIST("TimeAssocWithReassocMC/hReassocCheckVtxType"), static_cast<int>(ReassocCheckVtxType::kIsRecGoodMatchVsBadVtxTrue));
          } else {
            registry.fill(HIST("TimeAssocWithReassocMC/hReassocCheckVtxType"), static_cast<int>(ReassocCheckVtxType::kIsRecGoodMatchVsBadVtxFalse));
          }
        }
      }
    }
    LOG(info) << "No MC: " << nNoMC;
  }
  PROCESS_SWITCH(QaReassocMFTPbPb, processTimeAssocWithReassocMC, "process MC: check MFT tracks in reassociation with compatible collisions", false);

  Partition<MftTracksWCollsMC> tracksInAcc = (aod::fwdtrack::eta < trackCuts.maxEta) && (aod::fwdtrack::eta > trackCuts.minEta);

  void processAssocMC(MftTracksWCollsMC const& tracks,
                      aod::McParticles const& /*particles*/,
                      CollsGenCentFT0CExtra const& collisions,
                      aod::McCollisions const& mcCollisions)
  {
    buildLookupTable(collisions, mcCollisions);

    for (const auto& collision : collisions) {
      if (!isGoodEvent<true>(collision)) {
        continue;
      }
      // Select collisions with the largest number of contributors
      if (gConf.cfgRemoveSplitVertex && (!setRecCollSel.contains(collision.globalIndex()))) {
        continue;
      }
      if (!collision.has_mcCollision()) {
        continue;
      }
      auto mcCollision = collision.mcCollision_as<aod::McCollisions>();
      if (eventCuts.useZVtxCutMC && (std::abs(mcCollision.posZ()) >= eventCuts.maxZvtx)) {
        continue;
      }

      int nTrk = 0, nFakeTrk = 0, nGoodTrk = 0;
      for (const auto& track : tracks) {
        uint index = uint(track.collisionId() >= 0);
        if (track.has_mcParticle() && track.mcMask() == trackCuts.selMcMask) {
          // auto particle = track.mcParticle_as<aod::McParticles>();
          const auto& particle = track.mcParticle();
          auto trkCollId = track.has_collision() ? track.collisionId() : -1;
          auto ids = track.compatibleCollIds();
          registry.fill(HIST("TrackToColl/histTrackNumColls"), ids.size());
          registry.fill(HIST("TrackToColl/hAmbTrackType"), static_cast<int>(AmbTrkType::kAll));

          if (ids.empty() || (ids.size() == 1 && trkCollId == ids[0])) {
            if (ids.empty()) {
              registry.fill(HIST("TrackToColl/hAmbTrackType"), static_cast<int>(AmbTrkType::kOrphan));
            }
            if (ids.size() == 1 && trkCollId == ids[0]) {
              registry.fill(HIST("TrackToColl/hAmbTrackType"), static_cast<int>(AmbTrkType::kNonAmb));
            }
          }
          if (gConf.cfgRemoveOrphanTracks && ids.empty()) {
            continue;
          }
          if (gConf.cfgRemoveTrivialAssoc) {
            if (ids.empty() || (ids.size() == 1 && trkCollId == ids[0])) {
              registry.fill(HIST("TrackToColl/histNonAmbTrackNumColls"), ids.size());
              continue;
            }
          }
          if (gConf.cfgRemoveAmbiguousTracks && (track.compatibleCollIds().size() != 1)) {
            continue;
          }
          nTrk++;
          if ((particle.mcCollisionId() != collision.mcCollision().globalIndex())) {
            nFakeTrk++;
            continue;
          }
          if (collision.mcCollisionId() == particle.mcCollisionId()) {
            nGoodTrk++;
          }
          if (!ids.empty()) {
            if (ids.size() == 1) {
              if (trkCollId != ids[0]) {
                registry.fill(HIST("TrackToColl/hAmbTrackType"), static_cast<int>(AmbTrkType::kAmb));
              } else {
                registry.fill(HIST("TrackToColl/hAmbTrackType"), static_cast<int>(AmbTrkType::kNonAmbSame));
              }
            } else {
              registry.fill(HIST("TrackToColl/hAmbTrackType"), static_cast<int>(AmbTrkType::kAmbGt1));
            }
          }

          bool isAmbiguous = (ids.size() > 1);
          if (isAmbiguous) {
            registry.fill(HIST("TrackToColl/histAmbTrackNumColls"), ids.size());
            std::vector<float> ambVtxZ{};
            for (const auto& collIdx : ids) {
              const auto& ambColl = collisions.rawIteratorAt(collIdx);
              ambVtxZ.push_back(ambColl.posZ());
            }
            if (!ambVtxZ.empty()) {
              registry.fill(HIST("TrackToColl/histAmbTrackZvtxRMS"), computeRMS(ambVtxZ));
            }
          }

          float deltaX = -999.f;
          float deltaY = -999.f;
          float deltaZ = -999.f;
          if (index != 0u) {
            const auto& collision = track.collision_as<CollsGenCentFT0CExtra>();
            deltaX = collision.posX() - mcCollision.posX();
            deltaY = collision.posY() - mcCollision.posY();
            deltaZ = collision.posZ() - mcCollision.posZ();
            if (collision.has_mcCollision() && collision.mcCollisionId() == particle.mcCollisionId()) {
              hCollAssoc[index + 1]->Fill(track.pt(), track.eta(), deltaX, deltaY, deltaZ);
            } else {
              if (isAmbiguous) {
                for (const auto& collIdx : ids) {
                  auto ambColl = collisions.rawIteratorAt(collIdx);
                  if (ambColl.has_mcCollision() && ambColl.mcCollisionId() == particle.mcCollisionId()) {
                    hCollAssoc[index + 2]->Fill(track.pt(), track.eta(), deltaX, deltaY, deltaZ);
                    // hCollAssoc[index + 2]->Fill(track.pt(), track.eta(), ambColl.posX() - mcCollision.posX(), ambColl.posY() - mcCollision.posY(), ambColl.posZ() - mcCollision.posZ());
                    break;
                  }
                }
              }
            }
            hCollAssoc[index]->Fill(track.pt(), track.eta(), deltaX, deltaY, deltaZ);
          }
        } else {
          hCollAssoc[index]->Fill(track.pt(), track.eta(), -999.f, -999.f, -999.f);
        }
      }
      float frac = (nTrk > 0) ? static_cast<float>(nGoodTrk) / nTrk : -1.;
      registry.fill(HIST("TrackToColl/histFracGoodTracks"), frac);
      float fracFake = (nTrk > 0) ? static_cast<float>(nFakeTrk) / nTrk : -1.;
      registry.fill(HIST("TrackToColl/histFracTracksFakeMcColl"), fracFake);
    }
  }

  PROCESS_SWITCH(QaReassocMFTPbPb, processAssocMC, "Process collision-association information, requires extra table from TrackToCollisionAssociation task (fillTableOfCollIdsPerTrack=true)", false);

  void processReAssoc3dMC(BestTracks3dWCollsMC const& besttracks,
                          MftTracksWCollsMC const& /*tracks*/,
                          aod::McParticles const& /*particles*/,
                          CollsGenCentFT0CExtra const& collisions,
                          aod::McCollisions const& mcCollisions)
  {
    float cGen = -1;
    float crecMin = 105.f;
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
    float occGen = -1.;
    for (const auto& collision : collisions) {
      if (isGoodEvent<false>(collision)) {
        float o = getOccupancy(collision, eventCuts.occupancyEstimator);
        if (o > occGen) {
          occGen = o;
        }
      }
    }

    registry.fill(HIST("ReAssocMC/EvtGenRecReassoc"), 1., cGen, occGen);

    buildLookupTable(collisions, mcCollisions);

    int nNoMC{0};
    for (const auto& collision : collisions) {
      auto occ = getOccupancy(collision, eventCuts.occupancyEstimator);
      float crec = getRecoCent(collision);
      registry.fill(HIST("ReAssocMC/EvtGenRecReassoc"), 2., crec, occ);
      registry.fill(HIST("ReAssocMC/hReAssocMCEventStatus"), static_cast<int>(ReAssocMCEventStatus::kEvtReAsAll), crec, occ);
      if (!isGoodEvent<true>(collision)) {
        continue;
      }
      registry.fill(HIST("ReAssocMC/EvtGenRecReassoc"), 3., crec, occ);
      registry.fill(HIST("ReAssocMC/hReAssocMCEventStatus"), static_cast<int>(ReAssocMCEventStatus::kEvtReAsSelected), crec, occ);

      if (!collision.has_mcCollision()) {
        continue;
      }

      int64_t recCollId = collision.globalIndex();
      auto itMC = mapMcCollIdPerRecColl.find(recCollId);
      if (itMC == mapMcCollIdPerRecColl.end()) {
        nNoMC++;
        continue;
      }

      registry.fill(HIST("ReAssocMC/hReAssocMCEventStatus"), static_cast<int>(ReAssocMCEventStatus::kEvtReAsHasMcColl), crec, occ);
      auto mcColl = collision.template mcCollision_as<aod::McCollisions>();
      if (gConf.cfgRemoveSplitVertex && (!setRecCollSel.contains(collision.globalIndex()))) {
        continue;
      }
      registry.fill(HIST("ReAssocMC/hReAssocMCEventStatus"), static_cast<int>(ReAssocMCEventStatus::kEvtReAsSplitVtxRemoved), crec, occ);
      if (eventCuts.useZVtxCutMC && (std::abs(mcColl.posZ()) >= eventCuts.maxZvtx)) {
        continue;
      }
      registry.fill(HIST("ReAssocMC/hReAssocMCEventStatus"), static_cast<int>(ReAssocMCEventStatus::kEvtReAsZVtxCutMC), crec, occ);

      for (auto const& atrack : besttracks) {
        registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kTrkReAssocAll), crec, occ);
        if (!isBestTrackSelected<true>(atrack)) {
          continue;
        }
        registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kTrkBestSel), crec, occ);

        const float bestDcaZ = getDCAz(atrack);
        auto itrack = atrack.template mfttrack_as<MftTracksWCollsMC>();
        if (!isTrackSelected<true>(itrack)) {
          continue;
        }
        registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kTrkSel), crec, occ);
        if (!itrack.has_collision()) {
          continue;
        }
        registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kTrkHasColl), crec, occ);
        if (gConf.cfgRemoveReassigned) {
          if (itrack.collisionId() != atrack.bestCollisionId()) {
            continue;
          }
        }
        registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kTrkReassignedRemoved), crec, occ);

        auto ids = atrack.compatibleCollIds();
        if (gConf.cfgRemoveTrivialAssoc) {
          if (ids.empty() || (ids.size() == 1 && itrack.collisionId() == ids[0])) {
            continue;
          }
        }
        if (gConf.cfgRemoveAmbiguousTracks && (ids.size() != 1)) {
          continue;
        }

        hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkRec)]->Fill(itrack.pt(), itrack.eta(), 0., 0., 0., crec, occ);
        hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkRec)]->Fill(itrack.pt(), itrack.eta(), atrack.bestDCAXY(), bestDcaZ, 0., 0., crec, occ);

        if (itrack.collisionId() >= 0 && itrack.has_mcParticle() && itrack.mcMask() == trackCuts.selMcMask) {

          auto particle = itrack.template mcParticle_as<aod::McParticles>();
          if (!isChrgParticle(particle.pdgCode())) {
            continue;
          }
          if (gConf.cfgUseParticleSel && !isParticleSelected(particle)) {
            continue;
          }
          // auto collision = itrack.template collision_as<aod::McCollisionLabels>();
          if (eventCuts.useZDiffCut) {
            if (std::abs(collision.posZ() - atrack.mcParticle().mcCollision().posZ()) > eventCuts.maxZvtxDiff) {
              continue;
            }
          }

          float deltaX = -999.f;
          float deltaY = -999.f;
          float deltaZ = -999.f;

          auto vtxXtruth = atrack.mcParticle().mcCollision().posX();
          auto vtxYtruth = atrack.mcParticle().mcCollision().posY();
          auto vtxZtruth = atrack.mcParticle().mcCollision().posZ();

          const int bestRecColl = atrack.bestCollisionId();
          auto itMapVtxXrec = mapVtxXrec.find(bestRecColl);
          auto itMapVtxYrec = mapVtxYrec.find(bestRecColl);
          auto itMapVtxZrec = mapVtxZrec.find(bestRecColl);
          auto itMapMcCollIdPerRecColl = mapMcCollIdPerRecColl.find(bestRecColl);

          if (itMapVtxXrec == mapVtxXrec.end()) {
            continue;
          }
          if (itMapVtxYrec == mapVtxYrec.end()) {
            continue;
          }
          if (itMapVtxZrec == mapVtxZrec.end()) {
            continue;
          }
          if (itMapMcCollIdPerRecColl == mapMcCollIdPerRecColl.end()) {
            continue;
          }
          const float vtxXbest = itMapVtxXrec->second;
          const float vtxYbest = itMapVtxYrec->second;
          const float vtxZbest = itMapVtxZrec->second;
          // LOGP(info, "\t ---> \t .... \t vtxZrec: {} - collision.posZ(): {}", vtxZrec, collision.posZ());
          const float mcCollIdRec = itMapMcCollIdPerRecColl->second;
          // LOGP(info, "\t ---> \t .... \t mcCollIdRec: {} - bestMCCol: {}", mcCollIdRec, bestMCCol);

          deltaX = vtxXbest - vtxXtruth;
          deltaY = vtxYbest - vtxYtruth;
          deltaZ = vtxZbest - vtxZtruth;

          const auto dcaXtruth(particle.vx() - particle.mcCollision().posX());
          const auto dcaYtruth(particle.vy() - particle.mcCollision().posY());
          const auto dcaZtruth(particle.vz() - particle.mcCollision().posZ());
          auto dcaXYtruth = std::sqrt(dcaXtruth * dcaXtruth + dcaYtruth * dcaYtruth);

          registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kTrkHasMcPart), crec, occ);

          if (!ids.empty()) {
            if (collision.has_mcCollision() && mcCollIdRec == particle.mcCollisionId()) { // good coll
              if (isRecInCompColl(atrack)) {                                              // coll vertex is among compatible colls
                registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kTrkIdGt0), crec, occ);
                hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkIdGt0)]->Fill(itrack.pt(), itrack.eta(), deltaX, deltaY, deltaZ, crec, occ);
                hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkIdGt0)]->Fill(itrack.pt(), itrack.eta(), atrack.bestDCAXY(), bestDcaZ, dcaXYtruth, dcaZtruth, crec, occ);
              }
            }
            if (ids.size() == 1) {
              registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kTrkNonAmb), crec, occ);
              hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmb)]->Fill(itrack.pt(), itrack.eta(), deltaX, deltaY, deltaZ, crec, occ);
              hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmb)]->Fill(itrack.pt(), itrack.eta(), atrack.bestDCAXY(), bestDcaZ, dcaXYtruth, dcaZtruth, crec, occ);
              if (collision.mcCollisionId() == particle.mcCollisionId()) {
                registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kTrkNonAmbGood), crec, occ);
                hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbGood)]->Fill(itrack.pt(), itrack.eta(), deltaX, deltaY, deltaZ, crec, occ);
                hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbGood)]->Fill(itrack.pt(), itrack.eta(), atrack.bestDCAXY(), bestDcaZ, dcaXYtruth, dcaZtruth, crec, occ);
              } else {
                registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kTrkNonAmbBad), crec, occ);
                hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbBad)]->Fill(itrack.pt(), itrack.eta(), deltaX, deltaY, deltaZ, crec, occ);
                hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbBad)]->Fill(itrack.pt(), itrack.eta(), atrack.bestDCAXY(), bestDcaZ, dcaXYtruth, dcaZtruth, crec, occ);
              }
              if (itrack.collisionId() == ids[0]) { // non ambiguous
                registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kTrkNonAmbID), crec, occ);
                hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbID)]->Fill(itrack.pt(), itrack.eta(), deltaX, deltaY, deltaZ, crec, occ);
                hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbID)]->Fill(itrack.pt(), itrack.eta(), atrack.bestDCAXY(), bestDcaZ, dcaXYtruth, dcaZtruth, crec, occ);
                if (collision.mcCollisionId() == particle.mcCollisionId()) {
                  registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kTrkNonAmbIDGood), crec, occ);
                  hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbIDGood)]->Fill(itrack.pt(), itrack.eta(), deltaX, deltaY, deltaZ, crec, occ);
                  hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbIDGood)]->Fill(itrack.pt(), itrack.eta(), atrack.bestDCAXY(), bestDcaZ, dcaXYtruth, dcaZtruth, crec, occ);
                } else {
                  registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kTrkNonAmbIDBad), crec, occ);
                  hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbIDBad)]->Fill(itrack.pt(), itrack.eta(), deltaX, deltaY, deltaZ, crec, occ);
                  hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbIDBad)]->Fill(itrack.pt(), itrack.eta(), atrack.bestDCAXY(), bestDcaZ, dcaXYtruth, dcaZtruth, crec, occ);
                }
              } else if (itrack.collisionId() != ids[0]) { // ambiguous extra
                registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kTrkAmbID), crec, occ);
                hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkAmbID)]->Fill(itrack.pt(), itrack.eta(), deltaX, deltaY, deltaZ, crec, occ);
                hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkAmbID)]->Fill(itrack.pt(), itrack.eta(), atrack.bestDCAXY(), bestDcaZ, dcaXYtruth, dcaZtruth, crec, occ);
                if (collision.mcCollisionId() == particle.mcCollisionId()) {
                  registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kTrkAmbIDGood), crec, occ);
                  hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkAmbIDGood)]->Fill(itrack.pt(), itrack.eta(), deltaX, deltaY, deltaZ, crec, occ);
                  hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkAmbIDGood)]->Fill(itrack.pt(), itrack.eta(), atrack.bestDCAXY(), bestDcaZ, dcaXYtruth, dcaZtruth, crec, occ);
                } else {
                  registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kTrkAmbIDBad), crec, occ);
                  hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkAmbIDBad)]->Fill(itrack.pt(), itrack.eta(), deltaX, deltaY, deltaZ, crec, occ);
                  hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkAmbIDBad)]->Fill(itrack.pt(), itrack.eta(), atrack.bestDCAXY(), bestDcaZ, dcaXYtruth, dcaZtruth, crec, occ);
                }
              } else { // non ambiguous (extra)
                registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kTrkNonAmbIDExtra), crec, occ);
                hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbIDExtra)]->Fill(itrack.pt(), itrack.eta(), deltaX, deltaY, deltaZ, crec, occ);
                hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbIDExtra)]->Fill(itrack.pt(), itrack.eta(), atrack.bestDCAXY(), bestDcaZ, dcaXYtruth, dcaZtruth, crec, occ);
                if (collision.mcCollisionId() == particle.mcCollisionId()) {
                  registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kTrkNonAmbIDExtraGood), crec, occ);
                  hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbIDExtraGood)]->Fill(itrack.pt(), itrack.eta(), deltaX, deltaY, deltaZ, crec, occ);
                  hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbIDExtraGood)]->Fill(itrack.pt(), itrack.eta(), atrack.bestDCAXY(), bestDcaZ, dcaXYtruth, dcaZtruth, crec, occ);
                } else {
                  registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kTrkNonAmbIDExtraBad), crec, occ);
                  hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbIDExtraBad)]->Fill(itrack.pt(), itrack.eta(), deltaX, deltaY, deltaZ, crec, occ);
                  hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbIDExtraBad)]->Fill(itrack.pt(), itrack.eta(), atrack.bestDCAXY(), bestDcaZ, dcaXYtruth, dcaZtruth, crec, occ);
                }
              }
            } else {                                                  // ambiguous
              if (itrack.collisionId() != atrack.bestCollisionId()) { // re-associated
                registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kReAssoc), crec, occ);
                hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kReAssoc)]->Fill(itrack.pt(), itrack.eta(), deltaX, deltaY, deltaZ, crec, occ);
                hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kReAssoc)]->Fill(itrack.pt(), itrack.eta(), atrack.bestDCAXY(), bestDcaZ, dcaXYtruth, dcaZtruth, crec, occ);
                if (collision.has_mcCollision() && mcCollIdRec == particle.mcCollisionId()) { // good coll
                  registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kReAssocGood), crec, occ);
                  hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kReAssocGood)]->Fill(itrack.pt(), itrack.eta(), deltaX, deltaY, deltaZ, crec, occ);
                  hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kReAssocGood)]->Fill(itrack.pt(), itrack.eta(), atrack.bestDCAXY(), bestDcaZ, dcaXYtruth, dcaZtruth, crec, occ);
                  if (isRecInCompColl(atrack)) { // coll vertex is among compatible colls
                    registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kReAssocGoodIsCompTrue), crec, occ);
                    hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kReAssocGoodIsCompTrue)]->Fill(itrack.pt(), itrack.eta(), deltaX, deltaY, deltaZ, crec, occ);
                    hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kReAssocGoodIsCompTrue)]->Fill(itrack.pt(), itrack.eta(), atrack.bestDCAXY(), bestDcaZ, dcaXYtruth, dcaZtruth, crec, occ);
                  } else {
                    registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kReAssocGoodIsCompFalse), crec, occ);
                    hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kReAssocGoodIsCompFalse)]->Fill(itrack.pt(), itrack.eta(), deltaX, deltaY, deltaZ, crec, occ);
                    hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kReAssocGoodIsCompFalse)]->Fill(itrack.pt(), itrack.eta(), atrack.bestDCAXY(), bestDcaZ, dcaXYtruth, dcaZtruth, crec, occ);
                  }
                } else {
                  registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kReAssocBad), crec, occ);
                  hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kReAssocBad)]->Fill(itrack.pt(), itrack.eta(), deltaX, deltaY, deltaZ, crec, occ);
                  hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kReAssocBad)]->Fill(itrack.pt(), itrack.eta(), atrack.bestDCAXY(), bestDcaZ, dcaXYtruth, dcaZtruth, crec, occ);
                  if (isRecInCompColl(atrack)) { // coll vertex is among compatible colls
                    registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kReAssocBadIsCompTrue), crec, occ);
                    hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kReAssocBadIsCompTrue)]->Fill(itrack.pt(), itrack.eta(), deltaX, deltaY, deltaZ, crec, occ);
                    hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kReAssocBadIsCompTrue)]->Fill(itrack.pt(), itrack.eta(), atrack.bestDCAXY(), bestDcaZ, dcaXYtruth, dcaZtruth, crec, occ);
                  } else {
                    registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kReAssocBadIsCompFalse), crec, occ);
                    hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kReAssocBadIsCompFalse)]->Fill(itrack.pt(), itrack.eta(), deltaX, deltaY, deltaZ, crec, occ);
                    hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kReAssocBadIsCompFalse)]->Fill(itrack.pt(), itrack.eta(), atrack.bestDCAXY(), bestDcaZ, dcaXYtruth, dcaZtruth, crec, occ);
                  }
                }
              } else { // associated
                registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kAssoc), crec, occ);
                hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kAssoc)]->Fill(itrack.pt(), itrack.eta(), deltaX, deltaY, deltaZ, crec, occ);
                hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kAssoc)]->Fill(itrack.pt(), itrack.eta(), atrack.bestDCAXY(), bestDcaZ, dcaXYtruth, dcaZtruth, crec, occ);
                if (collision.has_mcCollision() && mcCollIdRec == particle.mcCollisionId()) { // good coll
                  registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kAssocGood), crec, occ);
                  hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kAssocGood)]->Fill(itrack.pt(), itrack.eta(), deltaX, deltaY, deltaZ, crec, occ);
                  hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kAssocGood)]->Fill(itrack.pt(), itrack.eta(), atrack.bestDCAXY(), bestDcaZ, dcaXYtruth, dcaZtruth, crec, occ);
                  if (isRecInCompColl(atrack)) { // coll vertex is among compatible colls
                    registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kAssocGoodIsCompTrue), crec, occ);
                    hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kAssocGoodIsCompTrue)]->Fill(itrack.pt(), itrack.eta(), deltaX, deltaY, deltaZ, crec, occ);
                    hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kAssocGoodIsCompTrue)]->Fill(itrack.pt(), itrack.eta(), atrack.bestDCAXY(), bestDcaZ, dcaXYtruth, dcaZtruth, crec, occ);
                  } else {
                    registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kAssocGoodIsCompFalse), crec, occ);
                    hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kAssocGoodIsCompFalse)]->Fill(itrack.pt(), itrack.eta(), deltaX, deltaY, deltaZ, crec, occ);
                    hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kAssocGoodIsCompFalse)]->Fill(itrack.pt(), itrack.eta(), atrack.bestDCAXY(), bestDcaZ, dcaXYtruth, dcaZtruth, crec, occ);
                  }
                } else {
                  registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kAssocBad), crec, occ);
                  hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kAssocBad)]->Fill(itrack.pt(), itrack.eta(), deltaX, deltaY, deltaZ, crec, occ);
                  hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kAssocBad)]->Fill(itrack.pt(), itrack.eta(), atrack.bestDCAXY(), bestDcaZ, dcaXYtruth, dcaZtruth, crec, occ);
                  if (isRecInCompColl(atrack)) { // coll vertex is among compatible colls
                    registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kAssocBadIsCompTrue), crec, occ);
                    hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kAssocBadIsCompTrue)]->Fill(itrack.pt(), itrack.eta(), deltaX, deltaY, deltaZ, crec, occ);
                    hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kAssocBadIsCompTrue)]->Fill(itrack.pt(), itrack.eta(), atrack.bestDCAXY(), bestDcaZ, dcaXYtruth, dcaZtruth, crec, occ);
                  } else {
                    registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kAssocBadIsCompFalse), crec, occ);
                    hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kAssocBadIsCompFalse)]->Fill(itrack.pt(), itrack.eta(), deltaX, deltaY, deltaZ, crec, occ);
                    hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kAssocBadIsCompFalse)]->Fill(itrack.pt(), itrack.eta(), atrack.bestDCAXY(), bestDcaZ, dcaXYtruth, dcaZtruth, crec, occ);
                  }
                }
              }
            }
          }
        }
      }
    }
    LOG(info) << "No MC: " << nNoMC;
  }
  PROCESS_SWITCH(QaReassocMFTPbPb, processReAssoc3dMC, "Process re-association information based on BestCollisionsFwd3d (3D) table", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<QaReassocMFTPbPb>(cfgc)};
}
