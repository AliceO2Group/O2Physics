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

#include "Functions.h"
#include "Index.h"
#include "bestCollisionTable.h"

#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/fwdtrackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/McCollisionExtra.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/GeomConstants.h"
#include "CommonConstants/MathConstants.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/Propagator.h"
#include "Field/MagneticField.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/RuntimeError.h"
#include "Framework/runDataProcessing.h"
#include "MathUtils/Utils.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "ReconstructionDataFormats/TrackFwd.h"

#include "Math/MatrixFunctions.h"
#include "Math/SMatrix.h"
#include "TGeoGlobalMagField.h"
#include "TMCProcess.h"
#include "TPDGCode.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>

using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<double, 5>;

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
auto static constexpr Czero = 0.f;
auto static constexpr Cninety = 90.f;
auto static constexpr ConeHeighty = 180.f;
auto static constexpr CminAccFT0A = 3.5f;
auto static constexpr CmaxAccFT0A = 4.9f;
auto static constexpr CminAccFT0C = -3.3f;
auto static constexpr CmaxAccFT0C = -2.1f;

constexpr int CevtSel = 15;
constexpr int CtrkSel = 7;
constexpr int CtrkTrkBestSel = 6;
constexpr int CambTrkType = 7;
constexpr int CselAmbTrkTypeAssocFlag = 16;
constexpr int CtrackToCollEvtType = 5;
constexpr int CreassocVtxType = 18;
constexpr int CevtReAsReAssocMCEventStatus = 5;
constexpr int CreAssocMCTrackStatus = 29;

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
  kSel = 0,
  kSelGoodVtxTrue,
  kSelGoodVtxBad,
  kSelNonAmbAll,
  kSelNonAmbGoodVtxTrue,
  kSelNonAmbGoodVtxBad,
  kSelNonAmbSameAll,
  kSelNonAmbSameGoodVtxTrue,
  kSelNonAmbSameGoodVtxBad,
  kSelAmbAll,
  kSelAmbGoodVtxTrue,
  kSelAmbGoodVtxBad,
  kSelAmbGt1All,
  kSelAmbGt1GoodVtxTrue,
  kSelAmbGt1GoodVtxBad,
  kSelOrphanNull,
  nSelAmbTrkTypeAssocFlag
};

enum class TrackToCollEvtType {
  kAllRecColl = 0,
  kIsInelGt0wMft,
  kEvtSel,
  kBestCollIdx,
  kIsMcColl,
  nTrackToCollEvtType
};

enum class VertexStatusMC {
  kNull = 0,
  kGood,
  kBad
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
  kTrkNonAmbAll,
  kTrkNonAmbGood,
  kTrkNonAmbBad,
  kTrkAmbAll,
  kTrkAmbGood,
  kTrkAmbBad,
  kTrkNonAmbAllE,
  kTrkNonAmbGoodE,
  kTrkNonAmbBadE,
  kAssoc,
  kAssocGood,
  kAssocGoodIsCompTrue,
  kAssocGoodIsCompFalse,
  kAssocBad,
  kAssocBadIsCompTrue,
  kAssocBadIsCompFalse,
  kReAssoc,
  kReAssocGood,
  kReAssocGoodIsCompTrue,
  kReAssocGoodIsCompFalse,
  kReAssocBad,
  kReAssocBadIsCompTrue,
  kReAssocBadIsCompFalse,
  nReAssocMCTrackStatus
};

enum class HistStatusReAssocVtx {
  kTrkNonAmbAll = 0,
  kTrkNonAmbGood,
  kTrkNonAmbBad,
  kTrkAmbAll,
  kTrkAmbGood,
  kTrkAmbBad,
  kTrkNonAmbAllE,
  kTrkNonAmbGoodE,
  kTrkNonAmbBadE,
  kAssoc,
  kAssocGood,
  kAssocGoodIsCompTrue,
  kAssocGoodIsCompFalse,
  kAssocBad,
  kAssocBadIsCompTrue,
  kAssocBadIsCompFalse,
  kReAssoc,
  kReAssocGood,
  kReAssocGoodIsCompTrue,
  kReAssocGoodIsCompFalse,
  kReAssocBad,
  kReAssocBadIsCompTrue,
  kReAssocBadIsCompFalse,
  nHistStatusReAssocVtx
};

std::unordered_map<int64_t, float> mapVtxXrec;
std::unordered_map<int64_t, float> mapVtxYrec;
std::unordered_map<int64_t, float> mapVtxZrec;
std::unordered_map<int64_t, float> mapVtxXgen;
std::unordered_map<int64_t, float> mapVtxYgen;
std::unordered_map<int64_t, float> mapVtxZgen;
std::unordered_map<int64_t, float> mapMcCollIdPerRecColl;

struct DndetaMFTPbPb {
  SliceCache cache;

  std::array<std::shared_ptr<THnSparse>, 4> hCollAssoc;
  std::array<std::shared_ptr<THnSparse>, 23> hReAssocVtxRes;
  std::array<std::shared_ptr<THnSparse>, 23> hReAssocDCA;
  // std::array<std::shared_ptr<THnSparse>, 23> hReAssocDCAPrim;
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

  HistogramRegistry registryAlign{"registryAlign", {}};
  std::array<std::string, 4> mftQuadrant = {"Q0", "Q1", "Q2", "Q3"};
  std::array<std::array<std::array<std::unordered_map<std::string, o2::framework::HistPtr>, 3>, 4>, 2> hAlignment;

  struct : ConfigurableGroup {
    Configurable<bool> cfgDoIR{"cfgDoIR", false, "Flag to retrieve Interaction rate from CCDB"};
    Configurable<bool> cfgUseIRCut{"cfgUseIRCut", false, "Flag to cut on IR rate"};
    Configurable<bool> cfgIRCrashOnNull{"cfgIRCrashOnNull", false, "Flag to avoid CTP RateFetcher crash"};
    Configurable<std::string> cfgIRSource{"cfgIRSource", "ZNC hadronic", "Estimator of the interaction rate (Pb-Pb: ZNC hadronic)"};
    Configurable<bool> cfgUseTrackSel{"cfgUseTrackSel", false, "Flag to apply track selection"};
    Configurable<bool> cfgUseParticleSel{"cfgUseParticleSel", false, "Flag to apply particle selection"};
    Configurable<bool> cfgUsePrimaries{"cfgUsePrimaries", false, "Select primary particles"};
    Configurable<bool> cfgRemoveTrivialAssoc{"cfgRemoveTrivialAssoc", false, "Skip trivial associations"};
    Configurable<bool> cfgRemoveAmbiguousTracks{"cfgRemoveAmbiguousTracks", false, "Remove ambiguous tracks"};
    Configurable<bool> cfgRemoveOrphanTracks{"cfgRemoveOrphanTracks", true, "Remove orphan tracks"};
    Configurable<bool> cfgRemoveReassigned{"cfgRemoveReassigned", false, "Remove reassgined tracks"};
    Configurable<bool> cfgRemoveSplitVertex{"cfgRemoveSplitVertex", true, "Remove split vertices"};
    Configurable<bool> cfgUseTrackParExtra{"cfgUseTrackParExtra", false, "Use table with refitted track parameters"};
    Configurable<bool> cfgUseInelgt0{"cfgUseInelgt0", false, "Use INEL > 0 condition"};
    Configurable<bool> cfgUseInelgt0wMFT{"cfgUseInelgt0wMFT", false, "Use INEL > 0 condition with MFT acceptance"};
    Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
    Configurable<uint> cfgDCAtype{"cfgDCAtype", 2, "DCA coordinate type [0: DCA-X, 1: DCA-Y, 2: DCA-XY]"};
  } gConf;

  struct : ConfigurableGroup {
    ConfigurableAxis interactionRateBins{"interactionRateBins", {500, 0, 50}, "Binning for the interaction rate (kHz)"};
    ConfigurableAxis occupancyBins{"occupancyBins", {VARIABLE_WIDTH, 0.0f, 250.0f, 500.0f, 750.0f, 1000.0f, 1500.0f, 2000.0f, 3000.0f, 4500.0f, 6000.0f, 8000.0f, 10000.0f, 50000.0f}, "Occupancy"};
    ConfigurableAxis centralityBins{"centralityBins", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100}, "Centrality"};
    ConfigurableAxis irBins{"irBins", {500, 0, 50}, "Interaction rate (kHz)"};
    ConfigurableAxis pvBins{"pvBins", {501, -0.5, 500.5}, ""};
    ConfigurableAxis fv0aMultBins{"fv0aMultBins", {501, -0.5, 500.5}, ""};
    ConfigurableAxis ft0aMultBins{"ft0aMultBins", {501, -0.5, 500.5}, ""};
    ConfigurableAxis ft0cMultBins{"ft0cMultBins", {501, -0.5, 500.5}, ""};
    ConfigurableAxis ptBins{"ptBins", {101, -0.5, 10.5}, "pT binning (GeV/c)"};
    ConfigurableAxis multBins{"multBins", {701, -0.5, 700.5}, "Multiplicity binning"};
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
    Configurable<bool> useChi2Cut{"useChi2Cut", false, "use track chi2 cut"};
    Configurable<float> maxChi2NCl{"maxChi2NCl", 1000.f, "maximum chi2 per MFT clusters"};
    Configurable<bool> usePtCut{"usePtCut", false, "use track pT cut"};
    Configurable<float> minPt{"minPt", 0., "minimum pT of the MFT tracks"};
    Configurable<bool> requireCA{"requireCA", false, "Use Cellular Automaton track-finding algorithm"};
    Configurable<float> maxDCAxy{"maxDCAxy", 0.01f, "Cut on dca XY"};
    Configurable<float> maxDCAz{"maxDCAz", 0.01f, "Cut on dca Z"};
  } trackCuts;

  struct : ConfigurableGroup {
    Configurable<float> maxZvtx{"maxZvtx", 20.0f, "maximum cut on z-vtx (cm)"};
    Configurable<float> minZvtx{"minZvtx", -20.0f, "minimum cut on z-vtx (cm)"};
    Configurable<bool> useZDiffCut{"useZDiffCut", false, "use Zvtx reco-mc diff. cut"};
    Configurable<float> maxZvtxDiff{"maxZvtxDiff", 1.0f, "max allowed Z vtx difference for reconstruced collisions (cm)"};
    Configurable<bool> useZVtxCutMC{"useZVtxCutMC", false, "use Zvtx cut in MC"};
    Configurable<bool> requireIsGoodZvtxFT0VsPV{"requireIsGoodZvtxFT0VsPV", true, "require events with PV position along z consistent (within 1 cm) between PV reconstructed using tracks and PV using FT0 A-C time difference"};
    Configurable<bool> requireRejectSameBunchPileup{"requireRejectSameBunchPileup", true, "reject collisions in case of pileup with another collision in the same foundBC"};
    Configurable<bool> requireNoCollInTimeRangeStrict{"requireNoCollInTimeRangeStrict", false, " requireNoCollInTimeRangeStrict"};
    Configurable<bool> requireNoCollInRofStrict{"requireNoCollInRofStrict", true, "requireNoCollInRofStrict"};
    Configurable<bool> requireNoCollInRofStandard{"requireNoCollInRofStandard", false, "requireNoCollInRofStandard"};
    Configurable<bool> requireNoHighMultCollInPrevRof{"requireNoHighMultCollInPrevRof", false, "requireNoHighMultCollInPrevRof"};
    Configurable<bool> requireNoCollInTimeRangeStd{"requireNoCollInTimeRangeStd", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 10 microseconds"};
    Configurable<bool> requireNoCollInTimeRangeNarrow{"requireNoCollInTimeRangeNarrow", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 10 microseconds"};
    Configurable<uint> occupancyEstimator{"occupancyEstimator", 1, "Occupancy estimator: 1 = trackOccupancyInTimeRange, 2 = ft0cOccupancyInTimeRange"};
    Configurable<float> minOccupancy{"minOccupancy", -1, "minimum occupancy from neighbouring collisions"};
    Configurable<float> maxOccupancy{"maxOccupancy", -1, "maximum occupancy from neighbouring collisions"};
    Configurable<float> minIR{"minIR", -1, "minimum IR (kHz) collisions"};
    Configurable<float> maxIR{"maxIR", -1, "maximum IR (kHz) collisions"};
  } eventCuts;

  Service<o2::framework::O2DatabasePDG> pdg;
  Service<ccdb::BasicCCDBManager> ccdb;
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
  TH2* gCurrentHadronicRate;
  RCTFlagsChecker rctChecker;
  RCTFlagsChecker rctCheckerExtra{kFT0Bad, kITSBad, kTPCBadTracking, kMFTBad};

  float bZ = 0;                                          // Magnetic field for MFT
  static constexpr double CcenterMFT[3] = {0, 0, -61.4}; // Field at center of MFT
  float mZShift = 0;                                     // z-vertex shift

  o2::parameters::GRPMagField* grpmag = nullptr;

  std::vector<int> ambiguousTrkIds;
  std::vector<int> reassignedTrkIds;
  std::vector<int> ambiguousTrkIdsMC;
  std::vector<int> reassignedTrkIdsMC;

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

    rctChecker.init(rctCuts.cfgEvtRCTFlagCheckerLabel, rctCuts.cfgEvtRCTFlagCheckerZDCCheck, rctCuts.cfgEvtRCTFlagCheckerLimitAcceptAsBad);

    ccdb->setURL(ccdbUrl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(ccdbNoLaterThan.value);
    ccdb->setFatalWhenNull(false);

    if (static_cast<int>(doprocessDataInclusive) +
          static_cast<int>(doprocessDatawBestTracksInclusive) >
        1) {
      LOGP(fatal,
           "Either processDataInclusive OR "
           "processDatawBestTracksInclusive should be enabled!");
    }
    if ((static_cast<int>(doprocessDataCentFT0C) +
           static_cast<int>(doprocessDatawBestTracksCentFT0C) >
         1) ||
        (static_cast<int>(doprocessDataCentFT0CVariant1) +
           static_cast<int>(doprocessDatawBestTracksCentFT0CVariant1) >
         1) ||
        (static_cast<int>(doprocessDataCentFT0M) +
           static_cast<int>(doprocessDatawBestTracksCentFT0M) >
         1) ||
        (static_cast<int>(doprocessDataCentNGlobal) +
           static_cast<int>(doprocessDatawBestTracksCentNGlobal) >
         1) ||
        (static_cast<int>(doprocessDataCentMFT) +
           static_cast<int>(doprocessDatawBestTracksCentMFT) >
         1)) {
      LOGP(fatal,
           "Either processDataCent[ESTIMATOR] OR "
           "processDatawBestTracksCent[ESTIMATOR] should "
           "be enabled!");
    }
    if (static_cast<int>(doprocessMCInclusive) +
          static_cast<int>(doprocessMCwBestTracksInclusive) >
        1) {
      LOGP(fatal,
           "Either processMCInclusive OR processMCwBestTracksInclusive "
           "should be enabled!");
    }
    if ((static_cast<int>(doprocessMCCentFT0C) +
           static_cast<int>(doprocessMCwBestTracksCentFT0C) >
         1) ||
        (static_cast<int>(doprocessMCCentFT0CVariant1) +
           static_cast<int>(doprocessMCwBestTracksCentFT0CVariant1) >
         1) ||
        (static_cast<int>(doprocessMCCentFT0M) +
           static_cast<int>(doprocessMCwBestTracksCentFT0M) >
         1) ||
        (static_cast<int>(doprocessMCCentNGlobal) +
           static_cast<int>(doprocessMCwBestTracksCentNGlobal) >
         1) ||
        (static_cast<int>(doprocessMCCentMFT) +
           static_cast<int>(doprocessMCwBestTracksCentMFT) >
         1)) {
      LOGP(fatal,
           "Either processMCCent[ESTIMATOR] OR "
           "processMCwBestTracksCent[ESTIMATOR] should "
           "be enabled!");
    }

    // registry.add("Events/hRCTSel", "Event accepted if RCT not selected;RCT Status;", {HistType::kTH1F, {{2, 0.5, 2.5}}});
    // auto hrctSel = registry.get<TH1>(HIST("Events/hRCTSel"));
    // auto* x = hrctSel->GetXaxis();
    // x->SetBinLabel(1, "All");
    // x->SetBinLabel(2, "kFT0Bad");

    registry.add("Events/hEvtSel", "Number of events; Cut; #Evt Passed Cut", {HistType::kTH1F, {{static_cast<int>(EvtSel::nEvtSel), -0.5, +static_cast<int>(EvtSel::nEvtSel) - 0.5}}});
    std::string labelEvtSel[CevtSel];
    labelEvtSel[static_cast<int>(EvtSel::evtAll)] = "All coll.";
    labelEvtSel[static_cast<int>(EvtSel::evtSel)] = "Sel 8";
    labelEvtSel[static_cast<int>(EvtSel::evtIsGoodZvtx)] = "kIsGoodZvtxFT0vsPV";
    labelEvtSel[static_cast<int>(EvtSel::evtNoSameBunchPileup)] = "NoSameBunchPileup";
    labelEvtSel[static_cast<int>(EvtSel::evtZvtxCut)] = "Z-vtx cut";
    labelEvtSel[static_cast<int>(EvtSel::evtNoCollInTimeRangeStd)] = "kNoCollInTimeRangeStd";
    labelEvtSel[static_cast<int>(EvtSel::evtNoCollInTimeRangeNarrow)] = "kNoCollInTimeRangeNarrow";
    labelEvtSel[static_cast<int>(EvtSel::evtNoCollInTimeRangeStrict)] = "kNoCollInTimeRangeStrict";
    labelEvtSel[static_cast<int>(EvtSel::evtNoCollInRofStrict)] = "kNoCollInRofStrict";
    labelEvtSel[static_cast<int>(EvtSel::evtNoCollInRofStandard)] = "kNoCollInRofStandard";
    labelEvtSel[static_cast<int>(EvtSel::evtNoHighMultCollInPrevRof)] = "kNoHighMultCollInPrevRof";
    labelEvtSel[static_cast<int>(EvtSel::evtBelowMinOccup)] = "Below min occup.";
    labelEvtSel[static_cast<int>(EvtSel::evtAboveMaxOccup)] = "Above max occup.";
    labelEvtSel[static_cast<int>(EvtSel::evtRCTFlagChecker)] = "RCT Flag Checker";
    labelEvtSel[static_cast<int>(EvtSel::evtRCTFlagCheckerExtra)] = "RCT Flag Checker Extra";
    registry.get<TH1>(HIST("Events/hEvtSel"))->SetMinimum(0.1);
    for (int iBin = 0; iBin < static_cast<int>(EvtSel::nEvtSel); iBin++) {
      registry.get<TH1>(HIST("Events/hEvtSel"))->GetXaxis()->SetBinLabel(iBin + 1, labelEvtSel[iBin].data());
    }

    registry.add("Tracks/hBestTrkSel", "Number of best tracks; Cut; #Tracks Passed Cut", {HistType::kTH1F, {{static_cast<int>(TrkTrkBestSel::nTrkTrkBestSel), -0.5, +static_cast<int>(TrkTrkBestSel::nTrkTrkBestSel) - 0.5}}});
    std::string labelTrkTrkBestSel[CtrkTrkBestSel];
    labelTrkTrkBestSel[static_cast<int>(TrkTrkBestSel::trkTrkBestSelAll)] = "All";
    labelTrkTrkBestSel[static_cast<int>(TrkTrkBestSel::trkTrkBestSelCollID)] = "Assigned (ID>=0)";
    labelTrkTrkBestSel[static_cast<int>(TrkTrkBestSel::trkTrkBestSelOrphan)] = "No orphans";
    labelTrkTrkBestSel[static_cast<int>(TrkTrkBestSel::trkTrkBestSelDCAxyCut)] = "DCA xy cut";
    labelTrkTrkBestSel[static_cast<int>(TrkTrkBestSel::trkTrkBestSelDCAzCut)] = "DCA z cut";
    labelTrkTrkBestSel[static_cast<int>(TrkTrkBestSel::trkTrkBestSelNumReassoc)] = "#Reassoc";
    registry.get<TH1>(HIST("Tracks/hBestTrkSel"))->SetMinimum(0.1);
    for (int iBin = 0; iBin < static_cast<int>(TrkTrkBestSel::nTrkTrkBestSel); iBin++) {
      registry.get<TH1>(HIST("Tracks/hBestTrkSel"))->GetXaxis()->SetBinLabel(iBin + 1, labelTrkTrkBestSel[iBin].data());
    }

    registry.add("Tracks/hTrkSel", "Number of tracks; Cut; #Tracks Passed Cut", {HistType::kTH1F, {{static_cast<int>(TrkSel::nTrkSel), -0.5, +static_cast<int>(TrkSel::nTrkSel) - 0.5}}});
    std::string labelTrkSel[CtrkSel];
    labelTrkSel[static_cast<int>(TrkSel::trkSelAll)] = "All";
    labelTrkSel[static_cast<int>(TrkSel::trkSelNCls)] = "Ncls";
    labelTrkSel[static_cast<int>(TrkSel::trkSelChi2Ncl)] = "Chi2";
    labelTrkSel[static_cast<int>(TrkSel::trkSelEta)] = "Eta";
    labelTrkSel[static_cast<int>(TrkSel::trkSelPhiCut)] = "Phi cut";
    labelTrkSel[static_cast<int>(TrkSel::trkSelPt)] = "Pt";
    labelTrkSel[static_cast<int>(TrkSel::trkSelCA)] = "CA";
    registry.get<TH1>(HIST("Tracks/hTrkSel"))->SetMinimum(0.1);
    for (int iBin = 0; iBin < static_cast<int>(TrkSel::nTrkSel); iBin++) {
      registry.get<TH1>(HIST("Tracks/hTrkSel"))->GetXaxis()->SetBinLabel(iBin + 1, labelTrkSel[iBin].data());
    }

    auto hBcSel = registry.add<TH1>("hBcSel", "hBcSel", HistType::kTH1F,
                                    {{3, -0.5f, +2.5f}});
    hBcSel->GetXaxis()->SetBinLabel(1, "Good BCs");
    hBcSel->GetXaxis()->SetBinLabel(2, "BCs with collisions");
    hBcSel->GetXaxis()->SetBinLabel(3, "BCs with pile-up/splitting");

    if (doprocessDataInclusive || doprocessDatawBestTracksInclusive ||
        doprocessMCInclusive || doprocessMCwBestTracksInclusive) {
      registry.add({"Events/Selection",
                    ";status;occupancy",
                    {HistType::kTH2F, {{2, 0.5, 2.5}, occupancyAxis}}});
      auto hstat = registry.get<TH2>(HIST("Events/Selection"));
      auto* x = hstat->GetXaxis();
      x->SetBinLabel(1, "All");
      x->SetBinLabel(2, "Selected");

      registry.add("Events/hInteractionRate", "; occupancy; IR (kHz)", kTH2F, {occupancyAxis, irAxis});

      registry.add({"Events/NtrkZvtx",
                    "; N_{trk}; Z_{vtx} (cm); occupancy",
                    {HistType::kTHnSparseF, {multAxis, zAxis, occupancyAxis}}});
      registry.add({"Tracks/EtaZvtx",
                    "; #eta; Z_{vtx} (cm); occupancy",
                    {HistType::kTHnSparseF, {etaAxis, zAxis, occupancyAxis}}});
      registry.add(
        {"Tracks/PhiEta",
         "; #varphi; #eta; occupancy",
         {HistType::kTHnSparseF, {phiAxis, etaAxis, occupancyAxis}}});

      qaregistry.add(
        {"Tracks/Chi2Eta",
         "; #chi^{2}; #it{#eta}; occupancy",
         {HistType::kTHnSparseF, {chiSqAxis, etaAxis, occupancyAxis}}});
      qaregistry.add({"Tracks/Chi2",
                      "; #chi^{2};",
                      {HistType::kTH2F, {chiSqAxis, occupancyAxis}}});
      qaregistry.add(
        {"Tracks/NclustersEta",
         "; nClusters; #eta; occupancy",
         {HistType::kTHnSparseF, {nclsAxis, etaAxis, occupancyAxis}}});
      qaregistry.add({"Tracks/NchSel",
                      "; N_{ch}; occupancy",
                      {HistType::kTH2F, {multAxis, occupancyAxis}}});
      qaregistry.add({"Tracks/NchBestSel",
                      "; N_{ch}; occupancy",
                      {HistType::kTH2F, {multAxis, occupancyAxis}}});

      qaregistry.add({"Tracks/TanLambda", "; TanLambda; centrality; occupancy", {HistType::kTH2F, {tanLambdaAxis, occupancyAxis}}});
      qaregistry.add({"Tracks/InvQPt", "; InvQPt; centrality; occupancy", {HistType::kTH2F, {invQPtAxis, occupancyAxis}}});
      qaregistry.add({"Tracks/Eta", "; #eta; centrality; occupancy", {HistType::kTH2F, {etaAxis, occupancyAxis}}});
      qaregistry.add({"Tracks/Phi", "; #varphi; centrality; occupancy", {HistType::kTH2F, {phiAxis, occupancyAxis}}});

      if (doprocessDatawBestTracksInclusive) {
        registry.add(
          {"Events/NtrkZvtxBest",
           "; N_{trk}; Z_{vtx} (cm); occupancy",
           {HistType::kTHnSparseF, {multAxis, zAxis, occupancyAxis}}});
        registry.add(
          {"Tracks/EtaZvtxBest",
           "; #eta; Z_{vtx} (cm); occupancy",
           {HistType::kTHnSparseF, {etaAxis, zAxis, occupancyAxis}}});
        registry.add(
          {"Tracks/PhiEtaBest",
           "; #varphi; #eta; occupancy",
           {HistType::kTHnSparseF, {phiAxis, etaAxis, occupancyAxis}}});
        qaregistry.add(
          {"Tracks/NclustersEtaBest",
           "; nClusters; #eta; occupancy",
           {HistType::kTHnSparseF, {nclsAxis, etaAxis, occupancyAxis}}});
        qaregistry.add(
          {"Tracks/DCA3d",
           "; p_{T} (GeV/c); #eta; DCA_{XY} (cm); DCA_{Z} (cm); occupancy",
           {HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, occupancyAxis}}});
        qaregistry.add(
          {"Tracks/ReTracksEtaZvtx",
           "; #eta; #it{z}_{vtx} (cm); occupancy",
           {HistType::kTHnSparseF, {etaAxis, zAxis, occupancyAxis}}});
        qaregistry.add(
          {"Tracks/ReTracksPhiEta",
           "; #varphi; #eta; occupancy",
           {HistType::kTHnSparseF, {phiAxis, etaAxis, occupancyAxis}}});
        qaregistry.add(
          {"Tracks/OrigTracksEtaZvtx",
           "; #eta; #it{z}_{vtx} (cm); occupancy",
           {HistType::kTHnSparseF, {etaAxis, zAxis, occupancyAxis}}});
        qaregistry.add(
          {"Tracks/OrigTracksPhiEta",
           "; #varphi; #eta; occupancy",
           {HistType::kTHnSparseF, {phiAxis, etaAxis, occupancyAxis}}});
        qaregistry.add(
          {"Tracks/RestTracksEtaZvtx",
           "; #eta; #it{z}_{vtx} (cm); occupancy",
           {HistType::kTHnSparseF, {etaAxis, zAxis, occupancyAxis}}});
        qaregistry.add(
          {"Tracks/RestTracksPhiEta",
           "; #varphi; #eta; occupancy",
           {HistType::kTHnSparseF, {phiAxis, etaAxis, occupancyAxis}}});
        qaregistry.add({"Tracks/TrackAmbDegree",
                        "; N_{coll}^{comp}; occupancy",
                        {HistType::kTH2F, {{51, -0.5, 50.5}, occupancyAxis}}});
        qaregistry.add({"Tracks/TanLambdaExtra", "; TanLambda; occupancy", {HistType::kTH2F, {tanLambdaAxis, occupancyAxis}}});
        qaregistry.add({"Tracks/InvQPtExtra", "; InvQPt; occupancy", {HistType::kTH2F, {invQPtAxis, occupancyAxis}}});
        qaregistry.add({"Tracks/EtaExtra", "; #eta; occupancy", {HistType::kTH2F, {etaAxis, occupancyAxis}}});
        qaregistry.add({"Tracks/PhiExtra", "; #varphi; occupancy", {HistType::kTH2F, {phiAxis, occupancyAxis}}});
      }
    }

    if (doprocessDataCentFT0C || doprocessDatawBestTracksCentFT0C ||
        doprocessMCCentFT0C || doprocessMCwBestTracksCentFT0C ||
        doprocessDataCentFT0CVariant1 ||
        doprocessDatawBestTracksCentFT0CVariant1 ||
        doprocessMCCentFT0CVariant1 || doprocessMCwBestTracksCentFT0CVariant1 ||
        doprocessDataCentFT0M || doprocessDatawBestTracksCentFT0M ||
        doprocessMCCentFT0M || doprocessMCwBestTracksCentFT0M ||
        doprocessDataCentNGlobal || doprocessDatawBestTracksCentNGlobal ||
        doprocessMCCentNGlobal || doprocessMCwBestTracksCentNGlobal ||
        doprocessDataCentMFT || doprocessDatawBestTracksCentMFT ||
        doprocessMCCentMFT || doprocessMCwBestTracksCentMFT) {
      registry.add({"Events/Centrality/Selection",
                    ";status;centrality;occupancy",
                    {HistType::kTHnSparseF,
                     {{2, 0.5, 2.5}, centralityAxis, occupancyAxis}}});
      auto hstat = registry.get<THnSparse>(HIST("Events/Centrality/Selection"));
      hstat->GetAxis(0)->SetBinLabel(1, "All");
      hstat->GetAxis(0)->SetBinLabel(2, "Selected");

      registry.add("Events/Centrality/hInteractionRate", "; centrality; occupancy; IR (kHz)", kTHnSparseF, {centralityAxis, occupancyAxis, irAxis});

      qaregistry.add({"Events/Centrality/hCent",
                      "; centrality; occupancy",
                      {HistType::kTH2F, {centralityAxis, occupancyAxis}},
                      true});
      qaregistry.add(
        {"Events/Centrality/hZvtxCent",
         "; Z_{vtx} (cm); centrality; occupancy",
         {HistType::kTHnSparseF, {zAxis, centralityAxis, occupancyAxis}}});
      registry.add({"Events/Centrality/NtrkZvtx",
                    "; N_{trk}; Z_{vtx} (cm); centrality; occupancy",
                    {HistType::kTHnSparseF,
                     {multAxis, zAxis, centralityAxis, occupancyAxis}}});
      registry.add({"Tracks/Centrality/EtaZvtx",
                    "; #eta; Z_{vtx} (cm); centrality; occupancy",
                    {HistType::kTHnSparseF,
                     {etaAxis, zAxis, centralityAxis, occupancyAxis}}});
      registry.add({"Tracks/Centrality/PhiEta",
                    "; #varphi; #eta; centrality; occupancy",
                    {HistType::kTHnSparseF,
                     {phiAxis, etaAxis, centralityAxis, occupancyAxis}}});

      qaregistry.add(
        {"Tracks/Centrality/NchSel",
         "; N_{ch}; centrality; occupancy",
         {HistType::kTHnSparseF, {multAxis, centralityAxis, occupancyAxis}}});
      qaregistry.add(
        {"Tracks/Centrality/NchBestSel",
         "; N_{ch}; centrality; occupancy",
         {HistType::kTHnSparseF, {multAxis, centralityAxis, occupancyAxis}}});
      qaregistry.add(
        {"Tracks/Centrality/Chi2Eta",
         "; #chi^{2}; #it{#eta}; centrality; occupancy",
         {HistType::kTHnSparseF,
          {chiSqAxis, etaAxis, centralityAxis, occupancyAxis}}});
      qaregistry.add({"Tracks/Centrality/Chi2",
                      "; #chi^{2}; centrality; occupancy",
                      {HistType::kTHnSparseF,
                       {chiSqAxis, centralityAxis, occupancyAxis}}});
      qaregistry.add({"Tracks/Centrality/NclustersEta",
                      "; nClusters; #eta; centrality; occupancy",
                      {HistType::kTHnSparseF,
                       {nclsAxis, etaAxis, centralityAxis, occupancyAxis}}});

      qaregistry.add({"Tracks/Centrality/TanLambda", "; TanLambda; centrality; occupancy", {HistType::kTHnSparseF, {tanLambdaAxis, centralityAxis, occupancyAxis}}});
      qaregistry.add({"Tracks/Centrality/InvQPt", "; InvQPt; centrality; occupancy", {HistType::kTHnSparseF, {invQPtAxis, centralityAxis, occupancyAxis}}});
      qaregistry.add({"Tracks/Centrality/Eta", "; #eta; centrality; occupancy", {HistType::kTHnSparseF, {etaAxis, centralityAxis, occupancyAxis}}});
      qaregistry.add({"Tracks/Centrality/Phi", "; #varphi; centrality; occupancy", {HistType::kTHnSparseF, {phiAxis, centralityAxis, occupancyAxis}}});

      if (doprocessDatawBestTracksCentFT0C ||
          doprocessDatawBestTracksCentFT0CVariant1 ||
          doprocessDatawBestTracksCentFT0M ||
          doprocessDatawBestTracksCentNGlobal ||
          doprocessDatawBestTracksCentMFT) {
        registry.add({"Events/Centrality/NtrkZvtxBest",
                      "; N_{trk}; Z_{vtx} (cm); centrality; occupancy",
                      {HistType::kTHnSparseF,
                       {multAxis, zAxis, centralityAxis, occupancyAxis}}});
        registry.add({"Tracks/Centrality/EtaZvtxBest",
                      "; #eta; Z_{vtx} (cm); centrality; occupancy",
                      {HistType::kTHnSparseF,
                       {etaAxis, zAxis, centralityAxis, occupancyAxis}}});
        registry.add({"Tracks/Centrality/PhiEtaBest",
                      "; #varphi; #eta; centrality; occupancy",
                      {HistType::kTHnSparseF,
                       {phiAxis, etaAxis, centralityAxis, occupancyAxis}}});
        qaregistry.add(
          {"Tracks/Centrality/NclustersEtaBest",
           "; nClusters; #eta; centrality; occupancy",
           {HistType::kTHnSparseF,
            {nclsAxis, etaAxis, centralityAxis, occupancyAxis}}});
        qaregistry.add({"Tracks/Centrality/TrackAmbDegree",
                        "; N_{coll}^{comp}; centrality; occupancy",
                        {HistType::kTHnSparseF,
                         {{51, -0.5, 50.5}, centralityAxis, occupancyAxis}}});
        qaregistry.add({"Tracks/Centrality/TanLambdaExtra", "; TanLambda; centrality; occupancy", {HistType::kTHnSparseF, {tanLambdaAxis, centralityAxis, occupancyAxis}}});
        qaregistry.add({"Tracks/Centrality/InvQPtExtra", "; InvQPt; centrality; occupancy", {HistType::kTHnSparseF, {invQPtAxis, centralityAxis, occupancyAxis}}});
        qaregistry.add({"Tracks/Centrality/EtaExtra", "; #eta; centrality; occupancy", {HistType::kTHnSparseF, {etaAxis, centralityAxis, occupancyAxis}}});
        qaregistry.add({"Tracks/Centrality/PhiExtra", "; #varphi; centrality; occupancy", {HistType::kTHnSparseF, {phiAxis, centralityAxis, occupancyAxis}}});
        qaregistry.add(
          {"Tracks/Centrality/DCA3d",
           "; p_{T} (GeV/c); #eta; DCA_{XY} (cm); DCA_{Z} (cm); centrality; occupancy",
           {HistType::kTHnSparseF,
            {ptAxis, etaAxis, dcaxyAxis, dcazAxis, centralityAxis, occupancyAxis}}});
        qaregistry.add({"Tracks/Centrality/ReTracksEtaZvtx",
                        "; #eta; #it{z}_{vtx} (cm); occupancy",
                        {HistType::kTHnSparseF,
                         {etaAxis, zAxis, centralityAxis, occupancyAxis}}});
        qaregistry.add({"Tracks/Centrality/ReTracksPhiEta",
                        "; #varphi; #eta; occupancy",
                        {HistType::kTHnSparseF,
                         {phiAxis, etaAxis, centralityAxis, occupancyAxis}}});
        qaregistry.add({"Tracks/Centrality/OrigTracksEtaZvtx",
                        "; #eta; #it{z}_{vtx} (cm); occupancy",
                        {HistType::kTHnSparseF,
                         {etaAxis, zAxis, centralityAxis, occupancyAxis}}});
        qaregistry.add({"Tracks/Centrality/OrigTracksPhiEta",
                        "; #varphi; #eta; occupancy",
                        {HistType::kTHnSparseF,
                         {phiAxis, etaAxis, centralityAxis, occupancyAxis}}});
        qaregistry.add({"Tracks/Centrality/RestTracksEtaZvtx",
                        "; #eta; #it{z}_{vtx} (cm); occupancy",
                        {HistType::kTHnSparseF,
                         {etaAxis, zAxis, centralityAxis, occupancyAxis}}});
        qaregistry.add({"Tracks/Centrality/RestTracksPhiEta",
                        "; #varphi; #eta; occupancy",
                        {HistType::kTHnSparseF,
                         {phiAxis, etaAxis, centralityAxis, occupancyAxis}}});
      }
    }

    if (doprocessMCInclusive || doprocessMCwBestTracksInclusive) {
      registry.add({"Events/EvtEffGen",
                    ";status;occupancy",
                    {HistType::kTH2F, {{3, 0.5, 3.5}, occupancyAxis}}});
      auto heff = registry.get<TH2>(HIST("Events/EvtEffGen"));
      auto* h = heff->GetXaxis();
      h->SetBinLabel(1, "All reconstructed");
      h->SetBinLabel(2, "Selected reconstructed");
      h->SetBinLabel(3, "All generated");

      registry.add({"Events/NtrkZvtxGen_t",
                    "; N_{trk}; Z_{vtx} (cm);",
                    {HistType::kTH2F, {multAxis, zAxis}}});
      registry.add({"Events/NtrkZvtxGen",
                    "; N_{trk}; Z_{vtx} (cm); occupancy",
                    {HistType::kTHnSparseF, {multAxis, zAxis, occupancyAxis}}});
      registry.add({"Tracks/EtaZvtxGen",
                    "; #eta; Z_{vtx} (cm); occupancy",
                    {HistType::kTHnSparseF, {etaAxis, zAxis, occupancyAxis}}});
      registry.add(
        {"Tracks/PhiEtaGen",
         "; #varphi; #eta;",
         {HistType::kTHnSparseF, {phiAxis, etaAxis, occupancyAxis}}});
      registry.add({"Tracks/EtaZvtxGen_t",
                    "; #eta; Z_{vtx} (cm);",
                    {HistType::kTH2F, {etaAxis, zAxis}}});
      registry.add({"Tracks/PhiEtaGen_t",
                    "; #varphi; #eta;",
                    {HistType::kTH2F, {phiAxis, etaAxis}}});
      qaregistry.add({"Events/NotFoundEventZvtx",
                      "; #it{z}_{vtx} (cm)",
                      {HistType::kTH1F, {zAxis}}});
      qaregistry.add({"Events/ZvtxDiff",
                      "; Z_{rec} - Z_{gen} (cm)",
                      {HistType::kTH1F, {deltaZAxis}}});
      qaregistry.add({"Events/SplitMult",
                      "; N_{gen}; #it{z}_{vtx} (cm)",
                      {HistType::kTH2F, {multAxis, zAxis}}});
    }

    if (doprocessMCCentFT0C || doprocessMCwBestTracksCentFT0C ||
        doprocessMCCentFT0CVariant1 || doprocessMCwBestTracksCentFT0CVariant1 ||
        doprocessMCCentFT0M || doprocessMCwBestTracksCentFT0M ||
        doprocessMCCentNGlobal || doprocessMCwBestTracksCentNGlobal ||
        doprocessMCCentMFT || doprocessMCwBestTracksCentMFT) {
      registry.add({"Events/Centrality/EvtEffGen",
                    ";status;centrality;occupancy",
                    {HistType::kTHnSparseF,
                     {{3, 0.5, 3.5}, centralityAxis, occupancyAxis}}});
      auto heff = registry.get<THnSparse>(HIST("Events/Centrality/EvtEffGen"));
      heff->GetAxis(0)->SetBinLabel(1, "All reconstructed");
      heff->GetAxis(0)->SetBinLabel(2, "Selected reconstructed");
      heff->GetAxis(0)->SetBinLabel(3, "All generated");

      registry.add(
        {"Events/Centrality/NtrkZvtxGen_t",
         "; N_{trk}; Z_{vtx} (cm); centrality",
         {HistType::kTHnSparseF, {multAxis, zAxis, centralityAxis}}});
      registry.add({"Events/Centrality/NtrkZvtxGen",
                    "; N_{trk}; Z_{vtx} (cm); centrality; occupancy",
                    {HistType::kTHnSparseF,
                     {multAxis, zAxis, centralityAxis, occupancyAxis}}});
      registry.add({"Events/Centrality/hRecCent",
                    "; centrality; occupancy",
                    {HistType::kTH2F, {centralityAxis, occupancyAxis}}});
      registry.add(
        {"Events/Centrality/hRecZvtxCent",
         "; Z_{vtx} (cm); centrality; occupancy",
         {HistType::kTHnSparseF, {zAxis, centralityAxis, occupancyAxis}}});
      registry.add({"Tracks/Centrality/EtaZvtxGen",
                    "; #eta; Z_{vtx} (cm); centrality; occupancy",
                    {HistType::kTHnSparseF,
                     {etaAxis, zAxis, centralityAxis, occupancyAxis}}});
      registry.add({"Tracks/Centrality/PhiEtaGen",
                    "; #varphi; #eta; centrality; occupancy",
                    {HistType::kTHnSparseF,
                     {phiAxis, etaAxis, centralityAxis, occupancyAxis}}});
      registry.add({"Tracks/Centrality/EtaZvtxGen_t",
                    "; #eta; Z_{vtx} (cm); centrality",
                    {HistType::kTHnSparseF, {etaAxis, zAxis, centralityAxis}}});
      registry.add(
        {"Tracks/Centrality/PhiEtaGen_t",
         "; #varphi; #eta; centrality",
         {HistType::kTHnSparseF, {phiAxis, etaAxis, centralityAxis}}});
      qaregistry.add({"Events/Centrality/NotFoundEventZvtx",
                      "; #it{z}_{vtx} (cm); centrality",
                      {HistType::kTH2F, {zAxis, centralityAxis}}});
      qaregistry.add({"Events/Centrality/ZvtxDiff",
                      "; Z_{rec} - Z_{gen} (cm); centrality",
                      {HistType::kTH2F, {deltaZAxis, centralityAxis}}});
      qaregistry.add(
        {"Events/Centrality/SplitMult",
         "; N_{gen}; #it{z}_{vtx} (cm); centrality",
         {HistType::kTHnSparseF, {multAxis, zAxis, centralityAxis}}});
    }

    if (doprocessTrkEffIdxBestInlusive) {
      qaregistry.add({"Tracks/hPtEtaEffGenFakeBest",
                      "; p_{T} (GeV/c); #eta",
                      {HistType::kTH2F,
                       {ptAxis, etaAxis}}});
      qaregistry.add({"Tracks/hPtEtaEffGenBest",
                      "; p_{T} (GeV/c); #eta",
                      {HistType::kTH2F,
                       {ptAxis, etaAxis}}});
      qaregistry.add({"Tracks/hPtEtaEffPrimBest",
                      "; p_{T} (GeV/c); #eta",
                      {HistType::kTH2F,
                       {ptAxis, etaAxis}}});
      qaregistry.add({"Tracks/hPtEtaEffSecBest",
                      "; p_{T} (GeV/c); #eta",
                      {HistType::kTH2F,
                       {ptAxis, etaAxis}}});
      qaregistry.add({"Tracks/hPtEtaEffGenDuplBest",
                      "; p_{T} (GeV/c); #eta",
                      {HistType::kTH2F,
                       {ptAxis, etaAxis}}});
      qaregistry.add({"Tracks/hPtEtaEffDuplBest",
                      "; p_{T} (GeV/c); #eta",
                      {HistType::kTH2F,
                       {ptAxis, etaAxis}}});
      qaregistry.add({"Tracks/NmftTrkPerPartBest",
                      "; #it{N}_{mft tracks per particle}",
                      {HistType::kTH1F, {multAxis}}});
    }

    if (doprocessTrkEffIdxBestCentFT0C) {
      qaregistry.add(
        {"Tracks/Centrality/hPtEtaEffGenFakeBest",
         "; p_{T} (GeV/c); #eta; centrality",
         {HistType::kTHnSparseF,
          {ptAxis, etaAxis, centralityAxis}}});
      qaregistry.add(
        {"Tracks/Centrality/hPtEtaEffGenBest",
         "; p_{T} (GeV/c); #eta; centrality",
         {HistType::kTHnSparseF,
          {ptAxis, etaAxis, centralityAxis}}});
      qaregistry.add(
        {"Tracks/Centrality/hPtEtaEffPrimBest",
         "; p_{T} (GeV/c); #eta; centrality",
         {HistType::kTHnSparseF,
          {ptAxis, etaAxis, centralityAxis}}});
      qaregistry.add(
        {"Tracks/Centrality/hPtEtaEffSecBest",
         "; p_{T} (GeV/c); #eta; Z_{vtx} (cm)",
         {HistType::kTHnSparseF,
          {ptAxis, etaAxis, centralityAxis}}});
      qaregistry.add(
        {"Tracks/Centrality/hPtEtaEffGenDuplBest",
         "; p_{T} (GeV/c); #eta; centrality",
         {HistType::kTHnSparseF,
          {ptAxis, etaAxis, centralityAxis}}});
      qaregistry.add(
        {"Tracks/Centrality/hPtEtaEffDuplBest",
         "; p_{T} (GeV/c); #eta; centrality",
         {HistType::kTHnSparseF,
          {ptAxis, etaAxis, centralityAxis}}});
      qaregistry.add(
        {"Tracks/Centrality/NmftTrkPerPartBest",
         "; #it{N}_{mft tracks per particle}; centrality",
         {HistType::kTHnSparseF, {multAxis, centralityAxis}}});
    }

    if (doprocessTrkEffIdxInlusive) {
      qaregistry.add({"Tracks/hPtEtaEffGen",
                      "; p_{T} (GeV/c); #eta; occupancy",
                      {HistType::kTHnSparseF,
                       {ptAxis, etaAxis, occupancyAxis}}});
      qaregistry.add({"Tracks/hPtEtaEffPrim",
                      "; p_{T} (GeV/c); #eta; occupancy",
                      {HistType::kTHnSparseF,
                       {ptAxis, etaAxis, occupancyAxis}}});
      qaregistry.add({"Tracks/hPtEtaEffSec",
                      "; p_{T} (GeV/c); #eta; occupancy",
                      {HistType::kTHnSparseF,
                       {ptAxis, etaAxis, occupancyAxis}}});
      qaregistry.add({"Tracks/hPtEtaEffGenDupl",
                      "; p_{T} (GeV/c); #eta; occupancy",
                      {HistType::kTHnSparseF,
                       {ptAxis, etaAxis, occupancyAxis}}});
      qaregistry.add({"Tracks/hPtEtaEffDupl",
                      "; p_{T} (GeV/c); #eta; occupancy",
                      {HistType::kTHnSparseF,
                       {ptAxis, etaAxis, occupancyAxis}}});
      qaregistry.add({"Tracks/NmftTrkPerPart",
                      "; #it{N}_{mft tracks per particle}; occupancy",
                      {HistType::kTH2F, {multAxis, occupancyAxis}}});
    }

    if (doprocessTrkEffIdxCentFT0C) {
      qaregistry.add(
        {"Tracks/Centrality/hPtEtaEffGen",
         "; p_{T} (GeV/c); #eta; centrality; "
         "occupancy",
         {HistType::kTHnSparseF,
          {ptAxis, etaAxis, centralityAxis, occupancyAxis}}});
      qaregistry.add(
        {"Tracks/Centrality/hPtEtaEffPrim",
         "; p_{T} (GeV/c); #eta; centrality; "
         "occupancy",
         {HistType::kTHnSparseF,
          {ptAxis, etaAxis, centralityAxis, occupancyAxis}}});
      qaregistry.add(
        {"Tracks/Centrality/hPtEtaEffSec",
         "; p_{T} (GeV/c); #eta; Z_{vtx} (cm); centrality; "
         "occupancy",
         {HistType::kTHnSparseF,
          {ptAxis, etaAxis, centralityAxis, occupancyAxis}}});
      qaregistry.add(
        {"Tracks/Centrality/hPtEtaEffGenDupl",
         "; p_{T} (GeV/c); #eta; centrality; "
         "occupancy",
         {HistType::kTHnSparseF,
          {ptAxis, etaAxis, centralityAxis, occupancyAxis}}});
      qaregistry.add(
        {"Tracks/Centrality/hPtEtaEffDupl",
         "; p_{T} (GeV/c); #eta; centrality; "
         "occupancy",
         {HistType::kTHnSparseF,
          {ptAxis, etaAxis, centralityAxis, occupancyAxis}}});
      qaregistry.add(
        {"Tracks/Centrality/NmftTrkPerPart",
         "; #it{N}_{mft tracks per particle}; centrality; occupancy",
         {HistType::kTHnSparseF, {multAxis, centralityAxis, occupancyAxis}}});
    }

    if (doprocessTrkEffBestInclusive) {
      qaregistry.add({"Tracks/hPtPhiEtaZvtxEffBestGen",
                      "; p_{T} (GeV/c); #varphi; #eta; Z_{vtx} (cm); occupancy",
                      {HistType::kTHnSparseF,
                       {ptAxis, phiAxis, etaAxis, zAxis, occupancyAxis}}});
      qaregistry.add({"Tracks/hPtPhiEtaZvtxEffBestRec",
                      "; p_{T} (GeV/c); #varphi; #eta; Z_{vtx} (cm); occupancy",
                      {HistType::kTHnSparseF,
                       {ptAxis, phiAxis, etaAxis, zAxis, occupancyAxis}}});
      qaregistry.add({"Tracks/hPtEffBestFakeRec",
                      " ; p_{T} (GeV/c); occupancy",
                      {HistType::kTHnSparseF, {ptAxis, phiAxis, etaAxis, zAxis, occupancyAxis}}});
    }

    if (doprocessTrkEffBestCentFT0C) {
      qaregistry.add(
        {"Tracks/Centrality/hPtPhiEtaZvtxEffBestGen",
         "; p_{T} (GeV/c); #varphi; #eta; Z_{vtx} (cm); centrality; "
         "occupancy",
         {HistType::kTHnSparseF,
          {ptAxis, phiAxis, etaAxis, zAxis, centralityAxis, occupancyAxis}}});
      qaregistry.add(
        {"Tracks/Centrality/hPtPhiEtaZvtxEffBestRec",
         "; p_{T} (GeV/c); #varphi; #eta; Z_{vtx} (cm); centrality; "
         "occupancy",
         {HistType::kTHnSparseF,
          {ptAxis, phiAxis, etaAxis, zAxis, centralityAxis, occupancyAxis}}});
      qaregistry.add(
        {"Tracks/Centrality/hPtEffBestFakeRec",
         "; p_{T} (GeV/c); centrality; occupancy",
         {HistType::kTHnSparseF, {ptAxis, phiAxis, etaAxis, zAxis, centralityAxis, occupancyAxis}}});
    }

    if (doprocessMcQAInclusive) {
      qaregistry.add({"Events/hRecPerGenColls",
                      "; #it{N}_{reco collisions} / #it{N}_{gen collisions};",
                      {HistType::kTH2F, {{200, 0., 2.}, occupancyAxis}}});
      qaregistry.add({"Tracks/hNmftTrks",
                      "; #it{N}_{mft tracks};",
                      {HistType::kTH2F, {{200, -0.5, 200.}, occupancyAxis}}});
      qaregistry.add({"Tracks/hFracAmbiguousMftTrks",
                      "; #it{N}_{ambiguous tracks} / #it{N}_{tracks};",
                      {HistType::kTH2F, {{100, 0., 1.}, occupancyAxis}}});
    }

    if (doprocessMcQACentFT0C) {
      qaregistry.add(
        {"Events/Centrality/hRecPerGenColls",
         "; #it{N}_{reco collisions} / #it{N}_{gen collisions}; centrality",
         {HistType::kTHnSparseF,
          {{200, 0., 2.}, centralityAxis, occupancyAxis}}});
      qaregistry.add({"Tracks/Centrality/hNmftTrks",
                      "; #it{N}_{mft tracks}; centrality",
                      {HistType::kTHnSparseF,
                       {{200, -0.5, 200.}, centralityAxis, occupancyAxis}}});
      qaregistry.add(
        {"Tracks/Centrality/hFracAmbiguousMftTrks",
         "; #it{N}_{ambiguous tracks} / #it{N}_{tracks}; centrality",
         {HistType::kTHnSparseF,
          {{100, 0., 1.}, centralityAxis, occupancyAxis}}});
    }

    if (doprocessCheckAmbiguousMftTracks) {
      qaregistry.add({"Tracks/hMftTracksAmbDegree",
                      " ; N_{coll}^{comp}",
                      {HistType::kTH1F, {{41, -0.5, 40.5}}}});
      qaregistry.add({"Tracks/hMftTracksAmbDegreeWithTrivial",
                      " ; N_{coll}^{comp}",
                      {HistType::kTH1F, {{41, -0.5, 40.5}}}});

      qaregistry.add("Tracks/hAmbTrackType", "hAmbTrackType", {HistType::kTH1F, {{static_cast<int>(AmbTrkType::nAmbTrkType), -0.5, +static_cast<int>(AmbTrkType::nAmbTrkType) - 0.5}}});
      std::string labelAmbiguity[CambTrkType];
      labelAmbiguity[static_cast<int>(AmbTrkType::kAll)] = "all";
      labelAmbiguity[static_cast<int>(AmbTrkType::kOrphan)] = "orphan";
      labelAmbiguity[static_cast<int>(AmbTrkType::kOrphanNull)] = "orphanNull";
      labelAmbiguity[static_cast<int>(AmbTrkType::kNonAmb)] = "nonAmb";
      labelAmbiguity[static_cast<int>(AmbTrkType::kNonAmbSame)] = "nonAmbSame";
      labelAmbiguity[static_cast<int>(AmbTrkType::kAmb)] = "Amb";
      labelAmbiguity[static_cast<int>(AmbTrkType::kAmbGt1)] = "AmbGt1";
      qaregistry.get<TH1>(HIST("Tracks/hAmbTrackType"))->SetMinimum(0.1);
      for (int iBin = 0; iBin < static_cast<int>(AmbTrkType::nAmbTrkType); iBin++) {
        qaregistry.get<TH1>(HIST("Tracks/hAmbTrackType"))->GetXaxis()->SetBinLabel(iBin + 1, labelAmbiguity[iBin].data());
      }
    }

    if (doprocessAssocMC) {
      registry.add("TrackToColl/hAmbTrackType", "hAmbTrackType", {HistType::kTH1F, {{static_cast<int>(AmbTrkType::nAmbTrkType), -0.5, +static_cast<int>(AmbTrkType::nAmbTrkType) - 0.5}}});
      std::string labelAmbiguity[CambTrkType];
      labelAmbiguity[static_cast<int>(AmbTrkType::kAll)] = "all";
      labelAmbiguity[static_cast<int>(AmbTrkType::kOrphan)] = "orphan";
      labelAmbiguity[static_cast<int>(AmbTrkType::kOrphanNull)] = "orphanNull";
      labelAmbiguity[static_cast<int>(AmbTrkType::kNonAmb)] = "nonAmb";
      labelAmbiguity[static_cast<int>(AmbTrkType::kNonAmbSame)] = "nonAmbSame";
      labelAmbiguity[static_cast<int>(AmbTrkType::kAmb)] = "Amb";
      labelAmbiguity[static_cast<int>(AmbTrkType::kAmbGt1)] = "AmbGt1";
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

      registry.add("ReAssocMC/hNReAssocRecoColls", "Number of times generated collisions are reconstructed; N; Counts", HistType::kTH1F, {{10, -0.5, 9.5}});

      registry.add("ReAssocMC/hReAssocMCEventStatus", ";status", {HistType::kTH1F, {{static_cast<int>(ReAssocMCEventStatus::nEvtReAsReAssocMCEventStatus), -0.5, +static_cast<int>(ReAssocMCEventStatus::nEvtReAsReAssocMCEventStatus) - 0.5}}});
      std::string labelReAssocMCEventStatus[CevtReAsReAssocMCEventStatus];
      labelReAssocMCEventStatus[static_cast<int>(ReAssocMCEventStatus::kEvtReAsAll)] = "All";
      labelReAssocMCEventStatus[static_cast<int>(ReAssocMCEventStatus::kEvtReAsSelected)] = "Selected";
      labelReAssocMCEventStatus[static_cast<int>(ReAssocMCEventStatus::kEvtReAsHasMcColl)] = "Has Mc Coll";
      labelReAssocMCEventStatus[static_cast<int>(ReAssocMCEventStatus::kEvtReAsSplitVtxRemoved)] = "Split Vtx Removed";
      labelReAssocMCEventStatus[static_cast<int>(ReAssocMCEventStatus::kEvtReAsZVtxCutMC)] = "Vtx-z cut MC";
      registry.get<TH1>(HIST("ReAssocMC/hReAssocMCEventStatus"))->SetMinimum(0.1);
      for (int iBin = 0; iBin < static_cast<int>(ReAssocMCEventStatus::nEvtReAsReAssocMCEventStatus); iBin++) {
        registry.get<TH1>(HIST("ReAssocMC/hReAssocMCEventStatus"))->GetXaxis()->SetBinLabel(iBin + 1, labelReAssocMCEventStatus[iBin].data());
      }

      registry.add("ReAssocMC/hReAssocMCTrackStatus", ";status", {HistType::kTH1F, {{static_cast<int>(ReAssocMCTrackStatus::nReAssocMCTrackStatus), -0.5, +static_cast<int>(ReAssocMCTrackStatus::nReAssocMCTrackStatus) - 0.5}}});
      std::string labelReAssocMCTrackStatus[CreAssocMCTrackStatus];
      labelReAssocMCTrackStatus[static_cast<int>(ReAssocMCTrackStatus::kTrkReAssocAll)] = "All";
      labelReAssocMCTrackStatus[static_cast<int>(ReAssocMCTrackStatus::kTrkBestSel)] = "Best sel";
      labelReAssocMCTrackStatus[static_cast<int>(ReAssocMCTrackStatus::kTrkSel)] = "Trk sel";
      labelReAssocMCTrackStatus[static_cast<int>(ReAssocMCTrackStatus::kTrkHasColl)] = "Has coll";
      labelReAssocMCTrackStatus[static_cast<int>(ReAssocMCTrackStatus::kTrkReassignedRemoved)] = "Reas rm";
      labelReAssocMCTrackStatus[static_cast<int>(ReAssocMCTrackStatus::kTrkHasMcPart)] = "Has part";
      labelReAssocMCTrackStatus[static_cast<int>(ReAssocMCTrackStatus::kTrkNonAmbAll)] = "Non-amb";
      labelReAssocMCTrackStatus[static_cast<int>(ReAssocMCTrackStatus::kTrkNonAmbGood)] = "Non-amb good coll.";
      labelReAssocMCTrackStatus[static_cast<int>(ReAssocMCTrackStatus::kTrkNonAmbBad)] = "Non-amb bad coll.";
      labelReAssocMCTrackStatus[static_cast<int>(ReAssocMCTrackStatus::kTrkAmbAll)] = "Amb";
      labelReAssocMCTrackStatus[static_cast<int>(ReAssocMCTrackStatus::kTrkAmbGood)] = "Amb good coll.";
      labelReAssocMCTrackStatus[static_cast<int>(ReAssocMCTrackStatus::kTrkAmbBad)] = "Amb bad coll.";
      labelReAssocMCTrackStatus[static_cast<int>(ReAssocMCTrackStatus::kTrkNonAmbAllE)] = "Non-amb (ex)";
      labelReAssocMCTrackStatus[static_cast<int>(ReAssocMCTrackStatus::kTrkNonAmbGoodE)] = "Non-amb good coll. (ex)";
      labelReAssocMCTrackStatus[static_cast<int>(ReAssocMCTrackStatus::kTrkNonAmbBadE)] = "Non-amb bad coll. (ex)";
      labelReAssocMCTrackStatus[static_cast<int>(ReAssocMCTrackStatus::kAssoc)] = "Assoc (gt1 amb)";
      labelReAssocMCTrackStatus[static_cast<int>(ReAssocMCTrackStatus::kAssocGood)] = "Assoc good";
      labelReAssocMCTrackStatus[static_cast<int>(ReAssocMCTrackStatus::kAssocGoodIsCompTrue)] = "Assoc good Comp True";
      labelReAssocMCTrackStatus[static_cast<int>(ReAssocMCTrackStatus::kAssocGoodIsCompFalse)] = "Assoc good Comp False";
      labelReAssocMCTrackStatus[static_cast<int>(ReAssocMCTrackStatus::kAssocBad)] = "Assoc bad";
      labelReAssocMCTrackStatus[static_cast<int>(ReAssocMCTrackStatus::kAssocBadIsCompTrue)] = "Assoc bad Comp True";
      labelReAssocMCTrackStatus[static_cast<int>(ReAssocMCTrackStatus::kAssocBadIsCompFalse)] = "Assoc bad Comp False";
      labelReAssocMCTrackStatus[static_cast<int>(ReAssocMCTrackStatus::kReAssoc)] = "ReAssoc";
      labelReAssocMCTrackStatus[static_cast<int>(ReAssocMCTrackStatus::kReAssocGood)] = "ReAssoc good";
      labelReAssocMCTrackStatus[static_cast<int>(ReAssocMCTrackStatus::kReAssocGoodIsCompTrue)] = "ReAssoc good Comp True";
      labelReAssocMCTrackStatus[static_cast<int>(ReAssocMCTrackStatus::kReAssocGoodIsCompFalse)] = "ReAssoc good Comp False";
      labelReAssocMCTrackStatus[static_cast<int>(ReAssocMCTrackStatus::kReAssocBad)] = "ReAssoc bad";
      labelReAssocMCTrackStatus[static_cast<int>(ReAssocMCTrackStatus::kReAssocBadIsCompTrue)] = "ReAssoc bad Comp True";
      labelReAssocMCTrackStatus[static_cast<int>(ReAssocMCTrackStatus::kReAssocBadIsCompFalse)] = "ReAssoc bad Comp False";
      registry.get<TH1>(HIST("ReAssocMC/hReAssocMCTrackStatus"))->SetMinimum(0.1);
      for (int iBin = 0; iBin < static_cast<int>(ReAssocMCTrackStatus::nReAssocMCTrackStatus); iBin++) {
        registry.get<TH1>(HIST("ReAssocMC/hReAssocMCTrackStatus"))->GetXaxis()->SetBinLabel(iBin + 1, labelReAssocMCTrackStatus[iBin].data());
      }

      // Vertex resolution
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbAll)] = registry.add<THnSparse>("ReAssocMC/hVtxResNonAmbAll", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbGood)] = registry.add<THnSparse>("ReAssocMC/hVtxResNonAmbGood", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbBad)] = registry.add<THnSparse>("ReAssocMC/hVtxResNonAmbBad", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkAmbAll)] = registry.add<THnSparse>("ReAssocMC/hVtxResAmbAll", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkAmbGood)] = registry.add<THnSparse>("ReAssocMC/hVtxResAmbGood", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkAmbBad)] = registry.add<THnSparse>("ReAssocMC/hVtxResAmbBad", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbAllE)] = registry.add<THnSparse>("ReAssocMC/hVtxResNonAmbAllE", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbGoodE)] = registry.add<THnSparse>("ReAssocMC/hVtxResNonAmbGoodE", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbBadE)] = registry.add<THnSparse>("ReAssocMC/hVtxResNonAmbBadE", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kAssoc)] = registry.add<THnSparse>("ReAssocMC/hVtxResAssoc", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kAssocGood)] = registry.add<THnSparse>("ReAssocMC/hVtxResAssocGood", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kAssocGoodIsCompTrue)] = registry.add<THnSparse>("ReAssocMC/hVtxResAssocGoodIsCompTrue", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kAssocGoodIsCompFalse)] = registry.add<THnSparse>("ReAssocMC/hVtxResAssocGoodIsCompFalse", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kAssocBad)] = registry.add<THnSparse>("ReAssocMC/hVtxResAssocBad", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kAssocBadIsCompTrue)] = registry.add<THnSparse>("ReAssocMC/hVtxResAssocBadIsCompTrue", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kAssocBadIsCompFalse)] = registry.add<THnSparse>("ReAssocMC/hVtxResAssocBadIsCompFalse", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kReAssoc)] = registry.add<THnSparse>("ReAssocMC/hVtxResReAssoc", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kReAssocGood)] = registry.add<THnSparse>("ReAssocMC/hVtxResReAssocGood", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kReAssocGoodIsCompTrue)] = registry.add<THnSparse>("ReAssocMC/hVtxResReAssocGoodIsCompTrue", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kReAssocGoodIsCompFalse)] = registry.add<THnSparse>("ReAssocMC/hVtxResReAssocGoodIsCompFalse", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kReAssocBad)] = registry.add<THnSparse>("ReAssocMC/hVtxResReAssocBad", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kReAssocBadIsCompTrue)] = registry.add<THnSparse>("ReAssocMC/hVtxResReAssocBadIsCompTrue", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis});
      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kReAssocBadIsCompFalse)] = registry.add<THnSparse>("ReAssocMC/hVtxResReAssocBadIsCompFalse", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis, deltaZAxis, deltaZAxis});

      // DCA
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbAll)] = registry.add<THnSparse>("ReAssocMC/hDCANonAmbAll", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbGood)] = registry.add<THnSparse>("ReAssocMC/hDCANonAmbGood", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbBad)] = registry.add<THnSparse>("ReAssocMC/hDCANonAmbBad", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkAmbAll)] = registry.add<THnSparse>("ReAssocMC/hDCAAmbAll", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkAmbGood)] = registry.add<THnSparse>("ReAssocMC/hDCAAmbGood", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkAmbBad)] = registry.add<THnSparse>("ReAssocMC/hDCAAmbBad", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbAllE)] = registry.add<THnSparse>("ReAssocMC/hDCANonAmbAllE", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbGoodE)] = registry.add<THnSparse>("ReAssocMC/hDCANonAmbGoodE", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbBadE)] = registry.add<THnSparse>("ReAssocMC/hDCANonAmbBadE", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kAssoc)] = registry.add<THnSparse>("ReAssocMC/hDCAAssoc", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kAssocGood)] = registry.add<THnSparse>("ReAssocMC/hDCAAssocGood", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kAssocGoodIsCompTrue)] = registry.add<THnSparse>("ReAssocMC/hDCAAssocGoodIsCompTrue", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kAssocGoodIsCompFalse)] = registry.add<THnSparse>("ReAssocMC/hDCAAssocGoodIsCompFalse", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kAssocBad)] = registry.add<THnSparse>("ReAssocMC/hDCAAssocBad", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kAssocBadIsCompTrue)] = registry.add<THnSparse>("ReAssocMC/hDCAAssocBadIsCompTrue", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kAssocBadIsCompFalse)] = registry.add<THnSparse>("ReAssocMC/hDCAAssocBadIsCompFalse", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kReAssoc)] = registry.add<THnSparse>("ReAssocMC/hDCAReAssoc", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kReAssocGood)] = registry.add<THnSparse>("ReAssocMC/hDCAReAssocGood", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kReAssocGoodIsCompTrue)] = registry.add<THnSparse>("ReAssocMC/hDCAReAssocGoodIsCompTrue", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kReAssocGoodIsCompFalse)] = registry.add<THnSparse>("ReAssocMC/hDCAReAssocGoodIsCompFalse", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kReAssocBad)] = registry.add<THnSparse>("ReAssocMC/hDCAReAssocBad", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kReAssocBadIsCompTrue)] = registry.add<THnSparse>("ReAssocMC/hDCAReAssocBadIsCompTrue", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis});
      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kReAssocBadIsCompFalse)] = registry.add<THnSparse>("ReAssocMC/hDCAReAssocBadIsCompFalse", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY} (cm)^{reco};  DCA_{Z} (cm)^{reco}; DCA_{XY} (cm);  DCA_{Z} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, dcaxyAxis, dcazAxis, dcaxyAxis, dcazAxis});
    }

    if (doprocessEfficiencyInclusive) {
      qaregistry.add({"Tracks/hEffRec",
                      "; p_{T} (GeV/c); #varphi; #eta; Z_{vtx} (cm); occupancy",
                      {HistType::kTHnSparseF,
                       {ptAxis, phiAxis, etaAxis, zAxis, occupancyAxis}}});
      qaregistry.add({"Tracks/hEffFake",
                      "; p_{T} (GeV/c); #varphi; #eta; Z_{vtx} (cm); occupancy",
                      {HistType::kTHnSparseF,
                       {ptAxis, phiAxis, etaAxis, zAxis, occupancyAxis}}});
    }

    if (doprocessEfficiencyCentFT0C) {
      qaregistry.add({"Tracks/Centrality/hEffRec",
                      "; p_{T} (GeV/c); #varphi; #eta; Z_{vtx} (cm); centrality; occupancy",
                      {HistType::kTHnSparseF,
                       {ptAxis, phiAxis, etaAxis, zAxis, centralityAxis, occupancyAxis}}});
      qaregistry.add({"Tracks/Centrality/hEffFake",
                      "; p_{T} (GeV/c); #varphi; #eta; Z_{vtx} (cm); centrality; occupancy",
                      {HistType::kTHnSparseF,
                       {ptAxis, phiAxis, etaAxis, zAxis, centralityAxis, occupancyAxis}}});
    }

    if (doprocessCorrelationwBestTracksInclusive) {
      qaregistry.add("Events/hMultMFTvsFT0A", "MultMFT_vs_FT0A", {HistType::kTH2F, {multAxis, multFT0aAxis}});
      qaregistry.add("Events/hMultMFTvsFT0C", "MultMFT_vs_FT0C", {HistType::kTH2F, {multAxis, multFT0cAxis}});
      qaregistry.add("Events/hNPVtracksVsFT0C", "NPVtracks_vs_FT0C", {HistType::kTH2F, {pvAxis, multFT0cAxis}});
      qaregistry.add("Events/hMultMFTvsFV0A", "MultMFT_vs_FV0A", {HistType::kTH2F, {multAxis, multFV0aAxis}});
      qaregistry.add("Events/hNPVtracksVsMultMFT", "NPVtracks_vs_MultMFT", {HistType::kTH2F, {pvAxis, multAxis}});
    }

    if (doprocessSecondariesMCInlcusive || doprocessSecondariesMCCentFT0C) {
      auto hNevt = registry.add<TH1>("Events/hNGenRecColls", "Number of generated and reconstructed MC collisions", HistType::kTH1F, {{3, 0.5, 3.5}});
      hNevt->GetXaxis()->SetBinLabel(1, "Reconstructed collisions");
      hNevt->GetXaxis()->SetBinLabel(2, "Generated collisions");

      if (doprocessSecondariesMCInlcusive) {
        qaregistry.add("Events/hTrackToCollEvtType", ";status", {HistType::kTH1F, {{static_cast<int>(TrackToCollEvtType::nTrackToCollEvtType), -0.5, +static_cast<int>(TrackToCollEvtType::nTrackToCollEvtType) - 0.5}}});
        std::string labelTrkToCollEvt[CtrackToCollEvtType];
        labelTrkToCollEvt[static_cast<int>(TrackToCollEvtType::kAllRecColl)] = "all rec";
        labelTrkToCollEvt[static_cast<int>(TrackToCollEvtType::kIsInelGt0wMft)] = "inel>0 mft";
        labelTrkToCollEvt[static_cast<int>(TrackToCollEvtType::kEvtSel)] = "evt sel";
        labelTrkToCollEvt[static_cast<int>(TrackToCollEvtType::kBestCollIdx)] = "bestColl Idx (numContrib)";
        labelTrkToCollEvt[static_cast<int>(TrackToCollEvtType::kIsMcColl)] = "is mc";
        qaregistry.get<TH1>(HIST("Events/hTrackToCollEvtType"))->SetMinimum(0.1);
        for (int iBin = 0; iBin < static_cast<int>(TrackToCollEvtType::nTrackToCollEvtType); iBin++) {
          qaregistry.get<TH1>(HIST("Events/hTrackToCollEvtType"))->GetXaxis()->SetBinLabel(iBin + 1, labelTrkToCollEvt[iBin].data());
        }
        // registry.add({"Events/EvtGenRec", ";status", {HistType::kTH1F, {{3, 0.5, 3.5}}}});
        // auto heff = registry.get<TH1>(HIST("Events/EvtGenRec"));
        // auto* h = heff->GetXaxis();
        // h->SetBinLabel(1, "All generated");
        // h->SetBinLabel(2, "All reconstructed");
        // h->SetBinLabel(3, "Selected reconstructed");

        qaregistry.add("TrkCompColls/hAmbTrackType", ";status", {HistType::kTH1F, {{static_cast<int>(AmbTrkType::nAmbTrkType), -0.5, +static_cast<int>(AmbTrkType::nAmbTrkType) - 0.5}}});
        std::string labelAmbiguity[CambTrkType];
        labelAmbiguity[static_cast<int>(AmbTrkType::kAll)] = "all";
        labelAmbiguity[static_cast<int>(AmbTrkType::kOrphan)] = "orphan";
        labelAmbiguity[static_cast<int>(AmbTrkType::kOrphanNull)] = "orphanNull";
        labelAmbiguity[static_cast<int>(AmbTrkType::kNonAmb)] = "nonAmb";
        labelAmbiguity[static_cast<int>(AmbTrkType::kNonAmbSame)] = "nonAmbSame";
        labelAmbiguity[static_cast<int>(AmbTrkType::kAmb)] = "Amb";
        labelAmbiguity[static_cast<int>(AmbTrkType::kAmbGt1)] = "AmbGt1";
        qaregistry.get<TH1>(HIST("TrkCompColls/hAmbTrackType"))->SetMinimum(0.1);
        for (int iBin = 0; iBin < static_cast<int>(AmbTrkType::nAmbTrkType); iBin++) {
          qaregistry.get<TH1>(HIST("TrkCompColls/hAmbTrackType"))->GetXaxis()->SetBinLabel(iBin + 1, labelAmbiguity[iBin].data());
        }

        registry.add({"Tracks/THnRecAmb", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm)", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis}}});
        registry.add({"Tracks/THnRec", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm)", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis}}});
        registry.add({"Tracks/THnRecNonAmb", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm)", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis}}});
        registry.add({"Tracks/THnRecAmbRest", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm)", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis}}});
        registry.add({"Tracks/THnGenPrim", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm)", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis}}});
        registry.add({"Tracks/THnGenSec", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm)", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis}}});
        registry.add({"Tracks/THnGenSecWeak", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm)", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis}}});
        registry.add({"Tracks/THnGenSecMat", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm)", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis}}});
        registry.add({"Tracks/THnGenPrimAmb", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm)", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis}}});
        registry.add({"Tracks/THnGenSecAmb", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm)", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis}}});
        registry.add({"Tracks/THnGenSecWeakAmb", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm)", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis}}});
        registry.add({"Tracks/THnGenSecMatAmb", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm)", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis}}});
      }
      if (doprocessSecondariesMCCentFT0C) {
        qaregistry.add("Events/Centrality/hTrackToCollEvtType", ";status;centrality", {HistType::kTH2F, {{static_cast<int>(TrackToCollEvtType::nTrackToCollEvtType), -0.5, +static_cast<int>(TrackToCollEvtType::nTrackToCollEvtType) - 0.5}, centralityAxis}});
        std::string labelTrkToCollEvt[CtrackToCollEvtType];
        labelTrkToCollEvt[static_cast<int>(TrackToCollEvtType::kAllRecColl)] = "all rec";
        labelTrkToCollEvt[static_cast<int>(TrackToCollEvtType::kIsInelGt0wMft)] = "inel>0 mft";
        labelTrkToCollEvt[static_cast<int>(TrackToCollEvtType::kEvtSel)] = "evt sel";
        labelTrkToCollEvt[static_cast<int>(TrackToCollEvtType::kBestCollIdx)] = "bestColl Idx (numContrib)";
        labelTrkToCollEvt[static_cast<int>(TrackToCollEvtType::kIsMcColl)] = "is mc";
        qaregistry.get<TH1>(HIST("Events/Centrality/hTrackToCollEvtType"))->SetMinimum(0.1);
        for (int iBin = 0; iBin < static_cast<int>(TrackToCollEvtType::nTrackToCollEvtType); iBin++) {
          qaregistry.get<TH2>(HIST("Events/Centrality/hTrackToCollEvtType"))->GetXaxis()->SetBinLabel(iBin + 1, labelTrkToCollEvt[iBin].data());
        }
        // registry.add({"Events/Centrality/EvtGenRec", ";status;centrality", {HistType::kTH2F, {{3, 0.5, 3.5}, centralityAxis}}});
        // auto heff = registry.get<TH2>(HIST("Events/Centrality/EvtGenRec"));
        // auto* h = heff->GetXaxis();
        // h->SetBinLabel(1, "All generated");
        // h->SetBinLabel(2, "All reconstructed");
        // h->SetBinLabel(3, "Selected reconstructed");

        qaregistry.add("TrkCompColls/Centrality/hAmbTrackType", ";status;centrality", {HistType::kTH2F, {{static_cast<int>(AmbTrkType::nAmbTrkType), -0.5, +static_cast<int>(AmbTrkType::nAmbTrkType) - 0.5}, centralityAxis}});
        std::string labelAmbiguity[CambTrkType];
        labelAmbiguity[static_cast<int>(AmbTrkType::kAll)] = "all";
        labelAmbiguity[static_cast<int>(AmbTrkType::kOrphan)] = "orphan";
        labelAmbiguity[static_cast<int>(AmbTrkType::kOrphanNull)] = "orphanNull";
        labelAmbiguity[static_cast<int>(AmbTrkType::kNonAmb)] = "nonAmb";
        labelAmbiguity[static_cast<int>(AmbTrkType::kNonAmbSame)] = "nonAmbSame";
        labelAmbiguity[static_cast<int>(AmbTrkType::kAmb)] = "Amb";
        labelAmbiguity[static_cast<int>(AmbTrkType::kAmbGt1)] = "AmbGt1";
        qaregistry.get<TH2>(HIST("TrkCompColls/Centrality/hAmbTrackType"))->SetMinimum(0.1);
        for (int iBin = 0; iBin < static_cast<int>(AmbTrkType::nAmbTrkType); iBin++) {
          qaregistry.get<TH2>(HIST("TrkCompColls/Centrality/hAmbTrackType"))->GetXaxis()->SetBinLabel(iBin + 1, labelAmbiguity[iBin].data());
        }

        registry.add({"Tracks/Centrality/THnRecAmb", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis, centralityAxis}}});
        registry.add({"Tracks/Centrality/THnRec", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis, centralityAxis}}});
        registry.add({"Tracks/Centrality/THnRecNonAmb", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis, centralityAxis}}});
        registry.add({"Tracks/Centrality/THnRecAmbRest", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis, centralityAxis}}});
        registry.add({"Tracks/Centrality/THnGenPrim", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis, centralityAxis}}});
        registry.add({"Tracks/Centrality/THnGenSec", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis, centralityAxis}}});
        registry.add({"Tracks/Centrality/THnGenSecWeak", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis, centralityAxis}}});
        registry.add({"Tracks/Centrality/THnGenSecMat", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis, centralityAxis}}});
        registry.add({"Tracks/Centrality/THnGenPrimAmb", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis, centralityAxis}}});
        registry.add({"Tracks/Centrality/THnGenSecAmb", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis, centralityAxis}}});
        registry.add({"Tracks/Centrality/THnGenSecWeakAmb", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis, centralityAxis}}});
        registry.add({"Tracks/Centrality/THnGenSecMatAmb", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis, centralityAxis}}});
      }
    }

    if (doprocessAlignmentInclusive) {
      for (size_t j = 0; j < mftQuadrant.size(); j++) {
        const auto& quadrant = mftQuadrant[j];
        std::string hPath = std::string("Alignment/") + quadrant + "/";
        hAlignment[0][j][0]["DCA_x_vs_z"] = registryAlign.add((hPath + "DCA_x_vs_z").c_str(), std::format("DCA(x) vs. z - {}", quadrant).c_str(), {HistType::kTH2F, {zAxis, dcaxyAxis}});
        hAlignment[0][j][0]["DCA_y_vs_z"] = registryAlign.add((hPath + "DCA_y_vs_z").c_str(), std::format("DCA(y) vs. z - {}", quadrant).c_str(), {HistType::kTH2F, {zAxis, dcaxyAxis}});
      }
    }

    if (doprocessReassocEfficiency) {
      registry.add("Events/hNReEffColls", "Number of times generated collisions are reconstructed; N; Counts", HistType::kTH1F, {{10, -0.5, 9.5}});
      registry.add("Events/hNchTVX", "; status;", HistType::kTH2F, {multAxis, {2, 0, 2}});

      registry.add("Events/Centrality/ReEffStatus", ";status;centrality", HistType::kTH2F, {{9, 0.5, 9.5}, centralityAxis});
      auto hstat = registry.get<TH2>(HIST("Events/Centrality/ReEffStatus"));
      hstat->GetXaxis()->SetBinLabel(1, "All tracks");
      hstat->GetXaxis()->SetBinLabel(2, "Ambiguous tracks");
      hstat->GetXaxis()->SetBinLabel(3, "Reassigned tracks");
      hstat->GetXaxis()->SetBinLabel(4, "Extra tracks");
      hstat->GetXaxis()->SetBinLabel(5, "Originally correctly reassgined");
      hstat->GetXaxis()->SetBinLabel(6, "Correctly reassigned");
      hstat->GetXaxis()->SetBinLabel(7, "Not reassigned (reassigned)");
      hstat->GetXaxis()->SetBinLabel(8, "Not reassigned");
      hstat->GetXaxis()->SetBinLabel(9, "Correctly assigned true");

      registry.add({"AmbTracks/hVtxzMCrec", " ; Z_{vtx} (cm)", {HistType::kTH1F, {zAxis}}});
      registry.add({"AmbTracks/DCAXY", " ; DCA_{XY} (cm)", {HistType::kTH1F, {dcaxyAxis}}});
      registry.add({"AmbTracks/DCAZ", " ; DCA_{Z} (cm)", {HistType::kTH1F, {dcazAxis}}});
      registry.add({"AmbTracks/DCAXYBest", " ; DCA_{XY} (cm)", {HistType::kTH1F, {dcaxyAxis}}});
      registry.add({"AmbTracks/DCAZBest", " ; DCA_{Z} (cm)", {HistType::kTH1F, {dcazAxis}}});
      registry.add({"AmbTracks/DCAXYBestPrim", " ; DCA_{XY} (cm)", {HistType::kTH1F, {dcaxyAxis}}});
      registry.add({"AmbTracks/DCAZBestPrim", " ; DCA_{Z} (cm)", {HistType::kTH1F, {dcazAxis}}});
      registry.add({"AmbTracks/DCAXYBestTrue", " ; DCA_{XY} (cm)", {HistType::kTH1F, {dcaxyAxis}}});
      registry.add({"AmbTracks/DCAZBestTrue", " ; DCA_{Z} (cm)", {HistType::kTH1F, {dcazAxis}}});
      registry.add({"AmbTracks/DCAXYBestFalse", " ; DCA_{XY} (cm)", {HistType::kTH1F, {dcaxyAxis}}});
      registry.add({"AmbTracks/DCAZBestFalse", " ; DCA_{Z} (cm)", {HistType::kTH1F, {dcazAxis}}});
      registry.add({"AmbTracks/DCAXYBestTrueOrigAssoc", " ; DCA_{XY} (cm)", {HistType::kTH1F, {dcaxyAxis}}});
      registry.add({"AmbTracks/DCAXYBestTrueOrigReAssoc", " ; DCA_{XY} (cm)", {HistType::kTH1F, {dcaxyAxis}}});
      registry.add({"AmbTracks/DCAZBestTrueOrigAssoc", " ; DCA_{Z} (cm)", {HistType::kTH1F, {dcazAxis}}});
      registry.add({"AmbTracks/DCAZBestTrueOrigReAssoc", " ; DCA_{Z} (cm)", {HistType::kTH1F, {dcazAxis}}});

      registry.add({"AmbTracks/Centrality/THnDCAxyBestGen", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm); DCA_{XY} (cm);  DCA_{Z} (cm)", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis, dcaxyAxis, dcazAxis, centralityAxis}}});
      registry.add({"AmbTracks/Centrality/THnDCAxyBestGenPrim", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm); DCA_{XY} (cm);  DCA_{Z} (cm)", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis, dcaxyAxis, dcazAxis, centralityAxis}}});
      registry.add({"AmbTracks/Centrality/THnDCAxyBestTrue", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm); DCA_{XY} (cm);  DCA_{Z} (cm)", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis, dcaxyAxis, dcazAxis, centralityAxis}}});
      registry.add({"AmbTracks/Centrality/THnDCAxyBestFalse", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm); DCA_{XY} (cm);  DCA_{Z} (cm)", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis, dcaxyAxis, dcazAxis, centralityAxis}}});
      registry.add({"AmbTracks/Centrality/THnDCAxyBestTrueOrigAssoc", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm); DCA_{XY} (cm);  DCA_{Z} (cm)", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis, dcaxyAxis, dcazAxis, centralityAxis}}});
      registry.add({"AmbTracks/Centrality/THnDCAxyBestTrueOrigReAssoc", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm); DCA_{XY} (cm);  DCA_{Z} (cm)", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis, dcaxyAxis, dcazAxis, centralityAxis}}});

      registry.add({"AmbTracks/BestGenDxyz", "; #Delta X (cm); #Delta Y (cm); #Delta Z (cm); DCA_{XY} (cm);  DCA_{Z} (cm)", {HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis, dcaxyAxis, dcazAxis}}});
      registry.add({"AmbTracks/BestPrimDxyz", "; #Delta X (cm); #Delta Y (cm); #Delta Z (cm); DCA_{XY} (cm);  DCA_{Z} (cm)", {HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis, dcaxyAxis, dcazAxis}}});
      registry.add({"AmbTracks/BestTrueDxyz", "; #Delta X (cm); #Delta Y (cm); #Delta Z (cm); DCA_{XY} (cm);  DCA_{Z} (cm)", {HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis, dcaxyAxis, dcazAxis}}});
      registry.add({"AmbTracks/BestFalseDxyz", "; #Delta X (cm); #Delta Y (cm); #Delta Z (cm); DCA_{XY} (cm);  DCA_{Z} (cm)", {HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis, dcaxyAxis, dcazAxis}}});
      registry.add({"AmbTracks/BestTrueOrigReAssocDxyz", "; #Delta X (cm); #Delta Y (cm); #Delta Z (cm); DCA_{XY} (cm);  DCA_{Z} (cm)", {HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis, dcaxyAxis, dcazAxis}}});
      registry.add({"AmbTracks/BestTrueOrigAssocDxyz", "; #Delta X (cm); #Delta Y (cm); #Delta Z (cm); DCA_{XY} (cm);  DCA_{Z} (cm)", {HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis, dcaxyAxis, dcazAxis}}});
    }

    if (doprocessEventAndSignalLossCentFT0C) {
      auto hNevt = registry.add<TH1>("Events/Centrality/hNRecCollsSigEvtLoss", "Number of reconstructed MC collisions", HistType::kTH1F, {{1, 0.5, 1.5}});
      hNevt->GetXaxis()->SetBinLabel(1, "Reconstructed collisions");

      registry.add("Events/Centrality/EvtSigLossStatus", ";status;centrality", {HistType::kTH2F, {{3, 0.5, 3.5}, centralityAxis}});
      auto hstat = registry.get<TH2>(HIST("Events/Centrality/EvtSigLossStatus"));
      hstat->GetXaxis()->SetBinLabel(1, "All MC gen events");
      hstat->GetXaxis()->SetBinLabel(2, "MC gen events with rec event with event selection");
      hstat->GetXaxis()->SetBinLabel(3, "MC gen events with no rec events");

      registry.add("Events/Centrality/hMultGenVsCent", "event mult MC gen", {HistType::kTH2F, {centralityAxis, multFT0cAxis}});
      registry.add("Events/Centrality/hMultGenVsCentNParticlesEta05", "event mult MC gen", {HistType::kTH2F, {centralityAxis, multAxis}});
      registry.add("Events/Centrality/hMultGenVsCentNParticlesEtaMFT", "event mult MC gen", {HistType::kTH2F, {centralityAxis, multAxis}});
      registry.add("Events/Centrality/hMultGenVsCentRec", "event mult MC gen vs centrality", {HistType::kTH2F, {centralityAxis, multFT0cAxis}});
      registry.add("Events/Centrality/hMultGenVsCentRecNParticlesEta05", "event mult MC gen vs centrality", {HistType::kTH2F, {centralityAxis, multAxis}});
      registry.add("Events/Centrality/hMultGenVsCentRecNParticlesEtaMFT", "event mult MC gen vs centrality", {HistType::kTH2F, {centralityAxis, multAxis}});

      registry.add({"Tracks/Centrality/EtaCentVsMultGen_t", "; #eta; centrality; mult gen", {HistType::kTHnSparseF, {etaAxis, centralityAxis, multFT0cAxis}}});
      registry.add({"Tracks/Centrality/EtaCentVsMultGen", "; #eta; centrality; mult gen", {HistType::kTHnSparseF, {etaAxis, centralityAxis, multFT0cAxis}}});

      registry.add({"Tracks/Centrality/EtaGen_t", "; #eta; centrality; occupancy", {HistType::kTH2F, {etaAxis, centralityAxis}}});
      registry.add({"Tracks/Centrality/EtaGen", "; #eta; centrality; occupancy", {HistType::kTH2F, {etaAxis, centralityAxis}}});
    }

    if (doprocessTimeAssocMC) {

      registry.add("TimeAssocMC/hTimeAssocMCEventStatus", ";status", {HistType::kTH1F, {{static_cast<int>(ReAssocMCEventStatus::nEvtReAsReAssocMCEventStatus), -0.5, +static_cast<int>(ReAssocMCEventStatus::nEvtReAsReAssocMCEventStatus) - 0.5}}});
      std::string labelTimeAssocMCEventStatus[CevtReAsReAssocMCEventStatus];
      labelTimeAssocMCEventStatus[static_cast<int>(ReAssocMCEventStatus::kEvtReAsAll)] = "All";
      labelTimeAssocMCEventStatus[static_cast<int>(ReAssocMCEventStatus::kEvtReAsSelected)] = "Selected";
      labelTimeAssocMCEventStatus[static_cast<int>(ReAssocMCEventStatus::kEvtReAsHasMcColl)] = "Has Mc Coll";
      labelTimeAssocMCEventStatus[static_cast<int>(ReAssocMCEventStatus::kEvtReAsSplitVtxRemoved)] = "Split Vtx Removed";
      labelTimeAssocMCEventStatus[static_cast<int>(ReAssocMCEventStatus::kEvtReAsZVtxCutMC)] = "Vtx-z cut MC";
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
      registry.add({"TimeAssocMC/hVTXkSelGoodVtxBad", "; #Delta X (cm); #Delta Y (cm); #Delta Z (cm)", {HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis}}});
      registry.add({"TimeAssocMC/hVTXkSelNonAmbAll", "; #Delta X (cm); #Delta Y (cm); #Delta Z (cm)", {HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis}}});
      registry.add({"TimeAssocMC/hVTXkSelNonAmbGoodVtxTrue", "; #Delta X (cm); #Delta Y (cm); #Delta Z (cm)", {HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis}}});
      registry.add({"TimeAssocMC/hVTXkSelNonAmbGoodVtxBad", "; #Delta X (cm); #Delta Y (cm); #Delta Z (cm)", {HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis}}});
      registry.add({"TimeAssocMC/hVTXkSelNonAmbSameAll", "; #Delta X (cm); #Delta Y (cm); #Delta Z (cm)", {HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis}}});
      registry.add({"TimeAssocMC/hVTXkSelNonAmbSameGoodVtxTrue", "; #Delta X (cm); #Delta Y (cm); #Delta Z (cm)", {HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis}}});
      registry.add({"TimeAssocMC/hVTXkSelNonAmbSameGoodVtxBad", "; #Delta X (cm); #Delta Y (cm); #Delta Z (cm)", {HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis}}});
      registry.add({"TimeAssocMC/hVTXkSelAmbAll", "; #Delta X (cm); #Delta Y (cm); #Delta Z (cm)", {HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis}}});
      registry.add({"TimeAssocMC/hVTXkSelAmbGoodVtxTrue", "; #Delta X (cm); #Delta Y (cm); #Delta Z (cm)", {HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis}}});
      registry.add({"TimeAssocMC/hVTXkSelAmbGoodVtxBad", "; #Delta X (cm); #Delta Y (cm); #Delta Z (cm)", {HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis}}});
      registry.add({"TimeAssocMC/hVTXkSelAmbGt1All", "; #Delta X (cm); #Delta Y (cm); #Delta Z (cm)", {HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis}}});
      registry.add({"TimeAssocMC/hVTXkSelAmbGt1GoodVtxTrue", "; #Delta X (cm); #Delta Y (cm); #Delta Z (cm)", {HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis}}});
      registry.add({"TimeAssocMC/hVTXkSelAmbGt1GoodVtxBad", "; #Delta X (cm); #Delta Y (cm); #Delta Z (cm)", {HistType::kTHnSparseF, {deltaZAxis, deltaZAxis, deltaZAxis}}});

      registry.add("TimeAssocMC/hTimeAssocCheckVtxType", ";status", {HistType::kTH1F, {{static_cast<int>(ReassocCheckVtxType::nReassocVtxType), -0.5, +static_cast<int>(ReassocCheckVtxType::nReassocVtxType) - 0.5}}});
      std::string labelReAssocVtxType[CreassocVtxType];
      labelReAssocVtxType[static_cast<int>(ReassocCheckVtxType::kIsTrueVtxAllTrue)] = "kIsTrueVtxAll=True";
      labelReAssocVtxType[static_cast<int>(ReassocCheckVtxType::kIsTrueVtxAllFalse)] = "kIsTrueVtxAll=False";
      labelReAssocVtxType[static_cast<int>(ReassocCheckVtxType::kIsTrueVtxVsGoodVtxTrue)] = "kIsTrueVtxVsGoodVtx=True";
      labelReAssocVtxType[static_cast<int>(ReassocCheckVtxType::kIsTrueVtxVsGoodVtxFalse)] = "kIsTrueVtxVsGoodVtx=False";
      labelReAssocVtxType[static_cast<int>(ReassocCheckVtxType::kIsTrueVtxVsBadVtxTrue)] = "kIsTrueVtxVsBadVtx=True";
      labelReAssocVtxType[static_cast<int>(ReassocCheckVtxType::kIsTrueVtxVsBadVtxFalse)] = "kIsTrueVtxVsBadVtx=False";
      registry.get<TH1>(HIST("TimeAssocMC/hTimeAssocCheckVtxType"))->SetMinimum(0.1);
      for (int iBin = 0; iBin < static_cast<int>(ReassocCheckVtxType::nReassocVtxType); iBin++) {
        registry.get<TH1>(HIST("TimeAssocMC/hTimeAssocCheckVtxType"))->GetXaxis()->SetBinLabel(iBin + 1, labelReAssocVtxType[iBin].data());
      }

      registry.add("TimeAssocMC/hAmbTrackType", ";status", {HistType::kTH1F, {{static_cast<int>(AmbTrkType::nAmbTrkType), -0.5, +static_cast<int>(AmbTrkType::nAmbTrkType) - 0.5}}});
      std::string labelAmbiguity[CambTrkType];
      labelAmbiguity[static_cast<int>(AmbTrkType::kAll)] = "all";
      labelAmbiguity[static_cast<int>(AmbTrkType::kOrphan)] = "orphan";
      labelAmbiguity[static_cast<int>(AmbTrkType::kOrphanNull)] = "orphanNull";
      labelAmbiguity[static_cast<int>(AmbTrkType::kNonAmb)] = "nonAmb";
      labelAmbiguity[static_cast<int>(AmbTrkType::kNonAmbSame)] = "nonAmbSame";
      labelAmbiguity[static_cast<int>(AmbTrkType::kAmb)] = "Amb";
      labelAmbiguity[static_cast<int>(AmbTrkType::kAmbGt1)] = "AmbGt1";
      registry.get<TH1>(HIST("TimeAssocMC/hAmbTrackType"))->SetMinimum(0.1);
      for (int iBin = 0; iBin < static_cast<int>(AmbTrkType::nAmbTrkType); iBin++) {
        registry.get<TH1>(HIST("TimeAssocMC/hAmbTrackType"))->GetXaxis()->SetBinLabel(iBin + 1, labelAmbiguity[iBin].data());
      }

      registry.add("TimeAssocMC/hAmbTrkTypeAssocFlag", ";status", {HistType::kTH1F, {{static_cast<int>(AmbTrkTypeAssocFlag::nSelAmbTrkTypeAssocFlag), -0.5, +static_cast<int>(AmbTrkTypeAssocFlag::nSelAmbTrkTypeAssocFlag) - 0.5}}});
      std::string lAmbTrackType[CselAmbTrkTypeAssocFlag];
      lAmbTrackType[static_cast<int>(AmbTrkTypeAssocFlag::kSel)] = "all sel";
      lAmbTrackType[static_cast<int>(AmbTrkTypeAssocFlag::kSelGoodVtxTrue)] = "all good vtx";
      lAmbTrackType[static_cast<int>(AmbTrkTypeAssocFlag::kSelGoodVtxBad)] = "all bad vtx";
      lAmbTrackType[static_cast<int>(AmbTrkTypeAssocFlag::kSelNonAmbAll)] = "non-amb";
      lAmbTrackType[static_cast<int>(AmbTrkTypeAssocFlag::kSelNonAmbGoodVtxTrue)] = "non-amb good vtx";
      lAmbTrackType[static_cast<int>(AmbTrkTypeAssocFlag::kSelNonAmbGoodVtxBad)] = "non-amb bad vtx";
      lAmbTrackType[static_cast<int>(AmbTrkTypeAssocFlag::kSelNonAmbSameAll)] = "non-amb (same)";
      lAmbTrackType[static_cast<int>(AmbTrkTypeAssocFlag::kSelNonAmbSameGoodVtxTrue)] = "non-amb (same) good vtx";
      lAmbTrackType[static_cast<int>(AmbTrkTypeAssocFlag::kSelNonAmbSameGoodVtxBad)] = "non-amb (same) bad vtx";
      lAmbTrackType[static_cast<int>(AmbTrkTypeAssocFlag::kSelAmbAll)] = "amb";
      lAmbTrackType[static_cast<int>(AmbTrkTypeAssocFlag::kSelAmbGoodVtxTrue)] = "amb good vtx";
      lAmbTrackType[static_cast<int>(AmbTrkTypeAssocFlag::kSelAmbGoodVtxBad)] = "amb bad vtx";
      lAmbTrackType[static_cast<int>(AmbTrkTypeAssocFlag::kSelAmbGt1All)] = "ambGt1";
      lAmbTrackType[static_cast<int>(AmbTrkTypeAssocFlag::kSelAmbGt1GoodVtxTrue)] = "ambGt1 good vtx";
      lAmbTrackType[static_cast<int>(AmbTrkTypeAssocFlag::kSelAmbGt1GoodVtxBad)] = "ambGt1 bad vtx";
      lAmbTrackType[static_cast<int>(AmbTrkTypeAssocFlag::kSelOrphanNull)] = "orhpan null";
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
      std::string labelReAssocVtxType[CreassocVtxType];
      labelReAssocVtxType[static_cast<int>(ReassocCheckVtxType::kIsTrueVtxAllTrue)] = "kIsTrueVtxAll=True";
      labelReAssocVtxType[static_cast<int>(ReassocCheckVtxType::kIsTrueVtxAllFalse)] = "kIsTrueVtxAll=False";
      labelReAssocVtxType[static_cast<int>(ReassocCheckVtxType::kIsRecGoodAllTrue)] = "IsRecGoodAll=True";
      labelReAssocVtxType[static_cast<int>(ReassocCheckVtxType::kIsRecGoodAllFalse)] = "kIsRecGoodAll=False";
      labelReAssocVtxType[static_cast<int>(ReassocCheckVtxType::kIsRecGoodMatchAllTrue)] = "kIsRecGoodMatchAll=True";
      labelReAssocVtxType[static_cast<int>(ReassocCheckVtxType::kIsRecGoodMatchAllFalse)] = "kIsRecGoodMatchAll=False";
      //
      labelReAssocVtxType[static_cast<int>(ReassocCheckVtxType::kIsTrueVtxVsGoodVtxTrue)] = "kIsTrueVtxVsGoodVtx=True";
      labelReAssocVtxType[static_cast<int>(ReassocCheckVtxType::kIsTrueVtxVsGoodVtxFalse)] = "kIsTrueVtxVsGoodVtx=False";
      labelReAssocVtxType[static_cast<int>(ReassocCheckVtxType::kIsRecGoodVsGoodVtxTrue)] = "kIsRecGoodVsGoodVtx=True";
      labelReAssocVtxType[static_cast<int>(ReassocCheckVtxType::kIsRecGoodVsGoodVtxFalse)] = "kIsRecGoodVsGoodVtx=False";
      labelReAssocVtxType[static_cast<int>(ReassocCheckVtxType::kIsRecGoodMatchVsGoodVtxTrue)] = "kIsRecGoodMatchVsGoodVtx=True";
      labelReAssocVtxType[static_cast<int>(ReassocCheckVtxType::kIsRecGoodMatchVsGoodVtxFalse)] = "kIsRecGoodMatchVsGoodVtx=False";
      //
      labelReAssocVtxType[static_cast<int>(ReassocCheckVtxType::kIsTrueVtxVsBadVtxTrue)] = "kIsTrueVtxVsBadVtx=True";
      labelReAssocVtxType[static_cast<int>(ReassocCheckVtxType::kIsTrueVtxVsBadVtxFalse)] = "kIsTrueVtxVsBadVtx=False";
      labelReAssocVtxType[static_cast<int>(ReassocCheckVtxType::kIsRecGoodVsBadVtxTrue)] = "kIsRecGoodVsBadVtx=True";
      labelReAssocVtxType[static_cast<int>(ReassocCheckVtxType::kIsRecGoodVsBadVtxFalse)] = "kIsRecGoodVsBadVtx=False";
      labelReAssocVtxType[static_cast<int>(ReassocCheckVtxType::kIsRecGoodMatchVsBadVtxTrue)] = "kIsRecGoodMatchVsBadVtx=True";
      labelReAssocVtxType[static_cast<int>(ReassocCheckVtxType::kIsRecGoodMatchVsBadVtxFalse)] = "kIsRecGoodMatchVsBadVtx=False";
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
  }

  /// Filters - tracks
  Filter filtTrkEta = (aod::fwdtrack::eta < trackCuts.maxEta) &&
                      (aod::fwdtrack::eta > trackCuts.minEta);
  Filter filtATrackID = (aod::fwdtrack::bestCollisionId >= CintZero);
  Filter filtATrackDCAxy = (nabs(aod::fwdtrack::bestDCAXY) < trackCuts.maxDCAxy);
  Filter filtATrackDCAz = (nabs(aod::fwdtrack::bestDCAZ) < trackCuts.maxDCAz);

  /// Filters - mc particles
  Filter primaries = (aod::mcparticle::flags & (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary) == (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary && (aod::mcparticle::eta < trackCuts.maxEta) && (aod::mcparticle::eta > trackCuts.minEta);

  /// Joined tables
  using FullBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;
  using CollBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
  using ExtBCs = soa::Join<aod::BCs, aod::Timestamps, aod::MatchedBCCollisionsSparseMulti>;

  /// Collisions
  using Colls = soa::Join<aod::Collisions, aod::EvSels>;
  using Coll = Colls::iterator;
  using CollsCentFT0C = soa::Join<aod::Collisions, aod::CentFT0Cs, aod::EvSels>;
  using CollsCentFT0CVariant1 =
    soa::Join<aod::Collisions, aod::CentFT0CVariant1s, aod::EvSels>;
  using CollsCentFT0M = soa::Join<aod::Collisions, aod::CentFT0Ms, aod::EvSels>;
  using CollsCentNGlobal =
    soa::Join<aod::Collisions, aod::CentNGlobals, aod::EvSels>;
  using CollsCentMFT = soa::Join<aod::Collisions, aod::CentMFTs, aod::EvSels>;
  using CollCentFT0C = CollsCentFT0C::iterator;
  using CollsGenCentFT0C = soa::Join<aod::McCollisionLabels, aod::Collisions,
                                     aod::CentFT0Cs, aod::EvSels>;
  using CollisionsWithMCLabels = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
  using CollGenCent = CollsGenCentFT0C::iterator;
  using CollsCorr = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::PVMults, aod::CentFT0Cs, aod::CentFV0As, aod::CentFT0CVariant1s, aod::CentFT0Ms, aod::CentNGlobals, aod::CentMFTs>;
  using CollsMCExtra = soa::Join<aod::McCollisions, aod::McCollsExtra>;
  using CollsMCExtraMult = soa::Join<aod::McCollisions, aod::MultMCExtras, aod::McCollsExtra>;

  /// Tracks
  using MFTTracksLabeled = soa::Join<aod::MFTTracks, aod::McMFTTrackLabels>;
  using MftTracksWColls = soa::Join<aod::MFTTracks, aod::MFTTrkCompColls>;
  using MftTracksWCollsMC = soa::Join<aod::MFTTracks, aod::MFTTrkCompColls, aod::McMFTTrackLabels>;
  using BestTracksMC = soa::Join<aod::MFTTracks, aod::BestCollisionsFwd3d, aod::McMFTTrackLabels>;
  using BestTracks3dWCollsMC = soa::Join<aod::MFTTracks, aod::MFTTrkCompColls, aod::BestCollisionsFwd3d, aod::McMFTTrackLabels>;
  using BestTracks2dWCollsMC = soa::Join<aod::MFTTracks, aod::MFTTrkCompColls, aod::BestCollisionsFwd, aod::McMFTTrackLabels>;

  /// Filtered tables
  using FiltMftTracks = soa::Filtered<aod::MFTTracks>;
  using FiltMcMftTracks = soa::Filtered<MFTTracksLabeled>;
  using FiltBestTracks = soa::Filtered<aod::BestCollisionsFwd3d>;
  using FiltMcBestTracks = soa::Filtered<BestTracksMC>;

  using FiltParticles = soa::Filtered<aod::McParticles>;

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

    o2::field::MagneticField* field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
    bZ = field->getBz(CcenterMFT);
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
    if (std::abs(besttrack.bestDCAZ()) >= trackCuts.maxDCAz) {
      return false;
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
      if (mftChi2NCl > trackCuts.maxChi2NCl)
        return false;
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
           (phi < ((PIHalf - 0.1) * PI) + trackCuts.phiCut)))
        return false;
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

  template <typename C, bool fillHis = false, typename T>
  int countTracks(T const& tracks, float z, float c, float occ)
  {
    auto nTrk = 0;
    for (auto const& track : tracks) {
      if (fillHis) {
        if constexpr (has_reco_cent<C>) {
          qaregistry.fill(HIST("Tracks/Centrality/Chi2Eta"), track.chi2(), track.eta(), c, occ);
          qaregistry.fill(HIST("Tracks/Centrality/Chi2"), track.chi2(), c, occ);
          qaregistry.fill(HIST("Tracks/Centrality/NclustersEta"), track.nClusters(), track.eta(), c, occ);
        } else {
          qaregistry.fill(HIST("Tracks/Chi2Eta"), track.chi2(), track.eta(), occ);
          qaregistry.fill(HIST("Tracks/Chi2"), track.chi2(), occ);
          qaregistry.fill(HIST("Tracks/NclustersEta"), track.nClusters(), track.eta(), occ);
        }
      }
      if (!isTrackSelected(track)) {
        continue;
      }
      if (fillHis) {
        float phi = track.phi();
        o2::math_utils::bringTo02Pi(phi);
        if (phi < Czero || TwoPI < phi) {
          continue;
        }
        if constexpr (has_reco_cent<C>) {
          registry.fill(HIST("Tracks/Centrality/EtaZvtx"), track.eta(), z, c, occ);
          registry.fill(HIST("Tracks/Centrality/PhiEta"), phi, track.eta(), c, occ);
          qaregistry.fill(HIST("Tracks/Centrality/TanLambda"), track.tgl(), c, occ);
          qaregistry.fill(HIST("Tracks/Centrality/InvQPt"), track.signed1Pt(), c, occ);
          qaregistry.fill(HIST("Tracks/Centrality/Eta"), track.eta(), c, occ);
          qaregistry.fill(HIST("Tracks/Centrality/Phi"), phi, c, occ);
        } else {
          registry.fill(HIST("Tracks/EtaZvtx"), track.eta(), z, occ);
          registry.fill(HIST("Tracks/PhiEta"), phi, track.eta(), occ);
          qaregistry.fill(HIST("Tracks/TanLambda"), track.tgl(), c, occ);
          qaregistry.fill(HIST("Tracks/InvQPt"), track.signed1Pt(), c, occ);
          qaregistry.fill(HIST("Tracks/Eta"), track.eta(), c, occ);
          qaregistry.fill(HIST("Tracks/Phi"), phi, c, occ);
        }
      }
      ++nTrk;
    }
    if (fillHis) {
      if constexpr (has_reco_cent<C>) {
        qaregistry.fill(HIST("Tracks/Centrality/NchSel"), nTrk, c, occ);
      } else {
        qaregistry.fill(HIST("Tracks/NchSel"), nTrk, occ);
      }
    }
    return nTrk;
  }

  template <typename C, bool fillHis = false, typename B>
  void countBestTracksExtra(B const& besttracksExtra, float c, float occ)
  {
    for (auto const& etrack : besttracksExtra) {
      if (fillHis) {
        float phi = etrack.phis();
        o2::math_utils::bringTo02Pi(phi);
        if (phi < Czero || TwoPI < phi) {
          continue;
        }
        if constexpr (has_reco_cent<C>) {
          if (gConf.cfgUseTrackParExtra) {
            qaregistry.fill(HIST("Tracks/Centrality/TanLambdaExtra"), etrack.tgl(), c, occ);
            qaregistry.fill(HIST("Tracks/Centrality/InvQPtExtra"), etrack.signed1Pt(), c, occ);
            qaregistry.fill(HIST("Tracks/Centrality/EtaExtra"), etrack.etas(), c, occ);
            qaregistry.fill(HIST("Tracks/Centrality/PhiExtra"), phi, c, occ);
          }
        } else {
          if (gConf.cfgUseTrackParExtra) {
            qaregistry.fill(HIST("Tracks/TanLambdaExtra"), etrack.tgl(), occ);
            qaregistry.fill(HIST("Tracks/InvQPtExtra"), etrack.signed1Pt(), occ);
            qaregistry.fill(HIST("Tracks/EtaExtra"), etrack.etas(), occ);
            qaregistry.fill(HIST("Tracks/PhiExtra"), phi, occ);
          }
        }
      }
    }
  }

  template <typename C, bool fillHis = false, typename T, typename B>
  int countBestTracks(T const& tracks, B const& besttracks, float z,
                      float c, float occ)
  {
    auto nATrk = 0;
    ambiguousTrkIds.reserve(besttracks.size());
    reassignedTrkIds.reserve(besttracks.size());
    for (auto const& atrack : besttracks) {
      if (!isBestTrackSelected(atrack)) {
        continue;
      }
      auto itrack = atrack.template mfttrack_as<T>();
      if (!isTrackSelected(itrack)) {
        continue;
      }
      ambiguousTrkIds.emplace_back(atrack.mfttrackId());
      ++nATrk;
      if (fillHis) {
        float phi = itrack.phi();
        o2::math_utils::bringTo02Pi(phi);
        if (phi < Czero || TwoPI < phi) {
          continue;
        }
        if constexpr (has_reco_cent<C>) {
          registry.fill(HIST("Tracks/Centrality/EtaZvtxBest"), itrack.eta(), z, c, occ);
          registry.fill(HIST("Tracks/Centrality/PhiEtaBest"), phi, itrack.eta(), c, occ);
          qaregistry.fill(HIST("Tracks/Centrality/DCA3d"), itrack.pt(), itrack.eta(), atrack.bestDCAXY(), atrack.bestDCAZ(), c, occ);
          qaregistry.fill(HIST("Tracks/Centrality/NclustersEtaBest"), itrack.nClusters(), itrack.eta(), c, occ);
          qaregistry.fill(HIST("Tracks/Centrality/TrackAmbDegree"), atrack.ambDegree(), c, occ);
        } else {
          registry.fill(HIST("Tracks/EtaZvtxBest"), itrack.eta(), z, occ);
          registry.fill(HIST("Tracks/PhiEtaBest"), phi, itrack.eta(), occ);
          qaregistry.fill(HIST("Tracks/DCA3d"), itrack.pt(), itrack.eta(), atrack.bestDCAXY(), atrack.bestDCAZ(), occ);
          qaregistry.fill(HIST("Tracks/NclustersEtaBest"), itrack.nClusters(), itrack.eta(), occ);
          qaregistry.fill(HIST("Tracks/TrackAmbDegree"), atrack.ambDegree(), occ);
        }
      }

      if (itrack.has_collision() && itrack.collisionId() != atrack.bestCollisionId()) {
        reassignedTrkIds.emplace_back(atrack.mfttrackId());
        if (fillHis) {
          registry.fill(HIST("Tracks/hBestTrkSel"), static_cast<int>(TrkTrkBestSel::trkTrkBestSelNumReassoc));
          float phi = itrack.phi();
          o2::math_utils::bringTo02Pi(phi);
          if (phi < Czero || TwoPI < phi) {
            continue;
          }
          if constexpr (has_reco_cent<C>) {
            qaregistry.fill(HIST("Tracks/Centrality/ReTracksEtaZvtx"), itrack.eta(), itrack.template collision_as<C>().posZ(), c, occ);
            qaregistry.fill(HIST("Tracks/Centrality/ReTracksPhiEta"), phi, itrack.eta(), c, occ);
          } else {
            qaregistry.fill(HIST("Tracks/ReTracksEtaZvtx"), itrack.eta(), itrack.template collision_as<C>().posZ(), occ);
            qaregistry.fill(HIST("Tracks/ReTracksPhiEta"), phi, itrack.eta(), occ);
          }
        }
      }
    }

    for (auto const& track : tracks) {
      if (!isTrackSelected(track)) {
        continue;
      }
      float phi = track.phi();
      o2::math_utils::bringTo02Pi(phi);
      if (phi < Czero || TwoPI < phi) {
        continue;
      }
      if (fillHis) {
        if constexpr (has_reco_cent<C>) {
          qaregistry.fill(HIST("Tracks/Centrality/OrigTracksEtaZvtx"), track.eta(), z, c, occ);
          qaregistry.fill(HIST("Tracks/Centrality/OrigTracksPhiEta"), phi, track.eta(), c, occ);
        } else {
          qaregistry.fill(HIST("Tracks/OrigTracksEtaZvtx"), track.eta(), z, occ);
          qaregistry.fill(HIST("Tracks/OrigTracksPhiEta"), phi, track.eta(), occ);
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
          qaregistry.fill(HIST("Tracks/Centrality/RestTracksEtaZvtx"), track.eta(), z, c, occ);
          qaregistry.fill(HIST("Tracks/Centrality/RestTracksPhiEta"), phi, track.eta(), c, occ);
          // registry.fill(HIST("Tracks/Centrality/EtaZvtxBest"), track.eta(), z, c, occ);
          // registry.fill(HIST("Tracks/Centrality/PhiEtaBest"), phi, track.eta(), c, occ);
          // qaregistry.fill(HIST("Tracks/Centrality/NclustersEtaBest"), track.nClusters(), track.eta(), c, occ);
        } else {
          qaregistry.fill(HIST("Tracks/RestTracksEtaZvtx"), track.eta(), z, occ);
          qaregistry.fill(HIST("Tracks/RestTracksPhiEta"), phi, track.eta(), occ);
          // registry.fill(HIST("Tracks/EtaZvtxBest"), track.eta(), z, occ);
          // registry.fill(HIST("Tracks/PhiEtaBest"), phi, track.eta(), occ);
          // qaregistry.fill(HIST("Tracks/NclustersEtaBest"), track.nClusters(), track.eta(), occ);
        }
      }
    }
    if (fillHis) {
      if constexpr (has_reco_cent<C>) {
        qaregistry.fill(HIST("Tracks/Centrality/NchBestSel"), nATrk, c, occ);
      } else {
        qaregistry.fill(HIST("Tracks/NchBestSel"), nATrk, occ);
      }
    }
    ambiguousTrkIds.clear();
    ambiguousTrkIds.shrink_to_fit();
    reassignedTrkIds.clear();
    reassignedTrkIds.shrink_to_fit();
    return nATrk;
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
      registry.fill(HIST("Events/hNchTVX"), nChrgMc, 0.5);
      return false;
    }
    registry.fill(HIST("Events/hNchTVX"), nChrgMc, 1.5);

    if (nChrgMc == CintZero) {
      return false;
    }

    return true;
  }

  template <typename P>
  bool isParticleSelected(P const& particle)
  {
    if (gConf.cfgUsePrimaries && !particle.isPhysicalPrimary()) {
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
           (phi < ((PIHalf - 0.1) * PI) + trackCuts.phiCut)))
        return false;
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

  void initHadronicRate(CollBCs::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    mRunNumber = bc.runNumber();
    if (gHadronicRate.find(mRunNumber) == gHadronicRate.end()) {
      auto runDuration = ccdb->getRunDuration(mRunNumber);
      mSOR = runDuration.first;
      mMinSeconds = std::floor(mSOR * 1.e-3);               /// round tsSOR to the highest integer lower than tsSOR
      float maxSec = std::ceil(runDuration.second * 1.e-3); /// round tsEOR to the lowest integer higher than tsEOR
      const AxisSpec axisSeconds{static_cast<int>((maxSec - mMinSeconds) / 20.f), 0, maxSec - mMinSeconds, "Seconds since SOR"};
      int hadronicRateBins = static_cast<int>(eventCuts.maxIR - eventCuts.minIR);
      gHadronicRate[mRunNumber] = registry.add<TH2>(Form("HadronicRate/%i", mRunNumber), ";Time since SOR (s);Hadronic rate (kHz)", kTH2D, {axisSeconds, {hadronicRateBins, eventCuts.minIR, eventCuts.maxIR}}).get();
    }
    gCurrentHadronicRate = gHadronicRate[mRunNumber];
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
  void fillHistMC(P const& particles, float c, float occ, float zvtx,
                  bool const gtZeroColl)
  {
    for (auto const& particle : particles) {
      if (!isChrgParticle(particle.pdgCode())) {
        continue;
      }
      if (gConf.cfgUseParticleSel && !isParticleSelected(particle)) {
        continue;
      }

      float phi = particle.phi();
      o2::math_utils::bringTo02Pi(phi);
      if (phi < Czero || TwoPI < phi) {
        continue;
      }
      if constexpr (isCent) {
        registry.fill(HIST("Tracks/Centrality/EtaZvtxGen_t"), particle.eta(),
                      zvtx, c);
        registry.fill(HIST("Tracks/Centrality/PhiEtaGen_t"), phi,
                      particle.eta(), c);
      } else {
        registry.fill(HIST("Tracks/EtaZvtxGen_t"), particle.eta(), zvtx);
        registry.fill(HIST("Tracks/PhiEtaGen_t"), phi, particle.eta());
      }

      if (gtZeroColl) {
        float phi = particle.phi();
        o2::math_utils::bringTo02Pi(phi);
        if (phi < Czero || TwoPI < phi) {
          continue;
        }
        if constexpr (isCent) {
          registry.fill(HIST("Tracks/Centrality/EtaZvtxGen"), particle.eta(),
                        zvtx, c, occ);
          registry.fill(HIST("Tracks/Centrality/PhiEtaGen"), phi,
                        particle.eta(), c, occ);
        } else {
          registry.fill(HIST("Tracks/EtaZvtxGen"), particle.eta(), zvtx, occ);
          registry.fill(HIST("Tracks/PhiEtaGen"), phi, particle.eta(), occ);
        }
      }
    }
  }

  /// @brief process function for general event statistics
  void processTagging(FullBCs const& bcs, CollsCentFT0C const& collisions)
  {
    std::vector<typename std::decay_t<decltype(collisions)>::iterator> cols;
    for (auto const& bc : bcs) {
      if ((bc.selection_bit(aod::evsel::kIsBBT0A) &&
           bc.selection_bit(aod::evsel::kIsBBT0C)) != 0) {
        registry.fill(HIST("hBcSel"), 0);
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
          registry.fill(HIST("hBcSel"), 1);
          if (cols.size() > 1) {
            registry.fill(HIST("hBcSel"), 2);
          }
        }
      }
    }
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processTagging, "Collect event sample stats",
                 true);

  /// @brief process function for counting tracks
  template <typename C>
  void processData(typename C::iterator const& collision,
                   FiltMftTracks const& tracks, CollBCs const& /*bcs*/)
  {
    auto occ = getOccupancy(collision, eventCuts.occupancyEstimator);
    float c = getRecoCent(collision);
    auto bc = collision.template foundBC_as<CollBCs>();
    if constexpr (has_reco_cent<C>) {
      registry.fill(HIST("Events/Centrality/Selection"), 1., c, occ);
    } else {
      registry.fill(HIST("Events/Selection"), 1., occ);
    }

    // const bool isFT0Bad = bc.rct_bit(kFT0Bad);
    // registry.fill(HIST("Events/hRCTSel"), 1.0);
    // if (!isFT0Bad) {
    //   registry.fill(HIST("Events/hRCTSel"), 2.0);
    // }

    if (!isGoodEvent<true>(collision)) {
      return;
    }

    if (gConf.cfgDoIR) {
      initHadronicRate(bc);
      float ir = !gConf.cfgIRSource.value.empty() ? rateFetcher.fetch(ccdb.service, bc.timestamp(), bc.runNumber(), gConf.cfgIRSource, gConf.cfgIRCrashOnNull) * 1.e-3 : -1;
      if constexpr (has_reco_cent<C>) {
        registry.fill(HIST("Events/Centrality/hInteractionRate"), c, occ, ir);
      } else {
        registry.fill(HIST("Events/hInteractionRate"), occ, ir);
      }
      float seconds = bc.timestamp() * 1.e-3 - mMinSeconds;
      if (gConf.cfgUseIRCut && (ir < eventCuts.minIR || ir > eventCuts.maxIR)) { // cut on hadronic rate
        return;
      }
      gCurrentHadronicRate->Fill(seconds, ir);
    }

    auto z = collision.posZ();
    if constexpr (has_reco_cent<C>) {
      registry.fill(HIST("Events/Centrality/Selection"), 2., c, occ);
    } else {
      registry.fill(HIST("Events/Selection"), 2., occ);
    }

    auto nTrk = countTracks<C, true>(tracks, z, c, occ);

    if constexpr (has_reco_cent<C>) {
      registry.fill(HIST("Events/Centrality/NtrkZvtx"), nTrk, z, c, occ);
    } else {
      registry.fill(HIST("Events/NtrkZvtx"), nTrk, z, occ);
    }
  }

  /// @brief process function for counting tracks (based on BestCollisionsFwd3d
  /// table)
  template <typename C>
  void processDatawBestTracks(
    typename C::iterator const& collision, FiltMftTracks const& tracks,
    soa::SmallGroups<aod::BestCollisionsFwd3d> const& besttracks, aod::BestCollisionsFwd3dExtra const& besttracksExtra, CollBCs const& /*bcs*/)
  {
    auto occ = getOccupancy(collision, eventCuts.occupancyEstimator);
    float c = getRecoCent(collision);
    auto bc = collision.template foundBC_as<CollBCs>();
    if constexpr (has_reco_cent<C>) {
      registry.fill(HIST("Events/Centrality/Selection"), 1., c, occ);
    } else {
      registry.fill(HIST("Events/Selection"), 1., occ);
    }

    if (!isGoodEvent<true>(collision)) {
      return;
    }

    if (gConf.cfgDoIR) {
      initHadronicRate(bc);
      float ir = !gConf.cfgIRSource.value.empty() ? rateFetcher.fetch(ccdb.service, bc.timestamp(), bc.runNumber(), gConf.cfgIRSource, gConf.cfgIRCrashOnNull) * 1.e-3 : -1;
      if constexpr (has_reco_cent<C>) {
        registry.fill(HIST("Events/Centrality/hInteractionRate"), c, occ, ir);
      } else {
        registry.fill(HIST("Events/hInteractionRate"), occ, ir);
      }
      float seconds = bc.timestamp() * 1.e-3 - mMinSeconds;
      if (gConf.cfgUseIRCut && (ir < eventCuts.minIR || ir > eventCuts.maxIR)) { // cut on hadronic rate
        return;
      }
      gCurrentHadronicRate->Fill(seconds, ir);
    }

    auto z = collision.posZ();
    if constexpr (has_reco_cent<C>) {
      registry.fill(HIST("Events/Centrality/Selection"), 2., c, occ);
      qaregistry.fill(HIST("Events/Centrality/hZvtxCent"), z, c, occ);
      qaregistry.fill(HIST("Events/Centrality/hCent"), c, occ);
    } else {
      registry.fill(HIST("Events/Selection"), 2., occ);
    }

    auto nBestTrks = countBestTracks<C, true>(tracks, besttracks, z, c, occ);
    countBestTracksExtra<C, true>(besttracksExtra, c, occ);

    if constexpr (has_reco_cent<C>) {
      registry.fill(HIST("Events/Centrality/NtrkZvtxBest"), nBestTrks, z, c,
                    occ);
    } else {
      registry.fill(HIST("Events/NtrkZvtxBest"), nBestTrks, z, occ);
    }
  }

  void processDataInclusive(Colls::iterator const& collision,
                            FiltMftTracks const& tracks, CollBCs const& bcs)
  {
    processData<Colls>(collision, tracks, bcs);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processDataInclusive,
                 "Count tracks (inclusive)", false);

  void processDataCentFT0C(CollsCentFT0C::iterator const& collision,
                           FiltMftTracks const& tracks, CollBCs const& bcs)
  {
    processData<CollsCentFT0C>(collision, tracks, bcs);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processDataCentFT0C,
                 "Count tracks in FT0C centrality bins", false);

  void
    processDataCentFT0CVariant1(CollsCentFT0CVariant1::iterator const& collision,
                                FiltMftTracks const& tracks, CollBCs const& bcs)
  {
    processData<CollsCentFT0CVariant1>(collision, tracks, bcs);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processDataCentFT0CVariant1,
                 "Count tracks in FT0CVariant1 centrality bins", false);

  void processDataCentFT0M(CollsCentFT0M::iterator const& collision,
                           FiltMftTracks const& tracks, CollBCs const& bcs)
  {
    processData<CollsCentFT0M>(collision, tracks, bcs);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processDataCentFT0M,
                 "Count tracks in FT0M centrality bins", false);

  void processDataCentNGlobal(CollsCentNGlobal::iterator const& collision,
                              FiltMftTracks const& tracks, CollBCs const& bcs)
  {
    processData<CollsCentNGlobal>(collision, tracks, bcs);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processDataCentNGlobal,
                 "Count tracks in NGlobal centrality bins", false);

  void processDataCentMFT(CollsCentMFT::iterator const& collision,
                          FiltMftTracks const& tracks, CollBCs const& bcs)
  {
    processData<CollsCentMFT>(collision, tracks, bcs);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processDataCentMFT,
                 "Count tracks in MFT centrality bins", false);

  void processDatawBestTracksInclusive(
    Colls::iterator const& collision, FiltMftTracks const& tracks,
    soa::SmallGroups<aod::BestCollisionsFwd3d> const& besttracks, aod::BestCollisionsFwd3dExtra const& besttracksExtra, CollBCs const& bcs)
  {
    processDatawBestTracks<Colls>(collision, tracks, besttracks, besttracksExtra, bcs);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processDatawBestTracksInclusive,
                 "Count tracks based on BestCollisionsFwd3d table (inclusive)",
                 false);

  void processDatawBestTracksCentFT0C(
    CollsCentFT0C::iterator const& collision, FiltMftTracks const& tracks,
    soa::SmallGroups<aod::BestCollisionsFwd3d> const& besttracks, aod::BestCollisionsFwd3dExtra const& besttracksExtra, CollBCs const& bcs)
  {
    processDatawBestTracks<CollsCentFT0C>(collision, tracks, besttracks, besttracksExtra, bcs);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processDatawBestTracksCentFT0C,
                 "Count tracks in FT0C centrality bins based on BestCollisionsFwd3d table",
                 false);

  void processDatawBestTracksCentFT0CVariant1(
    CollsCentFT0CVariant1::iterator const& collision,
    FiltMftTracks const& tracks,
    soa::SmallGroups<aod::BestCollisionsFwd3d> const& besttracks, aod::BestCollisionsFwd3dExtra const& besttracksExtra, CollBCs const& bcs)
  {
    processDatawBestTracks<CollsCentFT0CVariant1>(collision, tracks, besttracks, besttracksExtra, bcs);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processDatawBestTracksCentFT0CVariant1,
                 "Count tracks in FT0CVariant1 centrality bins based on "
                 "BestCollisionsFwd3d table",
                 false);

  void processDatawBestTracksCentFT0M(
    CollsCentFT0M::iterator const& collision, FiltMftTracks const& tracks,
    soa::SmallGroups<aod::BestCollisionsFwd3d> const& besttracks, aod::BestCollisionsFwd3dExtra const& besttracksExtra, CollBCs const& bcs)
  {
    processDatawBestTracks<CollsCentFT0M>(collision, tracks, besttracks, besttracksExtra, bcs);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processDatawBestTracksCentFT0M,
                 "Count tracks in FT0M centrality bins based on BestCollisionsFwd3d table",
                 false);

  void processDatawBestTracksCentNGlobal(
    CollsCentNGlobal::iterator const& collision, FiltMftTracks const& tracks,
    soa::SmallGroups<aod::BestCollisionsFwd3d> const& besttracks, aod::BestCollisionsFwd3dExtra const& besttracksExtra, CollBCs const& bcs)
  {
    processDatawBestTracks<CollsCentNGlobal>(collision, tracks, besttracks, besttracksExtra, bcs);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processDatawBestTracksCentNGlobal,
                 "Count tracks in NGlobal centrality bins based on "
                 "BestCollisionsFwd3d table",
                 false);

  void processDatawBestTracksCentMFT(
    CollsCentMFT::iterator const& collision, FiltMftTracks const& tracks,
    soa::SmallGroups<aod::BestCollisionsFwd3d> const& besttracks, aod::BestCollisionsFwd3dExtra const& besttracksExtra, CollBCs const& bcs)
  {
    processDatawBestTracks<CollsCentMFT>(collision, tracks, besttracks, besttracksExtra, bcs);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processDatawBestTracksCentMFT,
                 "Count tracks in MFT centrality bins based on BestCollisionsFwd3d table",
                 false);

  Preslice<FiltMcMftTracks> perCol = o2::aod::fwdtrack::collisionId;
  PresliceUnsorted<CollsGenCentFT0C> recColPerMcCol =
    aod::mccollisionlabel::mcCollisionId;
  Partition<FiltParticles> mcSample = (aod::mcparticle::eta < trackCuts.maxEta) && (aod::mcparticle::eta > trackCuts.minEta);

  /// @brief process template function to run on MC gen
  template <typename MC, typename C>
  void processMC(
    typename MC::iterator const& mcCollision,
    soa::SmallGroups<soa::Join<C, aod::McCollisionLabels>> const& collisions,
    FiltParticles const& particles, FiltMcMftTracks const& tracks)
  {
    bool gtZeroColl = false;
    int gtOneColl = 0;

    float cGen = -1;
    if constexpr (has_reco_cent<C>) {
      float crecMin = 105.f;
      for (const auto& collision : collisions) {
        if (isGoodEvent<false>(collision)) {
          float c = getRecoCent(collision);
          if (c < crecMin) {
            crecMin = c;
          }
        }
      }
      if (cGen < 0)
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

    for (auto const& collision : collisions) {
      float occrec = getOccupancy(collision, eventCuts.occupancyEstimator);
      float crec = getRecoCent(collision);

      if constexpr (has_reco_cent<C>) {
        registry.fill(HIST("Events/Centrality/EvtEffGen"), 1., crec, occrec);
      } else {
        registry.fill(HIST("Events/EvtEffGen"), 1., occrec);
      }

      if (isGoodEvent<false>(collision)) {
        gtZeroColl = true;
        ++gtOneColl;
        auto z = collision.posZ();

        if constexpr (has_reco_cent<C>) {
          registry.fill(HIST("Events/Centrality/EvtEffGen"), 2., crec, occrec);
          registry.fill(HIST("Events/Centrality/hRecCent"), crec, occrec);
          registry.fill(HIST("Events/Centrality/hRecZvtxCent"), z, crec,
                        occrec);
        } else {
          registry.fill(HIST("Events/EvtEffGen"), 2., occrec);
        }

        auto perColSample = tracks.sliceBy(perCol, collision.globalIndex());
        auto nTrkRec = countTracks<C, true>(perColSample, z, crec, occrec);

        if constexpr (has_reco_cent<C>) {
          qaregistry.fill(HIST("Events/Centrality/ZvtxDiff"),
                          collision.posZ() - mcCollision.posZ(), crec);
        } else {
          qaregistry.fill(HIST("Events/ZvtxDiff"),
                          collision.posZ() - mcCollision.posZ());
        }

        if (eventCuts.useZDiffCut) {
          if (std::abs(collision.posZ() - mcCollision.posZ()) >
              eventCuts.maxZvtxDiff) {
            continue;
          }
        }

        if constexpr (has_reco_cent<C>) {
          registry.fill(HIST("Events/Centrality/NtrkZvtxGen"), nTrkRec,
                        collision.posZ(), crec, occrec);
        } else {
          registry.fill(HIST("Events/NtrkZvtxGen"), nTrkRec, collision.posZ(),
                        occrec);
        }
      }
    }

    if constexpr (has_reco_cent<C>) {
      registry.fill(HIST("Events/Centrality/EvtEffGen"), 3., cGen, occGen);
    } else {
      registry.fill(HIST("Events/EvtEffGen"), 3., occGen);
    }

    auto perCollMCsample = mcSample->sliceByCached(
      aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);
    auto nchrg = countPart(perCollMCsample);
    auto zvtxMC = mcCollision.posZ();

    if (gtOneColl > 1) {
      if constexpr (has_reco_cent<C>) {
        qaregistry.fill(HIST("Events/Centrality/SplitMult"), nchrg, zvtxMC, cGen);
      } else {
        qaregistry.fill(HIST("Events/SplitMult"), nchrg, zvtxMC);
      }
    }

    auto nCharged = countPart(particles);
    if constexpr (has_reco_cent<C>) {
      registry.fill(HIST("Events/Centrality/NtrkZvtxGen_t"), nCharged, zvtxMC,
                    cGen);
    } else {
      registry.fill(HIST("Events/NtrkZvtxGen_t"), nCharged, zvtxMC);
    }

    fillHistMC<has_reco_cent<C>>(particles, cGen, occGen, zvtxMC, gtZeroColl);

    if (collisions.size() == 0) {
      if constexpr (has_reco_cent<C>) {
        qaregistry.fill(HIST("Events/Centrality/NotFoundEventZvtx"),
                        mcCollision.posZ(), cGen);
      } else {
        qaregistry.fill(HIST("Events/NotFoundEventZvtx"), mcCollision.posZ());
      }
    }
  }

  void processMCInclusive(
    aod::McCollisions::iterator const& mccollision,
    soa::SmallGroups<soa::Join<Colls, aod::McCollisionLabels>> const& collisions,
    FiltParticles const& particles, FiltMcMftTracks const& tracks)
  {
    processMC<aod::McCollisions, Colls>(mccollision, collisions, particles,
                                        tracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processMCInclusive,
                 "Count MC particles (inclusive)", false);

  void processMCCentFT0C(
    aod::McCollisions::iterator const& mccollision,
    soa::SmallGroups<soa::Join<CollsCentFT0C, aod::McCollisionLabels>> const& collisions,
    FiltParticles const& particles, FiltMcMftTracks const& tracks)
  {
    processMC<aod::McCollisions, CollsCentFT0C>(mccollision, collisions,
                                                particles, tracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processMCCentFT0C,
                 "Count MC particles in FT0C centrality bins", false);

  void processMCCentFT0CVariant1(
    aod::McCollisions::iterator const& mccollision,
    soa::SmallGroups<soa::Join<CollsCentFT0CVariant1,
                               aod::McCollisionLabels>> const& collisions,
    FiltParticles const& particles, FiltMcMftTracks const& tracks)
  {
    processMC<aod::McCollisions, CollsCentFT0CVariant1>(mccollision, collisions,
                                                        particles, tracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processMCCentFT0CVariant1,
                 "Count MC particles in FT0CVariant1 centrality bins", false);

  void processMCCentFT0M(
    aod::McCollisions::iterator const& mccollision,
    soa::SmallGroups<soa::Join<CollsCentFT0M, aod::McCollisionLabels>> const& collisions,
    FiltParticles const& particles, FiltMcMftTracks const& tracks)
  {
    processMC<aod::McCollisions, CollsCentFT0M>(mccollision, collisions,
                                                particles, tracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processMCCentFT0M,
                 "Count MC particles in FT0M centrality bins", false);

  void processMCCentNGlobal(
    aod::McCollisions::iterator const& mccollision,
    soa::SmallGroups<soa::Join<CollsCentNGlobal,
                               aod::McCollisionLabels>> const& collisions,
    FiltParticles const& particles, FiltMcMftTracks const& tracks)
  {
    processMC<aod::McCollisions, CollsCentNGlobal>(mccollision, collisions,
                                                   particles, tracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processMCCentNGlobal,
                 "Count MC particles in NGlobal centrality bins", false);

  void processMCCentMFT(
    aod::McCollisions::iterator const& mccollision,
    soa::SmallGroups<soa::Join<CollsCentMFT, aod::McCollisionLabels>> const& collisions,
    FiltParticles const& particles, FiltMcMftTracks const& tracks)
  {
    processMC<aod::McCollisions, CollsCentMFT>(mccollision, collisions,
                                               particles, tracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processMCCentMFT,
                 "Count MC particles in MFT centrality bins", false);

  PresliceUnsorted<aod::BestCollisionsFwd3d> perColU =
    aod::fwdtrack::bestCollisionId;

  /// @brief process template function to run on MC truth using
  /// aod::BestCollisionsFwd3d tracks
  template <typename MC, typename C>
  void processMCwBestTracks(
    typename MC::iterator const& mcCollision,
    soa::SmallGroups<soa::Join<C, aod::McCollisionLabels>> const& collisions,
    FiltParticles const& particles, FiltMcMftTracks const& tracks,
    FiltBestTracks const& besttracks)
  {
    bool gtZeroColl = false;
    float cGen = -1;
    if constexpr (has_reco_cent<C>) {
      float crecMin = 105.f;
      for (const auto& collision : collisions) {
        if (isGoodEvent<false>(collision)) {
          float c = getRecoCent(collision);
          if (c < crecMin) {
            crecMin = c;
          }
        }
      }
      if (cGen < 0)
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

    for (auto const& collision : collisions) {
      auto occrec = getOccupancy(collision, eventCuts.occupancyEstimator);
      float crec = getRecoCent(collision);

      if constexpr (has_reco_cent<C>) {
        registry.fill(HIST("Events/Centrality/EvtEffGen"), 1., crec, occrec);
      } else {
        registry.fill(HIST("Events/EvtEffGen"), 1., occrec);
      }

      if (isGoodEvent<false>(collision)) {
        gtZeroColl = true;
        auto z = collision.posZ();

        if constexpr (has_reco_cent<C>) {
          registry.fill(HIST("Events/Centrality/EvtEffGen"), 2., crec, occrec);
        } else {
          registry.fill(HIST("Events/EvtEffGen"), 2., occrec);
        }

        auto perCollisionSample =
          tracks.sliceBy(perCol, collision.globalIndex());
        auto perCollisionASample =
          besttracks.sliceBy(perColU, collision.globalIndex());
        auto nTrkRec = countBestTracks<C, false>(
          perCollisionSample, perCollisionASample, z, crec,
          collision.trackOccupancyInTimeRange());

        if constexpr (has_reco_cent<C>) {
          registry.fill(HIST("Events/Centrality/NtrkZvtxGen"), nTrkRec, z,
                        crec, occrec);
        } else {
          registry.fill(HIST("Events/NtrkZvtxGen"), nTrkRec, z, occrec);
        }
      }
    }

    if constexpr (has_reco_cent<C>) {
      registry.fill(HIST("Events/Centrality/EvtEffGen"), 3., cGen, occGen);
    } else {
      registry.fill(HIST("Events/EvtEffGen"), 3., occGen);
    }

    auto zvtxMC = mcCollision.posZ();
    auto nCharged = countPart(particles);
    if constexpr (has_reco_cent<C>) {
      registry.fill(HIST("Events/Centrality/NtrkZvtxGen_t"), nCharged, zvtxMC,
                    cGen);
    } else {
      registry.fill(HIST("Events/NtrkZvtxGen_t"), nCharged, zvtxMC);
    }

    fillHistMC<has_reco_cent<C>>(particles, cGen, occGen, zvtxMC, gtZeroColl);

    if (collisions.size() == 0) {
      if constexpr (has_reco_cent<C>) {
        qaregistry.fill(HIST("Events/Centrality/NotFoundEventZvtx"),
                        mcCollision.posZ(), cGen);
      } else {
        qaregistry.fill(HIST("Events/NotFoundEventZvtx"), mcCollision.posZ());
      }
    }
  }

  void processMCwBestTracksInclusive(
    aod::McCollisions::iterator const& mccollision,
    soa::SmallGroups<soa::Join<Colls, aod::McCollisionLabels>> const& collisions,
    FiltParticles const& particles, FiltMcMftTracks const& tracks,
    FiltBestTracks const& besttracks)
  {
    processMCwBestTracks<aod::McCollisions, Colls>(
      mccollision, collisions, particles, tracks, besttracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processMCwBestTracksInclusive,
                 "Count MC particles using aod::BestCollisionsFwd3d (inclusive)",
                 false);

  void processMCwBestTracksCentFT0C(
    aod::McCollisions::iterator const& mccollision,
    soa::SmallGroups<soa::Join<CollsCentFT0C, aod::McCollisionLabels>> const& collisions,
    FiltParticles const& particles, FiltMcMftTracks const& tracks,
    FiltBestTracks const& besttracks)
  {
    processMCwBestTracks<aod::McCollisions, CollsCentFT0C>(
      mccollision, collisions, particles, tracks, besttracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processMCwBestTracksCentFT0C,
                 "Count MC particles in FT0C centrality bins using aod::BestCollisionsFwd3d",
                 false);

  void processMCwBestTracksCentFT0CVariant1(
    aod::McCollisions::iterator const& mccollision,
    soa::SmallGroups<soa::Join<CollsCentFT0CVariant1,
                               aod::McCollisionLabels>> const& collisions,
    FiltParticles const& particles, FiltMcMftTracks const& tracks,
    FiltBestTracks const& besttracks)
  {
    processMCwBestTracks<aod::McCollisions, CollsCentFT0CVariant1>(
      mccollision, collisions, particles, tracks, besttracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processMCwBestTracksCentFT0CVariant1,
                 "Count MC particles in FT0CVariant1 centrality bins using "
                 "aod::BestCollisionsFwd3d",
                 false);

  void processMCwBestTracksCentFT0M(
    aod::McCollisions::iterator const& mccollision,
    soa::SmallGroups<soa::Join<CollsCentFT0M, aod::McCollisionLabels>> const& collisions,
    FiltParticles const& particles, FiltMcMftTracks const& tracks,
    FiltBestTracks const& besttracks)
  {
    processMCwBestTracks<aod::McCollisions, CollsCentFT0M>(
      mccollision, collisions, particles, tracks, besttracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processMCwBestTracksCentFT0M,
                 "Count MC particles in FT0M centrality bins using aod::BestCollisionsFwd3d",
                 false);

  void processMCwBestTracksCentNGlobal(
    aod::McCollisions::iterator const& mccollision,
    soa::SmallGroups<soa::Join<CollsCentNGlobal,
                               aod::McCollisionLabels>> const& collisions,
    FiltParticles const& particles, FiltMcMftTracks const& tracks,
    FiltBestTracks const& besttracks)
  {
    processMCwBestTracks<aod::McCollisions, CollsCentNGlobal>(
      mccollision, collisions, particles, tracks, besttracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processMCwBestTracksCentNGlobal,
                 "Count MC particles in NGlobal centrality bins using "
                 "aod::BestCollisionsFwd3d",
                 false);

  void processMCwBestTracksCentMFT(
    aod::McCollisions::iterator const& mccollision,
    soa::SmallGroups<soa::Join<CollsCentMFT, aod::McCollisionLabels>> const& collisions,
    FiltParticles const& particles, FiltMcMftTracks const& tracks,
    FiltBestTracks const& besttracks)
  {
    processMCwBestTracks<aod::McCollisions, CollsCentMFT>(
      mccollision, collisions, particles, tracks, besttracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processMCwBestTracksCentMFT,
                 "Count MC particles in MFT centrality bins using aod::BestCollisionsFwd3d",
                 false);

  using ParticlesI = soa::Join<aod::McParticles, aod::ParticlesToMftTracks>;
  Partition<ParticlesI> primariesI = (aod::mcparticle::flags & (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary) == (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary && (aod::mcparticle::eta < trackCuts.maxEta) && (aod::mcparticle::eta > trackCuts.minEta);

  Preslice<aod::McParticles> perMCCol = aod::mcparticle::mcCollisionId;

  template <typename C, typename MC>
  void processTrkEffIdxBest(
    typename soa::Join<C, aod::McCollisionLabels> const& collisions,
    MC const& /*mccollisions*/, ParticlesI const& /*particles*/,
    BestTracksMC const& atracks)
  {
    for (auto const& collision : collisions) {
      if (!isGoodEvent<false>(collision)) {
        continue;
      }
      if (!collision.has_mcCollision()) {
        continue;
      }

      float crec = getRecoCent(collision);
      auto mcCollision = collision.mcCollision();

      if (eventCuts.useZDiffCut) {
        if (std::abs(collision.posZ() - mcCollision.posZ()) > eventCuts.maxZvtxDiff) {
          continue;
        }
      }

      auto partsPerCol = primariesI->sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);
      partsPerCol.bindExternalIndices(&atracks);
      for (auto const& particle : partsPerCol) {
        if (!isChrgParticle(particle.pdgCode())) {
          continue;
        }
        if (gConf.cfgUseParticleSel && !isParticleSelected(particle)) {
          continue;
        }

        // MC gen
        if constexpr (has_reco_cent<C>) {
          if (particle.eta() > trackCuts.minEta && particle.eta() < trackCuts.maxEta) {
            if (std::abs(mcCollision.posZ()) < eventCuts.maxZvtx) {
              qaregistry.fill(HIST("Tracks/Centrality/hPtEtaEffGenBest"), particle.pt(), particle.eta(), crec);
            }
          }
        } else {
          if (particle.eta() > trackCuts.minEta && particle.eta() < trackCuts.maxEta) {
            if (std::abs(mcCollision.posZ()) < eventCuts.maxZvtx) {
              qaregistry.fill(HIST("Tracks/hPtEtaEffGenBest"), particle.pt(), particle.eta());
            }
          }
        }
        // MC rec
        if (particle.has_mfttracks()) {
          auto iscounted = false;
          auto ncnt = 0;
          auto relatedTracks = particle.template mfttracks_as<BestTracksMC>();
          for (auto const& track : relatedTracks) {
            if (!isBestTrackSelected<false>(track)) {
              continue;
            }
            ++ncnt;

            if constexpr (has_reco_cent<C>) {
              if (track.eta() > trackCuts.minEta && track.eta() < trackCuts.maxEta) {
                if (!iscounted) { // primaries
                  if (std::abs(mcCollision.posZ()) < eventCuts.maxZvtx) {
                    qaregistry.fill(HIST("Tracks/Centrality/hPtEtaEffPrimBest"), particle.pt(), particle.eta(), crec);
                  }
                  iscounted = true;
                }
              }
              if (ncnt > 1) { // secondaries
                if (track.eta() > trackCuts.minEta && track.eta() < trackCuts.maxEta) {
                  qaregistry.fill(HIST("Tracks/Centrality/hPtEtaEffSecBest"), particle.pt(), particle.eta(), crec);
                }
              }
            } else {
              if (track.eta() > trackCuts.minEta && track.eta() < trackCuts.maxEta) {
                if (!iscounted) { // primaries
                  if (std::abs(mcCollision.posZ()) < eventCuts.maxZvtx) {
                    qaregistry.fill(HIST("Tracks/hPtEtaEffPrimBest"), particle.pt(), particle.eta());
                  }
                  iscounted = true;
                }
              }
              if (ncnt > 1) { // secondaries
                if (track.eta() > trackCuts.minEta && track.eta() < trackCuts.maxEta) {
                  qaregistry.fill(HIST("Tracks/hPtEtaEffSecBest"), particle.pt(), particle.eta());
                }
              }
            }

            if constexpr (has_reco_cent<C>) {
              qaregistry.fill(HIST("Tracks/Centrality/NmftTrkPerPartBest"), ncnt, crec);
            } else {
              qaregistry.fill(HIST("Tracks/NmftTrkPerPartBest"), ncnt);
            }

            if (relatedTracks.size() > 1) { // duplicates
              if constexpr (has_reco_cent<C>) {
                qaregistry.fill(HIST("Tracks/Centrality/hPtEtaEffGenDuplBest"), particle.pt(), particle.eta(), crec);
                for (auto const& track : relatedTracks) {
                  qaregistry.fill(HIST("Tracks/Centrality/hPtEtaEffDuplBest"), track.pt(), track.eta(), crec);
                }
              } else {
                qaregistry.fill(HIST("Tracks/hPtEtaEffGenDuplBest"), particle.pt(), particle.eta());
                for (auto const& track : relatedTracks) {
                  qaregistry.fill(HIST("Tracks/hPtEtaEffDuplBest"), track.pt(), track.eta());
                }
              }
            }
          }
        } else {
          // MC FAKES
          if constexpr (has_reco_cent<C>) {
            if (particle.eta() > trackCuts.minEta && particle.eta() < trackCuts.maxEta) {
              if (std::abs(mcCollision.posZ()) < eventCuts.maxZvtx) {
                qaregistry.fill(HIST("Tracks/Centrality/hPtEtaEffGenFakeBest"), particle.pt(), particle.eta(), crec);
              }
            }
          } else {
            if (particle.eta() > trackCuts.minEta && particle.eta() < trackCuts.maxEta) {
              if (std::abs(mcCollision.posZ()) < eventCuts.maxZvtx) {
                qaregistry.fill(HIST("Tracks/hPtEtaEffGenFakeBest"), particle.pt(), particle.eta());
              }
            }
          }
        }
      }
    }
  }

  void processTrkEffIdxBestInlusive(
    soa::Join<Colls, aod::McCollisionLabels> const& collisions,
    aod::McCollisions const& mccollisions, ParticlesI const& particles,
    BestTracksMC const& atracks)
  {
    processTrkEffIdxBest<Colls, aod::McCollisions>(collisions, mccollisions, particles, atracks);
  }
  PROCESS_SWITCH(DndetaMFTPbPb, processTrkEffIdxBestInlusive, "Process tracking efficiency best (inclusive, indexed)", false);

  void processTrkEffIdxBestCentFT0C(
    soa::Join<CollsCentFT0C, aod::McCollisionLabels> const& collisions,
    aod::McCollisions const& mccollisions, ParticlesI const& particles,
    BestTracksMC const& atracks)
  {
    processTrkEffIdxBest<CollsCentFT0C, aod::McCollisions>(collisions, mccollisions, particles, atracks);
  }
  PROCESS_SWITCH(DndetaMFTPbPb, processTrkEffIdxBestCentFT0C, "Process tracking efficiency best (in FT0C centrality bins, indexed)", false);

  /// @brief process template function to calculate tracking efficiency (indexed
  /// as particle-to-MFT-tracks)
  template <typename C, typename MC>
  void processTrkEffIdx(
    typename soa::Join<C, aod::McCollisionLabels> const& collisions,
    MC const& /*mccollisions*/, ParticlesI const& /*particles*/,
    MFTTracksLabeled const& tracks)
  {
    for (auto const& collision : collisions) {
      if (!isGoodEvent<false>(collision)) {
        continue;
      }
      if (!collision.has_mcCollision()) {
        continue;
      }

      float crec = getRecoCent(collision);
      auto occrec = getOccupancy(collision, eventCuts.occupancyEstimator);
      auto mcCollision = collision.mcCollision();

      auto partsPerCol = primariesI->sliceByCached(
        aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);
      partsPerCol.bindExternalIndices(&tracks);

      for (auto const& particle : partsPerCol) {
        if (!isChrgParticle(particle.pdgCode())) {
          continue;
        }
        if (gConf.cfgUseParticleSel && !isParticleSelected(particle)) {
          continue;
        }

        // MC gen
        if constexpr (has_reco_cent<C>) {
          if (particle.eta() > trackCuts.minEta && particle.eta() < trackCuts.maxEta) {
            if (std::abs(mcCollision.posZ()) < eventCuts.maxZvtx) {
              qaregistry.fill(HIST("Tracks/Centrality/hPtEtaEffGen"), particle.pt(), particle.eta(), crec, occrec);
            }
          }
        } else {
          if (particle.eta() > trackCuts.minEta && particle.eta() < trackCuts.maxEta) {
            if (std::abs(mcCollision.posZ()) < eventCuts.maxZvtx) {
              qaregistry.fill(HIST("Tracks/hPtEtaEffGen"), particle.pt(), particle.eta(), occrec);
            }
          }
        }
        // MC rec
        if (particle.has_mfttracks()) {
          auto iscounted = false;
          auto ncnt = 0;
          auto relatedTracks = particle.template mfttracks_as<MFTTracksLabeled>();
          for (auto const& track : relatedTracks) {
            ++ncnt;
            if constexpr (has_reco_cent<C>) {
              if (track.eta() > trackCuts.minEta && track.eta() < trackCuts.maxEta) {
                if (!iscounted) { // primaries
                  if (std::abs(mcCollision.posZ()) < eventCuts.maxZvtx) {
                    qaregistry.fill(HIST("Tracks/Centrality/hPtEtaEffPrim"), particle.pt(), particle.eta(), crec, occrec);
                  }
                  iscounted = true;
                }
              }
              if (ncnt > 1) { // secondaries
                if (track.eta() > trackCuts.minEta && track.eta() < trackCuts.maxEta) {
                  qaregistry.fill(HIST("Tracks/Centrality/hPtEtaEffSec"), particle.pt(), particle.eta(), crec, occrec);
                }
              }
            } else {
              if (track.eta() > trackCuts.minEta && track.eta() < trackCuts.maxEta) {
                if (!iscounted) { // primaries
                  if (std::abs(mcCollision.posZ()) < eventCuts.maxZvtx) {
                    qaregistry.fill(HIST("Tracks/hPtEtaEffPrim"), particle.pt(), particle.eta(), occrec);
                  }
                  iscounted = true;
                }
              }
              if (ncnt > 1) { // secondaries
                if (track.eta() > trackCuts.minEta && track.eta() < trackCuts.maxEta) {
                  qaregistry.fill(HIST("Tracks/hPtEtaEffSec"), particle.pt(), particle.eta(), occrec);
                }
              }
            }
          }

          if constexpr (has_reco_cent<C>) {
            qaregistry.fill(HIST("Tracks/Centrality/NmftTrkPerPart"), ncnt, crec, occrec);
          } else {
            qaregistry.fill(HIST("Tracks/NmftTrkPerPart"), ncnt, occrec);
          }

          if (relatedTracks.size() > 1) { // duplicates
            if constexpr (has_reco_cent<C>) {
              qaregistry.fill(HIST("Tracks/Centrality/hPtEtaEffGenDupl"), particle.pt(), particle.eta(), crec, occrec);
              for (auto const& track : relatedTracks) {
                qaregistry.fill(HIST("Tracks/Centrality/hPtEtaEffDupl"), track.pt(), track.eta(), crec, occrec);
              }
            } else {
              qaregistry.fill(HIST("Tracks/hPtEtaEffGenDupl"), particle.pt(), particle.eta(), occrec);
              for (auto const& track : relatedTracks) {
                qaregistry.fill(HIST("Tracks/hPtEtaEffDupl"), track.pt(), track.eta(), occrec);
              }
            }
          }
        }
      }
    }
  }

  void processTrkEffIdxInlusive(
    soa::Join<Colls, aod::McCollisionLabels> const& collisions,
    aod::McCollisions const& mccollisions, ParticlesI const& particles,
    MFTTracksLabeled const& tracks)
  {
    processTrkEffIdx<Colls, aod::McCollisions>(collisions, mccollisions,
                                               particles, tracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processTrkEffIdxInlusive,
                 "Process tracking efficiency (inclusive, indexed)", false);

  void processTrkEffIdxCentFT0C(
    soa::Join<CollsCentFT0C, aod::McCollisionLabels> const& collisions,
    aod::McCollisions const& mccollisions, ParticlesI const& particles,
    MFTTracksLabeled const& tracks)
  {
    processTrkEffIdx<CollsCentFT0C, aod::McCollisions>(collisions, mccollisions,
                                                       particles, tracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processTrkEffIdxCentFT0C,
                 "Process tracking efficiency (in FT0C centrality bins, indexed)", false);

  /// @brief process function to calculate tracking efficiency (indexed) based
  /// on BestCollisionsFwd3d in FT0C bins
  template <typename C, typename MC>
  void processTrkEffBest(
    typename soa::Join<C, aod::McCollisionLabels>::iterator const& collision,
    MC const& /*mccollisions*/, FiltParticles const& particles,
    FiltMcMftTracks const& tracks,
    soa::SmallGroups<aod::BestCollisionsFwd3d> const& besttracks)
  {
    if (!isGoodEvent<false>(collision)) {
      return;
    }
    if (!collision.has_mcCollision()) {
      return;
    }

    float crec = getRecoCent(collision);
    auto occrec = getOccupancy(collision, eventCuts.occupancyEstimator);
    auto mcCollision = collision.mcCollision();
    auto partsPerCol = particles.sliceByCached(
      aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);

    for (auto const& particle : partsPerCol) {
      if (!isChrgParticle(particle.pdgCode())) {
        continue;
      }
      if (gConf.cfgUseParticleSel && !isParticleSelected(particle)) {
        continue;
      }
      if constexpr (has_reco_cent<C>) {
        if (particle.eta() > trackCuts.minEta && particle.eta() < trackCuts.maxEta) {
          qaregistry.fill(HIST("Tracks/Centrality/hPtPhiEtaZvtxEffBestGen"), particle.pt(), particle.phi(), particle.eta(), mcCollision.posZ(), crec, occrec);
        }
      } else {
        if (particle.eta() > trackCuts.minEta && particle.eta() < trackCuts.maxEta) {
          qaregistry.fill(HIST("Tracks/hPtPhiEtaZvtxEffBestGen"), particle.pt(), particle.phi(), particle.eta(), mcCollision.posZ(), occrec);
        }
      }
    }

    ambiguousTrkIdsMC.reserve(besttracks.size());
    reassignedTrkIdsMC.reserve(besttracks.size());

    for (auto const& track : besttracks) {
      ambiguousTrkIdsMC.emplace_back(track.mfttrackId());
      if (!isBestTrackSelected<false>(track)) {
        continue;
      }
      auto itrack = track.mfttrack_as<FiltMcMftTracks>();
      if (itrack.collisionId() != track.bestCollisionId()) {
        reassignedTrkIdsMC.emplace_back(track.mfttrackId());
      }
      if (!isTrackSelected<false>(itrack)) {
        continue;
      }
      if (itrack.has_mcParticle()) {
        auto particle = itrack.mcParticle_as<FiltParticles>();
        if (itrack.eta() > trackCuts.minEta && itrack.eta() < trackCuts.maxEta) {
          if constexpr (has_reco_cent<C>) {
            qaregistry.fill(HIST("Tracks/Centrality/hPtPhiEtaZvtxEffBestRec"),
                            particle.pt(), particle.phi(), particle.eta(),
                            mcCollision.posZ(), crec, occrec);
          } else {
            qaregistry.fill(HIST("Tracks/hPtPhiEtaZvtxEffBestRec"), particle.pt(),
                            particle.phi(), particle.eta(), mcCollision.posZ(),
                            occrec);
          }
        }
      } else {
        if constexpr (has_reco_cent<C>) {
          qaregistry.fill(HIST("Tracks/Centrality/hPtEffBestFakeRec"), itrack.pt(), itrack.phi(), itrack.eta(), mcCollision.posZ(), crec, occrec);
        } else {
          qaregistry.fill(HIST("Tracks/hPtEffBestFakeRec"), itrack.pt(), itrack.phi(), itrack.eta(), mcCollision.posZ(), occrec);
        }
      }
    }

    for (auto const& track : tracks) {
      if (std::find(ambiguousTrkIdsMC.begin(), ambiguousTrkIdsMC.end(), track.globalIndex()) != ambiguousTrkIdsMC.end()) {
        continue;
      }
      if (std::find(reassignedTrkIdsMC.begin(), reassignedTrkIdsMC.end(), track.globalIndex()) != reassignedTrkIdsMC.end()) {
        continue;
      }
      if (!isTrackSelected<false>(track)) {
        continue;
      }
      if (track.has_mcParticle()) {
        auto particle = track.mcParticle_as<FiltParticles>();
        if (track.eta() > trackCuts.minEta && track.eta() < trackCuts.maxEta) {
          if constexpr (has_reco_cent<C>) {
            qaregistry.fill(HIST("Tracks/Centrality/hPtPhiEtaZvtxEffBestRec"),
                            particle.pt(), particle.phi(), particle.eta(),
                            mcCollision.posZ(), crec, occrec);
          } else {
            qaregistry.fill(HIST("Tracks/hPtPhiEtaZvtxEffBestRec"), particle.pt(),
                            particle.phi(), particle.eta(), mcCollision.posZ(),
                            occrec);
          }
        }
      } else {
        if constexpr (has_reco_cent<C>) {
          qaregistry.fill(HIST("Tracks/Centrality/hPtEffBestFakeRec"), track.pt(), track.phi(), track.eta(), mcCollision.posZ(), crec, occrec);
        } else {
          qaregistry.fill(HIST("Tracks/hPtEffBestFakeRec"), track.pt(), track.phi(), track.eta(), mcCollision.posZ(), occrec);
        }
      }
    }
    ambiguousTrkIdsMC.clear();
    ambiguousTrkIdsMC.shrink_to_fit();
    reassignedTrkIdsMC.clear();
    reassignedTrkIdsMC.shrink_to_fit();
  }

  void processTrkEffBestInclusive(
    soa::Join<Colls, aod::McCollisionLabels>::iterator const& collision,
    aod::McCollisions const& mccollisions, FiltParticles const& particles,
    FiltMcMftTracks const& tracks,
    soa::SmallGroups<aod::BestCollisionsFwd3d> const& besttracks)
  {
    processTrkEffBest<Colls, aod::McCollisions>(collision, mccollisions,
                                                particles, tracks, besttracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processTrkEffBestInclusive,
                 "Process tracking efficiency (inclusive, based on BestCollisionsFwd3d)",
                 false);

  void processTrkEffBestCentFT0C(
    soa::Join<CollsCentFT0C, aod::McCollisionLabels>::iterator const& collision,
    aod::McCollisions const& mccollisions, FiltParticles const& particles,
    FiltMcMftTracks const& tracks,
    soa::SmallGroups<aod::BestCollisionsFwd3d> const& besttracks)
  {
    processTrkEffBest<CollsCentFT0C, aod::McCollisions>(
      collision, mccollisions, particles, tracks, besttracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processTrkEffBestCentFT0C,
                 "Process tracking efficiency (in FT0 centrality bins, based "
                 "on BestCollisionsFwd3d)",
                 false);

  Preslice<FiltMcMftTracks> filtMcTrkperCol = o2::aod::fwdtrack::collisionId;

  /// @brief process function to calculate MC efficiency and fraction of fake
  /// tracks
  template <typename C, typename MC>
  void processEfficiency(
    typename soa::Join<C, aod::McCollisionLabels> const& collisions,
    MC const& /*mccollisions*/, FiltParticles const& /*particles*/,
    FiltMcMftTracks const& tracks)
  {
    for (auto const& collision : collisions) {
      if (!isGoodEvent<false>(collision)) {
        continue;
      }

      float crec = getRecoCent(collision);
      auto occrec = getOccupancy(collision, eventCuts.occupancyEstimator);
      auto mcCollision = collision.mcCollision();
      auto perColTrks =
        tracks.sliceBy(filtMcTrkperCol, collision.globalIndex());

      for (auto const& track : perColTrks) {
        if (!isTrackSelected<false>(track)) {
          continue;
        }
        if (track.has_mcParticle()) {
          auto particle = track.template mcParticle_as<FiltParticles>();
          if constexpr (has_reco_cent<C>) {
            qaregistry.fill(HIST("Tracks/Centrality/hEffRec"), particle.pt(),
                            particle.phi(), particle.eta(), mcCollision.posZ(),
                            crec, occrec);
          } else {
            qaregistry.fill(HIST("Tracks/hEffRec"), particle.pt(),
                            particle.phi(), particle.eta(), mcCollision.posZ(),
                            crec, occrec);
          }
        } else {
          if constexpr (has_reco_cent<C>) {
            qaregistry.fill(HIST("Tracks/Centrality/hEffFake"), track.pt(),
                            track.phi(), track.eta(), mcCollision.posZ(), crec,
                            occrec);
          } else {
            qaregistry.fill(HIST("Tracks/hEffFake"), track.pt(), track.phi(),
                            track.eta(), mcCollision.posZ(), crec, occrec);
          }
        }
      }
    }
  }

  void processEfficiencyInclusive(
    soa::Join<Colls, aod::McCollisionLabels> const& collisions,
    aod::McCollisions const& mccollisions, FiltParticles const& particles,
    FiltMcMftTracks const& tracks)
  {
    processEfficiency<Colls, aod::McCollisions>(collisions, mccollisions,
                                                particles, tracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processEfficiencyInclusive,
                 "Process efficiencies (inclusive)", false);

  void processEfficiencyCentFT0C(
    soa::Join<CollsCentFT0C, aod::McCollisionLabels> const& collisions,
    aod::McCollisions const& mccollisions, FiltParticles const& particles,
    FiltMcMftTracks const& tracks)
  {
    processEfficiency<CollsCentFT0C, aod::McCollisions>(
      collisions, mccollisions, particles, tracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processEfficiencyCentFT0C,
                 "Process efficiencies in FT0C centrality bins", false);

  template <typename C>
  void processCorrelationwBestTracks(typename C::iterator const& collision, FiltMftTracks const& /*tracks*/, soa::SmallGroups<aod::BestCollisionsFwd3d> const& besttracks)
  {
    if (!isGoodEvent<false>(collision)) {
      return;
    }

    auto nBestTrks = 0;
    for (auto const& atrack : besttracks) {
      if (gConf.cfgUseTrackSel && !isBestTrackSelected<false>(atrack)) {
        continue;
      }
      auto itrack = atrack.template mfttrack_as<FiltMftTracks>();
      if (itrack.eta() < trackCuts.minEta || itrack.eta() > trackCuts.maxEta) {
        continue;
      }
      if (gConf.cfgUseTrackSel && !isTrackSelected<false>(itrack)) {
        continue;
      }
      nBestTrks++;
    }
    qaregistry.fill(HIST("Events/hMultMFTvsFT0A"), nBestTrks, collision.multFT0A());
    qaregistry.fill(HIST("Events/hMultMFTvsFT0C"), nBestTrks, collision.multFT0C());
    qaregistry.fill(HIST("Events/hNPVtracksVsFT0C"), collision.multNTracksPV(), collision.multFT0C());
    qaregistry.fill(HIST("Events/hMultMFTvsFV0A"), nBestTrks, collision.multFV0A());
    qaregistry.fill(HIST("Events/hNPVtracksVsMultMFT"), collision.multNTracksPV(), nBestTrks);
  }

  void processCorrelationwBestTracksInclusive(CollsCorr::iterator const& collision, FiltMftTracks const& tracks, soa::SmallGroups<aod::BestCollisionsFwd3d> const& besttracks)
  {
    processCorrelationwBestTracks<CollsCorr>(collision, tracks, besttracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processCorrelationwBestTracksInclusive, "Do correlation study based on BestCollisionsFwd3d table", false);

  int getQuadrantPhi(float phi)
  {
    if (phi >= Czero && phi < Cninety) {
      return 0;
    }
    if (phi >= Cninety && phi <= ConeHeighty) {
      return 1;
    }
    if (phi >= -ConeHeighty && phi < -Cninety) {
      return 2;
    }
    if (phi >= -Cninety && phi < Czero) {
      return 3;
    }
    return -1;
  }

  template <typename T>
  int getQuadrantTrack(T const& track)
  {
    float phi = static_cast<float>(track.phi()) * ConeHeighty / PI;
    return getQuadrantPhi(phi);
  }

  /// @brief process function to check MFT alignment (based on BestCollisionsFwd3d table)
  template <typename C>
  void processAlignment(typename C::iterator const& collision,
                        FiltMftTracks const& /*tracks*/,
                        soa::SmallGroups<aod::BestCollisionsFwd3d> const& besttracks,
                        CollBCs const& /*bcs*/
  )
  {
    auto bc = collision.template foundBC_as<CollBCs>();
    if (!isGoodEvent<true>(collision)) {
      return;
    }

    if (gConf.cfgDoIR) {
      initHadronicRate(bc);
      float ir = !gConf.cfgIRSource.value.empty() ? rateFetcher.fetch(ccdb.service, bc.timestamp(), bc.runNumber(), gConf.cfgIRSource, gConf.cfgIRCrashOnNull) * 1.e-3 : -1;
      float seconds = bc.timestamp() * 1.e-3 - mMinSeconds;
      if (gConf.cfgUseIRCut && (ir < eventCuts.minIR || ir > eventCuts.maxIR)) { // cut on hadronic rate
        return;
      }
      gCurrentHadronicRate->Fill(seconds, ir);
    }

    auto z = collision.posZ();
    for (auto const& atrack : besttracks) {
      if (!isBestTrackSelected(atrack)) {
        continue;
      }
      auto itrack = atrack.template mfttrack_as<FiltMftTracks>();
      if (!isTrackSelected(itrack)) {
        continue;
      }

      int quadrant = getQuadrantTrack(itrack);
      if (quadrant < 0) {
        continue;
      }
      std::get<std::shared_ptr<TH2>>(hAlignment[0][quadrant][0]["DCA_x_vs_z"])->Fill(z, atrack.bestDCAXY());
    }
  }

  void processAlignmentInclusive(Colls::iterator const& collision,
                                 FiltMftTracks const& tracks,
                                 soa::SmallGroups<aod::BestCollisionsFwd3d> const& besttracks,
                                 CollBCs const& bcs)
  {
    processAlignment<Colls>(collision, tracks, besttracks, bcs);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processAlignmentInclusive, "Check MFT alignment using tracks based on BestCollisionsFwd3d table (inclusive)", false);

  /// @brief process function to calculate signal loss based on MC
  void processEventAndSignalLossCentFT0C(CollsMCExtraMult::iterator const& mcCollision,
                                         soa::SmallGroups<soa::Join<CollsCentFT0C, aod::McCollisionLabels>> const& collisions,
                                         FiltParticles const& particles)
  {
    registry.fill(HIST("Events/Centrality/hNRecCollsSigEvtLoss"), 1.f, collisions.size());

    if (gConf.cfgUseInelgt0 && !mcCollision.isInelGt0()) {
      return;
    }
    if (eventCuts.useZVtxCutMC && (std::abs(mcCollision.posZ()) >= eventCuts.maxZvtx)) {
      return;
    }

    bool gtZeroColl = false;
    auto maxNcontributors = -1;
    auto centrality = -1;
    for (auto const& collision : collisions) {
      if (!isGoodEvent<false>(collision)) {
        continue;
      }
      if (std::abs(collision.posZ()) >= eventCuts.maxZvtx) {
        continue;
      }
      if (maxNcontributors < collision.numContrib()) {
        maxNcontributors = collision.numContrib();
        centrality = getRecoCent(collision);
      }
      gtZeroColl = true;
    }

    auto perCollMCsample = mcSample->sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);
    auto multMCNParticlesEtaMFT = countPart(perCollMCsample);

    registry.fill(HIST("Events/Centrality/EvtSigLossStatus"), 1., centrality);
    registry.fill(HIST("Events/Centrality/hMultGenVsCent"), centrality, mcCollision.multMCFT0C());
    registry.fill(HIST("Events/Centrality/hMultGenVsCentNParticlesEta05"), centrality, mcCollision.multMCNParticlesEta05());
    registry.fill(HIST("Events/Centrality/hMultGenVsCentNParticlesEtaMFT"), centrality, multMCNParticlesEtaMFT);

    if (collisions.size() == 0) {
      registry.fill(HIST("Events/Centrality/EvtSigLossStatus"), 3., centrality);
    }

    if (gtZeroColl) {
      registry.fill(HIST("Events/Centrality/EvtSigLossStatus"), 2., centrality);
      registry.fill(HIST("Events/Centrality/hMultGenVsCentRec"), centrality, mcCollision.multMCFT0C());
      registry.fill(HIST("Events/Centrality/hMultGenVsCentRecNParticlesEta05"), centrality, mcCollision.multMCNParticlesEta05());
      registry.fill(HIST("Events/Centrality/hMultGenVsCentRecNParticlesEtaMFT"), centrality, multMCNParticlesEtaMFT);
    }

    for (auto const& particle : particles) {
      if (!isChrgParticle(particle.pdgCode())) {
        continue;
      }
      if (gConf.cfgUseParticleSel && !isParticleSelected(particle)) {
        continue;
      }

      float phi = particle.phi();
      o2::math_utils::bringTo02Pi(phi);
      if (phi < Czero || TwoPI < phi) {
        continue;
      }

      registry.fill(HIST("Tracks/Centrality/EtaCentVsMultGen_t"), particle.eta(), centrality, mcCollision.multMCFT0C());
      registry.fill(HIST("Tracks/Centrality/EtaGen_t"), particle.eta(), centrality);

      if (gtZeroColl) {
        float phi = particle.phi();
        o2::math_utils::bringTo02Pi(phi);
        if (phi < Czero || TwoPI < phi) {
          continue;
        }
        registry.fill(HIST("Tracks/Centrality/EtaCentVsMultGen"), particle.eta(), centrality, mcCollision.multMCFT0C());
        registry.fill(HIST("Tracks/Centrality/EtaGen"), particle.eta(), centrality);
      }
    }
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processEventAndSignalLossCentFT0C, "Signal/event loss based on MC (in FT0C centrality bins)", false);

  Preslice<FiltMftTracks> filtTrkperCol = o2::aod::fwdtrack::collisionId;

  /// @brief process template function for MC QA checks
  template <typename C>
  void processMcQA(typename soa::Join<C, aod::McCollisionLabels> const& collisions,
                   MFTTracksLabeled const& tracks,
                   aod::AmbiguousMFTTracks const& atracks,
                   aod::McCollisions const& mcCollisions,
                   FiltParticles const& /*particles*/)
  {
    for (const auto& collision : collisions) {
      float crec = getRecoCent(collision);
      auto occrec = getOccupancy(collision, eventCuts.occupancyEstimator);

      if constexpr (has_reco_cent<C>) {
        qaregistry.fill(HIST("Events/Centrality/hRecPerGenColls"),
                        static_cast<float>(collisions.size()) /
                          mcCollisions.size(),
                        crec, occrec);
      } else {
        qaregistry.fill(HIST("Events/hRecPerGenColls"),
                        static_cast<float>(collisions.size()) /
                          mcCollisions.size(),
                        occrec);
      }

      if (!isGoodEvent<false>(collision)) {
        return;
      }

      auto trkPerColl = tracks.sliceBy(filtTrkperCol, collision.globalIndex());
      uint ntracks{0u}, nAtracks{0u};
      for (const auto& track : trkPerColl) {
        ntracks++;
        for (const auto& atrack : atracks) {
          if (atrack.mfttrackId() == track.globalIndex()) {
            nAtracks++;
            break;
          }
        }
      }
      if constexpr (has_reco_cent<C>) {
        qaregistry.fill(HIST("Tracks/Centrality/hNmftTrks"), ntracks, crec,
                        occrec);
        qaregistry.fill(HIST("Tracks/Centrality/hFracAmbiguousMftTrks"),
                        static_cast<float>(nAtracks) / ntracks, crec, occrec);
      } else {
        qaregistry.fill(HIST("Tracks/hNmftTrks"), ntracks, occrec);
        qaregistry.fill(HIST("Tracks/hFracAmbiguousMftTrks"),
                        static_cast<float>(nAtracks) / ntracks, occrec);
      }
    }
  }

  void processMcQAInclusive(
    soa::Join<Colls, aod::McCollisionLabels> const& collisions,
    MFTTracksLabeled const& tracks, aod::AmbiguousMFTTracks const& atracks,
    aod::McCollisions const& mcCollisions, FiltParticles const& particles)
  {
    processMcQA<Colls>(collisions, tracks, atracks, mcCollisions, particles);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processMcQAInclusive,
                 "Process MC QA checks (inclusive)", false);

  void processMcQACentFT0C(
    soa::Join<CollsCentFT0C, aod::McCollisionLabels> const& collisions,
    MFTTracksLabeled const& tracks, aod::AmbiguousMFTTracks const& atracks,
    aod::McCollisions const& mcCollisions, FiltParticles const& particles)
  {
    processMcQA<CollsCentFT0C>(collisions, tracks, atracks, mcCollisions,
                               particles);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processMcQACentFT0C,
                 "Process MC QA checks (in FT0 centrality bins)", false);

  /// @brief process function to check ambiguous tracks
  void processCheckAmbiguousMftTracks(aod::Collisions const&, MftTracksWColls const& tracks)
  {
    for (auto const& track : tracks) {
      auto trkCollId = track.has_collision() ? track.collisionId() : -1;
      auto ids = track.compatibleCollIds();
      qaregistry.fill(HIST("Tracks/hAmbTrackType"), static_cast<int>(AmbTrkType::kAll));
      if (ids.empty() || (ids.size() == 1 && trkCollId == ids[0])) {
        qaregistry.fill(HIST("Tracks/hMftTracksAmbDegreeWithTrivial"), track.compatibleCollIds().size());
        if (ids.empty()) {
          qaregistry.fill(HIST("Tracks/hAmbTrackType"), static_cast<int>(AmbTrkType::kOrphan));
        }
        if (ids.size() == 1 && trkCollId == ids[0]) {
          qaregistry.fill(HIST("Tracks/hAmbTrackType"), static_cast<int>(AmbTrkType::kNonAmb));
        }
        continue;
      }
      qaregistry.fill(HIST("Tracks/hMftTracksAmbDegree"), track.compatibleCollIds().size());

      if (track.compatibleCollIds().size() > 0) {
        if (track.compatibleCollIds().size() == 1) {
          if (track.collisionId() != track.compatibleCollIds()[0]) {
            qaregistry.fill(HIST("Tracks/hAmbTrackType"), static_cast<int>(AmbTrkType::kAmb));
          } else {
            qaregistry.fill(HIST("Tracks/hAmbTrackType"), static_cast<int>(AmbTrkType::kNonAmbSame));
          }
        } else {
          qaregistry.fill(HIST("Tracks/hAmbTrackType"), static_cast<int>(AmbTrkType::kAmbGt1));
        }
      }
    }
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processCheckAmbiguousMftTracks, "Process checks for Ambiguous MFT tracks (inclusive)", false);

  // Preslice<aod::McParticles> perColMc = aod::mcparticle::mcCollisionId;
  void processTimeAssocMC(CollsMCExtra const& mcCollisions,
                          CollisionsWithMCLabels const& collisions,
                          MftTracksWCollsMC const& tracks,
                          aod::McParticles const& /*particles*/,
                          aod::McCollisionLabels const& labels)
  {
    const auto& nRecoColls = collisions.size();
    LOG(info) << "reconstructed collisions: " << nRecoColls;
    const auto& nMcColls = mcCollisions.size();
    LOG(info) << "MC collisions: " << nMcColls;
    const auto& nLabels = labels.size();
    LOG(info) << "collision labels: " << nLabels;

    std::unordered_map<int64_t, int64_t> mapRecToMc;
    mapRecToMc.clear();
    mapRecToMc.reserve(nRecoColls);

    // std::unordered_map<int64_t, float> mapVtxXrec;
    mapVtxXrec.clear();
    mapVtxXrec.reserve(nRecoColls);
    // std::unordered_map<int64_t, float> mapVtxYrec;
    mapVtxYrec.clear();
    mapVtxYrec.reserve(nRecoColls);
    // std::unordered_map<int64_t, float> mapVtxZrec;
    mapVtxZrec.clear();
    mapVtxZrec.reserve(nRecoColls);

    if (nRecoColls <= CintZero) {
      return;
    }

    auto maxNcontributors = -1;
    auto bestCollIndex = -1;
    for (auto const& collision : collisions) {
      if (maxNcontributors < collision.numContrib()) {
        maxNcontributors = collision.numContrib();
        bestCollIndex = collision.globalIndex();
        mapVtxXrec.emplace(collision.globalIndex(), collision.posX());
        mapVtxYrec.emplace(collision.globalIndex(), collision.posY());
        mapVtxZrec.emplace(collision.globalIndex(), collision.posZ());
        mapRecToMc.emplace(collision.globalIndex(), collision.mcCollisionId());
      }
    }
    LOG(info) << "mapRecToMc size: " << mapRecToMc.size();
    LOG(info) << "mapVtxXrec size: " << mapVtxXrec.size();

    std::unordered_map<int64_t, float> mapVtxXgen;
    mapVtxXgen.clear();
    mapVtxXgen.reserve(nMcColls);
    std::unordered_map<int64_t, float> mapVtxYgen;
    mapVtxYgen.reserve(nMcColls);
    mapVtxYgen.clear();
    std::unordered_map<int64_t, float> mapVtxZgen;
    mapVtxZgen.clear();
    mapVtxZgen.reserve(nMcColls);

    for (const auto& mcCollision : mcCollisions) {
      mapVtxXgen.emplace(mcCollision.globalIndex(), mcCollision.posX());
      mapVtxYgen.emplace(mcCollision.globalIndex(), mcCollision.posY());
      mapVtxZgen.emplace(mcCollision.globalIndex(), mcCollision.posZ());
    }
    LOG(info) << "mapVtxXgen size: " << mapVtxXgen.size();

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
      registry.fill(HIST("TimeAssocMC/hTimeAssocMCEventStatus"), static_cast<int>(ReAssocMCEventStatus::kEvtReAsHasMcColl));

      int64_t recCollId = collision.globalIndex();
      auto itMC = mapRecToMc.find(recCollId);
      if (itMC == mapRecToMc.end()) {
        nNoMC++;
        LOGP(debug, "collison {} has no MC coll", recCollId);
        continue;
      }
      auto mcCollision = collision.mcCollision_as<CollsMCExtra>();
      if (gConf.cfgRemoveSplitVertex && (bestCollIndex != collision.globalIndex())) {
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
        if (trkCollId != recCollId) { // check if track is associated to rec coll
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
        if (ids.size() > 0) {
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

        if (track.collisionId() >= 0 && track.has_mcParticle() && track.mcMask() == 0) {
          auto itMCTrk = mapRecToMc.find(trkCollId);
          const auto& mcPart = track.mcParticle();
          if (!isChrgParticle(mcPart.pdgCode())) {
            continue;
          }
          if (gConf.cfgUseParticleSel && !isParticleSelected(mcPart)) {
            continue;
          }
          int64_t mcPartId = mcPart.mcCollisionId();

          // check if rec vertex is available in MC collisions
          for (const auto& mcTrkId : mapRecToMc) {
            if (mcTrkId.second == mcPartId) {
              isTrueVtx = true;
              break;
            }
          }

          // check if there is good or bad collision
          if (itMCTrk != mapRecToMc.end()) {
            int mcTrkCollId = itMCTrk->second;
            if (mcPartId == mcTrkCollId) { // particle.mcCollisionId == collision.mcCollisionId -> good vtx
              vtxFlag = static_cast<int>(VertexStatusMC::kGood);
            } else { // wrong vtx
              vtxFlag = static_cast<int>(VertexStatusMC::kBad);
            }
          }

          if (mapVtxXrec.find(trkCollId) == mapVtxXrec.end()) {
            continue;
          }
          if (mapVtxYrec.find(trkCollId) == mapVtxYrec.end()) {
            continue;
          }
          if (mapVtxZrec.find(trkCollId) == mapVtxZrec.end()) {
            continue;
          }
          if (mapRecToMc.find(trkCollId) == mapRecToMc.end()) {
            continue;
          }
          vtxX = mapVtxXrec.find(trkCollId)->second;
          vtxY = mapVtxYrec.find(trkCollId)->second;
          vtxZ = mapVtxZrec.find(trkCollId)->second;
          // LOGP(info, "\t ---> \t .... \t vtxZrec: {} - collision.posZ(): {}", vtxZrec, collision.posZ());
          int64_t mcCollIdRec = mapRecToMc.find(trkCollId)->second;
          // int64_t mcCollId = itMC->second;
          // LOGP(info, "\t ---> \t .... \t mcCollIdRec: {} - mcCollId: {} - bestMCCol: {}", mcCollIdRec, mcCollId, bestMCCol);
          if (mapVtxXgen.find(mcCollIdRec) == mapVtxXgen.end()) {
            continue;
          }
          if (mapVtxYgen.find(mcCollIdRec) == mapVtxYgen.end()) {
            continue;
          }
          if (mapVtxZgen.find(mcCollIdRec) == mapVtxZgen.end()) {
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
        }

        registry.fill(HIST("TimeAssocMC/VtxStatus"), vtxFlag);
        registry.fill(HIST("TimeAssocMC/hVertexResV1"), deltaXv1, deltaYv1, deltaZv1);
        registry.fill(HIST("TimeAssocMC/hVertexResV2"), deltaXv2, deltaYv2, deltaZv2);

        registry.fill(HIST("TimeAssocMC/hAmbTrkTypeAssocFlag"), static_cast<int>(AmbTrkTypeAssocFlag::kSel));

        if (isTrueVtx) {
          registry.fill(HIST("TimeAssocMC/hTimeAssocCheckVtxType"), static_cast<int>(ReassocCheckVtxType::kIsTrueVtxAllTrue));
        } else {
          registry.fill(HIST("TimeAssocMC/hTimeAssocCheckVtxType"), static_cast<int>(ReassocCheckVtxType::kIsTrueVtxAllFalse));
        }
        if (vtxFlag == static_cast<int>(VertexStatusMC::kGood)) {
          if (isTrueVtx) {
            registry.fill(HIST("TimeAssocMC/hVTXkSelGoodVtxTrue"), deltaXv2, deltaYv2, deltaZv2);
            registry.fill(HIST("TimeAssocMC/hAmbTrkTypeAssocFlag"), static_cast<int>(AmbTrkTypeAssocFlag::kSelGoodVtxTrue));
            registry.fill(HIST("TimeAssocMC/hTimeAssocCheckVtxType"), static_cast<int>(ReassocCheckVtxType::kIsTrueVtxVsGoodVtxTrue));
          } else {
            registry.fill(HIST("TimeAssocMC/hTimeAssocCheckVtxType"), static_cast<int>(ReassocCheckVtxType::kIsTrueVtxVsGoodVtxFalse));
          }
        }
        if (vtxFlag == static_cast<int>(VertexStatusMC::kBad)) {
          if (isTrueVtx) {
            registry.fill(HIST("TimeAssocMC/hVTXkSelGoodVtxBad"), deltaXv2, deltaYv2, deltaZv2);
            registry.fill(HIST("TimeAssocMC/hAmbTrkTypeAssocFlag"), static_cast<int>(AmbTrkTypeAssocFlag::kSelGoodVtxBad));
            registry.fill(HIST("TimeAssocMC/hTimeAssocCheckVtxType"), static_cast<int>(ReassocCheckVtxType::kIsTrueVtxVsBadVtxTrue));
          } else {
            registry.fill(HIST("TimeAssocMC/hTimeAssocCheckVtxType"), static_cast<int>(ReassocCheckVtxType::kIsTrueVtxVsBadVtxFalse));
          }
        }

        if (ids.size() > 0) {
          if (ids.size() == 1) {
            if (trkCollId == ids[0]) { // non ambiguous
              registry.fill(HIST("TimeAssocMC/hVTXkSelNonAmbAll"), deltaXv2, deltaYv2, deltaZv2);
              registry.fill(HIST("TimeAssocMC/hAmbTrkTypeAssocFlag"), static_cast<int>(AmbTrkTypeAssocFlag::kSelNonAmbAll));
              if (vtxFlag == static_cast<int>(VertexStatusMC::kGood)) {
                if (isTrueVtx) {
                  registry.fill(HIST("TimeAssocMC/hVTXkSelNonAmbGoodVtxTrue"), deltaXv2, deltaYv2, deltaZv2);
                  registry.fill(HIST("TimeAssocMC/hAmbTrkTypeAssocFlag"), static_cast<int>(AmbTrkTypeAssocFlag::kSelNonAmbGoodVtxTrue));
                }
              }
              if (vtxFlag == static_cast<int>(VertexStatusMC::kBad)) {
                if (isTrueVtx) {
                  registry.fill(HIST("TimeAssocMC/hVTXkSelNonAmbGoodVtxBad"), deltaXv2, deltaYv2, deltaZv2);
                  registry.fill(HIST("TimeAssocMC/hAmbTrkTypeAssocFlag"), static_cast<int>(AmbTrkTypeAssocFlag::kSelNonAmbGoodVtxBad));
                }
              }
            } else if (trkCollId != ids[0]) {
              registry.fill(HIST("TimeAssocMC/hVTXkSelAmbAll"), deltaXv2, deltaYv2, deltaZv2);
              registry.fill(HIST("TimeAssocMC/hAmbTrkTypeAssocFlag"), static_cast<int>(AmbTrkTypeAssocFlag::kSelAmbAll));
              if (vtxFlag == static_cast<int>(VertexStatusMC::kGood)) {
                if (isTrueVtx) {
                  registry.fill(HIST("TimeAssocMC/hVTXkSelAmbGoodVtxTrue"), deltaXv2, deltaYv2, deltaZv2);
                  registry.fill(HIST("TimeAssocMC/hAmbTrkTypeAssocFlag"), static_cast<int>(AmbTrkTypeAssocFlag::kSelAmbGoodVtxTrue));
                }
              }
              if (vtxFlag == static_cast<int>(VertexStatusMC::kBad)) {
                if (isTrueVtx) {
                  registry.fill(HIST("TimeAssocMC/hVTXkSelAmbGoodVtxBad"), deltaXv2, deltaYv2, deltaZv2);
                  registry.fill(HIST("TimeAssocMC/hAmbTrkTypeAssocFlag"), static_cast<int>(AmbTrkTypeAssocFlag::kSelAmbGoodVtxBad));
                }
              }
            } else { // non ambiguous (extra)
              registry.fill(HIST("TimeAssocMC/hVTXkSelNonAmbSameAll"), deltaXv2, deltaYv2, deltaZv2);
              registry.fill(HIST("TimeAssocMC/hAmbTrkTypeAssocFlag"), static_cast<int>(AmbTrkTypeAssocFlag::kSelNonAmbSameAll));
              if (vtxFlag == static_cast<int>(VertexStatusMC::kGood)) {
                if (isTrueVtx) {
                  registry.fill(HIST("TimeAssocMC/hVTXkSelNonAmbSameGoodVtxTrue"), deltaXv2, deltaYv2, deltaZv2);
                  registry.fill(HIST("TimeAssocMC/hAmbTrkTypeAssocFlag"), static_cast<int>(AmbTrkTypeAssocFlag::kSelNonAmbSameGoodVtxTrue));
                }
              }
              if (vtxFlag == static_cast<int>(VertexStatusMC::kBad)) {
                if (isTrueVtx) {
                  registry.fill(HIST("TimeAssocMC/hVTXkSelNonAmbSameGoodVtxBad"), deltaXv2, deltaYv2, deltaZv2);
                  registry.fill(HIST("TimeAssocMC/hAmbTrkTypeAssocFlag"), static_cast<int>(AmbTrkTypeAssocFlag::kSelNonAmbSameGoodVtxBad));
                }
              }
            }
          } else { // ambiguous
            registry.fill(HIST("TimeAssocMC/hVTXkSelAmbGt1All"), deltaXv2, deltaYv2, deltaZv2);
            registry.fill(HIST("TimeAssocMC/hAmbTrkTypeAssocFlag"), static_cast<int>(AmbTrkTypeAssocFlag::kSelAmbGt1All));
            if (vtxFlag == static_cast<int>(VertexStatusMC::kGood)) {
              if (isTrueVtx) {
                registry.fill(HIST("TimeAssocMC/hVTXkSelAmbGt1GoodVtxTrue"), deltaXv2, deltaYv2, deltaZv2);
                registry.fill(HIST("TimeAssocMC/hAmbTrkTypeAssocFlag"), static_cast<int>(AmbTrkTypeAssocFlag::kSelAmbGt1GoodVtxTrue));
              }
            }
            if (vtxFlag == static_cast<int>(VertexStatusMC::kBad)) {
              if (isTrueVtx) {
                registry.fill(HIST("TimeAssocMC/hVTXkSelAmbGt1GoodVtxBad"), deltaXv2, deltaYv2, deltaZv2);
                registry.fill(HIST("TimeAssocMC/hAmbTrkTypeAssocFlag"), static_cast<int>(AmbTrkTypeAssocFlag::kSelAmbGt1GoodVtxBad));
              }
            }
          }
        } else {
          registry.fill(HIST("TimeAssocMC/hAmbTrkTypeAssocFlag"), static_cast<int>(AmbTrkTypeAssocFlag::kSelOrphanNull));
        }
      }
    }
    LOG(info) << "No MC: " << nNoMC;
  }
  PROCESS_SWITCH(DndetaMFTPbPb, processTimeAssocMC, "process MC: check MFT tracks in compatible collisions (time-associasted)", false);

  void processTimeAssocWithReassocMC(CollsMCExtra const& mcCollisions,
                                     CollisionsWithMCLabels const& collisions,
                                     FiltMftTracks const& /*tracks*/,
                                     BestTracks3dWCollsMC const& besttracks,
                                     aod::McParticles const& /*particles*/,
                                     aod::McCollisionLabels const& labels,
                                     ExtBCs const& bcs)
  {
    if (bcs.size() == 0) {
      return;
    }
    auto bc = bcs.begin();
    initCCDB(bc);

    const auto& nRecoColls = collisions.size();
    LOG(info) << "reconstructed collisions: " << nRecoColls;
    const auto& nMcColls = mcCollisions.size();
    LOG(info) << "MC collisions: " << nMcColls;
    // LOGP(debug, "MC col {} has {} reco cols", mcCollision.globalIndex(), collisions.size());
    const auto& nLabels = labels.size();
    LOG(info) << "collision labels: " << nLabels;

    std::unordered_map<int64_t, int64_t> mapRecToMc;
    mapRecToMc.reserve(nRecoColls);
    std::unordered_map<int64_t, int64_t> mapMcToRec;
    mapMcToRec.reserve(nRecoColls);
    // std::unordered_map<int64_t, float> mapVtxXrec;
    mapVtxXrec.reserve(nRecoColls);
    // std::unordered_map<int64_t, float> mapVtxYrec;
    mapVtxYrec.reserve(nRecoColls);
    // std::unordered_map<int64_t, float> mapVtxZrec;
    mapVtxZrec.reserve(nRecoColls);

    if (nRecoColls <= CintZero) {
      return;
    }

    auto maxNcontributors = -1;
    auto bestCollIndex = -1;
    for (auto const& collision : collisions) {
      if (maxNcontributors < collision.numContrib()) {
        maxNcontributors = collision.numContrib();
        bestCollIndex = collision.globalIndex();
        mapVtxXrec.emplace(collision.globalIndex(), collision.posX());
        mapVtxYrec.emplace(collision.globalIndex(), collision.posY());
        mapVtxZrec.emplace(collision.globalIndex(), collision.posZ());
        mapRecToMc.emplace(collision.globalIndex(), collision.mcCollisionId());
        mapMcToRec.emplace(collision.mcCollisionId(), collision.globalIndex());
      }
    }
    LOG(info) << "mapRecToMc size: " << mapRecToMc.size();
    LOG(info) << "mapVtxXrec size: " << mapVtxXrec.size();

    std::unordered_map<int64_t, float> mapVtxXgen;
    mapVtxXgen.reserve(nMcColls);
    std::unordered_map<int64_t, float> mapVtxYgen;
    mapVtxYgen.reserve(nMcColls);
    std::unordered_map<int64_t, float> mapVtxZgen;
    mapVtxZgen.reserve(nMcColls);

    for (const auto& mcCollision : mcCollisions) {
      mapVtxXgen.emplace(mcCollision.globalIndex(), mcCollision.posX());
      mapVtxYgen.emplace(mcCollision.globalIndex(), mcCollision.posY());
      mapVtxZgen.emplace(mcCollision.globalIndex(), mcCollision.posZ());
    }
    LOG(info) << "mapVtxXgen size: " << mapVtxXgen.size();

    int nNoMC{0};
    for (const auto& collision : collisions) {
      int64_t recCollId = collision.globalIndex();
      auto itMC = mapRecToMc.find(recCollId);
      if (itMC == mapRecToMc.end()) {
        nNoMC++;
        LOGP(debug, "collison {} has no MC coll", recCollId);
        continue;
      }
      if (!collision.has_mcCollision()) {
        continue;
      }
      auto mcCollision = collision.mcCollision_as<CollsMCExtra>();
      // if (collision.globalIndex() != mcCollision.bestCollisionIndex()) {
      //   continue;
      // }
      if (gConf.cfgRemoveSplitVertex && (bestCollIndex != collision.globalIndex())) {
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
        // auto itrack = atrack.mfttrack_as<FiltMftTracks>();
        // if (!isTrackSelected(itrack)) {
        // continue;
        // }
        if (gConf.cfgUseTrackSel && !isTrackSelected<true>(atrack)) {
          continue;
        }
        float phi = atrack.phi();
        o2::math_utils::bringTo02Pi(phi);
        if (phi < Czero || TwoPI < phi) {
          continue;
        }
        auto trkBestCollId = atrack.has_collision() ? atrack.bestCollisionId() : -1;
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

        if (trkBestCollId >= 0 && atrack.has_mcParticle() && atrack.mcMask() == 0) {
          auto itRecToMc = mapRecToMc.find(trkBestCollId); // try mfttrackId ???
          int64_t mcPartId = atrack.mcParticle().mcCollisionId();
          auto const& idCompColl = atrack.compatibleCollIds();

          // check if rec vertex is available in MC collisions
          for (const auto& mcTrkId : mapRecToMc) {
            if (mcTrkId.second == mcPartId) {
              isTrueVtx = true;
              break;
            }
          }

          // check good rec vertex of corresponding mc coll is available in compatible rec coll
          if (!idCompColl.empty()) {
            for (auto const& id : idCompColl) {
              auto itMcCollId = mapRecToMc.find(id);
              if (itMcCollId != mapRecToMc.end()) {
                if (itMcCollId->second == mcPartId) {
                  isRecGoodMatch = true;
                  break;
                }
              }
            }
          }

          // check if there is good or bad collision
          if (itRecToMc != mapRecToMc.end()) {
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
          if (mapVtxXrec.find(trkBestCollId) == mapVtxXrec.end()) {
            continue;
          }
          if (mapVtxYrec.find(trkBestCollId) == mapVtxYrec.end()) {
            continue;
          }
          if (mapVtxZrec.find(trkBestCollId) == mapVtxZrec.end()) {
            continue;
          }
          if (mapRecToMc.find(trkBestCollId) == mapRecToMc.end()) {
            continue;
          }
          vtxXbest = mapVtxXrec.find(trkBestCollId)->second;
          vtxYbest = mapVtxYrec.find(trkBestCollId)->second;
          vtxZbest = mapVtxZrec.find(trkBestCollId)->second;
          // LOGP(info, "\t ---> \t .... \t vtxZrec: {} - collision.posZ(): {}", vtxZrec, collision.posZ());
          int64_t mcCollIdRec = mapRecToMc.find(trkBestCollId)->second;
          // int64_t mcCollId = itMC->second;
          // LOGP(info, "\t ---> \t .... \t mcCollIdRec: {} - mcCollId: {} - bestMCCol: {}", mcCollIdRec, mcCollId, bestMCCol);
          if (mapVtxXgen.find(mcCollIdRec) == mapVtxXgen.end()) {
            continue;
          }
          if (mapVtxYgen.find(mcCollIdRec) == mapVtxYgen.end()) {
            continue;
          }
          if (mapVtxZgen.find(mcCollIdRec) == mapVtxZgen.end()) {
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
        }

        if (vtxFlag == static_cast<int>(VertexStatusMC::kBad)) {
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
  PROCESS_SWITCH(DndetaMFTPbPb, processTimeAssocWithReassocMC, "process MC: check MFT tracks in reassociation with compatible collisions", false);

  Partition<MftTracksWCollsMC> tracksInAcc = (aod::fwdtrack::eta < trackCuts.maxEta) && (aod::fwdtrack::eta > trackCuts.minEta);

  void processAssocMC(CollisionsWithMCLabels const& collisions,
                      MftTracksWCollsMC const& tracks,
                      aod::McParticles const& /*particles*/,
                      aod::McCollisions const& /*mccollisions*/
  )
  {
    const auto& nRecoColls = collisions.size();
    // Generated evets with >= 1 reco collisions
    if (nRecoColls > CintZero) {
      auto maxNcontributors = -1;
      auto bestCollIndex = -1;
      for (const auto& collision : collisions) {
        if (maxNcontributors < collision.numContrib()) {
          maxNcontributors = collision.numContrib();
          bestCollIndex = collision.globalIndex();
        }
      }
      for (const auto& collision : collisions) {
        if (!isGoodEvent<true>(collision)) {
          continue;
        }
        // Select collisions with the largest number of contributors
        if (gConf.cfgRemoveSplitVertex && (bestCollIndex != collision.globalIndex())) {
          continue;
        }
        if (!collision.has_mcCollision()) {
          continue;
        }
        auto mcCollision = collision.template mcCollision_as<aod::McCollisions>();

        if (eventCuts.useZVtxCutMC && (std::abs(mcCollision.posZ()) >= eventCuts.maxZvtx)) {
          continue;
        }

        int nTrk = 0, nFakeTrk = 0, nGoodTrk = 0;
        for (const auto& track : tracks) {
          uint index = uint(track.collisionId() >= 0);
          if (track.has_mcParticle() && track.mcMask() == 0) {
            // auto particle = track.mcParticle_as<FiltParticles>();
            const auto& particle = track.mcParticle();
            auto trkCollId = track.has_collision() ? track.collisionId() : -1;
            auto ids = track.compatibleCollIds();
            qaregistry.fill(HIST("TrackToColl/histTrackNumColls"), ids.size());
            qaregistry.fill(HIST("TrackToColl/hAmbTrackType"), static_cast<int>(AmbTrkType::kAll));

            if (ids.empty() || (ids.size() == 1 && trkCollId == ids[0])) {
              if (ids.empty()) {
                qaregistry.fill(HIST("TrackToColl/hAmbTrackType"), static_cast<int>(AmbTrkType::kOrphan));
              }
              if (ids.size() == 1 && trkCollId == ids[0]) {
                qaregistry.fill(HIST("TrackToColl/hAmbTrackType"), static_cast<int>(AmbTrkType::kNonAmb));
              }
            }
            if (gConf.cfgRemoveOrphanTracks && ids.empty()) {
              continue;
            }
            if (gConf.cfgRemoveTrivialAssoc) {
              if (ids.empty() || (ids.size() == 1 && trkCollId == ids[0])) {
                qaregistry.fill(HIST("TrackToColl/histNonAmbTrackNumColls"), ids.size());
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
            if (ids.size() > 0) {
              if (ids.size() == 1) {
                if (trkCollId != ids[0]) {
                  qaregistry.fill(HIST("TrackToColl/hAmbTrackType"), static_cast<int>(AmbTrkType::kAmb));
                } else {
                  qaregistry.fill(HIST("TrackToColl/hAmbTrackType"), static_cast<int>(AmbTrkType::kNonAmbSame));
                }
              } else {
                qaregistry.fill(HIST("TrackToColl/hAmbTrackType"), static_cast<int>(AmbTrkType::kAmbGt1));
              }
            }

            bool isAmbiguous = (ids.size() > 1);
            if (isAmbiguous) {
              qaregistry.fill(HIST("TrackToColl/histAmbTrackNumColls"), ids.size());
              std::vector<float> ambVtxZ{};
              for (const auto& collIdx : ids) {
                const auto& ambColl = collisions.rawIteratorAt(collIdx);
                ambVtxZ.push_back(ambColl.posZ());
              }
              if (!ambVtxZ.empty()) {
                qaregistry.fill(HIST("TrackToColl/histAmbTrackZvtxRMS"), computeRMS(ambVtxZ));
              }
            }

            float deltaX = -999.f;
            float deltaY = -999.f;
            float deltaZ = -999.f;
            if (index) {
              const auto& collision = track.template collision_as<CollisionsWithMCLabels>();
              const auto& mcCollision = particle.template mcCollision_as<aod::McCollisions>();
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
        qaregistry.fill(HIST("TrackToColl/histFracGoodTracks"), frac);
        float fracFake = (nTrk > 0) ? static_cast<float>(nFakeTrk) / nTrk : -1.;
        qaregistry.fill(HIST("TrackToColl/histFracTracksFakeMcColl"), fracFake);
      }
    }
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processAssocMC, "Process collision-association information, requires extra table from TrackToCollisionAssociation task (fillTableOfCollIdsPerTrack=true)", false);

  template <typename B>
  void processReAssocMC(CollisionsWithMCLabels const& collisions,
                        B const& besttracks,
                        FiltMcMftTracks const& /*tracks*/,
                        aod::McCollisions const& /*mcCollisions*/,
                        aod::McParticles const& /*particles*/
  )
  {
    const auto& nRecoColls = collisions.size();
    registry.fill(HIST("ReAssocMC/hNReAssocRecoColls"), 1.f, nRecoColls);

    if (nRecoColls > CintZero) {

      mapVtxXrec.clear();
      mapVtxYrec.clear();
      mapVtxZrec.clear();
      mapMcCollIdPerRecColl.clear();
      mapVtxXrec.reserve(collisions.size());
      mapVtxYrec.reserve(collisions.size());
      mapVtxZrec.reserve(collisions.size());
      mapMcCollIdPerRecColl.reserve(collisions.size());

      auto maxNcontributors = -1;
      auto bestCollIndex = -1;
      for (auto const& collision : collisions) {
        if (maxNcontributors < collision.numContrib()) {
          maxNcontributors = collision.numContrib();
          bestCollIndex = collision.globalIndex();
          mapVtxXrec.emplace(collision.globalIndex(), collision.posX());
          mapVtxYrec.emplace(collision.globalIndex(), collision.posY());
          mapVtxZrec.emplace(collision.globalIndex(), collision.posZ());
          mapMcCollIdPerRecColl.emplace(collision.globalIndex(), collision.mcCollisionId());
        }
      }

      int nNoMC{0};
      for (const auto& collision : collisions) {
        registry.fill(HIST("ReAssocMC/hReAssocMCEventStatus"), static_cast<int>(ReAssocMCEventStatus::kEvtReAsAll));
        if (!isGoodEvent<true>(collision)) {
          continue;
        }
        registry.fill(HIST("ReAssocMC/hReAssocMCEventStatus"), static_cast<int>(ReAssocMCEventStatus::kEvtReAsSelected));
        if (!collision.has_mcCollision()) {
          continue;
        }
        registry.fill(HIST("ReAssocMC/hReAssocMCEventStatus"), static_cast<int>(ReAssocMCEventStatus::kEvtReAsHasMcColl));

        int64_t recCollId = collision.globalIndex();
        auto itMC = mapMcCollIdPerRecColl.find(recCollId);
        if (itMC == mapMcCollIdPerRecColl.end()) {
          nNoMC++;
          continue;
        }

        auto mcColl = collision.template mcCollision_as<aod::McCollisions>();
        // Select collisions with the largest number of contributors
        if (gConf.cfgRemoveSplitVertex && (bestCollIndex != collision.globalIndex())) {
          continue;
        }
        registry.fill(HIST("ReAssocMC/hReAssocMCEventStatus"), static_cast<int>(ReAssocMCEventStatus::kEvtReAsSplitVtxRemoved));
        if (eventCuts.useZVtxCutMC && (std::abs(mcColl.posZ()) >= eventCuts.maxZvtx)) {
          continue;
        }
        registry.fill(HIST("ReAssocMC/hReAssocMCEventStatus"), static_cast<int>(ReAssocMCEventStatus::kEvtReAsZVtxCutMC));

        for (auto const& atrack : besttracks) {
          registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kTrkReAssocAll));
          if (!isBestTrackSelected<true>(atrack)) {
            continue;
          }
          registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kTrkBestSel));

          auto itrack = atrack.template mfttrack_as<FiltMcMftTracks>();

          if (!isTrackSelected<true>(itrack)) {
            continue;
          }
          registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kTrkSel));
          float phi = itrack.phi();
          o2::math_utils::bringTo02Pi(phi);
          if (phi < Czero || TwoPI < phi) {
            continue;
          }
          if (!itrack.has_collision()) {
            continue;
          }
          registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kTrkHasColl));
          if (gConf.cfgRemoveReassigned) {
            if (itrack.collisionId() != atrack.bestCollisionId()) {
              continue;
            }
            registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kTrkReassignedRemoved));
          }

          if (itrack.collisionId() >= 0 && itrack.has_mcParticle() && itrack.mcMask() == 0) {

            auto particle = itrack.template mcParticle_as<aod::McParticles>();
            if (!isChrgParticle(particle.pdgCode())) {
              continue;
            }
            if (gConf.cfgUseParticleSel && !isParticleSelected(particle)) {
              continue;
            }
            auto collision = itrack.template collision_as<CollisionsWithMCLabels>();

            float deltaX = -999.f;
            float deltaY = -999.f;
            float deltaZ = -999.f;

            if (mapVtxXrec.find(atrack.bestCollisionId()) == mapVtxXrec.end()) {
              continue;
            }
            if (mapVtxYrec.find(atrack.bestCollisionId()) == mapVtxYrec.end()) {
              continue;
            }
            if (mapVtxZrec.find(atrack.bestCollisionId()) == mapVtxZrec.end()) {
              continue;
            }
            if (mapMcCollIdPerRecColl.find(atrack.bestCollisionId()) == mapMcCollIdPerRecColl.end()) {
              continue;
            }
            const float vtxXbest = mapVtxXrec.find(atrack.bestCollisionId())->second;
            const float vtxYbest = mapVtxYrec.find(atrack.bestCollisionId())->second;
            const float vtxZbest = mapVtxZrec.find(atrack.bestCollisionId())->second;
            // LOGP(info, "\t ---> \t .... \t vtxZrec: {} - collision.posZ(): {}", vtxZrec, collision.posZ());
            const float mcCollIdRec = mapMcCollIdPerRecColl.find(atrack.bestCollisionId())->second;
            // LOGP(info, "\t ---> \t .... \t mcCollIdRec: {} - bestMCCol: {}", mcCollIdRec, bestMCCol);

            auto vtxXtruth = atrack.mcParticle().mcCollision().posX();
            auto vtxYtruth = atrack.mcParticle().mcCollision().posY();
            auto vtxZtruth = atrack.mcParticle().mcCollision().posZ();

            deltaX = vtxXbest - vtxXtruth;
            deltaY = vtxYbest - vtxYtruth;
            deltaZ = vtxZbest - vtxZtruth;

            const auto dcaXtruth(particle.vx() - particle.mcCollision().posX());
            const auto dcaYtruth(particle.vy() - particle.mcCollision().posY());
            const auto dcaZtruth(particle.vz() - particle.mcCollision().posZ());
            auto dcaXYtruth = std::sqrt(dcaXtruth * dcaXtruth + dcaYtruth * dcaYtruth);

            // check good rec vertex of corresponding mc coll is available in compatible rec coll
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
                      break;
                    }
                  }
                }
              }
            }

            registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kTrkHasMcPart));

            if (ids.size() > 0) {
              if (ids.size() == 1) {
                if (itrack.collisionId() == ids[0]) { // non ambiguous
                  registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kTrkNonAmbAll));
                  hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbAll)]->Fill(particle.pt(), particle.eta(), deltaX, deltaY, deltaZ);
                  hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbAll)]->Fill(particle.pt(), particle.eta(), 0., 0., dcaXYtruth, dcaZtruth);
                  if (collision.mcCollisionId() == particle.mcCollisionId()) {
                    registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kTrkNonAmbGood));
                    hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbGood)]->Fill(particle.pt(), particle.eta(), deltaX, deltaY, deltaZ);
                    hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbGood)]->Fill(particle.pt(), particle.eta(), 0., 0., dcaXYtruth, dcaZtruth);
                  } else {
                    registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kTrkNonAmbBad));
                    hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbBad)]->Fill(particle.pt(), particle.eta(), deltaX, deltaY, deltaZ);
                    hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbBad)]->Fill(particle.pt(), particle.eta(), 0., 0., dcaXYtruth, dcaZtruth);
                  }
                } else if (itrack.collisionId() != ids[0]) { // ambiguous extra
                  registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kTrkAmbAll));
                  hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkAmbAll)]->Fill(particle.pt(), particle.eta(), deltaX, deltaY, deltaZ);
                  hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkAmbAll)]->Fill(particle.pt(), particle.eta(), atrack.bestDCAXY(), atrack.bestDCAZ(), dcaXYtruth, dcaZtruth);

                  if (collision.mcCollisionId() == particle.mcCollisionId()) {
                    registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kTrkAmbGood));
                    hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkAmbGood)]->Fill(particle.pt(), particle.eta(), deltaX, deltaY, deltaZ);
                    hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkAmbGood)]->Fill(particle.pt(), particle.eta(), atrack.bestDCAXY(), atrack.bestDCAZ(), dcaXYtruth, dcaZtruth);
                  } else {
                    registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kTrkAmbBad));
                    hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkAmbBad)]->Fill(particle.pt(), particle.eta(), deltaX, deltaY, deltaZ);
                    hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkAmbBad)]->Fill(particle.pt(), particle.eta(), atrack.bestDCAXY(), atrack.bestDCAZ(), dcaXYtruth, dcaZtruth);
                  }
                } else { // non ambiguous (extra)
                  registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kTrkNonAmbAllE));
                  hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbAllE)]->Fill(particle.pt(), particle.eta(), deltaX, deltaY, deltaZ);
                  hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbAllE)]->Fill(particle.pt(), particle.eta(), 0., 0., dcaXYtruth, dcaZtruth);

                  if (collision.mcCollisionId() == particle.mcCollisionId()) {
                    registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kTrkNonAmbGoodE));
                    hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbGoodE)]->Fill(particle.pt(), particle.eta(), deltaX, deltaY, deltaZ);
                    hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbGoodE)]->Fill(particle.pt(), particle.eta(), 0., 0., dcaXYtruth, dcaZtruth);
                  } else {
                    registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kTrkNonAmbBadE));
                    hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbBadE)]->Fill(particle.pt(), particle.eta(), deltaX, deltaY, deltaZ);
                    hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kTrkNonAmbBadE)]->Fill(particle.pt(), particle.eta(), 0., 0., dcaXYtruth, dcaZtruth);
                  }
                }
              } else { // ambiguous

                if (itrack.collisionId() == atrack.bestCollisionId()) { // associated
                  registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kAssoc));
                  hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kAssoc)]->Fill(particle.pt(), particle.eta(), deltaX, deltaY, deltaZ);
                  hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kAssoc)]->Fill(particle.pt(), particle.eta(), atrack.bestDCAXY(), atrack.bestDCAZ(), dcaXYtruth, dcaZtruth);

                  if (collision.has_mcCollision() && mcCollIdRec == particle.mcCollisionId()) { // good coll
                    registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kAssocGood));
                    hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kAssocGood)]->Fill(particle.pt(), particle.eta(), deltaX, deltaY, deltaZ);
                    hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kAssocGood)]->Fill(particle.pt(), particle.eta(), atrack.bestDCAXY(), atrack.bestDCAZ(), dcaXYtruth, dcaZtruth);

                    if (isInCoColl) { // coll vertex is among compatible colls
                      registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kAssocGoodIsCompTrue));
                      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kAssocGoodIsCompTrue)]->Fill(particle.pt(), particle.eta(), deltaX, deltaY, deltaZ);
                      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kAssocGoodIsCompTrue)]->Fill(particle.pt(), particle.eta(), atrack.bestDCAXY(), atrack.bestDCAZ(), dcaXYtruth, dcaZtruth);
                    } else {
                      registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kAssocGoodIsCompFalse));
                      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kAssocGoodIsCompFalse)]->Fill(particle.pt(), particle.eta(), deltaX, deltaY, deltaZ);
                      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kAssocGoodIsCompFalse)]->Fill(particle.pt(), particle.eta(), atrack.bestDCAXY(), atrack.bestDCAZ(), dcaXYtruth, dcaZtruth);
                    }
                  } else {
                    registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kAssocBad));
                    hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kAssocBad)]->Fill(particle.pt(), particle.eta(), deltaX, deltaY, deltaZ);
                    hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kAssocBad)]->Fill(particle.pt(), particle.eta(), atrack.bestDCAXY(), atrack.bestDCAZ(), dcaXYtruth, dcaZtruth);

                    if (isInCoColl) { // coll vertex is among compatible colls
                      registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kAssocBadIsCompTrue));
                      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kAssocBadIsCompTrue)]->Fill(particle.pt(), particle.eta(), deltaX, deltaY, deltaZ);
                      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kAssocBadIsCompTrue)]->Fill(particle.pt(), particle.eta(), atrack.bestDCAXY(), atrack.bestDCAZ(), dcaXYtruth, dcaZtruth);
                    } else {
                      registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kAssocBadIsCompFalse));
                      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kAssocBadIsCompFalse)]->Fill(particle.pt(), particle.eta(), deltaX, deltaY, deltaZ);
                      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kAssocBadIsCompFalse)]->Fill(particle.pt(), particle.eta(), atrack.bestDCAXY(), atrack.bestDCAZ(), dcaXYtruth, dcaZtruth);
                    }
                  }
                } else { // reassociated
                  registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kReAssoc));
                  if (collision.has_mcCollision() && mcCollIdRec == particle.mcCollisionId()) { // good coll
                    registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kReAssocGood));
                    hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kReAssocGood)]->Fill(particle.pt(), particle.eta(), deltaX, deltaY, deltaZ);
                    hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kReAssocGood)]->Fill(particle.pt(), particle.eta(), atrack.bestDCAXY(), atrack.bestDCAZ(), dcaXYtruth, dcaZtruth);

                    if (isInCoColl) { // coll vertex is among compatible colls
                      registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kReAssocGoodIsCompTrue));
                      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kReAssocGoodIsCompTrue)]->Fill(particle.pt(), particle.eta(), deltaX, deltaY, deltaZ);
                      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kReAssocGoodIsCompTrue)]->Fill(particle.pt(), particle.eta(), atrack.bestDCAXY(), atrack.bestDCAZ(), dcaXYtruth, dcaZtruth);
                    } else {
                      registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kReAssocGoodIsCompFalse));
                      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kReAssocGoodIsCompFalse)]->Fill(particle.pt(), particle.eta(), deltaX, deltaY, deltaZ);
                      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kReAssocGoodIsCompFalse)]->Fill(particle.pt(), particle.eta(), atrack.bestDCAXY(), atrack.bestDCAZ(), dcaXYtruth, dcaZtruth);
                    }
                  } else {
                    registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kReAssocBad));
                    hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kReAssocBad)]->Fill(particle.pt(), particle.eta(), deltaX, deltaY, deltaZ);
                    hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kReAssocBad)]->Fill(particle.pt(), particle.eta(), atrack.bestDCAXY(), atrack.bestDCAZ(), dcaXYtruth, dcaZtruth);

                    if (isInCoColl) { // coll vertex is among compatible colls
                      registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kReAssocBadIsCompTrue));
                      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kReAssocBadIsCompTrue)]->Fill(particle.pt(), particle.eta(), deltaX, deltaY, deltaZ);
                      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kReAssocBadIsCompTrue)]->Fill(particle.pt(), particle.eta(), atrack.bestDCAXY(), atrack.bestDCAZ(), dcaXYtruth, dcaZtruth);
                    } else {
                      registry.fill(HIST("ReAssocMC/hReAssocMCTrackStatus"), static_cast<int>(ReAssocMCTrackStatus::kReAssocBadIsCompFalse));
                      hReAssocVtxRes[static_cast<int>(HistStatusReAssocVtx::kReAssocBadIsCompFalse)]->Fill(particle.pt(), particle.eta(), deltaX, deltaY, deltaZ);
                      hReAssocDCA[static_cast<int>(HistStatusReAssocVtx::kReAssocBadIsCompFalse)]->Fill(particle.pt(), particle.eta(), atrack.bestDCAXY(), atrack.bestDCAZ(), dcaXYtruth, dcaZtruth);
                    }
                  }
                }
              }
            }
          }
        }
      }
      mapVtxXrec.clear();
      mapVtxYrec.clear();
      mapVtxZrec.clear();
      mapMcCollIdPerRecColl.clear();
      LOG(info) << "No MC: " << nNoMC;
    }
  }

  void processReAssoc3dMC(CollisionsWithMCLabels const& collisions,
                          BestTracks3dWCollsMC const& besttracks,
                          FiltMcMftTracks const& tracks,
                          aod::McCollisions const& mcCollisions,
                          aod::McParticles const& particles)
  {
    processReAssocMC<BestTracks3dWCollsMC>(collisions, besttracks, tracks, mcCollisions, particles);
  }
  PROCESS_SWITCH(DndetaMFTPbPb, processReAssoc3dMC, "Process re-association information based on BestCollisionsFwd3d table", false);

  Preslice<MftTracksWCollsMC> compCollPerCol = o2::aod::fwdtrack::collisionId;

  /// @brief process function to check reassociation efficiency
  void processReassocEfficiency(CollsMCExtraMult::iterator const& mcCollision,
                                soa::SmallGroups<CollsGenCentFT0C> const& collisions,
                                MftTracksWCollsMC const& tracks,
                                aod::McParticles const& particles,
                                ExtBCs const& bcs)
  {
    if (bcs.size() == 0) {
      return;
    }
    auto bc = bcs.begin();
    initCCDB(bc);

    // At least one generated primary in MFT acceptance + TVX triggered collisions
    if (gConf.cfgUseInelgt0wMFT && !isInelGt0wMft(particles)) {
      return;
    }
    if (eventCuts.useZVtxCutMC && (std::abs(mcCollision.posZ()) >= eventCuts.maxZvtx)) {
      return;
    }

    const auto& nRecoColls = collisions.size();
    registry.fill(HIST("Events/hNReEffColls"), 1.f, nRecoColls);

    // Generated evets with >= 1 reco collisions
    if (nRecoColls > CintZero) {

      auto maxNcontributors = -1;
      auto bestCollIndex = -1;
      auto crec = -1;
      for (auto const& collision : collisions) {
        if (!isGoodEvent<false>(collision)) {
          continue;
        }
        if (maxNcontributors < collision.numContrib()) {
          maxNcontributors = collision.numContrib();
          bestCollIndex = collision.globalIndex();
          crec = getRecoCent(collision);
        }
      }

      registry.fill(HIST("AmbTracks/hVtxzMCrec"), mcCollision.posZ());

      for (const auto& collision : collisions) {
        if (!isGoodEvent<false>(collision)) {
          continue;
        }
        // Select collisions with the largest number of contributors
        if (gConf.cfgRemoveSplitVertex && (bestCollIndex != collision.globalIndex())) {
          continue;
        }
        registry.fill(HIST("Events/Centrality/ReEffStatus"), 1, crec);
        if (!collision.has_mcCollision()) {
          continue;
        }

        std::array<double, 3> dcaInfOrig;
        std::array<double, 2> dcaInfo;
        double bestDCA[2];
        auto trkPerColl = tracks.sliceBy(compCollPerCol, collision.globalIndex());

        for (auto const& track : trkPerColl) {
          dcaInfOrig[0] = 999.f; // original DCAx from propagation
          dcaInfOrig[1] = 999.f; // original DCAy from propagation
          dcaInfOrig[2] = 999.f; // original DCAz from propagation
          dcaInfo[0] = 999.f;    // calcualted DCAxy
          dcaInfo[1] = 999.f;    // calculated DCAz - same as original
          bestDCA[0] = 999.f;    // minimal DCAxy
          bestDCA[1] = 999.f;    // minimal DCAz

          if (!isTrackSelected<false>(track)) {
            continue;
          }

          auto bestCol = track.has_collision() ? track.collisionId() : -1;
          uint index = uint(track.collisionId() >= 0);
          auto ids = track.compatibleCollIds();
          bool isAmbiguous = (ids.size() > 1);

          if (gConf.cfgRemoveReassigned) {
            if (ids.empty() || (ids.size() == 1 && bestCol == ids[0])) {
              continue;
            }
          }
          if (gConf.cfgRemoveOrphanTracks && ids.empty()) {
            continue;
          }
          if (!track.has_mcParticle()) {
            continue;
          }
          auto particle = track.mcParticle_as<aod::McParticles>();
          if (!isChrgParticle(particle.pdgCode())) {
            continue;
          }
          if (particle.eta() <= trackCuts.minEta || particle.eta() >= trackCuts.maxEta) {
            continue;
          }
          if (gConf.cfgUseParticleSel && !isParticleSelected(particle)) {
            continue;
          }
          int bestMCCol = -1;
          o2::track::TrackParCovFwd trackPar = o2::aod::fwdtrackutils::getTrackParCovFwdShift(track, mZShift);

          if (index) {
            auto mcCollision = particle.mcCollision_as<CollsMCExtraMult>();
            if (eventCuts.useZDiffCut) {
              if (std::abs(collision.posZ() - mcCollision.posZ()) > eventCuts.maxZvtxDiff) {
                continue;
              }
            }

            if (isAmbiguous) {
              for (const auto& collIdx : track.compatibleCollIds()) {
                auto ambColl = collisions.rawIteratorAt(collIdx);
                if (!ambColl.has_mcCollision()) {
                  continue;
                }

                trackPar.propagateToDCAhelix(bZ, {ambColl.posX(), ambColl.posY(), ambColl.posZ()}, dcaInfOrig);

                if (gConf.cfgDCAtype == 0) {
                  dcaInfo[0] = dcaInfOrig[0];
                } else if (gConf.cfgDCAtype == 1) {
                  dcaInfo[0] = dcaInfOrig[1];
                } else if (gConf.cfgDCAtype == 2) {
                  dcaInfo[0] = std::sqrt(dcaInfOrig[0] * dcaInfOrig[0] + dcaInfOrig[1] * dcaInfOrig[1]);
                }
                dcaInfo[1] = dcaInfOrig[2];

                if ((std::abs(dcaInfo[0]) < std::abs(bestDCA[0])) && (std::abs(dcaInfo[1]) < std::abs(bestDCA[1]))) {
                  bestCol = ambColl.globalIndex();
                  bestMCCol = ambColl.mcCollisionId();
                  bestDCA[0] = dcaInfo[0];
                  bestDCA[1] = dcaInfo[1];
                }

                // LOGP(info, "-> track {}: {}", track.globalIndex(), dcaInfo[0]);
                registry.fill(HIST("AmbTracks/DCAXY"), dcaInfo[0]);
                registry.fill(HIST("AmbTracks/DCAZ"), dcaInfo[1]);
              }

              registry.fill(HIST("AmbTracks/DCAXYBest"), bestDCA[0]);
              registry.fill(HIST("AmbTracks/DCAZBest"), bestDCA[1]);
              registry.fill(HIST("AmbTracks/Centrality/THnDCAxyBestGen"), particle.pt(), particle.eta(), mcCollision.posZ(), bestDCA[0], bestDCA[1], crec);
              registry.fill(HIST("AmbTracks/BestGenDxyz"), collision.posX() - mcCollision.posX(), collision.posY() - mcCollision.posY(), collision.posZ() - mcCollision.posZ(), bestDCA[0], bestDCA[1]);

              if (particle.isPhysicalPrimary()) {
                registry.fill(HIST("AmbTracks/DCAXYBestPrim"), bestDCA[0]);
                registry.fill(HIST("AmbTracks/DCAZBestPrim"), bestDCA[1]);
                registry.fill(HIST("AmbTracks/Centrality/THnDCAxyBestGenPrim"), particle.pt(), particle.eta(), mcCollision.posZ(), bestDCA[0], bestDCA[1], crec);
                registry.fill(HIST("AmbTracks/BestPrimDxyz"), collision.posX() - mcCollision.posX(), collision.posY() - mcCollision.posY(), collision.posZ() - mcCollision.posZ(), bestDCA[0], bestDCA[1]);
              }

              auto mcCollID = particle.mcCollisionId();
              registry.fill(HIST("Events/Centrality/ReEffStatus"), 2, crec);

              if (!track.has_collision()) {
                registry.fill(HIST("Events/Centrality/ReEffStatus"), 4, crec);
                if (bestMCCol == mcCollID) // correctly assigned to bestCol
                {
                  registry.fill(HIST("Events/Centrality/ReEffStatus"), 6, crec);
                  registry.fill(HIST("AmbTracks/DCAXYBestTrue"), bestDCA[0]);
                  registry.fill(HIST("AmbTracks/DCAZBestTrue"), bestDCA[1]);
                  registry.fill(HIST("AmbTracks/Centrality/THnDCAxyBestTrue"), particle.pt(), particle.eta(), mcCollision.posZ(), bestDCA[0], bestDCA[1], crec);
                  registry.fill(HIST("AmbTracks/BestTrueDxyz"), collision.posX() - mcCollision.posX(), collision.posY() - mcCollision.posY(), collision.posZ() - mcCollision.posZ(), bestDCA[0], bestDCA[1]);
                } else {
                  registry.fill(HIST("AmbTracks/DCAXYBestFalse"), bestDCA[0]);
                  registry.fill(HIST("AmbTracks/DCAZBestFalse"), bestDCA[1]);
                  registry.fill(HIST("AmbTracks/Centrality/THnDCAxyBestFalse"), particle.pt(), particle.eta(), mcCollision.posZ(), bestDCA[0], bestDCA[1], crec);
                  registry.fill(HIST("AmbTracks/BestFalseDxyz"), collision.posX() - mcCollision.posX(), collision.posY() - mcCollision.posY(), collision.posZ() - mcCollision.posZ(), bestDCA[0], bestDCA[1]);
                }
              } else if (track.collisionId() != bestCol) { // reassgined
                auto collOrig = collisions.rawIteratorAt(track.collisionId());
                registry.fill(HIST("Events/Centrality/ReEffStatus"), 3, crec);
                if (bestMCCol == mcCollID) // correctly reassigned
                {
                  registry.fill(HIST("Events/Centrality/ReEffStatus"), 6, crec);
                  registry.fill(HIST("AmbTracks/DCAXYBestTrue"), bestDCA[0]);
                  registry.fill(HIST("AmbTracks/DCAZBestTrue"), bestDCA[1]);
                  registry.fill(HIST("AmbTracks/Centrality/THnDCAxyBestTrue"), particle.pt(), particle.eta(), mcCollision.posZ(), bestDCA[0], bestDCA[1], crec);
                  registry.fill(HIST("AmbTracks/BestTrueDxyz"), collision.posX() - mcCollision.posX(), collision.posY() - mcCollision.posY(), collision.posZ() - mcCollision.posZ(), bestDCA[0], bestDCA[1]);
                } else { // uncorrectly reassigned
                  registry.fill(HIST("AmbTracks/DCAXYBestFalse"), bestDCA[0]);
                  registry.fill(HIST("AmbTracks/DCAZBestFalse"), bestDCA[1]);
                  registry.fill(HIST("AmbTracks/Centrality/THnDCAxyBestFalse"), particle.pt(), particle.eta(), mcCollision.posZ(), bestDCA[0], bestDCA[1], crec);
                  registry.fill(HIST("AmbTracks/BestFalseDxyz"), collision.posX() - mcCollision.posX(), collision.posY() - mcCollision.posY(), collision.posZ() - mcCollision.posZ(), bestDCA[0], bestDCA[1]);
                }
                if (collOrig.mcCollisionId() == mcCollID) { // initially correctly assigned
                  registry.fill(HIST("Events/Centrality/ReEffStatus"), 5, crec);
                  registry.fill(HIST("AmbTracks/DCAXYBestTrueOrigReAssoc"), bestDCA[0]);
                  registry.fill(HIST("AmbTracks/DCAZBestTrueOrigReAssoc"), bestDCA[1]);
                  registry.fill(HIST("AmbTracks/Centrality/THnDCAxyBestTrueOrigReAssoc"), particle.pt(), particle.eta(), mcCollision.posZ(), bestDCA[0], bestDCA[1], crec);
                  registry.fill(HIST("AmbTracks/BestTrueOrigReAssocDxyz"), collision.posX() - mcCollision.posX(), collision.posY() - mcCollision.posY(), collision.posZ() - mcCollision.posZ(), bestDCA[0], bestDCA[1]);
                }
              } else { // not reassigned - the track has a collision and track.collisionId() == bestCol
                if (track.collisionId() != bestCol) {
                  registry.fill(HIST("Events/Centrality/ReEffStatus"), 7, crec);
                  // LOGP(info, "-> track id {}: bestCollid {}", track.collisionId(), bestCol);
                }
                registry.fill(HIST("Events/Centrality/ReEffStatus"), 8, crec);
                if (bestMCCol == mcCollID) // correctly assigned
                {
                  registry.fill(HIST("Events/Centrality/ReEffStatus"), 9, crec);
                  registry.fill(HIST("AmbTracks/DCAXYBestTrueOrigAssoc"), bestDCA[0]);
                  registry.fill(HIST("AmbTracks/DCAZBestTrueOrigAssoc"), bestDCA[1]);
                  registry.fill(HIST("AmbTracks/Centrality/THnDCAxyBestTrueOrigAssoc"), particle.pt(), particle.eta(), mcCollision.posZ(), bestDCA[0], bestDCA[1], crec);
                  registry.fill(HIST("AmbTracks/BestTrueOrigAssocDxyz"), collision.posX() - mcCollision.posX(), collision.posY() - mcCollision.posY(), collision.posZ() - mcCollision.posZ(), bestDCA[0], bestDCA[1]);
                }
              }
            }
          }
        }
      }
    }
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processReassocEfficiency, "Process collision-reassociation efficiency based on track propagation DCA 3D (in FT0C centrality bins)", false);

  Preslice<MftTracksWCollsMC> mftTrkCompCollperCol = o2::aod::fwdtrack::collisionId;

  /// @brief process template function for DCA MC checks
  template <typename C, typename MC>
  void processSecondariesMC(typename soa::Join<C, aod::McCollisionLabels> const& collisions,
                            MftTracksWCollsMC const& tracks,
                            MC const& mcCollisions,
                            aod::McParticles const& particles,
                            ExtBCs const& bcs)
  {
    registry.fill(HIST("Events/hNGenRecColls"), 1.f, collisions.size());
    registry.fill(HIST("Events/hNGenRecColls"), 2.f, mcCollisions.size());

    if (bcs.size() == 0) {
      return;
    }
    auto bc = bcs.begin();
    initCCDB(bc);

    // Generated evets with >= 1 reco collisions
    const auto& nRecoColls = collisions.size();
    if (nRecoColls > CintZero) {
      auto maxNcontributors = -1;
      auto bestCollIndex = -1;
      auto crec = -1;
      for (auto const& collision : collisions) {
        if (!isGoodEvent<false>(collision)) {
          continue;
        }
        if (maxNcontributors < collision.numContrib()) {
          maxNcontributors = collision.numContrib();
          bestCollIndex = collision.globalIndex();
          if constexpr (has_reco_cent<C>) {
            crec = getRecoCent(collision);
          }
        }
      }

      for (const auto& collision : collisions) {
        if constexpr (has_reco_cent<C>) {
          qaregistry.fill(HIST("Events/Centrality/hTrackToCollEvtType"), static_cast<int>(TrackToCollEvtType::kAllRecColl), crec);
        } else {
          qaregistry.fill(HIST("Events/hTrackToCollEvtType"), static_cast<int>(TrackToCollEvtType::kAllRecColl));
        }
        // At least one generated primary in MFT acceptance + TVX triggered collisions
        if (gConf.cfgUseInelgt0wMFT && !isInelGt0wMft(particles)) {
          return;
        }
        if constexpr (has_reco_cent<C>) {
          qaregistry.fill(HIST("Events/Centrality/hTrackToCollEvtType"), static_cast<int>(TrackToCollEvtType::kIsInelGt0wMft), crec);
        } else {
          qaregistry.fill(HIST("Events/hTrackToCollEvtType"), static_cast<int>(TrackToCollEvtType::kIsInelGt0wMft));
        }
        if (!isGoodEvent<true>(collision)) {
          continue;
        }
        if constexpr (has_reco_cent<C>) {
          qaregistry.fill(HIST("Events/Centrality/hTrackToCollEvtType"), static_cast<int>(TrackToCollEvtType::kEvtSel), crec);
        } else {
          qaregistry.fill(HIST("Events/hTrackToCollEvtType"), static_cast<int>(TrackToCollEvtType::kEvtSel));
        }
        // Select collisions with the largest number of contributors
        if (gConf.cfgRemoveSplitVertex && (bestCollIndex != collision.globalIndex())) {
          continue;
        }
        if constexpr (has_reco_cent<C>) {
          qaregistry.fill(HIST("Events/Centrality/hTrackToCollEvtType"), static_cast<int>(TrackToCollEvtType::kBestCollIdx), crec);
        } else {
          qaregistry.fill(HIST("Events/hTrackToCollEvtType"), static_cast<int>(TrackToCollEvtType::kBestCollIdx));
        }
        if (!collision.has_mcCollision()) {
          continue;
        }
        if constexpr (has_reco_cent<C>) {
          qaregistry.fill(HIST("Events/Centrality/hTrackToCollEvtType"), static_cast<int>(TrackToCollEvtType::kIsMcColl), crec);
        } else {
          qaregistry.fill(HIST("Events/hTrackToCollEvtType"), static_cast<int>(TrackToCollEvtType::kIsMcColl));
        }

        auto trkPerColl = tracks.sliceBy(mftTrkCompCollperCol, collision.globalIndex());
        for (auto const& track : trkPerColl) {
          if (!isTrackSelected<false>(track)) {
            continue;
          }
          if (!track.has_collision()) {
            continue;
          }
          auto trkCollId = track.has_collision() ? track.collisionId() : -1;
          auto ids = track.compatibleCollIds();
          bool isAmbiguous = (ids.size() > 1);

          if constexpr (has_reco_cent<C>) {
            qaregistry.fill(HIST("TrkCompColls/Centrality/hAmbTrackType"), static_cast<int>(AmbTrkType::kAll), crec);
          } else {
            qaregistry.fill(HIST("TrkCompColls/hAmbTrackType"), static_cast<int>(AmbTrkType::kAll));
          }
          if (ids.empty()) {
            if constexpr (has_reco_cent<C>) {
              qaregistry.fill(HIST("TrkCompColls/Centrality/hAmbTrackType"), static_cast<int>(AmbTrkType::kOrphan), crec);
            } else {
              qaregistry.fill(HIST("TrkCompColls/hAmbTrackType"), static_cast<int>(AmbTrkType::kOrphan));
            }
          }
          if (ids.size() > 0) {
            if (ids.size() == 1) {
              if (trkCollId == ids[0]) {
                if constexpr (has_reco_cent<C>) {
                  qaregistry.fill(HIST("TrkCompColls/Centrality/hAmbTrackType"), static_cast<int>(AmbTrkType::kNonAmb), crec);
                } else {
                  qaregistry.fill(HIST("TrkCompColls/hAmbTrackType"), static_cast<int>(AmbTrkType::kNonAmb));
                }
              } else if (trkCollId != ids[0]) {
                if constexpr (has_reco_cent<C>) {
                  qaregistry.fill(HIST("TrkCompColls/Centrality/hAmbTrackType"), static_cast<int>(AmbTrkType::kAmb), crec);
                } else {
                  qaregistry.fill(HIST("TrkCompColls/hAmbTrackType"), static_cast<int>(AmbTrkType::kAmb));
                }
              } else {
                if constexpr (has_reco_cent<C>) {
                  qaregistry.fill(HIST("TrkCompColls/Centrality/hAmbTrackType"), static_cast<int>(AmbTrkType::kNonAmbSame), crec);
                } else {
                  qaregistry.fill(HIST("TrkCompColls/hAmbTrackType"), static_cast<int>(AmbTrkType::kNonAmbSame));
                }
              }
            } else {
              if constexpr (has_reco_cent<C>) {
                qaregistry.fill(HIST("TrkCompColls/Centrality/hAmbTrackType"), static_cast<int>(AmbTrkType::kAmbGt1), crec);
              } else {
                qaregistry.fill(HIST("TrkCompColls/hAmbTrackType"), static_cast<int>(AmbTrkType::kAmbGt1));
              }
            }
          } else {
            if constexpr (has_reco_cent<C>) {
              qaregistry.fill(HIST("TrkCompColls/Centrality/hAmbTrackType"), static_cast<int>(AmbTrkType::kOrphanNull), crec);
            } else {
              qaregistry.fill(HIST("TrkCompColls/hAmbTrackType"), static_cast<int>(AmbTrkType::kOrphanNull));
            }
          }

          if (gConf.cfgRemoveTrivialAssoc) {
            if (ids.empty() || (ids.size() == 1 && trkCollId == ids[0])) {
              qaregistry.fill(HIST("TrkCompColls/Centrality/histNonAmbTrackNumColls"), ids.size());
              continue;
            }
          }

          if (isAmbiguous) {
            if constexpr (has_reco_cent<C>) {
              registry.fill(HIST("Tracks/Centrality/THnRecAmb"), track.pt(), track.eta(), collision.posZ(), crec);
            } else {
              registry.fill(HIST("Tracks/THnRecAmb"), track.pt(), track.eta(), collision.posZ());
            }
          } else {
            if constexpr (has_reco_cent<C>) {
              registry.fill(HIST("Tracks/Centrality/THnRec"), track.pt(), track.eta(), collision.posZ(), crec);
            } else {
              registry.fill(HIST("Tracks/THnRec"), track.pt(), track.eta(), collision.posZ());
            }
            if (trkCollId == ids[0]) {
              if constexpr (has_reco_cent<C>) {
                registry.fill(HIST("Tracks/Centrality/THnRecNonAmb"), track.pt(), track.eta(), collision.posZ(), crec);
              } else {
                registry.fill(HIST("Tracks/THnRecNonAmb"), track.pt(), track.eta(), collision.posZ());
              }
            }
            if (trkCollId != ids[0]) {
              if constexpr (has_reco_cent<C>) {
                registry.fill(HIST("Tracks/Centrality/THnRecAmbRest"), track.pt(), track.eta(), collision.posZ(), crec);
              } else {
                registry.fill(HIST("Tracks/THnRecAmbRest"), track.pt(), track.eta(), collision.posZ());
              }
            }
          }

          // remove orphan tracks
          if (ids.empty()) {
            continue;
          }
          // remove fake tracks
          if (!track.has_mcParticle() || track.mcMask() != 0) {
            continue;
          }
          // assigned collision index
          uint index = uint(track.collisionId() >= 0);
          auto particle = track.template mcParticle_as<aod::McParticles>();
          if (!isChrgParticle(particle.pdgCode())) {
            continue;
          }
          if (particle.eta() <= trackCuts.minEta || particle.eta() >= trackCuts.maxEta) {
            continue;
          }
          if (gConf.cfgUseParticleSel && !isParticleSelected(particle)) {
            continue;
          }

          if (index) {
            auto mcCollision = particle.template mcCollision_as<aod::McCollisions>();
            if (eventCuts.useZVtxCutMC && (std::abs(mcCollision.posZ()) >= eventCuts.maxZvtx)) {
              continue;
            }
            if (eventCuts.useZDiffCut) {
              if (std::abs(collision.posZ() - mcCollision.posZ()) > eventCuts.maxZvtxDiff) {
                continue;
              }
            }

            if (collision.has_mcCollision() && collision.mcCollisionId() == particle.mcCollisionId()) {
              if (!particle.isPhysicalPrimary()) { // Secondaries (weak decays and material)
                if constexpr (has_reco_cent<C>) {
                  registry.fill(HIST("Tracks/Centrality/THnGenSec"), particle.pt(), particle.eta(), particle.mcCollision().posZ(), crec);
                } else {
                  registry.fill(HIST("Tracks/THnGenSec"), particle.pt(), particle.eta(), particle.mcCollision().posZ());
                }
                if (particle.getProcess() == TMCProcess::kPDecay) { // Particles from decay
                  if constexpr (has_reco_cent<C>) {
                    registry.fill(HIST("Tracks/Centrality/THnGenSecWeak"), particle.pt(), particle.eta(), particle.mcCollision().posZ(), crec);
                  } else {
                    registry.fill(HIST("Tracks/THnGenSecWeak"), particle.pt(), particle.eta(), particle.mcCollision().posZ());
                  }
                } else { // Particles from the material
                  if constexpr (has_reco_cent<C>) {
                    registry.fill(HIST("Tracks/Centrality/THnGenSecMat"), particle.pt(), particle.eta(), particle.mcCollision().posZ(), crec);
                  } else {
                    registry.fill(HIST("Tracks/THnGenSecMat"), particle.pt(), particle.eta(), particle.mcCollision().posZ());
                  }
                }
              } else { // Primaries
                if constexpr (has_reco_cent<C>) {
                  registry.fill(HIST("Tracks/Centrality/THnGenPrim"), particle.pt(), particle.eta(), particle.mcCollision().posZ(), crec);
                } else {
                  registry.fill(HIST("Tracks/THnGenPrim"), particle.pt(), particle.eta(), particle.mcCollision().posZ());
                }
              }
            } else {
              if (isAmbiguous) {
                for (const auto& collIdx : track.compatibleCollIds()) {
                  auto ambColl = collisions.rawIteratorAt(collIdx);
                  if (ambColl.has_mcCollision() && ambColl.mcCollisionId() == particle.mcCollisionId()) {
                    if (!particle.isPhysicalPrimary()) { // Secondaries (weak decays and material)
                      if constexpr (has_reco_cent<C>) {
                        registry.fill(HIST("Tracks/Centrality/THnGenSecAmb"), particle.pt(), particle.eta(), particle.mcCollision().posZ(), crec);
                      } else {
                        registry.fill(HIST("Tracks/THnGenSecAmb"), particle.pt(), particle.eta(), particle.mcCollision().posZ());
                      }
                      if (particle.getProcess() == TMCProcess::kPDecay) { // Particles from decay
                        if constexpr (has_reco_cent<C>) {
                          registry.fill(HIST("Tracks/Centrality/THnGenSecWeakAmb"), particle.pt(), particle.eta(), particle.mcCollision().posZ(), crec);
                        } else {
                          registry.fill(HIST("Tracks/THnGenSecWeakAmb"), particle.pt(), particle.eta(), particle.mcCollision().posZ());
                        }
                      } else { // Particles from the material
                        if constexpr (has_reco_cent<C>) {
                          registry.fill(HIST("Tracks/Centrality/THnGenSecMatAmb"), particle.pt(), particle.eta(), particle.mcCollision().posZ(), crec);
                        } else {
                          registry.fill(HIST("Tracks/THnGenSecMatAmb"), particle.pt(), particle.eta(), particle.mcCollision().posZ());
                        }
                      }
                    } else { // Primaries
                      if constexpr (has_reco_cent<C>) {
                        registry.fill(HIST("Tracks/Centrality/THnGenPrimAmb"), particle.pt(), particle.eta(), particle.mcCollision().posZ(), crec);
                      } else {
                        registry.fill(HIST("Tracks/THnGenPrimAmb"), particle.pt(), particle.eta(), particle.mcCollision().posZ());
                      }
                    }
                    break;
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  void processSecondariesMCInlcusive(soa::Join<Colls, aod::McCollisionLabels> const& collisions,
                                     MftTracksWCollsMC const& tracks,
                                     aod::McCollisions const& mccollisions,
                                     aod::McParticles const& particles,
                                     ExtBCs const& bcs)
  {
    processSecondariesMC<Colls, aod::McCollisions>(collisions, tracks, mccollisions, particles, bcs);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processSecondariesMCInlcusive, "Process secondaries checks (Inclusive)", false);

  void processSecondariesMCCentFT0C(soa::Join<CollsCentFT0C, aod::McCollisionLabels> const& collisions,
                                    MftTracksWCollsMC const& tracks,
                                    aod::McCollisions const& mccollisions,
                                    aod::McParticles const& particles,
                                    ExtBCs const& bcs)
  {
    processSecondariesMC<CollsCentFT0C, aod::McCollisions>(collisions, tracks, mccollisions, particles, bcs);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processSecondariesMCCentFT0C, "Process secondaries checks (in FT0C centrality bins)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<DndetaMFTPbPb>(cfgc)};
}
