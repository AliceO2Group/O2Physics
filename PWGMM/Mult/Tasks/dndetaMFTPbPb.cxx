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
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/MathConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/RuntimeError.h"
#include "Framework/runDataProcessing.h"
#include "MathUtils/Utils.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"

#include "TMCProcess.h"
#include "TPDGCode.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>
#include <string>
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

auto static constexpr kMinCharge = 3.f;
auto static constexpr kIntZero = 0;
auto static constexpr kZero = 0.f;

enum TrkSel {
  trkSelAll,
  trkSelNCls,
  trkSelChi2Ncl,
  trkSelEta,
  trkSelPhiCut,
  trkSelPt,
  trkSelCA,
  nTrkSel
};

enum TrkBestSel {
  trkBestSelAll,
  trkBestSelCollID,
  trkBestSelDCAxyCut,
  trkBestSelDCAzCut,
  trkBestSelNumReassoc,
  nTrkBestSel
};

enum AmbTrkType {
  kNonAmb = 0,
  kOrphan = 1,
  kNonAmbSame = 2,
  kAmb = 3,
  kAmbGt1 = 4,
  nAmbTrkType
};

struct DndetaMFTPbPb {
  SliceCache cache;

  std::array<std::shared_ptr<THnSparse>, 4> hCollAssoc;
  std::array<std::shared_ptr<THnSparse>, 4> hReAssoc;
  std::array<std::shared_ptr<THnSparse>, 6> hDCAMc;

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

  Configurable<bool> cfgDoIR{"cfgDoIR", false, "Flag to retrieve Interaction rate from CCDB"};
  Configurable<bool> cfgUseIRCut{"cfgUseIRCut", false, "Flag to cut on IR rate"};
  Configurable<bool> cfgIRCrashOnNull{"cfgIRCrashOnNull", false, "Flag to avoid CTP RateFetcher crash"};
  Configurable<std::string> cfgIRSource{"cfgIRSource", "ZNC hadronic", "Estimator of the interaction rate (Pb-Pb: ZNC hadronic)"};
  Configurable<bool> cfgUseTrackSel{"cfgUseTrackSel", false, "Flag to apply track selection"};
  Configurable<bool> cfgUseParticleSel{"cfgUseParticleSel", false, "Flag to apply particle selection"};
  Configurable<bool> cfgRemoveReassigned{"cfgRemoveReassigned", false, "Remove reassgined tracks"};

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
  } binOpt;

  struct : ConfigurableGroup {
    Configurable<bool> requireRCTFlagChecker{"requireRCTFlagChecker", false, "Check event quality in run condition table"};
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
    Configurable<bool> requireIsGoodZvtxFT0VsPV{"requireIsGoodZvtxFT0VsPV", true, "require events with PV position along z consistent (within 1 cm) between PV reconstructed using tracks and PV using FT0 A-C time difference"};
    Configurable<bool> requireRejectSameBunchPileup{"requireRejectSameBunchPileup", true, "reject collisions in case of pileup with another collision in the same foundBC"};
    Configurable<bool> requireNoCollInTimeRangeStrict{"requireNoCollInTimeRangeStrict", true, " requireNoCollInTimeRangeStrict"};
    Configurable<bool> requireNoCollInRofStrict{"requireNoCollInRofStrict", false, "requireNoCollInRofStrict"};
    Configurable<bool> requireNoCollInRofStandard{"requireNoCollInRofStandard", false, "requireNoCollInRofStandard"};
    Configurable<bool> requireNoHighMultCollInPrevRof{"requireNoHighMultCollInPrevRof", false, "requireNoHighMultCollInPrevRof"};
    Configurable<bool> requireNoCollInTimeRangeStd{"requireNoCollInTimeRangeStd", true, "reject collisions corrupted by the cannibalism, with other collisions within +/- 10 microseconds"};
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

  int mRunNumber{-1};
  uint64_t mSOR{0};
  float mMinSeconds{-1.};
  std::unordered_map<int, TH2*> gHadronicRate;
  ctpRateFetcher rateFetcher;
  TH2* gCurrentHadronicRate;
  RCTFlagsChecker rctChecker;

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

    auto hev = registry.add<TH1>("Events/hEvtSel", "hEvtSel", HistType::kTH1F,
                                 {{15, -0.5f, +14.5f}});
    hev->GetXaxis()->SetBinLabel(1, "All collisions");
    hev->GetXaxis()->SetBinLabel(2, "Ev. sel.");
    hev->GetXaxis()->SetBinLabel(3, "kIsGoodZvtxFT0vsPV");
    hev->GetXaxis()->SetBinLabel(4, "NoSameBunchPileup");
    hev->GetXaxis()->SetBinLabel(5, "Z-vtx cut");
    hev->GetXaxis()->SetBinLabel(6, "kNoCollInTimeRangeStd");
    hev->GetXaxis()->SetBinLabel(7, "kNoCollInTimeRangeNarrow");
    hev->GetXaxis()->SetBinLabel(8, "kNoCollInTimeRangeStrict");
    hev->GetXaxis()->SetBinLabel(9, "kNoCollInRofStrict");
    hev->GetXaxis()->SetBinLabel(10, "kNoCollInRofStandard");
    hev->GetXaxis()->SetBinLabel(11, "kNoHighMultCollInPrevRof");
    hev->GetXaxis()->SetBinLabel(12, "Below min occup.");
    hev->GetXaxis()->SetBinLabel(13, "Above max occup.");
    hev->GetXaxis()->SetBinLabel(14, "RCT Flag Checker");

    registry.add("Tracks/hBestTrkSel", "Number of best tracks; Cut; #Tracks Passed Cut", {HistType::kTH1F, {{nTrkBestSel, -0.5, +nTrkBestSel - 0.5}}});
    registry.get<TH1>(HIST("Tracks/hBestTrkSel"))->GetXaxis()->SetBinLabel(trkBestSelAll + 1, "All");
    registry.get<TH1>(HIST("Tracks/hBestTrkSel"))->GetXaxis()->SetBinLabel(trkBestSelCollID + 1, "Assigned (ID>=0)");
    registry.get<TH1>(HIST("Tracks/hBestTrkSel"))->GetXaxis()->SetBinLabel(trkBestSelDCAxyCut + 1, "DCA xy cut");
    registry.get<TH1>(HIST("Tracks/hBestTrkSel"))->GetXaxis()->SetBinLabel(trkBestSelDCAzCut + 1, "DCA z cut");
    registry.get<TH1>(HIST("Tracks/hBestTrkSel"))->GetXaxis()->SetBinLabel(trkBestSelNumReassoc + 1, "Reassociated");

    registry.add("Tracks/hTrkSel", "Number of tracks; Cut; #Tracks Passed Cut", {HistType::kTH1F, {{nTrkSel, -0.5, +nTrkSel - 0.5}}});
    registry.get<TH1>(HIST("Tracks/hTrkSel"))->GetXaxis()->SetBinLabel(trkSelAll + 1, "All");
    registry.get<TH1>(HIST("Tracks/hTrkSel"))->GetXaxis()->SetBinLabel(trkSelNCls + 1, "Ncl cut");
    registry.get<TH1>(HIST("Tracks/hTrkSel"))->GetXaxis()->SetBinLabel(trkSelChi2Ncl + 1, "#chi^{2}/Ncl cut");
    registry.get<TH1>(HIST("Tracks/hTrkSel"))->GetXaxis()->SetBinLabel(trkSelEta + 1, "#eta cut");
    registry.get<TH1>(HIST("Tracks/hTrkSel"))->GetXaxis()->SetBinLabel(trkSelPhiCut + 1, "#varphi cut");
    registry.get<TH1>(HIST("Tracks/hTrkSel"))->GetXaxis()->SetBinLabel(trkSelPt + 1, "#it{p}_{T} cut");
    registry.get<TH1>(HIST("Tracks/hTrkSel"))->GetXaxis()->SetBinLabel(trkSelCA + 1, "Tracking algorithm (CA)");

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
      // qaregistry.add({"Tracks/hAmbTrackType",
      //                 " ; Ambiguous track type",
      //                 {HistType::kTH1F, {{5, -0.5, 4.5}}}});

      qaregistry.add("Tracks/hAmbTrackType", "hAmbTrackType", {HistType::kTH1F, {{AmbTrkType::nAmbTrkType, -0.5, +AmbTrkType::nAmbTrkType - 0.5}}});
      std::string labelAmbiguity[AmbTrkType::nAmbTrkType];
      labelAmbiguity[AmbTrkType::kOrphan] = "orphan";
      labelAmbiguity[AmbTrkType::kNonAmb] = "nonAmbiguous";
      labelAmbiguity[AmbTrkType::kNonAmbSame] = "trkInCollTabHasSameAssoc";
      labelAmbiguity[AmbTrkType::kAmb] = "trkInCollTabHasDiffAssoc";
      labelAmbiguity[AmbTrkType::kAmbGt1] = "trkInCollTabHasGt1Assoc";
      qaregistry.get<TH1>(HIST("Tracks/hAmbTrackType"))->SetMinimum(0.1);

      for (int iBin = 0; iBin < AmbTrkType::nAmbTrkType; iBin++) {
        qaregistry.get<TH1>(HIST("Tracks/hAmbTrackType"))->GetXaxis()->SetBinLabel(iBin + 1, labelAmbiguity[iBin].data());
      }
    }

    if (doprocessCollAssocMC) {
      // tracks not associated to any collision
      hCollAssoc[0] = qaregistry.add<THnSparse>("TrackToColl/hNonAssocTracks", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis});
      // tracks associasted to a collision
      hCollAssoc[1] = qaregistry.add<THnSparse>("TrackToColl/hAssocTracks", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis});
      // tracks associated to the correct collision considering only first reco collision (based on the MC collision index)
      hCollAssoc[2] = qaregistry.add<THnSparse>("TrackToColl/hGoodAssocTracks", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis});
      // tracks associated to the correct collision considering all ambiguous reco collisions (based on the MC collision index)
      hCollAssoc[3] = qaregistry.add<THnSparse>("TrackToColl/hGoodAssocTracksAmb", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis});
      qaregistry.add("TrackToColl/histFracTracksFakeMcColl", "Fraction of tracks originating from fake collision; fraction; entries", {HistType::kTH1F, {{101, 0., 1.01}}});
      qaregistry.add("TrackToColl/histFracGoodTracks", "Fraction of tracks originating from the correct collision; fraction; entries", {HistType::kTH1F, {{101, 0., 1.01}}});
      qaregistry.add("TrackToColl/histAmbTrackNumColls", "Number of collisions associated to an ambiguous track; no. collisions; entries", {HistType::kTH1F, {{30, -0.5, 29.5}}});
    }

    if (doprocessReAssocMC) {
      // tracks not associated to any collision
      hReAssoc[0] = qaregistry.add<THnSparse>("ReAssoc/hNonAssocTracks", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis});
      // tracks associasted to a collision
      hReAssoc[1] = qaregistry.add<THnSparse>("ReAssoc/hAssocTracks", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis});
      // tracks associated to the correct collision considering only first reco collision (based on the MC collision index)
      hReAssoc[2] = qaregistry.add<THnSparse>("ReAssoc/hGoodAssocTracks", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis});
      // tracks associated to the correct collision considering all ambiguous reco collisions (based on the MC collision index)
      hReAssoc[3] = qaregistry.add<THnSparse>("ReAssoc/hGoodAssocTracksAmb", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {ptAxis, etaAxis, deltaZAxis});
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
        registry.add({"Events/EvtGenRec", ";status", {HistType::kTH1F, {{3, 0.5, 3.5}}}});
        auto heff = registry.get<TH1>(HIST("Events/EvtGenRec"));
        auto* h = heff->GetXaxis();
        h->SetBinLabel(1, "All generated");
        h->SetBinLabel(2, "All reconstructed");
        h->SetBinLabel(3, "Selected reconstructed");
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
        registry.add({"Events/Centrality/EvtGenRec", ";status;centrality", {HistType::kTH2F, {{3, 0.5, 3.5}, centralityAxis}}});
        auto heff = registry.get<TH2>(HIST("Events/Centrality/EvtGenRec"));
        auto* h = heff->GetXaxis();
        h->SetBinLabel(1, "All generated");
        h->SetBinLabel(2, "All reconstructed");
        h->SetBinLabel(3, "Selected reconstructed");
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

    if (doprocessDCAReassocMcInclusive || doprocessDCAReassocMcCentFT0C) {
      auto hNevt = registry.add<TH1>("Events/hNGenRecCollsReassoc", "Number of generated and reconstructed MC collisions", HistType::kTH1F, {{3, 0.5, 3.5}});
      hNevt->GetXaxis()->SetBinLabel(1, "Reconstructed collisions");
      hNevt->GetXaxis()->SetBinLabel(2, "Generated collisions");
      if (doprocessDCAReassocMcInclusive) {
        registry.add({"Events/EvtGenRecReassoc", ";status", {HistType::kTH2F, {{3, 0.5, 3.5}, occupancyAxis}}});
        auto heff = registry.get<TH2>(HIST("Events/EvtGenRecReassoc"));
        auto* h = heff->GetXaxis();
        h->SetBinLabel(1, "All generated");
        h->SetBinLabel(2, "All reconstructed");
        h->SetBinLabel(3, "Selected reconstructed");
        registry.add({"Tracks/THnDCAxyBestRec", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm); DCA_{XY} (cm);  DCA_{Z} (cm)", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis, dcaxyAxis, dcazAxis}}});
        registry.add({"Tracks/THnDCAxyBestRecFake", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm); DCA_{XY} (cm);  DCA_{Z} (cm)", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis, dcaxyAxis, dcazAxis}}});
        registry.add({"Tracks/THnDCAxyBestGenPrim", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm); DCA_{XY} (cm);  DCA_{Z} (cm)", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis, dcaxyAxis, dcazAxis}}});
        registry.add({"Tracks/THnDCAxyBestGenTruthPrim", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm); DCA_{XY} (cm);  DCA_{Z} (cm)", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis, dcaxyAxis, dcazAxis}}});
        registry.add({"Tracks/THnDCAxyBestGenPrimWrongColl", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm); DCA_{XY} (cm);  DCA_{Z} (cm)", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis, dcaxyAxis, dcazAxis}}});
        registry.add({"Tracks/THnDCAxyBestGenTruthPrimWrongColl", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm); DCA_{XY} (cm);  DCA_{Z} (cm)", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis, dcaxyAxis, dcazAxis}}});
        registry.add({"Tracks/THnDCAxyBestGenSec", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm); DCA_{XY} (cm);  DCA_{Z} (cm)", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis, dcaxyAxis, dcazAxis}}});
        registry.add({"Tracks/THnDCAxyBestGenTruthSec", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm); DCA_{XY} (cm);  DCA_{Z} (cm)", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis, dcaxyAxis, dcazAxis}}});
        registry.add({"Tracks/THnDCAxyBestGenSecWrongColl", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm); DCA_{XY} (cm);  DCA_{Z} (cm)", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis, dcaxyAxis, dcazAxis}}});
        registry.add({"Tracks/THnDCAxyBestGenTruthSecWrongColl", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm); DCA_{XY} (cm);  DCA_{Z} (cm)", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis, dcaxyAxis, dcazAxis}}});
        registry.add({"Tracks/THnDCAxyBestGenSecWeak", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm); DCA_{XY} (cm);  DCA_{Z} (cm)", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis, dcaxyAxis, dcazAxis}}});
        registry.add({"Tracks/THnDCAxyBestGenSecMat", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm); DCA_{XY} (cm);  DCA_{Z} (cm)", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis, dcaxyAxis, dcazAxis}}});
      }
      if (doprocessDCAReassocMcCentFT0C) {
        registry.add({"Events/Centrality/EvtGenRecReassoc", ";status;centrality", {HistType::kTHnSparseF, {{3, 0.5, 3.5}, centralityAxis, occupancyAxis}}});
        auto heff = registry.get<THnSparse>(HIST("Events/Centrality/EvtGenRecReassoc"));
        heff->GetAxis(0)->SetBinLabel(1, "All generated");
        heff->GetAxis(0)->SetBinLabel(2, "All reconstructed");
        heff->GetAxis(0)->SetBinLabel(3, "Selected reconstructed");
        registry.add({"Tracks/Centrality/THnDCAxyBestRec", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm); DCA_{XY} (cm);  DCA_{Z} (cm); centrality", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis, dcaxyAxis, dcazAxis, centralityAxis}}});
        registry.add({"Tracks/Centrality/THnDCAxyBestRecFake", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm); DCA_{XY} (cm);  DCA_{Z} (cm); centrality", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis, dcaxyAxis, dcazAxis, centralityAxis}}});
        registry.add({"Tracks/Centrality/THnDCAxyBestGenPrim", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm); DCA_{XY} (cm);  DCA_{Z} (cm); centrality", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis, dcaxyAxis, dcazAxis, centralityAxis}}});
        registry.add({"Tracks/Centrality/THnDCAxyBestGenTruthPrim", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm); DCA_{XY} (cm);  DCA_{Z} (cm); centrality", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis, dcaxyAxis, dcazAxis, centralityAxis}}});
        registry.add({"Tracks/Centrality/THnDCAxyBestGenPrimWrongColl", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm); DCA_{XY} (cm);  DCA_{Z} (cm); centrality", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis, dcaxyAxis, dcazAxis, centralityAxis}}});
        registry.add({"Tracks/Centrality/THnDCAxyBestGenTruthPrimWrongColl", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm); DCA_{XY} (cm);  DCA_{Z} (cm); centrality", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis, dcaxyAxis, dcazAxis, centralityAxis}}});
        registry.add({"Tracks/Centrality/THnDCAxyBestGenSec", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm); DCA_{XY} (cm);  DCA_{Z} (cm); centrality", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis, dcaxyAxis, dcazAxis, centralityAxis}}});
        registry.add({"Tracks/Centrality/THnDCAxyBestGenTruthSec", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm); DCA_{XY} (cm);  DCA_{Z} (cm); centrality", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis, dcaxyAxis, dcazAxis, centralityAxis}}});
        registry.add({"Tracks/Centrality/THnDCAxyBestGenSecWrongColl", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm); DCA_{XY} (cm);  DCA_{Z} (cm); centrality", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis, dcaxyAxis, dcazAxis, centralityAxis}}});
        registry.add({"Tracks/Centrality/THnDCAxyBestGenTruthSecWrongColl", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm); DCA_{XY} (cm);  DCA_{Z} (cm); centrality", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis, dcaxyAxis, dcazAxis, centralityAxis}}});
        registry.add({"Tracks/Centrality/THnDCAxyBestGenSecWeak", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm); DCA_{XY} (cm);  DCA_{Z} (cm); centrality", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis, dcaxyAxis, dcazAxis, centralityAxis}}});
        registry.add({"Tracks/Centrality/THnDCAxyBestGenSecMat", ";  p_{T} (GeV/c); #eta; Z_{vtx} (cm); DCA_{XY} (cm);  DCA_{Z} (cm); centrality", {HistType::kTHnSparseF, {ptAxis, etaAxis, zAxis, dcaxyAxis, dcazAxis, centralityAxis}}});
      }
    }
  }

  /// Filters - tracks
  Filter filtTrkEta = (aod::fwdtrack::eta < trackCuts.maxEta) &&
                      (aod::fwdtrack::eta > trackCuts.minEta);
  Filter filtATrackID = (aod::fwdtrack::bestCollisionId >= kIntZero);
  Filter filtATrackDCAxy = (nabs(aod::fwdtrack::bestDCAXY) < trackCuts.maxDCAxy);
  Filter filtATrackDCAz = (nabs(aod::fwdtrack::bestDCAZ) < trackCuts.maxDCAz);

  /// Filters - mc particles
  Filter primaries = (aod::mcparticle::flags & (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary) == (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary && (aod::mcparticle::eta < trackCuts.maxEta) && (aod::mcparticle::eta > trackCuts.minEta);

  /// Joined tables
  using FullBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;
  using CollBCs = soa::Join<aod::BCsWithTimestamps, aod::Run3MatchedToBCSparse>;
  // Collisions
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
  // Tracks
  using MFTTracksLabeled = soa::Join<aod::MFTTracks, aod::McMFTTrackLabels>;
  using MftTracksWColls = soa::Join<aod::MFTTracks, aod::MFTTrkCompColls>;
  using MftTracksWCollsMC = soa::Join<aod::MFTTracks, aod::MFTTrkCompColls, aod::McMFTTrackLabels>;

  using BestTracksMC = soa::Join<aod::MFTTracks, aod::BestCollisionsFwd3d, aod::McMFTTrackLabels>;

  /// Filtered tables
  using FiltMftTracks = soa::Filtered<aod::MFTTracks>;
  using FiltMcMftTracks = soa::Filtered<MFTTracksLabeled>;
  using FiltBestTracks = soa::Filtered<aod::BestCollisionsFwd3d>;
  using FiltMcBestTracks = soa::Filtered<BestTracksMC>;

  using FiltParticles = soa::Filtered<aod::McParticles>;

  template <bool fillHis = true, typename B>
  bool isBestTrackSelected(const B& besttrack)
  {
    if constexpr (fillHis) {
      registry.fill(HIST("Tracks/hBestTrkSel"), trkBestSelAll);
    }
    if (besttrack.bestCollisionId() < kIntZero) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Tracks/hBestTrkSel"), trkBestSelCollID);
    }
    if (std::abs(besttrack.bestDCAXY()) >= trackCuts.maxDCAxy) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Tracks/hBestTrkSel"), trkBestSelDCAxyCut);
    }
    if (std::abs(besttrack.bestDCAZ()) >= trackCuts.maxDCAxy) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Tracks/hBestTrkSel"), trkBestSelDCAzCut);
    }
    return true;
  }

  template <bool fillHis = true, typename T>
  bool isTrackSelected(const T& track)
  {
    if constexpr (fillHis) {
      registry.fill(HIST("Tracks/hTrkSel"), trkSelAll);
    }
    if (track.nClusters() < trackCuts.minNclusterMft) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Tracks/hTrkSel"), trkSelNCls);
    }
    if (trackCuts.useChi2Cut) {
      float nclMft = std::max(2.0f * track.nClusters() - 5.0f, 1.0f);
      float mftChi2NCl = track.chi2() / nclMft;
      if (mftChi2NCl > trackCuts.maxChi2NCl)
        return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Tracks/hTrkSel"), trkSelChi2Ncl);
    }
    if (track.eta() < trackCuts.minEta || track.eta() > trackCuts.maxEta) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Tracks/hTrkSel"), trkSelEta);
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
      registry.fill(HIST("Tracks/hTrkSel"), trkSelPhiCut);
    }
    if (trackCuts.usePtCut && track.pt() < trackCuts.minPt) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Tracks/hTrkSel"), trkSelPt);
    }
    if (trackCuts.requireCA && !track.isCA()) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Tracks/hTrkSel"), trkSelCA);
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
        if (phi < kZero || TwoPI < phi) {
          continue;
        }
        if constexpr (has_reco_cent<C>) {
          registry.fill(HIST("Tracks/Centrality/EtaZvtx"), track.eta(), z, c, occ);
          registry.fill(HIST("Tracks/Centrality/PhiEta"), phi, track.eta(), c, occ);
        } else {
          registry.fill(HIST("Tracks/EtaZvtx"), track.eta(), z, occ);
          registry.fill(HIST("Tracks/PhiEta"), phi, track.eta(), occ);
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
        if (phi < kZero || TwoPI < phi) {
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
          registry.fill(HIST("Tracks/hBestTrkSel"), trkBestSelNumReassoc);
          float phi = itrack.phi();
          o2::math_utils::bringTo02Pi(phi);
          if (phi < kZero || TwoPI < phi) {
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
      if (phi < kZero || TwoPI < phi) {
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
      if (cfgUseParticleSel && !isParticleSelected(particle)) {
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
      registry.fill(HIST("Events/hEvtSel"), 0);
    }
    if (!collision.sel8()) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), 1);
    }
    if (eventCuts.requireIsGoodZvtxFT0VsPV && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), 2);
    }
    if (eventCuts.requireRejectSameBunchPileup && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), 3);
    }
    if (collision.posZ() <= eventCuts.minZvtx || collision.posZ() >= eventCuts.maxZvtx) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), 4);
    }
    if (eventCuts.requireNoCollInTimeRangeStd &&
        !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), 5);
    }
    if (eventCuts.requireNoCollInTimeRangeNarrow &&
        !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow)) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), 6);
    }
    if (eventCuts.requireNoCollInTimeRangeStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict)) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), 7);
    }
    if (eventCuts.requireNoCollInRofStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStrict)) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), 8);
    }
    if (eventCuts.requireNoCollInRofStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), 9);
    }
    if (eventCuts.requireNoHighMultCollInPrevRof && !collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), 10);
    }
    if (eventCuts.minOccupancy >= 0 &&
        getOccupancy(collision, eventCuts.occupancyEstimator) <
          eventCuts.minOccupancy) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), 11);
    }
    if (eventCuts.maxOccupancy >= 0 &&
        getOccupancy(collision, eventCuts.occupancyEstimator) >
          eventCuts.maxOccupancy) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), 12);
    }

    if (rctCuts.requireRCTFlagChecker && !rctChecker(collision)) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), 13);
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
    return std::abs(charge) >= kMinCharge;
  }

  template <bool isCent, typename P>
  void fillHistMC(P const& particles, float c, float occ, float zvtx,
                  bool const gtZeroColl)
  {
    for (auto const& particle : particles) {
      if (!isChrgParticle(particle.pdgCode())) {
        continue;
      }
      if (cfgUseParticleSel && !isParticleSelected(particle)) {
        continue;
      }

      float phi = particle.phi();
      o2::math_utils::bringTo02Pi(phi);
      if (phi < kZero || TwoPI < phi) {
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
        if (phi < kZero || TwoPI < phi) {
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

    if (!isGoodEvent<true>(collision)) {
      return;
    }

    if (cfgDoIR) {
      initHadronicRate(bc);
      float ir = !cfgIRSource.value.empty() ? rateFetcher.fetch(ccdb.service, bc.timestamp(), bc.runNumber(), cfgIRSource, cfgIRCrashOnNull) * 1.e-3 : -1;
      if constexpr (has_reco_cent<C>) {
        registry.fill(HIST("Events/Centrality/hInteractionRate"), c, occ, ir);
      } else {
        registry.fill(HIST("Events/hInteractionRate"), occ, ir);
      }
      float seconds = bc.timestamp() * 1.e-3 - mMinSeconds;
      if (cfgUseIRCut && (ir < eventCuts.minIR || ir > eventCuts.maxIR)) { // cut on hadronic rate
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
    soa::SmallGroups<aod::BestCollisionsFwd3d> const& besttracks, CollBCs const& /*bcs*/)
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

    if (cfgDoIR) {
      initHadronicRate(bc);
      float ir = !cfgIRSource.value.empty() ? rateFetcher.fetch(ccdb.service, bc.timestamp(), bc.runNumber(), cfgIRSource, cfgIRCrashOnNull) * 1.e-3 : -1;
      if constexpr (has_reco_cent<C>) {
        registry.fill(HIST("Events/Centrality/hInteractionRate"), c, occ, ir);
      } else {
        registry.fill(HIST("Events/hInteractionRate"), occ, ir);
      }
      float seconds = bc.timestamp() * 1.e-3 - mMinSeconds;
      if (cfgUseIRCut && (ir < eventCuts.minIR || ir > eventCuts.maxIR)) { // cut on hadronic rate
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
    soa::SmallGroups<aod::BestCollisionsFwd3d> const& besttracks, CollBCs const& bcs)
  {
    processDatawBestTracks<Colls>(collision, tracks, besttracks, bcs);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processDatawBestTracksInclusive,
                 "Count tracks based on BestCollisionsFwd3d table (inclusive)",
                 false);

  void processDatawBestTracksCentFT0C(
    CollsCentFT0C::iterator const& collision, FiltMftTracks const& tracks,
    soa::SmallGroups<aod::BestCollisionsFwd3d> const& besttracks, CollBCs const& bcs)
  {
    processDatawBestTracks<CollsCentFT0C>(collision, tracks, besttracks, bcs);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processDatawBestTracksCentFT0C,
                 "Count tracks in FT0C centrality bins based on BestCollisionsFwd3d table",
                 false);

  void processDatawBestTracksCentFT0CVariant1(
    CollsCentFT0CVariant1::iterator const& collision,
    FiltMftTracks const& tracks,
    soa::SmallGroups<aod::BestCollisionsFwd3d> const& besttracks, CollBCs const& bcs)
  {
    processDatawBestTracks<CollsCentFT0CVariant1>(collision, tracks, besttracks, bcs);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processDatawBestTracksCentFT0CVariant1,
                 "Count tracks in FT0CVariant1 centrality bins based on "
                 "BestCollisionsFwd3d table",
                 false);

  void processDatawBestTracksCentFT0M(
    CollsCentFT0M::iterator const& collision, FiltMftTracks const& tracks,
    soa::SmallGroups<aod::BestCollisionsFwd3d> const& besttracks, CollBCs const& bcs)
  {
    processDatawBestTracks<CollsCentFT0M>(collision, tracks, besttracks, bcs);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processDatawBestTracksCentFT0M,
                 "Count tracks in FT0M centrality bins based on BestCollisionsFwd3d table",
                 false);

  void processDatawBestTracksCentNGlobal(
    CollsCentNGlobal::iterator const& collision, FiltMftTracks const& tracks,
    soa::SmallGroups<aod::BestCollisionsFwd3d> const& besttracks, CollBCs const& bcs)
  {
    processDatawBestTracks<CollsCentNGlobal>(collision, tracks, besttracks, bcs);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processDatawBestTracksCentNGlobal,
                 "Count tracks in NGlobal centrality bins based on "
                 "BestCollisionsFwd3d table",
                 false);

  void processDatawBestTracksCentMFT(
    CollsCentMFT::iterator const& collision, FiltMftTracks const& tracks,
    soa::SmallGroups<aod::BestCollisionsFwd3d> const& besttracks, CollBCs const& bcs)
  {
    processDatawBestTracks<CollsCentMFT>(collision, tracks, besttracks, bcs);
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
        if (cfgUseParticleSel && !isParticleSelected(particle)) {
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
        if (cfgUseParticleSel && !isParticleSelected(particle)) {
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
      if (cfgUseParticleSel && !isParticleSelected(particle)) {
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

  /// @brief process function to check ambiguous tracks
  void processCheckAmbiguousMftTracks(aod::Collisions const&, MftTracksWColls const& tracks)
  {
    for (auto const& track : tracks) {
      auto trkCollId = track.has_collision() ? track.collisionId() : -1;
      auto ids = track.compatibleCollIds();
      if (ids.empty() || (ids.size() == 1 && trkCollId == ids[0])) {
        qaregistry.fill(HIST("Tracks/hMftTracksAmbDegreeWithTrivial"), track.compatibleCollIds().size());
        if (ids.empty()) {
          qaregistry.fill(HIST("Tracks/hAmbTrackType"), AmbTrkType::kOrphan);
        }
        if (ids.size() == 1 && trkCollId == ids[0]) {
          qaregistry.fill(HIST("Tracks/hAmbTrackType"), AmbTrkType::kNonAmb);
        }
        continue;
      }
      qaregistry.fill(HIST("Tracks/hMftTracksAmbDegree"), track.compatibleCollIds().size());

      if (track.compatibleCollIds().size() > 0) {
        if (track.compatibleCollIds().size() == 1) {
          if (track.collisionId() != track.compatibleCollIds()[0]) {
            qaregistry.fill(HIST("Tracks/hAmbTrackType"), AmbTrkType::kAmb);
          } else {
            qaregistry.fill(HIST("Tracks/hAmbTrackType"), AmbTrkType::kNonAmbSame);
          }
        } else {
          qaregistry.fill(HIST("Tracks/hAmbTrackType"), AmbTrkType::kAmbGt1);
        }
      }
    }
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processCheckAmbiguousMftTracks, "Process checks for Ambiguous MFT tracks (inclusive)", false);

  Partition<MftTracksWCollsMC> tracksInAcc = (aod::fwdtrack::eta < trackCuts.maxEta) && (aod::fwdtrack::eta > trackCuts.minEta);

  template <typename C>
  void processCheckAssocMC(C const& collisions,
                           MftTracksWCollsMC const& tracks,
                           aod::McParticles const& /*particles*/,
                           aod::McCollisions const& /*mccollisions*/
  )
  {
    for (const auto& collision : collisions) {
      if (!collision.has_mcCollision()) {
        continue;
      }
      // auto mcCollision = collision.template mcCollision_as<aod::McCollisions>();
      auto tracksInColl = tracksInAcc->sliceByCached(aod::fwdtrack::collisionId, collision.globalIndex(), cache);
      int nTrk = 0, nFakeTrk = 0, nGoodTrk = 0;
      for (const auto& track : tracksInColl) {
        if (!track.has_mcParticle()) {
          continue;
        }
        nTrk++;
        auto particle = track.mcParticle();

        if ((particle.mcCollisionId() != collision.mcCollision().globalIndex())) {
          nFakeTrk++;
          continue;
        }
        if (collision.mcCollisionId() == particle.mcCollisionId()) {
          nGoodTrk++;
        }
      }
      float frac = (nTrk > 0) ? static_cast<float>(nGoodTrk) / nTrk : -1.;
      qaregistry.fill(HIST("TrackToColl/histFracGoodTracks"), frac);
      float fracFake = (nTrk > 0) ? static_cast<float>(nFakeTrk) / nTrk : -1.;
      qaregistry.fill(HIST("TrackToColl/histFracTracksFakeMcColl"), fracFake);
    }

    for (auto const& track : tracks) {
      uint index = uint(track.collisionId() >= 0);
      if (track.has_mcParticle()) {
        // auto particle = track.mcParticle_as<FiltParticles>();
        auto particle = track.mcParticle();
        bool isAmbiguous = (track.compatibleCollIds().size() != 1);
        if (isAmbiguous) {
          qaregistry.fill(HIST("TrackToColl/histAmbTrackNumColls"), track.compatibleCollIds().size());
        }
        float deltaZ = -999.f;
        if (index) {
          auto collision = track.collision_as<CollisionsWithMCLabels>();
          auto mcCollision = particle.mcCollision_as<aod::McCollisions>();
          deltaZ = collision.posZ() - mcCollision.posZ();
          if (collision.has_mcCollision() && collision.mcCollisionId() == particle.mcCollisionId()) {
            hCollAssoc[index + 1]->Fill(track.pt(), track.eta(), deltaZ);
          } else {
            if (isAmbiguous) {
              for (const auto& collIdx : track.compatibleCollIds()) {
                auto ambCollision = collisions.rawIteratorAt(collIdx);
                if (ambCollision.has_mcCollision() && ambCollision.mcCollisionId() == particle.mcCollisionId()) {
                  hCollAssoc[index + 2]->Fill(track.pt(), track.eta(), deltaZ);
                  break;
                }
              }
            }
          }
        }
        hCollAssoc[index]->Fill(track.pt(), track.eta(), deltaZ);
      } else {
        hCollAssoc[index]->Fill(track.pt(), track.eta(), -999.f);
      }
    }
  }

  void processCollAssocMC(CollisionsWithMCLabels const& collisions,
                          MftTracksWCollsMC const& tracks,
                          aod::McParticles const& particles,
                          aod::McCollisions const& mccollisions)
  {
    processCheckAssocMC<CollisionsWithMCLabels>(collisions, tracks, particles, mccollisions);
  }
  PROCESS_SWITCH(DndetaMFTPbPb, processCollAssocMC, "Process collision-association information, requires extra table from TrackToCollisionAssociation task (fillTableOfCollIdsPerTrack=true)", false);

  template <typename C>
  void processCheckReAssocMC(C const& /*collisions*/,
                             soa::SmallGroups<aod::BestCollisionsFwd3d> const& besttracks,
                             FiltMcMftTracks const& /*tracks*/,
                             FiltParticles const& /*particles*/,
                             aod::McCollisions const& /*mccollisions*/
  )
  {
    for (auto const& track : besttracks) {
      uint index = uint(track.bestCollisionId() >= 0); // assigned
      if (!isBestTrackSelected<false>(track)) {
        continue;
      }
      auto itrack = track.mfttrack_as<FiltMcMftTracks>();

      if (cfgRemoveReassigned) {
        if (itrack.collisionId() != track.bestCollisionId()) {
          continue;
        }
      }
      if (!isTrackSelected<false>(itrack)) {
        continue;
      }
      if (itrack.has_mcParticle()) {
        auto particle = itrack.mcParticle_as<FiltParticles>();

        float deltaZ = -999.f;
        if (index) {
          auto collision = itrack.collision_as<CollisionsWithMCLabels>();
          auto mcCollision = particle.mcCollision_as<aod::McCollisions>();
          deltaZ = collision.posZ() - mcCollision.posZ();
          if (collision.has_mcCollision() && collision.mcCollisionId() == particle.mcCollisionId()) {
            hReAssoc[index + 1]->Fill(itrack.pt(), itrack.eta(), deltaZ);
          } else {
            hReAssoc[index + 2]->Fill(itrack.pt(), itrack.eta(), deltaZ);
          }
        }
        hReAssoc[index]->Fill(itrack.pt(), itrack.eta(), deltaZ);
      } else {
        hReAssoc[index]->Fill(itrack.pt(), itrack.eta(), -999.f);
      }
    }
  }

  void processReAssocMC(CollisionsWithMCLabels const& collisions,
                        soa::SmallGroups<aod::BestCollisionsFwd3d> const& besttracks,
                        FiltMcMftTracks const& tracks,
                        FiltParticles const& particles,
                        aod::McCollisions const& mccollisions)
  {
    processCheckReAssocMC<CollisionsWithMCLabels>(collisions, besttracks, tracks, particles, mccollisions);
  }
  PROCESS_SWITCH(DndetaMFTPbPb, processReAssocMC, "Process re-association information based on BestCollisionsFwd3d table", false);

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

  Preslice<MftTracksWCollsMC> mftTrkCompCollperCol = o2::aod::fwdtrack::collisionId;

  /// @brief process template function for DCA MC checks
  template <typename C, typename MC>
  void processSecondariesMC(typename soa::Join<C, aod::McCollisionLabels> const& collisions,
                            MftTracksWCollsMC const& tracks,
                            MC const& mcCollisions,
                            aod::McParticles const& /*particles*/
  )
  {
    registry.fill(HIST("Events/hNGenRecColls"), 1.f, collisions.size());
    registry.fill(HIST("Events/hNGenRecColls"), 2.f, mcCollisions.size());

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

    if constexpr (has_reco_cent<C>) {
      registry.fill(HIST("Events/Centrality/EvtGenRec"), 1., cGen);
    } else {
      registry.fill(HIST("Events/EvtGenRec"), 1.);
    }

    for (const auto& collision : collisions) {
      float crec = getRecoCent(collision);

      if constexpr (has_reco_cent<C>) {
        registry.fill(HIST("Events/Centrality/EvtGenRec"), 2., crec);
      } else {
        registry.fill(HIST("Events/EvtGenRec"), 2.);
      }

      if (!isGoodEvent<false>(collision)) {
        continue;
      }

      if constexpr (has_reco_cent<C>) {
        registry.fill(HIST("Events/Centrality/EvtGenRec"), 3., crec);
      } else {
        registry.fill(HIST("Events/EvtGenRec"), 3.);
      }

      if (!collision.has_mcCollision()) {
        continue;
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
        bool isAmbiguous = (ids.size() != 1);

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

        uint index = uint(track.collisionId() >= 0);

        if (!track.has_mcParticle()) {
          LOGP(debug, "No MC particle for track, skip...");
          continue;
        }

        auto particle = track.template mcParticle_as<aod::McParticles>();
        if (!isChrgParticle(particle.pdgCode())) {
          continue;
        }
        if (particle.eta() <= trackCuts.minEta || particle.eta() >= trackCuts.maxEta) {
          continue;
        }
        if (cfgUseParticleSel && !isParticleSelected(particle)) {
          continue;
        }

        if (index) {
          auto mcCollision = particle.template mcCollision_as<aod::McCollisions>();
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
                auto ambCollision = collisions.rawIteratorAt(collIdx);
                if (ambCollision.has_mcCollision() && ambCollision.mcCollisionId() == particle.mcCollisionId()) {
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

  void processSecondariesMCInlcusive(soa::Join<Colls, aod::McCollisionLabels> const& collisions,
                                     MftTracksWCollsMC const& tracks,
                                     aod::McCollisions const& mccollisions,
                                     aod::McParticles const& particles)
  {
    processSecondariesMC<Colls, aod::McCollisions>(collisions, tracks, mccollisions, particles);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processSecondariesMCInlcusive, "Process secondaries checks (Inclusive)", false);

  void processSecondariesMCCentFT0C(soa::Join<CollsCentFT0C, aod::McCollisionLabels> const& collisions,
                                    MftTracksWCollsMC const& tracks,
                                    aod::McCollisions const& mccollisions,
                                    aod::McParticles const& particles)
  {
    processSecondariesMC<CollsCentFT0C, aod::McCollisions>(collisions, tracks, mccollisions, particles);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processSecondariesMCCentFT0C, "Process secondaries checks (in FT0C centrality bins)", false);

  template <typename C, typename MC>
  void processDCAReassocMc(typename soa::Join<C, aod::McCollisionLabels> const& collisions,
                           MC const& mcCollisions,
                           aod::McParticles const& /*particles*/,
                           BestTracksMC const& besttracks,
                           FiltMcMftTracks const& /*tracks*/
  )
  {
    registry.fill(HIST("Events/hNGenRecCollsReassoc"), 1.f, collisions.size());
    registry.fill(HIST("Events/hNGenRecCollsReassoc"), 2.f, mcCollisions.size());

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

    if constexpr (has_reco_cent<C>) {
      registry.fill(HIST("Events/Centrality/EvtGenRecReassoc"), 1., cGen, occGen);
    } else {
      registry.fill(HIST("Events/EvtGenRecReassoc"), 1., occGen);
    }

    for (const auto& collision : collisions) {
      auto occ = getOccupancy(collision, eventCuts.occupancyEstimator);
      float crec = getRecoCent(collision);

      if constexpr (has_reco_cent<C>) {
        registry.fill(HIST("Events/Centrality/EvtGenRecReassoc"), 2., crec, occ);
      } else {
        registry.fill(HIST("Events/EvtGenRecReassoc"), 2., occ);
      }

      if (!isGoodEvent<false>(collision)) {
        continue;
      }

      if constexpr (has_reco_cent<C>) {
        registry.fill(HIST("Events/Centrality/EvtGenRecReassoc"), 3., crec, occ);
      } else {
        registry.fill(HIST("Events/EvtGenRecReassoc"), 3., occ);
      }

      if (!collision.has_mcCollision()) {
        continue;
      }

      auto perCollisionASample = besttracks.sliceBy(perColU, collision.globalIndex());
      for (auto const& atrack : perCollisionASample) {
        if (!isBestTrackSelected<false>(atrack)) {
          continue;
        }
        auto itrack = atrack.template mfttrack_as<FiltMcMftTracks>();

        if (!isTrackSelected<false>(itrack)) {
          continue;
        }
        float phi = itrack.phi();
        o2::math_utils::bringTo02Pi(phi);
        if (phi < kZero || TwoPI < phi) {
          continue;
        }

        if (!itrack.has_collision()) {
          continue;
        }
        if (cfgRemoveReassigned) {
          if (itrack.collisionId() != atrack.bestCollisionId()) {
            continue;
          }
        }

        if constexpr (has_reco_cent<C>) {
          registry.fill(HIST("Tracks/Centrality/THnDCAxyBestRec"), itrack.pt(), itrack.eta(), collision.posZ(), atrack.bestDCAXY(), atrack.bestDCAZ(), crec);
        } else {
          registry.fill(HIST("Tracks/THnDCAxyBestRec"), itrack.pt(), itrack.eta(), collision.posZ(), atrack.bestDCAXY(), atrack.bestDCAZ());
        }

        if (itrack.has_mcParticle()) {
          auto particle = itrack.template mcParticle_as<aod::McParticles>();
          if (!isChrgParticle(particle.pdgCode())) {
            continue;
          }
          if (particle.eta() <= trackCuts.minEta || particle.eta() >= trackCuts.maxEta) {
            continue;
          }
          if (cfgUseParticleSel && !isParticleSelected(particle)) {
            continue;
          }

          const auto dcaXtruth(particle.vx() - particle.mcCollision().posX());
          const auto dcaYtruth(particle.vy() - particle.mcCollision().posY());
          const auto dcaZtruth(particle.vz() - particle.mcCollision().posZ());
          auto dcaXYtruth = std::sqrt(dcaXtruth * dcaXtruth + dcaYtruth * dcaYtruth);
          auto mcCollision = particle.template mcCollision_as<aod::McCollisions>();

          if (eventCuts.useZDiffCut) {
            if (std::abs(collision.posZ() - mcCollision.posZ()) > eventCuts.maxZvtxDiff) {
              continue;
            }
          }

          if (collision.has_mcCollision() && collision.mcCollisionId() == particle.mcCollisionId()) {
            if (!particle.isPhysicalPrimary()) { // Secondaries (weak decays and material)
              if constexpr (has_reco_cent<C>) {
                registry.fill(HIST("Tracks/Centrality/THnDCAxyBestGenSec"), particle.pt(), particle.eta(), mcCollision.posZ(), atrack.bestDCAXY(), atrack.bestDCAZ(), crec);
                registry.fill(HIST("Tracks/Centrality/THnDCAxyBestGenTruthSec"), particle.pt(), particle.eta(), mcCollision.posZ(), dcaXYtruth, dcaZtruth, crec);
              } else {
                registry.fill(HIST("Tracks/THnDCAxyBestGenSec"), particle.pt(), particle.eta(), mcCollision.posZ(), atrack.bestDCAXY(), atrack.bestDCAZ());
                registry.fill(HIST("Tracks/THnDCAxyBestGenTruthSec"), particle.pt(), particle.eta(), mcCollision.posZ(), dcaXYtruth, dcaZtruth);
              }
              if (particle.getProcess() == TMCProcess::kPDecay) { // Particles from decay
                if constexpr (has_reco_cent<C>) {
                  registry.fill(HIST("Tracks/Centrality/THnDCAxyBestGenSecWeak"), particle.pt(), particle.eta(), mcCollision.posZ(), atrack.bestDCAXY(), atrack.bestDCAZ(), crec);
                } else {
                  registry.fill(HIST("Tracks/THnDCAxyBestGenSecWeak"), particle.pt(), particle.eta(), mcCollision.posZ(), atrack.bestDCAXY(), atrack.bestDCAZ());
                }
              } else { // Particles from the material
                if constexpr (has_reco_cent<C>) {
                  registry.fill(HIST("Tracks/Centrality/THnDCAxyBestGenSecMat"), particle.pt(), particle.eta(), mcCollision.posZ(), atrack.bestDCAXY(), atrack.bestDCAZ(), crec);
                } else {
                  registry.fill(HIST("Tracks/THnDCAxyBestGenSecMat"), particle.pt(), particle.eta(), mcCollision.posZ(), atrack.bestDCAXY(), atrack.bestDCAZ());
                }
              }
            } else { // Primaries
              if constexpr (has_reco_cent<C>) {
                registry.fill(HIST("Tracks/Centrality/THnDCAxyBestGenPrim"), particle.pt(), particle.eta(), mcCollision.posZ(), atrack.bestDCAXY(), atrack.bestDCAZ(), crec);
                registry.fill(HIST("Tracks/Centrality/THnDCAxyBestGenTruthPrim"), particle.pt(), particle.eta(), mcCollision.posZ(), dcaXYtruth, dcaZtruth, crec);
              } else {
                registry.fill(HIST("Tracks/THnDCAxyBestGenPrim"), particle.pt(), particle.eta(), mcCollision.posZ(), atrack.bestDCAXY(), atrack.bestDCAZ());
                registry.fill(HIST("Tracks/THnDCAxyBestGenTruthPrim"), particle.pt(), particle.eta(), mcCollision.posZ(), dcaXYtruth, dcaZtruth);
              }
            }
          } else {                               // Wrong collision
            if (!particle.isPhysicalPrimary()) { // Secondaries (weak decays and material)
              if constexpr (has_reco_cent<C>) {
                registry.fill(HIST("Tracks/Centrality/THnDCAxyBestGenSecWrongColl"), particle.pt(), particle.eta(), mcCollision.posZ(), atrack.bestDCAXY(), atrack.bestDCAZ(), crec);
                registry.fill(HIST("Tracks/Centrality/THnDCAxyBestGenTruthSecWrongColl"), particle.pt(), particle.eta(), mcCollision.posZ(), dcaXYtruth, dcaZtruth, crec);
              } else {
                registry.fill(HIST("Tracks/THnDCAxyBestGenSecWrongColl"), particle.pt(), particle.eta(), mcCollision.posZ(), atrack.bestDCAXY(), atrack.bestDCAZ());
                registry.fill(HIST("Tracks/THnDCAxyBestGenTruthSecWrongColl"), particle.pt(), particle.eta(), mcCollision.posZ(), dcaXYtruth, dcaZtruth);
              }
            } else { // Primaries
              if constexpr (has_reco_cent<C>) {
                registry.fill(HIST("Tracks/Centrality/THnDCAxyBestGenPrimWrongColl"), particle.pt(), particle.eta(), mcCollision.posZ(), atrack.bestDCAXY(), atrack.bestDCAZ(), crec);
                registry.fill(HIST("Tracks/Centrality/THnDCAxyBestGenTruthPrimWrongColl"), particle.pt(), particle.eta(), mcCollision.posZ(), dcaXYtruth, dcaZtruth, crec);
              } else {
                registry.fill(HIST("Tracks/THnDCAxyBestGenPrimWrongColl"), particle.pt(), particle.eta(), mcCollision.posZ(), atrack.bestDCAXY(), atrack.bestDCAZ());
                registry.fill(HIST("Tracks/THnDCAxyBestGenTruthPrimWrongColl"), particle.pt(), particle.eta(), mcCollision.posZ(), dcaXYtruth, dcaZtruth);
              }
            }
          }
        } else {
          LOGP(debug, "No MC particle for ambiguous itrack, skip...");
          if constexpr (has_reco_cent<C>) {
            registry.fill(HIST("Tracks/Centrality/THnDCAxyBestRecFake"), itrack.pt(), itrack.eta(), collision.posZ(), atrack.bestDCAXY(), atrack.bestDCAZ(), crec);
          } else {
            registry.fill(HIST("Tracks/THnDCAxyBestRecFake"), itrack.pt(), itrack.eta(), collision.posZ(), atrack.bestDCAXY(), atrack.bestDCAZ());
          }
        }
      }
    }
  }

  void processDCAReassocMcInclusive(soa::Join<Colls, aod::McCollisionLabels> const& collisions,
                                    aod::McCollisions const& mccollisions,
                                    aod::McParticles const& particles,
                                    BestTracksMC const& besttracks,
                                    FiltMcMftTracks const& tracks)
  {
    processDCAReassocMc<Colls, aod::McCollisions>(collisions, mccollisions, particles, besttracks, tracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processDCAReassocMcInclusive, "Process MC DCA checks using re-association information based on BestCollisionsFwd3d table (Inclusive)", false);

  void processDCAReassocMcCentFT0C(soa::Join<CollsCentFT0C, aod::McCollisionLabels> const& collisions,
                                   aod::McCollisions const& mccollisions,
                                   aod::McParticles const& particles,
                                   BestTracksMC const& besttracks,
                                   FiltMcMftTracks const& tracks)
  {
    processDCAReassocMc<CollsCentFT0C, aod::McCollisions>(collisions, mccollisions, particles, besttracks, tracks);
  }

  PROCESS_SWITCH(DndetaMFTPbPb, processDCAReassocMcCentFT0C, "Process MC DCA checks using re-association information based on BestCollisionsFwd3d table (in FT0C centrality bins)", false);

  template <typename C>
  void processCorrelationwBestTracks(typename C::iterator const& collision, FiltMftTracks const& /*tracks*/, soa::SmallGroups<aod::BestCollisionsFwd3d> const& besttracks)
  {
    if (!isGoodEvent<false>(collision)) {
      return;
    }

    auto nBestTrks = 0;
    for (auto const& atrack : besttracks) {
      if (cfgUseTrackSel && !isBestTrackSelected<false>(atrack)) {
        continue;
      }
      auto itrack = atrack.template mfttrack_as<FiltMftTracks>();
      if (itrack.eta() < trackCuts.minEta || itrack.eta() > trackCuts.maxEta) {
        continue;
      }
      if (cfgUseTrackSel && !isTrackSelected<false>(itrack)) {
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
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<DndetaMFTPbPb>(cfgc)};
}
