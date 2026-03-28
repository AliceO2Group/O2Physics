// Copyright 2020-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file   PseudorapidityDensityMFT.cxx
/// \author Sarah Herrmann <sarah.herrmann@cern.ch>
/// \author Tulika Tripathy <tulika.tripathy@cern.ch>
/// \brief This code loops over MFT tracks and collisions and fills histograms
///        useful to compute dNdeta

#include "PWGMM/Mult/DataModel/bestCollisionTable.h"

#include "Common/Core/trackUtilities.h"
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
#include "ReconstructionDataFormats/TrackParametrization.h"

#include "Math/MatrixFunctions.h"
#include "Math/SMatrix.h"
#include "Math/SVector.h"
#include "TFile.h"
#include "TGeoGlobalMagField.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <unordered_set>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::track;
using namespace o2::aod::rctsel;

AxisSpec ptAxis = {2001, -0.005, 20.005};
AxisSpec deltazAxis = {6100, -6.1, 6.1};
AxisSpec zAxis = {3001, -30.1, 30.1};
AxisSpec phiAxis = {629, 0, o2::constants::math::TwoPI, "Rad", "phi axis"};
AxisSpec etaAxis = {18, -4.6, -1.};
// AxisSpec dcaXyAxis = {2000, -10, 10};
// AxisSpec dcaZAxis  = {2000, -10, 10};
// AxisSpec dcaXAxis  = {2000, -10, 10};
// AxisSpec dcaYAxis  = {2000, -10, 10};// previous AxisSpec dcaYAxis  = {2000, -10, 10};

AxisSpec dcaXyAxis = {6000, -30, 30};
AxisSpec dcaZAxis = {6000, -30, 30};
AxisSpec dcaXAxis = {6000, -30, 30};
AxisSpec dcaYAxis = {6000, -30, 30}; // previous AxisSpec dcaYAxis  = {2000, -10, 10};

// AxisSpec dcaXyAxis = {600, -0.15f, 0.15f};
// AxisSpec dcaZAxis  = {600, -0.15f, 0.15f};
// AxisSpec dcaXAxis  = {600, -0.15f, 0.15f};
// AxisSpec dcaYAxis  = {600, -0.15f, 0.15f};
// bin width 0.0005 cm: range [-30, 30] cm => 60/0.0005 = 120000 bins
// Keep bin width = 0.0005 cm (5 um): range [-1, 1] cm => 2.0/0.0005 = 4000 bins
// AxisSpec axisBinsDCA = {600, -0.15f, 0.15f, "#it{dca}_{xy} (cm)"};

AxisSpec centAxis = {{0, 10, 20, 30, 40, 50, 60, 70, 80, 100}};

// Vertex position axes (cm)
AxisSpec vxAxis = {200, -0.5, 0.5, "V_{x} (cm)"};
AxisSpec vyAxis = {200, -0.5, 0.5, "V_{y} (cm)"};
// Status axis for reco/truth (1=reco, 2=true)
AxisSpec recoTruthStatusAxis = {2, 0.5, 2.5, "status"};

// Delta-vertex axes (reco - true) in cm
AxisSpec deltaVxAxis = {400, -0.5, 0.5, "#DeltaV_{x} = V_{x}^{rec}-V_{x}^{true} (cm)"};
AxisSpec deltaVyAxis = {400, -0.5, 0.5, "#DeltaV_{y} = V_{y}^{rec}-V_{y}^{true} (cm)"};

static constexpr TrackSelectionFlags::flagtype TrackSelectionIts =
  TrackSelectionFlags::kITSNCls | TrackSelectionFlags::kITSChi2NDF |
  TrackSelectionFlags::kITSHits;

static constexpr TrackSelectionFlags::flagtype TrackSelectionTpc =
  TrackSelectionFlags::kTPCNCls |
  TrackSelectionFlags::kTPCCrossedRowsOverNCls |
  TrackSelectionFlags::kTPCChi2NDF;

static constexpr TrackSelectionFlags::flagtype TrackSelectionDca =
  TrackSelectionFlags::kDCAz | TrackSelectionFlags::kDCAxy;

// using MFTTracksLabeled = soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels>;
//  replace your alias with the extension included:
using FullBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;

// using MFTTracksLabeled =
//   soa::Join<o2::aod::MFTTracks,aod::BestCollisionsFwd3d,
//             /*aod::MFTTracks_001Extension, */ // exposes bestCollisionId, bestDCAXY, (and bestDCAZ if 3D)
//             aod::McMFTTrackLabels>;
using MFTTracksLabeled =
  soa::Join<o2::aod::MFTTracks, aod::MFTTrkCompColls, aod::BestCollisionsFwd3d,
            /*aod::MFTTracks_001Extension, */ // exposes bestCollisionId, bestDCAXY, (and bestDCAZ if 3D)
            aod::McMFTTrackLabels>;
using MFTTracksLabeled2d =
  soa::Join<o2::aod::MFTTracks, aod::BestCollisionsFwd,
            aod::McMFTTrackLabels>;

using MFTTracksLabeledOrg =
  soa::Join<o2::aod::MFTTracks, aod::MFTTrkCompColls,
            /*aod::MFTTracks_001Extension, */ // exposes bestCollisionId, bestDCAXY, (and bestDCAZ if 3D)
            aod::McMFTTrackLabels>;
// using McCollisionsWithExtra = soa::Join<aod::McCollisions, aod::McCollsExtra>;
using McCollisionsWithExtra = o2::soa::Join<o2::aod::McCollisions, o2::aod::McCollsExtra>;

using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<double, 5>;
struct PseudorapidityDensityMFT {
  SliceCache cache;
  Preslice<aod::MFTTracks> perCol = o2::aod::fwdtrack::collisionId;
  Preslice<aod::McParticles> perMcCol = aod::mcparticle::mcCollisionId;
  Preslice<aod::Tracks> perColCentral = aod::track::collisionId;

  Service<o2::framework::O2DatabasePDG> pdg;

  // --- CCDB magnetic field (needed for propagateToDCAhelix in this device) ---
  Service<o2::ccdb::BasicCCDBManager> ccdbMgr;
  Configurable<std::string> ccdburlMag{"ccdburlMag", "http://alice-ccdb.cern.ch",
                                       "CCDB url for GRPMagField"};
  Configurable<std::string> grpmagPathMag{"grpmagPathMag", "GLO/Config/GRPMagField",
                                          "CCDB path for GRPMagField"};

  int magRunNumber = -1;
  float bzMFT = 0.f;
  o2::parameters::GRPMagField* grpmag = nullptr;
  static constexpr double CenterMft[3] = {0., 0., -61.4};

  enum class GenRecoCutBin : int {
    AllRecoCollisions = 1,
    UseContBestCollisionIndex,
    HasMcCollision,
    IsTriggerTVX,
    NoTimeFrameBorder,
    NoITSROFrameBorder,
    NoSameBunchPileup,
    GoodZvtxFT0vsPV,
    NoCollInRofStandard,
    NoCollInRofStrict,
    NoCollInTimeRangeStandard,
    NoCollInTimeRangeStrict,
    NoHighMultCollInPrevRof,
    RctMFT,
    VzWindow,
    InelGt0
  };

  enum class TrackLabelSummaryBin : int {
    AllTracks = 1,
    NoMcLabel,
    FakeTracks,
    TrueTracks,
    PrimaryTracks,
    SecondaryTracks
  };

  enum class GenRecoTimeComTrackMode : int {
    AllNonOrphan = 0,
    NonOrphanNonAmbiguous,
    NonOrphanAmbiguous
  };

  enum class RightWrongBin : int {
    Right = 1,
    Wrong,
    Neither,
    Both
  };
  enum class EventSelectionBin : int {
    All = 1,
    Vz,
    VzItsRof,
    VzSelected,
    Sel8VzInelGt0,
    SelInelInelFwdGt0,
    Rejected,
    GoodBCs,
    BCsWithCollisions,
    BCsWithPileupSplitting,
    PerCollisionSampleGt0,
    MidtracksAndPerCollisionSampleGt0
  };
  enum class HashTableRowCountsBin : int {
    RowsSaved = 1,
    UniqueRecoColsSaved,
    UniqueBestRecoCols
  };
  enum class WrongVertexRecoExistsBin : int {
    RecoOfTrueExists = 1,
    RecoOfTrueMissing
  };
  enum class BoolBin : int {
    No = 0,
    Yes = 1
  };
  enum class SingleCountBin : int {
    Count = 1
  };
  enum class EventEfficiencyBin : int {
    Generated = 1,
    GeneratedInelGt0,
    Reconstructed,
    Selected,
    SelectedInelGt0
  };
  enum class CentralitySelectionBin : int {
    All = 1,
    Selected,
    Rejected
  };
  static constexpr float ForwardEtaMax = -2.0f;
  static constexpr float ForwardEtaMin = -3.9f;

  static constexpr float PhiVetoLow = 0.02f;
  static constexpr float PhiVetoPiMin = 3.10f;
  static constexpr float PhiVetoPiMax = 3.23f;
  static constexpr float PhiVetoHigh = 6.21f;

  static constexpr float NdfScale = 2.0f;
  static constexpr float NdfOffset = 5.0f;
  static constexpr float MinNdf = 1.0f;

  template <typename TrackT>
  static float getTrackNdf(TrackT const& track)
  {
    return std::max(NdfScale * track.nClusters() - NdfOffset, MinNdf);
  }
  static constexpr int NoCompatibleCollisions = 0;
  static constexpr int SingleCompatibleCollision = 1;

  static constexpr int OrphanAmbDegree = 0;
  static constexpr int NonAmbiguousAmbDegree = 1;

  static constexpr int ChargeUnitTimesThree = 3;

  void initMagField(FullBCs::iterator const& bc)
  {
    if (magRunNumber == bc.runNumber()) {
      return;
    }

    grpmag = ccdbMgr->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPathMag.value, bc.timestamp());
    if (!grpmag) {
      LOGF(warning, "GRPMagField not found in CCDB for ts=%lld", (long long)bc.timestamp());
      bzMFT = 0.f;
      magRunNumber = bc.runNumber();
      return;
    }

    // This sets TGeoGlobalMagField internally
    o2::base::Propagator::initFieldFromGRP(grpmag);
    magRunNumber = bc.runNumber();

    auto* field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
    if (field) {
      bzMFT = field->getBz(CenterMft);
      LOGP(info, "Initialized magnetic field for run {}: bzMFT={} kG", magRunNumber, bzMFT);
    } else {
      LOGF(warning, "TGeoGlobalMagField has no field even after initFieldFromGRP; bzMFT=0");
      bzMFT = 0.f;
    }
  }

  RCTFlagsChecker rctChecker{"CBT"};
  RCTFlagsChecker myChecker{kFT0Bad, kITSBad, kTPCBadTracking, kTPCBadPID, kMFTBad};
  Configurable<int> maxGenRecoEvents{"maxGenRecoEvents", -1,
                                     "Maximum number of MC collisions to process in processGenReco (-1 = all)"};

  int nProcessedGenReco = 0;

  Configurable<float> estimatorEta{"estimatorEta", 1.0,
                                   "eta range for INEL>0 sample definition"};

  Configurable<bool> useEvSel{"useEvSel", true, "use event selection"};
  Configurable<bool> useNoSameBunchPileup{"useNoSameBunchPileup", true, "reject collisions in case of pileup with another collision in the same foundBC"};
  Configurable<bool> useNoCollInRofStandard{"useNoCollInRofStandard", true, "Require evsel::kNoCollInRofStrict in processGenReco"};
  Configurable<bool> useNoCollInRofStrict{"useNoCollInRofStrict", true, "Require evsel::kNoCollInRofStrict in processGenReco"};
  Configurable<bool> useNoCollInTimeRangeStrict{"useNoCollInTimeRangeStrict", true, "Require evsel::kNoCollInTimeRangeStrict in processGenReco"};
  Configurable<bool> useNoCollInTimeRangeStandard{"useNoCollInTimeRangeStandard", true, "Require evsel::kNoCollInTimeRangeStandard in processGenReco"};
  Configurable<bool> useNoHighMultCollInPrevRof{"useNoHighMultCollInPrevRof", true, "Require evsel::kNoHighMultCollInPrevRof in processGenReco"};
  Configurable<bool> useGoodZvtxFT0vsPV{"useGoodZvtxFT0vsPV", true, "Require evsel::kIsGoodZvtxFT0vsPV in processGenReco"};
  Configurable<bool> useRctMFT{"useRctMFT", true, "Apply RCT runlist flags for MFT"};

  Configurable<bool> disableITSROFCut{"disableITSROFCut", false, "Disable ITS ROF cut for event selection"};

  ConfigurableAxis multBinning{"multBinning", {701, -0.5, 700.5}, ""};
  ConfigurableAxis etaBinning{"etaBinning", {36, -4.6, -1.}, ""};
  Configurable<bool> useZDiffCut{"useZDiffCut", true, "use Z difference cut"};
  Configurable<float> maxZDiff{"maxZDiff", 1.0f, "max allowed Z difference for reconstructed collisions (cm)"};

  Configurable<bool> usePhiCut{"usePhiCut", true, "use azimuthal angle cut"};
  Configurable<bool> useDCAxyCut{"useDCAxyCut", false, "use DCAxy cut"};
  Configurable<bool> useCont{"useCont", false, "No of contributors cut"};

  Configurable<bool> usePtCut{"usePtCut", false, "use Pt cut"};

  Configurable<bool> useDCAzCut{"useDCAzCut", false, "use DCAz cut"};

  Configurable<float> cfgPhiCut{"cfgPhiCut", 0.1f,
                                "Cut on azimuthal angle of MFT tracks"};
  Configurable<float> cfgPhiCut1{"cfgPhiCut1", 0.0f,
                                 "low Cut on azimuthal angle of MFT tracks"};
  Configurable<float> cfgPhiCut2{"cfgPhiCut2", 6.3f,
                                 "high Cut on azimuthal angle of MFT tracks"};
  Configurable<float> cfgVzCut1{"cfgVzCut1", -30.0f,
                                "Cut1 on vertex position of MFT tracks"};
  Configurable<float> cfgVzCut2{"cfgVzCut2", 30.0f,
                                "Cut2 on vertex position of MFT tracks"};
  Configurable<float> cfgnCluster{"cfgnCluster", 5.0f,
                                  "Cut on no of clusters per MFT track"};
  Configurable<float> cfgnEta1{"cfgnEta1", -4.5f,
                               "Cut on eta1"};
  Configurable<float> cfgnEta2{"cfgnEta2", -1.0f,
                               "Cut on eta1"};
  Configurable<float> cfgnPt{"cfgnPt", 10.0f,
                             "Cut on Pt"};
  Configurable<float> cfgChi2NDFMax{"cfgChi2NDFMax", 2000.0f, "Max allowed chi2/NDF for MFT tracks"};
  Configurable<float> maxDCAxy{"maxDCAxy", 2.0f, "Cut on dcaXY"};
  Configurable<float> maxDCAz{"maxDCAz", 2.0f, "Cut on dcaZ"};
  Configurable<bool> useLostByCutVeto{"useLostByCutVeto", true, "Reject tracks with lostNoRecoOfTrue or lostRecoExistsButNotCompatible"};

  Configurable<int> cfgGenRecoTimeComTrackMode{"cfgGenRecoTimeComTrackMode",
                                               static_cast<int>(GenRecoTimeComTrackMode::AllNonOrphan),
                                               "processGenRecoTimeCom track mode: AllNonOrphan=0, NonOrphanNonAmbiguous=1, NonOrphanAmbiguous=2"};

  HistogramRegistry registry{
    "registry",
    {
      {"TracksEtaZvtx",
       "; #eta; #it{z}_{vtx} (cm); tracks",
       {HistType::kTH2F, {etaBinning, zAxis}}}, //
      {"Tracks/EtaZvtx_gt0",
       "; #eta; #it{z}_{vtx} (cm); tracks",
       {HistType::kTH2F, {etaBinning, zAxis}}}, //
      {"TracksPhiEta",
       "; #varphi; #eta; tracks",
       {HistType::kTH2F, {phiAxis, etaBinning}}}, //
      {"TracksPhiZvtx",
       "; #varphi; #it{z}_{vtx} (cm); tracks",
       {HistType::kTH2F, {phiAxis, zAxis}}}, //
      {"TracksPtEta",
       " ; p_{T} (GeV/c); #eta",
       {HistType::kTH2F, {ptAxis, etaBinning}}}, //
      {"EventSelection",
       ";status;events",
       {HistType::kTH1F, {{15, 0.5, 15.5}}}},
      {"EventCounts",
       ";status;events",
       {HistType::kTH1F, {{2, 0.5, 2.5}}}},
    }};

  void init(InitContext&)
  {
    ccdbMgr->setURL(ccdburlMag.value); // or ccdburlMag.value (depending on your Configurable)
    ccdbMgr->setCaching(true);
    ccdbMgr->setLocalObjectValidityChecking();

    if (static_cast<int>(doprocessMult) +
          static_cast<int>(doprocessMultReassoc) +
          static_cast<int>(doprocessMultReassoc3d) +
          static_cast<int>(doprocessCountingCentrality) >
        1) {
      LOGP(fatal,
           "Exactly one process function between processMult, "
           "processMultReassoc, processMultReassoc3d and processCountingCentrality should be "
           "enabled!");
    }
    AxisSpec multAxis = {multBinning, "N_{trk}"};
    auto hstat = registry.get<TH1>(HIST("EventSelection"));
    auto* x = hstat->GetXaxis();
    x->SetBinLabel(static_cast<int>(EventSelectionBin::All), "All");
    x->SetBinLabel(static_cast<int>(EventSelectionBin::Vz), "Vz");
    x->SetBinLabel(static_cast<int>(EventSelectionBin::VzItsRof), "Vz+ITSRof");
    x->SetBinLabel(static_cast<int>(EventSelectionBin::VzSelected), "Vz+Selected");
    x->SetBinLabel(static_cast<int>(EventSelectionBin::Sel8VzInelGt0), "Sel8+Vz+INEL>0");
    x->SetBinLabel(static_cast<int>(EventSelectionBin::SelInelInelFwdGt0), "Sel INEL,INEL_fwd>0");
    x->SetBinLabel(static_cast<int>(EventSelectionBin::Rejected), "Rejected");
    x->SetBinLabel(static_cast<int>(EventSelectionBin::GoodBCs), "Good BCs");
    x->SetBinLabel(static_cast<int>(EventSelectionBin::BCsWithCollisions), "BCs with collisions");
    x->SetBinLabel(static_cast<int>(EventSelectionBin::BCsWithPileupSplitting), "BCs with pile-up/splitting");
    x->SetBinLabel(static_cast<int>(EventSelectionBin::PerCollisionSampleGt0), "percollisionSample>0");
    x->SetBinLabel(static_cast<int>(EventSelectionBin::MidtracksAndPerCollisionSampleGt0), "midtracks+percollisionSample>0");
    registry.add({"EventsNtrkZvtx",
                  "; N_{trk}; #it{z}_{vtx} (cm); events",
                  {HistType::kTH2F, {multAxis, zAxis}}});
    registry.add({"EventsNtrkZvtx_gt0",
                  "; N_{trk}; #it{z}_{vtx} (cm); events",
                  {HistType::kTH2F, {multAxis, zAxis}}});
    registry.add({"Tracks/2Danalysis/EventsNtrkZvtx_all",
                  "; N_{trk}; #it{z}_{vtx} (cm); events",
                  {HistType::kTH2F, {multAxis, zAxis}}});
    registry.add({"Tracks/2Danalysis/EventsNtrkZvtx_sel8",
                  "; N_{trk}; #it{z}_{vtx} (cm); events",
                  {HistType::kTH2F, {multAxis, zAxis}}});
    registry.add({"Tracks/2Danalysis/EventsNtrkZvtx_sel8_inelgt0",
                  "; N_{trk}; #it{z}_{vtx} (cm); events",
                  {HistType::kTH2F, {multAxis, zAxis}}});
    registry.add({"Tracks/2Danalysis/EventsNtrkZvtx_sel8_inelfwdgt0",
                  "; N_{trk}; #it{z}_{vtx} (cm); events",
                  {HistType::kTH2F, {multAxis, zAxis}}});
    registry.add({"Tracks/Control/DCAXY",
                  " ; DCA_{XY} (cm)",
                  {HistType::kTH1F, {dcaXyAxis}}});

    if (doprocessGenReco || doprocessGenRecoTimeCom) {
      registry.add({"EventsRecoCuts_GenReco",
                    ";cut;events",
                    {HistType::kTH1F, {{16, 0.5, 16.5}}}});
      {
        auto h = registry.get<TH1>(HIST("EventsRecoCuts_GenReco"));
        auto* x = h->GetXaxis();
        x->SetBinLabel(static_cast<int>(GenRecoCutBin::AllRecoCollisions), "All reco collisions (loop entry)");
        x->SetBinLabel(static_cast<int>(GenRecoCutBin::UseContBestCollisionIndex), "useContBestcollisionIndex");
        x->SetBinLabel(static_cast<int>(GenRecoCutBin::HasMcCollision), "has_mcCollision()");
        x->SetBinLabel(static_cast<int>(GenRecoCutBin::IsTriggerTVX), "kIsTriggerTVX (if useEvSel)");
        x->SetBinLabel(static_cast<int>(GenRecoCutBin::NoTimeFrameBorder), "kNoTimeFrameBorder (if useEvSel)");
        x->SetBinLabel(static_cast<int>(GenRecoCutBin::NoITSROFrameBorder), "kNoITSROFrameBorder (if useEvSel)");
        x->SetBinLabel(static_cast<int>(GenRecoCutBin::NoSameBunchPileup), "kNoSameBunchPileup");
        x->SetBinLabel(static_cast<int>(GenRecoCutBin::GoodZvtxFT0vsPV), "kIsGoodZvtxFT0vsPV");
        x->SetBinLabel(static_cast<int>(GenRecoCutBin::NoCollInRofStandard), "kNoCollInRofStandard (cfg)");
        x->SetBinLabel(static_cast<int>(GenRecoCutBin::NoCollInRofStrict), "kNoCollInRofStrict (cfg)");
        x->SetBinLabel(static_cast<int>(GenRecoCutBin::NoCollInTimeRangeStandard), "kNoCollInTimeRangeStandard (cfg)");
        x->SetBinLabel(static_cast<int>(GenRecoCutBin::NoCollInTimeRangeStrict), "kNoCollInTimeRangeStrict (cfg)");
        x->SetBinLabel(static_cast<int>(GenRecoCutBin::NoHighMultCollInPrevRof), "kNoHighMultCollInPrevRof (cfg)");
        x->SetBinLabel(static_cast<int>(GenRecoCutBin::RctMFT), "myChecker (cfg)");
        x->SetBinLabel(static_cast<int>(GenRecoCutBin::VzWindow), "Vz window");
        x->SetBinLabel(static_cast<int>(GenRecoCutBin::InelGt0), "INEL>0 (midtracks>0)");
        // x->SetBinLabel(11, "rctChecker (cfg)");

        // x->SetBinLabel(15, "Passed all event cuts");
      }

      AxisSpec trackLabelAxis = {6, 0.5, 6.5, "Track label category"};
      registry.add({"Purity/TrackLabelSummary",
                    ";Track label category;Counts",
                    {HistType::kTH1D, {trackLabelAxis}}});
      {
        auto h = registry.get<TH1>(HIST("Purity/TrackLabelSummary"));
        auto* x = h->GetXaxis();
        x->SetBinLabel(static_cast<int>(TrackLabelSummaryBin::AllTracks), "All tracks");
        x->SetBinLabel(static_cast<int>(TrackLabelSummaryBin::NoMcLabel), "No MC label");
        x->SetBinLabel(static_cast<int>(TrackLabelSummaryBin::FakeTracks), "Fake tracks");
        x->SetBinLabel(static_cast<int>(TrackLabelSummaryBin::TrueTracks), "True tracks");
        x->SetBinLabel(static_cast<int>(TrackLabelSummaryBin::PrimaryTracks), "Primary tracks");
        x->SetBinLabel(static_cast<int>(TrackLabelSummaryBin::SecondaryTracks), "Secondary tracks");
      }
      registry.add({"Purity/TrackEtaCategorySparse",
                    ";#eta_{all};#eta_{noMC};#eta_{fake};#eta_{true};#eta_{primary};#eta_{secondary}",
                    {HistType::kTHnSparseF, {etaBinning, etaBinning, etaBinning, etaBinning, etaBinning, etaBinning}}});

      registry.add({"Purity/WrongVertexRecoExists",
                    ";status;Counts",
                    {HistType::kTH1F, {{2, 0.5, 2.5}}}});
      {
        auto h = registry.get<TH1>(HIST("Purity/WrongVertexRecoExists"));
        if (h) {
          h->GetXaxis()->SetBinLabel(static_cast<int>(WrongVertexRecoExistsBin::RecoOfTrueExists), "Reco of true exists");
          h->GetXaxis()->SetBinLabel(static_cast<int>(WrongVertexRecoExistsBin::RecoOfTrueMissing), "Reco of true missing");
        }
      }

      registry.add({"Purity/HashTableRowCounts",
                    ";status;counts",
                    {HistType::kTH1F, {{3, 0.5, 3.5}}}});
      auto hHashTableRowCounts = registry.get<TH1>(HIST("Purity/HashTableRowCounts"));
      auto* xHash = hHashTableRowCounts->GetXaxis();
      xHash->SetBinLabel(static_cast<int>(HashTableRowCountsBin::RowsSaved), "rows saved");
      xHash->SetBinLabel(static_cast<int>(HashTableRowCountsBin::UniqueRecoColsSaved), "unique recoCol saved");
      xHash->SetBinLabel(static_cast<int>(HashTableRowCountsBin::UniqueBestRecoCols), "unique bestRecoCol");
      registry.add({"Purity/THnDCAChosenVsRight_Wrong",
                    ";#eta;DCA_{xy}^{chosen} (cm);DCA_{xy}^{calculated} (cm);DCA_{xy}^{right} (cm);DCA_{z}^{chosen} (cm);DCA_{z}^{calculated} (cm);DCA_{z}^{right} (cm)",
                    {HistType::kTHnSparseF, {etaBinning, dcaXyAxis, dcaXyAxis, dcaXyAxis, dcaZAxis, dcaZAxis, dcaZAxis}}});
      registry.add({"Purity/THnDCAChosenVsRight_Right",
                    ";#eta;DCA_{xy}^{chosen} (cm);DCA_{xy}^{calculated} (cm);DCA_{xy}^{right} (cm);DCA_{z}^{chosen} (cm);DCA_{z}^{calculated} (cm);DCA_{z}^{right} (cm)",
                    {HistType::kTHnSparseF, {etaBinning, dcaXyAxis, dcaXyAxis, dcaXyAxis, dcaZAxis, dcaZAxis, dcaZAxis}}});
      registry.add("Purity/RecoOfTrueExists",
                   "Any reco collision exists for track true MC collision id;exists (0=no,1=yes);tracks",
                   kTH1F, {{2, -0.5, 1.5}});
      registry.add("Purity/RecoOfTrueInCompatible",
                   "Reco collision(s) of true MC event present in track compatible collisions;inCompatible (0=no,1=yes);tracks",
                   kTH1F, {{2, -0.5, 1.5}});
      registry.add("Purity/RecoOfTrueExistsR",
                   "Any reco collision exists for track true MC collision id;exists (0=no,1=yes);tracks",
                   kTH1F, {{2, -0.5, 1.5}});
      registry.add("Purity/RecoOfTrueInCompatibleR",
                   "Reco collision(s) of true MC event present in track compatible collisions;inCompatible (0=no,1=yes);tracks",
                   kTH1F, {{2, -0.5, 1.5}});

      registry.add("Purity/RecoOfTrueExistsW",
                   "Any reco collision exists for track true MC collision id;exists (0=no,1=yes);tracks",
                   kTH1F, {{2, -0.5, 1.5}});
      registry.add("Purity/RecoOfTrueInCompatibleW",
                   "Reco collision(s) of true MC event present in track compatible collisions;inCompatible (0=no,1=yes);tracks",
                   kTH1F, {{2, -0.5, 1.5}});

      registry.add("Purity/hCorrectRecoIDinTheListR",
                   "Any reco collision exists for track true MC collision id;exists (0=no,1=yes);tracks",
                   kTH1F, {{2, -0.5, 1.5}});
      registry.add("Purity/hCorrectRecoIDinTheListW",
                   "Any reco collision exists for track true MC collision id;exists (0=no,1=yes);tracks",
                   kTH1F, {{2, -0.5, 1.5}});

      // P(Nch): number of selected MFT tracks per accepted reco collision (after all event+track cuts)
      // Tracks lost because of the OR cut ( !exists || !inCompatible )
      registry.add("Purity/LostByBoth",
                   "Tracks rejected by (!recoOfTrueExists || !recoOfTrueInCompatible);status (0=kept,1=lost);tracks",
                   kTH1F, {{2, -0.5, 1.5}});
      // Number of ITS-TPC contributors to the reconstructed collision (PV contributors)
      registry.add({"Purity/reco/CollisionNumContrib",
                    ";N_{contrib} ( PV contributors);collisions",
                    {HistType::kTH1F, {{3001, -0.5, 3000.5}}}});

      // Tracks that were WRONG (by your definition) BEFORE applying the cut
      registry.add("Purity/WrongBeforeRecoOfTrueCut",
                   "Tracks classified wrong BEFORE applying recoOfTrue cut;wrong (0=no,1=yes);tracks",
                   kTH1F, {{2, -0.5, 1.5}});

      // Optional but very useful: intersection (lost AND wrong)

      registry.add({"Purity/reco/PNchMFT_afterCuts",
                    ";N_{trk}^{MFT} (selected);events",
                    {HistType::kTH1F, {multAxis}}});
      registry.add({"Purity/DCAyVsDCAx_Right",
                    ";DCA_{x} (cm);DCA_{y} (cm)",
                    {HistType::kTH2F, {dcaXAxis, dcaYAxis}}});
      registry.add({"Purity/reco/woOrp/All",
                    ";bin;counts",
                    {HistType::kTH1F, {{1, 0.5, 1.5}}}});
      registry.add({"Purity/reco/woOrp/AllEta",
                    ";#eta;counts",
                    {HistType::kTH1F, {etaBinning}}});
      registry.add({"Purity/SelectedAfterDCAxy/PrimaryAll",
                    ";bin;counts",
                    {HistType::kTH1F, {{1, 0.5, 1.5}}}});
      registry.add({"Purity/SelectedAfterDCAxy/PrimaryAllEta",
                    ";#eta;counts",
                    {HistType::kTH1F, {etaBinning}}});
      registry.add({"TracksToPartPhiEta",
                    "; #varphi; #eta; tracks",
                    {HistType::kTH2F, {phiAxis, etaBinning}}});
      registry.add({"EventsReco",
                    "; status; events",
                    {HistType::kTH1F, {{7, 0.5, 7.5}}}});

      // Additional histograms for processGenReco (w/o orphan), grouped under Purity/reco/
      registry.add({"Purity/reco/weakStrange/SelectedTracksEta",
                    "; #eta; selected reco tracks from weak strange decays",
                    {HistType::kTH1F, {etaBinning}}});
      registry.add({"Purity/reco/weakStrange/SelectedTracksEtaZvtx",
                    "; #eta; #it{z}_{vtx} (cm); selected reco tracks from weak strange decays",
                    {HistType::kTH2F, {etaBinning, zAxis}}});
      registry.add({"Purity/reco/SelectedTracksEta",
                    "; #eta; selected reco tracks",
                    {HistType::kTH1F, {etaBinning}}});
      registry.add({"Purity/reco/woOrp/nTrk",
                    " ; N_{Trk}^{all}",
                    {HistType::kTH1F, {{701, -0.5, 700.5}}}});

      registry.add({"Purity/reco/woOrp/woOrpTracksEtaZvtx",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {etaBinning, zAxis}}});
      registry.add({"Purity/reco/woOrp/woOrpTracksPtZvtx",
                    "; p_{T} (GeV/c); #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {ptAxis, zAxis}}});
      registry.add({"Purity/reco/woOrp/woOrpPtZvtx_gt0",
                    "; p_{T} (GeV/c); #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {ptAxis, zAxis}}});

      registry.add({"Purity/reco/woOrp/woOrpEtaZvtx_gt0",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {etaBinning, zAxis}}});

      registry.add({"Purity/reco/woOrp/woOrpTracksDCAxyZvtx_gt0",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {dcaXyAxis, zAxis}}});

      registry.add({"Purity/reco/woOrp/woOrpTracksDCAzZvtx_gt0",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {dcaZAxis, zAxis}}});

      registry.add({"Purity/reco/woOrp/woOrpTracksPhiEta",
                    "; #varphi; #eta; tracks",
                    {HistType::kTH2F, {phiAxis, etaBinning}}});

      // --- Category-wise woOrp histograms (after track cuts & woOrp) ---
      // Fake tracks (no MC particle)
      registry.add({"Purity/reco/woOrp_fake/woOrpEtaZvtx_gt0",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {etaBinning, zAxis}}});
      registry.add({"Purity/reco/woOrp_fake/woOrpTracksEtaZvtx",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {etaBinning, zAxis}}});
      registry.add({"Purity/reco/woOrp_fake/woOrpTracksPtZvtx",
                    "; p_{T} (GeV/c); #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {ptAxis, zAxis}}});
      registry.add({"Purity/reco/woOrp_fake/woOrpPtZvtx_gt0",
                    "; p_{T} (GeV/c); #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {ptAxis, zAxis}}});
      registry.add({"Purity/reco/woOrp_fake/woOrpTracksPhiEta",
                    "; #varphi; #eta; tracks",
                    {HistType::kTH2F, {phiAxis, etaBinning}}});

      // Tracks that have an MC particle (matched)
      registry.add({"Purity/reco/woOrp_hasMC/woOrpEtaZvtx_gt0",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {etaBinning, zAxis}}});
      registry.add({"Purity/reco/woOrp_hasMC/woOrpTracksEtaZvtx",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {etaBinning, zAxis}}});
      registry.add({"Purity/reco/woOrp_hasMC/woOrpTracksPtZvtx",
                    "; p_{T} (GeV/c); #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {ptAxis, zAxis}}});
      registry.add({"Purity/reco/woOrp_hasMC/woOrpPtZvtx_gt0",
                    "; p_{T} (GeV/c); #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {ptAxis, zAxis}}});
      registry.add({"Purity/reco/woOrp_hasMC/woOrpTracksPhiEta",
                    "; #varphi; #eta; tracks",
                    {HistType::kTH2F, {phiAxis, etaBinning}}});

      // Secondary tracks (has MC but not primary)
      registry.add({"Purity/reco/woOrp_secondary/woOrpEtaZvtx_gt0",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {etaBinning, zAxis}}});
      registry.add({"Purity/reco/woOrp_secondary/woOrpTracksEtaZvtx",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {etaBinning, zAxis}}});
      registry.add({"Purity/reco/woOrp_secondary/woOrpTracksPtZvtx",
                    "; p_{T} (GeV/c); #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {ptAxis, zAxis}}});
      registry.add({"Purity/reco/woOrp_secondary/woOrpPtZvtx_gt0",
                    "; p_{T} (GeV/c); #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {ptAxis, zAxis}}});
      registry.add({"Purity/reco/woOrp_secondary/woOrpTracksPhiEta",
                    "; #varphi; #eta; tracks",
                    {HistType::kTH2F, {phiAxis, etaBinning}}});

      // Primary tracks (has MC and is physical primary)
      registry.add({"Purity/reco/woOrp_primary/woOrpEtaZvtx_gt0",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {etaBinning, zAxis}}});
      registry.add({"Purity/reco/woOrp_primary/woOrpTracksEtaZvtx",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {etaBinning, zAxis}}});
      registry.add({"Purity/reco/woOrp_primary/woOrpTracksPtZvtx",
                    "; p_{T} (GeV/c); #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {ptAxis, zAxis}}});
      registry.add({"Purity/reco/woOrp_primary/woOrpPtZvtx_gt0",
                    "; p_{T} (GeV/c); #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {ptAxis, zAxis}}});
      registry.add({"Purity/reco/woOrp_primary/woOrpTracksPhiEta",
                    "; #varphi; #eta; tracks",
                    {HistType::kTH2F, {phiAxis, etaBinning}}});

      // MC-side (generator truth) histograms for purity in processGenReco
      registry.add({"Purity/mc/PrimaryAll",
                    ";bin;counts",
                    {HistType::kTH1F, {{1, 0.5, 1.5}}}});
      registry.add({"Purity/mc/PrimaryAllEta",
                    ";#eta;counts",
                    {HistType::kTH1F, {etaBinning}}});
      registry.add({"Purity/mc/PrimaryTracksEtaZvtx",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {etaBinning, zAxis}}});
      registry.add({"Purity/mc/PrimaryTracksEtaZvtx_gt0",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {etaBinning, zAxis}}});
      registry.add({"Purity/mc/PrimaryTracksDCAxyZvtx_gt0",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {dcaXyAxis, zAxis}}});
      registry.add({"Purity/mc/PrimaryTracksDCAzZvtx_gt0",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {dcaZAxis, zAxis}}});

      registry.add({"Purity/mc/PrimaryTracksPtZvtx_gt0",
                    "; p_{T} (GeV/c); #it{z}_{vtx} (cm); primary tracks",
                    {HistType::kTH2F, {ptAxis, zAxis}}});
      registry.add({"Purity/mc/PrimaryTracksPhiEta",
                    "; #varphi; #eta; tracks",
                    {HistType::kTH2F, {phiAxis, etaBinning}}});
      // --- MC-side (NO-DCA) truth histograms for pre-DCA accounting ---
      registry.add({"Purity/mc_noDCA/PrimaryAll",
                    ";bin;counts",
                    {HistType::kTH1F, {{1, 0.5, 1.5}}}});
      registry.add({"Purity/mc_noDCA/PrimaryAllEta",
                    ";#eta;counts",
                    {HistType::kTH1F, {etaBinning}}});
      registry.add({"Purity/mc_noDCA/PrimaryTracksEtaZvtx",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {etaBinning, zAxis}}});
      registry.add({"Purity/mc_noDCA/PrimaryTracksEtaZvtx_gt0",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {etaBinning, zAxis}}});
      registry.add({"Purity/mc_noDCA/PrimaryTracksPtZvtx_gt0",
                    "; p_{T} (GeV/c); #it{z}_{vtx} (cm); primary tracks",
                    {HistType::kTH2F, {ptAxis, zAxis}}});
      registry.add({"Purity/mc_noDCA/PrimaryTracksPhiEta",
                    "; #varphi; #eta; tracks",
                    {HistType::kTH2F, {phiAxis, etaBinning}}});

      // Book keeping for pre-DCA selection counters (for symmetry with AfterDCAxy)
      registry.add({"Purity/SelectedBeforeDCAxy/PrimaryAll",
                    ";bin;counts",
                    {HistType::kTH1F, {{1, 0.5, 1.5}}}});
      registry.add({"Purity/SelectedBeforeDCAxy/PrimaryAllEta",
                    ";#eta;counts",
                    {HistType::kTH1F, {etaBinning}}});
      // --- Fake-track counters (reco side after DCA selections) ---
      registry.add({"Purity/Fakes/All",
                    ";bin;counts",
                    {HistType::kTH1F, {{1, 0.5, 1.5}}}});
      registry.add({"Purity/Fakes/AllEta",
                    ";#eta_{reco};counts",
                    {HistType::kTH1F, {etaBinning}}});
      registry.add({"Purity/Fakes/TracksEtaZvtx",
                    "; #eta_{reco}; #it{z}_{vtx}^{rec} (cm); tracks",
                    {HistType::kTH2F, {etaBinning, zAxis}}});
      registry.add({"Purity/Fakes/TracksPhiEta",
                    "; #varphi_{reco}; #eta_{reco}; tracks",
                    {HistType::kTH2F, {phiAxis, etaBinning}}});
      // Fake reason breakdown: 1-INEL>0 true (vz), 2-Reco sel (woOrp), 3-has_mcParticle, 4-isPrimary, 5-rightVertex

      // --- Purity calculation histograms (as profiles: purity = <isPrimary>) ---
      registry.add({"Purity/PurityOverall",
                    ";bin;purity",
                    {HistType::kTProfile, {{1, 0.5, 1.5}}}});
      registry.add({"Purity/PurityVsEta",
                    ";#eta;purity",
                    {HistType::kTProfile, {etaBinning}}});
      registry.add({"RightWrong",
                    ";category;counts",
                    {HistType::kTH1F, {{4, 0.5, 4.5}}}});
      auto hrw = registry.get<TH1>(HIST("RightWrong"));
      hrw->GetXaxis()->SetBinLabel(static_cast<int>(RightWrongBin::Right), "right");
      hrw->GetXaxis()->SetBinLabel(static_cast<int>(RightWrongBin::Wrong), "wrong");
      hrw->GetXaxis()->SetBinLabel(static_cast<int>(RightWrongBin::Neither), "neither");
      hrw->GetXaxis()->SetBinLabel(static_cast<int>(RightWrongBin::Both), "both");

      registry.add({"Purity/RightWrongLater",
                    ";category;counts",
                    {HistType::kTH1F, {{4, 0.5, 4.5}}}});
      auto hrw1 = registry.get<TH1>(HIST("Purity/RightWrongLater"));
      hrw1->GetXaxis()->SetBinLabel(static_cast<int>(RightWrongBin::Right), "right");
      hrw1->GetXaxis()->SetBinLabel(static_cast<int>(RightWrongBin::Wrong), "wrong");
      hrw1->GetXaxis()->SetBinLabel(static_cast<int>(RightWrongBin::Neither), "neither");
      hrw1->GetXaxis()->SetBinLabel(static_cast<int>(RightWrongBin::Both), "both");

      // Vertex-position difference for wrong-vertex associations (reco - true MC)
      registry.add({"deltaVZ_fromReco",
                    ";#Delta z_{vtx}^{reco-true} (cm);tracks",
                    {HistType::kTH1F, {zAxis}}});
      registry.add({"deltaVZ_fromTrue",
                    ";#Delta z_{vtx}^{reco-true} (cm);tracks",
                    {HistType::kTH1F, {zAxis}}});
      registry.add({"Purity/DeltaXWrong",
                    ";#Delta x_{vtx}^{reco-true} (cm);tracks",
                    {HistType::kTH1F, {dcaZAxis}}});
      registry.add({"Purity/DeltaYWrong",
                    ";#Delta y_{vtx}^{reco-true} (cm);tracks",
                    {HistType::kTH1F, {dcaZAxis}}});
      registry.add({"Purity/DeltaZWrong",
                    ";#Delta z_{vtx}^{reco-true} (cm);tracks",
                    {HistType::kTH1F, {dcaZAxis}}});
      registry.add({"Purity/xReco",
                    ";#Delta z_{vtx}^{reco} (cm);tracks",
                    {HistType::kTH1F, {zAxis}}});
      registry.add({"Purity/xTrue",
                    ";#Delta z_{vtx}^{true} (cm);tracks",
                    {HistType::kTH1F, {zAxis}}});
      registry.add({"Purity/yReco",
                    ";#Delta z_{vtx}^{reco} (cm);tracks",
                    {HistType::kTH1F, {zAxis}}});
      registry.add({"Purity/yTrue",
                    ";#Delta z_{vtx}^{true} (cm);tracks",
                    {HistType::kTH1F, {zAxis}}});
      registry.add({"Purity/zReco",
                    ";#Delta z_{vtx}^{reco} (cm);tracks",
                    {HistType::kTH1F, {zAxis}}});
      registry.add({"Purity/zTrue",
                    ";#Delta z_{vtx}^{true} (cm);tracks",
                    {HistType::kTH1F, {zAxis}}});

      // Vertex positions: store both reco and truth in one THnSparse with a status axis
      registry.add({"Purity/VtxXYZTruth",
                    "; V_{x} (cm); V_{y} (cm); V_{z} (cm)",
                    {HistType::kTHnSparseF, {vxAxis, vyAxis, zAxis}}});
      registry.add({"Purity/VtxXYZReco",
                    "; V_{x} (cm); V_{y} (cm); V_{z} (cm)",
                    {HistType::kTHnSparseF, {vxAxis, vyAxis, zAxis}}});

      // Delta vertex positions (reco - true)
      registry.add({"Purity/DeltaVtxXYZ",
                    "; #DeltaV_{x} (cm); #DeltaV_{y} (cm); #DeltaV_{z} (cm)",
                    {HistType::kTHnSparseF, {deltaVxAxis, deltaVyAxis, deltazAxis}}});

      registry.add({"Purity/DeltaXRight",
                    ";#Delta x_{vtx}^{reco-true} (cm) (right);tracks",
                    {HistType::kTH1F, {dcaZAxis}}});
      registry.add({"Purity/DeltaYRight",
                    ";#Delta y_{vtx}^{reco-true} (cm) (right);tracks",
                    {HistType::kTH1F, {dcaZAxis}}});
      registry.add({"Purity/DeltaZRight",
                    ";#Delta z_{vtx}^{reco-true} (cm) (right);tracks",
                    {HistType::kTH1F, {dcaZAxis}}});

      // 5D THnSparse histograms: (eta, DCAxy, DCAz, DCAx, DCAy)

      // 1) All / Primary / Secondary after track selection (no vertex classification)
      registry.add({"Purity/RecoSparseAll",
                    ";#eta;DCA_{xy} (cm);DCA_{z} (cm);DCA_{x} (cm);DCA_{y} (cm)",
                    {HistType::kTHnSparseF,
                     {etaBinning, dcaXyAxis, dcaZAxis, dcaXAxis, dcaYAxis}}});
      registry.add({"RecoSparseAllBest",
                    ";#eta;DCA_{xy} (cm);DCA_{z} (cm);DCA_{x} (cm);DCA_{y} (cm)",
                    {HistType::kTHnSparseF,
                     {etaBinning, dcaXyAxis, dcaZAxis, dcaXAxis, dcaYAxis}}});
      registry.add({"RecoSparseAllBestWrong",
                    ";#eta;DCA_{xy} (cm);DCA_{z} (cm);DCA_{x} (cm);DCA_{y} (cm)",
                    {HistType::kTHnSparseF,
                     {etaBinning, dcaXyAxis, dcaZAxis, dcaXAxis, dcaYAxis}}});

      registry.add({"Purity/RecoSparsePrimary",
                    ";#eta;DCA_{xy} (cm);DCA_{z} (cm);DCA_{x} (cm);DCA_{y} (cm)",
                    {HistType::kTHnSparseF,
                     {etaBinning, dcaXyAxis, dcaZAxis, dcaXAxis, dcaYAxis}}});

      registry.add({"Purity/RecoSparseSecondary",
                    ";#eta;DCA_{xy} (cm);DCA_{z} (cm);DCA_{x} (cm);DCA_{y} (cm)",
                    {HistType::kTHnSparseF,
                     {etaBinning, dcaXyAxis, dcaZAxis, dcaXAxis, dcaYAxis}}});

      // 2) Right-vertex: all / primary / secondary
      registry.add({"Purity/RecoSparseRightAll",
                    ";#eta;DCA_{xy} (cm);DCA_{z} (cm);DCA_{x} (cm);DCA_{y} (cm)",
                    {HistType::kTHnSparseF,
                     {etaBinning, dcaXyAxis, dcaZAxis, dcaXAxis, dcaYAxis}}});

      registry.add({"Purity/RecoSparseRightPrimary",
                    ";#eta;DCA_{xy} (cm);DCA_{z} (cm);DCA_{x} (cm);DCA_{y} (cm)",
                    {HistType::kTHnSparseF,
                     {etaBinning, dcaXyAxis, dcaZAxis, dcaXAxis, dcaYAxis}}});

      registry.add({"Purity/RecoSparseRightSecondary",
                    ";#eta;DCA_{xy} (cm);DCA_{z} (cm);DCA_{x} (cm);DCA_{y} (cm)",
                    {HistType::kTHnSparseF,
                     {etaBinning, dcaXyAxis, dcaZAxis, dcaXAxis, dcaYAxis}}});

      // 3) Wrong-vertex: all / primary / secondary
      registry.add({"Purity/RecoSparseWrongAll",
                    ";#eta;DCA_{xy} (cm);DCA_{z} (cm);DCA_{x} (cm);DCA_{y} (cm)",
                    {HistType::kTHnSparseF,
                     {etaBinning, dcaXyAxis, dcaZAxis, dcaXAxis, dcaYAxis}}});

      registry.add({"Purity/RecoSparseWrongPrimary",
                    ";#eta;DCA_{xy} (cm);DCA_{z} (cm);DCA_{x} (cm);DCA_{y} (cm)",
                    {HistType::kTHnSparseF,
                     {etaBinning, dcaXyAxis, dcaZAxis, dcaXAxis, dcaYAxis}}});

      registry.add({"Purity/RecoSparseWrongSecondary",
                    ";#eta;DCA_{xy} (cm);DCA_{z} (cm);DCA_{x} (cm);DCA_{y} (cm)",
                    {HistType::kTHnSparseF,
                     {etaBinning, dcaXyAxis, dcaZAxis, dcaXAxis, dcaYAxis}}});

      // 4) Generator-truth DCA 5D THnSparse: (eta_truth, DCAxy_truth, DCAz_truth, DCAx_truth, DCAy_truth)
      //    booked under Tracks/dca/Truth/
      registry.add({"Tracks/dca/Truth/THnDCAxyBestGenTruthAll",
                    ";#eta_{truth};DCA_{xy}^{truth} (cm);DCA_{z}^{truth} (cm);DCA_{x}^{truth} (cm);DCA_{y}^{truth} (cm)",
                    {HistType::kTHnSparseF,
                     {etaBinning, dcaXyAxis, dcaZAxis, dcaXAxis, dcaYAxis}}});
      registry.add({"Tracks/dca/Truth/THnDCAxyBestGenTruthPrim",
                    ";#eta_{truth};DCA_{xy}^{truth} (cm);DCA_{z}^{truth} (cm);DCA_{x}^{truth} (cm);DCA_{y}^{truth} (cm) (primary)",
                    {HistType::kTHnSparseF,
                     {etaBinning, dcaXyAxis, dcaZAxis, dcaXAxis, dcaYAxis}}});
      registry.add({"Tracks/dca/Truth/THnDCAxyBestGenTruthSec",
                    ";#eta_{truth};DCA_{xy}^{truth} (cm);DCA_{z}^{truth} (cm);DCA_{x}^{truth} (cm);DCA_{y}^{truth} (cm) (secondary)",
                    {HistType::kTHnSparseF,
                     {etaBinning, dcaXyAxis, dcaZAxis, dcaXAxis, dcaYAxis}}});

      // Right-vertex truth DCA
      registry.add({"Tracks/dca/Truth/THnDCAxyBestGenTruthRightAll",
                    ";#eta_{truth};DCA_{xy}^{truth} (cm);DCA_{z}^{truth} (cm);DCA_{x}^{truth} (cm);DCA_{y}^{truth} (cm) (right vertex)",
                    {HistType::kTHnSparseF,
                     {etaBinning, dcaXyAxis, dcaZAxis, dcaXAxis, dcaYAxis}}});
      registry.add({"Tracks/dca/Truth/THnDCAxyBestGenTruthRightPrim",
                    ";#eta_{truth};DCA_{xy}^{truth} (cm);DCA_{z}^{truth} (cm);DCA_{x}^{truth} (cm);DCA_{y}^{truth} (cm) (primary, right vertex)",
                    {HistType::kTHnSparseF,
                     {etaBinning, dcaXyAxis, dcaZAxis, dcaXAxis, dcaYAxis}}});
      registry.add({"Tracks/dca/Truth/THnDCAxyBestGenTruthRightSec",
                    ";#eta_{truth};DCA_{xy}^{truth} (cm);DCA_{z}^{truth} (cm);DCA_{x}^{truth} (cm);DCA_{y}^{truth} (cm) (secondary, right vertex)",
                    {HistType::kTHnSparseF,
                     {etaBinning, dcaXyAxis, dcaZAxis, dcaXAxis, dcaYAxis}}});

      // Wrong-vertex truth DCA
      registry.add({"Tracks/dca/Truth/THnDCAxyBestGenTruthWrongAll",
                    ";#eta_{truth};DCA_{xy}^{truth} (cm);DCA_{z}^{truth} (cm);DCA_{x}^{truth} (cm);DCA_{y}^{truth} (cm) (wrong vertex)",
                    {HistType::kTHnSparseF,
                     {etaBinning, dcaXyAxis, dcaZAxis, dcaXAxis, dcaYAxis}}});
      registry.add({"Tracks/dca/Truth/THnDCAxyBestGenTruthWrongPrim",
                    ";#eta_{truth};DCA_{xy}^{truth} (cm);DCA_{z}^{truth} (cm);DCA_{x}^{truth} (cm);DCA_{y}^{truth} (cm) (primary, wrong vertex)",
                    {HistType::kTHnSparseF,
                     {etaBinning, dcaXyAxis, dcaZAxis, dcaXAxis, dcaYAxis}}});
      registry.add({"Tracks/dca/Truth/THnDCAxyBestGenTruthWrongSec",
                    ";#eta_{truth};DCA_{xy}^{truth} (cm);DCA_{z}^{truth} (cm);DCA_{x}^{truth} (cm);DCA_{y}^{truth} (cm) (secondary, wrong vertex)",
                    {HistType::kTHnSparseF,
                     {etaBinning, dcaXyAxis, dcaZAxis, dcaXAxis, dcaYAxis}}});

      // 5) Delta-DCA THnSparse: (ΔDCAxy, ΔDCAz, ΔDCAx, ΔDCAy) for right/wrong vertex, all/primary/secondary
      registry.add({"Tracks/dca/Truth/THnDeltaDCARightAll",
                    ";#Delta DCA_{xy} (cm);#Delta DCA_{z} (cm);#Delta DCA_{x} (cm);#Delta DCA_{y} (cm)",
                    {HistType::kTHnSparseF,
                     {dcaXyAxis, dcaZAxis, dcaXAxis, dcaYAxis}}});
      registry.add({"Tracks/dca/Truth/THnDeltaDCARightPrim",
                    ";#Delta DCA_{xy} (cm);#Delta DCA_{z} (cm);#Delta DCA_{x} (cm);#Delta DCA_{y} (cm) (primary, right vertex)",
                    {HistType::kTHnSparseF,
                     {dcaXyAxis, dcaZAxis, dcaXAxis, dcaYAxis}}});
      registry.add({"Tracks/dca/Truth/THnDeltaDCARightSec",
                    ";#Delta DCA_{xy} (cm);#Delta DCA_{z} (cm);#Delta DCA_{x} (cm);#Delta DCA_{y} (cm) (secondary, right vertex)",
                    {HistType::kTHnSparseF,
                     {dcaXyAxis, dcaZAxis, dcaXAxis, dcaYAxis}}});

      registry.add({"Tracks/dca/Truth/THnDeltaDCAWrongAll",
                    ";#Delta DCA_{xy} (cm);#Delta DCA_{z} (cm);#Delta DCA_{x} (cm);#Delta DCA_{y} (cm) (wrong vertex)",
                    {HistType::kTHnSparseF,
                     {dcaXyAxis, dcaZAxis, dcaXAxis, dcaYAxis}}});
      registry.add({"Tracks/dca/Truth/THnDeltaDCAWrongPrim",
                    ";#Delta DCA_{xy} (cm);#Delta DCA_{z} (cm);#Delta DCA_{x} (cm);#Delta DCA_{y} (cm) (primary, wrong vertex)",
                    {HistType::kTHnSparseF,
                     {dcaXyAxis, dcaZAxis, dcaXAxis, dcaYAxis}}});
      registry.add({"Tracks/dca/Truth/THnDeltaDCAWrongSec",
                    ";#Delta DCA_{xy} (cm);#Delta DCA_{z} (cm);#Delta DCA_{x} (cm);#Delta DCA_{y} (cm) (secondary, wrong vertex)",
                    {HistType::kTHnSparseF,
                     {dcaXyAxis, dcaZAxis, dcaXAxis, dcaYAxis}}});

      registry.add({"Purity/RecoSparseAll_EventCount",
                    ";events;counts",
                    {HistType::kTH1F, {{1, 0.5, 1.5}}}});
      registry.add({"Purity/RecoSparseAll_EventCountBest",
                    ";events;counts",
                    {HistType::kTH1F, {{1, 0.5, 1.5}}}});
      registry.add({"Purity/RecoSparseRightAll_EventCount",
                    ";events;counts",
                    {HistType::kTH1F, {{1, 0.5, 1.5}}}});
      registry.add({"Purity/RecoSparseWrongAll_EventCount",
                    ";events;counts",
                    {HistType::kTH1F, {{1, 0.5, 1.5}}}});
      registry.add({"Purity/BestRecoColNotFound",
                    ";events;counts",
                    {HistType::kTH1F, {{1, 0.5, 1.5}}}});
    }

    if (doprocessGen) {
      registry.add({"EventsNtrkZvtxGen",
                    "; N_{trk}; #it{z}_{vtx} (cm); events",
                    {HistType::kTH2F, {multAxis, zAxis}}});
      registry.add({"EventsNtrkZvtxGen_t",
                    "; N_{trk}; #it{z}_{vtx} (cm); events",
                    {HistType::kTH2F, {multAxis, zAxis}}});
      registry.add({"EventsNtrkZvtxGen_gt0",
                    "; N_{trk}; #it{z}_{vtx} (cm); events",
                    {HistType::kTH2F, {multAxis, zAxis}}});
      registry.add({"EventsNtrkZvtxGen_gt0t",
                    "; N_{trk}; #it{z}_{vtx} (cm); events",
                    {HistType::kTH2F, {multAxis, zAxis}}});
      registry.add({"TracksEtaZvtxGen",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {etaBinning, zAxis}}});
      registry.add({"TracksEtaZvtxGen_t",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {etaBinning, zAxis}}});
      registry.add({"TracksEtaZvtxGen_gt0",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {etaBinning, zAxis}}});
      registry.add({"TracksEtaZvtxGen_gt0t",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {etaBinning, zAxis}}});
      registry.add({"TracksPhiEtaGen",
                    "; #varphi; #eta; tracks",
                    {HistType::kTH2F, {phiAxis, etaBinning}}});
      registry.add({"TracksPhiEtaGen_gt0",
                    "; #varphi; #eta; tracks",
                    {HistType::kTH2F, {phiAxis, etaBinning}}});
      registry.add({"TracksPhiEtaGen_gt0t",
                    "; #varphi; #eta; tracks",
                    {HistType::kTH2F, {phiAxis, etaBinning}}});
      registry.add({"TracksPhiZvtxGen",
                    "; #varphi; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {phiAxis, zAxis}}}); //
      registry.add({"TracksToPartPtEta",
                    " ; p_{T} (GeV/c); #eta",
                    {HistType::kTH2F, {ptAxis, etaBinning}}}); //
      registry.add({"TracksPtEtaGen",
                    " ; p_{T} (GeV/c); #eta",
                    {HistType::kTH2F, {ptAxis, etaBinning}}});
      registry.add({"TracksPtEtaGen_t",
                    " ; p_{T} (GeV/c); #eta",
                    {HistType::kTH2F, {ptAxis, etaBinning}}});
      registry.add({"NotFoundEventZvtx",
                    " ; #it{z}_{vtx} (cm)",
                    {HistType::kTH1F, {zAxis}}});
      registry.add({"EventsZposDiff",
                    " ; Z_{rec} - Z_{gen} (cm)",
                    {HistType::kTH1F, {deltazAxis}}});

      registry.add({"TracksEtaZvtxGen_gt0_primary",
                    "; #eta; #it{z}_{vtx} (cm); primary tracks",
                    {HistType::kTH2F, {etaBinning, zAxis}}});
      registry.add({"TracksPhiEtaGen_gt0_primary",
                    "; #varphi; #eta; primary tracks",
                    {HistType::kTH2F, {phiAxis, etaBinning}}});
      registry.add({"TracksPtZvtxGen_gt0_primary",
                    "; p_{T} (GeV/c); #it{z}_{vtx} (cm); primary tracks",
                    {HistType::kTH2F, {ptAxis, zAxis}}});
      registry.add({"EventsSplitMult", " ; N_{gen}", {HistType::kTH1F, {multAxis}}});
      auto heff = registry.get<TH1>(HIST("EventEfficiency"));
      x = heff->GetXaxis();
      x->SetBinLabel(static_cast<int>(EventEfficiencyBin::Generated), "Generated");
      x->SetBinLabel(static_cast<int>(EventEfficiencyBin::GeneratedInelGt0), "Generated INEL>0");
      x->SetBinLabel(static_cast<int>(EventEfficiencyBin::Reconstructed), "Reconstructed");
      x->SetBinLabel(static_cast<int>(EventEfficiencyBin::Selected), "Selected");
      x->SetBinLabel(static_cast<int>(EventEfficiencyBin::SelectedInelGt0), "Selected INEL>0");
    }

    if (doprocessMultReassoc || doprocessMultReassoc3d) {
      // Purity denominator histograms (reco side after DCAxy selection)

      registry.add({"Tracks/Control/DeltaZ",
                    " ; #it{z_{orig}}-#it{z_{reass}}",
                    {HistType::kTH1F, {zAxis}}});

      registry.add({"Tracks/Control/TrackAmbDegree",
                    " ; N_{coll}^{comp}",
                    {HistType::kTH1F, {{51, -0.5, 50.5}}}});
      registry.add({"Tracks/Control/TrackIsAmb",
                    " ; isAmbiguous",
                    {HistType::kTH1I, {{2, -0.5, 1.5}}}});

      registry.add({"Tracks/Control/ReassignedTracksEtaZvtx",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {etaBinning, zAxis}}});
      registry.add({"Tracks/Control/ReassignedTracksPhiEta",
                    "; #varphi; #eta; tracks",
                    {HistType::kTH2F, {phiAxis, etaBinning}}});
      registry.add({"Tracks/Control/ReassignedVertexCorr",
                    "; #it{z}_{vtx}^{orig} (cm); #it{z}_{vtx}^{re} (cm)",
                    {HistType::kTH2F, {zAxis, zAxis}}});

      registry.add({"Tracks/Control/notReassignedTracksEtaZvtx",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {etaBinning, zAxis}}});
      registry.add({"Tracks/Control/notReassignedTracksPhiEta",
                    "; #varphi; #eta; tracks",
                    {HistType::kTH2F, {phiAxis, etaBinning}}});
      registry.add({"Tracks/Control/notReassignedVertexCorr",
                    "; #it{z}_{vtx}^{orig} (cm); #it{z}_{vtx}^{re} (cm)",
                    {HistType::kTH2F, {zAxis, zAxis}}});
      registry.add({"Tracks/Control/Chi2NDF",
                    " ; #chi^{2}/ndf",
                    {HistType::kTH1F, {{5000, 0.0, 5000.0}}}});
      registry.add({"Tracks/Control/amb/AmbTracksEtaZvtx",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {etaBinning, zAxis}}}); //

      registry.add({"Tracks/Control/woOrp/nTrk",
                    " ; N_{Trk}^{all}",
                    {HistType::kTH1F, {{701, -0.5, 700.5}}}}); //
      registry.add({"Tracks/Control/amb/nTrkAmb",
                    " ; N_{Trk}^{amb}",
                    {HistType::kTH1F, {{701, -0.5, 700.5}}}}); //
      registry.add({"Tracks/Control/nonamb/nTrkNonAmb",
                    " ; N_{Trk}^{nonamb}",
                    {HistType::kTH1F, {{701, -0.5, 700.5}}}}); //

      registry.add({"Tracks/Control/amb/AmbTracksPhiEta",
                    "; #varphi; #eta; tracks",
                    {HistType::kTH2F, {phiAxis, etaBinning}}}); //
      registry.add({"Tracks/Control/amb/AmbVertexCorr",
                    "; #it{z}_{vtx}^{orig} (cm); #it{z}_{vtx}^{re} (cm)",
                    {HistType::kTH2F, {zAxis, zAxis}}}); //
      registry.add({"Tracks/Control/amb/EtaZvtxAmb_gt0",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {etaBinning, zAxis}}}); //
      registry.add({"Tracks/Control/amb/DCAxy_amb", " ; DCA_{xy} (cm) ambiguous",
                    //  {HistType::kTH1F,{{100000, 0.5, 100000.0}}}}); //
                    {HistType::kTH1F, {dcaXyAxis}}}); //

      registry.add({"Tracks/Control/nonamb/nonAmbTracksEtaZvtx",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {etaBinning, zAxis}}}); //

      registry.add({"Tracks/Control/nonamb/nonAmbTracksPhiEta",
                    "; #varphi; #eta; tracks",
                    {HistType::kTH2F, {phiAxis, etaBinning}}}); //
      registry.add({"Tracks/Control/nonamb/nonAmbVertexCorr",
                    "; #it{z}_{vtx}^{orig} (cm); #it{z}_{vtx}^{re} (cm)",
                    {HistType::kTH2F, {zAxis, zAxis}}}); //
      registry.add({"Tracks/Control/nonamb/EtaZvtxNonAmb_gt0",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {etaBinning, zAxis}}}); //
      registry.add({"Tracks/Control/nonamb/DCAxy_nonamb", " ; DCA_{xy}(cm) non-ambiguous",
                    //  {HistType::kTH1F,{{100000, 0.5, 100000.0}}}}); //
                    {HistType::kTH1F, {{dcaXyAxis}}}}); //

      registry.add({"Tracks/Control/woOrp/woOrpTracksEtaZvtx",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {etaBinning, zAxis}}}); //
      registry.add({"Tracks/Control/woOrp/woOrpEtaZvtx_gt0",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {etaBinning, zAxis}}}); //
      registry.add({"Tracks/2Danalysis/EtaZvtx",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {etaBinning, zAxis}}}); //
      registry.add({"Tracks/2Danalysis/EtaZvtx_sel8",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {etaBinning, zAxis}}}); //
      registry.add({"Tracks/2Danalysis/EtaZvtx_sel8_inelgt0",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {etaBinning, zAxis}}}); //
      registry.add({"Tracks/2Danalysis/EtaZvtx_sel8_inelfwdgt0",
                    "; #eta; #it{z}_{vtx} (cm); tracks",
                    {HistType::kTH2F, {etaBinning, zAxis}}}); //
      registry.add({"Tracks/Control/woOrp/woOrpTracksPhiEta",
                    "; #varphi; #eta; tracks",
                    {HistType::kTH2F, {phiAxis, etaBinning}}}); //
      registry.add({"Tracks/Control/woOrp/woOrpVertexCorr",
                    "; #it{z}_{vtx}^{orig} (cm); #it{z}_{vtx}^{re} (cm)",
                    {HistType::kTH2F, {zAxis, zAxis}}}); //
      registry.add({"Tracks/Control/woOrp/DCAxy_woOrp", " ; DCA_{xy}(cm) w/o orphan",
                    //  {HistType::kTH1F,{{100000, 0.5, 100000.0}}}}); //
                    {HistType::kTH1F, {{dcaXyAxis}}}}); //

      if (doprocessMultReassoc3d) {
        // DCAz histograms analogous to DCAxy, only for 3D reassociation
        registry.add({"Tracks/Control/DCAZ",
                      " ; DCA_{Z} (cm)",
                      {HistType::kTH1F, {dcaZAxis}}});
        registry.add({"Tracks/Control/amb/DCAz_amb",
                      " ; DCA_{z} (cm) ambiguous",
                      {HistType::kTH1F, {dcaZAxis}}});
        registry.add({"Tracks/Control/nonamb/DCAz_nonamb",
                      " ; DCA_{z}(cm) non-ambiguous",
                      {HistType::kTH1F, {dcaZAxis}}});
        registry.add({"Tracks/Control/woOrp/DCAz_woOrp",
                      " ; DCA_{z}(cm) w/o orphan",
                      {HistType::kTH1F, {dcaZAxis}}});
      }

      registry.add({"collisionID", " ; Collision ID",
                    //  {HistType::kTH1F,{{100000, 0.5, 100000.0}}}}); //
                    {HistType::kTH1F, {{100000, -50000.0, 50000.0}}}}); //
      registry.add({"collisionIDamb", " ; Collision ID amb",
                    //  {HistType::kTH1F,{{100000, 0.5, 100000.0}}}}); //
                    {HistType::kTH1F, {{100000, -50000.0, 50000.0}}}});                                                                                        //
      registry.add({"NonambEventCounts", " ; EventCounts Nonamb", {HistType::kTH1F, {{1, 0.5, 1.5}}}});                                                        //
      registry.add({"hNumCollisionsNonAmb_InelMFT", " ; Number of Collisions with Non-Ambiguous Tracks;Count;Frequency", {HistType::kTH1F, {{1, 0.5, 1.5}}}}); //
      registry.add({"hNumCollisionsAmb_InelMFT", " ; Number of Collisions with Non-Ambiguous Tracks;Count;Frequency", {HistType::kTH1F, {{1, 0.5, 1.5}}}});    //
      registry.add({"hNumCollisions_InelMFT", " ; Number of selected events with Inel>0 and MFT>0;Count;Frequency", {HistType::kTH1F, {{1, 0.5, 1.5}}}});      //
      registry.add({"hNumCollisions_Inel", " ; Number of selected events with Inel>0;Count;Frequency", {HistType::kTH1F, {{1, 0.5, 1.5}}}});                   //
      registry.add({"ambEventCounts", " ; EventCounts Nonamb", {HistType::kTH1F, {{1, 0.5, 1.5}}}});                                                           //
    }

    if (doprocessCountingCentrality) {
      registry.add({"Events/Centrality/Selection",
                    ";status;centrality;events",
                    {HistType::kTH2F, {{3, 0.5, 3.5}, centAxis}}});
      auto hstat = registry.get<TH2>(HIST("Events/Centrality/Selection"));
      auto* x = hstat->GetXaxis();
      x->SetBinLabel(static_cast<int>(CentralitySelectionBin::All), "All");
      x->SetBinLabel(static_cast<int>(CentralitySelectionBin::Selected), "Selected");
      x->SetBinLabel(static_cast<int>(CentralitySelectionBin::Rejected), "Rejected");

      registry.add({"Events/Centrality/NtrkZvtx",
                    "; N_{trk}; Z_{vtx} (cm); centrality",
                    {HistType::kTH3F, {multAxis, zAxis, centAxis}}});
      registry.add({"Tracks/Centrality/EtaZvtx",
                    "; #eta; Z_{vtx} (cm); centrality",
                    {HistType::kTH3F, {etaBinning, zAxis, centAxis}}});
      registry.add({"Tracks/Centrality/PhiEta",
                    "; #varphi; #eta; centrality",
                    {HistType::kTH3F, {phiAxis, etaBinning, centAxis}}});
      registry.add({"Tracks/Centrality/Control/PtEta",
                    " ; p_{T} (GeV/c); #eta; centrality",
                    {HistType::kTH3F, {ptAxis, etaBinning, centAxis}}});
      registry.add({"Tracks/Centrality/Control/DCAXYPt",
                    " ; p_{T} (GeV/c) ; DCA_{XY} (cm); centrality",
                    {HistType::kTH3F, {ptAxis, dcaXyAxis, centAxis}}});
      registry.add({"Tracks/Centrality/Control/ReassignedDCAXYPt",
                    " ; p_{T} (GeV/c) ; DCA_{XY} (cm); centrality",
                    {HistType::kTH3F, {ptAxis, dcaXyAxis, centAxis}}});
      registry.add({"Tracks/Centrality/Control/ExtraDCAXYPt",
                    " ; p_{T} (GeV/c) ; DCA_{XY} (cm); centrality",
                    {HistType::kTH3F, {ptAxis, dcaXyAxis, centAxis}}});
      registry.add({"Tracks/Centrality/Control/ExtraTracksEtaZvtx",
                    "; #eta; Z_{vtx} (cm); centrality",
                    {HistType::kTH3F, {etaBinning, zAxis, centAxis}}});
      registry.add({"Tracks/Centrality/Control/ExtraTracksPhiEta",
                    "; #varphi; #eta; centrality",
                    {HistType::kTH3F, {phiAxis, etaBinning, centAxis}}});
      registry.add({"Tracks/Centrality/Control/ReassignedTracksEtaZvtx",
                    "; #eta; Z_{vtx} (cm); centrality",
                    {HistType::kTH3F, {etaBinning, zAxis, centAxis}}});
      registry.add({"Tracks/Centrality/Control/ReassignedTracksPhiEta",
                    "; #varphi; #eta; centrality",
                    {HistType::kTH3F, {phiAxis, etaBinning, centAxis}}});
      registry.add({"Tracks/Centrality/Control/ReassignedVertexCorr",
                    "; Z_{vtx}^{orig} (cm); Z_{vtx}^{re} (cm); centrality",
                    {HistType::kTH3F, {zAxis, zAxis, centAxis}}});
    }

    if (doprocessGenCent) {
      registry.add({"Events/Centrality/EventEfficiency",
                    ";status;centrality;events",
                    {HistType::kTH2F, {{2, 0.5, 2.5}, centAxis}}});
      auto heff = registry.get<TH2>(HIST("Events/Centrality/EventEfficiency"));
      auto* x = heff->GetXaxis();
      x->SetBinLabel(1, "Generated");
      x->SetBinLabel(2, "Selected");

      registry.add("Events/Centrality/CentPercentileMCGen",
                   "CentPercentileMCGen", kTH1D, {centAxis}, false);
      registry.add({"Events/Centrality/NtrkZvtxGen",
                    "; N_{trk}; Z_{vtx} (cm); centrality",
                    {HistType::kTH3F, {multAxis, zAxis, centAxis}}});
      registry.add({"Events/Centrality/NtrkZvtxGen_t",
                    "; N_{trk}; Z_{vtx} (cm); centrality",
                    {HistType::kTH3F, {multAxis, zAxis, centAxis}}});
      registry.add({"Tracks/Centrality/EtaZvtxGen_t",
                    "; #eta; Z_{vtx} (cm); centrality",
                    {HistType::kTH3F, {etaBinning, zAxis, centAxis}}});
      registry.add({"Tracks/Centrality/EtaZvtxGen",
                    "; #eta; Z_{vtx} (cm); centrality",
                    {HistType::kTH3F, {etaBinning, zAxis, centAxis}}});
      registry.add({"Tracks/Centrality/PhiEtaGen",
                    "; #varphi; #eta; centrality",
                    {HistType::kTH3F, {phiAxis, etaBinning, centAxis}}});
    }
  }

  using FullBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;
  void processTagging(FullBCs const& bcs,
                      soa::Join<aod::Collisions, aod::EvSels> const& collisions)
  {

    std::vector<typename std::decay_t<decltype(collisions)>::iterator> cols;
    for (const auto& bc : bcs) {
      if (!useEvSel ||
          (useEvSel && ((bc.selection_bit(aod::evsel::kIsBBT0A) &&
                         bc.selection_bit(aod::evsel::kIsBBT0C)) != 0))) {
        registry.fill(HIST("EventSelection"), static_cast<int>(EventSelectionBin::GoodBCs));
        cols.clear();
        for (const auto& collision : collisions) {
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
          registry.fill(HIST("EventSelection"), static_cast<int>(EventSelectionBin::BCsWithCollisions));
          if (cols.size() > 1) {
            registry.fill(HIST("EventSelection"), static_cast<int>(EventSelectionBin::BCsWithPileupSplitting));
          }
        }
      }
    }
  }

  bool passGenRecoTrackMode(auto const& track) const
  {
    const auto compatibleColls = track.compatibleCollIds();
    const auto nCompatibleColls = compatibleColls.size();

    if (nCompatibleColls == NoCompatibleCollisions) {
      return false;
    }

    const auto mode =
      static_cast<GenRecoTimeComTrackMode>(cfgGenRecoTimeComTrackMode.value);

    switch (mode) {
      case GenRecoTimeComTrackMode::AllNonOrphan:
        return nCompatibleColls > NoCompatibleCollisions;
      case GenRecoTimeComTrackMode::NonOrphanNonAmbiguous:
        return nCompatibleColls == SingleCompatibleCollision;
      case GenRecoTimeComTrackMode::NonOrphanAmbiguous:
        return nCompatibleColls > SingleCompatibleCollision;
      default:
        return false;
    }
  }

  PROCESS_SWITCH(PseudorapidityDensityMFT, processTagging,
                 "Collect event sample stats", true);

  static constexpr float SampleEtaMax = -2.8f;
  static constexpr float SampleEtaMin = -3.2f;

  Partition<aod::MFTTracks> sample =
    (aod::fwdtrack::eta < SampleEtaMax) && (aod::fwdtrack::eta > SampleEtaMin);

  static constexpr float CentralEtaMax = 1.f;
  Partition<aod::Tracks> sampleCentral = (nabs(aod::track::eta) < CentralEtaMax);

  static constexpr int ValidBestCollisionIdMin = 0;
  static constexpr int InvalidCollisionId = -1;
  static constexpr float MaxBestDcaXy = 2.f;
  expressions::Filter atrackFilter =
    (aod::fwdtrack::bestCollisionId >= ValidBestCollisionIdMin) &&
    (aod::fwdtrack::eta < ForwardEtaMax) &&
    (aod::fwdtrack::eta > ForwardEtaMin) &&
    (nabs(aod::fwdtrack::bestDCAXY) <= MaxBestDcaXy);

  using CollwEv = soa::Join<aod::Collisions, aod::EvSels>;

  expressions::Filter trackSelectionCentral =
    ((aod::track::trackCutFlag & TrackSelectionIts) == TrackSelectionIts) &&
    ifnode((aod::track::v001::detectorMap & (uint8_t)o2::aod::track::TPC) ==
             (uint8_t)o2::aod::track::TPC,
           (aod::track::trackCutFlag & TrackSelectionTpc) ==
             TrackSelectionTpc,
           true) &&
    ((aod::track::trackCutFlag & TrackSelectionDca) == TrackSelectionDca) &&
    (nabs(aod::track::eta) < estimatorEta);

  using FiCentralTracks = soa::Filtered<
    soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection,
              aod::TracksDCA>>; // central tracks for INEL>0

  void processMult(CollwEv::iterator const& collision,
                   aod::MFTTracks const& tracks,
                   FiCentralTracks const& midtracks, aod::Tracks const&)
  {

    registry.fill(HIST("EventSelection"), static_cast<int>(EventSelectionBin::All));
    if (!useEvSel || (useEvSel && collision.sel8())) {
      registry.fill(HIST("EventSelection"), static_cast<int>(EventSelectionBin::VzSelected));
      auto z = collision.posZ();
      auto perCollisionSample = sampleCentral->sliceByCached(
        o2::aod::track::collisionId, collision.globalIndex(), cache);
      auto nTrk = perCollisionSample.size();

      registry.fill(HIST("EventsNtrkZvtx"), nTrk, z);

      if (midtracks.size() > 0) // INEL>0
      {
        registry.fill(HIST("EventSelection"), static_cast<int>(EventSelectionBin::Sel8VzInelGt0));
        registry.fill(HIST("EventsNtrkZvtx_gt0"), nTrk, z);
      }

      if (tracks.size() > 0) {
        for (const auto& track : tracks) {

          float phi = track.phi();
          o2::math_utils::bringTo02Pi(phi);

          if (usePhiCut) {
            if ((phi < cfgPhiCut) ||
                ((phi > o2::constants::math::PI - cfgPhiCut) && (phi < o2::constants::math::PI + cfgPhiCut)) ||
                (phi > o2::constants::math::TwoPI - cfgPhiCut) ||
                ((phi > ((o2::constants::math::PIHalf - 0.1) * o2::constants::math::PI) - cfgPhiCut) &&
                 (phi < ((o2::constants::math::PIHalf - 0.1) * o2::constants::math::PI) + cfgPhiCut)))
              continue;
          }

          registry.fill(HIST("TracksEtaZvtx"), track.eta(), z);
          if (midtracks.size() > 0) // INEL>0
          {
            registry.fill(HIST("Tracks/EtaZvtx_gt0"), track.eta(), z);
          }
          registry.fill(HIST("TracksPhiEta"), phi, track.eta());
          registry.fill(HIST("TracksPtEta"), track.pt(), track.eta());
          if ((track.eta() < ForwardEtaMax) && (track.eta() > ForwardEtaMin)) {
            registry.fill(HIST("TracksPhiZvtx"), phi, z);
          }
        }
      }

    } else {
      registry.fill(HIST("EventSelection"), static_cast<int>(EventSelectionBin::Rejected));
    }
  }

  PROCESS_SWITCH(PseudorapidityDensityMFT, processMult,
                 "Process reco or data info", true);
  // Common implementation for both BestCollisionsFwd and BestCollisionsFwd3d
  template <typename RetracksT>
  void processMultReassocCommon(CollwEv::iterator const& collision,
                                o2::aod::MFTTracks const&,
                                RetracksT const& retracks,
                                FiCentralTracks const& midtracks, aod::Tracks const&)
  {
    registry.fill(HIST("EventSelection"), static_cast<int>(EventSelectionBin::All));
    auto perCollisionSample = sampleCentral->sliceByCached(
      o2::aod::track::collisionId, collision.globalIndex(), cache);
    auto nTrk = perCollisionSample.size();
    auto z = collision.posZ();
    registry.fill(HIST("EventsNtrkZvtx"), nTrk, z);
    if ((z >= cfgVzCut1) && (z <= cfgVzCut2)) {
      registry.fill(HIST("Tracks/2Danalysis/EventsNtrkZvtx_all"), nTrk, z);
      registry.fill(HIST("EventSelection"), static_cast<int>(EventSelectionBin::Vz));
      for (const auto& retrack : retracks) {
        auto track = retrack.mfttrack();
        float ndf = getTrackNdf(track);
        float chi2ndf = track.chi2() / ndf;
        float phi = track.phi();
        o2::math_utils::bringTo02Pi(phi);
        if (usePhiCut) {
          if ((phi <= PhiVetoLow) ||
              ((phi >= PhiVetoPiMin) && (phi <= PhiVetoPiMax)) ||
              (phi >= PhiVetoHigh))
            continue;
        }
        float dcaXyCut = retrack.bestDCAXY();
        if (useDCAxyCut) {
          if (dcaXyCut > maxDCAxy)
            continue;
        }
        if constexpr (std::is_same_v<RetracksT, soa::SmallGroups<aod::BestCollisionsFwd3d>>) {
          float dcaZCut = retrack.bestDCAZ();
          if (useDCAzCut) {
            if (dcaZCut > maxDCAz)
              continue;
          }
        }
        if ((cfgnEta1 < track.eta()) && (track.eta() < cfgnEta2) && track.nClusters() >= cfgnCluster && retrack.ambDegree() > 0 && chi2ndf < cfgChi2NDFMax && (phi > cfgPhiCut1 && phi < cfgPhiCut2)) {
          registry.fill(HIST("Tracks/2Danalysis/EtaZvtx"), track.eta(), z);
        }
      }
      if (!disableITSROFCut && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
        return;
      }
      registry.fill(HIST("EventSelection"), static_cast<int>(EventSelectionBin::VzItsRof));
      if (!useEvSel || (useEvSel && collision.selection_bit(aod::evsel::kIsTriggerTVX) && collision.selection_bit(aod::evsel::kNoTimeFrameBorder) && collision.selection_bit(aod::evsel::kNoSameBunchPileup))) {
        registry.fill(HIST("EventSelection"), static_cast<int>(EventSelectionBin::VzSelected));
        registry.fill(HIST("Tracks/2Danalysis/EventsNtrkZvtx_sel8"), nTrk, z);
        std::unordered_set<int> uniqueEvents;
        std::unordered_set<int> uniqueEventsAmb;
        std::unordered_set<int> uniqueCollisions;
        std::unordered_set<int> uniqueCollisionsAmb;
        std::unordered_set<int> eventsInelMFT;
        std::unordered_set<int> eventsInel;
        if (midtracks.size() > 0) {
          registry.fill(HIST("EventSelection"), static_cast<int>(EventSelectionBin::Sel8VzInelGt0));
          registry.fill(HIST("EventsNtrkZvtx_gt0"), nTrk, z);
          registry.fill(HIST("Tracks/2Danalysis/EventsNtrkZvtx_sel8_inelgt0"), nTrk, z);
          eventsInel.insert(collision.globalIndex());
        }
        if (perCollisionSample.size() > 0) {
          registry.fill(HIST("EventSelection"), static_cast<int>(EventSelectionBin::PerCollisionSampleGt0));
        }
        if (midtracks.size() > 0 && perCollisionSample.size() > 0) {
          registry.fill(HIST("EventSelection"), static_cast<int>(EventSelectionBin::MidtracksAndPerCollisionSampleGt0));
        }
        int64_t i = 0.0, j = 0.0, k = 0.0;
        for (const auto& retrack : retracks) {
          auto track = retrack.mfttrack();
          float ndf = getTrackNdf(track);
          float chi2ndf = track.chi2() / ndf;
          float phi = track.phi();
          o2::math_utils::bringTo02Pi(phi);
          if (usePhiCut) {
            if ((phi <= PhiVetoLow) ||
                ((phi >= PhiVetoPiMin) && (phi <= PhiVetoPiMax)) ||
                (phi >= PhiVetoHigh))
              continue;
          }
          float dcaXyCut = retrack.bestDCAXY();
          if (useDCAxyCut) {
            if (dcaXyCut > maxDCAxy)
              continue;
          }
          if constexpr (std::is_same_v<RetracksT, soa::SmallGroups<aod::BestCollisionsFwd3d>>) {
            float dcaZCut = retrack.bestDCAZ();
            if (useDCAzCut) {
              if (dcaZCut > maxDCAz)
                continue;
            }
          }
          if ((cfgnEta1 < track.eta()) && (track.eta() < cfgnEta2) && track.nClusters() >= cfgnCluster && retrack.ambDegree() > 0 && chi2ndf < cfgChi2NDFMax && (phi > cfgPhiCut1 && phi < cfgPhiCut2)) {
            registry.fill(HIST("Tracks/Control/Chi2NDF"), chi2ndf);
            registry.fill(HIST("Tracks/2Danalysis/EtaZvtx_sel8"), track.eta(), z);
            if (midtracks.size() > 0 && retrack.ambDegree() > 0) {
              registry.fill(HIST("Tracks/2Danalysis/EtaZvtx_sel8_inelgt0"), track.eta(), z);
            }
          }
        }
        if (retracks.size() > 0) {
          registry.fill(HIST("EventSelection"), static_cast<int>(EventSelectionBin::SelInelInelFwdGt0));
          if (midtracks.size() > 0) {
            registry.fill(HIST("Tracks/2Danalysis/EventsNtrkZvtx_sel8_inelfwdgt0"), nTrk, z);
          }
          for (const auto& retrack : retracks) {
            auto track = retrack.mfttrack();
            float ndf = getTrackNdf(track);
            float chi2ndf = track.chi2() / ndf;
            float phi = track.phi();
            float dcaXyCut = retrack.bestDCAXY();
            o2::math_utils::bringTo02Pi(phi);
            // Declare dcaZCut only if needed below.
            if ((cfgnEta1 < track.eta()) && (track.eta() < cfgnEta2) && track.nClusters() >= cfgnCluster && chi2ndf < cfgChi2NDFMax && (phi > cfgPhiCut1 && phi < cfgPhiCut2)) {
              if (usePhiCut) {
                if ((phi <= PhiVetoLow) ||
                    ((phi >= PhiVetoPiMin) && (phi <= PhiVetoPiMax)) ||
                    (phi >= PhiVetoHigh))
                  continue;
              }
              if (useDCAxyCut) {
                if (dcaXyCut > maxDCAxy)
                  continue;
              }
              if constexpr (std::is_same_v<RetracksT, soa::SmallGroups<aod::BestCollisionsFwd3d>>) {
                float dcaZCut = retrack.bestDCAZ();
                if (useDCAzCut) {
                  if (dcaZCut > maxDCAz)
                    continue;
                }
              }

              registry.fill(HIST("TracksEtaZvtx"), track.eta(), z);
              if (midtracks.size() > 0 && retrack.ambDegree() > 0) {
                registry.fill(HIST("Tracks/EtaZvtx_gt0"), track.eta(), z);
                registry.fill(HIST("Tracks/2Danalysis/EtaZvtx_sel8_inelfwdgt0"), track.eta(), z);
                eventsInelMFT.insert(retrack.bestCollisionId());
              }
              if (retrack.ambDegree() != 0) {
                registry.fill(HIST("Tracks/Control/woOrp/woOrpEtaZvtx_gt0"), track.eta(), z);
                ++k;
              }
              float phi = track.phi();
              o2::math_utils::bringTo02Pi(phi);
              registry.fill(HIST("TracksPhiEta"), phi, track.eta());
              registry.fill(HIST("TracksPtEta"), track.pt(), track.eta());
              if ((track.eta() < ForwardEtaMax) && (track.eta() > ForwardEtaMin)) {
                registry.fill(HIST("TracksPhiZvtx"), phi, z);
              }
              if (track.collisionId() > InvalidCollisionId &&
                  retrack.ambDegree() == SingleCompatibleCollision) {
                registry.fill(HIST("collisionID"), track.collisionId());
              }
              if (track.collisionId() > InvalidCollisionId &&
                  retrack.ambDegree() > SingleCompatibleCollision) {
                registry.fill(HIST("collisionIDamb"), track.collisionId());
              }
              if (track.collisionId() != retrack.bestCollisionId()) {
                registry.fill(HIST("Tracks/Control/ReassignedTracksEtaZvtx"),
                              track.eta(), z);
                registry.fill(HIST("Tracks/Control/ReassignedTracksPhiEta"), phi,
                              track.eta());
                registry.fill(HIST("Tracks/Control/ReassignedVertexCorr"),
                              track.template collision_as<CollwEv>().posZ(), z);

                registry.fill(HIST("Tracks/Control/DeltaZ"),
                              track.template collision_as<CollwEv>().posZ() -
                                collision.posZ());
              }
              if (track.collisionId() == retrack.bestCollisionId()) {
                registry.fill(HIST("Tracks/Control/notReassignedTracksEtaZvtx"),
                              track.eta(), z);
                registry.fill(HIST("Tracks/Control/notReassignedTracksPhiEta"), phi,
                              track.eta());
                registry.fill(HIST("Tracks/Control/notReassignedVertexCorr"),
                              track.template collision_as<CollwEv>().posZ(), z);
              }

              registry.fill(HIST("Tracks/Control/TrackAmbDegree"),
                            retrack.ambDegree());
              registry.fill(HIST("Tracks/Control/DCAXY"), retrack.bestDCAXY());
              if constexpr (std::is_same_v<RetracksT, soa::SmallGroups<aod::BestCollisionsFwd3d>>) {
                registry.fill(HIST("Tracks/Control/DCAZ"), retrack.bestDCAZ());
              }
              int isAmbiguous = 0;

              if (retrack.ambDegree() > 1 && retrack.ambDegree() != 0) {
                isAmbiguous = 1;
                ++i;

                registry.fill(HIST("Tracks/Control/amb/EtaZvtxAmb_gt0"), track.eta(), z);

                registry.fill(HIST("Tracks/Control/amb/AmbTracksEtaZvtx"),
                              track.eta(), z);
                registry.fill(HIST("Tracks/Control/amb/AmbTracksPhiEta"), phi,
                              track.eta());
                registry.fill(HIST("Tracks/Control/amb/AmbVertexCorr"),
                              track.template collision_as<CollwEv>().posZ(), z);
                registry.fill(HIST("Tracks/Control/amb/DCAxy_amb"), retrack.bestDCAXY());
                if constexpr (std::is_same_v<RetracksT, soa::SmallGroups<aod::BestCollisionsFwd3d>>) {
                  registry.fill(HIST("Tracks/Control/amb/DCAz_amb"), retrack.bestDCAZ());
                }
                if (track.collisionId() == retrack.bestCollisionId()) {
                }
                uniqueEventsAmb.insert(retrack.bestCollisionId());
              }
              if (midtracks.size() > 0 && retrack.ambDegree() > 1 && retrack.ambDegree() != 0) {
                uniqueCollisionsAmb.insert(collision.globalIndex());
              }

              registry.fill(HIST("Tracks/Control/TrackIsAmb"), isAmbiguous);
              if (retrack.ambDegree() == 1 && retrack.ambDegree() != 0) {
                ++j;
                registry.fill(HIST("Tracks/Control/nonamb/EtaZvtxNonAmb_gt0"), track.eta(), z);
                registry.fill(HIST("Tracks/Control/nonamb/nonAmbTracksEtaZvtx"),
                              track.eta(), z);
                registry.fill(HIST("Tracks/Control/nonamb/nonAmbTracksPhiEta"), phi,
                              track.eta());
                registry.fill(HIST("Tracks/Control/nonamb/nonAmbVertexCorr"),
                              track.template collision_as<CollwEv>().posZ(), z);
                registry.fill(HIST("Tracks/Control/nonamb/DCAxy_nonamb"), retrack.bestDCAXY());
                if constexpr (std::is_same_v<RetracksT, soa::SmallGroups<aod::BestCollisionsFwd3d>>) {
                  registry.fill(HIST("Tracks/Control/nonamb/DCAz_nonamb"), retrack.bestDCAZ());
                }
                if (track.collisionId() == retrack.bestCollisionId()) {
                }
                uniqueEvents.insert(retrack.bestCollisionId());
              }
              if (midtracks.size() > 0 && retrack.ambDegree() == 1 && retrack.ambDegree() != 0) {
                uniqueCollisions.insert(collision.globalIndex());
              }
              if ((retrack.ambDegree() > 1) || (retrack.ambDegree() <= 1))

                if (retrack.ambDegree() != 0) {
                  registry.fill(HIST("Tracks/Control/woOrp/woOrpTracksEtaZvtx"),
                                track.eta(), z);
                  registry.fill(HIST("Tracks/Control/woOrp/woOrpTracksPhiEta"), phi,
                                track.eta());
                  registry.fill(HIST("Tracks/Control/woOrp/woOrpVertexCorr"),
                                track.template collision_as<CollwEv>().posZ(), z);
                  registry.fill(HIST("Tracks/Control/woOrp/DCAxy_woOrp"), retrack.bestDCAXY());
                  if constexpr (std::is_same_v<RetracksT, soa::SmallGroups<aod::BestCollisionsFwd3d>>) {
                    registry.fill(HIST("Tracks/Control/woOrp/DCAz_woOrp"), retrack.bestDCAZ());
                  }
                }
            }
          }
          registry.fill(HIST("ambEventCounts"), 1, uniqueEventsAmb.size());
          registry.fill(HIST("NonambEventCounts"), 1, uniqueEvents.size());
          registry.fill(HIST("hNumCollisionsNonAmb_InelMFT"), 1, uniqueCollisions.size());
          registry.fill(HIST("hNumCollisionsAmb_InelMFT"), 1, uniqueCollisionsAmb.size());
          registry.fill(HIST("hNumCollisions_InelMFT"), 1, eventsInelMFT.size());
        }
        registry.fill(HIST("Tracks/Control/amb/nTrkAmb"), i);
        registry.fill(HIST("Tracks/Control/nonamb/nTrkNonAmb"), j);
        registry.fill(HIST("Tracks/Control/woOrp/nTrk"), k);
        registry.fill(HIST("hNumCollisions_Inel"), 1, eventsInel.size());
      }
    } else {
      registry.fill(HIST("EventSelection"), static_cast<int>(EventSelectionBin::Rejected));
    }
  }

  void processMultReassoc(CollwEv::iterator const& collision,
                          o2::aod::MFTTracks const& mft,
                          soa::SmallGroups<aod::BestCollisionsFwd> const& retracks,
                          FiCentralTracks const& midtracks, aod::Tracks const& trk)
  {
    processMultReassocCommon(collision, mft, retracks, midtracks, trk);
  }

  void processMultReassoc3d(CollwEv::iterator const& collision,
                            o2::aod::MFTTracks const& mft,
                            soa::SmallGroups<aod::BestCollisionsFwd3d> const& retracks,
                            FiCentralTracks const& midtracks, aod::Tracks const& trk)
  {
    processMultReassocCommon(collision, mft, retracks, midtracks, trk);
  }
  PROCESS_SWITCH(PseudorapidityDensityMFT, processMultReassoc,
                 "Process reco or data info", false);

  PROCESS_SWITCH(PseudorapidityDensityMFT, processMultReassoc3d,
                 "Process reco or data info (3d)", false);

  using ExColsCent = soa::Join<aod::Collisions, aod::CentFT0Cs, aod::EvSels>;

  void processCountingCentrality(ExColsCent::iterator const& collision,
                                 aod::MFTTracks const& tracks)
  {
    auto c = collision.centFT0C();
    registry.fill(HIST("Events/Centrality/Selection"), 1., c);

    if (!useEvSel || collision.sel8()) {
      auto z = collision.posZ();
      registry.fill(HIST("Events/Centrality/Selection"), 2., c);
      auto perCollisionSample = sample->sliceByCached(
        o2::aod::fwdtrack::collisionId, collision.globalIndex(), cache);
      auto nTrk = perCollisionSample.size();

      registry.fill(HIST("Events/Centrality/NtrkZvtx"), nTrk, z, c);

      for (const auto& track : tracks) {

        float phi = track.phi();
        o2::math_utils::bringTo02Pi(phi);

        if (usePhiCut) {
          if ((phi < cfgPhiCut) ||
              ((phi > o2::constants::math::PI - cfgPhiCut) && (phi < o2::constants::math::PI + cfgPhiCut)) ||
              (phi > o2::constants::math::TwoPI - cfgPhiCut) ||
              ((phi > ((o2::constants::math::PIHalf - 0.1) * o2::constants::math::PI) - cfgPhiCut) &&
               (phi < ((o2::constants::math::PIHalf - 0.1) * o2::constants::math::PI) + cfgPhiCut)))
            continue;
        }

        registry.fill(HIST("Tracks/Centrality/EtaZvtx"), track.eta(), z, c);
        registry.fill(HIST("Tracks/Centrality/PhiEta"), phi, track.eta(), c);
      }

    } else {
      registry.fill(HIST("Events/Centrality/Selection"), 3.,
                    c); // rejected events
    }
  }

  PROCESS_SWITCH(PseudorapidityDensityMFT, processCountingCentrality,
                 "Count tracks in centrality bins", false);

  using Particles = soa::Filtered<aod::McParticles>;
  expressions::Filter primaries =
    (aod::mcparticle::flags &
     (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary) ==
    (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary;
  static constexpr float McSampleEtaMax = 1.1f;
  Partition<Particles> mcSample = nabs(aod::mcparticle::eta) < McSampleEtaMax;
  Partition<Particles> mcSampleCentral =
    nabs(aod::mcparticle::eta) < estimatorEta;

  void processGen(
    aod::McCollisions::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<aod::Collisions, aod::EvSels,
                                   aod::McCollisionLabels>> const& collisions,
    Particles const& particles, aod::MFTTracks const& /*tracks*/,
    FiCentralTracks const& midtracks)
  {
    registry.fill(HIST("EventEfficiency"), static_cast<int>(EventEfficiencyBin::Generated));

    auto perCollisionMCSample = mcSample->sliceByCached(
      aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);
    auto nCharged = 0;
    for (const auto& particle : perCollisionMCSample) {
      auto charge = 0.;
      auto p = pdg->GetParticle(particle.pdgCode());
      if (p != nullptr) {
        charge = p->Charge();
      }
      if (std::abs(charge) < ChargeUnitTimesThree) {
        continue;
      }
      nCharged++;
    }
    registry.fill(HIST("EventsNtrkZvtxGen_t"), nCharged, mcCollision.posZ());

    //--------for INEL>0
    auto perCollisionMCSampleCentral = mcSampleCentral->sliceByCached(
      aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);
    auto nChargedCentral = 0;
    for (const auto& particle : perCollisionMCSample) {
      auto charge = 0.;
      auto p = pdg->GetParticle(particle.pdgCode());
      if (p != nullptr) {
        charge = p->Charge();
      }
      if (std::abs(charge) < ChargeUnitTimesThree) {
        continue;
      }
      nChargedCentral++;
    }
    if ((mcCollision.posZ() >= cfgVzCut1) && (mcCollision.posZ() <= cfgVzCut2)) {
      if (nChargedCentral > 0) {
        registry.fill(HIST("EventEfficiency"), static_cast<int>(EventEfficiencyBin::GeneratedInelGt0));
        registry.fill(HIST("EventsNtrkZvtxGen_gt0t"), nCharged,
                      mcCollision.posZ());
      }
    }
    //-----------
    bool atLeastOne = false;
    bool atLeastOneGt0 = false;
    int moreThanOne = 0;

    LOGP(debug, "MC col {} has {} reco cols", mcCollision.globalIndex(),
         collisions.size());
    for (const auto& collision : collisions) {
      registry.fill(HIST("EventEfficiency"), static_cast<int>(EventEfficiencyBin::Reconstructed));
      if (!disableITSROFCut && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
        return;
      }
      if (!useEvSel || (useEvSel && collision.selection_bit(aod::evsel::kIsTriggerTVX) && collision.selection_bit(aod::evsel::kNoTimeFrameBorder) && collision.selection_bit(aod::evsel::kNoSameBunchPileup))) {
        atLeastOne = true;
        auto perCollisionSample = sample->sliceByCached(
          o2::aod::fwdtrack::collisionId, collision.globalIndex(), cache);

        registry.fill(HIST("EventEfficiency"), static_cast<int>(EventEfficiencyBin::Selected));
        auto perCollisionSampleCentral =
          midtracks.sliceBy(perColCentral, collision.globalIndex());
        if ((collision.posZ() >= cfgVzCut1) && (collision.posZ() <= cfgVzCut2) && (mcCollision.posZ() >= cfgVzCut1) && (mcCollision.posZ() <= cfgVzCut2)) {
          if (perCollisionSampleCentral.size() > 0) {
            registry.fill(HIST("EventEfficiency"), static_cast<int>(EventEfficiencyBin::SelectedInelGt0));
            atLeastOneGt0 = true;
            registry.fill(HIST("EventsNtrkZvtxGen_gt0"),
                          perCollisionSample.size(), collision.posZ());
          }

          registry.fill(HIST("EventsZposDiff"),
                        collision.posZ() - mcCollision.posZ());
          if (useZDiffCut) {
            if (std::abs(collision.posZ() - mcCollision.posZ()) > maxZDiff) {
              continue;
            }
          }
          registry.fill(HIST("EventsNtrkZvtxGen"), perCollisionSample.size(),
                        collision.posZ());
          ++moreThanOne;
        }
      }
    }
    if (collisions.size() == 0) {
      registry.fill(HIST("NotFoundEventZvtx"), mcCollision.posZ());
    }
    if (moreThanOne > 1) {
      registry.fill(HIST("EventsSplitMult"), nCharged);
    }
    if ((mcCollision.posZ() >= cfgVzCut1) && (mcCollision.posZ() <= cfgVzCut2)) {

      for (const auto& particle : particles) {
        auto p = pdg->GetParticle(particle.pdgCode());
        auto charge = 0;
        if (p != nullptr) {
          charge = static_cast<int>(p->Charge());
        }
        if (std::abs(charge) < ChargeUnitTimesThree) {
          continue;
        }
        float phi = particle.phi();
        o2::math_utils::bringTo02Pi(phi);
        float ptCut = particle.pt();

        if (usePhiCut) {
          if ((phi <= PhiVetoLow) ||
              ((phi >= PhiVetoPiMin) && (phi <= PhiVetoPiMax)) ||
              (phi >= PhiVetoHigh))
            continue;
        }
        if (usePtCut) {
          if (ptCut > cfgnPt)
            continue;
        }
        if (cfgnEta1 < particle.eta() && particle.eta() < cfgnEta2 && (phi > cfgPhiCut1 && phi < cfgPhiCut2)) {
          registry.fill(HIST("TracksEtaZvtxGen_t"), particle.eta(),
                        mcCollision.posZ());
          if (perCollisionMCSampleCentral.size() > 0) {
            registry.fill(HIST("TracksEtaZvtxGen_gt0t"), particle.eta(),
                          mcCollision.posZ());
            registry.fill(HIST("TracksPhiEtaGen_gt0t"), particle.phi(), particle.eta());
          }
          if (atLeastOne) {
            registry.fill(HIST("TracksEtaZvtxGen"), particle.eta(),
                          mcCollision.posZ());
            registry.fill(HIST("TracksPtEtaGen"), particle.pt(), particle.eta());
            if (atLeastOneGt0) {
              registry.fill(HIST("TracksEtaZvtxGen_gt0"), particle.eta(),
                            mcCollision.posZ());
              registry.fill(HIST("TracksPhiEtaGen_gt0"), particle.phi(), particle.eta());
              if (particle.isPhysicalPrimary()) {
                registry.fill(HIST("TracksEtaZvtxGen_gt0_primary"), particle.eta(), mcCollision.posZ());
                registry.fill(HIST("TracksPhiEtaGen_gt0_primary"), particle.phi(), particle.eta());
                registry.fill(HIST("TracksPtZvtxGen_gt0_primary"), particle.pt(), mcCollision.posZ());
              }
            }
          }

          registry.fill(HIST("TracksPhiEtaGen"), particle.phi(), particle.eta());
          registry.fill(HIST("TracksPhiZvtxGen"), particle.phi(),
                        mcCollision.posZ());
          registry.fill(HIST("TracksPtEtaGen_t"), particle.pt(), particle.eta());
        }
      }
    }
  }

  PROCESS_SWITCH(PseudorapidityDensityMFT, processGen,
                 "Process generator-level info", false);

  using ExColsGenCent =
    soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions,
                               aod::CentFT0Cs, aod::EvSels>>;

  void processGenCent(aod::McCollisions::iterator const& mcCollision,
                      ExColsGenCent const& collisions,
                      Particles const& particles,
                      MFTTracksLabeled const& /*tracks*/)
  {

    LOGP(debug, "MC col {} has {} reco cols", mcCollision.globalIndex(),
         collisions.size());

    float cGen = -1;
    bool atLeastOne = false;
    for (const auto& collision : collisions) {
      float cRec = -1;
      if constexpr (ExColsGenCent::template contains<aod::CentFT0Cs>()) {
        cRec = collision.centFT0C();
      }
      if (!useEvSel || (useEvSel && collision.sel8())) {
        if constexpr (ExColsGenCent::template contains<aod::CentFT0Cs>()) {
          if (!atLeastOne) {
            cGen = cRec;
          }
        }
        atLeastOne = true;

        registry.fill(HIST("Events/Centrality/EventEfficiency"), 2., cGen);
        registry.fill(HIST("Events/Centrality/CentPercentileMCGen"), cGen);

        auto perCollisionSample = sample->sliceByCached(
          o2::aod::fwdtrack::collisionId, collision.globalIndex(), cache);
        registry.fill(HIST("Events/Centrality/NtrkZvtxGen"),
                      perCollisionSample.size(), collision.posZ(), cGen);
      }
    }

    registry.fill(HIST("Events/Centrality/EventEfficiency"), 1., cGen);

    auto perCollisionMCSample = mcSample->sliceByCached(
      aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);
    auto nCharged = 0;

    for (const auto& particle : perCollisionMCSample) {
      auto p = pdg->GetParticle(particle.pdgCode());
      auto charge = 0;
      if (p != nullptr) {
        charge = static_cast<int>(p->Charge());
      }
      if (std::abs(charge) < ChargeUnitTimesThree) {
        continue;
      }
      nCharged++;
    }

    if constexpr (ExColsGenCent::template contains<aod::CentFT0Cs>()) {
      registry.fill(HIST("Events/Centrality/NtrkZvtxGen_t"), nCharged,
                    mcCollision.posZ(), cGen);
    }

    for (const auto& particle : particles) {
      auto p = pdg->GetParticle(particle.pdgCode());
      auto charge = 0;
      if (p != nullptr) {
        charge = static_cast<int>(p->Charge());
      }
      if (std::abs(charge) < ChargeUnitTimesThree) {
        continue;
      }

      if constexpr (ExColsGenCent::template contains<aod::CentFT0Cs>()) {
        registry.fill(HIST("Tracks/Centrality/EtaZvtxGen_t"), particle.eta(),
                      mcCollision.posZ(), cGen);
      }

      if (atLeastOne) {
        if constexpr (ExColsGenCent::template contains<aod::CentFT0Cs>()) {
          registry.fill(HIST("Tracks/Centrality/EtaZvtxGen"), particle.eta(),
                        mcCollision.posZ(), cGen);
          float phi = particle.phi();
          o2::math_utils::bringTo02Pi(phi);
          registry.fill(HIST("Tracks/Centrality/PhiEtaGen"), phi,
                        particle.eta(), cGen);
        }
      }
    }
  }

  PROCESS_SWITCH(PseudorapidityDensityMFT, processGenCent,
                 "Process generator-level info in centrality bins", false);

  void processGenPt(
    soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
    MFTTracksLabeled const& tracks, aod::McParticles const&)
  {
    if (!useEvSel || (useEvSel && collision.sel8())) {
      for (const auto& track : tracks) {
        if (!track.has_mcParticle()) {
          continue;
        }
        auto particle = track.mcParticle();
        if (!particle.isPhysicalPrimary()) {
          continue;
        }
        registry.fill(HIST("TracksToPartPtEta"), particle.pt(), particle.eta());
      }
    }
  }

  PROCESS_SWITCH(PseudorapidityDensityMFT, processGenPt,
                 "Process particle-level info of pt", false);

  void processGenRecoTimeCom(McCollisionsWithExtra::iterator const& mcCollision,
                             o2::soa::SmallGroups<soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>> const& collisions,
                             FullBCs const& bcs,
                             MFTTracksLabeledOrg const& tracks,
                             FiCentralTracks const& midtracks,
                             aod::McParticles const&)
  {
    const auto fillGenRecoCut = [&](GenRecoCutBin bin) {
      registry.fill(HIST("EventsRecoCuts_GenReco"), static_cast<int>(bin));
    };
    fillGenRecoCut(GenRecoCutBin::AllRecoCollisions);
    std::unordered_map<int, int> recoToMc;
    std::unordered_map<int, std::vector<int>> mcToReco; // MC collision id -> list of reco collision globalIndex

    for (const auto& collision : collisions) {
      int nSavedRows = 0;
      std::unordered_set<int> uniqueRecoColsSaved;
      int recoCol = collision.globalIndex(); // reconstructed vertex index
      int mcCol = collision.mcCollisionId(); // true MC collision index

      if (mcCol >= 0) {
        recoToMc[recoCol] = mcCol;
        mcToReco[mcCol].push_back(recoCol);

        ++nSavedRows;
        uniqueRecoColsSaved.insert(recoCol);
      }
      registry.fill(HIST("Purity/HashTableRowCounts"),
                    static_cast<int>(HashTableRowCountsBin::RowsSaved), nSavedRows);
      registry.fill(HIST("Purity/HashTableRowCounts"),
                    static_cast<int>(HashTableRowCountsBin::UniqueRecoColsSaved), uniqueRecoColsSaved.size());

      registry.fill(HIST("Purity/reco/CollisionNumContrib"), collision.numContrib());

      if (useCont && collision.globalIndex() != mcCollision.bestCollisionIndex()) {
        continue;
      }
      fillGenRecoCut(GenRecoCutBin::UseContBestCollisionIndex);

      if (!collision.has_mcCollision()) {
        LOGF(warning, "No MC collision found...");
        return;
      }
      fillGenRecoCut(GenRecoCutBin::HasMcCollision);

      auto countAndPassEvSelGenReco = [&](auto const& collision) {
        struct EvSelStep {
          bool enabled;
          decltype(aod::evsel::kIsTriggerTVX) bit;
          GenRecoCutBin bin;
        };

        const std::array<EvSelStep, 10> steps = {{
          {true, aod::evsel::kIsTriggerTVX, GenRecoCutBin::IsTriggerTVX},
          {true, aod::evsel::kNoTimeFrameBorder, GenRecoCutBin::NoTimeFrameBorder},
          {true, aod::evsel::kNoITSROFrameBorder, GenRecoCutBin::NoITSROFrameBorder},
          {useNoSameBunchPileup, aod::evsel::kNoSameBunchPileup, GenRecoCutBin::NoSameBunchPileup},
          {useGoodZvtxFT0vsPV, aod::evsel::kIsGoodZvtxFT0vsPV, GenRecoCutBin::GoodZvtxFT0vsPV},
          {useNoCollInRofStandard, aod::evsel::kNoCollInRofStandard, GenRecoCutBin::NoCollInRofStandard},
          {useNoCollInRofStrict, aod::evsel::kNoCollInRofStrict, GenRecoCutBin::NoCollInRofStrict},
          {useNoCollInTimeRangeStandard, aod::evsel::kNoCollInTimeRangeStandard, GenRecoCutBin::NoCollInTimeRangeStandard},
          {useNoCollInTimeRangeStrict, aod::evsel::kNoCollInTimeRangeStrict, GenRecoCutBin::NoCollInTimeRangeStrict},
          {useNoHighMultCollInPrevRof, aod::evsel::kNoHighMultCollInPrevRof, GenRecoCutBin::NoHighMultCollInPrevRof},
        }};

        if (!useEvSel) {
          for (const auto& step : steps) {
            fillGenRecoCut(step.bin);
          }
          fillGenRecoCut(GenRecoCutBin::RctMFT);
          return true;
        }

        for (const auto& step : steps) {
          if (!step.enabled) {
            fillGenRecoCut(step.bin);
            continue;
          }

          if (!collision.selection_bit(step.bit)) {
            return false;
          }
          fillGenRecoCut(step.bin);
        }

        if (useRctMFT && !myChecker(collision)) {
          return false;
        }
        fillGenRecoCut(GenRecoCutBin::RctMFT);

        return true;
      };

      if (!countAndPassEvSelGenReco(collision)) {
        continue;
      }

      const auto z = collision.posZ();
      if ((z < cfgVzCut1) || (z > cfgVzCut2)) {
        continue;
      }
      fillGenRecoCut(GenRecoCutBin::VzWindow);

      auto perCollisionSampleCentral = midtracks.sliceBy(perColCentral, collision.globalIndex());
      if (perCollisionSampleCentral.size() <= 0) {
        continue;
      }
      fillGenRecoCut(GenRecoCutBin::InelGt0);

      // constexpr uint8_t kFakeMcMask = 1u << 7;
      for (const auto& track : tracks) {
        float ndf = getTrackNdf(track);
        const float chi2ndf = track.chi2() / ndf;
        float phi = track.phi();
        // const float dcaXyCut = track.bestDCAXY();
        // const float dcaZCut = track.bestDCAZ();
        const float ptCut = track.pt();
        o2::math_utils::bringTo02Pi(phi);

        const bool failTrackCuts =
          track.nClusters() < cfgnCluster ||
          track.eta() <= cfgnEta1 ||
          track.eta() >= cfgnEta2 ||
          chi2ndf >= cfgChi2NDFMax ||
          phi <= cfgPhiCut1 ||
          phi >= cfgPhiCut2 ||
          (usePhiCut &&
           ((phi <= PhiVetoLow) ||
            ((phi >= PhiVetoPiMin) && (phi <= PhiVetoPiMax)) ||
            (phi >= PhiVetoHigh))) ||
          // (useDCAxyCut && dcaxyCut > maxDCAxy) ||
          // (useDCAzCut && std::abs(dcazCut) > maxDCAz) ||
          (usePtCut && ptCut > cfgnPt);

        if (failTrackCuts) {
          continue;
        }
        const bool hasMcLabel = track.has_mcParticle();
        const bool isFakeByLabel = hasMcLabel ? (track.mcMask() != 0) : false;
        const bool isTrueByLabel = hasMcLabel && !isFakeByLabel;
        const bool hasNoMcLabel = !hasMcLabel;
        const bool isPrimaryCharged = hasMcLabel && !isFakeByLabel && track.mcParticle().isPhysicalPrimary();
        const bool isSecondaryCharged = hasMcLabel && !isFakeByLabel && !track.mcParticle().isPhysicalPrimary();
        const float eta = track.eta();
        if (!passGenRecoTrackMode(track)) { // 0-> All nonorphans, 1->Non-Amb, 2->Amb
          continue;
        }
        int bin = static_cast<int>(RightWrongBin::Neither);
        bool recoOfTrueExists = false;

        if (isTrueByLabel) {
          int recoCol = track.collisionId();
          auto itRecoToMc = recoToMc.find(recoCol);
          const int mcOfTrack = track.mcParticle().mcCollisionId();

          // Check whether any reco vertex exists for the true MC collision of this track
          for (const auto& [recoId, mcId] : recoToMc) {
            if (mcId == mcOfTrack) {
              recoOfTrueExists = true;
              break;
            }
          }

          if (recoCol >= 0 && itRecoToMc != recoToMc.end()) {
            int mcFromReco = itRecoToMc->second;
            bin = (mcFromReco == mcOfTrack)
                    ? static_cast<int>(RightWrongBin::Right)
                    : static_cast<int>(RightWrongBin::Wrong);
          }
        }

        registry.fill(HIST("RightWrong"), bin);

        if (bin == static_cast<int>(RightWrongBin::Wrong)) {
          registry.fill(HIST("Purity/WrongVertexRecoExists"), recoOfTrueExists ? static_cast<int>(BoolBin::Yes) : static_cast<int>(BoolBin::No));
        }

        const auto fillTrackLabelSummary = [&](TrackLabelSummaryBin bin) {
          registry.fill(HIST("Purity/TrackLabelSummary"), static_cast<int>(bin));
        };
        const auto fillTrackEtaCategory = [&](TrackLabelSummaryBin bin) {
          constexpr float EtaSentinel = -999.f;

          float etaAll = EtaSentinel;
          float etaNoMc = EtaSentinel;
          float etaFake = EtaSentinel;
          float etaTrue = EtaSentinel;
          float etaPrimary = EtaSentinel;
          float etaSecondary = EtaSentinel;

          switch (bin) {
            case TrackLabelSummaryBin::AllTracks:
              etaAll = eta;
              break;
            case TrackLabelSummaryBin::NoMcLabel:
              etaNoMc = eta;
              break;
            case TrackLabelSummaryBin::FakeTracks:
              etaFake = eta;
              break;
            case TrackLabelSummaryBin::TrueTracks:
              etaTrue = eta;
              break;
            case TrackLabelSummaryBin::PrimaryTracks:
              etaPrimary = eta;
              break;
            case TrackLabelSummaryBin::SecondaryTracks:
              etaSecondary = eta;
              break;
          }

          registry.fill(HIST("Purity/TrackEtaCategorySparse"),
                        etaAll, etaNoMc, etaFake, etaTrue, etaPrimary, etaSecondary);
        };

        // registry.fill(HIST("Purity/TrackLabelStatus"), hasMcLabel ? 1.0 : 0.0);
        // registry.fill(HIST("Purity/TrackFakeStatus"), isFakeByLabel ? 1.0 : 0.0);

        fillTrackLabelSummary(TrackLabelSummaryBin::AllTracks);
        fillTrackEtaCategory(TrackLabelSummaryBin::AllTracks);

        if (hasNoMcLabel) {
          fillTrackLabelSummary(TrackLabelSummaryBin::NoMcLabel);
          fillTrackEtaCategory(TrackLabelSummaryBin::NoMcLabel);
          continue;
        }

        if (isFakeByLabel) {
          fillTrackLabelSummary(TrackLabelSummaryBin::FakeTracks);
          fillTrackEtaCategory(TrackLabelSummaryBin::FakeTracks);
          continue;
        }

        fillTrackLabelSummary(TrackLabelSummaryBin::TrueTracks);
        fillTrackEtaCategory(TrackLabelSummaryBin::TrueTracks);

        if (isPrimaryCharged) {
          fillTrackLabelSummary(TrackLabelSummaryBin::PrimaryTracks);
          fillTrackEtaCategory(TrackLabelSummaryBin::PrimaryTracks);
        }

        if (isSecondaryCharged) {
          fillTrackLabelSummary(TrackLabelSummaryBin::SecondaryTracks);
          fillTrackEtaCategory(TrackLabelSummaryBin::SecondaryTracks);
        }
        // registry.fill(HIST("RightWrong"), bin);

      } // Track loop 1
    } // Collision
  }
  PROCESS_SWITCH(PseudorapidityDensityMFT, processGenRecoTimeCom,
                 "Process for MC time compatible", false);

  // using MFTTracksLabeled = soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels>;
  // using MFTTracksLabeled =soa::Join<o2::aod::MFTTracks,aod::BestCollisionsFwd3d,aod::McMFTTrackLabels>;

  // aod::MFTTracks const& tracks
  // soa::Join<aod::McCollisions, aod::McCollsExtra>::iterator const& mcCollision
  // aod::McCollisions::iterator const& mcCollision
  // McCollisionsWithExtra::iterator const& mcCollision

  // void processMCeff(soa::Join<aod::McCollisions, aod::McCollsExtra>::iterator const& mcCollision, CollisionMCRecTable const& RecCols, TrackMCTrueTable const& GenParticles, FilTrackMCRecTable const& RecTracks)
  // soa::Join<aod::McCollisions, aod::McCollsExtra>::iterator const& mcCollision //This worked
  void processGenReco(McCollisionsWithExtra::iterator const& mcCollision,
                      o2::soa::SmallGroups<soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>> const& collisions,
                      FullBCs const& bcs,
                      MFTTracksLabeled const& tracks,
                      FiCentralTracks const& midtracks,
                      aod::McParticles const&)

  {

    const auto fillGenRecoCut = [&](GenRecoCutBin bin) {
      registry.fill(HIST("EventsRecoCuts_GenReco"), static_cast<int>(bin));
    };
    fillGenRecoCut(GenRecoCutBin::AllRecoCollisions);

    // if (maxGenRecoEvents >= 0 && nProcessedGenReco >= maxGenRecoEvents) {
    //   return;
    // }
    // ++nProcessedGenReco; // for HIR
    std::unordered_map<int, int> recoToMc;
    std::unordered_map<int, std::vector<int>> mcToReco; // MC collision id -> list of reco collision globalIndex
    std::unordered_map<float, float> recoToMcVZ;
    std::unordered_map<int, float> recoVtxX;
    std::unordered_map<int, float> recoVtxY;
    std::unordered_map<int, float> recoVtxZ;
    std::unordered_map<int, std::array<double, 3>> recoVtxByRecoId;
    std::unordered_map<int, std::array<double, 3>> recoVtxByMcId;

    recoVtxByRecoId.reserve(collisions.size());
    recoVtxByMcId.reserve(collisions.size());
    mcToReco.reserve(collisions.size());

    // --- Make sure magnetic field exists in THIS device before any propagation ---
    // IMPORTANT: calling collision.bc_as<FullBCs>() requires the BC table to be subscribed.
    // We subscribe by taking `FullBCs const& bcs` in the process signature and init once here.
    bool magInited = false;
    for (auto const& bc : bcs) {
      initMagField(bc);
      magInited = true;
      break; // once is enough (initMagField is internally guarded)
    }
    if (!magInited) {
      LOGF(fatal, "BC table is empty: cannot initialize magnetic field");
    }

    //_______________________________________________________________________________

    for (const auto& collision : collisions) {
      int nSavedRows = 0;
      std::unordered_set<int> uniqueRecoColsSaved;
      int recoCol = collision.globalIndex(); // reconstructed vertex index
      int mcCol = collision.mcCollisionId(); // true MC collision index

      if (mcCol >= 0) {
        recoToMc[recoCol] = mcCol;
        mcToReco[mcCol].push_back(recoCol);

        ++nSavedRows;
        uniqueRecoColsSaved.insert(recoCol);
      }
      registry.fill(HIST("Purity/HashTableRowCounts"),
                    static_cast<int>(HashTableRowCountsBin::RowsSaved), nSavedRows);
      registry.fill(HIST("Purity/HashTableRowCounts"),
                    static_cast<int>(HashTableRowCountsBin::UniqueRecoColsSaved), uniqueRecoColsSaved.size());

      recoVtxX[recoCol] = collision.posX();
      recoVtxY[recoCol] = collision.posY();
      recoVtxZ[recoCol] = collision.posZ();

      registry.fill(HIST("Purity/reco/CollisionNumContrib"), collision.numContrib());

      if (useCont && collision.globalIndex() != mcCollision.bestCollisionIndex()) {
        continue;
      }
      fillGenRecoCut(GenRecoCutBin::UseContBestCollisionIndex);

      if (!collision.has_mcCollision()) {
        LOGF(warning, "No MC collision found...");
        return;
      }
      fillGenRecoCut(GenRecoCutBin::HasMcCollision);

      auto countAndPassEvSelGenReco = [&](auto const& collision) {
        struct EvSelStep {
          bool enabled;
          decltype(aod::evsel::kIsTriggerTVX) bit;
          GenRecoCutBin bin;
        };

        const std::array<EvSelStep, 10> steps = {{
          {true, aod::evsel::kIsTriggerTVX, GenRecoCutBin::IsTriggerTVX},
          {true, aod::evsel::kNoTimeFrameBorder, GenRecoCutBin::NoTimeFrameBorder},
          {true, aod::evsel::kNoITSROFrameBorder, GenRecoCutBin::NoITSROFrameBorder},
          {useNoSameBunchPileup, aod::evsel::kNoSameBunchPileup, GenRecoCutBin::NoSameBunchPileup},
          {useGoodZvtxFT0vsPV, aod::evsel::kIsGoodZvtxFT0vsPV, GenRecoCutBin::GoodZvtxFT0vsPV},
          {useNoCollInRofStandard, aod::evsel::kNoCollInRofStandard, GenRecoCutBin::NoCollInRofStandard},
          {useNoCollInRofStrict, aod::evsel::kNoCollInRofStrict, GenRecoCutBin::NoCollInRofStrict},
          {useNoCollInTimeRangeStandard, aod::evsel::kNoCollInTimeRangeStandard, GenRecoCutBin::NoCollInTimeRangeStandard},
          {useNoCollInTimeRangeStrict, aod::evsel::kNoCollInTimeRangeStrict, GenRecoCutBin::NoCollInTimeRangeStrict},
          {useNoHighMultCollInPrevRof, aod::evsel::kNoHighMultCollInPrevRof, GenRecoCutBin::NoHighMultCollInPrevRof},
        }};

        if (!useEvSel) {
          for (const auto& step : steps) {
            fillGenRecoCut(step.bin);
          }
          fillGenRecoCut(GenRecoCutBin::RctMFT);
          return true;
        }

        for (const auto& step : steps) {
          if (!step.enabled) {
            fillGenRecoCut(step.bin);
            continue;
          }

          if (!collision.selection_bit(step.bit)) {
            return false;
          }
          fillGenRecoCut(step.bin);
        }

        if (useRctMFT && !myChecker(collision)) {
          return false;
        }
        fillGenRecoCut(GenRecoCutBin::RctMFT);

        return true;
      };

      if (!countAndPassEvSelGenReco(collision)) {
        continue;
      }

      const auto z = collision.posZ();
      if ((z < cfgVzCut1) || (z > cfgVzCut2)) {
        continue;
      }
      fillGenRecoCut(GenRecoCutBin::VzWindow);

      auto perCollisionSampleCentral = midtracks.sliceBy(perColCentral, collision.globalIndex());
      if (perCollisionSampleCentral.size() <= 0) {
        continue;
      }
      fillGenRecoCut(GenRecoCutBin::InelGt0);

      //________________________________________________________________________________

      registry.fill(HIST("Purity/xReco"), collision.posX());
      registry.fill(HIST("Purity/xTrue"), mcCollision.posX());
      registry.fill(HIST("Purity/yReco"), collision.posY());

      registry.fill(HIST("Purity/yTrue"), mcCollision.posY());
      registry.fill(HIST("Purity/zReco"), collision.posZ());
      registry.fill(HIST("Purity/zTrue"), mcCollision.posZ());

      // --- Vertex position THnSparse: status axis (1=reco, 2=true) ---
      registry.fill(HIST("Purity/VtxXYZTruth"), mcCollision.posX(), mcCollision.posY(), mcCollision.posZ());
      registry.fill(HIST("Purity/VtxXYZReco"), collision.posX(), collision.posY(), collision.posZ());

      // --- Delta vertex THnSparse (reco - true) ---
      registry.fill(HIST("Purity/DeltaVtxXYZ"),
                    collision.posX() - mcCollision.posX(), // Reco - Truth
                    collision.posY() - mcCollision.posY(),
                    collision.posZ() - mcCollision.posZ());

      recoVtxByRecoId[collision.globalIndex()] = {collision.posX(), collision.posY(), collision.posZ()};
      recoVtxByMcId[collision.mcCollisionId()] = {mcCollision.posX(), mcCollision.posY(), mcCollision.posZ()};

      int64_t woOrpCount = 0;
      bool filledRight = false;
      bool filledWrong = false;
      int nMftSelectedAfterCuts = 0;
      std::unordered_set<int> uniqueBestRecoCols;
      // auto particle = templatedTrack.template mcParticle_as<aod::McParticles>();

      if (tracks.size() > 0) {
        bool countedPrimary = false;
        for (const auto& track : tracks) { // track loop starts
          // All compatible collisions assigned by track-to-collision-associator
          if (!(midtracks.size() > 0))
            continue;
          // std::cout <<" midtracks.size() " <<midtracks.size()<< std::endl; // "    track.globalIndex()  "<<track.globalIndex() This is track id

          float ndf = getTrackNdf(track);
          float chi2ndf = track.chi2() / ndf;
          float phi = track.phi();
          float dcaXyCut = track.bestDCAXY(); // UNCOMMENT THIS for tracks with reassociation 3d
          float dcaZCut = track.bestDCAZ();   // UNCOMMENT THIS for tracks with reassociation 3d
          float ptCut = track.pt();

          o2::math_utils::bringTo02Pi(phi);
          const float etaReco = track.eta();
          const float dcaXYReco = dcaXyCut; // track.bestDCAXY()
          const float dcaZReco = dcaZCut;   // track.bestDCAZ()
          const float dcaXReco = dcaXYReco * std::cos(phi);
          const float dcaYReco = dcaXYReco * std::sin(phi);

          const bool failTrackCuts =
            track.nClusters() < cfgnCluster ||
            etaReco <= cfgnEta1 ||
            etaReco >= cfgnEta2 ||
            chi2ndf >= cfgChi2NDFMax ||
            phi <= cfgPhiCut1 ||
            phi >= cfgPhiCut2 ||
            (usePhiCut &&
             ((phi <= PhiVetoLow) ||
              ((phi >= PhiVetoPiMin) && (phi <= PhiVetoPiMax)) ||
              (phi >= PhiVetoHigh))) ||
            (useDCAxyCut && dcaXyCut > maxDCAxy) ||
            (useDCAzCut && std::abs(dcaZCut) > maxDCAz) ||
            (usePtCut && ptCut > cfgnPt);

          if (failTrackCuts) {
            continue;
          }

          const bool hasMcLabel = track.has_mcParticle();
          const bool isFakeByLabel = hasMcLabel ? (track.mcMask() != 0) : false;
          const bool isTrueByLabel = hasMcLabel && !isFakeByLabel;
          const bool isPrimaryCharged = hasMcLabel && !isFakeByLabel && track.mcParticle().isPhysicalPrimary();
          const bool isSecondaryCharged = hasMcLabel && !isFakeByLabel && !track.mcParticle().isPhysicalPrimary();
          const auto mcColObj = track.mcParticle().mcCollision_as<McCollisionsWithExtra>();
          const auto mcPart = track.mcParticle();

          if (!passGenRecoTrackMode(track)) { //
            continue;
          }
          int bin = static_cast<int>(RightWrongBin::Neither);
          bool recoOfTrueExists = false;
          bool recoOfTrueInCompatible = false;

          const int bestColID = track.bestCollisionId(); // same as track.collisionId();
          if (isTrueByLabel) {
            auto itRecoToMc = recoToMc.find(bestColID);
            const int mcOfTrack = track.mcParticle().mcCollisionId();
            const auto compatibleIds = track.compatibleCollIds();
            auto itRecoList = mcToReco.find(mcOfTrack);

            // 1) First check whether the correct reco collision of the true MC collision
            //    is present in the compatible-collision list.
            if (!compatibleIds.empty() && itRecoList != mcToReco.end() && !itRecoList->second.empty()) {
              for (const auto& trueRecoId : itRecoList->second) {
                for (const auto& compatibleId : compatibleIds) {
                  if (compatibleId == trueRecoId) {
                    recoOfTrueInCompatible = true;
                    break;
                  }
                }
                if (recoOfTrueInCompatible) {
                  break;
                }
              }
            }

            // 2) Then check whether any reco collision for the true MC collision exists at all.
            if (itRecoList != mcToReco.end() && !itRecoList->second.empty()) {
              recoOfTrueExists = true;
            }

            // 3) Finally classify the actually chosen reco collision as right/wrong.
            if (bestColID >= 0 && itRecoToMc != recoToMc.end()) {
              const int mcFromReco = itRecoToMc->second;
              bin = (mcFromReco == mcOfTrack)
                      ? static_cast<int>(RightWrongBin::Right)
                      : static_cast<int>(RightWrongBin::Wrong);
            }
          }

          registry.fill(HIST("RightWrong"), bin);
          registry.fill(HIST("Purity/RecoOfTrueExists"),
                        recoOfTrueExists ? static_cast<int>(BoolBin::Yes)
                                         : static_cast<int>(BoolBin::No));
          registry.fill(HIST("Purity/RecoOfTrueInCompatible"),
                        recoOfTrueInCompatible ? static_cast<int>(BoolBin::Yes)
                                               : static_cast<int>(BoolBin::No));

          if (bin == static_cast<int>(RightWrongBin::Wrong)) {
            registry.fill(HIST("Purity/WrongVertexRecoExists"),
                          recoOfTrueExists ? static_cast<int>(WrongVertexRecoExistsBin::RecoOfTrueExists)
                                           : static_cast<int>(WrongVertexRecoExistsBin::RecoOfTrueMissing));
            registry.fill(HIST("Purity/RecoOfTrueExistsW"),
                          recoOfTrueExists ? static_cast<int>(BoolBin::Yes)
                                           : static_cast<int>(BoolBin::No));
            registry.fill(HIST("Purity/RecoOfTrueInCompatibleW"),
                          recoOfTrueInCompatible ? static_cast<int>(BoolBin::Yes)
                                                 : static_cast<int>(BoolBin::No));
          }

          if (bin == static_cast<int>(RightWrongBin::Right)) {
            registry.fill(HIST("Purity/RecoOfTrueExistsR"),
                          recoOfTrueExists ? static_cast<int>(BoolBin::Yes)
                                           : static_cast<int>(BoolBin::No));
            registry.fill(HIST("Purity/RecoOfTrueInCompatibleR"),
                          recoOfTrueInCompatible ? static_cast<int>(BoolBin::Yes)
                                                 : static_cast<int>(BoolBin::No));
          }

          registry.fill(HIST("Purity/RecoSparseAll"),
                        etaReco, dcaXYReco, dcaZReco, dcaXReco, dcaYReco);

          if (isPrimaryCharged) {
            registry.fill(HIST("Purity/RecoSparsePrimary"),
                          etaReco, dcaXYReco, dcaZReco, dcaXReco, dcaYReco);
          } else {
            registry.fill(HIST("Purity/RecoSparseSecondary"),
                          etaReco, dcaXYReco, dcaZReco, dcaXReco, dcaYReco);
          }

          registry.fill(HIST("RecoSparseAllBest"),
                        etaReco, dcaXYReco, dcaZReco, dcaXReco, dcaYReco);
          if (bin == static_cast<int>(RightWrongBin::Wrong)) {
            float vzBest = 999.;
            float vzTrue = 999.;
            auto itVzBest = recoVtxZ.find(bestColID);
            if (itVzBest != recoVtxZ.end()) {
              vzBest = itVzBest->second;
            }
            auto itVzTrue = recoToMcVZ.find(vzBest);
            if (itVzTrue != recoToMcVZ.end()) {
              vzTrue = itVzBest->second;
            }
            double_t vztrueParticle = mcColObj.posZ();
            double_t diff1 = vzBest - vztrueParticle;
            double_t diff2 = vzBest - vzTrue;
            registry.fill(HIST("deltaVZ_fromReco"), diff1);
            registry.fill(HIST("deltaVZ_fromTrue"), diff2);
            registry.fill(HIST("RecoSparseAllBestWrong"),
                          etaReco, dcaXYReco, dcaZReco, dcaXReco, dcaYReco);
          }

          // Truth DCA components w.r.t. MC collision vertex
          const auto dcaXtruth = mcPart.vx() - mcColObj.posX(); // here mcColObj -> mcPart.mcCollision()
          const auto dcaYtruth = mcPart.vy() - mcColObj.posY();
          const auto dcaZtruth = mcPart.vz() - mcColObj.posZ();
          const auto dcaXYtruth = std::sqrt(dcaXtruth * dcaXtruth +
                                            dcaYtruth * dcaYtruth);

          const float etaTruth = mcPart.eta();
          const bool isPrimaryTruth = mcPart.isPhysicalPrimary();

          // Base truth histograms (independent of right/wrong vertex)
          registry.fill(HIST("Tracks/dca/Truth/THnDCAxyBestGenTruthAll"),
                        etaTruth, dcaXYtruth, dcaZtruth, dcaXtruth, dcaYtruth);
          if (isPrimaryTruth) {
            registry.fill(HIST("Tracks/dca/Truth/THnDCAxyBestGenTruthPrim"),
                          etaTruth, dcaXYtruth, dcaZtruth, dcaXtruth, dcaYtruth);
          } else {
            registry.fill(HIST("Tracks/dca/Truth/THnDCAxyBestGenTruthSec"),
                          etaTruth, dcaXYtruth, dcaZtruth, dcaXtruth, dcaYtruth);
          }

          registry.fill(HIST("Purity/reco/woOrp/woOrpTracksEtaZvtx"), track.eta(), z);
          registry.fill(HIST("Purity/reco/woOrp/woOrpTracksPtZvtx"), track.pt(), z);
          if (perCollisionSampleCentral.size() > 0) {
            registry.fill(HIST("Purity/reco/woOrp/woOrpEtaZvtx_gt0"), track.eta(), z);
            registry.fill(HIST("Purity/reco/woOrp/woOrpPtZvtx_gt0"), track.pt(), z);
            registry.fill(HIST("Purity/reco/woOrp/woOrpTracksDCAxyZvtx_gt0"), dcaXyCut, z);
            registry.fill(HIST("Purity/reco/woOrp/woOrpTracksDCAzZvtx_gt0"), dcaZCut, z);
          }
          registry.fill(HIST("Purity/reco/woOrp/woOrpTracksPhiEta"), phi, track.eta());

          if (isFakeByLabel) {
            //  std::cout << " track.eta() " <<track.eta()<<std::endl;
            registry.fill(HIST("Purity/reco/woOrp_fake/woOrpTracksEtaZvtx"), track.eta(), z);
            registry.fill(HIST("Purity/reco/woOrp_fake/woOrpTracksPtZvtx"), track.pt(), z);
            registry.fill(HIST("Purity/reco/woOrp_fake/woOrpTracksPhiEta"), phi, track.eta());
            if (perCollisionSampleCentral.size() > 0) {
              registry.fill(HIST("Purity/reco/woOrp_fake/woOrpEtaZvtx_gt0"), track.eta(), z);
              registry.fill(HIST("Purity/reco/woOrp_fake/woOrpPtZvtx_gt0"), track.pt(), z);
            }
          }
          if (isTrueByLabel) {
            // Has MC particle
            //  std::cout << " track.eta() has mc particle " <<track.pt()<<std::endl;
            registry.fill(HIST("Purity/reco/woOrp_hasMC/woOrpTracksEtaZvtx"), track.eta(), z);
            registry.fill(HIST("Purity/reco/woOrp_hasMC/woOrpTracksPtZvtx"), track.pt(), z);
            registry.fill(HIST("Purity/reco/woOrp_hasMC/woOrpTracksPhiEta"), phi, track.eta());
            if (perCollisionSampleCentral.size() > 0) {
              registry.fill(HIST("Purity/reco/woOrp_hasMC/woOrpEtaZvtx_gt0"), track.eta(), z);
              registry.fill(HIST("Purity/reco/woOrp_hasMC/woOrpPtZvtx_gt0"), track.pt(), z);
            }
          }
          if (isSecondaryCharged) {
            registry.fill(HIST("Purity/reco/woOrp_secondary/woOrpTracksEtaZvtx"), track.eta(), z);
            registry.fill(HIST("Purity/reco/woOrp_secondary/woOrpTracksPtZvtx"), track.pt(), z);
            registry.fill(HIST("Purity/reco/woOrp_secondary/woOrpTracksPhiEta"), phi, track.eta());
            if (perCollisionSampleCentral.size() > 0) {
              registry.fill(HIST("Purity/reco/woOrp_secondary/woOrpEtaZvtx_gt0"), track.eta(), z);
              registry.fill(HIST("Purity/reco/woOrp_secondary/woOrpPtZvtx_gt0"), track.pt(), z);
            }
          }
          if (isPrimaryCharged) {
            registry.fill(HIST("Purity/reco/woOrp_primary/woOrpTracksEtaZvtx"), track.eta(), z);
            registry.fill(HIST("Purity/reco/woOrp_primary/woOrpTracksPtZvtx"), track.pt(), z);
            registry.fill(HIST("Purity/reco/woOrp_primary/woOrpTracksPhiEta"), phi, track.eta());
            if (perCollisionSampleCentral.size() > 0) {
              registry.fill(HIST("Purity/reco/woOrp_primary/woOrpEtaZvtx_gt0"), track.eta(), z);
              registry.fill(HIST("Purity/reco/woOrp_primary/woOrpPtZvtx_gt0"), track.pt(), z);
            }
          }

          ++woOrpCount;
          // Category 2: all reco tracks after selections (woOrp)

          // --- Primary vs Fake accounting ---
          const float xTrue = mcColObj.posX();
          const float yTrue = mcColObj.posY();
          const float zTrue = mcColObj.posZ();

          std::array<double, 3> dcaInfOrig{999., 999., 999.};
          std::array<double, 2> dcaChosen{999., 999.};          // (DCAxy, DCAz) to chosen reco vertex
          std::array<double, 2> dcaRight{999., 999.};           // (DCAxy, DCAz) to truth MC vertex (reference)
          std::array<double, 3> dcaChosenXYZ{999., 999., 999.}; // (DCAx, DCAy, DCAz) to chosen reco vertex

          const double bZ = o2::base::Propagator::Instance()->getNominalBz();

          // Build forward track parameters once per track, then use a fresh copy for each vertex propagation
          std::vector<double> v1; // empty -> null cov
          SMatrix55 tcovs(v1.begin(), v1.end());
          SMatrix5 tpars(track.x(), track.y(), track.phi(), track.tgl(), track.signed1Pt());
          o2::track::TrackParCovFwd trackPar0{track.z(), tpars, tcovs, track.chi2()};

          auto trackPar = trackPar0; // copy
          dcaInfOrig = {999., 999., 999.};
          const std::array<double, 3> vtxChosen{collision.posX(), collision.posY(), collision.posZ()};
          trackPar.propagateToDCAhelix(bZ, vtxChosen, dcaInfOrig);
          dcaChosenXYZ = dcaInfOrig; // store (DCAx, DCAy, DCAz)
          dcaChosen[0] = std::sqrt(dcaInfOrig[0] * dcaInfOrig[0] + dcaInfOrig[1] * dcaInfOrig[1]);
          dcaChosen[1] = dcaInfOrig[2];

          dcaInfOrig = {999., 999., 999.};
          const std::array<double, 3> vtxTruth{xTrue, yTrue, zTrue};
          trackPar.propagateToDCAhelix(bZ, vtxTruth, dcaInfOrig);
          dcaRight[0] = std::sqrt(dcaInfOrig[0] * dcaInfOrig[0] + dcaInfOrig[1] * dcaInfOrig[1]);
          dcaRight[1] = dcaInfOrig[2];

          registry.fill(HIST("Purity/DCAyVsDCAx_Right"), dcaChosenXYZ[2], dcaChosenXYZ[1]);
          // Fill only for WRONG-vertex tracks (the diagnostic you want)
          if (bin == static_cast<int>(RightWrongBin::Wrong)) {
            //  std::cout <<"dcaxy choosen  " << dcaXyCut<<" dcaxy calculated  "<<dcaChosen[0]<<" dcaxy right "<<dcaRight[0]<< "dcaz choosen  " << dcaZCut<<" dcaz calculated  "<<dcaChosen[1]<<" dcaz right "<<dcaRight[1]<< std::endl; // "    track.globalIndex()  "<<track.globalIndex() This is track id

            registry.fill(HIST("Purity/THnDCAChosenVsRight_Wrong"),
                          etaReco, dcaXyCut, dcaChosen[0], dcaRight[0], dcaZCut, dcaChosen[1], dcaRight[1]);
          }
          if (bin == static_cast<int>(RightWrongBin::Right)) {
            // std::cout <<"dcaxy choosen  " << dcaXyCut<<" dcaxy calculated  "<<dcaChosen[0]<<" dcaxy right "<<dcaRight[0]<< "dcaz choosen  " << dcaZCut<<" dcaz calculated  "<<dcaChosen[1]<<" dcaz right "<<dcaRight[1]<< std::endl; // "    track.globalIndex()  "<<track.globalIndex() This is track id

            registry.fill(HIST("Purity/THnDCAChosenVsRight_Right"),
                          etaReco, dcaXyCut, dcaChosen[0], dcaRight[0], dcaZCut, dcaChosen[1], dcaRight[1]);

            ++nMftSelectedAfterCuts;
          }
          // --- end DCA-to-vertex diagnostic ---

          if (bin == static_cast<int>(RightWrongBin::Right)) {
            // Right-vertex: all / primary / secondary
            registry.fill(HIST("Purity/RecoSparseRightAll"),
                          etaReco, dcaXYReco, dcaZReco, dcaXReco, dcaYReco);
            if (!filledRight) {
              registry.fill(HIST("Purity/RecoSparseRightAll_EventCount"), static_cast<int>(SingleCountBin::Count));

              filledRight = true;
            }
            if (isPrimaryCharged) {
              registry.fill(HIST("Purity/RecoSparseRightPrimary"),
                            etaReco, dcaXYReco, dcaZReco, dcaXReco, dcaYReco);
            } else {
              registry.fill(HIST("Purity/RecoSparseRightSecondary"),
                            etaReco, dcaXYReco, dcaZReco, dcaXReco, dcaYReco);
            }
          } else if (bin == static_cast<int>(RightWrongBin::Wrong)) {
            // Wrong-vertex: all / primary / secondary
            registry.fill(HIST("Purity/RecoSparseWrongAll"),
                          etaReco, dcaXYReco, dcaZReco, dcaXReco, dcaYReco);
            if (!filledWrong) {
              registry.fill(HIST("Purity/RecoSparseWrongAll_EventCount"), static_cast<int>(SingleCountBin::Count));
              filledWrong = true;
            }
            if (isPrimaryCharged) {
              registry.fill(HIST("Purity/RecoSparseWrongPrimary"),
                            etaReco, dcaXYReco, dcaZReco, dcaXReco, dcaYReco);
            } else {
              registry.fill(HIST("Purity/RecoSparseWrongSecondary"),
                            etaReco, dcaXYReco, dcaZReco, dcaXReco, dcaYReco);
            }
          }
          // Reconstructed vertex position for bestbestColID, if available
          auto itVtxX = recoVtxX.find(bestColID);

          if (itVtxX != recoVtxX.end()) {

            const float xReco = itVtxX->second;
            const float yReco = recoVtxY[bestColID];
            const float zReco = recoVtxZ[bestColID];
            const bool recoVzIn = (zReco >= cfgVzCut1) && (zReco <= cfgVzCut2);
            const bool trueVzIn = (zTrue >= cfgVzCut1) && (zTrue <= cfgVzCut2);

            if (!(recoVzIn && trueVzIn)) {
              continue; // skip filling Delta* histos for this track
            }

            const float deltaXvtx = xReco - xTrue;
            const float deltaYvtx = yReco - yTrue;
            const float deltaZvtx = zReco - zTrue;

            // We are interested in how far the WRONG collisions are from the true one
            if (bin == static_cast<int>(RightWrongBin::Wrong)) {
              registry.fill(HIST("Purity/DeltaXWrong"), deltaXvtx);
              registry.fill(HIST("Purity/DeltaYWrong"), deltaYvtx);
              registry.fill(HIST("Purity/DeltaZWrong"), deltaZvtx);
            }
            if (bin == static_cast<int>(RightWrongBin::Right)) {
              registry.fill(HIST("Purity/DeltaXRight"), deltaXvtx);
              registry.fill(HIST("Purity/DeltaYRight"), deltaYvtx);
              registry.fill(HIST("Purity/DeltaZRight"), deltaZvtx);
            }
          }

          // --- Delta-DCA components: truth minus reco ---
          const float deltaDCAxy = dcaXYtruth - dcaXYReco;
          const float deltaDCAz = dcaZtruth - dcaZReco;
          const float deltaDCAx = dcaXtruth - dcaXReco;
          const float deltaDCAy = dcaYtruth - dcaYReco;
          if (bin == static_cast<int>(RightWrongBin::Right)) {
            // Right-vertex: all / primary / secondary
            registry.fill(HIST("Tracks/dca/Truth/THnDeltaDCARightAll"),
                          deltaDCAxy, deltaDCAz, deltaDCAx, deltaDCAy);
            if (isPrimaryCharged) {
              registry.fill(HIST("Tracks/dca/Truth/THnDeltaDCARightPrim"),
                            deltaDCAxy, deltaDCAz, deltaDCAx, deltaDCAy);
            } else {
              registry.fill(HIST("Tracks/dca/Truth/THnDeltaDCARightSec"),
                            deltaDCAxy, deltaDCAz, deltaDCAx, deltaDCAy);
            }
          } else {
            // Wrong-vertex: all / primary / secondary
            registry.fill(HIST("Tracks/dca/Truth/THnDeltaDCAWrongAll"),
                          deltaDCAxy, deltaDCAz, deltaDCAx, deltaDCAy);
            if (isPrimaryCharged) {
              registry.fill(HIST("Tracks/dca/Truth/THnDeltaDCAWrongPrim"),
                            deltaDCAxy, deltaDCAz, deltaDCAx, deltaDCAy);
            } else {
              registry.fill(HIST("Tracks/dca/Truth/THnDeltaDCAWrongSec"),
                            deltaDCAxy, deltaDCAz, deltaDCAx, deltaDCAy);
            }
          }

          // auto mcColObj = track.mcParticle().mcCollision_as<McCollisionsWithExtra>();
          // True primary match (purity numerator)
          registry.fill(HIST("Purity/mc/PrimaryAll"), static_cast<int>(SingleCountBin::Count));
          registry.fill(HIST("Purity/mc/PrimaryAllEta"), mcPart.eta());
          registry.fill(HIST("Purity/mc/PrimaryTracksEtaZvtx"), mcPart.eta(), mcCollision.posZ());
          if (perCollisionSampleCentral.size() > 0) {
            registry.fill(HIST("Purity/mc/PrimaryTracksEtaZvtx_gt0"), mcPart.eta(), mcCollision.posZ());
            registry.fill(HIST("Purity/mc/PrimaryTracksPtZvtx_gt0"), mcPart.pt(), mcCollision.posZ());
            registry.fill(HIST("Purity/mc/PrimaryTracksDCAxyZvtx_gt0"), dcaXyCut, mcCollision.posZ());
            registry.fill(HIST("Purity/mc/PrimaryTracksDCAzZvtx_gt0"), dcaZCut, mcCollision.posZ());
          }
          registry.fill(HIST("Purity/mc/PrimaryTracksPhiEta"), mcPart.phi(), mcPart.eta());
          registry.fill(HIST("Purity/SelectedAfterDCAxy/PrimaryAll"), static_cast<int>(SingleCountBin::Count));
          registry.fill(HIST("Purity/SelectedAfterDCAxy/PrimaryAllEta"), mcPart.eta());
          countedPrimary = true;

          // --- Purity profiles (mean of indicator gives purity) ---
          registry.fill(HIST("Purity/PurityOverall"),
                        static_cast<int>(SingleCountBin::Count),
                        countedPrimary ? static_cast<int>(BoolBin::Yes)
                                       : static_cast<int>(BoolBin::No));

          registry.fill(HIST("Purity/PurityVsEta"), track.eta(),
                        countedPrimary ? static_cast<int>(BoolBin::Yes)
                                       : static_cast<int>(BoolBin::No));
        } // track loop
        registry.fill(HIST("Purity/HashTableRowCounts"),
                      static_cast<int>(HashTableRowCountsBin::UniqueBestRecoCols), uniqueBestRecoCols.size());
      }
      registry.fill(HIST("Purity/reco/woOrp/nTrk"), woOrpCount);
      registry.fill(HIST("Purity/reco/PNchMFT_afterCuts"), nMftSelectedAfterCuts);
    } // collision
  }

  PROCESS_SWITCH(PseudorapidityDensityMFT, processGenReco,
                 "Process particle-level info of pt", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PseudorapidityDensityMFT>(cfgc)};
}
