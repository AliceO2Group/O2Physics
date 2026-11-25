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

/// \file emcTmMonitor.cxx
/// \brief Simple monitoring task for EMCal clusters
/// \author Marvin Hemmer <marvin.hemmer@cern.ch>
/// \since 24.02.2023
///
/// This task is meant to be used for monitoring the matching between global tracks and EMCal clusters
/// properties, such as:
/// - cluster energy over track momentum
/// - difference in eta
/// - difference in phi
/// The event selection (eventCuts) applies the standard cuts of the event selection task: sel8, z vertex, the evsel
/// flags (pile-up, ITS/TF frame borders, vertex quality, ITS layers), occupancy, RCT flags and a centrality or
/// multiplicity range with a configurable estimator. Optionally EMCal in the readout (Run 2: kINT7, Run 3: kTVXinEMC).
/// For pilot beam data, specific BC IDs can be vetoed (vetoBCID) or selected (selectBCID).
/// The MC process (processMCTruth) follows the leading MC particle of each cluster through the matching chain:
/// is its track reconstructed, extrapolated to the EMCal, a matching candidate (global track or V0 leg),
/// and matched to the cluster as primary or secondary track.

#include "PWGEM/PhotonMeson/DataModel/gammaTables.h" // for the V0 legs used in the secondary matching
#include "PWGJE/DataModel/EMCALClusters.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/LHCConstants.h>
#include <CommonConstants/MathConstants.h>
#include <EMCALBase/Geometry.h>
#include <EMCALBase/GeometryBase.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Concepts.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TAxis.h>
#include <TH2.h>
#include <THn.h>
#include <TMCProcess.h>
#include <TPDGCode.h>
#include <TParticlePDG.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace o2::framework;
using namespace o2::framework::expressions;
using CollisionEvSelIt = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels, o2::aod::CentFT0Cs, o2::aod::CentFT0Ms, o2::aod::CentFT0As, o2::aod::FT0Mults, o2::aod::PVMults>::iterator;
using BcEvSelIt = o2::soa::Join<o2::aod::BCs, o2::aod::BcSels>::iterator;
using SelectedClusters = o2::soa::Filtered<o2::aod::EMCALClusters>;
using TracksPID = o2::soa::Join<o2::aod::pidTPCFullEl, o2::aod::pidTPCFullPi, o2::aod::FullTracks, o2::aod::TrackSelection>;
using SelectedMcClusters = o2::soa::Filtered<o2::soa::Join<o2::aod::EMCALClusters, o2::aod::EMCALMCClusters>>;
using TracksMC = o2::soa::Join<o2::aod::FullTracks, o2::aod::TrackSelection, o2::aod::McTrackLabels>;

/// \brief Truth class of the leading MC particle of a cluster
enum class TmTruthClass : uint8_t {
  kConversionElectron = 0, // e+- with a photon as mother, produced in the transport (photon conversion)
  kPrimaryElectron,        // e+- produced by the generator (Dalitz, heavy flavour, ...)
  kOtherElectron,          // any other e+- (delta rays, Compton, weak decays in the transport, ...)
  kOtherCharged,           // charged particle that is not an e+- (hadrons, muons)
  kPhoton,
  kNeutral, // neutral particle that is not a photon
  kOther,   // no MC particle or unknown PDG code
  kNClasses
};

/// \brief What happened to the track of the leading MC particle of a cluster in the track matching
enum class TmTruthStatus : uint8_t {
  kNoRecoTrack = 0,     // no reconstructed track of this particle in this collision
  kNotExtrapolated,     // track(s) exist, but none was extrapolated to the EMCal (trackEta/PhiEmcal = -999)
  kNoMatchingCandidate, // extrapolated, but neither a global track nor a V0 leg: input to no matching
  kCandidateNotMatched, // matching candidate, but not matched to this cluster (dR > maxMatchingDistance or not among the closest tracks)
  kMatchedPrimary,      // matched to this cluster as primary (global) track
  kMatchedSecondary,    // matched to this cluster as secondary track (V0 leg)
  kNStatus
};

/// \brief Type of the reconstructed track
enum class TmTrackType : uint8_t {
  kNoTrack = 0,
  kGlobal, // passes isGlobalTrack()
  kItsTpc, // ITS and TPC, but not a global track
  kTpcOnly,
  kItsOnly,
  kOtherType,
  kNTypes
};

/// \brief Steps of the event selection, bins of the cut flow histogram (a disabled cut passes every event)
enum class EventCut : uint8_t {
  kAll = 0,
  kBC,
  kSel8,
  kEmcalReadout,
  kZVertex,
  kNoTimeFrameBorder,
  kNoITSROFrameBorder,
  kNoSameBunchPileup,
  kGoodZvtxFT0vsPV,
  kVertexITSTPC,
  kVertexTOFmatched,
  kNoCollInTimeRangeStandard,
  kNoCollInTimeRangeNarrow,
  kNoCollInRofStandard,
  kNoHighMultCollInPrevRof,
  kGoodITSLayer3,
  kGoodITSLayer0123,
  kGoodITSLayersAll,
  kTrackOccupancy,
  kFT0COccupancy,
  kRCT,
  kCentMult,
  kSelected,
  kNCuts
};

/// \brief Estimator of the centrality / multiplicity cut
enum class CentMultEstimator : int {
  kOff = 0,
  kCentFT0C,
  kCentFT0M,
  kCentFT0A,
  kMultFT0C,
  kMultFT0A,
  kMultFT0M,
  kMultNTracksPV,
  kNEstimators
};

struct EmcTmMonitor {
  HistogramRegistry mHistManager{"TrackMatchingMonitorHistograms", {}, OutputObjHandlingPolicy::AnalysisObject};
  o2::emcal::Geometry* mGeometry = nullptr;
  Service<o2::framework::O2DatabasePDG> mPdgDatabase{};

  Preslice<o2::aod::EMCALClusterCells> perCluster = o2::aod::emcalclustercell::emcalclusterId;
  Preslice<o2::aod::EMCALAmbiguousClusterCells> perClusterAmb = o2::aod::emcalclustercell::emcalambiguousclusterId;
  Preslice<o2::aod::EMCALMatchedTracks> perClusterMatchedTracks = o2::aod::emcalclustercell::emcalclusterId;
  Preslice<o2::aod::EMCMatchSecs> perClusterMatchedSecondaries = o2::aod::emcalclustercell::emcalclusterId;
  // configurable parameters
  Configurable<bool> doEventSel{"doEventSel", false, "deprecated, use eventCuts.requireEmcalReadout: demand kINT7 (Run 2) / kTVXinEMC (Run 3)"};
  Configurable<std::string> vetoBCID{"vetoBCID", "", "BC ID(s) to be excluded, this should be used as an alternative to the event selection"};
  Configurable<std::string> selectBCID{"selectBCID", "all", "BC ID(s) to be included, this should be used as an alternative to the event selection"};
  Configurable<double> vertexCut{"vertexCut", -1, "deprecated, use eventCuts.zVertexMax: if > 0, overrides it"};
  Configurable<int> clusterDefinition{"clusterDefinition", 10, "cluster definition to be selected, e.g. 10=kV3Default"};
  ConfigurableAxis clusterTimeBinning{"clusterTimeBinning", {1500, -600, 900}, ""};
  Configurable<bool> usePionRejection{"usePionRejection", false, "demand pion rection for electron signal with TPC PID"};
  Configurable<std::vector<float>> tpcNsigmaElectron{"tpcNsigmaElectron", {-1., +3.}, "TPC PID NSigma range for electron signal (first <= NSigma <= second)"};
  Configurable<std::vector<float>> tpcNsigmaBack{"tpcNsigmaBack", {-10., -4.}, "TPC PID NSigma range for electron background (first <= NSigma <= second)"};
  Configurable<std::vector<float>> tpcNsigmaPion{"tpcNsigmaPion", {-3., +3.}, "TPC PID NSigma range for pions for rejection if usePionRejection is enabled (first <= NSigma <= second)"};
  Configurable<float> minTime{"minTime", -25., "Minimum cluster time for time cut"};
  Configurable<float> maxTime{"maxTime", +20., "Maximum cluster time for time cut"};
  Configurable<float> minM02{"minM02", 0.1, "Minimum M02 for M02 cut"};
  Configurable<float> maxM02{"maxM02", 0.9, "Maximum M02 for M02 cut"};
  Configurable<float> maxM02HighPt{"maxM02HighPt", 0.6, "Maximum M02 for M02 cut for high pT"};
  Configurable<float> m02highPt{"m02highPt", 15., "pT threshold for maxM02HighPt cut. Set to negative value to disable it!"};
  Configurable<float> minDEta{"minDEta", 0.01, "Minimum dEta between track and cluster"};
  Configurable<float> minDPhi{"minDPhi", 0.01, "Minimum dPhi between track and cluster"};
  Configurable<std::vector<float>> eOverPRange{"eOverPRange", {0.9, 1.2}, "E/p range where one would search for electrons (first <= E/p <= second)"};
  Configurable<bool> hasTRD{"hasTRD", false, "demand track to have a hit in TRD."};
  Configurable<float> emcalRadius{"emcalRadius", 440.f, "MC truth: radius (cm) used to move the cluster eta from the vertex to the origin frame of the track extrapolation"};

  struct : ConfigurableGroup {
    std::string prefix = "eventCuts";
    Configurable<bool> requireSel8{"requireSel8", true, "require sel8"};
    Configurable<bool> requireEmcalReadout{"requireEmcalReadout", false, "require EMCal in the readout: kINT7 (Run 2) / kTVXinEMC (Run 3)"};
    Configurable<float> zVertexMax{"zVertexMax", 10.f, "maximum |z| of the primary vertex (cm), <= 0: off"};
    Configurable<bool> requireNoTimeFrameBorder{"requireNoTimeFrameBorder", false, "require kNoTimeFrameBorder"};
    Configurable<bool> requireNoITSROFrameBorder{"requireNoITSROFrameBorder", false, "require kNoITSROFrameBorder"};
    Configurable<bool> requireNoSameBunchPileup{"requireNoSameBunchPileup", false, "require kNoSameBunchPileup"};
    Configurable<bool> requireGoodZvtxFT0vsPV{"requireGoodZvtxFT0vsPV", false, "require kIsGoodZvtxFT0vsPV"};
    Configurable<bool> requireVertexITSTPC{"requireVertexITSTPC", false, "require kIsVertexITSTPC"};
    Configurable<bool> requireVertexTOFmatched{"requireVertexTOFmatched", false, "require kIsVertexTOFmatched"};
    Configurable<bool> requireNoCollInTimeRangeStandard{"requireNoCollInTimeRangeStandard", false, "require kNoCollInTimeRangeStandard"};
    Configurable<bool> requireNoCollInTimeRangeNarrow{"requireNoCollInTimeRangeNarrow", false, "require kNoCollInTimeRangeNarrow"};
    Configurable<bool> requireNoCollInRofStandard{"requireNoCollInRofStandard", false, "require kNoCollInRofStandard"};
    Configurable<bool> requireNoHighMultCollInPrevRof{"requireNoHighMultCollInPrevRof", false, "require kNoHighMultCollInPrevRof"};
    Configurable<bool> requireGoodITSLayer3{"requireGoodITSLayer3", false, "require kIsGoodITSLayer3"};
    Configurable<bool> requireGoodITSLayer0123{"requireGoodITSLayer0123", false, "require kIsGoodITSLayer0123"};
    Configurable<bool> requireGoodITSLayersAll{"requireGoodITSLayersAll", false, "require kIsGoodITSLayersAll"};
    Configurable<float> trackOccupancyMin{"trackOccupancyMin", -1.f, "minimum track occupancy (trackOccupancyInTimeRange), < 0: off"};
    Configurable<float> trackOccupancyMax{"trackOccupancyMax", -1.f, "maximum track occupancy (trackOccupancyInTimeRange), < 0: off"};
    Configurable<float> ft0cOccupancyMin{"ft0cOccupancyMin", -1.f, "minimum FT0C occupancy (ft0cOccupancyInTimeRange), < 0: off"};
    Configurable<float> ft0cOccupancyMax{"ft0cOccupancyMax", -1.f, "maximum FT0C occupancy (ft0cOccupancyInTimeRange), < 0: off"};
    Configurable<bool> applyRCT{"applyRCT", false, "apply the RCT flags of rctLabel"};
    Configurable<std::string> rctLabel{"rctLabel", "CBT_calo", "RCT selection label: CBT, CBT_hadronPID, CBT_electronPID, CBT_calo, CBT_muon, CBT_muon_glo"};
    Configurable<bool> rctCheckZDC{"rctCheckZDC", false, "RCT: also require good ZDC quality (Pb-Pb)"};
    Configurable<bool> rctTreatLimitedAcceptanceAsBad{"rctTreatLimitedAcceptanceAsBad", false, "RCT: reject runs flagged LimitedAcceptanceMCReproducible"};
    Configurable<int> centMultEstimator{"centMultEstimator", 0, "centrality / multiplicity cut: 0 = off, 1 = cent FT0C, 2 = cent FT0M, 3 = cent FT0A, 4 = mult FT0C, 5 = mult FT0A, 6 = mult FT0M, 7 = mult NTracksPV"};
    Configurable<float> centMultMin{"centMultMin", 0.f, "minimum centrality (%) / multiplicity"};
    Configurable<float> centMultMax{"centMultMax", 100.f, "maximum centrality (%) / multiplicity"};
  } eventCuts;

  std::vector<int> mVetoBCIDs;
  std::vector<int> mSelectBCIDs;
  o2::aod::rctsel::RCTFlagsChecker mRctChecker;

  /// \brief Create output histograms and initialize geometry
  void init(InitContext const&)
  {
    // load geometry just in case we need it
    mGeometry = o2::emcal::Geometry::GetInstanceFromRunNumber(300000);

    // create common axes
    LOG(info) << "Creating histograms";
    const AxisSpec bcAxis{3501, -0.5, 3500.5};
    const AxisSpec energyAxis{makeEnergyBinningAliPhysics(), "E_{clus} (GeV)"};
    const AxisSpec amplitudeAxisLarge{1000, 0., 100., "amplitudeLarge", "Amplitude (GeV)"};
    const AxisSpec dEtaAxis{100, -1.f * minDEta, minDEta, "d#it{#eta}"};
    const AxisSpec dPhiAxis{100, -1.f * minDPhi, minDPhi, "d#it{#varphi} (rad)"};
    const AxisSpec dRAxis{150, 0.0, 0.015, "d#it{R}"};
    const AxisSpec eoverpAxis{500, 0, 10, "#it{E}_{cluster}/#it{p}_{track}"};
    const AxisSpec nSigmaAxis{400, -10., +30., "N#sigma"};
    const AxisSpec trackptAxis{makePtBinning(), "#it{p}_{T,track}"};
    const AxisSpec trackpAxis{200, 0, 100, "#it{p}_{track}"};
    const AxisSpec clusterptAxis{makePtBinning(), "#it{p}_{T}"};
    const AxisSpec etaAxis{160, -0.8, 0.8, "#eta"};
    const AxisSpec phiAxis{72, 0, 2 * 3.14159, "#varphi (rad)"};
    const AxisSpec smAxis{20, -0.5, 19.5, "SM"};
    AxisSpec timeAxis{clusterTimeBinning, "t_{cl} (ns)"};

    const int maxMatched = 20; // maximum number of matched tracks, hardcoded in emcalCorrectionTask.cxx!
    const AxisSpec nmatchedtrack{maxMatched, -0.5, maxMatched + 0.5};

    // create histograms
    // event properties
    mHistManager.add("eventsAll", "Number of events", HistType::kTH1F, {{1, 0.5, 1.5}});
    mHistManager.add("eventsSelected", "Number of events", HistType::kTH1F, {{1, 0.5, 1.5}});
    mHistManager.add("eventVertexZAll", "z-vertex of event (all events)", HistType::kTH1F, {{200, -20, 20}});
    mHistManager.add("eventVertexZSelected", "z-vertex of event (selected events)", HistType::kTH1F, {{200, -20, 20}});
    auto hCutFlow = mHistManager.add<TH1>("Event/hCutFlow", "events after each step of the event selection (a disabled cut passes every event);;events", HistType::kTH1D, {{static_cast<int>(EventCut::kNCuts), -0.5, static_cast<int>(EventCut::kNCuts) - 0.5}});
    static constexpr std::array<std::string_view, static_cast<int>(EventCut::kNCuts)> CutLabels = {"all", "BC ID", "sel8", "EMCal readout", "z vertex", "NoTimeFrameBorder", "NoITSROFrameBorder", "NoSameBunchPileup", "GoodZvtxFT0vsPV", "VertexITSTPC", "VertexTOFmatched", "NoCollInTimeRangeStandard", "NoCollInTimeRangeNarrow", "NoCollInRofStandard", "NoHighMultCollInPrevRof", "GoodITSLayer3", "GoodITSLayer0123", "GoodITSLayersAll", "track occupancy", "FT0C occupancy", "RCT", "cent/mult", "selected"};
    for (std::size_t i = 0; i < CutLabels.size(); ++i) {
      hCutFlow->GetXaxis()->SetBinLabel(static_cast<int>(i) + 1, CutLabels[i].data());
    }
    const AxisSpec trackOccupancyAxis{200, 0., 20000., "track occupancy"};
    const AxisSpec ft0cOccupancyAxis{200, 0., 200000., "FT0C occupancy"};
    const bool isCentrality = eventCuts.centMultEstimator.value <= static_cast<int>(CentMultEstimator::kCentFT0A);
    const AxisSpec centMultAxis{isCentrality ? AxisSpec{110, -5., 105., "centrality (%)"} : AxisSpec{500, 0., 50000., "multiplicity"}};
    mHistManager.add("Event/hTrackOccupancyBefore", "track occupancy (before the event selection);;events", HistType::kTH1D, {trackOccupancyAxis});
    mHistManager.add("Event/hTrackOccupancyAfter", "track occupancy (after the event selection);;events", HistType::kTH1D, {trackOccupancyAxis});
    mHistManager.add("Event/hFT0COccupancyBefore", "FT0C occupancy (before the event selection);;events", HistType::kTH1D, {ft0cOccupancyAxis});
    mHistManager.add("Event/hFT0COccupancyAfter", "FT0C occupancy (after the event selection);;events", HistType::kTH1D, {ft0cOccupancyAxis});
    if (eventCuts.centMultEstimator.value != static_cast<int>(CentMultEstimator::kOff)) {
      mHistManager.add("Event/hCentMultBefore", "centrality / multiplicity of the cut (before the event selection);;events", HistType::kTH1D, {centMultAxis});
      mHistManager.add("Event/hCentMultAfter", "centrality / multiplicity of the cut (after the event selection);;events", HistType::kTH1D, {centMultAxis});
    }
    mHistManager.add("Event/hBCId", "BC ID of the selected events;BC ID;events", HistType::kTH1D, {{o2::constants::lhc::LHCMaxBunches, -0.5, o2::constants::lhc::LHCMaxBunches - 0.5}});

    if (eventCuts.centMultEstimator.value < static_cast<int>(CentMultEstimator::kOff) || eventCuts.centMultEstimator.value >= static_cast<int>(CentMultEstimator::kNEstimators)) {
      LOG(fatal) << "eventCuts.centMultEstimator " << eventCuts.centMultEstimator.value << " unknown, see the configurable description";
    }
    if (eventCuts.applyRCT.value) {
      mRctChecker.init(eventCuts.rctLabel.value, eventCuts.rctCheckZDC.value, eventCuts.rctTreatLimitedAcceptanceAsBad.value);
    }
    if (doEventSel.value) {
      LOG(warning) << "doEventSel is deprecated, use eventCuts.requireEmcalReadout. Requiring EMCal in the readout.";
    }
    if (vertexCut.value > 0.) {
      LOG(warning) << "vertexCut is deprecated, use eventCuts.zVertexMax. Using |z| < " << vertexCut.value << " cm instead of eventCuts.zVertexMax.";
    }

    // cluster properties (matched clusters)
    mHistManager.add("TrackEtaPhi", "#eta vs #varphi of all selected tracks", HistType::kTH2F, {etaAxis, phiAxis});                                                                             // eta vs phi of all selected tracks
    mHistManager.add("TrackEtaPhi_Neg", "#eta vs #varphi of all selected negative tracks", HistType::kTH2F, {etaAxis, phiAxis});                                                                // eta vs phi of all selected negative tracks
    mHistManager.add("TrackEtaPhi_Pos", "#eta vs #varphi of all selected positive tracks", HistType::kTH2F, {etaAxis, phiAxis});                                                                // eta vs phi of all selected positive tracks
    mHistManager.add("clusterEMatched", "Energy of cluster (with match)", HistType::kTH1F, {energyAxis});                                                                                       // energy of matched clusters
    mHistManager.add("MatchedTrackEtaPhi_BeforeCut", "#eta vs #varphi of all selected matched tracks before d#eta and #dvarphi cut", HistType::kTH2F, {etaAxis, phiAxis});                      // eta vs phi of all selected matched tracks before dEta and dPhi cut
    mHistManager.add("MatchedTrackEtaPhi_Neg_BeforeCut", "#eta vs #varphi of all selected negative matched tracks before d#eta and #dvarphi cut", HistType::kTH2F, {etaAxis, phiAxis});         // eta vs phi of all selected negative matched tracks before dEta and dPhi cut
    mHistManager.add("MatchedTrackEtaPhi_Pos_BeforeCut", "#eta vs #varphi of all selected positive matched tracks before d#eta and #dvarphi cut", HistType::kTH2F, {etaAxis, phiAxis});         // eta vs phi of all selected positive matched tracks before dEta and dPhi cut
    mHistManager.add("MatchedTrackEtaPhi", "#eta vs #varphi of all selected matched tracks", HistType::kTH2F, {etaAxis, phiAxis});                                                              // eta vs phi of all selected matched tracks
    mHistManager.add("MatchedTrackEtaPhi_Neg", "#eta vs #varphi of all selected negative matched tracks", HistType::kTH2F, {etaAxis, phiAxis});                                                 // eta vs phi of all selected negative matched tracks
    mHistManager.add("MatchedTrackEtaPhi_Pos", "#eta vs #varphi of all selected positive matched tracks", HistType::kTH2F, {etaAxis, phiAxis});                                                 // eta vs phi of all selected positive matched tracks
    mHistManager.add("clusterTM_dEtadPhi", "cluster trackmatching dEta/dPhi", HistType::kTH3F, {dEtaAxis, dPhiAxis, smAxis});                                                                   // dEta dPhi per SM
    mHistManager.add("clusterTM_dEtadPhi_ASide", "cluster trackmatching in A-Side dEta/dPhi", HistType::kTH3F, {dEtaAxis, dPhiAxis, smAxis});                                                   // dEta dPhi per SM in A-Aside
    mHistManager.add("clusterTM_dEtadPhi_CSide", "cluster trackmatching in C-Side tracks dEta/dPhi", HistType::kTH3F, {dEtaAxis, dPhiAxis, smAxis});                                            // dEta dPhi per SM in C-Side
    mHistManager.add("clusterTM_PosdEtadPhi", "cluster trackmatching positive tracks dEta/dPhi", HistType::kTH3F, {dEtaAxis, dPhiAxis, smAxis});                                                // dEta dPhi per SM positive track
    mHistManager.add("clusterTM_NegdEtadPhi", "cluster trackmatching negative tracks dEta/dPhi", HistType::kTH3F, {dEtaAxis, dPhiAxis, smAxis});                                                // dEta dPhi per SM negative track
    mHistManager.add("clusterTM_PosdEtadPhi_Pl0_75", "cluster trackmatching positive tracks, p < 0.75 dEta/dPhi", HistType::kTH3F, {dEtaAxis, dPhiAxis, smAxis});                               // dEta dPhi per SM positive track with p < 0.75 GeV/c
    mHistManager.add("clusterTM_NegdEtadPhi_Pl0_75", "cluster trackmatching negative tracks, p < 0.75 dEta/dPhi", HistType::kTH3F, {dEtaAxis, dPhiAxis, smAxis});                               // dEta dPhi per SM negative track with p < 0.75 GeV/c
    mHistManager.add("clusterTM_PosdEtadPhi_0_75leqPl1_25", "cluster trackmatching positive tracks, 0.75 <= p < 1.25 dEta/dPhi", HistType::kTH3F, {dEtaAxis, dPhiAxis, smAxis});                // dEta dPhi per SM positive track with 0.75 <= p < 1.25 GeV/c
    mHistManager.add("clusterTM_NegdEtadPhi_0_75leqPl1_25", "cluster trackmatching negative tracks, 0.75 <= p < 1.25 dEta/dPhi", HistType::kTH3F, {dEtaAxis, dPhiAxis, smAxis});                // dEta dPhi per SM negative track with 0.75 <= p < 1.25 GeV/c
    mHistManager.add("clusterTM_PosdEtadPhi_Pgeq1_25", "cluster trackmatching positive tracks, p >= 1.25 dEta/dPhi", HistType::kTH3F, {dEtaAxis, dPhiAxis, smAxis});                            // dEta dPhi per SM positive track with p >= 1.25 GeV/c
    mHistManager.add("clusterTM_NegdEtadPhi_Pgeq1_25", "cluster trackmatching negative tracks, p >= 1.25 dEta/dPhi", HistType::kTH3F, {dEtaAxis, dPhiAxis, smAxis});                            // dEta dPhi per SM negative track with p >= 1.25 GeV/c
    mHistManager.add("clusterTM_dEtaPt", "cluster trackmatching dEta/#it{p}_{T};d#it{#eta};#it{p}_{T} (GeV/#it{c})", HistType::kTH3F, {dEtaAxis, clusterptAxis, smAxis});                       // dEta vs pT per SM
    mHistManager.add("clusterTM_PosdPhiPt", "cluster trackmatching positive tracks dPhi/#it{p}_{T}", HistType::kTH3F, {dPhiAxis, clusterptAxis, smAxis});                                       // dPhi vs pT per SM positive track
    mHistManager.add("clusterTM_NegdPhiPt", "cluster trackmatching negative tracks dPh/#it{p}_{T}", HistType::kTH3F, {dPhiAxis, clusterptAxis, smAxis});                                        // dPhi vs pT per SM negative track
    mHistManager.add("clusterTM_dEtaTN", "cluster trackmatching dEta/TN;d#it{#eta};#it{N}_{matched tracks}", HistType::kTH2F, {dEtaAxis, nmatchedtrack});                                       // dEta compared to the Nth closest track
    mHistManager.add("clusterTM_dPhiTN", "cluster trackmatching dPhi/TN;d#it{#varphi} (rad);#it{N}_{matched tracks}", HistType::kTH2F, {dPhiAxis, nmatchedtrack});                              // dPhi compared to the Nth closest track
    mHistManager.add("clusterTM_dRTN", "cluster trackmatching dR/TN;d#it{R};#it{N}_{matched tracks}", HistType::kTH2F, {dRAxis, nmatchedtrack});                                                // dR compared to the Nth closest track
    mHistManager.add("clusterTM_NTrack", "cluster trackmatching NMatchedTracks", HistType::kTH1I, {nmatchedtrack});                                                                             // how many tracks are matched
    mHistManager.add("clusterTM_EoverP_E", "E/p ", HistType::kTH3F, {eoverpAxis, energyAxis, nmatchedtrack});                                                                                   // E/p vs p for the Nth closest track
    mHistManager.add("clusterTM_EoverP_Pt", "E/p vs track pT", HistType::kTH3F, {eoverpAxis, trackptAxis, nmatchedtrack});                                                                      // E/p vs track pT for the Nth closest track
    mHistManager.add("clusterTM_EvsP", "cluster E/track p", HistType::kTH3F, {energyAxis, trackpAxis, nmatchedtrack});                                                                          // E vs p for the Nth closest track
    mHistManager.add("clusterTM_EoverP_ep", "cluster E/electron p", HistType::kTH3F, {eoverpAxis, trackptAxis, nmatchedtrack});                                                                 // E over p vs track pT for the Nth closest electron/positron track
    mHistManager.add("clusterTM_EoverP_e", "cluster E/electron p", HistType::kTH3F, {eoverpAxis, trackptAxis, nmatchedtrack});                                                                  // E over p vs track pT for the Nth closest electron track
    mHistManager.add("clusterTM_EoverP_p", "cluster E/electron p", HistType::kTH3F, {eoverpAxis, trackptAxis, nmatchedtrack});                                                                  // E over p vs track pT for the Nth closest positron track
    mHistManager.add("clusterTM_EoverP_electron_ASide", "cluster E/electron p in A-Side", HistType::kTH3F, {eoverpAxis, trackptAxis, nmatchedtrack});                                           // E over p vs track pT for the Nth closest electron/positron track in A-Side
    mHistManager.add("clusterTM_EoverP_electron_CSide", "cluster E/electron p in C-Side", HistType::kTH3F, {eoverpAxis, trackptAxis, nmatchedtrack});                                           // E over p vs track pT for the Nth closest electron/positron track in C-Side
    mHistManager.add("clusterTM_NSigma_BeforeCut", "electron NSigma for matched tracks before d#eta and #dvarphi cut", HistType::kTH3F, {nSigmaAxis, trackptAxis, nmatchedtrack});              // NSigma_electron vs track pT for the Nth closest track before dEta and dPhi cut
    mHistManager.add("clusterTM_NSigma_neg_BeforeCut", "electron NSigma for matched negative tracks before d#eta and #dvarphi cut", HistType::kTH3F, {nSigmaAxis, trackptAxis, nmatchedtrack}); // NSigma_electron vs track pT for the Nth closest negative tracks before dEta and dPhi cut
    mHistManager.add("clusterTM_NSigma_pos_BeforeCut", "electron NSigma for matched positive tracks before d#eta and #dvarphi cut", HistType::kTH3F, {nSigmaAxis, trackptAxis, nmatchedtrack}); // NSigma_electron vs track pT for the Nth closest positive tracks before dEta and dPhi cut
    mHistManager.add("clusterTM_NSigma", "electron NSigma for matched track", HistType::kTH3F, {nSigmaAxis, trackptAxis, nmatchedtrack});                                                       // NSigma_electron vs track pT for the Nth closest track
    mHistManager.add("clusterTM_NSigma_neg", "electron NSigma for matched negative track", HistType::kTH3F, {nSigmaAxis, trackptAxis, nmatchedtrack});                                          // NSigma_electron vs track pT for the Nth closest negative tracks
    mHistManager.add("clusterTM_NSigma_pos", "electron NSigma for matched positive track", HistType::kTH3F, {nSigmaAxis, trackptAxis, nmatchedtrack});                                          // NSigma_electron vs track pT for the Nth closest positive tracks
    mHistManager.add("clusterTM_NSigma_cut", "electron NSigma for matched track with cuts", HistType::kTH3F, {nSigmaAxis, trackptAxis, nmatchedtrack});                                         // NSigma_electron vs track pT for the Nth closest track with cuts on E/p and cluster cuts
    mHistManager.add("clusterTM_NSigma_e", "NSigma electron", HistType::kTH3F, {nSigmaAxis, trackptAxis, nmatchedtrack});                                                                       // NSigma vs track pT for the Nth closest electron track
    mHistManager.add("clusterTM_NSigma_p", "NSigma positron", HistType::kTH3F, {nSigmaAxis, trackptAxis, nmatchedtrack});                                                                       // NSigma vs track pT for the Nth closest positron track
    mHistManager.add("clusterTM_NSigma_e_cut", "NSigma electron with E over p cut", HistType::kTH3F, {nSigmaAxis, trackptAxis, nmatchedtrack});                                                 // NSigma vs track pT for the Nth closest electron track with E/p cut
    mHistManager.add("clusterTM_NSigma_p_cut", "NSigma positron with E over p cut", HistType::kTH3F, {nSigmaAxis, trackptAxis, nmatchedtrack});                                                 // NSigma vs track pT for the Nth closest positron track with E/p cut
    mHistManager.add("TrackTM_NSigma", "electron NSigma for global tracks", HistType::kTH2F, {nSigmaAxis, trackptAxis});                                                                        // NSigma_electron vs track pT for global track
    mHistManager.add("TrackTM_NSigma_e", "NSigma e for global negative tracks", HistType::kTH2F, {nSigmaAxis, trackptAxis});                                                                    // NSigma vs track pT for negative global track
    mHistManager.add("TrackTM_NSigma_p", "NSigma e for global positive tracks", HistType::kTH2F, {nSigmaAxis, trackptAxis});                                                                    // NSigma vs track pT for positive global track
    mHistManager.add("clusterTM_NSigma_electron_ASide", "NSigma electron in A-Side", HistType::kTH3F, {nSigmaAxis, trackptAxis, nmatchedtrack});                                                // NSigma vs track pT for the Nth closest electron/positron track in A-Side
    mHistManager.add("clusterTM_NSigma_electron_CSide", "NSigma positron in C-Side", HistType::kTH3F, {nSigmaAxis, trackptAxis, nmatchedtrack});                                                // NSigma vs track pT for the Nth closest electron/positron track in C-Side
    mHistManager.add("clusterTM_EoverP_hadron", "cluster E/hadron p", HistType::kTH3F, {eoverpAxis, trackpAxis, nmatchedtrack});                                                                // E over p vs track pT for the Nth closest hadron track
    mHistManager.add("clusterTM_EoverP_hn", "cluster E/hadron p", HistType::kTH3F, {eoverpAxis, trackpAxis, nmatchedtrack});                                                                    // E over p vs track pT for the Nth closest negative hadron track
    mHistManager.add("clusterTM_EoverP_hp", "cluster E/hadron p", HistType::kTH3F, {eoverpAxis, trackpAxis, nmatchedtrack});                                                                    // E over p vs track pT for the Nth closest positive hadron track
    mHistManager.add("clusterTM_EoverP_hadron_ASide", "cluster E/hadron p in A-Side", HistType::kTH3F, {eoverpAxis, trackpAxis, nmatchedtrack});                                                // E over p vs track pT for the Nth closest hadron track in A-Side
    mHistManager.add("clusterTM_EoverP_hadron_CSide", "cluster E/hadron p in C-Side", HistType::kTH3F, {eoverpAxis, trackpAxis, nmatchedtrack});                                                // E over p vs track pT for the Nth closest hadron track in C-Side

    // MC truth of the track matching, only booked if processMCTruth is enabled
    if (doprocessMCTruth.value) {
      static constexpr std::array<std::string_view, static_cast<int>(TmTruthClass::kNClasses)> ClassLabels = {"conversion e^{#pm}", "primary e^{#pm}", "other e^{#pm}", "charged (non e)", "#gamma", "neutral (non #gamma)", "other"};
      static constexpr std::array<std::string_view, static_cast<int>(TmTruthStatus::kNStatus)> StatusLabels = {"no reco track", "not extrapolated", "no matching candidate", "candidate not matched", "matched primary", "matched secondary"};
      static constexpr std::array<std::string_view, static_cast<int>(TmTrackType::kNTypes)> TypeLabels = {"no track", "global", "ITS-TPC non-global", "TPC only", "ITS only", "other"};
      const AxisSpec classAxis{static_cast<int>(TmTruthClass::kNClasses), -0.5, static_cast<int>(TmTruthClass::kNClasses) - 0.5, "leading MC particle"};
      const AxisSpec statusAxis{static_cast<int>(TmTruthStatus::kNStatus), -0.5, static_cast<int>(TmTruthStatus::kNStatus) - 0.5, "track matching status"};
      const AxisSpec typeAxis{static_cast<int>(TmTrackType::kNTypes), -0.5, static_cast<int>(TmTrackType::kNTypes) - 0.5, "track type"};
      const AxisSpec v0LegAxis{2, -0.5, 1.5, "V0 leg"};
      const AxisSpec mcPtAxis{std::vector<double>{0., 0.25, 0.5, 0.75, 1., 1.5, 2., 3., 4., 6., 10., 20., 50.}, "#it{p}_{T,MC} (GeV/#it{c})"};
      // residuals: fine binning around 0, coarse binning up to MaxResidual, so the tails beyond the matching windows are visible
      constexpr double MaxResidual = 0.4;  // range of the residual axes
      constexpr double FineRange = 0.04;   // fine binning within +-FineRange
      constexpr double FineWidth = 0.002;  // bin width within +-FineRange
      constexpr double CoarseWidth = 0.01; // bin width outside +-FineRange
      const int nCoarse = static_cast<int>(std::lround((MaxResidual - FineRange) / CoarseWidth));
      const int nFine = static_cast<int>(std::lround((FineRange - (-FineRange)) / FineWidth));
      std::vector<double> residualBins;
      residualBins.reserve(nCoarse + nFine + nCoarse + 1); // negative coarse, fine, positive coarse edges incl. the upper end
      for (int i = 0; i < nCoarse; ++i) {
        residualBins.push_back(-MaxResidual + CoarseWidth * i);
      }
      for (int i = 0; i < nFine; ++i) {
        residualBins.push_back(-FineRange + FineWidth * i);
      }
      for (int i = 0; i <= nCoarse; ++i) {
        residualBins.push_back(FineRange + CoarseWidth * i);
      }
      const AxisSpec mcDEtaAxis{residualBins, "#Delta#eta (track - cluster)"};
      const AxisSpec mcDPhiAxis{residualBins, "#Delta#varphi (track - cluster) (rad)"};
      const AxisSpec radiusAxis{125, 0., 250., "#it{R}_{xy} of the conversion (cm)"};

      auto setLabels = [](TAxis* axis, auto const& labels) {
        for (std::size_t i = 0; i < labels.size(); ++i) {
          axis->SetBinLabel(static_cast<int>(i) + 1, labels[i].data());
        }
      };
      // histograms with labelled axes are TH2 or THn: the TH3 label merge of ROOT 6.36 is broken
      auto hLeadingClass = mHistManager.add<TH2>("MCTruth/hLeadingClass", "leading MC particle of the cluster", HistType::kTH2F, {classAxis, energyAxis});
      setLabels(hLeadingClass->GetXaxis(), ClassLabels);
      auto hStatus = mHistManager.add<THn>("MCTruth/hStatus", "track matching status of the leading MC particle", HistType::kTHnF, {statusAxis, classAxis, typeAxis, v0LegAxis, mcPtAxis});
      setLabels(hStatus->GetAxis(0), StatusLabels);
      setLabels(hStatus->GetAxis(1), ClassLabels);
      setLabels(hStatus->GetAxis(2), TypeLabels);
      for (const auto& name : {"MCTruth/hResidualEta", "MCTruth/hResidualPhi"}) {
        const bool isEta = std::string_view(name).ends_with("Eta");
        auto hRes = mHistManager.add<THn>(name, "residual of the extrapolated track of the leading MC particle (origin frame)", HistType::kTHnF, {isEta ? mcDEtaAxis : mcDPhiAxis, mcPtAxis, classAxis, typeAxis});
        setLabels(hRes->GetAxis(2), ClassLabels);
        setLabels(hRes->GetAxis(3), TypeLabels);
      }
      mHistManager.add("MCTruth/hResidualConversion", "conversion e^{#pm}: residual of the extrapolated track", HistType::kTH3F, {mcDEtaAxis, mcDPhiAxis, mcPtAxis});
      mHistManager.add("MCTruth/hResidualVsStored", "matched tracks: own minus stored residual;#Delta#eta_{own} - #Delta#eta_{stored};#Delta#varphi_{own} - #Delta#varphi_{stored}", HistType::kTH2F, {{200, -0.02, 0.02}, {200, -0.02, 0.02}});
      auto hConvRadius = mHistManager.add<THn>("MCTruth/hConvRadius", "conversion e^{#pm}: matching status vs conversion radius", HistType::kTHnF, {statusAxis, radiusAxis, mcPtAxis});
      setLabels(hConvRadius->GetAxis(0), StatusLabels);
      mHistManager.add("MCTruth/hEoverPConversion", "conversion e^{#pm}: #it{E}_{cls}/#it{p}_{track}", HistType::kTH2F, {eoverpAxis, trackptAxis});
    }

    if (!vetoBCID->empty()) {
      std::stringstream parser(vetoBCID.value);
      std::string token;
      int bcid = -1;
      while (std::getline(parser, token, ',')) {
        bcid = std::stoi(token);
        LOG(info) << "Veto BCID " << bcid;
        mVetoBCIDs.push_back(bcid);
      }
    }
    if (selectBCID.value != "all") {
      std::stringstream parser(selectBCID.value);
      std::string token;
      int bcid = -1;
      while (std::getline(parser, token, ',')) {
        bcid = std::stoi(token);
        LOG(info) << "Select BCID " << bcid;
        mSelectBCIDs.push_back(bcid);
      }
    }
  }

  // define cluster filter. It selects only those clusters which are of the type
  // sadly passing of the string at runtime is not possible for technical region so cluster definition is
  // an integer instead
  Filter clusterDefinitionSelection = (o2::aod::emcalcluster::definition == clusterDefinition) && (o2::aod::emcalcluster::time >= minTime) && (o2::aod::emcalcluster::time <= maxTime) && (o2::aod::emcalcluster::m02 > minM02) && (o2::aod::emcalcluster::m02 < maxM02);

  /// \brief Centrality / multiplicity of the collision for the estimator of eventCuts.centMultEstimator
  static float getCentMult(CollisionEvSelIt const& collision, const CentMultEstimator estimator)
  {
    switch (estimator) {
      case CentMultEstimator::kCentFT0C:
        return collision.centFT0C();
      case CentMultEstimator::kCentFT0M:
        return collision.centFT0M();
      case CentMultEstimator::kCentFT0A:
        return collision.centFT0A();
      case CentMultEstimator::kMultFT0C:
        return collision.multFT0C();
      case CentMultEstimator::kMultFT0A:
        return collision.multFT0A();
      case CentMultEstimator::kMultFT0M:
        return collision.multFT0M();
      case CentMultEstimator::kMultNTracksPV:
        return collision.multNTracksPV();
      default:
        return 0.f;
    }
  }

  /// \brief Whether a value lies in [min, max], a negative limit is off
  static bool isInRange(const float value, const float min, const float max)
  {
    return (min < 0.f || value >= min) && (max < 0.f || value <= max);
  }

  /// \brief BC ID veto / selection of vetoBCID and selectBCID
  bool isBCSelected(const int bcId) const
  {
    if (std::find(mVetoBCIDs.begin(), mVetoBCIDs.end(), bcId) != mVetoBCIDs.end()) {
      return false;
    }
    return mSelectBCIDs.empty() || std::find(mSelectBCIDs.begin(), mSelectBCIDs.end(), bcId) != mSelectBCIDs.end();
  }

  /// \brief Event selection of the monitor, see eventCuts. The cuts are applied in the order of EventCut.
  /// \param theCollision collision to check
  /// \param fillHistos fill the event histograms (only once if several process functions run)
  bool isEventSelected(CollisionEvSelIt const& theCollision, const bool fillHistos)
  {
    // fill the cut flow bin of a passed step
    auto passed = [&](const EventCut cut) {
      if (fillHistos) {
        mHistManager.fill(HIST("Event/hCutFlow"), static_cast<int>(cut));
      }
    };
    const auto estimator = static_cast<CentMultEstimator>(eventCuts.centMultEstimator.value);
    const bool applyCentMult = estimator != CentMultEstimator::kOff;
    if (fillHistos) {
      mHistManager.fill(HIST("eventsAll"), 1);
      mHistManager.fill(HIST("eventVertexZAll"), theCollision.posZ());
      mHistManager.fill(HIST("Event/hTrackOccupancyBefore"), theCollision.trackOccupancyInTimeRange());
      mHistManager.fill(HIST("Event/hFT0COccupancyBefore"), theCollision.ft0cOccupancyInTimeRange());
      if (applyCentMult) {
        mHistManager.fill(HIST("Event/hCentMultBefore"), getCentMult(theCollision, estimator));
      }
    }
    passed(EventCut::kAll);

    // BC ID veto / selection (pilot beams); needs aod::BCs in the process function, otherwise the bc() index is not bound
    const auto bc = theCollision.bc();
    const int bcId = static_cast<int>(bc.globalBC() % o2::constants::lhc::LHCMaxBunches);
    if (!isBCSelected(bcId)) {
      return false;
    }
    passed(EventCut::kBC);

    if (eventCuts.requireSel8.value && !theCollision.sel8()) {
      return false;
    }
    passed(EventCut::kSel8);

    // EMCal in the readout: kINT7 in Run 2, kTVXinEMC in Run 3
    if (eventCuts.requireEmcalReadout.value || doEventSel.value) {
      constexpr int BeginningRun3 = 300000;
      const bool hasEmcal = bc.runNumber() < BeginningRun3 ? theCollision.alias_bit(kINT7) : theCollision.alias_bit(kTVXinEMC);
      if (!hasEmcal) {
        return false;
      }
    }
    passed(EventCut::kEmcalReadout);

    const double zVertexMax = vertexCut.value > 0. ? vertexCut.value : static_cast<double>(eventCuts.zVertexMax.value);
    if (zVertexMax > 0. && std::abs(theCollision.posZ()) > zVertexMax) {
      return false;
    }
    passed(EventCut::kZVertex);

    // evsel flags, in the order of EventCut
    const std::array<std::pair<bool, o2::aod::evsel::EventSelectionFlags>, static_cast<int>(EventCut::kGoodITSLayersAll) - static_cast<int>(EventCut::kNoTimeFrameBorder) + 1> flags = {{
      {eventCuts.requireNoTimeFrameBorder.value, o2::aod::evsel::kNoTimeFrameBorder},
      {eventCuts.requireNoITSROFrameBorder.value, o2::aod::evsel::kNoITSROFrameBorder},
      {eventCuts.requireNoSameBunchPileup.value, o2::aod::evsel::kNoSameBunchPileup},
      {eventCuts.requireGoodZvtxFT0vsPV.value, o2::aod::evsel::kIsGoodZvtxFT0vsPV},
      {eventCuts.requireVertexITSTPC.value, o2::aod::evsel::kIsVertexITSTPC},
      {eventCuts.requireVertexTOFmatched.value, o2::aod::evsel::kIsVertexTOFmatched},
      {eventCuts.requireNoCollInTimeRangeStandard.value, o2::aod::evsel::kNoCollInTimeRangeStandard},
      {eventCuts.requireNoCollInTimeRangeNarrow.value, o2::aod::evsel::kNoCollInTimeRangeNarrow},
      {eventCuts.requireNoCollInRofStandard.value, o2::aod::evsel::kNoCollInRofStandard},
      {eventCuts.requireNoHighMultCollInPrevRof.value, o2::aod::evsel::kNoHighMultCollInPrevRof},
      {eventCuts.requireGoodITSLayer3.value, o2::aod::evsel::kIsGoodITSLayer3},
      {eventCuts.requireGoodITSLayer0123.value, o2::aod::evsel::kIsGoodITSLayer0123},
      {eventCuts.requireGoodITSLayersAll.value, o2::aod::evsel::kIsGoodITSLayersAll},
    }};
    for (std::size_t i = 0; i < flags.size(); ++i) {
      if (flags[i].first && !theCollision.selection_bit(flags[i].second)) {
        return false;
      }
      passed(static_cast<EventCut>(static_cast<int>(EventCut::kNoTimeFrameBorder) + static_cast<int>(i)));
    }

    // occupancy, track and FT0C based, each with its own range
    if (!isInRange(static_cast<float>(theCollision.trackOccupancyInTimeRange()), eventCuts.trackOccupancyMin.value, eventCuts.trackOccupancyMax.value)) {
      return false;
    }
    passed(EventCut::kTrackOccupancy);
    if (!isInRange(theCollision.ft0cOccupancyInTimeRange(), eventCuts.ft0cOccupancyMin.value, eventCuts.ft0cOccupancyMax.value)) {
      return false;
    }
    passed(EventCut::kFT0COccupancy);

    if (eventCuts.applyRCT.value && !mRctChecker.checkTable(theCollision)) {
      return false;
    }
    passed(EventCut::kRCT);

    if (applyCentMult) {
      const float centMult = getCentMult(theCollision, estimator);
      if (centMult < eventCuts.centMultMin.value || centMult > eventCuts.centMultMax.value) {
        return false;
      }
    }
    passed(EventCut::kCentMult);

    passed(EventCut::kSelected);
    if (fillHistos) {
      mHistManager.fill(HIST("eventsSelected"), 1);
      mHistManager.fill(HIST("eventVertexZSelected"), theCollision.posZ());
      mHistManager.fill(HIST("Event/hTrackOccupancyAfter"), theCollision.trackOccupancyInTimeRange());
      mHistManager.fill(HIST("Event/hFT0COccupancyAfter"), theCollision.ft0cOccupancyInTimeRange());
      if (applyCentMult) {
        mHistManager.fill(HIST("Event/hCentMultAfter"), getCentMult(theCollision, estimator));
      }
      mHistManager.fill(HIST("Event/hBCId"), bcId);
    }
    return true;
  }

  /// \brief Process EMCAL clusters that are matched to a collisions
  void processCollisions(CollisionEvSelIt const& theCollision, o2::aod::BCs const&, SelectedClusters const& clusters, o2::aod::EMCALClusterCells const&, o2::aod::Calos const&, o2::aod::EMCALMatchedTracks const& matchedtracks, TracksPID const& alltracks)
  {
    if (!isEventSelected(theCollision, true)) {
      return;
    }
    float lowP = 0.75f;
    float midP = 1.25f;

    for (const auto& alltrack : alltracks) {
      double trackEta = alltrack.trackEtaEmcal();
      double trackPhi = alltrack.trackPhiEmcal();
      if (alltrack.isGlobalTrack()) { // NSigma of all global tracks without matching
        mHistManager.fill(HIST("TrackTM_NSigma"), alltrack.tpcNSigmaEl(), alltrack.pt());
        mHistManager.fill(HIST("TrackEtaPhi"), trackEta, trackPhi);
        if (alltrack.sign() == 1) { // NSigma of all global positive tracks without matching
          mHistManager.fill(HIST("TrackTM_NSigma_p"), alltrack.tpcNSigmaEl(), alltrack.pt());
          mHistManager.fill(HIST("TrackEtaPhi_Pos"), trackEta, trackPhi);
        } else if (alltrack.sign() == -1) { // NSigma of all global negative tracks without matching
          mHistManager.fill(HIST("TrackTM_NSigma_e"), alltrack.tpcNSigmaEl(), alltrack.pt());
          mHistManager.fill(HIST("TrackEtaPhi_Neg"), trackEta, trackPhi);
        }
      }
    }

    LOG(debug) << "bunch crossing ID" << theCollision.bcId();
    // loop over all clusters from accepted collision
    for (const auto& cluster : clusters) {
      // fill histograms of cluster properties
      // in this implementation the cluster properties are directly
      // loaded from the flat table, in the future one should
      // consider using the AnalysisCluster object to work with
      // after loading.
      LOG(debug) << "Cluster energy: " << cluster.energy();
      LOG(debug) << "Cluster time: " << cluster.time();
      LOG(debug) << "Cluster M02: " << cluster.m02();
      LOG(debug) << "Cluster energy: " << cluster.energy();
      LOG(debug) << "Cluster index: " << cluster.index();
      LOG(debug) << "ncells in cluster: " << cluster.nCells();
      LOG(debug) << "real cluster id" << cluster.id();

      // Example of loop over tracks matched to cluster within dR=0.4, where only the
      // 20 most closest tracks are stored. If needed the number of tracks can be later
      // increased in the correction framework. Access to all track properties via match.track()
      // If you want to access tracks from a joined table you can access all track properties via
      // match.track_as<globTracks>() with
      // using globTracks = o2::soa::Join<o2::aod::Tracks, o2::aod::TrackSelection>;
      // In this example the counter t is just used to only look at the closest match
      float dEta{0}, dPhi{0}, pT{0}, absP{0}, trackEta{0}, trackPhi{0}, nSigmaEl{0};
      int supermoduleID = -1;
      try {
        supermoduleID = mGeometry->SuperModuleNumberFromEtaPhi(cluster.eta(), cluster.phi());
      } catch (o2::emcal::InvalidPositionException const& e) {
        // Imprecision of the position at the sector boundaries, mostly due to
        // vertex imprecision. Skip these clusters for the now.
        continue;
      }

      pT = cluster.energy() / std::cosh(cluster.eta());
      if (m02highPt > 0. && cluster.m02() >= maxM02HighPt && pT >= m02highPt) { // high pT M02 cut
        continue;
      }
      auto tracksofcluster = matchedtracks.sliceBy(perClusterMatchedTracks, cluster.globalIndex());
      int t = 0;
      for (const auto& match : tracksofcluster) {
        nSigmaEl = match.track_as<TracksPID>().tpcNSigmaEl();
        // exmple of how to access any property of the matched tracks (tracks are sorted by how close they are to cluster)
        LOG(debug) << "Pt of match" << match.track_as<TracksPID>().pt();
        absP = std::abs(match.track_as<TracksPID>().p());
        // only consider closest match
        trackEta = match.track_as<TracksPID>().trackEtaEmcal();
        trackPhi = match.track_as<TracksPID>().trackPhiEmcal();
        dPhi = match.deltaPhi();
        dEta = match.deltaEta();
        mHistManager.fill(HIST("MatchedTrackEtaPhi_BeforeCut"), trackEta, trackPhi);
        mHistManager.fill(HIST("clusterTM_NSigma_BeforeCut"), nSigmaEl, match.track_as<TracksPID>().pt(), t);
        if (match.track_as<TracksPID>().sign() == -1) {
          mHistManager.fill(HIST("MatchedTrackEtaPhi_Neg_BeforeCut"), trackEta, trackPhi);
          mHistManager.fill(HIST("clusterTM_NSigma_neg_BeforeCut"), nSigmaEl, match.track_as<TracksPID>().pt(), t);
        } else if (match.track_as<TracksPID>().sign() == 1) {
          mHistManager.fill(HIST("MatchedTrackEtaPhi_Pos_BeforeCut"), trackEta, trackPhi);
          mHistManager.fill(HIST("clusterTM_NSigma_pos_BeforeCut"), nSigmaEl, match.track_as<TracksPID>().pt(), t);
        }
        if (std::abs(dEta) >= minDEta || std::abs(dPhi) >= minDPhi) { // dEta and dPhi cut
          continue;
        }
        if (hasTRD && !(match.track_as<TracksPID>().hasTRD())) { // request TRD hit cut
          continue;
        }
        // only fill these for the first matched track:
        if (t == 0) {
          mHistManager.fill(HIST("clusterTM_dEtadPhi"), dEta, dPhi, supermoduleID);
          mHistManager.fill(HIST("clusterTM_dEtaPt"), dEta, pT, supermoduleID);
        }
        double eOverP = cluster.energy() / absP;
        mHistManager.fill(HIST("clusterTM_dEtaTN"), dEta, t);
        mHistManager.fill(HIST("clusterTM_dPhiTN"), dPhi, t);
        mHistManager.fill(HIST("clusterTM_dRTN"), std::sqrt(dPhi * dPhi + dEta * dEta), t);
        mHistManager.fill(HIST("clusterTM_EoverP_E"), eOverP, cluster.energy(), t);
        mHistManager.fill(HIST("clusterTM_dEtadPhi"), dEta, dPhi, t);
        mHistManager.fill(HIST("clusterEMatched"), cluster.energy(), t);
        mHistManager.fill(HIST("clusterTM_EvsP"), cluster.energy(), absP, t);
        mHistManager.fill(HIST("clusterTM_EoverP_Pt"), eOverP, match.track_as<TracksPID>().pt(), t);
        mHistManager.fill(HIST("clusterTM_NSigma"), nSigmaEl, match.track_as<TracksPID>().pt(), t);
        mHistManager.fill(HIST("MatchedTrackEtaPhi"), trackEta, trackPhi);
        if (eOverPRange->at(0) <= eOverP && eOverP <= eOverPRange->at(1)) {
          mHistManager.fill(HIST("clusterTM_NSigma_cut"), nSigmaEl, match.track_as<TracksPID>().pt(), t);
        }
        // A- and C-side
        if (match.track_as<TracksPID>().eta() > 0.0 && t == 0) {
          mHistManager.fill(HIST("clusterTM_dEtadPhi_ASide"), dEta, dPhi, supermoduleID);
        } else if (match.track_as<TracksPID>().eta() < 0.0 && t == 0) {
          mHistManager.fill(HIST("clusterTM_dEtadPhi_CSide"), dEta, dPhi, supermoduleID);
        }
        // positive and negative tracks seperate, with three different track momentum ranges
        if (match.track_as<TracksPID>().sign() == 1) {
          mHistManager.fill(HIST("clusterTM_NSigma_pos"), nSigmaEl, match.track_as<TracksPID>().pt(), t);
          mHistManager.fill(HIST("MatchedTrackEtaPhi_Pos"), trackEta, trackPhi);
          if (t == 0) {
            mHistManager.fill(HIST("clusterTM_PosdPhiPt"), dPhi, pT, supermoduleID);
            mHistManager.fill(HIST("clusterTM_PosdEtadPhi"), dEta, dPhi, supermoduleID);
            if (absP < lowP) {
              mHistManager.fill(HIST("clusterTM_PosdEtadPhi_Pl0_75"), dEta, dPhi, supermoduleID);
            } else if (absP >= midP) {
              mHistManager.fill(HIST("clusterTM_PosdEtadPhi_Pgeq1_25"), dEta, dPhi, supermoduleID);
            } else {
              mHistManager.fill(HIST("clusterTM_PosdEtadPhi_0_75leqPl1_25"), dEta, dPhi, supermoduleID);
            }
          }
        } else if (match.track_as<TracksPID>().sign() == -1) {
          mHistManager.fill(HIST("clusterTM_NSigma_neg"), nSigmaEl, match.track_as<TracksPID>().pt(), t);
          mHistManager.fill(HIST("MatchedTrackEtaPhi_Neg"), trackEta, trackPhi);
          if (t == 0) {
            mHistManager.fill(HIST("clusterTM_NegdPhiPt"), dPhi, pT, supermoduleID);
            mHistManager.fill(HIST("clusterTM_NegdEtadPhi"), dEta, dPhi, supermoduleID);
            if (absP < lowP) {
              mHistManager.fill(HIST("clusterTM_NegdEtadPhi_Pl0_75"), dEta, dPhi, supermoduleID);
            } else if (absP >= midP) {
              mHistManager.fill(HIST("clusterTM_NegdEtadPhi_Pgeq1_25"), dEta, dPhi, supermoduleID);
            } else {
              mHistManager.fill(HIST("clusterTM_NegdEtadPhi_0_75leqPl1_25"), dEta, dPhi, supermoduleID);
            }
          }
        }
        if (tpcNsigmaElectron->at(0) <= nSigmaEl && nSigmaEl <= tpcNsigmaElectron->at(1)) {                                                                                   // E/p for e+/e-
          if (usePionRejection && (tpcNsigmaPion->at(0) <= match.track_as<TracksPID>().tpcNSigmaPi() || match.track_as<TracksPID>().tpcNSigmaPi() <= tpcNsigmaPion->at(1))) { // with pion rejection
            mHistManager.fill(HIST("clusterTM_EoverP_ep"), eOverP, match.track_as<TracksPID>().pt(), t);
            if (match.track_as<TracksPID>().eta() >= 0.) {
              mHistManager.fill(HIST("clusterTM_EoverP_electron_ASide"), eOverP, match.track_as<TracksPID>().pt(), t);
              mHistManager.fill(HIST("clusterTM_NSigma_electron_ASide"), nSigmaEl, match.track_as<TracksPID>().pt(), t);
            } else {
              mHistManager.fill(HIST("clusterTM_EoverP_electron_CSide"), eOverP, match.track_as<TracksPID>().pt(), t);
              mHistManager.fill(HIST("clusterTM_NSigma_electron_CSide"), nSigmaEl, match.track_as<TracksPID>().pt(), t);
            }
            if (match.track_as<TracksPID>().sign() == -1) {
              mHistManager.fill(HIST("clusterTM_EoverP_e"), eOverP, match.track_as<TracksPID>().pt(), t);
              mHistManager.fill(HIST("clusterTM_NSigma_e"), nSigmaEl, match.track_as<TracksPID>().pt(), t);
              if (eOverPRange->at(0) <= eOverP && eOverP <= eOverPRange->at(1)) {
                mHistManager.fill(HIST("clusterTM_NSigma_e_cut"), nSigmaEl, match.track_as<TracksPID>().pt(), t);
              }
            } else if (match.track_as<TracksPID>().sign() == +1) {
              mHistManager.fill(HIST("clusterTM_EoverP_p"), eOverP, match.track_as<TracksPID>().pt(), t);
              mHistManager.fill(HIST("clusterTM_NSigma_p"), nSigmaEl, match.track_as<TracksPID>().pt(), t);
              if (eOverPRange->at(0) <= eOverP && eOverP <= eOverPRange->at(1)) {
                mHistManager.fill(HIST("clusterTM_NSigma_p_cut"), nSigmaEl, match.track_as<TracksPID>().pt(), t);
              }
            }
          } else { // without pion rejection
            mHistManager.fill(HIST("clusterTM_EoverP_ep"), eOverP, match.track_as<TracksPID>().pt(), t);
            if (match.track_as<TracksPID>().eta() >= 0.) {
              mHistManager.fill(HIST("clusterTM_EoverP_electron_ASide"), eOverP, match.track_as<TracksPID>().pt(), t);
              mHistManager.fill(HIST("clusterTM_NSigma_electron_ASide"), nSigmaEl, match.track_as<TracksPID>().pt(), t);
            } else {
              mHistManager.fill(HIST("clusterTM_EoverP_electron_CSide"), eOverP, match.track_as<TracksPID>().pt(), t);
              mHistManager.fill(HIST("clusterTM_NSigma_electron_CSide"), nSigmaEl, match.track_as<TracksPID>().pt(), t);
            }
            if (match.track_as<TracksPID>().sign() == -1) {
              mHistManager.fill(HIST("clusterTM_EoverP_e"), eOverP, match.track_as<TracksPID>().pt(), t);
              mHistManager.fill(HIST("clusterTM_NSigma_e"), nSigmaEl, match.track_as<TracksPID>().pt(), t);
              if (eOverPRange->at(0) <= eOverP && eOverP <= eOverPRange->at(1)) {
                mHistManager.fill(HIST("clusterTM_NSigma_e_cut"), nSigmaEl, match.track_as<TracksPID>().pt(), t);
              }
            } else if (match.track_as<TracksPID>().sign() == +1) {
              mHistManager.fill(HIST("clusterTM_EoverP_p"), eOverP, match.track_as<TracksPID>().pt(), t);
              mHistManager.fill(HIST("clusterTM_NSigma_p"), nSigmaEl, match.track_as<TracksPID>().pt(), t);
              if (eOverPRange->at(0) <= eOverP && eOverP <= eOverPRange->at(1)) {
                mHistManager.fill(HIST("clusterTM_NSigma_p_cut"), nSigmaEl, match.track_as<TracksPID>().pt(), t);
              }
            }
          }
        } else if (tpcNsigmaBack->at(0) <= nSigmaEl && nSigmaEl <= tpcNsigmaBack->at(1)) { // E/p for hadrons / background
          mHistManager.fill(HIST("clusterTM_EoverP_hadron"), eOverP, match.track_as<TracksPID>().pt(), t);
          if (match.track_as<TracksPID>().eta() >= 0.) {
            mHistManager.fill(HIST("clusterTM_EoverP_hadron_ASide"), eOverP, match.track_as<TracksPID>().pt(), t);
          } else {
            mHistManager.fill(HIST("clusterTM_EoverP_hadron_CSide"), eOverP, match.track_as<TracksPID>().pt(), t);
          }
          if (match.track_as<TracksPID>().sign() == -1) {
            mHistManager.fill(HIST("clusterTM_EoverP_hn"), eOverP, match.track_as<TracksPID>().pt(), t);
          } else if (match.track_as<TracksPID>().sign() == +1) {
            mHistManager.fill(HIST("clusterTM_EoverP_hp"), eOverP, match.track_as<TracksPID>().pt(), t);
          }
        }
        t++;
      }
      mHistManager.fill(HIST("clusterTM_NTrack"), t);
    }
  }
  PROCESS_SWITCH(EmcTmMonitor, processCollisions, "Process clusters from collision", true);

  /// \brief Information of a reconstructed track needed for the MC truth of the matching
  struct McTrackInfo {
    int64_t id{-1};
    bool isGlobal{false};
    bool hasITS{false};
    bool hasTPC{false};
    float etaEmcal{-999.f};
    float phiEmcal{-999.f};
    float p{0.f};
    float pt{0.f};
  };

  static constexpr float NotOnEmcal = -900.f; // trackEta/PhiEmcal of tracks that were not extrapolated to the EMCal are -999

  /// \brief Truth class of an MC particle
  template <o2::soa::is_iterator TMcParticle, o2::soa::is_table TMcParticles>
  TmTruthClass truthClassOf(TMcParticle const& mc, TMcParticles const& mcParticles)
  {
    const int pdg = mc.pdgCode();
    if (std::abs(pdg) == PDG_t::kElectron) {
      if (mc.producedByGenerator()) {
        return TmTruthClass::kPrimaryElectron;
      }
      if (mc.getProcess() == TMCProcess::kPPair && mc.has_mothers() && mcParticles.iteratorAt(mc.mothersIds()[0]).pdgCode() == PDG_t::kGamma) {
        return TmTruthClass::kConversionElectron;
      }
      return TmTruthClass::kOtherElectron;
    }
    if (pdg == PDG_t::kGamma) {
      return TmTruthClass::kPhoton;
    }
    const TParticlePDG* info = mPdgDatabase->GetParticle(pdg);
    if (info == nullptr) {
      return TmTruthClass::kOther;
    }
    return info->Charge() != 0. ? TmTruthClass::kOtherCharged : TmTruthClass::kNeutral;
  }

  /// \brief Type of a reconstructed track
  static TmTrackType trackTypeOf(McTrackInfo const& info)
  {
    if (info.isGlobal) {
      return TmTrackType::kGlobal;
    }
    if (info.hasITS && info.hasTPC) {
      return TmTrackType::kItsTpc;
    }
    if (info.hasTPC) {
      return TmTrackType::kTpcOnly;
    }
    return info.hasITS ? TmTrackType::kItsOnly : TmTrackType::kOtherType;
  }

  /// \brief MC truth of the track matching: follow the leading MC particle of each cluster through the matching chain.
  ///        Needs emcal-correction-task with processMCWithSecondaries and the photon conversion builder (aod::V0Legs).
  void processMCTruth(CollisionEvSelIt const& theCollision, o2::aod::BCs const&, SelectedMcClusters const& clusters, o2::aod::EMCALMatchedTracks const& matchedTracks, o2::aod::EMCMatchSecs const& matchedSecondaries, TracksMC const& tracks, o2::aod::V0Legs const& v0Legs, o2::aod::McParticles const& mcParticles)
  {
    if (!isEventSelected(theCollision, !doprocessCollisions.value)) {
      return;
    }

    // reconstructed tracks of each MC particle in this collision, and the tracks used as V0 legs (secondary matching)
    std::unordered_map<int64_t, std::vector<McTrackInfo>> tracksOfMcParticle;
    for (const auto& track : tracks) {
      if (!track.has_mcParticle()) {
        continue;
      }
      tracksOfMcParticle[track.mcParticleId()].push_back({track.globalIndex(), track.isGlobalTrack(), track.hasITS(), track.hasTPC(), track.trackEtaEmcal(), track.trackPhiEmcal(), track.p(), track.pt()});
    }
    std::unordered_set<int64_t> v0LegTracks;
    for (const auto& leg : v0Legs) {
      v0LegTracks.insert(leg.trackId());
    }

    for (const auto& cluster : clusters) {
      try {
        mGeometry->SuperModuleNumberFromEtaPhi(cluster.eta(), cluster.phi());
      } catch (o2::emcal::InvalidPositionException const& e) {
        continue; // same as in processCollisions
      }
      const float clusterPt = cluster.energy() / std::cosh(cluster.eta());
      if (m02highPt > 0. && cluster.m02() >= maxM02HighPt && clusterPt >= m02highPt) {
        continue;
      }

      // leading contributor: the labels are ordered by energy fraction (ClusterLabel::orderLabels in the ClusterFactory)
      const auto mcIds = cluster.mcParticleIds();
      const int64_t leadId = mcIds.empty() ? -1 : mcIds[0];
      TmTruthClass truthClass = TmTruthClass::kOther;
      float mcPt = 0.f;
      float convRadius = -1.f;
      if (leadId >= 0) {
        const auto mc = mcParticles.iteratorAt(leadId);
        mcPt = mc.pt();
        truthClass = truthClassOf(mc, mcParticles);
        convRadius = std::hypot(mc.vx(), mc.vy());
      }
      mHistManager.fill(HIST("MCTruth/hLeadingClass"), static_cast<int>(truthClass), cluster.energy());
      const bool isConversion = truthClass == TmTruthClass::kConversionElectron;
      if (truthClass != TmTruthClass::kConversionElectron && truthClass != TmTruthClass::kPrimaryElectron && truthClass != TmTruthClass::kOtherElectron && truthClass != TmTruthClass::kOtherCharged) {
        continue; // neutral leading particle: no track to match
      }

      // best reconstructed track of the leading particle: extrapolated and matching candidate first, then highest pT
      const McTrackInfo* best = nullptr;
      int bestScore = -1;
      if (auto it = tracksOfMcParticle.find(leadId); it != tracksOfMcParticle.end()) {
        for (const auto& info : it->second) {
          const bool isExtrapolated = info.etaEmcal > NotOnEmcal && info.phiEmcal > NotOnEmcal;
          const bool isCandidate = info.isGlobal || v0LegTracks.contains(info.id);
          const int score = 2 * (isExtrapolated ? 1 : 0) + (isCandidate ? 1 : 0);
          if (score > bestScore || (score == bestScore && info.pt > best->pt)) {
            best = &info;
            bestScore = score;
          }
        }
      }

      TmTruthStatus status = TmTruthStatus::kNoRecoTrack;
      TmTrackType type = TmTrackType::kNoTrack;
      bool isV0Leg = false;
      if (best != nullptr) {
        type = trackTypeOf(*best);
        isV0Leg = v0LegTracks.contains(best->id);
        const bool isExtrapolated = best->etaEmcal > NotOnEmcal && best->phiEmcal > NotOnEmcal;

        // is the track among the primary or secondary matches of this cluster
        bool isMatchedPrimary = false, isMatchedSecondary = false;
        float storedDEta = 0.f, storedDPhi = 0.f;
        for (const auto& match : matchedTracks.sliceBy(perClusterMatchedTracks, cluster.globalIndex())) {
          if (match.trackId() == best->id) {
            isMatchedPrimary = true;
            storedDEta = match.deltaEta();
            storedDPhi = match.deltaPhi();
            break;
          }
        }
        if (!isMatchedPrimary) {
          for (const auto& match : matchedSecondaries.sliceBy(perClusterMatchedSecondaries, cluster.globalIndex())) {
            if (match.trackId() == best->id) {
              isMatchedSecondary = true;
              storedDEta = match.deltaEta();
              storedDPhi = match.deltaPhi();
              break;
            }
          }
        }
        if (isMatchedPrimary) {
          status = TmTruthStatus::kMatchedPrimary;
        } else if (isMatchedSecondary) {
          status = TmTruthStatus::kMatchedSecondary;
        } else if (!isExtrapolated) {
          status = TmTruthStatus::kNotExtrapolated;
        } else if (!best->isGlobal && !isV0Leg) {
          status = TmTruthStatus::kNoMatchingCandidate;
        } else {
          status = TmTruthStatus::kCandidateNotMatched;
        }

        if (isExtrapolated) {
          // the matching compares to the cluster position seen from the origin, the cluster table has eta seen from the vertex
          const float clusterEtaOrigin = std::asinh(std::sinh(cluster.eta()) + theCollision.posZ() / emcalRadius);
          const float dEta = best->etaEmcal - clusterEtaOrigin;
          const float dPhi = std::remainder(best->phiEmcal - cluster.phi(), o2::constants::math::TwoPI);
          mHistManager.fill(HIST("MCTruth/hResidualEta"), dEta, mcPt, static_cast<int>(truthClass), static_cast<int>(type));
          mHistManager.fill(HIST("MCTruth/hResidualPhi"), dPhi, mcPt, static_cast<int>(truthClass), static_cast<int>(type));
          if (isConversion) {
            mHistManager.fill(HIST("MCTruth/hResidualConversion"), dEta, dPhi, mcPt);
          }
          if (isMatchedPrimary || isMatchedSecondary) {
            mHistManager.fill(HIST("MCTruth/hResidualVsStored"), dEta - storedDEta, dPhi - storedDPhi);
          }
        }
        if (isConversion && best->p > 0.f) {
          mHistManager.fill(HIST("MCTruth/hEoverPConversion"), cluster.energy() / best->p, best->pt);
        }
      }
      mHistManager.fill(HIST("MCTruth/hStatus"), static_cast<int>(status), static_cast<int>(truthClass), static_cast<int>(type), static_cast<int>(isV0Leg), mcPt);
      if (isConversion) {
        mHistManager.fill(HIST("MCTruth/hConvRadius"), static_cast<int>(status), convRadius, mcPt);
      }
    }
  }
  PROCESS_SWITCH(EmcTmMonitor, processMCTruth, "MC truth of the track matching for clusters from charged particles", false);

  /// \brief Create binning for cluster energy axis (variable bin size)
  /// \return vector with bin limits
  std::vector<double> makeEnergyBinning() const
  {
    auto fillBinLimits = [](std::vector<double>& binlimits, double max, double binwidth) {
      auto current = *binlimits.rbegin();
      while (current < max) {
        current += binwidth;
        binlimits.emplace_back(current);
      }
    };
    std::vector<double> result = {0.};
    fillBinLimits(result, 2., 0.1);
    fillBinLimits(result, 5., 0.2);
    fillBinLimits(result, 10., 0.5);
    fillBinLimits(result, 20., 1.);
    fillBinLimits(result, 50., 2.);
    fillBinLimits(result, 100., 5.);
    fillBinLimits(result, 200., 10.);
    return result;
  }

  /// \brief Create binning for cluster energy axis (variable bin size)
  /// direct port from binning often used in AliPhysics for debugging
  /// \return vector with bin limits
  std::vector<double> makeEnergyBinningAliPhysics() const
  {

    std::vector<double> result;
    int nBinsClusterE = 235;
    for (int i = 0; i < nBinsClusterE + 1; i++) {
      if (i < 1) { // o2-linter: disable=magic-number (just numbers for binning)
        result.emplace_back(0.3 * i);
      } else if (i < 55) { // o2-linter: disable=magic-number (just numbers for binning)
        result.emplace_back(0.3 + 0.05 * (i - 1));
      } else if (i < 105) { // o2-linter: disable=magic-number (just numbers for binning)
        result.emplace_back(3. + 0.1 * (i - 55));
      } else if (i < 140) { // o2-linter: disable=magic-number (just numbers for binning)
        result.emplace_back(8. + 0.2 * (i - 105));
      } else if (i < 170) { // o2-linter: disable=magic-number (just numbers for binning)
        result.emplace_back(15. + 0.5 * (i - 140));
      } else if (i < 190) { // o2-linter: disable=magic-number (just numbers for binning)
        result.emplace_back(30. + 1.0 * (i - 170));
      } else if (i < 215) { // o2-linter: disable=magic-number (just numbers for binning)
        result.emplace_back(50. + 2.0 * (i - 190));
      } else if (i < 235) { // o2-linter: disable=magic-number (just numbers for binning)
        result.emplace_back(100. + 5.0 * (i - 215));
      } else {
        result.emplace_back(200. + 10.0 * (i - 235));
      }
    }
    return result;
  }

  /// \brief Create binning for pT axis (variable bin size)
  /// direct port from binning often used in AliPhysics for debugging
  /// \return vector with bin limits
  std::vector<double> makePtBinning() const
  {
    std::vector<double> result;
    result.reserve(1000);
    double epsilon = 1e-6;
    double valGammaPt = 0;
    for (int i = 0; i < 1000; ++i) { // o2-linter: disable=magic-number (just numbers for binning)
      result.push_back(valGammaPt);
      if (valGammaPt < 1.0 - epsilon) { // o2-linter: disable=magic-number (just numbers for binning)
        valGammaPt += 0.1;
      } else if (valGammaPt < 5 - epsilon) { // o2-linter: disable=magic-number (just numbers for binning)
        valGammaPt += 0.2;
      } else if (valGammaPt < 10 - epsilon) { // o2-linter: disable=magic-number (just numbers for binning)
        valGammaPt += 0.5;
      } else if (valGammaPt < 50 - epsilon) { // o2-linter: disable=magic-number (just numbers for binning)
        valGammaPt += 1;
      } else if (valGammaPt < 100 - epsilon) { // o2-linter: disable=magic-number (just numbers for binning)
        valGammaPt += 5;
      } else {
        break;
      }
    }
    return result;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& context)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<EmcTmMonitor>(context)};
  return workflow;
}
