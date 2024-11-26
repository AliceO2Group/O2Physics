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

#include <climits>
#include <cstdlib>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/HistogramRegistry.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "EMCALBase/Geometry.h"
#include "EMCALCalib/BadChannelMap.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "DataFormatsEMCAL/Cell.h"
#include "DataFormatsEMCAL/Constants.h"
#include "DataFormatsEMCAL/AnalysisCluster.h"

#include "CommonDataFormat/InteractionRecord.h"

// \struct TrackMatchingMonitor
/// \brief Simple monitoring task for EMCal clusters
/// \author Marvin Hemmer <marvin.hemmer@cern.ch>
/// \since 24.02.2023
///
/// This task is meant to be used for monitoring the matching between global tracks and EMCal clusters
/// properties, such as:
/// - cluster energy over track momentum
/// - difference in eta
/// - difference in phi
/// Simple event selection using the flag doEventSel is provided, which selects INT7 events if set to 1
/// For pilot beam data, instead of relying on the event selection, one can veto specific BC IDS using the flag
/// fDoVetoBCID.
using namespace o2::framework;
using namespace o2::framework::expressions;
using collisionEvSelIt = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>::iterator;
using bcEvSelIt = o2::soa::Join<o2::aod::BCs, o2::aod::BcSels>::iterator;
using selectedClusters = o2::soa::Filtered<o2::aod::EMCALClusters>;
using selectedAmbiguousClusters = o2::soa::Filtered<o2::aod::EMCALAmbiguousClusters>;
using tracksPID = o2::soa::Join<o2::aod::pidTPCFullEl, o2::aod::pidTPCFullPi, o2::aod::FullTracks, o2::aod::TrackSelection>;
struct TrackMatchingMonitor {
  HistogramRegistry mHistManager{"TrackMatchingMonitorHistograms", {}, OutputObjHandlingPolicy::AnalysisObject};
  o2::emcal::Geometry* mGeometry = nullptr;

  Preslice<o2::aod::EMCALClusterCells> perCluster = o2::aod::emcalclustercell::emcalclusterId;
  Preslice<o2::aod::EMCALAmbiguousClusterCells> perClusterAmb = o2::aod::emcalclustercell::emcalambiguousclusterId;
  Preslice<o2::aod::EMCALMatchedTracks> perClusterMatchedTracks = o2::aod::emcalclustercell::emcalclusterId;
  // configurable parameters
  // TODO adapt mDoEventSel switch to also allow selection of other triggers (e.g. EMC7)
  Configurable<bool> mDoEventSel{"doEventSel", 0, "demand kINT7"};
  Configurable<std::string> mVetoBCID{"vetoBCID", "", "BC ID(s) to be excluded, this should be used as an alternative to the event selection"};
  Configurable<std::string> mSelectBCID{"selectBCID", "all", "BC ID(s) to be included, this should be used as an alternative to the event selection"};
  Configurable<double> mVertexCut{"vertexCut", -1, "apply z-vertex cut with value in cm"};
  Configurable<int> mClusterDefinition{"clusterDefinition", 10, "cluster definition to be selected, e.g. 10=kV3Default"};
  ConfigurableAxis mClusterTimeBinning{"clustertime-binning", {1500, -600, 900}, ""};
  Configurable<bool> hasPropagatedTracks{"hasPropagatedTracks", false, "temporary flag, only set to true when running over data which has the tracks propagated to EMCal/PHOS!"};
  Configurable<bool> usePionRejection{"usePionRejection", false, "demand pion rection for electron signal with TPC PID"};
  Configurable<std::vector<float>> tpcNsigmaElectron{"tpcNsigmaElectron", {-1., +3.}, "TPC PID NSigma range for electron signal (first <= NSigma <= second)"};
  Configurable<std::vector<float>> tpcNsigmaBack{"tpcNsigmaBack", {-10., -4.}, "TPC PID NSigma range for electron background (first <= NSigma <= second)"};
  Configurable<std::vector<float>> tpcNsigmaPion{"tpcNsigmaPion", {-3., +3.}, "TPC PID NSigma range for pions for rejection if usePionRejection is enabled (first <= NSigma <= second)"};
  Configurable<float> minTime{"minTime", -25., "Minimum cluster time for time cut"};
  Configurable<float> maxTime{"maxTime", +20., "Maximum cluster time for time cut"};
  Configurable<float> minM02{"minM02", 0.1, "Minimum M02 for M02 cut"};
  Configurable<float> maxM02{"maxM02", 0.9, "Maximum M02 for M02 cut"};
  Configurable<float> maxM02HighPt{"maxM02HighPt", 0.6, "Maximum M02 for M02 cut for high pT"};
  Configurable<float> M02highPt{"M02highPt", 15., "pT threshold for maxM02HighPt cut. Set to negative value to disable it!"};
  Configurable<float> minDEta{"minDEta", 0.01, "Minimum dEta between track and cluster"};
  Configurable<float> minDPhi{"minDPhi", 0.01, "Minimum dPhi between track and cluster"};
  Configurable<std::vector<float>> eOverPRange{"eOverPRange", {0.9, 1.2}, "E/p range where one would search for electrons (first <= E/p <= second)"};
  Configurable<bool> hasTRD{"hasTRD", false, "demand track to have a hit in TRD."};

  std::vector<int> mVetoBCIDs;
  std::vector<int> mSelectBCIDs;

  /// \brief Create output histograms and initialize geometry
  void init(InitContext const&)
  {
    // create histograms
    using o2HistType = HistType;
    using o2Axis = AxisSpec;

    // load geometry just in case we need it
    mGeometry = o2::emcal::Geometry::GetInstanceFromRunNumber(300000);

    // create common axes
    LOG(info) << "Creating histograms";
    const o2Axis bcAxis{3501, -0.5, 3500.5};
    const o2Axis energyAxis{makeEnergyBinningAliPhysics(), "E_{clus} (GeV)"};
    const o2Axis amplitudeAxisLarge{1000, 0., 100., "amplitudeLarge", "Amplitude (GeV)"};
    const o2Axis dEtaAxis{100, -1.f * minDEta, minDEta, "d#it{#eta}"};
    const o2Axis dPhiAxis{100, -1.f * minDPhi, minDPhi, "d#it{#varphi} (rad)"};
    const o2Axis dRAxis{150, 0.0, 0.015, "d#it{R}"};
    const o2Axis eoverpAxis{500, 0, 10, "#it{E}_{cluster}/#it{p}_{track}"};
    const o2Axis nSigmaAxis{400, -10., +30., "N#sigma"};
    const o2Axis trackptAxis{makePtBinning(), "#it{p}_{T,track}"};
    const o2Axis trackpAxis{200, 0, 100, "#it{p}_{track}"};
    const o2Axis clusterptAxis{makePtBinning(), "#it{p}_{T}"};
    const o2Axis etaAxis{160, -0.8, 0.8, "#eta"};
    const o2Axis phiAxis{72, 0, 2 * 3.14159, "#varphi (rad)"};
    const o2Axis smAxis{20, -0.5, 19.5, "SM"};
    o2Axis timeAxis{mClusterTimeBinning, "t_{cl} (ns)"};

    int MaxMatched = 20; // maximum number of matched tracks, hardcoded in emcalCorrectionTask.cxx!
    const o2Axis nmatchedtrack{MaxMatched, -0.5, MaxMatched + 0.5};

    // event properties
    mHistManager.add("eventsAll", "Number of events", o2HistType::kTH1F, {{1, 0.5, 1.5}});
    mHistManager.add("eventsSelected", "Number of events", o2HistType::kTH1F, {{1, 0.5, 1.5}});
    mHistManager.add("eventVertexZAll", "z-vertex of event (all events)", o2HistType::kTH1F, {{200, -20, 20}});
    mHistManager.add("eventVertexZSelected", "z-vertex of event (selected events)", o2HistType::kTH1F, {{200, -20, 20}});

    // cluster properties (matched clusters)
    mHistManager.add("TrackEtaPhi", "#eta vs #varphi of all selected tracks", o2HistType::kTH2F, {etaAxis, phiAxis});                                                                             // eta vs phi of all selected tracks
    mHistManager.add("TrackEtaPhi_Neg", "#eta vs #varphi of all selected negative tracks", o2HistType::kTH2F, {etaAxis, phiAxis});                                                                // eta vs phi of all selected negative tracks
    mHistManager.add("TrackEtaPhi_Pos", "#eta vs #varphi of all selected positive tracks", o2HistType::kTH2F, {etaAxis, phiAxis});                                                                // eta vs phi of all selected positive tracks
    mHistManager.add("clusterEMatched", "Energy of cluster (with match)", o2HistType::kTH1F, {energyAxis});                                                                                       // energy of matched clusters
    mHistManager.add("MatchedTrackEtaPhi_BeforeCut", "#eta vs #varphi of all selected matched tracks before d#eta and #dvarphi cut", o2HistType::kTH2F, {etaAxis, phiAxis});                      // eta vs phi of all selected matched tracks before dEta and dPhi cut
    mHistManager.add("MatchedTrackEtaPhi_Neg_BeforeCut", "#eta vs #varphi of all selected negative matched tracks before d#eta and #dvarphi cut", o2HistType::kTH2F, {etaAxis, phiAxis});         // eta vs phi of all selected negative matched tracks before dEta and dPhi cut
    mHistManager.add("MatchedTrackEtaPhi_Pos_BeforeCut", "#eta vs #varphi of all selected positive matched tracks before d#eta and #dvarphi cut", o2HistType::kTH2F, {etaAxis, phiAxis});         // eta vs phi of all selected positive matched tracks before dEta and dPhi cut
    mHistManager.add("MatchedTrackEtaPhi", "#eta vs #varphi of all selected matched tracks", o2HistType::kTH2F, {etaAxis, phiAxis});                                                              // eta vs phi of all selected matched tracks
    mHistManager.add("MatchedTrackEtaPhi_Neg", "#eta vs #varphi of all selected negative matched tracks", o2HistType::kTH2F, {etaAxis, phiAxis});                                                 // eta vs phi of all selected negative matched tracks
    mHistManager.add("MatchedTrackEtaPhi_Pos", "#eta vs #varphi of all selected positive matched tracks", o2HistType::kTH2F, {etaAxis, phiAxis});                                                 // eta vs phi of all selected positive matched tracks
    mHistManager.add("clusterTM_dEtadPhi", "cluster trackmatching dEta/dPhi", o2HistType::kTH3F, {dEtaAxis, dPhiAxis, smAxis});                                                                   // dEta dPhi per SM
    mHistManager.add("clusterTM_dEtadPhi_ASide", "cluster trackmatching in A-Side dEta/dPhi", o2HistType::kTH3F, {dEtaAxis, dPhiAxis, smAxis});                                                   // dEta dPhi per SM in A-Aside
    mHistManager.add("clusterTM_dEtadPhi_CSide", "cluster trackmatching in C-Side tracks dEta/dPhi", o2HistType::kTH3F, {dEtaAxis, dPhiAxis, smAxis});                                            // dEta dPhi per SM in C-Side
    mHistManager.add("clusterTM_PosdEtadPhi", "cluster trackmatching positive tracks dEta/dPhi", o2HistType::kTH3F, {dEtaAxis, dPhiAxis, smAxis});                                                // dEta dPhi per SM positive track
    mHistManager.add("clusterTM_NegdEtadPhi", "cluster trackmatching negative tracks dEta/dPhi", o2HistType::kTH3F, {dEtaAxis, dPhiAxis, smAxis});                                                // dEta dPhi per SM negative track
    mHistManager.add("clusterTM_PosdEtadPhi_Pl0_75", "cluster trackmatching positive tracks, p < 0.75 dEta/dPhi", o2HistType::kTH3F, {dEtaAxis, dPhiAxis, smAxis});                               // dEta dPhi per SM positive track with p < 0.75 GeV/c
    mHistManager.add("clusterTM_NegdEtadPhi_Pl0_75", "cluster trackmatching negative tracks, p < 0.75 dEta/dPhi", o2HistType::kTH3F, {dEtaAxis, dPhiAxis, smAxis});                               // dEta dPhi per SM negative track with p < 0.75 GeV/c
    mHistManager.add("clusterTM_PosdEtadPhi_0_75leqPl1_25", "cluster trackmatching positive tracks, 0.75 <= p < 1.25 dEta/dPhi", o2HistType::kTH3F, {dEtaAxis, dPhiAxis, smAxis});                // dEta dPhi per SM positive track with 0.75 <= p < 1.25 GeV/c
    mHistManager.add("clusterTM_NegdEtadPhi_0_75leqPl1_25", "cluster trackmatching negative tracks, 0.75 <= p < 1.25 dEta/dPhi", o2HistType::kTH3F, {dEtaAxis, dPhiAxis, smAxis});                // dEta dPhi per SM negative track with 0.75 <= p < 1.25 GeV/c
    mHistManager.add("clusterTM_PosdEtadPhi_Pgeq1_25", "cluster trackmatching positive tracks, p >= 1.25 dEta/dPhi", o2HistType::kTH3F, {dEtaAxis, dPhiAxis, smAxis});                            // dEta dPhi per SM positive track with p >= 1.25 GeV/c
    mHistManager.add("clusterTM_NegdEtadPhi_Pgeq1_25", "cluster trackmatching negative tracks, p >= 1.25 dEta/dPhi", o2HistType::kTH3F, {dEtaAxis, dPhiAxis, smAxis});                            // dEta dPhi per SM negative track with p >= 1.25 GeV/c
    mHistManager.add("clusterTM_dEtaPt", "cluster trackmatching dEta/#it{p}_{T};d#it{#eta};#it{p}_{T} (GeV/#it{c})", o2HistType::kTH3F, {dEtaAxis, clusterptAxis, smAxis});                       // dEta vs pT per SM
    mHistManager.add("clusterTM_PosdPhiPt", "cluster trackmatching positive tracks dPhi/#it{p}_{T}", o2HistType::kTH3F, {dPhiAxis, clusterptAxis, smAxis});                                       // dPhi vs pT per SM positive track
    mHistManager.add("clusterTM_NegdPhiPt", "cluster trackmatching negative tracks dPh/#it{p}_{T}", o2HistType::kTH3F, {dPhiAxis, clusterptAxis, smAxis});                                        // dPhi vs pT per SM negative track
    mHistManager.add("clusterTM_dEtaTN", "cluster trackmatching dEta/TN;d#it{#eta};#it{N}_{matched tracks}", o2HistType::kTH2F, {dEtaAxis, nmatchedtrack});                                       // dEta compared to the Nth closest track
    mHistManager.add("clusterTM_dPhiTN", "cluster trackmatching dPhi/TN;d#it{#varphi} (rad);#it{N}_{matched tracks}", o2HistType::kTH2F, {dPhiAxis, nmatchedtrack});                              // dPhi compared to the Nth closest track
    mHistManager.add("clusterTM_dRTN", "cluster trackmatching dR/TN;d#it{R};#it{N}_{matched tracks}", o2HistType::kTH2F, {dRAxis, nmatchedtrack});                                                // dR compared to the Nth closest track
    mHistManager.add("clusterTM_NTrack", "cluster trackmatching NMatchedTracks", o2HistType::kTH1I, {nmatchedtrack});                                                                             // how many tracks are matched
    mHistManager.add("clusterTM_EoverP_E", "E/p ", o2HistType::kTH3F, {eoverpAxis, energyAxis, nmatchedtrack});                                                                                   // E/p vs p for the Nth closest track
    mHistManager.add("clusterTM_EoverP_Pt", "E/p vs track pT", o2HistType::kTH3F, {eoverpAxis, trackptAxis, nmatchedtrack});                                                                      // E/p vs track pT for the Nth closest track
    mHistManager.add("clusterTM_EvsP", "cluster E/track p", o2HistType::kTH3F, {energyAxis, trackpAxis, nmatchedtrack});                                                                          // E vs p for the Nth closest track
    mHistManager.add("clusterTM_EoverP_ep", "cluster E/electron p", o2HistType::kTH3F, {eoverpAxis, trackptAxis, nmatchedtrack});                                                                 // E over p vs track pT for the Nth closest electron/positron track
    mHistManager.add("clusterTM_EoverP_e", "cluster E/electron p", o2HistType::kTH3F, {eoverpAxis, trackptAxis, nmatchedtrack});                                                                  // E over p vs track pT for the Nth closest electron track
    mHistManager.add("clusterTM_EoverP_p", "cluster E/electron p", o2HistType::kTH3F, {eoverpAxis, trackptAxis, nmatchedtrack});                                                                  // E over p vs track pT for the Nth closest positron track
    mHistManager.add("clusterTM_EoverP_electron_ASide", "cluster E/electron p in A-Side", o2HistType::kTH3F, {eoverpAxis, trackptAxis, nmatchedtrack});                                           // E over p vs track pT for the Nth closest electron/positron track in A-Side
    mHistManager.add("clusterTM_EoverP_electron_CSide", "cluster E/electron p in C-Side", o2HistType::kTH3F, {eoverpAxis, trackptAxis, nmatchedtrack});                                           // E over p vs track pT for the Nth closest electron/positron track in C-Side
    mHistManager.add("clusterTM_NSigma_BeforeCut", "electron NSigma for matched tracks before d#eta and #dvarphi cut", o2HistType::kTH3F, {nSigmaAxis, trackptAxis, nmatchedtrack});              // NSigma_electron vs track pT for the Nth closest track before dEta and dPhi cut
    mHistManager.add("clusterTM_NSigma_neg_BeforeCut", "electron NSigma for matched negative tracks before d#eta and #dvarphi cut", o2HistType::kTH3F, {nSigmaAxis, trackptAxis, nmatchedtrack}); // NSigma_electron vs track pT for the Nth closest negative tracks before dEta and dPhi cut
    mHistManager.add("clusterTM_NSigma_pos_BeforeCut", "electron NSigma for matched positive tracks before d#eta and #dvarphi cut", o2HistType::kTH3F, {nSigmaAxis, trackptAxis, nmatchedtrack}); // NSigma_electron vs track pT for the Nth closest positive tracks before dEta and dPhi cut
    mHistManager.add("clusterTM_NSigma", "electron NSigma for matched track", o2HistType::kTH3F, {nSigmaAxis, trackptAxis, nmatchedtrack});                                                       // NSigma_electron vs track pT for the Nth closest track
    mHistManager.add("clusterTM_NSigma_neg", "electron NSigma for matched negative track", o2HistType::kTH3F, {nSigmaAxis, trackptAxis, nmatchedtrack});                                          // NSigma_electron vs track pT for the Nth closest negative tracks
    mHistManager.add("clusterTM_NSigma_pos", "electron NSigma for matched positive track", o2HistType::kTH3F, {nSigmaAxis, trackptAxis, nmatchedtrack});                                          // NSigma_electron vs track pT for the Nth closest positive tracks
    mHistManager.add("clusterTM_NSigma_cut", "electron NSigma for matched track with cuts", o2HistType::kTH3F, {nSigmaAxis, trackptAxis, nmatchedtrack});                                         // NSigma_electron vs track pT for the Nth closest track with cuts on E/p and cluster cuts
    mHistManager.add("clusterTM_NSigma_e", "NSigma electron", o2HistType::kTH3F, {nSigmaAxis, trackptAxis, nmatchedtrack});                                                                       // NSigma vs track pT for the Nth closest electron track
    mHistManager.add("clusterTM_NSigma_p", "NSigma positron", o2HistType::kTH3F, {nSigmaAxis, trackptAxis, nmatchedtrack});                                                                       // NSigma vs track pT for the Nth closest positron track
    mHistManager.add("clusterTM_NSigma_e_cut", "NSigma electron with E over p cut", o2HistType::kTH3F, {nSigmaAxis, trackptAxis, nmatchedtrack});                                                 // NSigma vs track pT for the Nth closest electron track with E/p cut
    mHistManager.add("clusterTM_NSigma_p_cut", "NSigma positron with E over p cut", o2HistType::kTH3F, {nSigmaAxis, trackptAxis, nmatchedtrack});                                                 // NSigma vs track pT for the Nth closest positron track with E/p cut
    mHistManager.add("TrackTM_NSigma", "electron NSigma for global tracks", o2HistType::kTH2F, {nSigmaAxis, trackptAxis});                                                                        // NSigma_electron vs track pT for global track
    mHistManager.add("TrackTM_NSigma_e", "NSigma e for global negative tracks", o2HistType::kTH2F, {nSigmaAxis, trackptAxis});                                                                    // NSigma vs track pT for negative global track
    mHistManager.add("TrackTM_NSigma_p", "NSigma e for global positive tracks", o2HistType::kTH2F, {nSigmaAxis, trackptAxis});                                                                    // NSigma vs track pT for positive global track
    mHistManager.add("clusterTM_NSigma_electron_ASide", "NSigma electron in A-Side", o2HistType::kTH3F, {nSigmaAxis, trackptAxis, nmatchedtrack});                                                // NSigma vs track pT for the Nth closest electron/positron track in A-Side
    mHistManager.add("clusterTM_NSigma_electron_CSide", "NSigma positron in C-Side", o2HistType::kTH3F, {nSigmaAxis, trackptAxis, nmatchedtrack});                                                // NSigma vs track pT for the Nth closest electron/positron track in C-Side
    mHistManager.add("clusterTM_EoverP_hadron", "cluster E/hadron p", o2HistType::kTH3F, {eoverpAxis, trackpAxis, nmatchedtrack});                                                                // E over p vs track pT for the Nth closest hadron track
    mHistManager.add("clusterTM_EoverP_hn", "cluster E/hadron p", o2HistType::kTH3F, {eoverpAxis, trackpAxis, nmatchedtrack});                                                                    // E over p vs track pT for the Nth closest negative hadron track
    mHistManager.add("clusterTM_EoverP_hp", "cluster E/hadron p", o2HistType::kTH3F, {eoverpAxis, trackpAxis, nmatchedtrack});                                                                    // E over p vs track pT for the Nth closest positive hadron track
    mHistManager.add("clusterTM_EoverP_hadron_ASide", "cluster E/hadron p in A-Side", o2HistType::kTH3F, {eoverpAxis, trackpAxis, nmatchedtrack});                                                // E over p vs track pT for the Nth closest hadron track in A-Side
    mHistManager.add("clusterTM_EoverP_hadron_CSide", "cluster E/hadron p in C-Side", o2HistType::kTH3F, {eoverpAxis, trackpAxis, nmatchedtrack});                                                // E over p vs track pT for the Nth closest hadron track in C-Side

    if (mVetoBCID->length()) {
      std::stringstream parser(mVetoBCID.value);
      std::string token;
      int bcid;
      while (std::getline(parser, token, ',')) {
        bcid = std::stoi(token);
        LOG(info) << "Veto BCID " << bcid;
        mVetoBCIDs.push_back(bcid);
      }
    }
    if (mSelectBCID.value != "all") {
      std::stringstream parser(mSelectBCID.value);
      std::string token;
      int bcid;
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
  Filter clusterDefinitionSelection = (o2::aod::emcalcluster::definition == mClusterDefinition) && (o2::aod::emcalcluster::time >= minTime) && (o2::aod::emcalcluster::time <= maxTime) && (o2::aod::emcalcluster::m02 > minM02) && (o2::aod::emcalcluster::m02 < maxM02);

  /// \brief Process EMCAL clusters that are matched to a collisions
  void processCollisions(collisionEvSelIt const& theCollision, selectedClusters const& clusters, o2::aod::EMCALClusterCells const&, o2::aod::Calos const&, o2::aod::EMCALMatchedTracks const& matchedtracks, tracksPID const& alltracks)
  {
    mHistManager.fill(HIST("eventsAll"), 1);

    // do event selection if mDoEventSel is specified
    // currently the event selection is hard coded to kINT7
    // but other selections are possible that are defined in TriggerAliases.h
    bool isSelected = true;
    if (mDoEventSel) {
      if (theCollision.bc().runNumber() < 300000) {
        if (!theCollision.alias_bit(kINT7)) {
          isSelected = false;
        }
      } else {
        if (!theCollision.alias_bit(kTVXinEMC)) {
          isSelected = false;
        }
      }
    }
    if (!isSelected) {
      LOG(debug) << "Event not selected because it is not kINT7 or does not have EMCAL in readout, skipping";
      return;
    }
    mHistManager.fill(HIST("eventVertexZAll"), theCollision.posZ());
    if (mVertexCut > 0 && TMath::Abs(theCollision.posZ()) > mVertexCut) {
      LOG(debug) << "Event not selected because of z-vertex cut z= " << theCollision.posZ() << " > " << mVertexCut << " cm, skipping";
      return;
    }
    mHistManager.fill(HIST("eventsSelected"), 1);
    mHistManager.fill(HIST("eventVertexZSelected"), theCollision.posZ());

    for (const auto& alltrack : alltracks) {
      double trackEta, trackPhi;
      if (hasPropagatedTracks) { // only temporarily while not every data has the tracks propagated to EMCal/PHOS
        trackEta = alltrack.trackEtaEmcal();
        trackPhi = alltrack.trackPhiEmcal();
      } else {
        trackEta = alltrack.eta();
        trackPhi = alltrack.phi();
      }
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
      double dEta, dPhi, pT, abs_p, trackEta, trackPhi, NSigmaEl;
      int supermoduleID;
      try {
        supermoduleID = mGeometry->SuperModuleNumberFromEtaPhi(cluster.eta(), cluster.phi());
      } catch (o2::emcal::InvalidPositionException& e) {
        // Imprecision of the position at the sector boundaries, mostly due to
        // vertex imprecision. Skip these clusters for the now.
        continue;
      }

      pT = cluster.energy() / cosh(cluster.eta());
      if (M02highPt > 0. && cluster.m02() >= maxM02HighPt && pT >= M02highPt) { // high pT M02 cut
        continue;
      }
      auto tracksofcluster = matchedtracks.sliceBy(perClusterMatchedTracks, cluster.globalIndex());
      int t = 0;
      for (const auto& match : tracksofcluster) {
        NSigmaEl = match.track_as<tracksPID>().tpcNSigmaEl();
        // exmple of how to access any property of the matched tracks (tracks are sorted by how close they are to cluster)
        LOG(debug) << "Pt of match" << match.track_as<tracksPID>().pt();
        abs_p = abs(match.track_as<tracksPID>().p());
        // only consider closest match
        if (hasPropagatedTracks) { // only temporarily while not every data has the tracks propagated to EMCal/PHOS
          trackEta = match.track_as<tracksPID>().trackEtaEmcal();
          trackPhi = match.track_as<tracksPID>().trackPhiEmcal();
        } else {
          trackEta = match.track_as<tracksPID>().eta();
          trackPhi = match.track_as<tracksPID>().phi();
        }
        dPhi = trackPhi - cluster.phi();
        dEta = trackEta - cluster.eta();
        mHistManager.fill(HIST("MatchedTrackEtaPhi_BeforeCut"), trackEta, trackPhi);
        mHistManager.fill(HIST("clusterTM_NSigma_BeforeCut"), NSigmaEl, match.track_as<tracksPID>().pt(), t);
        if (match.track_as<tracksPID>().sign() == -1) {
          mHistManager.fill(HIST("MatchedTrackEtaPhi_Neg_BeforeCut"), trackEta, trackPhi);
          mHistManager.fill(HIST("clusterTM_NSigma_neg_BeforeCut"), NSigmaEl, match.track_as<tracksPID>().pt(), t);
        } else if (match.track_as<tracksPID>().sign() == 1) {
          mHistManager.fill(HIST("MatchedTrackEtaPhi_Pos_BeforeCut"), trackEta, trackPhi);
          mHistManager.fill(HIST("clusterTM_NSigma_pos_BeforeCut"), NSigmaEl, match.track_as<tracksPID>().pt(), t);
        }
        if (fabs(dEta) >= minDEta || fabs(dPhi) >= minDPhi) { // dEta and dPhi cut
          continue;
        }
        if (hasTRD && !(match.track_as<tracksPID>().hasTRD())) { // request TRD hit cut
          continue;
        }
        // only fill these for the first matched track:
        if (t == 0) {
          mHistManager.fill(HIST("clusterTM_dEtadPhi"), dEta, dPhi, supermoduleID);
          mHistManager.fill(HIST("clusterTM_dEtaPt"), dEta, pT, supermoduleID);
        }
        double eOverP = cluster.energy() / abs_p;
        mHistManager.fill(HIST("clusterTM_dEtaTN"), dEta, t);
        mHistManager.fill(HIST("clusterTM_dPhiTN"), dPhi, t);
        mHistManager.fill(HIST("clusterTM_dRTN"), std::sqrt(dPhi * dPhi + dEta * dEta), t);
        mHistManager.fill(HIST("clusterTM_EoverP_E"), eOverP, cluster.energy(), t);
        mHistManager.fill(HIST("clusterTM_dEtadPhi"), dEta, dPhi, t);
        mHistManager.fill(HIST("clusterEMatched"), cluster.energy(), t);
        mHistManager.fill(HIST("clusterTM_EvsP"), cluster.energy(), abs_p, t);
        mHistManager.fill(HIST("clusterTM_EoverP_Pt"), eOverP, match.track_as<tracksPID>().pt(), t);
        mHistManager.fill(HIST("clusterTM_NSigma"), NSigmaEl, match.track_as<tracksPID>().pt(), t);
        mHistManager.fill(HIST("MatchedTrackEtaPhi"), trackEta, trackPhi);
        if (eOverPRange->at(0) <= eOverP && eOverP <= eOverPRange->at(1)) {
          mHistManager.fill(HIST("clusterTM_NSigma_cut"), NSigmaEl, match.track_as<tracksPID>().pt(), t);
        }
        // A- and C-side
        if (match.track_as<tracksPID>().eta() > 0.0 && t == 0) {
          mHistManager.fill(HIST("clusterTM_dEtadPhi_ASide"), dEta, dPhi, supermoduleID);
        } else if (match.track_as<tracksPID>().eta() < 0.0 && t == 0) {
          mHistManager.fill(HIST("clusterTM_dEtadPhi_CSide"), dEta, dPhi, supermoduleID);
        }
        // positive and negative tracks seperate, with three different track momentum ranges
        if (match.track_as<tracksPID>().sign() == 1) {
          mHistManager.fill(HIST("clusterTM_NSigma_pos"), NSigmaEl, match.track_as<tracksPID>().pt(), t);
          mHistManager.fill(HIST("MatchedTrackEtaPhi_Pos"), trackEta, trackPhi);
          if (t == 0) {
            mHistManager.fill(HIST("clusterTM_PosdPhiPt"), dPhi, pT, supermoduleID);
            mHistManager.fill(HIST("clusterTM_PosdEtadPhi"), dEta, dPhi, supermoduleID);
            if (abs_p < 0.75) {
              mHistManager.fill(HIST("clusterTM_PosdEtadPhi_Pl0_75"), dEta, dPhi, supermoduleID);
            } else if (abs_p >= 1.25) {
              mHistManager.fill(HIST("clusterTM_PosdEtadPhi_Pgeq1_25"), dEta, dPhi, supermoduleID);
            } else {
              mHistManager.fill(HIST("clusterTM_PosdEtadPhi_0_75leqPl1_25"), dEta, dPhi, supermoduleID);
            }
          }
        } else if (match.track_as<tracksPID>().sign() == -1) {
          mHistManager.fill(HIST("clusterTM_NSigma_neg"), NSigmaEl, match.track_as<tracksPID>().pt(), t);
          mHistManager.fill(HIST("MatchedTrackEtaPhi_Neg"), trackEta, trackPhi);
          if (t == 0) {
            mHistManager.fill(HIST("clusterTM_NegdPhiPt"), dPhi, pT, supermoduleID);
            mHistManager.fill(HIST("clusterTM_NegdEtadPhi"), dEta, dPhi, supermoduleID);
            if (abs_p < 0.75) {
              mHistManager.fill(HIST("clusterTM_NegdEtadPhi_Pl0_75"), dEta, dPhi, supermoduleID);
            } else if (abs_p >= 1.25) {
              mHistManager.fill(HIST("clusterTM_NegdEtadPhi_Pgeq1_25"), dEta, dPhi, supermoduleID);
            } else {
              mHistManager.fill(HIST("clusterTM_NegdEtadPhi_0_75leqPl1_25"), dEta, dPhi, supermoduleID);
            }
          }
        }
        if (tpcNsigmaElectron->at(0) <= NSigmaEl && NSigmaEl <= tpcNsigmaElectron->at(1)) {                                                                                   // E/p for e+/e-
          if (usePionRejection && (tpcNsigmaPion->at(0) <= match.track_as<tracksPID>().tpcNSigmaPi() || match.track_as<tracksPID>().tpcNSigmaPi() <= tpcNsigmaPion->at(1))) { // with pion rejection
            mHistManager.fill(HIST("clusterTM_EoverP_ep"), eOverP, match.track_as<tracksPID>().pt(), t);
            if (match.track_as<tracksPID>().eta() >= 0.) {
              mHistManager.fill(HIST("clusterTM_EoverP_electron_ASide"), eOverP, match.track_as<tracksPID>().pt(), t);
              mHistManager.fill(HIST("clusterTM_NSigma_electron_ASide"), NSigmaEl, match.track_as<tracksPID>().pt(), t);
            } else {
              mHistManager.fill(HIST("clusterTM_EoverP_electron_CSide"), eOverP, match.track_as<tracksPID>().pt(), t);
              mHistManager.fill(HIST("clusterTM_NSigma_electron_CSide"), NSigmaEl, match.track_as<tracksPID>().pt(), t);
            }
            if (match.track_as<tracksPID>().sign() == -1) {
              mHistManager.fill(HIST("clusterTM_EoverP_e"), eOverP, match.track_as<tracksPID>().pt(), t);
              mHistManager.fill(HIST("clusterTM_NSigma_e"), NSigmaEl, match.track_as<tracksPID>().pt(), t);
              if (eOverPRange->at(0) <= eOverP && eOverP <= eOverPRange->at(1)) {
                mHistManager.fill(HIST("clusterTM_NSigma_e_cut"), NSigmaEl, match.track_as<tracksPID>().pt(), t);
              }
            } else if (match.track_as<tracksPID>().sign() == +1) {
              mHistManager.fill(HIST("clusterTM_EoverP_p"), eOverP, match.track_as<tracksPID>().pt(), t);
              mHistManager.fill(HIST("clusterTM_NSigma_p"), NSigmaEl, match.track_as<tracksPID>().pt(), t);
              if (eOverPRange->at(0) <= eOverP && eOverP <= eOverPRange->at(1)) {
                mHistManager.fill(HIST("clusterTM_NSigma_p_cut"), NSigmaEl, match.track_as<tracksPID>().pt(), t);
              }
            }
          } else { // without pion rejection
            mHistManager.fill(HIST("clusterTM_EoverP_ep"), eOverP, match.track_as<tracksPID>().pt(), t);
            if (match.track_as<tracksPID>().eta() >= 0.) {
              mHistManager.fill(HIST("clusterTM_EoverP_electron_ASide"), eOverP, match.track_as<tracksPID>().pt(), t);
              mHistManager.fill(HIST("clusterTM_NSigma_electron_ASide"), NSigmaEl, match.track_as<tracksPID>().pt(), t);
            } else {
              mHistManager.fill(HIST("clusterTM_EoverP_electron_CSide"), eOverP, match.track_as<tracksPID>().pt(), t);
              mHistManager.fill(HIST("clusterTM_NSigma_electron_CSide"), NSigmaEl, match.track_as<tracksPID>().pt(), t);
            }
            if (match.track_as<tracksPID>().sign() == -1) {
              mHistManager.fill(HIST("clusterTM_EoverP_e"), eOverP, match.track_as<tracksPID>().pt(), t);
              mHistManager.fill(HIST("clusterTM_NSigma_e"), NSigmaEl, match.track_as<tracksPID>().pt(), t);
              if (eOverPRange->at(0) <= eOverP && eOverP <= eOverPRange->at(1)) {
                mHistManager.fill(HIST("clusterTM_NSigma_e_cut"), NSigmaEl, match.track_as<tracksPID>().pt(), t);
              }
            } else if (match.track_as<tracksPID>().sign() == +1) {
              mHistManager.fill(HIST("clusterTM_EoverP_p"), eOverP, match.track_as<tracksPID>().pt(), t);
              mHistManager.fill(HIST("clusterTM_NSigma_p"), NSigmaEl, match.track_as<tracksPID>().pt(), t);
              if (eOverPRange->at(0) <= eOverP && eOverP <= eOverPRange->at(1)) {
                mHistManager.fill(HIST("clusterTM_NSigma_p_cut"), NSigmaEl, match.track_as<tracksPID>().pt(), t);
              }
            }
          }
        } else if (tpcNsigmaBack->at(0) <= NSigmaEl && NSigmaEl <= tpcNsigmaBack->at(1)) { // E/p for hadrons / background
          mHistManager.fill(HIST("clusterTM_EoverP_hadron"), eOverP, match.track_as<tracksPID>().pt(), t);
          if (match.track_as<tracksPID>().eta() >= 0.) {
            mHistManager.fill(HIST("clusterTM_EoverP_hadron_ASide"), eOverP, match.track_as<tracksPID>().pt(), t);
          } else {
            mHistManager.fill(HIST("clusterTM_EoverP_hadron_CSide"), eOverP, match.track_as<tracksPID>().pt(), t);
          }
          if (match.track_as<tracksPID>().sign() == -1) {
            mHistManager.fill(HIST("clusterTM_EoverP_hn"), eOverP, match.track_as<tracksPID>().pt(), t);
          } else if (match.track_as<tracksPID>().sign() == +1) {
            mHistManager.fill(HIST("clusterTM_EoverP_hp"), eOverP, match.track_as<tracksPID>().pt(), t);
          }
        }
        t++;
      }
      mHistManager.fill(HIST("clusterTM_NTrack"), t);
    }
  }
  PROCESS_SWITCH(TrackMatchingMonitor, processCollisions, "Process clusters from collision", true);

  /// \brief Create binning for cluster energy axis (variable bin size)
  /// \return vector with bin limits
  std::vector<double>
    makeEnergyBinning() const
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
    Int_t nBinsClusterE = 235;
    for (Int_t i = 0; i < nBinsClusterE + 1; i++) {
      if (i < 1)
        result.emplace_back(0.3 * i);
      else if (i < 55)
        result.emplace_back(0.3 + 0.05 * (i - 1));
      else if (i < 105)
        result.emplace_back(3. + 0.1 * (i - 55));
      else if (i < 140)
        result.emplace_back(8. + 0.2 * (i - 105));
      else if (i < 170)
        result.emplace_back(15. + 0.5 * (i - 140));
      else if (i < 190)
        result.emplace_back(30. + 1.0 * (i - 170));
      else if (i < 215)
        result.emplace_back(50. + 2.0 * (i - 190));
      else if (i < 235)
        result.emplace_back(100. + 5.0 * (i - 215));
      else if (i < 245)
        result.emplace_back(200. + 10.0 * (i - 235));
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
    for (int i = 0; i < 1000; ++i) {
      result.push_back(valGammaPt);
      if (valGammaPt < 1.0 - epsilon)
        valGammaPt += 0.1;
      else if (valGammaPt < 5 - epsilon)
        valGammaPt += 0.2;
      else if (valGammaPt < 10 - epsilon)
        valGammaPt += 0.5;
      else if (valGammaPt < 50 - epsilon)
        valGammaPt += 1;
      else if (valGammaPt < 100 - epsilon)
        valGammaPt += 5;
      else
        break;
    }
    return result;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<TrackMatchingMonitor>(cfgc, TaskName{"emc-tmmonitor"})};
  return workflow;
}
