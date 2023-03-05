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
using tracksPID = o2::soa::Join<o2::aod::pidTPCFullEl, o2::aod::pidTPCFullPi, o2::aod::FullTracks>;
struct TrackMatchingMonitor {
  HistogramRegistry mHistManager{"TrackMatchingMonitorHistograms", {}, OutputObjHandlingPolicy::QAObject};
  o2::emcal::Geometry* mGeometry = nullptr;

  Preslice<o2::aod::EMCALClusterCells> perCluster = o2::aod::emcalclustercell::emcalclusterId;
  Preslice<o2::aod::EMCALAmbiguousClusterCells> perClusterAmb = o2::aod::emcalclustercell::emcalclusterId;
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
  Configurable<std::vector<float>> tpcNsigmaElectron{"tpcNsigmaElectron", {-1., +3.}, "TPC PID NSigma range for electron signal"};
  Configurable<std::vector<float>> tpcNsigmaBack{"tpcNsigmaBack", {-10., -4.}, "TPC PID NSigma range for electron background"};
  Configurable<std::vector<float>> tpcNsigmaPion{"tpcNsigmaPion", {-3., +3.}, "TPC PID NSigma range for pions for rejection if usePionRejection is enabled"};

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
    const o2Axis dEtaAxis{100, -0.4, 0.4, "d#it{#eta}"};
    const o2Axis dPhiAxis{100, -0.4, 0.4, "d#it{#varphi} (rad)"};
    const o2Axis dRAxis{100, 0.0, 0.4, "d#it{R}"};
    const o2Axis eoverpAxis{500, 0, 10, "#it{E}_{cluster}/#it{p}_{track}"};
    const o2Axis trackptAxis{200, 0, 100, "#it{p}_{T,track}"};
    const o2Axis trackpAxis{200, 0, 100, "#it{p}_{track}"};
    const o2Axis clusterptAxis{200, 0, 100, "#it{p}_{T}"};
    o2Axis timeAxis{mClusterTimeBinning, "t_{cl} (ns)"};

    int MaxMatched = 20; // maximum number of matched tracks, hardcoded in emcalCorrectionTask.cxx!
    const o2Axis nmatchedtrack{MaxMatched, -0.5, MaxMatched + 0.5};

    // event properties
    mHistManager.add("eventsAll", "Number of events", o2HistType::kTH1F, {{1, 0.5, 1.5}});
    mHistManager.add("eventsSelected", "Number of events", o2HistType::kTH1F, {{1, 0.5, 1.5}});
    mHistManager.add("eventBCAll", "Bunch crossing ID of event (all events)", o2HistType::kTH1F, {bcAxis});
    mHistManager.add("eventBCSelected", "Bunch crossing ID of event (selected events)", o2HistType::kTH1F, {bcAxis});
    mHistManager.add("eventVertexZAll", "z-vertex of event (all events)", o2HistType::kTH1F, {{200, -20, 20}});
    mHistManager.add("eventVertexZSelected", "z-vertex of event (selected events)", o2HistType::kTH1F, {{200, -20, 20}});

    // cluster properties (matched clusters)
    mHistManager.add("clusterAmpFractionLeadingCell", "Fraction of energy in leading cell", o2HistType::kTH1F, {{100, 0, 1}});
    mHistManager.add("clusterEMatched", "Energy of cluster (with match)", o2HistType::kTH1F, {energyAxis});
    mHistManager.add("clusterTM_dEtadPhi", "cluster trackmatching dEta/dPhi", o2HistType::kTH3F, {dEtaAxis, dPhiAxis, nmatchedtrack});                                                    // dEta dPhi of the Nth clostest track
    mHistManager.add("clusterTM_dEtadPhi_ASide", "cluster trackmatching in A-Side dEta/dPhi", o2HistType::kTH3F, {dEtaAxis, dPhiAxis, nmatchedtrack});                                    // dEta dPhi of the Nth clostest track in A-Aside
    mHistManager.add("clusterTM_dEtadPhi_CSide", "cluster trackmatching in C-Side tracks dEta/dPhi", o2HistType::kTH3F, {dEtaAxis, dPhiAxis, nmatchedtrack});                             // dEta dPhi of the Nth clostest track in C-Side
    mHistManager.add("clusterTM_PosdEtadPhi", "cluster trackmatching positive tracks dEta/dPhi", o2HistType::kTH3F, {dEtaAxis, dPhiAxis, nmatchedtrack});                                 // dEta dPhi of the Nth clostest positive track
    mHistManager.add("clusterTM_NegdEtadPhi", "cluster trackmatching negative tracks dEta/dPhi", o2HistType::kTH3F, {dEtaAxis, dPhiAxis, nmatchedtrack});                                 // dEta dPhi of the Nth clostest negative track
    mHistManager.add("clusterTM_PosdEtadPhi_Pl0_75", "cluster trackmatching positive tracks, p < 0.75 dEta/dPhi", o2HistType::kTH3F, {dEtaAxis, dPhiAxis, nmatchedtrack});                // dEta dPhi of the Nth clostest positive track with p < 0.75 GeV/c
    mHistManager.add("clusterTM_NegdEtadPhi_Pl0_75", "cluster trackmatching negative tracks, p < 0.75 dEta/dPhi", o2HistType::kTH3F, {dEtaAxis, dPhiAxis, nmatchedtrack});                // dEta dPhi of the Nth clostest negative track with p < 0.75 GeV/c
    mHistManager.add("clusterTM_PosdEtadPhi_0_75leqPl1_25", "cluster trackmatching positive tracks, 0.75 <= p < 1.25 dEta/dPhi", o2HistType::kTH3F, {dEtaAxis, dPhiAxis, nmatchedtrack}); // dEta dPhi of the Nth clostest positive track with 0.75 <= p < 1.25 GeV/c
    mHistManager.add("clusterTM_NegdEtadPhi_0_75leqPl1_25", "cluster trackmatching negative tracks, 0.75 <= p < 1.25 dEta/dPhi", o2HistType::kTH3F, {dEtaAxis, dPhiAxis, nmatchedtrack}); // dEta dPhi of the Nth clostest negative track with 0.75 <= p < 1.25 GeV/c
    mHistManager.add("clusterTM_PosdEtadPhi_Pgeq1_25", "cluster trackmatching positive tracks, p >= 1.25 dEta/dPhi", o2HistType::kTH3F, {dEtaAxis, dPhiAxis, nmatchedtrack});             // dEta dPhi of the Nth clostest positive track with p >= 1.25 GeV/c
    mHistManager.add("clusterTM_NegdEtadPhi_Pgeq1_25", "cluster trackmatching negative tracks, p >= 1.25 dEta/dPhi", o2HistType::kTH3F, {dEtaAxis, dPhiAxis, nmatchedtrack});             // dEta dPhi of the Nth clostest negative track with p >= 1.25 GeV/c
    mHistManager.add("clusterTM_dEtaPt", "cluster trackmatching dEta/#it{p}_{T};d#it{#eta};#it{p}_{T} (GeV/#it{c})", o2HistType::kTH3F, {dEtaAxis, clusterptAxis, nmatchedtrack});        // dEta vs pT of the Nth clostest track
    mHistManager.add("clusterTM_PosdPhiPt", "cluster trackmatching positive tracks dPhi/#it{p}_{T}", o2HistType::kTH3F, {dPhiAxis, clusterptAxis, nmatchedtrack});                        // dPhi vs pT of the Nth clostest positive track
    mHistManager.add("clusterTM_NegdPhiPt", "cluster trackmatching negative tracks dPh/#it{p}_{T}", o2HistType::kTH3F, {dPhiAxis, clusterptAxis, nmatchedtrack});                         // dPhi vs pT of the Nth clostest negative track
    mHistManager.add("clusterTM_dEtaTN", "cluster trackmatching dEta/TN;d#it{#eta};#it{N}_{matched tracks}", o2HistType::kTH2F, {dEtaAxis, nmatchedtrack});                               // dEta compared to the Nth closest track
    mHistManager.add("clusterTM_dPhiTN", "cluster trackmatching dPhi/TN;d#it{#varphi} (rad);#it{N}_{matched tracks}", o2HistType::kTH2F, {dPhiAxis, nmatchedtrack});                      // dPhi compared to the Nth closest track
    mHistManager.add("clusterTM_dRTN", "cluster trackmatching dR/TN;d#it{R};#it{N}_{matched tracks}", o2HistType::kTH2F, {dRAxis, nmatchedtrack});                                        // dR compared to the Nth closest track
    mHistManager.add("clusterTM_NTrack", "cluster trackmatching NMatchedTracks", o2HistType::kTH1I, {nmatchedtrack});                                                                     // how many tracks are matched
    mHistManager.add("clusterTM_dEtaTNAli", "cluster trackmatching dEta/TN;d#it{#eta};#it{N}_{matched tracks}", o2HistType::kTH2F, {dEtaAxis, nmatchedtrack});                            // dEta compared to the Nth closest track with cuts from latest Pi0 Run2 analysis
    mHistManager.add("clusterTM_dPhiTNAli", "cluster trackmatching dPhi/TN;d#it{#varphi} (rad);#it{N}_{matched tracks}", o2HistType::kTH2F, {dPhiAxis, nmatchedtrack});                   // dPhi compared to the Nth closest track with cuts from latest Pi0 Run2 analysis
    mHistManager.add("clusterTM_dRTNAli", "cluster trackmatching dR/TN;d#it{R};#it{N}_{matched tracks}", o2HistType::kTH2F, {dRAxis, nmatchedtrack});                                     // dR compared to the Nth closest track with cuts from latest Pi0 Run2 analysis
    mHistManager.add("clusterTM_NTrackAli", "cluster trackmatching NMatchedTracks", o2HistType::kTH1I, {nmatchedtrack});                                                                  // how many tracks are matched with cuts from latest Pi0 Run2 analysis
    mHistManager.add("clusterTM_EoverP_E", "E/p ", o2HistType::kTH3F, {eoverpAxis, energyAxis, nmatchedtrack});                                                                           // E/p vs p for the Nth closest track
    mHistManager.add("clusterTM_EoverP_Pt", "E/p vs track pT", o2HistType::kTH3F, {eoverpAxis, trackptAxis, nmatchedtrack});                                                              // E/p vs track pT for the Nth closest track
    mHistManager.add("clusterTM_EvsP", "cluster E/track p", o2HistType::kTH3F, {energyAxis, trackpAxis, nmatchedtrack});                                                                  // E vs p for the Nth closest track
    mHistManager.add("clusterTM_EoverP_electron", "cluster E/electron p", o2HistType::kTH3F, {eoverpAxis, trackpAxis, nmatchedtrack});                                                    // E over p vs track pT for the Nth closest electron/positron track
    mHistManager.add("clusterTM_EoverP_hadron", "cluster E/hadron p)", o2HistType::kTH3F, {eoverpAxis, trackpAxis, nmatchedtrack});                                                       // E over p vs track pT for the Nth closest hadron track

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
  Filter clusterDefinitionSelection = (o2::aod::emcalcluster::definition == mClusterDefinition);

  /// \brief Process EMCAL clusters that are matched to a collisions
  void processCollisions(collisionEvSelIt const& theCollision, selectedClusters const& clusters, o2::aod::EMCALClusterCells const& emccluscells, o2::aod::Calos const& allcalos, o2::aod::EMCALMatchedTracks const& matchedtracks, tracksPID const& alltrack)
  {
    mHistManager.fill(HIST("eventsAll"), 1);

    // do event selection if mDoEventSel is specified
    // currently the event selection is hard coded to kINT7
    // but other selections are possible that are defined in TriggerAliases.h
    bool isSelected = true;
    if (mDoEventSel) {
      if (theCollision.bc().runNumber() < 300000) {
        if (!theCollision.alias()[kINT7]) {
          isSelected = false;
        }
      } else {
        if (!theCollision.alias()[kTVXinEMC]) {
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

    LOG(debug) << "bunch crossing ID" << theCollision.bcId();
    std::array<int, 20> numberOfClustersSM;
    std::fill(numberOfClustersSM.begin(), numberOfClustersSM.end(), 0);
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

      // example of loop over all cells of current cluster
      // cell.calo() allows to access the cell properties as defined in AnalysisDataModel
      // In this exammple, we loop over all cells and find the cell of maximum energy and plot the fraction
      // it carries of the whole cluster
      auto cellsofcluster = emccluscells.sliceBy(perCluster, cluster.globalIndex());
      double maxamp = 0;
      double ampfraction = 0;
      for (const auto& cell : cellsofcluster) {
        // example how to get any information of the cell associated with cluster
        LOG(debug) << "Cell ID:" << cell.calo().amplitude() << " Time " << cell.calo().time();
        if (cell.calo().amplitude() > maxamp) {
          maxamp = cell.calo().amplitude();
        }
      }
      ampfraction = maxamp / cluster.energy();
      mHistManager.fill(HIST("clusterAmpFractionLeadingCell"), ampfraction);

      // Example of loop over tracks matched to cluster within dR=0.4, where only the
      // 20 most closest tracks are stored. If needed the number of tracks can be later
      // increased in the correction framework. Access to all track properties via match.track()
      // If you want to access tracks from a joined table you can access all track properties via
      // match.track_as<globTracks>() with
      // using globTracks = o2::soa::Join<o2::aod::Tracks, o2::aod::TrackSelection>;
      // In this example the counter t is just used to only look at the closest match
      double dEta, dPhi, pT, abs_p;
      pT = cluster.energy() / cosh(cluster.eta());
      auto tracksofcluster = matchedtracks.sliceBy(perClusterMatchedTracks, cluster.globalIndex());
      int t = 0;
      int tAli = 0;
      for (const auto& match : tracksofcluster) {
        // exmple of how to access any property of the matched tracks (tracks are sorted by how close they are to cluster)
        LOG(debug) << "Pt of match" << match.track_as<tracksPID>().pt();
        abs_p = abs(match.track_as<tracksPID>().p());
        // only consider closest match
        if (hasPropagatedTracks) { // only temporarily while not every data has the tracks propagated to EMCal/PHOS
          dEta = match.track_as<tracksPID>().trackEtaEmcal() - cluster.eta();
          dPhi = match.track_as<tracksPID>().trackPhiEmcal() - cluster.phi();
        } else {
          dEta = match.track_as<tracksPID>().eta() - cluster.eta();
          dPhi = match.track_as<tracksPID>().phi() - cluster.phi();
        }
        mHistManager.fill(HIST("clusterTM_dEtaTN"), dEta, t);
        mHistManager.fill(HIST("clusterTM_dPhiTN"), dPhi, t);
        mHistManager.fill(HIST("clusterTM_dRTN"), std::sqrt(dPhi * dPhi + dEta * dEta), t);
        mHistManager.fill(HIST("clusterTM_EoverP_E"), cluster.energy() / abs_p, cluster.energy(), t);
        mHistManager.fill(HIST("clusterTM_dEtadPhi"), dEta, dPhi, t);
        mHistManager.fill(HIST("clusterEMatched"), cluster.energy(), t);
        mHistManager.fill(HIST("clusterTM_dEtaPt"), dEta, pT, t);
        mHistManager.fill(HIST("clusterTM_EvsP"), cluster.energy(), abs_p, t);
        mHistManager.fill(HIST("clusterTM_EoverP_Pt"), cluster.energy() / abs_p, match.track_as<tracksPID>().pt(), t);
        // A- and C-side
        if (match.track_as<tracksPID>().eta() > 0.0) {
          mHistManager.fill(HIST("clusterTM_dEtadPhi_ASide"), dEta, dPhi, t);
        } else if (match.track_as<tracksPID>().eta() < 0.0) {
          mHistManager.fill(HIST("clusterTM_dEtadPhi_CSide"), dEta, dPhi, t);
        }
        // positive and negative tracks seperate, with three different track momentum ranges
        if (match.track_as<tracksPID>().sign() == 1) {
          mHistManager.fill(HIST("clusterTM_PosdEtadPhi"), dEta, dPhi, t);
          mHistManager.fill(HIST("clusterTM_PosdPhiPt"), dPhi, pT, t);
          if (abs_p < 0.75) {
            mHistManager.fill(HIST("clusterTM_PosdEtadPhi_Pl0_75"), dEta, dPhi, t);
          } else if (abs_p >= 1.25) {
            mHistManager.fill(HIST("clusterTM_PosdEtadPhi_Pgeq1_25"), dEta, dPhi, t);
          } else {
            mHistManager.fill(HIST("clusterTM_PosdEtadPhi_0_75leqPl1_25"), dEta, dPhi, t);
          }
        } else {
          mHistManager.fill(HIST("clusterTM_NegdEtadPhi"), dEta, dPhi, t);
          mHistManager.fill(HIST("clusterTM_NegdPhiPt"), dPhi, pT, t);
          if (abs_p < 0.75) {
            mHistManager.fill(HIST("clusterTM_NegdEtadPhi_Pl0_75"), dEta, dPhi, t);
          } else if (abs_p >= 1.25) {
            mHistManager.fill(HIST("clusterTM_NegdEtadPhi_Pgeq1_25"), dEta, dPhi, t);
          } else {
            mHistManager.fill(HIST("clusterTM_NegdEtadPhi_0_75leqPl1_25"), dEta, dPhi, t);
          }
        }
        if (match.track_as<tracksPID>().tpcNSigmaEl() >= tpcNsigmaElectron->at(0) && match.track_as<tracksPID>().tpcNSigmaEl() <= tpcNsigmaElectron->at(1)) {                 // E/p for e+/e-
          if (usePionRejection && (match.track_as<tracksPID>().tpcNSigmaPi() <= tpcNsigmaPion->at(0) || match.track_as<tracksPID>().tpcNSigmaPi() >= tpcNsigmaPion->at(1))) { // with pion rejection
            mHistManager.fill(HIST("clusterTM_EoverP_electron"), cluster.energy() / abs_p, match.track_as<tracksPID>().pt(), t);
          } else {
            mHistManager.fill(HIST("clusterTM_EoverP_electron"), cluster.energy() / abs_p, match.track_as<tracksPID>().pt(), t);
          }
        } else if (match.track_as<tracksPID>().tpcNSigmaEl() >= tpcNsigmaBack->at(0) && match.track_as<tracksPID>().tpcNSigmaEl() >= tpcNsigmaBack->at(1)) { // E/p for hadrons / background
          mHistManager.fill(HIST("clusterTM_EoverP_hadron"), cluster.energy() / abs_p, match.track_as<tracksPID>().pt(), t);
        }
        if ((fabs(dEta) <= 0.01 + pow(match.track_as<tracksPID>().pt() + 4.07, -2.5)) &&
            (fabs(dPhi) <= 0.015 + pow(match.track_as<tracksPID>().pt() + 3.65, -2.)) &&
            cluster.energy() / abs_p < 1.75) {
          mHistManager.fill(HIST("clusterTM_dEtaTNAli"), dEta, tAli);
          mHistManager.fill(HIST("clusterTM_dPhiTNAli"), dPhi, tAli);
          mHistManager.fill(HIST("clusterTM_dRTNAli"), std::sqrt(dPhi * dPhi + dEta * dEta), tAli);
          tAli++;
        }
        t++;
      }
      mHistManager.fill(HIST("clusterTM_NTrackAli"), tAli);
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
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<TrackMatchingMonitor>(cfgc, TaskName{"emc-tmmonitor"})};
  return workflow;
}
