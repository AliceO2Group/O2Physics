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

/// \file electronSelectionWithTPCEMcal.cxx
/// \author Rashi Gupta <rashi.gupta@cern.ch>, IIT Indore
/// \author Ravindra Singh <ravindra.singh@cern.ch>, IIT Indore

#include "THnSparse.h"

#include "DataFormatsEMCAL/AnalysisCluster.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/PID/TPCPIDResponse.h"
#include "Common/Core/TrackSelection.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGHF/HFL/DataModel/ElectronSelectionTable.h"
#include "PWGJE/DataModel/EMCALClusters.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

auto massEl = o2::constants::physics::MassElectron;

struct HfElectronSelectionWithTPCEMcal {

  Produces<aod::HfSelEl> electronSel;
  // Configurables
  // Cluster information
  Configurable<bool> fillClusterInfo{"fillClusterInfo", true, "Fill histograms with EMCal cluster info before and after track match"};

  // Event Selection
  Configurable<float> mVertexCut{"mVertexCut", 10., "apply z-vertex cut with value in cm"};
  Configurable<bool> Isrun3{"Isrun3", true, "Data is from Run3 or Run2"};

  // Track selection
  Configurable<float> dcaXYTrackMax{"dcaXYTrackMax", 0.5f, "DCA XY cut"};
  Configurable<float> dcaZTrackMax{"dcaZTrackMax", 1.0f, "DCA Z cut"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.6f, "Eta range for electron tracks"};
  Configurable<float> etaTrackMin{"etaTrackMin", -0.6f, "Eta range for electron tracks"};
  Configurable<float> ptTrackMin{"ptTrackMin", 3.0f, "Transverse MOmentum range for electron tracks"};

  // EMcal and Dcal selection cut
  Configurable<float> etaTrackDCalNegativeMax{"etaTrackDCalNegativeMax", -0.22f, "Eta range for electron Dcal tracks"};
  Configurable<float> etaTrackDCalNegativeMin{"etaTrackDCalNegativeMin", -0.6f, "Eta range for electron tracks"};
  Configurable<float> etaTrackDCalPositiveMax{"etaTrackDCalPositiveMax", 0.6f, "Eta range for electron Dcal tracks"};
  Configurable<float> etaTrackDCalPositiveMin{"etaTrackDCalPositiveMin", 0.22f, "Eta range for electron tracks"};
  Configurable<float> phiTrackDCalMax{"phiTrackDCalMax", 3.3621f, "phi range for electron tracks associated Dcal"};
  Configurable<float> phiTrackDCalMin{"phiTrackDCalMin", 1.3955f, "phi range for electron tracks associated Dcal"};
  Configurable<float> phiTrackEMCalMax{"phiTrackEMCalMax", 5.708f, "phi range for electron tracks associated Emcal"};
  Configurable<float> phiTrackEMCalMin{"phiTrackEMCalMin", 4.5355f, "phi range for electron tracks associated Emcal"};

  // Track and Cluster matching cut
  Configurable<float> deltaEtaMatchMin{"deltaEtaMatchMin", 0.015f, "Eta distance of EMCAL cluster to its closest track"};
  Configurable<float> deltaPhiMatchMin{"deltaPhiMatchMin", 0.025f, "Phi distance of EMCAL cluster to its closest track"};
  Configurable<float> clusterTimeMax{"clusterTimeMax", 50.f, "Cluster time"};
  // Inclusive electron selection cut
  Configurable<float> eopElectronMin{"eopElectronMin", 0.8f, "Minimum E/p for electron tracks"};
  Configurable<float> eopElectronMax{"eopElectronMax", 1.2f, "Maximum E/p for electron tracks"};
  Configurable<float> m02ElectronMax{"m02ElectronMax", 0.9f, "max Electron M02"};
  Configurable<float> m02ElectronMin{"m02ElectronMin", 0.02f, "min Electron M02"};
  Configurable<float> m20ElectronMax{"m20ElectronMax", 1000.f, "max Electron M20"};
  Configurable<float> m20ElectronMin{"m20ElectronMin", 0.0f, "min Electron M20"};
  Configurable<float> tpcNsigmaElectronMin{"tpcNsigmaElectronMin", -0.5f, "min Electron TPCnsigma"};
  Configurable<float> tpcNsigmaElectronMax{"tpcNsigmaElectronMax", 3.0f, "max Electron TPCnsigma"};

  // Track selection cut for Mc Reco
  Configurable<float> mcRecDcaXYTrackMax{"mcRecDcaXYTrackMax", 0.5f, "MCReco DCA XY cut"};
  Configurable<float> mcRecDcaZTrackMax{"mcRecDcaZTrackMax", 1.0f, "MCReco DCA Z cut"};
  Configurable<float> mcRecEtaTrackMax{"mcRecEtaTrackMax", 0.6f, "MCReco Eta range for electron tracks"};
  Configurable<float> mcRecEtaTrackMin{"mcRecEtaTrackMin", -0.6f, "MCReco Eta range for electron tracks"};
  Configurable<float> mcRecPtTrackMin{"mcRecPtTrackMin", 3.0f, "McReco Transverse MOmentum range for electron tracks"};

  // EMcal and Dcal selection cut for Mc Reco
  Configurable<float> mcRecEtaTrackDCalNegativeMax{"mcRecEtaTrackDCalNegativeMax", -0.22f, "MCReco Eta range for electron Dcal tracks"};
  Configurable<float> mcRecEtaTrackDCalNegativeMin{"mcRecEtaTrackDCalNegativeMin", -0.6f, "MCReco Eta range for electron tracks"};
  Configurable<float> mcRecEtaTrackDCalPositiveMax{"mcRecEtaTrackDCalPositiveMax", 0.6f, "MCReco Eta range for electron Dcal tracks"};
  Configurable<float> mcRecEtaTrackDCalPositiveMin{"mcRecEtaTrackDCalPositiveMin", 0.22f, "MCReco Eta range for electron tracks"};
  Configurable<float> mcRecPhiTrackDCalMax{"mcRecPhiTrackDCalMax", 3.3621f, "McReco Phi range for electron tracks associated Dcal"};
  Configurable<float> mcRecPhiTrackDCalMin{"mcRecPhiTrackDCalMin", 1.3955f, "McReco Phi range for electron tracks associated Dcal"};
  Configurable<float> mcRecPhiTrackEMCalMax{"mcRecPhiTrackEMCalMax", 5.708f, "McReco Phi range for electron tracks associated Emcal"};
  Configurable<float> mcRecPhiTrackEMCalMin{"mcRecPhiTrackEMCalMin", 4.5355f, "McReco Phi range for electron tracks associated Emcal"};

  // Track and Cluster matching cut for Mc Reco
  Configurable<float> mcRecDeltaEtaMatchMin{"mcRecDeltaEtaMatchMin", 0.015f, "McReco Eta distance of EMCAL cluster to its closest track"};
  Configurable<float> mcRecDeltaPhiMatchMin{"mcRecDeltaPhiMatchMin", 0.025f, "McReco Phi distance of EMCAL cluster to its closest track"};
  Configurable<float> mcRecClusterTimeMax{"mcRecClusterTimeMax", 50.f, "McReco Cluster time"};

  // Inclusive electron selection cut for Mc Reco
  Configurable<float> mcRecM02ElectronMax{"mcRecM02ElectronMax", 0.9f, "MC Reco max Electron M02"};
  Configurable<float> mcRecM02ElectronMin{"mcRecM02ElectronMin", 0.02f, "MC Reco min Electron M02"};
  Configurable<float> mcRecM20ElectronMax{"mcRecM20ElectronMax", 1000.f, "MC Reco max Electron M20"};
  Configurable<float> mcRecM20ElectronMin{"mcRecM20ElectronMin", 0.0f, "MC Reco min Electron M20"};
  Configurable<float> mcRecTpcNsigmaElectronMin{"mcRecTpcNsigmaElectronMin", -0.5f, "MC Reco min Electron TPCnsigma"};
  Configurable<float> mcRecTpcNsigmaElectronMax{"mcRecTpcNsigmaElectronMax", 3.0f, "MC Reco max Electron TPCnsigma"};

  using CollisionTables = o2::soa::Filtered<o2::soa::Join<aod::Collisions, aod::Mults, aod::EvSels>>;
  using CollisionTable = CollisionTables::iterator;
  using TrackTables = o2::soa::Join<o2::aod::Tracks, o2::aod::TracksCov, o2::aod::TracksExtra, o2::aod::pidTPCFullEl, o2::aod::pidTOFFullEl, o2::aod::TracksDCA, o2::aod::TrackSelection, o2::aod::TrackSelectionExtension>;
  using TrackTable = TrackTables::iterator;

  using McCollisionTable = o2::soa::Filtered<o2::soa::Join<CollisionTables, aod::McCollisionLabels>>;
  using McCollisionTables = McCollisionTable::iterator;
  using McTrackTables = soa::Join<TrackTables, aod::McTrackLabels>;
  using McEMcalTable = soa::Join<o2::aod::EMCALClusters, aod::EMCALMCClusters>;

  Filter CollisionFilter = nabs(aod::collision::posZ) < mVertexCut && aod::collision::numContrib > (uint16_t)1;
  PresliceUnsorted<o2::aod::EMCALMatchedTracks> perClusterMatchedTracks = o2::aod::emcalmatchedtrack::trackId;

  HistogramConfigSpec hClusterInfoSpec{HistType::kTHnSparseD, {{300, 0.0, 30.0}, {100, -0.9, 0.9}, {200, 0, 6.3}, {50, 0, 50}, {1800, -900, 900}}};
  HistogramConfigSpec hDeltaPhiDeltaEtaClusterTrackSpecEnergy{HistType::kTHnSparseD, {{400, -0.2, 0.2}, {400, -0.2, 0.2}, {600, -300, 300}, {300, 0.0, 30.0}}};
  HistogramConfigSpec hPIDSpec{HistType::kTHnSparseD, {{60, 0, 3}, {500, 0.0, 50.0}, {500, 0., 50.}, {300, -15, 15}, {300, 0.0, 30.0}, {400, 0, 2}, {400, 0, 2}}};
  HistogramConfigSpec hTrackAllInfoSpec{HistType::kTHnSparseD, {{480, 0, 160}, {300, -15, 15}, {500, 0., 50.}, {500, 0., 50.}, {100, -1.5, 1.5}, {100, 0, 7}, {3, 0, 3}}};
  HistogramConfigSpec hTrackInfoSpec{HistType::kTHnSparseD, {{60, 0, 3}, {480, 0, 160}, {300, -15, 15}, {500, 0., 50.}, {500, 0., 50.}, {100, -1.5, 1.5}, {100, 0, 7}}};

  HistogramRegistry registry{
    "registry",
    {{"hNevents", "No of events", {HistType::kTH1F, {{3, 1, 4}}}},
     {"hZvertex", "z vertex", {HistType::kTH1F, {{100, -100, 100}}}},

     {"hTrackInformation", "Sparse TPC info; dE/dx;n#sigma;#it{p} (GeV#it{/c});#it{p}_{T} (GeV#it{/c});#eta;#varphi;passEMcal;", hTrackAllInfoSpec},
     {"hClusterInformationBefore", "Cluster Info before match; Energy (GeV);#eta;#varphi", hClusterInfoSpec},
     {"hClusterInformationAfter", "Cluster Info after match; Energy (GeV);#eta;#varphi", hClusterInfoSpec},
     {"hPIDafterMatch", "PID Info after match; E/P; dE/dx;n#sigma;#it{p} (GeV#it{/c});#it{p}_{T} (GeV#it{/c});#eta;#varphi;", hTrackInfoSpec},
     {"hEPRatioafterPID", "E/P Ratio after PID Cuts apply only trackwodca filter", {HistType::kTH2F, {{60, 0, 3}, {100, 0, 10}}}},
     {"hPIDafterPIDcuts", "PID Info after PID cuts; E/P; #it{p} (GeV#it{/c});#it{p}_{T} (GeV#it{/c});n_{#sigma}^{e};GeV;M02;M20", hPIDSpec},
     {"hClsTrkEtaPhiDiffTimeEnergy", "ClsTrkEtaPhiDiffTimeEnergy;#Delta#eta;#Delta#varphi;Sec; Energy (GeV)", hDeltaPhiDeltaEtaClusterTrackSpecEnergy}}};

  void init(InitContext&)
  {
    registry.get<THnSparse>(HIST("hTrackInformation"))->Sumw2();
    registry.get<THnSparse>(HIST("hClusterInformationBefore"))->Sumw2();
    registry.get<THnSparse>(HIST("hClusterInformationAfter"))->Sumw2();
    registry.get<THnSparse>(HIST("hPIDafterMatch"))->Sumw2();
    registry.get<THnSparse>(HIST("hPIDafterPIDcuts"))->Sumw2();
    registry.get<THnSparse>(HIST("hClsTrkEtaPhiDiffTimeEnergy"))->Sumw2();
  }

  // Track Selection Cut
  template <typename T>
  bool selTracks(T const& track)
  {
    if (!track.isGlobalTrackWoDCA())
      return false;
    if (std::abs(track.dcaXY()) > dcaXYTrackMax || std::abs(track.dcaZ()) > dcaZTrackMax)
      return false;
    if (track.eta() < etaTrackMin || track.eta() > etaTrackMax)
      return false;
    if ((track.phi() < phiTrackEMCalMin || track.phi() > phiTrackEMCalMax) && (track.phi() < phiTrackDCalMin || track.phi() > phiTrackDCalMax))
      return false;
    if (track.pt() < ptTrackMin)
      return false;
    return true;
  }

  // MC Reco Track Selection Cut
  template <typename T>
  bool mcSelTracks(T const& mcTrack)
  {
    if (!mcTrack.isGlobalTrackWoDCA())
      return false;
    if (std::abs(mcTrack.dcaXY()) > mcRecDcaXYTrackMax || std::abs(mcTrack.dcaZ()) > mcRecDcaZTrackMax)
      return false;
    if (mcTrack.eta() < mcRecEtaTrackMin || mcTrack.eta() > mcRecEtaTrackMax)
      return false;
    if ((mcTrack.phi() < mcRecPhiTrackEMCalMin || mcTrack.phi() > mcRecPhiTrackEMCalMax) && (mcTrack.phi() < mcRecPhiTrackDCalMin || mcTrack.phi() > mcRecPhiTrackDCalMax))
      return false;
    if (mcTrack.pt() < mcRecPtTrackMin)
      return false;
    return true;
  }

  // Electron Identification
  template <bool IsData, bool IsMC, typename TracksType, typename ClusterType, typename MatchType, typename CollisionType>
  void fillElectronTrack(CollisionType const& collision, TracksType const& tracks, ClusterType const& cluster, MatchType const& matchedTracks)
  {
    if (!(Isrun3 ? collision.sel8() : (collision.sel7() && collision.alias_bit(kINT7))))
      return;

    registry.fill(HIST("hNevents"), 1);
    registry.fill(HIST("hZvertex"), collision.posZ());

    /////////////////////////////////
    // cluster info before match ///
    ///////////////////////////////
    if (fillClusterInfo) {
      for (const auto& clusterBefore : cluster) {
        registry.fill(HIST("hClusterInformationBefore"), clusterBefore.energy(), clusterBefore.eta(), clusterBefore.phi(), clusterBefore.nCells(), clusterBefore.time());
      }
    }
    int passEMCal;
    float phiTrack = -999;
    float etaTrack = -999;
    float pTrack = -999;
    float ptTrack = -999;
    float dcaxyTrack = -999;
    float dcazTrack = -999;
    float tpcNsigmaTrack = -999;

    for (const auto& track : tracks) {

      phiTrack = track.phi();
      etaTrack = track.eta();
      pTrack = track.p();
      ptTrack = track.pt();
      dcaxyTrack = track.dcaXY();
      dcazTrack = track.dcaZ();
      tpcNsigmaTrack = track.tpcNSigmaEl();

      // Apply Track Selection

      if constexpr (IsData) {
        if (!selTracks(track))
          continue;
      }
      if constexpr (IsMC) {
        if (!mcSelTracks(track))
          continue;
      }
      passEMCal = 0;

      if ((phiTrack > phiTrackEMCalMin && phiTrack < phiTrackEMCalMax) && (etaTrack > etaTrackMin && etaTrack < etaTrackMax))
        passEMCal = 1; // EMcal acceptance passed
      if ((phiTrack > phiTrackDCalMin && phiTrack < phiTrackDCalMax) && ((etaTrack > etaTrackDCalPositiveMin && etaTrack < etaTrackDCalPositiveMax) || (etaTrack > etaTrackDCalNegativeMin && etaTrack < etaTrackDCalNegativeMax)))
        passEMCal = 2; // Dcal acceptance passed

      registry.fill(HIST("hTrackInformation"), track.tpcSignal(), tpcNsigmaTrack, pTrack, ptTrack, etaTrack, phiTrack, passEMCal); // track infor after filter bit

      auto tracksofcluster = matchedTracks.sliceBy(perClusterMatchedTracks, track.globalIndex());
      float phiMatchTrack = -999;
      float etaMatchTrack = -999;
      float pMatchTrack = -999;
      float ptMatchTrack = -999;
      float tpcNsigmaMatchTrack = -999;
      float phiMatchCluster = -999;
      float etaMatchCluster = -999;
      float eMatchCluster = -999;
      float m02MatchCluster = -999;
      float m20MatchCluster = -999;
      float timeCluster = -999;
      float cellCluster = -999;
      float deltaPhiMatch = -999.;
      float deltaEtaMatch = -999.;
      float eop = -999;
      bool isEMcal = 0;

      float trackRapidity = track.rapidity(massEl);

      for (const auto& ematchTrack : tracksofcluster) {

        auto matchTrack = ematchTrack.template track_as<TracksType>();

        auto cluster = ematchTrack.template emcalcluster_as<ClusterType>();

        phiMatchTrack = matchTrack.phi();
        etaMatchTrack = matchTrack.eta();
        pMatchTrack = matchTrack.p();
        ptMatchTrack = matchTrack.pt();
        tpcNsigmaMatchTrack = matchTrack.tpcNSigmaEl();
        phiMatchCluster = cluster.phi();
        etaMatchCluster = cluster.eta();
        eMatchCluster = cluster.energy();
        m02MatchCluster = cluster.m02();
        m20MatchCluster = cluster.m20();
        timeCluster = cluster.time();
        cellCluster = cluster.nCells();

        deltaPhiMatch = matchTrack.trackPhiEmcal() - phiMatchCluster;
        deltaEtaMatch = matchTrack.trackEtaEmcal() - etaMatchCluster;

        // Track and cluster Matching

        if constexpr (IsData) {
          if (std::abs(timeCluster) > clusterTimeMax)
            continue;
          if (std::abs(deltaPhiMatch) > deltaPhiMatchMin || std::abs(deltaEtaMatch) > deltaEtaMatchMin)
            continue;
        }
        if constexpr (IsMC) {
          if (std::abs(timeCluster) > mcRecClusterTimeMax)
            continue;
          if (std::abs(deltaPhiMatch) > mcRecDeltaPhiMatchMin || std::abs(deltaEtaMatch) > mcRecDeltaEtaMatchMin)
            continue;
        }

        registry.fill(HIST("hClsTrkEtaPhiDiffTimeEnergy"), deltaEtaMatch, deltaPhiMatch, cluster.time(), eMatchCluster);

        if (fillClusterInfo)
          registry.fill(HIST("hClusterInformationAfter"), eMatchCluster, etaMatchCluster, phiMatchCluster, cluster.nCells(), timeCluster);
        eop = eMatchCluster / pMatchTrack;
        registry.fill(HIST("hPIDafterMatch"), eop, matchTrack.tpcSignal(), tpcNsigmaMatchTrack, pMatchTrack, ptMatchTrack, etaMatchTrack, phiMatchTrack);

        // Apply Electron Identification cuts
        if constexpr (IsData) {
          if ((tpcNsigmaMatchTrack < tpcNsigmaElectronMin || tpcNsigmaMatchTrack > tpcNsigmaElectronMax) || (m02MatchCluster < m02ElectronMin || m02MatchCluster > m02ElectronMax) || (m20MatchCluster < m20ElectronMin || m20MatchCluster > m20ElectronMax))
            continue;
        }
        if constexpr (IsMC) {

          if ((tpcNsigmaMatchTrack < mcRecTpcNsigmaElectronMin || tpcNsigmaMatchTrack > mcRecTpcNsigmaElectronMax) || (m02MatchCluster < mcRecM02ElectronMin || m02MatchCluster > mcRecM02ElectronMax) || (m20MatchCluster < mcRecM20ElectronMin || m20MatchCluster > mcRecM20ElectronMax))
            continue;
        }

        registry.fill(HIST("hEPRatioafterPID"), eop, ptMatchTrack);
        if (eop < eopElectronMin || eop > eopElectronMax)
          continue;

        registry.fill(HIST("hPIDafterPIDcuts"), eop, pMatchTrack, ptMatchTrack, tpcNsigmaMatchTrack, eMatchCluster, m02MatchCluster, m20MatchCluster);
        isEMcal = 1;
        electronSel(matchTrack.collisionId(), matchTrack.globalIndex(), matchTrack.eta(), matchTrack.phi(), ptMatchTrack, pMatchTrack, trackRapidity, matchTrack.dcaXY(), matchTrack.dcaZ(), matchTrack.tpcNSigmaEl(), matchTrack.tofNSigmaEl(),
                    eMatchCluster, etaMatchCluster, phiMatchCluster, m02MatchCluster, m20MatchCluster, cellCluster, timeCluster, deltaEtaMatch, deltaPhiMatch, isEMcal);
      }

      /// Electron information without Emcal and use TPC and TOF
      if (isEMcal == 1)
        continue;
      electronSel(track.collisionId(), track.globalIndex(), etaTrack, phiTrack, ptTrack, pTrack, trackRapidity, dcaxyTrack, dcazTrack, track.tpcNSigmaEl(), track.tofNSigmaEl(),
                  eMatchCluster, etaMatchCluster, phiMatchCluster, m02MatchCluster, m20MatchCluster, cellCluster, timeCluster, deltaEtaMatch, deltaPhiMatch, isEMcal);
    }
  }

  void processData(CollisionTable const& collision, aod::EMCALClusters const& mAnalysisClusters, o2::aod::EMCALMatchedTracks const& matchedTracks, TrackTables const& tracks)
  {
    fillElectronTrack<true, false>(collision, tracks, mAnalysisClusters, matchedTracks);
  }
  PROCESS_SWITCH(HfElectronSelectionWithTPCEMcal, processData, "process Data info only", true);

  // group according to reconstructed Collisions
  void processMcRec(McCollisionTables const& mccollision, McTrackTables const& mctracks, McEMcalTable const& mcClusters, o2::aod::EMCALMatchedTracks const& matchedTracks)
  {
    fillElectronTrack<false, true>(mccollision, mctracks, mcClusters, matchedTracks);
  }
  PROCESS_SWITCH(HfElectronSelectionWithTPCEMcal, processMcRec, "Process MC Reco mode", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfElectronSelectionWithTPCEMcal>(cfgc)};
}
