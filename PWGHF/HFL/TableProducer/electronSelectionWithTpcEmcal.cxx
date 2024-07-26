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

/// \file electronSelectionWithTpcEmcal.cxx
/// \brief Task used to electron selection with tpc and emcal.
/// \author Rashi Gupta <rashi.gupta@cern.ch>, IIT Indore
/// \author Ravindra Singh <ravindra.singh@cern.ch>, IIT Indore

#include <boost/move/detail/meta_utils_core.hpp>
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

#include "PWGJE/DataModel/EMCALClusters.h"

#include "PWGHF/HFL/DataModel/ElectronSelectionTable.h"

using namespace o2;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct HfElectronSelectionWithTpcEmcal {

  Produces<aod::HfSelEl> electronSel;
  // Configurables
  // EMCal Cluster information
  Configurable<bool> fillEmcClusterInfo{"fillEmcClusterInfo", true, "Fill histograms with EMCal cluster info before and after track match"};

  // Event Selection
  Configurable<float> zPvPosMax{"zPvPosMax", 10., "Maximum z of the primary vertex (cm)"};
  Configurable<bool> isRun3{"isRun3", true, "Data is from Run3 or Run2"};

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

  // Track and  EMCal Cluster matching cut
  Configurable<float> deltaEtaMatchMin{"deltaEtaMatchMin", -0.013f, "Min Eta distance of EMCAL cluster to its closest track"};
  Configurable<float> deltaEtaMatchMax{"deltaEtaMatchMax", 0.0171f, "Max Eta distance of EMCAL cluster to its closest track"};
  Configurable<float> deltaPhiMatchMin{"deltaPhiMatchMin", -0.022f, "Min Phi distance of EMCAL cluster to its closest track"};
  Configurable<float> deltaPhiMatchMax{"deltaPhiMatchMax", 0.028f, "Max Phi distance of EMCAL cluster to its closest track"};
  Configurable<float> timeEmcClusterMax{"timeEmcClusterMax", 50.f, "EMCal Cluster time"};

  // Inclusive electron selection cut
  Configurable<float> eopElectronMin{"eopElectronMin", 0.8f, "Minimum E/p for electron tracks"};
  Configurable<float> eopElectronMax{"eopElectronMax", 1.2f, "Maximum E/p for electron tracks"};
  Configurable<float> m02EmcClusterElectronMax{"m02EmcClusterElectronMax", 0.9f, "max Electron  EMCal Cluster M02"};
  Configurable<float> m02EmcClusterElectronMin{"m02EmcClusterElectronMin", 0.02f, "min Electron  EMCal Cluster M02"};
  Configurable<float> m20EmcClusterElectronMax{"m20EmcClusterElectronMax", 1000.f, "max Electron  EMCal Cluster M20"};
  Configurable<float> m20EmcClusterElectronMin{"m20EmcClusterElectronMin", 0.0f, "min Electron  EMCal Cluster M20"};
  Configurable<float> tpcNsigmaElectronMin{"tpcNsigmaElectronMin", -0.5f, "min Electron TPCnsigma"};
  Configurable<float> tpcNsigmaElectronMax{"tpcNsigmaElectronMax", 3.0f, "max Electron TPCnsigma"};

  // Track and EMCal Cluster matching cut for Mc Reco
  Configurable<float> mcRecDeltaEtaMatchMin{"mcRecDeltaEtaMatchMin", -0.013f, "McReco Min Eta distance of EMCAL cluster to its closest track"};
  Configurable<float> mcRecDeltaEtaMatchMax{"mcRecDeltaEtaMatchMax", 0.0171f, "McReco Max Eta distance of EMCAL cluster to its closest track"};
  Configurable<float> mcRecDeltaPhiMatchMin{"mcRecDeltaPhiMatchMin", -0.022f, "McReco Min Phi distance of EMCAL cluster to its closest track"};
  Configurable<float> mcRecDeltaPhiMatchMax{"mcRecDeltaPhiMatchMax", 0.028f, "McReco Max Phi distance of EMCAL cluster to its closest track"};

  Configurable<float> mcRecTimeEmcClusterMax{"mcRecTimeEmcClusterMax", 50.f, "McReco EMCal Cluster time"};

  // Inclusive electron selection cut for Mc Reco
  Configurable<float> mcRecM02EmcClusterElectronMax{"mcRecM02EmcClusterElectronMax", 0.9f, "MC Reco max Electron EMCal Cluster M02"};
  Configurable<float> mcRecM02EmcClusterElectronMin{"mcRecM02EmcClusterElectronMin", 0.02f, "MC Reco min Electron  EMCal Cluster M02"};
  Configurable<float> mcRecM20EmcClusterElectronMax{"mcRecM20EmcClusterElectronMax", 1000.f, "MC Reco max Electron  EMCal Cluster M20"};
  Configurable<float> mcRecM20EmcClusterElectronMin{"mcRecM20EmcClusterElectronMin", 0.0f, "MC Reco min Electron   EMCal Cluster M20"};
  Configurable<float> mcRecTpcNsigmaElectronMin{"mcRecTpcNsigmaElectronMin", -0.5f, "MC Reco min Electron TPCnsigma"};
  Configurable<float> mcRecTpcNsigmaElectronMax{"mcRecTpcNsigmaElectronMax", 3.0f, "MC Reco max Electron TPCnsigma"};

  using TableCollisions = o2::soa::Filtered<o2::soa::Join<aod::Collisions, aod::Mults, aod::EvSels>>;
  using TableCollision = TableCollisions::iterator;
  using TableTracks = o2::soa::Join<o2::aod::Tracks, o2::aod::TracksCov, o2::aod::TracksExtra, o2::aod::pidTPCFullEl, o2::aod::pidTOFFullEl, o2::aod::TracksDCA, o2::aod::TrackSelection, o2::aod::TrackSelectionExtension>;

  using McTableCollisions = o2::soa::Filtered<o2::soa::Join<TableCollisions, aod::McCollisionLabels>>;
  using McTableCollision = McTableCollisions::iterator;
  using McTableTracks = soa::Join<TableTracks, aod::McTrackLabels>;
  using McTableEmcals = soa::Join<o2::aod::EMCALClusters, aod::EMCALMCClusters>;

  Filter CollisionFilter = nabs(aod::collision::posZ) < zPvPosMax && aod::collision::numContrib > (uint16_t)1;
  PresliceUnsorted<o2::aod::EMCALMatchedTracks> perClusterMatchedTracks = o2::aod::emcalmatchedtrack::trackId;

  HistogramConfigSpec hEmcClusterInfoSpec{HistType::kTHnSparseD, {{300, 0.0, 30.0}, {100, -0.9, 0.9}, {200, 0, 6.3}, {50, 0, 50}, {1800, -900, 900}}};
  HistogramConfigSpec hDeltaPhiDeltaEtaEmcClusterTrackSpecEnergy{HistType::kTHnSparseD, {{400, -0.2, 0.2}, {400, -0.2, 0.2}, {600, -300, 300}, {300, 0.0, 30.0}}};
  HistogramConfigSpec hPIDSpec{HistType::kTHnSparseD, {{60, 0, 3}, {500, 0.0, 50.0}, {500, 0., 50.}, {300, -15, 15}, {300, 0.0, 30.0}, {400, 0, 2}, {400, 0, 2}}};
  HistogramConfigSpec hTrackAllInfoSpec{HistType::kTHnSparseD, {{480, 0, 160}, {300, -15, 15}, {500, 0., 50.}, {500, 0., 50.}, {100, -1.5, 1.5}, {100, 0, 7}, {3, 0, 3}}};
  HistogramConfigSpec hTrackInfoSpec{HistType::kTHnSparseD, {{60, 0, 3}, {480, 0, 160}, {300, -15, 15}, {500, 0., 50.}, {500, 0., 50.}, {100, -1.5, 1.5}, {100, 0, 7}}};

  HistogramRegistry registry{
    "registry",
    {{"hNevents", "No of events", {HistType::kTH1F, {{3, 1, 4}}}},
     {"hZvertex", "z vertex", {HistType::kTH1F, {{100, -100, 100}}}},
     {"hTrackInformation", "Sparse TPC info; dE/dx;n#sigma;#it{p} (GeV#it{/c});#it{p}_{T} (GeV#it{/c});#eta;#varphi;passEMcal;", hTrackAllInfoSpec},
     {"hEmcClusterInformationBefore", "EMCal Cluster Info before match; Energy (GeV);#eta;#varphi", hEmcClusterInfoSpec},
     {"hEmcClusterInformationAfter", "EMCal Cluster Info after match; Energy (GeV);#eta;#varphi", hEmcClusterInfoSpec},
     {"hPIDafterMatch", "PID Info after match; E/P; dE/dx;n#sigma;#it{p} (GeV#it{/c});#it{p}_{T} (GeV#it{/c});#eta;#varphi;", hTrackInfoSpec},
     {"hEPRatioafterPID", "E/P Ratio after PID Cuts apply only trackwodca filter", {HistType::kTH2F, {{60, 0, 3}, {100, 0, 10}}}},
     {"hPIDafterPIDcuts", "PID Info after PID cuts; E/P; #it{p} (GeV#it{/c});#it{p}_{T} (GeV#it{/c});n_{#sigma}^{e};GeV;M02;M20", hPIDSpec},
     {"hEmcClsTrkEtaPhiDiffTimeEnergy", "EmcClsTrkEtaPhiDiffTimeEnergy;#Delta#eta;#Delta#varphi;Sec; Energy (GeV)", hDeltaPhiDeltaEtaEmcClusterTrackSpecEnergy}}};

  void init(InitContext&)
  {
    registry.get<THnSparse>(HIST("hTrackInformation"))->Sumw2();
    registry.get<THnSparse>(HIST("hEmcClusterInformationBefore"))->Sumw2();
    registry.get<THnSparse>(HIST("hEmcClusterInformationAfter"))->Sumw2();
    registry.get<THnSparse>(HIST("hPIDafterMatch"))->Sumw2();
    registry.get<THnSparse>(HIST("hPIDafterPIDcuts"))->Sumw2();
    registry.get<THnSparse>(HIST("hEmcClsTrkEtaPhiDiffTimeEnergy"))->Sumw2();
  }

  // Track Selection Cut
  template <typename T>
  bool selTracks(T const& track)
  {
    if (!track.isGlobalTrackWoDCA()) {
      return false;
    }
    if (std::abs(track.dcaXY()) > dcaXYTrackMax || std::abs(track.dcaZ()) > dcaZTrackMax) {
      return false;
    }
    if (track.eta() < etaTrackMin || track.eta() > etaTrackMax) {
      return false;
    }
    if ((track.phi() < phiTrackEMCalMin || track.phi() > phiTrackEMCalMax) && (track.phi() < phiTrackDCalMin || track.phi() > phiTrackDCalMax)) {
      return false;
    }
    if (track.pt() < ptTrackMin) {
      return false;
    }
    return true;
  }

  // Electron Identification
  template <bool isMc, typename TracksType, typename EmcClusterType, typename MatchType, typename CollisionType, typename ParticleType>
  void fillElectronTrack(CollisionType const& collision, TracksType const& tracks, EmcClusterType const& emcClusters, MatchType const& matchedTracks, ParticleType const& /*particlemc*/)
  {
    if (!(isRun3 ? collision.sel8() : (collision.sel7() && collision.alias_bit(kINT7))))
      return;

    registry.fill(HIST("hNevents"), 1);
    registry.fill(HIST("hZvertex"), collision.posZ());

    /////////////////////////////////
    // EMCal cluster info before match ///
    ///////////////////////////////
    if (fillEmcClusterInfo) {
      for (const auto& emcClusterBefore : emcClusters) {
        registry.fill(HIST("hEmcClusterInformationBefore"), emcClusterBefore.energy(), emcClusterBefore.eta(), emcClusterBefore.phi(), emcClusterBefore.nCells(), emcClusterBefore.time());
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
      if (!selTracks(track)) {
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
      float phiMatchEmcCluster = -999;
      float etaMatchEmcCluster = -999;
      float eMatchEmcCluster = -999;
      float m02MatchEmcCluster = -999;
      float m20MatchEmcCluster = -999;
      float timeEmcCluster = -999;
      float cellEmcCluster = -999;
      float deltaPhiMatch = -999.;
      float deltaEtaMatch = -999.;
      float eop = -999;
      bool isEMcal = false;

      float trackRapidity = track.rapidity(MassElectron);

      for (const auto& ematchTrack : tracksofcluster) {

        auto matchTrack = ematchTrack.template track_as<TracksType>();

        auto emcCluster = ematchTrack.template emcalcluster_as<EmcClusterType>();

        phiMatchTrack = matchTrack.phi();
        etaMatchTrack = matchTrack.eta();
        pMatchTrack = matchTrack.p();
        ptMatchTrack = matchTrack.pt();
        tpcNsigmaMatchTrack = matchTrack.tpcNSigmaEl();
        phiMatchEmcCluster = emcCluster.phi();
        etaMatchEmcCluster = emcCluster.eta();
        eMatchEmcCluster = emcCluster.energy();
        m02MatchEmcCluster = emcCluster.m02();
        m20MatchEmcCluster = emcCluster.m20();
        timeEmcCluster = emcCluster.time();
        cellEmcCluster = emcCluster.nCells();

        deltaPhiMatch = matchTrack.trackPhiEmcal() - phiMatchEmcCluster;
        deltaEtaMatch = matchTrack.trackEtaEmcal() - etaMatchEmcCluster;

        // Track and EMCal cluster Matching

        if constexpr (!isMc) {
          if (std::abs(timeEmcCluster) > timeEmcClusterMax) {
            continue;
          }
          if (deltaPhiMatch < deltaPhiMatchMin || deltaPhiMatch > deltaPhiMatchMax || deltaEtaMatch < deltaEtaMatchMin || deltaEtaMatch > deltaEtaMatchMax) {
            continue;
          }
        } else {
          if (std::abs(timeEmcCluster) > mcRecTimeEmcClusterMax) {
            continue;
          }
          if (deltaPhiMatch < mcRecDeltaPhiMatchMin || deltaPhiMatch > mcRecDeltaPhiMatchMax || deltaEtaMatch < mcRecDeltaEtaMatchMin || deltaEtaMatch > mcRecDeltaEtaMatchMax) {
            continue;
          }
        }

        registry.fill(HIST("hEmcClsTrkEtaPhiDiffTimeEnergy"), deltaEtaMatch, deltaPhiMatch, timeEmcCluster, eMatchEmcCluster);

        if (fillEmcClusterInfo)
          registry.fill(HIST("hEmcClusterInformationAfter"), eMatchEmcCluster, etaMatchEmcCluster, phiMatchEmcCluster, cellEmcCluster, timeEmcCluster);
        eop = eMatchEmcCluster / pMatchTrack;
        registry.fill(HIST("hPIDafterMatch"), eop, matchTrack.tpcSignal(), tpcNsigmaMatchTrack, pMatchTrack, ptMatchTrack, etaMatchTrack, phiMatchTrack);

        // Apply Electron Identification cuts
        if constexpr (!isMc) {
          if ((tpcNsigmaMatchTrack < tpcNsigmaElectronMin || tpcNsigmaMatchTrack > tpcNsigmaElectronMax) || (m02MatchEmcCluster < m02EmcClusterElectronMin || m02MatchEmcCluster > m02EmcClusterElectronMax) || (m20MatchEmcCluster < m20EmcClusterElectronMin || m20MatchEmcCluster > m20EmcClusterElectronMax)) {
            continue;
          }
        } else {
          if ((tpcNsigmaMatchTrack < mcRecTpcNsigmaElectronMin || tpcNsigmaMatchTrack > mcRecTpcNsigmaElectronMax) || (m02MatchEmcCluster < mcRecM02EmcClusterElectronMin || m02MatchEmcCluster > mcRecM02EmcClusterElectronMax) || (m20MatchEmcCluster < mcRecM20EmcClusterElectronMin || m20MatchEmcCluster > mcRecM20EmcClusterElectronMax)) {
            continue;
          }
        }

        registry.fill(HIST("hEPRatioafterPID"), eop, ptMatchTrack);
        if (eop < eopElectronMin || eop > eopElectronMax) {
          continue;
        }
        registry.fill(HIST("hPIDafterPIDcuts"), eop, pMatchTrack, ptMatchTrack, tpcNsigmaMatchTrack, eMatchEmcCluster, m02MatchEmcCluster, m20MatchEmcCluster);

        isEMcal = true;
        electronSel(matchTrack.collisionId(), matchTrack.globalIndex(), etaMatchTrack, phiMatchTrack, ptMatchTrack, pMatchTrack, trackRapidity, matchTrack.dcaXY(), matchTrack.dcaZ(), matchTrack.tpcNSigmaEl(), matchTrack.tofNSigmaEl(),
                    eMatchEmcCluster, etaMatchEmcCluster, phiMatchEmcCluster, m02MatchEmcCluster, m20MatchEmcCluster, cellEmcCluster, timeEmcCluster, deltaEtaMatch, deltaPhiMatch, isEMcal);
      }

      /// Electron information without Emcal and use TPC and TOF
      if (isEMcal) {
        continue;
      }
      electronSel(track.collisionId(), track.globalIndex(), etaTrack, phiTrack, ptTrack, pTrack, trackRapidity, dcaxyTrack, dcazTrack, track.tpcNSigmaEl(), track.tofNSigmaEl(),
                  eMatchEmcCluster, etaMatchEmcCluster, phiMatchEmcCluster, m02MatchEmcCluster, m20MatchEmcCluster, cellEmcCluster, timeEmcCluster, deltaEtaMatch, deltaPhiMatch, isEMcal);
    }
  }

  ///  Electron selection - for real data and data-like analysis
  void processData(TableCollision const& collision,
                   TableTracks const& tracks,
                   aod::EMCALClusters const& emcClusters,
                   o2::aod::EMCALMatchedTracks const& matchedTracks)
  {
    fillElectronTrack<false>(collision, tracks, emcClusters, matchedTracks, 0);
  }
  PROCESS_SWITCH(HfElectronSelectionWithTpcEmcal, processData, "process Data info only", true);

  ///  Electron selection - for MC reco-level analysis
  void processMcRec(McTableCollision const& mcCollision,
                    McTableTracks const& mcTracks,
                    McTableEmcals const& mcEmcClusters,
                    o2::aod::EMCALMatchedTracks const& matchedTracks,
                    aod::McParticles const& mcParticles)
  {
    fillElectronTrack<true>(mcCollision, mcTracks, mcEmcClusters, matchedTracks, mcParticles);
  }
  PROCESS_SWITCH(HfElectronSelectionWithTpcEmcal, processMcRec, "Process MC Reco mode", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfElectronSelectionWithTpcEmcal>(cfgc)};
}
