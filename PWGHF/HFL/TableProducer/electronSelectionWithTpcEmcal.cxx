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

const int etaAxisBins = 100;
const float trackEtaAxisMin = -1.5;
const float trackEtaAxisMax = 1.5;
const int phiAxisBins = 100;
const float trackPhiAxisMin = 0.;
const float trackPhiAxisMax = o2::constants::math::TwoPI;
const int passEMCalBins = 3;
const int passEMCalAxisMin = 0.;
const int passEMCalAxisMax = 3;
const int eopAxisBins = 60;
const float eopAxisMin = 0.;
const float eopAxisMax = 3.0;
const int pAxisBins = 500;
const float pAxisMin = 0.;
const float pAxisMax = 50.0;
const int m02AxisBins = 100;
const float m02AxisMin = 0.;
const float m02AxisMax = 2.0;
const int m20AxisBins = 100;
const float m20AxisMin = 0.;
const float m20AxisMax = 2.0;
const int nSigmaAxisBins = 300;
const float nSigmaAxisMin = -15.;
const float nSigmaAxisMax = 15.;
const int dEdxAxisBins = 480;
const float dEdxAxisMin = 0.;
const float dEdxAxisMax = 160.;
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

  Filter CollisionFilter = nabs(aod::collision::posZ) < zPvPosMax && aod::collision::numContrib > static_cast<uint16_t>(1);
  PresliceUnsorted<o2::aod::EMCALMatchedTracks> perClusterMatchedTracks = o2::aod::emcalmatchedtrack::trackId;

  HistogramConfigSpec hEmcClusterEnergySpec{HistType::kTH1F, {{300, 0.0, 30.0}}};
  HistogramConfigSpec hEmcClusterEtaPhiSpec{HistType::kTH2F, {{100, -0.9, 0.9}, {200, 0, 6.3}}};
  HistogramConfigSpec hEmcClusterEnergyCellSpec{HistType::kTH2F, {{400, 0.0, 30.0}, {50, 0, 50}}};
  HistogramConfigSpec hEmcClusterEnergyTimeSpec{HistType::kTH2F, {{300, 0.0, 30.0}, {1800, -900, 900}}};

  HistogramConfigSpec hDeltaPhiDeltaEtaEmcClusterTrackSpecEnergy{HistType::kTH3F, {{400, -0.2, 0.2}, {400, -0.2, 0.2}, {600, -300, 300}}};
  HistogramConfigSpec hAfterMatchEoPSigamSpec{HistType::kTHnSparseD, {{eopAxisBins, eopAxisMin, eopAxisMax}, {pAxisBins, pAxisMin, pAxisMax}, {nSigmaAxisBins, nSigmaAxisMin, nSigmaAxisMax}, {m02AxisBins, m02AxisMin, m02AxisMax}, {m20AxisBins, m20AxisMin, m20AxisMax}}};

  HistogramConfigSpec hTrackEnergyLossSpec{HistType::kTH3F, {{dEdxAxisBins, dEdxAxisMin, dEdxAxisMax}, {pAxisBins, pAxisMin, pAxisMax}, {passEMCalBins, passEMCalAxisMin, passEMCalAxisMax}}};

  HistogramConfigSpec hTracknSigmaSpec{HistType::kTH3F, {{nSigmaAxisBins, nSigmaAxisMin, nSigmaAxisMax}, {pAxisBins, pAxisMin, pAxisMax}, {passEMCalBins, passEMCalAxisMin, passEMCalAxisMax}}};

  HistogramRegistry registry{
    "registry",
    {{"hNevents", "No of events", {HistType::kTH1F, {{3, 1, 4}}}},
     {"hZvertex", "z vertex", {HistType::kTH1F, {{100, -100, 100}}}},
     {"hEmcClusterM02", "m02", {HistType::kTH1F, {{m02AxisBins, m02AxisMin, m02AxisMax}}}},
     {"hEmcClusterM20", "m20", {HistType::kTH1F, {{m20AxisBins, m20AxisMin, m20AxisMax}}}},
     {"hTrackEtaPhi", "TPC EtaPhi Info; #eta;#varphi;passEMcal;", {HistType::kTH3F, {{etaAxisBins, trackEtaAxisMin, trackEtaAxisMax}, {phiAxisBins, trackPhiAxisMin, trackPhiAxisMax}, {passEMCalBins, passEMCalAxisMin, passEMCalAxisMax}}}},
     {"hTrackEnergyLossVsP", " TPC Energy loss info vs P; dE/dx;#it{p} (GeV#it{/c});passEMcal;", hTrackEnergyLossSpec},
     {"hTrackEnergyLossVsPt", " TPC Energy loss info vs Pt; dE/dx;#it{p}_{T} (GeV#it{/c});passEMcal;", hTrackEnergyLossSpec},
     {"hTracknSigmaVsP", " TPC nSigma info vs P; n#sigma;#it{p} (GeV#it{/c});passEMcal;", hTracknSigmaSpec},
     {"hTracknSigmaVsPt", " TPC nSigma info vs Pt; n#sigma;#it{p}_{T} (GeV#it{/c});passEMcal;", hTracknSigmaSpec},
     {"hEmcClusterEnergy", "EMCal Cluster Info before match Energy; Energy (GeV)", hEmcClusterEnergySpec},
     {"hEmcClusterEtaPhi", "EMCal Cluster Info before match Eta  and Phi; #eta;#varphi;", hEmcClusterEtaPhiSpec},
     {"hEmcClusterEnergyCell", "EMCal Cluster Info before match Energy vs nCells; Energy (GeV);ncell;", hEmcClusterEnergyCellSpec},
     {"hEmcClusterEnergyTime", "EMCal Cluster Info before match Energy vs time; Energy (GeV); sec;", hEmcClusterEnergyTimeSpec},
     {"hEmcClusterAfterMatchEnergy", "EMCal Cluster Info After match Energy; Energy (GeV)", hEmcClusterEnergySpec},
     {"hEmcClusterAfterMatchEtaPhi", "EMCal Cluster Info After match Eta  and Phi; #eta;#varphi;", hEmcClusterEtaPhiSpec},
     {"hEmcClusterAfterMatchEnergyCells", "EMCal Cluster Info After match Energy vs nCells; Energy (GeV);ncell;", hEmcClusterEnergyCellSpec},
     {"hEmcClusterAfterMatchEnergyTime", "EMCal Cluster Info After match Energy vs time; Energy (GeV); sec;", hEmcClusterEnergyTimeSpec},

     {"hAfterMatchSigmaVsEoP", "PID Info after  match EoP vs Sigma ; E/P;#it{p}_{T} (GeV#it{/c});n#sigma; m02; m20;", hAfterMatchEoPSigamSpec},
     {"hAfterMatchEoPVsP", "PID Info after match  EoP vs P; E/P;#it{p} (GeV#it{/c});", {HistType::kTH2F, {{eopAxisBins, eopAxisMin, eopAxisMax}, {pAxisBins, pAxisMin, pAxisMax}}}},
     {"hAfterMatchSigmaVsP", "PID Info after match Sigma vs Momentum ; n#sigma; #it{p} (GeV#it{/c}; ", {HistType::kTH2F, {{nSigmaAxisBins, nSigmaAxisMin, nSigmaAxisMax}, {pAxisBins, pAxisMin, pAxisMax}}}},
     {"hAfterMatchEtaPhi", "PID Info after match Eta vs Phi ; #eta; #varphi; ", {HistType::kTH2F, {{etaAxisBins, trackEtaAxisMin, trackEtaAxisMax}, {phiAxisBins, trackPhiAxisMin, trackPhiAxisMax}}}},
     {"hAfterMatchEnergyLossVsP", "PID Info after match Energy loss info vs P ; dE/dx;#it{p} (GeV#it{/c});; ", {HistType::kTH2F, {{dEdxAxisBins, dEdxAxisMin, dEdxAxisMax}, {pAxisBins, pAxisMin, pAxisMax}}}},
     {"hAfterMatchEnergyLossVsPt", "PID Info after match Energy loss info vs Pt ;dE/dx;#it{p}_{T} (GeV#it{/c}); ", {HistType::kTH2F, {{dEdxAxisBins, dEdxAxisMin, dEdxAxisMax}, {pAxisBins, pAxisMin, pAxisMax}}}},

     {"hAfterPIDEtaPhi", "PID Info after PID Cuts Eta vs Phi ; #eta; #varphi; ", {HistType::kTH2F, {{etaAxisBins, trackEtaAxisMin, trackEtaAxisMax}, {phiAxisBins, trackPhiAxisMin, trackPhiAxisMax}}}},
     {"hEPRatioAfterPID", "E/P Ratio after PID Cuts apply only trackwodca filter", {HistType::kTH2F, {{pAxisBins, pAxisMin, pAxisMax}, {300, 0, 30}}}},
     {"hPIDAfterPIDCuts", "PID Info after PID cuts; E/P;#it{p}_{T} (GeV#it{/c});n#sigma;m02; m20;", hAfterMatchEoPSigamSpec},
     {"hEmcClsTrkEtaPhiDiffTimeEnergy", "EmcClsTrkEtaPhiDiffTimeEnergy;#Delta#eta;#Delta#varphi;Sec;", hDeltaPhiDeltaEtaEmcClusterTrackSpecEnergy}}};

  void init(o2::framework::InitContext&)
  {
    registry.get<THnSparse>(HIST("hAfterMatchSigmaVsEoP"))->Sumw2();
    registry.get<THnSparse>(HIST("hPIDAfterPIDCuts"))->Sumw2();
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
        registry.fill(HIST("hEmcClusterEnergy"), emcClusterBefore.energy());                                // track etaphi infor after filter bit
        registry.fill(HIST("hEmcClusterEtaPhi"), emcClusterBefore.eta(), emcClusterBefore.phi());           // track etaphi infor after filter bit
        registry.fill(HIST("hEmcClusterEnergyCell"), emcClusterBefore.energy(), emcClusterBefore.nCells()); // track etaphi infor after filter bit
        registry.fill(HIST("hEmcClusterEnergyTime"), emcClusterBefore.energy(), emcClusterBefore.time());   // track etaphi infor after filter bit
        registry.fill(HIST("hEmcClusterM02"), emcClusterBefore.m02());
        registry.fill(HIST("hEmcClusterM20"), emcClusterBefore.m20());
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
    int electronId = -999;
    for (const auto& track : tracks) {
      phiTrack = track.phi();
      etaTrack = track.eta();
      pTrack = track.p();
      ptTrack = track.pt();
      dcaxyTrack = track.dcaXY();
      dcazTrack = track.dcaZ();
      tpcNsigmaTrack = track.tpcNSigmaEl();
      electronId = track.globalIndex();
      // Apply Track Selection
      if (!selTracks(track)) {
        continue;
      }
      passEMCal = 0;

      if ((phiTrack > phiTrackEMCalMin && phiTrack < phiTrackEMCalMax) && (etaTrack > etaTrackMin && etaTrack < etaTrackMax))
        passEMCal = 1; // EMcal acceptance passed
      if ((phiTrack > phiTrackDCalMin && phiTrack < phiTrackDCalMax) && ((etaTrack > etaTrackDCalPositiveMin && etaTrack < etaTrackDCalPositiveMax) || (etaTrack > etaTrackDCalNegativeMin && etaTrack < etaTrackDCalNegativeMax)))
        passEMCal = 2; // Dcal acceptance passed

      registry.fill(HIST("hTrackEtaPhi"), etaTrack, phiTrack, passEMCal);                 // track etaphi infor after filter bit
      registry.fill(HIST("hTrackEnergyLossVsP"), track.tpcSignal(), pTrack, passEMCal);   // track etaphi infor after filter bit
      registry.fill(HIST("hTrackEnergyLossVsPt"), track.tpcSignal(), ptTrack, passEMCal); // track etaphi infor after filter bit
      registry.fill(HIST("hTracknSigmaVsP"), tpcNsigmaTrack, pTrack, passEMCal);          // track etaphi infor after filter bit
      registry.fill(HIST("hTracknSigmaVsPt"), tpcNsigmaTrack, ptTrack, passEMCal);        // track etaphi infor after filter bit

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

        registry.fill(HIST("hEmcClsTrkEtaPhiDiffTimeEnergy"), deltaEtaMatch, deltaPhiMatch, timeEmcCluster);

        if (fillEmcClusterInfo)
          registry.fill(HIST("hEmcClusterAfterMatchEnergy"), emcCluster.energy());                         // track etaphi infor after filter bit
        registry.fill(HIST("hEmcClusterAfterMatchEtaPhi"), emcCluster.eta(), emcCluster.phi());            // track etaphi infor after filter bit
        registry.fill(HIST("hEmcClusterAfterMatchEnergyCells"), emcCluster.energy(), emcCluster.nCells()); // track etaphi infor after filter bit
        registry.fill(HIST("hEmcClusterAfterMatchEnergyTime"), emcCluster.energy(), emcCluster.time());    // track etaphi infor after filter bit

        eop = eMatchEmcCluster / pMatchTrack;

        registry.fill(HIST("hAfterMatchSigmaVsEoP"), eop, ptMatchTrack, tpcNsigmaMatchTrack, m02MatchEmcCluster, m20MatchEmcCluster);
        registry.fill(HIST("hAfterMatchEoPVsP"), eop, pMatchTrack);
        registry.fill(HIST("hAfterMatchSigmaVsP"), tpcNsigmaMatchTrack, pMatchTrack);
        registry.fill(HIST("hAfterMatchEtaPhi"), etaMatchTrack, phiMatchTrack);
        registry.fill(HIST("hAfterMatchEnergyLossVsP"), matchTrack.tpcSignal(), pMatchTrack);
        registry.fill(HIST("hAfterMatchEnergyLossVsPt"), matchTrack.tpcSignal(), ptMatchTrack);

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

        registry.fill(HIST("hPIDAfterPIDCuts"), eop, ptMatchTrack, tpcNsigmaMatchTrack, m02MatchEmcCluster, m20MatchEmcCluster);
        registry.fill(HIST("hEPRatioAfterPID"), pMatchTrack, eMatchEmcCluster);
        registry.fill(HIST("hAfterPIDEtaPhi"), etaMatchTrack, phiMatchTrack);
        if (eop < eopElectronMin || eop > eopElectronMax) {
          continue;
        }

        isEMcal = true;
        // std::cout << " electron id  in selection" << electronId << std::endl; // FIXME: Use LOG
        electronSel(matchTrack.collisionId(), electronId, etaMatchTrack, phiMatchTrack, ptMatchTrack, pMatchTrack, trackRapidity, matchTrack.dcaXY(), matchTrack.dcaZ(), matchTrack.tpcNSigmaEl(), matchTrack.tofNSigmaEl(),
                    eMatchEmcCluster, etaMatchEmcCluster, phiMatchEmcCluster, m02MatchEmcCluster, m20MatchEmcCluster, cellEmcCluster, timeEmcCluster, deltaEtaMatch, deltaPhiMatch, isEMcal);
      }

      /// Electron information without Emcal and use TPC and TOF
      if (isEMcal) {
        continue;
      }
      electronSel(track.collisionId(), electronId, etaTrack, phiTrack, ptTrack, pTrack, trackRapidity, dcaxyTrack, dcazTrack, track.tpcNSigmaEl(), track.tofNSigmaEl(),
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
