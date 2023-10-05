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

/// \file taskCorrelationHfeHadrons.cxx
/// \author Rashi Gupta <rashi.gupta@cern.ch>, IIT Indore
/// \author Ravindra Singh <ravindra.singh@cern.ch>, IIT Indore

#include "THnSparse.h"

#include "DataFormatsEMCAL/AnalysisCluster.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/PID/TPCPIDResponse.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGJE/DataModel/EMCALClusters.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct HfTaskCorrelationHfeHadrons {
  // Configurables

  // Cluster information
  Configurable<bool> clusterInfo{"clusterInfo", true, "EMCal cluster info before and after track match"};
  // Run3 selection
  Configurable<bool> isRun3{"isRun3", true, "Data is from Run3 or Run2"};
  // Electron hadron correlation condition
  Configurable<bool> ptCondition{"ptCondition", true, "Electron pT should be greater than associate particle pT"};

  // Deltaphi binning
  Configurable<int> nBinsDeltaPhi{"nBinsDeltaPhi", 32, "Bins for #Delta#varphi bins"};

  // Track selection
  Configurable<float> dcaXYTrackMax{"dcaXYTrackMax", 0.5f, "DCA XY cut"};
  Configurable<float> dcaZTrackMax{"dcaZTrackMax", 1.0f, "DCA Z cut"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.6f, "Eta range for electron tracks"};
  Configurable<float> etaTrackMin{"etaTrackMin", -0.6f, "Eta range for electron tracks"};
  Configurable<float> ptTrackMin{"ptTrackMin", 2.0f, "Transverse MOmentum range for electron tracks"};

  // EMcal and Dcal selection cut
  Configurable<float> etaTrackDcalLeftMax{"etaTrackDcalLeftMax", -0.22f, "Eta range for electron Dcal tracks"};
  Configurable<float> etaTrackDcalLeftMin{"etaTrackDcalLeftMin", -0.6f, "Eta range for electron tracks"};
  Configurable<float> etaTrackDcalRightMax{"etaTrackDcalRightMax", 0.6f, "Eta range for electron Dcal tracks"};
  Configurable<float> etaTrackDcalRightMin{"etaTrackDcalRightMin", 0.22f, "Eta range for electron tracks"};
  Configurable<float> phiTrackDCalMax{"phiTrackDCalMax", 3.3621f, "phi range for electron tracks associated Dcal"};
  Configurable<float> phiTrackDcalMin{"phiTrackDcalMin", 1.3955f, "phi range for electron tracks associated Dcal"};
  Configurable<float> phiTrackEMcalMax{"phiTrackEMcalMax", 5.708f, "phi range for electron tracks associated Emcal"};
  Configurable<float> phiTrackEMcalMin{"phiTrackEMcalMin", 4.5355f, "phi range for electron tracks associated Emcal"};

  // Track and Cluster matching cut
  Configurable<float> deltaEtaMatchMin{"deltaEtaMatchMin", 0.01f, "Eta distance of EMCAL cluster to its closest track"};
  Configurable<float> deltaPhiMatchMin{"deltaPhiMatchMin", 0.01f, "Phi distance of EMCAL cluster to its closest track"};
  Configurable<float> timeClusterMax{"timeClusterMax", 100.f, "Cluster time"};

  // Inclusive electron selection cut
  Configurable<float> eopElMin{"eopElMin", 0.8f, "Minimum E/p for electron tracks"};
  Configurable<float> eopElMax{"eopElMax", 1.2f, "Maximum E/p for electron tracks"};
  Configurable<float> m02ElMax{"m02ElMax", 0.9f, "max Electron M02"};
  Configurable<float> m02ElMin{"m02ElMin", 0.02f, "min Electron M02"};
  Configurable<float> m20ElMax{"m20ElMax", 1000.f, "max Electron M20"};
  Configurable<float> m20ElMin{"m20ElMin", 0.0f, "min Electron M20"};
  Configurable<float> nSigmaTpcElMin{"nSigmaTpcElMin", -1.0f, "min Electron TPCnsigma"};
  Configurable<float> nSigmaTpcElMax{"nSigmaTpcElMax", 3.0f, "max Electron TPCnsigma"};

  // Hadron selection cut
  Configurable<float> etaTrackAssocMax{"etaTrackAssocMax", 0.8f, "Eta range for Associated or partner electron tracks"};
  Configurable<float> ptTrackAssocMin{"ptTrackAssocMin", 0.2f, "Minimum pT for associated track"};

  using CollisionTable = soa::Filtered<soa::Join<aod::Collisions, aod::Mults, aod::EvSels>>::iterator;
  using TrackTables = o2::soa::Join<o2::aod::Tracks, o2::aod::TracksExtra, o2::aod::pidTPCFullEl, o2::aod::TracksDCA, o2::aod::TrackSelection, o2::aod::TrackSelectionExtension>;
  using TrackTable = TrackTables::iterator;

  Filter CollisionFilter = nabs(aod::collision::posZ) < 10.f && aod::collision::numContrib > (uint16_t)1;

  PresliceUnsorted<o2::aod::EMCALMatchedTracks> perClusterMatchedTracks = o2::aod::emcalmatchedtrack::trackId;

  HistogramConfigSpec hClusterInfoSpec{HistType::kTHnSparseD, {{300, 0.0, 30.0}, {100, -0.9, 0.9}, {200, 0, 6.3}, {50, 0, 50}, {1800, -900, 900}}};
  HistogramConfigSpec hCorrelSpec{HistType::kTHnSparseD, {{30, 0., 30.}, {20, 0., 20.}, {nBinsDeltaPhi, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}, {50, -1.8, 1.8}}};
  HistogramConfigSpec hDeltaPhiDeltaEtaClusterTrackSpec{HistType::kTH3F, {{400, -0.2, 0.2}, {400, -0.2, 0.2}, {280, 0, 70}}};
  HistogramConfigSpec hPIDSpec{HistType::kTHnSparseD, {{500, 0.0, 50.0}, {500, 0., 50.}, {300, -15, 15}, {300, 0.0, 30.0}, {400, 0, 2}, {400, 0, 2}}};
  HistogramConfigSpec hTrackAllInfoSpec{HistType::kTHnSparseD, {{500, 0, 160}, {300, -15, 15}, {500, 0., 50.}, {500, 0., 50.}, {100, -1.5, 1.5}, {100, 0, 7}, {2000, -1, 1}, {2000, -1, 1}, {3, 0, 3}}};
  HistogramConfigSpec hTrackInfoSpec{HistType::kTHnSparseD, {{500, 0, 160}, {300, -15, 15}, {500, 0., 50.}, {500, 0., 50.}, {100, -1.5, 1.5}, {100, 0, 7}}};

  HistogramRegistry registry{
    "registry",
    {{"hNevents", "No of events", {HistType::kTH1F, {{3, 1, 4}}}},
     {"hZvertex", "z vertex", {HistType::kTH1F, {{100, -20, 20}}}},

     {"hHHCorrel", "Sparse for Delta phi and Delta eta with Hadron;p_{T}^{e} (GeV#it{/c});p_{T}^{h} (GeV#it{/c});#Delta#varphi;#Delta#eta;", hCorrelSpec},
     {"hInclusiveEHCorrel", "Sparse for Delta phi and Delta eta with Inclusive electron;p_{T}^{e} (GeV#it{/c});p_{T}^{h} (GeV#it{/c});#Delta#varphi;#Delta#eta;", hCorrelSpec},
     {"hTrackInformation", "Sparse TPC info; dE/dx;n#sigma;#it{p} (GeV#it{/c});#it{p}_{T} (GeV#it{/c});#eta;#varphi; DcaXY;Dcaz;passEMcal;", hTrackAllInfoSpec},
     {"hHadronInformation", "Sparse hadron info; dE/dx;n#sigma;#it{p} (GeV#it{/c});#it{p}_{T} (GeV#it{/c});#eta;#varphi;", hTrackInfoSpec},
     {"hClusterInformationBefore", "Cluster Info before match; Energy (GeV);#eta;#varphi", hClusterInfoSpec},
     {"hClusterInformationAfter", "Cluster Info after match; Energy (GeV);#eta;#varphi", hClusterInfoSpec},
     {"hPIDafterMatch", "PID Info after match;dE/dx;n#sigma;#it{p} (GeV#it{/c});#it{p}_{T} (GeV#it{/c});#eta;#varphi;", hTrackInfoSpec},
     {"hPIDafterPIDcuts", "PID Info after PID cuts; #it{p} (GeV#it{/c});#it{p}_{T} (GeV#it{/c});n_{#sigma}^{e};GeV;M02;M20", hPIDSpec},
     {"hClsTrkEtaPhiDiffTime", "ClsTrkEtaPhiDiffTime;#Delta#eta;#Delta#varphi;Sec;", hDeltaPhiDeltaEtaClusterTrackSpec}}};

  void init(InitContext&)
  {
    registry.get<THnSparse>(HIST("hInclusiveEHCorrel"))->Sumw2();
    registry.get<THnSparse>(HIST("hHHCorrel"))->Sumw2();
    registry.get<THnSparse>(HIST("hTrackInformation"))->Sumw2();
    registry.get<THnSparse>(HIST("hClusterInformationBefore"))->Sumw2();
    registry.get<THnSparse>(HIST("hClusterInformationAfter"))->Sumw2();
    registry.get<THnSparse>(HIST("hPIDafterMatch"))->Sumw2();
    registry.get<THnSparse>(HIST("hPIDafterPIDcuts"))->Sumw2();
    registry.get<THnSparse>(HIST("hHadronInformation"))->Sumw2();
  }

  // correlation function for electron hadron

  void correlation(TrackTable const& eTrack, TrackTables const& assoTrack, int sparseNo = -1)
  {
    if (sparseNo < 0 || sparseNo > 1) {
      LOGF(info, "Error: no matching function for call.");
      return;
    }
    std::shared_ptr<THnSparse> hEHCorrArray[2] = {
      registry.get<THnSparse>(HIST("hHHCorrel")),
      registry.get<THnSparse>(HIST("hInclusiveEHCorrel"))};

    // Construct Deta Phi between electrons and hadrons
    double deltaPhi = -999;
    double deltaEta = -999;
    double ptHadron = -999;
    double ptElectron = -999;
    double phiElectron = -999;
    double phiHadron = -999;
    double etaElectron = -999;
    double etaHadron = -999;
    for (const auto& hTrack : assoTrack) {
      if (hTrack.globalIndex() == eTrack.globalIndex())
        continue;
      ptHadron = hTrack.pt();
      ptElectron = eTrack.pt();
      phiElectron = eTrack.phi();
      phiHadron = hTrack.phi();
      etaElectron = eTrack.eta();
      etaHadron = hTrack.eta();

      // Apply Hadron cuts
      if (std::abs(etaHadron) > etaTrackAssocMax)
        continue;
      if (ptHadron < ptTrackAssocMin)
        continue;
      if (std::abs(hTrack.dcaXY()) > dcaXYTrackMax || std::abs(eTrack.dcaZ()) > dcaZTrackMax)
        continue;
      if (!hTrack.isGlobalTrackWoDCA())
        continue;

      registry.fill(HIST("hHadronInformation"), hTrack.tpcSignal(), hTrack.tpcNSigmaEl(), hTrack.p(), ptHadron, etaHadron, phiHadron);

      if (ptCondition && (ptElectron > ptHadron))
        continue;

      deltaPhi = RecoDecay::constrainAngle(phiElectron - phiHadron, -o2::constants::math::PIHalf);
      deltaEta = etaElectron - etaHadron;
      if (sparseNo >= 0)
        hEHCorrArray[sparseNo]->Fill(ptElectron, ptHadron, deltaPhi, deltaEta);
    }
  }

  void process(CollisionTable const& collision, aod::EMCALClusters const& mAnalysisClusters, o2::aod::EMCALMatchedTracks const& matchedTracks, TrackTables const& tracks)
  {
    if (!(isRun3 ? collision.sel8() : (collision.sel7() && collision.alias_bit(kINT7))))
      return;

    registry.fill(HIST("hNevents"), 1);
    registry.fill(HIST("hZvertex"), collision.posZ());

    /////////////////////////////////
    // cluster info before match ///
    ///////////////////////////////
    if (clusterInfo) {
      for (const auto& clusterbf : mAnalysisClusters) {
        registry.fill(HIST("hClusterInformationBefore"), clusterbf.energy(), clusterbf.eta(), clusterbf.phi(), clusterbf.nCells(), clusterbf.time());
      }
    }
    int passEmcal;
    double phiTrack = -999;
    double etaTrack = -999;
    double pTrack = -999;
    double ptTrack = -999;
    double dcaxyTrack = -999;
    double dcazTrack = -999;
    double nSigmaTpcTrack = -999;

    for (const auto& track : tracks) {

      phiTrack = track.phi();
      etaTrack = track.eta();
      pTrack = track.p();
      ptTrack = track.pt();
      dcaxyTrack = track.dcaXY();
      dcazTrack = track.dcaZ();
      nSigmaTpcTrack = track.tpcNSigmaEl();

      if (!track.isGlobalTrackWoDCA())
        continue;
      passEmcal = 0;
      if ((phiTrack > phiTrackEMcalMin && phiTrack < phiTrackEMcalMax) && (etaTrack > etaTrackMin && etaTrack < etaTrackMax))
        passEmcal = 1; // EMcal acceptance passed
      if ((phiTrack > phiTrackDcalMin && phiTrack < phiTrackDCalMax) && ((etaTrack > etaTrackDcalRightMin && etaTrack < etaTrackDcalRightMax) || (etaTrack > etaTrackDcalLeftMin && etaTrack < etaTrackDcalLeftMax)))
        passEmcal = 2;                                                                                                                                    // Dcal acceptance passed
      registry.fill(HIST("hTrackInformation"), track.tpcSignal(), nSigmaTpcTrack, pTrack, ptTrack, etaTrack, phiTrack, dcaxyTrack, dcazTrack, passEmcal); // track infor after filter bit

      // Apply Track cut
      if (ptTrack < ptTrackMin)
        continue;
      if (std::abs(dcaxyTrack) > dcaXYTrackMax || std::abs(dcazTrack) > dcaZTrackMax)
        continue;

      auto tracksofcluster = matchedTracks.sliceBy(perClusterMatchedTracks, track.globalIndex());
      double phiMatchTrack = -999;
      double etaMatchTrack = -999;
      double pMatchTrack = -999;
      double ptMatchTrack = -999;
      double nSigmaTpcMatchTrack = -999;
      double phiMatchCluster = -999;
      double etaMatchCluster = -999;
      double eMatchCluster = -999;
      double m02MatchCluster = -999;
      double m20MatchCluster = -999;
      double timeMatchCluster = -999;
      for (const auto& emtrack : tracksofcluster) {

        if (track.globalIndex() != emtrack.trackId())
          continue;
        auto mtrack = emtrack.track_as<TrackTables>();
        auto cluster = emtrack.emcalcluster_as<aod::EMCALClusters>();

        phiMatchTrack = mtrack.phi();
        etaMatchTrack = mtrack.eta();
        pMatchTrack = mtrack.p();
        ptMatchTrack = mtrack.pt();
        nSigmaTpcMatchTrack = mtrack.tpcNSigmaEl();
        phiMatchCluster = cluster.phi();
        etaMatchCluster = cluster.eta();
        eMatchCluster = cluster.energy();
        m02MatchCluster = cluster.m02();
        m20MatchCluster = cluster.m20();
        timeMatchCluster = cluster.time();
        correlation(mtrack, tracks, 0); //"0" stands for filling Di-hadron

        if (etaMatchTrack < etaTrackMin || etaMatchTrack > etaTrackMax)
          continue;
        double deltaPhiMatch = -999.;
        double deltaEtaMatch = -999.;

        deltaPhiMatch = mtrack.trackPhiEmcal() - phiMatchCluster;
        deltaEtaMatch = mtrack.trackEtaEmcal() - etaMatchCluster;

        registry.fill(HIST("hClsTrkEtaPhiDiffTime"), deltaEtaMatch, deltaPhiMatch, timeMatchCluster);

        // Track and cluster Matching
        if (std::abs(deltaPhiMatch) > deltaPhiMatchMin || std::abs(deltaEtaMatch) > deltaEtaMatchMin)
          continue;
        if (timeMatchCluster > timeClusterMax)
          continue;

        if (clusterInfo)
          registry.fill(HIST("hClusterInformationAfter"), eMatchCluster, etaMatchCluster, phiMatchCluster, cluster.nCells(), timeMatchCluster);

        registry.fill(HIST("hPIDafterMatch"), mtrack.tpcSignal(), nSigmaTpcMatchTrack, pMatchTrack, ptMatchTrack, etaMatchTrack, phiMatchTrack);

        double eop = eMatchCluster / pMatchTrack;

        // Apply Electron Identification cuts

        if ((nSigmaTpcMatchTrack < nSigmaTpcElMin || nSigmaTpcMatchTrack > nSigmaTpcElMax) || (eop < eopElMin || eop > eopElMax) || (m02MatchCluster < m02ElMin || m02MatchCluster > m02ElMax) || (m20MatchCluster < m20ElMin || m20MatchCluster > m20ElMax))
          continue;
        registry.fill(HIST("hPIDafterPIDcuts"), pMatchTrack, ptMatchTrack, nSigmaTpcMatchTrack, eMatchCluster, m02MatchCluster, m20MatchCluster);

        correlation(mtrack, tracks, 1); //"1" stands for filling Electron-hadron
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskCorrelationHfeHadrons>(cfgc)};
}
