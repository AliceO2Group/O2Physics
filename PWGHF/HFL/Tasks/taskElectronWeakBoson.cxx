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

/// \file taskElectronWeakBoson.cxx
/// \brief task for WeakBoson (W/Z) based on electron in mid-rapidity
/// \author S. Sakai & S. Ito (Univ. of Tsukuba)
#include <vector>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"

#include "EMCALBase/Geometry.h"
#include "EMCALCalib/BadChannelMap.h"

#include "DataFormatsEMCAL/Cell.h"
#include "DataFormatsEMCAL/Constants.h"
#include "DataFormatsEMCAL/AnalysisCluster.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"

#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGHF/Core/HfHelper.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct HfTaskElectronWeakBoson {

  // configurable parameters
  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pt registry"};
  Configurable<float> binPtmax{"binPtmax", 100.0, "maximum pt registry"};
  Configurable<int> nBinsE{"nBinsE", 100, "N bins in E registry"};
  Configurable<float> binEmax{"binEmax", 100.0, "maximum E registry"};

  Configurable<float> vtxZ{"vtxZ", 10.f, ""};

  Configurable<float> etaTrLow{"etaTrLow", -0.6f, "minimun track eta"};
  Configurable<float> etaTrUp{"etaTrUp", 0.6f, "maximum track eta"};
  Configurable<float> dcaxyMax{"dcaxyMax", 2.0f, "mximum DCA xy"};
  Configurable<float> chi2ItsMax{"chi2ItsMax", 15.0f, "its chi2 cut"};
  Configurable<float> ptMin{"ptMin", 3.0f, "minimum pT cut"};
  Configurable<float> chi2TpcMax{"chi2TpcMax", 4.0f, "tpc chi2 cut"};
  Configurable<float> nclItsMin{"nclItsMin", 2.0f, "its # of cluster cut"};
  Configurable<float> nclTpcMin{"nclTpcMin", 100.0f, "tpc # if cluster cut"};
  Configurable<float> nclcrossTpcMin{"nclcrossTpcMin", 100.0f, "tpc # of crossedRows cut"};
  Configurable<float> nsigTpcMinLose{"nsigTpcMinLose", -3.0, "tpc Nsig lose lower cut"};
  Configurable<float> nsigTpcMin{"nsigTpcMin", -1.0, "tpc Nsig lower cut"};
  Configurable<float> nsigTpcMax{"nsigTpcMax", 3.0, "tpc Nsig upper cut"};

  Configurable<float> phiEmcMin{"phiEmcMin", 1.39, "EMC phi acc min"};
  Configurable<float> phiEmcMax{"phiEmcMax", 3.36, "EMC phi acc max"};
  Configurable<int> clusterDefinition{"clusterDefinition", 10, "cluster definition to be selected, e.g. 10=kV3Default"};
  Configurable<float> timeEmcMin{"timeEmcMin", -25., "Minimum EMCcluster timing"};
  Configurable<float> timeEmcMax{"timeEmcMax", +20., "Maximum EMCcluster timing"};
  Configurable<float> m02Min{"m02Min", 0.1, "Minimum M02"};
  Configurable<float> m02Max{"m02Max", 0.9, "Maximum M02"};
  Configurable<float> rMatchMax{"rMatchMax", 0.05, "cluster - track matching cut"};
  Configurable<float> eopMin{"eopMin", 0.9, "Minimum eop"};
  Configurable<float> eopMax{"eopMax", 1.3, "Maximum eop"};

  Configurable<float> rIsolation{"rIsolation", 0.3, "cone radius for isolation cut"};
  Configurable<float> energyIsolationMax{"energyIsolationMax", 0.1, "isolation cut on energy"};
  Configurable<int> trackIsolationMax{"trackIsolationMax", 3, "Maximum number of tracks in isolation cone"};

  struct HfElectronCandidate {
    float pt, eta, phi, energy;
    int charge;
    HfElectronCandidate(float ptr, float e, float ph, float en, int ch)
      : pt(ptr), eta(e), phi(ph), energy(en), charge(ch) {}

    int sign() const { return charge; }
  };
  std::vector<HfElectronCandidate> selectedElectronsIso;
  std::vector<HfElectronCandidate> selectedElectronsAss;

  using SelectedClusters = o2::aod::EMCALClusters;
  // PbPb
  using TrackEle = o2::soa::Join<o2::aod::Tracks, o2::aod::FullTracks, o2::aod::TracksExtra, o2::aod::TracksDCA, o2::aod::TrackSelection, o2::aod::pidTPCFullEl>;

  // pp
  // using TrackEle = o2::soa::Filtered<o2::soa::Join<o2::aod::Tracks, o2::aod::FullTracks, o2::aod::TracksDCA, o2::aod::TrackSelection, o2::aod::pidTPCEl, o2::aod::pidTOFEl>>;

  // Filter
  Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < vtxZ);

  Filter etafilter = (aod::track::eta < etaTrUp) && (aod::track::eta > etaTrLow);
  Filter dcaxyfilter = (nabs(aod::track::dcaXY) < dcaxyMax);
  Filter filterGlobalTr = requireGlobalTrackInFilter();

  Filter clusterDefinitionSelection = (o2::aod::emcalcluster::definition == clusterDefinition) && (o2::aod::emcalcluster::time >= timeEmcMin) && (o2::aod::emcalcluster::time <= timeEmcMax) && (o2::aod::emcalcluster::m02 > m02Min) && (o2::aod::emcalcluster::m02 < m02Max);

  // Data Handling Objects
  Preslice<o2::aod::EMCALClusterCells> perCluster = o2::aod::emcalclustercell::emcalclusterId;
  Preslice<o2::aod::EMCALAmbiguousClusterCells> perClusterAmb = o2::aod::emcalclustercell::emcalambiguousclusterId;
  PresliceUnsorted<o2::aod::EMCALMatchedTracks> perClusterMatchedTracks = o2::aod::emcalmatchedtrack::trackId;

  // Histogram registry: an object to hold your registrygrams
  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {

    // define axes you want to use
    const AxisSpec axisZvtx{400, -20, 20, "Zvtx"};
    const AxisSpec axisCounter{1, 0, 1, "events"};
    const AxisSpec axisEta{200, -1.0, 1.0, "#eta"};
    const AxisSpec axisPt{nBinsPt, 0, binPtmax, "p_{T}"};
    const AxisSpec axisNsigma{100, -5, 5, "N#sigma"};
    const AxisSpec axisE{nBinsE, 0, binEmax, "Energy"};
    const AxisSpec axisM02{100, 0, 1, "M02"};
    const AxisSpec axisdPhi{200, -1, 1, "dPhi"};
    const AxisSpec axisdEta{200, -1, 1, "dEta"};
    const AxisSpec axisPhi{350, 0, 7, "Phi"};
    const AxisSpec axisEop{200, 0, 2, "Eop"};
    const AxisSpec axisChi2{500, 0.0, 50.0, "#chi^{2}"};
    const AxisSpec axisCluster{100, 0.0, 200.0, "counts"};
    const AxisSpec axisITSNCls{20, 0.0, 20, "counts"};
    const AxisSpec axisEMCtime{200, -100.0, 100, "EMC time"};
    const AxisSpec axisIsoEnergy{100, 0, 1, "Isolation energy(GeV/C)"};
    const AxisSpec axisIsoTrack{20, -0.5, 19.5, "Isolation Track"};
    const AxisSpec axisInvMassZ{200, 0, 200, "M_{ee} (GeV/c^{2})"};
    const AxisSpec axisInvMassDy{200, 0, 2, "M_{ee} (GeV/c^{2})"};

    // create registrygrams
    registry.add("hZvtx", "Z vertex", kTH1F, {axisZvtx});
    registry.add("hEventCounter", "hEventCounter", kTH1F, {axisCounter});
    registry.add("hITSchi2", "ITS #chi^{2}", kTH1F, {axisChi2});
    registry.add("hTPCchi2", "TPC #chi^{2}", kTH1F, {axisChi2});
    registry.add("hTPCnCls", "TPC NCls", kTH1F, {axisCluster});
    registry.add("hITSnCls", "ITS NCls", kTH1F, {axisITSNCls});
    registry.add("hTPCnClsCrossedRows", "TPC NClsCrossedRows", kTH1F, {axisCluster});
    registry.add("hEta", "track eta", kTH1F, {axisEta});
    registry.add("hPt", "track pt", kTH1F, {axisPt});
    registry.add("hTPCNsigma", "TPC electron Nsigma", kTH2F, {{axisPt}, {axisNsigma}});
    registry.add("hEnergy", "EMC cluster energy", kTH1F, {axisE});
    registry.add("hM02", "EMC M02", kTH2F, {{axisNsigma}, {axisM02}});
    registry.add("hM20", "EMC M20", kTH2F, {{axisNsigma}, {axisM02}});
    registry.add("hTrMatch", "Track EMC Match", kTH2F, {{axisdPhi}, {axisdEta}});
    registry.add("hTrMatch_mim", "Track EMC Match minimu minimumm", kTH2F, {{axisdPhi}, {axisdEta}});
    registry.add("hMatchPhi", "Match in Phi", kTH2F, {{axisPhi}, {axisPhi}});
    registry.add("hMatchEta", "Match in Eta", kTH2F, {{axisEta}, {axisEta}});
    registry.add("hEop", "energy momentum match", kTH2F, {{axisPt}, {axisEop}});
    registry.add("hEopIsolation", "energy momentum match after isolation", kTH2F, {{axisPt}, {axisEop}});
    registry.add("hEopIsolationTr", "energy momentum match after isolationTr", kTH2F, {{axisPt}, {axisEop}});
    registry.add("hEopNsigTPC", "Eop vs. Nsigma", kTH2F, {{axisNsigma}, {axisEop}});
    registry.add("hEMCtime", "EMC timing", kTH1F, {axisEMCtime});
    registry.add("hIsolationEnergy", "Isolation Energy", kTH2F, {{axisE}, {axisIsoEnergy}});
    registry.add("hIsolationTrack", "Isolation Track", kTH2F, {{axisE}, {axisIsoTrack}});
    registry.add("hInvMassZeeLs", "invariant mass for Z LS pair", kTH2F, {{axisPt}, {axisInvMassZ}});
    registry.add("hInvMassZeeUls", "invariant mass for Z ULS pair", kTH2F, {{axisPt}, {axisInvMassZ}});
    registry.add("hInvMassDyLs", "invariant mass for DY LS pair", kTH2F, {{axisPt}, {axisInvMassDy}});
    registry.add("hInvMassDyUls", "invariant mass for DY ULS pair", kTH2F, {{axisPt}, {axisInvMassDy}});
  }
  bool isIsolatedCluster(const o2::aod::EMCALCluster& cluster,
                         const SelectedClusters& clusters)
  {
    float energySum = 0.0;
    float isoEnergy = 10.0;
    float etaAssCluster = cluster.eta();
    float phiAssCluster = cluster.phi();

    for (const auto& associateCluster : clusters) {
      // Calculate angular distances
      double dEta = associateCluster.eta() - etaAssCluster;
      double dPhi = associateCluster.phi() - phiAssCluster;

      // Normalize φ difference
      dPhi = RecoDecay::constrainAngle(dPhi, -o2::constants::math::PI);

      // Calculate ΔR
      double deltaR = std::sqrt(dEta * dEta + dPhi * dPhi);

      // Sum energy within isolation cone
      if (deltaR < rIsolation) {
        energySum += associateCluster.energy();
      }
    }

    if (energySum > 0) {
      isoEnergy = energySum / cluster.energy() - 1.0;
    }

    registry.fill(HIST("hIsolationEnergy"), cluster.energy(), isoEnergy);

    return (isoEnergy < energyIsolationMax);
  }
  bool isIsolatedTrack(double etaEle,
                       double phiEle,
                       float ptEle,
                       TrackEle const& tracks)
  {
    int trackCount = 0;

    for (const auto& track : tracks) {
      // skip the reference track
      if (std::abs(track.pt() - ptEle) < 1e-4)
        continue;

      double dEta = track.eta() - etaEle;
      double dPhi = track.phi() - phiEle;
      dPhi = RecoDecay::constrainAngle(dPhi, -o2::constants::math::PI);

      double deltaR = std::sqrt(dEta * dEta + dPhi * dPhi);

      if (deltaR < rIsolation) {
        trackCount++;
      }
    }

    registry.fill(HIST("hIsolationTrack"), ptEle, trackCount);

    return (trackCount <= trackIsolationMax);
  }

  void process(soa::Filtered<aod::Collisions>::iterator const& collision,
               SelectedClusters const& emcClusters,
               TrackEle const& tracks,
               o2::aod::EMCALMatchedTracks const& matchedtracks)
  {
    registry.fill(HIST("hEventCounter"), 0.5);

    // LOGF(info, "Collision index : %d", collision.index());
    // LOGF(info, "Number of tracks: %d", tracks.size());
    // LOGF(info, "Number of clusters: %d", clusters.size());

    registry.fill(HIST("hZvtx"), collision.posZ());

    for (const auto& track : tracks) {

      if (std::abs(track.eta()) > etaTrUp)
        continue;
      if (track.tpcNClsCrossedRows() < nclcrossTpcMin)
        continue;
      if (std::abs(track.dcaXY()) > dcaxyMax)
        continue;
      if (track.itsChi2NCl() > chi2ItsMax)
        continue;
      if (track.tpcChi2NCl() > chi2TpcMax)
        continue;
      if (track.tpcNClsFound() < nclTpcMin)
        continue;
      if (track.itsNCls() < nclItsMin)
        continue;
      if (track.pt() < ptMin)
        continue;

      registry.fill(HIST("hEta"), track.eta());
      registry.fill(HIST("hITSchi2"), track.itsChi2NCl());
      registry.fill(HIST("hTPCchi2"), track.tpcChi2NCl());
      registry.fill(HIST("hTPCnCls"), track.tpcNClsFound());
      registry.fill(HIST("hITSnCls"), track.itsNCls());
      registry.fill(HIST("hTPCnClsCrossedRows"), track.tpcNClsCrossedRows());
      registry.fill(HIST("hPt"), track.pt());
      registry.fill(HIST("hTPCNsigma"), track.p(), track.tpcNSigmaEl());

      float energyTrk = 0.0;

      if (track.tpcNSigmaEl() > nsigTpcMinLose && track.tpcNSigmaEl() < nsigTpcMax) {
        selectedElectronsAss.emplace_back(
          track.pt(),
          track.eta(),
          track.phi(),
          energyTrk,
          track.sign());
      }

      // track - match

      //  continue;
      if (track.phi() < phiEmcMin || track.phi() > phiEmcMax)
        continue;
      auto tracksofcluster = matchedtracks.sliceBy(perClusterMatchedTracks, track.globalIndex());

      // LOGF(info, "Number of matched track: %d", tracksofcluster.size());

      double rMin = 999.9;
      double dPhiMin = 999.9;
      double dEtaMin = 999.9;
      bool isIsolated = false;
      bool isIsolatedTr = false;

      if (tracksofcluster.size()) {
        int nMatch = 0;
        for (const auto& match : tracksofcluster) {
          if (match.emcalcluster_as<SelectedClusters>().time() < timeEmcMin || match.emcalcluster_as<SelectedClusters>().time() > timeEmcMax)
            continue;
          if (match.emcalcluster_as<SelectedClusters>().m02() < m02Min || match.emcalcluster_as<SelectedClusters>().m02() > m02Max)
            continue;

          float m20Emc = match.emcalcluster_as<SelectedClusters>().m20();
          float m02Emc = match.emcalcluster_as<SelectedClusters>().m02();
          float energyEmc = match.emcalcluster_as<SelectedClusters>().energy();
          double phiEmc = match.emcalcluster_as<SelectedClusters>().phi();
          double etaEmc = match.emcalcluster_as<SelectedClusters>().eta();
          double timeEmc = match.emcalcluster_as<SelectedClusters>().time();
          // LOG(info) << "tr phi0 = " << match.track_as<TrackEle>().phi();
          // LOG(info) << "tr phi1 = " << track.phi();
          // LOG(info) << "emc phi = " << phiEmc;

          if (nMatch == 0) {
            double dEta = match.track_as<TrackEle>().trackEtaEmcal() - etaEmc;
            double dPhi = match.track_as<TrackEle>().trackPhiEmcal() - phiEmc;
            dPhi = RecoDecay::constrainAngle(dPhi, -o2::constants::math::PI);

            registry.fill(HIST("hMatchPhi"), phiEmc, match.track_as<TrackEle>().trackPhiEmcal());
            registry.fill(HIST("hMatchEta"), etaEmc, match.track_as<TrackEle>().trackEtaEmcal());

            double r = RecoDecay::sqrtSumOfSquares(dPhi, dEta);
            if (r < rMin) {
              rMin = r;
              dPhiMin = dPhi;
              dEtaMin = dEta;
            }
            registry.fill(HIST("hTrMatch"), dPhi, dEta);
            registry.fill(HIST("hEMCtime"), timeEmc);
            registry.fill(HIST("hEnergy"), energyEmc);

            if (r > rMatchMax)
              continue;

            const auto& cluster = match.emcalcluster_as<SelectedClusters>();

            double eop = energyEmc / match.track_as<TrackEle>().p();

            // LOG(info) << "E/p" << eop;
            registry.fill(HIST("hEopNsigTPC"), match.track_as<TrackEle>().tpcNSigmaEl(), eop);
            registry.fill(HIST("hM02"), match.track_as<TrackEle>().tpcNSigmaEl(), m02Emc);
            registry.fill(HIST("hM20"), match.track_as<TrackEle>().tpcNSigmaEl(), m20Emc);
            if (match.track_as<TrackEle>().tpcNSigmaEl() > nsigTpcMin && match.track_as<TrackEle>().tpcNSigmaEl() < nsigTpcMax) {
              registry.fill(HIST("hEop"), match.track_as<TrackEle>().pt(), eop);

              if (eop > eopMin && eop < eopMax) {
                isIsolated = isIsolatedCluster(cluster, emcClusters);
                isIsolatedTr = isIsolatedTrack(track.phi(), track.eta(), track.pt(), tracks);
              }

              if (isIsolated) {
                registry.fill(HIST("hEopIsolation"), match.track_as<TrackEle>().pt(), eop);

                selectedElectronsIso.emplace_back(
                  match.track_as<TrackEle>().pt(),
                  match.track_as<TrackEle>().eta(),
                  match.track_as<TrackEle>().phi(),
                  energyEmc,
                  match.track_as<TrackEle>().sign());
              }

              if (isIsolatedTr) {
                registry.fill(HIST("hEopIsolationTr"), match.track_as<TrackEle>().pt(), eop);
              }
            }
          }

          nMatch++;
        }
      }

      if (rMin < rMatchMax) {
        // LOG(info) << "R mim = " << rMin;
        registry.fill(HIST("hTrMatch_mim"), dPhiMin, dEtaMin);
      }

    } // end of track loop

    // calculate inv. mass
    if (selectedElectronsIso.size() > 1) {
      for (size_t i = 0; i < selectedElectronsIso.size(); i++) {
        const auto& e1 = selectedElectronsIso[i];
        for (size_t j = 0; j < selectedElectronsAss.size(); j++) {
          const auto& e2 = selectedElectronsAss[j];

          float ptIso = e1.pt;
          float ptAss = e2.pt;
          if (ptIso == ptAss)
            continue;
          auto arr1 = RecoDecayPtEtaPhi::pVector(e1.pt, e1.eta, e1.phi);
          auto arr2 = RecoDecayPtEtaPhi::pVector(e2.pt, e2.eta, e2.phi);
          double mass = RecoDecay::m(std::array{arr1, arr2}, std::array{o2::constants::physics::MassElectron, o2::constants::physics::MassElectron});
          if (e1.sign() * e2.sign() > 0) {
            registry.fill(HIST("hInvMassDyLs"), ptIso, mass);
          } else {
            registry.fill(HIST("hInvMassDyUls"), ptIso, mass);
          }

          if (ptAss < 20.0 && ptIso < 20.0)
            continue;

          if (e1.sign() * e2.sign() > 0) {
            registry.fill(HIST("hInvMassZeeLs"), ptIso, mass);
          } else {
            registry.fill(HIST("hInvMassZeeUls"), ptIso, mass);
          }
        }
      }
    } // end of inv. mass calculation
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskElectronWeakBoson>(cfgc)};
}
