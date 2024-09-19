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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"

#include "EMCALBase/Geometry.h"
#include "EMCALCalib/BadChannelMap.h"

#include "DataFormatsEMCAL/Cell.h"
#include "DataFormatsEMCAL/Constants.h"
#include "DataFormatsEMCAL/AnalysisCluster.h"

#include "PWGJE/DataModel/EMCALClusters.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct HfTaskElectronWeakBoson {

  using SelectedClusters = o2::aod::EMCALClusters;

  // PbPb
  using TrackEle = o2::soa::Join<o2::aod::Tracks, o2::aod::FullTracks, o2::aod::TracksExtra_001, o2::aod::TracksDCA, o2::aod::TrackSelection, o2::aod::pidTPCFullEl>;

  // pp
  // using TrackEle = o2::soa::Filtered<o2::soa::Join<o2::aod::Tracks, o2::aod::FullTracks, o2::aod::TracksDCA, o2::aod::TrackSelection, o2::aod::pidTPCEl, o2::aod::pidTOFEl>>;

  // configurable parameters
  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pt registry"};
  Configurable<float> BinPtmax{"BinPtmax", 100.0, "maximum pt registry"};
  Configurable<int> nBinsE{"nBinsE", 100, "N bins in E registry"};
  Configurable<float> BinEmax{"BinEmax", 100.0, "maximum E registry"};

  Configurable<float> vtxZ{"vtxZ", 10.f, ""};

  Configurable<float> etalowCut{"etalowCut", -0.6f, "eta lower cut"};
  Configurable<float> etaupCut{"etaupCut", 0.6f, "eta upper cut"};
  Configurable<float> dcaxyCut{"dcaxyCut", 2.0f, "dca xy cut"};
  Configurable<float> itschi2Cut{"itschi2Cut", 15.0f, "its chi2 cut"};
  Configurable<float> mimpTCut{"mimpTCut", 3.0f, "minimum pT cut"};
  Configurable<float> tpcchi2Cut{"tpcchi2Cut", 4.0f, "tpc chi2 cut"};
  Configurable<float> itsNclCut{"itsNclCut", 2.0f, "its # of cluster cut"};
  Configurable<float> tpcNclCut{"tpcNclCut", 100.0f, "tpc # if cluster cut"};
  Configurable<float> tpcNclCrCut{"tpcNclCrCut", 100.0f, "tpc # of crossedRows cut"};
  Configurable<float> tpcNsiglowCut{"tpcNsiglowCut", -1.0, "tpc Nsig lower cut"};
  Configurable<float> tpcNsigupCut{"tpcNsigupCut", 3.0, "tpc Nsig upper cut"};

  Configurable<float> emcaccPhimin{"emcaccPhimin", 1.39, "Maximum M20"};
  Configurable<float> emcaccPhimax{"emcaccPhimax", 3.36, "Maximum M20"};
  Configurable<int> ClusterDefinition{"ClusterDefinition", 10, "cluster definition to be selected, e.g. 10=kV3Default"};
  Configurable<float> minTime{"minTime", -25., "Minimum cluster time"};
  Configurable<float> maxTime{"maxTime", +20., "Maximum cluster time"};
  Configurable<float> minM02{"minM02", 0.1, "Minimum M02"};
  Configurable<float> maxM02{"maxM02", 0.9, "Maximum M02"};
  Configurable<float> minM20{"minM20", 0.1, "Minimum M20"};
  Configurable<float> maxM20{"maxM20", 0.6, "Maximum M20"};
  Configurable<float> MatchRCut{"MatchRCut", 0.1, "cluster - track matching cut"};

  // Filter
  Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < vtxZ);

  Filter etafilter = (aod::track::eta < etaupCut) && (aod::track::eta > etalowCut);
  Filter dcaxyfilter = (nabs(aod::track::dcaXY) < dcaxyCut);
  Filter filter_globalTr = requireGlobalTrackInFilter();

  Filter clusterDefinitionSelection = (o2::aod::emcalcluster::definition == ClusterDefinition) && (o2::aod::emcalcluster::time >= minTime) && (o2::aod::emcalcluster::time <= maxTime) && (o2::aod::emcalcluster::m02 > minM02) && (o2::aod::emcalcluster::m02 < maxM02);

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
    const AxisSpec axisPt{nBinsPt, 0, BinPtmax, "p_{T}"};
    const AxisSpec axisNsigma{100, -5, 5, "N#sigma"};
    const AxisSpec axisE{nBinsE, 0, BinEmax, "Energy"};
    const AxisSpec axisM02{100, 0, 1, "M02"};
    const AxisSpec axisdPhi{200, -1, 1, "dPhi"};
    const AxisSpec axisdEta{200, -1, 1, "dEta"};
    const AxisSpec axisPhi{350, 0, 7, "Phi"};
    const AxisSpec axisEop{200, 0, 2, "Eop"};
    const AxisSpec axisChi2{500, 0.0, 50.0, "#chi^{2}"};
    const AxisSpec axisCluster{100, 0.0, 200.0, "counts"};
    const AxisSpec axisITSNCls{20, 0.0, 20, "counts"};
    const AxisSpec axisEMCtime{200, -100.0, 100, "EMC time"};

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
    registry.add("hEopNsigTPC", "Eop vs. Nsigma", kTH2F, {{axisNsigma}, {axisEop}});
    registry.add("hEMCtime", "EMC timing", kTH1F, {axisEMCtime});
  }

  // void process(soa::Filtered<aod::Collisions>::iterator const& collision, SelectedClusters const& clusters, TrackEle const& tracks, o2::aod::EMCALMatchedTracks const& matchedtracks)
  void process(soa::Filtered<aod::Collisions>::iterator const& collision, SelectedClusters const& emcClusters, TrackEle const& tracks, o2::aod::EMCALMatchedTracks const& matchedtracks)
  // void process(soa::Filtered<aod::Collisions>::iterator const& collision, TrackEle const& tracks, o2::aod::EMCALMatchedTracks const& matchedtracks)
  {
    registry.fill(HIST("hEventCounter"), 0.5);

    // LOGF(info, "Collision index : %d", collision.index());
    // LOGF(info, "Number of tracks: %d", tracks.size());
    // LOGF(info, "Number of clusters: %d", clusters.size());

    registry.fill(HIST("hZvtx"), collision.posZ());

    for (const auto& track : tracks) {

      registry.fill(HIST("hEta"), track.eta());
      registry.fill(HIST("hITSchi2"), track.itsChi2NCl());
      registry.fill(HIST("hTPCchi2"), track.tpcChi2NCl());
      registry.fill(HIST("hTPCnCls"), track.tpcNClsFound());
      registry.fill(HIST("hITSnCls"), track.itsNCls());
      registry.fill(HIST("hTPCnClsCrossedRows"), track.tpcNClsCrossedRows());

      if (std::abs(track.eta()) > etaupCut)
        continue;
      if (track.tpcNClsCrossedRows() < tpcNclCrCut)
        continue;
      if (std::abs(track.dcaXY()) > dcaxyCut)
        continue;
      if (track.itsChi2NCl() > itschi2Cut)
        continue;
      if (track.tpcChi2NCl() > tpcchi2Cut)
        continue;
      if (track.tpcNClsFound() < tpcNclCut)
        continue;
      if (track.itsNCls() < itsNclCut)
        continue;
      if (track.pt() < mimpTCut)
        continue;
      registry.fill(HIST("hPt"), track.pt());
      registry.fill(HIST("hTPCNsigma"), track.p(), track.tpcNSigmaEl());

      // track - match

      //  continue;
      if (track.phi() < emcaccPhimin || track.phi() > emcaccPhimax)
        continue;
      auto tracksofcluster = matchedtracks.sliceBy(perClusterMatchedTracks, track.globalIndex());

      // LOGF(info, "Number of matched track: %d", tracksofcluster.size());

      double Rmim = 999.9;
      double dPhi_mim = 999.9;
      double dEta_mim = 999.9;

      if (tracksofcluster.size()) {
        int nmatch = 0;
        for (const auto& match : tracksofcluster) {
          if (match.emcalcluster_as<SelectedClusters>().time() < minTime || match.emcalcluster_as<SelectedClusters>().time() > maxTime)
            continue;
          if (match.emcalcluster_as<SelectedClusters>().m02() < minM02 || match.emcalcluster_as<SelectedClusters>().m02() > maxM02)
            continue;

          float emc_m20 = match.emcalcluster_as<SelectedClusters>().m20();
          float emc_m02 = match.emcalcluster_as<SelectedClusters>().m02();
          float emc_energy = match.emcalcluster_as<SelectedClusters>().energy();
          double emc_phi = match.emcalcluster_as<SelectedClusters>().phi();
          double emc_eta = match.emcalcluster_as<SelectedClusters>().eta();
          double emc_time = match.emcalcluster_as<SelectedClusters>().time();
          // LOG(info) << "tr phi0 = " << match.track_as<TrackEle>().phi();
          // LOG(info) << "tr phi1 = " << track.phi();
          // LOG(info) << "emc phi = " << emc_phi;
          if (nmatch == 0) {
            double dEta = match.track_as<TrackEle>().eta() - emc_eta;
            double dPhi = match.track_as<TrackEle>().phi() - emc_phi;
            if (dPhi > o2::constants::math::PI) {
              dPhi -= 2 * o2::constants::math::PI;
            } else if (dPhi < -o2::constants::math::PI) {
              dPhi += 2 * o2::constants::math::PI;
            }

            registry.fill(HIST("hMatchPhi"), emc_phi, match.track_as<TrackEle>().phi());
            registry.fill(HIST("hMatchEta"), emc_eta, match.track_as<TrackEle>().eta());

            double R = std::sqrt(std::pow(dPhi, 2) + std::pow(dEta, 2));
            if (R < Rmim) {
              Rmim = R;
              dPhi_mim = dPhi;
              dEta_mim = dEta;
            }
            registry.fill(HIST("hTrMatch"), dPhi, dEta);
            registry.fill(HIST("hEMCtime"), emc_time);
            registry.fill(HIST("hEnergy"), emc_energy);

            if (R < MatchRCut)
              continue;

            double eop = emc_energy / match.track_as<TrackEle>().p();
            // LOG(info) << "E/p" << eop;
            registry.fill(HIST("hEopNsigTPC"), match.track_as<TrackEle>().tpcNSigmaEl(), eop);
            registry.fill(HIST("hM02"), match.track_as<TrackEle>().tpcNSigmaEl(), emc_m02);
            registry.fill(HIST("hM20"), match.track_as<TrackEle>().tpcNSigmaEl(), emc_m20);
            if (match.track_as<TrackEle>().tpcNSigmaEl() > tpcNsiglowCut && match.track_as<TrackEle>().tpcNSigmaEl() < tpcNsigupCut) {
              registry.fill(HIST("hEop"), match.track_as<TrackEle>().pt(), eop);
            }
          }

          nmatch++;
        }
      }

      if (Rmim < MatchRCut) {
        // LOG(info) << "R mim = " << Rmim;
        registry.fill(HIST("hTrMatch_mim"), dPhi_mim, dEta_mim);
      }

    } // end of track loop
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfTaskElectronWeakBoson>(cfgc)};
}
