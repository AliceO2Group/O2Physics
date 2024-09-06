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
/// \file taskElectronWeakNoson.cxx
/// \briff task for WeakBoson (W/Z) based on electron in mid-rapidity
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

  using SelectedClusters = o2::soa::Filtered<o2::aod::EMCALClusters>;

  // PbPb
  using TrackEle = o2::soa::Filtered<o2::soa::Join<o2::aod::Tracks, o2::aod::FullTracks, o2::aod::TracksExtra_001, o2::aod::TracksDCA, o2::aod::TrackSelection, o2::aod::pidTPCFullEl>>;

  // pp
  // using TrackEle = o2::soa::Filtered<o2::soa::Join<o2::aod::Tracks, o2::aod::FullTracks, o2::aod::TracksDCA, o2::aod::TrackSelection, o2::aod::pidTPCEl, o2::aod::pidTOFEl>>;

  Preslice<o2::aod::EMCALClusterCells> perCluster = o2::aod::emcalclustercell::emcalclusterId;
  Preslice<o2::aod::EMCALAmbiguousClusterCells> perClusterAmb = o2::aod::emcalclustercell::emcalambiguousclusterId;
  PresliceUnsorted<o2::aod::EMCALMatchedTracks> perClusterMatchedTracks = o2::aod::emcalmatchedtrack::trackId;

  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pt histo"};
  Configurable<int> nBinsE{"nBinsE", 100, "N bins in E histo"};

  // event filter
  Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Configurable<float> vtxZ{"vtxZ", 10.f, ""};
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < vtxZ);

  // track cuts
  Configurable<float> etalow{"etalow", -0.6f, ""};
  Configurable<float> etaup{"etaup", 0.6f, ""};
  Configurable<float> dcaxy_cut{"dcaxy_cut", 2.0f, ""};
  Configurable<float> mimpT_cut{"mimpT_cut", 3.0f, "minimum pT cut"};
  Configurable<float> itschi2_cut{"itschi2_cut", 15.0f, "its chi2 cut"};
  Configurable<float> tpcchi2_cut{"tpcchi2_cut", 4.0f, "tpc chi2 cut"};
  Configurable<float> itsNcl_cut{"itsNcl_cut", 2.0f, "its # of cluster cut"};
  Configurable<float> tpcNcl_cut{"tpcNcl_cut", 100.0f, "tpc # if cluster cut"};
  Configurable<float> tpcNclCr_cut{"tpcNclCr_cut", 100.0f, "tpc # of crossedRows cut"};

  Filter filter_globalTr = requireGlobalTrackInFilter();
  Filter etafilter = (aod::track::eta < etaup) && (aod::track::eta > etalow);
  Filter dcaxyfilter = (nabs(aod::track::dcaXY) < dcaxy_cut);

  // cluster cut
  Configurable<int> mClusterDefinition{"mClusterDefinition", 10, "cluster definition to be selected, e.g. 10=kV3Default"};
  Configurable<float> minTime{"minTime", -25., "Minimum cluster time for time cut"};
  Configurable<float> maxTime{"maxTime", +20., "Maximum cluster time for time cut"};
  Configurable<float> minM02{"minM02", 0.1, "Minimum M02 for M02 cut"};
  Configurable<float> maxM02{"maxM02", 0.9, "Maximum M02 for M02 cut"};
  Configurable<float> minM20{"minM20", 0.1, "Minimum M20 for M20 cut"};
  Configurable<float> maxM20{"maxM20", 0.6, "Maximum M20 for M20 cut"};
  Configurable<float> MatchR_cut{"MatchR_cut", 0.1, "cluster - track matching cut"};

  Filter clusterDefinitionSelection = (o2::aod::emcalcluster::definition == mClusterDefinition) && (o2::aod::emcalcluster::time >= minTime) && (o2::aod::emcalcluster::time <= maxTime) && (o2::aod::emcalcluster::m02 > minM02) && (o2::aod::emcalcluster::m02 < maxM02);

  // Histogram registry: an object to hold your histograms
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {

    // define axes you want to use
    const AxisSpec axisZvtx{400, -20, 20, "Zvtx"};
    const AxisSpec axisCounter{1, 0, 1, "events"};
    const AxisSpec axisEta{200, -1.0, 1.0, "#eta"};
    const AxisSpec axisPt{nBinsPt, 0, 10, "p_{T}"};
    const AxisSpec axisNsigma{100, -5, 5, "N#sigma"};
    const AxisSpec axisE{nBinsE, 0, 10, "Energy"};
    const AxisSpec axisM02{100, 0, 1, "M02"};
    const AxisSpec axisdPhi{200, -1, 1, "dPhi"};
    const AxisSpec axisdEta{200, -1, 1, "dEta"};
    const AxisSpec axisPhi{350, 0, 7, "Phi"};
    const AxisSpec axisEop{200, 0, 2, "Eop"};
    const AxisSpec axisChi2{500, 0.0, 50.0, "#chi^{2}"};
    const AxisSpec axisCluster{100, 0.0, 200.0, "counts"};
    const AxisSpec axisITSNCls{20, 0.0, 20, "counts"};
    const AxisSpec axisEMCtime{200, -100.0, 100, "EMC time"};

    // create histograms
    histos.add("ZvtxHistogram", "ZvtxHistogram", kTH1F, {axisZvtx});
    histos.add("hEventCounter", "hEventCounter", kTH1F, {axisCounter});
    histos.add("ITS_Chi2_Hist", "ITS #chi^{2} Hist", kTH1F, {axisChi2});
    histos.add("TPC_Chi2_Hist", "TPC #chi^{2} Hist", kTH1F, {axisChi2});
    histos.add("TPC_NCls_Hist", "TPC_NCls_Hist", kTH1F, {axisCluster});
    histos.add("ITS_NCls_Hist", "ITS_NCls_Hist", kTH1F, {axisITSNCls});
    histos.add("TPC_NClsCrossedRows_Hist", "TPC_NClsCrossedRows_Hist", kTH1F, {axisCluster});
    histos.add("etaHistogram", "etaHistogram", kTH1F, {axisEta});
    histos.add("ptHistogram", "ptHistogram", kTH1F, {axisPt});
    histos.add("TPCElHistogram", "TPCElHistogram", kTH2F, {{axisPt}, {axisNsigma}});
    histos.add("EnergyHistogram", "EnergyHistogram", kTH1F, {axisE});
    histos.add("M02Histogram", "M02Histogram", kTH2F, {{axisNsigma}, {axisM02}});
    histos.add("M20Histogram", "M20Histogram", kTH2F, {{axisNsigma}, {axisM02}});
    histos.add("TrMatchHistogram", "TrMatchHistogram", kTH2F, {{axisdPhi}, {axisdEta}});
    histos.add("TrMatchHistogram_mim", "TrMatchHistogram_mim", kTH2F, {{axisdPhi}, {axisdEta}});
    histos.add("MatchPhiHistogram", "MatchPhiHistogram", kTH2F, {{axisPhi}, {axisPhi}});
    histos.add("MatchEtaHistogram", "MatchEtaHistogram", kTH2F, {{axisEta}, {axisEta}});
    histos.add("EopHistogram", "EopHistogram", kTH2F, {{axisPt}, {axisEop}});
    histos.add("EopNsigTPCHistogram", "EopNsigTPCHistogram", kTH2F, {{axisNsigma}, {axisEop}});
    histos.add("EMCtimeHistogram", "EMCtimeHistogram", kTH1F, {axisEMCtime});
  }

  // void process(soa::Filtered<aod::Collisions>::iterator const& collision, SelectedClusters const& clusters, TrackEle const& tracks, o2::aod::EMCALMatchedTracks const& matchedtracks)
  void process(soa::Filtered<aod::Collisions>::iterator const& collision, SelectedClusters const& emcClusters, TrackEle const& tracks, o2::aod::EMCALMatchedTracks const& matchedtracks)
  {
    histos.fill(HIST("hEventCounter"), 0.5);

    // LOGF(info, "Collision index : %d", collision.index());
    // LOGF(info, "Number of tracks: %d", tracks.size());
    // LOGF(info, "Number of clusters: %d", clusters.size());
    // LOGF(info, "Number of clusters: %d", emcClusters.size());

    histos.fill(HIST("ZvtxHistogram"), collision.posZ());

    for (const auto& track : tracks) {

      histos.fill(HIST("etaHistogram"), track.eta());
      histos.fill(HIST("ITS_Chi2_Hist"), track.itsChi2NCl());
      histos.fill(HIST("TPC_Chi2_Hist"), track.tpcChi2NCl());
      histos.fill(HIST("TPC_NCls_Hist"), track.tpcNClsFound());
      histos.fill(HIST("ITS_NCls_Hist"), track.itsNCls());
      histos.fill(HIST("TPC_NClsCrossedRows_Hist"), track.tpcNClsCrossedRows());

      if (std::abs(track.eta()) > etaup)
        continue;
      if (track.tpcNClsCrossedRows() < tpcNclCr_cut)
        continue;
      if (std::abs(track.dcaXY()) > dcaxy_cut)
        continue;
      if (track.itsChi2NCl() > itschi2_cut)
        continue;
      if (track.tpcChi2NCl() > tpcchi2_cut)
        continue;
      if (track.tpcNClsFound() < tpcNcl_cut)
        continue;
      if (track.itsNCls() < itsNcl_cut)
        continue;
      if (track.pt() < mimpT_cut)
        continue;
      histos.fill(HIST("ptHistogram"), track.pt());
      histos.fill(HIST("TPCElHistogram"), track.p(), track.tpcNSigmaEl());

      // track - match

      if (emcClusters.size() < 1)
        continue;
      if (track.phi() < 1.39 || track.phi() > 3.15)
        continue;
      auto tracksofcluster = matchedtracks.sliceBy(perClusterMatchedTracks, track.globalIndex());

      // LOGF(info, "Number of matched track: %d", tracksofcluster.size());

      double Rmim = 999.9;
      double dPhi_mim = 999.9;
      double dEta_mim = 999.9;

      if (tracksofcluster.size() > 0) {
        int nmatch = 0;
        for (const auto& match : tracksofcluster) {
          if (match.emcalcluster_as<SelectedClusters>().time() < minTime || match.emcalcluster_as<SelectedClusters>().time() > maxTime)
            continue;
          if (match.emcalcluster_as<SelectedClusters>().m02() < minM02 || match.emcalcluster_as<SelectedClusters>().m02() > maxM02)
            continue;

          double emc_m20 = match.emcalcluster_as<SelectedClusters>().m20();
          double emc_m02 = match.emcalcluster_as<SelectedClusters>().m02();
          double emc_energy = match.emcalcluster_as<SelectedClusters>().energy();
          double emc_phi = match.emcalcluster_as<SelectedClusters>().phi();
          double emc_eta = match.emcalcluster_as<SelectedClusters>().eta();
          double emc_time = match.emcalcluster_as<SelectedClusters>().time();
          // LOG(info) << "tr phi0 = " << match.track_as<TrackEle>().phi();
          // LOG(info) << "tr phi1 = " << track.phi();
          // LOG(info) << "emc phi = " << emc_phi;
          if (nmatch == 0) {
            double dPhi = match.track_as<TrackEle>().phi() - emc_phi;
            double dEta = match.track_as<TrackEle>().eta() - emc_eta;

            histos.fill(HIST("MatchPhiHistogram"), emc_phi, match.track_as<TrackEle>().phi());
            histos.fill(HIST("MatchEtaHistogram"), emc_eta, match.track_as<TrackEle>().eta());

            double R = sqrt(pow(dPhi, 2) + pow(dEta, 2));
            if (R < Rmim) {
              Rmim = R;
              dPhi_mim = dPhi;
              dEta_mim = dEta;
            }
            histos.fill(HIST("TrMatchHistogram"), dPhi, dEta);
            histos.fill(HIST("EMCtimeHistogram"), emc_time);

            if (R < MatchR_cut)
              continue;

            double eop = emc_energy / match.track_as<TrackEle>().p();
            // LOG(info) << "E/p" << eop;
            histos.fill(HIST("EopNsigTPCHistogram"), match.track_as<TrackEle>().tpcNSigmaEl(), eop);
            histos.fill(HIST("M02Histogram"), match.track_as<TrackEle>().tpcNSigmaEl(), emc_m02);
            histos.fill(HIST("M20Histogram"), match.track_as<TrackEle>().tpcNSigmaEl(), emc_m20);
            if (match.track_as<TrackEle>().tpcNSigmaEl() > -1.0 && match.track_as<TrackEle>().tpcNSigmaEl() < 3) {
              histos.fill(HIST("EopHistogram"), match.track_as<TrackEle>().pt(), eop);
            }
          }

          nmatch++;
        }
      }

      if (Rmim < 10.0) {
        // LOG(info) << "R mim = " << Rmim;
        histos.fill(HIST("TrMatchHistogram_mim"), dPhi_mim, dEta_mim);
      }

    } // end of track loop
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfTaskElectronWeakBoson>(cfgc)};
}
