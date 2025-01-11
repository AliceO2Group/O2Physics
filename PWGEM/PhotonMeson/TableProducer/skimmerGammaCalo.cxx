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

/// \brief skim cluster information to write photon cluster table in AO2D.root
/// dependencies: emcal-correction-task
/// \author marvin.hemmer@cern.ch

#include <algorithm>
#include <vector>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

#include "Common/Core/TableHelper.h"

// includes for the R recalculation
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "CCDB/BasicCCDBManager.h"

#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/emcalHistoDefinitions.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct skimmerGammaCalo {

  Preslice<o2::aod::EMCALClusterCells> CellperCluster = o2::aod::emcalclustercell::emcalclusterId;
  Preslice<o2::aod::EMCALMatchedTracks> MTperCluster = o2::aod::emcalclustercell::emcalclusterId;

  Produces<aod::SkimEMCClusters> tableGammaEMCReco;
  Produces<aod::EMCClusterMCLabels> tableEMCClusterMCLabels;
  Produces<aod::SkimEMCCells> tableCellEMCReco;
  Produces<aod::SkimEMCMTs> tableTrackEMCReco;

  // Configurable for filter/cuts
  Configurable<float> minTime{"minTime", -200., "Minimum cluster time for time cut"};
  Configurable<float> maxTime{"maxTime", +200., "Maximum cluster time for time cut"};
  Configurable<float> minM02{"minM02", 0.0, "Minimum M02 for M02 cut"};
  Configurable<float> maxM02{"maxM02", 1.0, "Maximum M02 for M02 cut"};
  Configurable<float> minE{"minE", 0.5, "Minimum energy for energy cut"};
  Configurable<float> maxdEta{"maxdEta", 0.1, "Set a maximum difference in eta for tracks and cluster to still count as matched"};
  Configurable<float> maxdPhi{"maxdPhi", 0.1, "Set a maximum difference in phi for tracks and cluster to still count as matched"};
  Configurable<bool> applyEveSel_at_skimming{"applyEveSel_at_skimming", false, "flag to apply minimal event selection at the skimming level"};
  Configurable<bool> inherit_from_emevent_photon{"inherit_from_emevent_photon", false, "flag to inherit task options from emevent-photon"};

  HistogramRegistry historeg{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  void init(o2::framework::InitContext& initContext)
  {
    historeg.add("hCaloClusterEIn", "hCaloClusterEIn", gHistoSpec_clusterE);
    historeg.add("hCaloClusterEOut", "hCaloClusterEOut", gHistoSpec_clusterE);
    historeg.add("hMTEtaPhi", "hMTEtaPhi", gHistoSpec_clusterTM_dEtadPhi);
    auto hCaloClusterFilter = historeg.add<TH1>("hCaloClusterFilter", "hCaloClusterFilter", kTH1I, {{5, 0, 5}});
    hCaloClusterFilter->GetXaxis()->SetBinLabel(1, "in");
    hCaloClusterFilter->GetXaxis()->SetBinLabel(2, "E cut");
    hCaloClusterFilter->GetXaxis()->SetBinLabel(3, "time cut");
    hCaloClusterFilter->GetXaxis()->SetBinLabel(4, "M02 cut");
    hCaloClusterFilter->GetXaxis()->SetBinLabel(5, "out");

    if (inherit_from_emevent_photon) {
      getTaskOptionValue(initContext, "create-emevent-photon", "applyEveSel_at_skimming", applyEveSel_at_skimming.value, true); // for EM users.
    }

    LOG(info) << "| Timing cut: " << minTime << " < t < " << maxTime;
    LOG(info) << "| M02 cut: " << minM02 << " < M02 < " << maxM02;
    LOG(info) << "| E cut: E > " << minE;
  }

  void processRec(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::EMCALClusters const& emcclusters, aod::EMCALClusterCells const& emcclustercells, aod::EMCALMatchedTracks const& emcmatchedtracks, aod::FullTracks const&)
  {
    if (applyEveSel_at_skimming && (!collision.selection_bit(o2::aod::evsel::kIsTriggerTVX) || !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder))) {
      return;
    }
    for (const auto& emccluster : emcclusters) {
      historeg.fill(HIST("hCaloClusterEIn"), emccluster.energy());
      historeg.fill(HIST("hCaloClusterFilter"), 0);

      // Energy cut
      if (emccluster.energy() < minE) {
        historeg.fill(HIST("hCaloClusterFilter"), 1);
        continue;
      }
      // timing cut
      if (emccluster.time() > maxTime || emccluster.time() < minTime) {
        historeg.fill(HIST("hCaloClusterFilter"), 2);
        continue;
      }
      // M02 cut
      if (emccluster.nCells() > 1 && (emccluster.m02() > maxM02 || emccluster.m02() < minM02)) {
        historeg.fill(HIST("hCaloClusterFilter"), 3);
        continue;
      }

      // Skimmed cell table
      auto groupedCells = emcclustercells.sliceBy(CellperCluster, emccluster.globalIndex());
      for (const auto& emcclustercell : groupedCells) {
        tableCellEMCReco(emcclustercell.emcalclusterId(), emcclustercell.caloId());
      }

      // Skimmed matched tracks table
      std::vector<int32_t> vTrackIds;
      std::vector<float> vEta;
      std::vector<float> vPhi;
      std::vector<float> vP;
      std::vector<float> vPt;
      auto groupedMTs = emcmatchedtracks.sliceBy(MTperCluster, emccluster.globalIndex());
      vTrackIds.reserve(groupedMTs.size());
      vEta.reserve(groupedMTs.size());
      vPhi.reserve(groupedMTs.size());
      vP.reserve(groupedMTs.size());
      vPt.reserve(groupedMTs.size());
      for (const auto& emcmatchedtrack : groupedMTs) {
        if (std::abs(emccluster.eta() - emcmatchedtrack.track_as<aod::FullTracks>().trackEtaEmcal()) >= maxdEta || std::abs(emccluster.phi() - emcmatchedtrack.track_as<aod::FullTracks>().trackPhiEmcal()) >= maxdPhi) {
          continue;
        }
        historeg.fill(HIST("hMTEtaPhi"), emccluster.eta() - emcmatchedtrack.track_as<aod::FullTracks>().trackEtaEmcal(), emccluster.phi() - emcmatchedtrack.track_as<aod::FullTracks>().trackPhiEmcal());
        vTrackIds.emplace_back(emcmatchedtrack.trackId());
        vEta.emplace_back(emcmatchedtrack.track_as<aod::FullTracks>().trackEtaEmcal());
        vPhi.emplace_back(emcmatchedtrack.track_as<aod::FullTracks>().trackPhiEmcal());
        vP.emplace_back(emcmatchedtrack.track_as<aod::FullTracks>().p());
        vPt.emplace_back(emcmatchedtrack.track_as<aod::FullTracks>().pt());
        tableTrackEMCReco(emcmatchedtrack.emcalclusterId(), emcmatchedtrack.track_as<aod::FullTracks>().trackEtaEmcal(), emcmatchedtrack.track_as<aod::FullTracks>().trackPhiEmcal(),
                          emcmatchedtrack.track_as<aod::FullTracks>().p(), emcmatchedtrack.track_as<aod::FullTracks>().pt());
      }

      historeg.fill(HIST("hCaloClusterEOut"), emccluster.energy());
      historeg.fill(HIST("hCaloClusterFilter"), 4);

      tableGammaEMCReco(emccluster.collisionId(), emccluster.definition(), emccluster.energy(), emccluster.eta(), emccluster.phi(), emccluster.m02(),
                        emccluster.nCells(), emccluster.time(), emccluster.isExotic(), vEta, vPhi, vP, vPt);
    }
  }
  void processMC(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Join<aod::EMCALClusters, aod::EMCALMCClusters> const& emcclusters, aod::McParticles const&)
  {
    if (applyEveSel_at_skimming && (!collision.selection_bit(o2::aod::evsel::kIsTriggerTVX) || !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder))) {
      return;
    }
    for (const auto& emccluster : emcclusters) {

      // Energy cut
      if (emccluster.energy() < minE) {
        continue;
      }
      // timing cut
      if (emccluster.time() > maxTime || emccluster.time() < minTime) {
        continue;
      }
      // M02 cut
      if (emccluster.nCells() > 1 && (emccluster.m02() > maxM02 || emccluster.m02() < minM02)) {
        continue;
      }
      std::vector<int32_t> mcLabels;
      for (size_t iCont = 0; iCont < emccluster.amplitudeA().size(); iCont++) {
        mcLabels.push_back(emccluster.mcParticleIds()[iCont]);
      }
      // LOGF(info, "---- New Cluster ---");
      // for (unsigned long int iCont = 0; iCont < mcLabels.size(); iCont++) {
      //   LOGF(info, "iCont = %d, mcParticle = %d, amplitudeA = %.5f", iCont, mcLabels.at(iCont), emccluster.amplitudeA()[iCont]);
      // }
      tableEMCClusterMCLabels(mcLabels);
      mcLabels.clear();
    }
  }
  PROCESS_SWITCH(skimmerGammaCalo, processRec, "process only reconstructed info", true);
  PROCESS_SWITCH(skimmerGammaCalo, processMC, "process MC info", false); // Run this in addition to processRec for MCs to copy the cluster mc labels from the EMCALMCClusters to the skimmed EMCClusterMCLabels table

  void processDummy(aod::Collision const&)
  {
    // do nothing
  }
  PROCESS_SWITCH(skimmerGammaCalo, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<skimmerGammaCalo>(cfgc, TaskName{"skimmer-gamma-calo"})};
  return workflow;
}
