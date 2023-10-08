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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

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
  Produces<aod::SkimEMCCells> tableCellEMCReco;
  Produces<aod::SkimEMCMTs> tableTrackEMCReco;

  // Configurable for filter/cuts
  Configurable<float> minTime{"minTime", -200., "Minimum cluster time for time cut"};
  Configurable<float> maxTime{"maxTime", +200., "Maximum cluster time for time cut"};
  Configurable<float> minM02{"minM02", 0.0, "Minimum M02 for M02 cut"};
  Configurable<float> maxM02{"maxM02", 1.0, "Maximum M02 for M02 cut"};
  Configurable<bool> hasPropagatedTracks{"hasPropagatedTracks", false, "temporary flag, only set to true when running over data which has the tracks propagated to EMCal/PHOS!"};

  HistogramRegistry historeg{
    "historeg",
    {},
    OutputObjHandlingPolicy::QAObject,
    false,
    true};

  void init(o2::framework::InitContext&)
  {
    historeg.add("hCaloClusterEIn", "hCaloClusterEIn", gHistoSpec_clusterE);
    historeg.add("hCaloClusterEOut", "hCaloClusterEOut", gHistoSpec_clusterE);
    historeg.add("hMTEtaPhi", "hMTEtaPhi", gHistoSpec_EtaPhi);
    auto hCaloClusterFilter = historeg.add<TH1>("hCaloClusterFilter", "hCaloClusterFilter", kTH1I, {{4, 0, 4}});
    hCaloClusterFilter->GetXaxis()->SetBinLabel(1, "in");
    hCaloClusterFilter->GetXaxis()->SetBinLabel(2, "time cut");
    hCaloClusterFilter->GetXaxis()->SetBinLabel(3, "M02 cut");
    hCaloClusterFilter->GetXaxis()->SetBinLabel(4, "out");

    LOG(info) << "| Timing cut: " << minTime << " < t < " << maxTime << std::endl;
    LOG(info) << "| M02 cut: " << minM02 << " < M02 < " << maxM02 << std::endl;
  }

  void processRec(aod::Collision const&, aod::EMCALClusters const& emcclusters, aod::EMCALClusterCells const& emcclustercells, aod::EMCALMatchedTracks const& emcmatchedtracks, aod::FullTracks const& tracks)
  {
    for (const auto& emccluster : emcclusters) {
      historeg.fill(HIST("hCaloClusterEIn"), emccluster.energy());
      historeg.fill(HIST("hCaloClusterFilter"), 0);

      // timing cut
      if (emccluster.time() > maxTime || emccluster.time() < minTime) {
        historeg.fill(HIST("hCaloClusterFilter"), 1);
        continue;
      }
      // M02 cut
      if (emccluster.m02() > maxM02 || emccluster.m02() < minM02) {
        historeg.fill(HIST("hCaloClusterFilter"), 2);
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
        if (hasPropagatedTracks) { // only temporarily while not every data has the tracks propagated to EMCal/PHOS
          historeg.fill(HIST("hMTEtaPhi"), emccluster.eta() - emcmatchedtrack.track_as<aod::FullTracks>().trackEtaEmcal(), emccluster.phi() - emcmatchedtrack.track_as<aod::FullTracks>().trackPhiEmcal());
          vTrackIds.emplace_back(emcmatchedtrack.trackId());
          vEta.emplace_back(emcmatchedtrack.track_as<aod::FullTracks>().trackEtaEmcal());
          vPhi.emplace_back(emcmatchedtrack.track_as<aod::FullTracks>().trackPhiEmcal());
          vP.emplace_back(emcmatchedtrack.track_as<aod::FullTracks>().p());
          vPt.emplace_back(emcmatchedtrack.track_as<aod::FullTracks>().pt());
          tableTrackEMCReco(emcmatchedtrack.emcalclusterId(), emcmatchedtrack.track_as<aod::FullTracks>().trackEtaEmcal(), emcmatchedtrack.track_as<aod::FullTracks>().trackPhiEmcal(),
                            emcmatchedtrack.track_as<aod::FullTracks>().p(), emcmatchedtrack.track_as<aod::FullTracks>().pt());
        } else {
          historeg.fill(HIST("hMTEtaPhi"), emccluster.eta() - emcmatchedtrack.track_as<aod::FullTracks>().eta(), emccluster.phi() - emcmatchedtrack.track_as<aod::FullTracks>().phi());
          vTrackIds.emplace_back(emcmatchedtrack.trackId());
          vEta.emplace_back(emcmatchedtrack.track_as<aod::FullTracks>().eta());
          vPhi.emplace_back(emcmatchedtrack.track_as<aod::FullTracks>().phi());
          vP.emplace_back(emcmatchedtrack.track_as<aod::FullTracks>().p());
          vPt.emplace_back(emcmatchedtrack.track_as<aod::FullTracks>().pt());
          tableTrackEMCReco(emcmatchedtrack.emcalclusterId(), emcmatchedtrack.track_as<aod::FullTracks>().eta(), emcmatchedtrack.track_as<aod::FullTracks>().phi(),
                            emcmatchedtrack.track_as<aod::FullTracks>().p(), emcmatchedtrack.track_as<aod::FullTracks>().pt());
        }
      }

      historeg.fill(HIST("hCaloClusterEOut"), emccluster.energy());
      historeg.fill(HIST("hCaloClusterFilter"), 3);

      tableGammaEMCReco(emccluster.collisionId(), emccluster.id(),
                        emccluster.energy(), emccluster.coreEnergy(), emccluster.eta(), emccluster.phi(), emccluster.m02(),
                        emccluster.m20(), emccluster.nCells(), emccluster.time(), emccluster.isExotic(), emccluster.distanceToBadChannel(), emccluster.nlm(),
                        emccluster.definition(), vTrackIds, vEta, vPhi, vP, vPt);
    }
  }
  PROCESS_SWITCH(skimmerGammaCalo, processRec, "process only reconstructed info", true);

  void processDummy(aod::Collision const& collision)
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
