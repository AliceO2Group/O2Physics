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

/// \file skimmerGammaCalo.cxx
/// \brief skim cluster information to write photon cluster table into derived AO2D.root
/// \author marvin.hemmer@cern.ch
/// dependencies: emcal-correction-task

#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/emcalHistoDefinitions.h"
#include "PWGJE/DataModel/EMCALClusters.h"

#include "Common/CCDB/TriggerAliases.h"
#include "Common/DataModel/EventSelection.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct SkimmerGammaCalo {

  Preslice<o2::aod::EMCALClusterCells> psCellperCluster = o2::aod::emcalclustercell::emcalclusterId;
  Preslice<o2::aod::EMCALMatchedTracks> psMTperCluster = o2::aod::emcalclustercell::emcalclusterId;
  Preslice<o2::aod::EMCMatchSecs> psMSperCluster = o2::aod::emcalclustercell::emcalclusterId;

  Produces<aod::SkimEMCClusters_001> tableGammaEMCReco;
  Produces<aod::EMCClusterMCLabels> tableEMCClusterMCLabels;
  Produces<aod::SkimEMCCells> tableCellEMCReco;

  // Configurable for filter/cuts
  Configurable<float> minTime{"minTime", -200., "Minimum cluster time for time cut"};
  Configurable<float> maxTime{"maxTime", +200., "Maximum cluster time for time cut"};
  Configurable<float> minM02{"minM02", 0.0, "Minimum M02 for M02 cut"};
  Configurable<float> maxM02{"maxM02", 1.0, "Maximum M02 for M02 cut"};
  Configurable<float> minE{"minE", 0.5, "Minimum energy for energy cut"};
  Configurable<bool> removeExotic{"removeExotic", false, "Flag to enable the removal of exotic clusters."};
  Configurable<std::vector<int>> clusterDefinitions{"clusterDefinitions", {0, 1, 2, 10, 11, 12, 13, 20, 21, 22, 30, 40, 41, 42, 43, 44, 45}, "Cluster definitions to be accepted (e.g. 13 for kV3MostSplitLowSeed)"};
  Configurable<float> maxdEta{"maxdEta", 0.1, "Set a maximum difference in eta for tracks and cluster to still count as matched"};
  Configurable<float> maxdPhi{"maxdPhi", 0.1, "Set a maximum difference in phi for tracks and cluster to still count as matched"};
  Configurable<float> maxEoverP{"maxEoverP", 1.5, "Set a maximum for cluster E / track p for track matching."};
  Configurable<float> maxdEtaSec{"maxdEtaSec", 0.1, "Set a maximum difference in eta for secondary tracks and cluster to still count as matched"};
  Configurable<float> maxdPhiSec{"maxdPhiSec", 0.1, "Set a maximum difference in phi for secondary tracks and cluster to still count as matched"};
  Configurable<bool> needEMCTrigger{"needEMCTrigger", false, "flag to only save events which have kTVXinEMC trigger bit. To reduce PbPb derived data size"};

  HistogramRegistry historeg{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  void init(o2::framework::InitContext&)
  {
    historeg.add("DefinitionIn", "Cluster definitions before cuts;#bf{Cluster definition};#bf{#it{N}_{clusters}}", HistType::kTH1F, {{51, -0.5, 50.5}});
    historeg.add("DefinitionOut", "Cluster definitions after cuts;#bf{Cluster definition};#bf{#it{N}_{clusters}}", HistType::kTH1F, {{51, -0.5, 50.5}});
    historeg.add("EIn", "Energy of clusters before cuts", gHistoSpecClusterE);
    historeg.add("EOut", "Energy of clusters after cuts", gHistoSpecClusterE);
    historeg.add("MTEtaPhiBeforeTM", "Eta phi of matched tracks before TM cuts", gHistoSpecClusterTMdEtadPhi);
    historeg.add("MTEtaPhiAfterTM", "Eta phi of matched tracks after TM cuts", gHistoSpecClusterTMdEtadPhi);
    historeg.add("MSTEtaPhiBeforeTM", "Eta phi of matched secondary tracks before TM cuts", gHistoSpecClusterTMdEtadPhi);
    historeg.add("MSTEtaPhiAfterTM", "Eta phi of matched secondary tracks after TM cuts", gHistoSpecClusterTMdEtadPhi);
    historeg.add("Eoverp", "E/p for cluster E and track p", gHistoSpecTMEoverP);
    historeg.add("M02In", "Shape of cluster before cuts;#bf{#it{M}_{02}};#bf{#it{N}_{clusters}}", HistType::kTH1F, {{200, 0, 2}});
    historeg.add("M02Out", "Shape of cluster after cuts;#bf{#it{M}_{02}};#bf{#it{N}_{clusters}}", HistType::kTH1F, {{200, 0, 2}});
    historeg.add("TimeIn", "Time of cluster before cuts;#bf{#it{t} (ns)};#bf{#it{N}_{clusters}}", HistType::kTH1F, {{200, -100, 100}});
    historeg.add("TimeOut", "Time of cluster after cuts;#bf{#it{t} (ns)};#bf{#it{N}_{clusters}}", HistType::kTH1F, {{200, -100, 100}});

    auto hCaloClusterFilter = historeg.add<TH1>("hCaloClusterFilter", "hCaloClusterFilter", kTH1I, {{7, 0, 7}});
    hCaloClusterFilter->GetXaxis()->SetBinLabel(1, "in");
    hCaloClusterFilter->GetXaxis()->SetBinLabel(2, "Definition cut");
    hCaloClusterFilter->GetXaxis()->SetBinLabel(3, "E cut");
    hCaloClusterFilter->GetXaxis()->SetBinLabel(4, "time cut");
    hCaloClusterFilter->GetXaxis()->SetBinLabel(5, "M02 cut");
    hCaloClusterFilter->GetXaxis()->SetBinLabel(6, "exotic cut");
    hCaloClusterFilter->GetXaxis()->SetBinLabel(7, "out");

    auto hCaloTrackFilter = historeg.add<TH1>("hCaloTrackFilter", "hCaloTrackFilter", kTH1I, {{4, 0, 4}});
    hCaloTrackFilter->GetXaxis()->SetBinLabel(1, "in");
    hCaloTrackFilter->GetXaxis()->SetBinLabel(2, "#Delta#eta #Delta#varphi");
    hCaloTrackFilter->GetXaxis()->SetBinLabel(3, "E/p cut");
    hCaloTrackFilter->GetXaxis()->SetBinLabel(4, "out");

    auto hCaloSecondaryTrackFilter = historeg.add<TH1>("hCaloSecondaryTrackFilter", "hCaloSecondaryTrackFilter", kTH1I, {{4, 0, 4}});
    hCaloSecondaryTrackFilter->GetXaxis()->SetBinLabel(1, "in");
    hCaloSecondaryTrackFilter->GetXaxis()->SetBinLabel(2, "#Delta#eta #Delta#varphi");
    hCaloSecondaryTrackFilter->GetXaxis()->SetBinLabel(3, "E/p cut");
    hCaloSecondaryTrackFilter->GetXaxis()->SetBinLabel(4, "out");

    LOG(info) << "| EMCal cluster cuts for skimming:";
    LOG(info) << "| Timing cut: " << minTime << " < t < " << maxTime;
    LOG(info) << "| M02 cut: " << minM02 << " < M02 < " << maxM02;
    LOG(info) << "| E cut: E > " << minE;
    LOG(info) << "| TM - dPhi cut: dPhi < " << maxdPhi;
    LOG(info) << "| TM - dEta cut: dEta < " << maxdEta;
    LOG(info) << "| TM - E/p cut: E/p < " << maxEoverP;
  }

  template <typename TSecondaries>
  static constexpr bool HasSecondaries = !std::is_same_v<TSecondaries, std::nullptr_t>;

  template <typename TCollision, typename TClusters, typename TClusterCells, typename TTracks, typename TMatchedTracks, typename TMatchedSecondaries = std::nullptr_t>
  void runAnalysis(TCollision const& collision, TClusters const& emcclusters, TClusterCells const& emcclustercells, TMatchedTracks const& emcmatchedtracks, TTracks const& /*tracks*/, TMatchedSecondaries const& secondaries = nullptr)
  {
    if (!collision.isSelected()) {
      return;
    }
    if (needEMCTrigger.value && !collision.alias_bit(kTVXinEMC)) {
      return;
    }

    for (const auto& emccluster : emcclusters) {
      historeg.fill(HIST("hCaloClusterFilter"), 0);

      historeg.fill(HIST("DefinitionIn"), emccluster.definition());
      historeg.fill(HIST("M02In"), emccluster.m02());
      historeg.fill(HIST("TimeIn"), emccluster.time());
      historeg.fill(HIST("EIn"), emccluster.energy());

      // Definition cut
      if (!(std::find(clusterDefinitions.value.begin(), clusterDefinitions.value.end(), emccluster.definition()) != clusterDefinitions.value.end())) {
        historeg.fill(HIST("hCaloClusterFilter"), 1);
        continue;
      }
      // Energy cut
      if (emccluster.energy() < minE) {
        historeg.fill(HIST("hCaloClusterFilter"), 2);
        continue;
      }
      // timing cut
      if (emccluster.time() > maxTime || emccluster.time() < minTime) {
        historeg.fill(HIST("hCaloClusterFilter"), 3);
        continue;
      }
      // M02 cut
      if (emccluster.nCells() > 1 && (emccluster.m02() > maxM02 || emccluster.m02() < minM02)) {
        historeg.fill(HIST("hCaloClusterFilter"), 4);
        continue;
      }
      if (removeExotic.value && emccluster.isExotic()) {
        historeg.fill(HIST("hCaloClusterFilter"), 5);
        continue;
      }
      historeg.fill(HIST("hCaloClusterFilter"), 6);

      // Skimmed cell table
      auto groupedCells = emcclustercells.sliceBy(psCellperCluster, emccluster.globalIndex());
      for (const auto& emcclustercell : groupedCells) {
        tableCellEMCReco(emcclustercell.emcalclusterId(), emcclustercell.caloId());
      }

      // Skimmed matched tracks table
      std::vector<float> vEta;
      std::vector<float> vPhi;
      std::vector<float> vP;
      std::vector<float> vPt;
      auto groupedMTs = emcmatchedtracks.sliceBy(psMTperCluster, emccluster.globalIndex());
      vEta.reserve(groupedMTs.size());
      vPhi.reserve(groupedMTs.size());
      vP.reserve(groupedMTs.size());
      vPt.reserve(groupedMTs.size());
      for (const auto& emcmatchedtrack : groupedMTs) {
        historeg.fill(HIST("hCaloTrackFilter"), 0);
        historeg.fill(HIST("MTEtaPhiBeforeTM"), emcmatchedtrack.deltaEta(), emcmatchedtrack.deltaPhi());
        if (std::abs(emcmatchedtrack.deltaEta()) >= maxdEta || std::abs(emcmatchedtrack.deltaPhi()) >= maxdPhi) {
          historeg.fill(HIST("hCaloTrackFilter"), 1);
          continue;
        }
        historeg.fill(HIST("Eoverp"), emccluster.energy(), emccluster.energy() / emcmatchedtrack.template track_as<aod::FullTracks>().p());
        if (emccluster.energy() / emcmatchedtrack.template track_as<aod::FullTracks>().p() > maxEoverP) {
          historeg.fill(HIST("hCaloTrackFilter"), 2);
          continue;
        }
        historeg.fill(HIST("hCaloTrackFilter"), 3);
        historeg.fill(HIST("MTEtaPhiAfterTM"), emcmatchedtrack.deltaEta(), emcmatchedtrack.deltaPhi());
        vEta.emplace_back(emcmatchedtrack.deltaEta());
        vPhi.emplace_back(emcmatchedtrack.deltaPhi());
        vP.emplace_back(emcmatchedtrack.template track_as<aod::FullTracks>().p());
        vPt.emplace_back(emcmatchedtrack.template track_as<aod::FullTracks>().pt());
      }

      std::vector<float> vEtaSecondaries = {};
      std::vector<float> vPhiSecondaries = {};
      std::vector<float> vPSecondaries = {};
      std::vector<float> vPtSecondaries = {};

      if constexpr (HasSecondaries<TMatchedSecondaries>) {
        auto groupedMatchedSecondaries = secondaries.sliceBy(psMSperCluster, emccluster.globalIndex());
        vEta.reserve(groupedMatchedSecondaries.size());
        vPhi.reserve(groupedMatchedSecondaries.size());
        vP.reserve(groupedMatchedSecondaries.size());
        vPt.reserve(groupedMatchedSecondaries.size());
        for (const auto& emcMatchedSecondary : groupedMatchedSecondaries) {
          historeg.fill(HIST("hCaloSecondaryTrackFilter"), 0);
          historeg.fill(HIST("MSTEtaPhiBeforeTM"), emcMatchedSecondary.deltaEta(), emcMatchedSecondary.deltaPhi());
          if (std::abs(emcMatchedSecondary.deltaEta()) >= maxdEtaSec || std::abs(emcMatchedSecondary.deltaPhi()) >= maxdPhiSec) {
            historeg.fill(HIST("hCaloSecondaryTrackFilter"), 1);
            continue;
          }
          historeg.fill(HIST("hCaloSecondaryTrackFilter"), 3);
          historeg.fill(HIST("MSTEtaPhiAfterTM"), emcMatchedSecondary.deltaEta(), emcMatchedSecondary.deltaPhi());
          vEta.emplace_back(emcMatchedSecondary.deltaEta());
          vPhi.emplace_back(emcMatchedSecondary.deltaPhi());
          vP.emplace_back(emcMatchedSecondary.template track_as<aod::FullTracks>().p());
          vPt.emplace_back(emcMatchedSecondary.template track_as<aod::FullTracks>().pt());
        }
      }

      historeg.fill(HIST("DefinitionOut"), emccluster.definition());
      historeg.fill(HIST("EOut"), emccluster.energy());
      historeg.fill(HIST("M02Out"), emccluster.m02());
      historeg.fill(HIST("TimeOut"), emccluster.time());

      tableGammaEMCReco(emccluster.collisionId(), emccluster.definition(), emccluster.energy(), emccluster.eta(), emccluster.phi(), emccluster.m02(),
                        emccluster.nCells(), emccluster.time(), emccluster.isExotic(), vPhi, vEta, vP, vPt, vPhiSecondaries, vEtaSecondaries, vPSecondaries, vPtSecondaries);
    }
  }

  void processRec(soa::Join<aod::Collisions, aod::EvSels, aod::EMEvSels>::iterator const& collision, aod::EMCALClusters const& emcclusters, aod::EMCALClusterCells const& emcclustercells, aod::EMCALMatchedTracks const& emcmatchedtracks, aod::FullTracks const& tracks)
  {
    runAnalysis(collision, emcclusters, emcclustercells, emcmatchedtracks, tracks);
  }
  PROCESS_SWITCH(SkimmerGammaCalo, processRec, "process only reconstructed info", true);

  void processRecWithSecondaries(soa::Join<aod::Collisions, aod::EvSels, aod::EMEvSels>::iterator const& collision, aod::EMCALClusters const& emcclusters, aod::EMCALClusterCells const& emcclustercells, aod::EMCALMatchedTracks const& emcmatchedtracks, aod::FullTracks const& tracks, aod::EMCMatchSecs const& emcmatchedsecondaries)
  {
    runAnalysis(collision, emcclusters, emcclustercells, emcmatchedtracks, tracks, emcmatchedsecondaries);
  }
  PROCESS_SWITCH(SkimmerGammaCalo, processRecWithSecondaries, "process reconstructed info with secondary track matching.", false);

  void processMC(soa::Join<aod::Collisions, aod::EvSels, aod::EMEvSels>::iterator const& collision, soa::Join<aod::EMCALClusters, aod::EMCALMCClusters> const& emcclusters, aod::McParticles const&)
  {
    if (!collision.isSelected()) {
      return;
    }

    if (needEMCTrigger.value && !collision.alias_bit(kTVXinEMC)) {
      return;
    }

    for (const auto& emccluster : emcclusters) {

      // Definition cut
      if (!(std::find(clusterDefinitions.value.begin(), clusterDefinitions.value.end(), emccluster.definition()) != clusterDefinitions.value.end())) {
        continue;
      }
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
  PROCESS_SWITCH(SkimmerGammaCalo, processMC, "process MC info", false); // Run this in addition to processRec for MCs to copy the cluster mc labels from the EMCALMCClusters to the skimmed EMCClusterMCLabels table

  void processDummy(aod::Collision const&)
  {
    // do nothing
  }
  PROCESS_SWITCH(SkimmerGammaCalo, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<SkimmerGammaCalo>(cfgc)};
  return workflow;
}
