// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file emcalQC.cxx
/// \brief Task that runs basic EMCal cluster QA for derived data in the EM format
/// \author Nicolas Strangmann (nicolas.strangmann@cern.ch) Goethe University Frankfurt
///

#include "PWGEM/PhotonMeson/Core/EMCPhotonCut.h"
#include "PWGEM/PhotonMeson/Core/EMPhotonEventCut.h"
#include "PWGEM/PhotonMeson/DataModel/GammaTablesRedux.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/ClusterHistograms.h"
#include "PWGEM/PhotonMeson/Utils/EventHistograms.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>

#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::photon;

using MyCollisions = soa::Join<aod::EMEvents_004, aod::EMEventsAlias, aod::EMEventsMult_000, aod::EMEventsCent_000, aod::EMEventsQvec_001, aod::EMEventsWeight>;
using MyCollision = MyCollisions::iterator;

using EMCalPhotons = soa::Join<aod::EMCEMEventIds, aod::MinClusters>;
using EMCalPhoton = EMCalPhotons::iterator;

struct EmcalQC {

  Configurable<bool> cfgDo2DQA{"cfgDo2DQA", true, "perform 2 dimensional cluster QA"};

  EMPhotonEventCut fEMEventCut;
  struct : ConfigurableGroup {
    std::string prefix = "eventcut_group";
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", 10.f, "max. Zvtx"};
    Configurable<bool> cfgRequireSel8{"cfgRequireSel8", true, "require sel8 in event cut"};
    Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "require FT0AND in event cut"};
    Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", false, "require No time frame border in event cut"};
    Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", false, "require no ITS readout frame border in event cut"};
    Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "require no same bunch pileup in event cut"};
    Configurable<bool> cfgRequireVertexITSTPC{"cfgRequireVertexITSTPC", false, "require Vertex ITSTPC in event cut"}; // ITS-TPC matched track contributes PV.
    Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", false, "require good Zvtx between FT0 vs. PV in event cut"};
    Configurable<bool> cfgRequireEMCReadoutInMB{"cfgRequireEMCReadoutInMB", false, "require the EMC to be read out in an MB collision (kTVXinEMC)"};
    Configurable<bool> cfgRequireEMCHardwareTriggered{"cfgRequireEMCHardwareTriggered", false, "require the EMC to be hardware triggered (kEMC7 or kDMC7)"};
    Configurable<int> cfgOccupancyMin{"cfgOccupancyMin", -1, "min. occupancy"};
    Configurable<int> cfgOccupancyMax{"cfgOccupancyMax", 1000000000, "max. occupancy"};
    Configurable<bool> onlyKeepWeightedEvents{"onlyKeepWeightedEvents", false, "flag to keep only weighted events (for JJ MCs) and remove all MB events (with weight = 1)"};
  } eventcuts;

  EMCPhotonCut fEMCCut;
  struct : ConfigurableGroup {
    std::string prefix = "emccut_group";
    Configurable<std::string> clusterDefinition{"clusterDefinition", "kV3Default", "Clusterizer to be selected, e.g. V3Default"};
    Configurable<float> minOpenAngle{"minOpenAngle", 0.0202, "apply min opening angle"};
    Configurable<float> minClusterTime{"minClusterTime", -20., "Minimum cluster time for EMCal time cut"};
    Configurable<float> maxClusterTime{"maxClusterTime", +25., "Maximum cluster time for EMCal time cut"};
    Configurable<float> minM02{"minM02", 0.1, "Minimum M02 for EMCal M02 cut"};
    Configurable<float> maxM02{"maxM02", 0.7, "Maximum M02 for EMCal M02 cut"};
    Configurable<float> minClusterE{"minClusterE", 0.7, "Minimum cluster energy for EMCal energy cut"};
    Configurable<int> minNCell{"minNCell", 1, "Minimum number of cells per cluster for EMCal NCell cut"};
    Configurable<std::vector<float>> cfgEMCTMEta{"cfgEMCTMEta", {0.01f, 4.07f, -2.5f}, "|eta| <= [0]+(pT+[1])^[2] for EMCal track matching"};
    Configurable<std::vector<float>> cfgEMCTMPhi{"cfgEMCTMPhi", {0.015f, 3.65f, -2.f}, "|phi| <= [0]+(pT+[1])^[2] for EMCal track matching"};
    Configurable<std::vector<float>> emcSecTMEta{"emcSecTMEta", {0.01f, 4.07f, -2.5f}, "|eta| <= [0]+(pT+[1])^[2] for EMCal track matching"};
    Configurable<std::vector<float>> emcSecTMPhi{"emcSecTMPhi", {0.015f, 3.65f, -2.f}, "|phi| <= [0]+(pT+[1])^[2] for EMCal track matching"};
    Configurable<float> tmEoverP{"tmEoverP", 1.75, "Minimum cluster energy over track momentum for EMCal track matching"};
    Configurable<bool> useExoticCut{"useExoticCut", true, "FLag to use the EMCal exotic cluster cut"};
    Configurable<bool> cfgEMCUseTM{"cfgEMCUseTM", false, "flag to use EMCal track matching cut or not"};
    Configurable<bool> emcUseSecondaryTM{"emcUseSecondaryTM", false, "flag to use EMCal secondary track matching cut or not"};
  } emccuts;

  void defineEMEventCut()
  {
    fEMEventCut = EMPhotonEventCut("fEMEventCut", "fEMEventCut");
    fEMEventCut.SetRequireSel8(eventcuts.cfgRequireSel8);
    fEMEventCut.SetRequireFT0AND(eventcuts.cfgRequireFT0AND);
    fEMEventCut.SetZvtxRange(-eventcuts.cfgZvtxMax, +eventcuts.cfgZvtxMax);
    fEMEventCut.SetRequireNoTFB(eventcuts.cfgRequireNoTFB);
    fEMEventCut.SetRequireNoITSROFB(eventcuts.cfgRequireNoITSROFB);
    fEMEventCut.SetRequireNoSameBunchPileup(eventcuts.cfgRequireNoSameBunchPileup);
    fEMEventCut.SetRequireVertexITSTPC(eventcuts.cfgRequireVertexITSTPC);
    fEMEventCut.SetRequireGoodZvtxFT0vsPV(eventcuts.cfgRequireGoodZvtxFT0vsPV);
    fEMEventCut.SetRequireEMCReadoutInMB(eventcuts.cfgRequireEMCReadoutInMB);
    fEMEventCut.SetRequireEMCHardwareTriggered(eventcuts.cfgRequireEMCHardwareTriggered);
  }

  void defineEMCCut()
  {
    fEMCCut.SetClusterizer(emccuts.clusterDefinition);
    fEMCCut.SetMinE(emccuts.minClusterE);
    fEMCCut.SetMinNCell(emccuts.minNCell);
    fEMCCut.SetM02Range(emccuts.minM02, emccuts.maxM02);
    fEMCCut.SetTimeRange(emccuts.minClusterTime, emccuts.maxClusterTime);

    fEMCCut.SetUseExoticCut(emccuts.useExoticCut);

    fEMCCut.SetTrackMatchingEtaParams(emccuts.cfgEMCTMEta->at(0), emccuts.cfgEMCTMEta->at(1), emccuts.cfgEMCTMEta->at(2));
    fEMCCut.SetTrackMatchingPhiParams(emccuts.cfgEMCTMPhi->at(0), emccuts.cfgEMCTMPhi->at(1), emccuts.cfgEMCTMPhi->at(2));
    fEMCCut.SetMinEoverP(emccuts.tmEoverP);
    fEMCCut.SetUseTM(emccuts.cfgEMCUseTM.value); // disables or enables TM

    fEMCCut.SetSecTrackMatchingEtaParams(emccuts.emcSecTMEta->at(0), emccuts.emcSecTMEta->at(1), emccuts.emcSecTMEta->at(2));
    fEMCCut.SetSecTrackMatchingPhiParams(emccuts.emcSecTMPhi->at(0), emccuts.emcSecTMPhi->at(1), emccuts.emcSecTMPhi->at(2));
    fEMCCut.SetUseSecondaryTM(emccuts.emcUseSecondaryTM.value); // disables or enables secondary TM
  }

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  SliceCache cache;

  PresliceOptional<MinMTracks> perEMCClusterMT = o2::aod::mintm::minClusterId;
  PresliceOptional<MinMSTracks> perEMCClusterMS = o2::aod::mintm::minClusterId;
  PresliceOptional<EMCalPhotons> perCollisionEMC = o2::aod::emccluster::emphotoneventId;

  void init(InitContext&)
  {

    defineEMCCut();
    defineEMEventCut();

    o2::aod::pwgem::photonmeson::utils::eventhistogram::addEventHistograms(&fRegistry);
    auto hEMCCollisionCounter = fRegistry.add<TH1>("Event/hEMCCollisionCounter", "Number of collisions after event cuts", HistType::kTH1D, {{7, 0.5, 7.5}}, false);
    hEMCCollisionCounter->GetXaxis()->SetBinLabel(1, "all");
    hEMCCollisionCounter->GetXaxis()->SetBinLabel(2, "+TVX");         // TVX
    hEMCCollisionCounter->GetXaxis()->SetBinLabel(3, "+|z|<10cm");    // TVX with z < 10cm
    hEMCCollisionCounter->GetXaxis()->SetBinLabel(4, "+Sel8");        // TVX with z < 10cm and Sel8
    hEMCCollisionCounter->GetXaxis()->SetBinLabel(5, "+Good z vtx");  // TVX with z < 10cm and Sel8 and good z xertex
    hEMCCollisionCounter->GetXaxis()->SetBinLabel(6, "+unique");      // TVX with z < 10cm and Sel8 and good z xertex and unique (only collision in the BC)
    hEMCCollisionCounter->GetXaxis()->SetBinLabel(7, "+EMC readout"); // TVX with z < 10cm and Sel8 and good z xertex and unique (only collision in the BC) and kTVXinEMC
    o2::aod::pwgem::photonmeson::utils::clusterhistogram::addClusterHistograms(&fRegistry, cfgDo2DQA);
  }

  bool isEventGood(MyCollision const& collision)
  {

    if (eventcuts.onlyKeepWeightedEvents && std::fabs(collision.weight() - 1.) < 1E-10) {
      return false;
    }

    fRegistry.fill(HIST("Event/hEMCCollisionCounter"), 1, collision.weight());
    if (!eventcuts.cfgRequireFT0AND || collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
      fRegistry.fill(HIST("Event/hEMCCollisionCounter"), 2, collision.weight());
      if (std::abs(collision.posZ()) < eventcuts.cfgZvtxMax) {
        fRegistry.fill(HIST("Event/hEMCCollisionCounter"), 3, collision.weight());
        if (!eventcuts.cfgRequireSel8 || collision.sel8()) {
          fRegistry.fill(HIST("Event/hEMCCollisionCounter"), 4, collision.weight());
          if (!eventcuts.cfgRequireGoodZvtxFT0vsPV || collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
            fRegistry.fill(HIST("Event/hEMCCollisionCounter"), 5, collision.weight());
            if (!eventcuts.cfgRequireNoSameBunchPileup || collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
              fRegistry.fill(HIST("Event/hEMCCollisionCounter"), 6, collision.weight());
              if (!eventcuts.cfgRequireEMCReadoutInMB || collision.alias_bit(kTVXinEMC)) {
                fRegistry.fill(HIST("Event/hEMCCollisionCounter"), 7, collision.weight());
              }
            }
          }
        }
      }
    }

    o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<0>(&fRegistry, collision, collision.weight());
    if (!fEMEventCut.IsSelected(collision)) {
      return false;
    }
    if (!(eventcuts.cfgOccupancyMin <= collision.trackOccupancyInTimeRange() && collision.trackOccupancyInTimeRange() < eventcuts.cfgOccupancyMax)) {
      return false;
    }

    o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<1>(&fRegistry, collision, collision.weight());
    fRegistry.fill(HIST("Event/before/hCollisionCounter"), 12.0, collision.weight()); // accepted
    fRegistry.fill(HIST("Event/after/hCollisionCounter"), 12.0, collision.weight());  // accepted
    return true;
  }

  template <o2::soa::is_table TClusters>
  void doClusterQA(MyCollision const& collision, TClusters const& clusters)
  {
    fRegistry.fill(HIST("Cluster/before/hNgamma"), clusters.size(), collision.weight());
    int ngAfter = 0;
    for (const auto& cluster : clusters) {
      // Fill the cluster properties before applying any cuts
      if (!fEMCCut.IsSelectedEMCal(EMCPhotonCut::EMCPhotonCuts::kDefinition, cluster)) {
        continue;
      }
      o2::aod::pwgem::photonmeson::utils::clusterhistogram::fillClusterHistograms<0>(&fRegistry, cluster, cfgDo2DQA, collision.weight());

      // Apply cuts one by one and fill in hClusterQualityCuts histogram
      fRegistry.fill(HIST("Cluster/hClusterQualityCuts"), 0., cluster.e(), collision.weight());

      // Define two boleans to see, whether the cluster "survives" the EMC cluster cuts to later check, whether the cuts in this task align with the ones in EMCPhotonCut.h:
      bool survivesIsSelectedEMCalCuts = true; // Survives "manual" cuts listed in this task
      bool survivesIsSelectedCuts = false;

      survivesIsSelectedCuts = fEMCCut.IsSelected(cluster);
      for (int icut = 0; icut < static_cast<int>(EMCPhotonCut::EMCPhotonCuts::kNCuts); icut++) { // Loop through different cut observables
        EMCPhotonCut::EMCPhotonCuts specificcut = static_cast<EMCPhotonCut::EMCPhotonCuts>(icut);
        if (!fEMCCut.IsSelectedEMCal(specificcut, cluster)) { // Check whether cluster passes this cluster requirement, if not, fill why in the next row
          fRegistry.fill(HIST("Cluster/hClusterQualityCuts"), icut + 1, cluster.e(), collision.weight());
          survivesIsSelectedEMCalCuts = false;
        }
      }

      if (survivesIsSelectedCuts != survivesIsSelectedEMCalCuts) {
        LOGF(info, "Cummulative application of IsSelectedEMCal cuts does not equal the IsSelected result");
        LOGF(info, "Cummulative application of IsSelectedEMCal resulted in %d", survivesIsSelectedEMCalCuts);
        LOGF(info, "IsSelected resulted in %d", survivesIsSelectedCuts);
      }

      if (survivesIsSelectedCuts) {
        o2::aod::pwgem::photonmeson::utils::clusterhistogram::fillClusterHistograms<1>(&fRegistry, cluster, cfgDo2DQA, collision.weight());
        fRegistry.fill(HIST("Cluster/hClusterQualityCuts"), static_cast<double>(EMCPhotonCut::EMCPhotonCuts::kNCuts) + 1., cluster.e(), collision.weight());
        ++ngAfter;
      }
    }
    fRegistry.fill(HIST("Cluster/after/hNgamma"), ngAfter, collision.weight());

    return;
  }

  template <o2::soa::is_table TClusters, o2::soa::is_table TMatchedTracks, o2::soa::is_table TMatchedSecondaries>
  void doClusterQA(MyCollision const& collision, TClusters const& clusters, TMatchedTracks const& primTracks, TMatchedSecondaries const& secTracks)
  {
    fRegistry.fill(HIST("Cluster/before/hNgamma"), clusters.size(), collision.weight());
    int ngAfter = 0;
    for (const auto& cluster : clusters) {
      // Fill the cluster properties before applying any cuts
      if (!fEMCCut.IsSelectedEMCal(EMCPhotonCut::EMCPhotonCuts::kDefinition, cluster)) {
        continue;
      }
      auto primTracksPerCluster = primTracks.sliceBy(perEMCClusterMT, cluster.globalIndex());
      auto secTracksPerCluster = secTracks.sliceBy(perEMCClusterMS, cluster.globalIndex());
      o2::aod::pwgem::photonmeson::utils::clusterhistogram::fillClusterHistograms<0>(&fRegistry, cluster, cfgDo2DQA, collision.weight(), primTracksPerCluster, secTracksPerCluster);

      // Apply cuts one by one and fill in hClusterQualityCuts histogram
      fRegistry.fill(HIST("Cluster/hClusterQualityCuts"), 0., cluster.e(), collision.weight());

      // Define two boleans to see, whether the cluster "survives" the EMC cluster cuts to later check, whether the cuts in this task align with the ones in EMCPhotonCut.h:
      bool survivesIsSelectedEMCalCuts = true; // Survives "manual" cuts listed in this task
      bool survivesIsSelectedCuts = false;

      survivesIsSelectedCuts = fEMCCut.IsSelected(cluster, primTracksPerCluster, secTracksPerCluster);
      for (int icut = 0; icut < static_cast<int>(EMCPhotonCut::EMCPhotonCuts::kNCuts); icut++) { // Loop through different cut observables
        EMCPhotonCut::EMCPhotonCuts specificcut = static_cast<EMCPhotonCut::EMCPhotonCuts>(icut);
        if (specificcut == EMCPhotonCut::EMCPhotonCuts::kTM || specificcut == EMCPhotonCut::EMCPhotonCuts::kSecondaryTM) {
          // will do track matching cuts extra later or never depending on chosen process function
          continue;
        }
        if (!fEMCCut.IsSelectedEMCal(specificcut, cluster)) { // Check whether cluster passes this cluster requirement, if not, fill why in the next row
          fRegistry.fill(HIST("Cluster/hClusterQualityCuts"), icut + 1, cluster.e(), collision.weight());
          survivesIsSelectedEMCalCuts = false;
        }
      }

      if (!fEMCCut.IsSelectedEMCal(EMCPhotonCut::EMCPhotonCuts::kTM, cluster, primTracksPerCluster)) { // Check whether cluster passes this cluster requirement, if not, fill why in the next row
        fRegistry.fill(HIST("Cluster/hClusterQualityCuts"), static_cast<double>(EMCPhotonCut::EMCPhotonCuts::kTM) + 1, cluster.e(), collision.weight());
        survivesIsSelectedEMCalCuts = false;
      }
      if (!fEMCCut.IsSelectedEMCal(EMCPhotonCut::EMCPhotonCuts::kSecondaryTM, cluster, secTracksPerCluster)) { // Check whether cluster passes this cluster requirement, if not, fill why in the next row
        fRegistry.fill(HIST("Cluster/hClusterQualityCuts"), static_cast<double>(EMCPhotonCut::EMCPhotonCuts::kSecondaryTM) + 1, cluster.e(), collision.weight());
        survivesIsSelectedEMCalCuts = false;
      }

      if (survivesIsSelectedCuts != survivesIsSelectedEMCalCuts) [[unlikely]] {
        LOGF(info, "Cummulative application of IsSelectedEMCal cuts does not equal the IsSelected result");
        LOGF(info, "Cummulative application of IsSelectedEMCal resulted in %d", survivesIsSelectedEMCalCuts);
        LOGF(info, "IsSelected resulted in %d", survivesIsSelectedCuts);
      }

      if (survivesIsSelectedCuts) {
        o2::aod::pwgem::photonmeson::utils::clusterhistogram::fillClusterHistograms<1>(&fRegistry, cluster, cfgDo2DQA, collision.weight(), primTracksPerCluster, secTracksPerCluster);
        fRegistry.fill(HIST("Cluster/hClusterQualityCuts"), static_cast<double>(EMCPhotonCut::EMCPhotonCuts::kNCuts) + 1., cluster.e(), collision.weight());
        ++ngAfter;
      }
    }
    fRegistry.fill(HIST("Cluster/after/hNgamma"), ngAfter, collision.weight());

    return;
  }

  void processQC(MyCollisions const& collisions, EMCalPhotons const& clusters)
  {
    for (const auto& collision : collisions) {

      if (!isEventGood(collision)) {
        continue;
      }

      auto clustersPerColl = clusters.sliceBy(perCollisionEMC, collision.collisionId());
      doClusterQA(collision, clustersPerColl);

    } // end of collision loop
  } // end of process

  void processQCTM(MyCollisions const& collisions, EMCalPhotons const& clusters, MinMTracks const& matchedPrims, MinMSTracks const& matchedSeconds)
  {
    for (const auto& collision : collisions) {

      if (!isEventGood(collision)) {
        continue;
      }

      auto clustersPerColl = clusters.sliceBy(perCollisionEMC, collision.collisionId());
      doClusterQA(collision, clustersPerColl, matchedPrims, matchedSeconds);

    } // end of collision loop
  } // end of process

  void processDummy(MyCollisions const&) {}

  PROCESS_SWITCH(EmcalQC, processQC, "run EMCal QC", false);
  PROCESS_SWITCH(EmcalQC, processQCTM, "run EMCal QC with track matching", false);
  PROCESS_SWITCH(EmcalQC, processDummy, "Dummy function", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<EmcalQC>(cfgc)};
}
