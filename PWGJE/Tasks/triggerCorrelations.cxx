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

// QA correlation task for jet trigger
//
// Authors: Nima Zardoshti

#include <string>

#include "EMCALBase/Geometry.h"
#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/HistogramRegistry.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/CCDB/TriggerAliases.h"

#include "Common/Core/RecoDecay.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/DataModel/Jet.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"

struct TriggerCorrelationsTask {

  HistogramRegistry registry{"registry",
                             {{"h3_jetcharged_pt_jetfull_pt_photon_pt_triggered", ";#it{p}_{T,jet}^{charged} (GeV/#it{c}); #it{p}_{T,jet}^{full} (GeV/#it{c}); #it{p}_{T,photon} (GeV/#it{c})", {HistType::kTH3F, {{201, -1.0, 200}, {201, -1.0, 200}, {201, -1.0, 200}}}, true},
                              {"h_collision_trigger_events", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}}, true}}};

  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<bool> doChargedJetTrigger{"doChargedJetTrigger", true, "add the charged jeet trigger to the QA"};
  Configurable<bool> doFullJetTrigger{"doFullJetTrigger", true, "add the full et trigger to the QA"};
  Configurable<bool> doGammaTrigger{"doGammaTrigger", true, "add the gamma trigger to the QA"};
  Configurable<float> jetsChargedR{"jetsChargedR", 0.6, "resolution parameter for triggered charged jets"};
  Configurable<float> jetsFullR{"jetsFullR", 0.2, "resolution parameter for triggered full jets"};
  Configurable<bool> rejectExoticClusters{"rejectExoticClusters", true, "reject exotic clusters"};
  Configurable<float> clusterTimeMin{"clusterTimeMin", -999, "minimum cluster time for gamma trigger"};
  Configurable<float> clusterTimeMax{"clusterTimeMax", 999, "maximum cluster time for gamma trigger"};
  Configurable<int> emcalTriggered{"emcalTriggered", -1, "-1 = min bias, 0 = L0"};
  Configurable<std::string> clusterDefinition{"clusterDefinition", "kV3Default", "cluster definition to be selected, e.g. V3Default"};

  void init(o2::framework::InitContext&)
  {
  }

  o2::aod::EMCALClusterDefinition clusterDef = o2::aod::emcalcluster::getClusterDefinitionFromString(clusterDefinition.value);
  Filter clusterDefinitionSelection = o2::aod::emcalcluster::definition == static_cast<int>(clusterDef);

  void processTriggeredCorrelations(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
                                    soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jetsCharged,
                                    soa::Join<aod::FullJets, aod::FullJetConstituents> const& jetsFull,
                                    soa::Filtered<o2::aod::EMCALClusters> const& clusters,
                                    soa::Join<aod::Tracks, aod::TrackSelection> const& tracks)
  {

    registry.fill(HIST("h_collision_trigger_events"), 0.5); // all events
    if (collision.posZ() < vertexZCut) {
      registry.fill(HIST("h_collision_trigger_events"), 1.5); // all events with z vertex cut
    }
    if (collision.sel8()) {
      registry.fill(HIST("h_collision_trigger_events"), 2.5); // events with sel8()
    }
    if (collision.alias_bit(kTVXinEMC)) {
      registry.fill(HIST("h_collision_trigger_events"), 3.5); // events with emcal bit
    }

    float jetsChargedLeadingPt = -1.0;
    float jetsFullLeadingPt = -1.0;
    float gammaLeadingPt = -1.0;
    if (collision.sel8()) {
      if (doChargedJetTrigger) {
        for (auto& jetCharged : jetsCharged) {
          if (jetCharged.r() == round(jetsChargedR * 100.0f)) {
            jetsChargedLeadingPt = jetCharged.pt();
            break;
          }
        }
      }
    }

    if ((emcalTriggered == -1 && collision.alias_bit(kTVXinEMC)) || (emcalTriggered == 0 && (collision.alias_bit(kEMC7) || collision.alias_bit(kDMC7)))) {
      if (doFullJetTrigger) {
        for (auto& jetFull : jetsFull) {
          if (jetFull.r() == round(jetsFullR * 100.0f)) {
            jetsFullLeadingPt = jetFull.pt();
            break;
          }
        }
      }

      if (doGammaTrigger) {
        for (const auto& cluster : clusters) {
          if (rejectExoticClusters && cluster.isExotic()) {
            continue;
          }
          if (cluster.time() < clusterTimeMin || cluster.time() > clusterTimeMax) {
            continue;
          }
          double gammaPt = cluster.energy() / std::cosh(cluster.eta());
          if (TVector2::Phi_0_2pi(cluster.phi()) < 4 && gammaPt > gammaLeadingPt) {
            gammaLeadingPt = gammaPt;
          }
        }
      }
    }

    float jetsChargedPtStart = 0.0;
    float jetsFullPtStart = 0.0;
    float gammaPtStart = 0.0;
    if (jetsChargedLeadingPt < 0.0) {
      jetsChargedPtStart = -1.0;
    }
    if (jetsFullLeadingPt < 0.0) {
      jetsFullPtStart = -1.0;
    }
    if (gammaLeadingPt < 0.0) {
      gammaPtStart = -1.0;
    }

    for (float jetChargedPt = jetsChargedPtStart; jetChargedPt <= jetsChargedLeadingPt; jetChargedPt += 1.0) {

      for (float jetFullPt = jetsFullPtStart; jetFullPt <= jetsFullLeadingPt; jetFullPt += 1.0) {

        for (float gammaPt = gammaPtStart; gammaPt <= gammaLeadingPt; gammaPt += 1.0) {

          registry.fill(HIST("h3_jetcharged_pt_jetfull_pt_photon_pt_triggered"), jetChargedPt, jetFullPt, gammaPt);
        }
      }
    }
  }
  PROCESS_SWITCH(TriggerCorrelationsTask, processTriggeredCorrelations, "QA for trigger correlations", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<TriggerCorrelationsTask>(cfgc, TaskName{"trigger-correlations"})};
}
