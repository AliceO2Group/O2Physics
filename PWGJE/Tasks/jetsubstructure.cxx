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

// jet analysis tasks (subscribing to jet finder task)
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>
//

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/O2DatabasePDGPlugin.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetSubstructure.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/FastJetUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"

struct JetSubstructureTask {
  Produces<aod::ChargedJetSubstructures> jetSubstructureDataTable;
  Produces<aod::ChargedMCDetectorLevelJetSubstructures> jetSubstructureMCDTable;
  Produces<aod::ChargedMCParticleLevelJetSubstructures> jetSubstructureMCPTable;
  Produces<aod::ChargedEventWiseSubtractedJetSubstructures> jetSubstructureDataSubTable;
  OutputObj<TH2F> hZg{"h_jet_zg_jet_pt"};
  OutputObj<TH2F> hRg{"h_jet_rg_jet_pt"};
  OutputObj<TH2F> hNsd{"h_jet_nsd_jet_pt"};

  Configurable<float> zCut{"zCut", 0.1, "soft drop z cut"};
  Configurable<float> beta{"beta", 0.0, "soft drop beta"};

  Service<o2::framework::O2DatabasePDG> pdg;
  std::vector<fastjet::PseudoJet> jetConstituents;
  std::vector<fastjet::PseudoJet> jetReclustered;
  JetFinder jetReclusterer;

  void init(InitContext const&)
  {
    hZg.setObject(new TH2F("h_jet_zg_jet_pt", ";z_{g}; #it{p}_{T,jet} (GeV/#it{c})",
                           10, 0.0, 0.5, 200, 0.0, 200.0));
    hRg.setObject(new TH2F("h_jet_rg_jet_pt", ";R_{g}; #it{p}_{T,jet} (GeV/#it{c})",
                           10, 0.0, 0.5, 200, 0.0, 200.0));
    hNsd.setObject(new TH2F("h_jet_nsd_jet_pt", ";n_{SD}; #it{p}_{T,jet} (GeV/#it{c})",
                            7, -0.5, 6.5, 200, 0.0, 200.0));

    jetReclusterer.isReclustering = true;
    jetReclusterer.algorithm = fastjet::JetAlgorithm::cambridge_algorithm;
  }

  template <typename T, typename U>
  void jetReclustering(T const& jet, U& outputTable)
  {
    jetReclustered.clear();
    fastjet::ClusterSequenceArea clusterSeq(jetReclusterer.findJets(jetConstituents, jetReclustered));
    jetReclustered = sorted_by_pt(jetReclustered);
    fastjet::PseudoJet daughterSubJet = jetReclustered[0];
    fastjet::PseudoJet parentSubJet1;
    fastjet::PseudoJet parentSubJet2;
    bool softDropped = false;
    auto nsd = 0.0;
    auto zg = -1.0;
    auto rg = -1.0;
    while (daughterSubJet.has_parents(parentSubJet1, parentSubJet2)) {
      if (parentSubJet1.perp() < parentSubJet2.perp()) {
        std::swap(parentSubJet1, parentSubJet2);
      }
      auto z = parentSubJet2.perp() / (parentSubJet1.perp() + parentSubJet2.perp());
      auto theta = parentSubJet1.delta_R(parentSubJet2);
      if (z >= zCut * TMath::Power(theta / (jet.r() / 100.f), beta)) {
        if (!softDropped) {
          zg = z;
          rg = theta;
          hZg->Fill(zg, jet.pt());
          hRg->Fill(rg, jet.pt());
          softDropped = true;
        }
        nsd++;
      }
      daughterSubJet = parentSubJet1;
    }
    hNsd->Fill(nsd, jet.pt());
    outputTable(zg, rg, nsd);
  }

  template <typename T, typename U, typename V>
  void analyseCharged(T const& jet, U const& tracks, V& outputTable)
  {
    jetConstituents.clear();
    for (auto& jetConstituent : jet.template tracks_as<U>()) {
      fastjetutilities::fillTracks(jetConstituent, jetConstituents, jetConstituent.globalIndex());
    }
    jetReclustering(jet, outputTable);
  }

  void processDummy(JetTracks const& tracks)
  {
  }
  PROCESS_SWITCH(JetSubstructureTask, processDummy, "Dummy process function turned on by default", true);

  void processChargedJetsData(soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>::iterator const& jet,
                              JetTracks const& tracks)
  {
    analyseCharged(jet, tracks, jetSubstructureDataTable);
  }
  PROCESS_SWITCH(JetSubstructureTask, processChargedJetsData, "charged jet substructure", false);

  void processChargedJetsEventWiseSubData(soa::Join<aod::ChargedEventWiseSubtractedJets, aod::ChargedEventWiseSubtractedJetConstituents>::iterator const& jet,
                                          JetTracksSub const& tracks)
  {
    analyseCharged(jet, tracks, jetSubstructureDataSubTable);
  }
  PROCESS_SWITCH(JetSubstructureTask, processChargedJetsEventWiseSubData, "eventwise-constituent subtracted charged jet substructure", false);

  void processChargedJetsMCD(typename soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>::iterator const& jet,
                             JetTracks const& tracks)
  {
    analyseCharged(jet, tracks, jetSubstructureMCDTable);
  }
  PROCESS_SWITCH(JetSubstructureTask, processChargedJetsMCD, "charged jet substructure", false);

  void processChargedJetsMCP(typename soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents>::iterator const& jet,
                             JetParticles const& particles)
  {
    jetConstituents.clear();
    for (auto& jetConstituent : jet.template tracks_as<JetParticles>()) {
      fastjetutilities::fillTracks(jetConstituent, jetConstituents, jetConstituent.globalIndex(), static_cast<int>(JetConstituentStatus::track), pdg->Mass(jetConstituent.pdgCode()));
    }
    jetReclustering(jet, jetSubstructureMCPTable);
  }
  PROCESS_SWITCH(JetSubstructureTask, processChargedJetsMCP, "charged jet substructure on MC particle level", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{

  return WorkflowSpec{adaptAnalysisTask<JetSubstructureTask>(
    cfgc, TaskName{"jet-substructure"})};
}
