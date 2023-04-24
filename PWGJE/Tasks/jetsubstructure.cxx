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
// Author: Nima Zardoshti
//

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "TDatabasePDG.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Common/Core/RecoDecay.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetSubstructure.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/FastJetUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"

template <typename SubstructureTable>
struct JetSubstructureTask {
  Produces<SubstructureTable> jetSubstructureTable;
  OutputObj<TH1F> hZg{"h_jet_zg"};
  OutputObj<TH1F> hRg{"h_jet_rg"};
  OutputObj<TH1F> hNsd{"h_jet_nsd"};

  Configurable<float> jetPtMin{"jetPtMin", 0.0, "minimum jet pT cut"};
  Configurable<float> zCut{"zCut", 0.1, "soft drop z cut"};
  Configurable<float> beta{"beta", 0.0, "soft drop beta"};
  Configurable<float> jetR{"jetR", 0.4, "jet resolution parameter"};
  Configurable<bool> doConstSub{"doConstSub", false, "do constituent subtraction"};

  std::vector<fastjet::PseudoJet> jetConstituents;
  std::vector<fastjet::PseudoJet> jetReclustered;
  JetFinder jetReclusterer;

  void init(InitContext const&)
  {
    hZg.setObject(new TH1F("h_jet_zg", "zg ;zg",
                           10, 0.0, 0.5));
    hRg.setObject(new TH1F("h_jet_rg", "rg ;rg",
                           10, 0.0, 0.5));
    hNsd.setObject(new TH1F("h_jet_nsd", "nsd ;nsd",
                            7, -0.5, 6.5));
    jetReclusterer.isReclustering = true;
    jetReclusterer.algorithm = fastjet::JetAlgorithm::cambridge_algorithm;
  }

  //Filter jetCuts = aod::jet::pt > f_jetPtMin; //how does this work?

  template <typename T>
  void jetReclustering(T const& jet)
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
      if (z >= zCut * TMath::Power(theta / jetR, beta)) {
        if (!softDropped) {
          zg = z;
          rg = theta;
          hZg->Fill(zg);
          hRg->Fill(rg);
          softDropped = true;
        }
        nsd++;
      }
      daughterSubJet = parentSubJet1;
    }
    hNsd->Fill(nsd);
    jetSubstructureTable(jet.globalIndex(), zg, rg, nsd);
  }

  void processDummy(aod::Tracks const& track)
  {
  }
  PROCESS_SWITCH(JetSubstructureTask, processDummy, "Dummy process function turned on by default", true);

  void processData(soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>::iterator const& jet,
                   aod::Tracks const& tracks,
                   aod::ChargedJetConstituentsSub const& constituentsSub)
  {
    jetConstituents.clear();

    if (doConstSub) {
      for (const auto& jetconstituentSub : constituentsSub) {
        FastJetUtilities::fillTracks(jetconstituentSub, jetConstituents);
      }
    } else {
      for (auto& jetConstituent : jet.tracks_as<aod::Tracks>()) {
        FastJetUtilities::fillTracks(jetConstituent, jetConstituents, jetConstituent.globalIndex());
      }
    }
    jetReclustering(jet);
  }
  PROCESS_SWITCH(JetSubstructureTask, processData, "jet substructure on data", false);

  void processMCD(soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>::iterator const& jet,
                  aod::Tracks const& tracks,
                  aod::ChargedJetConstituentsSub const& constituentsSub)
  {
    jetConstituents.clear();

    if (doConstSub) {
      for (const auto& jetconstituentSub : constituentsSub) {
        FastJetUtilities::fillTracks(jetconstituentSub, jetConstituents);
      }
    } else {
      for (auto& jetConstituent : jet.tracks_as<aod::Tracks>()) {
        FastJetUtilities::fillTracks(jetConstituent, jetConstituents, jetConstituent.globalIndex());
      }
    }
    jetReclustering(jet);
  }
  PROCESS_SWITCH(JetSubstructureTask, processMCD, "jet substructure on MC detector level", false);

  void processMCP(soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents>::iterator const& jet,
                  aod::McParticles const& particles)
  {
    jetConstituents.clear();
    for (auto& jetConstituent : jet.tracks_as<aod::McParticles>()) {
      FastJetUtilities::fillTracks(jetConstituent, jetConstituents, jetConstituent.globalIndex(), static_cast<int>(JetConstituentStatus::track), RecoDecay::getMassPDG(jetConstituent.pdgCode()));
    }
    jetReclustering(jet);
  }
  PROCESS_SWITCH(JetSubstructureTask, processMCP, "jet substructure on MC particle level", false);
};
using JetSubstructureDataLevel = JetSubstructureTask<o2::aod::ChargedJetSubstructure>;
using JetSubstructureMCDetectorLevel = JetSubstructureTask<o2::aod::ChargedMCDetectorLevelJetSubstructure>;
using JetSubstructureMCParticleLevel = JetSubstructureTask<o2::aod::ChargedMCParticleLevelJetSubstructure>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(adaptAnalysisTask<JetSubstructureDataLevel>(cfgc,
                                                                 SetDefaultProcesses{},
                                                                 TaskName{"jet-substructure-data"}));

  tasks.emplace_back(adaptAnalysisTask<JetSubstructureMCDetectorLevel>(cfgc,
                                                                       SetDefaultProcesses{},
                                                                       TaskName{"jet-substructure-mcd"}));

  tasks.emplace_back(adaptAnalysisTask<JetSubstructureMCParticleLevel>(cfgc,
                                                                       SetDefaultProcesses{},
                                                                       TaskName{"jet-substructure-mcp"}));

  return WorkflowSpec{tasks};
}
