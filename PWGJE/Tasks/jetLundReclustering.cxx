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

// \author
// Alice Caluisi   -   alice.caluisi@cern.ch
// \since September 2023

//
// Task performing jet reclustering and producing primary Lund Plane histograms
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "fastjet/contrib/LundGenerator.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "TDatabasePDG.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Common/Core/RecoDecay.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/DataModel/JetSubstructure.h"
#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"

using namespace o2;
using namespace o2::track;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace std;

struct JetLundReclustering {

  HistogramRegistry registry;

  std::vector<fastjet::PseudoJet> jetConstituents;
  std::vector<fastjet::PseudoJet> jetReclustered;
  JetFinder jetReclusterer;

  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<float> jetPtMin{"jetPtMin", 5.0, "minimum jet pT cut"};
  Configurable<float> jetR{"jetR", 0.4, "jet resolution parameter"};
  Configurable<float> jet_min_eta{"jet_min_eta", -0.5, "minimum jet eta"};
  Configurable<float> jet_max_eta{"jet_max_eta", 0.5, "maximum jet eta"};
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};

  std::vector<int> eventSelectionBits;

  void init(InitContext const&)
  {
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));

    registry.add("PrimaryLundPlane_kT", "Primary Lund 3D plane;ln(R/Delta);ln(k_{t}/GeV);{p}_{t}", {HistType::kTH3F, {{100, 0, 10}, {100, -10, 10}, {20, 0, 200}}});
    registry.add("PrimaryLundPlane_z", "Primary Lund 3D plane;ln(R/Delta);ln(1/z);{p}_{t}", {HistType::kTH3F, {{100, 0, 10}, {100, 0, 10}, {20, 0, 200}}});
    registry.add("jet_PtEtaPhi", "Correlation of jet #it{p}_{T}, #eta and #phi;#it{p}_{T,jet} (GeV/#it{c});#eta_{jet};#phi_{jet} [rad]", {HistType::kTH3F, {{100, 0, 200}, {180, -0.9, 0.9}, {180, 0., 2 * M_PI}}});

    jetReclusterer.isReclustering = true;
    jetReclusterer.algorithm = fastjet::JetAlgorithm::cambridge_algorithm;
  }

  Filter jetFilter = aod::jet::pt > jetPtMin&& aod::jet::r == nround(jetR.node() * 100.0f) && aod::jet::eta > jet_min_eta&& aod::jet::eta < jet_max_eta;
  Filter collisionFilter = nabs(aod::jcollision::posZ) < vertexZCut;

  // Reclustering function
  template <typename T>
  void jetReclustering(T const& jet)
  {
    jetReclustered.clear();
    fastjet::ClusterSequenceArea clusterSeq(jetReclusterer.findJets(jetConstituents, jetReclustered));
    jetReclustered = sorted_by_pt(jetReclustered);
    fastjet::PseudoJet pair = jetReclustered[0];
    fastjet::PseudoJet j1;
    fastjet::PseudoJet j2;
    while (pair.has_parents(j1, j2)) {
      if (j1.pt() < j2.pt()) {
        std::swap(j1, j2);
      }
      double deltaR = j1.delta_R(j2);
      double kt = j2.pt() * deltaR;
      double z = j2.pt() / (j1.pt() + j2.pt());
      double jetRadius = static_cast<double>(jet.r()) / 100.0;
      double coord1 = std::log(jetRadius / deltaR);
      double coord2 = std::log(kt);
      double coord3 = std::log(1 / z);
      registry.fill(HIST("PrimaryLundPlane_kT"), coord1, coord2, jet.pt());
      registry.fill(HIST("PrimaryLundPlane_z"), coord1, coord3, jet.pt());
      pair = j1;
    }
  }

  // Dummy process
  void processDummy(aod::JetCollisions const&)
  {
  }
  PROCESS_SWITCH(JetLundReclustering, processDummy, "Dummy process function, turned on by default", true);

  // Process function for charged jets
  void processChargedJets(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                          soa::Filtered<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>> const& jets,
                          aod::JetTracks const&)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    for (const auto& jet : jets) {
      registry.fill(HIST("jet_PtEtaPhi"), jet.pt(), jet.eta(), jet.phi());
      jetConstituents.clear();
      for (auto& jetConstituent : jet.tracks_as<aod::JetTracks>()) {
        fastjetutilities::fillTracks(jetConstituent, jetConstituents, jetConstituent.globalIndex());
      }
      // Perform jet reclustering
      jetReclustering(jet);
    }
  }
  PROCESS_SWITCH(JetLundReclustering, processChargedJets, "Process function for charged jets", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<JetLundReclustering>(cfgc, TaskName{"jet-lund-reclustering"})};
}
