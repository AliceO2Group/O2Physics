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
// Task performing Lund declustering and producing primary Lund plane
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "fastjet/contrib/LundGenerator.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/runDataProcessing.h"
#include "TDatabasePDG.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Common/Core/RecoDecay.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/FastJetUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace std;

struct JetLundPlane {

  HistogramRegistry histos{"LundHistos"};

  std::vector<fastjet::PseudoJet> jetConstituents;
  std::vector<fastjet::PseudoJet> jetReclustered;
  JetFinder jetReclusterer;

  void init(InitContext const&)
  {
    histos.add("PrimaryLundPlane3D", "Primary Lund 3D plane;ln(R/Delta);ln(k_{t}/GeV);{p}_{t}", {HistType::kTH3F, {{100, 0, 20}, {100, -10, 20}, {20, 0, 200}}});
    histos.add("PrimaryLundPlane2D", "Primary Lund 2D plane;ln(R/Delta);ln(k_{t}/GeV)", {HistType::kTH2F, {{100, 0, 20}, {100, -10, 20}}});
    histos.print();
    jetReclusterer.isReclustering = true;
    jetReclusterer.algorithm = fastjet::JetAlgorithm::cambridge_algorithm;
  }

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
        std::swap(j1,j2);
      }
      double deltaR = j1.delta_R(j2);
      double kt = j2.pt() * deltaR;
      double coord1 = std::log(jet.r()/deltaR);
      double coord2 = std::log(kt);
      histos.fill(HIST("PrimaryLundPlane3D"), coord1, coord2, jet.pt());
      histos.fill(HIST("PrimaryLundPlane2D"), coord1, coord2);
      pair = j1;
    }
  }

  void processDummy(aod::Tracks const& track)
  {
  }
  PROCESS_SWITCH(JetLundPlane, processDummy, "Dummy process function, turned on by default", true);
  
  void processChargedJets(soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>::iterator const& jet, aod::Tracks const& tracks)
  {
    jetConstituents.clear();
    for (auto& jetConstituent : jet.tracks_as<aod::Tracks>()) {
      FastJetUtilities::fillTracks(jetConstituent, jetConstituents, jetConstituent.globalIndex());
    }
    jetReclustering(jet);
  }
  PROCESS_SWITCH(JetLundPlane, processChargedJets, "Process function for charged jets", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<JetLundPlane>(cfgc,
                                                    SetDefaultProcesses{},
                                                    TaskName{"jet-lund-declustering"}));
  return WorkflowSpec{tasks};
}
