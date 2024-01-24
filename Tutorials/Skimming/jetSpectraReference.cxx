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

#include "TH1F.h"
#include "TTree.h"

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/HistogramRegistry.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/Core/JetFinder.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetSpectraReference {

  HistogramRegistry registry{
    "registry",
    {{"hJetPt", "Jet pT;Jet #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 100.}}}},
     {"hNJetConstituents", "Number of constituents;N;entries", {HistType::kTH1F, {{100, -0.5, 99.5}}}},
     {"hConstituentPt", "Constituent pT; Constituent #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 100.}}}}}};

  // Configurable<float> f_jetPtMin{"f_jetPtMin", 0.0, "minimum jet pT cut"};
  // Filter jetCuts = aod::jet::pt > f_jetPtMin; //how does this work?

  void process(soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>::iterator const& jet,
               aod::Tracks const& tracks)
  {
    registry.fill(HIST("hJetPt"), jet.pt());
    registry.fill(HIST("hNJetConstituents"), jet.tracks().size());
    for (auto& constituent : jet.tracks_as<aod::Tracks>()) {
      registry.fill(HIST("hConstituentPt"), constituent.pt());
    }
  }
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetSpectraReference>(cfgc, TaskName{"jetspectra-task-skim-reference"})};
}
