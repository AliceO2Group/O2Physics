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

// Author: Nima Zardoshti

#include <TMath.h>
#include <cmath>
#include <string>

#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetBkgSubUtils.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/DataModel/Jet.h"

#include "../filterTables.h"

#include "Framework/HistogramRegistry.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetHFFilterTask {

  HistogramRegistry registry;

  enum { kJetD0ChLowPt = 0,
         kJetD0ChHighPt = 1,
         kJetLcChLowPt = 2,
         kJetLcChHighPt = 3,
         kAllObjects = 4 };

  Produces<aod::JetHFFilters> tags;

  Configurable<float> jetD0ChLowPtThreshold{"jetD0ChLowPtThreshold", 0.0, "Threshold for charged D0 jet low pt trigger"};
  Configurable<float> jetD0ChHighPtThreshold{"jetD0ChHighPtThreshold", 0.0, "Threshold for charged D0 jet high pt trigger"};
  Configurable<float> jetD0ChR{"jetD0ChR", 0.6, "jet resolution parameter for charged D0 jet for low pt trigger"};
  Configurable<float> jetLcChLowPtThreshold{"jetLcChLowPtThreshold", 0.0, "Threshold for charged Lc jet low pt trigger"};
  Configurable<float> jetLcChHighPtThreshold{"jetLcChHighPtThreshold", 0.0, "Threshold for charged Lc jet high pt trigger"};
  Configurable<float> jetLcChR{"jetLcChR", 0.6, "jet resolution parameter for charged Lc jet for low pt trigger"};

  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};

  void init(o2::framework::InitContext&)
  {
    registry.add("h_d0jet_pt", "D^{0} - tagged jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}});
    registry.add("h_d0jet_pt_lowpt", "D^{0} - tagged jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}});
    registry.add("h_d0jet_pt_highpt", "D^{0} - tagged jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}});

    registry.add("h_lcjet_pt", "#Lambda^{+}_{c} - tagged jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}});
    registry.add("h_lcjet_pt_lowpt", "#Lambda^{+}_{c} - tagged jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}});
    registry.add("h_lcjet_pt_highpt", "#Lambda^{+}_{c} - tagged jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}});

    registry.add("h_collisions", "Collision ;entries", {HistType::kTH1F, {{5, 0.0, 5.0}}});
  }

  void processJets(soa::Join<aod::JetCollisions, aod::EvSels>::iterator const& /*collision*/, soa::Join<o2::aod::D0ChargedJets, o2::aod::D0ChargedJetConstituents> const& d0Jets, aod::CandidatesD0Data const& /*d0Candidates*/, soa::Join<o2::aod::LcChargedJets, o2::aod::LcChargedJetConstituents> const& lcJets, aod::CandidatesLcData const& /*lcCandidates*/, aod::JetTracks const& /*tracks*/)
  {
    registry.fill(HIST("h_collisions"), 0.5);
    bool keepEvent[kAllObjects]{false};
    for (auto const& d0Jet : d0Jets) {
      if (d0Jet.r() == round(jetD0ChR * 100.0f)) {
        registry.fill(HIST("h_d0jet_pt"), d0Jet.pt());
        if (d0Jet.pt() >= jetD0ChLowPtThreshold) {
          keepEvent[kJetD0ChLowPt] = true;
          registry.fill(HIST("h_d0jet_pt_lowpt"), d0Jet.pt());
          registry.fill(HIST("h_collisions"), 1.5);
        }
        if (d0Jet.pt() >= jetD0ChHighPtThreshold) {
          keepEvent[kJetD0ChHighPt] = true;
          registry.fill(HIST("h_d0jet_pt_highpt"), d0Jet.pt());
          registry.fill(HIST("h_collisions"), 2.5);
        }
      }
    }
    for (auto const& lcJet : lcJets) {
      if (lcJet.r() == round(jetLcChR * 100.0f)) {
        registry.fill(HIST("h_lcjet_pt"), lcJet.pt());
        if (lcJet.pt() >= jetLcChLowPtThreshold) {
          keepEvent[kJetLcChLowPt] = true;
          registry.fill(HIST("h_lcjet_pt_lowpt"), lcJet.pt());
          registry.fill(HIST("h_collisions"), 3.5);
        }
        if (lcJet.pt() >= jetLcChHighPtThreshold) {
          keepEvent[kJetLcChHighPt] = true;
          registry.fill(HIST("h_lcjet_pt_highpt"), lcJet.pt());
          registry.fill(HIST("h_collisions"), 4.5);
        }
      }
    }
    tags(keepEvent[kJetD0ChLowPt], keepEvent[kJetD0ChHighPt], keepEvent[kJetLcChLowPt], keepEvent[kJetLcChHighPt]);
  }
  PROCESS_SWITCH(JetHFFilterTask, processJets, "Do HF charged jet triggering", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{adaptAnalysisTask<JetHFFilterTask>(cfg)};
}
