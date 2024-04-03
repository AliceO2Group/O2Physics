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

static const std::vector<std::string> highPtObjectsNames = {"JetD0ChLowPt", "JetD0ChHighPt", "JetLcChLowPt", "JetLcChHighPt"};

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
  Configurable<bool> fillTHns{"fillTHns", true, "fill THn histograms"};

  Configurable<std::vector<double>> jetRadiiPlot{"jetRadiiPlot", std::vector<double>{0.2, 0.4, 0.6}, "jet resolution parameters"};

  void init(o2::framework::InitContext&)
  {
    auto jetRadiiPlotBins = (std::vector<double>)jetRadiiPlot;
    if (jetRadiiPlotBins.size() > 1) {
      jetRadiiPlotBins.push_back(jetRadiiPlotBins[jetRadiiPlotBins.size() - 1] + (TMath::Abs(jetRadiiPlotBins[jetRadiiPlotBins.size() - 1] - jetRadiiPlotBins[jetRadiiPlotBins.size() - 2])));
    } else {
      jetRadiiPlotBins.push_back(jetRadiiPlotBins[jetRadiiPlotBins.size() - 1] + 0.1);
    }

    registry.add("h_d0jet_pt", "D^{0} - tagged jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}});
    registry.add("h_d0jet_pt_lowpt", "D^{0} - tagged jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}});
    registry.add("h_d0jet_pt_highpt", "D^{0} - tagged jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}});

    registry.add("h_lcjet_pt", "#Lambda^{+}_{c} - tagged jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}});
    registry.add("h_lcjet_pt_lowpt", "#Lambda^{+}_{c} - tagged jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}});
    registry.add("h_lcjet_pt_highpt", "#Lambda^{+}_{c} - tagged jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}});

    // these might end up being too big
    registry.add("d0Thn", "Thn for D^{0}-tagged jets", {HistType::kTHnC, {{jetRadiiPlotBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}, {160, -1.0, 7.}, {200, 0., 200.}, {200, 0., 200.}, {100, -1.0, 1.0}, {160, -1.0, 7.}, {100, -1.0, 1.0}, {1700, 1.3, 3.0}}});
    registry.add("d0Thn_witheventcuts", "Thn for D^{0}-tagged jets with event cuts", {HistType::kTHnC, {{jetRadiiPlotBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}, {160, -1.0, 7.}, {200, 0., 200.}, {200, 0., 200.}, {100, -1.0, 1.0}, {160, -1.0, 7.}, {100, -1.0, 1.0}, {1700, 1.3, 3.0}}});
    registry.add("lcThn", "Thn for #Lambda^{+}_{c}-tagged jets", {HistType::kTHnC, {{jetRadiiPlotBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}, {160, -1.0, 7.}, {200, 0., 200.}, {200, 0., 200.}, {100, -1.0, 1.0}, {160, -1.0, 7.}, {100, -1.0, 1.0}, {1700, 1.3, 3.0}}});
    registry.add("lcThn_witheventcuts", "Thn for #Lambda^{+}_{c}-tagged jets with event cuts", {HistType::kTHnC, {{jetRadiiPlotBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}, {160, -1.0, 7.}, {200, 0., 200.}, {200, 0., 200.}, {100, -1.0, 1.0}, {160, -1.0, 7.}, {100, -1.0, 1.0}, {1700, 1.3, 3.0}}});
    // radius, JetPt, JetEta, Jet Phi, Jet Ntracks, HF pT, HF Eta, HF Phi, HF Y, HF Mass
  }

  void processJets(soa::Join<JetCollisions, aod::EvSels>::iterator const& collision, soa::Join<o2::aod::D0ChargedJets, o2::aod::D0ChargedJetConstituents> const& d0Jets, CandidatesD0Data const& d0Candidates, soa::Join<o2::aod::LcChargedJets, o2::aod::LcChargedJetConstituents> const& lcJets, CandidatesLcData const& lcCandidates, JetTracks const& tracks)
  {
    bool keepEvent[kAllObjects]{false};
    for (auto const& d0Jet : d0Jets) {
      if (fillTHns) {
        for (auto const& d0Candidate : d0Jet.hfcandidates_as<CandidatesD0Data>()) {
          registry.fill(HIST("d0Thn"), d0Jet.r() / 100.0, d0Jet.pt(), d0Jet.eta(), d0Jet.phi(), d0Jet.tracksIds().size() + d0Jet.hfcandidatesIds().size(), d0Candidate.pt(), d0Candidate.eta(), d0Candidate.phi(), d0Candidate.y(), d0Candidate.m());
          if (collision.posZ() < vertexZCut && collision.sel8() && collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
            registry.fill(HIST("d0Thn_witheventcuts"), d0Jet.r() / 100.0, d0Jet.pt(), d0Jet.eta(), d0Jet.phi(), d0Jet.tracksIds().size() + d0Jet.hfcandidatesIds().size(), d0Candidate.pt(), d0Candidate.eta(), d0Candidate.phi(), d0Candidate.y(), d0Candidate.m());
          }
          break;
        }
      }
      if (d0Jet.r() == round(jetD0ChR * 100.0f)) {
        registry.fill(HIST("h_d0jet_pt"), d0Jet.pt());
        if (d0Jet.pt() >= jetD0ChLowPtThreshold) {
          keepEvent[kJetD0ChLowPt] = true;
          registry.fill(HIST("h_d0jet_pt_lowpt"), d0Jet.pt());
        }
        if (d0Jet.pt() >= jetD0ChHighPtThreshold) {
          keepEvent[kJetD0ChHighPt] = true;
          registry.fill(HIST("h_d0jet_pt_highpt"), d0Jet.pt());
        }
      }
    }
    for (auto const& lcJet : lcJets) {
      if (fillTHns) {
        for (auto const& lcCandidate : lcJet.hfcandidates_as<CandidatesLcData>()) {
          registry.fill(HIST("lcThn"), lcJet.r() / 100.0, lcJet.pt(), lcJet.eta(), lcJet.phi(), lcJet.tracksIds().size() + lcJet.hfcandidatesIds().size(), lcCandidate.pt(), lcCandidate.eta(), lcCandidate.phi(), lcCandidate.y(), lcCandidate.m());
          if (collision.posZ() < vertexZCut && collision.sel8() && collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
            registry.fill(HIST("lcThn_witheventcuts"), lcJet.r() / 100.0, lcJet.pt(), lcJet.eta(), lcJet.phi(), lcJet.tracksIds().size() + lcJet.hfcandidatesIds().size(), lcCandidate.pt(), lcCandidate.eta(), lcCandidate.phi(), lcCandidate.y(), lcCandidate.m());
          }
          break;
        }
      }
      if (lcJet.r() == round(jetLcChR * 100.0f)) {
        registry.fill(HIST("h_lcjet_pt"), lcJet.pt());
        if (lcJet.pt() >= jetLcChLowPtThreshold) {
          keepEvent[kJetLcChLowPt] = true;
          registry.fill(HIST("h_lcjet_pt_lowpt"), lcJet.pt());
        }
        if (lcJet.pt() >= jetLcChHighPtThreshold) {
          keepEvent[kJetLcChHighPt] = true;
          registry.fill(HIST("h_lcjet_pt_highpt"), lcJet.pt());
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
