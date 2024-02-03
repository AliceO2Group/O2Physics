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

// h+jet analysis task
//
// Authors: Daniel Jones

#include <cmath>
#include <vector>

#include "TRandom3.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/HistogramRegistry.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/runDataProcessing.h"
#include "EventFiltering/filterTables.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/DataModel/Jet.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct hJetAnalysis {

  Configurable<float> pt_TTref_min{"pt_TTref_min", 5, "reference minimum trigger track pt"};
  Configurable<float> pt_TTref_max{"pt_TTref_max", 7, "reference maximum trigger track pt"};
  Configurable<float> pt_TTsig_min{"pt_TTsig_min", 20, "signal minimum trigger track pt"};
  Configurable<float> pt_TTsig_max{"pt_TTsig_max", 50, "signal maximum trigger track pt"};
  Configurable<float> frac_sig{"frac_sig", 0.5, "fraction of events to use for signal"};
  Configurable<float> jetR{"jetR", 0.4, "jet resolution parameter"};

  TRandom3* rand = new TRandom3(0);

  Filter jetCuts = aod::jet::r == nround(jetR.node() * 100.0f);

  HistogramRegistry registry{"registry",
                             {{"hNtrig", "number of triggers;trigger type;entries", {HistType::kTH1F, {{2, 0, 2}}}},
                              {"hPtTrack", "Track p_{T};p_{T};entries", {HistType::kTH1F, {{100, 0, 100}}}},
                              {"hEtaTrack", "Track #eta;#eta;entries", {HistType::kTH1F, {{20, -1, 1}}}},
                              {"hPhiTrack", "Track #phi;#phi;entries", {HistType::kTH1F, {{200, -M_PI, 2 * M_PI}}}},
                              {"hReferencePtDPhi", "jet p_{T} vs DPhi;p_{T,jet};#Delta#phi", {HistType::kTH2F, {{150, 0, 150}, {100, -M_PI, M_PI}}}},
                              {"hSignalPtDPhi", "jet p_{T} vs DPhi;p_{T,jet};#Delta#phi", {HistType::kTH2F, {{150, 0, 150}, {100, -M_PI, M_PI}}}},
                              {"hReferencePt", "jet p_{T};p_{T,jet};entries", {HistType::kTH1F, {{150, 0, 150}}}},
                              {"hSignalPt", "jet p_{T};p_{T,jet};entries", {HistType::kTH1F, {{150, 0, 150}}}},
                              {"hSignalLeadingTrack", "leading track p_{T};p_{T,jet};#Delta#phi;leading track p_{T}", {HistType::kTH3F, {{150, 0, 150}, {100, -M_PI, M_PI}, {150, 0, 150}}}},
                              {"hReferenceLeadingTrack", "leading track p_{T};p_{T,jet};#Delta#phi;leading track p_{T}", {HistType::kTH3F, {{150, 0, 150}, {100, -M_PI, M_PI}, {150, 0, 150}}}},
                              {"hJetSignalMultiplicity", "jet multiplicity;N_{jets};entries", {HistType::kTH1F, {{10, 0, 10}}}},
                              {"hJetReferenceMultiplicity", "jet multiplicity;N_{jets};entries", {HistType::kTH1F, {{10, 0, 10}}}},
                              {"hJetSignalConstituentMultiplicity", "jet constituent multiplicity;p_{T,jet};#Delta#phi;N_{constituents}", {HistType::kTH3F, {{150, 0, 150}, {100, -M_PI, M_PI}, {50, 0, 50}}}},
                              {"hJetReferenceConstituentMultiplicity", "jet constituent multiplicity;p_{T,jet};#Delta#phi;N_{constituents}", {HistType::kTH3F, {{150, 0, 150}, {100, -M_PI, M_PI}, {50, 0, 50}}}},
                              {"hJetSignalConstituentPt", "jet constituent p_{T};p_{T,jet};#Delta#phi;p_{T,constituent}", {HistType::kTH3F, {{150, 0, 150}, {100, -M_PI, M_PI}, {150, 0, 150}}}},
                              {"hJetReferenceConstituentPt", "jet constituent p_{T};p_{T,jet};#Delta#phi;p_{T,constituent}", {HistType::kTH3F, {{150, 0, 150}, {100, -M_PI, M_PI}, {150, 0, 150}}}},
                              {"hSigEventTriggers", "N_{triggers};events", {HistType::kTH1F, {{10, 0, 10}}}},
                              {"hRefEventTriggers", "N_{triggers};events", {HistType::kTH1F, {{10, 0, 10}}}},
                              {"hJetPt", "jet p_{T};p_{T,jet};entries", {HistType::kTH1F, {{150, 0, 150}}}},
                              {"hJetEta", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{20, -1, 1}}}},
                              {"hJetPhi", "jet #phi;#phi_{jet};entries", {HistType::kTH1F, {{200, -M_PI, 2 * M_PI}}}}}};

  void init(InitContext const&) {}

  float dPhi(float phi1, float phi2)
  {
    float dPhi = phi1 - phi2;
    if (dPhi < -M_PI)
      dPhi += 2 * M_PI;
    if (dPhi > M_PI)
      dPhi -= 2 * M_PI;
    return dPhi;
  }

  template <typename T, typename U>
  void fillHistograms(T const& jets, U const& tracks)
  {
    bool is_sig_col;
    std::vector<double> phi_TT_ar;
    double phi_TT;
    int trig_number;
    int n_TT = 0;
    double leadingPT = 0;

    float dice = rand->Rndm();
    if (dice < frac_sig)
      is_sig_col = true;
    else
      is_sig_col = false;

    for (auto& track : tracks) {
      if (is_sig_col && track.pt() < pt_TTsig_max && track.pt() > pt_TTsig_min) {
        phi_TT_ar.push_back(track.phi());
        n_TT++;
      }
      if (!is_sig_col && track.pt() < pt_TTref_max && track.pt() > pt_TTref_min) {
        phi_TT_ar.push_back(track.pt());
        n_TT++;
      }
      registry.fill(HIST("hPtTrack"), track.pt());
      registry.fill(HIST("hEtaTrack"), track.eta());
      registry.fill(HIST("hPhiTrack"), track.phi());
    }

    if (n_TT > 0) {
      trig_number = rand->Integer(n_TT);
      phi_TT = phi_TT_ar[trig_number];
      if (is_sig_col) {
        registry.fill(HIST("hNtrig"), 1.5);
        registry.fill(HIST("hJetSignalMultiplicity"), jets.size());
        registry.fill(HIST("hSigEventTriggers"), n_TT);
      }
      if (!is_sig_col) {
        registry.fill(HIST("hNtrig"), 0.5);
        registry.fill(HIST("hJetReferenceMultiplicity"), jets.size());
        registry.fill(HIST("hRefEventTriggers"), n_TT);
      }
    }

    for (auto& jet : jets) {
      registry.fill(HIST("hJetPt"), jet.pt());
      registry.fill(HIST("hJetEta"), jet.eta());
      registry.fill(HIST("hJetPhi"), jet.phi());
      if (n_TT > 0) {
        float dphi = dPhi(jet.phi(), phi_TT);
        if (is_sig_col) {
          registry.fill(HIST("hSignalPtDPhi"), jet.pt(), dphi);
          registry.fill(HIST("hSignalPt"), jet.pt());
          registry.fill(HIST("hJetSignalConstituentMultiplicity"), jet.pt(), dphi, jet.tracks().size());
          for (auto& constituent : jet.template tracks_as<U>()) {
            if (constituent.pt() > leadingPT) {
              leadingPT = constituent.pt();
            }
            registry.fill(HIST("hJetSignalConstituentPt"), jet.pt(), dphi, constituent.pt());
          }
          registry.fill(HIST("hSignalLeadingTrack"), jet.pt(), dphi, leadingPT);
        }
        if (!is_sig_col) {
          registry.fill(HIST("hReferencePtDPhi"), jet.pt(), dphi);
          registry.fill(HIST("hReferencePt"), jet.pt());
          registry.fill(HIST("hJetReferenceConstituentMultiplicity"), jet.pt(), dphi, jet.tracks().size());
          for (auto& constituent : jet.template tracks_as<U>()) {
            if (constituent.pt() > leadingPT) {
              leadingPT = constituent.pt();
            }
            registry.fill(HIST("hJetReferenceConstituentPt"), jet.pt(), dphi, constituent.pt());
          }
          registry.fill(HIST("hReferenceLeadingTrack"), jet.pt(), dphi, leadingPT);
        }
      }
    }
  }

  void processData(aod::JCollision const&, soa::Filtered<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>> const& jets, aod::JTracks const& tracks)
  {
    fillHistograms(jets, tracks);
  }
  PROCESS_SWITCH(hJetAnalysis, processData, "process data", true);

  void processMCD(aod::JCollision const&, soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>> const& jets, aod::JTracks const& tracks)
  {
    fillHistograms(jets, tracks);
  }
  PROCESS_SWITCH(hJetAnalysis, processMCD, "process MC detector level", false);

  void processMCP(aod::JMcCollision const&, soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents>> const& jets, aod::JMcParticles const& particles)
  {
    fillHistograms(jets, particles);
  }
  PROCESS_SWITCH(hJetAnalysis, processMCP, "process MC particle level", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<hJetAnalysis>(cfgc, TaskName{"hJetAnalysis"})}; }
