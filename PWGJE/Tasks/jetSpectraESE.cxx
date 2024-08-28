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

// jet analysis framework with ESE (19/08/2024)
//
/// \author Joachim Hansen <joachim.hansen@cern.ch>
//

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"

#include "PWGCF/Flow/DataModel/FlowESETable.h"
#include "PWGCF/Flow/Core/FFitWeights.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"

using ColWqVecFT0C = soa::Join<aod::Collisions, aod::CentFT0Cs, aod::qVecFT0Cs, aod::qPercentileFT0Cs>;
using JColwESE = soa::Join<aod::JCollisions, aod::CentFT0Cs, aod::qVecFT0Cs, aod::qPercentileFT0Cs>;

struct JetSpectraEseTask {

  HistogramRegistry registry{"registry",
                             {{"h_collisions", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}}},
                              {"h_jet_pt", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_jet_eta", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
                              {"h_jet_phi", "jet #phi;#phi_{jet};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
                              {"h_Cent_q2Percs", "qPercs;Centrality;", {HistType::kTH2F, {{100, 0, 100}, {250, 0, 5000}}}},
                              {"h_jet_pt_q2_0", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_jet_pt_q2_1", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_jet_pt_in", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_jet_pt_out", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}}

                             }};

  Configurable<float> jetPtMin{"jetPtMin", 5.0, "minimum jet pT cut"};
  Configurable<float> jetR{"jetR", 0.4, "jet resolution parameter"};
  Configurable<std::vector<float>> qLimits{"qLimits", {30, 70}, "jet resolution parameter"};

  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

  int eventSelection = -1;
  int trackSelection = -1;

  void init(o2::framework::InitContext&)
  {
    eventSelection = jetderiveddatautilities::initialiseEventSelection(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
  }

  Filter jetCuts = aod::jet::pt > jetPtMin&& aod::jet::r == nround(jetR.node() * 100.0f);

  int spCode(const int qPerc)
  {
    if (qPerc < 0)
      return -1;

    if (qPerc <= qLimits->at(0)) {
      return 0;
    } else if (qPerc >= qLimits->at(1)) {
      return 1;
    } else {
      return -1;
    }
  }
  std::string jetplane(float deltaPhi)
  {
    std::string s = "";
    if (deltaPhi < TMath::Pi() / 6.) {
      s = "in";
    } else if (deltaPhi > TMath::Pi() / 3.) {
      s = "out";
    }

    return s;
  }

  void processESEDataCharged(JColwESE::iterator const& collision, soa::Filtered<aod::ChargedJets> const& jets, aod::BCsWithTimestamps const&)
  {

    float vPsi2 = FFitWeights::EventPlane(collision, 2);

    auto qPerc = collision.qPERCFT0C();
    if (qPerc[0] > 0)
      registry.fill(HIST("h_Cent_q2Percs"), collision.centFT0C(), qPerc[0]);

    int code = spCode(qPerc[0]);

    registry.fill(HIST("h_collisions"), 0.5);
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }

    registry.fill(HIST("h_collisions"), 1.5);
    float plane{0};
    std::string pl{""};
    for (auto const& jet : jets) {

      plane = vPsi2 - jet.phi();
      pl = jetplane(plane);
      registry.fill(HIST("h_jet_pt"), jet.pt());

      if (code < 0)
        continue;

      if (code == 0)
        registry.fill(HIST("h_jet_pt_q2_0"), jet.pt());
      if (code == 1)
        registry.fill(HIST("h_jet_pt_q2_1"), jet.pt());
      if (pl == "in")
        registry.fill(HIST("h_jet_pt_in"), jet.pt());
      if (pl == "out")
        registry.fill(HIST("h_jet_pt_out"), jet.pt());
    }
  }
  PROCESS_SWITCH(JetSpectraEseTask, processESEDataCharged, "process self contained collisions", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<JetSpectraEseTask>(cfgc, TaskName{"jet-spectra-ese"})}; }
