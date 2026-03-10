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
//
// \file JetD0AngSubstructure.cxx
//
// \brief Analysis task for the reconstruction and study of charged jets
//        containing D_0 mesons in pp collisions.
// \inherited from D0 fragmentation and Ds
// \P. Dhankher

#include "PWGHF/Core/DecayChannels.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetSubtraction.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include <CommonConstants/MathConstants.h>
#include <Framework/AnalysisTask.h>
#include <Framework/ConfigContext.h>
#include <Framework/Configurable.h>
#include <Framework/DataProcessorSpec.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include "TVector3.h"

#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// Definition of a custom AOD table to store jet–D0 quantities
namespace o2::aod
{
namespace jet_obj
{
// Jet-related quantities
DECLARE_SOA_COLUMN(JetHfDist, jetHfDist, float);
DECLARE_SOA_COLUMN(JetPt, jetPt, float);
DECLARE_SOA_COLUMN(JetEta, jetEta, float);
DECLARE_SOA_COLUMN(JetPhi, jetPhi, float);
DECLARE_SOA_COLUMN(JetNConst, jetNConst, int);
DECLARE_SOA_COLUMN(JetAng, jetAng, float);
// D0 candidate quantities
DECLARE_SOA_COLUMN(HfPt, hfPt, float);
DECLARE_SOA_COLUMN(HfEta, hfEta, float);
DECLARE_SOA_COLUMN(HfPhi, hfPhi, float);
DECLARE_SOA_COLUMN(HfMass, hfMass, float);
DECLARE_SOA_COLUMN(HfY, hfY, float);
// ML scores
DECLARE_SOA_COLUMN(HfMlScore0, hfMlScore0, float);
DECLARE_SOA_COLUMN(HfMlScore1, hfMlScore1, float);
DECLARE_SOA_COLUMN(HfMlScore2, hfMlScore2, float);
} // namespace jet_obj
// AOD table definition
DECLARE_SOA_TABLE(JetObjTable, "AOD", "JETOBJTABLE",
                  jet_obj::JetHfDist,
                  jet_obj::JetPt,
                  jet_obj::JetEta,
                  jet_obj::JetPhi,
                  jet_obj::JetNConst,
                  jet_obj::JetAng,
                  jet_obj::HfPt,
                  jet_obj::HfEta,
                  jet_obj::HfPhi,
                  jet_obj::HfMass,
                  jet_obj::HfY,
                  jet_obj::HfMlScore0,
                  jet_obj::HfMlScore1,
                  jet_obj::HfMlScore2);
} // namespace o2::aod
struct JetD0AngSubstructure {
  /**
   * Histogram registry
   *
   * Contains:
   *  - Event and track histograms
   *  - Jet kinematic distributions
   *  - D0–jet substructure observables
   */
  HistogramRegistry registry{"registry",
                             {{"h_collision_counter", "# of collisions;", {HistType::kTH1F, {{2, 0., 2.}}}},
                              {"h_jet_counter", ";# of D^{0} jets;", {HistType::kTH1F, {{6, 0., 3.0}}}},
                              {"h_d0_jet_projection", ";z^{D^{0},jet}_{||};dN/dz^{D^{0},jet}_{||}", {HistType::kTH1F, {{1000, 0., 10.}}}},
                              {"h_d0_jet_distance_vs_projection", ";#DeltaR_{D^{0},jet};z^{D^{0},jet}_{||}", {HistType::kTH2F, {{1000, 0., 10.}, {1000, 0., 10.}}}},
                              {"h_d0_jet_distance", ";#DeltaR_{D^{0},jet};dN/d(#DeltaR)", {HistType::kTH1F, {{1000, 0., 10.}}}},
                              {"h_d0_jet_pt", ";p_{T,D^{0} jet};dN/dp_{T,D^{0} jet}", {HistType::kTH1F, {{200, 0., 10.}}}},
                              {"h_d0_jet_eta", ";#eta_{T,D^{0} jet};dN/d#eta_{D^{0} jet}", {HistType::kTH1F, {{250, -5., 5.}}}},
                              {"h_d0_jet_phi", ";#phi_{T,D^{0} jet};dN/d#phi_{D^{0} jet}", {HistType::kTH1F, {{250, -10., 10.}}}},
                              {"h_d0_mass", ";m_{D^{0}});dN/dm_{D^{0}}", {HistType::kTH1F, {{1000, 0., 10.}}}},
                              {"h_d0_eta", ";#eta_{D^{0}});dN/d#eta_{D_{}}", {HistType::kTH1F, {{250, -5., 5.}}}},
                              {"h_d0_phi", ";#phi_{D^{0}});dN/d#phi_{D^{0}}", {HistType::kTH1F, {{250, -10., 10.}}}},
                              {"h_d0_ang", ";#lambda_{#kappa}^{#alpha};counts", {HistType::kTH1F, {{100, 0., 1.}}}}}};

  // Configurables
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};

  Configurable<float> jetPtMin{"jetPtMin", 5.0, "minimum jet pT cut"};
  Configurable<float> jetR{"jetR", 0.4, "jet resolution parameter"};

  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
  Configurable<float> kappa{"kappa", 1.0, "angularity kappa"}; // to do: configurable from json
  Configurable<float> alpha{"alpha", 1.0, "angularity alpha"};

  std::vector<int> eventSelectionBits;
  int trackSelection = -1;

  // Output table producer
  Produces<aod::JetObjTable> ObjJetTable;

  float angularity;

  void init(o2::framework::InitContext&)
  {
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
  }

  template <typename T, typename U>
  void jetCalculateAngularity(T const& jet, U const& /*tracks*/)
  {
    angularity = 0.0;
    for (auto& constituent : jet.template tracks_as<U>()) {
      angularity += std::pow(constituent.pt(), kappa) * std::pow(jetutilities::deltaR(jet, constituent) / (jet.r() / 100.f), alpha);
    }
    angularity /= std::pow(jet.pt(), kappa);
  }

  void processDataChargedSubstructure(aod::JetCollision const& collision,
                                      soa::Join<aod::D0ChargedJets, aod::D0ChargedJetConstituents> const& jets,
                                      aod::CandidatesD0Data const&, aod::JetTracks const& tracks)
  {
    // apply event selection and fill histograms for sanity check
    registry.fill(HIST("h_collision_counter"), 0.5);

    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits) || !(std::abs(collision.posZ()) < vertexZCut)) {
      return;
    }
    registry.fill(HIST("h_collision_counter"), 1.5);

    // Loop over jets containing D0 candidates
    for (const auto& jet : jets) {
      // number of charged jets with D0
      registry.fill(HIST("h_jet_counter"), 0.5);
      // obtaining jet 3-vector
      TVector3 jetVector(jet.px(), jet.py(), jet.pz());

      // Loop over D0 candidates associated to the jet
      for (const auto& d0Candidate : jet.candidates_as<aod::CandidatesD0Data>()) {
        // obtaining jet 3-vector
        TVector3 d0Vector(d0Candidate.px(), d0Candidate.py(), d0Candidate.pz());

        // calculating fraction of the jet momentum carried by the D0 along the direction of the jet axis
        double zParallel = (jetVector * d0Vector) / (jetVector * jetVector);

        // calculating angular distance in eta-phi plane
        double axisDistance = jetutilities::deltaR(jet, d0Candidate);

        jetCalculateAngularity(jet, tracks);

        // filling histograms
        registry.fill(HIST("h_d0_jet_projection"), zParallel);
        registry.fill(HIST("h_d0_jet_distance_vs_projection"), axisDistance, zParallel);
        registry.fill(HIST("h_d0_jet_distance"), axisDistance);
        registry.fill(HIST("h_d0_jet_pt"), jet.pt());
        registry.fill(HIST("h_d0_jet_eta"), jet.eta());
        registry.fill(HIST("h_d0_jet_phi"), jet.phi());
        registry.fill(HIST("h_d0_mass"), d0Candidate.m());
        registry.fill(HIST("h_d0_eta"), d0Candidate.eta());
        registry.fill(HIST("h_d0_phi"), d0Candidate.phi());
        registry.fill(HIST("h_d0_ang"), angularity); // add more axis

        // filling table
        ObjJetTable(axisDistance,
                    jet.pt(), jet.eta(), jet.phi(), jet.tracks_as<aod::JetTracks>().size(), angularity,
                    d0Candidate.pt(), d0Candidate.eta(), d0Candidate.phi(), d0Candidate.m(), d0Candidate.y(), d0Candidate.mlScores()[0], d0Candidate.mlScores()[1], d0Candidate.mlScores()[2]);

        break; // get out of candidates' loop after first HF particle is found in jet
      } // end of D0 candidates loop

    } // end of jets loop

  } // end of process function
  PROCESS_SWITCH(JetD0AngSubstructure, processDataChargedSubstructure, "charged HF jet substructure", false);
};
// Workflow definition
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<JetD0AngSubstructure>(cfgc, TaskName{"jet-d0-ang-substructure"})}; }
