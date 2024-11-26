// Copyright 2019-2024 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \brief charm hadron hadronization task
/// \author Christian Reckziegel <christian.reckziegel@cern.ch>, Federal University of ABC
/// \since 15.03.2024
///
/// The task store data relevant to the calculation of hadronization observables radial
/// profile and/or jet momentum fraction for charmed hadrons
#include <vector>
#include <string>

#include "TVector3.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/runDataProcessing.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/HistogramRegistry.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetSubstructure.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetHFUtilities.h"
#include "PWGJE/Core/JetUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// calculate delta phi such that 0 < delta phi < pi
double deltaPhi(double phi1, double phi2)
{
  // Compute the absolute difference between phi1 and phi2
  double dphi = std::abs(phi1 - phi2);

  // Constrain angle between [min,min+2pi] = [-pi,-pi+2pi]
  dphi = RecoDecay::constrainAngle(dphi, -o2::constants::math::PI);

  // Return absolute value of distance
  return std::abs(dphi);
}

// creating table for storing distance data
namespace o2::aod
{
namespace jet_distance
{
DECLARE_SOA_COLUMN(JetHfDist, jethfdist, float);
DECLARE_SOA_COLUMN(JetPt, jetpt, float);
DECLARE_SOA_COLUMN(JetEta, jeteta, float);
DECLARE_SOA_COLUMN(JetPhi, jetphi, float);
DECLARE_SOA_COLUMN(JetNConst, jetnconst, float);
DECLARE_SOA_COLUMN(HfPt, hfpt, float);
DECLARE_SOA_COLUMN(HfEta, hfeta, float);
DECLARE_SOA_COLUMN(HfPhi, hfphi, float);
DECLARE_SOA_COLUMN(HfMass, hfmass, float);
DECLARE_SOA_COLUMN(HfY, hfy, float);
DECLARE_SOA_COLUMN(HfPrompt, hfPrompt, bool);
DECLARE_SOA_COLUMN(HfMatch, hfmatch, bool);
DECLARE_SOA_COLUMN(HfMlScore0, hfmlscore0, float);
DECLARE_SOA_COLUMN(HfMlScore1, hfmlscore1, float);
DECLARE_SOA_COLUMN(HfMlScore2, hfmlscore2, float);
DECLARE_SOA_COLUMN(HfMatchedFrom, hfmatchedfrom, int);
DECLARE_SOA_COLUMN(HfSelectedAs, hfselectedas, int);
DECLARE_SOA_COLUMN(MCJetHfDist, mcjethfdist, float);
DECLARE_SOA_COLUMN(MCJetPt, mcjetpt, float);
DECLARE_SOA_COLUMN(MCJetEta, mcjeteta, float);
DECLARE_SOA_COLUMN(MCJetPhi, mcjetphi, float);
DECLARE_SOA_COLUMN(MCJetNConst, mcjetnconst, float);
DECLARE_SOA_COLUMN(MCHfPt, mchfpt, float);
DECLARE_SOA_COLUMN(MCHfEta, mchfeta, float);
DECLARE_SOA_COLUMN(MCHfPhi, mchfphi, float);
DECLARE_SOA_COLUMN(MCHfY, mchfy, float);
DECLARE_SOA_COLUMN(MCHfPrompt, mchfPrompt, bool);
DECLARE_SOA_COLUMN(MCHfMatch, mchfmatch, bool);
} // namespace jet_distance
DECLARE_SOA_TABLE(JetDistanceTable, "AOD", "JETDISTTABLE",
                  jet_distance::JetHfDist,
                  jet_distance::JetPt,
                  jet_distance::JetEta,
                  jet_distance::JetPhi,
                  jet_distance::JetNConst,
                  jet_distance::HfPt,
                  jet_distance::HfEta,
                  jet_distance::HfPhi,
                  jet_distance::HfMass,
                  jet_distance::HfY,
                  jet_distance::HfMlScore0,
                  jet_distance::HfMlScore1,
                  jet_distance::HfMlScore2);
DECLARE_SOA_TABLE(MCPJetDistanceTable, "AOD", "MCPJETDISTTABLE",
                  jet_distance::MCJetHfDist,
                  jet_distance::MCJetPt,
                  jet_distance::MCJetEta,
                  jet_distance::MCJetPhi,
                  jet_distance::MCJetNConst,
                  jet_distance::MCHfPt,
                  jet_distance::MCHfEta,
                  jet_distance::MCHfPhi,
                  jet_distance::MCHfY,
                  jet_distance::MCHfPrompt,
                  jet_distance::MCHfMatch);
DECLARE_SOA_TABLE(MCDJetDistanceTable, "AOD", "MCDJETDISTTABLE",
                  jet_distance::JetHfDist,
                  jet_distance::JetPt,
                  jet_distance::JetEta,
                  jet_distance::JetPhi,
                  jet_distance::JetNConst,
                  jet_distance::HfPt,
                  jet_distance::HfEta,
                  jet_distance::HfPhi,
                  jet_distance::HfMass,
                  jet_distance::HfY,
                  jet_distance::HfPrompt,
                  jet_distance::HfMatch,
                  jet_distance::HfMlScore0,
                  jet_distance::HfMlScore1,
                  jet_distance::HfMlScore2,
                  jet_distance::HfMatchedFrom,
                  jet_distance::HfSelectedAs);
DECLARE_SOA_TABLE(MatchJetDistanceTable, "AOD", "MATCHTABLE",
                  jet_distance::MCJetHfDist,
                  jet_distance::MCJetPt,
                  jet_distance::MCJetEta,
                  jet_distance::MCJetPhi,
                  jet_distance::MCJetNConst,
                  jet_distance::MCHfPt,
                  jet_distance::MCHfEta,
                  jet_distance::MCHfPhi,
                  jet_distance::MCHfY,
                  jet_distance::MCHfPrompt,
                  jet_distance::JetHfDist,
                  jet_distance::JetPt,
                  jet_distance::JetEta,
                  jet_distance::JetPhi,
                  jet_distance::JetNConst,
                  jet_distance::HfPt,
                  jet_distance::HfEta,
                  jet_distance::HfPhi,
                  jet_distance::HfMass,
                  jet_distance::HfY,
                  jet_distance::HfPrompt,
                  jet_distance::HfMlScore0,
                  jet_distance::HfMlScore1,
                  jet_distance::HfMlScore2,
                  jet_distance::HfMatchedFrom,
                  jet_distance::HfSelectedAs);
} // namespace o2::aod

struct HfFragmentationFunctionTask {
  // producing new table
  Produces<aod::JetDistanceTable> distJetTable;
  Produces<aod::MCPJetDistanceTable> mcpdistJetTable;
  Produces<aod::MCDJetDistanceTable> mcddistJetTable;
  Produces<aod::MatchJetDistanceTable> matchJetTable;

  // Tables for MC jet matching
  using JetMCDTable = soa::Join<aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents, aod::D0ChargedMCDetectorLevelJetsMatchedToD0ChargedMCParticleLevelJets>;
  using JetMCPTable = soa::Join<aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents, aod::D0ChargedMCParticleLevelJetsMatchedToD0ChargedMCDetectorLevelJets>;

  // slices for accessing proper HF mcdjets collision associated to mccollisions
  PresliceUnsorted<aod::JetCollisionsMCD> CollisionsPerMCCollision = aod::jmccollisionlb::mcCollisionId;
  Preslice<JetMCDTable> D0MCDJetsPerCollision = aod::jet::collisionId;
  Preslice<JetMCPTable> D0MCPJetsPerMCCollision = aod::jet::mcCollisionId;

  // Histogram registry: an object to hold your histograms
  HistogramRegistry registry{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};

  int eventSelection = -1;

  void init(InitContext const&)
  {
    // initialise event selection:
    eventSelection = jetderiveddatautilities::initialiseEventSelection(static_cast<std::string>(eventSelections));

    // create histograms
    // collision system histograms
    std::vector<std::string> histLabels = {"mccollisions", "z_cut", "collisions", "sel8"};
    registry.add("h_collision_counter", ";# of collisions;", HistType::kTH1F, {{static_cast<int>(histLabels.size()), 0.0, static_cast<double>(histLabels.size())}});
    auto counter = registry.get<TH1>(HIST("h_collision_counter"));
    for (std::vector<std::string>::size_type iCounter = 0; iCounter < histLabels.size(); iCounter++) {
      counter->GetXaxis()->SetBinLabel(iCounter + 1, histLabels[iCounter].data());
    }
    registry.add("h_jet_counter", ";# of jets;", {HistType::kTH1F, {{6, 0., 3.0}}});
    auto jetCounter = registry.get<TH1>(HIST("h_jet_counter"));
    jetCounter->GetXaxis()->SetBinLabel(1, "particle level");
    jetCounter->GetXaxis()->SetBinLabel(2, "detector level");
    jetCounter->GetXaxis()->SetBinLabel(3, "particle matched jets");
    jetCounter->GetXaxis()->SetBinLabel(4, "detector matched jets");
    jetCounter->GetXaxis()->SetBinLabel(5, "mcd matched to mcp loop");
    jetCounter->GetXaxis()->SetBinLabel(6, "mcp matched to mcd loop");
    // D0 candidate histograms from data
    registry.add("h_d0_jet_projection", ";z^{D^{0},jet}_{||};dN/dz^{D^{0},jet}_{||}", {HistType::kTH1F, {{1000, 0., 10.}}});
    registry.add("h_d0_jet_distance_vs_projection", ";#DeltaR_{D^{0},jet};z^{D^{0},jet}_{||}", {HistType::kTH2F, {{1000, 0., 10.}, {1000, 0., 10.}}});
    registry.add("h_d0_jet_distance", ";#DeltaR_{D^{0},jet};dN/d(#DeltaR)", {HistType::kTH1F, {{1000, 0., 10.}}});
    registry.add("h_d0_jet_pt", ";p_{T,D^{0} jet};dN/dp_{T,D^{0} jet}", {HistType::kTH1F, {{200, 0., 10.}}});
    registry.add("h_d0_jet_eta", ";#eta_{T,D^{0} jet};dN/d#eta_{D^{0} jet}", {HistType::kTH1F, {{250, -5., 5.}}});
    registry.add("h_d0_jet_phi", ";#phi_{T,D^{0} jet};dN/d#phi_{D^{0} jet}", {HistType::kTH1F, {{250, -10., 10.}}});
    registry.add("h_d0_mass", ";m_{D^{0}} (GeV/c^{2});dN/dm_{D^{0}}", {HistType::kTH1F, {{1000, 0., 10.}}});
    registry.add("h_d0_eta", ";#eta_{D^{0}} (GeV/c^{2});dN/d#eta_{D^{0}}", {HistType::kTH1F, {{250, -5., 5.}}});
    registry.add("h_d0_phi", ";#phi_{D^{0}} (GeV/c^{2});dN/d#phi_{D^{0}}", {HistType::kTH1F, {{250, -10., 10.}}});
  }

  void processDummy(aod::TracksIU const&) {}
  PROCESS_SWITCH(HfFragmentationFunctionTask, processDummy, "Dummy process function turned on by default", true);

  void processDataChargedSubstructure(aod::JetCollision const& collision,
                                      soa::Join<aod::D0ChargedJets, aod::D0ChargedJetConstituents> const& jets,
                                      aod::CandidatesD0Data const&)
  {
    // apply event selection and fill histograms for sanity check
    registry.fill(HIST("h_collision_counter"), 2.0);
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection) || !(abs(collision.posZ()) < vertexZCut)) {
      return;
    }
    registry.fill(HIST("h_collision_counter"), 3.0);

    for (auto& jet : jets) {
      // fill jet counter histogram
      registry.fill(HIST("h_jet_counter"), 0.5);
      // obtaining jet 3-vector
      TVector3 jetVector(jet.px(), jet.py(), jet.pz());

      for (auto& d0Candidate : jet.candidates_as<aod::CandidatesD0Data>()) {

        // obtaining jet 3-vector
        TVector3 d0Vector(d0Candidate.px(), d0Candidate.py(), d0Candidate.pz());

        // calculating fraction of the jet momentum carried by the D0 along the direction of the jet axis
        double z_parallel = (jetVector * d0Vector) / (jetVector * jetVector);

        // calculating angular distance in eta-phi plane
        double axisDistance = jetutilities::deltaR(jet, d0Candidate);

        // filling histograms
        registry.fill(HIST("h_d0_jet_projection"), z_parallel);
        registry.fill(HIST("h_d0_jet_distance_vs_projection"), axisDistance, z_parallel);
        registry.fill(HIST("h_d0_jet_distance"), axisDistance);
        registry.fill(HIST("h_d0_jet_pt"), jet.pt());
        registry.fill(HIST("h_d0_jet_eta"), jet.eta());
        registry.fill(HIST("h_d0_jet_phi"), jet.phi());
        registry.fill(HIST("h_d0_mass"), d0Candidate.m());
        registry.fill(HIST("h_d0_eta"), d0Candidate.eta());
        registry.fill(HIST("h_d0_phi"), d0Candidate.phi());

        // filling table
        distJetTable(axisDistance,
                     jet.pt(), jet.eta(), jet.phi(), jet.tracks_as<aod::JetTracks>().size(),
                     d0Candidate.pt(), d0Candidate.eta(), d0Candidate.phi(), d0Candidate.m(), d0Candidate.y(), d0Candidate.mlScores()[0], d0Candidate.mlScores()[1], d0Candidate.mlScores()[2]);

        break; // get out of candidates' loop after first HF particle is found in jet
      } // end of D0 candidates loop

    } // end of jets loop

  } // end of process function
  PROCESS_SWITCH(HfFragmentationFunctionTask, processDataChargedSubstructure, "charged HF jet substructure", false);

  void processMcEfficiency(aod::JetMcCollisions const& mccollisions,
                           aod::JetCollisionsMCD const& collisions,
                           JetMCDTable const& mcdjets,
                           JetMCPTable const& mcpjets,
                           aod::CandidatesD0MCD const&,
                           aod::CandidatesD0MCP const&,
                           aod::JetTracks const&,
                           aod::JetParticles const&)
  {
    for (const auto& mccollision : mccollisions) {

      registry.fill(HIST("h_collision_counter"), 0.0);
      // skip collisions outside of |z| < vertexZCut
      if (abs(mccollision.posZ()) > vertexZCut) {
        continue;
      }
      registry.fill(HIST("h_collision_counter"), 1.0);

      // reconstructed collisions associated to same mccollision
      const auto collisionsPerMCCollision = collisions.sliceBy(CollisionsPerMCCollision, mccollision.globalIndex());
      for (const auto& collision : collisionsPerMCCollision) {

        registry.fill(HIST("h_collision_counter"), 2.0);
        if (!jetderiveddatautilities::selectCollision(collision, eventSelection) || !(abs(collision.posZ()) < vertexZCut)) {
          continue;
        }
        registry.fill(HIST("h_collision_counter"), 3.0);

        // d0 detector level jets associated to the current same collision
        const auto d0mcdJetsPerCollision = mcdjets.sliceBy(D0MCDJetsPerCollision, collision.globalIndex());
        for (const auto& mcdjet : d0mcdJetsPerCollision) {

          registry.fill(HIST("h_jet_counter"), 0.5);

          // obtain leading HF candidate in jet
          auto mcdd0cand = mcdjet.candidates_first_as<aod::CandidatesD0MCD>();

          if (mcdjet.has_matchedJetCand()) {
            registry.fill(HIST("h_jet_counter"), 1.5);
          }

          // reflection information for storage: +1 = D0, -1 = D0bar, 0 = neither
          int matchedFrom = 0;
          int decayChannel = 1 << aod::hf_cand_2prong::DecayType::D0ToPiK;

          if (mcdd0cand.flagMcMatchRec() == decayChannel) { // matched to D0 on truth level
            matchedFrom = 1;
          } else if (mcdd0cand.flagMcMatchRec() == -decayChannel) { // matched to D0bar on truth level
            matchedFrom = -1;
          } else { // matched to another kind of particle on truth level
            matchedFrom = 0;
          }

          // store data in MC detector level table
          mcddistJetTable(jetutilities::deltaR(mcdjet, mcdd0cand),
                          mcdjet.pt(), mcdjet.eta(), mcdjet.phi(), mcdjet.tracks_as<aod::JetTracks>().size(),                                                         // detector level jet
                          mcdd0cand.pt(), mcdd0cand.eta(), mcdd0cand.phi(), mcdd0cand.m(), mcdd0cand.y(), (mcdd0cand.originMcRec() == RecoDecay::OriginType::Prompt), // detector level D0 candidate
                          mcdjet.has_matchedJetCand(), mcdd0cand.mlScores()[0], mcdd0cand.mlScores()[1], mcdd0cand.mlScores()[2],                                     // ML scores for bkg, prompt and non-prompt
                          matchedFrom, mcdd0cand.candidateSelFlag());                                                                                                 // check whether detector level candidate is a reflection, CandidateSelFlag == 0 -> selected as D0, CandidateSelFlag == 1 -> selected as D0bar
        }
      }

      // d0 particle level jets associated to same mccollision
      const auto d0mcpJetsPerMCCollision = mcpjets.sliceBy(D0MCPJetsPerMCCollision, mccollision.globalIndex());
      for (const auto& mcpjet : d0mcpJetsPerMCCollision) {

        registry.fill(HIST("h_jet_counter"), 0.0);

        // obtain leading HF particle in jet
        auto mcpd0cand = mcpjet.candidates_first_as<aod::CandidatesD0MCP>();

        if (mcpjet.has_matchedJetCand()) {
          registry.fill(HIST("h_jet_counter"), 1.0);
        }

        // store data in MC detector level table (calculate angular distance in eta-phi plane on the fly)
        mcpdistJetTable(jetutilities::deltaR(mcpjet, mcpd0cand),
                        mcpjet.pt(), mcpjet.eta(), mcpjet.phi(), mcpjet.tracks_as<aod::JetParticles>().size(),                                       // particle level jet
                        mcpd0cand.pt(), mcpd0cand.eta(), mcpd0cand.phi(), mcpd0cand.y(), (mcpd0cand.originMcGen() == RecoDecay::OriginType::Prompt), // particle level D0
                        mcpjet.has_matchedJetCand());
      }
    }
  }
  PROCESS_SWITCH(HfFragmentationFunctionTask, processMcEfficiency, "non-matched and matched MC HF and jets", false);

  void processMcChargedMatched(aod::JetMcCollisions const& mccollisions,
                               aod::JetCollisionsMCD const& collisions,
                               JetMCDTable const& mcdjets,
                               JetMCPTable const&,
                               aod::CandidatesD0MCD const&,
                               aod::CandidatesD0MCP const&,
                               aod::JetTracks const&,
                               aod::JetParticles const&)
  {
    for (const auto& mccollision : mccollisions) {

      registry.fill(HIST("h_collision_counter"), 0.0);

      // skip collisions outside of |z| < vertexZCut
      if (abs(mccollision.posZ()) > vertexZCut) {
        continue;
      }
      registry.fill(HIST("h_collision_counter"), 1.0);

      // reconstructed collisions associated to same mccollision
      const auto collisionsPerMCCollision = collisions.sliceBy(CollisionsPerMCCollision, mccollision.globalIndex());
      for (const auto& collision : collisionsPerMCCollision) {

        registry.fill(HIST("h_collision_counter"), 2.0);
        if (!jetderiveddatautilities::selectCollision(collision, eventSelection) || !(abs(collision.posZ()) < vertexZCut)) {
          continue;
        }
        registry.fill(HIST("h_collision_counter"), 3.0);
        // d0 detector level jets associated to the current same collision
        const auto d0mcdJetsPerCollision = mcdjets.sliceBy(D0MCDJetsPerCollision, collision.globalIndex());
        for (const auto& mcdjet : d0mcdJetsPerCollision) {

          registry.fill(HIST("h_jet_counter"), 0.5);

          // comparison with fill on bin on 2.5 for sanity check
          if (mcdjet.has_matchedJetCand()) {
            registry.fill(HIST("h_jet_counter"), 1.5);
          }

          // obtain leading HF candidate in jet
          auto mcdd0cand = mcdjet.candidates_first_as<aod::CandidatesD0MCD>();

          // reflection information for storage: +1 = D0, -1 = D0bar, 0 = neither
          int matchedFrom = 0;
          int decayChannel = 1 << aod::hf_cand_2prong::DecayType::D0ToPiK;

          if (mcdd0cand.flagMcMatchRec() == decayChannel) { // matched to D0 on truth level
            matchedFrom = 1;
          } else if (mcdd0cand.flagMcMatchRec() == -decayChannel) { // matched to D0bar on truth level
            matchedFrom = -1;
          } else { // matched to another kind of particle on truth level
            matchedFrom = 0;
          }

          // loop through detector level matched to current particle level
          for (auto& mcpjet : mcdjet.matchedJetCand_as<JetMCPTable>()) {

            registry.fill(HIST("h_jet_counter"), 2.5);

            // obtain leading HF candidate in jet
            auto mcpd0cand = mcpjet.candidates_first_as<aod::CandidatesD0MCP>();

            // store matched particle and detector level data in one single table (calculate angular distance in eta-phi plane on the fly)
            matchJetTable(jetutilities::deltaR(mcpjet, mcpd0cand), mcpjet.pt(), mcpjet.eta(), mcpjet.phi(), mcpjet.tracks_as<aod::JetParticles>().size(),             // particle level jet
                          mcpd0cand.pt(), mcpd0cand.eta(), mcpd0cand.phi(), mcpd0cand.y(), (mcpd0cand.originMcGen() == RecoDecay::OriginType::Prompt),                // particle level D0
                          jetutilities::deltaR(mcdjet, mcdd0cand), mcdjet.pt(), mcdjet.eta(), mcdjet.phi(), mcdjet.tracks_as<aod::JetTracks>().size(),                // detector level jet
                          mcdd0cand.pt(), mcdd0cand.eta(), mcdd0cand.phi(), mcdd0cand.m(), mcdd0cand.y(), (mcdd0cand.originMcRec() == RecoDecay::OriginType::Prompt), // detector level D0
                          mcdd0cand.mlScores()[0], mcdd0cand.mlScores()[1], mcdd0cand.mlScores()[2],
                          matchedFrom, mcdd0cand.candidateSelFlag()); // check whether detector level candidate is a reflection, CandidateSelFlag == 0 -> selected as D0, CandidateSelFlag == 1 -> selected as D0bar
          }
        }
      }
    }
  }
  PROCESS_SWITCH(HfFragmentationFunctionTask, processMcChargedMatched, "matched MC HF and jets", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfFragmentationFunctionTask>(cfgc, TaskName{"jet-charm-hadronization"})};
}
