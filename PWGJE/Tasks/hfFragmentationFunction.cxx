// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file hfFragmentationFunction.cxx
/// \brief charm hadron hadronization task
/// \author Christian Reckziegel <christian.reckziegel@cern.ch>, Federal University of ABC
/// \since 15.03.2024
///
/// The task store data relevant to the calculation of hadronization observables radial
/// profile and/or jet momentum fraction for charmed hadrons

#include "PWGHF/Core/DecayChannels.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetHFUtilities.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TVector3.h>

#include <Rtypes.h>

#include <cstdlib>
#include <string>
#include <vector>

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

//
/// Collision counter selection indexes
///
/// The collision selection is done and stored in multiple steps, for later QA analysis.
/// In order not to hard code which bins should be filled throughout different process
/// function, this namespace with enums is create
namespace collisionSelections
{
enum CollisionSelectionStep {
  kMCCollisions = 0,                    ///< raw mccollisions with no selection, starts with 0
  kMCCollisionsZCut,                    ///< mccollisions with z vtx selection
  kMCCollisionsZCutSel8,                ///< mccollisions with z vtx and sel8 mc emulated selections
  kMCCollisionsZCutSel8HasCollisions,   ///< mccollisions with z vtx and sel8 mc emulated selections, with at least one reconstructed collisions
  kMCCollisionsZCutSel8SplitCollisions, ///< mccollisions with z vtx and sel8 mc emulated selections, with no split reconstructed collisions
  kRecoCollisions,                      ///< raw reconstructed collisions after previous mccollisions selection
  kRecoCollisionsZcut,                  ///< reconstructed collisions with z vtx selection after previous mccollisions selection
  kRecoCollisionsZcutSel8               ///< reconstructed collisions with z vtx and sel8 selections after previous mccollisions selection
};
}
// creating table for storing distance data
namespace o2::aod
{
namespace jet_distance
{
DECLARE_SOA_COLUMN(JetHfDist, jetHfDist, float);
DECLARE_SOA_COLUMN(JetPt, jetPt, float);
DECLARE_SOA_COLUMN(JetEta, jetEta, float);
DECLARE_SOA_COLUMN(JetPhi, jetPhi, float);
DECLARE_SOA_COLUMN(JetNConst, jetNConst, int);
DECLARE_SOA_COLUMN(HfPt, hfPt, float);
DECLARE_SOA_COLUMN(HfEta, hfEta, float);
DECLARE_SOA_COLUMN(HfPhi, hfPhi, float);
DECLARE_SOA_COLUMN(HfMass, hfMass, float);
DECLARE_SOA_COLUMN(HfY, hfY, float);
DECLARE_SOA_COLUMN(HfPrompt, hfPrompt, bool);
DECLARE_SOA_COLUMN(HfMatch, hfMatch, bool);
DECLARE_SOA_COLUMN(HfMlScore0, hfMlScore0, float);
DECLARE_SOA_COLUMN(HfMlScore1, hfMlScore1, float);
DECLARE_SOA_COLUMN(HfMlScore2, hfMlScore2, float);
DECLARE_SOA_COLUMN(HfMatchedFrom, hfMatchedFrom, int);
DECLARE_SOA_COLUMN(HfSelectedAs, hfSelectedAs, int);
DECLARE_SOA_COLUMN(McJetHfDist, mcJetHfDist, float);
DECLARE_SOA_COLUMN(McJetPt, mcJetPt, float);
DECLARE_SOA_COLUMN(McJetEta, mcJetEta, float);
DECLARE_SOA_COLUMN(McJetPhi, mcJetPhi, float);
DECLARE_SOA_COLUMN(McJetNConst, mcJetNConst, int);
DECLARE_SOA_COLUMN(McHfPt, mcHfPt, float);
DECLARE_SOA_COLUMN(McHfEta, mcHfEta, float);
DECLARE_SOA_COLUMN(McHfPhi, mcHfPhi, float);
DECLARE_SOA_COLUMN(McHfY, mcHfY, float);
DECLARE_SOA_COLUMN(McHfPrompt, mcHfPrompt, bool);
DECLARE_SOA_COLUMN(McHfMatch, mcHfMatch, bool);
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
                  jet_distance::McJetHfDist,
                  jet_distance::McJetPt,
                  jet_distance::McJetEta,
                  jet_distance::McJetPhi,
                  jet_distance::McJetNConst,
                  jet_distance::McHfPt,
                  jet_distance::McHfEta,
                  jet_distance::McHfPhi,
                  jet_distance::McHfY,
                  jet_distance::McHfPrompt,
                  jet_distance::McHfMatch);
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
                  jet_distance::McJetHfDist,
                  jet_distance::McJetPt,
                  jet_distance::McJetEta,
                  jet_distance::McJetPhi,
                  jet_distance::McJetNConst,
                  jet_distance::McHfPt,
                  jet_distance::McHfEta,
                  jet_distance::McHfPhi,
                  jet_distance::McHfY,
                  jet_distance::McHfPrompt,
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

struct HfFragmentationFunction {
  // producing new table
  Produces<aod::JetDistanceTable> distJetTable;
  Produces<aod::MCPJetDistanceTable> mcpdistJetTable;
  Produces<aod::MCDJetDistanceTable> mcddistJetTable;
  Produces<aod::MatchJetDistanceTable> matchJetTable;

  // Tables for MC jet matching
  using JetD0MCDTable = soa::Join<aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents, aod::D0ChargedMCDetectorLevelJetsMatchedToD0ChargedMCParticleLevelJets>;
  using JetD0MCPTable = soa::Join<aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents, aod::D0ChargedMCParticleLevelJetsMatchedToD0ChargedMCDetectorLevelJets>;
  using JetLcMCDTable = soa::Join<aod::LcChargedMCDetectorLevelJets, aod::LcChargedMCDetectorLevelJetConstituents, aod::LcChargedMCDetectorLevelJetsMatchedToLcChargedMCParticleLevelJets>;
  using JetLcMCPTable = soa::Join<aod::LcChargedMCParticleLevelJets, aod::LcChargedMCParticleLevelJetConstituents, aod::LcChargedMCParticleLevelJetsMatchedToLcChargedMCDetectorLevelJets>;

  // slices for accessing proper HF mcdjets collision associated to mccollisions
  PresliceUnsorted<aod::JetCollisionsMCD> collisionsPerMCCollisionPreslice = aod::jmccollisionlb::mcCollisionId;
  Preslice<JetD0MCDTable> d0MCDJetsPerCollisionPreslice = aod::jet::collisionId;
  Preslice<JetD0MCPTable> d0MCPJetsPerMCCollisionPreslice = aod::jet::mcCollisionId;
  Preslice<JetLcMCDTable> lcMCDJetsPerCollisionPreslice = aod::jet::collisionId;
  Preslice<JetLcMCPTable> lcMCPJetsPerMCCollisionPreslice = aod::jet::mcCollisionId;

  // Histogram registry: an object to hold your histograms
  HistogramRegistry registry{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<bool> applyMcEventSelection{"applyMcEventSelection", false, "Choose a boolean value"};
  Configurable<bool> applyRecoEventSelection{"applyRecoEventSelection", true, "Choose a boolean value"};
  Configurable<bool> rejectSplitCollisions{"rejectSplitCollisions", true, "reject generated events associated to more than one reconstructed collision"};

  std::vector<int> eventSelectionBits;

  void init(InitContext const&)
  {
    // initialise event selection:
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));

    // create histograms
    // collision counter histograms
    std::vector<std::string> histLabels = {"mccollisions", "mccollisions+z_cut", "mccollisions+z_{cut}+sel8", "mccollisions+z_{cut}+sel8+HasCollisions", "mccollisions+z_{cut}+sel8+NoSplitVtx", "collisions", "collisions+z_{cut}", "collisions+z_{cut}+sel8"};
    registry.add("h_collision_counter", ";# of collisions;", HistType::kTH1F, {{static_cast<int>(histLabels.size()), 0.0, static_cast<double>(histLabels.size())}});
    auto collCounter = registry.get<TH1>(HIST("h_collision_counter"));
    for (std::vector<std::string>::size_type iCounter = 0; iCounter < histLabels.size(); iCounter++) {
      collCounter->GetXaxis()->SetBinLabel(iCounter + 1, histLabels[iCounter].data());
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
    registry.add("h_hf_jet_projection", ";z^{HF,jet}_{||};dN/dz^{HF,jet}_{||}", {HistType::kTH1F, {{1000, 0., 10.}}});
    registry.add("h_hf_jet_distance_vs_projection", ";#DeltaR_{HF,jet};z^{HF,jet}_{||}", {HistType::kTH2F, {{1000, 0., 10.}, {1000, 0., 10.}}});
    registry.add("h_hf_jet_distance", ";#DeltaR_{HF,jet};dN/d(#DeltaR)", {HistType::kTH1F, {{1000, 0., 10.}}});
    registry.add("h_hf_jet_pt", ";p_{T,HF jet};dN/dp_{T,HF jet}", {HistType::kTH1F, {{200, 0., 10.}}});
    registry.add("h_hf_jet_eta", ";#eta_{T,HF jet};dN/d#eta_{HF jet}", {HistType::kTH1F, {{250, -5., 5.}}});
    registry.add("h_hf_jet_phi", ";#phi_{T,HF jet};dN/d#phi_{HF jet}", {HistType::kTH1F, {{250, -10., 10.}}});
    registry.add("h_hf_mass", ";m_{HF} (GeV/c^{2});dN/dm_{HF}", {HistType::kTH1F, {{1000, 0., 10.}}});
    registry.add("h_hf_eta", ";#eta_{HF} (GeV/c^{2});dN/d#eta_{HF}", {HistType::kTH1F, {{250, -5., 5.}}});
    registry.add("h_hf_phi", ";#phi_{HF} (GeV/c^{2});dN/d#phi_{HF}", {HistType::kTH1F, {{250, -10., 10.}}});
  }

  void processDummy(aod::TracksIU const&) {}
  PROCESS_SWITCH(HfFragmentationFunction, processDummy, "Dummy process function turned on by default", true);

  template <typename TJets, typename TCandidates>
  void analyzeData(aod::JetCollision const& collision,
                   TJets const& jets,
                   TCandidates const&,
                   aod::JetTracks const&)
  {
    // apply event selection and fill histograms for sanity check
    registry.fill(HIST("h_collision_counter"), collisionSelections::kRecoCollisions);
    if (applyRecoEventSelection && (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits) || !(std::abs(collision.posZ()) < vertexZCut))) {
      return;
    }
    registry.fill(HIST("h_collision_counter"), collisionSelections::kRecoCollisionsZcutSel8);

    for (const auto& jet : jets) {
      // fill jet counter histogram
      registry.fill(HIST("h_jet_counter"), 0.5);
      // obtaining jet 3-vector
      TVector3 jetVector(jet.px(), jet.py(), jet.pz());

      for (const auto& candidate : jet.template candidates_as<TCandidates>()) {

        // obtaining jet 3-vector
        TVector3 hadronMomentum(candidate.px(), candidate.py(), candidate.pz());

        // calculating fraction of the jet momentum carried by the HF hadron along the direction of the jet axis
        double zParallel = (jetVector * hadronMomentum) / (jetVector * jetVector);

        // calculating angular distance in eta-phi plane
        double axisDistance = jetutilities::deltaR(jet, candidate);

        // filling histograms
        registry.fill(HIST("h_hf_jet_projection"), zParallel);
        registry.fill(HIST("h_hf_jet_distance_vs_projection"), axisDistance, zParallel);
        registry.fill(HIST("h_hf_jet_distance"), axisDistance);
        registry.fill(HIST("h_hf_jet_pt"), jet.pt());
        registry.fill(HIST("h_hf_jet_eta"), jet.eta());
        registry.fill(HIST("h_hf_jet_phi"), jet.phi());
        registry.fill(HIST("h_hf_mass"), candidate.m());
        registry.fill(HIST("h_hf_eta"), candidate.eta());
        registry.fill(HIST("h_hf_phi"), candidate.phi());

        // filling table
        distJetTable(axisDistance,
                     jet.pt(), jet.eta(), jet.phi(), jet.template tracks_as<aod::JetTracks>().size() + jet.template candidates_as<TCandidates>().size(),
                     candidate.pt(), candidate.eta(), candidate.phi(), candidate.m(), candidate.y(), candidate.mlScores()[0], candidate.mlScores()[1], candidate.mlScores()[2]);

        break; // get out of candidates' loop after first HF particle is found in jet
      } // end of HF hadron candidates loop

    } // end of jets loop

  } // end of analyzeData function

  void processD0DataCharged(aod::JetCollision const& collision,
                            soa::Join<aod::D0ChargedJets, aod::D0ChargedJetConstituents> const& jets,
                            aod::CandidatesD0Data const& candidates,
                            aod::JetTracks const& jettracks)
  {
    analyzeData<soa::Join<aod::D0ChargedJets, aod::D0ChargedJetConstituents>, aod::CandidatesD0Data>(collision, jets, candidates, jettracks);
  }
  PROCESS_SWITCH(HfFragmentationFunction, processD0DataCharged, "Store kinematic charged D0 jet information from measured DATA", false);

  void processLcDataCharged(aod::JetCollision const& collision,
                            soa::Join<aod::LcChargedJets, aod::LcChargedJetConstituents> const& jets,
                            aod::CandidatesLcData const& candidates,
                            aod::JetTracks const& jettracks)
  {
    analyzeData<soa::Join<aod::LcChargedJets, aod::LcChargedJetConstituents>, aod::CandidatesLcData>(collision, jets, candidates, jettracks);
  }
  PROCESS_SWITCH(HfFragmentationFunction, processLcDataCharged, "Store kinematic charged Lc jet information from measured DATA", false);

  void processMcEfficiency(aod::JetMcCollisions const& mccollisions,
                           aod::JetCollisionsMCD const& collisions,
                           JetD0MCDTable const& mcdjets,
                           JetD0MCPTable const& mcpjets,
                           aod::CandidatesD0MCD const&,
                           aod::CandidatesD0MCP const&,
                           aod::JetTracks const&,
                           aod::JetParticles const&)
  {
    for (const auto& mccollision : mccollisions) {

      registry.fill(HIST("h_collision_counter"), collisionSelections::kMCCollisions);
      // skip collisions outside of |z| < vertexZCut
      if (applyMcEventSelection && (!jetderiveddatautilities::selectCollision(mccollision, eventSelectionBits) || !(std::abs(mccollision.posZ()) < vertexZCut))) {
        continue;
      }
      registry.fill(HIST("h_collision_counter"), collisionSelections::kMCCollisionsZCutSel8);

      // reconstructed collisions associated to same mccollision
      const auto collisionsPerMCCollision = collisions.sliceBy(collisionsPerMCCollisionPreslice, mccollision.globalIndex());
      for (const auto& collision : collisionsPerMCCollision) {

        registry.fill(HIST("h_collision_counter"), collisionSelections::kRecoCollisions);
        if (applyRecoEventSelection && (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits) || !(std::abs(collision.posZ()) < vertexZCut))) {
          continue;
        }
        registry.fill(HIST("h_collision_counter"), collisionSelections::kRecoCollisionsZcutSel8);

        // d0 detector level jets associated to the current same collision
        const auto d0mcdJetsPerCollision = mcdjets.sliceBy(d0MCDJetsPerCollisionPreslice, collision.globalIndex());
        for (const auto& mcdjet : d0mcdJetsPerCollision) {

          registry.fill(HIST("h_jet_counter"), 0.5);

          // obtain leading HF candidate in jet
          auto mcdd0cand = mcdjet.candidates_first_as<aod::CandidatesD0MCD>();

          if (mcdjet.has_matchedJetCand()) {
            registry.fill(HIST("h_jet_counter"), 1.5);
          }

          // reflection information for storage: D0 = +1, D0bar = -1, neither = 0
          int selectedAs = 0;

          // bitwise AND operation: Checks whether BIT(i) is set, regardless of other bits
          if (mcdd0cand.candidateSelFlag() & BIT(0)) { // CandidateSelFlag == BIT(0) -> selected as D0
            selectedAs = 1;
          } else if (mcdd0cand.candidateSelFlag() & BIT(1)) { // CandidateSelFlag == BIT(1) -> selected as D0bar
            selectedAs = -1;
          }

          // store data in MC detector level table
          mcddistJetTable(jetutilities::deltaR(mcdjet, mcdd0cand),
                          mcdjet.pt(), mcdjet.eta(), mcdjet.phi(), mcdjet.tracks_as<aod::JetTracks>().size() + mcdjet.candidates_as<aod::CandidatesD0MCD>().size(),   // detector level jet
                          mcdd0cand.pt(), mcdd0cand.eta(), mcdd0cand.phi(), mcdd0cand.m(), mcdd0cand.y(), (mcdd0cand.originMcRec() == RecoDecay::OriginType::Prompt), // detector level D0 candidate
                          mcdjet.has_matchedJetCand(), mcdd0cand.mlScores()[0], mcdd0cand.mlScores()[1], mcdd0cand.mlScores()[2],                                     // Machine Learning PID scores: background, prompt, non-prompt
                          static_cast<int>(mcdd0cand.flagMcMatchRec()), selectedAs);                                                                                  // +1/-1 = D0(bar)→Kπ, ±2..5 = other D0 channels, 0 = no match
        }
      }

      // d0 particle level jets associated to same mccollision
      const auto d0mcpJetsPerMCCollision = mcpjets.sliceBy(d0MCPJetsPerMCCollisionPreslice, mccollision.globalIndex());
      for (const auto& mcpjet : d0mcpJetsPerMCCollision) {

        registry.fill(HIST("h_jet_counter"), 0.0);

        // obtain leading HF particle in jet
        auto mcpd0cand = mcpjet.candidates_first_as<aod::CandidatesD0MCP>();

        if (mcpjet.has_matchedJetCand()) {
          registry.fill(HIST("h_jet_counter"), 1.0);
        }

        // store data in MC detector level table (calculate angular distance in eta-phi plane on the fly)
        mcpdistJetTable(jetutilities::deltaR(mcpjet, mcpd0cand),
                        mcpjet.pt(), mcpjet.eta(), mcpjet.phi(), mcpjet.tracks_as<aod::JetParticles>().size() + mcpjet.candidates_as<aod::CandidatesD0MCP>().size(), // particle level jet
                        mcpd0cand.pt(), mcpd0cand.eta(), mcpd0cand.phi(), mcpd0cand.y(), (mcpd0cand.originMcGen() == RecoDecay::OriginType::Prompt),                 // particle level D0
                        mcpjet.has_matchedJetCand());
      }
    }
  }
  PROCESS_SWITCH(HfFragmentationFunction, processMcEfficiency, "non-matched and matched MC HF and jets", false);

  template <typename TMCPJetsPerMCCollisionPreslice, typename TMCDJetsPerCollisionPreslice, typename TJetsMCP, typename TJetsMCD, typename TCandidatesMCP, typename TCandidatesMCD>
  void analyzeMC(TMCPJetsPerMCCollisionPreslice const& MCPJetsPerMCCollisionPreslice,
                 TMCDJetsPerCollisionPreslice const& MCDJetsPerCollisionPreslice,
                 aod::JetMcCollisions const& mccollisions,
                 aod::JetCollisionsMCD const& collisions,
                 TJetsMCP const& mcpjets,
                 TJetsMCD const& mcdjets,
                 TCandidatesMCP const&,
                 TCandidatesMCD const&,
                 aod::JetParticles const&,
                 aod::JetTracks const&)
  {
    for (const auto& mccollision : mccollisions) {

      // --- begin event selection
      registry.fill(HIST("h_collision_counter"), collisionSelections::kMCCollisions);
      // skip collisions outside of |z| < vertexZCut
      if (applyMcEventSelection && !(std::abs(mccollision.posZ()) < vertexZCut)) {
        continue;
      }
      registry.fill(HIST("h_collision_counter"), collisionSelections::kMCCollisionsZCut);
      if (applyMcEventSelection && !jetderiveddatautilities::selectCollision(mccollision, eventSelectionBits)) {
        continue;
      }
      registry.fill(HIST("h_collision_counter"), collisionSelections::kMCCollisionsZCutSel8);

      // reconstructed collisions associated to same mccollision
      const auto collisionsPerMCCollision = collisions.sliceBy(collisionsPerMCCollisionPreslice, mccollision.globalIndex());
      // only consider events with at least one reconstructed collision
      if (collisionsPerMCCollision.size() == 0) {
        continue;
      }
      // only consider events with no split vertices (one mccollision-to-one collision)
      registry.fill(HIST("h_collision_counter"), collisionSelections::kMCCollisionsZCutSel8HasCollisions);
      if (rejectSplitCollisions && collisionsPerMCCollision.size() > 1) {
        continue;
      }
      registry.fill(HIST("h_collision_counter"), collisionSelections::kMCCollisionsZCutSel8SplitCollisions);

      bool hasSelectedCollision = false;
      for (const auto& collision : collisionsPerMCCollision) {

        registry.fill(HIST("h_collision_counter"), collisionSelections::kRecoCollisions);
        if (applyRecoEventSelection && !(std::abs(collision.posZ()) < vertexZCut)) {
          continue;
        }
        registry.fill(HIST("h_collision_counter"), collisionSelections::kRecoCollisionsZcut);
        if (applyRecoEventSelection && !jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
          continue;
        }
        registry.fill(HIST("h_collision_counter"), collisionSelections::kRecoCollisionsZcutSel8);
        hasSelectedCollision = true;
      } // end of collisions loop

      if (!hasSelectedCollision) {
        continue;
      }

      // --- begin particle level jets storage
      // hf particle level jets associated to same mccollision
      const auto mcpJetsPerMCCollision = mcpjets.sliceBy(MCPJetsPerMCCollisionPreslice, mccollision.globalIndex());
      for (const auto& mcpjet : mcpJetsPerMCCollision) {

        registry.fill(HIST("h_jet_counter"), 0.0);

        // obtain leading HF particle in jet
        auto mcpcand = mcpjet.template candidates_first_as<TCandidatesMCP>();

        if (mcpjet.has_matchedJetCand()) {
          registry.fill(HIST("h_jet_counter"), 1.0);

          // loop over detector level matched to current particle level
          for (const auto& mcdjet : mcpjet.template matchedJetCand_as<TJetsMCD>()) {
            registry.fill(HIST("h_jet_counter"), 2.0);

            // obtain leading HF candidate in jet
            auto mcdcand = mcdjet.template candidates_first_as<TCandidatesMCD>();

            int selectedAs = 0;

            // bitwise AND operation: Checks whether BIT(i) is set, regardless of other bits
            if (mcdcand.candidateSelFlag() & BIT(0)) { // CandidateSelFlag == BIT(0) -> selected as HF
              selectedAs = 1;
            } else if (mcdcand.candidateSelFlag() & BIT(1)) { // CandidateSelFlag == BIT(1) -> selected as HFbar
              selectedAs = -1;
            }

            // store matched particle and detector level data in one single table (calculate angular distance in eta-phi plane on the fly)
            matchJetTable(jetutilities::deltaR(mcpjet, mcpcand), mcpjet.pt(), mcpjet.eta(), mcpjet.phi(), mcpjet.template tracks_as<aod::JetParticles>().size() + mcpjet.template candidates_as<TCandidatesMCP>().size(), // particle level jet
                          mcpcand.pt(), mcpcand.eta(), mcpcand.phi(), mcpcand.y(), (mcpcand.originMcGen() == RecoDecay::OriginType::Prompt),                                                                              // particle level HF
                          jetutilities::deltaR(mcdjet, mcdcand), mcdjet.pt(), mcdjet.eta(), mcdjet.phi(), mcdjet.template tracks_as<aod::JetTracks>().size() + mcdjet.template candidates_as<TCandidatesMCD>().size(),    // detector level jet
                          mcdcand.pt(), mcdcand.eta(), mcdcand.phi(), mcdcand.m(), mcdcand.y(), (mcdcand.originMcRec() == RecoDecay::OriginType::Prompt),                                                                 // detector level HF
                          mcdcand.mlScores()[0], mcdcand.mlScores()[1], mcdcand.mlScores()[2],                                                                                                                            // Machine Learning PID scores: background, prompt, non-prompt
                          static_cast<int>(mcdcand.flagMcMatchRec()), selectedAs);                                                                                                                                        // HF = +1, HFbar = -1, neither = 0
          }
        } else {
          // store matched particle and detector level data in one single table (calculate angular distance in eta-phi plane on the fly)
          matchJetTable(jetutilities::deltaR(mcpjet, mcpcand), mcpjet.pt(), mcpjet.eta(), mcpjet.phi(), mcpjet.template tracks_as<aod::JetParticles>().size() + mcpjet.template candidates_as<TCandidatesMCP>().size(), // particle level jet
                        mcpcand.pt(), mcpcand.eta(), mcpcand.phi(), mcpcand.y(), (mcpcand.originMcGen() == RecoDecay::OriginType::Prompt),                                                                              // particle level HF
                        -2, -2, -2, -2, -2,                                                                                                                                                                             // no detector-level jet found
                        -2, -2, -2, -2, -2, false,                                                                                                                                                                      // no detector-level jet found
                        -2, -2, -2,                                                                                                                                                                                     // no detector-level jet found
                        -2, -2);                                                                                                                                                                                        // no detector-level jet found
        }
      } // end of mcpjets loop

      // --- begin non-matched detector level jets storage (fake candidates and correlated background if present)
      // reconstructed collisions associated to same mccollision
      for (const auto& collision : collisionsPerMCCollision) {

        // Also apply the reconstructed level collisions selections
        if (applyRecoEventSelection && (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits) || !(std::abs(collision.posZ()) < vertexZCut))) {
          continue;
        }

        // d0 detector level jets associated to the current same collision
        const auto mcdJetsPerCollision = mcdjets.sliceBy(MCDJetsPerCollisionPreslice, collision.globalIndex());
        for (const auto& mcdjet : mcdJetsPerCollision) {

          registry.fill(HIST("h_jet_counter"), 0.5);

          // obtain leading HF candidate in jet
          auto mcdcand = mcdjet.template candidates_first_as<TCandidatesMCD>();

          if (mcdjet.has_matchedJetCand()) {
            registry.fill(HIST("h_jet_counter"), 1.5);
          } else { // store the detector level non-matched candidates

            // reflection information for storage: D0 = +1, D0bar = -1, neither = 0
            int selectedAs = 0;

            // bitwise AND operation: Checks whether BIT(i) is set, regardless of other bits
            if (mcdcand.candidateSelFlag() & BIT(0)) { // CandidateSelFlag == BIT(0) -> selected as D0
              selectedAs = 1;
            } else if (mcdcand.candidateSelFlag() & BIT(1)) { // CandidateSelFlag == BIT(1) -> selected as D0bar
              selectedAs = -1;
            }

            // store matched particle and detector level data in one single table (calculate angular distance in eta-phi plane on the fly)
            matchJetTable(-2, -2, -2, -2, -2,                                                                                                                                                                          // particle level jet
                          -2, -2, -2, -2, false,                                                                                                                                                                       // particle level HF
                          jetutilities::deltaR(mcdjet, mcdcand), mcdjet.pt(), mcdjet.eta(), mcdjet.phi(), mcdjet.template tracks_as<aod::JetTracks>().size() + mcdjet.template candidates_as<TCandidatesMCD>().size(), // detector level jet
                          mcdcand.pt(), mcdcand.eta(), mcdcand.phi(), mcdcand.m(), mcdcand.y(), (mcdcand.originMcRec() == RecoDecay::OriginType::Prompt),                                                              // detector level HF
                          mcdcand.mlScores()[0], mcdcand.mlScores()[1], mcdcand.mlScores()[2],                                                                                                                         // Machine Learning PID scores: background, prompt, non-prompt
                          static_cast<int>(mcdcand.flagMcMatchRec()), selectedAs);                                                                                                                                     // HF = +1, HFbar = -1, neither = 0
          }
        } // end of non-matched detector level jets loop
      } // end of collisions loop
    } // end of mccollisions loop
  } // end of analyzeMC function

  void processD0MC(aod::JetMcCollisions const& mccollisions,
                   aod::JetCollisionsMCD const& collisions,
                   JetD0MCPTable const& mcpjets,
                   JetD0MCDTable const& mcdjets,
                   aod::CandidatesD0MCP const& mcpcands,
                   aod::CandidatesD0MCD const& mcdcands,
                   aod::JetParticles const& jetparticles,
                   aod::JetTracks const& jettracks)
  {
    analyzeMC<Preslice<JetD0MCPTable>, Preslice<JetD0MCDTable>, JetD0MCPTable, JetD0MCDTable, aod::CandidatesD0MCP, aod::CandidatesD0MCD>(d0MCPJetsPerMCCollisionPreslice, d0MCDJetsPerCollisionPreslice, mccollisions, collisions, mcpjets, mcdjets, mcpcands, mcdcands, jetparticles, jettracks);
  }
  PROCESS_SWITCH(HfFragmentationFunction, processD0MC, "Store all simulated D0 jets information with matched candidate (if any found)", false);

  void processLcMC(aod::JetMcCollisions const& mccollisions,
                   aod::JetCollisionsMCD const& collisions,
                   JetLcMCPTable const& mcpjets,
                   JetLcMCDTable const& mcdjets,
                   aod::CandidatesLcMCP const& mcpcands,
                   aod::CandidatesLcMCD const& mcdcands,
                   aod::JetParticles const& jetparticles,
                   aod::JetTracks const& jettracks)
  {
    analyzeMC<Preslice<JetLcMCPTable>, Preslice<JetLcMCDTable>, JetLcMCPTable, JetLcMCDTable, aod::CandidatesLcMCP, aod::CandidatesLcMCD>(lcMCPJetsPerMCCollisionPreslice, lcMCDJetsPerCollisionPreslice, mccollisions, collisions, mcpjets, mcdjets, mcpcands, mcdcands, jetparticles, jettracks);
  }
  PROCESS_SWITCH(HfFragmentationFunction, processLcMC, "Store all simulated Lc jets information with matched candidate (if any found)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfFragmentationFunction>(cfgc)};
}
