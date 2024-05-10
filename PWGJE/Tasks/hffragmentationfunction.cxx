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

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// calculate delta phi such that 0 < delta phi < 2*pi
double DeltaPhi(double phi1, double phi2)
{
  // Compute the absolute difference between phi1 and phi2
  double dphi = std::abs(phi1 - phi2);
  if (dphi > M_PI) {
    // subtract 2pi if the difference if bigger than pi
    dphi = dphi - 2 * M_PI;
  }

  return dphi;
}

// creating table for storing distance data
namespace o2::aod
{
namespace DistanceSpace
{
DECLARE_SOA_COLUMN(JetHfDist, jethfdist, float);
DECLARE_SOA_COLUMN(JetPt, jetpt, float);
DECLARE_SOA_COLUMN(JetEta, jeteta, float);
DECLARE_SOA_COLUMN(JetPhi, jetphi, float);
DECLARE_SOA_COLUMN(HfPt, hfpt, float);
DECLARE_SOA_COLUMN(HfEta, hfeta, float);
DECLARE_SOA_COLUMN(HfPhi, hfphi, float);
DECLARE_SOA_COLUMN(HfMass, hfmass, float);
DECLARE_SOA_COLUMN(HfY, hfy, float);
} // namespace DistanceSpace
DECLARE_SOA_TABLE(JetDistanceTable, "AOD", "JETDISTTABLE",
                  DistanceSpace::JetHfDist,
                  DistanceSpace::JetPt,
                  DistanceSpace::JetEta,
                  DistanceSpace::JetPhi,
                  DistanceSpace::HfPt,
                  DistanceSpace::HfEta,
                  DistanceSpace::HfPhi,
                  DistanceSpace::HfMass,
                  DistanceSpace::HfY);
namespace MCPDistanceTable
{
DECLARE_SOA_COLUMN(JetHfDist, jethfdist, float);
DECLARE_SOA_COLUMN(JetPt, jetpt, float);
DECLARE_SOA_COLUMN(JetEta, jeteta, float);
DECLARE_SOA_COLUMN(JetPhi, jetphi, float);
DECLARE_SOA_COLUMN(HfPt, hfpt, float);
DECLARE_SOA_COLUMN(HfEta, hfeta, float);
DECLARE_SOA_COLUMN(HfPhi, hfphi, float);
DECLARE_SOA_COLUMN(HfY, hfy, float);
DECLARE_SOA_COLUMN(HfMatch, hfmatch, bool);
} // namespace MCPDistanceTable
DECLARE_SOA_TABLE(MCPJetDistanceTable, "AOD", "MCPJETDISTTABLE",
                  MCPDistanceTable::JetHfDist,
                  MCPDistanceTable::JetPt,
                  MCPDistanceTable::JetEta,
                  MCPDistanceTable::JetPhi,
                  MCPDistanceTable::HfPt,
                  MCPDistanceTable::HfEta,
                  MCPDistanceTable::HfPhi,
                  MCPDistanceTable::HfY,
                  MCPDistanceTable::HfMatch);
namespace MCDDistanceTable
{
DECLARE_SOA_COLUMN(JetHfDist, jethfdist, float);
DECLARE_SOA_COLUMN(JetPt, jetpt, float);
DECLARE_SOA_COLUMN(JetEta, jeteta, float);
DECLARE_SOA_COLUMN(JetPhi, jetphi, float);
DECLARE_SOA_COLUMN(HfPt, hfpt, float);
DECLARE_SOA_COLUMN(HfEta, hfeta, float);
DECLARE_SOA_COLUMN(HfPhi, hfphi, float);
DECLARE_SOA_COLUMN(HfMass, hfmass, float);
DECLARE_SOA_COLUMN(HfY, hfy, float);
DECLARE_SOA_COLUMN(HfMatch, hfmatch, bool);
} // namespace MCDDistanceTable
DECLARE_SOA_TABLE(MCDJetDistanceTable, "AOD", "MCDJETDISTTABLE",
                  MCDDistanceTable::JetHfDist,
                  MCDDistanceTable::JetPt,
                  MCDDistanceTable::JetEta,
                  MCDDistanceTable::JetPhi,
                  MCDDistanceTable::HfPt,
                  MCDDistanceTable::HfEta,
                  MCDDistanceTable::HfPhi,
                  MCDDistanceTable::HfMass,
                  MCDDistanceTable::HfY,
                  MCDDistanceTable::HfMatch);
} // namespace o2::aod

struct HfFragmentationFunctionTask {
  // producing new table
  Produces<aod::JetDistanceTable> distJetTable;
  Produces<aod::MCPJetDistanceTable> mcpdistJetTable;
  Produces<aod::MCDJetDistanceTable> mcddistJetTable;

  // Tables for MC jet matching
  using JetMCDTable = soa::Join<aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents, aod::D0ChargedMCDetectorLevelJetsMatchedToD0ChargedMCParticleLevelJets>;
  using JetMCPTable = soa::Join<aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents, aod::D0ChargedMCParticleLevelJetsMatchedToD0ChargedMCDetectorLevelJets>;

  // slices for accessing proper HF mcdjets collision associated to mccollisions
  PresliceUnsorted<JetCollisionsMCD> CollisionsPerMCCollision = aod::jmccollisionlb::mcCollisionId;
  Preslice<JetMCDTable> D0MCDJetsPerCollision = aod::jet::collisionId;
  Preslice<JetMCPTable> D0MCPJetsPerMCCollision = aod::jet::mcCollisionId;

  // Histogram registry: an object to hold your histograms
  HistogramRegistry registry{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {

    // create histograms
    // D0 candidate histograms from data
    registry.add("h_jet_counter", ";# jets;", {HistType::kTH1F, {{2, 0., 1.}}});
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

  void processDataChargedSubstructure(JetCollision const&,
                                      soa::Join<aod::D0ChargedJets, aod::D0ChargedJetConstituents> const& jets,
                                      JetTracks const&,
                                      CandidatesD0Data const&)
  {

    double axisDistance = 0;

    for (auto& jet : jets) {
      // fill jet counter histogram
      registry.fill(HIST("h_jet_counter"), 0.5);
      // obtaining jet 3-vector
      TVector3 jetVector(jet.px(), jet.py(), jet.pz());

      for (auto& d0Candidate : jet.hfcandidates_as<CandidatesD0Data>()) { // for jet constituents use -> auto& jetConstituent : jet.tracks_as<JetTracks>()

        // obtaining jet 3-vector
        TVector3 d0Vector(d0Candidate.px(), d0Candidate.py(), d0Candidate.pz());
        // calculating fraction of the jet momentum carried by the D0 along the direction of the jet axis
        double z_parallel = (jetVector * d0Vector) / (jetVector * jetVector);

        // calculating angular distance in eta-phi plane
        axisDistance = sqrt(pow(jet.eta() - d0Candidate.eta(), 2) + pow(DeltaPhi(jet.phi(), d0Candidate.phi()), 2));

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
        distJetTable(axisDistance, jet.pt(), jet.eta(), jet.phi(), d0Candidate.pt(), d0Candidate.eta(), d0Candidate.phi(), d0Candidate.m(), d0Candidate.y());

        break; // get out of candidates' loop after first HF particle is found in jet
      }        // end of D0 candidates loop

    } // end of jets loop

  } // end of process function
  PROCESS_SWITCH(HfFragmentationFunctionTask, processDataChargedSubstructure, "charged HF jet substructure", false);

  void processMcChargedMatched(JetMcCollisions const& mccollisions,
                               JetCollisionsMCD const& collisions,
                               JetMCDTable const& mcdjets, JetMCPTable const& mcpjets,
                               CandidatesD0MCD const& mcdd0cands, CandidatesD0MCP const& mcpd0cands,
                               JetTracks const& tracks, JetParticles const& particles)
  {
    double axisDistance = 0;

    for (auto& mccollision : mccollisions) {

      // reconstructed collisions associated to same mccollision
      const auto collisionsPerMCCollision = collisions.sliceBy(CollisionsPerMCCollision, mccollision.globalIndex());
      for (auto& collision : collisionsPerMCCollision) {
        // d0 detector level jets associated to the current same collision
        const auto d0mcdJetsPerCollision = mcdjets.sliceBy(D0MCDJetsPerCollision, collision.globalIndex());
        for (auto& mcdjet : d0mcdJetsPerCollision) {

          // obtain leading HF candidate in jet (is this correct?)
          auto mcdd0cand = mcdjet.template hfcandidates_first_as<CandidatesD0MCD>();

          // calculating angular distance in eta-phi plane
          axisDistance = sqrt(pow(mcdjet.eta() - mcdd0cand.eta(), 2) + pow(DeltaPhi(mcdjet.phi(), mcdd0cand.phi()), 2));

          // store data in MC detector level table
          mcddistJetTable(axisDistance, mcdjet.pt(), mcdjet.eta(), mcdjet.phi(), mcdd0cand.pt(), mcdd0cand.eta(), mcdd0cand.phi(), mcdd0cand.m(), mcdd0cand.y(), mcdjet.has_matchedJetCand());
        }
      }

      // d0 particle level jets associated to same mccollision
      const auto d0mcpJetsPerMCCollision = mcpjets.sliceBy(D0MCPJetsPerMCCollision, mccollision.globalIndex());
      for (auto& mcpjet : d0mcpJetsPerMCCollision) {

        // obtain leading HF particle in jet
        auto mcpd0cand = mcpjet.template hfcandidates_first_as<CandidatesD0MCP>();

        // calculating angular distance in eta-phi plane
        axisDistance = sqrt(pow(mcpjet.eta() - mcpd0cand.eta(), 2) + pow(DeltaPhi(mcpjet.phi(), mcpd0cand.phi()), 2));

        // store data in MC detector level table
        mcpdistJetTable(axisDistance, mcpjet.pt(), mcpjet.eta(), mcpjet.phi(), mcpd0cand.pt(), mcpd0cand.eta(), mcpd0cand.phi(), mcpd0cand.y(), mcpjet.has_matchedJetCand());
      }
    }
  }
  PROCESS_SWITCH(HfFragmentationFunctionTask, processMcChargedMatched, "matched MC HF jets", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfFragmentationFunctionTask>(cfgc, TaskName{"jet-charm-hadronization"})};
}
