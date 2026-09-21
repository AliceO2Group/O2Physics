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
/// \file jetDsSpecSubs.cxx
/// \brief Ds-tagged jet analysis with substructure histogram outputs
/// \author Monalisa Melo <monalisa.melo@cern.ch>, Universidade de São Paulo

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"

#include "Common/Core/RecoDecay.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/ConfigContext.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <THnSparse.h>
#include <TVector3.h>

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

consteval float getValFromBin(int bin)
{
  return static_cast<float>(bin) - 0.5f;
}

enum BinExpColCntr { AllCollisions = 1,
                     Sel8ZCut = 2 };

enum BinMCColCntr { AllMCCollisions = 1,
                    SelectedMCCollisions = 2,
                    AssociatedRecoCollisions = 3,
                    SelectedAssociatedRecoCollisions = 4
};

enum BinMCJetCntr {
  MCPJets = 1,
  MCPJetsWithMCDMatch = 2,
  MCDJets = 3,
  MCDJetsWithMCPMatch = 4,
  MatchedPairs = 5
};

//==============================================================
// AOD tables for Ds-tagged jet observables
//==============================================================
namespace o2::aod
{
namespace jet_ds_output
{
// Jet observables
DECLARE_SOA_COLUMN(JetPt, jetPt, float);
DECLARE_SOA_COLUMN(JetEta, jetEta, float);
DECLARE_SOA_COLUMN(JetPhi, jetPhi, float);
DECLARE_SOA_COLUMN(JetNConst, jetNConst, int);
// Ds observables
DECLARE_SOA_COLUMN(DsPt, dsPt, float);
DECLARE_SOA_COLUMN(DsEta, dsEta, float);
DECLARE_SOA_COLUMN(DsPhi, dsPhi, float);
DECLARE_SOA_COLUMN(DsMass, dsMass, float);
// Ds-jet observables
DECLARE_SOA_COLUMN(ZParallel, zParallel, float);
DECLARE_SOA_COLUMN(DeltaR, deltaR, float);
// MC information
DECLARE_SOA_COLUMN(DsOrigin, dsOrigin, int);
DECLARE_SOA_COLUMN(IsMatched, isMatched, bool);
// Matched MCP <-> MCD observables
// Particle level
DECLARE_SOA_COLUMN(McpJetPt, mcpJetPt, float);
DECLARE_SOA_COLUMN(McpJetEta, mcpJetEta, float);
DECLARE_SOA_COLUMN(McpJetPhi, mcpJetPhi, float);
DECLARE_SOA_COLUMN(McpJetNConst, mcpJetNConst, int);
DECLARE_SOA_COLUMN(McpDsPt, mcpDsPt, float);
DECLARE_SOA_COLUMN(McpDsEta, mcpDsEta, float);
DECLARE_SOA_COLUMN(McpDsPhi, mcpDsPhi, float);
DECLARE_SOA_COLUMN(McpZParallel, mcpZParallel, float);
DECLARE_SOA_COLUMN(McpDeltaR, mcpDeltaR, float);

DECLARE_SOA_COLUMN(McdJetPt, mcdJetPt, float);
DECLARE_SOA_COLUMN(McdJetEta, mcdJetEta, float);
DECLARE_SOA_COLUMN(McdJetPhi, mcdJetPhi, float);
DECLARE_SOA_COLUMN(McdJetNConst, mcdJetNConst, int);
DECLARE_SOA_COLUMN(McdDsPt, mcdDsPt, float);
DECLARE_SOA_COLUMN(McdDsEta, mcdDsEta, float);
DECLARE_SOA_COLUMN(McdDsPhi, mcdDsPhi, float);
DECLARE_SOA_COLUMN(McdDsMass, mcdDsMass, float);
DECLARE_SOA_COLUMN(McdZParallel, mcdZParallel, float);
DECLARE_SOA_COLUMN(McdDeltaR, mcdDeltaR, float);
// MC truth matching information
DECLARE_SOA_COLUMN(DsFlagMcMatchRec, dsFlagMcMatchRec, int);
DECLARE_SOA_COLUMN(DsFlagMcMatchGen, dsFlagMcMatchGen, int);
// MC truth information
DECLARE_SOA_COLUMN(MatchedDsOrigin, matchedDsOrigin, int);
} // namespace jet_ds_output
// Detector-level Ds-tagged jets
DECLARE_SOA_TABLE(DsMCDJetTable, "AOD", "DSMCDJET",
                  jet_ds_output::JetPt,
                  jet_ds_output::JetEta,
                  jet_ds_output::JetPhi,
                  jet_ds_output::JetNConst,
                  jet_ds_output::DsPt,
                  jet_ds_output::DsEta,
                  jet_ds_output::DsPhi,
                  jet_ds_output::DsMass,
                  jet_ds_output::ZParallel,
                  jet_ds_output::DeltaR,
                  jet_ds_output::DsOrigin,
                  jet_ds_output::DsFlagMcMatchRec,
                  jet_ds_output::IsMatched);
// Particle-level Ds-tagged jets
DECLARE_SOA_TABLE(DsMCPJetTable, "AOD", "DSMCPJET",
                  jet_ds_output::JetPt,
                  jet_ds_output::JetEta,
                  jet_ds_output::JetPhi,
                  jet_ds_output::JetNConst,
                  jet_ds_output::DsPt,
                  jet_ds_output::DsEta,
                  jet_ds_output::DsPhi,
                  jet_ds_output::ZParallel,
                  jet_ds_output::DeltaR,
                  jet_ds_output::DsOrigin,
                  jet_ds_output::DsFlagMcMatchGen,
                  jet_ds_output::IsMatched);
// Matched particle-level <-> detector-level Ds-tagged jet pairs
DECLARE_SOA_TABLE(DsMatchedJetTable, "AOD", "DSMATCHJET",
                  // Particle level
                  jet_ds_output::McpJetPt,
                  jet_ds_output::McpJetEta,
                  jet_ds_output::McpJetPhi,
                  jet_ds_output::McpJetNConst,
                  jet_ds_output::McpDsPt,
                  jet_ds_output::McpDsEta,
                  jet_ds_output::McpDsPhi,
                  jet_ds_output::McpZParallel,
                  jet_ds_output::McpDeltaR,
                  // Detector level
                  jet_ds_output::McdJetPt,
                  jet_ds_output::McdJetEta,
                  jet_ds_output::McdJetPhi,
                  jet_ds_output::McdJetNConst,
                  jet_ds_output::McdDsPt,
                  jet_ds_output::McdDsEta,
                  jet_ds_output::McdDsPhi,
                  jet_ds_output::McdDsMass,
                  jet_ds_output::McdZParallel,
                  jet_ds_output::McdDeltaR,
                  // Truth
                  jet_ds_output::MatchedDsOrigin);
} // namespace o2::aod

struct JetDsSpecSubs {

  //==================
  // Type definitions
  //==================

  Produces<aod::DsMCDJetTable> outputMCDJets;
  Produces<aod::DsMCPJetTable> outputMCPJets;
  Produces<aod::DsMatchedJetTable> outputMatchedJets;

  using DsCandidatesData = aod::CandidatesDsData;
  using DsCandidatesMCD = aod::CandidatesDsMCD;
  using DsCandidatesMCP = aod::CandidatesDsMCP;

  using DsDataJets = soa::Join<aod::DsChargedJets, aod::DsChargedJetConstituents>;
  using DsMCDJets = soa::Join<aod::DsChargedMCDetectorLevelJets, aod::DsChargedMCDetectorLevelJetConstituents, aod::DsChargedMCDetectorLevelJetsMatchedToDsChargedMCParticleLevelJets>;
  using DsMCPJets = soa::Join<aod::DsChargedMCParticleLevelJets, aod::DsChargedMCParticleLevelJetConstituents, aod::DsChargedMCParticleLevelJetsMatchedToDsChargedMCDetectorLevelJets>;

  using DsMCPJetsOnTheFly = soa::Join<aod::DsChargedMCParticleLevelJets, aod::DsChargedMCParticleLevelJetConstituents>;

  using DsDataJetsEWS = soa::Join<
    aod::DsChargedEventWiseSubtractedJets,
    aod::DsChargedEventWiseSubtractedJetConstituents>;

  // Inclusive charged jets
  using ChargedJets = soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>;
  using ChargedJetsEWS = soa::Join<aod::ChargedEventWiseSubtractedJets, aod::ChargedEventWiseSubtractedJetConstituents>;

  // Slices for access to proper HF MCD jet collision that is associated to MCCollision
  PresliceUnsorted<aod::JetCollisionsMCD> collisionsPerMCCollisionPreslice = aod::jmccollisionlb::mcCollisionId;
  Preslice<DsMCDJets> dsMCDJetsPerEXPCollisionPreslice = aod::jet::collisionId;
  // Preslice<DsMCDJetsEWS> dsMCDJetsEWSPerEXPCollisionPreslice = aod::jet::collisionId;
  Preslice<DsMCPJets> dsMCPJetsPerMCCollisionPreslice = aod::jet::mcCollisionId;
  // Preslice<DsMCPJetsEWS> dsMCPJetsEWSPerMCCollisionPreslice = aod::jet::mcCollisionId;

  // Event configurables
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> jetPtMin{"jetPtMin", 3.0, "minimum jet pT cut"};
  // Configurable<float> jetR{"jetR", 0.4, "jet resolution parameter"};

  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

  // Event-wise constituent subtraction jet tables
  // Configurable<float> centralityMin{"centralityMin", -999.f, "Minimum FT0M centrality"};
  // Configurable<float> centralityMax{"centralityMax", 999.f, "Maximum FT0M centrality"};

  // internals
  std::vector<int> eventSelectionBits;
  int trackSelection = -1;

  // Filters
  // Filter jetCuts = aod::jet::pt > jetPtMin&& aod::jet::r == nround(jetR.node() * 100.0f);
  Filter jetCuts = aod::jet::pt > jetPtMin;

  // Filtered jet tables
  using FilteredDsDataJets = soa::Filtered<DsDataJets>;
  using FilteredDsMCDJets = soa::Filtered<DsMCDJets>;
  using FilteredDsMCPJets = soa::Filtered<DsMCPJets>;

  using FilteredDsDataJetsEWS = soa::Filtered<DsDataJetsEWS>;
  // using FilteredDsMCDJetsEWS = soa::Filtered<DsMCDJetsEWS>;
  // using FilteredDsMCPJetsEWS = soa::Filtered<DsMCPJetsEWS>;

  // Filtered inclusive charged jets
  using FilteredChargedJets = soa::Filtered<ChargedJets>;
  using FilteredChargedJetsEWS = soa::Filtered<ChargedJetsEWS>;

  using FilteredDsMCPJetsOnTheFly = soa::Filtered<DsMCPJetsOnTheFly>;

  //=====================================================================================
  // Histogram definitions
  //=====================================================================================

  HistogramRegistry registry;

  void addInclusiveQAHistograms()
  {
    registry.add("h_collisions", "event status;event status;entries", {HistType::kTH1F, {{10, 0., 10.}}});

    // Track QA
    registry.add("h_track_pt", ";#it{p}_{T,track};entries", {HistType::kTH1F, {{200, 0., 200.}}});
    registry.add("h_track_eta", ";#eta_{track};entries", {HistType::kTH1F, {{100, -1., 1.}}});
    registry.add("h_track_phi", ";#varphi_{track};entries", {HistType::kTH1F, {{80, -1., 7.}}});

    // Inclusive charged-jet QA
    registry.add("h_njets_inclusive", ";Number of charged jets per event;Events", {HistType::kTH1F, {{100, 0., 100.}}});
    registry.add("h_jet_pt_inclusive", ";#it{p}_{T,jet} (GeV/#it{c});Entries", {HistType::kTH1F, {{200, 0., 200.}}});
    registry.add("h_jet_eta_inclusive", ";#eta_{jet};Entries", {HistType::kTH1F, {{100, -1., 1.}}});
    registry.add("h_jet_phi_inclusive", ";#varphi_{jet} (rad);Entries", {HistType::kTH1F, {{80, -1., 7.}}});
    registry.add("h_jet_mass_inclusive", ";#it{m}_{jet} (GeV/#it{c}^{2});Entries", {HistType::kTH1F, {{120, 0., 60.}}});
    registry.add("h_jet_nconst_inclusive", ";Jet constituents;Entries", {HistType::kTH1F, {{100, 0., 100.}}});
    registry.add("h_lambda11_inclusive", ";#lambda_{1}^{1};Entries", {HistType::kTH1F, {{100, 0., 1.}}});
    registry.add("h_lambda21_inclusive", ";#lambda_{2}^{1};Entries", {HistType::kTH1F, {{100, 0., 1.}}});
  }

  void addDataHistograms()
  {
    registry.add("h_event_counter_data", ";Selection step;Events", {HistType::kTH1F, {{3, 0.5, 3.5}}});

    // Jet QA
    registry.add("h_jet_pt_data", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}});
    registry.add("h_jet_eta_data", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1., 1.}}});
    registry.add("h_jet_phi_data", "jet #phi;#varphi_{jet};entries", {HistType::kTH1F, {{80, -1., 7.}}});

    // Ds QA
    registry.add("h_ds_mass_data", ";m_{D_{S}} (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{350, 1.6, 2.3}}});
    registry.add("h_ds_pt_data", ";#it{p}_{T,D_{S}} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 50.}}});
    registry.add("h_ds_eta_data", ";#eta_{D_{S}};entries", {HistType::kTH1F, {{100, -1., 1.}}});
    registry.add("h_ds_phi_data", ";#varphi_{D_{S}};entries", {HistType::kTH1F, {{80, -1., 7.}}});

    // Ds-tagged jet observables
    registry.add("h_ds_jet_projection_data", ";z^{D_{S},jet}_{||};entries", {HistType::kTH1F, {{200, 0., 1.}}});
    registry.add("h_ds_jet_distance_data", ";#DeltaR_{D_{S},jet};entries", {HistType::kTH1F, {{200, 0., 1.}}});
    registry.add("h_ds_jet_mass_data", ";m_{jet}^{ch} (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{300, 0., 25.}}});
    registry.add("h_ds_jet_lambda11_data", ";#lambda_{1}^{1};entries", {HistType::kTH1F, {{200, 0., 1.}}});
    registry.add("h_ds_jet_lambda21_data", ";#lambda_{2}^{1};entries", {HistType::kTH1F, {{200, 0., 1.}}});

    // Main data sparse
    registry.add("hSparse_ds_data", ";m_{D_{S}};#it{p}_{T,D_{S}};#it{p}_{T,jet};z^{D_{S},jet}_{||};#DeltaR_{D_{S},jet}", {HistType::kTHnSparseF, {{350, 1.6, 2.3}, {100, 0., 50.}, {100, 0., 100.}, {200, 0., 1.}, {200, 0., 1.}}});

    auto hEventCounter = registry.get<TH1>(HIST("h_event_counter_data"));
    hEventCounter->GetXaxis()->SetBinLabel(1, "Input collisions");
    hEventCounter->GetXaxis()->SetBinLabel(2, "Event selection");
    hEventCounter->GetXaxis()->SetBinLabel(3, "|z| < 10 cm");
  }

  void addMCEfficiencyHistograms()
  {
    // General MC counters
    registry.add("hMCJetCounter", "N_{jet};", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hMCColCounter", "N_{collisions};", {HistType::kTH1F, {{4, 0., 4.}}});
    // Detector-level: all MCD Ds-tagged jets
    // Detector-level jet QA
    registry.add("h_jet_pt_mcd", "detector-level jet pT;#it{p}_{T,jet}^{det} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 100.}}});
    registry.add("h_jet_eta_mcd", "detector-level jet #eta;#eta_{jet}^{det};entries", {HistType::kTH1F, {{100, -1., 1.}}});
    registry.add("h_jet_phi_mcd", "detector-level jet #varphi;#varphi_{jet}^{det};entries", {HistType::kTH1F, {{80, -1., 7.}}});
    // Detector-level Ds QA
    registry.add("h_ds_pt_mcd", ";#it{p}_{T,D_{S}}^{det} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 50.}}});
    registry.add("h_ds_eta_mcd", ";#eta_{D_{S}}^{det};entries", {HistType::kTH1F, {{60, -1., 1.}}});
    registry.add("h_ds_phi_mcd", ";#varphi_{D_{S}}^{det};entries", {HistType::kTH1F, {{80, -1., 7.}}});
    registry.add("h_ds_mass_mcd", ";m_{D_{S}}^{det} (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{350, 1.6, 2.3}}});
    // Detector-level sparse histograms
    registry.add(
      "hSparse_ds_mcd1",
      ";m_{D_{S}}^{rec};#it{p}_{T,D_{S}}^{det};#it{p}_{T,jet}^{det};z^{D_{S},jet}_{||,det};Origin(D_{S});Matching status",
      {HistType::kTHnSparseF,
       {{350, 1.6, 2.3},
        {100, 0., 50.},
        {100, 0., 100.},
        {200, 0., 1.},
        {3, -0.5, 2.5},
        {2, -0.5, 1.5}}});

    registry.add(
      "hSparse_ds_mcd2",
      ";#it{p}_{T,D_{S}}^{det};#it{p}_{T,jet}^{det};#DeltaR_{D_{S},jet}^{det}",
      {HistType::kTHnSparseF,
       {{100, 0., 50.},
        {100, 0., 100.},
        {200, 0., 1.}}});

    registry.add(
      "hSparse_ds_mcd3",
      ";#it{p}_{T,jet}^{det};z^{D_{S},jet}_{||,det};#DeltaR_{D_{S},jet}^{det}",
      {HistType::kTHnSparseF,
       {{100, 0., 100.},
        {200, 0., 1.},
        {200, 0., 1.}}});

    // Particle-level: all MCP Ds-tagged jets
    // Particle-level jet QA
    registry.add("h_jet_pt_mcp", "particle-level jet pT;#it{p}_{T,jet}^{part} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 100.}}});
    registry.add("h_jet_eta_mcp", "particle-level jet #eta;#eta_{jet}^{part};entries", {HistType::kTH1F, {{100, -1., 1.}}});
    registry.add("h_jet_phi_mcp", "particle-level jet #varphi;#varphi_{jet}^{part};entries", {HistType::kTH1F, {{80, -1., 7.}}});
    // Particle-level Ds QA
    registry.add("h_ds_pt_mcp", ";#it{p}_{T,D_{S}}^{part} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 50.}}});
    registry.add("h_ds_eta_mcp", ";#eta_{D_{S}}^{part};entries", {HistType::kTH1F, {{60, -1., 1.}}});
    registry.add("h_ds_phi_mcp", ";#varphi_{D_{S}}^{part};entries", {HistType::kTH1F, {{80, -1., 7.}}});

    // Particle-level sparse with origin and matching status
    registry.add(
      "hSparse_ds_mcp",
      ";#it{p}_{T,D_{S}}^{part};#it{p}_{T,jet}^{part};z^{D_{S},jet}_{||,part};#DeltaR_{D_{S},jet}^{part};Origin(D_{S});Matching status",
      {HistType::kTHnSparseF,
       {{100, 0., 50.},
        {100, 0., 100.},
        {200, 0., 1.},
        {200, 0., 1.},
        {3, -0.5, 2.5},
        {2, -0.5, 1.5}}});

    // Matched detector-level Ds-tagged jets
    // Matched detector-level jet QA
    registry.add("h_jet_pt_mcd_matched", "matched detector-level jet pT;#it{p}_{T,jet}^{det} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 100.}}});
    registry.add("h_jet_eta_mcd_matched", "matched detector-level jet #eta;#eta_{jet}^{det};entries", {HistType::kTH1F, {{100, -1., 1.}}});
    registry.add("h_jet_phi_mcd_matched", "matched detector-level jet #varphi;#varphi_{jet}^{det};entries", {HistType::kTH1F, {{80, -1., 7.}}});
    // Matched detector-level Ds QA
    registry.add("h_ds_pt_mcd_matched", ";#it{p}_{T,D_{S}}^{det} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 50.}}});
    registry.add("h_ds_eta_mcd_matched", ";#eta_{D_{S}}^{det};entries", {HistType::kTH1F, {{60, -1., 1.}}});
    registry.add("h_ds_phi_mcd_matched", ";#varphi_{D_{S}}^{det};entries", {HistType::kTH1F, {{80, -1., 7.}}});
    registry.add("h_ds_mass_mcd_matched", ";m_{D_{S}}^{det} (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{350, 1.6, 2.3}}});
    // Matched detector-level sparse histograms
    registry.add(
      "hSparse_ds_mcd_matched1",
      ";m_{D_{S}}^{rec};#it{p}_{T,D_{S}}^{det};#it{p}_{T,jet}^{det};z^{D_{S},jet}_{||,det};Origin(D_{S})",
      {HistType::kTHnSparseF,
       {{350, 1.6, 2.3},
        {100, 0., 50.},
        {100, 0., 100.},
        {200, 0., 1.},
        {3, -0.5, 2.5}}});

    registry.add(
      "hSparse_ds_mcd_matched2",
      ";#it{p}_{T,D_{S}}^{det};#it{p}_{T,jet}^{det};#DeltaR_{D_{S},jet}^{det}",
      {HistType::kTHnSparseF,
       {{100, 0., 50.},
        {100, 0., 100.},
        {200, 0., 1.}}});

    registry.add(
      "hSparse_ds_mcd_matched3",
      ";#it{p}_{T,jet}^{det};z^{D_{S},jet}_{||,det};#DeltaR_{D_{S},jet}^{det}",
      {HistType::kTHnSparseF,
       {{100, 0., 100.},
        {200, 0., 1.},
        {200, 0., 1.}}});

    // MCP <-> MCD matched response
    registry.add(
      "hSparseMatchedJets",
      "Matched Ds-tagged jets;#it{p}_{T,jet}^{part};#it{p}_{T,jet}^{det};#it{p}_{T,D_{S}}^{part};#it{p}_{T,D_{S}}^{det};z_{||}^{part};z_{||}^{det};Origin(D_{S});flagMcMatchRec",
      {HistType::kTHnSparseF,
       {{100, 0., 100.},
        {100, 0., 100.},
        {100, 0., 50.},
        {100, 0., 50.},
        {200, 0., 1.},
        {200, 0., 1.},
        {3, -0.5, 2.5},
        {21, -10.5, 10.5}}});

    // Counter labels
    auto mcCollisionCounter = registry.get<TH1>(HIST("hMCColCounter"));
    mcCollisionCounter->GetXaxis()->SetBinLabel(BinMCColCntr::AllMCCollisions, "All MC coll.");
    mcCollisionCounter->GetXaxis()->SetBinLabel(BinMCColCntr::SelectedMCCollisions, "Selected MC coll.");
    mcCollisionCounter->GetXaxis()->SetBinLabel(BinMCColCntr::AssociatedRecoCollisions, "MC-associated reco coll.");
    mcCollisionCounter->GetXaxis()->SetBinLabel(BinMCColCntr::SelectedAssociatedRecoCollisions, "Selected MC-associated reco coll.");

    auto jetCounter = registry.get<TH1>(HIST("hMCJetCounter"));
    jetCounter->GetXaxis()->SetBinLabel(BinMCJetCntr::MCPJets, "MCP Ds-jets");
    jetCounter->GetXaxis()->SetBinLabel(BinMCJetCntr::MCPJetsWithMCDMatch, "MCP jets w/ MCD match");
    jetCounter->GetXaxis()->SetBinLabel(BinMCJetCntr::MCDJets, "MCD Ds-jets");
    jetCounter->GetXaxis()->SetBinLabel(BinMCJetCntr::MCDJetsWithMCPMatch, "MCD jets w/ MCP match");
    jetCounter->GetXaxis()->SetBinLabel(BinMCJetCntr::MatchedPairs, "MCP-MCD pairs");

    // Matching-status labels
    auto hSparseMCD = registry.get<THnSparse>(HIST("hSparse_ds_mcd1"));
    hSparseMCD->GetAxis(5)->SetBinLabel(1, "Unmatched");
    hSparseMCD->GetAxis(5)->SetBinLabel(2, "Matched");

    auto hSparseMCP = registry.get<THnSparse>(HIST("hSparse_ds_mcp"));
    hSparseMCP->GetAxis(5)->SetBinLabel(1, "Unmatched");
    hSparseMCP->GetAxis(5)->SetBinLabel(2, "Matched");
  }

  void addMCPOnTheFlyHistograms()
  {
    registry.add("h_event_counter_mcp_on_the_fly", ";Selection step;MC collisions", {HistType::kTH1F, {{2, 0.5, 2.5}}});

    // Particle-level jet QA
    registry.add("h_jet_pt_mcp_on_the_fly", ";#it{p}_{T,jet}^{part} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 100.}}});
    registry.add("h_jet_eta_mcp_on_the_fly", ";#eta_{jet}^{part};entries", {HistType::kTH1F, {{100, -1., 1.}}});
    registry.add("h_jet_phi_mcp_on_the_fly", ";#varphi_{jet}^{part};entries", {HistType::kTH1F, {{80, -1., 7.}}});

    // Particle-level Ds QA
    registry.add("h_ds_pt_mcp_on_the_fly", ";#it{p}_{T,D_{S}}^{part} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 50.}}});
    registry.add("h_ds_eta_mcp_on_the_fly", ";#eta_{D_{S}}^{part};entries", {HistType::kTH1F, {{60, -1., 1.}}});
    registry.add("h_ds_phi_mcp_on_the_fly", ";#varphi_{D_{S}}^{part};entries", {HistType::kTH1F, {{80, -1., 7.}}});

    // Theoretical particle-level sparse
    registry.add("hSparse_ds_mcp_on_the_fly", ";#it{p}_{T,D_{S}}^{part};#it{p}_{T,jet}^{part};z^{D_{S},jet}_{||,part};#DeltaR_{D_{S},jet}^{part};Origin(D_{S})", {HistType::kTHnSparseF, {{100, 0., 50.}, {100, 0., 100.}, {200, 0., 1.}, {200, 0., 1.}, {2, -0.5, 1.5}}});

    auto hCounter = registry.get<TH1>(HIST("h_event_counter_mcp_on_the_fly"));
    hCounter->GetXaxis()->SetBinLabel(1, "Generated collisions");
    hCounter->GetXaxis()->SetBinLabel(2, "|z_{vtx}^{gen}| < cut");

    auto hSparse = registry.get<THnSparse>(HIST("hSparse_ds_mcp_on_the_fly"));
    hSparse->GetAxis(4)->SetBinLabel(1, "Prompt");
    hSparse->GetAxis(4)->SetBinLabel(2, "Non-prompt");
  }

  //========
  // INIT
  //========

  void init(InitContext const&)
  {
    // Initialise event and track selection criteria
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

    // Prevent simultaneous execution of processes using alternative jet collections
    if (doprocessCollisions && doprocessCollisionsEWS) {
      throw std::runtime_error("Enable either processCollisions or processCollisionsEWS, not both simultaneously");
    }

    if (doprocessDataChargedSubstructure && doprocessDataChargedSubstructureEWS) {
      throw std::runtime_error("Enable either processDataChargedSubstructure or processDataChargedSubstructureEWS, not both simultaneously");
    }
    // Determine which histogram groups are required
    const bool doInclusiveQA = doprocessCollisions || doprocessCollisionsEWS;
    const bool doData = doprocessDataChargedSubstructure || doprocessDataChargedSubstructureEWS;
    const bool doMCEfficiency = doprocessMonteCarloEfficiencyDs;
    const bool doMCPOnTheFly = doprocessMCPOnTheFly;

    // Register histograms only if the corresponding process is enabled
    if (doInclusiveQA) {
      addInclusiveQAHistograms();
    }

    if (doData) {
      addDataHistograms();
    }

    if (doMCEfficiency) {
      addMCEfficiencyHistograms();
    }

    if (doMCPOnTheFly) {
      addMCPOnTheFlyHistograms();
    }
  }
  //===============
  // Lambda compute
  //===============

  template <typename JET, typename TRACKS>
  float computeLambda(JET const& jet, TRACKS const& tracks, float alpha, float kappa)
  {
    if (jet.pt() <= 0.f) {
      return -1.f;
    }

    float sum = 0.f;
    for (auto const& trk : tracks) {
      const float dr = jetutilities::deltaR(jet, trk);
      sum += std::pow(trk.pt(), kappa) * std::pow(dr, alpha);
    }
    const float jetRadius = jet.r() / 100.f;

    const float denom = std::pow(jet.pt(), kappa) * std::pow(jetRadius, alpha);
    if (denom <= 0.f) {
      return -1.f;
    }
    return sum / denom;
  }

  //=================
  // Jet Mass compute
  //=================

  template <typename TRACKS>
  float computeJetMass(TRACKS const& tracks)
  {
    double sumPx = 0.0, sumPy = 0.0, sumPz = 0.0, sumE = 0.0;

    for (auto const& trk : tracks) {
      const double pt = trk.pt();
      const double phi = trk.phi();
      const double eta = trk.eta();

      const double px = pt * std::cos(phi);
      const double py = pt * std::sin(phi);
      const double pz = pt * std::sinh(eta);
      const double p = std::sqrt(px * px + py * py + pz * pz);

      sumPx += px;
      sumPy += py;
      sumPz += pz;
      sumE += p; // massless
    }

    const double m2 = sumE * sumE - (sumPx * sumPx + sumPy * sumPy + sumPz * sumPz);
    return (m2 > 0.0) ? static_cast<float>(std::sqrt(m2)) : 0.f;
  }

  //==============
  // Collision QA
  //==============

  // This function is shared by the pp and Pb–Pb (event-wise subtraction)
  template <typename JetTable>
  void fillInclusiveJetQA(aod::JetCollision const& collision,
                          aod::JetTracks const& tracks,
                          JetTable const& jets)
  {
    // Collision counter
    registry.fill(HIST("h_collisions"), getValFromBin(AllCollisions));

    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    registry.fill(HIST("h_collisions"), getValFromBin(Sel8ZCut));

    // Track QA
    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }

      registry.fill(HIST("h_track_pt"), track.pt());
      registry.fill(HIST("h_track_eta"), track.eta());
      registry.fill(HIST("h_track_phi"), track.phi());
    }

    // Inclusive charged-jet QA
    registry.fill(HIST("h_njets_inclusive"), jets.size());

    for (auto const& jet : jets) {

      registry.fill(HIST("h_jet_pt_inclusive"), jet.pt());
      registry.fill(HIST("h_jet_eta_inclusive"), jet.eta());
      registry.fill(HIST("h_jet_phi_inclusive"), jet.phi());

      auto constituents = jet.template tracks_as<aod::JetTracks>();

      registry.fill(HIST("h_jet_nconst_inclusive"), constituents.size());

      registry.fill(HIST("h_jet_mass_inclusive"), computeJetMass(constituents));
      registry.fill(HIST("h_lambda11_inclusive"), computeLambda(jet, constituents, 1.f, 1.f));
      registry.fill(HIST("h_lambda21_inclusive"), computeLambda(jet, constituents, 2.f, 1.f));
    }
  }
  // Inclusive charged-jet QA for pp collisions
  void processCollisions(aod::JetCollision const& collision,
                         aod::JetTracks const& tracks,
                         FilteredChargedJets const& jets)
  {
    fillInclusiveJetQA(collision, tracks, jets);
  }
  PROCESS_SWITCH(JetDsSpecSubs, processCollisions, "collision QA", false);
  // Inclusive charged-jet QA for Pb–Pb collisions using
  // event-wise constituent-subtracted jets
  void processCollisionsEWS(aod::JetCollision const& collision,
                            aod::JetTracks const& tracks,
                            FilteredChargedJetsEWS const& jets)
  {
    fillInclusiveJetQA(collision, tracks, jets);
  }
  PROCESS_SWITCH(JetDsSpecSubs, processCollisionsEWS, "collision QA EWS", false);

  //=====================================================================================
  // DATA function
  //=====================================================================================

  // Common implementation of the Ds-in-jet analysis for data
  // Shared by the pp and Pb–Pb (event-wise subtraction) workflows
  template <typename JetTable>
  void analyseDataDsJet(aod::JetCollision const& collision,
                        JetTable const& jets)
  {
    registry.fill(HIST("h_event_counter_data"), 1);

    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    registry.fill(HIST("h_event_counter_data"), 2);

    if (std::abs(collision.posZ()) >= vertexZCut) {
      return;
    }
    registry.fill(HIST("h_event_counter_data"), 3);

    for (const auto& jet : jets) {

      registry.fill(HIST("h_jet_pt_data"), jet.pt());
      registry.fill(HIST("h_jet_eta_data"), jet.eta());
      registry.fill(HIST("h_jet_phi_data"), jet.phi());

      auto jetTracks = jet.template tracks_as<aod::JetTracks>();

      const float lambda11 = computeLambda(jet, jetTracks, 1.f, 1.f);
      const float lambda21 = computeLambda(jet, jetTracks, 2.f, 1.f);

      const float mjet = computeJetMass(jetTracks);

      const TVector3 jetVector(jet.px(), jet.py(), jet.pz());

      // Loop over Ds candidates
      for (const auto& dsCandidate : jet.template candidates_as<aod::CandidatesDsData>()) {

        const TVector3 dsVector(dsCandidate.px(), dsCandidate.py(), dsCandidate.pz());

        // Axis distance Delta_R
        const float deltaR = jetutilities::deltaR(jet, dsCandidate);
        // zParallel defined as longitudinal momentum fraction along the jet axis
        const float zParallel = (jetVector * dsVector) / (jetVector * jetVector);

        // --- Ds-level observables ---
        registry.fill(HIST("h_ds_mass_data"), dsCandidate.m());
        registry.fill(HIST("h_ds_pt_data"), dsCandidate.pt());
        registry.fill(HIST("h_ds_eta_data"), dsCandidate.eta());
        registry.fill(HIST("h_ds_phi_data"), dsCandidate.phi());

        registry.fill(HIST("h_ds_jet_projection_data"), zParallel);
        registry.fill(HIST("h_ds_jet_distance_data"), deltaR);

        // Main THnSparse: invariant mass, pT, z, and DeltaR
        registry.fill(HIST("hSparse_ds_data"),
                      dsCandidate.m(),
                      dsCandidate.pt(),
                      jet.pt(),
                      zParallel,
                      deltaR);
      }

      if (!jet.template candidates_as<aod::CandidatesDsData>().empty()) {
        // Jet Mass
        registry.fill(HIST("h_ds_jet_mass_data"), mjet);

        // Jet angularity
        if (lambda11 >= 0.f) {
          registry.fill(HIST("h_ds_jet_lambda11_data"), lambda11);
        }
        if (lambda21 >= 0.f) {
          registry.fill(HIST("h_ds_jet_lambda21_data"), lambda21);
        }
      }
    }
  }

  //=====================================================================================
  // DATA process
  //=====================================================================================

  // Data analysis using standard charged jets (pp)
  void processDataChargedSubstructure(aod::JetCollision const& collision,
                                      FilteredDsDataJets const& jets,
                                      aod::CandidatesDsData const&,
                                      aod::JetTracks const&)
  {
    analyseDataDsJet(collision, jets);
  }
  PROCESS_SWITCH(JetDsSpecSubs, processDataChargedSubstructure, "Data charged jets", false);

  // Data analysis using event-wise constituent-subtracted charged jets (Pb–Pb)
  void processDataChargedSubstructureEWS(aod::JetCollision const& collision,
                                         FilteredDsDataJetsEWS const& jets,
                                         aod::CandidatesDsData const&,
                                         aod::JetTracks const&)
  {
    analyseDataDsJet(collision, jets);
  }
  PROCESS_SWITCH(JetDsSpecSubs, processDataChargedSubstructureEWS, "Data charged jets EWS", false);

  //==================
  // MC function
  //==================

  template <typename TMCPJetsPerMCCollisionPreslice,
            typename TJetsMCD,
            typename TJetsMCP,
            typename TCandidatesMCD,
            typename TCandidatesMCP>
  void analyseMonteCarloEfficiency(
    TMCPJetsPerMCCollisionPreslice const& MCPJetsPerMCCollisionPreslice,
    aod::JetMcCollisions const& mccollisions,
    aod::JetCollisionsMCD const& collisions,
    TJetsMCD const& mcdjets,
    TJetsMCP const& mcpjets,
    TCandidatesMCD const&,
    TCandidatesMCP const&,
    aod::JetTracks const&,
    aod::JetParticles const&)
  {
    // MC collisions and particle-level Ds-tagged jets
    for (const auto& mccollision : mccollisions) {
      // All MC collisions
      registry.fill(HIST("hMCColCounter"), getValFromBin(BinMCColCntr::AllMCCollisions));

      // MC collision selection
      if (!jetderiveddatautilities::selectCollision(mccollision, eventSelectionBits) ||
          !(std::abs(mccollision.posZ()) < vertexZCut)) {
        continue;
      }
      // Selected MC collisions
      registry.fill(HIST("hMCColCounter"), getValFromBin(BinMCColCntr::SelectedMCCollisions));

      // Particle-level jets belonging to the current MC collision
      const auto mcpJetsPerMCCollision = mcpjets.sliceBy(MCPJetsPerMCCollisionPreslice, mccollision.globalIndex());

      // Particle-level Ds-tagged jets
      for (const auto& mcpjet : mcpJetsPerMCCollision) {

        registry.fill(HIST("hMCJetCounter"), getValFromBin(BinMCJetCntr::MCPJets));

        // Obtain the Ds particle associated with the particle-level jet
        auto mcpcand = mcpjet.template candidates_first_as<TCandidatesMCP>();

        // Particle-level Ds-jet observables
        TVector3 jetVecMCP;
        jetVecMCP.SetPtEtaPhi(mcpjet.pt(), mcpjet.eta(), mcpjet.phi());

        TVector3 dsVecMCP;
        dsVecMCP.SetPtEtaPhi(mcpcand.pt(), mcpcand.eta(), mcpcand.phi());

        const float zParallelMCP = dsVecMCP.Dot(jetVecMCP) / jetVecMCP.Mag2();
        const float deltaRMCP = jetutilities::deltaR(mcpjet, mcpcand);

        const int originMCP = static_cast<int>(mcpcand.originMcGen());
        const int flagMcMatchGen = static_cast<int>(mcpcand.flagMcMatchGen());
        const bool isMatchedMCP = mcpjet.has_matchedJetCand();

        // Particle-level jet QA
        registry.fill(HIST("h_jet_pt_mcp"), mcpjet.pt());
        registry.fill(HIST("h_jet_eta_mcp"), mcpjet.eta());
        registry.fill(HIST("h_jet_phi_mcp"), mcpjet.phi());

        // Particle-level Ds QA
        registry.fill(HIST("h_ds_pt_mcp"), mcpcand.pt());
        registry.fill(HIST("h_ds_eta_mcp"), mcpcand.eta());
        registry.fill(HIST("h_ds_phi_mcp"), mcpcand.phi());

        // Particle-level sparse
        registry.fill(
          HIST("hSparse_ds_mcp"),
          mcpcand.pt(),
          mcpjet.pt(),
          zParallelMCP,
          deltaRMCP,
          originMCP,
          static_cast<float>(isMatchedMCP));

        // Store all particle-level Ds-tagged jets
        outputMCPJets(
          mcpjet.pt(),
          mcpjet.eta(),
          mcpjet.phi(),
          mcpjet.template tracks_as<aod::JetParticles>().size() +
            mcpjet.template candidates_as<TCandidatesMCP>().size(),
          mcpcand.pt(),
          mcpcand.eta(),
          mcpcand.phi(),
          zParallelMCP,
          deltaRMCP,
          originMCP,
          flagMcMatchGen,
          isMatchedMCP);

        // MCP -> MCD matched pairs
        if (isMatchedMCP) {

          registry.fill(HIST("hMCJetCounter"), getValFromBin(BinMCJetCntr::MCPJetsWithMCDMatch));

          for (const auto& mcdjet : mcpjet.template matchedJetCand_as<TJetsMCD>()) {

            const auto& collision = collisions.iteratorAt(mcdjet.collisionId());

            // Detector-level collision selection
            if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits) ||
                !(std::abs(collision.posZ()) < vertexZCut)) {
              continue;
            }

            // Obtain the reconstructed Ds candidate
            auto mcdcand = mcdjet.template candidates_first_as<TCandidatesMCD>();

            // Detector-level Ds-jet observables
            TVector3 jetVecMCD;
            jetVecMCD.SetPtEtaPhi(mcdjet.pt(), mcdjet.eta(), mcdjet.phi());

            TVector3 dsVecMCD;
            dsVecMCD.SetPtEtaPhi(mcdcand.pt(), mcdcand.eta(), mcdcand.phi());

            const float zParallelMCD = dsVecMCD.Dot(jetVecMCD) / jetVecMCD.Mag2();
            const float deltaRMCD = jetutilities::deltaR(mcdjet, mcdcand);

            const int originMCD = static_cast<int>(mcdcand.originMcRec());
            const int flagMcMatchRec = static_cast<int>(mcdcand.flagMcMatchRec());

            // Matched detector-level jet QA
            registry.fill(HIST("h_jet_pt_mcd_matched"), mcdjet.pt());
            registry.fill(HIST("h_jet_eta_mcd_matched"), mcdjet.eta());
            registry.fill(HIST("h_jet_phi_mcd_matched"), mcdjet.phi());

            // Matched detector-level Ds QA
            registry.fill(HIST("h_ds_pt_mcd_matched"), mcdcand.pt());
            registry.fill(HIST("h_ds_eta_mcd_matched"), mcdcand.eta());
            registry.fill(HIST("h_ds_phi_mcd_matched"), mcdcand.phi());
            registry.fill(HIST("h_ds_mass_mcd_matched"), mcdcand.m());

            // Matched detector-level sparse: mass, pT Ds, pT jet, z_parallel, origin
            registry.fill(
              HIST("hSparse_ds_mcd_matched1"),
              mcdcand.m(),
              mcdcand.pt(),
              mcdjet.pt(),
              zParallelMCD,
              originMCD);

            // Matched detector-level sparse: pT Ds, pT jet, DeltaR
            registry.fill(
              HIST("hSparse_ds_mcd_matched2"),
              mcdcand.pt(),
              mcdjet.pt(),
              deltaRMCD);

            // Matched detector-level sparse: pT jet, z_parallel, DeltaR
            registry.fill(
              HIST("hSparse_ds_mcd_matched3"),
              mcdjet.pt(),
              zParallelMCD,
              deltaRMCD);

            // MCP <-> MCD response
            registry.fill(
              HIST("hSparseMatchedJets"),
              mcpjet.pt(),
              mcdjet.pt(),
              mcpcand.pt(),
              mcdcand.pt(),
              zParallelMCP,
              zParallelMCD,
              originMCP,
              flagMcMatchRec);

            // One accepted MCP-MCD pair
            registry.fill(HIST("hMCJetCounter"), getValFromBin(BinMCJetCntr::MatchedPairs));

            // Store matched MCP <-> MCD pair
            outputMatchedJets(
              // Particle level
              mcpjet.pt(),
              mcpjet.eta(),
              mcpjet.phi(),
              mcpjet.template tracks_as<aod::JetParticles>().size() +
                mcpjet.template candidates_as<TCandidatesMCP>().size(),
              mcpcand.pt(),
              mcpcand.eta(),
              mcpcand.phi(),
              zParallelMCP,
              deltaRMCP,

              // Detector level
              mcdjet.pt(),
              mcdjet.eta(),
              mcdjet.phi(),
              mcdjet.template tracks_as<aod::JetTracks>().size() +
                mcdjet.template candidates_as<TCandidatesMCD>().size(),
              mcdcand.pt(),
              mcdcand.eta(),
              mcdcand.phi(),
              mcdcand.m(),
              zParallelMCD,
              deltaRMCD,

              // MC truth
              originMCP);

          } // end of matched MCD jets
        } // end of isMatchedMCP
      } // end of MCP jets
    } // end of MC collisions

    // Detector-level Ds-tagged jets
    // This loop is outside the MC-collision loop so that each
    // detector-level jet is processed only once.
    for (const auto& mcdjet : mcdjets) {

      const auto& collision = collisions.iteratorAt(mcdjet.collisionId());

      // Detector-level collision selection
      if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits) ||
          !(std::abs(collision.posZ()) < vertexZCut)) {
        continue;
      }
      registry.fill(HIST("hMCJetCounter"), getValFromBin(BinMCJetCntr::MCDJets));

      // Obtain the reconstructed Ds candidate
      auto mcdcand = mcdjet.template candidates_first_as<TCandidatesMCD>();

      // Detector-level Ds-jet observables
      TVector3 jetVecMCD;
      jetVecMCD.SetPtEtaPhi(mcdjet.pt(), mcdjet.eta(), mcdjet.phi());

      TVector3 dsVecMCD;
      dsVecMCD.SetPtEtaPhi(mcdcand.pt(), mcdcand.eta(), mcdcand.phi());

      const float zParallelMCD = dsVecMCD.Dot(jetVecMCD) / jetVecMCD.Mag2();
      const float deltaRMCD = jetutilities::deltaR(mcdjet, mcdcand);
      const int originMCD = static_cast<int>(mcdcand.originMcRec());
      const int flagMcMatchRec = static_cast<int>(mcdcand.flagMcMatchRec());
      const bool isMatchedMCD = mcdjet.has_matchedJetCand();

      if (isMatchedMCD) {
        registry.fill(HIST("hMCJetCounter"), getValFromBin(BinMCJetCntr::MCDJetsWithMCPMatch));
      }

      // Detector-level jet QA
      registry.fill(HIST("h_jet_pt_mcd"), mcdjet.pt());
      registry.fill(HIST("h_jet_eta_mcd"), mcdjet.eta());
      registry.fill(HIST("h_jet_phi_mcd"), mcdjet.phi());

      // Detector-level Ds QA
      registry.fill(HIST("h_ds_pt_mcd"), mcdcand.pt());
      registry.fill(HIST("h_ds_eta_mcd"), mcdcand.eta());
      registry.fill(HIST("h_ds_phi_mcd"), mcdcand.phi());
      registry.fill(HIST("h_ds_mass_mcd"), mcdcand.m());

      // Detector-level sparse: mass, pT Ds, pT jet, z_parallel, origin, matching status
      registry.fill(
        HIST("hSparse_ds_mcd1"),
        mcdcand.m(),
        mcdcand.pt(),
        mcdjet.pt(),
        zParallelMCD,
        originMCD,
        static_cast<float>(isMatchedMCD));

      // Detector-level sparse: pT Ds, pT jet, DeltaR
      registry.fill(
        HIST("hSparse_ds_mcd2"),
        mcdcand.pt(),
        mcdjet.pt(),
        deltaRMCD);

      // Detector-level sparse: pT jet, z_parallel, DeltaR
      registry.fill(
        HIST("hSparse_ds_mcd3"),
        mcdjet.pt(),
        zParallelMCD,
        deltaRMCD);

      // Store all detector-level Ds-tagged jets
      outputMCDJets(
        mcdjet.pt(),
        mcdjet.eta(),
        mcdjet.phi(),
        mcdjet.template tracks_as<aod::JetTracks>().size() +
          mcdjet.template candidates_as<TCandidatesMCD>().size(),
        mcdcand.pt(),
        mcdcand.eta(),
        mcdcand.phi(),
        mcdcand.m(),
        zParallelMCD,
        deltaRMCD,
        originMCD,
        flagMcMatchRec,
        isMatchedMCD);

    } // end of MCD jets
  } // end of analyseMonteCarloEfficiency

  void processMonteCarloEfficiencyDs(aod::JetMcCollisions const& mccollisions,
                                     aod::JetCollisionsMCD const& collisions,
                                     DsMCDJets const& mcdjets,
                                     DsMCPJets const& mcpjets,
                                     DsCandidatesMCD const& mcdDscand,
                                     DsCandidatesMCP const& mcpDscand,
                                     aod::JetTracks const& tracks,
                                     aod::JetParticles const& particles)
  {
    analyseMonteCarloEfficiency<Preslice<DsMCPJets>,
                                DsMCDJets,
                                DsMCPJets,
                                DsCandidatesMCD,
                                DsCandidatesMCP>(
      dsMCPJetsPerMCCollisionPreslice,
      mccollisions,
      collisions,
      mcdjets,
      mcpjets,
      mcdDscand,
      mcpDscand,
      tracks,
      particles);
  }
  PROCESS_SWITCH(JetDsSpecSubs, processMonteCarloEfficiencyDs, "Non-matched and matched MC Ds and jets", false);

  //=====================================================================================
  // MC particle-level process for on-the-fly simulations
  //=====================================================================================

  void processMCPOnTheFly(aod::JetMcCollision const& mccollision,
                          FilteredDsMCPJetsOnTheFly const& mcpjets,
                          DsCandidatesMCP const&)
  {
    // Count all generated MC collisions before the vertex selection
    registry.fill(HIST("h_event_counter_mcp_on_the_fly"), 1);

    // Apply the generated-vertex selection
    if (std::abs(mccollision.posZ()) >= vertexZCut) {
      return;
    }

    // Count generated MC collisions passing the vertex selection
    registry.fill(HIST("h_event_counter_mcp_on_the_fly"), 2);

    // Loop over particle-level Ds-tagged charged jets in the current MC collision
    for (const auto& mcpjet : mcpjets) {

      // Retrieve the leading generated Ds candidate associated with the jet
      const auto mcpCandidate = mcpjet.template candidates_first_as<DsCandidatesMCP>();

      // Classify the generated Ds origin: 0 = prompt, 1 = non-prompt
      const auto origin = mcpCandidate.originMcGen();

      // Keep only prompt and non-prompt generated Ds candidates.
      if (origin != RecoDecay::OriginType::Prompt && origin != RecoDecay::OriginType::NonPrompt) {
        continue;
      }
      // Encode prompt and non-prompt origins as 0 and 1, respectively.
      const int originMCP = origin == RecoDecay::OriginType::Prompt ? 0 : 1;

      // Build the particle-level jet and Ds momentum vectors
      const TVector3 jetVector(mcpjet.px(), mcpjet.py(), mcpjet.pz());
      const TVector3 candidateVector(mcpCandidate.px(), mcpCandidate.py(), mcpCandidate.pz());

      // Compute the longitudinal momentum fraction of the Ds along the jet axis
      const float zParallel = jetVector.Dot(candidateVector) / jetVector.Mag2();

      // Compute the angular distance between the Ds candidate and the jet axis
      const float deltaR = jetutilities::deltaR(mcpjet, mcpCandidate);

      // Fill particle-level jet distributions
      registry.fill(HIST("h_jet_pt_mcp_on_the_fly"), mcpjet.pt());
      registry.fill(HIST("h_jet_eta_mcp_on_the_fly"), mcpjet.eta());
      registry.fill(HIST("h_jet_phi_mcp_on_the_fly"), mcpjet.phi());

      // Fill particle-level Ds-candidate distributions
      registry.fill(HIST("h_ds_pt_mcp_on_the_fly"), mcpCandidate.pt());
      registry.fill(HIST("h_ds_eta_mcp_on_the_fly"), mcpCandidate.eta());
      registry.fill(HIST("h_ds_phi_mcp_on_the_fly"), mcpCandidate.phi());

      // Store the particle-level Ds-tagged jet observables and Ds origin
      registry.fill(HIST("hSparse_ds_mcp_on_the_fly"),
                    mcpCandidate.pt(),
                    mcpjet.pt(),
                    zParallel,
                    deltaR,
                    originMCP);
    }
  }
  PROCESS_SWITCH(JetDsSpecSubs, processMCPOnTheFly, "Process on-the-fly MC particle-level Ds-tagged jets", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<JetDsSpecSubs>(cfgc)};
}
