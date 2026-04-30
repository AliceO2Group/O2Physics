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
///
/// \file JetD0AngSubstructure.cxx
/// \brief Analysis task for the reconstruction and study of charged jets
/// containing D_0 mesons in pp collisions. The code is partially inherited from
/// hfFragmentationFunction.cxx.
///
/// \author P. Dhankher
/// \author L.J. Huisman

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetHFUtilities.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
//
#include "PWGHF/Core/DecayChannels.h"

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
#include <TVector3.h>

#include <Rtypes.h>

#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>
//
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// Definition of a custom AOD table to store jet–D0 quantities
namespace o2::aod
{
/*
// Jet-related quantities
*/

DECLARE_SOA_COLUMN(ExpJetHfDist, expJetHfDist, float);
DECLARE_SOA_COLUMN(ExpJetPt, expJetPt, float);
DECLARE_SOA_COLUMN(ExpJetEta, expJetEta, float);
DECLARE_SOA_COLUMN(ExpJetPhi, expJetPhi, float);
DECLARE_SOA_COLUMN(ExpJetNConst, expJetNConst, int);
DECLARE_SOA_COLUMN(ExpJetAng, expJetAng, float);
// D0 candidate quantities
DECLARE_SOA_COLUMN(ExpHfZParallel, expHfZParallel, float);

DECLARE_SOA_COLUMN(ExpHfPt, expHfPt, float);
DECLARE_SOA_COLUMN(ExpHfEta, expHfEta, float);
DECLARE_SOA_COLUMN(ExpHfPhi, expHfPhi, float);
DECLARE_SOA_COLUMN(ExpHfMass, expHfMass, float);
DECLARE_SOA_COLUMN(ExpHfY, expHfY, float);
// ML scores
DECLARE_SOA_COLUMN(ExpHfMlScore0, expHfMlScore0, float);
DECLARE_SOA_COLUMN(ExpHfMlScore1, expHfMlScore1, float);
DECLARE_SOA_COLUMN(ExpHfMlScore2, expHfMlScore2, float);

/*
// MCP quantities (Particle Level)
*/
// Jets
DECLARE_SOA_COLUMN(McpJetHfDist, mcpJetHfDist, float);
DECLARE_SOA_COLUMN(McpJetPt, mcpJetPt, float);
DECLARE_SOA_COLUMN(McpJetEta, mcpJetEta, float);
DECLARE_SOA_COLUMN(McpJetPhi, mcpJetPhi, float);
DECLARE_SOA_COLUMN(McpJetNConst, mcpJetNConst, float);
DECLARE_SOA_COLUMN(McpJetAng, mcpJetAng, float);
// D0 candidates (Heavy Flavour)
DECLARE_SOA_COLUMN(McpHfZParallel, mcpHfZParallel, float);

DECLARE_SOA_COLUMN(McpHfPt, mcpHfPt, float);
DECLARE_SOA_COLUMN(McpHfEta, mcpHfEta, float);
DECLARE_SOA_COLUMN(McpHfPhi, mcpHfPhi, float);
// DECLARE_SOA_COLUMN(HfMass,mcpHfMass, float);
DECLARE_SOA_COLUMN(McpHfY, mcpHfY, float);
DECLARE_SOA_COLUMN(McpHfPrompt, mcpHfPrompt, bool);
DECLARE_SOA_COLUMN(McpHfMatch, mcpHfMatch, bool);

/*
// MCD quantities (Detector Level)
*/
// Jets
DECLARE_SOA_COLUMN(McdJetHfDist, mcdJetHfDist, float);
DECLARE_SOA_COLUMN(McdJetPt, mcdJetPt, float);
DECLARE_SOA_COLUMN(McdJetEta, mcdJetEta, float);
DECLARE_SOA_COLUMN(McdJetPhi, mcdJetPhi, float);
DECLARE_SOA_COLUMN(McdJetNConst, mcdJetNConst, float);
DECLARE_SOA_COLUMN(McdJetAng, mcdJetAng, float);
// D0 candidates (Heavy Flavour)
DECLARE_SOA_COLUMN(McdHfZParallel, mcdHfZParallel, float);

DECLARE_SOA_COLUMN(McdHfPt, mcdHfPt, float);
DECLARE_SOA_COLUMN(McdHfEta, mcdHfEta, float);
DECLARE_SOA_COLUMN(McdHfPhi, mcdHfPhi, float);
DECLARE_SOA_COLUMN(McdHfMass, mcdHfMass, float);
DECLARE_SOA_COLUMN(McdHfY, mcdHfY, float);
DECLARE_SOA_COLUMN(McdHfPrompt, mcdHfPrompt, bool);
DECLARE_SOA_COLUMN(McdHfMatch, mcdHfMatch, bool);
// Other
DECLARE_SOA_COLUMN(McdHfMatchedFrom, mcdHfMatchedFrom, int);
DECLARE_SOA_COLUMN(McdHfSelectedAs, mcdHfSelectedAs, int);
// ML scores
DECLARE_SOA_COLUMN(McdHfMlScore0, mcdHfMlScore0, float);
DECLARE_SOA_COLUMN(McdHfMlScore1, mcdHfMlScore1, float);
DECLARE_SOA_COLUMN(McdHfMlScore2, mcdHfMlScore2, float);

/*
// AOD table definition
*/
DECLARE_SOA_TABLE(EXPJetObjTable, "AOD", "EXPJETOBJTABLE",
                  ExpJetHfDist,
                  ExpJetPt,
                  ExpJetEta,
                  ExpJetPhi,
                  ExpJetNConst,
                  ExpJetAng,
                  ExpHfZParallel,
                  ExpHfPt,
                  ExpHfEta,
                  ExpHfPhi,
                  ExpHfMass,
                  ExpHfY,
                  ExpHfMlScore0,
                  ExpHfMlScore1,
                  ExpHfMlScore2);
DECLARE_SOA_TABLE(MCPJetObjTable, "AOD", "MCPJETOBJTABLE",
                  McpJetHfDist,
                  McpJetPt,
                  McpJetEta,
                  McpJetPhi,
                  McpJetNConst,
                  McpJetAng,
                  McpHfZParallel,
                  McpHfPt,
                  McpHfEta,
                  McpHfPhi,
                  McpHfY,
                  McpHfPrompt,
                  McpHfMatch);
DECLARE_SOA_TABLE(MCDJetObjTable, "AOD", "MCDJETOBJTABLE",
                  McdJetHfDist,
                  McdJetPt,
                  McdJetEta,
                  McdJetPhi,
                  McdJetNConst,
                  McdJetAng,
                  McdHfZParallel,
                  McdHfPt,
                  McdHfEta,
                  McdHfPhi,
                  McdHfMass,
                  McdHfY,
                  McdHfPrompt,
                  McdHfMatch,
                  McdHfMlScore0,
                  McdHfMlScore1,
                  McdHfMlScore2,
                  McdHfMatchedFrom,
                  McdHfSelectedAs);
DECLARE_SOA_TABLE(MatchJetDistanceTable, "AOD", "MATCHTABLE",
                  McpJetHfDist,
                  McpJetPt,
                  McpJetEta,
                  McpJetPhi,
                  McpJetNConst,
                  McpJetAng,
                  McpHfZParallel,
                  McpHfPt,
                  McpHfEta,
                  McpHfPhi,
                  McpHfY,
                  McpHfPrompt,
                  McdJetHfDist,
                  McdJetPt,
                  McdJetEta,
                  McdJetPhi,
                  McdJetNConst,
                  McdJetAng,
                  McdHfZParallel,
                  McdHfPt,
                  McdHfEta,
                  McdHfPhi,
                  McdHfMass,
                  McdHfY,
                  McdHfPrompt,
                  McdHfMlScore0,
                  McdHfMlScore1,
                  McdHfMlScore2,
                  McdHfMatchedFrom,
                  McdHfSelectedAs);

} // namespace o2::aod

// Helps to avoid typos in histogram names when using the "HIST" macro
/*
hf(l): Heavy flavour
ex(p): Experimental data (detector level)
mcd: Monte Carlo data (detector level)
mcp: Monte Carlo data (particle level)
*/
namespace histnames
{
#define HNAME(Name) constexpr const char* Name = #Name;
/*
// Experimental Data (analyseDataChargedSubstructure)
*/
HNAME(ex_col);           // Collision Counter
HNAME(ex_jet);           // Jet Counter
HNAME(ex_jet_pt);        // Jet pT
HNAME(ex_jet_eta);       // Jet eta
HNAME(ex_jet_phi);       // Jet phi
HNAME(ex_jet_ang);       // Jet angularity
HNAME(ex_jet_proj);      // Projection of HF candidate momentum on jet axis
HNAME(ex_jet_dist);      // Jet–HF candidate angular distance
HNAME(ex_jet_dist_proj); // Jet–HF candidate distance vs projection
HNAME(ex_hfl_pt);        // HF candidate pT
HNAME(ex_hfl_mass);      // HF candidate mass
HNAME(ex_hfl_eta);       // HF candidate eta
HNAME(ex_hfl_phi);       // HF candidate phi
/*
// Monte Carlo Data Efficiency (analyseMonteCarloEfficiency)
*/
HNAME(mc_eff_col);          // Collision Counter
HNAME(mc_eff_jet);          // Jet Counter
HNAME(mc_eff_det_jet_pt);   // Detector level jet pT
HNAME(mc_eff_det_jet_eta);  // Detector level jet eta
HNAME(mc_eff_det_jet_phi);  // Detector level jet phi
HNAME(mc_eff_det_jet_ang);  // Detector level jet angularity
HNAME(mc_eff_det_hfl_pt);   // Detector level HF candidate pT
HNAME(mc_eff_det_hfl_mass); // Detector level HF candidate mass
HNAME(mc_eff_det_hfl_eta);  // Detector level HF candidate eta
HNAME(mc_eff_det_hfl_phi);  // Detector level HF candidate phi
HNAME(mc_eff_par_jet_pt);   // Particle level jet pT
HNAME(mc_eff_par_jet_eta);  // Particle level jet eta
HNAME(mc_eff_par_jet_phi);  // Particle level jet phi
HNAME(mc_eff_par_jet_ang);  // Particle level jet angularity
HNAME(mc_eff_par_hfl_pt);   // Particle level HF candidate pT
// HNAME(mc_eff_par_hfl_mass); // Particle level HF candidate mass -> PDG value
HNAME(mc_eff_par_hfl_eta); // Particle level HF candidate eta
HNAME(mc_eff_par_hfl_phi); // Particle level HF candidate phi
/*
// Monte Carlo Data Matching (analyseMonteCarlo)
*/
HNAME(mc_col);          // Collision Counter
HNAME(mc_jet);          // Jet Counter
HNAME(mc_det_jet_pt);   // Detector level jet pT
HNAME(mc_det_jet_eta);  // Detector level jet eta
HNAME(mc_det_jet_phi);  // Detector level jet phi
HNAME(mc_det_jet_ang);  // Detector level jet angularity
HNAME(mc_det_hfl_pt);   // Detector level HF candidate pT (matched)
HNAME(mc_det_hfl_mass); // Detector level HF candidate mass (matched)
HNAME(mc_det_hfl_eta);  // Detector level HF candidate eta (matched)
HNAME(mc_det_hfl_phi);  // Detector level HF candidate phi (matched)
// HNAME(mc_par_jet_pt);   // Particle level jet pT
// HNAME(mc_par_jet_eta);  // Particle level jet eta
// HNAME(mc_par_jet_phi);  // Particle level jet phi
// HNAME(mc_par_jet_ang);  // Particle level jet angularity
// HNAME(mc_par_hfl_pt);   // Particle level HF candidate pT
// HNAME(mc_par_hfl_mass); // Particle level HF candidate mass -> PDG value
// HNAME(mc_par_hfl_eta); // Particle level HF candidate eta
// HNAME(mc_par_hfl_phi); // Particle level HF candidate phi

#undef HNAME
} // namespace histnames

consteval float getValFromBin(int bin)
{
  return static_cast<float>(bin) - 0.5f;
}

enum BIN_EX_COLCNTR { AllCollisions = 1,
                      Sel8ZCut = 2 };

enum BIN_EX_JETCNTR { ChargedJets = 1 };
enum BIN_MC_COLCNTR { All = 1,
                      ZCut = 2,
                      Matched = 3,
                      MatchedSel8ZCut = 4 };

enum BIN_MC_JETCNTR { DetectorLevelJetInMCCollision = 1,
                      ParticleLevelJetInMCCollision = 2,
                      DetectorLevelJetWithMatchedCandidate = 3,
                      ParticleLevelJetWithMatchedCandidate = 4
};

struct JetD0AngSubstructure {

  // Output table producer
  Produces<aod::EXPJetObjTable> objJetTable;
  Produces<aod::MCDJetObjTable> mcdJetTable;
  Produces<aod::MCPJetObjTable> mcpJetTable;
  Produces<aod::MatchJetDistanceTable> matchJetTable;
  //
  using JetChargedTableD0 = soa::Join<aod::D0ChargedJets, aod::D0ChargedJetConstituents>;
  // MC Matching Tables
  using JetD0MCDTable = soa::Join<aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents, aod::D0ChargedMCDetectorLevelJetsMatchedToD0ChargedMCParticleLevelJets>;
  using JetD0MCPTable = soa::Join<aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents, aod::D0ChargedMCParticleLevelJetsMatchedToD0ChargedMCDetectorLevelJets>;
  // Type Aliases for the constituent jets in the table (to be used for
  // angularity calculation)
  using JetD0MCDTableConstituent = JetD0MCDTable::iterator;
  using JetD0MCPTableConstituent = JetD0MCPTable::iterator;
  // Slices for access to proper HF MCD jet collision that is associated to
  // MCCollision
  PresliceUnsorted<aod::JetCollisionsMCD> collisionsPerMCCollisionPreslice = aod::jmccollisionlb::mcCollisionId;
  Preslice<JetD0MCDTable> d0MCDJetsPerEXPCollisionPreslice = aod::jet::collisionId;
  Preslice<JetD0MCPTable> d0MCPJetsPerMCCollisionPreslice = aod::jet::mcCollisionId;

  HistogramRegistry registry{"registry", {}};

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

  void init(o2::framework::InitContext&)
  {
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
    addHistograms();
  }

  void addHistograms()
  {
    /*
    // Experimental Data Histograms
    */
    // Definitions
    registry.add(histnames::ex_col, "N_{coll};", {HistType::kTH1F, {{2, 0., 2.}}});
    registry.add(histnames::ex_jet, "N_{jet};", {HistType::kTH1F, {{1, 0., 1.0}}});
    registry.add(histnames::ex_jet_pt, ";p_{T,jet};dN/dp_{T,jet}", {HistType::kTH1F, {{1000, 0., 50.}}});
    registry.add(histnames::ex_jet_eta, ";#eta_{jet};dN/d#eta_{jet}", {HistType::kTH1F, {{250, -5., 5.}}});
    registry.add(histnames::ex_jet_phi, ";#phi_{jet};dN/d#phi_{jet}", {HistType::kTH1F, {{250, -10., 10.}}});
    registry.add(histnames::ex_jet_ang, ";angularity_{jet};dN/d(angularity)", {HistType::kTH1F, {{1000, 0., 10.}}});
    registry.add(histnames::ex_jet_dist, ";#DeltaR_{D^{0},jet};dN/d(#DeltaR)", {HistType::kTH1F, {{1000, 0., 10.}}});
    registry.add(histnames::ex_jet_dist_proj, ";#DeltaR_{D^{0},jet};z^{D^{0},jet}_{||}", {HistType::kTH2F, {{1000, 0., 10.}, {1000, 0., 10.}}});
    registry.add(histnames::ex_jet_proj, ";z^{D^{0},jet}_{||};dN/dz^{D^{0},jet}_{||}", {HistType::kTH1F, {{1000, 0., 10.}}});
    registry.add(histnames::ex_hfl_pt, ";p_{T,D^{0}};dN/dp_{T,D^{0}}", {HistType::kTH1F, {{200, 0., 10.}}});
    registry.add(histnames::ex_hfl_mass, ";m_{D^{0}};dN/dm_{D^{0}}", {HistType::kTH1F, {{1000, 0., 10.}}});
    registry.add(histnames::ex_hfl_eta, ";#eta_{D^{0}};dN/d#eta_{D^{0}}", {HistType::kTH1F, {{250, -5., 5.}}});
    registry.add(histnames::ex_hfl_phi, ";#phi_{D^{0}};dN/d#phi_{D^{0}}", {HistType::kTH1F, {{250, -10., 10.}}});
    // Labels
    auto expCollisionCounter = registry.get<TH1>(HIST(histnames::ex_col));
    expCollisionCounter->GetXaxis()->SetBinLabel(BIN_EX_COLCNTR::AllCollisions, "all");
    expCollisionCounter->GetXaxis()->SetBinLabel(BIN_EX_COLCNTR::Sel8ZCut, "sel8 + zcut");

    auto expJetCounter = registry.get<TH1>(HIST(histnames::ex_jet));
    expJetCounter->GetXaxis()->SetBinLabel(BIN_EX_JETCNTR::ChargedJets, "Charged jets with D0");
    /*
    // Monte Carlo Data Efficiency Histograms
    */
    // Definitions
    registry.add(histnames::mc_eff_col, "N_{coll};", {HistType::kTH1F, {{4, 0., 4.0}}});
    registry.add(histnames::mc_eff_jet, "N_{jet};", {HistType::kTH1F, {{4, 0., 4.0}}});
    registry.add(histnames::mc_eff_det_jet_pt, ";p_{T,det jet};dN/dp_{T,det jet}", {HistType::kTH1F, {{1000, 0., 50.}}});
    registry.add(histnames::mc_eff_det_jet_eta, ";#eta_{det jet};dN/d#eta_{det jet}", {HistType::kTH1F, {{250, -5., 5.}}});
    registry.add(histnames::mc_eff_det_jet_phi, ";#phi_{det jet};dN/d#phi_{det jet}", {HistType::kTH1F, {{250, -10., 10.}}});
    registry.add(histnames::mc_eff_det_jet_ang, ";angularity_{det jet};dN/d(angularity)", {HistType::kTH1F, {{1000, 0., 10.}}});
    registry.add(histnames::mc_eff_det_hfl_pt, ";p_{T,det D^{0}};dN/dp_{T,det D^{0}}", {HistType::kTH1F, {{200, 0., 10.}}});
    registry.add(histnames::mc_eff_det_hfl_mass, ";m_{det D^{0}};dN/dm_{det D^{0}}", {HistType::kTH1F, {{1000, 0., 10.}}});
    registry.add(histnames::mc_eff_det_hfl_eta, ";#eta_{det D^{0}};dN/d#eta_{det D^{0}}", {HistType::kTH1F, {{250, -5., 5.}}});
    registry.add(histnames::mc_eff_det_hfl_phi, ";#phi_{det D^{0}};dN/d#phi_{det D^{0}}", {HistType::kTH1F, {{250, -10., 10.}}});
    // Labels
    auto mcCollisionCounter = registry.get<TH1>(HIST(histnames::mc_eff_col));
    mcCollisionCounter->GetXaxis()->SetBinLabel(BIN_MC_COLCNTR::All, "mccollisions");
    mcCollisionCounter->GetXaxis()->SetBinLabel(BIN_MC_COLCNTR::ZCut, "z_cut");
    mcCollisionCounter->GetXaxis()->SetBinLabel(BIN_MC_COLCNTR::Matched, "collisions");
    mcCollisionCounter->GetXaxis()->SetBinLabel(BIN_MC_COLCNTR::MatchedSel8ZCut, "sel8");

    auto jetCounter = registry.get<TH1>(HIST(histnames::mc_eff_jet));
    jetCounter->GetXaxis()->SetBinLabel(BIN_MC_JETCNTR::ParticleLevelJetInMCCollision, "particle level");
    jetCounter->GetXaxis()->SetBinLabel(BIN_MC_JETCNTR::DetectorLevelJetInMCCollision, "detector level");
    jetCounter->GetXaxis()->SetBinLabel(BIN_MC_JETCNTR::DetectorLevelJetWithMatchedCandidate, "particle matched jets");
    jetCounter->GetXaxis()->SetBinLabel(BIN_MC_JETCNTR::ParticleLevelJetWithMatchedCandidate, "detector matched jets");
    /*
    // Monte Carlo Data Histograms
    */
    registry.add(histnames::mc_col, "N_{coll};", {HistType::kTH1F, {{4, 0., 4.0}}});
    registry.add(histnames::mc_jet, "N_{jet};", {HistType::kTH1F, {{4, 0., 4.0}}});
    registry.add(histnames::mc_det_jet_pt, ";p_{T,det jet};dN/dp_{T,det jet}", {HistType::kTH1F, {{1000, 0., 50.}}});
    registry.add(histnames::mc_det_jet_eta, ";#eta_{det jet};dN/d#eta_{det jet}", {HistType::kTH1F, {{250, -5., 5.}}});
    registry.add(histnames::mc_det_jet_phi, ";#phi_{det jet};dN/d#phi_{det jet}", {HistType::kTH1F, {{250, -10., 10.}}});
    registry.add(histnames::mc_det_jet_ang, ";angularity_{det jet};dN/d(angularity)", {HistType::kTH1F, {{1000, 0., 10.}}});
    registry.add(histnames::mc_det_hfl_pt, ";p_{T,det D^{0}};dN/dp_{T,det D^{0}}", {HistType::kTH1F, {{200, 0., 10.}}});
    registry.add(histnames::mc_det_hfl_mass, ";m_{det D^{0}};dN/dm_{det D^{0}}", {HistType::kTH1F, {{1000, 0., 10.}}});
    registry.add(histnames::mc_det_hfl_eta, ";#eta_{det D^{0}};dN/d#eta_{det D^{0}}", {HistType::kTH1F, {{250, -5., 5.}}});
    registry.add(histnames::mc_det_hfl_phi, ";#phi_{det D^{0}};dN/d#phi_{det D^{0}}", {HistType::kTH1F, {{250, -10., 10.}}});
    // Labels
    auto mcCollisionCounter2 = registry.get<TH1>(HIST(histnames::mc_col));
    mcCollisionCounter2->GetXaxis()->SetBinLabel(BIN_MC_COLCNTR::All, "mccollisions");
    mcCollisionCounter2->GetXaxis()->SetBinLabel(BIN_MC_COLCNTR::ZCut, "z_cut");
    mcCollisionCounter2->GetXaxis()->SetBinLabel(BIN_MC_COLCNTR::Matched, "collisions");
    mcCollisionCounter2->GetXaxis()->SetBinLabel(BIN_MC_COLCNTR::MatchedSel8ZCut, "sel8");

    auto jetCounter2 = registry.get<TH1>(HIST(histnames::mc_jet));
    jetCounter2->GetXaxis()->SetBinLabel(BIN_MC_JETCNTR::DetectorLevelJetInMCCollision, "detector level");
    jetCounter2->GetXaxis()->SetBinLabel(BIN_MC_JETCNTR::ParticleLevelJetInMCCollision, "particle level");
    jetCounter2->GetXaxis()->SetBinLabel(BIN_MC_JETCNTR::DetectorLevelJetWithMatchedCandidate, "particle matched jets");
    jetCounter2->GetXaxis()->SetBinLabel(BIN_MC_JETCNTR::ParticleLevelJetWithMatchedCandidate, "detector matched jets");
  };

  template <typename T, typename U>
  float jetCalculateAngularityEXP(T const& jet, U const& /*tracks*/)
  {
    float tAngularity = 0.0;
    for (const auto& constituent : jet.template tracks_as<U>()) {
      tAngularity += std::pow(constituent.pt(), kappa) * std::pow(jetutilities::deltaR(jet, constituent) / (jet.r() / 100.f), alpha);
    }
    tAngularity /= std::pow(jet.pt(), kappa);
    return tAngularity;
  }

  template <typename JetTableMCDConstituent>
  float jetCalculateAngularityMCD(JetTableMCDConstituent const& jet, aod::JetTracks const& tracks)
  {
    float a = 0.f;
    for (const auto& id : jet.tracksIds()) {
      const auto trk = tracks.iteratorAt(id);
      a += std::pow(trk.pt(), kappa) * std::pow(jetutilities::deltaR(jet, trk) / (jet.r() / 100.f), alpha);
    }
    return a / std::pow(jet.pt(), kappa);
  }

  template <typename JetTableMCPConstituent>
  float jetCalculateAngularityMCP(JetTableMCPConstituent const& jet, aod::JetParticles const& particles)
  {
    float a = 0.f;
    for (const auto& id : jet.tracksIds()) {
      const auto p = particles.iteratorAt(id);
      a += std::pow(p.pt(), kappa) * std::pow(jetutilities::deltaR(jet, p) / (jet.r() / 100.f), alpha);
    }
    return a / std::pow(jet.pt(), kappa);
  }

  template <typename JetChargedTable,
            typename CandidatesTable>
  void analyseDataChargedSubstructure(aod::JetCollision const& collision,
                                      JetChargedTable const& jets,
                                      CandidatesTable const& /*candidates*/,
                                      aod::JetTracks const& tracks)
  {
    // apply event selection and fill histograms for sanity check
    registry.fill(HIST(histnames::ex_col), getValFromBin(BIN_EX_COLCNTR::AllCollisions));
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits) || !(std::abs(collision.posZ()) < vertexZCut)) {
      return;
    }
    registry.fill(HIST(histnames::ex_col), getValFromBin(BIN_EX_COLCNTR::Sel8ZCut));

    // Loop over jets containing D0 candidates
    for (const auto& jet : jets) {
      // number of charged jets with D0
      registry.fill(HIST(histnames::ex_jet), getValFromBin(BIN_EX_JETCNTR::ChargedJets));
      // obtaining jet 3-vector
      TVector3 jetVector(jet.px(), jet.py(), jet.pz());

      // Loop over D0 candidates associated to the jet
      for (const auto& d0Candidate : jet.template candidates_as<CandidatesTable>()) {
        // obtaining jet 3-vector
        TVector3 d0Vector(d0Candidate.px(), d0Candidate.py(), d0Candidate.pz());

        // calculating fraction of the jet momentum carried by the D0 along the
        // direction of the jet axis
        double zParallel = (jetVector * d0Vector) / (jetVector * jetVector);

        // calculating angular distance in eta-phi plane
        double axisDistance = jetutilities::deltaR(jet, d0Candidate);

        float angularity = jetCalculateAngularityEXP(jet, tracks);

        // filling histograms
        registry.fill(HIST(histnames::ex_jet_proj), zParallel);
        registry.fill(HIST(histnames::ex_jet_dist_proj), axisDistance, zParallel);
        registry.fill(HIST(histnames::ex_jet_dist), axisDistance);
        registry.fill(HIST(histnames::ex_jet_pt), jet.pt());
        registry.fill(HIST(histnames::ex_jet_eta), jet.eta());
        registry.fill(HIST(histnames::ex_jet_phi), jet.phi());
        registry.fill(HIST(histnames::ex_jet_ang), angularity);
        registry.fill(HIST(histnames::ex_hfl_pt), d0Candidate.pt());
        registry.fill(HIST(histnames::ex_hfl_mass), d0Candidate.m());
        registry.fill(HIST(histnames::ex_hfl_eta), d0Candidate.eta());
        registry.fill(HIST(histnames::ex_hfl_phi), d0Candidate.phi()); // add more axis

        // filling table
        objJetTable(axisDistance,
                    jet.pt(),
                    jet.eta(),
                    jet.phi(),
                    jet.template tracks_as<aod::JetTracks>().size(),
                    angularity,
                    zParallel,
                    d0Candidate.pt(),
                    d0Candidate.eta(),
                    d0Candidate.phi(),
                    d0Candidate.m(),
                    d0Candidate.y(),
                    d0Candidate.mlScores()[0],
                    d0Candidate.mlScores()[1],
                    d0Candidate.mlScores()[2]);

        break; // get out of candidates' loop after first HF particle is found
               // in jet
      } // end of D0 candidates loop

    } // end of jets loop

  } // end of process function

  template <typename MCDJetsPerMCCollissionPreslice,
            typename MCPJetsPerMCCollissionPreslice,
            typename JetTableMCD,
            typename JetTableMCP,
            typename CandidatesMCD,
            typename CandidatesMCP>
  void analyseMonteCarloEfficiency(MCDJetsPerMCCollissionPreslice const& jetmcdpreslice,
                                   MCPJetsPerMCCollissionPreslice const& jetmcppreslice,
                                   aod::JetMcCollisions const& mccollisions,
                                   aod::JetCollisionsMCD const& collisions,
                                   JetTableMCD const& mcdjets,
                                   JetTableMCP const& mcpjets,
                                   CandidatesMCD const& /*mcdCandidates*/,
                                   CandidatesMCP const& /*mcpCandidates*/,
                                   aod::JetTracks const& tracks,
                                   aod::JetParticles const& particles)
  {
    for (const auto& mccollision : mccollisions) {

      registry.fill(HIST(histnames::mc_eff_col), getValFromBin(BIN_MC_COLCNTR::All));
      // skip collisions outside of |z| < vertexZCut
      if (std::abs(mccollision.posZ()) > vertexZCut) {
        continue;
      }
      registry.fill(HIST(histnames::mc_eff_col), getValFromBin(BIN_MC_COLCNTR::ZCut));

      // reconstructed collisions associated to same mccollision
      const auto collisionsPerMCCollision = collisions.sliceBy(collisionsPerMCCollisionPreslice, mccollision.globalIndex());
      for (const auto& collision : collisionsPerMCCollision) {

        registry.fill(HIST(histnames::mc_eff_col), getValFromBin(BIN_MC_COLCNTR::Matched));
        if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits) ||
            !(std::abs(collision.posZ()) < vertexZCut)) {
          continue;
        }
        registry.fill(HIST(histnames::mc_eff_col), getValFromBin(BIN_MC_COLCNTR::MatchedSel8ZCut));

        // d0 detector level jets associated to the current same collision
        const auto d0mcdJetsPerCollision = mcdjets.sliceBy(jetmcdpreslice, collision.globalIndex());
        for (const auto& mcdjet : d0mcdJetsPerCollision) {

          registry.fill(
            HIST(histnames::mc_eff_jet), getValFromBin(BIN_MC_JETCNTR::DetectorLevelJetInMCCollision));

          // obtain leading HF candidate in jet
          auto mcdd0cand = mcdjet.template candidates_first_as<CandidatesMCD>();

          if (mcdjet.has_matchedJetCand()) {
            registry.fill(
              HIST(histnames::mc_eff_jet), getValFromBin(BIN_MC_JETCNTR::DetectorLevelJetWithMatchedCandidate));
          }

          // reflection information for storage: D0 = +1, D0bar = -1, neither =
          // 0
          int matchedFrom = 0;
          int decayChannel = o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK;
          int selectedAs = 0;

          if (mcdd0cand.flagMcMatchRec() == decayChannel) { // matched to D0 on truth level
            matchedFrom = 1;
          } else if (mcdd0cand.flagMcMatchRec() == -decayChannel) { // matched to D0bar on truth level
            matchedFrom = -1;
          }
          // bitwise AND operation: Checks whether BIT(i) is set, regardless of
          // other bits
          if (mcdd0cand.candidateSelFlag() & BIT(0)) { // CandidateSelFlag == BIT(0) -> selected as D0
            selectedAs = 1;
          } else if (mcdd0cand.candidateSelFlag() & BIT(1)) { // CandidateSelFlag == BIT(1) -> selected as D0bar
            selectedAs = -1;
          }

          TVector3 mcdjetvector(mcdjet.px(), mcdjet.py(), mcdjet.pz());
          TVector3 mcdcandvector(mcdd0cand.px(), mcdd0cand.py(), mcdd0cand.pz());

          float mcdzparallel = (mcdjetvector * mcdcandvector) / (mcdjetvector * mcdjetvector);

          float angularity = jetCalculateAngularityMCD(mcdjet, tracks);
          registry.fill(HIST(histnames::mc_eff_det_jet_pt), mcdjet.pt());
          registry.fill(HIST(histnames::mc_eff_det_jet_eta), mcdjet.eta());
          registry.fill(HIST(histnames::mc_eff_det_jet_phi), mcdjet.phi());
          registry.fill(HIST(histnames::mc_eff_det_jet_ang), angularity);
          // Particle Histgrams
          registry.fill(HIST(histnames::mc_eff_det_hfl_pt), mcdd0cand.pt());
          registry.fill(HIST(histnames::mc_eff_det_hfl_mass), mcdd0cand.m());
          registry.fill(HIST(histnames::mc_eff_det_hfl_eta), mcdd0cand.eta());
          registry.fill(HIST(histnames::mc_eff_det_hfl_phi), mcdd0cand.phi());

          mcdJetTable(
            jetutilities::deltaR(mcdjet, mcdd0cand),
            mcdjet.pt(),
            mcdjet.eta(),
            mcdjet.phi(),
            mcdjet.template tracks_as<aod::JetTracks>().size(), // detector level jet
            angularity,
            mcdzparallel,
            mcdd0cand.pt(),
            mcdd0cand.eta(), mcdd0cand.phi(),
            mcdd0cand.m(), mcdd0cand.y(),
            (mcdd0cand.originMcRec() == RecoDecay::OriginType::Prompt), // detector level D0 candidate
            mcdjet.has_matchedJetCand(),
            mcdd0cand.mlScores()[0],
            mcdd0cand.mlScores()[1],
            mcdd0cand.mlScores()[2], // // Machine Learning PID scores: background, prompt, non-prompt

            matchedFrom,
            selectedAs); // D0 = +1, D0bar = -1, neither = 0
        }
      } // end of reconstructed collisions loop (detector level mc collisions
        // matching with real collisions)

      // d0 particle level jets associated to same mccollision
      const auto d0mcpJetsPerMCCollision = mcpjets.sliceBy(jetmcppreslice, mccollision.globalIndex());
      for (const auto& mcpjet : d0mcpJetsPerMCCollision) {

        registry.fill(HIST(histnames::mc_eff_jet), getValFromBin(BIN_MC_JETCNTR::ParticleLevelJetInMCCollision));

        // obtain leading HF particle in jet
        auto mcpd0cand = mcpjet.template candidates_first_as<CandidatesMCP>();

        if (mcpjet.has_matchedJetCand()) {
          registry.fill(HIST(histnames::mc_eff_jet), getValFromBin(BIN_MC_JETCNTR::ParticleLevelJetWithMatchedCandidate));
        }
        TVector3 mcpjetvector(mcpjet.px(), mcpjet.py(), mcpjet.pz());
        TVector3 mcpcandvector(mcpd0cand.px(), mcpd0cand.py(), mcpd0cand.pz());

        float mcpzparallel = (mcpjetvector * mcpcandvector) / (mcpjetvector * mcpjetvector);

        float angularity = jetCalculateAngularityMCP(mcpjet, particles);

        registry.fill(HIST(histnames::mc_eff_par_jet_pt), mcpjet.pt());
        registry.fill(HIST(histnames::mc_eff_par_jet_eta), mcpjet.eta());
        registry.fill(HIST(histnames::mc_eff_par_jet_phi), mcpjet.phi());
        registry.fill(HIST(histnames::mc_eff_par_jet_ang), angularity);
        // Particle Histgrams
        registry.fill(HIST(histnames::mc_eff_par_hfl_pt), mcpd0cand.pt());

        registry.fill(HIST(histnames::mc_eff_par_hfl_eta), mcpd0cand.eta());
        registry.fill(HIST(histnames::mc_eff_par_hfl_phi), mcpd0cand.phi());
        // store data in MC detector level table (calculate angular distance in
        // eta-phi plane on the fly)
        mcpJetTable(jetutilities::deltaR(mcpjet, mcpd0cand),
                    mcpjet.pt(),
                    mcpjet.eta(),
                    mcpjet.phi(),
                    mcpjet.template tracks_as<aod::JetParticles>().size(), // particle level jet
                    angularity,
                    mcpzparallel,
                    mcpd0cand.pt(),
                    mcpd0cand.eta(),
                    mcpd0cand.phi(),
                    // mcpd0cand.m(),
                    mcpd0cand.y(),
                    (mcpd0cand.originMcGen() == RecoDecay::OriginType::Prompt), // particle level D0
                    mcpjet.has_matchedJetCand());
      } // End of particle level jets loop, related to detector level collisions
        // and jets.
    }
  }

  template <typename MCPJetsPerMCCollissionPreslice,
            typename JetTableMCD,
            typename JetTableMCP,
            typename CandidatesMCD,
            typename CandidatesMCP>
  void analyseMonteCarlo(MCPJetsPerMCCollissionPreslice jetmcpreslice,
                         aod::JetMcCollisions const& mccollisions,
                         aod::JetCollisionsMCD const& collisions,
                         JetTableMCD const& /*mcdjets*/,
                         JetTableMCP const& mcpjets,
                         CandidatesMCD const& /*mcdCandidates*/,
                         CandidatesMCP const& /*mcpCandidates*/,
                         aod::JetTracks const& jettracks,
                         aod::JetParticles const& jetparticles)
  {
    for (const auto& mccollision : mccollisions) {
      registry.fill(HIST(histnames::mc_col), getValFromBin(BIN_MC_COLCNTR::All));
      // skip collisions outside of |z| < vertexZCut
      if (std::abs(mccollision.posZ()) > vertexZCut) {
        continue;
      }
      registry.fill(HIST(histnames::mc_col), getValFromBin(BIN_MC_COLCNTR::ZCut));

      // hf particle level jets associated to same mccollision
      const JetD0MCPTable mcpJetsPerMCCollision = mcpjets.sliceBy(jetmcpreslice, mccollision.globalIndex());
      for (const auto& mcpjet : mcpJetsPerMCCollision) {

        registry.fill(HIST(histnames::mc_jet), getValFromBin(BIN_MC_JETCNTR::ParticleLevelJetInMCCollision));

        // obtain leading HF particle in jet
        auto mcpcand = mcpjet.template candidates_first_as<CandidatesMCP>();

        TVector3 mcpjetvector(mcpjet.px(), mcpjet.py(), mcpjet.pz());
        TVector3 mcpcandvector(mcpcand.px(), mcpcand.py(), mcpcand.pz());
        float mcpzparallel = (mcpjetvector * mcpcandvector) / (mcpjetvector * mcpjetvector);

        if (mcpjet.has_matchedJetCand()) {
          registry.fill(HIST(histnames::mc_jet), getValFromBin(BIN_MC_JETCNTR::ParticleLevelJetWithMatchedCandidate));

          // loop over detector level matched to current particle level
          for (const auto& mcdjet : mcpjet.template matchedJetCand_as<JetTableMCD>()) {
            registry.fill(HIST(histnames::mc_jet), getValFromBin(BIN_MC_JETCNTR::DetectorLevelJetWithMatchedCandidate));

            // apply collision sel8 selection on detector level jet's collision
            const auto& collision = collisions.iteratorAt(mcdjet.collisionId());
            registry.fill(HIST(histnames::mc_col), getValFromBin(BIN_MC_COLCNTR::Matched));
            if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits) || !(std::abs(collision.posZ()) < vertexZCut)) {
              continue;
            }
            registry.fill(HIST(histnames::mc_col), getValFromBin(BIN_MC_COLCNTR::MatchedSel8ZCut));

            // obtain leading HF candidate in jet
            auto mcdcand = mcdjet.template candidates_first_as<CandidatesMCD>();

            // reflection information for storage: HF = +1, HFbar = -1, neither
            // = 0
            int matchedFrom = 0;
            int decayChannel = 0;
            if (jethfutilities::isD0Table<CandidatesMCD>()) {
              decayChannel = o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK;
            } else if (jethfutilities::isLcTable<CandidatesMCD>()) {
              decayChannel = o2::hf_decay::hf_cand_3prong::DecayChannelMain::LcToPKPi;
            }
            int selectedAs = 0;

            if (mcdcand.flagMcMatchRec() == decayChannel) { // matched to HF on truth level
              matchedFrom = 1;
            } else if (mcdcand.flagMcMatchRec() == -decayChannel) { // matched to HFbar on truth level
              matchedFrom = -1;
            }
            // bitwise AND operation: Checks whether BIT(i) is set, regardless
            // of other bits
            if (mcdcand.candidateSelFlag() & BIT(0)) { // CandidateSelFlag == BIT(0) -> selected as HF
              selectedAs = 1;
            } else if (mcdcand.candidateSelFlag() & BIT(1)) { // CandidateSelFlag == BIT(1) -> selected as
                                                              // HFbar
              selectedAs = -1;
            }

            TVector3 mcdjetvector(mcdjet.px(), mcdjet.py(), mcdjet.pz());
            TVector3 mcdcandvector(mcdcand.px(), mcdcand.py(), mcdcand.pz());

            float mcdzparallel = (mcdjetvector * mcdcandvector) / (mcdjetvector * mcdjetvector);

            float mcpAngularity = jetCalculateAngularityMCP(mcpjet, jetparticles);
            float mcdAngularity = jetCalculateAngularityMCD(mcdjet, jettracks);

            // filling histograms
            // registry.fill(HIST(histnames::h_exp_d0_jet_projection),
            // zParallel);
            // registry.fill(HIST(histnames::h_exp_d0_jet_distance_vs_projection),
            // axisDistance, zParallel);
            // registry.fill(HIST(histnames::h_exp_d0_jet_distance),
            // axisDistance); registry.fill(HIST(histnames::h_exp_d0_jet_pt),
            // jet.pt()); registry.fill(HIST(histnames::h_exp_d0_jet_eta),
            // jet.eta()); registry.fill(HIST(histnames::h_exp_d0_jet_phi),
            // jet.phi()); registry.fill(HIST(histnames::h_exp_d0_jet_ang),
            // angularity); Jet Histograms

            registry.fill(HIST(histnames::mc_det_jet_pt), mcdjet.pt());
            registry.fill(HIST(histnames::mc_det_jet_eta), mcdjet.eta());
            registry.fill(HIST(histnames::mc_det_jet_phi), mcdjet.phi());
            registry.fill(HIST(histnames::mc_det_jet_ang), mcdAngularity);
            // Particle Histgrams
            registry.fill(HIST(histnames::mc_det_hfl_pt), mcdcand.pt());
            registry.fill(HIST(histnames::mc_det_hfl_mass), mcdcand.m());
            registry.fill(HIST(histnames::mc_det_hfl_eta), mcdcand.eta());
            registry.fill(HIST(histnames::mc_det_hfl_phi), mcdcand.phi()); // add more axis

            // store matched particle and detector level data in one single
            // table (calculate angular distance in eta-phi plane on the fly)
            matchJetTable(jetutilities::deltaR(mcpjet, mcpcand),
                          mcpjet.pt(),
                          mcpjet.eta(),
                          mcpjet.phi(),
                          mcpjet.template tracks_as<aod::JetParticles>().size(), // particle level jet
                          mcpAngularity,
                          mcpzparallel,
                          mcpcand.pt(),
                          mcpcand.eta(),
                          mcpcand.phi(),
                          mcpcand.y(),
                          (mcpcand.originMcGen() == RecoDecay::OriginType::Prompt), // particle level HF
                          jetutilities::deltaR(mcdjet, mcdcand),
                          mcdjet.pt(),
                          mcdjet.eta(),
                          mcdjet.phi(),
                          mcdjet.template tracks_as<aod::JetTracks>().size(), // detector level jet
                          mcdAngularity,
                          mcdzparallel,
                          mcdcand.pt(),
                          mcdcand.eta(),
                          mcdcand.phi(),
                          mcdcand.m(),
                          mcdcand.y(),
                          (mcdcand.originMcRec() == RecoDecay::OriginType::Prompt), // detector level HF
                          mcdcand.mlScores()[0],
                          mcdcand.mlScores()[1],
                          mcdcand.mlScores()[2], // Machine Learning PID scores: background, prompt, non-prompt
                          matchedFrom,
                          selectedAs); // HF = +1, HFbar = -1, neither = 0
          }
        } else {
          // store matched particle and detector level data in one single table
          // (calculate angular distance in eta-phi plane on the fly)
          float mcpAngularity = jetCalculateAngularityMCP(mcpjet, jetparticles);
          // float mcpAngularity = 0.;
          matchJetTable(jetutilities::deltaR(mcpjet, mcpcand),
                        mcpjet.pt(),
                        mcpjet.eta(),
                        mcpjet.phi(),
                        mcpjet.template tracks_as<aod::JetParticles>().size(), // particle level jet
                        mcpAngularity,
                        mcpzparallel,
                        mcpcand.pt(),
                        mcpcand.eta(),
                        mcpcand.phi(),
                        mcpcand.y(),
                        (mcpcand.originMcGen() == RecoDecay::OriginType::Prompt), // particle level HF
                        -2,
                        -2,
                        -2,
                        -2,
                        -2,
                        -2, // detector level jet
                        -2,
                        -2,
                        -2,
                        -2,
                        -2,
                        -2,
                        -2, // detector level HF
                        -2,
                        -2,
                        -2, // Machine Learning PID scores: background, prompt, non-prompt
                        -2,
                        -2); // HF = +1, HFbar = -1, neither = 0
        }
      } // end of mcpjets loop
    } // end of mccollisions loop
  };

  void processataChargedSubstructureD0(aod::JetCollision const& collision,
                                       JetChargedTableD0 const& jets,
                                       aod::CandidatesD0Data const& candidates,
                                       aod::JetTracks const& tracks)
  {
    analyseDataChargedSubstructure<JetChargedTableD0, aod::CandidatesD0Data>(collision, jets, candidates, tracks);
  }

  PROCESS_SWITCH(JetD0AngSubstructure, processataChargedSubstructureD0, "charged HF jet substructure", false);

  void processMonteCarloEfficiencyD0(aod::JetMcCollisions const& mccollisions,
                                     aod::JetCollisionsMCD const& collisions,
                                     JetD0MCDTable const& mcdjets,
                                     JetD0MCPTable const& mcpjets,
                                     aod::CandidatesD0MCD const& mcdCandidates,
                                     aod::CandidatesD0MCP const& mcpCandidates,
                                     aod::JetTracks const& jettracks,
                                     aod::JetParticles const& jetparticles)
  {
    analyseMonteCarloEfficiency<Preslice<JetD0MCDTable>,
                                Preslice<JetD0MCPTable>,
                                JetD0MCDTable,
                                JetD0MCPTable,
                                aod::CandidatesD0MCD,
                                aod::CandidatesD0MCP>(d0MCDJetsPerEXPCollisionPreslice,
                                                      d0MCPJetsPerMCCollisionPreslice,
                                                      mccollisions,
                                                      collisions,
                                                      mcdjets,
                                                      mcpjets,
                                                      mcdCandidates,
                                                      mcpCandidates,
                                                      jettracks,
                                                      jetparticles);
  }

  PROCESS_SWITCH(JetD0AngSubstructure, processMonteCarloEfficiencyD0, "non-matched and matched MC D0 and jets", false);

  void processMonteCarloD0(aod::JetMcCollisions const& mccollisions,
                           aod::JetCollisionsMCD const& collisions,
                           JetD0MCDTable const& mcdjets,
                           JetD0MCPTable const& mcpjets,
                           aod::CandidatesD0MCD const& mcdCandidates,
                           aod::CandidatesD0MCP const& mcpCandidates,
                           aod::JetTracks const& jettracks,
                           aod::JetParticles const& jetparticles)
  {
    analyseMonteCarlo<Preslice<JetD0MCPTable>,
                      JetD0MCDTable,
                      JetD0MCPTable,
                      aod::CandidatesD0MCD,
                      aod::CandidatesD0MCP>(d0MCPJetsPerMCCollisionPreslice,
                                            mccollisions,
                                            collisions,
                                            mcdjets,
                                            mcpjets,
                                            mcdCandidates,
                                            mcpCandidates,
                                            jettracks,
                                            jetparticles);
  }

  PROCESS_SWITCH(JetD0AngSubstructure, processMonteCarloD0, "Store all simulated D0 jets information with matched candidate (if any found)", false);
};
// Workflow definition
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<JetD0AngSubstructure>(cfgc)};
}
