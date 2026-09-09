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
// Measure heavy flavour angularity(D0 & Lambda_c)

/// \author Rajdeep Nandi
/// \author Preeti Dhankher

#include "PWGJE/Core/JetCandidateUtilities.h"
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

#include <array>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace
{
constexpr int kOriginMcPrompt = 1;    // MC origin flag: prompt
constexpr int kOriginMcNonPrompt = 2; // MC origin flag: non-prompt
constexpr int8_t kCandSelFlagD0 = 0;
} // namespace

namespace o2::aod
{

// Jet-side columns (shared by all three tables)
DECLARE_SOA_COLUMN(JetHfDist, jetHfDist, float);
DECLARE_SOA_COLUMN(JetPt, jetPt, float);
DECLARE_SOA_COLUMN(JetEta, jetEta, float);
DECLARE_SOA_COLUMN(JetPhi, jetPhi, float);
DECLARE_SOA_COLUMN(JetNConst, jetNConst, int);
DECLARE_SOA_COLUMN(JetAng, jetAng, float);
DECLARE_SOA_COLUMN(JetGirth, jetGirth, float);
DECLARE_SOA_COLUMN(JetMass, jetMass, float);

//  D0-candidate columns (shared by all three tables)
DECLARE_SOA_COLUMN(HfZParallel, hfZParallel, float);
DECLARE_SOA_COLUMN(HfPt, hfPt, float);
DECLARE_SOA_COLUMN(HfEta, hfEta, float);
DECLARE_SOA_COLUMN(HfPhi, hfPhi, float);
DECLARE_SOA_COLUMN(HfMass, hfMass, float);
DECLARE_SOA_COLUMN(HfY, hfY, float);

//  ML scores (one entry per BDT class) — DATA & MCD only
DECLARE_SOA_COLUMN(HfMlScore0, hfMlScore0, float);
DECLARE_SOA_COLUMN(HfMlScore1, hfMlScore1, float);
DECLARE_SOA_COLUMN(HfMlScore2, hfMlScore2, float);

//  MC-truth columns — MCD & MCP only
DECLARE_SOA_COLUMN(HfFlagMcMatch, hfFlagMcMatch, int8_t);
DECLARE_SOA_COLUMN(HfOriginMc, hfOriginMc, int8_t);

DECLARE_SOA_COLUMN(JetMatchGeoStatus, jetMatchGeoStatus, int8_t);
DECLARE_SOA_COLUMN(JetMatchCandStatus, jetMatchCandStatus, int8_t);
DECLARE_SOA_COLUMN(JetMatchClean, jetMatchClean, int8_t);
DECLARE_SOA_COLUMN(JetMatchedPt, jetMatchedPt, float);
DECLARE_SOA_COLUMN(JetMatchedEta, jetMatchedEta, float);
DECLARE_SOA_COLUMN(JetMatchedPhi, jetMatchedPhi, float);
DECLARE_SOA_COLUMN(JetMatchedDR, jetMatchedDR, float);
DECLARE_SOA_COLUMN(JetMatchedD0Pt, jetMatchedD0Pt, float);

DECLARE_SOA_TABLE(D0JetObjTable, "AOD", "D0JETOBJTABLE",
                  JetHfDist,
                  JetPt,
                  JetEta,
                  JetPhi,
                  JetNConst,
                  JetAng,
                  JetGirth,
                  JetMass,
                  HfZParallel,
                  HfPt,
                  HfEta,
                  HfPhi,
                  HfMass,
                  HfY,
                  HfMlScore0,
                  HfMlScore1,
                  HfMlScore2);

DECLARE_SOA_TABLE(D0JetMCDObjTable, "AOD", "D0JETMCDOBJ",
                  JetHfDist,
                  JetPt,
                  JetEta,
                  JetPhi,
                  JetNConst,
                  JetAng,
                  JetGirth,
                  JetMass,
                  HfZParallel,
                  HfPt,
                  HfEta,
                  HfPhi,
                  HfMass,
                  HfY,
                  HfMlScore0,
                  HfMlScore1,
                  HfMlScore2,
                  HfFlagMcMatch,
                  HfOriginMc);

DECLARE_SOA_TABLE(D0JetMCPObjTable, "AOD", "D0JETMCPOBJ",
                  JetHfDist,
                  JetPt,
                  JetEta,
                  JetPhi,
                  JetNConst,
                  JetAng,
                  JetGirth,
                  JetMass,
                  HfZParallel,
                  HfPt,
                  HfEta,
                  HfPhi,
                  HfY,
                  HfFlagMcMatch,
                  HfOriginMc);

DECLARE_SOA_TABLE(D0JetMCDObjMatchedTable, "AOD", "D0JMCDOBJMATCH",
                  JetHfDist,
                  JetPt,
                  JetEta,
                  JetPhi,
                  JetNConst,
                  JetAng,
                  JetGirth,
                  JetMass,
                  HfZParallel,
                  HfPt,
                  HfEta,
                  HfPhi,
                  HfMass,
                  HfY,
                  HfMlScore0,
                  HfMlScore1,
                  HfMlScore2,
                  HfFlagMcMatch,
                  HfOriginMc,
                  JetMatchGeoStatus,
                  JetMatchCandStatus,
                  JetMatchClean,
                  JetMatchedPt,
                  JetMatchedEta,
                  JetMatchedPhi,
                  JetMatchedDR,
                  JetMatchedD0Pt);

} // namespace o2::aod

consteval float getValFromBin(int bin)
{
  return static_cast<float>(bin) - 0.5f;
}

/// Collision-counter bins
enum BinCollCntr {
  AllCollisions = 1,
  Sel8ZCut = 2
};

/// Jet-counter bins
enum BinJetCntr {
  ChargedJets = 1,
  D0TaggedJets = 2
};

/// MC-efficiency collision-counter bins (jet-matching based efficiency)
enum BinMCColCntr {
  AllMCCollisions = 1,
  SelectedMCCollisions = 2,
  AssociatedRecoCollisions = 3,
  SelectedAssociatedRecoCollisions = 4
};

/// MC-efficiency jet-counter bins (jet-matching based efficiency)
enum BinMCJetCntr {
  DetectorLevelJetInMCCollision = 1,
  ParticleLevelJetInMCCollision = 2,
  DetectorLevelJetWithMatchedCandidate = 3,
  ParticleLevelJetWithMatchedCandidate = 4
};

struct JetHFAngularityTask {

  // DATA

  using D0CandidatesData = aod::CandidatesD0Data;

  using D0DataJets = soa::Join<aod::D0ChargedJets,
                               aod::D0ChargedJetConstituents>;

  // MC-DETECTOR LEVEL(MCD)

  using D0CandidatesMCD = aod::CandidatesD0MCD;

  using D0MCDJets = soa::Join<aod::D0ChargedMCDetectorLevelJets,
                              aod::D0ChargedMCDetectorLevelJetConstituents>;

  using D0MCDJetsMatched = soa::Join<aod::D0ChargedMCDetectorLevelJets,
                                     aod::D0ChargedMCDetectorLevelJetConstituents,
                                     aod::D0ChargedMCDetectorLevelJetsMatchedToD0ChargedMCParticleLevelJets>;

  // MC PARTICLE LEVEL (MCP)

  using D0CandidatesMCP = aod::CandidatesD0MCP;

  using D0MCPJets = soa::Join<aod::D0ChargedMCParticleLevelJets,
                              aod::D0ChargedMCParticleLevelJetConstituents>;

  using D0MCPJetsMatched = soa::Join<aod::D0ChargedMCParticleLevelJets,
                                     aod::D0ChargedMCParticleLevelJetConstituents,
                                     aod::D0ChargedMCParticleLevelJetsMatchedToD0ChargedMCDetectorLevelJets>;

  // Configurables

  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range (cm)"};
  Configurable<float> jetPtMin{"jetPtMin", 5.0f, "Minimum jet pT (GeV/c)"};
  Configurable<float> jetR{"jetR", 0.4f, "Jet resolution parameter R"};

  Configurable<std::string> eventSelections{"eventSelections", "sel8",
                                            "Event selection string"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks",
                                            "Track selection string"};

  Produces<aod::D0JetObjTable> objJetTable;
  Produces<aod::D0JetMCDObjTable> objJetMCDTable;
  Produces<aod::D0JetMCPObjTable> objJetMCPTable;

  Produces<aod::D0JetMCDObjMatchedTable> objJetMCDMatchedTable;

  std::vector<int> eventSelectionBits;
  int trackSelection = -1;

  float massD0MCP = -1.f;

  // Slices for the jet-matching based MC efficiency calculation

  // Reconstructed collisions associated to a given JetMcCollision
  PresliceUnsorted<aod::JetCollisionsMCD> collisionsPerMCCollisionPreslice = aod::jmccollisionlb::mcCollisionId;

  // D0-tagged (matched) MCD jets per reconstructed collision
  Preslice<D0MCDJetsMatched> d0MCDJetsPerCollisionPreslice = aod::jet::collisionId;

  // D0-tagged (matched) MCP jets per MC collision
  Preslice<D0MCPJetsMatched> d0MCPJetsPerMCCollisionPreslice = aod::jet::mcCollisionId;

  // Filters

  Filter jetCutsPt = aod::jet::pt > jetPtMin;
  Filter jetCutsR = aod::jet::r == nround(jetR.node() * 100.0f);
  Filter collisionFilter = nabs(aod::jcollision::posZ) < vertexZCut;

  // Histograms

  HistogramRegistry registry{"registry", {

                                           //                 DATA

                                           // ----- Collision QA -----
                                           {"h_collisions", "Event status;status;entries", {HistType::kTH1F, {{4, 0., 4.}}}},

                                           // ----- Track QA -----
                                           {"h_track_pt", ";#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                                           {"h_track_eta", ";#eta_{track};entries", {HistType::kTH1F, {{100, -1., 1.}}}},
                                           {"h_track_phi", ";#varphi_{track};entries", {HistType::kTH1F, {{80, -1., 7.}}}},

                                           // ----- Jet QA -----
                                           {"h_jet_counter", ";type;counts", {HistType::kTH1F, {{4, 0., 4.}}}},

                                           {"h_jet_pt", ";#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                                           {"h_jet_eta", ";#eta_{jet};entries", {HistType::kTH1F, {{100, -1., 1.}}}},
                                           {"h_jet_phi", ";#phi_{jet};entries", {HistType::kTH1F, {{80, -1., 7.}}}},

                                           // ----- D0 candidate properties -----
                                           {"h_d0_mass", ";m_{K#pi} (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{300, 1.7, 2.0}}}},
                                           {"h_d0_pt", ";#it{p}_{T,D^{0}} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 100.}}}},
                                           {"h_d0_eta", ";#eta_{D^{0}};entries", {HistType::kTH1F, {{100, -1., 1.}}}},
                                           {"h_d0_phi", ";#phi_{D^{0}};entries", {HistType::kTH1F, {{100, -1., 7.}}}},
                                           {"h_d0_y", ";y_{D^{0}};entries", {HistType::kTH1F, {{100, -2., 2.}}}},

                                           // ----- D0-in-jet observables -----
                                           {"h_d0_jet_projection", ";z_{||}^{D^{0},jet};entries", {HistType::kTH1F, {{200, 0., 2.}}}},

                                           {"h_d0_jet_distance", ";#DeltaR_{D^{0},jet};entries", {HistType::kTH1F, {{200, 0., 1.}}}},

                                           // ----- Jet substructure (D0-tagged jets only) -----
                                           {"h_d0_jet_lambda11", ";#lambda_{1}^{1} (angularity);entries", {HistType::kTH1F, {{200, 0., 1.}}}},

                                           {"h_d0_jet_lambda12", ";#lambda_{2}^{1} (girth);entries", {HistType::kTH1F, {{200, 0., 1.}}}},

                                           {"h_d0_jet_mass", ";m_{jet} (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{200, 0., 50.}}}},

                                           // ----- Jet constituent multiplicity & ML scores
                                           {"h_jet_nconst", ";N_{const};entries", {HistType::kTH1F, {{100, -0.5, 99.5}}}},
                                           {"h_d0_mlscore0", ";ML score 0 (bkg);entries", {HistType::kTH1F, {{100, 0., 1.}}}},
                                           {"h_d0_mlscore1", ";ML score 1 (prompt);entries", {HistType::kTH1F, {{100, 0., 1.}}}},
                                           {"h_d0_mlscore2", ";ML score 2 (non-prompt);entries", {HistType::kTH1F, {{100, 0., 1.}}}},

                                           // ----- Sparse: (m_D0, pT_D0, pT_jet, z_parallel, DeltaR) -----
                                           {"hSparse_d0", ";m_{D^{0}};#it{p}_{T,D^{0}};#it{p}_{T,jet};z_{||};#DeltaR;#lambda_{1}^{1};#lambda_{2}^{1}", {HistType::kTHnSparseF, {{300, 1.7, 2.0}, {200, 0., 100.}, {200, 0., 100.}, {200, 0., 2.}, {200, 0., 1.}, {50, 0., 1.}, {50, 0., 1.}}}},

                                           //               MC DETECTOR LEVEL (MCD)

                                           {"h_collisions_mcd", "MCD event status;status;entries", {HistType::kTH1F, {{2, 0., 2.}}}},

                                           {"h_jet_counter_mcd", ";type;counts", {HistType::kTH1F, {{2, 0., 2.}}}},
                                           {"h_jet_pt_mcd", ";#it{p}_{T,jet}^{det} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                                           {"h_jet_eta_mcd", ";#eta_{jet}^{det};entries", {HistType::kTH1F, {{100, -1., 1.}}}},
                                           {"h_jet_phi_mcd", ";#phi_{jet}^{det};entries", {HistType::kTH1F, {{80, -1., 7.}}}},

                                           {"h_d0_mass_mcd", ";m_{K#pi}^{det} (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{300, 1.7, 2.0}}}},
                                           {"h_d0_pt_mcd", ";#it{p}_{T,D^{0}}^{det} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 100.}}}},
                                           {"h_d0_eta_mcd", ";#eta_{D^{0}}^{det};entries", {HistType::kTH1F, {{100, -1., 1.}}}},
                                           {"h_d0_phi_mcd", ";#phi_{D^{0}}^{det};entries", {HistType::kTH1F, {{100, -1., 7.}}}},
                                           {"h_d0_y_mcd", ";y_{D^{0}}^{det};entries", {HistType::kTH1F, {{100, -2., 2.}}}},

                                           {"h_d0_jet_projection_mcd", ";z_{||}^{D^{0},jet,det};entries", {HistType::kTH1F, {{200, 0., 2.}}}},
                                           {"h_d0_jet_distance_mcd", ";#DeltaR_{D^{0},jet}^{det};entries", {HistType::kTH1F, {{200, 0., 1.}}}},

                                           {"h_d0_jet_lambda11_mcd", ";#lambda_{1}^{1,det};entries", {HistType::kTH1F, {{200, 0., 1.}}}},
                                           {"h_d0_jet_lambda12_mcd", ";#lambda_{2}^{1,det};entries", {HistType::kTH1F, {{200, 0., 1.}}}},
                                           {"h_d0_jet_mass_mcd", ";m_{jet}^{det} (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{200, 0., 50.}}}},

                                           {"h_d0_origin_mcd", ";origin (0=none,1=prompt,2=non-prompt);entries", {HistType::kTH1F, {{3, 0., 3.}}}},

                                           {"hSparse_d0_mcd", ";m_{D^{0}}^{det};#it{p}_{T,D^{0}}^{det};#it{p}_{T,jet}^{det};z_{||}^{det};#DeltaR^{det};#lambda_{1}^{1,det};#lambda_{2}^{1,det}", {HistType::kTHnSparseF, {{300, 1.7, 2.0}, {200, 0., 100.}, {200, 0., 100.}, {200, 0., 2.}, {200, 0., 1.}, {50, 0., 1.}, {50, 0., 1.}}}},

                                           {"hSparse_d0_mcd_calib", ";m_{D^{0}}^{det};#it{p}_{T,D^{0}}^{det};#it{p}_{T,jet}^{det};z_{||}^{det};Origin(D^{0});#lambda_{1}^{1,det};#lambda_{2}^{1,det}", {HistType::kTHnSparseF, {{300, 1.7, 2.0}, {200, 0., 100.}, {200, 0., 100.}, {200, 0., 2.}, {3, -0.5, 2.5}, {50, 0., 1.}, {50, 0., 1.}}}},

                                           {"h_d0_mass_mcd_signal", ";m_{K#pi}^{det} (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{300, 1.7, 2.0}}}},
                                           {"h_d0_mass_mcd_reflection", ";m_{K#pi}^{det} (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{300, 1.7, 2.0}}}},

                                           {"hSparse_d0_mcd_reflection", ";m_{D^{0}}^{det};#it{p}_{T,D^{0}}^{det};#it{p}_{T,jet}^{det};z_{||}^{det};Origin(D^{0});#lambda_{1}^{1,det};#lambda_{2}^{1,det}", {HistType::kTHnSparseF, {{300, 1.7, 2.0}, {200, 0., 100.}, {200, 0., 100.}, {200, 0., 2.}, {3, -0.5, 2.5}, {50, 0., 1.}, {50, 0., 1.}}}},

                                           //                    MC PARTICLE LEVEL (MCP)

                                           {"h_collisions_mcp", "MCP event status;status;entries", {HistType::kTH1F, {{2, 0., 2.}}}},

                                           {"h_jet_counter_mcp", ";type;counts", {HistType::kTH1F, {{2, 0., 2.}}}},
                                           {"h_jet_pt_mcp", ";#it{p}_{T,jet}^{part} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                                           {"h_jet_eta_mcp", ";#eta_{jet}^{part};entries", {HistType::kTH1F, {{100, -1., 1.}}}},
                                           {"h_jet_phi_mcp", ";#phi_{jet}^{part};entries", {HistType::kTH1F, {{80, -1., 7.}}}},

                                           {"h_d0_pt_mcp", ";#it{p}_{T,D^{0}}^{part} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 100.}}}},
                                           {"h_d0_eta_mcp", ";#eta_{D^{0}}^{part};entries", {HistType::kTH1F, {{100, -1., 1.}}}},
                                           {"h_d0_phi_mcp", ";#phi_{D^{0}}^{part};entries", {HistType::kTH1F, {{100, -1., 7.}}}},
                                           {"h_d0_y_mcp", ";y_{D^{0}}^{part};entries", {HistType::kTH1F, {{100, -2., 2.}}}},

                                           {"h_d0_jet_projection_mcp", ";z_{||}^{D^{0},jet,part};entries", {HistType::kTH1F, {{200, 0., 2.}}}},
                                           {"h_d0_jet_distance_mcp", ";#DeltaR_{D^{0},jet}^{part};entries", {HistType::kTH1F, {{200, 0., 1.}}}},

                                           {"h_d0_jet_lambda11_mcp", ";#lambda_{1}^{1,part};entries", {HistType::kTH1F, {{200, 0., 1.}}}},
                                           {"h_d0_jet_lambda12_mcp", ";#lambda_{2}^{1,part};entries", {HistType::kTH1F, {{200, 0., 1.}}}},
                                           {"h_d0_jet_mass_mcp", ";m_{jet}^{part} (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{200, 0., 50.}}}},

                                           {"h_d0_origin_mcp", ";origin (0=none,1=prompt,2=non-prompt);entries", {HistType::kTH1F, {{3, 0., 3.}}}},

                                           //           JET-TO-JET MATCHING (MCD <-> MCP)

                                           {"h_jet_matching_geo_status_mcd", ";status (0=unmatched,1=matched);entries", {HistType::kTH1F, {{2, 0., 2.}}}},
                                           {"h_jet_matching_cand_status_mcd", ";status (0=unmatched,1=matched);entries", {HistType::kTH1F, {{2, 0., 2.}}}},
                                           {"h_jet_matching_clean_mcd", ";status (0=not clean,1=clean);entries", {HistType::kTH1F, {{2, 0., 2.}}}},

                                           {"h_jet_matching_geocand_disagree_mcd", ";geo & cand agree? (0=no,1=yes);entries", {HistType::kTH1F, {{2, 0., 2.}}}},

                                           {"h_jet_matching_dr_mcd", ";#DeltaR(jet^{det},jet^{part});entries", {HistType::kTH1F, {{100, 0., 0.5}}}},

                                           {"h_jet_matching_ngeo_mcd", ";number of geo-matched candidates;entries", {HistType::kTH1F, {{6, 0., 6.}}}},

                                           {"h_jet_matching_ncand_mcd", ";number of cand-matched candidates;entries", {HistType::kTH1F, {{6, 0., 6.}}}},

                                           {"h_jet_matching_dr_mcd_allcand", ";#DeltaR(jet^{det},jet^{part}) (all geo candidates);entries", {HistType::kTH1F, {{100, 0., 0.5}}}},

                                           {"h_jet_pt_response_matrix_allcand", ";#it{p}_{T,jet}^{det} (GeV/#it{c});#it{p}_{T,jet}^{part} (GeV/#it{c}) (all geo candidates)", {HistType::kTH2F, {{200, 0., 200.}, {200, 0., 200.}}}},

                                           // Detector-level vs matched particle-level jet pT

                                           {"h_jet_pt_response_matrix", ";#it{p}_{T,jet}^{det} (GeV/#it{c});#it{p}_{T,jet}^{part} (GeV/#it{c})", {HistType::kTH2F, {{200, 0., 200.}, {200, 0., 200.}}}},
                                           {"h_jet_pt_response_matrix_pr", ";#it{p}_{T,jet}^{det} (GeV/#it{c});#it{p}_{T,jet}^{part} (GeV/#it{c}) (prompt)", {HistType::kTH2F, {{200, 0., 200.}, {200, 0., 200.}}}},
                                           {"h_jet_pt_response_matrix_np", ";#it{p}_{T,jet}^{det} (GeV/#it{c});#it{p}_{T,jet}^{part} (GeV/#it{c}) (non-prompt)", {HistType::kTH2F, {{200, 0., 200.}, {200, 0., 200.}}}},

                                           // Axes: pT_jet^det, pT_D0^det, pT_jet^part, pT_D0^part
                                           {"h_response_d0pt_pr", ";#it{p}_{T,jet}^{det} (GeV/#it{c});#it{p}_{T,D^{0}}^{det} (GeV/#it{c});#it{p}_{T,jet}^{part} (GeV/#it{c});#it{p}_{T,D^{0}}^{part} (GeV/#it{c})", {HistType::kTHnSparseF, {{40, 0., 200.}, {50, 0., 50.}, {40, 0., 200.}, {50, 0., 50.}}}},
                                           {"h_response_d0pt_np", ";#it{p}_{T,jet}^{det} (GeV/#it{c});#it{p}_{T,D^{0}}^{det} (GeV/#it{c});#it{p}_{T,jet}^{part} (GeV/#it{c});#it{p}_{T,D^{0}}^{part} (GeV/#it{c})", {HistType::kTHnSparseF, {{40, 0., 200.}, {50, 0., 50.}, {40, 0., 200.}, {50, 0., 50.}}}},

                                           {"h_kineff_gen_num_pr", ";#it{p}_{T,jet}^{part} (GeV/#it{c});#it{p}_{T,D^{0}}^{part} (GeV/#it{c})", {HistType::kTH2F, {{40, 0., 200.}, {50, 0., 50.}}}},
                                           {"h_kineff_gen_den_pr", ";#it{p}_{T,jet}^{part} (GeV/#it{c});#it{p}_{T,D^{0}}^{part} (GeV/#it{c})", {HistType::kTH2F, {{40, 0., 200.}, {50, 0., 50.}}}},
                                           {"h_kineff_gen_num_np", ";#it{p}_{T,jet}^{part} (GeV/#it{c});#it{p}_{T,D^{0}}^{part} (GeV/#it{c})", {HistType::kTH2F, {{40, 0., 200.}, {50, 0., 50.}}}},
                                           {"h_kineff_gen_den_np", ";#it{p}_{T,jet}^{part} (GeV/#it{c});#it{p}_{T,D^{0}}^{part} (GeV/#it{c})", {HistType::kTH2F, {{40, 0., 200.}, {50, 0., 50.}}}},

                                           {"h_kineff_det_num_pr", ";#it{p}_{T,jet}^{det} (GeV/#it{c});#it{p}_{T,D^{0}}^{det} (GeV/#it{c})", {HistType::kTH2F, {{40, 0., 200.}, {50, 0., 50.}}}},
                                           {"h_kineff_det_den_pr", ";#it{p}_{T,jet}^{det} (GeV/#it{c});#it{p}_{T,D^{0}}^{det} (GeV/#it{c})", {HistType::kTH2F, {{40, 0., 200.}, {50, 0., 50.}}}},
                                           {"h_kineff_det_num_np", ";#it{p}_{T,jet}^{det} (GeV/#it{c});#it{p}_{T,D^{0}}^{det} (GeV/#it{c})", {HistType::kTH2F, {{40, 0., 200.}, {50, 0., 50.}}}},
                                           {"h_kineff_det_den_np", ";#it{p}_{T,jet}^{det} (GeV/#it{c});#it{p}_{T,D^{0}}^{det} (GeV/#it{c})", {HistType::kTH2F, {{40, 0., 200.}, {50, 0., 50.}}}},
                                           //             ANGULARITY RESPONSE MATRIX

                                           {"h_response_angularity_pr", ";#it{p}_{T,jet}^{det} (GeV/#it{c});#lambda_{1}^{1,det};#it{p}_{T,jet}^{part} (GeV/#it{c});#lambda_{1}^{1,part}", {HistType::kTHnSparseF, {{40, 0., 200.}, {50, 0., 1.}, {40, 0., 200.}, {50, 0., 1.}}}},
                                           {"h_response_angularity_np", ";#it{p}_{T,jet}^{det} (GeV/#it{c});#lambda_{1}^{1,det};#it{p}_{T,jet}^{part} (GeV/#it{c});#lambda_{1}^{1,part}", {HistType::kTHnSparseF, {{40, 0., 200.}, {50, 0., 1.}, {40, 0., 200.}, {50, 0., 1.}}}},
                                           {"h_response_girth_pr", ";#it{p}_{T,jet}^{det} (GeV/#it{c});#lambda_{2}^{1,det};#it{p}_{T,jet}^{part} (GeV/#it{c});#lambda_{2}^{1,part}", {HistType::kTHnSparseF, {{40, 0., 200.}, {50, 0., 1.}, {40, 0., 200.}, {50, 0., 1.}}}},
                                           {"h_response_girth_np", ";#it{p}_{T,jet}^{det} (GeV/#it{c});#lambda_{2}^{1,det};#it{p}_{T,jet}^{part} (GeV/#it{c});#lambda_{2}^{1,part}", {HistType::kTHnSparseF, {{40, 0., 200.}, {50, 0., 1.}, {40, 0., 200.}, {50, 0., 1.}}}},

                                           {"h_eff_run2_num_pr", ";#it{p}_{T,D^{0}}^{det} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 50.}}}},
                                           {"h_eff_run2_num_np", ";#it{p}_{T,D^{0}}^{det} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 50.}}}},
                                           {"h_eff_run2_den_pr", ";#it{p}_{T,D^{0}}^{part} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 50.}}}},
                                           {"h_eff_run2_den_np", ";#it{p}_{T,D^{0}}^{part} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 50.}}}},

                                           //        MC EFFICIENCY (jet-to-jet matching based)

                                           {"MC_CollisionCounter", "N_{collisions};", {HistType::kTH1F, {{4, 0., 4.}}}},
                                           {"MC_JetCounter", "N_{jet};", {HistType::kTH1F, {{4, 0., 4.}}}},

                                           {"hSparse_d0_mcd_eff", ";m_{D^{0}}^{det};#it{p}_{T,D^{0}}^{det};#it{p}_{T,jet}^{det};z_{||}^{det};Origin(D^{0});Matching status;#lambda_{1}^{1,det};#lambda_{2}^{1,det};IsReflection", {HistType::kTHnSparseF, {{300, 1.7, 2.0}, {200, 0., 100.}, {200, 0., 100.}, {200, 0., 2.}, {3, -0.5, 2.5}, {2, -0.5, 1.5}, {50, 0., 1.}, {50, 0., 1.}, {2, -0.5, 1.5}}}},

                                           {"hSparse_d0_mcp_eff", ";#it{p}_{T,D^{0}}^{part};#it{p}_{T,jet}^{part};z_{||}^{part};Origin(D^{0});Matching status;#lambda_{1}^{1,part};#lambda_{2}^{1,part}", {HistType::kTHnSparseF, {{200, 0., 100.}, {200, 0., 100.}, {200, 0., 2.}, {3, -0.5, 2.5}, {2, -0.5, 1.5}, {50, 0., 1.}, {50, 0., 1.}}}},

                                         }};

  // init()

  void init(InitContext const&)
  {
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(eventSelections.value);
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(trackSelections.value);

    massD0MCP = jetcandidateutilities::getTablePDGMass<D0CandidatesMCP>();

    // ---- DATA bin labels ----
    auto hColl = registry.get<TH1>(HIST("h_collisions"));
    hColl->GetXaxis()->SetBinLabel(AllCollisions, "all");
    hColl->GetXaxis()->SetBinLabel(Sel8ZCut, "sel8 + |z| cut");
    hColl->GetXaxis()->SetBinLabel(3, "has D0 jet");
    hColl->GetXaxis()->SetBinLabel(4, "rejected");

    auto hJet = registry.get<TH1>(HIST("h_jet_counter"));
    hJet->GetXaxis()->SetBinLabel(ChargedJets, "All charged jets");
    hJet->GetXaxis()->SetBinLabel(D0TaggedJets, "D0-tagged jets");

    // ---- MCD bin labels ----
    auto hCollMCD = registry.get<TH1>(HIST("h_collisions_mcd"));
    hCollMCD->GetXaxis()->SetBinLabel(1, "selected collisions");
    hCollMCD->GetXaxis()->SetBinLabel(2, "has D0 jet");

    auto hJetMCD = registry.get<TH1>(HIST("h_jet_counter_mcd"));
    hJetMCD->GetXaxis()->SetBinLabel(ChargedJets, "All charged jets");
    hJetMCD->GetXaxis()->SetBinLabel(D0TaggedJets, "D0-tagged jets");

    auto hOriginMCD = registry.get<TH1>(HIST("h_d0_origin_mcd"));
    hOriginMCD->GetXaxis()->SetBinLabel(1, "none/bkg");
    hOriginMCD->GetXaxis()->SetBinLabel(2, "prompt");
    hOriginMCD->GetXaxis()->SetBinLabel(3, "non-prompt");

    // ---- MCP bin labels ----
    auto hCollMCP = registry.get<TH1>(HIST("h_collisions_mcp"));
    hCollMCP->GetXaxis()->SetBinLabel(1, "all MC collisions");
    hCollMCP->GetXaxis()->SetBinLabel(2, "has D0 jet");

    auto hJetMCP = registry.get<TH1>(HIST("h_jet_counter_mcp"));
    hJetMCP->GetXaxis()->SetBinLabel(ChargedJets, "All charged jets");
    hJetMCP->GetXaxis()->SetBinLabel(D0TaggedJets, "D0-tagged jets");

    auto hOriginMCP = registry.get<TH1>(HIST("h_d0_origin_mcp"));
    hOriginMCP->GetXaxis()->SetBinLabel(1, "none/bkg");
    hOriginMCP->GetXaxis()->SetBinLabel(2, "prompt");
    hOriginMCP->GetXaxis()->SetBinLabel(3, "non-prompt");

    // ---- Jet-matching bin labels  ----
    auto hMatchGeoMCD = registry.get<TH1>(HIST("h_jet_matching_geo_status_mcd"));
    hMatchGeoMCD->GetXaxis()->SetBinLabel(1, "unmatched");
    hMatchGeoMCD->GetXaxis()->SetBinLabel(2, "matched");

    auto hMatchCandMCD = registry.get<TH1>(HIST("h_jet_matching_cand_status_mcd"));
    hMatchCandMCD->GetXaxis()->SetBinLabel(1, "unmatched");
    hMatchCandMCD->GetXaxis()->SetBinLabel(2, "matched");

    auto hMatchCleanMCD = registry.get<TH1>(HIST("h_jet_matching_clean_mcd"));
    hMatchCleanMCD->GetXaxis()->SetBinLabel(1, "not clean");
    hMatchCleanMCD->GetXaxis()->SetBinLabel(2, "clean");

    auto hGeoCandDisagreeMCD = registry.get<TH1>(HIST("h_jet_matching_geocand_disagree_mcd"));
    hGeoCandDisagreeMCD->GetXaxis()->SetBinLabel(1, "disagree");
    hGeoCandDisagreeMCD->GetXaxis()->SetBinLabel(2, "agree");

    // ---- MC efficiency (jet-to-jet matching based) bin labels ----
    auto mcCollisionCounter = registry.get<TH1>(HIST("MC_CollisionCounter"));
    mcCollisionCounter->GetXaxis()->SetBinLabel(BinMCColCntr::AllMCCollisions, "All MC collisions");
    mcCollisionCounter->GetXaxis()->SetBinLabel(BinMCColCntr::SelectedMCCollisions, "Selected MC collisions");
    mcCollisionCounter->GetXaxis()->SetBinLabel(BinMCColCntr::AssociatedRecoCollisions, "MC-Associated reco collisions");
    mcCollisionCounter->GetXaxis()->SetBinLabel(BinMCColCntr::SelectedAssociatedRecoCollisions, "Selected MC-Associated reco collisions");

    auto jetCounter = registry.get<TH1>(HIST("MC_JetCounter"));
    jetCounter->GetXaxis()->SetBinLabel(BinMCJetCntr::DetectorLevelJetInMCCollision, "Detector-level D0 jets");
    jetCounter->GetXaxis()->SetBinLabel(BinMCJetCntr::ParticleLevelJetInMCCollision, "Particle-level D0 jets");
    jetCounter->GetXaxis()->SetBinLabel(BinMCJetCntr::DetectorLevelJetWithMatchedCandidate, "Detector jets matched to particle");
    jetCounter->GetXaxis()->SetBinLabel(BinMCJetCntr::ParticleLevelJetWithMatchedCandidate, "Particle jets matched to detector");

    auto hSparseMCDEff = registry.get<THnSparse>(HIST("hSparse_d0_mcd_eff"));
    hSparseMCDEff->GetAxis(4)->SetBinLabel(1, "none/bkg");
    hSparseMCDEff->GetAxis(4)->SetBinLabel(2, "prompt");
    hSparseMCDEff->GetAxis(4)->SetBinLabel(3, "non-prompt");
    hSparseMCDEff->GetAxis(5)->SetBinLabel(1, "Unmatched");
    hSparseMCDEff->GetAxis(5)->SetBinLabel(2, "Matched");
    hSparseMCDEff->GetAxis(8)->SetBinLabel(1, "Not reflection");
    hSparseMCDEff->GetAxis(8)->SetBinLabel(2, "Reflection");

    auto hSparseMCPEff = registry.get<THnSparse>(HIST("hSparse_d0_mcp_eff"));
    hSparseMCPEff->GetAxis(3)->SetBinLabel(1, "none/bkg");
    hSparseMCPEff->GetAxis(3)->SetBinLabel(2, "prompt");
    hSparseMCPEff->GetAxis(3)->SetBinLabel(3, "non-prompt");
    hSparseMCPEff->GetAxis(4)->SetBinLabel(1, "Unmatched");
    hSparseMCPEff->GetAxis(4)->SetBinLabel(2, "Matched");

    // ---- Reflection sparse bin labels ----
    auto hSparseMCDRefl = registry.get<THnSparse>(HIST("hSparse_d0_mcd_reflection"));
    hSparseMCDRefl->GetAxis(4)->SetBinLabel(1, "none/bkg");
    hSparseMCDRefl->GetAxis(4)->SetBinLabel(2, "prompt");
    hSparseMCDRefl->GetAxis(4)->SetBinLabel(3, "non-prompt");
  }

  // Helper: reflection identification (MCD only)

  template <typename CAND>
  bool isReflectedCandidate(CAND const& cand)
  {
    const int8_t flagMcMatchRec = cand.flagMcMatchRec();
    if (flagMcMatchRec == 0) {
      return false; // not truth-matched -> background, not a reflection
    }
    const bool trueIsD0 = flagMcMatchRec > 0;                                       // sign encodes true species
    const bool thisRowIsD0Hypothesis = (cand.candidateSelFlag() == kCandSelFlagD0); // which hypothesis this row uses
    return thisRowIsD0Hypothesis != trueIsD0;
  }

  // Helper: angularity

  template <typename JET, typename TRACKS, typename CANDIDATES>
  float computeLambda(JET const& jet, TRACKS const& tracks, CANDIDATES const& candidates,
                      float alpha, float kappa)
  {
    if (jet.pt() <= 0.f) {
      return -1.f;
    }

    float sum = 0.f;
    for (auto const& trk : tracks) {
      const float dr = jetutilities::deltaR(jet, trk);
      sum += std::pow(trk.pt(), kappa) * std::pow(dr, alpha);
    }
    for (auto const& cand : candidates) {
      const float dr = jetutilities::deltaR(jet, cand);
      sum += std::pow(cand.pt(), kappa) * std::pow(dr, alpha);
    }

    const float jetRadius = jet.r() / 100.f; // stored as R×100 in the table
    const float denom = std::pow(jet.pt(), kappa) * std::pow(jetRadius, alpha);
    if (denom <= 0.f) {
      return -1.f;
    }
    return sum / denom;
  }

  // Helper: jet invariant mass

  template <typename TRACKS, typename CANDIDATES>
  float computeJetMass(TRACKS const& tracks,
                       CANDIDATES const& candidates,
                       double candMass = -1.)
  {
    std::array<double, 3> momTotal{0., 0., 0.};
    double energyTot = 0.;

    for (auto const& trk : tracks) {
      momTotal[0] += trk.px();
      momTotal[1] += trk.py();
      momTotal[2] += trk.pz();
      energyTot += RecoDecay::p(trk.px(), trk.py(), trk.pz()); // massless approximation for ordinary tracks
    }

    for (auto const& cand : candidates) {
      double m = 0;
      if (candMass > 0.) {
        m = candMass;
      } else if constexpr (requires { cand.m(); }) {
        m = cand.m(); // DATA / MCD: reconstructed invariant mass
      } else {
        m = 0.;
      }

      momTotal[0] += cand.px();
      momTotal[1] += cand.py();
      momTotal[2] += cand.pz();
      energyTot += RecoDecay::e(cand.px(), cand.py(), cand.pz(), m);
    }

    const double mass2 = RecoDecay::m2(momTotal, energyTot);
    return (mass2 > 0.) ? static_cast<float>(std::sqrt(mass2)) : 0.f;
  }

  // Process: collision QA (DATA)

  void processCollisions(aod::JetCollision const& collision,
                         aod::JetTracks const& tracks)
  {
    registry.fill(HIST("h_collisions"), getValFromBin(AllCollisions));

    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      registry.fill(HIST("h_collisions"), 3.5f); // rejected
      return;
    }
    registry.fill(HIST("h_collisions"), getValFromBin(Sel8ZCut));

    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
      registry.fill(HIST("h_track_pt"), track.pt());
      registry.fill(HIST("h_track_eta"), track.eta());
      registry.fill(HIST("h_track_phi"), track.phi());
    }
  }

  PROCESS_SWITCH(JetHFAngularityTask, processCollisions,
                 "Collision and track QA", true);

  // Process: DATA jet analysis

  void processDataChargedSubstructure(
    aod::JetCollision const& collision,
    soa::Filtered<D0DataJets> const& jets,
    D0CandidatesData const& /*candidates*/,
    aod::JetTracks const&)
  {
    // ---- Event selection ----
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }

    bool collisionHasD0Jet = false;

    for (const auto& jet : jets) {

      registry.fill(HIST("h_jet_counter"), getValFromBin(ChargedJets));
      registry.fill(HIST("h_jet_pt"), jet.pt());
      registry.fill(HIST("h_jet_eta"), jet.eta());
      registry.fill(HIST("h_jet_phi"), jet.phi());

      // ---- Jet-level quantities ----

      auto jetTracks = jet.template tracks_as<aod::JetTracks>();
      auto jetCandidates = jet.template candidates_as<D0CandidatesData>();

      const int nConst = static_cast<int>(jetTracks.size()) + static_cast<int>(jetCandidates.size());
      const float angularity = computeLambda(jet, jetTracks, jetCandidates, 1.f, 1.f); // λ_1^1
      const float girth = computeLambda(jet, jetTracks, jetCandidates, 2.f, 1.f);      // λ_2^1
      const float mjet = computeJetMass(jetTracks, jetCandidates);

      const std::array<double, 3> jetMom{jet.px(), jet.py(), jet.pz()};

      bool hasD0 = false;

      // ---- D0 candidate loop ----
      for (const auto& d0Candidate : jetCandidates) {

        hasD0 = true;

        const std::array<double, 3> d0Mom{d0Candidate.px(), d0Candidate.py(), d0Candidate.pz()};

        // Longitudinal momentum fraction: z_|| = (p_D0 . p_jet) / |p_jet|^2
        const float zParallel = RecoDecay::dotProd(d0Mom, jetMom) / RecoDecay::mag2(jetMom);
        // Angular separation between D0 and jet axis
        const float axisDistance = jetutilities::deltaR(jet, d0Candidate);

        // ---- Histograms ----
        registry.fill(HIST("h_d0_mass"), d0Candidate.m());
        registry.fill(HIST("h_d0_pt"), d0Candidate.pt());
        registry.fill(HIST("h_d0_eta"), d0Candidate.eta());
        registry.fill(HIST("h_d0_phi"), d0Candidate.phi());
        registry.fill(HIST("h_d0_y"), d0Candidate.y());
        registry.fill(HIST("h_d0_jet_projection"), zParallel);
        registry.fill(HIST("h_d0_jet_distance"), axisDistance);

        registry.fill(HIST("h_jet_nconst"), nConst);
        registry.fill(HIST("h_d0_mlscore0"), d0Candidate.mlScores()[0]);
        registry.fill(HIST("h_d0_mlscore1"), d0Candidate.mlScores()[1]);
        registry.fill(HIST("h_d0_mlscore2"), d0Candidate.mlScores()[2]);

        registry.fill(HIST("hSparse_d0"),
                      d0Candidate.m(), d0Candidate.pt(),
                      jet.pt(), zParallel, axisDistance,
                      angularity, girth);

        objJetTable(axisDistance,
                    jet.pt(),
                    jet.eta(),
                    jet.phi(),
                    nConst,
                    angularity,
                    girth,
                    mjet,
                    zParallel,
                    d0Candidate.pt(),
                    d0Candidate.eta(),
                    d0Candidate.phi(),
                    d0Candidate.m(),
                    d0Candidate.y(),
                    d0Candidate.mlScores()[0],
                    d0Candidate.mlScores()[1],
                    d0Candidate.mlScores()[2]);

      } // end D0 loop

      if (hasD0) {
        registry.fill(HIST("h_jet_counter"), getValFromBin(D0TaggedJets));
        registry.fill(HIST("h_d0_jet_lambda11"), angularity);
        registry.fill(HIST("h_d0_jet_lambda12"), girth);
        registry.fill(HIST("h_d0_jet_mass"), mjet);
        collisionHasD0Jet = true;
      }

    } // end jet loop

    if (collisionHasD0Jet) {
      registry.fill(HIST("h_collisions"), 2.5f); // has D0 jet
    }
  }

  PROCESS_SWITCH(JetHFAngularityTask, processDataChargedSubstructure,
                 "DATA: D0-tagged charged jet substructure", true);

  // Process: MC DETECTOR-LEVEL (MCD)

  void processMCDChargedSubstructure(
    aod::JetCollision const& collision,
    soa::Filtered<D0MCDJets> const& mcdjets,
    D0CandidatesMCD const& /*candidates*/,
    aod::JetTracks const&)
  {
    // ---- Event selection (same sel8 + |z| ) ----
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }

    registry.fill(HIST("h_collisions_mcd"), 0.5f); // selected collisions

    bool collisionHasD0Jet = false;

    for (const auto& jet : mcdjets) {

      registry.fill(HIST("h_jet_counter_mcd"), getValFromBin(ChargedJets));
      registry.fill(HIST("h_jet_pt_mcd"), jet.pt());
      registry.fill(HIST("h_jet_eta_mcd"), jet.eta());
      registry.fill(HIST("h_jet_phi_mcd"), jet.phi());

      // ---- Jet-level quantities (detector level) ----
      auto jetTracks = jet.template tracks_as<aod::JetTracks>();
      auto jetCandidates = jet.template candidates_as<D0CandidatesMCD>();

      const int nConst = static_cast<int>(jetTracks.size()) + static_cast<int>(jetCandidates.size());
      const float angularity = computeLambda(jet, jetTracks, jetCandidates, 1.f, 1.f); // λ_1^1
      const float girth = computeLambda(jet, jetTracks, jetCandidates, 2.f, 1.f);      // λ_2^1
      const float mjet = computeJetMass(jetTracks, jetCandidates);

      const std::array<double, 3> jetMom{jet.px(), jet.py(), jet.pz()};

      bool hasD0 = false;

      //  D0 candidate loop (with MC truth-matching info)
      for (const auto& d0Candidate : jetCandidates) {

        hasD0 = true;

        const std::array<double, 3> d0Mom{d0Candidate.px(), d0Candidate.py(), d0Candidate.pz()};

        const float zParallel = RecoDecay::dotProd(d0Mom, jetMom) / RecoDecay::mag2(jetMom);
        const float axisDistance = jetutilities::deltaR(jet, d0Candidate);

        const int8_t flagMcMatch = d0Candidate.flagMcMatchRec();
        const int8_t originMc = d0Candidate.originMcRec(); // 0 none, 1 prompt, 2 non-prompt

        // ---- Reflection identification ----
        const bool isReflection = isReflectedCandidate(d0Candidate);
        const bool isTrueSignal = (flagMcMatch != 0) && !isReflection;

        if (isReflection) {
          registry.fill(HIST("h_d0_mass_mcd_reflection"), d0Candidate.m());
          registry.fill(HIST("hSparse_d0_mcd_reflection"),
                        d0Candidate.m(), d0Candidate.pt(), jet.pt(), zParallel, originMc,
                        angularity, girth);
        } else if (isTrueSignal) {
          registry.fill(HIST("h_d0_mass_mcd_signal"), d0Candidate.m());
          registry.fill(HIST("hSparse_d0_mcd_calib"),
                        d0Candidate.m(), d0Candidate.pt(), jet.pt(), zParallel, originMc,
                        angularity, girth);
        }

        // ---- Histograms ----
        registry.fill(HIST("h_d0_mass_mcd"), d0Candidate.m());
        registry.fill(HIST("h_d0_pt_mcd"), d0Candidate.pt());
        registry.fill(HIST("h_d0_eta_mcd"), d0Candidate.eta());
        registry.fill(HIST("h_d0_phi_mcd"), d0Candidate.phi());
        registry.fill(HIST("h_d0_y_mcd"), d0Candidate.y());
        registry.fill(HIST("h_d0_jet_projection_mcd"), zParallel);
        registry.fill(HIST("h_d0_jet_distance_mcd"), axisDistance);
        registry.fill(HIST("h_d0_origin_mcd"), static_cast<float>(originMc));
        registry.fill(HIST("hSparse_d0_mcd"),
                      d0Candidate.m(), d0Candidate.pt(),
                      jet.pt(), zParallel, axisDistance,
                      angularity, girth);

        if (originMc == kOriginMcPrompt) {
          registry.fill(HIST("h_eff_run2_num_pr"), d0Candidate.pt());
        } else if (originMc == kOriginMcNonPrompt) {
          registry.fill(HIST("h_eff_run2_num_np"), d0Candidate.pt());
        }

        objJetMCDTable(axisDistance,
                       jet.pt(),
                       jet.eta(),
                       jet.phi(),
                       nConst,
                       angularity,
                       girth,
                       mjet,
                       zParallel,
                       d0Candidate.pt(),
                       d0Candidate.eta(),
                       d0Candidate.phi(),
                       d0Candidate.m(),
                       d0Candidate.y(),
                       d0Candidate.mlScores()[0],
                       d0Candidate.mlScores()[1],
                       d0Candidate.mlScores()[2],
                       flagMcMatch,
                       originMc);

      } // end D0 loop

      if (hasD0) {
        registry.fill(HIST("h_jet_counter_mcd"), getValFromBin(D0TaggedJets));
        registry.fill(HIST("h_d0_jet_lambda11_mcd"), angularity);
        registry.fill(HIST("h_d0_jet_lambda12_mcd"), girth);
        registry.fill(HIST("h_d0_jet_mass_mcd"), mjet);
        collisionHasD0Jet = true;
      }

    } // end jet loop

    if (collisionHasD0Jet) {
      registry.fill(HIST("h_collisions_mcd"), 1.5f); // has D0 jet
    }
  }

  PROCESS_SWITCH(JetHFAngularityTask, processMCDChargedSubstructure,
                 "MC DETECTOR LEVEL: D0-tagged charged jet substructure", false);

  // Process: MC PARTICLE-LEVEL (MCP)

  void processMCPChargedSubstructure(
    aod::JetMcCollision const& /*mcCollision*/,
    soa::Filtered<D0MCPJets> const& mcpjets,
    D0CandidatesMCP const& /*candidates*/,
    aod::JetParticles const&)
  {

    registry.fill(HIST("h_collisions_mcp"), 0.5f); // all MC collisions

    bool collisionHasD0Jet = false;

    for (const auto& jet : mcpjets) {

      registry.fill(HIST("h_jet_counter_mcp"), getValFromBin(ChargedJets));
      registry.fill(HIST("h_jet_pt_mcp"), jet.pt());
      registry.fill(HIST("h_jet_eta_mcp"), jet.eta());
      registry.fill(HIST("h_jet_phi_mcp"), jet.phi());

      //  Jet-level quantities (particle / generator level)
      auto jetParticles = jet.template tracks_as<aod::JetParticles>();
      auto jetCandidates = jet.template candidates_as<D0CandidatesMCP>();

      const int nConst = static_cast<int>(jetParticles.size()) + static_cast<int>(jetCandidates.size());
      const float angularity = computeLambda(jet, jetParticles, jetCandidates, 1.f, 1.f); // λ_1^1
      const float girth = computeLambda(jet, jetParticles, jetCandidates, 2.f, 1.f);      // λ_2^1
      const float mjet = computeJetMass(jetParticles, jetCandidates, massD0MCP);

      const std::array<double, 3> jetMom{jet.px(), jet.py(), jet.pz()};

      bool hasD0 = false;

      //  Generated D0 loop (with MC truth-matching info)
      for (const auto& d0Particle : jetCandidates) {

        hasD0 = true;

        const std::array<double, 3> d0Mom{d0Particle.px(), d0Particle.py(), d0Particle.pz()};

        const float zParallel = RecoDecay::dotProd(d0Mom, jetMom) / RecoDecay::mag2(jetMom);
        const float axisDistance = jetutilities::deltaR(jet, d0Particle);

        const int8_t flagMcMatch = d0Particle.flagMcMatchGen();
        const int8_t originMc = d0Particle.originMcGen(); // 0 none, 1 prompt, 2 non-prompt

        // ---- Histograms ----
        registry.fill(HIST("h_d0_pt_mcp"), d0Particle.pt());
        registry.fill(HIST("h_d0_eta_mcp"), d0Particle.eta());
        registry.fill(HIST("h_d0_phi_mcp"), d0Particle.phi());
        registry.fill(HIST("h_d0_y_mcp"), d0Particle.y());
        registry.fill(HIST("h_d0_jet_projection_mcp"), zParallel);
        registry.fill(HIST("h_d0_jet_distance_mcp"), axisDistance);
        registry.fill(HIST("h_d0_origin_mcp"), static_cast<float>(originMc));

        if (originMc == kOriginMcPrompt) {
          registry.fill(HIST("h_eff_run2_den_pr"), d0Particle.pt());
        } else if (originMc == kOriginMcNonPrompt) {
          registry.fill(HIST("h_eff_run2_den_np"), d0Particle.pt());
        }
        objJetMCPTable(axisDistance,
                       jet.pt(),
                       jet.eta(),
                       jet.phi(),
                       nConst,
                       angularity,
                       girth,
                       mjet,
                       zParallel,
                       d0Particle.pt(),
                       d0Particle.eta(),
                       d0Particle.phi(),
                       d0Particle.y(),
                       flagMcMatch,
                       originMc);

      } // end D0 loop

      if (hasD0) {
        registry.fill(HIST("h_jet_counter_mcp"), getValFromBin(D0TaggedJets));
        registry.fill(HIST("h_d0_jet_lambda11_mcp"), angularity);
        registry.fill(HIST("h_d0_jet_lambda12_mcp"), girth);
        registry.fill(HIST("h_d0_jet_mass_mcp"), mjet);
        collisionHasD0Jet = true;
      }

    } // end jet loop

    if (collisionHasD0Jet) {
      registry.fill(HIST("h_collisions_mcp"), 1.5f); // has D0 jet
    }
  }

  PROCESS_SWITCH(JetHFAngularityTask, processMCPChargedSubstructure,
                 "MC PARTICLE LEVEL: D0-tagged charged jet substructure", false);

  void processMCDChargedSubstructureMatched(aod::JetMcCollisions const& mccollisions,
                                            aod::JetCollisionsMCD const& collisions,
                                            D0MCDJetsMatched const& mcdjets,
                                            D0MCPJetsMatched const& mcpjets,
                                            D0CandidatesMCD const&,
                                            D0CandidatesMCP const&,
                                            aod::JetTracks const&,
                                            aod::JetParticles const&)
  {

    const int jetRBin = static_cast<int>(std::lround(static_cast<float>(jetR) * 100.0f));

    for (const auto& mccollision : mccollisions) {

      registry.fill(HIST("MC_CollisionCounter"), getValFromBin(BinMCColCntr::AllMCCollisions));

      if (std::abs(mccollision.posZ()) > vertexZCut) {
        continue;
      }
      registry.fill(HIST("MC_CollisionCounter"), getValFromBin(BinMCColCntr::SelectedMCCollisions));

      // const float mcWeight = mccollision.weight();

      const auto collisionsPerMCCollision = collisions.sliceBy(collisionsPerMCCollisionPreslice, mccollision.globalIndex());
      for (const auto& collision : collisionsPerMCCollision) {

        registry.fill(HIST("MC_CollisionCounter"), getValFromBin(BinMCColCntr::AssociatedRecoCollisions));

        if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits) ||
            !(std::abs(collision.posZ()) < vertexZCut)) {
          continue;
        }
        registry.fill(HIST("MC_CollisionCounter"), getValFromBin(BinMCColCntr::SelectedAssociatedRecoCollisions));

        const auto d0mcdJetsPerCollision = mcdjets.sliceBy(d0MCDJetsPerCollisionPreslice, collision.globalIndex());
        for (const auto& mcdjet : d0mcdJetsPerCollision) {
          if (mcdjet.pt() <= jetPtMin || mcdjet.r() != jetRBin) {
            continue;
          }
          registry.fill(HIST("MC_JetCounter"), getValFromBin(BinMCJetCntr::DetectorLevelJetInMCCollision));

          auto jetTracks = mcdjet.template tracks_as<aod::JetTracks>();
          auto jetCandidates = mcdjet.template candidates_as<D0CandidatesMCD>();

          const int nConst = static_cast<int>(jetTracks.size()) + static_cast<int>(jetCandidates.size());
          const float angularity = computeLambda(mcdjet, jetTracks, jetCandidates, 1.f, 1.f); // λ_1^1
          const float girth = computeLambda(mcdjet, jetTracks, jetCandidates, 2.f, 1.f);      // λ_2^1
          const float mjet = computeJetMass(jetTracks, jetCandidates);

          const std::array<double, 3> jetMomMCD{mcdjet.px(), mcdjet.py(), mcdjet.pz()};

          bool isGeoMatched = false;
          float matchedPt = -1.f;
          float matchedEta = -999.f;
          float matchedPhi = -999.f;
          float matchedDR = -1.f;
          float matchedAngularity = -1.f;
          float matchedGirth = -1.f;
          float matchedD0Pt = -1.f; // NEW
          int geoMatchedGlobalIndex = -1;
          int nGeoMatches = 0;

          if (mcdjet.has_matchedJetGeo()) {
            for (const auto& mcpjet : mcdjet.template matchedJetGeo_as<D0MCPJetsMatched>()) {
              ++nGeoMatches;

              const float dR = jetutilities::deltaR(mcdjet, mcpjet);
              registry.fill(HIST("h_jet_matching_dr_mcd_allcand"), dR);
              registry.fill(HIST("h_jet_pt_response_matrix_allcand"), mcdjet.pt(), mcpjet.pt());

              if (!isGeoMatched) {
                isGeoMatched = true;
                matchedPt = mcpjet.pt();
                matchedEta = mcpjet.eta();
                matchedPhi = mcpjet.phi();
                matchedDR = dR;
                geoMatchedGlobalIndex = mcpjet.globalIndex();

                auto mcpjetParticles = mcpjet.template tracks_as<aod::JetParticles>();
                auto mcpjetCandidates = mcpjet.template candidates_as<D0CandidatesMCP>();
                matchedAngularity = computeLambda(mcpjet, mcpjetParticles, mcpjetCandidates, 1.f, 1.f); // λ_1^1
                matchedGirth = computeLambda(mcpjet, mcpjetParticles, mcpjetCandidates, 2.f, 1.f);      // λ_2^1
                if (!mcpjetCandidates.empty()) {
                  auto matchedMcpD0 = mcpjet.template candidates_first_as<D0CandidatesMCP>();
                  matchedD0Pt = matchedMcpD0.pt();
                }
              }
            }
          }
          registry.fill(HIST("h_jet_matching_ngeo_mcd"), nGeoMatches);

          bool isCandMatched = false;
          int candMatchedGlobalIndex = -1;
          int nCandMatches = 0;

          if (mcdjet.has_matchedJetCand()) {
            for (const auto& mcpjet : mcdjet.template matchedJetCand_as<D0MCPJetsMatched>()) {
              ++nCandMatches;
              if (!isCandMatched) {
                isCandMatched = true;
                candMatchedGlobalIndex = mcpjet.globalIndex();
              }
            }
          }
          registry.fill(HIST("h_jet_matching_ncand_mcd"), nCandMatches);

          // "Clean" match: both criteria fire AND agree on the same target jet.
          const bool isCleanMatched = isGeoMatched && isCandMatched &&
                                      (geoMatchedGlobalIndex == candMatchedGlobalIndex);

          registry.fill(HIST("h_jet_matching_geo_status_mcd"), isGeoMatched ? 1.5f : 0.5f);
          registry.fill(HIST("h_jet_matching_cand_status_mcd"), isCandMatched ? 1.5f : 0.5f);
          registry.fill(HIST("h_jet_matching_clean_mcd"), isCleanMatched ? 1.5f : 0.5f);

          if (isGeoMatched) {
            const bool agree = isCandMatched && (geoMatchedGlobalIndex == candMatchedGlobalIndex);
            registry.fill(HIST("h_jet_matching_geocand_disagree_mcd"), agree ? 1.5f : 0.5f);
          }

          if (isCleanMatched) {
            registry.fill(HIST("h_jet_pt_response_matrix"), mcdjet.pt(), matchedPt);
            registry.fill(HIST("h_jet_matching_dr_mcd"), matchedDR);
          }

          const auto geoStatus = static_cast<int8_t>(isGeoMatched);
          const auto candStatus = static_cast<int8_t>(isCandMatched);
          const auto cleanStatus = static_cast<int8_t>(isCleanMatched);

          const int isCleanMatchedMCD = isCleanMatched ? 1 : 0;

          if (isCleanMatched) {
            registry.fill(HIST("MC_JetCounter"), getValFromBin(BinMCJetCntr::DetectorLevelJetWithMatchedCandidate));
          }

          for (const auto& mcdD0cand : jetCandidates) {

            const std::array<double, 3> candMomMCD{mcdD0cand.px(), mcdD0cand.py(), mcdD0cand.pz()};
            const float zParallelMCD = RecoDecay::dotProd(candMomMCD, jetMomMCD) / RecoDecay::mag2(jetMomMCD);
            const float axisDistance = jetutilities::deltaR(mcdjet, mcdD0cand);

            const int8_t flagMcMatch = mcdD0cand.flagMcMatchRec();
            const int originMcd = mcdD0cand.originMcRec(); // 0 none/bkg, 1 prompt, 2 non-prompt

            // Reflection identification: candidateSelFlag() hypothesis vs. true species
            const bool isReflection = isReflectedCandidate(mcdD0cand);
            const int isReflectionFlag = isReflection ? 1 : 0;
            const bool isTrueSignal = (flagMcMatch != 0) && !isReflection;

            if (isReflection) {
              registry.fill(HIST("h_d0_mass_mcd_reflection"), mcdD0cand.m());
              registry.fill(HIST("hSparse_d0_mcd_reflection"),
                            mcdD0cand.m(), mcdD0cand.pt(), mcdjet.pt(), zParallelMCD, originMcd,
                            angularity, girth);
            } else if (isTrueSignal) {
              registry.fill(HIST("h_d0_mass_mcd_signal"), mcdD0cand.m());
              registry.fill(HIST("hSparse_d0_mcd_calib"),
                            mcdD0cand.m(), mcdD0cand.pt(), mcdjet.pt(), zParallelMCD, originMcd,
                            angularity, girth);
            }

            registry.fill(HIST("hSparse_d0_mcd_eff"),
                          mcdD0cand.m(),
                          mcdD0cand.pt(),
                          mcdjet.pt(),
                          zParallelMCD,
                          originMcd,
                          isCleanMatchedMCD,
                          angularity,
                          girth,
                          isReflectionFlag);

            if (originMcd == kOriginMcPrompt) {
              registry.fill(HIST("h_eff_run2_num_pr"), mcdD0cand.pt());
              registry.fill(HIST("h_kineff_det_den_pr"), mcdjet.pt(), mcdD0cand.pt());
              if (isGeoMatched) {
                registry.fill(HIST("h_kineff_det_num_pr"), mcdjet.pt(), mcdD0cand.pt());
              }
            } else if (originMcd == kOriginMcNonPrompt) {
              registry.fill(HIST("h_eff_run2_num_np"), mcdD0cand.pt());
              registry.fill(HIST("h_kineff_det_den_np"), mcdjet.pt(), mcdD0cand.pt());
              if (isGeoMatched) {
                registry.fill(HIST("h_kineff_det_num_np"), mcdjet.pt(), mcdD0cand.pt());
              }
            }

            if (isCleanMatched) {
              if (originMcd == kOriginMcPrompt) {
                registry.fill(HIST("h_jet_pt_response_matrix_pr"), mcdjet.pt(), matchedPt);
                registry.fill(HIST("h_response_angularity_pr"), mcdjet.pt(), angularity, matchedPt, matchedAngularity);
                registry.fill(HIST("h_response_girth_pr"), mcdjet.pt(), girth, matchedPt, matchedGirth);
                registry.fill(HIST("h_response_d0pt_pr"), mcdjet.pt(), mcdD0cand.pt(), matchedPt, matchedD0Pt);
              } else if (originMcd == kOriginMcNonPrompt) {
                registry.fill(HIST("h_jet_pt_response_matrix_np"), mcdjet.pt(), matchedPt);
                registry.fill(HIST("h_response_angularity_np"), mcdjet.pt(), angularity, matchedPt, matchedAngularity);
                registry.fill(HIST("h_response_girth_np"), mcdjet.pt(), girth, matchedPt, matchedGirth);
                registry.fill(HIST("h_response_d0pt_np"), mcdjet.pt(), mcdD0cand.pt(), matchedPt, matchedD0Pt);
              }
            }

            objJetMCDMatchedTable(axisDistance,
                                  mcdjet.pt(),
                                  mcdjet.eta(),
                                  mcdjet.phi(),
                                  nConst,
                                  angularity,
                                  girth,
                                  mjet,
                                  zParallelMCD,
                                  mcdD0cand.pt(),
                                  mcdD0cand.eta(),
                                  mcdD0cand.phi(),
                                  mcdD0cand.m(),
                                  mcdD0cand.y(),
                                  mcdD0cand.mlScores()[0],
                                  mcdD0cand.mlScores()[1],
                                  mcdD0cand.mlScores()[2],
                                  flagMcMatch,
                                  originMcd,
                                  geoStatus,
                                  candStatus,
                                  cleanStatus,
                                  matchedPt,
                                  matchedEta,
                                  matchedPhi,
                                  matchedDR,
                                  matchedD0Pt);
          } // end D0 candidate loop
        } // end MCD jet loop
      } // end reconstructed-collision loop

      const auto d0mcpJetsPerMCCollision = mcpjets.sliceBy(d0MCPJetsPerMCCollisionPreslice, mccollision.globalIndex());
      for (const auto& mcpjet : d0mcpJetsPerMCCollision) {

        if (mcpjet.pt() <= jetPtMin || mcpjet.r() != jetRBin) {
          continue;
        }

        registry.fill(HIST("MC_JetCounter"), getValFromBin(BinMCJetCntr::ParticleLevelJetInMCCollision));

        auto mcpD0cand = mcpjet.template candidates_first_as<D0CandidatesMCP>();

        const int originMcp = mcpD0cand.originMcGen(); // 0 none/bkg, 1 prompt, 2 non-prompt

        bool mcpIsCleanMatched = false;
        if (mcpjet.has_matchedJetGeo() && mcpjet.has_matchedJetCand()) {
          int mcpGeoTarget = -1;
          for (const auto& mcdjetGeo : mcpjet.template matchedJetGeo_as<D0MCDJetsMatched>()) {
            mcpGeoTarget = mcdjetGeo.globalIndex();
            break; // first geo match, mirrors the MCD-side "first match wins" convention
          }
          int mcpCandTarget = -1;
          for (const auto& mcdjetCand : mcpjet.template matchedJetCand_as<D0MCDJetsMatched>()) {
            mcpCandTarget = mcdjetCand.globalIndex();
            break;
          }
          mcpIsCleanMatched = (mcpGeoTarget == mcpCandTarget) && (mcpGeoTarget != -1);
        }
        const int isCleanMatchedMCP = mcpIsCleanMatched ? 1 : 0;

        if (mcpIsCleanMatched) {
          registry.fill(HIST("MC_JetCounter"), getValFromBin(BinMCJetCntr::ParticleLevelJetWithMatchedCandidate));
        }

        const std::array<double, 3> jetMomMCP{mcpjet.px(), mcpjet.py(), mcpjet.pz()};
        const std::array<double, 3> candMomMCP{mcpD0cand.px(), mcpD0cand.py(), mcpD0cand.pz()};
        const float zParallelMCP = RecoDecay::dotProd(candMomMCP, jetMomMCP) / RecoDecay::mag2(jetMomMCP);

        auto mcpJetParticles = mcpjet.template tracks_as<aod::JetParticles>();
        auto mcpJetCandidates = mcpjet.template candidates_as<D0CandidatesMCP>();
        const float angularityMCP = computeLambda(mcpjet, mcpJetParticles, mcpJetCandidates, 1.f, 1.f); // λ_1^1
        const float girthMCP = computeLambda(mcpjet, mcpJetParticles, mcpJetCandidates, 2.f, 1.f);      // λ_2^1

        registry.fill(HIST("hSparse_d0_mcp_eff"),
                      mcpD0cand.pt(),
                      mcpjet.pt(),
                      zParallelMCP,
                      originMcp,
                      isCleanMatchedMCP,
                      angularityMCP,
                      girthMCP);

        const bool mcpIsGeoMatched = mcpjet.has_matchedJetGeo();
        if (originMcp == kOriginMcPrompt) {
          registry.fill(HIST("h_eff_run2_den_pr"), mcpD0cand.pt());
          registry.fill(HIST("h_kineff_gen_den_pr"), mcpjet.pt(), mcpD0cand.pt());
          if (mcpIsGeoMatched) {
            registry.fill(HIST("h_kineff_gen_num_pr"), mcpjet.pt(), mcpD0cand.pt());
          }
        } else if (originMcp == kOriginMcNonPrompt) {
          registry.fill(HIST("h_eff_run2_den_np"), mcpD0cand.pt());
          registry.fill(HIST("h_kineff_gen_den_np"), mcpjet.pt(), mcpD0cand.pt());
          if (mcpIsGeoMatched) {
            registry.fill(HIST("h_kineff_gen_num_np"), mcpjet.pt(), mcpD0cand.pt());
          }
        }
      } // end MCP jet loop
    } // end MC collision loop
  }

  PROCESS_SWITCH(JetHFAngularityTask, processMCDChargedSubstructureMatched,
                 "MC: unified efficiency + jet matching + response matrix + reflection (needs jet-matching workflow)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetHFAngularityTask>(cfgc)};
}
