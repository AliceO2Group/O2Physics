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

// QA task for MC-based outliers
//
/// \author Jaime Norman <jaime.norman@cern.ch>

#include "RecoDecay.h"

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetSubtraction.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include <Framework/Configurable.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TMath.h>
#include <TMathBase.h>
#include <TRandom3.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <type_traits>
#include <vector>

#include <math.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetOutlierQATask {

  HistogramRegistry registry;

  using JetParticlesWithOriginal = soa::Join<aod::JetParticles, aod::JMcParticlePIs>;

  Preslice<aod::JetTracks> perColTrack = aod::jtrack::collisionId;
  Preslice<soa::Join<aod::JetMcCollisions, aod::JMcCollisionPIs>> perColParticle = aod::jmccollision::mcCollisionId;
  Preslice<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetEventWeights>> perColJets = aod::jet::collisionId;
  Preslice<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets, aod::ChargedMCDetectorLevelJetEventWeights>> perColJetsMatched = aod::jet::collisionId;

  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> centralityMin{"centralityMin", -999.0, "minimum centrality"};
  Configurable<float> centralityMax{"centralityMax", 999.0, "maximum centrality"};
  Configurable<std::vector<double>> jetRadii{"jetRadii", std::vector<double>{0.4}, "jet resolution parameters"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum eta acceptance for tracks"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum eta acceptance for tracks"};
  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum pT acceptance for tracks"};
  Configurable<float> trackPtMax{"trackPtMax", 100.0, "maximum pT acceptance for tracks"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
  Configurable<float> pTHatMaxMCD{"pTHatMaxMCD", 999.0, "maximum fraction of hard scattering for jet acceptance in detector MC"};
  Configurable<float> pTHatMaxMCP{"pTHatMaxMCP", 999.0, "maximum fraction of hard scattering for jet acceptance in particle MC"};
  Configurable<float> pTHatExponent{"pTHatExponent", 6.0, "exponent of the event weight for the calculation of pTHat"};
  Configurable<float> pTHatMaxMCDOutlier{"pTHatMaxMCDOutlier", 3.0, "maximum fraction of hard scattering for jet acceptance in detector MC for outlier studies"};
  Configurable<double> jetPtMax{"jetPtMax", 200., "set jet pT bin max"};
  Configurable<float> jetEtaMin{"jetEtaMin", -99.0, "minimum jet pseudorapidity"};
  Configurable<float> jetEtaMax{"jetEtaMax", 99.0, "maximum jet pseudorapidity"};
  Configurable<float> jetAreaFractionMin{"jetAreaFractionMin", -99.0, "used to make a cut on the jet areas"};
  Configurable<float> leadingConstituentPtMin{"leadingConstituentPtMin", -99.0, "minimum pT selection on jet constituent"};
  Configurable<float> leadingConstituentPtMax{"leadingConstituentPtMax", 9999.0, "maximum pT selection on jet constituent"};
  Configurable<int> trackOccupancyInTimeRangeMax{"trackOccupancyInTimeRangeMax", 999999, "maximum occupancy of tracks in neighbouring collisions in a given time range; only applied to reconstructed collisions (data and mcd jets), not mc collisions (mcp jets)"};
  Configurable<int> trackOccupancyInTimeRangeMin{"trackOccupancyInTimeRangeMin", -999999, "minimum occupancy of tracks in neighbouring collisions in a given time range; only applied to reconstructed collisions (data and mcd jets), not mc collisions (mcp jets)"};
  Configurable<float> splitCollisionsDeltaZ{"splitCollisionsDeltaZ", 0.2, "threshold in delta z to assign as split collision"};
  Configurable<float> splitCollisionsDeltaZPart{"splitCollisionsDeltaZPart", 0.2, "threshold in delta z to assign as split collision particle level"};
  Configurable<unsigned int> splitCollisionsDeltaBC{"splitCollisionsDeltaBC", 5, "threshold in BC to assign as split collision"};
  Configurable<int> mergeCollisionsDeltaMin{"mergeCollisionsDeltaMin", -10, "number of prior collisions to search for close Z position"};
  Configurable<int> mergeCollisionsDeltaMax{"mergeCollisionsDeltaMax", 10, "number of following collisions to search for close Z position"};
  Configurable<int> maxNTracksJJdifferent{"maxNTracksJJdifferent", 10, "maximum number of tracks from different JJ collision to be considered for track rejection"};

  std::map<uint64_t, std::vector<int64_t>> fBCCollMap; // key: global BC, value: vector of reduced event global indices

  std::vector<int> eventSelectionBits;
  int trackSelection = -1;

  std::vector<double> jetPtBins;
  std::vector<double> jetPtBinsRhoAreaSub;

  void init(o2::framework::InitContext&)
  {
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

    auto jetPtTemp = 0.0;
    jetPtBins.push_back(jetPtTemp);
    while (jetPtTemp < jetPtMax) {
      if (jetPtTemp < 100.0) {
        jetPtTemp += 1.0;
        jetPtBins.push_back(jetPtTemp);
      } else if (jetPtTemp < 200.0) {
        jetPtTemp += 5.0;
        jetPtBins.push_back(jetPtTemp);
      } else {
        jetPtTemp += 10.0;
        jetPtBins.push_back(jetPtTemp);
      }
    }
    std::vector<double> pThatBinning = {0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0,
                                        3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 20.0, 30.0};

    AxisSpec jetPtAxis = {jetPtBins, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec pThatAxis = {pThatBinning, "p_{T} / #hat{p_{T}}"};

    if (doprocessJetsAmbiguous) {
      // outliers
      registry.add("h3_pthat_jet_pt_jet_ntracks_outliers", "#it#hat{p} (GeV/#it{c});#it{p}_{T,jet} (GeV/#it{c});N_{jet tracks}", {HistType::kTH3F, {{100, 0, 200}, jetPtAxis, {200, -0.5, 199.5}}});
      registry.add("h3_jet_pt_track_dEta_track_dPhi_outliers", "#it{p}_{T,jet} (GeV/#it{c});#Delta#eta;#Delta#phi", {HistType::kTH3F, {jetPtAxis, {100, -0.5, 0.5}, {100, -0.5, 0.5}}});
      registry.add("h3_pthat_jet_pt_jet_ntracks_all", "it#hat{p} (GeV/#it{c});#it{p}_{T,jet} (GeV/#it{c});N_{jet tracks}", {HistType::kTH3F, {{100, 0, 200}, jetPtAxis, {200, -0.5, 199.5}}});
      registry.add("h3_jet_pt_track_dEta_track_dPhi_all", "#it{p}_{T,jet} (GeV/#it{c});#Delta#eta;#Delta#phi", {HistType::kTH3F, {jetPtAxis, {100, -0.5, 0.5}, {100, -0.5, 0.5}}});

      registry.add("h3_pthat_jet_pt_jet_ntracks_outliers_noweight", "#it#hat{p} (GeV/#it{c});#it{p}_{T,jet} (GeV/#it{c});N_{jet tracks}", {HistType::kTH3F, {{100, 0, 200}, jetPtAxis, {200, -0.5, 199.5}}});
      registry.add("h3_jet_pt_track_dEta_track_dPhi_outliers_noweight", "#it{p}_{T,jet} (GeV/#it{c});#Delta#eta;#Delta#phi", {HistType::kTH3F, {jetPtAxis, {100, -0.5, 0.5}, {100, -0.5, 0.5}}});
      registry.add("h2_jet_pt_Dz_outliers_noweight", "#it{p}_{T,jet} (GeV/#it{c});z = #it{p}_{T,track} / #it{p}_{T,jet}", {HistType::kTH2F, {jetPtAxis, {30, 0, 1}}});
      registry.add("h3_pthat_jet_pt_jet_ntracks_all_noweight", "it#hat{p} (GeV/#it{c});#it{p}_{T,jet} (GeV/#it{c});N_{jet tracks}", {HistType::kTH3F, {{100, 0, 200}, jetPtAxis, {200, -0.5, 199.5}}});
      registry.add("h3_jet_pt_track_dEta_track_dPhi_all_noweight", "#it{p}_{T,jet} (GeV/#it{c});#Delta#eta;#Delta#phi", {HistType::kTH3F, {jetPtAxis, {100, -0.5, 0.5}, {100, -0.5, 0.5}}});
      registry.add("h2_jet_pt_Dz_all_noweight", "#it{p}_{T,jet} (GeV/#it{c});z = #it{p}_{T,track} / #it{p}_{T,jet}", {HistType::kTH2F, {jetPtAxis, {30, 0, 1}}});

      registry.add("h3_jet_pt_track_pt_pt_hat_ambiguous", "ambiguous;#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,track} (GeV/#it{c}); #hat{#it{p}_{T}} (GeV/#it{c})", {HistType::kTH3F, {{150, 0, 300}, {100, 0, 100}, {200, 0, 600}}});
      registry.add("h3_jet_pt_track_pt_pt_hat_unambiguous", "matched;#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,track} (GeV/#it{c}); #hat{#it{p}_{T}} (GeV/#it{c})", {HistType::kTH3F, {{150, 0, 300}, {100, 0, 100}, {200, 0, 600}}});
      registry.add("h3_jet_pt_frac_pt_ambiguous_pt_hat", "fraction pT;#it{p}_{T,jet} (GeV/#it{c});fraction of #it{p}_{T,track} unmatched; #hat{#it{p}_{T}} (GeV/#it{c})", {HistType::kTH3F, {{150, 0, 300}, {40, 0, 1.1}, {200, 0, 600}}});
      registry.add("h3_jet_pt_frac_constituents_ambiguous_pt_hat", "fraction const;#it{p}_{T,jet} (GeV/#it{c});fraction of constituents matched; #hat{#it{p}_{T}} (GeV/#it{c})", {HistType::kTH3F, {{150, 0, 300}, {40, 0, 1.1}, {200, 0, 600}}});
      registry.add("h3_jet_pt_track_pt_pt_hat_no_particle", "no matching particle;#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,track} (GeV/#it{c}); #hat{#it{p}_{T}} (GeV/#it{c})", {HistType::kTH3F, {{150, 0, 300}, {100, 0, 100}, {200, 0, 600}}});
      registry.add("h3_jet_pt_track_pt_pt_hat_with_particle", "with matching particle;#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,track} (GeV/#it{c}); #hat{#it{p}_{T}} (GeV/#it{c})", {HistType::kTH3F, {{150, 0, 300}, {100, 0, 100}, {200, 0, 600}}});
      registry.add("h3_jet_pt_frac_pt_unmatched_particle_pt_hat", "fraction pT;#it{p}_{T,jet} (GeV/#it{c});fraction of #it{p}_{T,track} unmatched; #hat{#it{p}_{T}} (GeV/#it{c})", {HistType::kTH3F, {{150, 0, 300}, {40, 0, 1.1}, {200, 0, 600}}});
      registry.add("h3_jet_pt_frac_constituents_unmatched_particle_pt_hat", "fraction const;#it{p}_{T,jet} (GeV/#it{c});fraction of constituents unmatched; #hat{#it{p}_{T}} (GeV/#it{c})", {HistType::kTH3F, {{150, 0, 300}, {40, 0, 1.1}, {200, 0, 600}}});

      // collision-level properties
      registry.add("h_collision_multiplicity_outlier", "multiplicity", {HistType::kTH1F, {{200, 0, 200}}});
      registry.add("h_collision_pTHat_multiplicity_outlier", "multiplicity", {HistType::kTH2F, {{100, 0, 500}, {200, 0, 200}}});
      registry.add("h_collision_multiplicity_all", "multiplicity", {HistType::kTH1F, {{200, 0, 200}}});
      registry.add("h_collision_pTHat_multiplicity_all", "multiplicity", {HistType::kTH2F, {{100, 0, 500}, {200, 0, 200}}});
      registry.add("h_collision_trackOccupancyInTimeRange_outlier", "track occupancy", {HistType::kTH1F, {{1000, 0, 10000}}});
      registry.add("h_collision_trackOccupancyInTimeRange_all", "track occupancy", {HistType::kTH1F, {{1000, 0, 10000}}});
      registry.add("h_collision_multFV0A_outlier", "mult V0A", {HistType::kTH1F, {{1000, 0, 1000}}});
      registry.add("h_collision_multFV0A_all", "mult V0A", {HistType::kTH1F, {{1000, 0, 1000}}});
      registry.add("h_collision_multFV0C_outlier", "mult V0A", {HistType::kTH1F, {{1000, 0, 1000}}});
      registry.add("h_collision_multFV0C_all", "mult V0A", {HistType::kTH1F, {{1000, 0, 1000}}});
      registry.add("h_collision_multFV0M_outlier", "mult V0A", {HistType::kTH1F, {{1000, 0, 1000}}});
      registry.add("h_collision_multFV0M_all", "mult V0A", {HistType::kTH1F, {{1000, 0, 1000}}});
    }
    if (doprocessCollisionsBC) {
      // delta Z checks
      registry.add("h_DeltaZ_InBunch", "Delta Z between two events in bunch", {HistType::kTH1F, {{1200, -30, 30}}});
      registry.add("h_DeltaZ_Z1_InBunch", "Delta Z between two events in bunch vs Z1", {HistType::kTH2F, {{1200, 30, 30}, {400, -12, 12}}});
      registry.add("h_Z1_Z2_InBunch", "Z1 vs Z2 between two events in bunch", {HistType::kTH2F, {{400, -12, 12}, {400, -12, 12}}});
      registry.add("h_DeltaZ_OutOfBunch", "Delta Z between two events out of bunch", {HistType::kTH1F, {{1200, -30, 30}}});
      registry.add("h_DeltaZ_Z1_OutOfBunch", "Delta Z between two events out of bunch vs Z1", {HistType::kTH2F, {{1200, 30, 30}, {400, -12, 12}}});
      registry.add("h_Z1_Z2_OutOfBunch", "Z1 vs Z2 between two events out of bunch", {HistType::kTH2F, {{400, -12, 12}, {400, -12, 12}}});
      registry.add("h_DeltaZ_InBunch_JJ", "Delta Z between two events in bunch", {HistType::kTH1F, {{1200, -30, 30}}});
      registry.add("h_DeltaZ_Z1_InBunch_JJ", "Delta Z between two events in bunch vs Z1", {HistType::kTH2F, {{1200, 30, 30}, {400, -12, 12}}});
      registry.add("h_Z1_Z2_InBunch_JJ", "Z1 vs Z2 between two events in bunch", {HistType::kTH2F, {{400, -12, 12}, {400, -12, 12}}});
      registry.add("h_DeltaZ_OutOfBunch_JJ", "Delta Z between two events out of bunch", {HistType::kTH1F, {{1200, -30, 30}}});
      registry.add("h_DeltaZ_Z1_OutOfBunch_JJ", "Delta Z between two events out of bunch vs Z1", {HistType::kTH2F, {{1200, 30, 30}, {400, -12, 12}}});
      registry.add("h_Z1_Z2_OutOfBunch_JJ", "Z1 vs Z2 between two events out of bunch", {HistType::kTH2F, {{400, -12, 12}, {400, -12, 12}}});
      registry.add("h_Z", "Delta Z between two events", {HistType::kTH1F, {{400, -12, 12}}});
    }
    if (doprocessTracksBC) {
      // track checks based on z position
      registry.add("h_Z_resolution", "Z resolution", {HistType::kTH1F, {{200, -0.05, 0.05}}});
      registry.add("h_Z_resolution_wide", "Z resolution", {HistType::kTH1F, {{200, -1, 1}}});
      registry.add("h_Z_reco_rejected", "Z position reconstructed rejected", {HistType::kTH1F, {{400, -12, 12}}});
      registry.add("h_Z_true_rejected", "Z position particle rejected", {HistType::kTH1F, {{400, -12, 12}}});
      registry.add("h_Z_reco_accepted", "Z position reconstructed accepted", {HistType::kTH1F, {{400, -12, 12}}});
      registry.add("h_Z_true_accepted", "Z position particle accepted", {HistType::kTH1F, {{400, -12, 12}}});

      registry.add("h_track_pt", "track pt;p_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0, 200}}});
      registry.add("h_track_eta", "track eta;#eta_{track};entries", {HistType::kTH1F, {{100, -5, 5}}});
      registry.add("h_track_phi", "track phi;#varphi_{track} (rad);entries", {HistType::kTH1F, {{160, -1.0, 7.0}}});
      registry.add("h_track_pt_eta", "track pt vs eta;p_{T,track} (GeV/#it{c});#eta_{track};entries", {HistType::kTH2F, {{200, 0, 200}, {100, -5, 5}}});
      registry.add("h_track_pt_phi", "track pt vs phi;p_{T,track} (GeV/#it{c});#varphi_{track} (rad);entries", {HistType::kTH2F, {{200, 0, 200}, {160, -1.0, 7.0}}});
      registry.add("h_pt_hard_track_pt", "Tracks vs pThard;#frac{p_{T}}{#hat{p}};p_{T}", {HistType::kTH2F, {pThatAxis, {200, 0, 200}}});

      registry.add("h_track_pt_accepted", "track pt accepted;p_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0, 200}}});
      registry.add("h_track_eta_accepted", "track eta accepted;#eta_{track};entries", {HistType::kTH1F, {{100, -5, 5}}});
      registry.add("h_track_phi_accepted", "track phi accepted;#varphi_{track} (rad);entries", {HistType::kTH1F, {{160, -1.0, 7.0}}});
      registry.add("h_track_pt_eta_accepted", "track pt vs eta accepted;p_{T,track} (GeV/#it{c});#eta_{track};entries", {HistType::kTH2F, {{200, 0, 200}, {100, -5, 5}}});
      registry.add("h_track_pt_phi_accepted", "track pt vs phi accepted;p_{T,track} (GeV/#it{c});#varphi_{track} (rad);entries", {HistType::kTH2F, {{200, 0, 200}, {160, -1.0, 7.0}}});
      registry.add("h_pt_hard_track_pt_accepted", "Tracks vs pThard;#frac{p_{T}}{#hat{p}};p_{T}", {HistType::kTH2F, {pThatAxis, {200, 0, 200}}});
      registry.add("h_track_pt_accepted_no_JJ_outlier", "track pT with no JJ outlier", {HistType::kTH1F, {{200, 0, 200}}});
      registry.add("h_track_pt_eta_accepted_no_JJ_outlier", "track pT vs eta with no JJ outlier;p_{T,track} (GeV/#it{c});#eta_{track};entries", {HistType::kTH2F, {{200, 0, 200}, {100, -5, 5}}});
      registry.add("h_track_pt_phi_accepted_no_JJ_outlier", "track pT vs phi with no JJ outlier;p_{T,track} (GeV/#it{c});#varphi_{track} (rad);entries", {HistType::kTH2F, {{200, 0, 200}, {160, -1.0, 7.0}}});
      registry.add("h_pt_hard_track_pt_accepted_no_JJ_outlier", "Tracks vs pThard;#frac{p_{T}}{#hat{p}};p_{T}", {HistType::kTH2F, {pThatAxis, {200, 0, 200}}});
      registry.add("h_track_pt_with_JJ_outlier", "track pT with JJ outlier", {HistType::kTH1F, {{200, 0, 200}}});
      registry.add("h_track_pt_eta_with_JJ_outlier", "track pT vs eta with JJ outlier;p_{T,track} (GeV/#it{c});#eta_{track};entries", {HistType::kTH2F, {{200, 0, 200}, {100, -5, 5}}});
      registry.add("h_track_pt_phi_with_JJ_outlier", "track pT vs phi with JJ outlier;p_{T,track} (GeV/#it{c});#varphi_{track} (rad);entries", {HistType::kTH2F, {{200, 0, 200}, {160, -1.0, 7.0}}});
      registry.add("h_pt_hard_track_pt_with_JJ_outlier", "Tracks vs pThard;#frac{p_{T}}{#hat{p}};p_{T}", {HistType::kTH2F, {pThatAxis, {200, 0, 200}}});
      registry.add("h_track_pt_with_MB_outlier", "track pT with MB outlier", {HistType::kTH1F, {{200, 0, 200}}});
      registry.add("h_track_pt_eta_with_MB_outlier", "track pT vs eta with MB outlier;p_{T,track} (GeV/#it{c});#eta_{track};entries", {HistType::kTH2F, {{200, 0, 200}, {100, -5, 5}}});
      registry.add("h_track_pt_phi_with_MB_outlier", "track pT vs phi with MB outlier;p_{T,track} (GeV/#it{c});#varphi_{track} (rad);entries", {HistType::kTH2F, {{200, 0, 200}, {160, -1.0, 7.0}}});
      registry.add("h_pt_hard_track_pt_with_MB_outlier", "Tracks vs pThard;#frac{p_{T}}{#hat{p}};p_{T}", {HistType::kTH2F, {pThatAxis, {200, 0, 200}}});

      registry.add("h_track_pt_no_collision", "track pt no collision;p_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0, 200}}});
      registry.add("h_track_pt_collision", "track pt collision;p_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0, 200}}});
      registry.add("h2_track_pt_pt_hat_no_particle", "track pt vs pt hat no particle;p_{T,track} (GeV/#it{c});#hat{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{200, 0, 200}, {600, 0, 600}}});
      registry.add("h2_track_pt_pt_hat_particle", "track pt vs pt hat particle;p_{T,track} (GeV/#it{c});#hat{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{200, 0, 200}, {600, 0, 600}}});

      registry.add("h_track_pt_outlier", "weight track pt", {HistType::kTH1F, {{200, 0, 200}}});
      registry.add("h2_pt_hat_track_pt", "track; #hat{#it{p}_{T}} (GeV/#it{c});#it{p}_{T,track} (GeV/#it{c})", {HistType::kTH2F, {{600, 0, 600}, {200, 0, 200}}});
      registry.add("h2_pt_hat_track_pt_outlier", "track; #hat{#it{p}_{T}} (GeV/#it{c});#it{p}_{T,track} (GeV/#it{c})", {HistType::kTH2F, {{600, 0, 600}, {200, 0, 200}}});
      registry.add("h2_neighbour_pt_hat_outlier", "neighbour; distance from collision; #hat{#it{p}_{T}} (GeV/#it{c})", {HistType::kTH2F, {{15, -7.5, 7.5}, {600, 0, 600}}});
      registry.add("h2_neighbour_track_pt_outlier", "neighbour; distance from collision; #it{p}_{T,track} (GeV/#it{c})", {HistType::kTH2F, {{15, -7.5, 7.5}, {200, 0, 100}}});
      registry.add("h2_neighbour_pt_hat_all", "neighbour; distance from collision; #hat{#it{p}_{T}} (GeV/#it{c})", {HistType::kTH2F, {{15, -7.5, 7.5}, {600, 0, 600}}});
      registry.add("h2_neighbour_track_pt_all", "neighbour; distance from collision; #it{p}_{T,track} (GeV/#it{c})", {HistType::kTH2F, {{15, -7.5, 7.5}, {200, 0, 100}}});

      registry.add("h_track_pt_outlier_same_collision", "weight track pt same collision", {HistType::kTH1F, {{200, 0, 200}}});
      registry.add("h_track_pt_outlier_different_collision_JJ", "weight track pt different jet-jet collision", {HistType::kTH1F, {{200, 0, 200}}});
      registry.add("h_track_pt_outlier_different_collision_MB", "weight track pt different MB collision", {HistType::kTH1F, {{200, 0, 200}}});

      registry.add("h2_outlier_event_Ntracks_different_selected_JJ", "number of selected tracks from different jet-jet events", {HistType::kTH2F, {{600, 0, 600}, {200, 0, 200}}});
      registry.add("h2_outlier_event_Ntracks_different_selected_MB", "number of selected tracks from different MB events", {HistType::kTH2F, {{600, 0, 600}, {200, 0, 200}}});
      registry.add("h2_outlier_event_Ntracks_same_selected_JJ", "number of selected tracks from same jet-jet events", {HistType::kTH2F, {{600, 0, 600}, {200, 0, 200}}});
      registry.add("h2_outlier_event_tracks_frac_different_JJ", "fraction of tracks from different jet-jet events", {HistType::kTH2F, {{600, 0, 600}, {100, 0, 1}}});
      registry.add("h2_outlier_event_tracks_frac_different_MB", "fraction of tracks from different MB events", {HistType::kTH2F, {{600, 0, 600}, {100, 0, 1}}});
      registry.add("h2_outlier_event_tracks_frac_different_selected_JJ", "fraction of selected tracks from different jet-jet events", {HistType::kTH2F, {{600, 0, 600}, {100, 0, 1}}});
      registry.add("h2_outlier_event_tracks_frac_different_selected_MB", "fraction of selected tracks from different MB events", {HistType::kTH2F, {{600, 0, 600}, {100, 0, 1}}});
      registry.add("h2_outlier_collision_ID_difference", "difference in collision ID between outlier collision and analysed collision", {HistType::kTH2F, {{600, 0, 600}, {200, -100, 100}}});
      registry.add("h_DeltaZ_Outlier", "Delta Z between outlier collision and analysed collision", {HistType::kTH1F, {{1200, -30, 30}}});
      registry.add("h2_DeltaZ_Outlier_difference", "Delta Z between outlier collision and analysed collision vs difference in collision ID", {HistType::kTH2F, {{1200, -30, 30}, {200, -100, 100}}});

      registry.add("h_track_pt_same_collision", "track pt from same collision or different MB collision;p_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0, 200}}});
      registry.add("h_track_pt_eta_same_collision", "track pt vs eta from same collision or different MB collision;p_{T,track} (GeV/#it{c});#eta_{track};entries", {HistType::kTH2F, {{200, 0, 200}, {100, -5, 5}}});
      registry.add("h_track_pt_phi_same_collision", "track pt vs phi from same collision or different MB collision;p_{T,track} (GeV/#it{c});#varphi_{track} (rad);entries", {HistType::kTH2F, {{200, 0, 200}, {160, -1.0, 7.0}}});
      registry.add("h2_collision_ID_difference_same_collision", "difference in collision ID between outlier collision and analysed collision", {HistType::kTH2F, {{600, 0, 600}, {200, -100, 100}}});
      registry.add("h_track_pt_no_JJ_different", "track pt from same collision or different MB collision;p_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0, 200}}});
      registry.add("h_track_pt_eta_no_JJ_different", "track pt vs eta from same collision or different MB collision;p_{T,track} (GeV/#it{c});#eta_{track};entries", {HistType::kTH2F, {{200, 0, 200}, {100, -5, 5}}});
      registry.add("h_track_pt_phi_no_JJ_different", "track pt vs phi from same collision or different MB collision;p_{T,track} (GeV/#it{c});#varphi_{track} (rad);entries", {HistType::kTH2F, {{200, 0, 200}, {160, -1.0, 7.0}}});
      registry.add("h2_collision_ID_difference_no_JJ_different", "difference in collision ID between outlier collision and analysed collision", {HistType::kTH2F, {{600, 0, 600}, {200, -100, 100}}});
    }
  }

  Filter trackCuts = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax);
  Filter eventCuts = (nabs(aod::jcollision::posZ) < vertexZCut && aod::jcollision::centFT0M >= centralityMin && aod::jcollision::centFT0M < centralityMax);

  template <typename T, typename U>
  bool isAcceptedJet(U const& jet)
  {

    if (jetAreaFractionMin > -98.0) {
      if (jet.area() < jetAreaFractionMin * M_PI * (jet.r() / 100.0) * (jet.r() / 100.0)) {
        return false;
      }
    }
    bool checkConstituentPt = true;
    bool checkConstituentMinPt = (leadingConstituentPtMin > -98.0);
    bool checkConstituentMaxPt = (leadingConstituentPtMax < 9998.0);
    if (!checkConstituentMinPt && !checkConstituentMaxPt) {
      checkConstituentPt = false;
    }

    if (checkConstituentPt) {
      bool isMinLeadingConstituent = !checkConstituentMinPt;
      bool isMaxLeadingConstituent = true;

      for (const auto& constituent : jet.template tracks_as<T>()) {
        double pt = constituent.pt();

        if (checkConstituentMinPt && pt >= leadingConstituentPtMin) {
          isMinLeadingConstituent = true;
        }
        if (checkConstituentMaxPt && pt > leadingConstituentPtMax) {
          isMaxLeadingConstituent = false;
        }
      }
      return isMinLeadingConstituent && isMaxLeadingConstituent;
    }

    return true;
  }

  void fillHistogramsAmbiguous(soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetEventWeights>::iterator const& jet,
                               float weight,
                               aod::AmbiguousTracks const& tracksAmbiguous)
  {

    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    if (jet.pt() > pTHatMaxMCD * pTHat) {
      return;
    }
    // outlier ID checks
    if (jet.pt() > pTHatMaxMCDOutlier * pTHat) {
      registry.fill(HIST("h3_pthat_jet_pt_jet_ntracks_outliers"), pTHat, jet.pt(), jet.tracksIds().size(), weight);
      registry.fill(HIST("h3_pthat_jet_pt_jet_ntracks_outliers_noweight"), pTHat, jet.pt(), jet.tracksIds().size());
    } else {
      registry.fill(HIST("h3_pthat_jet_pt_jet_ntracks_all"), pTHat, jet.pt(), jet.tracksIds().size(), weight);
      registry.fill(HIST("h3_pthat_jet_pt_jet_ntracks_all_noweight"), pTHat, jet.pt(), jet.tracksIds().size());
    }
    for (auto& constituent : jet.template tracks_as<aod::JetTracksMCD>()) {
      if (jet.pt() > pTHatMaxMCDOutlier * pTHat) {
        registry.fill(HIST("h3_jet_pt_track_dEta_track_dPhi_outliers"), jet.pt(), jet.eta() - constituent.eta(), jet.phi() - constituent.phi(), weight);
        registry.fill(HIST("h3_jet_pt_track_dEta_track_dPhi_outliers_noweight"), jet.pt(), jet.eta() - constituent.eta(), jet.phi() - constituent.phi());
        registry.fill(HIST("h2_jet_pt_Dz_outliers_noweight"), jet.pt(), constituent.pt() / jet.pt());
      } else {
        registry.fill(HIST("h3_jet_pt_track_dEta_track_dPhi_all"), jet.pt(), jet.eta() - constituent.eta(), jet.phi() - constituent.phi(), weight);
        registry.fill(HIST("h3_jet_pt_track_dEta_track_dPhi_all_noweight"), jet.pt(), jet.eta() - constituent.eta(), jet.phi() - constituent.phi());
        registry.fill(HIST("h2_jet_pt_Dz_all_noweight"), jet.pt(), constituent.pt() / jet.pt());
      }
    }

    // ambiguous/unmatched track checks
    auto iterAmbiguous = tracksAmbiguous.begin();
    int nAmbTracks = 0;
    int nUnmatchedTracks = 0;
    double pt_total = 0;
    double pt_amb = 0;
    double pt_unmatched = 0;
    for (auto& constituent : jet.template tracks_as<aod::JetTracksMCD>()) {
      pt_total += constituent.pt();
      bool has_MCparticle = constituent.has_mcParticle();
      if (!has_MCparticle) {
        // LOG(info) << "constituent NO MC PARTICLE: track.index()=" << constituent.index() << " track.globalIndex()=" << constituent.globalIndex();
        registry.fill(HIST("h3_jet_pt_track_pt_pt_hat_no_particle"), jet.pt(), constituent.pt(), pTHat, weight);
        pt_unmatched += constituent.pt();
        nUnmatchedTracks++;
      } else {
        registry.fill(HIST("h3_jet_pt_track_pt_pt_hat_with_particle"), jet.pt(), constituent.pt(), pTHat, weight);
      }

      bool goFillHisto = (iterAmbiguous != tracksAmbiguous.end());
      if (goFillHisto) {
        while (constituent.globalIndex() > iterAmbiguous.trackId()) {
          iterAmbiguous++;
          if (iterAmbiguous == tracksAmbiguous.end()) { /// all ambiguous tracks found
            goFillHisto = false;
            break;
          }
        }
      }
      if (goFillHisto) {
        if (constituent.globalIndex() == iterAmbiguous.trackId()) {
          nAmbTracks++;
          registry.fill(HIST("h3_jet_pt_track_pt_pt_hat_ambiguous"), jet.pt(), constituent.pt(), pTHat, weight);
          pt_amb += constituent.pt();
          nAmbTracks++;
        } else {
          registry.fill(HIST("h3_jet_pt_track_pt_pt_hat_unambiguous"), jet.pt(), constituent.pt(), pTHat, weight);
        }
      }
    }
    registry.fill(HIST("h3_jet_pt_frac_pt_ambiguous_pt_hat"), jet.pt(), pt_amb / pt_total, pTHat, weight);
    registry.fill(HIST("h3_jet_pt_frac_constituents_ambiguous_pt_hat"), jet.pt(), double(nAmbTracks) / double(jet.template tracks_as<aod::JetTracksMCD>().size()), pTHat, weight);

    registry.fill(HIST("h3_jet_pt_frac_pt_unmatched_particle_pt_hat"), jet.pt(), pt_unmatched / pt_total, pTHat, weight);
    registry.fill(HIST("h3_jet_pt_frac_constituents_unmatched_particle_pt_hat"), jet.pt(), double(nUnmatchedTracks) / double(jet.template tracks_as<aod::JetTracksMCD>().size()), pTHat, weight);
  }

  void processJetsAmbiguous(soa::Filtered<soa::Join<aod::JetCollisions, aod::JMcCollisionLbs>>::iterator const& collision,
                            aod::JetMcCollisions const&,
                            soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetEventWeights> const& jets,
                            aod::JetTracksMCD const&,
                            const aod::AmbiguousTracks& tracksAmbiguous)
  {
    //
    // jet-based outlier checks based on ambiguous tracks
    //
    if (collision.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      return;
    }
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    bool isOutlierEvent = false;
    int nTracksJet = 0;
    float pTHat = collision.mcCollision().ptHard();
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(jet)) {
        continue;
      }
      fillHistogramsAmbiguous(jet, jet.eventWeight(), tracksAmbiguous);
      nTracksJet += jet.tracksIds().size();
      if (jet.pt() > pTHatMaxMCDOutlier * pTHat) {
        isOutlierEvent = true;
      }
    }
    if (isOutlierEvent) {
      registry.fill(HIST("h_collision_pTHat_multiplicity_outlier"), pTHat, nTracksJet, collision.weight());
      registry.fill(HIST("h_collision_multiplicity_outlier"), nTracksJet, collision.weight());
      registry.fill(HIST("h_collision_trackOccupancyInTimeRange_outlier"), collision.trackOccupancyInTimeRange(), collision.weight());
      registry.fill(HIST("h_collision_multFV0A_outlier"), collision.multFV0A(), collision.weight());
      registry.fill(HIST("h_collision_multFV0C_outlier"), collision.multFV0C(), collision.weight());
      registry.fill(HIST("h_collision_multFV0M_outlier"), collision.multFV0M(), collision.weight());
    } else {
      registry.fill(HIST("h_collision_multiplicity_all"), nTracksJet, collision.weight());
      registry.fill(HIST("h_collision_pTHat_multiplicity_all"), pTHat, nTracksJet, collision.weight());
      registry.fill(HIST("h_collision_trackOccupancyInTimeRange_all"), collision.trackOccupancyInTimeRange(), collision.weight());
      registry.fill(HIST("h_collision_multFV0A_all"), collision.multFV0A(), collision.weight());
      registry.fill(HIST("h_collision_multFV0C_all"), collision.multFV0C(), collision.weight());
      registry.fill(HIST("h_collision_multFV0M_all"), collision.multFV0M(), collision.weight());
    }
  }
  PROCESS_SWITCH(JetOutlierQATask, processJetsAmbiguous, "jet finder QA mcd with weighted events", false);

  void processCollisionsBC(soa::Join<aod::JetCollisions, aod::JMcCollisionLbs, aod::JCollisionBCs> const& collisions,
                           aod::JetMcCollisions const&)
  {
    //
    // collision-based outlier checks based on BC and z position
    // based on 2-event correlation checks in PWGDQ/Tasks/tableReader_withAssoc.cxx
    //

    fBCCollMap.clear();
    for (auto const& collision : collisions) {
      // Fill the BC map of events
      if (fBCCollMap.find(collision.bcId()) == fBCCollMap.end()) {
        std::vector<int64_t> evIndices = {collision.globalIndex()};
        fBCCollMap[collision.bcId()] = evIndices;
      } else {
        auto& evIndices = fBCCollMap[collision.bcId()];
        evIndices.push_back(collision.globalIndex());
      }
      registry.fill(HIST("h_Z"), collision.posZ());
    }

    // Create a map for collisions which are candidate of being split
    // key: event global index, value: whether pileup event is a possible splitting
    // (not used ATM, but could be used to flag events in the future)
    std::map<int64_t, bool> collisionSplittingMap;

    // loop over the BC map, get the collision vectors and make in-bunch and out of bunch 2-event correlations
    for (auto bc1It = fBCCollMap.begin(); bc1It != fBCCollMap.end(); ++bc1It) {
      uint64_t bc1 = bc1It->first;
      auto bc1Events = bc1It->second;

      // same bunch event correlations, if more than 1 collisions in this bunch
      if (bc1Events.size() > 1) {
        for (auto ev1It = bc1Events.begin(); ev1It != bc1Events.end(); ++ev1It) {
          auto ev1 = collisions.rawIteratorAt(*ev1It);
          for (auto ev2It = std::next(ev1It); ev2It != bc1Events.end(); ++ev2It) {
            auto ev2 = collisions.rawIteratorAt(*ev2It);
            // compute 2-event quantities and mark the candidate split collisions
            float deltaZ = ev1.posZ() - ev2.posZ();
            if (TMath::Abs(deltaZ) < splitCollisionsDeltaZ) { // this is a possible collision split
              collisionSplittingMap[*ev1It] = true;
              collisionSplittingMap[*ev2It] = true;
            }
            registry.fill(HIST("h_DeltaZ_InBunch"), deltaZ);
            registry.fill(HIST("h_DeltaZ_Z1_InBunch"), deltaZ, ev1.posZ());
            registry.fill(HIST("h_Z1_Z2_InBunch"), ev1.posZ(), ev2.posZ());
            if (ev1.subGeneratorId() != jetderiveddatautilities::JCollisionSubGeneratorId::mbGap &&
                ev2.subGeneratorId() != jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) { // both are non-gap events
              registry.fill(HIST("h_DeltaZ_InBunch_JJ"), deltaZ);
              registry.fill(HIST("h_DeltaZ_Z1_InBunch_JJ"), deltaZ, ev1.posZ());
              registry.fill(HIST("h_Z1_Z2_InBunch_JJ"), ev1.posZ(), ev2.posZ());
            }
          } // end second event loop
        } // end first event loop
      } // end if BC1 events > 1

      // loop over the following BCs in the TF
      for (auto bc2It = std::next(bc1It); bc2It != fBCCollMap.end(); ++bc2It) {
        uint64_t bc2 = bc2It->first;
        if ((bc2 > bc1 ? bc2 - bc1 : bc1 - bc2) > splitCollisionsDeltaBC) {
          break;
        }
        auto bc2Events = bc2It->second;

        // loop over events in the first BC
        for (auto ev1It : bc1Events) {
          auto ev1 = collisions.rawIteratorAt(ev1It);
          // loop over events in the second BC
          for (auto ev2It : bc2Events) {
            auto ev2 = collisions.rawIteratorAt(ev2It);
            // compute 2-event quantities and mark the candidate split collisions
            float deltaZ = ev1.posZ() - ev2.posZ();
            if (TMath::Abs(deltaZ) < splitCollisionsDeltaZ) { // this is a possible collision split
              collisionSplittingMap[ev1It] = true;
              collisionSplittingMap[ev2It] = true;
            }
            registry.fill(HIST("h_DeltaZ_OutOfBunch"), deltaZ);
            registry.fill(HIST("h_DeltaZ_Z1_OutOfBunch"), deltaZ, ev1.posZ());
            registry.fill(HIST("h_Z1_Z2_OutOfBunch"), ev1.posZ(), ev2.posZ());
            if (ev1.subGeneratorId() != jetderiveddatautilities::JCollisionSubGeneratorId::mbGap &&
                ev2.subGeneratorId() != jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) { // both are non-gap events
              registry.fill(HIST("h_DeltaZ_OutOfBunch_JJ"), deltaZ);
              registry.fill(HIST("h_DeltaZ_Z1_OutOfBunch_JJ"), deltaZ, ev1.posZ());
              registry.fill(HIST("h_Z1_Z2_OutOfBunch_JJ"), ev1.posZ(), ev2.posZ());
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(JetOutlierQATask, processCollisionsBC, "jet finder QA outliers", false);

  void processTracksBC(soa::Filtered<soa::Join<aod::JetCollisions, aod::JMcCollisionLbs, aod::JCollisionBCs>> const& collisions,
                       soa::Join<aod::JetMcCollisions, aod::JMcCollisionPIs> const&,
                       aod::JetMcCollisions const& collisionsMC,
                       aod::JetTracksMCD const& tracks,
                       JetParticlesWithOriginal const&)
  {
    //
    // track-based outlier checks
    //

    // first check for collisions occuring close by in time and z in MC
    std::set<int> closeByCollisionIDs;
    for (auto const& collisionMC : collisionsMC) {
      if (collisionMC.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
        continue;
      }
      float posZtrue = collisionMC.posZ();
      for (auto const& collisionCloseMC : collisionsMC) { // check for closeby collisions in MC
        int diffColl = collisionCloseMC.globalIndex() - collisionMC.globalIndex();
        if (diffColl >= mergeCollisionsDeltaMin && diffColl <= mergeCollisionsDeltaMax) { // check if n collisions prior or after
          if (collisionCloseMC.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
            continue;
          }
          if (diffColl == 0) {
            continue;
          }
          if (TMath::Abs(collisionCloseMC.posZ() - posZtrue) < splitCollisionsDeltaZPart) {
            closeByCollisionIDs.insert(collisionMC.globalIndex()); // Save the ID of the close-by collision
            break;                                                 // closeby collision in MC, don't use this event
          }
        }
      }
    }
    // now make reconstructed-level checks
    for (auto const& collision : collisions) { // loop over reconstructed collisions
      if (collision.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
        continue;
      }
      if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
        continue;
      }
      if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
        continue;
      }
      float weight = collision.weight();
      float pTHat = collision.mcCollision().ptHard();
      bool isOutlierEventDifferentJJCollision = false;
      bool isOutlierEventDifferentMBCollision = false;

      const auto tracksColl = tracks.sliceBy(perColTrack, collision.globalIndex());
      // fill track histograms for all collisions
      for (auto const& track : tracksColl) {
        if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
          continue;
        }
        registry.fill(HIST("h_track_pt"), track.pt(), weight);
        registry.fill(HIST("h_track_eta"), track.eta(), weight);
        registry.fill(HIST("h_track_phi"), track.phi(), weight);
        registry.fill(HIST("h_track_pt_eta"), track.pt(), track.eta(), weight);
        registry.fill(HIST("h_track_pt_phi"), track.pt(), track.phi(), weight);
        registry.fill(HIST("h_pt_hard_track_pt"), pTHat != 0.0 ? track.pt() / pTHat : 0.0, track.pt(), weight);
        // checks on track distributions with/without collision and with/without MC particle
        if (!track.has_collision() || track.collisionId() != collision.globalIndex()) {
          registry.fill(HIST("h_track_pt_no_collision"), track.pt(), weight);
        } else {
          registry.fill(HIST("h_track_pt_collision"), track.pt(), weight);
        }
        bool has_MCparticle = track.has_mcParticle();
        if (!has_MCparticle) {
          registry.fill(HIST("h2_track_pt_pt_hat_no_particle"), track.pt(), collision.mcCollision().ptHard(), weight);
        } else {
          registry.fill(HIST("h2_track_pt_pt_hat_particle"), track.pt(), collision.mcCollision().ptHard(), weight);
        }
        // check outlier tracks and neighbouring collisions
        registry.fill(HIST("h2_pt_hat_track_pt"), pTHat, track.pt());
        // get MC info about track and collision
        auto mcParticleOutlier = track.mcParticle_as<JetParticlesWithOriginal>();
        auto collisionMCOutlier = collisionsMC.sliceBy(perColParticle, mcParticleOutlier.mcCollisionId());
        if (collisionMCOutlier.size() != 1) {
          LOG(info) << "size of collision outlier not expected";
          return;
        }
        int mcCollisionIDcoll = collision.mcCollisionId(); // Get the corresponding MC collision ID from the reco collision
        int mcCollisionIDOutlier = mcParticleOutlier.mcCollisionId();
        int subGenIDOutlier = collisionMCOutlier.begin().subGeneratorId();
        int outlierCollisionIDDifference = mcCollisionIDOutlier - mcCollisionIDcoll;

        int nMBdifferent = 0;
        int nMBdifferentSelected = 0;
        int nJJdifferent = 0;
        int nJJdifferentSelected = 0;
        int nJJsame = 0;
        int nJJsameSelected = 0;
        // ID outlier based on track pT relative to pTHat
        if (track.pt() > pTHatMaxMCDOutlier * pTHat) { // high weight outlier track
          registry.fill(HIST("h_track_pt_outlier"), track.pt());
          registry.fill(HIST("h2_pt_hat_track_pt_outlier"), pTHat, track.pt());
          for (auto const& collisionOutlier : collisions) { // find collisions closeby
            int diffColl = collision.globalIndex() - collisionOutlier.globalIndex();
            if (abs(diffColl) < 6) {
              float eventWeightOutlier = collisionOutlier.mcCollision().weight();
              double pTHatOutlier = collisionOutlier.mcCollision().ptHard();
              registry.fill(HIST("h2_neighbour_pt_hat_outlier"), float(diffColl + 0.1), pTHatOutlier, eventWeightOutlier);
              registry.fill(HIST("h2_neighbour_track_pt_outlier"), float(diffColl + 0.1), track.pt(), eventWeightOutlier);
            }
          }
          // now match tracks to their MC particle, check the MC collision ID of this particle, and
          // check what fraction of tracks in this event are associated to this MC collision
          // LOG(info) << "--- Loop over tracks in outlier event with pT/pThat = " << track.pt() / pTHat << "---";
          // LOG(info) << "N tracks in outlier event = " << tracksColl.size() << " pTHat = " << pTHat << " collisionID = " << collision.globalIndex() << " mcCollisionID = " << collision.mcCollisionId();
          for (auto const& trackOutlier : tracksColl) {
            if (!trackOutlier.has_mcParticle()) {
              continue;
            }
            bool isTrackSelected = false;
            if (jetderiveddatautilities::selectTrack(trackOutlier, trackSelection)) {
              isTrackSelected = true;
            }
            auto mcParticle = trackOutlier.mcParticle_as<JetParticlesWithOriginal>();
            auto collisionMC = collisionsMC.sliceBy(perColParticle, mcParticle.mcCollisionId());
            if (collisionMC.size() == 0) {
              LOG(info) << "no collision found for mcCollisionID = " << mcParticle.mcCollisionId();
              continue;
            }
            int mcCollisionIDtrack = mcParticle.mcCollisionId(); // Get the corresponding MC collision ID
            int subGenID = collisionMC.begin().subGeneratorId();
            if (mcCollisionIDtrack == mcCollisionIDcoll) {
              nJJsame++;
              if (isTrackSelected) {
                registry.fill(HIST("h_track_pt_outlier_same_collision"), trackOutlier.pt());
                nJJsameSelected++;
              }
            } else {
              if (subGenID == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) { // MB-gap
                nMBdifferent++;
                if (isTrackSelected) {
                  registry.fill(HIST("h_track_pt_outlier_different_collision_MB"), trackOutlier.pt());
                  nMBdifferentSelected++;
                }
              } else { // jet-jet
                nJJdifferent++;
                if (isTrackSelected) {
                  registry.fill(HIST("h_track_pt_outlier_different_collision_JJ"), trackOutlier.pt());
                  nJJdifferentSelected++;
                }
              }
            }
          }
          // LOG(info) << "nJJsame = " << nJJsame << " nJJdifferent = " << nJJdifferent << " nMBdifferent = " << nMBdifferent;
          // LOG(info) << "nJJsameSelected = " << nJJsameSelected << " nJJdifferentSelected = " << nJJdifferentSelected << " nMBdifferentSelected = " << nMBdifferentSelected;
          registry.fill(HIST("h2_outlier_event_Ntracks_different_selected_JJ"), pTHat, nJJdifferentSelected);
          registry.fill(HIST("h2_outlier_event_Ntracks_different_selected_MB"), pTHat, nMBdifferentSelected);
          registry.fill(HIST("h2_outlier_event_Ntracks_same_selected_JJ"), pTHat, nJJsameSelected);
          registry.fill(HIST("h2_outlier_event_tracks_frac_different_selected_MB"), pTHat, float(nMBdifferentSelected) / float(nJJdifferentSelected + nJJsameSelected + nMBdifferentSelected));
          registry.fill(HIST("h2_outlier_event_tracks_frac_different_JJ"), pTHat, float(nJJdifferent) / float(nJJdifferent + nJJsame + nMBdifferent));
          registry.fill(HIST("h2_outlier_event_tracks_frac_different_MB"), pTHat, float(nMBdifferent) / float(nJJdifferent + nJJsame + nMBdifferent));
          registry.fill(HIST("h2_outlier_event_tracks_frac_different_selected_JJ"), pTHat, float(nJJdifferentSelected) / float(nJJdifferentSelected + nJJsameSelected + nMBdifferentSelected));
          registry.fill(HIST("h2_outlier_event_tracks_frac_different_selected_MB"), pTHat, float(nMBdifferentSelected) / float(nJJdifferentSelected + nJJsameSelected + nMBdifferentSelected));
          // now check where outlier comes from
          // LOG(info) <<"outlier comes from " << (mcCollisionIDOutlier == mcCollisionIDcoll ? "same" : "different") << " event which is a " << (subGenIDOutlier == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap ? " MB-gap" : " jet-jet") << " collision with mcCollisionID = " << mcCollisionIDOutlier;
          registry.fill(HIST("h2_outlier_collision_ID_difference"), pTHat, float(outlierCollisionIDDifference));
          // if outlier comes from different collision, check which type and set flags
          if (mcCollisionIDOutlier != mcCollisionIDcoll && subGenIDOutlier != jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
            isOutlierEventDifferentJJCollision = true;
            float deltaZ = collisionMCOutlier.begin().posZ() - collision.mcCollision().posZ();
            registry.fill(HIST("h_DeltaZ_Outlier"), deltaZ);
            registry.fill(HIST("h2_DeltaZ_Outlier_difference"), deltaZ, float(outlierCollisionIDDifference));
          } else if (mcCollisionIDOutlier != mcCollisionIDcoll && subGenIDOutlier == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
            isOutlierEventDifferentMBCollision = true;
          }
        }
        // fill for tracks from same collision or different MB collision in collisions that likely aren't fully merged
        if (nJJdifferentSelected < maxNTracksJJdifferent &&
            (subGenIDOutlier == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap ||
             mcCollisionIDOutlier == mcCollisionIDcoll)) {
          registry.fill(HIST("h_track_pt_same_collision"), track.pt(), weight);
          registry.fill(HIST("h_track_pt_eta_same_collision"), track.pt(), track.eta(), weight);
          registry.fill(HIST("h_track_pt_phi_same_collision"), track.pt(), track.phi(), weight);
          registry.fill(HIST("h2_collision_ID_difference_same_collision"), pTHat, float(outlierCollisionIDDifference));
        }
        // fill tracks for events which have no JJ outlier tracks from different events
        if (nJJdifferentSelected == 0) {
          registry.fill(HIST("h_track_pt_no_JJ_different"), track.pt(), weight);
          registry.fill(HIST("h_track_pt_eta_no_JJ_different"), track.pt(), track.eta(), weight);
          registry.fill(HIST("h_track_pt_phi_no_JJ_different"), track.pt(), track.phi(), weight);
          registry.fill(HIST("h2_collision_ID_difference_no_JJ_different"), pTHat, float(outlierCollisionIDDifference));
        }
        // collision checks for all tracks
        for (auto const& collisionOutlier : collisions) { // find collisions closeby
          float eventWeightOutlier = collisionOutlier.mcCollision().weight();
          double pTHatOutlier = collisionOutlier.mcCollision().ptHard();
          int diffColl = collision.globalIndex() - collisionOutlier.globalIndex();

          if (abs(diffColl) < 6) {
            // LOG(info) << "pThat = " << pTHat << "pThat neighbour = "<<pTHatOutlier;
            registry.fill(HIST("h2_neighbour_pt_hat_all"), float(diffColl + 0.1), pTHatOutlier, eventWeightOutlier);
            registry.fill(HIST("h2_neighbour_track_pt_all"), float(diffColl + 0.1), track.pt(), eventWeightOutlier);
          }
        }
      }

      // now check for close collisions in Z and fill track histograms excluding these events
      float posZtrue = collision.mcCollision().posZ();
      float posZrec = collision.posZ();
      registry.fill(HIST("h_Z_resolution"), posZtrue - posZrec);
      registry.fill(HIST("h_Z_resolution_wide"), posZtrue - posZrec);
      int mcCollisionID = collision.mcCollisionId(); // Get the corresponding MC collision ID
      if (closeByCollisionIDs.find(mcCollisionID) != closeByCollisionIDs.end()) {
        registry.fill(HIST("h_Z_reco_rejected"), posZrec);
        registry.fill(HIST("h_Z_true_rejected"), posZtrue);
        continue; // Skip this reconstructed collision as it corresponds to a close-by MC collision
      }
      registry.fill(HIST("h_Z_reco_accepted"), posZrec);
      registry.fill(HIST("h_Z_true_accepted"), posZtrue);

      // fill track histograms for accepted collisions
      for (auto const& track : tracksColl) {
        if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
          continue;
        }
        registry.fill(HIST("h_track_pt_accepted"), track.pt(), weight);
        registry.fill(HIST("h_track_eta_accepted"), track.eta(), weight);
        registry.fill(HIST("h_track_phi_accepted"), track.phi(), weight);
        registry.fill(HIST("h_track_pt_eta_accepted"), track.pt(), track.eta(), weight);
        registry.fill(HIST("h_track_pt_phi_accepted"), track.pt(), track.phi(), weight);
        registry.fill(HIST("h_pt_hard_track_pt_accepted"), pTHat != 0.0 ? track.pt() / pTHat : 0.0, track.pt(), weight);

        if (!isOutlierEventDifferentJJCollision) {
          registry.fill(HIST("h_track_pt_accepted_no_JJ_outlier"), track.pt(), weight);
          registry.fill(HIST("h_track_pt_eta_accepted_no_JJ_outlier"), track.pt(), track.eta(), weight);
          registry.fill(HIST("h_track_pt_phi_accepted_no_JJ_outlier"), track.pt(), track.phi(), weight);
          registry.fill(HIST("h_pt_hard_track_pt_accepted_no_JJ_outlier"), pTHat != 0.0 ? track.pt() / pTHat : 0.0, track.pt(), weight);
        } else {
          registry.fill(HIST("h_track_pt_with_JJ_outlier"), track.pt(), weight);
          registry.fill(HIST("h_track_pt_eta_with_JJ_outlier"), track.pt(), track.eta(), weight);
          registry.fill(HIST("h_track_pt_phi_with_JJ_outlier"), track.pt(), track.phi(), weight);
          registry.fill(HIST("h_pt_hard_track_pt_with_JJ_outlier"), pTHat != 0.0 ? track.pt() / pTHat : 0.0, track.pt(), weight);
        }
        if (isOutlierEventDifferentMBCollision) {
          registry.fill(HIST("h_track_pt_with_MB_outlier"), track.pt(), weight);
          registry.fill(HIST("h_track_pt_eta_with_MB_outlier"), track.pt(), track.eta(), weight);
          registry.fill(HIST("h_track_pt_phi_with_MB_outlier"), track.pt(), track.phi(), weight);
          registry.fill(HIST("h_pt_hard_track_pt_with_MB_outlier"), pTHat != 0.0 ? track.pt() / pTHat : 0.0, track.pt(), weight);
        }
      }
    }
  }
  PROCESS_SWITCH(JetOutlierQATask, processTracksBC, "jet finder QA outliers", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<JetOutlierQATask>(cfgc, TaskName{"jet-outlier-qa"})}; // o2-linter: disable=name/o2-task,name/workflow-file
}
