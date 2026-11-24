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
/// \file neutralMesonsInFullJets.cxx
/// \brief Task for neutral meson measurement in full jets
/// \author Jan Scharf <jan.scharf@cern.ch>

/// TO DO : Add Neutral Mesons measurement to the current full jet measuremtn

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/runDataProcessing.h>

#include <TError.h>

#include <fairlogger/Logger.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <unordered_set>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct neutralMesonInFullJets {
  HistogramRegistry registry;
  // define different enum classes to avoid "magic numbers"
  // DataType used in eventSelection and saving of histograms, and more
  // MatchingScheme defines the typ of matching between MCP and MCD jets
  // ResponseType for saving response histograms
  // JetAcceptanceCut for customization of the eta and phi jet cut
  enum class DataType {
    Data,
    MCP,
    MCD
  };
  enum class MatchingScheme {
    Geo,
    Pt,
    GeoAndPt
  };
  enum class ResponseType {
    Response,
    Fake,
    Miss
  };
  enum class JetAccpetanceCut {
    FullR,
    HalfR
  };

  // CONFIGURABLES
  // General event Configurables
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> pTHatMax{"pTHatMax", 2.5, "maximum fraction of hard scattering for jet acceptance in MC"};
  Configurable<float> pTHatExponent{"pTHatExponent", 4.0, "exponent of the event weight for the calculation of pTHeventSelectionBitsat"}; // 6 for MB MC and 4 for JJ MC
  Configurable<float> pTHatMin{"pTHatAbsoluteMin", 0.0, "minimum value of pTHat"};
  Configurable<bool> doMBGapTrigger{"doMBGapTrigger", true, "set to true only when using MB-Gap Trigger JJ MC to reject MB events at the collision and track level"};
  Configurable<bool> applyRCTSelections{"applyRCTSelections", true, "decide to apply RCT selections"};
  Configurable<std::string> eventSelections{"eventSelections", "selMC", "choose event selection"}; // selMC for MC and sel8 for data
  Configurable<int> matchingScheme{"matchingScheme", 2, "Matching scheme: 0 = Geo, 1 = Pt, 2 = GeoAndPt"};

  // Jet Configurables
  Configurable<float> jetPtMin{"jetPtMin", 5.0, "minimum jet pT cut"};
  Configurable<int> jetAccpetanceCut{"jetAcceptanceCut", 0, "Jet acceptance cut: 0 = FullR, 1 = HalfR"};
  Configurable<float> jetEtaMax{"jetEtaMax", 0.7, "Min jet eta based on EMCAL acceptance"};
  Configurable<float> jetPhiMin{"jetPhiMin", 1.3962634, "Min jet phi based on EMCAL acceptance"};
  Configurable<float> jetPhiMax{"jetPhiMax", 3.2836100, "Max jet phi based on EMCAL acceptance"};

  // Track Configurables

  // Cluster or EMCAL Configurables

  std::vector<int> eventSelectionBits;

  void init(InitContext const&)
  {
    // Create all histograms to safe later
    // processDataJets:
    registry.add("h_data_jet_pt", "data jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 300.}}});
    registry.add("h_data_jet_eta", "data jet #eta; #eta; entries", {HistType::kTH1F, {{100, -0.8, 0.8}}});
    registry.add("h_data_jet_phi", "data jet #phi; #phi; entries", {HistType::kTH1F, {{100, 0, 2 * o2::constants::math::PI}}});
    registry.add("h_data_jet_eta_phi", "data jet #eta #phi; #eta; #phi; entries", {HistType::kTH2F, {{100, -0.8, 0.8}, {100, 0, 2 * o2::constants::math::PI}}});

    // processMCPJets
    registry.add("h_mcp_jet_pt", "MCP jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 300.}}});
    registry.add("h_mcp_jet_eta", "MCP jet #eta; #eta; entries", {HistType::kTH1F, {{100, -0.8, 0.8}}});
    registry.add("h_mcp_jet_phi", "MCP jet #phi; #phi; entries", {HistType::kTH1F, {{100, 0, 2 * o2::constants::math::PI}}});
    registry.add("h_mcp_jet_eta_phi", "MCP jet #eta #phi; #eta; #phi; entries", {HistType::kTH2F, {{100, -0.8, 0.8}, {100, 0, 2 * o2::constants::math::PI}}});

    // processMCDJets
    registry.add("h_mcd_jet_pt", "MCD jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0., 300.}}});
    registry.add("h_mcd_jet_eta", "MCD jet #eta; #eta; entries", {HistType::kTH1F, {{100, -0.8, 0.8}}});
    registry.add("h_mcd_jet_phi", "MCD jet #phi; #phi; entries", {HistType::kTH1F, {{100, 0, 2 * o2::constants::math::PI}}});
    registry.add("h_mcd_jet_eta_phi", "MCD jet #eta #phi; #eta; #phi; entries", {HistType::kTH2F, {{100, -0.8, 0.8}, {100, 0, 2 * o2::constants::math::PI}}});

    // Matching histograms using Geo "0" matching
    registry.add("h_jet_response_pt_geomatching", "Matched jets using geometrical matching; #p_{T, MCD}; #p_{T, MCP}; entries", {HistType::kTH2F, {{300, 0., 300.}, {300, 0., 300.}}});
    registry.add("h_jet_reponse_pt_fake_geomatching", "Fake jets using geometrical matching; #p_{T, MCP}; entries", {HistType::kTH1F, {{300, 0., 300.}}});
    registry.add("h_jet_reponse_pt_miss_geomatching", "Miss jets using geometrical matching; #p_{T, MCD}; entries", {HistType::kTH1F, {{300, 0., 300.}}});
    registry.add("h_jet_response_phi_geomatching", "Matched jets using geometrical matching; #phi_{MCD}; #phi_{MCP}; entries", {HistType::kTH2F, {{100, 0, 2 * o2::constants::math::PI}, {100, 0, 2 * o2::constants::math::PI}}});
    registry.add("h_jet_reponse_phi_fake_geomatching", "Fake jets using geometrical matching; #phi_{MCP}; entries", {HistType::kTH1F, {{100, 0, 2 * o2::constants::math::PI}}});
    registry.add("h_jet_reponse_phi_miss_geomatching", "Miss jets using geometrical matching; #phi_{MCD}; entries", {HistType::kTH1F, {{100, 0, 2 * o2::constants::math::PI}}});
    registry.add("h_jet_response_eta_geomatching", "Matched jets using geometrical matching; #eta_{MCD}; #eta_{MCP}; entries", {HistType::kTH2F, {{100, -0.8, 0.8}, {100, -0.8, 0.8}}});
    registry.add("h_jet_reponse_eta_fake_geomatching", "Fake jets using geometrical matching; #eta_{MCP}; entries", {HistType::kTH1F, {{100, -0.8, 0.8}}});
    registry.add("h_jet_reponse_eta_miss_geomatching", "Miss jets using geometrical matching; #eta_{MCD}; entries", {HistType::kTH1F, {{100, -0.8, 0.8}}});

    // Matching histograms using Pt "1" matching
    registry.add("h_jet_response_pt_ptmatching", "Matched jets using pT-fraction matching; #p_{T, MCD}; #p_{T, MCP}; entries", {HistType::kTH2F, {{300, 0., 300.}, {300, 0., 300.}}});
    registry.add("h_jet_reponse_pt_fake_ptmatching", "Fake jets using pT-fraction matching; #p_{T, MCP}; entries", {HistType::kTH1F, {{300, 0., 300.}}});
    registry.add("h_jet_reponse_pt_miss_ptmatching", "Miss jets using pT-fraction matching; #p_{T, MCD}; entries", {HistType::kTH1F, {{300, 0., 300.}}});
    registry.add("h_jet_response_phi_ptmatching", "Matched jets using pT-fraction matching; #phi_{MCD}; #phi_{MCP}; entries", {HistType::kTH2F, {{100, 0, 2 * o2::constants::math::PI}, {100, 0, 2 * o2::constants::math::PI}}});
    registry.add("h_jet_reponse_phi_fake_ptmatching", "Fake jets using pT-fraction matching; #phi_{MCP}; entries", {HistType::kTH1F, {{100, 0, 2 * o2::constants::math::PI}}});
    registry.add("h_jet_reponse_phi_miss_ptmatching", "Miss jets using pT-fraction matching; #phi_{MCD}; entries", {HistType::kTH1F, {{100, 0, 2 * o2::constants::math::PI}}});
    registry.add("h_jet_response_eta_ptmatching", "Matched jets using pT-fraction matching; #eta_{MCD}; #eta_{MCP}; entries", {HistType::kTH2F, {{100, -0.8, 0.8}, {100, -0.8, 0.8}}});
    registry.add("h_jet_reponse_eta_fake_ptmatching", "Fake jets using pT-fraction matching; #eta_{MCP}; entries", {HistType::kTH1F, {{100, -0.8, 0.8}}});
    registry.add("h_jet_reponse_eta_miss_ptmatching", "Miss jets using pT-fraction matching; #eta_{MCD}; entries", {HistType::kTH1F, {{100, -0.8, 0.8}}});

    // Matching histograms using Geo and Pt "2" matching
    registry.add("h_jet_response_pt_geoptmatching", "Matched jets using geometrical AND pT-fraction matching; #p_{T, MCD}; #p_{T, MCP}; entries", {HistType::kTH2F, {{300, 0., 300.}, {300, 0., 300.}}});
    registry.add("h_jet_reponse_pt_fake_geoptmatching", "Fake jets using geometrical AND pT-fraction matching; #p_{T, MCP}; entries", {HistType::kTH1F, {{300, 0., 300.}}});
    registry.add("h_jet_reponse_pt_miss_geoptmatching", "Miss jets using geometrical AND pT-fraction matching; #p_{T, MCD}; entries", {HistType::kTH1F, {{300, 0., 300.}}});
    registry.add("h_jet_response_phi_geoptmatching", "Matched jets using geometrical AND pT-fraction matching; #phi_{MCD}; #phi_{MCP}; entries", {HistType::kTH2F, {{100, 0, 2 * o2::constants::math::PI}, {100, 0, 2 * o2::constants::math::PI}}});
    registry.add("h_jet_reponse_phi_fake_geoptmatching", "Fake jets using geometrical AND pT-fraction matching; #phi_{MCP}; entries", {HistType::kTH1F, {{100, 0, 2 * o2::constants::math::PI}}});
    registry.add("h_jet_reponse_phi_miss_geoptmatching", "Miss jets using geometrical AND pT-fraction matching; #phi_{MCD}; entries", {HistType::kTH1F, {{100, 0, 2 * o2::constants::math::PI}}});
    registry.add("h_jet_response_eta_geoptmatching", "Matched jets using geometrical AND pT-fraction matching; #eta_{MCD}; #eta_{MCP}; entries", {HistType::kTH2F, {{100, -0.8, 0.8}, {100, -0.8, 0.8}}});
    registry.add("h_jet_reponse_eta_fake_geoptmatching", "Fake jets using geometrical AND pT-fraction matching; #eta_{MCP}; entries", {HistType::kTH1F, {{100, -0.8, 0.8}}});
    registry.add("h_jet_reponse_eta_miss_geoptmatching", "Miss jets using geometrical AND pT-fraction matching; #eta_{MCD}; entries", {HistType::kTH1F, {{100, -0.8, 0.8}}});

    // QA
    registry.add("h_mcd_collision_counter", "MCD events ; event status; entries", {HistType::kTH1F, {{12, 0.0, 12.0}}});
    registry.add("h_mcp_collision_counter", "MCP events ; event status; entries", {HistType::kTH1F, {{12, 0.0, 12.0}}});
    registry.add("h_data_jet_acceptance_counter", "Data jet acceptance Counter; Data jets; jet status; entries", {HistType::kTH1F, {{2, 0.0, 2.0}}});
    registry.add("h_mcd_jet_acceptance_counter", "MCD jet acceptance Counter; MCD jets; jet status; entries", {HistType::kTH1F, {{2, 0.0, 2.0}}});
    registry.add("h_mcp_jet_acceptance_counter", "MCP jet acceptance Counter; MCP jets; jet status; entries", {HistType::kTH1F, {{2, 0.0, 2.0}}});
  }

  // TABLES
  // Matching tables in both directions
  using JetTableMCDMatchedJoined = soa::Join<aod::FullMCDetectorLevelJets, aod::FullMCDetectorLevelJetsMatchedToFullMCParticleLevelJets>;
  using JetTableMCPMatchedJoined = soa::Join<aod::FullMCParticleLevelJets, aod::FullMCParticleLevelJetsMatchedToFullMCDetectorLevelJets>;
  // Tables with added EMCAL information, especially for acceptance
  using EMCCollisionsData = o2::soa::Join<aod::JetCollisions, aod::JEMCCollisionLbs>;   // JetCollisions with EMCAL Collision Labels
  using EMCCollisionsMCD = o2::soa::Join<aod::JetCollisionsMCD, aod::JEMCCollisionLbs>; // where, JetCollisionsMCD = JetCollisions+JMcCollisionLbs

  // Filters
  Filter jetCuts = aod::jet::pt > jetPtMin;
  Filter collisionFilter = nabs(aod::jcollision::posZ) < vertexZCut; // Filter is not sufficient and is performed again manually in the event Selection

  PresliceOptional<JetTableMCPMatchedJoined> JetMCPPerMcCollision = aod::jet::mcCollisionId;                      // used in MC matching
  PresliceUnsortedOptional<aod::JetCollisionsMCD> CollisionsPerMCPCollision = aod::jmccollisionlb::mcCollisionId; // ensure a for a MCP collision a matched MCD collision is accepted

  // HELPER FUNCTIONS
  // fill jet pt, eta and phi histograms for DataType Data, MCP and MCD
  template <DataType T, typename Jet>
  void fillJetHistograms(Jet const& jet,
                         float weight = 1.f)
  {
    if constexpr (T == DataType::Data) {
      registry.fill(HIST("h_data_jet_pt"), jet.pt(), weight);
      registry.fill(HIST("h_data_jet_eta"), jet.eta(), weight);
      registry.fill(HIST("h_data_jet_phi"), jet.phi(), weight);
      registry.fill(HIST("h_data_jet_eta_phi"), jet.eta(), jet.phi(), weight);
    } else if constexpr (T == DataType::MCP) {
      registry.fill(HIST("h_mcp_jet_pt"), jet.pt(), weight);
      registry.fill(HIST("h_mcp_jet_eta"), jet.eta(), weight);
      registry.fill(HIST("h_mcp_jet_phi"), jet.phi(), weight);
      registry.fill(HIST("h_mcp_jet_eta_phi"), jet.eta(), jet.phi(), weight);
    } else if constexpr (T == DataType::MCD) {
      registry.fill(HIST("h_mcd_jet_pt"), jet.pt(), weight);
      registry.fill(HIST("h_mcd_jet_eta"), jet.eta(), weight);
      registry.fill(HIST("h_mcd_jet_phi"), jet.phi(), weight);
      registry.fill(HIST("h_mcd_jet_eta_phi"), jet.eta(), jet.phi(), weight);
    }
  }

  // Helper function to fill collision counter hisograms for QA
  template <DataType T>
  void fillCollisionCounter(float step,
                            float weight)
  {
    if constexpr (T == DataType::MCD) {
      registry.fill(HIST("h_mcd_collision_counter"), step, weight);
    } else if constexpr (T == DataType::MCP) {
      registry.fill(HIST("h_mcp_collision_counter"), step, weight);
    }
  }

  // Jet-level acceptance cut on eta and phi, based on the acceptance of the EMCAL
  template <typename T>
  bool isJetInAcceptance(T const& jet) const
  {
    const double R = static_cast<double>(jet.r()) / 100.0;
    double radiusFraction = 1.0;

    // Different cut types, could be customized furthur and more methods can be added
    // radius Fraction sets the range in which jet is accepted, for FullR eta leads to -0.3 to 0.3
    switch (static_cast<JetAccpetanceCut>(jetAccpetanceCut.value)) {
      case JetAccpetanceCut::FullR:
        radiusFraction = 1.0;
        break;
      case JetAccpetanceCut::HalfR:
        radiusFraction = 0.5;
        break;
      default:
        LOGF(fatal, "Invalud jet acceptance cut selection: %d", jetAccpetanceCut.value);
    }
    const double margin = radiusFraction * R;

    // Eta acceptance
    const double jetFidEtaMax = jetEtaMax.value - margin;
    if (std::abs(jet.eta()) > jetFidEtaMax) {
      return false;
    }
    // Phi acceptance
    double jetFidPhiMin = jetPhiMin.value + margin;
    double jetFidPhiMax = jetPhiMax.value - margin;

    // Normalize phi to [0, 2pi)
    auto norm = [&](double a) {
      while (a < 0) {
        a += o2::constants::math::TwoPI;
      }
      while (a >= o2::constants::math::TwoPI) {
        a -= o2::constants::math::TwoPI;
      }
      return a;
    };

    // Normalize using the above function
    const double phi = norm(jet.phi());
    jetFidPhiMin = norm(jetFidPhiMin);
    jetFidPhiMax = norm(jetFidPhiMax);

    // cut on Phi
    if (jetFidPhiMin <= jetFidPhiMax) {
      return phi >= jetFidPhiMin && phi <= jetFidPhiMax;
    }

    // Wrap around case, the method should be fixed either by cutting the residual fraction or the jet radius
    LOGF(info, "Feducial window is wraped around (PhiMax < PhiMin).");
    return phi >= jetFidPhiMin || phi <= jetFidPhiMax;
  }

  // Helper function to check if for a MCP collision a matched MCD collision exisits that is accepted by the event Selection
  template <typename ParticleCollision, typename DetectorCollision>
  bool hasAcceptedDetectorCollision(ParticleCollision const& partcollisions,
                                    DetectorCollision const& detcollisions,
                                    float weight)
  {
    // To select MCP collision, it is ensured that it can be matched to accepted MCD collision, otherwise this collision is not usefull for matching and further study
    // In other words, the MCP collision needs to be measureable

    // Find all matched MCD collision to the specific MCP collision (no direct matching between MCP->MCD, therefore more possible)
    auto collisionspermcpjet = detcollisions.sliceBy(CollisionsPerMCPCollision, partcollisions.globalIndex());
    if (collisionspermcpjet.size() < 1) // at least one match is required
    {
      fillCollisionCounter<DataType::MCP>(7.5f, weight); // collisions that are cut, because they have no matched MCD collisions
      return false;
    }
    for (auto const& collision : collisionspermcpjet) {
      // detecor-level collision should also be selected based on the event Selection (not all cuts like vertex Z cut require, because that is already ensured by cutting on MCP collision)
      if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, doMBGapTrigger, applyRCTSelections, "CBT_calo")) {
        continue;
      }
      if (!collision.isAmbiguous() && jetderiveddatautilities::eventEMCAL(collision) && collision.alias_bit(kTVXinEMC)) {
        fillCollisionCounter<DataType::MCP>(8.5f, weight); // collisions that are cut, because the matched MCD collisions are not accepted by the standard event selection or are unmeasurebale in EMCAL
        return true;
      }
    }
    return false;
  }

  // Event-level selection for Data, MCP and MCD
  template <DataType T, typename Collision, typename Jets, typename DetectorCollisions = std::nullptr_t>
  bool eventSelection(Collision const& collision,
                      Jets const& jets,
                      float weight = 1.f,
                      float pTHat = 0.f,
                      DetectorCollisions const* detectorCollisions = nullptr)
  {
    // 0: All collisions
    fillCollisionCounter<T>(0.5, weight); // all collisions

    // 1: Additional Vertx Z selection (Filter is not sufficient)
    if (std::abs(collision.posZ()) > vertexZCut) {
      return false;
    }
    fillCollisionCounter<T>(1.5, weight); // collisions after vertex z cut

    // 2: MC pT-Hat outlier rejection
    if constexpr (T != DataType::Data) {
      if (pTHat <= 0.f || pTHat < pTHatMin) {
        fillCollisionCounter<T>(2.5f, weight); // collisions cut by the pT-Hat outlier rejection
        return false;
      }
      for (auto const& jet : jets) {
        if (jet.pt() > pTHatMax * pTHat) {
          fillCollisionCounter<T>(2.5f, weight); // collisions cut by the pT-Hat outlier rejection
          return false;
        }
      }
    }

    // not available for particle-level MC, workaround below
    if constexpr (T != DataType::MCP) {
      // 3: Standard event selection
      if (!jetderiveddatautilities::selectCollision(
            collision,
            eventSelectionBits,
            doMBGapTrigger,
            applyRCTSelections, "CBT_calo")) {
        fillCollisionCounter<T>(3.5f, weight); // collisions cut by the standard event selection (also MB Gap)
        return false;
      }

      // 4: EMCAL event selection
      if (collision.isAmbiguous() || !jetderiveddatautilities::eventEMCAL(collision) || !collision.alias_bit(kTVXinEMC)) {
        fillCollisionCounter<T>(4.5f, weight); // collisions cut by not measured by EMCAL
        return false;
      }
    }
    // workaournd for MCP
    else if constexpr (T == DataType::MCP) {
      // 3 and 4: MB Gap Selection + EMCAL selection using associated MCD collisions
      if (doMBGapTrigger && collision.getSubGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
        fillCollisionCounter<DataType::MCP>(6.5f, weight); // collisions cut by MB Gap
        return false;
      }
      if (detectorCollisions == nullptr || !hasAcceptedDetectorCollision(collision, *detectorCollisions, weight)) {
        return false;
      }
    }

    return true;
  }

  // retrieve matching mcpjet for a mcdjet based on the selected MatchingScheme, Geo, Pt or Geo and Pt
  template <MatchingScheme Scheme, typename MCDJet, typename Callback>
  void forEachAcceptedMatch(MCDJet const& mcdjet, Callback&& callback)
  {
    // Geo matching on the delta R closeness of the jets, Pt matching on the pt fraction of the matched constituents in the jet itself
    // Pt matching can be directional and needs to be accepted in both directions MCD->MCP and MCP->MCD (Geo is not directional MCD<->MCP)

    if constexpr (Scheme == MatchingScheme::Geo) {
      if (mcdjet.has_matchedJetGeo()) {
        for (auto const& mcpjet : mcdjet.template matchedJetGeo_as<JetTableMCPMatchedJoined>()) {
          callback(mcpjet);
        }
      }
    }

    else if constexpr (Scheme == MatchingScheme::Pt) {
      if (mcdjet.has_matchedJetPt()) {
        for (auto const& mcpjet : mcdjet.template matchedJetPt_as<JetTableMCPMatchedJoined>()) {
          // first loop for MCD->MCP matching
          bool hasReversePtMatch = false;
          if (mcpjet.has_matchedJetPt()) {
            for (auto const& matchedMCD : mcpjet.template matchedJetPt_as<JetTableMCDMatchedJoined>()) {
              // ensuring MCP->MCD matching
              if (matchedMCD.globalIndex() == mcdjet.globalIndex()) {
                hasReversePtMatch = true;
                break;
              }
            }
          }

          if (hasReversePtMatch) {
            callback(mcpjet);
          }
        }
      }
    }

    else if constexpr (Scheme == MatchingScheme::GeoAndPt) {
      std::unordered_set<int64_t> ptMatchedMCPJetIds;
      if (mcdjet.has_matchedJetPt()) {
        for (auto const& mcpjet : mcdjet.template matchedJetPt_as<JetTableMCPMatchedJoined>()) {
          ptMatchedMCPJetIds.insert(mcpjet.globalIndex());
        }
      }
      if (mcdjet.has_matchedJetGeo()) {
        for (auto const& mcpjet : mcdjet.template matchedJetGeo_as<JetTableMCPMatchedJoined>()) {
          if (ptMatchedMCPJetIds.find(mcpjet.globalIndex()) == ptMatchedMCPJetIds.end()) {
            continue;
          }
          bool hasReversePtMatch = false;
          if (mcpjet.has_matchedJetPt()) {
            for (auto const& matchedMCD : mcpjet.template matchedJetPt_as<JetTableMCDMatchedJoined>()) {
              if (matchedMCD.globalIndex() == mcdjet.globalIndex()) {
                hasReversePtMatch = true;
                break;
              }
            }
          }
          if (hasReversePtMatch) {
            callback(mcpjet);
          }
        }
      }
    }
  }

  // Helper function fill the histograms of the matching part 1 only response
  template <MatchingScheme Scheme, ResponseType Type, typename MCDJet, typename MCPJet>
  void fillMatchingResponse(MCDJet const& mcdjet, MCPJet const& mcpjet, float weight)
  {
    static_assert(Type == ResponseType::Response, "This is only for response histograms");
    if constexpr (Scheme == MatchingScheme::Geo) {
      registry.fill(HIST("h_jet_response_pt_geomatching"), mcdjet.pt(), mcpjet.pt(), weight);
      registry.fill(HIST("h_jet_response_phi_geomatching"), mcdjet.phi(), mcpjet.phi(), weight);
      registry.fill(HIST("h_jet_response_eta_geomatching"), mcdjet.eta(), mcpjet.eta(), weight);
    } else if constexpr (Scheme == MatchingScheme::Pt) {
      registry.fill(HIST("h_jet_response_pt_ptmatching"), mcdjet.pt(), mcpjet.pt(), weight);
      registry.fill(HIST("h_jet_response_phi_ptmatching"), mcdjet.phi(), mcpjet.phi(), weight);
      registry.fill(HIST("h_jet_response_eta_ptmatching"), mcdjet.eta(), mcpjet.eta(), weight);
    } else if constexpr (Scheme == MatchingScheme::GeoAndPt) {
      registry.fill(HIST("h_jet_response_pt_geoptmatching"), mcdjet.pt(), mcpjet.pt(), weight);
      registry.fill(HIST("h_jet_response_phi_geoptmatching"), mcdjet.phi(), mcpjet.phi(), weight);
      registry.fill(HIST("h_jet_response_eta_geoptmatching"), mcdjet.eta(), mcpjet.eta(), weight);
    }
  }

  // Helper function fill the histograms of the matching part 2 for miss and fake
  template <MatchingScheme Scheme, ResponseType Type, typename Jet>
  void fillMatchingResponse(Jet const& jet, float weight)
  {
    static_assert(Type != ResponseType::Response,
                  "Response histograms require both jets");

    if constexpr (Type == ResponseType::Fake) {

      if constexpr (Scheme == MatchingScheme::Geo) {
        registry.fill(HIST("h_jet_reponse_pt_fake_geomatching"), jet.pt(), weight);
        registry.fill(HIST("h_jet_reponse_phi_fake_geomatching"), jet.phi(), weight);
        registry.fill(HIST("h_jet_reponse_eta_fake_geomatching"), jet.eta(), weight);

      } else if constexpr (Scheme == MatchingScheme::Pt) {
        registry.fill(HIST("h_jet_reponse_pt_fake_ptmatching"), jet.pt(), weight);
        registry.fill(HIST("h_jet_reponse_phi_fake_ptmatching"), jet.phi(), weight);
        registry.fill(HIST("h_jet_reponse_eta_fake_ptmatching"), jet.eta(), weight);

      } else if constexpr (Scheme == MatchingScheme::GeoAndPt) {
        registry.fill(HIST("h_jet_reponse_pt_fake_geoptmatching"), jet.pt(), weight);
        registry.fill(HIST("h_jet_reponse_phi_fake_geoptmatching"), jet.phi(), weight);
        registry.fill(HIST("h_jet_reponse_eta_fake_geoptmatching"), jet.eta(), weight);
      }

    } else if constexpr (Type == ResponseType::Miss) {

      if constexpr (Scheme == MatchingScheme::Geo) {
        registry.fill(HIST("h_jet_reponse_pt_miss_geomatching"), jet.pt(), weight);
        registry.fill(HIST("h_jet_reponse_phi_miss_geomatching"), jet.phi(), weight);
        registry.fill(HIST("h_jet_reponse_eta_miss_geomatching"), jet.eta(), weight);

      } else if constexpr (Scheme == MatchingScheme::Pt) {
        registry.fill(HIST("h_jet_reponse_pt_miss_ptmatching"), jet.pt(), weight);
        registry.fill(HIST("h_jet_reponse_phi_miss_ptmatching"), jet.phi(), weight);
        registry.fill(HIST("h_jet_reponse_eta_miss_ptmatching"), jet.eta(), weight);

      } else if constexpr (Scheme == MatchingScheme::GeoAndPt) {
        registry.fill(HIST("h_jet_reponse_pt_miss_geoptmatching"), jet.pt(), weight);
        registry.fill(HIST("h_jet_reponse_phi_miss_geoptmatching"), jet.phi(), weight);
        registry.fill(HIST("h_jet_reponse_eta_miss_geoptmatching"), jet.eta(), weight);
      }
    }
  }

  // run the matching and the mcdjet loop, also aplly the event and jet selection
  template <MatchingScheme Scheme>
  void runJetMatching(
    EMCCollisionsMCD::iterator const& collision,
    JetTableMCDMatchedJoined const& mcdjets,
    JetTableMCPMatchedJoined const& mcpjets)
  {
    const float weight = collision.mcCollision().weight();
    const float pTHat = 10. / std::pow(weight, 1.0 / pTHatExponent);

    // MCD event selection
    if (!eventSelection<DataType::MCD>(collision, mcdjets, weight, pTHat)) {
      return;
    }

    // Retrieve MCP jets associated with this MC collision
    auto mcpJetsPerMcCollision = mcpjets.sliceBy(JetMCPPerMcCollision, collision.mcCollisionId());
    std::unordered_set<int64_t> matchedMCPJetIds;

    for (const auto& mcdjet : mcdjets) {
      bool hasAcceptedMatch = false;

      // Detector-level acceptance
      if (!isJetInAcceptance(mcdjet)) {
        continue;
      }

      forEachAcceptedMatch<Scheme>(mcdjet, [&](auto const& mcpjet) { // Jet Matching
        // Particle-level acceptance
        if (!isJetInAcceptance(mcpjet)) {
          return;
        }
        hasAcceptedMatch = true;
        matchedMCPJetIds.insert(mcpjet.globalIndex());
        fillMatchingResponse<Scheme, ResponseType::Response>(mcdjet, mcpjet, weight);
      });

      // Fake (Unmatched but accepted MCD jets)
      if (!hasAcceptedMatch) {
        fillMatchingResponse<Scheme, ResponseType::Fake>(mcdjet, weight);
      }
    }
    // Miss (Unmatched but accepted MCP jets)
    for (const auto& mcpjet : mcpJetsPerMcCollision) {

      // Particle-level acceptance
      if (!isJetInAcceptance(mcpjet)) {
        continue;
      }
      if (matchedMCPJetIds.find(mcpjet.globalIndex()) == matchedMCPJetIds.end()) {
        fillMatchingResponse<Scheme, ResponseType::Miss>(mcpjet, weight);
      }
    }
  }

  // PROCESSES

  // Dummy for testing
  void processDummy(aod::JetMcCollision const&)
  {
  }
  PROCESS_SWITCH(neutralMesonInFullJets, processDummy, "dummy process", false)

  // Measuring Jet pt, eta, phi for jets in data
  void processDataJets(soa::Filtered<EMCCollisionsData>::iterator const& collision,
                       soa::Filtered<aod::FullJets> const& jets)
  {
    if (!eventSelection<DataType::Data>(collision, jets)) {
      return;
    }

    for (auto const& jet : jets) {
      registry.fill(HIST("h_data_jet_acceptance_counter"), 0.5); // all jets in accepted collisions
      if (!isJetInAcceptance(jet)) {
        continue;
      }
      registry.fill(HIST("h_data_jet_acceptance_counter"), 1.5); // only accepted jets
      fillJetHistograms<DataType::Data>(jet);
    }
  }
  PROCESS_SWITCH(neutralMesonInFullJets, processDataJets, "Full jet measurement in data", false)

  // Measuring Jet pt, eta, phi for jets in MC detector-level
  void processMCDJets(soa::Filtered<EMCCollisionsMCD>::iterator const& collision,
                      soa::Filtered<aod::FullMCDetectorLevelJets> const& jets)
  {
    const float weight = collision.weight();
    const float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));

    if (!eventSelection<DataType::MCD>(
          collision, jets, weight, pTHat)) {
      return;
    }

    for (auto const& jet : jets) {
      registry.fill(HIST("h_mcd_jet_acceptance_counter"), 0.5); // all jets in accepted collisions
      if (!isJetInAcceptance(jet)) {
        continue;
      }
      registry.fill(HIST("h_mcd_jet_acceptance_counter"), 1.5); // only accepted jets
      fillJetHistograms<DataType::MCD>(jet, weight);
    }
  }
  PROCESS_SWITCH(neutralMesonInFullJets, processMCDJets, "Full detector-level MC jet measurement", false)

  // Measuring Jet pt, eta, phi for jets in MC particle-level
  void processMCPJets(aod::JetMcCollisions::iterator const& collision,
                      soa::Filtered<aod::FullMCParticleLevelJets> const& jets,
                      soa::SmallGroups<EMCCollisionsMCD> const& detectorCollisions)
  {
    // TO DO: Do I need to add a at least one jet in collision cut or is this already taken care of?
    const float weight = collision.weight();
    const float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));

    if (!eventSelection<DataType::MCP>(
          collision, jets, weight, pTHat, &detectorCollisions)) {
      return;
    }

    for (auto const& jet : jets) {
      registry.fill(HIST("h_mcp_jet_acceptance_counter"), 0.5); // all jets in accepted collisions
      if (!isJetInAcceptance(jet)) {
        continue;
      }
      registry.fill(HIST("h_mcp_jet_acceptance_counter"), 1.5); // only accepted jets
      fillJetHistograms<DataType::MCP>(jet, weight);
    }
  }
  PROCESS_SWITCH(neutralMesonInFullJets, processMCPJets, "Full particle-level MC jet measurement", false)

  // Matching MC detector and particle-level jets to calculate the response, fake and miss histogram for later use in unfolding
  void processJetMatching(EMCCollisionsMCD::iterator const& collision, // removed soa::Filtered here, was an issue in the matching
                          JetTableMCDMatchedJoined const& mcdjets,
                          aod::JetMcCollisions const&,
                          JetTableMCPMatchedJoined const& mcpjets)
  {
    switch (matchingScheme.value) {
      case 0:
        runJetMatching<MatchingScheme::Geo>(
          collision, mcdjets, mcpjets);
        break;
      case 1:
        runJetMatching<MatchingScheme::Pt>(
          collision, mcdjets, mcpjets);
        break;
      case 2:
        runJetMatching<MatchingScheme::GeoAndPt>(
          collision, mcdjets, mcpjets);
        break;
      default:
        LOGF(fatal, "Invalud matching scheme: %d", matchingScheme.value);
    }
  }
  PROCESS_SWITCH(neutralMesonInFullJets, processJetMatching, "Jet matching for MC-based unfolding", false)
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<neutralMesonInFullJets>(cfgc)};
}
