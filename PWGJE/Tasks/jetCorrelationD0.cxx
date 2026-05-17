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
/// \file jetCorrelationD0.cxx
/// \brief Task for analysing D0 triggered jet events.
/// \author Matthew Ockleton matthew.ockleton@cern.ch, University of Liverpool

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetHFUtilities.h"
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

#include <cstdlib>
#include <string>
#include <vector>

#include <math.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{

namespace d0collisionInfo
{
DECLARE_SOA_COLUMN(PosZ, posZ, float);
} // namespace d0collisionInfo

DECLARE_SOA_TABLE(CollisionTables, "AOD", "COLLINFOTABLE",
                  o2::soa::Index<>,
                  d0collisionInfo::PosZ);

DECLARE_SOA_TABLE(McCollisionTables, "AOD", "MCCOLLINFOTABLE",
                  o2::soa::Index<>,
                  d0collisionInfo::PosZ);

DECLARE_SOA_TABLE(MatchCollTables, "AOD", "MATCHCOLLTABLE",
                  o2::soa::Index<>,
                  d0collisionInfo::PosZ);

namespace collisionInfo
{
DECLARE_SOA_INDEX_COLUMN_CUSTOM(CollisionTable, collisionTable, "COLLINFOTABLES");
DECLARE_SOA_INDEX_COLUMN_CUSTOM(McCollisionTable, mcCollisionTable, "MCCOLLINFOTABLES");
DECLARE_SOA_INDEX_COLUMN_CUSTOM(MatchCollTable, matchCollTable, "MATCHCOLLTABLES");
} // namespace collisionInfo
namespace d0Info
{
// D0
DECLARE_SOA_COLUMN(D0PromptBDT, d0PromptBDT, float);
DECLARE_SOA_COLUMN(D0NonPromptBDT, d0NonPromptBDT, float);
DECLARE_SOA_COLUMN(D0BkgBDT, d0BkgBDT, float);
DECLARE_SOA_COLUMN(D0M, d0M, float);
DECLARE_SOA_COLUMN(D0Pt, d0Pt, float);
DECLARE_SOA_COLUMN(D0Eta, d0Eta, float);
DECLARE_SOA_COLUMN(D0Phi, d0Phi, float);
DECLARE_SOA_COLUMN(D0Y, d0Y, float);
DECLARE_SOA_COLUMN(D0McOrigin, d0McOrigin, float);
DECLARE_SOA_COLUMN(D0MD, d0MD, float);
DECLARE_SOA_COLUMN(D0PtD, d0PtD, float);
DECLARE_SOA_COLUMN(D0EtaD, d0EtaD, float);
DECLARE_SOA_COLUMN(D0PhiD, d0PhiD, float);
DECLARE_SOA_COLUMN(D0Reflection, d0Reflection, int);
} // namespace d0Info

DECLARE_SOA_TABLE(D0Tables, "AOD", "D0TABLE",
                  o2::soa::Index<>,
                  collisionInfo::CollisionTableId,
                  d0Info::D0PromptBDT,
                  d0Info::D0NonPromptBDT,
                  d0Info::D0BkgBDT,
                  d0Info::D0M,
                  d0Info::D0Pt,
                  d0Info::D0Eta,
                  d0Info::D0Phi,
                  d0Info::D0Y);

DECLARE_SOA_TABLE(D0McPTables, "AOD", "D0MCPTABLE",
                  o2::soa::Index<>,
                  collisionInfo::McCollisionTableId,
                  d0Info::D0McOrigin,
                  d0Info::D0Pt,
                  d0Info::D0Eta,
                  d0Info::D0Phi,
                  d0Info::D0Y);

namespace jetInfo
{
// D0 tables
DECLARE_SOA_INDEX_COLUMN(D0Table, d0Table);
DECLARE_SOA_INDEX_COLUMN(D0McPTable, d0McPTable);
// Jet
DECLARE_SOA_COLUMN(JetPt, jetPt, float);
DECLARE_SOA_COLUMN(JetEta, jetEta, float);
DECLARE_SOA_COLUMN(JetPhi, jetPhi, float);
DECLARE_SOA_COLUMN(PJetPt, pJetPt, float);
DECLARE_SOA_COLUMN(PJetEta, pJetEta, float);
DECLARE_SOA_COLUMN(PJetPhi, pJetPhi, float);
// D0-jet
DECLARE_SOA_COLUMN(D0JetDeltaPhi, d0JetDeltaPhi, float);
DECLARE_SOA_COLUMN(D0JetDeltaPhiP, d0JetDeltaPhiP, float);
} // namespace jetInfo

DECLARE_SOA_TABLE_STAGED(JetTables, "JETTABLE",
                         o2::soa::Index<>,
                         collisionInfo::CollisionTableId,
                         jetInfo::D0TableId,
                         jetInfo::JetPt,
                         jetInfo::JetEta,
                         jetInfo::JetPhi,
                         jetInfo::D0JetDeltaPhi);

DECLARE_SOA_TABLE_STAGED(JetMcPTables, "JETMCPTABLE",
                         o2::soa::Index<>,
                         collisionInfo::McCollisionTableId,
                         jetInfo::D0McPTableId,
                         jetInfo::JetPt,
                         jetInfo::JetEta,
                         jetInfo::JetPhi,
                         jetInfo::D0JetDeltaPhiP);

DECLARE_SOA_TABLE_STAGED(JetMatchedTables, "JETMATCHEDTABLE",
                         o2::soa::Index<>,
                         collisionInfo::MatchCollTableId,
                         jetInfo::JetPt,
                         jetInfo::JetEta,
                         jetInfo::JetPhi,
                         jetInfo::PJetPt,
                         jetInfo::PJetEta,
                         jetInfo::PJetPhi,
                         jetInfo::D0JetDeltaPhi,
                         jetInfo::D0JetDeltaPhiP);

} // namespace o2::aod

struct JetCorrelationD0 {
  // Define new table
  Produces<aod::CollisionTables> tableCollision;
  Produces<aod::MatchCollTables> tableMatchedCollision;
  Produces<aod::McCollisionTables> tableMcCollision;
  Produces<aod::D0Tables> tableD0;
  Produces<aod::D0McPTables> tableD0McParticle;
  Produces<aod::JetTables> tableJet;
  Produces<aod::JetMcPTables> tableJetMcParticle;
  Produces<aod::JetMatchedTables> tableJetMatched;

  // Configurables
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<bool> skipMBGapEvents{"skipMBGapEvents", false, "decide to run over MB gap events or not"};
  Configurable<bool> applyRCTSelections{"applyRCTSelections", true, "decide to apply RCT selections"};
  Configurable<float> jetPtCutMin{"jetPtCutMin", 5.0, "minimum value of jet pt"};
  Configurable<float> d0PtCutMin{"d0PtCutMin", 1.0, "minimum value of d0 pt"};
  Configurable<float> jetMcPtCutMin{"jetMcPtCutMin", 3.0, "minimum value of jet pt particle level"};
  Configurable<float> d0McPtCutMin{"d0McPtCutMin", 0.5, "minimum value of d0 pt particle level"};
  Configurable<float> vertexZCut{"vertexZCut", 10.0, "Accepted z-vertex range"};
  Configurable<float> pTHatExponent{"pTHatExponent", 6.0, "exponent of the event weight for the calculation of pTHat"};
  Configurable<float> pTHatMaxMcD{"pTHatMaxMcD", 999.0, "maximum fraction of hard scattering for jet acceptance in detector MC"};
  Configurable<float> pTHatMaxMcP{"pTHatMaxMcP", 999.0, "maximum fraction of hard scattering for jet acceptance in particle MC"};
  Configurable<float> pTHatAbsoluteMin{"pTHatAbsoluteMin", -99.0, "minimum value of pTHat"};

  // Filters
  Filter eventCuts = (nabs(aod::jcollision::posZ) < vertexZCut);
  std::vector<int> eventSelectionBits;

  // Histograms
  HistogramRegistry registry{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  template <typename T, typename U>
  void fillD0Histograms(T const& d0, U const& scores)
  {
    registry.fill(HIST("hD0MlBkg"), scores[0]);
    registry.fill(HIST("hD0MlNonPrompt"), scores[1]);
    registry.fill(HIST("hD0MlPrompt"), scores[2]);

    registry.fill(HIST("hD0Pt"), d0.pt());
    registry.fill(HIST("hD0M"), d0.m());
    registry.fill(HIST("hD0Eta"), d0.eta());
    registry.fill(HIST("hD0Phi"), d0.phi());
  }
  template <typename T, typename U>
  void fillJetHistograms(T const& jet, U const& dPhi)
  {
    registry.fill(HIST("hJetPt"), jet.pt());
    registry.fill(HIST("hJetEta"), jet.eta());
    registry.fill(HIST("hJetPhi"), jet.phi());
    registry.fill(HIST("hJet3D"), jet.pt(), jet.eta(), jet.phi());
    registry.fill(HIST("h_Jet_D0_Jet_dPhi"), dPhi);
    registry.fill(HIST("h_Jet_pT_D0_Jet_dPhi"), jet.pt(), dPhi);
  }

  template <typename T>
  bool applyCollisionSelections(T const& collision)
  {
    registry.fill(HIST("hCollisions"), 0.5); // All collisions
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents, applyRCTSelections)) {
      return false;
    }
    registry.fill(HIST("hCollisions"), 1.5); // Selected collisions
    registry.fill(HIST("hZvtxSelected"), collision.posZ());
    return true;
  }
  void init(InitContext const&)
  {
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    // General Axes
    AxisSpec axisEta = {100, -1.0, 1.0, "#eta"};
    AxisSpec axisPhi = {100, 0.0, o2::constants::math::TwoPI, "#phi"};
    AxisSpec axisInvMass = {500, 0, 10, "M (GeV/c)"};

    // General Histograms
    registry.add("hCollisions", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
    registry.add("hZvtxSelected", "Z vertex position;Z_{vtx};entries", {HistType::kTH1F, {{80, -20, 20}}});

    // D0 Histograms
    registry.add("hD0MlPrompt", "D0 ML Prompt Scores", {HistType::kTH1F, {{100, -1.0, 2.0}}});
    registry.add("hD0MlNonPrompt", "D0 ML NonPrompt Scores", {HistType::kTH1F, {{100, -1.0, 2.0}}});
    registry.add("hD0MlBkg", "D0 ML Background Scores", {HistType::kTH1F, {{100, -1.0, 2.0}}});

    registry.add("hD0Pt", "D^{0} p_{T};p_{T}^{D^{0}} (GeV/c);entries", {HistType::kTH1F, {{500, -100, 400, "p_{T}^{D^{0}} (GeV/c)"}}});
    registry.add("hD0M", "D^{0} Mass;M_{#pi K} (GeV/c);entries", HistType::kTH1F, {axisInvMass});
    registry.add("hD0Eta", "D^{0} #eta ;#eta_{D^{0}};entries", HistType::kTH1F, {axisEta});
    registry.add("hD0Phi", "D^{0} #phi ;#phi_{D^{0}};entries", HistType::kTH1F, {axisPhi});

    // Jet Histograms
    registry.add("hJetPt", "jet p_{T};p_{T,jet};entries", {HistType::kTH1F, {{500, -100, 400}}});
    registry.add("hJetEta", "jet #eta;#eta_{jet};entries", HistType::kTH1F, {axisEta});
    registry.add("hJetPhi", "jet #phi;#phi_{jet};entries", HistType::kTH1F, {axisPhi});
    registry.add("hJet3D", "3D jet distribution;p_{T};#eta;#phi", {HistType::kTH3F, {{500, -100, 400}, {100, -1.0, 1.0}, {100, 0.0, o2::constants::math::TwoPI}}});
    registry.add("h_Jet_D0_Jet_dPhi", "#Delta #phi _{D^{0}, jet}", kTH1F, {{100, 0, o2::constants::math::TwoPI}});
    registry.add("h_Jet_pT_D0_Jet_dPhi", "p_{T, jet} vs #Delta #phi _{D^{0}, jet}", kTH2F, {{100, 0, 100}, {100, 0, o2::constants::math::TwoPI}});

    // Matching histograms
    registry.add("hPtMatched", "p_{T} matching;p_{T,det};p_{T,part}", {HistType::kTH2F, {{500, -100, 400}, {400, 0, 400}}});
    registry.add("hPtMatched1d", "p_{T} matching 1d;p_{T,part}", {HistType::kTH1F, {{400, 0, 400}}});
    registry.add("hPhiMatched", "#phi matching;#phi_{det};#phi_{part}", {HistType::kTH2F, {{100, 0.0, o2::constants::math::TwoPI}, {100, 0.0, o2::constants::math::TwoPI}}});
    registry.add("hEtaMatched", "#eta matching;#eta_{det};#eta_{part}", {HistType::kTH2F, {{100, -1, 1}, {100, -1, 1}}});
    registry.add("hPtResolution", "p_{T} resolution;p_{T,part};Relative Resolution", {HistType::kTH2F, {{400, 0, 400}, {1000, -5.0, 5.0}}});
    registry.add("hPhiResolution", "#phi resolution;#p_{T,part};Resolution", {HistType::kTH2F, {{400, 0, 400}, {1000, -7.0, 7.0}}});
    registry.add("hEtaResolution", "#eta resolution;#p_{T,part};Resolution", {HistType::kTH2F, {{400, 0, 400}, {1000, -1.0, 1.0}}});
  }
  void processData(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                   aod::CandidatesD0Data const& d0Candidates,
                   soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets)
  {
    if (!applyCollisionSelections(collision)) {
      return;
    }
    tableCollision(collision.posZ());
    for (const auto& d0Candidate : d0Candidates) {
      if (d0Candidate.pt() < d0PtCutMin) { // once settled on a mlcut, then add the lower bound of the systematics as a cut here
        continue;
      }
      const auto scores = d0Candidate.mlScores();
      fillD0Histograms(d0Candidate, scores);
      tableD0(tableCollision.lastIndex(),
              scores[2],
              scores[1],
              scores[0],
              d0Candidate.m(),
              d0Candidate.pt(),
              d0Candidate.eta(),
              d0Candidate.phi(),
              d0Candidate.y());
      for (const auto& jet : jets) {
        if (jet.pt() < jetPtCutMin) {
          continue;
        }
        float dPhi = RecoDecay::constrainAngle(jet.phi() - d0Candidate.phi(), -o2::constants::math::PI);
        if (std::abs(dPhi - o2::constants::math::PI) > (o2::constants::math::PI / 2)) {
          continue;
        }
        fillJetHistograms(jet, dPhi);
        tableJet(tableCollision.lastIndex(),
                 tableD0.lastIndex(),
                 jet.pt(),
                 jet.eta(),
                 jet.phi(),
                 dPhi);
      }
    }
  }
  PROCESS_SWITCH(JetCorrelationD0, processData, "charged particle level jet analysis", true);

  void processMcDetector(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                         aod::CandidatesD0MCD const& d0Candidates,
                         soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents> const& jets)
  {
    if (!applyCollisionSelections(collision)) {
      return;
    }
    tableCollision(collision.posZ());
    for (const auto& d0Candidate : d0Candidates) {
      if (d0Candidate.pt() < d0PtCutMin) { // once settled on a mlcut, then add the lower biund of the systematics as a cut here
        continue;
      }
      const auto scores = d0Candidate.mlScores();
      fillD0Histograms(d0Candidate, scores);
      tableD0(tableCollision.lastIndex(), // might want to add some more detector level D0 quantities like prompt or non prompt info
              scores[2],
              scores[1],
              scores[0],
              d0Candidate.m(),
              d0Candidate.pt(),
              d0Candidate.eta(),
              d0Candidate.phi(),
              d0Candidate.y());
      for (const auto& jet : jets) {
        if (jet.pt() < jetPtCutMin) {
          continue;
        }
        float dPhi = RecoDecay::constrainAngle(jet.phi() - d0Candidate.phi(), -o2::constants::math::PI);
        if (std::abs(dPhi - o2::constants::math::PI) > (o2::constants::math::PI / 2)) {
          continue;
        }
        fillJetHistograms(jet, dPhi);
        tableJet(tableCollision.lastIndex(),
                 tableD0.lastIndex(),
                 jet.pt(),
                 jet.eta(),
                 jet.phi(),
                 dPhi);
      }
    }
  }
  PROCESS_SWITCH(JetCorrelationD0, processMcDetector, "charged detector level jet analysis", false);

  void processMcParticle(aod::JetMcCollision const& collision,
                         aod::CandidatesD0MCP const& d0McPCandidates,
                         soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents> const& jets)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents, applyRCTSelections)) {
      return;
    }
    tableMcCollision(collision.posZ());
    for (const auto& d0McPCandidate : d0McPCandidates) {
      if (d0McPCandidate.pt() < d0McPtCutMin) {
        continue;
      }
      tableD0McParticle(tableMcCollision.lastIndex(),
                        d0McPCandidate.originMcGen(),
                        d0McPCandidate.pt(),
                        d0McPCandidate.eta(),
                        d0McPCandidate.phi(),
                        d0McPCandidate.y());

      for (const auto& jet : jets) {
        if (jet.pt() < jetMcPtCutMin) {
          continue;
        }
        float dPhi = RecoDecay::constrainAngle(jet.phi() - d0McPCandidate.phi(), -o2::constants::math::PI);
        if (std::abs(dPhi - o2::constants::math::PI) > (o2::constants::math::PI / 2)) {
          continue;
        }
        fillJetHistograms(jet, dPhi);
        tableJetMcParticle(tableMcCollision.lastIndex(),
                           tableD0McParticle.lastIndex(),
                           jet.pt(),
                           jet.eta(),
                           jet.phi(),
                           dPhi);
      }
    }
  }
  PROCESS_SWITCH(JetCorrelationD0, processMcParticle, "charged MC Particle jets", false);

  void processMcMatched(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                        aod::CandidatesD0MCD const& d0Candidates,
                        aod::JetTracksMCD const& tracks,
                        aod::JetParticles const& particles,
                        soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets> const& McDJets,
                        aod::ChargedMCParticleLevelJets const&)
  {
    if (!applyCollisionSelections(collision)) {
      return;
    }
    tableMatchedCollision(collision.posZ());
    for (const auto& d0Candidate : d0Candidates) {
      if (d0Candidate.pt() < d0PtCutMin) { // once settled on a mlcut, then add the lower bound of the systematics as a cut here
        continue;
      }
      bool isMatched = false;
      const auto& d0Particle = jethfutilities::matchedHFParticle(d0Candidate, tracks, particles, isMatched);
      if (!isMatched) {
        continue;
      }
      for (const auto& McDJet : McDJets) {
        if (McDJet.pt() < jetPtCutMin) {
          continue;
        }
        float dPhiD = RecoDecay::constrainAngle(McDJet.phi() - d0Candidate.phi(), -o2::constants::math::PI);
        if (std::abs(dPhiD - o2::constants::math::PI) > (o2::constants::math::PI / 2)) {
          continue;
        }
        if (McDJet.has_matchedJetGeo()) { // geometric matching
          for (auto const& McPJet : McDJet.template matchedJetGeo_as<aod::ChargedMCParticleLevelJets>()) {
            float dPhiP = RecoDecay::constrainAngle(McPJet.phi() - d0Particle.phi(), -o2::constants::math::PI);
            // if (std::abs(dPhiP - o2::constants::math::PI) > (o2::constants::math::PI / 2)) {
            //   continue;
            // }
            tableJetMatched(tableMatchedCollision.lastIndex(),
                            McDJet.pt(),
                            McDJet.eta(),
                            McDJet.phi(),
                            McPJet.pt(),
                            McPJet.eta(),
                            McPJet.phi(),
                            dPhiD,
                            dPhiP);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(JetCorrelationD0, processMcMatched, "process matching of particle level jets to detector level jets", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetCorrelationD0>(cfgc)};
}
