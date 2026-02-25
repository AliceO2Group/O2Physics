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
/// \author Matthew Ockleton <matthew.ockleton@cern.ch>, University of Liverpool

#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"

#include "Common/Core/RecoDecay.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/Logger.h"
#include "Framework/runDataProcessing.h"
#include <Framework/ASoA.h>
#include <Framework/HistogramSpec.h>

#include <fairlogger/Logger.h>

#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
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
DECLARE_SOA_COLUMN(D0McOrigin, d0McOrigin, float);
DECLARE_SOA_COLUMN(D0MD, d0MD, float);
DECLARE_SOA_COLUMN(D0PtD, d0PtD, float);
DECLARE_SOA_COLUMN(D0EtaD, d0EtaD, float);
DECLARE_SOA_COLUMN(D0PhiD, d0PhiD, float);
DECLARE_SOA_COLUMN(D0Reflection, d0Reflection, int);
} // namespace d0Info

namespace jetInfo
{
// Jet
DECLARE_SOA_COLUMN(JetPt, jetPt, float);
DECLARE_SOA_COLUMN(JetEta, jetEta, float);
DECLARE_SOA_COLUMN(JetPhi, jetPhi, float);
// D0-jet
DECLARE_SOA_COLUMN(D0JetDeltaPhi, d0JetDeltaPhi, float);
} // namespace jetInfo

DECLARE_SOA_TABLE(CollisionTables, "AOD", "COLLISIONINFOTABLE",
                  o2::soa::Index<>,
                  collision::PosZ);

namespace indexColumns
{
DECLARE_SOA_INDEX_COLUMN(CollisionTable, collisionTable);
} // namespace o2::indexColumns

DECLARE_SOA_TABLE(D0DataTables, "AOD", "D0DATATABLE",
                  o2::soa::Index<>,
                  indexColumns::CollisionTableId,
                  d0Info::D0PromptBDT,
                  d0Info::D0NonPromptBDT,
                  d0Info::D0BkgBDT,
                  d0Info::D0M,
                  d0Info::D0Pt,
                  d0Info::D0Eta,
                  d0Info::D0Phi);

DECLARE_SOA_TABLE(D0McPTables, "AOD", "D0MCPARTICLELEVELTABLE",
                  o2::soa::Index<>,
                  indexColumns::CollisionTableId,
                  d0Info::D0McOrigin,
                  d0Info::D0Pt,
                  d0Info::D0Eta,
                  d0Info::D0Phi);

DECLARE_SOA_TABLE(D0McMatchedTables, "AOD", "D0MCMATCHEDTABLE",
                  o2::soa::Index<>,
                  indexColumns::CollisionTableId,
                  d0Info::D0Pt,
                  d0Info::D0Eta,
                  d0Info::D0Phi,
                  d0Info::D0McOrigin,
                  d0Info::D0Reflection);

namespace indexColumns
{
  DECLARE_SOA_INDEX_COLUMN(D0DataTable, d0Data);
  DECLARE_SOA_INDEX_COLUMN(D0McPTable, d0MCP);
  DECLARE_SOA_INDEX_COLUMN(D0McMatchedTable, d0MCMatched);
} // namespace o2::indexColumns

DECLARE_SOA_TABLE_STAGED(JetDataTables, "JETDATATABLE",
                         o2::soa::Index<>,
                         indexColumns::CollisionTableId,
                         indexColumns::D0DataTableId,
                         jetInfo::JetPt,
                         jetInfo::JetEta,
                         jetInfo::JetPhi,
                         jetInfo::D0JetDeltaPhi);

DECLARE_SOA_TABLE_STAGED(JetMCPTables, "JETMCPARTICLELEVELTABLE",
                         o2::soa::Index<>,
                         indexColumns::CollisionTableId,
                         indexColumns::D0McPTableId,
                         jetInfo::JetPt,
                         jetInfo::JetEta,
                         jetInfo::JetPhi,
                         jetInfo::D0JetDeltaPhi);

DECLARE_SOA_TABLE_STAGED(JetMCMatchedTables, "JETMCMATCHEDTABLE",
                         o2::soa::Index<>,
                         indexColumns::CollisionTableId,
                         indexColumns::D0McMatchedTableId,
                         jetInfo::JetPt,
                         jetInfo::JetEta,
                         jetInfo::JetPhi,
                         jetInfo::D0JetDeltaPhi);
} // namespace o2::aod

struct JetCorrelationD0 {
  // Define new table
  Produces<aod::CollisionTables> tableCollision;
  Produces<aod::D0DataTables> tableD0;
  Produces<aod::D0McPTables> tableD0MCParticle;
  Produces<aod::D0McMatchedTables> tableD0MCMatched;
  Produces<aod::JetDataTables> tableJet;
  Produces<aod::JetMCPTables> tableJetMCParticle;
  Produces<aod::JetMCMatchedTables> tableJetMCMatched;

  using TracksSelQuality = soa::Join<aod::TracksExtra, aod::TracksWMc>;

  // Configurables
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  // Configurable<std::string> triggerMasks{"triggerMasks", "", "possible JE Trigger masks: fJetChLowPt,fJetChHighPt,fTrackLowPt,fTrackHighPt,fJetD0ChLowPt,fJetD0ChHighPt,fJetLcChLowPt,fJetLcChHighPt,fEMCALReadout,fJetFullHighPt,fJetFullLowPt,fJetNeutralHighPt,fJetNeutralLowPt,fGammaVeryHighPtEMCAL,fGammaVeryHighPtDCAL,fGammaHighPtEMCAL,fGammaHighPtDCAL,fGammaLowPtEMCAL,fGammaLowPtDCAL,fGammaVeryLowPtEMCAL,fGammaVeryLowPtDCAL"};
  Configurable<float> jetPtCutMin{"jetPtCutMin", 5.0, "minimum value of jet pt"};
  Configurable<float> d0PtCutMin{"d0PtCutMin", 3.0, "minimum value of d0 pt"};
  Configurable<bool> doSumw{"doSumw", false, "enable sumw2 for weighted histograms"};
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> pTHatExponent{"pTHatExponent", 6.0, "exponent of the event weight for the calculation of pTHat"};
  Configurable<float> pTHatMaxMCD{"pTHatMaxMCD", 999.0, "maximum fraction of hard scattering for jet acceptance in detector MC"};
  Configurable<float> pTHatMaxMCP{"pTHatMaxMCP", 999.0, "maximum fraction of hard scattering for jet acceptance in particle MC"};
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
  void fillJetHistograms(T const& jet, U const& dphi)
  {
    registry.fill(HIST("hJetPt"), jet.pt());
    registry.fill(HIST("hJetEta"), jet.eta());
    registry.fill(HIST("hJetPhi"), jet.phi());
    registry.fill(HIST("hJet3D"), jet.pt(), jet.eta(), jet.phi());
    registry.fill(HIST("h_Jet_pT_D0_Jet_dPhi"), jet.pt(), dphi);
  }

  template <typename T, typename U>
  void applyCollisionSelections(T const& collision, U const& eventSelectionBits)
  {
    registry.fill(HIST("hCollisions"), 0.5); // All collisions
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    registry.fill(HIST("hCollisions"), 1.5); // Selected collisions
    registry.fill(HIST("hZvtxSelected"), collision.posZ());
  }

  template <typename T, typename U>
  // Jetbase is an MCD jet. We then loop through jettagv(MCP jets) to test if they match
  // void fillMatchedHistograms(T const& jetBase, float weight = 1.0) // float leadingTrackPtBase,
  void fillMatchedHistograms(T const& jetsBase, U const&, float weight = 1.0, float rho = 0.0)
  {
    for (const auto& jetBase : jetsBase) {
      if (jetBase.has_matchedJetGeo()) { // geometric matching
        for (auto& jetTag : jetBase.template matchedJetGeo_as<std::decay_t<U>>()) {
          registry.fill(HIST("hPtMatched"), jetBase.pt() - (rho * jetBase.area()), jetTag.pt(), weight);
          registry.fill(HIST("hPtMatched1d"), jetTag.pt(), weight);
          registry.fill(HIST("hPhiMatched"), jetBase.phi(), jetTag.phi(), weight);
          registry.fill(HIST("hEtaMatched"), jetBase.eta(), jetTag.eta(), weight);
          registry.fill(HIST("hPtResolution"), jetTag.pt(), (jetTag.pt() - (jetBase.pt() - (rho * jetBase.area()))) / jetTag.pt(), weight);
          registry.fill(HIST("hPhiResolution"), jetTag.pt(), jetTag.phi() - jetBase.phi(), weight);
          registry.fill(HIST("hEtaResolution"), jetTag.pt(), jetTag.eta() - jetBase.eta(), weight);
        }
      }
    }
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
    registry.add("h_Jet_pT_D0_Jet_dPhi", "p_{T, jet} vs #Delta #phi _{D^{0}, jet}", kTH2F, {{100, 0, 100}, {100, 0, o2::constants::math::TwoPI}});

    // Matching histograms
    registry.add("hPtMatched", "p_{T} matching;p_{T,det};p_{T,part}", {HistType::kTH2F, {{500, -100, 400}, {400, 0, 400}}}, doSumw);
    registry.add("hPtMatched1d", "p_{T} matching 1d;p_{T,part}", {HistType::kTH1F, {{400, 0, 400}}}, doSumw);
    registry.add("hPhiMatched", "#phi matching;#phi_{det};#phi_{part}", {HistType::kTH2F, {{100, 0.0, o2::constants::math::TwoPI}, {100, 0.0, o2::constants::math::TwoPI}}}, doSumw);
    registry.add("hEtaMatched", "#eta matching;#eta_{det};#eta_{part}", {HistType::kTH2F, {{100, -1, 1}, {100, -1, 1}}}, doSumw);
    registry.add("hPtResolution", "p_{T} resolution;p_{T,part};Relative Resolution", {HistType::kTH2F, {{400, 0, 400}, {1000, -5.0, 5.0}}}, doSumw);
    registry.add("hPhiResolution", "#phi resolution;#p_{T,part};Resolution", {HistType::kTH2F, {{400, 0, 400}, {1000, -7.0, 7.0}}}, doSumw);
    registry.add("hEtaResolution", "#eta resolution;#p_{T,part};Resolution", {HistType::kTH2F, {{400, 0, 400}, {1000, -1.0, 1.0}}}, doSumw);
  }

  void processData(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                   aod::CandidatesD0Data const& d0Candidates,
                   soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets)
  {
    applyCollisionSelections(collision, eventSelectionBits);
    tableCollision(collision.posZ());

    for (const auto& d0Candidate : d0Candidates) {
      const auto scores = d0Candidate.mlScores();
      if (d0Candidate.pt() < d0PtCutMin) {
        continue;
      }
      fillD0Histograms(d0Candidate, scores);
      tableD0(tableCollision.lastIndex(),
              scores[2],
              scores[1],
              scores[0],
              d0Candidate.m(),
              d0Candidate.pt(),
              d0Candidate.eta(),
              d0Candidate.phi());
      for (const auto& jet : jets) {
        if (jet.pt() < jetPtCutMin) {
          continue;
        }
        float dphi = RecoDecay::constrainAngle(jet.phi() - d0Candidate.phi());
        if (abs(dphi - M_PI) > (M_PI / 2)) {
          continue;
        }
        fillJetHistograms(jet, dphi);
        tableJet(tableCollision.lastIndex(),
                 tableD0.lastIndex(),
                 jet.pt(),
                 jet.eta(),
                 jet.phi(),
                 dphi);
      }
    }
  }
  PROCESS_SWITCH(JetCorrelationD0, processData, "charged particle level jet analysis", true);

  void processMCParticle(aod::JetMcCollision const& collision,
                         aod::CandidatesD0MCP const& d0MCPCandidates,
                         soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents>> const& jets)
  {
    // applyCollisionSelections(collision, eventSelectionBits);
    // Need to select on pz
    tableCollision(collision.posZ());

    for (const auto& d0MCPCandidate : d0MCPCandidates) {
      if (d0MCPCandidate.pt() < d0PtCutMin) {
        continue;
      }
      tableD0MCParticle(tableCollision.lastIndex(),
                        d0MCPCandidate.originMcGen(),
                        d0MCPCandidate.pt(),
                        d0MCPCandidate.eta(),
                        d0MCPCandidate.phi());

      for (const auto& jet : jets) {
        if (jet.pt() < jetPtCutMin) {
          continue;
        }
        float dphi = RecoDecay::constrainAngle(jet.phi() - d0MCPCandidate.phi());
        if (abs(dphi - M_PI) > (M_PI / 2)) {
          continue;
        }
        fillJetHistograms(jet, dphi);
        tableJetMCParticle(tableCollision.lastIndex(),
                           tableD0MCParticle.lastIndex(),
                           jet.pt(),
                           jet.eta(),
                           jet.phi(),
                           dphi);
      }
    }
  }
  PROCESS_SWITCH(JetCorrelationD0, processMCParticle, "process MC Particle jets", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetCorrelationD0>(cfgc)};
}
