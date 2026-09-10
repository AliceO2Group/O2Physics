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
/// \file jetCorrelationDMeson.cxx

/// \brief Task for analysing D-meson (D0, D*) triggered jet events

/// \author Matthew Ockleton matthew.ockleton@cern.ch, University of Liverpool
/// \author Mokshi Vaid  mokshi.vaid@cern.ch, University of Jammu

#include "PWGHF/Core/DecayChannels.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetHFUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <Rtypes.h>

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <type_traits>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace hf2Prong = o2::hf_decay::hf_cand_2prong;
namespace hfDstar = o2::hf_decay::hf_cand_dstar;

namespace o2::aod
{

namespace dcollisionInfo
{
DECLARE_SOA_COLUMN(PosZ, posZ, float);
} // namespace dcollisionInfo

DECLARE_SOA_TABLE(CollisionTables, "AOD", "COLLINFOTABLE",
                  o2::soa::Index<>,
                  dcollisionInfo::PosZ);

DECLARE_SOA_TABLE(McCollisionTables, "AOD", "MCCOLLINFOTABLE",
                  o2::soa::Index<>,
                  dcollisionInfo::PosZ);

DECLARE_SOA_TABLE(MatchCollTables, "AOD", "MATCHCOLLTABLE",
                  o2::soa::Index<>,
                  dcollisionInfo::PosZ);

namespace collisionInfo
{
DECLARE_SOA_INDEX_COLUMN_CUSTOM(CollisionTable, collisionTable, "COLLINFOTABLES");
DECLARE_SOA_INDEX_COLUMN_CUSTOM(McCollisionTable, mcCollisionTable, "MCCOLLINFOTABLES");
DECLARE_SOA_INDEX_COLUMN_CUSTOM(MatchCollTable, matchCollTable, "MATCHCOLLTABLES");
} // namespace collisionInfo

namespace dInfo
{
DECLARE_SOA_COLUMN(DPromptBDT, dPromptBDT, float);
DECLARE_SOA_COLUMN(DNonPromptBDT, dNonPromptBDT, float);
DECLARE_SOA_COLUMN(DBkgBDT, dBkgBDT, float);
DECLARE_SOA_COLUMN(DM, dM, float);
DECLARE_SOA_COLUMN(DDeltaM, dDeltaM, float); // -1 for D0 ("not applicable"); M(D*)-M(D0) for D*
DECLARE_SOA_COLUMN(DPt, dPt, float);
DECLARE_SOA_COLUMN(DEta, dEta, float);
DECLARE_SOA_COLUMN(DPhi, dPhi, float);
DECLARE_SOA_COLUMN(DY, dY, float);
DECLARE_SOA_COLUMN(DMcOrigin, dMcOrigin, float);
DECLARE_SOA_COLUMN(DCategory, dCategory, int);
DECLARE_SOA_COLUMN(DDecayChannel, dDecayChannel, int8_t);
} // namespace dInfo

DECLARE_SOA_TABLE(DTables, "AOD", "DTABLE",
                  o2::soa::Index<>,
                  collisionInfo::CollisionTableId,
                  dInfo::DPromptBDT,
                  dInfo::DNonPromptBDT,
                  dInfo::DBkgBDT,
                  dInfo::DM,
                  dInfo::DDeltaM,
                  dInfo::DPt,
                  dInfo::DEta,
                  dInfo::DPhi,
                  dInfo::DY);

// MC-detector level: shared table, category column
DECLARE_SOA_TABLE(DMcDTables, "AOD", "DMCDTABLE",
                  o2::soa::Index<>,
                  collisionInfo::McCollisionTableId,
                  dInfo::DPromptBDT,
                  dInfo::DNonPromptBDT,
                  dInfo::DBkgBDT,
                  dInfo::DM,
                  dInfo::DDeltaM,
                  dInfo::DPt,
                  dInfo::DEta,
                  dInfo::DPhi,
                  dInfo::DY,
                  dInfo::DCategory);

// MC-particle level
DECLARE_SOA_TABLE(DMcPTables, "AOD", "DMCPTABLE",
                  o2::soa::Index<>,
                  collisionInfo::McCollisionTableId,
                  dInfo::DMcOrigin,
                  dInfo::DPt,
                  dInfo::DEta,
                  dInfo::DPhi,
                  dInfo::DY,
                  dInfo::DDecayChannel);

namespace jetInfo
{
// Shared D-meson table indices (were D0Table/DstarTable etc. separately)
DECLARE_SOA_INDEX_COLUMN(DTable, dTable);
DECLARE_SOA_INDEX_COLUMN(DMcDTable, dMcDTable);
DECLARE_SOA_INDEX_COLUMN(DMcPTable, dMcPTable);
// Jet
DECLARE_SOA_COLUMN(JetPt, jetPt, float);
DECLARE_SOA_COLUMN(JetEta, jetEta, float);
DECLARE_SOA_COLUMN(JetPhi, jetPhi, float);
DECLARE_SOA_COLUMN(PJetPt, pJetPt, float);
DECLARE_SOA_COLUMN(PJetEta, pJetEta, float);
DECLARE_SOA_COLUMN(PJetPhi, pJetPhi, float);
// D-jet
DECLARE_SOA_COLUMN(DJetDeltaPhi, dJetDeltaPhi, float);
DECLARE_SOA_COLUMN(DJetDeltaPhiP, dJetDeltaPhiP, float);
} // namespace jetInfo

DECLARE_SOA_TABLE_STAGED(JetTables, "JETTABLE",
                         o2::soa::Index<>,
                         collisionInfo::CollisionTableId,
                         jetInfo::DTableId,
                         jetInfo::JetPt,
                         jetInfo::JetEta,
                         jetInfo::JetPhi,
                         jetInfo::DJetDeltaPhi);

DECLARE_SOA_TABLE_STAGED(JetMcDTables, "JETMCDTABLE",
                         o2::soa::Index<>,
                         collisionInfo::CollisionTableId,
                         jetInfo::DMcDTableId,
                         jetInfo::JetPt,
                         jetInfo::JetEta,
                         jetInfo::JetPhi,
                         jetInfo::DJetDeltaPhi);

DECLARE_SOA_TABLE_STAGED(JetMcPTables, "JETMCPTABLE",
                         o2::soa::Index<>,
                         collisionInfo::McCollisionTableId,
                         jetInfo::DMcPTableId,
                         jetInfo::JetPt,
                         jetInfo::JetEta,
                         jetInfo::JetPhi,
                         jetInfo::DJetDeltaPhiP);

DECLARE_SOA_TABLE_STAGED(JetMatchedTables, "JETMATCHEDTABLE",
                         o2::soa::Index<>,
                         collisionInfo::MatchCollTableId,
                         jetInfo::JetPt,
                         jetInfo::JetEta,
                         jetInfo::JetPhi,
                         jetInfo::PJetPt,
                         jetInfo::PJetEta,
                         jetInfo::PJetPhi,
                         jetInfo::DJetDeltaPhi,
                         jetInfo::DJetDeltaPhiP);

} // namespace o2::aod

struct JetCorrelationDmeson {
  // Define new tables
  Produces<aod::CollisionTables> tableCollision;
  Produces<aod::McCollisionTables> tableMcCollision;
  Produces<aod::MatchCollTables> tableMatchedCollision;
  Produces<aod::DTables> tableD;
  Produces<aod::DMcDTables> tableDMcDetector;
  Produces<aod::DMcPTables> tableDMcParticle;
  Produces<aod::JetTables> tableJet;
  Produces<aod::JetMcDTables> tableJetMcDetector;
  Produces<aod::JetMcPTables> tableJetMcParticle;
  Produces<aod::JetMatchedTables> tableJetMatched;

  // Configurables
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<bool> skipMBGapEvents{"skipMBGapEvents", false, "decide to run over MB gap events or not"};
  Configurable<bool> applyRCTSelections{"applyRCTSelections", true, "decide to apply RCT selections"};
  Configurable<float> jetPtCutMin{"jetPtCutMin", 5.0, "minimum value of jet pt"};
  Configurable<float> dPtCutMin{"dPtCutMin", 1.0, "minimum value of D-meson pt"};
  Configurable<float> jetMcPtCutMin{"jetMcPtCutMin", 3.0, "minimum value of jet pt particle level"};
  Configurable<float> dMcPtCutMin{"dMcPtCutMin", 0.5, "minimum value of D-meson pt particle level"};
  Configurable<float> vertexZCut{"vertexZCut", 10.0, "Accepted z-vertex range"};

  // Filters
  Filter eventCuts = (nabs(aod::jcollision::posZ) < vertexZCut);
  std::vector<int> eventSelectionBits;

  // Histograms
  HistogramRegistry registry{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // D0 and D* get their OWN distinct numeric values (not shared), since the
  // classification is tied to species-specific decay-channel information
  enum DCategory : int {
    Undefined = -1, // no truth match / unclassified
    // D0
    Signal = 0,     // correctly identified D0 and (bar), π+ K− + cc
    Reflection = 1, // true D0(bar) reconstructed with swapped mass hypothesis
    CorrBkg1 = 2,   // correlated background: π+ K− π0
    CorrBkg2 = 3,   // correlated background: π+ π−
    CorrBkg3 = 4,   // correlated background: π+ π− π0
    CorrBkg4 = 5,   // correlated background: K+ K−
    // D*
    DstarSignal = 6,     // correctly identified D* and (bar), (D0 → π+ K−) π+ + cc
    DstarReflection = 7, // true D*(bar) reconstructed with swapped mass hypothesis
    DstarCorrBkg = 8     // correlated background: (D0 → π+ K− π0) π+
  };

  template <typename T, typename U>
  void fillDHistograms(T const& dCandidate, U const& scores)
  {
    registry.fill(HIST("hDMlBkg"), scores[0]);
    registry.fill(HIST("hDMlNonPrompt"), scores[1]);
    registry.fill(HIST("hDMlPrompt"), scores[2]);

    registry.fill(HIST("hDPt"), dCandidate.pt());
    registry.fill(HIST("hDM"), dCandidate.m());
    registry.fill(HIST("hDEta"), dCandidate.eta());
    registry.fill(HIST("hDPhi"), dCandidate.phi());
    // (D0 candidates have no invMassCharm()).
    if constexpr (std::is_same_v<T, aod::CandidatesDstarData::iterator> ||
                  std::is_same_v<T, aod::CandidatesDstarMCD::iterator>) {
      registry.fill(HIST("hDDeltaM"), dCandidate.m() - dCandidate.invMassCharm());
    }
  }

  template <typename T, typename U>
  void fillJetHistograms(T const& jet, U const& dPhi)
  {
    registry.fill(HIST("hJetPt"), jet.pt());
    registry.fill(HIST("hJetEta"), jet.eta());
    registry.fill(HIST("hJetPhi"), jet.phi());
    registry.fill(HIST("hJet3D"), jet.pt(), jet.eta(), jet.phi());
    registry.fill(HIST("h_Jet_D_Jet_dPhi"), dPhi);
    registry.fill(HIST("h_Jet_pT_D_Jet_dPhi"), jet.pt(), dPhi);
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

    // D-meson Histograms (shared)
    registry.add("hDMlPrompt", "D ML Prompt Scores", {HistType::kTH1F, {{100, -1.0, 2.0}}});
    registry.add("hDMlNonPrompt", "D ML NonPrompt Scores", {HistType::kTH1F, {{100, -1.0, 2.0}}});
    registry.add("hDMlBkg", "D ML Background Scores", {HistType::kTH1F, {{100, -1.0, 2.0}}});

    registry.add("hDPt", "D p_{T};p_{T}^{D} (GeV/c);entries", {HistType::kTH1F, {{500, -100, 400, "p_{T}^{D} (GeV/c)"}}});
    registry.add("hDM", "D Mass;M (GeV/c);entries", HistType::kTH1F, {axisInvMass});
    registry.add("hDEta", "D #eta ;#eta_{D};entries", HistType::kTH1F, {axisEta});
    registry.add("hDPhi", "D #phi ;#phi_{D};entries", HistType::kTH1F, {axisPhi});

    // Delta-m: only ever filled for D*
    registry.add("hDDeltaM", "D* - D0 Mass Difference;#Delta M (GeV/c);entries", {HistType::kTH1F, {{300, 0.13, 0.17}}});

    // Jet Histograms
    registry.add("hJetPt", "jet p_{T};p_{T,jet};entries", {HistType::kTH1F, {{500, -100, 400}}});
    registry.add("hJetEta", "jet #eta;#eta_{jet};entries", HistType::kTH1F, {axisEta});
    registry.add("hJetPhi", "jet #phi;#phi_{jet};entries", HistType::kTH1F, {axisPhi});
    registry.add("hJet3D", "3D jet distribution;p_{T};#eta;#phi", {HistType::kTH3F, {{500, -100, 400}, {100, -1.0, 1.0}, {100, 0.0, o2::constants::math::TwoPI}}});
    registry.add("h_Jet_D_Jet_dPhi", "#Delta #phi _{D, jet}", kTH1F, {{100, 0, o2::constants::math::TwoPI}});
    registry.add("h_Jet_pT_D_Jet_dPhi", "p_{T, jet} vs #Delta #phi _{D, jet}", kTH2F, {{100, 0, 100}, {100, 0, o2::constants::math::TwoPI}});

    // Matching histograms
    registry.add("hPtMatched", "p_{T} matching;p_{T,det};p_{T,part}", {HistType::kTH2F, {{500, -100, 400}, {400, 0, 400}}});
    registry.add("hPtMatched1d", "p_{T} matching 1d;p_{T,part}", {HistType::kTH1F, {{400, 0, 400}}});
    registry.add("hPhiMatched", "#phi matching;#phi_{det};#phi_{part}", {HistType::kTH2F, {{100, 0.0, o2::constants::math::TwoPI}, {100, 0.0, o2::constants::math::TwoPI}}});
    registry.add("hEtaMatched", "#eta matching;#eta_{det};#eta_{part}", {HistType::kTH2F, {{100, -1, 1}, {100, -1, 1}}});
    registry.add("hPtResolution", "p_{T} resolution;p_{T,part};Relative Resolution", {HistType::kTH2F, {{400, 0, 400}, {1000, -5.0, 5.0}}});
    registry.add("hPhiResolution", "#phi resolution;#p_{T,part};Resolution", {HistType::kTH2F, {{400, 0, 400}, {1000, -7.0, 7.0}}});
    registry.add("hEtaResolution", "#eta resolution;#p_{T,part};Resolution", {HistType::kTH2F, {{400, 0, 400}, {1000, -1.0, 1.0}}});
  }

  // numeric values per species (see DCategory above) so a bare category
  // value is unambiguous even though it's stored in one shared column.

  template <typename U>
  int classifyDCandidate(int8_t dDecayChannel, int matchedFrom, int selectedAs)
  {
    int category = DCategory::Undefined;
    if constexpr (std::is_same_v<U, aod::CandidatesD0MCD>) {
      if ((std::abs(dDecayChannel) == hf2Prong::DecayChannelMain::D0ToPiK) && (matchedFrom != 0) && (selectedAs == matchedFrom)) {
        category = DCategory::Signal;
      } else if ((std::abs(dDecayChannel) == hf2Prong::DecayChannelMain::D0ToPiK) && (matchedFrom != 0) && (selectedAs == -1 * matchedFrom)) {
        category = DCategory::Reflection;
      } else if (std::abs(dDecayChannel) == hf2Prong::DecayChannelMain::D0ToPiKPi0) {
        category = DCategory::CorrBkg1;
      } else if (std::abs(dDecayChannel) == hf2Prong::DecayChannelMain::D0ToPiPi) {
        category = DCategory::CorrBkg2;
      } else if (std::abs(dDecayChannel) == hf2Prong::DecayChannelMain::D0ToPiPiPi0) {
        category = DCategory::CorrBkg3;
      } else if (std::abs(dDecayChannel) == hf2Prong::DecayChannelMain::D0ToKK) {
        category = DCategory::CorrBkg4;
      }
    } else if constexpr (std::is_same_v<U, aod::CandidatesDstarMCD>) {
      if ((std::abs(dDecayChannel) == hfDstar::DecayChannelMain::DstarToPiKPi) && (matchedFrom != 0) && (selectedAs == matchedFrom)) {
        category = DCategory::DstarSignal;
      } else if ((std::abs(dDecayChannel) == hfDstar::DecayChannelMain::DstarToPiKPi) && (matchedFrom != 0) && (selectedAs == -1 * matchedFrom)) {
        category = DCategory::DstarReflection;
      } else if (std::abs(dDecayChannel) == hfDstar::DecayChannelMain::DstarToPiKPiPi0) {
        category = DCategory::DstarCorrBkg;
      }
    }
    return category;
  }

  template <typename T, typename U, typename V>
  void doAnalysisData(T const& collision, U const& dCandidates, V const& jets)
  {
    if (!applyCollisionSelections(collision)) {
      return;
    }
    tableCollision(collision.posZ());
    for (const auto& dCandidate : dCandidates) {
      if (dCandidate.pt() < dPtCutMin) { // once settled on a mlcut, then add the lower bound of the systematics as a cut here
        continue;
      }
      const auto scores = dCandidate.mlScores();
      fillDHistograms(dCandidate, scores);

      float deltaM = -1.0; // sentinel: not applicable for D0
      if constexpr (std::is_same_v<U, aod::CandidatesDstarData>) {
        deltaM = dCandidate.m() - dCandidate.invMassCharm();
      }

      tableD(tableCollision.lastIndex(),
             scores[2],
             scores[1],
             scores[0],
             dCandidate.m(),
             deltaM,
             dCandidate.pt(),
             dCandidate.eta(),
             dCandidate.phi(),
             dCandidate.y());

      for (const auto& jet : jets) {
        if (jet.pt() < jetPtCutMin) {
          continue;
        }
        float dPhi = RecoDecay::constrainAngle(jet.phi() - dCandidate.phi(), -o2::constants::math::PI);
        if (std::abs(dPhi - o2::constants::math::PI) > (o2::constants::math::PI / 2)) {
          continue;
        }
        fillJetHistograms(jet, dPhi);
        tableJet(tableCollision.lastIndex(),
                 tableD.lastIndex(),
                 jet.pt(),
                 jet.eta(),
                 jet.phi(),
                 dPhi);
      }
    }
  }

  void processD0Data(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                     aod::CandidatesD0Data const& d0Candidates,
                     soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets)
  {
    doAnalysisData(collision, d0Candidates, jets);
  }
  PROCESS_SWITCH(JetCorrelationDmeson, processD0Data, "charged particle level jet analysis for D0", true);

  void processDstarData(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                        aod::CandidatesDstarData const& dstarCandidates,
                        soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets)
  {
    doAnalysisData(collision, dstarCandidates, jets);
  }
  PROCESS_SWITCH(JetCorrelationDmeson, processDstarData, "charged particle level jet analysis for D*", true);

  template <typename T, typename U, typename V>
  void doAnalysisMcDetector(T const& collision, U const& dCandidates, V const& jets)
  {
    if (!applyCollisionSelections(collision)) {
      return;
    }
    tableCollision(collision.posZ());
    for (const auto& dCandidate : dCandidates) {
      if (dCandidate.pt() < dPtCutMin) {
        continue;
      }
      const auto scores = dCandidate.mlScores();
      fillDHistograms(dCandidate, scores);

      int8_t dDecayChannel = dCandidate.flagMcMatchRec();
      int matchedFrom = 0;
      int selectedAs = 0;
      if (dDecayChannel > 0) {
        matchedFrom = 1;
      } else if (dDecayChannel < 0) {
        matchedFrom = -1;
      }
      if ((dCandidate.candidateSelFlag() & BIT(0)) != 0) {
        selectedAs = 1;
      } else if ((dCandidate.candidateSelFlag() & BIT(1)) != 0) {
        selectedAs = -1;
      }
      int category = classifyDCandidate<U>(dDecayChannel, matchedFrom, selectedAs);

      float deltaM = -1.0; // sentinel: not applicable for D0
      if constexpr (std::is_same_v<U, aod::CandidatesDstarMCD>) {
        deltaM = dCandidate.m() - dCandidate.invMassCharm();
      }

      tableDMcDetector(tableCollision.lastIndex(),
                       scores[2],
                       scores[1],
                       scores[0],
                       dCandidate.m(),
                       deltaM,
                       dCandidate.pt(),
                       dCandidate.eta(),
                       dCandidate.phi(),
                       dCandidate.y(),
                       category);
      for (const auto& jet : jets) {
        if (jet.pt() < jetPtCutMin) {
          continue;
        }
        float dPhi = RecoDecay::constrainAngle(jet.phi() - dCandidate.phi(), -o2::constants::math::PI);
        if (std::abs(dPhi - o2::constants::math::PI) > (o2::constants::math::PI / 2)) {
          continue;
        }
        fillJetHistograms(jet, dPhi);
        tableJetMcDetector(tableCollision.lastIndex(),
                           tableDMcDetector.lastIndex(),
                           jet.pt(),
                           jet.eta(),
                           jet.phi(),
                           dPhi);
      }
    }
  }

  void processD0McDetector(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                           aod::CandidatesD0MCD const& d0Candidates,
                           soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents> const& jets)
  {
    doAnalysisMcDetector(collision, d0Candidates, jets);
  }
  PROCESS_SWITCH(JetCorrelationDmeson, processD0McDetector, "charged detector level jet analysis for D0", false);

  void processDstarMcDetector(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                              aod::CandidatesDstarMCD const& dstarCandidates,
                              soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents> const& jets)
  {
    doAnalysisMcDetector(collision, dstarCandidates, jets);
  }
  PROCESS_SWITCH(JetCorrelationDmeson, processDstarMcDetector, "charged detector level jet analysis for D*", false);

  // Shared MC-particle logic.
  template <typename T, typename U, typename V>
  void doAnalysisMcParticle(T const& collision, U const& dMcPCandidates, V const& jets)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents, applyRCTSelections)) {
      return;
    }
    tableMcCollision(collision.posZ());
    for (const auto& dMcPCandidate : dMcPCandidates) {
      if (dMcPCandidate.pt() < dMcPtCutMin) {
        continue;
      }
      tableDMcParticle(tableMcCollision.lastIndex(),
                       dMcPCandidate.originMcGen(),
                       dMcPCandidate.pt(),
                       dMcPCandidate.eta(),
                       dMcPCandidate.phi(),
                       dMcPCandidate.y(),
                       dMcPCandidate.flagMcMatchGen());

      for (const auto& jet : jets) {
        if (jet.pt() < jetMcPtCutMin) {
          continue;
        }
        float dPhi = RecoDecay::constrainAngle(jet.phi() - dMcPCandidate.phi(), -o2::constants::math::PI);
        if (std::abs(dPhi - o2::constants::math::PI) > (o2::constants::math::PI / 2)) {
          continue;
        }
        fillJetHistograms(jet, dPhi);
        tableJetMcParticle(tableMcCollision.lastIndex(),
                           tableDMcParticle.lastIndex(),
                           jet.pt(),
                           jet.eta(),
                           jet.phi(),
                           dPhi);
      }
    }
  }

  void processD0McParticle(aod::JetMcCollision const& collision,
                           aod::CandidatesD0MCP const& d0McPCandidates,
                           soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents> const& jets)
  {
    doAnalysisMcParticle(collision, d0McPCandidates, jets);
  }
  PROCESS_SWITCH(JetCorrelationDmeson, processD0McParticle, "charged MC Particle jets for D0", false);

  void processDstarMcParticle(aod::JetMcCollision const& collision,
                              aod::CandidatesDstarMCP const& dstarMcPCandidates,
                              soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents> const& jets)
  {
    doAnalysisMcParticle(collision, dstarMcPCandidates, jets);
  }
  PROCESS_SWITCH(JetCorrelationDmeson, processDstarMcParticle, "charged MC Particle jets for D*", false);

  template <typename T, typename U, typename V, typename W, typename X, typename Y>
  void doAnalysisMcMatched(T const& collision, U const& dCandidates, V const& tracks,
                           W const& particles, X const& McDJets, Y const&)
  {
    if (!applyCollisionSelections(collision)) {
      return;
    }
    tableMatchedCollision(collision.posZ());
    for (const auto& dCandidate : dCandidates) {
      if (dCandidate.pt() < dPtCutMin) {
        continue;
      }
      const auto& dParticle = jethfutilities::matchedHFParticle(dCandidate, tracks, particles);
      for (const auto& McDJet : McDJets) {
        if (McDJet.pt() < jetPtCutMin) {
          continue;
        }
        float dPhiD = RecoDecay::constrainAngle(McDJet.phi() - dCandidate.phi(), -o2::constants::math::PI);
        if (std::abs(dPhiD - o2::constants::math::PI) > (o2::constants::math::PI / 2)) {
          continue;
        }
        if (McDJet.has_matchedJetGeo()) {
          for (auto const& McPJet : McDJet.template matchedJetGeo_as<Y>()) {
            float dPhiP = RecoDecay::constrainAngle(McPJet.phi() - dParticle.phi(), -o2::constants::math::PI);
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

  void processD0McMatched(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                          aod::CandidatesD0MCD const& d0Candidates,
                          aod::JetTracksMCD const& tracks,
                          aod::JetParticles const& particles,
                          soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets> const& McDJets,
                          aod::ChargedMCParticleLevelJets const& McPJets)
  {
    doAnalysisMcMatched(collision, d0Candidates, tracks, particles, McDJets, McPJets);
  }
  PROCESS_SWITCH(JetCorrelationDmeson, processD0McMatched, "process matching of particle level jets to detector level jets for D0", false);

  void processDstarMcMatched(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                             aod::CandidatesDstarMCD const& dstarCandidates,
                             aod::JetTracksMCD const& tracks,
                             aod::JetParticles const& particles,
                             soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets> const& McDJets,
                             aod::ChargedMCParticleLevelJets const& McPJets)
  {
    doAnalysisMcMatched(collision, dstarCandidates, tracks, particles, McDJets, McPJets);
  }
  PROCESS_SWITCH(JetCorrelationDmeson, processDstarMcMatched, "process matching of particle level jets to detector level jets for D*", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetCorrelationDmeson>(cfgc)};
}
