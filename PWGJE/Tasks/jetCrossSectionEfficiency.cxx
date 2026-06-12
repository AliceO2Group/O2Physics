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
// Step-by-step event selection efficiency for MC particle-level jets.
// Two ordering variants exposed as separate process functions:
//   - CollRecoFirst: reco collision required first, BC bits read from the reco-coll EvSel bitmask.
//   - BcBitsFirst:   BC bits read from the MC truth BC, then truth-side SBP, then reco at the end.
//
/// \author Joonsuk Bae <joonsuk.bae@cern.ch>

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/DataModel/Jet.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/RCTSelectionFlags.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <cmath>
#include <cstdint>
#include <string>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetCrossSectionEfficiency {

  HistogramRegistry registry;

  Configurable<std::string> eventSelections{"eventSelections", "sel8", "selTVX | selMC | selMCFull | sel8 | sel8Full"};
  Configurable<bool> skipMBGapEvents{"skipMBGapEvents", false, "reject min-bias gap events from hybrid MB+JJ MC productions"};
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> centralityMin{"centralityMin", -999.0f, "minimum centrality"};
  Configurable<float> centralityMax{"centralityMax", 999.0f, "maximum centrality"};
  Configurable<bool> checkCentFT0M{"checkCentFT0M", false, "0: centFT0C as default, 1: use centFT0M estimator"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9f, "minimum eta acceptance for tracks"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9f, "maximum eta acceptance for tracks"};
  Configurable<float> selectedJetsRadius{"selectedJetsRadius", 0.4f, "resolution parameter for histograms without radius"};
  Configurable<float> jetPtMin{"jetPtMin", 0.0f, "minimum jet pT"};
  Configurable<double> jetPtMax{"jetPtMax", 200.0, "set jet pT bin max"};
  Configurable<float> jetEtaMin{"jetEtaMin", -0.5f, "minimum jet pseudorapidity"};
  Configurable<float> jetEtaMax{"jetEtaMax", 0.5f, "maximum jet pseudorapidity"};

  Configurable<float> jetAreaFractionMin{"jetAreaFractionMin", -99.0f, "used to make a cut on the jet areas"};
  Configurable<float> leadingConstituentPtMinMCP{"leadingConstituentPtMinMCP", -99.0f, "minimum pT selection on MCP jet constituent"};
  Configurable<float> leadingConstituentPtMaxMCP{"leadingConstituentPtMaxMCP", 9999.0f, "maximum pT selection on MCP jet constituent"};

  Configurable<int> trackOccupancyInTimeRangeMax{"trackOccupancyInTimeRangeMax", 999999, "maximum track occupancy of tracks in neighbouring collisions in a given time range"};
  Configurable<int> trackOccupancyInTimeRangeMin{"trackOccupancyInTimeRangeMin", -999999, "minimum track occupancy of tracks in neighbouring collisions in a given time range"};
  Configurable<int> acceptSplitCollisions{"acceptSplitCollisions", 0, "0: only look at mcCollisions that are not split; 1: accept split mcCollisions, 2: accept split mcCollisions but only look at the first reco collision associated with it"};
  Configurable<float> pTHatMaxMCP{"pTHatMaxMCP", 999.0f, "maximum jet pT in units of pTHat for outlier rejection at MC particle level"};
  Configurable<float> pTHatAbsoluteMin{"pTHatAbsoluteMin", -99.0f, "minimum pTHat (drops events with pTHat below this)"};
  Configurable<float> pTHatExponent{"pTHatExponent", 6.0f, "exponent in the back-calculation of pTHat from the event weight"};
  Configurable<float> simPtRef{"simPtRef", 10.0f, "reference pT for the back-calculation of pTHat from the event weight"};
  Configurable<float> ptHardCalcMethodSwitch{"ptHardCalcMethodSwitch", 999.0f, "sentinel and threshold for the stored-ptHard branch (set <=0 to force weight-derived)"};
  Configurable<bool> applyRCT{"applyRCT", true, "apply RCT_pass step in the cascade (false: force the RCT step to pass)"};
  Configurable<std::string> rctSelectionsLabel{"rctSelectionsLabel", "CBT_hadronPID", "RCT selection preset name (see RCTSelectionFlags)"};

  o2::aod::rctsel::RCTFlagsChecker rctChecker;
  uint64_t rctMask = 0;

  bool applyTFB = true;
  bool applyROFB = true;
  bool applySBP = true;

  enum AcceptSplitCollisionsOptions {
    NonSplitOnly = 0,
    SplitOkCheckAnyAssocColl,      // 1
    SplitOkCheckFirstAssocCollOnly // 2
  };

  static constexpr float configSwitchLow = -98.0f;
  static constexpr float configSwitchHigh = 9998.0f;
  static constexpr float kBrokenPtHardSentinel = 1.0f;

  // CollRecoFirst: reco collision required first; BC bits read from the reco-coll EvSels.
  static constexpr int kBinCRF_Inel = 1;
  static constexpr int kBinCRF_RctPass = 2;
  static constexpr int kBinCRF_HasColl = 3;
  static constexpr int kBinCRF_Zreco = 4;
  static constexpr int kBinCRF_NoSplit = 5;
  static constexpr int kBinCRF_TVX = 6;
  static constexpr int kBinCRF_TFB = 7;
  static constexpr int kBinCRF_ROFB = 8;
  static constexpr int kBinCRF_SBP = 9;
  static constexpr int kBinCRF_N = 9;

  // BcBitsFirst: BC bits read from the MC truth BC; SBP from a Preslice count
  // (exactly one MC collision per truth BC) so it works before requiring reco.
  static constexpr int kBinBBF_Inel = 1;
  static constexpr int kBinBBF_RctPass = 2;
  static constexpr int kBinBBF_TVX = 3;
  static constexpr int kBinBBF_TFB = 4;
  static constexpr int kBinBBF_ROFB = 5;
  static constexpr int kBinBBF_TruthSBP = 6;
  static constexpr int kBinBBF_HasColl = 7;
  static constexpr int kBinBBF_Zreco = 8;
  static constexpr int kBinBBF_NoSplit = 9;
  static constexpr int kBinBBF_N = 9;

  Preslice<aod::JetMcCollisions> mcCollsPerBC = aod::jmccollision::bcId;

  void init(InitContext&)
  {
    if (!(acceptSplitCollisions == NonSplitOnly || acceptSplitCollisions == SplitOkCheckAnyAssocColl || acceptSplitCollisions == SplitOkCheckFirstAssocCollOnly)) {
      LOGF(fatal, "Configurable acceptSplitCollisions has wrong input value; stopping workflow");
    }

    rctChecker.init(static_cast<std::string>(rctSelectionsLabel));
    rctMask = rctChecker.value();

    const std::string es = eventSelections;
    if (es == "selTVX") {
      applyTFB = false;
      applyROFB = false;
      applySBP = false;
    } else if (es == "selMC") {
      applyTFB = true;
      applyROFB = false;
      applySBP = false;
    } else if (es == "selMCFull") {
      applyTFB = true;
      applyROFB = false;
      applySBP = true;
    } else if (es == "sel8") {
      applyTFB = true;
      applyROFB = true;
      applySBP = false;
    } else if (es == "sel8Full") {
      applyTFB = true;
      applyROFB = true;
      applySBP = true;
    } else {
      LOGF(fatal, "Configurable eventSelections=%s not supported; use selTVX, selMC, selMCFull, sel8, or sel8Full", es.c_str());
    }

    AxisSpec jetPtAxis = {200, 0., jetPtMax, "#it{p}_{T} (GeV/#it{c})"};

    if (doprocessCrossSectionEfficiency) {
      AxisSpec axCRF = {kBinCRF_N, 0.5, static_cast<double>(kBinCRF_N) + 0.5, "event selection (CollRecoFirst)"};
      registry.add("h2_jet_pt_part_eventselection_collRecoFirst",
                   "part jet pT vs event selection (CollRecoFirst);#it{p}_{T,jet}^{part} (GeV/#it{c});event selection;counts",
                   {HistType::kTH2F, {jetPtAxis, axCRF}});
      auto hCRF2 = registry.get<TH2>(HIST("h2_jet_pt_part_eventselection_collRecoFirst"));
      hCRF2->GetYaxis()->SetBinLabel(kBinCRF_Inel, "INEL");
      hCRF2->GetYaxis()->SetBinLabel(kBinCRF_RctPass, "+RCT_pass");
      hCRF2->GetYaxis()->SetBinLabel(kBinCRF_HasColl, "+hasRecoColl");
      hCRF2->GetYaxis()->SetBinLabel(kBinCRF_Zreco, "+|zReco|<10");
      hCRF2->GetYaxis()->SetBinLabel(kBinCRF_NoSplit, "+noSplit");
      hCRF2->GetYaxis()->SetBinLabel(kBinCRF_TVX, "+kTVX");
      hCRF2->GetYaxis()->SetBinLabel(kBinCRF_TFB, "+kNoTFB");
      hCRF2->GetYaxis()->SetBinLabel(kBinCRF_ROFB, "+kNoITSROFB");
      hCRF2->GetYaxis()->SetBinLabel(kBinCRF_SBP, "+kNoSBP");

      registry.add("h_mccollisions_eventselection_collRecoFirst",
                   "number of mc events vs event selection (CollRecoFirst);event selection;entries",
                   {HistType::kTH1F, {{kBinCRF_N, 0.5, static_cast<double>(kBinCRF_N) + 0.5}}});
      auto hCRF1 = registry.get<TH1>(HIST("h_mccollisions_eventselection_collRecoFirst"));
      hCRF1->GetXaxis()->SetBinLabel(kBinCRF_Inel, "INEL");
      hCRF1->GetXaxis()->SetBinLabel(kBinCRF_RctPass, "+RCT_pass");
      hCRF1->GetXaxis()->SetBinLabel(kBinCRF_HasColl, "+hasRecoColl");
      hCRF1->GetXaxis()->SetBinLabel(kBinCRF_Zreco, "+|zReco|<10");
      hCRF1->GetXaxis()->SetBinLabel(kBinCRF_NoSplit, "+noSplit");
      hCRF1->GetXaxis()->SetBinLabel(kBinCRF_TVX, "+kTVX");
      hCRF1->GetXaxis()->SetBinLabel(kBinCRF_TFB, "+kNoTFB");
      hCRF1->GetXaxis()->SetBinLabel(kBinCRF_ROFB, "+kNoITSROFB");
      hCRF1->GetXaxis()->SetBinLabel(kBinCRF_SBP, "+kNoSBP");
    }

    if (doprocessCrossSectionEfficiencyBcBitsFirst) {
      AxisSpec axBBF = {kBinBBF_N, 0.5, static_cast<double>(kBinBBF_N) + 0.5, "event selection (BcBitsFirst)"};
      registry.add("h2_jet_pt_part_eventselection_bcBitsFirst",
                   "part jet pT vs event selection (BcBitsFirst);#it{p}_{T,jet}^{part} (GeV/#it{c});event selection;counts",
                   {HistType::kTH2F, {jetPtAxis, axBBF}});
      auto hBBF2 = registry.get<TH2>(HIST("h2_jet_pt_part_eventselection_bcBitsFirst"));
      hBBF2->GetYaxis()->SetBinLabel(kBinBBF_Inel, "INEL");
      hBBF2->GetYaxis()->SetBinLabel(kBinBBF_RctPass, "+RCT_pass");
      hBBF2->GetYaxis()->SetBinLabel(kBinBBF_TVX, "+kTVX(truth)");
      hBBF2->GetYaxis()->SetBinLabel(kBinBBF_TFB, "+kNoTFB(truth)");
      hBBF2->GetYaxis()->SetBinLabel(kBinBBF_ROFB, "+kNoITSROFB(truth)");
      hBBF2->GetYaxis()->SetBinLabel(kBinBBF_TruthSBP, "+kNoSBP(truth)");
      hBBF2->GetYaxis()->SetBinLabel(kBinBBF_HasColl, "+hasColl");
      hBBF2->GetYaxis()->SetBinLabel(kBinBBF_Zreco, "+|zReco|<10");
      hBBF2->GetYaxis()->SetBinLabel(kBinBBF_NoSplit, "+noSplit");

      registry.add("h_mccollisions_eventselection_bcBitsFirst",
                   "number of mc events vs event selection (BcBitsFirst);event selection;entries",
                   {HistType::kTH1F, {{kBinBBF_N, 0.5, static_cast<double>(kBinBBF_N) + 0.5}}});
      auto hBBF1 = registry.get<TH1>(HIST("h_mccollisions_eventselection_bcBitsFirst"));
      hBBF1->GetXaxis()->SetBinLabel(kBinBBF_Inel, "INEL");
      hBBF1->GetXaxis()->SetBinLabel(kBinBBF_RctPass, "+RCT_pass");
      hBBF1->GetXaxis()->SetBinLabel(kBinBBF_TVX, "+kTVX(truth)");
      hBBF1->GetXaxis()->SetBinLabel(kBinBBF_TFB, "+kNoTFB(truth)");
      hBBF1->GetXaxis()->SetBinLabel(kBinBBF_ROFB, "+kNoITSROFB(truth)");
      hBBF1->GetXaxis()->SetBinLabel(kBinBBF_TruthSBP, "+kNoSBP(truth)");
      hBBF1->GetXaxis()->SetBinLabel(kBinBBF_HasColl, "+hasColl");
      hBBF1->GetXaxis()->SetBinLabel(kBinBBF_Zreco, "+|zReco|<10");
      hBBF1->GetXaxis()->SetBinLabel(kBinBBF_NoSplit, "+noSplit");
    }
  }

  template <typename TMcCollision>
  float computePtHat(TMcCollision const& mccollision)
  {
    float ptHardFromMc = ptHardCalcMethodSwitch;
    float storedPtHard = mccollision.ptHard();
    if (storedPtHard > kBrokenPtHardSentinel && storedPtHard < ptHardCalcMethodSwitch) {
      ptHardFromMc = storedPtHard;
    }
    float weight = mccollision.weight();
    return ptHardFromMc < ptHardCalcMethodSwitch
             ? ptHardFromMc
             : simPtRef / std::pow(weight, 1.0f / pTHatExponent);
  }

  template <typename TTracks, typename TJets>
  bool isAcceptedJet(TJets const& jet)
  {
    if (jetAreaFractionMin > configSwitchLow) {
      if (jet.area() < jetAreaFractionMin * o2::constants::math::PI * (jet.r() / 100.0) * (jet.r() / 100.0)) {
        return false;
      }
    }
    bool checkConstituentMinPt = (leadingConstituentPtMinMCP > configSwitchLow);
    bool checkConstituentMaxPt = (leadingConstituentPtMaxMCP < configSwitchHigh);
    bool checkConstituentPt = checkConstituentMinPt || checkConstituentMaxPt;

    if (checkConstituentPt) {
      bool isMinLeadingConstituent = !checkConstituentMinPt;
      bool isMaxLeadingConstituent = true;

      for (const auto& constituent : jet.template tracks_as<TTracks>()) {
        double constituentPt = constituent.pt();

        if (checkConstituentMinPt && constituentPt >= leadingConstituentPtMinMCP) {
          isMinLeadingConstituent = true;
        }
        if (checkConstituentMaxPt && constituentPt > leadingConstituentPtMaxMCP) {
          isMaxLeadingConstituent = false;
        }
      }
      return isMinLeadingConstituent && isMaxLeadingConstituent;
    }
    return true;
  }

  void processCrossSectionEfficiency(aod::JetMcCollisions::iterator const& mccollision,
                                     soa::SmallGroups<aod::JetCollisionsMCD> const& collisions,
                                     soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents> const& jets,
                                     aod::JetParticles const&,
                                     aod::JBCs const&)
  {
    if (skipMBGapEvents && mccollision.getSubGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      return;
    }

    bool hasRecoColl = (collisions.size() >= 1);
    bool passesZvtxCutReco = false;
    if (hasRecoColl) {
      if (acceptSplitCollisions == SplitOkCheckFirstAssocCollOnly) {
        auto const& col = collisions.begin();
        passesZvtxCutReco = (std::abs(col.posZ()) <= vertexZCut);
      } else {
        for (auto const& col : collisions) {
          if (std::abs(col.posZ()) <= vertexZCut) {
            passesZvtxCutReco = true;
            break;
          }
        }
      }
    }
    bool noSplitPass = (acceptSplitCollisions == NonSplitOnly) ? (collisions.size() == 1) : true;

    bool passesTVX = false, passesNoTFB = false, passesNoITSROFB = false, passesNoSBP = false;
    if (hasRecoColl) {
      auto const& col = collisions.begin();
      auto evSel = col.eventSel();
      passesTVX = (evSel & (1u << jetderiveddatautilities::JCollisionSel::selTVX)) != 0u;
      passesNoTFB = (evSel & (1u << jetderiveddatautilities::JCollisionSel::selNoTimeFrameBorder)) != 0u;
      passesNoITSROFB = (evSel & (1u << jetderiveddatautilities::JCollisionSel::selNoITSROFrameBorder)) != 0u;
      passesNoSBP = (evSel & (1u << jetderiveddatautilities::JCollisionSel::selNoSameBunchPileup)) != 0u;
    }

    bool passesRct = applyRCT ? (mccollision.bc_as<aod::JBCs>().rct_raw() & rctMask) == 0 : true;
    bool pass[kBinCRF_N + 1] = {false, true, passesRct, hasRecoColl, passesZvtxCutReco,
                                hasRecoColl && noSplitPass, passesTVX,
                                applyTFB ? passesNoTFB : true,
                                applyROFB ? passesNoITSROFB : true,
                                applySBP ? passesNoSBP : true};

    // Unified weight handling: MB MC -> weight=1 (no-op), JJ MC -> per-event sigma fraction.
    float weight = mccollision.weight();

    int sMax = 0;
    for (int s = kBinCRF_Inel; s <= kBinCRF_N; ++s) {
      if (!pass[s])
        break;
      registry.fill(HIST("h_mccollisions_eventselection_collRecoFirst"), static_cast<double>(s), weight);
      sMax = s;
    }
    if (sMax == 0) {
      return;
    }

    float pTHat = computePtHat(mccollision);
    if (pTHat < pTHatAbsoluteMin) {
      return;
    }

    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax) ||
          jet.pt() < jetPtMin || jet.pt() > pTHatMaxMCP * pTHat ||
          !isAcceptedJet<aod::JetParticles>(jet)) {
        continue;
      }
      for (int s = kBinCRF_Inel; s <= sMax; ++s) {
        registry.fill(HIST("h2_jet_pt_part_eventselection_collRecoFirst"), jet.pt(), static_cast<double>(s), weight);
      }
    }
  }
  PROCESS_SWITCH(JetCrossSectionEfficiency, processCrossSectionEfficiency,
                 "Cascade efficiency with collision-reco required first, BC bits read from the reco-coll EvSels bitmask", true);

  void processCrossSectionEfficiencyBcBitsFirst(aod::JetMcCollisions::iterator const& mccollision,
                                                soa::SmallGroups<aod::JetCollisionsMCD> const& collisions,
                                                soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents> const& jets,
                                                aod::JetParticles const&,
                                                aod::JBCs const&,
                                                aod::JetMcCollisions const& allMcCollisions)
  {
    if (skipMBGapEvents && mccollision.getSubGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      return;
    }

    auto truthBC = mccollision.bc_as<aod::JBCs>();
    bool passesTVXTruth = truthBC.selection_bit(aod::evsel::kIsTriggerTVX);
    bool passesNoTFBTruth = truthBC.selection_bit(aod::evsel::kNoTimeFrameBorder);
    bool passesNoITSROFBTruth = truthBC.selection_bit(aod::evsel::kNoITSROFrameBorder);

    auto sameBC = allMcCollisions.sliceBy(mcCollsPerBC, mccollision.bcId());
    bool truthNoSBP = (sameBC.size() == 1);

    bool hasRecoColl = (collisions.size() >= 1);
    bool passesZvtxCutReco = false;
    if (hasRecoColl) {
      if (acceptSplitCollisions == SplitOkCheckFirstAssocCollOnly) {
        auto const& col = collisions.begin();
        passesZvtxCutReco = (std::abs(col.posZ()) <= vertexZCut);
      } else {
        for (auto const& col : collisions) {
          if (std::abs(col.posZ()) <= vertexZCut) {
            passesZvtxCutReco = true;
            break;
          }
        }
      }
    }
    bool noSplitPass = (acceptSplitCollisions == NonSplitOnly) ? (collisions.size() == 1) : true;

    bool passesRct = applyRCT ? (truthBC.rct_raw() & rctMask) == 0 : true;
    bool pass[kBinBBF_N + 1] = {false, true, passesRct, passesTVXTruth,
                                applyTFB ? passesNoTFBTruth : true,
                                applyROFB ? passesNoITSROFBTruth : true,
                                applySBP ? truthNoSBP : true,
                                hasRecoColl, passesZvtxCutReco, hasRecoColl && noSplitPass};

    float weight = mccollision.weight();

    int sMax = 0;
    for (int s = kBinBBF_Inel; s <= kBinBBF_N; ++s) {
      if (!pass[s])
        break;
      registry.fill(HIST("h_mccollisions_eventselection_bcBitsFirst"), static_cast<double>(s), weight);
      sMax = s;
    }
    if (sMax == 0) {
      return;
    }

    float pTHat = computePtHat(mccollision);
    if (pTHat < pTHatAbsoluteMin) {
      return;
    }

    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax) ||
          jet.pt() < jetPtMin || jet.pt() > pTHatMaxMCP * pTHat ||
          !isAcceptedJet<aod::JetParticles>(jet)) {
        continue;
      }
      for (int s = kBinBBF_Inel; s <= sMax; ++s) {
        registry.fill(HIST("h2_jet_pt_part_eventselection_bcBitsFirst"), jet.pt(), static_cast<double>(s), weight);
      }
    }
  }
  PROCESS_SWITCH(JetCrossSectionEfficiency, processCrossSectionEfficiencyBcBitsFirst,
                 "Cascade efficiency with BC bits evaluated on MC truth BC before collision reco", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<JetCrossSectionEfficiency>(cfgc)};
}
