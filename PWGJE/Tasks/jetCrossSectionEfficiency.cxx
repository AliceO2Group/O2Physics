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
// Step-by-step event selection efficiency for MC particle-level jets
//
/// \author Joonsuk Bae <joonsuk.bae@cern.ch>

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/DataModel/Jet.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include <Framework/Configurable.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <algorithm>
#include <cmath>
#include <set>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetCrossSectionEfficiency {

  HistogramRegistry registry;

  Configurable<std::string> eventSelections{"eventSelections", "selTVX", "choose event selection (e.g., sel8, selMC, sel8Full, selTVX, etc.)"};
  Configurable<bool> skipMBGapEvents{"skipMBGapEvents", false, "flag to choose to reject min. bias gap events"};
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

  Configurable<int> trackOccupancyInTimeRangeMax{"trackOccupancyInTimeRangeMax", 999999, "maximum track occupancy of tracks in neighbouring collisions in a given time range; only applied to reconstructed collisions"};
  Configurable<int> trackOccupancyInTimeRangeMin{"trackOccupancyInTimeRangeMin", -999999, "minimum track occupancy of tracks in neighbouring collisions in a given time range; only applied to reconstructed collisions"};
  Configurable<int> acceptSplitCollisions{"acceptSplitCollisions", 0, "0: only look at mcCollisions that are not split; 1: accept split mcCollisions, 2: accept split mcCollisions but only look at the first reco collision associated with it"};
  Configurable<float> pTHatMaxMCP{"pTHatMaxMCP", 999.0f, "maximum fraction of hard scattering for jet acceptance in particle MC"};
  Configurable<float> pTHatExponent{"pTHatExponent", 6.0f, "exponent of the event weight for the calculation of pTHat"};
  Configurable<float> pTHatAbsoluteMin{"pTHatAbsoluteMin", -99.0f, "minimum value of pTHat"};

  enum AcceptSplitCollisionsOptions {
    NonSplitOnly = 0,
    SplitOkCheckAnyAssocColl,      // 1
    SplitOkCheckFirstAssocCollOnly // 2
  };

  static constexpr float configSwitchLow = -98.0f;
  static constexpr float configSwitchHigh = 9998.0f;

  std::vector<int> eventSelectionBits;

  struct StageInfo {
    int bit;
    int bin;
  };
  std::vector<StageInfo> activeStages;

  int binInel = 1;
  int binNoReco = 2;
  int binSplit = 3;
  int binZvtx = -1;
  int binCentrality = -1;

  const char* getStageLabel(int bit)
  {
    switch (bit) {
      case jetderiveddatautilities::JCollisionSel::selTVX:
        return "kTVX";
      case jetderiveddatautilities::JCollisionSel::selNoTimeFrameBorder:
        return "kTFBorder";
      case jetderiveddatautilities::JCollisionSel::selNoITSROFrameBorder:
        return "kITSROFBorder";
      case jetderiveddatautilities::JCollisionSel::selNoSameBunchPileup:
        return "NoSameBunchPileup";
      case jetderiveddatautilities::JCollisionSel::selIsGoodZvtxFT0vsPV:
        return "IsGoodZvtxFT0vsPV";
      case jetderiveddatautilities::JCollisionSel::selNoCollInTimeRangeStandard:
        return "NoCollInTimeRangeStd";
      case jetderiveddatautilities::JCollisionSel::selNoCollInRofStandard:
        return "NoCollInRofStd";
      default:
        return nullptr;
    }
  }

  void init(InitContext&)
  {
    if (!(acceptSplitCollisions == NonSplitOnly || acceptSplitCollisions == SplitOkCheckAnyAssocColl || acceptSplitCollisions == SplitOkCheckFirstAssocCollOnly)) {
      LOGF(fatal, "Configurable acceptSplitCollisions has wrong input value; stopping workflow");
    }

    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));

    auto hasBit = [&](int bit) {
      return std::find(eventSelectionBits.begin(), eventSelectionBits.end(), bit) != eventSelectionBits.end();
    };

    activeStages.clear();
    int nextEventSelectionBin = 4;
    std::set<int> addedBits;

    auto addStageIfNotAlreadyAdded = [&](int bit) {
      if (addedBits.find(bit) != addedBits.end()) {
        return;
      }
      addedBits.insert(bit);
      activeStages.push_back({bit, nextEventSelectionBin++});
    };

    // sel8, sel7, selKINT7 are composite selections that require multiple sub-bits.
    // initialiseEventSelectionBits returns them as single bits, but for step-by-step QA
    // we need to expand them into their constituent bits (kTVX, kTFBorder, kITSROFBorder).
    // Note: selMC, selMCFull, etc. are already expanded by initialiseEventSelectionBits.
    std::vector<std::pair<int, std::vector<int>>> compositeSelections = {
      {jetderiveddatautilities::JCollisionSel::sel8, {jetderiveddatautilities::JCollisionSel::selTVX, jetderiveddatautilities::JCollisionSel::selNoTimeFrameBorder, jetderiveddatautilities::JCollisionSel::selNoITSROFrameBorder}},
      {jetderiveddatautilities::JCollisionSel::sel7, {jetderiveddatautilities::JCollisionSel::selTVX, jetderiveddatautilities::JCollisionSel::selNoTimeFrameBorder, jetderiveddatautilities::JCollisionSel::selNoITSROFrameBorder}},
      {jetderiveddatautilities::JCollisionSel::selKINT7, {jetderiveddatautilities::JCollisionSel::selTVX, jetderiveddatautilities::JCollisionSel::selNoTimeFrameBorder, jetderiveddatautilities::JCollisionSel::selNoITSROFrameBorder}}};

    for (auto const& [compositeBit, subBits] : compositeSelections) {
      if (hasBit(compositeBit)) {
        for (auto subBit : subBits) {
          addStageIfNotAlreadyAdded(subBit);
        }
      }
    }

    std::vector<int> allPossibleBits = {
      jetderiveddatautilities::JCollisionSel::selTVX,
      jetderiveddatautilities::JCollisionSel::selNoTimeFrameBorder,
      jetderiveddatautilities::JCollisionSel::selNoITSROFrameBorder,
      jetderiveddatautilities::JCollisionSel::selNoSameBunchPileup,
      jetderiveddatautilities::JCollisionSel::selIsGoodZvtxFT0vsPV,
      jetderiveddatautilities::JCollisionSel::selNoCollInTimeRangeStandard,
      jetderiveddatautilities::JCollisionSel::selNoCollInRofStandard};

    for (auto bit : allPossibleBits) {
      if (hasBit(bit)) {
        addStageIfNotAlreadyAdded(bit);
      }
    }

    binZvtx = nextEventSelectionBin++;
    binCentrality = nextEventSelectionBin++;

    AxisSpec jetPtAxis = {200, 0., jetPtMax, "#it{p}_{T} (GeV/#it{c})"};
    int totalEventSelectionBins = nextEventSelectionBin - 1;
    AxisSpec eventSelectionAxis = {totalEventSelectionBins, 0.5, static_cast<double>(nextEventSelectionBin) - 0.5, "event selection"};

    if (doprocessCrossSectionEfficiency || doprocessCrossSectionEfficiencyWeighted) {
      registry.add("h2_jet_pt_part_eventselection",
                   "part jet pT vs event selection;#it{p}_{T,jet}^{part} (GeV/#it{c});event selection;counts",
                   {HistType::kTH2F, {jetPtAxis, eventSelectionAxis}});
      auto histJetPtVsEventSelection = registry.get<TH2>(HIST("h2_jet_pt_part_eventselection"));
      histJetPtVsEventSelection->GetYaxis()->SetBinLabel(binInel, "INEL");
      histJetPtVsEventSelection->GetYaxis()->SetBinLabel(binNoReco, "noRecoColl");
      histJetPtVsEventSelection->GetYaxis()->SetBinLabel(binSplit, "splitColl");
      for (auto const& stage : activeStages) {
        const char* label = getStageLabel(stage.bit);
        if (label) {
          histJetPtVsEventSelection->GetYaxis()->SetBinLabel(stage.bin, label);
        }
      }
      histJetPtVsEventSelection->GetYaxis()->SetBinLabel(binZvtx, "zvtx");
      histJetPtVsEventSelection->GetYaxis()->SetBinLabel(binCentrality, "centralitycut");
    }

    if (doprocessCrossSectionEfficiency) {
      registry.add("h_mccollisions_eventselection",
                   "number of mc events vs event selection;event selection;entries",
                   {HistType::kTH1F, {{totalEventSelectionBins, 0.5, static_cast<double>(nextEventSelectionBin) - 0.5}}});
      auto histMcCollisionsEventSelection = registry.get<TH1>(HIST("h_mccollisions_eventselection"));
      histMcCollisionsEventSelection->GetXaxis()->SetBinLabel(binInel, "INEL");
      histMcCollisionsEventSelection->GetXaxis()->SetBinLabel(binNoReco, "noRecoColl");
      histMcCollisionsEventSelection->GetXaxis()->SetBinLabel(binSplit, "splitColl");
      for (auto const& stage : activeStages) {
        const char* label = getStageLabel(stage.bit);
        if (label) {
          histMcCollisionsEventSelection->GetXaxis()->SetBinLabel(stage.bin, label);
        }
      }
      histMcCollisionsEventSelection->GetXaxis()->SetBinLabel(binZvtx, "zvtx");
      histMcCollisionsEventSelection->GetXaxis()->SetBinLabel(binCentrality, "centralitycut");
    }

    if (doprocessCrossSectionEfficiencyWeighted) {
      registry.add("h_mccollisions_eventselection_weighted",
                   "number of weighted mc events vs event selection;event selection;weighted entries",
                   {HistType::kTH1F, {{totalEventSelectionBins, 0.5, static_cast<double>(nextEventSelectionBin) - 0.5}}});
      auto histMcCollisionsEventSelectionWeighted = registry.get<TH1>(HIST("h_mccollisions_eventselection_weighted"));
      histMcCollisionsEventSelectionWeighted->GetXaxis()->SetBinLabel(binInel, "INEL");
      histMcCollisionsEventSelectionWeighted->GetXaxis()->SetBinLabel(binNoReco, "noRecoColl");
      histMcCollisionsEventSelectionWeighted->GetXaxis()->SetBinLabel(binSplit, "splitColl");
      for (auto const& stage : activeStages) {
        const char* label = getStageLabel(stage.bit);
        if (label) {
          histMcCollisionsEventSelectionWeighted->GetXaxis()->SetBinLabel(stage.bin, label);
        }
      }
      histMcCollisionsEventSelectionWeighted->GetXaxis()->SetBinLabel(binZvtx, "zvtx");
      histMcCollisionsEventSelectionWeighted->GetXaxis()->SetBinLabel(binCentrality, "centralitycut");
    }
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
                                     aod::JetParticles const&)
  {
    bool hasRecoColl = (collisions.size() >= 1);
    bool passesZvtxCut = (std::abs(mccollision.posZ()) <= vertexZCut);
    bool passesSplitCollCut = !(acceptSplitCollisions == NonSplitOnly && collisions.size() > 1);

    if (hasRecoColl) {
      bool occupancyIsGood = false;
      if (acceptSplitCollisions == SplitOkCheckFirstAssocCollOnly) {
        auto const& collision = collisions.begin();
        if ((trackOccupancyInTimeRangeMin < collision.trackOccupancyInTimeRange()) && (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMax)) {
          occupancyIsGood = true;
        }
      } else {
        for (auto const& collision : collisions) {
          if ((trackOccupancyInTimeRangeMin < collision.trackOccupancyInTimeRange()) && (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMax)) {
            occupancyIsGood = true;
            break;
          }
        }
      }
      if (!occupancyIsGood) {
        return;
      }
    }

    std::vector<bool> stagePassed(activeStages.size(), false);
    bool hasCustomEventSel = false;
    bool centralityIsGood = false;
    float centrality = mccollision.centFT0M();
    if (acceptSplitCollisions == SplitOkCheckFirstAssocCollOnly) {
      if (hasRecoColl) {
        auto const& collision = collisions.begin();
        for (std::size_t i = 0; i < activeStages.size(); ++i) {
          if (collision.eventSel() & (1 << activeStages[i].bit)) {
            stagePassed[i] = true;
          }
        }
        if (jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents)) {
          hasCustomEventSel = true;
        }
      }
      if ((centralityMin < centrality) && (centrality < centralityMax)) {
        centralityIsGood = true;
      }
    } else {
      for (auto const& collision : collisions) {
        for (std::size_t i = 0; i < activeStages.size(); ++i) {
          if (collision.eventSel() & (1 << activeStages[i].bit)) {
            stagePassed[i] = true;
          }
        }
        if (jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents)) {
          hasCustomEventSel = true;
        }
        float centrality = -1.0;
        checkCentFT0M ? centrality = collision.centFT0M() : centrality = collision.centFT0C();
        if ((centralityMin < centrality) && (centrality < centralityMax)) {
          centralityIsGood = true;
        }
      }
    }

    registry.fill(HIST("h_mccollisions_eventselection"), static_cast<double>(binInel));
    if (!hasRecoColl) {
    } else {
      registry.fill(HIST("h_mccollisions_eventselection"), static_cast<double>(binNoReco));
      if (!passesSplitCollCut) {
        goto endEventCounter;
      }
      registry.fill(HIST("h_mccollisions_eventselection"), static_cast<double>(binSplit));

      for (std::size_t i = 0; i < activeStages.size(); ++i) {
        if (!stagePassed[i]) {
          goto endEventCounter;
        }
        registry.fill(HIST("h_mccollisions_eventselection"), static_cast<double>(activeStages[i].bin));
      }

      if (!hasCustomEventSel || !passesZvtxCut) {
        goto endEventCounter;
      }
      registry.fill(HIST("h_mccollisions_eventselection"), static_cast<double>(binZvtx));

      if (!centralityIsGood) {
        goto endEventCounter;
      }
      registry.fill(HIST("h_mccollisions_eventselection"), static_cast<double>(binCentrality));
    }

  endEventCounter:

    for (auto const& jet : jets) {
      if (jet.r() != round(selectedJetsRadius * 100.0f)) {
        continue;
      }
      if (jet.pt() < jetPtMin) {
        continue;
      }
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetParticles>(jet)) {
        continue;
      }
      registry.fill(HIST("h2_jet_pt_part_eventselection"), jet.pt(), static_cast<double>(binInel));

      if (!hasRecoColl) {
        continue;
      }
      registry.fill(HIST("h2_jet_pt_part_eventselection"), jet.pt(), static_cast<double>(binNoReco));

      if (!passesSplitCollCut) {
        continue;
      }
      registry.fill(HIST("h2_jet_pt_part_eventselection"), jet.pt(), static_cast<double>(binSplit));

      for (std::size_t i = 0; i < activeStages.size(); ++i) {
        if (!stagePassed[i]) {
          goto nextJetUnweighted;
        }
        registry.fill(HIST("h2_jet_pt_part_eventselection"), jet.pt(), static_cast<double>(activeStages[i].bin));
      }

      if (!hasCustomEventSel || !passesZvtxCut) {
        goto nextJetUnweighted;
      }
      registry.fill(HIST("h2_jet_pt_part_eventselection"), jet.pt(), static_cast<double>(binZvtx));

      if (!centralityIsGood) {
        goto nextJetUnweighted;
      }
      registry.fill(HIST("h2_jet_pt_part_eventselection"), jet.pt(), static_cast<double>(binCentrality));

    nextJetUnweighted:;
    }
  }
  PROCESS_SWITCH(JetCrossSectionEfficiency, processCrossSectionEfficiency, "jet spectra QC for MC particle level with step-by-step cuts", false);

  void processCrossSectionEfficiencyWeighted(aod::JetMcCollisions::iterator const& mccollision,
                                             soa::SmallGroups<aod::JetCollisionsMCD> const& collisions,
                                             soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents> const& jets,
                                             aod::JetParticles const&)
  {
    bool hasRecoColl = (collisions.size() >= 1);
    bool passesZvtxCut = (std::abs(mccollision.posZ()) <= vertexZCut);
    bool passesSplitCollCut = !(acceptSplitCollisions == NonSplitOnly && collisions.size() > 1);

    if (hasRecoColl) {
      bool occupancyIsGood = false;
      if (acceptSplitCollisions == SplitOkCheckFirstAssocCollOnly) {
        auto const& collision = collisions.begin();
        if ((trackOccupancyInTimeRangeMin < collision.trackOccupancyInTimeRange()) && (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMax)) {
          occupancyIsGood = true;
        }
      } else {
        for (auto const& collision : collisions) {
          if ((trackOccupancyInTimeRangeMin < collision.trackOccupancyInTimeRange()) && (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMax)) {
            occupancyIsGood = true;
            break;
          }
        }
      }
      if (!occupancyIsGood) {
        return;
      }
    }

    std::vector<bool> stagePassed(activeStages.size(), false);
    bool hasCustomEventSel = false;
    bool centralityIsGood = false;
    float centrality = mccollision.centFT0M();
    if (acceptSplitCollisions == SplitOkCheckFirstAssocCollOnly) {
      if (hasRecoColl) {
        auto const& collision = collisions.begin();
        for (std::size_t i = 0; i < activeStages.size(); ++i) {
          if (collision.eventSel() & (1 << activeStages[i].bit)) {
            stagePassed[i] = true;
          }
        }
        if (jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents)) {
          hasCustomEventSel = true;
        }
      }
      if ((centralityMin < centrality) && (centrality < centralityMax)) {
        centralityIsGood = true;
      }
    } else {
      for (auto const& collision : collisions) {
        for (std::size_t i = 0; i < activeStages.size(); ++i) {
          if (collision.eventSel() & (1 << activeStages[i].bit)) {
            stagePassed[i] = true;
          }
        }
        if (jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents)) {
          hasCustomEventSel = true;
        }
        float centrality = -1.0;
        checkCentFT0M ? centrality = collision.centFT0M() : centrality = collision.centFT0C();
        if ((centralityMin < centrality) && (centrality < centralityMax)) {
          centralityIsGood = true;
        }
      }
    }

    float eventWeight = mccollision.weight();

    registry.fill(HIST("h_mccollisions_eventselection_weighted"), static_cast<double>(binInel), eventWeight);
    if (!hasRecoColl) {
    } else {
      registry.fill(HIST("h_mccollisions_eventselection_weighted"), static_cast<double>(binNoReco), eventWeight);
      if (!passesSplitCollCut) {
        goto endEventCounterWeighted;
      }
      registry.fill(HIST("h_mccollisions_eventselection_weighted"), static_cast<double>(binSplit), eventWeight);

      for (std::size_t i = 0; i < activeStages.size(); ++i) {
        if (!stagePassed[i]) {
          goto endEventCounterWeighted;
        }
        registry.fill(HIST("h_mccollisions_eventselection_weighted"), static_cast<double>(activeStages[i].bin), eventWeight);
      }

      if (!hasCustomEventSel || !passesZvtxCut) {
        goto endEventCounterWeighted;
      }
      registry.fill(HIST("h_mccollisions_eventselection_weighted"), static_cast<double>(binZvtx), eventWeight);

      if (!centralityIsGood) {
        goto endEventCounterWeighted;
      }
      registry.fill(HIST("h_mccollisions_eventselection_weighted"), static_cast<double>(binCentrality), eventWeight);
    }

  endEventCounterWeighted:

    float pTHat = 10.0f / (std::pow(eventWeight, 1.0f / pTHatExponent));

    for (auto const& jet : jets) {
      if (jet.r() != round(selectedJetsRadius * 100.0f)) {
        continue;
      }
      if (jet.pt() < jetPtMin) {
        continue;
      }
      if (jet.pt() > pTHatMaxMCP * pTHat || pTHat < pTHatAbsoluteMin) {
        continue;
      }
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetParticles>(jet)) {
        continue;
      }
      registry.fill(HIST("h2_jet_pt_part_eventselection"), jet.pt(), static_cast<double>(binInel), eventWeight);

      if (!hasRecoColl) {
        continue;
      }
      registry.fill(HIST("h2_jet_pt_part_eventselection"), jet.pt(), static_cast<double>(binNoReco), eventWeight);

      if (!passesSplitCollCut) {
        continue;
      }
      registry.fill(HIST("h2_jet_pt_part_eventselection"), jet.pt(), static_cast<double>(binSplit), eventWeight);

      for (std::size_t i = 0; i < activeStages.size(); ++i) {
        if (!stagePassed[i]) {
          goto nextJetWeighted;
        }
        registry.fill(HIST("h2_jet_pt_part_eventselection"), jet.pt(), static_cast<double>(activeStages[i].bin), eventWeight);
      }

      if (!hasCustomEventSel || !passesZvtxCut) {
        goto nextJetWeighted;
      }
      registry.fill(HIST("h2_jet_pt_part_eventselection"), jet.pt(), static_cast<double>(binZvtx), eventWeight);

      if (!centralityIsGood) {
        goto nextJetWeighted;
      }
      registry.fill(HIST("h2_jet_pt_part_eventselection"), jet.pt(), static_cast<double>(binCentrality), eventWeight);

    nextJetWeighted:;
    }
  }
  PROCESS_SWITCH(JetCrossSectionEfficiency, processCrossSectionEfficiencyWeighted, "jet spectra QC for MC particle level with step-by-step cuts (weighted)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<JetCrossSectionEfficiency>(cfgc)};
}
