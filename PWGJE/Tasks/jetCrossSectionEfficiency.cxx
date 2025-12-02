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
  Configurable<double> jetPtMax{"jetPtMax", 200.0, "set jet pT bin max"};
  Configurable<float> jetEtaMin{"jetEtaMin", -0.7f, "minimum jet pseudorapidity"};
  Configurable<float> jetEtaMax{"jetEtaMax", 0.7f, "maximum jet pseudorapidity"};

  Configurable<float> jetAreaFractionMin{"jetAreaFractionMin", -99.0f, "used to make a cut on the jet areas"};
  Configurable<float> leadingConstituentPtMin{"leadingConstituentPtMin", -99.0f, "minimum pT selection on jet constituent"};
  Configurable<float> leadingConstituentPtMax{"leadingConstituentPtMax", 9999.0f, "maximum pT selection on jet constituent"};

  Configurable<int> trackOccupancyInTimeRangeMax{"trackOccupancyInTimeRangeMax", 999999, "maximum track occupancy of tracks in neighbouring collisions in a given time range; only applied to reconstructed collisions"};
  Configurable<int> trackOccupancyInTimeRangeMin{"trackOccupancyInTimeRangeMin", -999999, "minimum track occupancy of tracks in neighbouring collisions in a given time range; only applied to reconstructed collisions"};
  Configurable<int> acceptSplitCollisions{"acceptSplitCollisions", 0, "0: only look at mcCollisions that are not split; 1: accept split mcCollisions, 2: accept split mcCollisions but only look at the first reco collision associated with it"};

  enum AcceptSplitCollisionsOptions {
    NonSplitOnly = 0,
    SplitOkCheckAnyAssocColl,      // 1
    SplitOkCheckFirstAssocCollOnly // 2
  };

  float configSwitchLow = -98.0f;
  float configSwitchHigh = 9998.0f;

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
  int binOccupancy = -1;

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
    int nextBin = 4;
    std::set<int> addedBits;

    auto addStageIfNotAlreadyAdded = [&](int bit) {
      if (addedBits.find(bit) != addedBits.end()) {
        return;
      }
      addedBits.insert(bit);
      activeStages.push_back({bit, nextBin++});
    };

    // sel8, sel7, selKINT7 are composite selections that require multiple sub-bits.
    // initialiseEventSelectionBits returns them as single bits, but for step-by-step QA
    // we need to expand them into their constituent bits (kTVX, kTFBorder, kITSROFBorder).
    // Note: selMC, selMCFull, etc. are already expanded by initialiseEventSelectionBits.
    std::vector<std::pair<int, std::vector<int>>> compositeSelections = {
      {jetderiveddatautilities::JCollisionSel::sel8, {jetderiveddatautilities::JCollisionSel::selTVX, jetderiveddatautilities::JCollisionSel::selNoTimeFrameBorder, jetderiveddatautilities::JCollisionSel::selNoITSROFrameBorder}},
      {jetderiveddatautilities::JCollisionSel::sel7, {jetderiveddatautilities::JCollisionSel::selTVX, jetderiveddatautilities::JCollisionSel::selNoTimeFrameBorder, jetderiveddatautilities::JCollisionSel::selNoITSROFrameBorder}},
      {jetderiveddatautilities::JCollisionSel::selKINT7, {jetderiveddatautilities::JCollisionSel::selTVX, jetderiveddatautilities::JCollisionSel::selNoTimeFrameBorder, jetderiveddatautilities::JCollisionSel::selNoITSROFrameBorder}}
    };

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
      jetderiveddatautilities::JCollisionSel::selNoCollInRofStandard
    };

    for (auto bit : allPossibleBits) {
      if (hasBit(bit)) {
        addStageIfNotAlreadyAdded(bit);
      }
    }

    binZvtx = nextBin++;
    binCentrality = nextBin++;
    binOccupancy = nextBin++;

    AxisSpec jetPtAxis = {200, 0., jetPtMax, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec eventSelectionAxis = {nextBin - 1, 0.5, static_cast<double>(nextBin) - 0.5, "event selection"};

    registry.add("h2_jet_pt_part_eventselection",
                 "part jet pT vs event selection;#it{p}_{T,jet}^{part} (GeV/#it{c});event selection;counts",
                 {HistType::kTH2F, {jetPtAxis, eventSelectionAxis}});
    auto h2 = registry.get<TH2>(HIST("h2_jet_pt_part_eventselection"));
    h2->GetYaxis()->SetBinLabel(binInel, "INEL");
    h2->GetYaxis()->SetBinLabel(binNoReco, "noRecoColl");
    h2->GetYaxis()->SetBinLabel(binSplit, "splitColl");
    for (auto const& stage : activeStages) {
      const char* label = getStageLabel(stage.bit);
      if (label) {
        h2->GetYaxis()->SetBinLabel(stage.bin, label);
      }
    }
    h2->GetYaxis()->SetBinLabel(binZvtx, "zvtx");
    h2->GetYaxis()->SetBinLabel(binCentrality, "centralitycut");
    h2->GetYaxis()->SetBinLabel(binOccupancy, "occupancycut");

    registry.add("h_mccollisions_eventselection",
                 "number of mc events vs event selection;event selection;entries",
                 {HistType::kTH1F, {{nextBin - 1, 0.5, static_cast<double>(nextBin) - 0.5}}});
    auto h1 = registry.get<TH1>(HIST("h_mccollisions_eventselection"));
    h1->GetXaxis()->SetBinLabel(binInel, "INEL");
    h1->GetXaxis()->SetBinLabel(binNoReco, "noRecoColl");
    h1->GetXaxis()->SetBinLabel(binSplit, "splitColl");
    for (auto const& stage : activeStages) {
      const char* label = getStageLabel(stage.bit);
      if (label) {
        h1->GetXaxis()->SetBinLabel(stage.bin, label);
      }
    }
    h1->GetXaxis()->SetBinLabel(binZvtx, "zvtx");
    h1->GetXaxis()->SetBinLabel(binCentrality, "centralitycut");
    h1->GetXaxis()->SetBinLabel(binOccupancy, "occupancycut");

    registry.add("h_mccollisions_eventselection_weighted",
                 "number of weighted mc events vs event selection;event selection;weighted entries",
                 {HistType::kTH1F, {{nextBin - 1, 0.5, static_cast<double>(nextBin) - 0.5}}});
    auto h1w = registry.get<TH1>(HIST("h_mccollisions_eventselection_weighted"));
    h1w->GetXaxis()->SetBinLabel(binInel, "INEL");
    h1w->GetXaxis()->SetBinLabel(binNoReco, "noRecoColl");
    h1w->GetXaxis()->SetBinLabel(binSplit, "splitColl");
    for (auto const& stage : activeStages) {
      const char* label = getStageLabel(stage.bit);
      if (label) {
        h1w->GetXaxis()->SetBinLabel(stage.bin, label);
      }
    }
    h1w->GetXaxis()->SetBinLabel(binZvtx, "zvtx");
    h1w->GetXaxis()->SetBinLabel(binCentrality, "centralitycut");
    h1w->GetXaxis()->SetBinLabel(binOccupancy, "occupancycut");
  }

  template <typename TTracks, typename TJets>
  bool isAcceptedJet(TJets const& jet)
  {
    if (jetAreaFractionMin > configSwitchLow) {
      if (jet.area() < jetAreaFractionMin * o2::constants::math::PI * (jet.r() / 100.0) * (jet.r() / 100.0)) {
        return false;
      }
    }
    bool checkConstituentMinPt = (leadingConstituentPtMin > configSwitchLow);
    bool checkConstituentMaxPt = (leadingConstituentPtMax < configSwitchHigh);
    bool checkConstituentPt = checkConstituentMinPt || checkConstituentMaxPt;

    if (checkConstituentPt) {
      bool isMinLeadingConstituent = !checkConstituentMinPt;
      bool isMaxLeadingConstituent = true;

      for (const auto& constituent : jet.template tracks_as<TTracks>()) {
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

  void processCrossSectionEfficiency(aod::JetMcCollisions::iterator const& mccollision,
                                     soa::SmallGroups<aod::JetCollisionsMCD> const& collisions,
                                     soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents> const& jets,
                                     aod::JetParticles const&)
  {
    bool hasRecoColl = (collisions.size() >= 1);
    bool passesZvtxCut = (std::abs(mccollision.posZ()) <= vertexZCut);
    bool passesSplitCollCut = !(acceptSplitCollisions == NonSplitOnly && collisions.size() > 1);

    std::vector<bool> stagePassed(activeStages.size(), false);
    bool hasCustomEventSel = false;
    bool centralityIsGood = false;
    bool occupancyIsGood = false;
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
        if ((trackOccupancyInTimeRangeMin < collision.trackOccupancyInTimeRange()) && (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMax)) {
          occupancyIsGood = true;
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
        if ((trackOccupancyInTimeRangeMin < collision.trackOccupancyInTimeRange()) && (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMax)) {
          occupancyIsGood = true;
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

      if (!occupancyIsGood) {
        goto endEventCounter;
      }
      registry.fill(HIST("h_mccollisions_eventselection"), static_cast<double>(binOccupancy));
    }

  endEventCounter:

    for (auto const& jet : jets) {
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

      if (!occupancyIsGood) {
        goto nextJetUnweighted;
      }
      registry.fill(HIST("h2_jet_pt_part_eventselection"), jet.pt(), static_cast<double>(binOccupancy));

    nextJetUnweighted:
      ;
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

    std::vector<bool> stagePassed(activeStages.size(), false);
    bool hasCustomEventSel = false;
    bool centralityIsGood = false;
    bool occupancyIsGood = false;
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
        if ((trackOccupancyInTimeRangeMin < collision.trackOccupancyInTimeRange()) && (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMax)) {
          occupancyIsGood = true;
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
        if ((trackOccupancyInTimeRangeMin < collision.trackOccupancyInTimeRange()) && (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMax)) {
          occupancyIsGood = true;
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

      if (!occupancyIsGood) {
        goto endEventCounterWeighted;
      }
      registry.fill(HIST("h_mccollisions_eventselection_weighted"), static_cast<double>(binOccupancy), eventWeight);
    }

  endEventCounterWeighted:

    for (auto const& jet : jets) {
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

      if (!occupancyIsGood) {
        goto nextJetWeighted;
      }
      registry.fill(HIST("h2_jet_pt_part_eventselection"), jet.pt(), static_cast<double>(binOccupancy), eventWeight);

    nextJetWeighted:
      ;
    }
  }
  PROCESS_SWITCH(JetCrossSectionEfficiency, processCrossSectionEfficiencyWeighted, "jet spectra QC for MC particle level with step-by-step cuts (weighted)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<JetCrossSectionEfficiency>(cfgc)};
}