// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

/// \file phiflowder.cxx
/// \brief Derived task for same-event and mixed-event Phi reconstruction
///
/// \author Prottay Das <prottay.das@cern.ch>

#include "PWGLF/DataModel/LFPhiFlowTables.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <Math/Vector4D.h>

#include <cstdint>
#include <map>
#include <string>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct phiflowder {

  // Event selection
  Configurable<float> centMin{"centMin", 0.f, "Minimum centrality"};
  Configurable<float> centMax{"centMax", 80.f, "Maximum centrality"};

  // PID selection from stored bitmask
  Configurable<bool> usePID{"usePID", false, "Enable PID selection using stored bitmask"};
  Configurable<int> pidChoice{"pidChoice", 2, "PID choice: 1, 2, 3, 4"};
  Configurable<bool> pidExclusive{"pidExclusive", false, "If true require only the chosen PID bit"};

  // Event mixing
  Configurable<int> nEvtMixing{"nEvtMixing", 5, "Number of events to mix"};
  ConfigurableAxis cfgVtxBins{"cfgVtxBins", {5, -10.0, 10.0}, "Mixing bins in z vertex"};
  ConfigurableAxis cfgMultBins{"cfgMultBins", {8, 0.0, 80.0}, "Mixing bins in centrality"};

  // Histogram axes
  ConfigurableAxis configAxisInvMass{"configAxisInvMass", {120, 0.9, 1.2}, "#it{M}_{K^{+}K^{-}} (GeV/#it{c}^{2})"};
  ConfigurableAxis configAxisPt{"configAxisPt", {100, 0.0, 10.0}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis configAxisCentrality{"configAxisCentrality", {8, 0.0, 80.0}, "Centrality (%)"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // PID bits: same convention as producer
  static constexpr uint8_t kPID1 = 1u << 0;
  static constexpr uint8_t kPID2 = 1u << 1;
  static constexpr uint8_t kPID3 = 1u << 2;
  static constexpr uint8_t kPID4 = 1u << 3;

  struct KaonDaughter {
    float px;
    float py;
    float pz;
    uint8_t pidMask;
  };

  void init(o2::framework::InitContext&)
  {

    histos.add("hCentrality", "Centrality distribution", kTH1F, {configAxisCentrality});

    histos.add("hSparseSame", "Same-event sparse", kTHnSparseF,
               {configAxisInvMass, configAxisPt, configAxisCentrality});
    histos.add("hSparseMixed", "Mixed-event sparse", kTHnSparseF,
               {configAxisInvMass, configAxisPt, configAxisCentrality});
  }

  uint8_t requiredPidBit() const
  {
    switch (pidChoice.value) {
      case 1:
        return kPID1;
      case 2:
        return kPID2;
      case 3:
        return kPID3;
      case 4:
        return kPID4;
      default:
        LOGF(warn, "pidChoice=%d is invalid. Using PID2.", pidChoice.value);
        return kPID2;
    }
  }

  bool passDaughterPID(uint8_t mask) const
  {
    if (!usePID.value) {
      return true;
    }

    const uint8_t bit = requiredPidBit();

    if (pidExclusive.value) {
      return mask == bit;
    }

    return (mask & bit) != 0;
  }

  template <typename TPhi>
  bool passPairPID(const TPhi& phiCand) const
  {
    if (!usePID.value) {
      return true;
    }

    const uint16_t pairMask = phiCand.kaonPidMask();

    const uint8_t kPlusMask = pairMask & 0xFF;
    const uint8_t kMinusMask = (pairMask >> 8) & 0xFF;

    return passDaughterPID(kPlusMask) && passDaughterPID(kMinusMask);
  }

  void collectUniqueKPlusForEvent(const aod::KaonTracks& kaontracks,
                                  int64_t eventId,
                                  std::map<int64_t, KaonDaughter>& uniqueKPlus)
  {
    uniqueKPlus.clear();

    for (const auto& phiCand : kaontracks) {
      if (phiCand.kaonkaoneventId() != eventId) {
        continue;
      }

      const uint16_t pairMask = phiCand.kaonPidMask();
      const uint8_t kPlusMask = pairMask & 0xFF;

      if (!passDaughterPID(kPlusMask)) {
        continue;
      }

      const int64_t kPlusId = phiCand.kaonIndex1();

      if (uniqueKPlus.find(kPlusId) == uniqueKPlus.end()) {
        uniqueKPlus[kPlusId] = {
          phiCand.d1Px(),
          phiCand.d1Py(),
          phiCand.d1Pz(),
          kPlusMask};
      }
    }
  }

  void collectUniqueKMinusForEvent(const aod::KaonTracks& kaontracks,
                                   int64_t eventId,
                                   std::map<int64_t, KaonDaughter>& uniqueKMinus)
  {
    uniqueKMinus.clear();

    for (const auto& phiCand : kaontracks) {
      if (phiCand.kaonkaoneventId() != eventId) {
        continue;
      }

      const uint16_t pairMask = phiCand.kaonPidMask();
      const uint8_t kMinusMask = (pairMask >> 8) & 0xFF;

      if (!passDaughterPID(kMinusMask)) {
        continue;
      }

      const int64_t kMinusId = phiCand.kaonIndex2();

      if (uniqueKMinus.find(kMinusId) == uniqueKMinus.end()) {
        uniqueKMinus[kMinusId] = {
          phiCand.d2Px(),
          phiCand.d2Py(),
          phiCand.d2Pz(),
          kMinusMask};
      }
    }
  }

  // Filter reduced events
  Filter centralityFilter = (aod::kaonkaonevent::cent >= centMin &&
                             aod::kaonkaonevent::cent < centMax);
  using EventCandidates = soa::Filtered<aod::KaonkaonEvents>;

  void processSameData(EventCandidates::iterator const& collision,
                       aod::KaonTracks const& kaontracks)
  {
    const float centrality = collision.cent();
    const int64_t eventId = collision.globalIndex();

    histos.fill(HIST("hCentrality"), centrality);

    for (const auto& phiCand : kaontracks) {
      if (phiCand.kaonkaoneventId() != eventId) {
        continue;
      }

      if (!passPairPID(phiCand)) {
        continue;
      }

      ROOT::Math::PxPyPzMVector kPlus(
        phiCand.d1Px(),
        phiCand.d1Py(),
        phiCand.d1Pz(),
        o2::constants::physics::MassKPlus);

      ROOT::Math::PxPyPzMVector kMinus(
        phiCand.d2Px(),
        phiCand.d2Py(),
        phiCand.d2Pz(),
        o2::constants::physics::MassKPlus);

      auto phi = kPlus + kMinus;

      const float phiMass = phi.M();
      const float phiPt = phi.Pt();

      histos.fill(HIST("hSparseSame"), phiMass, phiPt, centrality);
    }
  }

  PROCESS_SWITCH(phiflowder, processSameData, "Process same-event phi candidates", true);

  using BinningType = ColumnBinningPolicy<aod::kaonkaonevent::Posz,
                                          aod::kaonkaonevent::Cent>;

  void processMixedData(EventCandidates const& collisions,
                        aod::KaonTracks const& kaontracks)
  {

    BinningType colBinning{{cfgVtxBins, cfgMultBins}, false};

    std::map<int64_t, KaonDaughter> kPlusEvent1;
    std::map<int64_t, KaonDaughter> kMinusEvent2;
    // int nMixedEventPairs = 0;

    for (const auto& [collision1, collision2] :
         selfCombinations(colBinning, nEvtMixing.value, -1, collisions, collisions)) {

      if (collision1.globalIndex() == collision2.globalIndex()) {
        continue;
      }
      // nMixedEventPairs++;

      const float centrality = collision1.cent();

      collectUniqueKPlusForEvent(kaontracks,
                                 collision1.globalIndex(),
                                 kPlusEvent1);

      collectUniqueKMinusForEvent(kaontracks,
                                  collision2.globalIndex(),
                                  kMinusEvent2);

      // histos.fill(HIST("hNplusUnique"), kPlusEvent1.size());
      // histos.fill(HIST("hNminusUnique"), kMinusEvent2.size());

      // const int nMixedKPairs = static_cast<int>(kPlusEvent1.size() * kMinusEvent2.size());

      // histos.fill(HIST("hNmixedKPairsPerEventPair"), nMixedKPairs);

      for (const auto& [kPlusId, kPlus] : kPlusEvent1) {
        ROOT::Math::PxPyPzMVector kPlusVec(
          kPlus.px,
          kPlus.py,
          kPlus.pz,
          o2::constants::physics::MassKPlus);

        for (const auto& [kMinusId, kMinus] : kMinusEvent2) {
          ROOT::Math::PxPyPzMVector kMinusVec(
            kMinus.px,
            kMinus.py,
            kMinus.pz,
            o2::constants::physics::MassKPlus);

          auto phiMix = kPlusVec + kMinusVec;

          histos.fill(HIST("hSparseMixed"), phiMix.M(), phiMix.Pt(), centrality);
        }
      }
    }
    // histos.fill(HIST("hNMixEventPairs"), nMixedEventPairs);
  }
  PROCESS_SWITCH(phiflowder, processMixedData, "Process mixed-event phi candidates", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<phiflowder>(cfgc)};
}
