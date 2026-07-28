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

/// \file cksspinalignder.cxx
/// \brief Derived task for charged KStar spin alignment using reduced K0s and pion tables
///
/// \author prottay.das@cern.ch

#include "PWGLF/DataModel/LFCKSSpinalignmentTables.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <Math/GenVector/Boost.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct cksspinalignder {

  // Event selection
  Configurable<float> centMin{"centMin", 0.0f, "Minimum centrality"};
  Configurable<float> centMax{"centMax", 80.0f, "Maximum centrality"};

  // K0s selection from stored reduced table
  Configurable<bool> applyK0sSelection{"applyK0sSelection", true, "Apply additional K0s selection in derived task"};
  Configurable<float> cosPA{"cosPA", 0.995f, "Minimum V0 cosine pointing angle"};
  Configurable<float> radiusMin{"radiusMin", 1.2f, "Minimum V0 radius"};
  Configurable<float> radiusMax{"radiusMax", 100.0f, "Maximum V0 radius"};
  Configurable<float> dcaPion{"dcaPion", 0.1f, "Minimum V0 daughter DCA to PV"};
  Configurable<float> dcaDaughters{"dcaDaughters", 1.0f, "Maximum DCA between V0 daughters"};
  Configurable<float> k0sPtMin{"k0sPtMin", 0.0f, "Minimum K0s pT"};
  Configurable<float> k0sPtMax{"k0sPtMax", 10.0f, "Maximum K0s pT"};
  Configurable<float> k0sEtaMax{"k0sEtaMax", 0.8f, "Maximum K0s eta"};
  Configurable<float> k0sMassMin{"k0sMassMin", 0.45f, "Minimum K0s mass"};
  Configurable<float> k0sMassMax{"k0sMassMax", 0.55f, "Maximum K0s mass"};

  // Bachelor pion PID
  Configurable<bool> usePID{"usePID", true, "Apply bachelor pion PID using stored bitmask"};
  Configurable<int> pidChoice{"pidChoice", 4, "PID choice: 1, 2, 3, or 4"};
  Configurable<bool> pidExclusive{"pidExclusive", false, "If true require only selected PID bit. If false require selected bit among possible bits"};

  // Charge selection
  Configurable<int> chargeSelection{"chargeSelection", 0, "Bachelor pion charge selection: 0 = both, +1 = K*+, -1 = K*-"};

  // Event-plane choice
  Configurable<int> epChoice{"epChoice", 0, "Event plane: 0 = FT0C, 1 = FT0A, 2 = TPC"};

  // Event mixing
  Configurable<int> nEvtMixing{"nEvtMixing", 5, "Number of mixed events"};
  ConfigurableAxis cfgVtxBins{"cfgVtxBins", {10, -10.0f, 10.0f}, "Mixing bins in z vertex"};
  ConfigurableAxis cfgCentBins{"cfgCentBins", {8, 0.0f, 80.0f}, "Mixing bins in centrality"};

  // Sparse axes
  ConfigurableAxis configThnAxisInvMass{"configThnAxisInvMass", {150, 0.75f, 1.05f}, "#it{M}_{K^{0}_{S}#pi} (GeV/#it{c}^{2})"};
  ConfigurableAxis configThnAxisPt{"configThnAxisPt", {100, 0.0f, 10.0f}, "#it{p}_{T}^{K^{0}_{S}#pi} (GeV/#it{c})"};
  ConfigurableAxis configThnAxisSA{"configThnAxisSA", {20, -1.0f, 1.0f}, "cos#it{#theta}^{*}"};
  ConfigurableAxis configThnAxisCentrality{"configThnAxisCentrality", {8, 0.0f, 80.0f}, "Centrality (%)"};
  ConfigurableAxis configThnAxisCharge{"configThnAxisCharge", {2, -1.5f, 1.5f}, "Bachelor pion charge"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    histos.add("hCentrality", "Centrality distribution;Centrality (%);Events", kTH1F, {{80, 0.0f, 80.0f}});
    histos.add("hKShortMass", "K^{0}_{S} mass;M_{#pi^{+}#pi^{-}} (GeV/#it{c}^{2});Counts", kTH1F, {{200, 0.4f, 0.6f}});

    histos.add("hSparsesame",
               "Same-event K^{0}_{S}#pi;M_{K^{0}_{S}#pi};p_{T};cos#theta^{*};centrality;charge",
               kTHnSparseF,
               {configThnAxisInvMass, configThnAxisPt, configThnAxisSA, configThnAxisCentrality, configThnAxisCharge});

    histos.add("hSparsemix",
               "Mixed-event K^{0}_{S}#pi;M_{K^{0}_{S}#pi};p_{T};cos#theta^{*};centrality;charge",
               kTHnSparseF,
               {configThnAxisInvMass, configThnAxisPt, configThnAxisSA, configThnAxisCentrality, configThnAxisCharge});
  }

  // PID bit mapping from producer
  static constexpr uint8_t kPID1 = 1u << 0;
  static constexpr uint8_t kPID2 = 1u << 1;
  static constexpr uint8_t kPID3 = 1u << 2;
  static constexpr uint8_t kPID4 = 1u << 3;

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
        LOGF(warn, "pidChoice=%d is invalid. Using PID4.", pidChoice.value);
        return kPID4;
    }
  }

  template <typename TPion>
  bool passSelectedPID(const TPion& pion) const
  {
    if (!usePID.value) {
      return true;
    }

    const uint8_t mask = pion.pionPidMask();
    const uint8_t bit = requiredPidBit();

    if (pidExclusive.value) {
      return mask == bit;
    }

    return (mask & bit) != 0;
  }

  template <typename TPion>
  bool passChargeSelection(const TPion& pion) const
  {
    const int charge = static_cast<int>(pion.charge());

    if (charge == 0) {
      return false;
    }

    if (chargeSelection.value == 0) {
      return true;
    }

    return charge == chargeSelection.value;
  }

  template <typename TK0s>
  bool selectionV0(const TK0s& candidate) const
  {
    if (!applyK0sSelection.value) {
      return true;
    }

    ROOT::Math::PxPyPzMVector k0s(candidate.kShortPx(),
                                  candidate.kShortPy(),
                                  candidate.kShortPz(),
                                  candidate.kShortMass());

    const float k0sPt = k0s.Pt();
    const float k0sEta = k0s.Eta();

    if (std::abs(k0sEta) > k0sEtaMax.value) {
      return false;
    }

    if (candidate.kShortMass() < k0sMassMin.value ||
        candidate.kShortMass() > k0sMassMax.value) {
      return false;
    }

    if (candidate.v0Cospa() < cosPA.value) {
      return false;
    }

    if (candidate.v0Radius() < radiusMin.value ||
        candidate.v0Radius() > radiusMax.value) {
      return false;
    }

    if (candidate.dcaBetweenDaughter() > dcaDaughters.value) {
      return false;
    }

    if (candidate.dcaPositive() < dcaPion.value ||
        candidate.dcaNegative() < dcaPion.value) {
      return false;
    }

    if (k0sPt < k0sPtMin.value ||
        k0sPt > k0sPtMax.value) {
      return false;
    }

    return true;
  }

  template <typename TK0s, typename TCollision>
  bool k0sBelongsToEvent(const TK0s& k0s, const TCollision& collision) const
  {
    return k0s.kshortpioneventId() == collision.globalIndex();
  }

  template <typename TPion, typename TCollision>
  bool pionBelongsToEvent(const TPion& pion, const TCollision& collision) const
  {
    return pion.kshortpioneventId() == collision.globalIndex();
  }

  template <typename TCollision>
  float getEventPlane(const TCollision& collision) const
  {
    switch (epChoice.value) {
      case 0:
        return collision.psiFT0C();
      case 1:
        return collision.psiFT0A();
      case 2:
        return collision.psiTPC();
      default:
        LOGF(warn, "epChoice=%d is invalid. Using FT0C.", epChoice.value);
        return collision.psiFT0C();
    }
  }

  bool computeCosThetaStar(const ROOT::Math::PxPyPzMVector& k0s,
                           const ROOT::Math::PxPyPzMVector& kstar,
                           float psi,
                           float& cosThetaStar) const
  {
    ROOT::Math::Boost boost{kstar.BoostToCM()};
    ROOT::Math::PxPyPzMVector k0sCM = boost(k0s);
    ROOT::Math::XYZVector k0sCMVec = k0sCM.Vect();

    ROOT::Math::XYZVector eventPlaneVector(std::sin(2.0f * psi),
                                           -std::cos(2.0f * psi),
                                           0.0f);

    const double denom =
      std::sqrt(k0sCMVec.Mag2()) *
      std::sqrt(eventPlaneVector.Mag2());

    if (denom <= 0.0 || !std::isfinite(denom)) {
      return false;
    }

    cosThetaStar =
      static_cast<float>(eventPlaneVector.Dot(k0sCMVec) / denom);

    cosThetaStar = std::clamp(cosThetaStar, -1.0f, 1.0f);

    return true;
  }

  template <typename TK0s, typename TPion, typename TCollision>
  bool fillKstarPair(const TK0s& k0sCand,
                     const TPion& pionCand,
                     const TCollision& collision,
                     bool isMixed)
  {
    if (!selectionV0(k0sCand)) {
      return false;
    }

    if (!passSelectedPID(pionCand)) {
      return false;
    }

    if (!passChargeSelection(pionCand)) {
      return false;
    }

    // Same-event self-combination rejection:
    // bachelor pion must not be one of the K0s daughters.
    if (!isMixed) {
      if (pionCand.pionBachIndex() == k0sCand.pionIndex1() ||
          pionCand.pionBachIndex() == k0sCand.pionIndex2()) {
        return false;
      }
    }

    ROOT::Math::PxPyPzMVector k0s(k0sCand.kShortPx(),
                                  k0sCand.kShortPy(),
                                  k0sCand.kShortPz(),
                                  k0sCand.kShortMass());

    ROOT::Math::PxPyPzMVector pion(pionCand.pionBachPx(),
                                   pionCand.pionBachPy(),
                                   pionCand.pionBachPz(),
                                   o2::constants::physics::MassPionCharged);

    const auto kstar = k0s + pion;

    const float mass = kstar.M();
    const float pt = kstar.Pt();
    const float centrality = collision.cent();
    const float charge = static_cast<float>(pionCand.charge());
    const float psi = getEventPlane(collision);

    float cosThetaStar = 0.0f;

    if (!computeCosThetaStar(k0s, kstar, psi, cosThetaStar)) {
      return false;
    }

    if (isMixed) {
      histos.fill(HIST("hSparsemix"), mass, pt, cosThetaStar, centrality, charge);
    } else {
      histos.fill(HIST("hSparsesame"), mass, pt, cosThetaStar, centrality, charge);
    }

    return true;
  }

  Filter centralityFilter =
    (aod::kshortpionevent::cent >= centMin) &&
    (aod::kshortpionevent::cent < centMax);

  using EventCandidates = soa::Filtered<aod::KShortpionEvents>;

  void processSameData(EventCandidates::iterator const& collision,
                       aod::KShortTracks const& k0sTracks,
                       aod::PionTracks const& pionTracks)
  {
    const float centrality = collision.cent();

    histos.fill(HIST("hCentrality"), centrality);

    for (const auto& k0sCand : k0sTracks) {
      if (!k0sBelongsToEvent(k0sCand, collision)) {
        continue;
      }

      if (!selectionV0(k0sCand)) {
        continue;
      }

      histos.fill(HIST("hKShortMass"), k0sCand.kShortMass());

      for (const auto& pionCand : pionTracks) {
        if (!pionBelongsToEvent(pionCand, collision)) {
          continue;
        }

        fillKstarPair(k0sCand, pionCand, collision, false);
      }
    }
  }

  PROCESS_SWITCH(cksspinalignder, processSameData, "Process same-event charged KStar", true);

  SliceCache cache;
  using BinningType = ColumnBinningPolicy<aod::kshortpionevent::Posz, aod::kshortpionevent::Cent>;

  template <typename TK0sGroup, typename TPionGroup, typename TCollision>
  void fillMixedPairs(const TK0sGroup& k0sGroup,
                      const TPionGroup& pionGroup,
                      const TCollision& referenceCollision)
  {
    for (const auto& k0sCand : k0sGroup) {
      for (const auto& pionCand : pionGroup) {
        fillKstarPair(k0sCand, pionCand, referenceCollision, true);
      }
    }
  }

  void processMixedData(EventCandidates const& collisions,
                        aod::KShortTracks const& k0sTracks,
                        aod::PionTracks const& pionTracks)
  {
    BinningType colBinning{{cfgVtxBins, cfgCentBins}, true};

    for (const auto& [collision1, collision2] :
         selfCombinations(colBinning, nEvtMixing.value, -1, collisions, collisions)) {

      if (collision1.globalIndex() == collision2.globalIndex()) {
        continue;
      }

      // Direction 1: K0s from event 1 + pion from event 2.
      for (const auto& k0sCand : k0sTracks) {
        if (!k0sBelongsToEvent(k0sCand, collision1)) {
          continue;
        }

        for (const auto& pionCand : pionTracks) {
          if (!pionBelongsToEvent(pionCand, collision2)) {
            continue;
          }

          fillKstarPair(k0sCand, pionCand, collision1, true);
        }
      }

      // Direction 2: K0s from event 2 + pion from event 1.
      for (const auto& k0sCand : k0sTracks) {
        if (!k0sBelongsToEvent(k0sCand, collision2)) {
          continue;
        }

        for (const auto& pionCand : pionTracks) {
          if (!pionBelongsToEvent(pionCand, collision1)) {
            continue;
          }

          fillKstarPair(k0sCand, pionCand, collision2, true);
        }
      }
    }
  }

  PROCESS_SWITCH(cksspinalignder, processMixedData, "Process mixed-event charged KStar", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<cksspinalignder>(cfgc)};
}
