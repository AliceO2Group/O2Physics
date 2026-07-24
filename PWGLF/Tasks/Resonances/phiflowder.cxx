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

/// \file phiflowder.cxx
/// \brief Derived task for same-event and mixed-event Phi reconstruction for v1
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

  Configurable<float> centMin{"centMin", 0.f, "Minimum centrality"};
  Configurable<float> centMax{"centMax", 80.f, "Maximum centrality"};

  Configurable<bool> useStoredPID{"useStoredPID", false, "Apply PID-bit selection using stored kaon PID mask"};
  Configurable<int> pidChoice{"pidChoice", 2, "PID choice: 1, 2, 3, or 4"};
  Configurable<bool> pidExclusive{"pidExclusive", false, "Require exactly the selected PID bit"};

  Configurable<bool> applyDeepAngle{"applyDeepAngle", false, "Apply minimum opening-angle cut for K+K- pairs"};
  Configurable<float> deepAngleCut{"deepAngleCut", 0.04f, "Minimum K+K- opening angle"};

  Configurable<float> massMin{"massMin", 0.9f, "Minimum K+K- mass to fill"};
  Configurable<float> massMax{"massMax", 1.2f, "Maximum K+K- mass to fill"};

  ConfigurableAxis axisInvMass{"axisInvMass", {120, 0.9f, 1.2f}, "#it{M}_{K^{+}K^{-}} (GeV/#it{c}^{2})"};
  ConfigurableAxis axisPhiPt{"axisPhiPt", {100, 0.f, 10.f}, "#it{p}_{T}^{K^{+}K^{-}} (GeV/#it{c})"};
  ConfigurableAxis axisCent{"axisCent", {80, 0.f, 80.f}, "Centrality (%)"};
  ConfigurableAxis axisV1{"axisV1", {1000, -1.0f, 1.0f}, "v1"};
  ConfigurableAxis axisEta{"axisEta", {8, -0.8f, 0.8f}, "Eta"};
  ConfigurableAxis axisNPairs{"axisNPairs", {500, 0.f, 5000.f}, "Number of K^{+}K^{-} pairs per event"};

  Configurable<int> nEvtMixing{"nEvtMixing", 5, "Number of events to mix"};
  ConfigurableAxis cfgVtxBins{"cfgVtxBins", {10, -10, 10}, "Mixing bins - z-vertex"};
  ConfigurableAxis cfgCentBins{"cfgCentBins", {8, 0.0, 80}, "Mixing bins - centrality"};
  /*Configurable<float> ptMix{"ptMix", 0.2f, "ME: pT bin width"};
  Configurable<float> etaMix{"etaMix", 0.2f, "ME: eta bin width"};
  Configurable<float> phiMix{"phiMix", 0.3f, "ME: phi bin width"};*/

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  static constexpr uint8_t kPID1 = 1u << 0;
  static constexpr uint8_t kPID2 = 1u << 1;
  static constexpr uint8_t kPID3 = 1u << 2;
  static constexpr uint8_t kPID4 = 1u << 3;

  SliceCache cache;

  void init(InitContext const&)
  {
    histos.add("hCentrality", "Centrality distribution", kTH1F, {axisCent});
    histos.add("hSparseSame", "Same-event K^{+}K^{-};M;pT;centrality", kTHnSparseF, {axisInvMass, axisPhiPt, axisCent, axisEta, axisV1});
    histos.add("hSparseMix", "Mixed-event K^{+}K^{-};M;pT;centrality", kTHnSparseF, {axisInvMass, axisPhiPt, axisCent, axisEta, axisV1});
    histos.add("hpQxytpvscent", "hpQxytpvscent", HistType::kTHnSparseF, {axisCent, axisV1}, true);
    histos.add("hpQxpvscent", "hpQxpvscent", HistType::kTHnSparseF, {axisCent, axisV1}, true);
    histos.add("hpQxtvscent", "hpQxtvscent", HistType::kTHnSparseF, {axisCent, axisV1}, true);
    histos.add("hpQypvscent", "hpQypvscent", HistType::kTHnSparseF, {axisCent, axisV1}, true);
    histos.add("hpQytvscent", "hpQytvscent", HistType::kTHnSparseF, {axisCent, axisV1}, true);
    histos.add("hMixpairs", "hMixpairs", HistType::kTHnSparseF, {axisNPairs}, true);
  }

  uint8_t getRequiredPidBit() const
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
        LOGF(warn, "Invalid pidChoice=%d. PID2 will be used.", pidChoice.value);
        return kPID2;
    }
  }

  template <typename T>
  bool passStoredPID(const T& kaon) const
  {
    if (!useStoredPID.value) {
      return true;
    }

    const uint16_t mask = kaon.kaonPidMask();
    const uint8_t requiredBit = getRequiredPidBit();

    if (pidExclusive.value) {
      return mask == requiredBit;
    }

    return (mask & requiredBit) != 0;
  }

  bool passDeepAngle(float px1, float py1, float pz1, float px2, float py2, float pz2) const
  {
    if (!applyDeepAngle.value) {
      return true;
    }

    const double p1 = std::sqrt(px1 * px1 + py1 * py1 + pz1 * pz1);
    const double p2 = std::sqrt(px2 * px2 + py2 * py2 + pz2 * pz2);

    if (p1 <= 0.0 || p2 <= 0.0) {
      return false;
    }

    const double cosAngle = std::clamp((px1 * px2 + py1 * py2 + pz1 * pz2) / (p1 * p2), -1.0, 1.0);
    const double angle = std::acos(cosAngle);

    return angle >= deepAngleCut.value;
  }

  template <typename TPosKaons, typename TNegKaons>
  int fillMixedPairs(const TPosKaons& posGroup,
                     const TNegKaons& negGroup,
                     float centrality, float qxZDCA, float qxZDCC, float qyZDCA, float qyZDCC)
  {
    int nMixedPairs = 0;

    for (const auto& kPlus : posGroup) {
      if (!passStoredPID(kPlus)) {
        continue;
      }

      ROOT::Math::PxPyPzMVector plusVec(
        kPlus.px(),
        kPlus.py(),
        kPlus.pz(),
        o2::constants::physics::MassKPlus);

      for (const auto& kMinus : negGroup) {
        if (!passStoredPID(kMinus)) {
          continue;
        }

        if (!passDeepAngle(kPlus.px(),
                           kPlus.py(),
                           kPlus.pz(),
                           kMinus.px(),
                           kMinus.py(),
                           kMinus.pz())) {
          continue;
        }

        ROOT::Math::PxPyPzMVector minusVec(
          kMinus.px(),
          kMinus.py(),
          kMinus.pz(),
          o2::constants::physics::MassKPlus);

        const auto phiCandidate = plusVec + minusVec;

        const float phiMass = phiCandidate.M();
        const float phiPt = phiCandidate.Pt();
        float etaCand = phiCandidate.Eta();
        double phiCand = phiCandidate.Phi();
        if (phiCand < 0.0) {
          phiCand += 2.0 * TMath::Pi();
        }
        double const cosNPhi = std::cos(phiCand);
        double const sinNPhi = std::sin(phiCand);

        auto ux = cosNPhi; // real part of candidate q vector
        auto uy = sinNPhi; // imaginary part of candidate q vector
        auto oddv1 = ux * (qxZDCA - qxZDCC) + uy * (qyZDCA - qyZDCC);
        // auto evenv1 = ux * (qxZDCA + qxZDCC) + uy * (qyZDCA + qyZDCC);

        if (phiMass < massMin.value || phiMass > massMax.value) {
          continue;
        }

        histos.fill(HIST("hSparseMix"), phiMass, phiPt, centrality, etaCand, oddv1);

        nMixedPairs++;
      }
    }

    return nMixedPairs;
  }

  Filter centralityFilter = (aod::kaonevent::cent >= centMin) && (aod::kaonevent::cent < centMax);

  using EventCandidates = soa::Filtered<aod::KaonEvents>;

  Partition<aod::KaonTracks> posKaons = aod::kaonpair::charge > int8_t{0};
  Partition<aod::KaonTracks> negKaons = aod::kaonpair::charge < int8_t{0};

  void processSameData(EventCandidates::iterator const& collision, aod::KaonTracks const& /*kaontracks*/)
  {
    const float centrality = collision.cent();
    auto qxZDCA = collision.qxA();
    auto qyZDCA = collision.qyA();
    auto qxZDCC = collision.qxC();
    auto qyZDCC = collision.qyC();

    auto qxtQxp = qxZDCC * qxZDCA;
    auto qytQyp = qyZDCC * qyZDCA;
    auto qxytp = qxtQxp + qytQyp;
    // auto qxpQyt = qxZDCA * qyZDCC;
    // auto qxtQyp = qxZDCC * qyZDCA;

    histos.fill(HIST("hCentrality"), centrality);
    histos.fill(HIST("hpQxytpvscent"), centrality, qxytp);
    histos.fill(HIST("hpQxpvscent"), centrality, qxZDCA);
    histos.fill(HIST("hpQxtvscent"), centrality, qxZDCC);
    histos.fill(HIST("hpQypvscent"), centrality, qyZDCA);
    histos.fill(HIST("hpQytvscent"), centrality, qyZDCC);

    auto posThisColl = posKaons->sliceByCached(aod::kaonpair::kaoneventId, collision.globalIndex(), cache);
    auto negThisColl = negKaons->sliceByCached(aod::kaonpair::kaoneventId, collision.globalIndex(), cache);
    /*
    int nPlus = 0;
    int nMinus = 0;
    int nPhiSame = 0;

    for (const auto& kPlus : posThisColl) {
      if (passStoredPID(kPlus)) {
        nPlus++;
      }
    }

    for (const auto& kMinus : negThisColl) {
      if (passStoredPID(kMinus)) {
        nMinus++;
      }
    }

    histos.fill(HIST("hNplusPerEvent"), nPlus);
    histos.fill(HIST("hNminusPerEvent"), nMinus);
    */
    for (const auto& kPlus : posThisColl) {
      if (!passStoredPID(kPlus)) {
        continue;
      }

      ROOT::Math::PxPyPzMVector plusVec(kPlus.px(), kPlus.py(), kPlus.pz(), o2::constants::physics::MassKPlus);

      for (const auto& kMinus : negThisColl) {
        if (!passStoredPID(kMinus)) {
          continue;
        }

        if (!passDeepAngle(kPlus.px(), kPlus.py(), kPlus.pz(), kMinus.px(), kMinus.py(), kMinus.pz())) {
          continue;
        }

        ROOT::Math::PxPyPzMVector minusVec(kMinus.px(), kMinus.py(), kMinus.pz(), o2::constants::physics::MassKPlus);

        const auto phiCandidate = plusVec + minusVec;

        const float phiMass = phiCandidate.M();
        const float phiPt = phiCandidate.Pt();
        const float etaCand = phiCandidate.Eta();
        double phiCand = phiCandidate.Phi();
        if (phiCand < 0.0) {
          phiCand += 2.0 * TMath::Pi();
        }
        double const cosNPhi = std::cos(phiCand);
        double const sinNPhi = std::sin(phiCand);

        auto ux = cosNPhi; // real part of candidate q vector
        auto uy = sinNPhi; // imaginary part of candidate q vector
        auto oddv1 = ux * (qxZDCA - qxZDCC) + uy * (qyZDCA - qyZDCC);
        // auto evenv1 = ux * (qxZDCA + qxZDCC) + uy * (qyZDCA + qyZDCC);

        if (phiMass < massMin || phiMass > massMax) {
          continue;
        }

        histos.fill(HIST("hSparseSame"), phiMass, phiPt, centrality, etaCand, oddv1);
        // nPhiSame++;
      }
    }
  }

  PROCESS_SWITCH(phiflowder, processSameData, "Process same-event K+K- pairs", true);

  void processLikeSign(EventCandidates::iterator const& collision, aod::KaonTracks const& /*kaontracks*/)
  {
    const float centrality = collision.cent();

    auto posThisColl = posKaons->sliceByCached(aod::kaonpair::kaoneventId, collision.globalIndex(), cache);
    auto negThisColl = negKaons->sliceByCached(aod::kaonpair::kaoneventId, collision.globalIndex(), cache);
    /*
    int nPlus = 0;
    int nMinus = 0;
    int nPhiLike = 0;

    for (const auto& kPlus : posThisColl) {
      if (passStoredPID(kPlus)) {
        nPlus++;
      }
    }

    for (const auto& kMinus : negThisColl) {
      if (passStoredPID(kMinus)) {
        nMinus++;
      }
    }

    histos.fill(HIST("hNplusLikePerEvent"), nPlus);
    histos.fill(HIST("hNminusLikePerEvent"), nMinus);
    */

    // K+K+ pairs
    for (const auto& kPlus1 : posThisColl) {
      if (!passStoredPID(kPlus1)) {
        continue;
      }

      ROOT::Math::PxPyPzMVector plusVec1(kPlus1.px(), kPlus1.py(), kPlus1.pz(), o2::constants::physics::MassKPlus);

      for (const auto& kPlus2 : posThisColl) {
        if (!passStoredPID(kPlus2)) {
          continue;
        }

        // Reject self-pairs and count every unique K+K+ pair only once.
        if (kPlus2.kaonIndex() <= kPlus1.kaonIndex()) {
          continue;
        }

        if (!passDeepAngle(kPlus1.px(), kPlus1.py(), kPlus1.pz(), kPlus2.px(), kPlus2.py(), kPlus2.pz())) {
          continue;
        }

        ROOT::Math::PxPyPzMVector plusVec2(kPlus2.px(), kPlus2.py(), kPlus2.pz(), o2::constants::physics::MassKPlus);

        const auto likeCandidate = plusVec1 + plusVec2;
        const float likeMass = likeCandidate.M();
        const float likePt = likeCandidate.Pt();

        if (likeMass < massMin.value || likeMass > massMax.value) {
          continue;
        }
        histos.fill(HIST("hSparseLike"), likeMass, likePt, centrality);
        // nPhiLike++;
      }
    }

    // K-K- pairs
    for (const auto& kMinus1 : negThisColl) {
      if (!passStoredPID(kMinus1)) {
        continue;
      }

      ROOT::Math::PxPyPzMVector minusVec1(kMinus1.px(), kMinus1.py(), kMinus1.pz(), o2::constants::physics::MassKPlus);

      for (const auto& kMinus2 : negThisColl) {
        if (!passStoredPID(kMinus2)) {
          continue;
        }

        // Reject self-pairs and count every unique K-K- pair only once.
        if (kMinus2.kaonIndex() <= kMinus1.kaonIndex()) {
          continue;
        }

        if (!passDeepAngle(kMinus1.px(), kMinus1.py(), kMinus1.pz(), kMinus2.px(), kMinus2.py(), kMinus2.pz())) {
          continue;
        }

        ROOT::Math::PxPyPzMVector minusVec2(kMinus2.px(), kMinus2.py(), kMinus2.pz(), o2::constants::physics::MassKPlus);

        const auto likeCandidate = minusVec1 + minusVec2;
        const float likeMass = likeCandidate.M();
        const float likePt = likeCandidate.Pt();

        if (likeMass < massMin.value || likeMass > massMax.value) {
          continue;
        }

        histos.fill(HIST("hSparseLike"), likeMass, likePt, centrality);
        // nPhiLike++;
      }
    }
  }

  PROCESS_SWITCH(phiflowder, processLikeSign, "Process like-sign K+K+ and K-K- pairs", true);

  // Processing Event Mixing
  using BinningType = ColumnBinningPolicy<aod::kaonevent::Posz, aod::kaonevent::Cent>;
  BinningType colBinning{{cfgVtxBins, cfgCentBins}, true};

  void processMixedData(EventCandidates const& collisions,
                        aod::KaonTracks const& /*kaontracks*/)
  {
    BinningType colBinning{{cfgVtxBins, cfgCentBins}, true};
    for (const auto& [collision1, collision2] : selfCombinations(colBinning, nEvtMixing.value, -1, collisions, collisions)) {

      if (collision1.globalIndex() == collision2.globalIndex()) {
        continue;
      }

      const float centrality1 = collision1.cent();
      const float centrality2 = collision2.cent();

      auto qxZDCA = collision1.qxA();
      auto qyZDCA = collision1.qyA();
      auto qxZDCC = collision1.qxC();
      auto qyZDCC = collision1.qyC();

      // K+ from event 1 and K- from event 2
      auto posGroup1 = posKaons->sliceByCached(
        aod::kaonpair::kaoneventId,
        collision1.globalIndex(),
        cache);

      auto negGroup2 = negKaons->sliceByCached(
        aod::kaonpair::kaoneventId,
        collision2.globalIndex(),
        cache);

      // K+ from event 2 and K- from event 1
      auto posGroup2 = posKaons->sliceByCached(
        aod::kaonpair::kaoneventId,
        collision2.globalIndex(),
        cache);

      auto negGroup1 = negKaons->sliceByCached(
        aod::kaonpair::kaoneventId,
        collision1.globalIndex(),
        cache);

      int nMixedPairs = 0;

      nMixedPairs += fillMixedPairs(posGroup1,
                                    negGroup2,
                                    centrality1, qxZDCA, qxZDCC, qyZDCA, qyZDCC);

      nMixedPairs += fillMixedPairs(posGroup2,
                                    negGroup1,
                                    centrality2, qxZDCA, qxZDCC, qyZDCA, qyZDCC);

      histos.fill(HIST("hMixpairs"), nMixedPairs);
    }
  }

  PROCESS_SWITCH(phiflowder, processMixedData, "Process mixed-event K+K- pairs", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<phiflowder>(cfgc)};
}
