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
#include <Framework/ASoA.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/SliceCache.h>
#include <Framework/runDataProcessing.h>

#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
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
  ConfigurableAxis axisNKaons{"axisNKaons", {300, 0.f, 300.f}, "Number of stored kaons per event"};
  ConfigurableAxis axisNPairs{"axisNPairs", {500, 0.f, 5000.f}, "Number of K^{+}K^{-} pairs per event"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  static constexpr uint8_t kPID1 = 1u << 0;
  static constexpr uint8_t kPID2 = 1u << 1;
  static constexpr uint8_t kPID3 = 1u << 2;
  static constexpr uint8_t kPID4 = 1u << 3;

  SliceCache cache;

  void init(InitContext const&)
  {
    histos.add("hCentrality", "Centrality distribution", kTH1F, {axisCent});
    histos.add("hPhiMassSame", "Same-event K^{+}K^{-} invariant mass", kTH1F, {axisInvMass});
    histos.add("hPhiPtSame", "Same-event K^{+}K^{-} pT", kTH1F, {axisPhiPt});
    histos.add("hNplusPerEvent", "Number of stored K^{+} per event", kTH1F, {axisNKaons});
    histos.add("hNminusPerEvent", "Number of stored K^{-} per event", kTH1F, {axisNKaons});
    histos.add("hNphiSamePerEvent", "Number of same-event K^{+}K^{-} pairs per event", kTH1F, {axisNPairs});
    histos.add("hSparseSame", "Same-event K^{+}K^{-};M;pT;centrality", kTHnSparseF, {axisInvMass, axisPhiPt, axisCent});
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

  Filter centralityFilter = (aod::kaonevent::cent >= centMin) && (aod::kaonevent::cent < centMax);

  using EventCandidates = soa::Filtered<aod::KaonEvents>;

  Partition<aod::KaonTracks> posKaons = aod::kaonpair::charge > int8_t{0};
  Partition<aod::KaonTracks> negKaons = aod::kaonpair::charge < int8_t{0};

  void processSameData(EventCandidates::iterator const& collision, aod::KaonTracks const& /*kaontracks*/)
  {
    const float centrality = collision.cent();

    histos.fill(HIST("hCentrality"), centrality);

    auto posThisColl = posKaons->sliceByCached(aod::kaonpair::kaoneventId, collision.globalIndex(), cache);
    auto negThisColl = negKaons->sliceByCached(aod::kaonpair::kaoneventId, collision.globalIndex(), cache);

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

        if (phiMass < massMin || phiMass > massMax) {
          continue;
        }

        histos.fill(HIST("hPhiMassSame"), phiMass);
        histos.fill(HIST("hPhiPtSame"), phiPt);
        histos.fill(HIST("hSparseSame"), phiMass, phiPt, centrality);

        nPhiSame++;
      }
    }

    histos.fill(HIST("hNphiSamePerEvent"), nPhiSame);
  }

  PROCESS_SWITCH(phiflowder, processSameData, "Process same-event K+K- pairs", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<phiflowder>(cfgc)};
}
