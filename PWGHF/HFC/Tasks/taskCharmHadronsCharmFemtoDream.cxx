// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file taskCharmHadronsCharmFemtoDream.cxx
/// \brief Build DD pair femtoscopy distributions
/// \author Biao Zhang, Heidelberg University, biao.zhang@cern.ch

#include "PWGCF/DataModel/FemtoDerived.h"
#include "PWGCF/FemtoDream/Core/femtoDreamMath.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/runDataProcessing.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace o2::analysis::femtoDream;

struct HfTaskCharmHadronsCharmFemtoDream {
  // Like/unlike sign refers to charm flavour; D0 is electrically neutral.
  // Keep charge-conjugate channels separate so they can also be combined offline.
  enum PairChannel : int {
    kD0D0LikeSign = 0,
    kD0BarD0BarLikeSign = 1,
    kD0D0BarUnlikeSign = 2,
    kD0DstarPlusLikeSign = 3,
    kD0BarDstarMinusLikeSign = 4,
    kD0DstarMinusUnlikeSign = 5,
    kD0BarDstarPlusUnlikeSign = 6,
    kNPairChannels = 7
  };

  template <bool IsDstar, typename D0Row, typename OtherRow>
  static bool sharesDaughter(D0Row const& a, OtherRow const& b)
  {
    for (auto id : std::array{a.prong0Id(), a.prong1Id()}) {
      if (id < 0) {
        continue;
      }
      if (id == b.prong0Id() || id == b.prong1Id()) {
        return true;
      }
      if constexpr (IsDstar) {
        if (id == b.prong2Id()) {
          return true;
        }
      }
    }
    return false;
  }

  // a is D0; b is either D0 or Dstar. D0 charge means flavour, not electric
  // charge.
  template <bool IsDstar, typename D0Row, typename OtherRow>
  static PairChannel pairChannel(D0Row const& a, OtherRow const& b)
  {
    if constexpr (!IsDstar) {
      return a.charge() != b.charge() ? kD0D0BarUnlikeSign : (a.charge() > 0 ? kD0D0LikeSign : kD0BarD0BarLikeSign);
    }
    return a.charge() > 0 ? (b.charge() > 0 ? kD0DstarPlusLikeSign : kD0DstarMinusUnlikeSign)
                          : (b.charge() < 0 ? kD0BarDstarMinusLikeSign : kD0BarDstarPlusUnlikeSign);
  }

  // Canonical ordering keeps mass/pt axes independent of event ordering for D0D0.
  template <typename FirstRow, typename SecondRow>
  static bool reverseD0Order(FirstRow const& a, SecondRow const& b)
  {
    if (a.charge() != b.charge()) {
      return a.charge() < b.charge(); // D0 first for D0-D0bar
    }
    return a.pt() < b.pt(); // leading pT first for equal flavours
  }

  Configurable<float> ptMinD0{"ptMinD0", 0.f, "Minimum D0 pT"};
  Configurable<float> ptMaxD0{"ptMaxD0", 36.f, "Maximum D0 pT"};
  Configurable<float> ptMinDstar{"ptMinDstar", 0.f, "Minimum Dstar pT"};
  Configurable<float> ptMaxDstar{"ptMaxDstar", 36.f, "Maximum Dstar pT"};
  Configurable<float> etaMax{"etaMax", 0.8f, "Maximum absolute candidate eta"};
  Configurable<float> massMinD0{"massMinD0", 1.7f, "Minimum D0 mass (keep sidebands)"};
  Configurable<float> massMaxD0{"massMaxD0", 2.0f, "Maximum D0 mass"};
  Configurable<float> deltaMassMin{"deltaMassMin", 0.139f, "Minimum Dstar-D0 mass difference"};
  Configurable<float> deltaMassMax{"deltaMassMax", 0.17f, "Maximum Dstar-D0 mass difference"};
  Configurable<float> daughterMassMin{"daughterMassMin", 1.80f,"Minimum Dstar daughter D0 mass"};
  Configurable<float> daughterMassMax{"daughterMassMax", 1.93f,"Maximum Dstar daughter D0 mass"};
  Configurable<bool> useMl{"useMl", false, "Require valid ML scores and apply score cuts"};
  Configurable<float> maxBkgD0{"maxBkgD0", 1.f, "Maximum D0 background score"};
  Configurable<float> minPromptD0{"minPromptD0", 0.f, "Minimum D0 prompt score"};
  Configurable<float> maxBkgDstar{"maxBkgDstar", 1.f,"Maximum Dstar background score"};
  Configurable<float> minPromptDstar{"minPromptDstar", 0.f, "Minimum Dstar prompt score"};
  Configurable<int> charmHadCandSel{"charmHadCandSel", 1, "Minimum reduced charm candidate selection flag"};
  struct : ConfigurableGroup {
    std::string prefix = "eventSel";
    Configurable<bool> useCentrality{"useCentrality", false, "Apply percentile selection (requires a centrality-enabled producer)"};
    Configurable<int> multMin{"multMin", 0, "Minimum multiplicity (MultNtr)"};
    Configurable<int> multMax{"multMax", 99999, "Maximum multiplicity (MultNtr)"};
    Configurable<float> multPercentileMin{"multPercentileMin", 0.f, "Minimum multiplicity percentile"};
    Configurable<float> multPercentileMax{"multPercentileMax", 100.f, "Maximum multiplicity percentile"};
  } eventSel;

  struct : ConfigurableGroup {
    Configurable<bool> doMixEvent{"doMixEvent", true, "Enable mixed events"};
    Configurable<int> mixingBinPolicy{"mixingBinPolicy", 0, "0: multiplicity, 1: percentile, 2: both"};
    Configurable<int> mixingDepth{"mixingDepth", 5, "Number of neighbours for O2 event mixing"};
  } mixSetting;
  ConfigurableAxis mixingBinMult{"mixingBinMult", {VARIABLE_WIDTH, 0.f, 20.f, 60.f, 200.f}, "Mixing bins - multiplicity"};
  ConfigurableAxis mixingBinMultPercentile{"mixingBinMultPercentile", {VARIABLE_WIDTH, 0.f, 100.f}, "Mixing bins - multiplicity percentile"};
  ConfigurableAxis mixingBinVztx{"mixingBinVztx", {VARIABLE_WIDTH, -10.f, -4.f, 0.f, 4.f, 10.f}, "Mixing bins - z-vertex"};
  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultNtr> colBinningMult{{mixingBinVztx, mixingBinMult}, true};
  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultV0M> colBinningMultPercentile{{mixingBinVztx, mixingBinMultPercentile}, true};
  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultNtr, aod::femtodreamcollision::MultV0M> colBinningMultMultPercentile{{mixingBinVztx, mixingBinMult, mixingBinMultPercentile}, true};
  aod::femtodreamcollision::BitMaskType bitMask = 1 << 0;
  using FilteredCollisions = soa::Filtered<soa::Join<aod::FDCollisions, aod::FDColMasks>>;
  using FilteredCharmCand2Prongs = soa::Filtered<aod::FDHfCand2Prong>;
  using FilteredCharmCandDstars = soa::Filtered<aod::FDHfCandDstar>;

  Filter eventMultiplicity = aod::femtodreamcollision::multNtr >= eventSel.multMin && aod::femtodreamcollision::multNtr <= eventSel.multMax;
  Filter eventMultiplicityPercentile = ifnode(eventSel.useCentrality, aod::femtodreamcollision::multV0M >= eventSel.multPercentileMin && aod::femtodreamcollision::multV0M <= eventSel.multPercentileMax, Node{LiteralNode{true}});
  Filter hfCandSelFilter = aod::fdhf::candidateSelFlag >= charmHadCandSel;

  Partition<FilteredCharmCand2Prongs> partitionCharmHadron2Prong = ifnode(useMl, aod::fdhf::bdtBkg <= maxBkgD0 && aod::fdhf::bdtPrompt >= minPromptD0, Node{LiteralNode{true}});
  Partition<FilteredCharmCandDstars> partitionCharmHadronDstar = ifnode(useMl, aod::fdhf::bdtBkg <= maxBkgDstar && aod::fdhf::bdtPrompt >= minPromptDstar, Node{LiteralNode{true}});

  // Full-table process inputs do not automatically register grouping caches.
  // Register both keys explicitly for the sliceByCached calls below.
  Preslice<FilteredCharmCand2Prongs> perCollisionD0 = aod::femtodreamparticle::fdCollisionId;
  Preslice<FilteredCharmCandDstars> perCollisionDstar = aod::femtodreamparticle::fdCollisionId;
  SliceCache cache;
  HistogramRegistry registry{
    "Results",
    {},
    OutputObjHandlingPolicy::AnalysisObject};
  Produces<aod::FDHfCharm2Prong> rowFemtoResultCharm2Prong;
  Produces<aod::FDHfCharmDstar> rowFemtoResultCharmDstar;
  Produces<aod::FDHfColl> rowFemtoResultColl;

  void init(InitContext const&)
  {
    if (static_cast<int>(doprocessD0D0) + doprocessD0Dstar != 1) {
      LOGP(fatal, "Enable exactly one charm-charm analysis process");
    }
    if (mixSetting.mixingDepth < 0 || mixSetting.mixingBinPolicy < 0 || mixSetting.mixingBinPolicy > 2 || ptMinD0 < 0 || ptMinD0 >= ptMaxD0 ||
        ptMinDstar < 0 || ptMinDstar >= ptMaxDstar || etaMax <= 0 ||
        massMinD0 >= massMaxD0 || deltaMassMin >= deltaMassMax ||
        daughterMassMin >= daughterMassMax || charmHadCandSel < 1 ||
        eventSel.multMin > eventSel.multMax || eventSel.multPercentileMin > eventSel.multPercentileMax) {
      LOGP(fatal, "Invalid charm-charm selection or mixing configuration");
    }
    colBinningMult = {{mixingBinVztx, mixingBinMult}, true};
    colBinningMultPercentile = {{mixingBinVztx, mixingBinMultPercentile}, true};
    colBinningMultMultPercentile = {{mixingBinVztx, mixingBinMult, mixingBinMultPercentile}, true};
    const AxisSpec kstar{400, 0., 2., "k* (GeV/c)"};
    const AxisSpec massD0{300, massMinD0.value, massMaxD0.value,
                          "M(Kpi) (GeV/c2)"};
    const AxisSpec deltaMass{310, deltaMassMin.value, deltaMassMax.value,
                             "Delta M (GeV/c2)"};
    const AxisSpec pt{72, ptMinD0.value, ptMaxD0.value, "D0 pT (GeV/c)"};
    const AxisSpec ptStar{72, ptMinDstar.value, ptMaxDstar.value, "Dstar pT (GeV/c)"};
    const AxisSpec mult{mixingBinMult, "NTracksPV"};
    const AxisSpec channel{
      static_cast<int>(kNPairChannels), -0.5, static_cast<double>(kNPairChannels) - 0.5,
      "0=D0D0 LS,1=D0barD0bar LS,2=D0D0bar US,3=D0D*+ LS,4=D0barD*- LS,5=D0D*- US,6=D0barD*+ US"};
    registry.add("SE/D0D0", "Same event", kTHnSparseF,
                 {kstar, massD0, massD0, pt, pt, mult, channel});
    registry.add("ME/D0D0", "Mixed event", kTHnSparseF,
                 {kstar, massD0, massD0, pt, pt, mult, channel});
    registry.add("SE/D0Dstar", "Same event", kTHnSparseF,
                 {kstar, massD0, deltaMass, pt, ptStar, mult, channel});
    registry.add("ME/D0Dstar", "Mixed event", kTHnSparseF,
                 {kstar, massD0, deltaMass, pt, ptStar, mult, channel});
    registry.add("QA/events", ";0=in mixing range,1=outside;events", kTH1F,
                 {{2, -0.5, 1.5}});
    registry.add("MixingQA/hSECollisionBins", ";mixing bin;events", kTH1F, {{1000, -0.5, 999.5}});
    registry.add("MixingQA/hMECollisionBinsD0D0", ";mixing bin;event pairs", kTH1F, {{1000, -0.5, 999.5}});
    registry.add("MixingQA/hMECollisionBinsD0Dstar", ";mixing bin;event pairs", kTH1F, {{1000, -0.5, 999.5}});
    registry.add("QA/sharedDaughters", ";0=D0D0,1=D0Dstar;rejected pairs",
                 kTH1F, {{2, -0.5, 1.5}});
    registry.add("QA/d0MassPt", ";M(Kpi);pT", kTH2F, {massD0, pt});
    registry.add("QA/dstarMassPt", ";Delta M;pT", kTH2F, {deltaMass, ptStar});
  }

  template <typename Candidates>
  void validateMlScores(Candidates const& candidates)
  {
    if (!useMl) {
      return;
    }
    for (auto const& row : candidates) {
      if (!std::isfinite(row.bdtBkg()) || !std::isfinite(row.bdtPrompt()) ||
          row.bdtBkg() < 0.f || row.bdtPrompt() < 0.f) {
        LOGP(fatal, "ML cuts requested on missing/invalid scores: use an ML producer process");
      }
    }
  }

  template <bool IsDstar, bool FillQa = false, typename Row>
  bool select(Row const& row, float& mass)
  {
    if (std::abs(row.charge()) != 1 || !std::isfinite(row.pt()) ||
        !std::isfinite(row.eta()) || !std::isfinite(row.phi()) ||
        std::abs(row.eta()) >= etaMax ||
        row.pt() < (IsDstar ? ptMinDstar.value : ptMinD0.value) ||
        row.pt() >= (IsDstar ? ptMaxDstar.value : ptMaxD0.value)) {
      return false;
    }
    const std::array<double, 2> masses =
      row.charge() > 0 ? std::array{MassPiPlus, MassKPlus}
                       : std::array{MassKPlus, MassPiPlus};
    if constexpr (IsDstar) {
      const float daughterMass = row.mDaughD0(masses);
      if (!std::isfinite(daughterMass) || daughterMass < daughterMassMin ||
          daughterMass >= daughterMassMax) {
        return false;
      }
      mass = row.m(std::array{masses[0], masses[1], MassPiPlus}) - daughterMass;
    } else {
      mass = row.m(masses);
    }
    if (!std::isfinite(mass) ||
        mass < (IsDstar ? deltaMassMin.value : massMinD0.value) ||
        mass >= (IsDstar ? deltaMassMax.value : massMaxD0.value)) {
      return false;
    }
    if constexpr (IsDstar) {
      if constexpr (FillQa) {
        registry.fill(HIST("QA/dstarMassPt"), mass, row.pt());
      }
    } else {
      if constexpr (FillQa) {
        registry.fill(HIST("QA/d0MassPt"), mass, row.pt());
      }
    }
    return true;
  }

  template <bool Mixed, bool IsDstar, typename FirstRow, typename SecondRow>
  void fillPair(FirstRow const& first, SecondRow const& second,
                float firstMass, float secondMass, float mult)
  {
    if constexpr (!Mixed) {
      if (sharesDaughter<IsDstar>(first, second)) {
        registry.fill(HIST("QA/sharedDaughters"), IsDstar ? 1 : 0);
        return;
      }
    }
    const bool reverse = !IsDstar && reverseD0Order(first, second);
    const float mass1 = reverse ? secondMass : firstMass;
    const float mass2 = reverse ? firstMass : secondMass;
    const float pt1 = reverse ? second.pt() : first.pt();
    const float pt2 = reverse ? first.pt() : second.pt();
    // Read kinematics directly from the candidate rows. Only D0D0 can reverse;
    // its two nominal parent masses are equal. Never use Dstar delta mass here.
    const float kstar = reverse
                          ? FemtoDreamMath::getkstar(second, MassD0, first, MassD0)
                          : FemtoDreamMath::getkstar(first, MassD0, second, IsDstar ? MassDStar : MassD0);
    if (!std::isfinite(kstar)) {
      return;
    }
    const int channel = pairChannel<IsDstar>(first, second);
    if constexpr (IsDstar && Mixed) {
      registry.fill(HIST("ME/D0Dstar"), kstar, mass1, mass2, pt1, pt2, mult, channel);
    } else if constexpr (IsDstar) {
      registry.fill(HIST("SE/D0Dstar"), kstar, mass1, mass2, pt1, pt2, mult, channel);
    } else if constexpr (Mixed) {
      registry.fill(HIST("ME/D0D0"), kstar, mass1, mass2, pt1, pt2, mult, channel);
    } else {
      registry.fill(HIST("SE/D0D0"), kstar, mass1, mass2, pt1, pt2, mult, channel);
    }
  }

  template <typename Policy, typename Collision>
  int mixingBin(Policy const& policy, Collision col)
  {
    return policy.getBin(policy.getBinningValues(col));
  }

  template <bool IsDstar, typename D0Table, typename OtherTable, typename Policy>
  void doMixedEvent(FilteredCollisions const& cols, D0Table const& d0Rows,
                    OtherTable const& otherRows, Policy const& policy)
  {
    Partition<FilteredCollisions> partitionMaskedCol1 = (aod::femtodreamcollision::bitmaskTrackOne & bitMask) == bitMask;
    // TrackOne denotes D0 presence, TrackTwo denotes Dstar presence.
    Partition<FilteredCollisions> partitionMaskedCol2 = (aod::femtodreamcollision::bitmaskTrackOne & bitMask) == bitMask;
    Partition<FilteredCollisions> partitionMaskedColDstar = (aod::femtodreamcollision::bitmaskTrackTwo & bitMask) == bitMask;
    partitionMaskedCol1.bindTable(cols);
    if constexpr (IsDstar) {
      partitionMaskedColDstar.bindTable(cols);
    } else {
      partitionMaskedCol2.bindTable(cols);
    }
    auto const& secondCollisions = IsDstar ? *partitionMaskedColDstar.mFiltered : *partitionMaskedCol2.mFiltered;
    for (auto const& [collision1, collision2] : combinations(soa::CombinationsBlockFullIndexPolicy(policy, mixSetting.mixingDepth, -1, *partitionMaskedCol1.mFiltered, secondCollisions))) {
      if (collision1.globalIndex() == collision2.globalIndex()) {
        continue;
      }
      if constexpr (!IsDstar) {
        // Full policy returns both orientations for identical inputs. Count once.
        if (collision1.globalIndex() > collision2.globalIndex()) {
          continue;
        }
      }
      if (!std::isfinite(collision1.magField()) || collision1.magField() != collision2.magField()) {
        continue;
      }
      if constexpr (IsDstar) {
        registry.fill(HIST("MixingQA/hMECollisionBinsD0Dstar"), mixingBin(policy, collision1));
      } else {
        registry.fill(HIST("MixingQA/hMECollisionBinsD0D0"), mixingBin(policy, collision1));
      }
      auto first = d0Rows.sliceByCached(aod::femtodreamparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto second = otherRows.sliceByCached(aod::femtodreamparticle::fdCollisionId, collision2.globalIndex(), cache);
      // Species roles remain fixed; the full policy supplies the reverse event
      // orientation when eligible. Do not add a second manual reverse loop.
      for (auto const& [row1, row2] : combinations(soa::CombinationsFullIndexPolicy(first, second))) {
        float mass1{}, mass2{};
        if (select<false>(row1, mass1) && select<IsDstar>(row2, mass2)) {
          fillPair<true, IsDstar>(row1, row2, mass1, mass2, collision1.multNtr());
        }
      }
    }
  }

  // Write candidates once per collision, independently of SE/ME pair counts.
  // Reuse the same output tables and field conventions as the track/V0 tasks.
  template <bool WithDstar, typename CollType, typename D0Slice, typename DstarSlice>
  void fillTables(CollType const& col, D0Slice const& d0s, DstarSlice const& dstars)
  {
    int64_t timeStamp = -1;
    bool hasCandidate = false;
    auto recordTimeStamp = [&](int64_t value) {
      if (hasCandidate && value != timeStamp) {
        LOGP(fatal, "Inconsistent charm timestamps in one reduced collision");
      }
      timeStamp = value;
      hasCandidate = true;
    };
    for (auto const& part : d0s) {
      float mass{};
      if (!select<false, true>(part, mass)) {
        continue;
      }
      recordTimeStamp(part.timeStamp());
      rowFemtoResultCharm2Prong(
        col.globalIndex(), timeStamp, mass,
        part.pt(), part.eta(), part.phi(), part.prong0Id(), part.prong1Id(),
        part.charge(), part.bdtBkg(), part.bdtPrompt(), part.bdtFD(), 0, 0);
    }
    if constexpr (WithDstar) {
      for (auto const& part : dstars) {
        float mass{};
        if (!select<true, true>(part, mass)) {
          continue;
        }
        recordTimeStamp(part.timeStamp());
        const std::array<double, 2> daughterMasses = part.charge() > 0 ? std::array{MassPiPlus, MassKPlus} : std::array{MassKPlus, MassPiPlus};
        // CharmM is delta mass, CharmDaughM is M(Kpi), as in D+Track.
        rowFemtoResultCharmDstar(
          col.globalIndex(), timeStamp, mass,
          part.mDaughD0(daughterMasses), part.pt(), part.eta(), part.phi(),
          part.prong0Id(), part.prong1Id(), part.prong2Id(), part.charge(),
          part.bdtBkg(), part.bdtPrompt(), part.bdtFD(), 0, 0);
      }
    }
    // Retain D0-only/Dstar-only events even when no clean SE pair exists.
    // MC fields above are zero for these data processes, as in existing tasks.
    if (hasCandidate) {
      rowFemtoResultColl(col.globalIndex(), timeStamp, col.posZ(), col.multNtr());
    }
  }

  template <bool WithDstar, typename D0Table, typename DstarTable, typename Policy>
  void runWithPolicy(FilteredCollisions const& cols, D0Table const& d0Rows,
                     DstarTable const& dstarRows, Policy const& policy)
  {
    for (auto const& col : cols) {
      auto d0s = d0Rows.sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
      // Export before online mixing-bin cuts, so offline mixing can change bins.
      if constexpr (WithDstar) {
        auto dstars = dstarRows.sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
        fillTables<true>(col, d0s, dstars);
      } else {
        fillTables<false>(col, d0s, d0s);
      }
      const int bin = mixingBin(policy, col);
      if (bin < 0 || !std::isfinite(col.magField())) {
        registry.fill(HIST("QA/events"), 1);
        continue;
      }
      registry.fill(HIST("QA/events"), 0);
      registry.fill(HIST("MixingQA/hSECollisionBins"), bin);
      for (auto const& [row1, row2] : combinations(soa::CombinationsStrictlyUpperIndexPolicy(d0s, d0s))) {
        float mass1{}, mass2{};
        if (select<false>(row1, mass1) && select<false>(row2, mass2)) {
          fillPair<false, false>(row1, row2, mass1, mass2, col.multNtr());
        }
      }
      if constexpr (WithDstar) {
        auto dstars = dstarRows.sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
        for (auto const& [row1, row2] : combinations(soa::CombinationsFullIndexPolicy(d0s, dstars))) {
          float mass1{}, mass2{};
          if (select<false>(row1, mass1) && select<true>(row2, mass2)) {
            fillPair<false, true>(row1, row2, mass1, mass2, col.multNtr());
          }
        }
      }
    }
    if (mixSetting.doMixEvent && mixSetting.mixingDepth > 0) {
      doMixedEvent<false>(cols, d0Rows, d0Rows, policy);
      if constexpr (WithDstar) {
        doMixedEvent<true>(cols, d0Rows, dstarRows, policy);
      }
    }
  }

  template <bool WithDstar, typename DstarTable>
  void run(FilteredCollisions const& cols, FilteredCharmCand2Prongs const& d0Rows, DstarTable const& dstarRows)
  {
    validateMlScores(d0Rows);
    partitionCharmHadron2Prong.bindTable(d0Rows);
    if constexpr (WithDstar) {
      validateMlScores(dstarRows);
      partitionCharmHadronDstar.bindTable(dstarRows);
    }
    switch (mixSetting.mixingBinPolicy) {
      case aod::femtodreamcollision::kMult:
        runWithPolicy<WithDstar>(cols, partitionCharmHadron2Prong, partitionCharmHadronDstar, colBinningMult);
        break;
      case aod::femtodreamcollision::kMultPercentile:
        runWithPolicy<WithDstar>(cols, partitionCharmHadron2Prong, partitionCharmHadronDstar, colBinningMultPercentile);
        break;
      case aod::femtodreamcollision::kMultMultPercentile:
        runWithPolicy<WithDstar>(cols, partitionCharmHadron2Prong, partitionCharmHadronDstar, colBinningMultMultPercentile);
        break;
      default:
        LOGP(fatal, "Invalid mixing binning policy");
    }
  }

  void processD0D0(FilteredCollisions const& cols, FilteredCharmCand2Prongs const& d0s)
  {
    run<false>(cols, d0s, d0s);
  }
  PROCESS_SWITCH(HfTaskCharmHadronsCharmFemtoDream, processD0D0, "D0D0 data, SE and ME", true);

  void processD0Dstar(FilteredCollisions const& cols, FilteredCharmCand2Prongs const& d0s,
                      FilteredCharmCandDstars const& dstars)
  {
    run<true>(cols, d0s, dstars);
  }
  PROCESS_SWITCH(HfTaskCharmHadronsCharmFemtoDream, processD0Dstar, "D0D0 and D0Dstar data, SE and ME", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfTaskCharmHadronsCharmFemtoDream>(cfgc)};
}
