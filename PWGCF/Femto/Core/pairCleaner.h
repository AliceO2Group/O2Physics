// Copyright 2019-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file pairCleaner.h
/// \brief pair cleaner class
/// \author anton.riedel@tum.de, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTO_CORE_PAIRCLEANER_H_
#define PWGCF_FEMTO_CORE_PAIRCLEANER_H_

#include "PWGCF/Femto/Core/histManager.h"
#include "PWGCF/Femto/Core/modes.h"

#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/Logger.h>

#include <TH1.h>

#include <Rtypes.h>

#include <array>
#include <map>
#include <string>
#include <string_view>
#include <vector>

namespace o2::analysis::femto::paircleaner
{

// enum for pair cleaner histograms
enum PairCleanerHist {
  // standard 1D
  kStats,
  kKinematic, // kinematic variable of the pair/triplet
  kPairCleanerHistLast
};

// bin indices (0-based fill values) for the kStats histogram
enum PairCleanerStatsBin {
  kAllAnalyzed = 0,
  kAllBlocked = 1,
  kPairCleanerStatsBinLast = 2,
};

struct ConfPairCleanerBinning : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("PairCleanerBinning");
  o2::framework::ConfigurableAxis kinematic{"kinematic", {{100, 0.f, 1.f}}, "kinematic"};
};

constexpr std::array<histmanager::HistInfo<PairCleanerHist>, kPairCleanerHistLast>
  HistTable = {
    {
      // 1D
      {kStats, o2::framework::HistType::kTH1F, "hStats", "Stats; ; Entries"},
      {kKinematic, o2::framework::HistType::kTH1F, "hKinematic", "Kinematic of blocked pairs; Kinematic variable (GeV/#it{c}); Entries"},
    }};

template <typename T>
auto makePairCleanerHistSpecMap(const T& confPairCleaner)
{
  return std::map<PairCleanerHist, std::vector<o2::framework::AxisSpec>>{
    {kKinematic, {confPairCleaner.kinematic}},
  };
};

constexpr char PrefixPairCleanerTrackTrackSe[] = "TrackTrackCleaner/SE/";
constexpr char PrefixPairCleanerTrackTrackMe[] = "TrackTrackCleaner/ME/";
constexpr char PrefixPairCleanerV0V0Se[] = "V0V0Cleaner/SE/";
constexpr char PrefixPairCleanerV0V0Me[] = "V0V0Cleaner/ME/";
constexpr char PrefixPairCleanerD0D0Se[] = "D0D0Cleaner/SE/";
constexpr char PrefixPairCleanerD0D0Me[] = "D0D0Cleaner/ME/";
constexpr char PrefixPairCleanerTrackV0Se[] = "TrackV0Cleaner/SE/";
constexpr char PrefixPairCleanerTrackV0Me[] = "TrackV0Cleaner/ME/";
constexpr char PrefixPairCleanerTrackD0Se[] = "TrackD0Cleaner/SE/";
constexpr char PrefixPairCleanerTrackD0Me[] = "TrackD0Cleaner/ME/";
constexpr char PrefixPairCleanerTrackLcSe[] = "TrackLcCleaner/SE/";
constexpr char PrefixPairCleanerTrackLcMe[] = "TrackLcCleaner/ME/";
constexpr char PrefixPairCleanerTrackResonanceSe[] = "TrackResonanceCleaner/SE/";
constexpr char PrefixPairCleanerTrackResonanceMe[] = "TrackResonanceCleaner/ME/";
constexpr char PrefixPairCleanerV0ResonanceSe[] = "V0ResonanceCleaner/SE/";
constexpr char PrefixPairCleanerV0ResonanceMe[] = "V0ResonanceCleaner/ME/";
constexpr char PrefixPairCleanerTrackKinkSe[] = "TrackKinkCleaner/SE/";
constexpr char PrefixPairCleanerTrackKinkMe[] = "TrackKinkCleaner/ME/";
constexpr char PrefixPairCleanerTrackCascadeSe[] = "TrackCascadeCleaner/SE/";
constexpr char PrefixPairCleanerTrackCascadeMe[] = "TrackCascadeCleaner/ME/";
constexpr char PrefixPairCleanerMcParticleMcParticleSe[] = "McParticleMcParticleCleaner/SE/";
constexpr char PrefixPairCleanerMcParticleMcParticleMe[] = "McParticleMcParticleCleaner/ME/";

template <auto& prefix>
class BasePairCleaner
{
 public:
  BasePairCleaner() = default;
  virtual ~BasePairCleaner() = default;

  template <modes::Mode mode, typename T1, typename T2>
  void init(o2::framework::HistogramRegistry* registry,
            T1 const& specs,
            T2 const& pairCuts) // pair cuts configurable is defined in pairHistManager.h
  {
    if constexpr (modes::isFlagSet(mode, modes::Mode::kMc)) {
      mMixPairsWithCommonAncestor = pairCuts.mixOnlyCommonAncestor.value;
      mMixPairsWithNonCommonAncestor = pairCuts.mixOnlyNonCommonAncestor.value;
      mUseMotherAsAncestor = pairCuts.useMotherAsAncestor.value;
      if (mMixPairsWithCommonAncestor && mMixPairsWithNonCommonAncestor) {
        LOG(fatal) << "Both mixing with common and non-common ancestor is activated. Breaking...";
      }
    }

    mHistogramRegistry = registry;
    int nStatBins = static_cast<int>(kPairCleanerStatsBinLast);
    const o2::framework::AxisSpec axisStats = {nStatBins, -0.5, nStatBins - 0.5};
    mHistogramRegistry->add(std::string(prefix) + getHistNameV2(kStats, HistTable), getHistDesc(kStats, HistTable), getHistType(kStats, HistTable), {axisStats});
    mHistogramRegistry->get<TH1>(HIST(prefix) + HIST(histmanager::getHistName(kStats, HistTable)))->GetXaxis()->SetBinLabel(1, "All analyzed reco pairs");
    mHistogramRegistry->get<TH1>(HIST(prefix) + HIST(histmanager::getHistName(kStats, HistTable)))->GetXaxis()->SetBinLabel(2, "All blocked reco pairs");
    mHistogramRegistry->add(std::string(prefix) + getHistNameV2(kKinematic, HistTable), getHistDesc(kKinematic, HistTable), getHistType(kKinematic, HistTable), {specs.at(kKinematic)});
  }

 protected:
  template <typename T1, typename T2>
  bool isCleanParticlePair(T1 const& particle1, T2 const& particle2) const
  {
    return particle1.globalIndex() != particle2.globalIndex();
  };

  // const: only mutates *mHistogramRegistry through the stored pointer, not the pointer itself,
  // so it can be called from the const isCleanPair() overloads below.
  void fillAll() const
  {
    mHistogramRegistry->fill(HIST(prefix) + HIST(getHistName(kStats, HistTable)), kAllAnalyzed);
  }

  void fillBlocked(float kinematic) const
  {
    mHistogramRegistry->fill(HIST(prefix) + HIST(getHistName(kStats, HistTable)), kAllBlocked);
    mHistogramRegistry->fill(HIST(prefix) + HIST(getHistName(kKinematic, HistTable)), kinematic);
  }

  // mc only
  // ancestry is checked either with the partonic mother or with the first ancestor (i.e. the direct mother), depending on mUseMotherAsAncestor
  template <typename T1, typename T2, typename T3>
  bool mcPairHasCommonAncestor(T1 const& particle1, T2 const& particle2, T3 const& partonicMothers) const
  {
    if (mUseMotherAsAncestor) {
      return this->mcPairHasCommonMother(particle1, particle2);
    }
    return this->mcPairHasCommonPartonicMother(particle1, particle2, partonicMothers);
  };

  template <typename T1, typename T2, typename T3>
  bool mcPairHasNonCommonAncestor(T1 const& particle1, T2 const& particle2, T3 const& partonicMothers) const
  {
    if (mUseMotherAsAncestor) {
      // if one of the two particles has no associated mother, we cannot know if they have a common anchestor, so we break out with false
      if (!particle1.has_fMcMother() || !particle2.has_fMcMother()) {
        return false;
      }
      return !this->mcPairHasCommonMother(particle1, particle2);
    }
    // if one of the two particles has no associated partonic mother, we cannot know if they have a common anchestor, so we break out with false
    if (!particle1.has_fMcPartMoth() || !particle2.has_fMcPartMoth()) {
      return false;
    }
    return !this->mcPairHasCommonPartonicMother(particle1, particle2, partonicMothers);
  };

  // reco + mc
  template <typename T1, typename T2, typename T3, typename T4>
  bool pairHasCommonAncestor(T1 const& particle1, T2 const& particle2, T3 const& /*mcparticles*/, T4 const& partonicMothers) const
  {
    // if one of the two particles has no associated mc particle, we cannot know if they have a common anchestor, so we break out with false
    if (!particle1.has_fMcParticle() || !particle2.has_fMcParticle()) {
      return false;
    }

    // get mc particles
    auto mcParticle1 = particle1.template fMcParticle_as<T3>();
    auto mcParticle2 = particle2.template fMcParticle_as<T3>();

    return this->mcPairHasCommonAncestor(mcParticle1, mcParticle2, partonicMothers);
  };

  template <typename T1, typename T2, typename T3, typename T4>
  bool pairHasNonCommonAncestor(T1 const& particle1, T2 const& particle2, T3 const& /*mcparticles*/, T4 const& partonicMothers) const
  {
    // if one of the two particles has no associated mc particle, we cannot know if they have a common anchestor, so we break out with false
    if (!particle1.has_fMcParticle() || !particle2.has_fMcParticle()) {
      return false;
    }

    // get mc particles
    auto mcParticle1 = particle1.template fMcParticle_as<T3>();
    auto mcParticle2 = particle2.template fMcParticle_as<T3>();

    return this->mcPairHasNonCommonAncestor(mcParticle1, mcParticle2, partonicMothers);
  };

  bool mMixPairsWithCommonAncestor = false;
  bool mMixPairsWithNonCommonAncestor = false;

 private:
  // require both particles to originate from the same partonic mother
  template <typename T1, typename T2, typename T3>
  bool mcPairHasCommonPartonicMother(T1 const& particle1, T2 const& particle2, T3 const& /*partonicMothers*/) const
  {
    // if one of the two particles has no associated partonic mother, we cannot know if they have a common anchestor, so we break out with false
    if (!particle1.has_fMcPartMoth() || !particle2.has_fMcPartMoth()) {
      return false;
    }

    // get partonic mothers
    auto partonicMother1 = particle1.template fMcPartMoth_as<T3>();
    auto partonicMother2 = particle2.template fMcPartMoth_as<T3>();

    return partonicMother1.globalIndex() == partonicMother2.globalIndex();
  };

  // require both particles to have the same first ancestor, i.e. the same direct mother
  // there is exactly one row in the mother table per generated mother, so comparing the indices is sufficient
  template <typename T1, typename T2>
  bool mcPairHasCommonMother(T1 const& particle1, T2 const& particle2) const
  {
    // if one of the two particles has no associated mother, we cannot know if they have a common anchestor, so we break out with false
    if (!particle1.has_fMcMother() || !particle2.has_fMcMother()) {
      return false;
    }

    return particle1.fMcMotherId() == particle2.fMcMotherId();
  };

  o2::framework::HistogramRegistry* mHistogramRegistry = nullptr;
  bool mUseMotherAsAncestor = false;
};

template <auto& prefix>
class TrackTrackPairCleaner : public BasePairCleaner<prefix>
{
 public:
  TrackTrackPairCleaner() = default;

  template <typename T1, typename T2, typename T3, typename T4>
  bool isCleanPair(T1 const& track1, T2 const& track2, T3 const& /*trackTable*/, T4 const& pairHistManager) const
  {
    this->fillAll();
    bool isClean = this->isCleanParticlePair(track1, track2);
    if (!isClean) {
      this->fillBlocked(pairHistManager.getKinematic());
    }
    return isClean;
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  bool isCleanPair(T1 const& track1, T2 const& track2, T3 const& trackTable, T4 const& mcParticles, T5 const& partonicMothers, T6 const& pairHistManager) const
  {
    if (!this->isCleanPair(track1, track2, trackTable, pairHistManager)) {
      return false;
    }
    // pair is clean
    // no check if we require common or non-common ancestry
    if (this->mMixPairsWithCommonAncestor) {
      return this->pairHasCommonAncestor(track1, track2, mcParticles, partonicMothers);
    }
    if (this->mMixPairsWithNonCommonAncestor) {
      return this->pairHasNonCommonAncestor(track1, track2, mcParticles, partonicMothers);
    }
    return true;
  }
};

template <auto& prefix>
class V0V0PairCleaner : public BasePairCleaner<prefix> // also works for particles decaying into a positive and negative daughter, like resonances
{
 public:
  V0V0PairCleaner() = default;

  template <typename T1, typename T2, typename T3, typename T4>
  bool isCleanPair(T1 const& v01, T2 const& v02, T3 const& trackTable, T4 const& pairHistManager) const
  {
    this->fillAll();
    auto posDaughter1 = trackTable.rawIteratorAt(v01.posDauId() - trackTable.offset());
    auto negDaughter1 = trackTable.rawIteratorAt(v01.negDauId() - trackTable.offset());
    auto posDaughter2 = trackTable.rawIteratorAt(v02.posDauId() - trackTable.offset());
    auto negDaughter2 = trackTable.rawIteratorAt(v02.negDauId() - trackTable.offset());
    bool isClean = this->isCleanParticlePair(v01, v02) &&
                   this->isCleanParticlePair(posDaughter1, posDaughter2) && this->isCleanParticlePair(negDaughter1, negDaughter2) &&
                   this->isCleanParticlePair(posDaughter1, negDaughter2) && this->isCleanParticlePair(negDaughter1, posDaughter2);
    if (!isClean) {
      this->fillBlocked(pairHistManager.getKinematic());
    }
    return isClean;
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  bool isCleanPair(T1 const& v01, T2 const& v02, T3 const& trackTable, T4 const& mcParticles, T5 const& partonicMothers, T6 const& pairHistManager) const
  {
    if (!this->isCleanPair(v01, v02, trackTable, pairHistManager)) {
      return false;
    }
    // pair is clean
    // no check if we require common or non-common ancestry
    if (this->mMixPairsWithCommonAncestor) {
      return this->pairHasCommonAncestor(v01, v02, mcParticles, partonicMothers);
    }
    if (this->mMixPairsWithNonCommonAncestor) {
      return this->pairHasNonCommonAncestor(v01, v02, mcParticles, partonicMothers);
    }
    return true;
  }
};

template <auto& prefix>
class TrackV0PairCleaner : public BasePairCleaner<prefix> // also works for particles decaying into a positive and negative daughter, like resonances
{
 public:
  TrackV0PairCleaner() = default;

  template <typename T1, typename T2, typename T3, typename T4>
  bool isCleanPair(T1 const& track, T2 const& v0, T3 const& trackTable, T4 const& pairHistManager) const
  {
    this->fillAll();
    auto posDaughter = trackTable.rawIteratorAt(v0.posDauId() - trackTable.offset());
    auto negDaughter = trackTable.rawIteratorAt(v0.negDauId() - trackTable.offset());
    bool isClean = this->isCleanParticlePair(posDaughter, track) && this->isCleanParticlePair(negDaughter, track);
    if (!isClean) {
      this->fillBlocked(pairHistManager.getKinematic());
    }
    return isClean;
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  bool isCleanPair(T1 const& track1, T2 const& v0, T3 const& trackTable, T4 const& mcParticles, T5 const& partonicMothers, T6 const& pairHistManager) const
  {
    if (!this->isCleanPair(track1, v0, trackTable, pairHistManager)) {
      return false;
    }
    // pair is clean
    // now check if we require common or non-common ancestry
    if (this->mMixPairsWithCommonAncestor) {
      return this->pairHasCommonAncestor(track1, v0, mcParticles, partonicMothers);
    }
    if (this->mMixPairsWithNonCommonAncestor) {
      return this->pairHasNonCommonAncestor(track1, v0, mcParticles, partonicMothers);
    }
    return true;
  }
};

template <auto& prefix>
class TrackKinkPairCleaner : public BasePairCleaner<prefix>
{
 public:
  TrackKinkPairCleaner() = default;

  template <typename T1, typename T2, typename T3, typename T4>
  bool isCleanPair(T1 const& track, T2 const& kink, T3 const& trackTable, T4 const& pairHistManager) const
  {
    this->fillAll();
    auto chaDaughter = trackTable.rawIteratorAt(kink.chaDauId() - trackTable.offset());
    bool isClean = this->isCleanParticlePair(chaDaughter, track);
    if (!isClean) {
      this->fillBlocked(pairHistManager.getKinematic());
    }
    return isClean;
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  bool isCleanPair(T1 const& track1, T2 const& kink, T3 const& trackTable, T4 const& mcParticles, T5 const& partonicMothers, T6 const& pairHistManager) const
  {
    if (!this->isCleanPair(track1, kink, trackTable, pairHistManager)) {
      return false;
    }
    // pair is clean
    // now check if we require common or non-common ancestry
    if (this->mMixPairsWithCommonAncestor) {
      return this->pairHasCommonAncestor(track1, kink, mcParticles, partonicMothers);
    }
    if (this->mMixPairsWithNonCommonAncestor) {
      return this->pairHasNonCommonAncestor(track1, kink, mcParticles, partonicMothers);
    }
    return true;
  }
};

template <auto& prefix>
class TrackCascadePairCleaner : public BasePairCleaner<prefix>
{
 public:
  TrackCascadePairCleaner() = default;

  template <typename T1, typename T2, typename T3, typename T4>
  bool isCleanPair(T1 const& track, T2 const& cascade, T3 const& trackTable, T4 const& pairHistManager) const
  {
    this->fillAll();
    auto bachelor = trackTable.rawIteratorAt(cascade.bachelorId() - trackTable.offset());
    auto posDaughter = trackTable.rawIteratorAt(cascade.posDauId() - trackTable.offset());
    auto negDaughter = trackTable.rawIteratorAt(cascade.negDauId() - trackTable.offset());
    bool isClean = this->isCleanParticlePair(bachelor, track) && this->isCleanParticlePair(posDaughter, track) && this->isCleanParticlePair(negDaughter, track);
    if (!isClean) {
      this->fillBlocked(pairHistManager.getKinematic());
    }
    return isClean;
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  bool isCleanPair(T1 const& track1, T2 const& cascade, T3 const& trackTable, T4 const& mcParticles, T5 const& partonicMothers, T6 const& pairHistManager) const
  {
    if (!this->isCleanPair(track1, cascade, trackTable, pairHistManager)) {
      return false;
    }
    // pair is clean
    // now check if we require common or non-common ancestry
    if (this->mMixPairsWithCommonAncestor) {
      return this->pairHasCommonAncestor(track1, cascade, mcParticles, partonicMothers);
    }
    if (this->mMixPairsWithNonCommonAncestor) {
      return this->pairHasNonCommonAncestor(track1, cascade, mcParticles, partonicMothers);
    }
    return true;
  }
};

template <auto& prefix>
class TrackLcPairCleaner : public BasePairCleaner<prefix>
{
 public:
  TrackLcPairCleaner() = default;

  template <typename T1, typename T2, typename T3, typename T4>
  bool isCleanPair(T1 const& track, T2 const& lc, T3 const& trackTable, T4 const& pairHistManager) const
  {
    this->fillAll();
    auto prong0 = trackTable.rawIteratorAt(lc.prong0DauId() - trackTable.offset());
    auto prong1 = trackTable.rawIteratorAt(lc.prong1DauId() - trackTable.offset());
    auto prong2 = trackTable.rawIteratorAt(lc.prong2DauId() - trackTable.offset());
    bool isClean = this->isCleanParticlePair(prong0, track) && this->isCleanParticlePair(prong1, track) && this->isCleanParticlePair(prong2, track);
    if (!isClean) {
      this->fillBlocked(pairHistManager.getKinematic());
    }
    return isClean;
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  bool isCleanPair(T1 const& track1, T2 const& lc, T3 const& trackTable, T4 const& mcParticles, T5 const& partonicMothers, T6 const& pairHistManager) const
  {
    if (!this->isCleanPair(track1, lc, trackTable, pairHistManager)) {
      return false;
    }
    // pair is clean
    // now check if we require common or non-common ancestry
    if (this->mMixPairsWithCommonAncestor) {
      return this->pairHasCommonAncestor(track1, lc, mcParticles, partonicMothers);
    }
    if (this->mMixPairsWithNonCommonAncestor) {
      return this->pairHasNonCommonAncestor(track1, lc, mcParticles, partonicMothers);
    }
    return true;
  }
};

template <auto& prefix>
class McParticleMcParticlePairCleaner : public BasePairCleaner<prefix>
{
 public:
  McParticleMcParticlePairCleaner() = default;

  template <typename T1, typename T2, typename T3>
  bool isCleanPair(T1 const& particle1, T2 const& particle2, T3 const& pairHistManager) const
  {
    this->fillAll();
    bool isClean = this->isCleanParticlePair(particle1, particle2);
    if (!isClean) {
      this->fillBlocked(pairHistManager.getKinematic());
    }
    return isClean;
  }

  template <typename T1, typename T2, typename T3, typename T4>
  bool isCleanPair(T1 const& particle1, T2 const& particle2, T3 const& partonicMothers, T4 const& pairHistManager) const
  {
    if (!this->isCleanPair(particle1, particle2, pairHistManager)) {
      return false;
    }
    // pair is clean
    // no check if we require common or non-common ancestry
    if (this->mMixPairsWithCommonAncestor) {
      return this->mcPairHasCommonAncestor(particle1, particle2, partonicMothers);
    }
    if (this->mMixPairsWithNonCommonAncestor) {
      return this->mcPairHasNonCommonAncestor(particle1, particle2, partonicMothers);
    }
    return true;
  }
};

} // namespace o2::analysis::femto::paircleaner

#endif // PWGCF_FEMTO_CORE_PAIRCLEANER_H_
