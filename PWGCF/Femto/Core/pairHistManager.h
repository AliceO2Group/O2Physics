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

/// \file pairHistManager.h
/// \brief histogram manager for pair tasks
/// \author anton.riedel@tum.de, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTO_CORE_PAIRHISTMANAGER_H_
#define PWGCF_FEMTO_CORE_PAIRHISTMANAGER_H_

#include "PWGCF/Femto/Core/femtoUtils.h"
#include "PWGCF/Femto/Core/histManager.h"
#include "PWGCF/Femto/Core/modes.h"

#include "Framework/Configurable.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/HistogramSpec.h"

#include <Math/GenVector/Boost.h>
#include <Math/Vector4D.h>

#include <array>
#include <cmath>
#include <map>
#include <string>
#include <string_view>
#include <vector>

namespace o2::analysis::femto
{
namespace pairhistmanager
{
// enum for pair histograms
enum PairHist {
  // kinemtics
  kKstar,
  kKt,
  kMt,
  // 2d
  kPt1VsPt2,
  kPt1VsKstar,
  kPt2VsKstar,
  kPt1VsKt,
  kPt2VsKt,
  kPt1VsMt,
  kPt2VsMt,
  kPairHistogramLast
  //  more dimensions
};

enum MixingPoliciy {
  kVtxMult,
  kVtxCent,
  kVtxMultCent,
  kMixingPolicyLast
};

// Mixing configurables
struct ConfMixing : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("Mixing");
  o2::framework::ConfigurableAxis multBins{"multBins", {o2::framework::VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f}, "Mixing bins - multiplicity"};
  o2::framework::ConfigurableAxis centBins{"centBins", {o2::framework::VARIABLE_WIDTH, 0.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f}, "Mixing bins - centrality"};
  o2::framework::ConfigurableAxis vtxBins{"vtxBins", {o2::framework::VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  o2::framework::Configurable<int> depth{"depth", 5, "Number of events for mixing"};
  o2::framework::Configurable<int> policy{"policy", 0, "Binning policy for mixing (alywas in combination with z-vertex) -> 0: multiplicity, -> 1: centrality, -> 2: both"};
  o2::framework::Configurable<bool> sameSpecies{"sameSpecies", false, "Enable if partilce 1 and particle 2 are the same"};
  o2::framework::Configurable<int> seed{"seed", -1, "Seed to randomize particle 1 and particle 2 (if they are identical). Set to negative value to deactivate. Set to 0 to generate unique seed in time."};
};

struct ConfPairBinning : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("PairBinning");
  o2::framework::ConfigurableAxis kstar{"kstar", {{600, 0, 6}}, "kstar"};
  o2::framework::ConfigurableAxis kt{"kt", {{600, 0, 6}}, "kt"};
  o2::framework::ConfigurableAxis mt{"mt", {{500, 0.8, 5.8}}, "mt"};
};

// the enum gives the correct index in the array
constexpr std::array<histmanager::HistInfo<PairHist>, kPairHistogramLast> HistTable = {
  {{kKstar, o2::framework::kTH1F, "hKstar", "k*; k* (GeV/#it{c}); Entries"},
   {kKt, o2::framework::kTH1F, "hKt", "transverse momentum; k_{T} (GeV/#it{c}); Entries"},
   {kMt, o2::framework::kTH1F, "hMt", "transverse mass; m_{T} (GeV/#it{c}^{2}); Entries"},
   {kPt1VsPt2, o2::framework::kTH2F, "hPt1VsPt2", "track1 p_{T} vs track2 p_{T}; track1 p_T (GeV/#it{c}); track2 p_{T} (GeV/#it{c})"},
   {kPt1VsKstar, o2::framework::kTH2F, "hPt1VsKstar", "p_{T,1} vs k*; p_{T,2} (GeV/#it{c}); k* (GeV/#it{c})"},
   {kPt2VsKstar, o2::framework::kTH2F, "hPt2VsKstar", "p_{T,2} vs k*; p_{T,2} (GeV/#it{c}); k* (GeV/#it{c})"},
   {kPt1VsKt, o2::framework::kTH2F, "hPt1VsKt", "p_{T,1} vs k_{T}; p_{T,1} (GeV/#it{c}); k_{T} (GeV/#it{c})"},
   {kPt2VsKt, o2::framework::kTH2F, "hPt2VsKt", "p_{T,2} vs k_{T}; p_{T,2} (GeV/#it{c}); k_{T} (GeV/#it{c})"},
   {kPt1VsMt, o2::framework::kTH2F, "hPt1VsMt", "p_{T,1} vs m_{T}; p_{T,1} (GeV/#it{c}); m_{T} (GeV/#it{c})"},
   {kPt2VsMt, o2::framework::kTH2F, "hPt2VsMt", "p_{T,2} vs m_{T}; p_{T,2} (GeV/#it{c}); m_{T} (GeV/#it{c})"}}};

template <typename T1, typename T2, typename T3>
auto makePairHistSpecMap(const T1& confPairBinning, const T2& confObject1Binning, const T3& confObject2Binning)
{
  return std::map<PairHist, std::vector<framework::AxisSpec>>{
    {kKstar, {confPairBinning.kstar}},
    {kKt, {confPairBinning.kt}},
    {kMt, {confPairBinning.mt}},
    {kPt1VsPt2, {confObject1Binning.pt, confObject2Binning.pt}},
    {kPt1VsKstar, {confObject1Binning.pt, confPairBinning.kstar}},
    {kPt2VsKstar, {confObject2Binning.pt, confPairBinning.kstar}},
    {kPt1VsKt, {confObject1Binning.pt, confPairBinning.kt}},
    {kPt2VsKt, {confObject2Binning.pt, confPairBinning.kt}},
    {kPt1VsMt, {confObject1Binning.pt, confPairBinning.mt}},
    {kPt2VsMt, {confObject2Binning.pt, confPairBinning.mt}}};
};

constexpr char PrefixTrackTrackSe[] = "TrackTrack/SE/";
constexpr char PrefixTrackTrackMe[] = "TrackTrack/ME/";

constexpr char PrefixTrackV0Se[] = "TrackV0/SE/";
constexpr char PrefixTrackV0Me[] = "TrackV0/ME/";

constexpr char PrefixV0V0Se[] = "V0V0/SE/";
constexpr char PrefixV0V0Me[] = "V0V0/ME/";

constexpr char PrefixTrackResonanceSe[] = "TrackResonance/SE/";
constexpr char PrefixTrackResonanceMe[] = "TrackResonance/ME/";

constexpr char PrefixTrackCascadeSe[] = "TrackCascade/SE/";
constexpr char PrefixTrackCascadeMe[] = "TrackCascade/ME/";

constexpr char PrefixTrackKinkSe[] = "TrackKink/SE/";
constexpr char PrefixTrackKinkMe[] = "TrackKink/ME/";

constexpr std::string_view AnalysisDir = "Analysis/";
constexpr std::string_view QaDir = "QA/";

/// \class FemtoDreamEventHisto
/// \brief Class for histogramming event properties
// template <femtomodes::Mode mode>
template <const char* prefix, modes::Mode mode>
class PairHistManager
{
 public:
  /// Destructor
  virtual ~PairHistManager() = default;
  void init(o2::framework::HistogramRegistry* registry, std::map<PairHist, std::vector<o2::framework::AxisSpec>> Specs)
  {
    mHistogramRegistry = registry;

    if constexpr (isFlagSet(mode, modes::Mode::kAnalysis)) {
      std::string analysisDir = std::string(prefix) + std::string(AnalysisDir);

      mHistogramRegistry->add(analysisDir + getHistNameV2(kKstar, HistTable), getHistDesc(kKstar, HistTable), getHistType(kKstar, HistTable), {Specs[kKstar]});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kMt, HistTable), getHistDesc(kMt, HistTable), getHistType(kMt, HistTable), {Specs[kMt]});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kPt1VsPt2, HistTable), getHistDesc(kPt1VsPt2, HistTable), getHistType(kPt1VsPt2, HistTable), {Specs[kPt1VsPt2]});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kPt1VsKstar, HistTable), getHistDesc(kPt1VsKstar, HistTable), getHistType(kPt1VsKstar, HistTable), {Specs[kPt1VsKstar]});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kPt2VsKstar, HistTable), getHistDesc(kPt2VsKstar, HistTable), getHistType(kPt2VsKstar, HistTable), {Specs[kPt2VsKstar]});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kPt1VsKt, HistTable), getHistDesc(kPt1VsKt, HistTable), getHistType(kPt1VsKt, HistTable), {Specs[kPt1VsKt]});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kPt2VsKt, HistTable), getHistDesc(kPt2VsKt, HistTable), getHistType(kPt2VsKt, HistTable), {Specs[kPt2VsKt]});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kPt1VsMt, HistTable), getHistDesc(kPt1VsMt, HistTable), getHistType(kPt1VsMt, HistTable), {Specs[kPt1VsMt]});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kPt2VsMt, HistTable), getHistDesc(kPt2VsMt, HistTable), getHistType(kPt2VsMt, HistTable), {Specs[kPt2VsMt]});
    }

    // if constexpr (isFlagSet(mode, modes::Mode::kQA)) {
    // std::string qaDir = std::string(prefix) + std::string(QaDir);
    // }
  }

  void setMass(int PdgParticle1, int PdgParticle2)
  {
    mMass1 = o2::analysis::femto::utils::getMass(PdgParticle1);
    mMass2 = o2::analysis::femto::utils::getMass(PdgParticle2);
  }
  void setCharge(int chargeAbsParticle1, int chargeAbsParticle2)
  {
    // the pt stored is actually pt/z, so in case of particles with z > 1, we have to rescale the pt (this is so far only for He3 the case)
    // similarly, for neutral particles, no reason to rescale so we just set absolute charge to 1
    mAbsCharge1 = chargeAbsParticle1;
    mAbsCharge2 = chargeAbsParticle2;
  }

  template <typename T1, typename T2>
  void setPair(const T1& particle1, const T2& particle2)
  {
    // pt in track stable is stored from 1/signedPt from the original track table
    // in case of He with Z=2, we have to rescale the pt with the absolute charge
    mParticle1 = ROOT::Math::PtEtaPhiMVector{mAbsCharge1 * particle1.pt(), particle1.eta(), particle1.phi(), mMass1};
    mParticle2 = ROOT::Math::PtEtaPhiMVector{mAbsCharge2 * particle2.pt(), particle2.eta(), particle2.phi(), mMass2};
    auto partSum = mParticle1 + mParticle2;

    // set kT
    mKt = partSum.Pt() / 2.f;

    // set mT
    float averageMass = (mMass1 + mMass2) / 2.f;
    mMt = std::hypot(mKt, averageMass);

    // Boost Track1 to the pair rest frame and calculate k*
    auto track1 = ROOT::Math::PxPyPzEVector(mParticle1);
    ROOT::Math::Boost boostPrf(partSum.BoostToCM());
    mKstar = boostPrf(track1).P();
  }

  void fill()
  {
    if constexpr (isFlagSet(mode, modes::Mode::kAnalysis)) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstar, HistTable)), mKstar);
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kMt, HistTable)), mMt);
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kPt1VsPt2, HistTable)), mParticle1.Pt(), mParticle2.Pt());
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kPt1VsKstar, HistTable)), mParticle1.Pt(), mKstar);
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kPt1VsMt, HistTable)), mParticle1.Pt(), mMt);
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kPt1VsKt, HistTable)), mParticle1.Pt(), mKt);
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kPt2VsKstar, HistTable)), mParticle2.Pt(), mKstar);
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kPt2VsMt, HistTable)), mParticle2.Pt(), mMt);
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kPt2VsKt, HistTable)), mParticle2.Pt(), mKt);
    }

    // if constexpr (isFlagSet(mode, modes::Mode::kQA)) {
    //  mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(getHistName(kPtVsDcaz, HistTable)), track.pt(), track.dcaZ());
    // }
  }

 private:
  o2::framework::HistogramRegistry* mHistogramRegistry = nullptr;
  float mMass1 = 0.f;
  float mMass2 = 0.f;
  float mAbsCharge1 = 1.f;
  float mAbsCharge2 = 1.f;
  ROOT::Math::PtEtaPhiMVector mParticle1{};
  ROOT::Math::PtEtaPhiMVector mParticle2{};
  float mKstar = 0.f;
  float mKt = 0.f;
  float mMt = 0.f;
};

}; // namespace pairhistmanager
}; // namespace o2::analysis::femto
#endif // PWGCF_FEMTO_CORE_PAIRHISTMANAGER_H_
