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
  kKstarVsKt,
  kKstarVsMt,
  kKstarVsMult,
  kKstarVsCent,

  //  higher dimensions
  kKstarVsMtVsMult,
  kKstarVsMtVsMultVsCent,
  kKstarVsMtVsPt1VsP2VsMult,
  kKstarVsMtVsPt1VsP2VsMultVsCent,

  kPairHistogramLast
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
  o2::framework::Configurable<bool> plot1D{"plot1D", true, "Enable 1D histograms"};
  o2::framework::Configurable<bool> plot2D{"plot2D", true, "Enable 2D histograms"};
  o2::framework::Configurable<bool> plotKstarVsMtVsMult{"plotKstarVsMtVsMult", false, "Enable 3D histogram (Kstar Vs Mt Vs Mult)"};
  o2::framework::Configurable<bool> plotKstarVsMtVsMultVsCent{"plotKstarVsMtVsMultVsCent", false, "Enable 4D histogram (Kstar Vs Mt Vs Mult Vs Cent)"};
  o2::framework::Configurable<bool> plotKstarVsMtVsPt1VsP2VsMult{"plotKstarVsMtVsPt1VsP2VsMult", false, "Enable 5D histogram (Kstar Vs Mt Vs Pt1 Vs Pt2 Vs Mult)"};
  o2::framework::Configurable<bool> plotKstarVsMtVsPt1VsP2VsMultVsCent{"plotKstarVsMtVsPt1VsP2VsMultVsCent", false, "Enable 6D histogram (Kstar Vs Mt Vs Pt1 Vs Pt2 Vs Mult Vs Cent)"};
  o2::framework::ConfigurableAxis kstar{"kstar", {{600, 0, 6}}, "kstar"};
  o2::framework::ConfigurableAxis kt{"kt", {{600, 0, 6}}, "kt"};
  o2::framework::ConfigurableAxis mt{"mt", {{500, 0.8, 5.8}}, "mt"};
  o2::framework::ConfigurableAxis multiplicity{"multiplicity", {{50, 0, 200}}, "multiplicity"};
  o2::framework::ConfigurableAxis centrality{"centrality", {{10, 0, 100}}, "centrality (mult. percentile)"};
};

// the enum gives the correct index in the array
constexpr std::array<histmanager::HistInfo<PairHist>, kPairHistogramLast> HistTable = {
  {
    // 1D
    {kKstar, o2::framework::kTH1F, "hKstar", "k*; k* (GeV/#it{c}); Entries"},
    {kKt, o2::framework::kTH1F, "hKt", "transverse momentum; k_{T} (GeV/#it{c}); Entries"},
    {kMt, o2::framework::kTH1F, "hMt", "transverse mass; m_{T} (GeV/#it{c}^{2}); Entries"},
    // 2D
    {kPt1VsPt2, o2::framework::kTH2F, "hPt1VsPt2", "track1 p_{T} vs track2 p_{T}; track1 p_T (GeV/#it{c}); track2 p_{T} (GeV/#it{c})"},
    {kPt1VsKstar, o2::framework::kTH2F, "hPt1VsKstar", "p_{T,1} vs k*; p_{T,2} (GeV/#it{c}); k* (GeV/#it{c})"},
    {kPt2VsKstar, o2::framework::kTH2F, "hPt2VsKstar", "p_{T,2} vs k*; p_{T,2} (GeV/#it{c}); k* (GeV/#it{c})"},
    {kPt1VsKt, o2::framework::kTH2F, "hPt1VsKt", "p_{T,1} vs k_{T}; p_{T,1} (GeV/#it{c}); k_{T} (GeV/#it{c})"},
    {kPt2VsKt, o2::framework::kTH2F, "hPt2VsKt", "p_{T,2} vs k_{T}; p_{T,2} (GeV/#it{c}); k_{T} (GeV/#it{c})"},
    {kPt1VsMt, o2::framework::kTH2F, "hPt1VsMt", "p_{T,1} vs m_{T}; p_{T,1} (GeV/#it{c}); m_{T} (GeV/#it{c}^{2})"},
    {kPt2VsMt, o2::framework::kTH2F, "hPt2VsMt", "p_{T,2} vs m_{T}; p_{T,2} (GeV/#it{c}); m_{T} (GeV/#it{c}^{2})"},
    {kKstarVsKt, o2::framework::kTH2F, "hKstarVsKt", "k* vs k_{T}; k* (GeV/#it{c}); k_{T} (GeV/#it{c})"},
    {kKstarVsMt, o2::framework::kTH2F, "hKstarVsMt", "k* vs m_{T}; k* (GeV/#it{c}); m_{T} (GeV/#it{c}^{2})"},
    {kKstarVsCent, o2::framework::kTH2F, "hKstarVsCent", "k* vs Centrality (Mult. Percentile); k* (GeV/#it{c}); Centrality (%)"},
    {kKstarVsMult, o2::framework::kTH2F, "hKstarVsMult", "k* vs Multiplicity; k* (GeV/#it{c}); Multiplicity"},

    // n-D
    {kKstarVsMtVsMult, o2::framework::kTHnSparseF, "hKstarVsMtVsMult", "k* vs m_{T} vs multiplicity; k* (GeV/#it{c}); m_{T} (GeV/#it{c}^{2}); Multiplicity"},
    {kKstarVsMtVsMultVsCent, o2::framework::kTHnSparseF, "hKstarVsMtVsMultVsCent", "k* vs m_{T} vs multiplicity vs centrality; k* (GeV/#it{c}); m_{T} (GeV/#it{c}^{2}); Multiplicity; Centrality (%)"},
    {kKstarVsMtVsPt1VsP2VsMult, o2::framework::kTHnSparseF, "hKstarVsMtVsPt1VsPt2VsMult", "k* vs m_{T} vs p_{T,1} vs p_{T,2} vs multiplicity; k* (GeV/#it{c}); m_{T} (GeV/#it{c}^{2}); p_{T,1} (GeV/#it{c}); p_{T,2} (GeV/#it{c}); Multiplicity"},
    {kKstarVsMtVsPt1VsP2VsMultVsCent, o2::framework::kTHnSparseF, "hKstarVsMtVsPt1VsPt2VsMultVsCent", "k* vs m_{T} vs p_{T,1} vs p_{T,2} vs multiplicity vs centrality; k* (GeV/#it{c}); m_{T} (GeV/#it{c}^{2}); p_{T,1} (GeV/#it{c}); p_{T,2} (GeV/#it{c}); Multiplicity; Centrality"},

  }};

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
    {kPt2VsMt, {confObject2Binning.pt, confPairBinning.mt}},
    {kKstarVsKt, {confPairBinning.kstar, confPairBinning.kt}},
    {kKstarVsMt, {confPairBinning.kstar, confPairBinning.mt}},
    {kKstarVsMult, {confPairBinning.kstar, confPairBinning.multiplicity}},
    {kKstarVsCent, {confPairBinning.kstar, confPairBinning.centrality}},

    {kKstarVsMtVsMult, {confPairBinning.kstar, confPairBinning.mt, confPairBinning.multiplicity}},
    {kKstarVsMtVsMultVsCent, {confPairBinning.kstar, confPairBinning.mt, confPairBinning.multiplicity, confPairBinning.centrality}},
    {kKstarVsMtVsPt1VsP2VsMult, {confPairBinning.kstar, confPairBinning.mt, confObject1Binning.pt, confObject2Binning.pt, confPairBinning.multiplicity}},
    {kKstarVsMtVsPt1VsP2VsMultVsCent, {confPairBinning.kstar, confPairBinning.mt, confObject1Binning.pt, confObject2Binning.pt, confPairBinning.multiplicity, confPairBinning.centrality}},
  };
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
  PairHistManager() = default;
  ~PairHistManager() = default;

  template <typename T>
  void init(o2::framework::HistogramRegistry* registry, std::map<PairHist, std::vector<o2::framework::AxisSpec>> const& Specs, T const& ConfPairBinning)
  {
    mHistogramRegistry = registry;

    mPlot1d = ConfPairBinning.plot1D.value;
    mPlot2d = ConfPairBinning.plot2D.value;
    mPlotKstarVsMtVsMult = ConfPairBinning.plotKstarVsMtVsMult.value;
    mPlotKstarVsMtVsMultVsCent = ConfPairBinning.plotKstarVsMtVsMultVsCent.value;
    mPlotKstarVsMtVsPt1VsP2VsMult = ConfPairBinning.plotKstarVsMtVsPt1VsP2VsMult.value;
    mPlotKstarVsMtVsPt1VsP2VsMultVsCent = ConfPairBinning.plotKstarVsMtVsPt1VsP2VsMultVsCent.value;

    if constexpr (isFlagSet(mode, modes::Mode::kAnalysis)) {
      std::string analysisDir = std::string(prefix) + std::string(AnalysisDir);
      if (mPlot1d) {
        mHistogramRegistry->add(analysisDir + getHistNameV2(kKstar, HistTable), getHistDesc(kKstar, HistTable), getHistType(kKstar, HistTable), {Specs.at(kKstar)});
        mHistogramRegistry->add(analysisDir + getHistNameV2(kKt, HistTable), getHistDesc(kKt, HistTable), getHistType(kKt, HistTable), {Specs.at(kKt)});
        mHistogramRegistry->add(analysisDir + getHistNameV2(kMt, HistTable), getHistDesc(kMt, HistTable), getHistType(kMt, HistTable), {Specs.at(kMt)});
      }
      if (mPlot2d) {
        mHistogramRegistry->add(analysisDir + getHistNameV2(kPt1VsPt2, HistTable), getHistDesc(kPt1VsPt2, HistTable), getHistType(kPt1VsPt2, HistTable), {Specs.at(kPt1VsPt2)});
        mHistogramRegistry->add(analysisDir + getHistNameV2(kPt1VsKstar, HistTable), getHistDesc(kPt1VsKstar, HistTable), getHistType(kPt1VsKstar, HistTable), {Specs.at(kPt1VsKstar)});
        mHistogramRegistry->add(analysisDir + getHistNameV2(kPt2VsKstar, HistTable), getHistDesc(kPt2VsKstar, HistTable), getHistType(kPt2VsKstar, HistTable), {Specs.at(kPt2VsKstar)});
        mHistogramRegistry->add(analysisDir + getHistNameV2(kPt1VsKt, HistTable), getHistDesc(kPt1VsKt, HistTable), getHistType(kPt1VsKt, HistTable), {Specs.at(kPt1VsKt)});
        mHistogramRegistry->add(analysisDir + getHistNameV2(kPt2VsKt, HistTable), getHistDesc(kPt2VsKt, HistTable), getHistType(kPt2VsKt, HistTable), {Specs.at(kPt2VsKt)});
        mHistogramRegistry->add(analysisDir + getHistNameV2(kPt1VsMt, HistTable), getHistDesc(kPt1VsMt, HistTable), getHistType(kPt1VsMt, HistTable), {Specs.at(kPt1VsMt)});
        mHistogramRegistry->add(analysisDir + getHistNameV2(kPt2VsMt, HistTable), getHistDesc(kPt2VsMt, HistTable), getHistType(kPt2VsMt, HistTable), {Specs.at(kPt2VsMt)});

        mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsKt, HistTable), getHistDesc(kKstarVsKt, HistTable), getHistType(kKstarVsKt, HistTable), {Specs.at(kKstarVsKt)});
        mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsMt, HistTable), getHistDesc(kKstarVsMt, HistTable), getHistType(kKstarVsMt, HistTable), {Specs.at(kKstarVsMt)});
        mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsMult, HistTable), getHistDesc(kKstarVsMult, HistTable), getHistType(kKstarVsMult, HistTable), {Specs.at(kKstarVsMult)});
        mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsCent, HistTable), getHistDesc(kKstarVsCent, HistTable), getHistType(kKstarVsCent, HistTable), {Specs.at(kKstarVsCent)});
      }

      if (mPlotKstarVsMtVsMult) {
        mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsMtVsMult, HistTable), getHistDesc(kKstarVsMtVsMult, HistTable), getHistType(kKstarVsMtVsMult, HistTable), {Specs.at(kKstarVsMtVsMult)});
      }
      if (mPlotKstarVsMtVsMultVsCent) {
        mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsMtVsMultVsCent, HistTable), getHistDesc(kKstarVsMtVsMultVsCent, HistTable), getHistType(kKstarVsMtVsMultVsCent, HistTable), {Specs.at(kKstarVsMtVsMultVsCent)});
      }
      if (mPlotKstarVsMtVsPt1VsP2VsMult) {
        mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsMtVsPt1VsP2VsMult, HistTable), getHistDesc(kKstarVsMtVsPt1VsP2VsMult, HistTable), getHistType(kKstarVsMtVsPt1VsP2VsMult, HistTable), {Specs.at(kKstarVsMtVsPt1VsP2VsMult)});
      }
      if (mPlotKstarVsMtVsPt1VsP2VsMultVsCent) {
        mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsMtVsPt1VsP2VsMultVsCent, HistTable), getHistDesc(kKstarVsMtVsPt1VsP2VsMultVsCent, HistTable), getHistType(kKstarVsMtVsPt1VsP2VsMultVsCent, HistTable), {Specs.at(kKstarVsMtVsPt1VsP2VsMultVsCent)});
      }
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
    // the pt stored is actually as pt/z for tracks, so in case of particles with z > 1, we have to rescale the pt (this is so far only for He3 the case)
    // similarly, for neutral particles, no reason to rescale so we just set absolute charge to 1
    mAbsCharge1 = std::abs(chargeAbsParticle1);
    mAbsCharge2 = std::abs(chargeAbsParticle2);
  }

  template <typename T1, typename T2, typename T3>
  void setPair(const T1& particle1, const T2& particle2, const T3& col)
  {
    // set collision properties
    mMult = col.mult();
    mCent = col.cent();

    // pt in track table is calculated from 1/signedPt from the original track table
    // in case of He with Z=2, we have to rescale the pt with the absolute charge
    mParticle1 = ROOT::Math::PtEtaPhiMVector{mAbsCharge1 * particle1.pt(), particle1.eta(), particle1.phi(), mMass1};
    mParticle2 = ROOT::Math::PtEtaPhiMVector{mAbsCharge2 * particle2.pt(), particle2.eta(), particle2.phi(), mMass2};
    auto partSum = mParticle1 + mParticle2;

    // set kT
    mKt = partSum.Pt() / 2.f;

    // set mT
    float averageMass = (mMass1 + mMass2) / 2.f;
    mMt = std::hypot(mKt, averageMass);

    // Boost particle to the pair rest frame (Prf) and calculate k* (would be equivalent using particle 2)
    // make a copy of particle 1
    auto particle1Prf = ROOT::Math::PtEtaPhiMVector(mParticle1);
    // get lorentz boost into pair rest frame
    ROOT::Math::Boost boostPrf(partSum.BoostToCM());
    // boost particle 1 into pair rest frame and calculate its momentum, which has the same value as k*
    mKstar = boostPrf(particle1Prf).P();
  }

  void fill()
  {
    if constexpr (isFlagSet(mode, modes::Mode::kAnalysis)) {
      if (mPlot1d) {
        mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstar, HistTable)), mKstar);
        mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kMt, HistTable)), mMt);
        mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKt, HistTable)), mKt);
      }
      if (mPlot2d) {
        mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kPt1VsPt2, HistTable)), mParticle1.Pt(), mParticle2.Pt());
        mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kPt1VsKstar, HistTable)), mParticle1.Pt(), mKstar);
        mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kPt1VsMt, HistTable)), mParticle1.Pt(), mMt);
        mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kPt1VsKt, HistTable)), mParticle1.Pt(), mKt);
        mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kPt2VsKstar, HistTable)), mParticle2.Pt(), mKstar);
        mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kPt2VsMt, HistTable)), mParticle2.Pt(), mMt);
        mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kPt2VsKt, HistTable)), mParticle2.Pt(), mKt);
        mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsKt, HistTable)), mKstar, mKt);
        mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsMt, HistTable)), mKstar, mMt);
        mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsMult, HistTable)), mKstar, mMult);
        mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsCent, HistTable)), mKstar, mCent);
      }
      if (mPlotKstarVsMtVsMult) {
        mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsMtVsMult, HistTable)), mKstar, mMt, mMult);
      }
      if (mPlotKstarVsMtVsMultVsCent) {
        mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsMtVsMultVsCent, HistTable)), mKstar, mMt, mMult, mCent);
      }
      if (mPlotKstarVsMtVsPt1VsP2VsMult) {
        mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsMtVsPt1VsP2VsMult, HistTable)), mKstar, mMt, mParticle1.Pt(), mParticle2.pt(), mMult);
      }
      if (mPlotKstarVsMtVsPt1VsP2VsMultVsCent) {
        mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsMtVsPt1VsP2VsMultVsCent, HistTable)), mKstar, mMt, mParticle1.Pt(), mParticle2.pt(), mMult, mCent);
      }
    }

    // if constexpr (isFlagSet(mode, modes::Mode::kQA)) {
    //  mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(getHistName(kPtVsDcaz, HistTable)), track.pt(), track.dcaZ());
    // }
  }

 private:
  o2::framework::HistogramRegistry* mHistogramRegistry = nullptr;
  float mMass1 = 0.f;
  float mMass2 = 0.f;
  int mAbsCharge1 = 1;
  int mAbsCharge2 = 1;
  ROOT::Math::PtEtaPhiMVector mParticle1{};
  ROOT::Math::PtEtaPhiMVector mParticle2{};
  float mKstar = 0.f;
  float mKt = 0.f;
  float mMt = 0.f;
  float mMult = 0.f;
  float mCent = 0.f;

  // flags
  bool mPlot1d = true;
  bool mPlot2d = true;
  bool mPlotKstarVsMtVsMult = false;
  bool mPlotKstarVsMtVsMultVsCent = false;
  bool mPlotKstarVsMtVsPt1VsP2VsMult = false;
  bool mPlotKstarVsMtVsPt1VsP2VsMultVsCent = false;
};

}; // namespace pairhistmanager
}; // namespace o2::analysis::femto
#endif // PWGCF_FEMTO_CORE_PAIRHISTMANAGER_H_
