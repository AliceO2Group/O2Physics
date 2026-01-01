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
  // standard 1D
  kKstar,
  kKt,
  kMt,
  // standard 2D
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
  // 2D with mass
  kKstarVsMass1,
  kKstarVsMass2,
  kMass1VsMass2,
  // higher dimensions
  kKstarVsMtVsMult,
  kKstarVsMtVsMultVsCent,
  kKstarVsMtVsPt1VsPt2VsMult,
  kKstarVsMtVsPt1VsPt2VsMultVsCent,
  // higher dimensions with mass
  kKstarVsMass1VsMass2,
  kKstarVsMass1VsMult,
  kKstarVsMass2VsMult,
  kKstarVsMass1VsMass2VsMult,
  // mc
  kTrueKstarVsKstar,

  kPairHistogramLast
};

enum MixingPolicy {
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
  o2::framework::Configurable<bool> sameSpecies{"sameSpecies", false, "Enable if particle 1 and particle 2 are the same"};
  o2::framework::Configurable<int> seed{"seed", -1, "Seed to randomize particle 1 and particle 2 (if they are identical). Set to negative value to deactivate. Set to 0 to generate unique seed in time."};
};

struct ConfPairBinning : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("PairBinning");
  o2::framework::Configurable<bool> plot1D{"plot1D", true, "Enable 1D histograms"};
  o2::framework::Configurable<bool> plot2D{"plot2D", true, "Enable 2D histograms"};
  o2::framework::Configurable<bool> plotKstarVsMtVsMult{"plotKstarVsMtVsMult", false, "Enable 3D histogram (Kstar Vs Mt Vs Mult)"};
  o2::framework::Configurable<bool> plotKstarVsMtVsMultVsCent{"plotKstarVsMtVsMultVsCent", false, "Enable 4D histogram (Kstar Vs Mt Vs Mult Vs Cent)"};
  o2::framework::Configurable<bool> plotKstarVsMtVsPt1VsPt2VsMult{"plotKstarVsMtVsPt1VsPt2VsMult", false, "Enable 5D histogram (Kstar Vs Mt Vs Pt1 Vs Pt2 Vs Mult)"};
  o2::framework::Configurable<bool> plotKstarVsMtVsPt1VsPt2VsMultVsCent{"plotKstarVsMtVsPt1VsPt2VsMultVsCent", false, "Enable 6D histogram (Kstar Vs Mt Vs Pt1 Vs Pt2 Vs Mult Vs Cent)"};
  o2::framework::Configurable<bool> plotKstarVsMass1VsMass2{"plotKstarVsMass1VsMass2", false, "Enable 3D histogram (Kstar Vs Mass1 Vs Mass2)"};
  o2::framework::Configurable<bool> plotKstarVsMass1VsMult{"plotKstarVsMass1VsMult", false, "Enable 3D histogram (Kstar Vs Mass1 Vs Mult)"};
  o2::framework::Configurable<bool> plotKstarVsMass2VsMult{"plotKstarVsMass2VsMult", false, "Enable 3D histogram (Kstar Vs Mass2 Vs Mult)"};
  o2::framework::Configurable<bool> plotKstarVsMass1VsMass2VsMult{"plotKstarVsMass1VsMass2VsMult", false, "Enable 4D histogram (Kstar Vs Mass1 Vs Mass2 Vs Mult)"};
  o2::framework::ConfigurableAxis kstar{"kstar", {{600, 0, 6}}, "kstar"};
  o2::framework::ConfigurableAxis kt{"kt", {{600, 0, 6}}, "kt"};
  o2::framework::ConfigurableAxis mt{"mt", {{500, 0.8, 5.8}}, "mt"};
  o2::framework::ConfigurableAxis multiplicity{"multiplicity", {{50, 0, 200}}, "multiplicity"};
  o2::framework::ConfigurableAxis centrality{"centrality", {{10, 0, 100}}, "centrality (mult. percentile)"};
  o2::framework::ConfigurableAxis pt1{"pt1", {{100, 0, 6}}, "Pt binning for particle 1"};
  o2::framework::ConfigurableAxis pt2{"pt2", {{100, 0, 6}}, "Pt binning for particle 2"};
  o2::framework::ConfigurableAxis mass1{"mass1", {{100, 0, 2}}, "Mass binning for particle 1 (if particle has mass getter)"};
  o2::framework::ConfigurableAxis mass2{"mass2", {{100, 0, 2}}, "Mass binning for particle 2 (if particle has mass getter)"};
  o2::framework::Configurable<int> transverseMassType{"transverseMassType", static_cast<int>(modes::TransverseMassType::kAveragePdgMass), "Type of transverse mass (0-> Average Pdg Mass, 1-> Reduced Pdg Mass, 2-> Mt from combined 4 vector)"};
};

struct ConfPairCuts : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("PairCuts");
  o2::framework::Configurable<float> kstarMax{"kstarMax", -1, "Maximal kstar (set to -1 to deactivate)"};
  o2::framework::Configurable<float> kstarMin{"kstarMin", -1, "Minimal kstar (set to -1 to deactivate)"};
  o2::framework::Configurable<float> ktMax{"ktMax", -1, "Maximal kt (set to -1 to deactivate)"};
  o2::framework::Configurable<float> ktMin{"ktMin", -1, "Minimal kt (set to -1 to deactivate)"};
  o2::framework::Configurable<float> mtMax{"mtMax", -1, "Maximal mt (set to -1 to deactivate)"};
  o2::framework::Configurable<float> mtMin{"mtMin", -1, "Minimal mt (set to -1 to deactivate)"};
  o2::framework::Configurable<bool> commonAnchestor{"commonAnchestor", false, "Pair has common anchestor"};
  o2::framework::Configurable<bool> nonCommonAnchestor{"nonCommonAnchestor", false, "Pair has non-common anchestor"};
};

// the enum gives the correct index in the array
constexpr std::array<histmanager::HistInfo<PairHist>, kPairHistogramLast>
  HistTable = {
    {
      // 1D
      {kKstar, o2::framework::kTH1F, "hKstar", "k*; k* (GeV/#it{c}); Entries"},
      {kKt, o2::framework::kTH1F, "hKt", "transverse momentum; k_{T} (GeV/#it{c}); Entries"},
      {kMt, o2::framework::kTH1F, "hMt", "transverse mass; m_{T} (GeV/#it{c}^{2}); Entries"},
      // 2D
      {kPt1VsPt2, o2::framework::kTH2F, "hPt1VsPt2", " p_{T,1} vs p_{T,2}; p_{T,1} (GeV/#it{c}); p_{T,2} (GeV/#it{c})"},
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
      // 2D with mass
      {kKstarVsMass1, o2::framework::kTH2F, "hKstarVsMass1", "k* vs m_{1}; k* (GeV/#it{c}); m_{1} (GeV/#it{c}^{2})"},
      {kKstarVsMass2, o2::framework::kTH2F, "hKstarVsMass2", "k* vs m_{2}; k* (GeV/#it{c}); m_{2} (GeV/#it{c}^{2})"},
      {kMass1VsMass2, o2::framework::kTH2F, "hMass1VsMass2", "m_{1} vs m_{2}; m_{1} (GeV/#it{c}^{2}); m_{2} (GeV/#it{c}^{2})"},
      // n-D
      {kKstarVsMtVsMult, o2::framework::kTHnSparseF, "hKstarVsMtVsMult", "k* vs m_{T} vs multiplicity; k* (GeV/#it{c}); m_{T} (GeV/#it{c}^{2}); Multiplicity"},
      {kKstarVsMtVsMultVsCent, o2::framework::kTHnSparseF, "hKstarVsMtVsMultVsCent", "k* vs m_{T} vs multiplicity vs centrality; k* (GeV/#it{c}); m_{T} (GeV/#it{c}^{2}); Multiplicity; Centrality (%)"},
      {kKstarVsMtVsPt1VsPt2VsMult, o2::framework::kTHnSparseF, "hKstarVsMtVsPt1VsPt2VsMult", "k* vs m_{T} vs p_{T,1} vs p_{T,2} vs multiplicity; k* (GeV/#it{c}); m_{T} (GeV/#it{c}^{2}); p_{T,1} (GeV/#it{c}); p_{T,2} (GeV/#it{c}); Multiplicity"},
      {kKstarVsMtVsPt1VsPt2VsMultVsCent, o2::framework::kTHnSparseF, "hKstarVsMtVsPt1VsPt2VsMultVsCent", "k* vs m_{T} vs p_{T,1} vs p_{T,2} vs multiplicity vs centrality; k* (GeV/#it{c}); m_{T} (GeV/#it{c}^{2}); p_{T,1} (GeV/#it{c}); p_{T,2} (GeV/#it{c}); Multiplicity; Centrality"},
      // n-D with mass
      {kKstarVsMass1VsMass2, o2::framework::kTHnSparseF, "hKstarVsMass1VsMass2", "k* vs m_{1} vs m_{2}; k* (GeV/#it{c}); m_{1} (GeV/#it{c}^{2}); m_{2} (GeV/#it{c}^{2})"},
      {kKstarVsMass1VsMult, o2::framework::kTHnSparseF, "hKstarVsMass1VsMult", "k* vs m_{1} vs multiplicity; k* (GeV/#it{c}); m_{1} (GeV/#it{c}^{2}); Multiplicity"},
      {kKstarVsMass2VsMult, o2::framework::kTHnSparseF, "hKstarVsMass2VsMult", "k* vs m_{2} vs multiplicity; k* (GeV/#it{c}); m_{2} (GeV/#it{c}^{2}); Multiplicity"},
      {kKstarVsMass1VsMass2VsMult, o2::framework::kTHnSparseF, "hKstarVsMass1VsMass2VsMult", "k* vs m_{1} vs m_{2} vs multiplicity; k* (GeV/#it{c}); m_{1} (GeV/#it{c}^{2}); m_{2} (GeV/#it{c}^{2}); Multiplicity"},
      {kTrueKstarVsKstar, o2::framework::kTH2F, "hTrueKstarVsKstar", "k*_{True} vs k*; k*_{True} (GeV/#it{c});  k* (GeV/#it{c})"},
    }};

#define PAIR_HIST_ANALYSIS_MAP(conf)                                                                                   \
  {kKstar, {conf.kstar}},                                                                                              \
    {kKt, {conf.kt}},                                                                                                  \
    {kMt, {conf.mt}},                                                                                                  \
    {kPt1VsPt2, {conf.pt1, conf.pt2}},                                                                                 \
    {kPt1VsKstar, {conf.pt1, conf.kstar}},                                                                             \
    {kPt2VsKstar, {conf.pt2, conf.kstar}},                                                                             \
    {kPt1VsKt, {conf.pt1, conf.kt}},                                                                                   \
    {kPt2VsKt, {conf.pt2, conf.kt}},                                                                                   \
    {kPt1VsMt, {conf.pt1, conf.mt}},                                                                                   \
    {kPt2VsMt, {conf.pt2, conf.mt}},                                                                                   \
    {kKstarVsKt, {conf.kstar, conf.kt}},                                                                               \
    {kKstarVsMt, {conf.kstar, conf.mt}},                                                                               \
    {kKstarVsMult, {conf.kstar, conf.multiplicity}},                                                                   \
    {kKstarVsCent, {conf.kstar, conf.centrality}},                                                                     \
    {kKstarVsMass1, {conf.kstar, conf.mass1}},                                                                         \
    {kKstarVsMass2, {conf.kstar, conf.mass2}},                                                                         \
    {kMass1VsMass2, {conf.mass1, conf.mass2}},                                                                         \
    {kKstarVsMtVsMult, {conf.kstar, conf.mt, conf.multiplicity}},                                                      \
    {kKstarVsMtVsMultVsCent, {conf.kstar, conf.mt, conf.multiplicity, conf.centrality}},                               \
    {kKstarVsMtVsPt1VsPt2VsMult, {conf.kstar, conf.mt, conf.pt1, conf.pt2, conf.multiplicity}},                        \
    {kKstarVsMtVsPt1VsPt2VsMultVsCent, {conf.kstar, conf.mt, conf.pt1, conf.pt2, conf.multiplicity, conf.centrality}}, \
    {kKstarVsMass1VsMass2, {conf.kstar, conf.mass1, conf.mass2}},                                                      \
    {kKstarVsMass1VsMult, {conf.kstar, conf.mass1, conf.multiplicity}},                                                \
    {kKstarVsMass2VsMult, {conf.kstar, conf.mass2, conf.multiplicity}},                                                \
    {kKstarVsMass1VsMass2VsMult, {conf.kstar, conf.mass1, conf.mass2, conf.multiplicity}},

#define PAIR_HIST_MC_MAP(conf) \
  {kTrueKstarVsKstar, {conf.kstar, conf.kstar}},

template <typename T>
auto makePairHistSpecMap(const T& confPairBinning)
{
  return std::map<PairHist, std::vector<framework::AxisSpec>>{
    PAIR_HIST_ANALYSIS_MAP(confPairBinning)};
};

template <typename T>
auto makePairHistMcSpecMap(const T& confPairBinning)
{
  return std::map<PairHist, std::vector<framework::AxisSpec>>{
    PAIR_HIST_ANALYSIS_MAP(confPairBinning)
      PAIR_HIST_MC_MAP(confPairBinning)};
};

#undef PAIR_HIST_ANALYSIS_MAP
#undef PAIR_HIST_MC_MAP

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
constexpr std::string_view McDir = "MC/";

template <const char* prefix,
          modes::Particle particleType1,
          modes::Particle particleType2>
class PairHistManager
{
 public:
  PairHistManager() = default;
  ~PairHistManager() = default;

  template <modes::Mode mode, typename T1, typename T2>
  void init(o2::framework::HistogramRegistry* registry,
            std::map<PairHist, std::vector<o2::framework::AxisSpec>> const& Specs,
            T1 const& ConfPairBinning,
            T2 const& ConfPairCuts)
  {
    mHistogramRegistry = registry;

    // flags for histograms
    mPlot1d = ConfPairBinning.plot1D.value;
    mPlot2d = ConfPairBinning.plot2D.value;
    mPlotKstarVsMtVsMult = ConfPairBinning.plotKstarVsMtVsMult.value;
    mPlotKstarVsMtVsMultVsCent = ConfPairBinning.plotKstarVsMtVsMultVsCent.value;
    mPlotKstarVsMtVsPt1VsP2VsMult = ConfPairBinning.plotKstarVsMtVsPt1VsPt2VsMult.value;
    mPlotKstarVsMtVsPt1VsP2VsMultVsCent = ConfPairBinning.plotKstarVsMtVsPt1VsPt2VsMultVsCent.value;

    mPlotKstarVsMass1VsMass2 = ConfPairBinning.plotKstarVsMass1VsMass2.value;
    mPlotKstarVsMass1VsMult = ConfPairBinning.plotKstarVsMass1VsMult.value;
    mPlotKstarVsMass2VsMult = ConfPairBinning.plotKstarVsMass2VsMult.value;
    mPlotKstarVsMass1VsMass2VsMult = ConfPairBinning.plotKstarVsMass1VsMass2VsMult.value;

    // transverse mass type
    mMtType = static_cast<modes::TransverseMassType>(ConfPairBinning.transverseMassType.value);

    // values for cuts
    mKstarMin = ConfPairCuts.kstarMin.value;
    mKstarMax = ConfPairCuts.kstarMax.value;
    mKtMin = ConfPairCuts.ktMin.value;
    mKtMax = ConfPairCuts.ktMax.value;
    mMtMin = ConfPairCuts.mtMin.value;
    mMtMax = ConfPairCuts.mtMax.value;

    if constexpr (isFlagSet(mode, modes::Mode::kAnalysis)) {
      initAnalysis(Specs);
    }

    // if constexpr (isFlagSet(mode, modes::Mode::kQA)) {
    // std::string qaDir = std::string(prefix) + std::string(QaDir);
    // }
  }

  void setMass(int PdgParticle1, int PdgParticle2)
  {
    mPdgMass1 = o2::analysis::femto::utils::getMass(PdgParticle1);
    mPdgMass2 = o2::analysis::femto::utils::getMass(PdgParticle2);
    mAverageMass = (mPdgMass1 + mPdgMass2) / 2.f;
    mReducedMass = 2.f * (mPdgMass1 * mPdgMass2) / (mPdgMass1 + mPdgMass2);
  }
  void setCharge(int chargeAbsParticle1, int chargeAbsParticle2)
  {
    // the pt stored is actually as pt/z for tracks, so in case of particles with z > 1, we have to rescale the pt (this is so far only for He3 the case)
    mAbsCharge1 = std::abs(chargeAbsParticle1);
    mAbsCharge2 = std::abs(chargeAbsParticle2);
  }

  template <typename T1, typename T2>
  void setPair(const T1& particle1, const T2& particle2)
  {
    // pt in track table is calculated from 1/signedPt from the original track table
    // in case of He with Z=2, we have to rescale the pt with the absolute charge
    mParticle1 = ROOT::Math::PtEtaPhiMVector{mAbsCharge1 * particle1.pt(), particle1.eta(), particle1.phi(), mPdgMass1};
    mParticle2 = ROOT::Math::PtEtaPhiMVector{mAbsCharge2 * particle2.pt(), particle2.eta(), particle2.phi(), mPdgMass2};
    auto partSum = mParticle1 + mParticle2;

    // set kT
    mKt = partSum.Pt() / 2.f;

    // set mT
    computeMt(partSum);

    // Boost particle to the pair rest frame (Prf) and calculate k* (would be equivalent using particle 2)
    // make a copy of particle 1
    auto particle1Prf = ROOT::Math::PtEtaPhiMVector(mParticle1);
    // get lorentz boost into pair rest frame
    ROOT::Math::Boost boostPrf(partSum.BoostToCM());
    // boost particle 1 into pair rest frame and calculate its momentum, which has the same value as k*
    mKstar = boostPrf(particle1Prf).P();

    // if one of the particles has a mass getter, we cache the value for the filling later
    if constexpr (modes::hasMass(particleType1)) {
      mMass1 = particle1.mass();
    }
    if constexpr (modes::hasMass(particleType2)) {
      mMass2 = particle2.mass();
    }
  }

  template <typename T1, typename T2, typename T3>
  void setPair(const T1& particle1, const T2& particle2, const T3& col)
  {
    mMult = col.mult();
    mCent = col.cent();
    setPair(particle1, particle2);
  }

  template <typename T1, typename T2, typename T3, typename T4>
  void setPair(const T1& particle1, const T2& particle2, const T3& col1, const T4& col2)
  {
    mMult = 0.5f * (col1.mult() + col2.mult()); // if mixing with multiplicity, should be in the same mixing bin
    mCent = 0.5f * (col1.cent() + col2.cent()); // if mixing with centrality, should be in the same mixing bin
    setPair(particle1, particle2);
  }

  bool checkPairCuts() const
  {
    return (!(mKstarMin > 0.f) || mKstar > mKstarMin) &&
           (!(mKstarMax > 0.f) || mKstar < mKstarMax) &&
           (!(mKtMin > 0.f) || mKt > mKtMin) &&
           (!(mKtMax > 0.f) || mKt < mKtMax) &&
           (!(mMtMin > 0.f) || mMt > mMtMin) &&
           (!(mMtMax > 0.f) || mMt < mMtMax);
  }

  template <modes::Mode mode>
  void fill()
  {
    if constexpr (isFlagSet(mode, modes::Mode::kAnalysis)) {
      fillAnalysis();
    }
  }

  float getKstar() const { return mKstar; }

 private:
  void initAnalysis(std::map<PairHist, std::vector<o2::framework::AxisSpec>> const& Specs)
  {
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

      // special care for mass plots since not all particles have "mass"
      if constexpr (modes::hasMass(particleType1)) {
        mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsMass1, HistTable), getHistDesc(kKstarVsMass1, HistTable), getHistType(kKstarVsMass1, HistTable), {Specs.at(kKstarVsMass1)});
      }
      if constexpr (modes::hasMass(particleType2)) {
        mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsMass2, HistTable), getHistDesc(kKstarVsMass2, HistTable), getHistType(kKstarVsMass2, HistTable), {Specs.at(kKstarVsMass2)});
      }
      if constexpr (modes::hasMass(particleType1) && modes::hasMass(particleType2)) {
        mHistogramRegistry->add(analysisDir + getHistNameV2(kMass1VsMass2, HistTable), getHistDesc(kMass1VsMass2, HistTable), getHistType(kMass1VsMass2, HistTable), {Specs.at(kMass1VsMass2)});
      }
    }

    if (mPlotKstarVsMtVsMult) {
      mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsMtVsMult, HistTable), getHistDesc(kKstarVsMtVsMult, HistTable), getHistType(kKstarVsMtVsMult, HistTable), {Specs.at(kKstarVsMtVsMult)});
    }
    if (mPlotKstarVsMtVsMultVsCent) {
      mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsMtVsMultVsCent, HistTable), getHistDesc(kKstarVsMtVsMultVsCent, HistTable), getHistType(kKstarVsMtVsMultVsCent, HistTable), {Specs.at(kKstarVsMtVsMultVsCent)});
    }
    if (mPlotKstarVsMtVsPt1VsP2VsMult) {
      mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsMtVsPt1VsPt2VsMult, HistTable), getHistDesc(kKstarVsMtVsPt1VsPt2VsMult, HistTable), getHistType(kKstarVsMtVsPt1VsPt2VsMult, HistTable), {Specs.at(kKstarVsMtVsPt1VsPt2VsMult)});
    }
    if (mPlotKstarVsMtVsPt1VsP2VsMultVsCent) {
      mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsMtVsPt1VsPt2VsMultVsCent, HistTable), getHistDesc(kKstarVsMtVsPt1VsPt2VsMultVsCent, HistTable), getHistType(kKstarVsMtVsPt1VsPt2VsMultVsCent, HistTable), {Specs.at(kKstarVsMtVsPt1VsPt2VsMultVsCent)});
    }

    // again special care for particles with "mass"
    if constexpr (modes::hasMass(particleType1)) {
      if (mPlotKstarVsMass1VsMult) {
        mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsMass1VsMult, HistTable), getHistDesc(kKstarVsMass1VsMult, HistTable), getHistType(kKstarVsMass1VsMult, HistTable), {Specs.at(kKstarVsMass1VsMult)});
      }
    }
    if constexpr (modes::hasMass(particleType2)) {
      if (mPlotKstarVsMass2VsMult) {
        mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsMass2VsMult, HistTable), getHistDesc(kKstarVsMass2VsMult, HistTable), getHistType(kKstarVsMass2VsMult, HistTable), {Specs.at(kKstarVsMass2VsMult)});
      }
    }
    if constexpr (modes::hasMass(particleType1) && modes::hasMass(particleType2)) {
      if (mPlotKstarVsMass1VsMass2) {
        mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsMass1VsMass2, HistTable), getHistDesc(kKstarVsMass1VsMass2, HistTable), getHistType(kKstarVsMass1VsMass2, HistTable), {Specs.at(kKstarVsMass1VsMass2)});
      }
      if (mPlotKstarVsMass1VsMass2VsMult) {
        mHistogramRegistry->add(analysisDir + getHistNameV2(kKstarVsMass1VsMass2VsMult, HistTable), getHistDesc(kKstarVsMass1VsMass2VsMult, HistTable), getHistType(kKstarVsMass1VsMass2VsMult, HistTable), {Specs.at(kKstarVsMass1VsMass2VsMult)});
      }
    }
  }

  void fillAnalysis()
  {
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

      // // special care for mass plots since not all particles have "mass"
      if constexpr (modes::hasMass(particleType1)) {
        mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsMass1, HistTable)), mKstar, mMass1);
      }
      if constexpr (modes::hasMass(particleType2)) {
        mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsMass2, HistTable)), mKstar, mMass2);
      }
      if constexpr (modes::hasMass(particleType1) && modes::hasMass(particleType2)) {
        mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kMass1VsMass2, HistTable)), mMass1, mMass2);
      }
    }

    // n-D histograms are only filled if enabled
    if (mPlotKstarVsMtVsMult) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsMtVsMult, HistTable)), mKstar, mMt, mMult);
    }
    if (mPlotKstarVsMtVsMultVsCent) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsMtVsMultVsCent, HistTable)), mKstar, mMt, mMult, mCent);
    }
    if (mPlotKstarVsMtVsPt1VsP2VsMult) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsMtVsPt1VsPt2VsMult, HistTable)), mKstar, mMt, mParticle1.Pt(), mParticle2.pt(), mMult);
    }
    if (mPlotKstarVsMtVsPt1VsP2VsMultVsCent) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsMtVsPt1VsPt2VsMultVsCent, HistTable)), mKstar, mMt, mParticle1.Pt(), mParticle2.pt(), mMult, mCent);
    }

    // again special care for particles with "mass"
    if constexpr (modes::hasMass(particleType1)) {
      if (mPlotKstarVsMass1VsMult) {
        mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsMass1VsMult, HistTable)), mKstar, mMass1, mMult);
      }
    }
    if constexpr (modes::hasMass(particleType2)) {
      if (mPlotKstarVsMass2VsMult) {
        mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsMass2VsMult, HistTable)), mKstar, mMass2, mMult);
      }
    }
    if constexpr (modes::hasMass(particleType1) && modes::hasMass(particleType2)) {
      if (mPlotKstarVsMass1VsMass2) {
        mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsMass1VsMass2, HistTable)), mKstar, mMass1, mMass2);
      }
      if (mPlotKstarVsMass1VsMass2VsMult) {
        mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kKstarVsMass1VsMass2VsMult, HistTable)), mKstar, mMass1, mMass2, mMult);
      }
    }
  }

  void computeMt(ROOT::Math::PtEtaPhiMVector const& PairMomentum)
  {
    switch (mMtType) {
      case modes::TransverseMassType::kAveragePdgMass:
        mMt = std::hypot(mKt, mAverageMass);
        break;
      case modes::TransverseMassType::kReducedPdgMass:
        mMt = std::hypot(mKt, mReducedMass);
        break;
      case modes::TransverseMassType::kMt4Vector:
        mMt = PairMomentum.Mt() / 2.f;
        break;
      default:
        mMt = std::hypot(mKt, mAverageMass);
    }
  }

  o2::framework::HistogramRegistry* mHistogramRegistry = nullptr;
  float mPdgMass1 = 0.f;
  float mPdgMass2 = 0.f;

  modes::TransverseMassType mMtType = modes::TransverseMassType::kAveragePdgMass;
  float mAverageMass = 0.f;
  float mReducedMass = 0.f;

  int mAbsCharge1 = 1;
  int mAbsCharge2 = 1;
  ROOT::Math::PtEtaPhiMVector mParticle1{};
  ROOT::Math::PtEtaPhiMVector mParticle2{};
  float mMass1 = 0.f;
  float mMass2 = 0.f;
  float mKstar = 0.f;
  float mKt = 0.f;
  float mMt = 0.f;
  float mMult = 0.f;
  float mCent = 0.f;

  // cuts
  float mKstarMin = -1.f;
  float mKstarMax = -1.f;
  float mKtMin = -1.f;
  float mKtMax = -1.f;
  float mMtMin = -1.f;
  float mMtMax = -1.f;

  // flags
  bool mPlot1d = true;
  bool mPlot2d = true;
  bool mPlotKstarVsMtVsMult = false;
  bool mPlotKstarVsMtVsMultVsCent = false;
  bool mPlotKstarVsMtVsPt1VsP2VsMult = false;
  bool mPlotKstarVsMtVsPt1VsP2VsMultVsCent = false;

  bool mPlotKstarVsMass1VsMass2 = false;
  bool mPlotKstarVsMass1VsMult = false;
  bool mPlotKstarVsMass2VsMult = false;
  bool mPlotKstarVsMass1VsMass2VsMult = false;
};

}; // namespace pairhistmanager
}; // namespace o2::analysis::femto
#endif // PWGCF_FEMTO_CORE_PAIRHISTMANAGER_H_
