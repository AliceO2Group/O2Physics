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

/// \file tripletHistManager.h
/// \brief histogram manager for triplet tasks
/// \author anton.riedel@tum.de, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTO_CORE_TRIPLETHISTMANAGER_H_
#define PWGCF_FEMTO_CORE_TRIPLETHISTMANAGER_H_

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
namespace triplethistmanager
{
// enum for pair histograms
enum TripletHist {
  // standard 1D
  kQ3,
  kMt, // averate mt of all pairs in the triplet
  // standard 2D
  kPt1VsQ3,
  kPt2VsQ3,
  kPt3VsQ3,
  kQ3VsMt,
  kQ3VsMult,
  kQ3VsCent,
  // 2D with mass
  // kQ3VsMass1,
  // kQ3VsMass2,
  // kQ3VsMass3,
  // higher dimensions
  kPt1VsPt2VsPt3,
  kQ3VsPt1VsPt2VsPt3,
  kQ3VsMtVsMult,
  kQ3VsMtVsMultVsCent,
  kQ3VsMtVsPt1VsPt2VsPt3VsMult,
  kQ3VsMtVsPt1VsPt2VsPt3VsMultVsCent,
  // higher dimensions with mass
  // mc
  kTrueQ3VsQ3,
  kTrueMtVsMt,
  kTrueMultVsMult,
  kTrueCentVsCent,

  kTripletHistogramLast
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
  o2::framework::Configurable<int> seed{"seed", -1, "Seed to randomize particle 1/2/3 (if they are identical). Set to negative value to deactivate. Set to 0 to generate unique seed in time."};
  o2::framework::Configurable<bool> particle123AreSameSpecies{"particle123AreSameSpecies", false, "Particle 1,2 and 3 are of the same species"};
  o2::framework::Configurable<bool> particle12AreSameSpecies{"particle12AreSameSpecies", false, "Particle 1 and 2 are of the same species"};
};

struct ConfTripletBinning : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("TripletBinning");
  o2::framework::Configurable<bool> plot1D{"plot1D", true, "Enable 1D histograms"};
  o2::framework::Configurable<bool> plot2D{"plot2D", true, "Enable 2D histograms"};
  o2::framework::Configurable<bool> plotPt1VsPt2VsPt3{"plotPt1VsPt2VsPt3", false, "Enable 3D histogram (Pt1 vs Pt2 vs Pt3)"};
  o2::framework::Configurable<bool> plotQ3VsPt1VsPt2VsPt3{"plotQ3VsPt1VsPt2VsPt3", false, "Enable 4D histogram (Q3 vs Pt1 vs Pt2 vs Pt3)"};
  o2::framework::Configurable<bool> plotQ3VsMtVsMult{"plotQ3VsMtVsMult", false, "Enable 3D histogram (Q3 vs Mt vs Mult)"};
  o2::framework::Configurable<bool> plotQ3VsMtVsMultVsCent{"plotQ3VsMtVsMultVsCent", false, "Enable 4D histogram (Q3 vs Mt vs Mult Vs Cent)"};
  o2::framework::Configurable<bool> plotQ3VsMtVsPt1VsPt2VsPt3VsMult{"plotQ3VsMtVsPt1VsPt2VsPt3VsMult", false, "Enable 6D histogram (Q3 vs Mt Vs Pt1 vs Pt2 vs Pt3 vs Mult)"};
  o2::framework::Configurable<bool> plotQ3VsMtVsPt1VsPt2VsPt3VsMultVsCent{"plotQ3VsMtVsPt1VsPt2VsPt3VsMultVsCent", false, "Enable 3D histogram (Q3 vs Mt Vs Pt1 vs Pt2 vs Pt3 vs Mult vs Cent)"};
  o2::framework::ConfigurableAxis q3{"q3", {{600, 0, 6}}, "q3"};
  o2::framework::ConfigurableAxis mt{"mt", {{500, 0.8, 5.8}}, "mt"};
  o2::framework::ConfigurableAxis multiplicity{"multiplicity", {{50, 0, 200}}, "multiplicity"};
  o2::framework::ConfigurableAxis centrality{"centrality", {{10, 0, 100}}, "centrality (mult. percentile)"};
  o2::framework::ConfigurableAxis pt1{"pt1", {{100, 0, 6}}, "Pt binning for particle 1"};
  o2::framework::ConfigurableAxis pt2{"pt2", {{100, 0, 6}}, "Pt binning for particle 2"};
  o2::framework::ConfigurableAxis pt3{"pt3", {{100, 0, 6}}, "Pt binning for particle 3"};
  o2::framework::ConfigurableAxis mass1{"mass1", {{100, 0, 2}}, "Mass binning for particle 1 (if particle has mass getter)"};
  o2::framework::ConfigurableAxis mass2{"mass2", {{100, 0, 2}}, "Mass binning for particle 2 (if particle has mass getter)"};
  o2::framework::ConfigurableAxis mass3{"mass3", {{100, 0, 2}}, "Mass binning for particle 3 (if particle has mass getter)"};
  o2::framework::Configurable<int> transverseMassType{"transverseMassType", static_cast<int>(modes::TransverseMassType::kAveragePdgMass), "Type of transverse mass (0-> Average Pdg Mass, 1-> Reduced Pdg Mass, 2-> Mt from combined 4 vector)"};
};

struct ConfTripletCuts : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("TripletCuts");
  o2::framework::Configurable<float> q3Max{"q3Max", -1, "Maximal kstar (set to -1 to deactivate)"};
  o2::framework::Configurable<float> q3Min{"q3Min", -1, "Minimal kstar (set to -1 to deactivate)"};
  o2::framework::Configurable<float> mtMax{"mtMax", -1, "Maximal mt (set to -1 to deactivate)"};
  o2::framework::Configurable<float> mtMin{"mtMin", -1, "Minimal mt (set to -1 to deactivate)"};
  o2::framework::Configurable<bool> mixOnlyCommonAncestor{"mixOnlyCommonAncestor", false, "Require pair to have common anchestor (in the same event)"};
  o2::framework::Configurable<bool> mixOnlyNonCommonAncestor{"mixOnlyNonCommonAncestor", false, "Require pair to have non-common anchestor (in the same event)"};
};

// the enum gives the correct index in the array
constexpr std::array<histmanager::HistInfo<TripletHist>, kTripletHistogramLast>
  HistTable = {
    {
      // 1D
      {kQ3, o2::framework::kTH1F, "hQ3", "Q_{3}; Q_{3} (GeV/#it{c}); Entries"},
      {kMt, o2::framework::kTH1F, "hMt", "transverse mass; m_{T} (GeV/#it{c}^{2}); Entries"},
      // 2D
      {kPt1VsQ3, o2::framework::kTH2F, "hPt1VsQ3", " p_{T,1} vs Q_{3}; p_{T,1} (GeV/#it{c}); Q_{3} (GeV/#it{c})"},
      {kPt2VsQ3, o2::framework::kTH2F, "hPt2VsQ3", "p_{T,2} vs Q_{3}; p_{T,2} (GeV/#it{c}); Q_{3} (GeV/#it{c})"},
      {kPt3VsQ3, o2::framework::kTH2F, "hPt3VsQ3", "p_{T,3} vs Q_{3}; p_{T,2} (GeV/#it{c}); Q_{3} (GeV/#it{c})"},
      {kQ3VsMt, o2::framework::kTH2F, "hQ3VsMt", "Q_{3} vs m_{T}; Q_{3} (GeV/#it{c}); m_{T} (GeV/#it{c}^{2})"},
      {kQ3VsMult, o2::framework::kTH2F, "hQ3VsMult", "Q_{3} vs Multiplicity; Q_{3} (GeV/#it{c}); Multiplicity"},
      {kQ3VsCent, o2::framework::kTH2F, "hQ3VsCent", "Q_{3} vs Centrality; Q_{3} (GeV/#it{c}); Centrality"},
      // n-D
      {kPt1VsPt2VsPt3, o2::framework::kTHnSparseF, "hPt1VsPt2VsPt3", "k* vs m_{T} vs multiplicity; k* (GeV/#it{c}); m_{T} (GeV/#it{c}^{2}); Multiplicity"},
      {kQ3VsPt1VsPt2VsPt3, o2::framework::kTHnSparseF, "hQ3VsPt1VsPt2VsPt3", "Q_{3} vs p_{T,1} vs p_{T,2} vs p_{T,3}; Q_{3} (GeV/#it{c}); p_{T,1} (GeV/#it{c}) ; p_{T,2} (GeV/#it{c}); p_{T,3} (GeV/#it{c})"},
      {kQ3VsMtVsMult, o2::framework::kTHnSparseF, "hQ3VsMtVsMult", "Q_{3} vs m_{T} vs multiplicity; Q_{3} (GeV/#it{c}); m_{T} (GeV/#it{c}^{2}); Multiplicity"},
      {kQ3VsMtVsMultVsCent, o2::framework::kTHnSparseF, "hQ3VsMtVsMultVsCent", "Q_{3} vs m_{T} vs multiplicity vs centrality; Q_{3} (GeV/#it{c}); m_{T} (GeV/#it{c}^{2}); Multiplicity; Centrality"},
      // n-D with mass
      {kQ3VsMtVsPt1VsPt2VsPt3VsMult, o2::framework::kTHnSparseF, "hQ3VsMtVsPt1VsPt2VsPt3VsMult", "Q_{3} vs m_{T} vs p_{T,1} vs p_{T,2} vs p_{T,3} vs multiplicity; Q_{3} (GeV/#it{c}); m_{T} (GeV/#it{c}^{2}); p_{T,1} (GeV/#it{c}) ; p_{T,2} (GeV/#it{c}); p_{T,3} (GeV/#it{c}) ; Multiplicity"},
      {kQ3VsMtVsPt1VsPt2VsPt3VsMultVsCent, o2::framework::kTHnSparseF, "hQ3VsMtVsPt1VsPt2VsPt3VsMultVsCent", "Q_{3} vs m_{T} vs p_{T,1} vs p_{T,2} vs p_{T,3} vs multiplicity vs centrality; Q_{3} (GeV/#it{c}); m_{T} (GeV/#it{c}^{2}); p_{T,1} (GeV/#it{c}) ; p_{T,2} (GeV/#it{c}); p_{T,3} (GeV/#it{c}) vs multiplicity vs centrality"},
      {kTrueQ3VsQ3, o2::framework::kTH2F, "hTrueQ3VsQ3", "Q_{3,True} vs Q_{3}; Q_{3,true} (GeV/#it{c});  Q_{3} (GeV/#it{c})"},
      {kTrueMtVsMt, o2::framework::kTH2F, "hTrueMtVsMt", "m_{T,True} vs m_{T}; m_{T,True} (GeV/#it{c}^{2}); m_{T} (GeV/#it{c}^{2})"},
      {kTrueMultVsMult, o2::framework::kTH2F, "hTrueMultVsMult", "Multiplicity_{True} vs Multiplicity; Multiplicity_{True} ;  Multiplicity"},
      {kTrueCentVsCent, o2::framework::kTH2F, "hTrueCentVsCent", "Centrality_{True} vs Centrality; Centrality_{True} (%); Centrality (%)"},
    }};

#define TRIPLET_HIST_ANALYSIS_MAP(conf)                                                                  \
  {kQ3, {conf.q3}},                                                                                      \
    {kMt, {conf.mt}},                                                                                    \
    {kPt1VsQ3, {conf.pt1, conf.q3}},                                                                     \
    {kPt2VsQ3, {conf.pt2, conf.q3}},                                                                     \
    {kPt3VsQ3, {conf.pt3, conf.q3}},                                                                     \
    {kQ3VsMt, {conf.q3, conf.mt}},                                                                       \
    {kQ3VsMult, {conf.q3, conf.multiplicity}},                                                           \
    {kQ3VsCent, {conf.q3, conf.centrality}},                                                             \
    {kPt1VsPt2VsPt3, {conf.pt1, conf.pt2, conf.pt3}},                                                    \
    {kQ3VsPt1VsPt2VsPt3, {conf.q3, conf.pt1, conf.pt2, conf.pt3}},                                       \
    {kQ3VsMtVsMult, {conf.q3, conf.mt, conf.multiplicity}},                                              \
    {kQ3VsMtVsMultVsCent, {conf.q3, conf.mt, conf.multiplicity, conf.centrality}},                       \
    {kQ3VsMtVsPt1VsPt2VsPt3VsMult, {conf.q3, conf.mt, conf.pt1, conf.pt2, conf.pt3, conf.multiplicity}}, \
    {kQ3VsMtVsPt1VsPt2VsPt3VsMultVsCent, {conf.q3, conf.mt, conf.pt1, conf.pt2, conf.pt3, conf.multiplicity, conf.centrality}},

#define TRIPLET_HIST_MC_MAP(conf)                              \
  {kTrueQ3VsQ3, {conf.q3, conf.q3}},                           \
    {kTrueMtVsMt, {conf.mt, conf.mt}},                         \
    {kTrueMultVsMult, {conf.multiplicity, conf.multiplicity}}, \
    {kTrueCentVsCent, {conf.centrality, conf.centrality}},

template <typename T>
auto makeTripletHistSpecMap(const T& confPairBinning)
{
  return std::map<TripletHist, std::vector<framework::AxisSpec>>{
    TRIPLET_HIST_ANALYSIS_MAP(confPairBinning)};
};

template <typename T>
auto makeTripletMcHistSpecMap(const T& confPairBinning)
{
  return std::map<TripletHist, std::vector<framework::AxisSpec>>{
    TRIPLET_HIST_ANALYSIS_MAP(confPairBinning)
      TRIPLET_HIST_MC_MAP(confPairBinning)};
};

#undef TRIPLET_HIST_ANALYSIS_MAP
#undef TRIPLET_HIST_MC_MAP

constexpr char PrefixTrackTrackTrackSe[] = "TrackTrackTrack/SE/";
constexpr char PrefixTrackTrackTrackMe[] = "TrackTrackTrack/ME/";

constexpr std::string_view AnalysisDir = "Analysis/";
constexpr std::string_view QaDir = "QA/";
constexpr std::string_view McDir = "MC/";

template <const char* prefix,
          modes::Particle particleType1,
          modes::Particle particleType2,
          modes::Particle particleType3>
class TripletHistManager
{
 public:
  TripletHistManager() = default;
  ~TripletHistManager() = default;

  template <modes::Mode mode, typename T1, typename T2>
  void init(o2::framework::HistogramRegistry* registry,
            std::map<TripletHist, std::vector<o2::framework::AxisSpec>> const& Specs,
            T1 const& ConfTripletBinning,
            T2 const& ConfTripletCuts)
  {
    mHistogramRegistry = registry;

    // flags for histograms
    mPlot1d = ConfTripletBinning.plot1D.value;
    mPlot2d = ConfTripletBinning.plot2D.value;

    // transverse mass type
    mMtType = static_cast<modes::TransverseMassType>(ConfTripletBinning.transverseMassType.value);

    // values for cuts
    mQ3Min = ConfTripletCuts.q3Min.value;
    mQ3Max = ConfTripletCuts.q3Max.value;
    mMtMin = ConfTripletCuts.mtMin.value;
    mMtMax = ConfTripletCuts.mtMax.value;

    mPlotPt1VsPt2VsPt3 = ConfTripletBinning.plotPt1VsPt2VsPt3.value;
    mPlotQ3VsPt1VsPt2VsPt3 = ConfTripletBinning.plotQ3VsPt1VsPt2VsPt3.value;
    mPlotQ3VsMtVsMult = ConfTripletBinning.plotQ3VsMtVsMult.value;
    mPlotQ3VsMtVsMultVsCent = ConfTripletBinning.plotQ3VsMtVsMultVsCent.value;
    mPlotQ3VsMtVsPt1VsPt2VsPt3VsMult = ConfTripletBinning.plotQ3VsMtVsPt1VsPt2VsPt3VsMult.value;
    mPlotQ3VsMtVsPt1VsPt2VsPt3VsMultVsCent = ConfTripletBinning.plotQ3VsMtVsPt1VsPt2VsPt3VsMultVsCent.value;

    if constexpr (isFlagSet(mode, modes::Mode::kAnalysis)) {
      initAnalysis(Specs);
    }

    if constexpr (isFlagSet(mode, modes::Mode::kMc)) {
      initMc(Specs);
    }
  }

  void setMass(int PdgParticle1, int PdgParticle2, int PdgParticle3)
  {
    mPdgMass1 = o2::analysis::femto::utils::getMass(PdgParticle1);
    mPdgMass2 = o2::analysis::femto::utils::getMass(PdgParticle2);
    mPdgMass3 = o2::analysis::femto::utils::getMass(PdgParticle3);
  }
  void setCharge(int chargeAbsParticle1, int chargeAbsParticle2, int chargeAbsParticle3)
  {
    // the pt stored is actually as pt/z for tracks, so in case of particles with z > 1, we have to rescale the pt (this is so far only for He3 the case)
    mAbsCharge1 = std::abs(chargeAbsParticle1);
    mAbsCharge2 = std::abs(chargeAbsParticle2);
    mAbsCharge3 = std::abs(chargeAbsParticle3);
  }

  template <typename T1, typename T2, typename T3>
  void setTriplet(const T1& particle1, const T2& particle2, const T3& particle3)
  {
    // pt in track table is calculated from 1/signedPt from the original track table
    // in case of He with Z=2, we have to rescale the pt with the absolute charge
    mParticle1 = ROOT::Math::PtEtaPhiMVector(mAbsCharge1 * particle1.pt(), particle1.eta(), particle1.phi(), mPdgMass1);
    mParticle2 = ROOT::Math::PtEtaPhiMVector(mAbsCharge2 * particle2.pt(), particle2.eta(), particle2.phi(), mPdgMass2);
    mParticle3 = ROOT::Math::PtEtaPhiMVector(mAbsCharge3 * particle3.pt(), particle3.eta(), particle3.phi(), mPdgMass3);

    // set mT
    mMt = getMt(mParticle1, mPdgMass1, mParticle2, mPdgMass2, mParticle3, mPdgMass3);
    // set Q3
    mQ3 = getQ3(mParticle1, mPdgMass1, mParticle2, mPdgMass2, mParticle3, mPdgMass3);

    // if one of the particles has a mass getter, we cache the value for the filling later
    if constexpr (modes::hasMass(particleType1)) {
      mMass1 = particle1.mass();
    }
    if constexpr (modes::hasMass(particleType2)) {
      mMass2 = particle2.mass();
    }
    if constexpr (modes::hasMass(particleType3)) {
      mMass3 = particle2.mass();
    }
  }

  template <typename T1, typename T2, typename T3, typename T4>
  void setTriplet(T1 const& particle1, T2 const& particle2, T3 const& particle3, T4 const& col)
  {
    setTriplet(particle1, particle2, particle3);
    mMult = col.mult();
    mCent = col.cent();
  }

  template <typename T1, typename T2, typename T3, typename T4>
  void setTriplet(T1 const& particle1, T2 const& particle2, T3 const& particle3, T4 const& col1, T4 const& col2, T4 const& col3)
  {
    setTriplet(particle1, particle2, particle3);
    mMult = (col1.mult() + col2.mult() + col3.mult()) / 3.f; // if mixing with multiplicity, should be in the same mixing bin
    mCent = (col1.cent() + col2.cent() + col3.cent()) / 3.f; // if mixing with centrality, should be in the same mixing bin
  }

  template <typename T1, typename T2, typename T3, typename T4>
  void setTripletMc(T1 const& particle1, T2 const& particle2, T3 const& particle3, const T4& /*mcParticles*/)
  {
    if (!particle1.has_fMcParticle() || !particle2.has_fMcParticle() || !particle3.has_fMcParticle()) {
      mHasMcTriplet = false;
      return;
    }
    mHasMcTriplet = true;

    auto mcParticle1 = particle1.template fMcParticle_as<T4>();
    auto mcParticle2 = particle2.template fMcParticle_as<T4>();
    auto mcParticle3 = particle3.template fMcParticle_as<T4>();

    mTrueParticle1 = ROOT::Math::PtEtaPhiMVector(mAbsCharge1 * mcParticle1.pt(), mcParticle1.eta(), mcParticle1.phi(), mPdgMass1);
    mTrueParticle2 = ROOT::Math::PtEtaPhiMVector(mAbsCharge2 * mcParticle2.pt(), mcParticle2.eta(), mcParticle2.phi(), mPdgMass2);
    mTrueParticle3 = ROOT::Math::PtEtaPhiMVector(mAbsCharge3 * mcParticle3.pt(), mcParticle3.eta(), mcParticle3.phi(), mPdgMass3);

    // set true mT
    mTrueMt = getMt(mTrueParticle1, mPdgMass1, mTrueParticle2, mPdgMass2, mTrueParticle3, mPdgMass3);
    // set true Q3
    mTrueQ3 = getQ3(mTrueParticle1, mPdgMass1, mTrueParticle2, mPdgMass2, mTrueParticle3, mPdgMass3);
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  void setTripletMc(T1 const& particle1, T2 const& particle2, T3 const& particle3, T4 const& mcParticles, T5 const& col, const T6& /*mcCols*/)
  {
    setTriplet(particle1, particle2, particle3, col);
    setTripletMc(particle1, particle2, particle3, mcParticles);
    if (!col.has_fMcCol()) {
      mHasMcCol = false;
      return;
    }
    mHasMcCol = true;
    auto mcCol = col.template fMcCol_as<T6>();
    mTrueMult = mcCol.mult();
    mTrueCent = mcCol.cent();
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  void setTripletMc(T1 const& particle1, T2 const& particle2, T3 const& particle3, T4 const& mcParticles, T5 const& col1, T5 const& col2, T5 const& col3, const T6& /*mcCols*/)
  {
    setTriplet(particle1, particle2, particle3, col1, col2, col3);
    setTripletMc(particle1, particle2, particle3, mcParticles);
    if (!col1.has_fMcCol() || !col2.has_fMcCol() || !col2.has_fMcCol()) {
      mHasMcCol = false;
      return;
    }
    mHasMcCol = true;
    auto mcCol1 = col1.template fMcCol_as<T6>();
    auto mcCol2 = col2.template fMcCol_as<T6>();
    auto mcCol3 = col2.template fMcCol_as<T6>();
    mTrueMult = (mcCol1.mult() + mcCol2.mult() + mcCol3.mult()) / 3.f;
    mTrueCent = (mcCol1.cent() + mcCol2.cent() + mcCol3.cent()) / 3.f;
  }

  bool checkTripletCuts() const
  {
    return (!(mQ3Min > 0.f) || mQ3 > mQ3Min) &&
           (!(mQ3Max > 0.f) || mQ3 < mQ3Max) &&
           (!(mMtMin > 0.f) || mMt > mMtMin) &&
           (!(mMtMax > 0.f) || mMt < mMtMax);
  }

  template <modes::Mode mode>
  void fill()
  {
    if constexpr (isFlagSet(mode, modes::Mode::kAnalysis)) {
      fillAnalysis();
    }
    if constexpr (isFlagSet(mode, modes::Mode::kMc)) {
      fillMc();
    }
  }

  float getQ3() const { return mQ3; }

 private:
  template <typename T>
  T getqij(const T& vecparti, const T& vecpartj)
  {
    T trackSum = vecparti + vecpartj;
    T trackDifference = vecparti - vecpartj;
    float scaling = trackDifference.Dot(trackSum) / trackSum.Dot(trackSum);
    return trackDifference - scaling * trackSum;
  }

  template <typename T>
  float getQ3(const T& part1, const float mass1, const T& part2, const float mass2, const T& part3, const float mass3)
  {
    ROOT::Math::PtEtaPhiMVector vecpart1(part1.pt(), part1.eta(), part1.phi(), mass1);
    ROOT::Math::PtEtaPhiMVector vecpart2(part2.pt(), part2.eta(), part2.phi(), mass2);
    ROOT::Math::PtEtaPhiMVector vecpart3(part3.pt(), part3.eta(), part3.phi(), mass3);
    auto q12 = getqij(vecpart1, vecpart2);
    auto q23 = getqij(vecpart2, vecpart3);
    auto q31 = getqij(vecpart3, vecpart1);
    float q = q12.M2() + q23.M2() + q31.M2();
    return std::sqrt(-q);
  }

  void initAnalysis(std::map<TripletHist, std::vector<o2::framework::AxisSpec>> const& Specs)
  {
    std::string analysisDir = std::string(prefix) + std::string(AnalysisDir);
    if (mPlot1d) {
      mHistogramRegistry->add(analysisDir + getHistNameV2(kQ3, HistTable), getHistDesc(kQ3, HistTable), getHistType(kQ3, HistTable), {Specs.at(kQ3)});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kMt, HistTable), getHistDesc(kMt, HistTable), getHistType(kMt, HistTable), {Specs.at(kMt)});
    }
    if (mPlot2d) {
      mHistogramRegistry->add(analysisDir + getHistNameV2(kPt1VsQ3, HistTable), getHistDesc(kPt1VsQ3, HistTable), getHistType(kPt1VsQ3, HistTable), {Specs.at(kPt1VsQ3)});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kPt2VsQ3, HistTable), getHistDesc(kPt2VsQ3, HistTable), getHistType(kPt2VsQ3, HistTable), {Specs.at(kPt2VsQ3)});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kPt3VsQ3, HistTable), getHistDesc(kPt3VsQ3, HistTable), getHistType(kPt3VsQ3, HistTable), {Specs.at(kPt3VsQ3)});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kQ3VsMt, HistTable), getHistDesc(kQ3VsMt, HistTable), getHistType(kQ3VsMt, HistTable), {Specs.at(kQ3VsMt)});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kQ3VsMult, HistTable), getHistDesc(kQ3VsMult, HistTable), getHistType(kQ3VsMult, HistTable), {Specs.at(kQ3VsMult)});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kQ3VsCent, HistTable), getHistDesc(kQ3VsCent, HistTable), getHistType(kQ3VsCent, HistTable), {Specs.at(kQ3VsCent)});
      // TODO: implement histograms differential im mass of the particles
    }

    if (mPlotPt1VsPt2VsPt3) {
      mHistogramRegistry->add(analysisDir + getHistNameV2(kPt1VsPt2VsPt3, HistTable), getHistDesc(kPt1VsPt2VsPt3, HistTable), getHistType(kPt1VsPt2VsPt3, HistTable), {Specs.at(kPt1VsPt2VsPt3)});
    }
    if (mPlotQ3VsPt1VsPt2VsPt3) {
      mHistogramRegistry->add(analysisDir + getHistNameV2(kQ3VsPt1VsPt2VsPt3, HistTable), getHistDesc(kQ3VsPt1VsPt2VsPt3, HistTable), getHistType(kQ3VsPt1VsPt2VsPt3, HistTable), {Specs.at(kQ3VsPt1VsPt2VsPt3)});
    }
    if (mPlotQ3VsMtVsMult) {
      mHistogramRegistry->add(analysisDir + getHistNameV2(kQ3VsMtVsMult, HistTable), getHistDesc(kQ3VsMtVsMult, HistTable), getHistType(kQ3VsMtVsMult, HistTable), {Specs.at(kQ3VsMtVsMult)});
    }
    if (mPlotQ3VsMtVsMultVsCent) {
      mHistogramRegistry->add(analysisDir + getHistNameV2(kQ3VsMtVsMultVsCent, HistTable), getHistDesc(kQ3VsMtVsMultVsCent, HistTable), getHistType(kQ3VsMtVsMultVsCent, HistTable), {Specs.at(kQ3VsMtVsMultVsCent)});
    }
    if (mPlotQ3VsMtVsPt1VsPt2VsPt3VsMult) {
      mHistogramRegistry->add(analysisDir + getHistNameV2(kQ3VsMtVsPt1VsPt2VsPt3VsMult, HistTable), getHistDesc(kQ3VsMtVsPt1VsPt2VsPt3VsMult, HistTable), getHistType(kQ3VsMtVsPt1VsPt2VsPt3VsMult, HistTable), {Specs.at(kQ3VsMtVsPt1VsPt2VsPt3VsMult)});
    }
    if (mPlotQ3VsMtVsPt1VsPt2VsPt3VsMultVsCent) {
      mHistogramRegistry->add(analysisDir + getHistNameV2(kQ3VsMtVsPt1VsPt2VsPt3VsMultVsCent, HistTable), getHistDesc(kQ3VsMtVsPt1VsPt2VsPt3VsMultVsCent, HistTable), getHistType(kQ3VsMtVsPt1VsPt2VsPt3VsMultVsCent, HistTable), {Specs.at(kQ3VsMtVsPt1VsPt2VsPt3VsMultVsCent)});
    }
  }

  void initMc(std::map<TripletHist, std::vector<o2::framework::AxisSpec>> const& Specs)
  {
    std::string mcDir = std::string(prefix) + std::string(McDir);
    mHistogramRegistry->add(mcDir + getHistNameV2(kTrueQ3VsQ3, HistTable), getHistDesc(kTrueQ3VsQ3, HistTable), getHistType(kTrueQ3VsQ3, HistTable), {Specs.at(kTrueQ3VsQ3)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kTrueMtVsMt, HistTable), getHistDesc(kTrueMtVsMt, HistTable), getHistType(kTrueMtVsMt, HistTable), {Specs.at(kTrueMtVsMt)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kTrueMultVsMult, HistTable), getHistDesc(kTrueMultVsMult, HistTable), getHistType(kTrueMultVsMult, HistTable), {Specs.at(kTrueMultVsMult)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kTrueCentVsCent, HistTable), getHistDesc(kTrueCentVsCent, HistTable), getHistType(kTrueCentVsCent, HistTable), {Specs.at(kTrueCentVsCent)});
  }

  void fillAnalysis()
  {
    if (mPlot1d) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kQ3, HistTable)), mQ3);
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kMt, HistTable)), mMt);
    }
    if (mPlot2d) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kPt1VsQ3, HistTable)), mParticle1.Pt(), mQ3);
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kPt2VsQ3, HistTable)), mParticle2.Pt(), mQ3);
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kPt3VsQ3, HistTable)), mParticle3.Pt(), mQ3);
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kQ3VsMt, HistTable)), mQ3, mMt);
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kQ3VsMult, HistTable)), mQ3, mMult);
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kQ3VsCent, HistTable)), mQ3, mCent);
    }

    if (mPlotPt1VsPt2VsPt3) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kPt1VsPt2VsPt3, HistTable)), mParticle1.Pt(), mParticle2.Pt(), mParticle3.Pt());
    }
    if (mPlotQ3VsPt1VsPt2VsPt3) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kQ3VsPt1VsPt2VsPt3, HistTable)), mQ3, mParticle1.Pt(), mParticle2.Pt(), mParticle3.Pt());
    }
    if (mPlotQ3VsMtVsMult) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kQ3VsMtVsMult, HistTable)), mQ3, mMt, mMult);
    }
    if (mPlotQ3VsMtVsMultVsCent) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kQ3VsMtVsMultVsCent, HistTable)), mQ3, mMt, mMult, mCent);
    }
    if (mPlotQ3VsMtVsPt1VsPt2VsPt3VsMult) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kQ3VsMtVsPt1VsPt2VsPt3VsMult, HistTable)), mQ3, mMt, mParticle1.Pt(), mParticle2.Pt(), mParticle3.Pt(), mMult);
    }
    if (mPlotQ3VsMtVsPt1VsPt2VsPt3VsMultVsCent) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kQ3VsMtVsPt1VsPt2VsPt3VsMultVsCent, HistTable)), mQ3, mMt, mParticle1.Pt(), mParticle2.Pt(), mParticle3.Pt(), mMult, mCent);
    }
  }

  void fillMc()
  {
    if (mHasMcTriplet) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kTrueQ3VsQ3, HistTable)), mTrueQ3, mQ3);
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kTrueMtVsMt, HistTable)), mTrueMt, mMt);
    }
    if (mHasMcCol) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kTrueMultVsMult, HistTable)), mTrueMult, mMult);
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kTrueCentVsCent, HistTable)), mTrueCent, mCent);
    }
  }

  float getMt(ROOT::Math::PtEtaPhiMVector const& part1,
              float pdgMass1,
              ROOT::Math::PtEtaPhiMVector const& part2,
              float pdgMass2)
  {
    auto sum = part1 + part2;
    float mt = 0;
    float averageMass = 0;
    float reducedMass = 0;
    switch (mMtType) {
      case modes::TransverseMassType::kAveragePdgMass:
        averageMass = 0.5 * (pdgMass1 + pdgMass2);
        mt = std::hypot(0.5 * sum.Pt(), averageMass);
        break;
      case modes::TransverseMassType::kReducedPdgMass:
        reducedMass = 2.f * (pdgMass1 * pdgMass2 / (pdgMass1 + pdgMass2));
        mt = std::hypot(0.5 * sum.Pt(), reducedMass);
        break;
      case modes::TransverseMassType::kMt4Vector:
        mt = 0.5 * sum.Mt();
        break;
      default:
        LOG(fatal) << "Invalid transverse mass type, breaking...";
    }
    return mt;
  }

  float getMt(ROOT::Math::PtEtaPhiMVector const& part1,
              float mass1,
              ROOT::Math::PtEtaPhiMVector const& part2,
              float mass2,
              ROOT::Math::PtEtaPhiMVector const& part3,
              float mass3)
  {
    return (getMt(part1, mass1, part2, mass2) + getMt(part2, mass2, part3, mass3) + getMt(part1, mass1, part3, mass3)) / 3.;
  }

  o2::framework::HistogramRegistry* mHistogramRegistry = nullptr;
  float mPdgMass1 = 0.f;
  float mPdgMass2 = 0.f;
  float mPdgMass3 = 0.f;

  modes::TransverseMassType mMtType = modes::TransverseMassType::kAveragePdgMass;

  int mAbsCharge1 = 1;
  int mAbsCharge2 = 1;
  int mAbsCharge3 = 1;
  ROOT::Math::PtEtaPhiMVector mParticle1{};
  ROOT::Math::PtEtaPhiMVector mParticle2{};
  ROOT::Math::PtEtaPhiMVector mParticle3{};
  float mMass1 = 0.f;
  float mMass2 = 0.f;
  float mMass3 = 0.f;
  float mQ3 = 0.f;
  float mMt = 0.f;
  float mMult = 0.f;
  float mCent = 0.f;

  // mc
  ROOT::Math::PtEtaPhiMVector mTrueParticle1{};
  ROOT::Math::PtEtaPhiMVector mTrueParticle2{};
  ROOT::Math::PtEtaPhiMVector mTrueParticle3{};
  float mTrueQ3 = 0.f;
  float mTrueMt = 0.f;
  float mTrueMult = 0.f;
  float mTrueCent = 0.f;

  // cuts
  bool mHasMcTriplet = false;
  bool mHasMcCol = false;
  float mQ3Min = -1.f;
  float mQ3Max = -1.f;
  float mKtMin = -1.f;
  float mKtMax = -1.f;
  float mMtMin = -1.f;
  float mMtMax = -1.f;

  // flags
  bool mPlot1d = true;
  bool mPlot2d = true;

  bool mPlotPt1VsPt2VsPt3 = false;
  bool mPlotQ3VsPt1VsPt2VsPt3 = false;
  bool mPlotQ3VsMtVsMult = false;
  bool mPlotQ3VsMtVsMultVsCent = false;
  bool mPlotQ3VsMtVsPt1VsPt2VsPt3VsMult = false;
  bool mPlotQ3VsMtVsPt1VsPt2VsPt3VsMultVsCent = false;
};

}; // namespace triplethistmanager
}; // namespace o2::analysis::femto
#endif // PWGCF_FEMTO_CORE_TRIPLETHISTMANAGER_H_
