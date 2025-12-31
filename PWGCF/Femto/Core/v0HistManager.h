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

/// \file v0HistManager.h
/// \brief histogram manager for vzero histograms
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#ifndef PWGCF_FEMTO_CORE_V0HISTMANAGER_H_
#define PWGCF_FEMTO_CORE_V0HISTMANAGER_H_

#include "PWGCF/Femto/Core/histManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/trackHistManager.h"

#include "CommonConstants/MathConstants.h"
#include "Framework/Configurable.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/HistogramSpec.h"

#include <array>
#include <map>
#include <string>
#include <string_view>
#include <vector>

namespace o2::analysis::femto
{
namespace v0histmanager
{
// enum for track histograms
enum V0Hist {
  // analysis
  kPt,
  kEta,
  kPhi,
  kMass,
  kSign,
  // qa variables
  kMassLambda,
  kMassAntiLambda,
  kMassK0short,
  kCosPa,
  kDecayDauDca,
  kDecayVtxX,
  kDecayVtxY,
  kDecayVtxZ,
  kDecayVtx,
  kTransRadius,
  // 2d qa
  kPtVsEta,
  kPtVsPhi,
  kPhiVsEta,
  kPtVsCosPa,
  kPtVsLambdaMass,
  kPtVsAntiLambdaMass,
  kPtVsK0shortMass,
  kLambdaMassVsAntiLambdaMass,
  kK0shortMassVsLambdaMass,
  kK0shortMassVsAntiLambdaMass,
  // mc
  kOrigin,
  kPdg,
  kPdgMother,
  kPdgPartonicMother,
  kTruePt,
  kTrueEta,
  kTruePhi,
  // histograms for fraction estimation of v0s
  kNoMcParticle,
  kPrimary,
  kFromWrongCollision,
  kFromMaterial,
  kSecondary1,
  kSecondary2,
  kSecondary3,
  kSecondaryOther,

  kV0HistLast
};

constexpr std::size_t MaxSecondary = 3;

#define V0_DEFAULT_BINNING(defaultMassMin, defaultMassMax)                                         \
  o2::framework::ConfigurableAxis pt{"pt", {{600, 0, 6}}, "Pt"};                                   \
  o2::framework::ConfigurableAxis eta{"eta", {{300, -1.5, 1.5}}, "Eta"};                           \
  o2::framework::ConfigurableAxis phi{"phi", {{720, 0, 1.f * o2::constants::math::TwoPI}}, "Phi"}; \
  o2::framework::ConfigurableAxis mass{"mass", {{200, defaultMassMin, defaultMassMax}}, "Mass"};   \
  o2::framework::ConfigurableAxis sign{"sign", {{3, -1.5, 1.5}}, "Sign"};

template <const char* Prefix>
struct ConfLambdaBinning : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  V0_DEFAULT_BINNING(1.0, 1.2)
};
template <const char* Prefix>
struct ConfK0shortBinning : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  V0_DEFAULT_BINNING(0.475, 0.515)
};
#undef V0_DEFAULT_BINNING

constexpr const char PrefixLambdaBinning1[] = "LambdaBinning1";
using ConfLambdaBinning1 = ConfLambdaBinning<PrefixLambdaBinning1>;
constexpr const char PrefixK0shortBinning1[] = "K0shortBinning1";
using ConfK0shortBinning1 = ConfK0shortBinning<PrefixK0shortBinning1>;

template <const char* Prefix>
struct ConfV0QaBinning : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  o2::framework::Configurable<bool> plot2d{"plot2d", true, "Generate various 2D QA plots"};
  o2::framework::ConfigurableAxis cosPa{"cosPa", {{100, 0.9, 1}}, "Cosine of poiting angle"};
  o2::framework::ConfigurableAxis dauDcaAtDecay{"dauDcaAtDecay", {{150, 0, 1.5}}, "Daughter DCA at decay vertex"};
  o2::framework::ConfigurableAxis decayVertex{"decayVertex", {{100, 0, 100}}, "Decay vertex"};
  o2::framework::ConfigurableAxis transRadius{"transRadius", {{100, 0, 100}}, "Transverse radius"};
  o2::framework::ConfigurableAxis massLambda{"massLambda", {{200, 1, 1.2}}, "mass for antiparticle hypothesis"};
  o2::framework::ConfigurableAxis massAntiLambda{"massAntiLambda", {{100, 1, 1.2}}, "mass for antiparticle hypothesis"};
  o2::framework::ConfigurableAxis massK0short{"massK0short", {{200, 0.45, 0.55}}, "Mass for k0short hypothesis"};
};

constexpr const char PrefixLambdaQaBinning1[] = "LambdaQaBinning1";
using ConfLambdaQaBinning1 = ConfV0QaBinning<PrefixLambdaQaBinning1>;

constexpr const char PrefixK0shortQaBinning1[] = "K0shortQaBinning1";
using ConfK0shortQaBinning1 = ConfV0QaBinning<PrefixK0shortQaBinning1>;

template <const char* Prefix>
struct ConfV0McBinning : o2::framework::ConfigurableGroup {
  std::string prefix = std::string(Prefix);
  o2::framework::ConfigurableAxis pdgCodes{"pdgCodes", {{8001, -4000.5, 4000.5}}, "PDG codes of selected V0s"};
  o2::framework::ConfigurableAxis statusCode{"statusCode", {{21, -0.5, 20.5}}, "Status codes (i.e. Origin)"};
  o2::framework::Configurable<bool> plotOrigins{"plotOrigins", true, "Plot pt vs cosPa for different particle origins"};
  o2::framework::Configurable<std::vector<int>> pdgCodesForMothersOfSecondary{"pdgCodesForMothersOfSecondary", {3312, 3334}, "PDG codes of mothers of secondaries (Max 3 will be considered)"};
  o2::framework::ConfigurableAxis pt{"pt", {{150, 0, 3}}, "Pt"};
  o2::framework::ConfigurableAxis cosPa{"cosPa", {{100, 0.8, 1}}, "Cosine pointing angle"};
};

constexpr const char PrefixLambdaMcBinning1[] = "LambdaMcBinning1";
using ConfLambdaMcBinning1 = ConfV0McBinning<PrefixLambdaMcBinning1>;

constexpr const char PrefixK0shortMcBinning1[] = "K0shortMcBinning1";
using ConfK0shortMcBinning1 = ConfV0McBinning<PrefixK0shortMcBinning1>;

// must be in sync with enum TrackVariables
// the enum gives the correct index in the array
constexpr std::array<histmanager::HistInfo<V0Hist>, kV0HistLast> HistTable = {
  {{kPt, o2::framework::kTH1F, "hPt", "Transverse Momentum; p_{T} (GeV/#it{c}); Entries"},
   {kEta, o2::framework::kTH1F, "hEta", "Pseudorapdity; #eta; Entries"},
   {kPhi, o2::framework::kTH1F, "hPhi", "Azimuthal angle; #varphi; Entries"},
   {kMass, o2::framework::kTH1F, "hMass", "Invariant Mass; m_{Inv} (GeV/#it{c}^{2}); Entries"},
   {kSign, o2::framework::kTH1F, "hSign", "Sign (-1 -> antiparticle, 0 -> self conjugate, +1 -> particle); sign; Entries"},
   {kMassLambda, o2::framework::kTH1F, "hMassLambda", "#Lambda mass; m_{p#pi^{-}} (GeV/#it{c}^{2}); Entries"},
   {kMassAntiLambda, o2::framework::kTH1F, "hMassAntiLambda", "#bar{#Lambda} mass; m_{#bar{p}#pi^{+}} (GeV/#it{c}^{2}); Entries"},
   {kMassK0short, o2::framework::kTH1F, "hMassK0short", "K^{0}_{s} mass; m_{#pi^{+}#pi^{-}} (GeV/#it{c}^{2}); Entries"},
   {kCosPa, o2::framework::kTH1F, "hCosPa", "Cosine of pointing angle; cos(#alpha); Entries"},
   {kDecayDauDca, o2::framework::kTH1F, "hDauDca", "Daughter DCA at decay vertex ; DCA_{Decay vertex} (cm); Entries"},
   {kDecayVtxX, o2::framework::kTH1F, "hDecayVtxX", "X coordinate of decay vertex ; DV_{X} (cm); Entries"},
   {kDecayVtxY, o2::framework::kTH1F, "hDecayVtxY", "Y coordinate of decay vertex ; DV_{Y} (cm); Entries"},
   {kDecayVtxZ, o2::framework::kTH1F, "hDecayVtxZ", "Z coordinate of decay vertex ; DV_{Z} (cm); Entries"},
   {kDecayVtx, o2::framework::kTH1F, "hDecayVtx", "Distance of decay vertex from primary vertex ; DV (cm); Entries"},
   {kTransRadius, o2::framework::kTH1F, "hTransRadius", "Transverse radius ; r_{xy} (cm); Entries"},
   {kPtVsEta, o2::framework::kTH2F, "hPtVsEta", "p_{T} vs #eta; p_{T} (GeV/#it{c}) ; #eta"},
   {kPtVsPhi, o2::framework::kTH2F, "hPtVsPhi", "p_{T} vs #varphi; p_{T} (GeV/#it{c}) ; #varphi"},
   {kPhiVsEta, o2::framework::kTH2F, "hPhiVsEta", "#varphi vs #eta; #varphi ; #eta"},
   {kPtVsCosPa, o2::framework::kTH2F, "hPtVsCosPa", "Cosine of poiting angle vs p_{T}; cos(#alpha); p_{T} (GeV/#it{c})"},
   {kPtVsLambdaMass, o2::framework::kTH2F, "hPtVsLambdaMass", "p_{T} vs #Lambda mass; p_{T} (GeV/#it{c}); m_{p#pi^{-}} (GeV/#it{c}^{2})"},
   {kPtVsAntiLambdaMass, o2::framework::kTH2F, "hPtVsAntiLambdaMass", "p_{T} vs #bar{#Lambda} mass; p_{T} (GeV/#it{c}); m_{#bar{p}#pi^{+}} (GeV/#it{c}^{2})"},
   {kPtVsK0shortMass, o2::framework::kTH2F, "hPtVsK0shortMass", "p_{T} vs K^{0}_{S} mass; p_{T} (GeV/#it{c}); m_{#pi^{+}#pi^{-}} (GeV/#it{c}^{2})"},
   {kK0shortMassVsLambdaMass, o2::framework::kTH2F, "hK0shortMassVsLambdaMass", " K^{0}_{S} mass vs #Lambda mass; m_{#pi^{+}#pi^{-}} (GeV/#it{c}^{2}); m_{p#pi^{-}} (GeV/#it{c}^{2})"},
   {kK0shortMassVsAntiLambdaMass, o2::framework::kTH2F, "hK0shortMassVsAntiLambdaMass", "K^{0}_{S} mass vs #bar{#Lambda} mass; m_{#pi^{+}#pi^{-}} (GeV/#it{c}^{2}); m_{#bar{p}#pi^{+}} (GeV/#it{c}^{2})"},
   {kLambdaMassVsAntiLambdaMass, o2::framework::kTH2F, "hLambdaMassVsAntiLambdaMass", "#Lambda mass vs #bar{#Lambda}; m_{p#pi^{-}} (GeV/#it{c}^{2}); m_{#bar{p}#pi^{+}} (GeV/#it{c}^{2})"},
   {kOrigin, o2::framework::kTH1F, "hOrigin", "Status Codes (=Origin); Status Code; Entries"},
   {kPdg, o2::framework::kTH1F, "hPdg", "PDG Codes of reconstructed v0; PDG Code; Entries"},
   {kPdgMother, o2::framework::kTH1F, "hPdgMother", "PDG Codes of mother of reconstructed v0; PDG Code; Entries"},
   {kPdgPartonicMother, o2::framework::kTH1F, "hPdgPartonicMother", "PDG Codes of partonic mother of reconstructed v0; PDG Code; Entries"},
   {kTruePt, o2::framework::kTH1F, "hTruePt", "True transverse momentum; p_{T} (GeV/#it{c}); Entries"},
   {kTrueEta, o2::framework::kTH1F, "hTrueEta", "True pseudorapdity; #eta; Entries"},
   {kTruePhi, o2::framework::kTH1F, "hTruePhi", "True azimuthal angle; #varphi; Entries"},
   {kNoMcParticle, o2::framework::kTH2F, "hNoMcParticle", "Wrongly reconstructed particles; p_{T} (GeV/#it{c}); cos(#alpha)"},
   {kPrimary, o2::framework::kTH2F, "hPrimary", "Primary particles; p_{T} (GeV/#it{c}); cos(#alpha)"},
   {kFromWrongCollision, o2::framework::kTH2F, "hFromWrongCollision", "Particles associated to wrong collision; p_{T} (GeV/#it{c}); cos(#alpha)"},
   {kFromMaterial, o2::framework::kTH2F, "hFromMaterial", "Particles from material; p_{T} (GeV/#it{c}); cos(#alpha)"},
   {kSecondary1, o2::framework::kTH2F, "hFromSecondary1", "Particles from secondary decay; p_{T} (GeV/#it{c}); cos(#alpha)"},
   {kSecondary2, o2::framework::kTH2F, "hFromSecondary2", "Particles from seconary decay; p_{T} (GeV/#it{c}); cos(#alpha)"},
   {kSecondary3, o2::framework::kTH2F, "hFromSecondary3", "Particles from seconary decay; p_{T} (GeV/#it{c}); cos(#alpha)"},
   {kSecondaryOther, o2::framework::kTH2F, "hFromSecondaryOther", "Particles from every other seconary decay; p_{T} (GeV/#it{c}); cos(#alpha)"}},
};

#define V0_HIST_ANALYSIS_MAP(conf) \
  {kPt, {conf.pt}},                \
    {kEta, {conf.eta}},            \
    {kPhi, {conf.phi}},            \
    {kMass, {conf.mass}},          \
    {kSign, {conf.sign}},

#define V0_HIST_QA_MAP(confAnalysis, confQa)                                   \
  {kCosPa, {confQa.cosPa}},                                                    \
    {kDecayDauDca, {confQa.dauDcaAtDecay}},                                    \
    {kDecayVtxX, {confQa.decayVertex}},                                        \
    {kDecayVtxY, {confQa.decayVertex}},                                        \
    {kDecayVtxZ, {confQa.decayVertex}},                                        \
    {kDecayVtx, {confQa.decayVertex}},                                         \
    {kTransRadius, {confQa.transRadius}},                                      \
    {kPtVsEta, {confAnalysis.pt, confAnalysis.eta}},                           \
    {kPtVsPhi, {confAnalysis.pt, confAnalysis.phi}},                           \
    {kPhiVsEta, {confAnalysis.phi, confAnalysis.eta}},                         \
    {kPtVsCosPa, {confAnalysis.pt, confQa.cosPa}},                             \
    {kMassLambda, {confQa.massLambda}},                                        \
    {kMassAntiLambda, {confQa.massAntiLambda}},                                \
    {kMassK0short, {confQa.massK0short}},                                      \
    {kPtVsLambdaMass, {confAnalysis.pt, confQa.massLambda}},                   \
    {kPtVsAntiLambdaMass, {confAnalysis.pt, confQa.massAntiLambda}},           \
    {kPtVsK0shortMass, {confAnalysis.pt, confQa.massK0short}},                 \
    {kLambdaMassVsAntiLambdaMass, {confQa.massLambda, confQa.massAntiLambda}}, \
    {kK0shortMassVsLambdaMass, {confQa.massK0short, confQa.massLambda}},       \
    {kK0shortMassVsAntiLambdaMass, {confQa.massK0short, confQa.massAntiLambda}},

#define V0_HIST_MC_MAP(confAnalysis, confMc)          \
  {kTruePt, {confAnalysis.pt}},                       \
    {kTrueEta, {confAnalysis.eta}},                   \
    {kTruePhi, {confAnalysis.phi}},                   \
    {kOrigin, {confMc.statusCode}},                   \
    {kPdg, {confMc.pdgCodes}},                        \
    {kPdgMother, {confMc.pdgCodes}},                  \
    {kPdgPartonicMother, {confMc.pdgCodes}},          \
    {kNoMcParticle, {confMc.pt, confMc.cosPa}},       \
    {kPrimary, {confMc.pt, confMc.cosPa}},            \
    {kFromWrongCollision, {confMc.pt, confMc.cosPa}}, \
    {kFromMaterial, {confMc.pt, confMc.cosPa}},       \
    {kSecondary1, {confMc.pt, confMc.cosPa}},         \
    {kSecondary2, {confMc.pt, confMc.cosPa}},         \
    {kSecondary3, {confMc.pt, confMc.cosPa}},         \
    {kSecondaryOther, {confMc.pt, confMc.cosPa}},

template <typename T>
auto makeV0HistSpecMap(const T& confBinningAnalysis)
{
  return std::map<V0Hist, std::vector<framework::AxisSpec>>{
    V0_HIST_ANALYSIS_MAP(confBinningAnalysis)};
}

template <typename T1, typename T2>
std::map<V0Hist, std::vector<framework::AxisSpec>> makeV0QaHistSpecMap(T1 const& confBinningAnalysis, T2 const& confBinningQa)
{
  return std::map<V0Hist, std::vector<framework::AxisSpec>>{
    V0_HIST_ANALYSIS_MAP(confBinningAnalysis)
      V0_HIST_QA_MAP(confBinningAnalysis, confBinningQa)};
}

template <typename T1, typename T2, typename T3>
std::map<V0Hist, std::vector<framework::AxisSpec>> makeV0McQaHistSpecMap(T1 const& confBinningAnalysis, T2 const& confBinningQa, T3 const& confBinningMc)
{
  return std::map<V0Hist, std::vector<framework::AxisSpec>>{
    V0_HIST_ANALYSIS_MAP(confBinningAnalysis)
      V0_HIST_QA_MAP(confBinningAnalysis, confBinningQa)
        V0_HIST_MC_MAP(confBinningAnalysis, confBinningMc)};
}

#undef V0_HIST_ANALYSIS_MAP
#undef V0_HIST_QA_MAP
#undef V0_HIST_MC_MAP

constexpr char PrefixLambdaQa[] = "LambdaQA/";
constexpr char PrefixLambda1[] = "Lambda1/";
constexpr char PrefixLambda2[] = "Lambda2/";
constexpr char PrefixK0shortQa[] = "K0shortQa/";
constexpr char PrefixK0short1[] = "K0short1/";
constexpr char PrefixK0short2[] = "K0short2/";

constexpr char PrefixLambdaCascade[] = "LambdaCascadeQa/";

constexpr std::string_view AnalysisDir = "Kinematics/";
constexpr std::string_view QaDir = "QA/";
constexpr std::string_view McDir = "MC/";

constexpr int AbsChargeDaughters = 1;

/// \class FemtoDreamEventHisto
/// \brief Class for histogramming event properties
// template <femtomodes::Mode mode>
template <const char* v0Prefix,
          const char* posDauPrefix,
          const char* negDauPrefix,
          modes::V0 v0>
class V0HistManager
{
 public:
  V0HistManager() = default;
  ~V0HistManager() = default;

  // init for analysis
  template <modes::Mode mode>
  void init(o2::framework::HistogramRegistry* registry,
            std::map<V0Hist, std::vector<o2::framework::AxisSpec>> const& V0Specs,
            std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> const& PosDauSpecs,
            std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> const& NegDauSpecs)
  {
    mHistogramRegistry = registry;
    mPosDauManager.template init<mode>(registry, PosDauSpecs, AbsChargeDaughters);
    mNegDauManager.template init<mode>(registry, NegDauSpecs, AbsChargeDaughters);

    if constexpr (modes::isFlagSet(mode, modes::Mode::kAnalysis)) {
      initAnalysis(V0Specs);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kQa)) {
      initQa(V0Specs);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kMc)) {
      initMc(V0Specs);
    }
  }

  // init for qa
  template <modes::Mode mode, typename T1, typename T2, typename T3>
  void init(o2::framework::HistogramRegistry* registry,
            std::map<V0Hist, std::vector<o2::framework::AxisSpec>> const& V0Specs,
            T1 const& V0ConfBinningQa,
            std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> const& PosDauSpecs,
            T2 const& PosDauConfBinningQa,
            std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> const& NegDauSpecs,
            T3 const& NegDauConfBinningQa)
  {
    mHistogramRegistry = registry;
    mPosDauManager.template init<mode>(registry, PosDauSpecs, PosDauConfBinningQa, AbsChargeDaughters);
    mNegDauManager.template init<mode>(registry, NegDauSpecs, NegDauConfBinningQa, AbsChargeDaughters);
    this->enableOptionalHistograms(V0ConfBinningQa);

    if constexpr (modes::isFlagSet(mode, modes::Mode::kAnalysis)) {
      initAnalysis(V0Specs);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kQa)) {
      initQa(V0Specs);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kMc)) {
      initMc(V0Specs);
    }
  }

  // init for qa + mc
  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  void init(o2::framework::HistogramRegistry* registry,
            std::map<V0Hist, std::vector<o2::framework::AxisSpec>> const& V0Specs,
            T1 const& V0ConfBinningQa,
            T2 const& V0ConfBinningMc,
            std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> const& PosDauSpecs,
            T3 const& PosDauConfBinningQa,
            T4 const& PosDauConfBinningMc,
            std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> const& NegDauSpecs,
            T5 const& NegDauConfBinningQa,
            T6 const& NegDauConfBinningMc)
  {
    mHistogramRegistry = registry;
    mPosDauManager.template init<mode>(registry, PosDauSpecs, PosDauConfBinningQa, PosDauConfBinningMc, AbsChargeDaughters);
    mNegDauManager.template init<mode>(registry, NegDauSpecs, NegDauConfBinningQa, NegDauConfBinningMc, AbsChargeDaughters);
    this->enableOptionalHistograms(V0ConfBinningQa, V0ConfBinningMc);

    if constexpr (modes::isFlagSet(mode, modes::Mode::kAnalysis)) {
      initAnalysis(V0Specs);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kQa)) {
      initQa(V0Specs);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kMc)) {
      initMc(V0Specs);
    }
  }

  template <modes::Mode mode, typename T1, typename T2>
  void fill(T1 const& v0candidate, T2 const& tracks)
  {
    auto posDaughter = tracks.rawIteratorAt(v0candidate.posDauId() - tracks.offset());
    mPosDauManager.template fill<mode>(posDaughter, tracks);
    auto negDaughter = tracks.rawIteratorAt(v0candidate.negDauId() - tracks.offset());
    mNegDauManager.template fill<mode>(negDaughter, tracks);

    if constexpr (modes::isFlagSet(mode, modes::Mode::kAnalysis)) {
      fillAnalysis(v0candidate);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kQa)) {
      fillQa(v0candidate);
    }
  }

  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5>
  void fill(T1 const& v0candidate, T2 const& tracks, T3 const& mcParticles, T4 const& mcMothers, T5 const& mcPartonicMothers)
  {
    auto posDaughter = tracks.rawIteratorAt(v0candidate.posDauId() - tracks.offset());
    mPosDauManager.template fill<mode>(posDaughter, tracks, mcParticles, mcMothers, mcPartonicMothers);
    auto negDaughter = tracks.rawIteratorAt(v0candidate.negDauId() - tracks.offset());
    mNegDauManager.template fill<mode>(negDaughter, tracks, mcParticles, mcMothers, mcPartonicMothers);
    if constexpr (modes::isFlagSet(mode, modes::Mode::kAnalysis)) {
      fillAnalysis(v0candidate);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kQa)) {
      fillQa(v0candidate);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kMc)) {
      fillMc(v0candidate, mcParticles, mcMothers, mcPartonicMothers);
    }
  }

 private:
  // for qa
  template <typename T1>
  void enableOptionalHistograms(T1 const& V0ConfBinningQa)
  {
    mPlot2d = V0ConfBinningQa.plot2d.value;
  }

  // for qa and mc
  template <typename T1, typename T2>
  void enableOptionalHistograms(T1 const& V0ConfBinningQa, T2 const& V0ConfBinningMc)
  {
    this->enableOptionalHistograms(V0ConfBinningQa);

    mPlotOrigins = V0ConfBinningMc.plotOrigins.value;
    mPlotNSecondaries = V0ConfBinningMc.pdgCodesForMothersOfSecondary.value.size();

    for (std::size_t i = 0; i < MaxSecondary; i++) {
      if (i < V0ConfBinningMc.pdgCodesForMothersOfSecondary.value.size()) {
        mPdgCodesSecondaryMother.at(i) = std::abs(V0ConfBinningMc.pdgCodesForMothersOfSecondary.value.at(i));
      } else {
        mPdgCodesSecondaryMother.at(i) = 0;
      }
    }
  }

  void initAnalysis(std::map<V0Hist, std::vector<o2::framework::AxisSpec>> const& V0Specs)
  {
    std::string analysisDir = std::string(v0Prefix) + std::string(AnalysisDir);
    mHistogramRegistry->add(analysisDir + getHistNameV2(kPt, HistTable), getHistDesc(kPt, HistTable), getHistType(kPt, HistTable), {V0Specs.at(kPt)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kEta, HistTable), getHistDesc(kEta, HistTable), getHistType(kEta, HistTable), {V0Specs.at(kEta)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kPhi, HistTable), getHistDesc(kPhi, HistTable), getHistType(kPhi, HistTable), {V0Specs.at(kPhi)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kMass, HistTable), getHistDesc(kMass, HistTable), getHistType(kMass, HistTable), {V0Specs.at(kMass)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kSign, HistTable), getHistDesc(kSign, HistTable), getHistType(kSign, HistTable), {V0Specs.at(kSign)});
  }

  void initQa(std::map<V0Hist, std::vector<o2::framework::AxisSpec>> const& V0Specs)
  {
    std::string qaDir = std::string(v0Prefix) + std::string(QaDir);

    mHistogramRegistry->add(qaDir + getHistNameV2(kCosPa, HistTable), getHistDesc(kCosPa, HistTable), getHistType(kCosPa, HistTable), {V0Specs.at(kCosPa)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kDecayDauDca, HistTable), getHistDesc(kDecayDauDca, HistTable), getHistType(kDecayDauDca, HistTable), {V0Specs.at(kDecayDauDca)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kDecayVtxX, HistTable), getHistDesc(kDecayVtxX, HistTable), getHistType(kDecayVtxX, HistTable), {V0Specs.at(kDecayVtxX)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kDecayVtxY, HistTable), getHistDesc(kDecayVtxY, HistTable), getHistType(kDecayVtxY, HistTable), {V0Specs.at(kDecayVtxY)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kDecayVtxZ, HistTable), getHistDesc(kDecayVtxZ, HistTable), getHistType(kDecayVtxZ, HistTable), {V0Specs.at(kDecayVtxZ)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kDecayVtx, HistTable), getHistDesc(kDecayVtx, HistTable), getHistType(kDecayVtx, HistTable), {V0Specs.at(kDecayVtx)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kTransRadius, HistTable), getHistDesc(kTransRadius, HistTable), getHistType(kTransRadius, HistTable), {V0Specs.at(kTransRadius)});

    if (mPlot2d) {
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsEta, HistTable), getHistDesc(kPtVsEta, HistTable), getHistType(kPtVsEta, HistTable), {V0Specs.at(kPtVsEta)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsPhi, HistTable), getHistDesc(kPtVsPhi, HistTable), getHistType(kPtVsPhi, HistTable), {V0Specs.at(kPtVsPhi)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPhiVsEta, HistTable), getHistDesc(kPhiVsEta, HistTable), getHistType(kPhiVsEta, HistTable), {V0Specs.at(kPhiVsEta)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsCosPa, HistTable), getHistDesc(kPtVsCosPa, HistTable), getHistType(kPtVsCosPa, HistTable), {V0Specs.at(kPtVsCosPa)});

      mHistogramRegistry->add(qaDir + getHistNameV2(kMassLambda, HistTable), getHistDesc(kMassLambda, HistTable), getHistType(kMassLambda, HistTable), {V0Specs.at(kMassLambda)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kMassAntiLambda, HistTable), getHistDesc(kMassAntiLambda, HistTable), getHistType(kMassAntiLambda, HistTable), {V0Specs.at(kMassAntiLambda)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kMassK0short, HistTable), getHistDesc(kMassK0short, HistTable), getHistType(kMassK0short, HistTable), {V0Specs.at(kMassK0short)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsLambdaMass, HistTable), getHistDesc(kPtVsLambdaMass, HistTable), getHistType(kPtVsLambdaMass, HistTable), {V0Specs.at(kPtVsLambdaMass)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsAntiLambdaMass, HistTable), getHistDesc(kPtVsAntiLambdaMass, HistTable), getHistType(kPtVsAntiLambdaMass, HistTable), {V0Specs.at(kPtVsAntiLambdaMass)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsK0shortMass, HistTable), getHistDesc(kPtVsK0shortMass, HistTable), getHistType(kPtVsK0shortMass, HistTable), {V0Specs.at(kPtVsK0shortMass)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kLambdaMassVsAntiLambdaMass, HistTable), getHistDesc(kLambdaMassVsAntiLambdaMass, HistTable), getHistType(kLambdaMassVsAntiLambdaMass, HistTable), {V0Specs.at(kLambdaMassVsAntiLambdaMass)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kK0shortMassVsLambdaMass, HistTable), getHistDesc(kK0shortMassVsLambdaMass, HistTable), getHistType(kK0shortMassVsLambdaMass, HistTable), {V0Specs.at(kK0shortMassVsLambdaMass)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kK0shortMassVsAntiLambdaMass, HistTable), getHistDesc(kK0shortMassVsAntiLambdaMass, HistTable), getHistType(kK0shortMassVsAntiLambdaMass, HistTable), {V0Specs.at(kK0shortMassVsAntiLambdaMass)});
    }
  }

  void initMc(std::map<V0Hist, std::vector<o2::framework::AxisSpec>> const& V0Specs)
  {
    std::string mcDir = std::string(v0Prefix) + std::string(McDir);
    mHistogramRegistry->add(mcDir + getHistNameV2(kTruePt, HistTable), getHistDesc(kTruePt, HistTable), getHistType(kTruePt, HistTable), {V0Specs.at(kTruePt)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kTrueEta, HistTable), getHistDesc(kTrueEta, HistTable), getHistType(kTrueEta, HistTable), {V0Specs.at(kTrueEta)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kTruePhi, HistTable), getHistDesc(kTruePhi, HistTable), getHistType(kTruePhi, HistTable), {V0Specs.at(kTruePhi)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kOrigin, HistTable), getHistDesc(kOrigin, HistTable), getHistType(kOrigin, HistTable), {V0Specs.at(kOrigin)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kPdg, HistTable), getHistDesc(kPdg, HistTable), getHistType(kPdg, HistTable), {V0Specs.at(kPdg)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kPdgMother, HistTable), getHistDesc(kPdgMother, HistTable), getHistType(kPdgMother, HistTable), {V0Specs.at(kPdgMother)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kPdgPartonicMother, HistTable), getHistDesc(kPdgPartonicMother, HistTable), getHistType(kPdgPartonicMother, HistTable), {V0Specs.at(kPdgPartonicMother)});

    if (mPlotOrigins) {
      mHistogramRegistry->add(mcDir + getHistNameV2(kNoMcParticle, HistTable), getHistDesc(kNoMcParticle, HistTable), getHistType(kNoMcParticle, HistTable), {V0Specs.at(kNoMcParticle)});
      mHistogramRegistry->add(mcDir + getHistNameV2(kPrimary, HistTable), getHistDesc(kPrimary, HistTable), getHistType(kPrimary, HistTable), {V0Specs.at(kPrimary)});
      mHistogramRegistry->add(mcDir + getHistNameV2(kFromWrongCollision, HistTable), getHistDesc(kFromWrongCollision, HistTable), getHistType(kFromWrongCollision, HistTable), {V0Specs.at(kFromWrongCollision)});
      mHistogramRegistry->add(mcDir + getHistNameV2(kFromMaterial, HistTable), getHistDesc(kFromMaterial, HistTable), getHistType(kFromMaterial, HistTable), {V0Specs.at(kFromMaterial)});

      if (mPlotNSecondaries >= histmanager::kSecondaryPlotLevel1) {
        mHistogramRegistry->add(mcDir + getHistNameV2(kSecondary1, HistTable), getHistDesc(kSecondary1, HistTable), getHistType(kSecondary1, HistTable), {V0Specs.at(kSecondary1)});
      }
      if (mPlotNSecondaries >= histmanager::kSecondaryPlotLevel2) {
        mHistogramRegistry->add(mcDir + getHistNameV2(kSecondary2, HistTable), getHistDesc(kSecondary2, HistTable), getHistType(kSecondary2, HistTable), {V0Specs.at(kSecondary2)});
      }
      if (mPlotNSecondaries >= histmanager::kSecondaryPlotLevel3) {
        mHistogramRegistry->add(mcDir + getHistNameV2(kSecondary3, HistTable), getHistDesc(kSecondary3, HistTable), getHistType(kSecondary3, HistTable), {V0Specs.at(kSecondary3)});
      }
      mHistogramRegistry->add(mcDir + getHistNameV2(kSecondaryOther, HistTable), getHistDesc(kSecondaryOther, HistTable), getHistType(kSecondaryOther, HistTable), {V0Specs.at(kSecondaryOther)});
    }
  }

  template <typename T>
  void fillAnalysis(T const& v0candidate)
  {
    mHistogramRegistry->fill(HIST(v0Prefix) + HIST(AnalysisDir) + HIST(getHistName(kPt, HistTable)), v0candidate.pt());
    mHistogramRegistry->fill(HIST(v0Prefix) + HIST(AnalysisDir) + HIST(getHistName(kEta, HistTable)), v0candidate.eta());
    mHistogramRegistry->fill(HIST(v0Prefix) + HIST(AnalysisDir) + HIST(getHistName(kPhi, HistTable)), v0candidate.phi());
    mHistogramRegistry->fill(HIST(v0Prefix) + HIST(AnalysisDir) + HIST(getHistName(kMass, HistTable)), v0candidate.mass());

    if constexpr (modes::isEqual(v0, modes::V0::kLambda) || modes::isEqual(v0, modes::V0::kAntiLambda)) {
      mHistogramRegistry->fill(HIST(v0Prefix) + HIST(AnalysisDir) + HIST(getHistName(kSign, HistTable)), v0candidate.sign());
    }
    if constexpr (modes::isEqual(v0, modes::V0::kK0short)) {
      mHistogramRegistry->fill(HIST(v0Prefix) + HIST(AnalysisDir) + HIST(getHistName(kSign, HistTable)), 0);
    }
  }

  template <typename T>
  void fillQa(T const& v0candidate)
  {
    mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kCosPa, HistTable)), v0candidate.cosPa());
    mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kDecayDauDca, HistTable)), v0candidate.dauDca());
    mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kDecayVtxX, HistTable)), v0candidate.decayVtxX());
    mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kDecayVtxY, HistTable)), v0candidate.decayVtxY());
    mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kDecayVtxZ, HistTable)), v0candidate.decayVtxZ());
    mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kDecayVtx, HistTable)), v0candidate.decayVtx());
    mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kTransRadius, HistTable)), v0candidate.transRadius());

    if (mPlot2d) {
      mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kPtVsEta, HistTable)), v0candidate.pt(), v0candidate.eta());
      mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kPtVsPhi, HistTable)), v0candidate.pt(), v0candidate.phi());
      mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kPhiVsEta, HistTable)), v0candidate.phi(), v0candidate.eta());
      mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kPtVsCosPa, HistTable)), v0candidate.pt(), v0candidate.cosPa());

      if constexpr (modes::isEqual(v0, modes::V0::kLambda) || modes::isEqual(v0, modes::V0::kAntiLambda)) {
        float massLambda, massAntiLambda;
        if (v0candidate.sign() > 0) {
          massLambda = v0candidate.mass();
          massAntiLambda = v0candidate.massAnti();
        } else {
          massLambda = v0candidate.massAnti();
          massAntiLambda = v0candidate.mass();
        }
        mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kMassLambda, HistTable)), massLambda);
        mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kMassAntiLambda, HistTable)), massAntiLambda);
        mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kMassK0short, HistTable)), v0candidate.massK0short());
        mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kPtVsLambdaMass, HistTable)), v0candidate.pt(), massLambda);
        mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kPtVsAntiLambdaMass, HistTable)), v0candidate.pt(), massAntiLambda);
        mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kPtVsK0shortMass, HistTable)), v0candidate.pt(), v0candidate.massK0short());
        mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kLambdaMassVsAntiLambdaMass, HistTable)), massLambda, massAntiLambda);
        mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kK0shortMassVsLambdaMass, HistTable)), v0candidate.massK0short(), massLambda);
        mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kK0shortMassVsAntiLambdaMass, HistTable)), v0candidate.massK0short(), massAntiLambda);
      }

      if constexpr (modes::isEqual(v0, modes::V0::kK0short)) {
        mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kMassLambda, HistTable)), v0candidate.massLambda());
        mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kMassAntiLambda, HistTable)), v0candidate.massAntiLambda());
        mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kMassK0short, HistTable)), v0candidate.mass());
        mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kPtVsLambdaMass, HistTable)), v0candidate.pt(), v0candidate.massLambda());
        mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kPtVsAntiLambdaMass, HistTable)), v0candidate.pt(), v0candidate.massAntiLambda());
        mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kPtVsK0shortMass, HistTable)), v0candidate.pt(), v0candidate.mass());
        mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kLambdaMassVsAntiLambdaMass, HistTable)), v0candidate.massLambda(), v0candidate.massAntiLambda());
        mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kK0shortMassVsLambdaMass, HistTable)), v0candidate.mass(), v0candidate.massLambda());
        mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kK0shortMassVsAntiLambdaMass, HistTable)), v0candidate.mass(), v0candidate.massAntiLambda());
      }
    }
  }

  template <typename T1, typename T2, typename T3, typename T4>
  void fillMc(T1 const& v0candidate, T2 const& /*mcParticles*/, T3 const& /*mcMothers*/, T4 const& /*mcPartonicMothers*/)
  {
    // No MC Particle
    if (!v0candidate.has_fMcParticle()) {
      mHistogramRegistry->fill(HIST(v0Prefix) + HIST(McDir) + HIST(getHistName(kPdg, HistTable)), 0);
      mHistogramRegistry->fill(HIST(v0Prefix) + HIST(McDir) + HIST(getHistName(kOrigin, HistTable)), static_cast<float>(modes::McOrigin::kNoMcParticle));
      if (mPlotOrigins) {
        mHistogramRegistry->fill(HIST(v0Prefix) + HIST(McDir) + HIST(getHistName(kNoMcParticle, HistTable)), v0candidate.pt(), v0candidate.cosPa());
      }
      return;
    }

    // Retrieve MC particle
    auto mcParticle = v0candidate.template fMcParticle_as<T2>();
    mHistogramRegistry->fill(HIST(v0Prefix) + HIST(McDir) + HIST(getHistName(kTruePt, HistTable)), mcParticle.pt());
    mHistogramRegistry->fill(HIST(v0Prefix) + HIST(McDir) + HIST(getHistName(kTrueEta, HistTable)), mcParticle.eta());
    mHistogramRegistry->fill(HIST(v0Prefix) + HIST(McDir) + HIST(getHistName(kTruePhi, HistTable)), mcParticle.phi());
    mHistogramRegistry->fill(HIST(v0Prefix) + HIST(McDir) + HIST(getHistName(kOrigin, HistTable)), mcParticle.origin());
    mHistogramRegistry->fill(HIST(v0Prefix) + HIST(McDir) + HIST(getHistName(kPdg, HistTable)), mcParticle.pdgCode());

    // Mother PDG
    if (v0candidate.has_fMcMother()) {
      auto mother = v0candidate.template fMcMother_as<T3>();
      mHistogramRegistry->fill(HIST(v0Prefix) + HIST(McDir) + HIST(getHistName(kPdgMother, HistTable)), mother.pdgCode());
    } else {
      mHistogramRegistry->fill(HIST(v0Prefix) + HIST(McDir) + HIST(getHistName(kPdgMother, HistTable)), 0);
    }

    // Partonic Mother PDG
    if (v0candidate.has_fMcPartMoth()) {
      auto partonicMother = v0candidate.template fMcPartMoth_as<T4>();
      mHistogramRegistry->fill(HIST(v0Prefix) + HIST(McDir) + HIST(getHistName(kPdgPartonicMother, HistTable)), partonicMother.pdgCode());
    } else {
      mHistogramRegistry->fill(HIST(v0Prefix) + HIST(McDir) + HIST(getHistName(kPdgPartonicMother, HistTable)), 0);
    }

    // Plot origins
    if (mPlotOrigins) {
      switch (static_cast<modes::McOrigin>(mcParticle.origin())) {
        case modes::McOrigin::kPhysicalPrimary:
          mHistogramRegistry->fill(HIST(v0Prefix) + HIST(McDir) + HIST(getHistName(kPrimary, HistTable)), v0candidate.pt(), v0candidate.cosPa());
          break;
        case modes::McOrigin::kFromWrongCollision:
          mHistogramRegistry->fill(HIST(v0Prefix) + HIST(McDir) + HIST(getHistName(kFromWrongCollision, HistTable)), v0candidate.pt(), v0candidate.cosPa());
          break;
        case modes::McOrigin::kFromMaterial:
          mHistogramRegistry->fill(HIST(v0Prefix) + HIST(McDir) + HIST(getHistName(kFromMaterial, HistTable)), v0candidate.pt(), v0candidate.cosPa());
          break;
        case modes::McOrigin::kFromSecondaryDecay:
          if (v0candidate.has_fMcMother()) {
            auto mother = v0candidate.template fMcMother_as<T3>();
            int motherPdgCode = std::abs(mother.pdgCode());
            // Switch on PDG of the mother
            if (mPlotNSecondaries >= histmanager::kSecondaryPlotLevel1 && motherPdgCode == mPdgCodesSecondaryMother[0]) {
              mHistogramRegistry->fill(HIST(v0Prefix) + HIST(McDir) + HIST(getHistName(kSecondary1, HistTable)), v0candidate.pt(), v0candidate.cosPa());
            } else if (mPlotNSecondaries >= histmanager::kSecondaryPlotLevel2 && motherPdgCode == mPdgCodesSecondaryMother[1]) {
              mHistogramRegistry->fill(HIST(v0Prefix) + HIST(McDir) + HIST(getHistName(kSecondary2, HistTable)), v0candidate.pt(), v0candidate.cosPa());
            } else if (mPlotNSecondaries >= histmanager::kSecondaryPlotLevel3 && motherPdgCode == mPdgCodesSecondaryMother[2]) {
              mHistogramRegistry->fill(HIST(v0Prefix) + HIST(McDir) + HIST(getHistName(kSecondary3, HistTable)), v0candidate.pt(), v0candidate.cosPa());
            } else {
              mHistogramRegistry->fill(HIST(v0Prefix) + HIST(McDir) + HIST(getHistName(kSecondaryOther, HistTable)), v0candidate.pt(), v0candidate.cosPa());
            }
          }
          break;
        default:
          // Unknown origin → safely ignore
          break;
      }
    }
  }

  o2::framework::HistogramRegistry* mHistogramRegistry = nullptr;
  bool mPlot2d = false;
  bool mPlotOrigins = false;
  int mPlotNSecondaries = 0;
  std::array<int, MaxSecondary> mPdgCodesSecondaryMother = {0};
  trackhistmanager::TrackHistManager<posDauPrefix> mPosDauManager;
  trackhistmanager::TrackHistManager<negDauPrefix> mNegDauManager;
};
}; // namespace v0histmanager
}; // namespace o2::analysis::femto
#endif // PWGCF_FEMTO_CORE_V0HISTMANAGER_H_
