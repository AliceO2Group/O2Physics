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
  kTruePtVsPt,
  kTrueEtaVsEta,
  kTruePhiVsPhi,
  // histograms for fraction estimation of v0s
  kNoMcParticle,
  kPrimary,
  kFromWrongCollision,
  kFromMaterial,
  kMissidentified,
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
  o2::framework::ConfigurableAxis sign{"sign", {{3, -1.5, 1.5}}, "Sign"};                          \
  o2::framework::ConfigurableAxis pdgCodes{"pdgCodes", {{8001, -4000.5, 4000.5}}, "MC ONLY: PDG codes of reconstructed V0s"};

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
  o2::framework::Configurable<bool> plotOrigins{"plotOrigins", true, "MC ONLY: Plot pt vs cosPa for different particle origins"};
  o2::framework::Configurable<std::vector<int>> pdgCodesForMothersOfSecondary{"pdgCodesForMothersOfSecondary", {3312, 3334}, "MC ONLY: PDG codes of mothers of secondaries (Max 3 will be considered)"};
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
   {kTruePtVsPt, o2::framework::kTH2F, "hTruePtVsPt", "True transverse momentum vs transverse momentum; p_{T,True} (GeV/#it{c}); p_{T,True} (GeV/#it{c})"},
   {kTrueEtaVsEta, o2::framework::kTH2F, "hTrueEtaVsEta", "True pseudorapdity vs pseudorapdity; #eta_{True}; #eta"},
   {kTruePhiVsPhi, o2::framework::kTH2F, "hTruePhiVsPhi", "True azimuthal angle vs azimuthal angle; #varphi_{True}; #varphi"},
   {kNoMcParticle, o2::framework::kTH2F, "hNoMcParticle", "Wrongly reconstructed particles; p_{T} (GeV/#it{c}); cos(#alpha)"},
   {kPrimary, o2::framework::kTH2F, "hPrimary", "Primary particles; p_{T} (GeV/#it{c}); cos(#alpha)"},
   {kFromWrongCollision, o2::framework::kTH2F, "hFromWrongCollision", "Particles associated to wrong collision; p_{T} (GeV/#it{c}); cos(#alpha)"},
   {kFromMaterial, o2::framework::kTH2F, "hFromMaterial", "Particles from material; p_{T} (GeV/#it{c}); cos(#alpha)"},
   {kMissidentified, o2::framework::kTH2F, "hMissidentified", "Missidentified particles (fake/wrong PDG code); p_{T} (GeV/#it{c}); cos(#alpha)"},
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

#define V0_HIST_MC_MAP(conf)               \
  {kTruePtVsPt, {conf.pt, conf.pt}},       \
    {kTrueEtaVsEta, {conf.eta, conf.eta}}, \
    {kTruePhiVsPhi, {conf.phi, conf.phi}}, \
    {kPdg, {conf.pdgCodes}},               \
    {kPdgMother, {conf.pdgCodes}},         \
    {kPdgPartonicMother, {conf.pdgCodes}},

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

#define V0_HIST_MC_QA_MAP(confAnalysis, confQa)             \
  {kNoMcParticle, {confAnalysis.pt, confQa.cosPa}},         \
    {kPrimary, {confAnalysis.pt, confQa.cosPa}},            \
    {kFromWrongCollision, {confAnalysis.pt, confQa.cosPa}}, \
    {kFromMaterial, {confAnalysis.pt, confQa.cosPa}},       \
    {kMissidentified, {confAnalysis.pt, confQa.cosPa}},     \
    {kSecondary1, {confAnalysis.pt, confQa.cosPa}},         \
    {kSecondary2, {confAnalysis.pt, confQa.cosPa}},         \
    {kSecondary3, {confAnalysis.pt, confQa.cosPa}},         \
    {kSecondaryOther, {confAnalysis.pt, confQa.cosPa}},

template <typename T>
auto makeV0HistSpecMap(const T& confBinningAnalysis)
{
  return std::map<V0Hist, std::vector<framework::AxisSpec>>{
    V0_HIST_ANALYSIS_MAP(confBinningAnalysis)};
}

template <typename T>
auto makeV0McHistSpecMap(const T& confBinningAnalysis)
{
  return std::map<V0Hist, std::vector<framework::AxisSpec>>{
    V0_HIST_ANALYSIS_MAP(confBinningAnalysis)
      V0_HIST_MC_MAP(confBinningAnalysis)};
}

template <typename T1, typename T2>
std::map<V0Hist, std::vector<framework::AxisSpec>> makeV0QaHistSpecMap(T1 const& confBinningAnalysis, T2 const& confBinningQa)
{
  return std::map<V0Hist, std::vector<framework::AxisSpec>>{
    V0_HIST_ANALYSIS_MAP(confBinningAnalysis)
      V0_HIST_QA_MAP(confBinningAnalysis, confBinningQa)};
}

template <typename T1, typename T2>
std::map<V0Hist, std::vector<framework::AxisSpec>> makeV0McQaHistSpecMap(T1 const& confBinningAnalysis, T2 const& confBinningQa)
{
  return std::map<V0Hist, std::vector<framework::AxisSpec>>{
    V0_HIST_ANALYSIS_MAP(confBinningAnalysis)
      V0_HIST_QA_MAP(confBinningAnalysis, confBinningQa)
        V0_HIST_MC_MAP(confBinningAnalysis)
          V0_HIST_MC_QA_MAP(confBinningAnalysis, confBinningQa)};
}

#undef V0_HIST_ANALYSIS_MAP
#undef V0_HIST_QA_MAP
#undef V0_HIST_MC_MAP
#undef V0_HIST_MC_QA_MAP

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

  // init for analysis and mc
  template <modes::Mode mode, typename T>
  void init(o2::framework::HistogramRegistry* registry,
            std::map<V0Hist, std::vector<o2::framework::AxisSpec>> const& V0Specs,
            T const& ConfV0Selection,
            std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> const& PosDauSpecs,
            std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> const& NegDauSpecs)
  {
    mHistogramRegistry = registry;
    mPdgCode = std::abs(ConfV0Selection.pdgCodeAbs.value);

    int posDauPdgCodeAbs = 0;
    int negDauPdgCodeAbs = 0;
    const int absCharge = 1;
    const int signPlus = 1;
    const int signMinus = -1;

    if (mPdgCode == PDG_t::kK0Short) {
      posDauPdgCodeAbs = std::abs(PDG_t::kPiPlus);
      negDauPdgCodeAbs = std::abs(PDG_t::kPiMinus);
    } else if (mPdgCode == PDG_t::kLambda0) {
      if (ConfV0Selection.sign.value > 0) {
        posDauPdgCodeAbs = std::abs(PDG_t::kProton);
        negDauPdgCodeAbs = std::abs(PDG_t::kPiMinus);
      } else {
        mPdgCode = -1 * mPdgCode; // switch sign for antilambda
        posDauPdgCodeAbs = std::abs(PDG_t::kPiPlus);
        negDauPdgCodeAbs = std::abs(PDG_t::kProtonBar);
      }
    } else {
      LOG(fatal) << "PDG code for V0 has to be either Lambda or K0short";
    }

    mPosDauManager.template init<mode>(registry, PosDauSpecs, absCharge, signPlus, posDauPdgCodeAbs);
    mNegDauManager.template init<mode>(registry, NegDauSpecs, absCharge, signMinus, negDauPdgCodeAbs);

    if constexpr (modes::isFlagSet(mode, modes::Mode::kAnalysis)) {
      this->initAnalysis(V0Specs);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kMc)) {
      this->initMc(V0Specs);
    }
  }

  // init for analysis and qa and and mc
  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4>
  void init(o2::framework::HistogramRegistry* registry,
            std::map<V0Hist, std::vector<o2::framework::AxisSpec>> const& V0Specs,
            T1 const& ConfV0Selection,
            T2 const& ConfV0BinningQa,
            std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> const& PosDauSpecs,
            T3 const& ConfPosDauBinningQa,
            std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> const& NegDauSpecs,
            T4 const& ConfNegDauBinningQa)
  {
    mHistogramRegistry = registry;
    mPdgCode = std::abs(ConfV0Selection.pdgCodeAbs.value);
    this->enableOptionalHistograms(ConfV0BinningQa);

    int posDauPdgCode = 0;
    int negDauPdgCode = 0;
    const int absCharge = 1;
    const int signPlus = 1;
    const int signMinus = -1;

    if (mPdgCode == PDG_t::kK0Short) {
      posDauPdgCode = std::abs(PDG_t::kPiPlus);
      negDauPdgCode = std::abs(PDG_t::kPiMinus);
    } else if (mPdgCode == PDG_t::kLambda0) {
      if (ConfV0Selection.sign.value > 0) {
        posDauPdgCode = std::abs(PDG_t::kProton);
        negDauPdgCode = std::abs(PDG_t::kPiMinus);
      } else {
        mPdgCode = -1 * mPdgCode; // set pdg code for antilambda
        posDauPdgCode = std::abs(PDG_t::kPiPlus);
        negDauPdgCode = std::abs(PDG_t::kProtonBar);
      }
    } else {
      LOG(fatal) << "PDG code for V0 has to be either Lambda or K0short";
    }

    mPosDauManager.template init<mode>(registry, PosDauSpecs, absCharge, signPlus, posDauPdgCode, ConfPosDauBinningQa);
    mNegDauManager.template init<mode>(registry, NegDauSpecs, absCharge, signMinus, negDauPdgCode, ConfNegDauBinningQa);

    if constexpr (modes::isFlagSet(mode, modes::Mode::kAnalysis)) {
      this->initAnalysis(V0Specs);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kQa)) {
      this->initQa(V0Specs);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kMc)) {
      this->initMc(V0Specs);
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
      this->fillAnalysis(v0candidate);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kQa)) {
      this->fillQa(v0candidate);
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
      this->fillAnalysis(v0candidate);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kQa)) {
      this->fillQa(v0candidate);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kMc)) {
      this->template fillMc<mode>(v0candidate, mcParticles, mcMothers, mcPartonicMothers);
    }
  }

 private:
  template <typename T1>
  void enableOptionalHistograms(T1 const& V0ConfBinningQa)
  {
    mPlot2d = V0ConfBinningQa.plot2d.value;
    mPlotOrigins = V0ConfBinningQa.plotOrigins.value;
    mPlotNSecondaries = V0ConfBinningQa.pdgCodesForMothersOfSecondary.value.size();
    for (std::size_t i = 0; i < MaxSecondary; i++) {
      if (i < V0ConfBinningQa.pdgCodesForMothersOfSecondary.value.size()) {
        mPdgCodesSecondaryMother.at(i) = std::abs(V0ConfBinningQa.pdgCodesForMothersOfSecondary.value.at(i));
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
    mHistogramRegistry->add(mcDir + getHistNameV2(kTruePtVsPt, HistTable), getHistDesc(kTruePtVsPt, HistTable), getHistType(kTruePtVsPt, HistTable), {V0Specs.at(kTruePtVsPt)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kTrueEtaVsEta, HistTable), getHistDesc(kTrueEtaVsEta, HistTable), getHistType(kTrueEtaVsEta, HistTable), {V0Specs.at(kTrueEtaVsEta)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kTruePhiVsPhi, HistTable), getHistDesc(kTruePhiVsPhi, HistTable), getHistType(kTruePhiVsPhi, HistTable), {V0Specs.at(kTruePhiVsPhi)});

    // mc origin can be configured here
    const framework::AxisSpec axisOrigin = {static_cast<int>(modes::McOrigin::kMcOriginLast), -0.5, static_cast<double>(modes::McOrigin::kMcOriginLast) - 0.5};
    mHistogramRegistry->add(mcDir + getHistNameV2(kOrigin, HistTable), getHistDesc(kOrigin, HistTable), getHistType(kOrigin, HistTable), {axisOrigin});
    mHistogramRegistry->get<TH1>(HIST(v0Prefix) + HIST(McDir) + HIST(histmanager::getHistName(kOrigin, HistTable)))->GetXaxis()->SetBinLabel(1 + static_cast<int>(modes::McOrigin::kNoMcParticle), modes::mcOriginToString(modes::McOrigin::kNoMcParticle));
    mHistogramRegistry->get<TH1>(HIST(v0Prefix) + HIST(McDir) + HIST(histmanager::getHistName(kOrigin, HistTable)))->GetXaxis()->SetBinLabel(1 + static_cast<int>(modes::McOrigin::kFromWrongCollision), modes::mcOriginToString(modes::McOrigin::kFromWrongCollision));
    mHistogramRegistry->get<TH1>(HIST(v0Prefix) + HIST(McDir) + HIST(histmanager::getHistName(kOrigin, HistTable)))->GetXaxis()->SetBinLabel(1 + static_cast<int>(modes::McOrigin::kPhysicalPrimary), modes::mcOriginToString(modes::McOrigin::kPhysicalPrimary));
    mHistogramRegistry->get<TH1>(HIST(v0Prefix) + HIST(McDir) + HIST(histmanager::getHistName(kOrigin, HistTable)))->GetXaxis()->SetBinLabel(1 + static_cast<int>(modes::McOrigin::kFromSecondaryDecay), modes::mcOriginToString(modes::McOrigin::kFromSecondaryDecay));
    mHistogramRegistry->get<TH1>(HIST(v0Prefix) + HIST(McDir) + HIST(histmanager::getHistName(kOrigin, HistTable)))->GetXaxis()->SetBinLabel(1 + static_cast<int>(modes::McOrigin::kFromMaterial), modes::mcOriginToString(modes::McOrigin::kFromMaterial));
    mHistogramRegistry->get<TH1>(HIST(v0Prefix) + HIST(McDir) + HIST(histmanager::getHistName(kOrigin, HistTable)))->GetXaxis()->SetBinLabel(1 + static_cast<int>(modes::McOrigin::kMissidentified), modes::mcOriginToString(modes::McOrigin::kMissidentified));

    mHistogramRegistry->add(mcDir + getHistNameV2(kPdg, HistTable), getHistDesc(kPdg, HistTable), getHistType(kPdg, HistTable), {V0Specs.at(kPdg)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kPdgMother, HistTable), getHistDesc(kPdgMother, HistTable), getHistType(kPdgMother, HistTable), {V0Specs.at(kPdgMother)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kPdgPartonicMother, HistTable), getHistDesc(kPdgPartonicMother, HistTable), getHistType(kPdgPartonicMother, HistTable), {V0Specs.at(kPdgPartonicMother)});

    if (mPlotOrigins) {
      mHistogramRegistry->add(mcDir + getHistNameV2(kNoMcParticle, HistTable), getHistDesc(kNoMcParticle, HistTable), getHistType(kNoMcParticle, HistTable), {V0Specs.at(kNoMcParticle)});
      mHistogramRegistry->add(mcDir + getHistNameV2(kPrimary, HistTable), getHistDesc(kPrimary, HistTable), getHistType(kPrimary, HistTable), {V0Specs.at(kPrimary)});
      mHistogramRegistry->add(mcDir + getHistNameV2(kFromWrongCollision, HistTable), getHistDesc(kFromWrongCollision, HistTable), getHistType(kFromWrongCollision, HistTable), {V0Specs.at(kFromWrongCollision)});
      mHistogramRegistry->add(mcDir + getHistNameV2(kFromMaterial, HistTable), getHistDesc(kFromMaterial, HistTable), getHistType(kFromMaterial, HistTable), {V0Specs.at(kFromMaterial)});
      mHistogramRegistry->add(mcDir + getHistNameV2(kMissidentified, HistTable), getHistDesc(kMissidentified, HistTable), getHistType(kMissidentified, HistTable), {V0Specs.at(kMissidentified)});

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

  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4>
  void fillMc(T1 const& v0Candidate, T2 const& /*mcParticles*/, T3 const& /*mcMothers*/, T4 const& /*mcPartonicMothers*/)
  {
    // No MC Particle
    if (!v0Candidate.has_fMcParticle()) {
      mHistogramRegistry->fill(HIST(v0Prefix) + HIST(McDir) + HIST(getHistName(kPdg, HistTable)), 0);
      mHistogramRegistry->fill(HIST(v0Prefix) + HIST(McDir) + HIST(getHistName(kOrigin, HistTable)), static_cast<float>(modes::McOrigin::kNoMcParticle));
      if constexpr (modes::isFlagSet(mode, modes::Mode::kQa)) {
        if (mPlotOrigins) {
          mHistogramRegistry->fill(HIST(v0Prefix) + HIST(McDir) + HIST(getHistName(kNoMcParticle, HistTable)), v0Candidate.pt(), v0Candidate.cosPa());
        }
      }
      return;
    }

    // Retrieve MC particle
    auto mcParticle = v0Candidate.template fMcParticle_as<T2>();

    // missidentifed particles are special case
    // whether a particle is missidentfied or not cannot be known by the producer so we check it here
    bool isMissidentified = mcParticle.pdgCode() != mPdgCode;

    mHistogramRegistry->fill(HIST(v0Prefix) + HIST(McDir) + HIST(getHistName(kTruePtVsPt, HistTable)), mcParticle.pt(), v0Candidate.pt());
    mHistogramRegistry->fill(HIST(v0Prefix) + HIST(McDir) + HIST(getHistName(kTrueEtaVsEta, HistTable)), mcParticle.eta(), v0Candidate.eta());
    mHistogramRegistry->fill(HIST(v0Prefix) + HIST(McDir) + HIST(getHistName(kTruePhiVsPhi, HistTable)), mcParticle.phi(), v0Candidate.phi());
    if (isMissidentified) {
      mHistogramRegistry->fill(HIST(v0Prefix) + HIST(McDir) + HIST(getHistName(kOrigin, HistTable)), static_cast<int>(modes::McOrigin::kMissidentified));
    } else {
      mHistogramRegistry->fill(HIST(v0Prefix) + HIST(McDir) + HIST(getHistName(kOrigin, HistTable)), mcParticle.origin());
    }
    mHistogramRegistry->fill(HIST(v0Prefix) + HIST(McDir) + HIST(getHistName(kPdg, HistTable)), mcParticle.pdgCode());

    // get mother
    if (v0Candidate.has_fMcMother()) {
      auto mother = v0Candidate.template fMcMother_as<T3>();
      mHistogramRegistry->fill(HIST(v0Prefix) + HIST(McDir) + HIST(getHistName(kPdgMother, HistTable)), mother.pdgCode());
    } else {
      mHistogramRegistry->fill(HIST(v0Prefix) + HIST(McDir) + HIST(getHistName(kPdgMother, HistTable)), 0);
    }

    // get partonic mother
    if (v0Candidate.has_fMcPartMoth()) {
      auto partonicMother = v0Candidate.template fMcPartMoth_as<T4>();
      mHistogramRegistry->fill(HIST(v0Prefix) + HIST(McDir) + HIST(getHistName(kPdgPartonicMother, HistTable)), partonicMother.pdgCode());
    } else {
      mHistogramRegistry->fill(HIST(v0Prefix) + HIST(McDir) + HIST(getHistName(kPdgPartonicMother, HistTable)), 0);
    }

    if constexpr (modes::isFlagSet(mode, modes::Mode::kQa)) {
      if (mPlotOrigins) {
        // check first if particle is missidentified
        if (isMissidentified) {
          // if it is, we fill it as such
          mHistogramRegistry->fill(HIST(v0Prefix) + HIST(McDir) + HIST(getHistName(kMissidentified, HistTable)), v0Candidate.pt(), v0Candidate.cosPa());
        } else {
          // if not, we fill it acccoridng to its origin
          switch (static_cast<modes::McOrigin>(mcParticle.origin())) {
            case modes::McOrigin::kPhysicalPrimary:
              mHistogramRegistry->fill(HIST(v0Prefix) + HIST(McDir) + HIST(getHistName(kPrimary, HistTable)), v0Candidate.pt(), v0Candidate.cosPa());
              break;
            case modes::McOrigin::kFromWrongCollision:
              mHistogramRegistry->fill(HIST(v0Prefix) + HIST(McDir) + HIST(getHistName(kFromWrongCollision, HistTable)), v0Candidate.pt(), v0Candidate.cosPa());
              break;
            case modes::McOrigin::kFromMaterial:
              mHistogramRegistry->fill(HIST(v0Prefix) + HIST(McDir) + HIST(getHistName(kFromMaterial, HistTable)), v0Candidate.pt(), v0Candidate.cosPa());
              break;
            case modes::McOrigin::kFromSecondaryDecay:
              if (v0Candidate.has_fMcMother()) {
                auto mother = v0Candidate.template fMcMother_as<T3>();
                int motherPdgCode = std::abs(mother.pdgCode());
                // Switch on PDG of the mother
                if (mPlotNSecondaries >= histmanager::kSecondaryPlotLevel1 && motherPdgCode == mPdgCodesSecondaryMother[0]) {
                  mHistogramRegistry->fill(HIST(v0Prefix) + HIST(McDir) + HIST(getHistName(kSecondary1, HistTable)), v0Candidate.pt(), v0Candidate.cosPa());
                } else if (mPlotNSecondaries >= histmanager::kSecondaryPlotLevel2 && motherPdgCode == mPdgCodesSecondaryMother[1]) {
                  mHistogramRegistry->fill(HIST(v0Prefix) + HIST(McDir) + HIST(getHistName(kSecondary2, HistTable)), v0Candidate.pt(), v0Candidate.cosPa());
                } else if (mPlotNSecondaries >= histmanager::kSecondaryPlotLevel3 && motherPdgCode == mPdgCodesSecondaryMother[2]) {
                  mHistogramRegistry->fill(HIST(v0Prefix) + HIST(McDir) + HIST(getHistName(kSecondary3, HistTable)), v0Candidate.pt(), v0Candidate.cosPa());
                } else {
                  mHistogramRegistry->fill(HIST(v0Prefix) + HIST(McDir) + HIST(getHistName(kSecondaryOther, HistTable)), v0Candidate.pt(), v0Candidate.cosPa());
                }
              }
              break;
            default:
              LOG(warn) << "Encounted partilce with unknown origin!";
              break;
          }
        }
      }
    }
  }

  o2::framework::HistogramRegistry* mHistogramRegistry = nullptr;
  int mPdgCode = 0;
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
