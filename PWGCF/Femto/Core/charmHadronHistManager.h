// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file charmHadronHistManager.h
/// \brief histogram manager for charm hadron histograms
/// \author Igor Ptak, WUT, igor.tomasz.ptak@cern.ch

#ifndef PWGCF_FEMTO_CORE_CHARMHADRONHISTMANAGER_H_
#define PWGCF_FEMTO_CORE_CHARMHADRONHISTMANAGER_H_

#include "PWGCF/Femto/Core/histManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/trackHistManager.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/Logger.h>

#include <TH1.h>
#include <TPDGCode.h>

#include <array>
#include <map>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace o2::analysis::femto::charmhadronhistmanager
{
enum CharmHadronHist {
  kPt,
  kEta,
  kPhi,
  kMass,
  kSign,
  kPtVsMass,
  kPtVsEta,
  kPtVsPhi,
  kPhiVsEta,
  kMassD0,
  kMassD0bar,
  kMlBkg,
  kMlPrompt,
  kMlNonPrompt,
  kCpa,
  kCpaXY,
  kDecayLength,
  kDecayLengthXY,
  kImpactParameterProduct,
  kCosThetaStar,

  // mc
  kOrigin,
  kPdg,
  kPdgMother,
  kTruePtVsPt,
  kTrueEtaVsEta,
  kTruePhiVsPhi,
  kPtVsOrigin,

  kCharmHadronHistLast
};

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define CHARMHADRON_DEFAULT_BINNING(defaultMassMin, defaultMassMax)                                                                            \
  o2::framework::ConfigurableAxis pt{"pt", {{600, 0, 6}}, "Pt"};                                                                               \
  o2::framework::ConfigurableAxis eta{"eta", {{300, -1.5, 1.5}}, "Eta"};                                                                       \
  o2::framework::ConfigurableAxis phi{"phi", {{720, 0, 1.f * o2::constants::math::TwoPI}}, "Phi"};                                             \
  o2::framework::ConfigurableAxis mass{"mass", {{200, (defaultMassMin), (defaultMassMax)}}, "Mass"};                                           \
  o2::framework::ConfigurableAxis sign{"sign", {{3, -1.5, 1.5}}, "Sign"};                                                                      \
  o2::framework::ConfigurableAxis charmHadrons{"charmHadrons", {{8001, -4000.5, 4000.5}}, "MC ONLY: CharmHadrons codes of reconstructed D0s"}; \
  o2::framework::ConfigurableAxis pt2d{"pt2d", {{240, 0, 6}}, "Pt for 2D QA"};                                                                 \
  o2::framework::ConfigurableAxis eta2d{"eta2d", {{200, -1.5, 1.5}}, "Eta for 2D QA"};                                                         \
  o2::framework::ConfigurableAxis phi2d{"phi2d", {{200, 0, 1.f * o2::constants::math::TwoPI}}, "Phi for 2D QA"};

template <auto& Prefix>
struct ConfD0Binning : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  CHARMHADRON_DEFAULT_BINNING(1.7, 2.0)
};

#undef CHARMHADRON_DEFAULT_BINNING

constexpr const char PrefixD0Binning1[] = "D0Binning1";
constexpr const char PrefixD0Binning2[] = "D0Binning2";
using ConfD0Binning1 = ConfD0Binning<PrefixD0Binning1>;
using ConfD0Binning2 = ConfD0Binning<PrefixD0Binning2>;

template <auto& Prefix>
struct ConfD0QaBinning : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  o2::framework::Configurable<bool> plotTopology{"plotTopology", true, "Generate topological QA plots (cpa, decayLength, impactParameterProduct, cosThetaStar)"};
  o2::framework::ConfigurableAxis massD0{"massD0", {{200, 1.7, 2.0}}, "Mass of candidates selected as D0"};
  o2::framework::ConfigurableAxis massD0bar{"massD0bar", {{200, 1.7, 2.0}}, "Mass of candidates selected as D0bar"};
  o2::framework::ConfigurableAxis mlScore{"mlScore", {{100, 0.f, 1.f}}, "BDT ML score (bkg/prompt/non-prompt)"};
  o2::framework::ConfigurableAxis cpa{"cpa", {{100, 0.9f, 1.f}}, "Cosine of pointing angle"};
  o2::framework::ConfigurableAxis decayLength{"decayLength", {{200, 0.f, 0.2f}}, "Decay length (cm)"};
  o2::framework::ConfigurableAxis impactParameterProduct{"impactParameterProduct", {{200, -0.001f, 0.001f}}, "Product of daughter impact parameters (cm^2)"};
  o2::framework::ConfigurableAxis cosThetaStar{"cosThetaStar", {{100, -1.f, 1.f}}, "Cosine of decay angle in D0 rest frame"};
};

constexpr const char PrefixD0QaBinning1[] = "D0QaBinning1";
using ConfD0QaBinning1 = ConfD0QaBinning<PrefixD0QaBinning1>;

// must be in sync with enum CharmHadronHist
// the enum gives the correct index in the array
constexpr std::array<histmanager::HistInfo<CharmHadronHist>, kCharmHadronHistLast> HistTable = {
  {{kPt, o2::framework::HistType::kTH1F, "hPt", "Transverse Momentum; p_{T} (GeV/#it{c}); Entries"},
   {kEta, o2::framework::HistType::kTH1F, "hEta", "Pseudorapdity; #eta; Entries"},
   {kPhi, o2::framework::HistType::kTH1F, "hPhi", "Azimuthal angle; #varphi; Entries"},
   {kMass, o2::framework::HistType::kTH1F, "hMass", "Invariant Mass; m_{Inv} (GeV/#it{c}^{2}); Entries"},
   {kSign, o2::framework::HistType::kTH1F, "hSign", "Sign (-1 -> D0bar, +1 -> D0); sign; Entries"},
   {kPtVsMass, o2::framework::HistType::kTH2F, "hPtVsMass", "Transverse momentum vs invariant mass; p_{T} (GeV/#it{c}); m_{Inv} (GeV/#it{c}^{2})"},
   {kPtVsEta, o2::framework::HistType::kTH2F, "hPtVsEta", "p_{T} vs #eta; p_{T} (GeV/#it{c}); #eta"},
   {kPtVsPhi, o2::framework::HistType::kTH2F, "hPtVsPhi", "p_{T} vs #varphi; p_{T} (GeV/#it{c}); #varphi"},
   {kPhiVsEta, o2::framework::HistType::kTH2F, "hPhiVsEta", "#varphi vs #eta; #varphi; #eta"},
   {kMassD0, o2::framework::HistType::kTH1F, "hMassD0", "Invariant mass of D^{0} candidates; m_{Inv} (GeV/#it{c}^{2}); Entries"},
   {kMassD0bar, o2::framework::HistType::kTH1F, "hMassD0bar", "Invariant mass of #bar{D}^{0} candidates; m_{Inv} (GeV/#it{c}^{2}); Entries"},
   {kMlBkg, o2::framework::HistType::kTH1F, "hMlBkg", "BDT background score; ML score (bkg); Entries"},
   {kMlPrompt, o2::framework::HistType::kTH1F, "hMlPrompt", "BDT prompt score; ML score (prompt); Entries"},
   {kMlNonPrompt, o2::framework::HistType::kTH1F, "hMlNonPrompt", "BDT non-prompt score; ML score (non-prompt); Entries"},
   {kCpa, o2::framework::HistType::kTH1F, "hCpa", "Cosine of pointing angle; cos(#alpha); Entries"},
   {kCpaXY, o2::framework::HistType::kTH1F, "hCpaXY", "Cosine of pointing angle (XY); cos(#alpha)_{XY}; Entries"},
   {kDecayLength, o2::framework::HistType::kTH1F, "hDecayLength", "Decay length; L (cm); Entries"},
   {kDecayLengthXY, o2::framework::HistType::kTH1F, "hDecayLengthXY", "Decay length (XY); L_{XY} (cm); Entries"},
   {kImpactParameterProduct, o2::framework::HistType::kTH1F, "hImpactParameterProduct", "Product of daughter impact parameters; d_{0}^{K} #times d_{0}^{#pi} (cm^{2}); Entries"},
   {kCosThetaStar, o2::framework::HistType::kTH1F, "hCosThetaStar", "Cosine of decay angle in D0 rest frame; cos(#theta*); Entries"},
   {kOrigin, o2::framework::HistType::kTH1F, "hOrigin", "MC origin (prompt / non-prompt); origin; Entries"},
   {kPdg, o2::framework::HistType::kTH1F, "hPdg", "PDG code of matched generated particle; PDG code; Entries"},
   {kPdgMother, o2::framework::HistType::kTH1F, "hPdgMother", "PDG code of the mother; PDG code; Entries"},
   {kTruePtVsPt, o2::framework::HistType::kTH2F, "hTruePtVsPt", "True vs reconstructed p_{T}; p_{T,true} (GeV/#it{c}); p_{T} (GeV/#it{c})"},
   {kTrueEtaVsEta, o2::framework::HistType::kTH2F, "hTrueEtaVsEta", "True vs reconstructed #eta; #eta_{true}; #eta"},
   {kTruePhiVsPhi, o2::framework::HistType::kTH2F, "hTruePhiVsPhi", "True vs reconstructed #varphi; #varphi_{true}; #varphi"},
   {kPtVsOrigin, o2::framework::HistType::kTH2F, "hPtVsOrigin", "p_{T} vs MC origin; p_{T} (GeV/#it{c}); origin"}},
};

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define CHARMHADRON_HIST_ANALYSIS_MAP(conf)  \
  {kPt, {(conf).pt}},                        \
    {kEta, {(conf).eta}},                    \
    {kPhi, {(conf).phi}},                    \
    {kMass, {(conf).mass}},                  \
    {kSign, {(conf).sign}},                  \
    {kPtVsMass, {(conf).pt, (conf).mass}},   \
    {kPtVsEta, {(conf).pt2d, (conf).eta2d}}, \
    {kPtVsPhi, {(conf).pt2d, (conf).phi2d}}, \
    {kPhiVsEta, {(conf).phi2d, (conf).eta2d}},

template <typename T>
auto makeD0HistSpecMap(const T& confBinningAnalysis)
{
  return std::map<CharmHadronHist, std::vector<o2::framework::AxisSpec>>{
    CHARMHADRON_HIST_ANALYSIS_MAP(confBinningAnalysis)};
}

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define CHARMHADRON_HIST_MC_MAP(conf)          \
  {kPdg, {(conf).charmHadrons}},               \
    {kPdgMother, {(conf).charmHadrons}},       \
    {kTruePtVsPt, {(conf).pt, (conf).pt}},     \
    {kTrueEtaVsEta, {(conf).eta, (conf).eta}}, \
    {kTruePhiVsPhi, {(conf).phi, (conf).phi}}, \
    {kPtVsOrigin, {(conf).pt2d}},

template <typename T>
auto makeD0McHistSpecMap(const T& confBinningAnalysis)
{
  return std::map<CharmHadronHist, std::vector<o2::framework::AxisSpec>>{
    CHARMHADRON_HIST_ANALYSIS_MAP(confBinningAnalysis)
      CHARMHADRON_HIST_MC_MAP(confBinningAnalysis)};
}

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define CHARMHADRON_HIST_QA_MAP(conf)                           \
  {kMassD0, {(conf).massD0}},                                   \
    {kMassD0bar, {(conf).massD0bar}},                           \
    {kMlBkg, {(conf).mlScore}},                                 \
    {kMlPrompt, {(conf).mlScore}},                              \
    {kMlNonPrompt, {(conf).mlScore}},                           \
    {kCpa, {(conf).cpa}},                                       \
    {kCpaXY, {(conf).cpa}},                                     \
    {kDecayLength, {(conf).decayLength}},                       \
    {kDecayLengthXY, {(conf).decayLength}},                     \
    {kImpactParameterProduct, {(conf).impactParameterProduct}}, \
    {kCosThetaStar, {(conf).cosThetaStar}},

template <typename T>
auto makeD0QaHistSpecMap(const T& confBinningQa)
{
  return std::map<CharmHadronHist, std::vector<o2::framework::AxisSpec>>{
    CHARMHADRON_HIST_QA_MAP(confBinningQa)};
}

template <typename T1, typename T2>
auto makeD0McQaHistSpecMap(const T1& confBinningAnalysis, const T2& confBinningQa)
{
  return std::map<CharmHadronHist, std::vector<o2::framework::AxisSpec>>{
    CHARMHADRON_HIST_ANALYSIS_MAP(confBinningAnalysis)
      CHARMHADRON_HIST_QA_MAP(confBinningQa)
        CHARMHADRON_HIST_MC_MAP(confBinningAnalysis)};
}

#undef CHARMHADRON_HIST_ANALYSIS_MAP
#undef CHARMHADRON_HIST_MC_MAP
#undef CHARMHADRON_HIST_QA_MAP

// prefixes for the output directories in the histogram registry
constexpr char PrefixD01[] = "D01/";
constexpr char PrefixD02[] = "D02/";
constexpr char PrefixD0Qa[] = "D0QA/";
constexpr std::string_view AnalysisDir = "Analysis/";
constexpr std::string_view McDir = "MC/";
constexpr std::string_view QaDir = "QA/";

/// \class CharmHadronHistManager
/// \brief Class for histogramming charm hadron properties
template <auto& charmHadronPrefix,
          auto& prong0Prefix,
          auto& prong1Prefix,
          modes::CharmHadron charmHadron>
class CharmHadronHistManager
{
 public:
  CharmHadronHistManager() = default;
  ~CharmHadronHistManager() = default;

  // init for analysis
  template <modes::Mode mode, typename T>
  void init(o2::framework::HistogramRegistry* registry,
            std::map<CharmHadronHist, std::vector<o2::framework::AxisSpec>> const& CharmHadronSpecs,
            T const& ConfCharmHadronSelection,
            std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> const& Prong0Specs,
            std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> const& Prong1Specs)
  {
    mHistogramRegistry = registry;
    mPdgCode = std::abs(ConfCharmHadronSelection.pdgCodeAbs.value);

    auto [prong0PdgCodeAbs, prong1PdgCodeAbs] = this->resolveProngPdgCodes(ConfCharmHadronSelection.sign.value);

    mProng0Manager.template init<mode>(registry, Prong0Specs, AbsCharge, SignPlus, prong0PdgCodeAbs);
    mProng1Manager.template init<mode>(registry, Prong1Specs, AbsCharge, SignMinus, prong1PdgCodeAbs);

    if constexpr (modes::isFlagSet(mode, modes::Mode::kReco)) {
      this->initAnalysis(CharmHadronSpecs);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kMc)) {
      this->initMc(CharmHadronSpecs);
    }
  }

  // init for analysis and qa
  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4>
  void init(o2::framework::HistogramRegistry* registry,
            std::map<CharmHadronHist, std::vector<o2::framework::AxisSpec>> const& CharmHadronSpecs,
            std::map<CharmHadronHist, std::vector<o2::framework::AxisSpec>> const& CharmHadronQaSpecs,
            T1 const& ConfCharmHadronSelection,
            T2 const& ConfCharmHadronQaBinning,
            std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> const& Prong0Specs,
            T3 const& ConfProng0BinningQa,
            std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> const& Prong1Specs,
            T4 const& ConfProng1BinningQa)
  {
    mHistogramRegistry = registry;
    mPdgCode = std::abs(ConfCharmHadronSelection.pdgCodeAbs.value);
    this->enableOptionalHistograms(ConfCharmHadronQaBinning);

    auto [prong0PdgCodeAbs, prong1PdgCodeAbs] = this->resolveProngPdgCodes(ConfCharmHadronSelection.sign.value);

    mProng0Manager.template init<mode>(registry, Prong0Specs, AbsCharge, SignPlus, prong0PdgCodeAbs, ConfProng0BinningQa);
    mProng1Manager.template init<mode>(registry, Prong1Specs, AbsCharge, SignMinus, prong1PdgCodeAbs, ConfProng1BinningQa);

    if constexpr (modes::isFlagSet(mode, modes::Mode::kReco)) {
      this->initAnalysis(CharmHadronSpecs);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kQa)) {
      this->initQa(CharmHadronQaSpecs);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kMc)) {
      this->initMc(CharmHadronSpecs);
    }
  }

  template <modes::Mode mode, typename T1, typename T2>
  void fill(T1 const& charmHadronCandidate, T2 const& tracks)
  {
    auto prong0 = tracks.rawIteratorAt(charmHadronCandidate.posDauId() - tracks.offset());
    mProng0Manager.template fill<mode>(prong0, tracks);
    auto prong1 = tracks.rawIteratorAt(charmHadronCandidate.negDauId() - tracks.offset());
    mProng1Manager.template fill<mode>(prong1, tracks);

    if constexpr (modes::isFlagSet(mode, modes::Mode::kReco)) {
      this->fillAnalysis(charmHadronCandidate);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kQa)) {
      this->fillQa(charmHadronCandidate);
    }
  }

  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  void fill(T1 const& charmHadronCandidate, T2 const& tracks, T3 const& col, T4 const& mcParticles, T5 const& mcMothers, T6 const& mcPartonicMothers)
  {
    auto prong0 = tracks.rawIteratorAt(charmHadronCandidate.posDauId() - tracks.offset());
    mProng0Manager.template fill<mode>(prong0, tracks, col, mcParticles, mcMothers, mcPartonicMothers);
    auto prong1 = tracks.rawIteratorAt(charmHadronCandidate.negDauId() - tracks.offset());
    mProng1Manager.template fill<mode>(prong1, tracks, col, mcParticles, mcMothers, mcPartonicMothers);

    if constexpr (modes::isFlagSet(mode, modes::Mode::kReco)) {
      this->fillAnalysis(charmHadronCandidate);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kQa)) {
      this->fillQa(charmHadronCandidate);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kMc)) {
      this->fillMc(charmHadronCandidate, mcParticles, mcMothers, mcPartonicMothers);
    }
  }

 private:
  static constexpr int AbsCharge = 1;
  static constexpr int SignPlus = 1;
  static constexpr int SignMinus = -1;

  // charge of the prongs is fixed by the PWGHF reconstruction: prong0 is the positive
  // daughter, prong1 the negative one. The D0/D0bar hypothesis only swaps which prong is
  // the pion and which is the kaon, not their charge.
  std::pair<int, int> resolveProngPdgCodes(int sign)
  {
    if (mPdgCode != o2::constants::physics::Pdg::kD0) {
      LOG(fatal) << "PDG code for charm hadron has to be D0";
    }
    if (sign > 0) {
      // D0 -> pi+ K-
      return {std::abs(PDG_t::kPiPlus), std::abs(PDG_t::kKMinus)};
    }
    // D0bar -> K+ pi-
    mPdgCode = -1 * mPdgCode; // switch sign for D0bar
    return {std::abs(PDG_t::kKPlus), std::abs(PDG_t::kPiMinus)};
  }

  template <typename T>
  void enableOptionalHistograms(T const& ConfCharmHadronQaBinning)
  {
    mPlotTopology = ConfCharmHadronQaBinning.plotTopology.value;
  }

  void initAnalysis(std::map<CharmHadronHist, std::vector<o2::framework::AxisSpec>> const& CharmHadronSpecs)
  {
    std::string analysisDir = std::string(charmHadronPrefix) + std::string(AnalysisDir);
    mHistogramRegistry->add(analysisDir + getHistNameV2(kPt, HistTable), getHistDesc(kPt, HistTable), getHistType(kPt, HistTable), {CharmHadronSpecs.at(kPt)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kEta, HistTable), getHistDesc(kEta, HistTable), getHistType(kEta, HistTable), {CharmHadronSpecs.at(kEta)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kPhi, HistTable), getHistDesc(kPhi, HistTable), getHistType(kPhi, HistTable), {CharmHadronSpecs.at(kPhi)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kMass, HistTable), getHistDesc(kMass, HistTable), getHistType(kMass, HistTable), {CharmHadronSpecs.at(kMass)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kSign, HistTable), getHistDesc(kSign, HistTable), getHistType(kSign, HistTable), {CharmHadronSpecs.at(kSign)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kPtVsMass, HistTable), getHistDesc(kPtVsMass, HistTable), getHistType(kPtVsMass, HistTable), {CharmHadronSpecs.at(kPtVsMass)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kPtVsEta, HistTable), getHistDesc(kPtVsEta, HistTable), getHistType(kPtVsEta, HistTable), {CharmHadronSpecs.at(kPtVsEta)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kPtVsPhi, HistTable), getHistDesc(kPtVsPhi, HistTable), getHistType(kPtVsPhi, HistTable), {CharmHadronSpecs.at(kPtVsPhi)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kPhiVsEta, HistTable), getHistDesc(kPhiVsEta, HistTable), getHistType(kPhiVsEta, HistTable), {CharmHadronSpecs.at(kPhiVsEta)});
  }

  template <typename T>
  void fillAnalysis(T const& charmHadronCandidate)
  {
    mHistogramRegistry->fill(HIST(charmHadronPrefix) + HIST(AnalysisDir) + HIST(getHistName(kPt, HistTable)), charmHadronCandidate.pt());
    mHistogramRegistry->fill(HIST(charmHadronPrefix) + HIST(AnalysisDir) + HIST(getHistName(kEta, HistTable)), charmHadronCandidate.eta());
    mHistogramRegistry->fill(HIST(charmHadronPrefix) + HIST(AnalysisDir) + HIST(getHistName(kPhi, HistTable)), charmHadronCandidate.phi());
    mHistogramRegistry->fill(HIST(charmHadronPrefix) + HIST(AnalysisDir) + HIST(getHistName(kMass, HistTable)), charmHadronCandidate.mass());
    mHistogramRegistry->fill(HIST(charmHadronPrefix) + HIST(AnalysisDir) + HIST(getHistName(kSign, HistTable)), charmHadronCandidate.sign());
    mHistogramRegistry->fill(HIST(charmHadronPrefix) + HIST(AnalysisDir) + HIST(getHistName(kPtVsMass, HistTable)), charmHadronCandidate.pt(), charmHadronCandidate.mass());
    mHistogramRegistry->fill(HIST(charmHadronPrefix) + HIST(AnalysisDir) + HIST(getHistName(kPtVsEta, HistTable)), charmHadronCandidate.pt(), charmHadronCandidate.eta());
    mHistogramRegistry->fill(HIST(charmHadronPrefix) + HIST(AnalysisDir) + HIST(getHistName(kPtVsPhi, HistTable)), charmHadronCandidate.pt(), charmHadronCandidate.phi());
    mHistogramRegistry->fill(HIST(charmHadronPrefix) + HIST(AnalysisDir) + HIST(getHistName(kPhiVsEta, HistTable)), charmHadronCandidate.phi(), charmHadronCandidate.eta());
  }

  void initQa(std::map<CharmHadronHist, std::vector<o2::framework::AxisSpec>> const& CharmHadronQaSpecs)
  {
    std::string qaDir = std::string(charmHadronPrefix) + std::string(QaDir);
    mHistogramRegistry->add(qaDir + getHistNameV2(kMassD0, HistTable), getHistDesc(kMassD0, HistTable), getHistType(kMassD0, HistTable), {CharmHadronQaSpecs.at(kMassD0)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kMassD0bar, HistTable), getHistDesc(kMassD0bar, HistTable), getHistType(kMassD0bar, HistTable), {CharmHadronQaSpecs.at(kMassD0bar)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kMlBkg, HistTable), getHistDesc(kMlBkg, HistTable), getHistType(kMlBkg, HistTable), {CharmHadronQaSpecs.at(kMlBkg)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kMlPrompt, HistTable), getHistDesc(kMlPrompt, HistTable), getHistType(kMlPrompt, HistTable), {CharmHadronQaSpecs.at(kMlPrompt)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kMlNonPrompt, HistTable), getHistDesc(kMlNonPrompt, HistTable), getHistType(kMlNonPrompt, HistTable), {CharmHadronQaSpecs.at(kMlNonPrompt)});

    if (mPlotTopology) {
      mHistogramRegistry->add(qaDir + getHistNameV2(kCpa, HistTable), getHistDesc(kCpa, HistTable), getHistType(kCpa, HistTable), {CharmHadronQaSpecs.at(kCpa)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kCpaXY, HistTable), getHistDesc(kCpaXY, HistTable), getHistType(kCpaXY, HistTable), {CharmHadronQaSpecs.at(kCpaXY)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kDecayLength, HistTable), getHistDesc(kDecayLength, HistTable), getHistType(kDecayLength, HistTable), {CharmHadronQaSpecs.at(kDecayLength)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kDecayLengthXY, HistTable), getHistDesc(kDecayLengthXY, HistTable), getHistType(kDecayLengthXY, HistTable), {CharmHadronQaSpecs.at(kDecayLengthXY)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kImpactParameterProduct, HistTable), getHistDesc(kImpactParameterProduct, HistTable), getHistType(kImpactParameterProduct, HistTable), {CharmHadronQaSpecs.at(kImpactParameterProduct)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kCosThetaStar, HistTable), getHistDesc(kCosThetaStar, HistTable), getHistType(kCosThetaStar, HistTable), {CharmHadronQaSpecs.at(kCosThetaStar)});
    }
  }

  void initMc(std::map<CharmHadronHist, std::vector<o2::framework::AxisSpec>> const& CharmHadronSpecs)
  {
    std::string mcDir = std::string(charmHadronPrefix) + std::string(McDir);

    const o2::framework::AxisSpec axisOrigin = {static_cast<int>(modes::McOrigin::kMcOriginLast), -0.5, static_cast<double>(modes::McOrigin::kMcOriginLast) - 0.5};
    mHistogramRegistry->add(mcDir + getHistNameV2(kOrigin, HistTable), getHistDesc(kOrigin, HistTable), getHistType(kOrigin, HistTable), {axisOrigin});
    mHistogramRegistry->get<TH1>(HIST(charmHadronPrefix) + HIST(McDir) + HIST(histmanager::getHistName(kOrigin, HistTable)))->GetXaxis()->SetBinLabel(1 + static_cast<int>(modes::McOrigin::kPrompt), modes::mcOriginToString(modes::McOrigin::kPrompt));
    mHistogramRegistry->get<TH1>(HIST(charmHadronPrefix) + HIST(McDir) + HIST(histmanager::getHistName(kOrigin, HistTable)))->GetXaxis()->SetBinLabel(1 + static_cast<int>(modes::McOrigin::kNonPrompt), modes::mcOriginToString(modes::McOrigin::kNonPrompt));

    mHistogramRegistry->add(mcDir + getHistNameV2(kPdg, HistTable), getHistDesc(kPdg, HistTable), getHistType(kPdg, HistTable), {CharmHadronSpecs.at(kPdg)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kPdgMother, HistTable), getHistDesc(kPdgMother, HistTable), getHistType(kPdgMother, HistTable), {CharmHadronSpecs.at(kPdgMother)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kTruePtVsPt, HistTable), getHistDesc(kTruePtVsPt, HistTable), getHistType(kTruePtVsPt, HistTable), {CharmHadronSpecs.at(kTruePtVsPt)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kTrueEtaVsEta, HistTable), getHistDesc(kTrueEtaVsEta, HistTable), getHistType(kTrueEtaVsEta, HistTable), {CharmHadronSpecs.at(kTrueEtaVsEta)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kTruePhiVsPhi, HistTable), getHistDesc(kTruePhiVsPhi, HistTable), getHistType(kTruePhiVsPhi, HistTable), {CharmHadronSpecs.at(kTruePhiVsPhi)});
    std::vector<o2::framework::AxisSpec> ptVsOriginAxes = CharmHadronSpecs.at(kPtVsOrigin);
    ptVsOriginAxes.push_back(axisOrigin);
    mHistogramRegistry->add(mcDir + getHistNameV2(kPtVsOrigin, HistTable), getHistDesc(kPtVsOrigin, HistTable), getHistType(kPtVsOrigin, HistTable), ptVsOriginAxes);
  }

  template <typename T>
  void fillQa(T const& charmHadronCandidate)
  {
    // invariant mass split by hypothesis: D0 (sign > 0) or D0bar (sign < 0)
    if (charmHadronCandidate.sign() > 0) {
      mHistogramRegistry->fill(HIST(charmHadronPrefix) + HIST(QaDir) + HIST(getHistName(kMassD0, HistTable)), charmHadronCandidate.mass());
    } else {
      mHistogramRegistry->fill(HIST(charmHadronPrefix) + HIST(QaDir) + HIST(getHistName(kMassD0bar, HistTable)), charmHadronCandidate.mass());
    }

    // BDT scores of the accepted hypothesis: D0 (sign > 0) or D0bar (sign < 0)
    float mlBkg = 0.f;
    float mlPrompt = 0.f;
    float mlNonPrompt = 0.f;
    if (charmHadronCandidate.sign() > 0) {
      mlBkg = charmHadronCandidate.mlProbD0Bkg();
      mlPrompt = charmHadronCandidate.mlProbD0Prompt();
      mlNonPrompt = charmHadronCandidate.mlProbD0NonPrompt();
    } else {
      mlBkg = charmHadronCandidate.mlProbD0barBkg();
      mlPrompt = charmHadronCandidate.mlProbD0barPrompt();
      mlNonPrompt = charmHadronCandidate.mlProbD0barNonPrompt();
    }
    mHistogramRegistry->fill(HIST(charmHadronPrefix) + HIST(QaDir) + HIST(getHistName(kMlBkg, HistTable)), mlBkg);
    mHistogramRegistry->fill(HIST(charmHadronPrefix) + HIST(QaDir) + HIST(getHistName(kMlPrompt, HistTable)), mlPrompt);
    mHistogramRegistry->fill(HIST(charmHadronPrefix) + HIST(QaDir) + HIST(getHistName(kMlNonPrompt, HistTable)), mlNonPrompt);

    // topological variables PWGHF selects on
    if (mPlotTopology) {
      mHistogramRegistry->fill(HIST(charmHadronPrefix) + HIST(QaDir) + HIST(getHistName(kCpa, HistTable)), charmHadronCandidate.cpa());
      mHistogramRegistry->fill(HIST(charmHadronPrefix) + HIST(QaDir) + HIST(getHistName(kCpaXY, HistTable)), charmHadronCandidate.cpaXY());
      mHistogramRegistry->fill(HIST(charmHadronPrefix) + HIST(QaDir) + HIST(getHistName(kDecayLength, HistTable)), charmHadronCandidate.decayLength());
      mHistogramRegistry->fill(HIST(charmHadronPrefix) + HIST(QaDir) + HIST(getHistName(kDecayLengthXY, HistTable)), charmHadronCandidate.decayLengthXY());
      mHistogramRegistry->fill(HIST(charmHadronPrefix) + HIST(QaDir) + HIST(getHistName(kImpactParameterProduct, HistTable)), charmHadronCandidate.impactParameterProduct());
      mHistogramRegistry->fill(HIST(charmHadronPrefix) + HIST(QaDir) + HIST(getHistName(kCosThetaStar, HistTable)), charmHadronCandidate.cosThetaStar());
    }
  }

  template <typename T1, typename T2, typename T3, typename T4>
  void fillMc(T1 const& charmHadronCandidate, T2 const& /*mcParticles*/, T3 const& /*mcMothers*/, T4 const& /*mcPartonicMothers*/)
  {
    // no matched generated particle: reconstructed, but not a true D0
    if (!charmHadronCandidate.has_fMcParticle()) {
      mHistogramRegistry->fill(HIST(charmHadronPrefix) + HIST(McDir) + HIST(getHistName(kPdg, HistTable)), 0);
      mHistogramRegistry->fill(HIST(charmHadronPrefix) + HIST(McDir) + HIST(getHistName(kOrigin, HistTable)), static_cast<float>(modes::McOrigin::kNoMcParticle));
      return;
    }

    auto mcParticle = charmHadronCandidate.template fMcParticle_as<T2>();

    // resolution: generated vs reconstructed kinematics
    mHistogramRegistry->fill(HIST(charmHadronPrefix) + HIST(McDir) + HIST(getHistName(kTruePtVsPt, HistTable)), mcParticle.pt(), charmHadronCandidate.pt());
    mHistogramRegistry->fill(HIST(charmHadronPrefix) + HIST(McDir) + HIST(getHistName(kTrueEtaVsEta, HistTable)), mcParticle.eta(), charmHadronCandidate.eta());
    mHistogramRegistry->fill(HIST(charmHadronPrefix) + HIST(McDir) + HIST(getHistName(kTruePhiVsPhi, HistTable)), mcParticle.phi(), charmHadronCandidate.phi());

    // origin is resolved to prompt / non-prompt in the mc builder
    mHistogramRegistry->fill(HIST(charmHadronPrefix) + HIST(McDir) + HIST(getHistName(kOrigin, HistTable)), mcParticle.origin());
    mHistogramRegistry->fill(HIST(charmHadronPrefix) + HIST(McDir) + HIST(getHistName(kPtVsOrigin, HistTable)), charmHadronCandidate.pt(), mcParticle.origin());
    mHistogramRegistry->fill(HIST(charmHadronPrefix) + HIST(McDir) + HIST(getHistName(kPdg, HistTable)), mcParticle.pdgCode());

    // pdg of the mother, i.e. the beauty hadron for non-prompt D0s
    if (mcParticle.has_fMcMother()) {
      auto mother = mcParticle.template fMcMother_as<T3>();
      mHistogramRegistry->fill(HIST(charmHadronPrefix) + HIST(McDir) + HIST(getHistName(kPdgMother, HistTable)), mother.pdgCode());
    } else {
      mHistogramRegistry->fill(HIST(charmHadronPrefix) + HIST(McDir) + HIST(getHistName(kPdgMother, HistTable)), 0);
    }
  }

  o2::framework::HistogramRegistry* mHistogramRegistry = nullptr;
  int mPdgCode = 0;
  bool mPlotTopology = true;

  trackhistmanager::TrackHistManager<prong0Prefix> mProng0Manager;
  trackhistmanager::TrackHistManager<prong1Prefix> mProng1Manager;
};
}; // namespace o2::analysis::femto::charmhadronhistmanager
#endif // PWGCF_FEMTO_CORE_CHARMHADRONHISTMANAGER_H_
