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

/// \file sigma0HistManager.h
/// \brief histogram manager for Sigma0 histograms
/// \author Yakiv Paroviak, CTU in Prague, yakiv.paroviak@cern.ch

#ifndef PWGCF_FEMTO_CORE_SIGMA0HISTMANAGER_H_
#define PWGCF_FEMTO_CORE_SIGMA0HISTMANAGER_H_

#include "PWGCF/Femto/Core/histManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/trackHistManager.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/Logger.h>

#include <TH1.h>
#include <TPDGCode.h>

#include <array>
#include <cstdlib>
#include <map>
#include <string>
#include <string_view>
#include <vector>

namespace o2::analysis::femto::sigma0histmanager
{
// enum for Sigma0 histograms
enum Sigma0Hist {
  // analysis
  kPt,
  kEta,
  kPhi,
  kMass,
  kSign,
  kPtVsMass,
  // qa variables (from FSigma0Extras)
  kLambdaMass,
  kLambdaPt,
  kLambdaCosPa,
  kLambdaDauDca,
  kLambdaTransRadius,
  kPhotonPt,
  kPhotonCosPa,
  kPhotonTransRadius,
  kPhotonChi2Ndf,
  kOpeningAngle,
  // 2d qa
  kPtVsEta,
  kPtVsPhi,
  kPhiVsEta,
  kPtVsLambdaMass,
  kPtVsOpeningAngle,
  kLambdaPtVsPhotonPt,
  // mc
  kOrigin,
  kPdg,
  kPdgMother,
  kPdgPartonicMother,
  kTruePtVsPt,
  kTrueEtaVsEta,
  kTruePhiVsPhi,
  kSigma0HistLast
};

template <auto& Prefix>
struct ConfSigma0Binning : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  o2::framework::ConfigurableAxis pt{"pt", {{600, 0, 6}}, "Pt"};
  o2::framework::ConfigurableAxis eta{"eta", {{300, -1.5, 1.5}}, "Eta"};
  o2::framework::ConfigurableAxis phi{"phi", {{720, 0, 1.f * o2::constants::math::TwoPI}}, "Phi"};
  o2::framework::ConfigurableAxis mass{"mass", {{200, 1.1, 1.3}}, "Mass"};
  o2::framework::ConfigurableAxis sign{"sign", {{3, -1.5, 1.5}}, "Sign"};
  o2::framework::ConfigurableAxis pdgCodes{"pdgCodes", {{8001, -4000.5, 4000.5}}, "MC ONLY: PDG codes of reconstructed Sigma0s"};
};

constexpr const char PrefixSigma0Binning[] = "Sigma0Binning";
using ConfSigma0BinningDefault = ConfSigma0Binning<PrefixSigma0Binning>;

template <auto& Prefix>
struct ConfSigma0QaBinning : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  o2::framework::Configurable<bool> plot2d{"plot2d", true, "Enable 2d QA histograms"};
  o2::framework::ConfigurableAxis lambdaMass{"lambdaMass", {{200, 1., 1.2}}, "Mass of daughter Lambda"};
  o2::framework::ConfigurableAxis lambdaPt{"lambdaPt", {{600, 0, 6}}, "Pt of daughter Lambda"};
  o2::framework::ConfigurableAxis lambdaCosPa{"lambdaCosPa", {{100, 0.9, 1}}, "Cosine of pointing angle of daughter Lambda"};
  o2::framework::ConfigurableAxis lambdaDauDca{"lambdaDauDca", {{150, 0, 1.5}}, "DCA of Lambda daughters at Lambda decay vertex"};
  o2::framework::ConfigurableAxis lambdaTransRadius{"lambdaTransRadius", {{100, 0, 100}}, "Transverse radius of daughter Lambda decay vertex"};
  o2::framework::ConfigurableAxis photonPt{"photonPt", {{600, 0, 6}}, "Pt of daughter photon"};
  o2::framework::ConfigurableAxis photonCosPa{"photonCosPa", {{100, 0.9, 1}}, "Cosine of pointing angle of daughter photon"};
  o2::framework::ConfigurableAxis photonTransRadius{"photonTransRadius", {{180, 0, 90}}, "Transverse radius of photon conversion point"};
  o2::framework::ConfigurableAxis photonChi2Ndf{"photonChi2Ndf", {{100, 0, 50}}, "Chi2/NDF of photon conversion vertex"};
  o2::framework::ConfigurableAxis openingAngle{"openingAngle", {{180, 0, 1.f * o2::constants::math::PI}}, "Opening angle between Lambda and photon"};
};

constexpr const char PrefixSigma0QaBinning[] = "Sigma0QaBinning";
using ConfSigma0QaBinningDefault = ConfSigma0QaBinning<PrefixSigma0QaBinning>;

// must be in sync with enum Sigma0Hist
// the enum gives the correct index in the array
constexpr std::array<histmanager::HistInfo<Sigma0Hist>, kSigma0HistLast> HistTable = {
  {{kPt, o2::framework::HistType::kTH1F, "hPt", "Transverse Momentum; p_{T} (GeV/#it{c}); Entries"},
   {kEta, o2::framework::HistType::kTH1F, "hEta", "Pseudorapdity; #eta; Entries"},
   {kPhi, o2::framework::HistType::kTH1F, "hPhi", "Azimuthal angle; #varphi; Entries"},
   {kMass, o2::framework::HistType::kTH1F, "hMass", "Invariant Mass; m_{#Lambda#gamma} (GeV/#it{c}^{2}); Entries"},
   {kSign, o2::framework::HistType::kTH1F, "hSign", "Sign (-1 -> antiparticle, +1 -> particle); sign; Entries"},
   {kPtVsMass, o2::framework::HistType::kTH2F, "hPtVsMass", "Transverse momentum vs invariant mass; p_{T} (GeV/#it{c}); m_{#Lambda#gamma} (GeV/#it{c}^{2})"},
   {kLambdaMass, o2::framework::HistType::kTH1F, "hLambdaMass", "Invariant mass of daughter #Lambda; m_{#Lambda dau} (GeV/#it{c}^{2}); Entries"},
   {kLambdaPt, o2::framework::HistType::kTH1F, "hLambdaPt", "Transverse momentum of daughter #Lambda; p_{T,#Lambda dau} (GeV/#it{c}); Entries"},
   {kLambdaCosPa, o2::framework::HistType::kTH1F, "hLambdaCosPa", "Cosine of pointing angle of daughter #Lambda; cos_{#Lambda dau}(#alpha); Entries"},
   {kLambdaDauDca, o2::framework::HistType::kTH1F, "hLambdaDauDca", "Daughter DCA at #Lambda decay vertex; DCA_{#Lambda dau decay vertex} (cm); Entries"},
   {kLambdaTransRadius, o2::framework::HistType::kTH1F, "hLambdaTransRadius", "Transverse radius of daughter #Lambda; r_{xy,#Lambda dau} (cm); Entries"},
   {kPhotonPt, o2::framework::HistType::kTH1F, "hPhotonPt", "Transverse momentum of daughter #gamma; p_{T,#gamma dau} (GeV/#it{c}); Entries"},
   {kPhotonCosPa, o2::framework::HistType::kTH1F, "hPhotonCosPa", "Cosine of pointing angle of daughter #gamma; cos_{#gamma dau}(#alpha); Entries"},
   {kPhotonTransRadius, o2::framework::HistType::kTH1F, "hPhotonTransRadius", "Transverse radius of #gamma conversion point; r_{xy,#gamma dau} (cm); Entries"},
   {kPhotonChi2Ndf, o2::framework::HistType::kTH1F, "hPhotonChi2Ndf", "#chi^{2}/NDF of #gamma conversion vertex; #chi^{2}/NDF; Entries"},
   {kOpeningAngle, o2::framework::HistType::kTH1F, "hOpeningAngle", "Opening angle between #Lambda and #gamma; #theta_{#Lambda#gamma}; Entries"},
   {kPtVsEta, o2::framework::HistType::kTH2F, "hPtVsEta", "p_{T} vs #eta; p_{T} (GeV/#it{c}) ; #eta"},
   {kPtVsPhi, o2::framework::HistType::kTH2F, "hPtVsPhi", "p_{T} vs #varphi; p_{T} (GeV/#it{c}) ; #varphi"},
   {kPhiVsEta, o2::framework::HistType::kTH2F, "hPhiVsEta", "#varphi vs #eta; #varphi ; #eta"},
   {kPtVsLambdaMass, o2::framework::HistType::kTH2F, "hPtVsLambdaMass", "p_{T} vs mass of daughter #Lambda; p_{T} (GeV/#it{c}); m_{#Lambda dau} (GeV/#it{c}^{2})"},
   {kPtVsOpeningAngle, o2::framework::HistType::kTH2F, "hPtVsOpeningAngle", "p_{T} vs opening angle; p_{T} (GeV/#it{c}); #theta_{#Lambda#gamma}"},
   {kLambdaPtVsPhotonPt, o2::framework::HistType::kTH2F, "hLambdaPtVsPhotonPt", "p_{T} of daughter #Lambda vs p_{T} of daughter #gamma; p_{T,#Lambda dau} (GeV/#it{c}); p_{T,#gamma dau} (GeV/#it{c})"},
   {kOrigin, o2::framework::HistType::kTH1F, "hOrigin", "Status Codes (=Origin); Status Code; Entries"},
   {kPdg, o2::framework::HistType::kTH1F, "hPdg", "PDG Codes of reconstructed #Sigma^{0}; PDG Code; Entries"},
   {kPdgMother, o2::framework::HistType::kTH1F, "hPdgMother", "PDG Codes of mother of reconstructed #Sigma^{0}; PDG Code; Entries"},
   {kPdgPartonicMother, o2::framework::HistType::kTH1F, "hPdgPartonicMother", "PDG Codes of partonic mother of reconstructed #Sigma^{0}; PDG Code; Entries"},
   {kTruePtVsPt, o2::framework::HistType::kTH2F, "hTruePtVsPt", "True transverse momentum vs transverse momentum; p_{T,True} (GeV/#it{c}); p_{T} (GeV/#it{c})"},
   {kTrueEtaVsEta, o2::framework::HistType::kTH2F, "hTrueEtaVsEta", "True pseudorapdity vs pseudorapdity; #eta_{True}; #eta"},
   {kTruePhiVsPhi, o2::framework::HistType::kTH2F, "hTruePhiVsPhi", "True azimuthal angle vs azimuthal angle; #varphi_{True}; #varphi"}},
};

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define SIGMA0_HIST_ANALYSIS_MAP(conf) \
  {kPt, {(conf).pt}},                  \
    {kEta, {(conf).eta}},              \
    {kPhi, {(conf).phi}},              \
    {kMass, {(conf).mass}},            \
    {kSign, {(conf).sign}},            \
    {kPtVsMass, {(conf).pt, (conf).mass}},

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define SIGMA0_HIST_MC_MAP(conf)               \
  {kTruePtVsPt, {(conf).pt, (conf).pt}},       \
    {kTrueEtaVsEta, {(conf).eta, (conf).eta}}, \
    {kTruePhiVsPhi, {(conf).phi, (conf).phi}}, \
    {kPdg, {(conf).pdgCodes}},                 \
    {kPdgMother, {(conf).pdgCodes}},           \
    {kPdgPartonicMother, {(conf).pdgCodes}},

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define SIGMA0_HIST_QA_MAP(confAnalysis, confQa)                     \
  {kLambdaMass, {(confQa).lambdaMass}},                              \
    {kLambdaPt, {(confQa).lambdaPt}},                                \
    {kLambdaCosPa, {(confQa).lambdaCosPa}},                          \
    {kLambdaDauDca, {(confQa).lambdaDauDca}},                        \
    {kLambdaTransRadius, {(confQa).lambdaTransRadius}},              \
    {kPhotonPt, {(confQa).photonPt}},                                \
    {kPhotonCosPa, {(confQa).photonCosPa}},                          \
    {kPhotonTransRadius, {(confQa).photonTransRadius}},              \
    {kPhotonChi2Ndf, {(confQa).photonChi2Ndf}},                      \
    {kOpeningAngle, {(confQa).openingAngle}},                        \
    {kPtVsEta, {(confAnalysis).pt, (confAnalysis).eta}},             \
    {kPtVsPhi, {(confAnalysis).pt, (confAnalysis).phi}},             \
    {kPhiVsEta, {(confAnalysis).phi, (confAnalysis).eta}},           \
    {kPtVsLambdaMass, {(confAnalysis).pt, (confQa).lambdaMass}},     \
    {kPtVsOpeningAngle, {(confAnalysis).pt, (confQa).openingAngle}}, \
    {kLambdaPtVsPhotonPt, {(confQa).lambdaPt, (confQa).photonPt}},

template <typename T>
auto makeSigma0HistSpecMap(const T& confBinningAnalysis)
{
  return std::map<Sigma0Hist, std::vector<o2::framework::AxisSpec>>{
    SIGMA0_HIST_ANALYSIS_MAP(confBinningAnalysis)};
}

template <typename T>
auto makeSigma0McHistSpecMap(const T& confBinningAnalysis)
{
  return std::map<Sigma0Hist, std::vector<o2::framework::AxisSpec>>{
    SIGMA0_HIST_ANALYSIS_MAP(confBinningAnalysis)
      SIGMA0_HIST_MC_MAP(confBinningAnalysis)};
}

template <typename T1, typename T2>
std::map<Sigma0Hist, std::vector<o2::framework::AxisSpec>> makeSigma0QaHistSpecMap(T1 const& confBinningAnalysis, T2 const& confBinningQa)
{
  return std::map<Sigma0Hist, std::vector<o2::framework::AxisSpec>>{
    SIGMA0_HIST_ANALYSIS_MAP(confBinningAnalysis)
      SIGMA0_HIST_QA_MAP(confBinningAnalysis, confBinningQa)};
}

template <typename T1, typename T2>
std::map<Sigma0Hist, std::vector<o2::framework::AxisSpec>> makeSigma0McQaHistSpecMap(T1 const& confBinningAnalysis, T2 const& confBinningQa)
{
  return std::map<Sigma0Hist, std::vector<o2::framework::AxisSpec>>{
    SIGMA0_HIST_ANALYSIS_MAP(confBinningAnalysis)
      SIGMA0_HIST_QA_MAP(confBinningAnalysis, confBinningQa)
        SIGMA0_HIST_MC_MAP(confBinningAnalysis)};
}

#undef SIGMA0_HIST_ANALYSIS_MAP
#undef SIGMA0_HIST_MC_MAP
#undef SIGMA0_HIST_QA_MAP

constexpr char PrefixSigma0Qa[] = "Sigma0QA/";
constexpr char PrefixSigma0Se[] = "Sigma0/SE/";
constexpr char PrefixSigma0Me[] = "Sigma0/ME/";

constexpr std::string_view AnalysisDir = "Analysis/";
constexpr std::string_view QaDir = "QA/";
constexpr std::string_view McDir = "MC/";

/// \class Sigma0HistManager
/// \brief Class for Sigma0 and its dau track histogramming
///
/// The Sigma0 tables store only the indices of the Lambda daughter tracks (and of the photon), so the
/// daughter tracks are histogrammed with track histogram managers. QA of the Lambda and photon
/// daughters is read from FSigma0Extras, which needs to be joined to the Sigma0 table for QA modes.
template <auto& sigma0Prefix,
          auto& posDauPrefix,
          auto& negDauPrefix>
class Sigma0HistManager
{
 public:
  Sigma0HistManager() = default;
  ~Sigma0HistManager() = default;

  template <modes::Mode mode, typename T>
  void init(o2::framework::HistogramRegistry* registry,
            std::map<Sigma0Hist, std::vector<o2::framework::AxisSpec>> const& sigma0Specs,
            T const& ConfSigma0Selection,
            std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> const& PosDauSpecs,
            std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> const& NegDauSpecs)
  {
    mHistogramRegistry = registry;
    this->setPdgCodes(ConfSigma0Selection);

    mPosDauManager.template init<mode>(registry, PosDauSpecs, AbsCharge, SignPlus, mPosDauPdgCodeAbs);
    mNegDauManager.template init<mode>(registry, NegDauSpecs, AbsCharge, SignMinus, mNegDauPdgCodeAbs);

    if constexpr (modes::isFlagSet(mode, modes::Mode::kReco)) {
      this->initAnalysis(sigma0Specs);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kMc)) {
      this->initMc(sigma0Specs);
    }
  }

  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4>
  void init(o2::framework::HistogramRegistry* registry,
            std::map<Sigma0Hist, std::vector<o2::framework::AxisSpec>> const& sigma0Specs,
            T1 const& ConfSigma0Selection,
            T2 const& ConfSigma0BinningQa,
            std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> const& PosDauSpecs,
            T3 const& ConfPosDauQaBinning,
            std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> const& NegDauSpecs,
            T4 const& ConfNegDauQaBinning)
  {
    mHistogramRegistry = registry;
    this->setPdgCodes(ConfSigma0Selection);
    mPlot2d = ConfSigma0BinningQa.plot2d.value;

    mPosDauManager.template init<mode>(registry, PosDauSpecs, AbsCharge, SignPlus, mPosDauPdgCodeAbs, ConfPosDauQaBinning);
    mNegDauManager.template init<mode>(registry, NegDauSpecs, AbsCharge, SignMinus, mNegDauPdgCodeAbs, ConfNegDauQaBinning);

    if constexpr (modes::isFlagSet(mode, modes::Mode::kReco)) {
      this->initAnalysis(sigma0Specs);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kQa)) {
      this->initQa(sigma0Specs);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kMc)) {
      this->initMc(sigma0Specs);
    }
  }

  template <modes::Mode mode, typename T1, typename T2>
  void fill(T1 const& sigma0Candidate, T2 const& tracks)
  {
    auto posDaughter = tracks.rawIteratorAt(sigma0Candidate.posDauId() - tracks.offset());
    mPosDauManager.template fill<mode>(posDaughter, tracks);
    auto negDaughter = tracks.rawIteratorAt(sigma0Candidate.negDauId() - tracks.offset());
    mNegDauManager.template fill<mode>(negDaughter, tracks);

    if constexpr (modes::isFlagSet(mode, modes::Mode::kReco)) {
      this->fillAnalysis(sigma0Candidate);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kQa)) {
      this->fillQa(sigma0Candidate);
    }
  }

  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  void fill(T1 const& sigma0Candidate, T2 const& tracks, T3 const& col, T4 const& mcParticles, T5 const& mcMothers, T6 const& mcPartonicMothers)
  {
    auto posDaughter = tracks.rawIteratorAt(sigma0Candidate.posDauId() - tracks.offset());
    mPosDauManager.template fill<mode>(posDaughter, tracks, col, mcParticles, mcMothers, mcPartonicMothers);
    auto negDaughter = tracks.rawIteratorAt(sigma0Candidate.negDauId() - tracks.offset());
    mNegDauManager.template fill<mode>(negDaughter, tracks, col, mcParticles, mcMothers, mcPartonicMothers);

    if constexpr (modes::isFlagSet(mode, modes::Mode::kReco)) {
      this->fillAnalysis(sigma0Candidate);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kQa)) {
      this->fillQa(sigma0Candidate);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kMc)) {
      this->template fillMc<T4, T5, T6>(sigma0Candidate, col);
    }
  }

 private:
  static constexpr int AbsCharge = 1;
  static constexpr int SignPlus = 1;
  static constexpr int SignMinus = -1;

  // Sigma0 -> Lambda (p pi-) gamma; AntiSigma0 -> AntiLambda (pbar pi+) gamma
  template <typename T>
  void setPdgCodes(T const& ConfSigma0Selection)
  {
    mPdgCode = std::abs(ConfSigma0Selection.pdgCodeAbs.value);
    if (mPdgCode != PDG_t::kSigma0) {
      LOG(fatal) << "PDG code for Sigma0 has to be " << PDG_t::kSigma0;
    }
    if (ConfSigma0Selection.sign.value > 0) {
      mPosDauPdgCodeAbs = std::abs(PDG_t::kProton);
      mNegDauPdgCodeAbs = std::abs(PDG_t::kPiMinus);
    } else {
      mPdgCode = -1 * mPdgCode; // switch sign for AntiSigma0
      mPosDauPdgCodeAbs = std::abs(PDG_t::kPiPlus);
      mNegDauPdgCodeAbs = std::abs(PDG_t::kProtonBar);
    }
  }

  void initAnalysis(std::map<Sigma0Hist, std::vector<o2::framework::AxisSpec>> const& sigma0Specs)
  {
    std::string analysisDir = std::string(sigma0Prefix) + std::string(AnalysisDir);
    mHistogramRegistry->add(analysisDir + getHistNameV2(kPt, HistTable), getHistDesc(kPt, HistTable), getHistType(kPt, HistTable), {sigma0Specs.at(kPt)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kEta, HistTable), getHistDesc(kEta, HistTable), getHistType(kEta, HistTable), {sigma0Specs.at(kEta)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kPhi, HistTable), getHistDesc(kPhi, HistTable), getHistType(kPhi, HistTable), {sigma0Specs.at(kPhi)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kMass, HistTable), getHistDesc(kMass, HistTable), getHistType(kMass, HistTable), {sigma0Specs.at(kMass)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kSign, HistTable), getHistDesc(kSign, HistTable), getHistType(kSign, HistTable), {sigma0Specs.at(kSign)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kPtVsMass, HistTable), getHistDesc(kPtVsMass, HistTable), getHistType(kPtVsMass, HistTable), {sigma0Specs.at(kPtVsMass)});
  }

  void initQa(std::map<Sigma0Hist, std::vector<o2::framework::AxisSpec>> const& sigma0Specs)
  {
    std::string qaDir = std::string(sigma0Prefix) + std::string(QaDir);
    mHistogramRegistry->add(qaDir + getHistNameV2(kLambdaMass, HistTable), getHistDesc(kLambdaMass, HistTable), getHistType(kLambdaMass, HistTable), {sigma0Specs.at(kLambdaMass)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kLambdaPt, HistTable), getHistDesc(kLambdaPt, HistTable), getHistType(kLambdaPt, HistTable), {sigma0Specs.at(kLambdaPt)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kLambdaCosPa, HistTable), getHistDesc(kLambdaCosPa, HistTable), getHistType(kLambdaCosPa, HistTable), {sigma0Specs.at(kLambdaCosPa)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kLambdaDauDca, HistTable), getHistDesc(kLambdaDauDca, HistTable), getHistType(kLambdaDauDca, HistTable), {sigma0Specs.at(kLambdaDauDca)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kLambdaTransRadius, HistTable), getHistDesc(kLambdaTransRadius, HistTable), getHistType(kLambdaTransRadius, HistTable), {sigma0Specs.at(kLambdaTransRadius)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kPhotonPt, HistTable), getHistDesc(kPhotonPt, HistTable), getHistType(kPhotonPt, HistTable), {sigma0Specs.at(kPhotonPt)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kPhotonCosPa, HistTable), getHistDesc(kPhotonCosPa, HistTable), getHistType(kPhotonCosPa, HistTable), {sigma0Specs.at(kPhotonCosPa)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kPhotonTransRadius, HistTable), getHistDesc(kPhotonTransRadius, HistTable), getHistType(kPhotonTransRadius, HistTable), {sigma0Specs.at(kPhotonTransRadius)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kPhotonChi2Ndf, HistTable), getHistDesc(kPhotonChi2Ndf, HistTable), getHistType(kPhotonChi2Ndf, HistTable), {sigma0Specs.at(kPhotonChi2Ndf)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kOpeningAngle, HistTable), getHistDesc(kOpeningAngle, HistTable), getHistType(kOpeningAngle, HistTable), {sigma0Specs.at(kOpeningAngle)});

    if (mPlot2d) {
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsEta, HistTable), getHistDesc(kPtVsEta, HistTable), getHistType(kPtVsEta, HistTable), {sigma0Specs.at(kPtVsEta)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsPhi, HistTable), getHistDesc(kPtVsPhi, HistTable), getHistType(kPtVsPhi, HistTable), {sigma0Specs.at(kPtVsPhi)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPhiVsEta, HistTable), getHistDesc(kPhiVsEta, HistTable), getHistType(kPhiVsEta, HistTable), {sigma0Specs.at(kPhiVsEta)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsLambdaMass, HistTable), getHistDesc(kPtVsLambdaMass, HistTable), getHistType(kPtVsLambdaMass, HistTable), {sigma0Specs.at(kPtVsLambdaMass)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsOpeningAngle, HistTable), getHistDesc(kPtVsOpeningAngle, HistTable), getHistType(kPtVsOpeningAngle, HistTable), {sigma0Specs.at(kPtVsOpeningAngle)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kLambdaPtVsPhotonPt, HistTable), getHistDesc(kLambdaPtVsPhotonPt, HistTable), getHistType(kLambdaPtVsPhotonPt, HistTable), {sigma0Specs.at(kLambdaPtVsPhotonPt)});
    }
  }

  void initMc(std::map<Sigma0Hist, std::vector<o2::framework::AxisSpec>> const& sigma0Specs)
  {
    std::string mcDir = std::string(sigma0Prefix) + std::string(McDir);
    mHistogramRegistry->add(mcDir + getHistNameV2(kTruePtVsPt, HistTable), getHistDesc(kTruePtVsPt, HistTable), getHistType(kTruePtVsPt, HistTable), {sigma0Specs.at(kTruePtVsPt)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kTrueEtaVsEta, HistTable), getHistDesc(kTrueEtaVsEta, HistTable), getHistType(kTrueEtaVsEta, HistTable), {sigma0Specs.at(kTrueEtaVsEta)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kTruePhiVsPhi, HistTable), getHistDesc(kTruePhiVsPhi, HistTable), getHistType(kTruePhiVsPhi, HistTable), {sigma0Specs.at(kTruePhiVsPhi)});

    // mc origin can be configured here
    const o2::framework::AxisSpec axisOrigin = {static_cast<int>(modes::McOrigin::kMcOriginLast), -0.5, static_cast<double>(modes::McOrigin::kMcOriginLast) - 0.5};
    mHistogramRegistry->add(mcDir + getHistNameV2(kOrigin, HistTable), getHistDesc(kOrigin, HistTable), getHistType(kOrigin, HistTable), {axisOrigin});
    mHistogramRegistry->get<TH1>(HIST(sigma0Prefix) + HIST(McDir) + HIST(histmanager::getHistName(kOrigin, HistTable)))->GetXaxis()->SetBinLabel(1 + static_cast<int>(modes::McOrigin::kNoMcParticle), modes::mcOriginToString(modes::McOrigin::kNoMcParticle));
    mHistogramRegistry->get<TH1>(HIST(sigma0Prefix) + HIST(McDir) + HIST(histmanager::getHistName(kOrigin, HistTable)))->GetXaxis()->SetBinLabel(1 + static_cast<int>(modes::McOrigin::kFromWrongCollision), modes::mcOriginToString(modes::McOrigin::kFromWrongCollision));
    mHistogramRegistry->get<TH1>(HIST(sigma0Prefix) + HIST(McDir) + HIST(histmanager::getHistName(kOrigin, HistTable)))->GetXaxis()->SetBinLabel(1 + static_cast<int>(modes::McOrigin::kPhysicalPrimary), modes::mcOriginToString(modes::McOrigin::kPhysicalPrimary));
    mHistogramRegistry->get<TH1>(HIST(sigma0Prefix) + HIST(McDir) + HIST(histmanager::getHistName(kOrigin, HistTable)))->GetXaxis()->SetBinLabel(1 + static_cast<int>(modes::McOrigin::kFromSecondaryDecay), modes::mcOriginToString(modes::McOrigin::kFromSecondaryDecay));
    mHistogramRegistry->get<TH1>(HIST(sigma0Prefix) + HIST(McDir) + HIST(histmanager::getHistName(kOrigin, HistTable)))->GetXaxis()->SetBinLabel(1 + static_cast<int>(modes::McOrigin::kFromMaterial), modes::mcOriginToString(modes::McOrigin::kFromMaterial));
    mHistogramRegistry->get<TH1>(HIST(sigma0Prefix) + HIST(McDir) + HIST(histmanager::getHistName(kOrigin, HistTable)))->GetXaxis()->SetBinLabel(1 + static_cast<int>(modes::McOrigin::kMissidentified), modes::mcOriginToString(modes::McOrigin::kMissidentified));

    mHistogramRegistry->add(mcDir + getHistNameV2(kPdg, HistTable), getHistDesc(kPdg, HistTable), getHistType(kPdg, HistTable), {sigma0Specs.at(kPdg)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kPdgMother, HistTable), getHistDesc(kPdgMother, HistTable), getHistType(kPdgMother, HistTable), {sigma0Specs.at(kPdgMother)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kPdgPartonicMother, HistTable), getHistDesc(kPdgPartonicMother, HistTable), getHistType(kPdgPartonicMother, HistTable), {sigma0Specs.at(kPdgPartonicMother)});
  }

  template <typename T>
  void fillAnalysis(T const& sigma0Candidate)
  {
    mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(AnalysisDir) + HIST(getHistName(kPt, HistTable)), sigma0Candidate.pt());
    mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(AnalysisDir) + HIST(getHistName(kEta, HistTable)), sigma0Candidate.eta());
    mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(AnalysisDir) + HIST(getHistName(kPhi, HistTable)), sigma0Candidate.phi());
    mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(AnalysisDir) + HIST(getHistName(kMass, HistTable)), sigma0Candidate.mass());
    mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(AnalysisDir) + HIST(getHistName(kSign, HistTable)), sigma0Candidate.sign());
    mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(AnalysisDir) + HIST(getHistName(kPtVsMass, HistTable)), sigma0Candidate.pt(), sigma0Candidate.mass());
  }

  // needs FSigma0Extras joined to the Sigma0 table
  template <typename T>
  void fillQa(T const& sigma0Candidate)
  {
    mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(QaDir) + HIST(getHistName(kLambdaMass, HistTable)), sigma0Candidate.lambdaMass());
    mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(QaDir) + HIST(getHistName(kLambdaPt, HistTable)), sigma0Candidate.lambdaPt());
    mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(QaDir) + HIST(getHistName(kLambdaCosPa, HistTable)), sigma0Candidate.lambdaCosPa());
    mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(QaDir) + HIST(getHistName(kLambdaDauDca, HistTable)), sigma0Candidate.lambdaDauDca());
    mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(QaDir) + HIST(getHistName(kLambdaTransRadius, HistTable)), sigma0Candidate.lambdaTransRadius());
    mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(QaDir) + HIST(getHistName(kPhotonPt, HistTable)), sigma0Candidate.photonPt());
    mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(QaDir) + HIST(getHistName(kPhotonCosPa, HistTable)), sigma0Candidate.photonCosPa());
    mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(QaDir) + HIST(getHistName(kPhotonTransRadius, HistTable)), sigma0Candidate.photonTransRadius());
    mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(QaDir) + HIST(getHistName(kPhotonChi2Ndf, HistTable)), sigma0Candidate.photonChi2Ndf());
    mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(QaDir) + HIST(getHistName(kOpeningAngle, HistTable)), sigma0Candidate.openingAngle());

    if (mPlot2d) {
      mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(QaDir) + HIST(getHistName(kPtVsEta, HistTable)), sigma0Candidate.pt(), sigma0Candidate.eta());
      mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(QaDir) + HIST(getHistName(kPtVsPhi, HistTable)), sigma0Candidate.pt(), sigma0Candidate.phi());
      mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(QaDir) + HIST(getHistName(kPhiVsEta, HistTable)), sigma0Candidate.phi(), sigma0Candidate.eta());
      mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(QaDir) + HIST(getHistName(kPtVsLambdaMass, HistTable)), sigma0Candidate.pt(), sigma0Candidate.lambdaMass());
      mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(QaDir) + HIST(getHistName(kPtVsOpeningAngle, HistTable)), sigma0Candidate.pt(), sigma0Candidate.openingAngle());
      mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(QaDir) + HIST(getHistName(kLambdaPtVsPhotonPt, HistTable)), sigma0Candidate.lambdaPt(), sigma0Candidate.photonPt());
    }
  }

  template <typename TMcParticles, typename TMcMothers, typename TMcPartonicMothers, typename T1, typename T2>
  void fillMc(T1 const& sigma0Candidate, T2 const& col)
  {
    // No MC Particle
    if (!sigma0Candidate.has_fMcParticle()) {
      mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(McDir) + HIST(getHistName(kPdg, HistTable)), 0);
      mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(McDir) + HIST(getHistName(kOrigin, HistTable)), static_cast<float>(modes::McOrigin::kNoMcParticle));
      return;
    }

    // Retrieve MC particle
    auto mcParticle = sigma0Candidate.template fMcParticle_as<TMcParticles>();

    // whether a particle is associated to a wrong collision or not cannot be known by the producer so we check it here
    bool fromWrongCollision = mcParticle.fMcColId() != col.fMcColId();

    // whether a particle is missidentfied or not cannot be known by the producer so we check it here
    bool isMissidentified = mcParticle.pdgCode() != mPdgCode;

    mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(McDir) + HIST(getHistName(kTruePtVsPt, HistTable)), mcParticle.pt(), sigma0Candidate.pt());
    mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(McDir) + HIST(getHistName(kTrueEtaVsEta, HistTable)), mcParticle.eta(), sigma0Candidate.eta());
    mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(McDir) + HIST(getHistName(kTruePhiVsPhi, HistTable)), mcParticle.phi(), sigma0Candidate.phi());
    if (fromWrongCollision) {
      mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(McDir) + HIST(getHistName(kOrigin, HistTable)), static_cast<float>(modes::McOrigin::kFromWrongCollision));
    } else if (isMissidentified) {
      mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(McDir) + HIST(getHistName(kOrigin, HistTable)), static_cast<float>(modes::McOrigin::kMissidentified));
    } else {
      mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(McDir) + HIST(getHistName(kOrigin, HistTable)), mcParticle.origin());
    }
    mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(McDir) + HIST(getHistName(kPdg, HistTable)), mcParticle.pdgCode());

    // get mother
    if (mcParticle.has_fMcMother()) {
      auto mother = mcParticle.template fMcMother_as<TMcMothers>();
      mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(McDir) + HIST(getHistName(kPdgMother, HistTable)), mother.pdgCode());
    } else {
      mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(McDir) + HIST(getHistName(kPdgMother, HistTable)), 0);
    }

    // get partonic mother
    if (mcParticle.has_fMcPartMoth()) {
      auto partonicMother = mcParticle.template fMcPartMoth_as<TMcPartonicMothers>();
      mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(McDir) + HIST(getHistName(kPdgPartonicMother, HistTable)), partonicMother.pdgCode());
    } else {
      mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(McDir) + HIST(getHistName(kPdgPartonicMother, HistTable)), 0);
    }
  }

  o2::framework::HistogramRegistry* mHistogramRegistry = nullptr;
  int mPdgCode = 0;
  int mPosDauPdgCodeAbs = 0;
  int mNegDauPdgCodeAbs = 0;
  bool mPlot2d = false;

  trackhistmanager::TrackHistManager<posDauPrefix> mPosDauManager;
  trackhistmanager::TrackHistManager<negDauPrefix> mNegDauManager;
};
} // namespace o2::analysis::femto::sigma0histmanager
#endif // PWGCF_FEMTO_CORE_SIGMA0HISTMANAGER_H_
