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
/// \brief histogram manager for sigm histograms
/// \author Yakiv Paroviak, CTU in Prague, yakiv.paroviak@cern.ch

#ifndef PWGCF_FEMTO_CORE_SIGMA0HISTMANAGER_H_
#define PWGCF_FEMTO_CORE_SIGMA0HISTMANAGER_H_

#include "PWGCF/Femto/Core/histManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/trackHistManager.h"
#include "PWGCF/Femto/Core/v0HistManager.h"
#include "PWGCF/Femto/Core/photonHistManager.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/Logger.h>

#include <TH1.h>
#include <TPDGCode.h>

#include <array>
#include <cstddef>
#include <map>
#include <string>
#include <string_view>
#include <vector>

namespace o2::analysis::femto::sigma0histmanager
{
// enum for histograms
enum Sigma0Hist {
  // analysis
  kPt,
  kEta,
  kPhi,
  kMass,
  kSign,
  kPtVsMass,
  // qa
  kMassSigma0,
  kMassAntiSigma0,
  kCosPa,
  kDecayVtxX,
  kDecayVtxY,
  kDecayVtxZ,
  kDecayVtx,
  // 2d qa
  kPtVsEta,
  kPtVsPhi,
  kPhiVsEta,
  kPtVsCosPa,
  kPtVsSigma0Mass,
  kPtVsAntiSigma0Mass,
  kSigma0MassVsAntiSigma0Mass,
  // mc
  kOrigin,
  kPdg,
  kTruePtVsPt,
  kTrueEtaVsEta,
  kTruePhiVsPhi,
  // histograms for fraction estimation of sigma0s
  kNoMcParticle,
  kPrimary,

  kSigma0HistLast
};

constexpr std::size_t MaxSecondary = 3;

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define SIGMA0_DEFAULT_BINNING(defaultMassMin, defaultMassMax)                                       \
  o2::framework::ConfigurableAxis pt{"pt", {{600, 0, 6}}, "Pt"};                                     \
  o2::framework::ConfigurableAxis eta{"eta", {{300, -1.5, 1.5}}, "Eta"};                             \
  o2::framework::ConfigurableAxis phi{"phi", {{720, 0, 1.f * o2::constants::math::TwoPI}}, "Phi"};   \
  o2::framework::ConfigurableAxis mass{"mass", {{200, (defaultMassMin), (defaultMassMax)}}, "Mass"}; \
  o2::framework::ConfigurableAxis sign{"sign", {{3, -1.5, 1.5}}, "Sign"};                            \
  o2::framework::ConfigurableAxis pdgCodes{"pdgCodes", {{8001, -4000.5, 4000.5}}, "MC ONLY: PDG codes of reconstructed Sigma0s"};

template <auto& Prefix>
struct ConfSigma0Binning : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  SIGMA0_DEFAULT_BINNING(1.07, 1.32) // rough estimate of 5-sigma mass window for Sigma0
};

#undef SIGMA0_DEFAULT_BINNING

constexpr const char PrefixSigma0Binning[] = "Sigma0Binning";
using ConfSigma0BinningDefault = ConfSigma0Binning<PrefixSigma0Binning>;

template <auto& Prefix>
struct ConfSigma0QaBinning : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  o2::framework::Configurable<bool> plot2d{"plot2d", true, "Generate various 2D QA plots"};
  o2::framework::ConfigurableAxis cosPa{"cosPa", {{100, 0.9, 1}}, "Cosine of pointing angle"};
  o2::framework::ConfigurableAxis decayVertex{"decayVertex", {{100, 0, 100}}, "Decay vertex"};
  o2::framework::ConfigurableAxis massSigma0{"massSigma0", {{200, 1, 1.2}}, "mass for particle hypothesis"};
  o2::framework::ConfigurableAxis massAntiSigma0{"massAntiSigma0", {{100, 1, 1.2}}, "mass for antiparticle hypothesis"};
};

constexpr const char PrefixSigma0QaBinning[] = "Sigma0QaBinning";
using ConfSigma0QaBinningDefault = ConfSigma0QaBinning<PrefixSigma0QaBinning>;

// must be in sync with enum
constexpr std::array<histmanager::HistInfo<Sigma0Hist>, kSigma0HistLast> HistTable = {
  {{kPt, o2::framework::HistType::kTH1F, "hPt", "Transverse Momentum; p_{T} (GeV/#it{c}); Entries"},
   {kEta, o2::framework::HistType::kTH1F, "hEta", "Pseudorapdity; #eta; Entries"},
   {kPhi, o2::framework::HistType::kTH1F, "hPhi", "Azimuthal angle; #varphi; Entries"},
   {kMass, o2::framework::HistType::kTH1F, "hMass", "Invariant Mass; m_{Inv} (GeV/#it{c}^{2}); Entries"},
   {kSign, o2::framework::HistType::kTH1F, "hSign", "Sign (-1 -> antiparticle, 0 -> self conjugate, +1 -> particle); sign; Entries"},
   {kPtVsMass, o2::framework::HistType::kTH2F, "hPtVsMass", "Transverse momentum vs invariant mass; p_{T} (GeV/#it{c}); m_{Inv} (GeV/#it{c}^{2})"},
   {kMassSigma0, o2::framework::HistType::kTH1F, "hMassSigma0", "#Sigma0 mass; m_{p#pi^{-}} (GeV/#it{c}^{2}); Entries"},
   {kMassAntiSigma0, o2::framework::HistType::kTH1F, "hMassAntiSigma0", "#bar{#Sigma0} mass; m_{#bar{p}#pi^{+}} (GeV/#it{c}^{2}); Entries"},
   {kCosPa, o2::framework::HistType::kTH1F, "hCosPa", "Cosine of pointing angle; cos(#alpha); Entries"},
   {kDecayVtxX, o2::framework::HistType::kTH1F, "hDecayVtxX", "X coordinate of decay vertex ; DV_{X} (cm); Entries"},
   {kDecayVtxY, o2::framework::HistType::kTH1F, "hDecayVtxY", "Y coordinate of decay vertex ; DV_{Y} (cm); Entries"},
   {kDecayVtxZ, o2::framework::HistType::kTH1F, "hDecayVtxZ", "Z coordinate of decay vertex ; DV_{Z} (cm); Entries"},
   {kDecayVtx, o2::framework::HistType::kTH1F, "hDecayVtx", "Distance of decay vertex from primary vertex ; DV (cm); Entries"},
   {kPtVsEta, o2::framework::HistType::kTH2F, "hPtVsEta", "p_{T} vs #eta; p_{T} (GeV/#it{c}) ; #eta"},
   {kPtVsPhi, o2::framework::HistType::kTH2F, "hPtVsPhi", "p_{T} vs #varphi; p_{T} (GeV/#it{c}) ; #varphi"},
   {kPhiVsEta, o2::framework::HistType::kTH2F, "hPhiVsEta", "#varphi vs #eta; #varphi ; #eta"},
   {kPtVsCosPa, o2::framework::HistType::kTH2F, "hPtVsCosPa", "p_{T} vs Cosine of pointing angle;p_{T} (GeV/#it{c}); cos(#alpha);"},
   {kPtVsSigma0Mass, o2::framework::HistType::kTH2F, "hPtVsSigma0Mass", "p_{T} vs #Sigma0 mass; p_{T} (GeV/#it{c}); m_{p#pi^{-}} (GeV/#it{c}^{2})"},
   {kPtVsAntiSigma0Mass, o2::framework::HistType::kTH2F, "hPtVsAntiSigma0Mass", "p_{T} vs #bar{#Sigma0} mass; p_{T} (GeV/#it{c}); m_{#bar{p}#pi^{+}} (GeV/#it{c}^{2})"},
   {kSigma0MassVsAntiSigma0Mass, o2::framework::HistType::kTH2F, "hSigma0MassVsAntiSigma0Mass", "#Sigma0 mass vs #bar{#Sigma0}; m_{p#pi^{-}} (GeV/#it{c}^{2}); m_{#bar{p}#pi^{+}} (GeV/#it{c}^{2})"},
   {kOrigin, o2::framework::HistType::kTH1F, "hOrigin", "Status Codes (=Origin); Status Code; Entries"},
   {kPdg, o2::framework::HistType::kTH1F, "hPdg", "PDG Codes of reconstructed sigma0; PDG Code; Entries"},
   {kTruePtVsPt, o2::framework::HistType::kTH2F, "hTruePtVsPt", "True transverse momentum vs transverse momentum; p_{T,True} (GeV/#it{c}); p_{T,True} (GeV/#it{c})"},
   {kTrueEtaVsEta, o2::framework::HistType::kTH2F, "hTrueEtaVsEta", "True pseudorapdity vs pseudorapdity; #eta_{True}; #eta"},
   {kTruePhiVsPhi, o2::framework::HistType::kTH2F, "hTruePhiVsPhi", "True azimuthal angle vs azimuthal angle; #varphi_{True}; #varphi"},
   {kNoMcParticle, o2::framework::HistType::kTH2F, "hNoMcParticle", "Wrongly reconstructed particles; p_{T} (GeV/#it{c}); cos(#alpha)"},
   {kPrimary, o2::framework::HistType::kTH2F, "hPrimary", "Primary particles; p_{T} (GeV/#it{c}); cos(#alpha)"}},
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
    {kPdg, {(conf).pdgCodes}},                 

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define SIGMA0_HIST_QA_MAP(confAnalysis, confQa)                                     \
  {kCosPa, {(confQa).cosPa}},                                                        \
    {kDecayVtxX, {(confQa).decayVertex}},                                            \
    {kDecayVtxY, {(confQa).decayVertex}},                                            \
    {kDecayVtxZ, {(confQa).decayVertex}},                                            \
    {kDecayVtx, {(confQa).decayVertex}},                                             \
    {kPtVsEta, {(confAnalysis).pt, (confAnalysis).eta}},                             \
    {kPtVsPhi, {(confAnalysis).pt, (confAnalysis).phi}},                             \
    {kPhiVsEta, {(confAnalysis).phi, (confAnalysis).eta}},                           \
    {kPtVsCosPa, {(confAnalysis).pt, (confQa).cosPa}},                               \
    {kMassSigma0, {(confQa).massSigma0}},                                            \
    {kMassAntiSigma0, {(confQa).massAntiSigma0}},                                    \
    {kPtVsSigma0Mass, {(confAnalysis).pt, (confQa).massSigma0}},                     \
    {kPtVsAntiSigma0Mass, {(confAnalysis).pt, (confQa).massAntiSigma0}},             \
    {kSigma0MassVsAntiSigma0Mass, {(confQa).massSigma0, (confQa).massAntiSigma0}},

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define SIGMA0_HIST_MC_QA_MAP(confAnalysis, confQa)             \
  {kNoMcParticle, {(confAnalysis).pt, (confQa).cosPa}},         \
    {kPrimary, {(confAnalysis).pt, (confQa).cosPa}},            

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
        SIGMA0_HIST_MC_MAP(confBinningAnalysis)
          SIGMA0_HIST_MC_QA_MAP(confBinningAnalysis, confBinningQa)};
}

#undef SIGMA0_HIST_ANALYSIS_MAP
#undef SIGMA0_HIST_QA_MAP
#undef SIGMA0_HIST_MC_MAP
#undef SIGMA0_HIST_MC_QA_MAP

constexpr char PrefixSigma0Qa[] = "Sigma0QA/";
constexpr char PrefixSigma0[] = "Sigma0/";
constexpr char PrefixSigma0Cascade[] = "Sigma0CascadeQa/";

constexpr std::string_view AnalysisDir = "Analysis/";
constexpr std::string_view QaDir = "QA/";
constexpr std::string_view McDir = "MC/";

/// \class Sigma0HistManager
/// \brief Class for histogramming sigma0
template <auto& sigma0Prefix>
class Sigma0HistManager
{
 public:
  Sigma0HistManager() = default;
  ~Sigma0HistManager() = default;

    // minimal init
  template <modes::Mode mode, typename T1>
  void init(o2::framework::HistogramRegistry* registry,
            std::map<Sigma0Hist, std::vector<o2::framework::AxisSpec>> const& Sigma0Specs,
            T1 const& ConfSigma0Selection)
  {
    mHistogramRegistry = registry;
    mPdgCode = std::abs(ConfSigma0Selection.pdgCodeAbs.value);
    if (mPdgCode != PDG_t::kSigma0) {
      LOG(fatal) << "PDG code for Sigma0 has to be Sigma0";
    }
    if (ConfSigma0Selection.sign.value <= 0) {
      mPdgCode = -mPdgCode;
    }

    if constexpr (modes::isFlagSet(mode, modes::Mode::kReco)) {
      this->initAnalysis(Sigma0Specs);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kMc)) {
      this->initMc(Sigma0Specs);
    }
  }

  // init with qa
  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4>
  void init(o2::framework::HistogramRegistry* registry,
            std::map<Sigma0Hist, std::vector<o2::framework::AxisSpec>> const& Sigma0Specs,
            T1 const& ConfSigma0Selection,
            T2 const& ConfSigma0BinningQa)
  {
    this->template init<mode>(registry, Sigma0Specs, ConfSigma0Selection);
    this->enableOptionalHistograms(ConfSigma0BinningQa);
    if constexpr (modes::isFlagSet(mode, modes::Mode::kQa)) {
      this->initQa(Sigma0Specs);
    }
  }
    
  // init with qa + lambda and photon managers
  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4>
  void init(o2::framework::HistogramRegistry* registry,
            std::map<Sigma0Hist, std::vector<o2::framework::AxisSpec>> const& Sigma0Specs,
            T1 const& ConfSigma0Selection,
            T2 const& ConfSigma0BinningQa,
            std::map<o2::analysis::femto::v0histmanager::V0Hist, std::vector<o2::framework::AxisSpec>> const& LambdaSpecs,
            T3 const& ConfLambdaSelection,
            std::map<o2::analysis::femto::trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> const& PosDauSpecs,
            std::map<o2::analysis::femto::trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> const& NegDauSpecs,
            std::map<o2::analysis::femto::photonhistmanager::PhotonHist, std::vector<o2::framework::AxisSpec>> const& PhotonSpecs)
  {
    this->template init<mode>(registry, Sigma0Specs, ConfSigma0Selection, ConfSigma0BinningQa);
    mLambdaManager.template init<mode>(registry, LambdaSpecs, ConfLambdaSelection, PosDauSpecs, NegDauSpecs);
    mPhotonManager.template init<mode>(registry, PhotonSpecs);
  }

  // fill only sigma0
  template <modes::Mode mode, typename T1>
  void fill(T1 const& sigma0candidate)
  {
    if constexpr (modes::isFlagSet(mode, modes::Mode::kReco)) {
       this->fillAnalysis(sigma0candidate);
     }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kQa)) {
       this->fillQa(sigma0candidate);
    }
  }

    // wrapper for tripletProcessHelpers
  template <modes::Mode mode, typename T1, typename T2>
  void fill(T1 const& sigma0candidate, T2 const& tracks)
  {
       this->fillAnalysis(sigma0candidate);
  }

  // fill sigma0 + lambda + photon
  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4>
  void fill(T1 const& sigma0candidate, T2 const& v0table, T3 const& photontable, T4 const& tracks)
  {
    this->template fill<mode>(sigma0candidate);
    auto lambdau = v0table.rawIteratorAt(sigma0candidate.lambdaId() - v0table.offset());
    mLambdaManager.template fill<modes::Mode::kReco>(lambdau, tracks); // do not fill qa for daughters
    auto photondau = photontable.rawIteratorAt(sigma0candidate.photonId() - photontable.offset());
    mPhotonManager.template fill<modes::Mode::kReco>(photondau);
  }
  // fill mc
  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  void fill(T1 const& sigma0candidate, T2 const& col, T3 const& mcParticles, T4 const& mcMothers, T5 const& mcPartonicMothers)
  {
    this->template fill<mode>(sigma0candidate);
    if constexpr (modes::isFlagSet(mode, modes::Mode::kMc)) {
      this->template fillMc<mode>(sigma0candidate, col, mcParticles, mcMothers, mcPartonicMothers);
    }
  }
  // wrapper for tripletProcessHelpers
  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  void fill(T1 const& sigma0candidate, T2 const& trackTable, T3 const& col, T4 const& mcParticles, T5 const& mcMothers, T6 const& mcPartonicMothers)
  {
      this->template fillMc<mode>(sigma0candidate, col, mcParticles, mcMothers, mcPartonicMothers);
  }
  // fill mc + lambda + photon
  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
  void fill(T1 const& sigma0candidate, T2 const& col, T3 const& mcParticles, T4 const& mcMothers, T5 const& mcPartonicMothers, T6 const& v0table, T7 const& tracks)
  {
    this->template fillMc<mode>(sigma0candidate, col, mcParticles, mcMothers, mcPartonicMothers);
    auto lambdau = v0table.rawIteratorAt(sigma0candidate.lambdaId() - v0table.offset());
    mLambdaManager.template fill<modes::Mode::kReco_Mc>(lambdau, tracks, col, mcParticles, mcMothers, mcPartonicMothers);
    // photon MC manager is not implemented yet
    // auto photondau = photontable.rawIteratorAt(sigma0candidate.photonId() - photontable.offset());
    // mPhotonManager.template fill<mode>(photondau, col, mcParticles, mcMothers, mcPartonicMothers);
  }

 private:
  template <typename T1>
  void enableOptionalHistograms(T1 const& Sigma0ConfBinningQa)
  {
    mPlot2d = Sigma0ConfBinningQa.plot2d.value;
    mPlotOrigins = Sigma0ConfBinningQa.plotOrigins.value;
    mPlotNSecondaries = Sigma0ConfBinningQa.pdgCodesForMothersOfSecondary.value.size();
    for (std::size_t i = 0; i < MaxSecondary; i++) {
      if (i < Sigma0ConfBinningQa.pdgCodesForMothersOfSecondary.value.size()) {
        mPdgCodesSecondaryMother.at(i) = std::abs(Sigma0ConfBinningQa.pdgCodesForMothersOfSecondary.value.at(i));
      } else {
        mPdgCodesSecondaryMother.at(i) = 0;
      }
    }
  }

  void initAnalysis(std::map<Sigma0Hist, std::vector<o2::framework::AxisSpec>> const& Sigma0Specs)
  {
    std::string analysisDir = std::string(sigma0Prefix) + std::string(AnalysisDir);
    mHistogramRegistry->add(analysisDir + getHistNameV2(kPt, HistTable), getHistDesc(kPt, HistTable), getHistType(kPt, HistTable), {Sigma0Specs.at(kPt)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kEta, HistTable), getHistDesc(kEta, HistTable), getHistType(kEta, HistTable), {Sigma0Specs.at(kEta)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kPhi, HistTable), getHistDesc(kPhi, HistTable), getHistType(kPhi, HistTable), {Sigma0Specs.at(kPhi)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kMass, HistTable), getHistDesc(kMass, HistTable), getHistType(kMass, HistTable), {Sigma0Specs.at(kMass)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kSign, HistTable), getHistDesc(kSign, HistTable), getHistType(kSign, HistTable), {Sigma0Specs.at(kSign)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kPtVsMass, HistTable), getHistDesc(kPtVsMass, HistTable), getHistType(kPtVsMass, HistTable), {Sigma0Specs.at(kPtVsMass)});
  }

  void initQa(std::map<Sigma0Hist, std::vector<o2::framework::AxisSpec>> const& Sigma0Specs)
  {
    std::string qaDir = std::string(sigma0Prefix) + std::string(QaDir);

    mHistogramRegistry->add(qaDir + getHistNameV2(kCosPa, HistTable), getHistDesc(kCosPa, HistTable), getHistType(kCosPa, HistTable), {Sigma0Specs.at(kCosPa)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kDecayVtxX, HistTable), getHistDesc(kDecayVtxX, HistTable), getHistType(kDecayVtxX, HistTable), {Sigma0Specs.at(kDecayVtxX)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kDecayVtxY, HistTable), getHistDesc(kDecayVtxY, HistTable), getHistType(kDecayVtxY, HistTable), {Sigma0Specs.at(kDecayVtxY)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kDecayVtxZ, HistTable), getHistDesc(kDecayVtxZ, HistTable), getHistType(kDecayVtxZ, HistTable), {Sigma0Specs.at(kDecayVtxZ)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kDecayVtx, HistTable), getHistDesc(kDecayVtx, HistTable), getHistType(kDecayVtx, HistTable), {Sigma0Specs.at(kDecayVtx)});

    if (mPlot2d) {
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsEta, HistTable), getHistDesc(kPtVsEta, HistTable), getHistType(kPtVsEta, HistTable), {Sigma0Specs.at(kPtVsEta)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsPhi, HistTable), getHistDesc(kPtVsPhi, HistTable), getHistType(kPtVsPhi, HistTable), {Sigma0Specs.at(kPtVsPhi)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPhiVsEta, HistTable), getHistDesc(kPhiVsEta, HistTable), getHistType(kPhiVsEta, HistTable), {Sigma0Specs.at(kPhiVsEta)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsCosPa, HistTable), getHistDesc(kPtVsCosPa, HistTable), getHistType(kPtVsCosPa, HistTable), {Sigma0Specs.at(kPtVsCosPa)});

      mHistogramRegistry->add(qaDir + getHistNameV2(kMassSigma0, HistTable), getHistDesc(kMassSigma0, HistTable), getHistType(kMassSigma0, HistTable), {Sigma0Specs.at(kMassSigma0)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kMassAntiSigma0, HistTable), getHistDesc(kMassAntiSigma0, HistTable), getHistType(kMassAntiSigma0, HistTable), {Sigma0Specs.at(kMassAntiSigma0)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsSigma0Mass, HistTable), getHistDesc(kPtVsSigma0Mass, HistTable), getHistType(kPtVsSigma0Mass, HistTable), {Sigma0Specs.at(kPtVsSigma0Mass)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsAntiSigma0Mass, HistTable), getHistDesc(kPtVsAntiSigma0Mass, HistTable), getHistType(kPtVsAntiSigma0Mass, HistTable), {Sigma0Specs.at(kPtVsAntiSigma0Mass)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kSigma0MassVsAntiSigma0Mass, HistTable), getHistDesc(kSigma0MassVsAntiSigma0Mass, HistTable), getHistType(kSigma0MassVsAntiSigma0Mass, HistTable), {Sigma0Specs.at(kSigma0MassVsAntiSigma0Mass)});
    }
  }

  void initMc(std::map<Sigma0Hist, std::vector<o2::framework::AxisSpec>> const& Sigma0Specs)
  {
    std::string mcDir = std::string(sigma0Prefix) + std::string(McDir);
    mHistogramRegistry->add(mcDir + getHistNameV2(kTruePtVsPt, HistTable), getHistDesc(kTruePtVsPt, HistTable), getHistType(kTruePtVsPt, HistTable), {Sigma0Specs.at(kTruePtVsPt)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kTrueEtaVsEta, HistTable), getHistDesc(kTrueEtaVsEta, HistTable), getHistType(kTrueEtaVsEta, HistTable), {Sigma0Specs.at(kTrueEtaVsEta)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kTruePhiVsPhi, HistTable), getHistDesc(kTruePhiVsPhi, HistTable), getHistType(kTruePhiVsPhi, HistTable), {Sigma0Specs.at(kTruePhiVsPhi)});

    // mc origin can be configured here
    const o2::framework::AxisSpec axisOrigin = {static_cast<int>(modes::McOrigin::kMcOriginLast), -0.5, static_cast<double>(modes::McOrigin::kMcOriginLast) - 0.5};
    mHistogramRegistry->add(mcDir + getHistNameV2(kOrigin, HistTable), getHistDesc(kOrigin, HistTable), getHistType(kOrigin, HistTable), {axisOrigin});
    mHistogramRegistry->get<TH1>(HIST(sigma0Prefix) + HIST(McDir) + HIST(histmanager::getHistName(kOrigin, HistTable)))->GetXaxis()->SetBinLabel(1 + static_cast<int>(modes::McOrigin::kNoMcParticle), modes::mcOriginToString(modes::McOrigin::kNoMcParticle));
    mHistogramRegistry->get<TH1>(HIST(sigma0Prefix) + HIST(McDir) + HIST(histmanager::getHistName(kOrigin, HistTable)))->GetXaxis()->SetBinLabel(1 + static_cast<int>(modes::McOrigin::kFromWrongCollision), modes::mcOriginToString(modes::McOrigin::kFromWrongCollision));
    mHistogramRegistry->get<TH1>(HIST(sigma0Prefix) + HIST(McDir) + HIST(histmanager::getHistName(kOrigin, HistTable)))->GetXaxis()->SetBinLabel(1 + static_cast<int>(modes::McOrigin::kPhysicalPrimary), modes::mcOriginToString(modes::McOrigin::kPhysicalPrimary));

    mHistogramRegistry->add(mcDir + getHistNameV2(kPdg, HistTable), getHistDesc(kPdg, HistTable), getHistType(kPdg, HistTable), {Sigma0Specs.at(kPdg)});

    if (mPlotOrigins) {
      mHistogramRegistry->add(mcDir + getHistNameV2(kNoMcParticle, HistTable), getHistDesc(kNoMcParticle, HistTable), getHistType(kNoMcParticle, HistTable), {Sigma0Specs.at(kNoMcParticle)});
      mHistogramRegistry->add(mcDir + getHistNameV2(kPrimary, HistTable), getHistDesc(kPrimary, HistTable), getHistType(kPrimary, HistTable), {Sigma0Specs.at(kPrimary)});
    }
  }

  template <typename T>
  void fillAnalysis(T const& sigma0candidate)
  {
    mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(AnalysisDir) + HIST(getHistName(kPt, HistTable)), sigma0candidate.pt());
    mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(AnalysisDir) + HIST(getHistName(kEta, HistTable)), sigma0candidate.eta());
    mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(AnalysisDir) + HIST(getHistName(kPhi, HistTable)), sigma0candidate.phi());
    mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(AnalysisDir) + HIST(getHistName(kMass, HistTable)), sigma0candidate.mass());
    mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(AnalysisDir) + HIST(getHistName(kPtVsMass, HistTable)), sigma0candidate.pt(), sigma0candidate.mass());
    mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(AnalysisDir) + HIST(getHistName(kSign, HistTable)), sigma0candidate.sign());
  }

  template <typename T1>
  void fillQa(T1 const& sigma0candidate)
  {
    mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(QaDir) + HIST(getHistName(kCosPa, HistTable)), sigma0candidate.cosPa());
    mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(QaDir) + HIST(getHistName(kDecayVtxX, HistTable)), sigma0candidate.decayVtxX());
    mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(QaDir) + HIST(getHistName(kDecayVtxY, HistTable)), sigma0candidate.decayVtxY());
    mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(QaDir) + HIST(getHistName(kDecayVtxZ, HistTable)), sigma0candidate.decayVtxZ());
    mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(QaDir) + HIST(getHistName(kDecayVtx, HistTable)), sigma0candidate.decayVtx());

    float massSigma0 = 0;
    float massAntiSigma0 = 0;

    mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(QaDir) + HIST(getHistName(kMassSigma0, HistTable)), massSigma0);
    mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(QaDir) + HIST(getHistName(kMassAntiSigma0, HistTable)), massAntiSigma0);

    if (mPlot2d) {
      mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(QaDir) + HIST(getHistName(kPtVsEta, HistTable)), sigma0candidate.pt(), sigma0candidate.eta());
      mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(QaDir) + HIST(getHistName(kPtVsPhi, HistTable)), sigma0candidate.pt(), sigma0candidate.phi());
      mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(QaDir) + HIST(getHistName(kPhiVsEta, HistTable)), sigma0candidate.phi(), sigma0candidate.eta());
      mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(QaDir) + HIST(getHistName(kPtVsCosPa, HistTable)), sigma0candidate.pt(), sigma0candidate.cosPa());

      mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(QaDir) + HIST(getHistName(kPtVsSigma0Mass, HistTable)), sigma0candidate.pt(), massSigma0);
      mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(QaDir) + HIST(getHistName(kPtVsAntiSigma0Mass, HistTable)), sigma0candidate.pt(), massAntiSigma0);
      mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(QaDir) + HIST(getHistName(kSigma0MassVsAntiSigma0Mass, HistTable)), massSigma0, massAntiSigma0);
    }
  }

  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5>
  void fillMc(T1 const& sigma0Candidate, T2 const& col, T3 const& /*mcParticles*/, T4 const& /*mcMothers*/, T5 const& /*mcPartonicMothers*/)
  {
    // No MC Particle
    if (!sigma0Candidate.has_fMcParticle()) {
      mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(McDir) + HIST(getHistName(kPdg, HistTable)), 0);
      mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(McDir) + HIST(getHistName(kOrigin, HistTable)), static_cast<float>(modes::McOrigin::kNoMcParticle));
      if constexpr (modes::isFlagSet(mode, modes::Mode::kQa)) {
        if (mPlotOrigins) {
          mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(McDir) + HIST(getHistName(kNoMcParticle, HistTable)), sigma0Candidate.pt(), sigma0Candidate.cosPa());
        }
      }
      return;
    }

    // Retrieve MC particle
    auto mcParticle = sigma0Candidate.template fMcParticle_as<T3>();

    // whether a particle is associated to a wrong collision or not cannot be known by the producer so we check it here
    bool fromWrongCollision = mcParticle.fMcColId() != col.fMcColId();

    mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(McDir) + HIST(getHistName(kTruePtVsPt, HistTable)), mcParticle.pt(), sigma0Candidate.pt());
    mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(McDir) + HIST(getHistName(kTrueEtaVsEta, HistTable)), mcParticle.eta(), sigma0Candidate.eta());
    mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(McDir) + HIST(getHistName(kTruePhiVsPhi, HistTable)), mcParticle.phi(), sigma0Candidate.phi());
    if (!fromWrongCollision) {
      mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(McDir) + HIST(getHistName(kOrigin, HistTable)), mcParticle.origin());
    }
    mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(McDir) + HIST(getHistName(kPdg, HistTable)), mcParticle.pdgCode());

    if constexpr (modes::isFlagSet(mode, modes::Mode::kQa)) {
      if (mPlotOrigins) {
        // check first if particle is missidentified
        if (!fromWrongCollision) {
          // if not, we fill it acccoridng to its origin
          switch (static_cast<modes::McOrigin>(mcParticle.origin())) {
            case modes::McOrigin::kPhysicalPrimary:
              mHistogramRegistry->fill(HIST(sigma0Prefix) + HIST(McDir) + HIST(getHistName(kPrimary, HistTable)), sigma0Candidate.pt(), sigma0Candidate.cosPa());
              break;
            default:
              LOG(warn) << "Encountered particle with unknown origin!";
              break;
          }
        }
      }
    }
  }

  o2::framework::HistogramRegistry* mHistogramRegistry = nullptr;
  int mPdgCode = 0;
  bool mPlotOrigins = false;
  int mPlotNSecondaries = 0;
  std::array<int, MaxSecondary> mPdgCodesSecondaryMother = {0};
  bool mPlot2d = false;

  inline static constexpr std::string_view PrefixLambda = "Sigma0/DauLambda/";
  inline static constexpr std::string_view PrefixPos = "Sigma0/DauLambda/PosDau/";
  inline static constexpr std::string_view PrefixNeg = "Sigma0/DauLambda/NegDau/";
  inline static constexpr std::string_view PrefixPhoton = "Sigma0/DauPhoton/";
  v0histmanager::V0HistManager<PrefixLambda, PrefixPos, PrefixNeg, modes::V0::kLambda> mLambdaManager;
  photonhistmanager::PhotonHistManager<PrefixPhoton> mPhotonManager;
};
}; // namespace o2::analysis::femto::sigma0histmanager
#endif // PWGCF_FEMTO_CORE_Sigma0HISTMANAGER_H_
