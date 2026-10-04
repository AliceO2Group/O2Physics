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

/// \file photonHistManager.h
/// \brief histogram manager for photon (PCM) histograms
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTO_CORE_PHOTONHISTMANAGER_H_
#define PWGCF_FEMTO_CORE_PHOTONHISTMANAGER_H_

#include "PWGCF/Femto/Core/histManager.h"
#include "PWGCF/Femto/Core/modes.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>

#include <map>
#include <string>
#include <string_view>
#include <vector>

namespace o2::analysis::femto::photonhistmanager
{
// enum for photon (PCM) histograms
enum PhotonHist {
  // analysis
  kPt,
  kEta,
  kPhi,
  // qa variables
  kCosPa,
  kPa,
  kDcaToPvXY,
  kDcaToPvZ,
  kChi2Ndf,
  kV0Radius,
  kDecayVtxX,
  kDecayVtxY,
  kDecayVtxZ,
  kDecayVtx,
  kPosDauTpcNSigmaEl,
  kNegDauTpcNSigmaEl,
  kPosDauPt,
  kNegDauPt,
  // 2d qa
  kPtVsEta,
  kPtVsPhi,
  kPhiVsEta,
  kPtVsCosPa,
  kPosDauVsNegDauTpcNSigmaEl,
  kPosDauTpcSignalVsP,
  kNegDauTpcSignalVsP,

  kPhotonHistLast
};

struct ConfPhotonBinning : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("PhotonBinning");
  o2::framework::ConfigurableAxis pt{"pt", {{600, 0, 6}}, "Pt"};
  o2::framework::ConfigurableAxis eta{"eta", {{300, -1.5, 1.5}}, "Eta"};
  o2::framework::ConfigurableAxis phi{"phi", {{720, 0, 1.f * o2::constants::math::TwoPI}}, "Phi"};
};

struct ConfPhotonQaBinning : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("PhotonQaBinning");
  o2::framework::Configurable<bool> plot2d{"plot2d", true, "Generate various 2D QA plots"};
  o2::framework::ConfigurableAxis cosPa{"cosPa", {{100, 0.8, 1}}, "Cosine of pointing angle"};
  o2::framework::ConfigurableAxis pa{"pa", {{180, 0, 1.f * o2::constants::math::PI}}, "Pointing angle"};
  o2::framework::ConfigurableAxis dcaToPv{"dcaToPv", {{200, -1, 1}}, "DCA of the photon to the PV (cm)"};
  o2::framework::ConfigurableAxis chi2Ndf{"chi2Ndf", {{100, 0, 50}}, "Chi2/NDF of the KF conversion vertex"};
  o2::framework::ConfigurableAxis v0Radius{"v0Radius", {{200, 0, 200}}, "Transverse radius of the conversion point (cm)"};
  o2::framework::ConfigurableAxis decayVertex{"decayVertex", {{200, 0, 200}}, "Distance of the conversion point from the nominal interaction point (cm)"};
  o2::framework::ConfigurableAxis tpcNSigmaEl{"tpcNSigmaEl", {{200, -10, 10}}, "TPC n sigma electron for daughters"};
  o2::framework::ConfigurableAxis dauPt{"dauPt", {{300, 0, 3}}, "Pt of the daughters"};
  o2::framework::ConfigurableAxis dauTpcInnerParam{"dauTpcInnerParam", {{200, 0, 4}}, "Momentum of the daughters at the inner wall of the TPC (GeV/#it{c})"};
  o2::framework::ConfigurableAxis dauTpcSignal{"dauTpcSignal", {{200, 0, 200}}, "TPC dE/dx of the daughters"};
};

// must be in sync with enum PhotonHist
constexpr std::array<histmanager::HistInfo<PhotonHist>, kPhotonHistLast> HistTable = {
  {{kPt, o2::framework::HistType::kTH1F, "hPt", "Transverse Momentum; p_{T} (GeV/#it{c}); Entries"},
   {kEta, o2::framework::HistType::kTH1F, "hEta", "Pseudorapdity; #eta; Entries"},
   {kPhi, o2::framework::HistType::kTH1F, "hPhi", "Azimuthal angle; #varphi; Entries"},
   {kCosPa, o2::framework::HistType::kTH1F, "hCosPa", "Cosine of pointing angle; cos(#alpha); Entries"},
   {kPa, o2::framework::HistType::kTH1F, "hPa", "Pointing angle; #alpha; Entries"},
   {kDcaToPvXY, o2::framework::HistType::kTH1F, "hDcaToPvXY", "DCA_{xy} to PV; DCA_{xy} (cm); Entries"},
   {kDcaToPvZ, o2::framework::HistType::kTH1F, "hDcaToPvZ", "DCA_{z} to PV; DCA_{z} (cm); Entries"},
   {kChi2Ndf, o2::framework::HistType::kTH1F, "hChi2Ndf", "Chi2/NDF of KF conversion vertex; #chi^{2}/NDF; Entries"},
   {kV0Radius, o2::framework::HistType::kTH1F, "hV0Radius", "Transverse radius; r_{xy} (cm); Entries"},
   {kDecayVtxX, o2::framework::HistType::kTH1F, "hDecayVtxX", "X coordinate of conversion point; DV_{X} (cm); Entries"},
   {kDecayVtxY, o2::framework::HistType::kTH1F, "hDecayVtxY", "Y coordinate of conversion point; DV_{Y} (cm); Entries"},
   {kDecayVtxZ, o2::framework::HistType::kTH1F, "hDecayVtxZ", "Z coordinate of conversion point; DV_{Z} (cm); Entries"},
   {kDecayVtx, o2::framework::HistType::kTH1F, "hDecayVtx", "Distance of conversion point from primary vertex; DV (cm); Entries"},
   {kPosDauTpcNSigmaEl, o2::framework::HistType::kTH1F, "hPosDauTpcNSigmaEl", "TPC n#sigma_{e} of positive daughter; n#sigma_{TPC,e}; Entries"},
   {kNegDauTpcNSigmaEl, o2::framework::HistType::kTH1F, "hNegDauTpcNSigmaEl", "TPC n#sigma_{e} of negative daughter; n#sigma_{TPC,e}; Entries"},
   {kPosDauPt, o2::framework::HistType::kTH1F, "hPosDauPt", "p_{T} of positive daughter; p_{T} (GeV/#it{c}); Entries"},
   {kNegDauPt, o2::framework::HistType::kTH1F, "hNegDauPt", "p_{T} of negative daughter; p_{T} (GeV/#it{c}); Entries"},
   {kPtVsEta, o2::framework::HistType::kTH2F, "hPtVsEta", "p_{T} vs #eta; p_{T} (GeV/#it{c}) ; #eta"},
   {kPtVsPhi, o2::framework::HistType::kTH2F, "hPtVsPhi", "p_{T} vs #varphi; p_{T} (GeV/#it{c}) ; #varphi"},
   {kPhiVsEta, o2::framework::HistType::kTH2F, "hPhiVsEta", "#varphi vs #eta; #varphi ; #eta"},
   {kPtVsCosPa, o2::framework::HistType::kTH2F, "hPtVsCosPa", "p_{T} vs cosine of pointing angle; p_{T} (GeV/#it{c}); cos(#alpha)"},
   {kPosDauVsNegDauTpcNSigmaEl, o2::framework::HistType::kTH2F, "hPosDauVsNegDauTpcNSigmaEl", "TPC n#sigma_{e} of positive vs negative daughter; n#sigma_{TPC,e} (pos.); n#sigma_{TPC,e} (neg.)"},
   {kPosDauTpcSignalVsP, o2::framework::HistType::kTH2F, "hPosDauTpcSignalVsP", "TPC dE/dx vs p of positive daughter; p_{TPC} (GeV/#it{c}); dE/dx"},
   {kNegDauTpcSignalVsP, o2::framework::HistType::kTH2F, "hNegDauTpcSignalVsP", "TPC dE/dx vs p of negative daughter; p_{TPC} (GeV/#it{c}); dE/dx"}},
};

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define PHOTON_HIST_ANALYSIS_MAP(conf) \
  {kPt, {(conf).pt}},                  \
    {kEta, {(conf).eta}},              \
    {kPhi, {(conf).phi}},

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define PHOTON_HIST_QA_MAP(confAnalysis, confQa)                                \
  {kCosPa, {(confQa).cosPa}},                                                   \
    {kPa, {(confQa).pa}},                                                      \
    {kDcaToPvXY, {(confQa).dcaToPv}},                                          \
    {kDcaToPvZ, {(confQa).dcaToPv}},                                           \
    {kChi2Ndf, {(confQa).chi2Ndf}},                                            \
    {kV0Radius, {(confQa).v0Radius}},                                         \
    {kDecayVtxX, {(confQa).decayVertex}},                                      \
    {kDecayVtxY, {(confQa).decayVertex}},                                      \
    {kDecayVtxZ, {(confQa).decayVertex}},                                      \
    {kDecayVtx, {(confQa).decayVertex}},                                       \
    {kPosDauTpcNSigmaEl, {(confQa).tpcNSigmaEl}},                              \
    {kNegDauTpcNSigmaEl, {(confQa).tpcNSigmaEl}},                              \
    {kPosDauPt, {(confQa).dauPt}},                                            \
    {kNegDauPt, {(confQa).dauPt}},                                            \
    {kPtVsEta, {(confAnalysis).pt, (confAnalysis).eta}},                       \
    {kPtVsPhi, {(confAnalysis).pt, (confAnalysis).phi}},                       \
    {kPhiVsEta, {(confAnalysis).phi, (confAnalysis).eta}},                     \
    {kPtVsCosPa, {(confAnalysis).pt, (confQa).cosPa}},                         \
    {kPosDauVsNegDauTpcNSigmaEl, {(confQa).tpcNSigmaEl, (confQa).tpcNSigmaEl}}, \
    {kPosDauTpcSignalVsP, {(confQa).dauTpcInnerParam, (confQa).dauTpcSignal}}, \
    {kNegDauTpcSignalVsP, {(confQa).dauTpcInnerParam, (confQa).dauTpcSignal}},

template <typename T>
auto makePhotonHistSpecMap(const T& confBinningAnalysis)
{
  return std::map<PhotonHist, std::vector<o2::framework::AxisSpec>>{
    PHOTON_HIST_ANALYSIS_MAP(confBinningAnalysis)};
}

template <typename T1, typename T2>
std::map<PhotonHist, std::vector<o2::framework::AxisSpec>> makePhotonQaHistSpecMap(T1 const& confBinningAnalysis, T2 const& confBinningQa)
{
  return std::map<PhotonHist, std::vector<o2::framework::AxisSpec>>{
    PHOTON_HIST_ANALYSIS_MAP(confBinningAnalysis)
      PHOTON_HIST_QA_MAP(confBinningAnalysis, confBinningQa)};
}

#undef PHOTON_HIST_ANALYSIS_MAP
#undef PHOTON_HIST_QA_MAP

constexpr char PrefixPhotonQa[] = "PhotonQA/";
constexpr char PrefixTwoPhotonResonanceDau1Qa[] = "TwoPhotonResonanceDau1Qa/";
constexpr char PrefixTwoPhotonResonanceDau2Qa[] = "TwoPhotonResonanceDau2Qa/";

constexpr std::string_view AnalysisDir = "Analysis/";
constexpr std::string_view QaDir = "QA/";

/// \class PhotonHistManager
/// \brief Histogram manager for PCM photon candidates (reco + QA only, no MC)
template <auto& Prefix>
class PhotonHistManager
{
 public:
  PhotonHistManager() = default;
  ~PhotonHistManager() = default;

  // init for analysis only
  template <modes::Mode mode>
  void init(o2::framework::HistogramRegistry* registry, std::map<PhotonHist, std::vector<o2::framework::AxisSpec>> const& specs)
  {
    mHistogramRegistry = registry;
    if constexpr (modes::isFlagSet(mode, modes::Mode::kReco)) {
      this->initAnalysis(specs);
    }
  }

  // init for analysis + qa
  template <modes::Mode mode, typename T>
  void init(o2::framework::HistogramRegistry* registry, std::map<PhotonHist, std::vector<o2::framework::AxisSpec>> const& specs, T const& confQaBinning)
  {
    mHistogramRegistry = registry;
    mPlot2d = confQaBinning.plot2d.value;
    if constexpr (modes::isFlagSet(mode, modes::Mode::kReco)) {
      this->initAnalysis(specs);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kQa)) {
      this->initQa(specs);
    }
  }

  template <modes::Mode mode, typename T>
  void fill(T const& photon)
  {
    if constexpr (modes::isFlagSet(mode, modes::Mode::kReco)) {
      this->fillAnalysis(photon);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kQa)) {
      this->fillQa(photon);
    }
  }

 private:
  void initAnalysis(std::map<PhotonHist, std::vector<o2::framework::AxisSpec>> const& specs)
  {
    std::string analysisDir = std::string(Prefix) + std::string(AnalysisDir);
    mHistogramRegistry->add(analysisDir + getHistNameV2(kPt, HistTable), getHistDesc(kPt, HistTable), getHistType(kPt, HistTable), {specs.at(kPt)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kEta, HistTable), getHistDesc(kEta, HistTable), getHistType(kEta, HistTable), {specs.at(kEta)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kPhi, HistTable), getHistDesc(kPhi, HistTable), getHistType(kPhi, HistTable), {specs.at(kPhi)});
  }

  void initQa(std::map<PhotonHist, std::vector<o2::framework::AxisSpec>> const& specs)
  {
    std::string qaDir = std::string(Prefix) + std::string(QaDir);
    mHistogramRegistry->add(qaDir + getHistNameV2(kCosPa, HistTable), getHistDesc(kCosPa, HistTable), getHistType(kCosPa, HistTable), {specs.at(kCosPa)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kPa, HistTable), getHistDesc(kPa, HistTable), getHistType(kPa, HistTable), {specs.at(kPa)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kDcaToPvXY, HistTable), getHistDesc(kDcaToPvXY, HistTable), getHistType(kDcaToPvXY, HistTable), {specs.at(kDcaToPvXY)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kDcaToPvZ, HistTable), getHistDesc(kDcaToPvZ, HistTable), getHistType(kDcaToPvZ, HistTable), {specs.at(kDcaToPvZ)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kChi2Ndf, HistTable), getHistDesc(kChi2Ndf, HistTable), getHistType(kChi2Ndf, HistTable), {specs.at(kChi2Ndf)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kV0Radius, HistTable), getHistDesc(kV0Radius, HistTable), getHistType(kV0Radius, HistTable), {specs.at(kV0Radius)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kDecayVtxX, HistTable), getHistDesc(kDecayVtxX, HistTable), getHistType(kDecayVtxX, HistTable), {specs.at(kDecayVtxX)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kDecayVtxY, HistTable), getHistDesc(kDecayVtxY, HistTable), getHistType(kDecayVtxY, HistTable), {specs.at(kDecayVtxY)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kDecayVtxZ, HistTable), getHistDesc(kDecayVtxZ, HistTable), getHistType(kDecayVtxZ, HistTable), {specs.at(kDecayVtxZ)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kDecayVtx, HistTable), getHistDesc(kDecayVtx, HistTable), getHistType(kDecayVtx, HistTable), {specs.at(kDecayVtx)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kPosDauTpcNSigmaEl, HistTable), getHistDesc(kPosDauTpcNSigmaEl, HistTable), getHistType(kPosDauTpcNSigmaEl, HistTable), {specs.at(kPosDauTpcNSigmaEl)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kNegDauTpcNSigmaEl, HistTable), getHistDesc(kNegDauTpcNSigmaEl, HistTable), getHistType(kNegDauTpcNSigmaEl, HistTable), {specs.at(kNegDauTpcNSigmaEl)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kPosDauPt, HistTable), getHistDesc(kPosDauPt, HistTable), getHistType(kPosDauPt, HistTable), {specs.at(kPosDauPt)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kNegDauPt, HistTable), getHistDesc(kNegDauPt, HistTable), getHistType(kNegDauPt, HistTable), {specs.at(kNegDauPt)});

    if (mPlot2d) {
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsEta, HistTable), getHistDesc(kPtVsEta, HistTable), getHistType(kPtVsEta, HistTable), {specs.at(kPtVsEta)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsPhi, HistTable), getHistDesc(kPtVsPhi, HistTable), getHistType(kPtVsPhi, HistTable), {specs.at(kPtVsPhi)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPhiVsEta, HistTable), getHistDesc(kPhiVsEta, HistTable), getHistType(kPhiVsEta, HistTable), {specs.at(kPhiVsEta)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsCosPa, HistTable), getHistDesc(kPtVsCosPa, HistTable), getHistType(kPtVsCosPa, HistTable), {specs.at(kPtVsCosPa)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPosDauVsNegDauTpcNSigmaEl, HistTable), getHistDesc(kPosDauVsNegDauTpcNSigmaEl, HistTable), getHistType(kPosDauVsNegDauTpcNSigmaEl, HistTable), {specs.at(kPosDauVsNegDauTpcNSigmaEl)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPosDauTpcSignalVsP, HistTable), getHistDesc(kPosDauTpcSignalVsP, HistTable), getHistType(kPosDauTpcSignalVsP, HistTable), {specs.at(kPosDauTpcSignalVsP)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kNegDauTpcSignalVsP, HistTable), getHistDesc(kNegDauTpcSignalVsP, HistTable), getHistType(kNegDauTpcSignalVsP, HistTable), {specs.at(kNegDauTpcSignalVsP)});
    }
  }

  template <typename T>
  void fillAnalysis(T const& photon)
  {
    mHistogramRegistry->fill(HIST(Prefix) + HIST(AnalysisDir) + HIST(getHistName(kPt, HistTable)), photon.pt());
    mHistogramRegistry->fill(HIST(Prefix) + HIST(AnalysisDir) + HIST(getHistName(kEta, HistTable)), photon.eta());
    mHistogramRegistry->fill(HIST(Prefix) + HIST(AnalysisDir) + HIST(getHistName(kPhi, HistTable)), photon.phi());
  }

  template <typename T>
  void fillQa(T const& photon)
  {
    mHistogramRegistry->fill(HIST(Prefix) + HIST(QaDir) + HIST(getHistName(kCosPa, HistTable)), photon.cosPa());
    mHistogramRegistry->fill(HIST(Prefix) + HIST(QaDir) + HIST(getHistName(kPa, HistTable)), photon.pa());
    mHistogramRegistry->fill(HIST(Prefix) + HIST(QaDir) + HIST(getHistName(kDcaToPvXY, HistTable)), photon.dcaToPvXY());
    mHistogramRegistry->fill(HIST(Prefix) + HIST(QaDir) + HIST(getHistName(kDcaToPvZ, HistTable)), photon.dcaToPvZ());
    mHistogramRegistry->fill(HIST(Prefix) + HIST(QaDir) + HIST(getHistName(kChi2Ndf, HistTable)), photon.chi2Ndf());
    mHistogramRegistry->fill(HIST(Prefix) + HIST(QaDir) + HIST(getHistName(kV0Radius, HistTable)), photon.v0Radius());
    mHistogramRegistry->fill(HIST(Prefix) + HIST(QaDir) + HIST(getHistName(kDecayVtxX, HistTable)), photon.decayVtxX());
    mHistogramRegistry->fill(HIST(Prefix) + HIST(QaDir) + HIST(getHistName(kDecayVtxY, HistTable)), photon.decayVtxY());
    mHistogramRegistry->fill(HIST(Prefix) + HIST(QaDir) + HIST(getHistName(kDecayVtxZ, HistTable)), photon.decayVtxZ());
    mHistogramRegistry->fill(HIST(Prefix) + HIST(QaDir) + HIST(getHistName(kDecayVtx, HistTable)), photon.decayVtx());
    mHistogramRegistry->fill(HIST(Prefix) + HIST(QaDir) + HIST(getHistName(kPosDauTpcNSigmaEl, HistTable)), photon.posDauTpcNSigmaEl());
    mHistogramRegistry->fill(HIST(Prefix) + HIST(QaDir) + HIST(getHistName(kNegDauTpcNSigmaEl, HistTable)), photon.negDauTpcNSigmaEl());
    mHistogramRegistry->fill(HIST(Prefix) + HIST(QaDir) + HIST(getHistName(kPosDauPt, HistTable)), photon.posDauPt());
    mHistogramRegistry->fill(HIST(Prefix) + HIST(QaDir) + HIST(getHistName(kNegDauPt, HistTable)), photon.negDauPt());

    if (mPlot2d) {
      mHistogramRegistry->fill(HIST(Prefix) + HIST(QaDir) + HIST(getHistName(kPtVsEta, HistTable)), photon.pt(), photon.eta());
      mHistogramRegistry->fill(HIST(Prefix) + HIST(QaDir) + HIST(getHistName(kPtVsPhi, HistTable)), photon.pt(), photon.phi());
      mHistogramRegistry->fill(HIST(Prefix) + HIST(QaDir) + HIST(getHistName(kPhiVsEta, HistTable)), photon.phi(), photon.eta());
      mHistogramRegistry->fill(HIST(Prefix) + HIST(QaDir) + HIST(getHistName(kPtVsCosPa, HistTable)), photon.pt(), photon.cosPa());
      mHistogramRegistry->fill(HIST(Prefix) + HIST(QaDir) + HIST(getHistName(kPosDauVsNegDauTpcNSigmaEl, HistTable)), photon.posDauTpcNSigmaEl(), photon.negDauTpcNSigmaEl());
      mHistogramRegistry->fill(HIST(Prefix) + HIST(QaDir) + HIST(getHistName(kPosDauTpcSignalVsP, HistTable)), photon.posDauTpcInnerParam(), photon.posDauTpcSignal());
      mHistogramRegistry->fill(HIST(Prefix) + HIST(QaDir) + HIST(getHistName(kNegDauTpcSignalVsP, HistTable)), photon.negDauTpcInnerParam(), photon.negDauTpcSignal());
    }
  }

  o2::framework::HistogramRegistry* mHistogramRegistry = nullptr;
  bool mPlot2d = false;
};
} // namespace o2::analysis::femto::photonhistmanager
#endif // PWGCF_FEMTO_CORE_PHOTONHISTMANAGER_H_
