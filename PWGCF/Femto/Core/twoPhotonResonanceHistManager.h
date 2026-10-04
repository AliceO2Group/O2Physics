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

/// \file twoPhotonResonanceHistManager.h
/// \brief histogram manager for resonances built from two photons (pi0, eta, ...)
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTO_CORE_TWOPHOTONRESONANCEHISTMANAGER_H_
#define PWGCF_FEMTO_CORE_TWOPHOTONRESONANCEHISTMANAGER_H_

#include "PWGCF/Femto/Core/histManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/photonHistManager.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>

#include <array>
#include <map>
#include <string>
#include <string_view>
#include <vector>

namespace o2::analysis::femto::twophotonresonancehistmanager
{
enum TwoPhotonResonanceHist {
  // analysis
  kPt,
  kEta,
  kPhi,
  kMass,
  kPtVsMass,
  // 2d qa
  kPtVsEta,
  kPtVsPhi,
  kPhiVsEta,
  kTwoPhotonResonanceHistLast
};

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define TWOPHOTONRESONANCE_DEFAULT_BINNING(defaultMassMin, defaultMassMax)                          \
  o2::framework::ConfigurableAxis pt{"pt", {{600, 0, 6}}, "Pt"};                                    \
  o2::framework::ConfigurableAxis eta{"eta", {{300, -1.5, 1.5}}, "Eta"};                            \
  o2::framework::ConfigurableAxis phi{"phi", {{720, 0, 1.f * o2::constants::math::TwoPI}}, "Phi"};  \
  o2::framework::ConfigurableAxis mass{"mass", {{200, (defaultMassMin), (defaultMassMax)}}, "Mass"};

struct ConfPi0Binning : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("Pi0Binning");
  TWOPHOTONRESONANCE_DEFAULT_BINNING(0.f, 0.3f)
};

struct ConfEtaBinning : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("EtaBinning");
  TWOPHOTONRESONANCE_DEFAULT_BINNING(0.2f, 0.9f)
};
#undef TWOPHOTONRESONANCE_DEFAULT_BINNING

constexpr std::array<histmanager::HistInfo<TwoPhotonResonanceHist>, kTwoPhotonResonanceHistLast> HistTable = {
  {{kPt, o2::framework::HistType::kTH1F, "hPt", "Transverse Momentum; p_{T} (GeV/#it{c}); Entries"},
   {kEta, o2::framework::HistType::kTH1F, "hEta", "Pseudorapdity; #eta; Entries"},
   {kPhi, o2::framework::HistType::kTH1F, "hPhi", "Azimuthal angle; #varphi; Entries"},
   {kMass, o2::framework::HistType::kTH1F, "hMass", "Invariant mass; m_{#gamma#gamma} (GeV/#it{c}^{2}); Entries"},
   {kPtVsMass, o2::framework::HistType::kTH2F, "hPtVsMass", "p_{T} vs invariant mass; p_{T} (GeV/#it{c}); m_{#gamma#gamma} (GeV/#it{c}^{2})"},
   {kPtVsEta, o2::framework::HistType::kTH2F, "hPtVsEta", "p_{T} vs #eta; p_{T} (GeV/#it{c}) ; #eta"},
   {kPtVsPhi, o2::framework::HistType::kTH2F, "hPtVsPhi", "p_{T} vs #varphi;p_{T} (GeV/#it{c});#varphi"},
   {kPhiVsEta, o2::framework::HistType::kTH2F, "hPhiVsEta", "#varphi vs #eta; #varphi ; #eta"}},
};

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define TWOPHOTONRESONANCE_HIST_ANALYSIS_MAP(conf) \
  {kPt, {(conf).pt}},                              \
    {kEta, {(conf).eta}},                          \
    {kPhi, {(conf).phi}},                          \
    {kMass, {(conf).mass}},                        \
    {kPtVsMass, {(conf).pt, (conf).mass}},

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define TWOPHOTONRESONANCE_HIST_QA_MAP(conf)  \
  {kPtVsEta, {(conf).pt, (conf).eta}},        \
    {kPtVsPhi, {(conf).pt, (conf).phi}},      \
    {kPhiVsEta, {(conf).phi, (conf).eta}},    \
    {kPtVsMass, {(conf).pt, (conf).mass}},

template <typename T>
std::map<TwoPhotonResonanceHist, std::vector<o2::framework::AxisSpec>> makeTwoPhotonResonanceHistSpecMap(const T& confBinningAnalysis)
{
  return std::map<TwoPhotonResonanceHist, std::vector<o2::framework::AxisSpec>>{
    TWOPHOTONRESONANCE_HIST_ANALYSIS_MAP(confBinningAnalysis)};
}

template <typename T>
auto makeTwoPhotonResonanceQaHistSpecMap(const T& confBinningAnalysis)
{
  return std::map<TwoPhotonResonanceHist, std::vector<o2::framework::AxisSpec>>{
    TWOPHOTONRESONANCE_HIST_ANALYSIS_MAP(confBinningAnalysis)
      TWOPHOTONRESONANCE_HIST_QA_MAP(confBinningAnalysis)};
}

#undef TWOPHOTONRESONANCE_HIST_ANALYSIS_MAP
#undef TWOPHOTONRESONANCE_HIST_QA_MAP

constexpr char PrefixPi0[] = "Pi0/";
constexpr char PrefixEta[] = "Eta/";

constexpr std::string_view AnalysisDir = "Analysis/";
constexpr std::string_view QaDir = "QA/";

/// \class TwoPhotonResonanceHistManager
/// \brief Histogram manager for resonances built from two PCM photons (reco + QA only, no MC).
///        Both daughters share the same PhotonHistManager-based sub-manager treatment: there is no
///        pos/neg distinction, nor a PDG-dependent daughter lookup the way TwoTrackResonanceHistManager
///        needs (photon daughters are always massless, regardless of the mother species).
template <auto& resoPrefix,
          auto& dau1Prefix,
          auto& dau2Prefix,
          modes::TwoPhotonResonance reso>
class TwoPhotonResonanceHistManager
{
 public:
  TwoPhotonResonanceHistManager() = default;
  ~TwoPhotonResonanceHistManager() = default;

  template <modes::Mode mode, typename T1>
  void init(o2::framework::HistogramRegistry* registry,
            std::map<TwoPhotonResonanceHist, std::vector<o2::framework::AxisSpec>> const& ResoSpecs,
            std::map<photonhistmanager::PhotonHist, std::vector<o2::framework::AxisSpec>> const& DauSpecs,
            T1 const& ConfDauQaBinning)
  {
    mHistogramRegistry = registry;
    mDau1Manager.template init<mode>(registry, DauSpecs, ConfDauQaBinning);
    mDau2Manager.template init<mode>(registry, DauSpecs, ConfDauQaBinning);
    if constexpr (modes::isFlagSet(mode, modes::Mode::kReco)) {
      initAnalysis(ResoSpecs);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kQa)) {
      initQa(ResoSpecs);
    }
  }

  template <modes::Mode mode, typename T1, typename T2>
  void fill(T1 const& resonance, T2 const& photons)
  {
    auto dau1 = photons.rawIteratorAt(resonance.dau1Id() - photons.offset());
    mDau1Manager.template fill<mode>(dau1);
    auto dau2 = photons.rawIteratorAt(resonance.dau2Id() - photons.offset());
    mDau2Manager.template fill<mode>(dau2);
    if constexpr (modes::isFlagSet(mode, modes::Mode::kReco)) {
      fillAnalysis(resonance);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kQa)) {
      fillQa(resonance);
    }
  }

 private:
  void initAnalysis(std::map<TwoPhotonResonanceHist, std::vector<o2::framework::AxisSpec>> const& ResoSpecs)
  {
    std::string analysisDir = std::string(resoPrefix) + std::string(AnalysisDir);
    mHistogramRegistry->add(analysisDir + getHistNameV2(kPt, HistTable), getHistDesc(kPt, HistTable), getHistType(kPt, HistTable), {ResoSpecs.at(kPt)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kEta, HistTable), getHistDesc(kEta, HistTable), getHistType(kEta, HistTable), {ResoSpecs.at(kEta)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kPhi, HistTable), getHistDesc(kPhi, HistTable), getHistType(kPhi, HistTable), {ResoSpecs.at(kPhi)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kMass, HistTable), getHistDesc(kMass, HistTable), getHistType(kMass, HistTable), {ResoSpecs.at(kMass)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kPtVsMass, HistTable), getHistDesc(kPtVsMass, HistTable), getHistType(kPtVsMass, HistTable), {ResoSpecs.at(kPtVsMass)});
  }

  void initQa(std::map<TwoPhotonResonanceHist, std::vector<o2::framework::AxisSpec>> const& ResoSpecs)
  {
    std::string qaDir = std::string(resoPrefix) + std::string(QaDir);
    mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsEta, HistTable), getHistDesc(kPtVsEta, HistTable), getHistType(kPtVsEta, HistTable), {ResoSpecs.at(kPtVsEta)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsPhi, HistTable), getHistDesc(kPtVsPhi, HistTable), getHistType(kPtVsPhi, HistTable), {ResoSpecs.at(kPtVsPhi)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kPhiVsEta, HistTable), getHistDesc(kPhiVsEta, HistTable), getHistType(kPhiVsEta, HistTable), {ResoSpecs.at(kPhiVsEta)});
  }

  template <typename T>
  void fillAnalysis(T const& resonance)
  {
    mHistogramRegistry->fill(HIST(resoPrefix) + HIST(AnalysisDir) + HIST(getHistName(kPt, HistTable)), resonance.pt());
    mHistogramRegistry->fill(HIST(resoPrefix) + HIST(AnalysisDir) + HIST(getHistName(kEta, HistTable)), resonance.eta());
    mHistogramRegistry->fill(HIST(resoPrefix) + HIST(AnalysisDir) + HIST(getHistName(kPhi, HistTable)), resonance.phi());
    mHistogramRegistry->fill(HIST(resoPrefix) + HIST(AnalysisDir) + HIST(getHistName(kMass, HistTable)), resonance.mass());
    mHistogramRegistry->fill(HIST(resoPrefix) + HIST(AnalysisDir) + HIST(getHistName(kPtVsMass, HistTable)), resonance.pt(), resonance.mass());
  }

  template <typename T>
  void fillQa(T const& resonance)
  {
    mHistogramRegistry->fill(HIST(resoPrefix) + HIST(QaDir) + HIST(getHistName(kPtVsEta, HistTable)), resonance.pt(), resonance.eta());
    mHistogramRegistry->fill(HIST(resoPrefix) + HIST(QaDir) + HIST(getHistName(kPtVsPhi, HistTable)), resonance.pt(), resonance.phi());
    mHistogramRegistry->fill(HIST(resoPrefix) + HIST(QaDir) + HIST(getHistName(kPhiVsEta, HistTable)), resonance.phi(), resonance.eta());
  }

  o2::framework::HistogramRegistry* mHistogramRegistry = nullptr;
  photonhistmanager::PhotonHistManager<dau1Prefix> mDau1Manager;
  photonhistmanager::PhotonHistManager<dau2Prefix> mDau2Manager;
};
} // namespace o2::analysis::femto::twophotonresonancehistmanager
#endif // PWGCF_FEMTO_CORE_TWOPHOTONRESONANCEHISTMANAGER_H_
