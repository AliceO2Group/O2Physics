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

/// \file trackTrackTwoPhotonResonanceHistManager.h
/// \brief histogram manager for resonances built from two tracks and a two-photon resonance
///        (e.g. omega(782) -> pi+ pi- pi0)
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTO_CORE_TRACKTRACKTWOPHOTONRESONANCEHISTMANAGER_H_
#define PWGCF_FEMTO_CORE_TRACKTRACKTWOPHOTONRESONANCEHISTMANAGER_H_

#include "PWGCF/Femto/Core/histManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/photonHistManager.h"
#include "PWGCF/Femto/Core/trackHistManager.h"
#include "PWGCF/Femto/Core/twoPhotonResonanceHistManager.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>

#include <TPDGCode.h>

#include <array>
#include <map>
#include <string>
#include <string_view>
#include <vector>

namespace o2::analysis::femto::trackstwophotonresonancehistmanager
{
enum TrackTrackTwoPhotonResonanceHist {
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
  kTrackTrackTwoPhotonResonanceHistLast
};

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define TRACKTRACKTWOPHOTONRESONANCE_DEFAULT_BINNING(defaultMassMin, defaultMassMax)               \
  o2::framework::ConfigurableAxis pt{"pt", {{600, 0, 6}}, "Pt"};                                   \
  o2::framework::ConfigurableAxis eta{"eta", {{300, -1.5, 1.5}}, "Eta"};                           \
  o2::framework::ConfigurableAxis phi{"phi", {{720, 0, 1.f * o2::constants::math::TwoPI}}, "Phi"}; \
  o2::framework::ConfigurableAxis mass{"mass", {{200, (defaultMassMin), (defaultMassMax)}}, "Mass"};

struct ConfOmegaMesonBinning : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("OmegaMesonBinning");
  TRACKTRACKTWOPHOTONRESONANCE_DEFAULT_BINNING(0.6f, 0.95f)
};
#undef TRACKTRACKTWOPHOTONRESONANCE_DEFAULT_BINNING

constexpr std::array<histmanager::HistInfo<TrackTrackTwoPhotonResonanceHist>, kTrackTrackTwoPhotonResonanceHistLast> HistTable = {
  {{kPt, o2::framework::HistType::kTH1F, "hPt", "Transverse Momentum; p_{T} (GeV/#it{c}); Entries"},
   {kEta, o2::framework::HistType::kTH1F, "hEta", "Pseudorapdity; #eta; Entries"},
   {kPhi, o2::framework::HistType::kTH1F, "hPhi", "Azimuthal angle; #varphi; Entries"},
   {kMass, o2::framework::HistType::kTH1F, "hMass", "Invariant mass; m (GeV/#it{c}^{2}); Entries"},
   {kPtVsMass, o2::framework::HistType::kTH2F, "hPtVsMass", "p_{T} vs invariant mass; p_{T} (GeV/#it{c}); m (GeV/#it{c}^{2})"},
   {kPtVsEta, o2::framework::HistType::kTH2F, "hPtVsEta", "p_{T} vs #eta; p_{T} (GeV/#it{c}) ; #eta"},
   {kPtVsPhi, o2::framework::HistType::kTH2F, "hPtVsPhi", "p_{T} vs #varphi;p_{T} (GeV/#it{c});#varphi"},
   {kPhiVsEta, o2::framework::HistType::kTH2F, "hPhiVsEta", "#varphi vs #eta; #varphi ; #eta"}}};

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define TRACKTRACKTWOPHOTONRESONANCE_HIST_ANALYSIS_MAP(conf) \
  {kPt, {(conf).pt}},                                        \
    {kEta, {(conf).eta}},                                    \
    {kPhi, {(conf).phi}},                                    \
    {kMass, {(conf).mass}},                                  \
    {kPtVsMass, {(conf).pt, (conf).mass}},

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define TRACKTRACKTWOPHOTONRESONANCE_HIST_QA_MAP(conf) \
  {kPtVsEta, {(conf).pt, (conf).eta}},                 \
    {kPtVsPhi, {(conf).pt, (conf).phi}},               \
    {kPhiVsEta, {(conf).phi, (conf).eta}},

template <typename T>
auto makeTrackTrackTwoPhotonResonanceHistSpecMap(const T& confBinningAnalysis)
{
  return std::map<TrackTrackTwoPhotonResonanceHist, std::vector<o2::framework::AxisSpec>>{
    TRACKTRACKTWOPHOTONRESONANCE_HIST_ANALYSIS_MAP(confBinningAnalysis)};
}

template <typename T>
auto makeTrackTrackTwoPhotonResonanceQaHistSpecMap(const T& confBinningAnalysis)
{
  return std::map<TrackTrackTwoPhotonResonanceHist, std::vector<o2::framework::AxisSpec>>{
    TRACKTRACKTWOPHOTONRESONANCE_HIST_ANALYSIS_MAP(confBinningAnalysis)
      TRACKTRACKTWOPHOTONRESONANCE_HIST_QA_MAP(confBinningAnalysis)};
}

#undef TRACKTRACKTWOPHOTONRESONANCE_HIST_ANALYSIS_MAP
#undef TRACKTRACKTWOPHOTONRESONANCE_HIST_QA_MAP

constexpr char PrefixOmegaMeson[] = "OmegaMeson/";
constexpr char PrefixOmegaMesonNeutralDau[] = "OmegaMesonNeutralDau/";
constexpr char PrefixOmegaMesonNeutralDauPhotonDau1[] = "OmegaMesonNeutralDauPhotonDau1Qa/";
constexpr char PrefixOmegaMesonNeutralDauPhotonDau2[] = "OmegaMesonNeutralDauPhotonDau2Qa/";

constexpr std::string_view AnalysisDir = "Analysis/";
constexpr std::string_view QaDir = "QA/";

/// \class TrackTrackTwoPhotonResonanceHistManager
/// \brief Histogram manager for resonances built from two tracks and a two-photon resonance (e.g.
///        omega(782) -> pi+ pi- pi0). The charged daughters reuse the regular TrackHistManager;
///        the embedded two-photon-resonance daughter reuses TwoPhotonResonanceHistManager as-is
///        (same reco+QA histograms, including its own photon-leg QA, that the dedicated Pi0/Eta
///        QA task gets) by also threading the FPhotons table through to it.
template <auto& resoPrefix,
          auto& posDauPrefix,
          auto& negDauPrefix,
          auto& neutralDauPrefix,
          auto& neutralDauPhotonDau1Prefix,
          auto& neutralDauPhotonDau2Prefix,
          modes::TrackTrackTwoPhotonResonance reso,
          modes::TwoPhotonResonance neutralDauReso>
class TrackTrackTwoPhotonResonanceHistManager
{
 public:
  TrackTrackTwoPhotonResonanceHistManager() = default;
  ~TrackTrackTwoPhotonResonanceHistManager() = default;

  template <modes::Mode mode, typename T1, typename T2, typename T3>
  void init(o2::framework::HistogramRegistry* registry,
            std::map<TrackTrackTwoPhotonResonanceHist, std::vector<o2::framework::AxisSpec>> const& ResoSpecs,
            std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> const& PosDauSpecs,
            std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> const& NegDauSpecs,
            std::map<twophotonresonancehistmanager::TwoPhotonResonanceHist, std::vector<o2::framework::AxisSpec>> const& NeutralDauSpecs,
            std::map<photonhistmanager::PhotonHist, std::vector<o2::framework::AxisSpec>> const& NeutralDauPhotonSpecs,
            T1 const& ConfPosDauQaBinning,
            T2 const& ConfNegDauQaBinning,
            T3 const& ConfNeutralDauPhotonQaBinning)
  {
    mHistogramRegistry = registry;

    // omega(782) -> pi+ pi- pi0: both charged daughters are pions
    mPosDauManager.template init<mode>(registry, PosDauSpecs, 1, 1, kPiPlus, ConfPosDauQaBinning);
    mNegDauManager.template init<mode>(registry, NegDauSpecs, 1, -1, kPiPlus, ConfNegDauQaBinning);

    // the embedded pi0/eta candidate is a fully-fledged femto table with its own mask/extras
    // columns, so unlike the plain-track daughters above, it can just follow the outer mode.
    mNeutralDauManager.template init<mode>(registry, NeutralDauSpecs, NeutralDauPhotonSpecs, ConfNeutralDauPhotonQaBinning);

    if constexpr (modes::isFlagSet(mode, modes::Mode::kReco)) {
      initAnalysis(ResoSpecs);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kQa)) {
      initQa(ResoSpecs);
    }
  }

  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4>
  void fill(T1 const& resonance, T2 const& tracks, T3 const& neutralDaus, T4 const& photons)
  {
    auto posDau = tracks.rawIteratorAt(resonance.posDauId() - tracks.offset());
    mPosDauManager.template fill<mode>(posDau, tracks);
    auto negDau = tracks.rawIteratorAt(resonance.negDauId() - tracks.offset());
    mNegDauManager.template fill<mode>(negDau, tracks);
    auto neutralDau = neutralDaus.rawIteratorAt(resonance.neutralDauId() - neutralDaus.offset());
    mNeutralDauManager.template fill<mode>(neutralDau, photons);

    if constexpr (modes::isFlagSet(mode, modes::Mode::kReco)) {
      fillAnalysis(resonance);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kQa)) {
      fillQa(resonance);
    }
  }

 private:
  void initAnalysis(std::map<TrackTrackTwoPhotonResonanceHist, std::vector<o2::framework::AxisSpec>> const& ResoSpecs)
  {
    std::string analysisDir = std::string(resoPrefix) + std::string(AnalysisDir);
    mHistogramRegistry->add(analysisDir + getHistNameV2(kPt, HistTable), getHistDesc(kPt, HistTable), getHistType(kPt, HistTable), {ResoSpecs.at(kPt)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kEta, HistTable), getHistDesc(kEta, HistTable), getHistType(kEta, HistTable), {ResoSpecs.at(kEta)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kPhi, HistTable), getHistDesc(kPhi, HistTable), getHistType(kPhi, HistTable), {ResoSpecs.at(kPhi)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kMass, HistTable), getHistDesc(kMass, HistTable), getHistType(kMass, HistTable), {ResoSpecs.at(kMass)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kPtVsMass, HistTable), getHistDesc(kPtVsMass, HistTable), getHistType(kPtVsMass, HistTable), {ResoSpecs.at(kPtVsMass)});
  }

  void initQa(std::map<TrackTrackTwoPhotonResonanceHist, std::vector<o2::framework::AxisSpec>> const& ResoSpecs)
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
  trackhistmanager::TrackHistManager<posDauPrefix> mPosDauManager;
  trackhistmanager::TrackHistManager<negDauPrefix> mNegDauManager;
  twophotonresonancehistmanager::TwoPhotonResonanceHistManager<neutralDauPrefix, neutralDauPhotonDau1Prefix, neutralDauPhotonDau2Prefix, neutralDauReso> mNeutralDauManager;
};
} // namespace o2::analysis::femto::trackstwophotonresonancehistmanager
#endif // PWGCF_FEMTO_CORE_TRACKTRACKTWOPHOTONRESONANCEHISTMANAGER_H_
