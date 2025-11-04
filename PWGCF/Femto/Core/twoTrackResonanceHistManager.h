// Copyright 2019-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file twoTrackResonanceHistManager.h
/// \brief histogram manager for two track resonances
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#ifndef PWGCF_FEMTO_CORE_TWOTRACKRESONANCEHISTMANAGER_H_
#define PWGCF_FEMTO_CORE_TWOTRACKRESONANCEHISTMANAGER_H_

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
namespace twotrackresonancehistmanager
{
// enum for track histograms
enum TwoTrackResonanceHist {
  // analysis
  kPt,
  kEta,
  kPhi,
  kMass,
  kSign,
  // 2d qa
  kPtVsEta,
  kPtVsPhi,
  kPhiVsEta,
  kPtVsMass,
  kTwoTrackResonanceHistLast
};

#define TWOTRACKRESONANCE_DEFAULT_BINNING(defaultMassMin, defaultMassMax)                          \
  o2::framework::ConfigurableAxis pt{"pt", {{600, 0, 6}}, "Pt"};                                   \
  o2::framework::ConfigurableAxis eta{"eta", {{300, -1.5, 1.5}}, "Eta"};                           \
  o2::framework::ConfigurableAxis phi{"phi", {{720, 0, 1.f * o2::constants::math::TwoPI}}, "Phi"}; \
  o2::framework::ConfigurableAxis mass{"mass", {{200, defaultMassMin, defaultMassMax}}, "Mass"};   \
  o2::framework::ConfigurableAxis sign{"sign", {{3, -1.5, 1.5}}, "Sign"};

struct ConfPhiBinning : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("PhiBinning");
  TWOTRACKRESONANCE_DEFAULT_BINNING(0.8f, 1.2f)
};

struct ConfRho0Binning : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("Rho0Binning");
  TWOTRACKRESONANCE_DEFAULT_BINNING(0.5f, 1.f)
};

struct ConfKstar0Binning : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("Kstar0Binning");
  TWOTRACKRESONANCE_DEFAULT_BINNING(0.6f, 1.f)
};
#undef TWOTRACKRESONANCE_DEFAULT_BINNING

constexpr std::array<histmanager::HistInfo<TwoTrackResonanceHist>, kTwoTrackResonanceHistLast> HistTable = {
  {{kPt, o2::framework::kTH1F, "hPt", "Transverse Momentum; p_{T} (GeV/#it{c}); Entries"},
   {kEta, o2::framework::kTH1F, "hEta", "Pseudorapdity; #eta; Entries"},
   {kPhi, o2::framework::kTH1F, "hPhi", "Azimuthal angle; #varphi; Entries"},
   {kMass, o2::framework::kTH1F, "hMass", "Invariant mass; m (GeV/#it{c}^{2}); Entries"},
   {kSign, o2::framework::kTH1F, "hSign", "Sign (-1 -> antiparticle, 0 -> self conjugate, +1 -> particle); sign; Entries"},
   {kPtVsEta, o2::framework::kTH2F, "hPtVsEta", "p_{T} vs #eta; p_{T} (GeV/#it{c}) ; #eta"},
   {kPtVsPhi, o2::framework::kTH2F, "hPtVsPhi", "p_{T} vs #varphi;p_{T} (GeV/#it{c});#varphi"},
   {kPhiVsEta, o2::framework::kTH2F, "hPhiVsEta", "#varphi vs #eta; #varphi ; #eta"},
   {kPtVsMass, o2::framework::kTH2F, "hPtVsMass", "p_{T} vs invariant mass; p_{T} (GeV/#it{c}); m (GeV/#it{c}^{2})"}}};

template <typename T>
std::map<TwoTrackResonanceHist, std::vector<framework::AxisSpec>> makeTwoTrackResonanceHistSpecMap(const T& confBinningAnalysis)
{
  return std::map<TwoTrackResonanceHist, std::vector<framework::AxisSpec>>{
    {kPt, {confBinningAnalysis.pt}},
    {kEta, {confBinningAnalysis.eta}},
    {kPhi, {confBinningAnalysis.phi}},
    {kMass, {confBinningAnalysis.mass}},
    {kSign, {confBinningAnalysis.sign}}};
};

template <typename T>
auto makeTwoTrackResonanceQaHistSpecMap(const T& confBinningAnalysis)
{
  return std::map<TwoTrackResonanceHist, std::vector<framework::AxisSpec>>{
    {kPt, {confBinningAnalysis.pt}},
    {kEta, {confBinningAnalysis.eta}},
    {kPhi, {confBinningAnalysis.phi}},
    {kMass, {confBinningAnalysis.mass}},
    {kSign, {confBinningAnalysis.sign}},
    {kPtVsEta, {confBinningAnalysis.pt, confBinningAnalysis.eta}},
    {kPtVsPhi, {confBinningAnalysis.pt, confBinningAnalysis.phi}},
    {kPhiVsEta, {confBinningAnalysis.phi, confBinningAnalysis.eta}},
    {kPtVsMass, {confBinningAnalysis.pt, confBinningAnalysis.mass}}};
};

constexpr char PrefixRho[] = "Rho0/";
constexpr char PrefixPhi[] = "Phi/";
constexpr char PrefixKstar[] = "Kstar0/";

constexpr std::string_view AnalysisDir = "Kinematics/";
constexpr std::string_view QaDir = "QA/";

template <const char* resoPrefix, const char* posDauPrefix, const char* negDauPrefix, modes::Mode mode, modes::TwoTrackResonance reso>
class TwoTrackResonanceHistManager
{
 public:
  TwoTrackResonanceHistManager() = default;
  ~TwoTrackResonanceHistManager() = default;

  void init(o2::framework::HistogramRegistry* registry,
            std::map<TwoTrackResonanceHist, std::vector<o2::framework::AxisSpec>> const& ResoSpecs,
            std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> const& PosDauSpecs,
            std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> const& NegDauSpecs)
  {
    mHistogramRegistry = registry;
    mPosDauManager.init(registry, PosDauSpecs);
    mNegDauManager.init(registry, NegDauSpecs);
    if constexpr (modes::isFlagSet(mode, modes::Mode::kAnalysis)) {
      initAnalysis(ResoSpecs);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kQa)) {
      initQa(ResoSpecs);
    }
  }

  template <typename T1, typename T2>
  void enableOptionalHistograms(T1 const& PosDauConfBinningQa, T2 const& NegDauConfBinningQa)
  {
    mPosDauManager.enableOptionalHistograms(PosDauConfBinningQa);
    mNegDauManager.enableOptionalHistograms(NegDauConfBinningQa);
  }

  template <typename T1, typename T2>
  void init(o2::framework::HistogramRegistry* registry,
            std::map<TwoTrackResonanceHist, std::vector<o2::framework::AxisSpec>> ResoSpecs,
            std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> PosDauSpecs,
            T1 const& PosDauConfBinningQa,
            std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> NegDauSpecs,
            T2 const& NegDauConfBinningQa)
  {
    enableOptionalHistograms(PosDauConfBinningQa, NegDauConfBinningQa);
    init(registry, ResoSpecs, PosDauSpecs, NegDauSpecs);
  }

  template <typename T1, typename T2>
  void fill(T1 const& resonance, T2 const& tracks)
  {
    // this used to work, still under investigation
    // auto posDaughter = resonance.template posDau_as<T2>();
    // auto negDaughter = resonance.template negDau_as<T2>();
    auto posDaughter = tracks.rawIteratorAt(resonance.posDauId() - tracks.offset());
    mPosDauManager.fill(posDaughter, tracks);
    auto negDaughter = tracks.rawIteratorAt(resonance.negDauId() - tracks.offset());
    mNegDauManager.fill(negDaughter, tracks);
    if constexpr (modes::isFlagSet(mode, modes::Mode::kAnalysis)) {
      fillAnalysis(resonance);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kQa)) {
      fillQa(resonance);
    }
  }

 private:
  void initAnalysis(std::map<TwoTrackResonanceHist, std::vector<o2::framework::AxisSpec>> const& ResoSpecs)
  {
    std::string analysisDir = std::string(resoPrefix) + std::string(AnalysisDir);
    mHistogramRegistry->add(analysisDir + getHistNameV2(kPt, HistTable), getHistDesc(kPt, HistTable), getHistType(kPt, HistTable), {ResoSpecs.at(kPt)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kEta, HistTable), getHistDesc(kEta, HistTable), getHistType(kEta, HistTable), {ResoSpecs.at(kEta)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kPhi, HistTable), getHistDesc(kPhi, HistTable), getHistType(kPhi, HistTable), {ResoSpecs.at(kPhi)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kMass, HistTable), getHistDesc(kMass, HistTable), getHistType(kMass, HistTable), {ResoSpecs.at(kMass)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kSign, HistTable), getHistDesc(kSign, HistTable), getHistType(kSign, HistTable), {ResoSpecs.at(kSign)});
  }
  void initQa(std::map<TwoTrackResonanceHist, std::vector<o2::framework::AxisSpec>> const& ResoSpecs)
  {
    std::string qaDir = std::string(resoPrefix) + std::string(QaDir);
    mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsEta, HistTable), getHistDesc(kPtVsEta, HistTable), getHistType(kPtVsEta, HistTable), {ResoSpecs.at(kPtVsEta)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsPhi, HistTable), getHistDesc(kPtVsPhi, HistTable), getHistType(kPtVsPhi, HistTable), {ResoSpecs.at(kPtVsPhi)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kPhiVsEta, HistTable), getHistDesc(kPhiVsEta, HistTable), getHistType(kPhiVsEta, HistTable), {ResoSpecs.at(kPhiVsEta)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsMass, HistTable), getHistDesc(kPtVsMass, HistTable), getHistType(kPtVsMass, HistTable), {ResoSpecs.at(kPtVsMass)});
  }

  template <typename T>
  void fillAnalysis(T const& resonance)
  {
    mHistogramRegistry->fill(HIST(resoPrefix) + HIST(AnalysisDir) + HIST(getHistName(kPt, HistTable)), resonance.pt());
    mHistogramRegistry->fill(HIST(resoPrefix) + HIST(AnalysisDir) + HIST(getHistName(kEta, HistTable)), resonance.eta());
    mHistogramRegistry->fill(HIST(resoPrefix) + HIST(AnalysisDir) + HIST(getHistName(kPhi, HistTable)), resonance.phi());
    mHistogramRegistry->fill(HIST(resoPrefix) + HIST(AnalysisDir) + HIST(getHistName(kMass, HistTable)), resonance.mass());
    if constexpr (modes::isEqual(reso, modes::TwoTrackResonance::kPhi) || modes::isEqual(reso, modes::TwoTrackResonance::kRho0)) {
      mHistogramRegistry->fill(HIST(resoPrefix) + HIST(AnalysisDir) + HIST(getHistName(kSign, HistTable)), 0);
    }
    if constexpr (modes::isEqual(reso, modes::TwoTrackResonance::kKstar0) || modes::isEqual(reso, modes::TwoTrackResonance::kKstar0Bar)) {
      mHistogramRegistry->fill(HIST(resoPrefix) + HIST(AnalysisDir) + HIST(getHistName(kSign, HistTable)), resonance.sign());
    }
  }

  template <typename T>
  void fillQa(T const& resonance)
  {
    mHistogramRegistry->fill(HIST(resoPrefix) + HIST(QaDir) + HIST(getHistName(kPtVsEta, HistTable)), resonance.pt(), resonance.eta());
    mHistogramRegistry->fill(HIST(resoPrefix) + HIST(QaDir) + HIST(getHistName(kPtVsPhi, HistTable)), resonance.pt(), resonance.phi());
    mHistogramRegistry->fill(HIST(resoPrefix) + HIST(QaDir) + HIST(getHistName(kPhiVsEta, HistTable)), resonance.phi(), resonance.eta());
    mHistogramRegistry->fill(HIST(resoPrefix) + HIST(QaDir) + HIST(getHistName(kPtVsMass, HistTable)), resonance.pt(), resonance.mass());
  }

  o2::framework::HistogramRegistry* mHistogramRegistry = nullptr;
  trackhistmanager::TrackHistManager<posDauPrefix, mode> mPosDauManager;
  trackhistmanager::TrackHistManager<negDauPrefix, mode> mNegDauManager;
};
}; // namespace twotrackresonancehistmanager
}; // namespace o2::analysis::femto
#endif // PWGCF_FEMTO_CORE_TWOTRACKRESONANCEHISTMANAGER_H_
