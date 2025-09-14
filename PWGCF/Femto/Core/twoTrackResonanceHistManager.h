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
    {kSign, {confBinningAnalysis.sing}}};
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

/// \class FemtoDreamEventHisto
/// \brief Class for histogramming event properties
// template <femtomodes::Mode mode>
template <const char* resoPrefix, const char* posDauPrefix, const char* negDauPrefix, modes::Mode mode, modes::TwoTrackResonance reso>
class TwoTrackResonanceHistManager
{
 public:
  /// Destructor
  virtual ~TwoTrackResonanceHistManager() = default;
  /// Initializes histograms for the task
  /// \param registry Histogram registry to be passed
  ///
  void init(o2::framework::HistogramRegistry* registry,
            std::map<TwoTrackResonanceHist, std::vector<o2::framework::AxisSpec>> ResoSpecs,
            std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> PosDauSpecs,
            std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> NegDauSpecs)
  {
    mHistogramRegistry = registry;
    mPosDauManager.init(registry, PosDauSpecs);
    mNegDauManager.init(registry, NegDauSpecs);

    if constexpr (modes::isFlagSet(mode, modes::Mode::kAnalysis)) {
      std::string analysisDir = std::string(resoPrefix) + std::string(AnalysisDir);

      mHistogramRegistry->add(analysisDir + GetHistNamev2(kPt, HistTable), GetHistDesc(kPt, HistTable), GetHistType(kPt, HistTable), {ResoSpecs[kPt]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kEta, HistTable), GetHistDesc(kEta, HistTable), GetHistType(kEta, HistTable), {ResoSpecs[kEta]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kPhi, HistTable), GetHistDesc(kPhi, HistTable), GetHistType(kPhi, HistTable), {ResoSpecs[kPhi]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kMass, HistTable), GetHistDesc(kMass, HistTable), GetHistType(kMass, HistTable), {ResoSpecs[kMass]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kSign, HistTable), GetHistDesc(kSign, HistTable), GetHistType(kSign, HistTable), {ResoSpecs[kSign]});
    }

    if constexpr (modes::isFlagSet(mode, modes::Mode::kQa)) {
      std::string qaDir = std::string(resoPrefix) + std::string(QaDir);

      mHistogramRegistry->add(qaDir + GetHistNamev2(kPtVsEta, HistTable), GetHistDesc(kPtVsEta, HistTable), GetHistType(kPtVsEta, HistTable), {ResoSpecs[kPtVsEta]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPtVsPhi, HistTable), GetHistDesc(kPtVsPhi, HistTable), GetHistType(kPtVsPhi, HistTable), {ResoSpecs[kPtVsPhi]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPhiVsEta, HistTable), GetHistDesc(kPhiVsEta, HistTable), GetHistType(kPhiVsEta, HistTable), {ResoSpecs[kPhiVsEta]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPtVsMass, HistTable), GetHistDesc(kPtVsMass, HistTable), GetHistType(kPtVsMass, HistTable), {ResoSpecs[kPtVsMass]});
    }
  }

  template <typename T1, typename T2>
  void fill(T1 const& resonance, T2 const& /*tracks*/)
  {
    auto posDaughter = resonance.template posDau_as<T2>();
    mPosDauManager.fill(posDaughter);
    auto negDaughter = resonance.template negDau_as<T2>();
    mNegDauManager.fill(negDaughter);

    if constexpr (modes::isFlagSet(mode, modes::Mode::kAnalysis)) {
      mHistogramRegistry->fill(HIST(resoPrefix) + HIST(AnalysisDir) + HIST(GetHistName(kPt, HistTable)), resonance.pt());
      mHistogramRegistry->fill(HIST(resoPrefix) + HIST(AnalysisDir) + HIST(GetHistName(kEta, HistTable)), resonance.eta());
      mHistogramRegistry->fill(HIST(resoPrefix) + HIST(AnalysisDir) + HIST(GetHistName(kPhi, HistTable)), resonance.phi());
      mHistogramRegistry->fill(HIST(resoPrefix) + HIST(AnalysisDir) + HIST(GetHistName(kMass, HistTable)), resonance.mass());

      if constexpr (modes::isEqual(reso, modes::TwoTrackResonance::kPhi) || modes::isEqual(reso, modes::TwoTrackResonance::kRho0)) {
        mHistogramRegistry->fill(HIST(resoPrefix) + HIST(AnalysisDir) + HIST(GetHistName(kSign, HistTable)), 0);
      }
      if constexpr (modes::isEqual(reso, modes::TwoTrackResonance::kKstar0) || modes::isEqual(reso, modes::TwoTrackResonance::kKstar0Bar)) {
        mHistogramRegistry->fill(HIST(resoPrefix) + HIST(AnalysisDir) + HIST(GetHistName(kSign, HistTable)), resonance.sign());
      }
    }

    if constexpr (modes::isFlagSet(mode, modes::Mode::kQa)) {
      mHistogramRegistry->fill(HIST(resoPrefix) + HIST(QaDir) + HIST(GetHistName(kPtVsEta, HistTable)), resonance.pt(), resonance.eta());
      mHistogramRegistry->fill(HIST(resoPrefix) + HIST(QaDir) + HIST(GetHistName(kPtVsPhi, HistTable)), resonance.pt(), resonance.phi());
      mHistogramRegistry->fill(HIST(resoPrefix) + HIST(QaDir) + HIST(GetHistName(kPhiVsEta, HistTable)), resonance.phi(), resonance.eta());
      mHistogramRegistry->fill(HIST(resoPrefix) + HIST(QaDir) + HIST(GetHistName(kPtVsMass, HistTable)), resonance.pt(), resonance.mass());
    }
  }

 private:
  o2::framework::HistogramRegistry* mHistogramRegistry = nullptr;

  trackhistmanager::TrackHistManager<posDauPrefix, mode> mPosDauManager;
  trackhistmanager::TrackHistManager<negDauPrefix, mode> mNegDauManager;
};
}; // namespace twotrackresonancehistmanager
}; // namespace o2::analysis::femto
#endif // PWGCF_FEMTO_CORE_TWOTRACKRESONANCEHISTMANAGER_H_
