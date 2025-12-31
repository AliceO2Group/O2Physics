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

/// \file kinkHistManager.h
/// \brief histogram manager for kink histograms
/// \author Anton Riedel, TU München, anton.riedel@cern.ch
/// \author Henrik Fribert, TU München, henrik.fribert@cern.ch

#ifndef PWGCF_FEMTO_CORE_KINKHISTMANAGER_H_
#define PWGCF_FEMTO_CORE_KINKHISTMANAGER_H_

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
namespace kinkhistmanager
{
// enum for kink histograms
enum KinkHist {
  // analysis
  kPt,
  kEta,
  kPhi,
  kMass,
  kSign,
  // qa variables
  kKinkAngle,
  kDcaMothToPV,
  kDcaDaugToPV,
  kDecayVtxX,
  kDecayVtxY,
  kDecayVtxZ,
  kDecayVtx,
  kTransRadius,
  // 2d qa
  kPtVsEta,
  kPtVsPhi,
  kPhiVsEta,
  kPtVsKinkAngle,
  kPtVsDecayRadius,
  kKinkHistLast
};

#define KINK_DEFAULT_BINNING(defaultMassMin, defaultMassMax)                                       \
  o2::framework::ConfigurableAxis pt{"pt", {{600, 0, 6}}, "Pt"};                                   \
  o2::framework::ConfigurableAxis eta{"eta", {{300, -1.5, 1.5}}, "Eta"};                           \
  o2::framework::ConfigurableAxis phi{"phi", {{720, 0, 1.f * o2::constants::math::TwoPI}}, "Phi"}; \
  o2::framework::ConfigurableAxis mass{"mass", {{200, defaultMassMin, defaultMassMax}}, "Mass"};   \
  o2::framework::ConfigurableAxis sign{"sign", {{3, -1.5, 1.5}}, "Sign"};

template <const char* Prefix>
struct ConfSigmaBinning : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  KINK_DEFAULT_BINNING(1.1, 1.3)
};
template <const char* Prefix>
struct ConfSigmaPlusBinning : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  KINK_DEFAULT_BINNING(1.1, 1.3)
};
#undef KINK_DEFAULT_BINNING

constexpr const char PrefixSigmaBinning1[] = "SigmaBinning1";
using ConfSigmaBinning1 = ConfSigmaBinning<PrefixSigmaBinning1>;

constexpr const char PrefixSigmaPlusBinning1[] = "SigmaPlusBinning1";
using ConfSigmaPlusBinning1 = ConfSigmaPlusBinning<PrefixSigmaPlusBinning1>;

template <const char* Prefix>
struct ConfKinkQaBinning : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  o2::framework::Configurable<bool> plot2d{"plot2d", true, "Enable 2d QA h histograms"};
  o2::framework::ConfigurableAxis kinkAngle{"kinkAngle", {{100, 0, 3.15}}, "Kink Angle (rad)"};
  o2::framework::ConfigurableAxis dcaMothToPV{"dcaMothToPV", {{150, 0, 1.5}}, "Mother DCA to PV (cm)"};
  o2::framework::ConfigurableAxis dcaDaugToPV{"dcaDaugToPV", {{1000, 0, 100}}, "Daughter DCA to PV (cm)"};
  o2::framework::ConfigurableAxis decayVertex{"decayVertex", {{100, 0, 100}}, "Decay vertex position (cm)"};
  o2::framework::ConfigurableAxis transRadius{"transRadius", {{100, 0, 100}}, "Transverse radius (cm)"};
};

constexpr const char PrefixSigmaQaBinning1[] = "SigmaQaBinning1";
using ConfSigmaQaBinning1 = ConfKinkQaBinning<PrefixSigmaQaBinning1>;

constexpr const char PrefixSigmaPlusQaBinning1[] = "SigmaPlusQaBinning1";
using ConfSigmaPlusQaBinning1 = ConfKinkQaBinning<PrefixSigmaPlusQaBinning1>;

// must be in sync with enum KinkHist
// the enum gives the correct index in the array
constexpr std::array<histmanager::HistInfo<KinkHist>, kKinkHistLast> HistTable = {
  {{kPt, o2::framework::kTH1F, "hPt", "Transverse Momentum; p_{T} (GeV/#it{c}); Entries"},
   {kEta, o2::framework::kTH1F, "hEta", "Pseudorapidity; #eta; Entries"},
   {kPhi, o2::framework::kTH1F, "hPhi", "Azimuthal angle; #varphi; Entries"},
   {kMass, o2::framework::kTH1F, "hMass", "Invariant Mass; m_{Inv} (GeV/#it{c}^{2}); Entries"},
   {kSign, o2::framework::kTH1F, "hSign", "Sign; sign; Entries"},
   {kKinkAngle, o2::framework::kTH1F, "hKinkAngle", "Kink Angle; Angle (rad); Entries"},
   {kDcaMothToPV, o2::framework::kTH1F, "hDcaMothToPV", "Mother DCA to PV; DCA (cm); Entries"},
   {kDcaDaugToPV, o2::framework::kTH1F, "hDcaDaugToPV", "Daughter DCA to PV; DCA (cm); Entries"},
   {kDecayVtxX, o2::framework::kTH1F, "hDecayVtxX", "Decay Vertex X; x (cm); Entries"},
   {kDecayVtxY, o2::framework::kTH1F, "hDecayVtxY", "Decay Vertex Y; y (cm); Entries"},
   {kDecayVtxZ, o2::framework::kTH1F, "hDecayVtxZ", "Decay Vertex Z; z (cm); Entries"},
   {kDecayVtx, o2::framework::kTH1F, "hDecayVtx", "Decay Distance from PV; r (cm); Entries"},
   {kTransRadius, o2::framework::kTH1F, "hTransRadius", "Transverse Decay Radius; r_{xy} (cm); Entries"},
   {kPtVsEta, o2::framework::kTH2F, "hPtVsEta", "p_{T} vs #eta; p_{T} (GeV/#it{c}); #eta"},
   {kPtVsPhi, o2::framework::kTH2F, "hPtVsPhi", "p_{T} vs #varphi; p_{T} (GeV/#it{c}); #varphi"},
   {kPhiVsEta, o2::framework::kTH2F, "hPhiVsEta", "#varphi vs #eta; #varphi; #eta"},
   {kPtVsKinkAngle, o2::framework::kTH2F, "hPtVsKinkAngle", "p_{T} vs kink angle; p_{T} (GeV/#it{c}); kink angle (rad)"},
   {kPtVsDecayRadius, o2::framework::kTH2F, "hPtVsDecayRadius", "p_{T} vs transverse decay radius; p_{T} (GeV/#it{c}); r_{xy} (cm)"}}};

template <typename T>
auto makeKinkHistSpecMap(const T& confBinningAnalysis)
{
  return std::map<KinkHist, std::vector<framework::AxisSpec>>{
    {kPt, {confBinningAnalysis.pt}},
    {kEta, {confBinningAnalysis.eta}},
    {kPhi, {confBinningAnalysis.phi}},
    {kMass, {confBinningAnalysis.mass}},
    {kSign, {confBinningAnalysis.sign}}};
}

template <typename T1, typename T2>
std::map<KinkHist, std::vector<framework::AxisSpec>> makeKinkQaHistSpecMap(T1 const& confBinningAnalysis, T2 const& confBinningQa)
{
  return std::map<KinkHist, std::vector<framework::AxisSpec>>{
    {kPt, {confBinningAnalysis.pt}},
    {kEta, {confBinningAnalysis.eta}},
    {kPhi, {confBinningAnalysis.phi}},
    {kMass, {confBinningAnalysis.mass}},
    {kSign, {confBinningAnalysis.sign}},
    {kKinkAngle, {confBinningQa.kinkAngle}},
    {kDcaMothToPV, {confBinningQa.dcaMothToPV}},
    {kDcaDaugToPV, {confBinningQa.dcaDaugToPV}},
    {kDecayVtxX, {confBinningQa.decayVertex}},
    {kDecayVtxY, {confBinningQa.decayVertex}},
    {kDecayVtxZ, {confBinningQa.decayVertex}},
    {kDecayVtx, {confBinningQa.decayVertex}},
    {kTransRadius, {confBinningQa.transRadius}},
    {kPtVsEta, {confBinningAnalysis.pt, confBinningAnalysis.eta}},
    {kPtVsPhi, {confBinningAnalysis.pt, confBinningAnalysis.phi}},
    {kPhiVsEta, {confBinningAnalysis.phi, confBinningAnalysis.eta}},
    {kPtVsKinkAngle, {confBinningAnalysis.pt, confBinningQa.kinkAngle}},
    {kPtVsDecayRadius, {confBinningAnalysis.pt, confBinningQa.transRadius}}};
}

constexpr char PrefixSigmaQa[] = "SigmaQA/";
constexpr char PrefixSigma1[] = "Sigma1/";
constexpr char PrefixSigma2[] = "Sigma2/";
constexpr char PrefixSigmaPlusQa[] = "SigmaPlusQA/";
constexpr char PrefixSigmaPlus1[] = "SigmaPlus1/";
constexpr char PrefixSigmaPlus2[] = "SigmaPlus2/";

constexpr std::string_view AnalysisDir = "Kinematics/";
constexpr std::string_view QaDir = "QA/";

constexpr int AbsChargeDaughters = 1;

/// \class KinkHistManager
/// \brief Class for histogramming event properties
// template <femtomodes::Mode mode>
template <const char* kinkPrefix,
          const char* chaDauPrefix,
          modes::Kink kink>
class KinkHistManager
{
 public:
  KinkHistManager() = default;
  ~KinkHistManager() = default;

  template <modes::Mode mode>
  void init(o2::framework::HistogramRegistry* registry,
            std::map<KinkHist, std::vector<o2::framework::AxisSpec>> const& KinkSpecs,
            std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> const& ChaDauSpecs)
  {
    mHistogramRegistry = registry;
    mChaDauManager.template init<mode>(registry, ChaDauSpecs, AbsChargeDaughters);
    if constexpr (isFlagSet(mode, modes::Mode::kAnalysis)) {
      initAnalysis(KinkSpecs);
    }
    if constexpr (isFlagSet(mode, modes::Mode::kQa)) {
      initQa(KinkSpecs);
    }
  }

  template <modes::Mode mode, typename T1, typename T2>
  void init(o2::framework::HistogramRegistry* registry,
            std::map<KinkHist, std::vector<o2::framework::AxisSpec>> const& KinkSpecs,
            T1 const& KinkConfBinningQa,
            std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> const& ChaDauSpecs,
            T2 const& ChaDauConfBinningQa)
  {
    mHistogramRegistry = registry;
    mChaDauManager.template init<mode>(registry, ChaDauSpecs, ChaDauConfBinningQa, AbsChargeDaughters);
    this->enableOptionalHistograms(KinkConfBinningQa);
    if constexpr (isFlagSet(mode, modes::Mode::kAnalysis)) {
      initAnalysis(KinkSpecs);
    }
    if constexpr (isFlagSet(mode, modes::Mode::kQa)) {
      initQa(KinkSpecs);
    }
  }

  template <modes::Mode mode, typename T1, typename T2>
  void fill(T1 const& kinkcandidate, T2 const& tracks)
  {
    // this used to work, still under investigation
    // auto chaDaughter = kinkcandidate.template chaDau_as<T2>();
    auto chaDaughter = tracks.rawIteratorAt(kinkcandidate.chaDauId() - tracks.offset());
    mChaDauManager.template fill<mode>(chaDaughter, tracks);
    if constexpr (isFlagSet(mode, modes::Mode::kAnalysis)) {
      fillAnalysis(kinkcandidate);
    }
    if constexpr (isFlagSet(mode, modes::Mode::kQa)) {
      fillQa(kinkcandidate);
    }
  }

 private:
  template <typename T1>
  void enableOptionalHistograms(T1 const& KinkConfBinningQa)
  {
    mPlot2d = KinkConfBinningQa.plot2d.value;
  }

  void initAnalysis(std::map<KinkHist, std::vector<o2::framework::AxisSpec>> const& KinkSpecs)
  {
    std::string analysisDir = std::string(kinkPrefix) + std::string(AnalysisDir);
    mHistogramRegistry->add(analysisDir + getHistNameV2(kPt, HistTable), getHistDesc(kPt, HistTable), getHistType(kPt, HistTable), {KinkSpecs.at(kPt)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kEta, HistTable), getHistDesc(kEta, HistTable), getHistType(kEta, HistTable), {KinkSpecs.at(kEta)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kPhi, HistTable), getHistDesc(kPhi, HistTable), getHistType(kPhi, HistTable), {KinkSpecs.at(kPhi)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kMass, HistTable), getHistDesc(kMass, HistTable), getHistType(kMass, HistTable), {KinkSpecs.at(kMass)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kSign, HistTable), getHistDesc(kSign, HistTable), getHistType(kSign, HistTable), {KinkSpecs.at(kSign)});
  }

  void initQa(std::map<KinkHist, std::vector<o2::framework::AxisSpec>> const& KinkSpecs)
  {
    std::string qaDir = std::string(kinkPrefix) + std::string(QaDir);
    // Kink-specific QA histograms
    mHistogramRegistry->add(qaDir + getHistNameV2(kKinkAngle, HistTable), getHistDesc(kKinkAngle, HistTable), getHistType(kKinkAngle, HistTable), {KinkSpecs.at(kKinkAngle)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kDcaMothToPV, HistTable), getHistDesc(kDcaMothToPV, HistTable), getHistType(kDcaMothToPV, HistTable), {KinkSpecs.at(kDcaMothToPV)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kDcaDaugToPV, HistTable), getHistDesc(kDcaDaugToPV, HistTable), getHistType(kDcaDaugToPV, HistTable), {KinkSpecs.at(kDcaDaugToPV)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kDecayVtxX, HistTable), getHistDesc(kDecayVtxX, HistTable), getHistType(kDecayVtxX, HistTable), {KinkSpecs.at(kDecayVtxX)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kDecayVtxY, HistTable), getHistDesc(kDecayVtxY, HistTable), getHistType(kDecayVtxY, HistTable), {KinkSpecs.at(kDecayVtxY)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kDecayVtxZ, HistTable), getHistDesc(kDecayVtxZ, HistTable), getHistType(kDecayVtxZ, HistTable), {KinkSpecs.at(kDecayVtxZ)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kDecayVtx, HistTable), getHistDesc(kDecayVtx, HistTable), getHistType(kDecayVtx, HistTable), {KinkSpecs.at(kDecayVtx)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kTransRadius, HistTable), getHistDesc(kTransRadius, HistTable), getHistType(kTransRadius, HistTable), {KinkSpecs.at(kTransRadius)});
    if (mPlot2d) {
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsEta, HistTable), getHistDesc(kPtVsEta, HistTable), getHistType(kPtVsEta, HistTable), {KinkSpecs.at(kPtVsEta)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsPhi, HistTable), getHistDesc(kPtVsPhi, HistTable), getHistType(kPtVsPhi, HistTable), {KinkSpecs.at(kPtVsPhi)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPhiVsEta, HistTable), getHistDesc(kPhiVsEta, HistTable), getHistType(kPhiVsEta, HistTable), {KinkSpecs.at(kPhiVsEta)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsKinkAngle, HistTable), getHistDesc(kPtVsKinkAngle, HistTable), getHistType(kPtVsKinkAngle, HistTable), {KinkSpecs.at(kPtVsKinkAngle)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsDecayRadius, HistTable), getHistDesc(kPtVsDecayRadius, HistTable), getHistType(kPtVsDecayRadius, HistTable), {KinkSpecs.at(kPtVsDecayRadius)});
    }
  }

  /// Fill histograms for kink candidates
  /// \param kinkcandidate Kink candidate to fill histograms for
  template <typename T>
  void fillAnalysis(T const& kinkcandidate)
  {
    mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(AnalysisDir) + HIST(getHistName(kPt, HistTable)), kinkcandidate.pt());
    mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(AnalysisDir) + HIST(getHistName(kEta, HistTable)), kinkcandidate.eta());
    mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(AnalysisDir) + HIST(getHistName(kPhi, HistTable)), kinkcandidate.phi());
    mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(AnalysisDir) + HIST(getHistName(kMass, HistTable)), kinkcandidate.mass());

    if constexpr (isEqual(kink, modes::Kink::kSigma) || isEqual(kink, modes::Kink::kSigmaPlus)) {
      mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(AnalysisDir) + HIST(getHistName(kSign, HistTable)), kinkcandidate.sign());
    }
  }

  template <typename T>
  void fillQa(T const& kinkcandidate)
  {
    // Kink-specific QA histograms
    mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(QaDir) + HIST(getHistName(kKinkAngle, HistTable)), kinkcandidate.kinkAngle());
    mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(QaDir) + HIST(getHistName(kDcaMothToPV, HistTable)), kinkcandidate.dcaMothToPV());
    mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(QaDir) + HIST(getHistName(kDcaDaugToPV, HistTable)), kinkcandidate.dcaDaugToPV());
    mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(QaDir) + HIST(getHistName(kDecayVtxX, HistTable)), kinkcandidate.decayVtxX());
    mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(QaDir) + HIST(getHistName(kDecayVtxY, HistTable)), kinkcandidate.decayVtxY());
    mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(QaDir) + HIST(getHistName(kDecayVtxZ, HistTable)), kinkcandidate.decayVtxZ());
    // Calculate decay distance from PV
    float decayDistance = std::hypot(kinkcandidate.decayVtxX(), kinkcandidate.decayVtxY(), kinkcandidate.decayVtxZ());
    mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(QaDir) + HIST(getHistName(kDecayVtx, HistTable)), decayDistance);
    mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(QaDir) + HIST(getHistName(kTransRadius, HistTable)), kinkcandidate.transRadius());
    if (mPlot2d) {
      mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(QaDir) + HIST(getHistName(kPtVsEta, HistTable)), kinkcandidate.pt(), kinkcandidate.eta());
      mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(QaDir) + HIST(getHistName(kPtVsPhi, HistTable)), kinkcandidate.pt(), kinkcandidate.phi());
      mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(QaDir) + HIST(getHistName(kPhiVsEta, HistTable)), kinkcandidate.phi(), kinkcandidate.eta());
      mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(QaDir) + HIST(getHistName(kPtVsKinkAngle, HistTable)), kinkcandidate.pt(), kinkcandidate.kinkAngle());
      mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(QaDir) + HIST(getHistName(kPtVsDecayRadius, HistTable)), kinkcandidate.pt(), kinkcandidate.transRadius());
    }
  }

  o2::framework::HistogramRegistry* mHistogramRegistry = nullptr;
  trackhistmanager::TrackHistManager<chaDauPrefix> mChaDauManager;
  bool mPlot2d = true;
};
}; // namespace kinkhistmanager
}; // namespace o2::analysis::femto
#endif // PWGCF_FEMTO_CORE_KINKHISTMANAGER_H_
