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
#undef KINK_DEFAULT_BINNING

constexpr const char PrefixSigmaBinning1[] = "SigmaBinning1";
using ConfSigmaBinning1 = ConfSigmaBinning<PrefixSigmaBinning1>;

template <const char* Prefix>
struct ConfKinkQaBinning : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  o2::framework::ConfigurableAxis kinkAngle{"kinkAngle", {{100, 0, 3.15}}, "Kink Angle (rad)"};
  o2::framework::ConfigurableAxis dcaMothToPV{"dcaMothToPV", {{150, 0, 1.5}}, "Mother DCA to PV (cm)"};
  o2::framework::ConfigurableAxis dcaDaugToPV{"dcaDaugToPV", {{1000, 0, 100}}, "Daughter DCA to PV (cm)"};
  o2::framework::ConfigurableAxis decayVertex{"decayVertex", {{100, 0, 100}}, "Decay vertex position (cm)"};
  o2::framework::ConfigurableAxis transRadius{"transRadius", {{100, 0, 100}}, "Transverse radius (cm)"};
};

constexpr const char PrefixSigmaQaBinning1[] = "SigmaQaBinning1";
using ConfSigmaQaBinning1 = ConfKinkQaBinning<PrefixSigmaQaBinning1>;

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

constexpr std::string_view AnalysisDir = "Kinematics/";
constexpr std::string_view QaDir = "QA/";

/// \class KinkHistManager
/// \brief Class for histogramming event properties
// template <femtomodes::Mode mode>
template <const char* kinkPrefix,
          const char* chaDauPrefix,
          modes::Mode mode,
          modes::Kink kink>
class KinkHistManager
{
 public:
  /// Destructor
  virtual ~KinkHistManager() = default;

  /// Initializes histograms for the task
  /// \param registry Histogram registry to be passed
  ///
  void init(o2::framework::HistogramRegistry* registry,
            std::map<KinkHist, std::vector<o2::framework::AxisSpec>> KinkSpecs,
            std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> ChaDauSpecs)
  {
    mHistogramRegistry = registry;
    mChaDauManager.init(registry, ChaDauSpecs);

    if constexpr (isFlagSet(mode, modes::Mode::kAnalysis)) {
      std::string analysisDir = std::string(kinkPrefix) + std::string(AnalysisDir);
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kPt, HistTable), GetHistDesc(kPt, HistTable), GetHistType(kPt, HistTable), {KinkSpecs[kPt]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kEta, HistTable), GetHistDesc(kEta, HistTable), GetHistType(kEta, HistTable), {KinkSpecs[kEta]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kPhi, HistTable), GetHistDesc(kPhi, HistTable), GetHistType(kPhi, HistTable), {KinkSpecs[kPhi]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kMass, HistTable), GetHistDesc(kMass, HistTable), GetHistType(kMass, HistTable), {KinkSpecs[kMass]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kSign, HistTable), GetHistDesc(kSign, HistTable), GetHistType(kSign, HistTable), {KinkSpecs[kSign]});
    }

    if constexpr (isFlagSet(mode, modes::Mode::kQa)) {
      std::string qaDir = std::string(kinkPrefix) + std::string(QaDir);

      // Basic kinematic histograms
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPt, HistTable), GetHistDesc(kPt, HistTable), GetHistType(kPt, HistTable), {KinkSpecs[kPt]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kEta, HistTable), GetHistDesc(kEta, HistTable), GetHistType(kEta, HistTable), {KinkSpecs[kEta]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPhi, HistTable), GetHistDesc(kPhi, HistTable), GetHistType(kPhi, HistTable), {KinkSpecs[kPhi]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kMass, HistTable), GetHistDesc(kMass, HistTable), GetHistType(kMass, HistTable), {KinkSpecs[kMass]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kSign, HistTable), GetHistDesc(kSign, HistTable), GetHistType(kSign, HistTable), {KinkSpecs[kSign]});

      // Kink-specific QA histograms
      mHistogramRegistry->add(qaDir + GetHistNamev2(kKinkAngle, HistTable), GetHistDesc(kKinkAngle, HistTable), GetHistType(kKinkAngle, HistTable), {KinkSpecs[kKinkAngle]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kDcaMothToPV, HistTable), GetHistDesc(kDcaMothToPV, HistTable), GetHistType(kDcaMothToPV, HistTable), {KinkSpecs[kDcaMothToPV]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kDcaDaugToPV, HistTable), GetHistDesc(kDcaDaugToPV, HistTable), GetHistType(kDcaDaugToPV, HistTable), {KinkSpecs[kDcaDaugToPV]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kDecayVtxX, HistTable), GetHistDesc(kDecayVtxX, HistTable), GetHistType(kDecayVtxX, HistTable), {KinkSpecs[kDecayVtxX]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kDecayVtxY, HistTable), GetHistDesc(kDecayVtxY, HistTable), GetHistType(kDecayVtxY, HistTable), {KinkSpecs[kDecayVtxY]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kDecayVtxZ, HistTable), GetHistDesc(kDecayVtxZ, HistTable), GetHistType(kDecayVtxZ, HistTable), {KinkSpecs[kDecayVtxZ]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kDecayVtx, HistTable), GetHistDesc(kDecayVtx, HistTable), GetHistType(kDecayVtx, HistTable), {KinkSpecs[kDecayVtx]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kTransRadius, HistTable), GetHistDesc(kTransRadius, HistTable), GetHistType(kTransRadius, HistTable), {KinkSpecs[kTransRadius]});

      // 2D QA histograms
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPtVsEta, HistTable), GetHistDesc(kPtVsEta, HistTable), GetHistType(kPtVsEta, HistTable), {KinkSpecs[kPtVsEta]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPtVsPhi, HistTable), GetHistDesc(kPtVsPhi, HistTable), GetHistType(kPtVsPhi, HistTable), {KinkSpecs[kPtVsPhi]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPhiVsEta, HistTable), GetHistDesc(kPhiVsEta, HistTable), GetHistType(kPhiVsEta, HistTable), {KinkSpecs[kPhiVsEta]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPtVsKinkAngle, HistTable), GetHistDesc(kPtVsKinkAngle, HistTable), GetHistType(kPtVsKinkAngle, HistTable), {KinkSpecs[kPtVsKinkAngle]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPtVsDecayRadius, HistTable), GetHistDesc(kPtVsDecayRadius, HistTable), GetHistType(kPtVsDecayRadius, HistTable), {KinkSpecs[kPtVsDecayRadius]});
    }
  }

  /// Fill histograms for kink candidates
  /// \param kinkcandidate Kink candidate to fill histograms for
  template <typename T>
  void fill(T const& kinkcandidate)
  {
    if constexpr (isFlagSet(mode, modes::Mode::kAnalysis)) {
      mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(AnalysisDir) + HIST(GetHistName(kPt, HistTable)), kinkcandidate.pt());
      mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(AnalysisDir) + HIST(GetHistName(kEta, HistTable)), kinkcandidate.eta());
      mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(AnalysisDir) + HIST(GetHistName(kPhi, HistTable)), kinkcandidate.phi());
      mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(AnalysisDir) + HIST(GetHistName(kMass, HistTable)), kinkcandidate.mass());

      if constexpr (isEqual(kink, modes::Kink::kSigma)) {
        mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(AnalysisDir) + HIST(GetHistName(kSign, HistTable)), kinkcandidate.sign());
      }
    }

    if constexpr (isFlagSet(mode, modes::Mode::kQa)) {
      // Basic kinematic histograms
      mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(QaDir) + HIST(GetHistName(kPt, HistTable)), kinkcandidate.pt());
      mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(QaDir) + HIST(GetHistName(kEta, HistTable)), kinkcandidate.eta());
      mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(QaDir) + HIST(GetHistName(kPhi, HistTable)), kinkcandidate.phi());
      mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(QaDir) + HIST(GetHistName(kMass, HistTable)), kinkcandidate.mass());

      if constexpr (isEqual(kink, modes::Kink::kSigma)) {
        mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(QaDir) + HIST(GetHistName(kSign, HistTable)), kinkcandidate.sign());
      }

      // Kink-specific QA histograms
      mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(QaDir) + HIST(GetHistName(kKinkAngle, HistTable)), kinkcandidate.kinkAngle());
      mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(QaDir) + HIST(GetHistName(kDcaMothToPV, HistTable)), kinkcandidate.dcaMothToPV());
      mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(QaDir) + HIST(GetHistName(kDcaDaugToPV, HistTable)), kinkcandidate.dcaDaugToPV());
      mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(QaDir) + HIST(GetHistName(kDecayVtxX, HistTable)), kinkcandidate.decayVtxX());
      mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(QaDir) + HIST(GetHistName(kDecayVtxY, HistTable)), kinkcandidate.decayVtxY());
      mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(QaDir) + HIST(GetHistName(kDecayVtxZ, HistTable)), kinkcandidate.decayVtxZ());

      // Calculate decay distance from PV
      float decayDistance = std::sqrt(kinkcandidate.decayVtxX() * kinkcandidate.decayVtxX() +
                                      kinkcandidate.decayVtxY() * kinkcandidate.decayVtxY() +
                                      kinkcandidate.decayVtxZ() * kinkcandidate.decayVtxZ());
      mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(QaDir) + HIST(GetHistName(kDecayVtx, HistTable)), decayDistance);
      mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(QaDir) + HIST(GetHistName(kTransRadius, HistTable)), kinkcandidate.transRadius());

      // 2D QA histograms
      mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(QaDir) + HIST(GetHistName(kPtVsEta, HistTable)), kinkcandidate.pt(), kinkcandidate.eta());
      mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(QaDir) + HIST(GetHistName(kPtVsPhi, HistTable)), kinkcandidate.pt(), kinkcandidate.phi());
      mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(QaDir) + HIST(GetHistName(kPhiVsEta, HistTable)), kinkcandidate.phi(), kinkcandidate.eta());
      mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(QaDir) + HIST(GetHistName(kPtVsKinkAngle, HistTable)), kinkcandidate.pt(), kinkcandidate.kinkAngle());
      mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(QaDir) + HIST(GetHistName(kPtVsDecayRadius, HistTable)), kinkcandidate.pt(), kinkcandidate.transRadius());
    }
  }

  /// Fill histograms for kink candidates - overload with track table argument
  /// \param kinkcandidate Kink candidate to fill histograms for
  /// \param tracks Track table for daughter access
  template <typename T1, typename T2>
  void fill(T1 const& kinkcandidate, T2 const& /*tracks*/)
  {
    auto chaDaughter = kinkcandidate.template chaDau_as<T2>();
    mChaDauManager.fill(chaDaughter);
    fill(kinkcandidate);
  }

 private:
  o2::framework::HistogramRegistry* mHistogramRegistry;
  trackhistmanager::TrackHistManager<chaDauPrefix, mode> mChaDauManager;
};
}; // namespace kinkhistmanager
}; // namespace o2::analysis::femto
#endif // PWGCF_FEMTO_CORE_KINKHISTMANAGER_H_
