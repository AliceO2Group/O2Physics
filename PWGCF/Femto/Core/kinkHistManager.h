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
  // mc
  kOrigin,
  kPdg,
  kPdgMother,
  kPdgPartonicMother,
  kTruePt,
  kTrueEta,
  kTruePhi,
  // histograms for fraction estimation of kinks
  kNoMcParticle,
  kPrimary,
  kFromWrongCollision,
  kFromMaterial,
  kSecondary1,
  kSecondary2,
  kSecondary3,
  kSecondaryOther,

  kKinkHistLast
};

constexpr std::size_t MaxSecondary = 3;

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

template <const char* Prefix>
struct ConfKinkMcBinning : o2::framework::ConfigurableGroup {
  std::string prefix = std::string(Prefix);
  o2::framework::ConfigurableAxis pdgCodes{"pdgCodes", {{8001, -4000.5, 4000.5}}, "PDG codes of selected V0s"};
  o2::framework::ConfigurableAxis statusCode{"statusCode", {{21, -0.5, 20.5}}, "Status codes (i.e. Origin)"};
  o2::framework::Configurable<bool> plotOrigins{"plotOrigins", true, "Plot pt vs cosPa for different particle origins"};
  o2::framework::Configurable<std::vector<int>> pdgCodesForMothersOfSecondary{"pdgCodesForMothersOfSecondary", {3312, 3334}, "PDG codes of mothers of secondaries (Max 3 will be considered)"};
  o2::framework::ConfigurableAxis pt{"pt", {{150, 0, 3}}, "Pt"};
  o2::framework::ConfigurableAxis kinkAngle{"kinkAngle", {{315, 0, 3.15}}, "Kink angle"};
};

constexpr const char PrefixSigmaMcBinning1[] = "SigmaMcBinning1";
using ConfSigmaMcBinning = ConfKinkMcBinning<PrefixSigmaMcBinning1>;

constexpr const char PrefixSigmaPlusMcBinning1[] = "SigmaPlusMcBinning1";
using ConfSigmaPlusMcBinning = ConfKinkMcBinning<PrefixSigmaPlusMcBinning1>;

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
   {kPtVsDecayRadius, o2::framework::kTH2F, "hPtVsDecayRadius", "p_{T} vs transverse decay radius; p_{T} (GeV/#it{c}); r_{xy} (cm)"},
   {kOrigin, o2::framework::kTH1F, "hOrigin", "Status Codes (=Origin); Status Code; Entries"},
   {kPdg, o2::framework::kTH1F, "hPdg", "PDG Codes of reconstructed kinks; PDG Code; Entries"},
   {kPdgMother, o2::framework::kTH1F, "hPdgMother", "PDG Codes of mother of reconstructed kink; PDG Code; Entries"},
   {kPdgPartonicMother, o2::framework::kTH1F, "hPdgPartonicMother", "PDG Codes of partonic mother reconstructed knik; PDG Code; Entries"},
   {kTruePt, o2::framework::kTH1F, "hTruePt", "True transverse momentum; p_{T} (GeV/#it{c}); Entries"},
   {kTrueEta, o2::framework::kTH1F, "hTrueEta", "True pseudorapdity; #eta; Entries"},
   {kTruePhi, o2::framework::kTH1F, "hTruePhi", "True azimuthal angle; #varphi; Entries"},
   {kNoMcParticle, o2::framework::kTH2F, "hNoMcParticle", "Wrongly reconstructed particles; p_{T} (GeV/#it{c}); cos(#alpha)"},
   {kPrimary, o2::framework::kTH2F, "hPrimary", "Primary particles; p_{T} (GeV/#it{c}); cos(#alpha)"},
   {kFromWrongCollision, o2::framework::kTH2F, "hFromWrongCollision", "Particles associated to wrong collision; p_{T} (GeV/#it{c}); cos(#alpha)"},
   {kFromMaterial, o2::framework::kTH2F, "hFromMaterial", "Particles from material; p_{T} (GeV/#it{c}); cos(#alpha)"},
   {kSecondary1, o2::framework::kTH2F, "hFromSecondary1", "Particles from secondary decay; p_{T} (GeV/#it{c}); cos(#alpha)"},
   {kSecondary2, o2::framework::kTH2F, "hFromSecondary2", "Particles from seconary decay; p_{T} (GeV/#it{c}); cos(#alpha)"},
   {kSecondary3, o2::framework::kTH2F, "hFromSecondary3", "Particles from seconary decay; p_{T} (GeV/#it{c}); cos(#alpha)"},
   {kSecondaryOther, o2::framework::kTH2F, "hFromSecondaryOther", "Particles from every other seconary decay; p_{T} (GeV/#it{c}); cos(#alpha)"}}};

#define KINK_HIST_ANALYSIS_MAP(conf) \
  {kPt, {conf.pt}},                  \
    {kEta, {conf.eta}},              \
    {kPhi, {conf.phi}},              \
    {kMass, {conf.mass}},            \
    {kSign, {conf.sign}},

#define KINK_HIST_QA_MAP(confAnalysis, confQa)             \
  {kKinkAngle, {confQa.kinkAngle}},                        \
    {kDcaMothToPV, {confQa.dcaMothToPV}},                  \
    {kDcaDaugToPV, {confQa.dcaDaugToPV}},                  \
    {kDecayVtxX, {confQa.decayVertex}},                    \
    {kDecayVtxY, {confQa.decayVertex}},                    \
    {kDecayVtxZ, {confQa.decayVertex}},                    \
    {kDecayVtx, {confQa.decayVertex}},                     \
    {kTransRadius, {confQa.transRadius}},                  \
    {kPtVsEta, {confAnalysis.pt, confAnalysis.eta}},       \
    {kPtVsPhi, {confAnalysis.pt, confAnalysis.phi}},       \
    {kPhiVsEta, {confAnalysis.phi, confAnalysis.eta}},     \
    {kPtVsKinkAngle, {confAnalysis.pt, confQa.kinkAngle}}, \
    {kPtVsDecayRadius, {confAnalysis.pt, confQa.transRadius}},

#define KINK_HIST_MC_MAP(confAnalysis, confMc)            \
  {kTruePt, {confAnalysis.pt}},                           \
    {kTrueEta, {confAnalysis.eta}},                       \
    {kTruePhi, {confAnalysis.phi}},                       \
    {kOrigin, {confMc.statusCode}},                       \
    {kPdg, {confMc.pdgCodes}},                            \
    {kPdgMother, {confMc.pdgCodes}},                      \
    {kPdgPartonicMother, {confMc.pdgCodes}},              \
    {kNoMcParticle, {confMc.pt, confMc.kinkAngle}},       \
    {kPrimary, {confMc.pt, confMc.kinkAngle}},            \
    {kFromWrongCollision, {confMc.pt, confMc.kinkAngle}}, \
    {kFromMaterial, {confMc.pt, confMc.kinkAngle}},       \
    {kSecondary1, {confMc.pt, confMc.kinkAngle}},         \
    {kSecondary2, {confMc.pt, confMc.kinkAngle}},         \
    {kSecondary3, {confMc.pt, confMc.kinkAngle}},         \
    {kSecondaryOther, {confMc.pt, confMc.kinkAngle}},

template <typename T>
auto makeKinkHistSpecMap(const T& confBinningAnalysis)
{
  return std::map<KinkHist, std::vector<framework::AxisSpec>>{
    KINK_HIST_ANALYSIS_MAP(confBinningAnalysis)};
}

template <typename T1, typename T2>
std::map<KinkHist, std::vector<framework::AxisSpec>> makeKinkQaHistSpecMap(T1 const& confBinningAnalysis, T2 const& confBinningQa)
{
  return std::map<KinkHist, std::vector<framework::AxisSpec>>{
    KINK_HIST_ANALYSIS_MAP(confBinningAnalysis)
      KINK_HIST_QA_MAP(confBinningAnalysis, confBinningQa)};
}

template <typename T1, typename T2, typename T3>
std::map<KinkHist, std::vector<framework::AxisSpec>> makeKinkMcQaHistSpecMap(T1 const& confBinningAnalysis, T2 const& confBinningQa, T3 const& confBinningMc)
{
  return std::map<KinkHist, std::vector<framework::AxisSpec>>{
    KINK_HIST_ANALYSIS_MAP(confBinningAnalysis)
      KINK_HIST_QA_MAP(confBinningAnalysis, confBinningQa)
        KINK_HIST_MC_MAP(confBinningAnalysis, confBinningMc)};
}

#undef KINK_HIST_ANALYSIS_MAP
#undef KINK_HIST_QA_MAP
#undef KINK_HIST_MC_MAP

constexpr char PrefixSigmaQa[] = "SigmaQA/";
constexpr char PrefixSigma1[] = "Sigma1/";
constexpr char PrefixSigma2[] = "Sigma2/";
constexpr char PrefixSigmaPlusQa[] = "SigmaPlusQA/";
constexpr char PrefixSigmaPlus1[] = "SigmaPlus1/";
constexpr char PrefixSigmaPlus2[] = "SigmaPlus2/";

constexpr std::string_view AnalysisDir = "Kinematics/";
constexpr std::string_view QaDir = "QA/";
constexpr std::string_view McDir = "MC/";

constexpr int AbsChargeDaughters = 1;

template <const char* kinkPrefix,
          const char* chaDauPrefix,
          modes::Kink kink>
class KinkHistManager
{
 public:
  KinkHistManager() = default;
  ~KinkHistManager() = default;

  // init for analysis
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
    if constexpr (isFlagSet(mode, modes::Mode::kMc)) {
      initMc(KinkSpecs);
    }
  }

  // init for qa
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
    if constexpr (isFlagSet(mode, modes::Mode::kMc)) {
      initMc(KinkSpecs);
    }
  }

  // init for qa + mc
  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4>
  void init(o2::framework::HistogramRegistry* registry,
            std::map<KinkHist, std::vector<o2::framework::AxisSpec>> const& KinkSpecs,
            T1 const& KinkConfBinningQa,
            T2 const& KinkConfBinningMc,
            std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> const& ChaDauSpecs,
            T3 const& ChaDauConfBinningQa,
            T4 const& ChaDauConfBinningMc)
  {
    mHistogramRegistry = registry;
    mChaDauManager.template init<mode>(registry, ChaDauSpecs, ChaDauConfBinningQa, ChaDauConfBinningMc, AbsChargeDaughters);
    this->enableOptionalHistograms(KinkConfBinningQa, KinkConfBinningMc);
    if constexpr (isFlagSet(mode, modes::Mode::kAnalysis)) {
      initAnalysis(KinkSpecs);
    }
    if constexpr (isFlagSet(mode, modes::Mode::kQa)) {
      initQa(KinkSpecs);
    }
    if constexpr (isFlagSet(mode, modes::Mode::kMc)) {
      initMc(KinkSpecs);
    }
  }

  template <modes::Mode mode, typename T1, typename T2>
  void fill(T1 const& kinkCandidate, T2 const& tracks)
  {
    // this used to work, still under investigation
    // auto chaDaughter = kinkcandidate.template chaDau_as<T2>();
    auto chaDaughter = tracks.rawIteratorAt(kinkCandidate.chaDauId() - tracks.offset());
    mChaDauManager.template fill<mode>(chaDaughter, tracks);
    if constexpr (isFlagSet(mode, modes::Mode::kAnalysis)) {
      fillAnalysis(kinkCandidate);
    }
    if constexpr (isFlagSet(mode, modes::Mode::kQa)) {
      fillQa(kinkCandidate);
    }
  }

  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5>
  void fill(T1 const& kinkCandidate, T2 const& tracks, T3 const& mcParticles, T4 const& mcMothers, T5 const& mcPartonicMothers)
  {
    auto chaDaughter = tracks.rawIteratorAt(kinkCandidate.chaDauId() - tracks.offset());
    mChaDauManager.template fill<mode>(chaDaughter, tracks, mcParticles, mcMothers, mcPartonicMothers);
    if constexpr (modes::isFlagSet(mode, modes::Mode::kAnalysis)) {
      fillAnalysis(kinkCandidate);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kQa)) {
      fillQa(kinkCandidate);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kMc)) {
      fillMc(kinkCandidate, mcParticles, mcMothers, mcPartonicMothers);
    }
  }

 private:
  template <typename T1>
  void enableOptionalHistograms(T1 const& KinkConfBinningQa)
  {
    mPlot2d = KinkConfBinningQa.plot2d.value;
  }

  // for qa and mc
  template <typename T1, typename T2>
  void enableOptionalHistograms(T1 const& V0ConfBinningQa, T2 const& V0ConfBinningMc)
  {
    this->enableOptionalHistograms(V0ConfBinningQa);

    mPlotOrigins = V0ConfBinningMc.plotOrigins.value;
    mPlotNSecondaries = V0ConfBinningMc.pdgCodesForMothersOfSecondary.value.size();

    for (std::size_t i = 0; i < MaxSecondary; i++) {
      if (i < V0ConfBinningMc.pdgCodesForMothersOfSecondary.value.size()) {
        mPdgCodesSecondaryMother.at(i) = std::abs(V0ConfBinningMc.pdgCodesForMothersOfSecondary.value.at(i));
      } else {
        mPdgCodesSecondaryMother.at(i) = 0;
      }
    }
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

  void initMc(std::map<KinkHist, std::vector<o2::framework::AxisSpec>> const& KinkSpecs)
  {
    std::string mcDir = std::string(kinkPrefix) + std::string(McDir);
    mHistogramRegistry->add(mcDir + getHistNameV2(kTruePt, HistTable), getHistDesc(kTruePt, HistTable), getHistType(kTruePt, HistTable), {KinkSpecs.at(kTruePt)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kTrueEta, HistTable), getHistDesc(kTrueEta, HistTable), getHistType(kTrueEta, HistTable), {KinkSpecs.at(kTrueEta)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kTruePhi, HistTable), getHistDesc(kTruePhi, HistTable), getHistType(kTruePhi, HistTable), {KinkSpecs.at(kTruePhi)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kOrigin, HistTable), getHistDesc(kOrigin, HistTable), getHistType(kOrigin, HistTable), {KinkSpecs.at(kOrigin)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kPdg, HistTable), getHistDesc(kPdg, HistTable), getHistType(kPdg, HistTable), {KinkSpecs.at(kPdg)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kPdgMother, HistTable), getHistDesc(kPdgMother, HistTable), getHistType(kPdgMother, HistTable), {KinkSpecs.at(kPdgMother)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kPdgPartonicMother, HistTable), getHistDesc(kPdgPartonicMother, HistTable), getHistType(kPdgPartonicMother, HistTable), {KinkSpecs.at(kPdgPartonicMother)});

    if (mPlotOrigins) {
      mHistogramRegistry->add(mcDir + getHistNameV2(kNoMcParticle, HistTable), getHistDesc(kNoMcParticle, HistTable), getHistType(kNoMcParticle, HistTable), {KinkSpecs.at(kNoMcParticle)});
      mHistogramRegistry->add(mcDir + getHistNameV2(kPrimary, HistTable), getHistDesc(kPrimary, HistTable), getHistType(kPrimary, HistTable), {KinkSpecs.at(kPrimary)});
      mHistogramRegistry->add(mcDir + getHistNameV2(kFromWrongCollision, HistTable), getHistDesc(kFromWrongCollision, HistTable), getHistType(kFromWrongCollision, HistTable), {KinkSpecs.at(kFromWrongCollision)});
      mHistogramRegistry->add(mcDir + getHistNameV2(kFromMaterial, HistTable), getHistDesc(kFromMaterial, HistTable), getHistType(kFromMaterial, HistTable), {KinkSpecs.at(kFromMaterial)});

      if (mPlotNSecondaries >= histmanager::kSecondaryPlotLevel1) {
        mHistogramRegistry->add(mcDir + getHistNameV2(kSecondary1, HistTable), getHistDesc(kSecondary1, HistTable), getHistType(kSecondary1, HistTable), {KinkSpecs.at(kSecondary1)});
      }
      if (mPlotNSecondaries >= histmanager::kSecondaryPlotLevel2) {
        mHistogramRegistry->add(mcDir + getHistNameV2(kSecondary2, HistTable), getHistDesc(kSecondary2, HistTable), getHistType(kSecondary2, HistTable), {KinkSpecs.at(kSecondary2)});
      }
      if (mPlotNSecondaries >= histmanager::kSecondaryPlotLevel3) {
        mHistogramRegistry->add(mcDir + getHistNameV2(kSecondary3, HistTable), getHistDesc(kSecondary3, HistTable), getHistType(kSecondary3, HistTable), {KinkSpecs.at(kSecondary3)});
      }
      mHistogramRegistry->add(mcDir + getHistNameV2(kSecondaryOther, HistTable), getHistDesc(kSecondaryOther, HistTable), getHistType(kSecondaryOther, HistTable), {KinkSpecs.at(kSecondaryOther)});
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

  template <typename T1, typename T2, typename T3, typename T4>
  void fillMc(T1 const& kinkCandidate, T2 const& /*mcParticles*/, T3 const& /*mcMothers*/, T4 const& /*mcPartonicMothers*/)
  {
    // No MC Particle
    if (!kinkCandidate.has_fMcParticle()) {
      mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(McDir) + HIST(getHistName(kPdg, HistTable)), 0);
      mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(McDir) + HIST(getHistName(kOrigin, HistTable)), static_cast<float>(modes::McOrigin::kNoMcParticle));
      if (mPlotOrigins) {
        mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(McDir) + HIST(getHistName(kNoMcParticle, HistTable)), kinkCandidate.pt(), kinkCandidate.kinkAngle());
      }
      return;
    }

    // Retrieve MC particle
    auto mcParticle = kinkCandidate.template fMcParticle_as<T2>();
    mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(McDir) + HIST(getHistName(kTruePt, HistTable)), mcParticle.pt());
    mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(McDir) + HIST(getHistName(kTrueEta, HistTable)), mcParticle.eta());
    mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(McDir) + HIST(getHistName(kTruePhi, HistTable)), mcParticle.phi());
    mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(McDir) + HIST(getHistName(kOrigin, HistTable)), mcParticle.origin());
    mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(McDir) + HIST(getHistName(kPdg, HistTable)), mcParticle.pdgCode());

    // Mother PDG
    if (kinkCandidate.has_fMcMother()) {
      auto mother = kinkCandidate.template fMcMother_as<T3>();
      mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(McDir) + HIST(getHistName(kPdgMother, HistTable)), mother.pdgCode());
    } else {
      mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(McDir) + HIST(getHistName(kPdgMother, HistTable)), 0);
    }

    // Partonic Mother PDG
    if (kinkCandidate.has_fMcPartMoth()) {
      auto partonicMother = kinkCandidate.template fMcPartMoth_as<T4>();
      mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(McDir) + HIST(getHistName(kPdgPartonicMother, HistTable)), partonicMother.pdgCode());
    } else {
      mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(McDir) + HIST(getHistName(kPdgPartonicMother, HistTable)), 0);
    }

    // Plot origins
    if (mPlotOrigins) {
      switch (static_cast<modes::McOrigin>(mcParticle.origin())) {
        case modes::McOrigin::kPhysicalPrimary:
          mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(McDir) + HIST(getHistName(kPrimary, HistTable)), kinkCandidate.pt(), kinkCandidate.kinkAngle());
          break;
        case modes::McOrigin::kFromWrongCollision:
          mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(McDir) + HIST(getHistName(kFromWrongCollision, HistTable)), kinkCandidate.pt(), kinkCandidate.kinkAngle());
          break;
        case modes::McOrigin::kFromMaterial:
          mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(McDir) + HIST(getHistName(kFromMaterial, HistTable)), kinkCandidate.pt(), kinkCandidate.kinkAngle());
          break;
        case modes::McOrigin::kFromSecondaryDecay:
          if (kinkCandidate.has_fMcMother()) {
            auto mother = kinkCandidate.template fMcMother_as<T3>();
            int motherPdgCode = std::abs(mother.pdgCode());
            // Switch on PDG of the mother
            if (mPlotNSecondaries >= histmanager::kSecondaryPlotLevel1 && motherPdgCode == mPdgCodesSecondaryMother[0]) {
              mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(McDir) + HIST(getHistName(kSecondary1, HistTable)), kinkCandidate.pt(), kinkCandidate.kinkAngle());
            } else if (mPlotNSecondaries >= histmanager::kSecondaryPlotLevel2 && motherPdgCode == mPdgCodesSecondaryMother[1]) {
              mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(McDir) + HIST(getHistName(kSecondary2, HistTable)), kinkCandidate.pt(), kinkCandidate.kinkAngle());
            } else if (mPlotNSecondaries >= histmanager::kSecondaryPlotLevel3 && motherPdgCode == mPdgCodesSecondaryMother[2]) {
              mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(McDir) + HIST(getHistName(kSecondary3, HistTable)), kinkCandidate.pt(), kinkCandidate.kinkAngle());
            } else {
              mHistogramRegistry->fill(HIST(kinkPrefix) + HIST(McDir) + HIST(getHistName(kSecondaryOther, HistTable)), kinkCandidate.pt(), kinkCandidate.kinkAngle());
            }
          }
          break;
        default:
          // Unknown origin → safely ignore
          break;
      }
    }
  }

  o2::framework::HistogramRegistry* mHistogramRegistry = nullptr;
  trackhistmanager::TrackHistManager<chaDauPrefix> mChaDauManager;
  bool mPlot2d = true;
  bool mPlotOrigins = false;
  int mPlotNSecondaries = 0;
  std::array<int, MaxSecondary> mPdgCodesSecondaryMother = {0};
};
}; // namespace kinkhistmanager
}; // namespace o2::analysis::femto
#endif // PWGCF_FEMTO_CORE_KINKHISTMANAGER_H_
