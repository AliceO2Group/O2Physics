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

/// \file v0HistManager.h
/// \brief histogram manager for vzero histograms
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#ifndef PWGCF_FEMTO_CORE_V0HISTMANAGER_H_
#define PWGCF_FEMTO_CORE_V0HISTMANAGER_H_

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
namespace v0histmanager
{
// enum for track histograms
enum V0Hist {
  // analysis
  kPt,
  kEta,
  kPhi,
  kMass,
  kSign,
  // qa variables
  kMassLambda,
  kMassAntiLambda,
  kMassK0short,
  kCosPa,
  kDecayDauDca,
  kDecayVtxX,
  kDecayVtxY,
  kDecayVtxZ,
  kDecayVtx,
  kTransRadius,
  // 2d qa
  kPtVsEta,
  kPtVsPhi,
  kPhiVsEta,
  kPtVsCosPa,
  kPtVsLambdaMass,
  kPtVsAntiLambdaMass,
  kPtVsK0shortMass,
  kLambdaMassVsAntiLambdaMass,
  kK0shortMassVsLambdaMass,
  kK0shortMassVsAntiLambdaMass,
  kV0HistLast
};

#define V0_DEFAULT_BINNING(defaultMassMin, defaultMassMax)                                         \
  o2::framework::ConfigurableAxis pt{"pt", {{600, 0, 6}}, "Pt"};                                   \
  o2::framework::ConfigurableAxis eta{"eta", {{300, -1.5, 1.5}}, "Eta"};                           \
  o2::framework::ConfigurableAxis phi{"phi", {{720, 0, 1.f * o2::constants::math::TwoPI}}, "Phi"}; \
  o2::framework::ConfigurableAxis mass{"mass", {{200, defaultMassMin, defaultMassMax}}, "Mass"};   \
  o2::framework::ConfigurableAxis sign{"sign", {{3, -1.5, 1.5}}, "Sign"};

template <const char* Prefix>
struct ConfLambdaBinning : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  V0_DEFAULT_BINNING(1.0, 1.2)
};
template <const char* Prefix>
struct ConfK0shortBinning : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  V0_DEFAULT_BINNING(0.475, 0.515)
};
#undef V0_DEFAULT_BINNING

constexpr const char PrefixLambdaBinning1[] = "LambdaBinning1";
using ConfLambdaBinning1 = ConfLambdaBinning<PrefixLambdaBinning1>;
constexpr const char PrefixK0shortBinning1[] = "K0shortBinning1";
using ConfK0shortBinning1 = ConfK0shortBinning<PrefixK0shortBinning1>;

template <const char* Prefix>
struct ConfV0QaBinning : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  o2::framework::ConfigurableAxis cosPa{"cosPa", {{100, 0.9, 1}}, "Cosine of poiting angle"};
  o2::framework::ConfigurableAxis dauDcaAtDecay{"dauDcaAtDecay", {{150, 0, 1.5}}, "Daughter DCA at decay vertex"};
  o2::framework::ConfigurableAxis decayVertex{"decayVertex", {{100, 0, 100}}, "Decay vertex"};
  o2::framework::ConfigurableAxis transRadius{"transRadius", {{100, 0, 100}}, "Transverse radius"};
  o2::framework::ConfigurableAxis massLambda{"massLambda", {{200, 1, 1.2}}, "mass for antiparticle hypothesis"};
  o2::framework::ConfigurableAxis massAntiLambda{"massAntiLambda", {{100, 1, 1.2}}, "mass for antiparticle hypothesis"};
  o2::framework::ConfigurableAxis massK0short{"massK0short", {{200, 0.45, 0.55}}, "Mass for k0short hypothesis"};
};

constexpr const char PrefixLambdaQaBinning1[] = "LambdaQaBinning1";
using ConfLambdaQaBinning1 = ConfV0QaBinning<PrefixLambdaQaBinning1>;

constexpr const char PrefixK0shortQaBinning1[] = "K0shortQaBinning1";
using ConfK0shortQaBinning1 = ConfV0QaBinning<PrefixK0shortQaBinning1>;

// must be in sync with enum TrackVariables
// the enum gives the correct index in the array
constexpr std::array<histmanager::HistInfo<V0Hist>, kV0HistLast> HistTable = {
  {{kPt, o2::framework::kTH1F, "hPt", "Transverse Momentum; p_{T} (GeV/#it{c}); Entries"},
   {kEta, o2::framework::kTH1F, "hEta", "Pseudorapdity; #eta; Entries"},
   {kPhi, o2::framework::kTH1F, "hPhi", "Azimuthal angle; #varphi; Entries"},
   {kMass, o2::framework::kTH1F, "hMass", "Invariant Mass; m_{Inv} (GeV/#it{c}^{2}); Entries"},
   {kSign, o2::framework::kTH1F, "hSign", "Sign (-1 -> antiparticle, 0 -> self conjugate, +1 -> particle); sign; Entries"},
   {kMassLambda, o2::framework::kTH1F, "hMassLambda", "#Lambda mass; m_{p#pi^{-}} (GeV/#it{c}^{2}); Entries"},
   {kMassAntiLambda, o2::framework::kTH1F, "hMassAntiLambda", "#bar{#Lambda} mass; m_{#bar{p}#pi^{+}} (GeV/#it{c}^{2}); Entries"},
   {kMassK0short, o2::framework::kTH1F, "hMassK0short", "K^{0}_{s} mass; m_{#pi^{+}#pi^{-}} (GeV/#it{c}^{2}); Entries"},
   {kCosPa, o2::framework::kTH1F, "hCosPa", "Cosine of pointing angle; coa(#alpha); Entries"},
   {kDecayDauDca, o2::framework::kTH1F, "hDauDca", "Daughter DCA at decay vertex ; DCA_{Decay vertex} (cm); Entries"},
   {kDecayVtxX, o2::framework::kTH1F, "hDecayVtxX", "X coordinate of decay vertex ; DV_{X} (cm); Entries"},
   {kDecayVtxY, o2::framework::kTH1F, "hDecayVtxY", "Y coordinate of decay vertex ; DV_{Y} (cm); Entries"},
   {kDecayVtxZ, o2::framework::kTH1F, "hDecayVtxZ", "Z coordinate of decay vertex ; DV_{Z} (cm); Entries"},
   {kDecayVtx, o2::framework::kTH1F, "hDecayVtx", "Distance of decay vertex from primary vertex ; DV (cm); Entries"},
   {kTransRadius, o2::framework::kTH1F, "hTransRadius", "Transverse radius ; r_{xy} (cm); Entries"},
   {kPtVsEta, o2::framework::kTH2F, "hPtVsEta", "p_{T} vs #eta; p_{T} (GeV/#it{c}) ; #eta"},
   {kPtVsPhi, o2::framework::kTH2F, "hPtVsPhi", "p_{T} vs #varphi; p_{T} (GeV/#it{c}) ; #varphi"},
   {kPhiVsEta, o2::framework::kTH2F, "hPhiVsEta", "#varphi vs #eta; #varphi ; #eta"},
   {kPtVsCosPa, o2::framework::kTH2F, "hPtVsCosPa", "Cosine of poiting angle vs p_{T}; cos(#alpha); p_{T} (GeV/#it{c})"},
   {kPtVsLambdaMass, o2::framework::kTH2F, "hPtVsLambdaMass", "p_{T} vs #Lambda mass; p_{T} (GeV/#it{c}); m_{p#pi^{-}} (GeV/#it{c}^{2})"},
   {kPtVsAntiLambdaMass, o2::framework::kTH2F, "hPtVsAntiLambdaMass", "p_{T} vs #bar{#Lambda} mass; p_{T} (GeV/#it{c}); m_{#bar{p}#pi^{+}} (GeV/#it{c}^{2})"},
   {kPtVsK0shortMass, o2::framework::kTH2F, "hPtVsK0shortMass", "p_{T} vs K^{0}_{S} mass; p_{T} (GeV/#it{c}); m_{#pi^{+}#pi^{-}} (GeV/#it{c}^{2})"},
   {kK0shortMassVsLambdaMass, o2::framework::kTH2F, "hK0shortMassVsLambdaMass", " K^{0}_{S} mass vs #Lambda mass; m_{#pi^{+}#pi^{-}} (GeV/#it{c}^{2}); m_{p#pi^{-}} (GeV/#it{c}^{2})"},
   {kK0shortMassVsAntiLambdaMass, o2::framework::kTH2F, "hK0shortMassVsAntiLambdaMass", "K^{0}_{S} mass vs #bar{#Lambda} mass; m_{#pi^{+}#pi^{-}} (GeV/#it{c}^{2}); m_{#bar{p}#pi^{+}} (GeV/#it{c}^{2})"},
   {kLambdaMassVsAntiLambdaMass, o2::framework::kTH2F, "hLambdaMassVsAntiLambdaMass", "#Lambda mass vs #bar{#Lambda}; m_{p#pi^{-}} (GeV/#it{c}^{2}); m_{#bar{p}#pi^{+}} (GeV/#it{c}^{2})"}}};

template <typename T>
auto makeV0HistSpecMap(const T& confBinningAnalysis)
{
  return std::map<V0Hist, std::vector<framework::AxisSpec>>{
    {kPt, {confBinningAnalysis.pt}},
    {kEta, {confBinningAnalysis.eta}},
    {kPhi, {confBinningAnalysis.phi}},
    {kMass, {confBinningAnalysis.mass}},
    {kSign, {confBinningAnalysis.sign}}};
}

template <typename T1, typename T2>
std::map<V0Hist, std::vector<framework::AxisSpec>> makeV0QaHistSpecMap(T1 const& confBinningAnalysis, T2 const& confBinningQa)
{
  return std::map<V0Hist, std::vector<framework::AxisSpec>>{
    {kPt, {confBinningAnalysis.pt}},
    {kEta, {confBinningAnalysis.eta}},
    {kPhi, {confBinningAnalysis.phi}},
    {kMass, {confBinningAnalysis.mass}},
    {kSign, {confBinningAnalysis.sign}},
    {kCosPa, {confBinningQa.cosPa}},
    {kDecayDauDca, {confBinningQa.dauDcaAtDecay}},
    {kDecayVtxX, {confBinningQa.decayVertex}},
    {kDecayVtxY, {confBinningQa.decayVertex}},
    {kDecayVtxZ, {confBinningQa.decayVertex}},
    {kDecayVtx, {confBinningQa.decayVertex}},
    {kTransRadius, {confBinningQa.transRadius}},
    {kPtVsEta, {confBinningAnalysis.pt, confBinningAnalysis.eta}},
    {kPtVsPhi, {confBinningAnalysis.pt, confBinningAnalysis.phi}},
    {kPhiVsEta, {confBinningAnalysis.phi, confBinningAnalysis.eta}},
    {kPtVsCosPa, {confBinningAnalysis.pt, confBinningQa.cosPa}},
    {kMassLambda, {confBinningQa.massLambda}},
    {kMassAntiLambda, {confBinningQa.massAntiLambda}},
    {kMassK0short, {confBinningQa.massK0short}},
    {kPtVsLambdaMass, {confBinningAnalysis.pt, confBinningQa.massLambda}},
    {kPtVsAntiLambdaMass, {confBinningAnalysis.pt, confBinningQa.massAntiLambda}},
    {kPtVsK0shortMass, {confBinningAnalysis.pt, confBinningQa.massK0short}},
    {kLambdaMassVsAntiLambdaMass, {confBinningQa.massLambda, confBinningQa.massAntiLambda}},
    {kK0shortMassVsLambdaMass, {confBinningQa.massK0short, confBinningQa.massLambda}},
    {kK0shortMassVsAntiLambdaMass, {confBinningQa.massK0short, confBinningQa.massAntiLambda}}};
};

constexpr char PrefixLambdaQa[] = "LambdaQA/";
constexpr char PrefixLambda1[] = "Lambda1/";
constexpr char PrefixLambda2[] = "Lambda2/";
constexpr char PrefixK0shortQa[] = "K0shortQa/";
constexpr char PrefixK0short1[] = "K0short1/";
constexpr char PrefixK0short2[] = "K0short2/";

constexpr char PrefixLambdaCascade[] = "LambdaCascadeQa/";

constexpr std::string_view AnalysisDir = "Kinematics/";
constexpr std::string_view QaDir = "QA/";

/// \class FemtoDreamEventHisto
/// \brief Class for histogramming event properties
// template <femtomodes::Mode mode>
template <const char* v0Prefix,
          const char* posDauPrefix,
          const char* negDauPrefix,
          modes::Mode mode,
          modes::V0 v0>
class V0HistManager
{
 public:
  /// Destructor
  virtual ~V0HistManager() = default;
  /// Initializes histograms for the task
  /// \param registry Histogram registry to be passed
  ///
  void init(o2::framework::HistogramRegistry* registry,
            std::map<V0Hist, std::vector<o2::framework::AxisSpec>> V0Specs,
            std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> PosDauSpecs,
            std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> NegDauSpecs)
  {
    mHistogramRegistry = registry;
    mPosDauManager.init(registry, PosDauSpecs);
    mNegDauManager.init(registry, NegDauSpecs);

    if constexpr (modes::isFlagSet(mode, modes::Mode::kAnalysis)) {
      std::string analysisDir = std::string(v0Prefix) + std::string(AnalysisDir);
      mHistogramRegistry->add(analysisDir + getHistNameV2(kPt, HistTable), getHistDesc(kPt, HistTable), getHistType(kPt, HistTable), {V0Specs[kPt]});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kEta, HistTable), getHistDesc(kEta, HistTable), getHistType(kEta, HistTable), {V0Specs[kEta]});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kPhi, HistTable), getHistDesc(kPhi, HistTable), getHistType(kPhi, HistTable), {V0Specs[kPhi]});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kMass, HistTable), getHistDesc(kMass, HistTable), getHistType(kMass, HistTable), {V0Specs[kMass]});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kSign, HistTable), getHistDesc(kSign, HistTable), getHistType(kSign, HistTable), {V0Specs[kSign]});
    }

    if constexpr (modes::isFlagSet(mode, modes::Mode::kQa)) {
      std::string qaDir = std::string(v0Prefix) + std::string(QaDir);

      mHistogramRegistry->add(qaDir + getHistNameV2(kCosPa, HistTable), getHistDesc(kCosPa, HistTable), getHistType(kCosPa, HistTable), {V0Specs[kCosPa]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kDecayDauDca, HistTable), getHistDesc(kDecayDauDca, HistTable), getHistType(kDecayDauDca, HistTable), {V0Specs[kDecayDauDca]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kDecayVtxX, HistTable), getHistDesc(kDecayVtxX, HistTable), getHistType(kDecayVtxX, HistTable), {V0Specs[kDecayVtxX]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kDecayVtxY, HistTable), getHistDesc(kDecayVtxY, HistTable), getHistType(kDecayVtxY, HistTable), {V0Specs[kDecayVtxY]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kDecayVtxZ, HistTable), getHistDesc(kDecayVtxZ, HistTable), getHistType(kDecayVtxZ, HistTable), {V0Specs[kDecayVtxZ]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kDecayVtx, HistTable), getHistDesc(kDecayVtx, HistTable), getHistType(kDecayVtx, HistTable), {V0Specs[kDecayVtx]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kTransRadius, HistTable), getHistDesc(kTransRadius, HistTable), getHistType(kTransRadius, HistTable), {V0Specs[kTransRadius]});

      // qa 2d
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsEta, HistTable), getHistDesc(kPtVsEta, HistTable), getHistType(kPtVsEta, HistTable), {V0Specs[kPtVsEta]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsPhi, HistTable), getHistDesc(kPtVsPhi, HistTable), getHistType(kPtVsPhi, HistTable), {V0Specs[kPtVsPhi]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPhiVsEta, HistTable), getHistDesc(kPhiVsEta, HistTable), getHistType(kPhiVsEta, HistTable), {V0Specs[kPhiVsEta]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsCosPa, HistTable), getHistDesc(kPtVsCosPa, HistTable), getHistType(kPtVsCosPa, HistTable), {V0Specs[kPtVsCosPa]});

      mHistogramRegistry->add(qaDir + getHistNameV2(kMassLambda, HistTable), getHistDesc(kMassLambda, HistTable), getHistType(kMassLambda, HistTable), {V0Specs[kMassLambda]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kMassAntiLambda, HistTable), getHistDesc(kMassAntiLambda, HistTable), getHistType(kMassAntiLambda, HistTable), {V0Specs[kMassAntiLambda]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kMassK0short, HistTable), getHistDesc(kMassK0short, HistTable), getHistType(kMassK0short, HistTable), {V0Specs[kMassK0short]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsLambdaMass, HistTable), getHistDesc(kPtVsLambdaMass, HistTable), getHistType(kPtVsLambdaMass, HistTable), {V0Specs[kPtVsLambdaMass]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsAntiLambdaMass, HistTable), getHistDesc(kPtVsAntiLambdaMass, HistTable), getHistType(kPtVsAntiLambdaMass, HistTable), {V0Specs[kPtVsAntiLambdaMass]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsK0shortMass, HistTable), getHistDesc(kPtVsK0shortMass, HistTable), getHistType(kPtVsK0shortMass, HistTable), {V0Specs[kPtVsK0shortMass]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kLambdaMassVsAntiLambdaMass, HistTable), getHistDesc(kLambdaMassVsAntiLambdaMass, HistTable), getHistType(kLambdaMassVsAntiLambdaMass, HistTable), {V0Specs[kLambdaMassVsAntiLambdaMass]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kK0shortMassVsLambdaMass, HistTable), getHistDesc(kK0shortMassVsLambdaMass, HistTable), getHistType(kK0shortMassVsLambdaMass, HistTable), {V0Specs[kK0shortMassVsLambdaMass]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kK0shortMassVsAntiLambdaMass, HistTable), getHistDesc(kK0shortMassVsAntiLambdaMass, HistTable), getHistType(kK0shortMassVsAntiLambdaMass, HistTable), {V0Specs[kK0shortMassVsAntiLambdaMass]});
    }
  }

  template <typename T1, typename T2>
  void fill(T1 const& v0candidate, T2 const& tracks)
  {
    auto posDaughter = v0candidate.template posDau_as<T2>();
    mPosDauManager.fill(posDaughter, tracks);
    auto negDaughter = v0candidate.template negDau_as<T2>();
    mNegDauManager.fill(negDaughter, tracks);

    if constexpr (modes::isFlagSet(mode, modes::Mode::kAnalysis)) {
      mHistogramRegistry->fill(HIST(v0Prefix) + HIST(AnalysisDir) + HIST(getHistName(kPt, HistTable)), v0candidate.pt());
      mHistogramRegistry->fill(HIST(v0Prefix) + HIST(AnalysisDir) + HIST(getHistName(kEta, HistTable)), v0candidate.eta());
      mHistogramRegistry->fill(HIST(v0Prefix) + HIST(AnalysisDir) + HIST(getHistName(kPhi, HistTable)), v0candidate.phi());
      mHistogramRegistry->fill(HIST(v0Prefix) + HIST(AnalysisDir) + HIST(getHistName(kMass, HistTable)), v0candidate.mass());

      if constexpr (modes::isEqual(v0, modes::V0::kLambda) || modes::isEqual(v0, modes::V0::kAntiLambda)) {
        mHistogramRegistry->fill(HIST(v0Prefix) + HIST(AnalysisDir) + HIST(getHistName(kSign, HistTable)), v0candidate.sign());
      }
      if constexpr (modes::isEqual(v0, modes::V0::kK0short)) {
        mHistogramRegistry->fill(HIST(v0Prefix) + HIST(AnalysisDir) + HIST(getHistName(kSign, HistTable)), 0);
      }
    }

    if constexpr (modes::isFlagSet(mode, modes::Mode::kQa)) {
      mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kCosPa, HistTable)), v0candidate.cosPa());
      mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kDecayDauDca, HistTable)), v0candidate.dauDca());
      mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kDecayVtxX, HistTable)), v0candidate.decayVtxX());
      mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kDecayVtxY, HistTable)), v0candidate.decayVtxY());
      mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kDecayVtxZ, HistTable)), v0candidate.decayVtxZ());
      mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kDecayVtx, HistTable)), v0candidate.decayVtx());
      mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kTransRadius, HistTable)), v0candidate.transRadius());
      mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kPtVsEta, HistTable)), v0candidate.pt(), v0candidate.eta());
      mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kPtVsPhi, HistTable)), v0candidate.pt(), v0candidate.phi());
      mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kPhiVsEta, HistTable)), v0candidate.phi(), v0candidate.eta());
      mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kPtVsCosPa, HistTable)), v0candidate.pt(), v0candidate.cosPa());

      if constexpr (modes::isEqual(v0, modes::V0::kLambda) || modes::isEqual(v0, modes::V0::kAntiLambda)) {
        float massLambda, massAntiLambda;
        if (v0candidate.sign() > 0) {
          massLambda = v0candidate.mass();
          massAntiLambda = v0candidate.massAnti();
        } else {
          massLambda = v0candidate.massAnti();
          massAntiLambda = v0candidate.mass();
        }
        mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kMassLambda, HistTable)), massLambda);
        mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kMassAntiLambda, HistTable)), massAntiLambda);
        mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kMassK0short, HistTable)), v0candidate.massK0short());
        mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kPtVsLambdaMass, HistTable)), v0candidate.pt(), massLambda);
        mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kPtVsAntiLambdaMass, HistTable)), v0candidate.pt(), massAntiLambda);
        mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kPtVsK0shortMass, HistTable)), v0candidate.pt(), v0candidate.massK0short());
        mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kLambdaMassVsAntiLambdaMass, HistTable)), massLambda, massAntiLambda);
        mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kK0shortMassVsLambdaMass, HistTable)), v0candidate.massK0short(), massLambda);
        mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kK0shortMassVsAntiLambdaMass, HistTable)), v0candidate.massK0short(), massAntiLambda);
      }

      if constexpr (modes::isEqual(v0, modes::V0::kK0short)) {
        mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kMassLambda, HistTable)), v0candidate.massLambda());
        mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kMassAntiLambda, HistTable)), v0candidate.massAntiLambda());
        mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kMassK0short, HistTable)), v0candidate.mass());
        mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kPtVsLambdaMass, HistTable)), v0candidate.pt(), v0candidate.massLambda());
        mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kPtVsAntiLambdaMass, HistTable)), v0candidate.pt(), v0candidate.massAntiLambda());
        mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kPtVsK0shortMass, HistTable)), v0candidate.pt(), v0candidate.mass());
        mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kLambdaMassVsAntiLambdaMass, HistTable)), v0candidate.massLambda(), v0candidate.massAntiLambda());
        mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kK0shortMassVsLambdaMass, HistTable)), v0candidate.mass(), v0candidate.massLambda());
        mHistogramRegistry->fill(HIST(v0Prefix) + HIST(QaDir) + HIST(getHistName(kK0shortMassVsAntiLambdaMass, HistTable)), v0candidate.mass(), v0candidate.massAntiLambda());
      }
    }
  }

 private:
  o2::framework::HistogramRegistry* mHistogramRegistry = nullptr;

  trackhistmanager::TrackHistManager<posDauPrefix, mode> mPosDauManager;
  trackhistmanager::TrackHistManager<negDauPrefix, mode> mNegDauManager;
};
}; // namespace v0histmanager
}; // namespace o2::analysis::femto
#endif // PWGCF_FEMTO_CORE_V0HISTMANAGER_H_
