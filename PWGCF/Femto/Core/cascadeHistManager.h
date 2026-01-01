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

/// \file cascadeHistManager.h
/// \brief histogram manager for cascade histograms
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#ifndef PWGCF_FEMTO_CORE_CASCADEHISTMANAGER_H_
#define PWGCF_FEMTO_CORE_CASCADEHISTMANAGER_H_

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
namespace cascadehistmanager
{
// enum for track histograms
enum CascadeHist {
  // analysis
  kPt,
  kEta,
  kPhi,
  kMass,
  kSign,
  // qa variables
  kMassXi,
  kMassOmega,
  kCosPa,
  kDecayDauDca,
  kTransRadius,
  // 2d qa
  kPtVsEta,
  kPtVsPhi,
  kPhiVsEta,
  kPtVsCosPa,
  kPtVsMassXi,
  kPtVsMassOmega,
  kMassXiVsMassOmega,
  kCascadeHistLast
};

template <const char* Prefix>
struct ConfCascadeBinning : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  o2::framework::ConfigurableAxis pt{"pt", {{600, 0, 6}}, "Pt"};
  o2::framework::ConfigurableAxis eta{"eta", {{300, -1.5, 1.5}}, "Eta"};
  o2::framework::ConfigurableAxis phi{"phi", {{720, 0, 1.f * o2::constants::math::TwoPI}}, "Phi"};
  o2::framework::ConfigurableAxis mass{"mass", {{1000, 1.f, 2.f}}, "Mass"};
  o2::framework::ConfigurableAxis sign{"sign", {{3, -1.5, 1.5}}, "Sign"};
};

constexpr const char PrefixXiBinning[] = "XiBinning";
using ConfXiBinning = ConfCascadeBinning<PrefixXiBinning>;
constexpr const char PrefixOmegaBinning[] = "OmegaBinning";
using ConfOmegaBinning = ConfCascadeBinning<PrefixOmegaBinning>;

template <const char* Prefix>
struct ConfCascadeQaBinning : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  o2::framework::Configurable<bool> plot2d{"plot2d", true, "Enable 2d Qa histograms"};
  o2::framework::ConfigurableAxis cosPa{"cosPa", {{100, 0.9, 1}}, "Cosine of poiting angle"};
  o2::framework::ConfigurableAxis dauDcaAtDecay{"dauDcaAtDecay", {{150, 0, 1.5}}, "Daughter DCA at decay vertex"};
  o2::framework::ConfigurableAxis transRadius{"transRadius", {{100, 0, 100}}, "Transverse radius"};
  o2::framework::ConfigurableAxis massXi{"massXi", {{400, 1.2f, 1.6f}}, "mass for antiparticle hypothesis"};
  o2::framework::ConfigurableAxis massOmega{"massOmega", {{400, 1.4f, 1.8f}}, "mass for antiparticle hypothesis"};
};

constexpr const char PrefixXiQaBinning[] = "XiQaBinning";
using ConfXiQaBinning = ConfCascadeQaBinning<PrefixXiQaBinning>;

constexpr const char PrefixOmegatQaBinning[] = "OmegaQaBinning";
using ConfOmegaQaBinning = ConfCascadeQaBinning<PrefixOmegatQaBinning>;

// must be in sync with enum TrackVariables
// the enum gives the correct index in the array
constexpr std::array<histmanager::HistInfo<CascadeHist>, kCascadeHistLast> HistTable = {
  {{kPt, o2::framework::kTH1F, "hPt", "Transverse Momentum; p_{T} (GeV/#it{c}); Entries"},
   {kEta, o2::framework::kTH1F, "hEta", "Pseudorapdity; #eta; Entries"},
   {kPhi, o2::framework::kTH1F, "hPhi", "Azimuthal angle; #varphi; Entries"},
   {kMass, o2::framework::kTH1F, "hMass", "Invariant Mass; m_{Inv} (GeV/#it{c}^{2}); Entries"},
   {kSign, o2::framework::kTH1F, "hSign", "Sign (-1 -> antiparticle, 0 -> self conjugate, +1 -> particle); sign; Entries"},
   {kMassXi, o2::framework::kTH1F, "hMassXi", "Mass #Xi; m_{#Lambda#pi} (GeV/#it{c}^{2}); Entries"},
   {kMassOmega, o2::framework::kTH1F, "hMassOmega", "mass #Omega; m_{#LambdaK} (GeV/#it{c}^{2}); Entries"},
   {kCosPa, o2::framework::kTH1F, "hCosPa", "Cosine of pointing angle; coa(#alpha); Entries"},
   {kDecayDauDca, o2::framework::kTH1F, "hDauDca", "Daughter DCA at decay vertex ; DCA_{Decay vertex} (cm); Entries"},
   {kTransRadius, o2::framework::kTH1F, "hTransRadius", "Transverse radius ; r_{xy} (cm); Entries"},
   {kPtVsEta, o2::framework::kTH2F, "hPtVsEta", "p_{T} vs #eta; p_{T} (GeV/#it{c}) ; #eta"},
   {kPtVsPhi, o2::framework::kTH2F, "hPtVsPhi", "p_{T} vs #varphi; p_{T} (GeV/#it{c}) ; #varphi"},
   {kPhiVsEta, o2::framework::kTH2F, "hPhiVsEta", "#varphi vs #eta; #varphi ; #eta"},
   {kPtVsCosPa, o2::framework::kTH2F, "hPtVsCosPa", "Cosine of poiting angle vs p_{T}; cos(#alpha); p_{T} (GeV/#it{c})"},
   {kPtVsMassXi, o2::framework::kTH2F, "hPtVsMassXi", "p_{T} vs mass #Xi; p_{T} (GeV/#it{c}); m_{#Lambda#pi} (GeV/#it{c}^{2})"},
   {kPtVsMassOmega, o2::framework::kTH2F, "hPtVsMassOmega", "p_{T} vs mass #Omega; p_{T} (GeV/#it{c}); m_{#LambdaK} (GeV/#it{c}^{2})"},
   {kMassXiVsMassOmega, o2::framework::kTH2F, "hMassXiVsMassOmega", "mass #Xi vs mass #Omega; m_{#Lambda#pi} (GeV/#it{c}^{2}); m_{#LambdaK} (GeV/#it{c}^{2})"}}};

template <typename T>
auto makeCascadeHistSpecMap(const T& confBinningAnalysis)
{
  return std::map<CascadeHist, std::vector<framework::AxisSpec>>{
    {kPt, {confBinningAnalysis.pt}},
    {kEta, {confBinningAnalysis.eta}},
    {kPhi, {confBinningAnalysis.phi}},
    {kMass, {confBinningAnalysis.mass}},
    {kSign, {confBinningAnalysis.sign}}};
}

template <typename T1, typename T2>
std::map<CascadeHist, std::vector<framework::AxisSpec>> makeCascadeQaHistSpecMap(T1 const& confBinningAnalysis, T2 const& confBinningQa)
{
  return std::map<CascadeHist, std::vector<framework::AxisSpec>>{
    {kPt, {confBinningAnalysis.pt}},
    {kEta, {confBinningAnalysis.eta}},
    {kPhi, {confBinningAnalysis.phi}},
    {kMass, {confBinningAnalysis.mass}},
    {kSign, {confBinningAnalysis.sign}},
    {kCosPa, {confBinningQa.cosPa}},
    {kDecayDauDca, {confBinningQa.dauDcaAtDecay}},
    {kTransRadius, {confBinningQa.transRadius}},
    {kPtVsEta, {confBinningAnalysis.pt, confBinningAnalysis.eta}},
    {kPtVsPhi, {confBinningAnalysis.pt, confBinningAnalysis.phi}},
    {kPhiVsEta, {confBinningAnalysis.phi, confBinningAnalysis.eta}},
    {kPtVsCosPa, {confBinningAnalysis.pt, confBinningQa.cosPa}},
    {kMassXi, {confBinningQa.massXi}},
    {kMassOmega, {confBinningQa.massOmega}},
    {kPtVsMassXi, {confBinningAnalysis.pt, confBinningQa.massXi}},
    {kPtVsMassOmega, {confBinningAnalysis.pt, confBinningQa.massOmega}},
    {kMassXiVsMassOmega, {confBinningQa.massXi, confBinningQa.massOmega}}};
};

constexpr char PrefixXiQa[] = "XiQA/";
constexpr char PrefixXi[] = "Xi/";
constexpr char PrefixOmegaQa[] = "OmegaQa/";
constexpr char PrefixOmega[] = "Omega/";

constexpr char PrefixLambdaCascade[] = "LambdaCascadeQa/";

constexpr std::string_view AnalysisDir = "Kinematics/";
constexpr std::string_view QaDir = "QA/";

constexpr int AbsChargeDaughters = 1;

/// \class FemtoDreamEventHisto
/// \brief Class for histogramming event properties
// template <femtomodes::Mode mode>
template <const char* cascadePrefix,
          const char* bachelorPrefix,
          const char* posDauPrefix,
          const char* negDauPrefix,
          modes::Cascade cascade>
class CascadeHistManager
{
 public:
  CascadeHistManager() = default;
  ~CascadeHistManager() = default;

  template <modes::Mode mode, typename T>
  void init(o2::framework::HistogramRegistry* registry,
            std::map<CascadeHist, std::vector<o2::framework::AxisSpec>> const& cascadeSpecs,
            T const& ConfCascadeSelection,
            std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> const& BachelorSpecs,
            std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> const& PosDauSpecs,
            std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> const& NegDauSpecs)
  {
    mHistogramRegistry = registry;
    mPdgCode = std::abs(ConfCascadeSelection.pdgCodeAbs.value);

    int bachelorPdgCodeAbs = 0;
    int posDauPdgCodeAbs = 0;
    int negDauPdgCodeAbs = 0;
    const int absCharge = 1;
    int signBachelor = 0;
    const int signPlus = 1;
    const int signMinus = -1;

    if (mPdgCode == PDG_t::kXiMinus) {
      if (ConfCascadeSelection.sign.value < 0) {
        bachelorPdgCodeAbs = std::abs(PDG_t::kPiMinus);
        signBachelor = -1;
        posDauPdgCodeAbs = std::abs(PDG_t::kProton);
        negDauPdgCodeAbs = std::abs(PDG_t::kPiMinus);
      } else {
        mPdgCode = -1 * mPdgCode; // Xi+ has negative pdg code
        bachelorPdgCodeAbs = std::abs(PDG_t::kPiPlus);
        signBachelor = 1;
        posDauPdgCodeAbs = std::abs(PDG_t::kPiPlus);
        negDauPdgCodeAbs = std::abs(PDG_t::kProtonBar);
      }
    } else if (mPdgCode == PDG_t::kOmegaMinus) {
      if (ConfCascadeSelection.sign.value < 0) {
        bachelorPdgCodeAbs = std::abs(PDG_t::kKMinus);
        signBachelor = -1;
        posDauPdgCodeAbs = std::abs(PDG_t::kProton);
        negDauPdgCodeAbs = std::abs(PDG_t::kPiMinus);
      } else {
        mPdgCode = -1 * mPdgCode; // Omega+ has negative pdg code
        bachelorPdgCodeAbs = std::abs(PDG_t::kKPlus);
        signBachelor = 1;
        posDauPdgCodeAbs = std::abs(PDG_t::kPiPlus);
        negDauPdgCodeAbs = std::abs(PDG_t::kProtonBar);
      }
    } else {
      LOG(fatal) << "PDG code for Cascade has to be either Xi or Omega";
    }

    mBachelorManager.template init<mode>(registry, BachelorSpecs, absCharge, signBachelor, bachelorPdgCodeAbs);
    mPosDauManager.template init<mode>(registry, PosDauSpecs, absCharge, signPlus, posDauPdgCodeAbs);
    mNegDauManager.template init<mode>(registry, NegDauSpecs, absCharge, signMinus, negDauPdgCodeAbs);

    if constexpr (modes::isFlagSet(mode, modes::Mode::kAnalysis)) {
      initAnalysis(cascadeSpecs);
    }
  }

  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5>
  void init(o2::framework::HistogramRegistry* registry,
            std::map<CascadeHist, std::vector<o2::framework::AxisSpec>> const& cascadeSpecs,
            T1 const& ConfCascadeSelection,
            T2 const& ConfCascadeBinningQa,
            std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> const& BachelorSpecs,
            T3 const& ConfBachelorQaBinning,
            std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> const& PosDauSpecs,
            T4& ConfPosDauQaBinning,
            std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> const& NegDauSpecs,
            T5& ConfNegDauQaBinning)
  {
    mHistogramRegistry = registry;
    mPdgCode = std::abs(ConfCascadeSelection.pdgCodeAbs.value);
    this->enableOptionalHistograms(ConfCascadeBinningQa);

    int bachelorPdgCodeAbs = 0;
    int posDauPdgCodeAbs = 0;
    int negDauPdgCodeAbs = 0;
    const int absCharge = 1;
    int signBachelor = 0;
    const int signPlus = 1;
    const int signMinus = -1;

    if (mPdgCode == PDG_t::kXiMinus) {
      if (ConfCascadeSelection.sign.value < 0) {
        bachelorPdgCodeAbs = std::abs(PDG_t::kPiMinus);
        signBachelor = -1;
        posDauPdgCodeAbs = std::abs(PDG_t::kProton);
        negDauPdgCodeAbs = std::abs(PDG_t::kPiMinus);
      } else {
        mPdgCode = -1 * mPdgCode; // Xi+ has negative pdg code
        bachelorPdgCodeAbs = std::abs(PDG_t::kPiPlus);
        signBachelor = 1;
        posDauPdgCodeAbs = std::abs(PDG_t::kPiPlus);
        negDauPdgCodeAbs = std::abs(PDG_t::kProtonBar);
      }
    } else if (mPdgCode == PDG_t::kOmegaMinus) {
      if (ConfCascadeSelection.sign.value < 0) {
        bachelorPdgCodeAbs = std::abs(PDG_t::kKMinus);
        signBachelor = -1;
        posDauPdgCodeAbs = std::abs(PDG_t::kProton);
        negDauPdgCodeAbs = std::abs(PDG_t::kPiMinus);
      } else {
        mPdgCode = -1 * mPdgCode; // Omega+ has negative pdg code
        bachelorPdgCodeAbs = std::abs(PDG_t::kKPlus);
        signBachelor = 1;
        posDauPdgCodeAbs = std::abs(PDG_t::kPiPlus);
        negDauPdgCodeAbs = std::abs(PDG_t::kProtonBar);
      }
    } else {
      LOG(fatal) << "PDG code for Cascade has to be either Xi or Omega";
    }

    mBachelorManager.template init<mode>(registry, BachelorSpecs, absCharge, signBachelor, bachelorPdgCodeAbs, ConfBachelorQaBinning);
    mPosDauManager.template init<mode>(registry, PosDauSpecs, absCharge, signPlus, posDauPdgCodeAbs, ConfPosDauQaBinning);
    mNegDauManager.template init<mode>(registry, NegDauSpecs, absCharge, signMinus, negDauPdgCodeAbs, ConfNegDauQaBinning);

    if constexpr (modes::isFlagSet(mode, modes::Mode::kAnalysis)) {
      initAnalysis(cascadeSpecs);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kQa)) {
      initQa(cascadeSpecs);
    }
  }

  template <modes::Mode mode, typename T1, typename T2>
  void fill(T1 const& cascadeCandidate, T2 const& tracks)
  {
    // this used to work, still under investigation
    // auto bachelor = cascadeCandidate.template bachelor_as<T2>();
    // auto posDaughter = cascadeCandidate.template posDau_as<T2>();
    // auto negDaughter = cascadeCandidate.template negDau_as<T2>();
    auto posDaughter = tracks.rawIteratorAt(cascadeCandidate.posDauId() - tracks.offset());
    mPosDauManager.template fill<mode>(posDaughter, tracks);
    auto negDaughter = tracks.rawIteratorAt(cascadeCandidate.negDauId() - tracks.offset());
    mNegDauManager.template fill<mode>(negDaughter, tracks);
    auto bachelor = tracks.rawIteratorAt(cascadeCandidate.bachelorId() - tracks.offset());
    mBachelorManager.template fill<mode>(bachelor, tracks);

    if constexpr (modes::isFlagSet(mode, modes::Mode::kAnalysis)) {
      fillAnalysis(cascadeCandidate);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kQa)) {
      fillQa(cascadeCandidate);
    }
  }

 private:
  template <typename T>
  void enableOptionalHistograms(T const& CascadeConfBinningQa)
  {
    mPlot2d = CascadeConfBinningQa.plot2d.value;
  }

  void initAnalysis(std::map<CascadeHist, std::vector<o2::framework::AxisSpec>> const& cascadeSpecs)
  {
    std::string analysisDir = std::string(cascadePrefix) + std::string(AnalysisDir);
    mHistogramRegistry->add(analysisDir + getHistNameV2(kPt, HistTable), getHistDesc(kPt, HistTable), getHistType(kPt, HistTable), {cascadeSpecs.at(kPt)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kEta, HistTable), getHistDesc(kEta, HistTable), getHistType(kEta, HistTable), {cascadeSpecs.at(kEta)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kPhi, HistTable), getHistDesc(kPhi, HistTable), getHistType(kPhi, HistTable), {cascadeSpecs.at(kPhi)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kMass, HistTable), getHistDesc(kMass, HistTable), getHistType(kMass, HistTable), {cascadeSpecs.at(kMass)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kSign, HistTable), getHistDesc(kSign, HistTable), getHistType(kSign, HistTable), {cascadeSpecs.at(kSign)});
  }

  void initQa(std::map<CascadeHist, std::vector<o2::framework::AxisSpec>> const& cascadeSpecs)
  {
    std::string qaDir = std::string(cascadePrefix) + std::string(QaDir);
    mHistogramRegistry->add(qaDir + getHistNameV2(kCosPa, HistTable), getHistDesc(kCosPa, HistTable), getHistType(kCosPa, HistTable), {cascadeSpecs.at(kCosPa)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kDecayDauDca, HistTable), getHistDesc(kDecayDauDca, HistTable), getHistType(kDecayDauDca, HistTable), {cascadeSpecs.at(kDecayDauDca)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kTransRadius, HistTable), getHistDesc(kTransRadius, HistTable), getHistType(kTransRadius, HistTable), {cascadeSpecs.at(kTransRadius)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kMassXi, HistTable), getHistDesc(kMassXi, HistTable), getHistType(kMassXi, HistTable), {cascadeSpecs.at(kMassXi)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kMassOmega, HistTable), getHistDesc(kMassOmega, HistTable), getHistType(kMassOmega, HistTable), {cascadeSpecs.at(kMassOmega)});

    if (mPlot2d) {
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsEta, HistTable), getHistDesc(kPtVsEta, HistTable), getHistType(kPtVsEta, HistTable), {cascadeSpecs.at(kPtVsEta)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsPhi, HistTable), getHistDesc(kPtVsPhi, HistTable), getHistType(kPtVsPhi, HistTable), {cascadeSpecs.at(kPtVsPhi)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPhiVsEta, HistTable), getHistDesc(kPhiVsEta, HistTable), getHistType(kPhiVsEta, HistTable), {cascadeSpecs.at(kPhiVsEta)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsCosPa, HistTable), getHistDesc(kPtVsCosPa, HistTable), getHistType(kPtVsCosPa, HistTable), {cascadeSpecs.at(kPtVsCosPa)});

      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsMassXi, HistTable), getHistDesc(kPtVsMassXi, HistTable), getHistType(kPtVsMassXi, HistTable), {cascadeSpecs.at(kPtVsMassXi)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsMassOmega, HistTable), getHistDesc(kPtVsMassOmega, HistTable), getHistType(kPtVsMassOmega, HistTable), {cascadeSpecs.at(kPtVsMassOmega)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kMassXiVsMassOmega, HistTable), getHistDesc(kMassXiVsMassOmega, HistTable), getHistType(kMassXiVsMassOmega, HistTable), {cascadeSpecs.at(kMassXiVsMassOmega)});
    }
  }

  template <typename T>
  void fillAnalysis(T const& cascadeCandidate)
  {
    mHistogramRegistry->fill(HIST(cascadePrefix) + HIST(AnalysisDir) + HIST(getHistName(kPt, HistTable)), cascadeCandidate.pt());
    mHistogramRegistry->fill(HIST(cascadePrefix) + HIST(AnalysisDir) + HIST(getHistName(kEta, HistTable)), cascadeCandidate.eta());
    mHistogramRegistry->fill(HIST(cascadePrefix) + HIST(AnalysisDir) + HIST(getHistName(kPhi, HistTable)), cascadeCandidate.phi());
    mHistogramRegistry->fill(HIST(cascadePrefix) + HIST(AnalysisDir) + HIST(getHistName(kMass, HistTable)), cascadeCandidate.mass());
    mHistogramRegistry->fill(HIST(cascadePrefix) + HIST(AnalysisDir) + HIST(getHistName(kSign, HistTable)), cascadeCandidate.sign());
  }

  template <typename T>
  void fillQa(T const& cascadeCandidate)
  {
    float massXi = 0.f;
    float massOmega = 0.f;
    if constexpr (modes::isEqual(cascade, modes::Cascade::kXi)) {
      massXi = cascadeCandidate.mass();
      massOmega = cascadeCandidate.massOmega();
    }
    if constexpr (modes::isEqual(cascade, modes::Cascade::kOmega)) {
      massXi = cascadeCandidate.massXi();
      massOmega = cascadeCandidate.mass();
    }
    mHistogramRegistry->fill(HIST(cascadePrefix) + HIST(QaDir) + HIST(getHistName(kCosPa, HistTable)), cascadeCandidate.cascadeCosPa());
    mHistogramRegistry->fill(HIST(cascadePrefix) + HIST(QaDir) + HIST(getHistName(kDecayDauDca, HistTable)), cascadeCandidate.cascadeDauDca());
    mHistogramRegistry->fill(HIST(cascadePrefix) + HIST(QaDir) + HIST(getHistName(kTransRadius, HistTable)), cascadeCandidate.cascadeTransRadius());
    mHistogramRegistry->fill(HIST(cascadePrefix) + HIST(QaDir) + HIST(getHistName(kMassXi, HistTable)), massXi);
    mHistogramRegistry->fill(HIST(cascadePrefix) + HIST(QaDir) + HIST(getHistName(kMassOmega, HistTable)), massOmega);

    if (mPlot2d) {
      mHistogramRegistry->fill(HIST(cascadePrefix) + HIST(QaDir) + HIST(getHistName(kPtVsEta, HistTable)), cascadeCandidate.pt(), cascadeCandidate.eta());
      mHistogramRegistry->fill(HIST(cascadePrefix) + HIST(QaDir) + HIST(getHistName(kPtVsPhi, HistTable)), cascadeCandidate.pt(), cascadeCandidate.phi());
      mHistogramRegistry->fill(HIST(cascadePrefix) + HIST(QaDir) + HIST(getHistName(kPhiVsEta, HistTable)), cascadeCandidate.phi(), cascadeCandidate.eta());
      mHistogramRegistry->fill(HIST(cascadePrefix) + HIST(QaDir) + HIST(getHistName(kPtVsCosPa, HistTable)), cascadeCandidate.pt(), cascadeCandidate.cascadeCosPa());
      mHistogramRegistry->fill(HIST(cascadePrefix) + HIST(QaDir) + HIST(getHistName(kPtVsMassXi, HistTable)), cascadeCandidate.pt(), massXi);
      mHistogramRegistry->fill(HIST(cascadePrefix) + HIST(QaDir) + HIST(getHistName(kPtVsMassOmega, HistTable)), cascadeCandidate.pt(), massOmega);
      mHistogramRegistry->fill(HIST(cascadePrefix) + HIST(QaDir) + HIST(getHistName(kMassXiVsMassOmega, HistTable)), massXi, massOmega);
    }
  }

  o2::framework::HistogramRegistry* mHistogramRegistry = nullptr;
  bool mPlot2d = true;
  trackhistmanager::TrackHistManager<bachelorPrefix> mBachelorManager;
  trackhistmanager::TrackHistManager<posDauPrefix> mPosDauManager;
  trackhistmanager::TrackHistManager<negDauPrefix> mNegDauManager;
  int mPdgCode = 0;
};
}; // namespace cascadehistmanager
}; // namespace o2::analysis::femto
#endif // PWGCF_FEMTO_CORE_CASCADEHISTMANAGER_H_
