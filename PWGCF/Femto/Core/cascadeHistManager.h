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
  kLambdaCosPa,
  kLambdaDauDca,
  kLambdaTransRadius,
  kLambdaDcaToPv,
  // 2d qa
  kPtVsEta,
  kPtVsPhi,
  kPhiVsEta,
  kPtVsCosPa,
  kPtVsMassXi,
  kPtVsMassOmega,
  kMassXiVsMassOmega,
  // mc
  kOrigin,
  kPdg,
  kPdgMother,
  kPdgPartonicMother,
  kTruePtVsPt,
  kTrueEtaVsEta,
  kTruePhiVsPhi,
  // histograms for fraction estimation of v0s
  kNoMcParticle,
  kPrimary,
  kFromWrongCollision,
  kFromMaterial,
  kMissidentified,
  kSecondary1,
  kSecondary2,
  kSecondary3,
  kSecondaryOther,
  kCascadeHistLast
};

constexpr std::size_t MaxSecondary = 3;

template <const char* Prefix>
struct ConfCascadeBinning : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  o2::framework::ConfigurableAxis pt{"pt", {{600, 0, 6}}, "Pt"};
  o2::framework::ConfigurableAxis eta{"eta", {{300, -1.5, 1.5}}, "Eta"};
  o2::framework::ConfigurableAxis phi{"phi", {{720, 0, 1.f * o2::constants::math::TwoPI}}, "Phi"};
  o2::framework::ConfigurableAxis mass{"mass", {{1000, 1.f, 2.f}}, "Mass"};
  o2::framework::ConfigurableAxis sign{"sign", {{3, -1.5, 1.5}}, "Sign"};
  o2::framework::ConfigurableAxis pdgCodes{"pdgCodes", {{8001, -4000.5, 4000.5}}, "MC ONLY: PDG codes of reconstructed cascades"};
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
  o2::framework::ConfigurableAxis lambdaCosPa{"lambdaCosPa", {{100, 0.9, 1}}, "Cosine of poiting angle of daughter lambda"};
  o2::framework::ConfigurableAxis lambdaDauDca{"lambdaDauDca", {{150, 0, 1.0}}, "DCA of lambda daughters at lambda decay vertex"};
  o2::framework::ConfigurableAxis lambdaTransRadius{"lambdaTransRadius", {{100, 0, 100}}, "DCA of lambda daughters at lambda decay vertex"};
  o2::framework::ConfigurableAxis lambdaDcaToPv{"lambdaDcaToPv", {{100, 0, 200}}, "DCA of lambda daughter from primary vertex"};
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
   {kCosPa, o2::framework::kTH1F, "hCosPa", "Cosine of pointing angle; cos(#alpha); Entries"},
   {kDecayDauDca, o2::framework::kTH1F, "hDauDca", "Daughter DCA at decay vertex ; DCA_{Decay vertex} (cm); Entries"},
   {kTransRadius, o2::framework::kTH1F, "hTransRadius", "Transverse radius ; r_{xy} (cm); Entries"},
   {kLambdaCosPa, o2::framework::kTH1F, "hLambdaCosPa", "Cosine of poiting angle of daughter lambda ; cos_{#Lambda}(#alpha); Entries"},
   {kLambdaDauDca, o2::framework::kTH1F, "hLambdaDauDca", "Daughter DCA at #Lambda decay vertex ; DCA_{#Lambda decay vertex} (cm); Entries"},
   {kLambdaTransRadius, o2::framework::kTH1F, "hLambdaTransRadius", "Transverse radius of daughter #Lambda ; r_{xy} (cm); Entries"},
   {kLambdaDcaToPv, o2::framework::kTH1F, "hLambdaDcaToPv", "DCA to primary vertex of daughter #Lambda ; DCA (cm); Entries"},
   {kPtVsEta, o2::framework::kTH2F, "hPtVsEta", "p_{T} vs #eta; p_{T} (GeV/#it{c}) ; #eta"},
   {kPtVsPhi, o2::framework::kTH2F, "hPtVsPhi", "p_{T} vs #varphi; p_{T} (GeV/#it{c}) ; #varphi"},
   {kPhiVsEta, o2::framework::kTH2F, "hPhiVsEta", "#varphi vs #eta; #varphi ; #eta"},
   {kPtVsCosPa, o2::framework::kTH2F, "hPtVsCosPa", "Cosine of poiting angle vs p_{T}; cos(#alpha); p_{T} (GeV/#it{c})"},
   {kPtVsMassXi, o2::framework::kTH2F, "hPtVsMassXi", "p_{T} vs mass #Xi; p_{T} (GeV/#it{c}); m_{#Lambda#pi} (GeV/#it{c}^{2})"},
   {kPtVsMassOmega, o2::framework::kTH2F, "hPtVsMassOmega", "p_{T} vs mass #Omega; p_{T} (GeV/#it{c}); m_{#LambdaK} (GeV/#it{c}^{2})"},
   {kMassXiVsMassOmega, o2::framework::kTH2F, "hMassXiVsMassOmega", "mass #Xi vs mass #Omega; m_{#Lambda#pi} (GeV/#it{c}^{2}); m_{#LambdaK} (GeV/#it{c}^{2})"},
   {kOrigin, o2::framework::kTH1F, "hOrigin", "Status Codes (=Origin); Status Code; Entries"},
   {kPdg, o2::framework::kTH1F, "hPdg", "PDG Codes of reconstructed v0; PDG Code; Entries"},
   {kPdgMother, o2::framework::kTH1F, "hPdgMother", "PDG Codes of mother of reconstructed v0; PDG Code; Entries"},
   {kPdgPartonicMother, o2::framework::kTH1F, "hPdgPartonicMother", "PDG Codes of partonic mother of reconstructed v0; PDG Code; Entries"},
   {kTruePtVsPt, o2::framework::kTH2F, "hTruePtVsPt", "True transverse momentum vs transverse momentum; p_{T,True} (GeV/#it{c}); p_{T,True} (GeV/#it{c})"},
   {kTrueEtaVsEta, o2::framework::kTH2F, "hTrueEtaVsEta", "True pseudorapdity vs pseudorapdity; #eta_{True}; #eta"},
   {kTruePhiVsPhi, o2::framework::kTH2F, "hTruePhiVsPhi", "True azimuthal angle vs azimuthal angle; #varphi_{True}; #varphi"},
   {kNoMcParticle, o2::framework::kTH2F, "hNoMcParticle", "Wrongly reconstructed particles; p_{T} (GeV/#it{c}); cos(#alpha)"},
   {kPrimary, o2::framework::kTH2F, "hPrimary", "Primary particles; p_{T} (GeV/#it{c}); cos(#alpha)"},
   {kFromWrongCollision, o2::framework::kTH2F, "hFromWrongCollision", "Particles associated to wrong collision; p_{T} (GeV/#it{c}); cos(#alpha)"},
   {kFromMaterial, o2::framework::kTH2F, "hFromMaterial", "Particles from material; p_{T} (GeV/#it{c}); cos(#alpha)"},
   {kMissidentified, o2::framework::kTH2F, "hMissidentified", "Missidentified particles (fake/wrong PDG code); p_{T} (GeV/#it{c}); cos(#alpha)"},
   {kSecondary1, o2::framework::kTH2F, "hFromSecondary1", "Particles from secondary decay; p_{T} (GeV/#it{c}); cos(#alpha)"},
   {kSecondary2, o2::framework::kTH2F, "hFromSecondary2", "Particles from seconary decay; p_{T} (GeV/#it{c}); cos(#alpha)"},
   {kSecondary3, o2::framework::kTH2F, "hFromSecondary3", "Particles from seconary decay; p_{T} (GeV/#it{c}); cos(#alpha)"},
   {kSecondaryOther, o2::framework::kTH2F, "hFromSecondaryOther", "Particles from every other seconary decay; p_{T} (GeV/#it{c}); cos(#alpha)"}},
};

#define CASCADE_HIST_ANALYSIS_MAP(conf) \
  {kPt, {conf.pt}},                     \
    {kEta, {conf.eta}},                 \
    {kPhi, {conf.phi}},                 \
    {kMass, {conf.mass}},               \
    {kSign, {conf.sign}},

#define CASCADE_HIST_MC_MAP(conf)          \
  {kTruePtVsPt, {conf.pt, conf.pt}},       \
    {kTrueEtaVsEta, {conf.eta, conf.eta}}, \
    {kTruePhiVsPhi, {conf.phi, conf.phi}}, \
    {kPdg, {conf.pdgCodes}},               \
    {kPdgMother, {conf.pdgCodes}},         \
    {kPdgPartonicMother, {conf.pdgCodes}},

#define CASCADE_HIST_QA_MAP(confAnalysis, confQa)          \
  {kPt, {confAnalysis.pt}},                                \
    {kEta, {confAnalysis.eta}},                            \
    {kPhi, {confAnalysis.phi}},                            \
    {kMass, {confAnalysis.mass}},                          \
    {kSign, {confAnalysis.sign}},                          \
    {kCosPa, {confQa.cosPa}},                              \
    {kDecayDauDca, {confQa.dauDcaAtDecay}},                \
    {kTransRadius, {confQa.transRadius}},                  \
    {kLambdaCosPa, {confQa.lambdaCosPa}},                  \
    {kLambdaDauDca, {confQa.lambdaDauDca}},                \
    {kLambdaTransRadius, {confQa.lambdaTransRadius}},      \
    {kLambdaDcaToPv, {confQa.lambdaDcaToPv}},              \
    {kPtVsEta, {confAnalysis.pt, confAnalysis.eta}},       \
    {kPtVsPhi, {confAnalysis.pt, confAnalysis.phi}},       \
    {kPhiVsEta, {confAnalysis.phi, confAnalysis.eta}},     \
    {kPtVsCosPa, {confAnalysis.pt, confQa.cosPa}},         \
    {kMassXi, {confQa.massXi}},                            \
    {kMassOmega, {confQa.massOmega}},                      \
    {kPtVsMassXi, {confAnalysis.pt, confQa.massXi}},       \
    {kPtVsMassOmega, {confAnalysis.pt, confQa.massOmega}}, \
    {kMassXiVsMassOmega, {confQa.massXi, confQa.massOmega}},

#define CASCADE_HIST_MC_QA_MAP(confAnalysis, confQa)        \
  {kNoMcParticle, {confAnalysis.pt, confQa.cosPa}},         \
    {kPrimary, {confAnalysis.pt, confQa.cosPa}},            \
    {kFromWrongCollision, {confAnalysis.pt, confQa.cosPa}}, \
    {kFromMaterial, {confAnalysis.pt, confQa.cosPa}},       \
    {kMissidentified, {confAnalysis.pt, confQa.cosPa}},     \
    {kSecondary1, {confAnalysis.pt, confQa.cosPa}},         \
    {kSecondary2, {confAnalysis.pt, confQa.cosPa}},         \
    {kSecondary3, {confAnalysis.pt, confQa.cosPa}},         \
    {kSecondaryOther, {confAnalysis.pt, confQa.cosPa}},

template <typename T>
auto makeCascadeHistSpecMap(const T& confBinningAnalysis)
{
  return std::map<CascadeHist, std::vector<o2::framework::AxisSpec>>{
    CASCADE_HIST_ANALYSIS_MAP(confBinningAnalysis)};
}

template <typename T>
auto makeCascadeMcHistSpecMap(const T& confBinningAnalysis)
{
  return std::map<CascadeHist, std::vector<o2::framework::AxisSpec>>{
    CASCADE_HIST_ANALYSIS_MAP(confBinningAnalysis)
      CASCADE_HIST_MC_MAP(confBinningAnalysis)};
}

template <typename T1, typename T2>
std::map<CascadeHist, std::vector<o2::framework::AxisSpec>> makeCascadeQaHistSpecMap(T1 const& confBinningAnalysis, T2 const& confBinningQa)
{
  return std::map<CascadeHist, std::vector<o2::framework::AxisSpec>>{
    CASCADE_HIST_ANALYSIS_MAP(confBinningAnalysis)
      CASCADE_HIST_QA_MAP(confBinningAnalysis, confBinningQa)};
}

template <typename T1, typename T2>
std::map<CascadeHist, std::vector<o2::framework::AxisSpec>> makeCascadeMcQaHistSpecMap(T1 const& confBinningAnalysis, T2 const& confBinningQa)
{
  return std::map<CascadeHist, std::vector<o2::framework::AxisSpec>>{
    CASCADE_HIST_ANALYSIS_MAP(confBinningAnalysis)
      CASCADE_HIST_QA_MAP(confBinningAnalysis, confBinningQa)
        CASCADE_HIST_MC_MAP(confBinningAnalysis)
          CASCADE_HIST_MC_QA_MAP(confBinningAnalysis, confBinningQa)};
}

#undef CASCADE_HIST_ANALYSIS_MAP
#undef CASCADE_HIST_MC_MAP
#undef CASCADE_HIST_QA_MAP
#undef CASCADE_HIST_MC_QA_MAP

constexpr char PrefixXiQa[] = "XiQA/";
constexpr char PrefixXi[] = "Xi/";
constexpr char PrefixOmegaQa[] = "OmegaQa/";
constexpr char PrefixOmega[] = "Omega/";

constexpr char PrefixLambdaCascade[] = "LambdaCascadeQa/";

constexpr std::string_view AnalysisDir = "Kinematics/";
constexpr std::string_view QaDir = "QA/";
constexpr std::string_view McDir = "MC/";

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
        mPdgCode = -1 * mPdgCode; // Xi+ has negative pdg code because it is the anti particle
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
        mPdgCode = -1 * mPdgCode; // Omega+ has negative pdg code because it is the anti particle
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
    if constexpr (modes::isFlagSet(mode, modes::Mode::kMc)) {
      this->initMc(cascadeSpecs);
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
    if constexpr (modes::isFlagSet(mode, modes::Mode::kMc)) {
      this->initMc(cascadeSpecs);
    }
  }

  template <modes::Mode mode, typename T1, typename T2>
  void fill(T1 const& cascadeCandidate, T2 const& tracks)
  {
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

  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5>
  void fill(T1 const& cascadeCandidate, T2 const& tracks, T3 const& mcParticles, T4 const& mcMothers, T5 const& mcPartonicMothers)
  {
    auto posDaughter = tracks.rawIteratorAt(cascadeCandidate.posDauId() - tracks.offset());
    mPosDauManager.template fill<mode>(posDaughter, tracks, mcParticles, mcMothers, mcPartonicMothers);
    auto negDaughter = tracks.rawIteratorAt(cascadeCandidate.negDauId() - tracks.offset());
    mNegDauManager.template fill<mode>(negDaughter, tracks, mcParticles, mcMothers, mcPartonicMothers);
    auto bachelor = tracks.rawIteratorAt(cascadeCandidate.bachelorId() - tracks.offset());
    mBachelorManager.template fill<mode>(bachelor, tracks, mcParticles, mcMothers, mcPartonicMothers);

    if constexpr (modes::isFlagSet(mode, modes::Mode::kAnalysis)) {
      this->fillAnalysis(cascadeCandidate);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kQa)) {
      this->fillQa(cascadeCandidate);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kMc)) {
      this->template fillMc<mode>(cascadeCandidate, mcParticles, mcMothers, mcPartonicMothers);
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
    mHistogramRegistry->add(qaDir + getHistNameV2(kLambdaCosPa, HistTable), getHistDesc(kLambdaCosPa, HistTable), getHistType(kLambdaCosPa, HistTable), {cascadeSpecs.at(kLambdaCosPa)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kLambdaDauDca, HistTable), getHistDesc(kLambdaDauDca, HistTable), getHistType(kLambdaDauDca, HistTable), {cascadeSpecs.at(kLambdaDauDca)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kLambdaTransRadius, HistTable), getHistDesc(kLambdaTransRadius, HistTable), getHistType(kLambdaTransRadius, HistTable), {cascadeSpecs.at(kLambdaTransRadius)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kLambdaDcaToPv, HistTable), getHistDesc(kLambdaDcaToPv, HistTable), getHistType(kLambdaDcaToPv, HistTable), {cascadeSpecs.at(kLambdaDcaToPv)});
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

  void initMc(std::map<CascadeHist, std::vector<o2::framework::AxisSpec>> const& V0Specs)
  {
    std::string mcDir = std::string(cascadePrefix) + std::string(McDir);
    mHistogramRegistry->add(mcDir + getHistNameV2(kTruePtVsPt, HistTable), getHistDesc(kTruePtVsPt, HistTable), getHistType(kTruePtVsPt, HistTable), {V0Specs.at(kTruePtVsPt)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kTrueEtaVsEta, HistTable), getHistDesc(kTrueEtaVsEta, HistTable), getHistType(kTrueEtaVsEta, HistTable), {V0Specs.at(kTrueEtaVsEta)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kTruePhiVsPhi, HistTable), getHistDesc(kTruePhiVsPhi, HistTable), getHistType(kTruePhiVsPhi, HistTable), {V0Specs.at(kTruePhiVsPhi)});

    // mc origin can be configured here
    const o2::framework::AxisSpec axisOrigin = {static_cast<int>(modes::McOrigin::kMcOriginLast), -0.5, static_cast<double>(modes::McOrigin::kMcOriginLast) - 0.5};
    mHistogramRegistry->add(mcDir + getHistNameV2(kOrigin, HistTable), getHistDesc(kOrigin, HistTable), getHistType(kOrigin, HistTable), {axisOrigin});
    mHistogramRegistry->get<TH1>(HIST(cascadePrefix) + HIST(McDir) + HIST(histmanager::getHistName(kOrigin, HistTable)))->GetXaxis()->SetBinLabel(1 + static_cast<int>(modes::McOrigin::kNoMcParticle), modes::mcOriginToString(modes::McOrigin::kNoMcParticle));
    mHistogramRegistry->get<TH1>(HIST(cascadePrefix) + HIST(McDir) + HIST(histmanager::getHistName(kOrigin, HistTable)))->GetXaxis()->SetBinLabel(1 + static_cast<int>(modes::McOrigin::kFromWrongCollision), modes::mcOriginToString(modes::McOrigin::kFromWrongCollision));
    mHistogramRegistry->get<TH1>(HIST(cascadePrefix) + HIST(McDir) + HIST(histmanager::getHistName(kOrigin, HistTable)))->GetXaxis()->SetBinLabel(1 + static_cast<int>(modes::McOrigin::kPhysicalPrimary), modes::mcOriginToString(modes::McOrigin::kPhysicalPrimary));
    mHistogramRegistry->get<TH1>(HIST(cascadePrefix) + HIST(McDir) + HIST(histmanager::getHistName(kOrigin, HistTable)))->GetXaxis()->SetBinLabel(1 + static_cast<int>(modes::McOrigin::kFromSecondaryDecay), modes::mcOriginToString(modes::McOrigin::kFromSecondaryDecay));
    mHistogramRegistry->get<TH1>(HIST(cascadePrefix) + HIST(McDir) + HIST(histmanager::getHistName(kOrigin, HistTable)))->GetXaxis()->SetBinLabel(1 + static_cast<int>(modes::McOrigin::kFromMaterial), modes::mcOriginToString(modes::McOrigin::kFromMaterial));
    mHistogramRegistry->get<TH1>(HIST(cascadePrefix) + HIST(McDir) + HIST(histmanager::getHistName(kOrigin, HistTable)))->GetXaxis()->SetBinLabel(1 + static_cast<int>(modes::McOrigin::kMissidentified), modes::mcOriginToString(modes::McOrigin::kMissidentified));

    mHistogramRegistry->add(mcDir + getHistNameV2(kPdg, HistTable), getHistDesc(kPdg, HistTable), getHistType(kPdg, HistTable), {V0Specs.at(kPdg)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kPdgMother, HistTable), getHistDesc(kPdgMother, HistTable), getHistType(kPdgMother, HistTable), {V0Specs.at(kPdgMother)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kPdgPartonicMother, HistTable), getHistDesc(kPdgPartonicMother, HistTable), getHistType(kPdgPartonicMother, HistTable), {V0Specs.at(kPdgPartonicMother)});

    if (mPlotOrigins) {
      mHistogramRegistry->add(mcDir + getHistNameV2(kNoMcParticle, HistTable), getHistDesc(kNoMcParticle, HistTable), getHistType(kNoMcParticle, HistTable), {V0Specs.at(kNoMcParticle)});
      mHistogramRegistry->add(mcDir + getHistNameV2(kPrimary, HistTable), getHistDesc(kPrimary, HistTable), getHistType(kPrimary, HistTable), {V0Specs.at(kPrimary)});
      mHistogramRegistry->add(mcDir + getHistNameV2(kFromWrongCollision, HistTable), getHistDesc(kFromWrongCollision, HistTable), getHistType(kFromWrongCollision, HistTable), {V0Specs.at(kFromWrongCollision)});
      mHistogramRegistry->add(mcDir + getHistNameV2(kFromMaterial, HistTable), getHistDesc(kFromMaterial, HistTable), getHistType(kFromMaterial, HistTable), {V0Specs.at(kFromMaterial)});
      mHistogramRegistry->add(mcDir + getHistNameV2(kMissidentified, HistTable), getHistDesc(kMissidentified, HistTable), getHistType(kMissidentified, HistTable), {V0Specs.at(kMissidentified)});

      if (mPlotNSecondaries >= histmanager::kSecondaryPlotLevel1) {
        mHistogramRegistry->add(mcDir + getHistNameV2(kSecondary1, HistTable), getHistDesc(kSecondary1, HistTable), getHistType(kSecondary1, HistTable), {V0Specs.at(kSecondary1)});
      }
      if (mPlotNSecondaries >= histmanager::kSecondaryPlotLevel2) {
        mHistogramRegistry->add(mcDir + getHistNameV2(kSecondary2, HistTable), getHistDesc(kSecondary2, HistTable), getHistType(kSecondary2, HistTable), {V0Specs.at(kSecondary2)});
      }
      if (mPlotNSecondaries >= histmanager::kSecondaryPlotLevel3) {
        mHistogramRegistry->add(mcDir + getHistNameV2(kSecondary3, HistTable), getHistDesc(kSecondary3, HistTable), getHistType(kSecondary3, HistTable), {V0Specs.at(kSecondary3)});
      }
      mHistogramRegistry->add(mcDir + getHistNameV2(kSecondaryOther, HistTable), getHistDesc(kSecondaryOther, HistTable), getHistType(kSecondaryOther, HistTable), {V0Specs.at(kSecondaryOther)});
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
    mHistogramRegistry->fill(HIST(cascadePrefix) + HIST(QaDir) + HIST(getHistName(kLambdaCosPa, HistTable)), cascadeCandidate.lambdaCosPa());
    mHistogramRegistry->fill(HIST(cascadePrefix) + HIST(QaDir) + HIST(getHistName(kLambdaDauDca, HistTable)), cascadeCandidate.lambdaDauDca());
    mHistogramRegistry->fill(HIST(cascadePrefix) + HIST(QaDir) + HIST(getHistName(kLambdaTransRadius, HistTable)), cascadeCandidate.lambdaTransRadius());
    mHistogramRegistry->fill(HIST(cascadePrefix) + HIST(QaDir) + HIST(getHistName(kLambdaDcaToPv, HistTable)), cascadeCandidate.lambdaDcaToPv());
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

  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4>
  void fillMc(T1 const& cascadeCandidate, T2 const& /*mcParticles*/, T3 const& /*mcMothers*/, T4 const& /*mcPartonicMothers*/)
  {
    // No MC Particle
    if (!cascadeCandidate.has_fMcParticle()) {
      mHistogramRegistry->fill(HIST(cascadePrefix) + HIST(McDir) + HIST(getHistName(kPdg, HistTable)), 0);
      mHistogramRegistry->fill(HIST(cascadePrefix) + HIST(McDir) + HIST(getHistName(kOrigin, HistTable)), static_cast<float>(modes::McOrigin::kNoMcParticle));
      if constexpr (modes::isFlagSet(mode, modes::Mode::kQa)) {
        if (mPlotOrigins) {
          mHistogramRegistry->fill(HIST(cascadePrefix) + HIST(McDir) + HIST(getHistName(kNoMcParticle, HistTable)), cascadeCandidate.pt(), cascadeCandidate.cascadeCosPa());
        }
      }
      return;
    }

    // Retrieve MC particle
    auto mcParticle = cascadeCandidate.template fMcParticle_as<T2>();

    // missidentifed particles are special case
    // whether a particle is missidentfied or not cannot be known by the producer so we check it here
    bool isMissidentified = mcParticle.pdgCode() != mPdgCode;

    mHistogramRegistry->fill(HIST(cascadePrefix) + HIST(McDir) + HIST(getHistName(kTruePtVsPt, HistTable)), mcParticle.pt(), cascadeCandidate.pt());
    mHistogramRegistry->fill(HIST(cascadePrefix) + HIST(McDir) + HIST(getHistName(kTrueEtaVsEta, HistTable)), mcParticle.eta(), cascadeCandidate.eta());
    mHistogramRegistry->fill(HIST(cascadePrefix) + HIST(McDir) + HIST(getHistName(kTruePhiVsPhi, HistTable)), mcParticle.phi(), cascadeCandidate.phi());
    if (isMissidentified) {
      mHistogramRegistry->fill(HIST(cascadePrefix) + HIST(McDir) + HIST(getHistName(kOrigin, HistTable)), static_cast<int>(modes::McOrigin::kMissidentified));
    } else {
      mHistogramRegistry->fill(HIST(cascadePrefix) + HIST(McDir) + HIST(getHistName(kOrigin, HistTable)), mcParticle.origin());
    }
    mHistogramRegistry->fill(HIST(cascadePrefix) + HIST(McDir) + HIST(getHistName(kPdg, HistTable)), mcParticle.pdgCode());

    // get mother
    if (cascadeCandidate.has_fMcMother()) {
      auto mother = cascadeCandidate.template fMcMother_as<T3>();
      mHistogramRegistry->fill(HIST(cascadePrefix) + HIST(McDir) + HIST(getHistName(kPdgMother, HistTable)), mother.pdgCode());
    } else {
      mHistogramRegistry->fill(HIST(cascadePrefix) + HIST(McDir) + HIST(getHistName(kPdgMother, HistTable)), 0);
    }

    // get partonic mother
    if (cascadeCandidate.has_fMcPartMoth()) {
      auto partonicMother = cascadeCandidate.template fMcPartMoth_as<T4>();
      mHistogramRegistry->fill(HIST(cascadePrefix) + HIST(McDir) + HIST(getHistName(kPdgPartonicMother, HistTable)), partonicMother.pdgCode());
    } else {
      mHistogramRegistry->fill(HIST(cascadePrefix) + HIST(McDir) + HIST(getHistName(kPdgPartonicMother, HistTable)), 0);
    }

    if constexpr (modes::isFlagSet(mode, modes::Mode::kQa)) {
      if (mPlotOrigins) {
        // check first if particle is missidentified
        if (isMissidentified) {
          // if it is, we fill it as such
          mHistogramRegistry->fill(HIST(cascadePrefix) + HIST(McDir) + HIST(getHistName(kMissidentified, HistTable)), cascadeCandidate.pt(), cascadeCandidate.cascadeCosPa());
        } else {
          // if not, we fill it acccoridng to its origin
          switch (static_cast<modes::McOrigin>(mcParticle.origin())) {
            case modes::McOrigin::kPhysicalPrimary:
              mHistogramRegistry->fill(HIST(cascadePrefix) + HIST(McDir) + HIST(getHistName(kPrimary, HistTable)), cascadeCandidate.pt(), cascadeCandidate.cascadeCosPa());
              break;
            case modes::McOrigin::kFromWrongCollision:
              mHistogramRegistry->fill(HIST(cascadePrefix) + HIST(McDir) + HIST(getHistName(kFromWrongCollision, HistTable)), cascadeCandidate.pt(), cascadeCandidate.cascadeCosPa());
              break;
            case modes::McOrigin::kFromMaterial:
              mHistogramRegistry->fill(HIST(cascadePrefix) + HIST(McDir) + HIST(getHistName(kFromMaterial, HistTable)), cascadeCandidate.pt(), cascadeCandidate.cascadeCosPa());
              break;
            case modes::McOrigin::kFromSecondaryDecay:
              if (cascadeCandidate.has_fMcMother()) {
                auto mother = cascadeCandidate.template fMcMother_as<T3>();
                int motherPdgCode = std::abs(mother.pdgCode());
                // Switch on PDG of the mother
                if (mPlotNSecondaries >= histmanager::kSecondaryPlotLevel1 && motherPdgCode == mPdgCodesSecondaryMother[0]) {
                  mHistogramRegistry->fill(HIST(cascadePrefix) + HIST(McDir) + HIST(getHistName(kSecondary1, HistTable)), cascadeCandidate.pt(), cascadeCandidate.cascadeCosPa());
                } else if (mPlotNSecondaries >= histmanager::kSecondaryPlotLevel2 && motherPdgCode == mPdgCodesSecondaryMother[1]) {
                  mHistogramRegistry->fill(HIST(cascadePrefix) + HIST(McDir) + HIST(getHistName(kSecondary2, HistTable)), cascadeCandidate.pt(), cascadeCandidate.cascadeCosPa());
                } else if (mPlotNSecondaries >= histmanager::kSecondaryPlotLevel3 && motherPdgCode == mPdgCodesSecondaryMother[2]) {
                  mHistogramRegistry->fill(HIST(cascadePrefix) + HIST(McDir) + HIST(getHistName(kSecondary3, HistTable)), cascadeCandidate.pt(), cascadeCandidate.cascadeCosPa());
                } else {
                  mHistogramRegistry->fill(HIST(cascadePrefix) + HIST(McDir) + HIST(getHistName(kSecondaryOther, HistTable)), cascadeCandidate.pt(), cascadeCandidate.cascadeCosPa());
                }
              }
              break;
            default:
              LOG(warn) << "Encounted partilce with unknown origin!";
              break;
          }
        }
      }
    }
  }

  o2::framework::HistogramRegistry* mHistogramRegistry = nullptr;
  int mPdgCode = 0;
  bool mPlotOrigins = false;
  int mPlotNSecondaries = 0;
  std::array<int, MaxSecondary> mPdgCodesSecondaryMother = {0};
  bool mPlot2d = false;

  trackhistmanager::TrackHistManager<bachelorPrefix> mBachelorManager;
  trackhistmanager::TrackHistManager<posDauPrefix> mPosDauManager;
  trackhistmanager::TrackHistManager<negDauPrefix> mNegDauManager;
};
}; // namespace cascadehistmanager
}; // namespace o2::analysis::femto
#endif // PWGCF_FEMTO_CORE_CASCADEHISTMANAGER_H_
