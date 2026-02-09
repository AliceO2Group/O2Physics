// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file HFFilterHelpers.h
/// \brief Header file with definition of variables, methods, and tables used in the HFFilter.cxx task
///
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN
/// \author Marcel Lesch <marcel.lesch@tum.de>, TUM
/// \author Alexandre Bigot <alexandre.bigot@cern.ch>, Strasbourg University
/// \author Biao Zhang <biao.zhang@cern.ch>, CCNU
/// \author Federica Zanone <federica.zanone@cern.ch>, Heidelberg University
/// \author Antonio Palasciano <antonio.palasciano@cern.ch>, INFN Bari

#ifndef EVENTFILTERING_PWGHF_HFFILTERHELPERS_H_
#define EVENTFILTERING_PWGHF_HFFILTERHELPERS_H_

#include "EventFiltering/filterTables.h"
//
#include "PWGHF/Core/SelectorCuts.h"
//
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DCAFitter/DCAFitterN.h>
#include <MathUtils/BetheBlochAleph.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/Array2D.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/Logger.h>

#include <Math/GenVector/Boost.h>
#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TAxis.h>
#include <TH1.h>
#include <TH3.h>

#include <Rtypes.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <tuple>
#include <vector>

namespace o2::aod
{

namespace hffilters
{

enum HfTriggers {
  kHighPt2P = 0,
  kHighPt3P,
  kBeauty3P,
  kBeauty4P,
  kFemto2P,
  kFemto3P,
  kDoubleCharm2P,
  kDoubleCharm3P,
  kDoubleCharmMix,
  kV0Charm2P,
  kV0Charm3P,
  kCharmBarToXiBach,
  kSigmaCPPK,
  kSigmaC0K0,
  kPhotonCharm2P,
  kPhotonCharm3P,
  kSingleCharm2P,
  kSingleCharm3P,
  kSingleNonPromptCharm2P,
  kSingleNonPromptCharm3P,
  kCharmBarToXi2Bach,
  kPrCharm2P,
  kBtoJPsiKa,
  kBtoJPsiKstar,
  kBtoJPsiPhi,
  kBtoJPsiPrKa,
  kBtoJPsiPi,
  kSigmaCPr,
  kNtriggersHF
};

enum charmParticles {
  kD0 = 0,
  kDplus,
  kDs,
  kLc,
  kXic,
  kNCharmParticles
};

enum beautyParticles {
  kBplus = 0,
  kB0toDStar,
  kBc,
  kB0,
  kBs,
  kLb,
  kXib,
  kNBeautyParticles
};

enum beautyToJPsiParticles {
  kBplusToJPsi = 0,
  kB0ToJPsi,
  kBsToJPsi,
  kLbToJPsi,
  kBcToJPsi,
  kNBeautyParticlesToJPsi
};

enum bachelorTrackSelection {
  kRejected = 0,
  kSoftPion,
  kForBeauty,
  kSoftPionForBeauty,
  kPionForCharmBaryon,
  kKaonForCharmBaryon,
  kProtonForCharmBaryon,
  kSoftPionForSigmaC
};

enum PIDSpecies {
  kEl = 0,
  kPi,
  kAntiPi,
  kKa,
  kAntiKa,
  kPr,
  kAntiPr,
  kDe,
  kAntiDe
};

enum trackSpecies {
  kProtonForFemto,
  kDeuteronForFemto,
  kProtonForScPrCorr
};

enum V0Species {
  kPhoton = 0,
  kK0S,
  kLambda,
  kAntiLambda,
  kNV0
};

enum HfVtxStage : uint8_t {
  Skimmed = 0,
  BeautyVertex,
  CharmHadPiSelected,
  kNHfVtxStage
};

// Helper struct to pass V0 informations
struct V0Cand {
  std::array<float, 3> mom;
  std::array<float, 3> vtx;
  std::array<float, 21> cov;
  float etaPos;
  float etaNeg;
  float ptPos;
  float ptNeg;
  float pinTpcPos;
  float pinTpcNeg;
  float nClsFoundTpcPos;
  float nClsFoundTpcNeg;
  float nClsCrossedRowsTpcPos;
  float nClsCrossedRowsTpcNeg;
  float crossedRowsOverFindableClsTpcPos;
  float crossedRowsOverFindableClsTpcNeg;
  float signalTpcPos;
  float signalTpcNeg;
  float v0cosPA;
  float dcav0topv;
  float dcaV0daughters;
  float dcapostopv;
  float dcanegtopv;
  float alpha;
  float qtarm;
  float v0radius;
  float mK0Short;
  float mLambda;
  float mAntiLambda;
  float nSigmaPrTpcPos;
  float nSigmaPrTofPos;
  float nSigmaPrTpcNeg;
  float nSigmaPrTofNeg;
  float nSigmaPiTpcPos;
  float nSigmaPiTofPos;
  float nSigmaPiTpcNeg;
  float nSigmaPiTofNeg;
  bool hasTofPos;
  bool hasTofNeg;
};

// Helper struct to pass Cascade informations
struct CascCand {
  std::array<float, 3> mom;
  std::array<float, 3> vtx;
  std::array<float, 21> cov;
  V0Cand v0;
  float ptBach;
  float etaBach;
  float pinTpcBach;
  float nClsFoundTpcBach;
  float nClsCrossedRowsTpcBach;
  float crossedRowsOverFindableClsTpcBach;
  float signalTpcBach;
  float pt;
  float casccosPA;
  float cascradius;
  float dcaXYCascToPV;
  float dcacascdaughters;
  float mXi;
  float mOmega;
  float nSigmaPiTpcBach;
  float nSigmaPiTofBach;
  bool hasTofBach;
  int sign;
};

static const std::array<std::string, kNCharmParticles> charmParticleNames{"D0", "Dplus", "Ds", "Lc", "Xic"};
static const int nTotBeautyParts = static_cast<int>(kNBeautyParticles) + static_cast<int>(kNBeautyParticlesToJPsi);
static const std::array<std::string, nTotBeautyParts> beautyParticleNames{"Bplus", "B0toDStar", "Bc", "B0", "Bs", "Lb", "Xib", "BplusToJPsi", "B0ToJPsi", "BsToJPsi", "LbToJPsi", "BcToJPsi"};
static const std::array<int, kNCharmParticles> pdgCodesCharm{421, 411, 431, 4122, 4232};
static const std::array<std::string, 2> eventTitles = {"all", "rejected"};
static const std::vector<std::string> hfTriggerNames{filtering::HfHighPt2P::columnLabel(), filtering::HfHighPt3P::columnLabel(), filtering::HfBeauty3P::columnLabel(), filtering::HfBeauty4P::columnLabel(), filtering::HfFemto2P::columnLabel(), filtering::HfFemto3P::columnLabel(), filtering::HfDoubleCharm2P::columnLabel(), filtering::HfDoubleCharm3P::columnLabel(), filtering::HfDoubleCharmMix::columnLabel(), filtering::HfV0Charm2P::columnLabel(), filtering::HfV0Charm3P::columnLabel(), filtering::HfCharmBarToXiBach::columnLabel(), filtering::HfSigmaCPPK::columnLabel(), filtering::HfSigmaC0K0::columnLabel(), filtering::HfPhotonCharm2P::columnLabel(), filtering::HfPhotonCharm3P::columnLabel(), filtering::HfSingleCharm2P::columnLabel(), filtering::HfSingleCharm3P::columnLabel(), filtering::HfSingleNonPromptCharm2P::columnLabel(), filtering::HfSingleNonPromptCharm3P::columnLabel(), filtering::HfCharmBarToXi2Bach::columnLabel(), filtering::HfPrCharm2P::columnLabel(), filtering::HfBtoJPsiKa::columnLabel(), filtering::HfBtoJPsiKstar::columnLabel(), filtering::HfBtoJPsiPhi::columnLabel(), filtering::HfBtoJPsiPrKa::columnLabel(), filtering::HfBtoJPsiPi::columnLabel(), filtering::HfSigmaCPr::columnLabel()};

static const std::array<std::string, kNV0> v0Labels{"#gamma", "K_{S}^{0}", "#Lambda", "#bar{#Lambda}"};
static const std::array<std::string, kNV0> v0Names{"Photon", "K0S", "Lambda", "AntiLambda"};

static const std::tuple pdgCharmDaughters{
  std::array{-321, 211},        // D0
  std::array{-321, 211, 211},   // Dplus
  std::array{321, -321, 211},   // Ds
  std::array{2212, -321, 211},  // Lc
  std::array{2212, -321, 211}}; // Xic

constexpr float massPi = o2::constants::physics::MassPiPlus;
constexpr float massKa = o2::constants::physics::MassKPlus;
constexpr float massProton = o2::constants::physics::MassProton;
constexpr float massMu = o2::constants::physics::MassMuon;
constexpr float massDeuteron = o2::constants::physics::MassDeuteron;
constexpr float massGamma = o2::constants::physics::MassGamma;
constexpr float massK0S = o2::constants::physics::MassK0Short;
constexpr float massLambda = o2::constants::physics::MassLambda0;
constexpr float massXi = o2::constants::physics::MassXiMinus;
constexpr float massPhi = o2::constants::physics::MassPhi;
constexpr float massD0 = o2::constants::physics::MassD0;
constexpr float massDPlus = o2::constants::physics::MassDPlus;
constexpr float massDs = o2::constants::physics::MassDS;
constexpr float massLc = o2::constants::physics::MassLambdaCPlus;
constexpr float massXic = o2::constants::physics::MassXiCPlus;
constexpr float massDStar = o2::constants::physics::MassDStar;
constexpr float massBPlus = o2::constants::physics::MassBPlus;
constexpr float massB0 = o2::constants::physics::MassB0;
constexpr float massBs = o2::constants::physics::MassBS;
constexpr float massLb = o2::constants::physics::MassLambdaB0;
constexpr float massXib = o2::constants::physics::MassXiB0;
constexpr float massBc = o2::constants::physics::MassBCPlus;
constexpr float massSigmaCPlusPlus = o2::constants::physics::MassSigmaCPlusPlus;
constexpr float massSigmaC0 = o2::constants::physics::MassSigmaC0;
constexpr float massK0Star892 = o2::constants::physics::MassK0Star892;
constexpr float massJPsi = o2::constants::physics::MassJPsi;

static const o2::framework::AxisSpec ptAxis{50, 0.f, 50.f};
static const o2::framework::AxisSpec pAxis{50, 0.f, 10.f};
static const o2::framework::AxisSpec kstarAxis{200, 0.f, 2.f};
static const o2::framework::AxisSpec etaAxis{30, -1.5f, 1.5f};
static const o2::framework::AxisSpec nSigmaAxis{100, -10.f, 10.f};
static const o2::framework::AxisSpec alphaAxis{100, -1.f, 1.f};
static const o2::framework::AxisSpec qtAxis{100, 0.f, 0.25f};
static const o2::framework::AxisSpec bdtAxis{100, 0.f, 1.f};
static const o2::framework::AxisSpec phiAxis{36, 0., o2::constants::math::TwoPI};
static const std::array<o2::framework::AxisSpec, kNCharmParticles + 24> massAxisC = {o2::framework::AxisSpec{250, 1.65f, 2.15f}, o2::framework::AxisSpec{250, 1.65f, 2.15f}, o2::framework::AxisSpec{250, 1.75f, 2.25f}, o2::framework::AxisSpec{250, 2.05f, 2.55f}, o2::framework::AxisSpec{250, 2.25f, 2.75f}, o2::framework::AxisSpec{200, 0.139f, 0.159f}, o2::framework::AxisSpec{250, 0.f, 0.25f}, o2::framework::AxisSpec{250, 0.f, 0.25f}, o2::framework::AxisSpec{200, 0.48f, 0.88f}, o2::framework::AxisSpec{200, 0.48f, 0.88f}, o2::framework::AxisSpec{200, 1.1f, 1.4f}, o2::framework::AxisSpec{200, 1.1f, 1.4f}, o2::framework::AxisSpec{200, 1.1f, 1.4f}, o2::framework::AxisSpec{200, 1.1f, 1.4f}, o2::framework::AxisSpec{170, 0.13f, 0.3f}, o2::framework::AxisSpec{170, 0.13f, 0.3f}, o2::framework::AxisSpec{200, 0.4f, 0.8f}, o2::framework::AxisSpec{200, 0.4f, 0.8f}, o2::framework::AxisSpec{200, 0.4f, 0.8f}, o2::framework::AxisSpec{200, 0.4f, 0.8f}, o2::framework::AxisSpec{350, 2.3f, 3.0f}, o2::framework::AxisSpec{350, 2.3f, 3.0f}, o2::framework::AxisSpec{350, 2.3f, 3.0f}, o2::framework::AxisSpec{240, 2.4f, 3.6f}, o2::framework::AxisSpec{300, 0.7f, 1.3f}, o2::framework::AxisSpec{300, 0.7f, 1.3f}, o2::framework::AxisSpec{300, 0.7f, 1.3f}, o2::framework::AxisSpec{300, 0.7f, 1.3f}, o2::framework::AxisSpec{300, 0.14f, 0.26f}};
static const std::array<o2::framework::AxisSpec, nTotBeautyParts> massAxisB = {o2::framework::AxisSpec{500, 4.2f, 6.2f}, o2::framework::AxisSpec{500, 4.2f, 6.2f}, o2::framework::AxisSpec{500, 5.4f, 7.4f}, o2::framework::AxisSpec{500, 4.2f, 6.2f}, o2::framework::AxisSpec{500, 4.4f, 6.4f}, o2::framework::AxisSpec{400, 5.0f, 6.6f}, o2::framework::AxisSpec{500, 4.2f, 6.2f}, o2::framework::AxisSpec{500, 4.2f, 6.2f}, o2::framework::AxisSpec{500, 4.2f, 6.2f}, o2::framework::AxisSpec{500, 4.2f, 6.2f}, o2::framework::AxisSpec{400, 5.0f, 6.6f}, o2::framework::AxisSpec{240, 5.8f, 7.0f}};

// default values for configurables
// channels to trigger on for femto
constexpr int activeFemtoChannels[2][5] = {{1, 1, 1, 1, 0},  // pD0, pD+, pDs, pLc, pXic
                                           {0, 0, 0, 1, 0}}; // only for deLc
static const std::vector<std::string> labelsColumnsFemtoChannels = {"DZero", "DPlus", "Ds", "Lc", "Xic"};
static const std::vector<std::string> labelsRowsFemtoChannels = {"protonCharmFemto", "deuteronCharmFemto"};
constexpr float cutsPtThresholdsForFemto[1][2] = {{8., 1.4}}; // proton, deuteron
static const std::vector<std::string> labelsColumnsPtThresholdsForFemto = {"Proton", "Deuteron"};

// min and max pT for all tracks combined  (except for V0 and cascades)
constexpr float cutsPt[2][10] = {{1., 0.1, 0.8, 0.5, 0.1, 0.2, 0.4, 0.5, 0.3, 0.3},
                                 {100000., 100000., 5., 100000., 100000., 100000., 100000., 100000., 100000., 100000.}}; // beauty, D*, femto, SigmaC, Xic*+ -> SigmaC++K-, beauty to JPsi, Lc*->D0p
static const std::vector<std::string> labelsColumnsCutsPt = {"Beauty", "DstarPlus", "PrForFemto", "CharmBaryon", "SoftPiSigmaC", "SoftKaonXicResoToSigmaC", "DeForFemto", "BeautyToJPsi", "PrForLcReso", "PrForThetaC"};
static const std::vector<std::string> labelsRowsCutsPt = {"Minimum", "Maximum"};

// PID cuts
constexpr float cutsNsigma[4][9] = {
  {3., 3., 3., 5., 3., 3., 5., 3., 3.},               // TPC proton from Lc, pi/K from D0, K from 3-prong, femto selected proton, pi/K from Xic/Omegac, K from Xic*->SigmaC-Kaon, femto selected deuteron, K/p from beauty->JPsiX, proton from SigmaC-Pr correaltion
  {3., 3., 3., 2.5, 3., 3., 5., 3., 3.},              // TOF proton from Lc, pi/K from D0, K from 3-prong, femto selected proton, pi/K from Xic/Omegac, K from Xic*->SigmaC-Kaon, femto selected deuteron, K/p from beauty->JPsiX, proton from SigmaC-Pr correaltion
  {999., 999., 999., 2.5, 999., 999., 5., 999., 3.},  // Sum in quadrature of TPC and TOF (used only for femto selected proton and deuteron for pT < 4 GeV/c)
  {999., 999., 999., 999., 999., 999., -4., 999., 999.} // ITS used only for femto selected deuteron for less than pt threshold
};
static const std::vector<std::string> labelsColumnsNsigma = {"PrFromLc", "PiKaFromDZero", "KaFrom3Prong", "PrForFemto", "PiKaFromCharmBaryon", "SoftKaonFromXicResoToSigmaC", "DeForFemto", "KaPrFromBeautyToJPsi", "PrFromSigmaCPr"};
static const std::vector<std::string> labelsRowsNsigma = {"TPC", "TOF", "Comb", "ITS"};

// track cut
constexpr float cutsTrackQuality[2][7] = {{0., 0., 0., 999., 999., 0., 0.},
                                          {90, 80, 0.83, 160., 1., 5., 0.}};
static const std::vector<std::string> labelsColumnsTrackQuality = {"minTpcCluster", "minTpcRow", "minTpcCrossedOverFound", "maxTpcShared", "maxTpcFracShared", "minItsCluster", "minItsIbCluster"};

// high pt
constexpr float cutsHighPtThresholds[1][2] = {{8., 8.}}; // 2-prongs, 3-prongs
static const std::vector<std::string> labelsColumnsHighPtThresholds = {"2Prongs", "3Prongs"};

namespace hf_trigger_cuts_presel_beauty
{
static constexpr int nBinsPt = 2;
static constexpr int nCutVars = 4;
static constexpr int nCutVarsBtoJPsi = 6;
// default values for the pT bin edges (can be used to configure histogram axis)
// common for any beauty candidate
constexpr double binsPt[nBinsPt + 1] = {
  0.,
  5.,
  1000.0};
auto vecBinsPt = std::vector<double>{binsPt, binsPt + nBinsPt + 1};
// default values for the cuts
constexpr double cuts[nBinsPt][nCutVars] = {{0.4, -1, -1, 10.},  /* 0 < pt < 5 */
                                            {0.4, -1, -1, 10.}}; /* 5 < pt < 1000 */

constexpr double cutsBtoJPsi[nBinsPt][nCutVarsBtoJPsi] = {{1., 0.6, 0.9, 0.02, 0.02, 0.1},  /* 0 < pt < 5 */
                                                          {1., 0.8, 0.9, 0.02, 0.02, 0.1}}; /* 5 < pt < 1000 */

// row labels
static const std::vector<std::string> labelsPt{};
// column labels
static const std::vector<std::string> labelsColumnsTopolBeauty = {"DeltaMassB", "minCPA", "minDecayLength", "maxImpParProd"};
static const std::vector<std::string> labelsColumnsCutsBeautyToJPsi = {"minPtMuon", "DeltaMassB", "minCPA", "minDecayLength", "DeltaMassKK", "DeltaMassKPi"};

} // namespace hf_trigger_cuts_presel_beauty

// double charm
constexpr int activeDoubleCharmChannels[2][3] = {{1, 1, 1}, {1, 1, 0}}; // kDoubleCharm2P, kDoubleCharm3P, kDoubleCharmMix (second column to keep non-prompt)
static const std::vector<std::string> labelsColumnsDoubleCharmChannels = {"DoubleCharm2Prong", "DoubleCharm3Prong", "DoubleCharmMix"};
static const std::vector<std::string> labelsRowsDoubleCharmChannels = {"", "KeepNonprompt"};

// charm resonances
constexpr float cutsCharmReso[4][14] = {{0.0, 0.0, 0.0, 0.0, 0.4, 0., 0.0, 0.00, 0.21, 0.21, 0.0, 0.7, 0.7, 0.15},
                                        {0.155, 0.3, 0.3, 0.88, 0.88, 1.35, 0.18, 0.18, 0.25, 0.25, 0.8, 1.3, 1.3, 0.19},
                                        {0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 6.0, 0.0, 6.0, 0.0, 0.0, 0.0, 5.0},
                                        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}}; // D*+, D*0, Ds*0, Ds1+, Ds2*+, Xic*->D, SigmaC0, SigmaC++, SigmaC(2520)0, SigmaC(2520)++, Xic*->SigmaC, Lc*->D0P, Lc*->D*+P
static const std::vector<std::string> labelsColumnsDeltaMassCharmReso = {"DstarPlus", "DstarZero", "DsStarZero", "Ds1Plus", "Ds2StarPlus", "XicResoToD", "SigmaC0", "SigmaCPlusPlus", "SigmaC02520", "SigmaCPlusPlus2520", "XicResoToSigmaC", "LcResoToD0Pr", "ThetaC", "SigmaCPr"};
static const std::vector<std::string> labelsRowsDeltaMassCharmReso = {"deltaMassMin", "deltaMassMax", "ptMin", "ptMinCharmDaugh"};
// V0s for charm resonances
constexpr float cutsV0s[1][6] = {{0.85, 0.97, 0.5, 4., 0.02, 0.01}}; // cosPaGamma, cosPaK0sLambda, radiusK0sLambda, nSigmaPrLambda, deltaMassK0S, deltaMassLambda
static const std::vector<std::string> labelsColumnsV0s = {"CosPaGamma", "CosPaK0sLambda", "RadiusK0sLambda", "NSigmaPrLambda", "DeltaMassK0s", "DeltaMassLambda"};

// cascades for Xi + bachelor triggers
constexpr float cutsCascades[1][8] = {{0.2, 1., 0.01, 0.01, 0.99, 0.99, 0.3, 3.}}; // ptXiBachelor, deltaMassXi, deltaMassLambda, cosPaXi, cosPaLambda, DCAxyXi, nSigmaPid
static const std::vector<std::string> labelsColumnsCascades = {"PtBachelor", "PtXi", "DeltaMassXi", "DeltaMassLambda", "CosPAXi", "CosPaLambda", "DCAxyXi", "NsigmaPid"};
constexpr float cutsCharmBaryons[1][15] = {{5., 5., 1000., 2.35, 2.60, 2.35, 3., 3., 2.7, -2., -2., 1.e6, 1.e6, -1., -1.}}; // MinPtXiPi, MinPtXiKa, MinPtXiPiPi, MinMassXiPi, MinMassXiKa, MinMassXiPiPi, MaxMassXiPi, MaxMassXiKa, MaxMassXiPiPi, CosPaXiBach, CosPaXiBachBach, Chi2PcaXiBach, Chi2PcaXiBachBach, DecLenXiBach, DecLenBachBach
static const std::vector<std::string> labelsColumnsCharmBarCuts = {"MinPtXiPi", "MinPtXiKa", "MinPtXiPiPi", "MinMassXiPi", "MinMassXiKa", "MinMassXiPiPi", "MaxMassXiPi", "MaxMassXiKa", "MaxMassXiPiPi", "CosPaXiBach", "CosPaXiBachBach", "Chi2PcaXiBach", "Chi2PcaXiBachBach", "DecLenXiBach", "DecLenBachBach"};

//proton for SigmaC-pr trigger
constexpr float cutsSigmaCPrDefault[3][1] = {
  {0.399},   // ptPrMin
  {4.501},   // ptPrMax
  {1.0}   // ptTOFThreshold
};
static const std::vector<std::string> labelsRowsSigmaCPr = {
  "ptPrMin",
  "ptPrMax",
  "ptTOFThreshold"
};
static const std::vector<std::string> labelsColumnsSigmaCPr = {
  "SigmaCPr"
};

constexpr int requireStrangenessTrackedXi[1][2] = {{1, 0}};
static const std::vector<std::string> labelsColumnsCharmBaryons = {"CharmBarToXiBach", "CharmBarToXiBachBach"};

// dummy array
static const std::vector<std::string> labelsEmpty{};
static constexpr double cutsTrackDummy[o2::analysis::hf_cuts_single_track::NBinsPtTrack][o2::analysis::hf_cuts_single_track::NCutVarsTrack] = {{0., 10.}, {0., 10.}, {0., 10.}, {0., 10.}, {0., 10.}, {0., 10.}};
o2::framework::LabeledArray<double> cutsSingleTrackDummy{cutsTrackDummy[0], o2::analysis::hf_cuts_single_track::NBinsPtTrack, o2::analysis::hf_cuts_single_track::NCutVarsTrack, o2::analysis::hf_cuts_single_track::labelsPtTrack, o2::analysis::hf_cuts_single_track::labelsCutVarTrack};

// manual downscale factors for tests
constexpr double defDownscaleFactors[kNtriggersHF][1] = {{1.1}, {1.1}, {1.1}, {1.1}, {1.1}, {1.1}, {1.1}, {1.1}, {1.1}, {1.1}, {1.1}, {1.1}, {1.1}, {1.1}, {1.1}, {1.1}, {1.1}, {1.1}, {1.1}, {1.1}, {1.1}, {1.1}, {1.1}, {1.1}, {1.1}, {1.1}, {1.1}}; // one for each trigger
static const std::vector<std::string> labelsDownscaleFactor = {"Downscale factor"};

// Main helper class

class HfFilterHelper
{
 public:
  /// Default constructor
  HfFilterHelper() = default;

  // setters
  void setHighPtTriggerThresholds(float threshold2Prongs, float threshold3Prongs)
  {
    mPtThresholdHighPt2Prongs = threshold2Prongs;
    mPtThresholdHighPt3Prongs = threshold3Prongs;
  }
  void setPtTriggerThresholdsForFemto(float thresholdProtons, float thresholdDeuterons)
  {
    mPtThresholdProtonForFemto = thresholdProtons;
    mPtThresholdDeuteronForFemto = thresholdDeuterons;
  }
  void setForceTofForFemto(bool forceTofProtons, bool forceTofDeuterons)
  {
    mForceTofProtonForFemto = forceTofProtons;
    mForceTofDeuteronForFemto = forceTofDeuterons;
  }
  void setPtBinsSingleTracks(const std::vector<double>& ptBins) { mPtBinsTracks = ptBins; }
  void setPtBinsBeautyHadrons(const std::vector<double>& ptBins) { mPtBinsBeautyHadrons = ptBins; }
  void setCutsSingleTrackBeauty(const o2::framework::LabeledArray<double>& cutsSingleTrack3P, const o2::framework::LabeledArray<double>& cutsSingleTrack4P, const o2::framework::LabeledArray<double>& cutsSingleToJPsi)
  {
    mCutsSingleTrackBeauty3Prong = cutsSingleTrack3P;
    mCutsSingleTrackBeauty4Prong = cutsSingleTrack4P;
    mCutsSingleTrackBeautyToJPsi = cutsSingleToJPsi;
  }
  void setCutsBhadrons(const o2::framework::LabeledArray<double>& cutsBplus, const o2::framework::LabeledArray<double>& cutsB0toDstar, const o2::framework::LabeledArray<double>& cutsBc, const o2::framework::LabeledArray<double>& cutsB0, const o2::framework::LabeledArray<double>& cutsBs, const o2::framework::LabeledArray<double>& cutsLb, const o2::framework::LabeledArray<double>& cutsXib)
  {
    mCutsBhad[kBplus] = cutsBplus;
    mCutsBhad[kB0toDStar] = cutsB0toDstar;
    mCutsBhad[kBc] = cutsBc;
    mCutsBhad[kB0] = cutsB0;
    mCutsBhad[kBs] = cutsBs;
    mCutsBhad[kLb] = cutsLb;
    mCutsBhad[kXib] = cutsXib;
  }
  void setCutsBtoJPsi(const o2::framework::LabeledArray<double>& cuts)
  {
    mCutsBhadToJPsi = cuts;
  }
  void setPtLimitsProtonForFemto(float minPt, float maxPt)
  {
    mPtMinProtonForFemto = minPt;
    mPtMaxProtonForFemto = maxPt;
  }
  void setPtLimitsDeuteronForFemto(float minPt, float maxPt)
  {
    mPtMinDeuteronForFemto = minPt;
    mPtMaxDeuteronForFemto = maxPt;
  }
  void setPtLimitsBeautyBachelor(float minPt, float maxPt, float minPtBtoJPsiBach, float maxPtBtoJPsiBach)
  {
    mPtMinBeautyBachelor = minPt;
    mPtMaxBeautyBachelor = maxPt;
    mPtMinBeautyToJPsiBachelor = minPtBtoJPsiBach;
    mPtMaxBeautyToJPsiBachelor = maxPtBtoJPsiBach;
  }
  void setPtLimitsDstarSoftPion(float minPt, float maxPt)
  {
    mPtMinSoftPionForDstar = minPt;
    mPtMaxSoftPionForDstar = maxPt;
  }
  void setPtRangeSoftPiSigmaC(float minPt, float maxPt)
  {
    mPtMinSoftPionForSigmaC = minPt;
    mPtMaxSoftPionForSigmaC = maxPt;
  }
  void setPtDeltaMassRangeSigmaC(float minDeltaMassSigmaCZero, float maxDeltaMassSigmaCZero, float minDeltaMassSigmaCPlusPlus, float maxDeltaMassSigmaCPlusPlus, float minDeltaMassSigmaC2520Zero, float maxDeltaMassSigmaC2520Zero, float minDeltaMassSigmaC2520PlusPlus, float maxDeltaMassSigmaC2520PlusPlus, float minPtSigmaCZero, float minPtSigmaCPlusPlus, float minPtSigmaC2520Zero, float minPtSigmaC2520PlusPlus)
  {
    mDeltaMassMinSigmaCZero = minDeltaMassSigmaCZero;
    mDeltaMassMaxSigmaCZero = maxDeltaMassSigmaCZero;
    mDeltaMassMinSigmaC2520Zero = minDeltaMassSigmaC2520Zero;
    mDeltaMassMaxSigmaC2520Zero = maxDeltaMassSigmaC2520Zero;
    mDeltaMassMinSigmaCPlusPlus = minDeltaMassSigmaCPlusPlus;
    mDeltaMassMaxSigmaCPlusPlus = maxDeltaMassSigmaCPlusPlus;
    mDeltaMassMinSigmaC2520PlusPlus = minDeltaMassSigmaC2520PlusPlus;
    mDeltaMassMaxSigmaC2520PlusPlus = maxDeltaMassSigmaC2520PlusPlus;
    mPtMinSigmaCZero = minPtSigmaCZero;
    mPtMinSigmaC2520Zero = minPtSigmaC2520Zero;
    mPtMinSigmaCPlusPlus = minPtSigmaCPlusPlus;
    mPtMinSigmaC2520PlusPlus = minPtSigmaC2520PlusPlus;
  }
    void setParSigmaCPr(float minDeltaMassSigmaC, float maxDeltaMassSigmaC, float minPtSigmaC, float minPtProton, float maxPtProton, float minPtForTOF, bool forceTOF)
  {
    mMinDeltaMassScSigmaCPr = minDeltaMassSigmaC;
    mMaxDeltaMassScSigmaCPr = maxDeltaMassSigmaC;

    mMinPtScSigmaPr = minPtSigmaC;

    mMinPtPrSigmaCPr = minPtProton;
    mMaxPtPrSigmaCPr = maxPtProton;

    mForceTOFForPrSigmaCPr = forceTOF;
    mThresholdPtTOFForPrSigmaCPr = minPtForTOF;
  }
  void setPtRangeSoftKaonXicResoToSigmaC(float minPt, float maxPt)
  {
    mPtMinSoftKaonForXicResoToSigmaC = minPt;
    mPtMaxSoftKaonForXicResoToSigmaC = maxPt;
  }
  void setPtLimitsCharmBaryonBachelor(float minPt, float maxPt)
  {
    mPtMinCharmBaryonBachelor = minPt;
    mPtMaxCharmBaryonBachelor = maxPt;
  }
  void setPtLimitsLcResonanceBachelor(float minPt, float maxPt)
  {
    mPtMinLcResonanceBachelor = minPt;
    mPtMaxLcResonanceBachelor = maxPt;
  }
  void setPtLimitsThetaCBachelor(float minPt, float maxPt)
  {
    mPtMinThetaCBachelor = minPt;
    mPtMaxThetaCBachelor = maxPt;
  }

  void setNsigmaProtonCutsForFemto(std::array<float, 4> nSigmaCuts) { mNSigmaPrCutsForFemto = nSigmaCuts; }
  void setNsigmaDeuteronCutsForFemto(std::array<float, 4> nSigmaCuts) { mNSigmaDeCutsForFemto = nSigmaCuts;}
  void setNsigmaProtonCutsForSigmaCPr(std::array<float, 4> nSigmaCuts) { mNSigmaPrCutsForSigmaCPr = nSigmaCuts; }

  void setDeuteronTrackSelectionForFemto(float minTpcCluster, float minTpcRow, float minTpcCrossedOverFound, float maxTpcShared, float maxTpcFracShared, float minItsCluster, float minItsIbCluster)
  {
    mMinTpcCluster = minTpcCluster;
    mMinTpcRow = minTpcRow;
    mMinTpcCrossedOverFound = minTpcCrossedOverFound;
    mMaxTpcShared = maxTpcShared;
    mMaxTpcFracShared = maxTpcFracShared;
    mMinItsCluster = minItsCluster;
    mMinItsIbCluster = minItsIbCluster;
  }

  void setNsigmaProtonCutsForCharmBaryons(float nSigmaTpc, float nSigmaTof)
  {
    mNSigmaTpcPrCutForCharmBaryons = nSigmaTpc;
    mNSigmaTofPrCutForCharmBaryons = nSigmaTof;
  }
  void setNsigmaKaonCutsFor3Prongs(float nSigmaTpc, float nSigmaTof)
  {
    mNSigmaTpcKaCutFor3Prongs = nSigmaTpc;
    mNSigmaTofKaCutFor3Prongs = nSigmaTof;
  }
  void setNsigmaPionKaonCutsForDzero(float nSigmaTpc, float nSigmaTof)
  {
    mNSigmaTpcPiKaCutForDzero = nSigmaTpc;
    mNSigmaTofPiKaCutForDzero = nSigmaTof;
  }
  void setNsigmaKaonProtonCutsForBeautyToJPsi(float nSigmaTpc, float nSigmaTof)
  {
    mNSigmaTpcPrKaCutForBeautyToJPsi = nSigmaTpc;
    mNSigmaTofPrKaCutForBeautyToJPsi = nSigmaTof;
  }
  void setV0Selections(float minGammaCosPa, float minK0sLambdaCosPa, float minK0sLambdaRadius, float nSigmaPrFromLambda, float deltaMassK0s, float deltaMassLambda)
  {
    mMinGammaCosinePa = minGammaCosPa;
    mMinK0sLambdaCosinePa = minK0sLambdaCosPa;
    mMinK0sLambdaRadius = minK0sLambdaRadius;
    mMaxNsigmaPrForLambda = nSigmaPrFromLambda;
    mDeltaMassK0s = deltaMassK0s;
    mDeltaMassLambda = deltaMassLambda;
  }
  void setXiSelections(float minPtXiBachelor, float minPtXi, float deltaMassXi, float deltaMassLambda, float cosPaXi, float cosPaLambdaFromXi, float maxDcaxyXi, float nSigma)
  {
    mMinPtXiBachelor = minPtXiBachelor;
    mMinPtXi = minPtXi;
    mDeltaMassXi = deltaMassXi;
    mDeltaMassLambdaFromXi = deltaMassLambda;
    mCosPaXi = cosPaXi;
    mCosPaLambdaFromXi = cosPaLambdaFromXi;
    mMaxDcaXyXi = maxDcaxyXi;
    mMaxNsigmaXiDau = nSigma;
  }
  void setCutsSingleTrackCharmBaryonBachelor(const o2::framework::LabeledArray<double>& cutsSingleTrack) { mCutsSingleTrackCharmBaryonBachelor = cutsSingleTrack; }
  void setNsigmaPiCutsForCharmBaryonBachelor(float nSigmaTpc, float nSigmaTof)
  {
    mNSigmaTpcPiCharmBaryonBachelor = nSigmaTpc;
    mNSigmaTofPiCharmBaryonBachelor = nSigmaTof;
  }
  void setNsigmaTpcKaonFromXicResoToSigmaC(float nSigmaTpc, float nSigmaTof)
  {
    mNSigmaTpcKaonFromXicResoToSigmaC = nSigmaTpc;
    mNSigmaTofKaonFromXicResoToSigmaC = nSigmaTof;
  }

  void setXiBachelorSelections(float ptMinXiPi, float ptMinXiKa, float ptMinXiPiPi, float massMinXiPi, float massMinXiKa, float massMinXiPiPi, float massMaxXiPi, float massMaxXiKa, float massMaxXiPiPi, float cosPaMinXiBach, float cosPaMinXiBachBach, float chi2PcaMaxXiBach, float chi2PcaMaxXiBachBach, float decLenMinXiBach, float decLenMinXiBachBach)
  {
    mPtMinXiBach[0] = ptMinXiPi;
    mPtMinXiBach[1] = ptMinXiKa;
    mPtMinXiBach[2] = ptMinXiPiPi;
    mMassMinXiBach[0] = massMinXiPi;
    mMassMinXiBach[1] = massMinXiKa;
    mMassMinXiBach[2] = massMinXiPiPi;
    mMassMaxXiBach[0] = massMaxXiPi;
    mMassMaxXiBach[1] = massMaxXiKa;
    mMassMaxXiBach[2] = massMaxXiPiPi;
    mCosPaMinXiBach[0] = cosPaMinXiBach;
    mCosPaMinXiBach[1] = cosPaMinXiBachBach;
    mChi2PcaMaxXiBach[0] = chi2PcaMaxXiBach;
    mChi2PcaMaxXiBach[1] = chi2PcaMaxXiBachBach;
    mDecLenMinXiBach[0] = decLenMinXiBach;
    mDecLenMinXiBach[1] = decLenMinXiBachBach;
  }

  void setTpcPidCalibrationOption(int opt) { mTpcPidCalibrationOption = opt; }

  void setMassResolParametrisation(const std::string& recoPass)
  {
    if (recoPass == "2023_pass3") {
      mSigmaPars2Prongs[0] = 0.01424f;
      mSigmaPars2Prongs[1] = 0.00178f;
      mDeltaMassPars2Prongs[0] = -0.0025f;
      mDeltaMassPars2Prongs[1] = 0.0001f;
      mSigmaPars3Prongs[0] = 0.00796f;
      mSigmaPars3Prongs[1] = 0.00176f;
      mDeltaMassPars3Prongs[0] = -0.0025f;
      mDeltaMassPars3Prongs[1] = 0.0001f;
    } else if (recoPass == "2025_pass1") {
      mSigmaPars2Prongs[0] = 0.01424f;
      mSigmaPars2Prongs[1] = 0.00178f;
      mDeltaMassPars2Prongs[0] = -0.013f;
      mDeltaMassPars2Prongs[1] = 0.00029f;
      mSigmaPars3Prongs[0] = 0.00796f;
      mSigmaPars3Prongs[1] = 0.00176f;
      mDeltaMassPars3Prongs[0] = -0.013f;
      mDeltaMassPars3Prongs[1] = 0.00029f;
    } else {
      LOGP(fatal, "Mass resolution parametrisation {} not supported! Please set 2023_pass3", recoPass.data());
    }
  }

  void setNumSigmaForDeltaMassCharmHadCut(float nSigma) { mNumSigmaDeltaMassCharmHad = nSigma; }

  void setPreselDsToKKPi(const std::vector<double>& ptBins, const o2::framework::LabeledArray<double>& preselections)
  {
    mPtBinsPreselDsToKKPi = ptBins;
    mPreselDsToKKPi = preselections;
  }

  // helper functions for selections
  template <typename T>
  bool isSelectedHighPt2Prong(const T& pt);
  template <typename T>
  bool isSelectedHighPt3Prong(const T& pt);
  template <o2::aod::hffilters::HfTriggers whichTrigger, typename T, typename T1, typename T2>
  int16_t isSelectedTrackForSoftPionOrBeauty(const T& track, const T1& trackPar, const T2& dca);
  template <typename T1, typename T2, typename H2>
  bool isSelectedTrack4Corr(const T1& track, const T2& trackPar, const int& activateQA, H2 hTPCPID, H2 hTOFPID, const int& trackSpecies);
  template <typename T>
  int8_t isDzeroPreselected(const T& trackPos, const T& trackNeg);
  template <typename T>
  int8_t isDplusPreselected(const T& trackOppositeCharge);
  template <typename P, typename T>
  int8_t isDsPreselected(const P& pTrackSameChargeFirst, const P& pTrackSameChargeSecond, const P& pTrackOppositeCharge, const T& trackOppositeCharge);
  template <typename T>
  int8_t isCharmBaryonPreselected(const T& trackSameChargeFirst, const T& trackSameChargeSecond, const T& trackOppositeCharge);
  template <typename T, typename H2>
  int8_t isSelectedD0InMassRange(const T& pTrackPos, const T& pTrackNeg, const float& ptD, int8_t isSelected, const int& activateQA, H2 hMassVsPt);
  template <typename T, typename H2>
  int8_t isSelectedDplusInMassRange(const T& pTrackSameChargeFirst, const T& pTrackSameChargeSecond, const T& pTrackOppositeCharge, const float& ptD, const int& activateQA, H2 hMassVsPt);
  template <typename T, typename H2>
  int8_t isSelectedDsInMassRange(const T& pTrackSameChargeFirst, const T& pTrackSameChargeSecond, const T& pTrackOppositeCharge, const float& ptD, int8_t isSelected, const int& activateQA, H2 hMassVsPt);
  template <typename T, typename H2>
  int8_t isSelectedLcInMassRange(const T& pTrackSameChargeFirst, const T& pTrackSameChargeSecond, const T& pTrackOppositeCharge, const float& ptLc, const int8_t isSelected, const int& activateQA, H2 hMassVsPt);
  template <int charge, typename T, typename H2>
  int8_t isSelectedSigmaCInDeltaMassRange(const T& pTrackSameChargeFirst, const T& pTrackSameChargeSecond, const T& pTrackOppositeCharge, const T& pTrackSoftPi, const float ptSigmaC, const int8_t isSelectedLc, H2 hMassVsPt, const int& activateQA);
  template <typename T, typename H2>
  int8_t isSelectedXicInMassRange(const T& pTrackSameChargeFirst, const T& pTrackSameChargeSecond, const T& pTrackOppositeCharge, const float& ptXic, const int8_t isSelected, const int& activateQA, H2 hMassVsPt);
  template <typename V0, typename H2>
  int8_t isSelectedV0(const V0& v0, const int& activateQA, H2 hV0Selected, std::array<H2, 4>& hArmPod);
  template <typename Photon, typename T, typename H2>
  bool isSelectedPhoton(const Photon& photon, const std::array<T, 2>& dauTracks, const int& activateQA, H2 hV0Selected, std::array<H2, 4>& hArmPod);
  template <typename Casc>
  bool isSelectedCascade(const Casc& casc);
  template <typename T, typename T2>
  int16_t isSelectedBachelorForCharmBaryon(const T& track, const T2& dca);
  template <bool is4beauty = false, typename T>
  bool isSelectedProton4CharmOrBeautyBaryons(const T& track);
  template <typename T, typename U>
  int8_t isBDTSelected(const T& scores, const U& thresholdBDTScores);
  template <bool isKaonTrack, typename T>
  bool isSelectedKaonFromXicResoToSigmaC(const T& track);
  template <typename T1, typename T2, typename T3, typename T4>
  bool isSelectedBhadron(T1 const& pVecTrack0, T1 const& pVecTrack1, T2 const& dcaTrack0, T2 const& dcaTrack1, const T3& primVtx, const T4& secVtx, const int whichB);
  template <typename T1, typename T2>
  bool isSelectedBhadronInMassRange(T1 const& ptCand, T2 const& massCand, const int whichB);
  template <typename T1, typename T2, typename T3>
  bool isSelectedBzeroToDstar(T1 const& pVecTrack0, T1 const& pVecTrack1, T1 const& pVecTrack2, const T2& primVtx, const T3& secVtx);
  template <int Nprongs, typename T1, typename T2, typename T3, typename T4, typename H2>
  int8_t isSelectedBhadronToJPsi(std::array<T1, Nprongs> pVecDauTracks, std::array<T2, Nprongs - 2> tracksDauNoMu, const T3& primVtx, const T4& secVtx, const int& activateQA, std::array<H2, nTotBeautyParts>& hMassVsPt);
  template <typename T1>
  bool isCharmHadronMassInSbRegions(T1 const& massHypo1, T1 const& massHypo2, const float& lowLimitSB, const float& upLimitSB);
  template <typename T, typename C, typename H2>
  bool isSelectedXiBach(T const& trackParCasc, T const& trackParBachelor, int8_t isSelBachelor, C const& collision, o2::vertexing::DCAFitterN<2>& dcaFitter, const int& activateQA, H2 hMassVsPtXiPi, H2 hMassVsPtXiKa);
  template <int Nprongs, typename T, typename C, typename H2>
  bool isSelectedXiBachBach(T const& trackParCasc, std::array<T, 2> const& trackParBachelor, C const& collision, o2::vertexing::DCAFitterN<Nprongs>& dcaFitter, const int& activateQA, H2 hMassVsPtXiPiPi);
  template <bool is4ThetaC = false, typename T>
  bool isSelectedProtonFromLcResoOrThetaC(const T& track);
  // helpers
  template <typename T>
  T computeRelativeMomentum(const std::array<T, 3>& pTrack, const std::array<T, 3>& CharmCandMomentum, const T& CharmMass);
  template <typename T>
  int computeNumberOfCandidates(const std::vector<std::vector<T>>& indices);
  template <typename T1>
  int setVtxConfiguration(T1& vertexer, bool useAbsDCA);
  template <typename V, typename T, typename C>
  bool buildV0(V const& v0Indices, T const& tracks, C const& collision, o2::vertexing::DCAFitterN<2>& dcaFitter, const std::vector<int>& vetoedTrackIds, V0Cand& v0Cand);
  template <typename Casc, typename T, typename C, typename V>
  bool buildCascade(Casc const& cascIndices, V const& v0Indices, T const& tracks, C const& collision, o2::vertexing::DCAFitterN<2>& dcaFitter, const std::vector<int>& vetoedTrackIds, CascCand& cascCand);

  // PID
  void setValuesBB(o2::ccdb::CcdbApi& ccdbApi, aod::BCsWithTimestamps::iterator const& bunchCrossing, const std::array<std::string, 8>& ccdbPaths);
  void setTpcRecalibMaps(o2::framework::Service<o2::ccdb::BasicCCDBManager> const& ccdb, aod::BCsWithTimestamps::iterator const& bunchCrossing, const std::string& ccdbPath);

 private:
  // selections
  template <bool is4beauty = false, typename T>
  bool isSelectedKaon4Charm3ProngOrBeautyToJPsi(const T& track);

  // PID
  float getTPCSplineCalib(const float tpcPin, const float dEdx, const int& pidSpecies);
  template <typename T>
  float getTPCSplineCalib(const T& track, const int& pidSpecies);
  float getTPCPostCalib(const float tpcPin, const float tpcNCls, const float eta, const float tpcNSigma, const int& pidSpecies);
  template <typename T>
  float getTPCPostCalib(const T& track, const int& pidSpecies);

  // helpers
  template <typename T1, typename T2>
  int findBin(T1 const& binsPt, T2 value);
  template <typename T>
  std::array<T, 2> alphaAndQtAP(std::array<T, 3> const& momPos, std::array<T, 3> const& momNeg);

  // selections
  std::vector<double> mPtBinsTracks{};                                            // vector of pT bins for single track cuts
  std::vector<double> mPtBinsBeautyHadrons{};                                     // vector of pT bins for beauty hadron candidates
  o2::framework::LabeledArray<double> mCutsSingleTrackBeauty3Prong{};             // dca selections for the 3-prong b-hadron pion daughter
  o2::framework::LabeledArray<double> mCutsSingleTrackBeauty4Prong{};             // dca selections for the 4-prong b-hadron pion daughter
  o2::framework::LabeledArray<double> mCutsSingleTrackBeautyToJPsi{};             // dca selections for the b-hadron -> JPsi X daughters (not the muons)
  float mPtMinSoftPionForDstar{0.1};                                              // minimum pt for the D*+ soft pion
  float mPtMinSoftPionForSigmaC{0.1};                                             // minimum pt for the Σ0,++ soft pion
  float mPtMaxSoftPionForSigmaC{10000.f};                                         // maximum pt for the Σ0,++ soft pion
  float mPtMinSoftKaonForXicResoToSigmaC{0.1};                                    // minimum pt for the soft kaon of Xic* to SigmaC-Kaon
  float mPtMaxSoftKaonForXicResoToSigmaC{10000.f};                                // maximum pt for the soft kaon of Xic* to SigmaC-Kaon
  float mPtMinBeautyBachelor{0.5};                                                // minimum pt for the b-hadron pion daughter
  float mPtMinBeautyToJPsiBachelor{0.5};                                          // minimum pt for the b-hadron -> JPsi X daughters (not the muons)
  float mPtMinProtonForFemto{0.8};                                                // minimum pt for the proton for femto
  float mPtMinDeuteronForFemto{0.8};                                              // minimum pt for the deuteron for femto
  float mPtMinCharmBaryonBachelor{0.5};                                           // minimum pt for the bachelor pion from Xic/Omegac decays
  float mPtMinLcResonanceBachelor{0.3};                                           // minimum pt for the bachelor proton from Lc resonance decays
  float mPtMinThetaCBachelor{0.3};                                                // minimum pt for the bachelor proton from ThetaC decays
  float mPtMaxSoftPionForDstar{2.};                                               // maximum pt for the D*+ soft pion
  float mPtMaxBeautyBachelor{100000.};                                            // maximum pt for the b-hadron pion daughter
  float mPtMaxBeautyToJPsiBachelor{100000.};                                      // maximum pt for the b-hadron -> JPsi X daughters (not the muons)
  float mPtMaxProtonForFemto{5.0};                                                // maximum pt for the proton for femto
  float mPtMaxDeuteronForFemto{5.0};                                              // maximum pt for the deuteron for femto
  float mPtMaxCharmBaryonBachelor{100000.};                                       // maximum pt for the bachelor pion from Xic/Omegac decays
  float mPtMaxLcResonanceBachelor{100000.};                                       // maximum pt for the bachelor proton from Lc resonance decays
  float mPtMaxThetaCBachelor{100000.};                                            // maximum pt for the bachelor proton from ThetaC decays
  float mPtThresholdProtonForFemto{8.};                                           // pt threshold to change strategy for proton PID for femto
  float mPtThresholdDeuteronForFemto{1.4};                                        // pt threshold to change strategy for deuteron PID for femto
  float mPtMinSigmaCZero{0.f};                                                    // pt min SigmaC0 candidate
  float mPtMinSigmaC2520Zero{0.f};                                                // pt min SigmaC(2520)0 candidate
  float mPtMinSigmaCPlusPlus{0.f};                                                // pt min SigmaC++ candidate
  float mPtMinSigmaC2520PlusPlus{0.f};                                            // pt min SigmaC(2520)++ candidate
  std::array<float, 4> mNSigmaPrCutsForFemto{3., 3., 3., -4.};                    // cut values for Nsigma TPC, TOF, combined, ITS for femto protons
  std::array<float, 4> mNSigmaDeCutsForFemto{3., 3., 3., -4.};                    // cut values for Nsigma TPC, TOF, combined, ITS for femto deuterons
  std::array<float, 4> mNSigmaPrCutsForSigmaCPr{3., 3., 3., -4.};                 // cut values for Nsigma TPC, TOF, combined, ITS for proton in Sc-p correaltion
  float mNSigmaTpcPrCutForCharmBaryons{3.};                                       // maximum Nsigma TPC for protons in Lc and Xic decays
  float mNSigmaTofPrCutForCharmBaryons{3.};                                       // maximum Nsigma TOF for protons in Lc and Xic decays
  float mNSigmaTpcKaCutFor3Prongs{3.};                                            // maximum Nsigma TPC for kaons in 3-prong decays
  float mNSigmaTofKaCutFor3Prongs{3.};                                            // maximum Nsigma TOF for kaons in 3-prong decays
  float mNSigmaTpcPiKaCutForDzero{3.};                                            // maximum Nsigma TPC for pions/kaons in D0 decays
  float mNSigmaTofPiKaCutForDzero{3.};                                            // maximum Nsigma TOF for pions/kaons in D0 decays
  float mNSigmaTpcPrKaCutForBeautyToJPsi{3.};                                     // maximum Nsigma TPC for kaons and protons in B->JPsiX decays
  float mNSigmaTofPrKaCutForBeautyToJPsi{3.};                                     // maximum Nsigma TPC for kaons and protons in B->JPsiX decays
  float mDeltaMassMinSigmaCZero{0.155};                                           // minimum delta mass M(pKpipi)-M(pKpi) of SigmaC0 candidates
  float mDeltaMassMaxSigmaCZero{0.18};                                            // maximum delta mass M(pKpipi)-M(pKpi) of SigmaC0 candidates
  float mDeltaMassMinSigmaC2520Zero{0.2};                                         // minimum delta mass M(pKpipi)-M(pKpi) of SigmaC(2520)0 candidates
  float mDeltaMassMaxSigmaC2520Zero{0.26};                                        // maximum delta mass M(pKpipi)-M(pKpi) of SigmaC(2520)0 candidates
  float mDeltaMassMinSigmaCPlusPlus{0.155};                                       // minimum delta mass M(pKpipi)-M(pKpi) of SigmaC++ candidates
  float mDeltaMassMaxSigmaCPlusPlus{0.18};                                        // maximum delta mass M(pKpipi)-M(pKpi) of SigmaC++ candidates
  float mDeltaMassMinSigmaC2520PlusPlus{0.2};                                     // minimum delta mass M(pKpipi)-M(pKpi) of SigmaC(2520)++ candidates
  float mDeltaMassMaxSigmaC2520PlusPlus{0.26};                                    // maximum delta mass M(pKpipi)-M(pKpi) of SigmaC(2520)++ candidates
  float mMinGammaCosinePa{0.85};                                                  // minimum cosp for gammas
  float mMinK0sLambdaCosinePa{0.97};                                              // minimum cosp for K0S and Lambda in charm excited decays
  float mMinK0sLambdaRadius{0.5};                                                 // minimum radius for K0S and Lambda in charm excited decays
  float mMaxNsigmaPrForLambda{4.};                                                // maximum Nsigma TPC and TOF for protons in Lambda decays
  float mDeltaMassK0s{0.02};                                                      // delta mass cut for K0S in charm excited decays
  float mDeltaMassLambda{0.01};                                                   // delta mass cut for Lambda in charm excited decays
  float mMinPtXiBachelor{0.1};                                                    // minimum pt for Xi bachelor in Xic/Omegac decays
  float mMinPtXi{1.};                                                             // minimum pt for Xi in Xic/Omegac decays
  float mDeltaMassXi{0.01};                                                       // delta mass cut for Xi in Xic/Omegac decays
  float mDeltaMassLambdaFromXi{0.01};                                             // delta mass cut for Lambda <- Xi in Xic/Omegac decays
  float mCosPaXi{0.99};                                                           // minimum cosp for Xi in Xic/Omegac decays
  float mCosPaLambdaFromXi{0.99};                                                 // minimum cosp for Xi in Xic/Omegac decays
  float mMaxDcaXyXi{0.3};                                                         // maximum dca for Xi in Xic/Omegac decays
  float mMaxNsigmaXiDau{3.};                                                      // maximum Nsigma TPC and TOF for Xi daughter tracks
  o2::framework::LabeledArray<double> mCutsSingleTrackCharmBaryonBachelor{};      // dca selections for the bachelor pion from Xic/Omegac decays
  float mNSigmaTpcPiCharmBaryonBachelor{3.};                                      // maximum Nsigma TPC for pions in Xic/Omegac decays
  float mNSigmaTofPiCharmBaryonBachelor{3.};                                      // maximum Nsigma TOF for pions in Xic/Omegac decays
  float mNumSigmaDeltaMassCharmHad{2.5};                                          // number of sigmas for delta mass cut for charm hadrons in B and charm excited decays
  std::array<float, 2> mSigmaPars2Prongs{};                                       // parameters (intercept, slope) for parametrisation of mass sigma vs pT for 2-prongs
  std::array<float, 2> mDeltaMassPars2Prongs{};                                   // parameters (intercept, slope) for parametrisation of mass delta wrt PDG vs pT for 2-prongs
  std::array<float, 2> mSigmaPars3Prongs{};                                       // parameters (intercept, slope) for parametrisation of mass sigma vs pT for 3-prongs
  std::array<float, 2> mDeltaMassPars3Prongs{};                                   // parameters (intercept, slope) for parametrisation of mass delta wrt PDG vs pT for 3-prongs
  float mPtThresholdHighPt2Prongs{8.};                                            // threshold for high pT triggers for 2-prongs
  float mPtThresholdHighPt3Prongs{8.};                                            // threshold for high pT triggers for 3-prongs
  float mNSigmaTpcKaonFromXicResoToSigmaC{3.};                                    // maximum Nsigma TPC for kaons in Xic*->SigmaC-Kaon
  float mNSigmaTofKaonFromXicResoToSigmaC{3.};                                    // maximum Nsigma TOF for kaons in Xic*->SigmaC-Kaon
  bool mForceTofProtonForFemto = true;                                            // flag to force TOF PID for protons
  bool mForceTofDeuteronForFemto = false;                                         // flag to force TOF PID for deuterons
  std::array<float, 3> mPtMinXiBach{5., 5., 5.};                                  // minimum pT for XiBachelor candidates
  std::array<float, 3> mMassMinXiBach{2.35, 2.6, 2.35};                           // minimum invariant-mass for XiBachelor candidates
  std::array<float, 3> mMassMaxXiBach{3.0, 3.0, 2.7};                             // maximum invariant-mass for XiBachelor candidates
  std::array<float, 2> mCosPaMinXiBach{-2.f, -2.f};                               // minimum cosine of pointing angle for XiBachelor candidates
  std::array<float, 2> mChi2PcaMaxXiBach{1.e6, 1.e6};                             // minimum chi2PCA for XiBachelor candidates
  std::array<float, 2> mDecLenMinXiBach{-2.f, -2.f};                              // minimum decay length for XiBachelor candidates
  std::array<o2::framework::LabeledArray<double>, kNBeautyParticles> mCutsBhad{}; // selections for B-hadron candidates (DeltaMass, CPA, DecayLength, ImpactParameterProduct)
  o2::framework::LabeledArray<double> mCutsBhadToJPsi{};                          // selections for B->JPsi candidates (PtMinMu, DeltaMass, CPA, DecayLength)
  float mMinTpcCluster{90.};                                                      // Minimum number of TPC clusters required on a track
  float mMinTpcRow{80.};                                                          // Minimum number of TPC rows (pad rows) traversed by the track
  float mMinTpcCrossedOverFound{0.83};                                            // Minimum ratio of crossed TPC rows over findable clusters
  float mMaxTpcShared{160.};                                                      // Maximum allowed number of shared TPC clusters between tracks
  float mMaxTpcFracShared{1.};                                                    // Maximum allowed fraction of shared TPC clusters relative to total clusters
  float mMinItsCluster{1.};                                                       // Minimum required number of ITS clusters
  float mMinItsIbCluster{1.};                                                     // Minimum required number of ITS clusters for IB
// SigmaC–p (ScPr) trigger
float mMinDeltaMassScSigmaCPr{0.15f};                                             // min Delta mass (SigmaC) for SigmaC-Proton trigger
float mMaxDeltaMassScSigmaCPr{0.19f};                                             // max Delta mass (SigmaC) for SigmaC-Proton trigger
float mMinPtScSigmaPr{4.99f};                                                     // min pT(SigmaC) for SigmaC-Proton trigger
float mMinPtPrSigmaCPr{0.399f};                                                    // min pT(proton) for SigmaC-Proton trigger
float mMaxPtPrSigmaCPr{4.501f};                                                    // max pT(proton) for SigmaC-Proton trigger
bool  mForceTOFForPrSigmaCPr{true};                                               // force TOF for protonfor SigmaC-Proton trigger
float mThresholdPtTOFForPrSigmaCPr{1.0f};                                          // pT threshold above which TOF is required for SigmaC-Proton trigger

  // PID recalibrations
  int mTpcPidCalibrationOption{0};                          // Option for TPC PID calibration (0 -> AO2D, 1 -> postcalibrations, 2 -> alternative bethe bloch parametrisation)
  std::array<TH3F*, 8> mHistMapPiPrKaDe{};                  // Map for TPC PID postcalibrations for pions, kaon, protons and deuterons
  std::array<std::vector<double>, 8> mBetheBlochPiKaPrDe{}; // Bethe-Bloch parametrisations for pions, antipions, kaons, antikaons, protons, antiprotons, deuterons, antideuterons in TPC
  // Ds cuts from track-index-skim-creator
  std::vector<double> mPtBinsPreselDsToKKPi{};           // pT bins for pre-selections for Ds from track-index-skim-creator
  o2::framework::LabeledArray<double> mPreselDsToKKPi{}; // pre-selections for Ds from track-index-skim-creator
};

/// Selection of high-pt 2-prong candidates
/// \param pt is the pt of the 2-prong candidate
template <typename T>
inline bool HfFilterHelper::isSelectedHighPt2Prong(const T& pt)
{
  if (pt < mPtThresholdHighPt2Prongs) {
    return false;
  }
  return true;
}

/// Selection of high-pt 3-prong candidates
/// \param pt is the pt of the 3-prong candidate
template <typename T>
inline bool HfFilterHelper::isSelectedHighPt3Prong(const T& pt)
{
  if (pt < mPtThresholdHighPt3Prongs) {
    return false;
  }
  return true;
}

/// Single-track cuts for bachelor track of beauty candidates
/// \param track is a track parameter
/// \param trackPar is a track parameter
/// \param dca is the 2d array with dcaXY and dcaZ of the track
/// \return a flag that encodes the selection for soft pions BIT(kSoftPion), tracks for beauty BIT(kForBeauty), or soft pions for beauty BIT(kSoftPionForBeauty)
template <o2::aod::hffilters::HfTriggers whichTrigger, typename T, typename T1, typename T2>
inline int16_t HfFilterHelper::isSelectedTrackForSoftPionOrBeauty(const T& track, const T1& trackPar, const T2& dca)
{

  int16_t retValue{BIT(kSoftPion) | BIT(kForBeauty) | BIT(kSoftPionForBeauty) | BIT(kSoftPionForSigmaC)};

  if (!track.isGlobalTrackWoDCA()) {
    return kRejected;
  }

  auto pT = trackPar.getPt();
  auto pTBinTrack = findBin(mPtBinsTracks, pT);
  if (pTBinTrack == -1) {
    return kRejected;
  }

  // D*+ soft pion pt cut
  // We can keep ot for all triggers (SigmaC ones included), assuming that the D* soft pion is the softest
  if (pT < mPtMinSoftPionForDstar) { // soft pion min pT cut should be less stringent than usual tracks
    return kRejected;
  }

  if (std::fabs(trackPar.getEta()) > 0.8) {
    return kRejected;
  }

  if (std::fabs(dca[1]) > 2.f) {
    return kRejected;
  }

  if constexpr (whichTrigger == kSigmaCPPK || whichTrigger == kSigmaC0K0 || whichTrigger == kSigmaCPr) {

    // SigmaC0,++ soft pion pt cut
    if (pT < mPtMinSoftPionForSigmaC || pT > mPtMaxSoftPionForSigmaC) {
      return kRejected;
    }

    // We do not need any further selection for SigmaC soft-pi
    // The current track is a good SigmaC soft-pi candidate
    return retValue;
  }

  if (pT > mPtMaxSoftPionForDstar) {
    CLRBIT(retValue, kSoftPion);
    CLRBIT(retValue, kSoftPionForBeauty);
  }

  // below only regular beauty tracks, not required for soft pions
  float ptMin{-1.f}, ptMax{1000.f};
  if constexpr (whichTrigger == kBeauty3P || whichTrigger == kBeauty4P) {
    ptMin = mPtMinBeautyBachelor;
    ptMax = mPtMaxBeautyBachelor;
  } else if constexpr (whichTrigger == kBtoJPsiKa || whichTrigger == kBtoJPsiPi || whichTrigger == kBtoJPsiKstar || whichTrigger == kBtoJPsiPhi || whichTrigger == kBtoJPsiPrKa) {
    ptMin = mPtMinBeautyToJPsiBachelor;
    ptMax = mPtMaxBeautyToJPsiBachelor;
  }

  if (pT < ptMin || pT > ptMax) {
    CLRBIT(retValue, kForBeauty);
  }

  float minDca{1000.f}, maxDca{0.f};
  if constexpr (whichTrigger == kBeauty3P) {
    minDca = mCutsSingleTrackBeauty3Prong.get(pTBinTrack, 0u);
    maxDca = mCutsSingleTrackBeauty3Prong.get(pTBinTrack, 1u);
  } else if constexpr (whichTrigger == kBeauty4P) {
    minDca = mCutsSingleTrackBeauty4Prong.get(pTBinTrack, 0u);
    maxDca = mCutsSingleTrackBeauty4Prong.get(pTBinTrack, 1u);
  } else if constexpr (whichTrigger == kBtoJPsiKa || whichTrigger == kBtoJPsiPi || whichTrigger == kBtoJPsiKstar || whichTrigger == kBtoJPsiPhi || whichTrigger == kBtoJPsiPrKa) {
    minDca = mCutsSingleTrackBeautyToJPsi.get(pTBinTrack, 0u);
    maxDca = mCutsSingleTrackBeautyToJPsi.get(pTBinTrack, 1u);
  }

  if (std::fabs(dca[0]) < minDca) { // minimum DCAxy
    CLRBIT(retValue, kForBeauty);
    CLRBIT(retValue, kSoftPionForBeauty);
  }
  if (std::fabs(dca[0]) > maxDca) { // maximum DCAxy
    CLRBIT(retValue, kForBeauty);
    CLRBIT(retValue, kSoftPionForBeauty);
  }

  return retValue;
}

/// Basic selection of proton or deuteron candidates
/// \param track is a track
/// \param trackPar is a track parameter
/// \param activateQA flag to activate the filling of QA histos
/// \param hProtonTPCPID histo with NsigmaTPC vs. p
/// \param hProtonTOFPID histo with NsigmaTOF vs. p
/// \param trackSpecies flag to choose proton or deuteron
/// \return true if track passes all cuts
template <typename T1, typename T2, typename H2>
inline bool HfFilterHelper::isSelectedTrack4Corr(const T1& track, const T2& trackPar, const int& activateQA, H2 hTPCPID, H2 hTOFPID, const int& trackSpecies)
{
  float pt = trackPar.getPt();
  float ptMin, ptMax, ptThresholdPidStrategy;
  std::array<float, 4> nSigmaCuts;
  bool forceTof = false; // flag to force TOF PID

  // Assign particle-specific parameters
  switch (trackSpecies) {
    case kProtonForFemto:
      ptMin = mPtMinProtonForFemto;
      ptMax = mPtMaxProtonForFemto;
      nSigmaCuts = mNSigmaPrCutsForFemto;
      forceTof = mForceTofProtonForFemto;
      ptThresholdPidStrategy = mPtThresholdProtonForFemto;
      break;
    case kDeuteronForFemto:
      ptMin = mPtMinDeuteronForFemto;
      ptMax = mPtMaxDeuteronForFemto;
      nSigmaCuts = mNSigmaDeCutsForFemto;
      forceTof = mForceTofDeuteronForFemto;
      ptThresholdPidStrategy = mPtThresholdDeuteronForFemto;
      break;
    case kProtonForScPrCorr:
      ptMin = mMinPtPrSigmaCPr;
      ptMax = mMaxPtPrSigmaCPr;
      nSigmaCuts = mNSigmaPrCutsForSigmaCPr;
      forceTof = mForceTOFForPrSigmaCPr;
      ptThresholdPidStrategy = mThresholdPtTOFForPrSigmaCPr;
      break;
    default:
      return false; // Unknown particle type
  }

  // Common selection criteria
  if (pt < ptMin || pt > ptMax) {
    return false;
  }

  if (std::fabs(trackPar.getEta()) > 0.8) {
    return false;
  }

  if (!track.isGlobalTrack()) {
    return false; // use only global tracks
  }
  // PID evaluation
  float NSigmaITS = (trackSpecies == kDeuteronForFemto) ? track.itsNSigmaDe() : track.itsNSigmaPr(); // only used for deuteron
  float NSigmaTPC = (trackSpecies == kDeuteronForFemto) ? track.tpcNSigmaDe() : track.tpcNSigmaPr();
  float NSigmaTOF = (trackSpecies == kDeuteronForFemto) ? track.tofNSigmaDe() : track.tofNSigmaPr();
  if (!forceTof && !track.hasTOF()) {
    NSigmaTOF = 0.; // always accepted
  }

  // Apply TPC PID post-calibration(only available for proton, dummy for deuteron)
  if (mTpcPidCalibrationOption == 1) {
    NSigmaTPC = getTPCPostCalib(track, trackSpecies == kDeuteronForFemto ? kDe : kPr);
  } else if (mTpcPidCalibrationOption == 2) {
    if (track.sign() > 0) {
      NSigmaTPC = getTPCSplineCalib(track, trackSpecies == kDeuteronForFemto ? kDe : kPr);
    } else {
      NSigmaTPC = getTPCSplineCalib(track, trackSpecies == kDeuteronForFemto ? kAntiDe : kAntiPr);
    }
  }

  float NSigma = std::sqrt(NSigmaTPC * NSigmaTPC + NSigmaTOF * NSigmaTOF);
  float momentum = track.p();
  if (trackSpecies == kProtonForFemto || trackSpecies == kProtonForScPrCorr) {
    if (momentum <= ptThresholdPidStrategy) {
      if (NSigma > nSigmaCuts[2]) {
        return false;
      }
    } else {
      if (std::fabs(NSigmaTPC) > nSigmaCuts[0] || std::fabs(NSigmaTOF) > nSigmaCuts[1]) {
        return false;
      }
    }
  }
  // For deuterons: Determine whether to apply TOF based on pt threshold
  if (trackSpecies == kDeuteronForFemto) {

    if (track.tpcNClsFound() < mMinTpcCluster) {
      return false;
    }
    if (track.tpcNClsCrossedRows() < mMinTpcRow) {
      return false;
    }
    if (track.tpcCrossedRowsOverFindableCls() < mMinTpcCrossedOverFound) {
      return false;
    }
    if (track.tpcNClsShared() > mMaxTpcShared) {
      return false;
    }
    if (track.tpcFractionSharedCls() > mMaxTpcFracShared) {
      return false;
    }
    if (track.itsNCls() < mMinItsCluster) {
      return false;
    }
    if (track.itsNClsInnerBarrel() < mMinItsIbCluster) {
      return false;
    }

    // Apply different PID strategy in different pt range
    // one side selection only
    if (momentum <= ptThresholdPidStrategy) {
      if (std::fabs(NSigmaTPC) > nSigmaCuts[0] || NSigmaITS < -nSigmaCuts[3]) { // Use TPC and ITS below the threshold, NSigmaITS for deuteron with a lower limit
        return false;
      }
    } else {
      if (NSigma > nSigmaCuts[2]) { // Use combined TPC and TOF above the threshold
        return false;
      }
    }
  }

  if (activateQA > 1) {
    hTPCPID->Fill(track.p(), NSigmaTPC);
    if ((!forceTof || track.hasTOF())) {
      if (trackSpecies == kProtonForFemto)
        hTOFPID->Fill(momentum, NSigmaTOF);
      else if (trackSpecies == kDeuteronForFemto && momentum > ptThresholdPidStrategy)
        hTOFPID->Fill(momentum, NSigmaTOF);
    }
  }

  return true;
}

/// Basic additional selection of D+ candidates
/// \param trackOppositeCharge is the opposite charge track
/// \param mNSigmaTpcKaCutFor3Prongs max NsigmaTPC for kaon candidates
/// \param mNSigmaTofKaCutFor3Prongs max NsigmaTOF for kaon candidates
/// \return BIT(0) for Kpipi
template <typename T>
inline int8_t HfFilterHelper::isDplusPreselected(const T& trackOppositeCharge)
{
  int8_t retValue = 0;

  // check PID of opposite charge track
  if (!isSelectedKaon4Charm3ProngOrBeautyToJPsi(trackOppositeCharge)) {
    return retValue;
  }

  retValue |= BIT(0);
  return retValue;
}

/// Basic additional selection of Ds candidates
/// \param pTrackSameChargeFirst is the first same-charge track momentum
/// \param pTrackSameChargeSecond is the second same-charge track momentum
/// \param pTrackOppositeCharge is the opposite charge track momentum
/// \param trackOppositeCharge is the opposite charge track
/// \return BIT(0) for KKpi, BIT(1) for piKK
template <typename P, typename T>
inline int8_t HfFilterHelper::isDsPreselected(const P& pTrackSameChargeFirst, const P& pTrackSameChargeSecond, const P& pTrackOppositeCharge, const T& trackOppositeCharge)
{
  int8_t retValue = 0;

  // check PID of opposite charge track
  if (!isSelectedKaon4Charm3ProngOrBeautyToJPsi(trackOppositeCharge)) {
    return retValue;
  }

  // check delta-mass for phi resonance
  auto ptDs = RecoDecay::pt(pTrackSameChargeFirst, pTrackSameChargeSecond, pTrackOppositeCharge);
  auto ptBinDs = findBin(mPtBinsPreselDsToKKPi, ptDs);
  if (ptBinDs == -1) {
    return retValue;
  }

  auto invMassKKFirst = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackOppositeCharge}, std::array{massKa, massKa});
  auto invMassKKSecond = RecoDecay::m(std::array{pTrackSameChargeSecond, pTrackOppositeCharge}, std::array{massKa, massKa});

  float cutValueMassKK = mPreselDsToKKPi.get(ptBinDs, 4u);
  if (std::fabs(invMassKKFirst - massPhi) < cutValueMassKK) {
    retValue |= BIT(0);
  }
  if (std::fabs(invMassKKSecond - massPhi) < cutValueMassKK) {
    retValue |= BIT(1);
  }

  return retValue;
}

/// Basic additional selection of Lc->pKpi and Xic->pKpi candidates
/// \param trackSameChargeFirst is the first same-charge track
/// \param trackSameChargeSecond is the second same-charge track
/// \param trackOppositeCharge is the opposite charge track
/// \return BIT(0) for pKpi, BIT(1) for piKp
template <typename T>
inline int8_t HfFilterHelper::isCharmBaryonPreselected(const T& trackSameChargeFirst, const T& trackSameChargeSecond, const T& trackOppositeCharge)
{
  int8_t retValue = 0;
  // check PID of opposite charge track
  if (!isSelectedKaon4Charm3ProngOrBeautyToJPsi(trackOppositeCharge)) {
    return retValue;
  }
  if (isSelectedProton4CharmOrBeautyBaryons(trackSameChargeFirst)) {
    retValue |= BIT(0);
  }
  if (isSelectedProton4CharmOrBeautyBaryons(trackSameChargeSecond)) {
    retValue |= BIT(1);
  }

  return retValue;
}

/// Basic additional selection of D0 candidates
/// \param trackPos is the positive track
/// \param trackNeg is the negative track
/// \return BIT(0) for D0, BIT(1) for D0bar
template <typename T>
inline int8_t HfFilterHelper::isDzeroPreselected(const T& trackPos, const T& trackNeg)
{
  int8_t retValue = 0;

  float NSigmaPiTPCPos = trackPos.tpcNSigmaPi();
  float NSigmaPiTOFPos = trackPos.tofNSigmaPi();
  float NSigmaKaTPCPos = trackPos.tpcNSigmaKa();
  float NSigmaKaTOFPos = trackPos.tofNSigmaKa();

  float NSigmaPiTPCNeg = trackNeg.tpcNSigmaPi();
  float NSigmaPiTOFNeg = trackNeg.tofNSigmaPi();
  float NSigmaKaTPCNeg = trackNeg.tpcNSigmaKa();
  float NSigmaKaTOFNeg = trackNeg.tofNSigmaKa();

  if (mTpcPidCalibrationOption == 1) {
    NSigmaPiTPCPos = getTPCPostCalib(trackPos, kPi);
    NSigmaPiTPCNeg = getTPCPostCalib(trackNeg, kPi);
    NSigmaKaTPCPos = getTPCPostCalib(trackPos, kKa);
    NSigmaKaTPCNeg = getTPCPostCalib(trackNeg, kKa);
  } else if (mTpcPidCalibrationOption == 2) {
    NSigmaPiTPCPos = getTPCSplineCalib(trackPos, kPi);
    NSigmaPiTPCNeg = getTPCSplineCalib(trackNeg, kAntiPi);
    NSigmaKaTPCPos = getTPCSplineCalib(trackPos, kKa);
    NSigmaKaTPCNeg = getTPCSplineCalib(trackNeg, kAntiKa);
  }

  if ((std::fabs(NSigmaPiTPCPos) <= mNSigmaTpcPiKaCutForDzero && (!trackPos.hasTOF() || std::fabs(NSigmaPiTOFPos) <= mNSigmaTofPiKaCutForDzero)) && (std::fabs(NSigmaKaTPCNeg) <= mNSigmaTpcPiKaCutForDzero && (!trackNeg.hasTOF() || std::fabs(NSigmaKaTOFNeg) <= mNSigmaTofPiKaCutForDzero))) {
    retValue |= BIT(0);
  }
  if ((std::fabs(NSigmaPiTPCNeg) <= mNSigmaTpcPiKaCutForDzero && (!trackNeg.hasTOF() || std::fabs(NSigmaPiTOFNeg) <= mNSigmaTofPiKaCutForDzero)) && (std::fabs(NSigmaKaTPCPos) <= mNSigmaTpcPiKaCutForDzero && (!trackPos.hasTOF() || std::fabs(NSigmaKaTOFPos) <= mNSigmaTofPiKaCutForDzero))) {
    retValue |= BIT(1);
  }

  return retValue;
}

/// Mass selection of D0 candidates to build Bplus candidates
/// \param pTrackPos is the positive track momentum
/// \param pTrackNeg is the negative track momentum
/// \param ptD is the pt of the D0 meson candidate
/// \param isSelected is the flag containing the selection tag for the D0 candidate
/// \param activateQA flag to activate the filling of QA histos
/// \param hMassVsPt histo with invariant mass vs pt
/// \return 1 for D0, 2 for D0bar, 3 for both
template <typename T, typename H2>
inline int8_t HfFilterHelper::isSelectedD0InMassRange(const T& pTrackPos, const T& pTrackNeg, const float& ptD, int8_t isSelected, const int& activateQA, H2 hMassVsPt)
{
  float peakMean = (ptD < 10) ? ((massD0 + mDeltaMassPars2Prongs[0]) + mDeltaMassPars2Prongs[1] * ptD) : massD0;
  float peakWidth = mSigmaPars2Prongs[0] + mSigmaPars2Prongs[1] * ptD;

  int8_t retValue = 0;
  if (TESTBIT(isSelected, 0)) {
    auto invMassD0 = RecoDecay::m(std::array{pTrackPos, pTrackNeg}, std::array{massPi, massKa});
    if (activateQA) {
      hMassVsPt->Fill(ptD, invMassD0);
    }
    if (std::fabs(invMassD0 - peakMean) < mNumSigmaDeltaMassCharmHad * peakWidth || ptD > mPtThresholdHighPt2Prongs) {
      retValue |= BIT(0);
    }
  }
  if (TESTBIT(isSelected, 1)) {
    auto invMassD0bar = RecoDecay::m(std::array{pTrackPos, pTrackNeg}, std::array{massKa, massPi});
    if (activateQA) {
      hMassVsPt->Fill(ptD, invMassD0bar);
    }
    if (std::fabs(invMassD0bar - peakMean) < mNumSigmaDeltaMassCharmHad * peakWidth || ptD > mPtThresholdHighPt2Prongs) {
      retValue |= BIT(1);
    }
  }

  return retValue;
}

/// Mass selection of D+ candidates to build B0 candidates
/// \param pTrackSameChargeFirst is the first same-charge track momentum
/// \param pTrackSameChargeFirst is the second same-charge track momentum
/// \param pTrackSameChargeFirst is the opposite charge track momentum
/// \param ptD is the pt of the D+ meson candidate
/// \param activateQA flag to activate the filling of QA histos
/// \param hMassVsPt histo with invariant mass vs pt
/// \return BIT(0) (==1) for D+, 0 otherwise
template <typename T, typename H2>
inline int8_t HfFilterHelper::isSelectedDplusInMassRange(const T& pTrackSameChargeFirst, const T& pTrackSameChargeSecond, const T& pTrackOppositeCharge, const float& ptD, const int& activateQA, H2 hMassVsPt)
{
  float peakMean = (ptD < 10) ? ((massDPlus + mDeltaMassPars3Prongs[0]) + mDeltaMassPars3Prongs[1] * ptD) : massDPlus;
  float peakWidth = mSigmaPars3Prongs[0] + mSigmaPars3Prongs[1] * ptD;

  auto invMassDplus = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackSameChargeSecond, pTrackOppositeCharge}, std::array{massPi, massPi, massKa});
  if (activateQA) {
    hMassVsPt->Fill(ptD, invMassDplus);
  }

  if (std::fabs(invMassDplus - peakMean) > mNumSigmaDeltaMassCharmHad * peakWidth && ptD < mPtThresholdHighPt3Prongs) {
    return 0;
  }

  return BIT(0);
}

/// Mass selection of of Ds candidates to build Bs candidates
/// \param pTrackSameChargeFirst is the first same-charge track momentum
/// \param pTrackSameChargeFirst is the second same-charge track momentum
/// \param pTrackSameChargeFirst is the opposite charge track momentum
/// \param ptD is the pt of the Ds meson candidate
/// \param isSelected is the flag containing the selection tag for the Ds candidate
/// \param activateQA flag to activate the filling of QA histos
/// \param hMassVsPt histo with invariant mass vs pt
/// \return BIT(0) for KKpi, BIT(1) for piKK, BIT(2) for phipi, BIT(3) for piphi
template <typename T, typename H2>
inline int8_t HfFilterHelper::isSelectedDsInMassRange(const T& pTrackSameChargeFirst, const T& pTrackSameChargeSecond, const T& pTrackOppositeCharge, const float& ptD, int8_t isSelected, const int& activateQA, H2 hMassVsPt)
{
  float peakMean = (ptD < 10) ? ((massDs + mDeltaMassPars3Prongs[0]) + mDeltaMassPars3Prongs[1] * ptD) : massDs;
  float peakWidth = mSigmaPars3Prongs[0] + mSigmaPars3Prongs[1] * ptD;

  int8_t retValue = 0;
  if (TESTBIT(isSelected, 0)) {
    auto invMassDsToKKPi = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackOppositeCharge, pTrackSameChargeSecond}, std::array{massKa, massKa, massPi});
    if (activateQA) {
      hMassVsPt->Fill(ptD, invMassDsToKKPi);
    }
    if (std::fabs(invMassDsToKKPi - peakMean) < peakWidth * mNumSigmaDeltaMassCharmHad || ptD > mPtThresholdHighPt3Prongs) {
      retValue |= BIT(0);
    }
  }
  if (TESTBIT(isSelected, 1)) {
    auto invMassDsToPiKK = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackOppositeCharge, pTrackSameChargeSecond}, std::array{massPi, massKa, massKa});
    if (activateQA) {
      hMassVsPt->Fill(ptD, invMassDsToPiKK);
    }
    if (std::fabs(invMassDsToPiKK - peakMean) < peakWidth * mNumSigmaDeltaMassCharmHad || ptD > mPtThresholdHighPt3Prongs) {
      retValue |= BIT(1);
    }
  }

  return retValue;
}

/// Mass selection of Lc candidates to build Lb candidates
/// \param pTrackSameChargeFirst is the first same-charge track momentum
/// \param pTrackSameChargeSecond is the second same-charge track momentum
/// \param pTrackOppositeCharge is the opposite charge track momentum
/// \param ptLc is the pt of the D0 meson candidate
/// \param isSelected is the flag containing the selection tag for the D0 candidate
/// \param activateQA flag to activate the filling of QA histos
/// \param hMassVsPt histo with invariant mass vs pt
/// \return BIT(0) for pKpi with mass cut, BIT(1) for piKp with mass cut
template <typename T, typename H2>
inline int8_t HfFilterHelper::isSelectedLcInMassRange(const T& pTrackSameChargeFirst, const T& pTrackSameChargeSecond, const T& pTrackOppositeCharge, const float& ptLc, const int8_t isSelected, const int& activateQA, H2 hMassVsPt)
{
  float peakMean = (ptLc < 10) ? ((massLc + mDeltaMassPars3Prongs[0]) + mDeltaMassPars3Prongs[1] * ptLc) : massLc;
  float peakWidth = mSigmaPars3Prongs[0] + mSigmaPars3Prongs[1] * ptLc;

  int8_t retValue = 0;
  if (TESTBIT(isSelected, 0)) {
    auto invMassLcToPKPi = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackOppositeCharge, pTrackSameChargeSecond}, std::array{massProton, massKa, massPi});
    if (activateQA) {
      hMassVsPt->Fill(ptLc, invMassLcToPKPi);
    }
    if (std::fabs(invMassLcToPKPi - peakMean) < peakWidth * mNumSigmaDeltaMassCharmHad || ptLc > mPtThresholdHighPt3Prongs) {
      retValue |= BIT(0);
    }
  }
  if (TESTBIT(isSelected, 1)) {
    auto invMassLcToPiKP = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackOppositeCharge, pTrackSameChargeSecond}, std::array{massPi, massKa, massProton});
    if (activateQA) {
      hMassVsPt->Fill(ptLc, invMassLcToPiKP);
    }
    if (std::fabs(invMassLcToPiKP - peakMean) < peakWidth * mNumSigmaDeltaMassCharmHad || ptLc > mPtThresholdHighPt3Prongs) {
      retValue |= BIT(1);
    }
  }

  return retValue;
}

/// Delta mass selection on SigmaC candidates
template <int charge, typename T, typename H2>
inline int8_t HfFilterHelper::isSelectedSigmaCInDeltaMassRange(const T& pTrackSameChargeFirst, const T& pTrackSameChargeSecond, const T& pTrackOppositeCharge, const T& pTrackSoftPi, const float ptSigmaC, const int8_t isSelectedLc, H2 hMassVsPt, const int& activateQA)
{
  int8_t retValue = 0;
  if (TESTBIT(isSelectedLc, 0)) {
    /// Lc->pKpi case
    auto invMassLcToPKPi = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackOppositeCharge, pTrackSameChargeSecond}, std::array{massProton, massKa, massPi});
    std::array<float, 4> massDausSigmaCToLcPKPi{massProton, massKa, massPi, massPi};
    float invMassSigmaCToLcPKPi = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackOppositeCharge, pTrackSameChargeSecond, pTrackSoftPi}, massDausSigmaCToLcPKPi);
    float deltaMassPKPi = invMassSigmaCToLcPKPi - invMassLcToPKPi;
    bool isSigmaC2520{false};
    bool isSigmaC2455{false};
    if constexpr (charge == 0) {
      isSigmaC2455 = (mDeltaMassMinSigmaCZero < deltaMassPKPi && deltaMassPKPi < mDeltaMassMaxSigmaCZero && ptSigmaC > mPtMinSigmaCZero);
      isSigmaC2520 = (mDeltaMassMinSigmaC2520Zero < deltaMassPKPi && deltaMassPKPi < mDeltaMassMaxSigmaC2520Zero && ptSigmaC > mPtMinSigmaC2520Zero);
    } else if constexpr (charge == 2) {
      isSigmaC2455 = (mDeltaMassMinSigmaCPlusPlus < deltaMassPKPi && deltaMassPKPi < mDeltaMassMaxSigmaCPlusPlus && ptSigmaC > mPtMinSigmaCPlusPlus);
      isSigmaC2520 = (mDeltaMassMinSigmaC2520PlusPlus < deltaMassPKPi && deltaMassPKPi < mDeltaMassMaxSigmaC2520PlusPlus && ptSigmaC > mPtMinSigmaC2520PlusPlus);
    } else if constexpr (charge == -1){
      if (deltaMassPKPi > mMinDeltaMassScSigmaCPr && deltaMassPKPi < mMaxDeltaMassScSigmaCPr && ptSigmaC > mMinPtScSigmaPr){ // sigmaC charge independent for SigmaCPr
      SETBIT(retValue, 4); // SigmaCPr bit
      if (activateQA) {
        hMassVsPt->Fill(ptSigmaC, deltaMassPKPi);
      }
      }
    }

    if (isSigmaC2455 || isSigmaC2520) {
      retValue |= BIT(0);
      if (isSigmaC2455) {
        SETBIT(retValue, 2);
      }
      if (isSigmaC2520) {
        SETBIT(retValue, 3);
      }
      /// QA plot
      if (activateQA) {
        hMassVsPt->Fill(ptSigmaC, deltaMassPKPi);
      }
    }
  }
  
  if (TESTBIT(isSelectedLc, 1)) {
    /// Lc->piKp case
    auto invMassLcToPiKP = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackOppositeCharge, pTrackSameChargeSecond}, std::array{massPi, massKa, massProton});
    std::array<float, 4> massDausSigmaCToLcPiKP{massPi, massKa, massProton, massPi};
    float invMassSigmaCToLcPiKP = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackOppositeCharge, pTrackSameChargeSecond, pTrackSoftPi}, massDausSigmaCToLcPiKP);
    float deltaMassPiKP = invMassSigmaCToLcPiKP - invMassLcToPiKP;
    bool isSigmaC2520{false};
    bool isSigmaC2455{false};
    if constexpr (charge == 0) {
      isSigmaC2455 = (mDeltaMassMinSigmaCZero < deltaMassPiKP && deltaMassPiKP < mDeltaMassMaxSigmaCZero && ptSigmaC > mPtMinSigmaCZero);
      isSigmaC2520 = (mDeltaMassMinSigmaC2520Zero < deltaMassPiKP && deltaMassPiKP < mDeltaMassMaxSigmaC2520Zero && ptSigmaC > mPtMinSigmaC2520Zero);
    } else if constexpr (charge == 2) {
      isSigmaC2455 = (mDeltaMassMinSigmaCPlusPlus < deltaMassPiKP && deltaMassPiKP < mDeltaMassMaxSigmaCPlusPlus && ptSigmaC > mPtMinSigmaCPlusPlus);
      isSigmaC2520 = (mDeltaMassMinSigmaC2520PlusPlus < deltaMassPiKP && deltaMassPiKP < mDeltaMassMaxSigmaC2520PlusPlus && ptSigmaC > mPtMinSigmaC2520PlusPlus);
    } else if constexpr (charge == -1){
      if (deltaMassPiKP > mMinDeltaMassScSigmaCPr && deltaMassPiKP < mMaxDeltaMassScSigmaCPr && ptSigmaC > mMinPtScSigmaPr){ // sigmaC charge independent for SigmaCPr
      SETBIT(retValue, 4); // SigmaCPr bit
      if (activateQA) {
        hMassVsPt->Fill(ptSigmaC, deltaMassPiKP);
      }
      }
    }
    if (isSigmaC2455 || isSigmaC2520) {
      retValue |= BIT(1);
      if (isSigmaC2455) {
        SETBIT(retValue, 2);
      }
      if (isSigmaC2520) {
        SETBIT(retValue, 3);
      }
      /// QA plot
      if (activateQA) {
        hMassVsPt->Fill(ptSigmaC, deltaMassPiKP);
      }
    }
  }
  /// TODO: add QA plot

  return retValue;
}

/// Mass selection of Xic candidates to build Xib candidates
/// \param pTrackSameChargeFirst is the first same-charge track momentum
/// \param pTrackSameChargeSecond is the second same-charge track momentum
/// \param pTrackOppositeCharge is the opposite charge track momentum
/// \param ptXic is the pt of the Xic baryon candidate
/// \param isSelected is the flag containing the selection tag for the D0 candidate
/// \param activateQA flag to activate the filling of QA histos
/// \param hMassVsPt histo with invariant mass vs pt
/// \return BIT(0) for pKpi with mass cut, BIT(1) for piKp with mass cut
template <typename T, typename H2>
inline int8_t HfFilterHelper::isSelectedXicInMassRange(const T& pTrackSameChargeFirst, const T& pTrackSameChargeSecond, const T& pTrackOppositeCharge, const float& ptXic, const int8_t isSelected, const int& activateQA, H2 hMassVsPt)
{
  float peakMean = (ptXic < 10) ? ((massLc + mDeltaMassPars3Prongs[0]) + mDeltaMassPars3Prongs[1] * ptXic) : massXic;
  float peakWidth = mSigmaPars3Prongs[0] + mSigmaPars3Prongs[1] * massXic;

  int8_t retValue = 0;
  if (TESTBIT(isSelected, 0)) {
    auto invMassXicToPKPi = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackOppositeCharge, pTrackSameChargeSecond}, std::array{massProton, massKa, massPi});
    if (activateQA) {
      hMassVsPt->Fill(ptXic, invMassXicToPKPi);
    }
    if (std::fabs(invMassXicToPKPi - peakMean) < peakWidth * mNumSigmaDeltaMassCharmHad || ptXic > mPtThresholdHighPt3Prongs) {
      retValue |= BIT(0);
    }
  }
  if (TESTBIT(isSelected, 1)) {
    auto invMassXicToPiKP = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackOppositeCharge, pTrackSameChargeSecond}, std::array{massPi, massKa, massProton});
    if (activateQA) {
      hMassVsPt->Fill(ptXic, invMassXicToPiKP);
    }
    if (std::fabs(invMassXicToPiKP - peakMean) < peakWidth * mNumSigmaDeltaMassCharmHad || ptXic > mPtThresholdHighPt3Prongs) {
      retValue |= BIT(1);
    }
  }

  return retValue;
}

/// Basic selection of V0 candidates
/// \param v0 is the v0 candidate
/// \param activateQA flag to fill QA histos
/// \param hV0Selected is the pointer to the QA histo for selected V0S
/// \param hArmPod is the pointer to an array of QA histo AP plot after selection
/// \return an integer passes all cuts
template <typename V0, typename H2>
inline int8_t HfFilterHelper::isSelectedV0(const V0& v0, const int& activateQA, H2 hV0Selected, std::array<H2, 4>& hArmPod)
{
  int8_t isSelected{BIT(kK0S) | BIT(kLambda) | BIT(kAntiLambda)};

  if (activateQA > 1) {
    for (int iV0{kK0S}; iV0 < kNV0; ++iV0) {
      hV0Selected->Fill(0., iV0);
    }
  }

  // eta of daughters
  if (std::fabs(v0.etaPos) > 1. || std::fabs(v0.etaNeg) > 1.) { // cut all V0 daughters with |eta| > 1.
    if (activateQA > 1) {
      for (int iV0{kK0S}; iV0 < kNV0; ++iV0) {
        hV0Selected->Fill(1., iV0);
      }
    }
    return kRejected;
  }

  // V0 radius
  if (v0.v0radius < mMinK0sLambdaRadius) {
    for (int iV0{kK0S}; iV0 < kNV0; ++iV0) {
      CLRBIT(isSelected, iV0);
      if (activateQA > 1) {
        hV0Selected->Fill(2., iV0);
      }
    }
  }

  for (int iV0{kK0S}; iV0 < kNV0; ++iV0) {
    if (TESTBIT(isSelected, iV0) && v0.v0cosPA < mMinK0sLambdaCosinePa) {
      CLRBIT(isSelected, iV0);
      if (activateQA > 1) {
        hV0Selected->Fill(3., iV0);
      }
    }
  }

  // armenteros-podolanski / mass
  if (TESTBIT(isSelected, kK0S) && std::fabs(v0.mK0Short - massK0S) > mDeltaMassK0s) {
    CLRBIT(isSelected, kK0S);
    if (activateQA > 1) {
      hV0Selected->Fill(4., kK0S);
    }
  }
  if (TESTBIT(isSelected, kLambda) && std::fabs(v0.mLambda - massLambda) > mDeltaMassLambda) {
    CLRBIT(isSelected, kLambda);
    if (activateQA > 1) {
      hV0Selected->Fill(4., kLambda);
    }
  }
  if (TESTBIT(isSelected, kAntiLambda) && std::fabs(v0.mAntiLambda - massLambda) > mDeltaMassLambda) {
    CLRBIT(isSelected, kAntiLambda);
    if (activateQA > 1) {
      hV0Selected->Fill(4., kAntiLambda);
    }
  }

  // DCA V0 and V0 daughters
  for (int iV0{kK0S}; iV0 < kNV0; ++iV0) {
    if (TESTBIT(isSelected, iV0) && v0.dcav0topv > 0.1f) { // we want only primary V0s
      CLRBIT(isSelected, iV0);
      if (activateQA > 1) {
        hV0Selected->Fill(5., iV0);
      }
    }
    if (TESTBIT(isSelected, iV0) && (v0.dcaV0daughters > 1.f || std::fabs(v0.dcapostopv) < 0.05f || std::fabs(v0.dcanegtopv) < 0.05f)) {
      CLRBIT(isSelected, iV0);
      if (activateQA > 1) {
        hV0Selected->Fill(6., iV0);
      }
    }
  }

  // PID (Lambda/AntiLambda only)
  float nSigmaPrTpc[2] = {v0.nSigmaPrTpcPos, v0.nSigmaPrTpcNeg};
  float nSigmaPrTof[2] = {v0.nSigmaPrTofPos, v0.nSigmaPrTofNeg};
  float pInTpc[2] = {v0.pinTpcPos, v0.pinTpcNeg};
  if (mTpcPidCalibrationOption == 1) {
    float nClsTpc[2] = {v0.nClsFoundTpcPos, v0.nClsFoundTpcNeg};
    float etaDaus[2] = {v0.etaPos, v0.etaNeg};
    for (int iDau{0}; iDau < 2; ++iDau) {
      nSigmaPrTpc[iDau] = getTPCPostCalib(pInTpc[iDau], nClsTpc[iDau], etaDaus[iDau], nSigmaPrTpc[iDau], kPr);
    }
  } else if (mTpcPidCalibrationOption == 2) {
    float signalTpc[2] = {v0.signalTpcPos, v0.signalTpcNeg};
    for (int iDau{0}; iDau < 2; ++iDau) {
      nSigmaPrTpc[iDau] = getTPCSplineCalib(pInTpc[iDau], signalTpc[iDau], (iDau == 0) ? kPr : kAntiPr);
    }
  }

  if (TESTBIT(isSelected, kLambda) && (std::fabs(nSigmaPrTpc[0]) > mMaxNsigmaPrForLambda || (v0.hasTofPos && std::fabs(nSigmaPrTof[0]) > mMaxNsigmaPrForLambda))) {
    CLRBIT(isSelected, kLambda);
    if (activateQA > 1) {
      hV0Selected->Fill(7., kLambda);
    }
  }
  if (TESTBIT(isSelected, kAntiLambda) && (std::fabs(nSigmaPrTpc[1]) > mMaxNsigmaPrForLambda || (v0.hasTofNeg && std::fabs(nSigmaPrTof[1]) > mMaxNsigmaPrForLambda))) {
    CLRBIT(isSelected, kAntiLambda);
    if (activateQA > 1) {
      hV0Selected->Fill(7., kAntiLambda);
    }
  }

  if (activateQA) {
    for (int iV0{kK0S}; iV0 < kNV0; ++iV0) {
      if (TESTBIT(isSelected, iV0)) {
        hArmPod[iV0]->Fill(v0.alpha, v0.qtarm);
        if (activateQA > 1) {
          hV0Selected->Fill(8., iV0);
        }
      }
    }
  }

  return isSelected;
}

/// Basic selection of photon candidates
/// \param photon is the photon candidate
/// \param dauTracks is a 2-element array with positive and negative V0 daughter tracks
/// \param activateQA flag to fill QA histos
/// \param hV0Selected is the pointer to the QA histo for selected V0s
/// \param hArmPod is the pointer to an array of QA histo AP plot after selection
/// \return an integer passes all cuts
template <typename Photon, typename T, typename H2>
inline bool HfFilterHelper::isSelectedPhoton(const Photon& photon, const std::array<T, 2>& dauTracks, const int& activateQA, H2 hV0Selected, std::array<H2, 4>& hArmPod)
{

  if (activateQA > 1) {
    hV0Selected->Fill(0., kPhoton);
  }

  // eta of daughters
  if (std::fabs(dauTracks[0].eta()) > 1. || std::fabs(dauTracks[1].eta()) > 1.) { // cut all V0 daughters with |eta| > 1.
    if (activateQA > 1) {
      hV0Selected->Fill(1., kPhoton);
    }
    return false;
  }

  // radius
  if (photon.v0radius() < 0. || photon.v0radius() > 180.) {
    if (activateQA > 1) {
      hV0Selected->Fill(2., kPhoton);
    }
    return false;
  }

  // cosine of pointing angle
  if (photon.cospa() < mMinGammaCosinePa) {
    if (activateQA > 1) {
      hV0Selected->Fill(3., kPhoton);
    }
    return false;
  }

  if (activateQA) {
    hArmPod[kPhoton]->Fill(photon.alpha(), photon.qtarm());
    if (activateQA > 1) {
      hV0Selected->Fill(8., kPhoton);
    }
  }

  return true;
}

/// Basic selection of cascade candidates
/// \param casc is the cascade candidate
/// \return true if cascade passes all cuts
template <typename Casc>
inline bool HfFilterHelper::isSelectedCascade(const Casc& casc)
{

  // Xi min pT
  if (casc.pt < mMinPtXi) {
    return false;
  }

  // eta of daughters
  if (std::fabs(casc.v0.etaPos) > 1. || std::fabs(casc.v0.etaNeg) > 1. || std::fabs(casc.etaBach) > 1.) { // cut all V0 daughters with |eta| > 1.
    return false;
  }

  // V0 radius
  if (casc.v0.v0radius < 1.2) {
    return false;
  }

  // cascade radius
  if (casc.cascradius < 0.6) {
    return false;
  }

  // V0 cosp
  if (casc.v0.v0cosPA < mCosPaLambdaFromXi) {
    return false;
  }

  // cascade cosp
  if (casc.casccosPA < mCosPaXi) {
    return false;
  }

  // cascade DCAxy to PV
  if (std::fabs(casc.dcaXYCascToPV) > mMaxDcaXyXi) {
    return false;
  }

  // Xi bachelor min pT
  if (casc.ptBach < mMinPtXiBachelor) {
    return false;
  }

  // dau dca
  if (std::fabs(casc.v0.dcaV0daughters) > 1.f || std::fabs(casc.dcacascdaughters) > 1.f) {
    return false;
  }

  // cascade mass
  if (std::fabs(casc.mXi - massXi) > mDeltaMassXi) {
    return false;
  }

  // V0 mass
  if ((casc.sign < 0 && std::fabs(casc.v0.mLambda - massLambda) > mDeltaMassLambdaFromXi) || (casc.sign > 0 && std::fabs(casc.v0.mAntiLambda - massLambda) > mDeltaMassLambdaFromXi)) {
    return false;
  }

  // PID
  float nSigmaPrTpc[3] = {-999, casc.v0.nSigmaPrTpcPos, casc.v0.nSigmaPrTpcNeg};
  float nSigmaPrTof[3] = {-999., casc.v0.nSigmaPrTofPos, casc.v0.nSigmaPrTofNeg};
  float nSigmaPiTpc[3] = {casc.nSigmaPiTpcBach, casc.v0.nSigmaPiTpcPos, casc.v0.nSigmaPiTpcNeg};
  float nSigmaPiTof[3] = {casc.nSigmaPiTofBach, casc.v0.nSigmaPiTofPos, casc.v0.nSigmaPiTofNeg};
  float pInTpc[3] = {casc.pinTpcBach, casc.v0.pinTpcPos, casc.v0.pinTpcNeg};
  float nClsTpc[3] = {casc.nClsFoundTpcBach, casc.v0.nClsFoundTpcPos, casc.v0.nClsFoundTpcNeg};
  float nCrossedRowsTpc[3] = {casc.nClsCrossedRowsTpcBach, casc.v0.nClsCrossedRowsTpcPos, casc.v0.nClsCrossedRowsTpcNeg};
  float crossedRowsOverFindableClsTpc[3] = {casc.crossedRowsOverFindableClsTpcBach, casc.v0.crossedRowsOverFindableClsTpcPos, casc.v0.crossedRowsOverFindableClsTpcNeg};
  if (mTpcPidCalibrationOption == 1) {
    float etaDaus[3] = {casc.etaBach, casc.v0.etaPos, casc.v0.etaNeg};
    for (int iDau{0}; iDau < 3; ++iDau) {
      nSigmaPiTpc[iDau] = getTPCPostCalib(pInTpc[iDau], nClsTpc[iDau], etaDaus[iDau], nSigmaPrTpc[iDau], kPi);
      if (iDau == 0) {
        continue;
      }
      nSigmaPrTpc[iDau] = getTPCPostCalib(pInTpc[iDau], nClsTpc[iDau], etaDaus[iDau], nSigmaPrTpc[iDau], kPr);
    }
  } else if (mTpcPidCalibrationOption == 2) {
    float signalTpc[3] = {casc.signalTpcBach, casc.v0.signalTpcPos, casc.v0.signalTpcNeg};
    for (int iDau{0}; iDau < 3; ++iDau) {
      nSigmaPiTpc[iDau] = getTPCSplineCalib(pInTpc[iDau], signalTpc[iDau], (iDau == 0) ? kPi : kAntiPi);
      if (iDau == 0) {
        continue;
      }
      nSigmaPrTpc[iDau] = getTPCSplineCalib(pInTpc[iDau], signalTpc[iDau], kAntiPr);
    }
  }

  // PID to V0 tracks
  if (casc.sign < 0) { // Xi-
    if (std::fabs(nSigmaPrTpc[1]) > mMaxNsigmaXiDau && (casc.v0.hasTofPos && std::fabs(nSigmaPrTof[1]) > mMaxNsigmaXiDau)) {
      return false;
    }
    if (std::fabs(nSigmaPiTpc[2]) > mMaxNsigmaXiDau && (casc.v0.hasTofNeg && std::fabs(nSigmaPiTof[2]) > mMaxNsigmaXiDau)) {
      return false;
    }
  } else if (casc.sign > 0) { // Xi+
    if (std::fabs(nSigmaPrTpc[2]) > mMaxNsigmaXiDau && (casc.v0.hasTofNeg && std::fabs(nSigmaPrTof[2]) > mMaxNsigmaXiDau)) {
      return false;
    }
    if (std::fabs(nSigmaPiTpc[1]) > mMaxNsigmaXiDau && (casc.v0.hasTofPos && std::fabs(nSigmaPiTof[1]) > mMaxNsigmaXiDau)) {
      return false;
    }
  }

  // bachelor PID
  if (std::fabs(nSigmaPiTpc[0]) > mMaxNsigmaXiDau && (casc.hasTofBach && std::fabs(nSigmaPiTof[0]) > mMaxNsigmaXiDau)) {
    return false;
  }

  // additional track cuts
  for (int iTrack{0}; iTrack < 3; ++iTrack) {
    //  TPC clusters selections
    if (nClsTpc[iTrack] < 70) { // TODO: put me as a configurable please
      return false;
    }
    if (nCrossedRowsTpc[iTrack] < 70) {
      return false;
    }
    if (crossedRowsOverFindableClsTpc[iTrack] < 0.8) {
      return false;
    }
  }

  return true;
}

/// Single-track cuts for bachelor track of charm baryon candidates
/// \param track is a track
/// \param dca is the 2d array with dcaXY and dcaZ of the track
/// \return 0 if rejected, or a bitmap that contains the information whether it is selected as pion and/or kaon
template <typename T, typename T2>
inline int16_t HfFilterHelper::isSelectedBachelorForCharmBaryon(const T& track, const T2& dca)
{
  int16_t retValue{BIT(kPionForCharmBaryon) | BIT(kKaonForCharmBaryon)};

  if (!track.isGlobalTrackWoDCA()) {
    return kRejected;
  }

  float pt = track.pt();
  if (pt < mPtMinCharmBaryonBachelor || pt > mPtMaxCharmBaryonBachelor) {
    return kRejected;
  }

  auto pTBinTrack = findBin(mPtBinsTracks, pt);
  if (pTBinTrack == -1) {
    return kRejected;
  }

  if (std::fabs(dca[0]) < mCutsSingleTrackCharmBaryonBachelor.get(pTBinTrack, 0u)) {
    return kRejected; // minimum DCAxy
  }
  if (std::fabs(dca[0]) > mCutsSingleTrackCharmBaryonBachelor.get(pTBinTrack, 1u)) {
    return kRejected; // maximum DCAxy
  }

  if (std::fabs(dca[1]) > 2.f) {
    return kRejected; // maximum DCAz
  }

  if (track.tpcNClsFound() < 70) {
    return kRejected;
  }

  if (track.itsNCls() < 3) {
    return kRejected;
  }

  float nSigmaPiTpc = track.tpcNSigmaPi();
  float nSigmaKaTpc = track.tpcNSigmaKa();
  float nSigmaPiTof = track.tofNSigmaPi();
  float nSigmaKaTof = track.tofNSigmaKa();
  if (mTpcPidCalibrationOption == 1) {
    nSigmaPiTpc = getTPCPostCalib(track, kPi);
    nSigmaKaTpc = getTPCPostCalib(track, kKa);
  } else if (mTpcPidCalibrationOption == 2) {
    nSigmaPiTpc = getTPCSplineCalib(track, (track.sign() > 0) ? kPi : kAntiPi);
    nSigmaKaTpc = getTPCSplineCalib(track, (track.sign() > 0) ? kKa : kAntiKa);
  }

  if ((track.hasTPC() && std::fabs(nSigmaPiTpc) > mNSigmaTpcPiCharmBaryonBachelor) && (track.hasTOF() && std::fabs(nSigmaPiTof) > mNSigmaTofPiCharmBaryonBachelor)) {
    CLRBIT(retValue, kPionForCharmBaryon);
  }
  if ((track.hasTPC() && std::fabs(nSigmaKaTpc) > mNSigmaTpcPiCharmBaryonBachelor) && (track.hasTOF() && std::fabs(nSigmaKaTof) > mNSigmaTofPiCharmBaryonBachelor)) {
    CLRBIT(retValue, kKaonForCharmBaryon);
  }

  return retValue;
}

/// BDT selections
/// \param scores is a 3-element array with BDT out scores
/// \param thresholdBDTScores is the LabelledArray containing the BDT cut values
/// \return 0 if rejected, otherwise bitmap with BIT(RecoDecay::OriginType::Prompt) and/or BIT(RecoDecay::OriginType::NonPrompt) on
template <typename T, typename U>
inline int8_t HfFilterHelper::isBDTSelected(const T& scores, const U& thresholdBDTScores)
{
  int8_t retValue = 0;
  if (scores.size() < 3) {
    return retValue;
  }

  if (scores[0] > thresholdBDTScores.get(0u, 0u)) {
    return retValue;
  }
  retValue |= BIT(RecoDecay::OriginType::None); // signal, but not yet tagged as prompt or nonprompt
  if (scores[1] > thresholdBDTScores.get(0u, 1u)) {
    retValue |= BIT(RecoDecay::OriginType::Prompt);
  }
  if (scores[2] > thresholdBDTScores.get(0u, 2u)) {
    retValue |= BIT(RecoDecay::OriginType::NonPrompt);
  }

  return retValue;
}

/// Computation of the relative momentum between particle pairs
/// \param pTrack is the track momentum array
/// \param ProtonMass is the mass of a proton
/// \param CharmCandMomentum is the three momentum of a charm candidate
/// \param CharmMass is the mass of the charm hadron
/// \return relative momentum of pair
template <typename T>
inline T HfFilterHelper::computeRelativeMomentum(const std::array<T, 3>& pTrack, const std::array<T, 3>& CharmCandMomentum, const T& CharmMass)
{
  ROOT::Math::PxPyPzMVector part1(pTrack[0], pTrack[1], pTrack[2], massProton);
  ROOT::Math::PxPyPzMVector part2(CharmCandMomentum[0], CharmCandMomentum[1], CharmCandMomentum[2], CharmMass);

  ROOT::Math::PxPyPzMVector trackSum = part1 + part2;
  ROOT::Math::Boost boostv12{trackSum.BoostToCM()};
  ROOT::Math::PxPyPzMVector part1CM = boostv12(part1);
  ROOT::Math::PxPyPzMVector part2CM = boostv12(part2);
  ROOT::Math::PxPyPzMVector trackRelK = part1CM - part2CM;

  T kStar = 0.5 * trackRelK.P();
  return kStar;
} // float computeRelativeMomentum(const T& track, const std::array<float, 3>& CharmCandMomentum, const float& CharmMass)

/// Computation of the number of candidates in an event that do not share daughter tracks
/// \return 0 or 1 in case of less than 2 independent candidates in a single event, 2 otherwise
template <typename T>
inline int HfFilterHelper::computeNumberOfCandidates(const std::vector<std::vector<T>>& indices)
{
  if (indices.size() < 2) {
    return indices.size();
  }

  std::vector<int> numIndependentCand{};
  for (auto iCand{0u}; iCand < indices.size(); ++iCand) {
    int nIndependent = 0;
    for (auto iCandSecond{0u}; iCandSecond < indices.size(); ++iCandSecond) {
      if (iCand == iCandSecond) {
        continue;
      } else {
        bool hasOverlap = false;
        for (auto idxFirst{0u}; idxFirst < indices[iCand].size(); ++idxFirst) {
          for (auto idxSecond{0u}; idxSecond < indices[iCandSecond].size(); ++idxSecond) {
            if (indices[iCand][idxFirst] == indices[iCandSecond][idxSecond]) {
              hasOverlap = true;
              break;
            }
          }
        }
        if (!hasOverlap) {
          nIndependent++;
        }
      }
    }
    numIndependentCand.push_back(nIndependent);
  }
  std::sort(numIndependentCand.begin(), numIndependentCand.end());

  if (numIndependentCand.back() == 0) {
    return numIndependentCand.back();
  }

  return 2;
}

/// PID postcalibrations

/// load the TPC spline from the CCDB
/// \param ccdbApi is Api for CCDB
/// \param bunchCrossing is the timestamp of bunchcrossing for the run number
/// \param ccdbPaths  are the paths on CCDB for pions, antipions, kaons, antikaons, protons, antiprotons
inline void HfFilterHelper::setValuesBB(o2::ccdb::CcdbApi& ccdbApi, aod::BCsWithTimestamps::iterator const& bunchCrossing, const std::array<std::string, 8>& ccdbPaths)
{
  for (int iSpecie{0u}; iSpecie < 8; ++iSpecie) {
    std::map<std::string, std::string> metadata;
    auto hSpline = ccdbApi.retrieveFromTFileAny<TH1F>(ccdbPaths[iSpecie], metadata, bunchCrossing.timestamp());

    if (!hSpline) {
      LOG(fatal) << "File from CCDB in path " << ccdbPaths[iSpecie] << " was not found for run " << bunchCrossing.runNumber();
    }

    TAxis* axis = hSpline->GetXaxis();
    mBetheBlochPiKaPrDe[iSpecie] = {static_cast<double>(hSpline->GetBinContent(axis->FindBin("bb1"))),
                                    static_cast<double>(hSpline->GetBinContent(axis->FindBin("bb2"))),
                                    static_cast<double>(hSpline->GetBinContent(axis->FindBin("bb3"))),
                                    static_cast<double>(hSpline->GetBinContent(axis->FindBin("bb4"))),
                                    static_cast<double>(hSpline->GetBinContent(axis->FindBin("bb5"))),
                                    static_cast<double>(hSpline->GetBinContent(axis->FindBin("Resolution")))};
  }
}

/// load the TPC PID recalibration maps from the CCDB
/// \param ccdb is the CCDB object
/// \param bunchCrossing is the timestamp of bunchcrossing for the run number
/// \param ccdbPath is the path on CCDB for postcalibrations
inline void HfFilterHelper::setTpcRecalibMaps(o2::framework::Service<o2::ccdb::BasicCCDBManager> const& ccdb, aod::BCsWithTimestamps::iterator const& bunchCrossing, const std::string& ccdbPath)
{
  auto calibList = ccdb->getForTimeStamp<TList>(ccdbPath, bunchCrossing.timestamp());
  if (!calibList) {
    LOG(fatal) << "Can not find the TPC Post Calibration object!";
  }
  std::array<std::string, 8> mapNames = {"mean_map_pion", "sigma_map_pion", "mean_map_kaon", "sigma_map_kaon", "mean_map_proton", "sigma_map_proton", "mean_map_deuteron", "sigma_map_deuteron"};

  for (size_t iMap = 0; iMap < mapNames.size(); iMap++) {
    mHistMapPiPrKaDe[iMap] = nullptr;
  }

  for (size_t iMap = 0; iMap < mapNames.size(); iMap++) {

    mHistMapPiPrKaDe[iMap] = reinterpret_cast<TH3F*>(calibList->FindObject(mapNames[iMap].data()));
    if (!mHistMapPiPrKaDe[iMap]) {
      LOG(fatal) << "Cannot find histogram: " << mapNames[iMap].data();
      return;
    }
  }
}

/// Basic selection of proton candidates for Lc
/// \param track is a track
/// \param nsigmaTPCProton max NsigmaTPC for proton candidates
/// \param nsigmaTOFProton max NsigmaTOF for proton candidates
/// \return true if track passes all cuts
template <bool is4beauty, typename T>
inline bool HfFilterHelper::isSelectedProton4CharmOrBeautyBaryons(const T& track)
{
  float NSigmaTPC = track.tpcNSigmaPr();
  float NSigmaTOF = track.tofNSigmaPr();

  if (mTpcPidCalibrationOption == 1) {
    NSigmaTPC = getTPCPostCalib(track, kPr);
  } else if (mTpcPidCalibrationOption == 2) {
    if (track.sign() > 0) {
      NSigmaTPC = getTPCSplineCalib(track, kPr);
    } else {
      NSigmaTPC = getTPCSplineCalib(track, kAntiPr);
    }
  }

  if constexpr (is4beauty) {
    if (std::fabs(NSigmaTPC) > mNSigmaTpcPrKaCutForBeautyToJPsi) {
      return false;
    }
    if (track.hasTOF() && std::fabs(NSigmaTOF) > mNSigmaTofPrKaCutForBeautyToJPsi) {
      return false;
    }
  } else {
    if (std::fabs(NSigmaTPC) > mNSigmaTpcPrCutForCharmBaryons) {
      return false;
    }
    if (track.hasTOF() && std::fabs(NSigmaTOF) > mNSigmaTofPrCutForCharmBaryons) {
      return false;
    }
  }

  return true;
}

/// Basic selection of kaon candidates for kaons from Xic*->SigmaC-Kaon
/// \param isKaonTrack true if we are using a K+- track, false if we are using a K0s (V0)
/// \param track is a track
/// \return true if track passes all cuts
template <bool isKaonTrack, typename T>
inline bool HfFilterHelper::isSelectedKaonFromXicResoToSigmaC(const T& track)
{

  // pt selections
  float pt = track.pt();
  if (pt < mPtMinSoftKaonForXicResoToSigmaC || pt > mPtMaxSoftKaonForXicResoToSigmaC) {
    return false;
  }

  if constexpr (isKaonTrack) {
    /// if the kaon is a track, and not a K0s (V0), check the PID as well
    return isSelectedKaon4Charm3ProngOrBeautyToJPsi(track);
  }

  return true;
}

/// Basic selection of kaon candidates for charm candidates
/// \param track is a track
/// \return true if track passes all cuts
template <bool is4beauty, typename T>
inline bool HfFilterHelper::isSelectedKaon4Charm3ProngOrBeautyToJPsi(const T& track)
{
  float NSigmaTPC = track.tpcNSigmaKa();
  float NSigmaTOF = track.tofNSigmaKa();

  if (mTpcPidCalibrationOption == 1) {
    NSigmaTPC = getTPCPostCalib(track, kKa);
  } else if (mTpcPidCalibrationOption == 2) {
    if (track.sign() > 0) {
      NSigmaTPC = getTPCSplineCalib(track, kKa);
    } else {
      NSigmaTPC = getTPCSplineCalib(track, kAntiKa);
    }
  }

  if constexpr (is4beauty) {
    if (std::fabs(NSigmaTPC) > mNSigmaTpcPrKaCutForBeautyToJPsi) {
      return false;
    }
    if (track.hasTOF() && std::fabs(NSigmaTOF) > mNSigmaTofPrKaCutForBeautyToJPsi) {
      return false;
    }
  } else {
    if (std::fabs(NSigmaTPC) > mNSigmaTpcKaCutFor3Prongs) {
      return false;
    }
    if (track.hasTOF() && std::fabs(NSigmaTOF) > mNSigmaTofKaCutFor3Prongs) {
      return false;
    }
  }

  return true;
}

/// Basic selection of proton candidates forLc and ThetaC decays
/// \param track is a track
/// \return true if track passes all cuts
template <bool is4ThetaC, typename T>
inline bool HfFilterHelper::isSelectedProtonFromLcResoOrThetaC(const T& track)
{

  // pt selections
  float pt = track.pt();
  if constexpr (is4ThetaC) {
    if (pt < mPtMinThetaCBachelor || pt > mPtMaxThetaCBachelor) {
      return false;
    }
  } else {
    if (pt < mPtMinLcResonanceBachelor || pt > mPtMaxLcResonanceBachelor) {
      return false;
    }
  }

  return true;
}

/// Method to perform selections for B+ candidates after vertex reconstruction
/// \param pVecTrack0 is the array for the candidate D daughter momentum after reconstruction of secondary vertex
/// \param pVecTrack1 is the array for the candidate bachelor pion momentum after reconstruction of secondary vertex
/// \param dcaTrack0  is the dca of the D daughter track
/// \param dcaTrack1  is the dca of the pion daughter track
/// \param primVtx is the primary vertex
/// \param secVtx is the secondary vertex
/// \param whichB is the B-hadron species
/// \return true if the beauty candidate passes all cuts
template <typename T1, typename T2, typename T3, typename T4>
inline bool HfFilterHelper::isSelectedBhadron(T1 const& pVecTrack0, T1 const& pVecTrack1, T2 const& dcaTrack0, T2 const& dcaTrack1, const T3& primVtx, const T4& secVtx, const int whichB)
{
  if (whichB == kB0toDStar) {
    LOGP(fatal, "Wrong function used for selection of B0 -> D*pi, please use isSelectedBzeroToDstar");
  }

  auto pVecB = RecoDecay::pVec(pVecTrack0, pVecTrack1);
  auto pTB = RecoDecay::pt(pVecB);
  auto binPtB = findBin(mPtBinsBeautyHadrons, pTB);
  if (binPtB == -1) {
    return false;
  }
  auto cpa = RecoDecay::cpa(primVtx, secVtx, pVecB);
  auto decayLength = RecoDecay::distance(primVtx, secVtx);
  auto impactParameterProduct = dcaTrack0[0] * dcaTrack1[0];

  if (cpa < mCutsBhad[whichB].get(binPtB, 1u)) {
    return false;
  }
  if (decayLength < mCutsBhad[whichB].get(binPtB, 2u)) {
    return false;
  }
  if (impactParameterProduct > mCutsBhad[whichB].get(binPtB, 3u)) {
    return false;
  }

  return true;
}

/// Method to perform selections for B+ candidates after vertex reconstruction
/// \param pVecTrack0 is the array for the candidate D daughter momentum after reconstruction of secondary vertex
/// \param pVecTrack1 is the array for the soft pion momentum after reconstruction of secondary vertex
/// \param pVecTrack2 is the array for the candidate bachelor pion momentum after reconstruction of secondary vertex
/// \param primVtx is the primary vertex
/// \param secVtx is the secondary vertex
/// \return true if the beauty candidate passes all cuts
template <typename T1, typename T2, typename T3>
inline bool HfFilterHelper::isSelectedBzeroToDstar(T1 const& pVecTrack0, T1 const& pVecTrack1, T1 const& pVecTrack2, const T2& primVtx, const T3& secVtx)
{
  auto pVecB = RecoDecay::pVec(pVecTrack0, pVecTrack1, pVecTrack2);
  auto pTB = RecoDecay::pt(pVecB);
  auto binPtB = findBin(mPtBinsBeautyHadrons, pTB);
  if (binPtB == -1) {
    return false;
  }
  auto cpa = RecoDecay::cpa(primVtx, secVtx, pVecB);
  auto decayLength = RecoDecay::distance(primVtx, secVtx);

  if (cpa < mCutsBhad[kB0toDStar].get(binPtB, 1u)) {
    return false;
  }
  if (decayLength < mCutsBhad[kB0toDStar].get(binPtB, 2u)) {
    return false;
  }

  return true;
}

/// Method to perform selections for B+ candidates after vertex reconstruction
/// \param ptCand is the pT of the beauty candidate
/// \param massCand is the mass of the beauty candidate
/// \param whichB is the B-hadron species
/// \return true if the beauty candidate passes all cuts
template <typename T1, typename T2>
inline bool HfFilterHelper::isSelectedBhadronInMassRange(T1 const& ptCand, T2 const& massCand, const int whichB)
{
  auto binPtB = findBin(mPtBinsBeautyHadrons, ptCand);
  if (binPtB == -1) {
    return false;
  }

  float massBhad{-1};
  switch (whichB) {
    case kBplus: {
      massBhad = massBPlus;
      break;
    }
    case kB0toDStar: {
      massBhad = massB0;
      break;
    }
    case kB0: {
      massBhad = massB0;
      break;
    }
    case kBs: {
      massBhad = massBs;
      break;
    }
    case kBc: {
      massBhad = massBc;
      break;
    }
    case kLb: {
      massBhad = massLb;
      break;
    }
    case kXib: {
      massBhad = massXib;
      break;
    }
  }

  if (std::fabs(massCand - massBhad) > mCutsBhad[whichB].get(binPtB, 0u)) {
    return false;
  }

  return true;
}

/// Method to perform selections for B -> JPsiX candidates after vertex reconstruction
/// \param pVecDauTracks is the array of momentum vectors of all daughter tracks
/// \param tracksDauNoMu is the array of tracks for the daughters that are no muons
/// \param primVtx is the primary vertex
/// \param secVtx is the secondary vertex
/// \param activateQA is the flag to enable the
/// \param hMassVsPt is the array of histograms for QA
/// \return true if the beauty candidate passes all cuts
template <int Nprongs, typename T1, typename T2, typename T3, typename T4, typename H2>
inline int8_t HfFilterHelper::isSelectedBhadronToJPsi(std::array<T1, Nprongs> pVecDauTracks, std::array<T2, Nprongs - 2> tracksDauNoMu, const T3& primVtx, const T4& secVtx, const int& activateQA, std::array<H2, nTotBeautyParts>& hMassVsPt)
{
  int8_t isSelected{0};

  auto pVecJPsi = RecoDecay::pVec(pVecDauTracks[0], pVecDauTracks[1]);
  const int offset = static_cast<int>(kNBeautyParticles);

  if constexpr (Nprongs == 3) {
    auto pVecBhad = RecoDecay::pVec(pVecDauTracks[0], pVecDauTracks[1], pVecDauTracks[2]);
    auto ptBhad = RecoDecay::pt(pVecBhad);
    auto binPtB = findBin(mPtBinsBeautyHadrons, ptBhad);
    if (binPtB == -1) {
      return isSelected;
    }
    auto ptMu1 = RecoDecay::pt(pVecDauTracks[0]);
    auto ptMu2 = RecoDecay::pt(pVecDauTracks[1]);
    if (ptMu1 < mCutsBhadToJPsi.get(binPtB, 0u) || ptMu2 < mCutsBhadToJPsi.get(binPtB, 0u)) {
      return isSelected;
    }

    if (RecoDecay::cpa(primVtx, secVtx, pVecBhad) < mCutsBhadToJPsi.get(binPtB, 2u)) {
      return isSelected;
    }

    if (RecoDecay::distance(primVtx, secVtx) < mCutsBhadToJPsi.get(binPtB, 3u)) {
      return isSelected;
    }

    if (isSelectedKaon4Charm3ProngOrBeautyToJPsi<true>(tracksDauNoMu[0])) {
      auto massJPsiKa = RecoDecay::m(std::array{pVecJPsi, pVecDauTracks[2]}, std::array{massJPsi, massKa});
      if (std::fabs(massJPsiKa - massBPlus) < mCutsBhadToJPsi.get(binPtB, 1u)) {
        SETBIT(isSelected, kBplusToJPsi);
        if (activateQA) {
          hMassVsPt[offset + kBplusToJPsi]->Fill(ptBhad, massJPsiKa);
        }
      }
    }
    auto massJPsiPi = RecoDecay::m(std::array{pVecJPsi, pVecDauTracks[2]}, std::array{massJPsi, massPi});
    if (std::fabs(massJPsiPi - massBc) < mCutsBhadToJPsi.get(binPtB, 1u)) {
      SETBIT(isSelected, kBcToJPsi);
      if (activateQA) {
        hMassVsPt[offset + kBcToJPsi]->Fill(ptBhad, massJPsiPi);
      }
    }
  } else if constexpr (Nprongs == 4) {
    auto pVecBhad = RecoDecay::pVec(pVecDauTracks[0], pVecDauTracks[1], pVecDauTracks[2], pVecDauTracks[3]);
    auto ptBhad = RecoDecay::pt(pVecBhad);
    auto binPtB = findBin(mPtBinsBeautyHadrons, ptBhad);
    if (binPtB == -1) {
      return isSelected;
    }
    auto ptMu1 = RecoDecay::pt(pVecDauTracks[0]);
    auto ptMu2 = RecoDecay::pt(pVecDauTracks[1]);
    if (ptMu1 < mCutsBhadToJPsi.get(binPtB, 0u) || ptMu2 < mCutsBhadToJPsi.get(binPtB, 0u)) {
      return isSelected;
    }

    if (RecoDecay::cpa(primVtx, secVtx, pVecBhad) < mCutsBhadToJPsi.get(binPtB, 2u)) {
      return isSelected;
    }

    if (RecoDecay::distance(primVtx, secVtx) < mCutsBhadToJPsi.get(binPtB, 3u)) {
      return isSelected;
    }

    bool isFirstKaon = isSelectedKaon4Charm3ProngOrBeautyToJPsi<true>(tracksDauNoMu[0]);
    bool isSeconKaon = isSelectedKaon4Charm3ProngOrBeautyToJPsi<true>(tracksDauNoMu[1]);
    bool isFirstProton = isSelectedProton4CharmOrBeautyBaryons<true>(tracksDauNoMu[0]);
    bool isSecondProton = isSelectedProton4CharmOrBeautyBaryons<true>(tracksDauNoMu[1]);
    auto massKaKa = RecoDecay::m(std::array{pVecDauTracks[2], pVecDauTracks[3]}, std::array{massKa, massKa});
    if (isFirstKaon && isSeconKaon) {
      if (std::fabs(massKaKa - massPhi) < mCutsBhadToJPsi.get(binPtB, 4u)) {
        auto massJPsiKaKa = RecoDecay::m(std::array{pVecJPsi, pVecDauTracks[2], pVecDauTracks[3]}, std::array{massJPsi, massKa, massKa});
        if (std::fabs(massJPsiKaKa - massBs) < mCutsBhadToJPsi.get(binPtB, 1u)) {
          SETBIT(isSelected, kBsToJPsi);
          if (activateQA) {
            hMassVsPt[offset + kBsToJPsi]->Fill(ptBhad, massJPsiKaKa);
          }
        }
      }
    }
    if (isFirstKaon) {
      auto massKaPi = RecoDecay::m(std::array{pVecDauTracks[2], pVecDauTracks[3]}, std::array{massKa, massPi});
      if (std::fabs(massKaPi - massK0Star892) < mCutsBhadToJPsi.get(binPtB, 5u)) {
        auto massJPsiKaPi = RecoDecay::m(std::array{pVecJPsi, pVecDauTracks[2], pVecDauTracks[3]}, std::array{massJPsi, massKa, massPi});
        if (std::fabs(massJPsiKaPi - massB0) < mCutsBhadToJPsi.get(binPtB, 1u)) {
          SETBIT(isSelected, kB0ToJPsi);
          if (activateQA) {
            hMassVsPt[offset + kB0ToJPsi]->Fill(ptBhad, massJPsiKaPi);
          }
        }
      }
    }
    if (isSeconKaon) {
      auto massPiKa = RecoDecay::m(std::array{pVecDauTracks[2], pVecDauTracks[3]}, std::array{massPi, massKa});
      if (std::fabs(massPiKa - massK0Star892) < mCutsBhadToJPsi.get(binPtB, 5u)) {
        auto massJPsiPiKa = RecoDecay::m(std::array{pVecJPsi, pVecDauTracks[2], pVecDauTracks[3]}, std::array{massJPsi, massPi, massKa});
        if (std::fabs(massJPsiPiKa - massB0) < mCutsBhadToJPsi.get(binPtB, 1u)) {
          SETBIT(isSelected, kB0ToJPsi);
          if (activateQA) {
            hMassVsPt[offset + kB0ToJPsi]->Fill(ptBhad, massJPsiPiKa);
          }
        }
      }
    }
    if (isFirstProton && isSeconKaon) {
      auto massLbToJPsiPrKa = RecoDecay::m(std::array{pVecDauTracks[0], pVecDauTracks[1], pVecDauTracks[2], pVecDauTracks[3]}, std::array{massMu, massMu, massProton, massKa});
      if (std::fabs(massLbToJPsiPrKa - massLb) < mCutsBhadToJPsi.get(binPtB, 1u)) {
        SETBIT(isSelected, kLbToJPsi);
        if (activateQA) {
          hMassVsPt[offset + kLbToJPsi]->Fill(ptBhad, massLbToJPsiPrKa);
        }
      }
    }
    if (isFirstKaon && isSecondProton) {
      auto massLbToJPsiKaPr = RecoDecay::m(std::array{pVecDauTracks[0], pVecDauTracks[1], pVecDauTracks[2], pVecDauTracks[3]}, std::array{massMu, massMu, massKa, massProton});
      if (std::fabs(massLbToJPsiKaPr - massLb) < mCutsBhadToJPsi.get(binPtB, 1u)) {
        SETBIT(isSelected, kLbToJPsi);
        if (activateQA) {
          hMassVsPt[offset + kLbToJPsi]->Fill(ptBhad, massLbToJPsiKaPr);
        }
      }
    }
  }

  return isSelected;
}

/// Method to check if charm candidates has mass between sideband limits
/// \param massHypo1 is the array for the candidate D daughter momentum after reconstruction of secondary vertex
/// \param massHypo2 is the array for the candidate bachelor pion momentum after reconstruction of secondary vertex
/// \param lowLimitSB  is the dca of the D daughter track
/// \param upLimitSB  is the dca of the pion daughter track
/// \return true if the candidate passes the mass selection.
template <typename T1>
inline bool HfFilterHelper::isCharmHadronMassInSbRegions(T1 const& massHypo1, T1 const& massHypo2, const float& lowLimitSB, const float& upLimitSB)
{

  if ((massHypo1 < lowLimitSB || massHypo1 > upLimitSB) && (massHypo2 < lowLimitSB || massHypo2 > upLimitSB)) {
    return false;
  }

  return true;
}

/// Method to check if charm candidates has mass between sideband limits
/// \param trackParCasc is the cascade track parametrisation
/// \param trackParBachelor is the bachelor track parametrisation
/// \param isSelBachelor flag for bachelor selection (Pi/Ka)
/// \param collision is the collision containing the candidate
/// \param dcaFitter is the DCAFitter
/// \param activateQA is the flag to activate the QA
/// \param hMassVsPtXiPi is the 2D histogram with pT vs mass(XiPi)
/// \param hMassVsPtXiKa is the 2D histogram with pT vs mass(XiKa)
template <typename T, typename C, typename H2>
inline bool HfFilterHelper::isSelectedXiBach(T const& trackParCasc, T const& trackParBachelor, int8_t isSelBachelor, C const& collision, o2::vertexing::DCAFitterN<2>& dcaFitter, const int& activateQA, H2 hMassVsPtXiPi, H2 hMassVsPtXiKa)
{
  bool isSelectedXiPi{false}, isSelectedXiKa{false};

  // compute pT
  std::array<float, 3> pVecBachelor{}, pVecCascade{};
  getPxPyPz(trackParBachelor, pVecBachelor);
  getPxPyPz(trackParCasc, pVecCascade);
  auto ptXiBach = RecoDecay::pt(RecoDecay::pVec(pVecCascade, pVecBachelor));

  // compute first mass hypo
  float massXiPi{0.f};
  if (TESTBIT(isSelBachelor, kPionForCharmBaryon)) {
    massXiPi = RecoDecay::m(std::array{pVecCascade, pVecBachelor}, std::array{massXi, massPi});
    if (ptXiBach >= mPtMinXiBach[0] && massXiPi >= mMassMinXiBach[0] && massXiPi <= mMassMaxXiBach[0]) {
      isSelectedXiPi = true;
    }
  }

  // compute second mass hypo
  float massXiKa{0.f};
  if (TESTBIT(isSelBachelor, kKaonForCharmBaryon)) {
    massXiKa = RecoDecay::m(std::array{pVecCascade, pVecBachelor}, std::array{massXi, massKa});
    if (ptXiBach >= mPtMinXiBach[1] && massXiKa >= mMassMinXiBach[1] && massXiKa <= mMassMaxXiBach[1]) {
      isSelectedXiKa = true;
    }
  }

  bool isSelected = isSelectedXiPi || isSelectedXiKa;

  if (isSelected && (mCosPaMinXiBach[0] > -1.f || mChi2PcaMaxXiBach[0] < 1.e6 || mDecLenMinXiBach[0] > 0.)) { // if selected by pT and mass, check topology if applicable
    int nCand = 0;
    try {
      nCand = dcaFitter.process(trackParCasc, trackParBachelor);
    } catch (...) {
      LOG(error) << "Exception caught in DCA fitter process call for Xi + bachelor!";
      return false;
    }
    if (nCand == 0) {
      return false;
    }

    const auto& vtx = dcaFitter.getPCACandidate();
    const auto& trackCascProp = dcaFitter.getTrack(0);
    const auto& trackBachProp = dcaFitter.getTrack(1);
    std::array<float, 3> momCasc{}, momBach{};
    trackCascProp.getPxPyPzGlo(momCasc);
    trackBachProp.getPxPyPzGlo(momBach);
    auto momXiBach = RecoDecay::pVec(momCasc, momBach);

    if (dcaFitter.getChi2AtPCACandidate() > mChi2PcaMaxXiBach[0]) {
      return false;
    }

    std::array<float, 3> primVtx = {collision.posX(), collision.posY(), collision.posZ()};
    if (RecoDecay::cpa(primVtx, std::array{vtx[0], vtx[1], vtx[2]}, momXiBach) < mCosPaMinXiBach[0]) {
      return false;
    }

    if (RecoDecay::distance(primVtx, vtx) < mDecLenMinXiBach[0]) {
      return false;
    }

    if (activateQA) {
      if (isSelectedXiPi) {
        hMassVsPtXiPi->Fill(ptXiBach, massXiPi);
      }
      if (isSelectedXiKa) {
        hMassVsPtXiKa->Fill(ptXiBach, massXiKa);
      }
    }
  }

  return isSelected;
}

/// Method to check if charm candidates has mass between sideband limits
/// \param trackParCasc is the cascade track parametrisation
/// \param trackParBachelor is the array with two bachelor track parametrisations
/// \param collision is the collision containing the candidate
/// \param dcaFitter is the DCAFitter
/// \param activateQA is the flag to activate the QA
/// \param hMassVsPtXiPiPi is the 2D histogram with pT vs mass(XiPiPi)
template <int Nprongs, typename T, typename C, typename H2>
inline bool HfFilterHelper::isSelectedXiBachBach(T const& trackParCasc, std::array<T, 2> const& trackParBachelor, C const& collision, o2::vertexing::DCAFitterN<Nprongs>& dcaFitter, const int& activateQA, H2 hMassVsPtXiPiPi)
{
  // compute pT
  std::array<float, 3> pVecBachelorFirst{}, pVecBachelorSecond{}, pVecCascade{};
  getPxPyPz(trackParBachelor[0], pVecBachelorFirst);
  getPxPyPz(trackParBachelor[1], pVecBachelorSecond);
  getPxPyPz(trackParCasc, pVecCascade);
  auto ptXiBachBach = RecoDecay::pt(RecoDecay::pVec(pVecCascade, pVecBachelorFirst, pVecBachelorSecond));
  if (ptXiBachBach < mPtMinXiBach[2]) {
    return false;
  }

  // compute mass
  float massXiPiPi = RecoDecay::m(std::array{pVecCascade, pVecBachelorFirst, pVecBachelorSecond}, std::array{massXi, massPi, massPi});
  if (massXiPiPi < mMassMinXiBach[2] || massXiPiPi > mMassMaxXiBach[2]) {
    return false;
  }

  if (mCosPaMinXiBach[1] > -1.f || mChi2PcaMaxXiBach[1] < 1.e6 || mDecLenMinXiBach[0] > 0.) { // check topology if applicable
    int nCand = 0;
    if constexpr (Nprongs == 3) {
      try {
        nCand = dcaFitter.process(trackParBachelor[0], trackParBachelor[1], trackParCasc);
      } catch (...) {
        LOG(error) << "Exception caught in DCA fitter process call for Xi + bachelor + bachelor!";
        return false;
      }
      if (nCand == 0) {
        return false;
      }
    } else if constexpr (Nprongs == 2) {
      try {
        nCand = dcaFitter.process(trackParBachelor[0], trackParBachelor[1]);
      } catch (...) {
        LOG(error) << "Exception caught in DCA fitter process call for Xi + bachelor + bachelor!";
        return false;
      }
      if (nCand == 0) {
        return false;
      }
    }

    const auto& vtx = dcaFitter.getPCACandidate();

    std::array<float, 3> momCasc{pVecCascade}, momBachFirst{}, momBachSecond{};
    const auto& trackBachFirstProp = dcaFitter.getTrack(0);
    const auto& trackBachSecondProp = dcaFitter.getTrack(1);
    trackBachFirstProp.getPxPyPzGlo(momBachFirst);
    trackBachSecondProp.getPxPyPzGlo(momBachSecond);
    if constexpr (Nprongs == 3) {
      const auto& trackCascProp = dcaFitter.getTrack(2);
      trackCascProp.getPxPyPzGlo(momCasc);
    }
    auto momXiBachBach = RecoDecay::pVec(momCasc, momBachFirst, momBachSecond);

    if (dcaFitter.getChi2AtPCACandidate() > mChi2PcaMaxXiBach[1]) {
      return false;
    }

    std::array<float, 3> primVtx = {collision.posX(), collision.posY(), collision.posZ()};
    if (RecoDecay::cpa(primVtx, std::array{vtx[0], vtx[1], vtx[2]}, momXiBachBach) < mCosPaMinXiBach[1]) {
      return false;
    }

    if (RecoDecay::distance(primVtx, vtx) < mDecLenMinXiBach[1]) {
      return false;
    }

    if (activateQA) {
      hMassVsPtXiPiPi->Fill(ptXiBachBach, massXiPiPi);
    }
  }

  return true;
}

/// Update the TPC PID based on the spline of particles
/// \param track is a track parameter
/// \param pidSpecies is the particle species to be considered
/// \return updated nsigma value for TPC PID
template <typename T>
inline float HfFilterHelper::getTPCSplineCalib(const T& track, const int& pidSpecies)
{
  float tpcPin = track.tpcInnerParam();
  float dEdx = track.tpcSignal();

  return getTPCSplineCalib(tpcPin, dEdx, pidSpecies);
}

/// Update the TPC PID baesd on the spline of particles
/// \param tpcPin is the TPC momentum at innermost update
/// \param dEdx is the TPC dEdx
/// \param pidSpecies is the particle species to be considered
/// \return updated nsigma value for TPC PID
inline float HfFilterHelper::getTPCSplineCalib(const float tpcPin, const float dEdx, const int& pidSpecies)
{
  float mMassPar{0.};
  if (pidSpecies == kPi || pidSpecies == kAntiPi) {
    mMassPar = massPi;
  } else if (pidSpecies == kKa || pidSpecies == kAntiKa) {
    mMassPar = massKa;
  } else if (pidSpecies == kPr || pidSpecies == kAntiPr) {
    mMassPar = massProton;
  } else if (pidSpecies == kDe || pidSpecies == kAntiDe) {
    mMassPar = massDeuteron;
  } else {
    LOGP(fatal, "TPC recalibrated Nsigma requested for unknown particle species, return 999");
    return 999.;
  }

  auto bgScaling = 1 / mMassPar;
  double expBethe = common::BetheBlochAleph(static_cast<double>(tpcPin * bgScaling), mBetheBlochPiKaPrDe[pidSpecies][0], mBetheBlochPiKaPrDe[pidSpecies][1], mBetheBlochPiKaPrDe[pidSpecies][2], mBetheBlochPiKaPrDe[pidSpecies][3], mBetheBlochPiKaPrDe[pidSpecies][4]);
  double expSigma = expBethe * mBetheBlochPiKaPrDe[pidSpecies][5];
  return static_cast<float>((dEdx - expBethe) / expSigma);
}

/// compute TPC postcalibrated nsigma based on calibration histograms from CCDB
/// \param tpcPin is the TPC momentum at innermost update
/// \param tpcNCls is the number of found TPC clusters
/// \param eta is the pseudorapidity
/// \param tpcNSigma is the original Nsigma
/// \param pidSpecies is the PID species
/// \return the corrected Nsigma value for the PID species
inline float HfFilterHelper::getTPCPostCalib(const float tpcPin, const float tpcNCls, const float eta, const float tpcNSigma, const int& pidSpecies)
{
  int iHist{0};
  if (pidSpecies == kPi) {
    iHist = 0;
  } else if (pidSpecies == kKa) {
    iHist = 2;
  } else if (pidSpecies == kPr) {
    iHist = 4;
  } else {
    LOG(fatal) << "Wrong PID Species be selected, please check!";
  }

  if (!mHistMapPiPrKaDe[iHist] || !mHistMapPiPrKaDe[iHist + 1]) {
    LOGP(warn, "Postcalibration TPC PID histograms not set. Use default Nsigma values.");
  }

  auto binTPCNCls = mHistMapPiPrKaDe[iHist]->GetXaxis()->FindBin(tpcNCls);
  binTPCNCls = (binTPCNCls == 0 ? 1 : binTPCNCls);
  binTPCNCls = std::min(mHistMapPiPrKaDe[iHist]->GetXaxis()->GetNbins(), binTPCNCls);
  auto binPin = mHistMapPiPrKaDe[iHist]->GetYaxis()->FindBin(tpcPin);
  binPin = (binPin == 0 ? 1 : binPin);
  binPin = std::min(mHistMapPiPrKaDe[iHist]->GetYaxis()->GetNbins(), binPin);
  auto binEta = mHistMapPiPrKaDe[iHist]->GetZaxis()->FindBin(eta);
  binEta = (binEta == 0 ? 1 : binEta);
  binEta = std::min(mHistMapPiPrKaDe[iHist]->GetZaxis()->GetNbins(), binEta);

  auto mean = mHistMapPiPrKaDe[iHist]->GetBinContent(binTPCNCls, binPin, binEta);
  auto width = mHistMapPiPrKaDe[iHist + 1]->GetBinContent(binTPCNCls, binPin, binEta);

  return (tpcNSigma - mean) / width;
}

/// compute TPC postcalibrated nsigma based on calibration histograms from CCDB
/// \param track is the track
/// \param pidSpecies is the PID species
/// \return the corrected Nsigma value for the PID species
template <typename T>
inline float HfFilterHelper::getTPCPostCalib(const T& track, const int& pidSpecies)
{
  float tpcNCls = track.tpcNClsFound();
  float tpcPin = track.tpcInnerParam();
  float eta = track.eta();
  float tpcNSigma{-999.};

  if (pidSpecies == kPi) {
    tpcNSigma = track.tpcNSigmaPi();
  } else if (pidSpecies == kKa) {
    tpcNSigma = track.tpcNSigmaKa();
  } else if (pidSpecies == kPr) {
    tpcNSigma = track.tpcNSigmaPr();
  } else {
    LOG(fatal) << "Wrong PID Species be selected, please check!";
  }

  return getTPCPostCalib(tpcPin, tpcNCls, eta, tpcNSigma, pidSpecies);
}

/// Finds pT bin in an array.
/// \param bins  array of pT bins
/// \param value  pT
/// \return index of the pT bin
/// \note Accounts for the offset so that pt bin array can be used to also configure a histogram axis.
template <typename T1, typename T2>
inline int HfFilterHelper::findBin(T1 const& binsPt, T2 value)
{
  if (value < binsPt.front()) {
    return -1;
  }
  if (value >= binsPt.back()) {
    return -1;
  }
  return std::distance(binsPt.begin(), std::upper_bound(binsPt.begin(), binsPt.end(), value)) - 1;
}

/// Set vertxing configuration
/// \param vertexer o2::vertexing::DCAFitterN<N> object
template <typename T1>
inline int HfFilterHelper::setVtxConfiguration(T1& vertexer, bool useAbsDCA)
{
  // Fitter initialisation
  vertexer.setPropagateToPCA(true);
  vertexer.setMaxR(200.);
  vertexer.setMaxDZIni(1.e9);
  vertexer.setMaxDXYIni(4.);
  vertexer.setMinParamChange(1.e-3);
  vertexer.setMinRelChi2Change(0.9);
  vertexer.setMaxChi2(0.9);
  vertexer.setUseAbsDCA(useAbsDCA);
  vertexer.setWeightedFinalPCA(false);
  return 1;
}

/// Utility to compute AP alpha and qt
/// \param momPos momentum array of positive daughter
/// \param momNeg momentum array of negative daughter
template <typename T>
inline std::array<T, 2> HfFilterHelper::alphaAndQtAP(std::array<T, 3> const& momPos, std::array<T, 3> const& momNeg)
{
  float momTot = RecoDecay::p(momPos[0] + momNeg[0], momPos[1] + momNeg[1], momPos[2] + momNeg[2]);
  float lQlNeg = RecoDecay::dotProd(momNeg, std::array{momPos[0] + momNeg[0], momPos[1] + momNeg[1], momPos[2] + momNeg[2]}) / momTot;
  float lQlPos = RecoDecay::dotProd(momPos, std::array{momPos[0] + momNeg[0], momPos[1] + momNeg[1], momPos[2] + momNeg[2]}) / momTot;
  float alpha = (lQlPos - lQlNeg) / (lQlPos + lQlNeg);
  float qtarm = std::sqrt(RecoDecay::p2(momNeg) - lQlNeg * lQlNeg);

  std::array<float, 2> alphaAndQt = {alpha, qtarm};
  return alphaAndQt;
}

/// build V0 candidate from table with track indices
/// \param v0Indices V0 candidate from AO2D table (track indices)
/// \param tracks track table
/// \param collision collision
/// \param dcaFitter DCA fitter to be used
/// \param vetoedTrackIds vector with forbidden track indices, if any
template <typename V, typename T, typename C>
inline bool HfFilterHelper::buildV0(V const& v0Indices, T const& tracks, C const& collision, o2::vertexing::DCAFitterN<2>& dcaFitter, const std::vector<int>& vetoedTrackIds, V0Cand& v0Cand)
{
  auto trackPos = tracks.rawIteratorAt(v0Indices.posTrackId());
  auto trackNeg = tracks.rawIteratorAt(v0Indices.negTrackId());

  // minimal track cuts
  if (!trackPos.hasTPC() || !trackNeg.hasTPC()) {
    return false;
  }

  if (trackPos.tpcNClsCrossedRows() < 50 || trackNeg.tpcNClsCrossedRows() < 50) {
    return false;
  }

  if (std::find(vetoedTrackIds.begin(), vetoedTrackIds.end(), trackPos.globalIndex()) != vetoedTrackIds.end() || std::find(vetoedTrackIds.begin(), vetoedTrackIds.end(), trackNeg.globalIndex()) != vetoedTrackIds.end()) {
    return false;
  }

  auto trackParCovPos = getTrackParCov(trackPos);
  auto trackParCovNeg = getTrackParCov(trackNeg);
  std::array<float, 3> primVtx = {collision.posX(), collision.posY(), collision.posZ()};
  std::array<float, 2> dcaInfoPos, dcaInfoNeg;
  o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParCovPos, 2.f, dcaFitter.getMatCorrType(), &dcaInfoPos);
  o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParCovNeg, 2.f, dcaFitter.getMatCorrType(), &dcaInfoNeg);

  // reconstruct vertex
  int nCand = 0;
  try {
    nCand = dcaFitter.process(trackParCovPos, trackParCovNeg);
  } catch (...) {
    LOG(error) << "Exception caught in DCA fitter process call in V0!";
    return false;
  }
  if (nCand == 0) {
    return false;
  }

  // compute candidate momentum from tracks propagated to decay vertex
  auto& trackPosProp = dcaFitter.getTrack(0);
  auto& trackNegProp = dcaFitter.getTrack(1);
  std::array<float, 3> momPos{}, momNeg{};
  trackPosProp.getPxPyPzGlo(momPos);
  trackNegProp.getPxPyPzGlo(momNeg);
  v0Cand.mom = RecoDecay::pVec(momPos, momNeg);

  // fill V0 quantities
  v0Cand.dcapostopv = dcaInfoPos[0];
  v0Cand.dcanegtopv = dcaInfoNeg[0];
  v0Cand.ptPos = RecoDecay::pt(momPos);
  v0Cand.ptNeg = RecoDecay::pt(momNeg);
  v0Cand.pinTpcPos = trackPos.tpcInnerParam();
  v0Cand.pinTpcNeg = trackNeg.tpcInnerParam();
  v0Cand.nClsFoundTpcPos = trackPos.tpcNClsFound();
  v0Cand.nClsFoundTpcNeg = trackNeg.tpcNClsFound();
  v0Cand.nClsCrossedRowsTpcPos = trackPos.tpcNClsCrossedRows();
  v0Cand.nClsCrossedRowsTpcNeg = trackNeg.tpcNClsCrossedRows();
  v0Cand.crossedRowsOverFindableClsTpcPos = trackPos.tpcCrossedRowsOverFindableCls();
  v0Cand.crossedRowsOverFindableClsTpcNeg = trackNeg.tpcCrossedRowsOverFindableCls();
  v0Cand.signalTpcPos = trackPos.tpcSignal();
  v0Cand.signalTpcNeg = trackNeg.tpcSignal();
  v0Cand.etaPos = RecoDecay::eta(momPos);
  v0Cand.etaNeg = RecoDecay::eta(momNeg);
  v0Cand.dcaV0daughters = std::sqrt(dcaFitter.getChi2AtPCACandidate());

  const auto& vtx = dcaFitter.getPCACandidate();
  for (int iCoord{0}; iCoord < 3; ++iCoord) {
    v0Cand.vtx[iCoord] = vtx[iCoord];
  }
  auto covVtxV = dcaFitter.calcPCACovMatrix(0);
  v0Cand.cov = {};
  v0Cand.cov[0] = covVtxV(0, 0);
  v0Cand.cov[1] = covVtxV(1, 0);
  v0Cand.cov[2] = covVtxV(1, 1);
  v0Cand.cov[3] = covVtxV(2, 0);
  v0Cand.cov[4] = covVtxV(2, 1);
  v0Cand.cov[5] = covVtxV(2, 2);
  std::array<float, 21> covTpositive = {0.};
  std::array<float, 21> covTnegative = {0.};
  trackPosProp.getCovXYZPxPyPzGlo(covTpositive);
  trackNegProp.getCovXYZPxPyPzGlo(covTnegative);
  constexpr int MomInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
  for (int iCoord{0}; iCoord < 6; ++iCoord) {
    v0Cand.cov[MomInd[iCoord]] = covTpositive[MomInd[iCoord]] + covTnegative[MomInd[iCoord]];
  }
  v0Cand.v0radius = std::hypot(vtx[0], vtx[1]);
  v0Cand.v0cosPA = RecoDecay::cpa(primVtx, vtx, v0Cand.mom);

  auto trackParV0 = dcaFitter.createParentTrackParCov();
  trackParV0.setAbsCharge(0);
  trackParV0.setPID(o2::track::PID::K0);
  std::array<float, 2> dcaInfoV0;
  o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParV0, 2.f, dcaFitter.getMatCorrType(), &dcaInfoV0);
  v0Cand.dcav0topv = dcaInfoV0[0];

  v0Cand.mK0Short = RecoDecay::m(std::array{momPos, momNeg}, std::array{massPi, massPi});
  v0Cand.mLambda = RecoDecay::m(std::array{momPos, momNeg}, std::array{massProton, massPi});
  v0Cand.mAntiLambda = RecoDecay::m(std::array{momPos, momNeg}, std::array{massPi, massProton});

  auto alphaAndQt = alphaAndQtAP(momPos, momNeg);
  v0Cand.alpha = alphaAndQt[0];
  v0Cand.qtarm = alphaAndQt[1];

  v0Cand.hasTofPos = trackPos.hasTOF();
  v0Cand.hasTofNeg = trackNeg.hasTOF();
  v0Cand.nSigmaPrTpcPos = trackPos.tpcNSigmaPr();
  v0Cand.nSigmaPrTofPos = trackPos.tofNSigmaPr();
  v0Cand.nSigmaPrTpcNeg = trackNeg.tpcNSigmaPr();
  v0Cand.nSigmaPrTofNeg = trackNeg.tofNSigmaPr();
  v0Cand.nSigmaPiTpcPos = trackPos.tpcNSigmaPi();
  v0Cand.nSigmaPiTofPos = trackPos.tofNSigmaPi();
  v0Cand.nSigmaPiTpcNeg = trackNeg.tpcNSigmaPi();
  v0Cand.nSigmaPiTofNeg = trackNeg.tofNSigmaPi();

  return true;
}

/// build cascade candidate from table with track indices
/// \param cascIndices cascade candidate from AO2D table (track indices)
/// \param v0Indices V0 candidate from AO2D table (track indices)
/// \param tracks track table
/// \param collision collision
/// \param dcaFitter DCA fitter to be used
/// \param vetoedTrackIds vector with forbidden track indices, if any
template <typename Casc, typename T, typename C, typename V>
inline bool HfFilterHelper::buildCascade(Casc const& cascIndices, V const& v0Indices, T const& tracks, C const& collision, o2::vertexing::DCAFitterN<2>& dcaFitter, const std::vector<int>& vetoedTrackIds, CascCand& cascCand)
{
  auto v0 = v0Indices.rawIteratorAt(cascIndices.v0Id());
  auto trackBachelor = tracks.rawIteratorAt(cascIndices.bachelorId());

  // minimal track cuts
  if (!trackBachelor.hasTPC()) {
    return false;
  }

  if (trackBachelor.tpcNClsCrossedRows() < 50) {
    return false;
  }

  if (std::find(vetoedTrackIds.begin(), vetoedTrackIds.end(), trackBachelor.globalIndex()) != vetoedTrackIds.end()) {
    return false;
  }

  std::array<float, 2> dcaInfoBach;
  auto bachTrackParCov = getTrackParCov(trackBachelor);
  o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, bachTrackParCov, 2.f, dcaFitter.getMatCorrType(), &dcaInfoBach);

  // first we build V0 candidate
  V0Cand v0Cand;
  if (!buildV0(v0, tracks, collision, dcaFitter, vetoedTrackIds, v0Cand)) {
    return false;
  }

  // Set up covariance matrices (should in fact be optional)
  auto v0TrackParCov = o2::track::TrackParCov({v0Cand.vtx[0], v0Cand.vtx[1], v0Cand.vtx[2]}, {v0Cand.mom[0], v0Cand.mom[1], v0Cand.mom[2]}, v0Cand.cov, 0, true);
  v0TrackParCov.setAbsCharge(0);
  v0TrackParCov.setPID(o2::track::PID::Lambda);

  int nCand = 0;
  try {
    nCand = dcaFitter.process(v0TrackParCov, bachTrackParCov);
  } catch (...) {
    LOG(error) << "Exception caught in DCA fitter process call in Xi!";
    return false;
  }
  if (nCand == 0) {
    return false;
  }

  // compute candidate momentum from tracks propagated to decay vertex
  auto& trackV0Prop = dcaFitter.getTrack(0);
  auto& trackBachProp = dcaFitter.getTrack(1);
  std::array<float, 3> momV0{}, momBach{};
  trackV0Prop.getPxPyPzGlo(momV0);
  trackBachProp.getPxPyPzGlo(momBach);
  cascCand.mom = RecoDecay::pVec(momV0, momBach);
  cascCand.sign = trackBachelor.sign();

  cascCand.v0 = v0Cand;
  cascCand.ptBach = RecoDecay::pt(momBach);
  cascCand.etaBach = RecoDecay::eta(momBach);
  cascCand.pinTpcBach = trackBachelor.tpcInnerParam();
  cascCand.nClsFoundTpcBach = trackBachelor.tpcNClsFound();
  cascCand.nClsCrossedRowsTpcBach = trackBachelor.tpcNClsCrossedRows();
  cascCand.crossedRowsOverFindableClsTpcBach = trackBachelor.tpcCrossedRowsOverFindableCls();
  cascCand.signalTpcBach = trackBachelor.tpcSignal();
  cascCand.pt = RecoDecay::pt(cascCand.mom);

  std::array<float, 3> primVtx = {collision.posX(), collision.posY(), collision.posZ()};
  const auto& vtx = dcaFitter.getPCACandidate();
  for (int iCoord{0}; iCoord < 3; ++iCoord) {
    cascCand.vtx[iCoord] = vtx[iCoord];
  }
  auto covVtxV = dcaFitter.calcPCACovMatrix(0);
  cascCand.cov = {};
  cascCand.cov[0] = covVtxV(0, 0);
  cascCand.cov[1] = covVtxV(1, 0);
  cascCand.cov[2] = covVtxV(1, 1);
  cascCand.cov[3] = covVtxV(2, 0);
  cascCand.cov[4] = covVtxV(2, 1);
  cascCand.cov[5] = covVtxV(2, 2);
  std::array<float, 21> covTv0 = {0.};
  std::array<float, 21> covTbachelor = {0.};
  trackV0Prop.getCovXYZPxPyPzGlo(covTv0);
  trackBachProp.getCovXYZPxPyPzGlo(covTbachelor);
  constexpr int MomInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
  for (int iCoord{0}; iCoord < 6; ++iCoord) {
    cascCand.cov[MomInd[iCoord]] = covTv0[MomInd[iCoord]] + covTbachelor[MomInd[iCoord]];
  }

  cascCand.cascradius = std::hypot(vtx[0], vtx[1]);
  cascCand.casccosPA = RecoDecay::cpa(primVtx, vtx, cascCand.mom);

  auto trackParCasc = dcaFitter.createParentTrackParCov();
  trackParCasc.setAbsCharge(1);
  trackParCasc.setPID(o2::track::PID::XiMinus);
  std::array<float, 2> dcaInfoCasc;
  o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParCasc, 2.f, dcaFitter.getMatCorrType(), &dcaInfoCasc);
  cascCand.dcaXYCascToPV = dcaInfoCasc[0];
  cascCand.dcacascdaughters = std::sqrt(dcaFitter.getChi2AtPCACandidate());
  cascCand.mXi = RecoDecay::m(std::array{momBach, momV0}, std::array{massPi, massLambda});
  cascCand.mOmega = RecoDecay::m(std::array{momBach, momV0}, std::array{massKa, massLambda});

  cascCand.hasTofBach = trackBachelor.hasTOF();
  cascCand.nSigmaPiTpcBach = trackBachelor.tpcNSigmaPi();
  cascCand.nSigmaPiTofBach = trackBachelor.tofNSigmaPi();

  return true;
}

} // namespace hffilters

/// definition of tables
namespace hftraining
{
DECLARE_SOA_COLUMN(InvMassD0, invMassD0, float);                 //!
DECLARE_SOA_COLUMN(InvMassD0bar, invMassD0bar, float);           //!
DECLARE_SOA_COLUMN(InvMassDplus, invMassDplus, float);           //!
DECLARE_SOA_COLUMN(InvMassDsToKKPi, invMassDsToKKPi, float);     //!
DECLARE_SOA_COLUMN(InvMassDsToPiKK, invMassDsToPiKK, float);     //!
DECLARE_SOA_COLUMN(InvMassLcToPKPi, invMassLcToPKPi, float);     //!
DECLARE_SOA_COLUMN(InvMassLcToPiKP, invMassLcToPiKP, float);     //!
DECLARE_SOA_COLUMN(InvMassXicToPKPi, invMassXicToPKPi, float);   //!
DECLARE_SOA_COLUMN(InvMassXicToPiKP, invMassXicToPiKP, float);   //!
DECLARE_SOA_COLUMN(PT2Prong, pT2Prong, float);                   //!
DECLARE_SOA_COLUMN(PT3Prong, pT3Prong, float);                   //!
DECLARE_SOA_COLUMN(DeltaMassKKFirst, deltaMassKKFirst, float);   //!
DECLARE_SOA_COLUMN(DeltaMassKKSecond, deltaMassKKSecond, float); //!
DECLARE_SOA_COLUMN(PT1, pT1, float);                             //!
DECLARE_SOA_COLUMN(DCAPrimXY1, dcaPrimXY1, float);               //!
DECLARE_SOA_COLUMN(DCAPrimZ1, dcaPrimZ1, float);                 //!
DECLARE_SOA_COLUMN(NsigmaPiTPC1, nsigmaPiTPC1, float);           //!
DECLARE_SOA_COLUMN(NsigmaKaTPC1, nsigmaKaTPC1, float);           //!
DECLARE_SOA_COLUMN(NsigmaPrTPC1, nsigmaPrTPC1, float);           //!
DECLARE_SOA_COLUMN(NsigmaPiTOF1, nsigmaPiTOF1, float);           //!
DECLARE_SOA_COLUMN(NsigmaKaTOF1, nsigmaKaTOF1, float);           //!
DECLARE_SOA_COLUMN(NsigmaPrTOF1, nsigmaPrTOF1, float);           //!
DECLARE_SOA_COLUMN(PT2, pT2, float);                             //!
DECLARE_SOA_COLUMN(DCAPrimXY2, dcaPrimXY2, float);               //!
DECLARE_SOA_COLUMN(DCAPrimZ2, dcaPrimZ2, float);                 //!
DECLARE_SOA_COLUMN(NsigmaPiTPC2, nsigmaPiTPC2, float);           //!
DECLARE_SOA_COLUMN(NsigmaKaTPC2, nsigmaKaTPC2, float);           //!
DECLARE_SOA_COLUMN(NsigmaPrTPC2, nsigmaPrTPC2, float);           //!
DECLARE_SOA_COLUMN(NsigmaPiTOF2, nsigmaPiTOF2, float);           //!
DECLARE_SOA_COLUMN(NsigmaKaTOF2, nsigmaKaTOF2, float);           //!
DECLARE_SOA_COLUMN(NsigmaPrTOF2, nsigmaPrTOF2, float);           //!
DECLARE_SOA_COLUMN(PT3, pT3, float);                             //!
DECLARE_SOA_COLUMN(DCAPrimXY3, dcaPrimXY3, float);               //!
DECLARE_SOA_COLUMN(DCAPrimZ3, dcaPrimZ3, float);                 //!
DECLARE_SOA_COLUMN(NsigmaPiTPC3, nsigmaPiTPC3, float);           //!
DECLARE_SOA_COLUMN(NsigmaKaTPC3, nsigmaKaTPC3, float);           //!
DECLARE_SOA_COLUMN(NsigmaPrTPC3, nsigmaPrTPC3, float);           //!
DECLARE_SOA_COLUMN(NsigmaPiTOF3, nsigmaPiTOF3, float);           //!
DECLARE_SOA_COLUMN(NsigmaKaTOF3, nsigmaKaTOF3, float);           //!
DECLARE_SOA_COLUMN(NsigmaPrTOF3, nsigmaPrTOF3, float);           //!
DECLARE_SOA_COLUMN(FlagOrigin, flagOrigin, int8_t);              //!
DECLARE_SOA_COLUMN(Channel, channel, int8_t);                    //!
DECLARE_SOA_COLUMN(HFSelBit, hfselbit, int8_t);                  //!
DECLARE_SOA_COLUMN(IsInCorrectColl, isInCorrectColl, bool);      //!
} // namespace hftraining

DECLARE_SOA_TABLE(HFTrigTrain2P, "AOD", "HFTRIGTRAIN2P", //!
                  hftraining::InvMassD0,
                  hftraining::InvMassD0bar,
                  hftraining::PT2Prong,
                  hftraining::PT1,
                  hftraining::DCAPrimXY1,
                  hftraining::DCAPrimZ1,
                  hftraining::NsigmaPiTPC1,
                  hftraining::NsigmaKaTPC1,
                  hftraining::NsigmaPiTOF1,
                  hftraining::NsigmaKaTOF1,
                  hftraining::PT2,
                  hftraining::DCAPrimXY2,
                  hftraining::DCAPrimZ2,
                  hftraining::NsigmaPiTPC2,
                  hftraining::NsigmaKaTPC2,
                  hftraining::NsigmaPiTOF2,
                  hftraining::NsigmaKaTOF2,
                  hftraining::FlagOrigin,
                  hftraining::IsInCorrectColl);
DECLARE_SOA_TABLE(HFTrigTrain3P, "AOD", "HFTRIGTRAIN3P", //!
                  hftraining::InvMassDplus,
                  hftraining::InvMassDsToKKPi,
                  hftraining::InvMassDsToPiKK,
                  hftraining::InvMassLcToPKPi,
                  hftraining::InvMassLcToPiKP,
                  hftraining::InvMassXicToPKPi,
                  hftraining::InvMassXicToPiKP,
                  hftraining::PT3Prong,
                  hftraining::DeltaMassKKFirst,
                  hftraining::DeltaMassKKSecond,
                  hftraining::PT1,
                  hftraining::DCAPrimXY1,
                  hftraining::DCAPrimZ1,
                  hftraining::NsigmaPiTPC1,
                  hftraining::NsigmaKaTPC1,
                  hftraining::NsigmaPrTPC1,
                  hftraining::NsigmaPiTOF1,
                  hftraining::NsigmaKaTOF1,
                  hftraining::NsigmaPrTOF1,
                  hftraining::PT2,
                  hftraining::DCAPrimXY2,
                  hftraining::DCAPrimZ2,
                  hftraining::NsigmaPiTPC2,
                  hftraining::NsigmaKaTPC2,
                  hftraining::NsigmaPrTPC2,
                  hftraining::NsigmaPiTOF2,
                  hftraining::NsigmaKaTOF2,
                  hftraining::NsigmaPrTOF2,
                  hftraining::PT3,
                  hftraining::DCAPrimXY3,
                  hftraining::DCAPrimZ3,
                  hftraining::NsigmaPiTPC3,
                  hftraining::NsigmaKaTPC3,
                  hftraining::NsigmaPrTPC3,
                  hftraining::NsigmaPiTOF3,
                  hftraining::NsigmaKaTOF3,
                  hftraining::NsigmaPrTOF3,
                  hftraining::FlagOrigin,
                  hftraining::Channel,
                  hftraining::HFSelBit,
                  hftraining::IsInCorrectColl);

namespace hfoptimisationTree
{
DECLARE_SOA_COLUMN(CollisionIndex, collisionIndex, int); //!
DECLARE_SOA_COLUMN(ParticleID, particleID, int);         //!
DECLARE_SOA_COLUMN(Pt, pt, float);                       //!
DECLARE_SOA_COLUMN(BkgBDT, bkgBDT, float);               //!
DECLARE_SOA_COLUMN(PromptBDT, promptBDT, float);         //!
DECLARE_SOA_COLUMN(NonpromptBDT, nonpromptBDT, float);   //!
DECLARE_SOA_COLUMN(DCAXY, dcaXY, float);                 //!
DECLARE_SOA_COLUMN(KStar, kStar, float);                 //!
DECLARE_SOA_COLUMN(NsigmaPrTPC, nsigmaPrTPC, float);     //!
DECLARE_SOA_COLUMN(NsigmaPrTOF, nsigmaPrTOF, float);     //!
DECLARE_SOA_COLUMN(NsigmaDeTPC, nsigmaDeTPC, float);     //!
DECLARE_SOA_COLUMN(NsigmaDeTOF, nsigmaDeTOF, float);     //!
} // namespace hfoptimisationTree

DECLARE_SOA_TABLE(HFOptimisationTreeBeauty, "AOD", "HFOPTIMTREEB", //!
                  hfoptimisationTree::CollisionIndex,
                  hfoptimisationTree::ParticleID,
                  hfoptimisationTree::Pt,
                  hfoptimisationTree::BkgBDT,
                  hfoptimisationTree::PromptBDT,
                  hfoptimisationTree::NonpromptBDT,
                  hfoptimisationTree::DCAXY);
DECLARE_SOA_TABLE(HFOptimisationTreeCharm, "AOD", "HFOPTIMTREEC", //!
                  hfoptimisationTree::CollisionIndex,
                  hfoptimisationTree::ParticleID,
                  hfoptimisationTree::Pt,
                  hfoptimisationTree::BkgBDT,
                  hfoptimisationTree::PromptBDT,
                  hfoptimisationTree::NonpromptBDT);
DECLARE_SOA_TABLE(HFOptimisationTreeFemto, "AOD", "HFOPTIMTREEF", //!
                  hfoptimisationTree::CollisionIndex,
                  hfoptimisationTree::ParticleID,
                  hfoptimisationTree::Pt,
                  hfoptimisationTree::BkgBDT,
                  hfoptimisationTree::PromptBDT,
                  hfoptimisationTree::NonpromptBDT,
                  hfoptimisationTree::KStar,
                  hfoptimisationTree::NsigmaPrTPC,
                  hfoptimisationTree::NsigmaPrTOF,
                  hfoptimisationTree::NsigmaDeTPC,
                  hfoptimisationTree::NsigmaDeTOF);
DECLARE_SOA_TABLE(HFOptimisationTreeCollisions, "AOD", "HFOPTIMTREECOLL", //!
                  hfoptimisationTree::CollisionIndex)
} // namespace o2::aod

#endif // EVENTFILTERING_PWGHF_HFFILTERHELPERS_H_
