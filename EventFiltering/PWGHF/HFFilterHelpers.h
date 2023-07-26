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
// O2 includes

/// \file HFFilterHelpers.h
/// \brief Header file with definition of variables, methods, and tables used in the HFFilter.cxx task
///
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN
/// \author Marcel Lesch <marcel.lesch@tum.de>, TUM
/// \author Alexandre Bigot <alexandre.bigot@cern.ch>, Strasbourg University
/// \author Biao Zhang <biao.zhang@cern.ch>, CCNU
/// \author Federica Zanone <federica.zanone@cern.ch>, Heidelberg University

#ifndef EVENTFILTERING_PWGHF_HFFILTERHELPERS_H_
#define EVENTFILTERING_PWGHF_HFFILTERHELPERS_H_

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "Math/GenVector/Boost.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

#include "CCDB/CcdbApi.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/DataTypes.h"
#include "Framework/HistogramRegistry.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

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
  kB0,
  kBs,
  kLb,
  kXib,
  kNBeautyParticles
};

enum bachelorTrackSelection {
  kRejected = 0,
  kSoftPion,
  kForBeauty,
  kSoftPionForBeauty,
  kPionForCharmBaryon,
  kKaonForCharmBaryon
};

enum PIDSpecies {
  kEl = 0,
  kKa,
  kPi,
  kPr
};

enum V0Species {
  kPhoton = 0,
  kK0S,
  kLambda,
  kAntiLambda,
  kNV0
};

static const std::array<std::string, kNCharmParticles> charmParticleNames{"D0", "Dplus", "Ds", "Lc", "Xic"};
static const std::array<std::string, kNBeautyParticles> beautyParticleNames{"Bplus", "B0toDStar", "B0", "Bs", "Lb", "Xib"};
static const std::array<int, kNCharmParticles> pdgCodesCharm{421, 411, 431, 4122, 4232};
static const std::array<std::string, kNtriggersHF + 2> eventTitles = {"all", "rejected", "w/ high-#it{p}_{T} 2p charm", "w/ high-#it{p}_{T} 3p charm", "w/ 3p beauty", "w/ 4p beauty", "w/ 2p femto", "w/ 3p femto", "w/ 2p double charm", "w/ 3p double charm", "w/ 2p and 3p double charm", "w/ 2p + V0", "w/ 3p + V0", "w/ charm baryon"};
static const std::array<std::string, kNtriggersHF> hfTriggerNames{"kHighPt2P", "kHighPt3P", "kBeauty3P", "kBeauty4P", "kFemto2P", "kFemto3P", "kDoubleCharm2P", "kDoubleCharm3P", "kDoubleCharmMix", "kV0Charm2P", "kV0Charm3P", "kCharmBarToXiBach"};
static const std::array<std::string, kNV0> v0Labels{"#gamma", "K_{S}^{0}", "#Lambda", "#bar{#Lambda}"};
static const std::array<std::string, kNV0> v0Names{"Photon", "K0S", "Lambda", "AntiLambda"};

static const std::tuple pdgCharmDaughters{
  std::array{-321, 211},        // D0
  std::array{-321, 211, 211},   // Dplus
  std::array{321, -321, 211},   // Ds
  std::array{2212, -321, 211},  // Lc
  std::array{2212, -321, 211}}; // Xic

static const float massPi = 0.13957;
static const float massK = 0.493677;
static const float massProton = 0.938272;
static const float massPhi = 1.019455;
static const float massD0 = 1.86484;
static const float massDPlus = 1.86962;
static const float massDs = 1.9685;
static const float massLc = 2.28646;
static const float massXic = 2.4679;
static const float massDStar = 2.01027;
static const float massBPlus = 5.27915;
static const float massB0 = 5.27953;
static const float massBs = 5.3663;
static const float massLb = 5.6202;
static const float massXib = 5.7924;
static const float massGamma = 0.;
static const float massK0S = 0.497614;
static const float massLambda = 1.11568;
static const float massXi = 1.32171;

static const o2::framework::AxisSpec ptAxis{50, 0.f, 50.f};
static const o2::framework::AxisSpec pAxis{50, 0.f, 10.f};
static const o2::framework::AxisSpec kstarAxis{100, 0.f, 1.f};
static const o2::framework::AxisSpec etaAxis{30, -1.5f, 1.5f};
static const o2::framework::AxisSpec nSigmaAxis{100, -10.f, 10.f};
static const o2::framework::AxisSpec alphaAxis{100, -1.f, 1.f};
static const o2::framework::AxisSpec qtAxis{100, 0.f, 0.25f};
static const o2::framework::AxisSpec bdtAxis{100, 0.f, 1.f};
static const o2::framework::AxisSpec phiAxis{36, 0., o2::constants::math::TwoPI};
static const std::array<o2::framework::AxisSpec, kNCharmParticles + 8> massAxisC = {o2::framework::AxisSpec{100, 1.65f, 2.05f}, o2::framework::AxisSpec{100, 1.65f, 2.05f}, o2::framework::AxisSpec{100, 1.75f, 2.15f}, o2::framework::AxisSpec{100, 2.05f, 2.45f}, o2::framework::AxisSpec{100, 2.25f, 2.65f}, o2::framework::AxisSpec{100, 0.139f, 0.159f}, o2::framework::AxisSpec{100, 0.f, 0.25f}, o2::framework::AxisSpec{100, 0.f, 0.25f}, o2::framework::AxisSpec{100, 0.48f, 0.88f}, o2::framework::AxisSpec{100, 0.48f, 0.88f}, o2::framework::AxisSpec{100, 1.1f, 1.4f}, o2::framework::AxisSpec{100, 2.3f, 2.9f}, o2::framework::AxisSpec{100, 2.3f, 2.9f}};
static const std::array<o2::framework::AxisSpec, kNBeautyParticles> massAxisB = {o2::framework::AxisSpec{240, 4.8f, 6.0f}, o2::framework::AxisSpec{240, 4.8f, 6.0f}, o2::framework::AxisSpec{240, 4.8f, 6.0f}, o2::framework::AxisSpec{240, 4.8f, 6.0f}, o2::framework::AxisSpec{240, 5.0f, 6.2f}, o2::framework::AxisSpec{240, 5.0f, 6.2f}};

// default values for configurables
// channels to trigger on for femto
constexpr int activeFemtoChannels[1][5] = {{1, 1, 1, 1, 0}}; // pD0, pD+, pDs, pLc, pXic
static const std::vector<std::string> labelsColumnsFemtoChannels = {"protonDZero", "protonDPlus", "protonDs", "protonLc", "protonXic"};

// min pT for all tracks combined  (except for V0 and cascades)
constexpr float cutsMinPt[1][4] = {{0.5, 0.1, 0.8, 0.5}}; // beauty, D*, femto, charm baryons
static const std::vector<std::string> labelsColumnsMinPt = {"Beauty", "DstarPlus", "Femto", "CharmBaryon"};

// min pT for all tracks combined  (except for V0 and cascades)
constexpr float cutsNsigma[3][5] = {{3., 3., 3., 5., 3.},           // TPC proton from Lc, pi/K from D0, K from 3-prong, femto, pi/K from Xic/Omegac
                                    {3., 3., 3., 2.5, 3.},          // TOF proton from Lc, pi/K from D0, K from 3-prong, femto, pi/K from Xic/Omegac
                                    {999., 999., 999., 2.5, 999.}}; // Sum in quadrature of TPC and TOF (used only for femto for pT < 4 GeV/c)
static const std::vector<std::string> labelsColumnsNsigma = {"PrFromLc", "PiKaFromDZero", "KaFrom3Prong", "Femto", "PiKaFromCharmBaryon"};
static const std::vector<std::string> labelsRowsNsigma = {"TPC", "TOF", "Comb"};

// high pt
constexpr float cutsHighPtThresholds[1][2] = {{8., 8.}}; // 2-prongs, 3-prongs
static const std::vector<std::string> labelsColumnsHighPtThresholds = {"2Prongs", "3Prongs"};

// beauty
constexpr float cutsDeltaMassB[1][kNBeautyParticles + 1] = {{0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.04}}; // B+, B0, B0toDstar, Bs, Lb, Xib, charm daughter
static const std::vector<std::string> labelsColumnsDeltaMassB = {"Bplus", "BZero", "BZeroToDstar", "Bs", "Lb", "Xib", "CharmDau"};

// double charm
constexpr int activeDoubleCharmChannels[1][3] = {{1, 1, 1}}; // kDoubleCharm2P, kDoubleCharm3P, kDoubleCharmMix
static const std::vector<std::string> labelsColumnsDoubleCharmChannels = {"DoubleCharm2Prong", "DoubleCharm3Prong", "DoubleCharmMix"};

// charm resonances
constexpr float cutsMassCharmReso[1][6] = {{0.01, 0.3, 0.3, 0.88, 0.88, 1.4}}; // D*+, D*0, Ds*0, Ds1+, Ds2*+, Xic*
static const std::vector<std::string> labelsColumnsDeltaMasseCharmReso = {"DstarPlus", "DstarZero", "DsStarZero", "Ds1Plus", "Ds2StarPlus", "XicStar"};
// V0s for charm resonances
constexpr float cutsV0s[1][6] = {{0.85, 0.97, 0.5, 4., 0.02, 0.01}}; // cosPaGamma, cosPaK0sLambda, radiusK0sLambda, nSigmaPrLambda, deltaMassK0S, deltaMassLambda
static const std::vector<std::string> labelsColumnsV0s = {"CosPaGamma", "CosPaK0sLambda", "RadiusK0sLambda", "NSigmaPrLambda", "DeltaMassK0s", "DeltaMassLambda"};

// cascades for Xi + bachelor triggers
constexpr float cutsCascades[1][7] = {{0.2, 0.01, 0.01, 0.99, 0.99, 0.3, 3.}}; // ptXiBachelor, deltaMassXi, deltaMassLambda, cosPaXi, cosPaLambda, DCAxyXi, nSigmaPid
static const std::vector<std::string> labelsColumnsCascades = {"PtBachelor", "DeltaMassXi", "DeltaMassLambda", "CosPAXi", "CosPaLambda", "DCAxyXi", "NsigmaPid"};
constexpr float cutsCharmBaryons[1][4] = {{3., 3., 2.35, 2.60}}; // MinPtXiPi, MinPtXiKa, MinMassXiPi, MinMassXiKa
static const std::vector<std::string> labelsColumnsCharmBaryons = {"MinPtXiPi", "MinPtXiKa", "MinMassXiPi", "MinMassXiKa"};

// dummy array
static const std::vector<std::string> labelsEmpty{};
static constexpr double cutsTrackDummy[hf_cuts_single_track::nBinsPtTrack][hf_cuts_single_track::nCutVarsTrack] = {{0., 10.}, {0., 10.}, {0., 10.}, {0., 10.}, {0., 10.}, {0., 10.}};
o2::framework::LabeledArray<double> cutsSingleTrackDummy{cutsTrackDummy[0], hf_cuts_single_track::nBinsPtTrack, hf_cuts_single_track::nCutVarsTrack, hf_cuts_single_track::labelsPtTrack, hf_cuts_single_track::labelsCutVarTrack};

// Main helper class

class HfFilterHelper
{
 public:
  /// Default constructor
  HfFilterHelper() = default;

  // selections
  template <typename T, typename T1, typename T2, typename T3, typename T4>
  int8_t isSelectedTrackForSoftPionOrBeauty(const T track, const T1& trackPar, const T2& dca, const float& pTMinSoftPion, const float& pTMinBeautyBachelor, const T3& pTBinsTrack, const T4& cutsSingleTrackBeauty)
  template <typename T1, typename T2, typename H2, typename H3>
  bool isSelectedProton4Femto(const T1& track, const T2& trackPar, const float& femtoMinProtonPt, const float& femtoMaxNsigmaProton, const int setTPCCalib, H3 hMapProton, const std::array<std::vector<double>, 2>& hSplineProton, const int& activateQA, H2 hProtonTPCPID, H2 hProtonTOFPID)
  template <typename T, typename H3>
  bool isSelectedProton4CharmBaryons(const T& track, const float& nsigmaTPCProtonLc, const float& nsigmaTOFProtonLc, const int setTPCCalib, H3 hMapProton, const std::array<std::vector<double>, 2>& hSplineProton);
  template <typename T, typename H3>
  int8_t isDzeroPreselected(const T& trackPos, const T& trackNeg, const float& nsigmaTPCPionKaonDzero, const float& nsigmaTOFPionKaonDzero, const int& setTPCCalib, H3 hMapPion, const std::array<std::vector<double>, 2>& hSplinePion, const std::array<std::vector<double>, 2>& hSplineKaon)
  template <typename T, typename H3>
  int8_t isDplusPreselected(const T& trackOppositeCharge, const float& nsigmaTPCKaon3Prong, const float& nsigmaTOFKaon3Prong, const int setTPCCalib, H3 hMapKaon, std::array<std::vector<double>, 2>& hSplineKaon)
  template <typename P, typename T, typename H3>
  int8_t isDsPreselected(const P& pTrackSameChargeFirst, const P& pTrackSameChargeSecond, const P& pTrackOppositeCharge, const T& trackOppositeCharge, const float& nsigmaTPCKaon3Prong, const float& nsigmaTOFKaon3Prong, const int setTPCCalib, H3 hMapKaon, const std::array<std::vector<double>, 2>& hSplineKaon)
  template <typename T, typename H3>
  int8_t isCharmBaryonPreselected(const T& trackSameChargeFirst, const T& trackSameChargeSecond, const T& trackOppositeCharge, const float& nsigmaTPCProtonLc, const float& nsigmaTOFProtonLc, const float& nsigmaTPCKaon3Prong, const float& nsigmaTOFKaon3Prong, const int setTPCCalib, H3 hMapProton, const std::array<std::vector<double>, 2>& hSplineProton, H3 hMapKaon, const std::array<std::vector<double>, 2>& hSplineKaon)
  template <typename T, typename H2>
  int8_t isSelectedD0InMassRange(const T& pTrackPos, const T& pTrackNeg, const float& ptD, const float& phiD, int8_t isSelected, const float& deltaMassCharmHadronForBeauty, const int& activateQA, H2 hMassVsPt)
  template <typename T, typename H2>
  int8_t isSelectedDplusInMassRange(const T& pTrackSameChargeFirst, const T& pTrackSameChargeSecond, const T& pTrackOppositeCharge, const float& ptD, const float& phiD, const float& deltaMassCharmHadronForBeauty, const int& activateQA, H2 hMassVsPt)
  template <typename T, typename H2>
  int8_t isSelectedDsInMassRange(const T& pTrackSameChargeFirst, const T& pTrackSameChargeSecond, const T& pTrackOppositeCharge, const float& ptD, const float& phiD, int8_t isSelected, const float& deltaMassCharmHadronForBeauty, const int& activateQA, H2 hMassVsPt)
  template <typename T, typename H2>
  int8_t isSelectedLcInMassRange(const T& pTrackSameChargeFirst, const T& pTrackSameChargeSecond, const T& pTrackOppositeCharge, const float& ptLc, const float& phiLc, const int8_t isSelected, const float& deltaMassCharmHadronForBeauty, const int& activateQA, H2 hMassVsPt)
  template <typename T, typename H2>
  int8_t isSelectedXicInMassRange(const T& pTrackSameChargeFirst, const T& pTrackSameChargeSecond, const T& pTrackOppositeCharge, const float& ptXic, const float& phiXic, const int8_t isSelected, const float& deltaMassCharmHadronForBeauty, const int& activateQA, H2 hMassVsPt)
  template <typename V0, typename Coll, typename T, typename H2, typename H3>
  int8_t isSelectedV0(const V0& v0, const array<T, 2>& dauTracks, const Coll& collision, const float& minGammaCosinePa, const float& minV0CosinePa, const float& minV0Radius, const float& maxNsigmaPrForLambda, const float& deltaMassK0s, const float& deltaMassLambda, const int& setTPCCalib, H3 hMapProton, const std::array<std::vector<double>, 2>& hSplineProton, const int& activateQA, H2 hV0Selected, std::array<H2, 4>& hArmPod)
  template <typename Casc, typename V0, typename T, typename Coll, typename H3>
  bool isSelectedCascade(const Casc& casc, const V0& v0, const array<T, 3>& dauTracks, const Coll& collision, const float& minPtXiBachelor, const float& deltaMassXi, const float& deltaMassLambda, const float& cosPAXi, const float& cosPALambda, const float& DCAxyXi, const float& maxNsigma, const int& setTPCCalib, H3 hMapPion, H3 hMapProton, const std::array<std::vector<double>, 2>& hSplinePion, const std::array<std::vector<double>, 2>& hSplineProton)
  template <typename T, typename T2, typename T3, typename T4, typename H3>
  int8_t isSelectedBachelorForCharmBaryon(const T& track, const T2& dca, const float& minPt, const T3& pTBinsTrack, const T4& cutsSingleTrack, const float& maxNsigmaTPC, const float& maxNsigmaTOF, const int& setTPCCalib, H3 hMapPion, const std::array<std::vector<double>, 2>& hSplinePion, const std::array<std::vector<double>, 2>& hSplineKaon)
  template <typename T, typename U>
  int8_t isBDTSelected(const T& scores, const U& thresholdBDTScores)

  // helpers
  template <typename T>
  T computeRelativeMomentum(const std::array<T, 3>& pTrack, const std::array<T, 3>& CharmCandMomentum, const T& CharmMass)
  template <typename T>
  int computeNumberOfCandidates(std::vector<std::vector<T>> indices)

  // PID
  std::vector<double> setValuesBB(o2::ccdb::CcdbApi& ccdbApi, aod::BCsWithTimestamps::iterator const& bunchCrossing, const std::string ccdbPath);

  // ML
  Ort::Experimental::Session* initONNXSession(std::string& onnxFile, std::string partName, Ort::Env& env, Ort::SessionOptions& sessionOpt, std::vector<std::vector<int64_t>>& inputShapes, int& dataType, bool loadModelsFromCCDB, o2::ccdb::CcdbApi& ccdbApi, std::string mlModelPathCCDB, int64_t timestampCCDB)
  template <typename T>
  std::array<T, 3> predictONNX(std::vector<T>& inputFeatures, std::shared_ptr<Ort::Experimental::Session>& session, std::vector<std::vector<int64_t>>& inputShapes)

 private:
  // selections
  template <typename T, typename H3>
  bool isSelectedKaon4Charm3Prong(const T& track, const float& nsigmaTPCKaon3Prong, const float& nsigmaTOFKaon3Prong, const int setTPCCalib, H3 hMapKaon, const std::array<std::vector<double>, 2>& hSplineKaon)

  // PID
  template <typename T>
  double getTPCSplineCalib(const T& track, const float mMassPar, const std::vector<double> paraBetheBloch);
  template <typename T, typename H3>
  float getTPCPostCalib(const array<H3, 2>& hCalibMap, const T& track, const int pidSpecies)

  std::array<TH3F*, 2> mHistMapPion{nullptr, nullptr};    // Map for TPC PID postcalibrations for pions
  std::array<TH3F*, 2> mHistMapProton{nullptr, nullptr};  // Map for TPC PID postcalibrations for protons
  std::array<std::vector<double>, 2> mBetheBlochPion{};   // Bethe-Bloch parametrisations for pions and antipions in TPC
  std::array<std::vector<double>, 2> mBetheBlochKaon{};   // Bethe-Bloch parametrisations for kaons and antikaons in TPC
  std::array<std::vector<double>, 2> mBetheBlochProton{}; // Bethe-Bloch parametrisations for proton and antiprotons in TPC
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
                  hfoptimisationTree::NsigmaPrTOF);
DECLARE_SOA_TABLE(HFOptimisationTreeCollisions, "AOD", "HFOPTIMTREECOLL", //!
                  hfoptimisationTree::CollisionIndex)
} // namespace o2::aod

#endif // EVENTFILTERING_PWGHF_HFFILTERHELPERS_H_
