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
#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/MathConstants.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/DataTypes.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/AnalysisHelpers.h"
#include "Framework/O2DatabasePDGPlugin.h"

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
  kPi,
  kAntiPi,
  kKa,
  kAntiKa,
  kPr,
  kAntiPr
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

constexpr float massPi = o2::analysis::pdg::MassPiPlus;
constexpr float massKa = o2::analysis::pdg::MassKPlus;
constexpr float massProton = o2::analysis::pdg::MassProton;
constexpr float massGamma = o2::analysis::pdg::MassGamma;
constexpr float massK0S = o2::analysis::pdg::MassK0Short;
constexpr float massLambda = o2::analysis::pdg::MassLambda0;
constexpr float massXi = o2::analysis::pdg::MassXiMinus;
constexpr float massPhi = o2::analysis::pdg::MassPhi;
constexpr float massD0 = o2::analysis::pdg::MassD0;
constexpr float massDPlus = o2::analysis::pdg::MassDPlus;
constexpr float massDs = o2::analysis::pdg::MassDS;
constexpr float massLc = o2::analysis::pdg::MassLambdaCPlus;
constexpr float massXic = o2::analysis::pdg::MassXiCPlus;
constexpr float massDStar = o2::analysis::pdg::MassDStar;
constexpr float massBPlus = o2::analysis::pdg::MassBPlus;
constexpr float massB0 = o2::analysis::pdg::MassB0;
constexpr float massBs = o2::analysis::pdg::MassBS;
constexpr float massLb = o2::analysis::pdg::MassLambdaB0;
constexpr float massXib = o2::analysis::pdg::MassXiB0;

static const o2::framework::AxisSpec ptAxis{50, 0.f, 50.f};
static const o2::framework::AxisSpec pAxis{50, 0.f, 10.f};
static const o2::framework::AxisSpec kstarAxis{100, 0.f, 1.f};
static const o2::framework::AxisSpec etaAxis{30, -1.5f, 1.5f};
static const o2::framework::AxisSpec nSigmaAxis{100, -10.f, 10.f};
static const o2::framework::AxisSpec alphaAxis{100, -1.f, 1.f};
static const o2::framework::AxisSpec qtAxis{100, 0.f, 0.25f};
static const o2::framework::AxisSpec bdtAxis{100, 0.f, 1.f};
static const o2::framework::AxisSpec phiAxis{36, 0., o2::constants::math::TwoPI};
static const std::array<o2::framework::AxisSpec, kNCharmParticles + 8> massAxisC = {o2::framework::AxisSpec{100, 1.65f, 2.05f}, o2::framework::AxisSpec{100, 1.65f, 2.05f}, o2::framework::AxisSpec{100, 1.75f, 2.15f}, o2::framework::AxisSpec{100, 2.05f, 2.45f}, o2::framework::AxisSpec{100, 2.25f, 2.65f}, o2::framework::AxisSpec{100, 0.139f, 0.159f}, o2::framework::AxisSpec{100, 0.f, 0.25f}, o2::framework::AxisSpec{100, 0.f, 0.25f}, o2::framework::AxisSpec{200, 0.48f, 0.88f}, o2::framework::AxisSpec{200, 0.48f, 0.88f}, o2::framework::AxisSpec{100, 1.1f, 1.4f}, o2::framework::AxisSpec{100, 2.3f, 2.9f}, o2::framework::AxisSpec{100, 2.3f, 2.9f}};
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
static constexpr double cutsTrackDummy[o2::analysis::hf_cuts_single_track::nBinsPtTrack][o2::analysis::hf_cuts_single_track::nCutVarsTrack] = {{0., 10.}, {0., 10.}, {0., 10.}, {0., 10.}, {0., 10.}, {0., 10.}};
o2::framework::LabeledArray<double> cutsSingleTrackDummy{cutsTrackDummy[0], o2::analysis::hf_cuts_single_track::nBinsPtTrack, o2::analysis::hf_cuts_single_track::nCutVarsTrack, o2::analysis::hf_cuts_single_track::labelsPtTrack, o2::analysis::hf_cuts_single_track::labelsCutVarTrack};

// Main helper class

class HfFilterHelper
{
 public:
  /// Default constructor
  HfFilterHelper() = default;

  // setters
  void setPtBinsSingleTracks(std::vector<double> ptBins) { mPtBinsTracks = ptBins; }
  void setMinPtBeautyBachelor(float minPt) { mPtMinBeautyBachelor = minPt; }
  void setMinPtDstarSoftPion(float minPt) { mPtMinSoftPionForDstar = minPt; }
  void setCutsSingleTrackBeauty(o2::framework::LabeledArray<double> cutsSingleTrack3P, o2::framework::LabeledArray<double> cutsSingleTrack4P)
  {
    mCutsSingleTrackBeauty3Prong = cutsSingleTrack3P;
    mCutsSingleTrackBeauty4Prong = cutsSingleTrack4P;
  }
  void setMinPtProtonForFemto(float minPt) { mPtMinProtonForFemto = minPt; }
  void setPtThresholdPidStrategyForFemto(float ptThreshold) { mPtThresholdPidStrategyForFemto = ptThreshold; }
  void setNsigmaProtonCutsForFemto(std::array<float, 3> nSigmaCuts) { mNSigmaPrCutsForFemto = nSigmaCuts; }
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
  void setDeltaMassCharmHadForBeauty(float delta) { mDeltaMassCharmHadForBeauty = delta; }
  void setV0Selections(float minGammaCosPa, float minK0sLambdaCosPa, float minK0sLambdaRadius, float nSigmaPrFromLambda, float deltaMassK0s, float deltaMassLambda)
  {
    mMinGammaCosinePa = minGammaCosPa;
    mMinK0sLambdaCosinePa = minK0sLambdaCosPa;
    mMinK0sLambdaRadius = minK0sLambdaRadius;
    mMaxNsigmaPrForLambda = nSigmaPrFromLambda;
    mDeltaMassK0s = deltaMassK0s;
    mDeltaMassLambda = deltaMassLambda;
  }
  void setXiSelections(float minPtXiBachelor, float deltaMassXi, float deltaMassLambda, float cosPaXi, float cosPaLambdaFromXi, float maxDcaxyXi, float nSigma)
  {
    mMinPtXiBachelor = minPtXiBachelor;
    mDeltaMassXi = deltaMassXi;
    mDeltaMassLambdaFromXi = deltaMassLambda;
    mCosPaXi = cosPaXi;
    mCosPaLambdaFromXi = cosPaLambdaFromXi;
    mMaxDcaXyXi = maxDcaxyXi;
    mMaxNsigmaXiDau = nSigma;
  }
  void setCutsSingleTrackCharmBaryonBachelor(o2::framework::LabeledArray<double> cutsSingleTrack) { mCutsSingleTrackCharmBaryonBachelor = cutsSingleTrack; }
  void setMinPtCharmBaryonBachelor(float minPt) { mPtMinCharmBaryonBachelor = minPt; }
  void setNsigmaPiCutsForCharmBaryonBachelor(float nSigmaTpc, float nSigmaTof)
  {
    mNSigmaTpcPiCharmBaryonBachelor = nSigmaTpc;
    mNSigmaTofPiCharmBaryonBachelor = nSigmaTof;
  }

  void setTpcPidCalibrationOption(int opt) { mTpcPidCalibrationOption = opt; }

  // helper functions for selections
  template <typename T, typename T1, typename T2>
  int8_t isSelectedTrackForSoftPionOrBeauty(const T track, const T1& trackPar, const T2& dca, const int& whichTrigger);
  template <typename T1, typename T2, typename H2>
  bool isSelectedProton4Femto(const T1& track, const T2& trackPar, const int& activateQA, H2 hProtonTPCPID, H2 hProtonTOFPID);
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
  template <typename T, typename H2>
  int8_t isSelectedXicInMassRange(const T& pTrackSameChargeFirst, const T& pTrackSameChargeSecond, const T& pTrackOppositeCharge, const float& ptXic, const int8_t isSelected, const int& activateQA, H2 hMassVsPt);
  template <typename V0, typename Coll, typename T, typename H2>
  int8_t isSelectedV0(const V0& v0, const std::array<T, 2>& dauTracks, const Coll& collision, const int& activateQA, H2 hV0Selected, std::array<H2, 4>& hArmPod);
  template <typename Casc, typename V0, typename T, typename Coll>
  bool isSelectedCascade(const Casc& casc, const V0& v0, const std::array<T, 3>& dauTracks, const Coll& collision);
  template <typename T, typename T2>
  int8_t isSelectedBachelorForCharmBaryon(const T& track, const T2& dca);
  template <typename T, typename U>
  int8_t isBDTSelected(const T& scores, const U& thresholdBDTScores);

  // helpers
  template <typename T>
  T computeRelativeMomentum(const std::array<T, 3>& pTrack, const std::array<T, 3>& CharmCandMomentum, const T& CharmMass);
  template <typename T>
  int computeNumberOfCandidates(std::vector<std::vector<T>> indices);

  // PID
  void setValuesBB(o2::ccdb::CcdbApi& ccdbApi, aod::BCsWithTimestamps::iterator const& bunchCrossing, const std::array<std::string, 6>& ccdbPaths);
  void setTpcRecalibMaps(o2::framework::Service<o2::ccdb::BasicCCDBManager> const& ccdb, aod::BCsWithTimestamps::iterator const& bunchCrossing, const std::string& ccdbPath);

  // ML
  Ort::Experimental::Session* initONNXSession(std::string& onnxFile, std::string partName, Ort::Env& env, Ort::SessionOptions& sessionOpt, std::vector<std::vector<int64_t>>& inputShapes, int& dataType, bool loadModelsFromCCDB, o2::ccdb::CcdbApi& ccdbApi, std::string mlModelPathCCDB, int64_t timestampCCDB);
  template <typename T>
  std::array<T, 3> predictONNX(std::vector<T>& inputFeatures, std::shared_ptr<Ort::Experimental::Session>& session, std::vector<std::vector<int64_t>>& inputShapes);

 private:
  // selections
  template <typename T>
  bool isSelectedKaon4Charm3Prong(const T& track, const float& nsigmaTPCKaon3Prong, const float& nsigmaTOFKaon3Prong);
  template <typename T>
  bool isSelectedProton4CharmBaryons(const T& track, const float& nsigmaTPCProton, const float& nsigmaTOFProton);

  // PID
  template <typename T>
  double getTPCSplineCalib(const T& track, const int& pidSpecies);
  template <typename T>
  float getTPCPostCalib(const T& track, const int& pidSpecies);

  // helpers
  template <typename T1, typename T2>
  int findBin(T1 const& binsPt, T2 value);

  // selections
  std::vector<double> mPtBinsTracks{};                                       // vector of pT bins for single track cuts
  o2::framework::LabeledArray<double> mCutsSingleTrackBeauty3Prong{};        // dca selections for the 3-prong b-hadron pion daughter
  o2::framework::LabeledArray<double> mCutsSingleTrackBeauty4Prong{};        // dca selections for the 4-prong b-hadron pion daughter
  float mPtMinSoftPionForDstar{0.1};                                         // minimum pt for the D*+ soft pion
  float mPtMinBeautyBachelor{0.5};                                           // minimum pt for the b-hadron pion daughter
  float mPtMinProtonForFemto{0.8};                                           // minimum pt for the proton for femto
  float mPtThresholdPidStrategyForFemto{8.};                                 // pt threshold to change strategy for proton PID for femto
  std::array<float, 3> mNSigmaPrCutsForFemto{3., 3., 3.};                    // cut values for Nsigma TPC, TOF, combined for femto protons
  float mNSigmaTpcPrCutForCharmBaryons{3.};                                  // maximum Nsigma TPC for protons in Lc and Xic decays
  float mNSigmaTofPrCutForCharmBaryons{3.};                                  // maximum Nsigma TOF for protons in Lc and Xic decays
  float mNSigmaTpcKaCutFor3Prongs{3.};                                       // maximum Nsigma TPC for kaons in 3-prong decays
  float mNSigmaTofKaCutFor3Prongs{3.};                                       // maximum Nsigma TOF for kaons in 3-prong decays
  float mNSigmaTpcPiKaCutForDzero{3.};                                       // maximum Nsigma TPC for pions/kaons in D0 decays
  float mNSigmaTofPiKaCutForDzero{3.};                                       // maximum Nsigma TOF for pions/kaons in D0 decays
  float mDeltaMassCharmHadForBeauty{0.04};                                   // delta mass cut for charm hadrons in B and charm excited decays
  float mMinGammaCosinePa{0.85};                                             // minimum cosp for gammas
  float mMinK0sLambdaCosinePa{0.97};                                         // minimum cosp for K0S and Lambda in charm excited decays
  float mMinK0sLambdaRadius{0.5};                                            // minimum radius for K0S and Lambda in charm excited decays
  float mMaxNsigmaPrForLambda{4.};                                           // maximum Nsigma TPC and TOF for protons in Lambda decays
  float mDeltaMassK0s{0.02};                                                 // delta mass cut for K0S in charm excited decays
  float mDeltaMassLambda{0.01};                                              // delta mass cut for Lambda in charm excited decays
  float mMinPtXiBachelor{0.1};                                               // minimum pt for Xi bachelor
  float mDeltaMassXi{0.01};                                                  // delta mass cut for Xi in Xic/Omegac decays
  float mDeltaMassLambdaFromXi{0.01};                                        // delta mass cut for Lambda <- Xi in Xic/Omegac decays
  float mCosPaXi{0.99};                                                      // minimum cosp for Xi in Xic/Omegac decays
  float mCosPaLambdaFromXi{0.99};                                            // minimum cosp for Xi in Xic/Omegac decays
  float mMaxDcaXyXi{0.3};                                                    // maximum dca for Xi in Xic/Omegac decays
  float mMaxNsigmaXiDau{3.};                                                 // maximum Nsigma TPC and TOF for Xi daughter tracks
  float mPtMinCharmBaryonBachelor{0.5};                                      // minimum pt for the bachelor pion from Xic/Omegac decays
  o2::framework::LabeledArray<double> mCutsSingleTrackCharmBaryonBachelor{}; // dca selections for the bachelor pion from Xic/Omegac decays
  float mNSigmaTpcPiCharmBaryonBachelor{3.};                                 // maximum Nsigma TPC for pions in Xic/Omegac decays
  float mNSigmaTofPiCharmBaryonBachelor{3.};                                 // maximum Nsigma TOF for pions in Xic/Omegac decays

  // PID recalibrations
  int mTpcPidCalibrationOption{0};                        // Option for TPC PID calibration (0 -> AO2D, 1 -> postcalibrations, 2 -> alternative bethe bloch parametrisation)
  std::array<TH3F*, 4> mHistMapPiPr{};                    // Map for TPC PID postcalibrations for pions and protons
  std::array<std::vector<double>, 6> mBetheBlochPiKaPr{}; // Bethe-Bloch parametrisations for pions, antipions, kaons, antikaons, protons, antiprotons in TPC
};

/// Single-track cuts for bachelor track of beauty candidates
/// \param track is a track parameter
/// \param trackPar is a track parameter
/// \param dca is the 2d array with dcaXY and dcaZ of the track
/// \return a flag that encodes the selection for soft pions BIT(kSoftPion), tracks for beauty BIT(kForBeauty), or soft pions for beauty BIT(kSoftPionForBeauty)
template <typename T, typename T1, typename T2>
inline int8_t HfFilterHelper::isSelectedTrackForSoftPionOrBeauty(const T track, const T1& trackPar, const T2& dca, const int& whichTrigger)
{

  int8_t retValue{BIT(kSoftPion) | BIT(kForBeauty) | BIT(kSoftPionForBeauty)};

  if (!track.isGlobalTrackWoDCA()) {
    return kRejected;
  }

  auto pT = trackPar.getPt();
  auto pTBinTrack = findBin(mPtBinsTracks, pT);
  if (pTBinTrack == -1) {
    return kRejected;
  }

  if (pT < mPtMinSoftPionForDstar) { // soft pion should be less stringent than usual tracks
    return kRejected;
  }

  if (std::fabs(trackPar.getEta()) > 0.8) {
    return kRejected;
  }

  if (std::fabs(dca[1]) > 2.f) {
    return kRejected;
  }

  // below only regular beauty tracks, not required for soft pions
  if (pT < mPtMinBeautyBachelor) {
    CLRBIT(retValue, kForBeauty);
  }

  float minDca = 1000.f;
  float maxDca = 0.f;
  if (whichTrigger == kBeauty3P) {
    minDca = mCutsSingleTrackBeauty3Prong.get(pTBinTrack, 0u);
    maxDca = mCutsSingleTrackBeauty3Prong.get(pTBinTrack, 1u);
  } else if (whichTrigger == kBeauty4P) {
    minDca = mCutsSingleTrackBeauty4Prong.get(pTBinTrack, 0u);
    maxDca = mCutsSingleTrackBeauty4Prong.get(pTBinTrack, 1u);
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

/// Basic selection of proton candidates
/// \param track is a track
/// \param trackPar is a track parameter
/// \param activateQA flag to activate the filling of QA histos
/// \param hProtonTPCPID histo with NsigmaTPC vs. p
/// \param hProtonTOFPID histo with NsigmaTOF vs. p
/// \return true if track passes all cuts
template <typename T1, typename T2, typename H2>
inline bool HfFilterHelper::isSelectedProton4Femto(const T1& track, const T2& trackPar, const int& activateQA, H2 hProtonTPCPID, H2 hProtonTOFPID)
{
  if (trackPar.getPt() < mPtMinProtonForFemto) {
    return false;
  }

  if (std::fabs(trackPar.getEta()) > 0.8) {
    return false;
  }

  if (!track.isGlobalTrack()) {
    return false; // use only global tracks
  }

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

  float NSigma = std::sqrt(NSigmaTPC * NSigmaTPC + NSigmaTOF * NSigmaTOF);

  if (trackPar.getPt() <= mPtThresholdPidStrategyForFemto) {
    if (NSigma > mNSigmaPrCutsForFemto[2]) {
      return false;
    }
  } else {
    if (std::fabs(NSigmaTPC) > mNSigmaPrCutsForFemto[0] || std::fabs(NSigmaTOF) > mNSigmaPrCutsForFemto[1]) {
      return false;
    }
  }

  if (activateQA > 1) {
    hProtonTPCPID->Fill(track.p(), NSigmaTPC);
    hProtonTOFPID->Fill(track.p(), NSigmaTOF);
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
  if (!isSelectedKaon4Charm3Prong(trackOppositeCharge, mNSigmaTpcKaCutFor3Prongs, mNSigmaTofKaCutFor3Prongs)) {
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
  if (!isSelectedKaon4Charm3Prong(trackOppositeCharge, mNSigmaTpcKaCutFor3Prongs, mNSigmaTofKaCutFor3Prongs)) {
    return retValue;
  }

  // check delta-mass for phi resonance
  auto invMassKKFirst = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackOppositeCharge}, std::array{massKa, massKa});
  auto invMassKKSecond = RecoDecay::m(std::array{pTrackSameChargeSecond, pTrackOppositeCharge}, std::array{massKa, massKa});

  if (std::fabs(invMassKKFirst - massPhi) < 0.02) {
    retValue |= BIT(0);
  }
  if (std::fabs(invMassKKSecond - massPhi) < 0.02) {
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
  if (!isSelectedKaon4Charm3Prong(trackOppositeCharge, mNSigmaTpcKaCutFor3Prongs, mNSigmaTofKaCutFor3Prongs)) {
    return retValue;
  }
  if (isSelectedProton4CharmBaryons(trackSameChargeFirst, mNSigmaTpcPrCutForCharmBaryons, mNSigmaTofPrCutForCharmBaryons)) {
    retValue |= BIT(0);
  }
  if (isSelectedProton4CharmBaryons(trackSameChargeSecond, mNSigmaTpcPrCutForCharmBaryons, mNSigmaTofPrCutForCharmBaryons)) {
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
  int8_t retValue = 0;
  if (TESTBIT(isSelected, 0)) {
    auto invMassD0 = RecoDecay::m(std::array{pTrackPos, pTrackNeg}, std::array{massPi, massKa});
    if (activateQA) {
      hMassVsPt->Fill(ptD, invMassD0);
    }
    if (std::fabs(invMassD0 - massD0) < mDeltaMassCharmHadForBeauty || ptD > 10) {
      retValue |= BIT(0);
    }
  }
  if (TESTBIT(isSelected, 1)) {
    auto invMassD0bar = RecoDecay::m(std::array{pTrackPos, pTrackNeg}, std::array{massKa, massPi});
    if (activateQA) {
      hMassVsPt->Fill(ptD, invMassD0bar);
    }
    if (std::fabs(invMassD0bar - massD0) < mDeltaMassCharmHadForBeauty || ptD > 10) {
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
  auto invMassDplus = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackSameChargeSecond, pTrackOppositeCharge}, std::array{massPi, massPi, massKa});
  if (activateQA) {
    hMassVsPt->Fill(ptD, invMassDplus);
  }

  if (std::fabs(invMassDplus - massDPlus) > mDeltaMassCharmHadForBeauty && ptD > 0) {
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
  int8_t retValue = 0;
  if (TESTBIT(isSelected, 0)) {
    auto invMassDsToKKPi = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackOppositeCharge, pTrackSameChargeSecond}, std::array{massKa, massKa, massPi});
    if (activateQA) {
      hMassVsPt->Fill(ptD, invMassDsToKKPi);
    }
    if (std::fabs(invMassDsToKKPi - massDs) < mDeltaMassCharmHadForBeauty || ptD > 10) {
      retValue |= BIT(0);
    }
  }
  if (TESTBIT(isSelected, 1)) {
    auto invMassDsToPiKK = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackOppositeCharge, pTrackSameChargeSecond}, std::array{massPi, massKa, massKa});
    if (activateQA) {
      hMassVsPt->Fill(ptD, invMassDsToPiKK);
    }
    if (std::fabs(invMassDsToPiKK - massDs) < mDeltaMassCharmHadForBeauty || ptD > 10) {
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
  int8_t retValue = 0;
  if (TESTBIT(isSelected, 0)) {
    auto invMassLcToPKPi = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackOppositeCharge, pTrackSameChargeSecond}, std::array{massProton, massKa, massPi});
    if (activateQA) {
      hMassVsPt->Fill(ptLc, invMassLcToPKPi);
    }
    if (std::fabs(invMassLcToPKPi - massLc) < mDeltaMassCharmHadForBeauty || ptLc > 10) {
      retValue |= BIT(0);
    }
  }
  if (TESTBIT(isSelected, 1)) {
    auto invMassLcToPiKP = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackOppositeCharge, pTrackSameChargeSecond}, std::array{massPi, massKa, massProton});
    if (activateQA) {
      hMassVsPt->Fill(ptLc, invMassLcToPiKP);
    }
    if (std::fabs(invMassLcToPiKP - massLc) < mDeltaMassCharmHadForBeauty || ptLc > 10) {
      retValue |= BIT(1);
    }
  }

  return retValue;
}

/// Mass selection of Xic candidates to build Lb candidates
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
  int8_t retValue = 0;
  if (TESTBIT(isSelected, 0)) {
    auto invMassXicToPKPi = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackOppositeCharge, pTrackSameChargeSecond}, std::array{massProton, massKa, massPi});
    if (activateQA) {
      hMassVsPt->Fill(ptXic, invMassXicToPKPi);
    }
    if (std::fabs(invMassXicToPKPi - massXic) < mDeltaMassCharmHadForBeauty || ptXic > 10) {
      retValue |= BIT(0);
    }
  }
  if (TESTBIT(isSelected, 1)) {
    auto invMassXicToPiKP = RecoDecay::m(std::array{pTrackSameChargeFirst, pTrackOppositeCharge, pTrackSameChargeSecond}, std::array{massPi, massKa, massProton});
    if (activateQA) {
      hMassVsPt->Fill(ptXic, invMassXicToPiKP);
    }
    if (std::fabs(invMassXicToPiKP - massXic) < mDeltaMassCharmHadForBeauty || ptXic > 10) {
      retValue |= BIT(1);
    }
  }

  return retValue;
}

/// Basic selection of V0 candidates
/// \param v0 is the v0 candidate
/// \param dauTracks is a 2-element array with positive and negative V0 daughter tracks
/// \param collision is the current collision
/// \param activateQA flag to fill QA histos
/// \param hV0Selected is the pointer to the QA histo for selected gammas
/// \param hArmPod is the pointer to an array of QA histo AP plot before selection
/// \return an integer passes all cuts
template <typename V0, typename Coll, typename T, typename H2>
inline int8_t HfFilterHelper::isSelectedV0(const V0& v0, const std::array<T, 2>& dauTracks, const Coll& collision, const int& activateQA, H2 hV0Selected, std::array<H2, 4>& hArmPod)
{
  int8_t isSelected{BIT(kPhoton) | BIT(kK0S) | BIT(kLambda) | BIT(kAntiLambda)};

  if (activateQA > 1) {
    for (int iV0{kPhoton}; iV0 < kNV0; ++iV0) {
      hV0Selected->Fill(0., iV0);
    }
  }

  // eta of daughters
  if (std::fabs(dauTracks[0].eta()) > 1. || std::fabs(dauTracks[1].eta()) > 1.) { // cut all V0 daughters with |eta| > 1.
    if (activateQA > 1) {
      for (int iV0{kPhoton}; iV0 < kNV0; ++iV0) {
        hV0Selected->Fill(1., iV0);
      }
    }
    return kRejected;
  }

  // V0 radius
  if (v0.v0radius() < 0. || v0.v0radius() > 180.) {
    CLRBIT(isSelected, kPhoton);
    if (activateQA > 1) {
      hV0Selected->Fill(2., kPhoton);
    }
  }
  if (v0.v0radius() < mMinK0sLambdaRadius) {
    for (int iV0{kK0S}; iV0 < kNV0; ++iV0) {
      CLRBIT(isSelected, iV0);
      if (activateQA > 1) {
        hV0Selected->Fill(2., iV0);
      }
    }
  }

  auto v0CosinePa = v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ());
  // cosine of pointing angle
  if (TESTBIT(isSelected, kPhoton) && v0CosinePa < mMinGammaCosinePa) {
    CLRBIT(isSelected, kPhoton);
    if (activateQA > 1) {
      hV0Selected->Fill(3., kPhoton);
    }
  }
  for (int iV0{kK0S}; iV0 < kNV0; ++iV0) {
    if (TESTBIT(isSelected, iV0) && v0CosinePa < mMinK0sLambdaCosinePa) {
      CLRBIT(isSelected, iV0);
      if (activateQA > 1) {
        hV0Selected->Fill(3., iV0);
      }
    }
  }

  // armenteros-podolanski / mass
  if (TESTBIT(isSelected, kPhoton) && (std::pow(v0.alpha() / 0.95, 2) + std::pow(v0.qtarm() / 0.05, 2)) >= 1) {
    CLRBIT(isSelected, kPhoton);
    if (activateQA > 1) {
      hV0Selected->Fill(4., kPhoton);
    }
  }
  if (TESTBIT(isSelected, kK0S) && std::fabs(v0.mK0Short() - massK0S) > mDeltaMassK0s) {
    CLRBIT(isSelected, kK0S);
    if (activateQA > 1) {
      hV0Selected->Fill(4., kK0S);
    }
  }
  if (TESTBIT(isSelected, kLambda) && std::fabs(v0.mLambda() - massLambda) > mDeltaMassLambda) {
    CLRBIT(isSelected, kLambda);
    if (activateQA > 1) {
      hV0Selected->Fill(4., kLambda);
    }
  }
  if (TESTBIT(isSelected, kAntiLambda) && std::fabs(v0.mAntiLambda() - massLambda) > mDeltaMassLambda) {
    CLRBIT(isSelected, kAntiLambda);
    if (activateQA > 1) {
      hV0Selected->Fill(4., kAntiLambda);
    }
  }

  // DCA V0 and V0 daughters
  for (int iV0{kK0S}; iV0 < kNV0; ++iV0) {
    if (TESTBIT(isSelected, iV0) && v0.dcav0topv(collision.posX(), collision.posY(), collision.posZ()) > 0.1f) { // we want only primary V0s
      CLRBIT(isSelected, iV0);
      if (activateQA > 1) {
        hV0Selected->Fill(5., iV0);
      }
    }
    if (TESTBIT(isSelected, iV0) && (v0.dcaV0daughters() > 1.f || std::fabs(v0.dcapostopv()) < 0.05f || std::fabs(v0.dcanegtopv()) < 0.05f)) {
      CLRBIT(isSelected, iV0);
      if (activateQA > 1) {
        hV0Selected->Fill(6., iV0);
      }
    }
  }

  // psi pair (photon only)
  if (TESTBIT(isSelected, kPhoton) && std::fabs(v0.psipair()) > 0.1) {
    CLRBIT(isSelected, kPhoton);
    if (activateQA > 1) {
      hV0Selected->Fill(7., kPhoton);
    }
  }

  // PID (Lambda/AntiLambda only)
  float nSigmaPrTpc[2] = {dauTracks[0].tpcNSigmaPr(), dauTracks[1].tpcNSigmaPr()};
  float nSigmaPrTof[2] = {dauTracks[0].tofNSigmaPr(), dauTracks[1].tofNSigmaPr()};
  if (mTpcPidCalibrationOption == 1) {
    for (int iDau{0}; iDau < 2; ++iDau) {
      nSigmaPrTpc[iDau] = getTPCPostCalib(dauTracks[iDau], kPr);
    }
  } else if (mTpcPidCalibrationOption == 2) {
    for (int iDau{0}; iDau < 2; ++iDau) {
      nSigmaPrTpc[iDau] = getTPCSplineCalib(dauTracks[iDau], (iDau == 0) ? kPr : kAntiPr);
    }
  }

  if (TESTBIT(isSelected, kLambda) && ((dauTracks[0].hasTPC() && std::fabs(nSigmaPrTpc[0]) > mMaxNsigmaPrForLambda) || (dauTracks[0].hasTOF() && std::fabs(nSigmaPrTof[0]) > mMaxNsigmaPrForLambda))) {
    CLRBIT(isSelected, kLambda);
    if (activateQA > 1) {
      hV0Selected->Fill(8., kLambda);
    }
  }
  if (TESTBIT(isSelected, kAntiLambda) && ((dauTracks[1].hasTPC() && std::fabs(nSigmaPrTpc[1]) > mMaxNsigmaPrForLambda) || (dauTracks[1].hasTOF() && std::fabs(nSigmaPrTof[1]) > mMaxNsigmaPrForLambda))) {
    CLRBIT(isSelected, kAntiLambda);
    if (activateQA > 1) {
      hV0Selected->Fill(8., kAntiLambda);
    }
  }

  if (activateQA) {
    for (int iV0{kPhoton}; iV0 < kNV0; ++iV0) {
      if (TESTBIT(isSelected, iV0)) {
        hArmPod[iV0]->Fill(v0.alpha(), v0.qtarm());
        if (activateQA > 1) {
          hV0Selected->Fill(9., iV0);
        }
      }
    }
  }

  return isSelected;
}

/// Basic selection of cascade candidates
/// \param casc is the cascade candidate
/// \param v0 is the cascade daughter
/// \param dauTracks is a 3-element array with bachelor, positive and negative V0 daughter tracks
/// \param collision is the collision
/// \return true if cascade passes all cuts
template <typename Casc, typename V0, typename T, typename Coll>
inline bool HfFilterHelper::isSelectedCascade(const Casc& casc, const V0& v0, const std::array<T, 3>& dauTracks, const Coll& collision)
{
  // eta of daughters
  if (std::fabs(dauTracks[0].eta()) > 1. || std::fabs(dauTracks[1].eta()) > 1. || std::fabs(dauTracks[2].eta()) > 1.) { // cut all V0 daughters with |eta| > 1.
    return false;
  }

  // V0 radius
  if (v0.v0radius() < 1.2) {
    return false;
  }

  // cascade radius
  if (casc.cascradius() < 0.6) {
    return false;
  }

  // V0 cosp
  if (casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < mCosPaLambdaFromXi) {
    return false;
  }

  // cascade cosp
  if (casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) < mCosPaXi) {
    return false;
  }

  // cascade DCAxy to PV
  if (std::fabs(casc.dcaXYCascToPV()) > mMaxDcaXyXi) {
    return false;
  }

  // Xi bachelor min pT
  if (dauTracks[0].pt() < mMinPtXiBachelor) {
    return false;
  }

  // dau dca
  if (std::fabs(casc.dcaV0daughters()) > 1.f || std::fabs(casc.dcacascdaughters()) > 1.f) {
    return false;
  }

  // cascade mass
  if (std::fabs(casc.mXi() - massXi) > mDeltaMassXi) {
    return false;
  }

  // V0 mass
  if (std::fabs(casc.mLambda() - massLambda) > mDeltaMassLambdaFromXi) {
    return false;
  }

  // PID
  float nSigmaPrTpc[3] = {-999., dauTracks[1].tpcNSigmaPr(), dauTracks[2].tpcNSigmaPr()};
  float nSigmaPrTof[3] = {-999., dauTracks[1].tofNSigmaPr(), dauTracks[2].tofNSigmaPr()};
  float nSigmaPiTpc[3] = {dauTracks[0].tpcNSigmaPi(), dauTracks[1].tpcNSigmaPi(), dauTracks[2].tpcNSigmaPi()};
  float nSigmaPiTof[3] = {dauTracks[0].tofNSigmaPi(), dauTracks[1].tofNSigmaPi(), dauTracks[2].tofNSigmaPi()};
  if (mTpcPidCalibrationOption == 1) {
    for (int iDau{0}; iDau < 3; ++iDau) {
      nSigmaPiTpc[iDau] = getTPCPostCalib(dauTracks[iDau], kPi);
      if (iDau == 0) {
        continue;
      }
      nSigmaPrTpc[iDau] = getTPCPostCalib(dauTracks[iDau], kPr);
    }
  } else if (mTpcPidCalibrationOption == 2) {
    for (int iDau{0}; iDau < 3; ++iDau) {
      nSigmaPiTpc[iDau] = getTPCSplineCalib(dauTracks[iDau], (dauTracks[iDau].sign() > 0) ? kPi : kAntiPi);
      if (iDau == 0) {
        continue;
      }
      nSigmaPrTpc[iDau] = getTPCSplineCalib(dauTracks[iDau], (dauTracks[iDau].sign() > 0) ? kPr : kAntiPr);
    }
  }

  // PID to V0 tracks
  if (dauTracks[0].sign() < 0) { // Xi-
    if ((dauTracks[1].hasTPC() && std::fabs(nSigmaPrTpc[1]) > mMaxNsigmaXiDau) && (dauTracks[1].hasTOF() && std::fabs(nSigmaPrTof[1]) > mMaxNsigmaXiDau)) {
      return false;
    }
    if ((dauTracks[2].hasTPC() && std::fabs(nSigmaPiTpc[2]) > mMaxNsigmaXiDau) && (dauTracks[2].hasTOF() && std::fabs(nSigmaPiTof[2]) > mMaxNsigmaXiDau)) {
      return false;
    }
  } else if (dauTracks[0].sign() > 0) { // Xi+
    if ((dauTracks[2].hasTPC() && std::fabs(nSigmaPrTpc[2]) > mMaxNsigmaXiDau) && (dauTracks[2].hasTOF() && std::fabs(nSigmaPrTof[2]) > mMaxNsigmaXiDau)) {
      return false;
    }
    if ((dauTracks[1].hasTPC() && std::fabs(nSigmaPiTpc[1]) > mMaxNsigmaXiDau) && (dauTracks[1].hasTOF() && std::fabs(nSigmaPiTof[1]) > mMaxNsigmaXiDau)) {
      return false;
    }
  }

  // bachelor PID
  if ((dauTracks[0].hasTPC() && std::fabs(nSigmaPiTpc[0]) > mMaxNsigmaXiDau) && (dauTracks[0].hasTOF() && std::fabs(nSigmaPiTof[0]) > mMaxNsigmaXiDau)) {
    return false;
  }

  // additional track cuts
  for (const auto& dauTrack : dauTracks) {
    //  TPC clusters selections
    if (dauTrack.tpcNClsFound() < 70) { // TODO: put me as a configurable please
      return false;
    }
    if (dauTrack.tpcNClsCrossedRows() < 70) {
      return false;
    }
    if (dauTrack.tpcCrossedRowsOverFindableCls() < 0.8) {
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
inline int8_t HfFilterHelper::isSelectedBachelorForCharmBaryon(const T& track, const T2& dca)
{
  int8_t retValue{BIT(kPionForCharmBaryon) | BIT(kKaonForCharmBaryon)};

  if (!track.isGlobalTrackWoDCA()) {
    return kRejected;
  }

  if (track.pt() < mPtMinCharmBaryonBachelor) {
    return kRejected;
  }

  auto pTBinTrack = findBin(mPtBinsTracks, track.pt());
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
inline int HfFilterHelper::computeNumberOfCandidates(std::vector<std::vector<T>> indices)
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

/// ML helper methods

/// Iinitialisation of ONNX session
/// \param onnxFile is the onnx file name
/// \param partName is the particle name
/// \param env is the ONNX environment
/// \param sessionOpt is the ONNX session options
/// \param inputShapes is the input shape
/// \param dataType is the data type (1=float, 11=double)
/// \param loadModelsFromCCDB is the flag to decide whether the ONNX file is read from CCDB or not
/// \param ccdbApi is the CCDB API
/// \param mlModelPathCCDB is the model path in CCDB
/// \param timestampCCDB is the CCDB timestamp
/// \return the pointer to the ONNX Ort::Experimental::Session
inline Ort::Experimental::Session* HfFilterHelper::initONNXSession(std::string& onnxFile, std::string partName, Ort::Env& env, Ort::SessionOptions& sessionOpt, std::vector<std::vector<int64_t>>& inputShapes, int& dataType, bool loadModelsFromCCDB, o2::ccdb::CcdbApi& ccdbApi, std::string mlModelPathCCDB, int64_t timestampCCDB)
{
  // hard coded, we do not let the user change this
  sessionOpt.SetIntraOpNumThreads(1);
  sessionOpt.SetInterOpNumThreads(1);
  Ort::Experimental::Session* session = nullptr;

  std::map<std::string, std::string> metadata;
  bool retrieveSuccess = true;
  if (loadModelsFromCCDB) {
    retrieveSuccess = ccdbApi.retrieveBlob(mlModelPathCCDB + partName, ".", metadata, timestampCCDB, false, onnxFile);
  }
  if (retrieveSuccess) {
    session = new Ort::Experimental::Session{env, onnxFile, sessionOpt};
    inputShapes = session->GetInputShapes();
    if (inputShapes[0][0] < 0) {
      LOGF(warning, Form("Model for %s with negative input shape likely because converted with hummingbird, setting it to 1.", partName.data()));
      inputShapes[0][0] = 1;
    }

    Ort::TypeInfo typeInfo = session->GetInputTypeInfo(0);
    auto tensorInfo = typeInfo.GetTensorTypeAndShapeInfo();
    dataType = tensorInfo.GetElementType();
  } else {
    LOG(fatal) << "Error encountered while fetching/loading the ML model from CCDB! Maybe the ML model doesn't exist yet for this runnumber/timestamp?";
  }

  return session;
}

/// Iinitialisation of ONNX session
/// \param inputFeatures is the vector with input features
/// \param session is the ONNX Ort::Experimental::Session
/// \param inputShapes is the input shape
/// \return the array with the three output scores
template <typename T>
inline std::array<T, 3> HfFilterHelper::predictONNX(std::vector<T>& inputFeatures, std::shared_ptr<Ort::Experimental::Session>& session, std::vector<std::vector<int64_t>>& inputShapes)
{
  std::array<T, 3> scores{-1., 2., 2.};
  std::vector<Ort::Value> inputTensor{};
  inputTensor.push_back(Ort::Experimental::Value::CreateTensor<T>(inputFeatures.data(), inputFeatures.size(), inputShapes[0]));

  // double-check the dimensions of the input tensor
  if (inputTensor[0].GetTensorTypeAndShapeInfo().GetShape()[0] > 0) { // vectorial models can have negative shape if the shape is unknown
    assert(inputTensor[0].IsTensor() && inputTensor[0].GetTensorTypeAndShapeInfo().GetShape() == inputShapes[0]);
  }
  try {
    auto outputTensor = session->Run(session->GetInputNames(), inputTensor, session->GetOutputNames());
    assert(outputTensor.size() == session->GetOutputNames().size() && outputTensor[1].IsTensor());
    auto typeInfo = outputTensor[1].GetTensorTypeAndShapeInfo();
    assert(typeInfo.GetElementCount() == 3); // we need multiclass
    scores[0] = outputTensor[1].GetTensorMutableData<T>()[0];
    scores[1] = outputTensor[1].GetTensorMutableData<T>()[1];
    scores[2] = outputTensor[1].GetTensorMutableData<T>()[2];
  } catch (const Ort::Exception& exception) {
    LOG(error) << "Error running model inference: " << exception.what();
  }

  return scores;
}

/// PID postcalibrations

/// load the TPC spline from the CCDB
/// \param ccdbApi is Api for CCDB
/// \param bunchCrossing is the timestamp of bunchcrossing for the run number
/// \param ccdbPaths  are the paths on CCDB for pions, antipions, kaons, antikaons, protons, antiprotons
inline void HfFilterHelper::setValuesBB(o2::ccdb::CcdbApi& ccdbApi, aod::BCsWithTimestamps::iterator const& bunchCrossing, const std::array<std::string, 6>& ccdbPaths)
{
  for (int iSpecie{0u}; iSpecie < 6; ++iSpecie) {
    std::map<std::string, std::string> metadata;
    auto hSpline = ccdbApi.retrieveFromTFileAny<TH1F>(ccdbPaths[iSpecie], metadata, bunchCrossing.timestamp());

    if (!hSpline) {
      LOG(fatal) << "File from CCDB in path " << ccdbPaths[iSpecie] << " was not found for run " << bunchCrossing.runNumber();
    }

    TAxis* axis = hSpline->GetXaxis();
    mBetheBlochPiKaPr[iSpecie] = {static_cast<double>(hSpline->GetBinContent(axis->FindBin("bb1"))),
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
  mHistMapPiPr[0] = nullptr;
  mHistMapPiPr[1] = nullptr;
  mHistMapPiPr[2] = nullptr;
  mHistMapPiPr[3] = nullptr;

  mHistMapPiPr[0] = reinterpret_cast<TH3F*>(calibList->FindObject("mean_map_pion"));
  mHistMapPiPr[1] = reinterpret_cast<TH3F*>(calibList->FindObject("sigma_map_pion"));
  mHistMapPiPr[2] = reinterpret_cast<TH3F*>(calibList->FindObject("mean_map_proton"));
  mHistMapPiPr[3] = reinterpret_cast<TH3F*>(calibList->FindObject("sigma_map_proton"));

  if (!mHistMapPiPr[0] || !mHistMapPiPr[1] || !mHistMapPiPr[0] || !mHistMapPiPr[1]) {
    LOG(fatal) << "Can not find histograms!";
  }
}

/// Basic selection of proton candidates for Lc
/// \param track is a track
/// \param nsigmaTPCProton max NsigmaTPC for proton candidates
/// \param nsigmaTOFProton max NsigmaTOF for proton candidates
/// \return true if track passes all cuts
template <typename T>
inline bool HfFilterHelper::isSelectedProton4CharmBaryons(const T& track, const float& nsigmaTPCProton, const float& nsigmaTOFProton)
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

  if (std::fabs(NSigmaTPC) > nsigmaTPCProton) {
    return false;
  }
  if (track.hasTOF() && std::fabs(NSigmaTOF) > nsigmaTOFProton) {
    return false;
  }

  return true;
}

/// Basic selection of kaon candidates for charm candidates
/// \param track is a track
/// \param nsigmaTPCKaon max NsigmaTPC for kaon candidates
/// \param nsigmaTOFKaon max NsigmaTOF for kaon candidates
/// \return true if track passes all cuts
template <typename T>
inline bool HfFilterHelper::isSelectedKaon4Charm3Prong(const T& track, const float& nsigmaTPCKaon, const float& nsigmaTOFKaon)
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

  if (std::fabs(NSigmaTPC) > mNSigmaTpcKaCutFor3Prongs) {
    return false;
  }
  if (track.hasTOF() && std::fabs(NSigmaTOF) > mNSigmaTofKaCutFor3Prongs) {
    return false;
  }

  return true;
}

/// Update the TPC PID baesd on the spline of particles
/// \param track is a track parameter
/// \param pidSpecies is the particle species to be considered
/// \return updated nsigma value for TPC PID
template <typename T>
inline double HfFilterHelper::getTPCSplineCalib(const T& track, const int& pidSpecies)
{
  float mMassPar{0.};
  if (pidSpecies == kPi || pidSpecies == kAntiPi) {
    mMassPar = massPi;
  } else if (pidSpecies == kKa || pidSpecies == kAntiKa) {
    mMassPar = massKa;
  } else if (pidSpecies == kPr || pidSpecies == kAntiPr) {
    mMassPar = massProton;
  } else {
    LOGP(fatal, "TPC recalibrated Nsigma requested for unknown particle species, return 999");
    return 999.;
  }

  auto bgScaling = 1 / mMassPar;
  double expBethe = tpc::BetheBlochAleph(static_cast<double>(track.tpcInnerParam() * bgScaling), mBetheBlochPiKaPr[pidSpecies][0], mBetheBlochPiKaPr[pidSpecies][1], mBetheBlochPiKaPr[pidSpecies][2], mBetheBlochPiKaPr[pidSpecies][3], mBetheBlochPiKaPr[pidSpecies][4]);
  double expSigma = expBethe * mBetheBlochPiKaPr[pidSpecies][5];
  return static_cast<float>((track.tpcSignal() - expBethe) / expSigma);
}

/// compute TPC postcalibrated nsigma based on calibration histograms from CCDB
/// \param hCalibMean calibration histograms of mean from CCDB
/// \param hCalibSigma calibration histograms of sigma from CCDB
/// \param track is the track
/// \param pidSpecies is the PID species
/// \return the corrected Nsigma value for the PID species
template <typename T>
inline float HfFilterHelper::getTPCPostCalib(const T& track, const int& pidSpecies)
{
  float tpcNCls = track.tpcNClsFound();
  float tpcPin = track.tpcInnerParam();
  float eta = track.eta();
  float tpcNSigma{0.};
  int iHist{0};

  if (pidSpecies == kKa) {
    tpcNSigma = track.tpcNSigmaKa();
    iHist = 0; // same as pions
  } else if (pidSpecies == kPi) {
    tpcNSigma = track.tpcNSigmaPi();
    iHist = 0;
  } else if (pidSpecies == kPr) {
    tpcNSigma = track.tpcNSigmaPr();
    iHist = 2;
  } else {
    LOG(fatal) << "Wrong PID Species be selected, please check!";
  }
  if (!mHistMapPiPr[iHist] || !mHistMapPiPr[iHist + 1]) {
    LOGP(warn, "Postcalibration TPC PID histograms not set. Use default Nsigma values.");
  }

  auto binTPCNCls = mHistMapPiPr[iHist]->GetXaxis()->FindBin(tpcNCls);
  binTPCNCls = (binTPCNCls == 0 ? 1 : binTPCNCls);
  binTPCNCls = std::min(mHistMapPiPr[iHist]->GetXaxis()->GetNbins(), binTPCNCls);
  auto binPin = mHistMapPiPr[iHist]->GetYaxis()->FindBin(tpcPin);
  binPin = (binPin == 0 ? 1 : binPin);
  binPin = std::min(mHistMapPiPr[iHist]->GetYaxis()->GetNbins(), binPin);
  auto binEta = mHistMapPiPr[iHist]->GetZaxis()->FindBin(eta);
  binEta = (binEta == 0 ? 1 : binEta);
  binEta = std::min(mHistMapPiPr[iHist]->GetZaxis()->GetNbins(), binEta);

  auto mean = mHistMapPiPr[iHist]->GetBinContent(binTPCNCls, binPin, binEta);
  auto width = mHistMapPiPr[iHist + 1]->GetBinContent(binTPCNCls, binPin, binEta);

  return (tpcNSigma - mean) / width;
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
