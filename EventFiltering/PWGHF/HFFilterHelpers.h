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
#include "EventFiltering/filterTables.h"

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
  kKaonForCharmBaryon,
  kSoftPionForSigmaC
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
static const std::array<std::string, 2> eventTitles = {"all", "rejected"};
static const std::array<std::string, kNtriggersHF> hfTriggerNames{filtering::HfHighPt2P::columnLabel(), filtering::HfHighPt3P::columnLabel(), filtering::HfBeauty3P::columnLabel(), filtering::HfBeauty4P::columnLabel(), filtering::HfFemto2P::columnLabel(), filtering::HfFemto3P::columnLabel(), filtering::HfDoubleCharm2P::columnLabel(), filtering::HfDoubleCharm3P::columnLabel(), filtering::HfDoubleCharmMix::columnLabel(), filtering::HfV0Charm2P::columnLabel(), filtering::HfV0Charm3P::columnLabel(), filtering::HfCharmBarToXiBach::columnLabel(), filtering::HfSigmaCPPK::columnLabel(), filtering::HfSigmaC0K0::columnLabel(), filtering::HfPhotonCharm2P::columnLabel(), filtering::HfPhotonCharm3P::columnLabel()};

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
constexpr float massSigmaCPlusPlus = o2::constants::physics::MassSigmaCPlusPlus;
constexpr float massSigmaC0 = o2::constants::physics::MassSigmaC0;

static const o2::framework::AxisSpec ptAxis{50, 0.f, 50.f};
static const o2::framework::AxisSpec pAxis{50, 0.f, 10.f};
static const o2::framework::AxisSpec kstarAxis{100, 0.f, 1.f};
static const o2::framework::AxisSpec etaAxis{30, -1.5f, 1.5f};
static const o2::framework::AxisSpec nSigmaAxis{100, -10.f, 10.f};
static const o2::framework::AxisSpec alphaAxis{100, -1.f, 1.f};
static const o2::framework::AxisSpec qtAxis{100, 0.f, 0.25f};
static const o2::framework::AxisSpec bdtAxis{100, 0.f, 1.f};
static const o2::framework::AxisSpec phiAxis{36, 0., o2::constants::math::TwoPI};
static const std::array<o2::framework::AxisSpec, kNCharmParticles + 17> massAxisC = {o2::framework::AxisSpec{100, 1.65f, 2.05f}, o2::framework::AxisSpec{100, 1.65f, 2.05f}, o2::framework::AxisSpec{100, 1.75f, 2.15f}, o2::framework::AxisSpec{100, 2.05f, 2.45f}, o2::framework::AxisSpec{100, 2.25f, 2.65f}, o2::framework::AxisSpec{100, 0.139f, 0.159f}, o2::framework::AxisSpec{100, 0.f, 0.25f}, o2::framework::AxisSpec{100, 0.f, 0.25f}, o2::framework::AxisSpec{200, 0.48f, 0.88f}, o2::framework::AxisSpec{200, 0.48f, 0.88f}, o2::framework::AxisSpec{100, 1.1f, 1.4f}, o2::framework::AxisSpec{100, 1.1f, 1.4f}, o2::framework::AxisSpec{100, 1.1f, 1.4f}, o2::framework::AxisSpec{100, 1.1f, 1.4f}, o2::framework::AxisSpec{170, 0.13f, 0.3f}, o2::framework::AxisSpec{170, 0.13f, 0.3f}, o2::framework::AxisSpec{200, 0.4f, 0.8f}, o2::framework::AxisSpec{200, 0.4f, 0.8f}, o2::framework::AxisSpec{200, 0.4f, 0.8f}, o2::framework::AxisSpec{200, 0.4f, 0.8f}, o2::framework::AxisSpec{100, 2.3f, 2.9f}, o2::framework::AxisSpec{100, 2.3f, 2.9f}};
static const std::array<o2::framework::AxisSpec, kNBeautyParticles> massAxisB = {o2::framework::AxisSpec{240, 4.8f, 6.0f}, o2::framework::AxisSpec{240, 4.8f, 6.0f}, o2::framework::AxisSpec{240, 4.8f, 6.0f}, o2::framework::AxisSpec{240, 4.8f, 6.0f}, o2::framework::AxisSpec{240, 5.0f, 6.2f}, o2::framework::AxisSpec{240, 5.0f, 6.2f}};

// default values for configurables
// channels to trigger on for femto
constexpr int activeFemtoChannels[1][5] = {{1, 1, 1, 1, 0}}; // pD0, pD+, pDs, pLc, pXic
static const std::vector<std::string> labelsColumnsFemtoChannels = {"protonDZero", "protonDPlus", "protonDs", "protonLc", "protonXic"};

// min and max pT for all tracks combined  (except for V0 and cascades)
constexpr float cutsPt[2][6] = {{1., 0.1, 0.8, 0.5, 0.1, 0.2},
                                {100000., 100000., 5., 100000., 100000., 100000.}}; // beauty, D*, femto, SigmaC, Xic*+ -> SigmaC++K-
static const std::vector<std::string> labelsColumnsCutsPt = {"Beauty", "DstarPlus", "Femto", "CharmBaryon", "SoftPiSigmaC", "SoftKaonXicResoToSigmaC"};
static const std::vector<std::string> labelsRowsCutsPt = {"Minimum", "Maximum"};

// PID cuts
constexpr float cutsNsigma[3][6] = {{3., 3., 3., 5., 3., 3.},             // TPC proton from Lc, pi/K from D0, K from 3-prong, femto, pi/K from Xic/Omegac, K from Xic*->SigmaC-Kaon
                                    {3., 3., 3., 2.5, 3., 3.},            // TOF proton from Lc, pi/K from D0, K from 3-prong, femto, pi/K from Xic/Omegac, K from Xic*->SigmaC-Kaon
                                    {999., 999., 999., 2.5, 999., 999.}}; // Sum in quadrature of TPC and TOF (used only for femto for pT < 4 GeV/c)
static const std::vector<std::string> labelsColumnsNsigma = {"PrFromLc", "PiKaFromDZero", "KaFrom3Prong", "Femto", "PiKaFromCharmBaryon", "SoftKaonFromXicResoToSigmaC"};
static const std::vector<std::string> labelsRowsNsigma = {"TPC", "TOF", "Comb"};

// high pt
constexpr float cutsHighPtThresholds[1][2] = {{8., 8.}}; // 2-prongs, 3-prongs
static const std::vector<std::string> labelsColumnsHighPtThresholds = {"2Prongs", "3Prongs"};

// beauty
constexpr float cutsDeltaMassB[1][kNBeautyParticles] = {{0.4, 0.4, 0.4, 0.4, 0.4, 0.4}}; // B+, B0, B0toDstar, Bs, Lb, Xib
static const std::vector<std::string> labelsColumnsDeltaMassB = {"Bplus", "BZero", "BZeroToDstar", "Bs", "Lb", "Xib"};

// double charm
constexpr int activeDoubleCharmChannels[1][3] = {{1, 1, 1}}; // kDoubleCharm2P, kDoubleCharm3P, kDoubleCharmMix
static const std::vector<std::string> labelsColumnsDoubleCharmChannels = {"DoubleCharm2Prong", "DoubleCharm3Prong", "DoubleCharmMix"};

// charm resonances
constexpr float cutsCharmReso[3][11] = {{0.0, 0.0, 0.0, 0.0, 0.4, 0., 0.0, 0.00, 0.21, 0.21, 0.0},
                                        {0.155, 0.3, 0.3, 0.88, 0.88, 1.35, 0.18, 0.18, 0.25, 0.25, 0.8},
                                        {0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 6.0, 0.0, 6.0, 0.0}}; // D*+, D*0, Ds*0, Ds1+, Ds2*+, Xic*->D, SigmaC0, SigmaC++, SigmaC(2520)0, SigmaC(2520)++, Xic*->SigmaC
static const std::vector<std::string> labelsColumnsDeltaMassCharmReso = {"DstarPlus", "DstarZero", "DsStarZero", "Ds1Plus", "Ds2StarPlus", "XicResoToD", "SigmaC0", "SigmaCPlusPlus", "SigmaC02520", "SigmaCPlusPlus2520", "XicResoToSigmaC"};
static const std::vector<std::string> labelsRowsDeltaMassCharmReso = {"deltaMassMin", "deltaMassMax", "ptMin"};
// V0s for charm resonances
constexpr float cutsV0s[1][6] = {{0.85, 0.97, 0.5, 4., 0.02, 0.01}}; // cosPaGamma, cosPaK0sLambda, radiusK0sLambda, nSigmaPrLambda, deltaMassK0S, deltaMassLambda
static const std::vector<std::string> labelsColumnsV0s = {"CosPaGamma", "CosPaK0sLambda", "RadiusK0sLambda", "NSigmaPrLambda", "DeltaMassK0s", "DeltaMassLambda"};

// cascades for Xi + bachelor triggers
constexpr float cutsCascades[1][8] = {{0.2, 1., 0.01, 0.01, 0.99, 0.99, 0.3, 3.}}; // ptXiBachelor, deltaMassXi, deltaMassLambda, cosPaXi, cosPaLambda, DCAxyXi, nSigmaPid
static const std::vector<std::string> labelsColumnsCascades = {"PtBachelor", "PtXi", "DeltaMassXi", "DeltaMassLambda", "CosPAXi", "CosPaLambda", "DCAxyXi", "NsigmaPid"};
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
  void setHighPtTriggerThresholds(float threshold2Prongs, float threshold3Prongs)
  {
    mPtThresholdHighPt2Prongs = threshold2Prongs;
    mPtThresholdHighPt3Prongs = threshold3Prongs;
  }
  void setPtBinsSingleTracks(std::vector<double> ptBins) { mPtBinsTracks = ptBins; }
  void setCutsSingleTrackBeauty(o2::framework::LabeledArray<double> cutsSingleTrack3P, o2::framework::LabeledArray<double> cutsSingleTrack4P)
  {
    mCutsSingleTrackBeauty3Prong = cutsSingleTrack3P;
    mCutsSingleTrackBeauty4Prong = cutsSingleTrack4P;
  }
  void setPtLimitsProtonForFemto(float minPt, float maxPt)
  {
    mPtMinProtonForFemto = minPt;
    mPtMaxProtonForFemto = maxPt;
  }
  void setPtLimitsBeautyBachelor(float minPt, float maxPt)
  {
    mPtMinBeautyBachelor = minPt;
    mPtMaxBeautyBachelor = maxPt;
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
  void setCutsSingleTrackCharmBaryonBachelor(o2::framework::LabeledArray<double> cutsSingleTrack) { mCutsSingleTrackCharmBaryonBachelor = cutsSingleTrack; }
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

  void setTpcPidCalibrationOption(int opt) { mTpcPidCalibrationOption = opt; }

  void setMassResolParametrisation(std::string recoPass)
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
    } else {
      LOGP(fatal, "Mass resolution parametrisation {} not supported! Please set 2023_pass3", recoPass.data());
    }
  }

  void setNumSigmaForDeltaMassCharmHadCut(float nSigma) { mNumSigmaDeltaMassCharmHad = nSigma; }

  // helper functions for selections
  template <typename T>
  bool isSelectedHighPt2Prong(const T& pt);
  template <typename T>
  bool isSelectedHighPt3Prong(const T& pt);
  template <typename T, typename T1, typename T2>
  int8_t isSelectedTrackForSoftPionOrBeauty(const T& track, const T1& trackPar, const T2& dca, const int& whichTrigger);
  template <typename T1, typename T2, typename H2>
  bool isSelectedProton4Femto(const T1& track, const T2& trackPar, const int& activateQA, H2 hProtonTPCPID, H2 hProtonTOFPID, bool forceTof);
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
  template <typename V0, typename Coll, typename T, typename H2>
  int8_t isSelectedV0(const V0& v0, const std::array<T, 2>& dauTracks, const Coll& collision, const int& activateQA, H2 hV0Selected, std::array<H2, 4>& hArmPod);
  template <typename Photon, typename T, typename H2>
  inline bool isSelectedPhoton(const Photon& photon, const std::array<T, 2>& dauTracks, const int& activateQA, H2 hV0Selected, std::array<H2, 4>& hArmPod);
  template <typename Casc, typename T, typename Coll>
  bool isSelectedCascade(const Casc& casc, const std::array<T, 3>& dauTracks, const Coll& collision);
  template <typename T, typename T2>
  int8_t isSelectedBachelorForCharmBaryon(const T& track, const T2& dca);
  template <typename T, typename U>
  int8_t isBDTSelected(const T& scores, const U& thresholdBDTScores);
  template <bool isKaonTrack, typename T>
  bool isSelectedKaonFromXicResoToSigmaC(const T& track);

  // helpers
  template <typename T>
  T computeRelativeMomentum(const std::array<T, 3>& pTrack, const std::array<T, 3>& CharmCandMomentum, const T& CharmMass);
  template <typename T>
  int computeNumberOfCandidates(std::vector<std::vector<T>> indices);

  // PID
  void setValuesBB(o2::ccdb::CcdbApi& ccdbApi, aod::BCsWithTimestamps::iterator const& bunchCrossing, const std::array<std::string, 6>& ccdbPaths);
  void setTpcRecalibMaps(o2::framework::Service<o2::ccdb::BasicCCDBManager> const& ccdb, aod::BCsWithTimestamps::iterator const& bunchCrossing, const std::string& ccdbPath);

 private:
  // selections
  template <typename T>
  bool isSelectedKaon4Charm3Prong(const T& track);
  template <typename T>
  bool isSelectedProton4CharmBaryons(const T& track);

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
  float mPtMinSoftPionForSigmaC{0.1};                                        // minimum pt for the Σ0,++ soft pion
  float mPtMaxSoftPionForSigmaC{10000.f};                                    // maximum pt for the Σ0,++ soft pion
  float mPtMinSoftKaonForXicResoToSigmaC{0.1};                               // minimum pt for the soft kaon of Xic* to SigmaC-Kaon
  float mPtMaxSoftKaonForXicResoToSigmaC{10000.f};                           // maximum pt for the soft kaon of Xic* to SigmaC-Kaon
  float mPtMinBeautyBachelor{0.5};                                           // minimum pt for the b-hadron pion daughter
  float mPtMinProtonForFemto{0.8};                                           // minimum pt for the proton for femto
  float mPtMinCharmBaryonBachelor{0.5};                                      // minimum pt for the bachelor pion from Xic/Omegac decays
  float mPtMaxSoftPionForDstar{2.};                                          // maximum pt for the D*+ soft pion
  float mPtMaxBeautyBachelor{100000.};                                       // maximum pt for the b-hadron pion daughter
  float mPtMaxProtonForFemto{5.0};                                           // maximum pt for the proton for femto
  float mPtMaxCharmBaryonBachelor{100000.};                                  // maximum pt for the bachelor pion from Xic/Omegac decays
  float mPtThresholdPidStrategyForFemto{8.};                                 // pt threshold to change strategy for proton PID for femto
  float mPtMinSigmaCZero{0.f};                                               // pt min SigmaC0 candidate
  float mPtMinSigmaC2520Zero{0.f};                                           // pt min SigmaC(2520)0 candidate
  float mPtMinSigmaCPlusPlus{0.f};                                           // pt min SigmaC++ candidate
  float mPtMinSigmaC2520PlusPlus{0.f};                                       // pt min SigmaC(2520)++ candidate
  std::array<float, 3> mNSigmaPrCutsForFemto{3., 3., 3.};                    // cut values for Nsigma TPC, TOF, combined for femto protons
  float mNSigmaTpcPrCutForCharmBaryons{3.};                                  // maximum Nsigma TPC for protons in Lc and Xic decays
  float mNSigmaTofPrCutForCharmBaryons{3.};                                  // maximum Nsigma TOF for protons in Lc and Xic decays
  float mNSigmaTpcKaCutFor3Prongs{3.};                                       // maximum Nsigma TPC for kaons in 3-prong decays
  float mNSigmaTofKaCutFor3Prongs{3.};                                       // maximum Nsigma TOF for kaons in 3-prong decays
  float mNSigmaTpcPiKaCutForDzero{3.};                                       // maximum Nsigma TPC for pions/kaons in D0 decays
  float mNSigmaTofPiKaCutForDzero{3.};                                       // maximum Nsigma TOF for pions/kaons in D0 decays
  float mDeltaMassMinSigmaCZero{0.155};                                      // minimum delta mass M(pKpipi)-M(pKpi) of SigmaC0 candidates
  float mDeltaMassMaxSigmaCZero{0.18};                                       // maximum delta mass M(pKpipi)-M(pKpi) of SigmaC0 candidates
  float mDeltaMassMinSigmaC2520Zero{0.2};                                    // minimum delta mass M(pKpipi)-M(pKpi) of SigmaC(2520)0 candidates
  float mDeltaMassMaxSigmaC2520Zero{0.26};                                   // maximum delta mass M(pKpipi)-M(pKpi) of SigmaC(2520)0 candidates
  float mDeltaMassMinSigmaCPlusPlus{0.155};                                  // minimum delta mass M(pKpipi)-M(pKpi) of SigmaC++ candidates
  float mDeltaMassMaxSigmaCPlusPlus{0.18};                                   // maximum delta mass M(pKpipi)-M(pKpi) of SigmaC++ candidates
  float mDeltaMassMinSigmaC2520PlusPlus{0.2};                                // minimum delta mass M(pKpipi)-M(pKpi) of SigmaC(2520)++ candidates
  float mDeltaMassMaxSigmaC2520PlusPlus{0.26};                               // maximum delta mass M(pKpipi)-M(pKpi) of SigmaC(2520)++ candidates
  float mMinGammaCosinePa{0.85};                                             // minimum cosp for gammas
  float mMinK0sLambdaCosinePa{0.97};                                         // minimum cosp for K0S and Lambda in charm excited decays
  float mMinK0sLambdaRadius{0.5};                                            // minimum radius for K0S and Lambda in charm excited decays
  float mMaxNsigmaPrForLambda{4.};                                           // maximum Nsigma TPC and TOF for protons in Lambda decays
  float mDeltaMassK0s{0.02};                                                 // delta mass cut for K0S in charm excited decays
  float mDeltaMassLambda{0.01};                                              // delta mass cut for Lambda in charm excited decays
  float mMinPtXiBachelor{0.1};                                               // minimum pt for Xi bachelor in Xic/Omegac decays
  float mMinPtXi{1.};                                                        // minimum pt for Xi in Xic/Omegac decays
  float mDeltaMassXi{0.01};                                                  // delta mass cut for Xi in Xic/Omegac decays
  float mDeltaMassLambdaFromXi{0.01};                                        // delta mass cut for Lambda <- Xi in Xic/Omegac decays
  float mCosPaXi{0.99};                                                      // minimum cosp for Xi in Xic/Omegac decays
  float mCosPaLambdaFromXi{0.99};                                            // minimum cosp for Xi in Xic/Omegac decays
  float mMaxDcaXyXi{0.3};                                                    // maximum dca for Xi in Xic/Omegac decays
  float mMaxNsigmaXiDau{3.};                                                 // maximum Nsigma TPC and TOF for Xi daughter tracks
  o2::framework::LabeledArray<double> mCutsSingleTrackCharmBaryonBachelor{}; // dca selections for the bachelor pion from Xic/Omegac decays
  float mNSigmaTpcPiCharmBaryonBachelor{3.};                                 // maximum Nsigma TPC for pions in Xic/Omegac decays
  float mNSigmaTofPiCharmBaryonBachelor{3.};                                 // maximum Nsigma TOF for pions in Xic/Omegac decays
  float mNumSigmaDeltaMassCharmHad{2.5};                                     // number of sigmas for delta mass cut for charm hadrons in B and charm excited decays
  std::array<float, 2> mSigmaPars2Prongs{};                                  // parameters (intercept, slope) for parametrisation of mass sigma vs pT for 2-prongs
  std::array<float, 2> mDeltaMassPars2Prongs{};                              // parameters (intercept, slope) for parametrisation of mass delta wrt PDG vs pT for 2-prongs
  std::array<float, 2> mSigmaPars3Prongs{};                                  // parameters (intercept, slope) for parametrisation of mass sigma vs pT for 3-prongs
  std::array<float, 2> mDeltaMassPars3Prongs{};                              // parameters (intercept, slope) for parametrisation of mass delta wrt PDG vs pT for 3-prongs
  float mPtThresholdHighPt2Prongs{8.};                                       // threshold for high pT triggers for 2-prongs
  float mPtThresholdHighPt3Prongs{8.};                                       // threshold for high pT triggers for 3-prongs
  float mNSigmaTpcKaonFromXicResoToSigmaC{3.};                               // maximum Nsigma TPC for kaons in Xic*->SigmaC-Kaon
  float mNSigmaTofKaonFromXicResoToSigmaC{3.};                               // maximum Nsigma TOF for kaons in Xic*->SigmaC-Kaon

  // PID recalibrations
  int mTpcPidCalibrationOption{0};                        // Option for TPC PID calibration (0 -> AO2D, 1 -> postcalibrations, 2 -> alternative bethe bloch parametrisation)
  std::array<TH3F*, 6> mHistMapPiPrKa{};                  // Map for TPC PID postcalibrations for pions, kaon and protons
  std::array<std::vector<double>, 6> mBetheBlochPiKaPr{}; // Bethe-Bloch parametrisations for pions, antipions, kaons, antikaons, protons, antiprotons in TPC
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
template <typename T, typename T1, typename T2>
inline int8_t HfFilterHelper::isSelectedTrackForSoftPionOrBeauty(const T& track, const T1& trackPar, const T2& dca, const int& whichTrigger)
{

  int8_t retValue{BIT(kSoftPion) | BIT(kForBeauty) | BIT(kSoftPionForBeauty) | BIT(kSoftPionForSigmaC)};

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

  if (whichTrigger == kSigmaCPPK || whichTrigger == kSigmaC0K0) {

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
  if (pT < mPtMinBeautyBachelor || pT > mPtMaxBeautyBachelor) {
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
/// \param forceTof flag to force TOF PID
/// \return true if track passes all cuts
template <typename T1, typename T2, typename H2>
inline bool HfFilterHelper::isSelectedProton4Femto(const T1& track, const T2& trackPar, const int& activateQA, H2 hProtonTPCPID, H2 hProtonTOFPID, bool forceTof)
{
  float pt = trackPar.getPt();
  if (pt < mPtMinProtonForFemto || pt > mPtMaxProtonForFemto) {
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
  if (!forceTof && !track.hasTOF()) {
    NSigmaTOF = 0.; // always accepted
  }

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
    if (forceTof || track.hasTOF()) {
      hProtonTOFPID->Fill(track.p(), NSigmaTOF);
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
  if (!isSelectedKaon4Charm3Prong(trackOppositeCharge)) {
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
  if (!isSelectedKaon4Charm3Prong(trackOppositeCharge)) {
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
  if (!isSelectedKaon4Charm3Prong(trackOppositeCharge)) {
    return retValue;
  }
  if (isSelectedProton4CharmBaryons(trackSameChargeFirst)) {
    retValue |= BIT(0);
  }
  if (isSelectedProton4CharmBaryons(trackSameChargeSecond)) {
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
/// \param dauTracks is a 2-element array with positive and negative V0 daughter tracks
/// \param collision is the current collision
/// \param activateQA flag to fill QA histos
/// \param hV0Selected is the pointer to the QA histo for selected V0S
/// \param hArmPod is the pointer to an array of QA histo AP plot after selection
/// \return an integer passes all cuts
template <typename V0, typename Coll, typename T, typename H2>
inline int8_t HfFilterHelper::isSelectedV0(const V0& v0, const std::array<T, 2>& dauTracks, const Coll& /*collision*/, const int& activateQA, H2 hV0Selected, std::array<H2, 4>& hArmPod)
{
  int8_t isSelected{BIT(kK0S) | BIT(kLambda) | BIT(kAntiLambda)};

  if (activateQA > 1) {
    for (int iV0{kK0S}; iV0 < kNV0; ++iV0) {
      hV0Selected->Fill(0., iV0);
    }
  }

  // eta of daughters
  if (std::fabs(dauTracks[0].eta()) > 1. || std::fabs(dauTracks[1].eta()) > 1.) { // cut all V0 daughters with |eta| > 1.
    if (activateQA > 1) {
      for (int iV0{kK0S}; iV0 < kNV0; ++iV0) {
        hV0Selected->Fill(1., iV0);
      }
    }
    return kRejected;
  }

  // V0 radius
  if (v0.v0radius() < mMinK0sLambdaRadius) {
    for (int iV0{kK0S}; iV0 < kNV0; ++iV0) {
      CLRBIT(isSelected, iV0);
      if (activateQA > 1) {
        hV0Selected->Fill(2., iV0);
      }
    }
  }

  auto v0CosinePa = v0.v0cosPA();
  for (int iV0{kK0S}; iV0 < kNV0; ++iV0) {
    if (TESTBIT(isSelected, iV0) && v0CosinePa < mMinK0sLambdaCosinePa) {
      CLRBIT(isSelected, iV0);
      if (activateQA > 1) {
        hV0Selected->Fill(3., iV0);
      }
    }
  }

  // armenteros-podolanski / mass
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
    if (TESTBIT(isSelected, iV0) && v0.dcav0topv() > 0.1f) { // we want only primary V0s
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
      hV0Selected->Fill(7., kLambda);
    }
  }
  if (TESTBIT(isSelected, kAntiLambda) && ((dauTracks[1].hasTPC() && std::fabs(nSigmaPrTpc[1]) > mMaxNsigmaPrForLambda) || (dauTracks[1].hasTOF() && std::fabs(nSigmaPrTof[1]) > mMaxNsigmaPrForLambda))) {
    CLRBIT(isSelected, kAntiLambda);
    if (activateQA > 1) {
      hV0Selected->Fill(7., kAntiLambda);
    }
  }

  if (activateQA) {
    for (int iV0{kK0S}; iV0 < kNV0; ++iV0) {
      if (TESTBIT(isSelected, iV0)) {
        hArmPod[iV0]->Fill(v0.alpha(), v0.qtarm());
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
/// \param dauTracks is a 3-element array with bachelor, positive and negative V0 daughter tracks
/// \param collision is the collision
/// \return true if cascade passes all cuts
template <typename Casc, typename T, typename Coll>
inline bool HfFilterHelper::isSelectedCascade(const Casc& casc, const std::array<T, 3>& dauTracks, const Coll& collision)
{

  // Xi min pT
  if (casc.pt() < mMinPtXi) {
    return false;
  }

  // eta of daughters
  if (std::fabs(dauTracks[0].eta()) > 1. || std::fabs(dauTracks[1].eta()) > 1. || std::fabs(dauTracks[2].eta()) > 1.) { // cut all V0 daughters with |eta| > 1.
    return false;
  }

  // V0 radius
  if (casc.v0radius() < 1.2) {
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
  std::array<std::string, 6> mapNames = {"mean_map_pion", "sigma_map_pion", "mean_map_kaon", "sigma_map_kaon", "mean_map_proton", "sigma_map_proton"};

  for (size_t iMap = 0; iMap < mapNames.size(); iMap++) {
    mHistMapPiPrKa[iMap] = nullptr;
  }

  for (size_t iMap = 0; iMap < mapNames.size(); iMap++) {

    mHistMapPiPrKa[iMap] = reinterpret_cast<TH3F*>(calibList->FindObject(mapNames[iMap].data()));
    if (!mHistMapPiPrKa[iMap]) {
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
template <typename T>
inline bool HfFilterHelper::isSelectedProton4CharmBaryons(const T& track)
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

  if (std::fabs(NSigmaTPC) > mNSigmaTpcPrCutForCharmBaryons) {
    return false;
  }
  if (track.hasTOF() && std::fabs(NSigmaTOF) > mNSigmaTofPrCutForCharmBaryons) {
    return false;
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
    return isSelectedKaon4Charm3Prong(track);
  }

  return true;
}

/// Basic selection of kaon candidates for charm candidates
/// \param track is a track
/// \return true if track passes all cuts
template <typename T>
inline bool HfFilterHelper::isSelectedKaon4Charm3Prong(const T& track)
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

  if (pidSpecies == kPi) {
    tpcNSigma = track.tpcNSigmaPi();
    iHist = 0;
  } else if (pidSpecies == kKa) {
    tpcNSigma = track.tpcNSigmaKa();
    iHist = 2;
  } else if (pidSpecies == kPr) {
    tpcNSigma = track.tpcNSigmaPr();
    iHist = 4;
  } else {
    LOG(fatal) << "Wrong PID Species be selected, please check!";
  }
  if (!mHistMapPiPrKa[iHist] || !mHistMapPiPrKa[iHist + 1]) {
    LOGP(warn, "Postcalibration TPC PID histograms not set. Use default Nsigma values.");
  }

  auto binTPCNCls = mHistMapPiPrKa[iHist]->GetXaxis()->FindBin(tpcNCls);
  binTPCNCls = (binTPCNCls == 0 ? 1 : binTPCNCls);
  binTPCNCls = std::min(mHistMapPiPrKa[iHist]->GetXaxis()->GetNbins(), binTPCNCls);
  auto binPin = mHistMapPiPrKa[iHist]->GetYaxis()->FindBin(tpcPin);
  binPin = (binPin == 0 ? 1 : binPin);
  binPin = std::min(mHistMapPiPrKa[iHist]->GetYaxis()->GetNbins(), binPin);
  auto binEta = mHistMapPiPrKa[iHist]->GetZaxis()->FindBin(eta);
  binEta = (binEta == 0 ? 1 : binEta);
  binEta = std::min(mHistMapPiPrKa[iHist]->GetZaxis()->GetNbins(), binEta);

  auto mean = mHistMapPiPrKa[iHist]->GetBinContent(binTPCNCls, binPin, binEta);
  auto width = mHistMapPiPrKa[iHist + 1]->GetBinContent(binTPCNCls, binPin, binEta);

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
