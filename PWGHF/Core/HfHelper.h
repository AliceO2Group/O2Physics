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

/// \file HfHelper.h
/// \brief Class with helper functions for HF analyses
///
/// \author Vít Kučera <vit.kucera@cern.ch>, Inha University

#ifndef PWGHF_CORE_HFHELPER_H_
#define PWGHF_CORE_HFHELPER_H_

#include "PWGHF/Utils/utilsAnalysis.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelectorPID.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>

#include <Math/GenVector/Boost.h>
#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>

#include <array>
#include <cmath>
#include <vector>

template <typename T>
concept IsB0ToDstarPiChannel = requires(T candidate) {
  candidate.prongD0Id();
};

class HfHelper
{
 public:
  /// Default constructor
  HfHelper() = default;

  /// Default destructor
  ~HfHelper() = default;

  // 2-prong

  // D0(bar) → π± K∓

  template <typename T>
  static auto ctD0(const T& candidate)
  {
    return candidate.ct(o2::constants::physics::MassD0);
  }

  template <typename T>
  static auto yD0(const T& candidate)
  {
    return candidate.y(o2::constants::physics::MassD0);
  }

  template <typename T>
  static auto eD0(const T& candidate)
  {
    return candidate.e(o2::constants::physics::MassD0);
  }

  template <typename T>
  static auto invMassD0ToPiK(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassPiPlus, o2::constants::physics::MassKPlus});
  }

  template <typename T>
  static auto invMassD0barToKPi(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassKPlus, o2::constants::physics::MassPiPlus});
  }

  template <typename T>
  static auto cosThetaStarD0(const T& candidate)
  {
    return candidate.cosThetaStar(std::array{o2::constants::physics::MassPiPlus, o2::constants::physics::MassKPlus}, o2::constants::physics::MassD0, 1);
  }

  template <typename T>
  static auto cosThetaStarD0bar(const T& candidate)
  {
    return candidate.cosThetaStar(std::array{o2::constants::physics::MassKPlus, o2::constants::physics::MassPiPlus}, o2::constants::physics::MassD0, 0);
  }

  // J/ψ

  template <typename T>
  static auto ctJpsi(const T& candidate)
  {
    return candidate.ct(o2::constants::physics::MassJPsi);
  }

  template <typename T>
  static auto yJpsi(const T& candidate)
  {
    return candidate.y(o2::constants::physics::MassJPsi);
  }

  template <typename T>
  static auto eJpsi(const T& candidate)
  {
    return candidate.e(o2::constants::physics::MassJPsi);
  }

  // J/ψ → e+ e−
  template <typename T>
  static auto invMassJpsiToEE(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassElectron, o2::constants::physics::MassElectron});
  }
  // J/ψ → μ+ μ−

  template <typename T>
  static auto invMassJpsiToMuMu(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassMuonPlus, o2::constants::physics::MassMuonMinus});
  }

  // hf_cand_casc

  template <typename T>
  static auto invMassLcToK0sP(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassProton, o2::constants::physics::MassK0Short}); // first daughter is bachelor
  }

  template <typename T>
  static auto invMassGammaToEE(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassElectron, o2::constants::physics::MassElectron});
  }

  template <typename T>
  static auto ctV0K0s(const T& candidate)
  {
    return candidate.ctV0(o2::constants::physics::MassK0Short);
  }

  template <typename T>
  static auto ctV0Lambda(const T& candidate)
  {
    return candidate.ctV0(o2::constants::physics::MassLambda0);
  }

  // B± → D0bar(D0) π±

  template <typename T>
  static auto ctBplus(const T& candidate)
  {
    return candidate.ct(o2::constants::physics::MassBPlus);
  }

  template <typename T>
  static auto yBplus(const T& candidate)
  {
    return candidate.y(o2::constants::physics::MassBPlus);
  }

  template <typename T>
  static auto eBplus(const T& candidate)
  {
    return candidate.e(o2::constants::physics::MassBPlus);
  }

  template <typename T>
  static auto invMassBplusToD0Pi(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassD0, o2::constants::physics::MassPiPlus});
  }

  template <typename T>
  static auto invMassBplusToJpsiK(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassMuon, o2::constants::physics::MassMuon, o2::constants::physics::MassKPlus});
  }

  template <typename T>
  static auto cosThetaStarBplus(const T& candidate)
  {
    return candidate.cosThetaStar(std::array{o2::constants::physics::MassD0, o2::constants::physics::MassPiPlus}, o2::constants::physics::MassBPlus, 1);
  }

  // 3-prong

  // D± → π± K∓ π±

  template <typename T>
  static auto ctDplus(const T& candidate)
  {
    return candidate.ct(o2::constants::physics::MassDPlus);
  }

  template <typename T>
  static auto yDplus(const T& candidate)
  {
    return candidate.y(o2::constants::physics::MassDPlus);
  }

  template <typename T>
  static auto eDplus(const T& candidate)
  {
    return candidate.e(o2::constants::physics::MassDPlus);
  }

  template <typename T>
  static auto invMassDplusToPiKPi(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassPiPlus, o2::constants::physics::MassKPlus, o2::constants::physics::MassPiPlus});
  }

  template <typename T>
  static auto invMassDplusToPiKPi(const T& pVec0, const T& pVec1, const T& pVec2)
  {
    return RecoDecay::m(std::array{pVec0, pVec1, pVec2},
                        std::array{o2::constants::physics::MassPiPlus, o2::constants::physics::MassKPlus, o2::constants::physics::MassPiPlus});
  }

  // Ds± → K± K∓ π±

  template <typename T>
  static auto ctDs(const T& candidate)
  {
    return candidate.ct(o2::constants::physics::MassDS);
  }

  template <typename T>
  static auto yDs(const T& candidate)
  {
    return candidate.y(o2::constants::physics::MassDS);
  }

  template <typename T>
  static auto eDs(const T& candidate)
  {
    return candidate.e(o2::constants::physics::MassDS);
  }

  template <typename T>
  static auto invMassDsToKKPi(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassKPlus, o2::constants::physics::MassKPlus, o2::constants::physics::MassPiPlus});
  }

  template <typename T>
  static auto invMassDsToPiKK(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassPiPlus, o2::constants::physics::MassKPlus, o2::constants::physics::MassKPlus});
  }

  template <typename T>
  static auto massKKPairDsToKKPi(const T& candidate)
  {
    return RecoDecay::m(std::array{candidate.pVectorProng0(), candidate.pVectorProng1()}, std::array{o2::constants::physics::MassKPlus, o2::constants::physics::MassKPlus});
  }

  template <typename T>
  static auto massKKPairDsToPiKK(const T& candidate)
  {
    return RecoDecay::m(std::array{candidate.pVectorProng1(), candidate.pVectorProng2()}, std::array{o2::constants::physics::MassKPlus, o2::constants::physics::MassKPlus});
  }

  template <typename T>
  static auto deltaMassPhiDsToKKPi(const T& candidate)
  {
    return std::abs(massKKPairDsToKKPi(candidate) - o2::constants::physics::MassPhi);
  }

  template <typename T>
  static auto deltaMassPhiDsToPiKK(const T& candidate)
  {
    return std::abs(massKKPairDsToPiKK(candidate) - o2::constants::physics::MassPhi);
  }

  /// Calculate the cosine of the angle between the pion and the opposite sign kaon in the phi rest frame
  /// \param candidate Ds candidate from aod::HfCand3Prong table
  /// \param option mass hypothesis considered: 0 = KKPi, 1 = PiKK
  /// \return cosine of pion-kaon angle in the phi rest frame
  template <typename T>
  static auto cosPiKPhiRestFrame(const T& candidate, int option)
  {
    // Ported from AliAODRecoDecayHF3Prong::CosPiKPhiRFrame
    std::array<float, 3> momPi{};
    std::array<float, 3> momK1{};
    std::array<float, 3> momK2{};

    if (option == 0) { // KKPi
      momPi = candidate.pVectorProng2();
      momK1 = candidate.pVectorProng1();
      momK2 = candidate.pVectorProng0();
    } else { // PiKK
      momPi = candidate.pVectorProng0();
      momK1 = candidate.pVectorProng1();
      momK2 = candidate.pVectorProng2();
    }

    ROOT::Math::PxPyPzMVector const vecPi(momPi[0], momPi[1], momPi[2], o2::constants::physics::MassPiPlus);
    ROOT::Math::PxPyPzMVector const vecK1(momK1[0], momK1[1], momK1[2], o2::constants::physics::MassKPlus);
    ROOT::Math::PxPyPzMVector const vecK2(momK2[0], momK2[1], momK2[2], o2::constants::physics::MassKPlus);
    ROOT::Math::PxPyPzMVector const vecPhi = vecK1 + vecK2;

    ROOT::Math::Boost const boostToPhiRestFrame(vecPhi.BoostToCM());
    auto momPiPhiRestFrame = boostToPhiRestFrame(vecPi).Vect();
    auto momK1PhiRestFrame = boostToPhiRestFrame(vecK1).Vect();

    return momPiPhiRestFrame.Dot(momK1PhiRestFrame) / std::sqrt(momPiPhiRestFrame.Mag2() * momK1PhiRestFrame.Mag2());
  }

  template <typename T>
  static auto cos3PiKDsToKKPi(const T& candidate)
  {
    auto cosPiK = cosPiKPhiRestFrame(candidate, 0);
    return cosPiK * cosPiK * cosPiK;
  }

  template <typename T>
  static auto absCos3PiKDsToKKPi(const T& candidate)
  {
    return std::abs(cos3PiKDsToKKPi(candidate));
  }

  template <typename T>
  static auto cos3PiKDsToPiKK(const T& candidate)
  {
    auto cosPiK = cosPiKPhiRestFrame(candidate, 1);
    return cosPiK * cosPiK * cosPiK;
  }

  template <typename T>
  static auto absCos3PiKDsToPiKK(const T& candidate)
  {
    return std::abs(cos3PiKDsToPiKK(candidate));
  }

  // Λc± → p± K∓ π±

  template <typename T>
  static auto ctLc(const T& candidate)
  {
    return candidate.ct(o2::constants::physics::MassLambdaCPlus);
  }

  template <typename T>
  static auto yLc(const T& candidate)
  {
    return candidate.y(o2::constants::physics::MassLambdaCPlus);
  }

  template <typename T>
  static auto eLc(const T& candidate)
  {
    return candidate.e(o2::constants::physics::MassLambdaCPlus);
  }

  template <typename T>
  static auto invMassLcToPKPi(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassProton, o2::constants::physics::MassKPlus, o2::constants::physics::MassPiPlus});
  }

  template <typename T>
  static auto invMassLcToPiKP(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassPiPlus, o2::constants::physics::MassKPlus, o2::constants::physics::MassProton});
  }

  template <typename T>
  static auto invMassKPiPairLcToPKPi(const T& candidate)
  {
    return RecoDecay::m(std::array{candidate.pVectorProng1(), candidate.pVectorProng2()}, std::array{o2::constants::physics::MassKPlus, o2::constants::physics::MassPiPlus});
  }

  template <typename T>
  static auto invMassKPiPairLcToPiKP(const T& candidate)
  {
    return RecoDecay::m(std::array{candidate.pVectorProng1(), candidate.pVectorProng0()}, std::array{o2::constants::physics::MassKPlus, o2::constants::physics::MassPiPlus});
  }

  template <typename T>
  static auto invMassPKPairLcToPKPi(const T& candidate)
  {
    return RecoDecay::m(std::array{candidate.pVectorProng0(), candidate.pVectorProng1()}, std::array{o2::constants::physics::MassProton, o2::constants::physics::MassKPlus});
  }

  template <typename T>
  static auto invMassPKPairLcToPiKP(const T& candidate)
  {
    return RecoDecay::m(std::array{candidate.pVectorProng2(), candidate.pVectorProng1()}, std::array{o2::constants::physics::MassProton, o2::constants::physics::MassKPlus});
  }

  template <typename T>
  static auto invMassPPiPairLcToPKPi(const T& candidate)
  {
    return RecoDecay::m(std::array{candidate.pVectorProng0(), candidate.pVectorProng2()}, std::array{o2::constants::physics::MassProton, o2::constants::physics::MassPiPlus});
  }

  template <typename T>
  static auto invMassPPiPairLcToPiKP(const T& candidate)
  {
    return RecoDecay::m(std::array{candidate.pVectorProng2(), candidate.pVectorProng0()}, std::array{o2::constants::physics::MassProton, o2::constants::physics::MassPiPlus});
  }

  // Cd± → De± K∓ π±

  template <typename T>
  static auto invMassCdToDeKPi(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassDeuteron, o2::constants::physics::MassKPlus, o2::constants::physics::MassPiPlus});
  }

  template <typename T>
  static auto invMassCdToPiKDe(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassPiPlus, o2::constants::physics::MassKPlus, o2::constants::physics::MassDeuteron});
  }

  // Ξc± → p± K∓ π±

  template <typename T>
  static auto ctXic(const T& candidate)
  {
    return candidate.ct(o2::constants::physics::MassXiCPlus);
  }

  template <typename T>
  static auto yXic(const T& candidate)
  {
    return candidate.y(o2::constants::physics::MassXiCPlus);
  }

  template <typename T>
  static auto eXic(const T& candidate)
  {
    return candidate.e(o2::constants::physics::MassXiCPlus);
  }

  template <typename T>
  static auto invMassXicToPKPi(const T& candidate)
  {
    return invMassLcToPKPi(candidate);
  }

  template <typename T>
  static auto invMassXicToPiKP(const T& candidate)
  {
    return invMassLcToPiKP(candidate);
  }

  // hf_cand_casc_lf_2prong

  template <typename T>
  static auto invMassXiczeroToXiPi(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassXiMinus, o2::constants::physics::MassPiPlus});
  }

  template <typename T>
  static auto invMassOmegaczeroToOmegaPi(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassOmegaMinus, o2::constants::physics::MassPiPlus});
  }

  // hf_cand_casc_lf_3prong

  template <typename T>
  static auto invMassXicplusToXiPiPi(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassXiMinus, o2::constants::physics::MassPiPlus, o2::constants::physics::MassPiPlus});
  }

  // hf_cand_x

  // X → Jpsi π+ π-
  template <typename T>
  static auto ctX(const T& candidate)
  {
    return candidate.ct(o2::constants::physics::MassX3872);
  }

  template <typename T>
  static auto yX(const T& candidate)
  {
    return candidate.y(o2::constants::physics::MassX3872);
  }

  template <typename T>
  static auto eX(const T& candidate)
  {
    return candidate.e(o2::constants::physics::MassX3872);
  }

  template <typename T>
  static auto invMassXToJpsiPiPi(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassJPsi, o2::constants::physics::MassPiPlus, o2::constants::physics::MassPiPlus});
  }

  /// Difference between the X mass and the sum of the J/psi and di-pion masses
  template <typename T>
  static auto qX(const T& candidate)
  {
    auto piVec1 = std::array{candidate.pxProng1(), candidate.pyProng1(), candidate.pzProng1()};
    auto piVec2 = std::array{candidate.pxProng2(), candidate.pyProng2(), candidate.pzProng2()};

    auto arrayMomenta = std::array{piVec1, piVec2};
    double massPiPi = RecoDecay::m(arrayMomenta, std::array{o2::constants::physics::MassPiPlus, o2::constants::physics::MassPiPlus});

    double massJpsiPiPi = invMassXToJpsiPiPi(candidate);
    return std::abs(massJpsiPiPi - o2::constants::physics::MassJPsi - massPiPi);
  }

  /// Angular difference between the J/psi and the pion
  template <typename T>
  static auto dRX(const T& candidate, int numPi)
  {
    double etaJpsi = RecoDecay::eta(std::array{candidate.pxProng0(), candidate.pyProng0(), candidate.pzProng0()});
    double phiJpsi = RecoDecay::phi(candidate.pxProng0(), candidate.pyProng0());

    double etaPi, phiPi;

    if (numPi <= 1) {
      etaPi = RecoDecay::eta(std::array{candidate.pxProng1(), candidate.pyProng1(), candidate.pzProng1()});
      phiPi = RecoDecay::phi(candidate.pxProng1(), candidate.pyProng1());
    } else {
      etaPi = RecoDecay::eta(std::array{candidate.pxProng2(), candidate.pyProng2(), candidate.pzProng2()});
      phiPi = RecoDecay::phi(candidate.pxProng2(), candidate.pyProng2());
    }

    double const deltaEta = etaJpsi - etaPi;
    double const deltaPhi = RecoDecay::constrainAngle(phiJpsi - phiPi, -o2::constants::math::PI);

    return RecoDecay::sqrtSumOfSquares(deltaEta, deltaPhi);
  }

  /// Difference in pT between the two pions
  template <typename T>
  static auto balancePtPionsX(const T& candidate)
  {
    double ptPi1 = RecoDecay::pt(candidate.pxProng1(), candidate.pyProng1());
    double ptPi2 = RecoDecay::pt(candidate.pxProng2(), candidate.pyProng2());
    return std::abs(ptPi1 - ptPi2) / (ptPi1 + ptPi2);
  }

  // Ξcc±± → p± K∓ π± π±

  template <typename T>
  static auto ctXicc(const T& candidate)
  {
    return candidate.ct(o2::constants::physics::MassXiCCPlusPlus);
  }

  template <typename T>
  static auto yXicc(const T& candidate)
  {
    return candidate.y(o2::constants::physics::MassXiCCPlusPlus);
  }

  template <typename T>
  static auto eXicc(const T& candidate)
  {
    return candidate.e(o2::constants::physics::MassXiCCPlusPlus);
  }

  template <typename T>
  static auto invMassXiccToXicPi(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassXiCPlus, o2::constants::physics::MassPiPlus});
  }

  // chic → Jpsi gamma

  template <typename T>
  static auto ctChic(const T& candidate)
  {
    return candidate.ct(o2::constants::physics::MassChiC1);
  }

  template <typename T>
  static auto yChic(const T& candidate)
  {
    return candidate.y(o2::constants::physics::MassChiC1);
  }

  template <typename T>
  static auto eChic(const T& candidate)
  {
    return candidate.e(o2::constants::physics::MassChiC1);
  }
  template <typename T>
  static auto invMassChicToJpsiGamma(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassJPsi, 0.});
  }

  // Λb → Λc+ π- → p K- π+ π-

  template <typename T>
  static auto ctLb(const T& candidate)
  {
    return candidate.ct(o2::constants::physics::MassLambdaB0);
  }

  template <typename T>
  static auto yLb(const T& candidate)
  {
    return candidate.y(o2::constants::physics::MassLambdaB0);
  }

  template <typename T>
  static auto eLb(const T& candidate)
  {
    return candidate.e(o2::constants::physics::MassLambdaB0);
  }
  template <typename T>
  static auto invMassLbToLcPi(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassLambdaCPlus, o2::constants::physics::MassPiPlus});
  }

  // B0(B0bar) → D∓ π±

  template <typename T>
  static auto ctB0(const T& candidate)
  {
    return candidate.ct(o2::constants::physics::MassB0);
  }

  template <typename T>
  static auto yB0(const T& candidate)
  {
    return candidate.y(o2::constants::physics::MassB0);
  }

  template <typename T>
  static auto eB0(const T& candidate)
  {
    return candidate.e(o2::constants::physics::MassB0);
  }

  template <typename T>
  static auto invMassB0ToDPi(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassDMinus, o2::constants::physics::MassPiPlus});
  }

  template <IsB0ToDstarPiChannel T>
  static auto invMassB0ToDPi(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassD0, o2::constants::physics::MassPiPlus, o2::constants::physics::MassPiPlus});
  }

  template <typename T>
  static auto cosThetaStarB0(const T& candidate)
  {
    return candidate.cosThetaStar(std::array{o2::constants::physics::MassDMinus, o2::constants::physics::MassPiPlus}, o2::constants::physics::MassB0, 1);
  }

  // Bs(bar) → Ds∓ π±

  template <typename T>
  static auto ctBs(const T& candidate)
  {
    return candidate.ct(o2::constants::physics::MassBS);
  }

  template <typename T>
  static auto yBs(const T& candidate)
  {
    return candidate.y(o2::constants::physics::MassBS);
  }

  template <typename T>
  static auto eBs(const T& candidate)
  {
    return candidate.e(o2::constants::physics::MassBS);
  }

  template <typename T>
  static auto invMassBsToDsPi(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassDSBar, o2::constants::physics::MassPiPlus});
  }

  template <typename T>
  static auto invMassBsToJpsiPhi(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassMuon, o2::constants::physics::MassMuon, o2::constants::physics::MassKPlus, o2::constants::physics::MassKPlus});
  }

  template <typename T>
  static auto cosThetaStarBs(const T& candidate)
  {
    return candidate.cosThetaStar(std::array{o2::constants::physics::MassDSBar, o2::constants::physics::MassPiPlus}, o2::constants::physics::MassBS, 1);
  }

  /// Σc0,++ → Λc+(→pK-π+) π-,+

  /// @brief Sc inv. mass using reco mass for Lc in pKpi and PDG mass for pion
  template <typename T, typename U>
  static auto invMassScRecoLcToPKPi(const T& candidateSc, const U& candidateLc)
  {
    return candidateSc.m(std::array{static_cast<double>(invMassLcToPKPi(candidateLc)), o2::constants::physics::MassPiPlus});
  }

  /// @brief Sc inv. mass using reco mass for Lc in piKp and PDG mass for pion
  template <typename T, typename U>
  static auto invMassScRecoLcToPiKP(const T& candidateSc, const U& candidateLc)
  {
    return candidateSc.m(std::array{static_cast<double>(invMassLcToPiKP(candidateLc)), o2::constants::physics::MassPiPlus});
  }

  template <typename T>
  static auto ySc0(const T& candidate)
  {
    return candidate.y(o2::constants::physics::MassSigmaC0);
  }

  template <typename T>
  static auto yScPlusPlus(const T& candidate)
  {
    return candidate.y(o2::constants::physics::MassSigmaCPlusPlus);
  }

  /// Σc0,++ → Λc+(→K0sP) π-,+
  /// @brief Sc inv. mass using reco mass for Lc in K0sP and PDG mass for pion
  template <typename T, typename U>
  static auto invMassScRecoLcToK0sP(const T& candidateSc, const U& candidateLc)
  {
    return candidateSc.m(std::array{static_cast<double>(invMassLcToK0sP(candidateLc)), o2::constants::physics::MassPiMinus});
  }

  /// Apply topological cuts as defined in SelectorCuts.h
  /// \param candB0 B0 candidate
  /// \param cuts B0 candidate selection per pT bin"
  /// \param binsPt pT bin limits
  /// \return true if candidate passes all selections
  template <typename T1, typename T2, typename T3>
  static bool selectionB0ToDPiTopol(const T1& candB0, const T2& cuts, const T3& binsPt)
  {
    auto ptCandB0 = candB0.pt();
    auto ptD = RecoDecay::pt(candB0.pxProng0(), candB0.pyProng0());
    auto ptPi = RecoDecay::pt(candB0.pxProng1(), candB0.pyProng1());

    int pTBin = o2::analysis::findBin(binsPt, ptCandB0);
    if (pTBin == -1) {
      // LOGF(info, "B0 topol selection failed at getpTBin");
      return false;
    }

    // // check that the candidate pT is within the analysis range
    // if (ptCandB0 < ptCandMin || ptCandB0 >= ptCandMax) {
    //   return false;
    // }

    // B0 mass cut
    if (std::abs(invMassB0ToDPi(candB0) - o2::constants::physics::MassB0) > cuts->get(pTBin, "m")) {
      // Printf("B0 topol selection failed at mass diff check");
      return false;
    }

    // pion pt
    if (ptPi < cuts->get(pTBin, "pT Pi")) {
      return false;
    }

    // D- pt
    if (ptD < cuts->get(pTBin, "pT D")) {
      return false;
    }

    /*
    // D mass cut | already applied in candidateSelectorDplusToPiKPi.cxx
    if (std::abs(invMassDplusToPiKPi(hfCandD) - o2::constants::physics::MassDMinus) > cuts->get(pTBin, "DeltaMD")) {
      return false;
    }
    */

    // B0 Decay length
    if (candB0.decayLength() < cuts->get(pTBin, "B0 decLen")) {
      return false;
    }

    // B0 Decay length XY
    if (candB0.decayLengthXY() < cuts->get(pTBin, "B0 decLenXY")) {
      return false;
    }

    // B0 chi2PCA cut
    if (candB0.chi2PCA() > cuts->get(pTBin, "Chi2PCA")) {
      return false;
    }

    // B0 CPA cut
    if (candB0.cpa() < cuts->get(pTBin, "CPA")) {
      return false;
    }

    // d0 of pi
    if (std::abs(candB0.impactParameter1()) < cuts->get(pTBin, "d0 Pi")) {
      return false;
    }

    // d0 of D
    if (std::abs(candB0.impactParameter0()) < cuts->get(pTBin, "d0 D")) {
      return false;
    }

    return true;
  }

  /// Apply PID selection
  /// \param pidTrackPi PID status of trackPi (prong1 of B0 candidate)
  /// \param acceptPIDNotApplicable switch to accept Status::NotApplicable
  /// \return true if prong1 of B0 candidate passes all selections
  static bool selectionB0ToDPiPid(const int pidTrackPi, const bool acceptPIDNotApplicable)
  {
    if (!acceptPIDNotApplicable && pidTrackPi != TrackSelectorPID::Accepted) {
      return false;
    }
    if (acceptPIDNotApplicable && pidTrackPi == TrackSelectorPID::Rejected) {
      return false;
    }

    return true;
  }

  // Apply topological cuts as defined in SelectorCuts.h
  /// \param candBp B+ candidate
  /// \param cuts B+ candidate selection per pT bin"
  /// \param binsPt pT bin limits
  /// \return true if candidate passes all selections
  template <typename T1, typename T2, typename T3>
  static bool selectionBplusToD0PiTopol(const T1& candBp, const T2& cuts, const T3& binsPt)
  {
    auto ptCandBp = candBp.pt();
    auto ptPi = RecoDecay::pt(candBp.pxProng1(), candBp.pyProng1());

    int pTBin = o2::analysis::findBin(binsPt, ptCandBp);
    if (pTBin == -1) {
      return false;
    }

    // B+ mass cut
    if (std::abs(invMassBplusToD0Pi(candBp) - o2::constants::physics::MassBPlus) > cuts->get(pTBin, "m")) {
      return false;
    }

    // pion pt
    if (ptPi < cuts->get(pTBin, "pT Pi")) {
      return false;
    }

    // d0(D0)xd0(pi)
    if (candBp.impactParameterProduct() > cuts->get(pTBin, "Imp. Par. Product")) {
      return false;
    }

    // B Decay length
    if (candBp.decayLength() < cuts->get(pTBin, "B decLen")) {
      return false;
    }

    // B Decay length XY
    if (candBp.decayLengthXY() < cuts->get(pTBin, "B decLenXY")) {
      return false;
    }

    // B+ CPA cut
    if (candBp.cpa() < cuts->get(pTBin, "CPA")) {
      return false;
    }

    // d0 of pi
    if (std::abs(candBp.impactParameter1()) < cuts->get(pTBin, "d0 Pi")) {
      return false;
    }

    // d0 of D
    if (std::abs(candBp.impactParameter0()) < cuts->get(pTBin, "d0 D0")) {
      return false;
    }

    return true;
  }

  /// Apply PID selection
  /// \param pidTrackPi PID status of trackPi (prong1 of B+ candidate)
  /// \param acceptPIDNotApplicable switch to accept Status::NotApplicable
  /// \return true if prong1 of B+ candidate passes all selections
  static bool selectionBplusToD0PiPid(const int pidTrackPi, const bool acceptPIDNotApplicable)
  {
    if (!acceptPIDNotApplicable && pidTrackPi != TrackSelectorPID::Accepted) {
      return false;
    }
    if (acceptPIDNotApplicable && pidTrackPi == TrackSelectorPID::Rejected) {
      return false;
    }

    return true;
  }

  // Apply topological cuts as defined in SelectorCuts.h
  /// \param candBp B+ candidate
  /// \param cuts B+ candidate selection per pT bin
  /// \param binsPt pT bin limits
  /// \return true if candidate passes all selections
  template <typename T1, typename T2, typename T3>
  static bool selectionBplusToJpsiKTopol(const T1& candBp, const T2& cuts, const T3& binsPt)
  {
    auto ptCandBp = candBp.pt();
    auto mCandBp = invMassBplusToJpsiK(candBp);
    auto ptJpsi = RecoDecay::pt(candBp.pxProng0(), candBp.pyProng0());
    auto ptKa = RecoDecay::pt(candBp.pxProng1(), candBp.pyProng1());
    auto candJpsi = candBp.jpsi();
    float pseudoPropDecLen = candBp.decayLengthXY() * mCandBp / ptCandBp;

    int binPt = o2::analysis::findBin(binsPt, ptCandBp);
    if (binPt == -1) {
      return false;
    }

    // B+ mass cut
    if (std::abs(mCandBp - o2::constants::physics::MassBPlus) > cuts->get(binPt, "m")) {
      return false;
    }

    // kaon pt
    if (ptKa < cuts->get(binPt, "pT K")) {
      return false;
    }

    // J/Psi pt
    if (ptJpsi < cuts->get(binPt, "pT J/Psi")) {
      return false;
    }

    // J/Psi mass
    if (std::abs(candJpsi.m() - o2::constants::physics::MassJPsi) < cuts->get(binPt, "DeltaM J/Psi")) {
      return false;
    }

    // d0(J/Psi)xd0(K)
    if (candBp.impactParameterProduct() > cuts->get(binPt, "B Imp. Par. Product")) {
      return false;
    }

    // B+ Decay length
    if (candBp.decayLength() < cuts->get(binPt, "B decLen")) {
      return false;
    }

    // B+ Decay length XY
    if (candBp.decayLengthXY() < cuts->get(binPt, "B decLenXY")) {
      return false;
    }

    // B+ CPA cut
    if (candBp.cpa() < cuts->get(binPt, "CPA")) {
      return false;
    }

    // B+ CPAXY cut
    if (candBp.cpaXY() < cuts->get(binPt, "CPAXY")) {
      return false;
    }

    // d0 of K
    if (std::abs(candBp.impactParameter1()) < cuts->get(binPt, "d0 K")) {
      return false;
    }

    // d0 of J/Psi
    if (std::abs(candBp.impactParameter0()) < cuts->get(binPt, "d0 J/Psi")) {
      return false;
    }

    // B pseudoproper decay length
    if (pseudoPropDecLen < cuts->get(binPt, "B pseudoprop. decLen")) {
      return false;
    }

    return true;
  }

  /// Apply PID selection
  /// \param pidTrackKa PID status of trackKa (prong1 of B+ candidate)
  /// \param acceptPIDNotApplicable switch to accept Status::NotApplicable
  /// \return true if prong1 of B+ candidate passes all selections
  static bool selectionBplusToJpsiKPid(const int pidTrackKa, const bool acceptPIDNotApplicable)
  {
    if (!acceptPIDNotApplicable && pidTrackKa != TrackSelectorPID::Accepted) {
      return false;
    }
    if (acceptPIDNotApplicable && pidTrackKa == TrackSelectorPID::Rejected) {
      return false;
    }

    return true;
  }

  /// Apply topological cuts as defined in SelectorCuts.h
  /// \param candBs Bs candidate
  /// \param cuts Bs candidate selections
  /// \param binsPt pT bin limits
  /// \return true if candidate passes all selections
  template <typename T1, typename T2, typename T3>
  static bool selectionBsToDsPiTopol(const T1& candBs, const T2& cuts, const T3& binsPt)
  {
    auto ptCandBs = candBs.pt();
    auto ptDs = RecoDecay::pt(candBs.pxProng0(), candBs.pyProng0());
    auto ptPi = RecoDecay::pt(candBs.pxProng1(), candBs.pyProng1());

    int pTBin = o2::analysis::findBin(binsPt, ptCandBs);
    if (pTBin == -1) {
      return false;
    }

    // Bs mass cut
    if (std::abs(invMassBsToDsPi(candBs) - o2::constants::physics::MassBS) > cuts->get(pTBin, "m")) {
      return false;
    }

    // pion pt
    if (ptPi < cuts->get(pTBin, "pT Pi")) {
      return false;
    }

    // Ds pt
    if (ptDs < cuts->get(pTBin, "pT Ds")) {
      return false;
    }

    // Bs Decay length
    if (candBs.decayLength() < cuts->get(pTBin, "Bs decLen")) {
      return false;
    }

    // Bs Decay length XY
    if (candBs.decayLengthXY() < cuts->get(pTBin, "Bs decLenXY")) {
      return false;
    }

    // Bs chi2PCA cut
    if (candBs.chi2PCA() > cuts->get(pTBin, "Chi2PCA")) {
      return false;
    }

    // Bs CPA cut
    if (candBs.cpa() < cuts->get(pTBin, "CPA")) {
      return false;
    }

    // d0 of pi
    if (std::abs(candBs.impactParameter1()) < cuts->get(pTBin, "d0 Pi")) {
      return false;
    }

    // d0 of Ds
    if (std::abs(candBs.impactParameter0()) < cuts->get(pTBin, "d0 Ds")) {
      return false;
    }

    // d0(Ds)xd0(pi)
    if (candBs.impactParameterProduct() > cuts->get(pTBin, "Imp. Par. Product")) {
      return false;
    }

    return true;
  }

  /// Apply PID selection
  /// \param pidTrackPi PID status of trackPi (prong1 of Bs candidate)
  /// \param acceptPIDNotApplicable switch to accept Status::NotApplicable
  /// \return true if prong1 of Bs candidate passes all selections
  static bool selectionBsToDsPiPid(const int pidTrackPi, const bool acceptPIDNotApplicable)
  {
    if (!acceptPIDNotApplicable && pidTrackPi != TrackSelectorPID::Accepted) {
      return false;
    }
    if (acceptPIDNotApplicable && pidTrackPi == TrackSelectorPID::Rejected) {
      return false;
    }

    return true;
  }

  // Apply topological cuts as defined in SelectorCuts.h
  /// \param candBs Bs candidate
  /// \param candKa0 kaon candidate 0 (phi daughter)
  /// \param candKa1 kaon candidate 1 (phi daughter)
  /// \param cuts Bs candidate selection per pT bin
  /// \param binsPt pT bin limits
  /// \return true if candidate passes all selections
  template <typename T1, typename T2, typename T3, typename T4, typename T5>
  static bool selectionBsToJpsiPhiTopol(const T1& candBs, const T2& candKa0, const T3& candKa1, const T4& cuts, const T5& binsPt)
  {
    auto ptCandBs = candBs.pt();
    auto mCandBs = invMassBsToJpsiPhi(candBs);
    std::array<float, 3> pVecKa0 = candKa0.pVector();
    std::array<float, 3> pVecKa1 = candKa1.pVector();
    auto mCandPhi = RecoDecay::m(std::array{pVecKa0, pVecKa1}, std::array{o2::constants::physics::MassKPlus, o2::constants::physics::MassKPlus});
    auto ptJpsi = RecoDecay::pt(candBs.pxProng0(), candBs.pyProng0());
    auto candJpsi = candBs.jpsi();
    float pseudoPropDecLen = candBs.decayLengthXY() * mCandBs / ptCandBs;

    int binPt = o2::analysis::findBin(binsPt, ptCandBs);
    if (binPt == -1) {
      return false;
    }

    // Bs mass cut
    if (std::abs(mCandBs - o2::constants::physics::MassBPlus) > cuts->get(binPt, "m")) {
      return false;
    }

    // kaon pt
    if (candKa0.pt() < cuts->get(binPt, "pT K") &&
        candKa1.pt() < cuts->get(binPt, "pT K")) {
      return false;
    }

    // J/Psi pt
    if (ptJpsi < cuts->get(binPt, "pT J/Psi")) {
      return false;
    }

    // phi mass
    if (std::abs(mCandPhi - o2::constants::physics::MassPhi) < cuts->get(binPt, "DeltaM phi")) {
      return false;
    }

    // J/Psi mass
    if (std::abs(candJpsi.m() - o2::constants::physics::MassJPsi) < cuts->get(binPt, "DeltaM J/Psi")) {
      return false;
    }

    // d0(J/Psi)xd0(phi)
    if (candBs.impactParameterProduct() > cuts->get(binPt, "B Imp. Par. Product")) {
      return false;
    }

    // Bs Decay length
    if (candBs.decayLength() < cuts->get(binPt, "B decLen")) {
      return false;
    }

    // Bs Decay length XY
    if (candBs.decayLengthXY() < cuts->get(binPt, "B decLenXY")) {
      return false;
    }

    // Bs CPA cut
    if (candBs.cpa() < cuts->get(binPt, "CPA")) {
      return false;
    }

    // Bs CPAXY cut
    if (candBs.cpaXY() < cuts->get(binPt, "CPAXY")) {
      return false;
    }

    // d0 of phi
    if (std::abs(candBs.impactParameter1()) < cuts->get(binPt, "d0 phi")) {
      return false;
    }

    // d0 of J/Psi
    if (std::abs(candBs.impactParameter0()) < cuts->get(binPt, "d0 J/Psi")) {
      return false;
    }

    // B pseudoproper decay length
    if (pseudoPropDecLen < cuts->get(binPt, "B pseudoprop. decLen")) {
      return false;
    }

    return true;
  }

  /// Apply PID selection
  /// \param pidTrackKa PID status of trackKa (prong1 of B+ candidate)
  /// \param acceptPIDNotApplicable switch to accept Status::NotApplicable
  /// \return true if prong1 of B+ candidate passes all selections
  static bool selectionBsToJpsiPhiPid(const int pidTrackKa, const bool acceptPIDNotApplicable)
  {
    if (!acceptPIDNotApplicable && pidTrackKa != TrackSelectorPID::Accepted) {
      return false;
    }
    if (acceptPIDNotApplicable && pidTrackKa == TrackSelectorPID::Rejected) {
      return false;
    }

    return true;
  }

  /// Apply topological cuts as defined in SelectorCuts.h
  /// \param candLb Lb candidate
  /// \param cuts Lb candidate selection per pT bin"
  /// \param binsPt pT bin limits
  /// \return true if candidate passes all selections
  template <typename T1, typename T2, typename T3>
  static bool selectionLbToLcPiTopol(const T1& candLb, const T2& cuts, const T3& binsPt)
  {
    auto ptCandLb = candLb.pt();
    auto ptLc = candLb.ptProng0();
    auto ptPi = candLb.ptProng1();

    int pTBin = o2::analysis::findBin(binsPt, ptCandLb);
    if (pTBin == -1) {
      return false;
    }

    // Lb mass cut
    if (std::abs(invMassLbToLcPi(candLb) - o2::constants::physics::MassLambdaB0) > cuts->get(pTBin, "m")) {
      return false;
    }

    // pion pt
    if (ptPi < cuts->get(pTBin, "pT Pi")) {
      return false;
    }

    // Lc pt
    if (ptLc < cuts->get(pTBin, "pT Lc+")) {
      return false;
    }

    // Lb Decay length
    if (candLb.decayLength() < cuts->get(pTBin, "Lb decLen")) {
      return false;
    }

    // Lb Decay length XY
    if (candLb.decayLengthXY() < cuts->get(pTBin, "Lb decLenXY")) {
      return false;
    }

    // Lb chi2PCA cut
    if (candLb.chi2PCA() > cuts->get(pTBin, "Chi2PCA")) {
      return false;
    }

    // Lb CPA cut
    if (candLb.cpa() < cuts->get(pTBin, "CPA")) {
      return false;
    }

    // d0 of pi
    if (std::abs(candLb.impactParameter1()) < cuts->get(pTBin, "d0 Pi")) {
      return false;
    }

    // d0 of Lc
    if (std::abs(candLb.impactParameter0()) < cuts->get(pTBin, "d0 Lc+")) {
      return false;
    }

    return true;
  }

  /// \param pidTrackPi PID status of trackPi (prong1 of Lb candidate)
  /// \param acceptPIDNotApplicable switch to accept Status::NotApplicable
  /// \return true if prong1 of Lb candidate passes all selections
  static bool selectionLbToLcPiPid(const int pidTrackPi, const bool acceptPIDNotApplicable)
  {
    if (!acceptPIDNotApplicable && pidTrackPi != TrackSelectorPID::Accepted) {
      return false;
    }
    if (acceptPIDNotApplicable && pidTrackPi == TrackSelectorPID::Rejected) {
      return false;
    }

    return true;
  }

  /// Apply selection on ML scores for charm-hadron daughter in b-hadron decays (common for all the beauty channels)
  /// \param cuts ML score selection per bin of charm-hadron pT
  /// \param binsPtC pT bin limits of charm hadron
  /// \param mlScores vector with ml scores of charm hadron (position 0:bkg 1:prompt 2:nonprompt)
  /// \return true if b-hadron candidate passes all selections
  template <typename T1, typename T2>
  static bool applySelectionDmesMlScoresForB(const T1& cuts, const T2& binsPtC, float ptC, const std::vector<float>& mlScores)
  {
    int pTBin = o2::analysis::findBin(binsPtC, ptC);
    if (pTBin == -1) {
      return false;
    }

    if (mlScores[0] > cuts->get(pTBin, "ML score charm bkg")) {
      return false;
    }

    if (mlScores[1] > cuts->get(pTBin, "ML score charm prompt")) { // we want non-prompt for beauty
      return false;
    }

    if (mlScores[2] < cuts->get(pTBin, "ML score charm nonprompt")) { // we want non-prompt for beauty
      return false;
    }

    return true;
  }

  /// Apply selection on ML scores for charm-hadron daughter in b-hadron decays (could be common for all the beauty channels)
  /// \param candB b-hadron candidates
  /// \param cuts ML score selection per bin of charm-hadron pT
  /// \param binsPtC pT bin limits of charm hadron
  /// \return true if b-hadron candidate passes all selections
  template <typename T1, typename T2, typename T3>
  static bool selectionDmesMlScoresForB(const T1& candD, const T2& cuts, const T3& binsPtC, const std::vector<float>& mlScores)
  {
    return applySelectionDmesMlScoresForB(cuts, binsPtC, candD.pt(), mlScores);
  }

  /// Apply selection on ML scores for charm-hadron daughter in b-hadron decays in reduced format (common for all the beauty channels)
  /// \param candB b-hadron candidates
  /// \param cuts ML score selection per bin of charm-hadron pT
  /// \param binsPtC pT bin limits of charm hadron
  /// \return true if b-hadron candidate passes all selections
  template <typename T1, typename T2, typename T3>
  static bool selectionDmesMlScoresForBReduced(const T1& candB, const T2& cuts, const T3& binsPtC)
  {
    std::vector<float> mlScores;
    mlScores.push_back(candB.prong0MlScoreBkg());
    mlScores.push_back(candB.prong0MlScorePrompt());
    mlScores.push_back(candB.prong0MlScoreNonprompt()); // we want non-prompt for beauty
    return applySelectionDmesMlScoresForB(cuts, binsPtC, RecoDecay::pt(candB.pxProng0(), candB.pyProng0()), mlScores);
  }
};

#endif // PWGHF_CORE_HFHELPER_H_
