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

#include <vector>

#include "Math/GenVector/Boost.h"
#include "Math/Vector4D.h"
#include <TPDGCode.h>

#include "CommonConstants/PhysicsConstants.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelectorPID.h"

#include "PWGHF/Utils/utilsAnalysis.h"

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
  auto ctD0(const T& candidate)
  {
    return candidate.ct(o2::constants::physics::MassD0);
  }

  template <typename T>
  auto yD0(const T& candidate)
  {
    return candidate.y(o2::constants::physics::MassD0);
  }

  template <typename T>
  auto eD0(const T& candidate)
  {
    return candidate.e(o2::constants::physics::MassD0);
  }

  template <typename T>
  auto invMassD0ToPiK(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassPiPlus, o2::constants::physics::MassKPlus});
  }

  template <typename T>
  auto invMassD0barToKPi(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassKPlus, o2::constants::physics::MassPiPlus});
  }

  template <typename T>
  auto cosThetaStarD0(const T& candidate)
  {
    return candidate.cosThetaStar(std::array{o2::constants::physics::MassPiPlus, o2::constants::physics::MassKPlus}, o2::constants::physics::MassD0, 1);
  }

  template <typename T>
  auto cosThetaStarD0bar(const T& candidate)
  {
    return candidate.cosThetaStar(std::array{o2::constants::physics::MassKPlus, o2::constants::physics::MassPiPlus}, o2::constants::physics::MassD0, 0);
  }

  // J/ψ

  template <typename T>
  auto ctJpsi(const T& candidate)
  {
    return candidate.ct(o2::constants::physics::MassJPsi);
  }

  template <typename T>
  auto yJpsi(const T& candidate)
  {
    return candidate.y(o2::constants::physics::MassJPsi);
  }

  template <typename T>
  auto eJpsi(const T& candidate)
  {
    return candidate.e(o2::constants::physics::MassJPsi);
  }

  // J/ψ → e+ e−
  template <typename T>
  auto invMassJpsiToEE(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassElectron, o2::constants::physics::MassElectron});
  }
  // J/ψ → μ+ μ−

  template <typename T>
  auto invMassJpsiToMuMu(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassMuonPlus, o2::constants::physics::MassMuonMinus});
  }

  // hf_cand_casc

  template <typename T>
  auto invMassLcToK0sP(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassProton, o2::constants::physics::MassK0Short}); // first daughter is bachelor
  }

  template <typename T>
  auto invMassGammaToEE(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassElectron, o2::constants::physics::MassElectron});
  }

  template <typename T>
  auto ctV0K0s(const T& candidate)
  {
    return candidate.ctV0(o2::constants::physics::MassK0Short);
  }

  template <typename T>
  auto ctV0Lambda(const T& candidate)
  {
    return candidate.ctV0(o2::constants::physics::MassLambda0);
  }

  // B± → D0bar(D0) π±

  template <typename T>
  auto ctBplus(const T& candidate)
  {
    return candidate.ct(o2::constants::physics::MassBPlus);
  }

  template <typename T>
  auto yBplus(const T& candidate)
  {
    return candidate.y(o2::constants::physics::MassBPlus);
  }

  template <typename T>
  auto eBplus(const T& candidate)
  {
    return candidate.e(o2::constants::physics::MassBPlus);
  }

  template <typename T>
  auto invMassBplusToD0Pi(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassD0, o2::constants::physics::MassPiPlus});
  }

  template <typename T>
  auto cosThetaStarBplus(const T& candidate)
  {
    return candidate.cosThetaStar(std::array{o2::constants::physics::MassD0, o2::constants::physics::MassPiPlus}, o2::constants::physics::MassBPlus, 1);
  }

  // 3-prong

  // D± → π± K∓ π±

  template <typename T>
  auto ctDplus(const T& candidate)
  {
    return candidate.ct(o2::constants::physics::MassDPlus);
  }

  template <typename T>
  auto yDplus(const T& candidate)
  {
    return candidate.y(o2::constants::physics::MassDPlus);
  }

  template <typename T>
  auto eDplus(const T& candidate)
  {
    return candidate.e(o2::constants::physics::MassDPlus);
  }

  template <typename T>
  auto invMassDplusToPiKPi(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassPiPlus, o2::constants::physics::MassKPlus, o2::constants::physics::MassPiPlus});
  }

  template <typename T>
  auto invMassDplusToPiKPi(const T& pVec0, const T& pVec1, const T& pVec2)
  {
    return RecoDecay::m(std::array{pVec0, pVec1, pVec2},
                        std::array{o2::constants::physics::MassPiPlus, o2::constants::physics::MassKPlus, o2::constants::physics::MassPiPlus});
  }

  // Ds± → K± K∓ π±

  template <typename T>
  auto ctDs(const T& candidate)
  {
    return candidate.ct(o2::constants::physics::MassDS);
  }

  template <typename T>
  auto yDs(const T& candidate)
  {
    return candidate.y(o2::constants::physics::MassDS);
  }

  template <typename T>
  auto eDs(const T& candidate)
  {
    return candidate.e(o2::constants::physics::MassDS);
  }

  template <typename T>
  auto invMassDsToKKPi(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassKPlus, o2::constants::physics::MassKPlus, o2::constants::physics::MassPiPlus});
  }

  template <typename T>
  auto invMassDsToPiKK(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassPiPlus, o2::constants::physics::MassKPlus, o2::constants::physics::MassKPlus});
  }

  template <typename T>
  auto massKKPairDsToKKPi(const T& candidate)
  {
    return RecoDecay::m(std::array{candidate.pVectorProng0(), candidate.pVectorProng1()}, std::array{o2::constants::physics::MassKPlus, o2::constants::physics::MassKPlus});
  }

  template <typename T>
  auto massKKPairDsToPiKK(const T& candidate)
  {
    return RecoDecay::m(std::array{candidate.pVectorProng1(), candidate.pVectorProng2()}, std::array{o2::constants::physics::MassKPlus, o2::constants::physics::MassKPlus});
  }

  template <typename T>
  auto deltaMassPhiDsToKKPi(const T& candidate)
  {
    return std::abs(massKKPairDsToKKPi(candidate) - o2::constants::physics::MassPhi);
  }

  template <typename T>
  auto deltaMassPhiDsToPiKK(const T& candidate)
  {
    return std::abs(massKKPairDsToPiKK(candidate) - o2::constants::physics::MassPhi);
  }

  /// Calculate the cosine of the angle between the pion and the opposite sign kaon in the phi rest frame
  /// \param candidate Ds candidate from aod::HfCand3Prong table
  /// \param option mass hypothesis considered: 0 = KKPi, 1 = PiKK
  /// \return cosine of pion-kaon angle in the phi rest frame
  template <typename T>
  auto cosPiKPhiRestFrame(const T& candidate, int option)
  {
    // Ported from AliAODRecoDecayHF3Prong::CosPiKPhiRFrame
    std::array<float, 3> momPi;
    std::array<float, 3> momK1;
    std::array<float, 3> momK2;

    if (option == 0) { // KKPi
      momPi = candidate.pVectorProng2();
      momK1 = candidate.pVectorProng1();
      momK2 = candidate.pVectorProng0();
    } else { // PiKK
      momPi = candidate.pVectorProng0();
      momK1 = candidate.pVectorProng1();
      momK2 = candidate.pVectorProng2();
    }

    ROOT::Math::PxPyPzMVector vecPi(momPi[0], momPi[1], momPi[2], o2::constants::physics::MassPiPlus);
    ROOT::Math::PxPyPzMVector vecK1(momK1[0], momK1[1], momK1[2], o2::constants::physics::MassKPlus);
    ROOT::Math::PxPyPzMVector vecK2(momK2[0], momK2[1], momK2[2], o2::constants::physics::MassKPlus);
    ROOT::Math::PxPyPzMVector vecPhi = vecK1 + vecK2;

    ROOT::Math::Boost boostToPhiRestFrame(vecPhi.BoostToCM());
    auto momPiPhiRestFrame = boostToPhiRestFrame(vecPi).Vect();
    auto momK1PhiRestFrame = boostToPhiRestFrame(vecK1).Vect();

    return momPiPhiRestFrame.Dot(momK1PhiRestFrame) / std::sqrt(momPiPhiRestFrame.Mag2() * momK1PhiRestFrame.Mag2());
  }

  template <typename T>
  auto cos3PiKDsToKKPi(const T& candidate)
  {
    auto cosPiK = cosPiKPhiRestFrame(candidate, 0);
    return cosPiK * cosPiK * cosPiK;
  }

  template <typename T>
  auto absCos3PiKDsToKKPi(const T& candidate)
  {
    return std::abs(cos3PiKDsToKKPi(candidate));
  }

  template <typename T>
  auto cos3PiKDsToPiKK(const T& candidate)
  {
    auto cosPiK = cosPiKPhiRestFrame(candidate, 1);
    return cosPiK * cosPiK * cosPiK;
  }

  template <typename T>
  auto absCos3PiKDsToPiKK(const T& candidate)
  {
    return std::abs(cos3PiKDsToPiKK(candidate));
  }

  // Λc± → p± K∓ π±

  template <typename T>
  auto ctLc(const T& candidate)
  {
    return candidate.ct(o2::constants::physics::MassLambdaCPlus);
  }

  template <typename T>
  auto yLc(const T& candidate)
  {
    return candidate.y(o2::constants::physics::MassLambdaCPlus);
  }

  template <typename T>
  auto eLc(const T& candidate)
  {
    return candidate.e(o2::constants::physics::MassLambdaCPlus);
  }

  template <typename T>
  auto invMassLcToPKPi(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassProton, o2::constants::physics::MassKPlus, o2::constants::physics::MassPiPlus});
  }

  template <typename T>
  auto invMassLcToPiKP(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassPiPlus, o2::constants::physics::MassKPlus, o2::constants::physics::MassProton});
  }

  template <typename T>
  auto invMassKPiPairLcToPKPi(const T& candidate)
  {
    return RecoDecay::m(std::array{candidate.pVectorProng1(), candidate.pVectorProng2()}, std::array{o2::constants::physics::MassKPlus, o2::constants::physics::MassPiPlus});
  }

  template <typename T>
  auto invMassKPiPairLcToPiKP(const T& candidate)
  {
    return RecoDecay::m(std::array{candidate.pVectorProng1(), candidate.pVectorProng0()}, std::array{o2::constants::physics::MassKPlus, o2::constants::physics::MassPiPlus});
  }

  template <typename T>
  auto invMassPKPairLcToPKPi(const T& candidate)
  {
    return RecoDecay::m(std::array{candidate.pVectorProng0(), candidate.pVectorProng1()}, std::array{o2::constants::physics::MassProton, o2::constants::physics::MassKPlus});
  }

  template <typename T>
  auto invMassPKPairLcToPiKP(const T& candidate)
  {
    return RecoDecay::m(std::array{candidate.pVectorProng2(), candidate.pVectorProng1()}, std::array{o2::constants::physics::MassProton, o2::constants::physics::MassKPlus});
  }

  template <typename T>
  auto invMassPPiPairLcToPKPi(const T& candidate)
  {
    return RecoDecay::m(std::array{candidate.pVectorProng0(), candidate.pVectorProng2()}, std::array{o2::constants::physics::MassProton, o2::constants::physics::MassPiPlus});
  }

  template <typename T>
  auto invMassPPiPairLcToPiKP(const T& candidate)
  {
    return RecoDecay::m(std::array{candidate.pVectorProng2(), candidate.pVectorProng0()}, std::array{o2::constants::physics::MassProton, o2::constants::physics::MassPiPlus});
  }

  // Ξc± → p± K∓ π±

  template <typename T>
  auto ctXic(const T& candidate)
  {
    return candidate.ct(o2::constants::physics::MassXiCPlus);
  }

  template <typename T>
  auto yXic(const T& candidate)
  {
    return candidate.y(o2::constants::physics::MassXiCPlus);
  }

  template <typename T>
  auto eXic(const T& candidate)
  {
    return candidate.e(o2::constants::physics::MassXiCPlus);
  }

  template <typename T>
  auto invMassXicToPKPi(const T& candidate)
  {
    return invMassLcToPKPi(candidate);
  }

  template <typename T>
  auto invMassXicToPiKP(const T& candidate)
  {
    return invMassLcToPiKP(candidate);
  }

  // hf_cand_casc_lf_2prong

  template <typename T>
  auto invMassXiczeroToXiPi(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassXiMinus, o2::constants::physics::MassPiPlus});
  }

  template <typename T>
  auto invMassOmegaczeroToOmegaPi(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassOmegaMinus, o2::constants::physics::MassPiPlus});
  }

  // hf_cand_casc_lf_3prong

  template <typename T>
  auto invMassXicplusToXiPiPi(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassXiMinus, o2::constants::physics::MassPiPlus, o2::constants::physics::MassPiPlus});
  }

  // hf_cand_x

  // X → Jpsi π+ π-
  template <typename T>
  auto ctX(const T& candidate)
  {
    return candidate.ct(o2::constants::physics::MassX3872);
  }

  template <typename T>
  auto yX(const T& candidate)
  {
    return candidate.y(o2::constants::physics::MassX3872);
  }

  template <typename T>
  auto eX(const T& candidate)
  {
    return candidate.e(o2::constants::physics::MassX3872);
  }

  template <typename T>
  auto invMassXToJpsiPiPi(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassJPsi, o2::constants::physics::MassPiPlus, o2::constants::physics::MassPiPlus});
  }

  /// Difference between the X mass and the sum of the J/psi and di-pion masses
  template <typename T>
  auto qX(const T& candidate)
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
  auto dRX(const T& candidate, int numPi)
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

    double deltaEta = etaJpsi - etaPi;
    double deltaPhi = RecoDecay::constrainAngle(phiJpsi - phiPi, -o2::constants::math::PI);

    return RecoDecay::sqrtSumOfSquares(deltaEta, deltaPhi);
  }

  /// Difference in pT between the two pions
  template <typename T>
  auto balancePtPionsX(const T& candidate)
  {
    double ptPi1 = RecoDecay::pt(candidate.pxProng1(), candidate.pyProng1());
    double ptPi2 = RecoDecay::pt(candidate.pxProng2(), candidate.pyProng2());
    return std::abs(ptPi1 - ptPi2) / (ptPi1 + ptPi2);
  }

  // Ξcc±± → p± K∓ π± π±

  template <typename T>
  auto ctXicc(const T& candidate)
  {
    return candidate.ct(o2::constants::physics::MassXiCCPlusPlus);
  }

  template <typename T>
  auto yXicc(const T& candidate)
  {
    return candidate.y(o2::constants::physics::MassXiCCPlusPlus);
  }

  template <typename T>
  auto eXicc(const T& candidate)
  {
    return candidate.e(o2::constants::physics::MassXiCCPlusPlus);
  }

  template <typename T>
  auto invMassXiccToXicPi(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassXiCPlus, o2::constants::physics::MassPiPlus});
  }

  // chic → Jpsi gamma

  template <typename T>
  auto ctChic(const T& candidate)
  {
    return candidate.ct(o2::constants::physics::MassChiC1);
  }

  template <typename T>
  auto yChic(const T& candidate)
  {
    return candidate.y(o2::constants::physics::MassChiC1);
  }

  template <typename T>
  auto eChic(const T& candidate)
  {
    return candidate.e(o2::constants::physics::MassChiC1);
  }
  template <typename T>
  auto invMassChicToJpsiGamma(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassJPsi, 0.});
  }

  // Λb → Λc+ π- → p K- π+ π-

  template <typename T>
  auto ctLb(const T& candidate)
  {
    return candidate.ct(o2::constants::physics::MassLambdaB0);
  }

  template <typename T>
  auto yLb(const T& candidate)
  {
    return candidate.y(o2::constants::physics::MassLambdaB0);
  }

  template <typename T>
  auto eLb(const T& candidate)
  {
    return candidate.e(o2::constants::physics::MassLambdaB0);
  }
  template <typename T>
  auto invMassLbToLcPi(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassLambdaCPlus, o2::constants::physics::MassPiPlus});
  }

  // B0(B0bar) → D∓ π±

  template <typename T>
  auto ctB0(const T& candidate)
  {
    return candidate.ct(o2::constants::physics::MassB0);
  }

  template <typename T>
  auto yB0(const T& candidate)
  {
    return candidate.y(o2::constants::physics::MassB0);
  }

  template <typename T>
  auto eB0(const T& candidate)
  {
    return candidate.e(o2::constants::physics::MassB0);
  }

  template <typename T>
  auto invMassB0ToDPi(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassDMinus, o2::constants::physics::MassPiPlus});
  }

  template <typename T>
  auto cosThetaStarB0(const T& candidate)
  {
    return candidate.cosThetaStar(std::array{o2::constants::physics::MassDMinus, o2::constants::physics::MassPiPlus}, o2::constants::physics::MassB0, 1);
  }

  // Bs(bar) → Ds∓ π±

  template <typename T>
  auto ctBs(const T& candidate)
  {
    return candidate.ct(o2::constants::physics::MassBS);
  }

  template <typename T>
  auto yBs(const T& candidate)
  {
    return candidate.y(o2::constants::physics::MassBS);
  }

  template <typename T>
  auto eBs(const T& candidate)
  {
    return candidate.e(o2::constants::physics::MassBS);
  }

  template <typename T>
  auto invMassBsToDsPi(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassDSBar, o2::constants::physics::MassPiPlus});
  }

  template <typename T>
  auto cosThetaStarBs(const T& candidate)
  {
    return candidate.cosThetaStar(std::array{o2::constants::physics::MassDSBar, o2::constants::physics::MassPiPlus}, o2::constants::physics::MassBS, 1);
  }

  /// Σc0,++ → Λc+(→pK-π+) π-,+

  /// @brief Sc inv. mass using reco mass for Lc in pKpi and PDG mass for pion
  template <typename T, typename U>
  auto invMassScRecoLcToPKPi(const T& candidateSc, const U& candidateLc)
  {
    return candidateSc.m(std::array{static_cast<double>(invMassLcToPKPi(candidateLc)), o2::constants::physics::MassPiPlus});
  }

  /// @brief Sc inv. mass using reco mass for Lc in piKp and PDG mass for pion
  template <typename T, typename U>
  auto invMassScRecoLcToPiKP(const T& candidateSc, const U& candidateLc)
  {
    return candidateSc.m(std::array{static_cast<double>(invMassLcToPiKP(candidateLc)), o2::constants::physics::MassPiPlus});
  }

  template <typename T>
  auto ySc0(const T& candidate)
  {
    return candidate.y(o2::constants::physics::MassSigmaC0);
  }

  template <typename T>
  auto yScPlusPlus(const T& candidate)
  {
    return candidate.y(o2::constants::physics::MassSigmaCPlusPlus);
  }

  /// Σc0,++ → Λc+(→K0sP) π-,+
  /// @brief Sc inv. mass using reco mass for Lc in K0sP and PDG mass for pion
  template <typename T, typename U>
  auto invMassScRecoLcToK0sP(const T& candidateSc, const U& candidateLc)
  {
    return candidateSc.m(std::array{static_cast<double>(invMassLcToK0sP(candidateLc)), o2::constants::physics::MassPiMinus});
  }

  /// Apply topological cuts as defined in SelectorCuts.h
  /// \param candB0 B0 candidate
  /// \param cuts B0 candidate selection per pT bin"
  /// \param binsPt pT bin limits
  /// \return true if candidate passes all selections
  template <typename T1, typename T2, typename T3>
  bool selectionB0ToDPiTopol(const T1& candB0, const T2& cuts, const T3& binsPt)
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
  template <typename T1 = int, typename T2 = bool>
  bool selectionB0ToDPiPid(const T1& pidTrackPi, const T2& acceptPIDNotApplicable)
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
  bool selectionBplusToD0PiTopol(const T1& candBp, const T2& cuts, const T3& binsPt)
  {
    auto ptcandBp = candBp.pt();
    auto ptPi = RecoDecay::pt(candBp.pxProng1(), candBp.pyProng1());

    int pTBin = o2::analysis::findBin(binsPt, ptcandBp);
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
  template <typename T1 = int, typename T2 = bool>
  bool selectionBplusToD0PiPid(const T1& pidTrackPi, const T2& acceptPIDNotApplicable)
  {
    if (!acceptPIDNotApplicable && pidTrackPi != TrackSelectorPID::Accepted) {
      return false;
    }
    if (acceptPIDNotApplicable && pidTrackPi == TrackSelectorPID::Rejected) {
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
  bool selectionBsToDsPiTopol(const T1& candBs, const T2& cuts, const T3& binsPt)
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
  template <typename T1 = int, typename T2 = bool>
  bool selectionBsToDsPiPid(const T1& pidTrackPi, const T2& acceptPIDNotApplicable)
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
  bool applySelectionDmesMlScoresForB(const T1& cuts, const T2& binsPtC, float ptC, std::vector<float> mlScores)
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
  bool selectionDmesMlScoresForB(const T1& candD, const T2& cuts, const T3& binsPtC, const std::vector<float>& mlScores)
  {
    return applySelectionDmesMlScoresForB(cuts, binsPtC, candD.pt(), mlScores);
  }

  /// Apply selection on ML scores for charm-hadron daughter in b-hadron decays in reduced format (common for all the beauty channels)
  /// \param candB b-hadron candidates
  /// \param cuts ML score selection per bin of charm-hadron pT
  /// \param binsPtC pT bin limits of charm hadron
  /// \return true if b-hadron candidate passes all selections
  template <typename T1, typename T2, typename T3>
  bool selectionDmesMlScoresForBReduced(const T1& candB, const T2& cuts, const T3& binsPtC)
  {
    std::vector<float> mlScores;
    mlScores.push_back(candB.prong0MlScoreBkg());
    mlScores.push_back(candB.prong0MlScorePrompt());
    mlScores.push_back(candB.prong0MlScoreNonprompt()); // we want non-prompt for beauty
    return applySelectionDmesMlScoresForB(cuts, binsPtC, RecoDecay::pt(candB.pxProng0(), candB.pyProng0()), mlScores);
  }

 private:
};

#endif // PWGHF_CORE_HFHELPER_H_
