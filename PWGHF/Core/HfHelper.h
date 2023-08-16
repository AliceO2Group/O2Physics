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

#include <Math/GenVector/Boost.h>
#include <Math/Vector4D.h>
#include <TPDGCode.h>

#include "Framework/O2DatabasePDGPlugin.h"

#include "Common/Core/RecoDecay.h"

#include "PWGHF/Core/SelectorCuts.h"

class HfHelper
{
 public:
  /// Default constructor
  HfHelper() = default;

  /// Default destructor
  ~HfHelper() = default;

  o2::framework::Service<o2::framework::O2DatabasePDG> pdg;

  auto mass(int pdgCode)
  {
    return pdg->Mass(pdgCode);
  }

  // 2-prong

  // D0(bar) → π± K∓

  template <typename T>
  auto ctD0(const T& candidate)
  {
    return candidate.ct(mass(o2::analysis::pdg::kD0));
  }

  template <typename T>
  auto yD0(const T& candidate)
  {
    return candidate.y(mass(o2::analysis::pdg::kD0));
  }

  template <typename T>
  auto eD0(const T& candidate)
  {
    return candidate.e(mass(o2::analysis::pdg::kD0));
  }

  template <typename T>
  auto invMassD0ToPiK(const T& candidate)
  {
    return candidate.m(std::array{mass(kPiPlus), mass(kKPlus)});
  }

  template <typename T>
  auto invMassD0barToKPi(const T& candidate)
  {
    return candidate.m(std::array{mass(kKPlus), mass(kPiPlus)});
  }

  template <typename T>
  auto cosThetaStarD0(const T& candidate)
  {
    return candidate.cosThetaStar(std::array{mass(kPiPlus), mass(kKPlus)}, mass(o2::analysis::pdg::kD0), 1);
  }

  template <typename T>
  auto cosThetaStarD0bar(const T& candidate)
  {
    return candidate.cosThetaStar(std::array{mass(kKPlus), mass(kPiPlus)}, mass(o2::analysis::pdg::kD0), 0);
  }

  // J/ψ

  template <typename T>
  auto ctJpsi(const T& candidate)
  {
    return candidate.ct(mass(o2::analysis::pdg::kJPsi));
  }

  template <typename T>
  auto yJpsi(const T& candidate)
  {
    return candidate.y(mass(o2::analysis::pdg::kJPsi));
  }

  template <typename T>
  auto eJpsi(const T& candidate)
  {
    return candidate.e(mass(o2::analysis::pdg::kJPsi));
  }

  // J/ψ → e+ e−
  template <typename T>
  auto invMassJpsiToEE(const T& candidate)
  {
    return candidate.m(std::array{mass(kElectron), mass(kElectron)});
  }
  // J/ψ → μ+ μ−

  template <typename T>
  auto invMassJpsiToMuMu(const T& candidate)
  {
    return candidate.m(std::array{mass(kMuonPlus), mass(kMuonMinus)});
  }

  // hf_cand_casc

  template <typename T>
  auto invMassLcToK0sP(const T& candidate)
  {
    return candidate.m(std::array{mass(kProton), mass(kK0Short)}); // first daughter is bachelor
  }

  template <typename T>
  auto invMassGammaToEE(const T& candidate)
  {
    return candidate.m(std::array{mass(kElectron), mass(kElectron)});
  }

  template <typename T>
  auto ctV0K0s(const T& candidate)
  {
    return candidate.ctV0(mass(kK0Short));
  }

  template <typename T>
  auto ctV0Lambda(const T& candidate)
  {
    return candidate.ctV0(mass(kLambda0));
  }

  // B± → D0bar(D0) π±

  template <typename T>
  auto ctBplus(const T& candidate)
  {
    return candidate.ct(mass(o2::analysis::pdg::kBPlus));
  }

  template <typename T>
  auto yBplus(const T& candidate)
  {
    return candidate.y(mass(o2::analysis::pdg::kBPlus));
  }

  template <typename T>
  auto eBplus(const T& candidate)
  {
    return candidate.e(mass(o2::analysis::pdg::kBPlus));
  }

  template <typename T>
  auto invMassBplusToD0Pi(const T& candidate)
  {
    return candidate.m(std::array{mass(o2::analysis::pdg::kD0), mass(kPiPlus)});
  }

  template <typename T>
  auto cosThetaStarBplus(const T& candidate)
  {
    return candidate.cosThetaStar(std::array{mass(o2::analysis::pdg::kD0), mass(kPiPlus)}, mass(o2::analysis::pdg::kBPlus), 1);
  }

  // 3-prong

  // D± → π± K∓ π±

  template <typename T>
  auto ctDplus(const T& candidate)
  {
    return candidate.ct(mass(o2::analysis::pdg::kDPlus));
  }

  template <typename T>
  auto yDplus(const T& candidate)
  {
    return candidate.y(mass(o2::analysis::pdg::kDPlus));
  }

  template <typename T>
  auto eDplus(const T& candidate)
  {
    return candidate.e(mass(o2::analysis::pdg::kDPlus));
  }

  template <typename T>
  auto invMassDplusToPiKPi(const T& candidate)
  {
    return candidate.m(std::array{mass(kPiPlus), mass(kKPlus), mass(kPiPlus)});
  }

  // Ds± → K± K∓ π±

  template <typename T>
  auto ctDs(const T& candidate)
  {
    return candidate.ct(mass(o2::analysis::pdg::kDS));
  }

  template <typename T>
  auto yDs(const T& candidate)
  {
    return candidate.y(mass(o2::analysis::pdg::kDS));
  }

  template <typename T>
  auto eDs(const T& candidate)
  {
    return candidate.e(mass(o2::analysis::pdg::kDS));
  }

  template <typename T>
  auto invMassDsToKKPi(const T& candidate)
  {
    return candidate.m(std::array{mass(kKPlus), mass(kKPlus), mass(kPiPlus)});
  }

  template <typename T>
  auto invMassDsToPiKK(const T& candidate)
  {
    return candidate.m(std::array{mass(kPiPlus), mass(kKPlus), mass(kKPlus)});
  }

  template <typename T>
  auto deltaMassPhiDsToKKPi(const T& candidate)
  {
    double invMassKKpair = RecoDecay::m(std::array{candidate.pVectorProng0(), candidate.pVectorProng1()}, std::array{mass(kKPlus), mass(kKPlus)});
    return std::abs(invMassKKpair - mass(o2::analysis::pdg::kPhi));
  }

  template <typename T>
  auto deltaMassPhiDsToPiKK(const T& candidate)
  {
    double invMassKKpair = RecoDecay::m(std::array{candidate.pVectorProng1(), candidate.pVectorProng2()}, std::array{mass(kKPlus), mass(kKPlus)});
    return std::abs(invMassKKpair - mass(o2::analysis::pdg::kPhi));
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

    ROOT::Math::PxPyPzMVector vecPi(momPi[0], momPi[1], momPi[2], mass(kPiPlus));
    ROOT::Math::PxPyPzMVector vecK1(momK1[0], momK1[1], momK1[2], mass(kKPlus));
    ROOT::Math::PxPyPzMVector vecK2(momK2[0], momK2[1], momK2[2], mass(kKPlus));
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
  auto cos3PiKDsToPiKK(const T& candidate)
  {
    auto cosPiK = cosPiKPhiRestFrame(candidate, 1);
    return cosPiK * cosPiK * cosPiK;
  }

  // Λc± → p± K∓ π±

  template <typename T>
  auto ctLc(const T& candidate)
  {
    return candidate.ct(mass(o2::analysis::pdg::kLambdaCPlus));
  }

  template <typename T>
  auto yLc(const T& candidate)
  {
    return candidate.y(mass(o2::analysis::pdg::kLambdaCPlus));
  }

  template <typename T>
  auto eLc(const T& candidate)
  {
    return candidate.e(mass(o2::analysis::pdg::kLambdaCPlus));
  }

  template <typename T>
  auto invMassLcToPKPi(const T& candidate)
  {
    return candidate.m(std::array{mass(kProton), mass(kKPlus), mass(kPiPlus)});
  }

  template <typename T>
  auto invMassLcToPiKP(const T& candidate)
  {
    return candidate.m(std::array{mass(kPiPlus), mass(kKPlus), mass(kProton)});
  }

  // Ξc± → p± K∓ π±

  template <typename T>
  auto ctXic(const T& candidate)
  {
    return candidate.ct(mass(o2::analysis::pdg::kXiCPlus));
  }

  template <typename T>
  auto yXic(const T& candidate)
  {
    return candidate.y(mass(o2::analysis::pdg::kXiCPlus));
  }

  template <typename T>
  auto eXic(const T& candidate)
  {
    return candidate.e(mass(o2::analysis::pdg::kXiCPlus));
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
    return candidate.m(std::array{mass(kXiMinus), mass(kPiPlus)});
  }

  template <typename T>
  auto invMassOmegaczeroToOmegaPi(const T& candidate)
  {
    return candidate.m(std::array{mass(kOmegaMinus), mass(kPiPlus)});
  }

  // hf_cand_casc_lf_3prong

  template <typename T>
  auto invMassXicplusToXiPiPi(const T& candidate)
  {
    return candidate.m(std::array{mass(kXiMinus), mass(kPiPlus), mass(kPiPlus)});
  }

  // hf_cand_x

  // X → Jpsi π+ π-
  // TODO: add pdg code for X (9920443), temporarily hardcode mass here:
  float massX = 3.872; // replace this with: "mass(9920443)" when pdg is added
  template <typename T>
  auto ctX(const T& candidate)
  {
    return candidate.ct(massX);
  }

  template <typename T>
  auto yX(const T& candidate)
  {
    return candidate.y(massX);
  }

  template <typename T>
  auto eX(const T& candidate)
  {
    return candidate.e(massX);
  }

  template <typename T>
  auto invMassXToJpsiPiPi(const T& candidate)
  {
    return candidate.m(std::array{mass(443), mass(kPiPlus), mass(kPiPlus)});
  }

  /// Difference between the X mass and the sum of the J/psi and di-pion masses
  template <typename T>
  auto qX(const T& candidate)
  {
    auto piVec1 = std::array{candidate.pxProng1(), candidate.pyProng1(), candidate.pzProng1()};
    auto piVec2 = std::array{candidate.pxProng2(), candidate.pyProng2(), candidate.pzProng2()};
    double massPi = mass(kPiPlus);

    auto arrayMomenta = std::array{piVec1, piVec2};
    double massPiPi = RecoDecay::m(arrayMomenta, std::array{massPi, massPi});

    // PDG mass, as reported in CMS paper https://arxiv.org/pdf/1302.3968.pdf
    double massJpsi = mass(o2::analysis::pdg::kJPsi);

    double massX = invMassXToJpsiPiPi(candidate);
    return std::abs(massX - massJpsi - massPiPi);
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
    return candidate.ct(mass(o2::analysis::pdg::kXiCCPlusPlus));
  }

  template <typename T>
  auto yXicc(const T& candidate)
  {
    return candidate.y(mass(o2::analysis::pdg::kXiCCPlusPlus));
  }

  template <typename T>
  auto eXicc(const T& candidate)
  {
    return candidate.e(mass(o2::analysis::pdg::kXiCCPlusPlus));
  }

  template <typename T>
  auto invMassXiccToXicPi(const T& candidate)
  {
    return candidate.m(std::array{mass(o2::analysis::pdg::kXiCPlus), mass(kPiPlus)});
  }

  // chic → Jpsi gamma

  template <typename T>
  auto ctChic(const T& candidate)
  {
    return candidate.ct(mass(o2::analysis::pdg::kChiC1));
  }

  template <typename T>
  auto yChic(const T& candidate)
  {
    return candidate.y(mass(o2::analysis::pdg::kChiC1));
  }

  template <typename T>
  auto eChic(const T& candidate)
  {
    return candidate.e(mass(o2::analysis::pdg::kChiC1));
  }
  template <typename T>
  auto invMassChicToJpsiGamma(const T& candidate)
  {
    return candidate.m(std::array{mass(o2::analysis::pdg::kJPsi), 0.});
  }

  // Λb → Λc+ π- → p K- π+ π-

  template <typename T>
  auto ctLb(const T& candidate)
  {
    return candidate.ct(mass(o2::analysis::pdg::kLambdaB0));
  }

  template <typename T>
  auto yLb(const T& candidate)
  {
    return candidate.y(mass(o2::analysis::pdg::kLambdaB0));
  }

  template <typename T>
  auto eLb(const T& candidate)
  {
    return candidate.e(mass(o2::analysis::pdg::kLambdaB0));
  }
  template <typename T>
  auto invMassLbToLcPi(const T& candidate)
  {
    return candidate.m(std::array{mass(o2::analysis::pdg::kLambdaCPlus), mass(kPiPlus)});
  }

  // B0(B0bar) → D∓ π±

  template <typename T>
  auto ctB0(const T& candidate)
  {
    return candidate.ct(mass(o2::analysis::pdg::kB0));
  }

  template <typename T>
  auto yB0(const T& candidate)
  {
    return candidate.y(mass(o2::analysis::pdg::kB0));
  }

  template <typename T>
  auto eB0(const T& candidate)
  {
    return candidate.e(mass(o2::analysis::pdg::kB0));
  }

  template <typename T>
  auto invMassB0ToDPi(const T& candidate)
  {
    return candidate.m(std::array{mass(o2::analysis::pdg::kDMinus), mass(kPiPlus)});
  }

  template <typename T>
  auto cosThetaStarB0(const T& candidate)
  {
    return candidate.cosThetaStar(std::array{mass(o2::analysis::pdg::kDMinus), mass(kPiPlus)}, mass(o2::analysis::pdg::kB0), 1);
  }

  // Bs(bar) → Ds∓ π±

  template <typename T>
  auto ctBs(const T& candidate)
  {
    return candidate.ct(mass(o2::analysis::pdg::kBS));
  }

  template <typename T>
  auto yBs(const T& candidate)
  {
    return candidate.y(mass(o2::analysis::pdg::kBS));
  }

  template <typename T>
  auto eBs(const T& candidate)
  {
    return candidate.e(mass(o2::analysis::pdg::kBS));
  }

  template <typename T>
  auto invMassBsToDsPi(const T& candidate)
  {
    return candidate.m(std::array{mass(o2::analysis::pdg::kDSBar), mass(kPiPlus)});
  }

  template <typename T>
  auto cosThetaStarBs(const T& candidate)
  {
    return candidate.cosThetaStar(std::array{mass(o2::analysis::pdg::kDSBar), mass(kPiPlus)}, mass(o2::analysis::pdg::kBS), 1);
  }

  /// Σc0,++ → Λc+(→pK-π+) π-,+

  /// @brief Sc inv. mass using reco mass for Lc in pKpi and PDG mass for pion
  template <typename T, typename U>
  auto invMassScRecoLcToPKPi(const T& candidateSc, const U& candidateLc)
  {
    return candidateSc.m(std::array{static_cast<double>(invMassLcToPKPi(candidateLc)), mass(kPiPlus)});
  }

  /// @brief Sc inv. mass using reco mass for Lc in piKp and PDG mass for pion
  template <typename T, typename U>
  auto invMassScRecoLcToPiKP(const T& candidateSc, const U& candidateLc)
  {
    return candidateSc.m(std::array{static_cast<double>(invMassLcToPiKP(candidateLc)), mass(kPiPlus)});
  }

  template <typename T>
  auto ySc0(const T& candidate)
  {
    return candidate.y(mass(o2::analysis::pdg::kSigmaC0));
  }

  template <typename T>
  auto yScPlusPlus(const T& candidate)
  {
    return candidate.y(mass(o2::analysis::pdg::kSigmaCPlusPlus));
  }

 private:
};

#endif // PWGHF_CORE_HFHELPER_H_
