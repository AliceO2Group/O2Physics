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

/// \file HFFilterTables.h
/// \brief Header file with definition of variables, methods, and tables used in the HFFilter.cxx task
///
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN
/// \author Marcel Lesch <marcel.lesch@tum.de>, TUM
/// \author Alexandre Bigot <alexandre.bigot@cern.ch>, Strasbourg University
/// \author Biao Zhang <biao.zhang@cern.ch>, CCNU

#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Boost.h"

#include "Framework/AnalysisDataModel.h"
#include "Common/Core/RecoDecay.h"

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

enum beautyTrackSelection {
  kRejected = 0,
  kSoftPion,
  kRegular
};

static const std::array<std::string, kNtriggersHF> HfTriggerNames{"highPt", "beauty", "femto", "doubleCharm"};
static const std::array<std::string, kNCharmParticles> charmParticleNames{"D0", "Dplus", "Ds", "Lc", "Xic"};
static const std::array<std::string, kNBeautyParticles> beautyParticleNames{"Bplus", "B0toDStar", "B0", "Bs", "Lb", "Xib"};
static const std::array<int, kNCharmParticles> pdgCodesCharm{421, 411, 431, 4122, 4232};

static const std::tuple pdgCharmDaughters{
  std::array{-321, 211},        // D0
  std::array{-321, 211, 211},   // Dplus
  std::array{321, -321, 211},   // Ds
  std::array{2212, -321, 211},  // Lc
  std::array{2212, -321, 211}}; // Xic

static const float massPi = RecoDecay::getMassPDG(211);
static const float massK = RecoDecay::getMassPDG(321);
static const float massProton = RecoDecay::getMassPDG(2212);
static const float massPhi = RecoDecay::getMassPDG(333);
static const float massD0 = RecoDecay::getMassPDG(421);
static const float massDPlus = RecoDecay::getMassPDG(411);
static const float massDs = RecoDecay::getMassPDG(431);
static const float massLc = RecoDecay::getMassPDG(4122);
static const float massXic = RecoDecay::getMassPDG(4232);
static const float massDStar = RecoDecay::getMassPDG(413);
static const float massBPlus = RecoDecay::getMassPDG(511);
static const float massB0 = RecoDecay::getMassPDG(521);
static const float massBs = RecoDecay::getMassPDG(531);
static const float massLb = RecoDecay::getMassPDG(5122);
static const float massXib = RecoDecay::getMassPDG(5232);

/// Computation of the relative momentum between particle pairs
/// \param track is a track
/// \param ProtonMass is the mass of a proton
/// \param CharmCandMomentum is the three momentum of a charm candidate
/// \param CharmMass is the mass of the charm hadron
/// \return relative momentum of pair
template <typename T> // template <typename T, typename C>
float computeRelativeMomentum(const T& track, const std::array<float, 3>& CharmCandMomentum, const float& CharmMass)
{
  ROOT::Math::PxPyPzMVector part1(track.px(), track.py(), track.pz(), massProton);
  ROOT::Math::PxPyPzMVector part2(CharmCandMomentum[0], CharmCandMomentum[1], CharmCandMomentum[2], CharmMass);

  ROOT::Math::PxPyPzMVector trackSum = part1 + part2;
  ROOT::Math::Boost boostv12{trackSum.BoostToCM()};
  ROOT::Math::PxPyPzMVector part1CM = boostv12(part1);
  ROOT::Math::PxPyPzMVector part2CM = boostv12(part2);
  ROOT::Math::PxPyPzMVector trackRelK = part1CM - part2CM;

  float kStar = 0.5 * trackRelK.P();
  return kStar;
} // float computeRelativeMomentum(const T& track, const std::array<float, 3>& CharmCandMomentum, const float& CharmMass)

} // namespace hffilters

// definition of tables

namespace extra2Prong
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
} // namespace extra2Prong
namespace extra3Prong
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
} // namespace extra3Prong
DECLARE_SOA_TABLE(Colls2Prong, "AOD", "COLLSID2P", o2::aod::extra2Prong::CollisionId);
DECLARE_SOA_TABLE(Colls3Prong, "AOD", "COLLSID3P", o2::aod::extra3Prong::CollisionId);

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
                  hftraining::FlagOrigin);
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
                  hftraining::HFSelBit);

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
} // namespace o2::aod