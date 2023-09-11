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
///
/// \brief
/// \author Paul Buehler, paul.buehler@oeaw.ac.at
/// \since  01.10.2021

#ifndef PWGUD_CORE_UDHELPERS_H_
#define PWGUD_CORE_UDHELPERS_H_

#include <vector>
#include <bitset>
#include "TLorentzVector.h"
#include "Framework/Logger.h"
#include "DataFormatsFT0/Digit.h"
#include "CommonConstants/LHCConstants.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "PWGUD/Core/UPCHelpers.h"
#include "PWGUD/Core/DGCutparHolder.h"

using namespace o2;
using namespace o2::framework;

// namespace with helpers for UD framework
namespace udhelpers
{

// -----------------------------------------------------------------------------
// return net charge of PV tracks
template <bool onlyPV, typename std::enable_if<onlyPV>::type* = nullptr, typename TCs>
int8_t netCharge(TCs tracks)
{
  int8_t nch = 0;
  for (auto track : tracks) {
    if (track.isPVContributor()) {
      nch += track.sign();
    }
  }
  return nch;
}

// -----------------------------------------------------------------------------
// return net charge of tracks
template <bool onlyPV, typename std::enable_if<!onlyPV>::type* = nullptr, typename TCs>
int8_t netCharge(TCs tracks)
{
  int8_t nch = 0;
  for (auto track : tracks) {
    nch += track.sign();
  }
  return nch;
}

// -----------------------------------------------------------------------------
// return fraction of PV tracks with a TOF hit
template <bool onlyPV, typename std::enable_if<onlyPV>::type* = nullptr, typename TCs>
float rPVtrwTOF(TCs tracks, int nPVTracks)
{
  float rpvrwTOF = 0.;
  for (auto& track : tracks) {
    if (track.isPVContributor() && track.hasTOF()) {
      rpvrwTOF += 1.;
    }
  }
  if (nPVTracks > 0) {
    rpvrwTOF /= nPVTracks;
  }
  return rpvrwTOF;
}

// return fraction of tracks with a TOF hit
template <bool onlyPV, typename std::enable_if<!onlyPV>::type* = nullptr, typename TCs>
float rPVtrwTOF(TCs tracks, int nPVTracks)
{
  float rpvrwTOF = 0.;
  for (auto& track : tracks) {
    if (track.hasTOF()) {
      rpvrwTOF += 1.;
    }
  }
  if (nPVTracks > 0) {
    rpvrwTOF /= nPVTracks;
  }
  return rpvrwTOF;
}

// -----------------------------------------------------------------------------
// The associations between collsisions and BCs can be ambiguous.
// By default a collision is associated with the BC closest in time.
// Any BC falling within a BC window of meanBC +- deltaBC could potentially be the
// true BC.
//
template <typename T>
T compatibleBCs(uint64_t meanBC, int deltaBC, T const& bcs);

template <typename I, typename T>
T compatibleBCs(I& bcIter, uint64_t meanBC, int deltaBC, T const& bcs);

// In this variant of compatibleBCs the bcIter is ideally placed within
// [minBC, maxBC], but it does not need to be. The range is given by meanBC +- delatBC.
template <typename I, typename T>
T compatibleBCs(I& bcIter, uint64_t meanBC, int deltaBC, T const& bcs)
{
  // range of BCs to consider
  uint64_t minBC = (uint64_t)deltaBC < meanBC ? meanBC - (uint64_t)deltaBC : 0;
  uint64_t maxBC = meanBC + (uint64_t)deltaBC;
  LOGF(debug, "  minBC %d maxBC %d bcIterator %d (%d)", minBC, maxBC, bcIter.globalBC(), bcIter.globalIndex());

  // check [min,max]BC to overlap with [bcs.iteratorAt([0,bcs.size() - 1])
  if (maxBC < bcs.iteratorAt(0).globalBC() || minBC > bcs.iteratorAt(bcs.size() - 1).globalBC()) {
    LOGF(info, "<compatibleBCs> No overlap of [%d, %d] and [%d, %d]", minBC, maxBC, bcs.iteratorAt(0).globalBC(), bcs.iteratorAt(bcs.size() - 1).globalBC());
    return T{{bcs.asArrowTable()->Slice(0, 0)}, (uint64_t)0};
  }

  // find slice of BCs table with BC in [minBC, maxBC]
  int moveCount = 0;
  int64_t minBCId = bcIter.globalIndex();
  int64_t maxBCId = bcIter.globalIndex();

  // lower limit
  if (bcIter.globalBC() < minBC) {
    while (bcIter != bcs.end() && bcIter.globalBC() < minBC) {
      ++bcIter;
      ++moveCount;
      minBCId = bcIter.globalIndex();
    }
  } else {
    while (bcIter.globalIndex() > 0 && bcIter.globalBC() >= minBC) {
      minBCId = bcIter.globalIndex();
      --bcIter;
      --moveCount;
    }
  }

  // upper limit limit
  if (bcIter.globalBC() < maxBC) {
    while (bcIter != bcs.end() && bcIter.globalBC() <= maxBC) {
      maxBCId = bcIter.globalIndex();
      ++bcIter;
      ++moveCount;
    }

  } else {
    while (bcIter.globalIndex() > 0 && bcIter.globalBC() > maxBC) {
      --bcIter;
      --moveCount;
      maxBCId = bcIter.globalIndex();
    }
  }
  LOGF(debug, "  BC range: %d - %d", minBCId, maxBCId);

  // reset bcIter
  bcIter.moveByIndex(-moveCount);

  // create bc slice
  T slice{{bcs.asArrowTable()->Slice(minBCId, maxBCId - minBCId + 1)}, (uint64_t)minBCId};
  bcs.copyIndexBindings(slice);
  LOGF(debug, "  size of slice %d", slice.size());
  return slice;
}

// In this variant of compatibleBCs the range of compatible BCs is calculated from the
// collision time and the time resolution dt. Typically the range is +- 4*dt.
template <typename C, typename T>
T compatibleBCs(C const& collision, int ndt, T const& bcs, int nMinBCs = 7)
{
  LOGF(debug, "Collision time / resolution [ns]: %f / %f", collision.collisionTime(), collision.collisionTimeRes());

  // return if collisions has no associated BC
  if (!collision.has_foundBC() || ndt < 0) {
    return T{{bcs.asArrowTable()->Slice(0, 0)}, (uint64_t)0};
  }

  // get associated BC
  auto bcIter = collision.template foundBC_as<T>();

  // due to the filling scheme the most probable BC may not be the one estimated from the collision time
  uint64_t mostProbableBC = bcIter.globalBC();
  uint64_t meanBC = mostProbableBC + std::lround(collision.collisionTime() / o2::constants::lhc::LHCBunchSpacingNS);

  // enforce minimum number for deltaBC
  int deltaBC = std::ceil(collision.collisionTimeRes() / o2::constants::lhc::LHCBunchSpacingNS * ndt);
  if (deltaBC < nMinBCs) {
    deltaBC = nMinBCs;
  }
  LOGF(debug, "BC %d,  deltaBC %d", bcIter.globalIndex(), deltaBC);

  return compatibleBCs(bcIter, meanBC, deltaBC, bcs);
}

// In this variant of compatibleBCs the range of compatible BCs is defined by meanBC +- deltaBC.
template <typename T>
T compatibleBCs(uint64_t meanBC, int deltaBC, T const& bcs)
{
  // find BC with globalBC ~ meanBC
  uint64_t ind = (uint64_t)(bcs.size() / 2);
  auto bcIter = bcs.iteratorAt(ind);

  return compatibleBCs(bcIter, meanBC, deltaBC, bcs);
}

// -----------------------------------------------------------------------------
// Same as above but for collisions with MC information
template <typename F, typename T>
T MCcompatibleBCs(F const& collision, int ndt, T const& bcs, int nMinBCs = 7)
{
  LOGF(debug, "Collision time / resolution [ns]: %f / %f", collision.collisionTime(), collision.collisionTimeRes());

  // return if collisions has no associated BC
  if (!collision.has_foundBC()) {
    LOGF(info, "Collision %i - no BC found!", collision.globalIndex());
    return T{{bcs.asArrowTable()->Slice(0, 0)}, (uint64_t)0};
  }

  // get associated BC
  auto bcIter = collision.template foundBC_as<T>();

  // due to the filling scheme the most probable BC may not be the one estimated from the collision time
  uint64_t mostProbableBC = bcIter.globalBC();
  uint64_t meanBC = mostProbableBC + std::lround(collision.collisionTime() / o2::constants::lhc::LHCBunchSpacingNS);

  // enforce minimum number for deltaBC
  int deltaBC = std::ceil(collision.collisionTimeRes() / o2::constants::lhc::LHCBunchSpacingNS * ndt);
  if (deltaBC < nMinBCs) {
    deltaBC = nMinBCs;
  }

  return compatibleBCs(bcIter, meanBC, deltaBC, bcs);
}

// -----------------------------------------------------------------------------
// function to check if track provides good PID information
// Checks the nSigma for any particle assumption to be within limits.
template <typename TC>
bool hasGoodPID(DGCutparHolder diffCuts, TC track)
{
  // El, Mu, Pi, Ka, and Pr are considered
  // at least one nSigma must be within set limits
  LOGF(debug, "TPC PID %f / %f / %f / %f / %f",
       track.tpcNSigmaEl(),
       track.tpcNSigmaMu(),
       track.tpcNSigmaPi(),
       track.tpcNSigmaKa(),
       track.tpcNSigmaPr());
  if (TMath::Abs(track.tpcNSigmaEl()) < diffCuts.maxNSigmaTPC()) {
    return true;
  }
  if (TMath::Abs(track.tpcNSigmaMu()) < diffCuts.maxNSigmaTPC()) {
    return true;
  }
  if (TMath::Abs(track.tpcNSigmaPi()) < diffCuts.maxNSigmaTPC()) {
    return true;
  }
  if (TMath::Abs(track.tpcNSigmaKa()) < diffCuts.maxNSigmaTPC()) {
    return true;
  }
  if (TMath::Abs(track.tpcNSigmaPr()) < diffCuts.maxNSigmaTPC()) {
    return true;
  }

  if (track.hasTOF()) {
    LOGF(debug, "TOF PID %f / %f / %f / %f / %f",
         track.tofNSigmaEl(),
         track.tofNSigmaMu(),
         track.tofNSigmaPi(),
         track.tofNSigmaKa(),
         track.tofNSigmaPr());
    if (TMath::Abs(track.tofNSigmaEl()) < diffCuts.maxNSigmaTOF()) {
      return true;
    }
    if (TMath::Abs(track.tofNSigmaMu()) < diffCuts.maxNSigmaTOF()) {
      return true;
    }
    if (TMath::Abs(track.tofNSigmaPi()) < diffCuts.maxNSigmaTOF()) {
      return true;
    }
    if (TMath::Abs(track.tofNSigmaKa()) < diffCuts.maxNSigmaTOF()) {
      return true;
    }
    if (TMath::Abs(track.tofNSigmaPr()) < diffCuts.maxNSigmaTOF()) {
      return true;
    }
  }
  return false;
}

// -----------------------------------------------------------------------------
float FV0AmplitudeA(aod::FV0A&& fv0)
{
  float totAmplitude = 0;
  for (auto amp : fv0.amplitude()) {
    totAmplitude += amp;
  }

  return totAmplitude;
}

// -----------------------------------------------------------------------------
float FT0AmplitudeA(aod::FT0&& ft0)
{
  float totAmplitude = 0;
  for (auto amp : ft0.amplitudeA()) {
    totAmplitude += amp;
  }

  return totAmplitude;
}

// -----------------------------------------------------------------------------
float FT0AmplitudeC(aod::FT0&& ft0)
{
  float totAmplitude = 0;
  for (auto amp : ft0.amplitudeC()) {
    totAmplitude += amp;
  }

  return totAmplitude;
}

// -----------------------------------------------------------------------------
int16_t FDDAmplitudeA(aod::FDD&& fdd)
{
  int16_t totAmplitude = 0;
  for (auto amp : fdd.chargeA()) {
    totAmplitude += amp;
  }

  return totAmplitude;
}

// -----------------------------------------------------------------------------
int16_t FDDAmplitudeC(aod::FDD&& fdd)
{
  int16_t totAmplitude = 0;
  for (auto amp : fdd.chargeC()) {
    totAmplitude += amp;
  }

  return totAmplitude;
}

// -----------------------------------------------------------------------------
template <typename T>
bool cleanFV0(T& bc, float maxFITtime, float limitA)
{
  if (bc.has_foundFV0()) {
    bool ota = std::abs(bc.foundFV0().time()) <= maxFITtime;
    bool oma = FV0AmplitudeA(bc.foundFV0()) <= limitA;
    return ota && oma;
  } else {
    return true;
  }
}

// -----------------------------------------------------------------------------
template <typename T>
bool cleanFT0(T& bc, float maxFITtime, float limitA, float limitC)
{
  if (bc.has_foundFT0()) {

    // check times and amplitudes
    bool ota = std::abs(bc.foundFT0().timeA()) <= maxFITtime;
    bool otc = std::abs(bc.foundFT0().timeC()) <= maxFITtime;
    bool oma = FT0AmplitudeA(bc.foundFT0()) <= limitA;
    bool omc = FT0AmplitudeC(bc.foundFT0()) <= limitC;

    // compare decisions with FT0 trigger decisions
    std::bitset<8> triggers = bc.foundFT0().triggerMask();
    bool ora = !triggers[o2::ft0::Triggers::bitA];
    bool orc = !triggers[o2::ft0::Triggers::bitC];
    LOGF(debug, "ota %f otc %f ora/FT0AmplitudeA %d/%d orc/FT0AmplitudeC %d/%d", bc.foundFT0().timeA(), bc.foundFT0().timeC(), ora, oma, orc, omc);

    return ota && oma && otc && omc;
  } else {
    return true;
  }
}

// -----------------------------------------------------------------------------
template <typename T>
bool cleanFDD(T& bc, float maxFITtime, float limitA, float limitC)
{
  if (bc.has_foundFDD()) {
    bool ota = std::abs(bc.foundFDD().timeA()) <= maxFITtime;
    bool otc = std::abs(bc.foundFDD().timeC()) <= maxFITtime;
    bool oma = FDDAmplitudeA(bc.foundFDD()) <= limitA;
    bool omc = FDDAmplitudeC(bc.foundFDD()) <= limitC;
    return ota && oma && otc && omc;
  } else {
    return true;
  }
}

// -----------------------------------------------------------------------------
// FIT amplitude limits
//  lims[0]: FV0A
//  lims[1]: FT0A
//  lims[2]: FT0C
//  lims[3]: FDDA
//  lims[4]: FDDC

template <typename T>
bool cleanFIT(T& bc, float maxFITtime, std::vector<float> lims)
{
  return cleanFV0(bc, maxFITtime, lims[0]) &&
         cleanFT0(bc, maxFITtime, lims[1], lims[2]) &&
         cleanFDD(bc, maxFITtime, lims[3], lims[4]);
}
template <typename T>
bool cleanFITCollision(T& col, float maxFITtime, std::vector<float> lims)
{
  bool isCleanFV0 = true;
  if (col.has_foundFV0()) {
    isCleanFV0 = (std::abs(col.foundFV0().time()) <= maxFITtime) && (FV0AmplitudeA(col.foundFV0()) < lims[0]);
  }
  bool isCleanFT0 = true;
  if (col.has_foundFT0()) {
    isCleanFT0 = (std::abs(col.foundFT0().timeA()) <= maxFITtime) && (FT0AmplitudeA(col.foundFT0()) < lims[1]) &&
                 (std::abs(col.foundFT0().timeC()) <= maxFITtime) && (FT0AmplitudeC(col.foundFT0()) < lims[2]);
  }
  bool isCleanFDD = true;
  if (col.has_foundFDD()) {
    isCleanFDD = (std::abs(col.foundFDD().timeA()) <= maxFITtime) && (FDDAmplitudeA(col.foundFDD()) < lims[3]) &&
                 (std::abs(col.foundFDD().timeC()) <= maxFITtime) && (FDDAmplitudeC(col.foundFDD()) < lims[4]);
  }
  return (isCleanFV0 && isCleanFT0 && isCleanFDD);
}
// -----------------------------------------------------------------------------
// fill BB and BG information into FITInfo
template <typename BCR>
void fillBGBBFlags(upchelpers::FITInfo& info, uint64_t const& minbc, BCR const& bcrange)
{
  for (auto const& bc2u : bcrange) {

    // 0 <= bit <= 31
    auto bit = bc2u.globalBC() - minbc;
    if (!bc2u.selection_bit(o2::aod::evsel::kNoBGT0A))
      SETBIT(info.BGFT0Apf, bit);
    if (!bc2u.selection_bit(o2::aod::evsel::kNoBGT0C))
      SETBIT(info.BGFT0Cpf, bit);
    if (bc2u.selection_bit(o2::aod::evsel::kIsBBT0A))
      SETBIT(info.BBFT0Apf, bit);
    if (bc2u.selection_bit(o2::aod::evsel::kIsBBT0C))
      SETBIT(info.BBFT0Cpf, bit);
    if (!bc2u.selection_bit(o2::aod::evsel::kNoBGV0A))
      SETBIT(info.BGFV0Apf, bit);
    if (bc2u.selection_bit(o2::aod::evsel::kIsBBV0A))
      SETBIT(info.BBFV0Apf, bit);
    if (!bc2u.selection_bit(o2::aod::evsel::kNoBGFDA))
      SETBIT(info.BGFDDApf, bit);
    if (!bc2u.selection_bit(o2::aod::evsel::kNoBGFDC))
      SETBIT(info.BGFDDCpf, bit);
    if (bc2u.selection_bit(o2::aod::evsel::kIsBBFDA))
      SETBIT(info.BBFDDApf, bit);
    if (bc2u.selection_bit(o2::aod::evsel::kIsBBFDC))
      SETBIT(info.BBFDDCpf, bit);
  }
}

// -----------------------------------------------------------------------------
// extract FIT information
template <typename B>
void getFITinfo(upchelpers::FITInfo info, uint64_t const& bcnum, B const& bcs, aod::FT0s const& ft0s, aod::FV0As const& fv0as, aod::FDDs const& fdds)
{
  // find bc with globalBC = bcnum
  Partition<B> selbc = aod::bc::globalBC == bcnum;
  selbc.bindTable(bcs);

  // if BC exists then update FIT information for this BC
  if (selbc.size() > 0) {
    auto bc = selbc.begin();

    // FT0
    if (bc.has_foundFT0()) {
      auto ft0 = ft0s.iteratorAt(bc.foundFT0Id());
      info.timeFT0A = ft0.timeA();
      info.timeFT0C = ft0.timeC();
      const auto& ampsA = ft0.amplitudeA();
      const auto& ampsC = ft0.amplitudeC();
      info.ampFT0A = 0.;
      for (auto amp : ampsA) {
        info.ampFT0A += amp;
      }
      info.ampFT0C = 0.;
      for (auto amp : ampsC) {
        info.ampFT0C += amp;
      }
      info.triggerMaskFT0 = ft0.triggerMask();
    }

    // FV0A
    if (bc.has_foundFV0()) {
      auto fv0a = fv0as.iteratorAt(bc.foundFV0Id());
      info.timeFV0A = fv0a.time();
      const auto& amps = fv0a.amplitude();
      info.ampFV0A = 0.;
      for (auto amp : amps) {
        info.ampFV0A += amp;
      }
      info.triggerMaskFV0A = fv0a.triggerMask();
    }

    // FDD
    if (bc.has_foundFDD()) {
      auto fdd = fdds.iteratorAt(bc.foundFDDId());
      info.timeFDDA = fdd.timeA();
      info.timeFDDC = fdd.timeC();
      const auto& ampsA = fdd.chargeA();
      const auto& ampsC = fdd.chargeC();
      info.ampFDDA = 0.;
      for (auto amp : ampsA) {
        info.ampFDDA += amp;
      }
      info.ampFDDC = 0.;
      for (auto amp : ampsC) {
        info.ampFDDC += amp;
      }
      info.triggerMaskFDD = fdd.triggerMask();
    }
  }

  // fill BG and BB flags
  auto minbc = bcnum - 16;
  auto maxbc = bcnum + 15;
  Partition<B> bcrange = aod::bc::globalBC >= minbc && aod::bc::globalBC <= maxbc;
  bcrange.bindTable(bcs);
  fillBGBBFlags(info, minbc, bcrange);
}

// -----------------------------------------------------------------------------
template <typename T>
bool cleanZDC(T const& bc, aod::Zdcs& zdcs, std::vector<float>& lims, SliceCache& cache)
{
  const auto& ZdcBC = zdcs.sliceByCached(aod::zdc::bcId, bc.globalIndex(), cache);
  return (ZdcBC.size() == 0);
}

// -----------------------------------------------------------------------------
template <typename T>
bool cleanCalo(T const& bc, aod::Calos& calos, std::vector<float>& lims, SliceCache& cache)
{
  const auto& CaloBC = calos.sliceByCached(aod::calo::bcId, bc.globalIndex(), cache);
  return (CaloBC.size() == 0);
}

// -----------------------------------------------------------------------------
// check if all tracks come from same MCCollision
template <typename T>
int64_t sameMCCollision(T tracks, aod::McCollisions mccols, aod::McParticles mcparts)
{
  int64_t colID = -1;
  for (auto const& track : tracks) {
    if (track.has_mcParticle()) {
      auto mcpart = track.mcParticle();
      if (mcpart.has_mcCollision()) {
        auto mccol = mcpart.mcCollision();
        if (colID < 0) {
          colID = mccol.globalIndex();
        } else {
          if (colID != mccol.globalIndex()) {
            return (int64_t)-1;
          }
        }
      } else {
        return (int64_t)-1;
      }
    } else {
      return (int64_t)-1;
    }
  }

  return colID;
}

// -----------------------------------------------------------------------------
// In PYTHIA a central diffractive produced (CD) particle has the ID
// 9900110. Check the particles of a MC event to contain a CD particle.
template <typename T>
bool isPythiaCDE(T MCparts)
{
  for (auto mcpart : MCparts) {
    if (mcpart.pdgCode() == 9900110) {
      return true;
    }
  }
  return false;
}

// -----------------------------------------------------------------------------
// In pp events produced with GRANIITTI the stack starts with
// 22212/22212/99/22212/2212/99/90
template <typename T>
bool isGraniittiCDE(T MCparts)
{

  for (auto MCpart : MCparts) {
    LOGF(debug, " MCpart.pdgCode() %d", MCpart.pdgCode());
  }
  LOGF(debug, "");

  if (MCparts.size() < 7) {
    return false;
  } else {
    if (MCparts.iteratorAt(0).pdgCode() != 2212)
      return false;
    if (MCparts.iteratorAt(1).pdgCode() != 2212)
      return false;
    if (MCparts.iteratorAt(2).pdgCode() != 99)
      return false;
    if (MCparts.iteratorAt(3).pdgCode() != 2212)
      return false;
    if (MCparts.iteratorAt(4).pdgCode() != 2212)
      return false;
    if (MCparts.iteratorAt(5).pdgCode() != 99)
      return false;
    if (MCparts.iteratorAt(6).pdgCode() != 90)
      return false;
  }

  return true;
}

// -----------------------------------------------------------------------------
// Invariant mass of GRANIITTI generated event
template <typename T>
TLorentzVector ivmGraniittiCDE(T MCparts)
{
  TLorentzVector ivm = TLorentzVector(0., 0., 0., 0.);

  // is this a GRANIITTI generated event?
  if (isGraniittiCDE(MCparts)) {
    TLorentzVector lvtmp;

    for (int ii = 7; ii < MCparts.size(); ii++) {
      auto mcPart = MCparts.iteratorAt(ii);
      LOGF(debug, "  part %d / %d", mcPart.pdgCode(), mcPart.getGenStatusCode());
      if (mcPart.getGenStatusCode() == 0) {
        lvtmp.SetXYZT(mcPart.px(), mcPart.py(), mcPart.pz(), mcPart.e());
        ivm += lvtmp;
      }
    }
    LOGF(debug, "");
  }

  return ivm;
}

// -----------------------------------------------------------------------------

} // namespace udhelpers

#endif // PWGUD_CORE_UDHELPERS_H_
