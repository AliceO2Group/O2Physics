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
#include "DataFormatsFIT/Triggers.h"
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
  for (const auto& track : tracks) {
    if (track.isPVContributor()) {
      nch += track.sign();
    }
  }
  return nch;
}

// return net charge of tracks
template <bool onlyPV, typename std::enable_if<!onlyPV>::type* = nullptr, typename TCs>
int8_t netCharge(TCs tracks)
{
  int8_t nch = 0;
  for (const auto& track : tracks) {
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
  for (const auto& track : tracks) {
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
  for (const auto& track : tracks) {
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
T compatibleBCs(uint64_t const& meanBC, int const& deltaBC, T const& bcs);

template <typename B, typename T>
T compatibleBCs(B const& bc, uint64_t const& meanBC, int const& deltaBC, T const& bcs);

// In this variant of compatibleBCs the bcIter is ideally placed within
// [minBC, maxBC], but it does not need to be. The range is given by meanBC +- delatBC.
template <typename B, typename T>
T compatibleBCs(B const& bc, uint64_t const& meanBC, int const& deltaBC, T const& bcs)
{
  // get BCs iterator
  auto bcIter = bcs.iteratorAt(bc.globalIndex());

  // range of BCs to consider
  uint64_t minBC = static_cast<uint64_t>(deltaBC) < meanBC ? meanBC - static_cast<uint64_t>(deltaBC) : 0;
  uint64_t maxBC = meanBC + static_cast<uint64_t>(deltaBC);
  LOGF(debug, "  minBC %d maxBC %d bcIterator %d (%d) #BCs %d", minBC, maxBC, bcIter.globalBC(), bcIter.globalIndex(), bcs.size());

  // check [min,max]BC to overlap with [bcs.iteratorAt([0,bcs.size() - 1])
  if (maxBC < bcs.iteratorAt(0).globalBC() || minBC > bcs.iteratorAt(bcs.size() - 1).globalBC()) {
    LOGF(debug, "<compatibleBCs> No overlap of [%d, %d] and [%d, %d]", minBC, maxBC, bcs.iteratorAt(0).globalBC(), bcs.iteratorAt(bcs.size() - 1).globalBC());
    return T{{bcs.asArrowTable()->Slice(0, 0)}, static_cast<uint64_t>(0)};
  }

  // find slice of BCs table with BC in [minBC, maxBC]
  int64_t minBCId = bcIter.globalIndex();
  int64_t maxBCId = bcIter.globalIndex();

  // lower limit
  if (bcIter.globalBC() < minBC) {
    while (bcIter != bcs.end() && bcIter.globalBC() < minBC) {
      ++bcIter;
      minBCId = bcIter.globalIndex();
    }
  } else {
    while (bcIter.globalIndex() > 0 && bcIter.globalBC() >= minBC) {
      minBCId = bcIter.globalIndex();
      --bcIter;
    }
  }

  // upper limit limit
  if (bcIter.globalBC() < maxBC) {
    while (bcIter != bcs.end() && bcIter.globalBC() <= maxBC) {
      maxBCId = bcIter.globalIndex();
      ++bcIter;
    }
  } else {
    while (bcIter.globalIndex() > 0 && bcIter.globalBC() > maxBC) {
      --bcIter;
      maxBCId = bcIter.globalIndex();
    }
  }

  // create bc slice
  T bcslice{{bcs.asArrowTable()->Slice(minBCId, maxBCId - minBCId + 1)}, static_cast<uint64_t>(minBCId)};
  bcs.copyIndexBindings(bcslice);
  LOGF(debug, "  size of slice %d", bcslice.size());
  return bcslice;
}

// In this variant of compatibleBCs the range of compatible BCs is calculated from the
// collision time and the time resolution dt. Typically the range is +- 4*dt.
template <typename C, typename T>
T compatibleBCs(C const& collision, int ndt, T const& bcs, int nMinBCs = 7)
{
  LOGF(debug, "Collision time / resolution [ns]: %f / %f", collision.collisionTime(), collision.collisionTimeRes());

  // return if collisions has no associated BC
  if (!collision.has_foundBC() || ndt < 0) {
    return T{{bcs.asArrowTable()->Slice(0, 0)}, static_cast<uint64_t>(0)};
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
T compatibleBCs(uint64_t const& meanBC, int const& deltaBC, T const& bcs)
{
  // find BC with globalBC ~ meanBC
  uint64_t ind = static_cast<uint64_t>(bcs.size() / 2);
  auto bcIter = bcs.iteratorAt(ind);

  return compatibleBCs(bcIter, meanBC, deltaBC, bcs);
}

// -----------------------------------------------------------------------------
// Same as above but for collisions with MC information
template <typename F, typename T>
T MCcompatibleBCs(F const& collision, int const& ndt, T const& bcs, int const& nMinBCs = 7)
{
  LOGF(debug, "Collision time / resolution [ns]: %f / %f", collision.collisionTime(), collision.collisionTimeRes());

  // return if collisions has no associated BC
  if (!collision.has_foundBC()) {
    LOGF(debug, "Collision %i - no BC found!", collision.globalIndex());
    return T{{bcs.asArrowTable()->Slice(0, 0)}, static_cast<uint64_t>(0)};
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
  if (std::abs(track.tpcNSigmaEl()) < diffCuts.maxNSigmaTPC()) {
    return true;
  }
  if (std::abs(track.tpcNSigmaMu()) < diffCuts.maxNSigmaTPC()) {
    return true;
  }
  if (std::abs(track.tpcNSigmaPi()) < diffCuts.maxNSigmaTPC()) {
    return true;
  }
  if (std::abs(track.tpcNSigmaKa()) < diffCuts.maxNSigmaTPC()) {
    return true;
  }
  if (std::abs(track.tpcNSigmaPr()) < diffCuts.maxNSigmaTPC()) {
    return true;
  }

  if (track.hasTOF()) {
    LOGF(debug, "TOF PID %f / %f / %f / %f / %f",
         track.tofNSigmaEl(),
         track.tofNSigmaMu(),
         track.tofNSigmaPi(),
         track.tofNSigmaKa(),
         track.tofNSigmaPr());
    if (std::abs(track.tofNSigmaEl()) < diffCuts.maxNSigmaTOF()) {
      return true;
    }
    if (std::abs(track.tofNSigmaMu()) < diffCuts.maxNSigmaTOF()) {
      return true;
    }
    if (std::abs(track.tofNSigmaPi()) < diffCuts.maxNSigmaTOF()) {
      return true;
    }
    if (std::abs(track.tofNSigmaKa()) < diffCuts.maxNSigmaTOF()) {
      return true;
    }
    if (std::abs(track.tofNSigmaPr()) < diffCuts.maxNSigmaTOF()) {
      return true;
    }
  }
  return false;
}

// -----------------------------------------------------------------------------
template <typename TFV0>
float FV0AmplitudeA(TFV0 fv0)
{
  const auto& ampsA = fv0.amplitude();
  return std::accumulate(ampsA.begin(), ampsA.end(), 0.f);
}

// -----------------------------------------------------------------------------
template <typename TFT0>
float FT0AmplitudeA(TFT0 ft0)
{
  const auto& ampsA = ft0.amplitudeA();
  return std::accumulate(ampsA.begin(), ampsA.end(), 0.f);
}

// -----------------------------------------------------------------------------
template <typename TFT0>
float FT0AmplitudeC(TFT0 ft0)
{
  const auto& ampsC = ft0.amplitudeC();
  return std::accumulate(ampsC.begin(), ampsC.end(), 0.f);
}

// -----------------------------------------------------------------------------
template <typename TFDD>
float FDDAmplitudeA(TFDD fdd)
{
  std::vector<int16_t> ampsA(fdd.chargeA(), fdd.chargeA() + 8);
  return std::accumulate(ampsA.begin(), ampsA.end(), 0);
}

// -----------------------------------------------------------------------------
template <typename TFDD>
float FDDAmplitudeC(TFDD fdd)
{
  std::vector<int16_t> ampsC(fdd.chargeC(), fdd.chargeC() + 8);
  return std::accumulate(ampsC.begin(), ampsC.end(), 0);
}

// -----------------------------------------------------------------------------
template <typename T>
bool cleanFV0(T& bc, float maxFITtime, float limitA)
{
  if (bc.has_foundFV0() && limitA >= 0.) {
    bool ota = std::abs(bc.foundFV0().time()) <= maxFITtime;
    bool oma = FV0AmplitudeA(bc.foundFV0()) <= limitA;
    return ota && oma;
  } else {
    return true;
  }
}

// -----------------------------------------------------------------------------
template <typename T>
bool cleanFT0A(T& bc, float maxFITtime, float limitA)
{
  if (bc.has_foundFT0() && limitA >= 0.) {
    bool ota = std::abs(bc.foundFT0().timeA()) <= maxFITtime;
    bool oma = FT0AmplitudeA(bc.foundFT0()) <= limitA;
    return ota && oma;
  } else {
    return true;
  }
}

// -----------------------------------------------------------------------------
template <typename T>
bool cleanFT0C(T& bc, float maxFITtime, float limitC)
{
  if (bc.has_foundFT0() && limitC >= 0.) {
    bool otc = std::abs(bc.foundFT0().timeC()) <= maxFITtime;
    bool omc = FT0AmplitudeC(bc.foundFT0()) <= limitC;
    return otc && omc;
  } else {
    return true;
  }
}

// -----------------------------------------------------------------------------
template <typename T>
bool cleanFDDA(T& bc, float maxFITtime, float limitA)
{
  if (bc.has_foundFDD() && limitA >= 0.) {
    bool ota = std::abs(bc.foundFDD().timeA()) <= maxFITtime;
    bool oma = FDDAmplitudeA(bc.foundFDD()) <= limitA;
    return ota && oma;
  } else {
    return true;
  }
}

// -----------------------------------------------------------------------------
template <typename T>
bool cleanFDDC(T& bc, float maxFITtime, float limitC)
{
  if (bc.has_foundFDD() && limitC >= 0.) {
    bool otc = std::abs(bc.foundFDD().timeC()) <= maxFITtime;
    bool omc = FDDAmplitudeC(bc.foundFDD()) <= limitC;
    return otc && omc;
  } else {
    return true;
  }
}

// -----------------------------------------------------------------------------
template <typename T>
bool cleanFT0(T& bc, float maxFITtime, float limitA, float limitC)
{
  return cleanFT0A(bc, maxFITtime, limitA) &&
         cleanFT0C(bc, maxFITtime, limitC);
}

// -----------------------------------------------------------------------------
template <typename T>
bool cleanFDD(T& bc, float maxFITtime, float limitA, float limitC)
{
  return cleanFDDA(bc, maxFITtime, limitA) &&
         cleanFDDC(bc, maxFITtime, limitC);
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
    isCleanFV0 = lims[0] < 0. ? true : (std::abs(col.foundFV0().time()) <= maxFITtime) && (FV0AmplitudeA(col.foundFV0()) < lims[0]);
  }
  bool isCleanFT0 = true;
  if (col.has_foundFT0()) {
    isCleanFT0 = (lims[1] < 0. ? true : (std::abs(col.foundFT0().timeA()) <= maxFITtime) && (FT0AmplitudeA(col.foundFT0()) < lims[1])) &&
                 (lims[2] < 0. ? true : (std::abs(col.foundFT0().timeC()) <= maxFITtime) && (FT0AmplitudeC(col.foundFT0()) < lims[2]));
  }
  bool isCleanFDD = true;
  if (col.has_foundFDD()) {
    isCleanFDD = (lims[3] < 0. ? true : (std::abs(col.foundFDD().timeA()) <= maxFITtime) && (FDDAmplitudeA(col.foundFDD()) < lims[3])) &&
                 (lims[4] < 0. ? true : (std::abs(col.foundFDD().timeC()) <= maxFITtime) && (FDDAmplitudeC(col.foundFDD()) < lims[4]));
  }
  return (isCleanFV0 && isCleanFT0 && isCleanFDD);
}

// -----------------------------------------------------------------------------
template <typename T>
bool cleanFITA(T& bc, float maxFITtime, std::vector<float> lims)
{
  return cleanFV0(bc, maxFITtime, lims[0]) &&
         cleanFT0A(bc, maxFITtime, lims[1]) &&
         cleanFDDA(bc, maxFITtime, lims[3]);
}

// -----------------------------------------------------------------------------
template <typename T>
bool cleanFITC(T& bc, float maxFITtime, std::vector<float> lims)
{
  return cleanFT0C(bc, maxFITtime, lims[2]) &&
         cleanFDDC(bc, maxFITtime, lims[4]);
}

// -----------------------------------------------------------------------------
template <typename T>
bool TVX(T& bc)
{
  bool tvx = false;
  if (bc.has_foundFT0()) {
    auto ft0 = bc.foundFT0();
    tvx = TESTBIT(ft0.triggerMask(), o2::fit::Triggers::bitVertex);
  }
  return tvx;
}

// -----------------------------------------------------------------------------
template <typename T>
bool TSC(T& bc)
{
  bool tsc = false;
  if (bc.has_foundFT0()) {
    auto ft0 = bc.foundFT0();
    tsc = TESTBIT(ft0.triggerMask(), o2::fit::Triggers::bitSCen);
  }
  return tsc;
}

// -----------------------------------------------------------------------------
template <typename T>
bool TCE(T& bc)
{
  bool tce = false;
  if (bc.has_foundFT0()) {
    auto ft0 = bc.foundFT0();
    tce = TESTBIT(ft0.triggerMask(), o2::fit::Triggers::bitCen);
  }
  return tce;
}

// -----------------------------------------------------------------------------
template <typename T>
bool TOR(T& bc, float maxFITtime, std::vector<float> lims)
{
  auto torA = !cleanFT0A(bc, maxFITtime, lims[1]);
  auto torC = !cleanFT0C(bc, maxFITtime, lims[2]);
  return torA || torC;
}

// -----------------------------------------------------------------------------
// returns true if veto is active
// return false if veto is not active
template <typename T>
bool FITveto(T const& bc, DGCutparHolder const& diffCuts)
{
  // return if FIT veto is found in bc
  // Double Gap (DG) condition
  // 4 types of vetoes:
  //  0 TVX
  //  1 TSC
  //  2 TCE
  //  3 TOR
  if (diffCuts.withTVX()) {
    return TVX(bc);
  }
  if (diffCuts.withTSC()) {
    return TSC(bc);
  }
  if (diffCuts.withTCE()) {
    return TCE(bc);
  }
  if (diffCuts.withTOR()) {
    return !cleanFIT(bc, diffCuts.maxFITtime(), diffCuts.FITAmpLimits());
  }
  return false;
}

// -----------------------------------------------------------------------------

template <typename T>
bool cutNoTimeFrameBorder(T const& coll)
// Reject collisions close to TF borders due to incomplete TPC drift volume.
// https://its.cern.ch/jira/browse/O2-4623
// Return true when event is good.
{
  return coll.selection_bit(o2::aod::evsel::kNoTimeFrameBorder);
}

// -----------------------------------------------------------------------------

template <typename T>
bool cutNoSameBunchPileup(T const& coll)
// Rejects collisions which are associated with the same "found-by-T0" bunch crossing.
// Could be partially due to the pileup with another collision in the same foundBC.
// See more in slides 12-14 of https://indico.cern.ch/event/1396220/#1-event-selection-with-its-rof.
// Return true when event is good.
{
  return coll.selection_bit(o2::aod::evsel::kNoSameBunchPileup);
}

// -----------------------------------------------------------------------------

template <typename T>
bool cutNoITSROFrameBorder(T const& coll)
// Reject events affected by the ITS ROF border.
// https://its.cern.ch/jira/browse/O2-4309
// Return true when event is good.
{
  return coll.selection_bit(o2::aod::evsel::kNoITSROFrameBorder);
}

// -----------------------------------------------------------------------------

template <typename T>
bool cutIsGoodZvtxFT0vsPV(T const& coll)
// Removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference.
// The large vertexZ difference can be due to the in-bunch pileup or wrong BC assigned to a collision.
// Return true when event is good.
{
  return coll.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV);
}

// -----------------------------------------------------------------------------

template <typename T>
bool cutIsVertexITSTPC(T const& coll)
// Selects collisions with at least one ITS-TPC PV track, and thus rejects vertices built from ITS-only tracks.
// Has an effect only on the pp data, in Pb-Pb ITS-only vertices are already rejected by default.
// Return true when event is good.
{
  return coll.selection_bit(o2::aod::evsel::kIsVertexITSTPC);
}

// -----------------------------------------------------------------------------

template <typename T>
bool cutIsVertexTRDmatched(T const& coll)
// Selects collisions with at least one TRD PV track.
// Return true when event is good.
{
  return coll.selection_bit(o2::aod::evsel::kIsVertexTRDmatched);
}

// -----------------------------------------------------------------------------

template <typename T>
bool cutIsVertexTOFmatched(T const& coll)
// Selects collisions with at least one TOF PV track.
// Return true when event is good.
{
  return coll.selection_bit(o2::aod::evsel::kIsVertexTOFmatched);
}

// -----------------------------------------------------------------------------

template <typename T>
bool goodCollision(T const& coll, DGCutparHolder const& diffCuts)
// Return true if collision is accepted according to user-chosen rules from event selection task
{
  bool accepted = true;
  std::vector<int> sels = diffCuts.collisionSel();
  if (sels[0])
    accepted = accepted && cutNoTimeFrameBorder(coll);
  if (sels[1])
    accepted = accepted && cutNoSameBunchPileup(coll);
  if (sels[2])
    accepted = accepted && cutNoITSROFrameBorder(coll);
  if (sels[3])
    accepted = accepted && cutIsGoodZvtxFT0vsPV(coll);
  if (sels[4])
    accepted = accepted && cutIsVertexITSTPC(coll);
  if (sels[5])
    accepted = accepted && cutIsVertexTRDmatched(coll);
  if (sels[6])
    accepted = accepted && cutIsVertexTOFmatched(coll);

  return accepted;
}

// -----------------------------------------------------------------------------
// fill BB and BG information into FITInfo
template <typename BCR>
void fillBGBBFlags(upchelpers::FITInfo& info, uint64_t const& minbc, BCR const& bcrange)
{
  for (auto const& bc : bcrange) {

    // 0 <= bit <= 31
    auto bit = bc.globalBC() - minbc;
    if (bit < 0 || bit > 31)
      continue;

    if (!bc.selection_bit(o2::aod::evsel::kNoBGT0A))
      SETBIT(info.BGFT0Apf, bit);
    if (!bc.selection_bit(o2::aod::evsel::kNoBGT0C))
      SETBIT(info.BGFT0Cpf, bit);
    if (bc.selection_bit(o2::aod::evsel::kIsBBT0A))
      SETBIT(info.BBFT0Apf, bit);
    if (bc.selection_bit(o2::aod::evsel::kIsBBT0C))
      SETBIT(info.BBFT0Cpf, bit);
    if (!bc.selection_bit(o2::aod::evsel::kNoBGV0A))
      SETBIT(info.BGFV0Apf, bit);
    if (bc.selection_bit(o2::aod::evsel::kIsBBV0A))
      SETBIT(info.BBFV0Apf, bit);
    if (!bc.selection_bit(o2::aod::evsel::kNoBGFDA))
      SETBIT(info.BGFDDApf, bit);
    if (!bc.selection_bit(o2::aod::evsel::kNoBGFDC))
      SETBIT(info.BGFDDCpf, bit);
    if (bc.selection_bit(o2::aod::evsel::kIsBBFDA))
      SETBIT(info.BBFDDApf, bit);
    if (bc.selection_bit(o2::aod::evsel::kIsBBFDC))
      SETBIT(info.BBFDDCpf, bit);
  }
}

// -----------------------------------------------------------------------------
// extract FIT information
template <typename BC, typename BCS>
void getFITinfo(upchelpers::FITInfo& info, BC& bc, BCS const& bcs, aod::FT0s const& ft0s, aod::FV0As const& fv0as, aod::FDDs const& fdds)
{
  // FV0A
  if (bc.has_foundFV0()) {
    auto fv0 = fv0as.iteratorAt(bc.foundFV0Id());
    info.timeFV0A = fv0.time();
    info.ampFV0A = FV0AmplitudeA(fv0);
    info.triggerMaskFV0A = fv0.triggerMask();
  }

  // FT0
  if (bc.has_foundFT0()) {
    auto ft0 = ft0s.iteratorAt(bc.foundFT0Id());
    info.timeFT0A = ft0.timeA();
    info.timeFT0C = ft0.timeC();
    info.ampFT0A = FT0AmplitudeA(ft0);
    info.ampFT0C = FT0AmplitudeC(ft0);
    info.triggerMaskFT0 = ft0.triggerMask();
  }

  // FDD
  if (bc.has_foundFDD()) {
    auto fdd = fdds.iteratorAt(bc.foundFDDId());
    info.timeFDDA = fdd.timeA();
    info.timeFDDC = fdd.timeC();
    info.ampFDDA = FDDAmplitudeA(fdd);
    info.ampFDDC = FDDAmplitudeC(fdd);
    info.triggerMaskFDD = fdd.triggerMask();
  }

  // fill BG and BB flags
  auto bcnum = bc.globalBC();
  auto bcrange = compatibleBCs(bc, bcnum, 16, bcs);
  LOGF(debug, "size of bcrange %d", bcrange.size());
  fillBGBBFlags(info, bcnum - 16, bcrange);
}

// -----------------------------------------------------------------------------
template <typename T>
bool cleanZDC(T const& bc, aod::Zdcs& zdcs, std::vector<float>& /*lims*/, SliceCache& cache)
{
  const auto& ZdcBC = zdcs.sliceByCached(aod::zdc::bcId, bc.globalIndex(), cache);
  return (ZdcBC.size() == 0);
}

// -----------------------------------------------------------------------------
template <typename T>
bool cleanCalo(T const& bc, aod::Calos& calos, std::vector<float>& /*lims*/, SliceCache& cache)
{
  const auto& CaloBC = calos.sliceByCached(aod::calo::bcId, bc.globalIndex(), cache);
  return (CaloBC.size() == 0);
}

// -----------------------------------------------------------------------------
// check if all tracks come from same MCCollision
template <typename T>
int64_t sameMCCollision(T tracks, aod::McCollisions, aod::McParticles)
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
            return static_cast<int64_t>(-1);
          }
        }
      } else {
        return static_cast<int64_t>(-1);
      }
    } else {
      return static_cast<int64_t>(-1);
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
  for (const auto& mcpart : MCparts) {
    if (mcpart.pdgCode() == 9900110) {
      return true;
    }
  }
  return false;
}

// -----------------------------------------------------------------------------
// In J/Psi -> mu+ + mu- events generated with STARlight the stack starts with
// 443013, 13, -13 or 443013, -13, 13
template <typename T>
bool isSTARLightJPsimumu(T MCparts)
{
  if (MCparts.size() < 3) {
    return false;
  } else {
    if (MCparts.iteratorAt(0).pdgCode() != 443013)
      return false;
    if (std::abs(MCparts.iteratorAt(1).pdgCode()) != 13)
      return false;
    if (MCparts.iteratorAt(2).pdgCode() != -MCparts.iteratorAt(1).pdgCode())
      return false;
  }
  return true;
}

// -----------------------------------------------------------------------------
// PbPb di electron production
// [15, 11, 13], [15, 11, 13]
template <typename T>
bool isUpcgen(T MCparts)
{
  if (MCparts.size() < 4) {
    return false;
  } else {
    auto pid1 = std::abs(MCparts.iteratorAt(0).pdgCode());
    auto pid2 = std::abs(MCparts.iteratorAt(1).pdgCode());
    if (pid1 != 11 && pid1 != 13 && pid1 != 15)
      return false;
    if (pid2 != 11 && pid2 != 13 && pid2 != 15)
      return false;
  }
  return true;
}

// -----------------------------------------------------------------------------
// In pp events produced with GRANIITTI the stack starts with
// 22212/22212/22212/2212/[211,321,]/[211,321,]
template <typename T>
bool isGraniittiCDE(T MCparts)
{

  // for (auto MCpart : MCparts) {
  //   LOGF(info, " MCpart.pdgCode() %d", MCpart.pdgCode());
  // }
  // LOGF(debug, "");

  if (MCparts.size() < 6) {
    return false;
  } else {
    if (MCparts.iteratorAt(0).pdgCode() != 2212)
      return false;
    if (MCparts.iteratorAt(1).pdgCode() != 2212)
      return false;
    if (MCparts.iteratorAt(2).pdgCode() != 2212)
      return false;
    if (MCparts.iteratorAt(3).pdgCode() != 2212)
      return false;
  }

  return true;
}

// -----------------------------------------------------------------------------
// function to select MC events of interest
template <typename T>
int isOfInterest(T MCparts)
{

  // PYTHIA CDE
  if (isPythiaCDE(MCparts)) {
    return 1;
  }

  // GRANIITTI CDE
  if (isGraniittiCDE(MCparts)) {
    return 2;
  }

  // STARLIGHT J/Psi -> mu+ + mu-
  if (isSTARLightJPsimumu(MCparts)) {
    return 3;
  }

  // Upcgen
  if (isUpcgen(MCparts)) {
    return 4;
  }

  return 0;
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
