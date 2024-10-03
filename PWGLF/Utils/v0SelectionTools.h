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

#ifndef PWGLF_UTILS_V0SELECTIONTOOLS_H_
#define PWGLF_UTILS_V0SELECTIONTOOLS_H_

#include "v0SelectionBits.h"
#include "v0SelectionGroup.h"

// simple checkers, but ensure 64 bit integers
#define bitset(var, nbit) ((var) |= (static_cast<uint64_t>(1) << static_cast<uint64_t>(nbit)))
#define bitcheck(var, nbit) ((var) & (static_cast<uint64_t>(1) << static_cast<uint64_t>(nbit)))

namespace v0data
{

// utility method to calculate a selection map for the V0s
template <typename TV0, typename TTrack, typename TCollision>
uint64_t computeReconstructionBitmap(TV0 v0, TTrack posTrackExtra, TTrack negTrackExtra, TCollision collision, float rapidityLambda, float rapidityK0Short, const v0SelectionGroup& v0sels)
// precalculate this information so that a check is one mask operation, not many
{
  uint64_t bitMap = 0;
  // Base topological variables
  if (v0.v0radius() > v0sels.getv0radius())
    bitset(bitMap, v0data::selRadius);
  if (v0.v0radius() < v0sels.getv0radiusMax())
    bitset(bitMap, v0data::selRadiusMax);
  if (TMath::Abs(v0.dcapostopv()) > v0sels.getdcapostopv())
    bitset(bitMap, v0data::selDCAPosToPV);
  if (TMath::Abs(v0.dcanegtopv()) > v0sels.getdcanegtopv())
    bitset(bitMap, v0data::selDCANegToPV);
  if (v0.v0cosPA() > v0sels.getv0cospa())
    bitset(bitMap, v0data::selCosPA);
  if (v0.dcaV0daughters() < v0sels.getdcav0dau())
    bitset(bitMap, v0data::selDCAV0Dau);

  // rapidity
  if (TMath::Abs(rapidityLambda) < v0sels.getRapidityCut())
    bitset(bitMap, v0data::selLambdaRapidity);
  if (TMath::Abs(rapidityK0Short) < v0sels.getRapidityCut())
    bitset(bitMap, v0data::selK0ShortRapidity);

  // ITS quality flags
  if (posTrackExtra.itsNCls() >= v0sels.getminITSclusters())
    bitset(bitMap, v0data::selPosGoodITSTrack);
  if (negTrackExtra.itsNCls() >= v0sels.getminITSclusters())
    bitset(bitMap, v0data::selNegGoodITSTrack);

  // TPC quality flags
  if (posTrackExtra.tpcCrossedRows() >= v0sels.getminTPCrows())
    bitset(bitMap, v0data::selPosGoodTPCTrack);
  if (negTrackExtra.tpcCrossedRows() >= v0sels.getminTPCrows())
    bitset(bitMap, v0data::selNegGoodTPCTrack);

  // TPC PID
  if (fabs(posTrackExtra.tpcNSigmaPi()) < v0sels.getTpcPidNsigmaCut())
    bitset(bitMap, v0data::selTPCPIDPositivePion);
  if (fabs(posTrackExtra.tpcNSigmaPr()) < v0sels.getTpcPidNsigmaCut())
    bitset(bitMap, v0data::selTPCPIDPositiveProton);
  if (fabs(negTrackExtra.tpcNSigmaPi()) < v0sels.getTpcPidNsigmaCut())
    bitset(bitMap, v0data::selTPCPIDNegativePion);
  if (fabs(negTrackExtra.tpcNSigmaPr()) < v0sels.getTpcPidNsigmaCut())
    bitset(bitMap, v0data::selTPCPIDNegativeProton);

  // TOF PID in DeltaT (deprecated, kept for compatibility)
  // Positive track
  if (fabs(v0.posTOFDeltaTLaPr()) < v0sels.getmaxDeltaTimeProton())
    bitset(bitMap, v0data::selTOFDeltaTPositiveProtonLambda);
  if (fabs(v0.posTOFDeltaTLaPi()) < v0sels.getmaxDeltaTimePion())
    bitset(bitMap, v0data::selTOFDeltaTPositivePionLambda);
  if (fabs(v0.posTOFDeltaTK0Pi()) < v0sels.getmaxDeltaTimePion())
    bitset(bitMap, v0data::selTOFDeltaTPositivePionK0Short);
  // Negative track
  if (fabs(v0.negTOFDeltaTLaPr()) < v0sels.getmaxDeltaTimeProton())
    bitset(bitMap, v0data::selTOFDeltaTNegativeProtonLambda);
  if (fabs(v0.negTOFDeltaTLaPi()) < v0sels.getmaxDeltaTimePion())
    bitset(bitMap, v0data::selTOFDeltaTNegativePionLambda);
  if (fabs(v0.negTOFDeltaTK0Pi()) < v0sels.getmaxDeltaTimePion())
    bitset(bitMap, v0data::selTOFDeltaTNegativePionK0Short);

  // TOF PID in NSigma
  // Positive track
  if (fabs(v0.tofNSigmaLaPr()) < v0sels.getTofPidNsigmaCutLaPr())
    bitset(bitMap, v0data::selTOFNSigmaPositiveProtonLambda);
  if (fabs(v0.tofNSigmaALaPi()) < v0sels.getTofPidNsigmaCutLaPi())
    bitset(bitMap, v0data::selTOFNSigmaPositivePionLambda);
  if (fabs(v0.tofNSigmaK0PiPlus()) < v0sels.getTofPidNsigmaCutK0Pi())
    bitset(bitMap, v0data::selTOFNSigmaPositivePionK0Short);
  // Negative track
  if (fabs(v0.tofNSigmaALaPr()) < v0sels.getTofPidNsigmaCutLaPr())
    bitset(bitMap, v0data::selTOFNSigmaNegativeProtonLambda);
  if (fabs(v0.tofNSigmaLaPi()) < v0sels.getTofPidNsigmaCutLaPi())
    bitset(bitMap, v0data::selTOFNSigmaNegativePionLambda);
  if (fabs(v0.tofNSigmaK0PiMinus()) < v0sels.getTofPidNsigmaCutK0Pi())
    bitset(bitMap, v0data::selTOFNSigmaNegativePionK0Short);

  // ITS only tag
  if (posTrackExtra.tpcCrossedRows() < 1)
    bitset(bitMap, v0data::selPosItsOnly);
  if (negTrackExtra.tpcCrossedRows() < 1)
    bitset(bitMap, v0data::selNegItsOnly);

  // TPC only tag
  if (posTrackExtra.detectorMap() != o2::aod::track::TPC)
    bitset(bitMap, v0data::selPosNotTPCOnly);
  if (negTrackExtra.detectorMap() != o2::aod::track::TPC)
    bitset(bitMap, v0data::selNegNotTPCOnly);

  // proper lifetime
  if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0 < v0sels.getlifetimeCutLambda())
    bitset(bitMap, v0data::selLambdaCTau);
  if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short < v0sels.getlifetimeCutK0Short())
    bitset(bitMap, v0data::selK0ShortCTau);

  // armenteros
  if (v0.qtarm() * v0sels.getarmPodCut() > TMath::Abs(v0.alpha()) || v0sels.getarmPodCut() < 1e-4)
    bitset(bitMap, v0data::selK0ShortArmenteros);

  return bitMap;
}
} // namespace v0data

#endif // PWGLF_UTILS_V0SELECTIONTOOLS_H_
