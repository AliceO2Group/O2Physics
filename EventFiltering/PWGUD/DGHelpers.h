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

#ifndef O2_ANALYSIS_DGHELPERS_
#define O2_ANALYSIS_DGHELPERS_

#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "Framework/DataTypes.h"
#include "CommonConstants/LHCConstants.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "EventFiltering/PWGUD/DGCutparHolder.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

float FV0AmplitudeA(aod::FV0A&& fv0);
template <typename T>
bool cleanFV0(T& bc, float limitA);

float FT0AmplitudeA(aod::FT0&& ft0);
float FT0AmplitudeC(aod::FT0&& ft0);
template <typename T>
bool cleanFT0(T& bc, float limitA, float limitC);

int16_t FDDAmplitudeA(aod::FDD&& fdd);
int16_t FDDAmplitudeC(aod::FDD&& fdd);
template <typename T>
bool cleanFDD(T& bc, float limitA, float limitC);

template <typename T>
bool cleanFIT(T& bc, std::vector<float> lims);
template <typename T>
bool cleanFITCollision(T& col, std::vector<float> lims);

template <typename T>
bool cleanZDC(T& bc, aod::Zdcs& zdcs, std::vector<float>& lims);

template <typename T>
bool cleanCalo(T& bc, aod::Calos& calos, std::vector<float>& lims);

template <typename TC>
bool hasGoodPID(DGCutparHolder diffCuts, TC track);

template <typename T>
T compatibleBCs(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, int ndt, T const& bcs, int nMinBCs = 7);

// -----------------------------------------------------------------------------
// add here Selectors for different types of diffractive events
// Selector for Double Gap events
struct DGSelector {
 public:
  DGSelector()
  {
    fPDG = TDatabasePDG::Instance();
  }

  // Function to check if collisions passes filter
  template <typename CC, typename BCs, typename TCs, typename FWs>
  int IsSelected(DGCutparHolder diffCuts, CC const& collision, BCs& bcRange, TCs& tracks, FWs& fwdtracks)
  {
    LOGF(debug, "Collision %f", collision.collisionTime());
    LOGF(debug, "Number of close BCs: %i", bcRange.size());

    // check that there are no FIT signals in any of the compatible BCs
    // Double Gap (DG) condition
    auto lims = diffCuts.FITAmpLimits();

    for (auto& bc : bcRange) {
      LOGF(debug, "Amplitudes FV0A %f FT0 %f / %f FDD %i / %i",
           bc.has_foundFV0() ? FV0AmplitudeA(bc.foundFV0()) : -1.,
           bc.has_foundFT0() ? FT0AmplitudeA(bc.foundFT0()) : -1.,
           bc.has_foundFT0() ? FT0AmplitudeC(bc.foundFT0()) : -1.,
           bc.has_foundFDD() ? FDDAmplitudeA(bc.foundFDD()) : -1,
           bc.has_foundFDD() ? FDDAmplitudeC(bc.foundFDD()) : -1);
      LOGF(debug, "  clean FV0A %i FT0 %i FDD %i", cleanFV0(bc, lims[0]), cleanFT0(bc, lims[1], lims[2]), cleanFDD(bc, lims[3], lims[4]));

      if (!cleanFIT(bc, diffCuts.FITAmpLimits())) {
        return 1;
      }
    }

    // no activity in muon arm
    LOGF(debug, "Muons %i", fwdtracks.size());
    for (auto& muon : fwdtracks) {
      LOGF(debug, "  %i / %f / %f / %f", muon.trackType(), muon.eta(), muon.pt(), muon.p());
    }
    if (fwdtracks.size() > 0) {
      return 2;
    }

    // no global tracks which are not vtx tracks
    // no vtx tracks which are not global tracks
    auto rgtrwTOF = 0.;
    for (auto& track : tracks) {
      if (track.isGlobalTrack() && !track.isPVContributor()) {
        return 3;
      }
      if (diffCuts.globalTracksOnly() && !track.isGlobalTrack() && track.isPVContributor()) {
        return 4;
      }

      // update fraction of PV tracks with TOF hit
      if (track.isPVContributor() && track.hasTOF()) {
        rgtrwTOF += 1.;
      }
    }
    if (collision.numContrib() > 0) {
      rgtrwTOF /= collision.numContrib();
    }
    if (rgtrwTOF < diffCuts.minRgtrwTOF()) {
      return 5;
    }

    // number of vertex tracks
    if (collision.numContrib() < diffCuts.minNTracks() || collision.numContrib() > diffCuts.maxNTracks()) {
      return 6;
    }

    // PID, pt, and eta of tracks, invariant mass, and net charge
    // consider only vertex tracks

    // which particle hypothesis?
    auto mass2Use = 0.;
    TParticlePDG* pdgparticle = fPDG->GetParticle(diffCuts.pidHypothesis());
    if (pdgparticle != nullptr) {
      mass2Use = pdgparticle->Mass();
    }

    auto netCharge = 0;
    auto lvtmp = TLorentzVector();
    auto ivm = TLorentzVector();
    for (auto& track : tracks) {
      if (track.isPVContributor()) {

        // PID
        if (!hasGoodPID(diffCuts, track)) {
          return 7;
        }

        // pt
        lvtmp.SetXYZM(track.px(), track.py(), track.pz(), mass2Use);
        if (lvtmp.Perp() < diffCuts.minPt() || lvtmp.Perp() > diffCuts.maxPt()) {
          return 8;
        }

        // eta
        if (lvtmp.Eta() < diffCuts.minEta() || lvtmp.Eta() > diffCuts.maxEta()) {
          return 9;
        }
        netCharge += track.sign();
        ivm += lvtmp;
      }
    }

    // net charge
    if (netCharge < diffCuts.minNetCharge() || netCharge > diffCuts.maxNetCharge()) {
      return 10;
    }
    // invariant mass
    if (ivm.M() < diffCuts.minIVM() || ivm.M() > diffCuts.maxIVM()) {
      return 11;
    }

    // if we arrive here then the event is good!
    return 0;
  };

  template <typename BCs, typename TCs>
  int IsSelected(DGCutparHolder diffCuts, BCs& bc, TCs& tracks)
  {
    // check that there are no FIT signals in bc
    // Double Gap (DG) condition
    if (!cleanFIT(bc, diffCuts.FITAmpLimits())) {
      return 1;
    }

    // number of tracks
    if ((int)tracks.size() < diffCuts.minNTracks() || (int)tracks.size() > diffCuts.maxNTracks()) {
      return 6;
    }

    // PID, pt, and eta of tracks, invariant mass, and net charge
    // which particle hypothesis?
    auto mass2Use = 0.;
    TParticlePDG* pdgparticle = fPDG->GetParticle(diffCuts.pidHypothesis());
    if (pdgparticle != nullptr) {
      mass2Use = pdgparticle->Mass();
    }

    auto netCharge = 0;
    auto lvtmp = TLorentzVector();
    auto ivm = TLorentzVector();
    for (auto& track : tracks) {
      // PID
      if (!hasGoodPID(diffCuts, track)) {
        return 7;
      }

      // pt
      lvtmp.SetXYZM(track.px(), track.py(), track.pz(), mass2Use);
      if (lvtmp.Perp() < diffCuts.minPt() || lvtmp.Perp() > diffCuts.maxPt()) {
        return 8;
      }

      // eta
      if (lvtmp.Eta() < diffCuts.minEta() || lvtmp.Eta() > diffCuts.maxEta()) {
        return 9;
      }
      netCharge += track.sign();
      ivm += lvtmp;
    }

    // net charge
    if (netCharge < diffCuts.minNetCharge() || netCharge > diffCuts.maxNetCharge()) {
      return 10;
    }
    // invariant mass
    if (ivm.M() < diffCuts.minIVM() || ivm.M() > diffCuts.maxIVM()) {
      return 11;
    }

    // if we arrive here then the event is good!
    return 0;
  };

 private:
  TDatabasePDG* fPDG;
};

// -----------------------------------------------------------------------------
// The associations between collsisions and BCs can be ambiguous.
// By default a collision is associated with the BC closest in time.
// The collision time t_coll is determined by the tracks which are used to
// reconstruct the vertex. t_coll has an uncertainty dt_coll.
// Any BC with a BC time t_BC falling within a time window of +- ndt*dt_coll
// around t_coll could potentially be the true BC. ndt is typically 4. The
// total width of the time window is required to be at least 2*nMinBCs* LHCBunchSpacingNS.

template <typename T>
T compatibleBCs(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, int ndt, T const& bcs, int nMinBCs)
{
  LOGF(debug, "Collision time / resolution [ns]: %f / %f", collision.collisionTime(), collision.collisionTimeRes());

  // return if collisions has no associated BC
  if (!collision.has_foundBC()) {
    return T{{bcs.asArrowTable()->Slice(0, 0)}, (uint64_t)0};
  }

  // get associated BC
  auto bcIter = collision.foundBC_as<T>();

  // due to the filling scheme the most probable BC may not be the one estimated from the collision time
  uint64_t mostProbableBC = bcIter.globalBC();
  uint64_t meanBC = mostProbableBC + std::lround(collision.collisionTime() / o2::constants::lhc::LHCBunchSpacingNS);

  // enforce minimum number for deltaBC
  int deltaBC = std::ceil(collision.collisionTimeRes() / o2::constants::lhc::LHCBunchSpacingNS * ndt);
  if (deltaBC < nMinBCs) {
    deltaBC = nMinBCs;
  }

  int64_t minBC = meanBC - deltaBC;
  uint64_t maxBC = meanBC + deltaBC;
  if (minBC < 0) {
    minBC = 0;
  }

  // find slice of BCs table with BC in [minBC, maxBC]
  int64_t maxBCId = bcIter.globalIndex();
  int moveCount = 0; // optimize to avoid to re-create the iterator
  while (bcIter != bcs.end() && bcIter.globalBC() <= maxBC && (int64_t)bcIter.globalBC() >= minBC) {
    LOGF(debug, "Table id %d BC %llu", bcIter.globalIndex(), bcIter.globalBC());
    maxBCId = bcIter.globalIndex();
    ++bcIter;
    ++moveCount;
  }

  bcIter.moveByIndex(-moveCount); // Move back to original position
  int64_t minBCId = collision.bcId();
  while (bcIter != bcs.begin() && bcIter.globalBC() <= maxBC && (int64_t)bcIter.globalBC() >= minBC) {
    LOGF(debug, "Table id %d BC %llu", bcIter.globalIndex(), bcIter.globalBC());
    minBCId = bcIter.globalIndex();
    --bcIter;
  }

  LOGF(debug, "  BC range: %i (%d) - %i (%d)", minBC, minBCId, maxBC, maxBCId);

  T slice{{bcs.asArrowTable()->Slice(minBCId, maxBCId - minBCId + 1)}, (uint64_t)minBCId};
  bcs.copyIndexBindings(slice);
  return slice;
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
bool cleanFV0(T& bc, float limitA)
{
  if (bc.has_foundFV0()) {
    return (FV0AmplitudeA(bc.foundFV0()) < limitA);
  } else {
    return true;
  }
}

// -----------------------------------------------------------------------------
template <typename T>
bool cleanFT0(T& bc, float limitA, float limitC)
{
  if (bc.has_foundFT0()) {
    return (FT0AmplitudeA(bc.foundFT0()) < limitA) && (FT0AmplitudeC(bc.foundFT0()) < limitC);
  } else {
    return true;
  }
}

// -----------------------------------------------------------------------------
template <typename T>
bool cleanFDD(T& bc, float limitA, float limitC)
{
  if (bc.has_foundFDD()) {
    return (FDDAmplitudeA(bc.foundFDD()) < limitA) && (FDDAmplitudeC(bc.foundFDD()) < limitC);
  } else {
    return true;
  }
}

// -----------------------------------------------------------------------------
template <typename T>
bool cleanFIT(T& bc, std::vector<float> lims)
{
  return cleanFV0(bc, lims[0]) && cleanFT0(bc, lims[1], lims[2]) && cleanFDD(bc, lims[3], lims[4]);
}
template <typename T>
bool cleanFITCollision(T& col, std::vector<float> lims)
{
  bool isCleanFV0 = true;
  if (col.has_foundFV0()) {
    isCleanFV0 = (FV0AmplitudeA(col.foundFV0()) < lims[0]);
  }
  bool isCleanFT0 = true;
  if (col.has_foundFT0()) {
    isCleanFT0 = (FT0AmplitudeA(col.foundFT0()) < lims[1]) && (FT0AmplitudeC(col.foundFT0()) < lims[2]);
  }
  bool isCleanFDD = true;
  if (col.has_foundFDD()) {
    isCleanFDD = (FDDAmplitudeA(col.foundFDD()) < lims[3]) && (FDDAmplitudeC(col.foundFDD()) < lims[4]);
  }
  return (isCleanFV0 && isCleanFT0 && isCleanFDD);
}
// -----------------------------------------------------------------------------
template <typename T>
bool cleanZDC(T& bc, aod::Zdcs& zdcs, std::vector<float>& lims)
{
  const auto& ZdcBC = zdcs.sliceByCached(aod::zdc::bcId, bc.globalIndex());
  return (ZdcBC.size() == 0);
}

// -----------------------------------------------------------------------------
template <typename T>
bool cleanCalo(T& bc, aod::Calos& calos, std::vector<float>& lims)
{
  const auto& CaloBC = calos.sliceByCached(aod::calo::bcId, bc.globalIndex());
  return (CaloBC.size() == 0);
}

// -----------------------------------------------------------------------------
#endif // O2_ANALYSIS_DGHELPERS_
