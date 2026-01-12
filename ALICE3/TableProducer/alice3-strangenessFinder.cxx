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

/// \file alice3-strangenessFinder.cxx
///
/// \brief finding of V0 and cascade candidates for ALICE 3
///
/// This task finds and build condidates for strange hadrons (K0s, Lambda, AntiLambda, Xi-, Xi+, Omega-, Omega+)
/// using the output of the on-the-fly tracker.
///
/// \author Lucia Anna Tarasovičová, Pavol Jozef Šafárik University (SK)
///

#include "ALICE3/DataModel/OTFPIDTrk.h"
#include "ALICE3/DataModel/OTFRICH.h"
#include "ALICE3/DataModel/OTFStrangeness.h"
#include "ALICE3/DataModel/OTFTOF.h"
#include "ALICE3/DataModel/tracksAlice3.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "DCAFitter/DCAFitterN.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"
#include <Framework/AnalysisHelpers.h>

#include <cstdlib>

using namespace o2;
// using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::constants::physics;

using Alice3TracksWPid = soa::Join<aod::Tracks, aod::TracksCov, aod::McTrackLabels, aod::TracksDCA, aod::UpgradeTrkPids, aod::UpgradeTofs, aod::UpgradeRichs>;
using Alice3Tracks = soa::Join<aod::Tracks, aod::TracksCov, aod::McTrackLabels, aod::TracksDCA, aod::TracksCovExtension>;

struct alice3strangenessFinder {
  SliceCache cache;

  Produces<aod::V0CandidateIndices> v0CandidateIndices; // contains V0 candidate indices
  Produces<aod::V0CandidateCores> v0CandidateCores;     // contains V0 candidate core information

  Configurable<bool> mcSameMotherCheck{"mcSameMotherCheck", true, "check if tracks come from the same MC mother"};

  Configurable<float> nSigmaTOF{"nSigmaTOF", 5.0f, "Nsigma for TOF PID (if enabled)"};
  Configurable<float> dcaXYconstant{"dcaXYconstant", -1.0f, "[0] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> dcaXYpTdep{"dcaXYpTdep", 0.0, "[1] in |DCAxy| > [0]+[1]/pT"};

  o2::vertexing::DCAFitterN<2> fitter;
  o2::vertexing::DCAFitterN<3> fitter3;

  // partitions for D mesons
  Partition<Alice3Tracks> positiveSecondaryTracks =
    aod::track::signed1Pt > 0.0f && nabs(aod::track::dcaXY) > dcaXYconstant + dcaXYpTdep* nabs(aod::track::signed1Pt);
  Partition<Alice3Tracks> negativeSecondaryTracks =
    aod::track::signed1Pt < 0.0f && nabs(aod::track::dcaXY) > dcaXYconstant + dcaXYpTdep* nabs(aod::track::signed1Pt);
  // Partition<Alice3TracksWPid> negativeSecondaryPions = nabs(aod::upgrade_tof::nSigmaPionInnerTOF) < nSigmaTOF && nabs(aod::upgrade_tof::nSigmaPionOuterTOF) < nSigmaTOF && aod::track::signed1Pt < 0.0f && nabs(aod::track::dcaXY) > dcaXYconstant + dcaXYpTdep* nabs(aod::track::signed1Pt);
  // Partition<Alice3TracksWPid> positiveSecondaryPions = nabs(aod::upgrade_tof::nSigmaPionInnerTOF) < nSigmaTOF && nabs(aod::upgrade_tof::nSigmaPionOuterTOF) < nSigmaTOF && aod::track::signed1Pt > 0.0f && nabs(aod::track::dcaXY) > dcaXYconstant + dcaXYpTdep* nabs(aod::track::signed1Pt);
  // Partition<Alice3TracksWPid> secondaryProtons = nabs(aod::upgrade_tof::nSigmaProtonInnerTOF) < nSigmaTOF && nabs(aod::upgrade_tof::nSigmaProtonOuterTOF) < nSigmaTOF && aod::track::signed1Pt > 0.0f && nabs(aod::track::dcaXY) > dcaXYconstant + dcaXYpTdep* nabs(aod::track::signed1Pt);
  // Partition<Alice3TracksWPid> secondaryAntiProtons = nabs(aod::upgrade_tof::nSigmaProtonInnerTOF) < nSigmaTOF && nabs(aod::upgrade_tof::nSigmaProtonOuterTOF) < nSigmaTOF && aod::track::signed1Pt < 0.0f && nabs(aod::track::dcaXY) > dcaXYconstant + dcaXYpTdep* nabs(aod::track::signed1Pt);

  struct {
    float dcaDau;
    std::array<float, 3> posSV;
    std::array<float, 3> P;
    std::array<float, 3> Pdaug; // positive track
    std::array<float, 3> Ndaug; // negative track
    float cosPA;
    float dcaToPV;
  } v0cand;

  void init(InitContext&)
  {
    // Initialization code here
  }
  /// function to check if tracks have the same mother in MC
  template <typename TTrackType>
  bool checkSameMother(TTrackType const& track1, TTrackType const& track2)
  {
    bool sameMother = false;
    if (!track1.has_mcParticle() || !track2.has_mcParticle())
      return sameMother;
    auto mcParticle1 = track1.template mcParticle_as<aod::McParticles>();
    auto mcParticle2 = track2.template mcParticle_as<aod::McParticles>();
    if (!mcParticle1.has_mothers() || !mcParticle2.has_mothers())
      return sameMother;
    for (auto& mcParticleMother1 : mcParticle1.template mothers_as<aod::McParticles>()) {
      for (auto& mcParticleMother2 : mcParticle2.template mothers_as<aod::McParticles>()) {
        if (mcParticleMother1.globalIndex() == mcParticleMother2.globalIndex()) {
          sameMother = true;
        }
      }
    }
    return sameMother;
  }

  template <typename TTrackType>
  bool buildDecayCandidateTwoBody(TTrackType const& posTrack, TTrackType const& negTrack)
  {
    o2::track::TrackParCov posTrackCov = getTrackParCov(posTrack);
    o2::track::TrackParCov negTrackCov = getTrackParCov(negTrack);

    //}-{}-{}-{}-{}-{}-{}-{}-{}-{}
    // Move close to minima
    int nCand = 0;
    try {
      nCand = fitter.process(posTrackCov, negTrackCov);
    } catch (...) {
      return false;
    }
    if (nCand == 0) {
      return false;
    }
    //}-{}-{}-{}-{}-{}-{}-{}-{}-{}

    posTrackCov = fitter.getTrack(0);
    negTrackCov = fitter.getTrack(1);
    std::array<float, 3> posP;
    std::array<float, 3> negP;
    posTrackCov.getPxPyPzGlo(posP);
    negTrackCov.getPxPyPzGlo(negP);
    v0cand.dcaDau = TMath::Sqrt(fitter.getChi2AtPCACandidate());
    v0cand.Pdaug[0] = posP[0];
    v0cand.Pdaug[1] = posP[1];
    v0cand.Pdaug[2] = posP[2];
    v0cand.Ndaug[0] = negP[0];
    v0cand.Ndaug[1] = negP[1];
    v0cand.Ndaug[2] = negP[2];
    v0cand.P[0] = posP[0] + negP[0];
    v0cand.P[1] = posP[1] + negP[1];
    v0cand.P[2] = posP[2] + negP[2];
    const auto posSV = fitter.getPCACandidate();
    v0cand.posSV[0] = posSV[0];
    v0cand.posSV[1] = posSV[1];
    v0cand.posSV[2] = posSV[2];

    return true;
  }
  float CalculateDCAStraightToPV(float X, float Y, float Z, float Px, float Py, float Pz, float pvX, float pvY, float pvZ)
  {
    return std::sqrt((std::pow((pvY - Y) * Pz - (pvZ - Z) * Py, 2) + std::pow((pvX - X) * Pz - (pvZ - Z) * Px, 2) + std::pow((pvX - X) * Py - (pvY - Y) * Px, 2)) / (Px * Px + Py * Py + Pz * Pz));
  }
  void processFindV0CandidateNoPid(aod::Collision const& collision, Alice3Tracks const&, aod::McParticles const&)
  {
    auto negativeSecondaryTracksGrouped = negativeSecondaryTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto positiveSecondaryTracksGrouped = positiveSecondaryTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    for (auto const& posTrack : positiveSecondaryTracksGrouped) {
      for (auto const& negTrack : negativeSecondaryTracksGrouped) {

        if (mcSameMotherCheck && !checkSameMother(posTrack, negTrack))
          continue;
        if (!buildDecayCandidateTwoBody(posTrack, negTrack))
          continue;
        v0cand.cosPA = RecoDecay::cpa(std::array{collision.posX(), collision.posY(), collision.posZ()}, std::array{v0cand.posSV[0], v0cand.posSV[1], v0cand.posSV[2]}, std::array{v0cand.P[0], v0cand.P[1], v0cand.P[2]});
        v0cand.dcaToPV = CalculateDCAStraightToPV(
          v0cand.posSV[0], v0cand.posSV[1], v0cand.posSV[2],
          v0cand.P[0], v0cand.P[1], v0cand.P[2],
          collision.posX(), collision.posY(), collision.posZ());
        v0CandidateIndices(collision.globalIndex(),
                           posTrack.globalIndex(),
                           negTrack.globalIndex());
        v0CandidateCores(
          v0cand.posSV[0], v0cand.posSV[1], v0cand.posSV[2],
          v0cand.Pdaug[0], v0cand.Pdaug[1], v0cand.Pdaug[2],
          v0cand.Ndaug[0], v0cand.Ndaug[1], v0cand.Ndaug[2],
          v0cand.dcaDau, posTrack.dcaXY(), negTrack.dcaXY(),
          v0cand.cosPA, v0cand.dcaToPV);
      }
    }
  }
  //    void processFindV0CandidateWithPid(aod::Collision const& collision, aod::McParticles const& mcParticles, Alice3TracksWPid const&)
  //     {
  //         auto negativeSecondaryPionsGrouped = negativeSecondaryPions->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
  //         auto positiveSecondaryPionsGrouped = positiveSecondaryPions->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
  //         auto secondaryProtonsGrouped = secondaryProtons->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
  //         auto secondaryAntiProtonsGrouped = secondaryAntiProtons->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
  //     }
  PROCESS_SWITCH(alice3strangenessFinder, processFindV0CandidateNoPid, "find V0 without PID", true);
  // PROCESS_SWITCH(alice3strangenessFinder, processFindV0CandidateWithPid, "find V0 with PID", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<alice3strangenessFinder>(cfgc)};
}
