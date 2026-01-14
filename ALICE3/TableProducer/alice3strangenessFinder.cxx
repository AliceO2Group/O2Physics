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

/// \file alice3strangenessFinder.cxx
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
#include <Framework/Configurable.h>

#include <cstdlib>

using namespace o2;
// using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::constants::physics;

using Alice3TracksWPid = soa::Join<aod::Tracks, aod::TracksCov, aod::McTrackLabels, aod::TracksDCA, aod::UpgradeTrkPids, aod::UpgradeTofs, aod::UpgradeRichs>;
using Alice3Tracks = soa::Join<aod::StoredTracks, aod::StoredTracksCov, aod::McTrackLabels, aod::TracksDCA, aod::TracksCovExtension, aod::TracksAlice3>;

struct Alice3strangenessFinder {
  SliceCache cache;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Produces<aod::V0CandidateIndices> v0CandidateIndices; // contains V0 candidate indices
  Produces<aod::V0CandidateCores> v0CandidateCores;     // contains V0 candidate core information

  Configurable<float> nSigmaTOF{"nSigmaTOF", 5.0f, "Nsigma for TOF PID (if enabled)"};
  Configurable<float> dcaXYconstant{"dcaXYconstant", -1.0f, "[0] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> dcaXYpTdep{"dcaXYpTdep", 0.0, "[1] in |DCAxy| > [0]+[1]/pT"};

  // Vertexing
  Configurable<bool> propagateToPCA{"propagateToPCA", false, "create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCA{"useAbsDCA", true, "Minimise abs. distance rather than chi2"};
  Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  Configurable<double> maxR{"maxR", 150., "reject PCA's above this radius"};
  Configurable<double> maxDZIni{"maxDZIni", 5, "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> maxDXYIni{"maxDXYIni", 4, "reject (if>0) PCA candidate if tracks DXY exceeds threshold"};
  Configurable<double> maxVtxChi2{"maxVtxChi2", 2, "reject (if>0) vtx. chi2 above this value"};
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations is chi2/chi2old > this"};
  // Operation and minimisation criteria
  Configurable<float> magneticField{"magneticField", 20.0f, "Magnetic field (in kilogauss)"};
  Configurable<bool> doDCAplotsD{"doDCAplotsD", true, "do daughter prong DCA plots for D mesons"};
  Configurable<bool> doDCAplots3Prong{"doDCAplots3Prong", true, "do daughter prong DCA plots for Lc baryons"};
  Configurable<bool> doTopoPlotsForSAndB{"doTopoPlotsForSAndB", true, "do topological variable distributions for S and B separately"};
  Configurable<float> dcaDaughtersSelection{"dcaDaughtersSelection", 1000.0f, "DCA between daughters (cm)"};
  Configurable<bool> mcSameMotherCheck{"mcSameMotherCheck", true, "check if tracks come from the same MC mother"};
  // propagation options
  Configurable<bool> usePropagator{"usePropagator", false, "use external propagator"};
  Configurable<bool> refitWithMatCorr{"refitWithMatCorr", false, "refit V0 applying material corrections"};
  Configurable<bool> useCollinearV0{"useCollinearV0", true, "use collinear approximation for V0 fitting"};

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
    std::array<float, 3> pV0;
    std::array<float, 3> pPos; // positive track
    std::array<float, 3> pNeg; // negative track
    float cosPA;
    float dcaToPV;
  } v0cand;

  void init(InitContext&)
  {
    // Initialization code here
    fitter.setBz(magneticField);
    fitter.setUseAbsDCA(useAbsDCA);
    fitter.setPropagateToPCA(propagateToPCA);
    fitter.setMaxR(maxR);
    fitter.setMinParamChange(minParamChange);
    fitter.setMinRelChi2Change(minRelChi2Change);
    fitter.setMaxDZIni(maxDZIni);
    fitter.setMaxDXYIni(maxDXYIni);
    fitter.setMaxChi2(maxVtxChi2);
    fitter.setUsePropagator(usePropagator);
    fitter.setRefitWithMatCorr(refitWithMatCorr);
    fitter.setCollinear(useCollinearV0);
    fitter.setMatCorrType(o2::base::Propagator::MatCorrType::USEMatCorrNONE);

    histos.add("hV0Counter", "", kTH1D, {{4, 0, 4}}); // For QA reasons, counting found V0, 0: K0s, 1: Lambda, 2:AntiLambda, 3: wrongly identified V0
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
    if (mcParticle2.globalIndex() == mcParticle1.globalIndex()) { // for the V0 daughters we store the mc label of the mother particle in the daughter tracks
      sameMother = true;
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
    if (!fitter.isPropagateTracksToVertexDone() && !fitter.propagateTracksToVertex()) {
      LOG(debug) << "RejProp failed";
      return false;
    }

    posTrackCov = fitter.getTrack(0);
    negTrackCov = fitter.getTrack(1);
    std::array<float, 3> posP;
    std::array<float, 3> negP;
    posTrackCov.getPxPyPzGlo(posP);
    negTrackCov.getPxPyPzGlo(negP);
    v0cand.dcaDau = std::sqrt(fitter.getChi2AtPCACandidate());
    v0cand.pPos[0] = posP[0];
    v0cand.pPos[1] = posP[1];
    v0cand.pPos[2] = posP[2];
    v0cand.pNeg[0] = negP[0];
    v0cand.pNeg[1] = negP[1];
    v0cand.pNeg[2] = negP[2];
    v0cand.pV0[0] = posP[0] + negP[0];
    v0cand.pV0[1] = posP[1] + negP[1];
    v0cand.pV0[2] = posP[2] + negP[2];
    const auto posSV = fitter.getPCACandidatePos();
    v0cand.posSV[0] = posSV[0];
    v0cand.posSV[1] = posSV[1];
    v0cand.posSV[2] = posSV[2];

    return true;
  }
  float calculateDCAStraightToPV(float X, float Y, float Z, float Px, float Py, float Pz, float pvX, float pvY, float pvZ)
  {
    return std::sqrt((std::pow((pvY - Y) * Pz - (pvZ - Z) * Py, 2) + std::pow((pvX - X) * Pz - (pvZ - Z) * Px, 2) + std::pow((pvX - X) * Py - (pvY - Y) * Px, 2)) / (Px * Px + Py * Py + Pz * Pz));
  }
  void processFindV0CandidateNoPid(aod::Collision const& collision, Alice3Tracks const&, aod::McParticles const&)
  {
    auto negativeSecondaryTracksGrouped = negativeSecondaryTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto positiveSecondaryTracksGrouped = positiveSecondaryTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    for (auto const& posTrack : positiveSecondaryTracksGrouped) {
      if (!posTrack.isReconstructed()) {
        continue; // no ghost tracks
      }
      for (auto const& negTrack : negativeSecondaryTracksGrouped) {
        if (!negTrack.isReconstructed()) {
          continue; // no ghost tracks
        }
        auto mcParticle1 = posTrack.template mcParticle_as<aod::McParticles>();

        if (mcSameMotherCheck && !checkSameMother(posTrack, negTrack))
          continue;
        if (!buildDecayCandidateTwoBody(posTrack, negTrack))
          continue;
        v0cand.cosPA = RecoDecay::cpa(std::array{collision.posX(), collision.posY(), collision.posZ()}, std::array{v0cand.posSV[0], v0cand.posSV[1], v0cand.posSV[2]}, std::array{v0cand.pV0[0], v0cand.pV0[1], v0cand.pV0[2]});
        v0cand.dcaToPV = calculateDCAStraightToPV(
          v0cand.posSV[0], v0cand.posSV[1], v0cand.posSV[2],
          v0cand.pV0[0], v0cand.pV0[1], v0cand.pV0[2],
          collision.posX(), collision.posY(), collision.posZ());
        v0CandidateIndices(collision.globalIndex(),
                           posTrack.globalIndex(),
                           negTrack.globalIndex(),
                           mcParticle1.globalIndex());
        v0CandidateCores(
          v0cand.posSV[0], v0cand.posSV[1], v0cand.posSV[2],
          v0cand.pPos[0], v0cand.pPos[1], v0cand.pPos[2],
          v0cand.pNeg[0], v0cand.pNeg[1], v0cand.pNeg[2],
          v0cand.dcaDau, posTrack.dcaXY(), negTrack.dcaXY(),
          v0cand.cosPA, v0cand.dcaToPV);
        if (mcParticle1.pdgCode() == kK0Short) {
          histos.fill(HIST("hV0Counter"), 0.5);
        } else if (mcParticle1.pdgCode() == kLambda0) {
          histos.fill(HIST("hV0Counter"), 1.5);
        } else if (mcParticle1.pdgCode() == kLambda0Bar) {
          histos.fill(HIST("hV0Counter"), 2.5);
        } else {
          histos.fill(HIST("hV0Counter"), 3.5);
        }
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
  PROCESS_SWITCH(Alice3strangenessFinder, processFindV0CandidateNoPid, "find V0 without PID", true);
  // PROCESS_SWITCH(alice3strangenessFinder, processFindV0CandidateWithPid, "find V0 with PID", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<Alice3strangenessFinder>(cfgc)};
}
