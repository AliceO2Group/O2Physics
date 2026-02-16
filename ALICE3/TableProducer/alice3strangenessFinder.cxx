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

#include "PWGLF/DataModel/LFStrangenessTables.h"

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
  Produces<aod::StoredCascCores> tableCascadeCores;
  Produces<aod::CascIndices> tableCascadeIndices;

  Configurable<float> nSigmaTOF{"nSigmaTOF", 5.0f, "Nsigma for TOF PID (if enabled)"};
  Configurable<float> dcaXYconstant{"dcaXYconstant", -1.0f, "[0] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> dcaXYpTdep{"dcaXYpTdep", 0.0, "[1] in |DCAxy| > [0]+[1]/pT"};

  Configurable<float> bachMinConstDCAxy{"bachMinConstDCAxy", -1.0f, "[0] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> bachMinPtDepDCAxy{"bachMinPtDepDCAxy", 0.0, "[1] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> bachMinConstDCAz{"bachMinConstDCAz", -1.0f, "[0] in |DCAz| > [0]+[1]/pT"};
  Configurable<float> bachMinPtDepDCAz{"bachMinPtDepDCAz", 0.0, "[1] in |DCAz| > [0]+[1]/pT"};

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
  Partition<Alice3Tracks> bachelorTracks =
    nabs(aod::track::dcaXY) > bachMinConstDCAxy + bachMinPtDepDCAxy* nabs(aod::track::signed1Pt) && nabs(aod::track::dcaZ) > bachMinConstDCAz + bachMinPtDepDCAz* nabs(aod::track::signed1Pt);

  // Partition<Alice3TracksWPid> negativeSecondaryPions = nabs(aod::upgrade_tof::nSigmaPionInnerTOF) < nSigmaTOF && nabs(aod::upgrade_tof::nSigmaPionOuterTOF) < nSigmaTOF && aod::track::signed1Pt < 0.0f && nabs(aod::track::dcaXY) > dcaXYconstant + dcaXYpTdep* nabs(aod::track::signed1Pt);
  // Partition<Alice3TracksWPid> positiveSecondaryPions = nabs(aod::upgrade_tof::nSigmaPionInnerTOF) < nSigmaTOF && nabs(aod::upgrade_tof::nSigmaPionOuterTOF) < nSigmaTOF && aod::track::signed1Pt > 0.0f && nabs(aod::track::dcaXY) > dcaXYconstant + dcaXYpTdep* nabs(aod::track::signed1Pt);
  // Partition<Alice3TracksWPid> secondaryProtons = nabs(aod::upgrade_tof::nSigmaProtonInnerTOF) < nSigmaTOF && nabs(aod::upgrade_tof::nSigmaProtonOuterTOF) < nSigmaTOF && aod::track::signed1Pt > 0.0f && nabs(aod::track::dcaXY) > dcaXYconstant + dcaXYpTdep* nabs(aod::track::signed1Pt);
  // Partition<Alice3TracksWPid> secondaryAntiProtons = nabs(aod::upgrade_tof::nSigmaProtonInnerTOF) < nSigmaTOF && nabs(aod::upgrade_tof::nSigmaProtonOuterTOF) < nSigmaTOF && aod::track::signed1Pt < 0.0f && nabs(aod::track::dcaXY) > dcaXYconstant + dcaXYpTdep* nabs(aod::track::signed1Pt);

  struct Candidate {
    // decay properties
    float dcaDau{};
    float eta{};
    std::array<float, 3> p{};
    std::array<float, 3> posSV{};
    std::array<float, 3> pDau0{};
    std::array<float, 3> pDau1{};
    std::array<float, o2::track::kLabCovMatSize> parentTrackCovMatrix{};
    float cosPA{};
    float dcaToPV{};
  };

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

    auto hV0Counter = histos.add<TH1>("hV0Counter", "hV0Counter", kTH1D, {{4, 0, 4}});
    hV0Counter->GetXaxis()->SetBinLabel(1, "K0S");
    hV0Counter->GetXaxis()->SetBinLabel(2, "Lambda");
    hV0Counter->GetXaxis()->SetBinLabel(3, "AntiLambda");
    hV0Counter->GetXaxis()->SetBinLabel(4, "Misidentified");

    auto hCascadeCounter = histos.add<TH1>("hCascadeCounter", "hCascadeCounter", kTH1D, {{5, 0, 5}});
    hCascadeCounter->GetXaxis()->SetBinLabel(1, "Xi");
    hCascadeCounter->GetXaxis()->SetBinLabel(2, "AntiXi");
    hCascadeCounter->GetXaxis()->SetBinLabel(3, "Omega");
    hCascadeCounter->GetXaxis()->SetBinLabel(4, "AntiOmega");
    hCascadeCounter->GetXaxis()->SetBinLabel(5, "Misidentified");
  }

  float calculateDCAStraightToPV(float X, float Y, float Z, float Px, float Py, float Pz, float pvX, float pvY, float pvZ)
  {
    return std::sqrt((std::pow((pvY - Y) * Pz - (pvZ - Z) * Py, 2) + std::pow((pvX - X) * Pz - (pvZ - Z) * Px, 2) + std::pow((pvX - X) * Py - (pvY - Y) * Px, 2)) / (Px * Px + Py * Py + Pz * Pz));
  }

  /// function to check if tracks have the same mother in MC
  template <typename TTrackType>
  bool checkSameMother(TTrackType const& track1, TTrackType const& track2)
  {
    bool returnValue = false;
    // Association check
    if (track1.has_mcParticle() && track2.has_mcParticle()) {
      auto mcParticle1 = track1.template mcParticle_as<aod::McParticles>();
      auto mcParticle2 = track2.template mcParticle_as<aod::McParticles>();
      if (mcParticle1.has_mothers() && mcParticle2.has_mothers()) {
        for (const auto& mcParticleMother1 : mcParticle1.template mothers_as<aod::McParticles>()) {
          for (const auto& mcParticleMother2 : mcParticle2.template mothers_as<aod::McParticles>()) {
            if (mcParticleMother1.globalIndex() == mcParticleMother2.globalIndex()) {
              returnValue = true;
            }
          }
        }
      }
    } // end association check
    return returnValue;
  }

  template <typename TTrackType>
  bool buildDecayCandidateTwoBody(TTrackType const& t0, TTrackType const& t1, std::array<float, 3> vtx, Candidate& thisCandidate)
  {
    //}-{}-{}-{}-{}-{}-{}-{}-{}-{}
    // Move close to minima
    int nCand = 0;
    try {
      nCand = fitter.process(t0, t1);
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

    o2::track::TrackParCov t0New = fitter.getTrack(0);
    o2::track::TrackParCov t1New = fitter.getTrack(1);
    t0New.getPxPyPzGlo(thisCandidate.pDau0);
    t1New.getPxPyPzGlo(thisCandidate.pDau1);

    thisCandidate.dcaDau = std::sqrt(fitter.getChi2AtPCACandidate());
    thisCandidate.p[0] = thisCandidate.pDau0[0] + thisCandidate.pDau1[0];
    thisCandidate.p[1] = thisCandidate.pDau0[1] + thisCandidate.pDau1[1];
    thisCandidate.p[2] = thisCandidate.pDau0[2] + thisCandidate.pDau1[2];

    const auto posSV = fitter.getPCACandidatePos();
    thisCandidate.posSV[0] = posSV[0];
    thisCandidate.posSV[1] = posSV[1];
    thisCandidate.posSV[2] = posSV[2];

    std::array<float, o2::track::kLabCovMatSize> covA = {0};
    std::array<float, o2::track::kLabCovMatSize> covB = {0};
    fitter.getTrack(0).getCovXYZPxPyPzGlo(covA);
    fitter.getTrack(1).getCovXYZPxPyPzGlo(covB);

    static constexpr std::array<int, 6> MomentumIndices = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
    for (size_t i = 0; i < MomentumIndices.size(); i++) {
      int j = MomentumIndices[i];
      thisCandidate.parentTrackCovMatrix[j] = covA[j] + covB[j];
    }

    auto covVtx = fitter.calcPCACovMatrix();
    thisCandidate.parentTrackCovMatrix[0] = covVtx(0, 0);
    thisCandidate.parentTrackCovMatrix[1] = covVtx(1, 0);
    thisCandidate.parentTrackCovMatrix[2] = covVtx(1, 1);
    thisCandidate.parentTrackCovMatrix[3] = covVtx(2, 0);
    thisCandidate.parentTrackCovMatrix[4] = covVtx(2, 1);
    thisCandidate.parentTrackCovMatrix[5] = covVtx(2, 2);

    thisCandidate.eta = RecoDecay::eta(std::array{thisCandidate.p[0], thisCandidate.p[1], thisCandidate.p[2]});
    thisCandidate.cosPA = RecoDecay::cpa(vtx, std::array{thisCandidate.posSV[0], thisCandidate.posSV[1], thisCandidate.posSV[2]},
                                         std::array{thisCandidate.p[0], thisCandidate.p[1], thisCandidate.p[2]});
    thisCandidate.dcaToPV = calculateDCAStraightToPV(thisCandidate.posSV[0], thisCandidate.posSV[1], thisCandidate.posSV[2],
                                                     thisCandidate.p[0], thisCandidate.p[1], thisCandidate.p[2],
                                                     vtx[0], vtx[1], vtx[2]);

    return true;
  }

  void processFindV0CandidateNoPid(aod::Collision const& collision, Alice3Tracks const&, aod::McParticles const&)
  {
    auto negativeSecondaryTracksGrouped = negativeSecondaryTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto positiveSecondaryTracksGrouped = positiveSecondaryTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto bachelorTracksGrouped = bachelorTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    const std::array<float, 3> vtx = {collision.posX(), collision.posY(), collision.posZ()};

    for (auto const& posTrack : positiveSecondaryTracksGrouped) {
      if (!posTrack.isReconstructed()) {
        continue; // no ghost tracks
      }

      o2::track::TrackParCov pos = getTrackParCov(posTrack);
      for (auto const& negTrack : negativeSecondaryTracksGrouped) {
        if (!negTrack.isReconstructed()) {
          continue; // no ghost tracks
        }

        if (mcSameMotherCheck && !checkSameMother(posTrack, negTrack)) {
          continue; // keep only if same mother
        }

        o2::track::TrackParCov neg = getTrackParCov(negTrack);
        Candidate v0cand;
        if (!buildDecayCandidateTwoBody(pos, neg, vtx, v0cand)) {
          continue; // failed at building candidate
        }

        auto mcParticle1 = posTrack.template mcParticle_as<aod::McParticles>();
        if (mcParticle1.pdgCode() == kK0Short) {
          histos.fill(HIST("hV0Counter"), 0.5);
        } else if (mcParticle1.pdgCode() == kLambda0) {
          histos.fill(HIST("hV0Counter"), 1.5);
        } else if (mcParticle1.pdgCode() == kLambda0Bar) {
          histos.fill(HIST("hV0Counter"), 2.5);
        } else {
          histos.fill(HIST("hV0Counter"), 3.5);
        }

        v0CandidateIndices(collision.globalIndex(),
                           posTrack.globalIndex(),
                           negTrack.globalIndex(),
                           mcParticle1.globalIndex());

        v0CandidateCores(v0cand.posSV[0], v0cand.posSV[1], v0cand.posSV[2],
                         v0cand.pDau0[0], v0cand.pDau0[1], v0cand.pDau0[2],
                         v0cand.pDau1[0], v0cand.pDau1[1], v0cand.pDau1[2],
                         v0cand.dcaDau, posTrack.dcaXY(), negTrack.dcaXY(),
                         v0cand.cosPA, v0cand.dcaToPV);

        o2::track::TrackParCov v0(v0cand.posSV, v0cand.p, v0cand.parentTrackCovMatrix, 0);
        for (const auto& bachTrack : bachelorTracksGrouped) {
          if (bachTrack.globalIndex() == posTrack.globalIndex() || bachTrack.globalIndex() == negTrack.globalIndex()) {
            continue; // avoid using any track that was already used
          }

          // TODO mc same mother check

          Candidate cascCand;
          o2::track::TrackParCov bach = getTrackParCov(bachTrack);
          if (!buildDecayCandidateTwoBody(v0, bach, vtx, cascCand)) {
            continue; // failed at building candidate
          }

          const float massXi = RecoDecay::m(std::array{std::array{cascCand.pDau0[0], cascCand.pDau0[1], cascCand.pDau0[2]},
                                                       std::array{cascCand.pDau1[0], cascCand.pDau1[1], cascCand.pDau1[2]}},
                                            std::array{o2::constants::physics::MassLambda, o2::constants::physics::MassPionCharged});

          const float massOm = RecoDecay::m(std::array{std::array{cascCand.pDau0[0], cascCand.pDau0[1], cascCand.pDau0[2]},
                                                       std::array{cascCand.pDau1[0], cascCand.pDau1[1], cascCand.pDau1[2]}},
                                            std::array{o2::constants::physics::MassLambda, o2::constants::physics::MassKaonCharged});

          tableCascadeIndices(0, // cascade index, dummy value
                              posTrack.globalIndex(),
                              negTrack.globalIndex(),
                              bachTrack.globalIndex(),
                              collision.globalIndex());

          const float dcaPosToPV = calculateDCAStraightToPV(posTrack.x(), posTrack.y(), posTrack.z(),
                                                            posTrack.px(), posTrack.py(), posTrack.pz(),
                                                            vtx[0], vtx[1], vtx[2]);

          const float dcaNegToPV = calculateDCAStraightToPV(negTrack.x(), negTrack.y(), negTrack.z(),
                                                            negTrack.px(), negTrack.py(), negTrack.pz(),
                                                            vtx[0], vtx[1], vtx[2]);

          const float dcaBachToPV = calculateDCAStraightToPV(bachTrack.x(), bachTrack.y(), bachTrack.z(),
                                                             bachTrack.px(), bachTrack.py(), bachTrack.pz(),
                                                             vtx[0], vtx[1], vtx[2]);

          tableCascadeCores(bachTrack.sign(), massXi, massOm,
                            cascCand.posSV[0], cascCand.posSV[1], cascCand.posSV[2],
                            v0cand.posSV[0], v0cand.posSV[1], v0cand.posSV[2],
                            v0cand.pDau0[0], v0cand.pDau0[1], v0cand.pDau0[2],
                            v0cand.pDau1[0], v0cand.pDau1[1], v0cand.pDau1[2],
                            cascCand.pDau1[0], cascCand.pDau1[1], cascCand.pDau1[2],
                            cascCand.p[0], cascCand.p[1], cascCand.p[2],
                            v0cand.dcaDau, cascCand.dcaDau,
                            dcaPosToPV, dcaNegToPV, dcaBachToPV,
                            cascCand.dcaToPV, cascCand.dcaToPV);

          auto mcParticle2 = bachTrack.template mcParticle_as<aod::McParticles>();
          if (mcParticle2.pdgCode() == PDG_t::kXiMinus) {
            histos.fill(HIST("hCascadeCounter"), 0.5);
          } else if (mcParticle2.pdgCode() == PDG_t::kXiPlusBar) {
            histos.fill(HIST("hCascadeCounter"), 1.5);
          } else if (mcParticle2.pdgCode() == PDG_t::kOmegaMinus) {
            histos.fill(HIST("hCascadeCounter"), 2.5);
          } else if (mcParticle2.pdgCode() == PDG_t::kOmegaPlusBar) {
            histos.fill(HIST("hCascadeCounter"), 3.5);
          } else {
            histos.fill(HIST("hCascadeCounter"), 4.5);
          }
        } // end bachTrack
      } // end negTrack
    } // end posTrack
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
