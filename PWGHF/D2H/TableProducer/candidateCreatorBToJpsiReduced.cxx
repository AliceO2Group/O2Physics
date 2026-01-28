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

/// \file candidateCreatorBToJpsiReduced.cxx
/// \brief Reconstruction of B->J/Psi hadron candidates
///
/// \author Fabrizio Chinu <fabrizio.chinu@cern.ch>, Università degli Studi and INFN Torino
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

#include "PWGHF/D2H/DataModel/ReducedDataModel.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/Utils/utilsTrkCandHf.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"

#include <CommonConstants/PhysicsConstants.h>
#include <DCAFitter/DCAFitterN.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/WorkflowSpec.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/DCA.h>
#include <ReconstructionDataFormats/Track.h>
#include <ReconstructionDataFormats/TrackParametrizationWithError.h>

#include <TH1.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <memory>
#include <stdexcept>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_trkcandsel;

enum DecayChannel : uint8_t {
  B0ToJpsiK0Star = 0,
  BplusToJpsiK,
  BsToJpsiPhi
};

/// Reconstruction of B+ candidates
struct HfCandidateCreatorBToJpsiReduced {
  Produces<aod::HfCandBpJPBase> rowCandidateBpBase;        // table defined in CandidateReconstructionTables.h
  Produces<aod::HfRedBplus2JpsiDaus> rowCandidateBpProngs; // table defined in ReducedDataModel.h
  Produces<aod::HfCandBsJPBase> rowCandidateBsBase;        // table defined in CandidateReconstructionTables.h
  Produces<aod::HfRedBs2JpsiDaus> rowCandidateBsProngs;    // table defined in ReducedDataModel.h

  // vertexing
  Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCA{"useAbsDCA", true, "Minimise abs. distance rather than chi2"};
  Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any B+ is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations is chi2/chi2old > this"};

  Configurable<bool> runJpsiToee{"runJpsiToee", false, "Run analysis for J/Psi to ee (debug)"};
  // selection
  Configurable<double> invMassWindowJpsiHadTolerance{"invMassWindowJpsiHadTolerance", 0.01, "invariant-mass window tolerance for J/Psi K pair preselections (GeV/c2)"};

  float myInvMassWindowJpsiK{1.}, myInvMassWindowJpsiPhi{1.}; // variable that will store the value of invMassWindowJpsiK (defined in dataCreatorJpsiKReduced.cxx)
  double massBplus{0.}, massBs{0.};
  double bz{0.};
  o2::vertexing::DCAFitterN<2> df2; // fitter for B vertex (2-prong vertex fitter)
  o2::vertexing::DCAFitterN<3> df3; // fitter for B vertex (3-prong vertex fitter)
  o2::vertexing::DCAFitterN<4> df4; // fitter for B vertex (4-prong vertex fitter)

  using HfRedCollisionsWithExtras = soa::Join<aod::HfRedCollisions, aod::HfRedCollExtras>;

  Preslice<soa::Join<aod::HfRedJpsis, aod::HfRedJpsiCov>> candsJpsiPerCollision = hf_track_index_reduced::hfRedCollisionId;
  Preslice<soa::Join<aod::HfRedBach0Bases, aod::HfRedBach0Cov>> tracksLf0PerCollision = hf_track_index_reduced::hfRedCollisionId;
  Preslice<soa::Join<aod::HfRedBach1Bases, aod::HfRedBach1Cov>> tracksLf1PerCollision = hf_track_index_reduced::hfRedCollisionId;

  std::shared_ptr<TH1> hCandidatesB, hCandidatesPhi;
  o2::base::Propagator::MatCorrType noMatCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    // invariant-mass window cut
    massBplus = o2::constants::physics::MassBPlus;
    massBs = o2::constants::physics::MassBS;

    // Initialize fitters
    df2.setPropagateToPCA(propagateToPCA);
    df2.setMaxR(maxR);
    df2.setMaxDZIni(maxDZIni);
    df2.setMinParamChange(minParamChange);
    df2.setMinRelChi2Change(minRelChi2Change);
    df2.setUseAbsDCA(useAbsDCA);
    df2.setWeightedFinalPCA(useWeightedFinalPCA);
    df2.setMatCorrType(noMatCorr);

    df3.setPropagateToPCA(propagateToPCA);
    df3.setMaxR(maxR);
    df3.setMaxDZIni(maxDZIni);
    df3.setMinParamChange(minParamChange);
    df3.setMinRelChi2Change(minRelChi2Change);
    df3.setUseAbsDCA(useAbsDCA);
    df3.setWeightedFinalPCA(useWeightedFinalPCA);
    df3.setMatCorrType(noMatCorr);

    df4.setPropagateToPCA(propagateToPCA);
    df4.setMaxR(maxR);
    df4.setMaxDZIni(maxDZIni);
    df4.setMinParamChange(minParamChange);
    df4.setMinRelChi2Change(minRelChi2Change);
    df4.setUseAbsDCA(useAbsDCA);
    df4.setWeightedFinalPCA(useWeightedFinalPCA);
    df4.setMatCorrType(noMatCorr);

    // histograms
    registry.add("hMassJpsi", "J/Psi mass;#it{M}_{#mu#mu} (GeV/#it{c}^{2});Counts", {HistType::kTH1F, {{600, 2.5, 3.7, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassBplusToJpsiK", "2-prong candidates;inv. mass (B^{+} #rightarrow #overline{D^{0}}#pi^{#plus} #rightarrow #pi^{#minus}K^{#plus}#pi^{#plus}) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 3., 8.}}});
    registry.add("hCovPVXX", "2-prong candidates;XX element of cov. matrix of prim. vtx. position (cm^{2});entries", {HistType::kTH1F, {{100, 0., 1.e-4}}});
    registry.add("hCovSVXX", "2-prong candidates;XX element of cov. matrix of sec. vtx. position (cm^{2});entries", {HistType::kTH1F, {{100, 0., 0.2}}});
    registry.add("hEvents", "Events;;entries", HistType::kTH1F, {{1, 0.5, 1.5}});

    /// candidate monitoring
    hCandidatesB = registry.add<TH1>("hFitCandidatesJpsi", "candidates counter", {HistType::kTH1D, {axisCands}});
    hCandidatesPhi = registry.add<TH1>("hFitCandidatesPhi", "candidates counter", {HistType::kTH1D, {axisCands}});
    setLabelHistoCands(hCandidatesB);
    setLabelHistoCands(hCandidatesPhi);
  }

  /// Main function to perform B+ candidate creation
  /// \param collision the collision
  /// \param candsJpsiThisColl J/Psi candidates in this collision
  /// \param tracksLfThisCollisionArr LF tracks in this collision
  /// \param invMass2JpsiHadMin minimum B invariant-mass
  /// \param invMass2JpsiHadMax maximum B invariant-mass
  template <uint8_t DecChannel, typename Cands, typename TTracks0, typename TTracks1, typename Coll>
  void runCandidateCreation(Coll const& collision,
                            Cands const& candsJpsiThisColl,
                            TTracks0 const& tracksLfDau0ThisCollision,
                            TTracks1 const& tracksLfDau1ThisCollision,
                            const float invMass2JpsiHadMin,
                            const float invMass2JpsiHadMax)
  {
    auto primaryVertex = getPrimaryVertex(collision);
    auto covMatrixPV = primaryVertex.getCov();

    // Set the magnetic field from ccdb
    bz = collision.bz();
    df2.setBz(bz);
    df3.setBz(bz);
    df4.setBz(bz);

    for (const auto& candJpsi : candsJpsiThisColl) {
      o2::track::TrackParametrizationWithError<float> trackPosParCov(
        candJpsi.xDauPos(), candJpsi.alphaDauPos(), {candJpsi.yDauPos(), candJpsi.zDauPos(), candJpsi.snpDauPos(), candJpsi.tglDauPos(), candJpsi.signed1PtDauPos()}, 1 /*Charge*/, 1 /*Muon*/);
      o2::track::TrackParametrizationWithError<float> trackNegParCov(
        candJpsi.xDauNeg(), candJpsi.alphaDauNeg(), {candJpsi.yDauNeg(), candJpsi.zDauNeg(), candJpsi.snpDauNeg(), candJpsi.tglDauNeg(), candJpsi.signed1PtDauNeg()}, -1 /*Charge*/, 1 /*Muon*/);

      // ---------------------------------
      // reconstruct J/Psi candidate
      o2::track::TrackParCov trackParCovJpsi{};
      std::array<float, 3> pVecJpsi{};
      registry.fill(HIST("hFitCandidatesJpsi"), SVFitting::BeforeFit);
      try {
        if (df2.process(trackPosParCov, trackNegParCov) == 0) {
          LOG(info) << "DCAFitterN failed to reconstruct J/Psi candidate, skipping the candidate.";
          continue;
        }
      } catch (const std::runtime_error& error) {
        LOG(info) << "Run time error found: " << error.what() << ". DCAFitterN cannot work, skipping the candidate.";
        registry.fill(HIST("hFitCandidatesJpsi"), SVFitting::Fail);
        continue;
      }
      registry.fill(HIST("hFitCandidatesJpsi"), SVFitting::FitOk);

      std::array<float, 3> pVecDauPos{candJpsi.pxDauPos(), candJpsi.pyDauPos(), candJpsi.pzDauPos()};
      std::array<float, 3> pVecDauNeg{candJpsi.pxDauNeg(), candJpsi.pyDauNeg(), candJpsi.pzDauNeg()};

      df2.getTrack(0).getPxPyPzGlo(pVecDauPos);
      df2.getTrack(1).getPxPyPzGlo(pVecDauNeg);
      pVecJpsi = RecoDecay::pVec(pVecDauPos, pVecDauNeg);
      trackParCovJpsi = df2.createParentTrackParCov();
      trackParCovJpsi.setAbsCharge(0); // to be sure

      float invMassJpsi{0.f};
      if (runJpsiToee) {
        invMassJpsi = RecoDecay::m2(std::array{pVecDauPos, pVecDauNeg}, std::array{o2::constants::physics::MassElectron, o2::constants::physics::MassElectron});
      } else {
        invMassJpsi = RecoDecay::m2(std::array{pVecDauPos, pVecDauNeg}, std::array{o2::constants::physics::MassMuon, o2::constants::physics::MassMuon});
      }
      invMassJpsi = std::sqrt(invMassJpsi);
      registry.fill(HIST("hMassJpsi"), invMassJpsi);

      for (const auto& trackLf0 : tracksLfDau0ThisCollision) {
        // this track is among daughters
        if (trackLf0.trackId() == candJpsi.prongPosId() || trackLf0.trackId() == candJpsi.prongNegId()) {
          continue;
        }
        auto trackParCovLf0 = getTrackParCov(trackLf0);
        std::array<float, 3> pVecTrackLf0{};
        if constexpr (DecChannel == DecayChannel::BplusToJpsiK) {
          // ---------------------------------
          // reconstruct the 3-prong B+ vertex
          hCandidatesB->Fill(SVFitting::BeforeFit);
          try {
            if (df3.process(trackPosParCov, trackNegParCov, trackParCovLf0) == 0) {
              continue;
            }
          } catch (const std::runtime_error& error) {
            LOG(info) << "Run time error found: " << error.what() << ". DCAFitterN cannot work, skipping the candidate.";
            hCandidatesB->Fill(SVFitting::Fail);
            continue;
          }
          hCandidatesB->Fill(SVFitting::FitOk);
          // JpsiK passed B+ reconstruction

          // calculate relevant properties
          const auto& secondaryVertexBplus = df3.getPCACandidate();
          auto chi2PCA = df3.getChi2AtPCACandidate();
          auto covMatrixPCA = df3.calcPCACovMatrixFlat();
          registry.fill(HIST("hCovSVXX"), covMatrixPCA[0]);
          registry.fill(HIST("hCovPVXX"), covMatrixPV[0]);

          // get JPsi daughters and K tracks (propagated to the B+ vertex if propagateToPCA==true)
          // track.getPxPyPzGlo(pVec) modifies pVec of track
          df3.getTrack(0).getPxPyPzGlo(pVecDauPos);   // momentum of positive Jpsi daughter at the B+ vertex
          df3.getTrack(1).getPxPyPzGlo(pVecDauNeg);   // momentum of negative Jpsi daughter at the B+ vertex
          df3.getTrack(2).getPxPyPzGlo(pVecTrackLf0); // momentum of K at the B+ vertex

          // compute invariant mass square and apply selection
          auto invMass2JpsiK = RecoDecay::m2(std::array{pVecDauPos, pVecDauNeg, pVecTrackLf0}, std::array{o2::constants::physics::MassMuon, o2::constants::physics::MassMuon, o2::constants::physics::MassKPlus});
          if ((invMass2JpsiK < invMass2JpsiHadMin) || (invMass2JpsiK > invMass2JpsiHadMax)) {
            continue;
          }
          registry.fill(HIST("hMassBplusToJpsiK"), std::sqrt(invMass2JpsiK));

          // compute impact parameters of Jpsi and K
          o2::dataformats::DCA dcaDauPos, dcaDauNeg, dcaKaon;
          trackPosParCov.propagateToDCA(primaryVertex, bz, &dcaDauPos);
          trackNegParCov.propagateToDCA(primaryVertex, bz, &dcaDauNeg);
          trackParCovLf0.propagateToDCA(primaryVertex, bz, &dcaKaon);

          // get uncertainty of the decay length
          double phi, theta;
          // getPointDirection modifies phi and theta
          getPointDirection(std::array{collision.posX(), collision.posY(), collision.posZ()}, secondaryVertexBplus, phi, theta);
          auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
          auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));

          // fill the candidate table for the B+ here:
          rowCandidateBpBase(collision.globalIndex(),
                             collision.posX(), collision.posY(), collision.posZ(),
                             secondaryVertexBplus[0], secondaryVertexBplus[1], secondaryVertexBplus[2],
                             errorDecayLength, errorDecayLengthXY,
                             chi2PCA,
                             pVecDauPos[0], pVecDauPos[1], pVecDauPos[2],
                             pVecDauNeg[0], pVecDauNeg[1], pVecDauNeg[2],
                             pVecTrackLf0[0], pVecTrackLf0[1], pVecTrackLf0[2],
                             dcaDauPos.getY(), dcaDauNeg.getY(), dcaKaon.getY(),
                             std::sqrt(dcaDauPos.getSigmaY2()), std::sqrt(dcaDauNeg.getSigmaY2()), std::sqrt(dcaKaon.getSigmaY2()));

          rowCandidateBpProngs(candJpsi.globalIndex(), trackLf0.globalIndex());
        } else if constexpr (DecChannel == DecayChannel::BsToJpsiPhi) {
          for (const auto& trackLf1 : tracksLfDau1ThisCollision) {
            // this track is among daughters
            if (trackLf1.trackId() == candJpsi.prongPosId() || trackLf1.trackId() == candJpsi.prongNegId()) {
              continue;
            }
            auto trackParCovLf1 = getTrackParCov(trackLf1);
            std::array<float, 3> pVecTrackLf1 = trackLf1.pVector();

            // ---------------------------------
            // reconstruct the 4-prong Bs vertex
            hCandidatesB->Fill(SVFitting::BeforeFit);
            try {
              if (df4.process(trackPosParCov, trackNegParCov, trackParCovLf0, trackParCovLf1) == 0) {
                continue;
              }
            } catch (const std::runtime_error& error) {
              LOG(info) << "Run time error found: " << error.what() << ". DCAFitterN cannot work, skipping the candidate.";
              hCandidatesB->Fill(SVFitting::Fail);
              continue;
            }
            hCandidatesB->Fill(SVFitting::FitOk);
            // passed Bs reconstruction

            // calculate relevant properties
            const auto& secondaryVertexBs = df4.getPCACandidate();
            auto chi2PCA = df4.getChi2AtPCACandidate();
            auto covMatrixPCA = df4.calcPCACovMatrixFlat();
            registry.fill(HIST("hCovSVXX"), covMatrixPCA[0]);
            registry.fill(HIST("hCovPVXX"), covMatrixPV[0]);

            // get JPsi daughters and K tracks (propagated to the Bs vertex if propagateToPCA==true)
            // track.getPxPyPzGlo(pVec) modifies pVec of track
            df4.getTrack(0).getPxPyPzGlo(pVecDauPos);   // momentum of Jpsi at the B+ vertex
            df4.getTrack(1).getPxPyPzGlo(pVecDauNeg);   // momentum of Jpsi at the B+ vertex
            df4.getTrack(2).getPxPyPzGlo(pVecTrackLf0); // momentum of K at the B+ vertex
            df4.getTrack(3).getPxPyPzGlo(pVecTrackLf1); // momentum of K at the B+ vertex

            // compute invariant mass square and apply selection
            auto invMass2JpsiPhi = RecoDecay::m2(std::array{pVecDauPos, pVecDauNeg, pVecTrackLf0, pVecTrackLf1}, std::array{o2::constants::physics::MassMuon, o2::constants::physics::MassMuon, o2::constants::physics::MassKPlus, o2::constants::physics::MassKPlus});
            if ((invMass2JpsiPhi < invMass2JpsiHadMin) || (invMass2JpsiPhi > invMass2JpsiHadMax)) {
              continue;
            }
            registry.fill(HIST("hMassBplusToJpsiK"), std::sqrt(invMass2JpsiPhi));

            // compute impact parameters of Jpsi and K
            o2::dataformats::DCA dcaDauPos, dcaDauNeg, dcaTrackLf0, dcaTrackLf1;
            trackPosParCov.propagateToDCA(primaryVertex, bz, &dcaDauPos);
            trackNegParCov.propagateToDCA(primaryVertex, bz, &dcaDauNeg);
            trackParCovLf0.propagateToDCA(primaryVertex, bz, &dcaTrackLf0);
            trackParCovLf1.propagateToDCA(primaryVertex, bz, &dcaTrackLf1);

            // get uncertainty of the decay length
            double phi, theta;
            // getPointDirection modifies phi and theta
            getPointDirection(std::array{collision.posX(), collision.posY(), collision.posZ()}, secondaryVertexBs, phi, theta);
            auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
            auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));

            // fill the candidate table for the Bs here:
            rowCandidateBsBase(collision.globalIndex(),
                               collision.posX(), collision.posY(), collision.posZ(),
                               secondaryVertexBs[0], secondaryVertexBs[1], secondaryVertexBs[2],
                               errorDecayLength, errorDecayLengthXY,
                               chi2PCA,
                               pVecDauPos[0], pVecDauPos[1], pVecDauPos[2],
                               pVecDauNeg[0], pVecDauNeg[1], pVecDauPos[2],
                               pVecTrackLf0[0], pVecTrackLf0[1], pVecTrackLf0[2],
                               pVecTrackLf1[0], pVecTrackLf1[1], pVecTrackLf0[2],
                               dcaDauPos.getY(), dcaDauNeg.getY(), dcaTrackLf0.getY(), dcaTrackLf1.getY(),
                               std::sqrt(dcaDauPos.getSigmaY2()), std::sqrt(dcaDauNeg.getSigmaY2()), std::sqrt(dcaTrackLf0.getSigmaY2()), std::sqrt(dcaTrackLf1.getSigmaY2()));
            rowCandidateBsProngs(candJpsi.globalIndex(), trackLf0.globalIndex(), trackLf1.globalIndex());
          }
        }
      } // TrackLf0 loop
    } // J/Psi loop
  } // end runCandidateCreation

  void processDataBplus(HfRedCollisionsWithExtras const& collisions,
                        soa::Join<aod::HfRedJpsis, aod::HfRedJpsiCov> const& candsJpsi,
                        soa::Join<aod::HfRedBach0Bases, aod::HfRedBach0Cov> const& tracksKaon,
                        aod::HfOrigColCounts const& collisionsCounter,
                        aod::HfCfgBpToJpsi const& configs)
  {
    // Jpsi K invariant-mass window cut
    for (const auto& config : configs) {
      myInvMassWindowJpsiK = config.myInvMassWindowJpsiK();
    }
    // invMassWindowJpsiHadTolerance is used to apply a slightly tighter cut than in JpsiK pair preselection
    // to avoid accepting JpsiK pairs that were not formed in JpsiK pair creator
    double const invMass2JpsiKMin = (massBplus - myInvMassWindowJpsiK + invMassWindowJpsiHadTolerance) * (massBplus - myInvMassWindowJpsiK + invMassWindowJpsiHadTolerance);
    double const invMass2JpsiKMax = (massBplus + myInvMassWindowJpsiK - invMassWindowJpsiHadTolerance) * (massBplus + myInvMassWindowJpsiK - invMassWindowJpsiHadTolerance);

    for (const auto& collisionCounter : collisionsCounter) {
      registry.fill(HIST("hEvents"), 1, collisionCounter.originalCollisionCount());
    }

    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsJpsiThisColl = candsJpsi.sliceBy(candsJpsiPerCollision, thisCollId);
      auto tracksKaonThisCollision = tracksKaon.sliceBy(tracksLf0PerCollision, thisCollId);
      runCandidateCreation<DecayChannel::BplusToJpsiK>(collision, candsJpsiThisColl, tracksKaonThisCollision, tracksKaonThisCollision, invMass2JpsiKMin, invMass2JpsiKMax);
    }
  } // processDataBplus
  PROCESS_SWITCH(HfCandidateCreatorBToJpsiReduced, processDataBplus, "Process data for B+", true);

  void processDataBs(HfRedCollisionsWithExtras const& collisions,
                     soa::Join<aod::HfRedJpsis, aod::HfRedJpsiCov> const& candsJpsi,
                     soa::Join<aod::HfRedBach0Bases, aod::HfRedBach0Cov> const& tracksLfDau0,
                     soa::Join<aod::HfRedBach1Bases, aod::HfRedBach1Cov> const& tracksLfDau1,
                     aod::HfOrigColCounts const& collisionsCounter,
                     aod::HfCfgBsToJpsis const& configs)
  {
    // Jpsi K invariant-mass window cut
    for (const auto& config : configs) {
      myInvMassWindowJpsiPhi = config.myInvMassWindowJpsiPhi();
    }
    // invMassWindowJpsiHadTolerance is used to apply a slightly tighter cut than in JpsiK pair preselection
    // to avoid accepting JpsiK pairs that were not formed in JpsiK pair creator
    double const invMass2JpsiKMin = (massBs - myInvMassWindowJpsiPhi + invMassWindowJpsiHadTolerance) * (massBs - myInvMassWindowJpsiPhi + invMassWindowJpsiHadTolerance);
    double const invMass2JpsiKMax = (massBs + myInvMassWindowJpsiPhi - invMassWindowJpsiHadTolerance) * (massBs + myInvMassWindowJpsiPhi - invMassWindowJpsiHadTolerance);

    for (const auto& collisionCounter : collisionsCounter) {
      registry.fill(HIST("hEvents"), 1, collisionCounter.originalCollisionCount());
    }

    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsJpsiThisColl = candsJpsi.sliceBy(candsJpsiPerCollision, thisCollId);
      auto tracksLf0ThisCollision = tracksLfDau0.sliceBy(tracksLf0PerCollision, thisCollId);
      auto tracksLf1ThisCollision = tracksLfDau1.sliceBy(tracksLf1PerCollision, thisCollId);
      runCandidateCreation<DecayChannel::BsToJpsiPhi>(collision, candsJpsiThisColl, tracksLf0ThisCollision, tracksLf1ThisCollision, invMass2JpsiKMin, invMass2JpsiKMax);
    }
  } // processDataBs
  PROCESS_SWITCH(HfCandidateCreatorBToJpsiReduced, processDataBs, "Process data for Bs", false);
}; // struct

/// Extends the table base with expression columns and performs MC matching.
struct HfCandidateCreatorBToJpsiReducedExpressions {
  Spawns<aod::HfCandBpJPExt> rowCandidateBPlus;
  Spawns<aod::HfCandBsJPExt> rowCandidateBs;
  Produces<aod::HfMcRecRedBps> rowBplusMcRec;
  Produces<aod::HfMcRecRedBss> rowBsMcRec;

  /// Fill candidate information at MC reconstruction level
  /// \param rowsJpsiHadMcRec MC reco information on Jpsi hadron pairs
  /// \param candsBIds prong global indices of B candidates
  template <uint8_t DecChannel, typename McRec, typename CCands>
  void fillMcRec(McRec const& rowsJpsiHadMcRec, CCands const& candsBIds)
  {
    for (const auto& candB : candsBIds) {
      bool filledMcInfo{false};
      if constexpr (DecChannel == DecayChannel::BplusToJpsiK) {
        for (const auto& rowJpsiHadMcRec : rowsJpsiHadMcRec) {
          if ((rowJpsiHadMcRec.jpsiId() != candB.jpsiId()) || (rowJpsiHadMcRec.bachKaId() != candB.bachKaId())) {
            continue;
          }
          rowBplusMcRec(rowJpsiHadMcRec.flagMcMatchRec(), rowJpsiHadMcRec.flagMcDecayChanRec(), rowJpsiHadMcRec.flagWrongCollision(), rowJpsiHadMcRec.debugMcRec(), rowJpsiHadMcRec.ptMother());
          filledMcInfo = true;
          break;
        }
        if (!filledMcInfo) { // protection to get same size tables in case something went wrong: we created a candidate that was not preselected in the Jpsi-K creator
          rowBplusMcRec(0, -1, -1, -1, -1.f);
        }
      } else if constexpr (DecChannel == DecayChannel::BsToJpsiPhi) {
        for (const auto& rowJpsiHadMcRec : rowsJpsiHadMcRec) {
          if ((rowJpsiHadMcRec.jpsiId() != candB.jpsiId()) || (rowJpsiHadMcRec.prong0PhiId() != candB.prong0PhiId()) || (rowJpsiHadMcRec.prong1PhiId() != candB.prong1PhiId())) {
            continue;
          }
          rowBsMcRec(rowJpsiHadMcRec.flagMcMatchRec(), rowJpsiHadMcRec.flagMcDecayChanRec(), rowJpsiHadMcRec.flagWrongCollision(), rowJpsiHadMcRec.debugMcRec(), rowJpsiHadMcRec.ptMother());
          filledMcInfo = true;
          break;
        }
        if (!filledMcInfo) { // protection to get same size tables in case something went wrong: we created a candidate that was not preselected in the Jpsi-K creator
          rowBsMcRec(0, -1, -1, -1, -1.f);
        }
      }
    }
  }

  void processMcBPlus(HfMcRecRedJPKs const& rowsJpsiKMcRec, HfRedBplus2JpsiDaus const& candsBplus)
  {
    fillMcRec<DecayChannel::BplusToJpsiK>(rowsJpsiKMcRec, candsBplus);
  }
  PROCESS_SWITCH(HfCandidateCreatorBToJpsiReducedExpressions, processMcBPlus, "Process MC", false);

  void processMcBs(HfMcRecRedJPPhis const& rowsJpsiPhiMcRec, HfRedBs2JpsiDaus const& bs)
  {
    fillMcRec<DecayChannel::BsToJpsiPhi>(rowsJpsiPhiMcRec, bs);
  }
  PROCESS_SWITCH(HfCandidateCreatorBToJpsiReducedExpressions, processMcBs, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCandidateCreatorBToJpsiReduced>(cfgc),
                      adaptAnalysisTask<HfCandidateCreatorBToJpsiReducedExpressions>(cfgc)};
}
