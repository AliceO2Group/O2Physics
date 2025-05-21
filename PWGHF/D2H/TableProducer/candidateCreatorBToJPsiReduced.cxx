
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

/// \file candidateCreatorBToJPsiReduced.cxx
/// \brief Reconstruction of B->J/Psi hadron candidates
///
/// \author Fabrizio Chinu <fabrizio.chinu@cern.ch>, Università degli Studi and INFN Torino
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

#include <memory>

#include "CommonConstants/PhysicsConstants.h"
#include "DCAFitter/DCAFitterN.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/V0.h"

#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/D2H/DataModel/ReducedDataModel.h"
#include "PWGHF/Utils/utilsTrkCandHf.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_trkcandsel;

/// Reconstruction of B+ candidates
struct HfCandidateCreatorBToJPsiReduced {
  Produces<aod::HfCandBplusBase> rowCandidateBase;       // table defined in CandidateReconstructionTables.h
  Produces<aod::HfRedBplus2JPsiDaus> rowCandidateProngs; // table defined in ReducedDataModel.h

  // vertexing
  Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCA{"useAbsDCA", true, "Minimise abs. distance rather than chi2"};
  Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any B+ is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations is chi2/chi2old > this"};

  Configurable<bool> runJPsiToee{"runJPsiToee", false, "Run analysis for J/Psi to ee (debug)"};
  // selection
  Configurable<double> invMassWindowJPsiKTolerance{"invMassWindowJPsiKTolerance", 0.01, "invariant-mass window tolerance for J/Psi K pair preselections (GeV/c2)"};

  float myInvMassWindowJPsiK{1.}; // variable that will store the value of invMassWindowJPsiK (defined in dataCreatorJPsiKReduced.cxx)
  double massBplus{0.};
  double bz{0.};
  o2::vertexing::DCAFitterN<2> df2; // fitter for B vertex (2-prong vertex fitter)

  using HfRedCollisionsWithExtras = soa::Join<aod::HfRedCollisions, aod::HfRedCollExtras>;

  Preslice<soa::Join<aod::HfRedJPsis, aod::HfRedJPsiCov>> candsJPsiPerCollision = hf_track_index_reduced::hfRedCollisionId;
  Preslice<soa::Join<aod::HfRedBachProng0Bases, aod::HfRedBachProng0Cov>> tracksKaonPerCollision = hf_track_index_reduced::hfRedCollisionId;

  std::shared_ptr<TH1> hCandidates;
  o2::base::Propagator::MatCorrType noMatCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    // invariant-mass window cut
    massBplus = o2::constants::physics::MassBPlus;

    // Initialize fitter
    df2.setPropagateToPCA(propagateToPCA);
    df2.setMaxR(maxR);
    df2.setMaxDZIni(maxDZIni);
    df2.setMinParamChange(minParamChange);
    df2.setMinRelChi2Change(minRelChi2Change);
    df2.setUseAbsDCA(useAbsDCA);
    df2.setWeightedFinalPCA(useWeightedFinalPCA);
    df2.setMatCorrType(noMatCorr);

    // histograms
    registry.add("hMassJPsi", "J/Psi mass;#it{M}_{#mu#mu} (GeV/#it{c}^{2});Counts", {HistType::kTH1F, {{600, 2.5, 3.7, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassBplusToJPsiK", "2-prong candidates;inv. mass (B^{+} #rightarrow #overline{D^{0}}#pi^{#plus} #rightarrow #pi^{#minus}K^{#plus}#pi^{#plus}) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 3., 8.}}});
    registry.add("hCovPVXX", "2-prong candidates;XX element of cov. matrix of prim. vtx. position (cm^{2});entries", {HistType::kTH1F, {{100, 0., 1.e-4}}});
    registry.add("hCovSVXX", "2-prong candidates;XX element of cov. matrix of sec. vtx. position (cm^{2});entries", {HistType::kTH1F, {{100, 0., 0.2}}});
    registry.add("hEvents", "Events;;entries", HistType::kTH1F, {{1, 0.5, 1.5}});

    /// candidate monitoring
    hCandidates = registry.add<TH1>("hFitCandidatesJPsi", "candidates counter", {HistType::kTH1D, {axisCands}});
    setLabelHistoCands(hCandidates);
  }

  /// Main function to perform B+ candidate creation
  /// \param collision the collision
  /// \param candsJPsiThisColl J/Psi candidates in this collision
  /// \param tracksKaonThisCollision kaon tracks in this collision
  /// \param invMass2JPsiKMin minimum B+ invariant-mass
  /// \param invMass2JPsiKMax maximum B+ invariant-mass
  template <typename Cands, typename Kaons, typename Coll>
  void runCandidateCreation(Coll const& collision,
                            Cands const& candsJPsiThisColl,
                            Kaons const& tracksKaonThisCollision,
                            const float& invMass2JPsiKMin,
                            const float& invMass2JPsiKMax)
  {
    auto primaryVertex = getPrimaryVertex(collision);
    auto covMatrixPV = primaryVertex.getCov();

    // Set the magnetic field from ccdb
    bz = collision.bz();
    df2.setBz(bz);

    for (const auto& candJPsi : candsJPsiThisColl) {
      o2::track::TrackParametrizationWithError<float> trackPosParCov(
        candJPsi.xDauPos(), candJPsi.alphaDauPos(), {candJPsi.yDauPos(), candJPsi.zDauPos(), candJPsi.snpDauPos(), candJPsi.tglDauPos(), candJPsi.signed1PtDauPos()}, 1 /*Charge*/, 1 /*Muon*/);
      o2::track::TrackParametrizationWithError<float> trackNegParCov(
        candJPsi.xDauNeg(), candJPsi.alphaDauNeg(), {candJPsi.yDauNeg(), candJPsi.zDauNeg(), candJPsi.snpDauNeg(), candJPsi.tglDauNeg(), candJPsi.signed1PtDauNeg()}, -1 /*Charge*/, 1 /*Muon*/);

      // ---------------------------------
      // reconstruct J/Psi candidate
      o2::track::TrackParCov trackParCovJPsi{};
      std::array<float, 3> pVecJPsi{};
      registry.fill(HIST("hFitCandidatesJPsi"), SVFitting::BeforeFit);
      try {
        if (df2.process(trackPosParCov, trackNegParCov) == 0) {
          LOG(info) << "DCAFitterN failed to reconstruct J/Psi candidate, skipping the candidate.";
          continue;
        }
      } catch (const std::runtime_error& error) {
        LOG(info) << "Run time error found: " << error.what() << ". DCAFitterN cannot work, skipping the candidate.";
        registry.fill(HIST("hFitCandidatesJPsi"), SVFitting::Fail);
        continue;
      }
      registry.fill(HIST("hFitCandidatesJPsi"), SVFitting::FitOk);

      std::array<float, 3> pVecDauPos{candJPsi.pxDauPos(), candJPsi.pyDauPos(), candJPsi.pzDauPos()};
      std::array<float, 3> pVecDauNeg{candJPsi.pxDauNeg(), candJPsi.pyDauNeg(), candJPsi.pzDauNeg()};
      LOG(info) << "pVecDauPos before: " << pVecDauPos[0] << " " << pVecDauPos[1] << " " << pVecDauPos[2];

      // df2.getTrack(0).getPxPyPzGlo(pVecDauPos);
      // df2.getTrack(1).getPxPyPzGlo(pVecDauNeg);
      LOG(info) << "pVecDauPos: " << pVecDauPos[0] << " " << pVecDauPos[1] << " " << pVecDauPos[2];
      LOG(info) << "pVecDauPos from table: " << candJPsi.pxDauPos() << " " << candJPsi.pyDauPos() << " " << candJPsi.pzDauPos();
      pVecJPsi = RecoDecay::pVec(pVecDauPos, pVecDauNeg);
      trackParCovJPsi = df2.createParentTrackParCov();
      trackParCovJPsi.setAbsCharge(0); // to be sure

      float invMassJPsi{0.f};
      if (runJPsiToee) {
        invMassJPsi = RecoDecay::m2(std::array{pVecDauPos, pVecDauNeg}, std::array{o2::constants::physics::MassElectron, o2::constants::physics::MassElectron});
      } else {
        invMassJPsi = RecoDecay::m2(std::array{pVecDauPos, pVecDauNeg}, std::array{o2::constants::physics::MassMuon, o2::constants::physics::MassMuon});
      }
      invMassJPsi = std::sqrt(invMassJPsi);
      registry.fill(HIST("hMassJPsi"), invMassJPsi);

      for (const auto& trackKaon : tracksKaonThisCollision) {
        // this track is among daughters
        if (trackKaon.trackId() == candJPsi.prongPosId() || trackKaon.trackId() == candJPsi.prongNegId()) {
          continue;
        }
        auto trackParCovKa = getTrackParCov(trackKaon);
        std::array<float, 3> pVecKaon = trackKaon.pVector();

        // compute invariant mass square and apply selection
        auto invMass2JPsiK = RecoDecay::m2(std::array{pVecJPsi, pVecKaon}, std::array{o2::constants::physics::MassJPsi, o2::constants::physics::MassKPlus});
        if ((invMass2JPsiK < invMass2JPsiKMin) || (invMass2JPsiK > invMass2JPsiKMax)) {
          continue;
        }
        // ---------------------------------
        // reconstruct the 2-prong B+ vertex
        hCandidates->Fill(SVFitting::BeforeFit);
        try {
          if (df2.process(trackParCovJPsi, trackParCovKa) == 0) {
            continue;
          }
        } catch (const std::runtime_error& error) {
          LOG(info) << "Run time error found: " << error.what() << ". DCAFitterN cannot work, skipping the candidate.";
          hCandidates->Fill(SVFitting::Fail);
          continue;
        }
        hCandidates->Fill(SVFitting::FitOk);
        // JPsiK passed B+ reconstruction

        // calculate relevant properties
        const auto& secondaryVertexBplus = df2.getPCACandidate();
        auto chi2PCA = df2.getChi2AtPCACandidate();
        auto covMatrixPCA = df2.calcPCACovMatrixFlat();
        registry.fill(HIST("hCovSVXX"), covMatrixPCA[0]);
        registry.fill(HIST("hCovPVXX"), covMatrixPV[0]);

        // propagate JPsi and K to the B+ vertex
        df2.propagateTracksToVertex();
        // track.getPxPyPzGlo(pVec) modifies pVec of track
        df2.getTrack(0).getPxPyPzGlo(pVecJPsi); // momentum of JPsi at the B+ vertex
        df2.getTrack(1).getPxPyPzGlo(pVecKaon); // momentum of K at the B+ vertex

        registry.fill(HIST("hMassBplusToJPsiK"), std::sqrt(invMass2JPsiK));

        // compute impact parameters of JPsi and K
        o2::dataformats::DCA dcaJPsi;
        o2::dataformats::DCA dcaKaon;
        trackParCovJPsi.propagateToDCA(primaryVertex, bz, &dcaJPsi);
        trackParCovKa.propagateToDCA(primaryVertex, bz, &dcaKaon);

        // get uncertainty of the decay length
        double phi, theta;
        // getPointDirection modifies phi and theta
        getPointDirection(std::array{collision.posX(), collision.posY(), collision.posZ()}, secondaryVertexBplus, phi, theta);
        auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
        auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));

        // fill the candidate table for the B+ here:
        rowCandidateBase(collision.globalIndex(),
                         collision.posX(), collision.posY(), collision.posZ(),
                         secondaryVertexBplus[0], secondaryVertexBplus[1], secondaryVertexBplus[2],
                         errorDecayLength, errorDecayLengthXY,
                         chi2PCA,
                         pVecJPsi[0], pVecJPsi[1], pVecJPsi[2],
                         pVecKaon[0], pVecKaon[1], pVecKaon[2],
                         dcaJPsi.getY(), dcaKaon.getY(),
                         std::sqrt(dcaJPsi.getSigmaY2()), std::sqrt(dcaKaon.getSigmaY2()));

        rowCandidateProngs(candJPsi.globalIndex(), trackKaon.globalIndex());

      } // kaon loop
    } // J/Psi loop
  } // end runCandidateCreation

  void processData(HfRedCollisionsWithExtras const& collisions,
                   soa::Join<aod::HfRedJPsis, aod::HfRedJPsiCov> const& candsJPsi,
                   soa::Join<aod::HfRedBachProng0Bases, aod::HfRedBachProng0Cov> const& tracksKaon,
                   aod::HfOrigColCounts const& collisionsCounter,
                   aod::HfCandBpToJPsiKConfigs const& configs)
  {
    // JPsi K invariant-mass window cut
    for (const auto& config : configs) {
      myInvMassWindowJPsiK = config.myInvMassWindowJPsiK();
    }
    // invMassWindowJPsiKTolerance is used to apply a slightly tighter cut than in JPsiK pair preselection
    // to avoid accepting JPsiK pairs that were not formed in JPsiK pair creator
    double invMass2JPsiKMin = (massBplus - myInvMassWindowJPsiK + invMassWindowJPsiKTolerance) * (massBplus - myInvMassWindowJPsiK + invMassWindowJPsiKTolerance);
    double invMass2JPsiKMax = (massBplus + myInvMassWindowJPsiK - invMassWindowJPsiKTolerance) * (massBplus + myInvMassWindowJPsiK - invMassWindowJPsiKTolerance);

    for (const auto& collisionCounter : collisionsCounter) {
      registry.fill(HIST("hEvents"), 1, collisionCounter.originalCollisionCount());
    }

    static int ncol = 0;

    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsJPsiThisColl = candsJPsi.sliceBy(candsJPsiPerCollision, thisCollId);
      auto tracksKaonThisCollision = tracksKaon.sliceBy(tracksKaonPerCollision, thisCollId);
      runCandidateCreation(collision, candsJPsiThisColl, tracksKaonThisCollision, invMass2JPsiKMin, invMass2JPsiKMax);
      if (ncol % 10000 == 0) {
        LOG(debug) << ncol << " collisions parsed";
      }
      ncol++;
    }
  } // processData

  PROCESS_SWITCH(HfCandidateCreatorBToJPsiReduced, processData, "Process data", true);
}; // struct

/// Extends the table base with expression columns and performs MC matching.
struct HfCandidateCreatorBplusToJPsiReducedExpressions {
  Spawns<aod::HfCandBplusExt> rowCandidateBPlus;
  Spawns<aod::HfRedBachProng0Ext> rowTracksExt;
  Produces<aod::HfMcRecRedBps> rowBplusMcRec;

  /// Fill candidate information at MC reconstruction level
  /// \param rowsJPsiKMcRec MC reco information on JPsiK pairs
  /// \param candsBplus prong global indices of B+ candidates
  template <typename McRec>
  void fillBplusMcRec(McRec const& rowsJPsiKMcRec, HfRedBplus2JPsiDaus const& candsBplus)
  {
    for (const auto& candBplus : candsBplus) {
      bool filledMcInfo{false};
      for (const auto& rowJPsiKMcRec : rowsJPsiKMcRec) {
        if ((rowJPsiKMcRec.prong0Id() != candBplus.jPsiId()) || (rowJPsiKMcRec.prong1Id() != candBplus.bachKaId())) {
          continue;
        }
        rowBplusMcRec(rowJPsiKMcRec.flagMcMatchRec(), rowJPsiKMcRec.flagWrongCollision(), rowJPsiKMcRec.debugMcRec(), rowJPsiKMcRec.ptMother());
        filledMcInfo = true;
        break;
      }
      if (!filledMcInfo) { // protection to get same size tables in case something went wrong: we created a candidate that was not preselected in the JPsi-K creator
        rowBplusMcRec(0, -1, -1, -1.f);
      }
    }
  }

  void processMc(HfMcRecRedD0Pis const& rowsJPsiKMcRec, HfRedBplus2JPsiDaus const& candsBplus)
  {
    fillBplusMcRec(rowsJPsiKMcRec, candsBplus);
  }
  PROCESS_SWITCH(HfCandidateCreatorBplusToJPsiReducedExpressions, processMc, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCandidateCreatorBToJPsiReduced>(cfgc),
                      adaptAnalysisTask<HfCandidateCreatorBplusToJPsiReducedExpressions>(cfgc)};
}
