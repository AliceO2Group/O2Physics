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

/// \file candidateCreatorBplusReduced.cxx
/// \brief Reconstruction of B+ candidates
///
/// \author Antonio Palasciano <antonio.palasciano@cern.ch>, Università degli Studi di Bari

#include "PWGHF/D2H/DataModel/ReducedDataModel.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/Utils/utilsTrkCandHf.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"

#include <CommonConstants/PhysicsConstants.h>
#include <DCAFitter/DCAFitterN.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/DCA.h>

#include <TH1.h>

#include <array>
#include <cmath>
#include <memory>
#include <numeric>
#include <stdexcept>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_trkcandsel;

/// Reconstruction of B+ candidates
struct HfCandidateCreatorBplusReduced {
  Produces<aod::HfCandBplusBase> rowCandidateBase;         // table defined in CandidateReconstructionTables.h
  Produces<aod::HfRedBplusProngs> rowCandidateProngs;      // table defined in ReducedDataModel.h
  Produces<aod::HfRedBplusD0Mls> rowCandidateDmesMlScores; // table defined in ReducedDataModel.h

  // vertexing
  Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCA{"useAbsDCA", true, "Minimise abs. distance rather than chi2"};
  Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any B+ is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations is chi2/chi2old > this"};
  // selection
  Configurable<double> invMassWindowD0PiTolerance{"invMassWindowD0PiTolerance", 0.01, "invariant-mass window tolerance for D0Pi pair preselections (GeV/c2)"};

  float myInvMassWindowD0Pi{1.}; // variable that will store the value of invMassWindowD0Pi (defined in dataCreatorD0PiReduced.cxx)
  double massPi{0.};
  double massD0{0.};
  double massBplus{0.};
  double bz{0.};
  o2::vertexing::DCAFitterN<2> df2; // fitter for B vertex (2-prong vertex fitter)

  using HfRedCollisionsWithExtras = soa::Join<aod::HfRedCollisions, aod::HfRedCollExtras>;

  Preslice<soa::Join<aod::HfRed2Prongs, aod::HfRed2ProngsCov>> candsDPerCollision = hf_track_index_reduced::hfRedCollisionId;
  Preslice<soa::Join<aod::HfRed2Prongs, aod::HfRed2ProngsCov, aod::HfRed2ProngsMl>> candsDWithMlPerCollision = hf_track_index_reduced::hfRedCollisionId;
  Preslice<soa::Join<aod::HfRedTracks, aod::HfRedTracksCov>> tracksPionPerCollision = hf_track_index_reduced::hfRedCollisionId;

  std::shared_ptr<TH1> hCandidates;
  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    std::array<bool, 2> doprocess{doprocessData, doprocessDataWithDmesMl};
    if ((std::accumulate(doprocess.begin(), doprocess.end(), 0)) != 1) {
      LOGP(fatal, "Only one process function for data should be enabled at a time.");
    }

    // invariant-mass window cut
    massPi = o2::constants::physics::MassPiPlus;
    massD0 = o2::constants::physics::MassD0;
    massBplus = o2::constants::physics::MassBPlus;

    // Initialize fitter
    df2.setPropagateToPCA(propagateToPCA);
    df2.setMaxR(maxR);
    df2.setMaxDZIni(maxDZIni);
    df2.setMinParamChange(minParamChange);
    df2.setMinRelChi2Change(minRelChi2Change);
    df2.setUseAbsDCA(useAbsDCA);
    df2.setWeightedFinalPCA(useWeightedFinalPCA);

    // histograms
    registry.add("hMassBplusToD0Pi", "2-prong candidates;inv. mass (B^{+} #rightarrow #overline{D^{0}}#pi^{#plus} #rightarrow #pi^{#minus}K^{#plus}#pi^{#plus}) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 3., 8.}}});
    registry.add("hCovPVXX", "2-prong candidates;XX element of cov. matrix of prim. vtx. position (cm^{2});entries", {HistType::kTH1F, {{100, 0., 1.e-4}}});
    registry.add("hCovSVXX", "2-prong candidates;XX element of cov. matrix of sec. vtx. position (cm^{2});entries", {HistType::kTH1F, {{100, 0., 0.2}}});
    registry.add("hEvents", "Events;;entries", HistType::kTH1F, {{1, 0.5, 1.5}});

    /// candidate monitoring
    hCandidates = registry.add<TH1>("hCandidates", "candidates counter", {HistType::kTH1D, {axisCands}});
    setLabelHistoCands(hCandidates);
  }

  /// Main function to perform B+ candidate creation
  /// \param withDmesMl is the flag to use the table with ML scores for the D0 daughter (only possible if present in the derived data)
  /// \param collision the collision
  /// \param candsDThisColl B+ candidates in this collision
  /// \param tracksPionThisCollision pion tracks in this collision
  /// \param invMass2D0PiMin minimum B+ invariant-mass
  /// \param invMass2D0PiMax maximum B+ invariant-mass
  template <bool WithDmesMl, typename Cands, typename Pions, typename Coll>
  void runCandidateCreation(Coll const& collision,
                            Cands const& candsDThisColl,
                            Pions const& tracksPionThisCollision,
                            const float& invMass2D0PiMin,
                            const float& invMass2D0PiMax)
  {
    auto primaryVertex = getPrimaryVertex(collision);
    auto covMatrixPV = primaryVertex.getCov();

    // Set the magnetic field from ccdb
    bz = collision.bz();
    df2.setBz(bz);

    for (const auto& candD0 : candsDThisColl) {
      auto trackParCovD = getTrackParCov(candD0);
      std::array<float, 3> pVecD0 = candD0.pVector();

      for (const auto& trackPion : tracksPionThisCollision) {
        auto trackParCovPi = getTrackParCov(trackPion);
        std::array<float, 3> pVecPion = trackPion.pVector();

        // compute invariant mass square and apply selection
        auto invMass2D0Pi = RecoDecay::m2(std::array{pVecD0, pVecPion}, std::array{massD0, massPi});
        if ((invMass2D0Pi < invMass2D0PiMin) || (invMass2D0Pi > invMass2D0PiMax)) {
          continue;
        }
        // ---------------------------------
        // reconstruct the 2-prong B+ vertex
        hCandidates->Fill(SVFitting::BeforeFit);
        try {
          if (df2.process(trackParCovD, trackParCovPi) == 0) {
            continue;
          }
        } catch (const std::runtime_error& error) {
          LOG(info) << "Run time error found: " << error.what() << ". DCAFitterN cannot work, skipping the candidate.";
          hCandidates->Fill(SVFitting::Fail);
          continue;
        }
        hCandidates->Fill(SVFitting::FitOk);
        // D0Pi passed B+ reconstruction

        // calculate relevant properties
        const auto& secondaryVertexBplus = df2.getPCACandidate();
        auto chi2PCA = df2.getChi2AtPCACandidate();
        auto covMatrixPCA = df2.calcPCACovMatrixFlat();
        registry.fill(HIST("hCovSVXX"), covMatrixPCA[0]);
        registry.fill(HIST("hCovPVXX"), covMatrixPV[0]);

        // propagate D0 and Pi to the B+ vertex
        df2.propagateTracksToVertex();
        // track.getPxPyPzGlo(pVec) modifies pVec of track
        df2.getTrack(0).getPxPyPzGlo(pVecD0);   // momentum of D0 at the B+ vertex
        df2.getTrack(1).getPxPyPzGlo(pVecPion); // momentum of Pi at the B+ vertex

        registry.fill(HIST("hMassBplusToD0Pi"), std::sqrt(invMass2D0Pi));

        // compute impact parameters of D0 and Pi
        o2::dataformats::DCA dcaD0;
        o2::dataformats::DCA dcaPion;
        trackParCovD.propagateToDCA(primaryVertex, bz, &dcaD0);
        trackParCovPi.propagateToDCA(primaryVertex, bz, &dcaPion);

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
                         pVecD0[0], pVecD0[1], pVecD0[2],
                         pVecPion[0], pVecPion[1], pVecPion[2],
                         dcaD0.getY(), dcaPion.getY(),
                         std::sqrt(dcaD0.getSigmaY2()), std::sqrt(dcaPion.getSigmaY2()));

        rowCandidateProngs(candD0.globalIndex(), trackPion.globalIndex());

        if constexpr (WithDmesMl) {
          if (trackPion.signed1Pt() < 0) {
            rowCandidateDmesMlScores(candD0.mlScoreBkgMassHypo0(), candD0.mlScorePromptMassHypo0(), candD0.mlScoreNonpromptMassHypo0());
          } else {
            rowCandidateDmesMlScores(candD0.mlScoreBkgMassHypo1(), candD0.mlScorePromptMassHypo1(), candD0.mlScoreNonpromptMassHypo1());
          }
        }
      } // pi loop
    } // D0 loop
  } // end runCandidateCreation

  void processData(HfRedCollisionsWithExtras const& collisions,
                   soa::Join<aod::HfRed2Prongs, aod::HfRed2ProngsCov> const& candsD,
                   soa::Join<aod::HfRedTrackBases, aod::HfRedTracksCov> const& tracksPion,
                   aod::HfOrigColCounts const& collisionsCounter,
                   aod::HfCandBpConfigs const& configs)
  {
    // D0Pi invariant-mass window cut
    for (const auto& config : configs) {
      myInvMassWindowD0Pi = config.myInvMassWindowD0Pi();
    }
    // invMassWindowD0PiTolerance is used to apply a slightly tighter cut than in D0Pi pair preselection
    // to avoid accepting D0Pi pairs that were not formed in D0Pi pair creator
    double const invMass2D0PiMin = (massBplus - myInvMassWindowD0Pi + invMassWindowD0PiTolerance) * (massBplus - myInvMassWindowD0Pi + invMassWindowD0PiTolerance);
    double const invMass2D0PiMax = (massBplus + myInvMassWindowD0Pi - invMassWindowD0PiTolerance) * (massBplus + myInvMassWindowD0Pi - invMassWindowD0PiTolerance);

    for (const auto& collisionCounter : collisionsCounter) {
      registry.fill(HIST("hEvents"), 1, collisionCounter.originalCollisionCount());
    }

    static int ncol = 0;

    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsD.sliceBy(candsDPerCollision, thisCollId);
      auto tracksPionThisCollision = tracksPion.sliceBy(tracksPionPerCollision, thisCollId);
      runCandidateCreation<false>(collision, candsDThisColl, tracksPionThisCollision, invMass2D0PiMin, invMass2D0PiMax);
      if (ncol % 10000 == 0) {
        LOG(debug) << ncol << " collisions parsed";
      }
      ncol++;
    }
  } // processData

  PROCESS_SWITCH(HfCandidateCreatorBplusReduced, processData, "Process data without any ML score", true);

  void processDataWithDmesMl(HfRedCollisionsWithExtras const& collisions,
                             soa::Join<aod::HfRed2Prongs, aod::HfRed2ProngsCov, aod::HfRed2ProngsMl> const& candsD,
                             soa::Join<aod::HfRedTrackBases, aod::HfRedTracksCov> const& tracksPion,
                             aod::HfOrigColCounts const& collisionsCounter,
                             aod::HfCandBpConfigs const& configs)
  {
    // D0Pi invariant-mass window cut
    for (const auto& config : configs) {
      myInvMassWindowD0Pi = config.myInvMassWindowD0Pi();
    }
    // invMassWindowD0PiTolerance is used to apply a slightly tighter cut than in D0Pi pair preselection
    // to avoid accepting D0Pi pairs that were not formed in D0Pi pair creator
    float const invMass2D0PiMin = (massBplus - myInvMassWindowD0Pi + invMassWindowD0PiTolerance) * (massBplus - myInvMassWindowD0Pi + invMassWindowD0PiTolerance);
    float const invMass2D0PiMax = (massBplus + myInvMassWindowD0Pi - invMassWindowD0PiTolerance) * (massBplus + myInvMassWindowD0Pi - invMassWindowD0PiTolerance);

    for (const auto& collisionCounter : collisionsCounter) {
      registry.fill(HIST("hEvents"), 1, collisionCounter.originalCollisionCount());
    }

    static int ncol = 0;
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsD.sliceBy(candsDPerCollision, thisCollId);
      auto tracksPionThisCollision = tracksPion.sliceBy(tracksPionPerCollision, thisCollId);
      runCandidateCreation<true>(collision, candsDThisColl, tracksPionThisCollision, invMass2D0PiMin, invMass2D0PiMax);
      if (ncol % 10000 == 0) {
        LOGP(debug, "collisions parsed {}", ncol);
      }
      ncol++;
    }
  } // processDataWithDmesMl

  PROCESS_SWITCH(HfCandidateCreatorBplusReduced, processDataWithDmesMl, "Process data with ML scores of D mesons", false);
}; // struct

/// Extends the table base with expression columns and performs MC matching.
struct HfCandidateCreatorBplusReducedExpressions {
  Spawns<aod::HfCandBplusExt> rowCandidateBPlus;
  Spawns<aod::HfRedTracksExt> rowTracksExt;
  Produces<aod::HfMcRecRedBps> rowBplusMcRec;
  Produces<aod::HfMcCheckBps> rowBplusMcCheck;

  /// Fill candidate information at MC reconstruction level
  /// \param checkDecayTypeMc
  /// \param rowsD0PiMcRec MC reco information on D0Pi pairs
  /// \param candsBplus prong global indices of B+ candidates
  template <bool CheckDecayTypeMc, typename McRec>
  void fillBplusMcRec(McRec const& rowsD0PiMcRec, HfRedBplusProngs const& candsBplus)
  {
    for (const auto& candBplus : candsBplus) {
      bool filledMcInfo{false};
      for (const auto& rowD0PiMcRec : rowsD0PiMcRec) {
        if ((rowD0PiMcRec.prong0Id() != candBplus.prong0Id()) || (rowD0PiMcRec.prong1Id() != candBplus.prong1Id())) {
          continue;
        }
        rowBplusMcRec(rowD0PiMcRec.flagMcMatchRec(), -1 /*channel*/, rowD0PiMcRec.flagWrongCollision(), rowD0PiMcRec.debugMcRec(), rowD0PiMcRec.ptMother());
        filledMcInfo = true;
        if constexpr (CheckDecayTypeMc) {
          rowBplusMcCheck(rowD0PiMcRec.pdgCodeBeautyMother(),
                          rowD0PiMcRec.pdgCodeCharmMother(),
                          rowD0PiMcRec.pdgCodeProng0(),
                          rowD0PiMcRec.pdgCodeProng1(),
                          rowD0PiMcRec.pdgCodeProng2());
        }
        break;
      }
      if (!filledMcInfo) { // protection to get same size tables in case something went wrong: we created a candidate that was not preselected in the D0-Pi creator
        rowBplusMcRec(0, -1, -1, -1, -1.f);
        if constexpr (CheckDecayTypeMc) {
          rowBplusMcCheck(-1, -1, -1, -1, -1);
        }
      }
    }
  }

  void processMc(HfMcRecRedD0Pis const& rowsD0PiMcRec, HfRedBplusProngs const& candsBplus)
  {
    fillBplusMcRec<false>(rowsD0PiMcRec, candsBplus);
  }
  PROCESS_SWITCH(HfCandidateCreatorBplusReducedExpressions, processMc, "Process MC", false);

  void processMcWithDecayTypeCheck(soa::Join<HfMcRecRedD0Pis, HfMcCheckD0Pis> const& rowsD0PiMcRec, HfRedBplusProngs const& candsBplus)
  {
    fillBplusMcRec<true>(rowsD0PiMcRec, candsBplus);
  }
  PROCESS_SWITCH(HfCandidateCreatorBplusReducedExpressions, processMcWithDecayTypeCheck, "Process MC with decay type checks", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCandidateCreatorBplusReduced>(cfgc),
                      adaptAnalysisTask<HfCandidateCreatorBplusReducedExpressions>(cfgc)};
}
