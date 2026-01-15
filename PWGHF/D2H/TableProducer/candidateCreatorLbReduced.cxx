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

/// \file candidateCreatorLbReduced.cxx
/// \brief Reconstruction of Lb candidates
///
/// \author Biao Zhang <biao.zhang@cern.ch>, Heidelberg University

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
#include <utility>

using namespace o2;
using namespace o2::aod;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_trkcandsel;

/// Reconstruction of Lb candidates
struct HfCandidateCreatorLbReduced {
  Produces<aod::HfCandLbBase> rowCandidateBase;       // table defined in CandidateReconstructionTables.h
  Produces<aod::HfRedLbProngs> rowCandidateProngs;    // table defined in ReducedDataModel.h
  Produces<aod::HfRedLbLcMls> rowCandidateLcMlScores; // table defined in ReducedDataModel.h

  // vertexing
  Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCA{"useAbsDCA", true, "Minimise abs. distance rather than chi2"};
  Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  Configurable<float> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<float> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<float> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any Lb is smaller than this"};
  Configurable<float> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations is chi2/chi2old > this"};
  // selection
  Configurable<float> invMassWindowLcPiTolerance{"invMassWindowLcPiTolerance", 0.01, "invariant-mass window tolerance for LcPi pair preselections (GeV/c2)"};

  float myInvMassWindowLcPi{1.}; // variable that will store the value of invMassWindowLcPi
  float bz{0.};

  o2::vertexing::DCAFitterN<2> df2; // fitter for B vertex (2-prong vertex fitter)

  using HfRedCollisionsWithExtras = soa::Join<aod::HfRedCollisions, aod::HfRedCollExtras>;

  Preslice<soa::Join<aod::HfRed3Prongs, aod::HfRed3ProngsCov>> candsLcPerCollision = hf_track_index_reduced::hfRedCollisionId;
  Preslice<soa::Join<aod::HfRed3Prongs, aod::HfRed3ProngsCov, aod::HfRed3ProngsMl>> candsDWithMlPerCollision = hf_track_index_reduced::hfRedCollisionId;
  Preslice<soa::Join<aod::HfRedTrackBases, aod::HfRedTracksCov>> tracksPionPerCollision = hf_track_index_reduced::hfRedCollisionId;

  std::shared_ptr<TH1> hCandidates;
  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    std::array<bool, 2> doprocess{doprocessData, doprocessDataWithLcMl};
    if ((std::accumulate(doprocess.begin(), doprocess.end(), 0)) != 1) {
      LOGP(fatal, "Only one process function for data should be enabled at a time.");
    }

    // Initialize fitter
    df2.setPropagateToPCA(propagateToPCA);
    df2.setMaxR(maxR);
    df2.setMaxDZIni(maxDZIni);
    df2.setMinParamChange(minParamChange);
    df2.setMinRelChi2Change(minRelChi2Change);
    df2.setUseAbsDCA(useAbsDCA);
    df2.setWeightedFinalPCA(useWeightedFinalPCA);

    // histograms
    registry.add("hMassLambdaB0ToLcPi", "2-prong candidates;inv. mass (#Lambda_{b}^{0} #rightarrow #Lambda_{c}^{#plus}#pi^{#minus} #rightarrow pK^{#minus}#pi^{#plus}#pi^{#minus}) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 3., 8.}}});
    registry.add("hCovPVXX", "2-prong candidates;XX element of cov. matrix of prim. vtx. position (cm^{2});entries", {HistType::kTH1F, {{100, 0., 1.e-4}}});
    registry.add("hCovSVXX", "2-prong candidates;XX element of cov. matrix of sec. vtx. position (cm^{2});entries", {HistType::kTH1F, {{100, 0., 0.2}}});
    registry.add("hEvents", "Events;;entries", HistType::kTH1F, {{1, 0.5, 1.5}});

    /// candidate monitoring
    hCandidates = registry.add<TH1>("hCandidates", "candidates counter", {HistType::kTH1D, {axisCands}});
    setLabelHistoCands(hCandidates);
  }

  template <typename Config>
  std::pair<float, float> computeInvMass2LcPiWindow(Config const& configs,
                                                    float invMassWindowLcPiTolerance)
  {

    myInvMassWindowLcPi = 0.0f;
    for (const auto& config : configs) {
      myInvMassWindowLcPi = config.myInvMassWindowLcPi();
    }

    float const deltaMin = MassLambdaB0 - myInvMassWindowLcPi + invMassWindowLcPiTolerance;
    float const deltaMax = MassLambdaB0 + myInvMassWindowLcPi - invMassWindowLcPiTolerance;

    float const invMass2LcPiMin = deltaMin * deltaMin;
    float const invMass2LcPiMax = deltaMax * deltaMax;

    return {invMass2LcPiMin, invMass2LcPiMax};
  }

  /// Main function to perform Lb candidate creation
  /// \param withLcMl is the flag to use the table with ML scores for the Lc daughter (only possible if present in the derived data)
  /// \param collision the collision
  /// \param candsLcThisColl Lc candidates in this collision
  /// \param tracksPionThisCollision pion tracks in this collision
  /// \param invMass2LcPiMin minimum Lb invariant-mass
  /// \param invMass2LcPiMax maximum Lb invariant-mass
  template <bool WithLcMl, typename Cands, typename Pions, typename Coll>
  void runCandidateCreation(Coll const& collision,
                            Cands const& candsLcThisColl,
                            Pions const& tracksPionThisCollision,
                            float invMass2LcPiMin,
                            float invMass2LcPiMax)
  {
    auto primaryVertex = getPrimaryVertex(collision);
    auto covMatrixPV = primaryVertex.getCov();

    // Set the magnetic field from ccdb
    bz = collision.bz();
    df2.setBz(bz);

    for (const auto& candLc : candsLcThisColl) {
      auto trackParCovD = getTrackParCov(candLc);
      std::array<float, 3> pVecLc = candLc.pVector();

      for (const auto& trackPion : tracksPionThisCollision) {
        // this track is among daughters
        if (trackPion.trackId() == candLc.prong0Id() || trackPion.trackId() == candLc.prong1Id() || trackPion.trackId() == candLc.prong2Id()) {
          continue;
        }

        auto trackParCovPi = getTrackParCov(trackPion);
        std::array<float, 3> pVecPion = trackPion.pVector();

        // compute invariant mass square and apply selection
        auto invMass2LcPi = RecoDecay::m2(std::array{pVecLc, pVecPion}, std::array{MassLambdaCPlus, MassPiPlus});
        if ((invMass2LcPi < invMass2LcPiMin) || (invMass2LcPi > invMass2LcPiMax)) {
          continue;
        }
        // ---------------------------------
        // reconstruct the 2-prong Lb vertex
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

        // LcPi passed Lb reconstruction

        // calculate relevant properties
        const auto& secondaryVertexLb = df2.getPCACandidate();
        auto chi2PCA = df2.getChi2AtPCACandidate();
        auto covMatrixPCA = df2.calcPCACovMatrixFlat();
        registry.fill(HIST("hCovSVXX"), covMatrixPCA[0]);
        registry.fill(HIST("hCovPVXX"), covMatrixPV[0]);

        // propagate Lc and Pi to the Lb vertex
        df2.propagateTracksToVertex();
        // track.getPxPyPzGlo(pVec) modifies pVec of track
        df2.getTrack(0).getPxPyPzGlo(pVecLc);   // momentum of Lc at the Lb vertex
        df2.getTrack(1).getPxPyPzGlo(pVecPion); // momentum of Pi at the Lb vertex

        registry.fill(HIST("hMassLambdaB0ToLcPi"), std::sqrt(invMass2LcPi));

        // compute impact parameters of D and Pi
        o2::dataformats::DCA dcaLc;
        o2::dataformats::DCA dcaPion;
        trackParCovD.propagateToDCA(primaryVertex, bz, &dcaLc);
        trackParCovPi.propagateToDCA(primaryVertex, bz, &dcaPion);

        // get uncertainty of the decay length
        float phi{}, theta{};
        // getPointDirection modifies phi and theta
        getPointDirection(std::array{collision.posX(), collision.posY(), collision.posZ()}, secondaryVertexLb, phi, theta);
        auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
        auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));

        // fill the candidate table for the Lb here:
        rowCandidateBase(collision.globalIndex(),
                         collision.posX(), collision.posY(), collision.posZ(),
                         secondaryVertexLb[0], secondaryVertexLb[1], secondaryVertexLb[2],
                         errorDecayLength, errorDecayLengthXY,
                         chi2PCA,
                         pVecLc[0], pVecLc[1], pVecLc[2],
                         pVecPion[0], pVecPion[1], pVecPion[2],
                         dcaLc.getY(), dcaPion.getY(),
                         std::sqrt(dcaLc.getSigmaY2()), std::sqrt(dcaPion.getSigmaY2()));

        rowCandidateProngs(candLc.globalIndex(), trackPion.globalIndex());

        if constexpr (WithLcMl) {
          if (candLc.invMassHypo0() > 0) {
            rowCandidateLcMlScores(candLc.mlScoreBkgMassHypo0(), candLc.mlScorePromptMassHypo0(), candLc.mlScoreNonpromptMassHypo0());
          } else {
            rowCandidateLcMlScores(candLc.mlScoreBkgMassHypo1(), candLc.mlScorePromptMassHypo1(), candLc.mlScoreNonpromptMassHypo1());
          }
        } // pi loop
      } // Lc loop
    }
  }

  void processData(HfRedCollisionsWithExtras const& collisions,
                   soa::Join<aod::HfRed3Prongs, aod::HfRed3ProngsCov> const& candsLc,
                   soa::Join<aod::HfRedTrackBases, aod::HfRedTracksCov> const& tracksPion,
                   aod::HfOrigColCounts const& collisionsCounter,
                   aod::HfCandLbConfigs const& configs)
  {
    // LcPi invariant-mass window cut
    // invMassWindowLcPiTolerance is used to apply a slightly tighter cut than in LcPi pair preselection
    // to avoid accepting LcPi pairs that were not formed in LcPi pair creator
    auto [invMass2LcPiMin, invMass2LcPiMax] = computeInvMass2LcPiWindow(configs, invMassWindowLcPiTolerance);

    for (const auto& collisionCounter : collisionsCounter) {
      registry.fill(HIST("hEvents"), 1, collisionCounter.originalCollisionCount());
    }

    static int ncol = 0;
    static constexpr int PrintFrequency = 10000;
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsLcThisColl = candsLc.sliceBy(candsLcPerCollision, thisCollId);
      auto tracksPionThisCollision = tracksPion.sliceBy(tracksPionPerCollision, thisCollId);
      runCandidateCreation<false>(collision, candsLcThisColl, tracksPionThisCollision, invMass2LcPiMin, invMass2LcPiMax);
      if (ncol % PrintFrequency == 0) {
        LOGP(debug, "collisions parsed {}", ncol);
      }
      ncol++;
    }
  } // processData

  PROCESS_SWITCH(HfCandidateCreatorLbReduced, processData, "Process data without any ML score", true);

  void processDataWithLcMl(HfRedCollisionsWithExtras const& collisions,
                           soa::Join<aod::HfRed3Prongs, aod::HfRed3ProngsCov, aod::HfRed3ProngsMl> const& candsLc,
                           soa::Join<aod::HfRedTrackBases, aod::HfRedTracksCov> const& tracksPion,
                           aod::HfOrigColCounts const& collisionsCounter,
                           aod::HfCandLbConfigs const& configs)
  {
    // LcPi invariant-mass window cut
    // invMassWindowLcPiTolerance is used to apply a slightly tighter cut than in LcPi pair preselection
    // to avoid accepting LcPi pairs that were not formed in LcPi pair creator
    auto [invMass2LcPiMin, invMass2LcPiMax] = computeInvMass2LcPiWindow(configs, invMassWindowLcPiTolerance);

    for (const auto& collisionCounter : collisionsCounter) {
      registry.fill(HIST("hEvents"), 1, collisionCounter.originalCollisionCount());
    }

    static int ncol = 0;
    static constexpr int PrintFrequency = 10000;
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsLcThisColl = candsLc.sliceBy(candsLcPerCollision, thisCollId);
      auto tracksPionThisCollision = tracksPion.sliceBy(tracksPionPerCollision, thisCollId);
      runCandidateCreation<true>(collision, candsLcThisColl, tracksPionThisCollision, invMass2LcPiMin, invMass2LcPiMax);
      if (ncol % PrintFrequency == 0) {
        LOGP(debug, "collisions parsed {}", ncol);
      }
      ncol++;
    }
  } // processDataWithLcMl

  PROCESS_SWITCH(HfCandidateCreatorLbReduced, processDataWithLcMl, "Process data with ML scores of Lc", false);
}; // struct

/// Extends the table base with expression columns and performs MC matching.
struct HfCandidateCreatorLbReducedExpressions {
  Spawns<aod::HfCandLbExt> rowCandidateLb;
  Spawns<aod::HfRedTracksExt> rowTracksExt;
  Produces<aod::HfMcRecRedLbs> rowLbMcRec;
  Produces<aod::HfMcCheckLbs> rowLbMcCheck;

  /// Fill candidate information at MC reconstruction level
  /// \param checkDecayTypeMc
  /// \param rowsLcPiMcRec MC reco information on LcPi pairs
  /// \param candsLb prong global indices of Lb candidates
  template <bool CheckDecayTypeMc, typename McRec>
  void fillLbMcRec(McRec const& rowsLcPiMcRec, HfRedLbProngs const& candsLb)
  {
    for (const auto& candLb : candsLb) {
      bool filledMcInfo{false};
      for (const auto& rowLcPiMcRec : rowsLcPiMcRec) {
        if ((rowLcPiMcRec.prong0Id() != candLb.prong0Id()) || (rowLcPiMcRec.prong1Id() != candLb.prong1Id())) {
          continue;
        }
        rowLbMcRec(rowLcPiMcRec.flagMcMatchRec(), rowLcPiMcRec.flagWrongCollision(), rowLcPiMcRec.debugMcRec(), rowLcPiMcRec.ptMother());
        filledMcInfo = true;
        if constexpr (CheckDecayTypeMc) {
          rowLbMcCheck(rowLcPiMcRec.pdgCodeBeautyMother(),
                       rowLcPiMcRec.pdgCodeCharmMother(),
                       rowLcPiMcRec.pdgCodeProng0(),
                       rowLcPiMcRec.pdgCodeProng1(),
                       rowLcPiMcRec.pdgCodeProng2(),
                       rowLcPiMcRec.pdgCodeProng3());
        }
        break;
      }
      if (!filledMcInfo) { // protection to get same size tables in case something went wrong: we created a candidate that was not preselected in the LcPi creator
        rowLbMcRec(0, -1, -1, -1.f);
        if constexpr (CheckDecayTypeMc) {
          rowLbMcCheck(-1, -1, -1, -1, -1, -1);
        }
      }
    }
  }

  void processMc(HfMcRecRedLcPis const& rowsLcPiMcRec, HfRedLbProngs const& candsLb)
  {
    fillLbMcRec<false>(rowsLcPiMcRec, candsLb);
  }
  PROCESS_SWITCH(HfCandidateCreatorLbReducedExpressions, processMc, "Process MC", false);

  void processMcWithDecayTypeCheck(soa::Join<HfMcRecRedLcPis, HfMcCheckLcPis> const& rowsLcPiMcRec, HfRedLbProngs const& candsLb)
  {
    fillLbMcRec<true>(rowsLcPiMcRec, candsLb);
  }
  PROCESS_SWITCH(HfCandidateCreatorLbReducedExpressions, processMcWithDecayTypeCheck, "Process MC with decay type checks", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCandidateCreatorLbReduced>(cfgc),
                      adaptAnalysisTask<HfCandidateCreatorLbReducedExpressions>(cfgc)};
}
