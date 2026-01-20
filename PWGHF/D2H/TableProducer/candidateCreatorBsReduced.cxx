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

/// \file candidateCreatorBsReduced.cxx
/// \brief Reconstruction of Bs candidates
///
/// \author Fabio Catalano <fabio.catalano@cern.ch>, CERN

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

/// Reconstruction of Bs candidates
struct HfCandidateCreatorBsReduced {
  Produces<aod::HfCandBsBase> rowCandidateBase;         // table defined in CandidateReconstructionTables.h
  Produces<aod::HfRedBsProngs> rowCandidateProngs;      // table defined in ReducedDataModel.h
  Produces<aod::HfRedBsDsMls> rowCandidateDmesMlScores; // table defined in ReducedDataModel.h

  // vertexing
  Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCA{"useAbsDCA", true, "Minimise abs. distance rather than chi2"};
  Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  Configurable<float> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<float> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<float> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any B0 is smaller than this"};
  Configurable<float> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations is chi2/chi2old > this"};
  // selection
  Configurable<float> invMassWindowDPiTolerance{"invMassWindowDPiTolerance", 0.01, "invariant-mass window tolerance for DsPi pair preselections (GeV/c2)"};

  float myInvMassWindowDPi{1.}; // variable that will store the value of invMassWindowCharmHadPi (defined in dataCreatorCharmHadPiReduced.cxx)
  float massPi{o2::constants::physics::MassPiPlus};
  float massD{o2::constants::physics::MassDS};
  float massB{o2::constants::physics::MassBS};
  float bz{0.};

  o2::vertexing::DCAFitterN<2> df2; // fitter for B vertex (2-prong vertex fitter)

  using HfRedCollisionsWithExtras = soa::Join<aod::HfRedCollisions, aod::HfRedCollExtras>;

  Preslice<soa::Join<aod::HfRed3Prongs, aod::HfRed3ProngsCov>> candsDPerCollision = hf_track_index_reduced::hfRedCollisionId;
  Preslice<soa::Join<aod::HfRed3Prongs, aod::HfRed3ProngsCov, aod::HfRed3ProngsMl>> candsDWithMlPerCollision = hf_track_index_reduced::hfRedCollisionId;
  Preslice<soa::Join<aod::HfRedTrackBases, aod::HfRedTracksCov>> tracksPionPerCollision = hf_track_index_reduced::hfRedCollisionId;

  std::shared_ptr<TH1> hCandidates;
  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    std::array<bool, 2> doprocess{doprocessData, doprocessDataWithDmesMl};
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
    registry.add("hMassBsToDsPi", "2-prong candidates;inv. mass (B^{0}_{s} #rightarrow D_{s}^{#minus}#pi^{#plus} #rightarrow #K^{#minus}K^{#plus}#pi^{#minus}#pi^{#plus}) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 3., 8.}}});
    registry.add("hCovPVXX", "2-prong candidates;XX element of cov. matrix of prim. vtx. position (cm^{2});entries", {HistType::kTH1F, {{100, 0., 1.e-4}}});
    registry.add("hCovSVXX", "2-prong candidates;XX element of cov. matrix of sec. vtx. position (cm^{2});entries", {HistType::kTH1F, {{100, 0., 0.2}}});
    registry.add("hEvents", "Events;;entries", HistType::kTH1F, {{1, 0.5, 1.5}});

    /// candidate monitoring
    hCandidates = registry.add<TH1>("hCandidates", "candidates counter", {HistType::kTH1D, {axisCands}});
    setLabelHistoCands(hCandidates);
  }

  /// Main function to perform Bs candidate creation
  /// \param withDmesMl is the flag to use the table with ML scores for the Ds- daughter (only possible if present in the derived data)
  /// \param collision the collision
  /// \param candsDThisColl Bs candidates in this collision
  /// \param tracksPionThisCollision pion tracks in this collision
  /// \param invMass2DPiMin minimum Bs invariant-mass
  /// \param invMass2DPiMax maximum Bs invariant-mass
  template <bool WithDmesMl, typename Cands, typename Pions, typename Coll>
  void runCandidateCreation(Coll const& collision,
                            Cands const& candsDThisColl,
                            Pions const& tracksPionThisCollision,
                            const float& invMass2DPiMin,
                            const float& invMass2DPiMax)
  {
    auto primaryVertex = getPrimaryVertex(collision);
    auto covMatrixPV = primaryVertex.getCov();

    // Set the magnetic field from ccdb
    bz = collision.bz();
    df2.setBz(bz);

    for (const auto& candD : candsDThisColl) {
      auto trackParCovD = getTrackParCov(candD);
      std::array<float, 3> pVecD = candD.pVector();

      for (const auto& trackPion : tracksPionThisCollision) {
        // this track is among daughters
        if (trackPion.trackId() == candD.prong0Id() || trackPion.trackId() == candD.prong1Id() || trackPion.trackId() == candD.prong2Id()) {
          continue;
        }

        auto trackParCovPi = getTrackParCov(trackPion);
        std::array<float, 3> pVecPion = trackPion.pVector();

        // compute invariant mass square and apply selection
        auto invMass2DPi = RecoDecay::m2(std::array{pVecD, pVecPion}, std::array{massD, massPi});
        if ((invMass2DPi < invMass2DPiMin) || (invMass2DPi > invMass2DPiMax)) {
          continue;
        }
        // ---------------------------------
        // reconstruct the 2-prong Bs vertex
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
        hCandidates->Fill(SVFitting::FitOk); // DsPi passed Bs reconstruction

        // calculate relevant properties
        const auto& secondaryVertexB = df2.getPCACandidate();
        auto chi2PCA = df2.getChi2AtPCACandidate();
        auto covMatrixPCA = df2.calcPCACovMatrixFlat();
        registry.fill(HIST("hCovSVXX"), covMatrixPCA[0]);
        registry.fill(HIST("hCovPVXX"), covMatrixPV[0]);

        // get Ds and Pi tracks (propagated to the B0 vertex if propagateToPCA==true)
        // track.getPxPyPzGlo(pVec) modifies pVec of track
        df2.getTrack(0).getPxPyPzGlo(pVecD);    // momentum of Ds at the Bs vertex
        df2.getTrack(1).getPxPyPzGlo(pVecPion); // momentum of Pi at the Bs vertex

        registry.fill(HIST("hMassBsToDsPi"), std::sqrt(invMass2DPi));

        // compute impact parameters of Ds and Pi
        o2::dataformats::DCA dcaD;
        o2::dataformats::DCA dcaPion;
        trackParCovD.propagateToDCA(primaryVertex, bz, &dcaD);
        trackParCovPi.propagateToDCA(primaryVertex, bz, &dcaPion);

        // get uncertainty of the decay length
        float phi, theta;
        // getPointDirection modifies phi and theta
        getPointDirection(std::array{collision.posX(), collision.posY(), collision.posZ()}, secondaryVertexB, phi, theta);
        auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
        auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));

        // fill the candidate table for the Bs here:
        rowCandidateBase(collision.globalIndex(),
                         collision.posX(), collision.posY(), collision.posZ(),
                         secondaryVertexB[0], secondaryVertexB[1], secondaryVertexB[2],
                         errorDecayLength, errorDecayLengthXY,
                         chi2PCA,
                         pVecD[0], pVecD[1], pVecD[2],
                         pVecPion[0], pVecPion[1], pVecPion[2],
                         dcaD.getY(), dcaPion.getY(),
                         std::sqrt(dcaD.getSigmaY2()), std::sqrt(dcaPion.getSigmaY2()));

        rowCandidateProngs(candD.globalIndex(), trackPion.globalIndex());

        if constexpr (WithDmesMl) {
          if (candD.invMassHypo0() > 0) {
            rowCandidateDmesMlScores(candD.mlScoreBkgMassHypo0(), candD.mlScorePromptMassHypo0(), candD.mlScoreNonpromptMassHypo0());
          } else {
            rowCandidateDmesMlScores(candD.mlScoreBkgMassHypo1(), candD.mlScorePromptMassHypo1(), candD.mlScoreNonpromptMassHypo1());
          }
          // TODO: here we are assuming that only one of the two hypotheses is filled, to be checked
        }
      } // pi loop
    } // D loop
  }

  void processData(HfRedCollisionsWithExtras const& collisions,
                   soa::Join<aod::HfRed3Prongs, aod::HfRed3ProngsCov> const& candsD,
                   soa::Join<aod::HfRedTrackBases, aod::HfRedTracksCov> const& tracksPion,
                   aod::HfOrigColCounts const& collisionsCounter,
                   aod::HfCandBsConfigs const& configs)
  {
    // DsPi invariant-mass window cut
    for (const auto& config : configs) {
      myInvMassWindowDPi = config.myInvMassWindowDPi();
    }
    // invMassWindowDPiTolerance is used to apply a slightly tighter cut than in DsPi pair preselection
    // to avoid accepting DsPi pairs that were not formed in DsPi pair creator
    float const invMass2DPiMin = (massB - myInvMassWindowDPi + invMassWindowDPiTolerance) * (massB - myInvMassWindowDPi + invMassWindowDPiTolerance);
    float const invMass2DPiMax = (massB + myInvMassWindowDPi - invMassWindowDPiTolerance) * (massB + myInvMassWindowDPi - invMassWindowDPiTolerance);

    for (const auto& collisionCounter : collisionsCounter) {
      registry.fill(HIST("hEvents"), 1, collisionCounter.originalCollisionCount());
    }

    static int ncol = 0;
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsD.sliceBy(candsDPerCollision, thisCollId);
      auto tracksPionThisCollision = tracksPion.sliceBy(tracksPionPerCollision, thisCollId);
      runCandidateCreation<false>(collision, candsDThisColl, tracksPionThisCollision, invMass2DPiMin, invMass2DPiMax);
      if (ncol % 10000 == 0) {
        LOGP(debug, "collisions parsed {}", ncol);
      }
      ncol++;
    }
  } // processData

  PROCESS_SWITCH(HfCandidateCreatorBsReduced, processData, "Process data without any ML score", true);

  void processDataWithDmesMl(HfRedCollisionsWithExtras const& collisions,
                             soa::Join<aod::HfRed3Prongs, aod::HfRed3ProngsCov, aod::HfRed3ProngsMl> const& candsD,
                             soa::Join<aod::HfRedTrackBases, aod::HfRedTracksCov> const& tracksPion,
                             aod::HfOrigColCounts const& collisionsCounter,
                             aod::HfCandBsConfigs const& configs)
  {
    // DPi invariant-mass window cut
    for (const auto& config : configs) {
      myInvMassWindowDPi = config.myInvMassWindowDPi();
    }
    // invMassWindowDPiTolerance is used to apply a slightly tighter cut than in DsPi pair preselection
    // to avoid accepting DPi pairs that were not formed in DsPi pair creator
    float const invMass2DPiMin = (massB - myInvMassWindowDPi + invMassWindowDPiTolerance) * (massB - myInvMassWindowDPi + invMassWindowDPiTolerance);
    float const invMass2DPiMax = (massB + myInvMassWindowDPi - invMassWindowDPiTolerance) * (massB + myInvMassWindowDPi - invMassWindowDPiTolerance);

    for (const auto& collisionCounter : collisionsCounter) {
      registry.fill(HIST("hEvents"), 1, collisionCounter.originalCollisionCount());
    }

    static int ncol = 0;
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsD.sliceBy(candsDPerCollision, thisCollId);
      auto tracksPionThisCollision = tracksPion.sliceBy(tracksPionPerCollision, thisCollId);
      runCandidateCreation<true>(collision, candsDThisColl, tracksPionThisCollision, invMass2DPiMin, invMass2DPiMax);
      if (ncol % 10000 == 0) {
        LOGP(debug, "collisions parsed {}", ncol);
      }
      ncol++;
    }
  } // processDataWithDmesMl

  PROCESS_SWITCH(HfCandidateCreatorBsReduced, processDataWithDmesMl, "Process data with ML scores of D mesons", false);
}; // struct

/// Extends the table base with expression columns and performs MC matching.
struct HfCandidateCreatorBsReducedExpressions {
  Spawns<aod::HfCandBsExt> rowCandidateBs;
  Spawns<aod::HfRedTracksExt> rowTracksExt;
  Produces<aod::HfMcRecRedBss> rowBsMcRec;
  Produces<aod::HfMcCheckBss> rowBsMcCheck;

  /// Fill candidate information at MC reconstruction level
  /// \param checkDecayTypeMc
  /// \param rowsDPiMcRec MC reco information on DsPi pairs
  /// \param candsB prong global indices of Bs candidates
  template <bool CheckDecayTypeMc, typename McRec>
  void fillBsMcRec(McRec const& rowsDPiMcRec, HfRedBsProngs const& candsB)
  {
    for (const auto& candB : candsB) {
      bool filledMcInfo{false};
      for (const auto& rowDPiMcRec : rowsDPiMcRec) {
        if ((rowDPiMcRec.prong0Id() != candB.prong0Id()) || (rowDPiMcRec.prong1Id() != candB.prong1Id())) {
          continue;
        }
        rowBsMcRec(rowDPiMcRec.flagMcMatchRec(), -1 /*channel*/, rowDPiMcRec.flagWrongCollision(), rowDPiMcRec.debugMcRec(), rowDPiMcRec.ptMother());
        filledMcInfo = true;
        if constexpr (CheckDecayTypeMc) {
          rowBsMcCheck(rowDPiMcRec.pdgCodeBeautyMother(),
                       rowDPiMcRec.pdgCodeCharmMother(),
                       rowDPiMcRec.pdgCodeProng0(),
                       rowDPiMcRec.pdgCodeProng1(),
                       rowDPiMcRec.pdgCodeProng2(),
                       rowDPiMcRec.pdgCodeProng3());
        }
        break;
      }
      if (!filledMcInfo) { // protection to get same size tables in case something went wrong: we created a candidate that was not preselected in the DsPi creator
        rowBsMcRec(0, -1, -1, -1, -1.f);
        if constexpr (CheckDecayTypeMc) {
          rowBsMcCheck(-1, -1, -1, -1, -1, -1);
        }
      }
    }
  }

  void processMc(HfMcRecRedDsPis const& rowsDPiMcRec, HfRedBsProngs const& candsB)
  {
    fillBsMcRec<false>(rowsDPiMcRec, candsB);
  }
  PROCESS_SWITCH(HfCandidateCreatorBsReducedExpressions, processMc, "Process MC", false);

  void processMcWithDecayTypeCheck(soa::Join<HfMcRecRedDsPis, HfMcCheckDsPis> const& rowsDPiMcRec, HfRedBsProngs const& candsB)
  {
    fillBsMcRec<true>(rowsDPiMcRec, candsB);
  }
  PROCESS_SWITCH(HfCandidateCreatorBsReducedExpressions, processMcWithDecayTypeCheck, "Process MC with decay type checks", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCandidateCreatorBsReduced>(cfgc),
                      adaptAnalysisTask<HfCandidateCreatorBsReducedExpressions>(cfgc)};
}
