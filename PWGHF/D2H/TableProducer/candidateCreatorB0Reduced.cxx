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

/// \file candidateCreatorB0Reduced.cxx
/// \brief Reconstruction of B0 candidates
///
/// \author Alexandre Bigot <alexandre.bigot@cern.ch>, IPHC Strasbourg
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

#include <memory>

#include "CommonConstants/PhysicsConstants.h"
#include "DCAFitter/DCAFitterN.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"

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

/// Reconstruction of B0 candidates
struct HfCandidateCreatorB0Reduced {
  Produces<aod::HfCandB0Base> rowCandidateBase;         // table defined in CandidateReconstructionTables.h
  Produces<aod::HfRedB0Prongs> rowCandidateProngs;      // table defined in ReducedDataModel.h
  Produces<aod::HfRedB0DpMls> rowCandidateDmesMlScores; // table defined in ReducedDataModel.h

  // vertexing
  Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCA{"useAbsDCA", true, "Minimise abs. distance rather than chi2"};
  Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  Configurable<float> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<float> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<float> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any B0 is smaller than this"};
  Configurable<float> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations is chi2/chi2old > this"};
  // selection
  Configurable<float> invMassWindowDPiTolerance{"invMassWindowDPiTolerance", 0.01, "invariant-mass window tolerance for DPi pair preselections (GeV/c2)"};

  float myInvMassWindowDPi{1.}; // variable that will store the value of invMassWindowDPi (defined in dataCreatorDplusPiReduced.cxx)
  float massPi{0.};
  float massD{0.};
  float massB0{0.};
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

    // invariant-mass window cut
    massPi = o2::constants::physics::MassPiPlus;
    massD = o2::constants::physics::MassDMinus;
    massB0 = o2::constants::physics::MassB0;

    // Initialize fitter
    df2.setPropagateToPCA(propagateToPCA);
    df2.setMaxR(maxR);
    df2.setMaxDZIni(maxDZIni);
    df2.setMinParamChange(minParamChange);
    df2.setMinRelChi2Change(minRelChi2Change);
    df2.setUseAbsDCA(useAbsDCA);
    df2.setWeightedFinalPCA(useWeightedFinalPCA);

    // histograms
    registry.add("hMassB0ToDPi", "2-prong candidates;inv. mass (B^{0} #rightarrow D^{#minus}#pi^{#plus} #rightarrow #pi^{#minus}K^{#plus}#pi^{#minus}#pi^{#plus}) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 3., 8.}}});
    registry.add("hCovPVXX", "2-prong candidates;XX element of cov. matrix of prim. vtx. position (cm^{2});entries", {HistType::kTH1F, {{100, 0., 1.e-4}}});
    registry.add("hCovSVXX", "2-prong candidates;XX element of cov. matrix of sec. vtx. position (cm^{2});entries", {HistType::kTH1F, {{100, 0., 0.2}}});
    registry.add("hEvents", "Events;;entries", HistType::kTH1F, {{1, 0.5, 1.5}});

    /// candidate monitoring
    hCandidates = registry.add<TH1>("hCandidates", "candidates counter", {HistType::kTH1D, {axisCands}});
    setLabelHistoCands(hCandidates);
  }

  /// Main function to perform B0 candidate creation
  /// \param withDmesMl is the flag to use the table with ML scores for the D- daughter (only possible if present in the derived data)
  /// \param collision the collision
  /// \param candsDThisColl B0 candidates in this collision
  /// \param tracksPionThisCollision pion tracks in this collision
  /// \param invMass2DPiMin minimum B0 invariant-mass
  /// \param invMass2DPiMax maximum B0 invariant-mass
  template <bool withDmesMl, typename Cands, typename Pions, typename Coll>
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
        // reconstruct the 2-prong B0 vertex
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

        // DPi passed B0 reconstruction

        // calculate relevant properties
        const auto& secondaryVertexB0 = df2.getPCACandidate();
        auto chi2PCA = df2.getChi2AtPCACandidate();
        auto covMatrixPCA = df2.calcPCACovMatrixFlat();
        registry.fill(HIST("hCovSVXX"), covMatrixPCA[0]);
        registry.fill(HIST("hCovPVXX"), covMatrixPV[0]);

        // propagate D and Pi to the B0 vertex
        df2.propagateTracksToVertex();
        // track.getPxPyPzGlo(pVec) modifies pVec of track
        df2.getTrack(0).getPxPyPzGlo(pVecD);    // momentum of D at the B0 vertex
        df2.getTrack(1).getPxPyPzGlo(pVecPion); // momentum of Pi at the B0 vertex

        registry.fill(HIST("hMassB0ToDPi"), std::sqrt(invMass2DPi));

        // compute impact parameters of D and Pi
        o2::dataformats::DCA dcaD;
        o2::dataformats::DCA dcaPion;
        trackParCovD.propagateToDCA(primaryVertex, bz, &dcaD);
        trackParCovPi.propagateToDCA(primaryVertex, bz, &dcaPion);

        // get uncertainty of the decay length
        float phi, theta;
        // getPointDirection modifies phi and theta
        getPointDirection(std::array{collision.posX(), collision.posY(), collision.posZ()}, secondaryVertexB0, phi, theta);
        auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
        auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));

        // fill the candidate table for the B0 here:
        rowCandidateBase(collision.globalIndex(),
                         collision.posX(), collision.posY(), collision.posZ(),
                         secondaryVertexB0[0], secondaryVertexB0[1], secondaryVertexB0[2],
                         errorDecayLength, errorDecayLengthXY,
                         chi2PCA,
                         pVecD[0], pVecD[1], pVecD[2],
                         pVecPion[0], pVecPion[1], pVecPion[2],
                         dcaD.getY(), dcaPion.getY(),
                         std::sqrt(dcaD.getSigmaY2()), std::sqrt(dcaPion.getSigmaY2()));

        rowCandidateProngs(candD.globalIndex(), trackPion.globalIndex());

        if constexpr (withDmesMl) {
          rowCandidateDmesMlScores(candD.mlScoreBkgMassHypo0(), candD.mlScorePromptMassHypo0(), candD.mlScoreNonpromptMassHypo0());
        }
      } // pi loop
    }   // D loop
  }

  void processData(HfRedCollisionsWithExtras const& collisions,
                   soa::Join<aod::HfRed3Prongs, aod::HfRed3ProngsCov> const& candsD,
                   soa::Join<aod::HfRedTrackBases, aod::HfRedTracksCov> const& tracksPion,
                   aod::HfOrigColCounts const& collisionsCounter,
                   aod::HfCandB0Configs const& configs)
  {
    // DPi invariant-mass window cut
    for (const auto& config : configs) {
      myInvMassWindowDPi = config.myInvMassWindowDPi();
    }
    // invMassWindowDPiTolerance is used to apply a slightly tighter cut than in DPi pair preselection
    // to avoid accepting DPi pairs that were not formed in DPi pair creator
    float invMass2DPiMin = (massB0 - myInvMassWindowDPi + invMassWindowDPiTolerance) * (massB0 - myInvMassWindowDPi + invMassWindowDPiTolerance);
    float invMass2DPiMax = (massB0 + myInvMassWindowDPi - invMassWindowDPiTolerance) * (massB0 + myInvMassWindowDPi - invMassWindowDPiTolerance);

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

  PROCESS_SWITCH(HfCandidateCreatorB0Reduced, processData, "Process data without any ML score", true);

  void processDataWithDmesMl(HfRedCollisionsWithExtras const& collisions,
                             soa::Join<aod::HfRed3Prongs, aod::HfRed3ProngsCov, aod::HfRed3ProngsMl> const& candsD,
                             soa::Join<aod::HfRedTrackBases, aod::HfRedTracksCov> const& tracksPion,
                             aod::HfOrigColCounts const& collisionsCounter,
                             aod::HfCandB0Configs const& configs)
  {
    // DPi invariant-mass window cut
    for (const auto& config : configs) {
      myInvMassWindowDPi = config.myInvMassWindowDPi();
    }
    // invMassWindowDPiTolerance is used to apply a slightly tighter cut than in DPi pair preselection
    // to avoid accepting DPi pairs that were not formed in DPi pair creator
    float invMass2DPiMin = (massB0 - myInvMassWindowDPi + invMassWindowDPiTolerance) * (massB0 - myInvMassWindowDPi + invMassWindowDPiTolerance);
    float invMass2DPiMax = (massB0 + myInvMassWindowDPi - invMassWindowDPiTolerance) * (massB0 + myInvMassWindowDPi - invMassWindowDPiTolerance);

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

  PROCESS_SWITCH(HfCandidateCreatorB0Reduced, processDataWithDmesMl, "Process data with ML scores of D mesons", false);
}; // struct

/// Extends the table base with expression columns and performs MC matching.
struct HfCandidateCreatorB0ReducedExpressions {
  Spawns<aod::HfCandB0Ext> rowCandidateB0;
  Spawns<aod::HfRedTracksExt> rowTracksExt;
  Produces<aod::HfMcRecRedB0s> rowB0McRec;
  Produces<aod::HfMcCheckB0s> rowB0McCheck;

  /// Fill candidate information at MC reconstruction level
  /// \param checkDecayTypeMc
  /// \param rowsDPiMcRec MC reco information on DPi pairs
  /// \param candsB0 prong global indices of B0 candidates
  template <bool checkDecayTypeMc, typename McRec>
  void fillB0McRec(McRec const& rowsDPiMcRec, HfRedB0Prongs const& candsB0)
  {
    for (const auto& candB0 : candsB0) {
      bool filledMcInfo{false};
      for (const auto& rowDPiMcRec : rowsDPiMcRec) {
        if ((rowDPiMcRec.prong0Id() != candB0.prong0Id()) || (rowDPiMcRec.prong1Id() != candB0.prong1Id())) {
          continue;
        }
        rowB0McRec(rowDPiMcRec.flagMcMatchRec(), rowDPiMcRec.flagWrongCollision(), rowDPiMcRec.debugMcRec(), rowDPiMcRec.ptMother());
        filledMcInfo = true;
        if constexpr (checkDecayTypeMc) {
          rowB0McCheck(rowDPiMcRec.pdgCodeBeautyMother(),
                       rowDPiMcRec.pdgCodeCharmMother(),
                       rowDPiMcRec.pdgCodeProng0(),
                       rowDPiMcRec.pdgCodeProng1(),
                       rowDPiMcRec.pdgCodeProng2(),
                       rowDPiMcRec.pdgCodeProng3());
        }
        break;
      }
      if (!filledMcInfo) { // protection to get same size tables in case something went wrong: we created a candidate that was not preselected in the D-Pi creator
        rowB0McRec(0, -1, -1, -1.f);
        if constexpr (checkDecayTypeMc) {
          rowB0McCheck(-1, -1, -1, -1, -1, -1);
        }
      }
    }
  }

  void processMc(HfMcRecRedDpPis const& rowsDPiMcRec, HfRedB0Prongs const& candsB0)
  {
    fillB0McRec<false>(rowsDPiMcRec, candsB0);
  }
  PROCESS_SWITCH(HfCandidateCreatorB0ReducedExpressions, processMc, "Process MC", false);

  void processMcWithDecayTypeCheck(soa::Join<HfMcRecRedDpPis, HfMcCheckDpPis> const& rowsDPiMcRec, HfRedB0Prongs const& candsB0)
  {
    fillB0McRec<true>(rowsDPiMcRec, candsB0);
  }
  PROCESS_SWITCH(HfCandidateCreatorB0ReducedExpressions, processMcWithDecayTypeCheck, "Process MC with decay type checks", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCandidateCreatorB0Reduced>(cfgc),
                      adaptAnalysisTask<HfCandidateCreatorB0ReducedExpressions>(cfgc)};
}
