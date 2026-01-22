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

/// Reconstruction of B0 candidates
struct HfCandidateCreatorB0Reduced {
  Produces<aod::HfCandB0Base> rowCandidateBase;              // table defined in CandidateReconstructionTables.h
  Produces<aod::HfCandB0DStar> rowCandidateB0DStar;          // table defined in CandidateReconstructionTables.h
  Produces<aod::HfRedB0Prongs> rowCandidateProngs;           // table defined in ReducedDataModel.h
  Produces<aod::HfRedB0ProngDStars> rowCandidateProngsDStar; // table defined in ReducedDataModel.h
  Produces<aod::HfRedB0DpMls> rowCandidateDmesMlScores;      // table defined in ReducedDataModel.h

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
  float bz{0.};

  o2::vertexing::DCAFitterN<2> df2; // fitter for B vertex (2-prong vertex fitter for DpPi candidates)
  o2::vertexing::DCAFitterN<3> df3; // fitter for B vertex (3-prong vertex fitter for DstarPi candidates)

  using HfRedCollisionsWithExtras = soa::Join<aod::HfRedCollisions, aod::HfRedCollExtras>;
  using HfSoftPiWCovAndPid = soa::Join<aod::HfRedSoftPiBases, aod::HfRedSoftPiCov, aod::HfRedSoftPiPid>;

  Preslice<soa::Join<aod::HfRed3Prongs, aod::HfRed3ProngsCov>> candsDplusPerCollision = hf_track_index_reduced::hfRedCollisionId;
  Preslice<soa::Join<aod::HfRed3Prongs, aod::HfRed3ProngsCov, aod::HfRed3ProngsMl>> candsDplusWithMlPerCollision = hf_track_index_reduced::hfRedCollisionId;
  Preslice<soa::Join<aod::HfRed2Prongs, aod::HfRed2ProngsCov>> candsDstarPerCollision = hf_track_index_reduced::hfRedCollisionId;
  Preslice<soa::Join<aod::HfRed2Prongs, aod::HfRed2ProngsCov, aod::HfRed3ProngsMl>> candsDstarWithMlPerCollision = hf_track_index_reduced::hfRedCollisionId;
  Preslice<soa::Join<aod::HfRedTrackBases, aod::HfRedTracksCov>> tracksPionPerCollision = hf_track_index_reduced::hfRedCollisionId;

  std::shared_ptr<TH1> hCandidates;
  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    std::array<bool, 4> doprocess{doprocessDataDplusPi, doprocessDataDplusPiWithDmesMl, doprocessDataDstarPi, doprocessDataDstarPiWithDmesMl};
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

    df3.setPropagateToPCA(propagateToPCA);
    df3.setMaxR(maxR);
    df3.setMaxDZIni(maxDZIni);
    df3.setMinParamChange(minParamChange);
    df3.setMinRelChi2Change(minRelChi2Change);
    df3.setUseAbsDCA(useAbsDCA);
    df3.setWeightedFinalPCA(useWeightedFinalPCA);

    // histograms
    if (doprocessDataDplusPi || doprocessDataDplusPiWithDmesMl) {
      registry.add("hMassB0ToDPi", "2-prong candidates;inv. mass (B^{0} #rightarrow D^{#minus}#pi^{#plus} #rightarrow #pi^{#minus}K^{#plus}#pi^{#minus}#pi^{#plus}) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{500, 3., 8.}}});
    } else if (doprocessDataDstarPi || doprocessDataDstarPiWithDmesMl) {
      registry.add("hMassB0ToDPi", "3-prong candidates;inv. mass (B^{0} #rightarrow D^{0}#pi^{#plus}#pi^{#plus} #rightarrow #pi^{#minus}K^{#plus}#pi^{#minus}#pi^{#plus}) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{500, 3., 8.}}});
    }
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
  template <bool WithDmesMl, typename Cands, typename Pions, typename SoftPions, typename Coll>
  void runCandidateCreationDStar(Coll const& collision,
                                 Cands const& candsDThisColl,
                                 SoftPions const& softPions,
                                 Pions const& tracksPionThisCollision,
                                 const float invMass2DPiMin,
                                 const float invMass2DPiMax)
  {
    auto primaryVertex = getPrimaryVertex(collision);
    auto covMatrixPV = primaryVertex.getCov();

    // Set the magnetic field from ccdb
    bz = collision.bz();
    df3.setBz(bz);

    for (const auto& candD : candsDThisColl) {
      auto trackParCovD = getTrackParCov(candD);
      std::array<float, 3> pVecD0 = candD.pVector();
      auto trackSoftPion = softPions.rawIteratorAt(candD.globalIndex());
      std::array<float, 3> pVecSoftPion = trackSoftPion.pVector();
      auto trackParCovSoftPi = getTrackParCov(trackSoftPion);

      for (const auto& trackBachPion : tracksPionThisCollision) {
        // this track is among daughters
        if (trackBachPion.trackId() == candD.prong0Id() || trackBachPion.trackId() == candD.prong1Id() || trackBachPion.trackId() == trackSoftPion.trackId()) {
          continue;
        }

        auto trackParCovPi = getTrackParCov(trackBachPion);
        std::array<float, 3> pVecBachPion = trackBachPion.pVector();

        // compute invariant mass square and apply selection
        auto invMass2DPi = RecoDecay::m2(std::array{pVecD0, pVecSoftPion, pVecBachPion}, std::array{o2::constants::physics::MassD0, o2::constants::physics::MassPiPlus, o2::constants::physics::MassPiPlus});
        if ((invMass2DPi < invMass2DPiMin) || (invMass2DPi > invMass2DPiMax)) {
          continue;
        }

        // ---------------------------------
        // reconstruct the 3-prong B0 vertex
        hCandidates->Fill(SVFitting::BeforeFit);
        try {
          if (df3.process(trackParCovD, trackParCovSoftPi, trackParCovPi) == 0) {
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
        const auto& secondaryVertexB0 = df3.getPCACandidate();
        auto chi2PCA = df3.getChi2AtPCACandidate();
        auto covMatrixPCA = df3.calcPCACovMatrixFlat();
        registry.fill(HIST("hCovSVXX"), covMatrixPCA[0]);
        registry.fill(HIST("hCovPVXX"), covMatrixPV[0]);

        // get D and Pi tracks (propagated to the B0 vertex if propagateToPCA==true)
        // track.getPxPyPzGlo(pVec) modifies pVec of track
        df3.getTrack(0).getPxPyPzGlo(pVecD0);       // momentum of D at the B0 vertex
        df3.getTrack(1).getPxPyPzGlo(pVecSoftPion); // momentum of SoftPi at the B0 vertex
        df3.getTrack(2).getPxPyPzGlo(pVecBachPion); // momentum of Pi at the B0 vertex

        registry.fill(HIST("hMassB0ToDPi"), std::sqrt(invMass2DPi));

        // compute impact parameters of D and Pi
        o2::dataformats::DCA dcaD;
        o2::dataformats::DCA dcaSoftPion;
        o2::dataformats::DCA dcaBachPion;
        trackParCovD.propagateToDCA(primaryVertex, bz, &dcaD);
        trackParCovSoftPi.propagateToDCA(primaryVertex, bz, &dcaSoftPion);
        trackParCovPi.propagateToDCA(primaryVertex, bz, &dcaBachPion);

        // get uncertainty of the decay length
        float phi, theta;
        // getPointDirection modifies phi and theta
        getPointDirection(std::array{collision.posX(), collision.posY(), collision.posZ()}, secondaryVertexB0, phi, theta);
        auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
        auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));

        // fill the candidate table for the B0 here:
        rowCandidateB0DStar(collision.globalIndex(),
                            collision.posX(), collision.posY(), collision.posZ(),
                            secondaryVertexB0[0], secondaryVertexB0[1], secondaryVertexB0[2],
                            errorDecayLength, errorDecayLengthXY,
                            chi2PCA,
                            pVecD0[0], pVecD0[1], pVecD0[2],
                            pVecSoftPion[0], pVecSoftPion[1], pVecSoftPion[2],
                            pVecBachPion[0], pVecBachPion[1], pVecBachPion[2],
                            dcaD.getY(), dcaSoftPion.getY(), dcaBachPion.getY(),
                            std::sqrt(dcaD.getSigmaY2()), std::sqrt(dcaSoftPion.getSigmaY2()), std::sqrt(dcaBachPion.getSigmaY2()));

        rowCandidateProngsDStar(candD.globalIndex(), trackBachPion.globalIndex(), trackSoftPion.globalIndex());

        if constexpr (WithDmesMl) {
          rowCandidateDmesMlScores(candD.mlScoreBkgMassHypo0(), candD.mlScorePromptMassHypo0(), candD.mlScoreNonpromptMassHypo0());
        }
      } // pi loop
    } // D loop
  }

  /// Main function to perform B0 candidate creation
  /// \param withDmesMl is the flag to use the table with ML scores for the D- daughter (only possible if present in the derived data)
  /// \param collision the collision
  /// \param candsDThisColl B0 candidates in this collision
  /// \param tracksPionThisCollision pion tracks in this collision
  /// \param invMass2DPiMin minimum B0 invariant-mass
  /// \param invMass2DPiMax maximum B0 invariant-mass
  template <bool WithDmesMl, typename Cands, typename Pions, typename Coll>
  void runCandidateCreation(Coll const& collision,
                            Cands const& candsDThisColl,
                            Pions const& tracksPionThisCollision,
                            const float invMass2DPiMin,
                            const float invMass2DPiMax)
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
        auto invMass2DPi = RecoDecay::m2(std::array{pVecD, pVecPion}, std::array{o2::constants::physics::MassDMinus, o2::constants::physics::MassPiPlus});
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

        // get D and Pi tracks (propagated to the B0 vertex if propagateToPCA==true)
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

        if constexpr (WithDmesMl) {
          rowCandidateDmesMlScores(candD.mlScoreBkgMassHypo0(), candD.mlScorePromptMassHypo0(), candD.mlScoreNonpromptMassHypo0());
        }
      } // pi loop
    } // D loop
  }

  void processDataDplusPi(HfRedCollisionsWithExtras const& collisions,
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
    float const invMass2DPiMin = (o2::constants::physics::MassB0 - myInvMassWindowDPi + invMassWindowDPiTolerance) * (o2::constants::physics::MassB0 - myInvMassWindowDPi + invMassWindowDPiTolerance);
    float const invMass2DPiMax = (o2::constants::physics::MassB0 + myInvMassWindowDPi - invMassWindowDPiTolerance) * (o2::constants::physics::MassB0 + myInvMassWindowDPi - invMassWindowDPiTolerance);

    for (const auto& collisionCounter : collisionsCounter) {
      registry.fill(HIST("hEvents"), 1, collisionCounter.originalCollisionCount());
    }

    static int ncol = 0;
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsD.sliceBy(candsDplusPerCollision, thisCollId);
      auto tracksPionThisCollision = tracksPion.sliceBy(tracksPionPerCollision, thisCollId);
      runCandidateCreation<false>(collision, candsDThisColl, tracksPionThisCollision, invMass2DPiMin, invMass2DPiMax);
      if (ncol % 10000 == 0) {
        LOGP(debug, "collisions parsed {}", ncol);
      }
      ncol++;
    }
  } // processDataDplusPi

  PROCESS_SWITCH(HfCandidateCreatorB0Reduced, processDataDplusPi, "Process data D-pi without any ML score", true);

  void processDataDplusPiWithDmesMl(HfRedCollisionsWithExtras const& collisions,
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
    float const invMass2DPiMin = (o2::constants::physics::MassB0 - myInvMassWindowDPi + invMassWindowDPiTolerance) * (o2::constants::physics::MassB0 - myInvMassWindowDPi + invMassWindowDPiTolerance);
    float const invMass2DPiMax = (o2::constants::physics::MassB0 + myInvMassWindowDPi - invMassWindowDPiTolerance) * (o2::constants::physics::MassB0 + myInvMassWindowDPi - invMassWindowDPiTolerance);

    for (const auto& collisionCounter : collisionsCounter) {
      registry.fill(HIST("hEvents"), 1, collisionCounter.originalCollisionCount());
    }

    static int ncol = 0;
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsD.sliceBy(candsDplusPerCollision, thisCollId);
      auto tracksPionThisCollision = tracksPion.sliceBy(tracksPionPerCollision, thisCollId);
      runCandidateCreation<true>(collision, candsDThisColl, tracksPionThisCollision, invMass2DPiMin, invMass2DPiMax);
      if (ncol % 10000 == 0) {
        LOGP(debug, "collisions parsed {}", ncol);
      }
      ncol++;
    }
  } // processDataDplusPiWithDmesMl

  PROCESS_SWITCH(HfCandidateCreatorB0Reduced, processDataDplusPiWithDmesMl, "Process data D-pi with ML scores of D mesons", false);

  void processDataDstarPi(HfRedCollisionsWithExtras const& collisions,
                          soa::Join<aod::HfRed2Prongs, aod::HfRed2ProngsCov> const& candsD,
                          HfSoftPiWCovAndPid const& softPions,
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
    float const invMass2DPiMin = (o2::constants::physics::MassB0 - myInvMassWindowDPi + invMassWindowDPiTolerance) * (o2::constants::physics::MassB0 - myInvMassWindowDPi + invMassWindowDPiTolerance);
    float const invMass2DPiMax = (o2::constants::physics::MassB0 + myInvMassWindowDPi - invMassWindowDPiTolerance) * (o2::constants::physics::MassB0 + myInvMassWindowDPi - invMassWindowDPiTolerance);

    for (const auto& collisionCounter : collisionsCounter) {
      registry.fill(HIST("hEvents"), 1, collisionCounter.originalCollisionCount());
    }

    static int ncol = 0;
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsD.sliceBy(candsDstarPerCollision, thisCollId);
      auto tracksPionThisCollision = tracksPion.sliceBy(tracksPionPerCollision, thisCollId);
      runCandidateCreationDStar<false>(collision, candsDThisColl, softPions, tracksPionThisCollision, invMass2DPiMin, invMass2DPiMax);
      if (ncol % 10000 == 0) {
        LOGP(debug, "collisions parsed {}", ncol);
      }
      ncol++;
    }
  } // processDataDstarPi

  PROCESS_SWITCH(HfCandidateCreatorB0Reduced, processDataDstarPi, "Process data D*pi without any ML score", false);

  void processDataDstarPiWithDmesMl(HfRedCollisionsWithExtras const& collisions,
                                    soa::Join<aod::HfRed2Prongs, aod::HfRed2ProngsCov, aod::HfRed3ProngsMl> const& candsD,
                                    HfSoftPiWCovAndPid const& softPions,
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
    float const invMass2DPiMin = (o2::constants::physics::MassB0 - myInvMassWindowDPi + invMassWindowDPiTolerance) * (o2::constants::physics::MassB0 - myInvMassWindowDPi + invMassWindowDPiTolerance);
    float const invMass2DPiMax = (o2::constants::physics::MassB0 + myInvMassWindowDPi - invMassWindowDPiTolerance) * (o2::constants::physics::MassB0 + myInvMassWindowDPi - invMassWindowDPiTolerance);

    for (const auto& collisionCounter : collisionsCounter) {
      registry.fill(HIST("hEvents"), 1, collisionCounter.originalCollisionCount());
    }

    static int ncol = 0;
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsD.sliceBy(candsDstarWithMlPerCollision, thisCollId);
      auto tracksPionThisCollision = tracksPion.sliceBy(tracksPionPerCollision, thisCollId);
      runCandidateCreationDStar<true>(collision, candsDThisColl, softPions, tracksPionThisCollision, invMass2DPiMin, invMass2DPiMax);
      if (ncol % 10000 == 0) {
        LOGP(debug, "collisions parsed {}", ncol);
      }
      ncol++;
    }
  } // processDataDstarPiWithDmesMl

  PROCESS_SWITCH(HfCandidateCreatorB0Reduced, processDataDstarPiWithDmesMl, "Process data D*pi with ML scores of D mesons", false);

}; // struct

/// Extends the table base with expression columns and performs MC matching.
struct HfCandidateCreatorB0ReducedExpressions {
  Spawns<aod::HfCandB0Ext> rowCandidateB0;
  Spawns<aod::HfCandB0DStExt> rowCandidateB0DSt;
  Spawns<aod::HfRedTracksExt> rowTracksExt;
  Produces<aod::HfMcRecRedB0s> rowB0McRec;
  Produces<aod::HfMcCheckB0s> rowB0McCheck;

  /// Fill candidate information at MC reconstruction level
  /// \param checkDecayTypeMc
  /// \param rowsDPiMcRec MC reco information on DPi pairs
  /// \param candsB0 prong global indices of B0 candidates
  template <bool CheckDecayTypeMc, typename McRec, typename B0Prongs>
  void fillB0McRec(McRec const& rowsDPiMcRec, B0Prongs const& candsB0)
  {
    for (const auto& candB0 : candsB0) {
      bool filledMcInfo{false};
      for (const auto& rowDPiMcRec : rowsDPiMcRec) {
        if constexpr (std::is_same_v<B0Prongs, aod::HfRedB0Prongs>) {
          if ((rowDPiMcRec.prong0Id() != candB0.prong0Id()) || (rowDPiMcRec.prong1Id() != candB0.prong1Id())) {
            continue;
          }
        } else if constexpr (std::is_same_v<B0Prongs, aod::HfRedB0ProngDStars>) { // No need to check ID of soft pion, it is the same as D0
          if ((rowDPiMcRec.prongD0Id() != candB0.prongD0Id()) || (rowDPiMcRec.prongBachPiId() != candB0.prongBachPiId())) {
            continue;
          }
        }
        rowB0McRec(rowDPiMcRec.flagMcMatchRec(), -1 /*channel*/, rowDPiMcRec.flagWrongCollision(), rowDPiMcRec.debugMcRec(), rowDPiMcRec.ptMother());
        filledMcInfo = true;
        if constexpr (CheckDecayTypeMc) {
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
        rowB0McRec(0, -1, -1, -1, -1.f);
        if constexpr (CheckDecayTypeMc) {
          rowB0McCheck(-1, -1, -1, -1, -1, -1);
        }
      }
    }
  }

  void processMcDplusPi(HfMcRecRedDpPis const& rowsDPiMcRec, HfRedB0Prongs const& candsB0)
  {
    fillB0McRec<false>(rowsDPiMcRec, candsB0);
  }
  PROCESS_SWITCH(HfCandidateCreatorB0ReducedExpressions, processMcDplusPi, "Process MC for DplusPi", false);

  void processMcDplusPiWithDecayTypeCheck(soa::Join<HfMcRecRedDpPis, HfMcCheckDpPis> const& rowsDPiMcRec, HfRedB0Prongs const& candsB0)
  {
    fillB0McRec<true>(rowsDPiMcRec, candsB0);
  }
  PROCESS_SWITCH(HfCandidateCreatorB0ReducedExpressions, processMcDplusPiWithDecayTypeCheck, "Process MC with decay type checks for DplusPi", false);

  void processMcDstarPi(HfMcRecRedDStarPis const& rowsDPiMcRec, HfRedB0ProngDStars const& candsB0)
  {
    fillB0McRec<false>(rowsDPiMcRec, candsB0);
  }
  PROCESS_SWITCH(HfCandidateCreatorB0ReducedExpressions, processMcDstarPi, "Process MC for DstarPi", false);

  void processMcDstarPiWithDecayTypeCheck(soa::Join<HfMcRecRedDStarPis, HfMcCheckDpPis> const& rowsDPiMcRec, HfRedB0ProngDStars const& candsB0)
  {
    fillB0McRec<true>(rowsDPiMcRec, candsB0);
  }
  PROCESS_SWITCH(HfCandidateCreatorB0ReducedExpressions, processMcDstarPiWithDecayTypeCheck, "Process MC with decay type checks for DstarPi", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCandidateCreatorB0Reduced>(cfgc),
                      adaptAnalysisTask<HfCandidateCreatorB0ReducedExpressions>(cfgc)};
}
