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

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// Reconstruction of B+ candidates
struct HfCandidateCreatorBplusReduced {
  Produces<aod::HfCandBplusBase> rowCandidateBase; // table defined in CandidateReconstructionTables.h
  Produces<aod::HfRedBplusProngs> rowCandidateProngs; // table defined in ReducedDataModel.h

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
  // variable that will store the value of invMassWindowD0Pi (defined in dataCreatorD0PiReduced.cxx)
  float myInvMassWindowD0Pi{1.};

  double massPi{0.};
  double massD0{0.};
  double massBplus{0.};
  double bz{0.};

  // Fitter for B vertex (2-prong vertex filter)
  o2::vertexing::DCAFitterN<2> df2;

  Preslice<soa::Join<aod::HfRed2Prongs, aod::HfRed2ProngsCov>> candsDPerCollision = hf_track_index_reduced::hfRedCollisionId;
  Preslice<soa::Join<aod::HfRedTracks, aod::HfRedTracksCov>> tracksPionPerCollision = hf_track_index_reduced::hfRedCollisionId;

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    // histograms
    registry.add("hMassBplusToD0Pi", "2-prong candidates;inv. mass (B^{+} #rightarrow #overline{D^{0}}#pi^{#plus} #rightarrow #pi^{#minus}K^{#plus}#pi^{#plus}) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 3., 8.}}});
    registry.add("hCovPVXX", "2-prong candidates;XX element of cov. matrix of prim. vtx. position (cm^{2});entries", {HistType::kTH1F, {{100, 0., 1.e-4}}});
    registry.add("hCovSVXX", "2-prong candidates;XX element of cov. matrix of sec. vtx. position (cm^{2});entries", {HistType::kTH1F, {{100, 0., 0.2}}});
    registry.add("hEvents", "Events;;entries", HistType::kTH1F, {{1, 0.5, 1.5}});

    // Initialize fitter
    df2.setPropagateToPCA(propagateToPCA);
    df2.setMaxR(maxR);
    df2.setMaxDZIni(maxDZIni);
    df2.setMinParamChange(minParamChange);
    df2.setMinRelChi2Change(minRelChi2Change);
    df2.setUseAbsDCA(useAbsDCA);
    df2.setWeightedFinalPCA(useWeightedFinalPCA);

    // invariant-mass window cut
    massPi = o2::analysis::pdg::MassPiPlus;
    massD0 = o2::analysis::pdg::MassD0;
    massBplus = o2::analysis::pdg::MassBPlus;
  }

  void process(aod::HfRedCollisions const& collisions,
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
    double invMass2D0PiMin = (massBplus - myInvMassWindowD0Pi + invMassWindowD0PiTolerance) * (massBplus - myInvMassWindowD0Pi + invMassWindowD0PiTolerance);
    double invMass2D0PiMax = (massBplus + myInvMassWindowD0Pi - invMassWindowD0PiTolerance) * (massBplus + myInvMassWindowD0Pi - invMassWindowD0PiTolerance);

    for (const auto& collisionCounter : collisionsCounter) {
      registry.fill(HIST("hEvents"), 1, collisionCounter.originalCollisionCount());
    }

    static int ncol = 0;

    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto primaryVertex = getPrimaryVertex(collision);
      auto covMatrixPV = primaryVertex.getCov();

      if (ncol % 10000 == 0) {
        LOG(debug) << ncol << " collisions parsed";
      }
      ncol++;

      // Set the magnetic field from ccdb
      bz = collision.bz();
      df2.setBz(bz);

      auto candsDThisColl = candsD.sliceBy(candsDPerCollision, thisCollId);
      for (const auto& candD0 : candsDThisColl) {
        auto trackParCovD = getTrackParCov(candD0);
        std::array<float, 3> pVecD0 = {candD0.px(), candD0.py(), candD0.pz()};

        auto tracksPionThisCollision = tracksPion.sliceBy(tracksPionPerCollision, thisCollId);
        for (const auto& trackPion : tracksPionThisCollision) {
          auto trackParCovPi = getTrackParCov(trackPion);
          std::array<float, 3> pVecPion = {trackPion.px(), trackPion.py(), trackPion.pz()};

          // compute invariant mass square and apply selection
          auto invMass2D0Pi = RecoDecay::m2(std::array{pVecD0, pVecPion}, std::array{massD0, massPi});
          if ((invMass2D0Pi < invMass2D0PiMin) || (invMass2D0Pi > invMass2D0PiMax)) {
            continue;
          }
          // ---------------------------------
          // reconstruct the 2-prong B+ vertex
          if (df2.process(trackParCovD, trackParCovPi) == 0) {
            continue;
          }
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

          int hfFlag = BIT(hf_cand_bplus::DecayType::BplusToD0Pi);

          // fill the candidate table for the B+ here:
          rowCandidateBase(thisCollId,
                           collision.posX(), collision.posY(), collision.posZ(),
                           secondaryVertexBplus[0], secondaryVertexBplus[1], secondaryVertexBplus[2],
                           errorDecayLength, errorDecayLengthXY,
                           chi2PCA,
                           pVecD0[0], pVecD0[1], pVecD0[2],
                           pVecPion[0], pVecPion[1], pVecPion[2],
                           dcaD0.getY(), dcaPion.getY(),
                           std::sqrt(dcaD0.getSigmaY2()), std::sqrt(dcaPion.getSigmaY2()),
                           hfFlag);

          rowCandidateProngs(candD0.globalIndex(), trackPion.globalIndex());
        } // pi loop
      }   // D0 loop
    }     // collision loop
  }       // process
};        // struct

/// Extends the table base with expression columns and performs MC matching.
struct HfCandidateCreatorBplusReducedExpressions {
  Spawns<aod::HfCandBplusExt> rowCandidateBPlus;
  Spawns<aod::HfRedTracksExt> rowTracksExt;
  Produces<aod::HfMcRecRedBps> rowBplusMcRec;

  void processMc(HfMcRecRedD0Pis const& rowsD0PiMcRec, HfRedBplusProngs const& candsBplus)
  {
    for (const auto& candBplus : candsBplus) {
      for (const auto& rowD0PiMcRec : rowsD0PiMcRec) {
        if ((rowD0PiMcRec.prong0Id() != candBplus.prong0Id()) || (rowD0PiMcRec.prong1Id() != candBplus.prong1Id())) {
          continue;
        }
        rowBplusMcRec(rowD0PiMcRec.flagMcMatchRec(), rowD0PiMcRec.ptMother());
      }
    }
  }
  PROCESS_SWITCH(HfCandidateCreatorBplusReducedExpressions, processMc, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCandidateCreatorBplusReduced>(cfgc),
                      adaptAnalysisTask<HfCandidateCreatorBplusReducedExpressions>(cfgc)};
}
