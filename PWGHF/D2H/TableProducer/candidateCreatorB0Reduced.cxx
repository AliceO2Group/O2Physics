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

#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "DCAFitter/DCAFitterN.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/V0.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/D2H/DataModel/ReducedDataModel.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// Reconstruction of B0 candidates
struct HfCandidateCreatorB0Reduced {
  Produces<aod::HfCandB0Base> rowCandidateBase; // table defined in CandidateReconstructionTables.h

  // vertexing
  Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCA{"useAbsDCA", true, "Minimise abs. distance rather than chi2"};
  Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any B0 is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations is chi2/chi2old > this"};
  // selection
  Configurable<double> invMassWindowB0{"invMassWindowB0", 0.3, "invariant-mass window for B0 candidates"};

  double massPi = RecoDecay::getMassPDG(kPiPlus);
  double massD = RecoDecay::getMassPDG(pdg::Code::kDMinus);
  double massB0 = RecoDecay::getMassPDG(pdg::Code::kB0);
  double massDPi{0.};
  double bz{0.};

  // Fitter for B vertex (2-prong vertex filter)
  o2::vertexing::DCAFitterN<2> df2;

  Preslice<aod::HfCand3ProngReduced> candsDPerCollision = hf_track_index_reduced::hfReducedCollisionId;
  Preslice<aod::HfTracksReduced> tracksPionPerCollision = hf_track_index_reduced::hfReducedCollisionId;

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    // histograms
    registry.add("hMassB0ToDPi", "2-prong candidates;inv. mass (B^{0} #rightarrow D^{#minus}#pi^{#plus} #rightarrow #pi^{#minus}K^{#plus}#pi^{#minus}#pi^{#plus}) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 3., 8.}}});
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
  }

  void process(aod::HfReducedCollisions const& collisions,
               aod::HfCand3ProngReduced const& candsD,
               aod::HfTracksReduced const& tracksPion,
               aod::HfOriginalCollisionsCounter const& collisionsCounter)
  {
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
      for (const auto& candD : candsDThisColl) {
        auto trackParCovD = getTrackParCov(candD);
        std::array<float, 3> pVecD = {candD.px(), candD.py(), candD.pz()};

        auto tracksPionThisCollision = tracksPion.sliceBy(tracksPionPerCollision, thisCollId);
        for (const auto& trackPion : tracksPionThisCollision) {
          auto trackParCovPi = getTrackParCov(trackPion);
          std::array<float, 3> pVecPion = {trackPion.px(), trackPion.py(), trackPion.pz()};

          // ---------------------------------
          // reconstruct the 2-prong B0 vertex
          if (df2.process(trackParCovD, trackParCovPi) == 0) {
            continue;
          }
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

          // compute invariant
          massDPi = RecoDecay::m(array{pVecD, pVecPion}, array{massD, massPi});

          if (std::abs(massDPi - massB0) > invMassWindowB0) {
            continue;
          }
          registry.fill(HIST("hMassB0ToDPi"), massDPi);

          // compute impact parameters of D and Pi
          o2::dataformats::DCA dcaD;
          o2::dataformats::DCA dcaPion;
          trackParCovD.propagateToDCA(primaryVertex, bz, &dcaD);
          trackParCovPi.propagateToDCA(primaryVertex, bz, &dcaPion);

          // get uncertainty of the decay length
          double phi, theta;
          // getPointDirection modifies phi and theta
          getPointDirection(array{collision.posX(), collision.posY(), collision.posZ()}, secondaryVertexB0, phi, theta);
          auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
          auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));

          int hfFlag = BIT(hf_cand_b0::DecayType::B0ToDPi);

          // fill the candidate table for the B0 here:
          rowCandidateBase(thisCollId,
                           collision.posX(), collision.posY(), collision.posZ(),
                           secondaryVertexB0[0], secondaryVertexB0[1], secondaryVertexB0[2],
                           errorDecayLength, errorDecayLengthXY,
                           chi2PCA,
                           pVecD[0], pVecD[1], pVecD[2],
                           pVecPion[0], pVecPion[1], pVecPion[2],
                           dcaD.getY(), dcaPion.getY(),
                           std::sqrt(dcaD.getSigmaY2()), std::sqrt(dcaPion.getSigmaY2()),
                           candD.globalIndex(), trackPion.globalIndex(),
                           hfFlag);
        } // pi loop
      }   // D loop
    }     // collision loop
  }       // process
};        // struct

/// Extends the table base with expression columns and performs MC matching.
struct HfCandidateCreatorB0ReducedExpressions {
  Spawns<aod::HfCandB0Ext> rowCandidateB0;
  Produces<aod::HfB0McRecReduced> rowB0McRec;

  void processMc(HfDPiMcRecReduced const& rowsDPiMcRec)
  {
    for (const auto& candB0 : *rowCandidateB0) {
      for (const auto& rowDPiMcRec : rowsDPiMcRec) {
        if ((rowDPiMcRec.prong0Id() != candB0.prong0Id()) || (rowDPiMcRec.prong1Id() != candB0.prong1Id())) {
          continue;
        }
        rowB0McRec(rowDPiMcRec.flagMcMatchRec(), rowDPiMcRec.originMcRec(), rowDPiMcRec.debugMcRec(), rowDPiMcRec.ptMother());
      }
    }
  }
  PROCESS_SWITCH(HfCandidateCreatorB0ReducedExpressions, processMc, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCandidateCreatorB0Reduced>(cfgc),
                      adaptAnalysisTask<HfCandidateCreatorB0ReducedExpressions>(cfgc)};
}
