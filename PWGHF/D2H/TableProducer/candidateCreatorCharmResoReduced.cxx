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

/// \file candidateCreatorCharmResoReduced.cxx
/// \brief Reconstruction of Resonance candidates
///
/// \author Luca Aglietta <luca.aglietta@cern.ch>, Università degli Studi di Torino
#include <algorithm>

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

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace full
{
// Candidate
DECLARE_SOA_COLUMN(M, m, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(P, p, float);
// Daughters: Prong0: D, Prong1: V0
DECLARE_SOA_COLUMN(PtProng0, ptProng0, float);
DECLARE_SOA_COLUMN(MProng0, mProng0, float);
DECLARE_SOA_COLUMN(PtProng1, ptProng1, float);
DECLARE_SOA_COLUMN(MProng1, mProng1, float);
} //namespace full
} //namespace aod

// put the arguments into the table
DECLARE_SOA_TABLE(HfCandCharmReso, "AOD", "HFCANDCHARMRESO",
                full::M,
                full::Pt,
                full::P,
                full::MProng0,
                full::PtProng0,
                full::MProng1,
                full::PtProng1,
                hf_cand::PxProng0,
                hf_cand::PyProng0,
                hf_cand::PzProng0,
                hf_cand::PxProng1,
                hf_cand::PyProng1,
                hf_cand::PzProng1
                );
/// Reconstruction of B+ candidates
struct HfCandidateCreatorBplusReduced {
//Produces: Tables with resonance info
Produces<aod::HfCandCharmReso> rowCandidateReso;    

// Configurables
Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
Configurable<double> invMassWindowDV0{"invMassWindowDV0", 0.5, "invariant-mass window for DV0 pair selection (GeV/c2)"};
  
//Preslicing of D candidatess and V0s based on collisionId
Preslice<aod::HfRedD> candsDPerCollision = hf_track_index_reduced::hfRedCollisionId;
Preslice<aod::HfRedVzeros> candsV0PerCollision = hf_track_index_reduced::hfRedCollisionId;

//Useful constants
double massK0{0.};
double massLambda{0.};
double massDplus{0.};
double massDstar{0.};
double massReso{0.};
double bz{0.};

//Histogram registry: if task make it with a THNsparse with all variables you want to save
HistogramRegistry registry{"registry"};

void init(InitContext const&)
{
    // histograms
    registry.add("hMassBplusToD0Pi", "2-prong candidates;inv. mass (B^{+} #rightarrow #overline{D^{0}}#pi^{#plus} #rightarrow #pi^{#minus}K^{#plus}#pi^{#plus}) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 3., 8.}}});
    registry.add("hEvents", "Events;;entries", HistType::kTH1F, {{1, 0.5, 1.5}});

    // invariant-mass window cut
    massK0=o2::constants::physics::MassK0Short;
    massLambda=o2::constants::physics::MassLambda;
    massDplus=o2::constants::physics::MassDPlus;
    massDstar=o2::constants::physics::MassDStar;
    massReso=2.537;
}

/// Pion selection (D0 Pi <-- B+)
  /// \param trackPion is a track with the pion hypothesis
  /// \param trackParCovPion is the track parametrisation of the pion
  /// \param dcaPion is the 2-D array with track DCAs of the pion
  /// \param track0 is prong0 of selected D0 candidate
  /// \param track1 is prong1 of selected D0 candidate
  /// \param candD0 is the D0 candidate
  /// \return true if trackPion passes all cuts
  template <typename T1, typename T2>
  bool isV0Selected(const T1& preselD, const T2& PreselV0)
  {
    std::array<int, 3> dDaughtersIDs = {preselD.prong0Id(), preselD.prong1Id(), preselD.prong2Id()};
    // reject VOs that share daughters with D
    if (std::find(dDaughtersIDs.begin(), dDaughtersIDs.end(), PreselV0.prong0Id()) != dDaughtersIDs.end() || std::find(dDaughtersIDs.begin(), dDaughtersIDs.end(), PreselV0.prong1Id()) != dDaughtersIDs.end() ) {
      return false;
    }
    return true;
  }

//in process function nested loop and reject resonance candidates with track duplicates
void process(aod::HfRedCollisions const& collisions,
            aod::HfRedD const& candsD,
            aod::HfRedVzeros const& candsV0,
            aod::HfOrigColCounts const& collisionsCounter,
            aod::HfCandBpConfigs const& configs)
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

        auto candsDThisColl = candsD.sliceBy(candsDPerCollision, thisCollId);
        for (const auto& candD : candsDThisColl) {
            std::array<float, 3> pVecD = {candD.px(), candD.py(), candD.pz()};

            auto candV0ThisCollision = candsV0.sliceBy(candsV0PerCollision, thisCollId);
            for (const auto& candV0 : candV0ThisCollision) {
            std::array<float, 3> pVecV0 = {candV0.px(), candV0.py(), candV0.pz()};

            // compute invariant mass square and apply selection
            auto invMass2DV0 = RecoDecay::m2(std::array{pVecD, pVecV0}, std::array{massDplus, massK0});
            auto invMass2DV0Min = (massReso - invMassWindowDV0)*(massReso - invMassWindowDV0)
            auto invMass2DV0Max = (massReso + invMassWindowDV0)*(massReso + invMassWindowDV0)
            if ((invMass2DV0 < invMass2DV0Min) || (invMass2DV0 > invMass2DV0Max)) {
                continue;
            }
            // ---------------------------------

            registry.fill(HIST("hMassBplusToD0Pi"), std::sqrt(invMass2D0Pi));

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


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
return WorkflowSpec{adaptAnalysisTask<HfCandidateCreatorBplusReduced>(cfgc),
                    adaptAnalysisTask<HfCandidateCreatorBplusReducedExpressions>(cfgc)};
}
