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

/// \file Lambdacascadecorrelation.cxx
/// \brief Correlation-balance functions of multistrange baryons
/// \author Oveis Sheibani <oveis.sheibani@cern.ch>
//
// o2-linter: disable=name/workflow-file (filename retained for back-compat with existing alienv install and user JSON configs; file contains three structs LambdaCascadeProducer / LambdaTracksExtProducer / LambdaXiCorrelation so no single struct name can be the file name)

#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/Utils/inelGt.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>
#include <TList.h>
#include <TObject.h>
#include <TPDGCode.h>
#include <TString.h>
#include <TTree.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <string_view>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace o2::constants::math;
using namespace o2::soa;

// use parameters + cov mat non-propagated, aux info + (extension propagated)
using FullTracksExt = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov>;
using FullTracksExtIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU>;
using FullTracksExtWithPID = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr>;
using FullTracksExtIUWithPID = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr>;

// [Phase 12d → reverted Phase 16b] Data-model declarations inlined back into
// the .cxx so the task is contained in a single source file (ALICE submission
// policy). Originally extracted to a sibling header in Phase 12d.
//
// Phase log:
//   Phase 8  — CascadeFlags::IsTrueCascade (MC-truth purity tag)
//   Phase 8  — LambdaTracks::DcaV0ToPV / V0Radius / Pos+NegItsNCls /
//              PassesPrimaryTopo
//   Phase 9  — CascadeFlags::IsItsTracked  (ITS strangeness tracking)
//   Phase 9  — LambdaTracks::Pos+NegItsClusterMap
//   Phase 10 — LambdaTracks::LProper, Pos+NegDcaXY
//   Phase 14 — CascadeFlags::CascCutBits, LambdaTracks::CutBits +
//              raw cut inputs (tpcNSigma, mK0Short, qtArm, alphaArm, cTau)

namespace o2::aod
{
namespace cascadeflags
{
DECLARE_SOA_COLUMN(IsSelected, isSelected, int); //~!
// [Phase 8] MC-truth purity flag for the cascade row.
// True iff the cascade has a matched MC particle whose pdgCode is
// ±3312 (Ξ) or ±3334 (Ω) AND that MC particle is physical-primary.
// False on data and for combinatorial fakes (typical failure: a primary
// Λ paired with a stray pion that happened to satisfy the cascade-vertex
// fit). Lets the correlator filter out non-Ξ/Ω cascades when MC truth
// is available (LabeledCascades branch).
DECLARE_SOA_COLUMN(IsTrueCascade, isTrueCascade, bool);
// [Phase 9] ITS-strangeness-tracking flag for the cascade row.
// True iff the cascade has a matching aod::AssignedTrackedCascades row,
// i.e. ITS reconstructed a track segment for the parent Ξ⁻/Ω⁻ in
// the IB layers (cτ_Ξ ≈ 4.91 cm — the Ξ traverses the inner barrel
// before decaying often enough to be tracked). This is the closest
// thing ALICE has to "ITS confirms this is a real cascade" on data.
// Available on Run 3 only; ITS strangeness tracking is not produced
// for Run 2 AODs.
DECLARE_SOA_COLUMN(IsItsTracked, isItsTracked, bool);
// [Phase 14] Per-cascade bitmask: each bit = "passed cut N". Lets the
// downstream tree consumer reconstruct any cascade selection offline
// without re-running the workflow. Bit semantics defined by enum
// CascCutBit below.
DECLARE_SOA_COLUMN(CascCutBits, cascCutBits, uint32_t);
} // namespace cascadeflags
DECLARE_SOA_TABLE(CascadeFlags, "AOD", "CASCADEFLAGS", //!
                  cascadeflags::IsSelected,
                  cascadeflags::IsTrueCascade,
                  cascadeflags::IsItsTracked,
                  cascadeflags::CascCutBits);
using CascDataExtSelected = soa::Join<CascDataExt, CascadeFlags>;
} // namespace o2::aod

namespace o2::aod
{
namespace lambdacollision
{
DECLARE_SOA_COLUMN(Cent, cent, float);
DECLARE_SOA_COLUMN(Mult, mult, float);
DECLARE_SOA_COLUMN(RefCollId, refCollId, int64_t);
} // namespace lambdacollision
DECLARE_SOA_TABLE(LambdaCollisions, "AOD", "LAMBDACOLS", o2::soa::Index<>,
                  lambdacollision::Cent,
                  lambdacollision::Mult,
                  lambdacollision::RefCollId,
                  aod::collision::PosX,
                  aod::collision::PosY,
                  aod::collision::PosZ);
using LambdaCollision = LambdaCollisions::iterator;

namespace lambdamcgencollision
{
DECLARE_SOA_COLUMN(RefMcCollId, refMcCollId, int64_t); // original McCollision global index
}
DECLARE_SOA_TABLE(LambdaMcGenCollisions, "AOD", "LMCGENCOLS", o2::soa::Index<>,
                  lambdacollision::Cent,
                  lambdacollision::Mult,
                  lambdamcgencollision::RefMcCollId,
                  o2::aod::mccollision::PosX,
                  o2::aod::mccollision::PosY,
                  o2::aod::mccollision::PosZ);
using LambdaMcGenCollision = LambdaMcGenCollisions::iterator;

namespace lambdatrack
{
DECLARE_SOA_INDEX_COLUMN(LambdaCollision, lambdaCollision);
DECLARE_SOA_COLUMN(Px, px, float);
DECLARE_SOA_COLUMN(Py, py, float);
DECLARE_SOA_COLUMN(Pz, pz, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Rap, rap, float);
DECLARE_SOA_COLUMN(Mass, mass, float);
DECLARE_SOA_COLUMN(PosTrackId, posTrackId, int64_t);
DECLARE_SOA_COLUMN(NegTrackId, negTrackId, int64_t);
DECLARE_SOA_COLUMN(CosPA, cosPA, float);
DECLARE_SOA_COLUMN(DcaDau, dcaDau, float);
DECLARE_SOA_COLUMN(V0Type, v0Type, int8_t);
DECLARE_SOA_COLUMN(V0PrmScd, v0PrmScd, int8_t);
DECLARE_SOA_COLUMN(CorrFact, corrFact, float);
// [Phase 7] Mother PDG of the truth-matched MC Λ. 0 on data (no truth).
// In MC: 0 if the standalone Λ has no MC particle; else the parent PDG.
// Lets downstream consumers tag feed-down (e.g., motherPdg == ±3312 for
// Ξ-feeddown, ±3334 for Ω-feeddown, ±3212 for Σ⁰, etc.) without rerunning
// MC truth-matching themselves.
DECLARE_SOA_COLUMN(MotherPdg, motherPdg, int);
// [Phase 8] Per-V0 topology snapshot stored on the row so:
//   (i)  the optional Λ TTree carries the variables needed to do a
//        pT-differential primary-fraction template fit on data,
//   (ii) the partition can gate on a topology-only "looks primary" flag
//        on data (where v0PrmScd is uninformative).
DECLARE_SOA_COLUMN(DcaV0ToPV, dcaV0ToPV, float);
DECLARE_SOA_COLUMN(V0Radius, v0Radius, float);
DECLARE_SOA_COLUMN(PosItsNCls, posItsNCls, int8_t);
DECLARE_SOA_COLUMN(NegItsNCls, negItsNCls, int8_t);
// [Phase 8] Topology-only "looks primary" flag, computed in the V0 loop
// from the primCfg bundle:
//    dcav0topv      < cPrimMaxDcaV0ToPv
// && v0cosPA(PV)    > cPrimMinV0CosPA
// && v0radius       < cPrimMaxV0Radius
// && |dcapostopv|, |dcanegtopv| < cPrimMaxDauDcaToPv
// && posItsNCls, negItsNCls    >= cPrimMinDauItsNCls
// On data this is the only handle on primary-Λ-ness, so the trigger
// partition AND's it. On MC the partition AND's it WITH the truth bit
// (v0PrmScd==kPrimary), and the Λ TTree exposes the components so an
// external macro can build templates per motherPdg and refine the cuts.
DECLARE_SOA_COLUMN(PassesPrimaryTopo, passesPrimaryTopo, bool);
// [Phase 9] ITS hit-map per V0 daughter (uint8_t bitmask, one bit per
// of the 7 ITS layers). Lets a downstream macro do a geometric
// consistency check: a daughter from a V0 vertex at radius r > r_layer
// CANNOT have a hit on that layer — if it does, the V0 is mislocated
// or the "daughter" is a primary track misassigned. Tighter than the
// itsNCls integer count, free at this stage (just propagated through).
DECLARE_SOA_COLUMN(PosItsClusterMap, posItsClusterMap, uint8_t);
DECLARE_SOA_COLUMN(NegItsClusterMap, negItsClusterMap, uint8_t);
// [Phase 10] Pseudo-proper transverse decay length L_proper = L_xy * M_Λ / pT.
// Computed treating the V0 vertex distance from PV as L_xy. For primary
// Λ this is the actual proper-time projection (exponential with cτ=7.89 cm);
// for feed-down Λ it OVER-estimates because L_xy includes the parent
// flight (e.g. Ξ cτ ≈ 4.91 cm). The discriminator is the upper tail.
DECLARE_SOA_COLUMN(LProper, lProper, float);
// [Phase 10] Per-daughter DCA-XY-to-PV (signed). Already-cut on at production
// (cMinDcaProtonToPV / cMinDcaPionToPV) but not previously exposed downstream.
// Needed for the optional MVA / template-fit refinement and for inspecting
// the low-pT regime where these tighten primary-Λ purity.
DECLARE_SOA_COLUMN(PosDcaXY, posDcaXY, float);
DECLARE_SOA_COLUMN(NegDcaXY, negDcaXY, float);
// [Phase 14] Per-V0 cut bitmask: each bit = "passed cut N". Computed
// unconditionally for every accepted V0 (and every V0 when the
// diagnostic mode cFillLambdaTreeAllCandidates is on). Lets a downstream
// macro reconstruct any cut combination offline without re-running.
// Bit semantics defined by enum LambdaCutBit below.
DECLARE_SOA_COLUMN(CutBits, cutBits, uint32_t);
// [Phase 14] Raw cut-input values exposed on the Λ row so post-hoc cuts
// can be tightened/loosened without re-running. tpcNSigma values are
// matched to the v0Type hypothesis (proton-leg vs pion-leg).
DECLARE_SOA_COLUMN(TpcNSigmaPosPr, tpcNSigmaPosPr, float);  // proton hypothesis on positive daughter
DECLARE_SOA_COLUMN(TpcNSigmaNegPi, tpcNSigmaNegPi, float);  // pion   hypothesis on negative daughter
DECLARE_SOA_COLUMN(TpcNSigmaPosPi, tpcNSigmaPosPi, float);  // pion   hypothesis on positive daughter
DECLARE_SOA_COLUMN(TpcNSigmaNegPr, tpcNSigmaNegPr, float);  // proton hypothesis on negative daughter
DECLARE_SOA_COLUMN(MK0Short, mK0Short, float);              // K0s mass hypothesis (for sideband studies)
DECLARE_SOA_COLUMN(QtArm, qtArm, float);                    // Armenteros qT
DECLARE_SOA_COLUMN(AlphaArm, alphaArm, float);              // Armenteros α
DECLARE_SOA_COLUMN(CTau, cTau, float);                      // proper time × c (ctau)
} // namespace lambdatrack
DECLARE_SOA_TABLE(LambdaTracks, "AOD", "LAMBDATRACKS", o2::soa::Index<>,
                  lambdatrack::LambdaCollisionId,
                  lambdatrack::Px,
                  lambdatrack::Py,
                  lambdatrack::Pz,
                  lambdatrack::Pt,
                  lambdatrack::Eta,
                  lambdatrack::Phi,
                  lambdatrack::Rap,
                  lambdatrack::Mass,
                  lambdatrack::PosTrackId,
                  lambdatrack::NegTrackId,
                  lambdatrack::CosPA,
                  lambdatrack::DcaDau,
                  lambdatrack::V0Type,
                  lambdatrack::V0PrmScd,
                  lambdatrack::CorrFact,
                  lambdatrack::MotherPdg,
                  lambdatrack::DcaV0ToPV,
                  lambdatrack::V0Radius,
                  lambdatrack::PosItsNCls,
                  lambdatrack::NegItsNCls,
                  lambdatrack::PassesPrimaryTopo,
                  lambdatrack::PosItsClusterMap,
                  lambdatrack::NegItsClusterMap,
                  lambdatrack::LProper,
                  lambdatrack::PosDcaXY,
                  lambdatrack::NegDcaXY,
                  lambdatrack::CutBits,
                  lambdatrack::TpcNSigmaPosPr,
                  lambdatrack::TpcNSigmaNegPi,
                  lambdatrack::TpcNSigmaPosPi,
                  lambdatrack::TpcNSigmaNegPr,
                  lambdatrack::MK0Short,
                  lambdatrack::QtArm,
                  lambdatrack::AlphaArm,
                  lambdatrack::CTau);
using LambdaTrack = LambdaTracks::iterator;

namespace lambdatrackext
{
DECLARE_SOA_COLUMN(LambdaSharingDaughter, lambdaSharingDaughter, bool);
DECLARE_SOA_COLUMN(LambdaSharingDauIds, lambdaSharingDauIds, std::vector<int64_t>);
DECLARE_SOA_COLUMN(TrueLambdaFlag, trueLambdaFlag, bool);
} // namespace lambdatrackext
DECLARE_SOA_TABLE(LambdaTracksExt, "AOD", "LAMBDATRACKSEXT",
                  lambdatrackext::LambdaSharingDaughter,
                  lambdatrackext::LambdaSharingDauIds,
                  lambdatrackext::TrueLambdaFlag);

using LambdaTrackExt = LambdaTracksExt::iterator;

namespace lambdamcgentrack
{
DECLARE_SOA_INDEX_COLUMN(LambdaMcGenCollision, lambdaMcGenCollision);
}
DECLARE_SOA_TABLE(LambdaMcGenTracks, "AOD", "LMCGENTRACKS", o2::soa::Index<>,
                  lambdamcgentrack::LambdaMcGenCollisionId,
                  o2::aod::mcparticle::Px,
                  o2::aod::mcparticle::Py,
                  o2::aod::mcparticle::Pz,
                  lambdatrack::Pt,
                  lambdatrack::Eta,
                  lambdatrack::Phi,
                  lambdatrack::Rap,
                  lambdatrack::Mass,
                  lambdatrack::PosTrackId,
                  lambdatrack::NegTrackId,
                  lambdatrack::V0Type,
                  lambdatrack::CosPA,
                  lambdatrack::DcaDau,
                  lambdatrack::V0PrmScd,
                  lambdatrack::CorrFact);
using LambdaMcGenTrack = LambdaMcGenTracks::iterator;

} // namespace o2::aod

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::PVMults>;
using MyCollisionsMult = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;
using MyCascades = soa::Filtered<aod::CascDataExtSelected>;
using LabeledCascades = soa::Join<aod::CascDataExt, aod::McCascLabels>;

// [Phase 12d] All Λ-side data-model declarations live in
// Lambdacascadecorrelation_DataModel.h, included above.

enum CollisionLabels {
  kTotColBeforeHasMcCollision = 1,
  kTotCol,
  kPassSelCol
};

enum TrackLabels {
  kTracksBeforeHasMcParticle = 1,
  kAllV0Tracks,
  kV0KShortMassRej,
  kNotLambdaNotAntiLambda,
  kV0IsBothLambdaAntiLambda,
  kNotLambdaAfterSel,
  kV0IsLambdaOrAntiLambda,
  kPassV0DauTrackSel,
  kPassV0KinCuts,
  kPassV0TopoSel,
  kAllSelPassed,
  kPrimaryLambda,
  kSecondaryLambda,
  kLambdaDauNotMcParticle,
  kLambdaNotPrPiMinus,
  kAntiLambdaNotAntiPrPiPlus,
  kPassTrueLambdaSel,
  kEffCorrPtCent,
  kEffCorrPtRapCent,
  kNoEffCorr,
  kPFCorrPtCent,
  kPFCorrPtRapCent,
  kNoPFCorr,
  kGenTotAccLambda,
  kGenLambdaNoDau,
  kGenLambdaToPrPi
};

enum CentEstType {
  kCentFT0M = 0,
  kCentFV0A
};

enum RunType {
  kRun3 = 0,
  kRun2
};

enum ParticleType {
  kLambda = 0,
  kAntiLambda
};

enum ParticlePairType {
  kLambdaAntiLambda = 0,
  kLambdaLambda,
  kAntiLambdaAntiLambda
};

enum ShareDauLambda {
  kUniqueLambda = 0,
  kLambdaShareDau
};

enum RecGenType {
  kRec = 0,
  kGen
};

enum DMCType {
  kData = 0,
  kMC
};

// [Phase 12a] Centralised constants — replace scattered magic literals.
// PDG codes are taken from ROOT's TPDGCode.h enum (PDG_t) so the o2 linter
// is happy and the source of truth is the upstream framework header.
namespace lcorr_const
{
// PDG codes for the cascades we tag/filter on. Drawn from ::kXiMinus etc.
// in <TPDGCode.h> so the magic-literal-detector won't complain.
constexpr int kLambdaPdg     = ::kLambda0;     // 3122, Λ (uds)
constexpr int kXiMinusPdg    = ::kXiMinus;     // 3312, Ξ⁻ (dss)
constexpr int kOmegaMinusPdg = ::kOmegaMinus;  // 3334, Ω⁻ (sss)
constexpr int kSigma0Pdg     = ::kSigma0;      // 3212, Σ⁰ (uds, EM-decay → Λγ)

// ITS Inner-Barrel layer mask (Layers 0,1,2 in itsClusterMap bitfield).
// Used by the Phase 10 "≥1 daughter has an IB hit" primary-Λ requirement.
constexpr uint8_t kItsIBMask = 0x07;

// [Phase 16j] Cascade species-selection flag values written into the
// cascadeflags::IsSelected column. Promoted to namespace scope so both
// LambdaCascadeProducer (sets them) and LambdaXiCorrelation (consumes
// them) can refer to the same names.
constexpr int kFlagRejected   = 0; // rejected by processCandidate cut chain
constexpr int kFlagXiOnly     = 1; // bachelor passes pion-PID only → Ξ-eligible
constexpr int kFlagXiAndOmega = 2; // bachelor passes both pion AND kaon PID → both
constexpr int kFlagOmegaOnly  = 3; // bachelor passes kaon-PID only → Ω-eligible

// [Phase 16j] cVetoMode values for the auto-correlation veto policy.
constexpr int kVetoModeOff    = 0; // no veto
constexpr int kVetoModeStrict = 1; // veto only when BOTH daughters shared
constexpr int kVetoModeLoose  = 2; // veto when EITHER daughter shared

// [Phase 16j] cItsTrackMode values.
constexpr int kItsTrackModeOff      = 0; // ignore ITS-tracking flag
constexpr int kItsTrackModeRequired = 1; // require isItsTracked == true
constexpr int kItsTrackModeRescue   = 2; // accept isItsTracked OR truth-match
} // namespace lcorr_const

// [Phase 14] Per-Λ cut-bit enum. Bit i in lambdatrack::cutBits = 1 iff the
// candidate passed cut i. Computed UNCONDITIONALLY for every V0 row
// emitted into the LambdaTracks table — i.e. early-return short-circuits
// in selV0Particle do NOT hide later-stage results. This lets an offline
// macro reconstruct any cut combination from the tree without re-running.
enum LambdaCutBit : uint32_t {
  kCutMassWindow      = 0,   // cMinV0Mass < mLambda < cMaxV0Mass (matched to v0Type)
  kCutDauPid          = 1,   // |nσ(p)| < cTpcNsigmaCut AND |nσ(π)| < cut
  kCutDauTrackQual    = 2,   // selTrack on both daughters
  kCutDauDcaToPV      = 3,   // proton-leg DCA > min AND pion-leg DCA > min
  kCutKinematic       = 4,   // pT and |y or η| within cMinV0Pt..cMaxV0Pt and cMaxV0Rap
  kCutDcaV0Dau        = 5,   // dcaV0daughters in [min,max]
  kCutDcaV0ToPV       = 6,   // dcav0topv in [min,max]
  kCutV0Radius        = 7,   // v0radius in [min,max]
  kCutCtau            = 8,   // ctau in [min,max]
  kCutCosPA           = 9,   // v0cosPA > cMinV0CosPA
  kCutK0sRej          = 10,  // K0s mass-window rejection passed (or flag off)
  kCutAmbiguousVeto   = 11,  // !hasAmbiguousDaughters or veto disabled
  kCutMcTrueLambda    = 12,  // MC: selTrueMcRecLambda (always true on data)
  kCutPhase10Prim     = 13,  // passesPrimaryTopo (set in V0 loop)
  kCutBitMax          = 14
};

// [Phase 14] Per-cascade cut-bit enum. Bit i in cascadeflags::cascCutBits
// = 1 iff the candidate passed cut i. Mirrors the structure of
// LambdaCutBit but for the cascade selection chain in processCandidate.
enum CascCutBit : uint32_t {
  kCascCutTpcRowsV0Dau   = 0,   // V0-daughter tpcNClsCrossedRows
  kCascCutTpcRowsBach    = 1,   // bachelor tpcNClsCrossedRows
  kCascCutItsClsV0Dau    = 2,   // V0-daughter itsNCls
  kCascCutItsClsBach     = 3,   // bachelor itsNCls
  kCascCutItsChi2        = 4,   // ITS chi²/cluster (V0 + bachelor)
  kCascCutTpcChi2        = 5,   // TPC chi²/cluster (V0 + bachelor)
  kCascCutCascPt         = 6,   // casc.pt() > cMinCascPt
  kCascCutTopology       = 7,   // v0radius/cascradius/cosPA/dcav0topv/v0mass
  kCascCutRadiusOrder    = 8,   // cascRadius < v0Radius (consistency)
  kCascCutTrackEta       = 9,   // |η| of V0 dau + bachelor
  kCascCutCascEta        = 10,  // |η| of cascade
  kCascCutTpcNSigPr      = 11,  // proton-leg TPC nσ
  kCascCutTpcNSigPi      = 12,  // pion-leg TPC nσ
  kCascCutBachPidXi      = 13,  // bachelor pion nσ (Ξ hypothesis)
  kCascCutBachPidOm      = 14,  // bachelor kaon nσ (Ω hypothesis)
  kCascCutCompetingMass  = 15,  // competing-mass cut (Ω vs Ξ)
  kCascCutMcTrueXiOmega  = 16,  // MC: pdgCode ∈ {±3312,±3334} && isPhysPrim
  kCascCutItsTracked     = 17,  // matched in aod::AssignedTrackedCascades
  kCascCutBitMax         = 18
};

enum CorrHistDim {
  OneDimCorr = 1,
  TwoDimCorr,
  ThreeDimCorr
};

enum PrmScdType {
  kPrimary = 0,
  kSecondary
};

enum PrmScdPairType {
  kPP = 0,
  kPS,
  kSP,
  kSS
};

// =============================================================================
// [P1][R1] Shared event-selection helpers
//
// Both LambdaTableProducer (selCollision) and CascadeSelector (eventSelection)
// historically applied their own, slightly different, sets of event cuts. That
// produced collision populations that disagreed in subtle ways (e.g. LTP did
// not veto same-bunch pileup by default, CSEL did) and broke the downstream
// Λ-Ξ correlator's normalisation.
//
// This namespace provides:
//   - a POD `EventCuts` describing every cut either producer might apply,
//   - a templated `applyEventSelection<RunType>(col, cuts, centValue, &reason)`
//     that evaluates only the cuts whose `use*` flag is true,
//   - a `logEventCuts(tag, cuts)` that prints the cut configuration in a
//     unified, single-line, grep-friendly format.
//
// Each producer builds an `EventCuts` from its own Configurables in init() and
// uses the shared function. Defaults are unchanged — nothing differs at
// runtime unless the user touches the configurables — but the duplicated logic
// is gone and the [EVENTSEL-LTP] / [EVENTSEL-CSEL] log lines now make the
// producer-mismatch visible side by side.
// =============================================================================
namespace lcorr_evsel
{

struct EventCuts {
  // Vertex Z. LTP's original semantics: posZ <= min || posZ >= max → reject
  // (open interval). CSEL's were |posZ| > max (closed). We adopt LTP's
  // convention here; CSEL fills minVtxZ=-maxVertexZ, maxVtxZ=+maxVertexZ
  // (the difference matters only for collisions exactly at the boundary,
  // which is statistically zero events).
  bool useVtxZ{false};
  float minVtxZ{-10.f}, maxVtxZ{10.f};

  // Trigger
  bool useSel8{false}; // Run3
  bool useInt7{false}; // Run2
  bool useSel7{false}; // Run2

  // Centrality range (LTP only by default)
  bool useCentRange{false};
  float minCent{0.f}, maxCent{100.f};

  // INEL multiplicity cut (CSEL only by default)
  // Keeps multNTracksPVeta1 > inelMin convention from CascadeSelector.
  bool useInel{false};
  int inelMin{0};

  // Run3 selection bits
  bool useTriggerTvx{false};
  bool useTfBorder{false};
  bool useItsRoBorder{false};
  bool useItsTpcVtx{false};
  bool useNoSameBunchPileup{false};
  bool useZVtxTimeDiff{false};
  bool useIsGoodITSLayers{false};
};

// Apply the configured cuts to a collision. Returns true if accepted; on
// rejection writes a static C-string naming the failing cut into *reason
// so callers can bin-fill, log, etc. centValue is supplied by the caller
// because LTP selects the centrality estimator before the call.
template <RunType run, typename C>
inline bool applyEventSelection(C const& col, EventCuts const& cuts,
                                float centValue, const char*& reason)
{
  reason = "ok";

  // Vertex Z (LTP open-interval semantics)
  if (cuts.useVtxZ) {
    if (col.posZ() <= cuts.minVtxZ || col.posZ() >= cuts.maxVtxZ) {
      reason = "VtxZ";
      return false;
    }
  }

  // Trigger — Run3
  if constexpr (run == kRun3) {
    if (cuts.useSel8 && !col.sel8()) {
      reason = "Sel8";
      return false;
    }
  } else {
    // Trigger — Run2
    if (cuts.useInt7 && !col.alias_bit(kINT7)) {
      reason = "Int7";
      return false;
    }
    if (cuts.useSel7 && !col.sel7()) {
      reason = "Sel7";
      return false;
    }
  }

  // Centrality range (estimator already chosen by caller → centValue passed in)
  if (cuts.useCentRange) {
    if (centValue <= cuts.minCent || centValue >= cuts.maxCent) {
      reason = "CentRange";
      return false;
    }
  }

  // Run3 selection bits + INEL
  if constexpr (run == kRun3) {
    if (cuts.useInel && col.multNTracksPVeta1() <= cuts.inelMin) {
      reason = "INEL";
      return false;
    }
    if (cuts.useTriggerTvx && !col.selection_bit(aod::evsel::kIsTriggerTVX)) {
      reason = "TriggerTVX";
      return false;
    }
    if (cuts.useTfBorder && !col.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      reason = "TFBorder";
      return false;
    }
    if (cuts.useItsRoBorder && !col.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      reason = "ITSROBorder";
      return false;
    }
    if (cuts.useItsTpcVtx && !col.selection_bit(aod::evsel::kIsVertexITSTPC)) {
      reason = "ItsTpcVtx";
      return false;
    }
    if (cuts.useNoSameBunchPileup && !col.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      reason = "NoSameBunchPileup";
      return false;
    }
    if (cuts.useZVtxTimeDiff && !col.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      reason = "ZVtxTimeDiff";
      return false;
    }
    if (cuts.useIsGoodITSLayers && !col.selection_bit(aod::evsel::kIsGoodITSLayersAll)) {
      reason = "GoodITSLayers";
      return false;
    }
  }

  return true;
}

// Single-line, grep-friendly summary of the cut configuration. tag is the
// producer label, e.g. "EVENTSEL-LTP" or "EVENTSEL-CSEL".
inline void logEventCuts(const char* tag, EventCuts const& cuts)
{
  LOGF(info,
       "[%s] vtxZ=%d[%.2f,%.2f] sel8=%d int7=%d sel7=%d "
       "centRange=%d[%.1f,%.1f] inel=%d>%d "
       "trigTvx=%d tfBorder=%d itsROBorder=%d itsTpcVtx=%d "
       "noBunchPileup=%d zVtxTimeDiff=%d goodITSLayers=%d",
       tag,
       (int)cuts.useVtxZ, (double)cuts.minVtxZ, (double)cuts.maxVtxZ,
       (int)cuts.useSel8, (int)cuts.useInt7, (int)cuts.useSel7,
       (int)cuts.useCentRange, (double)cuts.minCent, (double)cuts.maxCent,
       (int)cuts.useInel, cuts.inelMin,
       (int)cuts.useTriggerTvx, (int)cuts.useTfBorder,
       (int)cuts.useItsRoBorder, (int)cuts.useItsTpcVtx,
       (int)cuts.useNoSameBunchPileup, (int)cuts.useZVtxTimeDiff,
       (int)cuts.useIsGoodITSLayers);
}

} // namespace lcorr_evsel

// =============================================================================
// [Phase 4] LambdaCascadeProducer — merger of LambdaTableProducer (LTP) and
// the standalone CascadeSelector (CSEL).
//
// Why merged: event selection used to run twice (once in LTP, once in CSEL),
// with subtly different defaults, producing inconsistent collision populations
// in the downstream LambdaXiCorrelation. The merge guarantees a single
// `selCollision()` call per event drives BOTH the Lambda table production
// and the cascade flagging — they can no longer disagree.
//
// Output tables (unchanged from before):
//   - aod::LambdaCollisions, aod::LambdaTracks                     (Lambda)
//   - aod::LambdaMcGenCollisions, aod::LambdaMcGenTracks            (Lambda MC)
//   - aod::CascadeFlags                                             (Cascade)
//
// Histogram registries: kept distinct (`histos` for Lambda-side, `cascRegistry`
// for cascade-side) so every existing path in AnalysisResults.root stays
// byte-identical. Downstream macros / hyperloop wagons / plotting scripts
// require no changes.
// =============================================================================
struct LambdaCascadeProducer {

  // === Outputs ===
  // Lambda-side
  Produces<aod::LambdaCollisions> lambdaCollisionTable;
  Produces<aod::LambdaTracks> lambdaTrackTable;
  Produces<aod::LambdaMcGenCollisions> lambdaMCGenCollisionTable;
  Produces<aod::LambdaMcGenTracks> lambdaMCGenTrackTable;
  // Cascade-side
  Produces<aod::CascadeFlags> cascflags;

  // Collisions
  Configurable<int> cCentEstimator{"cCentEstimator", 0, "Centrality Estimator : 0-FT0M, 1-FV0A"};
  Configurable<float> cMinZVtx{"cMinZVtx", -10.0, "Min VtxZ cut"};
  Configurable<float> cMaxZVtx{"cMaxZVtx", 10.0, "Max VtxZ cut"};
  Configurable<float> cMinMult{"cMinMult", 0., "Minumum Multiplicity"};
  Configurable<float> cMaxMult{"cMaxMult", 100.0, "Maximum Multiplicity"};
  Configurable<bool> cSel8Trig{"cSel8Trig", true, "Sel8 (T0A + T0C) Selection Run3"};
  Configurable<bool> cInt7Trig{"cInt7Trig", false, "kINT7 MB Trigger"};
  Configurable<bool> cSel7Trig{"cSel7Trig", false, "Sel7 (V0A + V0C) Selection Run2"};
  Configurable<bool> cTriggerTvxSel{"cTriggerTvxSel", false, "Trigger Time and Vertex Selection"};
  Configurable<bool> cTFBorder{"cTFBorder", false, "Timeframe Border Selection"};
  Configurable<bool> cNoItsROBorder{"cNoItsROBorder", false, "No ITSRO Border Cut"};
  Configurable<bool> cItsTpcVtx{"cItsTpcVtx", false, "ITS+TPC Vertex Selection"};
  Configurable<bool> cPileupReject{"cPileupReject", false, "Pileup rejection"};
  Configurable<bool> cZVtxTimeDiff{"cZVtxTimeDiff", false, "z-vtx time diff selection"};
  Configurable<bool> cIsGoodITSLayers{"cIsGoodITSLayers", false, "Good ITS Layers All"};

  // Tracks
  Configurable<float> cTrackMinPt{"cTrackMinPt", 0.15, "p_{T} minimum"};
  Configurable<float> cTrackMaxPt{"cTrackMaxPt", 999.0, "p_{T} maximum"};
  Configurable<float> cTrackEtaCut{"cTrackEtaCut", 0.8, "Pseudorapidity cut"};
  Configurable<int> cMinTpcCrossedRows{"cMinTpcCrossedRows", 70, "TPC Min Crossed Rows"};
  Configurable<float> cMinTpcCROverCls{"cMinTpcCROverCls", 0.8, "Tpc Min Crossed Rows Over Findable Clusters"};
  Configurable<float> cMaxTpcSharedClusters{"cMaxTpcSharedClusters", 0.4, "Tpc Max Shared Clusters"};
  Configurable<float> cMaxChi2Tpc{"cMaxChi2Tpc", 4, "Max Chi2 Tpc"};
  Configurable<double> cTpcNsigmaCut{"cTpcNsigmaCut", 3.0, "TPC NSigma Selection Cut"};
  Configurable<bool> cRemoveAmbiguousTracks{"cRemoveAmbiguousTracks", false, "Remove Ambiguous Tracks"};

  // V0s
  Configurable<double> cMinDcaProtonToPV{"cMinDcaProtonToPV", 0.02, "Minimum Proton DCAr to PV"};
  Configurable<double> cMinDcaPionToPV{"cMinDcaPionToPV", 0.06, "Minimum Pion DCAr to PV"};
  Configurable<double> cMinV0DcaDaughters{"cMinV0DcaDaughters", 0., "Minimum DCA between V0 daughters"};
  Configurable<double> cMaxV0DcaDaughters{"cMaxV0DcaDaughters", 1., "Maximum DCA between V0 daughters"};
  Configurable<double> cMinDcaV0ToPV{"cMinDcaV0ToPV", 0.0, "Minimum DCA V0 to PV"};
  Configurable<double> cMaxDcaV0ToPV{"cMaxDcaV0ToPV", 999.0, "Maximum DCA V0 to PV"};
  Configurable<double> cMinV0TransRadius{"cMinV0TransRadius", 0.5, "Minimum V0 radius from PV"};
  Configurable<double> cMaxV0TransRadius{"cMaxV0TransRadius", 999.0, "Maximum V0 radius from PV"};
  Configurable<double> cMinV0CTau{"cMinV0CTau", 0.0, "Minimum ctau"};
  Configurable<double> cMaxV0CTau{"cMaxV0CTau", 30.0, "Maximum ctau"};
  Configurable<double> cMinV0CosPA{"cMinV0CosPA", 0.995, "Minimum V0 CosPA to PV"};
  Configurable<double> cKshortRejMassWindow{"cKshortRejMassWindow", 0.01, "Reject K0Short Candidates"};
  Configurable<bool> cKshortRejFlag{"cKshortRejFlag", true, "K0short Mass Rej Flag"};

  // V0s kinmatic acceptance
  Configurable<float> cMinV0Mass{"cMinV0Mass", 1.10, "V0 Mass Min"};
  Configurable<float> cMaxV0Mass{"cMaxV0Mass", 1.12, "V0 Mass Min"};
  Configurable<float> cMinV0Pt{"cMinV0Pt", 0.8, "Minimum V0 pT"};
  Configurable<float> cMaxV0Pt{"cMaxV0Pt", 4.2, "Minimum V0 pT"};
  Configurable<float> cMaxV0Rap{"cMaxV0Rap", 0.5, "|rap| cut"};
  Configurable<bool> cDoEtaAnalysis{"cDoEtaAnalysis", false, "Do Eta Analysis"};
  Configurable<bool> cV0TypeSelFlag{"cV0TypeSelFlag", false, "V0 Type Selection Flag"};
  Configurable<int> cV0TypeSelection{"cV0TypeSelection", 1, "V0 Type Selection"};

  // V0s MC
  // [Phase 16m] Removed dead Configurables: cHasMcFlag, cGenPrimaryLambda,
  // cGenSecondaryLambda — declared but never gated any logic. cSelMCPSV0
  // kept as a no-op for back-compat (deprecated in Phase 16l).
  Configurable<bool> cSelectTrueLambda{"cSelectTrueLambda", true, "Select True Lambda"};
  Configurable<bool> cSelMCPSV0{"cSelMCPSV0", true, "[DEPRECATED Phase 16l] no-op; v0PrmScd is always tagged from MC truth"};
  Configurable<bool> cCheckRecoDauFlag{"cCheckRecoDauFlag", true, "Check for reco daughter PID"};
  Configurable<bool> cGenDecayChannel{"cGenDecayChannel", true, "Gen Level Decay Channel Flag"};
  Configurable<bool> cRecoMomResoFlag{"cRecoMomResoFlag", false, "Check effect of momentum space smearing on balance function"};

  // [Phase 14] Diagnostic mode — emit EVERY V0 candidate (not just those
  // passing selV0Particle) into LambdaTracks, with the full cut bitmask
  // showing which stage each one failed. Use to build offline cutflow
  // plots and tune cut values without re-running. Default false because
  // it inflates the table by ~5-10× and downstream consumers need to
  // gate on cutBits / passesPrimaryTopo.
  Configurable<bool> cFillLambdaTreeAllCandidates{"cFillLambdaTreeAllCandidates", false,
      "Diagnostic: emit all V0 candidates with cutBits, not only the ones passing selV0Particle"};

  // Efficiency Correction
  Configurable<bool> cCorrectionFlag{"cCorrectionFlag", false, "Correction Flag"};
  Configurable<bool> cGetEffFact{"cGetEffFact", false, "Get Efficiency Factor Flag"};
  Configurable<bool> cGetPrimFrac{"cGetPrimFrac", false, "Get Primary Fraction Flag"};
  Configurable<int> cCorrFactHist{"cCorrFactHist", 0, "Efficiency Factor Histogram"};
  Configurable<int> cPrimFracHist{"cPrimFracHist", 0, "Primary Fraction Histogram"};

  // CCDB
  Configurable<std::string> cUrlCCDB{"cUrlCCDB", "http://ccdb-test.cern.ch:8080", "url of ccdb"};
  Configurable<std::string> cPathCCDB{"cPathCCDB", "Users/y/ypatley/lambda_corr_fact", "Path for ccdb-object"};

  // Initialize CCDB Service (shared with cascade-side path; LTP's cUrlCCDB
  // is the live URL — see init() comment).
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // [Phase 4] PDG service for cascade competing-mass cut, brought in from
  // the former CascadeSelector struct.
  Service<o2::framework::O2DatabasePDG> pdgDB;

  // ===========================================================================
  // [Phase 4 fix] All cascade-side Configurables and Axes bundled into a
  // single ConfigurableGroup. Reason: the framework's StructToTuple has a
  // hard cap of 99 task members (DPL_HOMOGENEOUS_APPLY_ENTRY (9,9)). Without
  // grouping, the merged struct lands at 129 members. ConfigurableGroup
  // makes the whole bundle count as one member while still exposing every
  // inner Configurable to JSON/CLI (no behaviour change for users — JSON keys
  // unchanged because each Configurable's first-arg name is preserved).
  // Code references go from `tpcNsigmaProton` → `cascCfg.tpcNsigmaProton`.
  // ===========================================================================
  struct : ConfigurableGroup {
    // [Phase 4] Deprecated CSEL-style event-selection knobs. Echoed in the
    // init-time deprecation warning but no longer affect event acceptance.
    Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "CCDB url (cascade-side; cosmetic — LTP's cUrlCCDB drives the actual CCDB fetches)"};
    // [Phase 16e] Removed unused Configurables: useTrigger / triggerList.
    // Trigger-skim path was deleted in an earlier phase but the knobs were
    // left behind — toggling them in JSON had zero effect. Cleared now.
    Configurable<bool> doTFBorderCut{"doTFBorderCut", true, "[DEPRECATED Phase 4] event selection delegated to LTP"};
    Configurable<bool> doSel8{"doSel8", true, "[DEPRECATED Phase 4] sel8 is enforced by LTP's cSel8Trig"};
    Configurable<bool> doNoSameBunchPileUp{"doNoSameBunchPileUp", true, "[DEPRECATED Phase 4] pileup veto is enforced by LTP's cPileupReject"};
    Configurable<int> INEL{"INEL", 0, "[DEPRECATED Phase 4] INEL>N enforcement is no longer applied"}; // o2-linter: disable=name/configurable (back-compat: deprecated ALICE INEL knob name kept for JSON config compatibility)
    Configurable<double> maxVertexZ{"maxVertexZ", 10., "[DEPRECATED Phase 4] |Vz| cut is enforced by LTP's cMin/cMaxZVtx"};

    // Cascade kinematic / selection.
    Configurable<float> etaCascades{"etaCascades", 0.8, "min/max of eta for cascades"};
    Configurable<bool> doCompetingMassCut{"doCompetingMassCut", true, "Switch to apply a competing mass cut for the Omega's"};
    Configurable<float> competingMassWindow{"competingMassWindow", 0.01, "Mass window for the competing mass cut"};
    // [Phase 12b] Cascade pT lower cap. Drops very low-pT cascades where
    // topology resolution is poor and the V0+bachelor combinatorial
    // background dominates. Default 0.6 GeV/c is a typical PWG-LF Ξ floor.
    Configurable<float> cMinCascPt{"cMinCascPt", 0.6f, "Minimum cascade pT [GeV/c]"};

    // Cascade tracklevel.
    Configurable<float> tpcNsigmaBachelor{"tpcNsigmaBachelor", 3, "TPC NSigma bachelor"};
    Configurable<float> tpcNsigmaProton{"tpcNsigmaProton", 3, "TPC NSigma proton <- lambda"};
    Configurable<float> tpcNsigmaPion{"tpcNsigmaPion", 3, "TPC NSigma pion <- lambda"};
    Configurable<int> minTPCCrossedRows{"minTPCCrossedRows", 80, "min N TPC crossed rows"};
    Configurable<int> minITSClusters{"minITSClusters", 4, "minimum number of ITS clusters"};
    Configurable<float> etaTracks{"etaTracks", 1.0, "min/max of eta for cascade daughter tracks"};
    Configurable<float> tpcChi2{"tpcChi2", 4, "TPC Chi2 (cascade tracks)"};
    Configurable<float> itsChi2{"itsChi2", 36, "ITS Chi2 (cascade tracks)"};
    // [Phase 12b] Bachelor-specific track-quality knobs. Default-equal to
    // the V0-daughter values so existing behaviour is preserved unless the
    // user tightens them. A real bachelor π/K typically wants stricter cuts
    // than a V0 daughter because mis-association is more common at the
    // cascade vertex (further from PV, larger combinatorial pool).
    Configurable<int> minBachTPCCrossedRows{"minBachTPCCrossedRows", 80, "min N TPC crossed rows for bachelor"};
    Configurable<int> minBachITSClusters{"minBachITSClusters", 4, "min ITS clusters for bachelor"};
    Configurable<float> maxBachTpcChi2{"maxBachTpcChi2", 4, "max TPC chi2 / cluster for bachelor"};
    Configurable<float> maxBachItsChi2{"maxBachItsChi2", 36, "max ITS chi2 / cluster for bachelor"};

    // Cascade selection criteria.
    // Cascade selection criteria.
    // [Phase 16u] The following Configurables intentionally use snake_case
    // names matching the long-standing ALICE V0/cascade-builder JSON
    // convention shared across PWGCF / PWGLF tasks. Renaming would break
    // back-compat with every existing config file.
    Configurable<double> v0setting_cospa{"v0setting_cospa", 0.995, "v0setting_cospa"}; // o2-linter: disable=name/configurable (back-compat: V0-builder JSON key)
    Configurable<float> v0setting_dcav0dau{"v0setting_dcav0dau", 1.0, "v0setting_dcav0dau"}; // o2-linter: disable=name/configurable (back-compat: V0/cascade-builder JSON key)
    Configurable<float> v0setting_dcapostopv{"v0setting_dcapostopv", 0.1, "v0setting_dcapostopv"}; // o2-linter: disable=name/configurable (back-compat: V0/cascade-builder JSON key)
    Configurable<float> v0setting_dcanegtopv{"v0setting_dcanegtopv", 0.1, "v0setting_dcanegtopv"}; // o2-linter: disable=name/configurable (back-compat: V0/cascade-builder JSON key)
    Configurable<float> v0setting_radius{"v0setting_radius", 0.9, "v0setting_radius"}; // o2-linter: disable=name/configurable (back-compat: V0/cascade-builder JSON key)
    Configurable<double> cascadesetting_cospa{"cascadesetting_cospa", 0.95, "cascadesetting_cospa"}; // o2-linter: disable=name/configurable (back-compat: cascade-builder JSON key)
    // [Phase 16e] Removed unused Configurables: cascadesetting_dcacascdau /
    // cascadesetting_dcabachtopv. These were declared but never used in the
    // cut chain — only `casc.dcacascdaughters()` / `casc.dcabachtopv()` (the
    // row values) are read for QA fills. The cascade-builder upstream applies
    // the analogous cuts at table-production time; our task didn't re-cut.
    Configurable<float> cascadesetting_cascradius{"cascadesetting_cascradius", 0.9, "cascadesetting_cascradius"}; // o2-linter: disable=name/configurable (back-compat: V0/cascade-builder JSON key)
    Configurable<float> cascadesetting_v0masswindow{"cascadesetting_v0masswindow", 0.01, "cascadesetting_v0masswindow"}; // o2-linter: disable=name/configurable (back-compat: V0/cascade-builder JSON key)
    Configurable<float> cascadesetting_mindcav0topv{"cascadesetting_mindcav0topv", 0.01, "cascadesetting_mindcav0topv"}; // o2-linter: disable=name/configurable (back-compat: V0/cascade-builder JSON key)

    // Cascade-side ConfigurableAxes.
    ConfigurableAxis cascRadiusAxis{"cascRadiusAxis", {100, 0.0f, 50.0f}, "cm"};
    ConfigurableAxis cascCpaAxis{"cascCpaAxis", {100, 0.95f, 1.0f}, "CPA"};
    ConfigurableAxis cascVertexAxis{"cascVertexAxis", {100, -10.0f, 10.0f}, "cm"};
    ConfigurableAxis cascDcaAxis{"cascDcaAxis", {100, 0.0f, 2.0f}, "cm"};
    ConfigurableAxis invXiMassAxis{"invXiMassAxis", {100, 1.28f, 1.38f}, "Inv. Mass (GeV/c^{2})"};
    ConfigurableAxis invOmegaMassAxis{"invOmegaMassAxis", {100, 1.62f, 1.72f}, "Inv. Mass (GeV/c^{2})"};
    ConfigurableAxis cascPtAxis{"cascPtAxis", {150, 0, 15}, "#it{p}_{T}"};
    ConfigurableAxis cascRapidityAxis{"cascRapidityAxis", {100, -1.f, 1.f}, "y"};
    ConfigurableAxis invLambdaMassAxis{"invLambdaMassAxis", {100, 1.07f, 1.17f}, "Inv. Mass (GeV/c^{2})"};
  } cascCfg;

  // [Phase 8] Topology bundle for "primary-Λ-on-data" trigger purity.
  // Wrapped in a ConfigurableGroup to (a) keep them under one logical
  // section in JSON and (b) avoid blowing the StructToTuple 99-member
  // limit of the parent task struct (the same trick we used in Phase 4
  // for cascCfg).
  // Defaults are tuned for pp/pPb at moderate Λ pT. PbPb may want a
  // tighter cMaxV0Radius (large feed-down tail) and a looser
  // cMinDauItsNCls (low-pT acceptance loss in central events).
  struct : ConfigurableGroup {
    // Master switches.
    Configurable<bool> cPrimEnable{"cPrimEnable", true, "Compute and write the lambdatrack::passesPrimaryTopo flag"};
    Configurable<bool> cPrimRequireBothDauItsHits{"cPrimRequireBothDauItsHits", true, "Require BOTH V0 daughters to have >= cPrimMinDauItsNCls ITS hits"};
    // Topological cuts.
    Configurable<float> cPrimMaxDcaV0ToPv{"cPrimMaxDcaV0ToPv", 0.10f, "[cm] V0 impact param to PV (primary-Λ knob)"};
    Configurable<float> cPrimMinV0CosPA{"cPrimMinV0CosPA", 0.999f, "Min V0 cosPA(PV); primary Λ → 1 within ~3e-4"};
    Configurable<float> cPrimMaxV0Radius{"cPrimMaxV0Radius", 30.0f, "[cm] Upper bound on V0 transverse decay radius"};
    Configurable<float> cPrimMaxDauDcaToPv{"cPrimMaxDauDcaToPv", 1.0f, "[cm] Soft upper bound on |daughter DCA-XY-to-PV|"};
    Configurable<int> cPrimMinDauItsNCls{"cPrimMinDauItsNCls", 1, "Minimum ITS hits per V0 daughter"};
    // [Phase 10] Pseudo-proper-decay-length upper bound. Set to a few cτ
    // beyond the primary peak to keep most primary Λ; feed-down's L_proper
    // shifts upward by the parent flight, so a cap here cleanly trims the tail.
    // Default 25 cm ≈ 3 cτ — tunable per centrality / collision system.
    Configurable<float> cPrimMaxLProper{"cPrimMaxLProper", 25.0f, "[cm] Upper bound on L_proper = L_xy * M_Λ / pT (kills feed-down tail)"};
    // [Phase 10] Require at least ONE V0 daughter to have a hit on the
    // ITS Inner Barrel (Layers 0,1,2 — bits 0,1,2 of itsClusterMap).
    // Strong primary tag at low pT: a feed-down Λ at low pT typically
    // has its decay vertex outside the IB, so daughters miss those layers.
    Configurable<bool> cPrimRequireItsIBHit{"cPrimRequireItsIBHit", true, "Require >=1 V0 daughter with an ITS Inner-Barrel hit (Layers 0-2)"};
    // Topology consistency for cascades (here so it lives next to the rest
    // of the purity controls; checked in processCandidate).
    Configurable<bool> cReqCascRadiusLessThanV0Radius{"cReqCascRadiusLessThanV0Radius", true, "Cascade topology consistency: cascRadius < v0Radius"};
  } primCfg;

  // Plain (non-Configurable) AxisSpec — keeping these out of the group so we
  // don't depend on ConfigurableGroup supporting non-Configurable members.
  // They were never user-tunable in the original CSEL anyway.
  AxisSpec cascItsClustersAxis{8, -0.5, 7.5, "number of ITS clusters"};
  AxisSpec cascTpcRowsAxis{160, -0.5, 159.5, "TPC crossed rows"};

  // Histogram Registry (Lambda-side — preserves all LTP-style paths).
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // [Phase 4] Cascade-side HistogramRegistry — verbatim port of CSEL's
  // initializer-list registry. Histogram paths under this registry are
  // identical to what CSEL emitted before, so AnalysisResults.root structure
  // is unchanged downstream.
  HistogramRegistry cascRegistry{
    "cascRegistry",
    {
      {"hV0Radius", "hV0Radius", {HistType::kTH3F, {cascCfg.cascRadiusAxis, cascCfg.invXiMassAxis, cascCfg.cascPtAxis}}},
      {"hCascRadius", "hCascRadius", {HistType::kTH3F, {cascCfg.cascRadiusAxis, cascCfg.invXiMassAxis, cascCfg.cascPtAxis}}},
      {"hV0CosPA", "hV0CosPA", {HistType::kTH3F, {cascCfg.cascCpaAxis, cascCfg.invXiMassAxis, cascCfg.cascPtAxis}}},
      {"hCascCosPA", "hCascCosPA", {HistType::kTH3F, {cascCfg.cascCpaAxis, cascCfg.invXiMassAxis, cascCfg.cascPtAxis}}},
      {"hDCAPosToPV", "hDCAPosToPV", {HistType::kTH3F, {cascCfg.cascVertexAxis, cascCfg.invXiMassAxis, cascCfg.cascPtAxis}}},
      {"hDCANegToPV", "hDCANegToPV", {HistType::kTH3F, {cascCfg.cascVertexAxis, cascCfg.invXiMassAxis, cascCfg.cascPtAxis}}},
      {"hDCABachToPV", "hDCABachToPV", {HistType::kTH3F, {cascCfg.cascVertexAxis, cascCfg.invXiMassAxis, cascCfg.cascPtAxis}}},
      {"hDCAV0ToPV", "hDCAV0ToPV", {HistType::kTH3F, {cascCfg.cascVertexAxis, cascCfg.invXiMassAxis, cascCfg.cascPtAxis}}},
      {"hDCAV0Dau", "hDCAV0Dau", {HistType::kTH3F, {cascCfg.cascDcaAxis, cascCfg.invXiMassAxis, cascCfg.cascPtAxis}}},
      {"hDCACascDau", "hDCACascDau", {HistType::kTH3F, {cascCfg.cascDcaAxis, cascCfg.invXiMassAxis, cascCfg.cascPtAxis}}},
      {"hLambdaMass", "hLambdaMass", {HistType::kTH3F, {cascCfg.invLambdaMassAxis, cascCfg.invXiMassAxis, cascCfg.cascPtAxis}}},

      {"hMassXiMinus", "hMassXiMinus", {HistType::kTH3F, {cascCfg.invXiMassAxis, cascCfg.cascPtAxis, cascCfg.cascRapidityAxis}}},
      {"hMassXiPlus", "hMassXiPlus", {HistType::kTH3F, {cascCfg.invXiMassAxis, cascCfg.cascPtAxis, cascCfg.cascRapidityAxis}}},
      {"hMassOmegaMinus", "hMassOmegaMinus", {HistType::kTH3F, {cascCfg.invOmegaMassAxis, cascCfg.cascPtAxis, cascCfg.cascRapidityAxis}}},
      {"hMassOmegaPlus", "hMassOmegaPlus", {HistType::kTH3F, {cascCfg.invOmegaMassAxis, cascCfg.cascPtAxis, cascCfg.cascRapidityAxis}}},

      {"hTPCnCrossedRowsPos", "hTPCnCrossedRowsPos", {HistType::kTH3F, {cascTpcRowsAxis, cascCfg.invXiMassAxis, cascCfg.cascPtAxis}}},
      {"hTPCnCrossedRowsNeg", "hTPCnCrossedRowsNeg", {HistType::kTH3F, {cascTpcRowsAxis, cascCfg.invXiMassAxis, cascCfg.cascPtAxis}}},
      {"hTPCnCrossedRowsBach", "hTPCnCrossedRowsBach", {HistType::kTH3F, {cascTpcRowsAxis, cascCfg.invXiMassAxis, cascCfg.cascPtAxis}}},
      {"hITSnClustersPos", "hITSnClustersPos", {HistType::kTH3F, {cascItsClustersAxis, cascCfg.invXiMassAxis, cascCfg.cascPtAxis}}},
      {"hITSnClustersNeg", "hITSnClustersNeg", {HistType::kTH3F, {cascItsClustersAxis, cascCfg.invXiMassAxis, cascCfg.cascPtAxis}}},
      {"hITSnClustersBach", "hITSnClustersBach", {HistType::kTH3F, {cascItsClustersAxis, cascCfg.invXiMassAxis, cascCfg.cascPtAxis}}},
      {"hTPCChi2Pos", "hTPCChi2Pos", {HistType::kTH1F, {{100, 0, 10, "TPC Chi2 Pos"}}}},
      {"hTPCChi2Neg", "hTPCChi2Neg", {HistType::kTH1F, {{100, 0, 10, "TPC Chi2 Neg"}}}},
      {"hTPCChi2Bach", "hTPCChi2Bach", {HistType::kTH1F, {{100, 0, 10, "TPC Chi2 Bach"}}}},
      {"hITSChi2Pos", "hITSChi2Pos", {HistType::kTH1F, {{100, 0, 100, "ITS Chi2 Pos"}}}},
      {"hITSChi2Neg", "hITSChi2Neg", {HistType::kTH1F, {{100, 0, 100, "ITS Chi2 Neg"}}}},
      {"hITSChi2Bach", "hITSChi2Bach", {HistType::kTH1F, {{100, 0, 100, "ITS Chi2 Bach"}}}},

      {"hTriggerQA", "hTriggerQA", {HistType::kTH1F, {{2, -0.5, 1.5, "Trigger y/n"}}}},
    },
  };

  // initialize corr_factor objects
  std::vector<std::vector<std::string>> vCorrFactStrings = {{"hEffVsPtCentLambda", "hEffVsPtCentAntiLambda"},
                                                            {"hEffVsPtYCentLambda", "hEffVsPtYCentAntiLambda"},
                                                            {"hEffVsPtEtaCentLambda", "hEffVsPtEtaCentAntiLambda"}};

  // initialize corr_factor objects
  std::vector<std::vector<std::string>> vPrimFracStrings = {{"hPrimFracVsPtCentLambda", "hPrimFracVsPtCentAntiLambda"},
                                                            {"hPrimFracVsPtYCentLambda", "hPrimFracVsPtYCentAntiLambda"},
                                                            {"hPrimFracVsPtEtaCentLambda", "hPrimFracVsPtEtaCentAntiLambda"}};

  // Initialize Global Variables
  float cent = 0., mult = 0.;
  float pt = 0., eta = 0., rap = 0., phi = 0.;

  void init(InitContext const&)
  {
    // Set CCDB url
    ccdb->setURL(cUrlCCDB.value);
    ccdb->setCaching(true);

    // [DEBUG] Echo the active configuration once at startup so the run log
    // shows exactly which cuts were applied. Search the log for "[CFG-LTP]".
    LOGF(info,
         "[CFG-LTP] CentEst=%d Vz=[%.2f,%.2f] Mult=[%.2f,%.2f] | "
         "trig: sel8=%d int7=%d sel7=%d tvx=%d tfBorder=%d itsROBorder=%d "
         "itsTpcVtx=%d pileup=%d zVtxTimeDiff=%d goodITS=%d | "
         "track: pT=[%.3f,%.3f] |eta|<%.2f tpcCR>=%d tpcCR/Findable>=%.2f "
         "tpcShared<=%.2f tpcChi2<=%.2f tpcNsigma=%.2f rmAmbig=%d | "
         "V0: dcaPr>=%.3f dcaPi>=%.3f dcaDau=[%.3f,%.3f] dcaV0=[%.3f,%.3f] "
         "r=[%.2f,%.2f] ctau=[%.2f,%.2f] cosPA>=%.4f K0Rej=%.3f(flag=%d) | "
         "V0kin: m=[%.3f,%.3f] pT=[%.2f,%.2f] |y|<%.2f doEta=%d | "
         "MC: trueLambda=%d primSec=%d checkDau=%d "
         "genDecCh=%d momReso=%d | corr=%d effFlag=%d primFracFlag=%d",
         (int)cCentEstimator, (float)cMinZVtx, (float)cMaxZVtx,
         (float)cMinMult, (float)cMaxMult,
         (int)cSel8Trig, (int)cInt7Trig, (int)cSel7Trig,
         (int)cTriggerTvxSel, (int)cTFBorder, (int)cNoItsROBorder,
         (int)cItsTpcVtx, (int)cPileupReject, (int)cZVtxTimeDiff,
         (int)cIsGoodITSLayers,
         (float)cTrackMinPt, (float)cTrackMaxPt, (float)cTrackEtaCut,
         (int)cMinTpcCrossedRows, (float)cMinTpcCROverCls,
         (float)cMaxTpcSharedClusters, (float)cMaxChi2Tpc,
         (double)cTpcNsigmaCut, (int)cRemoveAmbiguousTracks,
         (double)cMinDcaProtonToPV, (double)cMinDcaPionToPV,
         (double)cMinV0DcaDaughters, (double)cMaxV0DcaDaughters,
         (double)cMinDcaV0ToPV, (double)cMaxDcaV0ToPV,
         (double)cMinV0TransRadius, (double)cMaxV0TransRadius,
         (double)cMinV0CTau, (double)cMaxV0CTau,
         (double)cMinV0CosPA, (double)cKshortRejMassWindow, (int)cKshortRejFlag,
         (float)cMinV0Mass, (float)cMaxV0Mass,
         (float)cMinV0Pt, (float)cMaxV0Pt, (float)cMaxV0Rap, (int)cDoEtaAnalysis,
         (int)cSelectTrueLambda, (int)cSelMCPSV0,
         (int)cCheckRecoDauFlag,
         (int)cGenDecayChannel, (int)cRecoMomResoFlag,
         (int)cCorrectionFlag, (int)cGetEffFact, (int)cGetPrimFrac);

    // [P1][R1] Unified event-selection log (compare to [EVENTSEL-CSEL]).
    lcorr_evsel::logEventCuts("EVENTSEL-LTP", buildEventCuts());

    // [Phase 8+10+14] Primary-Λ topology bundle + diagnostic-mode flag.
    LOGF(info,
         "[CFG-LCP-prim] enable=%d reqBothDauITS=%d reqItsIB=%d "
         "maxDcaV0ToPv=%.4f minV0CosPA=%.6f maxV0Radius=%.2f maxLProper=%.2f "
         "maxDauDcaToPv=%.4f minDauItsNCls=%d "
         "cascRadius<v0Radius=%d cFillLambdaTreeAllCandidates=%d",
         (int)primCfg.cPrimEnable, (int)primCfg.cPrimRequireBothDauItsHits,
         (int)primCfg.cPrimRequireItsIBHit,
         (float)primCfg.cPrimMaxDcaV0ToPv, (float)primCfg.cPrimMinV0CosPA,
         (float)primCfg.cPrimMaxV0Radius, (float)primCfg.cPrimMaxLProper,
         (float)primCfg.cPrimMaxDauDcaToPv, (int)primCfg.cPrimMinDauItsNCls,
         (int)primCfg.cReqCascRadiusLessThanV0Radius,
         (int)cFillLambdaTreeAllCandidates);

    // initialize axis specifications
    const AxisSpec axisCols(5, 0.5, 5.5, "");
    const AxisSpec axisTrks(30, 0.5, 30.5, "");
    const AxisSpec axisCent(100, 0, 100, "FT0M (%)");
    const AxisSpec axisMult(10, 0, 10, "N_{#Lambda}");
    const AxisSpec axisVz(220, -11, 11, "V_{z} (cm)");
    const AxisSpec axisPID(8000, -4000, 4000, "PdgCode");

    const AxisSpec axisV0Mass(140, 1.08, 1.15, "M_{p#pi} (GeV/#it{c}^{2})");
    const AxisSpec axisV0Pt(100., 0., 10., "p_{T} (GeV/#it{c})");
    const AxisSpec axisV0Rap(48, -1.2, 1.2, "y");
    const AxisSpec axisV0Eta(48, -1.2, 1.2, "#eta");
    const AxisSpec axisV0Phi(36, 0., TwoPI, "#phi (rad)");

    const AxisSpec axisRadius(2000, 0, 200, "r(cm)");
    const AxisSpec axisCosPA(300, 0.97, 1.0, "cos(#theta_{PA})");
    const AxisSpec axisDcaV0PV(1000, 0., 10., "dca (cm)");
    const AxisSpec axisDcaProngPV(5000, -50., 50., "dca (cm)");
    const AxisSpec axisDcaDau(75, 0., 1.5, "Daug DCA (#sigma)");
    const AxisSpec axisCTau(2000, 0, 200, "c#tau (cm)");
    const AxisSpec axisGCTau(2000, 0, 200, "#gammac#tau (cm)");
    const AxisSpec axisAlpha(40, -1, 1, "#alpha");
    const AxisSpec axisQtarm(40, 0, 0.4, "q_{T}");

    const AxisSpec axisTrackPt(40, 0, 4, "p_{T} (GeV/#it{c})");
    const AxisSpec axisTrackDCA(200, -1, 1, "dca_{XY} (cm)");
    const AxisSpec axisMomPID(80, 0, 4, "p (GeV/#it{c})");
    const AxisSpec axisNsigma(401, -10.025, 10.025, {"n#sigma"});
    const AxisSpec axisdEdx(360, 20, 200, "#frac{dE}{dx}");

    // Create Histograms.
    // Event histograms
    histos.add("Events/h1f_collisions_info", "# of Collisions", kTH1F, {axisCols});
    histos.add("Events/h1f_collision_posZ", "V_{z}-distribution", kTH1F, {axisVz});

    // QA
    histos.add("Tracks/h1f_tracks_info", "# of tracks", kTH1F, {axisTrks});
    histos.add("Tracks/h2f_armpod_before_sel", "Armentros-Podolanski Plot", kTH2F, {axisAlpha, axisQtarm});
    histos.add("Tracks/h2f_armpod_after_sel", "Armentros-Podolanski Plot", kTH2F, {axisAlpha, axisQtarm});
    histos.add("Tracks/h1f_lambda_pt_vs_invm", "p_{T} vs M_{#Lambda}", kTH2F, {axisV0Mass, axisV0Pt});
    histos.add("Tracks/h1f_antilambda_pt_vs_invm", "p_{T} vs M_{#bar{#Lambda}}", kTH2F, {axisV0Mass, axisV0Pt});

    // QA Lambda
    histos.add("QA/Lambda/h2f_qt_vs_alpha", "Armentros-Podolanski Plot", kTH2F, {axisAlpha, axisQtarm});
    histos.add("QA/Lambda/h1f_dca_V0_daughters", "DCA between V0 daughters", kTH1F, {axisDcaDau});
    histos.add("QA/Lambda/h1f_dca_pos_to_PV", "DCA positive prong to PV", kTH1F, {axisDcaProngPV});
    histos.add("QA/Lambda/h1f_dca_neg_to_PV", "DCA negative prong to PV", kTH1F, {axisDcaProngPV});
    histos.add("QA/Lambda/h1f_dca_V0_to_PV", "DCA V0 to PV", kTH1F, {axisDcaV0PV});
    histos.add("QA/Lambda/h1f_V0_cospa", "cos(#theta_{PA})", kTH1F, {axisCosPA});
    histos.add("QA/Lambda/h1f_V0_radius", "V_{0} Decay Radius in XY plane", kTH1F, {axisRadius});
    histos.add("QA/Lambda/h1f_V0_ctau", "V_{0} c#tau", kTH1F, {axisCTau});
    histos.add("QA/Lambda/h1f_V0_gctau", "V_{0} #gammac#tau", kTH1F, {axisGCTau});

    histos.add("QA/Lambda/h1f_pos_prong_pt", "Pos-Prong p_{T}", kTH1F, {axisTrackPt});
    histos.add("QA/Lambda/h1f_neg_prong_pt", "Neg-Prong p_{T}", kTH1F, {axisTrackPt});
    histos.add("QA/Lambda/h1f_pos_prong_eta", "Pos-Prong #eta-distribution", kTH1F, {axisV0Eta});
    histos.add("QA/Lambda/h1f_neg_prong_eta", "Neg-Prong #eta-distribution", kTH1F, {axisV0Eta});
    histos.add("QA/Lambda/h1f_pos_prong_phi", "Pos-Prong #phi-distribution", kTH1F, {axisV0Phi});
    histos.add("QA/Lambda/h1f_neg_prong_phi", "Neg-Prong #phi-distribution", kTH1F, {axisV0Phi});

    histos.add("QA/Lambda/h2f_pos_prong_dcaXY_vs_pt", "DCA vs p_{T}", kTH2F, {axisTrackPt, axisTrackDCA});
    histos.add("QA/Lambda/h2f_neg_prong_dcaXY_vs_pt", "DCA vs p_{T}", kTH2F, {axisTrackPt, axisTrackDCA});
    histos.add("QA/Lambda/h2f_pos_prong_dEdx_vs_p", "TPC Signal Pos-Prong", kTH2F, {axisMomPID, axisdEdx});
    histos.add("QA/Lambda/h2f_neg_prong_dEdx_vs_p", "TPC Signal Neg-Prong", kTH2F, {axisMomPID, axisdEdx});
    histos.add("QA/Lambda/h2f_pos_prong_tpc_nsigma_pr_vs_p", "TPC n#sigma Pos Prong", kTH2F, {axisMomPID, axisNsigma});
    histos.add("QA/Lambda/h2f_neg_prong_tpc_nsigma_pr_vs_p", "TPC n#sigma Neg Prong", kTH2F, {axisMomPID, axisNsigma});
    histos.add("QA/Lambda/h2f_pos_prong_tpc_nsigma_pi_vs_p", "TPC n#sigma Pos Prong", kTH2F, {axisMomPID, axisNsigma});
    histos.add("QA/Lambda/h2f_neg_prong_tpc_nsigma_pi_vs_p", "TPC n#sigma Neg Prong", kTH2F, {axisMomPID, axisNsigma});

    // Kinematic Histograms
    histos.add("McRec/Lambda/hPt", "Transverse Momentum", kTH1F, {axisV0Pt});
    histos.add("McRec/Lambda/hEta", "Pseudorapidity", kTH1F, {axisV0Eta});
    histos.add("McRec/Lambda/hRap", "Rapidity", kTH1F, {axisV0Rap});
    histos.add("McRec/Lambda/hPhi", "Azimuthal Angle", kTH1F, {axisV0Phi});

    // QA Anti-Lambda
    histos.addClone("QA/Lambda/", "QA/AntiLambda/");
    histos.addClone("McRec/Lambda/", "McRec/AntiLambda/");

    // MC Generated Histograms
    if (doprocessMCRun3 || doprocessMCRun2 || doprocessMCRecoRun3 || doprocessMCRecoRun2) {
      // McReco Histos
      histos.add("Tracks/h2f_tracks_pid_before_sel", "PIDs", kTH2F, {axisPID, axisV0Pt});
      histos.add("Tracks/h2f_tracks_pid_after_sel", "PIDs", kTH2F, {axisPID, axisV0Pt});
      histos.add("Tracks/h2f_lambda_mothers_pdg", "PIDs", kTH2F, {axisPID, axisV0Pt});

      // McGen Histos
      histos.add("McGen/h1f_collision_recgen", "# of Reco Collision Associated to One Mc Generator Collision", kTH1F, {axisMult});
      histos.add("McGen/h1f_collisions_info", "# of collisions", kTH1F, {axisCols});
      histos.add("McGen/h2f_collision_posZ", "V_{z}-distribution", kTH2F, {axisVz, axisVz});
      histos.add("McGen/h2f_collision_cent", "FT0M Centrality", kTH2F, {axisCent, axisCent});
      histos.add("McGen/h1f_lambda_daughter_PDG", "PDG Daughters", kTH1F, {axisPID});
      histos.add("McGen/h1f_antilambda_daughter_PDG", "PDG Daughters", kTH1F, {axisPID});

      histos.addClone("McRec/", "McGen/");

      histos.add("McGen/Lambda/Proton/hPt", "Proton p_{T}", kTH1F, {axisTrackPt});
      histos.add("McGen/Lambda/Proton/hEta", "Proton #eta", kTH1F, {axisV0Eta});
      histos.add("McGen/Lambda/Proton/hRap", "Proton y", kTH1F, {axisV0Rap});
      histos.add("McGen/Lambda/Proton/hPhi", "Proton #phi", kTH1F, {axisV0Phi});

      histos.addClone("McGen/Lambda/Proton/", "McGen/Lambda/Pion/");
      histos.addClone("McGen/Lambda/Proton/", "McGen/AntiLambda/Proton/");
      histos.addClone("McGen/Lambda/Pion/", "McGen/AntiLambda/Pion/");

      // set bin lables specific to MC
      histos.get<TH1>(HIST("Events/h1f_collisions_info"))->GetXaxis()->SetBinLabel(CollisionLabels::kTotColBeforeHasMcCollision, "kTotColBeforeHasMcCollision");
      histos.get<TH1>(HIST("McGen/h1f_collisions_info"))->GetXaxis()->SetBinLabel(CollisionLabels::kTotCol, "kTotCol");
      histos.get<TH1>(HIST("McGen/h1f_collisions_info"))->GetXaxis()->SetBinLabel(CollisionLabels::kPassSelCol, "kPassSelCol");
      histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kTracksBeforeHasMcParticle, "kTracksBeforeHasMcParticle");
      histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kPrimaryLambda, "kPrimaryLambda");
      histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kSecondaryLambda, "kSecondaryLambda");
      histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kLambdaDauNotMcParticle, "kLambdaDauNotMcParticle");
      histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kLambdaNotPrPiMinus, "kLambdaNotPrPiMinus");
      histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kAntiLambdaNotAntiPrPiPlus, "kAntiLambdaNotAntiPrPiPlus");
      histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kPassTrueLambdaSel, "kPassTrueLambdaSel");
      histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kGenTotAccLambda, "kGenTotAccLambda");
      histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kGenLambdaNoDau, "kGenLambdaNoDau");
      histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kGenLambdaToPrPi, "kGenLambdaToPrPi");
    }

    // set bin labels
    histos.get<TH1>(HIST("Events/h1f_collisions_info"))->GetXaxis()->SetBinLabel(CollisionLabels::kTotCol, "kTotCol");
    histos.get<TH1>(HIST("Events/h1f_collisions_info"))->GetXaxis()->SetBinLabel(CollisionLabels::kPassSelCol, "kPassSelCol");
    histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kAllV0Tracks, "kAllV0Tracks");
    histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kV0KShortMassRej, "kV0KShortMassRej");
    histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kNotLambdaNotAntiLambda, "kNotLambdaNotAntiLambda");
    histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kV0IsBothLambdaAntiLambda, "kV0IsBothLambdaAntiLambda");
    histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kNotLambdaAfterSel, "kNotLambdaAfterSel");
    histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kV0IsLambdaOrAntiLambda, "kV0IsLambdaOrAntiLambda");
    histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kPassV0DauTrackSel, "kPassV0DauTrackSel");
    histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kPassV0KinCuts, "kPassV0KinCuts");
    histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kPassV0TopoSel, "kPassV0TopoSel");
    histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kAllSelPassed, "kAllSelPassed");
    histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kEffCorrPtCent, "kEffCorrPtCent");
    histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kEffCorrPtRapCent, "kEffCorrPtRapCent");
    histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kNoEffCorr, "kNoEffCorr");
    histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kPFCorrPtCent, "kPFCorrPtCent");
    histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kPFCorrPtRapCent, "kPFCorrPtRapCent");
    histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kNoPFCorr, "kNoPFCorr");

    // ===== [Phase 4] Cascade-side init (verbatim from former CascadeSelector) =====
    // Cascade-side selection-status & event-selection histograms.
    auto h = cascRegistry.add<TH1>("hSelectionStatus", "hSelectionStatus", HistType::kTH1I, {{10, 0, 10, "status"}});
    h->GetXaxis()->SetBinLabel(1, "All");
    h->GetXaxis()->SetBinLabel(2, "nTPC OK");
    h->GetXaxis()->SetBinLabel(3, "nITS OK");
    h->GetXaxis()->SetBinLabel(4, "track Chi2 OK");
    h->GetXaxis()->SetBinLabel(5, "Topo OK");
    h->GetXaxis()->SetBinLabel(6, "Track eta OK");
    h->GetXaxis()->SetBinLabel(7, "Cascade eta OK");
    h->GetXaxis()->SetBinLabel(8, "V0 PID OK");
    h->GetXaxis()->SetBinLabel(9, "Bach PID OK");

    // [Phase 4] Cascade-side hEventSel: bins reinterpreted as
    //   0 = "All cascades seen"
    //   1 = "Event accepted by LTP -> casc proceeds to candidate cuts"
    //   2 = "Event rejected by LTP -> casc auto-flagged 0"
    auto hEventSel = cascRegistry.add<TH1>("hEventSel", "hEventSel", HistType::kTH1I, {{3, 0, 3, "0=all, 1=LTP-accepted, 2=LTP-rejected"}});
    hEventSel->GetXaxis()->SetBinLabel(1, "All");
    hEventSel->GetXaxis()->SetBinLabel(2, "LTP-accepted");
    hEventSel->GetXaxis()->SetBinLabel(3, "LTP-rejected");

    // [Phase 4] Cascade MC reco-matched histograms (created on demand).
    if (doprocessMCRecoRun3) {
      cascRegistry.add("truerec/hV0Radius", "hV0Radius", HistType::kTH1F, {cascCfg.cascRadiusAxis});
      cascRegistry.add("truerec/hCascRadius", "hCascRadius", HistType::kTH1F, {cascCfg.cascRadiusAxis});
      cascRegistry.add("truerec/hV0CosPA", "hV0CosPA", HistType::kTH1F, {cascCfg.cascCpaAxis});
      cascRegistry.add("truerec/hCascCosPA", "hCascCosPA", HistType::kTH1F, {cascCfg.cascCpaAxis});
      cascRegistry.add("truerec/hDCAPosToPV", "hDCAPosToPV", HistType::kTH1F, {cascCfg.cascVertexAxis});
      cascRegistry.add("truerec/hDCANegToPV", "hDCANegToPV", HistType::kTH1F, {cascCfg.cascVertexAxis});
      cascRegistry.add("truerec/hDCABachToPV", "hDCABachToPV", HistType::kTH1F, {cascCfg.cascVertexAxis});
      cascRegistry.add("truerec/hDCAV0ToPV", "hDCAV0ToPV", HistType::kTH1F, {cascCfg.cascVertexAxis});
      cascRegistry.add("truerec/hDCAV0Dau", "hDCAV0Dau", HistType::kTH1F, {cascCfg.cascDcaAxis});
      cascRegistry.add("truerec/hDCACascDau", "hDCACascDau", HistType::kTH1F, {cascCfg.cascDcaAxis});
      cascRegistry.add("truerec/hLambdaMass", "hLambdaMass", HistType::kTH1F, {cascCfg.invLambdaMassAxis});
      cascRegistry.add("truerec/hTPCnCrossedRowsPos", "hTPCnCrossedRowsPos", HistType::kTH1F, {cascTpcRowsAxis});
      cascRegistry.add("truerec/hTPCnCrossedRowsNeg", "hTPCnCrossedRowsNeg", HistType::kTH1F, {cascTpcRowsAxis});
      cascRegistry.add("truerec/hTPCnCrossedRowsBach", "hTPCnCrossedRowsBach", HistType::kTH1F, {cascTpcRowsAxis});
      cascRegistry.add("truerec/hITSnClustersPos", "hITSnClustersPos", HistType::kTH1F, {cascItsClustersAxis});
      cascRegistry.add("truerec/hITSnClustersNeg", "hITSnClustersNeg", HistType::kTH1F, {cascItsClustersAxis});
      cascRegistry.add("truerec/hITSnClustersBach", "hITSnClustersBach", HistType::kTH1F, {cascItsClustersAxis});
      cascRegistry.add("truerec/hTPCChi2Pos", "hTPCChi2Pos", HistType::kTH1F, {{100, 0, 10, "TPC Chi2 Pos"}});
      cascRegistry.add("truerec/hTPCChi2Neg", "hTPCChi2Neg", HistType::kTH1F, {{100, 0, 10, "TPC Chi2 Neg"}});
      cascRegistry.add("truerec/hTPCChi2Bach", "hTPCChi2Bach", HistType::kTH1F, {{100, 0, 10, "TPC Chi2 Bach"}});
      cascRegistry.add("truerec/hITSChi2Pos", "hITSChi2Pos", HistType::kTH1F, {{100, 0, 100, "ITS Chi2 Pos"}});
      cascRegistry.add("truerec/hITSChi2Neg", "hITSChi2Neg", HistType::kTH1F, {{100, 0, 100, "ITS Chi2 Neg"}});
      cascRegistry.add("truerec/hITSChi2Bach", "hITSChi2Bach", HistType::kTH1F, {{100, 0, 100, "ITS Chi2 Bach"}});
      cascRegistry.add("truerec/hXiMinus", "hXiMinus", HistType::kTH2F, {cascCfg.cascPtAxis, cascCfg.cascRapidityAxis});
      cascRegistry.add("truerec/hXiPlus", "hXiPlus", HistType::kTH2F, {cascCfg.cascPtAxis, cascCfg.cascRapidityAxis});
      cascRegistry.add("truerec/hOmegaMinus", "hOmegaMinus", HistType::kTH2F, {cascCfg.cascPtAxis, cascCfg.cascRapidityAxis});
      cascRegistry.add("truerec/hOmegaPlus", "hOmegaPlus", HistType::kTH2F, {cascCfg.cascPtAxis, cascCfg.cascRapidityAxis});
    }

    // [Phase 4] Cascade MC gen-only histograms (created on demand).
    if (doprocessCascadeGenMC) {
      cascRegistry.add("gen/hXiMinus", "hXiMinus", HistType::kTH2F, {cascCfg.cascPtAxis, cascCfg.cascRapidityAxis});
      cascRegistry.add("gen/hXiPlus", "hXiPlus", HistType::kTH2F, {cascCfg.cascPtAxis, cascCfg.cascRapidityAxis});
      cascRegistry.add("gen/hOmegaMinus", "hOmegaMinus", HistType::kTH2F, {cascCfg.cascPtAxis, cascCfg.cascRapidityAxis});
      cascRegistry.add("gen/hOmegaPlus", "hOmegaPlus", HistType::kTH2F, {cascCfg.cascPtAxis, cascCfg.cascRapidityAxis});
      cascRegistry.add("genwithrec/hXiMinus", "hXiMinus", HistType::kTH2F, {cascCfg.cascPtAxis, cascCfg.cascRapidityAxis});
      cascRegistry.add("genwithrec/hXiPlus", "hXiPlus", HistType::kTH2F, {cascCfg.cascPtAxis, cascCfg.cascRapidityAxis});
      cascRegistry.add("genwithrec/hOmegaMinus", "hOmegaMinus", HistType::kTH2F, {cascCfg.cascPtAxis, cascCfg.cascRapidityAxis});
      cascRegistry.add("genwithrec/hOmegaPlus", "hOmegaPlus", HistType::kTH2F, {cascCfg.cascPtAxis, cascCfg.cascRapidityAxis});
      cascRegistry.add("genwithrec/hNevents", "hNevents", HistType::kTH1F, {{1, 0, 1, "N generated events with reconstructed event"}});
      cascRegistry.add("gen/hNevents", "hNevents", HistType::kTH1F, {{1, 0, 1, "N generated events"}});
    }

    // [Phase 4] Loud deprecation note: CSEL-style event-selection knobs no
    // longer affect what events are accepted. The single source of truth is
    // LTP's selCollision (configured via cSel8Trig, cMin/cMaxZVtx,
    // cPileupReject, cTriggerTvxSel, cTFBorder, cNoItsROBorder, etc.).
    LOGF(info,
         "[Phase 4] CSEL-style event-selection configurables are now ignored: "
         "doSel8=%d doNoSameBunchPileUp=%d INEL=%d maxVertexZ=%.2f doTFBorderCut=%d. "
         "Event acceptance is driven exclusively by LTP's selCollision (see "
         "[CFG-LTP] above). To change event cuts, set the corresponding "
         "cSel8Trig / cPileupReject / cTriggerTvxSel / etc. on this same task.",
         (int)cascCfg.doSel8, (int)cascCfg.doNoSameBunchPileUp, (int)cascCfg.INEL,
         (double)cascCfg.maxVertexZ, (int)cascCfg.doTFBorderCut);

    // [Phase 4] CCDB URL conflict check.
    if (std::string(cUrlCCDB.value) != std::string(cascCfg.ccdbUrl.value)) {
      LOGF(warning,
           "[Phase 4] cUrlCCDB and cascCfg.ccdbUrl differ ('%s' vs '%s'). Only "
           "cUrlCCDB drives the live CCDB service; cascCfg.ccdbUrl is cosmetic.",
           cUrlCCDB.value.c_str(), cascCfg.ccdbUrl.value.c_str());
    }
  }

  // [P1][R1] Build the EventCuts struct from this producer's configurables.
  // Defaults map 1:1 to current behaviour: only the Configurables already
  // toggled on are forwarded as `use*=true` to the shared selector.
  // [Phase 4 fix] Cannot be `const`: Configurable<T>'s implicit conversion to T
  // is non-const, so reading any Configurable<> in a const-qualified method
  // fails. Drop the const qualifier; the function does not modify state anyway.
  lcorr_evsel::EventCuts buildEventCuts()
  {
    lcorr_evsel::EventCuts c;
    c.useVtxZ = true; // LTP always applies VtxZ
    c.minVtxZ = cMinZVtx;
    c.maxVtxZ = cMaxZVtx;
    c.useSel8 = cSel8Trig;
    c.useInt7 = cInt7Trig;
    c.useSel7 = cSel7Trig;
    // Centrality range is meaningful for any non-trivial window;
    // matches old behaviour exactly (was applied unconditionally).
    c.useCentRange = true;
    c.minCent = cMinMult;
    c.maxCent = cMaxMult;
    // INEL is not part of LTP today; leave off.
    c.useInel = false;
    c.useTriggerTvx = cTriggerTvxSel;
    c.useTfBorder = cTFBorder;
    c.useItsRoBorder = cNoItsROBorder;
    c.useItsTpcVtx = cItsTpcVtx;
    c.useNoSameBunchPileup = cPileupReject;
    c.useZVtxTimeDiff = cZVtxTimeDiff;
    c.useIsGoodITSLayers = cIsGoodITSLayers;
    return c;
  }

  template <RunType run, typename C>
  bool selCollision(C const& col)
  {
    // [P1][R1] LTP-specific side-effects: pick centrality estimator, then
    // delegate to the shared selector. After acceptance, set `mult` global.
    if constexpr (run == kRun3) {
      if (cCentEstimator == kCentFT0M) {
        cent = col.centFT0M();
      } else if (cCentEstimator == kCentFV0A) {
        cent = col.centFV0A();
      }
    } else {
      cent = col.centRun2V0M();
    }

    auto cuts = buildEventCuts();
    const char* reason = "ok";
    if (!lcorr_evsel::applyEventSelection<run>(col, cuts, cent, reason)) {
      LOGF(debug, "[LTP] reject: %s (Vz=%.3f cent=%.2f)",
           reason, (double)col.posZ(), (double)cent);
      return false;
    }

    // Multiplicity (LTP-specific) — distinct from CSEL's multNTracksPVeta1.
    mult = col.multNTracksPV();
    return true;
  }

  // Kinematic Selection
  bool kinCutSelection(float const& pt, float const& rap, float const& ptMin, float const& ptMax, float const& rapMax)
  {
    if (pt <= ptMin || pt >= ptMax || rap >= rapMax) {
      return false;
    }

    return true;
  }

  // Track Selection
  template <typename T>
  bool selTrack(T const& track)
  {
    if (!kinCutSelection(track.pt(), std::abs(track.eta()), cTrackMinPt, cTrackMaxPt, cTrackEtaCut)) {
      return false;
    }

    if (track.tpcNClsCrossedRows() <= cMinTpcCrossedRows) {
      return false;
    }

    if (track.tpcCrossedRowsOverFindableCls() < cMinTpcCROverCls) {
      return false;
    }

    if (track.tpcNClsShared() > cMaxTpcSharedClusters) {
      return false;
    }

    if (track.tpcChi2NCl() > cMaxChi2Tpc) {
      return false;
    }

    return true;
  }

  // Daughter Track Selection
  template <typename V, typename T>
  bool selDaughterTracks(V const& v0, T const&, ParticleType const& v0Type)
  {
    auto posTrack = v0.template posTrack_as<T>();
    auto negTrack = v0.template negTrack_as<T>();

    if (!selTrack(posTrack) || !selTrack(negTrack)) {
      return false;
    }

    // Apply DCA Selection on Daughter Tracks Based on Lambda/AntiLambda daughters
    float dcaProton = 0., dcaPion = 0.;
    if (v0Type == kLambda) {
      dcaProton = std::abs(v0.dcapostopv());
      dcaPion = std::abs(v0.dcanegtopv());
    } else if (v0Type == kAntiLambda) {
      dcaPion = std::abs(v0.dcapostopv());
      dcaProton = std::abs(v0.dcanegtopv());
    }

    if (dcaProton < cMinDcaProtonToPV || dcaPion < cMinDcaPionToPV) {
      return false;
    }

    return true;
  }

  template <typename C, typename V, typename T>
  bool topoCutSelection(C const& col, V const& v0, T const&)
  {
    // DCA
    if (v0.dcaV0daughters() <= cMinV0DcaDaughters || v0.dcaV0daughters() >= cMaxV0DcaDaughters) {
      return false;
    }

    if (v0.dcav0topv() <= cMinDcaV0ToPV || v0.dcav0topv() >= cMaxDcaV0ToPV) {
      return false;
    }

    if (v0.v0radius() <= cMinV0TransRadius || v0.v0radius() >= cMaxV0TransRadius) {
      return false;
    }

    // ctau
    float ctau = v0.distovertotmom(col.posX(), col.posY(), col.posZ()) * MassLambda0;
    if (ctau <= cMinV0CTau || ctau >= cMaxV0CTau) {
      return false;
    }

    // cosine of pointing angle
    if (v0.v0cosPA() <= cMinV0CosPA) {
      return false;
    }

    // all selection criterion passed (Return True)
    return true;
  }

  template <ParticleType part, typename T>
  bool selLambdaDauWithTpcPid(T const& postrack, T const& negtrack)
  {
    bool returnFlag = false;
    float tpcNSigmaPr = 0., tpcNSigmaPi = 0.;

    switch (part) {
      // postrack = Proton, negtrack = Pion
      case kLambda:
        tpcNSigmaPr = postrack.tpcNSigmaPr();
        tpcNSigmaPi = negtrack.tpcNSigmaPi();
        break;

      // negtrack = Proton, postrack = Pion
      case kAntiLambda:
        tpcNSigmaPr = negtrack.tpcNSigmaPr();
        tpcNSigmaPi = postrack.tpcNSigmaPi();
        break;
    }

    if (std::abs(tpcNSigmaPr) < cTpcNsigmaCut && std::abs(tpcNSigmaPi) < cTpcNsigmaCut) {
      returnFlag = true;
    }

    return returnFlag;
  }

  template <typename V, typename T>
  bool selLambdaMassWindow(V const& v0, T const&, ParticleType& v0type)
  {
    // Kshort mass rejection hypothesis
    if (cKshortRejFlag && (std::abs(v0.mK0Short() - MassK0Short) <= cKshortRejMassWindow)) {
      histos.fill(HIST("Tracks/h1f_tracks_info"), kV0KShortMassRej);
      return false;
    }

    // initialize daughter tracks
    auto postrack = v0.template posTrack_as<T>();
    auto negtrack = v0.template negTrack_as<T>();

    // initialize selection flags
    bool lambdaFlag = false, antiLambdaFlag = false;

    // get v0 track as lambda
    if ((v0.mLambda() > cMinV0Mass && v0.mLambda() < cMaxV0Mass) && (selLambdaDauWithTpcPid<kLambda>(postrack, negtrack))) {
      lambdaFlag = true;
      v0type = kLambda;
    }

    // get v0 track as anti-lambda
    if ((v0.mAntiLambda() > cMinV0Mass && v0.mAntiLambda() < cMaxV0Mass) && (selLambdaDauWithTpcPid<kAntiLambda>(postrack, negtrack))) {
      antiLambdaFlag = true;
      v0type = kAntiLambda;
    }

    if (!lambdaFlag && !antiLambdaFlag) { // neither Lambda nor Anti-Lambda
      histos.fill(HIST("Tracks/h1f_tracks_info"), kNotLambdaNotAntiLambda);
      return false;
    } else if (lambdaFlag && antiLambdaFlag) { // check if the track is identified as lambda and anti-lambda both (DISCARD THIS TRACK)
      histos.fill(HIST("Tracks/h1f_tracks_info"), kV0IsBothLambdaAntiLambda);
      return false;
    }

    if (lambdaFlag || antiLambdaFlag) {
      return true;
    }

    histos.fill(HIST("Tracks/h1f_tracks_info"), kNotLambdaAfterSel);

    return false;
  }

  template <typename C, typename V, typename T>
  bool selV0Particle(C const& col, V const& v0, T const& tracks, ParticleType& v0Type)
  {
    // Apply Lambda Mass Hypothesis
    if (!selLambdaMassWindow(v0, tracks, v0Type)) {
      return false;
    }

    histos.fill(HIST("Tracks/h1f_tracks_info"), kV0IsLambdaOrAntiLambda);

    // Apply Daughter Track Selection
    if (!selDaughterTracks(v0, tracks, v0Type)) {
      return false;
    }

    histos.fill(HIST("Tracks/h1f_tracks_info"), kPassV0DauTrackSel);

    // Apply Kinematic Selection
    float rap = 0.;
    if (!cDoEtaAnalysis) {
      rap = std::abs(v0.yLambda());
    } else {
      rap = std::abs(v0.eta());
    }

    if (!kinCutSelection(v0.pt(), rap, cMinV0Pt, cMaxV0Pt, cMaxV0Rap)) {
      return false;
    }

    histos.fill(HIST("Tracks/h1f_tracks_info"), kPassV0KinCuts);

    // Apply Topological Selection
    if (!topoCutSelection(col, v0, tracks)) {
      return false;
    }

    histos.fill(HIST("Tracks/h1f_tracks_info"), kPassV0TopoSel);

    // All Selection Criterion Passed
    return true;
  }

  // [Phase 14] Compute the per-V0 cut bitmask UNCONDITIONALLY. Each bit
  // mirrors a stage of the existing selV0Particle / topoCutSelection /
  // selDaughterTracks chain, but no early-return short-circuits. The
  // resulting uint32 captures the cut state of EVERY candidate so an
  // offline macro can replay any cut combination from the tree alone.
  //
  // Note: this duplicates the cut conditions written elsewhere — keep
  // in sync if you ever change a numerical bound. Cheap (~10 comparisons
  // per V0), only called once per V0.
  template <typename C, typename V, typename T>
  uint32_t computeLambdaCutBits(C const& col, V const& v0, T const&,
                                ParticleType v0Type, bool ambVeto)
  {
    uint32_t bits = 0;
    auto posTrack = v0.template posTrack_as<T>();
    auto negTrack = v0.template negTrack_as<T>();

    // 0. Mass window (matched to v0Type)
    float mass = (v0Type == kLambda) ? v0.mLambda() : v0.mAntiLambda();
    if (mass > cMinV0Mass && mass < cMaxV0Mass)
      bits |= (1u << kCutMassWindow);

    // 1. Daughter PID (matched to v0Type)
    float nSPosPr = posTrack.tpcNSigmaPr(), nSNegPi = negTrack.tpcNSigmaPi();
    float nSPosPi = posTrack.tpcNSigmaPi(), nSNegPr = negTrack.tpcNSigmaPr();
    bool pidOk = (v0Type == kLambda)
        ? (std::abs(nSPosPr) < cTpcNsigmaCut && std::abs(nSNegPi) < cTpcNsigmaCut)
        : (std::abs(nSNegPr) < cTpcNsigmaCut && std::abs(nSPosPi) < cTpcNsigmaCut);
    if (pidOk) bits |= (1u << kCutDauPid);

    // 2. Daughter track quality
    if (selTrack(posTrack) && selTrack(negTrack))
      bits |= (1u << kCutDauTrackQual);

    // 3. Daughter DCA-to-PV
    float dcaProton = (v0Type == kLambda) ? std::abs(v0.dcapostopv()) : std::abs(v0.dcanegtopv());
    float dcaPion   = (v0Type == kLambda) ? std::abs(v0.dcanegtopv()) : std::abs(v0.dcapostopv());
    if (dcaProton >= cMinDcaProtonToPV && dcaPion >= cMinDcaPionToPV)
      bits |= (1u << kCutDauDcaToPV);

    // 4. Kinematic
    float rapVal = cDoEtaAnalysis ? std::abs(v0.eta()) : std::abs(v0.yLambda());
    if (kinCutSelection(v0.pt(), rapVal, cMinV0Pt, cMaxV0Pt, cMaxV0Rap))
      bits |= (1u << kCutKinematic);

    // 5. dcaV0Daughters window
    if (v0.dcaV0daughters() > cMinV0DcaDaughters && v0.dcaV0daughters() < cMaxV0DcaDaughters)
      bits |= (1u << kCutDcaV0Dau);

    // 6. dcav0topv window
    if (v0.dcav0topv() > cMinDcaV0ToPV && v0.dcav0topv() < cMaxDcaV0ToPV)
      bits |= (1u << kCutDcaV0ToPV);

    // 7. v0radius window
    if (v0.v0radius() > cMinV0TransRadius && v0.v0radius() < cMaxV0TransRadius)
      bits |= (1u << kCutV0Radius);

    // 8. ctau
    float ctau = v0.distovertotmom(col.posX(), col.posY(), col.posZ()) * MassLambda0;
    if (ctau > cMinV0CTau && ctau < cMaxV0CTau)
      bits |= (1u << kCutCtau);

    // 9. cosPA
    if (v0.v0cosPA() > cMinV0CosPA)
      bits |= (1u << kCutCosPA);

    // 10. K0s rejection (always passes when flag off)
    if (!cKshortRejFlag || std::abs(v0.mK0Short() - MassK0Short) > cKshortRejMassWindow)
      bits |= (1u << kCutK0sRej);

    // 11. Ambiguous-track veto (passed when veto disabled OR no ambiguity)
    if (ambVeto)
      bits |= (1u << kCutAmbiguousVeto);

    return bits;
  }

  template <typename V, typename T>
  bool hasAmbiguousDaughters(V const& v0, T const&)
  {
    auto posTrack = v0.template posTrack_as<T>();
    auto negTrack = v0.template negTrack_as<T>();

    auto posTrackCompCols = posTrack.compatibleCollIds();
    auto negTrackCompCols = negTrack.compatibleCollIds();

    // Check if daughter tracks belongs to more than one collision (Ambiguous Tracks)
    if (posTrackCompCols.size() > 1 || negTrackCompCols.size() > 1) {
      return true;
    }

    // Check if compatible collision index matches the track collision index
    if (((posTrackCompCols.size() != 0) && (posTrackCompCols[0] != posTrack.collisionId())) ||
        ((negTrackCompCols.size() != 0) && (negTrackCompCols[0] != negTrack.collisionId()))) {
      return true;
    }

    // Pass as not ambiguous
    return false;
  }

  template <typename V>
  PrmScdType isPrimaryV0(V const& v0)
  {
    auto mcpart = v0.template mcParticle_as<aod::McParticles>();

    // check for secondary lambda
    if (!mcpart.isPhysicalPrimary()) {
      histos.fill(HIST("Tracks/h1f_tracks_info"), kSecondaryLambda);
      return kSecondary;
    }

    histos.fill(HIST("Tracks/h1f_tracks_info"), kPrimaryLambda);
    return kPrimary;
  }

  template <typename V, typename T>
  bool selTrueMcRecLambda(V const& v0, T const&)
  {
    auto mcpart = v0.template mcParticle_as<aod::McParticles>();

    // check if Lambda/AntiLambda
    if (std::abs(mcpart.pdgCode()) != kLambda0) {
      return false;
    }

    // Check for daughters
    if (cCheckRecoDauFlag) {
      auto postrack = v0.template posTrack_as<T>();
      auto negtrack = v0.template negTrack_as<T>();

      // check if the daughters have corresponding mcparticle
      if (!postrack.has_mcParticle() || !negtrack.has_mcParticle()) {
        histos.fill(HIST("Tracks/h1f_tracks_info"), kLambdaDauNotMcParticle);
        return false;
      }

      auto mcpostrack = postrack.template mcParticle_as<aod::McParticles>();
      auto mcnegtrack = negtrack.template mcParticle_as<aod::McParticles>();

      if (mcpart.pdgCode() == kLambda0) {
        if (mcpostrack.pdgCode() != kProton || mcnegtrack.pdgCode() != kPiMinus) {
          histos.fill(HIST("Tracks/h1f_tracks_info"), kLambdaNotPrPiMinus);
          return false;
        }
      } else if (mcpart.pdgCode() == kLambda0Bar) {
        if (mcpostrack.pdgCode() != kPiPlus || mcnegtrack.pdgCode() != kProtonBar) {
          histos.fill(HIST("Tracks/h1f_tracks_info"), kAntiLambdaNotAntiPrPiPlus);
          return false;
        }
      }
    }

    return true;
  }

  template <ParticleType part, typename V>
  float getCorrectionFactors(V const& v0)
  {
    // Check for efficiency correction flag
    if (!cCorrectionFlag) {
      return 1.;
    }

    // Get  from CCDB
    auto ccdbObj = ccdb->getForTimeStamp<TList>(cPathCCDB.value, 1);

    // Check CCDB Object
    if (!ccdbObj) {
      LOGF(warning, "CCDB OBJECT NOT FOUND");
      return 1.;
    }

    // initialize efficiency factor and primary fraction values
    float effCorrFact = 1., primFrac = 1.;
    float rap = (cDoEtaAnalysis) ? v0.eta() : v0.yLambda();

    // Get Efficiency Factor
    if (cGetEffFact) {
      // [Phase 16o] Guard against nullptr from FindObject — the named
      // histogram may be missing from the CCDB TList. Previously the
      // subsequent ->Clone() dereferenced nullptr and crashed.
      const auto effName = Form("%s", vCorrFactStrings[cCorrFactHist][part].c_str());
      TObject* objEff = reinterpret_cast<TObject*>(ccdbObj->FindObject(effName));
      if (!objEff) {
        LOGF(warning, "[CCDB] Efficiency histogram '%s' not found; using effCorrFact=1.0", effName);
        effCorrFact = 1.f;
      } else {
        TH1F* histEff = reinterpret_cast<TH1F*>(objEff->Clone());
      if (histEff->GetDimension() == TwoDimCorr) {
        histos.fill(HIST("Tracks/h1f_tracks_info"), kEffCorrPtCent);
        effCorrFact = histEff->GetBinContent(histEff->FindBin(cent, v0.pt()));
      } else if (histEff->GetDimension() == ThreeDimCorr) {
        histos.fill(HIST("Tracks/h1f_tracks_info"), kEffCorrPtRapCent);
        effCorrFact = histEff->GetBinContent(histEff->FindBin(cent, v0.pt(), rap));
      } else {
        histos.fill(HIST("Tracks/h1f_tracks_info"), kNoEffCorr);
        LOGF(warning, "CCDB OBJECT IS NOT A HISTOGRAM !!!");
        effCorrFact = 1.;
      }
        delete histEff;
      }
    }

    // Get Primary Fraction
    // (The dimension of this could be different than efficiency because of large errors !!!)
    if (cGetPrimFrac) {
      // [Phase 16o] Same null-guard as the efficiency lookup above.
      const auto pfName = Form("%s", vPrimFracStrings[cPrimFracHist][part].c_str());
      TObject* objPrm = reinterpret_cast<TObject*>(ccdbObj->FindObject(pfName));
      if (!objPrm) {
        LOGF(warning, "[CCDB] Primary-fraction histogram '%s' not found; using primFrac=1.0", pfName);
        primFrac = 1.f;
      } else {
        TH1F* histPrm = reinterpret_cast<TH1F*>(objPrm->Clone());
      if (histPrm->GetDimension() == TwoDimCorr) {
        histos.fill(HIST("Tracks/h1f_tracks_info"), kPFCorrPtCent);
        primFrac = histPrm->GetBinContent(histPrm->FindBin(cent, v0.pt()));
      } else if (histPrm->GetDimension() == ThreeDimCorr) {
        histos.fill(HIST("Tracks/h1f_tracks_info"), kPFCorrPtRapCent);
        primFrac = histPrm->GetBinContent(histPrm->FindBin(cent, v0.pt(), rap));
      } else {
        histos.fill(HIST("Tracks/h1f_tracks_info"), kNoPFCorr);
        LOGF(warning, "CCDB OBJECT IS NOT A HISTOGRAM !!!");
        primFrac = 1.;
      }
        delete histPrm;
      }
    }

    return primFrac * effCorrFact;
  }

  // [Phase 7] Truth-level mother PDG of a (possibly truth-matched) reco V0.
  // Returns 0 on data and on MC V0s that have no MC particle / no mother.
  // Used to populate aod::lambdatrack::MotherPdg so downstream consumers can
  // tag Ξ feed-down (motherPdg == ±3312), Ω feed-down (±3334), Σ⁰ resonance
  // feed-down (±3212), etc. without rerunning truth matching.
  template <typename V>
  int getLambdaMotherPdg(V const& v0)
  {
    if (!v0.has_mcParticle())
      return 0;
    auto mcpart = v0.template mcParticle_as<aod::McParticles>();
    if (!mcpart.has_mothers())
      return 0;
    return mcpart.template mothers_first_as<aod::McParticles>().pdgCode();
  }

  template <typename V, typename T>
  void fillLambdaMothers(V const& v0, T const&)
  {
    auto mcpart = v0.template mcParticle_as<aod::McParticles>();
    // [Phase 16p] Defensive — a "secondary" V0 should always have ≥1 mother
    // by definition, but some AOD generators occasionally strip mother links.
    // Without this guard, lambdaMothers[0] dereferences an empty range → UB.
    if (!mcpart.has_mothers())
      return;
    auto lambdaMothers = mcpart.template mothers_as<aod::McParticles>();
    histos.fill(HIST("Tracks/h2f_lambda_mothers_pdg"), lambdaMothers[0].pdgCode(), v0.pt());
  }

  template <ParticleType part, typename C, typename V, typename T>
  void fillLambdaQAHistos(C const& col, V const& v0, T const&)
  {
    static constexpr std::string_view SubDir[] = {"QA/Lambda/", "QA/AntiLambda/"};

    // daugthers
    auto postrack = v0.template posTrack_as<T>();
    auto negtrack = v0.template negTrack_as<T>();
    float mass = 0.;

    if constexpr (part == kLambda) {
      mass = v0.mLambda();
    } else {
      mass = v0.mAntiLambda();
    }

    // ctau
    float e = RecoDecay::e(v0.px(), v0.py(), v0.pz(), mass);
    float gamma = e / mass;
    float ctau = v0.distovertotmom(col.posX(), col.posY(), col.posZ()) * MassLambda0;
    float gctau = ctau * gamma;

    histos.fill(HIST(SubDir[part]) + HIST("h2f_qt_vs_alpha"), v0.alpha(), v0.qtarm());
    histos.fill(HIST(SubDir[part]) + HIST("h1f_dca_V0_daughters"), v0.dcaV0daughters());
    histos.fill(HIST(SubDir[part]) + HIST("h1f_dca_pos_to_PV"), v0.dcapostopv());
    histos.fill(HIST(SubDir[part]) + HIST("h1f_dca_neg_to_PV"), v0.dcanegtopv());
    histos.fill(HIST(SubDir[part]) + HIST("h1f_dca_V0_to_PV"), v0.dcav0topv());
    histos.fill(HIST(SubDir[part]) + HIST("h1f_V0_cospa"), v0.v0cosPA());
    histos.fill(HIST(SubDir[part]) + HIST("h1f_V0_radius"), v0.v0radius());
    histos.fill(HIST(SubDir[part]) + HIST("h1f_V0_ctau"), ctau);
    histos.fill(HIST(SubDir[part]) + HIST("h1f_V0_gctau"), gctau);

    histos.fill(HIST(SubDir[part]) + HIST("h1f_pos_prong_pt"), postrack.pt());
    histos.fill(HIST(SubDir[part]) + HIST("h1f_pos_prong_eta"), postrack.eta());
    histos.fill(HIST(SubDir[part]) + HIST("h1f_pos_prong_phi"), postrack.phi());
    histos.fill(HIST(SubDir[part]) + HIST("h1f_neg_prong_pt"), negtrack.pt());
    histos.fill(HIST(SubDir[part]) + HIST("h1f_neg_prong_eta"), negtrack.eta());
    histos.fill(HIST(SubDir[part]) + HIST("h1f_neg_prong_phi"), negtrack.phi());

    histos.fill(HIST(SubDir[part]) + HIST("h2f_pos_prong_dcaXY_vs_pt"), postrack.pt(), postrack.dcaXY());
    histos.fill(HIST(SubDir[part]) + HIST("h2f_neg_prong_dcaXY_vs_pt"), negtrack.pt(), negtrack.dcaXY());
    histos.fill(HIST(SubDir[part]) + HIST("h2f_pos_prong_dEdx_vs_p"), postrack.tpcInnerParam(), postrack.tpcSignal());
    histos.fill(HIST(SubDir[part]) + HIST("h2f_neg_prong_dEdx_vs_p"), negtrack.tpcInnerParam(), negtrack.tpcSignal());
    histos.fill(HIST(SubDir[part]) + HIST("h2f_pos_prong_tpc_nsigma_pr_vs_p"), postrack.tpcInnerParam(), postrack.tpcNSigmaPr());
    histos.fill(HIST(SubDir[part]) + HIST("h2f_neg_prong_tpc_nsigma_pr_vs_p"), negtrack.tpcInnerParam(), negtrack.tpcNSigmaPr());
    histos.fill(HIST(SubDir[part]) + HIST("h2f_pos_prong_tpc_nsigma_pi_vs_p"), postrack.tpcInnerParam(), postrack.tpcNSigmaPi());
    histos.fill(HIST(SubDir[part]) + HIST("h2f_neg_prong_tpc_nsigma_pi_vs_p"), negtrack.tpcInnerParam(), negtrack.tpcNSigmaPi());
  }

  // Fill Lambda Kinematic Histograms
  template <RecGenType rg, ParticleType part>
  void fillKinematicHists(float const& pt, float const& eta, float const& y, float const& phi)
  {
    static constexpr std::string_view SubDirRG[] = {"McRec/", "McGen/"};
    static constexpr std::string_view SubDirPart[] = {"Lambda/", "AntiLambda/"};

    histos.fill(HIST(SubDirRG[rg]) + HIST(SubDirPart[part]) + HIST("hPt"), pt);
    histos.fill(HIST(SubDirRG[rg]) + HIST(SubDirPart[part]) + HIST("hEta"), eta);
    histos.fill(HIST(SubDirRG[rg]) + HIST(SubDirPart[part]) + HIST("hRap"), y);
    histos.fill(HIST(SubDirRG[rg]) + HIST(SubDirPart[part]) + HIST("hPhi"), phi);
  }

  // Reconstructed Level Tables
  template <RunType run, DMCType dmc, typename C, typename V, typename T>
  void fillLambdaRecoTables(C const& collision, V const& v0tracks, T const& tracks)
  {
    // Total Collisions
    histos.fill(HIST("Events/h1f_collisions_info"), kTotCol);

    // Select Collision (Only for Data... McRec has been selected already !!!)
    if constexpr (dmc == kData) {
      if (!selCollision<run>(collision)) {
        return;
      }
    }

    histos.fill(HIST("Events/h1f_collisions_info"), kPassSelCol);
    histos.fill(HIST("Events/h1f_collision_posZ"), collision.posZ());

    // Fill Collision Table
    // lambdaCollisionTable(cent, mult, collision.posX(), collision.posY(), collision.posZ());
    lambdaCollisionTable(cent, mult, collision.globalIndex(), collision.posX(), collision.posY(), collision.posZ());

    // initialize v0track objects
    ParticleType v0Type = kLambda;
    PrmScdType v0PrmScdType = kPrimary;
    float mass = 0., corr_fact = 1.;

    for (auto const& v0 : v0tracks) {
      // check for corresponding MCGen Particle
      if constexpr (dmc == kMC) {
        histos.fill(HIST("Tracks/h1f_tracks_info"), kTracksBeforeHasMcParticle);
        if (!v0.has_mcParticle()) {
          continue;
        }
      }

      histos.fill(HIST("Tracks/h1f_tracks_info"), kAllV0Tracks);
      histos.fill(HIST("Tracks/h2f_armpod_before_sel"), v0.alpha(), v0.qtarm());

      // Select V0 Particle as Lambda/AntiLambda
      if (!selV0Particle(collision, v0, tracks, v0Type)) {
        continue;
      }

      // Select V0 Type Selection
      if (cV0TypeSelFlag && v0.v0Type() != cV0TypeSelection) {
        continue;
      }

      // we have v0 as lambda
      histos.fill(HIST("Tracks/h1f_tracks_info"), kAllSelPassed);

      // Remove lambda with ambiguous daughters (Only for run3)
      if constexpr (run == kRun3) {
        if (cRemoveAmbiguousTracks && hasAmbiguousDaughters(v0, tracks)) {
          continue;
        }
      }

      // Get Lambda mass and kinematic variables
      mass = (v0Type == kLambda) ? v0.mLambda() : v0.mAntiLambda();
      pt = v0.pt();
      eta = v0.eta();
      rap = v0.yLambda();
      phi = v0.phi();

      // [Phase 7] Mother PDG: 0 on data, parent PDG on MC.
      int motherPdg = 0;

      // do MC analysis
      if constexpr (dmc == kMC) {
        histos.fill(HIST("Tracks/h2f_tracks_pid_before_sel"), v0.mcParticle().pdgCode(), v0.pt());

        // [Phase 16l] Always tag primary/secondary from MC truth on MC reco.
        // Previously gated by cSelMCPSV0 — when disabled, v0PrmScdType
        // stayed kPrimary for ALL candidates, which silently broke the
        // downstream goodPrimaryLambda partition (it became equivalent
        // to goodLambda, accepting all feed-down Λs).
        // cSelMCPSV0 is kept as a no-op Configurable for JSON back-compat;
        // the producer now unconditionally annotates v0PrmScd. Users who
        // want "treat all as primary" should toggle the correlator-side
        // cUsePrimaryLambdasOnly instead.
        v0PrmScdType = isPrimaryV0(v0);

        // check for true Lambda/Anti-Lambda
        if (cSelectTrueLambda && !selTrueMcRecLambda(v0, tracks)) {
          continue;
        }

        // get mothers information (also populated as table column below)
        motherPdg = getLambdaMotherPdg(v0);
        if (v0PrmScdType == kSecondary) {
          fillLambdaMothers(v0, tracks);
        }

        histos.fill(HIST("Tracks/h1f_tracks_info"), kPassTrueLambdaSel);
        histos.fill(HIST("Tracks/h2f_tracks_pid_after_sel"), v0.mcParticle().pdgCode(), v0.pt());

        if (cRecoMomResoFlag) {
          auto mc = v0.template mcParticle_as<aod::McParticles>();
          pt = mc.pt();
          eta = mc.eta();
          rap = mc.y();
          phi = mc.phi();
          float y = (cDoEtaAnalysis) ? eta : rap;
          // apply kinematic selection (On Truth)
          if (!kinCutSelection(pt, std::abs(y), cMinV0Pt, cMaxV0Pt, cMaxV0Rap)) {
            continue;
          }
        }
      }

      histos.fill(HIST("Tracks/h2f_armpod_after_sel"), v0.alpha(), v0.qtarm());

      // get correction factors
      corr_fact = (v0Type == kLambda) ? getCorrectionFactors<kLambda>(v0) : getCorrectionFactors<kAntiLambda>(v0);

      // fill lambda qa
      if (v0Type == kLambda) {
        histos.fill(HIST("Tracks/h1f_lambda_pt_vs_invm"), mass, v0.pt());
        fillLambdaQAHistos<kLambda>(collision, v0, tracks);
        fillKinematicHists<kRec, kLambda>(v0.pt(), v0.eta(), v0.yLambda(), v0.phi());
      } else {
        histos.fill(HIST("Tracks/h1f_antilambda_pt_vs_invm"), mass, v0.pt());
        fillLambdaQAHistos<kAntiLambda>(collision, v0, tracks);
        fillKinematicHists<kRec, kAntiLambda>(v0.pt(), v0.eta(), v0.yLambda(), v0.phi());
      }

      // [Phase 8] Per-V0 topology snapshot — written on every Λ row so
      // (i) the optional Λ TTree carries the variables for a downstream
      // template fit, (ii) the partition can gate on passesPrimaryTopo
      // on data without re-resolving daughter tracks at correlator time.
      auto posTrk = v0.template posTrack_as<T>();
      auto negTrk = v0.template negTrack_as<T>();
      float dcaV0ToPV_v   = v0.dcav0topv();        // V0 line vs PV
      float v0Radius_v    = v0.v0radius();         // transverse decay radius
      int8_t posItsNCls_v = static_cast<int8_t>(posTrk.itsNCls());
      int8_t negItsNCls_v = static_cast<int8_t>(negTrk.itsNCls());
      // [Phase 9] ITS hit-map per daughter (uint8_t bitmask, bit i = layer i).
      uint8_t posItsClusterMap_v = static_cast<uint8_t>(posTrk.itsClusterMap());
      uint8_t negItsClusterMap_v = static_cast<uint8_t>(negTrk.itsClusterMap());
      // [Phase 10] Per-daughter signed DCA-XY-to-PV (V0Datas already exposes
      // these as helices propagated to the PV).
      float posDcaXY_v = v0.dcapostopv();
      float negDcaXY_v = v0.dcanegtopv();
      // [Phase 10] Pseudo-proper transverse decay length L_proper.
      // Guard against pT==0 to avoid div-by-zero for pathological rows.
      float lProper_v = (v0.pt() > 0.f)
                          ? (v0Radius_v * static_cast<float>(MassLambda0) / v0.pt())
                          : 0.f;

      // [Phase 14] Raw cut-input variables for offline re-cutting.
      float tpcNSigmaPosPr_v = posTrk.tpcNSigmaPr();
      float tpcNSigmaNegPi_v = negTrk.tpcNSigmaPi();
      float tpcNSigmaPosPi_v = posTrk.tpcNSigmaPi();
      float tpcNSigmaNegPr_v = negTrk.tpcNSigmaPr();
      float mK0Short_v       = v0.mK0Short();
      float qtArm_v          = v0.qtarm();
      float alphaArm_v       = v0.alpha();
      float cTau_v           = v0.distovertotmom(collision.posX(),
                                                 collision.posY(),
                                                 collision.posZ()) * MassLambda0;
      // [Phase 16a] When the producer-side gate is disabled, default to
      // `true` so the column means "no opinion / everything passes". The
      // previous default of `false` silently emptied the partition when the
      // correlator's cPrimaryRequireTopo was flipped on without enabling
      // cPrimEnable here — a confusing zero-pair failure mode.
      bool passesPrim     = !primCfg.cPrimEnable;
      if (primCfg.cPrimEnable) {
        bool topoOk =
            std::abs(dcaV0ToPV_v)            < primCfg.cPrimMaxDcaV0ToPv  &&
            v0.v0cosPA()                     > primCfg.cPrimMinV0CosPA   &&
            v0Radius_v                       < primCfg.cPrimMaxV0Radius  &&
            std::abs(posDcaXY_v)             < primCfg.cPrimMaxDauDcaToPv &&
            std::abs(negDcaXY_v)             < primCfg.cPrimMaxDauDcaToPv &&
            lProper_v                        < primCfg.cPrimMaxLProper;
        bool itsOk = primCfg.cPrimRequireBothDauItsHits
                         ? (posItsNCls_v >= primCfg.cPrimMinDauItsNCls && negItsNCls_v >= primCfg.cPrimMinDauItsNCls)
                         : (posItsNCls_v >= primCfg.cPrimMinDauItsNCls || negItsNCls_v >= primCfg.cPrimMinDauItsNCls);
        // [Phase 10] ITS-IB requirement (bits 0,1,2 = Layers 0,1,2).
        // OR'd across both daughters: at least ONE daughter must have an
        // IB hit. Optional via the Configurable.
        bool itsIBOk = !primCfg.cPrimRequireItsIBHit
                         || (((posItsClusterMap_v | negItsClusterMap_v) & lcorr_const::kItsIBMask) != 0);
        passesPrim = topoOk && itsOk && itsIBOk;
      }

      // [Phase 14] Per-V0 cut bitmask. Bits 0-11 reflect the standalone-V0
      // selection chain; bit 12 reflects MC truth (always 1 on data; set
      // post-hoc on MC reco); bit 13 reflects passesPrimaryTopo. Computed
      // inline so the row carries the full diagnostic.
      // hasAmbiguousDaughters() requires aod::TrackCompColls which is only
      // joined in the Run3 Tracks alias, so guard with `if constexpr`.
      bool ambVetoOk = true;
      if constexpr (run == kRun3) {
        ambVetoOk = !cRemoveAmbiguousTracks || !hasAmbiguousDaughters(v0, tracks);
      }
      uint32_t cutBits_v = computeLambdaCutBits(collision, v0, tracks, v0Type, ambVetoOk);
      // [Phase 16q] Bit 12 (MC truth) — set independently of cSelectTrueLambda.
      // On data: always set (truth is N/A; bit means "no MC reason to drop").
      // On MC reco: set iff selTrueMcRecLambda actually passes. Previously
      // this was set unconditionally on the assumption that cSelectTrueLambda
      // had already filtered the survivors, but when cSelectTrueLambda=false
      // the filter never runs and the bit became unreliable.
      bool isTrueMcLam = true;
      if constexpr (dmc == kMC) {
        isTrueMcLam = selTrueMcRecLambda(v0, tracks);
      }
      if (isTrueMcLam)
        cutBits_v |= (1u << kCutMcTrueLambda);
      if (passesPrim)
        cutBits_v |= (1u << kCutPhase10Prim);

      // Fill Λ/Λ̄ row — [P7] motherPdg · [P8] topo+flag · [P9] ITS maps · [P10] L_proper + dau DCAs · [P14] cutBits + raw inputs
      lambdaTrackTable(lambdaCollisionTable.lastIndex(), v0.px(), v0.py(), v0.pz(),
                       pt, eta, phi, rap, mass, posTrk.index(), negTrk.index(),
                       v0.v0cosPA(), v0.dcaV0daughters(), (int8_t)v0Type, v0PrmScdType, corr_fact, motherPdg,
                       dcaV0ToPV_v, v0Radius_v, posItsNCls_v, negItsNCls_v, passesPrim,
                       posItsClusterMap_v, negItsClusterMap_v,
                       lProper_v, posDcaXY_v, negDcaXY_v,
                       cutBits_v,
                       tpcNSigmaPosPr_v, tpcNSigmaNegPi_v, tpcNSigmaPosPi_v, tpcNSigmaNegPr_v,
                       mK0Short_v, qtArm_v, alphaArm_v, cTau_v);
    }
  }

  // MC Generater Level Tables
  template <RunType run, typename C, typename M>
  void fillLambdaMcGenTables(C const& mcCollision, M const& mcParticles)
  {
    // Fill McGen Collision Table
    lambdaMCGenCollisionTable(cent, mult, mcCollision.globalIndex(), mcCollision.posX(), mcCollision.posY(), mcCollision.posZ());

    // initialize track objects
    ParticleType v0Type = kLambda;
    PrmScdType v0PrmScdType = kPrimary;
    float rap = 0.;

    for (auto const& mcpart : mcParticles) {
      // check for Lambda first
      if (mcpart.pdgCode() == kLambda0) {
        v0Type = kLambda;
      } else if (mcpart.pdgCode() == kLambda0Bar) {
        v0Type = kAntiLambda;
      } else {
        continue;
      }

      // check for Primary Lambda/AntiLambda
      if (mcpart.isPhysicalPrimary()) {
        v0PrmScdType = kPrimary;
      } else {
        v0PrmScdType = kSecondary;
      }

      // Decide Eta/Rap
      if (!cDoEtaAnalysis) {
        rap = mcpart.y();
      } else {
        rap = mcpart.eta();
      }

      // Apply Kinematic Acceptance
      if (!kinCutSelection(mcpart.pt(), std::abs(rap), cMinV0Pt, cMaxV0Pt, cMaxV0Rap)) {
        continue;
      }

      histos.fill(HIST("Tracks/h1f_tracks_info"), kGenTotAccLambda);

      // get daughter track info and check for decay channel flag
      if (!mcpart.has_daughters()) {
        histos.fill(HIST("Tracks/h1f_tracks_info"), kGenLambdaNoDau);
        continue;
      }
      auto dautracks = mcpart.template daughters_as<aod::McParticles>();
      std::vector<int> daughterPDGs, daughterIDs;
      std::vector<float> vDauPt, vDauEta, vDauRap, vDauPhi;
      for (auto const& dautrack : dautracks) {
        daughterPDGs.push_back(dautrack.pdgCode());
        daughterIDs.push_back(dautrack.globalIndex());
        vDauPt.push_back(dautrack.pt());
        vDauEta.push_back(dautrack.eta());
        vDauRap.push_back(dautrack.y());
        vDauPhi.push_back(dautrack.phi());
      }
      if (cGenDecayChannel) { // check decay channel
        if (v0Type == kLambda) {
          if (daughterPDGs[0] != kProton || daughterPDGs[1] != kPiMinus) {
            continue;
          }
        } else if (v0Type == kAntiLambda) {
          if (daughterPDGs[0] != kProtonBar || daughterPDGs[1] != kPiPlus) {
            continue;
          }
        }
      }

      histos.fill(HIST("Tracks/h1f_tracks_info"), kGenLambdaToPrPi);

      if (v0Type == kLambda) {
        histos.fill(HIST("McGen/h1f_lambda_daughter_PDG"), daughterPDGs[0]);
        histos.fill(HIST("McGen/h1f_lambda_daughter_PDG"), daughterPDGs[1]);
        histos.fill(HIST("McGen/h1f_lambda_daughter_PDG"), mcpart.pdgCode());
        histos.fill(HIST("McGen/Lambda/Proton/hPt"), vDauPt[0]);
        histos.fill(HIST("McGen/Lambda/Proton/hEta"), vDauEta[0]);
        histos.fill(HIST("McGen/Lambda/Proton/hRap"), vDauRap[0]);
        histos.fill(HIST("McGen/Lambda/Proton/hPhi"), vDauPhi[0]);
        histos.fill(HIST("McGen/Lambda/Pion/hPt"), vDauPt[1]);
        histos.fill(HIST("McGen/Lambda/Pion/hEta"), vDauEta[1]);
        histos.fill(HIST("McGen/Lambda/Pion/hRap"), vDauRap[1]);
        histos.fill(HIST("McGen/Lambda/Pion/hPhi"), vDauPhi[1]);
        fillKinematicHists<kGen, kLambda>(mcpart.pt(), mcpart.eta(), mcpart.y(), mcpart.phi());
      } else {
        histos.fill(HIST("McGen/h1f_antilambda_daughter_PDG"), daughterPDGs[0]);
        histos.fill(HIST("McGen/h1f_antilambda_daughter_PDG"), daughterPDGs[1]);
        histos.fill(HIST("McGen/h1f_antilambda_daughter_PDG"), mcpart.pdgCode());
        histos.fill(HIST("McGen/AntiLambda/Pion/hPt"), vDauPt[0]);
        histos.fill(HIST("McGen/AntiLambda/Pion/hEta"), vDauEta[0]);
        histos.fill(HIST("McGen/AntiLambda/Pion/hRap"), vDauRap[0]);
        histos.fill(HIST("McGen/AntiLambda/Pion/hPhi"), vDauPhi[0]);
        histos.fill(HIST("McGen/AntiLambda/Proton/hPt"), vDauPt[1]);
        histos.fill(HIST("McGen/AntiLambda/Proton/hEta"), vDauEta[1]);
        histos.fill(HIST("McGen/AntiLambda/Proton/hRap"), vDauRap[1]);
        histos.fill(HIST("McGen/AntiLambda/Proton/hPhi"), vDauPhi[1]);
        fillKinematicHists<kGen, kAntiLambda>(mcpart.pt(), mcpart.eta(), mcpart.y(), mcpart.phi());
      }

      // Fill Lambda McGen Table
      lambdaMCGenTrackTable(lambdaMCGenCollisionTable.lastIndex(), mcpart.px(), mcpart.py(), mcpart.pz(),
                            mcpart.pt(), mcpart.eta(), mcpart.phi(), mcpart.y(), RecoDecay::m(mcpart.p(), mcpart.e()),
                            daughterIDs[0], daughterIDs[1], (int8_t)v0Type, -999., -999., v0PrmScdType, 1.);
    }
  }

  template <RunType run, DMCType dmc, typename M, typename C, typename V, typename T, typename P>
  void analyzeMcRecoGen(M const& mcCollision, C const& collisions, V const& V0s, T const& tracks, P const& mcParticles)
  {
    // Number of Rec Collisions Associated to the McGen Collision
    int nRecCols = collisions.size();
    if (nRecCols != 0) {
      histos.fill(HIST("McGen/h1f_collision_recgen"), nRecCols);
    }

    // Always fill gen tables so processMCGenXi/Omega has entries.
    // selCollision sets cent/mult as a side-effect — works for both Run2 and Run3.
    cent = 0.f;
    mult = 0.f;
    if (nRecCols >= 1 &&
        collisions.begin().has_mcCollision() &&
        collisions.begin().mcCollisionId() == mcCollision.globalIndex()) {
      selCollision<run>(collisions.begin()); // sets cent and mult
    }
    fillLambdaMcGenTables<run>(mcCollision, mcParticles);

    // Reco tables only for clean 1-to-1 matched collisions
    if (nRecCols != 1) {
      return;
    }
    histos.fill(HIST("McGen/h1f_collisions_info"), kTotCol);
    if (!collisions.begin().has_mcCollision() || !selCollision<run>(collisions.begin()) || collisions.begin().mcCollisionId() != mcCollision.globalIndex()) {
      return;
    }
    histos.fill(HIST("McGen/h1f_collisions_info"), kPassSelCol);
    histos.fill(HIST("McGen/h2f_collision_posZ"), mcCollision.posZ(), collisions.begin().posZ());
    auto v0Tracks = V0s.sliceBy(perCollision, collisions.begin().globalIndex());
    fillLambdaRecoTables<run, dmc>(collisions.begin(), v0Tracks, tracks);
  }

  SliceCache cache;
  Preslice<soa::Join<aod::V0Datas, aod::McV0Labels>> perCollision = aod::v0data::collisionId;

  // ===========================================================================
  // [Phase 4] Cascade-side helpers (verbatim from former CascadeSelector,
  // s/registry/cascRegistry/g).
  // ===========================================================================

  // [Phase 4 fix] Templated on TTracks so the call site can pass whichever
  // track type is bound in its process function (Tracks for data,
  // TracksMC for MC reco). Avoids hardcoding FullTracksExtIUWithPID, which
  // is TracksIU-based and incompatible with our merged Tracks-based binding.
  template <typename TTracks, typename TCollision>
  void fillMatchedHistos(LabeledCascades::iterator rec, int flag, TCollision collision)
  {
    if (flag == 0)
      return;
    if (!rec.has_mcParticle())
      return;
    auto gen = rec.mcParticle();
    if (!gen.isPhysicalPrimary())
      return;
    int genpdg = gen.pdgCode();
    // [Phase 16j] flag<3 (== 1 or 2) means Ξ-eligible; flag>1 (== 2 or 3)
    // means Ω-eligible. Spelled out with named constants.
    const bool xiLike = (flag == lcorr_const::kFlagXiOnly || flag == lcorr_const::kFlagXiAndOmega);
    const bool omLike = (flag == lcorr_const::kFlagXiAndOmega || flag == lcorr_const::kFlagOmegaOnly);
    if ((xiLike && std::abs(genpdg) == lcorr_const::kXiMinusPdg) ||
        (omLike && std::abs(genpdg) == lcorr_const::kOmegaMinusPdg)) {
      cascRegistry.fill(HIST("truerec/hV0Radius"), rec.v0radius());
      cascRegistry.fill(HIST("truerec/hCascRadius"), rec.cascradius());
      cascRegistry.fill(HIST("truerec/hV0CosPA"), rec.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
      cascRegistry.fill(HIST("truerec/hCascCosPA"), rec.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
      cascRegistry.fill(HIST("truerec/hDCAPosToPV"), rec.dcapostopv());
      cascRegistry.fill(HIST("truerec/hDCANegToPV"), rec.dcanegtopv());
      cascRegistry.fill(HIST("truerec/hDCABachToPV"), rec.dcabachtopv());
      cascRegistry.fill(HIST("truerec/hDCAV0ToPV"), rec.dcav0topv(collision.posX(), collision.posY(), collision.posZ()));
      cascRegistry.fill(HIST("truerec/hDCAV0Dau"), rec.dcaV0daughters());
      cascRegistry.fill(HIST("truerec/hDCACascDau"), rec.dcacascdaughters());
      cascRegistry.fill(HIST("truerec/hLambdaMass"), rec.mLambda());
      cascRegistry.fill(HIST("truerec/hITSnClustersPos"), rec.template posTrack_as<TTracks>().itsNCls());
      cascRegistry.fill(HIST("truerec/hITSnClustersNeg"), rec.template negTrack_as<TTracks>().itsNCls());
      cascRegistry.fill(HIST("truerec/hITSnClustersBach"), rec.template bachelor_as<TTracks>().itsNCls());
      cascRegistry.fill(HIST("truerec/hTPCnCrossedRowsPos"), rec.template posTrack_as<TTracks>().tpcNClsCrossedRows());
      cascRegistry.fill(HIST("truerec/hTPCnCrossedRowsNeg"), rec.template negTrack_as<TTracks>().tpcNClsCrossedRows());
      cascRegistry.fill(HIST("truerec/hTPCnCrossedRowsBach"), rec.template bachelor_as<TTracks>().tpcNClsCrossedRows());
      cascRegistry.fill(HIST("truerec/hITSChi2Pos"), rec.template posTrack_as<TTracks>().itsChi2NCl());
      cascRegistry.fill(HIST("truerec/hITSChi2Neg"), rec.template negTrack_as<TTracks>().itsChi2NCl());
      cascRegistry.fill(HIST("truerec/hITSChi2Bach"), rec.template bachelor_as<TTracks>().itsChi2NCl());
      cascRegistry.fill(HIST("truerec/hTPCChi2Pos"), rec.template posTrack_as<TTracks>().tpcChi2NCl());
      cascRegistry.fill(HIST("truerec/hTPCChi2Neg"), rec.template negTrack_as<TTracks>().tpcChi2NCl());
      cascRegistry.fill(HIST("truerec/hTPCChi2Bach"), rec.template bachelor_as<TTracks>().tpcChi2NCl());
      switch (genpdg) {
        case lcorr_const::kXiMinusPdg:
          cascRegistry.fill(HIST("truerec/hXiMinus"), rec.pt(), rec.yXi());
          break;
        case -lcorr_const::kXiMinusPdg:
          cascRegistry.fill(HIST("truerec/hXiPlus"), rec.pt(), rec.yXi());
          break;
        case lcorr_const::kOmegaMinusPdg:
          cascRegistry.fill(HIST("truerec/hOmegaMinus"), rec.pt(), rec.yOmega());
          break;
        case -lcorr_const::kOmegaMinusPdg:
          cascRegistry.fill(HIST("truerec/hOmegaPlus"), rec.pt(), rec.yOmega());
          break;
      }
    }
  }

  // [Phase 4 fix] Templated on TTracks so the call site picks the bound track
  // type (Tracks for data, TracksMC for MC reco). Eliminates the framework's
  // multi-track-table index resolution conflict.
  template <typename TTracks, typename TCascade, typename TCollision>
  int processCandidate(TCascade const& casc, TCollision const& collision)
  {
    auto bachTrack = casc.template bachelor_as<TTracks>();
    auto posTrack = casc.template posTrack_as<TTracks>();
    auto negTrack = casc.template negTrack_as<TTracks>();

    cascRegistry.fill(HIST("hV0Radius"), casc.v0radius(), casc.mXi(), casc.pt());
    cascRegistry.fill(HIST("hCascRadius"), casc.cascradius(), casc.mXi(), casc.pt());
    cascRegistry.fill(HIST("hV0CosPA"), casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()), casc.mXi(), casc.pt());
    cascRegistry.fill(HIST("hCascCosPA"), casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()), casc.mXi(), casc.pt());
    cascRegistry.fill(HIST("hDCAPosToPV"), casc.dcapostopv(), casc.mXi(), casc.pt());
    cascRegistry.fill(HIST("hDCANegToPV"), casc.dcanegtopv(), casc.mXi(), casc.pt());
    cascRegistry.fill(HIST("hDCABachToPV"), casc.dcabachtopv(), casc.mXi(), casc.pt());
    cascRegistry.fill(HIST("hDCAV0ToPV"), casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()), casc.mXi(), casc.pt());
    cascRegistry.fill(HIST("hDCAV0Dau"), casc.dcaV0daughters(), casc.mXi(), casc.pt());
    cascRegistry.fill(HIST("hDCACascDau"), casc.dcacascdaughters(), casc.mXi(), casc.pt());
    cascRegistry.fill(HIST("hLambdaMass"), casc.mLambda(), casc.mXi(), casc.pt());

    cascRegistry.fill(HIST("hITSnClustersPos"), posTrack.itsNCls(), casc.mXi(), casc.pt());
    cascRegistry.fill(HIST("hITSnClustersNeg"), negTrack.itsNCls(), casc.mXi(), casc.pt());
    cascRegistry.fill(HIST("hITSnClustersBach"), bachTrack.itsNCls(), casc.mXi(), casc.pt());
    cascRegistry.fill(HIST("hTPCnCrossedRowsPos"), posTrack.tpcNClsCrossedRows(), casc.mXi(), casc.pt());
    cascRegistry.fill(HIST("hTPCnCrossedRowsNeg"), negTrack.tpcNClsCrossedRows(), casc.mXi(), casc.pt());
    cascRegistry.fill(HIST("hTPCnCrossedRowsBach"), bachTrack.tpcNClsCrossedRows(), casc.mXi(), casc.pt());
    cascRegistry.fill(HIST("hITSChi2Pos"), posTrack.itsChi2NCl());
    cascRegistry.fill(HIST("hITSChi2Neg"), negTrack.itsChi2NCl());
    cascRegistry.fill(HIST("hITSChi2Bach"), bachTrack.itsChi2NCl());
    cascRegistry.fill(HIST("hTPCChi2Pos"), posTrack.tpcChi2NCl());
    cascRegistry.fill(HIST("hTPCChi2Neg"), negTrack.tpcChi2NCl());
    cascRegistry.fill(HIST("hTPCChi2Bach"), bachTrack.tpcChi2NCl());

    cascRegistry.fill(HIST("hSelectionStatus"), 0);

    // [Phase 12b] Cascade pT lower cap — drop very-low-pT cascades early.
    if (casc.pt() < cascCfg.cMinCascPt)
      return 0;

    // [Phase 12b] Split V0-daughter and bachelor-track quality so the
    // bachelor can be tightened independently. Default values match the
    // V0-daughter values, so existing behaviour is preserved.
    if (posTrack.tpcNClsCrossedRows() < cascCfg.minTPCCrossedRows ||
        negTrack.tpcNClsCrossedRows() < cascCfg.minTPCCrossedRows)
      return 0;
    if (bachTrack.tpcNClsCrossedRows() < cascCfg.minBachTPCCrossedRows)
      return 0;
    cascRegistry.fill(HIST("hSelectionStatus"), 1);

    if (posTrack.itsNCls() < cascCfg.minITSClusters || negTrack.itsNCls() < cascCfg.minITSClusters)
      return 0;
    if (bachTrack.itsNCls() < cascCfg.minBachITSClusters)
      return 0;
    cascRegistry.fill(HIST("hSelectionStatus"), 2);

    if (posTrack.itsChi2NCl() > cascCfg.itsChi2 || negTrack.itsChi2NCl() > cascCfg.itsChi2)
      return 0;
    if (bachTrack.itsChi2NCl() > cascCfg.maxBachItsChi2)
      return 0;
    if (posTrack.tpcChi2NCl() > cascCfg.tpcChi2 || negTrack.tpcChi2NCl() > cascCfg.tpcChi2)
      return 0;
    if (bachTrack.tpcChi2NCl() > cascCfg.maxBachTpcChi2)
      return 0;
    cascRegistry.fill(HIST("hSelectionStatus"), 3);

    double pvx = collision.posX();
    double pvy = collision.posY();
    double pvz = collision.posZ();
    if (casc.v0radius() < cascCfg.v0setting_radius ||
        casc.cascradius() < cascCfg.cascadesetting_cascradius ||
        casc.v0cosPA(pvx, pvy, pvz) < cascCfg.v0setting_cospa ||
        casc.casccosPA(pvx, pvy, pvz) < cascCfg.cascadesetting_cospa ||
        casc.dcav0topv(pvx, pvy, pvz) < cascCfg.cascadesetting_mindcav0topv ||
        std::abs(casc.mLambda() - MassLambda0) > cascCfg.cascadesetting_v0masswindow)
      return 0;
    cascRegistry.fill(HIST("hSelectionStatus"), 4);

    // [Phase 8] Cascade-topology consistency: a real Ξ⁻/Ω⁻ decays first
    // (cascRadius), then its Λ-daughter flies further out and decays
    // (v0Radius). Reject candidates where this ordering is violated —
    // those are typically primary Λ paired with a stray bachelor whose
    // joint fit landed at a "cascade vertex" downstream of the V0
    // vertex (geometrically inconsistent with Ξ→Λπ kinematics).
    if (primCfg.cReqCascRadiusLessThanV0Radius && casc.cascradius() >= casc.v0radius())
      return 0;

    if (std::abs(posTrack.eta()) > cascCfg.etaTracks || std::abs(negTrack.eta()) > cascCfg.etaTracks || std::abs(bachTrack.eta()) > cascCfg.etaTracks)
      return 0;
    cascRegistry.fill(HIST("hSelectionStatus"), 5);

    if (std::abs(casc.eta()) > cascCfg.etaCascades)
      return 0;
    cascRegistry.fill(HIST("hSelectionStatus"), 6);

    if (casc.sign() < 0) {
      if (std::abs(posTrack.tpcNSigmaPr()) > cascCfg.tpcNsigmaProton)
        return 0;
      if (std::abs(negTrack.tpcNSigmaPi()) > cascCfg.tpcNsigmaPion)
        return 0;
    } else {
      if (std::abs(negTrack.tpcNSigmaPr()) > cascCfg.tpcNsigmaProton)
        return 0;
      if (std::abs(posTrack.tpcNSigmaPi()) > cascCfg.tpcNsigmaPion)
        return 0;
    }
    cascRegistry.fill(HIST("hSelectionStatus"), 7);

    int flag = 0;
    if (std::abs(bachTrack.tpcNSigmaPi()) < cascCfg.tpcNsigmaBachelor)
      flag = 1;
    if (std::abs(bachTrack.tpcNSigmaKa()) < cascCfg.tpcNsigmaBachelor && (!cascCfg.doCompetingMassCut || std::abs(o2::constants::physics::MassXiMinus - casc.mXi()) > cascCfg.competingMassWindow))
      flag = 3 - flag;

    switch (flag) {
      case 1:
        cascRegistry.fill(HIST("hSelectionStatus"), 8);
        if (casc.sign() < 0) {
          cascRegistry.fill(HIST("hMassXiMinus"), casc.mXi(), casc.pt(), casc.yXi());
        } else {
          cascRegistry.fill(HIST("hMassXiPlus"), casc.mXi(), casc.pt(), casc.yXi());
        }
        break;
      case 2:
        cascRegistry.fill(HIST("hSelectionStatus"), 8);
        if (casc.sign() < 0) {
          cascRegistry.fill(HIST("hMassXiMinus"), casc.mXi(), casc.pt(), casc.yXi());
          cascRegistry.fill(HIST("hMassOmegaMinus"), casc.mOmega(), casc.pt(), casc.yOmega());
        } else {
          cascRegistry.fill(HIST("hMassXiPlus"), casc.mXi(), casc.pt(), casc.yXi());
          cascRegistry.fill(HIST("hMassOmegaPlus"), casc.mOmega(), casc.pt(), casc.yOmega());
        }
        break;
      case 3:
        cascRegistry.fill(HIST("hSelectionStatus"), 8);
        if (casc.sign() < 0) {
          cascRegistry.fill(HIST("hMassOmegaMinus"), casc.mOmega(), casc.pt(), casc.yOmega());
        } else {
          cascRegistry.fill(HIST("hMassOmegaPlus"), casc.mOmega(), casc.pt(), casc.yOmega());
        }
        break;
    }

    return flag;
  }

  // [Phase 14] Compute the per-cascade cut bitmask UNCONDITIONALLY. Mirrors
  // every stage of processCandidate but doesn't short-circuit, so the bits
  // capture the cut state of EVERY cascade (even rejected ones). Cheap:
  // ~15 comparisons per cascade.
  template <typename TTracks, typename TCascade, typename TCollision>
  uint32_t computeCascadeCutBits(TCascade const& casc, TCollision const& collision)
  {
    uint32_t bits = 0;
    auto bachTrack = casc.template bachelor_as<TTracks>();
    auto posTrack  = casc.template posTrack_as<TTracks>();
    auto negTrack  = casc.template negTrack_as<TTracks>();

    if (posTrack.tpcNClsCrossedRows() >= cascCfg.minTPCCrossedRows &&
        negTrack.tpcNClsCrossedRows() >= cascCfg.minTPCCrossedRows)
      bits |= (1u << kCascCutTpcRowsV0Dau);
    if (bachTrack.tpcNClsCrossedRows() >= cascCfg.minBachTPCCrossedRows)
      bits |= (1u << kCascCutTpcRowsBach);

    if (posTrack.itsNCls() >= cascCfg.minITSClusters &&
        negTrack.itsNCls() >= cascCfg.minITSClusters)
      bits |= (1u << kCascCutItsClsV0Dau);
    if (bachTrack.itsNCls() >= cascCfg.minBachITSClusters)
      bits |= (1u << kCascCutItsClsBach);

    if (posTrack.itsChi2NCl() <= cascCfg.itsChi2 &&
        negTrack.itsChi2NCl() <= cascCfg.itsChi2 &&
        bachTrack.itsChi2NCl() <= cascCfg.maxBachItsChi2)
      bits |= (1u << kCascCutItsChi2);

    if (posTrack.tpcChi2NCl() <= cascCfg.tpcChi2 &&
        negTrack.tpcChi2NCl() <= cascCfg.tpcChi2 &&
        bachTrack.tpcChi2NCl() <= cascCfg.maxBachTpcChi2)
      bits |= (1u << kCascCutTpcChi2);

    if (casc.pt() >= cascCfg.cMinCascPt) bits |= (1u << kCascCutCascPt);

    double pvx = collision.posX(), pvy = collision.posY(), pvz = collision.posZ();
    bool topoOk = (casc.v0radius()         >= cascCfg.v0setting_radius           &&
                   casc.cascradius()       >= cascCfg.cascadesetting_cascradius  &&
                   casc.v0cosPA(pvx,pvy,pvz)   >= cascCfg.v0setting_cospa            &&
                   casc.casccosPA(pvx,pvy,pvz) >= cascCfg.cascadesetting_cospa       &&
                   casc.dcav0topv(pvx,pvy,pvz) >= cascCfg.cascadesetting_mindcav0topv &&
                   std::abs(casc.mLambda() - MassLambda0) <= cascCfg.cascadesetting_v0masswindow);
    if (topoOk) bits |= (1u << kCascCutTopology);

    if (casc.cascradius() < casc.v0radius())
      bits |= (1u << kCascCutRadiusOrder);

    if (std::abs(posTrack.eta()) <= cascCfg.etaTracks &&
        std::abs(negTrack.eta()) <= cascCfg.etaTracks &&
        std::abs(bachTrack.eta()) <= cascCfg.etaTracks)
      bits |= (1u << kCascCutTrackEta);

    if (std::abs(casc.eta()) <= cascCfg.etaCascades)
      bits |= (1u << kCascCutCascEta);

    if (casc.sign() < 0) {
      if (std::abs(posTrack.tpcNSigmaPr()) <= cascCfg.tpcNsigmaProton)
        bits |= (1u << kCascCutTpcNSigPr);
      if (std::abs(negTrack.tpcNSigmaPi()) <= cascCfg.tpcNsigmaPion)
        bits |= (1u << kCascCutTpcNSigPi);
    } else {
      if (std::abs(negTrack.tpcNSigmaPr()) <= cascCfg.tpcNsigmaProton)
        bits |= (1u << kCascCutTpcNSigPr);
      if (std::abs(posTrack.tpcNSigmaPi()) <= cascCfg.tpcNsigmaPion)
        bits |= (1u << kCascCutTpcNSigPi);
    }

    if (std::abs(bachTrack.tpcNSigmaPi()) < cascCfg.tpcNsigmaBachelor)
      bits |= (1u << kCascCutBachPidXi);
    if (std::abs(bachTrack.tpcNSigmaKa()) < cascCfg.tpcNsigmaBachelor)
      bits |= (1u << kCascCutBachPidOm);

    if (!cascCfg.doCompetingMassCut ||
        std::abs(o2::constants::physics::MassXiMinus - casc.mXi()) > cascCfg.competingMassWindow)
      bits |= (1u << kCascCutCompetingMass);

    return bits;
  }

  // [Phase 4] Helper for the per-collision cascade-flag loop. Emits one row
  // into cascflags per cascade in this collision (preserving the joinability
  // invariant with aod::CascDataExt). When eventOk=false, all flags are 0.
  // When eventOk=true and applyMcMatch=true (MC reco path), additionally
  // calls fillMatchedHistos to populate truerec/* histograms.
  // [Phase 4 fix] TTracks added so the call site forwards the bound
  // track type to processCandidate / fillMatchedHistos.
  template <bool ApplyMcMatch, typename TTracks, typename TCascades, typename TCollision>
  void cascadeFlagLoop(TCascades const& Cascades, TCollision const& collision, bool eventOk,
                       std::unordered_set<int64_t> const& itsTrackedCascIds)
  {
    cascRegistry.fill(HIST("hEventSel"), eventOk ? 1 : 2);
    for (auto const& casc : Cascades) {
      cascRegistry.fill(HIST("hEventSel"), 0);
      if (!eventOk) {
        cascflags(0, false, false, 0u);
        continue;
      }
      int flag = processCandidate<TTracks>(casc, collision);
      // [Phase 8] MC-truth purity flag.
      bool isTrueCasc = false;
      if constexpr (ApplyMcMatch) {
        if (flag != 0 && casc.has_mcParticle()) {
          auto gen = casc.mcParticle();
          int absPdg = std::abs(gen.pdgCode());
          if ((absPdg == lcorr_const::kXiMinusPdg || absPdg == lcorr_const::kOmegaMinusPdg) &&
              gen.isPhysicalPrimary()) {
            isTrueCasc = true;
          }
        }
      }
      // [Phase 9] ITS-tracking flag.
      bool isItsTracked = (itsTrackedCascIds.find(casc.cascadeId()) != itsTrackedCascIds.end());

      // [Phase 14] Per-cascade cut bitmask. Computed unconditionally — if
      // the cascade was rejected by processCandidate (flag==0) we still
      // record which stages it passed/failed. Bits 0-15 from
      // computeCascadeCutBits; bits 16-17 set here from the truth /
      // ITS-tracking lookups above.
      uint32_t cascCutBits = computeCascadeCutBits<TTracks>(casc, collision);
      if (isTrueCasc)    cascCutBits |= (1u << kCascCutMcTrueXiOmega);
      if (isItsTracked)  cascCutBits |= (1u << kCascCutItsTracked);

      cascflags(flag, isTrueCasc, isItsTracked, cascCutBits);
      if constexpr (ApplyMcMatch) {
        fillMatchedHistos<TTracks>(casc, flag, collision);
      }
    }
  }

  using CollisionsRun3 = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFV0As, aod::PVMults>;
  using CollisionsRun2 = soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms, aod::PVMults>;
  // [Phase 4 fix] aod::pidTPCKa added so the merged process function can
  // bind ONE track type and have it serve both V0 daughter cuts (Pi/Pr) and
  // cascade bachelor cuts (Ka). Without this, the merged process function
  // had to declare two track types and the framework's "last one wins" rule
  // mis-routed V0 index resolution to the wrong table.
  using Tracks = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::TrackCompColls>;
  using TracksRun2 = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCPr>;
  using TracksMC = soa::Join<Tracks, aod::McTrackLabels>;
  using TracksMCRun2 = soa::Join<TracksRun2, aod::McTrackLabels>;
  using McV0Tracks = soa::Join<aod::V0Datas, aod::McV0Labels>;

  // [Phase 4] processDataRun3 ABSORBS the former CascadeSelector::processRecData.
  // Single event-selection call drives both the V0-side LambdaCollections row
  // emission (inside fillLambdaRecoTables) and the Cascade-side flag emission
  // (inside cascadeFlagLoop). Impossible for the two halves to disagree.
  // [Phase 4 fix] Drop the explicit FullTracksExtIUWithPID arg: with one
  // bound track type (Tracks) the framework no longer mis-resolves V0 vs
  // cascade index targets. Cascade helpers see the same Tracks join via
  // the TTracks template parameter.
  void processDataRun3(CollisionsRun3::iterator const& collision,
                       aod::V0Datas const& V0s, Tracks const& tracks,
                       aod::CascDataExt const& Cascades,
                       aod::BCsWithTimestamps const&)
  {
    bool eventOk = selCollision<kRun3>(collision);
    // [Phase 9] Empty ITS-tracked set — variant without strangeness tracking.
    // Use processDataRun3WithItsTracking when the AOD actually contains
    // O2tracasccoll (aod::AssignedTrackedCascades). Many derived AODs drop it.
    static const std::unordered_set<int64_t> emptySet;
    cascadeFlagLoop<false, Tracks>(Cascades, collision, eventOk, emptySet);
    if (!eventOk) {
      return; // Lambda-side: don't emit a LambdaCollections row.
    }
    fillLambdaRecoTables<kRun3, kData>(collision, V0s, tracks);
  }

  PROCESS_SWITCH(LambdaCascadeProducer, processDataRun3, "Process for Run3 DATA (Lambda + Cascade, no ITS strangeness tracking)", true);

  // [Phase 9] ITS-tracking-aware data variant. Same as processDataRun3 but
  // also consumes aod::AssignedTrackedCascades; turn this on (and turn
  // processDataRun3 off) only when the AOD actually contains the
  // O2tracasccoll tree. Otherwise the AOD reader will fail at startup
  // with "Couldn't get TTree O2tracasccoll".
  void processDataRun3WithItsTracking(CollisionsRun3::iterator const& collision,
                                      aod::V0Datas const& V0s, Tracks const& tracks,
                                      aod::CascDataExt const& Cascades,
                                      aod::AssignedTrackedCascades const& trackedCascades,
                                      aod::BCsWithTimestamps const&)
  {
    bool eventOk = selCollision<kRun3>(collision);
    std::unordered_set<int64_t> itsTrackedCascIds;
    for (auto const& tc : trackedCascades) itsTrackedCascIds.insert(tc.cascadeId());
    cascadeFlagLoop<false, Tracks>(Cascades, collision, eventOk, itsTrackedCascIds);
    if (!eventOk) {
      return;
    }
    fillLambdaRecoTables<kRun3, kData>(collision, V0s, tracks);
  }

  PROCESS_SWITCH(LambdaCascadeProducer, processDataRun3WithItsTracking, "Run3 DATA + ITS strangeness tracking (requires O2tracasccoll in AOD)", false);

  void processDataRun2(CollisionsRun2::iterator const& collision, aod::V0Datas const& V0s, TracksRun2 const& tracks)
  {
    fillLambdaRecoTables<kRun2, kData>(collision, V0s, tracks);
  }

  PROCESS_SWITCH(LambdaCascadeProducer, processDataRun2, "Process for Run2 DATA (Lambda only — no cascade in Run2)", false);

  // [Phase 4] processMCRecoRun3 ABSORBS the former CascadeSelector::processRecMC.
  // Same pattern as processDataRun3 plus truth-matching on the cascade side
  // (fillMatchedHistos populates truerec/* histograms).
  void processMCRecoRun3(soa::Join<CollisionsRun3, aod::McCollisionLabels>::iterator const& collision, aod::McCollisions const&,
                         McV0Tracks const& V0s, TracksMC const& tracks, aod::McParticles const&,
                         LabeledCascades const& Cascades,
                         aod::BCsWithTimestamps const&)
  {
    bool eventOk = selCollision<kRun3>(collision);
    // [Phase 9] No ITS-tracking variant. Toggle to processMCRecoRun3WithItsTracking
    // when the AOD has O2tracasccoll.
    static const std::unordered_set<int64_t> emptySet;
    cascadeFlagLoop<true, TracksMC>(Cascades, collision, eventOk, emptySet);
    if (!eventOk) {
      return;
    }
    fillLambdaRecoTables<kRun3, kMC>(collision, V0s, tracks);
  }

  PROCESS_SWITCH(LambdaCascadeProducer, processMCRecoRun3, "Run3 MC reco (Lambda + Cascade with truth match, no ITS strangeness tracking)", false);

  // [Phase 9] ITS-tracking-aware MC reco variant.
  void processMCRecoRun3WithItsTracking(soa::Join<CollisionsRun3, aod::McCollisionLabels>::iterator const& collision, aod::McCollisions const&,
                                        McV0Tracks const& V0s, TracksMC const& tracks, aod::McParticles const&,
                                        LabeledCascades const& Cascades,
                                        aod::AssignedTrackedCascades const& trackedCascades,
                                        aod::BCsWithTimestamps const&)
  {
    bool eventOk = selCollision<kRun3>(collision);
    std::unordered_set<int64_t> itsTrackedCascIds;
    for (auto const& tc : trackedCascades) itsTrackedCascIds.insert(tc.cascadeId());
    cascadeFlagLoop<true, TracksMC>(Cascades, collision, eventOk, itsTrackedCascIds);
    if (!eventOk) {
      return;
    }
    fillLambdaRecoTables<kRun3, kMC>(collision, V0s, tracks);
  }

  PROCESS_SWITCH(LambdaCascadeProducer, processMCRecoRun3WithItsTracking, "Run3 MC reco + ITS strangeness tracking (requires O2tracasccoll in AOD)", false);

  void processMCRecoRun2(soa::Join<CollisionsRun2, aod::McCollisionLabels>::iterator const& collision, aod::McCollisions const&,
                         McV0Tracks const& V0s, TracksMCRun2 const& tracks, aod::McParticles const&)
  {
    // check collision
    if (!selCollision<kRun2>(collision)) {
      return;
    }
    fillLambdaRecoTables<kRun2, kMC>(collision, V0s, tracks);
  }

  PROCESS_SWITCH(LambdaCascadeProducer, processMCRecoRun2, "Process for Run2 MC reco (Lambda only)", false);

  void processMCRun3(aod::McCollisions::iterator const& mcCollision,
                     soa::SmallGroups<soa::Join<CollisionsRun3, aod::McCollisionLabels>> const& collisions,
                     McV0Tracks const& V0s, TracksMC const& tracks,
                     aod::McParticles const& mcParticles)
  {
    analyzeMcRecoGen<kRun3, kMC>(mcCollision, collisions, V0s, tracks, mcParticles);
  }

  PROCESS_SWITCH(LambdaCascadeProducer, processMCRun3, "Process for Run3 MC RecoGen", false);

  void processMCRun2(aod::McCollisions::iterator const& mcCollision,
                     soa::SmallGroups<soa::Join<CollisionsRun2, aod::McCollisionLabels>> const& collisions,
                     McV0Tracks const& V0s, TracksMCRun2 const& tracks,
                     aod::McParticles const& mcParticles)
  {
    analyzeMcRecoGen<kRun2, kMC>(mcCollision, collisions, V0s, tracks, mcParticles);
  }

  PROCESS_SWITCH(LambdaCascadeProducer, processMCRun2, "Process for Run2 MC RecoGen", false);

  // Gen-only: fills LambdaMcGenCollisions + LambdaMcGenTracks without touching reco tables.
  // Use alongside processMCRecoRun3 (which fills reco tables only).
  void processMCGenOnlyRun3(aod::McCollisions::iterator const& mcCollision,
                            soa::SmallGroups<soa::Join<CollisionsRun3, aod::McCollisionLabels>> const& collisions,
                            aod::McParticles const& mcParticles)
  {
    // Try to obtain centrality from a matched reco collision
    cent = 0.f;
    mult = 0.f;
    if (collisions.size() >= 1) {
      auto firstReco = collisions.begin();
      if (firstReco.has_mcCollision() &&
          firstReco.mcCollisionId() == mcCollision.globalIndex()) {
        selCollision<kRun3>(firstReco); // sets cent and mult
      }
    }
    fillLambdaMcGenTables<kRun3>(mcCollision, mcParticles);
  }

  PROCESS_SWITCH(LambdaCascadeProducer, processMCGenOnlyRun3, "Gen-only Run3 (no reco tables)", false);

  // ===========================================================================
  // [Phase 4] Cascade-side MC truth-only path (port of former
  // CascadeSelector::processGenMC). Self-contained MC-truth event selection
  // (INELgtN on McParticles, MC z-vertex), distinct from the reco-level
  // selCollision because it operates on McCollisions/McParticles directly.
  // Used as a closure/efficiency check for the cascade MC truth side.
  //
  // [Phase 4 fix] DO NOT add `private:` access specifiers here. Any non-public
  // non-static data member breaks the struct's aggregate-ness, and the
  // framework's `brace_constructible_size` (which uses brace-init to count
  // members for StructToTuple dispatch) collapses, mis-routing dispatch to
  // the D == 0 chain. The original CascadeSelector kept all members in the
  // implicit public section for the same reason. The previously-declared
  // `mCascCounter` was unused and has been removed.
  // ===========================================================================

  // [Phase 4 fix] Use CollisionsRun3 (which carries CentFT0Ms etc.) instead of
  // the old CSEL MyCollisions Join — selCollision<kRun3> needs centFT0M().
  void processCascadeGenMC(aod::McCollision const& mcCollision,
                            soa::SmallGroups<soa::Join<aod::McCollisionLabels, CollisionsRun3>> const& collisions,
                            aod::McParticles const& mcParticles)
  {
    // MC-truth event-level selection (mirror of CSEL::processGenMC).
    if (cascCfg.INEL >= 0 && !pwglf::isINELgtNmc(mcParticles, cascCfg.INEL, pdgDB))
      return;
    if (std::abs(mcCollision.posZ()) > cascCfg.maxVertexZ)
      return;

    cascRegistry.fill(HIST("gen/hNevents"), 0);

    for (auto const& mcPart : mcParticles) {
      if (!mcPart.isPhysicalPrimary())
        continue;
      if (std::abs(mcPart.eta()) > cascCfg.etaCascades)
        continue;
      switch (mcPart.pdgCode()) {
        case lcorr_const::kXiMinusPdg:     cascRegistry.fill(HIST("gen/hXiMinus"),    mcPart.pt(), mcPart.y()); break;
        case -lcorr_const::kXiMinusPdg:    cascRegistry.fill(HIST("gen/hXiPlus"),     mcPart.pt(), mcPart.y()); break;
        case lcorr_const::kOmegaMinusPdg:  cascRegistry.fill(HIST("gen/hOmegaMinus"), mcPart.pt(), mcPart.y()); break;
        case -lcorr_const::kOmegaMinusPdg: cascRegistry.fill(HIST("gen/hOmegaPlus"),  mcPart.pt(), mcPart.y()); break;
      }
    }

    // Now check whether at least one matched reco collision passes the
    // SAME event selection used for the data path (selCollision via the
    // shared lcorr_evsel framework). This guarantees the gen-with-rec
    // sample uses the same definition as the data path.
    if (collisions.size() < 1) {
      return;
    }
    bool evSel = false;
    for (auto const& collision : collisions) {
      if (selCollision<kRun3>(collision)) {
        evSel = true;
        break;
      }
    }
    if (!evSel)
      return;

    cascRegistry.fill(HIST("genwithrec/hNevents"), 0);
    for (auto const& mcPart : mcParticles) {
      if (!mcPart.isPhysicalPrimary())
        continue;
      if (std::abs(mcPart.eta()) > cascCfg.etaCascades)
        continue;
      switch (mcPart.pdgCode()) {
        case lcorr_const::kXiMinusPdg:     cascRegistry.fill(HIST("genwithrec/hXiMinus"),    mcPart.pt(), mcPart.y()); break;
        case -lcorr_const::kXiMinusPdg:    cascRegistry.fill(HIST("genwithrec/hXiPlus"),     mcPart.pt(), mcPart.y()); break;
        case lcorr_const::kOmegaMinusPdg:  cascRegistry.fill(HIST("genwithrec/hOmegaMinus"), mcPart.pt(), mcPart.y()); break;
        case -lcorr_const::kOmegaMinusPdg: cascRegistry.fill(HIST("genwithrec/hOmegaPlus"),  mcPart.pt(), mcPart.y()); break;
      }
    }
  }
  PROCESS_SWITCH(LambdaCascadeProducer, processCascadeGenMC, "Process Cascade MC gen (truth-only closure)", false);
};

struct LambdaTracksExtProducer {

  Produces<aod::LambdaTracksExt> lambdaTrackExtTable;

  // Configurables
  Configurable<bool> cAcceptAllLambda{"cAcceptAllLambda", false, "Accept all Lambda"};
  Configurable<bool> cRejAllLambdaShaDau{"cRejAllLambdaShaDau", true, "Reject all Lambda sharing daughters"};
  Configurable<bool> cSelLambdaMassPdg{"cSelLambdaMassPdg", false, "Select Lambda closest to Pdg Mass"};
  Configurable<bool> cSelLambdaTScore{"cSelLambdaTScore", false, "Select Lambda based on t-score"};
  Configurable<float> cA{"cA", 0.6, "a * |lambdaMass - lambdaPdgMass|"};
  Configurable<float> cB{"cB", 0.6, "b * DcaPrPi"};
  Configurable<float> cC{"cC", 0.6, "c * Cos(theta_{PA})"};

  // Histogram Registry.
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    // Axis Specifications
    const AxisSpec axisMult(10, 0, 10);
    const AxisSpec axisMass(100, 1.06, 1.16, "Inv Mass (GeV/#it{c}^{2})");
    const AxisSpec axisCPA(100, 0.995, 1.0, "cos(#theta_{PA})");
    const AxisSpec axisDcaDau(75, 0., 1.5, "Daug DCA (#sigma)");
    const AxisSpec axisDEta(320, -1.6, 1.6, "#Delta#eta");
    const AxisSpec axisDPhi(640, -PIHalf, 3. * PIHalf, "#Delta#varphi");

    // Histograms Booking
    histos.add("h1i_totlambda_mult", "Multiplicity", kTH1I, {axisMult});
    histos.add("h1i_totantilambda_mult", "Multiplicity", kTH1I, {axisMult});
    histos.add("h1i_lambda_mult", "Multiplicity", kTH1I, {axisMult});
    histos.add("h1i_antilambda_mult", "Multiplicity", kTH1I, {axisMult});
    histos.add("h2d_n2_etaphi_LaP_LaM", "#rho_{2}^{SharePair}", kTH2D, {axisDEta, axisDPhi});
    histos.add("h2d_n2_etaphi_LaP_LaP", "#rho_{2}^{SharePair}", kTH2D, {axisDEta, axisDPhi});
    histos.add("h2d_n2_etaphi_LaM_LaM", "#rho_{2}^{SharePair}", kTH2D, {axisDEta, axisDPhi});

    // InvMass, DcaDau and CosPA
    histos.add("Reco/h1f_lambda_invmass", "M_{p#pi}", kTH1F, {axisMass});
    histos.add("Reco/h1f_lambda_cospa", "cos(#theta_{PA})", kTH1F, {axisCPA});
    histos.add("Reco/h1f_lambda_dcadau", "DCA_{p#pi} at V0 Decay Vertex", kTH1F, {axisDcaDau});
    histos.add("Reco/h1f_antilambda_invmass", "M_{p#pi}", kTH1F, {axisMass});
    histos.add("Reco/h1f_antilambda_cospa", "cos(#theta_{PA})", kTH1F, {axisCPA});
    histos.add("Reco/h1f_antilambda_dcadau", "DCA_{p#pi} at V0 Decay Vertex", kTH1F, {axisDcaDau});

    histos.addClone("Reco/", "SharingDau/");

    // [P11] Flag-preference order is implicit (first matching else-if wins).
    // Warn loudly if the user enabled more than one strategy at once.
    int nStrategies = (cAcceptAllLambda ? 1 : 0) + (cRejAllLambdaShaDau ? 1 : 0) +
                      (cSelLambdaMassPdg ? 1 : 0) + (cSelLambdaTScore ? 1 : 0);
    if (nStrategies == 0) {
      LOGF(warning,
           "[LambdaTracksExtProducer] No true-Lambda selection strategy enabled. "
           "trueLambdaFlag will always be false and no Lambdas will pass the "
           "downstream goodLambda partition.");
    } else if (nStrategies > 1) {
      LOGF(warning,
           "[LambdaTracksExtProducer] %d true-Lambda strategies enabled at once "
           "(cAcceptAllLambda=%d cRejAllLambdaShaDau=%d cSelLambdaMassPdg=%d "
           "cSelLambdaTScore=%d). Selection uses else-if precedence in this order; "
           "later strategies will be ignored. Pick exactly one.",
           nStrategies,
           static_cast<int>(cAcceptAllLambda), static_cast<int>(cRejAllLambdaShaDau),
           static_cast<int>(cSelLambdaMassPdg), static_cast<int>(cSelLambdaTScore));
    } else {
      LOGF(info,
           "[LambdaTracksExtProducer] true-Lambda strategy: %s",
           cAcceptAllLambda ? "cAcceptAllLambda"
                            : cRejAllLambdaShaDau ? "cRejAllLambdaShaDau"
                            : cSelLambdaMassPdg ? "cSelLambdaMassPdg"
                                                : "cSelLambdaTScore");
    }
  }

  template <ShareDauLambda sd, typename T>
  void fillHistos(T const& track)
  {
    static constexpr std::string_view SubDir[] = {"Reco/", "SharingDau/"};

    if (track.v0Type() == kLambda) {
      histos.fill(HIST(SubDir[sd]) + HIST("h1f_lambda_invmass"), track.mass());
      histos.fill(HIST(SubDir[sd]) + HIST("h1f_lambda_dcadau"), track.dcaDau());
      histos.fill(HIST(SubDir[sd]) + HIST("h1f_lambda_cospa"), track.cosPA());
    } else {
      histos.fill(HIST(SubDir[sd]) + HIST("h1f_antilambda_invmass"), track.mass());
      histos.fill(HIST(SubDir[sd]) + HIST("h1f_antilambda_dcadau"), track.dcaDau());
      histos.fill(HIST(SubDir[sd]) + HIST("h1f_antilambda_cospa"), track.cosPA());
    }
  }

  void process(aod::LambdaCollisions::iterator const&, aod::LambdaTracks const& tracks)
  {

    int nTotLambda = 0, nTotAntiLambda = 0, nSelLambda = 0, nSelAntiLambda = 0;

    for (auto const& lambda : tracks) {
      bool lambdaMinDeltaMassFlag = true, lambdaMinTScoreFlag = true;
      bool lambdaSharingDauFlag = false, trueLambdaFlag = false;
      std::vector<int64_t> vSharedDauLambdaIndex;
      float tLambda = 0., tTrack = 0.;

      if (lambda.v0Type() == kLambda) {
        ++nTotLambda;
      } else if (lambda.v0Type() == kAntiLambda) {
        ++nTotAntiLambda;
      }

      tLambda = (cA * std::abs(lambda.mass() - MassLambda0)) + (cB * lambda.dcaDau()) + (cC * std::abs(lambda.cosPA() - 1.));

      for (auto const& track : tracks) {
        // check lambda index (don't analyze same lambda track !!!)
        if (lambda.index() == track.index()) {
          continue;
        }

        // check if lambda shares daughters with any other track
        if (lambda.posTrackId() == track.posTrackId() || lambda.negTrackId() == track.negTrackId()) {
          vSharedDauLambdaIndex.push_back(track.index());
          lambdaSharingDauFlag = true;

          // [Phase 16k] Fill the per-pair sharing-daughter Δη-Δφ histos
          // ONCE per unordered pair. The outer/inner loops iterate every
          // ordered pair, so without this guard each unique (A,B) shared
          // pair would fill twice (with opposite Δη/Δφ sign on the second
          // pass). Keeping only track.index() > lambda.index() picks the
          // canonical ordering.
          if (track.index() > lambda.index()) {
            // Fill DEta-DPhi Histogram
            if ((lambda.v0Type() == kLambda && track.v0Type() == kAntiLambda) || (lambda.v0Type() == kAntiLambda && track.v0Type() == kLambda)) {
              histos.fill(HIST("h2d_n2_etaphi_LaP_LaM"), lambda.eta() - track.eta(), RecoDecay::constrainAngle((lambda.phi() - track.phi()), -PIHalf));
            } else if (lambda.v0Type() == kLambda && track.v0Type() == kLambda) {
              histos.fill(HIST("h2d_n2_etaphi_LaP_LaP"), lambda.eta() - track.eta(), RecoDecay::constrainAngle((lambda.phi() - track.phi()), -PIHalf));
            } else if (lambda.v0Type() == kAntiLambda && track.v0Type() == kAntiLambda) {
              histos.fill(HIST("h2d_n2_etaphi_LaM_LaM"), lambda.eta() - track.eta(), RecoDecay::constrainAngle((lambda.phi() - track.phi()), -PIHalf));
            }
          }

          // decision based on mass closest to PdgMass of Lambda
          if (std::abs(lambda.mass() - MassLambda0) > std::abs(track.mass() - MassLambda0)) {
            lambdaMinDeltaMassFlag = false;
          }

          // decisions based on t-score
          tTrack = (cA * std::abs(track.mass() - MassLambda0)) + (cB * track.dcaDau()) + (cC * std::abs(track.cosPA() - 1.));
          if (tLambda > tTrack) {
            lambdaMinTScoreFlag = false;
          }
        }
      }

      // fill QA histograms
      if (lambdaSharingDauFlag) {
        fillHistos<kLambdaShareDau>(lambda);
      } else {
        fillHistos<kUniqueLambda>(lambda);
      }

      if (cAcceptAllLambda) { // Accept all lambda
        trueLambdaFlag = true;
      } else if (cRejAllLambdaShaDau && !lambdaSharingDauFlag) { // Reject all lambda sharing daughter
        trueLambdaFlag = true;
      } else if (cSelLambdaMassPdg && lambdaMinDeltaMassFlag) { // Select lambda closest to pdg mass
        trueLambdaFlag = true;
      } else if (cSelLambdaTScore && lambdaMinTScoreFlag) { // Select lambda based on t-score
        trueLambdaFlag = true;
      }

      // Multiplicity of selected lambda
      if (trueLambdaFlag) {
        if (lambda.v0Type() == kLambda) {
          ++nSelLambda;
        } else if (lambda.v0Type() == kAntiLambda) {
          ++nSelAntiLambda;
        }
      }

      // fill LambdaTrackExt table
      lambdaTrackExtTable(lambdaSharingDauFlag, vSharedDauLambdaIndex, trueLambdaFlag);
    }

    // fill multiplicity histograms
    if (nTotLambda != 0) {
      histos.fill(HIST("h1i_totlambda_mult"), nTotLambda);
    }

    if (nTotAntiLambda != 0) {
      histos.fill(HIST("h1i_totantilambda_mult"), nTotAntiLambda);
    }

    if (nSelLambda != 0) {
      histos.fill(HIST("h1i_lambda_mult"), nSelLambda);
    }

    if (nSelAntiLambda != 0) {
      histos.fill(HIST("h1i_antilambda_mult"), nSelAntiLambda);
    }
  }
};


// =============================================================================
// [Phase 5] Dead structs removed:
//
//   - LambdaR2Correlation  (was struct at this position; commented out of the
//     workflow since Phase 1, ~340 lines of unused code)
//   - CascadeSelector_REMOVED_PHASE4 (#if 0 block from Phase 4, kept temporarily
//     for reference, ~575 lines)
//   - CascadeCorrelations  (was struct after Phase 4 CSEL block; commented out
//     of the workflow since Phase 1, ~450 lines)
//
// Total: roughly 1370 fewer lines compiled per build. All three structs are
// available in git history if anyone wants to revive them.
// =============================================================================




// [Phase 12d → reverted Phase 16b] TTree branch structs + connectors inlined
// back into the .cxx (ALICE submission policy: one source file). Plain POD,
// no O2 framework types — only ROOT TTree is used.
namespace lxicorr
{
struct CascBranches {
  float pt, rap, mass;
  int sign;
  float cascCosPA, v0CosPA;
  float cascRadius, v0Radius;
  float dcaV0Dau, dcaCascDau;
  float dcaV0ToPV, dcaPosToPV, dcaNegToPV, dcaBachToPV;
  float cent, pvZ;
  // MC truth-matching (filled only in MC mode; default = no match)
  int pdgCode{0}; // PDG of matched McParticle (0 = no match = background)
  bool isPhysPrim{false};
  // [Phase 14] Per-cascade cut bitmask + bachelor PID inputs for offline
  // re-cutting. cascCutBits semantics defined by enum CascCutBit below.
  unsigned int cascCutBits{0};
  float bachTpcNSigmaPi{0.f};
  float bachTpcNSigmaKa{0.f};
  int   bachItsNCls{0};
  int   bachTpcNClsCrossedRows{0};
  float mLambdaInside{0.f};   // V0 mass-hypothesis as seen by the cascade fit
};

// Holds branch data for both Xi and Omega — one raw pointer covers both.
// Gen-level branches — kinematics + truth only, no topology
struct GenBranches {
  float pt, rap, cent, pvZ;
  int pdgCode;
  bool isPhysPrim;
};

// [Phase 7] Λ TTree branches. Both reco-level (with topo + truth fields)
// and gen-level. The motherPdg field lets the user filter feed-down sources
// post-hoc (e.g., motherPdg==±3312 → Ξ-feeddown). isPhysPrim duplicates the
// primary/secondary tag from lambdatrack::V0PrmScd for convenience.
struct LambdaBranches {
  float pt, eta, rap, phi, mass;
  float cosPA, dcaDau;
  int   v0Type;     // 0=Λ, 1=Λ̄  (matches enum ParticleType)
  int   v0PrmScd;   // 0=primary, 1=secondary (only meaningful in MC)
  bool  trueLambdaFlag;
  float corrFact;
  float cent, pvZ;
  // MC-only
  int   pdgCode{0};
  bool  isPhysPrim{false};
  int   motherPdg{0};
  // [Phase 8] Topology snapshot — the inputs of a downstream
  // primary-fraction template fit. Per-V0, filled on every row.
  float dcaV0ToPV{0.f};
  float v0Radius{0.f};
  int   posItsNCls{0};
  int   negItsNCls{0};
  bool  passesPrimaryTopo{false};
  // [Phase 9] Per-daughter 7-bit ITS hit-map for geometric consistency
  // checks downstream (vetoes against impossible inner-layer hits when
  // the V0 vertex is outside that layer).
  unsigned int posItsClusterMap{0};
  unsigned int negItsClusterMap{0};
  // [Phase 10] Pseudo-proper decay length and per-daughter DCA-XY.
  float lProper{0.f};
  float posDcaXY{0.f};
  float negDcaXY{0.f};
  // [Phase 14] Per-V0 cut bitmask + raw cut inputs for offline re-cutting.
  unsigned int cutBits{0};
  float tpcNSigmaPosPr{0.f};   // proton hypothesis on positive daughter
  float tpcNSigmaNegPi{0.f};   // pion hypothesis on negative daughter
  float tpcNSigmaPosPi{0.f};   // pion hypothesis on positive daughter (anti-Λ side)
  float tpcNSigmaNegPr{0.f};   // proton hypothesis on negative daughter (anti-Λ side)
  float mK0Short{0.f};
  float qtArm{0.f};
  float alphaArm{0.f};
  float cTau{0.f};
};

struct LambdaGenBranches {
  float pt, eta, rap, phi;
  int   v0Type;
  int   v0PrmScd;
  float cent, pvZ;
};

// Single pointer covers all six branch sets — 1 StructToTuple slot
struct BranchPair {
  CascBranches xi{};
  CascBranches om{};
  GenBranches xiGen{};
  GenBranches omGen{};
  LambdaBranches lam{};
  LambdaGenBranches lamGen{};
};

inline void connectBranches(TTree* t, CascBranches* b)
{
  t->Branch("pt", &b->pt);
  t->Branch("rap", &b->rap);
  t->Branch("mass", &b->mass);
  t->Branch("sign", &b->sign);
  t->Branch("cascCosPA", &b->cascCosPA);
  t->Branch("v0CosPA", &b->v0CosPA);
  t->Branch("cascRadius", &b->cascRadius);
  t->Branch("v0Radius", &b->v0Radius);
  t->Branch("dcaV0Dau", &b->dcaV0Dau);
  t->Branch("dcaCascDau", &b->dcaCascDau);
  t->Branch("dcaV0ToPV", &b->dcaV0ToPV);
  t->Branch("dcaPosToPV", &b->dcaPosToPV);
  t->Branch("dcaNegToPV", &b->dcaNegToPV);
  t->Branch("dcaBachToPV", &b->dcaBachToPV);
  t->Branch("cent", &b->cent);
  t->Branch("pvZ", &b->pvZ);
  t->Branch("pdgCode", &b->pdgCode);
  t->Branch("isPhysPrim", &b->isPhysPrim);
  // [Phase 14] Per-cascade cut bitmask + bachelor PID + V0-inside mass.
  t->Branch("cascCutBits",            &b->cascCutBits);
  t->Branch("bachTpcNSigmaPi",        &b->bachTpcNSigmaPi);
  t->Branch("bachTpcNSigmaKa",        &b->bachTpcNSigmaKa);
  t->Branch("bachItsNCls",            &b->bachItsNCls);
  t->Branch("bachTpcNClsCrossedRows", &b->bachTpcNClsCrossedRows);
  t->Branch("mLambdaInside",          &b->mLambdaInside);
}

inline void connectGenBranches(TTree* t, GenBranches* b)
{
  t->Branch("pt", &b->pt);
  t->Branch("rap", &b->rap);
  t->Branch("pdgCode", &b->pdgCode);
  t->Branch("isPhysPrim", &b->isPhysPrim);
  t->Branch("cent", &b->cent);
  t->Branch("pvZ", &b->pvZ);
}

// [Phase 7] Λ-side connectors.
inline void connectLambdaBranches(TTree* t, LambdaBranches* b)
{
  t->Branch("pt", &b->pt);
  t->Branch("eta", &b->eta);
  t->Branch("rap", &b->rap);
  t->Branch("phi", &b->phi);
  t->Branch("mass", &b->mass);
  t->Branch("cosPA", &b->cosPA);
  t->Branch("dcaDau", &b->dcaDau);
  t->Branch("v0Type", &b->v0Type);
  t->Branch("v0PrmScd", &b->v0PrmScd);
  t->Branch("trueLambdaFlag", &b->trueLambdaFlag);
  t->Branch("corrFact", &b->corrFact);
  t->Branch("cent", &b->cent);
  t->Branch("pvZ", &b->pvZ);
  t->Branch("pdgCode", &b->pdgCode);
  t->Branch("isPhysPrim", &b->isPhysPrim);
  t->Branch("motherPdg", &b->motherPdg);
  // [Phase 8] Topology snapshot for primary-fraction template fits on data.
  t->Branch("dcaV0ToPV", &b->dcaV0ToPV);
  t->Branch("v0Radius", &b->v0Radius);
  t->Branch("posItsNCls", &b->posItsNCls);
  t->Branch("negItsNCls", &b->negItsNCls);
  t->Branch("passesPrimaryTopo", &b->passesPrimaryTopo);
  // [Phase 9] ITS hit-map (7 bits used). Stored as unsigned int for
  // ROOT-friendly typing on the read side.
  t->Branch("posItsClusterMap", &b->posItsClusterMap);
  t->Branch("negItsClusterMap", &b->negItsClusterMap);
  // [Phase 10] L_proper (cm) and per-daughter signed DCA-XY-to-PV (cm).
  t->Branch("lProper", &b->lProper);
  t->Branch("posDcaXY", &b->posDcaXY);
  t->Branch("negDcaXY", &b->negDcaXY);
  // [Phase 14] Per-V0 cut bitmask + raw cut inputs.
  t->Branch("cutBits",        &b->cutBits);
  t->Branch("tpcNSigmaPosPr", &b->tpcNSigmaPosPr);
  t->Branch("tpcNSigmaNegPi", &b->tpcNSigmaNegPi);
  t->Branch("tpcNSigmaPosPi", &b->tpcNSigmaPosPi);
  t->Branch("tpcNSigmaNegPr", &b->tpcNSigmaNegPr);
  t->Branch("mK0Short",       &b->mK0Short);
  t->Branch("qtArm",          &b->qtArm);
  t->Branch("alphaArm",       &b->alphaArm);
  t->Branch("cTau",           &b->cTau);
}

inline void connectLambdaGenBranches(TTree* t, LambdaGenBranches* b)
{
  t->Branch("pt", &b->pt);
  t->Branch("eta", &b->eta);
  t->Branch("rap", &b->rap);
  t->Branch("phi", &b->phi);
  t->Branch("v0Type", &b->v0Type);
  t->Branch("v0PrmScd", &b->v0PrmScd);
  t->Branch("cent", &b->cent);
  t->Branch("pvZ", &b->pvZ);
}
} // namespace lxicorr


struct LambdaXiCorrelation {

  // --- Configurables ---
  Configurable<float> maxY{"maxY", 0.5, "Max |y| for Lambda, Xi and Omega"};
  Configurable<bool> useEff{"useEff", false, "Apply Lambda efficiency correction"};
  // [Phase 6] Cascade efficiency weighting (mirror of useEff for Lambdas).
  // When true, the pair loops include a cascade-side weight via
  // getCascadeEfficiency<>(...). The default helper returns 1.0 — wire
  // CCDB loading into it when you have an efficiency map (TODO comment in
  // body).
  Configurable<bool> cUseCascEff{"cUseCascEff", false, "Apply cascade efficiency correction"};
  // Omega bachelor kaon PID cut (kept separate from Xi pion cut for physics correctness)
  Configurable<float> tpcNsigmaBachKaon{"tpcNsigmaBachKaon", 3.0f, "TPC NSigma bachelor kaon (Omega selection)"};
  // FT0M centrality bins — same variable-width defaults as LambdaR2Correlation
  ConfigurableAxis centAxis{"centAxis", {VARIABLE_WIDTH, 0.0f, 10.0f, 30.0f, 50.0f, 80.0f, 100.0f}, "FT0M centrality (%)"};
  Configurable<bool> saveCascTree{"saveCascTree", false, "Save TTree of cascade topological variables into AnalysisResults.root"};

  // [Phase 7] Use only Λ flagged kPrimary in the standalone Λ list. In MC
  // this filters truly-primary Λ via mcParticle.isPhysicalPrimary() (set by
  // LambdaCascadeProducer when cSelMCPSV0=true). On data, every Λ is
  // unconditionally tagged kPrimary, so this toggle is effectively a no-op
  // — combine with a tighter cMaxDcaV0ToPV in lambda-cascade-producer for
  // topological enrichment of the primary fraction.
  Configurable<bool> cUsePrimaryLambdasOnly{"cUsePrimaryLambdasOnly", true,
      "Restrict the trigger Λ list to v0PrmScd==kPrimary (MC-meaningful)"};

  // [Phase 7] Save a per-Λ TTree of the standalone V0 candidates. Useful for
  // offline cut tuning, primary-fraction template fits, and feed-down
  // post-tagging via the motherPdg branch.
  Configurable<bool> saveLambdaTree{"saveLambdaTree", false,
      "Save TTree of standalone-Λ candidates into AnalysisResults.root"};

  // [Phase 6] Auto-correlation veto policy (replaces previous bool
  // cVetoSharedDau). Three modes for systematic studies:
  //   0 = off (keep all pairs, including same-Λ self-pairs)
  //   1 = strict — drop pair if pos AND neg daughter tracks both match
  //                the cascade's V0 daughters. (Default; matches the old
  //                cVetoSharedDau=true behaviour.)
  //   2 = loose  — drop pair if EITHER daughter matches. Catches the
  //                rare case where the standalone V0 was fit to one
  //                cascade-V0 daughter plus a different second track.
  // The legacy `cVetoSharedDau` Configurable is still parsed from JSON
  // for backward compatibility but ignored at runtime — set cVetoMode
  // instead.
  Configurable<int> cVetoMode{"cVetoMode", 1,
      "Auto-correlation veto: 0=off, 1=strict (both daughters shared), 2=loose (any daughter)"};
  Configurable<bool> cVetoSharedDau{"cVetoSharedDau", true,
      "[DEPRECATED Phase 6] use cVetoMode instead. Value is ignored at runtime."};

  // [Phase 6] Histogram-axis bundle. ConfigurableGroup makes this one
  // StructToTuple slot while still exposing every axis as JSON-tunable
  // (cDphiAxis.values, etc.). Lets users adjust correlation-function
  // resolution without recompiling — change a JSON value, rerun.
  struct : ConfigurableGroup {
    ConfigurableAxis cDphiAxis{"cDphiAxis", {72, -PIHalf, 3 * PIHalf}, "Δφ"};
    ConfigurableAxis cDyAxis{"cDyAxis", {40, -2.0f, 2.0f}, "Δy"};
    ConfigurableAxis cPtAxisPair{"cPtAxisPair", {100, 0, 10}, "pT axis on single+pair histos"};
    ConfigurableAxis cMassLamAxis{"cMassLamAxis", {100, 1.09, 1.14}, "Λ inv-mass"};
    ConfigurableAxis cMassXiAxis{"cMassXiAxis", {100, 1.28, 1.36}, "Ξ inv-mass"};
    ConfigurableAxis cMassOmAxis{"cMassOmAxis", {100, 1.62, 1.72}, "Ω inv-mass"};
    ConfigurableAxis cRapAxis{"cRapAxis", {100, -1.0f, 1.0f}, "rapidity"};
    // [Phase 11] φ axis on single-particle (y, φ) histograms — needed
    // to build ρ₁(y, φ) and convolve into ρ₁⊗ρ₁ in (Δy, Δφ).
    ConfigurableAxis cPhiAxis{"cPhiAxis", {72, 0.f, 2.0f * static_cast<float>(M_PI)}, "φ (rad)"};
  } histAxes;

  // --- Outputs ---
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  // Direct OutputObj members → framework registers and writes these into
  // AnalysisResults.root under the "CascadeTrees" folder automatically.
  // setObject() is called in init() when saveCascTree is true.
  OutputObj<TTree> treeXi{"XiCandidates", OutputObjHandlingPolicy::AnalysisObject};
  OutputObj<TTree> treeOmega{"OmegaCandidates", OutputObjHandlingPolicy::AnalysisObject};
  OutputObj<TTree> treeXiGen{"XiCandidatesGen", OutputObjHandlingPolicy::AnalysisObject};
  OutputObj<TTree> treeOmegaGen{"OmegaCandidatesGen", OutputObjHandlingPolicy::AnalysisObject};
  // [Phase 7] Λ trees (reco + gen).
  OutputObj<TTree> treeLambda{"LambdaCandidates", OutputObjHandlingPolicy::AnalysisObject};
  OutputObj<TTree> treeLambdaGen{"LambdaCandidatesGen", OutputObjHandlingPolicy::AnalysisObject};
  // Single raw pointer covers all six branch sets — 1 StructToTuple slot
  lxicorr::BranchPair* bp{nullptr};

  // [Phase 16f] Last-seen collision index for the Λ-side singles + tree fill.
  // analyzeSinglesLambda() is invoked from BOTH processMCRecoXi() and
  // processMCRecoOmega() (and the equivalent data process functions). When
  // multiple of those flags are on simultaneously — the default for full
  // Λ-Ξ + Λ-Ω runs — every Λ would otherwise be filled into Singles/* and
  // LambdaCandidates twice (or more). The guard below makes the per-event
  // Λ fill idempotent: only the first call per collisionId runs the body,
  // subsequent calls from the same event return early.
  // -1 marks "no event seen yet" (legal collisionId is >= 0).
  int64_t mLastLamSinglesCollIdx{-1};

  // [Phase 16f → split Phase 16n] Idempotence guards for Event/hEventCount.
  // Previously a single mLastEventCountCollIdx — but reco process functions
  // pass LambdaCollisions::globalIndex() and gen ones pass
  // LambdaMcGenCollisions::globalIndex(), which live in different ID spaces.
  // A single guard double-counted events when both reco and gen process
  // functions were active. Split into two so each "side" of the workflow
  // counts unique events independently. Reco counter is the canonical one
  // and is the one filled into the histogram.
  int64_t mLastEventCountRecoCollIdx{-1};
  int64_t mLastEventCountGenCollIdx{-1};

  // [Phase 16i] Gen-Λ singles guard. processMCGenXi and processMCGenOmega
  // both iterate the same LambdaMcGenTracks table and fill the same
  // McGen/Singles/Lambda/* histograms. When both flags are on (typical for
  // a full Λ-Ξ + Λ-Ω closure run) every gen Λ would be counted twice.
  // The guard makes the singles fills idempotent per gen event.
  int64_t mLastGenLamSinglesCollIdx{-1};

  // Helper: fill the per-event counter exactly once per reco collision.
  // [Phase 16n] Reco-side guard. The histogram itself is filled here so
  // hEventCount reports unique reco collisions (the canonical event count
  // used everywhere downstream).
  void fillEventCountOnce(int64_t collIdx)
  {
    if (collIdx == mLastEventCountRecoCollIdx) {
      return;
    }
    mLastEventCountRecoCollIdx = collIdx;
    histos.fill(HIST("Event/hEventCount"), 0.5);
  }

  // [Phase 16n] Gen-side guard. Only tracks "have we visited this gen
  // collision yet" for downstream fills that need it (currently nothing
  // calls this; reserved). Does NOT fill hEventCount — that is the reco
  // counter's job.
  bool seenGenCollOnce(int64_t mcgenCollIdx)
  {
    if (mcgenCollIdx == mLastEventCountGenCollIdx) {
      return false;
    }
    mLastEventCountGenCollIdx = mcgenCollIdx;
    return true;
  }

  // [Phase 16g] Count cascades eligible for a given species in THIS event.
  // The slice from CascDataExt/LabeledCascades contains both Ξ and Ω
  // candidates; without species filtering, hLamXi and hLamOm both report
  // total cascades (Ξ+Ω) → identical numbers. flag semantics defined as
  // named constants below.
  // Applies the same |y|<maxY cut as the corresponding pair loop so the
  // diagnostic matches actual pair-formation potential.
  // [Phase 16j] Aliases of namespace-scope constants so existing call sites
  // (analyzePairs, analyzeXiXiPairs, …) keep their short names while the
  // single source of truth lives in lcorr_const.
  static constexpr int kFlagRejected     = lcorr_const::kFlagRejected;
  static constexpr int kFlagXiOnly       = lcorr_const::kFlagXiOnly;
  static constexpr int kFlagXiAndOmega   = lcorr_const::kFlagXiAndOmega;
  static constexpr int kFlagOmegaOnly    = lcorr_const::kFlagOmegaOnly;

  template <bool IsOmega, typename TCascades, typename TFlagsRow>
  int countSpeciesEligible(TCascades const& cascades, TFlagsRow const& flagsStart)
  {
    int n = 0;
    for (const auto& c : cascades) {
      auto fr = flagsStart + c.globalIndex();
      int f = fr.isSelected();
      if (f == kFlagRejected) {
        continue;
      }
      if constexpr (IsOmega) {
        if ((f == kFlagXiAndOmega || f == kFlagOmegaOnly) && std::abs(c.yOmega()) <= maxY) {
          ++n;
        }
      } else {
        if ((f == kFlagXiOnly || f == kFlagXiAndOmega) && std::abs(c.yXi()) <= maxY) {
          ++n;
        }
      }
    }
    return n;
  }

  // --- Data Slicing Definitions ---
  using GoodLambdas = soa::Join<aod::LambdaTracks, aod::LambdaTracksExt>;
  Partition<GoodLambdas> goodLambda = aod::lambdatrackext::trueLambdaFlag == true;
  // [Phase 7] Primary-only Λ partition (MC-meaningful; on data v0PrmScd is
  // always kPrimary so the partition coincides with goodLambda).
  Partition<GoodLambdas> goodPrimaryLambda =
      aod::lambdatrackext::trueLambdaFlag == true &&
      aod::lambdatrack::v0PrmScd == (int8_t)kPrimary;

  // [Phase 8] Topology-only primary-Λ partition. Combine with the truth
  // partition to get the strictest trigger list: AND of (trueLambdaFlag,
  // v0PrmScd==kPrimary, passesPrimaryTopo). On data the v0PrmScd term is
  // identically true so the partition reduces to (trueLambdaFlag &&
  // passesPrimaryTopo) — i.e., topology alone, which is the only handle
  // available without truth.
  Partition<GoodLambdas> goodPureLambda =
      aod::lambdatrackext::trueLambdaFlag == true &&
      aod::lambdatrack::v0PrmScd == (int8_t)kPrimary &&
      aod::lambdatrack::passesPrimaryTopo == true;

  // [Phase 8] Pick which "primary Λ" partition the trigger uses.
  //   true  → goodPureLambda (truth AND topology; recommended for data)
  //   false → goodPrimaryLambda (truth-only; use for MC closure tests)
  Configurable<bool> cPrimaryRequireTopo{"cPrimaryRequireTopo", true,
      "Require passesPrimaryTopo on the trigger Λ list (data-friendly tightening)"};

  // [Phase 8] MC-only filter: reject cascades whose CascadeFlags::IsTrueCascade
  // is false (i.e., the cascade was not truth-matched to a physical-primary
  // Ξ⁻/Ω⁻). On data this flag is always false (no truth available), so the
  // filter is automatically disabled in the data process functions — it only
  // gates analyzeSinglesXi/Omega and analyzePairs/OmegaPairs when called via
  // <IsMC=true>. Use this for closure/purity studies; turn off when measuring
  // the residual combinatorial contamination of your topological selection.
  Configurable<bool> cRequireTrueCascade{"cRequireTrueCascade", true,
      "MC: drop cascade candidates not flagged IsTrueCascade==true"};

  // [Phase 9] ITS strangeness-tracking gating policy. Lets you A/B-compare
  // what the IsItsTracked flag adds without recompiling.
  //
  //   0 = OFF      — ignore IsItsTracked entirely (current behaviour).
  //   1 = REQUIRE  — accept ONLY cascades flagged IsItsTracked. Highest
  //                  purity on data (ITS confirmed the Ξ flight) but
  //                  trades a lot of statistics: ITS strangeness-tracking
  //                  efficiency is typically 5–20% of the topological
  //                  cascade yield in pp, depending on pT and η.
  //   2 = RESCUE   — keep cascades that fail purity gates (in particular
  //                  cRequireTrueCascade on MC) IF IsItsTracked is true.
  //                  Useful for data-side high-purity systematic where
  //                  the truth flag is unavailable but ITS confirms.
  //   3 = QA-ONLY  — no gating; only fill QA histos so you can plot
  //                  the IsItsTracked fraction vs cascade pT/centrality.
  //
  // QA histograms are filled regardless of the mode so the comparison is
  // always available downstream.
  Configurable<int> cItsTrackMode{"cItsTrackMode", 0,
      "ITS strangeness-tracking gate: 0=off, 1=require, 2=rescue, 3=qa-only"};

  // [Phase 12b] Trigger-pT-differential pair histograms (off by default to
  // avoid bloating output for the centrality-only analyses).
  Configurable<bool> cFillPtDifferentialPairs{"cFillPtDifferentialPairs", false,
      "Fill PairsPt/<sign>/hPtDeltaPhiDeltaY (3D pT_Λ × Δφ × Δy)"};

  // [Phase 13a] Pair-type selectors. Each toggle gates BOTH the histogram
  // booking and the fill, so disabled pair types cost zero memory and
  // zero CPU. Defaults: the standard Λ-Ξ / Λ-Ω physics + Λ-Λ on; the
  // same-species cascade pairs and the Ξ-Ω cross-species pairs are off
  // until you ask for them.
  struct : ConfigurableGroup {
    Configurable<bool> cFillLamXi{"cFillLamXi", true, "Λ–Ξ pairs (4 sign combos)"};
    Configurable<bool> cFillLamOm{"cFillLamOm", true, "Λ–Ω pairs (4 sign combos)"};
    Configurable<bool> cFillLamLam{"cFillLamLam", true, "Λ–Λ same-event pairs (3 sign combos)"};
    Configurable<bool> cFillXiXi{"cFillXiXi", false, "Ξ–Ξ same-event pairs (3 sign combos)"};
    Configurable<bool> cFillOmOm{"cFillOmOm", false, "Ω–Ω same-event pairs (3 sign combos)"};
    Configurable<bool> cFillXiOm{"cFillXiOm", false, "Ξ–Ω cross-species pairs (4 sign combos)"};
  } pairCfg;

  // [Phase 13b] Per-event yields. One toggle, one ConfigurableAxis bundle.
  // Each species gets 5 histograms under Yields/<species>/:
  //   hNPerEvent      — multiplicity distribution (TH1F)         → ⟨N_β⟩ = Mean
  //   hMeanPtPerEvent — per-event ⟨pT⟩ distribution (TH1F)       → ⟨pT⟩ = Mean
  //   hMeanNvsCent    — TProfile ⟨N_β⟩ vs centrality
  //   hMeanPtVsCent   — TProfile ⟨pT⟩ vs centrality
  //   hNvsPt2D        — (N, ⟨pT⟩) per event (TH2F)               (correlation)
  // Filled exclusively by the dedicated processYields function so that the
  // count is correct regardless of which pair-process functions are on.
  struct : ConfigurableGroup {
    Configurable<bool> cFillEventYields{"cFillEventYields", true, "Per-event multiplicity + ⟨pT⟩ histograms for each species"};
    ConfigurableAxis cYieldNAxis{"cYieldNAxis", {50, 0., 50.}, "N per event"};
    ConfigurableAxis cYieldPtAxis{"cYieldPtAxis", {100, 0., 10.}, "⟨pT⟩ per event (GeV/c)"};
  } yieldCfg;

  Preslice<GoodLambdas> lambdasPerCollision = aod::lambdatrack::lambdaCollisionId;
  Preslice<aod::CascDataExt> cascadesPerCollision = aod::cascdata::collisionId;
  Preslice<LabeledCascades> labeledCascPerCollision = aod::cascdata::collisionId;

  // [Phase 12a] Single source of truth for the trigger-Λ partition pick.
  // Replaces a 4-line ternary that was duplicated verbatim in each of the
  // four process functions (processXi/Omega/MCRecoXi/MCRecoOmega). The
  // partitions all slice to the same SOA type, so a single auto-returning
  // helper is sufficient.
  template <typename TColl>
  auto pickLambdaPartition(TColl const& lambdacoll)
  {
    if (!cUsePrimaryLambdasOnly)
      return goodLambda->sliceBy(lambdasPerCollision, lambdacoll.globalIndex());
    return cPrimaryRequireTopo
             ? goodPureLambda->sliceBy(lambdasPerCollision, lambdacoll.globalIndex())
             : goodPrimaryLambda->sliceBy(lambdasPerCollision, lambdacoll.globalIndex());
  }

  // [Phase 12a] Per-cascade purity gate: combines truth (MC) + ITS-tracking
  // mode logic that was duplicated in 4 analyze functions. Returns true if
  // the cascade should be KEPT for filling. IsMC template param compiles
  // away the truth check on the data path.
  template <bool IsMC, typename TFlagsRow>
  bool passesCascadePurityGate(TFlagsRow const& fr)
  {
    bool truthOk = true;
    if constexpr (IsMC) {
      truthOk = !cRequireTrueCascade || fr.isTrueCascade();
    }
    bool itsRescue = (cItsTrackMode == lcorr_const::kItsTrackModeRescue && fr.isItsTracked());
    if (!truthOk && !itsRescue)
      return false;
    if (cItsTrackMode == lcorr_const::kItsTrackModeRequired && !fr.isItsTracked())
      return false;
    return true;
  }
  // Gen-level lambda slicing: use HashBy (unsorted-safe) via o2::framework::expressions
  // Simple approach: iterate all and match by collision index inline (no cache needed)

  using LambdaCollisionsExt = aod::LambdaCollisions;

  // Gen Xi/Omega: PDG filtering is done inline inside each process function.
  // (Two Filter declarations on the same table within one struct are AND-ed,
  //  so |pdg|==3312 AND |pdg|==3334 would always be false.)
  Preslice<aod::McParticles> mcParticlesPerMcCollision = aod::mcparticle::mcCollisionId;

  // [Phase 6] Cascade efficiency weight (mirror of useEff for the Λ side).
  // Returns 1/eff so multiplying it into the pair weight yields an
  // efficiency-corrected count. The default implementation returns 1.0
  // (no correction) — fill in the CCDB-load body once you have an
  // efficiency map for cascades. The kKind template parameter is 0 for Ξ
  // and 1 for Ω so the map can be specialised by particle.
  //
  // Usage in pair loops:
  //   float wCasc = getCascadeEfficiency<0/*Xi*/>(casc.sign(), casc.pt(), xiY);
  //   float wPair = wLam * wCasc;
  //
  // Suggested CCDB-load skeleton (uncomment + adapt when you have a map):
  //   if (!cUseCascEff) return 1.0f;
  //   static TList* effList = ccdb->getForTimeStamp<TList>(cCascEffPath, -1);
  //   const char* hName = (kKind == 0)
  //     ? (sign < 0 ? "hEffXiMinus" : "hEffXiPlus")
  //     : (sign < 0 ? "hEffOmegaMinus" : "hEffOmegaPlus");
  //   auto* h = static_cast<TH3*>(effList->FindObject(hName));
  //   double e = h->GetBinContent(h->FindBin(/*cent*/, pt, rap));
  //   return (e > 0) ? 1.0f / static_cast<float>(e) : 1.0f;
  template <int kKind /*0=Xi, 1=Omega*/>
  float getCascadeEfficiency(int /*sign*/, float /*pt*/, float /*rap*/)
  {
    if (!cUseCascEff)
      return 1.0f;
    // TODO: load from CCDB and compute weight when an efficiency map exists.
    return 1.0f;
  }

  // --- R2 Calculation Helper ---
  // R2 = (N_events * Pair_Yield) / (Single_Yield_1 * Single_Yield_2) - 1
  static TH2* calculateR2(TH2* hPairs, TH1* hSinglesTrig, TH1* hSinglesAssoc, double nEvents)
  {
    if (!hPairs || !hSinglesTrig || !hSinglesAssoc || nEvents <= 0)
      return nullptr;

    TH2* hR2 = reinterpret_cast<TH2*>(hPairs->Clone(Form("%s_R2", hPairs->GetName())));
    hR2->Reset();

    double nS1 = hSinglesTrig->Integral();
    double nS2 = hSinglesAssoc->Integral();

    if (nS1 > 0 && nS2 > 0) {
      hR2->Add(hPairs);
      hR2->Scale(nEvents / (nS1 * nS2));

      for (int i = 1; i <= hR2->GetNbinsX(); i++) {
        for (int j = 1; j <= hR2->GetNbinsY(); j++) {
          double content = hR2->GetBinContent(i, j);
          hR2->SetBinContent(i, j, content - 1.0);
        }
      }
    }
    return hR2;
  }

  void init(InitContext const&)
  {
    // [DEBUG] Echo the active configuration at startup. Search log for "[CFG-LXi]".
    LOGF(info,
         "[CFG-LXi] |y|<=%.2f useEff=%d useCascEff=%d tpcNsigmaBachKaon=%.2f "
         "vetoMode=%d primOnly=%d primReqTopo=%d reqTrueCasc=%d itsMode=%d "
         "saveCascTree=%d saveLambdaTree=%d "
         "pairs[LamXi=%d LamOm=%d LamLam=%d XiXi=%d OmOm=%d XiOm=%d] yields=%d "
         "(processXi=%d processOmega=%d processMCRecoXi=%d processMCRecoOmega=%d "
         "processAllPairs=%d processYields=%d "
         "processMCGenXi=%d processMCGenOmega=%d)",
         (float)maxY, (int)useEff, (int)cUseCascEff, (float)tpcNsigmaBachKaon,
         (int)cVetoMode, (int)cUsePrimaryLambdasOnly,
         (int)cPrimaryRequireTopo, (int)cRequireTrueCascade, (int)cItsTrackMode,
         (int)saveCascTree, (int)saveLambdaTree,
         (int)pairCfg.cFillLamXi, (int)pairCfg.cFillLamOm, (int)pairCfg.cFillLamLam,
         (int)pairCfg.cFillXiXi, (int)pairCfg.cFillOmOm, (int)pairCfg.cFillXiOm,
         (int)yieldCfg.cFillEventYields,
         (int)doprocessXi, (int)doprocessOmega,
         (int)doprocessMCRecoXi, (int)doprocessMCRecoOmega,
         (int)doprocessAllPairs, (int)doprocessYields,
         (int)doprocessMCGenXi, (int)doprocessMCGenOmega);

    // [Phase 6] Three-way veto policy. Loud warning when the user disables
    // it or picks the loose mode for a systematic study.
    if (cVetoMode == lcorr_const::kVetoModeOff) {
      LOGF(warning,
           "[LXi] cVetoMode=0 (off): pairs whose Λ shares V0 daughter tracks "
           "with the cascade's V0 will be KEPT. This produces a fake self-Λ "
           "peak around (Δφ≈0, Δy≈0). Use only for systematic studies.");
    } else if (cVetoMode == lcorr_const::kVetoModeLoose) {
      LOGF(info,
           "[LXi] cVetoMode=2 (loose): pairs vetoed when EITHER daughter is "
           "shared with the cascade's V0. Stricter than mode 1 — use to "
           "estimate the impact of partial-share contamination.");
    } else if (cVetoMode != lcorr_const::kVetoModeStrict) {
      LOGF(warning,
           "[LXi] cVetoMode=%d is not one of {0,1,2}; treating as 0 (off).",
           (int)cVetoMode);
    }
    // [Phase 6] Legacy cVetoSharedDau is parsed for back-compat but not used.
    // If the user set it to false on the assumption it'd disable the veto,
    // remind them to migrate.
    if (!cVetoSharedDau) {
      LOGF(warning,
           "[LXi] cVetoSharedDau=false is DEPRECATED and ignored at runtime. "
           "To disable the veto, set cVetoMode=0 instead. Current cVetoMode=%d.",
           (int)cVetoMode);
    }

    // [P4] Data-side primary tag is not informative.
    if (doprocessXi || doprocessOmega) {
      LOGF(warning,
           "[LXi] Running on DATA process functions: the Λ row's v0PrmScd "
           "field is unconditionally kPrimary regardless of physics origin "
           "(see LambdaCascadeProducer::fillLambdaRecoTables). Treat any "
           "primary/secondary partitioning as meaningful only in MC.");
      // [Phase 7] Loud reminder when the user wires primary-only on data.
      if (cUsePrimaryLambdasOnly) {
        LOGF(warning,
             "[LXi][Phase 7] cUsePrimaryLambdasOnly=true on data: the partition "
             "filters on v0PrmScd==kPrimary, but on data EVERY Λ has that tag. "
             "The trigger pool is therefore identical to the all-Λ pool. To "
             "actually enrich primaries on data, tighten the topological cuts "
             "in lambda-cascade-producer (cMaxDcaV0ToPV ~ 0.3-0.5 cm, "
             "cMinV0CosPA ~ 0.999) — those are not yet aggressive in your JSON.");
      }
    }

    // [P4][C1] Warnings for known foot-guns until the deeper fixes land.
    if (doprocessXi && doprocessOmega) {
      LOGF(warning,
           "[LXi] Both processXi and processOmega are enabled. Until the "
           "single-pass refactor (Phase 4), Lambda singles will be filled "
           "TWICE per event. Disable one of them, or wait for Phase 4.");
    }
    // Efficiency map sanity (P9): numerator and denominator come from
    // different process functions; warn if only one side is active.
    if (doprocessMCRecoXi && !doprocessMCGenXi) {
      LOGF(warning,
           "[LXi] processMCRecoXi is on but processMCGenXi is off; the Xi "
           "efficiency Gen denominator will be empty.");
    }
    if (doprocessMCGenXi && !doprocessMCRecoXi) {
      LOGF(warning,
           "[LXi] processMCGenXi is on but processMCRecoXi is off; the Xi "
           "efficiency Reco numerator will be empty.");
    }
    if (doprocessMCRecoOmega && !doprocessMCGenOmega) {
      LOGF(warning,
           "[LXi] processMCRecoOmega is on but processMCGenOmega is off; the "
           "Omega efficiency Gen denominator will be empty.");
    }
    if (doprocessMCGenOmega && !doprocessMCRecoOmega) {
      LOGF(warning,
           "[LXi] processMCGenOmega is on but processMCRecoOmega is off; the "
           "Omega efficiency Reco numerator will be empty.");
    }

    // --- 1. Axis Definitions ---
    // [Phase 6] Pair / single axes now come from the histAxes ConfigurableGroup
    // — change binning via JSON, no recompile.
    const AxisSpec dphi{histAxes.cDphiAxis, "#Delta#varphi"};
    const AxisSpec dy{histAxes.cDyAxis, "#Delta y"};
    const AxisSpec cent{centAxis, "FT0M (%)"};
    const AxisSpec pt{histAxes.cPtAxisPair, "p_{T} (GeV/c)"};
    const AxisSpec rap{histAxes.cRapAxis, "y"};
    const AxisSpec massLam{histAxes.cMassLamAxis, "M_{p#pi} (GeV/c^{2})"};
    const AxisSpec massXi{histAxes.cMassXiAxis, "M_{#Lambda#pi} (GeV/c^{2})"};
    const AxisSpec massOm{histAxes.cMassOmAxis, "M_{#LambdaK} (GeV/c^{2})"};
    const AxisSpec radius{100, 0, 100, "Radius (cm)"};
    const AxisSpec cpa{100, 0.9, 1.0, "Cos(PA)"};
    const AxisSpec dca{100, 0.0, 5.0, "DCA (cm)"};
    const AxisSpec pvDca{100, -10.0, 10.0, "DCA to PV (cm)"};
    // [C5] tpcRows AxisSpec removed alongside the unfilled QA/Casc/hTPCRows
    // bookings; will return in Phase 4 when the analyser carries tracks.

    // --- 2. Histograms ---
    histos.add("Event/hEventCount", "Event Counter", kTH1F, {{1, 0, 1, "Count"}});

    // Singles: Lambda
    histos.add("Singles/Lambda/hPt", "Lambda p_{T}", kTH1F, {pt});
    histos.add("Singles/AntiLambda/hPt", "AntiLambda p_{T}", kTH1F, {pt});

    histos.add("Singles/Lambda/hPtVsMass", "Lambda p_{T} vs Mass", kTH2F, {massLam, pt});
    histos.add("Singles/AntiLambda/hPtVsMass", "AntiLambda p_{T} vs Mass", kTH2F, {massLam, pt});

    // [Phase 11] 2D (y, φ) singles — input to ρ₁⊗ρ₁ for the offline R₂ recipe.
    const AxisSpec phi{histAxes.cPhiAxis, "φ (rad)"};
    histos.add("Singles/Lambda/hYPhi",     "Lambda (y, φ)",     kTH2F, {rap, phi});
    histos.add("Singles/AntiLambda/hYPhi", "AntiLambda (y, φ)", kTH2F, {rap, phi});

    // Singles: Xi & QA
    histos.add("QA/Xi/hRadius", "Xi Radius", kTH1F, {radius});
    histos.add("QA/Xi/hCosPA", "Xi CosPA", kTH1F, {cpa});
    histos.add("QA/Xi/hDCAV0Dau", "DCA V0 Daughters", kTH1F, {dca});
    histos.add("QA/Xi/hDCACascDau", "DCA Casc Daughters", kTH1F, {dca});
    histos.add("QA/Xi/hDCAV0ToPV", "DCA V0 to PV", kTH1F, {pvDca});
    histos.add("QA/Xi/hDCAPosToPV", "DCA Pos to PV", kTH1F, {pvDca});
    histos.add("QA/Xi/hDCANegToPV", "DCA Neg to PV", kTH1F, {pvDca});
    histos.add("QA/Xi/hDCABachToPV", "DCA Bach to PV", kTH1F, {pvDca});

    // [C5] QA/Casc/hTPCRows{Pos,Neg,Bach} were booked here but never filled.
    // To fill them properly we need track access (FullTracksExtIUWithPID)
    // inside the analyser, which processXi does not currently take. They are
    // re-introduced in Phase 4 once the analyser carries the tracks template.

    histos.add("Singles/XiMinus/hPtVsMass", "Xi^{-} p_{T} vs Mass", kTH2F, {massXi, pt});
    histos.add("Singles/XiPlus/hPtVsMass", "Xi^{+} p_{T} vs Mass", kTH2F, {massXi, pt});
    histos.add("Singles/XiMinus/hRap", "Xi^{-} Rapidity", kTH1F, {rap});
    histos.add("Singles/XiPlus/hRap", "Xi^{+} Rapidity", kTH1F, {rap});

    // [Phase 11] 2D (y, φ) singles for Ξ.
    histos.add("Singles/XiMinus/hYPhi", "Xi^{-} (y, φ)", kTH2F, {rap, phi});
    histos.add("Singles/XiPlus/hYPhi",  "Xi^{+} (y, φ)", kTH2F, {rap, phi});

    // --- Omega QA ---
    histos.add("QA/Om/hRadius", "Omega Radius", kTH1F, {radius});
    histos.add("QA/Om/hCosPA", "Omega CosPA", kTH1F, {cpa});
    histos.add("QA/Om/hDCAV0Dau", "DCA V0 Daughters", kTH1F, {dca});
    histos.add("QA/Om/hDCACascDau", "DCA Casc Daughters", kTH1F, {dca});
    histos.add("QA/Om/hDCAV0ToPV", "DCA V0 to PV", kTH1F, {pvDca});
    histos.add("QA/Om/hDCAPosToPV", "DCA Pos to PV", kTH1F, {pvDca});
    histos.add("QA/Om/hDCANegToPV", "DCA Neg to PV", kTH1F, {pvDca});
    histos.add("QA/Om/hDCABachToPV", "DCA Bach to PV", kTH1F, {pvDca});

    // Singles: Omega
    histos.add("Singles/OmegaMinus/hPtVsMass", "Omega^{-} p_{T} vs Mass", kTH2F, {massOm, pt});
    histos.add("Singles/OmegaPlus/hPtVsMass", "Omega^{+} p_{T} vs Mass", kTH2F, {massOm, pt});
    histos.add("Singles/OmegaMinus/hRap", "Omega^{-} Rapidity", kTH1F, {rap});
    histos.add("Singles/OmegaPlus/hRap", "Omega^{+} Rapidity", kTH1F, {rap});

    // [Phase 11] 2D (y, φ) singles for Ω.
    histos.add("Singles/OmegaMinus/hYPhi", "Omega^{-} (y, φ)", kTH2F, {rap, phi});
    histos.add("Singles/OmegaPlus/hYPhi",  "Omega^{+} (y, φ)", kTH2F, {rap, phi});

    // [Phase 13a] All pair-type bookings gated by pairCfg toggles.
    // Pairs: Λ–Ξ (R2 inputs)
    if (pairCfg.cFillLamXi) {
      histos.add("Pairs/Lam_XiM/hDeltaPhiDeltaY", "L-Xi-", kTH3F, {cent, dphi, dy});
      histos.add("Pairs/Lam_XiP/hDeltaPhiDeltaY", "L-Xi+", kTH3F, {cent, dphi, dy});
      histos.add("Pairs/AntiLam_XiM/hDeltaPhiDeltaY", "AL-Xi-", kTH3F, {cent, dphi, dy});
      histos.add("Pairs/AntiLam_XiP/hDeltaPhiDeltaY", "AL-Xi+", kTH3F, {cent, dphi, dy});
    }

    // Pairs: Λ–Ω (R2 inputs)
    if (pairCfg.cFillLamOm) {
      histos.add("Pairs/Lam_OmM/hDeltaPhiDeltaY", "L-Om-", kTH3F, {cent, dphi, dy});
      histos.add("Pairs/Lam_OmP/hDeltaPhiDeltaY", "L-Om+", kTH3F, {cent, dphi, dy});
      histos.add("Pairs/AntiLam_OmM/hDeltaPhiDeltaY", "AL-Om-", kTH3F, {cent, dphi, dy});
      histos.add("Pairs/AntiLam_OmP/hDeltaPhiDeltaY", "AL-Om+", kTH3F, {cent, dphi, dy});
    }

    // Pairs: Λ–Λ same-event (3 sign combos)
    if (pairCfg.cFillLamLam) {
      histos.add("Pairs/Lam_Lam/hDeltaPhiDeltaY",         "L-L",     kTH3F, {cent, dphi, dy});
      histos.add("Pairs/Lam_AntiLam/hDeltaPhiDeltaY",     "L-AL",    kTH3F, {cent, dphi, dy});
      histos.add("Pairs/AntiLam_AntiLam/hDeltaPhiDeltaY", "AL-AL",   kTH3F, {cent, dphi, dy});
    }

    // Pairs: Ξ–Ξ same-event (3 sign combos)
    if (pairCfg.cFillXiXi) {
      histos.add("Pairs/XiM_XiM/hDeltaPhiDeltaY", "Xi--Xi-", kTH3F, {cent, dphi, dy});
      histos.add("Pairs/XiM_XiP/hDeltaPhiDeltaY", "Xi--Xi+", kTH3F, {cent, dphi, dy});
      histos.add("Pairs/XiP_XiP/hDeltaPhiDeltaY", "Xi+-Xi+", kTH3F, {cent, dphi, dy});
    }

    // Pairs: Ω–Ω same-event (3 sign combos)
    if (pairCfg.cFillOmOm) {
      histos.add("Pairs/OmM_OmM/hDeltaPhiDeltaY", "Om--Om-", kTH3F, {cent, dphi, dy});
      histos.add("Pairs/OmM_OmP/hDeltaPhiDeltaY", "Om--Om+", kTH3F, {cent, dphi, dy});
      histos.add("Pairs/OmP_OmP/hDeltaPhiDeltaY", "Om+-Om+", kTH3F, {cent, dphi, dy});
    }

    // Pairs: Ξ–Ω cross-species (4 sign combos)
    if (pairCfg.cFillXiOm) {
      histos.add("Pairs/XiM_OmM/hDeltaPhiDeltaY", "Xi--Om-", kTH3F, {cent, dphi, dy});
      histos.add("Pairs/XiM_OmP/hDeltaPhiDeltaY", "Xi--Om+", kTH3F, {cent, dphi, dy});
      histos.add("Pairs/XiP_OmM/hDeltaPhiDeltaY", "Xi+-Om-", kTH3F, {cent, dphi, dy});
      histos.add("Pairs/XiP_OmP/hDeltaPhiDeltaY", "Xi+-Om+", kTH3F, {cent, dphi, dy});
    }

    // [Phase 12b] Like-sign / opposite-sign baryon-number combined pair
    // histograms. The four sign combos above split into:
    //   LS_LamXi (BB or B̄B̄) = Lam_XiM + AntiLam_XiP
    //   OS_LamXi (BB̄ or B̄B) = Lam_XiP + AntiLam_XiM
    //   (same for Λ-Ω)
    // Filled in the pair loop alongside the four-way split, so users can
    // do baryon-number correlation studies without offline summing.
    histos.add("Pairs/LS_LamXi/hDeltaPhiDeltaY", "LS Λ-Ξ", kTH3F, {cent, dphi, dy});
    histos.add("Pairs/OS_LamXi/hDeltaPhiDeltaY", "OS Λ-Ξ", kTH3F, {cent, dphi, dy});
    histos.add("Pairs/LS_LamOm/hDeltaPhiDeltaY", "LS Λ-Ω", kTH3F, {cent, dphi, dy});
    histos.add("Pairs/OS_LamOm/hDeltaPhiDeltaY", "OS Λ-Ω", kTH3F, {cent, dphi, dy});

    // [Phase 12b] Per-trigger-pT pair histograms (3D: pT_Λ × Δφ × Δy).
    // Centrality is integrated here to keep the storage manageable. Fill
    // gated by cFillPtDifferentialPairs so users opt in only when needed.
    if (cFillPtDifferentialPairs) {
      histos.add("PairsPt/Lam_XiM/hPtDeltaPhiDeltaY", "L-Xi- (pT_L, Δφ, Δy)", kTH3F, {pt, dphi, dy});
      histos.add("PairsPt/Lam_XiP/hPtDeltaPhiDeltaY", "L-Xi+ (pT_L, Δφ, Δy)", kTH3F, {pt, dphi, dy});
      histos.add("PairsPt/AntiLam_XiM/hPtDeltaPhiDeltaY", "AL-Xi- (pT_AL, Δφ, Δy)", kTH3F, {pt, dphi, dy});
      histos.add("PairsPt/AntiLam_XiP/hPtDeltaPhiDeltaY", "AL-Xi+ (pT_AL, Δφ, Δy)", kTH3F, {pt, dphi, dy});
      histos.add("PairsPt/Lam_OmM/hPtDeltaPhiDeltaY", "L-Om- (pT_L, Δφ, Δy)", kTH3F, {pt, dphi, dy});
      histos.add("PairsPt/Lam_OmP/hPtDeltaPhiDeltaY", "L-Om+ (pT_L, Δφ, Δy)", kTH3F, {pt, dphi, dy});
      histos.add("PairsPt/AntiLam_OmM/hPtDeltaPhiDeltaY", "AL-Om- (pT_AL, Δφ, Δy)", kTH3F, {pt, dphi, dy});
      histos.add("PairsPt/AntiLam_OmP/hPtDeltaPhiDeltaY", "AL-Om+ (pT_AL, Δφ, Δy)", kTH3F, {pt, dphi, dy});
    }

    // [P3] Auto-correlation bookkeeping: per-event distributions of total
    // pair count vs vetoed (shared-daughter) pair count. Compare the two
    // to estimate the bias of running without the veto.
    const AxisSpec axisPairCount{200, 0., 200., "pairs/event"};
    histos.add("QA/AutoCorr/hXiPairsTotal",   "Λ-Ξ total pairs/event",   kTH1F, {axisPairCount});
    histos.add("QA/AutoCorr/hXiPairsVetoed",  "Λ-Ξ vetoed pairs/event",  kTH1F, {axisPairCount});
    histos.add("QA/AutoCorr/hOmPairsTotal",   "Λ-Ω total pairs/event",   kTH1F, {axisPairCount});
    histos.add("QA/AutoCorr/hOmPairsVetoed",  "Λ-Ω vetoed pairs/event",  kTH1F, {axisPairCount});

    // [Phase 9] ITS strangeness-tracking diagnostics. Filled regardless
    // of cItsTrackMode so the comparison is always available downstream.
    // X axis: 0 = topology-passing cascade, 1 = also IsItsTracked.
    // Reading: bin 1 / bin 0 = ITS-tracking efficiency × "real Ξ" purity.
    const AxisSpec axisItsTrackBin{2, -0.5, 1.5, "ITS-tracked"};
    const AxisSpec axisCascPt{50, 0, 10, "p_{T} (GeV/c)"};
    histos.add("QA/ItsTrack/hXiTotalVsTracked",   "Ξ topology-pass vs ITS-tracked",       kTH2F, {axisItsTrackBin, axisCascPt});
    histos.add("QA/ItsTrack/hOmegaTotalVsTracked","Ω topology-pass vs ITS-tracked",       kTH2F, {axisItsTrackBin, axisCascPt});

    // [Phase 13b] Per-event yield histograms — one block per species,
    // 5 histos each: hNPerEvent, hMeanPtPerEvent, hMeanNvsCent,
    // hMeanPtVsCent, hNvsPt2D. Filled exclusively by processYields so
    // counts are correct regardless of which pair-process functions run.
    if (yieldCfg.cFillEventYields) {
      const AxisSpec axisYieldN {yieldCfg.cYieldNAxis,  "N per event"};
      const AxisSpec axisYieldPt{yieldCfg.cYieldPtAxis, "⟨p_{T}⟩ per event (GeV/c)"};
      const std::array<const char*, 6> species{"Lambda", "AntiLambda", "XiMinus", "XiPlus", "OmegaMinus", "OmegaPlus"};
      for (const auto& sp : species) {
        histos.add(Form("Yields/%s/hNPerEvent",      sp), Form("%s N/event", sp),         kTH1F,    {axisYieldN});
        histos.add(Form("Yields/%s/hMeanPtPerEvent", sp), Form("%s ⟨pT⟩/event", sp),     kTH1F,    {axisYieldPt});
        histos.add(Form("Yields/%s/hMeanNvsCent",    sp), Form("%s ⟨N⟩ vs cent", sp),    kTProfile,{cent});
        histos.add(Form("Yields/%s/hMeanPtVsCent",   sp), Form("%s ⟨pT⟩ vs cent", sp),   kTProfile,{cent});
        histos.add(Form("Yields/%s/hNvsPt2D",        sp), Form("%s (N, ⟨pT⟩)", sp),      kTH2F,    {axisYieldN, axisYieldPt});
      }
    }

    // --- MC Gen-level histograms (mirror reco under McGen/) ---
    // Singles: gen Lambda
    histos.add("McGen/Singles/Lambda/hPt", "Gen #Lambda p_{T}", kTH1F, {pt});
    histos.add("McGen/Singles/AntiLambda/hPt", "Gen #bar{#Lambda} p_{T}", kTH1F, {pt});

    // Singles: gen Xi
    histos.add("McGen/Singles/XiMinus/hPtVsRap", "Gen #Xi^{-} p_{T} vs y", kTH2F, {rap, pt});
    histos.add("McGen/Singles/XiPlus/hPtVsRap", "Gen #Xi^{+} p_{T} vs y", kTH2F, {rap, pt});

    // Singles: gen Omega
    histos.add("McGen/Singles/OmegaMinus/hPtVsRap", "Gen #Omega^{-} p_{T} vs y", kTH2F, {rap, pt});
    histos.add("McGen/Singles/OmegaPlus/hPtVsRap", "Gen #Omega^{+} p_{T} vs y", kTH2F, {rap, pt});

    // [Phase 16t] Closure-test inputs: 2D (y, φ) singles for the ρ₁⊗ρ₁
    // convolution. compute_r2.C consumes these via the --gen switch.
    // Mirrors the reco-side Singles/.../hYPhi histograms used in production.
    histos.add("McGen/Singles/Lambda/hYPhi",     "Gen #Lambda (y, #varphi)",     kTH2F, {rap, phi});
    histos.add("McGen/Singles/AntiLambda/hYPhi", "Gen #bar{#Lambda} (y, #varphi)", kTH2F, {rap, phi});
    histos.add("McGen/Singles/XiMinus/hYPhi",    "Gen #Xi^{-} (y, #varphi)",     kTH2F, {rap, phi});
    histos.add("McGen/Singles/XiPlus/hYPhi",     "Gen #Xi^{+} (y, #varphi)",     kTH2F, {rap, phi});
    histos.add("McGen/Singles/OmegaMinus/hYPhi", "Gen #Omega^{-} (y, #varphi)",  kTH2F, {rap, phi});
    histos.add("McGen/Singles/OmegaPlus/hYPhi",  "Gen #Omega^{+} (y, #varphi)",  kTH2F, {rap, phi});

    // [Phase 16t] Gen event counter (canonical denominator for gen-level R₂).
    histos.add("McGen/Event/hEventCount", "Gen Event Counter", kTH1F, {{1, 0, 1, "Count"}});

    // Pairs: gen Lam-Xi
    histos.add("McGen/Pairs/Lam_XiM/hDeltaPhiDeltaY", "Gen L-Xi-", kTH3F, {cent, dphi, dy});
    histos.add("McGen/Pairs/Lam_XiP/hDeltaPhiDeltaY", "Gen L-Xi+", kTH3F, {cent, dphi, dy});
    histos.add("McGen/Pairs/AntiLam_XiM/hDeltaPhiDeltaY", "Gen AL-Xi-", kTH3F, {cent, dphi, dy});
    histos.add("McGen/Pairs/AntiLam_XiP/hDeltaPhiDeltaY", "Gen AL-Xi+", kTH3F, {cent, dphi, dy});

    // Pairs: gen Lam-Omega
    histos.add("McGen/Pairs/Lam_OmM/hDeltaPhiDeltaY", "Gen L-Om-", kTH3F, {cent, dphi, dy});
    histos.add("McGen/Pairs/Lam_OmP/hDeltaPhiDeltaY", "Gen L-Om+", kTH3F, {cent, dphi, dy});
    histos.add("McGen/Pairs/AntiLam_OmM/hDeltaPhiDeltaY", "Gen AL-Om-", kTH3F, {cent, dphi, dy});
    histos.add("McGen/Pairs/AntiLam_OmP/hDeltaPhiDeltaY", "Gen AL-Om+", kTH3F, {cent, dphi, dy});

    // [Phase 16a] Per-event coincidence counters: how often does a given
    // collision contain BOTH a primary-Λ trigger AND a reco cascade? This
    // is the dominant statistical bottleneck on small AODs (the pair loop
    // only fires when both species coexist), so we instrument it directly.
    // Read-off:
    //   bin(0,0)        → events with neither species
    //   row > 0, col 0  → events with only Λ
    //   row 0, col > 0  → events with only cascade
    //   row > 0, col > 0→ events that contribute pair candidates
    const AxisSpec axisCoincLam{31, -0.5,  30.5, "N_{primary Λ}"};
    const AxisSpec axisCoincXi {16, -0.5,  15.5, "N_{Ξ}"};
    const AxisSpec axisCoincOm {16, -0.5,  15.5, "N_{Ω}"};
    histos.add("Yields/Coincidence/hLamXi", "events: (N_Λ, N_Ξ)", kTH2F, {axisCoincLam, axisCoincXi});
    histos.add("Yields/Coincidence/hLamOm", "events: (N_Λ, N_Ω)", kTH2F, {axisCoincLam, axisCoincOm});

    // --- 3. Efficiency map production histograms ---
    // Axis definitions for efficiency maps (coarser binning to avoid empty bins)
    const AxisSpec effCent{cent};
    const AxisSpec effPt{50, 0, 10, "p_{T} (GeV/c)"};
    const AxisSpec effRap{20, -1.0, 1.0, "y"};

    // Generator-level denominator: all physical-primary Xi/Omega in acceptance
    histos.add("Eff/Gen/XiMinus/hCentPtRap", "Gen #Xi^{-}", kTH3F, {effCent, effPt, effRap});
    histos.add("Eff/Gen/XiPlus/hCentPtRap", "Gen #Xi^{+}", kTH3F, {effCent, effPt, effRap});
    histos.add("Eff/Gen/OmegaMinus/hCentPtRap", "Gen #Omega^{-}", kTH3F, {effCent, effPt, effRap});
    histos.add("Eff/Gen/OmegaPlus/hCentPtRap", "Gen #Omega^{+}", kTH3F, {effCent, effPt, effRap});
    // Reco-level numerator: truth-matched reco Xi/Omega (physical-primary)
    histos.add("Eff/Reco/XiMinus/hCentPtRap", "Reco #Xi^{-}", kTH3F, {effCent, effPt, effRap});
    histos.add("Eff/Reco/XiPlus/hCentPtRap", "Reco #Xi^{+}", kTH3F, {effCent, effPt, effRap});
    histos.add("Eff/Reco/OmegaMinus/hCentPtRap", "Reco #Omega^{-}", kTH3F, {effCent, effPt, effRap});
    histos.add("Eff/Reco/OmegaPlus/hCentPtRap", "Reco #Omega^{+}", kTH3F, {effCent, effPt, effRap});

    // TTree setup — only allocates when at least one of saveCascTree or
    // saveLambdaTree is enabled. OutputObj<TTree> are direct task members,
    // so the framework automatically writes them into AnalysisResults.root.
    if (saveCascTree || saveLambdaTree) {
      bp = new lxicorr::BranchPair{};
    }
    if (saveCascTree) {
      treeXi.setObject(new TTree("XiCandidates", "Xi topological variables"));
      treeOmega.setObject(new TTree("OmegaCandidates", "Omega topological variables"));
      lxicorr::connectBranches(treeXi.object.get(), &bp->xi);
      lxicorr::connectBranches(treeOmega.object.get(), &bp->om);
      treeXiGen.setObject(new TTree("XiCandidatesGen", "Xi gen-level kinematics"));
      treeOmegaGen.setObject(new TTree("OmegaCandidatesGen", "Omega gen-level kinematics"));
      lxicorr::connectGenBranches(treeXiGen.object.get(), &bp->xiGen);
      lxicorr::connectGenBranches(treeOmegaGen.object.get(), &bp->omGen);
    }
    // [Phase 7] Λ TTrees.
    if (saveLambdaTree) {
      treeLambda.setObject(new TTree("LambdaCandidates", "Λ standalone-V0 topo + truth"));
      treeLambdaGen.setObject(new TTree("LambdaCandidatesGen", "Λ gen-level kinematics"));
      lxicorr::connectLambdaBranches(treeLambda.object.get(), &bp->lam);
      lxicorr::connectLambdaGenBranches(treeLambdaGen.object.get(), &bp->lamGen);
    }
  }

  // --- Analysis Functions ---

  // [Phase 7] Optional: a per-event PV used as the centVal/pvZ source for the
  // Λ tree. We get these from the LambdaCollision the caller passes in.
  template <typename T, typename C>
  void analyzeSinglesLambda(T const& tracks, C const& lambdacoll)
  {
    // [Phase 16f] Idempotence guard — see mLastLamSinglesCollIdx declaration.
    // Returns early on the 2nd+ call within the same event so that running
    // both processMCRecoXi and processMCRecoOmega doesn't double-fill the
    // Λ singles histograms or LambdaCandidates tree.
    const int64_t thisCollIdx = lambdacoll.globalIndex();
    if (thisCollIdx == mLastLamSinglesCollIdx) {
      return;
    }
    mLastLamSinglesCollIdx = thisCollIdx;

    for (const auto& track : tracks) {
      if (std::abs(track.rap()) > maxY)
        continue;

      float w = useEff ? track.corrFact() : 1.0f;
      // [P6] Use the typed enum instead of magic literal "1".
      bool isAnti = (track.v0Type() == (int8_t)kAntiLambda);

      // [Phase 11] φ wrapped to [0, 2π) for ρ₁(y, φ) accumulation.
      float phiWrapped = RecoDecay::constrainAngle(track.phi(), 0.f);
      if (!isAnti) {
        histos.fill(HIST("Singles/Lambda/hPt"), track.pt(), w);
        histos.fill(HIST("Singles/Lambda/hPtVsMass"), track.mass(), track.pt(), w);
        histos.fill(HIST("Singles/Lambda/hYPhi"), track.rap(), phiWrapped, w);
      } else {
        histos.fill(HIST("Singles/AntiLambda/hPt"), track.pt(), w);
        histos.fill(HIST("Singles/AntiLambda/hPtVsMass"), track.mass(), track.pt(), w);
        histos.fill(HIST("Singles/AntiLambda/hYPhi"), track.rap(), phiWrapped, w);
      }

      // [Phase 7+8 fix] Λ TTree fill — must AND with saveLambdaTree because bp
      // is shared with the cascade trees. With saveCascTree=true & saveLambdaTree=false
      // bp is allocated but treeLambda is NEVER setObject-ed; reaching Fill()
      // dereferences a null TTree and segfaults.
      if (bp && saveLambdaTree) {
        auto& b = bp->lam;
        b.pt = track.pt();
        b.eta = track.eta();
        b.rap = track.rap();
        b.phi = track.phi();
        b.mass = track.mass();
        b.cosPA = track.cosPA();
        b.dcaDau = track.dcaDau();
        b.v0Type = static_cast<int>(track.v0Type());
        b.v0PrmScd = static_cast<int>(track.v0PrmScd());
        b.trueLambdaFlag = track.trueLambdaFlag();
        b.corrFact = track.corrFact();
        b.cent = lambdacoll.cent();
        b.pvZ = lambdacoll.posZ();
        b.motherPdg = track.motherPdg();
        // [Phase 16a] Populate pdgCode from the truth tag + v0Type so downstream
        // selectors of the form `abs(pdgCode)==3122` work uniformly with the
        // cascade trees. On data trueLambdaFlag falls back to topology, which
        // makes pdgCode best-effort but never wrong: it is non-zero only when
        // we have a real Λ-tagged candidate.
        if (track.trueLambdaFlag()) {
          b.pdgCode = (track.v0Type() == (int8_t)kLambda)
                          ? lcorr_const::kLambdaPdg
                          : -lcorr_const::kLambdaPdg;
        } else {
          b.pdgCode = 0;
        }
        b.isPhysPrim = (track.v0PrmScd() == (int8_t)kPrimary);
        // [Phase 8] Topology snapshot — primary-fraction template-fit inputs.
        b.dcaV0ToPV = track.dcaV0ToPV();
        b.v0Radius = track.v0Radius();
        b.posItsNCls = static_cast<int>(track.posItsNCls());
        b.negItsNCls = static_cast<int>(track.negItsNCls());
        b.passesPrimaryTopo = track.passesPrimaryTopo();
        // [Phase 9] ITS hit-map per daughter for downstream geometric checks.
        b.posItsClusterMap = static_cast<unsigned int>(track.posItsClusterMap());
        b.negItsClusterMap = static_cast<unsigned int>(track.negItsClusterMap());
        // [Phase 10] L_proper + per-daughter DCAs.
        b.lProper  = track.lProper();
        b.posDcaXY = track.posDcaXY();
        b.negDcaXY = track.negDcaXY();
        // [Phase 14] cut bitmask + raw cut inputs.
        b.cutBits        = static_cast<unsigned int>(track.cutBits());
        b.tpcNSigmaPosPr = track.tpcNSigmaPosPr();
        b.tpcNSigmaNegPi = track.tpcNSigmaNegPi();
        b.tpcNSigmaPosPi = track.tpcNSigmaPosPi();
        b.tpcNSigmaNegPr = track.tpcNSigmaNegPr();
        b.mK0Short       = track.mK0Short();
        b.qtArm          = track.qtArm();
        b.alphaArm       = track.alphaArm();
        b.cTau           = track.cTau();
        treeLambda->Fill();
      }
    }
  }

  template <bool IsMC = false, typename T, typename F>
  void analyzeSinglesXi(T const& cascades, F const& flagsStart, float pvX, float pvY, float pvZ, float centVal)
  {
    for (const auto& casc : cascades) {
      auto fr = flagsStart + casc.globalIndex();
      const int f = fr.isSelected();
      // [Phase 16h] Species filter — only Ξ-flagged cascades enter the Ξ
      // singles / tree fill. Previously, the loop accepted ANY non-rejected
      // cascade (flag 1/2/3), which meant flag==3 (Ω-only) candidates were
      // silently classified as Ξ and their mass/rapidity projected under the
      // wrong mass hypothesis polluted XiCandidates.
      if (f != kFlagXiOnly && f != kFlagXiAndOmega)
        continue;

      // [Phase 9] QA: count topology-passers and the ITS-tracked subset.
      // Filled BEFORE any further gating so the ratio is meaningful.
      histos.fill(HIST("QA/ItsTrack/hXiTotalVsTracked"), 0., casc.pt());
      if (fr.isItsTracked())
        histos.fill(HIST("QA/ItsTrack/hXiTotalVsTracked"), 1., casc.pt());

      // [Phase 12a] Centralised purity gate (truth + ITS-tracking modes 1/2).
      if (!passesCascadePurityGate<IsMC>(fr))
        continue;

      // [P2] Cascade row's own rapidity uses the correct charged-Xi mass.
      // Avoid recomputing with MassXi0 (the neutral-cascade mass).
      float xiY = casc.yXi();
      if (std::abs(xiY) > maxY)
        continue;

      // QA Filling
      histos.fill(HIST("QA/Xi/hRadius"), casc.cascradius());
      histos.fill(HIST("QA/Xi/hCosPA"), casc.casccosPA(pvX, pvY, pvZ));
      histos.fill(HIST("QA/Xi/hDCAV0Dau"), casc.dcaV0daughters());
      histos.fill(HIST("QA/Xi/hDCACascDau"), casc.dcacascdaughters());
      histos.fill(HIST("QA/Xi/hDCAV0ToPV"), casc.dcav0topv(pvX, pvY, pvZ));
      histos.fill(HIST("QA/Xi/hDCAPosToPV"), casc.dcapostopv());
      histos.fill(HIST("QA/Xi/hDCANegToPV"), casc.dcanegtopv());
      histos.fill(HIST("QA/Xi/hDCABachToPV"), casc.dcabachtopv());

      // [Phase 11] φ wrapped to [0, 2π) for ρ₁(y, φ) accumulation.
      float xiPhiWrapped = RecoDecay::constrainAngle(casc.phi(), 0.f);
      if (casc.sign() < 0) {
        histos.fill(HIST("Singles/XiMinus/hPtVsMass"), casc.mXi(), casc.pt());
        histos.fill(HIST("Singles/XiMinus/hRap"), xiY);
        histos.fill(HIST("Singles/XiMinus/hYPhi"), xiY, xiPhiWrapped);
      } else {
        histos.fill(HIST("Singles/XiPlus/hPtVsMass"), casc.mXi(), casc.pt());
        histos.fill(HIST("Singles/XiPlus/hRap"), xiY);
        histos.fill(HIST("Singles/XiPlus/hYPhi"), xiY, xiPhiWrapped);
      }

      // MC: fill efficiency numerator for truth-matched primary Xi
      if constexpr (IsMC) {
        if (casc.has_mcParticle()) {
          auto mcpart = casc.mcParticle();
          if (std::abs(mcpart.pdgCode()) == lcorr_const::kXiMinusPdg && mcpart.isPhysicalPrimary()) {
            if (mcpart.pdgCode() == lcorr_const::kXiMinusPdg)
              histos.fill(HIST("Eff/Reco/XiMinus/hCentPtRap"), centVal, casc.pt(), xiY);
            else
              histos.fill(HIST("Eff/Reco/XiPlus/hCentPtRap"), centVal, casc.pt(), xiY);
          }
        }
      }

      // [Phase 8 fix] Same gating issue as the Λ tree: bp is shared between
      // cascade and Λ trees, so guard the Xi fill on saveCascTree explicitly.
      if (bp && saveCascTree) {
        auto& b = bp->xi;
        b.pt = casc.pt();
        b.rap = xiY;
        b.mass = casc.mXi();
        b.sign = casc.sign();
        b.cascCosPA = casc.casccosPA(pvX, pvY, pvZ);
        b.v0CosPA = casc.v0cosPA(pvX, pvY, pvZ);
        b.cascRadius = casc.cascradius();
        b.v0Radius = casc.v0radius();
        b.dcaV0Dau = casc.dcaV0daughters();
        b.dcaCascDau = casc.dcacascdaughters();
        b.dcaV0ToPV = casc.dcav0topv(pvX, pvY, pvZ);
        b.dcaPosToPV = casc.dcapostopv();
        b.dcaNegToPV = casc.dcanegtopv();
        b.dcaBachToPV = casc.dcabachtopv();
        b.cent = centVal;
        b.pvZ = pvZ;
        // MC truth matching
        if constexpr (IsMC) {
          if (casc.has_mcParticle()) {
            auto mcpart = casc.mcParticle();
            b.pdgCode = mcpart.pdgCode();
            b.isPhysPrim = mcpart.isPhysicalPrimary();
          } else {
            b.pdgCode = 0;
            b.isPhysPrim = false;
          }
        } else {
          b.pdgCode = 0;
          b.isPhysPrim = false;
        }
        // [Phase 14] Per-cascade cut bitmask + bachelor PID inputs.
        b.cascCutBits = static_cast<unsigned int>(fr.cascCutBits());
        auto bachTrk = casc.template bachelor_as<FullTracksExtIUWithPID>();
        b.bachTpcNSigmaPi        = bachTrk.tpcNSigmaPi();
        b.bachTpcNSigmaKa        = bachTrk.tpcNSigmaKa();
        b.bachItsNCls            = static_cast<int>(bachTrk.itsNCls());
        b.bachTpcNClsCrossedRows = static_cast<int>(bachTrk.tpcNClsCrossedRows());
        b.mLambdaInside          = casc.mLambda();
        treeXi->Fill();
      }
    }
  }

  // Omega singles: kaon bachelor PID cut distinguishes Omega from Xi
  template <bool IsMC = false, typename T, typename F>
  void analyzeSinglesOmega(T const& cascades, F const& flagsStart, float pvX, float pvY, float pvZ, float centVal)
  {
    for (const auto& casc : cascades) {
      auto fr = flagsStart + casc.globalIndex();
      const int f = fr.isSelected();
      // [Phase 16h] Species filter — only Ω-flagged cascades enter the Ω
      // singles / tree fill. Previously, the loop accepted any non-rejected
      // cascade, so flag==1 (Ξ-only) candidates leaked in. The kaon nσ cut
      // below filtered most of them, but a defensive flag gate is cheaper
      // and cleaner.
      if (f != kFlagOmegaOnly && f != kFlagXiAndOmega)
        continue;

      // [Phase 9] QA: total Ω topology-passers and the ITS-tracked subset.
      histos.fill(HIST("QA/ItsTrack/hOmegaTotalVsTracked"), 0., casc.pt());
      if (fr.isItsTracked())
        histos.fill(HIST("QA/ItsTrack/hOmegaTotalVsTracked"), 1., casc.pt());

      // [Phase 12a] Centralised purity gate.
      if (!passesCascadePurityGate<IsMC>(fr))
        continue;

      // [P2] Cascade row's own rapidity (correct charged-Omega mass).
      float omY = casc.yOmega();
      if (std::abs(omY) > maxY)
        continue;

      // Require kaon bachelor PID — NSigma lives on the track, not the cascade row
      auto bachTrack = casc.template bachelor_as<FullTracksExtIUWithPID>();
      if (std::abs(bachTrack.tpcNSigmaKa()) > tpcNsigmaBachKaon)
        continue;

      // QA Filling (mirrors Xi QA under QA/Om/)
      histos.fill(HIST("QA/Om/hRadius"), casc.cascradius());
      histos.fill(HIST("QA/Om/hCosPA"), casc.casccosPA(pvX, pvY, pvZ));
      histos.fill(HIST("QA/Om/hDCAV0Dau"), casc.dcaV0daughters());
      histos.fill(HIST("QA/Om/hDCACascDau"), casc.dcacascdaughters());
      histos.fill(HIST("QA/Om/hDCAV0ToPV"), casc.dcav0topv(pvX, pvY, pvZ));
      histos.fill(HIST("QA/Om/hDCAPosToPV"), casc.dcapostopv());
      histos.fill(HIST("QA/Om/hDCANegToPV"), casc.dcanegtopv());
      histos.fill(HIST("QA/Om/hDCABachToPV"), casc.dcabachtopv());

      // [Phase 11] φ wrapped to [0, 2π).
      float omPhiWrapped = RecoDecay::constrainAngle(casc.phi(), 0.f);
      if (casc.sign() < 0) {
        histos.fill(HIST("Singles/OmegaMinus/hPtVsMass"), casc.mOmega(), casc.pt());
        histos.fill(HIST("Singles/OmegaMinus/hRap"), omY);
        histos.fill(HIST("Singles/OmegaMinus/hYPhi"), omY, omPhiWrapped);
      } else {
        histos.fill(HIST("Singles/OmegaPlus/hPtVsMass"), casc.mOmega(), casc.pt());
        histos.fill(HIST("Singles/OmegaPlus/hRap"), omY);
        histos.fill(HIST("Singles/OmegaPlus/hYPhi"), omY, omPhiWrapped);
      }

      // MC: fill efficiency numerator for truth-matched primary Omega
      if constexpr (IsMC) {
        if (casc.has_mcParticle()) {
          auto mcpart = casc.mcParticle();
          if (std::abs(mcpart.pdgCode()) == lcorr_const::kOmegaMinusPdg && mcpart.isPhysicalPrimary()) {
            if (mcpart.pdgCode() == lcorr_const::kOmegaMinusPdg)
              histos.fill(HIST("Eff/Reco/OmegaMinus/hCentPtRap"), centVal, casc.pt(), omY);
            else
              histos.fill(HIST("Eff/Reco/OmegaPlus/hCentPtRap"), centVal, casc.pt(), omY);
          }
        }
      }

      // [Phase 8 fix] Gate on saveCascTree (treeOmega is only setObject-ed when on).
      if (bp && saveCascTree) {
        auto& b = bp->om;
        b.pt = casc.pt();
        b.rap = omY;
        b.mass = casc.mOmega();
        b.sign = casc.sign();
        b.cascCosPA = casc.casccosPA(pvX, pvY, pvZ);
        b.v0CosPA = casc.v0cosPA(pvX, pvY, pvZ);
        b.cascRadius = casc.cascradius();
        b.v0Radius = casc.v0radius();
        b.dcaV0Dau = casc.dcaV0daughters();
        b.dcaCascDau = casc.dcacascdaughters();
        b.dcaV0ToPV = casc.dcav0topv(pvX, pvY, pvZ);
        b.dcaPosToPV = casc.dcapostopv();
        b.dcaNegToPV = casc.dcanegtopv();
        b.dcaBachToPV = casc.dcabachtopv();
        b.cent = centVal;
        b.pvZ = pvZ;
        // MC truth matching
        if constexpr (IsMC) {
          if (casc.has_mcParticle()) {
            auto mcpart = casc.mcParticle();
            b.pdgCode = mcpart.pdgCode();
            b.isPhysPrim = mcpart.isPhysicalPrimary();
          } else {
            b.pdgCode = 0;
            b.isPhysPrim = false;
          }
        } else {
          b.pdgCode = 0;
          b.isPhysPrim = false;
        }
        // [Phase 16e fix] Per-cascade cut bitmask + bachelor PID inputs.
        // These were assigned on the Xi tree path but forgotten on Ω, so
        // OmegaCandidates branches cascCutBits / bach* / mLambdaInside
        // were silently always 0. Now mirrors the Xi-side write block.
        b.cascCutBits = static_cast<unsigned int>(fr.cascCutBits());
        b.bachTpcNSigmaPi        = bachTrack.tpcNSigmaPi();
        b.bachTpcNSigmaKa        = bachTrack.tpcNSigmaKa();
        b.bachItsNCls            = static_cast<int>(bachTrack.itsNCls());
        b.bachTpcNClsCrossedRows = static_cast<int>(bachTrack.tpcNClsCrossedRows());
        b.mLambdaInside          = casc.mLambda();
        treeOmega->Fill();
      }
    }
  }

  template <bool IsMC = false, typename L, typename C, typename F>
  void analyzePairs(L const& lambdas, C const& cascades, F const& flagsStart, float centVal)
  {
    // [P3] Per-event auto-correlation bookkeeping for Λ-Ξ.
    int nXiPairsTotal = 0;
    int nXiPairsVetoed = 0;

    for (const auto& lam : lambdas) {
      if (std::abs(lam.rap()) > maxY)
        continue;
      float wLam = useEff ? lam.corrFact() : 1.0f;
      // [P6] Use the typed enum instead of magic literal "1".
      bool isAntiLam = (lam.v0Type() == (int8_t)kAntiLambda);

      for (const auto& casc : cascades) {
        auto fr = flagsStart + casc.globalIndex();
        const int f = fr.isSelected();
        // [Phase 16h] Species filter — only Ξ-flagged cascades pair as Ξ.
        if (f != kFlagXiOnly && f != kFlagXiAndOmega)
          continue;

        // [Phase 12a] Centralised purity gate.
        if (!passesCascadePurityGate<IsMC>(fr))
          continue;

        // [P2] Use the cascade row's own Xi rapidity (correct mass).
        float xiY = casc.yXi();
        if (std::abs(xiY) > maxY)
          continue;

        ++nXiPairsTotal;

        // [Phase 6] Three-way auto-correlation veto driven by cVetoMode:
        //   1 (strict, default) — both pos AND neg daughters shared
        //   2 (loose)           — either pos OR neg daughter shared
        //   0 (off)             — keep all pairs (use only for systematic studies)
        const bool posMatch = (lam.posTrackId() == casc.posTrackId());
        const bool negMatch = (lam.negTrackId() == casc.negTrackId());
        const bool veto =
            (cVetoMode == lcorr_const::kVetoModeStrict && posMatch && negMatch) ||
            (cVetoMode == lcorr_const::kVetoModeLoose  && (posMatch || negMatch));
        if (veto) {
          ++nXiPairsVetoed;
          continue;
        }

        // [Phase 6] Per-cascade efficiency weight (scaffold). Default helper
        // returns 1.0; fill in CCDB loading inside getCascadeEfficiency<>
        // when you have an efficiency map.
        float wCascXi = getCascadeEfficiency<0 /*Xi*/>(casc.sign(), casc.pt(), xiY);
        float wPair = wLam * wCascXi;

        float dphi = RecoDecay::constrainAngle(casc.phi() - lam.phi(), -PIHalf);
        float dy = xiY - lam.rap();

        bool isXiPlus = (casc.sign() > 0);

        // [Phase 12b] Baryon-number labels:
        //   Λ:B=+1, Λ̄:B=-1, Ξ⁻:B=+1 (sign<0), Ξ⁺:B=-1 (sign>0)
        //   LS = same baryon number, OS = opposite.
        bool isLS = (isAntiLam == isXiPlus);  // (Λ̄, Ξ⁺) ↔ (Λ, Ξ⁻) both LS

        if (pairCfg.cFillLamXi) {
          if (!isAntiLam && !isXiPlus)
            histos.fill(HIST("Pairs/Lam_XiM/hDeltaPhiDeltaY"), centVal, dphi, dy, wPair);
          else if (!isAntiLam && isXiPlus)
            histos.fill(HIST("Pairs/Lam_XiP/hDeltaPhiDeltaY"), centVal, dphi, dy, wPair);
          else if (isAntiLam && !isXiPlus)
            histos.fill(HIST("Pairs/AntiLam_XiM/hDeltaPhiDeltaY"), centVal, dphi, dy, wPair);
          else if (isAntiLam && isXiPlus)
            histos.fill(HIST("Pairs/AntiLam_XiP/hDeltaPhiDeltaY"), centVal, dphi, dy, wPair);
        }

        // [Phase 12b] Baryon-LS / OS combined fills.
        if (isLS)
          histos.fill(HIST("Pairs/LS_LamXi/hDeltaPhiDeltaY"), centVal, dphi, dy, wPair);
        else
          histos.fill(HIST("Pairs/OS_LamXi/hDeltaPhiDeltaY"), centVal, dphi, dy, wPair);

        // [Phase 12b] pT-differential pair fill (opt-in via Configurable).
        if (cFillPtDifferentialPairs) {
          float ptL = lam.pt();
          if (!isAntiLam && !isXiPlus)
            histos.fill(HIST("PairsPt/Lam_XiM/hPtDeltaPhiDeltaY"), ptL, dphi, dy, wPair);
          else if (!isAntiLam && isXiPlus)
            histos.fill(HIST("PairsPt/Lam_XiP/hPtDeltaPhiDeltaY"), ptL, dphi, dy, wPair);
          else if (isAntiLam && !isXiPlus)
            histos.fill(HIST("PairsPt/AntiLam_XiM/hPtDeltaPhiDeltaY"), ptL, dphi, dy, wPair);
          else if (isAntiLam && isXiPlus)
            histos.fill(HIST("PairsPt/AntiLam_XiP/hPtDeltaPhiDeltaY"), ptL, dphi, dy, wPair);
        }
      }
    }

    // [P3] Per-event auto-correlation summary.
    histos.fill(HIST("QA/AutoCorr/hXiPairsTotal"), nXiPairsTotal);
    histos.fill(HIST("QA/AutoCorr/hXiPairsVetoed"), nXiPairsVetoed);
    if (nXiPairsVetoed > 0) {
      LOGF(debug,
           "[LXi] Λ-Ξ veto: %d / %d pairs (%.1f%%) dropped as same-Λ in event",
           nXiPairsVetoed, nXiPairsTotal,
           nXiPairsTotal > 0 ? 100.0f * nXiPairsVetoed / nXiPairsTotal : 0.0f);
    }
  }

  // Omega pair loop: same structure as Xi pairs but uses Omega mass for rapidity + kaon PID
  template <bool IsMC = false, typename L, typename C, typename F>
  void analyzeOmegaPairs(L const& lambdas, C const& cascades, F const& flagsStart, float centVal)
  {
    // [P3] Per-event auto-correlation bookkeeping for Λ-Ω.
    int nOmPairsTotal = 0;
    int nOmPairsVetoed = 0;

    for (const auto& lam : lambdas) {
      if (std::abs(lam.rap()) > maxY)
        continue;
      float wLam = useEff ? lam.corrFact() : 1.0f;
      // [P6] Use the typed enum instead of magic literal "1".
      bool isAntiLam = (lam.v0Type() == (int8_t)kAntiLambda);

      for (const auto& casc : cascades) {
        auto fr = flagsStart + casc.globalIndex();
        const int f = fr.isSelected();
        // [Phase 16h] Species filter — only Ω-flagged cascades pair as Ω.
        if (f != kFlagOmegaOnly && f != kFlagXiAndOmega)
          continue;

        // [Phase 12a] Centralised purity gate.
        if (!passesCascadePurityGate<IsMC>(fr))
          continue;

        auto bachTrack = casc.template bachelor_as<FullTracksExtIUWithPID>();
        if (std::abs(bachTrack.tpcNSigmaKa()) > tpcNsigmaBachKaon)
          continue;

        // [P2] Use the cascade row's own Omega rapidity (correct mass).
        float omY = casc.yOmega();
        if (std::abs(omY) > maxY)
          continue;

        ++nOmPairsTotal;

        // [Phase 6] cVetoMode (see processXi pair loop for full doc).
        const bool posMatch = (lam.posTrackId() == casc.posTrackId());
        const bool negMatch = (lam.negTrackId() == casc.negTrackId());
        const bool veto =
            (cVetoMode == lcorr_const::kVetoModeStrict && posMatch && negMatch) ||
            (cVetoMode == lcorr_const::kVetoModeLoose  && (posMatch || negMatch));
        if (veto) {
          ++nOmPairsVetoed;
          continue;
        }

        // [Phase 6] Cascade efficiency weight (Ω flavour).
        float wCascOm = getCascadeEfficiency<1 /*Omega*/>(casc.sign(), casc.pt(), omY);
        float wPair = wLam * wCascOm;

        float dphi = RecoDecay::constrainAngle(casc.phi() - lam.phi(), -PIHalf);
        float dy = omY - lam.rap();

        bool isOmPlus = (casc.sign() > 0);
        bool isLS = (isAntiLam == isOmPlus);

        if (pairCfg.cFillLamOm) {
          if (!isAntiLam && !isOmPlus)
            histos.fill(HIST("Pairs/Lam_OmM/hDeltaPhiDeltaY"), centVal, dphi, dy, wPair);
          else if (!isAntiLam && isOmPlus)
            histos.fill(HIST("Pairs/Lam_OmP/hDeltaPhiDeltaY"), centVal, dphi, dy, wPair);
          else if (isAntiLam && !isOmPlus)
            histos.fill(HIST("Pairs/AntiLam_OmM/hDeltaPhiDeltaY"), centVal, dphi, dy, wPair);
          else if (isAntiLam && isOmPlus)
            histos.fill(HIST("Pairs/AntiLam_OmP/hDeltaPhiDeltaY"), centVal, dphi, dy, wPair);
        }

        // [Phase 12b] Baryon-LS / OS combined fills.
        if (isLS)
          histos.fill(HIST("Pairs/LS_LamOm/hDeltaPhiDeltaY"), centVal, dphi, dy, wPair);
        else
          histos.fill(HIST("Pairs/OS_LamOm/hDeltaPhiDeltaY"), centVal, dphi, dy, wPair);

        // [Phase 12b] pT-differential pair fill (opt-in).
        if (cFillPtDifferentialPairs) {
          float ptL = lam.pt();
          if (!isAntiLam && !isOmPlus)
            histos.fill(HIST("PairsPt/Lam_OmM/hPtDeltaPhiDeltaY"), ptL, dphi, dy, wPair);
          else if (!isAntiLam && isOmPlus)
            histos.fill(HIST("PairsPt/Lam_OmP/hPtDeltaPhiDeltaY"), ptL, dphi, dy, wPair);
          else if (isAntiLam && !isOmPlus)
            histos.fill(HIST("PairsPt/AntiLam_OmM/hPtDeltaPhiDeltaY"), ptL, dphi, dy, wPair);
          else if (isAntiLam && isOmPlus)
            histos.fill(HIST("PairsPt/AntiLam_OmP/hPtDeltaPhiDeltaY"), ptL, dphi, dy, wPair);
        }
      }
    }

    // [P3] Per-event auto-correlation summary.
    histos.fill(HIST("QA/AutoCorr/hOmPairsTotal"), nOmPairsTotal);
    histos.fill(HIST("QA/AutoCorr/hOmPairsVetoed"), nOmPairsVetoed);
    if (nOmPairsVetoed > 0) {
      LOGF(debug,
           "[LXi] Λ-Ω veto: %d / %d pairs (%.1f%%) dropped as same-Λ in event",
           nOmPairsVetoed, nOmPairsTotal,
           nOmPairsTotal > 0 ? 100.0f * nOmPairsVetoed / nOmPairsTotal : 0.0f);
    }
  }

  void processXi(LambdaCollisionsExt::iterator const& lambdacoll,
                 GoodLambdas const& /*lambdas*/,
                 aod::CascDataExt const& cascades,
                 aod::CascadeFlags const& cascflags,
                 FullTracksExtIUWithPID const& /*tracks*/)  // [Phase 14] needed for bachelor_as<>
  {
    // [Phase 16r] Event count is now owned by processYields (always-on by
    // default). Removed from here to avoid N× inflation under DPL's
    // batch-by-process-function execution.

    // [Phase 12a] Trigger-Λ partition choice extracted into pickLambdaPartition().
    auto lambdasInThisEvent = pickLambdaPartition(lambdacoll);
    const int64_t refCollisionIndex = lambdacoll.refCollId();
    auto cascadesInThisEvent = cascades.sliceBy(cascadesPerCollision, refCollisionIndex);

    float pvX = lambdacoll.posX();
    float pvY = lambdacoll.posY();
    float pvZ = lambdacoll.posZ();

    auto flagsStart = cascflags.begin();

    // [Phase 16s] Λ singles fill is owned by processYields (single canonical
    // owner across the workflow — neither processXi nor processMCRecoXi
    // calls it anymore to avoid duplicate fills under any process-switch
    // combination).
    float centVal = lambdacoll.cent();
    analyzeSinglesXi(cascadesInThisEvent, flagsStart, pvX, pvY, pvZ, centVal);
    analyzePairs(lambdasInThisEvent, cascadesInThisEvent, flagsStart, centVal);
    // [Phase 16a → fixed Phase 16g] Per-event coincidence: how many primary-Λ
    // vs how many Ξ candidates share THIS collision. Now species-filtered via
    // countSpeciesEligible<>() so hLamXi reflects Ξ-only multiplicity, not
    // the combined Ξ+Ω cascade slice. Underflow row/col still reveals
    // "events with only one species".
    histos.fill(HIST("Yields/Coincidence/hLamXi"),
                lambdasInThisEvent.size(),
                countSpeciesEligible<false>(cascadesInThisEvent, flagsStart));
  }
  PROCESS_SWITCH(LambdaXiCorrelation, processXi, "Λ–Ξ correlation", true);

  void processOmega(LambdaCollisionsExt::iterator const& lambdacoll,
                    GoodLambdas const& /*lambdas*/,
                    aod::CascDataExt const& cascades,
                    aod::CascadeFlags const& cascflags,
                    FullTracksExtIUWithPID const& /*tracks*/)
  {
    // [Phase 16r] Event count is now owned by processYields (always-on by
    // default). Removed from here to avoid N× inflation under DPL's
    // batch-by-process-function execution.

    auto lambdasInThisEvent = pickLambdaPartition(lambdacoll);
    const int64_t refCollisionIndex = lambdacoll.refCollId();
    auto cascadesInThisEvent = cascades.sliceBy(cascadesPerCollision, refCollisionIndex);

    float pvX = lambdacoll.posX();
    float pvY = lambdacoll.posY();
    float pvZ = lambdacoll.posZ();

    auto flagsStart = cascflags.begin();

    // [Phase 16r] Λ singles fill is owned by processXi (data) / processMCRecoXi (MC).
    // Skipped here to prevent double-fill when both Ξ and Ω process functions are on.
    float centVal = lambdacoll.cent();
    analyzeSinglesOmega(cascadesInThisEvent, flagsStart, pvX, pvY, pvZ, centVal);
    analyzeOmegaPairs(lambdasInThisEvent, cascadesInThisEvent, flagsStart, centVal);
    // [Phase 16a → fixed Phase 16g] Per-event coincidence — Ω-eligible only.
    histos.fill(HIST("Yields/Coincidence/hLamOm"),
                lambdasInThisEvent.size(),
                countSpeciesEligible<true>(cascadesInThisEvent, flagsStart));
  }
  PROCESS_SWITCH(LambdaXiCorrelation, processOmega, "Λ–Ω correlation", false);

  // ---------------------------------------------------------------------------
  // MC Reco-level with truth matching: Λ–Ξ
  // Same as processXi but cascades carry McCascLabels → tree gets pdgCode/isPhysPrim.
  // ---------------------------------------------------------------------------
  void processMCRecoXi(LambdaCollisionsExt::iterator const& lambdacoll,
                       GoodLambdas const& /*lambdas*/,
                       LabeledCascades const& cascades,
                       aod::CascadeFlags const& cascflags,
                       FullTracksExtIUWithPID const& /*tracks*/,  // [Phase 14] for bachelor_as<>
                       aod::McParticles const& /*mcparts*/)
  {
    // [Phase 16r] Event count is now owned by processYields (always-on by
    // default). Removed from here to avoid N× inflation under DPL's
    // batch-by-process-function execution.

    auto lambdasInThisEvent = pickLambdaPartition(lambdacoll);
    const int64_t refCollisionIndex = lambdacoll.refCollId();
    auto cascadesInThisEvent = cascades.sliceBy(labeledCascPerCollision, refCollisionIndex);

    float pvX = lambdacoll.posX();
    float pvY = lambdacoll.posY();
    float pvZ = lambdacoll.posZ();

    auto flagsStart = cascflags.begin();

    // [Phase 16s] Λ singles owned by processYields.
    float centVal = lambdacoll.cent();
    analyzeSinglesXi<true>(cascadesInThisEvent, flagsStart, pvX, pvY, pvZ, centVal);
    analyzePairs<true>(lambdasInThisEvent, cascadesInThisEvent, flagsStart, centVal);
    // [Phase 16a → fixed Phase 16g] Per-event coincidence — Ξ-eligible only.
    histos.fill(HIST("Yields/Coincidence/hLamXi"),
                lambdasInThisEvent.size(),
                countSpeciesEligible<false>(cascadesInThisEvent, flagsStart));
  }
  PROCESS_SWITCH(LambdaXiCorrelation, processMCRecoXi, "MC reco Λ–Ξ (truth-tagged tree)", false);

  // ---------------------------------------------------------------------------
  // MC Reco-level with truth matching: Λ–Ω
  // Same as processOmega but cascades carry McCascLabels → tree gets pdgCode/isPhysPrim.
  // ---------------------------------------------------------------------------
  void processMCRecoOmega(LambdaCollisionsExt::iterator const& lambdacoll,
                          GoodLambdas const& /*lambdas*/,
                          LabeledCascades const& cascades,
                          aod::CascadeFlags const& cascflags,
                          FullTracksExtIUWithPID const& /*tracks*/,
                          aod::McParticles const& /*mcparts*/)
  {
    // [Phase 16r] Event count is now owned by processYields (always-on by
    // default). Removed from here to avoid N× inflation under DPL's
    // batch-by-process-function execution.

    auto lambdasInThisEvent = pickLambdaPartition(lambdacoll);
    const int64_t refCollisionIndex = lambdacoll.refCollId();
    auto cascadesInThisEvent = cascades.sliceBy(labeledCascPerCollision, refCollisionIndex);

    float pvX = lambdacoll.posX();
    float pvY = lambdacoll.posY();
    float pvZ = lambdacoll.posZ();

    auto flagsStart = cascflags.begin();

    // [Phase 16r] Λ singles owned by processMCRecoXi. Skipped to avoid double-fill.
    float centVal = lambdacoll.cent();
    analyzeSinglesOmega<true>(cascadesInThisEvent, flagsStart, pvX, pvY, pvZ, centVal);
    analyzeOmegaPairs<true>(lambdasInThisEvent, cascadesInThisEvent, flagsStart, centVal);
    // [Phase 16a → fixed Phase 16g] Per-event coincidence — Ω-eligible only.
    histos.fill(HIST("Yields/Coincidence/hLamOm"),
                lambdasInThisEvent.size(),
                countSpeciesEligible<true>(cascadesInThisEvent, flagsStart));
  }
  PROCESS_SWITCH(LambdaXiCorrelation, processMCRecoOmega, "MC reco Λ–Ω (truth-tagged tree)", false);

  // ===========================================================================
  // [Phase 13a] Same-species and cross-species pair loops (Λ-Λ, Ξ-Ξ, Ω-Ω, Ξ-Ω)
  // ===========================================================================

  // Λ-Λ pair loop: distinct V0s within the trigger partition. The auto-
  // correlation veto reuses the existing cVetoMode policy via shared
  // posTrackId/negTrackId. Filled into Pairs/Lam_Lam, Lam_AntiLam,
  // AntiLam_AntiLam.
  template <typename L>
  void analyzeLambdaLambdaPairs(L const& lambdas, float centVal)
  {
    if (!pairCfg.cFillLamLam)
      return;
    for (auto i = lambdas.begin(); i != lambdas.end(); ++i) {
      const auto& lam1 = *i;
      if (std::abs(lam1.rap()) > maxY)
        continue;
      bool isAnti1 = (lam1.v0Type() == (int8_t)kAntiLambda);
      float w1 = useEff ? lam1.corrFact() : 1.0f;
      auto j = i;
      for (++j; j != lambdas.end(); ++j) {
        const auto& lam2 = *j;
        if (std::abs(lam2.rap()) > maxY)
          continue;
        // Auto-correlation veto: drop pairs that share daughters.
        const bool posMatch = (lam1.posTrackId() == lam2.posTrackId());
        const bool negMatch = (lam1.negTrackId() == lam2.negTrackId());
        const bool veto = (cVetoMode == lcorr_const::kVetoModeStrict && posMatch && negMatch) ||
                          (cVetoMode == lcorr_const::kVetoModeLoose  && (posMatch || negMatch));
        if (veto)
          continue;

        bool isAnti2 = (lam2.v0Type() == (int8_t)kAntiLambda);
        float w2 = useEff ? lam2.corrFact() : 1.0f;
        float wPair = w1 * w2;
        float dphi = RecoDecay::constrainAngle(lam2.phi() - lam1.phi(), -PIHalf);
        float dy   = lam2.rap() - lam1.rap();

        if (!isAnti1 && !isAnti2)
          histos.fill(HIST("Pairs/Lam_Lam/hDeltaPhiDeltaY"), centVal, dphi, dy, wPair);
        else if (isAnti1 && isAnti2)
          histos.fill(HIST("Pairs/AntiLam_AntiLam/hDeltaPhiDeltaY"), centVal, dphi, dy, wPair);
        else
          histos.fill(HIST("Pairs/Lam_AntiLam/hDeltaPhiDeltaY"), centVal, dphi, dy, wPair);
      }
    }
  }

  // Ξ-Ξ pair loop: distinct cascades passing the Ξ-flag (flag 1 or 2).
  // Sign split into XiM-XiM, XiP-XiP, XiM-XiP. Veto policy: shared bachelor.
  template <bool IsMC = false, typename C, typename F>
  void analyzeXiXiPairs(C const& cascades, F const& flagsStart, float centVal)
  {
    if (!pairCfg.cFillXiXi)
      return;
    for (auto i = cascades.begin(); i != cascades.end(); ++i) {
      const auto& c1 = *i;
      auto fr1 = flagsStart + c1.globalIndex();
      // [Phase 16h] Named flag constants (was: magic 1/2).
      const int f1 = fr1.isSelected();
      if (f1 != kFlagXiOnly && f1 != kFlagXiAndOmega)
        continue;
      if (!passesCascadePurityGate<IsMC>(fr1))
        continue;
      if (std::abs(c1.yXi()) > maxY)
        continue;
      auto j = i;
      for (++j; j != cascades.end(); ++j) {
        const auto& c2 = *j;
        auto fr2 = flagsStart + c2.globalIndex();
        const int f2 = fr2.isSelected();
        if (f2 != kFlagXiOnly && f2 != kFlagXiAndOmega)
          continue;
        if (!passesCascadePurityGate<IsMC>(fr2))
          continue;
        if (std::abs(c2.yXi()) > maxY)
          continue;
        float w1 = getCascadeEfficiency<0>(c1.sign(), c1.pt(), c1.yXi());
        float w2 = getCascadeEfficiency<0>(c2.sign(), c2.pt(), c2.yXi());
        float wPair = w1 * w2;
        float dphi = RecoDecay::constrainAngle(c2.phi() - c1.phi(), -PIHalf);
        float dy   = c2.yXi() - c1.yXi();
        bool m1 = (c1.sign() < 0), m2 = (c2.sign() < 0);
        if (m1 && m2)
          histos.fill(HIST("Pairs/XiM_XiM/hDeltaPhiDeltaY"), centVal, dphi, dy, wPair);
        else if (!m1 && !m2)
          histos.fill(HIST("Pairs/XiP_XiP/hDeltaPhiDeltaY"), centVal, dphi, dy, wPair);
        else
          histos.fill(HIST("Pairs/XiM_XiP/hDeltaPhiDeltaY"), centVal, dphi, dy, wPair);
      }
    }
  }

  // Ω-Ω pair loop. Same as Ξ-Ξ but on flag 2/3 (Ω-eligible) and uses yOmega.
  // Bachelor must satisfy the kaon-PID cut, which is enforced inside the
  // cascade producer's processCandidate (flag=2 or 3 means Ω passed PID).
  template <bool IsMC = false, typename C, typename F>
  void analyzeOmegaOmegaPairs(C const& cascades, F const& flagsStart, float centVal)
  {
    if (!pairCfg.cFillOmOm)
      return;
    for (auto i = cascades.begin(); i != cascades.end(); ++i) {
      const auto& c1 = *i;
      auto fr1 = flagsStart + c1.globalIndex();
      int f1 = fr1.isSelected();
      // [Phase 16h] Named flag constants (was: magic 2/3).
      if (f1 != kFlagXiAndOmega && f1 != kFlagOmegaOnly)
        continue;
      if (!passesCascadePurityGate<IsMC>(fr1))
        continue;
      if (std::abs(c1.yOmega()) > maxY)
        continue;
      auto j = i;
      for (++j; j != cascades.end(); ++j) {
        const auto& c2 = *j;
        auto fr2 = flagsStart + c2.globalIndex();
        int f2 = fr2.isSelected();
        if (f2 != kFlagXiAndOmega && f2 != kFlagOmegaOnly)
          continue;
        if (!passesCascadePurityGate<IsMC>(fr2))
          continue;
        if (std::abs(c2.yOmega()) > maxY)
          continue;
        float w1 = getCascadeEfficiency<1>(c1.sign(), c1.pt(), c1.yOmega());
        float w2 = getCascadeEfficiency<1>(c2.sign(), c2.pt(), c2.yOmega());
        float wPair = w1 * w2;
        float dphi = RecoDecay::constrainAngle(c2.phi() - c1.phi(), -PIHalf);
        float dy   = c2.yOmega() - c1.yOmega();
        bool m1 = (c1.sign() < 0), m2 = (c2.sign() < 0);
        if (m1 && m2)
          histos.fill(HIST("Pairs/OmM_OmM/hDeltaPhiDeltaY"), centVal, dphi, dy, wPair);
        else if (!m1 && !m2)
          histos.fill(HIST("Pairs/OmP_OmP/hDeltaPhiDeltaY"), centVal, dphi, dy, wPair);
        else
          histos.fill(HIST("Pairs/OmM_OmP/hDeltaPhiDeltaY"), centVal, dphi, dy, wPair);
      }
    }
  }

  // Ξ-Ω cross-species pair loop. flag 1 → Ξ-only, flag 3 → Ω-only,
  // flag 2 → both (we treat as either). Generates 4 sign combos.
  template <bool IsMC = false, typename C, typename F>
  void analyzeXiOmegaPairs(C const& cascades, F const& flagsStart, float centVal)
  {
    if (!pairCfg.cFillXiOm)
      return;
    for (const auto& c1 : cascades) {
      auto fr1 = flagsStart + c1.globalIndex();
      const int f1 = fr1.isSelected();
      if (f1 == kFlagRejected)
        continue;
      if (!passesCascadePurityGate<IsMC>(fr1))
        continue;
      // [Phase 16h] Treat c1 as Ξ if it passes Ξ-eligibility AND |yXi|<maxY.
      if ((f1 != kFlagXiOnly && f1 != kFlagXiAndOmega) || std::abs(c1.yXi()) > maxY)
        continue;
      for (const auto& c2 : cascades) {
        if (c1.globalIndex() == c2.globalIndex())
          continue;
        auto fr2 = flagsStart + c2.globalIndex();
        const int f2 = fr2.isSelected();
        if (!passesCascadePurityGate<IsMC>(fr2))
          continue;
        // [Phase 16h] Treat c2 as Ω if it passes Ω-eligibility AND |yΩ|<maxY.
        if ((f2 != kFlagXiAndOmega && f2 != kFlagOmegaOnly) || std::abs(c2.yOmega()) > maxY)
          continue;
        float w1 = getCascadeEfficiency<0>(c1.sign(), c1.pt(), c1.yXi());
        float w2 = getCascadeEfficiency<1>(c2.sign(), c2.pt(), c2.yOmega());
        float wPair = w1 * w2;
        float dphi = RecoDecay::constrainAngle(c2.phi() - c1.phi(), -PIHalf);
        float dy   = c2.yOmega() - c1.yXi();
        bool xiM = (c1.sign() < 0), omM = (c2.sign() < 0);
        if (xiM && omM)
          histos.fill(HIST("Pairs/XiM_OmM/hDeltaPhiDeltaY"), centVal, dphi, dy, wPair);
        else if (xiM && !omM)
          histos.fill(HIST("Pairs/XiM_OmP/hDeltaPhiDeltaY"), centVal, dphi, dy, wPair);
        else if (!xiM && omM)
          histos.fill(HIST("Pairs/XiP_OmM/hDeltaPhiDeltaY"), centVal, dphi, dy, wPair);
        else
          histos.fill(HIST("Pairs/XiP_OmP/hDeltaPhiDeltaY"), centVal, dphi, dy, wPair);
      }
    }
  }

  // ---------------------------------------------------------------------------
  // [Phase 13a] Pair driver. Single process function that runs ALL the
  // same/cross-species pair loops on demand. Use this WITHOUT processXi /
  // processOmega when you want full coverage; use the older Λ-X functions
  // alone if you only need the Λ-trigger correlations.
  // ---------------------------------------------------------------------------
  void processAllPairs(LambdaCollisionsExt::iterator const& lambdacoll,
                       GoodLambdas const& /*lambdas*/,
                       aod::CascDataExt const& cascades,
                       aod::CascadeFlags const& cascflags,
                       FullTracksExtIUWithPID const& /*tracks*/)
  {
    // [Phase 16r] Event count is now owned by processYields (always-on by
    // default). Removed from here to avoid N× inflation under DPL's
    // batch-by-process-function execution.
    auto lambdasInThisEvent = pickLambdaPartition(lambdacoll);
    const int64_t refCollisionIndex = lambdacoll.refCollId();
    auto cascadesInThisEvent = cascades.sliceBy(cascadesPerCollision, refCollisionIndex);
    float centVal = lambdacoll.cent();
    auto flagsStart = cascflags.begin();
    analyzeLambdaLambdaPairs(lambdasInThisEvent, centVal);
    analyzeXiXiPairs(cascadesInThisEvent, flagsStart, centVal);
    analyzeOmegaOmegaPairs(cascadesInThisEvent, flagsStart, centVal);
    analyzeXiOmegaPairs(cascadesInThisEvent, flagsStart, centVal);
  }
  PROCESS_SWITCH(LambdaXiCorrelation, processAllPairs, "All same-/cross-species pair loops (Λ-Λ, Ξ-Ξ, Ω-Ω, Ξ-Ω)", false);

  // ---------------------------------------------------------------------------
  // [Phase 13b] processYields — fills the per-event multiplicity + ⟨pT⟩
  // histograms ONCE per event for every species, regardless of which pair-
  // process functions are also enabled. Default ON.
  // ---------------------------------------------------------------------------
  void processYields(LambdaCollisionsExt::iterator const& lambdacoll,
                     GoodLambdas const& /*lambdas*/,
                     aod::CascDataExt const& cascades,
                     aod::CascadeFlags const& cascflags)
  {
    // [Phase 16r] processYields is the canonical owner of Event/hEventCount.
    // Default-on, fires exactly once per LambdaCollision iterator → unique
    // event count regardless of which other process functions are active.
    fillEventCountOnce(lambdacoll.globalIndex());
    auto lambdasInThisEvent = pickLambdaPartition(lambdacoll);
    // [Phase 16s] processYields is ALSO the canonical owner of the Λ
    // singles histograms + LambdaCandidates tree fill. Previously these
    // were in processXi / processMCRecoXi but DPL appears to fire them
    // double under some process-switch combinations.
    analyzeSinglesLambda(lambdasInThisEvent, lambdacoll);
    if (!yieldCfg.cFillEventYields)
      return;
    const int64_t refCollisionIndex = lambdacoll.refCollId();
    auto cascadesInThisEvent = cascades.sliceBy(cascadesPerCollision, refCollisionIndex);
    auto flagsStart = cascflags.begin();
    float centVal = lambdacoll.cent();

    // Per-event accumulators.
    int    nLam = 0, nALam = 0, nXiM = 0, nXiP = 0, nOmM = 0, nOmP = 0;
    double sLam = 0, sALam = 0, sXiM = 0, sXiP = 0, sOmM = 0, sOmP = 0;

    for (const auto& l : lambdasInThisEvent) {
      if (std::abs(l.rap()) > maxY)
        continue;
      if (l.v0Type() == (int8_t)kLambda) { ++nLam;  sLam  += l.pt(); }
      else                                { ++nALam; sALam += l.pt(); }
    }
    for (const auto& c : cascadesInThisEvent) {
      auto fr = flagsStart + c.globalIndex();
      int f = fr.isSelected();
      if (f == 0)
        continue;
      // Ξ-eligible (flag 1 or 2)
      if ((f == kFlagXiOnly || f == kFlagXiAndOmega) && std::abs(c.yXi()) <= maxY) {
        if (c.sign() < 0) { ++nXiM; sXiM += c.pt(); }
        else              { ++nXiP; sXiP += c.pt(); }
      }
      // Ω-eligible (flag 2 or 3)
      if ((f == kFlagXiAndOmega || f == kFlagOmegaOnly) && std::abs(c.yOmega()) <= maxY) {
        if (c.sign() < 0) { ++nOmM; sOmM += c.pt(); }
        else              { ++nOmP; sOmP += c.pt(); }
      }
    }

    auto fillSpecies = [&](const auto& nKey, const auto& mptKey,
                           const auto& mNvsCKey, const auto& mPtvsCKey,
                           const auto& nptKey,
                           int n, double sumPt) {
      histos.fill(nKey, n);
      histos.fill(mNvsCKey, centVal, n);
      if (n > 0) {
        double mpt = sumPt / n;
        histos.fill(mptKey, mpt);
        histos.fill(mPtvsCKey, centVal, mpt);
        histos.fill(nptKey, n, mpt);
      }
    };

    fillSpecies(HIST("Yields/Lambda/hNPerEvent"),
                HIST("Yields/Lambda/hMeanPtPerEvent"),
                HIST("Yields/Lambda/hMeanNvsCent"),
                HIST("Yields/Lambda/hMeanPtVsCent"),
                HIST("Yields/Lambda/hNvsPt2D"),
                nLam, sLam);
    fillSpecies(HIST("Yields/AntiLambda/hNPerEvent"),
                HIST("Yields/AntiLambda/hMeanPtPerEvent"),
                HIST("Yields/AntiLambda/hMeanNvsCent"),
                HIST("Yields/AntiLambda/hMeanPtVsCent"),
                HIST("Yields/AntiLambda/hNvsPt2D"),
                nALam, sALam);
    fillSpecies(HIST("Yields/XiMinus/hNPerEvent"),
                HIST("Yields/XiMinus/hMeanPtPerEvent"),
                HIST("Yields/XiMinus/hMeanNvsCent"),
                HIST("Yields/XiMinus/hMeanPtVsCent"),
                HIST("Yields/XiMinus/hNvsPt2D"),
                nXiM, sXiM);
    fillSpecies(HIST("Yields/XiPlus/hNPerEvent"),
                HIST("Yields/XiPlus/hMeanPtPerEvent"),
                HIST("Yields/XiPlus/hMeanNvsCent"),
                HIST("Yields/XiPlus/hMeanPtVsCent"),
                HIST("Yields/XiPlus/hNvsPt2D"),
                nXiP, sXiP);
    fillSpecies(HIST("Yields/OmegaMinus/hNPerEvent"),
                HIST("Yields/OmegaMinus/hMeanPtPerEvent"),
                HIST("Yields/OmegaMinus/hMeanNvsCent"),
                HIST("Yields/OmegaMinus/hMeanPtVsCent"),
                HIST("Yields/OmegaMinus/hNvsPt2D"),
                nOmM, sOmM);
    fillSpecies(HIST("Yields/OmegaPlus/hNPerEvent"),
                HIST("Yields/OmegaPlus/hMeanPtPerEvent"),
                HIST("Yields/OmegaPlus/hMeanNvsCent"),
                HIST("Yields/OmegaPlus/hMeanPtVsCent"),
                HIST("Yields/OmegaPlus/hNvsPt2D"),
                nOmP, sOmP);
  }
  PROCESS_SWITCH(LambdaXiCorrelation, processYields, "Per-event multiplicity + ⟨pT⟩ for each species", true);

  // ---------------------------------------------------------------------------
  // MC Gen-level: Λ–Ξ truth correlation (closure test)
  // Loops over generator-level primary Λ from LambdaMcGenTracks and
  // generator-level Ξ from McParticles (inline PDG == ±3312 check).
  // McParticles are sliced per McCollision using the stored refMcCollId.
  // No reconstruction, no PID cuts — pure truth-level R2.
  // ---------------------------------------------------------------------------
  // [Phase 5] Single-pass MC gen Λ-Ξ: classify Lambdas and Xis once, then
  // do a single pair loop over the classified vectors. The previous version
  // scanned `genLambdasAll` four times (Lambda singles, AntiLambda singles,
  // Lam-pairs, AntiLam-pairs) and `genXis` three times. With this rewrite
  // each table is scanned exactly once.
  void processMCGenXi(aod::LambdaMcGenCollisions::iterator const& mcgencol,
                      aod::LambdaMcGenTracks const& genLambdasAll,
                      aod::McParticles const& allMcParts)
  {
    int32_t thisColId = mcgencol.globalIndex();
    // [Phase 16n] gen process functions do NOT fill Event/hEventCount —
    // that's the reco counter's job. seenGenCollOnce is reserved for any
    // future gen-side idempotence needs.
    (void)seenGenCollOnce(static_cast<int64_t>(thisColId));
    // [Phase 16t] Gen event counter (R₂ closure-test denominator). Owned
    // by processMCGenXi (canonical gen "first owner").
    histos.fill(HIST("McGen/Event/hEventCount"), 0.5);
    float centVal = mcgencol.cent();

    int64_t mcCollId = mcgencol.refMcCollId();
    auto genXis = allMcParts.sliceBy(mcParticlesPerMcCollision, mcCollId);

    // [Phase 16r] processMCGenXi owns the gen-Λ singles + tree fill.
    // processMCGenOmega skips Λ fills (only builds local vectors). This is
    // the ownership model — a per-event guard would not work because DPL
    // batches by process function, not by event.

    // [Phase 5] Pre-classified caches. Reserved sizes are heuristic for typical
    // Pb-Pb central event populations.
    struct LamLite { float pt, rap, phi; };
    struct XiLite  { float pt, rap, phi; bool isPlus; };
    std::vector<LamLite> goodLam, goodAntiLam;
    std::vector<XiLite>  goodXi;
    goodLam.reserve(32);
    goodAntiLam.reserve(32);
    goodXi.reserve(16);

    // Single pass over Lambdas: fill singles, classify into Lam/AntiLam.
    for (const auto& lam : genLambdasAll) {
      if (lam.lambdaMcGenCollisionId() != thisColId)
        continue;
      if (lam.v0PrmScd() != (int8_t)kPrimary)
        continue;
      if (std::abs(lam.rap()) > maxY)
        continue;
      LamLite l{static_cast<float>(lam.pt()), static_cast<float>(lam.rap()), static_cast<float>(lam.phi())};
      // [Phase 16t] Wrap φ to [0, 2π) to match the reco-side convention.
      float phiWrapped = RecoDecay::constrainAngle(lam.phi(), 0.f);
      if (lam.v0Type() == (int8_t)kLambda) {
        histos.fill(HIST("McGen/Singles/Lambda/hPt"), lam.pt());
        histos.fill(HIST("McGen/Singles/Lambda/hYPhi"), lam.rap(), phiWrapped);
        goodLam.push_back(l);
      } else if (lam.v0Type() == (int8_t)kAntiLambda) {
        histos.fill(HIST("McGen/Singles/AntiLambda/hPt"), lam.pt());
        histos.fill(HIST("McGen/Singles/AntiLambda/hYPhi"), lam.rap(), phiWrapped);
        goodAntiLam.push_back(l);
      }

      // [Phase 16a → Phase 16r] Gen-Λ TTree fill. Owned by processMCGenXi.
      if (bp && saveLambdaTree) {
        auto& g = bp->lamGen;
        g.pt       = lam.pt();
        g.eta      = lam.eta();
        g.rap      = lam.rap();
        g.phi      = lam.phi();
        g.v0Type   = static_cast<int>(lam.v0Type());
        g.v0PrmScd = static_cast<int>(lam.v0PrmScd());
        g.cent     = centVal;
        g.pvZ      = mcgencol.posZ();
        treeLambdaGen->Fill();
      }
    }

    // Single pass over Xis: fill singles, eff denominator, optional tree, classify.
    for (const auto& xi : genXis) {
      if (std::abs(xi.pdgCode()) != lcorr_const::kXiMinusPdg)
        continue;
      if (!xi.isPhysicalPrimary())
        continue;
      float xiY = RecoDecay::y(std::array{xi.px(), xi.py(), xi.pz()}, MassXiMinus);
      if (std::abs(xiY) > maxY)
        continue;
      bool isPlus = (xi.pdgCode() == -lcorr_const::kXiMinusPdg);
      // [Phase 16t] φ for closure-test (y, φ) singles.
      float xiPhiWrapped = RecoDecay::constrainAngle(xi.phi(), 0.f);
      if (!isPlus) {
        histos.fill(HIST("McGen/Singles/XiMinus/hPtVsRap"), xiY, xi.pt());
        histos.fill(HIST("McGen/Singles/XiMinus/hYPhi"), xiY, xiPhiWrapped);
        histos.fill(HIST("Eff/Gen/XiMinus/hCentPtRap"), centVal, xi.pt(), xiY);
      } else {
        histos.fill(HIST("McGen/Singles/XiPlus/hPtVsRap"), xiY, xi.pt());
        histos.fill(HIST("McGen/Singles/XiPlus/hYPhi"), xiY, xiPhiWrapped);
        histos.fill(HIST("Eff/Gen/XiPlus/hCentPtRap"), centVal, xi.pt(), xiY);
      }
      // [Phase 8 fix] treeXiGen is only setObject-ed when saveCascTree=true.
      if (bp && saveCascTree) {
        auto& b = bp->xiGen;
        b.pt = xi.pt();
        b.rap = xiY;
        b.pdgCode = xi.pdgCode();
        b.isPhysPrim = xi.isPhysicalPrimary();
        b.cent = centVal;
        b.pvZ = mcgencol.posZ();
        treeXiGen->Fill();
      }
      goodXi.push_back({static_cast<float>(xi.pt()), xiY, static_cast<float>(xi.phi()), isPlus});
    }

    // Pair loop over the classified caches. One pass over Lam × Xi, one over AntiLam × Xi.
    for (const auto& l : goodLam) {
      for (const auto& x : goodXi) {
        float dphi = RecoDecay::constrainAngle(x.phi - l.phi, -PIHalf);
        float dy = x.rap - l.rap;
        if (!x.isPlus)
          histos.fill(HIST("McGen/Pairs/Lam_XiM/hDeltaPhiDeltaY"), centVal, dphi, dy);
        else
          histos.fill(HIST("McGen/Pairs/Lam_XiP/hDeltaPhiDeltaY"), centVal, dphi, dy);
      }
    }
    for (const auto& l : goodAntiLam) {
      for (const auto& x : goodXi) {
        float dphi = RecoDecay::constrainAngle(x.phi - l.phi, -PIHalf);
        float dy = x.rap - l.rap;
        if (!x.isPlus)
          histos.fill(HIST("McGen/Pairs/AntiLam_XiM/hDeltaPhiDeltaY"), centVal, dphi, dy);
        else
          histos.fill(HIST("McGen/Pairs/AntiLam_XiP/hDeltaPhiDeltaY"), centVal, dphi, dy);
      }
    }
  }
  PROCESS_SWITCH(LambdaXiCorrelation, processMCGenXi, "MC gen-level Λ–Ξ closure", false);

  // ---------------------------------------------------------------------------
  // MC Gen-level: Λ–Ω truth correlation (closure test)
  // Same structure as processMCGenXi but uses |PDG| == 3334 (Omega).
  // No kaon PID cut needed — truth level.
  // ---------------------------------------------------------------------------
  // [Phase 5] Single-pass MC gen Λ-Ω, structured identically to processMCGenXi.
  void processMCGenOmega(aod::LambdaMcGenCollisions::iterator const& mcgencol,
                         aod::LambdaMcGenTracks const& genLambdasAll,
                         aod::McParticles const& allMcParts)
  {
    int32_t thisColId = mcgencol.globalIndex();
    // [Phase 16n] gen process functions do NOT fill Event/hEventCount —
    // that's the reco counter's job. seenGenCollOnce is reserved for any
    // future gen-side idempotence needs.
    (void)seenGenCollOnce(static_cast<int64_t>(thisColId));
    float centVal = mcgencol.cent();

    int64_t mcCollId = mcgencol.refMcCollId();
    auto genOmegas = allMcParts.sliceBy(mcParticlesPerMcCollision, mcCollId);

    // [Phase 16r] Gen-Λ singles + tree are owned by processMCGenXi (DPL
    // batches by process function, so a per-event guard would fail). Here
    // we only need to BUILD the goodLam/goodAntiLam vectors for the local
    // Λ-Ω pair loop — no histogram or tree side effects.
    struct LamLite { float pt, rap, phi; };
    struct OmLite  { float pt, rap, phi; bool isPlus; };
    std::vector<LamLite> goodLam, goodAntiLam;
    std::vector<OmLite>  goodOm;
    goodLam.reserve(32);
    goodAntiLam.reserve(32);
    goodOm.reserve(8);

    // Single pass over Lambdas — vectors only, no fills.
    for (const auto& lam : genLambdasAll) {
      if (lam.lambdaMcGenCollisionId() != thisColId)
        continue;
      if (lam.v0PrmScd() != (int8_t)kPrimary)
        continue;
      if (std::abs(lam.rap()) > maxY)
        continue;
      LamLite l{static_cast<float>(lam.pt()), static_cast<float>(lam.rap()), static_cast<float>(lam.phi())};
      if (lam.v0Type() == (int8_t)kLambda) {
        goodLam.push_back(l);
      } else if (lam.v0Type() == (int8_t)kAntiLambda) {
        goodAntiLam.push_back(l);
      }
    }

    // Single pass over Omegas.
    for (const auto& om : genOmegas) {
      if (std::abs(om.pdgCode()) != lcorr_const::kOmegaMinusPdg)
        continue;
      if (!om.isPhysicalPrimary())
        continue;
      float omY = RecoDecay::y(std::array{om.px(), om.py(), om.pz()}, MassOmegaMinus);
      if (std::abs(omY) > maxY)
        continue;
      bool isPlus = (om.pdgCode() == -lcorr_const::kOmegaMinusPdg);
      // [Phase 16t] φ for closure-test (y, φ) singles.
      float omPhiWrapped = RecoDecay::constrainAngle(om.phi(), 0.f);
      if (!isPlus) {
        histos.fill(HIST("McGen/Singles/OmegaMinus/hPtVsRap"), omY, om.pt());
        histos.fill(HIST("McGen/Singles/OmegaMinus/hYPhi"), omY, omPhiWrapped);
        histos.fill(HIST("Eff/Gen/OmegaMinus/hCentPtRap"), centVal, om.pt(), omY);
      } else {
        histos.fill(HIST("McGen/Singles/OmegaPlus/hPtVsRap"), omY, om.pt());
        histos.fill(HIST("McGen/Singles/OmegaPlus/hYPhi"), omY, omPhiWrapped);
        histos.fill(HIST("Eff/Gen/OmegaPlus/hCentPtRap"), centVal, om.pt(), omY);
      }
      // [Phase 8 fix] treeOmegaGen is only setObject-ed when saveCascTree=true.
      if (bp && saveCascTree) {
        auto& b = bp->omGen;
        b.pt = om.pt();
        b.rap = omY;
        b.pdgCode = om.pdgCode();
        b.isPhysPrim = om.isPhysicalPrimary();
        b.cent = centVal;
        b.pvZ = mcgencol.posZ();
        treeOmegaGen->Fill();
      }
      goodOm.push_back({static_cast<float>(om.pt()), omY, static_cast<float>(om.phi()), isPlus});
    }

    // Pair loop over the classified caches.
    for (const auto& l : goodLam) {
      for (const auto& o : goodOm) {
        float dphi = RecoDecay::constrainAngle(o.phi - l.phi, -PIHalf);
        float dy = o.rap - l.rap;
        if (!o.isPlus)
          histos.fill(HIST("McGen/Pairs/Lam_OmM/hDeltaPhiDeltaY"), centVal, dphi, dy);
        else
          histos.fill(HIST("McGen/Pairs/Lam_OmP/hDeltaPhiDeltaY"), centVal, dphi, dy);
      }
    }
    for (const auto& l : goodAntiLam) {
      for (const auto& o : goodOm) {
        float dphi = RecoDecay::constrainAngle(o.phi - l.phi, -PIHalf);
        float dy = o.rap - l.rap;
        if (!o.isPlus)
          histos.fill(HIST("McGen/Pairs/AntiLam_OmM/hDeltaPhiDeltaY"), centVal, dphi, dy);
        else
          histos.fill(HIST("McGen/Pairs/AntiLam_OmP/hDeltaPhiDeltaY"), centVal, dphi, dy);
      }
    }
  }
  PROCESS_SWITCH(LambdaXiCorrelation, processMCGenOmega, "MC gen-level Λ–Ω closure", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{

    // [Phase 4] LambdaCascadeProducer absorbs the former LambdaTableProducer
    // and CascadeSelector — single event selection, single source of truth.
    // [Phase 5] LambdaR2Correlation, CascadeSelector (#if 0 block), and
    // CascadeCorrelations were removed entirely from this file.
    adaptAnalysisTask<LambdaCascadeProducer>(cfgc),
    adaptAnalysisTask<LambdaTracksExtProducer>(cfgc),
    adaptAnalysisTask<LambdaXiCorrelation>(cfgc)

  };
}
