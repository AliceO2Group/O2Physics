// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
//

/// \file photonconversionbuilder.ch
/// \brief this task produces photon data table with KFParticle.
/// \author Daiki Sekihata <daiki.sekihata@cern.ch>, Tokyo
#ifndef PWGEM_PHOTONMESON_TABLEPRODUCER_PHOTONCONVERSIONBUILDER_H_
#define PWGEM_PHOTONMESON_TABLEPRODUCER_PHOTONCONVERSIONBUILDER_H_

#ifndef HomogeneousField
#define HomogeneousField // o2-linter: disable=name/macro (name coming from KFParticle, not us, needed for KFParticle::SetField)
#endif

#include "PWGEM/Dilepton/Utils/PairUtilities.h"
#include "PWGEM/PhotonMeson/Core/EmMlResponsePCM.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCandidate.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/DataModel/EventTables.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/ConversionPointProfile.h"
#include "PWGEM/PhotonMeson/Utils/PCMUtilities.h"
#include "PWGEM/PhotonMeson/Utils/TrackSelection.h"
#include "PWGEM/PhotonMeson/Utils/V0DedupUtilities.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TPCVDriftManager.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Tools/KFparticle/KFUtilities.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/LHCConstants.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DataFormatsTPC/VDriftCorrFact.h>
#include <DetectorsBase/GeometryManager.h>
#include <DetectorsBase/MatLayerCylSet.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Concepts.h>
#include <Framework/Configurable.h>
#include <Framework/DataTypes.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>
#include <MathUtils/Primitive2D.h>
#include <ReconstructionDataFormats/HelixHelper.h>
#include <ReconstructionDataFormats/PID.h>
#include <ReconstructionDataFormats/Track.h>

#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TPDGCode.h>

#include <KFPTrack.h>
#include <KFPVertex.h>
#include <KFParticle.h>

#include <GPUROOTCartesianFwd.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <ratio>
#include <set>
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

// NOLINTBEGIN(google-global-names-in-headers)
using o2::framework::AxisSpec;
using o2::framework::Configurable;
using o2::framework::ConfigurableGroup;
using o2::framework::HistType;
using o2::framework::kTH1D;
using o2::framework::kTH1F;
using o2::framework::kTH2F;
using o2::framework::kTH3F;
using o2::framework::kTHnSparseF;
using o2::framework::Preslice;
using o2::framework::Produces;
using o2::pwgem::photonmeson::cpRadiusProfile;
using o2::pwgem::photonmeson::cpScanLeg;
using o2::pwgem::photonmeson::FinderCollision;
using o2::pwgem::photonmeson::FinderCuts;
using o2::pwgem::photonmeson::FinderLeg;
using o2::pwgem::photonmeson::FinderStep;
using o2::pwgem::photonmeson::LegScan;
using o2::pwgem::photonmeson::qInv;
using o2::pwgem::photonmeson::RadiusProfile;
using o2::pwgem::photonmeson::solveConflictGroup;
// NOLINTEND(google-global-names-in-headers)

using MyCollisions = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels, o2::aod::PMEvSels, o2::aod::CentFT0Ms, o2::aod::CentFT0As, o2::aod::CentFT0Cs>;
// using MyCollisionsWithSWT = soa::Join<MyCollisions, aod::EMSWTriggerBitsTMP>;
using MyCollisionsMC = o2::soa::Join<MyCollisions, o2::aod::McCollisionLabels>;
using MyBCs = o2::soa::Join<o2::aod::BCsWithTimestamps, o2::aod::PcmObjects>;

using MyTracksIU = o2::soa::Join<o2::aod::TracksIU, o2::aod::TracksExtra, o2::aod::TracksCovIU, o2::aod::pidTPCFullEl, o2::aod::pidTPCFullPi>;
using MyTracksIUMC = o2::soa::Join<MyTracksIU, o2::aod::McTrackLabels, o2::aod::mcTPCTuneOnData>;

enum MatCorrTypeEM {
  None = 0,
  TGeo = 1,
  LUT = 2
};

enum TrackPropMode {
  kProper = 0,
  kFast = 1,
  kBoth = 2
};

enum V0DeduplicationMode {
  Pairwise = 0,       // old algorithm, V0 with best PCA + cosPA at the same time in a group wins
  GreedyMatching = 1, // new algorithm, best score wins (score combines PCA and cosPA -- see detailed implementation in PCMUtilities.cxx)
  GroupMatching = 2,  // new algorithm, best combination of leg-disjoint V0s wins, not the single best V0
  KeepAll = 3         // no deduplication, keep all V0s
};

// Struct needed for deduplication mode 1. Used to keep track  of v0 properties and a score based on pca and cosPA
struct V0CandidateHelper {
  int v0ID = -1;
  int colID = -1;
  int posID = -1;
  int eleID = -1;
  float cosPA = -1.f;
  float pca = -1.f;
  float score = -1.f;
  float mee = 0.f; // e+e- mass at the secondary vertex (GeV/c^2)
  float pxAtPV = 0.f;
  float pyAtPV = 0.f;
  float pzAtPV = 0.f;
  float svX = 0.f;
  float svY = 0.f;
  float svZ = 0.f;
  float pxAtSV = 0.f; // momentum at the conversion point without the PV constraint (pointing residuals)
  float pyAtSV = 0.f;
  float pzAtSV = 0.f;
  float purity = -999.f;
  float cosPAxy = -2.f; // transverse pointing angle
  float cosPA3D = -2.f; // 3D pointing angle

  V0CandidateHelper() = default;

  V0CandidateHelper(int v0, int col, int pos, int ele, float cpa, float p, float s, float m = 0.f)
    : v0ID(v0), colID(col), posID(pos), eleID(ele), cosPA(cpa), pca(p), score(s), mee(m)
  {
  }

  bool operator==(const V0CandidateHelper& other) const { return score == other.score; }
  bool operator<(const V0CandidateHelper& other) const { return score < other.score; }
  bool operator>(const V0CandidateHelper& other) const { return score > other.score; }
};

// (v0.globalIndex(), collision.globalIndex(), pos.globalIndex(), ele.globalIndex())
// Identifies a V0 candidate by  collision, positive-track and negative track IDs, and the SVertexer V0 index (or -1 for finder and resolver photons)
struct CandKey {
  int64_t v0Id = -1; // SVertexer V0 index, -1 for finder and resolver photons
  int64_t collisionId = -1;
  int64_t posId = -1;
  int64_t eleId = -1;
  auto operator<=>(const CandKey&) const = default;
};

// Stores whether a trial fit succeded and the score of the candidate (for the ambiguity resolver, the score is used to decide which of two crossed candidates to keep)
struct PairFit {
  bool ok = false;
  float score = 999.f;
};

// options for fitting a leg pair that was not built by the S-Vertexer (finder, ambiguity resolver)
struct SourcePV {
  float x = 0.f, y = 0.f, z = 0.f, cXX = 0.f, cXY = 0.f, cXZ = 0.f, cYY = 0.f, cYZ = 0.f, cZZ = 0.f, c2 = 0.f;
  uint16_t nContrib = 0;
  [[nodiscard]] float posX() const { return x; }
  [[nodiscard]] float posY() const { return y; }
  [[nodiscard]] float posZ() const { return z; }
  [[nodiscard]] float covXX() const { return cXX; }
  [[nodiscard]] float covXY() const { return cXY; }
  [[nodiscard]] float covXZ() const { return cXZ; }
  [[nodiscard]] float covYY() const { return cYY; }
  [[nodiscard]] float covYZ() const { return cYZ; }
  [[nodiscard]] float covZZ() const { return cZZ; }
  [[nodiscard]] float chi2() const { return c2; }
  [[nodiscard]] uint16_t numContrib() const { return nContrib; }
};

struct LegDcaKey {
  int64_t trackId = -1, collisionId = -1;
  float x = 0.f, y = 0.f, z = 0.f;
  bool sideBased = false;
  bool operator==(LegDcaKey const& o) const { return trackId == o.trackId && collisionId == o.collisionId && x == o.x && y == o.y && z == o.z && sideBased == o.sideBased; }
};
struct LegDcaKeyHash {
  size_t operator()(LegDcaKey const& k) const
  {
    size_t h = std::hash<int64_t>{}(k.trackId) * 1000003u ^ std::hash<int64_t>{}(k.collisionId);
    h = h * 1000003u ^ std::hash<float>{}(k.x);
    h = h * 1000003u ^ std::hash<float>{}(k.y);
    return (h * 1000003u ^ std::hash<float>{}(k.z)) * 2u + (k.sideBased ? 1u : 0u);
  }
};
// Caches a propagated leg (electron PID) at its DCA to the vertex, for reuse in the ambiguity resolver and the mixed-event builder
struct LegAtDca {
  std::array<float, 2> dca{};
  o2::track::TrackParCov trk; // leg (electron PID) propagated to its DCA to the vertex
};

struct FitOptions {
  bool nonSVertexer = false;    // finder PCA cut for mixed legs, side-based TPC drift correction
  bool relaxAP = false;         // no Armenteros-Podolanski cut
  bool purityGate = false;      // finder: require purity score > finderMinScore
  float minPurity = -999.f;     // purity gate threshold for this fit; below -998 the configured finderMinScore
  PairFit* probe = nullptr;     // only return the dedup score, do not add a candidate
  const SourcePV* pv = nullptr; // run 2: fit against this vertex instead of the collision's (original pairings at their own source vertex)
};

// twin of PCMFreeLegs for the mixed instance of the builder
namespace o2::aod
{
DECLARE_SOA_TABLE(PCMFreeLegsMix, "AOD", "PCMFREELEGMIX", //!
                  o2::soa::Index<>, pcmfreeleg::CollisionId, pcmfreeleg::TrackId, v0leg::Sign,
                  v0legcov::X, v0legcov::Alpha, v0legcov::Y, v0legcov::Z, v0legcov::Snp, v0legcov::Tgl, v0legcov::Signed1Pt, v0legcov::Cov,
                  track::TPCNClsFindable, track::TPCNClsFindableMinusFound, track::TPCNClsFindableMinusCrossedRows, track::TPCNClsShared,
                  track::TPCChi2NCl, track::TPCInnerParam, track::TPCSignal, pidtpc::TPCNSigmaEl, pidtpc::TPCNSigmaPi,
                  track::ITSClusterSizes, track::ITSChi2NCl, track::DetectorMap);
DECLARE_SOA_TABLE(PCMFreeLegIdsMix, "AOD", "PCMFREELEGIDMIX", //! twin of PCMFreeLegIds for the mixed instance (never written)
                  o2::soa::Index<>, pcmfreeleg::CollisionId, pcmfreelegid::TrackIdDelta);
} // namespace o2::aod

// one tested pair of stored photons in the ambiguity resolver (also read by the MC diagnostics)
struct SwapRecord {
  CandKey keyI{}, keyJ{};
  float q = -1.f;
  bool complete = false;
  bool swapped = false;
  float gain = -999.f;
  int fail = 0; // bits: 1 stored I refit, 2 stored J refit, 4 crossed I prefilter, 8 crossed J prefilter, 16 crossed I fit, 32 crossed J fit
  int line = 0; // rejection line in fitV0Candidate of the first crossed fit that failed
};

// MC truth of the legs of a genuine photon, and the SVertexer and finder results for this leg pair
struct TrueLegTrace {
  int svLine = -1;     // SVertexer V0 of these legs: -1 none, 0 accepted as candidate, > 0 source line of the rejection in fitV0Candidate
  int finderStep = -9; // finder: -9 pair not tested (tgl window, leg not offered, finder off), -8 both legs with ITS, -7 pair already built, else evalFinderPair
  int finderLine = -1; // finder fit: -1 not fitted, 0 accepted, > 0 source line of the rejection in fitV0Candidate
};

struct DedupDiag {
  std::map<CandKey, std::vector<CandKey>> blockersByKey;            // rejected candidate -> candidates that took one of its legs
  std::set<CandKey> finderKeys;                                     // candidates built by the run-1 finder (the other v0Id = -1 ones come from the swap resolver)
  std::map<std::pair<int64_t, int64_t>, TrueLegTrace> trueLegTrace; // MC: leg pairs of genuine photons, SVertexer fit and finder result
  std::vector<V0CandidateHelper> snapshot;                          // all candidates, before deduplication
  std::vector<CandKey> stored;                                      // surviving candidates
  std::map<CandKey, std::pair<int, float>> lossMargin;
  std::vector<SwapRecord> swaps;
};

template <o2::soa::is_table TV0PhotonKF, o2::soa::is_table TV0legs, o2::soa::is_table TV0legsXYZ, o2::soa::is_table TV0legsDeDxMC, o2::soa::is_table TV0photonsphivpsi, o2::soa::is_table TV0legsCov, o2::soa::is_table TV0photonsKFChi2, bool kMixedInstance = false>
struct PhotonConversionBuilder {

  // using namespace o2;
  // using namespace o2::soa;
  // using namespace o2::framework;
  // using namespace o2::framework::expressions;
  // using namespace o2::constants::physics;
  // using namespace o2::pwgem::photonmeson;

  Produces<TV0PhotonKF> v0photonskf;
  Produces<TV0legs> v0legs;
  Produces<TV0legsXYZ> v0legsXYZ;
  Produces<TV0legsDeDxMC> v0legsDeDxMC;
  Produces<TV0photonsphivpsi> v0photonsphivpsi;
  Produces<TV0legsCov> v0legsCov;
  Produces<TV0photonsKFChi2> v0photonskfchi2; // KF chi2 of the conversion-vertex fit and to the PV
  Produces<std::conditional_t<kMixedInstance, o2::aod::PCMFreeLegsMix, o2::aod::PCMFreeLegs>> pcmfreelegs;
  Produces<std::conditional_t<kMixedInstance, o2::aod::PCMFreeLegIdsMix, o2::aod::PCMFreeLegIds>> pcmfreelegids; // linked derived data: free legs as track indices

  // Produces<aod::V0PhotonsKFCov> v0photonskfcov;
  // Produces<aod::EMEventsNgPCM> events_ngpcm;

  // CCDB options
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  // Operation and minimisation criteria
  Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};
  Configurable<int> useMatCorrType{"useMatCorrType", 0, "0: none, 1: TGeo, 2: LUT"};
  Configurable<int> modeTrackPropagation{"modeTrackPropagation", 0, "0: use real track propagation, including material, 1: use fast approximation using only geometry, 2: Use real track propagation and make comparison to fast propagation (only for debugging and testing)"};
  Configurable<float> propV0LegsRadius{"propV0LegsRadius", 60.f, "Radius to which the V0 legs are propagated to calculate psipair and phiV"};
  o2::framework::Configurable<bool> useSideBasedCorrection{"useSideBasedCorrection", false, "when moving TPC-only tracks, use TPC side flags (with legacy fallback) instead of tgl sign"};
  Configurable<bool> fillLegCov{"fillLegCov", false, "fill V0LegsCov (the table is written in any case, empty if false): the leg track parameters and covariance given to the V0 fit (TPC-only legs moved to the collision), one row per V0Legs row"};
  Configurable<bool> fillFreeLegIds{"fillFreeLegIds", false, "fill PCMFreeLegIds (the table is written in any case, empty if false): the free legs (as for fillFreeLegs) as AO2D (collision, track) indices only, for the mixer on linked derived data (the leg parameters are read from the parent AO2D); only selected collisions"};
  Configurable<bool> fillFreeLegs{"fillFreeLegs", false, "fill PCMFreeLegs (the table is written in any case, empty if false): tracks that pass the leg cuts (as the finder takes them) and are the leg of no stored photon, with the parameters the V0 fit would use (TPC-only legs moved to their own collision); pooled by the mixer with cfgPoolFreeLegs"};

  // deduplication settings
  struct : ConfigurableGroup {
    std::string prefix = "dedup_group";
    Configurable<int> deduplicationMode{"deduplicationMode", 0, "0: Pairwise deduplication, 1: Based on Greedy matching (best score wins), 2: Based on Group matching (crossed pairs are both kept), 3: Keep all V0s, our default in the config is mode 0, however if a wrong configuration is used, the default will be mode 1 (Greedy matching) to avoid crashes"};
    Configurable<float> deduplicationScoreWeight{"deduplicationScoreWeight", 0.5f, "0.:only pca goes into the score, 1: only cosPA goes itno score, any number in between is a mix of pca and cosPA"};
    Configurable<bool> cfgDedupTruthMaps{"cfgDedupTruthMaps", true, "fill the fine (dEta, dPhi, q) truth maps in MC"};
    Configurable<int> dedupMaxGroupSize{"dedupMaxGroupSize", 12, "group matching (mode 2): conflict groups up to this size are solved exactly, larger ones greedily (capped at 16)"};
    Configurable<float> groupMatchingScoreThreshold{"groupMatchingScoreThreshold", -1.f, "dedup mode 2, choice of the best leg-disjoint subset of a conflict group. <= 0 (master): most photons first, then lowest summed score. > 0: largest sum of max(threshold - score, 1e-3) over the photons, so a photon with a score near or above the threshold adds almost nothing and two poor candidates cannot outweigh one good one (e.g. 0.4)"};
    Configurable<float> gapWeight{"gapWeight", 0.f, "dedup score: weight of the helix-circle gap of the legs (0 = off)"};
    Configurable<float> gapScale{"gapScale", 1.f, "dedup score: scale of the circle gap (cm)"};
    Configurable<bool> useCosPA3D{"useCosPA3D", false, "use the 3D cosPA instead of the xy cosPA in V0PhotonCandidate"};
    Configurable<bool> resolveSplitLegs{"resolveSplitLegs", false, "after the dedup, remove photons that share a split track with a better photon - avoids duplicates for photon candidates"};
    Configurable<float> splitLegMaxDTgl{"splitLegMaxDTgl", 0.03f, "split legs: max |tgl1 - tgl2|"};
    Configurable<float> splitLegMaxChi2{"splitLegMaxChi2", 15.f, "split legs: max chi2 of two track pieces"};
    Configurable<bool> dedupHalfSwapRule{"dedupHalfSwapRule", false, "dedup mode 2, after the choice of a conflict group: a non-selected candidate whose two legs went to two different winners replaces both if it scores better than either and the two remaining legs of the winners form no candidate (two fakes on the legs of one photon cannot outweigh it through the score floor); same rule in run 2"};
    Configurable<bool> cfgExtendedDiag{"cfgExtendedDiag", false, "MC: truth of the stored photons, pairs and swaps"};
  } dedup;

  // recovery finder for not thinned data (step 2) and ambiguity resolver settings (step 5)
  struct : ConfigurableGroup {
    std::string prefix = "recovery_group";
    Configurable<bool> runFinder{"runFinder", false, "photon finder for free legs without ITS, candidates enter the dedup"};
    Configurable<float> finderMinScore{"finderMinScore", 0.f, "finder: min. purity score (log-likelihood ratio true / fake)"};
    Configurable<float> finderMaxDeltaTgl{"finderMaxDeltaTgl", 0.15f, "finder: max |tgl+ - tgl-|"};
    Configurable<float> finderMaxCircleGap{"finderMaxCircleGap", 5.f, "finder: max |circle gap| (cm)"};
    Configurable<float> finderMinCosXY{"finderMinCosXY", 0.9f, "finder: min. xy pointing at the touching point"};
    Configurable<float> finderMaxDeltaZ{"finderMaxDeltaZ", 3.f, "finder: max |z+ - z-| at the touching point (cm)"};
    Configurable<float> finderMaxDz{"finderMaxDz", 4.f, "finder: max |z of the photon line at the beam axis - z_PV| (cm)"};
    Configurable<float> finderRMargin{"finderRMargin", 5.f, "finder: margin below the min. conversion radius (cm)"};
    Configurable<float> finderTimeMargin{"finderTimeMargin", 100.f, "finder: margin of the leg time windows (ns)"};
    Configurable<float> finderTimeMarginFwd{"finderTimeMarginFwd", 3000.f, "finder: extra late margin for TPC-only legs (ns)"};
    Configurable<float> finderMaxPCAMixed{"finderMaxPCAMixed", -1.f, "finder: max PCA for an ITS and a non-ITS leg (cm), < 0 builder cuts"};
    Configurable<float> finderFreeCollinearLegs{"finderFreeCollinearLegs", 0.f, "finder: two candidates of one collision with q_inv below this (GeV/c) count as collinear; their legs stay available to the finder (leg pairs already built are not refitted). Recovers true photons next to a close partner photon. 0 = off"};
    Configurable<bool> finderUniformFit{"finderUniformFit", false, "finder candidates (run 1) and finder-origin originals (run 2) are fitted exactly like SVertexer candidates: no side-based drift correction, AP cut applied, no purity gate (symmetric contest between finder-built and SVertexer-built photons)"};
    Configurable<bool> svertexerLegRule{"svertexerLegRule", false, "run 1: an SVertexer V0 candidate with a leg without ITS is fitted like a finder candidate (no AP cut, purity gate; the TPC drift correction of moved legs is useSideBasedCorrection), so that run 1 and run 2 (uniformPairFit) apply one leg-based rule to every pairing"};
    Configurable<bool> finderFreeAllLegs{"finderFreeAllLegs", false, "finder: no leg blocking at all, legs of accepted candidates stay available (leg pairs already built are not refitted), so the candidate set does not depend on the partner photons"};
    Configurable<bool> runSwapResolver{"runSwapResolver", false, "ambiguity resolver: for close stored photon pairs, fit the crossed leg pairing (no AP cut) and keep the pairing with the lower dedup score"};
    Configurable<float> swapMaxQ{"swapMaxQ", 0.1f, "swap resolver: max q_inv of the two photons (GeV/c)"};
    Configurable<float> swapMaxConvDist{"swapMaxConvDist", 20.f, "swap resolver: max distance of the conversion points (cm)"};
  } recovery;

  // builder run 2 of the co-reconstructed mixed event (processMixed)
  struct : ConfigurableGroup {
    std::string prefix = "mixed_group";
    Configurable<bool> prefilterCrossed{"prefilterCrossed", true, "crossed leg pairs (legs of two different run-1 photons): fit only if |dtgl| < recovery finderMaxDeltaTgl and |circle gap| < finderMaxCircleGap (CPU; validation: off)"};
    Configurable<bool> skipSameSourcePairs{"skipSameSourcePairs", true, "mixed event: do not count pairs whose 4 legs all come from the same source event (with pooled conflict groups a group can give two photons of one partner event: a same-event pair inside the mixed event). Data-safe (mixer event counter), pair counting only"};
    Configurable<int> freeLegRuleMask{"freeLegRuleMask", 5, "with freeLegRun1Rules: which run-1 rules a crossed pairing with a pooled free leg gets: bit 0 AP cut, bit 1 purity gate (pairings with a leg without ITS), bit 2 finder fit options (side-based drift correction, finder PCA; pairings with a leg without ITS); 7 = all run-1 rules, 0 = the lenient crossed fit"};
    Configurable<float> freeLegMinPurity{"freeLegMinPurity", -999.f, "with freeLegRun1Rules and bit 1 of freeLegRuleMask: purity threshold for crossed pairings with a pooled free leg; below -998 recovery.finderMinScore (as run 1)"};
    Configurable<bool> mixedDiagnostics{"mixedDiagnostics", false, "MC diagnostics of run 2 (hCrossedFit, four-leg decisions, matched configurations, mixed pair classes, fake classes, hQSym, swap/brems pointing, alternative pairings): CPU on data, no effect on the stored photons"};
    Configurable<bool> diagFreeLegFits{"diagFreeLegFits", false, "diagnostic: every crossed pairing with a pooled free leg is fitted twice as a trial, with the lenient crossed options and with the run-1 rules, and Mixed/Decision/hFreeLegFit records both results (accepted or the rejection line), split by MC truth; the trial candidates are removed again, the real fit follows the configuration (CPU and QA side effect: V0 QA histograms count the trial fits)"};
    Configurable<bool> freeLegRun1Rules{"freeLegRun1Rules", true, "crossed pairings of a pooled free leg (sourceV0Id -2) with a stored leg are fitted as run 1 would have fitted them: with a leg without ITS as a finder candidate (side-based drift correction, AP cut, purity gate), with both legs with ITS as an SVertexer V0 (AP cut, no purity gate); stored x stored pairings keep crossedRelaxAP (overlay mixed event)"};
    Configurable<bool> freeFreeFinderOnly{"freeFreeFinderOnly", false, "with skipFreeFreeCrossed off: a crossed pairing of two pooled free legs is formed only if at least one leg has no ITS (the run-1 finder domain); free x free with both legs with ITS could only have been an SVertexer V0 in run 1 and is skipped (overlay mixed event)"};
    Configurable<bool> skipFreeFreeCrossed{"skipFreeFreeCrossed", true, "crossed pairings of two pooled free legs (sourceV0Id -2 on both) are not formed: run 1 rejected these legs (finder purity gate, AP cut), the run-2 crossed fit would accept them as fakes (overlay mixed event)"};
    Configurable<bool> skipSameSourceCrossed{"skipSameSourceCrossed", true, "crossed pairings of two legs from one source event are not formed (injection mixed event, mixer cfgInject: the current event's pairings were decided by run 1, only the injected legs form new candidates; the same-event artificial event is then exactly the run-1 output). Candidate formation and the hMatched configurations of the mixed event"};
    Configurable<bool> crossedRelaxAP{"crossedRelaxAP", true, "crossed leg pairs are fitted without the AP cut, like the originals without an SVertexer V0 (symmetric fit options, control)"};
    Configurable<bool> refitOriginalsAtSourceVertex{"refitOriginalsAtSourceVertex", true, "the original pairings are refitted against the primary vertex of their own source event (position and covariance carried in PCMMixLegs), so a mixed photon enters the decision and the output exactly as it left run 1; the crossed pairings use the vertex of the artificial event"};
    Configurable<bool> uniformPairFit{"uniformPairFit", false, "decide in run 2 only: every leg pairing of an artificial event (original or crossed, inside or across pooled conflict groups) passes the same prefilter and is fitted with the same options, chosen by the legs as the run-1 finder does (a pair with a leg without ITS gets the finder options: side-based drift correction, no AP cut, purity gate); for run 1 keep-all + mixer cfgPoolGroups, so that the same and the mixed event are decided by identical code on identical kinds of candidates"};
  } mixed;

  // event-level skim for derived data
  struct : ConfigurableGroup {
    std::string prefix = "skim_group";
    Configurable<int> minNv0PerCollision{"minNv0PerCollision", -1, "store V0 photons only for collisions with at least this many surviving candidates (event-level skim for derived data). negative = store all"};
    Configurable<float> minMinvPairPerCollision{"minMinvPairPerCollision", -1.f, "store a collision's photons only if at least one photon pair has qinv (minv) above this (GeV/c^2). negative = disabled"};
    Configurable<float> maxMinvPairPerCollision{"maxMinvPairPerCollision", -1.f, "store a collision's photons only if at least one photon pair has qinv (minv) below this (GeV/c^2). negative = disabled"};
  } skim;

  // single track cuts
  struct : ConfigurableGroup {
    std::string prefix = "trackcut_group";
    Configurable<int> min_ncluster_tpc{"min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<int> mincrossedrows{"mincrossedrows", 40, "min crossed rows"};
    Configurable<bool> moveTPCTracks{"moveTPCTracks", true, "Move TPC-only tracks under the collision assumption"};
    Configurable<bool> disableITSonlyTracks{"disableITSonlyTracks", false, "disable ITSonly tracks in V0 legs"};
    Configurable<bool> disableTPConlyTracks{"disableTPConlyTracks", false, "disable TPConly tracks in V0 legs"};
    Configurable<bool> requireITShit{"requireITShit", false, "require ITS hit to V0 legs"};
    Configurable<float> maxchi2tpc{"maxchi2tpc", 5.0, "max chi2/NclsTPC"}; // default 4.0 + 1.0
    Configurable<float> maxchi2its{"maxchi2its", 6.0, "max chi2/NclsITS"}; // default 5.0 + 1.0
    Configurable<float> maxpt_itsonly{"maxpt_itsonly", 0.15, "max pT for ITSonly tracks at SV"};
    Configurable<float> maxTPCNsigmaEl{"maxTPCNsigmaEl", 4.0, "max. TPC n sigma for electron"};
    Configurable<float> minTPCNsigmaEl{"minTPCNsigmaEl", -4.0, "min. TPC n sigma for electron"};
    Configurable<float> max_frac_shared_clusters_tpc{"max_frac_shared_clusters_tpc", 999.f, "max fraction of shared clusters in TPC"};
    Configurable<float> dcanegtopv{"dcanegtopv", 0.1, "DCA Neg To PV"};
    Configurable<float> dcapostopv{"dcapostopv", 0.1, "DCA Pos To PV"};
    Configurable<float> maxX{"maxX", 83.1, "max X for track IU"};
    Configurable<float> min_pt_trackiu{"min_pt_trackiu", 0.05, "min pT for trackiu"}; // this comes from online processing. pT of track seed is above 50 MeV/c in B = 0.5 T, 20 MeV/c in B = 0.2 T.
  } trackcuts;

  // v0 photon cuts
  struct : ConfigurableGroup {
    std::string prefix = "v0cut_group";
    Configurable<float> min_v0cospa_tpconly{"min_v0cospa_tpconly", 0.99, "min V0 CosPA to V0s with TPConly tracks"}; // double -> N.B. dcos(x)/dx = 0 at x=0)
    Configurable<float> min_v0cospa_its{"min_v0cospa_its", 0.99, "min V0 CosPA to V0s with ITs hits"};
    Configurable<float> max_dcav0dau_tpconly{"max_dcav0dau_tpconly", 3.0, "max distance btween 2 legs to V0s with TPConly tracks"};
    Configurable<float> max_dcav0dau_its{"max_dcav0dau_its", 0.5, "max distance btween 2 legs to V0s with ITS hits"};
    Configurable<float> max_dcav0dau_itsibss{"max_dcav0dau_itsibss", 1.0, "max distance btween 2 legs to V0s with ITS hits on ITSib SS"};
    Configurable<float> max_dcav0dau_tpc_inner_fc{"max_dcav0dau_tpc_inner_fc", 1.5, "max distance btween 2 legs to V0s with ITS hits on TPC inner FC"};
    Configurable<float> min_v0radius{"min_v0radius", 1.0, "min v0 radius"};
    Configurable<float> max_v0radius{"max_v0radius", 90.0, "max v0 radius"};
    Configurable<int> conversionPointMethod{"conversionPointMethod", 0, "0: helix touching point of the legs (Run 2), 1: radius-profiled pair fit through the material, falls back to 0 without a valid minimum"};
    Configurable<float> margin_r_its{"margin_r_its", 3.0, "margin for r cut in cm"};
    Configurable<float> margin_r_tpc{"margin_r_tpc", 7.0, "margin for r cut in cm"};
    Configurable<float> margin_r_itstpc_tpc{"margin_r_itstpc_tpc", 7.0, "margin for r cut in cm"};
    Configurable<float> margin_z{"margin_z", 7.0, "margin for z cut in cm"};
    Configurable<float> max_alpha_ap{"max_alpha_ap", 0.95, "max alpha for AP cut"};
    Configurable<float> max_qt_ap{"max_qt_ap", 0.01, "max qT for AP cut"};
    Configurable<float> min_pt_v0{"min_pt_v0", 0.1, "min pT for v0 photons at PV"};
    Configurable<float> max_pt_v0_itsonly{"max_pt_v0_itsonly", 0.3, "max pT for v0 photons wth 2 ITSonly tracks at PV"};
    Configurable<float> max_eta_v0{"max_eta_v0", 0.9, "max eta for v0 photons at PV"};
    Configurable<float> kfMassConstrain{"kfMassConstrain", -1.f, "mass constrain for the KFParticle mother particle"};
    Configurable<float> max_r_req_its{"max_r_req_its", 16.0, "max Rxy for V0 with ITS hits"};
    Configurable<float> min_r_tpconly{"min_r_tpconly", 32.0, "min Rxy for V0 with TPConly tracks"};
    Configurable<float> max_r_itsmft_ss{"max_r_itsmft_ss", 66.0, "max Rxy for ITS/MFT SS"};
    Configurable<float> max_dcatopv_xy_v0{"max_dcatopv_xy_v0", +1e+10, "max. DCAxy to PV for V0"};
    Configurable<float> max_dcatopv_z_v0{"max_dcatopv_z_v0", +1e+10, "max. DCAz to PV for V0"};
    Configurable<bool> reject_v0_on_itsib{"reject_v0_on_itsib", true, "flag to reject v0s on ITSib"};
  } v0cuts;

  // PCM ML inference
  struct : ConfigurableGroup {
    std::string prefix = "mlcut_group";
    Configurable<bool> applyPCMMl{"applyPCMMl", false, "Flag to apply ML selections"};
    Configurable<bool> use2DBinning{"use2DBinning", false, "Flag to enable/disable 2D binning for ML application"};
    Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};
    Configurable<int> nClassesPCMMl{"nClassesPCMMl", static_cast<int>(o2::analysis::em_cuts_ml::NCutScores), "Number of classes in ML model"};
    Configurable<int> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
    Configurable<int> centTypePCMMl{"centTypePCMMl", 2, "Centrality type for 2D ML application: FT0M:0, FT0A:1, FT0C:2"};
    Configurable<std::vector<int>> cutDirPCMMl{"cutDirPCMMl", std::vector<int>{o2::analysis::em_cuts_ml::vecCutDir}, "Whether to reject score values greater or smaller than the threshold"};
    Configurable<std::vector<std::string>> namesInputFeatures{"namesInputFeatures", std::vector<std::string>{"feature1", "feature2"}, "Names of ML model input features"};
    Configurable<std::vector<std::string>> modelPathsCCDB{"modelPathsCCDB", std::vector<std::string>{"path_ccdb/BDT_PCM/"}, "Paths of models on CCDB"};
    Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"ModelHandler_onnx_PCM.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
    Configurable<std::vector<std::string>> labelsBinsPCMMl{"labelsBinsPCMMl", std::vector<std::string>{"bin 0", "bin 1"}, "Labels for bins"};
    Configurable<std::vector<std::string>> labelsCutScoresPCMMl{"labelsCutScoresPCMMl", std::vector<std::string>{o2::analysis::em_cuts_ml::labelsCutScore}, "Labels for cut scores"};
    Configurable<std::vector<double>> binsPtPCMMl{"binsPtPCMMl", std::vector<double>{0.0, +1e+10}, "pT bin limits for ML application"};
    Configurable<std::vector<double>> binsCentPCMMl{"binsCentPCMMl", std::vector<double>{0.0, 100.0}, "Centrality bin limits for ML application"};
    Configurable<std::vector<double>> cutsPCMMlFlat{"cutsPCMMlFlat", {0.5}, "Flattened ML cuts: [bin0_score0, bin0_score1, ..., binN_scoreM]"};
  } mlcuts;

  o2::analysis::EmMlResponsePCM<float> emMlResponse;
  std::vector<float> outputML;
  V0PhotonCandidate v0photoncandidate;
  o2::ccdb::CcdbApi ccdbApi;

  int mRunNumber{};
  float d_bz{};
  float maxSnp{};  // max sine phi for propagation
  float maxStep{}; // max step size (cm) for propagation
  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb{};
  o2::base::MatLayerCylSet* lut = nullptr;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  o2::aod::common::TPCVDriftManager mVDriftMgr;
  uint64_t mLastCCDBTimestamp{0};
  std::set<CandKey> mRelaxedKeys;
  std::unordered_map<LegDcaKey, LegAtDca, LegDcaKeyHash> mLegDcaCache; // leg DCA to the vertex per DF (fitV0Candidate)
  std::map<CandKey, std::array<int, 4>> mGroupOfKey;
  std::vector<int> mGroupSizeOfCand;
  struct PermuteResult {
    int colID = -1;
    int nMembers = 0;
    bool changed = false;
    bool greedy = false;
  };
  std::vector<PermuteResult> mPermuteResults; // run 2, mode 2: one entry per conflict group with >= 2 candidates
  int mNAnomalyDump = 0;

  o2::framework::HistogramRegistry registry{
    "registry",
    {
      {"hCollisionCounter", "hCollisionCounter", {HistType::kTH1F, {{1, 0.5f, 1.5f}}}},
      {"hNFreeLegs", "free legs filled per DF (fillFreeLegs);legs;DFs", {HistType::kTH1F, {{200, 0.f, 20000.f}}}},
      {"V0/hAP", "Armenteros Podolanski;#alpha;q_{T} (GeV/c)", {HistType::kTH2F, {{200, -1.0f, 1.0f}, {250, 0, 0.25}}}},
      {"V0/hConversionPointXY", "conversion point in XY;X (cm);Y (cm)", {HistType::kTH2F, {{400, -100.0f, 100.0f}, {400, -100.f, 100.f}}}},
      {"V0/hConversionPointRZ", "conversion point in RZ;Z (cm);R_{xy} (cm)", {HistType::kTH2F, {{200, -100.0f, 100.0f}, {200, 0.f, 100.f}}}},
      {"V0/hPt", "pT of V0 at PV;p_{T,#gamma} (GeV/c)", {HistType::kTH1F, {{1000, 0.0f, 10.0f}}}},
      {"V0/hEtaPhi", "#eta vs. #varphi of V0 at PV;#varphi (rad.);#eta", {HistType::kTH2F, {{72, 0.0f, o2::constants::math::TwoPI}, {200, -1, +1}}}},
      {"V0/hCosPA", "cosine of pointing angle;cosine of pointing angle", {HistType::kTH1F, {{100, 0.99f, 1.f}}}},
      {"V0/hCosPA_Rxy", "cosine of pointing angle;r_{xy} (cm);cosine of pointing angle", {HistType::kTH2F, {{200, 0, 100}, {100, 0.99f, 1.f}}}},
      {"V0/hCosPAXY_Rxy", "cosine of pointing angle;r_{xy} (cm);cosine of pointing angle", {HistType::kTH2F, {{200, 0, 100}, {100, 0.99f, 1.f}}}},
      {"V0/hCosPARZ_Rxy", "cosine of pointing angle;r_{xy} (cm);cosine of pointing angle", {HistType::kTH2F, {{200, 0, 100}, {100, 0.99f, 1.f}}}},
      {"V0/hPCA", "distance between 2 legs at SV;PCA (cm)", {HistType::kTH1F, {{500, 0.0f, 5.f}}}},
      {"V0/hPCA_Rxy", "distance between 2 legs at SV;R_{xy} (cm);PCA (cm)", {HistType::kTH2F, {{200, 0, 100}, {500, 0.0f, 5.f}}}},
      {"V0/hPCA_CosPA", "distance between 2 legs at SV vs. cosPA;cosine of pointing angle;PCA (cm)", {HistType::kTH2F, {{100, 0.99, 1}, {500, 0.0f, 5.f}}}},
      {"V0/hDCAxyz", "DCA to PV;DCA_{xy} (cm);DCA_{z} (cm)", {HistType::kTH2F, {{200, -5.f, +5.f}, {200, -5.f, +5.f}}}},
      {"V0/hMeeSV_Rxy", "mee at SV vs. R_{xy};R_{xy} (cm);m_{ee} at SV (GeV/c^{2})", {HistType::kTH2F, {{200, 0.0f, 100.f}, {100, 0, 0.1f}}}},
      {"V0/hRxy_minX_ITSonly_ITSonly", "min trackiu X vs. R_{xy};trackiu X (cm);min trackiu X - R_{xy} (cm)", {HistType::kTH2F, {{100, 0.0f, 100.f}, {100, -50.0, 50.0f}}}},
      {"V0/hRxy_minX_ITSTPC_ITSTPC", "min trackiu X vs. R_{xy};trackiu X (cm);min trackiu X - R_{xy} (cm)", {HistType::kTH2F, {{100, 0.0f, 100.f}, {100, -50.0, 50.0f}}}},
      {"V0/hRxy_minX_ITSTPC_ITSonly", "min trackiu X vs. R_{xy};trackiu X (cm);min trackiu X - R_{xy} (cm)", {HistType::kTH2F, {{100, 0.0f, 100.f}, {100, -50.0, 50.0f}}}},
      {"V0/hRxy_minX_ITSTPC_TPC", "min trackiu X vs. R_{xy};trackiu X (cm);min trackiu X - R_{xy} (cm)", {HistType::kTH2F, {{100, 0.0f, 100.f}, {100, -50.0, 50.0f}}}},
      {"V0/hRxy_minX_TPC_TPC", "min trackiu X vs. R_{xy};trackiu X (cm);min trackiu X - R_{xy} (cm)", {HistType::kTH2F, {{100, 0.0f, 100.f}, {100, -50.0, 50.0f}}}},
      {"V0/hPCA_diffX", "PCA vs. trackiu X - R_{xy};distance btween 2 legs (cm);min trackiu X - R_{xy} (cm)", {HistType::kTH2F, {{500, 0.0f, 5.f}, {100, -50.0, 50.0f}}}},
      {"V0/hPhiVPsiPair", "phiV vs. psi pair;#psi_{pair} (rad.);#phi_{V} (rad.)", {HistType::kTH2F, {{500, -o2::constants::math::PI, o2::constants::math::PI}, {500, 0.0f, o2::constants::math::TwoPI}}}},
      {"V0Leg/hPt", "pT of leg at SV;p_{T,e} (GeV/c)", {HistType::kTH1F, {{1000, 0.0f, 10.0f}}}},
      {"V0Leg/hEtaPhi", "#eta vs. #varphi of leg at SV;#varphi (rad.);#eta", {HistType::kTH2F, {{72, 0.0f, o2::constants::math::TwoPI}, {200, -1, +1}}}},
      {"V0Leg/hRelDeltaPt", "pT resolution;p_{T} (GeV/c);#Deltap_{T}/p_{T}", {HistType::kTH2F, {{1000, 0.f, 10.f}, {100, 0, 1}}}},
      {"V0Leg/hDCAxyz", "DCA xy vs. z to PV;DCA_{xy} (cm);DCA_{z} (cm)", {HistType::kTH2F, {{200, -50.f, 50.f}, {200, -50.f, +50.f}}}},
      {"V0Leg/hdEdx_Pin", "TPC dE/dx vs. p_{in};p_{in} (GeV/c);TPC dE/dx", {HistType::kTH2F, {{1000, 0.f, 10.f}, {200, 0.f, 200.f}}}},
      {"V0Leg/hTPCNsigmaEl", "TPC dE/dx vs. p_{in};p_{in} (GeV/c);n #sigma_{e}^{TPC}", {HistType::kTH2F, {{1000, 0.f, 10.f}, {100, -5.f, +5.f}}}},
      {"V0Leg/hXZ", "track iu x vs. z;z (cm);x (cm)", {HistType::kTH2F, {{200, -100.f, 100.f}, {200, 0.f, 100.f}}}},
    }};

  void init(o2::framework::InitContext&)
  {
    mRunNumber = 0;
    d_bz = 0;
    maxSnp = 0.85f;  // could be changed later
    maxStep = 2.00f; // could be changed later

    static constexpr std::array<const char*, 4> DedupNames = {"pairwise", "greedy matching", "group matching", "keep all"};
    if (dedup.deduplicationMode < 0 || dedup.deduplicationMode >= static_cast<int>(DedupNames.size())) {
      LOG(fatal) << "unknown deduplicationMode " << dedup.deduplicationMode.value;
    }
    LOGF(info, "photon-conversion-builder: deduplicationMode = %d (%s), score weight = %.2f",
         dedup.deduplicationMode.value, DedupNames[dedup.deduplicationMode.value], dedup.deduplicationScoreWeight.value);
    LOGF(info, "photon-conversion-builder: minNv0PerCollision = %d, minMinvPairPerCollision = %f, maxMinvPairPerCollision = %f",
         skim.minNv0PerCollision.value, skim.minMinvPairPerCollision.value, skim.maxMinvPairPerCollision.value);

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    switch (static_cast<MatCorrTypeEM>(useMatCorrType.value)) {
      case MatCorrTypeEM::TGeo:
        LOGF(info, "TGeo correction requested, loading geometry");
        if (!o2::base::GeometryManager::isGeometryLoaded()) {
          ccdb->get<TGeoManager>(geoPath);
        }
        matCorr = o2::base::Propagator::MatCorrType::USEMatCorrTGeo;
        break;
      case MatCorrTypeEM::LUT:
        LOGF(info, "LUT correction requested, loading LUT");
        lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(lutPath));
        matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
        break;
      default:
        LOGF(info, "no correction requested, loading LUT by default!");
        lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(lutPath));
        matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
        break;
    }

    if (mlcuts.applyPCMMl) {
      if (mlcuts.use2DBinning) {
        int binsNPt = static_cast<int>(mlcuts.binsPtPCMMl->size()) - 1;
        int binsNCent = static_cast<int>(mlcuts.binsCentPCMMl->size()) - 1;
        int binsN = binsNPt * binsNCent;
        if (binsN * static_cast<int>(mlcuts.cutDirPCMMl->size()) != static_cast<int>(mlcuts.cutsPCMMlFlat->size())) {
          LOG(fatal) << "Mismatch in number of bins and cuts provided for 2D ML application: binsN * cutDirPCMMl: " << binsN * static_cast<int>(mlcuts.cutDirPCMMl->size()) << " bins vs. cutsPCMMlFlat: " << mlcuts.cutsPCMMlFlat->size() << " cuts";
        }
        if (binsN != static_cast<int>(mlcuts.onnxFileNames->size())) {
          LOG(fatal) << "Mismatch in number of bins and ONNX files provided for 2D ML application: binsN " << binsN << " bins vs. onnxFileNames: " << mlcuts.onnxFileNames->size() << " ONNX files";
        }
        if (binsN != static_cast<int>(mlcuts.labelsBinsPCMMl->size())) {
          LOG(fatal) << "Mismatch in number of bins and labels provided for 2D ML application: binsN:" << binsN << " bins vs. labelsBinsPCMMl: " << mlcuts.labelsBinsPCMMl->size() << " labels";
        }
        if (static_cast<int>(mlcuts.cutDirPCMMl->size()) != mlcuts.nClassesPCMMl) {
          LOG(fatal) << "Mismatch in number of classes and cut directions provided for 2D ML application: nClassesPCMMl: " << mlcuts.nClassesPCMMl << " classes vs. cutDirPCMMl: " << mlcuts.cutDirPCMMl->size() << " cut directions";
        }
        if (static_cast<int>(mlcuts.labelsCutScoresPCMMl->size()) != mlcuts.nClassesPCMMl) {
          LOG(fatal) << "Mismatch in number of labels for cut scores and number of classes provided for 2D ML application: nClassesPCMMl: " << mlcuts.nClassesPCMMl << " classes vs. labelsCutScoresPCMMl: " << mlcuts.labelsCutScoresPCMMl->size() << " labels";
        }
        o2::framework::LabeledArray<double> cutsPCMMl(mlcuts.cutsPCMMlFlat->data(), binsN, mlcuts.nClassesPCMMl, mlcuts.labelsBinsPCMMl, mlcuts.labelsCutScoresPCMMl);
        emMlResponse.configure2D(mlcuts.binsPtPCMMl, mlcuts.binsCentPCMMl, cutsPCMMl, mlcuts.cutDirPCMMl, mlcuts.nClassesPCMMl);
      } else {
        int binsNPt = static_cast<int>(mlcuts.binsPtPCMMl->size()) - 1;
        if (binsNPt * static_cast<int>(mlcuts.cutDirPCMMl->size()) != static_cast<int>(mlcuts.cutsPCMMlFlat->size())) {
          LOG(fatal) << "Mismatch in number of pT bins and cuts provided for ML application: binsNPt * cutDirPCMMl:" << binsNPt * mlcuts.cutDirPCMMl->size() << " bins vs. cutsPCMMlFlat: " << mlcuts.cutsPCMMlFlat->size() << " cuts";
        }
        if (binsNPt != static_cast<int>(mlcuts.onnxFileNames->size())) {
          LOG(fatal) << "Mismatch in number of pT bins and ONNX files provided for ML application: binsNPt " << binsNPt << " bins vs. onnxFileNames: " << mlcuts.onnxFileNames->size() << " ONNX files";
        }
        if (binsNPt != static_cast<int>(mlcuts.labelsBinsPCMMl->size())) {
          LOG(fatal) << "Mismatch in number of pT bins and labels provided for ML application: binsNPt:" << binsNPt << " bins vs. labelsBinsPCMMl: " << mlcuts.labelsBinsPCMMl->size() << " labels";
        }
        if (mlcuts.nClassesPCMMl != static_cast<int>(mlcuts.cutDirPCMMl->size())) {
          LOG(fatal) << "Mismatch in number of classes and cut directions provided for ML application: nClassesPCMMl: " << mlcuts.nClassesPCMMl << " classes vs. cutDirPCMMl: " << mlcuts.cutDirPCMMl->size() << " cut directions";
        }
        if (static_cast<int>(mlcuts.labelsCutScoresPCMMl->size()) != mlcuts.nClassesPCMMl) {
          LOG(fatal) << "Mismatch in number of labels for cut scores and number of classes provided for ML application: nClassesPCMMl:" << mlcuts.nClassesPCMMl << " classes vs. labelsCutScoresPCMMl: " << mlcuts.labelsCutScoresPCMMl->size() << " labels";
        }
        o2::framework::LabeledArray<double> cutsPCMMl(mlcuts.cutsPCMMlFlat->data(), binsNPt, mlcuts.nClassesPCMMl, mlcuts.labelsBinsPCMMl, mlcuts.labelsCutScoresPCMMl);
        emMlResponse.configure(mlcuts.binsPtPCMMl, cutsPCMMl, mlcuts.cutDirPCMMl, mlcuts.nClassesPCMMl);
      }
      if (mlcuts.loadModelsFromCCDB) {
        ccdbApi.init(ccdburl);
        emMlResponse.setModelPathsCCDB(mlcuts.onnxFileNames, ccdbApi, mlcuts.modelPathsCCDB, mlcuts.timestampCCDB);
      } else {
        emMlResponse.setModelPathsLocal(mlcuts.onnxFileNames);
      }
      emMlResponse.cacheInputFeaturesIndices(mlcuts.namesInputFeatures);
      emMlResponse.init();
      if (mlcuts.nClassesPCMMl == o2::analysis::NClassesML::kTwo) {
        registry.add("V0/hBDTBackgroundScoreBeforeCutVsPt", "BDT background score before cut vs pT; pT (GeV/c); BDT background score", {HistType::kTH2F, {{1000, 0.0f, 20.0f}, {1000, 0.0f, 1.0f}}});
        registry.add("V0/hBDTBackgroundScoreAfterCutVsPt", "BDT background score after cut vs pT; pT (GeV/c); BDT background score", {HistType::kTH2F, {{1000, 0.0f, 20.0f}, {1000, 0.0f, 1.0f}}});
        registry.add("V0/hBDTSignalScoreBeforeCutVsPt", "BDT signal score before cut vs pT; pT (GeV/c); BDT signal score", {HistType::kTH2F, {{1000, 0.0f, 20.0f}, {1000, 0.0f, 1.0f}}});
        registry.add("V0/hBDTSignalScoreAfterCutVsPt", "BDT signal score after cut vs pT; pT (GeV/c); BDT signal score", {HistType::kTH2F, {{1000, 0.0f, 20.0f}, {1000, 0.0f, 1.0f}}});
      } else if (mlcuts.nClassesPCMMl == o2::analysis::NClassesML::kThree) {
        registry.add("V0/hBDTBackgroundScoreBeforeCutVsPt", "BDT background score before cut vs pT; pT (GeV/c); BDT background score", {HistType::kTH2F, {{1000, 0.0f, 20.0f}, {1000, 0.0f, 1.0f}}});
        registry.add("V0/hBDTBackgroundScoreAfterCutVsPt", "BDT background score after cut vs pT; pT (GeV/c); BDT background score", {HistType::kTH2F, {{1000, 0.0f, 20.0f}, {1000, 0.0f, 1.0f}}});
        registry.add("V0/hBDTPrimaryPhotonScoreBeforeCutVsPt", "BDT primary photon score before cut vs pT; pT (GeV/c); BDT primary photon score", {HistType::kTH2F, {{1000, 0.0f, 20.0f}, {1000, 0.0f, 1.0f}}});
        registry.add("V0/hBDTPrimaryPhotonScoreAfterCutVsPt", "BDT primary photon score after cut vs pT; pT (GeV/c); BDT primary photon score", {HistType::kTH2F, {{1000, 0.0f, 20.0f}, {1000, 0.0f, 1.0f}}});
        registry.add("V0/hBDTSecondaryPhotonScoreBeforeCutVsPt", "BDT secondary photon score before cut vs pT; pT (GeV/c); BDT secondary photon score", {HistType::kTH2F, {{1000, 0.0f, 20.0f}, {1000, 0.0f, 1.0f}}});
        registry.add("V0/hBDTSecondaryPhotonScoreAfterCutVsPt", "BDT secondary photon score after cut vs pT; pT (GeV/c); BDT secondary photon score", {HistType::kTH2F, {{1000, 0.0f, 20.0f}, {1000, 0.0f, 1.0f}}});
      } else {
        registry.add("V0/hBDTScoreBeforeCutVsPt", "BDT score before cut vs pT; pT (GeV/c); BDT score", {HistType::kTH2F, {{1000, 0.0f, 20.0f}, {1000, 0.0f, 1.0f}}});
        registry.add("V0/hBDTScoreAfterCutVsPt", "BDT score after cut vs pT; pT (GeV/c); BDT score", {HistType::kTH2F, {{1000, 0.0f, 20.0f}, {1000, 0.0f, 1.0f}}});
      }
    }

    // Compare proper propagation and fast geometrical propagation
    if (modeTrackPropagation == TrackPropMode::kBoth) {
      registry.add("V0Leg/hDCAxyPropagationCompare", "Comparison of DCA_{xy} propagation;DCA_{xy} (cm) proper propagation; DCA_{xy} (cm) geom. propagation", {HistType::kTH2F, {{200, -10., 10.}, {200, -10., 10.}}});
      registry.add("V0Leg/hDCAzPropagationCompare", "Comparison of DCA_{z} propagation;DCA_{z} (cm) proper propagation; DCA_{z} (cm) geom. propagation", {HistType::kTH2F, {{200, -10., 10.}, {200, -10., 10.}}});
      registry.add("V0/hPhivPropagationCompare", "Comparison of #phi_{v};#phi_{v} proper propagation; #phi_{v} geom. propagation", {HistType::kTH2F, {{100, 0., 1.6}, {100, 0., 1.6}}});
      registry.add("V0/hPsiPairPropagationCompare", "Comparison of #Psi_{pair};#Psi_{pair} proper propagation; #Psi_{pair} geom. propagation", {HistType::kTH2F, {{100, 0., 1.6}, {100, 0., 1.6}}});
    }

    // Make sure the dedup.deduplicationScoreWeight is between 0 and 1
    if (dedup.deduplicationScoreWeight > 1.f) { // o2-linter: disable=magic-number (score has to be below unity)
      LOG(warning) << "deduplicationScoreWeight is larger than unity which is not allowed";
    }
    if (dedup.deduplicationScoreWeight < 0.f) { // o2-linter: disable=magic-number (score has to be above zero)
      LOG(warning) << "deduplicationScoreWeight is smaller than zero which is not allowed";
    }

    registry.add("V0/hSVFitReject", "SVertexer V0 fits (run 1, dedup pass): source line of the rejecting cut, 0 = accepted (CPU order of the cuts);line;V0s", kTH1F, {{2400, -0.5f, 2399.5f}});
    if (v0cuts.conversionPointMethod == 1) {
      registry.add("V0/hConvPointProfile", "conversionPointMethod 1;0 fallback to the touching point, 1 profile minimum used;candidates", kTH1F, {{2, -0.5f, 1.5f}});
      registry.add("V0/hConvPointProfileShift", "conversionPointMethod 1;R_{profile} - R_{touching point} (cm);candidates", kTH1F, {{200, -10.f, 10.f}});
    }
    if (recovery.runSwapResolver && dedup.deduplicationMode == V0DeduplicationMode::KeepAll) {
      LOG(fatal) << "runSwapResolver needs deduplicationMode 0, 1 or 2";
    }
    if (dedup.resolveSplitLegs) {
      registry.add("Dedup/hSplitLegRemoved", "split legs: photons removed per DF;photons;DFs", kTH1F, {{100, -0.5f, 99.5f}});
    }
    if (dedup.dedupHalfSwapRule) {
      registry.add("Dedup/hHalfSwapRule", "half-swap rule: two winners replaced by the candidate on their legs;replacements per DF;DFs", kTH1F, {{100, -0.5f, 99.5f}});
    }
    if (recovery.runFinder) {
      registry.add("Finder/hSteps", "finder;step passed;pairs", kTH1D, {{FinderStep::kNFindStep, -0.5f, static_cast<float>(FinderStep::kNFindStep) - 0.5f}});
      registry.add("Finder/hScore", "finder;purity score;pairs", kTH1F, {{100, -20.f, 10.f}});
      registry.add("Finder/hNAdded", "finder;photons added per DF;DFs", kTH1F, {{200, 0.f, 2000.f}});
    }
    if (recovery.runSwapResolver) {
      registry.add("Swap/hTested", "swap resolver;0 tested, 1 swapped;pairs", kTH1D, {{2, -0.5f, 1.5f}});
      registry.add("Swap/hQ", "swap resolver;q_{inv} (GeV/c);0 tested, 1 swapped", kTH2F, {{100, 0.f, 0.1f}, {2, -0.5f, 1.5f}});
    }

    if (doprocessMC) {
      const AxisSpec axQ{60, 0.f, 0.3f, "q_{inv}^{true} (GeV/c)"};
      const AxisSpec axDEta{80, -1.6f, 1.6f, "#Delta#eta_{#gamma#gamma}^{true}"};
      const AxisSpec axClass{3, -0.5f, 2.5f, "0 = true, 1 = cross-leg fake, 2 = other fake"};
      // lostByPartnerFake is a subset of lostByCrossFake, and a pair can appear in several of them.
      for (const auto& nm : {"before", "bothSurvive", "bothSurviveSameColl", "surviveAsCrossedPairing", "halfCrossedStored", "oneLost", "bothLost",
                             "lostByCrossFake", "lostByPartnerFake", "lostByOtherFake", "lostByTrue"}) {
        const std::string dir = std::string("MCDedup/Pairs/") + nm + "/";
        registry.add((dir + "hQ").c_str(), "truth photon pairs with >=1 true V0 candidate before deduplication;q_{inv}^{true} (GeV/c);pairs", kTH1F, {axQ}, true);
        registry.add((dir + "hDEta").c_str(), "truth photon pairs with >=1 true V0 candidate before deduplication;#Delta#eta^{true};pairs", kTH1F, {axDEta}, true);
      }
      registry.add("MCDedup/Candidates/hClassStage", "candidate composition;class;0 = before, 1 = stored", kTH2F, {axClass, {2, -0.5f, 1.5f}}, true);

      const AxisSpec axBlocker{4, -0.5f, 3.5f, "0 = true, 1 = cross-leg fake, 2 = other fake, 3 = no blocker (cut, not dedup)"};

      registry.add("MCDedup/Photons/hFate", "fate of the true photons;0 = alive, 1 = lost (blocked), 2 = lost (no blocker);photons", kTH1F, {{3, -0.5f, 2.5f}}, true);
      registry.add("MCDedup/Photons/hFateProps", "true photons with a candidate;0 alive, 1 lost to a cross-leg fake, 2 lost to another candidate, 3 lost without blocker;p_{T} (GeV/c);R_{xy} (cm);legs with ITS;pointing angle (rad);PCA (cm);other candidates sharing a leg", kTHnSparseF, {{4, -0.5f, 3.5f}, {20, 0.f, 2.f}, {18, 0.f, 90.f}, {3, -0.5f, 2.5f}, {16, 0.f, 0.08f}, {15, 0.f, 3.f}, {11, -0.5f, 10.5f}});
      registry.add("MCDedup/FourLeg/hRun1", "run 1: two true photon candidates whose crossed pairings are candidates too (MC);q_{inv} of the true photons (GeV/c);outcome: 0 both true stored, 1 both crossed stored (swap), 2 one true stored, 3 one crossed stored, 4 none, 5 other;isolated 4-candidate decision: 0 true, 1 crossed, 2 other;other candidates using one of the 4 legs;photons related (same mother, or mother e^{#pm});one true candidate is a resolver/finder photon (no SVertexer V0)", kTHnSparseF, {{10, 0.f, 0.05f}, {6, -0.5f, 5.5f}, {3, -0.5f, 2.5f}, {6, -0.5f, 5.5f}, {2, -0.5f, 1.5f}, {2, -0.5f, 1.5f}});
      registry.add("MCDedup/HalfSwap/hSource", "run 1, true photons with a candidate (MC);0 stored, 1 both legs in two other stored candidates (half-swap source), 2 one leg used, 3 no leg used, 4 both legs in one other stored candidate;score S of the best true candidate;objective: c(S_{w1}) + c(S_{w2}) - c(S_{true}), c = max(s_{0} - S, 0.001);(S_{w1} + S_{w2})/2 - S_{true};p_{T} of the true candidate (GeV/c)", kTHnSparseF, {{5, -0.5f, 4.5f}, {40, 0.f, 2.f}, {40, -0.5f, 0.5f}, {40, -1.f, 1.f}, {10, 0.f, 2.f}});
      registry.add("MCDedup/FourLeg/hRun1Origin", "run 1: as hRun1, with the origin of the candidates (MC);q_{inv} of the true photons (GeV/c);outcome: 0 both true stored, 1 swap, 2 one true stored, 3 one crossed stored, 4 none, 5 other;crossed candidates: 0 both SVertexer V0s, 1 one built by the swap resolver, 2 both built by the resolver;true candidates: 0 both SVertexer V0s, 1 one non-SVertexer, 2 both non-SVertexer;other candidates on the 4 legs (0, 1+);photons related", kTHnSparseF, {{10, 0.f, 0.05f}, {6, -0.5f, 5.5f}, {3, -0.5f, 2.5f}, {3, -0.5f, 2.5f}, {2, -0.5f, 1.5f}, {2, -0.5f, 1.5f}});
      registry.add("MCDedup/Inventory/hPhoton", "eligible true photons (MC);#it{p} (GeV/c);R_{conv} (cm);legs with ITS;fate: 0 no true candidate, 1 candidate lost, 2 stored, 3 lost but a leg in a stored crossed candidate;q_{true} to the closest other eligible photon (GeV/c)", kTHnSparseF, {{20, 0.f, 2.f}, {18, 0.f, 90.f}, {3, -0.5f, 2.5f}, {4, -0.5f, 3.5f}, {20, 0.f, 0.1f}});
      registry.add("MCDedup/Inventory/hPair", "eligible pairs of true photons (MC);q_{true} (GeV/c);|E_{1} - E_{2}| / (E_{1} + E_{2});k_{T} (GeV/c);fate: 0 both stored, 1 one stored, 2 none stored and no crossed stored, 3 both crossed stored (swap), 4 one crossed stored, 5 other;true candidates of the pair (0, 1, 2);crossed candidates of the pair (0, 1, 2)", kTHnSparseF, {{20, 0.f, 0.1f}, {10, 0.f, 1.f}, {6, 0.f, 1.5f}, {6, -0.5f, 5.5f}, {3, -0.5f, 2.5f}, {3, -0.5f, 2.5f}});
      registry.add("MCDedup/Inventory/hSwapOrigin", "eligible pairs of true photons with a stored crossed candidate (MC);q_{true} (GeV/c);fate: 0 both crossed stored (swap), 1 one crossed stored;true candidates of the pair (0, 1, 2);stored crossed: 0 SV+SV, 1 SV+finder, 2 finder+finder, 3 SV+resolver, 4 finder+resolver, 5 resolver+resolver;history of photon 1 (the higher code);history of photon 2", kTHnSparseF, {{20, 0.f, 0.1f}, {2, -0.5f, 1.5f}, {3, -0.5f, 2.5f}, {6, -0.5f, 5.5f}, {7, -0.5f, 6.5f}, {7, -0.5f, 6.5f}});
      registry.add("MCDedup/Inventory/hNeverDetail", "true photons without a true candidate in eligible pairs with q_{true} < 20 MeV (MC, entries = pairs);pair fate (as hPair);history (as hSwapOrigin);detail: SVertexer or finder fit rejection line (history 1, 5), finder step (history 2-4: -8 both ITS, -7 built, -9 not tested, else evalFinderPair)", kTH3F, {{6, -0.5f, 5.5f}, {7, -0.5f, 6.5f}, {2010, -10.5f, 1999.5f}});
      registry.add("MCDedup/Inventory/hSwapLoss", "same-event full swaps of two true photons (MC);q_{true} (GeV/c);true candidates of the pair (0, 1, 2);resolver: 0 crossed in different collisions, 1 not tested q_{rec} > swapMaxQ, 2 not tested conv. distance > swapMaxConvDist, 3 true pairing fails the prefilter, 4 true pairing rejected by a fit cut, 5 stored refit failed, 6 complete, crossed score better, 7 complete, gain > 0 but leg taken by another swap, 8 other;gain = S_{crossed} - S_{true} (resolver sign: S_{stored} - S_{alternative})", kTHnSparseF, {{20, 0.f, 0.1f}, {3, -0.5f, 2.5f}, {9, -0.5f, 8.5f}, {40, -2.f, 2.f}});
      registry.add("MCDedup/Inventory/hSwapLossCut", "same-event full swaps whose true pairing the resolver fit rejects (MC, q_{true} < 20 MeV);true candidates of the pair (0, 1, 2);rejection line in fitV0Candidate", kTH2F, {{3, -0.5f, 2.5f}, {2000, -0.5f, 1999.5f}});
      registry.add("MCDedup/Swap/hResolverTruth", "swap resolver, every tested pair of close stored photons (MC);true photons among the 2 stored;true photons among the 2 alternatives;outcome: 0 incomplete (prefilter), 1 incomplete (fit cut), 2 incomplete (stored refit), 3 complete kept, 4 swapped;gain S_{stored} - S_{alternative};q_{inv} of the stored pair (GeV/c)", kTHnSparseF, {{3, -0.5f, 2.5f}, {3, -0.5f, 2.5f}, {5, -0.5f, 4.5f}, {40, -2.f, 2.f}, {20, 0.f, 0.1f}});
      registry.add("MCDedup/Photons/hBlockerClass", "blockers of a killed true photon (ENTRIES = blockers, not photons);class of the blocker;blocker entries", kTH1F, {axBlocker}, true);
      registry.add("MCDedup/Photons/hKillMargin", "mode 0/1 only;S_{victim} - S_{winner};kills", kTH1F, {{200, 0.f, 1.f}}, true);
      registry.add("MCDedup/Photons/hVictimVsWinnerPCA", "PCA of the two;PCA_{victim} (cm);PCA_{winner} (cm)", kTH2F, {{60, 0.f, 3.f}, {60, 0.f, 3.f}}, true);
      // Reason for photon building failure
      // x = 0 -> an equally large solution containing it existed
      // x > 0 -> keeping it would have cost x candidates
      registry.add("MCDedup/Photons/hLostMargin", "why the true photon lost (mode 2);cardinality deficit;#Sigma S(best set with it) - #Sigma S(winner)", kTH2F, {{5, -0.5f, 4.5f}, {100, 0.f, 2.f}}, true);
      if (dedup.cfgDedupTruthMaps) {
        const AxisSpec axDEtaFine{640, -1.6f, 1.6f, "#Delta#eta^{true}"};
        const AxisSpec axDPhi{144, -o2::constants::math::PI, o2::constants::math::PI, "#Delta#varphi^{true} (rad)"};
        for (const auto& nm : {"before", "bothSurvive", "bothLost", "lostByPartnerFake"}) {
          registry.add((std::string("MCDedup/Pairs/") + nm + "/hMap").c_str(), "truth photon pairs", kTHnSparseF, {axDEtaFine, axDPhi, axQ}, true);
        }
      }
      if (dedup.cfgExtendedDiag) {
        const AxisSpec axPt{50, 0.f, 5.f, "p_{T,#gamma} (GeV/c)"};
        const AxisSpec axQrec{50, 0.f, 0.1f, "q_{inv} (GeV/c)"};
        registry.add("MCDedup/Stored/hClassSource", "stored photons;class;0 SVertexer, 1 finder or swap;p_{T,#gamma} (GeV/c)", kTH3F, {axClass, {2, -0.5f, 1.5f}, axPt});
        registry.add("MCDedup/Stored/hPairClass", "stored photon pairs, same collision;0 true-true, 1 crossed doublet, 2 other;q_{inv} (GeV/c)", kTH2F, {{3, -0.5f, 2.5f}, axQrec});
        registry.add("MCDedup/Stored/hSwapTrueCands", "stored crossed doublets (full swaps);true pairings among the dedup candidates (0, 1, 2);q_{inv} (GeV/c)", kTH2F, {{3, -0.5f, 2.5f}, axQrec});
        registry.add("MCDedup/Copies/hChoice", "true photons with a candidate in the dedup sample;collisions with a true candidate;0 = right coll., 1 = wrong (right available), 2 = wrong (no right copy), 3 = lost", kTH2F, {{5, 0.5f, 5.5f}, {4, -0.5f, 3.5f}});
        registry.add("MCDedup/Copies/hRule", "true photons, >= 2 copies, right copy available;0 = current dedup choice, 1 = smallest |#Deltaz| of the photon line, 2 = largest xy cos(PA) (master), 3 = largest 3D cos(PA);0 = wrong, 1 = right", kTH2F, {{4, -0.5f, 3.5f}, {2, -0.5f, 1.5f}});
        registry.add("MCDedup/Copies/hPAxy", "copies of true photons with >= 2 copies;0 = wrong coll., 1 = right coll.;transverse pointing angle (rad)", kTH2F, {{2, -0.5f, 1.5f}, {200, 0.f, 0.2f}});
        registry.add("MCDedup/Copies/hPA3D", "copies of true photons with >= 2 copies;0 = wrong coll., 1 = right coll.;3D pointing angle (rad)", kTH2F, {{2, -0.5f, 1.5f}, {200, 0.f, 0.2f}});
        registry.add("MCDedup/Swap/hCategory", "swap resolver;3 #times true stored + true crossed;0 kept, 1 swapped;q_{inv} (GeV/c)", kTH3F, {{9, -0.5f, 8.5f}, {2, -0.5f, 1.5f}, axQrec});
      }
    }

    if (doprocessMixed && doprocessMixedDerived) {
      LOG(fatal) << "enable processMixed (AO2D input) or processMixedDerived (derived data), not both";
    }
    if (doprocessMixed || doprocessMixedDerived) {
      if (doprocessRec || doprocessMC || doprocessRec_OnlyIfDielectron) {
        LOG(fatal) << "processMixed (builder run 2) has to run in its own instance (photon-conversion-builder-mixed) without processRec / processMC";
      }
      addMixedHistograms();
    }

    mVDriftMgr.setUseSideBasedCorrection(useSideBasedCorrection.value);
  }

  void addMixedHistograms()
  {
    const AxisSpec axType{2, -0.5f, 1.5f, "artificial event (0 same, 1 mixed)"};
    registry.add("Mixed/hEvents", "artificial events", kTH1F, {axType});
    registry.add("Mixed/hNPhotonsIn", "run-1 photons given to run 2;artificial event;photons per event", kTH2F, {axType, {201, -0.5f, 200.5f}});
    registry.add("Mixed/hCandidates", "candidate cut flow;artificial event;0 original pairs, 1 original refit ok, 2 crossed pairs, 3 crossed after pre-cut, 4 crossed fit ok", kTH2F, {axType, {5, -0.5f, 4.5f}});
    registry.add("Mixed/hStored", "stored photons;artificial event;0 original pairing, 1 crossed pairing", kTH2F, {axType, {2, -0.5f, 1.5f}});
    registry.add("Mixed/hNStoredVsIn", "photons per artificial event;artificial event;photons in (run 1);photons out (run 2)", kTH3F, {axType, {101, -0.5f, 100.5f}, {101, -0.5f, 100.5f}});
    registry.add("Mixed/hNOriginalLost", "run-1 photons not stored again;artificial event;lost originals per event", kTH2F, {axType, {21, -0.5f, 20.5f}});
    registry.add("Mixed/hNCrossedStored", "crossed photons stored;artificial event;crossed photons per event", kTH2F, {axType, {21, -0.5f, 20.5f}});
    registry.add("Mixed/hOrigRelDPt", "stored original pairing vs run 1;artificial event;p_{T}^{run 2} / p_{T}^{run 1} - 1", kTH2F, {axType, {400, -0.02f, 0.02f}});
    registry.add("Mixed/hOrigDPhi", "stored original pairing vs run 1;artificial event;#varphi^{run 2} - #varphi^{run 1} (rad)", kTH2F, {axType, {400, -0.02f, 0.02f}});
    registry.add("Mixed/hOrigDEta", "stored original pairing vs run 1;artificial event;#eta^{run 2} - #eta^{run 1}", kTH2F, {axType, {400, -0.02f, 0.02f}});
    registry.add("Mixed/hOrigExact", "stored original pairing vs run 1;artificial event;0 differs, 1 identical p (float)", kTH2F, {axType, {2, -0.5f, 1.5f}});
    registry.add("Mixed/hCrossedSourceQ", "stored crossed photons: the two run-1 photons of its legs;artificial event;q_{inv} (GeV/c)", kTH2F, {axType, {200, 0.f, 0.2f}});
    registry.add("Mixed/hCrossedSourceDist", "stored crossed photons: the two run-1 photons of its legs;artificial event;distance of the conversion points (cm)", kTH2F, {axType, {200, 0.f, 40.f}});
    registry.add("Mixed/hCrossedSourceEvents", "stored crossed photons;artificial event;legs from 0 the same event, 1 two events", kTH2F, {axType, {2, -0.5f, 1.5f}});
    registry.add("Mixed/hStoredTruth", "MC, stored photons;artificial event;0 original true, 1 original fake, 2 crossed true, 3 crossed fake, 4 no label", kTH2F, {axType, {5, -0.5f, 4.5f}});
    registry.add("Mixed/hTimeCandidates", "CPU of the candidate fits;artificial event;time per event (ms)", kTH2F, {axType, {1000, 0.f, 1000.f}});
    registry.add("Mixed/hTimeDF", "CPU of run 2 per DF;time (ms)", kTH1F, {{1000, 0.f, 100000.f}});
    const AxisSpec axPairTruth{5, -0.5f, 4.5f, "pair (0 all, 1 two true photons, 2 full swap of two true photons, 3 other with a fake, 4 no label)"};
    const AxisSpec axQ{500, 0.f, 0.5f, "q_{inv} (GeV/c)"};
    registry.add("Mixed/Pairs/hQBefore", "photon pairs before run 2 (run-1 photons);artificial event", kTH3F, {axType, axPairTruth, axQ});
    registry.add("Mixed/Pairs/hQAfter", "photon pairs after run 2;artificial event", kTH3F, {axType, axPairTruth, axQ});
    const AxisSpec axAlt{9, -0.5f, 8.5f, "3 #times alternatives passing the fit + (0 kept, 1 changed, 2 other)"};
    const AxisSpec axQLow{100, 0.f, 0.1f, "q_{inv} of the input pair (GeV/c)"};
    registry.add("Mixed/Pairs/hAltTrueTrue", "input pairs of two true photons;artificial event", kTH3F, {axType, axAlt, axQLow});
    registry.add("Mixed/Pairs/hFakeClass", "photon pairs by MC class (MC);artificial event;0 before / 1 after run 2;class (see mixFakeClass);q_{inv} (GeV/c);legs with ITS;crossed pairing passed the run-2 fit", kTHnSparseF, {axType, {2, -0.5f, 1.5f}, {18, -0.5f, 17.5f}, {200, 0.f, 0.2f}, {5, -0.5f, 4.5f}, {2, -0.5f, 1.5f}});
    const AxisSpec axDR{std::vector<double>{0., 0.5, 1., 2., 3., 5., 10., 20., 50., 200.}, "|R_{A} - R_{B}| (cm)"};
    registry.add("Mixed/Pairs/hFakeClassOrder", "photon pairs by MC class and conversion order (MC);artificial event;0 before / 1 after run 2;class (see mixFakeClass);q_{inv} (GeV/c);0 softer photon at the smaller radius (kept), 1 softer photon at the larger radius;|R_{A} - R_{B}| (cm)", kTHnSparseF, {axType, {2, -0.5f, 1.5f}, {18, -0.5f, 17.5f}, {40, 0.f, 0.2f}, {2, -0.5f, 1.5f}, axDR});
    registry.add("Mixed/Pairs/hAltSwap", "input pairs = full swap of two true photons;artificial event", kTH3F, {axType, axAlt, axQLow});
    registry.add("Mixed/Pairs/hFakeClassGeo", "photon pairs by MC class and geometry (MC);artificial event;0 before / 1 after run 2;class (see mixFakeClass);q_{inv} (GeV/c);crossed pairing passed the run-2 fit;same physical-primary ancestor;conversion-point distance (cm);log_{10} closest same-sign legs #sqrt{#Delta#varphi^{2}+#Deltatan#lambda^{2}};max TPC shared-cluster fraction", kTHnSparseF, {axType, {2, -0.5f, 1.5f}, {18, -0.5f, 17.5f}, {50, 0.f, 0.1f}, {2, -0.5f, 1.5f}, {2, -0.5f, 1.5f}, {40, 0.f, 40.f}, {30, -5.f, 1.f}, {10, 0.f, 1.f}});
    registry.add("Mixed/Pairs/hGroupCtx", "photon pairs after run 2, conflict group of the run-2 dedup (MC);artificial event;class (see mixFakeClass);q_{inv} (GeV/c);crossed pairing passed the run-2 fit;same conflict group;candidates in the group;legs in the group;photons stored from the group", kTHnSparseF, {axType, {18, -0.5f, 17.5f}, {20, 0.f, 0.1f}, {2, -0.5f, 1.5f}, {2, -0.5f, 1.5f}, {20, 0.5f, 20.5f}, {15, 3.5f, 33.5f}, {10, 0.5f, 10.5f}});
    registry.add("Mixed/Decision/hFourLeg", "pairs of run-1 photons whose two crossed pairings also passed the fit (MC), selected before the decision;artificial event;MC: 0 original assignment = two conversions, 1 crossed assignment = two conversions, 2 neither;isolated 4-candidate decision: 0 original, 1 crossed, 2 other;actual decision: 0 original stored, 1 crossed stored, 2 other;q_{inv} original assignment (GeV/c);q_{inv} crossed assignment (GeV/c);candidates in the conflict group;S_{crossed} - S_{original};an original was fitted without the AP cut", kTHnSparseF, {axType, {3, -0.5f, 2.5f}, {3, -0.5f, 2.5f}, {3, -0.5f, 2.5f}, {50, 0.f, 0.1f}, {50, 0.f, 0.1f}, {12, 3.5f, 15.5f}, {80, -1.f, 1.f}, {2, -0.5f, 1.5f}});
    registry.add("Mixed/Decision/hFourLegKin", "4-leg configurations with min(q_{orig}, q_{crossed}) < 50 MeV (MC);artificial event;MC: 0 original = two conversions, 1 crossed = two conversions, 2 neither;isolated decision: 0 original, 1 crossed, 2 other;#Sigma PCA_{crossed} - #Sigma PCA_{original} (cm);#Sigma #theta_{PA, crossed} - #Sigma #theta_{PA, original} (rad);|#DeltaR_{xy}| of the conversions (cm);|#Deltaz| of the conversions (cm);momentum asymmetry of the more similar same-sign legs;photons related by MC history (brems, siblings)", kTHnSparseF, {axType, {3, -0.5f, 2.5f}, {3, -0.5f, 2.5f}, {40, -2.f, 2.f}, {40, -0.1f, 0.1f}, {20, 0.f, 20.f}, {20, 0.f, 20.f}, {10, 0.f, 1.f}, {2, -0.5f, 1.5f}});
    {
      const AxisSpec axTruth2{2, -0.5f, 1.5f, "MC: 0 original = two conversions, 1 crossed = two conversions"};
      const AxisSpec axChoice{3, -0.5f, 2.5f, "isolated decision: 0 genuine assignment, 1 wrong assignment, 2 other"};
      const AxisSpec axRel2{2, -0.5f, 1.5f, "photons related by MC history"};
      registry.add("Mixed/Decision/hFourLegTruth", "4-leg configurations with a genuine assignment, min(q) < 50 MeV (MC);artificial event", kTHnSparseF, {axType, axTruth2, axChoice, {4, -0.5f, 3.5f, "MC relation: 0 none, 1 radiated by a leg, 2 siblings, 3 leg daughter of a leg"}, {10, 0.f, 0.05f, "q_{inv} of the genuine assignment (GeV/c)"}, {9, -0.5f, 8.5f, "3 #times (genuine photons with S #geq s_{0}) + wrong candidates with S #geq s_{0}"}, {24, -60.f, 60.f, "#Sigma#theta_{PA} wrong - genuine (mrad)"}, {20, -1.f, 1.f, "#Sigma PCA wrong - genuine (cm)"}, {std::vector<double>{0., 0.005, 0.05, 0.1, 0.2, 0.4, 1.}, "max TPC shared-cluster fraction of the 4 legs"}});
      registry.add("Mixed/Decision/hFourLegTrack", "4-leg configurations with a genuine assignment, min(q) < 50 MeV (MC);artificial event", kTHnSparseF, {axType, axTruth2, axChoice, axRel2, {15, -2.f, 1.5f, "log_{10} min distance in the TPC of the legs the wrong assignment pairs (cm)"}, {15, -2.f, 1.5f, "log_{10} min distance in the TPC of the same-sign legs (cm)"}, {5, -0.5f, 4.5f, "legs with ITS"}, {std::vector<double>{-10., 722., 2228., 1.e6}, "occupancy of the artificial event (terciles)"}});
      registry.add("Mixed/Decision/hFourLegResid", "4-leg configurations with a genuine assignment, min(q) < 50 MeV (MC), residuals signed toward the partner;artificial event", kTHnSparseF, {axType, axTruth2, axChoice, axRel2, {20, -0.01f, 0.01f, "#varphi residual photon 1 (rad)"}, {20, -0.01f, 0.01f, "#varphi residual photon 2 (rad)"}, {20, -0.01f, 0.01f, "#eta residual photon 1"}, {20, -0.01f, 0.01f, "#eta residual photon 2"}});
      registry.add("Mixed/Decision/hFourLegAlt", "4-leg configurations with a genuine assignment, min(q) < 50 MeV (MC), isolated decision with other scores;artificial event", kTHnSparseF, {axType, axTruth2, axRel2, {3, -0.5f, 2.5f, "baseline: 0 genuine, 1 wrong, 2 other"}, {3, -0.5f, 2.5f, "pointing only (w = 1)"}, {3, -0.5f, 2.5f, "PCA only (w = 0)"}, {3, -0.5f, 2.5f, "baseline score, raw-sum objective"}, {3, -0.5f, 2.5f, "baseline weight without gap term"}});
      registry.add("Mixed/Decision/hFourLegOrigin", "4-leg configurations with a genuine assignment, min(q) < 50 MeV (MC), common origin of the two genuine photons;artificial event", kTHnSparseF, {axType, axTruth2, axChoice, {4, -0.5f, 3.5f, "MC relation: 0 none, 1 radiated by a leg, 2 siblings, 3 leg daughter of a leg"}, {10, 0.f, 0.05f, "q_{inv} of the genuine assignment (GeV/c)"}, {5, -1.5f, 3.5f, "common origin: -1 n/a, 0 none, 1 same chain root, 2 same great-grandmother, 3 same grandmother"}, {8, -0.5f, 7.5f, "PDG of the chain root: 0 other, 1 pi0, 2 eta, 3 omega, 4 eta', 5 rho/K*, 6 K0, 7 baryon"}});
      registry.add("Mixed/Decision/hFourLegPoint", "4-leg configurations with a genuine assignment, min(q) < 50 MeV (MC), pointing to the PV split into transverse (#varphi) and longitudinal (#lambda) parts;artificial event", kTHnSparseF, {axType, axTruth2, axChoice, axRel2, {5, 0.f, 0.05f, "q_{inv} of the genuine assignment (GeV/c)"}, {24, -30.f, 30.f, "#Sigma|#Delta#varphi| wrong - genuine (mrad)"}, {24, -30.f, 30.f, "#Sigma|#Delta#lambda| wrong - genuine (mrad)"}, {20, -20.f, 20.f, "#Delta#lambda genuine photon 1 (mrad)"}, {20, -20.f, 20.f, "#Delta#lambda genuine photon 2 (mrad)"}, {20, -20.f, 20.f, "#Delta#varphi genuine photon 1 (mrad)"}, {20, -20.f, 20.f, "#Delta#varphi genuine photon 2 (mrad)"}});
    }
    registry.add("Mixed/Decision/hFreeLegFit", "crossed pairings with a pooled free leg (run 2, mixed_group.diagFreeLegFits): trial fits with both option sets (MC);artificial event;truth: 0 legs of two genuine photons with q_{true} < 50 MeV (swap ingredient), 1 other;lenient crossed fit: 0 accepted, else rejection line;run-1 rules: 0 accepted, else rejection line;legs with ITS;free legs in the pairing (1, 2)", kTHnSparseF, {axType, {2, -0.5f, 1.5f}, {2000, -0.5f, 1999.5f}, {2000, -0.5f, 1999.5f}, {3, -0.5f, 2.5f}, {2, 0.5f, 2.5f}});
    registry.add("Mixed/Decision/hCrossedFit", "crossed leg pairings in run 2: result of the fit (MC);artificial event;two genuine photons, q_{true} < 50 MeV;0 accepted, 1 prefilter (only pairs inside the tgl window are visited), else source line of the rejection in fitV0Candidate;legs with ITS;q_{true} of the two photons (GeV/c)", kTHnSparseF, {axType, {2, -0.5f, 1.5f}, {1600, -0.5f, 1599.5f}, {3, -0.5f, 2.5f}, {5, 0.f, 0.05f}});
    registry.add("Mixed/Decision/hMatched", "two genuine MC photons whose 4 legs are all in the artificial event (same definition for SE and ME, q from the MC photons);artificial event;q_{true} (GeV/c);run-1 input pairing: 0 both genuine photons, 1 swapped (both crossed), 2 other;genuine pairings fitted in run 2 (0-2);crossed pairings fitted in run 2 (0-2);run-2 outcome: 0 both genuine stored, 1 both crossed stored, 2 other;photons related (same mother, or one radiated by a leg of the other);PV contributors of the event (charged multiplicity);photons in the artificial event;occupancy of the event;occupancy of the source event - current event (signed, the photon with the larger difference; 0 for SE);TPC-only legs among the 4", kTHnSparseF, {axType, {10, 0.f, 0.05f}, {3, -0.5f, 2.5f}, {3, -0.5f, 2.5f}, {3, -0.5f, 2.5f}, {3, -0.5f, 2.5f}, {2, -0.5f, 1.5f}, {std::vector<double>{0., 50., 100., 200., 400., 800., 1200., 1600., 2000., 2500., 3000., 5000.}, "PV contributors"}, {12, 0.f, 60.f, "photons in the artificial event"}, {std::vector<double>{-10., 722., 2228., 1.e6}, "occupancy"}, {std::vector<double>{-20000., -3000., -2000., -1000., -500., -250., 250., 500., 1000., 2000., 3000., 20000.}, "occupancy partner - current"}, {5, -0.5f, 4.5f, "TPC-only legs"}});
    registry.add("Mixed/Decision/hMatchedCE", "two genuine MC photons whose 4 legs are all in the artificial event, central-electrode variables;artificial event;q_{true} (GeV/c);run-1 input pairing: 0 both genuine photons, 1 swapped, 2 other;run-2 outcome: 0 both genuine stored, 1 both crossed stored, 2 other;photons related;min |#eta| of the two MC photons;min |z| of the two true conversion points (cm);legs crossing the central electrode in the TPC (of 4);TPC-only legs", kTHnSparseF, {axType, {5, 0.f, 0.05f}, {3, -0.5f, 2.5f}, {3, -0.5f, 2.5f}, {2, -0.5f, 1.5f}, {std::vector<double>{0., 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6, 0.9, 2.}, "min |#eta|"}, {std::vector<double>{0., 1., 2., 5., 10., 20., 50., 300.}, "min |z_{conv}| (cm)"}, {5, -0.5f, 4.5f}, {5, -0.5f, 4.5f}});
    registry.add("Mixed/Decision/hMatchedLegs", "two genuine MC photons whose 4 legs are all in the artificial event, partner of each leg in the stored photons;artificial event;q_{true} (GeV/c);run-1 input pairing: 0 both genuine photons, 1 swapped, 2 other;run-2 outcome: 0 both genuine stored, 1 both crossed stored, 2 other;photons related;legs not stored;legs stored with their own genuine partner;legs stored with the crossed partner (other photon of the pair);legs stored with a leg of a third genuine photon (both legs in the event);legs stored with a leg of an incomplete MC photon (arrived inside a stored fake);legs stored with a leg without photon mother", kTHnSparseF, {axType, {5, 0.f, 0.05f}, {3, -0.5f, 2.5f}, {3, -0.5f, 2.5f}, {2, -0.5f, 1.5f}, {5, -0.5f, 4.5f}, {5, -0.5f, 4.5f}, {5, -0.5f, 4.5f}, {5, -0.5f, 4.5f}, {5, -0.5f, 4.5f}, {5, -0.5f, 4.5f}});
    registry.add("Mixed/Decision/hFourLegAll", "pairs of run-1 photons whose two crossed pairings are candidates INCLUDING those appended by the run-2 resolvers (MC);artificial event;MC: 0 original = two conversions, 1 crossed = two conversions, 2 neither;actual decision: 0 original stored, 1 crossed stored, 2 other;q_{inv} of the genuine (truth 0: original, 1: crossed) assignment (GeV/c);crossed candidates: 0 both from the run-2 fit, 1 one appended by a resolver, 2 both appended;original photons: 0 both SVertexer V0s in run 1, 1 one not, 2 both not", kTHnSparseF, {axType, {3, -0.5f, 2.5f}, {3, -0.5f, 2.5f}, {10, 0.f, 0.05f}, {3, -0.5f, 2.5f}, {3, -0.5f, 2.5f}});
    registry.add("Mixed/Decision/hPermute", "conflict groups of the run-2 dedup with >= 2 candidates;artificial event;candidates;reversed order selects other candidates;greedy (no enumeration)", kTHnSparseF, {axType, {30, 0.5f, 30.5f}, {2, -0.5f, 1.5f}, {2, -0.5f, 1.5f}});
    registry.add("Mixed/Pairs/hPairOrigin", "photon pairs after run 2 by MC class and common origin of the two MC photons (MC);artificial event;MC class (see mixFakeClass);q_{inv} (GeV/c);common origin: -1 same photon / no label, 0 none, 1 same chain root, 2 same great-grandmother, 3 same grandmother", kTHnSparseF, {axType, {18, -0.5f, 17.5f}, {50, 0.f, 0.1f}, {5, -1.5f, 3.5f}});
    registry.add("Mixed/Pairs/hSameSource", "mixed event after run 2: pairs whose 4 legs all come from one source event (same-event pairs inside the mixed event, from pooled groups);q_{inv} (GeV/c);MC class (see mixFakeClass)", kTH2F, {{200, 0.f, 0.2f}, {18, -0.5f, 17.5f}});
    registry.add("Mixed/Pairs/hQSym", "photon pairs after run 2, swap-symmetric counting: a stored pair whose two crossed pairings (same 4 legs) also passed the run-2 fit counts 1/2 at its q and 1/2 at q of the crossed pairing;artificial event;MC class of the stored pair (see mixFakeClass);q_{inv} (GeV/c)", kTH3F, {axType, {18, -0.5f, 17.5f}, {500, 0.f, 0.5f}});
    registry.add("Mixed/Pairs/hQClean", "photon pairs, all legs TPC shared fraction < 0.4;artificial event;0 before / 1 after run 2;q_{inv} (GeV/c)", kTH3F, {axType, {2, -0.5f, 1.5f}, {500, 0.f, 0.5f}});
    registry.add("Mixed/Pairs/hBrems", "photon pairs, bremsstrahlung pointing;artificial event;0 before / 1 after run 2;class (see mixFakeClass);MC relation (0 none, 1 radiated by a leg, 2 siblings, 3 leg daughter of a leg);q_{inv} (GeV/c);log_{10} d_{min} (cm);log_{10} angle to the leg (rad);radius of the closest point (cm);B converts further out;p_{T} of B (GeV/c);all legs TPC shared fraction < 0.4", kTHnSparseF, {axType, {2, -0.5f, 1.5f}, {18, -0.5f, 17.5f}, {4, -0.5f, 3.5f}, {50, 0.f, 0.1f}, {40, -3.f, 2.f}, {40, -4.f, 0.5f}, {45, 0.f, 90.f}, {2, -0.5f, 1.5f}, {20, 0.f, 1.f}, {2, -0.5f, 1.5f}});
    registry.add("Mixed/Pairs/hSwapDecision", "full swaps after run 2 (classes 14 and 17);artificial event;0 created in run 2, 1 carried from run 1, 2 one of each;0 true pairings scored better, 1 swap scored better, 2 a true pairing failed the fit;q_{inv} (GeV/c);all legs TPC shared fraction < 0.4", kTHnSparseF, {axType, {3, -0.5f, 2.5f}, {3, -0.5f, 2.5f}, {50, 0.f, 0.1f}, {2, -0.5f, 1.5f}});
  }

  void initCCDB(MyBCs::iterator const& bc)
  {
    initCCDBFor(bc.runNumber(), bc.timestamp(), bc.grpMagField());
  }

  // run number, timestamp and GRPMagField of a BC (AO2D) or of a PMEvents row with EmMagFields (derived data, builder run 2)
  template <typename TGrpMagField>
  void initCCDBFor(const int runNumber, const uint64_t timestamp, TGrpMagField const& grpMagField)
  {
    if (mRunNumber == runNumber) {
      return;
    }
    mRunNumber = runNumber;
    // In case override, don't proceed, please - no CCDB access required
    if (d_bz_input > -990) { // o2-linter: disable=magic-number (override value)
      d_bz = d_bz_input;
      o2::parameters::GRPMagField grpmag;
      if (std::fabs(d_bz) > 1e-5) {                   // o2-linter: disable=magic-number (override value)
        grpmag.setL3Current(30000.f / (d_bz / 5.0f)); // o2-linter: disable=magic-number (override value)
      }
      o2::base::Propagator::initFieldFromGRP(&grpmag);

      return;
    }

    auto run3grp_timestamp = timestamp;
    o2::parameters::GRPMagField grp = grpMagField; // copy: getNominalL3Field is not const
    o2::base::Propagator::initFieldFromGRP(&grp);
    // Fetch magnetic field from ccdb for current collision
    d_bz = grp.getNominalL3Field();
    LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";

    if (useMatCorrType == 2) { // o2-linter: disable=magic-number (material budget correction)
      // setMatLUT only after magfield has been initalized (setMatLUT has implicit and problematic init field call if not)
      o2::base::Propagator::Instance()->setMatLUT(lut);
    }
    /// Set magnetic field for KF vertexing
    const float magneticField = o2::base::Propagator::Instance()->getNominalBz();
    KFParticle::SetField(magneticField);

    mVDriftMgr.init(&ccdb->instance());
  }

  void updateCCDB(MyBCs::iterator const& bc)
  {
    updateCCDBAt(bc.timestamp());
  }

  void updateCCDBAt(const uint64_t timestamp)
  {
    mVDriftMgr.update(timestamp);
    mLastCCDBTimestamp = timestamp;
  }

  std::pair<int8_t, std::set<uint8_t>> its_ib_Requirement = {0, {0, 1, 2}}; // no hit on 3 ITS ib layers.
  template <bool isMC, typename TTrack>
  bool checkV0leg(TTrack const& track)
  {
    using namespace o2::pwgem::photonmeson;
    if constexpr (isMC) {
      if (!track.has_mcParticle()) {
        return false;
      }
    }

    if (trackcuts.disableITSonlyTracks && isITSonlyTrack(track)) {
      return false;
    }

    if (trackcuts.disableTPConlyTracks && isTPConlyTrack(track)) {
      return false;
    }

    if (trackcuts.requireITShit && !track.hasITS()) {
      return false;
    }

    if (track.x() > trackcuts.maxX) {
      return false;
    }

    if (!track.hasITS() && !track.hasTPC()) {
      return false;
    }

    if (track.hasITS() && !track.hasTPC() && (track.hasTRD() || track.hasTOF())) { // remove unrealistic track. this should not happen.
      return false;
    }

    if (track.hasTPC()) {
      if (track.tpcNClsFound() < trackcuts.min_ncluster_tpc) {
        return false;
      }
      if (track.tpcNClsCrossedRows() < trackcuts.mincrossedrows || track.tpcChi2NCl() > trackcuts.maxchi2tpc) {
        return false;
      }
      if (track.tpcFractionSharedCls() > trackcuts.max_frac_shared_clusters_tpc) {
        return false;
      }
      if (track.tpcNSigmaEl() < trackcuts.minTPCNsigmaEl || trackcuts.maxTPCNsigmaEl < track.tpcNSigmaEl()) {
        return false;
      }
    }

    if (track.hasITS()) {
      if (track.itsChi2NCl() > trackcuts.maxchi2its) {
        return false;
      }

      if (v0cuts.reject_v0_on_itsib) {
        auto hits_ib = std::count_if(its_ib_Requirement.second.begin(), its_ib_Requirement.second.end(), [&](auto&& requiredLayer) { return track.itsClusterMap() & (1 << requiredLayer); });
        bool its_ob_only = hits_ib <= its_ib_Requirement.first;
        if (!its_ob_only) {
          return false;
        }
      }
    }

    return true;
  }

  float cospaXY_KF(const KFParticle& kfp, const KFParticle& PV)
  {
    float lx = kfp.GetX() - PV.GetX(); // flight length X
    float ly = kfp.GetY() - PV.GetY(); // flight length Y

    float px = kfp.GetPx();
    float py = kfp.GetPy();
    float cospaXY = RecoDecay::dotProd(std::array{lx, ly}, std::array{px, py}) / (RecoDecay::sqrtSumOfSquares(lx, ly) * RecoDecay::sqrtSumOfSquares(px, py));
    if (cospaXY < -1.f) {
      return -1.f;
    }
    if (cospaXY > 1.f) {
      return 1.f;
    }
    return cospaXY;
  }

  float cospaRZ_KF(const KFParticle& kfp, const KFParticle& PV)
  {
    float lx = kfp.GetX() - PV.GetX();              // flight length X
    float ly = kfp.GetY() - PV.GetY();              // flight length Y
    float lz = kfp.GetZ() - PV.GetZ();              // flight length Z
    float lt = RecoDecay::sqrtSumOfSquares(lx, ly); // flight length R, i.e. transverse plane.

    float pt = RecoDecay::sqrtSumOfSquares(kfp.GetPx(), kfp.GetPy());
    float pz = kfp.GetPz();

    float cospaRZ = RecoDecay::dotProd(std::array{lt, lz}, std::array{pt, pz}) / (RecoDecay::sqrtSumOfSquares(lt, lz) * RecoDecay::sqrtSumOfSquares(pt, pz));
    if (cospaRZ < -1.f) {
      return -1.f;
    }
    if (cospaRZ > 1.f) {
      return 1.f;
    }
    return cospaRZ;
  }

  template <bool isMC, typename TTrack, typename TShiftedTrack>
  void fillTrackTable(TTrack const& track, TShiftedTrack const& shiftedtrack, const float dcaXY, const float dcaZ)
  {
    v0legs(track.collisionId(), track.globalIndex(), track.sign(),
           shiftedtrack.GetPx(), shiftedtrack.GetPy(), shiftedtrack.GetPz(), dcaXY, dcaZ,
           track.tpcNClsFindable(), track.tpcNClsFindableMinusFound(), track.tpcNClsFindableMinusCrossedRows(), track.tpcNClsShared(),
           track.tpcChi2NCl(), track.tpcInnerParam(), track.tpcSignal(),
           track.tpcNSigmaEl(), track.tpcNSigmaPi(),
           track.itsClusterSizes(), track.itsChi2NCl(), track.detectorMap());

    v0legsXYZ(shiftedtrack.GetX(), shiftedtrack.GetY(), shiftedtrack.GetZ());

    if constexpr (isMC) {
      v0legsDeDxMC(track.mcTunedTPCSignal());
    }
  }

  void fillLegCovTable(o2::track::TrackParCov const& trackParCov)
  {
    const auto cov = trackParCov.getCov();
    v0legsCov(trackParCov.getX(), trackParCov.getAlpha(), trackParCov.getY(), trackParCov.getZ(),
              trackParCov.getSnp(), trackParCov.getTgl(), trackParCov.getQ2Pt(), cov.data());
  }

  // free legs (fillFreeLegs): after the tables are written; a track that passes the leg cuts, has a collision and the TPC (as the
  // finder takes them) and is the leg of no stored photon; TPC-only legs moved to their own collision as in the fit
  template <bool isMC, class TBCs, class TCollisions, class TTracks>
  void writeFreeLegs(TCollisions const& collisions, TTracks const& tracks, std::set<int64_t> const& skimRejected)
  {
    std::set<int64_t> usedLegs;
    for (const auto& key : stored_fullv0Ids) {
      usedLegs.insert(key.posId);
      usedLegs.insert(key.eleId);
    }
    int nFree = 0;
    int64_t lastFreeLegTrack = 0; // PCMFreeLegIds: track index stored as the difference to the previous row of this DF
    for (const auto& track : tracks) {
      if (track.collisionId() < 0 || !track.hasTPC() || usedLegs.contains(track.globalIndex()) || skimRejected.contains(track.collisionId()) ||
          track.pt() < trackcuts.min_pt_trackiu || !checkV0leg<isMC>(track)) {
        continue;
      }
      const auto collision = collisions.rawIteratorAt(track.collisionId());
      if (!collision.isSelected()) { // the mixer only takes selected events (PMEvents): free legs of other collisions are never used
        continue;
      }
      auto par = getTrackParCov(track);
      if constexpr (!std::is_void_v<TBCs>) {
        if (trackcuts.moveTPCTracks && o2::pwgem::photonmeson::isTPConlyTrack(track) && !mVDriftMgr.moveTPCTrack<TBCs, TCollisions>(collision, track, par)) {
          continue;
        }
      }
      if (fillFreeLegIds) { // tracks are visited in increasing index order: store the difference to the previous free leg of this DF
        pcmfreelegids(track.collisionId(), static_cast<int>(track.globalIndex() - lastFreeLegTrack));
        lastFreeLegTrack = track.globalIndex();
      }
      if (!fillFreeLegs) {
        ++nFree;
        continue;
      }
      const auto cov = par.getCov();
      pcmfreelegs(track.collisionId(), track.globalIndex(), track.sign(),
                  par.getX(), par.getAlpha(), par.getY(), par.getZ(), par.getSnp(), par.getTgl(), par.getQ2Pt(), cov.data(),
                  track.tpcNClsFindable(), track.tpcNClsFindableMinusFound(), track.tpcNClsFindableMinusCrossedRows(), track.tpcNClsShared(),
                  track.tpcChi2NCl(), track.tpcInnerParam(), track.tpcSignal(), track.tpcNSigmaEl(), track.tpcNSigmaPi(),
                  track.itsClusterSizes(), track.itsChi2NCl(), track.detectorMap());
      ++nFree;
    }
    registry.fill(HIST("hNFreeLegs"), static_cast<float>(nFree));
  }

  template <bool isMC, class TBCs, class TCollisions, class TTracks, typename TV0>
  void fillV0Table(TV0 const& v0, const bool filltable)
  {
    const auto& pos = v0.template posTrack_as<TTracks>();
    const auto& ele = v0.template negTrack_as<TTracks>();
    const auto& collision = v0.template collision_as<TCollisions>();                            // collision where this v0 belongs to.
    const bool finderLike = recovery.svertexerLegRule.value && !(pos.hasITS() && ele.hasITS()); // one leg-based rule for every pairing
    fitV0Candidate<isMC, TBCs, TCollisions, TTracks>(collision, pos, ele, v0.globalIndex(), filltable,
                                                     FitOptions{.nonSVertexer = finderLike, .relaxAP = finderLike, .purityGate = finderLike && recovery.runFinder.value});
    if (!filltable) {
      registry.fill(HIST("V0/hSVFitReject"), static_cast<float>(mLastFitReject));
    }
    if constexpr (isMC) {
      if (auto* t = filltable ? nullptr : trueLegTrace(pos.globalIndex(), ele.globalIndex()); t != nullptr) {
        t->svLine = t->svLine == 0 ? 0 : mLastFitReject; // accepted in any collision wins
      }
    }
  }

  // the purity score (finder LLR tables) is only needed for the finder gate and the split-leg order
  bool usePurity() const { return recovery.runFinder.value || dedup.resolveSplitLegs.value; }

  template <typename TTrack>
  static float purityScore(TTrack const& pos, TTrack const& ele, const float cosPA3D, const float pca, const float mee, const float gap, const float psipair)
  {
    return o2::pwgem::photonmeson::finderPurityScore({cosPA3D, pca, mee, gap, std::fabs(pos.tgl() - ele.tgl()), psipair,
                                                      std::max(std::fabs(pos.tpcNSigmaEl()), std::fabs(ele.tpcNSigmaEl())),
                                                      std::min(std::fabs(pos.tpcNSigmaPi()), std::fabs(ele.tpcNSigmaPi()))});
  }

  // deduplication score, lower = better: master score (cosPA, PCA) + optional circle-gap term
  float dedupScore(const float cosPA, const float pca, const float gap)
  {
    float score = getScoreV0(cosPA, pca, dedup.deduplicationScoreWeight);
    if (dedup.gapWeight > 0.f) {
      score += dedup.gapWeight * std::min(std::fabs(gap) / dedup.gapScale.value, 3.f); // o2-linter: disable=magic-number (cap)
    }
    return score;
  }

  // fit of one leg pair with the V0 cuts
  template <bool isMC, class TBCs, class TCollisions, class TTracks, typename TCollision, typename TTrack>
  void fitV0Candidate(TCollision const& collision, TTrack const& pos, TTrack const& ele, const int64_t v0Id, const bool filltable, FitOptions const& opt = {})
  {
    using namespace o2::pwgem::photonmeson;
    mVDriftMgr.setUseSideBasedCorrection(opt.nonSVertexer || useSideBasedCorrection.value);
    mLastFitReject = 0;
    const float pvX = opt.pv ? opt.pv->posX() : collision.posX(), pvY = opt.pv ? opt.pv->posY() : collision.posY(), pvZ = opt.pv ? opt.pv->posZ() : collision.posZ();

    if (pos.sign() * ele.sign() > 0) { // reject same sign pair
      mLastFitReject = __LINE__;       // rejection point
      return;
    }

    if (pos.pt() < trackcuts.min_pt_trackiu || ele.pt() < trackcuts.min_pt_trackiu) {
      mLastFitReject = __LINE__; // rejection point
      return;
    }

    if (pos.globalIndex() == ele.globalIndex()) {
      mLastFitReject = __LINE__; // rejection point
      return;
    }

    if (isITSonlyTrack(pos) && !ele.hasITS()) {
      mLastFitReject = __LINE__; // rejection point
      return;
    }

    if (isITSonlyTrack(ele) && !pos.hasITS()) {
      mLastFitReject = __LINE__; // rejection point
      return;
    }

    if (!checkV0leg<isMC>(pos) || !checkV0leg<isMC>(ele)) {
      mLastFitReject = __LINE__; // rejection point
      return;
    }

    // LOGF(info, "v0.collisionId() = %d , v0.posTrackId() = %d , v0.negTrackId() = %d", v0.collisionId(), v0.posTrackId(), v0.negTrackId());

    // if(isTPConlyTrack(ele)){
    //   // LOGF(info, "TPConly: ele.globalIndex() = %d, ele.x() = %f, ele.y() = %f, ele.z() = %f, ele.tgl() = %f, ele.alpha() = %f, ele.snp() = %f, ele.signed1Pt() = %f", ele.globalIndex(), ele.x(), ele.y(), ele.z(), ele.tgl(), ele.alpha(), ele.snp(), ele.signed1Pt());
    //   // LOGF(info, "TPConly: ele.globalIndex() = %d, ele.cYY() = %f, ele.cZY() = %f, ele.cZZ() = %f, ele.cSnpY() = %f, ele.cSnpZ() = %f, ele.cSnpSnp() = %f, ele.cTglY() = %f, ele.cTglZ() = %f, ele.cTglSnp() = %f, ele.cTglTgl() = %f, ele.c1PtY() = %f, ele.c1PtZ() = %f, ele.c1PtSnp() = %f, ele.c1PtTgl() = %f, ele.c1Pt21Pt2() = %f", ele.globalIndex(), ele.cYY(), ele.cZY(), ele.cZZ(), ele.cSnpY(), ele.cSnpZ(), ele.cSnpSnp(), ele.cTglY(), ele.cTglZ(), ele.cTglSnp(), ele.cTglTgl(), ele.c1PtY(), ele.c1PtZ(), ele.c1PtSnp(), ele.c1PtTgl(), ele.c1Pt21Pt2());
    // }

    // Calculate DCA with respect to the collision associated to the v0, not individual tracks
    std::array<float, 2> dcaInfo{};

    auto pTrack = getTrackParCov(pos);
    bool failedMovePos = false;
    if constexpr (!std::is_void_v<TBCs>) { // TBCs = void: legs of builder run 2, already moved in run 1
      failedMovePos = trackcuts.moveTPCTracks && isTPConlyTrack(pos) && !mVDriftMgr.moveTPCTrack<TBCs, TCollisions>(collision, pos, pTrack);
    }
    if (failedMovePos) {
      if (opt.probe == nullptr) {
        LOGP(error, "failed correction for positive tpc track");
      }
      mLastFitReject = __LINE__; // rejection point (diagnostic)
      return;
    }
    const bool sideBased = opt.nonSVertexer || useSideBasedCorrection.value;
    const bool cacheDca = modeTrackPropagation == TrackPropMode::kProper;
    std::optional<o2::track::TrackParCov> posAtDca, eleAtDca; // legs at their DCA to the vertex (reused for phi_V / psi_pair)
    const LegDcaKey posKey{pos.globalIndex(), collision.globalIndex(), pvX, pvY, pvZ, sideBased};
    const LegAtDca* posCached = nullptr;
    if (cacheDca) {
      if (const auto it = mLegDcaCache.find(posKey); it != mLegDcaCache.end()) {
        posCached = &it->second;
      }
    }
    auto pTrackC = pTrack;
    o2::math_utils::Point3D<float> vtxPrim{pvX, pvY, pvZ};
    if (modeTrackPropagation != TrackPropMode::kProper) {
      dcaInfo = CalculateDCAFast(pTrackC, vtxPrim, d_bz);
    }
    pTrackC.setPID(o2::track::PID::Electron);

    std::array<float, 2> dcaInfoFast = dcaInfo;
    if (posCached != nullptr) {
      dcaInfo = posCached->dca;
      posAtDca = posCached->trk;
    } else if (modeTrackPropagation != TrackPropMode::kFast) {
      const bool propOk = o2::base::Propagator::Instance()->propagateToDCABxByBz({pvX, pvY, pvZ}, pTrackC, 2.f, matCorr, &dcaInfo);
      if (cacheDca && propOk) { // a failed propagation leaves dcaInfo as it was: not cached
        mLegDcaCache.emplace(posKey, LegAtDca{dcaInfo, pTrackC});
        posAtDca = pTrackC;
      }
      if (filltable && modeTrackPropagation == TrackPropMode::kBoth) {
        registry.fill(HIST("V0Leg/hDCAxyPropagationCompare"), dcaInfo[0], dcaInfoFast[0]);
        registry.fill(HIST("V0Leg/hDCAzPropagationCompare"), dcaInfo[1], dcaInfoFast[1]);
      }
    }
    auto posdcaXY = dcaInfo[0];
    auto posdcaZ = dcaInfo[1];

    auto nTrack = getTrackParCov(ele);
    bool failedMoveEle = false;
    if constexpr (!std::is_void_v<TBCs>) { // TBCs = void: legs of builder run 2, already moved in run 1
      failedMoveEle = trackcuts.moveTPCTracks && isTPConlyTrack(ele) && !mVDriftMgr.moveTPCTrack<TBCs, TCollisions>(collision, ele, nTrack);
    }
    if (failedMoveEle) {
      if (opt.probe == nullptr) {
        LOGP(error, "failed correction for negative tpc track");
      }
      mLastFitReject = __LINE__; // rejection point (diagnostic)
      return;
    }
    const LegDcaKey eleKey{ele.globalIndex(), collision.globalIndex(), pvX, pvY, pvZ, sideBased};
    const LegAtDca* eleCached = nullptr;
    if (cacheDca) {
      if (const auto it = mLegDcaCache.find(eleKey); it != mLegDcaCache.end()) {
        eleCached = &it->second;
      }
    }
    auto nTrackC = nTrack;
    if (modeTrackPropagation != TrackPropMode::kProper) {
      dcaInfo = CalculateDCAFast(nTrackC, vtxPrim, d_bz);
    }
    nTrackC.setPID(o2::track::PID::Electron);

    dcaInfoFast = dcaInfo;
    if (eleCached != nullptr) {
      dcaInfo = eleCached->dca;
      eleAtDca = eleCached->trk;
    } else if (modeTrackPropagation != TrackPropMode::kFast) {
      const bool propOk = o2::base::Propagator::Instance()->propagateToDCABxByBz({pvX, pvY, pvZ}, nTrackC, 2.f, matCorr, &dcaInfo);
      if (cacheDca && propOk) { // a failed propagation leaves dcaInfo as it was: not cached
        mLegDcaCache.emplace(eleKey, LegAtDca{dcaInfo, nTrackC});
        eleAtDca = nTrackC;
      }
      if (filltable && modeTrackPropagation == TrackPropMode::kBoth) {
        registry.fill(HIST("V0Leg/hDCAxyPropagationCompare"), dcaInfo[0], dcaInfoFast[0]);
        registry.fill(HIST("V0Leg/hDCAzPropagationCompare"), dcaInfo[1], dcaInfoFast[1]);
      }
    }
    auto eledcaXY = dcaInfo[0];
    auto eledcaZ = dcaInfo[1];

    if (std::fabs(posdcaXY) < trackcuts.dcapostopv || std::fabs(eledcaXY) < trackcuts.dcanegtopv) {
      mLastFitReject = __LINE__; // rejection point (diagnostic)
      return;
    }

    std::array<float, 3> xyz = {0.f, 0.f, 0.f};
    if (v0cuts.conversionPointMethod != 1) {
      // Calculate the conversion point as the weighted average of the two legs' closest approach points -- avoids expensice KFParticle construction for every candidate
      const float bzNominal = o2::base::Propagator::Instance()->getNominalBz();
      const o2::track::TrackAuxPar cP(pTrack, bzNominal), cN(nTrack, bzNominal);
      const float xC = (cP.xC * cN.rC + cN.xC * cP.rC) / (cP.rC + cN.rC);
      const float yC = (cP.yC * cN.rC + cN.yC * cP.rC) / (cP.rC + cN.rC);
      if (RecoDecay::sqrtSumOfSquares(xC, yC) > trackcuts.maxX + v0cuts.margin_r_tpc) {
        mLastFitReject = __LINE__; // rejection point
        return;
      }
    }
    Vtx_recalculationParCov(o2::base::Propagator::Instance(), pTrack, nTrack, xyz, matCorr);
    if (v0cuts.conversionPointMethod == 1) {
      LegScan scanP, scanN;
      cpScanLeg(pTrack, d_bz, matCorr, scanP);
      cpScanLeg(nTrack, d_bz, matCorr, scanN);
      RadiusProfile prof;
      cpRadiusProfile(scanP, scanN, pos.hasITS() ? 0.f : 0.25f, ele.hasITS() ? 0.f : 0.25f, prof); // o2-linter: disable=magic-number (0.5 cm extra z for legs without ITS)
      const bool useProfile = prof.iMin >= 0 && prof.status != 4;                                  // o2-linter: disable=magic-number (profile class: failed)
      if (useProfile) {
        registry.fill(HIST("V0/hConvPointProfileShift"), prof.rMin - std::hypot(xyz[0], xyz[1]));
        xyz = prof.xyz;
      }
      registry.fill(HIST("V0/hConvPointProfile"), useProfile ? 1.f : 0.f);
    }
    float rxy_tmp = RecoDecay::sqrtSumOfSquares(xyz[0], xyz[1]);
    if (rxy_tmp > trackcuts.maxX + v0cuts.margin_r_tpc) {
      mLastFitReject = __LINE__; // rejection point (diagnostic)
      return;
    }
    if (rxy_tmp < std::fabs(xyz[2]) * std::tan(2 * std::atan(std::exp(-v0cuts.max_eta_v0))) - v0cuts.margin_z) {
      mLastFitReject = __LINE__; // rejection point (diagnostic) // RZ line cut
      return;
    }

    KFPTrack kfp_track_pos = createKFPTrackFromTrackParCov(pTrack, pos.sign(), pos.tpcNClsFound(), pos.tpcChi2NCl());
    KFPTrack kfp_track_ele = createKFPTrackFromTrackParCov(nTrack, ele.sign(), ele.tpcNClsFound(), ele.tpcChi2NCl());
    KFParticle kfp_pos(kfp_track_pos, kPositron);
    KFParticle kfp_ele(kfp_track_ele, kElectron);
    std::array<const KFParticle*, 2> GammaDaughters{&kfp_pos, &kfp_ele};

    KFParticle gammaKF;
    gammaKF.SetConstructMethod(2);
    gammaKF.Construct(GammaDaughters.data(), 2);
    if (v0cuts.kfMassConstrain > -0.1) { // o2-linter: disable=magic-number (nonsensical default value)
      gammaKF.SetNonlinearMassConstraint(v0cuts.kfMassConstrain);
    }
    KFPVertex kfpVertex = opt.pv ? createKFPVertexFromCollision(*opt.pv) : createKFPVertexFromCollision(collision);
    KFParticle KFPV(kfpVertex);
    const float chi2NDFVtx = gammaKF.GetNDF() > 0 ? gammaKF.GetChi2() / gammaKF.GetNDF() : -1.f; // before the PV constraint
    const float chi2PV = gammaKF.GetDeviationFromVertex(KFPV);

    // Transport the gamma to the recalculated decay vertex
    KFParticle gammaKF_DecayVtx = gammaKF; // with respect to (0,0,0)
    gammaKF_DecayVtx.TransportToPoint(xyz.data());

    float cospa_kf = cpaFromKF(gammaKF_DecayVtx, KFPV);
    if (!ele.hasITS() && !pos.hasITS()) {
      if (cospa_kf < v0cuts.min_v0cospa_tpconly) {
        mLastFitReject = __LINE__; // rejection point
        return;
      }
    } else {
      if (cospa_kf < v0cuts.min_v0cospa_its) {
        mLastFitReject = __LINE__; // rejection point
        return;
      }
    }

    float rxy = RecoDecay::sqrtSumOfSquares(gammaKF_DecayVtx.GetX(), gammaKF_DecayVtx.GetY());
    if (rxy > trackcuts.maxX + v0cuts.margin_r_tpc) {
      mLastFitReject = __LINE__; // rejection point
      return;
    }
    if (rxy < std::fabs(gammaKF_DecayVtx.GetZ()) * std::tan(2 * std::atan(std::exp(-v0cuts.max_eta_v0))) - v0cuts.margin_z) {
      mLastFitReject = __LINE__; // rejection point  // RZ line cut
      return;
    }
    if (rxy < v0cuts.min_v0radius || v0cuts.max_v0radius < rxy) {
      mLastFitReject = __LINE__; // rejection point
      return;
    }

    if (!filltable && opt.probe == nullptr) {
      if (isITSTPCTrack(pos) && isITSTPCTrack(ele)) {
        registry.fill(HIST("V0/hRxy_minX_ITSTPC_ITSTPC"), std::min(pTrack.getX(), nTrack.getX()), std::min(pTrack.getX(), nTrack.getX()) - rxy); // trackiu.x() - rxy should be positive
      } else if (isITSonlyTrack(pos) && isITSonlyTrack(ele)) {
        registry.fill(HIST("V0/hRxy_minX_ITSonly_ITSonly"), std::min(pTrack.getX(), nTrack.getX()), std::min(pTrack.getX(), nTrack.getX()) - rxy); // trackiu.x() - rxy should be positive
      } else if ((isITSTPCTrack(pos) && isITSonlyTrack(ele)) || (isITSTPCTrack(ele) && isITSonlyTrack(pos))) {
        registry.fill(HIST("V0/hRxy_minX_ITSTPC_ITSonly"), std::min(pTrack.getX(), nTrack.getX()), std::min(pTrack.getX(), nTrack.getX()) - rxy); // trackiu.x() - rxy should be positive
      } else if (isITSTPCTrack(pos) && !ele.hasITS()) {
        registry.fill(HIST("V0/hRxy_minX_ITSTPC_TPC"), std::min(pTrack.getX(), 83.f), std::min(pTrack.getX(), 83.f) - rxy); // trackiu.x() - rxy should be positive
      } else if (isITSTPCTrack(ele) && !pos.hasITS()) {
        registry.fill(HIST("V0/hRxy_minX_ITSTPC_TPC"), std::min(nTrack.getX(), 83.f), std::min(nTrack.getX(), 83.f) - rxy); // trackiu.x() - rxy should be positive
      } else {
        registry.fill(HIST("V0/hRxy_minX_TPC_TPC"), std::min(83.f, 83.f), std::min(83.f, 83.f) - rxy); // trackiu.x() - rxy should be positive
      }
    }

    if (pos.hasITS() && ele.hasITS()) { // ITSonly-ITSonly, ITSTPC-ITSTPC, ITSTPC-ITSonly
      if (rxy > std::min(pTrack.getX(), nTrack.getX()) + v0cuts.margin_r_its) {
        mLastFitReject = __LINE__; // rejection point
        return;
      }
    } else if (!pos.hasITS() && ele.hasITS()) { // ITSTPC-TPC
      if (rxy > std::min(83.f, nTrack.getX()) + v0cuts.margin_r_itstpc_tpc) {
        mLastFitReject = __LINE__; // rejection point
        return;
      }
    } else if (pos.hasITS() && !ele.hasITS()) { // ITSTPC-TPC
      if (rxy > std::min(pTrack.getX(), 83.f) + v0cuts.margin_r_itstpc_tpc) {
        mLastFitReject = __LINE__; // rejection point
        return;
      }
    } else if (!pos.hasITS() && !ele.hasITS()) { // TPC-TPC
      if (rxy > std::min(83.f, 83.f) + v0cuts.margin_r_tpc) {
        mLastFitReject = __LINE__; // rejection point
        return;
      }
    }

    if ((!pos.hasITS() || !ele.hasITS()) && rxy < v0cuts.max_r_req_its) { // conversion points smaller than v0cuts.max_r_req_its have to be detected with ITS hits.
      mLastFitReject = __LINE__;                                          // rejection point
      return;
    }

    if ((!pos.hasITS() && !ele.hasITS()) && rxy < v0cuts.min_r_tpconly) { // TPConly tracks can detect conversion points larger than v0cuts.min_r_tpconly.
      mLastFitReject = __LINE__;                                          // rejection point
      return;
    }

    // Apply a topological constraint of the gamma to the PV. Parameters will be given at the primary vertex.
    KFParticle gammaKF_PV = gammaKF;
    gammaKF_PV.SetProductionVertex(KFPV);
    float v0pt = RecoDecay::sqrtSumOfSquares(gammaKF_PV.GetPx(), gammaKF_PV.GetPy());
    float v0eta = RecoDecay::eta(std::array{gammaKF_PV.GetPx(), gammaKF_PV.GetPy(), gammaKF_PV.GetPz()});
    float v0phi = RecoDecay::constrainAngle(RecoDecay::phi(gammaKF_PV.GetPx(), gammaKF_PV.GetPy()));

    // KFParticle gammaKF_DecayVtx2 = gammaKF;
    // gammaKF_DecayVtx2.SetProductionVertex(KFPV);
    // gammaKF_DecayVtx2.TransportToPoint(xyz);
    // LOGF(info, "gammaKF_PV.GetPx() = %f, gammaKF_DecayVtx.GetPx() = %f, gammaKF_DecayVtx2.GetPx() = %f", gammaKF_PV.GetPx(), gammaKF_DecayVtx.GetPx(), gammaKF_DecayVtx2.GetPx());
    // LOGF(info, "gammaKF_PV.GetPy() = %f, gammaKF_DecayVtx.GetPy() = %f, gammaKF_DecayVtx2.GetPy() = %f", gammaKF_PV.GetPy(), gammaKF_DecayVtx.GetPy(), gammaKF_DecayVtx2.GetPy());
    // LOGF(info, "gammaKF_PV.GetPz() = %f, gammaKF_DecayVtx.GetPz() = %f, gammaKF_DecayVtx2.GetPz() = %f", gammaKF_PV.GetPz(), gammaKF_DecayVtx.GetPz(), gammaKF_DecayVtx2.GetPz());

    if (std::fabs(v0eta) > v0cuts.max_eta_v0 || v0pt < v0cuts.min_pt_v0) {
      mLastFitReject = __LINE__; // rejection point
      return;
    }

    if (isITSonlyTrack(ele) && isITSonlyTrack(pos) && v0pt > v0cuts.max_pt_v0_itsonly) {
      mLastFitReject = __LINE__; // rejection point
      return;
    }

    // phi_V and psi_pair (4 material propagations) are only stored with the candidate: computed after the cuts (CPU; identical values)
    float phiv = 999.f;
    float psipair = 999.f;
    float phivFast = 999.f;
    float psipairFast = 999.f;
    float baseR = std::hypot(xyz[0], xyz[1]);
    // This method uses the track Helix instead of the full propagation.
    // Hence, it is only an approximation but much faster
    if (modeTrackPropagation != TrackPropMode::kProper) {

      o2::track::TrackAuxPar helixPosEle(nTrack, d_bz);
      o2::track::TrackAuxPar helixPosPos(pTrack, d_bz);

      float diffX = helixPosEle.xC - helixPosPos.xC;
      float diffY = helixPosEle.yC - helixPosPos.yC;
      auto phiHelix = RecoDecay::constrainAngle<float, float>(std::atan2(diffY, diffX) - o2::constants::math::PI / 2.);

      // Electron
      float arcLenghtEle = helixPosEle.rC * 0.9 > propV0LegsRadius ? std::asin(propV0LegsRadius / helixPosEle.rC) * helixPosEle.rC : o2::constants::math::PI / 2.2 * helixPosEle.rC; // This assumes that the photon momentum vector is a tangent of the circle // o2-linter: disable=magic-number (geometrical assumption for propagation)
      auto propTrackEle = getPropMomentumFromTrackHelix(arcLenghtEle, ele, helixPosEle, d_bz / 10., phiHelix - ele.phi());
      // Positron
      float arcLenghtPos = helixPosPos.rC * 0.9 > propV0LegsRadius ? std::asin(propV0LegsRadius / helixPosPos.rC) * helixPosPos.rC : o2::constants::math::PI / 2.2 * helixPosPos.rC; // This assumes that the photon momentum vector is a tangent of the circle // o2-linter: disable=magic-number (geometrical assumption for propagation)
      auto propTrackPos = getPropMomentumFromTrackHelix(arcLenghtPos, pos, helixPosPos, d_bz / 10., phiHelix - pos.phi());

      phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(propTrackPos[0], propTrackPos[1], propTrackPos[2], propTrackEle[0], propTrackEle[1], propTrackEle[2], pos.sign(), ele.sign(), d_bz);
      psipair = o2::aod::pwgem::dilepton::utils::pairutil::getPsiPair(propTrackPos[0], propTrackPos[1], propTrackPos[2], propTrackEle[0], propTrackEle[1], propTrackEle[2]);

      // Store values for later comparison
      phivFast = phiv;
      psipairFast = psipair;
    }
    // This uses the full propagation including material effects
    if (modeTrackPropagation != TrackPropMode::kFast) {
      std::array<float, 3> offsetsR = {propV0LegsRadius, 30.f, 10.f};
      bool pPropagatedSuccess = false;
      bool nPropagatedSuccess = false;
      auto pTrackProp = pTrack;
      auto nTrackProp = nTrack;
      for (const auto& offsetR : offsetsR) {
        // the legs at their DCA to the vertex: the same propagation as for the DCA cut above (reused when it succeeded there)
        if (posAtDca) {
          pTrackProp = *posAtDca;
        } else {
          pTrackProp = pTrack;
          pTrackProp.setPID(o2::track::PID::Electron);
          o2::base::Propagator::Instance()->propagateToDCABxByBz({pvX, pvY, pvZ}, pTrackProp, 2.f, matCorr, &dcaInfo);
        }
        if (eleAtDca) {
          nTrackProp = *eleAtDca;
        } else {
          nTrackProp = nTrack;
          nTrackProp.setPID(o2::track::PID::Electron);
          o2::base::Propagator::Instance()->propagateToDCABxByBz({pvX, pvY, pvZ}, nTrackProp, 2.f, matCorr, &dcaInfo);
        }

        pPropagatedSuccess = o2::base::Propagator::Instance()->propagateToR(pTrackProp, baseR + offsetR);
        nPropagatedSuccess = o2::base::Propagator::Instance()->propagateToR(nTrackProp, baseR + offsetR);

        if (pPropagatedSuccess && nPropagatedSuccess) {
          KFPTrack kfp_track_posProp = createKFPTrackFromTrackParCov(pTrackProp, pos.sign(), pos.tpcNClsFound(), pos.tpcChi2NCl());
          KFPTrack kfp_track_eleProp = createKFPTrackFromTrackParCov(nTrackProp, ele.sign(), ele.tpcNClsFound(), ele.tpcChi2NCl());
          phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(kfp_track_posProp.GetPx(), kfp_track_posProp.GetPy(), kfp_track_posProp.GetPz(), kfp_track_eleProp.GetPx(), kfp_track_eleProp.GetPy(), kfp_track_eleProp.GetPz(), pos.sign(), ele.sign(), d_bz);
          psipair = o2::aod::pwgem::dilepton::utils::pairutil::getPsiPair(kfp_track_posProp.GetPx(), kfp_track_posProp.GetPy(), kfp_track_posProp.GetPz(), kfp_track_eleProp.GetPx(), kfp_track_eleProp.GetPy(), kfp_track_eleProp.GetPz());
          break;
        }
        LOG(debug) << "Propagation to offset" << offsetR << " cm failed for " << (pPropagatedSuccess ? "negative" : "positive") << " track. Trying smaller offset.";
      }
      if (filltable && modeTrackPropagation == TrackPropMode::kBoth) {
        registry.fill(HIST("V0/hPhivPropagationCompare"), phiv, phivFast);
        registry.fill(HIST("V0/hPsiPairPropagationCompare"), psipair, psipairFast);
      }
    }

    if (phiv == 999.f || psipair == 999.f) { // o2-linter: disable=magic-number (nonsensical default value)
      LOG(debug) << "Propagation failed for all radii (" << propV0LegsRadius << ", 30, 10 cm). Using default values for phiv and psipair (999.f).";
    }
    KFParticle kfp_pos_DecayVtx = kfp_pos;         // Don't set Primary Vertex
    KFParticle kfp_ele_DecayVtx = kfp_ele;         // Don't set Primary Vertex
    kfp_pos_DecayVtx.TransportToPoint(xyz.data()); // Don't set Primary Vertex
    kfp_ele_DecayVtx.TransportToPoint(xyz.data()); // Don't set Primary Vertex

    float cospaXYKF = cospaXY_KF(gammaKF_DecayVtx, KFPV);
    float cospaRZKF = cospaRZ_KF(gammaKF_DecayVtx, KFPV);
    auto centType = static_cast<CentType>(mlcuts.centTypePCMMl.value);
    // cosPA slot of the candidate (dedup score, cosPA cut): xy as in master, or 3D with dedup_group.useCosPA3D
    v0photoncandidate.setPhotonCandidate(gammaKF_DecayVtx, gammaKF_PV, pos, kfp_pos_DecayVtx, ele, kfp_ele_DecayVtx, collision, dedup.useCosPA3D ? cospa_kf : cospaXYKF, cospaRZKF, cospaXYKF, psipair, phiv, centType, posdcaXY, eledcaXY, posdcaZ, eledcaZ);

    if (!ele.hasITS() && !pos.hasITS()) { // V0s with TPConly-TPConly
      if (v0cuts.max_r_itsmft_ss < rxy && rxy < trackcuts.maxX + v0cuts.margin_r_tpc) {
        if (v0photoncandidate.getPCA() > v0cuts.max_dcav0dau_tpc_inner_fc) {
          mLastFitReject = __LINE__; // rejection point
          return;
        }
      } else {
        if (v0photoncandidate.getPCA() > v0cuts.max_dcav0dau_tpconly) {
          mLastFitReject = __LINE__; // rejection point
          return;
        }
      }
    } else if (opt.nonSVertexer && recovery.finderMaxPCAMixed > 0.f && pos.hasITS() != ele.hasITS()) {
      if (v0photoncandidate.getPCA() > recovery.finderMaxPCAMixed) {
        mLastFitReject = __LINE__; // rejection point
        return;
      }
    } else { // V0s with ITS hits
      if (rxy < v0cuts.max_r_req_its) {
        if (v0photoncandidate.getPCA() > v0cuts.max_dcav0dau_itsibss) {
          mLastFitReject = __LINE__; // rejection point
          return;
        }
      } else {
        if (v0photoncandidate.getPCA() > v0cuts.max_dcav0dau_its) {
          mLastFitReject = __LINE__; // rejection point
          return;
        }
      }
    }

    if (isITSonlyTrack(pos) && v0photoncandidate.getPosPt() > trackcuts.maxpt_itsonly) {
      mLastFitReject = __LINE__; // rejection point
      return;
    }

    if (isITSonlyTrack(ele) && v0photoncandidate.getElePt() > trackcuts.maxpt_itsonly) {
      mLastFitReject = __LINE__; // rejection point
      return;
    }

    if (v0photoncandidate.getChi2NDF() > 6e+3) { // protection for uint16. // o2-linter: disable=magic-number (protection for uint16)
      mLastFitReject = __LINE__;                 // rejection point
      return;
    }

    if (std::fabs(v0photoncandidate.getDcaXYToPV()) > v0cuts.max_dcatopv_xy_v0 || std::fabs(v0photoncandidate.getDcaZToPV()) > v0cuts.max_dcatopv_z_v0) {
      mLastFitReject = __LINE__; // rejection point
      return;
    }

    if (!opt.relaxAP && !checkAP(v0photoncandidate.getAlpha(), v0photoncandidate.getQt(), v0cuts.max_alpha_ap, v0cuts.max_qt_ap)) { // store only photon conversions
      mLastFitReject = __LINE__;                                                                                                    // rejection point
      return;
    }

    ROOT::Math::PxPyPzMVector vposDedup(kfp_pos_DecayVtx.GetPx(), kfp_pos_DecayVtx.GetPy(), kfp_pos_DecayVtx.GetPz(), o2::constants::physics::MassElectron);
    ROOT::Math::PxPyPzMVector veleDedup(kfp_ele_DecayVtx.GetPx(), kfp_ele_DecayVtx.GetPy(), kfp_ele_DecayVtx.GetPz(), o2::constants::physics::MassElectron);
    const auto meeSV = static_cast<float>((vposDedup + veleDedup).M());

    //  passed all cuts: score and dedup candidate
    const float gap = circleGap(o2::track::TrackAuxPar(pTrack, d_bz), o2::track::TrackAuxPar(nTrack, d_bz));
    const float score = dedupScore(v0photoncandidate.getCosPA(), v0photoncandidate.getPCA(), gap);
    if (opt.probe != nullptr) {
      *opt.probe = {.ok = true, .score = score};
      mLastFitReject = __LINE__; // rejection point (diagnostic)
      return;
    }
    const float purity = usePurity() ? purityScore(pos, ele, cospa_kf, v0photoncandidate.getPCA(), meeSV, gap, v0photoncandidate.getPsiPair()) : -999.f;
    if (opt.purityGate) {
      registry.fill(HIST("Finder/hScore"), purity);
      if (purity <= (opt.minPurity > -998.f ? opt.minPurity : recovery.finderMinScore.value)) { // o2-linter: disable=magic-number (unset marker)
        mLastFitReject = __LINE__;                                                              // rejection point (diagnostic)
        return;
      }
    }

    pca_map[CandKey{v0Id, collision.globalIndex(), pos.globalIndex(), ele.globalIndex()}] = v0photoncandidate.getPCA();
    cospa_map[CandKey{v0Id, collision.globalIndex(), pos.globalIndex(), ele.globalIndex()}] = v0photoncandidate.getCosPA();

    V0CandidateHelper v0Helper(v0Id, collision.globalIndex(), pos.globalIndex(), ele.globalIndex(),
                               v0photoncandidate.getCosPA(), v0photoncandidate.getPCA(), score, meeSV);
    v0Helper.pxAtPV = gammaKF_PV.GetPx();
    v0Helper.pyAtPV = gammaKF_PV.GetPy();
    v0Helper.pzAtPV = gammaKF_PV.GetPz();
    v0Helper.svX = gammaKF_DecayVtx.GetX();
    v0Helper.svY = gammaKF_DecayVtx.GetY();
    v0Helper.svZ = gammaKF_DecayVtx.GetZ();
    v0Helper.pxAtSV = gammaKF_DecayVtx.GetPx();
    v0Helper.pyAtSV = gammaKF_DecayVtx.GetPy();
    v0Helper.pzAtSV = gammaKF_DecayVtx.GetPz();
    v0Helper.purity = purity;
    v0Helper.cosPAxy = cospaXYKF;
    v0Helper.cosPA3D = cospa_kf;
    vecV0Dedup.emplace_back(v0Helper);

    if (mlcuts.applyPCMMl) {
      bool isSelectedML = false;
      std::vector<float> mlInputFeatures = emMlResponse.getInputFeatures(v0photoncandidate, pos, ele);
      if (mlcuts.use2DBinning) {
        isSelectedML = emMlResponse.isSelectedMl(mlInputFeatures, v0photoncandidate.getPt(), v0photoncandidate.getCent(), outputML);
      } else {
        isSelectedML = emMlResponse.isSelectedMl(mlInputFeatures, v0photoncandidate.getPt(), outputML);
      }
      if (filltable) {
        if (mlcuts.nClassesPCMMl == o2::analysis::NClassesML::kTwo) {
          registry.fill(HIST("V0/hBDTBackgroundScoreBeforeCutVsPt"), v0photoncandidate.getPt(), outputML[0]);
          registry.fill(HIST("V0/hBDTSignalScoreBeforeCutVsPt"), v0photoncandidate.getPt(), outputML[1]);
        } else if (mlcuts.nClassesPCMMl == o2::analysis::NClassesML::kThree) {
          registry.fill(HIST("V0/hBDTPrimaryPhotonScoreBeforeCutVsPt"), v0photoncandidate.getPt(), outputML[0]);
          registry.fill(HIST("V0/hBDTSecondaryPhotonScoreBeforeCutVsPt"), v0photoncandidate.getPt(), outputML[1]);
          registry.fill(HIST("V0/hBDTBackgroundScoreBeforeCutVsPt"), v0photoncandidate.getPt(), outputML[2]);
        } else {
          registry.fill(HIST("V0/hBDTScoreBeforeCutVsPt"), v0photoncandidate.getPt(), outputML[0]);
        }
      }
      if (!isSelectedML) {
        mLastFitReject = __LINE__; // rejection point
        return;
      }
      if (filltable) {
        if (mlcuts.nClassesPCMMl == o2::analysis::NClassesML::kTwo) {
          registry.fill(HIST("V0/hBDTBackgroundScoreAfterCutVsPt"), v0photoncandidate.getPt(), outputML[0]);
          registry.fill(HIST("V0/hBDTSignalScoreAfterCutVsPt"), v0photoncandidate.getPt(), outputML[1]);
        } else if (mlcuts.nClassesPCMMl == o2::analysis::NClassesML::kThree) {
          registry.fill(HIST("V0/hBDTPrimaryPhotonScoreAfterCutVsPt"), v0photoncandidate.getPt(), outputML[0]);
          registry.fill(HIST("V0/hBDTSecondaryPhotonScoreAfterCutVsPt"), v0photoncandidate.getPt(), outputML[1]);
          registry.fill(HIST("V0/hBDTBackgroundScoreAfterCutVsPt"), v0photoncandidate.getPt(), outputML[2]);
        } else {
          registry.fill(HIST("V0/hBDTScoreAfterCutVsPt"), v0photoncandidate.getPt(), outputML[0]);
        }
      }
    }

    if (filltable) {
      registry.fill(HIST("V0/hAP"), v0photoncandidate.getAlpha(), v0photoncandidate.getQt());
      registry.fill(HIST("V0/hConversionPointXY"), gammaKF_DecayVtx.GetX(), gammaKF_DecayVtx.GetY());
      registry.fill(HIST("V0/hConversionPointRZ"), gammaKF_DecayVtx.GetZ(), rxy);
      registry.fill(HIST("V0/hPt"), v0photoncandidate.getPt());
      registry.fill(HIST("V0/hEtaPhi"), v0phi, v0eta);
      registry.fill(HIST("V0/hCosPA"), v0photoncandidate.getCosPA());
      registry.fill(HIST("V0/hCosPA_Rxy"), rxy, v0photoncandidate.getCosPA());
      registry.fill(HIST("V0/hPCA"), v0photoncandidate.getPCA());
      registry.fill(HIST("V0/hPCA_CosPA"), v0photoncandidate.getCosPA(), v0photoncandidate.getPCA());
      registry.fill(HIST("V0/hPCA_Rxy"), rxy, v0photoncandidate.getPCA());
      registry.fill(HIST("V0/hDCAxyz"), v0photoncandidate.getDcaXYToPV(), v0photoncandidate.getDcaZToPV());
      registry.fill(HIST("V0/hPCA_diffX"), v0photoncandidate.getPCA(), std::min(pTrack.getX(), nTrack.getX()) - rxy); // trackiu.x() - rxy should be positive
      registry.fill(HIST("V0/hPhiVPsiPair"), v0photoncandidate.getPsiPair(), v0photoncandidate.getPhiV());

      // LOGF(info, "cospa_kf = %f, cospaXY_kf = %f, cospaRZ_kf = %f", cospa_kf, cospaXY_kf, cospaRZ_kf);
      registry.fill(HIST("V0/hCosPAXY_Rxy"), rxy, cospaXYKF);
      registry.fill(HIST("V0/hCosPARZ_Rxy"), rxy, cospaRZKF);

      for (const auto& leg : {kfp_pos_DecayVtx, kfp_ele_DecayVtx}) {
        float legpt = RecoDecay::sqrtSumOfSquares(leg.GetPx(), leg.GetPy());
        float legeta = RecoDecay::eta(std::array{leg.GetPx(), leg.GetPy(), leg.GetPz()});
        float legphi = RecoDecay::constrainAngle(RecoDecay::phi(leg.GetPx(), leg.GetPy()));
        registry.fill(HIST("V0Leg/hPt"), legpt);
        registry.fill(HIST("V0Leg/hEtaPhi"), legphi, legeta);
      } // end of leg loop
      for (const auto& leg : {pos, ele}) {
        registry.fill(HIST("V0Leg/hdEdx_Pin"), leg.tpcInnerParam(), leg.tpcSignal());
        registry.fill(HIST("V0Leg/hTPCNsigmaEl"), leg.tpcInnerParam(), leg.tpcNSigmaEl());
      } // end of leg loop
      for (const auto& leg : {pTrack, nTrack}) {
        registry.fill(HIST("V0Leg/hXZ"), leg.getZ(), leg.getX());
        registry.fill(HIST("V0Leg/hRelDeltaPt"), leg.getPt(), leg.getPt() * std::sqrt(leg.getSigma1Pt2()));
      } // end of leg loop
      registry.fill(HIST("V0Leg/hDCAxyz"), posdcaXY, posdcaZ);
      registry.fill(HIST("V0Leg/hDCAxyz"), eledcaXY, eledcaZ);

      ROOT::Math::PxPyPzMVector vpos_sv(kfp_pos_DecayVtx.GetPx(), kfp_pos_DecayVtx.GetPy(), kfp_pos_DecayVtx.GetPz(), o2::constants::physics::MassElectron);
      ROOT::Math::PxPyPzMVector vele_sv(kfp_ele_DecayVtx.GetPx(), kfp_ele_DecayVtx.GetPy(), kfp_ele_DecayVtx.GetPz(), o2::constants::physics::MassElectron);
      ROOT::Math::PxPyPzMVector v0_sv = vpos_sv + vele_sv;
      registry.fill(HIST("V0/hMeeSV_Rxy"), rxy, v0_sv.M());

      v0photonskf(collision.globalIndex(), v0Id, v0legs.lastIndex() + 1, v0legs.lastIndex() + 2,
                  gammaKF_DecayVtx.GetX(), gammaKF_DecayVtx.GetY(), gammaKF_DecayVtx.GetZ(),
                  gammaKF_PV.GetPx(), gammaKF_PV.GetPy(), gammaKF_PV.GetPz(),
                  v0_sv.M(), v0photoncandidate.getDcaXYToPV(), v0photoncandidate.getDcaZToPV(),
                  cospa_kf, cospaXYKF, cospaRZKF,
                  v0photoncandidate.getPCA(), v0photoncandidate.getAlpha(), v0photoncandidate.getQt(), v0photoncandidate.getChi2NDF());
      v0photonsphivpsi(v0photoncandidate.getPhiV(), v0photoncandidate.getPsiPair());
      v0photonskfchi2(chi2NDFVtx, chi2PV);

      // v0photonskfcov(gammaKF_PV.GetCovariance(9), gammaKF_PV.GetCovariance(14), gammaKF_PV.GetCovariance(20), gammaKF_PV.GetCovariance(13), gammaKF_PV.GetCovariance(19), gammaKF_PV.GetCovariance(18));

      fillTrackTable<isMC>(pos, kfp_pos_DecayVtx, posdcaXY, posdcaZ); // positive leg first
      fillTrackTable<isMC>(ele, kfp_ele_DecayVtx, eledcaXY, eledcaZ); // negative leg second
      if (fillLegCov) {
        fillLegCovTable(pTrack);
        fillLegCovTable(nTrack);
      }
    } // end of fill table
  }

  //  bookkeeping of candidates and stored photons
  static CandKey keyOf(V0CandidateHelper const& c)
  {
    return CandKey{.v0Id = c.v0ID, .collisionId = c.colID, .posId = c.posID, .eleId = c.eleID};
  }

  // MC diagnostics: the trace of a leg pair of one genuine photon, nullptr for any other pair or without diagnostics
  TrueLegTrace* trueLegTrace(const int64_t posId, const int64_t eleId)
  {
    if (mDiag == nullptr) {
      return nullptr;
    }
    const auto ip = mLegMother.find(posId);
    const auto ie = mLegMother.find(eleId);
    if (ip == mLegMother.end() || ie == mLegMother.end() || ip->second != ie->second) {
      return nullptr;
    }
    return &mDiag->trueLegTrace[{posId, eleId}];
  }

  std::map<CandKey, int> candidateIndexByKey() const
  {
    std::map<CandKey, int> index;
    for (int i = 0; i < static_cast<int>(vecV0Dedup.size()); ++i) {
      index[keyOf(vecV0Dedup[i])] = i;
    }
    return index;
  }

  // replaces the stored photons after the dedup and keeps stored_v0Ids and nv0_map consistent
  void setStored(std::vector<CandKey> keys)
  {
    std::stable_sort(keys.begin(), keys.end(), [](const auto& a, const auto& b) { return a.collisionId < b.collisionId; });
    stored_fullv0Ids.swap(keys);
    stored_v0Ids.clear();
    for (auto& [col, n] : nv0_map) { // o2-linter: disable=const-ref-in-for-loop (counters are reset in place)
      n = 0;
    }
    for (const auto& key : stored_fullv0Ids) {
      stored_v0Ids.emplace_back(key.posId, key.eleId);
      nv0_map[key.collisionId]++;
    }
  }

  // finder leg: helix circle and the time window in which the track can belong to a collision (TPC-only legs: from the drift-time asymmetry)
  template <typename TBCs, typename TCollisions, typename TTrack>
  bool makeFinderLeg(TTrack const& track, TCollisions const& collisions, const int64_t bcRef, FinderLeg& leg)
  {
    using namespace o2::pwgem::photonmeson;
    if (track.collisionId() < 0 || !track.hasTPC()) {
      return false;
    }
    leg.id = track.globalIndex();
    leg.tgl = track.tgl();
    leg.noITS = !track.hasITS();
    leg.movable = isTPConlyTrack(track) && (track.flags() & o2::aod::track::TrackTimeAsym) != 0;
    const auto bc = collisions.rawIteratorAt(track.collisionId()).template bc_as<TBCs>();
    leg.t = static_cast<double>(static_cast<int64_t>(bc.globalBC()) - bcRef) * o2::constants::lhc::LHCBunchSpacingNS + track.trackTime();
    if (leg.movable) {
      const int side = tpcDriftSide(track);
      if (side == TPCSideUndetermined) {
        return false;
      }
      leg.side = side;
      leg.tMin = leg.t - track.tpcDeltaTBwd() - recovery.finderTimeMargin;
      leg.tMax = leg.t + track.tpcDeltaTFwd() + recovery.finderTimeMarginFwd;
    } else {
      constexpr float NSigmaTime = 4.f;
      const float sigma = std::max(track.trackTimeRes(), 1.f);
      const float halfWidth = (track.flags() & o2::aod::track::TrackTimeResIsRange) != 0 ? sigma : NSigmaTime * sigma;
      leg.side = 0;
      leg.tMin = leg.t - halfWidth - recovery.finderTimeMargin;
      leg.tMax = leg.t + halfWidth + recovery.finderTimeMargin;
    }
    leg.par = getTrackPar(track);
    leg.circle = o2::track::TrackAuxPar(leg.par, d_bz);
    return true;
  }

  // pairs the free legs of the whole DF (not used by any candidate unless finderFreeCollinearLegs, at least one leg without ITS); a pair that passes the geometry
  // (evalFinderPair) and matches a collision in z is fitted and becomes a dedup candidate with v0Id = -1
  template <bool isMC, class TBCs, class TCollisions, class TTracks>
  void runPhotonFinder(TCollisions const& collisions, TTracks const& tracks)
  {
    std::set<int64_t> usedLegs;
    std::set<std::pair<int64_t, int64_t>> builtPairs;
    // legs of accepted candidates are not free; with finderFreeCollinearLegs, legs of collinear candidate pairs stay free; with
    // finderFreeAllLegs no leg is blocked (only pairs already built are skipped)
    std::vector<bool> collinear(vecV0Dedup.size(), recovery.finderFreeAllLegs.value);
    if (recovery.finderFreeCollinearLegs > 0.f && !recovery.finderFreeAllLegs) {
      std::map<int, std::vector<size_t>> candsOfColl;
      for (size_t i = 0; i < vecV0Dedup.size(); ++i) {
        candsOfColl[vecV0Dedup[i].colID].push_back(i);
      }
      for (const auto& [colId, members] : candsOfColl) {
        for (size_t a = 0; a < members.size(); ++a) {
          for (size_t b = a + 1; b < members.size(); ++b) {
            const auto& ci = vecV0Dedup[members[a]];
            const auto& cj = vecV0Dedup[members[b]];
            if (ci.posID == cj.posID || ci.eleID == cj.eleID) {
              continue;
            }
            if (qInv(ci.pxAtPV, ci.pyAtPV, ci.pzAtPV, cj.pxAtPV, cj.pyAtPV, cj.pzAtPV) < recovery.finderFreeCollinearLegs) {
              collinear[members[a]] = true;
              collinear[members[b]] = true;
            }
          }
        }
      }
    }
    std::vector<char> inCandidate(tracks.size(), 0); // CPU: builtPairs is looked up only when both legs are in a candidate
    for (size_t i = 0; i < vecV0Dedup.size(); ++i) {
      const auto& c = vecV0Dedup[i];
      builtPairs.emplace(c.posID, c.eleID);
      if (c.posID >= 0 && c.posID < static_cast<int64_t>(inCandidate.size()) && c.eleID >= 0 && c.eleID < static_cast<int64_t>(inCandidate.size())) {
        inCandidate[c.posID] = inCandidate[c.eleID] = 1;
      }
      if (!collinear[i]) {
        usedLegs.insert(c.posID);
        usedLegs.insert(c.eleID);
      }
    }
    std::vector<FinderCollision> colls;
    int64_t bcRef = -1;
    for (const auto& col : collisions) {
      if constexpr (isMC) {
        if (!col.has_mcCollision()) {
          continue;
        }
      }
      if (!col.isSelected()) {
        continue;
      }
      const auto bc = col.template bc_as<TBCs>();
      if (bcRef < 0) {
        bcRef = static_cast<int64_t>(bc.globalBC());
      }
      colls.push_back({static_cast<double>(static_cast<int64_t>(bc.globalBC()) - bcRef) * o2::constants::lhc::LHCBunchSpacingNS + col.collisionTime(), col.globalIndex(), col.posX(), col.posY(), col.posZ()});
    }
    if (colls.empty()) {
      return;
    }
    std::sort(colls.begin(), colls.end(), [](const FinderCollision& a, const FinderCollision& b) { return a.t < b.t; });
    const auto* vd = ccdb->getForTimeStamp<o2::tpc::VDriftCorrFact>("TPC/Calib/VDriftTgl", mLastCCDBTimestamp);
    if (vd == nullptr) {
      LOG(warning) << "photon finder: no TPC drift velocity for timestamp " << mLastCCDBTimestamp;
      return;
    }
    const double vDrift = vd->refVDrift * vd->corrFact * 1e-3;
    float xB = 0.f, yB = 0.f;
    for (const auto& c : colls) {
      xB += c.x / static_cast<float>(colls.size());
      yB += c.y / static_cast<float>(colls.size());
    }

    std::vector<FinderLeg> posLegs, eleLegs;
    for (const auto& track : tracks) {
      if (usedLegs.contains(track.globalIndex()) || track.pt() < trackcuts.min_pt_trackiu || !checkV0leg<isMC>(track)) {
        continue;
      }
      FinderLeg leg;
      if (makeFinderLeg<TBCs>(track, collisions, bcRef, leg)) {
        (track.sign() > 0 ? posLegs : eleLegs).push_back(leg);
      }
    }
    std::sort(eleLegs.begin(), eleLegs.end(), [](const FinderLeg& a, const FinderLeg& b) { return a.tgl < b.tgl; });

    const FinderCuts cuts{recovery.finderMaxDeltaTgl, recovery.finderMaxCircleGap, v0cuts.min_r_tpconly, v0cuts.max_r_req_its, recovery.finderRMargin,
                          trackcuts.maxX + v0cuts.margin_r_tpc, recovery.finderMinCosXY, recovery.finderMaxDeltaZ, recovery.finderMaxDz, d_bz};
    const FitOptions finderFit{.nonSVertexer = !recovery.finderUniformFit.value, .purityGate = !recovery.finderUniformFit.value};
    const size_t nBefore = vecV0Dedup.size();
    std::array<double, FinderStep::kNFindStep> stepCount{}; // filled into Finder/hSteps once per DF (CPU: one fill per pair dominated on data)
    auto isBuilt = [&](int64_t a, int64_t b) { return inCandidate[a] && inCandidate[b] && builtPairs.contains({a, b}); };
    // one leg pair, in the order of the original loop (positive legs in track order, negative legs in tgl order)
    auto tryPair = [&](FinderLeg const& p, FinderLeg const& e) {
      TrueLegTrace* trace = nullptr;
      if constexpr (isMC) {
        trace = trueLegTrace(p.id, e.id);
      }
      if ((!p.noITS && !e.noITS) || isBuilt(p.id, e.id)) {
        if (trace != nullptr) {
          trace->finderStep = (!p.noITS && !e.noITS) ? -8 : -7; // o2-linter: disable=magic-number (trace codes)
        }
        return;
      }
      int64_t colId = -1;
      const int step = evalFinderPair(p, e, colls, cuts, xB, yB, vDrift, colId);
      if (step >= 0 && step < FinderStep::kNFindStep) {
        stepCount[step] += 1.;
      }
      if (trace != nullptr) {
        trace->finderStep = step;
      }
      if (step == FinderStep::kFindCollision) {
        fitV0Candidate<isMC, TBCs, TCollisions, TTracks>(collisions.rawIteratorAt(colId), tracks.rawIteratorAt(p.id), tracks.rawIteratorAt(e.id), -1, false, finderFit);
        if (trace != nullptr) {
          trace->finderLine = mLastFitReject;
        }
      }
    };
    const auto tglLess = [](const FinderLeg& c, float v) { return c.tgl < v; };
    // time-binning of the negative legs: the time window of a positive leg is intersected with the time windows of the negative legs in the same bin. Reason for implementation:  faster than checking all negative legs for each positive leg
    const bool timeBinned = (!isMC || mDiag == nullptr) && !eleLegs.empty();
    std::vector<std::vector<uint32_t>> eleInBin;
    double tBin0 = 0.;
    const double binWidth = 5000.; // ns
    if (timeBinned) {
      double tLoAll = eleLegs.front().tMin, tHiAll = eleLegs.front().tMax;
      for (const auto& e : eleLegs) {
        tLoAll = std::min(tLoAll, e.tMin);
        tHiAll = std::max(tHiAll, e.tMax);
      }
      tBin0 = tLoAll;
      eleInBin.resize(static_cast<size_t>((tHiAll - tLoAll) / binWidth) + 1);
      for (uint32_t k = 0; k < eleLegs.size(); ++k) { // k increasing: every bin stays sorted by index = tgl order
        const auto b0 = static_cast<size_t>((eleLegs[k].tMin - tBin0) / binWidth);
        const auto b1 = std::min(static_cast<size_t>((eleLegs[k].tMax - tBin0) / binWidth), eleInBin.size() - 1);
        for (size_t b = b0; b <= b1; ++b) {
          eleInBin[b].push_back(k);
        }
      }
    }
    std::vector<uint32_t> tried;
    for (const auto& p : posLegs) {
      const auto lo = static_cast<uint32_t>(std::lower_bound(eleLegs.begin(), eleLegs.end(), p.tgl - recovery.finderMaxDeltaTgl, tglLess) - eleLegs.begin());
      if (!timeBinned) {
        for (uint32_t k = lo; k < eleLegs.size() && eleLegs[k].tgl <= p.tgl + recovery.finderMaxDeltaTgl; ++k) {
          tryPair(p, eleLegs[k]);
        }
        continue;
      }
      const auto hi = static_cast<uint32_t>(std::upper_bound(eleLegs.begin(), eleLegs.end(), p.tgl + recovery.finderMaxDeltaTgl, [](float v, const FinderLeg& c) { return v < c.tgl; }) - eleLegs.begin());
      if (lo >= hi || p.tMax < tBin0) {
        continue;
      }
      tried.clear();
      const auto b0 = static_cast<size_t>(std::max(0., (p.tMin - tBin0) / binWidth));
      const auto b1 = std::min(static_cast<size_t>((p.tMax - tBin0) / binWidth), eleInBin.size() - 1);
      for (size_t b = b0; b <= b1 && b < eleInBin.size(); ++b) {
        const auto& v = eleInBin[b];
        for (auto it = std::lower_bound(v.begin(), v.end(), lo); it != v.end() && *it < hi; ++it) {
          tried.push_back(*it);
        }
      }
      if (b1 > b0) {
        std::sort(tried.begin(), tried.end());
        tried.erase(std::unique(tried.begin(), tried.end()), tried.end());
      }
      for (const auto k : tried) {
        tryPair(p, eleLegs[k]);
      }
    }
    for (int st = 0; st < FinderStep::kNFindStep; ++st) {
      if (stepCount[st] > 0.) {
        registry.fill(HIST("Finder/hSteps"), static_cast<float>(st), stepCount[st]);
      }
    }
    registry.fill(HIST("Finder/hNAdded"), static_cast<float>(vecV0Dedup.size() - nBefore));
    if (mDiag != nullptr) {
      for (size_t i = nBefore; i < vecV0Dedup.size(); ++i) {
        mDiag->finderKeys.insert(keyOf(vecV0Dedup[i]));
      }
    }
  }

  // two stored photons that use two pieces of one track (same tgl, compatible parameters): keep the better one
  template <class TTracks>
  void resolveSplitLegConflicts(TTracks const& tracks, DedupDiag* diag)
  {
    const auto indexOfKey = candidateIndexByKey();
    std::vector<std::pair<float, size_t>> order;
    order.reserve(stored_fullv0Ids.size());
    for (size_t k = 0; k < stored_fullv0Ids.size(); ++k) {
      const auto it = indexOfKey.find(stored_fullv0Ids[k]);
      order.emplace_back(it == indexOfKey.end() ? std::numeric_limits<float>::max() : (usePurity() ? -vecV0Dedup[it->second].purity : vecV0Dedup[it->second].score), k);
    }
    std::stable_sort(order.begin(), order.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
    auto samePiece = [&](int64_t idA, int64_t idB) {
      if (idA == idB) {
        return false;
      }
      const auto trA = tracks.rawIteratorAt(idA);
      const auto trB = tracks.rawIteratorAt(idB);
      if (std::fabs(trA.tgl() - trB.tgl()) > dedup.splitLegMaxDTgl) {
        return false;
      }
      const float chi2 = splitPieceChi2(trA, trB, d_bz);
      return chi2 >= 0.f && chi2 < dedup.splitLegMaxChi2;
    };
    // the keptPosOf and keptEleOf maps are keyed by collisionId (or 0 for mixed instances) and contain a map of tgl to the kept leg information
    struct KeptLeg {
      int64_t id;
      size_t owner;
      size_t order;
    };
    std::unordered_map<int64_t, std::multimap<float, KeptLeg>> keptPosOf, keptEleOf;
    std::unordered_map<int64_t, float> tglOf;
    auto tgl = [&](int64_t id) {
      const auto it = tglOf.find(id);
      return it != tglOf.end() ? it->second : (tglOf[id] = tracks.rawIteratorAt(id).tgl());
    };
    auto firstOwner = [&](std::multimap<float, KeptLeg> const& kept, int64_t id) {
      const float t = tgl(id);
      size_t bestOrder = std::numeric_limits<size_t>::max();
      int64_t owner = -1;
      for (auto it = kept.lower_bound(t - dedup.splitLegMaxDTgl); it != kept.end() && it->first <= t + dedup.splitLegMaxDTgl; ++it) {
        if (it->second.order < bestOrder && samePiece(id, it->second.id)) {
          bestOrder = it->second.order;
          owner = static_cast<int64_t>(it->second.owner);
        }
      }
      return owner;
    };
    std::vector<char> removed(stored_fullv0Ids.size(), 0);
    int nRemoved = 0;
    for (const auto& [rank, k] : order) {
      const auto& key = stored_fullv0Ids[k];
      auto& keptPos = keptPosOf[kMixedInstance ? static_cast<int64_t>(key.collisionId) : 0];
      auto& keptEle = keptEleOf[kMixedInstance ? static_cast<int64_t>(key.collisionId) : 0];
      int64_t owner = firstOwner(keptPos, key.posId);
      if (owner < 0) {
        owner = firstOwner(keptEle, key.eleId);
      }
      if (owner >= 0) {
        removed[k] = 1;
        ++nRemoved;
        if (diag != nullptr) {
          diag->blockersByKey[key].push_back(stored_fullv0Ids[owner]);
        }
        continue;
      }
      keptPos.emplace(tgl(key.posId), KeptLeg{key.posId, k, keptPos.size()});
      keptEle.emplace(tgl(key.eleId), KeptLeg{key.eleId, k, keptEle.size()});
    }
    if (nRemoved > 0) {
      std::vector<CandKey> kept;
      kept.reserve(stored_fullv0Ids.size() - nRemoved);
      for (size_t k = 0; k < stored_fullv0Ids.size(); ++k) {
        if (removed[k] == 0) {
          kept.push_back(stored_fullv0Ids[k]);
        }
      }
      setStored(std::move(kept));
    }
    registry.fill(HIST("Dedup/hSplitLegRemoved"), static_cast<float>(nRemoved));
  }

  // ambiguity resolver
  // two stored photons of one collision, close in q and in space: fit the crossed leg pairing; if both crossed photons pass and
  // their summed score is lower, swap
  template <bool isMC, class TBCs, class TCollisions, class TTracks>
  void runSwapResolver(TCollisions const& collisions, TTracks const& tracks, DedupDiag* diag)
  {
    const auto indexOfKey = candidateIndexByKey();
    std::map<int64_t, std::vector<int>> storedOfColl;
    for (const auto& key : stored_fullv0Ids) {
      if (const auto it = indexOfKey.find(key); it != indexOfKey.end()) {
        storedOfColl[key.collisionId].push_back(it->second);
      }
    }
    auto fitPair = [&](int64_t colId, int64_t posId, int64_t eleId, PairFit* probe, bool original = false) {
      const auto pos = tracks.rawIteratorAt(posId);
      SourcePV srcPv;
      const SourcePV* pv = nullptr;
      if constexpr (requires { pos.sourcePvX(); }) {
        if (original && mixed.refitOriginalsAtSourceVertex) {
          srcPv = SourcePV{pos.sourcePvX(), pos.sourcePvY(), pos.sourcePvZ(), pos.sourcePvCovXX(), pos.sourcePvCovXY(), pos.sourcePvCovXZ(), pos.sourcePvCovYY(), pos.sourcePvCovYZ(), pos.sourcePvCovZZ(), pos.sourcePvChi2(), pos.sourcePvNContrib()};
          pv = &srcPv;
        }
      }
      const auto ele = tracks.rawIteratorAt(eleId);
      const bool finderLike = mixed.uniformPairFit ? !(pos.hasITS() && ele.hasITS()) : !recovery.finderUniformFit.value; // uniformPairFit (run 2): options by the legs
      fitV0Candidate<isMC, TBCs, TCollisions, TTracks>(collisions.rawIteratorAt(colId), pos, ele, -1, false,
                                                       FitOptions{.nonSVertexer = finderLike, .relaxAP = finderLike, .purityGate = mixed.uniformPairFit && finderLike && recovery.runFinder.value, .probe = probe, .pv = pv});
    };
    auto probeFit = [&](int64_t colId, int64_t posId, int64_t eleId, bool original = false) {
      PairFit out;
      fitPair(colId, posId, eleId, &out, original);
      return out;
    };
    auto crossedPrefilter = [&](int64_t posId, int64_t eleId) {
      const auto pos = tracks.rawIteratorAt(posId);
      const auto ele = tracks.rawIteratorAt(eleId);
      if (std::fabs(pos.tgl() - ele.tgl()) > recovery.finderMaxDeltaTgl) {
        return false;
      }
      const o2::track::TrackAuxPar cp(getTrackPar(pos), d_bz);
      const o2::track::TrackAuxPar ce(getTrackPar(ele), d_bz);
      return std::fabs(circleGap(cp, ce)) < recovery.finderMaxCircleGap;
    };
    std::map<int, PairFit> storedFit;
    auto storedFitOf = [&](int idx) -> const PairFit& {
      auto it = storedFit.find(idx);
      if (it == storedFit.end()) {
        const auto& c = vecV0Dedup[idx];
        it = storedFit.emplace(idx, probeFit(c.colID, c.posID, c.eleID, true)).first;
      }
      return it->second;
    };

    //  test all close pairs of stored photons
    struct Test {
      SwapRecord rec;
      float gain = -999.f;
    };
    std::vector<SwapRecord> swaps;
    for (const auto& [colId, members] : storedOfColl) {
      std::vector<Test> tests;
      for (size_t a = 0; a < members.size(); ++a) {
        for (size_t b = a + 1; b < members.size(); ++b) {
          const auto& ci = vecV0Dedup[members[a]];
          const auto& cj = vecV0Dedup[members[b]];
          if (ci.posID == cj.posID || ci.eleID == cj.eleID) {
            continue;
          }
          const float q = qInv(ci.pxAtPV, ci.pyAtPV, ci.pzAtPV, cj.pxAtPV, cj.pyAtPV, cj.pzAtPV);
          if (q > recovery.swapMaxQ || std::hypot(ci.svX - cj.svX, ci.svY - cj.svY, ci.svZ - cj.svZ) > recovery.swapMaxConvDist) {
            continue;
          }
          const PairFit storedI = storedFitOf(members[a]);
          const PairFit storedJ = storedFitOf(members[b]);
          const bool preI = crossedPrefilter(ci.posID, cj.eleID);
          const PairFit crossedI = preI ? probeFit(colId, ci.posID, cj.eleID) : PairFit{};
          const int lineI = preI && !crossedI.ok ? mLastFitReject : 0;
          const bool preJ = crossedPrefilter(cj.posID, ci.eleID);
          const PairFit crossedJ = preJ ? probeFit(colId, cj.posID, ci.eleID) : PairFit{};
          const int lineJ = preJ && !crossedJ.ok ? mLastFitReject : 0;
          Test t;
          t.rec = {keyOf(ci), keyOf(cj), q, storedI.ok && storedJ.ok && crossedI.ok && crossedJ.ok, false};
          t.rec.fail = (storedI.ok ? 0 : 1) | (storedJ.ok ? 0 : 2) | (preI ? 0 : 4) | (preJ ? 0 : 8) | (preI && !crossedI.ok ? 16 : 0) | (preJ && !crossedJ.ok ? 32 : 0); // o2-linter: disable=magic-number (bit mask)
          t.rec.line = lineI > 0 ? lineI : lineJ;
          if (t.rec.complete) { // gain > 0: the crossed pairing has the lower (better) summed dedup score
            t.gain = (storedI.score + storedJ.score) - (crossedI.score + crossedJ.score);
            t.rec.gain = t.gain;
          }
          tests.push_back(t);
        }
      }
      std::vector<int> order(tests.size());
      std::iota(order.begin(), order.end(), 0);
      std::sort(order.begin(), order.end(), [&](int x, int y) { return tests[x].gain > tests[y].gain; });
      std::set<CandKey> used;
      for (const int& k : order) {
        auto& t = tests[k];
        if (t.rec.complete && t.gain > 0.f && !used.contains(t.rec.keyI) && !used.contains(t.rec.keyJ)) {
          t.rec.swapped = true;
          used.insert(t.rec.keyI);
          used.insert(t.rec.keyJ);
          swaps.push_back(t.rec);
        }
        registry.fill(HIST("Swap/hTested"), 0.f);
        registry.fill(HIST("Swap/hQ"), t.rec.q, t.rec.swapped ? 1.f : 0.f);
        if (diag != nullptr) {
          diag->swaps.push_back(t.rec);
        }
      }
    }
    if (swaps.empty()) {
      return;
    }

    //  apply swaps: the crossed photons become candidates and replace the stored pair
    const std::set<CandKey> storedBefore(stored_fullv0Ids.begin(), stored_fullv0Ids.end());
    std::set<CandKey> removedKeys;
    std::vector<CandKey> addedKeys;
    for (const auto& s : swaps) {
      const int64_t colId = s.keyI.collisionId;
      const int64_t posI = s.keyI.posId, eleI = s.keyI.eleId, posJ = s.keyJ.posId, eleJ = s.keyJ.eleId;
      const CandKey a1 = CandKey{.v0Id = -1, .collisionId = colId, .posId = posI, .eleId = eleJ};
      const CandKey a2 = CandKey{.v0Id = -1, .collisionId = colId, .posId = posJ, .eleId = eleI};
      if (storedBefore.contains(a1) || storedBefore.contains(a2)) {
        continue;
      }
      if (!pca_map.contains(a1)) {
        fitPair(colId, posI, eleJ, nullptr);
      }
      if (!pca_map.contains(a2)) {
        fitPair(colId, posJ, eleI, nullptr);
      }
      if (!pca_map.contains(a1) || !pca_map.contains(a2)) {
        continue;
      }
      removedKeys.insert(s.keyI);
      removedKeys.insert(s.keyJ);
      addedKeys.push_back(a1);
      addedKeys.push_back(a2);
      mRelaxedKeys.insert(a1); // refitted for the table without the AP cut, as in the resolver fit
      mRelaxedKeys.insert(a2);
      registry.fill(HIST("Swap/hTested"), 1.f);
      if (diag != nullptr) {
        diag->blockersByKey[s.keyI] = {a1, a2};
        diag->blockersByKey[s.keyJ] = {a1, a2};
      }
    }
    std::vector<CandKey> kept;
    kept.reserve(stored_fullv0Ids.size());
    for (const auto& key : stored_fullv0Ids) {
      if (!removedKeys.contains(key)) {
        kept.push_back(key);
      }
    }
    kept.insert(kept.end(), addedKeys.begin(), addedKeys.end());
    setStored(std::move(kept));
  }

  //  deduplication
  Preslice<o2::aod::V0s> perCollision = o2::aod::v0::collisionId;
  std::map<CandKey, float> pca_map;                      // (v0.globalIndex(), collision.globalIndex(), pos.globalIndex(), ele.globalIndex()) -> pca
  std::map<CandKey, float> cospa_map;                    // (v0.globalIndex(), collision.globalIndex(), pos.globalIndex(), ele.globalIndex()) -> cospa
  int mLastFitReject = 0;                                // source line of the last rejection in fitV0Candidate, 0 = accepted (diagnostic)
  DedupDiag* mDiag = nullptr;                            // run 1 with MC diagnostics: where the finder and fillV0Table record the genuine leg pairs
  std::unordered_map<int64_t, int64_t> mLegMother;       // MC: track -> first mother of its MC electron (genuine leg pairs share it)
  std::unordered_map<int64_t, float> mOccOfSourceEvent;  // run 2: mixer event uid -> occupancy (from its same-event artificial event)
  std::vector<V0CandidateHelper> vecV0Dedup;             // vector with all V0Candidates that is used to sort them by score (see struct V0CandidateHelper for more details of content)
  std::vector<std::pair<int64_t, int64_t>> stored_v0Ids; // (pos.globalIndex(), ele.globalIndex())
  std::vector<CandKey> stored_fullv0Ids;                 // (v0.globalIndex(), collision.globalIndex(), pos.globalIndex(), ele.globalIndex())
  std::unordered_map<int64_t, int> nv0_map;              // map collisionId -> nv0

  // mode 0: a candidate is kept if no other candidate sharing a leg has a smaller PCA and no copy in another collision has a larger cosPA
  void selectPairwise(DedupDiag* diag)
  {
    for (const auto& [key, value] : pca_map) {
      auto v0Id = key.v0Id;
      auto collisionId = key.collisionId;
      auto posId = key.posId;
      auto eleId = key.eleId;
      float v0pca = value;
      float cospa = cospa_map[key];
      bool is_closest_v0 = true;
      bool is_most_aligned_v0 = true;

      for (const auto& [key_tmp, value_tmp] : pca_map) {
        auto collisionId_tmp = key_tmp.collisionId;
        auto posId_tmp = key_tmp.posId;
        auto eleId_tmp = key_tmp.eleId;
        float v0pca_tmp = value_tmp;
        float cospa_tmp = cospa_map[key_tmp];

        if (key == key_tmp) { // skip exactly the same candidate
          continue;
        }

        if (collisionId != collisionId_tmp && eleId == eleId_tmp && posId == posId_tmp && cospa < cospa_tmp) {
          if (diag != nullptr) {
            diag->blockersByKey[key].push_back(key_tmp);
          }
          is_most_aligned_v0 = false;
          break;
        }
        if ((eleId == eleId_tmp || posId == posId_tmp) && v0pca > v0pca_tmp) {
          if (diag != nullptr) {
            diag->blockersByKey[key].push_back(key_tmp);
          }
          is_closest_v0 = false;
          break;
        }
      } // end of pca_map tmp loop

      bool is_stored = std::find(stored_v0Ids.begin(), stored_v0Ids.end(), std::make_pair(posId, eleId)) != stored_v0Ids.end();
      if (is_closest_v0 && is_most_aligned_v0 && !is_stored) {
        // auto v0 = v0s.rawIteratorAt(v0Id);
        // auto collision = collisions.rawIteratorAt(collisionId);
        // auto pos = tracks.rawIteratorAt(posId);
        // auto ele = tracks.rawIteratorAt(eleId);
        // LOGF(info, "!accept! | collision id = %d | v0id1 = %d , posid1 = %d , eleid1 = %d , pca1 = %f , cospa = %f", collisionId, v0Id, posId, eleId, v0pca, cospa);

        // fillV0Table<isMC, TCollisions, TTracks>(v0, true);
        stored_v0Ids.emplace_back(std::make_pair(posId, eleId));
        stored_fullv0Ids.push_back(CandKey{.v0Id = v0Id, .collisionId = collisionId, .posId = posId, .eleId = eleId});
        nv0_map[collisionId]++;
      }
    } // end of pca_map loop
    if (recovery.runFinder || recovery.runSwapResolver) {
      std::stable_sort(stored_fullv0Ids.begin(), stored_fullv0Ids.end(), [](const auto& a, const auto& b) { return a.collisionId < b.collisionId; });
    }
    // LOGF(info, "pca_map.size() = %d", pca_map.size());
  }

  // mode 3: every candidate is kept, including all duplicates
  void keepAll()
  {
    // Keep-all: Every candidate that survived the quality cuts is stored, including all duplicates
    stored_v0Ids.clear();
    stored_fullv0Ids.clear();
    nv0_map.clear();
    for (const auto& v0Cand : vecV0Dedup) {
      stored_v0Ids.emplace_back(v0Cand.posID, v0Cand.eleID);
      stored_fullv0Ids.push_back(keyOf(v0Cand));
      nv0_map[v0Cand.colID]++;
    }
    std::sort(stored_fullv0Ids.begin(), stored_fullv0Ids.end(),
              [](const auto& a, const auto& b) {
                return a.collisionId < b.collisionId;
              });
  }

  // mode 2: candidates sharing a leg form a conflict group; per group the best leg-disjoint subset is kept
  void selectByConflictGroup(DedupDiag* diag)
  {
    const int maxEnum = std::min(dedup.dedupMaxGroupSize.value, 16); // o2-linter: disable=magic-number (subset check of photons)
    const int nCand = static_cast<int>(vecV0Dedup.size());
    stored_v0Ids.clear();
    stored_fullv0Ids.clear();
    nv0_map.clear();
    mGroupOfKey.clear();
    mPermuteResults.clear();
    if constexpr (kMixedInstance) {
      mGroupSizeOfCand.assign(nCand, 0);
    }
    int nGreedyFallback = 0;

    //  conflict groups: candidates connected through a shared leg
    std::vector<int> parent(nCand);
    for (int i = 0; i < nCand; ++i) {
      parent[i] = i;
    }
    auto findRoot = [&parent](int x) {
      while (parent[x] != x) {
        parent[x] = parent[parent[x]]; // path halving
        x = parent[x];
      }
      return x;
    };
    std::unordered_map<int, int> firstCandOfTrack;
    for (int i = 0; i < nCand; ++i) {
      for (const int& trackId : {vecV0Dedup[i].posID, vecV0Dedup[i].eleID}) {
        auto [it, inserted] = firstCandOfTrack.try_emplace(trackId, i);
        if (!inserted) {
          parent[findRoot(i)] = findRoot(it->second);
        }
      }
    }
    std::map<int, std::vector<int>> groups;
    for (int i = 0; i < nCand; ++i) {
      groups[findRoot(i)].push_back(i);
    }

    //  decide each group
    for (auto& [root, members] : groups) { // o2-linter: disable=const-ref-in-for-loop (members is sorted in place below)
      const int nMembers = static_cast<int>(members.size());
      std::vector<char> selected(nMembers, 0);

      // Fix the order inside the group: best score first
      std::sort(members.begin(), members.end(), [this](int a, int b) {
        if (vecV0Dedup[a].score != vecV0Dedup[b].score) {
          return vecV0Dedup[a].score < vecV0Dedup[b].score;
        }
        if (vecV0Dedup[a].v0ID != vecV0Dedup[b].v0ID) {
          return vecV0Dedup[a].v0ID < vecV0Dedup[b].v0ID;
        }
        return std::tie(vecV0Dedup[a].posID, vecV0Dedup[a].eleID, vecV0Dedup[a].colID) < std::tie(vecV0Dedup[b].posID, vecV0Dedup[b].eleID, vecV0Dedup[b].colID);
      });

      if (nMembers == 1) {
        selected[0] = 1;
      } else if (nMembers <= maxEnum) {

        std::unordered_map<int, int> legBit;
        std::vector<uint64_t> legMask(nMembers, 0);
        for (int k = 0; k < nMembers; ++k) {
          for (const int& trackId : {vecV0Dedup[members[k]].posID, vecV0Dedup[members[k]].eleID}) {
            const int bit = legBit.try_emplace(trackId, static_cast<int>(legBit.size())).first->second;
            legMask[k] |= (1ull << bit); // at most 2 * 16 = 32 legs, fits into 64 bit
          }
        }
        std::vector<float> scores(nMembers);
        for (int k = 0; k < nMembers; ++k) {
          scores[k] = vecV0Dedup[members[k]].score;
        }
        const auto sol = solveConflictGroup(legMask, scores, dedup.groupMatchingScoreThreshold.value, diag != nullptr);
        for (int k = 0; k < nMembers; ++k) {
          selected[k] = ((sol.bestMask & (1u << k)) != 0) ? 1 : 0;
          if (diag != nullptr && selected[k] == 0 && sol.bestWithCount[k] >= 0) {
            const auto& cand = vecV0Dedup[members[k]];
            diag->lossMargin[keyOf(cand)] =
              std::make_pair(sol.bestCount - sol.bestWithCount[k], sol.bestWithScore[k] - sol.bestScore);
          }
        }
        if constexpr (kMixedInstance) { // permutation test: the same group in reversed order must select the same candidates
          const std::vector<uint64_t> legMaskRev(legMask.rbegin(), legMask.rend());
          const std::vector<float> scoresRev(scores.rbegin(), scores.rend());
          const auto solRev = solveConflictGroup(legMaskRev, scoresRev, dedup.groupMatchingScoreThreshold.value, false);
          bool changed = false;
          for (int k = 0; k < nMembers; ++k) {
            changed = changed || (((solRev.bestMask & (1u << (nMembers - 1 - k))) != 0) != (selected[k] != 0));
          }
          mPermuteResults.push_back({vecV0Dedup[members[0]].colID, nMembers, changed, false});
        }
      } else {
        // Too large to enumerate: greedy by score
        ++nGreedyFallback;
        if constexpr (kMixedInstance) {
          mPermuteResults.push_back({vecV0Dedup[members[0]].colID, nMembers, false, true});
        }
        std::vector<int> usedTracks;
        usedTracks.reserve(2 * nMembers);
        for (int k = 0; k < nMembers; ++k) {
          const auto& cand = vecV0Dedup[members[k]];
          if (std::find(usedTracks.begin(), usedTracks.end(), cand.posID) != usedTracks.end() ||
              std::find(usedTracks.begin(), usedTracks.end(), cand.eleID) != usedTracks.end()) {
            continue;
          }
          usedTracks.push_back(cand.posID);
          usedTracks.push_back(cand.eleID);
          selected[k] = 1;
        }
      }

      // run 1 (MC diagnostic of the swap decision): half-swap rule: if a candidate is worse than two others that share its legs, and the crossed candidate is not present, replace it
      if (dedup.dedupHalfSwapRule && nMembers >= 3) {
        int nReplaced = 0;
        for (int t = 0; t < nMembers; ++t) {
          if (selected[t] != 0) {
            continue;
          }
          const auto& cT = vecV0Dedup[members[t]];
          int w1 = -1, w2 = -1;
          for (int k = 0; k < nMembers; ++k) {
            if (selected[k] == 0 || k == t) {
              continue;
            }
            const auto& cW = vecV0Dedup[members[k]];
            if (cW.posID == cT.posID) {
              w1 = k;
            }
            if (cW.eleID == cT.eleID) {
              w2 = k;
            }
          }
          if (w1 < 0 || w2 < 0 || w1 == w2) {
            continue;
          }
          const auto& cW1 = vecV0Dedup[members[w1]];
          const auto& cW2 = vecV0Dedup[members[w2]];
          if (!(cT.score < std::min(cW1.score, cW2.score))) {
            continue;
          }
          bool complementIsCandidate = false;
          for (int k = 0; k < nMembers; ++k) {
            const auto& cC = vecV0Dedup[members[k]];
            if (cC.posID == cW2.posID && cC.eleID == cW1.eleID) {
              complementIsCandidate = true;
              break;
            }
          }
          if (complementIsCandidate) {
            continue;
          }
          selected[w1] = 0;
          selected[w2] = 0;
          selected[t] = 1;
          ++nReplaced;
        }
        registry.fill(HIST("Dedup/hHalfSwapRule"), static_cast<float>(nReplaced));
      }

      std::array<int, 4> groupInfo{root, nMembers, 0, static_cast<int>(std::count(selected.begin(), selected.end(), 1))};
      if constexpr (kMixedInstance) {
        std::set<int> legsOfGroup;
        for (const int& m : members) {
          legsOfGroup.insert(vecV0Dedup[m].posID);
          legsOfGroup.insert(vecV0Dedup[m].eleID);
          mGroupSizeOfCand[m] = nMembers;
        }
        groupInfo[2] = static_cast<int>(legsOfGroup.size());
      }
      //  store the winners, book the blockers of the losers
      for (int k = 0; k < nMembers; ++k) {
        const auto& cand = vecV0Dedup[members[k]];
        if (selected[k] != 0) {
          stored_v0Ids.emplace_back(cand.posID, cand.eleID);
          stored_fullv0Ids.push_back(keyOf(cand));
          nv0_map[cand.colID]++;
          if constexpr (kMixedInstance) {
            mGroupOfKey[keyOf(cand)] = groupInfo;
          }
          continue;
        }
        if (diag == nullptr) {
          continue;
        }
        // all winners that took one leg
        const CandKey victimKey = keyOf(cand);
        for (int j = 0; j < nMembers; ++j) {
          if (selected[j] == 0) {
            continue;
          }
          const auto& winner = vecV0Dedup[members[j]];
          if (winner.posID == cand.posID || winner.eleID == cand.eleID) {
            diag->blockersByKey[victimKey].push_back(keyOf(winner));
          }
        }
      }
    }

    if (nGreedyFallback > 0) {
      LOGF(warning, "group matching: %d conflict group(s) larger than dedupMaxGroupSize = %d were solved greedily", nGreedyFallback, maxEnum);
    }

    // the table expects candidates ordered by collision
    std::sort(stored_fullv0Ids.begin(), stored_fullv0Ids.end(),
              [](const auto& a, const auto& b) {
                return a.collisionId < b.collisionId;
              });
  }

  // mode 1: candidates sorted by score (lower = better), a candidate is kept if both legs are still free
  void selectGreedy(const size_t nTracks, DedupDiag* diag)
  {
    // Sort best candidates first, depending on score
    std::sort(vecV0Dedup.begin(), vecV0Dedup.end());
    // container to keep track of which electron and positron tracks have already been used
    std::vector<bool> usedLegs(nTracks, false);
    std::unordered_map<int, CandKey> ownerOfLeg; // diagnostics only
    stored_v0Ids.clear();
    stored_fullv0Ids.clear();
    nv0_map.clear();
    for (const auto& v0Cand : vecV0Dedup) {
      const CandKey candKey = keyOf(v0Cand);
      if (usedLegs[v0Cand.posID] || usedLegs[v0Cand.eleID]) {
        if (diag != nullptr) {
          // both legs, not just the first one found: a cross-leg fake is usually blocked by
          // two different winners, and that is exactly the case being measured here
          auto& blockers = diag->blockersByKey[candKey];
          for (const int& legId : {v0Cand.posID, v0Cand.eleID}) {
            const auto itOwner = ownerOfLeg.find(legId);
            if (itOwner == ownerOfLeg.end()) {
              continue;
            }
            if (std::find(blockers.begin(), blockers.end(), itOwner->second) == blockers.end()) {
              blockers.push_back(itOwner->second); // the same winner may hold both legs
            }
          }
        }
        continue;
      }
      usedLegs[v0Cand.posID] = true;
      usedLegs[v0Cand.eleID] = true;
      if (diag != nullptr) {
        ownerOfLeg[v0Cand.posID] = candKey;
        ownerOfLeg[v0Cand.eleID] = candKey;
      }
      stored_v0Ids.emplace_back(v0Cand.posID, v0Cand.eleID);
      stored_fullv0Ids.push_back(keyOf(v0Cand));
      nv0_map[v0Cand.colID]++;
    }

    // sort accepted candidates by collision ID to be compatible with table
    std::sort(stored_fullv0Ids.begin(), stored_fullv0Ids.end(),
              [](const auto& a, const auto& b) {
                return a.collisionId < b.collisionId;
              });
  }

  template <bool isMC, bool isTriggerAnalysis, bool enableFilter, typename TCollisions, typename TV0s, typename TTracks, typename TBCs>
  void build(TCollisions const& collisions, TV0s const& v0s, TTracks const& tracks, TBCs const&, DedupDiag* diag = nullptr)
  {
    mRelaxedKeys.clear();
    mLegDcaCache.clear();
    mDiag = diag;
    mLegMother.clear();
    if constexpr (isMC) {
      if (diag != nullptr) { // genuine leg pairs = opposite-sign MC electrons with the same first mother (MCDedup/Inventory/hSwapOrigin)
        for (const auto& track : tracks) {
          if (!track.has_mcParticle()) {
            continue;
          }
          const auto mc = track.template mcParticle_as<o2::aod::McParticles>();
          if (std::abs(mc.pdgCode()) == kElectron && mc.has_mothers()) {
            mLegMother[track.globalIndex()] = mc.mothersIds()[0];
          }
        }
      }
    }
    for (const auto& collision : collisions) {
      if constexpr (isMC) {
        if (!collision.has_mcCollision()) {
          continue;
        }
      }

      if (!collision.isSelected()) {
        continue;
      }

      // if constexpr (isTriggerAnalysis) {
      //   if (collision.triggerMask_raw() == 0) {
      //     continue;
      //   }
      // }

      nv0_map[collision.globalIndex()] = 0;

      const auto& bc = collision.template foundBC_as<MyBCs>();
      initCCDB(bc);
      registry.fill(HIST("hCollisionCounter"), 1);

      updateCCDB(bc); // delay update until is needed

      vecV0Dedup.reserve(30000); // rough estimate for number of V0s per DF
      const auto& v0s_per_coll = v0s.sliceBy(perCollision, collision.globalIndex());
      // LOGF(info, "n v0 = %d", v0s_per_coll.size());
      for (const auto& v0 : v0s_per_coll) {
        // LOGF(info, "collision.globalIndex() = %d, v0.globalIndex() = %d, v0.posTrackId() = %d, v0.negTrackId() = %d", collision.globalIndex(), v0.globalIndex(), v0.posTrackId() , v0.negTrackId());
        fillV0Table<isMC, TBCs, TCollisions, TTracks>(v0, false);
      } // end of v0 loop
    } // end of collision loop

    if (recovery.runFinder) {
      runPhotonFinder<isMC, TBCs, TCollisions, TTracks>(collisions, tracks);
    }

    stored_v0Ids.reserve(pca_map.size());     // number of photon candidates per DF
    stored_fullv0Ids.reserve(pca_map.size()); // number of photon candidates per DF

    switch (dedup.deduplicationMode.value) {
      case V0DeduplicationMode::Pairwise:
        selectPairwise(diag);
        break;
      case V0DeduplicationMode::KeepAll:
        keepAll();
        break;
      case V0DeduplicationMode::GroupMatching:
        selectByConflictGroup(diag);
        break;
      default:
        selectGreedy(tracks.size(), diag);
        break;
    }

    if (dedup.resolveSplitLegs) {
      resolveSplitLegConflicts(tracks, diag);
    }
    if (recovery.runSwapResolver) {
      runSwapResolver<isMC, TBCs, TCollisions, TTracks>(collisions, tracks, diag);
    }

    if (diag != nullptr) {
      diag->snapshot = vecV0Dedup;
      diag->stored = stored_fullv0Ids;
    }

    std::set<int64_t> skimRejected;
    if (skim.minNv0PerCollision > 0 || skim.maxMinvPairPerCollision > 0.f || skim.minMinvPairPerCollision > 0.f) { // o2-linter: disable=magic-number (check if configurable is enabled)
      const std::set<CandKey> storedSetSkim(stored_fullv0Ids.begin(), stored_fullv0Ids.end());
      std::unordered_map<int64_t, std::vector<const V0CandidateHelper*>> candsPerColl;
      for (const auto& c : vecV0Dedup) {
        if (storedSetSkim.contains(keyOf(c))) {
          candsPerColl[c.colID].push_back(&c);
        }
      }
      for (const auto& [colId, cands] : candsPerColl) {
        if (skim.minNv0PerCollision > 0 && static_cast<int>(cands.size()) < skim.minNv0PerCollision) {
          skimRejected.insert(colId);
          continue;
        }
        if (skim.maxMinvPairPerCollision > 0.f || skim.minMinvPairPerCollision > 0.f) { // o2-linter: disable=magic-number (check if configurable is enabled)
          bool hasPairOfInterest = false;
          const float minMinv = skim.minMinvPairPerCollision > 0.f ? skim.minMinvPairPerCollision.value : 0.f;
          const float maxMinv = skim.maxMinvPairPerCollision > 0.f ? skim.maxMinvPairPerCollision.value : std::numeric_limits<float>::max();
          for (size_t i = 0; i < cands.size() && !hasPairOfInterest; ++i) {
            for (size_t j = i + 1; j < cands.size(); ++j) {
              ROOT::Math::PxPyPzMVector v1(cands[i]->pxAtPV, cands[i]->pyAtPV, cands[i]->pzAtPV, 0.f);
              ROOT::Math::PxPyPzMVector v2(cands[j]->pxAtPV, cands[j]->pyAtPV, cands[j]->pzAtPV, 0.f);
              const double massPair = (v1 + v2).M();
              if (massPair > minMinv && massPair < maxMinv) {
                hasPairOfInterest = true;
                break;
              }
            }
          }
          if (!hasPairOfInterest) {
            skimRejected.insert(colId);
          }
        }
      }
    }

    for (const auto& fullv0Id : stored_fullv0Ids) {
      auto v0Id = fullv0Id.v0Id;
      const auto collisionId = fullv0Id.collisionId;
      if (skimRejected.contains(collisionId)) {
        continue; // collision failed the skim - do not store its photons/legs
      }
      // auto posId = fullv0Id.posId;
      // auto eleId = fullv0Id.eleId;
      // LOGF(info, "!accept! | collision id = %d | v0id = %d , posid = %d , eleid = %d", collisionId, v0Id, posId, eleId);

      if (v0Id >= 0) {
        auto v0 = v0s.rawIteratorAt(v0Id);
        if constexpr (enableFilter) {
          auto collision_tmp = v0.template collision_as<TCollisions>();                                               // collision where this v0 belongs.
          if (!(collision_tmp.neeuls() >= 1 || collision_tmp.neeuls() + nv0_map[collision_tmp.globalIndex()] >= 2)) { // o2-linter: disable=magic-number (nonsensical default value)
            continue;
          }
          // LOGF(info, "collision_tmp.globalIndex() = %d, collision_tmp.neeuls() = %d, nv0_map = %d", collision_tmp.globalIndex(), collision_tmp.neeuls(), nv0_map[collision_tmp.globalIndex()]);
        }

        fillV0Table<isMC, TBCs, TCollisions, TTracks>(v0, true);
        continue;
      }
      const auto collision = collisions.rawIteratorAt(collisionId);
      if constexpr (enableFilter) {
        if (!(collision.neeuls() >= 1 || collision.neeuls() + nv0_map[collisionId] >= 2)) { // o2-linter: disable=magic-number (nonsensical default value)
          continue;
        }
      }
      fitV0Candidate<isMC, TBCs, TCollisions, TTracks>(collision, tracks.rawIteratorAt(fullv0Id.posId), tracks.rawIteratorAt(fullv0Id.eleId), v0Id, true,
                                                       FitOptions{.nonSVertexer = !recovery.finderUniformFit.value, .relaxAP = mRelaxedKeys.contains(fullv0Id) && !recovery.finderUniformFit.value});
    } // end of fullv0Id loop
    if (fillFreeLegs || fillFreeLegIds) {
      writeFreeLegs<isMC, TBCs, TCollisions, TTracks>(collisions, tracks, skimRejected);
    }

    // for (const auto& collision : collisions) {
    //   if constexpr (isMC) {
    //     if (!collision.has_mcCollision()) {
    //       continue;
    //     }
    //   }
    //   // events_ngpcm(nv0_map[collision.globalIndex()]);
    // } // end of collision loop

    pca_map.clear();
    cospa_map.clear();
    nv0_map.clear();
    stored_v0Ids.clear();
    stored_v0Ids.shrink_to_fit();
    stored_fullv0Ids.clear();
    stored_fullv0Ids.shrink_to_fit();
    vecV0Dedup.clear();
  } // end of build

  // MC photon pair information for the deduplication diagnostics (run 1) and the mixer (run 2)
  struct MixPairPhoton {
    float px = 0.f, py = 0.f, pz = 0.f;
    std::pair<int64_t, int> posKey{-1, -1}, eleKey{-1, -1};
    bool posTrue = false, eleTrue = false;
    int64_t posRow = -1, eleRow = -1;      // PCMMixLegs rows
    bool original = true;                  // leg pairing of a run-1 photon (false: crossed pairing built in run 2)
    std::array<int, 4> group{-1, 0, 0, 0}; // run-2 conflict group (id, candidates, legs, stored photons), id -1: unknown
    float vx = 0.f, vy = 0.f, vz = 0.f;    // conversion point (run 1 before, run 2 after)
    struct LegMC {                         // MC of one leg for the fake-pair classes (ids within the DF of the source event)
      int64_t evt = -1;
      int pid = -1, mid = -1, gmid = -1, pdg = 0, mpdg = 0, trk = -1, anc = -1;
      int ggid = -1, root = -1, rootPdg = 0; // mother of the grandmother, root of the decay chain of the photon
      bool its = false;
      float phi = 0.f, tgl = 0.f, sharedFrac = 0.f; // track direction and TPC shared-cluster fraction (data-computable)
      float vx = 0.f, vy = 0.f, vz = 0.f;           // run-1 conversion point of the photon the leg belonged to
    };
    LegMC posMC, eleMC;
  };
  template <typename TLeg>
  static void setMixPairLegs(MixPairPhoton& p, TLeg const& pos, TLeg const& ele)
  {
    p.posKey = {pos.sourceEventUid(), pos.mcMotherId()};
    p.eleKey = {ele.sourceEventUid(), ele.mcMotherId()};
    p.posTrue = pos.mcMotherId() >= 0 && pos.mcMotherPdg() == kGamma;
    p.eleTrue = ele.mcMotherId() >= 0 && ele.mcMotherPdg() == kGamma;
    using LegMC = typename MixPairPhoton::LegMC;
    auto legMC = [](TLeg const& l) {
      return LegMC{l.sourceEventUid(), l.mcParticleId(), l.mcMotherId(), l.mcGrandMotherId(), l.mcPdg(), l.mcMotherPdg(), l.origTrackId(), l.mcAncestorId(),
                   l.mcGGMotherId(), l.mcRootId(), l.mcRootPdg(), static_cast<bool>(l.hasITS()), l.phi(), l.tgl(), l.tpcFractionSharedCls(), l.sourceVx(), l.sourceVy(), l.sourceVz()};
    };
    p.posMC = legMC(pos);
    p.eleMC = legMC(ele);
  }

  static int mixCommonOrigin(MixPairPhoton::LegMC const& x, MixPairPhoton::LegMC const& y)
  {
    if (x.mid < 0 || y.mid < 0 || (x.evt == y.evt && x.mid == y.mid)) {
      return -1;
    }
    if (x.evt != y.evt) {
      return 0;
    }
    if (x.gmid >= 0 && x.gmid == y.gmid) {
      return 3;
    }
    if (x.ggid >= 0 && x.ggid == y.ggid) {
      return 2;
    }
    return (x.root >= 0 && x.root == y.root) ? 1 : 0;
  }
  // common origin of the two MC photons whose legs form a photon pair (true pair, full or half swap): first leg vs the first leg of
  // another MC photon among the remaining three
  static int mixPairOrigin(MixPairPhoton const& a, MixPairPhoton const& b)
  {
    const auto& x = a.posMC;
    for (const auto* y : {&a.eleMC, &b.posMC, &b.eleMC}) {
      if (y->mid >= 0 && x.mid >= 0 && !(x.evt == y->evt && x.mid == y->mid)) {
        return mixCommonOrigin(x, *y);
      }
    }
    return -1;
  }

  static bool mixRelated(MixPairPhoton const& a, MixPairPhoton const& b)
  {
    using L = MixPairPhoton::LegMC;
    auto radiatedBy = [](L const& x, L const& y) { return x.evt == y.evt && x.gmid >= 0 && x.gmid == y.pid; };
    auto siblings = [](L const& x, L const& y) { return x.evt == y.evt && x.gmid >= 0 && x.gmid == y.gmid && x.mid != y.mid; };
    for (const auto* x : {&a.posMC, &a.eleMC}) {
      for (const auto* y : {&b.posMC, &b.eleMC}) {
        if (radiatedBy(*x, *y) || radiatedBy(*y, *x) || siblings(*x, *y)) {
          return true;
        }
      }
    }
    return false;
  }

  static int mixFakeClass(MixPairPhoton const& a, MixPairPhoton const& b)
  {
    using L = MixPairPhoton::LegMC;
    const std::array<const L*, 2> la{&a.posMC, &a.eleMC}, lb{&b.posMC, &b.eleMC};
    auto samePart = [](L const& x, L const& y) { return x.evt == y.evt && x.pid >= 0 && x.pid == y.pid; };
    auto sameMother = [](L const& x, L const& y) { return x.evt == y.evt && x.mid >= 0 && x.mid == y.mid; };
    auto daughterOf = [](L const& x, L const& y) { return x.evt == y.evt && x.mid >= 0 && x.mid == y.pid; };
    auto truePhoton = [&](MixPairPhoton const& p) { return p.posMC.mpdg == kGamma && p.eleMC.mpdg == kGamma && sameMother(p.posMC, p.eleMC); };
    for (const auto* x : la) {
      for (const auto* y : lb) {
        if (x->evt == y->evt && x->trk >= 0 && x->trk == y->trk) {
          return 12;
        }
      }
    }
    if (truePhoton(a) && truePhoton(b) && sameMother(a.posMC, b.posMC)) {
      return 1;
    }
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) {
        if (samePart(*la[i], *lb[j])) {
          return (la[1 - i]->mpdg == kGamma && lb[1 - j]->mpdg == kGamma) ? 2 : 3;
        }
      }
    }
    for (const auto* x : la) {
      for (const auto* y : lb) {
        if (x->mpdg == kGamma && sameMother(*x, *y)) {
          return 4;
        }
      }
    }
    for (const auto* x : la) {
      for (const auto* y : lb) {
        if (daughterOf(*x, *y) || daughterOf(*y, *x)) {
          return 5;
        }
      }
    }
    if (mixRelated(a, b)) {
      return 5;
    }
    for (const auto* x : la) {
      for (const auto* y : lb) {
        if (sameMother(*x, *y)) {
          return 6;
        }
      }
    }
    const bool aFake = !truePhoton(a), bFake = !truePhoton(b);
    if (aFake && bFake) {
      return 7;
    }
    const auto& f = aFake ? a : b;
    if (sameMother(f.posMC, f.eleMC)) {
      return 8;
    }
    if (std::abs(f.posMC.pdg) != kElectron || std::abs(f.eleMC.pdg) != kElectron) {
      return 9;
    }
    if (f.posMC.mpdg == kGamma && f.eleMC.mpdg == kGamma) {
      return 10;
    }
    return 11;
  }
  // two photons come from one shower: a leg of each has the same physical-primary ancestor (photons of one primary photon/e+-)
  static bool mixSameAncestor(MixPairPhoton const& a, MixPairPhoton const& b)
  {
    for (const auto* x : {&a.posMC, &a.eleMC}) {
      for (const auto* y : {&b.posMC, &b.eleMC}) {
        if (x->evt == y->evt && x->anc >= 0 && x->anc == y->anc) {
          return true;
        }
      }
    }
    return false;
  }

  static std::array<float, 3> mixPairGeometry(MixPairPhoton const& a, MixPairPhoton const& b)
  {
    using L = MixPairPhoton::LegMC;
    auto mid = [](L const& p, L const& e) { return std::array<float, 3>{0.5f * (p.vx + e.vx), 0.5f * (p.vy + e.vy), 0.5f * (p.vz + e.vz)}; };
    const auto ca = mid(a.posMC, a.eleMC), cb = mid(b.posMC, b.eleMC);
    auto legDist = [](L const& x, L const& y) { return std::hypot(RecoDecay::constrainAngle(x.phi - y.phi, -o2::constants::math::PI), x.tgl - y.tgl); };
    return {std::hypot(ca[0] - cb[0], ca[1] - cb[1], ca[2] - cb[2]), std::min(legDist(a.posMC, b.posMC), legDist(a.eleMC, b.eleMC)),
            std::max({a.posMC.sharedFrac, a.eleMC.sharedFrac, b.posMC.sharedFrac, b.eleMC.sharedFrac})};
  }
  // MC photon pair relationship: 0 unrelated, 1 radiated by the other, 2 siblings (same grandmother), 3 daughter of the other
  static int mixRelKind(MixPairPhoton const& a, MixPairPhoton const& b)
  {
    using L = MixPairPhoton::LegMC;
    auto radiatedBy = [](L const& x, L const& y) { return x.evt == y.evt && x.gmid >= 0 && x.gmid == y.pid; };
    auto siblings = [](L const& x, L const& y) { return x.evt == y.evt && x.gmid >= 0 && x.gmid == y.gmid && x.mid != y.mid; };
    auto daughterOf = [](L const& x, L const& y) { return x.evt == y.evt && x.mid >= 0 && x.mid == y.pid; };
    int kind = 0;
    for (const auto* x : {&a.posMC, &a.eleMC}) {
      for (const auto* y : {&b.posMC, &b.eleMC}) {
        if (radiatedBy(*x, *y) || radiatedBy(*y, *x)) {
          return 1;
        }
        if (daughterOf(*x, *y) || daughterOf(*y, *x)) {
          kind = 3;
        } else if (kind == 0 && siblings(*x, *y)) {
          kind = 2;
        }
      }
    }
    return kind;
  }
  struct BremsPointing {
    float dMin = 999.f, angle = 3.15f, rEmit = -1.f, outward = 0.f, ptB = 0.f;
  };
  template <typename TLegs>
  BremsPointing mixBremsPointing(MixPairPhoton const& a, MixPairPhoton const& b, TLegs const& legs) const
  {
    constexpr float StepCm = 0.25f, MaxArcCm = 100.f, MarginCm = 2.f;
    BremsPointing best;
    for (int order = 0; order < 2; ++order) {
      const auto& emitter = order == 0 ? a : b;
      const auto& brems = order == 0 ? b : a;
      const float pB = std::hypot(brems.px, brems.py, brems.pz);
      if (pB <= 0.f) {
        continue;
      }
      const std::array<float, 3> u{brems.px / pB, brems.py / pB, brems.pz / pB};
      const float rA = std::hypot(emitter.vx, emitter.vy), rB = std::hypot(brems.vx, brems.vy);
      for (const int64_t row : {emitter.posRow, emitter.eleRow}) {
        const auto leg = legs.rawIteratorAt(row);
        const o2::track::TrackPar par(leg.x(), leg.alpha(), std::array<float, 5>{leg.y(), leg.z(), leg.snp(), leg.tgl(), leg.signed1Pt()});
        o2::math_utils::CircleXYf_t circle;
        float sna = 0.f, csa = 0.f;
        par.getCircleParams(d_bz, circle, sna, csa);
        std::array<float, 3> x0{}, p0{};
        par.getXYZGlo(x0);
        par.getPxPyPzGlo(p0);
        if (circle.rC <= 0.f) {
          continue;
        }
        const float psi0 = std::atan2(x0[1] - circle.yC, x0[0] - circle.xC);
        const float sense = ((x0[0] - circle.xC) * p0[1] - (x0[1] - circle.yC) * p0[0]) > 0.f ? 1.f : -1.f; // direction of motion on the circle
        const float tgl = par.getTgl(), norm = std::sqrt(1.f + tgl * tgl);
        constexpr int NSteps = static_cast<int>((MaxArcCm + MaxArcCm) / StepCm);
        for (int iStep = 0; iStep <= NSteps; ++iStep) {
          const float s = -MaxArcCm + static_cast<float>(iStep) * StepCm;
          const float psi = psi0 + sense * s / circle.rC;
          const std::array<float, 3> h{circle.xC + circle.rC * std::cos(psi), circle.yC + circle.rC * std::sin(psi), x0[2] + tgl * s};
          const float rH = std::hypot(h[0], h[1]);
          if (rH < rA - MarginCm || rH > rB + MarginCm) {
            continue; // emission between the two conversions
          }
          const std::array<float, 3> d{h[0] - brems.vx, h[1] - brems.vy, h[2] - brems.vz};
          const float t = d[0] * u[0] + d[1] * u[1] + d[2] * u[2];
          if (t > 0.f) {
            continue; // in front of B's conversion point
          }
          const float dist = std::hypot(d[0] - t * u[0], d[1] - t * u[1], d[2] - t * u[2]);
          if (dist < best.dMin) {
            const std::array<float, 3> tan{-sense * std::sin(psi) / norm, sense * std::cos(psi) / norm, tgl / norm};
            best.dMin = dist;
            best.angle = std::acos(std::clamp(tan[0] * u[0] + tan[1] * u[1] + tan[2] * u[2], -1.f, 1.f));
            best.rEmit = rH;
            best.outward = rB > rA ? 1.f : 0.f;
            best.ptB = std::hypot(brems.px, brems.py);
          }
        }
      }
    }
    return best;
  }
  // skip pairs of photons that come from the same source event (MC truth) when mixing, if requested
  bool skipSameSourcePair(const float type, MixPairPhoton const& a, MixPairPhoton const& b)
  {
    if (type < 0.5f) {
      return false;
    }
    const int64_t e = a.posMC.evt;
    const bool same = e >= 0 && a.eleMC.evt == e && b.posMC.evt == e && b.eleMC.evt == e;
    return same && mixed.skipSameSourcePairs.value;
  }
  void fillSameSourcePairs(const float type, std::vector<MixPairPhoton> const& photons)
  {
    if (type < 0.5f) {
      return;
    }
    for (size_t i = 0; i < photons.size(); ++i) {
      for (size_t j = i + 1; j < photons.size(); ++j) {
        const auto& a = photons[i];
        const auto& b = photons[j];
        const int64_t e = a.posMC.evt;
        if (e >= 0 && a.eleMC.evt == e && b.posMC.evt == e && b.eleMC.evt == e) {
          registry.fill(HIST("Mixed/Pairs/hSameSource"), qInv(a.px, a.py, a.pz, b.px, b.py, b.pz), mixPairClass(a, b));
        }
      }
    }
  }
  void fillMixFakeClasses(const float type, const float stage, const int64_t colId, std::vector<MixPairPhoton> const& photons,
                          std::set<std::tuple<int64_t, int64_t, int64_t>> const& candLegs)
  {
    for (size_t i = 0; i < photons.size(); ++i) {
      for (size_t j = i + 1; j < photons.size(); ++j) {
        const auto& a = photons[i];
        const auto& b = photons[j];
        if (stage > 0.5f && skipSameSourcePair(type, a, b)) {
          continue;
        }
        const float q = qInv(a.px, a.py, a.pz, b.px, b.py, b.pz);
        if (q < 0.5f && mixPairGeometry(a, b)[2] < 0.4f) { // o2-linter: disable=magic-number (normalisation range, analysis leg cut)
          registry.fill(HIST("Mixed/Pairs/hQClean"), type, stage, q);
        }
        if (q >= 0.2f) { // o2-linter: disable=magic-number (diagnostic range)
          continue;
        }
        const auto nIts = static_cast<float>(a.posMC.its + a.eleMC.its + b.posMC.its + b.eleMC.its);
        const float amb = (candLegs.contains({colId, a.posRow, b.eleRow}) || candLegs.contains({colId, b.posRow, a.eleRow})) ? 1.f : 0.f;
        const float cls = mixPairClass(a, b);
        registry.fill(HIST("Mixed/Pairs/hFakeClass"), type, stage, 0.f, q, nIts, amb);
        registry.fill(HIST("Mixed/Pairs/hFakeClass"), type, stage, cls, q, nIts, amb);
        {
          const float rA = std::hypot(a.vx, a.vy), rB = std::hypot(b.vx, b.vy);
          const bool aSofter = std::hypot(a.px, a.py, a.pz) < std::hypot(b.px, b.py, b.pz);
          const float rSoft = aSofter ? rA : rB, rHard = aSofter ? rB : rA;
          registry.fill(HIST("Mixed/Pairs/hFakeClassOrder"), type, stage, cls, q, rSoft < rHard ? 0.f : 1.f, std::fabs(rA - rB));
        }
        if (q < 0.1f) { // o2-linter: disable=magic-number (diagnostic range)
          const auto g = mixPairGeometry(a, b);
          registry.fill(HIST("Mixed/Pairs/hFakeClassGeo"), type, stage, cls, q, amb, mixSameAncestor(a, b) ? 1.f : 0.f, g[0],
                        std::log10(std::max(g[1], 1e-6f)), g[2]); // o2-linter: disable=magic-number (log floor)
          if (stage > 0.5f) {
            registry.fill(HIST("Mixed/Pairs/hPairOrigin"), type, cls, q, static_cast<float>(mixPairOrigin(a, b)));
          }
          if (stage > 0.5f && a.group[0] >= 0 && b.group[0] >= 0) {
            // same group: its context; different groups: the larger one
            const bool sameGroup = a.group[0] == b.group[0];
            const auto& gr = (sameGroup || a.group[1] >= b.group[1]) ? a.group : b.group;
            registry.fill(HIST("Mixed/Pairs/hGroupCtx"), type, cls, q, amb, sameGroup ? 1.f : 0.f, static_cast<float>(std::min(gr[1], 20)),
                          static_cast<float>(std::min(gr[2], 33)), static_cast<float>(std::min(gr[3], 10))); // o2-linter: disable=magic-number (last bins)
          }
        }
      }
    }
  }
  // MC class of a photon pair (see fillMixFakeClasses)
  static float mixPairClass(MixPairPhoton const& a, MixPairPhoton const& b)
  {
    const int truth = mixPairTruth(a, b);
    if (truth == 3) {
      return static_cast<float>(mixFakeClass(a, b));
    }
    if (truth == 1 || truth == 2) {
      return (truth == 1 ? 13.f : 14.f) + (mixRelated(a, b) ? 3.f : 0.f);
    }
    return 15.f;
  }

  template <typename TLegs>
  void fillMixSwapBrems(const float type, const float stage, const int64_t colId, std::vector<MixPairPhoton> const& photons,
                        std::map<std::tuple<int64_t, int64_t, int64_t>, float> const& candScore, TLegs const& legs)
  {
    constexpr float MaxQ = 0.1f, MaxShared = 0.4f;
    auto score = [&](int64_t pos, int64_t ele) {
      const auto it = candScore.find({colId, pos, ele});
      return it == candScore.end() ? -1.f : it->second;
    };
    for (size_t i = 0; i < photons.size(); ++i) {
      for (size_t j = i + 1; j < photons.size(); ++j) {
        const auto& a = photons[i];
        const auto& b = photons[j];
        if (stage > 0.5f && skipSameSourcePair(type, a, b)) {
          continue;
        }
        const float q = qInv(a.px, a.py, a.pz, b.px, b.py, b.pz);
        if (q >= MaxQ) {
          continue;
        }
        const float cls = mixPairClass(a, b);
        const float clean = mixPairGeometry(a, b)[2] < MaxShared ? 1.f : 0.f;
        const auto bp = mixBremsPointing(a, b, legs);
        registry.fill(HIST("Mixed/Pairs/hBrems"), type, stage, cls, static_cast<float>(mixRelKind(a, b)), q,
                      std::log10(std::max(bp.dMin, 1e-3f)), std::log10(std::max(bp.angle, 1e-4f)), bp.rEmit, bp.outward, bp.ptB, clean); // o2-linter: disable=magic-number (log floors)
        if (stage < 0.5f || (cls != 14.f && cls != 17.f)) {                                                                              // o2-linter: disable=magic-number (swap classes)
          continue;
        }
        const float sa = score(a.posRow, a.eleRow), sb = score(b.posRow, b.eleRow);
        const float t1 = score(a.posRow, b.eleRow), t2 = score(b.posRow, a.eleRow);
        float decision = 2.f;
        if (t1 >= 0.f && t2 >= 0.f) {
          decision = (sa + sb < t1 + t2) ? 1.f : 0.f;
        }
        const float origin = (a.original && b.original) ? 1.f : ((!a.original && !b.original) ? 0.f : 2.f);
        registry.fill(HIST("Mixed/Pairs/hSwapDecision"), type, origin, decision, q, clean);
      }
    }
  }
  // truth of a photon pair: 1 two true photons, 2 the legs of two true photons crossed (full swap), 3 other with a fake, 4 no label
  static int mixPairTruth(MixPairPhoton const& a, MixPairPhoton const& b)
  {
    if (a.posKey.second < 0 || a.eleKey.second < 0 || b.posKey.second < 0 || b.eleKey.second < 0) {
      return 4;
    }
    const bool allTrue = a.posTrue && a.eleTrue && b.posTrue && b.eleTrue;
    if (allTrue && a.posKey == a.eleKey && b.posKey == b.eleKey && a.posKey != b.posKey) {
      return 1;
    }
    if (allTrue && a.posKey == b.eleKey && a.eleKey == b.posKey && a.posKey != a.eleKey) {
      return 2;
    }
    return 3;
  }
  template <bool after>
  void fillMixPairs(const float type, std::vector<MixPairPhoton> const& photons)
  {
    for (size_t i = 0; i < photons.size(); ++i) {
      for (size_t j = i + 1; j < photons.size(); ++j) {
        const auto& a = photons[i];
        const auto& b = photons[j];
        if (after && skipSameSourcePair(type, a, b)) {
          continue;
        }
        const float q = qInv(a.px, a.py, a.pz, b.px, b.py, b.pz);
        const auto truth = static_cast<float>(mixPairTruth(a, b));
        if constexpr (after) {
          registry.fill(HIST("Mixed/Pairs/hQAfter"), type, 0.f, q);
          registry.fill(HIST("Mixed/Pairs/hQAfter"), type, truth, q);
        } else {
          registry.fill(HIST("Mixed/Pairs/hQBefore"), type, 0.f, q);
          registry.fill(HIST("Mixed/Pairs/hQBefore"), type, truth, q);
        }
      }
    }
  }

  // fill the occupancy of the source events and the deduplication diagnostics for the mixed-event photon pairs
  template <typename TMixCollisions, typename TMixLegs>
  void fillMatchedConfigs(TMixCollisions const& collisions, TMixLegs const& legs)
  {
    for (const auto& collision : collisions) { // same-event artificial events carry the occupancy of their source event
      if (!collision.isMixed()) {
        mOccOfSourceEvent[collision.sourceEventUid()] = collision.occupancy();
      }
    }
    using Key3 = std::tuple<int64_t, int64_t, int64_t>;
    std::set<Key3> cands, stored;
    for (const auto& c : vecV0Dedup) {
      cands.insert(Key3{c.colID, c.posID, c.eleID});
    }
    std::unordered_map<int64_t, int64_t> storedPartner;
    for (const auto& k : stored_fullv0Ids) {
      stored.insert(Key3{k.collisionId, k.posId, k.eleId});
      storedPartner.try_emplace(k.posId, k.eleId);
      storedPartner.try_emplace(k.eleId, k.posId);
    }
    for (const auto& collision : collisions) {
      const int64_t colId = collision.globalIndex();
      const float type = collision.isMixed() ? 1.f : 0.f;
      // MC photons with both legs present: (source event, mother) -> (pos row, ele row)
      std::map<std::pair<int64_t, int>, std::array<int64_t, 2>> photonLegs;
      for (int i = collision.firstLeg(); i < collision.firstLeg() + collision.nLegs(); ++i) {
        const auto& l = mLegIt[i];
        if (l.mcMotherId() < 0 || l.mcMotherPdg() != kGamma) {
          continue;
        }
        auto [it, inserted] = photonLegs.try_emplace(std::make_pair(static_cast<int64_t>(l.sourceEventUid()), l.mcMotherId()), std::array<int64_t, 2>{-1, -1});
        it->second[l.sign() > 0 ? 0 : 1] = i;
      }
      std::vector<std::array<int64_t, 2>> photons;
      for (const auto& [k, rows] : photonLegs) {
        if (rows[0] >= 0 && rows[1] >= 0) {
          photons.push_back(rows);
        }
      }
      for (size_t i = 0; i < photons.size(); ++i) {
        for (size_t j = i + 1; j < photons.size(); ++j) {
          const auto& ap = mLegIt[photons[i][0]];
          const auto& am = mLegIt[photons[i][1]];
          const auto& bp = mLegIt[photons[j][0]];
          const auto& bm = mLegIt[photons[j][1]];
          // related: same mother (siblings) or one photon radiated by a leg particle of the other (bremsstrahlung), as mixRelated
          const bool sameEvt = ap.sourceEventUid() == bp.sourceEventUid();
          if (type > 0.5f && sameEvt && mixed.skipSameSourceCrossed) {
            continue; // injection: a configuration of two current-event photons is a same-event one
          }
          const float qTrue = qInv(ap.mcMotherPx(), ap.mcMotherPy(), ap.mcMotherPz(), bp.mcMotherPx(), bp.mcMotherPy(), bp.mcMotherPz());
          if (qTrue >= 0.05f) { // o2-linter: disable=magic-number (diagnostic range)
            continue;
          }
          const int gA = ap.mcGrandMotherId(), gB = bp.mcGrandMotherId();
          const bool related = sameEvt && ((gA >= 0 && gA == gB) || (gA >= 0 && (gA == bp.mcParticleId() || gA == bm.mcParticleId())) ||
                                           (gB >= 0 && (gB == ap.mcParticleId() || gB == am.mcParticleId())));
          const int64_t A = photons[i][0], a = photons[i][1], B = photons[j][0], b = photons[j][1]; // pos / ele rows
          float input = 2.f;
          if (ap.sourcePhotonUid() == am.sourcePhotonUid() && bp.sourcePhotonUid() == bm.sourcePhotonUid() && ap.sourcePhotonUid() != bp.sourcePhotonUid()) {
            input = 0.f;
          } else if (ap.sourcePhotonUid() == bm.sourcePhotonUid() && bp.sourcePhotonUid() == am.sourcePhotonUid() && ap.sourcePhotonUid() != bp.sourcePhotonUid()) {
            input = 1.f;
          }
          const int nGen = static_cast<int>(cands.contains(Key3{colId, A, a})) + static_cast<int>(cands.contains(Key3{colId, B, b}));
          const int nCross = static_cast<int>(cands.contains(Key3{colId, A, b})) + static_cast<int>(cands.contains(Key3{colId, B, a}));
          const bool sGen = stored.contains(Key3{colId, A, a}) && stored.contains(Key3{colId, B, b});
          const bool sCross = stored.contains(Key3{colId, A, b}) && stored.contains(Key3{colId, B, a});
          float dOcc = 0.f; // ME: occupancy of the photon's source event - current event, the photon with the larger |difference|
          if (type > 0.5f) {
            float best = 0.f;
            for (const int64_t uid : {static_cast<int64_t>(ap.sourceEventUid()), static_cast<int64_t>(bp.sourceEventUid())}) {
              const auto it = mOccOfSourceEvent.find(uid);
              const float d = it == mOccOfSourceEvent.end() ? -19999.f : it->second - collision.occupancy();
              if (std::fabs(d) > std::fabs(best)) {
                best = d;
              }
            }
            dOcc = best;
          }
          const int nTpcOnly = static_cast<int>(!ap.hasITS()) + static_cast<int>(!am.hasITS()) + static_cast<int>(!bp.hasITS()) + static_cast<int>(!bm.hasITS());
          registry.fill(HIST("Mixed/Decision/hMatched"), type, qTrue, input, static_cast<float>(nGen), static_cast<float>(nCross), sGen ? 0.f : (sCross ? 1.f : 2.f), related ? 1.f : 0.f, static_cast<float>(collision.numContrib()), static_cast<float>(collision.nLegs()) / 2.f, collision.occupancy(), std::clamp(dOcc, -19999.f, 19999.f), static_cast<float>(nTpcOnly));
          // central electrode: true conversion point (ME partners shifted by the alignment, < 2 cm) and the leg's tan lambda (unchanged by it)
          constexpr float kTpcRin = 85.f, kTpcRout = 247.f; // o2-linter: disable=magic-number (TPC drift volume, diagnostic)
          auto crossesCE = [&](auto const& l) {
            const float rConv = std::hypot(l.mcVx(), l.mcVy());
            if (rConv > kTpcRout) {
              return 0;
            }
            const float z1 = l.mcVz() + l.tgl() * (std::max(rConv, kTpcRin) - rConv), z2 = l.mcVz() + l.tgl() * (kTpcRout - rConv);
            return z1 * z2 <= 0.f ? 1 : 0;
          };
          auto absEta = [](float px, float py, float pz) { return std::fabs(std::asinh(pz / std::max(std::hypot(px, py), 1.e-6f))); };
          const float minEta = std::min(absEta(ap.mcMotherPx(), ap.mcMotherPy(), ap.mcMotherPz()), absEta(bp.mcMotherPx(), bp.mcMotherPy(), bp.mcMotherPz()));
          const float minZ = std::min(std::fabs(ap.mcVz()), std::fabs(bp.mcVz()));
          const int nCE = crossesCE(ap) + crossesCE(am) + crossesCE(bp) + crossesCE(bm);
          registry.fill(HIST("Mixed/Decision/hMatchedCE"), type, qTrue, input, sGen ? 0.f : (sCross ? 1.f : 2.f), related ? 1.f : 0.f, std::min(minEta, 1.99f), std::min(minZ, 299.f), static_cast<float>(nCE), static_cast<float>(nTpcOnly));
          // leg fates: 0 not stored, 1 own genuine partner, 2 crossed partner, 3 leg of a third complete photon, 4 leg of an incomplete MC photon, 5 no photon mother
          std::array<int, 6> nFate{};
          const std::array<int64_t, 4> rows{A, a, B, b};
          const std::array<int64_t, 4> own{a, A, b, B}, crossed{b, B, a, A};
          for (int k = 0; k < 4; ++k) {
            const auto it = storedPartner.find(rows[k]);
            if (it == storedPartner.end()) {
              ++nFate[0];
            } else if (it->second == own[k]) {
              ++nFate[1];
            } else if (it->second == crossed[k]) {
              ++nFate[2];
            } else {
              const auto p = legs.rawIteratorAt(it->second);
              if (p.mcMotherId() < 0 || p.mcMotherPdg() != kGamma) {
                ++nFate[5];
              } else {
                const auto pIt = photonLegs.find(std::make_pair(static_cast<int64_t>(p.sourceEventUid()), p.mcMotherId()));
                const bool complete = pIt != photonLegs.end() && pIt->second[0] >= 0 && pIt->second[1] >= 0;
                ++nFate[complete ? 3 : 4];
              }
            }
          }
          registry.fill(HIST("Mixed/Decision/hMatchedLegs"), type, qTrue, input, sGen ? 0.f : (sCross ? 1.f : 2.f), related ? 1.f : 0.f,
                        static_cast<float>(nFate[0]), static_cast<float>(nFate[1]), static_cast<float>(nFate[2]), static_cast<float>(nFate[3]),
                        static_cast<float>(nFate[4]), static_cast<float>(nFate[5]));
        }
      }
    }
  }

  // truth of a photon pair: 0 both genuine, 1 crossed, 2 other (at least one leg not from a true photon), for all candidates of a collision
  template <typename TMixCollisions, typename TMixLegs>
  void fillFourLegAll(TMixCollisions const& collisions, TMixLegs const& legs, const size_t nFit)
  {
    using Key3 = std::tuple<int64_t, int64_t, int64_t>;
    std::map<Key3, int> firstCand;
    std::map<int64_t, std::vector<int>> originalsOfColl;
    for (size_t i = 0; i < vecV0Dedup.size(); ++i) {
      const auto& c = vecV0Dedup[i];
      if (firstCand.emplace(Key3{c.colID, c.posID, c.eleID}, static_cast<int>(i)).second && c.v0ID >= 0 && i < nFit) {
        originalsOfColl[c.colID].push_back(static_cast<int>(i));
      }
    }
    std::set<Key3> storedSet;
    for (const auto& k : stored_fullv0Ids) {
      storedSet.insert(Key3{k.collisionId, k.posId, k.eleId});
    }
    auto truePhoton = [&](int64_t posRow, int64_t eleRow) {
      const auto p = legs.rawIteratorAt(posRow);
      const auto e = legs.rawIteratorAt(eleRow);
      return p.mcMotherId() >= 0 && p.sourceEventUid() == e.sourceEventUid() && p.mcMotherId() == e.mcMotherId() && p.mcMotherPdg() == kGamma;
    };
    for (const auto& [colId, orig] : originalsOfColl) {
      const float type = collisions.rawIteratorAt(colId).isMixed() ? 1.f : 0.f;
      for (size_t i = 0; i < orig.size(); ++i) {
        for (size_t j = i + 1; j < orig.size(); ++j) {
          const auto& a = vecV0Dedup[orig[i]];
          const auto& b = vecV0Dedup[orig[j]];
          const auto itX1 = firstCand.find(Key3{colId, a.posID, b.eleID});
          const auto itX2 = firstCand.find(Key3{colId, b.posID, a.eleID});
          if (itX1 == firstCand.end() || itX2 == firstCand.end()) {
            continue;
          }
          const auto& x1 = vecV0Dedup[itX1->second];
          const auto& x2 = vecV0Dedup[itX2->second];
          const auto ap = legs.rawIteratorAt(a.posID);
          const auto bp = legs.rawIteratorAt(b.posID);
          const bool distinctMothers = ap.sourceEventUid() != bp.sourceEventUid() || ap.mcMotherId() != bp.mcMotherId();
          const bool physO = truePhoton(a.posID, a.eleID) && truePhoton(b.posID, b.eleID) && distinctMothers;
          const bool physX = truePhoton(a.posID, b.eleID) && truePhoton(b.posID, a.eleID) && distinctMothers;
          const float truth = physO ? 0.f : (physX ? 1.f : 2.f);
          const float qGen = physX ? qInv(x1.pxAtPV, x1.pyAtPV, x1.pzAtPV, x2.pxAtPV, x2.pyAtPV, x2.pzAtPV) : qInv(a.pxAtPV, a.pyAtPV, a.pzAtPV, b.pxAtPV, b.pyAtPV, b.pzAtPV);
          if (qGen >= 0.05f) { // o2-linter: disable=magic-number (diagnostic range)
            continue;
          }
          const bool storedO = storedSet.contains(Key3{colId, a.posID, a.eleID}) && storedSet.contains(Key3{colId, b.posID, b.eleID});
          const bool storedX = storedSet.contains(Key3{colId, a.posID, b.eleID}) && storedSet.contains(Key3{colId, b.posID, a.eleID});
          const float actual = storedO ? 0.f : (storedX ? 1.f : 2.f);
          const int crossedOrigin = (static_cast<size_t>(itX1->second) >= nFit ? 1 : 0) + (static_cast<size_t>(itX2->second) >= nFit ? 1 : 0);
          const int origOrigin = (ap.sourceV0Id() < 0 ? 1 : 0) + (bp.sourceV0Id() < 0 ? 1 : 0);
          registry.fill(HIST("Mixed/Decision/hFourLegAll"), type, truth, actual, qGen, static_cast<float>(crossedOrigin), static_cast<float>(origOrigin));
        }
      }
    }
  }

  template <typename TMixCollisions, typename TMixLegs>
  void fillFourLegTruthAligned(TMixCollisions const& collisions, TMixLegs const& legs, const int64_t colId, const float type, const float truth,
                               const float isolated, const float actual, V0CandidateHelper const& a, V0CandidateHelper const& b,
                               V0CandidateHelper const& x1, V0CandidateHelper const& x2, const int groupSize)
  {
    const bool origTrue = truth < 0.5f;
    const auto& t1 = origTrue ? a : x1;
    const auto& t2 = origTrue ? b : x2;
    const auto& w1 = origTrue ? x1 : a;
    const auto& w2 = origTrue ? x2 : b;
    // decision in truth terms: 0 genuine assignment, 1 wrong assignment, 2 other
    auto inTruth = [&](float decision) { return decision > 1.5f ? 2.f : ((decision < 0.5f) == origTrue ? 0.f : 1.f); };
    const float choice = inTruth(isolated);
    MixPairPhoton p1, p2;
    setMixPairLegs(p1, legs.rawIteratorAt(t1.posID), legs.rawIteratorAt(t1.eleID));
    setMixPairLegs(p2, legs.rawIteratorAt(t2.posID), legs.rawIteratorAt(t2.eleID));
    const int relKind = mixRelKind(p1, p2);
    const auto rel = static_cast<float>(relKind);
    const float rel2 = relKind > 0 ? 1.f : 0.f;
    const float s0 = dedup.groupMatchingScoreThreshold.value;
    const float qTrue = qInv(t1.pxAtPV, t1.pyAtPV, t1.pzAtPV, t2.pxAtPV, t2.pyAtPV, t2.pzAtPV);
    const int nClipT = (t1.score >= s0 ? 1 : 0) + (t2.score >= s0 ? 1 : 0);
    const int nClipW = (w1.score >= s0 ? 1 : 0) + (w2.score >= s0 ? 1 : 0);
    auto angle = [](float cpa) { return std::acos(std::clamp(cpa, -1.f, 1.f)); };
    const float dTheta = 1e3f * ((angle(w1.cosPA) + angle(w2.cosPA)) - (angle(t1.cosPA) + angle(t2.cosPA))); // mrad
    const float dPCA = (w1.pca + w2.pca) - (t1.pca + t2.pca);
    const std::array<int64_t, 4> rows{t1.posID, t1.eleID, t2.posID, t2.eleID};
    float fshMax = 0.f;
    int nITS = 0;
    for (const auto& r : rows) {
      const auto l = legs.rawIteratorAt(r);
      fshMax = std::max(fshMax, static_cast<float>(l.tpcFractionSharedCls()));
      nITS += l.hasITS() ? 1 : 0;
    }
    registry.fill(HIST("Mixed/Decision/hFourLegTruth"), type, truth, choice, rel, qTrue, static_cast<float>(3 * nClipT + nClipW), dTheta, dPCA, fshMax);
    {
      const int origin = mixCommonOrigin(p1.posMC, p2.posMC);
      auto rootClass = [](int pdg) {
        switch (std::abs(pdg)) {
          case 111:
            return 1.f; // o2-linter: disable=pdg/explicit-code (diagnostic classes)
          case 221:
            return 2.f; // o2-linter: disable=pdg/explicit-code
          case 223:
            return 3.f; // o2-linter: disable=pdg/explicit-code
          case 331:
            return 4.f; // o2-linter: disable=pdg/explicit-code
          case 113:
          case 213:
          case 313:
          case 323:
            return 5.f; // o2-linter: disable=pdg/explicit-code
          case 310:
          case 130:
          case 311:
            return 6.f; // o2-linter: disable=pdg/explicit-code
          default:
            return std::abs(pdg) > 1000 ? 7.f : 0.f; // o2-linter: disable=magic-number (baryons)
        }
      };
      registry.fill(HIST("Mixed/Decision/hFourLegOrigin"), type, truth, choice, rel, qTrue, static_cast<float>(origin), origin > 0 ? rootClass(p1.posMC.rootPdg) : 0.f);
    }

    // TPC leg proximity: the minimum distance of the two legs of one photon to the two legs of the other photon, at 17 radii in the TPC
    constexpr int NRadii = 17;
    auto positions = [&](int64_t row) {
      std::array<std::array<float, 3>, NRadii> xyz{};
      std::array<bool, NRadii> ok{};
      const auto tp0 = getTrackPar(legs.rawIteratorAt(row));
      for (int k = 0; k < NRadii; ++k) {
        auto tp = tp0;
        float x = 0.f;
        ok[k] = tp.getXatLabR(85.f + 10.f * k, x, d_bz) && tp.propagateParamTo(x, d_bz); // o2-linter: disable=magic-number (TPC radii)
        if (ok[k]) {
          tp.getXYZGlo(xyz[k]);
        }
      }
      return std::make_pair(xyz, ok);
    };
    const auto pT1p = positions(t1.posID), pT1e = positions(t1.eleID), pT2p = positions(t2.posID), pT2e = positions(t2.eleID);
    auto minDist = [&](auto const& u, auto const& v) {
      float d = 999.f;
      for (int k = 0; k < NRadii; ++k) {
        if (u.second[k] && v.second[k]) {
          d = std::min(d, std::hypot(u.first[k][0] - v.first[k][0], u.first[k][1] - v.first[k][1], u.first[k][2] - v.first[k][2]));
        }
      }
      return d;
    };
    const float dOS = std::min(minDist(pT1p, pT2e), minDist(pT2p, pT1e)); // legs that the wrong assignment pairs
    const float dSS = std::min(minDist(pT1p, pT2p), minDist(pT1e, pT2e)); // control: same-sign legs
    const float occ = collisions.rawIteratorAt(colId).occupancy();
    registry.fill(HIST("Mixed/Decision/hFourLegTrack"), type, truth, choice, rel2, std::log10(std::max(dOS, 1e-3f)), std::log10(std::max(dSS, 1e-3f)),
                  static_cast<float>(nITS), occ); // o2-linter: disable=magic-number (log floor)

    auto mother = [&](int64_t row) {
      const auto l = legs.rawIteratorAt(row);
      return std::array<float, 3>{l.mcMotherPx(), l.mcMotherPy(), l.mcMotherPz()};
    };
    const auto m1 = mother(t1.posID), m2 = mother(t2.posID);
    auto phi = [](float px, float py) { return std::atan2(py, px); };
    auto eta = [](float px, float py, float pz) { return std::asinh(pz / std::max(std::hypot(px, py), 1e-6f)); };
    const float uPhi = RecoDecay::constrainAngle(phi(m2[0], m2[1]) - phi(m1[0], m1[1]), -o2::constants::math::PI) > 0.f ? 1.f : -1.f;
    const float uEta = eta(m2[0], m2[1], m2[2]) - eta(m1[0], m1[1], m1[2]) > 0.f ? 1.f : -1.f;
    const float rPhi1 = uPhi * RecoDecay::constrainAngle(phi(t1.pxAtPV, t1.pyAtPV) - phi(m1[0], m1[1]), -o2::constants::math::PI);
    const float rPhi2 = -uPhi * RecoDecay::constrainAngle(phi(t2.pxAtPV, t2.pyAtPV) - phi(m2[0], m2[1]), -o2::constants::math::PI);
    const float rEta1 = uEta * (eta(t1.pxAtPV, t1.pyAtPV, t1.pzAtPV) - eta(m1[0], m1[1], m1[2]));
    const float rEta2 = -uEta * (eta(t2.pxAtPV, t2.pyAtPV, t2.pzAtPV) - eta(m2[0], m2[1], m2[2]));
    registry.fill(HIST("Mixed/Decision/hFourLegResid"), type, truth, choice, rel2, rPhi1, rPhi2, rEta1, rEta2);

    const std::vector<uint64_t> legMask{0b0011, 0b1100, 0b1001, 0b0110};
    auto decide = [&](std::vector<float> const& sc, float thr) {
      const auto sol = solveConflictGroup(legMask, sc, thr, false);
      return inTruth(sol.bestMask == 0b0011 ? 0.f : (sol.bestMask == 0b1100 ? 1.f : 2.f)); // o2-linter: disable=magic-number (subset masks)
    };
    auto scores = [&](float w) {
      return std::vector<float>{getScoreV0(a.cosPA, a.pca, w), getScoreV0(b.cosPA, b.pca, w), getScoreV0(x1.cosPA, x1.pca, w), getScoreV0(x2.cosPA, x2.pca, w)};
    };
    const std::vector<float> baseScores{a.score, b.score, x1.score, x2.score};
    registry.fill(HIST("Mixed/Decision/hFourLegAlt"), type, truth, rel2, choice, decide(scores(1.f), s0), decide(scores(0.f), s0), decide(baseScores, -1.f),
                  decide(scores(dedup.deduplicationScoreWeight.value), s0));

    const auto col = collisions.rawIteratorAt(colId);
    auto pointing = [&](V0CandidateHelper const& c) {
      const float dx = c.svX - col.posX(), dy = c.svY - col.posY(), dz = c.svZ - col.posZ();
      const float dPhi = RecoDecay::constrainAngle(std::atan2(c.pyAtSV, c.pxAtSV) - std::atan2(dy, dx), -o2::constants::math::PI);
      const float dLam = std::atan2(c.pzAtSV, std::hypot(c.pxAtSV, c.pyAtSV)) - std::atan2(dz, std::hypot(dx, dy));
      return std::make_pair(1e3f * dPhi, 1e3f * dLam); // mrad
    };
    const auto pt1 = pointing(t1), pt2 = pointing(t2), pw1 = pointing(w1), pw2 = pointing(w2);
    const float dPhiSum = (std::fabs(pw1.first) + std::fabs(pw2.first)) - (std::fabs(pt1.first) + std::fabs(pt2.first));
    const float dLamSum = (std::fabs(pw1.second) + std::fabs(pw2.second)) - (std::fabs(pt1.second) + std::fabs(pt2.second));
    registry.fill(HIST("Mixed/Decision/hFourLegPoint"), type, truth, choice, rel2, qTrue, dPhiSum, dLamSum, pt1.second, pt2.second, pt1.first, pt2.first);

    // log dump of the few cases where the isolated decision is wrong, but the actual decision is correct and genuine
    if (type < 0.5f && origTrue && isolated > 0.5f && isolated < 1.5f && actual < 0.5f && mNAnomalyDump < 200) { // o2-linter: disable=magic-number (cap)
      ++mNAnomalyDump;
      LOGF(info, "[4leg-anom] col %ld group %d scores A %.4f B %.4f X1 %.4f X2 %.4f s0 %.2f q %.4f rel %d", colId, groupSize, a.score, b.score, x1.score, x2.score, s0, qTrue, relKind);
    }
  }

  template <typename TMixCollisions, typename TMixLegs>
  void fillFourLegDecisions(TMixCollisions const& collisions, TMixLegs const& legs, const size_t nFit)
  {
    fillFourLegAll(collisions, legs, nFit);
    using Key3 = std::tuple<int64_t, int64_t, int64_t>;
    std::map<Key3, int> firstCand; // first fit of a leg pair = its dedup candidate
    std::map<int64_t, std::vector<int>> originalsOfColl;
    for (size_t i = 0; i < nFit; ++i) {
      const auto& c = vecV0Dedup[i];
      if (firstCand.emplace(Key3{c.colID, c.posID, c.eleID}, static_cast<int>(i)).second && c.v0ID >= 0) {
        originalsOfColl[c.colID].push_back(static_cast<int>(i));
      }
    }
    std::set<Key3> storedSet;
    for (const auto& k : stored_fullv0Ids) {
      storedSet.insert(Key3{k.collisionId, k.posId, k.eleId});
    }
    auto truePhoton = [&](int64_t posRow, int64_t eleRow) {
      const auto p = legs.rawIteratorAt(posRow);
      const auto e = legs.rawIteratorAt(eleRow);
      return p.mcMotherId() >= 0 && p.sourceEventUid() == e.sourceEventUid() && p.mcMotherId() == e.mcMotherId() && p.mcMotherPdg() == kGamma;
    };
    const float s0 = dedup.groupMatchingScoreThreshold.value;
    for (const auto& [colId, orig] : originalsOfColl) {
      const float type = collisions.rawIteratorAt(colId).isMixed() ? 1.f : 0.f;
      for (size_t i = 0; i < orig.size(); ++i) {
        for (size_t j = i + 1; j < orig.size(); ++j) {
          const auto& a = vecV0Dedup[orig[i]];
          const auto& b = vecV0Dedup[orig[j]];
          const auto itX1 = firstCand.find(Key3{colId, a.posID, b.eleID});
          const auto itX2 = firstCand.find(Key3{colId, b.posID, a.eleID});
          if (itX1 == firstCand.end() || itX2 == firstCand.end()) {
            continue; // not both assignments available
          }
          const auto& x1 = vecV0Dedup[itX1->second];
          const auto& x2 = vecV0Dedup[itX2->second];
          const float qO = qInv(a.pxAtPV, a.pyAtPV, a.pzAtPV, b.pxAtPV, b.pyAtPV, b.pzAtPV);
          const float qX = qInv(x1.pxAtPV, x1.pyAtPV, x1.pzAtPV, x2.pxAtPV, x2.pyAtPV, x2.pzAtPV);
          if (std::min(qO, qX) >= 0.1f) { // o2-linter: disable=magic-number (diagnostic range)
            continue;
          }
          // the isolated decision is the dedup objective on only these four candidates, the actual decision is the stored set after the dedup and the resolvers
          const std::vector<uint64_t> legMask{0b0011, 0b1100, 0b1001, 0b0110};
          const std::vector<float> scores{a.score, b.score, x1.score, x2.score};
          const auto sol = solveConflictGroup(legMask, scores, s0, false);
          const float isolated = sol.bestMask == 0b0011 ? 0.f : (sol.bestMask == 0b1100 ? 1.f : 2.f); // o2-linter: disable=magic-number (subset masks)
          const bool storedO = storedSet.contains(Key3{colId, a.posID, a.eleID}) && storedSet.contains(Key3{colId, b.posID, b.eleID});
          const bool storedX = storedSet.contains(Key3{colId, a.posID, b.eleID}) && storedSet.contains(Key3{colId, b.posID, a.eleID});
          const float actual = storedO ? 0.f : (storedX ? 1.f : 2.f);
          const auto ap = legs.rawIteratorAt(a.posID);
          const auto bp = legs.rawIteratorAt(b.posID);
          const bool distinctMothers = ap.sourceEventUid() != bp.sourceEventUid() || ap.mcMotherId() != bp.mcMotherId();
          const bool physO = truePhoton(a.posID, a.eleID) && truePhoton(b.posID, b.eleID) && distinctMothers;
          const bool physX = truePhoton(a.posID, b.eleID) && truePhoton(b.posID, a.eleID) && distinctMothers;
          const float truth = physO ? 0.f : (physX ? 1.f : 2.f);
          const int groupSize = std::max({mGroupSizeOfCand[orig[i]], mGroupSizeOfCand[orig[j]], mGroupSizeOfCand[itX1->second], mGroupSizeOfCand[itX2->second]});
          // originals without an SVertexer V0 are fitted without the AP cut, the crossed pairings with it (asymmetric options)
          const bool relaxedOriginal = ap.sourceV0Id() < 0 || bp.sourceV0Id() < 0;
          registry.fill(HIST("Mixed/Decision/hFourLeg"), type, truth, isolated, actual, qO, qX, static_cast<float>(std::min(groupSize, 15)), // o2-linter: disable=magic-number (last bin)
                        (x1.score + x2.score) - (a.score + b.score), relaxedOriginal ? 1.f : 0.f);
          if (std::min(qO, qX) < 0.05f) { // o2-linter: disable=magic-number (diagnostic range)
            // what tips the 4-leg decision: score components, geometry and kinematics of the two conversions, MC relation
            auto angle = [](float cpa) { return std::acos(std::clamp(cpa, -1.f, 1.f)); };
            const float dPCA = (x1.pca + x2.pca) - (a.pca + b.pca);
            const float dAngle = (angle(x1.cosPA) + angle(x2.cosPA)) - (angle(a.cosPA) + angle(b.cosPA));
            const float dR = std::fabs(std::hypot(a.svX, a.svY) - std::hypot(b.svX, b.svY));
            const float dZ = std::fabs(a.svZ - b.svZ);
            auto pLeg = [&](int64_t row) {
              const auto l = legs.rawIteratorAt(row);
              return std::fabs(1.f / l.signed1Pt()) * std::sqrt(1.f + l.tgl() * l.tgl());
            };
            auto asym = [](float p1, float p2) { return std::fabs(p1 - p2) / std::max(p1 + p2, 1e-6f); };
            const float sameSignSim = std::min(asym(pLeg(a.posID), pLeg(b.posID)), asym(pLeg(a.eleID), pLeg(b.eleID)));
            MixPairPhoton pa, pb;
            setMixPairLegs(pa, legs.rawIteratorAt(a.posID), legs.rawIteratorAt(a.eleID));
            setMixPairLegs(pb, legs.rawIteratorAt(b.posID), legs.rawIteratorAt(b.eleID));
            registry.fill(HIST("Mixed/Decision/hFourLegKin"), type, truth, isolated, dPCA, dAngle, dR, dZ, sameSignSim, mixRelated(pa, pb) ? 1.f : 0.f);
          }
          if (truth < 1.5f && std::min(qO, qX) < 0.05f) { // o2-linter: disable=magic-number (genuine assignment exists; diagnostic range)
            fillFourLegTruthAligned(collisions, legs, colId, type, truth, isolated, actual, a, b, x1, x2, groupSize);
          }
        }
      }
    }
    for (const auto& r : mPermuteResults) {
      registry.fill(HIST("Mixed/Decision/hPermute"), collisions.rawIteratorAt(r.colID).isMixed() ? 1.f : 0.f, static_cast<float>(std::min(r.nMembers, 30)), // o2-linter: disable=magic-number (last bin)
                    r.changed ? 1.f : 0.f, r.greedy ? 1.f : 0.f);
    }
  }

  // run 2: build the mixed event, fill the crossed pairings, and fill the diagnostics of the run-1 photons (pooled free legs) before the dedup and resolvers
  template <typename TMixCollisions, typename TMixLegs>
  void buildMixed(TMixCollisions const& collisions, TMixLegs const& legs)
  {
    const auto tDF = std::chrono::steady_clock::now();
    mRelaxedKeys.clear();
    mLegDcaCache.clear();
    // build the leg iterator cache for the current event
    mLegIt.clear();
    mLegIt.reserve(legs.size());
    for (int64_t i = 0; i < static_cast<int64_t>(legs.size()); ++i) {
      mLegIt.push_back(legs.rawIteratorAt(i));
    }
    std::unordered_map<int64_t, std::vector<MixPairPhoton>> inputOfColl; // artificial event -> input photons (pair diagnostics)
    auto crossedPrefilter = [&](auto const& pos, auto const& ele) {
      if (std::fabs(pos.tgl() - ele.tgl()) > recovery.finderMaxDeltaTgl) {
        return false;
      }
      const o2::track::TrackAuxPar cp(getTrackPar(pos), d_bz);
      const o2::track::TrackAuxPar ce(getTrackPar(ele), d_bz);
      return std::fabs(circleGap(cp, ce)) < recovery.finderMaxCircleGap;
    };
    for (const auto& collision : collisions) {
      const float type = collision.isMixed() ? 1.f : 0.f;
      const auto t0 = std::chrono::steady_clock::now();
      registry.fill(HIST("Mixed/hEvents"), type);
      std::vector<int64_t> posLegs, eleLegs;
      for (int i = collision.firstLeg(); i < collision.firstLeg() + collision.nLegs(); ++i) {
        (mLegIt[i].sign() > 0 ? posLegs : eleLegs).push_back(i);
      }
      std::unordered_map<int64_t, o2::track::TrackAuxPar> circleOf;
      std::unordered_map<int64_t, int64_t> eleOfPhoton;
      std::unordered_map<int64_t, int> nEleOfSource;
      std::vector<std::pair<float, int64_t>> eleByTgl;
      eleByTgl.reserve(eleLegs.size());
      for (const auto& ie : eleLegs) {
        eleOfPhoton.emplace(mLegIt[ie].sourcePhotonUid(), ie);
        nEleOfSource[mLegIt[ie].sourceEventUid()]++;
        eleByTgl.emplace_back(mLegIt[ie].tgl(), ie);
      }
      std::sort(eleByTgl.begin(), eleByTgl.end());
      const bool tglWindow = mixed.prefilterCrossed && !mixed.uniformPairFit;
      if (mixed.prefilterCrossed) {
        for (const auto& i : posLegs) {
          circleOf.emplace(i, o2::track::TrackAuxPar(getTrackPar(mLegIt[i]), d_bz));
        }
        for (const auto& i : eleLegs) {
          circleOf.emplace(i, o2::track::TrackAuxPar(getTrackPar(mLegIt[i]), d_bz));
        }
      }
      auto crossedPrefilterCached = [&](int64_t ip, int64_t ie) {
        return std::fabs(mLegIt[ip].tgl() - mLegIt[ie].tgl()) <= recovery.finderMaxDeltaTgl && std::fabs(circleGap(circleOf.at(ip), circleOf.at(ie))) < recovery.finderMaxCircleGap;
      };
      std::vector<MixPairPhoton> before; // crossed pairings of this artificial event, before the dedup and resolvers
      for (const auto& ip : posLegs) {
        const auto& pos = mLegIt[ip];
        if (const auto itE = eleOfPhoton.find(pos.sourcePhotonUid()); itE != eleOfPhoton.end()) {
          const auto ie = itE->second;
          const auto& ele = mLegIt[ie];
          {
            MixPairPhoton p;
            p.px = pos.sourcePx();
            p.py = pos.sourcePy();
            p.pz = pos.sourcePz();
            setMixPairLegs(p, pos, ele);
            p.posRow = ip;
            p.eleRow = ie;
            p.vx = pos.sourceVx();
            p.vy = pos.sourceVy();
            p.vz = pos.sourceVz();
            before.push_back(p);
          }
        }
      }
      registry.fill(HIST("Mixed/hNPhotonsIn"), type, static_cast<float>(before.size())); // run-1 photons (pooled free legs are single legs)
      fillMixPairs<false>(type, before);
      inputOfColl[collision.globalIndex()] = before;
      std::vector<int64_t> visit; // negative legs paired with this positive leg, in the original (row) order
      for (const auto& ip : posLegs) {
        const auto& pos = mLegIt[ip];
        visit.clear();
        if (tglWindow) { // crossed pairings are prefiltered by tgl, then by circle gap
          const float tglP = pos.tgl();
          for (auto it = std::lower_bound(eleByTgl.begin(), eleByTgl.end(), std::make_pair(tglP - recovery.finderMaxDeltaTgl - 1e-4f, int64_t{-1})); it != eleByTgl.end() && it->first <= tglP + recovery.finderMaxDeltaTgl + 1e-4f; ++it) {
            visit.push_back(it->second);
          }
          if (const auto itE = eleOfPhoton.find(pos.sourcePhotonUid()); itE != eleOfPhoton.end()) {
            visit.push_back(itE->second);
          }
          std::sort(visit.begin(), visit.end());
          visit.erase(std::unique(visit.begin(), visit.end()), visit.end());
          // crossed pairings outside the window: counted in the cut flow as prefilter rejections (not split by MC in hCrossedFit)
          const bool hasPartner = eleOfPhoton.contains(pos.sourcePhotonUid());
          int nOutside = mixed.skipSameSourceCrossed ? static_cast<int>(eleLegs.size()) - nEleOfSource[pos.sourceEventUid()] : static_cast<int>(eleLegs.size()) - (hasPartner ? 1 : 0);
          for (const auto& ie : visit) {
            const auto& ele = mLegIt[ie];
            if (ele.sourcePhotonUid() != pos.sourcePhotonUid() && !(mixed.skipSameSourceCrossed && pos.sourceEventUid() == ele.sourceEventUid())) {
              nOutside--;
            }
          }
          if (nOutside > 0) {
            registry.fill(HIST("Mixed/hCandidates"), type, 2.f, static_cast<double>(nOutside));
          }
        } else {
          visit.assign(eleLegs.begin(), eleLegs.end());
        }
        for (const auto& ie : visit) {
          const auto& ele = mLegIt[ie];
          const size_t nBefore = vecV0Dedup.size();
          const bool originalPairing = pos.sourcePhotonUid() == ele.sourcePhotonUid();
          if (!originalPairing && mixed.skipSameSourceCrossed && pos.sourceEventUid() == ele.sourceEventUid()) {
            continue; // injection: same-source crossed pairings were decided by run 1
          }
          if (!originalPairing && mixed.skipFreeFreeCrossed && pos.sourceV0Id() == -2 && ele.sourceV0Id() == -2) { // o2-linter: disable=magic-number (pooled free leg)
            continue;                                                                                              // overlay: two free legs, rejected by run 1
          }
          if (!originalPairing && mixed.freeFreeFinderOnly && pos.sourceV0Id() == -2 && ele.sourceV0Id() == -2 && pos.hasITS() && ele.hasITS()) { // o2-linter: disable=magic-number (pooled free leg)
            continue;                                                                                                                             // overlay: two free legs with ITS, only an SVertexer V0 in run 1
          }
          if (mixed.uniformPairFit) {
            // the crossed pairings are fitted with the same options as the originals, but without the AP cut if they were not found by an SVertexer V0
            const bool finderLike = !(pos.hasITS() && ele.hasITS());
            registry.fill(HIST("Mixed/hCandidates"), type, originalPairing ? 0.f : 2.f);
            if (mixed.prefilterCrossed && !crossedPrefilter(pos, ele)) {
              continue;
            }
            if (!originalPairing) {
              registry.fill(HIST("Mixed/hCandidates"), type, 3.f);
            }
            const SourcePV srcPv{pos.sourcePvX(), pos.sourcePvY(), pos.sourcePvZ(), pos.sourcePvCovXX(), pos.sourcePvCovXY(), pos.sourcePvCovXZ(),
                                 pos.sourcePvCovYY(), pos.sourcePvCovYZ(), pos.sourcePvCovZZ(), pos.sourcePvChi2(), pos.sourcePvNContrib()};
            fitV0Candidate<false, void, TMixCollisions, TMixLegs>(collision, pos, ele, originalPairing ? ip : -1, false,
                                                                  FitOptions{.nonSVertexer = finderLike, .relaxAP = finderLike, .purityGate = finderLike && recovery.runFinder.value, .pv = (originalPairing && mixed.refitOriginalsAtSourceVertex) ? &srcPv : nullptr});
            if (vecV0Dedup.size() > nBefore) {
              registry.fill(HIST("Mixed/hCandidates"), type, originalPairing ? 1.f : 4.f);
            }
            continue;
          }
          if (originalPairing) { // the original pairings are fitted with the same options as in run 1, but without the AP cut if they were not found by an SVertexer V0
            registry.fill(HIST("Mixed/hCandidates"), type, 0.f);
            const SourcePV srcPv{pos.sourcePvX(), pos.sourcePvY(), pos.sourcePvZ(), pos.sourcePvCovXX(), pos.sourcePvCovXY(), pos.sourcePvCovXZ(),
                                 pos.sourcePvCovYY(), pos.sourcePvCovYZ(), pos.sourcePvCovZZ(), pos.sourcePvChi2(), pos.sourcePvNContrib()};
            fitV0Candidate<false, void, TMixCollisions, TMixLegs>(collision, pos, ele, ip, false,
                                                                  FitOptions{.nonSVertexer = pos.sourceV0Id() < 0 && !recovery.finderUniformFit.value, .relaxAP = pos.sourceV0Id() < 0 && !recovery.finderUniformFit.value, .pv = mixed.refitOriginalsAtSourceVertex ? &srcPv : nullptr});
            if (vecV0Dedup.size() > nBefore) {
              registry.fill(HIST("Mixed/hCandidates"), type, 1.f);
            }
            continue;
          }
          registry.fill(HIST("Mixed/hCandidates"), type, 2.f);
          // crossed pairings: the same prefilter and the same leg-based fit options for every pairing (original pairings keep v0Id = row of the positive leg)
          const bool twoPhotons = pos.mcMotherPdg() == kGamma && ele.mcMotherPdg() == kGamma && pos.mcMotherId() >= 0 && ele.mcMotherId() >= 0 &&
                                  (pos.sourceEventUid() != ele.sourceEventUid() || pos.mcMotherId() != ele.mcMotherId());
          const float qTwo = twoPhotons ? qInv(pos.mcMotherPx(), pos.mcMotherPy(), pos.mcMotherPz(), ele.mcMotherPx(), ele.mcMotherPy(), ele.mcMotherPz()) : 1.f;
          const float genuineClose = (twoPhotons && qTwo < 0.05f) ? 1.f : 0.f; // o2-linter: disable=magic-number (diagnostic range)
          const auto nITSx = static_cast<float>((pos.hasITS() ? 1 : 0) + (ele.hasITS() ? 1 : 0));
          if (mixed.prefilterCrossed && !crossedPrefilterCached(ip, ie)) {
            if (mixed.mixedDiagnostics) {
              registry.fill(HIST("Mixed/Decision/hCrossedFit"), type, genuineClose, 1.f, nITSx, std::min(qTwo, 0.0499f));
            }
            continue;
          }
          registry.fill(HIST("Mixed/hCandidates"), type, 3.f);
          const bool anyFree = pos.sourceV0Id() == -2 || ele.sourceV0Id() == -2; // o2-linter: disable=magic-number (pooled free leg)
          const bool freeLegPair = mixed.freeLegRun1Rules && anyFree;
          const bool finderLikeFree = freeLegPair && !(pos.hasITS() && ele.hasITS());
          if (mixed.diagFreeLegFits && anyFree) { // trial fits with both option sets, candidates removed again
            const bool finderLike = !(pos.hasITS() && ele.hasITS());
            auto trial = [&](FitOptions const& o) {
              const size_t n0 = vecV0Dedup.size();
              fitV0Candidate<false, void, TMixCollisions, TMixLegs>(collision, pos, ele, -1, false, o);
              // cppcheck-suppress knownConditionTrueFalse
              if (vecV0Dedup.size() > n0) {
                const auto k = keyOf(vecV0Dedup.back());
                pca_map.erase(k);
                cospa_map.erase(k);
                vecV0Dedup.pop_back();
                return 0;
              }
              return std::max(mLastFitReject, 2);
            };
            const int lenient = trial(FitOptions{.relaxAP = mixed.crossedRelaxAP.value});
            const int run1 = trial(FitOptions{.nonSVertexer = finderLike, .purityGate = finderLike && recovery.runFinder.value});
            const auto nFree = static_cast<float>((pos.sourceV0Id() == -2 ? 1 : 0) + (ele.sourceV0Id() == -2 ? 1 : 0)); // o2-linter: disable=magic-number (pooled free leg)
            registry.fill(HIST("Mixed/Decision/hFreeLegFit"), type, genuineClose > 0.5f ? 0.f : 1.f, static_cast<float>(std::min(lenient, 1999)), static_cast<float>(std::min(run1, 1999)), nITSx, nFree);
          }
          const int mask = mixed.freeLegRuleMask.value;
          const FitOptions freeOpt{.nonSVertexer = finderLikeFree && (mask & 4) != 0, .relaxAP = (mask & 1) == 0 && mixed.crossedRelaxAP.value, // o2-linter: disable=magic-number (bit mask)
                                   .purityGate = finderLikeFree && (mask & 2) != 0 && recovery.runFinder.value,
                                   .minPurity = mixed.freeLegMinPurity.value}; // o2-linter: disable=magic-number (bit mask)
          fitV0Candidate<false, void, TMixCollisions, TMixLegs>(collision, pos, ele, -1, false, freeLegPair ? freeOpt : FitOptions{.relaxAP = mixed.crossedRelaxAP.value});
          if (vecV0Dedup.size() > nBefore) {
            registry.fill(HIST("Mixed/hCandidates"), type, 4.f);
          }
          if (mixed.mixedDiagnostics) {
            registry.fill(HIST("Mixed/Decision/hCrossedFit"), type, genuineClose, vecV0Dedup.size() > nBefore ? 0.f : static_cast<float>(mLastFitReject), nITSx, std::min(qTwo, 0.0499f));
          }
        }
      }
      registry.fill(HIST("Mixed/hTimeCandidates"), type, std::chrono::duration<float, std::milli>(std::chrono::steady_clock::now() - t0).count());
    }

    const size_t nFit = vecV0Dedup.size(); // candidates of the run-2 fit (the resolvers may append more)
    switch (dedup.deduplicationMode.value) {
      case V0DeduplicationMode::Pairwise:
        selectPairwise(nullptr);
        break;
      case V0DeduplicationMode::KeepAll:
        keepAll();
        break;
      case V0DeduplicationMode::GroupMatching:
        selectByConflictGroup(nullptr);
        break;
      default:
        selectGreedy(legs.size(), nullptr);
        break;
    }
    if (dedup.resolveSplitLegs) {
      resolveSplitLegConflicts(legs, nullptr);
    }
    if (recovery.runSwapResolver) {
      runSwapResolver<false, void, TMixCollisions, TMixLegs>(collisions, legs, nullptr);
    }
    if (mixed.mixedDiagnostics && dedup.deduplicationMode.value == static_cast<int>(V0DeduplicationMode::GroupMatching)) {
      fillFourLegDecisions(collisions, legs, nFit);
    }
    if (mixed.mixedDiagnostics) {
      fillMatchedConfigs(collisions, legs);
    }

    std::unordered_map<int64_t, std::array<int, 2>> nStored;       // artificial event -> (originals, crossed)
    std::unordered_map<int64_t, std::vector<MixPairPhoton>> after; // artificial event -> stored photons
    for (const auto& key : stored_fullv0Ids) {
      const auto collision = collisions.rawIteratorAt(key.collisionId);
      const auto pos = legs.rawIteratorAt(key.posId);
      const auto ele = legs.rawIteratorAt(key.eleId);
      const bool original = key.v0Id >= 0;
      const float type = collision.isMixed() ? 1.f : 0.f;
      const auto nWritten = v0photonskf.lastIndex();
      const SourcePV srcPv{pos.sourcePvX(), pos.sourcePvY(), pos.sourcePvZ(), pos.sourcePvCovXX(), pos.sourcePvCovXY(), pos.sourcePvCovXZ(),
                           pos.sourcePvCovYY(), pos.sourcePvCovYZ(), pos.sourcePvCovZZ(), pos.sourcePvChi2(), pos.sourcePvNContrib()};
      const bool finderLike = mixed.uniformPairFit && !(pos.hasITS() && ele.hasITS()); // uniformPairFit: options by the legs, as in the candidate fit
      fitV0Candidate<false, void, TMixCollisions, TMixLegs>(collision, pos, ele, key.v0Id, true,
                                                            FitOptions{.nonSVertexer = mixed.uniformPairFit ? finderLike : (original && pos.sourceV0Id() < 0 && !recovery.finderUniformFit.value),
                                                                       .relaxAP = mRelaxedKeys.contains(key) || (mixed.uniformPairFit ? finderLike : (original ? (pos.sourceV0Id() < 0 && !recovery.finderUniformFit.value) : mixed.crossedRelaxAP.value)),
                                                                       .purityGate = finderLike && recovery.runFinder.value,
                                                                       .pv = (original && mixed.refitOriginalsAtSourceVertex) ? &srcPv : nullptr});
      if (v0photonskf.lastIndex() == nWritten) {
        continue; // rejected by the ML cut (not a dedup candidate then either)
      }
      const auto& c = vecV0Dedup.back();
      nStored[key.collisionId][original ? 0 : 1]++;
      registry.fill(HIST("Mixed/hStored"), type, original ? 0.f : 1.f);
      if (original) {
        const float pt1 = std::hypot(pos.sourcePx(), pos.sourcePy());
        const float pt2 = std::hypot(c.pxAtPV, c.pyAtPV);
        registry.fill(HIST("Mixed/hOrigRelDPt"), type, pt2 / pt1 - 1.f);
        registry.fill(HIST("Mixed/hOrigDPhi"), type, RecoDecay::constrainAngle(RecoDecay::phi(c.pxAtPV, c.pyAtPV) - RecoDecay::phi(pos.sourcePx(), pos.sourcePy()), -o2::constants::math::PI));
        registry.fill(HIST("Mixed/hOrigDEta"), type, RecoDecay::eta(std::array{c.pxAtPV, c.pyAtPV, c.pzAtPV}) - RecoDecay::eta(std::array{pos.sourcePx(), pos.sourcePy(), pos.sourcePz()}));
        registry.fill(HIST("Mixed/hOrigExact"), type, (c.pxAtPV == pos.sourcePx() && c.pyAtPV == pos.sourcePy() && c.pzAtPV == pos.sourcePz()) ? 1.f : 0.f);
      } else if (pos.sourceV0Id() != -2 && ele.sourceV0Id() != -2) { // -2: pooled free leg, no run-1 photon
        registry.fill(HIST("Mixed/hCrossedSourceQ"), type, qInv(pos.sourcePx(), pos.sourcePy(), pos.sourcePz(), ele.sourcePx(), ele.sourcePy(), ele.sourcePz()));
        registry.fill(HIST("Mixed/hCrossedSourceDist"), type, std::hypot(pos.sourceVx() - ele.sourceVx(), pos.sourceVy() - ele.sourceVy(), pos.sourceVz() - ele.sourceVz()));
        registry.fill(HIST("Mixed/hCrossedSourceEvents"), type, pos.sourceEventUid() == ele.sourceEventUid() ? 0.f : 1.f);
      }
      float truth = 4.f; // no label
      if (pos.mcMotherId() >= 0 && ele.mcMotherId() >= 0) {
        const bool truePhoton = pos.sourceEventUid() == ele.sourceEventUid() && pos.mcMotherId() == ele.mcMotherId() && pos.mcMotherPdg() == kGamma;
        truth = (original ? 0.f : 2.f) + (truePhoton ? 0.f : 1.f);
      }
      registry.fill(HIST("Mixed/hStoredTruth"), type, truth);
      MixPairPhoton stored;
      stored.px = c.pxAtPV;
      stored.py = c.pyAtPV;
      stored.pz = c.pzAtPV;
      setMixPairLegs(stored, pos, ele);
      stored.posRow = key.posId;
      stored.eleRow = key.eleId;
      stored.original = original;
      if (const auto itGroup = mGroupOfKey.find(key); itGroup != mGroupOfKey.end()) {
        stored.group = itGroup->second;
      }
      stored.vx = c.svX;
      stored.vy = c.svY;
      stored.vz = c.svZ;
      after[key.collisionId].push_back(stored);
    }
    for (const auto& collision : collisions) {
      const float type = collision.isMixed() ? 1.f : 0.f;
      const auto it = nStored.find(collision.globalIndex());
      const std::array<int, 2> n = it == nStored.end() ? std::array<int, 2>{0, 0} : it->second;
      const auto itIn = inputOfColl.find(collision.globalIndex());
      const int nIn = itIn == inputOfColl.end() ? 0 : static_cast<int>(itIn->second.size()); // run-1 photons (not the pooled free legs)
      registry.fill(HIST("Mixed/hNStoredVsIn"), type, nIn, n[0] + n[1]);
      registry.fill(HIST("Mixed/hNOriginalLost"), type, nIn - n[0]);
      registry.fill(HIST("Mixed/hNCrossedStored"), type, n[1]);
      if (const auto itA = after.find(collision.globalIndex()); mixed.mixedDiagnostics && itA != after.end()) {
        fillMixPairs<true>(type, itA->second);
        fillSameSourcePairs(type, itA->second);
      }
    }

    if (mixed.mixedDiagnostics) { // MC diagnostics of the pairings (cost on data)
      // alternative pairings of the input pairs: passed the run-2 fit (candidate) / stored
      std::set<std::tuple<int64_t, int64_t, int64_t>> candLegs, storedLegs;
      for (const auto& c : vecV0Dedup) {
        candLegs.emplace(c.colID, c.posID, c.eleID);
      }
      for (const auto& k : stored_fullv0Ids) {
        storedLegs.emplace(k.collisionId, k.posId, k.eleId);
      }
      for (const auto& [colId, input] : inputOfColl) {
        fillMixFakeClasses(collisions.rawIteratorAt(colId).isMixed() ? 1.f : 0.f, 0.f, colId, input, candLegs);
      }
      for (const auto& [colId, out] : after) {
        fillMixFakeClasses(collisions.rawIteratorAt(colId).isMixed() ? 1.f : 0.f, 1.f, colId, out, candLegs);
      }
      // swap-symmetric counting (data-computable): the decision between the two pairings of an ambiguous 4-leg configuration no longer
      // enters the pair count, in the same and the mixed event alike
      {
        std::map<std::tuple<int64_t, int64_t, int64_t>, std::array<float, 3>> candMom; // first fit of a leg pair = its dedup candidate
        for (const auto& c : vecV0Dedup) {
          candMom.emplace(std::make_tuple(static_cast<int64_t>(c.colID), static_cast<int64_t>(c.posID), static_cast<int64_t>(c.eleID)), std::array<float, 3>{c.pxAtPV, c.pyAtPV, c.pzAtPV});
        }
        for (const auto& [colId, out] : after) {
          const float type = collisions.rawIteratorAt(colId).isMixed() ? 1.f : 0.f;
          for (size_t i = 0; i < out.size(); ++i) {
            for (size_t j = i + 1; j < out.size(); ++j) {
              const auto& a = out[i];
              const auto& b = out[j];
              if (skipSameSourcePair(type, a, b)) {
                continue;
              }
              const float q = qInv(a.px, a.py, a.pz, b.px, b.py, b.pz);
              const auto itX1 = candMom.find(std::make_tuple(static_cast<int64_t>(colId), a.posRow, b.eleRow));
              const auto itX2 = candMom.find(std::make_tuple(static_cast<int64_t>(colId), b.posRow, a.eleRow));
              const bool ambiguous = itX1 != candMom.end() && itX2 != candMom.end();
              const float qX = ambiguous ? qInv(itX1->second[0], itX1->second[1], itX1->second[2], itX2->second[0], itX2->second[1], itX2->second[2]) : 999.f;
              if (q >= 0.5f && (!ambiguous || qX >= 0.5f)) { // o2-linter: disable=magic-number (histogram range)
                continue;
              }
              const float cls = mixPairClass(a, b);
              if (ambiguous) {
                registry.fill(HIST("Mixed/Pairs/hQSym"), type, cls, q, 0.5f);
                registry.fill(HIST("Mixed/Pairs/hQSym"), type, cls, qX, 0.5f);
              } else {
                registry.fill(HIST("Mixed/Pairs/hQSym"), type, cls, q);
              }
            }
          }
        }
      }
      std::map<std::tuple<int64_t, int64_t, int64_t>, float> candScore; // first fit of a leg pair = its dedup candidate
      for (const auto& c : vecV0Dedup) {
        candScore.emplace(std::make_tuple(static_cast<int64_t>(c.colID), static_cast<int64_t>(c.posID), static_cast<int64_t>(c.eleID)), c.score);
      }
      for (const auto& [colId, input] : inputOfColl) {
        fillMixSwapBrems(collisions.rawIteratorAt(colId).isMixed() ? 1.f : 0.f, 0.f, colId, input, candScore, legs);
      }
      for (const auto& [colId, out] : after) {
        fillMixSwapBrems(collisions.rawIteratorAt(colId).isMixed() ? 1.f : 0.f, 1.f, colId, out, candScore, legs);
      }
      for (const auto& [colId, input] : inputOfColl) {
        const float type = collisions.rawIteratorAt(colId).isMixed() ? 1.f : 0.f;
        for (size_t i = 0; i < input.size(); ++i) {
          for (size_t j = i + 1; j < input.size(); ++j) {
            const auto& a = input[i];
            const auto& b = input[j];
            const int truth = mixPairTruth(a, b);
            if (truth != 1 && truth != 2) {
              continue;
            }
            const float q = qInv(a.px, a.py, a.pz, b.px, b.py, b.pz);
            if (q >= 0.1f) { // o2-linter: disable=magic-number (diagnostic range)
              continue;
            }
            const auto alt1 = std::make_tuple(colId, a.posRow, b.eleRow);
            const auto alt2 = std::make_tuple(colId, b.posRow, a.eleRow);
            const int nAlt = static_cast<int>(candLegs.contains(alt1)) + static_cast<int>(candLegs.contains(alt2));
            int outcome = 2;
            if (storedLegs.contains({colId, a.posRow, a.eleRow}) && storedLegs.contains({colId, b.posRow, b.eleRow})) {
              outcome = 0;
            } else if (storedLegs.contains(alt1) && storedLegs.contains(alt2)) {
              outcome = 1;
            }
            if (truth == 1) {
              registry.fill(HIST("Mixed/Pairs/hAltTrueTrue"), type, static_cast<float>(3 * nAlt + outcome), q);
            } else {
              registry.fill(HIST("Mixed/Pairs/hAltSwap"), type, static_cast<float>(3 * nAlt + outcome), q);
            }
          }
        }
      }
    }

    pca_map.clear();
    cospa_map.clear();
    nv0_map.clear();
    stored_v0Ids.clear();
    stored_fullv0Ids.clear();
    vecV0Dedup.clear();
    registry.fill(HIST("Mixed/hTimeDF"), std::chrono::duration<float, std::milli>(std::chrono::steady_clock::now() - tDF).count());
  }

  o2::framework::expressions::Filter v0Filter = o2::aod::v0::v0Type > static_cast<uint8_t>(0);
  using FilteredV0s = o2::soa::Filtered<o2::aod::V0s>;

  enum PairFate {
    kBefore = 0,
    kBothSurvive,
    kOneLost,
    kBothLost,
    kLostByCrossFake,
    kLostByPartnerFake,
    kLostByOtherFake,
    kLostByTrue
  };

  void fillTruthPair(int fate, float q, float dEta)
  {
    switch (fate) {
      case kBefore:
        registry.fill(HIST("MCDedup/Pairs/before/hQ"), q);
        registry.fill(HIST("MCDedup/Pairs/before/hDEta"), dEta);
        break;
      case kBothSurvive:
        registry.fill(HIST("MCDedup/Pairs/bothSurvive/hQ"), q);
        registry.fill(HIST("MCDedup/Pairs/bothSurvive/hDEta"), dEta);
        break;
      case kOneLost:
        registry.fill(HIST("MCDedup/Pairs/oneLost/hQ"), q);
        registry.fill(HIST("MCDedup/Pairs/oneLost/hDEta"), dEta);
        break;
      case kBothLost:
        registry.fill(HIST("MCDedup/Pairs/bothLost/hQ"), q);
        registry.fill(HIST("MCDedup/Pairs/bothLost/hDEta"), dEta);
        break;
      case kLostByCrossFake:
        registry.fill(HIST("MCDedup/Pairs/lostByCrossFake/hQ"), q);
        registry.fill(HIST("MCDedup/Pairs/lostByCrossFake/hDEta"), dEta);
        break;
      case kLostByPartnerFake:
        registry.fill(HIST("MCDedup/Pairs/lostByPartnerFake/hQ"), q);
        registry.fill(HIST("MCDedup/Pairs/lostByPartnerFake/hDEta"), dEta);
        break;
      case kLostByOtherFake:
        registry.fill(HIST("MCDedup/Pairs/lostByOtherFake/hQ"), q);
        registry.fill(HIST("MCDedup/Pairs/lostByOtherFake/hDEta"), dEta);
        break;
      case kLostByTrue:
        registry.fill(HIST("MCDedup/Pairs/lostByTrue/hQ"), q);
        registry.fill(HIST("MCDedup/Pairs/lostByTrue/hDEta"), dEta);
        break;
      default:
        break;
    }
  }

  void fillTruthPairMap(int fate, float dEta, float dPhi, float q)
  {
    if (!dedup.cfgDedupTruthMaps) {
      return;
    }
    switch (fate) {
      case kBefore:
        registry.fill(HIST("MCDedup/Pairs/before/hMap"), dEta, dPhi, q);
        break;
      case kBothSurvive:
        registry.fill(HIST("MCDedup/Pairs/bothSurvive/hMap"), dEta, dPhi, q);
        break;
      case kBothLost:
        registry.fill(HIST("MCDedup/Pairs/bothLost/hMap"), dEta, dPhi, q);
        break;
      case kLostByPartnerFake:
        registry.fill(HIST("MCDedup/Pairs/lostByPartnerFake/hMap"), dEta, dPhi, q);
        break;
      default:
        break;
    }
  }

  struct PhotonFate {
    bool alive = false;
    int best = -1;                // best-score candidate
    std::vector<int> cands;       // all true candidates of this photon
    std::set<int64_t> aliveColls; // collisions with a surviving candidate
  };

  // which collision the copies of a true photon end up in
  template <o2::soa::is_table TCollisions, o2::soa::is_table TMCParticles>
  void fillCopyDiagnostics(TCollisions const& collisions, TMCParticles const& mcparticles, DedupDiag const& diag, std::map<int64_t, PhotonFate> const& fate)
  {
    auto isRightColl = [&](int64_t colId, int64_t mid) {
      const auto col = collisions.rawIteratorAt(colId);
      return col.has_mcCollision() && col.mcCollisionId() == mcparticles.iteratorAt(mid).mcCollisionId();
    };
    auto dzPointing = [&](V0CandidateHelper const& c) {
      const auto col = collisions.rawIteratorAt(c.colID);
      const float pt = std::hypot(c.pxAtPV, c.pyAtPV);
      return pt > 0.f ? std::fabs(c.svZ - std::hypot(c.svX - col.posX(), c.svY - col.posY()) * c.pzAtPV / pt - col.posZ()) : 999.f; // o2-linter: disable=magic-number (no direction)
    };
    for (const auto& [mid, f] : fate) {
      std::set<int64_t> copyColls;
      bool rightAvailable = false;
      for (const int& ci : f.cands) {
        copyColls.insert(diag.snapshot[ci].colID);
        rightAvailable = rightAvailable || isRightColl(diag.snapshot[ci].colID, mid);
      }
      int choice = 3; // o2-linter: disable=magic-number (bin 3 = lost)
      if (f.alive) {
        const auto& midRef = mid; // ordinary variable, not a structured binding
        const bool storedRight = std::any_of(f.aliveColls.begin(), f.aliveColls.end(), [&](int64_t c) { return isRightColl(c, midRef); });
        choice = storedRight ? 0 : (rightAvailable ? 1 : 2);
      }
      const int nCopies = static_cast<int>(copyColls.size());
      registry.fill(HIST("MCDedup/Copies/hChoice"), static_cast<float>(std::min(nCopies, 5)), static_cast<float>(choice)); // o2-linter: disable=magic-number (last bin)
      constexpr int MinCopies = 2;
      if (nCopies < MinCopies || !rightAvailable) {
        continue;
      }
      int bestDz = -1, bestPA = -1, bestPA3D = -1;
      for (const int& ci : f.cands) {
        const auto& c = diag.snapshot[ci];
        const float right = isRightColl(c.colID, mid) ? 1.f : 0.f;
        registry.fill(HIST("MCDedup/Copies/hPAxy"), right, std::acos(std::clamp(c.cosPAxy, -1.f, 1.f)));
        registry.fill(HIST("MCDedup/Copies/hPA3D"), right, std::acos(std::clamp(c.cosPA3D, -1.f, 1.f)));
        if (bestDz < 0 || dzPointing(c) < dzPointing(diag.snapshot[bestDz])) {
          bestDz = ci;
        }
        if (bestPA < 0 || c.cosPAxy > diag.snapshot[bestPA].cosPAxy) {
          bestPA = ci;
        }
        if (bestPA3D < 0 || c.cosPA3D > diag.snapshot[bestPA3D].cosPA3D) {
          bestPA3D = ci;
        }
      }
      if (f.alive) { // same photons for all rules
        registry.fill(HIST("MCDedup/Copies/hRule"), 0.f, choice == 0 ? 1.f : 0.f);
        registry.fill(HIST("MCDedup/Copies/hRule"), 1.f, isRightColl(diag.snapshot[bestDz].colID, mid) ? 1.f : 0.f);
        registry.fill(HIST("MCDedup/Copies/hRule"), 2.f, isRightColl(diag.snapshot[bestPA].colID, mid) ? 1.f : 0.f);
        registry.fill(HIST("MCDedup/Copies/hRule"), 3.f, isRightColl(diag.snapshot[bestPA3D].colID, mid) ? 1.f : 0.f);
      }
    }
  }

  template <typename TMCParticles>
  void fillRun1FourLeg(TMCParticles const& mcparticles, DedupDiag const& diag, std::vector<uint8_t> const& cls, std::vector<int64_t> const& motherPos,
                       std::set<CandKey> const& storedSet, std::map<int64_t, int>& nCandsOfTrack)
  {
    using Key3 = std::tuple<int64_t, int64_t, int64_t>;
    const int nCand = static_cast<int>(diag.snapshot.size());
    std::map<Key3, int> firstCand;                           // first fit of a leg pair = its candidate
    std::map<std::pair<int64_t, int64_t>, int> trueOfPhoton; // (collision, MC photon) -> first true candidate
    for (int i = 0; i < nCand; ++i) {
      const auto& c = diag.snapshot[i];
      firstCand.emplace(Key3{c.colID, c.posID, c.eleID}, i);
      if (cls[i] == kV0True) {
        trueOfPhoton.emplace(std::make_pair(static_cast<int64_t>(c.colID), motherPos[i]), i);
      }
    }
    std::map<int64_t, std::vector<std::pair<int64_t, int>>> photonsOfColl;
    for (const auto& [k, i] : trueOfPhoton) {
      photonsOfColl[k.first].emplace_back(k.second, i);
    }
    auto stored = [&](int i) { return storedSet.contains(keyOf(diag.snapshot[i])); };
    const std::vector<uint64_t> legMask{0b0011, 0b1100, 0b1001, 0b0110};
    for (const auto& [colId, phs] : photonsOfColl) {
      for (size_t i = 0; i < phs.size(); ++i) {
        for (size_t j = i + 1; j < phs.size(); ++j) {
          const auto& a = diag.snapshot[phs[i].second];
          const auto& b = diag.snapshot[phs[j].second];
          if (a.posID == b.posID || a.eleID == b.eleID) {
            continue; // two MC photons sharing a track (should not happen for true candidates)
          }
          const auto itX1 = firstCand.find(Key3{colId, a.posID, b.eleID});
          const auto itX2 = firstCand.find(Key3{colId, b.posID, a.eleID});
          if (itX1 == firstCand.end() || itX2 == firstCand.end()) {
            continue;
          }
          const auto mA = mcparticles.rawIteratorAt(phs[i].first);
          const auto mB = mcparticles.rawIteratorAt(phs[j].first);
          const float qTrue = qInv(mA.px(), mA.py(), mA.pz(), mB.px(), mB.py(), mB.pz());
          if (qTrue >= 0.05f) { // o2-linter: disable=magic-number (diagnostic range)
            continue;
          }
          const int iA = phs[i].second, iB = phs[j].second, iX1 = itX1->second, iX2 = itX2->second;
          const int nT = (stored(iA) ? 1 : 0) + (stored(iB) ? 1 : 0), nX = (stored(iX1) ? 1 : 0) + (stored(iX2) ? 1 : 0);
          float outcome = 5.f;
          if (nT == 2 && nX == 0) {
            outcome = 0.f;
          } else if (nT == 0 && nX == 2) {
            outcome = 1.f;
          } else if (nT == 1 && nX == 0) {
            outcome = 2.f;
          } else if (nT == 0 && nX == 1) {
            outcome = 3.f;
          } else if (nT == 0 && nX == 0) {
            outcome = 4.f;
          }
          const std::vector<float> scores{a.score, b.score, diag.snapshot[iX1].score, diag.snapshot[iX2].score};
          const auto sol = solveConflictGroup(legMask, scores, dedup.groupMatchingScoreThreshold.value, false);
          const float isolated = sol.bestMask == 0b0011 ? 0.f : (sol.bestMask == 0b1100 ? 1.f : 2.f);                               // o2-linter: disable=magic-number (subset masks)
          const int nThird = nCandsOfTrack[a.posID] + nCandsOfTrack[a.eleID] + nCandsOfTrack[b.posID] + nCandsOfTrack[b.eleID] - 8; // o2-linter: disable=magic-number (each leg is in 2 of the 4)
          bool related = false;
          if (mA.has_mothers() && mB.has_mothers()) {
            const auto gA = mcparticles.rawIteratorAt(mA.mothersIds()[0]);
            const auto gB = mcparticles.rawIteratorAt(mB.mothersIds()[0]);
            related = gA.globalIndex() == gB.globalIndex() || std::abs(gA.pdgCode()) == kElectron || std::abs(gB.pdgCode()) == kElectron;
          }
          const bool finder = a.v0ID < 0 || b.v0ID < 0;
          registry.fill(HIST("MCDedup/FourLeg/hRun1"), qTrue, outcome, isolated, static_cast<float>(std::clamp(nThird, 0, 5)), related ? 1.f : 0.f, finder ? 1.f : 0.f); // o2-linter: disable=magic-number (last bin)
          const int crossedOrigin = (diag.snapshot[iX1].v0ID < 0 ? 1 : 0) + (diag.snapshot[iX2].v0ID < 0 ? 1 : 0);                                                       // finder off: non-SVertexer = swap resolver
          const int trueOrigin = (a.v0ID < 0 ? 1 : 0) + (b.v0ID < 0 ? 1 : 0);
          registry.fill(HIST("MCDedup/FourLeg/hRun1Origin"), qTrue, outcome, static_cast<float>(crossedOrigin), static_cast<float>(trueOrigin), nThird > 0 ? 1.f : 0.f, related ? 1.f : 0.f);
        }
      }
    }
  }

  template <o2::soa::is_table TCollisions, o2::soa::is_table TTracks, o2::soa::is_table TMCParticles>
  void fillDedupTruthDiagnostics(TCollisions const& collisions, TTracks const& tracks, TMCParticles const& mcparticles, DedupDiag const& diag)
  {
    const int nCand = static_cast<int>(diag.snapshot.size());
    if (nCand == 0) {
      return;
    }
    const std::set<CandKey> storedSet(diag.stored.begin(), diag.stored.end());

    std::vector<CandKey> key(nCand);
    std::vector<uint8_t> cls(nCand, static_cast<uint8_t>(kV0OtherFake));
    std::vector<int64_t> motherPos(nCand, -1), motherEle(nCand, -1);
    std::map<CandKey, int> indexOfKey;
    for (int i = 0; i < nCand; ++i) {
      const auto& c = diag.snapshot[i];
      key[i] = keyOf(c);
      indexOfKey[key[i]] = i;
      cls[i] = static_cast<uint8_t>(classifyV0Truth(tracks.rawIteratorAt(c.posID), tracks.rawIteratorAt(c.eleID),
                                                    mcparticles, motherPos[i], motherEle[i]));
      registry.fill(HIST("MCDedup/Candidates/hClassStage"), static_cast<float>(cls[i]), 0.f);
      if (storedSet.contains(key[i]) > 0) {
        registry.fill(HIST("MCDedup/Candidates/hClassStage"), static_cast<float>(cls[i]), 1.f);
      }
    }

    std::map<int64_t, PhotonFate> fate;
    for (int i = 0; i < nCand; ++i) {
      if (cls[i] != kV0True) {
        continue;
      }
      auto& f = fate[motherPos[i]];
      f.cands.push_back(i);
      if (storedSet.contains(key[i]) > 0) {
        f.alive = true;
        f.aliveColls.insert(key[i].collisionId);
      }
      if (f.best < 0 || diag.snapshot[i].score < diag.snapshot[f.best].score) {
        f.best = i;
      }
    }

    std::map<int64_t, int> nCandsOfTrack;
    for (int i = 0; i < nCand; ++i) {
      nCandsOfTrack[diag.snapshot[i].posID]++;
      nCandsOfTrack[diag.snapshot[i].eleID]++;
    }
    auto fillFateProps = [&](int ci, float fateBin) {
      const auto& c = diag.snapshot[ci];
      const int nITS = (tracks.rawIteratorAt(c.posID).hasITS() ? 1 : 0) + (tracks.rawIteratorAt(c.eleID).hasITS() ? 1 : 0);
      const int nOther = nCandsOfTrack[c.posID] + nCandsOfTrack[c.eleID] - 2;
      registry.fill(HIST("MCDedup/Photons/hFateProps"), fateBin, std::hypot(c.pxAtPV, c.pyAtPV), std::hypot(c.svX, c.svY), static_cast<float>(nITS),
                    std::acos(std::clamp(c.cosPA, -1.f, 1.f)), c.pca, static_cast<float>(std::min(nOther, 10))); // o2-linter: disable=magic-number (last bin)
    };

    std::map<int64_t, std::vector<int>> blockersOfPhoton;
    for (const auto& [mid, f] : fate) {
      if (f.alive) {
        registry.fill(HIST("MCDedup/Photons/hFate"), 0.f);
        const auto itStored = std::find_if(f.cands.begin(), f.cands.end(), [&](int ci) { return storedSet.contains(key[ci]); });
        fillFateProps(itStored != f.cands.end() ? *itStored : f.best, 0.f);
        continue;
      }
      auto& blockers = blockersOfPhoton[mid];
      for (const int& ci : f.cands) {
        const auto itBlocker = diag.blockersByKey.find(key[ci]);
        if (itBlocker == diag.blockersByKey.end()) {
          continue;
        }
        for (const auto& bKey : itBlocker->second) {
          const auto itIndex = indexOfKey.find(bKey);
          if (itIndex == indexOfKey.end()) {
            continue;
          }
          if (std::find(blockers.begin(), blockers.end(), itIndex->second) == blockers.end()) {
            blockers.push_back(itIndex->second);
          }
        }
      }
      // Why was it not in the winning set? Take the most favourable of its candidates: the one
      // that could have been kept at the smallest cost in cardinality.
      int bestDeficit = -1;
      float bestExcess = 0.f;
      for (const int& ci : f.cands) {
        const auto itMargin = diag.lossMargin.find(key[ci]);
        if (itMargin == diag.lossMargin.end()) {
          continue;
        }
        if (bestDeficit < 0 || itMargin->second.first < bestDeficit ||
            (itMargin->second.first == bestDeficit && itMargin->second.second < bestExcess)) {
          bestDeficit = itMargin->second.first;
          bestExcess = itMargin->second.second;
        }
      }
      if (bestDeficit >= 0) {
        registry.fill(HIST("MCDedup/Photons/hLostMargin"), static_cast<float>(bestDeficit), bestExcess);
      }
      // hVictimBlockerClass below counts one entry per blocker and a victim can have two.
      registry.fill(HIST("MCDedup/Photons/hFate"), blockers.empty() ? 2.f : 1.f);
      const bool lostToCross = std::any_of(blockers.begin(), blockers.end(), [&](int w) { return cls[w] == kV0CrossLegFake; });
      if (f.best >= 0) {
        fillFateProps(f.best, blockers.empty() ? 3.f : (lostToCross ? 1.f : 2.f));
      }

      if (blockers.empty()) {
        registry.fill(HIST("MCDedup/Photons/hBlockerClass"), 3.f); // o2-linter: disable=magic-number (bin 3 = no blocker)
        continue;
      }
      for (const int& w : blockers) {
        registry.fill(HIST("MCDedup/Photons/hBlockerClass"), static_cast<float>(cls[w]));
      }
      // margin and PCA comparison relation the best candidate to its blocker

      if (f.best >= 0 && dedup.deduplicationMode.value < static_cast<int>(V0DeduplicationMode::GroupMatching)) {
        const int w = blockers.front();
        registry.fill(HIST("MCDedup/Photons/hKillMargin"), diag.snapshot[f.best].score - diag.snapshot[w].score);
        registry.fill(HIST("MCDedup/Photons/hVictimVsWinnerPCA"), diag.snapshot[f.best].pca, diag.snapshot[w].pca);
      }
    }

    fillRun1FourLeg(mcparticles, diag, cls, motherPos, storedSet, nCandsOfTrack);

    {
      std::map<std::pair<int64_t, int64_t>, int> storedOfTrack; // (collision, track) -> stored candidate
      for (int i = 0; i < nCand; ++i) {
        if (storedSet.contains(key[i])) {
          const auto& c = diag.snapshot[i];
          storedOfTrack[{c.colID, c.posID}] = i;
          storedOfTrack[{c.colID, c.eleID}] = i;
        }
      }
      const float s0 = dedup.groupMatchingScoreThreshold.value;
      auto contrib = [&](float sc) { return std::max(s0 - sc, 1e-3f); }; // o2-linter: disable=magic-number (objective floor)
      for (const auto& [mid, f] : fate) {
        if (f.best < 0) {
          continue;
        }
        const auto& t = diag.snapshot[f.best];
        const float pt = std::hypot(t.pxAtPV, t.pyAtPV);
        if (f.alive) {
          registry.fill(HIST("MCDedup/HalfSwap/hSource"), 0.f, t.score, 0.f, 0.f, pt);
          continue;
        }
        const auto ip = storedOfTrack.find({t.colID, t.posID});
        const auto ie = storedOfTrack.find({t.colID, t.eleID});
        const bool hp = ip != storedOfTrack.end(), he = ie != storedOfTrack.end();
        float cat = 3.f;
        if (hp && he) {
          cat = ip->second != ie->second ? 1.f : 4.f;
        } else if (hp || he) {
          cat = 2.f;
        }
        float dObj = 0.f, dRaw = 0.f;
        if (cat == 1.f) {
          const auto& w1 = diag.snapshot[ip->second];
          const auto& w2 = diag.snapshot[ie->second];
          dObj = contrib(w1.score) + contrib(w2.score) - contrib(t.score);
          dRaw = 0.5f * (w1.score + w2.score) - t.score;
        }
        registry.fill(HIST("MCDedup/HalfSwap/hSource"), cat, t.score, std::clamp(dObj, -0.499f, 0.499f), std::clamp(dRaw, -0.999f, 0.999f), pt);
      }
    }

    if (dedup.cfgExtendedDiag) {
      fillCopyDiagnostics(collisions, mcparticles, diag, fate);
    }

    // stored cross-leg fakes by collisions
    std::map<std::pair<int64_t, int64_t>, std::set<int64_t>> storedCrossColl;
    for (int i = 0; i < nCand; ++i) {
      if (cls[i] == kV0CrossLegFake && storedSet.contains(key[i])) {
        storedCrossColl[{motherPos[i], motherEle[i]}].insert(key[i].collisionId);
      }
    }

    //  truth pairs, grouped by MC collision
    std::map<int, std::vector<int64_t>> mothersByMcColl;
    for (const auto& [mid, f] : fate) {
      mothersByMcColl[mcparticles.iteratorAt(mid).mcCollisionId()].push_back(mid);
    }
    constexpr float MaxQTruth = 0.3f;
    for (const auto& [mcCol, mids] : mothersByMcColl) {
      for (size_t a = 0; a < mids.size(); ++a) {
        for (size_t b = a + 1; b < mids.size(); ++b) {
          const auto g1 = mcparticles.iteratorAt(mids[a]);
          const auto g2 = mcparticles.iteratorAt(mids[b]);
          const float e1 = std::hypot(g1.px(), g1.py(), g1.pz());
          const float e2 = std::hypot(g2.px(), g2.py(), g2.pz());
          const float qTrue = std::sqrt(std::max(0.f, 2.f * (e1 * e2 - g1.px() * g2.px() - g1.py() * g2.py() - g1.pz() * g2.pz())));
          if (qTrue > MaxQTruth) {
            continue;
          }
          const float dEta = g1.eta() - g2.eta();
          const float dPhi = RecoDecay::constrainAngle(static_cast<float>(g1.phi() - g2.phi()), -o2::constants::math::PI);
          fillTruthPair(kBefore, qTrue, dEta);
          fillTruthPairMap(kBefore, dEta, dPhi, qTrue);

          const auto& f1 = fate.at(mids[a]);
          const auto& f2 = fate.at(mids[b]);
          const int nLost = (f1.alive ? 0 : 1) + (f2.alive ? 0 : 1);
          if (nLost == 0) {
            fillTruthPair(kBothSurvive, qTrue, dEta);
            fillTruthPairMap(kBothSurvive, dEta, dPhi, qTrue);
            bool sameColl = false;
            for (const auto& c1 : f1.aliveColls) {
              sameColl = sameColl || (f2.aliveColls.count(c1) > 0);
            }
            if (sameColl) {
              registry.fill(HIST("MCDedup/Pairs/bothSurviveSameColl/hQ"), qTrue);
            }
            continue;
          }
          { // the pair stored as the crossed pairing (A+B-, B+A-) in one collision or one crossed pair only
            const auto itAB = storedCrossColl.find({mids[a], mids[b]});
            const auto itBA = storedCrossColl.find({mids[b], mids[a]});
            bool bothCrossed = false;
            if (itAB != storedCrossColl.end() && itBA != storedCrossColl.end()) {
              for (const auto& c : itAB->second) {
                bothCrossed = bothCrossed || itBA->second.contains(c);
              }
            }
            if (bothCrossed) {
              registry.fill(HIST("MCDedup/Pairs/surviveAsCrossedPairing/hQ"), qTrue);
            } else if (itAB != storedCrossColl.end() || itBA != storedCrossColl.end()) {
              registry.fill(HIST("MCDedup/Pairs/halfCrossedStored/hQ"), qTrue);
            }
          }
          if (nLost == 1) {
            fillTruthPair(kOneLost, qTrue, dEta);
          } else {
            fillTruthPair(kBothLost, qTrue, dEta);
            fillTruthPairMap(kBothLost, dEta, dPhi, qTrue);
          }

          // Reason for photon building failure
          bool byCross = false, byPartner = false, byOther = false, byTrue = false;
          auto attribute = [&](int64_t victim, int64_t partner) {
            if (fate.at(victim).alive) {
              return;
            }
            const auto itB = blockersOfPhoton.find(victim);
            if (itB == blockersOfPhoton.end()) {
              return;
            }
            for (const int& w : itB->second) {
              if (cls[w] == kV0CrossLegFake) {
                byCross = true;
                // partner fake = a cross-leg fake using a leg of the other photon
                if (motherPos[w] == partner || motherEle[w] == partner) {
                  byPartner = true;
                }
              } else if (cls[w] == kV0OtherFake) {
                byOther = true;
              } else if (cls[w] == kV0True) {
                byTrue = true; // real photon against real photon
              }
            }
          };
          attribute(mids[a], mids[b]);
          attribute(mids[b], mids[a]);
          if (byCross) {
            fillTruthPair(kLostByCrossFake, qTrue, dEta);
          }
          if (byPartner) {
            fillTruthPair(kLostByPartnerFake, qTrue, dEta);
            fillTruthPairMap(kLostByPartnerFake, dEta, dPhi, qTrue);
          }
          if (byOther) {
            fillTruthPair(kLostByOtherFake, qTrue, dEta);
          }
          if (byTrue) {
            fillTruthPair(kLostByTrue, qTrue, dEta);
          }
        }
      }
    }

    {
      struct Elig {
        int64_t pos = -1, ele = -1;
        int nIts = 0;
        float vx = 0.f, vy = 0.f;
      };
      std::map<int64_t, Elig> elig;
      for (const auto& track : tracks) {
        if (track.pt() < trackcuts.min_pt_trackiu || !checkV0leg<true>(track)) {
          continue;
        }
        const int64_t mid = getPhotonMotherId(track, mcparticles);
        if (mid < 0) {
          continue;
        }
        auto& e = elig[mid];
        if (track.sign() > 0) {
          e.pos = track.globalIndex();
        } else {
          e.ele = track.globalIndex();
        }
        if (track.hasITS()) {
          e.nIts++;
        }
        const auto mc = track.template mcParticle_as<TMCParticles>();
        e.vx = mc.vx();
        e.vy = mc.vy();
      }
      std::map<int64_t, int> crossedStoredOfPhoton;
      std::map<std::pair<int64_t, int64_t>, std::pair<int, int>> crossedOfPair;    // (offered, stored)
      std::map<std::pair<int64_t, int64_t>, std::vector<int>> crossedOriginOfPair; // stored crossed: 0 SVertexer, 1 finder, 2 swap resolver
      std::map<std::pair<int64_t, int64_t>, std::vector<int>> crossedIndexOfPair;  // stored crossed: candidate index
      std::map<std::pair<CandKey, CandKey>, const SwapRecord*> swapRecordOf;
      for (const auto& r : diag.swaps) {
        swapRecordOf[std::minmax(r.keyI, r.keyJ)] = &r;
      }
      for (int i = 0; i < nCand; ++i) {
        if (cls[i] != kV0CrossLegFake) {
          continue;
        }
        auto& cp = crossedOfPair[std::minmax(motherPos[i], motherEle[i])];
        cp.first++;
        if (storedSet.contains(key[i])) {
          cp.second++;
          crossedOriginOfPair[std::minmax(motherPos[i], motherEle[i])].push_back(key[i].v0Id >= 0 ? 0 : (diag.finderKeys.contains(key[i]) ? 1 : 2));
          crossedIndexOfPair[std::minmax(motherPos[i], motherEle[i])].push_back(i);
          crossedStoredOfPhoton[motherPos[i]]++;
          crossedStoredOfPhoton[motherEle[i]]++;
        }
      }
      std::map<int, std::vector<int64_t>> eligByMcColl;
      for (const auto& [mid, e] : elig) {
        if (e.pos >= 0 && e.ele >= 0) {
          eligByMcColl[mcparticles.iteratorAt(mid).mcCollisionId()].push_back(mid);
        }
      }
      auto history = [&](int64_t mid, int& detail) { // run-1 history of one true photon, see MCDedup/Inventory/hSwapOrigin
        detail = -10;
        if (fate.contains(mid)) {
          return 0;
        }
        const auto& e = elig.at(mid);
        const auto it = diag.trueLegTrace.find({e.pos, e.ele});
        const TrueLegTrace t = it == diag.trueLegTrace.end() ? TrueLegTrace{} : it->second;
        if (t.svLine > 0) {
          detail = t.svLine;
          return 1;
        }
        if (t.svLine == 0) {
          return 6; // o2-linter: disable=magic-number (history codes)
        }
        detail = t.finderStep;
        if (e.nIts == 2 || t.finderStep == -8) { // o2-linter: disable=magic-number (both legs with ITS)
          detail = -8;                           // o2-linter: disable=magic-number (trace code)
          return 2;
        }
        if (t.finderStep == FinderStep::kFindCollision) {
          detail = t.finderLine;
          return t.finderLine > 0 ? 5 : 6; // o2-linter: disable=magic-number (history codes)
        }
        return t.finderStep >= -1 ? 4 : 3; // o2-linter: disable=magic-number (history codes)
      };
      auto photonFate = [&](int64_t mid) { // 0 no true candidate, 1 candidate lost, 2 stored, 3 lost but a leg in a stored crossed candidate
        const auto it = fate.find(mid);
        const bool inCrossed = crossedStoredOfPhoton.contains(mid);
        if (it == fate.end()) {
          return inCrossed ? 3 : 0;
        }
        if (it->second.alive) {
          return 2;
        }
        return inCrossed ? 3 : 1;
      };
      for (const auto& [mcCol, mids] : eligByMcColl) {
        std::vector<float> closest(mids.size(), 0.0999f);
        for (size_t a = 0; a < mids.size(); ++a) {
          const auto ga = mcparticles.iteratorAt(mids[a]);
          for (size_t b = a + 1; b < mids.size(); ++b) {
            const auto gb = mcparticles.iteratorAt(mids[b]);
            const float qTrue = qInv(ga.px(), ga.py(), ga.pz(), gb.px(), gb.py(), gb.pz());
            closest[a] = std::min(closest[a], qTrue);
            closest[b] = std::min(closest[b], qTrue);
            if (qTrue >= 0.1f) {
              continue;
            }
            const float e1 = ga.p(), e2 = gb.p();
            const float asym = (e1 + e2) > 0.f ? std::fabs(e1 - e2) / (e1 + e2) : 0.f;
            const float kt = 0.5f * std::hypot(ga.px() + gb.px(), ga.py() + gb.py());
            const int fa = photonFate(mids[a]), fb = photonFate(mids[b]);
            const auto itx = crossedOfPair.find(std::minmax(mids[a], mids[b]));
            const int nCrossOff = itx == crossedOfPair.end() ? 0 : std::min(itx->second.first, 2);
            const int nCrossSt = itx == crossedOfPair.end() ? 0 : itx->second.second;
            const int nStoredTrue = (fa == 2 ? 1 : 0) + (fb == 2 ? 1 : 0);
            int pairFate = 5;
            if (nStoredTrue == 2) {
              pairFate = 0;
            } else if (nStoredTrue == 1 && nCrossSt == 0) {
              pairFate = 1;
            } else if (nStoredTrue == 0 && nCrossSt == 0) {
              pairFate = 2;
            } else if (nStoredTrue == 0 && nCrossSt == 2) {
              pairFate = 3;
            } else if (nCrossSt == 1) {
              pairFate = 4;
            }
            const int nTrueCand = (fate.contains(mids[a]) ? 1 : 0) + (fate.contains(mids[b]) ? 1 : 0);
            registry.fill(HIST("MCDedup/Inventory/hPair"), qTrue, std::min(asym, 0.999f), std::min(kt, 1.499f), static_cast<float>(pairFate), static_cast<float>(nTrueCand), static_cast<float>(nCrossOff));
            int detA = -10, detB = -10;
            int hA = history(mids[a], detA), hB = history(mids[b], detB);
            if (qTrue < 0.02f) { // o2-linter: disable=magic-number (close pairs)
              if (hA > 0) {
                registry.fill(HIST("MCDedup/Inventory/hNeverDetail"), static_cast<float>(pairFate), static_cast<float>(hA), static_cast<float>(detA));
              }
              if (hB > 0) {
                registry.fill(HIST("MCDedup/Inventory/hNeverDetail"), static_cast<float>(pairFate), static_cast<float>(hB), static_cast<float>(detB));
              }
            }
            if (pairFate == 3 || pairFate == 4) { // o2-linter: disable=magic-number (swap, one crossed stored)
              auto o = crossedOriginOfPair.at(std::minmax(mids[a], mids[b]));
              const int o1 = std::min(o[0], o.size() > 1 ? o[1] : o[0]), o2 = std::max(o[0], o.size() > 1 ? o[1] : o[0]);
              static constexpr std::array<std::array<int, 3>, 3> OriginCode{{{0, 1, 3}, {1, 2, 4}, {3, 4, 5}}};
              if (hA < hB) {
                std::swap(hA, hB);
              }
              registry.fill(HIST("MCDedup/Inventory/hSwapOrigin"), qTrue, pairFate == 3 ? 0.f : 1.f, static_cast<float>(nTrueCand), static_cast<float>(OriginCode[o1][o2]), static_cast<float>(hA), static_cast<float>(hB));
            }
            if (pairFate == 3) { // o2-linter: disable=magic-number (full swap): what the swap resolver did with the two stored crossed photons
              const auto& ix = crossedIndexOfPair.at(std::minmax(mids[a], mids[b]));
              const auto& x1 = diag.snapshot[ix[0]];
              const auto& x2 = diag.snapshot[ix.size() > 1 ? ix[1] : ix[0]];
              int outcome = 8; // o2-linter: disable=magic-number (outcome codes)
              float gain = -1.999f;
              int line = 0;
              const float qRec = qInv(x1.pxAtPV, x1.pyAtPV, x1.pzAtPV, x2.pxAtPV, x2.pyAtPV, x2.pzAtPV);
              if (x1.colID != x2.colID) {
                outcome = 0;
              } else if (!recovery.runSwapResolver.value) {
                outcome = 8; // o2-linter: disable=magic-number (outcome codes)
              } else if (qRec > recovery.swapMaxQ) {
                outcome = 1;
              } else if (std::hypot(x1.svX - x2.svX, x1.svY - x2.svY, x1.svZ - x2.svZ) > recovery.swapMaxConvDist) {
                outcome = 2;
              } else if (const auto it = swapRecordOf.find(std::minmax(keyOf(x1), keyOf(x2))); it != swapRecordOf.end()) {
                const auto& r = *it->second;
                if (r.complete) {
                  gain = std::clamp(-r.gain, -1.999f, 1.999f); // crossed (stored) minus true (alternative)
                  outcome = r.gain > 0.f ? 7 : 6;              // o2-linter: disable=magic-number (outcome codes)
                } else if ((r.fail & 3) != 0) {                // o2-linter: disable=magic-number (bit mask)
                  outcome = 5;                                 // o2-linter: disable=magic-number (outcome codes)
                } else if ((r.fail & 48) != 0) {               // o2-linter: disable=magic-number (bit mask)
                  outcome = 4;                                 // o2-linter: disable=magic-number (outcome codes)
                  line = r.line;
                } else {
                  outcome = 3;
                }
              }
              registry.fill(HIST("MCDedup/Inventory/hSwapLoss"), qTrue, static_cast<float>(nTrueCand), static_cast<float>(outcome), gain);
              if (outcome == 4 && qTrue < 0.02f) { // o2-linter: disable=magic-number (close pairs)
                registry.fill(HIST("MCDedup/Inventory/hSwapLossCut"), static_cast<float>(nTrueCand), static_cast<float>(line));
              }
            }
          }
        }
        for (size_t a = 0; a < mids.size(); ++a) {
          const auto ga = mcparticles.iteratorAt(mids[a]);
          const auto& e = elig.at(mids[a]);
          registry.fill(HIST("MCDedup/Inventory/hPhoton"), std::min(ga.p(), 1.999f), std::min(std::hypot(e.vx, e.vy), 89.9f), static_cast<float>(e.nIts), static_cast<float>(photonFate(mids[a])), closest[a]);
        }
      }
    }
  }

  template <o2::soa::is_table TTracks, o2::soa::is_table TMCParticles>
  void fillExtendedDiagnostics(TTracks const& tracks, TMCParticles const& mcparticles, DedupDiag const& diag)
  {
    struct Stored {
      CandKey key;
      int cls = 0;
      int64_t mPos = -1;
      int64_t mEle = -1;
      std::array<float, 3> p{};
    };
    std::map<CandKey, const V0CandidateHelper*> candOfKey;
    std::set<std::tuple<int64_t, int64_t, int64_t>> candLegs; // (collision, pos, ele) of all dedup candidates
    for (const auto& c : diag.snapshot) {
      candOfKey[keyOf(c)] = &c;
      candLegs.emplace(c.colID, c.posID, c.eleID);
    }
    std::map<int64_t, std::vector<Stored>> storedOfColl;
    for (const auto& key : diag.stored) {
      const auto it = candOfKey.find(key);
      if (it == candOfKey.end()) {
        continue;
      }
      Stored s{key, 0, -1, -1, {it->second->pxAtPV, it->second->pyAtPV, it->second->pzAtPV}};
      s.cls = classifyV0Truth(tracks.rawIteratorAt(key.posId), tracks.rawIteratorAt(key.eleId), mcparticles, s.mPos, s.mEle);
      const float pt = std::hypot(s.p[0], s.p[1]);
      registry.fill(HIST("MCDedup/Stored/hClassSource"), static_cast<float>(s.cls), key.v0Id < 0 ? 1.f : 0.f, pt);
      storedOfColl[key.collisionId].push_back(s);
    }
    constexpr float MaxQ = 0.1f;
    for (const auto& [col, photons] : storedOfColl) {
      for (size_t i = 0; i < photons.size(); ++i) {
        for (size_t j = i + 1; j < photons.size(); ++j) {
          const auto& a = photons[i];
          const auto& b = photons[j];
          const float q = qInv(a.p[0], a.p[1], a.p[2], b.p[0], b.p[1], b.p[2]);
          if (q > MaxQ) {
            continue;
          }
          int cat = 2;
          if (a.cls == kV0True && b.cls == kV0True && a.mPos != b.mPos) {
            cat = 0;
          } else if (a.cls == kV0CrossLegFake && b.cls == kV0CrossLegFake && a.mPos == b.mEle && a.mEle == b.mPos) {
            cat = 1;
            const int nTrue = static_cast<int>(candLegs.contains({col, a.key.posId, b.key.eleId})) + static_cast<int>(candLegs.contains({col, b.key.posId, a.key.eleId}));
            registry.fill(HIST("MCDedup/Stored/hSwapTrueCands"), static_cast<float>(nTrue), q);
          }
          registry.fill(HIST("MCDedup/Stored/hPairClass"), static_cast<float>(cat), q);
        }
      }
    }
    for (const auto& s : diag.swaps) {
      const std::array<std::pair<int64_t, int64_t>, 4> legs = {{{s.keyI.posId, s.keyI.eleId}, {s.keyJ.posId, s.keyJ.eleId}, {s.keyI.posId, s.keyJ.eleId}, {s.keyJ.posId, s.keyI.eleId}}};
      std::array<int, 4> isTrue{};
      for (int k = 0; k < 4; ++k) { // o2-linter: disable=magic-number (four leg pairs)
        int64_t mPos = -1, mEle = -1;
        isTrue[k] = classifyV0Truth(tracks.rawIteratorAt(legs[k].first), tracks.rawIteratorAt(legs[k].second), mcparticles, mPos, mEle) == kV0True ? 1 : 0;
      }
      registry.fill(HIST("MCDedup/Swap/hCategory"), static_cast<float>(3 * (isTrue[0] + isTrue[1]) + isTrue[2] + isTrue[3]), s.swapped ? 1.f : 0.f, s.q);
      int outcome = 3; // o2-linter: disable=magic-number (outcome codes)
      if (s.swapped) {
        outcome = 4; // o2-linter: disable=magic-number (outcome codes)
      } else if (!s.complete) {
        outcome = (s.fail & 3) != 0 ? 2 : ((s.fail & 48) != 0 ? 1 : 0); // o2-linter: disable=magic-number (bit mask)
      }
      registry.fill(HIST("MCDedup/Swap/hResolverTruth"), static_cast<float>(isTrue[0] + isTrue[1]), static_cast<float>(isTrue[2] + isTrue[3]), static_cast<float>(outcome), s.complete ? std::clamp(s.gain, -1.999f, 1.999f) : -1.999f, std::min(s.q, 0.0999f));
    }
  }

  // ===== process functions =====
  void processRec(MyCollisions const& collisions, FilteredV0s const& v0s, MyTracksIU const& tracks, MyBCs const& bcs)
  {
    build<false, false, false>(collisions, v0s, tracks, bcs);
  }
  PROCESS_SWITCH(PhotonConversionBuilder, processRec, "process reconstructed info for data", true);

  // void processRec_SWT(MyCollisionsWithSWT const& collisions, FilteredV0s const& v0s, MyTracksIU const& tracks, MyBCs const& bcs)
  // {
  //   build<false, true, false>(collisions, v0s, tracks, bcs);
  // }
  // PROCESS_SWITCH(PhotonConversionBuilder, processRec_SWT, "process reconstructed info for data", false);

  void processMC(MyCollisionsMC const& collisions, FilteredV0s const& v0s, MyTracksIUMC const& tracks,
                 MyBCs const& bcs, o2::aod::McParticles const& mcparticles)
  {
    DedupDiag diag;
    build<true, false, false>(collisions, v0s, tracks, bcs, &diag);
    fillDedupTruthDiagnostics(collisions, tracks, mcparticles, diag);
    if (dedup.cfgExtendedDiag) {
      fillExtendedDiagnostics(tracks, mcparticles, diag);
    }
  }
  PROCESS_SWITCH(PhotonConversionBuilder, processMC, "process reconstructed info for MC", false);

  void processRec_OnlyIfDielectron(o2::soa::Join<MyCollisions, o2::aod::EMEventsNee> const& collisions, FilteredV0s const& v0s, MyTracksIU const& tracks, MyBCs const& bcs)
  {
    build<false, false, true>(collisions, v0s, tracks, bcs);
  }
  PROCESS_SWITCH(PhotonConversionBuilder, processRec_OnlyIfDielectron, "process reconstructed info for data", false);

  using MixCollisionsIn = std::conditional_t<kMixedInstance, o2::aod::PCMMixCollisions, o2::aod::Collisions>;
  using MixLegsIn = std::conditional_t<kMixedInstance, o2::aod::PCMMixLegs, o2::aod::BCs>;
  std::vector<decltype(std::declval<MixLegsIn const&>().rawIteratorAt(0))> mLegIt; // run 2: one iterator per PCMMixLegs row of the DF (buildMixed)
  void processMixed(MixCollisionsIn const& collisions, MixLegsIn const& legs, MyBCs const& bcs)
  {
    if constexpr (kMixedInstance) {
      if (collisions.size() == 0 || bcs.size() == 0) {
        return;
      }
      const auto bc = bcs.begin(); // B field / material of this DF's run
      initCCDB(bc);
      updateCCDB(bc);
      buildMixed(collisions, legs);
    } else {
      LOG(fatal) << "processMixed only in photon-conversion-builder-mixed (o2-analysis-em-photon-conversion-builder-mixed)";
    }
  }
  PROCESS_SWITCH(PhotonConversionBuilder, processMixed, "builder run 2 of the co-reconstructed mixed event (input from photon-mixed-event-producer), own instance only", false);

  void processMixedDerived(MixCollisionsIn const& collisions, MixLegsIn const& legs)
  {
    if constexpr (kMixedInstance) {
      if (collisions.size() == 0) {
        return;
      }
      const auto event = collisions.begin(); // B field / material of this DF's run
      const auto* grp = ccdb->template getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", event.timestamp());
      if (grp == nullptr) {
        LOG(fatal) << "processMixedDerived: no GRPMagField for timestamp " << event.timestamp();
      }
      initCCDBFor(event.runNumber(), event.timestamp(), *grp);
      updateCCDBAt(event.timestamp());
      buildMixed(collisions, legs);
    } else {
      LOG(fatal) << "processMixedDerived only in photon-conversion-builder-mixed (o2-analysis-em-photon-conversion-builder-mixed)";
    }
  }
  PROCESS_SWITCH(PhotonConversionBuilder, processMixedDerived, "builder run 2 on derived data (input from photon-mixed-event-producer processDerived), own instance only", false);

  // void processRec_SWT_OnlyIfDielectron(soa::Join<MyCollisionsWithSWT, aod::EMEventsNee> const& collisions, FilteredV0s const& v0s, MyTracksIU const& tracks, MyBCs const& bcs)
  // {
  //   build<false, true, true>(collisions, v0s, tracks, bcs);
  // }
  // PROCESS_SWITCH(PhotonConversionBuilder, processRec_SWT_OnlyIfDielectron, "process reconstructed info for data", false);
};

#endif // PWGEM_PHOTONMESON_TABLEPRODUCER_PHOTONCONVERSIONBUILDER_H_
