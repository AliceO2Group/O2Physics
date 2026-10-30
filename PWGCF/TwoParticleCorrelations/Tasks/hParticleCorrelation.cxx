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
///
/// \file hParticleCorrelation.cxx
/// \brief Task for same- and mixed-event azimuthal angular correlations of trigger hadrons with charged hadrons, identified particles, and resonances
///
/// \author Rahul Verma (rahul.verma@iitb.ac.in) :: Durgesh Bhatt (durgesh.bhatt@cern.ch) :: Sadhana Dash (sadhana@phy.iitb.ac.in)

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/SliceCache.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics; // for constants

namespace o2::aod
{
namespace resonancecndt
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                         //! Candidate to collision
DECLARE_SOA_INDEX_COLUMN_FULL(PosTrack, posTrack, int, Tracks, "_Pos"); //! Positive track
DECLARE_SOA_INDEX_COLUMN_FULL(NegTrack, negTrack, int, Tracks, "_Neg"); //! Negative track

// Candidate kinematics
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);

DECLARE_SOA_COLUMN(Px, px, float);
DECLARE_SOA_COLUMN(Py, py, float);
DECLARE_SOA_COLUMN(Pz, pz, float);

// Invariant-mass hypotheses
DECLARE_SOA_COLUMN(MPhi1020, mPhi1020, float);
DECLARE_SOA_COLUMN(MKStar892, mKStar892, float);
DECLARE_SOA_COLUMN(MKStar892Bar, mKStar892Bar, float);
DECLARE_SOA_COLUMN(MLambda1520, mLambda1520, float);
DECLARE_SOA_COLUMN(MLambda1520Bar, mLambda1520Bar, float);

// Mass-region tags
DECLARE_SOA_COLUMN(Phi1020Tag, phi1020Tag, uint8_t);
DECLARE_SOA_COLUMN(KStar892Tag, kStar892Tag, uint8_t);
DECLARE_SOA_COLUMN(KStar892BarTag, kStar892BarTag, uint8_t);
DECLARE_SOA_COLUMN(Lambda1520Tag, lambda1520Tag, uint8_t);
DECLARE_SOA_COLUMN(Lambda1520BarTag, lambda1520BarTag, uint8_t);

} // namespace resonancecndt

DECLARE_SOA_TABLE(ResonanceCndts, "AOD", "RESONCNDTS", o2::soa::Index<>,
                  o2::aod::resonancecndt::CollisionId,
                  o2::aod::resonancecndt::PosTrackId,
                  o2::aod::resonancecndt::NegTrackId,

                  o2::aod::resonancecndt::Pt,
                  o2::aod::resonancecndt::Eta,
                  o2::aod::resonancecndt::Phi,

                  o2::aod::resonancecndt::Px,
                  o2::aod::resonancecndt::Py,
                  o2::aod::resonancecndt::Pz,

                  o2::aod::resonancecndt::MPhi1020,
                  o2::aod::resonancecndt::MKStar892,
                  o2::aod::resonancecndt::MKStar892Bar,
                  o2::aod::resonancecndt::MLambda1520,
                  o2::aod::resonancecndt::MLambda1520Bar,

                  o2::aod::resonancecndt::Phi1020Tag,
                  o2::aod::resonancecndt::KStar892Tag,
                  o2::aod::resonancecndt::KStar892BarTag,
                  o2::aod::resonancecndt::Lambda1520Tag,
                  o2::aod::resonancecndt::Lambda1520BarTag);

//______________________________________________________________________________
// Derived correlation tables
//______________________________________________________________________________

namespace hpcorrcollision
{
DECLARE_SOA_COLUMN(Centrality, centrality, float); //! Value of the mixing estimator chosen in the same-event task
DECLARE_SOA_COLUMN(CorrMask, corrMask, uint64_t);  //! Mixing pools this collision belongs to (CorrPresenceMask bits)
} // namespace hpcorrcollision

DECLARE_SOA_TABLE(HPCorrColls, "AOD", "HPCORRCOLL", o2::soa::Index<>,
                  o2::aod::collision::PosZ,
                  hpcorrcollision::Centrality,
                  hpcorrcollision::CorrMask);

//______________________________________________________________________________
// Trigger tracks
namespace hpcorrtrigger
{
DECLARE_SOA_INDEX_COLUMN(HPCorrColl, hpCorrColl); // o2-linter: disable=name/o2-column (Keep the derived collision-table relation)
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
} // namespace hpcorrtrigger

DECLARE_SOA_TABLE(HPCorrTriggers, "AOD", "HPCORRTRIG", o2::soa::Index<>,
                  hpcorrtrigger::HPCorrCollId,
                  hpcorrtrigger::Pt,
                  hpcorrtrigger::Eta,
                  hpcorrtrigger::Phi);

//______________________________________________________________________________
// Associated particles: tracks, identified hadrons and resonance candidates
namespace hpcorrassoc
{
DECLARE_SOA_INDEX_COLUMN(HPCorrColl, hpCorrColl); // o2-linter: disable=name/o2-column (Keep the derived collision-table relation)
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Sign, sign, int8_t);           //! Tracks: +1 (h+, pi+, K+, p) / -1 (h-, pi-, K-, pbar). Resonances: 0 (species, incl. anti-K*0 / anti-Lambda(1520), is in UsedMask)
DECLARE_SOA_COLUMN(UsedMask, usedMask, uint64_t); //! CorrPresenceMask bits of every correlation this particle was used in: 0-7 tracks (h, pi, K, p x low/high pT), 18-47 resonances (species x pT x mass window)
} // namespace hpcorrassoc

DECLARE_SOA_TABLE(HPCorrAssocs, "AOD", "HPCORRASSOC", o2::soa::Index<>,
                  hpcorrassoc::HPCorrCollId,
                  hpcorrassoc::Pt,
                  hpcorrassoc::Eta,
                  hpcorrassoc::Phi,
                  hpcorrassoc::Sign,
                  hpcorrassoc::UsedMask);
} // namespace o2::aod

// Delta phi in [-pi/2, 3pi/2)
float computeDeltaPhi(float particlePhi, float leadPhi)
{
  return RecoDecay::constrainAngle(particlePhi - leadPhi, -o2::constants::math::PIHalf);
}

template <typename T>
float computeRapidity(const T& track, const float mass)
{
  const float energy = RecoDecay::e(track.px(), track.py(), track.pz(), mass);
  return 0.5f * std::log((energy + track.pz()) / (energy - track.pz()));
}

template <typename T>
inline T getCfg(const auto& cfgLabelledArray, int32_t row, int32_t col)
{
  double val = cfgLabelledArray->get(row, col); // Always stored as double

  if constexpr (std::is_same_v<T, bool>) {
    return val != 0.0; // will return anything nonzero as true
  } else if constexpr (std::is_same_v<T, float>) {
    // cppcheck-suppress suspiciousFloatingPointCast
    return static_cast<float>(val);
  } else if constexpr (std::is_same_v<T, double>) {
    return val;
  } else if constexpr (std::is_same_v<T, int>) {
    return static_cast<int>(val);
  } else {
    static_assert(
      std::is_same_v<T, bool> || std::is_same_v<T, float> || std::is_same_v<T, double> || std::is_same_v<T, int>,
      "getCfg<T>() only supports T = bool, float, double or int");
  }
}

int binarySearchVector(const int64_t Key, const std::vector<int64_t>& List, int low, int high) // v2 : list by reference
{
  while (low <= high) {
    int mid = low + (high - low) / 2;
    if (Key == List[mid]) {
      return mid;
    }

    if (Key > List[mid]) {
      low = mid + 1; // If Key is greater, ignore left  half, update the low
    } else {
      high = mid - 1;
    } // If Key is smaller, ignore right half, update the high
  }
  return -1; // Element is not present
}

// do a fastest search in an sorted array. ==> Binary Search in an array.
template <typename T>
bool checkTrackInList(const T& track, const std::vector<int64_t>& ParticleList) // v2 : list by reference
{
  return static_cast<bool>(binarySearchVector(track.globalIndex(), ParticleList, 0, ParticleList.size() - 1) != -1);
}

enum CollisionRejectionTag {
  kCollAccepted = 0,

  // Basic event selection
  kCollRejectSel8,
  kCollRejectVertexZ,
  kCollRejectTriggerTVX,

  // TF / ROF borders
  kCollRejectITSROFrameBorder,
  kCollRejectTimeFrameBorder,

  // Vertex quality
  kCollRejectVertexITSTPC,
  kCollRejectGoodZvtxFT0vsPV,
  kCollRejectVertexTOFmatched,
  kCollRejectVertexTRDmatched,
  kCollRejectGoodITSLayersAll,

  // Pileup / neighbouring collisions
  kCollRejectSameBunchPileup,
  kCollRejectCollInTimeRangeStandard,
  kCollRejectCollInTimeRangeStrict,
  kCollRejectCollInTimeRangeNarrow,
  kCollRejectCollInRofStandard,
  kCollRejectCollInRofStrict,
  kCollRejectHighMultCollInPrevRof,

  // INEL classes
  kCollRejectINELgt0,
  kCollRejectINELgt1,

  // Multiplicity
  kCollRejectMultiplicityLow,
  kCollRejectMultiplicityHigh,

  // Centrality
  kCollRejectCentralityLow,
  kCollRejectCentralityHigh,

  // Occupancy
  kCollRejectOccupancyLow,
  kCollRejectOccupancyHigh,

  // Interaction rate
  kCollRejectInteractionRateLow,
  kCollRejectInteractionRateHigh,

  // RCT
  kCollRejectRCT,
  kCollRejectCorrelationRCT,

  kNCollisionRejectionTags
};

enum TrackRejectionTag {
  kTrackAccepted = 0,

  // Kinematics
  kTrackRejectPt,
  kTrackRejectEta,
  kTrackRejectMomentum,
  kTrackRejectCharge,

  // Standard O2 track selections
  kTrackRejectQualityTrack,
  kTrackRejectQualityTrackITS,
  kTrackRejectQualityTrackTPC,
  kTrackRejectPrimaryTrack,
  kTrackRejectInAcceptanceTrack,
  kTrackRejectGlobalTrack,
  kTrackRejectGlobalTrackWoDCA,
  kTrackRejectGlobalTrackWoPtEta,
  kTrackRejectGlobalTrackWoTPCCluster,
  kTrackRejectGlobalTrackWoDCATPCCluster,
  kTrackRejectPVContributor,

  // Detector presence
  kTrackRejectITS,
  kTrackRejectTPC,
  kTrackRejectTOF,
  kTrackRejectTRD,

  // TPC quality
  kTrackRejectTPCNClsFound,
  kTrackRejectTPCCrossedRows,
  kTrackRejectTPCCrossedRowsOverFindable,
  kTrackRejectTPCFoundOverFindable,
  kTrackRejectTPCFractionShared,
  kTrackRejectTPCChi2,

  // ITS quality
  kTrackRejectITSNCls,
  kTrackRejectITSNClsInnerBarrel,
  kTrackRejectITSChi2,

  // DCA
  kTrackRejectFixedDCAxy,
  kTrackRejectFixedDCAz,
  kTrackRejectPtDependentDCAxy,
  kTrackRejectPtDependentDCAz,

  kNTrackRejectionTags
};

enum QAStageEnum {
  kBeforeSelection = 0,
  kAfterSelFunction,
  kAfterSelection,
  kNQAStages
};

static constexpr std::array<std::string_view, kNQAStages> QAStageDire = {
  "BeforeSelection/",
  "AfterSelFunction/",
  "AfterSelection/"};

enum EventQABin {
  kEventAllFiltered = 0,
  kEventPassedSelCollision,
  kEventPassedMinFilteredTracks,
  kEventPassedMinSelectedTracks,
  kNEventQABins
};

enum TrackQABin {
  kTrackAllFiltered = 0,
  kTrackPassedSelectionTrack,
  kNTrackQABins
};

enum DataFrameQABin {
  kDFWithZeroFilteredColls = 0,
  kNDataFrameQABins
};

enum EventTypeEnum {
  kSameEvent = 0,
  kMixedEvent,
  kNEventTypes
};

static constexpr std::array<std::string_view, kNEventTypes> EventTypeDire = {
  "SE/",
  "ME/"};

enum PairTypeEnum {
  kHH = 0,
  kHId,
  kHPhi,
  kHKStar,
  kHKStarBar,
  kHLambda,
  kHLambdaBar,
  kPhiPhi,
  kNPairTypes
};

static constexpr std::array<std::string_view, kNPairTypes> PairTypeDire = {
  "hh/",
  "hIdentified/",
  "hPhi/",
  "hKStar/",
  "hKStarBar/",
  "hLambda/",
  "hLambdaBar/",
  "PhiPhi/"};

enum AnalysisTypeEnum {
  kQA = 0,
  kCorr,
  kNAnalysisTypes
};

static constexpr std::array<std::string_view, kNAnalysisTypes> AnalysisTypeDire = {
  "QA/",
  "Corr/"};

enum CorrRoleEnum {
  kTrigger = 0,
  kAssocLowPt,
  kAssocHighPt,
  kNCorrRoles
};

static constexpr std::array<std::string_view, kNCorrRoles> CorrRoleDire = {
  "Trigger/",
  "AssocLowPt/",
  "AssocHighPt/"};

enum PrimVtxMassRegionTag : uint8_t {
  kMassOutside = 0,
  kMassPeak = 1,
  kMassLSB = 2,
  kMassRSB = 3,
  kMassNone = 4,
  kNMassRegions = 5
};

static constexpr std::array<std::string_view, kNMassRegions> MassRegionDire = {
  "Outside/",
  "Peak/",
  "LSB/",
  "RSB/",
  ""};

enum DauTypeEnum {
  kPosDau = 0,
  kNegDau,
  kNDauTypes
};

static constexpr std::array<std::string_view, kNDauTypes> DauTypeDire = {
  "PosDau/",
  "NegDau/"};

// template <typename T>
// uint8_t getMassRegionTag(T mass, T lsbLow, T lsbUp, T peakLow, T peakUp, T rsbLow, T rsbUp)
uint8_t getMassRegionTag(float mass, float lsbLow, float lsbUp, float peakLow, float peakUp, float rsbLow, float rsbUp)
{
  if (lsbLow < mass && mass < lsbUp) {
    return kMassLSB;
  }
  if (peakLow < mass && mass < peakUp) {
    return kMassPeak;
  }
  if (rsbLow < mass && mass < rsbUp) {
    return kMassRSB;
  }
  return kMassOutside;
}

enum CorrCountTypeEnum {
  kCountHH = 0,
  kCountHPi,
  kCountHKa,
  kCountHPr,
  kCountHPhi,
  kCountHKStar,
  kCountHKStarBar,
  kCountHLambda,
  kCountHLambdaBar,
  kNCorrCountTypes
};

enum CorrCountPtEnum {
  kCountLowPt = 0,
  kCountHighPt,
  kNCorrCountPt
};

static constexpr int NCorrCountChannels = static_cast<int>(kNCorrCountTypes) * static_cast<int>(kNCorrCountPt);
static constexpr std::array<std::string_view, kNCorrCountTypes> CorrCountTypeName = {
  "hh",
  "hPi",
  "hKa",
  "hPr",
  "hPhi",
  "hKStar",
  "hKStarBar",
  "hLambda",
  "hLambdaBar"};

using CorrCountArray = std::array<std::array<uint64_t, kNCorrCountPt>, kNCorrCountTypes>;
enum ResoCorrCountTypeEnum {
  kResoCountHPhi = 0,
  kResoCountHKStar,
  kResoCountHKStarBar,
  kResoCountHLambda,
  kResoCountHLambdaBar,
  kNResoCorrCountTypes
};

using ResoRegionCountArray = std::array<std::array<std::array<uint64_t, kNMassRegions>, kNCorrCountPt>, kNResoCorrCountTypes>;

struct MixingBinStatus {
  uint64_t nCollisions = 0;         // nCollisions = total selected collisions entering this mixing bin in this dataframe
  CorrCountArray nCorrCollisions{}; // nCorrCollisions[hPhi][LowPt] = how many COLLISIONS in this mixing bin had >=1 hPhi-low correlation
};

enum PidEnum {
  kPi = 0,
  kKa,
  kPr,
  kEl,
  kMu,
  kDe,
  kNPid
};

static constexpr std::array<std::string_view, kNPid> PidTypeDire = {
  "Pi/",
  "Ka/",
  "Pr/",
  "El/",
  "Mu/",
  "De/"};

enum CutSettingEnum {
  kThrPforTOF = 0,
  kIdCutTypeLowP,
  kNSigmaTPCLowP,
  kNSigmaTOFLowP,
  kNSigmaRadLowP,
  kIdCutTypeHighP,
  kNSigmaTPCHighP,
  kNSigmaTOFHighP,
  kNSigmaRadHighP,
  kDoVetoOthers,
  kDoRelativeTPCcheck,
  kDoRelativeTOFcheck,
  kDoRelativeTPCTOFcheck,
  kNCutSettings
};

enum VetoSettingEnum {
  kDoVetoTPC = 0,
  kDoVetoTOF,
  kVetoTPC,
  kVetoTOF,
  kNVetoSettings
};

enum IdentificationType {
  kTPCidentified = 0,
  kTPCTOFidentified,
  kUnidentified,
  kNIdentificationTypes
};

enum TpcTofCutType {
  kCircularCut = 0,
  kRectangularCut,
  kEllipsoidalCut,
  kNCutTypes
};

enum MixingEstimatorEnum {
  kMixCentFT0C = 0,
  kMixCentFT0M,
  kMixCentFT0A,
  kMixCentFV0A,
  kNMixingEstimators
};

enum CorrPresenceMask : uint64_t {
  kMaskHHLowPt = 1ULL << 0,
  kMaskHHHighPt = 1ULL << 1,
  kMaskHPiLowPt = 1ULL << 2,
  kMaskHPiHighPt = 1ULL << 3,
  kMaskHKaLowPt = 1ULL << 4,
  kMaskHKaHighPt = 1ULL << 5,
  kMaskHPrLowPt = 1ULL << 6,
  kMaskHPrHighPt = 1ULL << 7,

  kMaskHPhiLowPt = 1ULL << 8,
  kMaskHPhiHighPt = 1ULL << 9,
  kMaskHKStarLowPt = 1ULL << 10,
  kMaskHKStarHighPt = 1ULL << 11,
  kMaskHKStarBarLowPt = 1ULL << 12,
  kMaskHKStarBarHighPt = 1ULL << 13,
  kMaskHLambdaLowPt = 1ULL << 14,
  kMaskHLambdaHighPt = 1ULL << 15,
  kMaskHLambdaBarLowPt = 1ULL << 16,
  kMaskHLambdaBarHighPt = 1ULL << 17,

  kMaskHPhiLowPtPeak = 1ULL << 18,
  kMaskHPhiLowPtLSB = 1ULL << 19,
  kMaskHPhiLowPtRSB = 1ULL << 20,
  kMaskHPhiHighPtPeak = 1ULL << 21,
  kMaskHPhiHighPtLSB = 1ULL << 22,
  kMaskHPhiHighPtRSB = 1ULL << 23,

  kMaskHKStarLowPtPeak = 1ULL << 24,
  kMaskHKStarLowPtLSB = 1ULL << 25,
  kMaskHKStarLowPtRSB = 1ULL << 26,
  kMaskHKStarHighPtPeak = 1ULL << 27,
  kMaskHKStarHighPtLSB = 1ULL << 28,
  kMaskHKStarHighPtRSB = 1ULL << 29,

  kMaskHKStarBarLowPtPeak = 1ULL << 30,
  kMaskHKStarBarLowPtLSB = 1ULL << 31,
  kMaskHKStarBarLowPtRSB = 1ULL << 32,
  kMaskHKStarBarHighPtPeak = 1ULL << 33,
  kMaskHKStarBarHighPtLSB = 1ULL << 34,
  kMaskHKStarBarHighPtRSB = 1ULL << 35,

  kMaskHLambdaLowPtPeak = 1ULL << 36,
  kMaskHLambdaLowPtLSB = 1ULL << 37,
  kMaskHLambdaLowPtRSB = 1ULL << 38,
  kMaskHLambdaHighPtPeak = 1ULL << 39,
  kMaskHLambdaHighPtLSB = 1ULL << 40,
  kMaskHLambdaHighPtRSB = 1ULL << 41,

  kMaskHLambdaBarLowPtPeak = 1ULL << 42,
  kMaskHLambdaBarLowPtLSB = 1ULL << 43,
  kMaskHLambdaBarLowPtRSB = 1ULL << 44,
  kMaskHLambdaBarHighPtPeak = 1ULL << 45,
  kMaskHLambdaBarHighPtLSB = 1ULL << 46,
  kMaskHLambdaBarHighPtRSB = 1ULL << 47,

  kMaskAll = (1ULL << 48) - 1ULL
};

constexpr uint64_t getCorrMaskBit(int corrType, int ptType)
{
  return 1ULL << (corrType * kNCorrCountPt + ptType);
}

constexpr uint64_t getResoRegionMaskBit(int resoType, int ptType, int massRegion)
{
  if (massRegion < kMassPeak || massRegion > kMassRSB) {
    return 0ULL;
  }
  return 1ULL << (18 + resoType * 6 + ptType * 3 + (massRegion - kMassPeak));
}

// Bits 18-47 (kMaskHPhiLowPtPeak ... kMaskHLambdaBarHighPtRSB): resonance (species x pT x mass-window) pools
static constexpr uint64_t MaskResoRegions = static_cast<uint64_t>(kMaskAll) & ~((1ULL << 18) - 1ULL);

static constexpr std::array<std::array<double, kNCutSettings>, kNPid> DefaultPIDcheckValues = {{
  {0.7, 0, 3.0, 3.0, 9.0, 0, 3.0, 3.0, 9.0, 0, 1, 0, 1}, // Pi
  {0.5, 0, 3.0, 3.0, 9.0, 0, 3.0, 3.0, 9.0, 0, 1, 0, 1}, // Ka
  {0.8, 0, 3.0, 3.0, 9.0, 0, 3.0, 3.0, 9.0, 0, 1, 0, 1}, // Pr
  {0.4, 0, 3.0, 3.0, 9.0, 0, 3.0, 3.0, 9.0, 0, 1, 0, 1}, // El
  {0.4, 0, 3.0, 3.0, 9.0, 0, 3.0, 3.0, 9.0, 0, 1, 0, 1}, // Mu
  {0.4, 0, 3.0, 3.0, 9.0, 0, 3.0, 3.0, 9.0, 0, 1, 0, 1}  // De
}};

static constexpr std::array<std::array<double, kNVetoSettings>, kNPid> DefaultPidVetoValues = {{
  {1, 1, 3.0, 3.0}, // Pi
  {1, 1, 3.0, 3.0}, // Ka
  {1, 1, 3.0, 3.0}, // Pr
  {0, 0, 3.0, 3.0}, // El
  {0, 0, 3.0, 3.0}, // Mu
  {0, 0, 3.0, 3.0}  // De
}};

//______________________________________________________________________________________________________________
// Common functions independent of configurables used across the tasks

template <int pidMode, typename T>
bool relativeIdOthersTPC(const T& track)
{
  float distSpecies = 1.e6f;

  if constexpr (pidMode == kPi) {
    distSpecies = std::fabs(track.tpcNSigmaPi());
  } else if constexpr (pidMode == kKa) {
    distSpecies = std::fabs(track.tpcNSigmaKa());
  } else if constexpr (pidMode == kPr) {
    distSpecies = std::fabs(track.tpcNSigmaPr());
  } else {
    // FIX ME! : Relative check will be updated in future, currently implemented for pion, kaon and proton only
    // else if constexpr (pidMode == kEl) {
    //   distSpecies = std::fabs(track.tpcNSigmaEl());
    // } else if constexpr (pidMode == kMu) {
    //   distSpecies = std::fabs(track.tpcNSigmaMu());
    // } else if constexpr (pidMode == kDe) {
    //   distSpecies = std::fabs(track.tpcNSigmaDe());
    // }

    return true;
  }

  if constexpr (pidMode != kPi) {
    if (std::fabs(track.tpcNSigmaPi()) <= distSpecies) {
      return false;
    }
  }

  if constexpr (pidMode != kKa) {
    if (std::fabs(track.tpcNSigmaKa()) <= distSpecies) {
      return false;
    }
  }

  if constexpr (pidMode != kPr) {
    if (std::fabs(track.tpcNSigmaPr()) <= distSpecies) {
      return false;
    }
  }

  // FIX ME! : Relative check will be updated in future, currently implemented for pion, kaon and proton only
  // if constexpr (pidMode != kEl)
  //   if (std::fabs(track.tpcNSigmaEl()) <= distSpecies)
  //     return false;

  // if constexpr (pidMode != kMu)
  //   if (std::fabs(track.tpcNSigmaMu()) <= distSpecies)
  //     return false;

  // if constexpr (pidMode != kDe)
  //   if (std::fabs(track.tpcNSigmaDe()) <= distSpecies)
  //     return false;

  return true;
}

template <int pidMode, typename T>
bool relativeIdOthersTOF(const T& track)
{
  float distSpecies = 1.e6f;

  if constexpr (pidMode == kPi) {
    distSpecies = std::fabs(track.tofNSigmaPi());
  } else if constexpr (pidMode == kKa) {
    distSpecies = std::fabs(track.tofNSigmaKa());
  } else if constexpr (pidMode == kPr) {
    distSpecies = std::fabs(track.tofNSigmaPr());

    // FIX ME! : Relative check will be updated in future, currently implemented for pion, kaon and proton only
    // } else if constexpr (pidMode == kEl) {
    //   distSpecies = std::fabs(track.tofNSigmaEl());
    // } else if constexpr (pidMode == kMu) {
    //   distSpecies = std::fabs(track.tofNSigmaMu());
    // } else if constexpr (pidMode == kDe) {
    //   distSpecies = std::fabs(track.tofNSigmaDe());
  } else {
    return true;
  }

  if constexpr (pidMode != kPi) {
    if (std::fabs(track.tofNSigmaPi()) <= distSpecies) {
      return false;
    }
  }

  if constexpr (pidMode != kKa) {
    if (std::fabs(track.tofNSigmaKa()) <= distSpecies) {
      return false;
    }
  }

  if constexpr (pidMode != kPr) {
    if (std::fabs(track.tofNSigmaPr()) <= distSpecies) {
      return false;
    }
  }
  // FIX ME! : Relative check will be updated in future, currently implemented for pion, kaon and proton only
  // if constexpr (pidMode != kEl)
  //   if (std::fabs(track.tofNSigmaEl()) <= distSpecies)
  //     return false;

  // if constexpr (pidMode != kMu)
  //   if (std::fabs(track.tofNSigmaMu()) <= distSpecies)
  //     return false;

  // if constexpr (pidMode != kDe)
  //   if (std::fabs(track.tofNSigmaDe()) <= distSpecies)
  //     return false;

  return true;
}

template <int pidMode, typename T>
bool relativeIdOthersTPCTOF(const T& track)
{
  float distSpeciesSq = 1.e6f;

  if constexpr (pidMode == kPi) {
    distSpeciesSq = track.tpcNSigmaPi() * track.tpcNSigmaPi() + track.tofNSigmaPi() * track.tofNSigmaPi();
  } else if constexpr (pidMode == kKa) {
    distSpeciesSq = track.tpcNSigmaKa() * track.tpcNSigmaKa() + track.tofNSigmaKa() * track.tofNSigmaKa();
  } else if constexpr (pidMode == kPr) {
    distSpeciesSq = track.tpcNSigmaPr() * track.tpcNSigmaPr() + track.tofNSigmaPr() * track.tofNSigmaPr();
    // FIX ME! : Relative check will be updated in future, currently implemented for pion, kaon and proton only
    // } else if constexpr (pidMode == kEl) {
    //   distSpeciesSq = track.tpcNSigmaEl() * track.tpcNSigmaEl() + track.tofNSigmaEl() * track.tofNSigmaEl();
    // } else if constexpr (pidMode == kMu) {
    //   distSpeciesSq = track.tpcNSigmaMu() * track.tpcNSigmaMu() + track.tofNSigmaMu() * track.tofNSigmaMu();
    // } else if constexpr (pidMode == kDe) {
    //   distSpeciesSq = track.tpcNSigmaDe() * track.tpcNSigmaDe() + track.tofNSigmaDe() * track.tofNSigmaDe();
  } else {
    return true;
  }

  if constexpr (pidMode != kPi) {
    const float distSq = track.tpcNSigmaPi() * track.tpcNSigmaPi() + track.tofNSigmaPi() * track.tofNSigmaPi();
    if (distSq <= distSpeciesSq) {
      return false;
    }
  }

  if constexpr (pidMode != kKa) {
    const float distSq = track.tpcNSigmaKa() * track.tpcNSigmaKa() + track.tofNSigmaKa() * track.tofNSigmaKa();
    if (distSq <= distSpeciesSq) {
      return false;
    }
  }

  if constexpr (pidMode != kPr) {
    const float distSq = track.tpcNSigmaPr() * track.tpcNSigmaPr() + track.tofNSigmaPr() * track.tofNSigmaPr();
    if (distSq <= distSpeciesSq) {
      return false;
    }
  }
  // FIX ME! : Relative check will be updated in future, currently implemented for pion, kaon and proton only
  // if constexpr (pidMode != kEl) {
  //   const float distSq = track.tpcNSigmaEl() * track.tpcNSigmaEl() + track.tofNSigmaEl() * track.tofNSigmaEl();
  //   if (distSq <= distSpeciesSq)
  //     return false;
  // }

  // if constexpr (pidMode != kMu) {
  //   const float distSq = track.tpcNSigmaMu() * track.tpcNSigmaMu() + track.tofNSigmaMu() * track.tofNSigmaMu();
  //   if (distSq <= distSpeciesSq)
  //     return false;
  // }

  // if constexpr (pidMode != kDe) {
  //   const float distSq = track.tpcNSigmaDe() * track.tpcNSigmaDe() + track.tofNSigmaDe() * track.tofNSigmaDe();
  //   if (distSq <= distSpeciesSq)
  //     return false;
  // }

  return true;
}

template <int pidMode, typename T>
bool selIdRectangularCut(const T& track, const float& nSigmaTPC, const float& nSigmaTOF)
{
  if constexpr (pidMode == kPi) {
    return std::fabs(track.tpcNSigmaPi()) < nSigmaTPC && std::fabs(track.tofNSigmaPi()) < nSigmaTOF;
  } else if constexpr (pidMode == kKa) {
    return std::fabs(track.tpcNSigmaKa()) < nSigmaTPC && std::fabs(track.tofNSigmaKa()) < nSigmaTOF;
  } else if constexpr (pidMode == kPr) {
    return std::fabs(track.tpcNSigmaPr()) < nSigmaTPC && std::fabs(track.tofNSigmaPr()) < nSigmaTOF;
  } else if constexpr (pidMode == kEl) {
    return std::fabs(track.tpcNSigmaEl()) < nSigmaTPC && std::fabs(track.tofNSigmaEl()) < nSigmaTOF;
  } else if constexpr (pidMode == kMu) {
    return std::fabs(track.tpcNSigmaMu()) < nSigmaTPC && std::fabs(track.tofNSigmaMu()) < nSigmaTOF;
  } else if constexpr (pidMode == kDe) {
    return std::fabs(track.tpcNSigmaDe()) < nSigmaTPC && std::fabs(track.tofNSigmaDe()) < nSigmaTOF;
  } else {
    return false;
  }
}

template <int pidMode, typename T>
bool selIdEllipsoidalCut(const T& track, const float& nSigmaTPC, const float& nSigmaTOF)
{
  if constexpr (pidMode == kPi) {
    float tpc = track.tpcNSigmaPi() / nSigmaTPC;
    float tof = track.tofNSigmaPi() / nSigmaTOF;
    return (tpc * tpc + tof * tof) < 1.0f;
  } else if constexpr (pidMode == kKa) {
    float tpc = track.tpcNSigmaKa() / nSigmaTPC;
    float tof = track.tofNSigmaKa() / nSigmaTOF;
    return (tpc * tpc + tof * tof) < 1.0f;
  } else if constexpr (pidMode == kPr) {
    float tpc = track.tpcNSigmaPr() / nSigmaTPC;
    float tof = track.tofNSigmaPr() / nSigmaTOF;
    return (tpc * tpc + tof * tof) < 1.0f;
  } else if constexpr (pidMode == kEl) {
    float tpc = track.tpcNSigmaEl() / nSigmaTPC;
    float tof = track.tofNSigmaEl() / nSigmaTOF;
    return (tpc * tpc + tof * tof) < 1.0f;
  } else if constexpr (pidMode == kMu) {
    float tpc = track.tpcNSigmaMu() / nSigmaTPC;
    float tof = track.tofNSigmaMu() / nSigmaTOF;
    return (tpc * tpc + tof * tof) < 1.0f;
  } else if constexpr (pidMode == kDe) {
    float tpc = track.tpcNSigmaDe() / nSigmaTPC;
    float tof = track.tofNSigmaDe() / nSigmaTOF;
    return (tpc * tpc + tof * tof) < 1.0f;
  } else {
    return false;
  }
}

// Circular cut : (tpcNSigma^2 + tofNSigma^2) < nSigmaSquaredRad ; for 3 sigma, nSigmaSquaredRad = 9.0
template <int pidMode, typename T>
constexpr bool selIdCircularCut(const T& track, const float& nSigmaSquaredRad)
{
  float tpc = 0.f;
  float tof = 0.f;

  if constexpr (pidMode == kPi) {
    tpc = track.tpcNSigmaPi();
    tof = track.tofNSigmaPi();
  } else if constexpr (pidMode == kKa) {
    tpc = track.tpcNSigmaKa();
    tof = track.tofNSigmaKa();
  } else if constexpr (pidMode == kPr) {
    tpc = track.tpcNSigmaPr();
    tof = track.tofNSigmaPr();
  } else if constexpr (pidMode == kEl) {
    tpc = track.tpcNSigmaEl();
    tof = track.tofNSigmaEl();
  } else if constexpr (pidMode == kMu) {
    tpc = track.tpcNSigmaMu();
    tof = track.tofNSigmaMu();
  } else if constexpr (pidMode == kDe) {
    tpc = track.tpcNSigmaDe();
    tof = track.tofNSigmaDe();
  } else {
    return false; // unknown pidMode
  }
  return (tpc * tpc + tof * tof) < nSigmaSquaredRad;
}

template <int pidMode, typename T>
bool applyVetoOthersTPC(const T& track, const std::array<bool, kNPid>& doVetoTPC, const std::array<float, kNPid>& vetoTPC)
{
  if constexpr (pidMode != kPi) {
    if (doVetoTPC[kPi] && std::fabs(track.tpcNSigmaPi()) < vetoTPC[kPi]) {
      return false;
    }
  }

  if constexpr (pidMode != kKa) {
    if (doVetoTPC[kKa] && std::fabs(track.tpcNSigmaKa()) < vetoTPC[kKa]) {
      return false;
    }
  }

  if constexpr (pidMode != kPr) {
    if (doVetoTPC[kPr] && std::fabs(track.tpcNSigmaPr()) < vetoTPC[kPr]) {
      return false;
    }
  }

  if constexpr (pidMode != kEl) {
    if (doVetoTPC[kEl] && std::fabs(track.tpcNSigmaEl()) < vetoTPC[kEl]) {
      return false;
    }
  }

  if constexpr (pidMode != kMu) {
    if (doVetoTPC[kMu] && std::fabs(track.tpcNSigmaMu()) < vetoTPC[kMu]) {
      return false;
    }
  }

  if constexpr (pidMode != kDe) {
    if (doVetoTPC[kDe] && std::fabs(track.tpcNSigmaDe()) < vetoTPC[kDe]) {
      return false;
    }
  }

  return true;
}

template <int pidMode, typename T>
bool applyVetoOthersTOF(const T& track, const std::array<bool, kNPid>& doVetoTOF, const std::array<float, kNPid>& vetoTOF)
{
  if constexpr (pidMode != kPi) {
    if (doVetoTOF[kPi] && std::fabs(track.tofNSigmaPi()) < vetoTOF[kPi]) {
      return false;
    }
  }

  if constexpr (pidMode != kKa) {
    if (doVetoTOF[kKa] && std::fabs(track.tofNSigmaKa()) < vetoTOF[kKa]) {
      return false;
    }
  }

  if constexpr (pidMode != kPr) {
    if (doVetoTOF[kPr] && std::fabs(track.tofNSigmaPr()) < vetoTOF[kPr]) {
      return false;
    }
  }

  if constexpr (pidMode != kEl) {
    if (doVetoTOF[kEl] && std::fabs(track.tofNSigmaEl()) < vetoTOF[kEl]) {
      return false;
    }
  }

  if constexpr (pidMode != kMu) {
    if (doVetoTOF[kMu] && std::fabs(track.tofNSigmaMu()) < vetoTOF[kMu]) {
      return false;
    }
  }

  if constexpr (pidMode != kDe) {
    if (doVetoTOF[kDe] && std::fabs(track.tofNSigmaDe()) < vetoTOF[kDe]) {
      return false;
    }
  }

  return true;
}

// Common Helper functions are over.
//__________________________________________________________________________________________________________________________

//__________________________________________________________________________________________________________________________
// Common Histograms booking functions - to be used inside
template <typename H>
void addBasicTrackQAHistos(H& histReg, const std::string& basePath, const AxisSpec& axisP, const AxisSpec& axisPt, const AxisSpec& axisEta, const AxisSpec& axisPhi, const AxisSpec& axisDcaXY, const AxisSpec& axisDcaZ, const AxisSpec& axisSign)
{
  histReg.add((basePath + "P").c_str(), "p;p (GeV/c);Counts", kTH1F, {axisP});
  histReg.add((basePath + "Pt").c_str(), "p_{T};p_{T} (GeV/c);Counts", kTH1F, {axisPt});
  histReg.add((basePath + "Eta").c_str(), "#eta;#eta;Counts", kTH1F, {axisEta});
  histReg.add((basePath + "Phi").c_str(), "#varphi;#varphi;Counts", kTH1F, {axisPhi});
  histReg.add((basePath + "DCAxy").c_str(), "DCA_{xy};DCA_{xy} (cm);Counts", kTH1F, {axisDcaXY});
  histReg.add((basePath + "DCAz").c_str(), "DCA_{z};DCA_{z} (cm);Counts", kTH1F, {axisDcaZ});
  histReg.add((basePath + "Sign").c_str(), "Sign;Sign;Counts", kTH1F, {axisSign});
  histReg.add((basePath + "DCAxyVsPt").c_str(), "DCA_{xy} vs p_{T};p_{T} (GeV/c);DCA_{xy} (cm)", kTH2F, {axisPt, axisDcaXY});
  histReg.add((basePath + "DCAzVsPt").c_str(), "DCA_{z} vs p_{T};p_{T} (GeV/c);DCA_{z} (cm)", kTH2F, {axisPt, axisDcaZ});
}

template <typename H>
void addIdentifiedQAHistos(H& histReg, const std::string& basePath, const AxisSpec& axisP, const AxisSpec& axisRapidity, const AxisSpec& axisTPCSignal, const AxisSpec& axisTOFBeta, const AxisSpec& axisTPCNSigma, const AxisSpec& axisTOFNSigma, const AxisSpec& axisIdMethod)
{
  histReg.add((basePath + "Rapidity").c_str(), "Rapidity;y;Counts", HistType::kTH1F, {axisRapidity});
  histReg.add((basePath + "TPCSignalVsP").c_str(), "TPC signal vs p;p (GeV/c);TPC signal", HistType::kTH2F, {axisP, axisTPCSignal});
  histReg.add((basePath + "TOFBetaVsP").c_str(), "TOF #beta vs p;p (GeV/c);#beta", HistType::kTH2F, {axisP, axisTOFBeta});
  histReg.add((basePath + "TPCNSigmaVsP").c_str(), "TPC n#sigma vs p;p (GeV/c);n#sigma_{TPC}", HistType::kTH2F, {axisP, axisTPCNSigma});
  histReg.add((basePath + "TOFNSigmaVsP").c_str(), "TOF n#sigma vs p;p (GeV/c);n#sigma_{TOF}", HistType::kTH2F, {axisP, axisTOFNSigma});
  histReg.add((basePath + "TPCNSigmaVsTOFNSigma").c_str(), "TPC n#sigma vs TOF n#sigma;n#sigma_{TPC};n#sigma_{TOF}", HistType::kTH2F, {axisTPCNSigma, axisTOFNSigma});
  histReg.add((basePath + "IdMethodVsP").c_str(), "Identification method vs p;p (GeV/c);ID method", HistType::kTH2F, {axisP, axisIdMethod});
}

// Common Histograms booking functions are over
//__________________________________________________________________________________________________________________________

struct HParticleCorrelation1ResonanceProducer {
  Produces<o2::aod::ResonanceCndts> resonanceCndts;
  HistogramRegistry resonanceQAPlots{"resonanceQAPlots", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  static constexpr float ResonancePidNSigmaMax = 3.5f;

  struct : ConfigurableGroup {
    Configurable<bool> c0001printDebugMessages{"c0001printDebugMessages", false, "Print debug messages"};
  } cfgDebug;

  struct : ConfigurableGroup {
    Configurable<bool> c0002requireSel8{"c0002requireSel8", true, "Require sel8 event selection"};
    Configurable<float> c0003cutZvertex{"c0003cutZvertex", 10.0f, "Maximum |z_{vtx}| (cm)"};
  } cfgEvent;

  struct : ConfigurableGroup {
    Configurable<float> c0004ptMin{"c0004ptMin", 0.15f, "Minimum track pT (GeV/c)"};
    Configurable<float> c0005ptMax{"c0005ptMax", 100.0f, "Maximum track pT (GeV/c)"};
    Configurable<float> c0006etaMax{"c0006etaMax", 0.8f, "Maximum |eta|"};

    Configurable<bool> c0007useFixedDCAxy{"c0007useFixedDCAxy", false, "Apply fixed DCAxy cut"};
    Configurable<float> c0008dcaXYMax{"c0008dcaXYMax", 0.1f, "Maximum |DCAxy| (cm)"};

    Configurable<bool> c0009useFixedDCAz{"c0009useFixedDCAz", false, "Apply fixed DCAz cut"};
    Configurable<float> c0010dcaZMax{"c0010dcaZMax", 0.2f, "Maximum |DCAz| (cm)"};
  } cfgTrackCuts;

  struct : ConfigurableGroup {
    Configurable<bool> c0011doPhi1020{"c0011doPhi1020", true, "Build Phi(1020) candidates"};
    Configurable<bool> c0012doKStar892{"c0012doKStar892", true, "Build K*(892)0 candidates"};
    Configurable<bool> c0013doKStar892Bar{"c0013doKStar892Bar", true, "Build anti-K*(892)0 candidates"};
    Configurable<bool> c0014doLambda1520{"c0014doLambda1520", true, "Build Lambda(1520) candidates"};
    Configurable<bool> c0015doLambda1520Bar{"c0015doLambda1520Bar", true, "Build anti-Lambda(1520) candidates"};
  } cfgResonances;

  struct : ConfigurableGroup {
    Configurable<float> c0016phi1020PeakLow{"c0016phi1020PeakLow", 1.013f, "Phi(1020) peak lower mass"};
    Configurable<float> c0017phi1020PeakUp{"c0017phi1020PeakUp", 1.026f, "Phi(1020) peak upper mass"};
    Configurable<float> c0018phi1020LSBLow{"c0018phi1020LSBLow", 0.995f, "Phi(1020) LSB lower mass"};
    Configurable<float> c0019phi1020LSBUp{"c0019phi1020LSBUp", 1.005f, "Phi(1020) LSB upper mass"};
    Configurable<float> c0020phi1020RSBLow{"c0020phi1020RSBLow", 1.040f, "Phi(1020) RSB lower mass"};
    Configurable<float> c0021phi1020RSBUp{"c0021phi1020RSBUp", 1.060f, "Phi(1020) RSB upper mass"};
  } cfgPhi1020Mass;

  struct : ConfigurableGroup {
    Configurable<float> c0022kstar892PeakLow{"c0022kstar892PeakLow", 0.846f, "K*(892)0 peak lower mass"};
    Configurable<float> c0023kstar892PeakUp{"c0023kstar892PeakUp", 0.946f, "K*(892)0 peak upper mass"};
    Configurable<float> c0024kstar892LSBLow{"c0024kstar892LSBLow", 0.700f, "K*(892)0 LSB lower mass"};
    Configurable<float> c0025kstar892LSBUp{"c0025kstar892LSBUp", 0.780f, "K*(892)0 LSB upper mass"};
    Configurable<float> c0026kstar892RSBLow{"c0026kstar892RSBLow", 1.010f, "K*(892)0 RSB lower mass"};
    Configurable<float> c0027kstar892RSBUp{"c0027kstar892RSBUp", 1.090f, "K*(892)0 RSB upper mass"};
  } cfgKstar892Mass;

  struct : ConfigurableGroup {
    Configurable<float> c0028lambda1520PeakLow{"c0028lambda1520PeakLow", 1.500f, "Lambda(1520) peak lower mass"};
    Configurable<float> c0029lambda1520PeakUp{"c0029lambda1520PeakUp", 1.540f, "Lambda(1520) peak upper mass"};
    Configurable<float> c0030lambda1520LSBLow{"c0030lambda1520LSBLow", 1.440f, "Lambda(1520) LSB lower mass"};
    Configurable<float> c0031lambda1520LSBUp{"c0031lambda1520LSBUp", 1.480f, "Lambda(1520) LSB upper mass"};
    Configurable<float> c0032lambda1520RSBLow{"c0032lambda1520RSBLow", 1.560f, "Lambda(1520) RSB lower mass"};
    Configurable<float> c0033lambda1520RSBUp{"c0033lambda1520RSBUp", 1.600f, "Lambda(1520) RSB upper mass"};
  } cfgLambda1520Mass;

  void init(InitContext const&)
  {
    const AxisSpec axisPhi1020Mass = {300, 0.98f, 1.08f, "m_{K^{+}K^{-}} (GeV/#it{c}^{2})"};
    const AxisSpec axisKStar892Mass = {300, 0.65f, 1.15f, "m_{K#pi} (GeV/#it{c}^{2})"};
    const AxisSpec axisLambda1520Mass = {300, 1.40f, 1.65f, "m_{pK} (GeV/#it{c}^{2})"};

    resonanceQAPlots.add("mPhi1020", "#phi(1020);m_{K^{+}K^{-}} (GeV/#it{c}^{2});Counts", kTH1F, {axisPhi1020Mass});
    resonanceQAPlots.add("mKStar892", "K^{*}(892)^{0};m_{K^{+}#pi^{-}} (GeV/#it{c}^{2});Counts", kTH1F, {axisKStar892Mass});
    resonanceQAPlots.add("mKStar892Bar", "#bar{K}^{*}(892)^{0};m_{#pi^{+}K^{-}} (GeV/#it{c}^{2});Counts", kTH1F, {axisKStar892Mass});
    resonanceQAPlots.add("mLambda1520", "#Lambda(1520);m_{pK^{-}} (GeV/#it{c}^{2});Counts", kTH1F, {axisLambda1520Mass});
    resonanceQAPlots.add("mLambda1520Bar", "#bar{#Lambda}(1520);m_{K^{+}#bar{p}} (GeV/#it{c}^{2});Counts", kTH1F, {axisLambda1520Mass});
  }

  template <typename T>
  bool selPion(const T& track)
  {
    if (std::abs(track.tpcNSigmaPi()) < ResonancePidNSigmaMax) {
      return true;
    }
    if (track.hasTOF()) {
      if (std::abs(track.tofNSigmaPi()) < ResonancePidNSigmaMax) {
        return true;
      }
    }
    return false;
  }

  template <typename T>
  bool selKaon(const T& track)
  {
    if (std::abs(track.tpcNSigmaKa()) < ResonancePidNSigmaMax) {
      return true;
    }
    if (track.hasTOF()) {
      if (std::abs(track.tofNSigmaKa()) < ResonancePidNSigmaMax) {
        return true;
      }
    }
    return false;
  }

  template <typename T>
  bool selProton(const T& track)
  {
    if (std::abs(track.tpcNSigmaPr()) < ResonancePidNSigmaMax) {
      return true;
    }
    if (track.hasTOF()) {
      if (std::abs(track.tofNSigmaPr()) < ResonancePidNSigmaMax) {
        return true;
      }
    }
    return false;
  }

  // Event Filter
  Filter eventFilter = (!cfgEvent.c0002requireSel8) || (o2::aod::evsel::sel8 == true);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < cfgEvent.c0003cutZvertex);

  // Track Filter
  Filter ptFilter = (o2::aod::track::pt > cfgTrackCuts.c0004ptMin) && (o2::aod::track::pt < cfgTrackCuts.c0005ptMax);
  Filter etaFilter = (nabs(o2::aod::track::eta) < cfgTrackCuts.c0006etaMax);
  Filter dcaFilter = ((!cfgTrackCuts.c0007useFixedDCAxy) || (nabs(o2::aod::track::dcaXY) < cfgTrackCuts.c0008dcaXYMax)) &&
                     ((!cfgTrackCuts.c0009useFixedDCAz) || (nabs(o2::aod::track::dcaZ) < cfgTrackCuts.c0010dcaZMax));
  using MyFilteredCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>;
  using MyFilteredTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::TOFSignal, aod::pidTOFbeta, aod::pidTOFmass, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullEl, aod::pidTPCFullDe, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullEl, aod::pidTOFFullDe>>;

  Preslice<MyFilteredTracks> tracksPerCollisionPreslice = o2::aod::track::collisionId;
  SliceCache cache;
  Partition<MyFilteredTracks> posTracks = aod::track::signed1Pt > 0.0f;
  Partition<MyFilteredTracks> negTracks = aod::track::signed1Pt < 0.0f;

  int dfNumber = 0;
  // void processNothing(aod::HPCorrColls const&)
  void processNothing(aod::Origins const& origins)
  {
    if (cfgDebug.c0001printDebugMessages) {
      LOG(info) << "DEBUG :: Process Nothing :: df_" << dfNumber << " :: origins = " << origins.size();
    }
    // Intentionally empty.
    // Keeps the task alive when running purely on derived data.
  }
  PROCESS_SWITCH(HParticleCorrelation1ResonanceProducer, processNothing, "Dummy process for derived-data analysis", true);

  void processSameEvent(MyFilteredCollisions const& collisions, MyFilteredTracks const& fullTracks, o2::aod::Origins const& /*Origins*/)
  {
    dfNumber++;
    if (cfgDebug.c0001printDebugMessages) {
      LOG(info) << "DEBUG :: df_" << dfNumber << " :: SE :: collisions = " << collisions.size() << " :: fullTracks = " << fullTracks.size();
    }

    float ePosPi = 0.0f;
    float ePosKa = 0.0f;
    float ePosPr = 0.0f;

    float eNegPi = 0.0f;
    float eNegKa = 0.0f;
    float eNegPr = 0.0f;

    float pt = 0.0f;
    float eta = 0.0f;
    float phi = 0.0f;
    float px = 0.0f;
    float py = 0.0f;
    float pz = 0.0f;
    float p = 0.0f;

    float mPhi1020 = -1.0f;
    float mKStar892 = -1.0f;
    float mKStar892Bar = -1.0f;
    float mLambda1520 = -1.0f;
    float mLambda1520Bar = -1.0f;

    uint8_t phi1020Tag = kMassOutside;
    uint8_t kStar892Tag = kMassOutside;
    uint8_t kStar892BarTag = kMassOutside;
    uint8_t lambda1520Tag = kMassOutside;
    uint8_t lambda1520BarTag = kMassOutside;

    bool posIsPi = false, posIsKa = false, posIsPr = false, negIsPi = false, negIsKa = false, negIsPr = false;
    bool fillTable = false;

    for (const auto& collision : collisions) {
      auto posTracksPerColl = posTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
      auto negTracksPerColl = negTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

      for (const auto& posTrack : posTracksPerColl) {
        posIsPi = selPion(posTrack);
        posIsKa = selKaon(posTrack);
        posIsPr = selProton(posTrack);
        if (!(posIsPi || posIsKa || posIsPr)) {
          continue;
        }
        // Precompute energies for fast speed
        ePosPi = RecoDecay::e(posTrack.px(), posTrack.py(), posTrack.pz(), MassPiPlus);
        ePosKa = RecoDecay::e(posTrack.px(), posTrack.py(), posTrack.pz(), MassKPlus);
        ePosPr = RecoDecay::e(posTrack.px(), posTrack.py(), posTrack.pz(), MassProton);
        for (const auto& negTrack : negTracksPerColl) {
          negIsPi = selPion(negTrack);
          negIsKa = selKaon(negTrack);
          negIsPr = selProton(negTrack);
          if (!(negIsPi || negIsKa || negIsPr)) {
            continue;
          }

          eNegPi = RecoDecay::e(negTrack.px(), negTrack.py(), negTrack.pz(), MassPiPlus);
          eNegKa = RecoDecay::e(negTrack.px(), negTrack.py(), negTrack.pz(), MassKPlus);
          eNegPr = RecoDecay::e(negTrack.px(), negTrack.py(), negTrack.pz(), MassProton);

          px = posTrack.px() + negTrack.px();
          py = posTrack.py() + negTrack.py();
          pz = posTrack.pz() + negTrack.pz();
          p = RecoDecay::p(px, py, pz); // definition = std::sqrt(px*px+py*py+pz*pz);

          fillTable = false;

          mPhi1020 = -1.0f;
          mKStar892 = -1.0f;
          mKStar892Bar = -1.0f;
          mLambda1520 = -1.0f;
          mLambda1520Bar = -1.0f;

          phi1020Tag = kMassOutside;
          kStar892Tag = kMassOutside;
          kStar892BarTag = kMassOutside;
          lambda1520Tag = kMassOutside;
          lambda1520BarTag = kMassOutside;

          // phi(1020) -> K+ + K-
          if (cfgResonances.c0011doPhi1020 && posIsKa && negIsKa) {
            mPhi1020 = RecoDecay::m(p, ePosKa + eNegKa);
            phi1020Tag = getMassRegionTag(mPhi1020, cfgPhi1020Mass.c0018phi1020LSBLow, cfgPhi1020Mass.c0019phi1020LSBUp, cfgPhi1020Mass.c0016phi1020PeakLow, cfgPhi1020Mass.c0017phi1020PeakUp, cfgPhi1020Mass.c0020phi1020RSBLow, cfgPhi1020Mass.c0021phi1020RSBUp);
            if (phi1020Tag != kMassOutside) {
              fillTable = true;
            }
            resonanceQAPlots.fill(HIST("mPhi1020"), mPhi1020);
          }

          // K(892)*   -> K+ + pi-
          if (cfgResonances.c0012doKStar892 && posIsKa && negIsPi) {
            mKStar892 = RecoDecay::m(p, ePosKa + eNegPi);
            kStar892Tag = getMassRegionTag(mKStar892, cfgKstar892Mass.c0024kstar892LSBLow, cfgKstar892Mass.c0025kstar892LSBUp, cfgKstar892Mass.c0022kstar892PeakLow, cfgKstar892Mass.c0023kstar892PeakUp, cfgKstar892Mass.c0026kstar892RSBLow, cfgKstar892Mass.c0027kstar892RSBUp);
            if (kStar892Tag != kMassOutside) {
              fillTable = true;
            }
            resonanceQAPlots.fill(HIST("mKStar892"), mKStar892);
          }

          // K(892)*Bar -> K- + pi+
          if (cfgResonances.c0013doKStar892Bar && posIsPi && negIsKa) {
            mKStar892Bar = RecoDecay::m(p, ePosPi + eNegKa);
            kStar892BarTag = getMassRegionTag(mKStar892Bar, cfgKstar892Mass.c0024kstar892LSBLow, cfgKstar892Mass.c0025kstar892LSBUp, cfgKstar892Mass.c0022kstar892PeakLow, cfgKstar892Mass.c0023kstar892PeakUp, cfgKstar892Mass.c0026kstar892RSBLow, cfgKstar892Mass.c0027kstar892RSBUp);
            if (kStar892BarTag != kMassOutside) {
              fillTable = true;
            }
            resonanceQAPlots.fill(HIST("mKStar892Bar"), mKStar892Bar);
          }

          // Λ(1520) -> P+ + Ka-
          if (cfgResonances.c0014doLambda1520 && posIsPr && negIsKa) {
            mLambda1520 = RecoDecay::m(p, ePosPr + eNegKa);
            lambda1520Tag = getMassRegionTag(mLambda1520, cfgLambda1520Mass.c0030lambda1520LSBLow, cfgLambda1520Mass.c0031lambda1520LSBUp, cfgLambda1520Mass.c0028lambda1520PeakLow, cfgLambda1520Mass.c0029lambda1520PeakUp, cfgLambda1520Mass.c0032lambda1520RSBLow, cfgLambda1520Mass.c0033lambda1520RSBUp);
            if (lambda1520Tag != kMassOutside) {
              fillTable = true;
            }
            resonanceQAPlots.fill(HIST("mLambda1520"), mLambda1520);
          }

          // Λ(1520)Bar -> PBar + Ka+
          if (cfgResonances.c0015doLambda1520Bar && posIsKa && negIsPr) {
            mLambda1520Bar = RecoDecay::m(p, ePosKa + eNegPr);
            lambda1520BarTag = getMassRegionTag(mLambda1520Bar, cfgLambda1520Mass.c0030lambda1520LSBLow, cfgLambda1520Mass.c0031lambda1520LSBUp, cfgLambda1520Mass.c0028lambda1520PeakLow, cfgLambda1520Mass.c0029lambda1520PeakUp, cfgLambda1520Mass.c0032lambda1520RSBLow, cfgLambda1520Mass.c0033lambda1520RSBUp);
            if (lambda1520BarTag != kMassOutside) {
              fillTable = true;
            }
            resonanceQAPlots.fill(HIST("mLambda1520Bar"), mLambda1520Bar);
          }

          if (!fillTable) {
            continue;
          }
          pt = RecoDecay::pt(px, py);                             // definition = std::sqrt(px*px+py*py);
          eta = RecoDecay::eta(std::array<float, 3>{px, py, pz}); // definition = 0.5 * std::log((p + pz) / (p - pz + 1e-12)); // 1e-12 protect divide-by-zero
          phi = RecoDecay::phi(px, py);                           // definition = std::atan2(py, px);

          resonanceCndts(collision.globalIndex(), posTrack.globalIndex(), negTrack.globalIndex(),
                         pt, eta, phi, px, py, pz,
                         mPhi1020, mKStar892, mKStar892Bar, mLambda1520, mLambda1520Bar,
                         phi1020Tag, kStar892Tag, kStar892BarTag, lambda1520Tag, lambda1520BarTag);
        } // Neg Tracks
      } // Pos Tracks
    } // collision Loop is over
  }
  PROCESS_SWITCH(HParticleCorrelation1ResonanceProducer, processSameEvent, "Process Same event", true);
};

struct HParticleCorrelation2SameEvent {

  Produces<aod::HPCorrColls> derivedCollisions;
  Produces<aod::HPCorrTriggers> derivedTriggers;
  Produces<aod::HPCorrAssocs> derivedAssocs;

  // Hisogram redistry:
  HistogramRegistry seEventQA{"seEventQA", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry seTrackQA{"seTrackQA", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry hhCorrelation{"hhCorrelation", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry hIdCorrelation{"hIdCorrelation", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry hPhiCorrelation{"hPhiCorrelation", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry hKStarCorrelation{"hKStarCorrelation", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry hKStarBarCorrelation{"hKStarBarCorrelation", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry hLambdaCorrelation{"hLambdaCorrelation", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry hLambdaBarCorrelation{"hLambdaBarCorrelation", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry resonanceMassSpectra{"resonanceMassSpectra", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry mixingQA{"mixingQA", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  struct : ConfigurableGroup {
    Configurable<bool> c1001printDebugMessages{"c1001printDebugMessages", false, "Print debug messages"};
  } cfgDebug;

  struct : ConfigurableGroup {
    // Basic event selection
    Configurable<float> c1002cutZvertex{"c1002cutZvertex", 10.f, "Maximum |z_{vtx}| (cm)"};        // Imp : All
    Configurable<bool> c1003requireSel8{"c1003requireSel8", true, "Require sel8 event selection"}; // Imp : All
    Configurable<bool> c1004requireTriggerTVX{"c1004requireTriggerTVX", false, "Require kIsTriggerTVX"};

    Configurable<int> c1005minNFilteredTracks{"c1005minNFilteredTracks", 3, "Minimum number of filtered tracks required per collision"};
    Configurable<int> c1006minNSelectedTracks{"c1006minNSelectedTracks", 3, "Minimum number of selected tracks required per collision"};

    // Time-frame / ROF border selections
    Configurable<bool> c1007requireNoITSROFrameBorder{"c1007requireNoITSROFrameBorder", false, "Require kNoITSROFrameBorder"};
    Configurable<bool> c1008requireNoTimeFrameBorder{"c1008requireNoTimeFrameBorder", false, "Require kNoTimeFrameBorder"};

    // Vertex quality
    Configurable<bool> c1009requireVertexITSTPC{"c1009requireVertexITSTPC", true, "Require kIsVertexITSTPC"};           // Imp : Light Ions + pp
    Configurable<bool> c1010requireGoodZvtxFT0vsPV{"c1010requireGoodZvtxFT0vsPV", false, "Require kIsGoodZvtxFT0vsPV"}; // Imp : Light Ions
    Configurable<bool> c1011requireVertexTOFmatched{"c1011requireVertexTOFmatched", false, "Require kIsVertexTOFmatched"};
    Configurable<bool> c1012requireVertexTRDmatched{"c1012requireVertexTRDmatched", false, "Require kIsVertexTRDmatched"};
    Configurable<bool> c1013requireGoodITSLayersAll{"c1013requireGoodITSLayersAll", false, "Require kIsGoodITSLayersAll"}; // Imp : Light Ions

    // Pileup / neighbouring-collision rejection
    Configurable<bool> c1014requireNoSameBunchPileup{"c1014requireNoSameBunchPileup", true, "Require kNoSameBunchPileup"}; // Imp : Light Ions + pp
    Configurable<bool> c1015requireNoCollInTimeRangeStandard{"c1015requireNoCollInTimeRangeStandard", false, "Require kNoCollInTimeRangeStandard"};
    Configurable<bool> c1016requireNoCollInTimeRangeStrict{"c1016requireNoCollInTimeRangeStrict", false, "Require kNoCollInTimeRangeStrict"};
    Configurable<bool> c1017requireNoCollInTimeRangeNarrow{"c1017requireNoCollInTimeRangeNarrow", false, "Require kNoCollInTimeRangeNarrow"};
    Configurable<bool> c1018requireNoCollInRofStandard{"c1018requireNoCollInRofStandard", false, "Require kNoCollInRofStandard"};
    Configurable<bool> c1019requireNoCollInRofStrict{"c1019requireNoCollInRofStrict", false, "Require kNoCollInRofStrict"};
    Configurable<bool> c1020requireNoHighMultCollInPrevRof{"c1020requireNoHighMultCollInPrevRof", false, "Require kNoHighMultCollInPrevRof"};

    // INEL event classes
    Configurable<bool> c1021requireINELgt0{"c1021requireINELgt0", false, "Require INEL > 0"};
    Configurable<bool> c1022requireINELgt1{"c1022requireINELgt1", false, "Require INEL > 1"};

    // Multiplicity selection
    Configurable<bool> c1023useMultiplicitySelection{"c1023useMultiplicitySelection", false, "Apply multiplicity selection"};
    Configurable<int> c1024multiplicityEstimator{"c1024multiplicityEstimator", 0, "0: multNTracksPV, 1: numContrib, 2: multFT0C, 3: multFT0M"};
    Configurable<float> c1025minMultiplicity{"c1025minMultiplicity", 0.f, "Minimum multiplicity"};
    Configurable<float> c1026maxMultiplicity{"c1026maxMultiplicity", 1.e9f, "Maximum multiplicity"};

    // Centrality selection
    Configurable<bool> c1027useCentralitySelection{"c1027useCentralitySelection", false, "Apply centrality selection"};
    Configurable<int> c1028centralityEstimator{"c1028centralityEstimator", 0, "0: centFT0C, 1: centFT0M, 2: centFT0A, 3: centFV0A"};
    Configurable<float> c1029minCentrality{"c1029minCentrality", 0.f, "Minimum centrality percentile"};
    Configurable<float> c1030maxCentrality{"c1030maxCentrality", 100.f, "Maximum centrality percentile"};

    // Occupancy selection
    Configurable<bool> c1031useOccupancySelection{"c1031useOccupancySelection", false, "Apply event occupancy selection"};
    Configurable<bool> c1032useFT0CbasedOccupancy{"c1032useFT0CbasedOccupancy", false, "Use FT0C occupancy instead of track occupancy"};
    Configurable<float> c1033minOccupancy{"c1033minOccupancy", -1.f, "Minimum occupancy"};
    Configurable<float> c1034maxOccupancy{"c1034maxOccupancy", 1.e9f, "Maximum occupancy"};

    // Interaction-rate selection
    Configurable<bool> c1035useInteractionRateSelection{"c1035useInteractionRateSelection", false, "Apply interaction-rate selection"};
    Configurable<float> c1036minInteractionRate{"c1036minInteractionRate", -1.f, "Minimum interaction rate"};
    Configurable<float> c1037maxInteractionRate{"c1037maxInteractionRate", 1.e9f, "Maximum interaction rate"};

    // RCT / detector-quality selection
    Configurable<bool> c1038requireRCTFlagChecker{"c1038requireRCTFlagChecker", false, "Apply Run Condition Table event-quality selection"};
    Configurable<bool> c1039requireCorrelationAnalysisRCTFlagChecker{"c1039requireCorrelationAnalysisRCTFlagChecker", false, "Apply correlation-analysis RCT selection"};
    Configurable<std::string> c1040rctFlagCheckerLabel{"c1040rctFlagCheckerLabel", "CBT_muon_global", "RCT flag checker label"};
    Configurable<bool> c1041rctCheckZDC{"c1041rctCheckZDC", false, "Include ZDC in RCT detector-quality check"};
    Configurable<bool> c1042rctTreatLimitedAcceptanceAsBad{"c1042rctTreatLimitedAcceptanceAsBad", false, "Treat limited detector acceptance as bad"};
  } cfgEvent;

  struct : ConfigurableGroup {
    // Kinematics
    Configurable<float> c1043ptMin{"c1043ptMin", 0.2f, "Minimum track pT (GeV/c)"};
    Configurable<float> c1044ptMax{"c1044ptMax", 1.e10f, "Maximum track pT (GeV/c)"};
    Configurable<float> c1045etaMax{"c1045etaMax", 0.8f, "Maximum |eta|"};

    // Standard O2 track-selection flags
    Configurable<bool> c1046requireGlobalTrack{"c1046requireGlobalTrack", false, "Require isGlobalTrack()"};
    Configurable<bool> c1047requireGlobalTrackWoDCA{"c1047requireGlobalTrackWoDCA", false, "Require isGlobalTrackWoDCA()"};
    Configurable<bool> c1048requirePVContributor{"c1048requirePVContributor", false, "Require track to be a PV contributor"};

    // Detector matching / presence
    Configurable<bool> c1049requireITS{"c1049requireITS", false, "Require ITS information"};
    Configurable<bool> c1050requireTPC{"c1050requireTPC", false, "Require TPC information"};
    Configurable<bool> c1051requireTOF{"c1051requireTOF", false, "Require TOF information"};
    Configurable<bool> c1052requireTRD{"c1052requireTRD", false, "Require TRD information"};

    // TPC quality
    Configurable<int> c1053tpcNClsFoundMin{"c1053tpcNClsFoundMin", 0, "Minimum number of found TPC clusters"};
    Configurable<int> c1054tpcNClsCrossedRowsMin{"c1054tpcNClsCrossedRowsMin", 80, "Minimum number of TPC crossed rows"};
    Configurable<float> c1055tpcCrossedRowsOverFindableMin{"c1055tpcCrossedRowsOverFindableMin", 0.f, "Minimum TPC crossed rows / findable clusters"};
    Configurable<float> c1056tpcFoundOverFindableMin{"c1056tpcFoundOverFindableMin", 0.f, "Minimum TPC found / findable clusters"};
    Configurable<float> c1057tpcFractionSharedMax{"c1057tpcFractionSharedMax", 1.f, "Maximum fraction of shared TPC clusters"};
    Configurable<float> c1058tpcChi2NClMin{"c1058tpcChi2NClMin", 0.f, "Minimum TPC chi2 per cluster"};
    Configurable<float> c1059tpcChi2NClMax{"c1059tpcChi2NClMax", 1.e10f, "Maximum TPC chi2 per cluster"};

    // ITS quality
    Configurable<int> c1060itsNClsMin{"c1060itsNClsMin", 0, "Minimum number of ITS clusters"};
    Configurable<int> c1061itsNClsMax{"c1061itsNClsMax", 7, "Maximum number of ITS clusters"};
    Configurable<int> c1062itsNClsInnerBarrelMin{"c1062itsNClsInnerBarrelMin", 0, "Minimum number of ITS inner-barrel clusters"};
    Configurable<float> c1063itsChi2NClMin{"c1063itsChi2NClMin", 0.f, "Minimum ITS chi2 per cluster"};
    Configurable<float> c1064itsChi2NClMax{"c1064itsChi2NClMax", 1.e10f, "Maximum ITS chi2 per cluster"};

    // Fixed DCA
    Configurable<bool> c1065useFixedDCAxy{"c1065useFixedDCAxy", true, "Apply fixed |DCAxy| cut"};
    Configurable<float> c1066dcaXYMax{"c1066dcaXYMax", 0.1f, "Maximum fixed |DCAxy| (cm)"};
    Configurable<bool> c1067useFixedDCAz{"c1067useFixedDCAz", true, "Apply fixed |DCAz| cut"};
    Configurable<float> c1068dcaZMax{"c1068dcaZMax", 0.2f, "Maximum fixed |DCAz| (cm)"};

    // pT-dependent DCAxy:
    // |DCAxy| < A + B * pT^C
    Configurable<bool> c1069usePtDependentDCAxy{"c1069usePtDependentDCAxy", false, "Apply pT-dependent DCAxy cut"};
    Configurable<float> c1070dcaXYPtA{"c1070dcaXYPtA", 0.0105f, "A coefficient of pT-dependent DCAxy cut"};
    Configurable<float> c1071dcaXYPtB{"c1071dcaXYPtB", 0.035f, "B coefficient of pT-dependent DCAxy cut"};
    Configurable<float> c1072dcaXYPtC{"c1072dcaXYPtC", 1.1f, "Power C of pT-dependent DCAxy cut"};

    // pT-dependent DCAz:
    // |DCAz| < A + B * pT^C
    Configurable<bool> c1073usePtDependentDCAz{"c1073usePtDependentDCAz", false, "Apply pT-dependent DCAz cut"};
    Configurable<float> c1074dcaZPtA{"c1074dcaZPtA", 0.1f, "A coefficient of pT-dependent DCAz cut"};
    Configurable<float> c1075dcaZPtB{"c1075dcaZPtB", 0.0f, "B coefficient of pT-dependent DCAz cut"};
    Configurable<float> c1076dcaZPtC{"c1076dcaZPtC", 0.0f, "Power C of pT-dependent DCAz cut"};
  } cfgTrackCuts;

  struct : ConfigurableGroup {
    Configurable<float> c1077triggerPtLow{"c1077triggerPtLow", 4.0f, "Minimum pT for trigger tracks"};
    Configurable<float> c1078triggerPtHigh{"c1078triggerPtHigh", 8.0f, "Maximum pT for trigger tracks"};
    Configurable<float> c1079assocPtLowMin{"c1079assocPtLowMin", 0.0f, "Minimum pT for low-pT associated tracks"};
    Configurable<float> c1080assocPtLowMax{"c1080assocPtLowMax", 2.0f, "Maximum pT for low-pT associated tracks"};
    Configurable<float> c1081assocPtHighMin{"c1081assocPtHighMin", 2.0f, "Minimum pT for high-pT associated tracks"};
    Configurable<float> c1082assocPtHighMax{"c1082assocPtHighMax", 4.0f, "Maximum pT for high-pT associated tracks"};
  } cfgPartitions;

  struct : ConfigurableGroup {
    Configurable<float> c1083resoRapidityMin{"c1083resoRapidityMin", -0.5f, "Minimum rapidity for resonance candidates"};
    Configurable<float> c1084resoRapidityMax{"c1084resoRapidityMax", 0.5f, "Maximum rapidity for resonance candidates"};

    Configurable<float> c1085phiPtLowMin{"c1085phiPtLowMin", 0.0f, "Minimum pT for low-pT Phi(1020)"};
    Configurable<float> c1086phiPtLowMax{"c1086phiPtLowMax", 2.0f, "Maximum pT for low-pT Phi(1020)"};
    Configurable<float> c1087phiPtHighMin{"c1087phiPtHighMin", 2.0f, "Minimum pT for high-pT Phi(1020)"};
    Configurable<float> c1088phiPtHighMax{"c1088phiPtHighMax", 4.0f, "Maximum pT for high-pT Phi(1020)"};

    Configurable<float> c1089kstarPtLowMin{"c1089kstarPtLowMin", 0.0f, "Minimum pT for low-pT K*(892)0"};
    Configurable<float> c1090kstarPtLowMax{"c1090kstarPtLowMax", 2.0f, "Maximum pT for low-pT K*(892)0"};
    Configurable<float> c1091kstarPtHighMin{"c1091kstarPtHighMin", 2.0f, "Minimum pT for high-pT K*(892)0"};
    Configurable<float> c1092kstarPtHighMax{"c1092kstarPtHighMax", 4.0f, "Maximum pT for high-pT K*(892)0"};

    Configurable<float> c1093kstarBarPtLowMin{"c1093kstarBarPtLowMin", 0.0f, "Minimum pT for low-pT anti-K*(892)0"};
    Configurable<float> c1094kstarBarPtLowMax{"c1094kstarBarPtLowMax", 2.0f, "Maximum pT for low-pT anti-K*(892)0"};
    Configurable<float> c1095kstarBarPtHighMin{"c1095kstarBarPtHighMin", 2.0f, "Minimum pT for high-pT anti-K*(892)0"};
    Configurable<float> c1096kstarBarPtHighMax{"c1096kstarBarPtHighMax", 4.0f, "Maximum pT for high-pT anti-K*(892)0"};

    Configurable<float> c1097lambda1520PtLowMin{"c1097lambda1520PtLowMin", 0.0f, "Minimum pT for low-pT Lambda(1520)"};
    Configurable<float> c1098lambda1520PtLowMax{"c1098lambda1520PtLowMax", 2.0f, "Maximum pT for low-pT Lambda(1520)"};
    Configurable<float> c1099lambda1520PtHighMin{"c1099lambda1520PtHighMin", 2.0f, "Minimum pT for high-pT Lambda(1520)"};
    Configurable<float> c1100lambda1520PtHighMax{"c1100lambda1520PtHighMax", 4.0f, "Maximum pT for high-pT Lambda(1520)"};

    Configurable<float> c1101lambda1520BarPtLowMin{"c1101lambda1520BarPtLowMin", 0.0f, "Minimum pT for low-pT anti-Lambda(1520)"};
    Configurable<float> c1102lambda1520BarPtLowMax{"c1102lambda1520BarPtLowMax", 2.0f, "Maximum pT for low-pT anti-Lambda(1520)"};
    Configurable<float> c1103lambda1520BarPtHighMin{"c1103lambda1520BarPtHighMin", 2.0f, "Minimum pT for high-pT anti-Lambda(1520)"};
    Configurable<float> c1104lambda1520BarPtHighMax{"c1104lambda1520BarPtHighMax", 4.0f, "Maximum pT for high-pT anti-Lambda(1520)"};
  } cfgResPartitions;

  struct : ConfigurableGroup {
    ConfigurableAxis c1105axisDeltaPhi{"c1105axisDeltaPhi", {80, -2.0f, 6.0f}, "#Delta#varphi"};
    ConfigurableAxis c1106axisDeltaEta{"c1106axisDeltaEta", {84, -2.1f, 2.1f}, "#Delta#eta"};

    ConfigurableAxis c1107axisCorrSparseTriggerPt{"c1107axisCorrSparseTriggerPt", {8, 4.0f, 8.0f}, "#it{p}_{T}^{trig} (GeV/#it{c})"};
    ConfigurableAxis c1108axisCorrSparseAssocPt{"c1108axisCorrSparseAssocPt", {8, 0.0f, 4.0f}, "#it{p}_{T}^{assoc} (GeV/#it{c})"};
    ConfigurableAxis c1109axisCorrSparseDeltaPhi{"c1109axisCorrSparseDeltaPhi", {32, -2.0f, 6.0f}, "#Delta#varphi"};
    ConfigurableAxis c1110axisCorrSparseDeltaEta{"c1110axisCorrSparseDeltaEta", {24, -2.1f, 2.1f}, "#Delta#eta"};
    ConfigurableAxis c1111axisCorrCentrality{"c1111axisCorrCentrality", {VARIABLE_WIDTH, 0.0, 1.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0}, "Centrality (%) axis of the correlation histograms (estimator chosen by c1138mixingEstimator); keep equal to c2008"};
  } cfgAxis;

  Configurable<bool> c1112fillDauQAOnce{"c1112fillDauQAOnce", false, "Fill each resonance daughter in QA only once per event"};
  Configurable<bool> c1113rejectResoWithAnyTriggerDaughter{"c1113rejectResoWithAnyTriggerDaughter", false, "Reject resonance if either daughter belongs to the event trigger population"};
  Configurable<bool> c1114requireSelectedTriggerForInvariantMass{"c1114requireSelectedTriggerForInvariantMass", true, "Require at least one selected trigger track before filling resonance invariant-mass spectra"};

  struct : ConfigurableGroup {
    Configurable<float> c1115phi1020PeakLow{"c1115phi1020PeakLow", 1.013f, "Phi(1020) peak lower mass"};
    Configurable<float> c1116phi1020PeakUp{"c1116phi1020PeakUp", 1.026f, "Phi(1020) peak upper mass"};
    Configurable<float> c1117phi1020LSBLow{"c1117phi1020LSBLow", 0.995f, "Phi(1020) LSB lower mass"};
    Configurable<float> c1118phi1020LSBUp{"c1118phi1020LSBUp", 1.005f, "Phi(1020) LSB upper mass"};
    Configurable<float> c1119phi1020RSBLow{"c1119phi1020RSBLow", 1.040f, "Phi(1020) RSB lower mass"};
    Configurable<float> c1120phi1020RSBUp{"c1120phi1020RSBUp", 1.060f, "Phi(1020) RSB upper mass"};
  } cfgPhi1020CorrMass;

  struct : ConfigurableGroup {
    Configurable<float> c1121kStar892PeakLow{"c1121kStar892PeakLow", 0.846f, "K*(892)0 correlation peak lower mass"};
    Configurable<float> c1122kStar892PeakUp{"c1122kStar892PeakUp", 0.946f, "K*(892)0 correlation peak upper mass"};
    Configurable<float> c1123kStar892LSBLow{"c1123kStar892LSBLow", 0.700f, "K*(892)0 correlation LSB lower mass"};
    Configurable<float> c1124kStar892LSBUp{"c1124kStar892LSBUp", 0.780f, "K*(892)0 correlation LSB upper mass"};
    Configurable<float> c1125kStar892RSBLow{"c1125kStar892RSBLow", 1.010f, "K*(892)0 correlation RSB lower mass"};
    Configurable<float> c1126kStar892RSBUp{"c1126kStar892RSBUp", 1.090f, "K*(892)0 correlation RSB upper mass"};
  } cfgKStar892CorrMass;

  struct : ConfigurableGroup {
    Configurable<float> c1127lambda1520PeakLow{"c1127lambda1520PeakLow", 1.500f, "Lambda(1520) correlation peak lower mass"};
    Configurable<float> c1128lambda1520PeakUp{"c1128lambda1520PeakUp", 1.540f, "Lambda(1520) correlation peak upper mass"};
    Configurable<float> c1129lambda1520LSBLow{"c1129lambda1520LSBLow", 1.440f, "Lambda(1520) correlation LSB lower mass"};
    Configurable<float> c1130lambda1520LSBUp{"c1130lambda1520LSBUp", 1.480f, "Lambda(1520) correlation LSB upper mass"};
    Configurable<float> c1131lambda1520RSBLow{"c1131lambda1520RSBLow", 1.560f, "Lambda(1520) correlation RSB lower mass"};
    Configurable<float> c1132lambda1520RSBUp{"c1132lambda1520RSBUp", 1.600f, "Lambda(1520) correlation RSB upper mass"};
  } cfgLambda1520CorrMass;

  struct : ConfigurableGroup {
    Configurable<uint64_t> c1133pairMask{"c1133pairMask", 0ULL, "Correlation mask to store; 0 = any non-zero correlation"};
    Configurable<int> c1134mixingBin{"c1134mixingBin", -1, "Mixing bin to store; -1 = all valid mixing bins"};
    Configurable<bool> c1135requireAllPairBits{"c1135requireAllPairBits", false, "Require all requested pair-mask bits instead of any requested bit"};
    Configurable<bool> c1136resetGlobalCountersPerDF{"c1136resetGlobalCountersPerDF", false, "Reset derived-data global counters at the beginning of every dataframe"};
  } cfgDerivedData;

  struct : ConfigurableGroup {
    Configurable<int> c1137nEvtMixing{"c1137nEvtMixing", 5, "Number of events to mix"};
    Configurable<int> c1138mixingEstimator{"c1138mixingEstimator", kMixCentFT0C, "Mixing estimator: 0=FT0C, 1=FT0M, 2=FT0A, 3=FV0A"};
    ConfigurableAxis c1139axisVtxMixing{"c1139axisVtxMixing", {VARIABLE_WIDTH, -10.0, -8.0, -6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0}, "Mixing bins - z vertex"};
    ConfigurableAxis c1140axisCentMixing{"c1140axisCentMixing", {VARIABLE_WIDTH, -1.0, 20.0, 50.0, 80.0, 101.0}, "Mixing bins - centrality"};
    ConfigurableAxis c1141axisMixingOccupancy{"c1141axisMixingOccupancy", {101, -0.5, 100.5}, "Number of collisions in mixing pool"};
  } cfgMixing;

  using BinningTypeVtxZFT0C = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
  using BinningTypeVtxZFT0M = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;
  using BinningTypeVtxZFT0A = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0A>;
  using BinningTypeVtxZFV0A = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFV0A>;

  BinningTypeVtxZFT0C colBinningFT0C{{cfgMixing.c1139axisVtxMixing, cfgMixing.c1140axisCentMixing}, true};
  BinningTypeVtxZFT0M colBinningFT0M{{cfgMixing.c1139axisVtxMixing, cfgMixing.c1140axisCentMixing}, true};
  BinningTypeVtxZFT0A colBinningFT0A{{cfgMixing.c1139axisVtxMixing, cfgMixing.c1140axisCentMixing}, true};
  BinningTypeVtxZFV0A colBinningFV0A{{cfgMixing.c1139axisVtxMixing, cfgMixing.c1140axisCentMixing}, true};

  int nMixBins = 0;

  template <typename C>
  int getMixingBin(const C& collision)
  {
    switch (cfgMixing.c1138mixingEstimator) {
      case kMixCentFT0C:
        return colBinningFT0C.getBin({collision.posZ(), collision.centFT0C()});
      case kMixCentFT0M:
        return colBinningFT0M.getBin({collision.posZ(), collision.centFT0M()});
      case kMixCentFT0A:
        return colBinningFT0A.getBin({collision.posZ(), collision.centFT0A()});
      case kMixCentFV0A:
        return colBinningFV0A.getBin({collision.posZ(), collision.centFV0A()});
      default:
        return -1;
    }
  }

  template <typename C>
  float getMixingEstimatorValue(const C& collision)
  {
    switch (cfgMixing.c1138mixingEstimator) {
      case kMixCentFT0C:
        return collision.centFT0C();
      case kMixCentFT0M:
        return collision.centFT0M();
      case kMixCentFT0A:
        return collision.centFT0A();
      case kMixCentFV0A:
        return collision.centFV0A();
      default:
        return -1.0f;
    }
  }

  struct : ConfigurableGroup {
    Configurable<LabeledArray<double>> c1142pidConfigSetting{"c1142pidConfigSetting", {DefaultPIDcheckValues[0].data(), kNPid, kNCutSettings, {"Pi", "Ka", "Pr", "El", "Mu", "De"}, {"ThrPforTOF", "IdCutTypeLowP", "NSigmaTPCLowP", "NSigmaTOFLowP", "NSigmaRadLowP", "IdCutTypeHighP", "NSigmaTPCHighP", "NSigmaTOFHighP", "NSigmaRadHighP", "doVetoOthers", "doRelativeTPCCheck", "doRelativeTOFcheck", "doRelativeTPCTOFcheck"}}, "Cut values for particle identification"};
    Configurable<LabeledArray<double>> c1143pidVetoSetting{"c1143pidVetoSetting", {DefaultPidVetoValues[0].data(), kNPid, kNVetoSettings, {"Pi", "Ka", "Pr", "El", "Mu", "De"}, {"doVetoTPC", "doVetoTOF", "vetoTPC", "vetoTOF"}}, "Veto cuts for particle identification"};
    Configurable<bool> c1144cfgId07CheckTofBeta{"c1144cfgId07CheckTofBeta", false, "Require beta > 0 for reliable TOF"};
  } cfgIdCut;

  void init(InitContext const&)
  {
    if (cfgDebug.c1001printDebugMessages) {
      LOGF(info, "Starting init");
    }

    if ((static_cast<uint64_t>(cfgDerivedData.c1133pairMask) & ~(0xFFULL | MaskResoRegions)) != 0ULL) {
      LOG(fatal) << "c1133pairMask = " << static_cast<uint64_t>(cfgDerivedData.c1133pairMask) << " contains bits that are not mixing pools; allowed are bits 0-7 (h-h, h-pi, h-K, h-p) and 18-47 (resonance Peak/LSB/RSB)";
    }
    // Axes
    AxisSpec axisVertexZ = {30, -15., 15., "vrtx_{Z} [cm]"};
    AxisSpec axisCentFT0C = {1200, -10.0, 110.0, "centFT0C(percentile)"};
    AxisSpec axisMult = {150, -1.0, 149.0};

    const AxisSpec axisEventVertexZ = {60, -15.0f, 15.0f, "z_{vtx} (cm)"};
    const AxisSpec axisEventCentrality = {102, -1.0f, 101.0f, "Centrality (%)"};
    const AxisSpec axisMixingBinStatus = {2, -0.5f, 1.5f, "Mixing-bin status"};

    const AxisSpec axisP = {200, 0.0f, 10.0f, "#it{p} (GeV/#it{c})"};
    const AxisSpec axisPt = {200, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisTPCInnerParam = {200, 0.0f, 10.0f, "#it{p}_{tpcInnerParam} (GeV/#it{c})"};
    const AxisSpec axisTOFExpMom = {200, 0.0f, 10.0f, "#it{p}_{tofExpMom} (GeV/#it{c})"};

    const AxisSpec axisEta = {100, -5.0f, 5.0f, "#eta"};
    const AxisSpec axisPhi = {110, -1.0f, 10.0f, "#phi (radians)"};
    const AxisSpec axisRapidity = {200, -5.0f, 5.0f, "Rapidity (y)"};
    const AxisSpec axisDcaXY = {100, -5.0f, 5.0f, "dcaXY"};
    const AxisSpec axisDcaZ = {100, -5.0f, 5.0f, "dcaZ"};
    const AxisSpec axisDcaXYwide = {2000, -100.0f, 100.0f, "dcaXY"};
    const AxisSpec axisDcaZwide = {2000, -100.0f, 100.0f, "dcaZ"};
    const AxisSpec axisSign = {10, -5.0f, 5.0f, "track.sign"};

    const AxisSpec axisTPCSignal = {100, -1.0f, 1000.0f, "tpcSignal"};
    const AxisSpec axisTOFBeta = {40, -2.0f, 2.0f, "tofBeta"};

    const AxisSpec axisTPCSignalFine = {10010, -1.0f, 1000.0f, "tpcSignal"};
    const AxisSpec axisTOFBetaFine = {400, -2.0f, 2.0f, "tofBeta"};

    const AxisSpec axisTPCNSigma = {200, -10.0f, 10.0f, "n#sigma_{TPC}"};
    const AxisSpec axisTOFNSigma = {200, -10.0f, 10.0f, "n#sigma_{TOF}"};
    const AxisSpec axisIdMethod = {2, -0.5f, 1.5f, "ID method (0=TPC, 1=TPC+TOF)"};

    AxisSpec axisDeltaPhiSpec{cfgAxis.c1105axisDeltaPhi, "#Delta#varphi"};
    AxisSpec axisDeltaEtaSpec{cfgAxis.c1106axisDeltaEta, "#Delta#eta"};

    AxisSpec axisCorrSparseTriggerPt{cfgAxis.c1107axisCorrSparseTriggerPt, "#it{p}_{T}^{trig} (GeV/#it{c})"};
    AxisSpec axisCorrSparseAssocPt{cfgAxis.c1108axisCorrSparseAssocPt, "#it{p}_{T}^{assoc} (GeV/#it{c})"};
    AxisSpec axisCorrSparseDeltaPhi{cfgAxis.c1109axisCorrSparseDeltaPhi, "#Delta#varphi"};
    AxisSpec axisCorrSparseDeltaEta{cfgAxis.c1110axisCorrSparseDeltaEta, "#Delta#eta"};
    AxisSpec axisCorrCentrality{cfgAxis.c1111axisCorrCentrality, "Centrality (%)"};

    const AxisSpec axisPhi1020Mass = {300, 0.98f, 1.08f, "m_{K^{+}K^{-}} (GeV/#it{c}^{2})"};
    const AxisSpec axisKStar892Mass = {300, 0.65f, 1.15f, "m_{K#pi} (GeV/#it{c}^{2})"};
    const AxisSpec axisLambda1520Mass = {300, 1.40f, 1.65f, "m_{pK} (GeV/#it{c}^{2})"};
    const AxisSpec axisResoMassCent = {102, -1.0f, 101.0f, "Centrality (%)"};

    AxisSpec axisVtxMixSpec{cfgMixing.c1139axisVtxMixing, "z_{vtx} (cm)"};
    AxisSpec axisCentMixSpec{cfgMixing.c1140axisCentMixing, "Centrality (%)"};
    AxisSpec axisPoolOccupancy{cfgMixing.c1141axisMixingOccupancy, "N eligible collisions"};
    nMixBins = axisVtxMixSpec.getNbins() * axisCentMixSpec.getNbins();

    // Re-build the binnings now that the configurables are set (as data members they hold the default axes)
    colBinningFT0C = BinningTypeVtxZFT0C{{cfgMixing.c1139axisVtxMixing, cfgMixing.c1140axisCentMixing}, true};
    colBinningFT0M = BinningTypeVtxZFT0M{{cfgMixing.c1139axisVtxMixing, cfgMixing.c1140axisCentMixing}, true};
    colBinningFT0A = BinningTypeVtxZFT0A{{cfgMixing.c1139axisVtxMixing, cfgMixing.c1140axisCentMixing}, true};
    colBinningFV0A = BinningTypeVtxZFV0A{{cfgMixing.c1139axisVtxMixing, cfgMixing.c1140axisCentMixing}, true};

    const AxisSpec axisMixBin{nMixBins, -0.5, static_cast<double>(nMixBins) - 0.5, "Mixing bin"};
    const AxisSpec axisCorrChannel{NCorrCountChannels, -0.5, static_cast<double>(NCorrCountChannels) - 0.5, "Correlation channel"};
    const AxisSpec axisReadyBins{nMixBins + 1, -0.5, static_cast<double>(nMixBins) + 0.5, "N ready mixing bins"};

    const int nEventQABins = static_cast<int>(kNEventQABins);
    const int nTrackQABins = static_cast<int>(kNTrackQABins);
    const int nDataFrameQABins = static_cast<int>(kNDataFrameQABins);
    const int nCollisionRejectionBins = static_cast<int>(kNCollisionRejectionTags);
    const int nTrackRejectionBins = static_cast<int>(kNTrackRejectionTags);

    auto addEventQAHistos = [&](auto& histReg, const std::string& basePath) {
      histReg.add((basePath + "VertexZ").c_str(), "Collision vertex;z_{vtx} (cm);Counts", HistType::kTH1F, {axisEventVertexZ});

      histReg.add((basePath + "CentFT0C").c_str(), "FT0C centrality;Centrality (%);Counts", HistType::kTH1F, {axisEventCentrality});
      histReg.add((basePath + "CentFT0M").c_str(), "FT0M centrality;Centrality (%);Counts", HistType::kTH1F, {axisEventCentrality});
      histReg.add((basePath + "CentFT0A").c_str(), "FT0A centrality;Centrality (%);Counts", HistType::kTH1F, {axisEventCentrality});
      histReg.add((basePath + "CentFV0A").c_str(), "FV0A centrality;Centrality (%);Counts", HistType::kTH1F, {axisEventCentrality});

      histReg.add((basePath + "VtxZVsCentFT0C").c_str(), "Vertex z vs FT0C centrality;z_{vtx} (cm);Centrality (%)", HistType::kTH2F, {axisEventVertexZ, axisEventCentrality});
      histReg.add((basePath + "VtxZVsCentFT0M").c_str(), "Vertex z vs FT0M centrality;z_{vtx} (cm);Centrality (%)", HistType::kTH2F, {axisEventVertexZ, axisEventCentrality});
      histReg.add((basePath + "VtxZVsCentFT0A").c_str(), "Vertex z vs FT0A centrality;z_{vtx} (cm);Centrality (%)", HistType::kTH2F, {axisEventVertexZ, axisEventCentrality});
      histReg.add((basePath + "VtxZVsCentFV0A").c_str(), "Vertex z vs FV0A centrality;z_{vtx} (cm);Centrality (%)", HistType::kTH2F, {axisEventVertexZ, axisEventCentrality});

      histReg.add((basePath + "CentFT0CVsCentFT0M").c_str(), "FT0C vs FT0M centrality;FT0C (%);FT0M (%)", HistType::kTH2F, {axisEventCentrality, axisEventCentrality});
      histReg.add((basePath + "CentFT0CVsCentFT0A").c_str(), "FT0C vs FT0A centrality;FT0C (%);FT0A (%)", HistType::kTH2F, {axisEventCentrality, axisEventCentrality});
      histReg.add((basePath + "CentFT0MVsCentFT0A").c_str(), "FT0M vs FT0A centrality;FT0M (%);FT0A (%)", HistType::kTH2F, {axisEventCentrality, axisEventCentrality});

      histReg.add((basePath + "CentFT0CVsCentFV0A").c_str(), "FT0C vs FV0A centrality;FT0C (%);FV0A (%)", HistType::kTH2F, {axisEventCentrality, axisEventCentrality});
      histReg.add((basePath + "CentFT0MVsCentFV0A").c_str(), "FT0M vs FV0A centrality;FT0M (%);FV0A (%)", HistType::kTH2F, {axisEventCentrality, axisEventCentrality});
      histReg.add((basePath + "CentFT0AVsCentFV0A").c_str(), "FT0A vs FV0A centrality;FT0A (%);FV0A (%)", HistType::kTH2F, {axisEventCentrality, axisEventCentrality});

      histReg.add((basePath + "MixingEstimator").c_str(), "Centrality estimator used for mixing;Centrality (%);Counts", HistType::kTH1F, {axisEventCentrality});
      histReg.add((basePath + "VtxZVsMixingEstimator").c_str(), "Physical mixing-bin map;z_{vtx} (cm);Mixing centrality (%)", HistType::kTH2F, {axisVtxMixSpec, axisCentMixSpec});
      histReg.add((basePath + "MixingBin").c_str(), "Mixing-bin index;Mixing bin;Counts", HistType::kTH1F, {axisMixBin});
      histReg.add((basePath + "MixingBinStatus").c_str(), "Mixing-bin validity;Status;Counts", HistType::kTH1F, {axisMixingBinStatus});
    };

    addEventQAHistos(seEventQA, "SE/Events/BeforeSelection/");
    seEventQA.addClone("SE/Events/BeforeSelection/", "SE/Events/AfterSelFunction/");
    seEventQA.addClone("SE/Events/BeforeSelection/", "SE/Events/AfterSelection/");

    auto hMixStatusBefore = seEventQA.get<TH1>(HIST("SE/Events/BeforeSelection/MixingBinStatus"));
    hMixStatusBefore->GetXaxis()->SetBinLabel(1, "Invalid");
    hMixStatusBefore->GetXaxis()->SetBinLabel(2, "Valid");

    auto hMixStatusAfter = seEventQA.get<TH1>(HIST("SE/Events/AfterSelection/MixingBinStatus"));
    hMixStatusAfter->GetXaxis()->SetBinLabel(1, "Invalid");
    hMixStatusAfter->GetXaxis()->SetBinLabel(2, "Valid");

    seEventQA.add("SE/Events/EventSelection", "Event selection QA", HistType::kTH1F, {{nEventQABins, -0.5, static_cast<double>(nEventQABins) - 0.5}});
    seTrackQA.add("SE/Tracks/TrackSelection", "Track selection QA", HistType::kTH1F, {{nTrackQABins, -0.5, static_cast<double>(nTrackQABins) - 0.5}});
    seEventQA.add("SE/Events/DataFrameQA", "Dataframe QA", HistType::kTH1F, {{nDataFrameQABins, -0.5, static_cast<double>(nDataFrameQABins) - 0.5}});
    seEventQA.add("SE/Events/CollisionRejection", "Collision rejection reason", HistType::kTH1F, {{nCollisionRejectionBins, -0.5, static_cast<double>(nCollisionRejectionBins) - 0.5}});
    seTrackQA.add("SE/Tracks/TrackRejection", "Track rejection reason", HistType::kTH1F, {{nTrackRejectionBins, -0.5, static_cast<double>(nTrackRejectionBins) - 0.5}});

    auto hEventQA = seEventQA.get<TH1>(HIST("SE/Events/EventSelection"));
    hEventQA->GetXaxis()->SetBinLabel(static_cast<int>(kEventAllFiltered) + 1, "Filtered collision");
    hEventQA->GetXaxis()->SetBinLabel(static_cast<int>(kEventPassedSelCollision) + 1, "Passed selCollision");
    hEventQA->GetXaxis()->SetBinLabel(static_cast<int>(kEventPassedMinFilteredTracks) + 1, "Passed min filtered tracks");
    hEventQA->GetXaxis()->SetBinLabel(static_cast<int>(kEventPassedMinSelectedTracks) + 1, "Passed min selected tracks");

    auto hTrackQA = seTrackQA.get<TH1>(HIST("SE/Tracks/TrackSelection"));
    hTrackQA->GetXaxis()->SetBinLabel(static_cast<int>(kTrackAllFiltered) + 1, "Filtered track");
    hTrackQA->GetXaxis()->SetBinLabel(static_cast<int>(kTrackPassedSelectionTrack) + 1, "Passed selectionTrack");

    auto hDataFrameQA = seEventQA.get<TH1>(HIST("SE/Events/DataFrameQA"));
    hDataFrameQA->GetXaxis()->SetBinLabel(static_cast<int>(kDFWithZeroFilteredColls) + 1, "Zero filtered collisions");

    auto hCollRej = seEventQA.get<TH1>(HIST("SE/Events/CollisionRejection"));
    hCollRej->GetXaxis()->SetBinLabel(kCollAccepted + 1, "Accepted");
    hCollRej->GetXaxis()->SetBinLabel(kCollRejectSel8 + 1, "sel8");
    hCollRej->GetXaxis()->SetBinLabel(kCollRejectVertexZ + 1, "Vertex Z");
    hCollRej->GetXaxis()->SetBinLabel(kCollRejectTriggerTVX + 1, "TVX");
    hCollRej->GetXaxis()->SetBinLabel(kCollRejectITSROFrameBorder + 1, "ITS ROF border");
    hCollRej->GetXaxis()->SetBinLabel(kCollRejectTimeFrameBorder + 1, "TF border");
    hCollRej->GetXaxis()->SetBinLabel(kCollRejectVertexITSTPC + 1, "ITS-TPC vertex");
    hCollRej->GetXaxis()->SetBinLabel(kCollRejectGoodZvtxFT0vsPV + 1, "FT0/PV z");
    hCollRej->GetXaxis()->SetBinLabel(kCollRejectVertexTOFmatched + 1, "TOF vertex");
    hCollRej->GetXaxis()->SetBinLabel(kCollRejectVertexTRDmatched + 1, "TRD vertex");
    hCollRej->GetXaxis()->SetBinLabel(kCollRejectGoodITSLayersAll + 1, "Good ITS layers");
    hCollRej->GetXaxis()->SetBinLabel(kCollRejectSameBunchPileup + 1, "Same bunch pileup");
    hCollRej->GetXaxis()->SetBinLabel(kCollRejectCollInTimeRangeStandard + 1, "Time range std");
    hCollRej->GetXaxis()->SetBinLabel(kCollRejectCollInTimeRangeStrict + 1, "Time range strict");
    hCollRej->GetXaxis()->SetBinLabel(kCollRejectCollInTimeRangeNarrow + 1, "Time range narrow");
    hCollRej->GetXaxis()->SetBinLabel(kCollRejectCollInRofStandard + 1, "ROF std");
    hCollRej->GetXaxis()->SetBinLabel(kCollRejectCollInRofStrict + 1, "ROF strict");
    hCollRej->GetXaxis()->SetBinLabel(kCollRejectHighMultCollInPrevRof + 1, "High mult prev ROF");
    hCollRej->GetXaxis()->SetBinLabel(kCollRejectINELgt0 + 1, "INEL > 0");
    hCollRej->GetXaxis()->SetBinLabel(kCollRejectINELgt1 + 1, "INEL > 1");

    auto hTrackRej = seTrackQA.get<TH1>(HIST("SE/Tracks/TrackRejection"));

    hTrackRej->GetXaxis()->SetBinLabel(kTrackAccepted + 1, "Accepted");
    hTrackRej->GetXaxis()->SetBinLabel(kTrackRejectPt + 1, "pT");
    hTrackRej->GetXaxis()->SetBinLabel(kTrackRejectEta + 1, "Eta");
    hTrackRej->GetXaxis()->SetBinLabel(kTrackRejectMomentum + 1, "Momentum");
    hTrackRej->GetXaxis()->SetBinLabel(kTrackRejectCharge + 1, "Charge");
    hTrackRej->GetXaxis()->SetBinLabel(kTrackRejectQualityTrack + 1, "Quality track");
    hTrackRej->GetXaxis()->SetBinLabel(kTrackRejectQualityTrackITS + 1, "Quality ITS");
    hTrackRej->GetXaxis()->SetBinLabel(kTrackRejectQualityTrackTPC + 1, "Quality TPC");
    hTrackRej->GetXaxis()->SetBinLabel(kTrackRejectPrimaryTrack + 1, "Primary track");
    hTrackRej->GetXaxis()->SetBinLabel(kTrackRejectInAcceptanceTrack + 1, "Acceptance track");
    hTrackRej->GetXaxis()->SetBinLabel(kTrackRejectGlobalTrack + 1, "Global track");
    hTrackRej->GetXaxis()->SetBinLabel(kTrackRejectGlobalTrackWoDCA + 1, "Global wo DCA");
    hTrackRej->GetXaxis()->SetBinLabel(kTrackRejectGlobalTrackWoPtEta + 1, "Global wo pT eta");
    hTrackRej->GetXaxis()->SetBinLabel(kTrackRejectGlobalTrackWoTPCCluster + 1, "Global wo TPC cluster");
    hTrackRej->GetXaxis()->SetBinLabel(kTrackRejectGlobalTrackWoDCATPCCluster + 1, "Global wo DCA/TPC");
    hTrackRej->GetXaxis()->SetBinLabel(kTrackRejectPVContributor + 1, "PV contributor");
    hTrackRej->GetXaxis()->SetBinLabel(kTrackRejectITS + 1, "No ITS");
    hTrackRej->GetXaxis()->SetBinLabel(kTrackRejectTPC + 1, "No TPC");
    hTrackRej->GetXaxis()->SetBinLabel(kTrackRejectTOF + 1, "No TOF");
    hTrackRej->GetXaxis()->SetBinLabel(kTrackRejectTRD + 1, "No TRD");
    hTrackRej->GetXaxis()->SetBinLabel(kTrackRejectTPCNClsFound + 1, "TPC NCls found");
    hTrackRej->GetXaxis()->SetBinLabel(kTrackRejectTPCCrossedRows + 1, "TPC crossed rows");
    hTrackRej->GetXaxis()->SetBinLabel(kTrackRejectTPCCrossedRowsOverFindable + 1, "TPC crossed/findable");
    hTrackRej->GetXaxis()->SetBinLabel(kTrackRejectTPCFoundOverFindable + 1, "TPC found/findable");
    hTrackRej->GetXaxis()->SetBinLabel(kTrackRejectTPCFractionShared + 1, "TPC shared fraction");
    hTrackRej->GetXaxis()->SetBinLabel(kTrackRejectTPCChi2 + 1, "TPC chi2");
    hTrackRej->GetXaxis()->SetBinLabel(kTrackRejectITSNCls + 1, "ITS NCls");
    hTrackRej->GetXaxis()->SetBinLabel(kTrackRejectITSNClsInnerBarrel + 1, "ITS IB NCls");
    hTrackRej->GetXaxis()->SetBinLabel(kTrackRejectITSChi2 + 1, "ITS chi2");
    hTrackRej->GetXaxis()->SetBinLabel(kTrackRejectFixedDCAxy + 1, "Fixed DCAxy");
    hTrackRej->GetXaxis()->SetBinLabel(kTrackRejectFixedDCAz + 1, "Fixed DCAz");
    hTrackRej->GetXaxis()->SetBinLabel(kTrackRejectPtDependentDCAxy + 1, "pT-dep DCAxy");
    hTrackRej->GetXaxis()->SetBinLabel(kTrackRejectPtDependentDCAz + 1, "pT-dep DCAz");

    //______________________________________________________________________________
    // Correlation histogram booking helpers
    auto addTrackQAHistos = [&](auto& histReg, const std::string& basePath) {
      addBasicTrackQAHistos(histReg, basePath, axisP, axisPt, axisEta, axisPhi, axisDcaXY, axisDcaZ, axisSign);
      histReg.add((basePath + "TPCSignalVsTPCInnerParam").c_str(), "TPC signal vs p_{TPC inner};p_{TPC inner} (GeV/c);TPC signal", kTH2F, {axisTPCInnerParam, axisTPCSignal});
      histReg.add((basePath + "TPCSignalVsP").c_str(), "TPC signal vs p;p (GeV/c);TPC signal", kTH2F, {axisP, axisTPCSignal});
      histReg.add((basePath + "TOFBetaVsTOFExpMom").c_str(), "TOF #beta vs p_{TOF};p_{TOF} (GeV/c);#beta", kTH2F, {axisTOFExpMom, axisTOFBeta});
      histReg.add((basePath + "TOFBetaVsP").c_str(), "TOF #beta vs p;p (GeV/c);#beta", kTH2F, {axisP, axisTOFBeta});
    };

    auto addFullTrackPIDQAHistos = [axisP, axisTPCNSigma, axisTOFNSigma](auto& histReg, const std::string& basePath) {
      histReg.add((basePath + "TPCNSigmaVsP").c_str(), "TPC n#sigma vs p;p (GeV/c);n#sigma_{TPC}", HistType::kTH2F, {axisP, axisTPCNSigma});
      histReg.add((basePath + "TOFNSigmaVsP").c_str(), "TOF n#sigma vs p;p (GeV/c);n#sigma_{TOF}", HistType::kTH2F, {axisP, axisTOFNSigma});
      histReg.add((basePath + "TPCNSigmaVsTOFNSigma").c_str(), "TPC n#sigma vs TOF n#sigma;n#sigma_{TPC};n#sigma_{TOF}", HistType::kTH2F, {axisTPCNSigma, axisTOFNSigma});
    };

    auto addFullTrackQAHistos = [&](auto& histReg, const std::string& basePath) {
      addTrackQAHistos(histReg, basePath);

      histReg.add((basePath + "TPCInnerParam").c_str(), "TPC inner momentum;p_{TPC inner} (GeV/c);Counts", HistType::kTH1F, {axisTPCInnerParam});
      histReg.add((basePath + "TOFExpMom").c_str(), "TOF expected momentum;p_{TOF} (GeV/c);Counts", HistType::kTH1F, {axisTOFExpMom});

      histReg.add((basePath + "DCAxyVsP").c_str(), "DCA_{xy} vs p;p (GeV/c);DCA_{xy} (cm)", HistType::kTH2F, {axisP, axisDcaXY});
      histReg.add((basePath + "DCAxyVsTPCInnerParam").c_str(), "DCA_{xy} vs p_{TPC inner};p_{TPC inner} (GeV/c);DCA_{xy} (cm)", HistType::kTH2F, {axisTPCInnerParam, axisDcaXY});
      histReg.add((basePath + "DCAxyVsTOFExpMom").c_str(), "DCA_{xy} vs p_{TOF};p_{TOF} (GeV/c);DCA_{xy} (cm)", HistType::kTH2F, {axisTOFExpMom, axisDcaXY});

      histReg.add((basePath + "DCAzVsP").c_str(), "DCA_{z} vs p;p (GeV/c);DCA_{z} (cm)", HistType::kTH2F, {axisP, axisDcaZ});
      histReg.add((basePath + "DCAzVsTPCInnerParam").c_str(), "DCA_{z} vs p_{TPC inner};p_{TPC inner} (GeV/c);DCA_{z} (cm)", HistType::kTH2F, {axisTPCInnerParam, axisDcaZ});
      histReg.add((basePath + "DCAzVsTOFExpMom").c_str(), "DCA_{z} vs p_{TOF};p_{TOF} (GeV/c);DCA_{z} (cm)", HistType::kTH2F, {axisTOFExpMom, axisDcaZ});

      histReg.add((basePath + "PVsPt").c_str(), "p vs p_{T};p (GeV/c);p_{T} (GeV/c)", HistType::kTH2F, {axisP, axisPt});
      histReg.add((basePath + "PVsTPCInnerParam").c_str(), "p vs p_{TPC inner};p (GeV/c);p_{TPC inner} (GeV/c)", HistType::kTH2F, {axisP, axisTPCInnerParam});
      histReg.add((basePath + "PVsTOFExpMom").c_str(), "p vs p_{TOF};p (GeV/c);p_{TOF} (GeV/c)", HistType::kTH2F, {axisP, axisTOFExpMom});

      histReg.add((basePath + "TPCSignalVsTOFExpMom").c_str(), "TPC signal vs p_{TOF};p_{TOF} (GeV/c);TPC signal", HistType::kTH2F, {axisTOFExpMom, axisTPCSignal});
      histReg.add((basePath + "TOFBetaVsTPCInnerParam").c_str(), "TOF #beta vs p_{TPC inner};p_{TPC inner} (GeV/c);#beta", HistType::kTH2F, {axisTPCInnerParam, axisTOFBeta});

      addFullTrackPIDQAHistos(histReg, basePath + "PID/Pi/");
      histReg.addClone(basePath + "PID/Pi/", basePath + "PID/Ka/");
      histReg.addClone(basePath + "PID/Pi/", basePath + "PID/Pr/");
      histReg.addClone(basePath + "PID/Pi/", basePath + "PID/El/");
      histReg.addClone(basePath + "PID/Pi/", basePath + "PID/Mu/");
      histReg.addClone(basePath + "PID/Pi/", basePath + "PID/De/");
    };
    addFullTrackQAHistos(seTrackQA, "SE/Tracks/BeforeSelection/");
    seTrackQA.addClone("SE/Tracks/BeforeSelection/", "SE/Tracks/AfterSelection/");

    auto addCorrHistos = [&](auto& histReg, const std::string& basePath) {
      histReg.add((basePath + "DeltaPhi").c_str(), "Trigger-associate;#Delta#varphi;Pairs", HistType::kTH1F, {axisDeltaPhiSpec});
      histReg.add((basePath + "DeltaEta").c_str(), "Trigger-associate;#Delta#eta;Pairs", HistType::kTH1F, {axisDeltaEtaSpec});
      histReg.add((basePath + "DeltaPhiDeltaEta").c_str(), "Trigger-associate;#Delta#varphi;#Delta#eta;Centrality (%)", HistType::kTH3F, {axisDeltaPhiSpec, axisDeltaEtaSpec, axisCorrCentrality});
      histReg.add((basePath + "CorrSparse").c_str(), "Correlation sparse;p_{T}^{trig};p_{T}^{assoc};#Delta#varphi;#Delta#eta;Centrality (%)", HistType::kTHnSparseF, {axisCorrSparseTriggerPt, axisCorrSparseAssocPt, axisCorrSparseDeltaPhi, axisCorrSparseDeltaEta, axisCorrCentrality});
    };

    auto addResoMassCorrHistos = [&](auto& histReg, const std::string& basePath, const AxisSpec& massAxis) {
      histReg.add((basePath + "DeltaPhiVsInvMass").c_str(), "#Delta#varphi vs invariant mass;Invariant mass (GeV/c^{2});#Delta#varphi;Centrality (%)", HistType::kTH3F, {massAxis, axisDeltaPhiSpec, axisCorrCentrality});
      histReg.add((basePath + "DeltaEtaVsInvMass").c_str(), "#Delta#eta vs invariant mass;Invariant mass (GeV/c^{2});#Delta#eta;Centrality (%)", HistType::kTH3F, {massAxis, axisDeltaEtaSpec, axisCorrCentrality});
    };

    auto addResoAssocQAHistos = [&](auto& histReg, const std::string& basePath, const AxisSpec& massAxis, const std::string& massHistTitle) {
      histReg.add((basePath + "P").c_str(), "Resonance;p (GeV/c);Counts", HistType::kTH1F, {axisP});
      histReg.add((basePath + "Pt").c_str(), "Resonance;p_{T} (GeV/c);Counts", HistType::kTH1F, {axisPt});
      histReg.add((basePath + "Eta").c_str(), "Resonance;#eta;Counts", HistType::kTH1F, {axisEta});
      histReg.add((basePath + "Phi").c_str(), "Resonance;#varphi;Counts", HistType::kTH1F, {axisPhi});
      histReg.add((basePath + "Rapidity").c_str(), "Resonance;y;Counts", HistType::kTH1F, {axisRapidity});
      histReg.add((basePath + "InvMass").c_str(), massHistTitle.c_str(), HistType::kTH1F, {massAxis});
      histReg.add((basePath + "InvMassVsPt").c_str(), "Invariant mass vs p_{T};p_{T} (GeV/c);Invariant mass (GeV/c^{2})", HistType::kTH2F, {axisPt, massAxis});
    };

    auto addResoDauQAHistos = [&](auto& histReg, const std::string& basePath) {
      addBasicTrackQAHistos(histReg, basePath, axisP, axisPt, axisEta, axisPhi, axisDcaXY, axisDcaZ, axisSign);
      addIdentifiedQAHistos(histReg, basePath, axisP, axisRapidity, axisTPCSignal, axisTOFBeta, axisTPCNSigma, axisTOFNSigma, axisIdMethod);
    };

    // h-h ______________________________________________________________________________
    addTrackQAHistos(hhCorrelation, "SE/hh/QA/Trigger/");
    hhCorrelation.addClone("SE/hh/QA/Trigger/", "SE/hh/QA/AssocLowPt/");
    hhCorrelation.addClone("SE/hh/QA/Trigger/", "SE/hh/QA/AssocHighPt/");
    addCorrHistos(hhCorrelation, "SE/hh/Corr/AssocLowPt/");
    hhCorrelation.addClone("SE/hh/Corr/AssocLowPt/", "SE/hh/Corr/AssocHighPt/");

    // h-Id ______________________________________________________________________________
    auto registerIdentifiedHistograms = [&](auto& histReg, const std::string& pairName) {
      const std::string base = "SE/" + pairName + "/";

      addTrackQAHistos(histReg, base + "QA/Trigger/");
      addBasicTrackQAHistos(histReg, base + "QA/AssocLowPt/Pi/", axisP, axisPt, axisEta, axisPhi, axisDcaXY, axisDcaZ, axisSign);
      addIdentifiedQAHistos(histReg, base + "QA/AssocLowPt/Pi/", axisP, axisRapidity, axisTPCSignal, axisTOFBeta, axisTPCNSigma, axisTOFNSigma, axisIdMethod);
      histReg.addClone(base + "QA/AssocLowPt/Pi/", base + "QA/AssocLowPt/Ka/");
      histReg.addClone(base + "QA/AssocLowPt/Pi/", base + "QA/AssocLowPt/Pr/");
      histReg.addClone(base + "QA/AssocLowPt/", base + "QA/AssocHighPt/");

      addCorrHistos(histReg, base + "Corr/AssocLowPt/Pi/");
      histReg.addClone(base + "Corr/AssocLowPt/Pi/", base + "Corr/AssocLowPt/Ka/");
      histReg.addClone(base + "Corr/AssocLowPt/Pi/", base + "Corr/AssocLowPt/Pr/");
      histReg.addClone(base + "Corr/AssocLowPt/", base + "Corr/AssocHighPt/");
    };
    registerIdentifiedHistograms(hIdCorrelation, "hIdentified");

    // h-Resonance ______________________________________________________________________________
    auto registerResonanceHistograms = [&](auto& histReg, const std::string& pairName, const AxisSpec& massAxis, const std::string& massHistTitle) {
      const std::string base = "SE/" + pairName + "/";
      const std::string qaPeak = base + "QA/AssocLowPt/Peak/";

      addTrackQAHistos(histReg, base + "QA/Trigger/");
      addResoAssocQAHistos(histReg, qaPeak + "Reso/", massAxis, massHistTitle);
      addResoDauQAHistos(histReg, qaPeak + "PosDau/");
      histReg.addClone(qaPeak + "PosDau/", qaPeak + "NegDau/");
      histReg.addClone(base + "QA/AssocLowPt/Peak/", base + "QA/AssocLowPt/LSB/");
      histReg.addClone(base + "QA/AssocLowPt/Peak/", base + "QA/AssocLowPt/RSB/");
      histReg.addClone(base + "QA/AssocLowPt/", base + "QA/AssocHighPt/");

      addCorrHistos(histReg, base + "Corr/AssocLowPt/Peak/");
      histReg.addClone(base + "Corr/AssocLowPt/Peak/", base + "Corr/AssocLowPt/LSB/");
      histReg.addClone(base + "Corr/AssocLowPt/Peak/", base + "Corr/AssocLowPt/RSB/");
      addResoMassCorrHistos(histReg, base + "Corr/AssocLowPt/MassDependent/", massAxis);
      histReg.addClone(base + "Corr/AssocLowPt/", base + "Corr/AssocHighPt/");
    };

    registerResonanceHistograms(hPhiCorrelation, "hPhi", axisPhi1020Mass, "#phi(1020);m_{K^{+}K^{-}} (GeV/#it{c}^{2});Counts");
    registerResonanceHistograms(hKStarCorrelation, "hKStar", axisKStar892Mass, "K^{*}(892)^{0};m_{K^{+}#pi^{-}} (GeV/#it{c}^{2});Counts");
    registerResonanceHistograms(hKStarBarCorrelation, "hKStarBar", axisKStar892Mass, "#bar{K}^{*}(892)^{0};m_{#pi^{+}K^{-}} (GeV/#it{c}^{2});Counts");
    registerResonanceHistograms(hLambdaCorrelation, "hLambda", axisLambda1520Mass, "#Lambda(1520);m_{pK^{-}} (GeV/#it{c}^{2});Counts");
    registerResonanceHistograms(hLambdaBarCorrelation, "hLambdaBar", axisLambda1520Mass, "#bar{#Lambda}(1520);m_{K^{+}#bar{p}} (GeV/#it{c}^{2});Counts");

    auto addInvariantMassSpectra = [&](const std::string& species, const AxisSpec& massAxis) {
      const std::string seBase = "SE/InvariantMass/" + species + "/";
      const std::string meBase = "ME/InvariantMass/" + species + "/";

      resonanceMassSpectra.add((seBase + "US").c_str(), "Unlike-sign;m_{inv} (GeV/#it{c}^{2});p_{T} (GeV/#it{c});Centrality (%)", kTHnSparseF, {massAxis, axisPt, axisResoMassCent});
      resonanceMassSpectra.add((seBase + "LS").c_str(), "Like-sign;m_{inv} (GeV/#it{c}^{2});p_{T} (GeV/#it{c});Centrality (%)", kTHnSparseF, {massAxis, axisPt, axisResoMassCent});
      resonanceMassSpectra.add((meBase + "US").c_str(), "Mixed-event;m_{inv} (GeV/#it{c}^{2});p_{T} (GeV/#it{c});Centrality (%)", kTHnSparseF, {massAxis, axisPt, axisResoMassCent});
    };

    addInvariantMassSpectra("Phi", axisPhi1020Mass);
    addInvariantMassSpectra("KStar", axisKStar892Mass);
    addInvariantMassSpectra("KStarBar", axisKStar892Mass);
    addInvariantMassSpectra("Lambda1520", axisLambda1520Mass);
    addInvariantMassSpectra("Lambda1520Bar", axisLambda1520Mass);

    //______________________________________________________________________________
    // Event-mixing QA
    mixingQA.add("Mixing/DataFramesProcessed", "Processed dataframes;Dataframe;Count", HistType::kTH1D, {{1, -0.5, 0.5}});
    mixingQA.add("Mixing/AllCollisionsPerBin", "Selected collisions per mixing bin;Mixing bin;Collisions", HistType::kTH1D, {axisMixBin});
    mixingQA.add("Mixing/AllCollisionsZCent", "Selected collisions in physical mixing pools;z_{vtx} (cm);Centrality (%)", HistType::kTH2D, {axisVtxMixSpec, axisCentMixSpec});
    mixingQA.add("Mixing/EligibleCollisionsPerBin", "Eligible collisions per mixing pool;Mixing bin;Correlation channel", HistType::kTH2D, {axisMixBin, axisCorrChannel});
    mixingQA.add("Mixing/EligibleCollisionsZCentChannel", "Eligible collisions in physical mixing pools;z_{vtx} (cm);Centrality (%);Correlation channel", HistType::kTH3D, {axisVtxMixSpec, axisCentMixSpec, axisCorrChannel});
    mixingQA.add("Mixing/PerDF/AllPoolOccupancy", "Mixing-pool occupancy per dataframe;Mixing bin;N collisions", HistType::kTH2D, {axisMixBin, axisPoolOccupancy});
    mixingQA.add("Mixing/PerDF/EligiblePoolOccupancy", "Eligible pool occupancy per dataframe;Mixing bin;Correlation channel;N eligible collisions", HistType::kTH3D, {axisMixBin, axisCorrChannel, axisPoolOccupancy});
    mixingQA.add("Mixing/PerDF/PoolOccupancyDistribution", "Distribution of mixing-pool occupancies;Correlation channel;N eligible collisions", HistType::kTH2D, {axisCorrChannel, axisPoolOccupancy});
    mixingQA.add("Mixing/PerDF/ReadyPoolMap", "Dataframes where pool was ready;Mixing bin;Correlation channel", HistType::kTH2D, {axisMixBin, axisCorrChannel});
    mixingQA.add("Mixing/PerDF/NReadyMixBins", "Number of ready mixing bins per dataframe;Correlation channel;N ready mixing bins", HistType::kTH2D, {axisCorrChannel, axisReadyBins});

    auto setCorrChannelLabels = [](TAxis* axis) {
      for (int iCorr = 0; iCorr < kNCorrCountTypes; ++iCorr) {
        axis->SetBinLabel(iCorr * kNCorrCountPt + kCountLowPt + 1, (std::string(CorrCountTypeName[iCorr]) + "_LowPt").c_str());
        axis->SetBinLabel(iCorr * kNCorrCountPt + kCountHighPt + 1, (std::string(CorrCountTypeName[iCorr]) + "_HighPt").c_str());
      }
    };
    setCorrChannelLabels(mixingQA.get<TH2>(HIST("Mixing/EligibleCollisionsPerBin"))->GetYaxis());
    setCorrChannelLabels(mixingQA.get<TH3>(HIST("Mixing/EligibleCollisionsZCentChannel"))->GetZaxis());
    setCorrChannelLabels(mixingQA.get<TH3>(HIST("Mixing/PerDF/EligiblePoolOccupancy"))->GetYaxis());
    setCorrChannelLabels(mixingQA.get<TH2>(HIST("Mixing/PerDF/PoolOccupancyDistribution"))->GetXaxis());
    setCorrChannelLabels(mixingQA.get<TH2>(HIST("Mixing/PerDF/ReadyPoolMap"))->GetYaxis());
    setCorrChannelLabels(mixingQA.get<TH2>(HIST("Mixing/PerDF/NReadyMixBins"))->GetXaxis());

    mixingQA.add("DerivedData/RowsWritten", "Rows written to the derived tables;Table;Rows", HistType::kTH1D, {{4, -0.5, 3.5}});
    mixingQA.get<TH1>(HIST("DerivedData/RowsWritten"))->GetXaxis()->SetBinLabel(1, "Collisions");
    mixingQA.get<TH1>(HIST("DerivedData/RowsWritten"))->GetXaxis()->SetBinLabel(2, "Triggers");
    mixingQA.get<TH1>(HIST("DerivedData/RowsWritten"))->GetXaxis()->SetBinLabel(3, "TrackAssocs");
    mixingQA.get<TH1>(HIST("DerivedData/RowsWritten"))->GetXaxis()->SetBinLabel(4, "ResoAssocs");

    if (cfgDebug.c1001printDebugMessages) {
      LOGF(info, "Finishing init");
    }
  } // Init function is over.

  //_______________________________ Particle Identification _______________________________
  //   p-dependent identification
  //     p <  ThrPforTOF : TPC+TOF (circular cut) if TOF is reliable, otherwise TPC only
  //     p >= ThrPforTOF : TPC+TOF (circular cut), TOF is required
  // Default cuts are 3 sigma (see DefaultPIDcheckValues).

  // If vetoIdOthers = true; it passed all veto checks
  // If vetoIdOthers = false; it failed veto check with some other particle
  template <int pidMode, typename T>
  bool vetoIdOthersTPC(const T& track)
  {
    // Static is only run once, ever.
    static const std::array<bool, kNPid> doVetoTPC = {
      getCfg<bool>(cfgIdCut.c1143pidVetoSetting, kPi, kDoVetoTPC),
      getCfg<bool>(cfgIdCut.c1143pidVetoSetting, kKa, kDoVetoTPC),
      getCfg<bool>(cfgIdCut.c1143pidVetoSetting, kPr, kDoVetoTPC),
      getCfg<bool>(cfgIdCut.c1143pidVetoSetting, kEl, kDoVetoTPC),
      getCfg<bool>(cfgIdCut.c1143pidVetoSetting, kMu, kDoVetoTPC),
      getCfg<bool>(cfgIdCut.c1143pidVetoSetting, kDe, kDoVetoTPC)};

    static const std::array<float, kNPid> vetoTPC = {
      getCfg<float>(cfgIdCut.c1143pidVetoSetting, kPi, kVetoTPC),
      getCfg<float>(cfgIdCut.c1143pidVetoSetting, kKa, kVetoTPC),
      getCfg<float>(cfgIdCut.c1143pidVetoSetting, kPr, kVetoTPC),
      getCfg<float>(cfgIdCut.c1143pidVetoSetting, kEl, kVetoTPC),
      getCfg<float>(cfgIdCut.c1143pidVetoSetting, kMu, kVetoTPC),
      getCfg<float>(cfgIdCut.c1143pidVetoSetting, kDe, kVetoTPC)};

    return applyVetoOthersTPC<pidMode>(track, doVetoTPC, vetoTPC);
  }

  template <int pidMode, typename T>
  bool vetoIdOthersTOF(const T& track)
  {
    // Only computed once
    static const std::array<bool, kNPid> doVetoTOF = {
      getCfg<bool>(cfgIdCut.c1143pidVetoSetting, kPi, kDoVetoTOF),
      getCfg<bool>(cfgIdCut.c1143pidVetoSetting, kKa, kDoVetoTOF),
      getCfg<bool>(cfgIdCut.c1143pidVetoSetting, kPr, kDoVetoTOF),
      getCfg<bool>(cfgIdCut.c1143pidVetoSetting, kEl, kDoVetoTOF),
      getCfg<bool>(cfgIdCut.c1143pidVetoSetting, kMu, kDoVetoTOF),
      getCfg<bool>(cfgIdCut.c1143pidVetoSetting, kDe, kDoVetoTOF)};

    static const std::array<float, kNPid> vetoTOF = {
      getCfg<float>(cfgIdCut.c1143pidVetoSetting, kPi, kVetoTOF),
      getCfg<float>(cfgIdCut.c1143pidVetoSetting, kKa, kVetoTOF),
      getCfg<float>(cfgIdCut.c1143pidVetoSetting, kPr, kVetoTOF),
      getCfg<float>(cfgIdCut.c1143pidVetoSetting, kEl, kVetoTOF),
      getCfg<float>(cfgIdCut.c1143pidVetoSetting, kMu, kVetoTOF),
      getCfg<float>(cfgIdCut.c1143pidVetoSetting, kDe, kVetoTOF)};

    return applyVetoOthersTOF<pidMode>(track, doVetoTOF, vetoTOF);
  }

  template <int pidMode, typename T>
  bool vetoIdOthersTPCTOF(const T& track)
  {
    // If either veto fails, reject
    return vetoIdOthersTPC<pidMode>(track) && vetoIdOthersTOF<pidMode>(track);
  }

  // Check if TOF is reliable
  template <typename T>
  inline bool checkReliableTOF(const T& track)
  {
    // which check makes the information of TOF relaiable? should track.beta() be checked? e.g.:
    if (cfgIdCut.c1144cfgId07CheckTofBeta) {
      return (track.hasTOF() && track.beta() > 0.0f);
    }
    return track.hasTOF();
  }

  template <int pidMode, typename T>
  bool idTPC(const T& track, const float& nSigmaTPC, float& nSigmaIdDistSq)
  {
    if constexpr (pidMode == kPi) {
      nSigmaIdDistSq = track.tpcNSigmaPi() * track.tpcNSigmaPi();
    } else if constexpr (pidMode == kKa) {
      nSigmaIdDistSq = track.tpcNSigmaKa() * track.tpcNSigmaKa();
    } else if constexpr (pidMode == kPr) {
      nSigmaIdDistSq = track.tpcNSigmaPr() * track.tpcNSigmaPr();
    } else if constexpr (pidMode == kEl) {
      nSigmaIdDistSq = track.tpcNSigmaEl() * track.tpcNSigmaEl();
    } else if constexpr (pidMode == kMu) {
      nSigmaIdDistSq = track.tpcNSigmaMu() * track.tpcNSigmaMu();
    } else if constexpr (pidMode == kDe) {
      nSigmaIdDistSq = track.tpcNSigmaDe() * track.tpcNSigmaDe();
    } else {
      nSigmaIdDistSq = 1000000;
    }

    static const bool doVetoOthers = getCfg<bool>(cfgIdCut.c1142pidConfigSetting, pidMode, kDoVetoOthers);
    if (doVetoOthers) {
      if (!vetoIdOthersTPC<pidMode>(track)) {
        // If vetoIdOthers = true; it passed all veto checks
        // If vetoIdOthers = false; it failed veto check with some other particle
        return false;
      }
    }

    static const bool doRelativeTPCcheck = getCfg<bool>(cfgIdCut.c1142pidConfigSetting, pidMode, kDoRelativeTPCcheck);
    if (doRelativeTPCcheck) {
      if (!relativeIdOthersTPC<pidMode>(track)) {
        // If relativeIdOthersTPC = true; particle has stronger nSigma compared to others
        // If relativeIdOthersTPC = false; some particle has stronger nSigma compared to it
        return false;
      }
    }

    if constexpr (pidMode == kPi) {
      return std::fabs(track.tpcNSigmaPi()) < nSigmaTPC;
    } else if constexpr (pidMode == kKa) {
      return std::fabs(track.tpcNSigmaKa()) < nSigmaTPC;
    } else if constexpr (pidMode == kPr) {
      return std::fabs(track.tpcNSigmaPr()) < nSigmaTPC;
    } else if constexpr (pidMode == kEl) {
      return std::fabs(track.tpcNSigmaEl()) < nSigmaTPC;
    } else if constexpr (pidMode == kMu) {
      return std::fabs(track.tpcNSigmaMu()) < nSigmaTPC;
    } else if constexpr (pidMode == kDe) {
      return std::fabs(track.tpcNSigmaDe()) < nSigmaTPC;
    } else {
      return false;
    }
  }

  template <int pidMode, typename T>
  bool idTPCTOF(const T& track, const int& pidCutType, const float& nSigmaTPC, const float& nSigmaTOF, const float& nSigmaSquaredRad, float& nSigmaIdDistSq)
  {
    if constexpr (pidMode == kPi) {
      nSigmaIdDistSq = track.tpcNSigmaPi() * track.tpcNSigmaPi() + track.tofNSigmaPi() * track.tofNSigmaPi();
    } else if constexpr (pidMode == kKa) {
      nSigmaIdDistSq = track.tpcNSigmaKa() * track.tpcNSigmaKa() + track.tofNSigmaKa() * track.tofNSigmaKa();
    } else if constexpr (pidMode == kPr) {
      nSigmaIdDistSq = track.tpcNSigmaPr() * track.tpcNSigmaPr() + track.tofNSigmaPr() * track.tofNSigmaPr();
    } else if constexpr (pidMode == kEl) {
      nSigmaIdDistSq = track.tpcNSigmaEl() * track.tpcNSigmaEl() + track.tofNSigmaEl() * track.tofNSigmaEl();
    } else if constexpr (pidMode == kMu) {
      nSigmaIdDistSq = track.tpcNSigmaMu() * track.tpcNSigmaMu() + track.tofNSigmaMu() * track.tofNSigmaMu();
    } else if constexpr (pidMode == kDe) {
      nSigmaIdDistSq = track.tpcNSigmaDe() * track.tpcNSigmaDe() + track.tofNSigmaDe() * track.tofNSigmaDe();
    } else {
      nSigmaIdDistSq = 1000000;
    }
    static const bool doVetoOthers = getCfg<bool>(cfgIdCut.c1142pidConfigSetting, pidMode, kDoVetoOthers);
    if (doVetoOthers) {
      if (!vetoIdOthersTPCTOF<pidMode>(track)) {
        // If vetoIdOthers = true; it passed all veto checks
        // If vetoIdOthers = false; it failed veto check with some other particle
        return false;
      }
    }

    static const bool doRelativeTOFcheck = getCfg<bool>(cfgIdCut.c1142pidConfigSetting, pidMode, kDoRelativeTOFcheck);
    if (doRelativeTOFcheck) {
      if (!relativeIdOthersTOF<pidMode>(track)) {
        // If relativeIdOthersTOF = true; particle has stronger nSigma compared to others
        // If relativeIdOthersTOF = false; some particle has stronger nSigma compared to it
        return false;
      }
    }

    static const bool doRelativeTPCTOFcheck = getCfg<bool>(cfgIdCut.c1142pidConfigSetting, pidMode, kDoRelativeTPCTOFcheck);
    if (doRelativeTPCTOFcheck) {
      if (!relativeIdOthersTPCTOF<pidMode>(track)) {
        // If relativeIdOthersTPCTOF = true; particle has stronger nSigma compared to others
        // If relativeIdOthersTPCTOF = false; some particle has stronger nSigma compared to it
        return false;
      }
    }

    switch (pidCutType) {
      case kCircularCut:
        return selIdCircularCut<pidMode>(track, nSigmaSquaredRad);
      case kRectangularCut:
        return selIdRectangularCut<pidMode>(track, nSigmaTPC, nSigmaTOF);
      case kEllipsoidalCut:
        return selIdEllipsoidalCut<pidMode>(track, nSigmaTPC, nSigmaTOF);
      default:
        return false;
    }
  }

  template <int pidMode, typename T>
  bool selPdependent(const T& track, int& IdMethod, float& nSigmaIdDistSq)
  {
    // Static cache inside function - initialized once on first call
    static const auto thrPforTOF = getCfg<float>(cfgIdCut.c1142pidConfigSetting, pidMode, kThrPforTOF);
    static const auto idCutTypeLowP = getCfg<int>(cfgIdCut.c1142pidConfigSetting, pidMode, kIdCutTypeLowP);
    static const auto nSigmaTPCLowP = getCfg<float>(cfgIdCut.c1142pidConfigSetting, pidMode, kNSigmaTPCLowP);
    static const auto nSigmaTOFLowP = getCfg<float>(cfgIdCut.c1142pidConfigSetting, pidMode, kNSigmaTOFLowP);
    static const auto nSigmaRadLowP = getCfg<float>(cfgIdCut.c1142pidConfigSetting, pidMode, kNSigmaRadLowP);
    static const auto idCutTypeHighP = getCfg<int>(cfgIdCut.c1142pidConfigSetting, pidMode, kIdCutTypeHighP);
    static const auto nSigmaTPCHighP = getCfg<float>(cfgIdCut.c1142pidConfigSetting, pidMode, kNSigmaTPCHighP);
    static const auto nSigmaTOFHighP = getCfg<float>(cfgIdCut.c1142pidConfigSetting, pidMode, kNSigmaTOFHighP);
    static const auto nSigmaRadHighP = getCfg<float>(cfgIdCut.c1142pidConfigSetting, pidMode, kNSigmaRadHighP);

    if (track.p() < thrPforTOF) {
      if (checkReliableTOF(track)) {
        if (idTPCTOF<pidMode>(track, idCutTypeLowP, nSigmaTPCLowP, nSigmaTOFLowP, nSigmaRadLowP, nSigmaIdDistSq)) {
          IdMethod = kTPCTOFidentified;
          return true;
        }
      } else {
        if (idTPC<pidMode>(track, nSigmaTPCLowP, nSigmaIdDistSq)) {
          IdMethod = kTPCidentified;
          return true;
        }
      }
    } else {
      if (checkReliableTOF(track)) {
        if (idTPCTOF<pidMode>(track, idCutTypeHighP, nSigmaTPCHighP, nSigmaTOFHighP, nSigmaRadHighP, nSigmaIdDistSq)) {
          IdMethod = kTPCTOFidentified;
          return true;
        }
      }
    }
    return false;
  }

  //______________________________Identification Functions________________________________________________________________
  // Pion
  template <typename T>
  bool selPion(const T& track, int& IdMethod, float& nSigmaIdDistSq)
  {
    return selPdependent<kPi>(track, IdMethod, nSigmaIdDistSq);
  }

  // Kaon
  template <typename T>
  bool selKaon(const T& track, int& IdMethod, float& nSigmaIdDistSq)
  {
    return selPdependent<kKa>(track, IdMethod, nSigmaIdDistSq);
  }

  // Proton
  template <typename T>
  bool selProton(const T& track, int& IdMethod, float& nSigmaIdDistSq)
  {
    return selPdependent<kPr>(track, IdMethod, nSigmaIdDistSq);
  }

  // Electron
  template <typename T>
  bool selElectron(const T& track, int& IdMethod, float& nSigmaIdDistSq)
  {
    return selPdependent<kEl>(track, IdMethod, nSigmaIdDistSq);
  }

  // Muon
  template <typename T>
  bool selMuon(const T& track, int& IdMethod, float& nSigmaIdDistSq)
  {
    return selPdependent<kMu>(track, IdMethod, nSigmaIdDistSq);
  }

  // Deuteron
  template <typename T>
  bool selDeuteron(const T& track, int& IdMethod, float& nSigmaIdDistSq)
  {
    return selPdependent<kDe>(track, IdMethod, nSigmaIdDistSq);
  }

  // Two-argument versions
  template <typename T>
  bool selPion(const T& track, int& IdMethod)
  {
    float nSigmaIdDistSq = -1.0f;
    return selPion(track, IdMethod, nSigmaIdDistSq);
  }

  template <typename T>
  bool selKaon(const T& track, int& IdMethod)
  {
    float nSigmaIdDistSq = -1.0f;
    return selKaon(track, IdMethod, nSigmaIdDistSq);
  }

  template <typename T>
  bool selProton(const T& track, int& IdMethod)
  {
    float nSigmaIdDistSq = -1.0f;
    return selProton(track, IdMethod, nSigmaIdDistSq);
  }

  template <typename T>
  bool selElectron(const T& track, int& IdMethod)
  {
    float nSigmaIdDistSq = -1.0f;
    return selElectron(track, IdMethod, nSigmaIdDistSq);
  }

  template <typename T>
  bool selMuon(const T& track, int& IdMethod)
  {
    float nSigmaIdDistSq = -1.0f;
    return selMuon(track, IdMethod, nSigmaIdDistSq);
  }

  template <typename T>
  bool selDeuteron(const T& track, int& IdMethod)
  {
    float nSigmaIdDistSq = -1.0f;
    return selDeuteron(track, IdMethod, nSigmaIdDistSq);
  }
  //

  template <typename T>
  CollisionRejectionTag selCollision(T const& collision)
  {
    // Basic event selection
    if (cfgEvent.c1003requireSel8 && !collision.sel8()) {
      return kCollRejectSel8;
    }

    if (std::abs(collision.posZ()) >= cfgEvent.c1002cutZvertex) {
      return kCollRejectVertexZ;
    }

    if (cfgEvent.c1004requireTriggerTVX && !collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
      return kCollRejectTriggerTVX;
    }

    // TF / ITS-ROF borders
    if (cfgEvent.c1007requireNoITSROFrameBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return kCollRejectITSROFrameBorder;
    }

    if (cfgEvent.c1008requireNoTimeFrameBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return kCollRejectTimeFrameBorder;
    }

    // Vertex quality
    if (cfgEvent.c1009requireVertexITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return kCollRejectVertexITSTPC;
    }

    if (cfgEvent.c1010requireGoodZvtxFT0vsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return kCollRejectGoodZvtxFT0vsPV;
    }

    if (cfgEvent.c1011requireVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
      return kCollRejectVertexTOFmatched;
    }

    if (cfgEvent.c1012requireVertexTRDmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTRDmatched)) {
      return kCollRejectVertexTRDmatched;
    }

    if (cfgEvent.c1013requireGoodITSLayersAll && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      return kCollRejectGoodITSLayersAll;
    }

    // Pileup / neighbouring-collision rejection
    if (cfgEvent.c1014requireNoSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return kCollRejectSameBunchPileup;
    }

    if (cfgEvent.c1015requireNoCollInTimeRangeStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return kCollRejectCollInTimeRangeStandard;
    }

    if (cfgEvent.c1016requireNoCollInTimeRangeStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict)) {
      return kCollRejectCollInTimeRangeStrict;
    }

    if (cfgEvent.c1017requireNoCollInTimeRangeNarrow && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow)) {
      return kCollRejectCollInTimeRangeNarrow;
    }

    if (cfgEvent.c1018requireNoCollInRofStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      return kCollRejectCollInRofStandard;
    }

    if (cfgEvent.c1019requireNoCollInRofStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStrict)) {
      return kCollRejectCollInRofStrict;
    }

    if (cfgEvent.c1020requireNoHighMultCollInPrevRof && !collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
      return kCollRejectHighMultCollInPrevRof;
    }

    // INEL classes
    if (cfgEvent.c1021requireINELgt0 && !collision.isInelGt0()) {
      return kCollRejectINELgt0;
    }

    if (cfgEvent.c1022requireINELgt1 && !collision.isInelGt1()) {
      return kCollRejectINELgt1;
    }

    return kCollAccepted;
  } // selCollision is Over.

  template <typename T>
  TrackRejectionTag selectionTrack(T const& track)
  {

    // Filter cuts again: resonance daughters are reached through an index into the full track table,
    // which does not re-apply the Filter (the resonance producer uses looser cuts).
    if (track.pt() <= cfgTrackCuts.c1043ptMin || track.pt() >= cfgTrackCuts.c1044ptMax) {
      return kTrackRejectPt;
    }
    if (std::abs(track.eta()) >= cfgTrackCuts.c1045etaMax) {
      return kTrackRejectEta;
    }
    if (cfgTrackCuts.c1065useFixedDCAxy && std::abs(track.dcaXY()) >= cfgTrackCuts.c1066dcaXYMax) {
      return kTrackRejectFixedDCAxy;
    }
    if (cfgTrackCuts.c1067useFixedDCAz && std::abs(track.dcaZ()) >= cfgTrackCuts.c1068dcaZMax) {
      return kTrackRejectFixedDCAz;
    }
    // Filter cuts are over

    // Standard O2 flags
    if (cfgTrackCuts.c1046requireGlobalTrack && !track.isGlobalTrack()) {
      return kTrackRejectGlobalTrack;
    }

    if (cfgTrackCuts.c1047requireGlobalTrackWoDCA && !track.isGlobalTrackWoDCA()) {
      return kTrackRejectGlobalTrackWoDCA;
    }

    if (cfgTrackCuts.c1048requirePVContributor && !track.isPVContributor()) {
      return kTrackRejectPVContributor;
    }

    // Detector presence
    if (cfgTrackCuts.c1049requireITS && !track.hasITS()) {
      return kTrackRejectITS;
    }

    if (cfgTrackCuts.c1050requireTPC && !track.hasTPC()) {
      return kTrackRejectTPC;
    }

    if (cfgTrackCuts.c1051requireTOF && !track.hasTOF()) {
      return kTrackRejectTOF;
    }

    if (cfgTrackCuts.c1052requireTRD && !track.hasTRD()) {
      return kTrackRejectTRD;
    }

    // TPC quality
    if (track.tpcNClsFound() < cfgTrackCuts.c1053tpcNClsFoundMin) {
      return kTrackRejectTPCNClsFound;
    }

    if (track.tpcNClsCrossedRows() < cfgTrackCuts.c1054tpcNClsCrossedRowsMin) {
      return kTrackRejectTPCCrossedRows;
    }

    if (track.tpcCrossedRowsOverFindableCls() < cfgTrackCuts.c1055tpcCrossedRowsOverFindableMin) {
      return kTrackRejectTPCCrossedRowsOverFindable;
    }

    if (track.tpcFoundOverFindableCls() < cfgTrackCuts.c1056tpcFoundOverFindableMin) {
      return kTrackRejectTPCFoundOverFindable;
    }

    if (track.tpcFractionSharedCls() > cfgTrackCuts.c1057tpcFractionSharedMax) {
      return kTrackRejectTPCFractionShared;
    }

    if (track.tpcChi2NCl() < cfgTrackCuts.c1058tpcChi2NClMin || track.tpcChi2NCl() > cfgTrackCuts.c1059tpcChi2NClMax) {
      return kTrackRejectTPCChi2;
    }

    // ITS quality
    if (track.itsNCls() < cfgTrackCuts.c1060itsNClsMin || track.itsNCls() > cfgTrackCuts.c1061itsNClsMax) {
      return kTrackRejectITSNCls;
    }

    if (track.itsNClsInnerBarrel() < cfgTrackCuts.c1062itsNClsInnerBarrelMin) {
      return kTrackRejectITSNClsInnerBarrel;
    }

    if (track.itsChi2NCl() < cfgTrackCuts.c1063itsChi2NClMin || track.itsChi2NCl() > cfgTrackCuts.c1064itsChi2NClMax) {
      return kTrackRejectITSChi2;
    }

    // pT-dependent DCA
    if (cfgTrackCuts.c1069usePtDependentDCAxy) {
      const float dcaXYMax = cfgTrackCuts.c1070dcaXYPtA + cfgTrackCuts.c1071dcaXYPtB / std::pow(track.pt(), cfgTrackCuts.c1072dcaXYPtC);
      if (std::abs(track.dcaXY()) > dcaXYMax) {
        return kTrackRejectPtDependentDCAxy;
      }
    }

    if (cfgTrackCuts.c1073usePtDependentDCAz) {
      const float dcaZMax = cfgTrackCuts.c1074dcaZPtA;
      if (std::abs(track.dcaZ()) > dcaZMax) {
        return kTrackRejectPtDependentDCAz;
      }
    }

    return kTrackAccepted;
  }

  template <int qaStage, typename C>
  void fillEventQA(const C& collision)
  {
    seEventQA.fill(HIST("SE/Events/") + HIST(QAStageDire[qaStage]) + HIST("VertexZ"), collision.posZ());

    seEventQA.fill(HIST("SE/Events/") + HIST(QAStageDire[qaStage]) + HIST("CentFT0C"), collision.centFT0C());
    seEventQA.fill(HIST("SE/Events/") + HIST(QAStageDire[qaStage]) + HIST("CentFT0M"), collision.centFT0M());
    seEventQA.fill(HIST("SE/Events/") + HIST(QAStageDire[qaStage]) + HIST("CentFT0A"), collision.centFT0A());
    seEventQA.fill(HIST("SE/Events/") + HIST(QAStageDire[qaStage]) + HIST("CentFV0A"), collision.centFV0A());

    seEventQA.fill(HIST("SE/Events/") + HIST(QAStageDire[qaStage]) + HIST("VtxZVsCentFT0C"), collision.posZ(), collision.centFT0C());
    seEventQA.fill(HIST("SE/Events/") + HIST(QAStageDire[qaStage]) + HIST("VtxZVsCentFT0M"), collision.posZ(), collision.centFT0M());
    seEventQA.fill(HIST("SE/Events/") + HIST(QAStageDire[qaStage]) + HIST("VtxZVsCentFT0A"), collision.posZ(), collision.centFT0A());
    seEventQA.fill(HIST("SE/Events/") + HIST(QAStageDire[qaStage]) + HIST("VtxZVsCentFV0A"), collision.posZ(), collision.centFV0A());

    seEventQA.fill(HIST("SE/Events/") + HIST(QAStageDire[qaStage]) + HIST("CentFT0CVsCentFT0M"), collision.centFT0C(), collision.centFT0M());
    seEventQA.fill(HIST("SE/Events/") + HIST(QAStageDire[qaStage]) + HIST("CentFT0CVsCentFT0A"), collision.centFT0C(), collision.centFT0A());
    seEventQA.fill(HIST("SE/Events/") + HIST(QAStageDire[qaStage]) + HIST("CentFT0MVsCentFT0A"), collision.centFT0M(), collision.centFT0A());

    seEventQA.fill(HIST("SE/Events/") + HIST(QAStageDire[qaStage]) + HIST("CentFT0CVsCentFV0A"), collision.centFT0C(), collision.centFV0A());
    seEventQA.fill(HIST("SE/Events/") + HIST(QAStageDire[qaStage]) + HIST("CentFT0MVsCentFV0A"), collision.centFT0M(), collision.centFV0A());
    seEventQA.fill(HIST("SE/Events/") + HIST(QAStageDire[qaStage]) + HIST("CentFT0AVsCentFV0A"), collision.centFT0A(), collision.centFV0A());

    const float mixingEstimatorValue = getMixingEstimatorValue(collision);
    const int mixingBin = getMixingBin(collision);
    const bool validMixingBin = mixingBin >= 0 && mixingBin < nMixBins;

    seEventQA.fill(HIST("SE/Events/") + HIST(QAStageDire[qaStage]) + HIST("MixingEstimator"), mixingEstimatorValue);
    seEventQA.fill(HIST("SE/Events/") + HIST(QAStageDire[qaStage]) + HIST("VtxZVsMixingEstimator"), collision.posZ(), mixingEstimatorValue);
    seEventQA.fill(HIST("SE/Events/") + HIST(QAStageDire[qaStage]) + HIST("MixingBinStatus"), validMixingBin ? 1.0f : 0.0f);

    if (validMixingBin) {
      seEventQA.fill(HIST("SE/Events/") + HIST(QAStageDire[qaStage]) + HIST("MixingBin"), mixingBin);
    }
  }

  template <int qaStage, int pidType, typename T>
  void fillFullTrackPIDQA(const T& track)
  {
    float tpcNSigma = 0.0f;
    float tofNSigma = 0.0f;

    if constexpr (pidType == kPi) {
      tpcNSigma = track.tpcNSigmaPi();
      tofNSigma = track.tofNSigmaPi();
    } else if constexpr (pidType == kKa) {
      tpcNSigma = track.tpcNSigmaKa();
      tofNSigma = track.tofNSigmaKa();
    } else if constexpr (pidType == kPr) {
      tpcNSigma = track.tpcNSigmaPr();
      tofNSigma = track.tofNSigmaPr();
    } else if constexpr (pidType == kEl) {
      tpcNSigma = track.tpcNSigmaEl();
      tofNSigma = track.tofNSigmaEl();
    } else if constexpr (pidType == kMu) {
      tpcNSigma = track.tpcNSigmaMu();
      tofNSigma = track.tofNSigmaMu();
    } else if constexpr (pidType == kDe) {
      tpcNSigma = track.tpcNSigmaDe();
      tofNSigma = track.tofNSigmaDe();
    }

    seTrackQA.fill(HIST("SE/Tracks/") + HIST(QAStageDire[qaStage]) + HIST("PID/") + HIST(PidTypeDire[pidType]) + HIST("TPCNSigmaVsP"), track.p(), tpcNSigma);
    seTrackQA.fill(HIST("SE/Tracks/") + HIST(QAStageDire[qaStage]) + HIST("PID/") + HIST(PidTypeDire[pidType]) + HIST("TOFNSigmaVsP"), track.p(), tofNSigma);
    seTrackQA.fill(HIST("SE/Tracks/") + HIST(QAStageDire[qaStage]) + HIST("PID/") + HIST(PidTypeDire[pidType]) + HIST("TPCNSigmaVsTOFNSigma"), tpcNSigma, tofNSigma);
  }

  template <int qaStage, typename T>
  void fillFullTrackQA(const T& track)
  {
    seTrackQA.fill(HIST("SE/Tracks/") + HIST(QAStageDire[qaStage]) + HIST("P"), track.p());
    seTrackQA.fill(HIST("SE/Tracks/") + HIST(QAStageDire[qaStage]) + HIST("Pt"), track.pt());
    seTrackQA.fill(HIST("SE/Tracks/") + HIST(QAStageDire[qaStage]) + HIST("TPCInnerParam"), track.tpcInnerParam());
    seTrackQA.fill(HIST("SE/Tracks/") + HIST(QAStageDire[qaStage]) + HIST("TOFExpMom"), track.tofExpMom());
    seTrackQA.fill(HIST("SE/Tracks/") + HIST(QAStageDire[qaStage]) + HIST("Eta"), track.eta());
    seTrackQA.fill(HIST("SE/Tracks/") + HIST(QAStageDire[qaStage]) + HIST("Phi"), track.phi());
    seTrackQA.fill(HIST("SE/Tracks/") + HIST(QAStageDire[qaStage]) + HIST("DCAxy"), track.dcaXY());
    seTrackQA.fill(HIST("SE/Tracks/") + HIST(QAStageDire[qaStage]) + HIST("DCAz"), track.dcaZ());
    seTrackQA.fill(HIST("SE/Tracks/") + HIST(QAStageDire[qaStage]) + HIST("Sign"), track.sign());

    seTrackQA.fill(HIST("SE/Tracks/") + HIST(QAStageDire[qaStage]) + HIST("DCAxyVsP"), track.p(), track.dcaXY());
    seTrackQA.fill(HIST("SE/Tracks/") + HIST(QAStageDire[qaStage]) + HIST("DCAxyVsPt"), track.pt(), track.dcaXY());
    seTrackQA.fill(HIST("SE/Tracks/") + HIST(QAStageDire[qaStage]) + HIST("DCAxyVsTPCInnerParam"), track.tpcInnerParam(), track.dcaXY());
    seTrackQA.fill(HIST("SE/Tracks/") + HIST(QAStageDire[qaStage]) + HIST("DCAxyVsTOFExpMom"), track.tofExpMom(), track.dcaXY());

    seTrackQA.fill(HIST("SE/Tracks/") + HIST(QAStageDire[qaStage]) + HIST("DCAzVsP"), track.p(), track.dcaZ());
    seTrackQA.fill(HIST("SE/Tracks/") + HIST(QAStageDire[qaStage]) + HIST("DCAzVsPt"), track.pt(), track.dcaZ());
    seTrackQA.fill(HIST("SE/Tracks/") + HIST(QAStageDire[qaStage]) + HIST("DCAzVsTPCInnerParam"), track.tpcInnerParam(), track.dcaZ());
    seTrackQA.fill(HIST("SE/Tracks/") + HIST(QAStageDire[qaStage]) + HIST("DCAzVsTOFExpMom"), track.tofExpMom(), track.dcaZ());

    seTrackQA.fill(HIST("SE/Tracks/") + HIST(QAStageDire[qaStage]) + HIST("PVsPt"), track.p(), track.pt());
    seTrackQA.fill(HIST("SE/Tracks/") + HIST(QAStageDire[qaStage]) + HIST("PVsTPCInnerParam"), track.p(), track.tpcInnerParam());
    seTrackQA.fill(HIST("SE/Tracks/") + HIST(QAStageDire[qaStage]) + HIST("PVsTOFExpMom"), track.p(), track.tofExpMom());

    seTrackQA.fill(HIST("SE/Tracks/") + HIST(QAStageDire[qaStage]) + HIST("TPCSignalVsP"), track.p(), track.tpcSignal());
    seTrackQA.fill(HIST("SE/Tracks/") + HIST(QAStageDire[qaStage]) + HIST("TPCSignalVsTPCInnerParam"), track.tpcInnerParam(), track.tpcSignal());
    seTrackQA.fill(HIST("SE/Tracks/") + HIST(QAStageDire[qaStage]) + HIST("TPCSignalVsTOFExpMom"), track.tofExpMom(), track.tpcSignal());

    seTrackQA.fill(HIST("SE/Tracks/") + HIST(QAStageDire[qaStage]) + HIST("TOFBetaVsP"), track.p(), track.beta());
    seTrackQA.fill(HIST("SE/Tracks/") + HIST(QAStageDire[qaStage]) + HIST("TOFBetaVsTPCInnerParam"), track.tpcInnerParam(), track.beta());
    seTrackQA.fill(HIST("SE/Tracks/") + HIST(QAStageDire[qaStage]) + HIST("TOFBetaVsTOFExpMom"), track.tofExpMom(), track.beta());

    fillFullTrackPIDQA<qaStage, kPi>(track);
    fillFullTrackPIDQA<qaStage, kKa>(track);
    fillFullTrackPIDQA<qaStage, kPr>(track);
    fillFullTrackPIDQA<qaStage, kEl>(track);
    fillFullTrackPIDQA<qaStage, kMu>(track);
    fillFullTrackPIDQA<qaStage, kDe>(track);
  }

  template <int eventType, int pairType, int analysisType, int roleType, typename H, typename T>
  void fillCorrParticleQA(H& histReg, const T& track)
  {
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[analysisType]) + HIST(CorrRoleDire[roleType]) + HIST("P"), track.p());
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[analysisType]) + HIST(CorrRoleDire[roleType]) + HIST("Pt"), track.pt());
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[analysisType]) + HIST(CorrRoleDire[roleType]) + HIST("Eta"), track.eta());
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[analysisType]) + HIST(CorrRoleDire[roleType]) + HIST("Phi"), track.phi());
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[analysisType]) + HIST(CorrRoleDire[roleType]) + HIST("DCAxy"), track.dcaXY());
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[analysisType]) + HIST(CorrRoleDire[roleType]) + HIST("DCAz"), track.dcaZ());
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[analysisType]) + HIST(CorrRoleDire[roleType]) + HIST("Sign"), track.sign());

    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[analysisType]) + HIST(CorrRoleDire[roleType]) + HIST("DCAxyVsPt"), track.pt(), track.dcaXY());
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[analysisType]) + HIST(CorrRoleDire[roleType]) + HIST("DCAzVsPt"), track.pt(), track.dcaZ());

    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[analysisType]) + HIST(CorrRoleDire[roleType]) + HIST("TPCSignalVsTPCInnerParam"), track.tpcInnerParam(), track.tpcSignal());
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[analysisType]) + HIST(CorrRoleDire[roleType]) + HIST("TPCSignalVsP"), track.p(), track.tpcSignal());

    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[analysisType]) + HIST(CorrRoleDire[roleType]) + HIST("TOFBetaVsTOFExpMom"), track.tofExpMom(), track.beta());
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[analysisType]) + HIST(CorrRoleDire[roleType]) + HIST("TOFBetaVsP"), track.p(), track.beta());
  }

  template <int eventType, int pidType, int roleType, typename H, typename T, typename A>
  void fillIdentifiedCorrelation(H& histReg, const T& trigger, const A& associate, const int idMethod)
  {
    const float dPhi = computeDeltaPhi(associate.phi(), trigger.phi());
    const float dEta = associate.eta() - trigger.eta();

    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[kHId]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(PidTypeDire[pidType]) + HIST("P"), associate.p());
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[kHId]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(PidTypeDire[pidType]) + HIST("Pt"), associate.pt());
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[kHId]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(PidTypeDire[pidType]) + HIST("Eta"), associate.eta());
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[kHId]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(PidTypeDire[pidType]) + HIST("Phi"), associate.phi());
    if constexpr (pidType == kPi) {
      histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[kHId]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(PidTypeDire[pidType]) + HIST("Rapidity"), computeRapidity(associate, MassPiPlus));
    } else if constexpr (pidType == kKa) {
      histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[kHId]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(PidTypeDire[pidType]) + HIST("Rapidity"), computeRapidity(associate, MassKPlus));
    } else if constexpr (pidType == kPr) {
      histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[kHId]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(PidTypeDire[pidType]) + HIST("Rapidity"), computeRapidity(associate, MassProton));
    }
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[kHId]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(PidTypeDire[pidType]) + HIST("DCAxy"), associate.dcaXY());
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[kHId]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(PidTypeDire[pidType]) + HIST("DCAz"), associate.dcaZ());
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[kHId]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(PidTypeDire[pidType]) + HIST("Sign"), associate.sign());

    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[kHId]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(PidTypeDire[pidType]) + HIST("DCAxyVsPt"), associate.pt(), associate.dcaXY());
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[kHId]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(PidTypeDire[pidType]) + HIST("DCAzVsPt"), associate.pt(), associate.dcaZ());

    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[kHId]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(PidTypeDire[pidType]) + HIST("TPCSignalVsP"), associate.p(), associate.tpcSignal());
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[kHId]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(PidTypeDire[pidType]) + HIST("TOFBetaVsP"), associate.p(), associate.beta());
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[kHId]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(PidTypeDire[pidType]) + HIST("IdMethodVsP"), associate.p(), idMethod);

    if constexpr (pidType == kPi) {
      histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[kHId]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(PidTypeDire[pidType]) + HIST("TPCNSigmaVsP"), associate.p(), associate.tpcNSigmaPi());
      histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[kHId]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(PidTypeDire[pidType]) + HIST("TOFNSigmaVsP"), associate.p(), associate.tofNSigmaPi());
      histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[kHId]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(PidTypeDire[pidType]) + HIST("TPCNSigmaVsTOFNSigma"), associate.tpcNSigmaPi(), associate.tofNSigmaPi());
    } else if constexpr (pidType == kKa) {
      histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[kHId]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(PidTypeDire[pidType]) + HIST("TPCNSigmaVsP"), associate.p(), associate.tpcNSigmaKa());
      histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[kHId]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(PidTypeDire[pidType]) + HIST("TOFNSigmaVsP"), associate.p(), associate.tofNSigmaKa());
      histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[kHId]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(PidTypeDire[pidType]) + HIST("TPCNSigmaVsTOFNSigma"), associate.tpcNSigmaKa(), associate.tofNSigmaKa());
    } else if constexpr (pidType == kPr) {
      histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[kHId]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(PidTypeDire[pidType]) + HIST("TPCNSigmaVsP"), associate.p(), associate.tpcNSigmaPr());
      histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[kHId]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(PidTypeDire[pidType]) + HIST("TOFNSigmaVsP"), associate.p(), associate.tofNSigmaPr());
      histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[kHId]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(PidTypeDire[pidType]) + HIST("TPCNSigmaVsTOFNSigma"), associate.tpcNSigmaPr(), associate.tofNSigmaPr());
    }

    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[kHId]) + HIST(AnalysisTypeDire[kCorr]) + HIST(CorrRoleDire[roleType]) + HIST(PidTypeDire[pidType]) + HIST("DeltaPhi"), dPhi);
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[kHId]) + HIST(AnalysisTypeDire[kCorr]) + HIST(CorrRoleDire[roleType]) + HIST(PidTypeDire[pidType]) + HIST("DeltaEta"), dEta);
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[kHId]) + HIST(AnalysisTypeDire[kCorr]) + HIST(CorrRoleDire[roleType]) + HIST(PidTypeDire[pidType]) + HIST("DeltaPhiDeltaEta"), dPhi, dEta, currentCentrality);
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[kHId]) + HIST(AnalysisTypeDire[kCorr]) + HIST(CorrRoleDire[roleType]) + HIST(PidTypeDire[pidType]) + HIST("CorrSparse"), trigger.pt(), associate.pt(), dPhi, dEta, currentCentrality);
  }

  //______________________________________________________________________________
  // Resonance candidate QA
  template <int eventType, int pairType, int roleType, int massRegion, typename H, typename A>
  void fillResoAssocQA(H& histReg, const A& associate, const float invMass)
  {
    const float p = std::sqrt(associate.px() * associate.px() + associate.py() * associate.py() + associate.pz() * associate.pz());
    const float energy = std::sqrt(p * p + invMass * invMass);
    const float rapidity = 0.5f * std::log((energy + associate.pz()) / (energy - associate.pz()));

    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(MassRegionDire[massRegion]) + HIST("Reso/") + HIST("P"), p);
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(MassRegionDire[massRegion]) + HIST("Reso/") + HIST("Pt"), associate.pt());
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(MassRegionDire[massRegion]) + HIST("Reso/") + HIST("Eta"), associate.eta());
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(MassRegionDire[massRegion]) + HIST("Reso/") + HIST("Phi"), associate.phi());
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(MassRegionDire[massRegion]) + HIST("Reso/") + HIST("Rapidity"), rapidity);
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(MassRegionDire[massRegion]) + HIST("Reso/") + HIST("InvMass"), invMass);
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(MassRegionDire[massRegion]) + HIST("Reso/") + HIST("InvMassVsPt"), associate.pt(), invMass);
  }

  //______________________________________________________________________________
  // Resonance daughter QA
  template <int eventType, int pairType, int roleType, int massRegion, int dauType, typename H, typename T>
  void fillResoDauQA(H& histReg, const T& track, const int idMethod)
  {
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(MassRegionDire[massRegion]) + HIST(DauTypeDire[dauType]) + HIST("P"), track.p());
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(MassRegionDire[massRegion]) + HIST(DauTypeDire[dauType]) + HIST("Pt"), track.pt());
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(MassRegionDire[massRegion]) + HIST(DauTypeDire[dauType]) + HIST("Eta"), track.eta());
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(MassRegionDire[massRegion]) + HIST(DauTypeDire[dauType]) + HIST("Phi"), track.phi());
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(MassRegionDire[massRegion]) + HIST(DauTypeDire[dauType]) + HIST("DCAxy"), track.dcaXY());
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(MassRegionDire[massRegion]) + HIST(DauTypeDire[dauType]) + HIST("DCAz"), track.dcaZ());
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(MassRegionDire[massRegion]) + HIST(DauTypeDire[dauType]) + HIST("Sign"), track.sign());

    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(MassRegionDire[massRegion]) + HIST(DauTypeDire[dauType]) + HIST("DCAxyVsPt"), track.pt(), track.dcaXY());
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(MassRegionDire[massRegion]) + HIST(DauTypeDire[dauType]) + HIST("DCAzVsPt"), track.pt(), track.dcaZ());

    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(MassRegionDire[massRegion]) + HIST(DauTypeDire[dauType]) + HIST("TPCSignalVsP"), track.p(), track.tpcSignal());
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(MassRegionDire[massRegion]) + HIST(DauTypeDire[dauType]) + HIST("TOFBetaVsP"), track.p(), track.beta());
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(MassRegionDire[massRegion]) + HIST(DauTypeDire[dauType]) + HIST("IdMethodVsP"), track.p(), idMethod);

    float tpcNSigma = 0.0f;
    float tofNSigma = 0.0f;
    float massHypothesis = 0.0f;

    if constexpr (pairType == kHPhi || (pairType == kHKStar && dauType == kPosDau) || (pairType == kHKStarBar && dauType == kNegDau) ||
                  (pairType == kHLambda && dauType == kNegDau) || (pairType == kHLambdaBar && dauType == kPosDau)) {
      tpcNSigma = track.tpcNSigmaKa();
      tofNSigma = track.tofNSigmaKa();
      massHypothesis = MassKPlus;
    } else if constexpr ((pairType == kHKStar && dauType == kNegDau) || (pairType == kHKStarBar && dauType == kPosDau)) {
      tpcNSigma = track.tpcNSigmaPi();
      tofNSigma = track.tofNSigmaPi();
      massHypothesis = MassPiPlus;
    } else if constexpr ((pairType == kHLambda && dauType == kPosDau) || (pairType == kHLambdaBar && dauType == kNegDau)) {
      tpcNSigma = track.tpcNSigmaPr();
      tofNSigma = track.tofNSigmaPr();
      massHypothesis = MassProton;
    }

    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(MassRegionDire[massRegion]) + HIST(DauTypeDire[dauType]) + HIST("Rapidity"), computeRapidity(track, massHypothesis));
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(MassRegionDire[massRegion]) + HIST(DauTypeDire[dauType]) + HIST("TPCNSigmaVsP"), track.p(), tpcNSigma);
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(MassRegionDire[massRegion]) + HIST(DauTypeDire[dauType]) + HIST("TOFNSigmaVsP"), track.p(), tofNSigma);
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(MassRegionDire[massRegion]) + HIST(DauTypeDire[dauType]) + HIST("TPCNSigmaVsTOFNSigma"), tpcNSigma, tofNSigma);
  }

  //______________________________________________________________________________
  // Generic h-h correlation
  template <int eventType, int pairType, int analysisType, int roleType, typename H, typename T, typename A>
  void fillTrigAssociateCorrelation(H& histReg, const T& trigger, const A& associate)
  {
    const float dPhi = computeDeltaPhi(associate.phi(), trigger.phi());
    const float dEta = associate.eta() - trigger.eta();

    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[analysisType]) + HIST(CorrRoleDire[roleType]) + HIST("DeltaPhi"), dPhi);
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[analysisType]) + HIST(CorrRoleDire[roleType]) + HIST("DeltaEta"), dEta);
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[analysisType]) + HIST(CorrRoleDire[roleType]) + HIST("DeltaPhiDeltaEta"), dPhi, dEta, currentCentrality);
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[analysisType]) + HIST(CorrRoleDire[roleType]) + HIST("CorrSparse"), trigger.pt(), associate.pt(), dPhi, dEta, currentCentrality);
  }

  //______________________________________________________________________________
  // Resonance correlation
  template <int eventType, int pairType, int roleType, int massRegion, typename H, typename T, typename A>
  void fillResoCorrelation(H& histReg, const T& trigger, const A& associate)
  {
    const float dPhi = computeDeltaPhi(associate.phi(), trigger.phi());
    const float dEta = associate.eta() - trigger.eta();

    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[kCorr]) + HIST(CorrRoleDire[roleType]) + HIST(MassRegionDire[massRegion]) + HIST("DeltaPhi"), dPhi);
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[kCorr]) + HIST(CorrRoleDire[roleType]) + HIST(MassRegionDire[massRegion]) + HIST("DeltaEta"), dEta);
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[kCorr]) + HIST(CorrRoleDire[roleType]) + HIST(MassRegionDire[massRegion]) + HIST("DeltaPhiDeltaEta"), dPhi, dEta, currentCentrality);
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[kCorr]) + HIST(CorrRoleDire[roleType]) + HIST(MassRegionDire[massRegion]) + HIST("CorrSparse"), trigger.pt(), associate.pt(), dPhi, dEta, currentCentrality);
  }

  // Resonance correlation as a function of invariant mass
  template <int eventType, int pairType, int roleType, typename H, typename T, typename A>
  void fillResoMassCorrelation(H& histReg, const T& trigger, const A& associate, const float invMass)
  {
    const float dPhi = computeDeltaPhi(associate.phi(), trigger.phi());
    const float dEta = associate.eta() - trigger.eta();

    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[kCorr]) + HIST(CorrRoleDire[roleType]) + HIST("MassDependent/") + HIST("DeltaPhiVsInvMass"), invMass, dPhi, currentCentrality);
    histReg.fill(HIST(EventTypeDire[eventType]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[kCorr]) + HIST(CorrRoleDire[roleType]) + HIST("MassDependent/") + HIST("DeltaEtaVsInvMass"), invMass, dEta, currentCentrality);
  }

  template <int eventType, int pairType, int roleType, int massRegion, typename H, typename T, typename A, typename P, typename N>
  void fillResoRegion(H& histReg, const T& trigger, const A& associate, const P& posDauTrack, const N& negDauTrack, const float invMass, const int posDauIdMethod, const int negDauIdMethod)
  {
    fillResoAssocQA<eventType, pairType, roleType, massRegion>(histReg, associate, invMass);
    fillResoDauQA<eventType, pairType, roleType, massRegion, kPosDau>(histReg, posDauTrack, posDauIdMethod);
    fillResoDauQA<eventType, pairType, roleType, massRegion, kNegDau>(histReg, negDauTrack, negDauIdMethod);
    fillResoCorrelation<eventType, pairType, roleType, massRegion>(histReg, trigger, associate);
  }

  template <int pairType, typename P, typename N>
  bool selectResonanceDaughters(const P& posDauTrack, const N& negDauTrack, int& posDauIdMethod, int& negDauIdMethod)
  {
    // Phi        -> K+ K-       -----------------------
    // K(892)*    -> K+ + pi-    -----------------------
    // K(892)*Bar -> K- + pi+    -----------------------
    // Λ(1520)    -> P+ + Ka-    -----------------------
    // Λ(1520)Bar -> PBar + Ka+  -----------------------

    posDauIdMethod = kUnidentified;
    negDauIdMethod = kUnidentified;
    if constexpr (pairType == kHPhi) {
      return selKaon(posDauTrack, posDauIdMethod) && selKaon(negDauTrack, negDauIdMethod);
    } else if constexpr (pairType == kHKStar) {
      return selKaon(posDauTrack, posDauIdMethod) && selPion(negDauTrack, negDauIdMethod);
    } else if constexpr (pairType == kHKStarBar) {
      return selPion(posDauTrack, posDauIdMethod) && selKaon(negDauTrack, negDauIdMethod);
    } else if constexpr (pairType == kHLambda) {
      return selProton(posDauTrack, posDauIdMethod) && selKaon(negDauTrack, negDauIdMethod);
    } else if constexpr (pairType == kHLambdaBar) {
      return selKaon(posDauTrack, posDauIdMethod) && selProton(negDauTrack, negDauIdMethod);
    }
    return false;
  }

  template <int pairType, typename A>
  uint8_t getResonanceMassRegion(const A& associate, float& invMass)
  {
    if constexpr (pairType == kHPhi) {
      invMass = associate.mPhi1020();
      return getMassRegionTag(invMass, cfgPhi1020CorrMass.c1117phi1020LSBLow, cfgPhi1020CorrMass.c1118phi1020LSBUp, cfgPhi1020CorrMass.c1115phi1020PeakLow, cfgPhi1020CorrMass.c1116phi1020PeakUp, cfgPhi1020CorrMass.c1119phi1020RSBLow, cfgPhi1020CorrMass.c1120phi1020RSBUp);
    } else if constexpr (pairType == kHKStar) {
      invMass = associate.mKStar892();
      return getMassRegionTag(invMass, cfgKStar892CorrMass.c1123kStar892LSBLow, cfgKStar892CorrMass.c1124kStar892LSBUp, cfgKStar892CorrMass.c1121kStar892PeakLow, cfgKStar892CorrMass.c1122kStar892PeakUp, cfgKStar892CorrMass.c1125kStar892RSBLow, cfgKStar892CorrMass.c1126kStar892RSBUp);
    } else if constexpr (pairType == kHKStarBar) {
      invMass = associate.mKStar892Bar();
      return getMassRegionTag(invMass, cfgKStar892CorrMass.c1123kStar892LSBLow, cfgKStar892CorrMass.c1124kStar892LSBUp, cfgKStar892CorrMass.c1121kStar892PeakLow, cfgKStar892CorrMass.c1122kStar892PeakUp, cfgKStar892CorrMass.c1125kStar892RSBLow, cfgKStar892CorrMass.c1126kStar892RSBUp);
    } else if constexpr (pairType == kHLambda) {
      invMass = associate.mLambda1520();
      return getMassRegionTag(invMass, cfgLambda1520CorrMass.c1129lambda1520LSBLow, cfgLambda1520CorrMass.c1130lambda1520LSBUp, cfgLambda1520CorrMass.c1127lambda1520PeakLow, cfgLambda1520CorrMass.c1128lambda1520PeakUp, cfgLambda1520CorrMass.c1131lambda1520RSBLow, cfgLambda1520CorrMass.c1132lambda1520RSBUp);
    } else if constexpr (pairType == kHLambdaBar) {
      invMass = associate.mLambda1520Bar();
      return getMassRegionTag(invMass, cfgLambda1520CorrMass.c1129lambda1520LSBLow, cfgLambda1520CorrMass.c1130lambda1520LSBUp, cfgLambda1520CorrMass.c1127lambda1520PeakLow, cfgLambda1520CorrMass.c1128lambda1520PeakUp, cfgLambda1520CorrMass.c1131lambda1520RSBLow, cfgLambda1520CorrMass.c1132lambda1520RSBUp);
    }

    invMass = -1.0f;
    return kMassOutside;
  }

  template <int pairType, int roleType, int pidType = -1>
  void incrementCorrCount(CorrCountArray& counts)
  {
    constexpr int PtType = roleType == kAssocLowPt ? kCountLowPt : kCountHighPt;

    if constexpr (pairType == kHH) {
      ++counts[kCountHH][PtType];
    } else if constexpr (pairType == kHId && pidType == kPi) {
      ++counts[kCountHPi][PtType];
    } else if constexpr (pairType == kHId && pidType == kKa) {
      ++counts[kCountHKa][PtType];
    } else if constexpr (pairType == kHId && pidType == kPr) {
      ++counts[kCountHPr][PtType];
    } else if constexpr (pairType == kHPhi) {
      ++counts[kCountHPhi][PtType];
    } else if constexpr (pairType == kHKStar) {
      ++counts[kCountHKStar][PtType];
    } else if constexpr (pairType == kHKStarBar) {
      ++counts[kCountHKStarBar][PtType];
    } else if constexpr (pairType == kHLambda) {
      ++counts[kCountHLambda][PtType];
    } else if constexpr (pairType == kHLambdaBar) {
      ++counts[kCountHLambdaBar][PtType];
    }
  }

  template <int pairType, int roleType>
  void incrementResoRegionCount(ResoRegionCountArray& counts, uint8_t massRegion)
  {
    constexpr int PtType = roleType == kAssocLowPt ? kCountLowPt : kCountHighPt;

    if (massRegion != kMassPeak && massRegion != kMassLSB && massRegion != kMassRSB) {
      return;
    }

    if constexpr (pairType == kHPhi) {
      ++counts[kResoCountHPhi][PtType][massRegion];
    } else if constexpr (pairType == kHKStar) {
      ++counts[kResoCountHKStar][PtType][massRegion];
    } else if constexpr (pairType == kHKStarBar) {
      ++counts[kResoCountHKStarBar][PtType][massRegion];
    } else if constexpr (pairType == kHLambda) {
      ++counts[kResoCountHLambda][PtType][massRegion];
    } else if constexpr (pairType == kHLambdaBar) {
      ++counts[kResoCountHLambdaBar][PtType][massRegion];
    }
  }

  template <int pairType, typename T1, typename T2>
  float getInvariantMass(const T1& dau1, const T2& dau2)
  {
    const float p = RecoDecay::p(dau1.px() + dau2.px(), dau1.py() + dau2.py(), dau1.pz() + dau2.pz());

    if constexpr (pairType == kHPhi) {
      return RecoDecay::m(p, RecoDecay::e(dau1.px(), dau1.py(), dau1.pz(), MassKPlus) + RecoDecay::e(dau2.px(), dau2.py(), dau2.pz(), MassKPlus));
    } else if constexpr (pairType == kHKStar) {
      return RecoDecay::m(p, RecoDecay::e(dau1.px(), dau1.py(), dau1.pz(), MassKPlus) + RecoDecay::e(dau2.px(), dau2.py(), dau2.pz(), MassPiPlus));
    } else if constexpr (pairType == kHKStarBar) {
      return RecoDecay::m(p, RecoDecay::e(dau1.px(), dau1.py(), dau1.pz(), MassPiPlus) + RecoDecay::e(dau2.px(), dau2.py(), dau2.pz(), MassKPlus));
    } else if constexpr (pairType == kHLambda) {
      return RecoDecay::m(p, RecoDecay::e(dau1.px(), dau1.py(), dau1.pz(), MassProton) + RecoDecay::e(dau2.px(), dau2.py(), dau2.pz(), MassKPlus));
    } else if constexpr (pairType == kHLambdaBar) {
      return RecoDecay::m(p, RecoDecay::e(dau1.px(), dau1.py(), dau1.pz(), MassKPlus) + RecoDecay::e(dau2.px(), dau2.py(), dau2.pz(), MassProton));
    }

    return -1.0f;
  }

  template <int pairType, bool isLikeSign, typename C, typename T1, typename T2>
  bool fillInvariantMassPair(const C& collision, const T1& dau1, const T2& dau2)
  {
    int dau1IdMethod = kUnidentified;
    int dau2IdMethod = kUnidentified;

    if (!selectResonanceDaughters<pairType>(dau1, dau2, dau1IdMethod, dau2IdMethod)) {
      return false;
    }

    const float mass = getInvariantMass<pairType>(dau1, dau2);
    const float pt = RecoDecay::pt(dau1.px() + dau2.px(), dau1.py() + dau2.py());
    const float pz = dau1.pz() + dau2.pz();
    const float energy = RecoDecay::e(dau1.px() + dau2.px(), dau1.py() + dau2.py(), pz, mass);
    const float rapidity = 0.5f * std::log((energy + pz) / (energy - pz));
    if (rapidity < cfgResPartitions.c1083resoRapidityMin || rapidity > cfgResPartitions.c1084resoRapidityMax) {
      return false;
    }

    float centrality = -1.0f;
    switch (cfgEvent.c1028centralityEstimator) {
      case 0:
        centrality = collision.centFT0C();
        break;
      case 1:
        centrality = collision.centFT0M();
        break;
      case 2:
        centrality = collision.centFT0A();
        break;
      case 3:
        centrality = collision.centFV0A();
        break;
      default:
        return false;
    }

    if constexpr (pairType == kHPhi && !isLikeSign) {
      resonanceMassSpectra.fill(HIST("SE/InvariantMass/Phi/US"), mass, pt, centrality);
    } else if constexpr (pairType == kHPhi && isLikeSign) {
      resonanceMassSpectra.fill(HIST("SE/InvariantMass/Phi/LS"), mass, pt, centrality);
    } else if constexpr (pairType == kHKStar && !isLikeSign) {
      resonanceMassSpectra.fill(HIST("SE/InvariantMass/KStar/US"), mass, pt, centrality);
    } else if constexpr (pairType == kHKStar && isLikeSign) {
      resonanceMassSpectra.fill(HIST("SE/InvariantMass/KStar/LS"), mass, pt, centrality);
    } else if constexpr (pairType == kHKStarBar && !isLikeSign) {
      resonanceMassSpectra.fill(HIST("SE/InvariantMass/KStarBar/US"), mass, pt, centrality);
    } else if constexpr (pairType == kHKStarBar && isLikeSign) {
      resonanceMassSpectra.fill(HIST("SE/InvariantMass/KStarBar/LS"), mass, pt, centrality);
    } else if constexpr (pairType == kHLambda && !isLikeSign) {
      resonanceMassSpectra.fill(HIST("SE/InvariantMass/Lambda1520/US"), mass, pt, centrality);
    } else if constexpr (pairType == kHLambda && isLikeSign) {
      resonanceMassSpectra.fill(HIST("SE/InvariantMass/Lambda1520/LS"), mass, pt, centrality);
    } else if constexpr (pairType == kHLambdaBar && !isLikeSign) {
      resonanceMassSpectra.fill(HIST("SE/InvariantMass/Lambda1520Bar/US"), mass, pt, centrality);
    } else if constexpr (pairType == kHLambdaBar && isLikeSign) {
      resonanceMassSpectra.fill(HIST("SE/InvariantMass/Lambda1520Bar/LS"), mass, pt, centrality);
    }

    return true;
  }

  template <typename C, typename P, typename N>
  void fillUnlikeSignInvariantMass(const C& collision, const P& posTracksPerColl, const N& negTracksPerColl)
  {
    for (const auto& posTrack : posTracksPerColl) {
      if (selectionTrack(posTrack) != kTrackAccepted) {
        continue;
      }

      for (const auto& negTrack : negTracksPerColl) {
        if (selectionTrack(negTrack) != kTrackAccepted) {
          continue;
        }

        fillInvariantMassPair<kHPhi, false>(collision, posTrack, negTrack);
        fillInvariantMassPair<kHKStar, false>(collision, posTrack, negTrack);
        fillInvariantMassPair<kHKStarBar, false>(collision, posTrack, negTrack);
        fillInvariantMassPair<kHLambda, false>(collision, posTrack, negTrack);
        fillInvariantMassPair<kHLambdaBar, false>(collision, posTrack, negTrack);
      }
    }
  }

  template <bool positiveSign, typename C, typename T>
  void fillLikeSignInvariantMass(const C& collision, const T& tracksPerColl)
  {
    for (const auto& track1 : tracksPerColl) {
      if (selectionTrack(track1) != kTrackAccepted) {
        continue;
      }

      for (const auto& track2 : tracksPerColl) {
        if (track2.globalIndex() <= track1.globalIndex()) {
          continue;
        }
        if (selectionTrack(track2) != kTrackAccepted) {
          continue;
        }

        fillInvariantMassPair<kHPhi, true>(collision, track1, track2);

        if constexpr (positiveSign) {
          if (!fillInvariantMassPair<kHKStar, true>(collision, track1, track2)) {
            fillInvariantMassPair<kHKStar, true>(collision, track2, track1);
          }
          if (!fillInvariantMassPair<kHLambda, true>(collision, track1, track2)) {
            fillInvariantMassPair<kHLambda, true>(collision, track2, track1);
          }
        } else {
          if (!fillInvariantMassPair<kHKStarBar, true>(collision, track1, track2)) {
            fillInvariantMassPair<kHKStarBar, true>(collision, track2, track1);
          }
          if (!fillInvariantMassPair<kHLambdaBar, true>(collision, track1, track2)) {
            fillInvariantMassPair<kHLambdaBar, true>(collision, track2, track1);
          }
        }
      }
    }
  }

  uint64_t getStoredCorrMask(uint64_t corrMask, int mixingBin)
  {
    if (mixingBin < 0 || mixingBin >= nMixBins) {
      return 0ULL;
    }

    if (cfgDerivedData.c1134mixingBin >= 0 && mixingBin != cfgDerivedData.c1134mixingBin) {
      return 0ULL;
    }

    const uint64_t requestedMask = cfgDerivedData.c1133pairMask;
    if (requestedMask == 0ULL) {
      return corrMask;
    }

    if (cfgDerivedData.c1135requireAllPairBits && (corrMask & requestedMask) != requestedMask) {
      return 0ULL;
    }

    // Keeping only the requested pools, so tracks / resonances of other channels are not written.
    // Resonance pools are the Peak/LSB/RSB bits (kMaskHPhiLowPtPeak ... kMaskHLambdaBarHighPtRSB).
    return corrMask & requestedMask;
  }

  //______________________________________________________________________________
  //   Tracks     : bit 0-7 = pools h, pi, K, p x low/high pT (same bit numbers as CorrPresenceMask),
  //                bit 15  = used as trigger
  //   Resonances : bit r   = CorrPresenceMask pool bit (18 + r), r = species * 6 + pT class * 3 + (mass window - 1)
  // Setting a bit again (e.g. for another trigger) changes nothing, so there are no duplicates,
  // and no sorting or searching is needed. Re-sized per collision; the vectors keep their capacity.
  static constexpr uint16_t TriggerUsedFlag = 1 << 15;
  std::vector<uint16_t> trackUsedFlags;
  std::vector<uint32_t> resoUsedFlags;
  int64_t firstTrackGI = 0;
  int64_t firstResoGI = 0;
  float currentCentrality = -1.0f; // mixing-estimator value of the current collision, used as centrality axis of the correlation histograms

  template <typename T>
  void markTrackUsed(const T& track, const uint16_t flag)
  {
    const int64_t i = track.globalIndex() - firstTrackGI;
    if (i >= 0 && i < static_cast<int64_t>(trackUsedFlags.size())) {
      trackUsedFlags[i] |= flag;
    }
  }

  template <typename R>
  void markResoUsed(const R& resonance, const uint32_t flag)
  {
    const int64_t i = resonance.globalIndex() - firstResoGI;
    if (i >= 0 && i < static_cast<int64_t>(resoUsedFlags.size())) {
      resoUsedFlags[i] |= flag;
    }
  }

  // Particle (+1) / antiparticle (-1) / self-conjugate (0) flag from the associate type
  static int8_t getResoSign(const int assocType)
  {
    if (assocType == kCountHKStar || assocType == kCountHLambda) {
      return 1;
    }
    if (assocType == kCountHKStarBar || assocType == kCountHLambdaBar) {
      return -1;
    }
    return 0; // phi(1020)
  }

  // 01-Start-obtaining trigger-LowPtAssociates correlation
  template <typename FilterTrkType, int eventType, int pairType, int roleType, typename H, typename T, typename A>
  void executeTrigAssocCorrelation(H& histReg, const T& trigger, const A& associates, const auto& triggerGIList, CorrCountArray& corrCountPerColl, ResoRegionCountArray& resoRegionCountsPerColl)
  {
    constexpr int PtType = roleType == kAssocLowPt ? kCountLowPt : kCountHighPt;
    int associateSelTag = kTrackAccepted;
    int assocPosDauSelTag = kTrackAccepted;
    int assocNegDauSelTag = kTrackAccepted;
    int posDauIdMethod = kUnidentified;
    int negDauIdMethod = kUnidentified;
    int pionIdMethod = kUnidentified;
    int kaonIdMethod = kUnidentified;
    int protonIdMethod = kUnidentified;

    for (const auto& associate : associates) {
      //__________________________________________________________________________
      // Hadron-hadron and hadron-Identified correlation
      if constexpr (pairType == kHH || pairType == kHId) {
        associateSelTag = selectionTrack(associate);
        if (associateSelTag != kTrackAccepted) {
          continue;
        }

        if constexpr (pairType == kHH) {
          fillCorrParticleQA<eventType, pairType, kQA, roleType>(histReg, associate);
          fillTrigAssociateCorrelation<eventType, pairType, kCorr, roleType>(histReg, trigger, associate);
          incrementCorrCount<pairType, roleType>(corrCountPerColl);
          markTrackUsed(associate, static_cast<uint16_t>(getCorrMaskBit(kCountHH, PtType)));
        }

        if constexpr (pairType == kHId) {
          pionIdMethod = kUnidentified;
          kaonIdMethod = kUnidentified;
          protonIdMethod = kUnidentified;

          if (selPion(associate, pionIdMethod)) {
            fillIdentifiedCorrelation<eventType, kPi, roleType>(histReg, trigger, associate, pionIdMethod);
            incrementCorrCount<pairType, roleType, kPi>(corrCountPerColl);
            markTrackUsed(associate, static_cast<uint16_t>(getCorrMaskBit(kCountHPi, PtType)));
          }
          if (selKaon(associate, kaonIdMethod)) {
            fillIdentifiedCorrelation<eventType, kKa, roleType>(histReg, trigger, associate, kaonIdMethod);
            incrementCorrCount<pairType, roleType, kKa>(corrCountPerColl);
            markTrackUsed(associate, static_cast<uint16_t>(getCorrMaskBit(kCountHKa, PtType)));
          }
          if (selProton(associate, protonIdMethod)) {
            fillIdentifiedCorrelation<eventType, kPr, roleType>(histReg, trigger, associate, protonIdMethod);
            incrementCorrCount<pairType, roleType, kPr>(corrCountPerColl);
            markTrackUsed(associate, static_cast<uint16_t>(getCorrMaskBit(kCountHPr, PtType)));
          }
        }
      }

      //__________________________________________________________________________
      // Hadron-Resonance correlation
      if constexpr (pairType == kHPhi || pairType == kHKStar || pairType == kHKStarBar || pairType == kHLambda || pairType == kHLambdaBar) {
        const auto& posDauTrack = associate.template posTrack_as<FilterTrkType>();
        const auto& negDauTrack = associate.template negTrack_as<FilterTrkType>();

        // Track Selection ------------------------
        assocPosDauSelTag = selectionTrack(posDauTrack);
        if (assocPosDauSelTag != kTrackAccepted) {
          continue;
        }

        assocNegDauSelTag = selectionTrack(negDauTrack);
        if (assocNegDauSelTag != kTrackAccepted) {
          continue;
        }

        if (!selectResonanceDaughters<pairType>(posDauTrack, negDauTrack, posDauIdMethod, negDauIdMethod)) {
          continue;
        }

        float invMass = -1.0f;
        const uint8_t massRegion = getResonanceMassRegion<pairType>(associate, invMass);
        const float rapidity = computeRapidity(associate, invMass);
        if (rapidity < cfgResPartitions.c1083resoRapidityMin || rapidity > cfgResPartitions.c1084resoRapidityMax) {
          continue;
        }

        // Trigger-resonance daughter overlap -----expensive check ----
        if (c1113rejectResoWithAnyTriggerDaughter) {
          // Strict mode: reject the whole resonance if either daughter belongs to the event trigger population.
          if (checkTrackInList(posDauTrack, triggerGIList) || checkTrackInList(negDauTrack, triggerGIList)) {
            continue;
          }
        } else {
          // Standard pair-wise autocorrelation removal: reject only if the CURRENT trigger is a daughter.
          if (trigger.globalIndex() == posDauTrack.globalIndex() || trigger.globalIndex() == negDauTrack.globalIndex()) {
            continue;
          }
        }

        fillResoMassCorrelation<eventType, pairType, roleType>(histReg, trigger, associate, invMass);

        if (massRegion == kMassOutside) {
          continue;
        }

        switch (massRegion) {
          case kMassPeak:
            fillResoRegion<eventType, pairType, roleType, kMassPeak>(histReg, trigger, associate, posDauTrack, negDauTrack, invMass, posDauIdMethod, negDauIdMethod);
            break;
          case kMassLSB:
            fillResoRegion<eventType, pairType, roleType, kMassLSB>(histReg, trigger, associate, posDauTrack, negDauTrack, invMass, posDauIdMethod, negDauIdMethod);
            break;
          case kMassRSB:
            fillResoRegion<eventType, pairType, roleType, kMassRSB>(histReg, trigger, associate, posDauTrack, negDauTrack, invMass, posDauIdMethod, negDauIdMethod);
            break;
          default:
            continue;
        }

        incrementCorrCount<pairType, roleType>(corrCountPerColl);
        incrementResoRegionCount<pairType, roleType>(resoRegionCountsPerColl, massRegion);
        markResoUsed(associate, static_cast<uint32_t>(getResoRegionMaskBit(pairType - kHPhi, PtType, massRegion) >> 18));
      }
    } // associate particle loop is over
  } // function definition is over

  template <typename FilterTrkType, int eventType, int pairType, typename H, typename T, typename A, typename B>
  void executeCorrelation(H& histReg, const T& triggers, const A& associatesLowPt, const B& associatesHighPt, CorrCountArray& corrCountPerColl, ResoRegionCountArray& resoRegionCountsPerColl)
  {
    // Trigger Selection
    // int nTrigger = 0;
    int triggerSelTag = kTrackAccepted;

    // Option 1 : Trigger tracks wont be accepted as daughter track of reconstructed resonance.
    // Option 2 :
    std::vector<int64_t> triggerGIList;
    std::vector<int64_t> posDauGIListPeak;
    std::vector<int64_t> negDauGIListPeak;
    std::vector<int64_t> posDauGIListLSB;
    std::vector<int64_t> negDauGIListLSB;
    std::vector<int64_t> posDauGIListRSSB;
    std::vector<int64_t> negDauGIListRSSB;

    if (c1113rejectResoWithAnyTriggerDaughter) {
      for (const auto& trigger : triggers) {
        if (selectionTrack(trigger) != kTrackAccepted) {
          continue;
        }
        triggerGIList.push_back(trigger.globalIndex());
      }
    }

    for (const auto& trigger : triggers) {
      triggerSelTag = selectionTrack(trigger);
      if (triggerSelTag != kTrackAccepted) {
        continue;
      }
      // fillTrigQA
      fillCorrParticleQA<eventType, pairType, kQA, kTrigger>(histReg, trigger);
      markTrackUsed(trigger, TriggerUsedFlag);
      // ++nTrigger;

      // triggerTrackIndexList.push_back(triggerTrack.globalIndex());
      executeTrigAssocCorrelation<FilterTrkType, eventType, pairType, kAssocLowPt>(histReg, trigger, associatesLowPt, triggerGIList, corrCountPerColl, resoRegionCountsPerColl);
      executeTrigAssocCorrelation<FilterTrkType, eventType, pairType, kAssocHighPt>(histReg, trigger, associatesHighPt, triggerGIList, corrCountPerColl, resoRegionCountsPerColl);
    } // trigger particle loop
  } // execute Correlation Function is over

  // Event Filter
  Filter eventFilter = (!cfgEvent.c1003requireSel8) || (o2::aod::evsel::sel8 == true);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < cfgEvent.c1002cutZvertex);

  // Track Filter
  Filter ptFilter = (o2::aod::track::pt > cfgTrackCuts.c1043ptMin) && (o2::aod::track::pt < cfgTrackCuts.c1044ptMax);
  Filter etaFilter = (nabs(o2::aod::track::eta) < cfgTrackCuts.c1045etaMax);
  Filter dcaFilter = ((!cfgTrackCuts.c1065useFixedDCAxy) || (nabs(o2::aod::track::dcaXY) < cfgTrackCuts.c1066dcaXYMax)) &&
                     ((!cfgTrackCuts.c1067useFixedDCAz) || (nabs(o2::aod::track::dcaZ) < cfgTrackCuts.c1068dcaZMax));
  using MyFilteredCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::CentFV0As, aod::Mults>>; // ,
  // using MyFilteredTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::TOFSignal, aod::pidTOFbeta, aod::pidTOFmass, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullEl, aod::pidTPCFullDe, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullEl, aod::pidTOFFullDe>>;
  using MyFilteredTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::TOFSignal, aod::pidTOFbeta, aod::pidTOFmass, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullDe, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullDe>>;
  // For manual sliceBy
  Preslice<MyFilteredTracks> tracksPerCollisionPreslice = o2::aod::track::collisionId;
  Preslice<aod::ResonanceCndts> resoPerCollisionPreslice = aod::resonancecndt::collisionId; // registers RESONCNDTS/fIndexCollisions in the SliceCache

  // // definition of partitions
  SliceCache cache;
  Partition<MyFilteredTracks> triggerTracks = cfgPartitions.c1077triggerPtLow < aod::track::pt && aod::track::pt < cfgPartitions.c1078triggerPtHigh;
  Partition<MyFilteredTracks> posTracks = aod::track::signed1Pt > 0.0f;
  Partition<MyFilteredTracks> negTracks = aod::track::signed1Pt < 0.0f;

  Partition<MyFilteredTracks> assocTracksLowPt = cfgPartitions.c1079assocPtLowMin < aod::track::pt && aod::track::pt < cfgPartitions.c1080assocPtLowMax;
  Partition<MyFilteredTracks> assocTracksHighPt = cfgPartitions.c1081assocPtHighMin < aod::track::pt && aod::track::pt < cfgPartitions.c1082assocPtHighMax;

  Partition<aod::ResonanceCndts> assocPhiLowPt = (cfgResPartitions.c1085phiPtLowMin < aod::resonancecndt::pt) && (aod::resonancecndt::pt < cfgResPartitions.c1086phiPtLowMax) && (aod::resonancecndt::phi1020Tag != static_cast<uint8_t>(kMassOutside));
  Partition<aod::ResonanceCndts> assocPhiHighPt = (cfgResPartitions.c1087phiPtHighMin < aod::resonancecndt::pt) && (aod::resonancecndt::pt < cfgResPartitions.c1088phiPtHighMax) && (aod::resonancecndt::phi1020Tag != static_cast<uint8_t>(kMassOutside));

  Partition<aod::ResonanceCndts> assocKStarLowPt = (cfgResPartitions.c1089kstarPtLowMin < aod::resonancecndt::pt) && (aod::resonancecndt::pt < cfgResPartitions.c1090kstarPtLowMax) && (aod::resonancecndt::kStar892Tag != static_cast<uint8_t>(kMassOutside));
  Partition<aod::ResonanceCndts> assocKStarHighPt = (cfgResPartitions.c1091kstarPtHighMin < aod::resonancecndt::pt) && (aod::resonancecndt::pt < cfgResPartitions.c1092kstarPtHighMax) && (aod::resonancecndt::kStar892Tag != static_cast<uint8_t>(kMassOutside));

  Partition<aod::ResonanceCndts> assocKStarBarLowPt = (cfgResPartitions.c1093kstarBarPtLowMin < aod::resonancecndt::pt) && (aod::resonancecndt::pt < cfgResPartitions.c1094kstarBarPtLowMax) && (aod::resonancecndt::kStar892BarTag != static_cast<uint8_t>(kMassOutside));
  Partition<aod::ResonanceCndts> assocKStarBarHighPt = (cfgResPartitions.c1095kstarBarPtHighMin < aod::resonancecndt::pt) && (aod::resonancecndt::pt < cfgResPartitions.c1096kstarBarPtHighMax) && (aod::resonancecndt::kStar892BarTag != static_cast<uint8_t>(kMassOutside));

  Partition<aod::ResonanceCndts> assocLambda1520LowPt = (cfgResPartitions.c1097lambda1520PtLowMin < aod::resonancecndt::pt) && (aod::resonancecndt::pt < cfgResPartitions.c1098lambda1520PtLowMax) && (aod::resonancecndt::lambda1520Tag != static_cast<uint8_t>(kMassOutside));
  Partition<aod::ResonanceCndts> assocLambda1520HighPt = (cfgResPartitions.c1099lambda1520PtHighMin < aod::resonancecndt::pt) && (aod::resonancecndt::pt < cfgResPartitions.c1100lambda1520PtHighMax) && (aod::resonancecndt::lambda1520Tag != static_cast<uint8_t>(kMassOutside));

  Partition<aod::ResonanceCndts> assocLambda1520BarLowPt = (cfgResPartitions.c1101lambda1520BarPtLowMin < aod::resonancecndt::pt) && (aod::resonancecndt::pt < cfgResPartitions.c1102lambda1520BarPtLowMax) && (aod::resonancecndt::lambda1520BarTag != static_cast<uint8_t>(kMassOutside));
  Partition<aod::ResonanceCndts> assocLambda1520BarHighPt = (cfgResPartitions.c1103lambda1520BarPtHighMin < aod::resonancecndt::pt) && (aod::resonancecndt::pt < cfgResPartitions.c1104lambda1520BarPtHighMax) && (aod::resonancecndt::lambda1520BarTag != static_cast<uint8_t>(kMassOutside));

  int dfNumber = 0;
  int64_t iCollGlobalCount = -1;

  void processNothing(aod::Origins const& origins)
  {
    if (cfgDebug.c1001printDebugMessages) {
      LOG(info) << "DEBUG :: Process Nothing :: df_" << dfNumber << " :: origins = " << origins.size();
    }
    // Intentionally empty.
    // Keeps the task alive when running purely on derived data.
  }
  PROCESS_SWITCH(HParticleCorrelation2SameEvent, processNothing, "Dummy process for derived-data analysis", true);

  void processSameEvent(MyFilteredCollisions const& collisions, MyFilteredTracks const& fullTracks, aod::ResonanceCndts const& resonanceCndts, o2::aod::Origins const& Origins, aod::BCsWithTimestamps const&)
  {
    dfNumber++;
    mixingQA.fill(HIST("Mixing/DataFramesProcessed"), 0.0);
    std::vector<MixingBinStatus> mixingBinStatusPerDF(nMixBins);
    uint64_t dataframeID = 0;

    if (cfgDerivedData.c1136resetGlobalCountersPerDF) {
      iCollGlobalCount = -1;
    }

    if (Origins.size() == 0) {
      dataframeID = std::numeric_limits<uint64_t>::max() - static_cast<uint64_t>(dfNumber);
      LOG(warn) << "Origins table is empty. Using synthetic dataframeID = " << dataframeID;
    } else {
      dataframeID = Origins.iteratorAt(0).dataframeID();
    }

    if (cfgDebug.c1001printDebugMessages) {
      LOG(info) << "DEBUG :: df_" << dfNumber << " :: Origins = " << Origins.size() << " :: dataframeID = " << dataframeID;
    }

    assocPhiLowPt.bindExternalIndices(&fullTracks);
    assocPhiHighPt.bindExternalIndices(&fullTracks);
    assocKStarLowPt.bindExternalIndices(&fullTracks);
    assocKStarHighPt.bindExternalIndices(&fullTracks);
    assocKStarBarLowPt.bindExternalIndices(&fullTracks);
    assocKStarBarHighPt.bindExternalIndices(&fullTracks);
    assocLambda1520LowPt.bindExternalIndices(&fullTracks);
    assocLambda1520HighPt.bindExternalIndices(&fullTracks);
    assocLambda1520BarLowPt.bindExternalIndices(&fullTracks);
    assocLambda1520BarHighPt.bindExternalIndices(&fullTracks);

    if (collisions.size() == 0) {
      seEventQA.fill(HIST("SE/Events/DataFrameQA"), kDFWithZeroFilteredColls);

      if (cfgDebug.c1001printDebugMessages) {
        LOG(info) << "DEBUG :: df_" << dfNumber << " :: SE :: collisions = 0 :: No filtered collisions found in this dataframe";
      }
    } else if (cfgDebug.c1001printDebugMessages) {
      auto bc = collisions.iteratorAt(0).bc_as<aod::BCsWithTimestamps>();
      int currentRunNumber = bc.runNumber();

      LOG(info) << "DEBUG :: df_" << dfNumber << " :: SE :: currentRunNumber = " << currentRunNumber << " :: " << collisions.iteratorAt(0).bc_as<aod::BCsWithTimestamps>().runNumber();
      LOG(info) << "DEBUG :: df_" << dfNumber << " :: SE :: collisions = " << collisions.size() << " :: fullTracks = " << fullTracks.size() << " :: resonanceCndts = " << resonanceCndts.size();
    }

    int nTrack = 0;
    CollisionRejectionTag collisionTag = kCollAccepted;
    TrackRejectionTag trackTag = kTrackAccepted;

    for (const auto& collision : collisions) { // CollisionLoop-Start
      seEventQA.fill(HIST("SE/Events/EventSelection"), static_cast<float>(kEventAllFiltered));
      fillEventQA<kBeforeSelection>(collision);

      collisionTag = selCollision(collision);
      seEventQA.fill(HIST("SE/Events/CollisionRejection"), static_cast<float>(collisionTag));
      if (collisionTag != kCollAccepted) {
        continue;
      }
      seEventQA.fill(HIST("SE/Events/EventSelection"), static_cast<float>(kEventPassedSelCollision));
      fillEventQA<kAfterSelFunction>(collision);
      const auto tracksPerCollision = fullTracks.sliceBy(tracksPerCollisionPreslice, collision.globalIndex());
      if (tracksPerCollision.size() < cfgEvent.c1005minNFilteredTracks) {
        continue;
      }
      seEventQA.fill(HIST("SE/Events/EventSelection"), static_cast<float>(kEventPassedMinFilteredTracks));

      nTrack = 0;
      for (const auto& track : tracksPerCollision) {
        seTrackQA.fill(HIST("SE/Tracks/TrackSelection"), static_cast<float>(kTrackAllFiltered));
        fillFullTrackQA<kBeforeSelection>(track);

        trackTag = selectionTrack(track);
        seTrackQA.fill(HIST("SE/Tracks/TrackRejection"), static_cast<float>(trackTag));

        if (trackTag != kTrackAccepted) {
          continue;
        }

        seTrackQA.fill(HIST("SE/Tracks/TrackSelection"), static_cast<float>(kTrackPassedSelectionTrack));
        fillFullTrackQA<kAfterSelection>(track);
        nTrack++;
      } // track loop is over.

      if (nTrack < cfgEvent.c1006minNSelectedTracks) {
        continue;
      }

      const auto resonancesPerCollision = resonanceCndts.sliceBy(resoPerCollisionPreslice, collision.globalIndex());
      trackUsedFlags.clear();
      resoUsedFlags.clear();
      if (tracksPerCollision.size() > 0) {
        firstTrackGI = tracksPerCollision.begin().globalIndex();
        trackUsedFlags.assign(tracksPerCollision.iteratorAt(tracksPerCollision.size() - 1).globalIndex() - firstTrackGI + 1, 0);
      }
      if (resonancesPerCollision.size() > 0) {
        firstResoGI = resonancesPerCollision.begin().globalIndex();
        resoUsedFlags.assign(resonancesPerCollision.iteratorAt(resonancesPerCollision.size() - 1).globalIndex() - firstResoGI + 1, 0);
      }

      fillEventQA<kAfterSelection>(collision);
      CorrCountArray corrCountsPerColl{};
      ResoRegionCountArray resoRegionCountsPerColl{};

      seEventQA.fill(HIST("SE/Events/EventSelection"), static_cast<float>(kEventPassedMinSelectedTracks));

      // Get the Partitions
      auto triggerTracksPerColl = triggerTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

      auto assocTrksPerCollLowPt = assocTracksLowPt->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
      auto assocTrksPerCollHighPt = assocTracksHighPt->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

      auto assocPhiPerCollLowPt = assocPhiLowPt->sliceByCached(aod::resonancecndt::collisionId, collision.globalIndex(), cache);
      auto assocPhiPerCollHighPt = assocPhiHighPt->sliceByCached(aod::resonancecndt::collisionId, collision.globalIndex(), cache);

      auto assocKStarPerCollLowPt = assocKStarLowPt->sliceByCached(aod::resonancecndt::collisionId, collision.globalIndex(), cache);
      auto assocKStarPerCollHighPt = assocKStarHighPt->sliceByCached(aod::resonancecndt::collisionId, collision.globalIndex(), cache);

      auto assocKStarBarPerCollLowPt = assocKStarBarLowPt->sliceByCached(aod::resonancecndt::collisionId, collision.globalIndex(), cache);
      auto assocKStarBarPerCollHighPt = assocKStarBarHighPt->sliceByCached(aod::resonancecndt::collisionId, collision.globalIndex(), cache);

      auto assocLambdaPerCollLowPt = assocLambda1520LowPt->sliceByCached(aod::resonancecndt::collisionId, collision.globalIndex(), cache);
      auto assocLambdaPerCollHighPt = assocLambda1520HighPt->sliceByCached(aod::resonancecndt::collisionId, collision.globalIndex(), cache);

      auto assocLambdaBarPerCollLowPt = assocLambda1520BarLowPt->sliceByCached(aod::resonancecndt::collisionId, collision.globalIndex(), cache);
      auto assocLambdaBarPerCollHighPt = assocLambda1520BarHighPt->sliceByCached(aod::resonancecndt::collisionId, collision.globalIndex(), cache);

      currentCentrality = getMixingEstimatorValue(collision);
      executeCorrelation<MyFilteredTracks, kSameEvent, kHH>(hhCorrelation, triggerTracksPerColl, assocTrksPerCollLowPt, assocTrksPerCollHighPt, corrCountsPerColl, resoRegionCountsPerColl);
      executeCorrelation<MyFilteredTracks, kSameEvent, kHId>(hIdCorrelation, triggerTracksPerColl, assocTrksPerCollLowPt, assocTrksPerCollHighPt, corrCountsPerColl, resoRegionCountsPerColl);
      executeCorrelation<MyFilteredTracks, kSameEvent, kHPhi>(hPhiCorrelation, triggerTracksPerColl, assocPhiPerCollLowPt, assocPhiPerCollHighPt, corrCountsPerColl, resoRegionCountsPerColl);
      executeCorrelation<MyFilteredTracks, kSameEvent, kHKStar>(hKStarCorrelation, triggerTracksPerColl, assocKStarPerCollLowPt, assocKStarPerCollHighPt, corrCountsPerColl, resoRegionCountsPerColl);
      executeCorrelation<MyFilteredTracks, kSameEvent, kHKStarBar>(hKStarBarCorrelation, triggerTracksPerColl, assocKStarBarPerCollLowPt, assocKStarBarPerCollHighPt, corrCountsPerColl, resoRegionCountsPerColl);
      executeCorrelation<MyFilteredTracks, kSameEvent, kHLambda>(hLambdaCorrelation, triggerTracksPerColl, assocLambdaPerCollLowPt, assocLambdaPerCollHighPt, corrCountsPerColl, resoRegionCountsPerColl);
      executeCorrelation<MyFilteredTracks, kSameEvent, kHLambdaBar>(hLambdaBarCorrelation, triggerTracksPerColl, assocLambdaBarPerCollLowPt, assocLambdaBarPerCollHighPt, corrCountsPerColl, resoRegionCountsPerColl);

      // Unlike Sign Signal
      // LikeSign Signal
      bool hasSelectedTrigger = false;
      if (c1114requireSelectedTriggerForInvariantMass) {
        for (const auto& trigger : triggerTracksPerColl) {
          if (selectionTrack(trigger) == kTrackAccepted) {
            hasSelectedTrigger = true;
            break;
          }
        }
      }

      if (!c1114requireSelectedTriggerForInvariantMass || hasSelectedTrigger) {
        auto posTracksPerColl = posTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
        auto negTracksPerColl = negTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

        fillUnlikeSignInvariantMass(collision, posTracksPerColl, negTracksPerColl);
        fillLikeSignInvariantMass<true>(collision, posTracksPerColl);
        fillLikeSignInvariantMass<false>(collision, negTracksPerColl);
      }

      //______________________________________________________________________________
      // Event-mixing information for this collision
      uint64_t corrPresenceMask = 0ULL;
      const int mixingBin = getMixingBin(collision);
      const float mixingEstimatorValue = getMixingEstimatorValue(collision);
      const bool validMixingBin = mixingBin >= 0 && mixingBin < nMixBins;

      if (validMixingBin) {
        ++mixingBinStatusPerDF[mixingBin].nCollisions;
        mixingQA.fill(HIST("Mixing/AllCollisionsPerBin"), mixingBin);
        mixingQA.fill(HIST("Mixing/AllCollisionsZCent"), collision.posZ(), mixingEstimatorValue);
      }

      for (int iCorr = 0; iCorr < kNCorrCountTypes; ++iCorr) {
        for (int iPt = 0; iPt < kNCorrCountPt; ++iPt) {
          if (corrCountsPerColl[iCorr][iPt] == 0) {
            continue;
          }

          const int corrChannel = iCorr * kNCorrCountPt + iPt;
          corrPresenceMask |= getCorrMaskBit(iCorr, iPt);

          if (validMixingBin) {
            ++mixingBinStatusPerDF[mixingBin].nCorrCollisions[iCorr][iPt];
            mixingQA.fill(HIST("Mixing/EligibleCollisionsPerBin"), mixingBin, corrChannel);
            mixingQA.fill(HIST("Mixing/EligibleCollisionsZCentChannel"), collision.posZ(), mixingEstimatorValue, corrChannel);
          }
        }
      }

      for (int iReso = 0; iReso < kNResoCorrCountTypes; ++iReso) {
        for (int iPt = 0; iPt < kNCorrCountPt; ++iPt) {
          for (int iMassRegion = kMassPeak; iMassRegion <= kMassRSB; ++iMassRegion) {
            if (resoRegionCountsPerColl[iReso][iPt][iMassRegion] == 0) {
              continue;
            }
            corrPresenceMask |= getResoRegionMaskBit(iReso, iPt, iMassRegion);
          }
        }
      }

      int nTriggersWritten = 0;
      int nTrackAssocsWritten = 0;
      int nResoAssocsWritten = 0;
      //______________________________________________________________________________
      // Derived tables: only the triggers and associates marked as used, and only for the stored pools.
      // Every particle is written once; its UsedMask holds all correlations it was used in.
      const uint64_t storedMask = getStoredCorrMask(corrPresenceMask, mixingBin);
      if (storedMask != 0ULL) {
        derivedCollisions(collision.posZ(), mixingEstimatorValue, storedMask);
        iCollGlobalCount++;

        // Tracks: triggers, and h / pi / K / p associates (track pools = CorrPresenceMask bits 0-7)
        for (const auto& track : tracksPerCollision) {
          const int64_t i = track.globalIndex() - firstTrackGI;
          if (i < 0 || i >= static_cast<int64_t>(trackUsedFlags.size()) || trackUsedFlags[i] == 0) {
            continue;
          }
          const uint16_t flags = trackUsedFlags[i];
          if ((flags & TriggerUsedFlag) != 0u) {
            nTriggersWritten++;
            derivedTriggers(iCollGlobalCount, track.pt(), track.eta(), track.phi());
          }
          const uint64_t usedMask = static_cast<uint64_t>(flags & 0xFFu) & storedMask;
          if (usedMask != 0ULL) {
            nTrackAssocsWritten++;
            derivedAssocs(iCollGlobalCount, track.pt(), track.eta(), track.phi(), static_cast<int8_t>(track.sign() > 0 ? 1 : -1), usedMask);
          }
        }

        // Resonance candidates (resonance pools = CorrPresenceMask bits 18-47; flag bit r = pool bit 18 + r)
        for (const auto& resonance : resonancesPerCollision) {
          const int64_t i = resonance.globalIndex() - firstResoGI;
          if (i < 0 || i >= static_cast<int64_t>(resoUsedFlags.size()) || resoUsedFlags[i] == 0) {
            continue;
          }
          const uint64_t usedMask = (static_cast<uint64_t>(resoUsedFlags[i]) << 18) & storedMask;
          if (usedMask != 0ULL) {
            nResoAssocsWritten++;
            derivedAssocs(iCollGlobalCount, resonance.pt(), resonance.eta(), resonance.phi(), static_cast<int8_t>(0), usedMask);
          }
        }

        mixingQA.fill(HIST("DerivedData/RowsWritten"), 0.0);
        mixingQA.fill(HIST("DerivedData/RowsWritten"), 1.0, nTriggersWritten);
        mixingQA.fill(HIST("DerivedData/RowsWritten"), 2.0, nTrackAssocsWritten);
        mixingQA.fill(HIST("DerivedData/RowsWritten"), 3.0, nResoAssocsWritten);
      } // Derived data creation block is over.
    } // collision loop - Over.

    //______________________________________________________________________________
    // Dataframe-level mixing-pool QA
    std::array<uint64_t, NCorrCountChannels> nReadyMixBins{};

    for (int iMixBin = 0; iMixBin < nMixBins; ++iMixBin) {
      const auto& binStatus = mixingBinStatusPerDF[iMixBin];
      mixingQA.fill(HIST("Mixing/PerDF/AllPoolOccupancy"), iMixBin, binStatus.nCollisions);

      for (int iCorr = 0; iCorr < kNCorrCountTypes; ++iCorr) {
        for (int iPt = 0; iPt < kNCorrCountPt; ++iPt) {
          const int corrChannel = iCorr * kNCorrCountPt + iPt;
          const uint64_t nEligibleCollisions = binStatus.nCorrCollisions[iCorr][iPt];

          mixingQA.fill(HIST("Mixing/PerDF/EligiblePoolOccupancy"), iMixBin, corrChannel, nEligibleCollisions);
          mixingQA.fill(HIST("Mixing/PerDF/PoolOccupancyDistribution"), corrChannel, nEligibleCollisions);

          if (nEligibleCollisions >= static_cast<uint64_t>(cfgMixing.c1137nEvtMixing)) {
            ++nReadyMixBins[corrChannel];
            mixingQA.fill(HIST("Mixing/PerDF/ReadyPoolMap"), iMixBin, corrChannel);
          }
        }
      }
    }

    for (int corrChannel = 0; corrChannel < NCorrCountChannels; ++corrChannel) {
      mixingQA.fill(HIST("Mixing/PerDF/NReadyMixBins"), corrChannel, nReadyMixBins[corrChannel]);
    }
  }
  PROCESS_SWITCH(HParticleCorrelation2SameEvent, processSameEvent, "Process Same event", true);
};

struct HParticleCorrelation3MixedEvent {

  HistogramRegistry hhCorrelation{"hhCorrelation", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry hIdCorrelation{"hIdCorrelation", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry hPhiCorrelation{"hPhiCorrelation", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry hKStarCorrelation{"hKStarCorrelation", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry hKStarBarCorrelation{"hKStarBarCorrelation", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry hLambdaCorrelation{"hLambdaCorrelation", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry hLambdaBarCorrelation{"hLambdaBarCorrelation", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry mixingOperationQA{"mixingOperationQA", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Mixing channels: 8 track channels (h-h, h-pi, h-K, h-p x low/high pT) followed by
  // 30 resonance channels (5 species x low/high pT x Peak/LSB/RSB)
  static constexpr int NMassWindows = 3;
  static constexpr int NTrackMixChannels = (kCountHPr + 1) * static_cast<int>(kNCorrCountPt);
  static constexpr int NMixOperationChannels = NTrackMixChannels + static_cast<int>(kNResoCorrCountTypes) * static_cast<int>(kNCorrCountPt) * NMassWindows;

  static constexpr std::array<std::string_view, NMixOperationChannels> MixOperationChannelName = {
    "hh_LowPt", "hh_HighPt", "hPi_LowPt", "hPi_HighPt", "hKa_LowPt", "hKa_HighPt", "hPr_LowPt", "hPr_HighPt",
    "hPhi_LowPt_Peak", "hPhi_LowPt_LSB", "hPhi_LowPt_RSB", "hPhi_HighPt_Peak", "hPhi_HighPt_LSB", "hPhi_HighPt_RSB",
    "hKStar_LowPt_Peak", "hKStar_LowPt_LSB", "hKStar_LowPt_RSB", "hKStar_HighPt_Peak", "hKStar_HighPt_LSB", "hKStar_HighPt_RSB",
    "hKStarBar_LowPt_Peak", "hKStarBar_LowPt_LSB", "hKStarBar_LowPt_RSB", "hKStarBar_HighPt_Peak", "hKStarBar_HighPt_LSB", "hKStarBar_HighPt_RSB",
    "hLambda_LowPt_Peak", "hLambda_LowPt_LSB", "hLambda_LowPt_RSB", "hLambda_HighPt_Peak", "hLambda_HighPt_LSB", "hLambda_HighPt_RSB",
    "hLambdaBar_LowPt_Peak", "hLambdaBar_LowPt_LSB", "hLambdaBar_LowPt_RSB", "hLambdaBar_HighPt_Peak", "hLambdaBar_HighPt_LSB", "hLambdaBar_HighPt_RSB"};

  struct : ConfigurableGroup {
    Configurable<bool> c2001printDebugMessages{"c2001printDebugMessages", false, "Print debug messages"};
  } cfgDebug;

  struct : ConfigurableGroup {
    ConfigurableAxis c2002axisDeltaPhi{"c2002axisDeltaPhi", {80, -2.0f, 6.0f}, "#Delta#varphi"};
    ConfigurableAxis c2003axisDeltaEta{"c2003axisDeltaEta", {84, -2.1f, 2.1f}, "#Delta#eta"};
    ConfigurableAxis c2004axisCorrSparseTriggerPt{"c2004axisCorrSparseTriggerPt", {8, 4.0f, 8.0f}, "#it{p}_{T}^{trig} (GeV/#it{c})"};
    ConfigurableAxis c2005axisCorrSparseAssocPt{"c2005axisCorrSparseAssocPt", {8, 0.0f, 4.0f}, "#it{p}_{T}^{assoc} (GeV/#it{c})"};
    ConfigurableAxis c2006axisCorrSparseDeltaPhi{"c2006axisCorrSparseDeltaPhi", {32, -2.0f, 6.0f}, "#Delta#varphi"};
    ConfigurableAxis c2007axisCorrSparseDeltaEta{"c2007axisCorrSparseDeltaEta", {24, -2.1f, 2.1f}, "#Delta#eta"};
    ConfigurableAxis c2008axisCorrCentrality{"c2008axisCorrCentrality", {VARIABLE_WIDTH, 0.0, 1.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0}, "Centrality (%) axis of the correlation histograms (estimator chosen by c1138mixingEstimator)"};
  } cfgAxis;

  // Selections (track cuts, PID, pT classes, mass windows, rapidity, estimator) are applied in the same-event
  // task and stored as tags. Only the mixing itself is configured here; keep these axes equal to c1139/c1140.
  struct : ConfigurableGroup {
    Configurable<int> c2009nEvtMixing{"c2009nEvtMixing", 5, "Number of partner events mixed with each trigger event"};
    ConfigurableAxis c2010axisVtxMixing{"c2010axisVtxMixing", {VARIABLE_WIDTH, -10.0, -8.0, -6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0}, "Mixing bins - z vertex"};
    ConfigurableAxis c2011axisCentMixing{"c2011axisCentMixing", {VARIABLE_WIDTH, -1.0, 20.0, 50.0, 80.0, 101.0}, "Mixing bins - centrality (estimator chosen by c1138mixingEstimator in the same-event task)"};
  } cfgMixing;

  // Re-built in init() with the configured axes (as a data member it holds the default axes)
  using BinningTypeMixing = ColumnBinningPolicy<aod::collision::PosZ, aod::hpcorrcollision::Centrality>;
  BinningTypeMixing colBinning{{cfgMixing.c2010axisVtxMixing, cfgMixing.c2011axisCentMixing}, true};
  int nMixBins = 0;

  void init(InitContext const&)
  {
    if (cfgDebug.c2001printDebugMessages) {
      LOGF(info, "Starting init");
    }

    const AxisSpec axisPt = {200, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisEta = {100, -5.0f, 5.0f, "#eta"};
    const AxisSpec axisPhi = {110, -1.0f, 10.0f, "#phi (radians)"};
    const AxisSpec axisSign = {3, -1.5f, 1.5f, "Sign (tracks: charge, resonances: 0)"};

    AxisSpec axisDeltaPhiSpec{cfgAxis.c2002axisDeltaPhi, "#Delta#varphi"};
    AxisSpec axisDeltaEtaSpec{cfgAxis.c2003axisDeltaEta, "#Delta#eta"};
    AxisSpec axisCorrSparseTriggerPt{cfgAxis.c2004axisCorrSparseTriggerPt, "#it{p}_{T}^{trig} (GeV/#it{c})"};
    AxisSpec axisCorrSparseAssocPt{cfgAxis.c2005axisCorrSparseAssocPt, "#it{p}_{T}^{assoc} (GeV/#it{c})"};
    AxisSpec axisCorrSparseDeltaPhi{cfgAxis.c2006axisCorrSparseDeltaPhi, "#Delta#varphi"};
    AxisSpec axisCorrSparseDeltaEta{cfgAxis.c2007axisCorrSparseDeltaEta, "#Delta#eta"};
    AxisSpec axisCorrCentrality{cfgAxis.c2008axisCorrCentrality, "Centrality (%)"};

    auto addAssocQAHistos = [&](auto& histReg, const std::string& basePath) {
      histReg.add((basePath + "Pt").c_str(), "Associate;p_{T} (GeV/c);Counts", HistType::kTH1F, {axisPt});
      histReg.add((basePath + "Eta").c_str(), "Associate;#eta;Counts", HistType::kTH1F, {axisEta});
      histReg.add((basePath + "Phi").c_str(), "Associate;#varphi;Counts", HistType::kTH1F, {axisPhi});
      histReg.add((basePath + "Sign").c_str(), "Associate;Sign;Counts", HistType::kTH1F, {axisSign});
    };

    auto addCorrHistos = [&](auto& histReg, const std::string& basePath) {
      histReg.add((basePath + "DeltaPhi").c_str(), "Trigger-associate;#Delta#varphi;Pairs", HistType::kTH1F, {axisDeltaPhiSpec});
      histReg.add((basePath + "DeltaEta").c_str(), "Trigger-associate;#Delta#eta;Pairs", HistType::kTH1F, {axisDeltaEtaSpec});
      histReg.add((basePath + "DeltaPhiDeltaEta").c_str(), "Trigger-associate;#Delta#varphi;#Delta#eta;Centrality (%)", HistType::kTH3F, {axisDeltaPhiSpec, axisDeltaEtaSpec, axisCorrCentrality});
      histReg.add((basePath + "CorrSparse").c_str(), "Correlation sparse;p_{T}^{trig};p_{T}^{assoc};#Delta#varphi;#Delta#eta;Centrality (%)", HistType::kTHnSparseF, {axisCorrSparseTriggerPt, axisCorrSparseAssocPt, axisCorrSparseDeltaPhi, axisCorrSparseDeltaEta, axisCorrCentrality});
    };

    // Histogram paths mirror the same-event task (SE/... -> ME/...)
    // h-h
    addAssocQAHistos(hhCorrelation, "ME/hh/QA/AssocLowPt/");
    hhCorrelation.addClone("ME/hh/QA/AssocLowPt/", "ME/hh/QA/AssocHighPt/");
    addCorrHistos(hhCorrelation, "ME/hh/Corr/AssocLowPt/");
    hhCorrelation.addClone("ME/hh/Corr/AssocLowPt/", "ME/hh/Corr/AssocHighPt/");

    // h-identified
    addAssocQAHistos(hIdCorrelation, "ME/hIdentified/QA/AssocLowPt/Pi/");
    hIdCorrelation.addClone("ME/hIdentified/QA/AssocLowPt/Pi/", "ME/hIdentified/QA/AssocLowPt/Ka/");
    hIdCorrelation.addClone("ME/hIdentified/QA/AssocLowPt/Pi/", "ME/hIdentified/QA/AssocLowPt/Pr/");
    hIdCorrelation.addClone("ME/hIdentified/QA/AssocLowPt/", "ME/hIdentified/QA/AssocHighPt/");
    addCorrHistos(hIdCorrelation, "ME/hIdentified/Corr/AssocLowPt/Pi/");
    hIdCorrelation.addClone("ME/hIdentified/Corr/AssocLowPt/Pi/", "ME/hIdentified/Corr/AssocLowPt/Ka/");
    hIdCorrelation.addClone("ME/hIdentified/Corr/AssocLowPt/Pi/", "ME/hIdentified/Corr/AssocLowPt/Pr/");
    hIdCorrelation.addClone("ME/hIdentified/Corr/AssocLowPt/", "ME/hIdentified/Corr/AssocHighPt/");

    // h-resonance
    auto registerResonanceHistograms = [&](auto& histReg, const std::string& pairName) {
      const std::string base = "ME/" + pairName + "/";
      addAssocQAHistos(histReg, base + "QA/AssocLowPt/Peak/Reso/");
      histReg.addClone(base + "QA/AssocLowPt/Peak/", base + "QA/AssocLowPt/LSB/");
      histReg.addClone(base + "QA/AssocLowPt/Peak/", base + "QA/AssocLowPt/RSB/");
      histReg.addClone(base + "QA/AssocLowPt/", base + "QA/AssocHighPt/");
      addCorrHistos(histReg, base + "Corr/AssocLowPt/Peak/");
      histReg.addClone(base + "Corr/AssocLowPt/Peak/", base + "Corr/AssocLowPt/LSB/");
      histReg.addClone(base + "Corr/AssocLowPt/Peak/", base + "Corr/AssocLowPt/RSB/");
      histReg.addClone(base + "Corr/AssocLowPt/", base + "Corr/AssocHighPt/");
    };
    registerResonanceHistograms(hPhiCorrelation, "hPhi");
    registerResonanceHistograms(hKStarCorrelation, "hKStar");
    registerResonanceHistograms(hKStarBarCorrelation, "hKStarBar");
    registerResonanceHistograms(hLambdaCorrelation, "hLambda");
    registerResonanceHistograms(hLambdaBarCorrelation, "hLambdaBar");

    //______________________________________________________________________________
    // Mixing binning (re-built here with the configured axes) and mixing QA
    AxisSpec axisVtxMix{cfgMixing.c2010axisVtxMixing, "z_{vtx} (cm)"};
    AxisSpec axisCentMix{cfgMixing.c2011axisCentMixing, "Centrality (%)"};
    nMixBins = axisVtxMix.getNbins() * axisCentMix.getNbins();
    colBinning = BinningTypeMixing{{cfgMixing.c2010axisVtxMixing, cfgMixing.c2011axisCentMixing}, true};

    const AxisSpec axisMixBin = {nMixBins, -0.5, static_cast<double>(nMixBins) - 0.5, "Mixing bin"};
    const AxisSpec axisMixChannel = {NMixOperationChannels, -0.5, static_cast<double>(NMixOperationChannels) - 0.5, "Mixing channel"};
    const AxisSpec axisMixPartners = {static_cast<int>(cfgMixing.c2009nEvtMixing) + 1, -0.5, static_cast<double>(cfgMixing.c2009nEvtMixing) + 0.5, "N mixed partner events"};
    const AxisSpec axisCollPerDF = {200, 0.0, 20000.0, "Derived collisions per dataframe"};
    const AxisSpec axisMixVtxZQA = {100, -10.0f, 10.0f, "z_{vtx} (cm)"};
    const AxisSpec axisMixDeltaVtxZQA = {100, -20.0f, 20.0f, "#Delta z_{vtx} (cm)"};
    const AxisSpec axisMixCentralityQA = {102, -1.0f, 101.0f, "Centrality (%)"};
    const AxisSpec axisMixDeltaCentralityQA = {204, -102.0f, 102.0f, "#Delta Centrality (%)"};

    mixingOperationQA.add("ME/Mixing/DataFramesProcessed", "Processed dataframes;Dataframe;Count", HistType::kTH1D, {{1, -0.5, 0.5}});
    mixingOperationQA.add("ME/Mixing/CollisionsPerDF", "Derived collisions per (merged) dataframe;N collisions;Dataframes", HistType::kTH1D, {axisCollPerDF});
    mixingOperationQA.add("ME/Mixing/CollisionsPerBinChannel", "Collisions in mixing pools;Mixing bin;Mixing channel", HistType::kTH2D, {axisMixBin, axisMixChannel});
    mixingOperationQA.add("ME/Mixing/CollisionPairsPerChannel", "Mixed collision pairs;Mixing channel;Mixed collision pairs", HistType::kTH1D, {axisMixChannel});
    mixingOperationQA.add("ME/Mixing/CollisionPairsPerBinChannel", "Mixed collision pairs;Mixing bin;Mixing channel", HistType::kTH2D, {axisMixBin, axisMixChannel});
    mixingOperationQA.add("ME/Mixing/PartnersPerTriggerEvent", "Mixed partners per trigger event;Mixing channel;N partner events", HistType::kTH2D, {axisMixChannel, axisMixPartners});
    mixingOperationQA.add("ME/Mixing/PartnersPerTriggerEventVsBin", "Mixed partners per trigger event;Mixing bin;Mixing channel;N partner events", HistType::kTH3D, {axisMixBin, axisMixChannel, axisMixPartners});
    mixingOperationQA.add("ME/Mixing/VtxZ1VsVtxZ2", "Mixed-event vertex positions;z_{vtx}^{trig} (cm);z_{vtx}^{assoc} (cm)", HistType::kTH2F, {axisMixVtxZQA, axisMixVtxZQA});
    mixingOperationQA.add("ME/Mixing/DeltaVtxZVsChannel", "Mixed-event #Delta z_{vtx};Mixing channel;z_{vtx}^{trig}-z_{vtx}^{assoc} (cm)", HistType::kTH2F, {axisMixChannel, axisMixDeltaVtxZQA});
    mixingOperationQA.add("ME/Mixing/Centrality1VsCentrality2", "Mixing estimator of mixed events;Estimator^{trig} (%);Estimator^{assoc} (%)", HistType::kTH2F, {axisMixCentralityQA, axisMixCentralityQA});
    mixingOperationQA.add("ME/Mixing/DeltaCentralityVsChannel", "Mixing-estimator difference;Mixing channel;Estimator^{trig}-Estimator^{assoc} (%)", HistType::kTH2F, {axisMixChannel, axisMixDeltaCentralityQA});

    auto setMixOperationChannelLabels = [](TAxis* axis) {
      for (int iChannel = 0; iChannel < NMixOperationChannels; ++iChannel) {
        axis->SetBinLabel(iChannel + 1, MixOperationChannelName[iChannel].data());
      }
    };
    setMixOperationChannelLabels(mixingOperationQA.get<TH2>(HIST("ME/Mixing/CollisionsPerBinChannel"))->GetYaxis());
    setMixOperationChannelLabels(mixingOperationQA.get<TH1>(HIST("ME/Mixing/CollisionPairsPerChannel"))->GetXaxis());
    setMixOperationChannelLabels(mixingOperationQA.get<TH2>(HIST("ME/Mixing/CollisionPairsPerBinChannel"))->GetYaxis());
    setMixOperationChannelLabels(mixingOperationQA.get<TH2>(HIST("ME/Mixing/PartnersPerTriggerEvent"))->GetXaxis());
    setMixOperationChannelLabels(mixingOperationQA.get<TH3>(HIST("ME/Mixing/PartnersPerTriggerEventVsBin"))->GetYaxis());
    setMixOperationChannelLabels(mixingOperationQA.get<TH2>(HIST("ME/Mixing/DeltaVtxZVsChannel"))->GetXaxis());
    setMixOperationChannelLabels(mixingOperationQA.get<TH2>(HIST("ME/Mixing/DeltaCentralityVsChannel"))->GetXaxis());

    if (cfgDebug.c2001printDebugMessages) {
      LOGF(info, "Finishing init");
    }
  }

  SliceCache cache;
  Preslice<aod::HPCorrTriggers> triggersPerCollPreslice = aod::hpcorrtrigger::hpCorrCollId;
  Preslice<aod::HPCorrAssocs> assocsPerCollPreslice = aod::hpcorrassoc::hpCorrCollId; // registers HPCORRASSOC/fIndexHPCorrColls in the SliceCache

  //______________________________________________________________________________
  // Associate partitions: one per mixing channel. Every stored particle carries UsedMask = all correlations
  // h-h, h-pi, h-K, h-p
  Partition<aod::HPCorrAssocs> assocHHLowPt = ncheckbit(aod::hpcorrassoc::usedMask, static_cast<uint64_t>(kMaskHHLowPt));
  Partition<aod::HPCorrAssocs> assocHHHighPt = ncheckbit(aod::hpcorrassoc::usedMask, static_cast<uint64_t>(kMaskHHHighPt));
  Partition<aod::HPCorrAssocs> assocHPiLowPt = ncheckbit(aod::hpcorrassoc::usedMask, static_cast<uint64_t>(kMaskHPiLowPt));
  Partition<aod::HPCorrAssocs> assocHPiHighPt = ncheckbit(aod::hpcorrassoc::usedMask, static_cast<uint64_t>(kMaskHPiHighPt));
  Partition<aod::HPCorrAssocs> assocHKaLowPt = ncheckbit(aod::hpcorrassoc::usedMask, static_cast<uint64_t>(kMaskHKaLowPt));
  Partition<aod::HPCorrAssocs> assocHKaHighPt = ncheckbit(aod::hpcorrassoc::usedMask, static_cast<uint64_t>(kMaskHKaHighPt));
  Partition<aod::HPCorrAssocs> assocHPrLowPt = ncheckbit(aod::hpcorrassoc::usedMask, static_cast<uint64_t>(kMaskHPrLowPt));
  Partition<aod::HPCorrAssocs> assocHPrHighPt = ncheckbit(aod::hpcorrassoc::usedMask, static_cast<uint64_t>(kMaskHPrHighPt));

  // h-Phi
  Partition<aod::HPCorrAssocs> assocHPhiLowPtPeak = ncheckbit(aod::hpcorrassoc::usedMask, static_cast<uint64_t>(kMaskHPhiLowPtPeak));
  Partition<aod::HPCorrAssocs> assocHPhiLowPtLSB = ncheckbit(aod::hpcorrassoc::usedMask, static_cast<uint64_t>(kMaskHPhiLowPtLSB));
  Partition<aod::HPCorrAssocs> assocHPhiLowPtRSB = ncheckbit(aod::hpcorrassoc::usedMask, static_cast<uint64_t>(kMaskHPhiLowPtRSB));
  Partition<aod::HPCorrAssocs> assocHPhiHighPtPeak = ncheckbit(aod::hpcorrassoc::usedMask, static_cast<uint64_t>(kMaskHPhiHighPtPeak));
  Partition<aod::HPCorrAssocs> assocHPhiHighPtLSB = ncheckbit(aod::hpcorrassoc::usedMask, static_cast<uint64_t>(kMaskHPhiHighPtLSB));
  Partition<aod::HPCorrAssocs> assocHPhiHighPtRSB = ncheckbit(aod::hpcorrassoc::usedMask, static_cast<uint64_t>(kMaskHPhiHighPtRSB));

  // h-KStar
  Partition<aod::HPCorrAssocs> assocHKStarLowPtPeak = ncheckbit(aod::hpcorrassoc::usedMask, static_cast<uint64_t>(kMaskHKStarLowPtPeak));
  Partition<aod::HPCorrAssocs> assocHKStarLowPtLSB = ncheckbit(aod::hpcorrassoc::usedMask, static_cast<uint64_t>(kMaskHKStarLowPtLSB));
  Partition<aod::HPCorrAssocs> assocHKStarLowPtRSB = ncheckbit(aod::hpcorrassoc::usedMask, static_cast<uint64_t>(kMaskHKStarLowPtRSB));
  Partition<aod::HPCorrAssocs> assocHKStarHighPtPeak = ncheckbit(aod::hpcorrassoc::usedMask, static_cast<uint64_t>(kMaskHKStarHighPtPeak));
  Partition<aod::HPCorrAssocs> assocHKStarHighPtLSB = ncheckbit(aod::hpcorrassoc::usedMask, static_cast<uint64_t>(kMaskHKStarHighPtLSB));
  Partition<aod::HPCorrAssocs> assocHKStarHighPtRSB = ncheckbit(aod::hpcorrassoc::usedMask, static_cast<uint64_t>(kMaskHKStarHighPtRSB));

  // h-KStarBar
  Partition<aod::HPCorrAssocs> assocHKStarBarLowPtPeak = ncheckbit(aod::hpcorrassoc::usedMask, static_cast<uint64_t>(kMaskHKStarBarLowPtPeak));
  Partition<aod::HPCorrAssocs> assocHKStarBarLowPtLSB = ncheckbit(aod::hpcorrassoc::usedMask, static_cast<uint64_t>(kMaskHKStarBarLowPtLSB));
  Partition<aod::HPCorrAssocs> assocHKStarBarLowPtRSB = ncheckbit(aod::hpcorrassoc::usedMask, static_cast<uint64_t>(kMaskHKStarBarLowPtRSB));
  Partition<aod::HPCorrAssocs> assocHKStarBarHighPtPeak = ncheckbit(aod::hpcorrassoc::usedMask, static_cast<uint64_t>(kMaskHKStarBarHighPtPeak));
  Partition<aod::HPCorrAssocs> assocHKStarBarHighPtLSB = ncheckbit(aod::hpcorrassoc::usedMask, static_cast<uint64_t>(kMaskHKStarBarHighPtLSB));
  Partition<aod::HPCorrAssocs> assocHKStarBarHighPtRSB = ncheckbit(aod::hpcorrassoc::usedMask, static_cast<uint64_t>(kMaskHKStarBarHighPtRSB));

  // h-Lambda
  Partition<aod::HPCorrAssocs> assocHLambdaLowPtPeak = ncheckbit(aod::hpcorrassoc::usedMask, static_cast<uint64_t>(kMaskHLambdaLowPtPeak));
  Partition<aod::HPCorrAssocs> assocHLambdaLowPtLSB = ncheckbit(aod::hpcorrassoc::usedMask, static_cast<uint64_t>(kMaskHLambdaLowPtLSB));
  Partition<aod::HPCorrAssocs> assocHLambdaLowPtRSB = ncheckbit(aod::hpcorrassoc::usedMask, static_cast<uint64_t>(kMaskHLambdaLowPtRSB));
  Partition<aod::HPCorrAssocs> assocHLambdaHighPtPeak = ncheckbit(aod::hpcorrassoc::usedMask, static_cast<uint64_t>(kMaskHLambdaHighPtPeak));
  Partition<aod::HPCorrAssocs> assocHLambdaHighPtLSB = ncheckbit(aod::hpcorrassoc::usedMask, static_cast<uint64_t>(kMaskHLambdaHighPtLSB));
  Partition<aod::HPCorrAssocs> assocHLambdaHighPtRSB = ncheckbit(aod::hpcorrassoc::usedMask, static_cast<uint64_t>(kMaskHLambdaHighPtRSB));

  // h-LambdaBar
  Partition<aod::HPCorrAssocs> assocHLambdaBarLowPtPeak = ncheckbit(aod::hpcorrassoc::usedMask, static_cast<uint64_t>(kMaskHLambdaBarLowPtPeak));
  Partition<aod::HPCorrAssocs> assocHLambdaBarLowPtLSB = ncheckbit(aod::hpcorrassoc::usedMask, static_cast<uint64_t>(kMaskHLambdaBarLowPtLSB));
  Partition<aod::HPCorrAssocs> assocHLambdaBarLowPtRSB = ncheckbit(aod::hpcorrassoc::usedMask, static_cast<uint64_t>(kMaskHLambdaBarLowPtRSB));
  Partition<aod::HPCorrAssocs> assocHLambdaBarHighPtPeak = ncheckbit(aod::hpcorrassoc::usedMask, static_cast<uint64_t>(kMaskHLambdaBarHighPtPeak));
  Partition<aod::HPCorrAssocs> assocHLambdaBarHighPtLSB = ncheckbit(aod::hpcorrassoc::usedMask, static_cast<uint64_t>(kMaskHLambdaBarHighPtLSB));
  Partition<aod::HPCorrAssocs> assocHLambdaBarHighPtRSB = ncheckbit(aod::hpcorrassoc::usedMask, static_cast<uint64_t>(kMaskHLambdaBarHighPtRSB));

  //______________________________________________________________________________
  // Collision pools: one per mixing channel. A collision is in the pool of a channel if its CorrMask has the
  // channel bit, i.e. it had at least one trigger-associate pair of this channel in the same event.
  // h-h, h-pi, h-K, h-p
  Partition<aod::HPCorrColls> collHHLowPt = ncheckbit(aod::hpcorrcollision::corrMask, static_cast<uint64_t>(kMaskHHLowPt));
  Partition<aod::HPCorrColls> collHHHighPt = ncheckbit(aod::hpcorrcollision::corrMask, static_cast<uint64_t>(kMaskHHHighPt));
  Partition<aod::HPCorrColls> collHPiLowPt = ncheckbit(aod::hpcorrcollision::corrMask, static_cast<uint64_t>(kMaskHPiLowPt));
  Partition<aod::HPCorrColls> collHPiHighPt = ncheckbit(aod::hpcorrcollision::corrMask, static_cast<uint64_t>(kMaskHPiHighPt));
  Partition<aod::HPCorrColls> collHKaLowPt = ncheckbit(aod::hpcorrcollision::corrMask, static_cast<uint64_t>(kMaskHKaLowPt));
  Partition<aod::HPCorrColls> collHKaHighPt = ncheckbit(aod::hpcorrcollision::corrMask, static_cast<uint64_t>(kMaskHKaHighPt));
  Partition<aod::HPCorrColls> collHPrLowPt = ncheckbit(aod::hpcorrcollision::corrMask, static_cast<uint64_t>(kMaskHPrLowPt));
  Partition<aod::HPCorrColls> collHPrHighPt = ncheckbit(aod::hpcorrcollision::corrMask, static_cast<uint64_t>(kMaskHPrHighPt));

  // h-Phi
  Partition<aod::HPCorrColls> collHPhiLowPtPeak = ncheckbit(aod::hpcorrcollision::corrMask, static_cast<uint64_t>(kMaskHPhiLowPtPeak));
  Partition<aod::HPCorrColls> collHPhiLowPtLSB = ncheckbit(aod::hpcorrcollision::corrMask, static_cast<uint64_t>(kMaskHPhiLowPtLSB));
  Partition<aod::HPCorrColls> collHPhiLowPtRSB = ncheckbit(aod::hpcorrcollision::corrMask, static_cast<uint64_t>(kMaskHPhiLowPtRSB));
  Partition<aod::HPCorrColls> collHPhiHighPtPeak = ncheckbit(aod::hpcorrcollision::corrMask, static_cast<uint64_t>(kMaskHPhiHighPtPeak));
  Partition<aod::HPCorrColls> collHPhiHighPtLSB = ncheckbit(aod::hpcorrcollision::corrMask, static_cast<uint64_t>(kMaskHPhiHighPtLSB));
  Partition<aod::HPCorrColls> collHPhiHighPtRSB = ncheckbit(aod::hpcorrcollision::corrMask, static_cast<uint64_t>(kMaskHPhiHighPtRSB));

  // h-KStar
  Partition<aod::HPCorrColls> collHKStarLowPtPeak = ncheckbit(aod::hpcorrcollision::corrMask, static_cast<uint64_t>(kMaskHKStarLowPtPeak));
  Partition<aod::HPCorrColls> collHKStarLowPtLSB = ncheckbit(aod::hpcorrcollision::corrMask, static_cast<uint64_t>(kMaskHKStarLowPtLSB));
  Partition<aod::HPCorrColls> collHKStarLowPtRSB = ncheckbit(aod::hpcorrcollision::corrMask, static_cast<uint64_t>(kMaskHKStarLowPtRSB));
  Partition<aod::HPCorrColls> collHKStarHighPtPeak = ncheckbit(aod::hpcorrcollision::corrMask, static_cast<uint64_t>(kMaskHKStarHighPtPeak));
  Partition<aod::HPCorrColls> collHKStarHighPtLSB = ncheckbit(aod::hpcorrcollision::corrMask, static_cast<uint64_t>(kMaskHKStarHighPtLSB));
  Partition<aod::HPCorrColls> collHKStarHighPtRSB = ncheckbit(aod::hpcorrcollision::corrMask, static_cast<uint64_t>(kMaskHKStarHighPtRSB));

  // h-KStarBar
  Partition<aod::HPCorrColls> collHKStarBarLowPtPeak = ncheckbit(aod::hpcorrcollision::corrMask, static_cast<uint64_t>(kMaskHKStarBarLowPtPeak));
  Partition<aod::HPCorrColls> collHKStarBarLowPtLSB = ncheckbit(aod::hpcorrcollision::corrMask, static_cast<uint64_t>(kMaskHKStarBarLowPtLSB));
  Partition<aod::HPCorrColls> collHKStarBarLowPtRSB = ncheckbit(aod::hpcorrcollision::corrMask, static_cast<uint64_t>(kMaskHKStarBarLowPtRSB));
  Partition<aod::HPCorrColls> collHKStarBarHighPtPeak = ncheckbit(aod::hpcorrcollision::corrMask, static_cast<uint64_t>(kMaskHKStarBarHighPtPeak));
  Partition<aod::HPCorrColls> collHKStarBarHighPtLSB = ncheckbit(aod::hpcorrcollision::corrMask, static_cast<uint64_t>(kMaskHKStarBarHighPtLSB));
  Partition<aod::HPCorrColls> collHKStarBarHighPtRSB = ncheckbit(aod::hpcorrcollision::corrMask, static_cast<uint64_t>(kMaskHKStarBarHighPtRSB));

  // h-Lambda
  Partition<aod::HPCorrColls> collHLambdaLowPtPeak = ncheckbit(aod::hpcorrcollision::corrMask, static_cast<uint64_t>(kMaskHLambdaLowPtPeak));
  Partition<aod::HPCorrColls> collHLambdaLowPtLSB = ncheckbit(aod::hpcorrcollision::corrMask, static_cast<uint64_t>(kMaskHLambdaLowPtLSB));
  Partition<aod::HPCorrColls> collHLambdaLowPtRSB = ncheckbit(aod::hpcorrcollision::corrMask, static_cast<uint64_t>(kMaskHLambdaLowPtRSB));
  Partition<aod::HPCorrColls> collHLambdaHighPtPeak = ncheckbit(aod::hpcorrcollision::corrMask, static_cast<uint64_t>(kMaskHLambdaHighPtPeak));
  Partition<aod::HPCorrColls> collHLambdaHighPtLSB = ncheckbit(aod::hpcorrcollision::corrMask, static_cast<uint64_t>(kMaskHLambdaHighPtLSB));
  Partition<aod::HPCorrColls> collHLambdaHighPtRSB = ncheckbit(aod::hpcorrcollision::corrMask, static_cast<uint64_t>(kMaskHLambdaHighPtRSB));

  // h-LambdaBar
  Partition<aod::HPCorrColls> collHLambdaBarLowPtPeak = ncheckbit(aod::hpcorrcollision::corrMask, static_cast<uint64_t>(kMaskHLambdaBarLowPtPeak));
  Partition<aod::HPCorrColls> collHLambdaBarLowPtLSB = ncheckbit(aod::hpcorrcollision::corrMask, static_cast<uint64_t>(kMaskHLambdaBarLowPtLSB));
  Partition<aod::HPCorrColls> collHLambdaBarLowPtRSB = ncheckbit(aod::hpcorrcollision::corrMask, static_cast<uint64_t>(kMaskHLambdaBarLowPtRSB));
  Partition<aod::HPCorrColls> collHLambdaBarHighPtPeak = ncheckbit(aod::hpcorrcollision::corrMask, static_cast<uint64_t>(kMaskHLambdaBarHighPtPeak));
  Partition<aod::HPCorrColls> collHLambdaBarHighPtLSB = ncheckbit(aod::hpcorrcollision::corrMask, static_cast<uint64_t>(kMaskHLambdaBarHighPtLSB));
  Partition<aod::HPCorrColls> collHLambdaBarHighPtRSB = ncheckbit(aod::hpcorrcollision::corrMask, static_cast<uint64_t>(kMaskHLambdaBarHighPtRSB));

  // Number of mixed partners of every collision as trigger event (QA only, re-used by all channels)
  std::vector<int> nPartnersPerColl;

  //______________________________________________________________________________
  // Channel bookkeeping: PairTypeEnum (+ PidEnum for identified hadrons) -> CorrCountTypeEnum,
  // which gives the CorrPresenceMask bit and the mixing-QA channel number
  template <int pairType, int pidType>
  static constexpr int getCorrType()
  {
    if constexpr (pairType == kHH) {
      return kCountHH;
    } else if constexpr (pairType == kHId) {
      return pidType == kPi ? kCountHPi : (pidType == kKa ? kCountHKa : kCountHPr);
    } else if constexpr (pairType == kHPhi) {
      return kCountHPhi;
    } else if constexpr (pairType == kHKStar) {
      return kCountHKStar;
    } else if constexpr (pairType == kHKStarBar) {
      return kCountHKStarBar;
    } else if constexpr (pairType == kHLambda) {
      return kCountHLambda;
    } else {
      return kCountHLambdaBar;
    }
  }

  //______________________________________________________________________________
  // Histogram filling for one mixed trigger-associate pair
  template <int pairType, int roleType, int pidType, int massRegion, typename H, typename T, typename A>
  void fillMixedPair(H& histReg, const T& trigger, const A& associate, const float centrality)
  {
    const float dPhi = computeDeltaPhi(associate.phi(), trigger.phi());
    const float dEta = associate.eta() - trigger.eta();

    if constexpr (pairType == kHH) {
      histReg.fill(HIST(EventTypeDire[kMixedEvent]) + HIST(PairTypeDire[kHH]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST("Pt"), associate.pt());
      histReg.fill(HIST(EventTypeDire[kMixedEvent]) + HIST(PairTypeDire[kHH]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST("Eta"), associate.eta());
      histReg.fill(HIST(EventTypeDire[kMixedEvent]) + HIST(PairTypeDire[kHH]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST("Phi"), associate.phi());
      histReg.fill(HIST(EventTypeDire[kMixedEvent]) + HIST(PairTypeDire[kHH]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST("Sign"), associate.sign());

      histReg.fill(HIST(EventTypeDire[kMixedEvent]) + HIST(PairTypeDire[kHH]) + HIST(AnalysisTypeDire[kCorr]) + HIST(CorrRoleDire[roleType]) + HIST("DeltaPhi"), dPhi);
      histReg.fill(HIST(EventTypeDire[kMixedEvent]) + HIST(PairTypeDire[kHH]) + HIST(AnalysisTypeDire[kCorr]) + HIST(CorrRoleDire[roleType]) + HIST("DeltaEta"), dEta);
      histReg.fill(HIST(EventTypeDire[kMixedEvent]) + HIST(PairTypeDire[kHH]) + HIST(AnalysisTypeDire[kCorr]) + HIST(CorrRoleDire[roleType]) + HIST("DeltaPhiDeltaEta"), dPhi, dEta, centrality);
      histReg.fill(HIST(EventTypeDire[kMixedEvent]) + HIST(PairTypeDire[kHH]) + HIST(AnalysisTypeDire[kCorr]) + HIST(CorrRoleDire[roleType]) + HIST("CorrSparse"), trigger.pt(), associate.pt(), dPhi, dEta, centrality);
    } else if constexpr (pairType == kHId) {
      histReg.fill(HIST(EventTypeDire[kMixedEvent]) + HIST(PairTypeDire[kHId]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(PidTypeDire[pidType]) + HIST("Pt"), associate.pt());
      histReg.fill(HIST(EventTypeDire[kMixedEvent]) + HIST(PairTypeDire[kHId]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(PidTypeDire[pidType]) + HIST("Eta"), associate.eta());
      histReg.fill(HIST(EventTypeDire[kMixedEvent]) + HIST(PairTypeDire[kHId]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(PidTypeDire[pidType]) + HIST("Phi"), associate.phi());
      histReg.fill(HIST(EventTypeDire[kMixedEvent]) + HIST(PairTypeDire[kHId]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(PidTypeDire[pidType]) + HIST("Sign"), associate.sign());

      histReg.fill(HIST(EventTypeDire[kMixedEvent]) + HIST(PairTypeDire[kHId]) + HIST(AnalysisTypeDire[kCorr]) + HIST(CorrRoleDire[roleType]) + HIST(PidTypeDire[pidType]) + HIST("DeltaPhi"), dPhi);
      histReg.fill(HIST(EventTypeDire[kMixedEvent]) + HIST(PairTypeDire[kHId]) + HIST(AnalysisTypeDire[kCorr]) + HIST(CorrRoleDire[roleType]) + HIST(PidTypeDire[pidType]) + HIST("DeltaEta"), dEta);
      histReg.fill(HIST(EventTypeDire[kMixedEvent]) + HIST(PairTypeDire[kHId]) + HIST(AnalysisTypeDire[kCorr]) + HIST(CorrRoleDire[roleType]) + HIST(PidTypeDire[pidType]) + HIST("DeltaPhiDeltaEta"), dPhi, dEta, centrality);
      histReg.fill(HIST(EventTypeDire[kMixedEvent]) + HIST(PairTypeDire[kHId]) + HIST(AnalysisTypeDire[kCorr]) + HIST(CorrRoleDire[roleType]) + HIST(PidTypeDire[pidType]) + HIST("CorrSparse"), trigger.pt(), associate.pt(), dPhi, dEta, centrality);
    } else {
      histReg.fill(HIST(EventTypeDire[kMixedEvent]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(MassRegionDire[massRegion]) + HIST("Reso/Pt"), associate.pt());
      histReg.fill(HIST(EventTypeDire[kMixedEvent]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(MassRegionDire[massRegion]) + HIST("Reso/Eta"), associate.eta());
      histReg.fill(HIST(EventTypeDire[kMixedEvent]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(MassRegionDire[massRegion]) + HIST("Reso/Phi"), associate.phi());
      histReg.fill(HIST(EventTypeDire[kMixedEvent]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[kQA]) + HIST(CorrRoleDire[roleType]) + HIST(MassRegionDire[massRegion]) + HIST("Reso/Sign"), associate.sign());

      histReg.fill(HIST(EventTypeDire[kMixedEvent]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[kCorr]) + HIST(CorrRoleDire[roleType]) + HIST(MassRegionDire[massRegion]) + HIST("DeltaPhi"), dPhi);
      histReg.fill(HIST(EventTypeDire[kMixedEvent]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[kCorr]) + HIST(CorrRoleDire[roleType]) + HIST(MassRegionDire[massRegion]) + HIST("DeltaEta"), dEta);
      histReg.fill(HIST(EventTypeDire[kMixedEvent]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[kCorr]) + HIST(CorrRoleDire[roleType]) + HIST(MassRegionDire[massRegion]) + HIST("DeltaPhiDeltaEta"), dPhi, dEta, centrality);
      histReg.fill(HIST(EventTypeDire[kMixedEvent]) + HIST(PairTypeDire[pairType]) + HIST(AnalysisTypeDire[kCorr]) + HIST(CorrRoleDire[roleType]) + HIST(MassRegionDire[massRegion]) + HIST("CorrSparse"), trigger.pt(), associate.pt(), dPhi, dEta, centrality);
    }
  }

  template <typename C>
  void fillMixingPairQA(const int mixChannel, const int mixBin, const C& collision1, const C& collision2)
  {
    mixingOperationQA.fill(HIST("ME/Mixing/CollisionPairsPerChannel"), mixChannel);
    mixingOperationQA.fill(HIST("ME/Mixing/CollisionPairsPerBinChannel"), mixBin, mixChannel);
    mixingOperationQA.fill(HIST("ME/Mixing/VtxZ1VsVtxZ2"), collision1.posZ(), collision2.posZ());
    mixingOperationQA.fill(HIST("ME/Mixing/DeltaVtxZVsChannel"), mixChannel, collision1.posZ() - collision2.posZ());
    mixingOperationQA.fill(HIST("ME/Mixing/Centrality1VsCentrality2"), collision1.centrality(), collision2.centrality());
    mixingOperationQA.fill(HIST("ME/Mixing/DeltaCentralityVsChannel"), mixChannel, collision1.centrality() - collision2.centrality());
  }

  //______________________________________________________________________________
  // Mix one channel with selfCombinations.
  // collisionPool  : collisions whose CorrMask has the channel bit (partition on HPCorrColls)
  // assocPartition : associates whose UsedMask has the channel bit (partition on HPCorrAssocs)
  // selfCombinations groups the pool by mixing bin (z-vertex x centrality) and pairs every collision with the
  // next c2009nEvtMixing collisions of the same bin. Each pair comes once, collision1 before collision2:
  // the triggers are taken from collision1 and the associates from collision2.
  template <int pairType, int roleType, int pidType, int massRegion, typename H, typename P, typename A>
  void mixChannel(H& histReg, P& collisionPool, A& assocPartition, aod::HPCorrTriggers const& triggers, const int64_t nCollisions)
  {
    constexpr int PtType = roleType == kAssocLowPt ? kCountLowPt : kCountHighPt;
    constexpr int CorrType = getCorrType<pairType, pidType>();
    constexpr bool IsReso = CorrType >= kCountHPhi;
    constexpr int MixChannel = IsReso ? NTrackMixChannels + (CorrType - kCountHPhi) * kNCorrCountPt * NMassWindows + PtType * NMassWindows + (massRegion - kMassPeak)
                                      : CorrType * kNCorrCountPt + PtType;
    static_assert(MixChannel >= 0 && MixChannel < NMixOperationChannels, "Invalid mixing channel");

    nPartnersPerColl.assign(nCollisions, 0);

    for (const auto& [collision1, collision2] : selfCombinations(colBinning, cfgMixing.c2009nEvtMixing, -1, collisionPool, collisionPool)) {
      const int mixBin1 = colBinning.getBin({collision1.posZ(), collision1.centrality()});
      const int mixBin2 = colBinning.getBin({collision2.posZ(), collision2.centrality()});
      if (mixBin1 != mixBin2) {
        LOG(fatal) << "EVENT MIXING ERROR :: collisions from different mixing bins :: " << mixBin1 << " :: " << mixBin2;
      }

      nPartnersPerColl[collision1.globalIndex()]++;
      fillMixingPairQA(MixChannel, mixBin1, collision1, collision2);

      const auto triggers1 = triggers.sliceBy(triggersPerCollPreslice, collision1.globalIndex());
      const auto associates2 = assocPartition->sliceByCached(aod::hpcorrassoc::hpCorrCollId, collision2.globalIndex(), cache);

      for (const auto& trigger : triggers1) {
        for (const auto& associate : associates2) {
          fillMixedPair<pairType, roleType, pidType, massRegion>(histReg, trigger, associate, collision1.centrality());
        }
      }
    }

    // Pool QA: every collision of the pool and the number of partners it was mixed with as trigger event
    for (const auto& collision : collisionPool) {
      const int mixBin = colBinning.getBin({collision.posZ(), collision.centrality()});
      if (mixBin < 0) {
        continue;
      }
      const int nPartners = nPartnersPerColl[collision.globalIndex()];
      mixingOperationQA.fill(HIST("ME/Mixing/CollisionsPerBinChannel"), mixBin, MixChannel);
      mixingOperationQA.fill(HIST("ME/Mixing/PartnersPerTriggerEvent"), MixChannel, nPartners);
      mixingOperationQA.fill(HIST("ME/Mixing/PartnersPerTriggerEventVsBin"), mixBin, MixChannel, nPartners);
    }
  }

  int dfNumber = 0;
  void processNothing(aod::Origins const& origins)
  {
    if (cfgDebug.c2001printDebugMessages) {
      LOG(info) << "DEBUG :: Process Nothing :: df_" << dfNumber << " :: origins = " << origins.size();
    }
    // Intentionally empty. Keeps the task alive when the mixing is not run.
  }
  PROCESS_SWITCH(HParticleCorrelation3MixedEvent, processNothing, "Dummy process for analysis", true);

  void processMixEventInDeriveData(aod::HPCorrColls const& collisions, aod::HPCorrTriggers const& triggers, aod::HPCorrAssocs const& assocs)
  {
    dfNumber++;
    const int64_t nCollisions = collisions.size();
    mixingOperationQA.fill(HIST("ME/Mixing/DataFramesProcessed"), 0.0);
    mixingOperationQA.fill(HIST("ME/Mixing/CollisionsPerDF"), nCollisions);

    if (cfgDebug.c2001printDebugMessages) {
      LOG(info) << "DEBUG :: df_" << dfNumber << " :: collisions = " << nCollisions << " :: triggers = " << triggers.size() << " :: assocs = " << assocs.size();
    }

    // Index sanity check: every trigger / associate must point to a collision of this (merged) dataframe
    int64_t nBadIndex = 0;
    for (const auto& trigger : triggers) {
      if (trigger.hpCorrCollId() < 0 || trigger.hpCorrCollId() >= nCollisions) {
        ++nBadIndex;
      }
    }
    for (const auto& assoc : assocs) {
      if (assoc.hpCorrCollId() < 0 || assoc.hpCorrCollId() >= nCollisions) {
        ++nBadIndex;
      }
    }
    if (nBadIndex > 0) {
      LOG(fatal) << "DERIVED DATA ERROR :: df_" << dfNumber << " :: " << nBadIndex << " trigger/associate rows point outside the collision table";
    }

    // h-h
    mixChannel<kHH, kAssocLowPt, -1, kMassNone>(hhCorrelation, collHHLowPt, assocHHLowPt, triggers, nCollisions);
    mixChannel<kHH, kAssocHighPt, -1, kMassNone>(hhCorrelation, collHHHighPt, assocHHHighPt, triggers, nCollisions);

    // h-pi, h-K, h-p
    mixChannel<kHId, kAssocLowPt, kPi, kMassNone>(hIdCorrelation, collHPiLowPt, assocHPiLowPt, triggers, nCollisions);
    mixChannel<kHId, kAssocHighPt, kPi, kMassNone>(hIdCorrelation, collHPiHighPt, assocHPiHighPt, triggers, nCollisions);
    mixChannel<kHId, kAssocLowPt, kKa, kMassNone>(hIdCorrelation, collHKaLowPt, assocHKaLowPt, triggers, nCollisions);
    mixChannel<kHId, kAssocHighPt, kKa, kMassNone>(hIdCorrelation, collHKaHighPt, assocHKaHighPt, triggers, nCollisions);
    mixChannel<kHId, kAssocLowPt, kPr, kMassNone>(hIdCorrelation, collHPrLowPt, assocHPrLowPt, triggers, nCollisions);
    mixChannel<kHId, kAssocHighPt, kPr, kMassNone>(hIdCorrelation, collHPrHighPt, assocHPrHighPt, triggers, nCollisions);

    // h-Phi
    mixChannel<kHPhi, kAssocLowPt, -1, kMassPeak>(hPhiCorrelation, collHPhiLowPtPeak, assocHPhiLowPtPeak, triggers, nCollisions);
    mixChannel<kHPhi, kAssocLowPt, -1, kMassLSB>(hPhiCorrelation, collHPhiLowPtLSB, assocHPhiLowPtLSB, triggers, nCollisions);
    mixChannel<kHPhi, kAssocLowPt, -1, kMassRSB>(hPhiCorrelation, collHPhiLowPtRSB, assocHPhiLowPtRSB, triggers, nCollisions);
    mixChannel<kHPhi, kAssocHighPt, -1, kMassPeak>(hPhiCorrelation, collHPhiHighPtPeak, assocHPhiHighPtPeak, triggers, nCollisions);
    mixChannel<kHPhi, kAssocHighPt, -1, kMassLSB>(hPhiCorrelation, collHPhiHighPtLSB, assocHPhiHighPtLSB, triggers, nCollisions);
    mixChannel<kHPhi, kAssocHighPt, -1, kMassRSB>(hPhiCorrelation, collHPhiHighPtRSB, assocHPhiHighPtRSB, triggers, nCollisions);

    // h-KStar
    mixChannel<kHKStar, kAssocLowPt, -1, kMassPeak>(hKStarCorrelation, collHKStarLowPtPeak, assocHKStarLowPtPeak, triggers, nCollisions);
    mixChannel<kHKStar, kAssocLowPt, -1, kMassLSB>(hKStarCorrelation, collHKStarLowPtLSB, assocHKStarLowPtLSB, triggers, nCollisions);
    mixChannel<kHKStar, kAssocLowPt, -1, kMassRSB>(hKStarCorrelation, collHKStarLowPtRSB, assocHKStarLowPtRSB, triggers, nCollisions);
    mixChannel<kHKStar, kAssocHighPt, -1, kMassPeak>(hKStarCorrelation, collHKStarHighPtPeak, assocHKStarHighPtPeak, triggers, nCollisions);
    mixChannel<kHKStar, kAssocHighPt, -1, kMassLSB>(hKStarCorrelation, collHKStarHighPtLSB, assocHKStarHighPtLSB, triggers, nCollisions);
    mixChannel<kHKStar, kAssocHighPt, -1, kMassRSB>(hKStarCorrelation, collHKStarHighPtRSB, assocHKStarHighPtRSB, triggers, nCollisions);

    // h-KStarBar
    mixChannel<kHKStarBar, kAssocLowPt, -1, kMassPeak>(hKStarBarCorrelation, collHKStarBarLowPtPeak, assocHKStarBarLowPtPeak, triggers, nCollisions);
    mixChannel<kHKStarBar, kAssocLowPt, -1, kMassLSB>(hKStarBarCorrelation, collHKStarBarLowPtLSB, assocHKStarBarLowPtLSB, triggers, nCollisions);
    mixChannel<kHKStarBar, kAssocLowPt, -1, kMassRSB>(hKStarBarCorrelation, collHKStarBarLowPtRSB, assocHKStarBarLowPtRSB, triggers, nCollisions);
    mixChannel<kHKStarBar, kAssocHighPt, -1, kMassPeak>(hKStarBarCorrelation, collHKStarBarHighPtPeak, assocHKStarBarHighPtPeak, triggers, nCollisions);
    mixChannel<kHKStarBar, kAssocHighPt, -1, kMassLSB>(hKStarBarCorrelation, collHKStarBarHighPtLSB, assocHKStarBarHighPtLSB, triggers, nCollisions);
    mixChannel<kHKStarBar, kAssocHighPt, -1, kMassRSB>(hKStarBarCorrelation, collHKStarBarHighPtRSB, assocHKStarBarHighPtRSB, triggers, nCollisions);

    // h-Lambda
    mixChannel<kHLambda, kAssocLowPt, -1, kMassPeak>(hLambdaCorrelation, collHLambdaLowPtPeak, assocHLambdaLowPtPeak, triggers, nCollisions);
    mixChannel<kHLambda, kAssocLowPt, -1, kMassLSB>(hLambdaCorrelation, collHLambdaLowPtLSB, assocHLambdaLowPtLSB, triggers, nCollisions);
    mixChannel<kHLambda, kAssocLowPt, -1, kMassRSB>(hLambdaCorrelation, collHLambdaLowPtRSB, assocHLambdaLowPtRSB, triggers, nCollisions);
    mixChannel<kHLambda, kAssocHighPt, -1, kMassPeak>(hLambdaCorrelation, collHLambdaHighPtPeak, assocHLambdaHighPtPeak, triggers, nCollisions);
    mixChannel<kHLambda, kAssocHighPt, -1, kMassLSB>(hLambdaCorrelation, collHLambdaHighPtLSB, assocHLambdaHighPtLSB, triggers, nCollisions);
    mixChannel<kHLambda, kAssocHighPt, -1, kMassRSB>(hLambdaCorrelation, collHLambdaHighPtRSB, assocHLambdaHighPtRSB, triggers, nCollisions);

    // h-LambdaBar
    mixChannel<kHLambdaBar, kAssocLowPt, -1, kMassPeak>(hLambdaBarCorrelation, collHLambdaBarLowPtPeak, assocHLambdaBarLowPtPeak, triggers, nCollisions);
    mixChannel<kHLambdaBar, kAssocLowPt, -1, kMassLSB>(hLambdaBarCorrelation, collHLambdaBarLowPtLSB, assocHLambdaBarLowPtLSB, triggers, nCollisions);
    mixChannel<kHLambdaBar, kAssocLowPt, -1, kMassRSB>(hLambdaBarCorrelation, collHLambdaBarLowPtRSB, assocHLambdaBarLowPtRSB, triggers, nCollisions);
    mixChannel<kHLambdaBar, kAssocHighPt, -1, kMassPeak>(hLambdaBarCorrelation, collHLambdaBarHighPtPeak, assocHLambdaBarHighPtPeak, triggers, nCollisions);
    mixChannel<kHLambdaBar, kAssocHighPt, -1, kMassLSB>(hLambdaBarCorrelation, collHLambdaBarHighPtLSB, assocHLambdaBarHighPtLSB, triggers, nCollisions);
    mixChannel<kHLambdaBar, kAssocHighPt, -1, kMassRSB>(hLambdaBarCorrelation, collHLambdaBarHighPtRSB, assocHLambdaBarHighPtRSB, triggers, nCollisions);
  }
  PROCESS_SWITCH(HParticleCorrelation3MixedEvent, processMixEventInDeriveData, "Mixed-event correlations from the minimal derived data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& context)
{
  return WorkflowSpec{
    adaptAnalysisTask<HParticleCorrelation1ResonanceProducer>(context),
    adaptAnalysisTask<HParticleCorrelation2SameEvent>(context),
    adaptAnalysisTask<HParticleCorrelation3MixedEvent>(context)};
}
