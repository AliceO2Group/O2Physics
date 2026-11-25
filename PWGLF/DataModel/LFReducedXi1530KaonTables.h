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

/// \file LFReducedXi1530KaonTables.h
/// \brief Reduced data model for Xi(1530)^0--K femtoscopy.
///
/// The Xi(1530) candidate stores two independent uint64_t selection words:
///   * TrackSelectionBits: daughter PID / track-quality / prompt-pion DCAz;
///   * TopologySelectionBits: Lambda/V0 and Xi/cascade topology, mass and
///     lifetime working points.
///
/// This split keeps the table compact while retaining both:
///   * the Xi(1530) analysis-note working points/systematic variations; and
///   * the pp Xi-candidate selections used in Run 3 strangeness production.
///
///
/// External-kaon DCA choices remain in a compact uint8_t bitmap; raw TPC/TOF
/// n-sigma values are retained for downstream PID studies.

#ifndef PWGLF_DATAMODEL_LFREDUCEDXI1530KAONTABLES_H_
#define PWGLF_DATAMODEL_LFREDUCEDXI1530KAONTABLES_H_

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>

#include <cmath>
#include <cstdint>

namespace o2::aod
{

namespace redxistark
{
// +1: Xi- pi+ US; -1: anti-Xi+ pi- US; +2/-2: corresponding LS controls.
enum XiStarChannel : int8_t {
  kXiStarUS = +1,
  kAntiXiStarUS = -1,
  kXiStarLS = +2,
  kAntiXiStarLS = -2
};

// -----------------------------------------------------------------------------
// TRACK/PID SELECTION WORD
// -----------------------------------------------------------------------------
// Bit positions only; numerical thresholds live in producer Configurables.
enum XiStarTrackSelBit : uint8_t {
  // Prompt Xi* pion, recipe A:
  // no TOF -> TPC only; TOF -> circular sqrt(TPC^2+TOF^2).
  kPiFirstPID3Sigma = 0,
  kPiFirstPID4Sigma,
  kPiFirstPID5Sigma,
  kPiFirstPID6Sigma,

  // Prompt Xi* pion, recipe B:
  // TPC-only below configurable pT threshold, circular above it.
  kPiFirstPtThresholdPID3Sigma,
  kPiFirstPtThresholdPID4Sigma,
  kPiFirstPtThresholdPID5Sigma,
  kPiFirstPtThresholdPID6Sigma,

  // Prompt Xi* pion, pure TPC-only PID.  These bits allow the analysis
  // to reproduce a TPC-only "all pion tracks" PID strategy downstream,
  // independent of whether the prompt pion also has TOF information.
  kPiFirstTPC3Sigma,
  kPiFirstTPC4Sigma,
  kPiFirstTPC5Sigma,
  kPiFirstTPC6Sigma,

  // Xi bachelor pion, availability-based hybrid PID.
  kXiBachelorPiHybridPID3Sigma,
  kXiBachelorPiHybridPID4Sigma,
  kXiBachelorPiHybridPID5Sigma,
  kXiBachelorPiHybridPID6Sigma,

  // Xi bachelor pion, TPC-only PID.  The 4.8-sigma point reproduces the
  // dedicated pp-Xi table.
  kXiBachelorPiTPC3Sigma,
  kXiBachelorPiTPC4Sigma,
  kXiBachelorPiTPC4p8Sigma,
  kXiBachelorPiTPC5Sigma,
  kXiBachelorPiTPC6Sigma,

  // Lambda daughter pion: TPC only.
  kLambdaPiTPC3Sigma,
  kLambdaPiTPC4Sigma,
  kLambdaPiTPC4p8Sigma,
  kLambdaPiTPC5Sigma,
  kLambdaPiTPC6Sigma,

  // Lambda daughter proton: TPC only.
  kLambdaPrTPC3Sigma,
  kLambdaPrTPC4Sigma,
  kLambdaPrTPC5Sigma,
  kLambdaPrTPC6Sigma,

  // Prompt-pion quality.
  kPiFirstTPCRows70,
  kPiFirstTPCRows80,
  kPiFirstTPCRows90,
  kPiFirstDCAz1cm,
  kPiFirstDCAzDefault,
  kPiFirstDCAz0p1cm,

  // Cascade-daughter quality.
  kAllCascDaughtersTPCRows50,
  kV0DaughtersTPCRows70,
  kV0DaughtersTPCRowsVar1,
  kV0DaughtersTPCRowsVar2,

  kNXiStarTrackSelBits
};
static_assert(kNXiStarTrackSelBits <= 64, "Xi(1530) track-selection bitmap exceeds uint64_t capacity");

constexpr uint64_t xiStarTrackSelMask(XiStarTrackSelBit bit)
{
  return uint64_t{1} << static_cast<uint8_t>(bit);
}

// -----------------------------------------------------------------------------
// TOPOLOGY SELECTION WORD
// -----------------------------------------------------------------------------
enum XiStarTopoSelBit : uint8_t {
  // V0 daughter DCA to PV.
  kV0PionDcaPV005 = 0, // > 0.05 cm
  kV0PionDcaPV006,     // > 0.06 cm
  kV0PionDcaPV010,     // > 0.10 cm
  kV0PionDcaPV020,     // > 0.20 cm (dedicated Lambda)

  kV0ProtonDcaPV005, // > 0.05 cm
  kV0ProtonDcaPV006, // > 0.06 cm
  kV0ProtonDcaPV007, // > 0.07 cm (dedicated Lambda)
  kV0ProtonDcaPV010, // > 0.10 cm

  // Keep the Xi(1530)-analysis convention requested by the analysis:
  // MINIMUM Lambda/V0 DCA to PV.
  kV0DcaPV000, // > 0.00 cm
  kV0DcaPV003, // > 0.03 cm
  kV0DcaPV010, // > 0.10 cm

  // DCA between V0 daughters / native fitter metric.
  kV0DcaDaughters1p0, // < 1.0
  kV0DcaDaughters0p5, // < 0.5
  kV0DcaDaughters0p1, // < 0.1

  // V0 pointing angle.
  kV0CosPA097,
  kV0CosPA098,
  kV0CosPA09876, // pp-Xi selection
  kV0CosPA099,
  kV0CosPA0995, // dedicated Lambda selection

  // Minimum V0 radius.
  kV0Radius0p9,
  kV0Radius1p01,
  kV0Radius1p2, // pp-Xi selection
  kV0Radius2p5,
  kV0Radius3p0, // dedicated Lambda selection

  // Lambda mass window.
  kLambdaMass10MeV,
  kLambdaMass8MeV,
  kLambdaMass6MeV,
  kLambdaMass11p6MeV, // pp-Xi selection

  // Lambda proper lifetime. Var1/Var2 are reserved for future systematic
  // values; their defaults can equal the nominal until values are supplied.
  kV0LifetimeDefault,
  kV0LifetimeVar1,
  kV0LifetimeVar2,

  // Xi bachelor DCA to PV.
  kCascBachelorDcaPV005,
  kCascBachelorDcaPV006,
  kCascBachelorDcaPV010,

  // DCA between bachelor and V0 / native cascade fitter metric.
  kCascDcaDaughters1p0,
  kCascDcaDaughters0p25,
  kCascDcaDaughters0p20,

  // Cascade pointing angle.
  kCascCosPA097,
  kCascCosPA098,
  kCascCosPA09947, // pp-Xi selection
  kCascCosPA0995,

  // Minimum cascade radius.
  kCascRadius0p9,
  kCascRadius1p0, // pp-Xi selection
  kCascRadius1p01,
  kCascRadius1p3,

  // Xi mass window used by the Xi(1530) analysis.
  kXiMass10MeV,
  kXiMass8MeV,
  kXiMass6MeV,

  // Xi proper-lifetime working points (mL/p in cm).
  // Numerical thresholds are Configurable in the producer.
  kXiLifetimeDefault, // default: 22.6 cm (~4.6*c*tau_Xi)
  kXiLifetimeVar1,    // default variation: 15 cm
  kXiLifetimeVar2,    // default variation: 12 cm

  // Additional pp-Xi cuts.
  kBachBaryonDCAxy0020, // > 0.020 cm
  kXiRapidity05,        // |y_Xi| < 0.5

  kNXiStarTopoSelBits
};
static_assert(kNXiStarTopoSelBits <= 64, "Xi(1530) topology-selection bitmap exceeds uint64_t capacity");

constexpr uint64_t xiStarTopoSelMask(XiStarTopoSelBit bit)
{
  return uint64_t{1} << static_cast<uint8_t>(bit);
}

} // namespace redxistark

// -----------------------------------------------------------------------------
// Event table
// -----------------------------------------------------------------------------
namespace redxistarkevent
{
DECLARE_SOA_COLUMN(FT0MPercentile, fT0MPercentile, float); //! FT0M multiplicity percentile (%)
DECLARE_SOA_COLUMN(Bz, bz, float);                         //! nominal L3 magnetic field along z (T)
} // namespace redxistarkevent

DECLARE_SOA_TABLE(RedXiStarKEvents, "AOD", "REDXISTARKEVENT",
                  o2::soa::Index<>,
                  collision::PosZ,
                  collision::NumContrib,
                  redxistarkevent::FT0MPercentile,
                  redxistarkevent::Bz);
using RedXiStarKEvent = RedXiStarKEvents::iterator;

// -----------------------------------------------------------------------------
// Xi(1530) candidate table
// -----------------------------------------------------------------------------
namespace redxistarcandidate
{
DECLARE_SOA_INDEX_COLUMN(RedXiStarKEvent, redXiStarKEvent);

DECLARE_SOA_COLUMN(Channel, channel, int8_t);
DECLARE_SOA_COLUMN(Px, px, float);
DECLARE_SOA_COLUMN(Py, py, float);
DECLARE_SOA_COLUMN(Pz, pz, float);
DECLARE_SOA_COLUMN(Mass, mass, float); //! reconstructed M(Xi pi)

DECLARE_SOA_COLUMN(TrackSelectionBits, trackSelectionBits, uint64_t);
DECLARE_SOA_COLUMN(TopologySelectionBits, topologySelectionBits, uint64_t);

// Prompt-pion kinematics are stored for downstream CPR QA.  The three
// displaced Xi daughters remain ID-only.
DECLARE_SOA_COLUMN(XiStarPionPt, xiStarPionPt, float);
DECLARE_SOA_COLUMN(XiStarPionEta, xiStarPionEta, float);
DECLARE_SOA_COLUMN(XiStarPionPhi, xiStarPionPhi, float);
DECLARE_SOA_COLUMN(XiStarPionIndex, xiStarPionIndex, int64_t);
DECLARE_SOA_COLUMN(XiBachelorIndex, xiBachelorIndex, int64_t);
DECLARE_SOA_COLUMN(V0PositiveIndex, v0PositiveIndex, int64_t);
DECLARE_SOA_COLUMN(V0NegativeIndex, v0NegativeIndex, int64_t);

DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, [](float px, float py) -> float {
  return std::sqrt(px * px + py * py);
});
DECLARE_SOA_DYNAMIC_COLUMN(P, p, [](float px, float py, float pz) -> float {
  return std::sqrt(px * px + py * py + pz * pz);
});
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, [](float px, float py, float pz) -> float {
  const float pt = std::sqrt(px * px + py * py);
  return pt > 0.f ? std::asinh(pz / pt) : 0.f;
});
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, [](float px, float py) -> float {
  return std::atan2(py, px);
});
} // namespace redxistarcandidate

DECLARE_SOA_TABLE(XiStarCandidates, "AOD", "REDXSTARCAND",
                  o2::soa::Index<>,
                  redxistarcandidate::RedXiStarKEventId,
                  redxistarcandidate::Channel,
                  redxistarcandidate::Px,
                  redxistarcandidate::Py,
                  redxistarcandidate::Pz,
                  redxistarcandidate::Mass,
                  redxistarcandidate::TrackSelectionBits,
                  redxistarcandidate::TopologySelectionBits,
                  redxistarcandidate::XiStarPionPt,
                  redxistarcandidate::XiStarPionEta,
                  redxistarcandidate::XiStarPionPhi,
                  redxistarcandidate::XiStarPionIndex,
                  redxistarcandidate::XiBachelorIndex,
                  redxistarcandidate::V0PositiveIndex,
                  redxistarcandidate::V0NegativeIndex,
                  redxistarcandidate::Pt<redxistarcandidate::Px, redxistarcandidate::Py>,
                  redxistarcandidate::P<redxistarcandidate::Px, redxistarcandidate::Py, redxistarcandidate::Pz>,
                  redxistarcandidate::Eta<redxistarcandidate::Px, redxistarcandidate::Py, redxistarcandidate::Pz>,
                  redxistarcandidate::Phi<redxistarcandidate::Px, redxistarcandidate::Py>);
using XiStarCandidate = XiStarCandidates::iterator;

// -----------------------------------------------------------------------------
// External primary-kaon table
// -----------------------------------------------------------------------------
namespace redxistarkaon
{
DECLARE_SOA_INDEX_COLUMN(RedXiStarKEvent, redXiStarKEvent);

enum KaonSelBit : uint8_t {
  kKaonDCAxyDefault = 0,
  kKaonDCAxyVar1,
  kKaonDCAxyVar2,
  kKaonDCAzDefault,
  kKaonDCAzVar1,
  kKaonDCAzVar2,
  kNKaonSelBits
};
static_assert(kNKaonSelBits <= 8, "Kaon selection bitmap exceeds uint8_t capacity");

constexpr uint8_t kaonSelMask(KaonSelBit bit)
{
  return uint8_t{1} << static_cast<uint8_t>(bit);
}

DECLARE_SOA_COLUMN(Charge, charge, int8_t);
DECLARE_SOA_COLUMN(Px, px, float);
DECLARE_SOA_COLUMN(Py, py, float);
DECLARE_SOA_COLUMN(Pz, pz, float);
DECLARE_SOA_COLUMN(TrackIndex, trackIndex, int64_t);
DECLARE_SOA_COLUMN(SelectionBits, selectionBits, uint8_t);

DECLARE_SOA_COLUMN(TPCNClsCrossedRows, tPCNClsCrossedRows, uint8_t);
DECLARE_SOA_COLUMN(TPCNSigmaKa, tPCNSigmaKa, float);
DECLARE_SOA_COLUMN(TOFNSigmaKa, tOFNSigmaKa, float);
DECLARE_SOA_COLUMN(HasTOF, hasTOF, bool);

DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, [](float px, float py) -> float {
  return std::sqrt(px * px + py * py);
});
DECLARE_SOA_DYNAMIC_COLUMN(P, p, [](float px, float py, float pz) -> float {
  return std::sqrt(px * px + py * py + pz * pz);
});
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, [](float px, float py, float pz) -> float {
  const float pt = std::sqrt(px * px + py * py);
  return pt > 0.f ? std::asinh(pz / pt) : 0.f;
});
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, [](float px, float py) -> float {
  return std::atan2(py, px);
});
} // namespace redxistarkaon

DECLARE_SOA_TABLE(KaonCandidates, "AOD", "REDXSTARKAON",
                  o2::soa::Index<>,
                  redxistarkaon::RedXiStarKEventId,
                  redxistarkaon::Charge,
                  redxistarkaon::Px,
                  redxistarkaon::Py,
                  redxistarkaon::Pz,
                  redxistarkaon::TrackIndex,
                  redxistarkaon::SelectionBits,
                  redxistarkaon::TPCNClsCrossedRows,
                  redxistarkaon::TPCNSigmaKa,
                  redxistarkaon::TOFNSigmaKa,
                  redxistarkaon::HasTOF,
                  redxistarkaon::Pt<redxistarkaon::Px, redxistarkaon::Py>,
                  redxistarkaon::P<redxistarkaon::Px, redxistarkaon::Py, redxistarkaon::Pz>,
                  redxistarkaon::Eta<redxistarkaon::Px, redxistarkaon::Py, redxistarkaon::Pz>,
                  redxistarkaon::Phi<redxistarkaon::Px, redxistarkaon::Py>);
using KaonCandidate = KaonCandidates::iterator;

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFREDUCEDXI1530KAONTABLES_H_
