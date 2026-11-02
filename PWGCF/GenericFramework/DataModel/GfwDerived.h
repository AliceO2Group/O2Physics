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

/// \file GfwDerived.h
/// \brief Derived collision and track tables for the GFW v02 flow analysis.
///
/// Selection columns are bit masks. A later analysis keeps a row when every bit
/// set in its mask is also set on the row (`ncheckbit` / `(value & mask) == mask`).
/// PID is stored only as bits. The nominal, tight, and loose windows are fixed
/// at production.
///
/// \author Maxim Virta, AU, maxim.virta@cern.ch

#ifndef PWGCF_GENERICFRAMEWORK_DATAMODEL_GFWDERIVED_H_
#define PWGCF_GENERICFRAMEWORK_DATAMODEL_GFWDERIVED_H_

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>

#include <Rtypes.h>

#include <cstdint>

namespace o2::aod
{
namespace gfwcollision
{
DECLARE_SOA_COLUMN(Centrality, centrality, float);     //! Centrality estimator chosen at production
DECLARE_SOA_COLUMN(Multiplicity, multiplicity, float); //! Global-track multiplicity, multNTracksGlobal
DECLARE_SOA_COLUMN(EventSel, eventSel, uint16_t);      //! Event-selection bitmask, see EventSelection

/// One bit per event-level decision. Bits record that the decision passed.
/// sel8 is required by the producer and is not stored here.
enum EventSelection : uint16_t {
  kNoSameBunchPileup = BIT(0),
  kIsGoodZvtxFT0vsPV = BIT(1),
  kNoCollInTimeRangeStandard = BIT(2),
  kNoCollInRofStandard = BIT(3),
  kNoHighMultCollInPrevRof = BIT(4),
  kIsVertexITSTPC = BIT(5),
  kIsGoodITSLayersAll = BIT(6),
  kNoTimeFrameBorder = BIT(7),
  kNoITSROFrameBorder = BIT(8),
  kOccupancy = BIT(9),    //! track occupancy is within the production threshold
  kNotTVXinTRD = BIT(10), //! collision is not a TRD-triggered TVX event
  kGoodVertex = BIT(11),  //! PV contributor and z-resolution quality
  kMultCorr = BIT(12)     //! multiplicity-correlation cut
};

inline constexpr uint16_t DefaultEventMask = 0;
} // namespace gfwcollision

DECLARE_SOA_TABLE(GfwCollisions, "AOD", "GFWCOLLISION", //!
                  o2::soa::Index<>,
                  bc::RunNumber,
                  collision::PosZ,
                  gfwcollision::Centrality,
                  gfwcollision::Multiplicity,
                  timestamp::Timestamp,
                  gfwcollision::EventSel);
using GfwCollision = GfwCollisions::iterator;

namespace gfwtrack
{
DECLARE_SOA_INDEX_COLUMN(GfwCollision, gfwCollision); //! Parent derived collision
DECLARE_SOA_COLUMN(Pt, pt, float);                    //! pT (GeV/c)
DECLARE_SOA_COLUMN(Eta, eta, float);                  //! Pseudorapidity
DECLARE_SOA_COLUMN(Phi, phi, float);                  //! Azimuthal angle
DECLARE_SOA_COLUMN(HasTOF, hasTOF, bool);             //! Track has a TOF match
DECLARE_SOA_COLUMN(TrackSel, trackSel, uint32_t);     //! Track and PID bitmask, see TrackSelection

/// Track-quality bits occupy 0-13. TPC PID bits occupy 14-22 and TOF PID bits occupy 23-31.
/// Global-track, ITS-cluster, and ITS-chi2 requirements are applied by the producer and are not stored.
/// DCA bits are |DCA| <= nSigma times a pT-dependent shape.
/// A PID bit is set when low * scale < nSigma < up * scale.
/// Above the TOF pT cut, every PID bit stays unset when the track has no TOF match.
/// Within each detector, PidWindowBits covers π/K/p, followed by the tight window and then the loose window.
enum TrackSelection : uint32_t {
  kDCAxy = BIT(0), //! |DCAxy| within the nominal nSigma times f_xy(pT)
  kDCAxyTight = BIT(1),
  kDCAxyLoose = BIT(2),
  kDCAz = BIT(3), //! |DCAz| within the nominal nSigma times f_z(pT)
  kDCAzTight = BIT(4),
  kDCAzLoose = BIT(5), // DCA bit selections

  kTPCXrows = BIT(6),
  kTPCXrowsTight = BIT(7),
  kTPCXrowsLoose = BIT(8),
  kTPCCls = BIT(9),
  kTPCClsTight = BIT(10),
  kTPCClsLoose = BIT(11),
  kChi2TPC = BIT(12),
  kChi2TPCLoose = BIT(13), // TPC bit selections

  // TPC PID selections: nominal, tight, then loose
  kTpcPi = BIT(14),
  kTpcKa = BIT(15),
  kTpcPr = BIT(16),
  kTpcPiTight = BIT(17),
  kTpcKaTight = BIT(18),
  kTpcPrTight = BIT(19),
  kTpcPiLoose = BIT(20),
  kTpcKaLoose = BIT(21),
  kTpcPrLoose = BIT(22),
  // TOF PID selections: nominal, tight, then loose
  kTofPi = BIT(23),
  kTofKa = BIT(24),
  kTofPr = BIT(25),
  kTofPiTight = BIT(26),
  kTofKaTight = BIT(27),
  kTofPrTight = BIT(28),
  kTofPiLoose = BIT(29),
  kTofKaLoose = BIT(30),
  kTofPrLoose = BIT(31)
};

/// Quality bits corresponding to the default flowGfwV02 track filter plus trackSelected().
inline constexpr uint32_t NominalMask = static_cast<uint32_t>(kDCAxy | kDCAz | kTPCXrows | kTPCCls | kChi2TPC);
/// Loosest quality bits. Production can optionally require these while still leaving the nominal cuts free.
inline constexpr uint32_t LooseSkimMask = static_cast<uint32_t>(kDCAxyLoose | kDCAzLoose | kTPCXrowsLoose | kTPCClsLoose | kChi2TPCLoose);
/// Tracks that enter the multiplicity-correlation count. Global track and ITS chi2 are already required by the producer.
inline constexpr uint32_t MultCorrTrackMask = static_cast<uint32_t>(kDCAz | kChi2TPC);

inline constexpr int PidWindowBits = 3;
inline constexpr int PidWindowNominal = 0;
inline constexpr int PidWindowTight = 1;
inline constexpr int PidWindowLoose = 2;
} // namespace gfwtrack

DECLARE_SOA_TABLE(GfwTracks, "AOD", "GFWTRACK", //!
                  o2::soa::Index<>,
                  gfwtrack::GfwCollisionId,
                  gfwtrack::Pt, gfwtrack::Eta, gfwtrack::Phi,
                  gfwtrack::HasTOF, gfwtrack::TrackSel);
using GfwTrack = GfwTracks::iterator;
} // namespace o2::aod

#endif // PWGCF_GENERICFRAMEWORK_DATAMODEL_GFWDERIVED_H_
