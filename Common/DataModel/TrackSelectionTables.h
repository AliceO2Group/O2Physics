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

#ifndef COMMON_DATAMODEL_TRACKSELECTIONTABLES_H_
#define COMMON_DATAMODEL_TRACKSELECTIONTABLES_H_

#include <Framework/AnalysisDataModel.h>

#include <cstdint>

namespace o2::aod
{
namespace fwdtrack
{
DECLARE_SOA_COLUMN(FwdDcaX, fwdDcaX, float); //! Impact parameter in X of forward track to the primary vertex
DECLARE_SOA_COLUMN(FwdDcaY, fwdDcaY, float); //! Impact parameter in Y of forward track to the primary vertex
} // namespace fwdtrack

namespace track
{
// Columns to store the DCA to the primary vertex
DECLARE_SOA_COLUMN(DcaXY, dcaXY, float);             //! Impact parameter in XY of the track to the primary vertex
DECLARE_SOA_COLUMN(DcaZ, dcaZ, float);               //! Impact parameter in Z of the track to the primary vertex
DECLARE_SOA_COLUMN(SigmaDcaXY2, sigmaDcaXY2, float); //! Impact parameter sigma^2 in XY of the track to the primary vertex
DECLARE_SOA_COLUMN(SigmaDcaZ2, sigmaDcaZ2, float);   //! Impact parameter sigma^2 in Z of the track to the primary vertex

struct TrackSelectionFlags {
 public:
  typedef uint16_t flagtype;
  // Single cut masks
  static constexpr flagtype kTrackType = 1 << 0;
  static constexpr flagtype kPtRange = 1 << 1;
  static constexpr flagtype kEtaRange = 1 << 2;
  static constexpr flagtype kTPCNCls = 1 << 3;
  static constexpr flagtype kTPCCrossedRows = 1 << 4;
  static constexpr flagtype kTPCCrossedRowsOverNCls = 1 << 5;
  static constexpr flagtype kTPCChi2NDF = 1 << 6;
  static constexpr flagtype kTPCRefit = 1 << 7;
  static constexpr flagtype kITSNCls = 1 << 8;
  static constexpr flagtype kITSChi2NDF = 1 << 9;
  static constexpr flagtype kITSRefit = 1 << 10;
  static constexpr flagtype kITSHits = 1 << 11;
  static constexpr flagtype kGoldenChi2 = 1 << 12;
  static constexpr flagtype kDCAxy = 1 << 13;
  static constexpr flagtype kDCAz = 1 << 14;
  // Combo masks
  static constexpr flagtype kQualityTracksITS = kTrackType | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits;
  static constexpr flagtype kQualityTracksTPC = kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit;
  static constexpr flagtype kQualityTracks = kTrackType | kQualityTracksITS | kQualityTracksTPC;
  static constexpr flagtype kQualityTracksWoTPCCluster = kQualityTracksITS | kTPCChi2NDF | kTPCRefit;
  static constexpr flagtype kPrimaryTracks = kGoldenChi2 | kDCAxy | kDCAz;
  static constexpr flagtype kInAcceptanceTracks = kPtRange | kEtaRange;
  static constexpr flagtype kGlobalTrack = kQualityTracks | kPrimaryTracks | kInAcceptanceTracks;
  static constexpr flagtype kGlobalTrackWoTPCCluster = kQualityTracksWoTPCCluster | kPrimaryTracks | kInAcceptanceTracks;
  static constexpr flagtype kGlobalTrackWoPtEta = kQualityTracks | kPrimaryTracks;
  static constexpr flagtype kGlobalTrackWoDCA = kQualityTracks | kInAcceptanceTracks;
  static constexpr flagtype kGlobalTrackWoDCAxy = kQualityTracks | kInAcceptanceTracks | kDCAz;
  static constexpr flagtype kGlobalTrackWoDCATPCCluster = kQualityTracksWoTPCCluster | kInAcceptanceTracks;

  /// @brief Function to check flag content
  /// @param flags bitmask contained in the track
  /// @param mask bitmask to check against the one contained in the track
  /// @return true if the check is successful
  static bool checkFlag(const TrackSelectionFlags::flagtype flags,
                        const TrackSelectionFlags::flagtype mask)
  {
    return (flags & mask) == mask;
  }
};

#define requireTrackCutInFilter(mask) ((o2::aod::track::trackCutFlag & o2::aod::track::mask) == o2::aod::track::mask)
#define requireQualityTracksInFilter() requireTrackCutInFilter(TrackSelectionFlags::kQualityTracks)
#define requireQualityTracksITSInFilter() requireTrackCutInFilter(TrackSelectionFlags::kQualityTracksITS)
#define requirePrimaryTracksInFilter() requireTrackCutInFilter(TrackSelectionFlags::kPrimaryTracks)
#define requireInAcceptanceTracksInFilter() requireTrackCutInFilter(TrackSelectionFlags::kInAcceptanceTracks)
#define requireGlobalTrackInFilter() requireTrackCutInFilter(TrackSelectionFlags::kGlobalTrack)
#define requireGlobalTrackWoTPCClusterInFilter() requireTrackCutInFilter(TrackSelectionFlags::kGlobalTrackWoTPCCluster)
#define requireGlobalTrackWoPtEtaInFilter() requireTrackCutInFilter(TrackSelectionFlags::kGlobalTrackWoPtEta)
#define requireGlobalTrackWoDCAInFilter() requireTrackCutInFilter(TrackSelectionFlags::kGlobalTrackWoDCA)
#define requireGlobalTrackWoDCATPCClusterInFilter() requireTrackCutInFilter(TrackSelectionFlags::kGlobalTrackWoDCATPCCluster)
#define requireTrackWithinBeamPipe (nabs(o2::aod::track::x) < o2::constants::geom::XBeamPipeOuterRef)

// Columns to store track filter decisions
DECLARE_SOA_COLUMN(IsGlobalTrackSDD, isGlobalTrackSDD, uint8_t);               //!
DECLARE_SOA_COLUMN(TrackCutFlag, trackCutFlag, TrackSelectionFlags::flagtype); //! Flag with the single cut passed flagged (general selection... stil being tuned)
DECLARE_SOA_COLUMN(TrackCutFlagFb1, trackCutFlagFb1, bool);                    //! Flag with the single cut passed flagged for the first selection criteria (as general but 1 point in ITS IB)
DECLARE_SOA_COLUMN(TrackCutFlagFb2, trackCutFlagFb2, bool);                    //! Flag with the single cut passed flagged for the second selection criteria (as general but 2 point2 in ITS IB)
DECLARE_SOA_COLUMN(TrackCutFlagFb3, trackCutFlagFb3, bool);                    //! Flag with the single cut passed flagged for the third selection criteria (HF-like: global w/o tight DCA selection)
DECLARE_SOA_COLUMN(TrackCutFlagFb4, trackCutFlagFb4, bool);                    //! Flag with the single cut passed flagged for the fourth selection criteria (nuclei)
DECLARE_SOA_COLUMN(TrackCutFlagFb5, trackCutFlagFb5, bool);                    //! Flag with the single cut passed flagged for the fith selection criteria (jet validation - reduced set of cuts)

#define DECLARE_DYN_TRKSEL_COLUMN(name, getter, mask) \
  DECLARE_SOA_DYNAMIC_COLUMN(name, getter, [](TrackSelectionFlags::flagtype flags) -> bool { return TrackSelectionFlags::checkFlag(flags, mask); });

// Single selection
DECLARE_SOA_DYNAMIC_COLUMN(CheckFlag, checkFlag,
                           [](TrackSelectionFlags::flagtype flags,
                              TrackSelectionFlags::flagtype mask) -> bool { return TrackSelectionFlags::checkFlag(flags, mask); }); //! Checks the single cut
// Combo selections
DECLARE_DYN_TRKSEL_COLUMN(IsQualityTrack, isQualityTrack, TrackSelectionFlags::kQualityTracks);                                          //! Passed the combined track cut: kQualityTracks
DECLARE_DYN_TRKSEL_COLUMN(IsQualityTrackITS, isQualityTrackITS, TrackSelectionFlags::kQualityTracksITS);                                 //! Passed the combined track cut: kQualityTracksITS
DECLARE_DYN_TRKSEL_COLUMN(IsQualityTrackTPC, isQualityTrackTPC, TrackSelectionFlags::kQualityTracksTPC);                                 //! Passed the combined track cut: kQualityTracksTPC
DECLARE_DYN_TRKSEL_COLUMN(IsPrimaryTrack, isPrimaryTrack, TrackSelectionFlags::kPrimaryTracks);                                          //! Passed the combined track cut: kPrimaryTracks
DECLARE_DYN_TRKSEL_COLUMN(IsInAcceptanceTrack, isInAcceptanceTrack, TrackSelectionFlags::kInAcceptanceTracks);                           //! Passed the combined track cut: kInAcceptanceTracks
DECLARE_DYN_TRKSEL_COLUMN(IsGlobalTrack, isGlobalTrack, TrackSelectionFlags::kGlobalTrack);                                              //! Passed the combined track cut: kGlobalTrack
DECLARE_DYN_TRKSEL_COLUMN(IsGlobalTrackWoTPCCluster, isGlobalTrackWoTPCCluster, TrackSelectionFlags::kGlobalTrackWoTPCCluster);          //! Passed the combined track cut: kGlobalTrackWoTPCCluster
DECLARE_DYN_TRKSEL_COLUMN(IsGlobalTrackWoPtEta, isGlobalTrackWoPtEta, TrackSelectionFlags::kGlobalTrackWoPtEta);                         //! Passed the combined track cut: kGlobalTrackWoPtEta
DECLARE_DYN_TRKSEL_COLUMN(IsGlobalTrackWoDCA, isGlobalTrackWoDCA, TrackSelectionFlags::kGlobalTrackWoDCA);                               //! Passed the combined track cut: kGlobalTrackWoDCA
DECLARE_DYN_TRKSEL_COLUMN(IsGlobalTrackWoDCATPCCluster, isGlobalTrackWoDCATPCCluster, TrackSelectionFlags::kGlobalTrackWoDCATPCCluster); //! Passed the combined track cut: kGlobalTrackWoDCATPCCluster

#undef DECLARE_DYN_TRKSEL_COLUMN

// Single selection
DECLARE_SOA_COLUMN(PassedTrackType, passedTrackType, bool);                           //! Passed the track cut: kTrackType
DECLARE_SOA_COLUMN(PassedPtRange, passedPtRange, bool);                               //! Passed the track cut: kPtRange
DECLARE_SOA_COLUMN(PassedEtaRange, passedEtaRange, bool);                             //! Passed the track cut: kEtaRange
DECLARE_SOA_COLUMN(PassedTPCNCls, passedTPCNCls, bool);                               //! Passed the track cut: kTPCNCls
DECLARE_SOA_COLUMN(PassedTPCCrossedRows, passedTPCCrossedRows, bool);                 //! Passed the track cut: kTPCCrossedRows
DECLARE_SOA_COLUMN(PassedTPCCrossedRowsOverNCls, passedTPCCrossedRowsOverNCls, bool); //! Passed the track cut: kTPCCrossedRowsOverNCls
DECLARE_SOA_COLUMN(PassedTPCChi2NDF, passedTPCChi2NDF, bool);                         //! Passed the track cut: kTPCChi2NDF
DECLARE_SOA_COLUMN(PassedTPCRefit, passedTPCRefit, bool);                             //! Passed the track cut: kTPCRefit
DECLARE_SOA_COLUMN(PassedITSNCls, passedITSNCls, bool);                               //! Passed the track cut: kITSNCls
DECLARE_SOA_COLUMN(PassedITSChi2NDF, passedITSChi2NDF, bool);                         //! Passed the track cut: kITSChi2NDF
DECLARE_SOA_COLUMN(PassedITSRefit, passedITSRefit, bool);                             //! Passed the track cut: kITSRefit
DECLARE_SOA_COLUMN(PassedITSHits, passedITSHits, bool);                               //! Passed the track cut: kITSHits
DECLARE_SOA_COLUMN(PassedGoldenChi2, passedGoldenChi2, bool);                         //! Passed the track cut: kGoldenChi2
DECLARE_SOA_COLUMN(PassedDCAxy, passedDCAxy, bool);                                   //! Passed the track cut: kDCAxy
DECLARE_SOA_COLUMN(PassedDCAz, passedDCAz, bool);                                     //! Passed the track cut: kDCAz
DECLARE_SOA_COLUMN(PassedITSHitsFB1, passedITSHitsFB1, bool);                         //! Passed the track cut: kITSHits defined for FB1
DECLARE_SOA_COLUMN(PassedITSHitsFB2, passedITSHitsFB2, bool);                         //! Passed the track cut: kITSHits defined for FB2

// Combo selections (not yet implemented, being thinking about it)
// DECLARE_SOA_DYNAMIC_COLUMN(IsQualityTrack, isQualityTrack,
//                          [](bool passedTrackType, bool passedTPCNCls, bool passedITSChi2NDF) -> bool { return (passedTrackType && passedTPCNCls && passedITSChi2NDF); }); //! Passed the combined track cut: kQualityTracks

} // namespace track

DECLARE_SOA_TABLE(TracksDCA, "AOD", "TRACKDCA", //! DCA information for the track
                  track::DcaXY,
                  track::DcaZ);
DECLARE_SOA_TABLE(TracksDCACov, "AOD", "TRACKDCACOV",
                  track::SigmaDcaXY2,
                  track::SigmaDcaZ2); //! DCA cov. matrix information for the track

DECLARE_SOA_TABLE(TrackSelection, "AOD", "TRACKSELECTION", //! Information on the track selection decision + split dynamic information
                  track::IsGlobalTrackSDD,
                  track::TrackCutFlag,
                  track::TrackCutFlagFb1,
                  track::TrackCutFlagFb2,
                  track::TrackCutFlagFb3,
                  track::TrackCutFlagFb4,
                  track::TrackCutFlagFb5,
                  track::IsQualityTrack<track::TrackCutFlag>,
                  track::IsQualityTrackITS<track::TrackCutFlag>,
                  track::IsQualityTrackTPC<track::TrackCutFlag>,
                  track::IsPrimaryTrack<track::TrackCutFlag>,
                  track::IsInAcceptanceTrack<track::TrackCutFlag>,
                  track::IsGlobalTrack<track::TrackCutFlag>,
                  track::IsGlobalTrackWoTPCCluster<track::TrackCutFlag>,
                  track::IsGlobalTrackWoPtEta<track::TrackCutFlag>,
                  track::IsGlobalTrackWoDCA<track::TrackCutFlag>,
                  track::IsGlobalTrackWoDCATPCCluster<track::TrackCutFlag>);

DECLARE_SOA_TABLE(TrackSelectionExtension, "AOD", "TRACKSELEXTRA", //! Information on the track selections set by each Filter Bit
                  track::PassedTrackType,
                  track::PassedPtRange,
                  track::PassedEtaRange,
                  track::PassedTPCNCls,
                  track::PassedTPCCrossedRows,
                  track::PassedTPCCrossedRowsOverNCls,
                  track::PassedTPCChi2NDF,
                  track::PassedTPCRefit,
                  track::PassedITSNCls,
                  track::PassedITSChi2NDF,
                  track::PassedITSRefit,
                  track::PassedITSHits,
                  track::PassedGoldenChi2,
                  track::PassedDCAxy,
                  track::PassedDCAz,
                  track::PassedITSHitsFB1,
                  track::PassedITSHitsFB2);

DECLARE_SOA_TABLE(FwdTracksDCA, "AOD", "FWDTRACKDCA", //! DCA information for the forward track
                  fwdtrack::FwdDcaX,
                  fwdtrack::FwdDcaY);

} // namespace o2::aod

#endif // COMMON_DATAMODEL_TRACKSELECTIONTABLES_H_
