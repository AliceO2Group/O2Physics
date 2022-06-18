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

#ifndef O2_ANALYSIS_TRACKSELECTIONTABLES_H_
#define O2_ANALYSIS_TRACKSELECTIONTABLES_H_

#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace track
{
// Columns to store the DCA to the primary vertex
DECLARE_SOA_COLUMN(DcaXY, dcaXY, float); //! Impact parameter in XY of the track to the primary vertex
DECLARE_SOA_COLUMN(DcaZ, dcaZ, float);   //! Impact parameter in Z of the track to the primary vertex

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
  static constexpr flagtype kQualityTracks = kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits;
  static constexpr flagtype kPrimaryTracks = kGoldenChi2 | kDCAxy | kDCAz;
  static constexpr flagtype kInAcceptanceTracks = kPtRange | kEtaRange;
  static constexpr flagtype kGlobalTrack = kQualityTracks | kPrimaryTracks | kInAcceptanceTracks;
  static constexpr flagtype kGlobalTrackWoPtEta = kQualityTracks | kPrimaryTracks;
};

#define requireTrackCutInFilter(mask) (aod::track::trackCutFlag & aod::track::mask) == aod::track::mask
#define requireGlobalTrackInFilter() requireTrackCutInFilter(TrackSelectionFlags::kGlobalTrack)
#define requireGlobalTrackWoPtEtaInFilter() requireTrackCutInFilter(TrackSelectionFlags::kGlobalTrackWoPtEta)
#define requireTrackWithinBeamPipe (nabs(aod::track::x) < o2::constants::geom::XBeamPipeOuterRef)

// Columns to store track filter decisions
DECLARE_SOA_COLUMN(IsGlobalTrackSDD, isGlobalTrackSDD, uint8_t);               //!
DECLARE_SOA_COLUMN(TrackCutFlag, trackCutFlag, TrackSelectionFlags::flagtype); //! Flag with the single cut passed flagged
#define DECLARE_DYN_TRKSEL_COLUMN(name, getter, mask) \
  DECLARE_SOA_DYNAMIC_COLUMN(name, getter, [](TrackSelectionFlags::flagtype flags) -> bool { return (flags & mask) == mask; });

DECLARE_DYN_TRKSEL_COLUMN(PassedTrackType, passedTrackType, TrackSelectionFlags::kTrackType);                                        //! Passed the track cut: kTrackType
DECLARE_DYN_TRKSEL_COLUMN(PassedPtRange, passedPtRange, TrackSelectionFlags::kPtRange);                                              //! Passed the track cut: kPtRange
DECLARE_DYN_TRKSEL_COLUMN(PassedEtaRange, passedEtaRange, TrackSelectionFlags::kEtaRange);                                           //! Passed the track cut: kEtaRange
DECLARE_DYN_TRKSEL_COLUMN(PassedTPCNCls, passedTPCNCls, TrackSelectionFlags::kTPCNCls);                                              //! Passed the track cut: kTPCNCls
DECLARE_DYN_TRKSEL_COLUMN(PassedTPCCrossedRows, passedTPCCrossedRows, TrackSelectionFlags::kTPCCrossedRows);                         //! Passed the track cut: kTPCCrossedRows
DECLARE_DYN_TRKSEL_COLUMN(PassedTPCCrossedRowsOverNCls, passedTPCCrossedRowsOverNCls, TrackSelectionFlags::kTPCCrossedRowsOverNCls); //! Passed the track cut: kTPCCrossedRowsOverNCls
DECLARE_DYN_TRKSEL_COLUMN(PassedTPCChi2NDF, passedTPCChi2NDF, TrackSelectionFlags::kTPCChi2NDF);                                     //! Passed the track cut: kTPCChi2NDF
DECLARE_DYN_TRKSEL_COLUMN(PassedTPCRefit, passedTPCRefit, TrackSelectionFlags::kTPCRefit);                                           //! Passed the track cut: kTPCRefit
DECLARE_DYN_TRKSEL_COLUMN(PassedITSNCls, passedITSNCls, TrackSelectionFlags::kITSNCls);                                              //! Passed the track cut: kITSNCls
DECLARE_DYN_TRKSEL_COLUMN(PassedITSChi2NDF, passedITSChi2NDF, TrackSelectionFlags::kITSChi2NDF);                                     //! Passed the track cut: kITSChi2NDF
DECLARE_DYN_TRKSEL_COLUMN(PassedITSRefit, passedITSRefit, TrackSelectionFlags::kITSRefit);                                           //! Passed the track cut: kITSRefit
DECLARE_DYN_TRKSEL_COLUMN(PassedITSHits, passedITSHits, TrackSelectionFlags::kITSHits);                                              //! Passed the track cut: kITSHits
DECLARE_DYN_TRKSEL_COLUMN(PassedGoldenChi2, passedGoldenChi2, TrackSelectionFlags::kGoldenChi2);                                     //! Passed the track cut: kGoldenChi2
DECLARE_DYN_TRKSEL_COLUMN(PassedDCAxy, passedDCAxy, TrackSelectionFlags::kDCAxy);                                                    //! Passed the track cut: kDCAxy
DECLARE_DYN_TRKSEL_COLUMN(PassedDCAz, passedDCAz, TrackSelectionFlags::kDCAz);                                                       //! Passed the track cut: kDCAz
DECLARE_DYN_TRKSEL_COLUMN(IsGlobalTrack, isGlobalTrack, TrackSelectionFlags::kGlobalTrack);                                          //! Passed the track cut: kGlobalTrack
#undef DECLARE_DYN_TRKSEL_COLUMN

} // namespace track
DECLARE_SOA_TABLE(TracksDCA, "AOD", "TRACKDCA", //! DCA information for the track
                  track::DcaXY,
                  track::DcaZ);

DECLARE_SOA_TABLE(TrackSelection, "AOD", "TRACKSELECTION", //! Information on the track selection decision + split dynamic information
                  track::IsGlobalTrackSDD,
                  track::TrackCutFlag,
                  track::PassedTrackType<track::TrackCutFlag>,
                  track::PassedPtRange<track::TrackCutFlag>,
                  track::PassedEtaRange<track::TrackCutFlag>,
                  track::PassedTPCNCls<track::TrackCutFlag>,
                  track::PassedTPCCrossedRows<track::TrackCutFlag>,
                  track::PassedTPCCrossedRowsOverNCls<track::TrackCutFlag>,
                  track::PassedTPCChi2NDF<track::TrackCutFlag>,
                  track::PassedTPCRefit<track::TrackCutFlag>,
                  track::PassedITSNCls<track::TrackCutFlag>,
                  track::PassedITSChi2NDF<track::TrackCutFlag>,
                  track::PassedITSRefit<track::TrackCutFlag>,
                  track::PassedITSHits<track::TrackCutFlag>,
                  track::PassedGoldenChi2<track::TrackCutFlag>,
                  track::PassedDCAxy<track::TrackCutFlag>,
                  track::PassedDCAz<track::TrackCutFlag>,
                  track::IsGlobalTrack<track::TrackCutFlag>);

} // namespace o2::aod

#endif // O2_ANALYSIS_TRACKSELECTIONTABLES_H_
