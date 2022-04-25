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
  static constexpr flagtype kGlobalTrack = kTrackType | kPtRange | kEtaRange | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits | kGoldenChi2 | kDCAxy | kDCAz;
};

// Columns to store track filter decisions
DECLARE_SOA_COLUMN(IsGlobalTrackSDD, isGlobalTrackSDD, uint8_t);               //!
DECLARE_SOA_COLUMN(TrackCutFlag, trackCutFlag, TrackSelectionFlags::flagtype); //! Flag with the single cut passed flagged
#define DECLARE_DYN_TRKSEL_COLUMN(name, getter, mask) \
  DECLARE_SOA_DYNAMIC_COLUMN(name, getter, [](TrackSelectionFlags::flagtype flags) -> bool { return (flags & mask) == mask; });

DECLARE_SOA_EXPRESSION_COLUMN(IsGlobalTrack, isGlobalTrack, bool, //! Flag for global tracks
                              (aod::track::trackCutFlag& TrackSelectionFlags::kGlobalTrack) == TrackSelectionFlags::kGlobalTrack);

DECLARE_DYN_TRKSEL_COLUMN(IsTrackType, isTrackType, TrackSelectionFlags::kTrackType);                                        //! Passed the track cut: kTrackType
DECLARE_DYN_TRKSEL_COLUMN(IsPtRange, isPtRange, TrackSelectionFlags::kPtRange);                                              //! Passed the track cut: kPtRange
DECLARE_DYN_TRKSEL_COLUMN(IsEtaRange, isEtaRange, TrackSelectionFlags::kEtaRange);                                           //! Passed the track cut: kEtaRange
DECLARE_DYN_TRKSEL_COLUMN(IsTPCNCls, isTPCNCls, TrackSelectionFlags::kTPCNCls);                                              //! Passed the track cut: kTPCNCls
DECLARE_DYN_TRKSEL_COLUMN(IsTPCCrossedRows, isTPCCrossedRows, TrackSelectionFlags::kTPCCrossedRows);                         //! Passed the track cut: kTPCCrossedRows
DECLARE_DYN_TRKSEL_COLUMN(IsTPCCrossedRowsOverNCls, isTPCCrossedRowsOverNCls, TrackSelectionFlags::kTPCCrossedRowsOverNCls); //! Passed the track cut: kTPCCrossedRowsOverNCls
DECLARE_DYN_TRKSEL_COLUMN(IsTPCChi2NDF, isTPCChi2NDF, TrackSelectionFlags::kTPCChi2NDF);                                     //! Passed the track cut: kTPCChi2NDF
DECLARE_DYN_TRKSEL_COLUMN(IsTPCRefit, isTPCRefit, TrackSelectionFlags::kTPCRefit);                                           //! Passed the track cut: kTPCRefit
DECLARE_DYN_TRKSEL_COLUMN(IsITSNCls, isITSNCls, TrackSelectionFlags::kITSNCls);                                              //! Passed the track cut: kITSNCls
DECLARE_DYN_TRKSEL_COLUMN(IsITSChi2NDF, isITSChi2NDF, TrackSelectionFlags::kITSChi2NDF);                                     //! Passed the track cut: kITSChi2NDF
DECLARE_DYN_TRKSEL_COLUMN(IsITSRefit, isITSRefit, TrackSelectionFlags::kITSRefit);                                           //! Passed the track cut: kITSRefit
DECLARE_DYN_TRKSEL_COLUMN(IsITSHits, isITSHits, TrackSelectionFlags::kITSHits);                                              //! Passed the track cut: kITSHits
DECLARE_DYN_TRKSEL_COLUMN(IsGoldenChi2, isGoldenChi2, TrackSelectionFlags::kGoldenChi2);                                     //! Passed the track cut: kGoldenChi2
DECLARE_DYN_TRKSEL_COLUMN(IsDCAxy, isDCAxy, TrackSelectionFlags::kDCAxy);                                                    //! Passed the track cut: kDCAxy
DECLARE_DYN_TRKSEL_COLUMN(IsDCAz, isDCAz, TrackSelectionFlags::kDCAz);                                                       //! Passed the track cut: kDCAz
#undef DECLARE_DYN_TRKSEL_COLUMN

} // namespace track
DECLARE_SOA_TABLE(TracksExtended, "AOD", "TRACKEXTENDED", //!
                  track::DcaXY,
                  track::DcaZ);

DECLARE_SOA_TABLE(TrackSelectionStore, "AOD", "TRACKSELECTION", //! Stored information on the track selection decision + split dynamic information
                  track::IsGlobalTrackSDD,
                  track::TrackCutFlag,
                  track::IsTrackType<track::TrackCutFlag>,
                  track::IsPtRange<track::TrackCutFlag>,
                  track::IsEtaRange<track::TrackCutFlag>,
                  track::IsTPCNCls<track::TrackCutFlag>,
                  track::IsTPCCrossedRows<track::TrackCutFlag>,
                  track::IsTPCCrossedRowsOverNCls<track::TrackCutFlag>,
                  track::IsTPCChi2NDF<track::TrackCutFlag>,
                  track::IsTPCRefit<track::TrackCutFlag>,
                  track::IsITSNCls<track::TrackCutFlag>,
                  track::IsITSChi2NDF<track::TrackCutFlag>,
                  track::IsITSRefit<track::TrackCutFlag>,
                  track::IsITSHits<track::TrackCutFlag>,
                  track::IsGoldenChi2<track::TrackCutFlag>,
                  track::IsDCAxy<track::TrackCutFlag>,
                  track::IsDCAz<track::TrackCutFlag>);

DECLARE_SOA_EXTENDED_TABLE(TrackSelection, TrackSelectionStore, "TRACKSELECTION", //! Split information on the track selection decision
                           track::IsGlobalTrack);

} // namespace o2::aod

#endif // O2_ANALYSIS_TRACKSELECTIONTABLES_H_
