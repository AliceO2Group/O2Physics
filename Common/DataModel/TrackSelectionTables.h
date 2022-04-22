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
DECLARE_SOA_DYNAMIC_COLUMN(IsGlobalTrack, isGlobalTrack, //! Flag for global tracks
                           [](TrackSelectionFlags::flagtype flags) -> bool { return (flags & TrackSelectionFlags::kGlobalTrack) == TrackSelectionFlags::kGlobalTrack; });
DECLARE_SOA_COLUMN(IsGlobalTrackSDD, isGlobalTrackSDD, uint8_t);               //!
DECLARE_SOA_COLUMN(TrackCutFlag, trackCutFlag, TrackSelectionFlags::flagtype); //! Flag with the single cut passed flagged
DECLARE_SOA_DYNAMIC_COLUMN(IsTrackType, isTrackType,                           //! Passed the track cut: kTrackType
                           [](TrackSelectionFlags::flagtype flags) -> bool { return (flags & TrackSelectionFlags::kTrackType) == TrackSelectionFlags::kTrackType; });
DECLARE_SOA_DYNAMIC_COLUMN(IsPtRange, isPtRange, //! Passed the track cut: kPtRange
                           [](TrackSelectionFlags::flagtype flags) -> bool { return (flags & TrackSelectionFlags::kPtRange) == TrackSelectionFlags::kPtRange; });
DECLARE_SOA_DYNAMIC_COLUMN(IsEtaRange, isEtaRange, //! Passed the track cut: kEtaRange
                           [](TrackSelectionFlags::flagtype flags) -> bool { return (flags & TrackSelectionFlags::kEtaRange) == TrackSelectionFlags::kEtaRange; });
DECLARE_SOA_DYNAMIC_COLUMN(IsTPCNCls, isTPCNCls, //! Passed the track cut: kTPCNCls
                           [](TrackSelectionFlags::flagtype flags) -> bool { return (flags & TrackSelectionFlags::kTPCNCls) == TrackSelectionFlags::kTPCNCls; });
DECLARE_SOA_DYNAMIC_COLUMN(IsTPCCrossedRows, isTPCCrossedRows, //! Passed the track cut: kTPCCrossedRows
                           [](TrackSelectionFlags::flagtype flags) -> bool { return (flags & TrackSelectionFlags::kTPCCrossedRows) == TrackSelectionFlags::kTPCCrossedRows; });
DECLARE_SOA_DYNAMIC_COLUMN(IsTPCCrossedRowsOverNCls, isTPCCrossedRowsOverNCls, //! Passed the track cut: kTPCCrossedRowsOverNCls
                           [](TrackSelectionFlags::flagtype flags) -> bool { return (flags & TrackSelectionFlags::kTPCCrossedRowsOverNCls) == TrackSelectionFlags::kTPCCrossedRowsOverNCls; });
DECLARE_SOA_DYNAMIC_COLUMN(IsTPCChi2NDF, isTPCChi2NDF, //! Passed the track cut: kTPCChi2NDF
                           [](TrackSelectionFlags::flagtype flags) -> bool { return (flags & TrackSelectionFlags::kTPCChi2NDF) == TrackSelectionFlags::kTPCChi2NDF; });
DECLARE_SOA_DYNAMIC_COLUMN(IsTPCRefit, isTPCRefit, //! Passed the track cut: kTPCRefit
                           [](TrackSelectionFlags::flagtype flags) -> bool { return (flags & TrackSelectionFlags::kTPCRefit) == TrackSelectionFlags::kTPCRefit; });
DECLARE_SOA_DYNAMIC_COLUMN(IsITSNCls, isITSNCls, //! Passed the track cut: kITSNCls
                           [](TrackSelectionFlags::flagtype flags) -> bool { return (flags & TrackSelectionFlags::kITSNCls) == TrackSelectionFlags::kITSNCls; });
DECLARE_SOA_DYNAMIC_COLUMN(IsITSChi2NDF, isITSChi2NDF, //! Passed the track cut: kITSChi2NDF
                           [](TrackSelectionFlags::flagtype flags) -> bool { return (flags & TrackSelectionFlags::kITSChi2NDF) == TrackSelectionFlags::kITSChi2NDF; });
DECLARE_SOA_DYNAMIC_COLUMN(IsITSRefit, isITSRefit, //! Passed the track cut: kITSRefit
                           [](TrackSelectionFlags::flagtype flags) -> bool { return (flags & TrackSelectionFlags::kITSRefit) == TrackSelectionFlags::kITSRefit; });
DECLARE_SOA_DYNAMIC_COLUMN(IsITSHits, isITSHits, //! Passed the track cut: kITSHits
                           [](TrackSelectionFlags::flagtype flags) -> bool { return (flags & TrackSelectionFlags::kITSHits) == TrackSelectionFlags::kITSHits; });
DECLARE_SOA_DYNAMIC_COLUMN(IsGoldenChi2, isGoldenChi2, //! Passed the track cut: kGoldenChi2
                           [](TrackSelectionFlags::flagtype flags) -> bool { return (flags & TrackSelectionFlags::kGoldenChi2) == TrackSelectionFlags::kGoldenChi2; });
DECLARE_SOA_DYNAMIC_COLUMN(IsDCAxy, isDCAxy, //! Passed the track cut: kDCAxy
                           [](TrackSelectionFlags::flagtype flags) -> bool { return (flags & TrackSelectionFlags::kDCAxy) == TrackSelectionFlags::kDCAxy; });
DECLARE_SOA_DYNAMIC_COLUMN(IsDCAz, isDCAz, //! Passed the track cut: kDCAz
                           [](TrackSelectionFlags::flagtype flags) -> bool { return (flags & TrackSelectionFlags::kDCAz) == TrackSelectionFlags::kDCAz; });

} // namespace track
DECLARE_SOA_TABLE(TracksExtended, "AOD", "TRACKEXTENDED", //!
                  track::DcaXY,
                  track::DcaZ);

DECLARE_SOA_TABLE(TrackSelection, "AOD", "TRACKSELECTION", //!
                  track::IsGlobalTrackSDD,
                  track::TrackCutFlag,
                  track::IsGlobalTrack<track::TrackCutFlag>,
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
} // namespace o2::aod

#endif // O2_ANALYSIS_TRACKSELECTIONTABLES_H_
