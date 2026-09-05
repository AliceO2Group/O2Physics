// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file TrackTunerCCDBObjects.h
/// \brief Declarative CCDB columns for the TrackTuner DCA / Q-over-pt calibrations.
///
/// The DCA calibration is published under a different path per data-taking period, so
/// the column declares a uniformity-value-to-path mapping rather than a single path:
/// the fetcher picks the entry whose run range contains the row's run number. This
/// replaces TrackTuner::getPathInputFileAutomaticFromCCDB(), whose run ranges these are.
/// A run matching no range is a fatal error, as it was before — there is deliberately no
/// fallback entry, since silently using another period's calibration is worse than stopping.
///
/// The ranges are the column's *default*; the whole mapping can be replaced at runtime
/// through the "ccdb:fTrackTunerDca" option, so adding a period need not be a code change.

#ifndef COMMON_DATAMODEL_TRACKTUNERCCDBOBJECTS_H_
#define COMMON_DATAMODEL_TRACKTUNERCCDBOBJECTS_H_

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>

#include <TList.h>

namespace o2::aod
{
namespace ccdbTrackTuner
{
DECLARE_SOA_CCDB_COLUMN(TrackTunerDca, trackTunerDca, TList, //!
                        "520259-529691=Users/m/mfaggin/test/inputsTrackTuner/pp2023/pass4/vsPhi;"
                        "534998-543113=Users/m/mfaggin/test/inputsTrackTuner/pp2023/pass4/vsPhi;"
                        "529397-529418=Users/m/mfaggin/test/inputsTrackTuner/PbPb2023/apass4/vsPhi;"
                        "543437-545367=Users/m/mfaggin/test/inputsTrackTuner/PbPb2023/apass4/vsPhi;"
                        "549559-558807=Users/m/mfaggin/test/inputsTrackTuner/pp2024/pass1_minBias/vsPhi;"
                        "564356-564445=Users/m/mfaggin/test/inputsTrackTuner/OO/LHC25ae;"
                        "564468-564472=Users/m/mfaggin/test/inputsTrackTuner/OO/LHC25af;"
                        "559348-559387=Users/m/mfaggin/test/inputsTrackTuner/pp2024/ppRef/polarity_positive;"
                        "559408-559456=Users/m/mfaggin/test/inputsTrackTuner/pp2024/ppRef/polarity_negative");

DECLARE_SOA_CCDB_COLUMN(TrackTunerQOverPt, trackTunerQOverPt, TList, //!
                        "Users/h/hsharma/qOverPtGraphs");
} // namespace ccdbTrackTuner

/// Uniform in the run number: one calibration per data-taking period, so the fetcher
/// resolves the path and queries once per distinct run rather than once per BC.
DECLARE_SOA_UNIFORM_TABLE(TrackTunerCCDBObjects, aod::Timestamps, o2::aod::timestamp::Timestamp,
                          aod::BCs, o2::aod::bc::RunNumber, 1, "TRKTUNERCCDB", //!
                          ccdbTrackTuner::TrackTunerDca, ccdbTrackTuner::TrackTunerQOverPt);
} // namespace o2::aod

#endif // COMMON_DATAMODEL_TRACKTUNERCCDBOBJECTS_H_
