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

/// \file AliasTables.h
/// \brief Table aliases
///
/// \author Gian Michele Innocenti <gian.michele.innocenti@cern.ch>, CERN
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN

#ifndef PWGHF_DATAMODEL_ALIASTABLES_H_
#define PWGHF_DATAMODEL_ALIASTABLES_H_

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <Framework/AnalysisDataModel.h>

namespace o2::aod
{
using BcFullInfos = soa::Join<aod::BCsWithTimestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;

using TracksWCov = soa::Join<Tracks, TracksCov>;
using TracksWDca = soa::Join<Tracks, TracksDCA>;
using TracksWExtra = soa::Join<Tracks, TracksExtra>;
using TracksWCovDca = soa::Join<Tracks, TracksCov, TracksDCA>;
using TracksWCovExtra = soa::Join<Tracks, TracksCov, TracksExtra>;
using TracksWDcaExtra = soa::Join<Tracks, TracksDCA, TracksExtra>;
using TracksWCovDcaExtra = soa::Join<Tracks, TracksCov, TracksDCA, TracksExtra>;

using TracksWMc = soa::Join<Tracks, McTrackLabels>;

using TracksPidEl = soa::Join<aod::pidTPCFullEl, aod::pidTOFFullEl>;
using TracksPidMu = soa::Join<aod::pidTPCFullMu, aod::pidTOFFullMu>;
using TracksPidPi = soa::Join<aod::pidTPCFullPi, aod::pidTOFFullPi>;
using TracksPidKa = soa::Join<aod::pidTPCFullKa, aod::pidTOFFullKa>;
using TracksPidPr = soa::Join<aod::pidTPCFullPr, aod::pidTOFFullPr>;
using TracksPidDe = soa::Join<aod::pidTPCFullDe, aod::pidTOFFullDe>;
using TracksPidTr = soa::Join<aod::pidTPCFullTr, aod::pidTOFFullTr>;
using TracksPidHe = soa::Join<aod::pidTPCFullHe, aod::pidTOFFullHe>;
using TracksPidAl = soa::Join<aod::pidTPCFullAl, aod::pidTOFFullAl>;

using TracksPidTinyEl = soa::Join<aod::pidTPCEl, aod::pidTOFEl>;
using TracksPidTinyMu = soa::Join<aod::pidTPCMu, aod::pidTOFMu>;
using TracksPidTinyPi = soa::Join<aod::pidTPCPi, aod::pidTOFPi>;
using TracksPidTinyKa = soa::Join<aod::pidTPCKa, aod::pidTOFKa>;
using TracksPidTinyPr = soa::Join<aod::pidTPCPr, aod::pidTOFPr>;
using TracksPidTinyDe = soa::Join<aod::pidTPCDe, aod::pidTOFDe>;
using TracksPidTinyTr = soa::Join<aod::pidTPCTr, aod::pidTOFTr>;
using TracksPidTinyHe = soa::Join<aod::pidTPCHe, aod::pidTOFHe>;
using TracksPidTinyAl = soa::Join<aod::pidTPCAl, aod::pidTOFAl>;
} // namespace o2::aod

#endif // PWGHF_DATAMODEL_ALIASTABLES_H_
