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

/// \file UDFwdTracksExtraConverter.cxx
/// \brief Converts UDFwdTracksExtra table from version 000 to 001

/// This task allows for the conversion of the UDFwdTracksExtra table from the version 000,
/// that is before the introduction of global tracks (and so contains only MCH-MID and
/// MCH only tracks), to the version 001, that includes global tracks
/// Global tracks introduced here: https://github.com/AliceO2Group/O2Physics/commit/72f673611ddcd7c39978787dbed9f77e6e7c3d6a)

/// executable name o2-analysis-ud-fwd-tracks-extra-converter

/// \author Andrea Giovanni Riffero <andrea.giovanni.riffero@cern.ch>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "PWGUD/DataModel/UDTables.h"

using namespace o2;
using namespace o2::framework;

// Converts UDFwdTracksExtra for version 000 to 001
// v_000 only MID-MCH and MCH only tracks
// v_001 global tracks
struct UDFwdTracksExtraConverter {
  Produces<o2::aod::UDFwdTracksExtra_001> udFwdTracksExtra_001;

  void process(o2::aod::UDFwdTracksExtra_000 const& tracks)
  {
    int trkType = 3; // trackType of MCH-MID tracks is 3

    for (const auto& track : tracks) {

      if (track.chi2MatchMCHMID() > 0)
        trkType = 3; // trackType of MCH-MID tracks is 3
      if (track.chi2MatchMCHMID() < 0)
        trkType = 4; // trackType of MCH only tracks is 4

      udFwdTracksExtra_001(trkType,
                           track.nClusters(),
                           track.pDca(),
                           track.rAtAbsorberEnd(),
                           track.chi2(),
                           track.chi2MatchMCHMID(),
                           0.0f, // dummy mchmftChi2, not available in version 000
                           track.mchBitMap(),
                           track.midBitMap(),
                           track.midBoards());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UDFwdTracksExtraConverter>(cfgc),
  };
}
