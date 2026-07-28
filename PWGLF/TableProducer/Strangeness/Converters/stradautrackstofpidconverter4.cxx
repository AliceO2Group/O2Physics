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
//
/// \file stradautrackstofpidconverter4.cxx
/// \brief Converts DauTrackTOFPID_002 into DauTrackTOFPID_003
///
/// \author David Dobrigkeit Chinellato <david.dobrigkeit.chinellato@cern.ch>, Austrian Academy of Sciences & MBI
/// \author Romain Schotter <romain.schotter@cern.ch>, Austrian Academy of Sciences & MBI

#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/runDataProcessing.h>

using namespace o2;
using namespace o2::framework;

// converts DauTrackTOFPIDs_002 to _003
struct stradautrackstofpidconverter4 {
  Produces<aod::DauTrackTOFPIDs_003> dautracktofpids;

  void process(aod::DauTrackTOFPIDs_002 const& dauTracks, soa::Join<aod::StraCollisions, aod::StraStamps> const& collisions)
  {
    // create new TOFPIDs
    dautracktofpids.reserve(dauTracks.size());
    for (const auto& dauTrack : dauTracks) {
      uint64_t bc = 0;
      if (dauTrack.straCollisionId() >= 0) {
        auto collision = collisions.rawIteratorAt(dauTrack.straCollisionId());
        bc = collision.globalBC();
      }
      dautracktofpids(
        bc,
        dauTrack.dauTrackExtraId(),
        dauTrack.tofSignal(),
        dauTrack.tofEvTime(),
        dauTrack.tofEvTimeErr(),
        dauTrack.length(),
        dauTrack.tofExpMom());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<stradautrackstofpidconverter4>(cfgc)};
}
