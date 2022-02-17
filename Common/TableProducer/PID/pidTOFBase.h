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

///
/// \file   pidTOFBase.h
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Base to build tasks for TOF PID tasks.
///

#include "Framework/AnalysisTask.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Common/Core/PID/PIDTOF.h"
#include "Common/Core/PID/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include <CCDB/BasicCCDBManager.h>
#include "TableHelper.h"

using namespace o2;
using namespace o2::track;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// Task to produce the TOF signal from the trackTime information
struct tofSignal {
  o2::framework::Produces<o2::aod::TOFSignal> table;
  bool enableTable = false;

  void init(o2::framework::InitContext& initContext)
  {
    // Checking that the table is requested in the workflow and enabling it
    enableTable = isTableRequiredInWorkflow(initContext, "TOFSignal");
  }
  using Trks = o2::soa::Join<aod::Tracks, aod::TracksExtra>;
  void process(Trks const& tracks)
  {
    if (!enableTable) {
      return;
    }
    table.reserve(tracks.size());
    for (auto& t : tracks) {
      table(o2::pid::tof::TOFSignal<Trks::iterator>::GetTOFSignal(t));
    }
  }
};

/// Selection criteria for tracks used for TOF event time
template <typename trackType>
bool filterForTOFEventTime(const trackType& tr)
{
  return (tr.hasTOF() && tr.p() > 0.5f && tr.p() < 2.f && tr.trackType() == o2::aod::track::TrackTypeEnum::Track);
} // accept all

/// Specialization of TOF event time maker
template <typename trackType,
          bool (*trackFilter)(const trackType&),
          template <typename T, o2::track::PID::ID> typename response,
          typename trackTypeContainer,
          typename responseParametersType>
o2::tof::eventTimeContainer evTimeMakerForTracks(const trackTypeContainer& tracks,
                                                 const responseParametersType& responseParameters)
{
  return o2::tof::evTimeMakerFromParam<trackTypeContainer, trackType, trackFilter, response, responseParametersType>(tracks, responseParameters);
}
