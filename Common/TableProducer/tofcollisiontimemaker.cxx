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
/// \file   tofcollisiontimemaker.cxx
/// \author Francesca Ercolessi francesca.ercolessi@cern.ch
/// \author Francesco Noferini francesco.noferini@cern.ch
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Task to produce the Event Time with the TOF using routines in the O2
///         This task will be merged in the workflow for the TOF-PID as it is one of its essential components
///

// O2 includes
#include "Framework/AnalysisTask.h"
#include "TOFReconstruction/EventTimeMaker.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "ReconstructionDataFormats/Track.h"
#include <CCDB/BasicCCDBManager.h>
#include "Common/Core/PID/PIDResponse.h"
#include "Common/Core/PID/PIDTOF.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::pid;
using namespace o2::framework::expressions;
using namespace o2::track;

namespace o2::aod
{
namespace tofeventtime
{
DECLARE_SOA_COLUMN(TOFEvTime, tofEvTime, float);       //! TOF event time
DECLARE_SOA_COLUMN(TOFEvTimeErr, tofEvTimeErr, float); //! TOF event time error
} // namespace tofeventtime

DECLARE_SOA_TABLE(TOFEventTime, "AOD", "TOFSignal", //! Table of the TOF event time
                  tofeventtime::TOFEvTime,
                  tofeventtime::TOFEvTimeErr);
} // namespace o2::aod

struct tofCollisionTimeMaker { /// Task that produces the TOF collision time
  Produces<o2::aod::TOFEventTime> table;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> paramfile{"param-file", "", "Path to the parametrization object, if emtpy the parametrization is not taken from file"};
  Configurable<std::string> sigmaname{"param-sigma", "TOFReso", "Name of the parametrization for the expected sigma, used in both file and CCDB mode"};
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPath{"ccdbPath", "Analysis/PID/TOF", "Path of the TOF parametrization on the CCDB"};
  Configurable<long> timestamp{"ccdb-timestamp", -1, "timestamp of the object"};
  Configurable<bool> addqa{"add-qa", 0, "Flag to add QA"};
  DetectorResponse response;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::QAObject};
  void init(o2::framework::InitContext& initContext)
  {
    // Getting the parametrization parameters
    ccdb->setURL(url.value);
    ccdb->setTimestamp(timestamp.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    // Not later than now objects
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    //
    const std::vector<float> p = {0.008, 0.008, 0.002, 40.0};
    response.SetParameters(DetectorResponse::kSigma, p);
    const std::string fname = paramfile.value;
    if (!fname.empty()) { // Loading the parametrization from file
      LOG(info) << "Loading exp. sigma parametrization from file" << fname << ", using param: " << sigmaname.value;
      response.LoadParamFromFile(fname.data(), sigmaname.value, DetectorResponse::kSigma);
    } else { // Loading it from CCDB
      std::string path = ccdbPath.value + "/" + sigmaname.value;
      LOG(info) << "Loading exp. sigma parametrization from CCDB, using path: " << path << " for timestamp " << timestamp.value;
      response.LoadParam(DetectorResponse::kSigma, ccdb->getForTimeStamp<Parametrization>(path, timestamp.value));
    }
    if (!addqa) {
      return;
    }
    histos.add("eventTime", "", kTH1F, {{1000, -10, 10}});
    histos.add("eventTimeErr", "", kTH1F, {{1000, -10, 10}});
  }
  using Trks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TOFSignal>;
  template <o2::track::PID::ID pid>
  using ResponseImplementation = o2::pid::tof::ExpTimes<Trks::iterator, pid>;
  void process(aod::Collision&, Trks const& tracks)
  {
    constexpr auto responsePi = ResponseImplementation<PID::Pion>();
    constexpr auto responseKa = ResponseImplementation<PID::Kaon>();
    constexpr auto responsePr = ResponseImplementation<PID::Proton>();

    table.reserve(tracks.size());
    auto evTime = o2::tof::evTimeMaker(tracks);
    for (auto& t : tracks) {
      // Expected times
      const float expPi = responsePi.GetExpectedSignal(t);
      const float expKa = responseKa.GetExpectedSignal(t);
      const float expPr = responsePr.GetExpectedSignal(t);
      // Expected resolutions
      const float expResoPi = responsePi.GetExpectedSigma(response, t);
      const float expResoKa = responseKa.GetExpectedSigma(response, t);
      const float expResoPr = responsePr.GetExpectedSigma(response, t);

      table(evTime.eventTime, evTime.eventTimeError);
      if (!addqa) {
        continue;
      }
      histos.fill(HIST("eventTime"), evTime.eventTime);
      histos.fill(HIST("eventTimeErr"), evTime.eventTimeError);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  auto workflow = WorkflowSpec{adaptAnalysisTask<tofCollisionTimeMaker>(cfgc)};
  return workflow;
}
