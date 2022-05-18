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
/// \file   pidTOF.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Task to produce PID tables for TOF split for each particle with only the Nsigma information.
///         The event time maker can be used to produce event TOF times.
///         Only the tables for the mass hypotheses requested are filled, the others are sent empty.
///         QA histograms for the TOF PID can be produced by adding `--add-qa 1` to the workflow
///

// O2 includes
#include "Framework/AnalysisTask.h"
#include "ReconstructionDataFormats/Track.h"
#include <CCDB/BasicCCDBManager.h>
#include "Common/Core/PID/PIDResponse.h"
#include "Common/Core/PID/PIDTOF.h"
#include "TableHelper.h"
#include "pidTOFBase.h"
#include "DPG/Tasks/qaPIDTOF.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::pid;
using namespace o2::framework::expressions;
using namespace o2::track;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  std::vector<ConfigParamSpec> options{{"add-qa", VariantType::Int, 0, {"Produce TOF PID QA histograms"}}};
  std::swap(workflowOptions, options);
}

#include "Framework/runDataProcessing.h"

/// Task to produce the response table
struct tofPid {
  // Tables to produce
  Produces<o2::aod::pidTOFEl> tablePIDEl;
  Produces<o2::aod::pidTOFMu> tablePIDMu;
  Produces<o2::aod::pidTOFPi> tablePIDPi;
  Produces<o2::aod::pidTOFKa> tablePIDKa;
  Produces<o2::aod::pidTOFPr> tablePIDPr;
  Produces<o2::aod::pidTOFDe> tablePIDDe;
  Produces<o2::aod::pidTOFTr> tablePIDTr;
  Produces<o2::aod::pidTOFHe> tablePIDHe;
  Produces<o2::aod::pidTOFAl> tablePIDAl;
  // Detector response and input parameters
  DetectorResponse response;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> paramfile{"param-file", "", "Path to the parametrization object, if emtpy the parametrization is not taken from file"};
  Configurable<std::string> sigmaname{"param-sigma", "TOFReso", "Name of the parametrization for the expected sigma, used in both file and CCDB mode"};
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPath{"ccdbPath", "Analysis/PID/TOF", "Path of the TOF parametrization on the CCDB"};
  Configurable<long> timestamp{"ccdb-timestamp", -1, "timestamp of the object"};
  // Configuration flags to include and exclude particle hypotheses
  Configurable<int> pidEl{"pid-el", -1, {"Produce PID information for the Electron mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidMu{"pid-mu", -1, {"Produce PID information for the Muon mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidPi{"pid-pi", -1, {"Produce PID information for the Pion mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidKa{"pid-ka", -1, {"Produce PID information for the Kaon mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidPr{"pid-pr", -1, {"Produce PID information for the Proton mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidDe{"pid-de", -1, {"Produce PID information for the Deuterons mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidTr{"pid-tr", -1, {"Produce PID information for the Triton mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidHe{"pid-he", -1, {"Produce PID information for the Helium3 mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidAl{"pid-al", -1, {"Produce PID information for the Alpha mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};

  void init(o2::framework::InitContext& initContext)
  {
    // Check that both processes are not enabled
    if (doprocessNoEvTime == true && doprocessEvTime == true) {
      LOGF(fatal, "Cannot enable processNoEvTime and processEvTime at the same time. Please choose one.");
    }
    // Checking the tables are requested in the workflow and enabling them
    auto enableFlag = [&](const std::string particle, Configurable<int>& flag) {
      const std::string table = "pidTOF" + particle;
      if (isTableRequiredInWorkflow(initContext, table)) {
        if (flag < 0) {
          flag.value = 1;
          LOG(info) << "Auto-enabling table: " + table;
        } else if (flag > 0) {
          flag.value = 1;
          LOG(info) << "Table enabled: " + table;
        } else {
          LOG(info) << "Table disabled: " + table;
        }
      }
    };

    enableFlag("El", pidEl);
    enableFlag("Mu", pidMu);
    enableFlag("Pi", pidPi);
    enableFlag("Ka", pidKa);
    enableFlag("Pr", pidPr);
    enableFlag("De", pidDe);
    enableFlag("Tr", pidTr);
    enableFlag("He", pidHe);
    enableFlag("Al", pidAl);

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
  }

  using TrksEvTime = soa::Join<aod::Tracks, aod::TracksExtra, aod::TOFSignal, aod::TOFEvTime>;
  template <o2::track::PID::ID pid>
  using ResponseImplementationEvTime = o2::pid::tof::ExpTimes<TrksEvTime::iterator, pid>;
  void processEvTime(TrksEvTime const& tracks, aod::Collisions const&)
  {
    constexpr auto responseEl = ResponseImplementationEvTime<PID::Electron>();
    constexpr auto responseMu = ResponseImplementationEvTime<PID::Muon>();
    constexpr auto responsePi = ResponseImplementationEvTime<PID::Pion>();
    constexpr auto responseKa = ResponseImplementationEvTime<PID::Kaon>();
    constexpr auto responsePr = ResponseImplementationEvTime<PID::Proton>();
    constexpr auto responseDe = ResponseImplementationEvTime<PID::Deuteron>();
    constexpr auto responseTr = ResponseImplementationEvTime<PID::Triton>();
    constexpr auto responseHe = ResponseImplementationEvTime<PID::Helium3>();
    constexpr auto responseAl = ResponseImplementationEvTime<PID::Alpha>();

    auto reserveTable = [&tracks](const Configurable<int>& flag, auto& table) {
      if (flag.value != 1) {
        return;
      }
      table.reserve(tracks.size());
    };

    reserveTable(pidEl, tablePIDEl);
    reserveTable(pidMu, tablePIDMu);
    reserveTable(pidPi, tablePIDPi);
    reserveTable(pidKa, tablePIDKa);
    reserveTable(pidPr, tablePIDPr);
    reserveTable(pidDe, tablePIDDe);
    reserveTable(pidTr, tablePIDTr);
    reserveTable(pidHe, tablePIDHe);
    reserveTable(pidAl, tablePIDAl);

    for (auto const& trk : tracks) { // Loop on collisions
      if (!trk.has_collision()) {    // Track was not assigned, cannot compute event time
        auto fillEmptyTable = [](const Configurable<int>& flag, auto& table) {
          if (flag.value != 1) {
            return;
          }
          aod::pidutils::packInTable<aod::pidtof_tiny::binning>(-999.f,
                                                                table);
        };

        fillEmptyTable(pidEl, tablePIDEl);
        fillEmptyTable(pidMu, tablePIDMu);
        fillEmptyTable(pidPi, tablePIDPi);
        fillEmptyTable(pidKa, tablePIDKa);
        fillEmptyTable(pidPr, tablePIDPr);
        fillEmptyTable(pidDe, tablePIDDe);
        fillEmptyTable(pidTr, tablePIDTr);
        fillEmptyTable(pidHe, tablePIDHe);
        fillEmptyTable(pidAl, tablePIDAl);

        continue;
      }

      // Check and fill enabled tables
      auto makeTable = [&trk, this](const Configurable<int>& flag, auto& table, const auto& responsePID) {
        if (flag.value == 1) {
          aod::pidutils::packInTable<aod::pidtof_tiny::binning>(responsePID.GetSeparation(response, trk, trk.tofEvTime(), trk.tofEvTimeErr()),
                                                                table);
        }
      };

      makeTable(pidEl, tablePIDEl, responseEl);
      makeTable(pidMu, tablePIDMu, responseMu);
      makeTable(pidPi, tablePIDPi, responsePi);
      makeTable(pidKa, tablePIDKa, responseKa);
      makeTable(pidPr, tablePIDPr, responsePr);
      makeTable(pidDe, tablePIDDe, responseDe);
      makeTable(pidTr, tablePIDTr, responseTr);
      makeTable(pidHe, tablePIDHe, responseHe);
      makeTable(pidAl, tablePIDAl, responseAl);
    }
  }

  PROCESS_SWITCH(tofPid, processEvTime, "Produce TOF response with TOF event time", false);

  using Trks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TOFSignal>;
  template <o2::track::PID::ID pid>
  using ResponseImplementation = o2::pid::tof::ExpTimes<Trks::iterator, pid>;
  void processNoEvTime(Trks const& tracks, aod::Collisions const&)
  {
    constexpr auto responseEl = ResponseImplementation<PID::Electron>();
    constexpr auto responseMu = ResponseImplementation<PID::Muon>();
    constexpr auto responsePi = ResponseImplementation<PID::Pion>();
    constexpr auto responseKa = ResponseImplementation<PID::Kaon>();
    constexpr auto responsePr = ResponseImplementation<PID::Proton>();
    constexpr auto responseDe = ResponseImplementation<PID::Deuteron>();
    constexpr auto responseTr = ResponseImplementation<PID::Triton>();
    constexpr auto responseHe = ResponseImplementation<PID::Helium3>();
    constexpr auto responseAl = ResponseImplementation<PID::Alpha>();

    // Check and fill enabled tables
    auto makeTable = [&tracks](const Configurable<int>& flag, auto& table, const DetectorResponse& response, const auto& responsePID) {
      if (flag.value == 1) {
        // Prepare memory for enabled tables
        table.reserve(tracks.size());
        for (auto const& trk : tracks) { // Loop on Tracks
          aod::pidutils::packInTable<aod::pidtof_tiny::binning>(responsePID.GetSeparation(response, trk),
                                                                table);
        }
      }
    };
    makeTable(pidEl, tablePIDEl, response, responseEl);
    makeTable(pidMu, tablePIDMu, response, responseMu);
    makeTable(pidPi, tablePIDPi, response, responsePi);
    makeTable(pidKa, tablePIDKa, response, responseKa);
    makeTable(pidPr, tablePIDPr, response, responsePr);
    makeTable(pidDe, tablePIDDe, response, responseDe);
    makeTable(pidTr, tablePIDTr, response, responseTr);
    makeTable(pidHe, tablePIDHe, response, responseHe);
    makeTable(pidAl, tablePIDAl, response, responseAl);
  }

  PROCESS_SWITCH(tofPid, processNoEvTime, "Produce TOF response without TOF event time, standard for Run 2", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  auto workflow = WorkflowSpec{adaptAnalysisTask<tofPid>(cfgc)};
  if (cfgc.options().get<int>("add-qa")) {
    workflow.push_back(adaptAnalysisTask<tofPidQa>(cfgc));
  }
  return workflow;
}
