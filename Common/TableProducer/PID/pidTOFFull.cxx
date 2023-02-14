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
/// \file   pidTOFFull.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Task to produce PID tables for TOF split for each particle.
///         Only the tables for the mass hypotheses requested are filled, the others are sent empty.
///

// O2 includes
#include <CCDB/BasicCCDBManager.h>
#include "TOFBase/EventTimeMaker.h"
#include "Framework/AnalysisTask.h"
#include "ReconstructionDataFormats/Track.h"

// O2Physics includes
#include "TableHelper.h"
#include "pidTOFBase.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::pid;
using namespace o2::framework::expressions;
using namespace o2::track;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  std::vector<ConfigParamSpec> options{{"add-qa", VariantType::Int, 0, {"Legacy. No effect."}}};
  std::swap(workflowOptions, options);
}

#include "Framework/runDataProcessing.h"

/// Task to produce the response table
struct tofPidFull {
  // Tables to produce
  Produces<o2::aod::pidTOFFullEl> tablePIDEl;
  Produces<o2::aod::pidTOFFullMu> tablePIDMu;
  Produces<o2::aod::pidTOFFullPi> tablePIDPi;
  Produces<o2::aod::pidTOFFullKa> tablePIDKa;
  Produces<o2::aod::pidTOFFullPr> tablePIDPr;
  Produces<o2::aod::pidTOFFullDe> tablePIDDe;
  Produces<o2::aod::pidTOFFullTr> tablePIDTr;
  Produces<o2::aod::pidTOFFullHe> tablePIDHe;
  Produces<o2::aod::pidTOFFullAl> tablePIDAl;
  // Detector response parameters
  o2::pid::tof::TOFResoParams mRespParams;
  o2::pid::tof::TOFResoParamsV2 mRespParamsV2;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> paramfile{"param-file", "", "Path to the parametrization object, if empty the parametrization is not taken from file"};
  Configurable<std::string> sigmaname{"param-sigma", "TOFResoParams", "Name of the parametrization for the expected sigma, used in both file and CCDB mode"};
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPath{"ccdbPath", "Analysis/PID/TOF", "Path of the TOF parametrization on the CCDB"};
  Configurable<std::string> passName{"passName", "", "Name of the pass inside of the CCDB parameter collection. If empty, the automatically deceted from metadata (to be implemented!!!)"};
  Configurable<int64_t> timestamp{"ccdb-timestamp", -1, "timestamp of the object"};

  Configurable<bool> enableTimeDependentResponse{"enableTimeDependentResponse", false, "Flag to use the collision timestamp to fetch the PID Response"};
  Configurable<bool> useParamCollection{"useParamCollection", false, "Flag to use the parameter collection instead of the legacy parameter distribution chain"};
  Configurable<bool> fatalOnPassNotAvailable{"fatalOnPassNotAvailable", true, "Flag to throw a fatal if the pass is not available in the retrieved CCDB object"};
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
  // Running variables
  std::string parametrizationPath = "";

  void init(o2::framework::InitContext& initContext)
  {
    if (doprocessWSlice == true && doprocessWoSlice == true && doprocessWoSliceDev == true) {
      LOGF(fatal, "Cannot enable processWoSlice and processWSlice and doprocessWoSliceDev at the same time. Please choose one.");
    }
    if (doprocessWSlice == false && doprocessWoSlice == false && doprocessWoSliceDev == false) {
      LOGF(fatal, "Cannot run without any of processWoSlice and processWSlice and doprocessWoSliceDev enabled. Please choose one.");
    }

    // Checking the tables are requested in the workflow and enabling them
    auto enableFlag = [&](const std::string particle, Configurable<int>& flag) {
      enableFlagIfTableRequired(initContext, "pidTOFFull" + particle, flag);
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
    const std::string fname = paramfile.value;
    if (!fname.empty()) { // Loading the parametrization from file
      LOG(info) << "Loading exp. sigma parametrization from file" << fname << ", using param: " << sigmaname.value;
      if (useParamCollection) {
        mRespParamsV2.loadParamFromFile(fname.data(), sigmaname.value);
      } else {
        mRespParams.LoadParamFromFile(fname.data(), sigmaname.value);
      }
    } else { // Loading it from CCDB
      parametrizationPath = ccdbPath.value + "/" + sigmaname.value;
      if (!enableTimeDependentResponse) {
        LOG(info) << "Loading exp. sigma parametrization from CCDB, using path: '" << parametrizationPath << "' for timestamp " << timestamp.value;
        if (useParamCollection) {
          // TODO: implement the automatic pass name detection from metadata
          if (passName.value == "") {
            passName.value = "unanchored"; // temporary default
            LOG(warning) << "Passed autodetect mode for pass, not implemented yet, waiting for metadata. Taking '" << passName.value << "'";
          }
          LOG(info) << "Using parameter collection, starting from pass '" << passName.value << "'";
          o2::tof::ParameterCollection* paramCollection = ccdb->getForTimeStamp<o2::tof::ParameterCollection>(parametrizationPath, timestamp.value);
          paramCollection->print();
          if (!paramCollection->retrieveParameters(mRespParamsV2, passName.value)) {
            if (fatalOnPassNotAvailable) {
              LOGF(fatal, "Pass '%s' not available in the retrieved CCDB object", passName.value.data());
            } else {
              LOGF(warning, "Pass '%s' not available in the retrieved CCDB object", passName.value.data());
            }
          }
          mRespParamsV2.print();
        } else {
          mRespParams.SetParameters(ccdb->getForTimeStamp<o2::pid::tof::TOFResoParams>(parametrizationPath, timestamp.value));
          mRespParams.Print();
        }
      }
    }
  }

  using Trks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TOFSignal, aod::TOFEvTime, aod::pidEvTimeFlags>;
  // Define slice per collision
  Preslice<Trks> perCollision = aod::track::collisionId;
  template <o2::track::PID::ID pid>
  using ResponseImplementation = o2::pid::tof::ExpTimes<Trks::iterator, pid>;
  void processWSlice(Trks const& tracks, aod::Collisions const&, aod::BCsWithTimestamps const&)
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

    int lastCollisionId = -1;          // Last collision ID analysed
    for (auto const& track : tracks) { // Loop on all tracks
      if (!track.has_collision()) {    // Track was not assigned, cannot compute NSigma (no event time) -> filling with empty table
        auto makeTableEmpty = [&](const Configurable<int>& flag, auto& table) {
          if (flag.value != 1) {
            return;
          }
          table(-999.f, -999.f);
        };

        makeTableEmpty(pidEl, tablePIDEl);
        makeTableEmpty(pidMu, tablePIDMu);
        makeTableEmpty(pidPi, tablePIDPi);
        makeTableEmpty(pidKa, tablePIDKa);
        makeTableEmpty(pidPr, tablePIDPr);
        makeTableEmpty(pidDe, tablePIDDe);
        makeTableEmpty(pidTr, tablePIDTr);
        makeTableEmpty(pidHe, tablePIDHe);
        makeTableEmpty(pidAl, tablePIDAl);

        continue;
      }

      if (track.collisionId() == lastCollisionId) { // Tracks from last collision already processed
        continue;
      }

      // Fill new table for the tracks in a collision
      lastCollisionId = track.collisionId(); // Cache last collision ID
      timestamp.value = track.collision().bc_as<aod::BCsWithTimestamps>().timestamp();
      if (enableTimeDependentResponse) {
        LOG(debug) << "Updating parametrization from path '" << parametrizationPath << "' and timestamp " << timestamp.value;
        if (useParamCollection) {
          if (!ccdb->getForTimeStamp<o2::tof::ParameterCollection>(parametrizationPath, timestamp.value)->retrieveParameters(mRespParamsV2, passName.value)) {
            if (fatalOnPassNotAvailable) {
              LOGF(fatal, "Pass '%s' not available in the retrieved CCDB object", passName.value.data());
            } else {
              LOGF(warning, "Pass '%s' not available in the retrieved CCDB object", passName.value.data());
            }
          }
        } else {
          mRespParams.SetParameters(ccdb->getForTimeStamp<o2::pid::tof::TOFResoParams>(parametrizationPath, timestamp));
        }
      }

      const auto& tracksInCollision = tracks.sliceBy(perCollision, lastCollisionId);
      for (auto const& trkInColl : tracksInCollision) { // Loop on tracks
        // Check and fill enabled tables
        auto makeTable = [&trkInColl, this](const Configurable<int>& flag, auto& table, const auto& responsePID) {
          if (flag.value != 1) {
            return;
          }
          if (useParamCollection) {
            table(responsePID.GetExpectedSigma(mRespParamsV2, trkInColl),
                  responsePID.GetSeparation(mRespParamsV2, trkInColl));
          } else {
            table(responsePID.GetExpectedSigma(mRespParams, trkInColl),
                  responsePID.GetSeparation(mRespParams, trkInColl));
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
  }
  PROCESS_SWITCH(tofPidFull, processWSlice, "Process with track slices", true);

  void processWoSlice(Trks const& tracks, aod::Collisions const&, aod::BCsWithTimestamps const&)
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

    int lastCollisionId = -1;          // Last collision ID analysed
    for (auto const& track : tracks) { // Loop on all tracks
      if (!track.has_collision()) {    // Track was not assigned, cannot compute NSigma (no event time) -> filling with empty table
        auto makeTableEmpty = [&](const Configurable<int>& flag, auto& table) {
          if (flag.value != 1) {
            return;
          }
          table(-999.f, -999.f);
        };

        makeTableEmpty(pidEl, tablePIDEl);
        makeTableEmpty(pidMu, tablePIDMu);
        makeTableEmpty(pidPi, tablePIDPi);
        makeTableEmpty(pidKa, tablePIDKa);
        makeTableEmpty(pidPr, tablePIDPr);
        makeTableEmpty(pidDe, tablePIDDe);
        makeTableEmpty(pidTr, tablePIDTr);
        makeTableEmpty(pidHe, tablePIDHe);
        makeTableEmpty(pidAl, tablePIDAl);

        continue;
      }

      if (enableTimeDependentResponse && (track.collisionId() != lastCollisionId)) { // Time dependent calib is enabled and this is a new collision
        lastCollisionId = track.collisionId();                                       // Cache last collision ID
        timestamp.value = track.collision().bc_as<aod::BCsWithTimestamps>().timestamp();
        LOG(debug) << "Updating parametrization from path '" << parametrizationPath << "' and timestamp " << timestamp.value;
        if (useParamCollection) {
          if (!ccdb->getForTimeStamp<o2::tof::ParameterCollection>(parametrizationPath, timestamp.value)->retrieveParameters(mRespParamsV2, passName.value)) {
            if (fatalOnPassNotAvailable) {
              LOGF(fatal, "Pass '%s' not available in the retrieved CCDB object", passName.value.data());
            } else {
              LOGF(warning, "Pass '%s' not available in the retrieved CCDB object", passName.value.data());
            }
          }
        } else {
          mRespParams.SetParameters(ccdb->getForTimeStamp<o2::pid::tof::TOFResoParams>(parametrizationPath, timestamp));
        }
      }

      // Check and fill enabled tables
      auto makeTable = [&track, this](const Configurable<int>& flag, auto& table, const auto& responsePID) {
        if (flag.value != 1) {
          return;
        }
        if (useParamCollection) {
          table(responsePID.GetExpectedSigma(mRespParamsV2, track),
                responsePID.GetSeparation(mRespParamsV2, track));
        } else {
          table(responsePID.GetExpectedSigma(mRespParams, track),
                responsePID.GetSeparation(mRespParams, track));
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
  PROCESS_SWITCH(tofPidFull, processWoSlice, "Process without track slices", false);

  void processWoSliceDev(Trks const& tracks, aod::Collisions const&, aod::BCsWithTimestamps const&)
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

#define doReserveTable(Particle)               \
  if (pid##Particle.value == 1) {              \
    tablePID##Particle.reserve(tracks.size()); \
  }

    doReserveTable(El);
    doReserveTable(Mu);
    doReserveTable(Pi);
    doReserveTable(Ka);
    doReserveTable(Pr);
    doReserveTable(De);
    doReserveTable(Tr);
    doReserveTable(He);
    doReserveTable(Al);

#undef doReserveTable

    int lastCollisionId = -1;          // Last collision ID analysed
    for (auto const& track : tracks) { // Loop on all tracks
      if (!track.has_collision()) {    // Track was not assigned, cannot compute NSigma (no event time) -> filling with empty table

#define doFillTableEmpty(Particle) \
  if (pid##Particle.value == 1) {  \
    tablePID##Particle(-999.f,     \
                       -999.f);    \
  }

        doFillTableEmpty(El);
        doFillTableEmpty(Mu);
        doFillTableEmpty(Pi);
        doFillTableEmpty(Ka);
        doFillTableEmpty(Pr);
        doFillTableEmpty(De);
        doFillTableEmpty(Tr);
        doFillTableEmpty(He);
        doFillTableEmpty(Al);

#undef doFillTableEmpty

        continue;
      }

      if (enableTimeDependentResponse && (track.collisionId() != lastCollisionId)) { // Time dependent calib is enabled and this is a new collision
        lastCollisionId = track.collisionId();                                       // Cache last collision ID
        timestamp.value = track.collision().bc_as<aod::BCsWithTimestamps>().timestamp();
        LOG(debug) << "Updating parametrization from path '" << parametrizationPath << "' and timestamp " << timestamp.value;
        mRespParams.SetParameters(ccdb->getForTimeStamp<o2::pid::tof::TOFResoParams>(parametrizationPath, timestamp));
      }

// Check and fill enabled tables
#define doFillTable(Particle)                                                   \
  if (pid##Particle.value == 1) {                                               \
    tablePID##Particle(response##Particle.GetExpectedSigma(mRespParams, track), \
                       response##Particle.GetSeparation(mRespParams, track));   \
  }

      doFillTable(El);
      doFillTable(Mu);
      doFillTable(Pi);
      doFillTable(Ka);
      doFillTable(Pr);
      doFillTable(De);
      doFillTable(Tr);
      doFillTable(He);
      doFillTable(Al);

#undef doFillTable
    }
  }
  PROCESS_SWITCH(tofPidFull, processWoSliceDev, "Process without track slices dev", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<tofPidFull>(cfgc)};
}
