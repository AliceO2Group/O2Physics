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
///         Only the tables for the mass hypotheses requested are filled, the others are sent empty.
///

#include "pidTOFBase.h"

#include "Common/Core/TableHelper.h"
#include "Common/DataModel/PIDResponseTOF.h"

#include <CCDB/BasicCCDBManager.h>
#include <DataFormatsTOF/ParameterContainers.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/Variant.h>
#include <PID/PIDTOF.h>
#include <ReconstructionDataFormats/PID.h>

#include <TGraph.h>
#include <TString.h>

#include <chrono>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

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

#include <Framework/runDataProcessing.h>

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
  // Detector response parameters
  o2::pid::tof::TOFResoParamsV2 mRespParamsV2;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<bool> inheritFromBaseTask{"inheritFromBaseTask", true, "Flag to iherit all common configurables from the TOF base task"};
  // CCDB configuration (inherited from TOF base task)
  Configurable<std::string> url{"ccdb-url", "", "url of the ccdb repository"};
  Configurable<int64_t> timestamp{"ccdb-timestamp", -1, "timestamp of the object"};
  // TOF Calib configuration (inherited from TOF base task)
  Configurable<std::string> paramFileName{"paramFileName", "", "Path to the parametrization object. If empty the parametrization is not taken from file"};
  Configurable<std::string> parametrizationPath{"parametrizationPath", "", "Path of the TOF parametrization on the CCDB or in the file, if the paramFileName is not empty"};
  Configurable<std::string> passName{"passName", "", "Name of the pass inside of the CCDB parameter collection. If empty, the automatically deceted from metadata (to be implemented!!!)"};
  Configurable<std::string> timeShiftCCDBPath{"timeShiftCCDBPath", "", "Path of the TOF time shift vs eta. If empty none is taken"};
  Configurable<bool> loadResponseFromCCDB{"loadResponseFromCCDB", false, "Flag to load the response from the CCDB"};
  Configurable<bool> enableTimeDependentResponse{"enableTimeDependentResponse", false, "Flag to use the collision timestamp to fetch the PID Response"};
  Configurable<bool> fatalOnPassNotAvailable{"fatalOnPassNotAvailable", true, "Flag to throw a fatal if the pass is not available in the retrieved CCDB object"};
  // Configuration flags to include and exclude particle hypotheses
  Configurable<LabeledArray<int>> enableParticle{"enableParticle",
                                                 {defaultParameters[0], nSpecies, nParameters, particleNames, parameterNames},
                                                 "Produce PID information for the various mass hypotheses. Values different than -1 override the automatic setup: the corresponding table can be set off (0) or on (1)"};

  // Running variables
  std::vector<int> mEnabledParticles; // Vector of enabled PID hypotheses to loop on when making tables
  int mLastCollisionId = -1;          // Last collision ID analysed
  void init(o2::framework::InitContext& initContext)
  {
    if (inheritFromBaseTask.value) { // Inheriting from base task
      if (!getTaskOptionValue(initContext, "tof-signal", "ccdb-url", url.value, true)) {
        LOG(fatal) << "Could not get ccdb-url from tof-signal task";
      }
      if (!getTaskOptionValue(initContext, "tof-signal", "ccdb-timestamp", timestamp.value, true)) {
        LOG(fatal) << "Could not get ccdb-timestamp from tof-signal task";
      }
      if (!getTaskOptionValue(initContext, "tof-event-time", "paramFileName", paramFileName.value, true)) {
        LOG(fatal) << "Could not get paramFileName from tof-event-time task";
      }
      if (!getTaskOptionValue(initContext, "tof-event-time", "parametrizationPath", parametrizationPath.value, true)) {
        LOG(fatal) << "Could not get parametrizationPath from tof-event-time task";
      }
      if (!getTaskOptionValue(initContext, "tof-event-time", "passName", passName.value, true)) {
        LOG(fatal) << "Could not get passName from tof-event-time task";
      }
      if (!getTaskOptionValue(initContext, "tof-signal", "timeShiftCCDBPath", timeShiftCCDBPath.value, true)) {
        LOG(fatal) << "Could not get timeShiftCCDBPath from tof-signal task";
      }
      if (!getTaskOptionValue(initContext, "tof-event-time", "loadResponseFromCCDB", loadResponseFromCCDB.value, true)) {
        LOG(fatal) << "Could not get loadResponseFromCCDB from tof-event-time task";
      }
      if (!getTaskOptionValue(initContext, "tof-event-time", "enableTimeDependentResponse", enableTimeDependentResponse.value, true)) {
        LOG(fatal) << "Could not get enableTimeDependentResponse from tof-event-time task";
      }
      if (!getTaskOptionValue(initContext, "tof-event-time", "fatalOnPassNotAvailable", fatalOnPassNotAvailable.value, true)) {
        LOG(fatal) << "Could not get fatalOnPassNotAvailable from tof-event-time task";
      }
    }
    if (doprocessWSlice == true && doprocessWoSlice == true) {
      LOGF(fatal, "Cannot enable processWoSlice and processWSlice at the same time. Please choose one.");
    }
    if (doprocessWSlice == false && doprocessWoSlice == false) {
      LOGF(fatal, "Cannot run without any of processWoSlice and processWSlice enabled. Please choose one.");
    }

    // Checking the tables are requested in the workflow and enabling them
    for (int i = 0; i < nSpecies; i++) {
      int f = enableParticle->get(particleNames[i].c_str(), "Enable");
      enableFlagIfTableRequired(initContext, "pidTOF" + particleNames[i], f);
      if (f == 1) {
        mEnabledParticles.push_back(i);
      }
    }
    // Printing enabled tables
    LOG(info) << "++ Enabled tables:";
    for (const int& i : mEnabledParticles) {
      LOG(info) << "++  pidTOF" << particleNames[i] << " is enabled";
    }

    // Getting the parametrization parameters
    ccdb->setURL(url.value);
    ccdb->setTimestamp(timestamp.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    // Not later than now objects
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    //

    // TODO: implement the automatic pass name detection from metadata
    if (passName.value == "") {
      passName.value = "unanchored"; // temporary default
      LOG(warning) << "Passed autodetect mode for pass, not implemented yet, waiting for metadata. Taking '" << passName.value << "'";
    }
    LOG(info) << "Using parameter collection, starting from pass '" << passName.value << "'";

    const std::string fname = paramFileName.value;
    if (!fname.empty()) { // Loading the parametrization from file
      LOG(info) << "Loading exp. sigma parametrization from file " << fname << ", using param: " << parametrizationPath.value;
      if (1) {
        o2::tof::ParameterCollection paramCollection;
        paramCollection.loadParamFromFile(fname, parametrizationPath.value);
        LOG(info) << "+++ Loaded parameter collection from file +++";
        if (!paramCollection.retrieveParameters(mRespParamsV2, passName.value)) {
          if (fatalOnPassNotAvailable) {
            LOGF(fatal, "Pass '%s' not available in the retrieved CCDB object", passName.value.data());
          } else {
            LOGF(warning, "Pass '%s' not available in the retrieved CCDB object", passName.value.data());
          }
        } else {
          mRespParamsV2.setShiftParameters(paramCollection.getPars(passName.value));
          mRespParamsV2.printShiftParameters();
        }
      } else {
        mRespParamsV2.loadParamFromFile(fname.data(), parametrizationPath.value);
      }
    } else if (loadResponseFromCCDB) { // Loading it from CCDB
      LOG(info) << "Loading exp. sigma parametrization from CCDB, using path: " << parametrizationPath.value << " for timestamp " << timestamp.value;
      o2::tof::ParameterCollection* paramCollection = ccdb->getForTimeStamp<o2::tof::ParameterCollection>(parametrizationPath.value, timestamp.value);
      paramCollection->print();
      if (!paramCollection->retrieveParameters(mRespParamsV2, passName.value)) { // Attempt at loading the parameters with the pass defined
        if (fatalOnPassNotAvailable) {
          LOGF(fatal, "Pass '%s' not available in the retrieved CCDB object", passName.value.data());
        } else {
          LOGF(warning, "Pass '%s' not available in the retrieved CCDB object", passName.value.data());
        }
      } else { // Pass is available, load non standard parameters
        mRespParamsV2.setShiftParameters(paramCollection->getPars(passName.value));
        mRespParamsV2.printShiftParameters();
      }
    }
    mRespParamsV2.print();
    if (timeShiftCCDBPath.value != "") {
      if (timeShiftCCDBPath.value.find(".root") != std::string::npos) {
        mRespParamsV2.setTimeShiftParameters(timeShiftCCDBPath.value, "gmean_Pos", true);
        mRespParamsV2.setTimeShiftParameters(timeShiftCCDBPath.value, "gmean_Neg", false);
      } else {
        mRespParamsV2.setTimeShiftParameters(ccdb->getForTimeStamp<TGraph>(Form("%s/pos", timeShiftCCDBPath.value.c_str()), timestamp.value), true);
        mRespParamsV2.setTimeShiftParameters(ccdb->getForTimeStamp<TGraph>(Form("%s/neg", timeShiftCCDBPath.value.c_str()), timestamp.value), false);
      }
    }
  }

  // Reserves an empty table for the given particle ID with size of the given track table
  void reserveTable(const int id, const int64_t& size)
  {
    switch (id) {
      case 0:
        tablePIDEl.reserve(size);
        break;
      case 1:
        tablePIDMu.reserve(size);
        break;
      case 2:
        tablePIDPi.reserve(size);
        break;
      case 3:
        tablePIDKa.reserve(size);
        break;
      case 4:
        tablePIDPr.reserve(size);
        break;
      case 5:
        tablePIDDe.reserve(size);
        break;
      case 6:
        tablePIDTr.reserve(size);
        break;
      case 7:
        tablePIDHe.reserve(size);
        break;
      case 8:
        tablePIDAl.reserve(size);
        break;
      default:
        LOG(fatal) << "Wrong particle ID in reserveTable()";
        break;
    }
  }

  // Makes the table empty for the given particle ID, filling it with dummy values
  void makeTableEmpty(const int id)
  {
    switch (id) {
      case 0:
        aod::pidtof_tiny::binning::packInTable(-999.f, tablePIDEl);
        break;
      case 1:
        aod::pidtof_tiny::binning::packInTable(-999.f, tablePIDMu);
        break;
      case 2:
        aod::pidtof_tiny::binning::packInTable(-999.f, tablePIDPi);
        break;
      case 3:
        aod::pidtof_tiny::binning::packInTable(-999.f, tablePIDKa);
        break;
      case 4:
        aod::pidtof_tiny::binning::packInTable(-999.f, tablePIDPr);
        break;
      case 5:
        aod::pidtof_tiny::binning::packInTable(-999.f, tablePIDDe);
        break;
      case 6:
        aod::pidtof_tiny::binning::packInTable(-999.f, tablePIDTr);
        break;
      case 7:
        aod::pidtof_tiny::binning::packInTable(-999.f, tablePIDHe);
        break;
      case 8:
        aod::pidtof_tiny::binning::packInTable(-999.f, tablePIDAl);
        break;
      default:
        LOG(fatal) << "Wrong particle ID in makeTableEmpty()";
        break;
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

    for (auto const& pidId : mEnabledParticles) {
      reserveTable(pidId, tracks.size());
    }

    int lastCollisionId = -1;          // Last collision ID analysed
    for (auto const& track : tracks) { // Loop on all tracks
      if (!track.has_collision()) {    // Track was not assigned, cannot compute NSigma (no event time) -> filling with empty table
        for (auto const& pidId : mEnabledParticles) {
          makeTableEmpty(pidId);
        }
        continue;
      }

      if (track.collisionId() == lastCollisionId) { // Tracks from last collision already processed
        continue;
      }

      // Fill new table for the tracks in a collision
      lastCollisionId = track.collisionId(); // Cache last collision ID
      timestamp.value = track.collision().bc_as<aod::BCsWithTimestamps>().timestamp();
      if (enableTimeDependentResponse) {
        LOG(debug) << "Updating parametrization from path '" << parametrizationPath.value << "' and timestamp " << timestamp.value;
        if (!ccdb->getForTimeStamp<o2::tof::ParameterCollection>(parametrizationPath.value, timestamp.value)->retrieveParameters(mRespParamsV2, passName.value)) {
          if (fatalOnPassNotAvailable) {
            LOGF(fatal, "Pass '%s' not available in the retrieved CCDB object", passName.value.data());
          } else {
            LOGF(warning, "Pass '%s' not available in the retrieved CCDB object", passName.value.data());
          }
        }
      }

      const auto& tracksInCollision = tracks.sliceBy(perCollision, lastCollisionId);
      for (auto const& trkInColl : tracksInCollision) { // Loop on tracks
        for (auto const& pidId : mEnabledParticles) {   // Loop on enabled particle hypotheses
          switch (pidId) {
            case 0:
              aod::pidtof_tiny::binning::packInTable(responseEl.GetSeparation(mRespParamsV2, trkInColl),
                                                     tablePIDEl);
              break;
            case 1:
              aod::pidtof_tiny::binning::packInTable(responseMu.GetSeparation(mRespParamsV2, trkInColl),
                                                     tablePIDMu);
              break;
            case 2:
              aod::pidtof_tiny::binning::packInTable(responsePi.GetSeparation(mRespParamsV2, trkInColl),
                                                     tablePIDPi);
              break;
            case 3:
              aod::pidtof_tiny::binning::packInTable(responseKa.GetSeparation(mRespParamsV2, trkInColl),
                                                     tablePIDKa);
              break;
            case 4:
              aod::pidtof_tiny::binning::packInTable(responsePr.GetSeparation(mRespParamsV2, trkInColl),
                                                     tablePIDPr);
              break;
            case 5:
              aod::pidtof_tiny::binning::packInTable(responseDe.GetSeparation(mRespParamsV2, trkInColl),
                                                     tablePIDDe);
              break;
            case 6:
              aod::pidtof_tiny::binning::packInTable(responseTr.GetSeparation(mRespParamsV2, trkInColl),
                                                     tablePIDTr);
              break;
            case 7:
              aod::pidtof_tiny::binning::packInTable(responseHe.GetSeparation(mRespParamsV2, trkInColl),
                                                     tablePIDHe);
              break;
            case 8:
              aod::pidtof_tiny::binning::packInTable(responseAl.GetSeparation(mRespParamsV2, trkInColl),
                                                     tablePIDAl);
              break;
            default:
              LOG(fatal) << "Wrong particle ID in processWSlice()";
              break;
          }
        }
      }
    }
  }
  PROCESS_SWITCH(tofPid, processWSlice, "Process with track slices", true);

  using TrksIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TOFSignal, aod::TOFEvTime, aod::pidEvTimeFlags>;
  template <o2::track::PID::ID pid>
  using ResponseImplementationIU = o2::pid::tof::ExpTimes<TrksIU::iterator, pid>;
  void processWoSlice(TrksIU const& tracks, aod::Collisions const&, aod::BCsWithTimestamps const&)
  {
    constexpr auto responseEl = ResponseImplementationIU<PID::Electron>();
    constexpr auto responseMu = ResponseImplementationIU<PID::Muon>();
    constexpr auto responsePi = ResponseImplementationIU<PID::Pion>();
    constexpr auto responseKa = ResponseImplementationIU<PID::Kaon>();
    constexpr auto responsePr = ResponseImplementationIU<PID::Proton>();
    constexpr auto responseDe = ResponseImplementationIU<PID::Deuteron>();
    constexpr auto responseTr = ResponseImplementationIU<PID::Triton>();
    constexpr auto responseHe = ResponseImplementationIU<PID::Helium3>();
    constexpr auto responseAl = ResponseImplementationIU<PID::Alpha>();

    for (auto const& pidId : mEnabledParticles) {
      reserveTable(pidId, tracks.size());
    }

    for (auto const& track : tracks) { // Loop on all tracks
      if (!track.has_collision()) {    // Track was not assigned, cannot compute NSigma (no event time) -> filling with empty table
        for (auto const& pidId : mEnabledParticles) {
          makeTableEmpty(pidId);
        }
        continue;
      }

      if (enableTimeDependentResponse && (track.collisionId() != mLastCollisionId)) { // Time dependent calib is enabled and this is a new collision
        mLastCollisionId = track.collisionId();                                       // Cache last collision ID
        timestamp.value = track.collision().bc_as<aod::BCsWithTimestamps>().timestamp();
        LOG(debug) << "Updating parametrization from path '" << parametrizationPath.value << "' and timestamp " << timestamp.value;
        if (!ccdb->getForTimeStamp<o2::tof::ParameterCollection>(parametrizationPath.value, timestamp.value)->retrieveParameters(mRespParamsV2, passName.value)) {
          if (fatalOnPassNotAvailable) {
            LOGF(fatal, "Pass '%s' not available in the retrieved CCDB object", passName.value.data());
          } else {
            LOGF(warning, "Pass '%s' not available in the retrieved CCDB object", passName.value.data());
          }
        }
      }

      for (auto const& pidId : mEnabledParticles) { // Loop on enabled particle hypotheses
        switch (pidId) {
          case 0:
            aod::pidtof_tiny::binning::packInTable(responseEl.GetSeparation(mRespParamsV2, track),
                                                   tablePIDEl);
            break;
          case 1:
            aod::pidtof_tiny::binning::packInTable(responseMu.GetSeparation(mRespParamsV2, track),
                                                   tablePIDMu);
            break;
          case 2:
            aod::pidtof_tiny::binning::packInTable(responsePi.GetSeparation(mRespParamsV2, track),
                                                   tablePIDPi);
            break;
          case 3:
            aod::pidtof_tiny::binning::packInTable(responseKa.GetSeparation(mRespParamsV2, track),
                                                   tablePIDKa);
            break;
          case 4:
            aod::pidtof_tiny::binning::packInTable(responsePr.GetSeparation(mRespParamsV2, track),
                                                   tablePIDPr);
            break;
          case 5:
            aod::pidtof_tiny::binning::packInTable(responseDe.GetSeparation(mRespParamsV2, track),
                                                   tablePIDDe);
            break;
          case 6:
            aod::pidtof_tiny::binning::packInTable(responseTr.GetSeparation(mRespParamsV2, track),
                                                   tablePIDTr);
            break;
          case 7:
            aod::pidtof_tiny::binning::packInTable(responseHe.GetSeparation(mRespParamsV2, track),
                                                   tablePIDHe);
            break;
          case 8:
            aod::pidtof_tiny::binning::packInTable(responseAl.GetSeparation(mRespParamsV2, track),
                                                   tablePIDAl);
            break;
          default:
            LOG(fatal) << "Wrong particle ID in processWoSlice()";
            break;
        }
      }
    }
  }
  PROCESS_SWITCH(tofPid, processWoSlice, "Process without track slices and on TrackIU (faster but only Run3)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<tofPid>(cfgc)}; }
