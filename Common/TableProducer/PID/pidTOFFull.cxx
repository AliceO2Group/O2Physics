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
  o2::pid::tof::TOFResoParamsV2 mRespParamsV2;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> paramFileName{"paramFileName", "", "Path to the parametrization object. If empty the parametrization is not taken from file"};
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> parametrizationPath{"parametrizationPath", "TOF/Calib/Params", "Path of the TOF parametrization on the CCDB or in the file, if the paramFileName is not empty"};
  Configurable<std::string> timeShiftCCDBPath{"timeShiftCCDBPath", "", "Path of the TOF time shift vs eta. If empty none is taken"}; // temporary for pass3
  Configurable<std::string> passName{"passName", "", "Name of the pass inside of the CCDB parameter collection. If empty, the automatically deceted from metadata (to be implemented!!!)"};
  Configurable<int64_t> timestamp{"ccdb-timestamp", -1, "timestamp of the object"};

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
    if (doprocessWSlice == true && doprocessWoSlice == true) {
      LOGF(fatal, "Cannot enable processWoSlice and processWSlice at the same time. Please choose one.");
    }
    if (doprocessWSlice == false && doprocessWoSlice == false) {
      LOGF(fatal, "Cannot run without any of processWoSlice and processWSlice enabled. Please choose one.");
    }

    // Checking the tables are requested in the workflow and enabling them
    for (int i = 0; i < nSpecies; i++) {
      int f = enableParticle->get(particleNames[i].c_str(), "Enable");
      enableFlagIfTableRequired(initContext, "pidTOFFull" + particleNames[i], f);
      if (f == 1) {
        mEnabledParticles.push_back(i);
      }
    }
    // Printing enabled tables
    LOG(info) << "++ Enabled tables:";
    for (const int& i : mEnabledParticles) {
      LOG(info) << "++  pidTOFFull" << particleNames[i] << " is enabled";
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
        tablePIDEl(-999.f, -999.f);
        break;
      case 1:
        tablePIDMu(-999.f, -999.f);
        break;
      case 2:
        tablePIDPi(-999.f, -999.f);
        break;
      case 3:
        tablePIDKa(-999.f, -999.f);
        break;
      case 4:
        tablePIDPr(-999.f, -999.f);
        break;
      case 5:
        tablePIDDe(-999.f, -999.f);
        break;
      case 6:
        tablePIDTr(-999.f, -999.f);
        break;
      case 7:
        tablePIDHe(-999.f, -999.f);
        break;
      case 8:
        tablePIDAl(-999.f, -999.f);
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
    float resolution = 1.f;            // Last resolution assigned
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
              resolution = responseEl.GetExpectedSigma(mRespParamsV2, trkInColl);
              tablePIDEl(resolution,
                         responseEl.GetSeparation(mRespParamsV2, trkInColl, resolution));
              break;
            case 1:
              resolution = responseMu.GetExpectedSigma(mRespParamsV2, trkInColl);
              tablePIDMu(resolution,
                         responseMu.GetSeparation(mRespParamsV2, trkInColl, resolution));
              break;
            case 2:
              resolution = responsePi.GetExpectedSigma(mRespParamsV2, trkInColl);
              tablePIDPi(resolution,
                         responsePi.GetSeparation(mRespParamsV2, trkInColl, resolution));
              break;
            case 3:
              resolution = responseKa.GetExpectedSigma(mRespParamsV2, trkInColl);
              tablePIDKa(resolution,
                         responseKa.GetSeparation(mRespParamsV2, trkInColl, resolution));
              break;
            case 4:
              resolution = responsePr.GetExpectedSigma(mRespParamsV2, trkInColl);
              tablePIDPr(resolution,
                         responsePr.GetSeparation(mRespParamsV2, trkInColl, resolution));
              break;
            case 5:
              resolution = responseDe.GetExpectedSigma(mRespParamsV2, trkInColl);
              tablePIDDe(resolution,
                         responseDe.GetSeparation(mRespParamsV2, trkInColl, resolution));
              break;
            case 6:
              resolution = responseTr.GetExpectedSigma(mRespParamsV2, trkInColl);
              tablePIDTr(resolution,
                         responseTr.GetSeparation(mRespParamsV2, trkInColl, resolution));
              break;
            case 7:
              resolution = responseHe.GetExpectedSigma(mRespParamsV2, trkInColl);
              tablePIDHe(resolution,
                         responseHe.GetSeparation(mRespParamsV2, trkInColl, resolution));
              break;
            case 8:
              resolution = responseAl.GetExpectedSigma(mRespParamsV2, trkInColl);
              tablePIDAl(resolution,
                         responseAl.GetSeparation(mRespParamsV2, trkInColl, resolution));
              break;
            default:
              LOG(fatal) << "Wrong particle ID in processWSlice()";
              break;
          }
        }
      }
    }
  }
  PROCESS_SWITCH(tofPidFull, processWSlice, "Process with track slices", true);

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
    float resolution = 1.f;            // Last resolution assigned
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
            resolution = responseEl.GetExpectedSigma(mRespParamsV2, track);
            tablePIDEl(resolution,
                       responseEl.GetSeparation(mRespParamsV2, track, resolution));
            break;
          case 1:
            resolution = responseMu.GetExpectedSigma(mRespParamsV2, track);
            tablePIDMu(resolution,
                       responseMu.GetSeparation(mRespParamsV2, track, resolution));
            break;
          case 2:
            resolution = responsePi.GetExpectedSigma(mRespParamsV2, track);
            tablePIDPi(resolution,
                       responsePi.GetSeparation(mRespParamsV2, track));
            break;
          case 3:
            resolution = responseKa.GetExpectedSigma(mRespParamsV2, track);
            tablePIDKa(resolution,
                       responseKa.GetSeparation(mRespParamsV2, track, resolution));
            break;
          case 4:
            resolution = responsePr.GetExpectedSigma(mRespParamsV2, track);
            tablePIDPr(resolution,
                       responsePr.GetSeparation(mRespParamsV2, track, resolution));
            break;
          case 5:
            resolution = responseDe.GetExpectedSigma(mRespParamsV2, track);
            tablePIDDe(resolution,
                       responseDe.GetSeparation(mRespParamsV2, track, resolution));
            break;
          case 6:
            resolution = responseTr.GetExpectedSigma(mRespParamsV2, track);
            tablePIDTr(resolution,
                       responseTr.GetSeparation(mRespParamsV2, track, resolution));
            break;
          case 7:
            resolution = responseHe.GetExpectedSigma(mRespParamsV2, track);
            tablePIDHe(resolution,
                       responseHe.GetSeparation(mRespParamsV2, track, resolution));
            break;
          case 8:
            resolution = responseAl.GetExpectedSigma(mRespParamsV2, track);
            tablePIDAl(resolution,
                       responseAl.GetSeparation(mRespParamsV2, track, resolution));
            break;
          default:
            LOG(fatal) << "Wrong particle ID in processWoSlice()";
            break;
        }
      }
    }
  }
  PROCESS_SWITCH(tofPidFull, processWoSlice, "Process without track slices and on TrackIU (faster but only Run3)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<tofPidFull>(cfgc)}; }
