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
/// \file   mcPidTof.cxx
/// \author Fabrizio Grosa fabrizio.grosa@cern.ch
/// \brief  Task to produce PID tables for TOF split for pi, K, p, copied from https://github.com/AliceO2Group/O2Physics/blob/master/Common/TableProducer/PID/pidTOFFull.cxx
///         In addition, it applies postcalibrations for MC.
///

// O2 includes
#include <CCDB/BasicCCDBManager.h>
#include "TOFBase/EventTimeMaker.h"
#include "Framework/AnalysisTask.h"
#include "ReconstructionDataFormats/Track.h"

// O2Physics includes
#include "TableHelper.h"
#include "Common/TableProducer/PID/pidTOFBase.h"

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
struct mcPidTof {
  // Tables to produce
  Produces<o2::aod::pidTOFFullPi> tablePIDPi;
  Produces<o2::aod::pidTOFFullKa> tablePIDKa;
  Produces<o2::aod::pidTOFFullPr> tablePIDPr;
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
  // postcalibrations to overcome MC FT0 timing issue
  std::map<int, TGraph*> gMcPostCalibMean{};
  std::map<int, TGraph*> gMcPostCalibSigma{};
  int currentRun{0};
  struct : ConfigurableGroup {
    std::string prefix = "mcRecalib";
    Configurable<bool> enable{"enable", false, "enable MC recalibration for Pi/Ka/Pr"};
    Configurable<std::string> ccdbPath{"ccdbPath", "Users/f/fgrosa/RecalibMcPidTOF/", "path for MC recalibration objects in CCDB"};
    Configurable<std::string> passName{"passName", "apass6", "reco pass of MC anchoring"};
  } mcRecalib;

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

    // Checking the tables are requested in the workflow and enabling them (only pi, K, p)
    std::array<int, 3> supportedSpecies = {2, 3, 4};
    for (auto iSpecie{0u}; iSpecie < supportedSpecies.size(); ++iSpecie) {
      int flag = -1;
      enableFlagIfTableRequired(initContext, "mcPidTof" + particleNames[supportedSpecies[iSpecie]], flag);
    }
    // Printing enabled tables
    LOG(info) << "++ Enabled tables:";
    for (const int& pidId : mEnabledParticles) {
      LOG(info) << "++  mcPidTof" << particleNames[pidId] << " is enabled";
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
      case 2:
        tablePIDPi.reserve(size);
        break;
      case 3:
        tablePIDKa.reserve(size);
        break;
      case 4:
        tablePIDPr.reserve(size);
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
      case 2:
        tablePIDPi(-999.f, -999.f);
        break;
      case 3:
        tablePIDKa(-999.f, -999.f);
        break;
      case 4:
        tablePIDPr(-999.f, -999.f);
        break;
      default:
        LOG(fatal) << "Wrong particle ID in makeTableEmpty()";
        break;
    }
  }

  /// Retrieve MC postcalibration objects from CCDB
  /// \param timestamp timestamp
  void retrieveMcPostCalibFromCcdb(int timestamp)
  {
    std::map<std::string, std::string> metadata;
    metadata["RecoPassName"] = mcRecalib.passName;
    auto calibList = ccdb->getSpecific<TList>(mcRecalib.ccdbPath, timestamp, metadata);
    for (auto const& pidId : mEnabledParticles) { // Loop on enabled particle hypotheses
      gMcPostCalibMean[pidId] = reinterpret_cast<TGraph*>(calibList->FindObject(Form("Mean%s", particleNames[pidId].data())));
      gMcPostCalibSigma[pidId] = reinterpret_cast<TGraph*>(calibList->FindObject(Form("Sigma%s", particleNames[pidId].data())));
    }
  }

  /// Apply MC postcalibrations
  /// \param pidId particle id
  /// \param pt track pT
  template <typename T>
  T applyMcRecalib(int pidId, T trackPt, T nSigma)
  {
    float shift{0.f}, scaleWidth{0.f};
    int nPoints = gMcPostCalibMean[pidId]->GetN();
    double ptMin = gMcPostCalibMean[pidId]->GetX()[0];
    double ptMax = gMcPostCalibMean[pidId]->GetX()[nPoints - 1];
    if (trackPt < ptMin) {
      shift = gMcPostCalibMean[pidId]->Eval(ptMin);
      scaleWidth = gMcPostCalibSigma[pidId]->Eval(ptMin);
    } else if (trackPt > ptMax) {
      shift = gMcPostCalibMean[pidId]->Eval(ptMax);
      scaleWidth = gMcPostCalibSigma[pidId]->Eval(ptMax);
    } else {
      shift = gMcPostCalibMean[pidId]->Eval(trackPt);
      scaleWidth = gMcPostCalibSigma[pidId]->Eval(trackPt);
    }

    T nSigmaCorr = (nSigma - shift) / scaleWidth;
    return nSigmaCorr;
  }

  using Trks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TOFSignal, aod::TOFEvTime, aod::pidEvTimeFlags>;
  // Define slice per collision
  Preslice<Trks> perCollision = aod::track::collisionId;
  template <o2::track::PID::ID pid>
  using ResponseImplementation = o2::pid::tof::ExpTimes<Trks::iterator, pid>;
  void processWSlice(Trks const& tracks, aod::Collisions const&, aod::BCsWithTimestamps const&)
  {
    constexpr auto responsePi = ResponseImplementation<PID::Pion>();
    constexpr auto responseKa = ResponseImplementation<PID::Kaon>();
    constexpr auto responsePr = ResponseImplementation<PID::Proton>();

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

      // in case of MC recalibrations, check if the objects from CCDB has to be updated
      if (mcRecalib.enable) {
        int runNumber = track.collision().bc_as<aod::BCsWithTimestamps>().runNumber();
        if (runNumber != currentRun) {
          // update postcalibration files
          retrieveMcPostCalibFromCcdb(timestamp);
        }
        currentRun = runNumber;
      }

      const auto& tracksInCollision = tracks.sliceBy(perCollision, lastCollisionId);
      for (auto const& trkInColl : tracksInCollision) { // Loop on tracks
        for (auto const& pidId : mEnabledParticles) {   // Loop on enabled particle hypotheses
          float nSigma{-999.f};
          switch (pidId) {
            case 2:
              resolution = responsePi.GetExpectedSigma(mRespParamsV2, trkInColl);
              nSigma = responsePi.GetSeparation(mRespParamsV2, trkInColl, resolution);
              if (mcRecalib.enable) {
                nSigma = applyMcRecalib(pidId, trkInColl.pt(), nSigma);
              }
              tablePIDPi(resolution, nSigma);
              break;
            case 3:
              resolution = responseKa.GetExpectedSigma(mRespParamsV2, trkInColl);
              nSigma = responseKa.GetSeparation(mRespParamsV2, trkInColl, resolution);
              if (mcRecalib.enable) {
                nSigma = applyMcRecalib(pidId, trkInColl.pt(), nSigma);
              }
              tablePIDKa(resolution, nSigma);
              break;
            case 4:
              resolution = responsePr.GetExpectedSigma(mRespParamsV2, trkInColl);
              nSigma = responsePr.GetSeparation(mRespParamsV2, trkInColl, resolution);
              if (mcRecalib.enable) {
                nSigma = applyMcRecalib(pidId, trkInColl.pt(), nSigma);
              }
              tablePIDPr(resolution, nSigma);
              break;
            default:
              LOG(fatal) << "Wrong particle ID in processWSlice()";
              break;
          }
        }
      }
    }
  }
  PROCESS_SWITCH(mcPidTof, processWSlice, "Process with track slices", true);

  using TrksIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TOFSignal, aod::TOFEvTime, aod::pidEvTimeFlags>;
  template <o2::track::PID::ID pid>
  using ResponseImplementationIU = o2::pid::tof::ExpTimes<TrksIU::iterator, pid>;
  void processWoSlice(TrksIU const& tracks, aod::Collisions const&, aod::BCsWithTimestamps const&)
  {
    constexpr auto responsePi = ResponseImplementationIU<PID::Pion>();
    constexpr auto responseKa = ResponseImplementationIU<PID::Kaon>();
    constexpr auto responsePr = ResponseImplementationIU<PID::Proton>();

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

      // in case of MC recalibrations, check if the objects from CCDB has to be updated
      if (mcRecalib.enable) {
        int runNumber = track.collision().bc_as<aod::BCsWithTimestamps>().runNumber();
        if (runNumber != currentRun) {
          // update postcalibration files
          retrieveMcPostCalibFromCcdb(timestamp);
        }
        currentRun = runNumber;
      }

      float nSigma{-999.f};
      for (auto const& pidId : mEnabledParticles) { // Loop on enabled particle hypotheses
        switch (pidId) {
          case 2:
            resolution = responsePi.GetExpectedSigma(mRespParamsV2, track);
            nSigma = responsePi.GetSeparation(mRespParamsV2, track, resolution);
            if (mcRecalib.enable) {
              nSigma = applyMcRecalib(pidId, track.pt(), nSigma);
            }
            tablePIDPi(resolution, nSigma);
            break;
          case 3:
            resolution = responseKa.GetExpectedSigma(mRespParamsV2, track);
            nSigma = responseKa.GetSeparation(mRespParamsV2, track, resolution);
            if (mcRecalib.enable) {
              nSigma = applyMcRecalib(pidId, track.pt(), nSigma);
            }
            tablePIDKa(resolution, nSigma);
            break;
          case 4:
            resolution = responsePr.GetExpectedSigma(mRespParamsV2, track);
            nSigma = responsePr.GetSeparation(mRespParamsV2, track, resolution);
            if (mcRecalib.enable) {
              nSigma = applyMcRecalib(pidId, track.pt(), nSigma);
            }
            tablePIDPr(resolution, nSigma);
            break;
          default:
            LOG(fatal) << "Wrong particle ID in processWoSlice()";
            break;
        }
      }
    }
  }
  PROCESS_SWITCH(mcPidTof, processWoSlice, "Process without track slices and on TrackIU (faster but only Run3)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<mcPidTof>(cfgc)}; }
