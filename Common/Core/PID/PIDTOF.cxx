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
/// \file   PIDTOF.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \since  02/07/2020
/// \brief  Implementation of the TOF detector response for PID
///

#include "PIDTOF.h"

// O2Physics includes
#include "Common/Core/CollisionTypeHelper.h"

// O2 includes
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "Framework/Logger.h"

#include <map>
#include <string>

namespace o2::pid::tof
{

void TOFResoParamsV3::setResolutionParametrizationRun2(std::unordered_map<std::string, float> const& pars)
{
  std::array<std::string, 13> paramNames{"TrkRes.Pi.P0", "TrkRes.Pi.P1", "TrkRes.Pi.P2", "TrkRes.Pi.P3", "time_resolution",
                                         "TrkRes.Ka.P0", "TrkRes.Ka.P1", "TrkRes.Ka.P2", "TrkRes.Ka.P3",
                                         "TrkRes.Pr.P0", "TrkRes.Pr.P1", "TrkRes.Pr.P2", "TrkRes.Pr.P3"};
  // Now we override the parametrization to use the Run 2 one
  for (int i = 0; i < 9; i++) {
    if (mResolution[i]) {
      delete mResolution[i];
    }
    mResolution[i] = new TF2(Form("tofResTrack.%s_Run2", particleNames[i]), "-10", 0., 20, -1, 1.); // With negative values the old one is used
  }
  // Print the map
  for (const auto& [key, value] : pars) {
    LOG(info) << "Key: " << key << " Value: " << value;
  }
  for (int i = 0; i < 13; i++) {
    setParameter(i, pars.at(paramNames[i]));
  }
}

// Time shift for post calibration to realign as a function of eta
void TOFResoParamsV3::setTimeShiftParameters(std::unordered_map<std::string, float> const& pars, const bool positive)
{
  std::string baseOpt = positive ? "TimeShift.Pos." : "TimeShift.Neg.";

  if (pars.count(baseOpt + "GetN") == 0) { // If the map does not contain the number of eta bins, we assume that no correction has to be applied
    return;
  }
  const int nPoints = static_cast<int>(pars.at(baseOpt + "GetN"));
  if (nPoints <= 0) {
    LOG(fatal) << "TOFResoParamsV3 shift: time must be positive";
  }
  TGraph graph;
  for (int i = 0; i < nPoints; ++i) {
    graph.AddPoint(pars.at(Form("TimeShift.eta%i", i)), pars.at(Form("TimeShift.cor%i", i)));
  }
  setTimeShiftParameters(&graph, positive);
}
void TOFResoParamsV3::setTimeShiftParameters(std::string const& filename, std::string const& objname, const bool positive)
{
  TFile f(filename.c_str(), "READ");
  if (f.IsOpen()) {
    if (positive) {
      f.GetObject(objname.c_str(), gPosEtaTimeCorr);
    } else {
      f.GetObject(objname.c_str(), gNegEtaTimeCorr);
    }
    f.Close();
  }
  LOG(info) << "Set the Time Shift parameters from file " << filename << " and object " << objname << " for " << (positive ? "positive" : "negative");
}
void TOFResoParamsV3::setTimeShiftParameters(TGraph* g, const bool positive)
{
  if (!g) {
    LOG(info) << "No Time Shift parameter is passed for " << (positive ? "positive" : "negative");
    return;
  }
  if (positive) {
    gPosEtaTimeCorr = g;
  } else {
    gNegEtaTimeCorr = g;
  }
  LOG(info) << "Set the Time Shift parameters from object " << g->GetName() << " " << g->GetTitle() << " for " << (positive ? "positive" : "negative");
}
float TOFResoParamsV3::getTimeShift(float eta, int16_t sign) const
{
  if (sign > 0) {
    if (!gPosEtaTimeCorr) {
      return 0.f;
    }
    return gPosEtaTimeCorr->Eval(eta);
  }
  if (!gNegEtaTimeCorr) {
    return 0.f;
  }
  return gNegEtaTimeCorr->Eval(eta);
}

} // namespace o2::pid::tof

#include "Framework/CommonServices.h"
#include "Framework/Plugins.h"
#include "Framework/ServiceHandle.h"
#include "Framework/ServiceSpec.h"

#include "TDatabasePDG.h"

using namespace o2::framework;

o2::pid::tof::TOFResoParamsV3 o2::pid::tof::TOFResponseImpl::parameters;
o2::common::core::MetadataHelper o2::pid::tof::TOFResponseImpl::metadataInfo;
bool o2::pid::tof::TOFResponseImpl::mIsInit = false;
int o2::pid::tof::TOFResponseImpl::mLastRunNumber = -1;

void o2::pid::tof::TOFResponseImpl::inheritFromBaseTask(o2::framework::InitContext& initContext, const std::string task)
{
  if (mIsInit) {
    LOG(fatal) << "TOFResponseImpl already initialized, cannot re-initialize";
  }
  getCfg(initContext, "ccdb-url", mUrl, task);
  getCfg(initContext, "ccdb-path-grplhcif", mPathGrpLhcIf, task);
  getCfg(initContext, "ccdb-timestamp", mTimestamp, task);
  getCfg(initContext, "timeShiftCCDBPathPos", mTimeShiftCCDBPathPos, task);
  getCfg(initContext, "timeShiftCCDBPathNeg", mTimeShiftCCDBPathNeg, task);
  getCfg(initContext, "timeShiftCCDBPathPosMC", mTimeShiftCCDBPathPosMC, task);
  getCfg(initContext, "timeShiftCCDBPathNegMC", mTimeShiftCCDBPathNegMC, task);
  getCfg(initContext, "paramFileName", mParamFileName, task);
  getCfg(initContext, "parametrizationPath", mParametrizationPath, task);
  getCfg(initContext, "reconstructionPass", mReconstructionPass, task);
  getCfg(initContext, "reconstructionPassDefault", mReconstructionPassDefault, task);
  getCfg(initContext, "fatalOnPassNotAvailable", mFatalOnPassNotAvailable, task);
  getCfg(initContext, "enableTimeDependentResponse", mEnableTimeDependentResponse, task);
  getCfg(initContext, "collisionSystem", mCollisionSystem, task);
  getCfg(initContext, "autoSetProcessFunctions", mAutoSetProcessFunctions, task);
}

void o2::pid::tof::TOFResponseImpl::initSetup(o2::ccdb::BasicCCDBManager* ccdb,
                                              o2::framework::InitContext& initContext)
{
  if (mIsInit) {
    LOG(fatal) << "TOFResponseImpl already initialized, cannot re-initialize";
  }

  if (!ccdb) {
    LOG(fatal) << "CCDB manager is not set, cannot initialize TOFResponseImpl";
  }
  inheritFromBaseTask(initContext); // Gets the configuration parameters from the base task (tof-signal)
  mCcdb = ccdb;                     // Set the CCDB manager
  mCcdb->setURL(mUrl);
  mCcdb->setTimestamp(mTimestamp);
  mCcdb->setCaching(true);
  mCcdb->setLocalObjectValidityChecking();
  // Not later than now objects
  mCcdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());

  mIsInit = true; // Set the initialization flag

  // Then the information about the metadata
  if (mReconstructionPass == "metadata") {
    LOG(info) << "Getting pass from metadata";
    if (metadataInfo.isMC()) {
      mReconstructionPass = metadataInfo.get("AnchorPassName");
    } else {
      mReconstructionPass = metadataInfo.get("RecoPassName");
    }
    LOG(info) << "Passed autodetect mode for pass. Taking '" << mReconstructionPass << "'";
  }
  LOG(info) << "Using parameter collection, starting from pass '" << mReconstructionPass << "'";

  if (!mParamFileName.empty()) { // Loading the parametrization from file
    LOG(info) << "Loading exp. sigma parametrization from file " << mParamFileName << ", using param: " << mParametrizationPath << " and pass " << mReconstructionPass;
    o2::tof::ParameterCollection paramCollection;
    paramCollection.loadParamFromFile(mParamFileName, mParametrizationPath);
    LOG(info) << "+++ Loaded parameter collection from file +++";
    if (!paramCollection.retrieveParameters(parameters, mReconstructionPass)) {
      if (mFatalOnPassNotAvailable) {
        LOG(fatal) << "Pass '" << mReconstructionPass << "' not available in the retrieved object from file";
      } else {
        LOG(warning) << "Pass '" << mReconstructionPass << "' not available in the retrieved object from file, fetching '" << mReconstructionPassDefault << "'";
        if (!paramCollection.retrieveParameters(parameters, mReconstructionPassDefault)) {
          paramCollection.print();
          LOG(fatal) << "Cannot get default pass for calibration " << mReconstructionPassDefault;
        } else {
          if (metadataInfo.isRun3()) {
            parameters.setResolutionParametrization(paramCollection.getPars(mReconstructionPassDefault));
          } else {
            parameters.setResolutionParametrizationRun2(paramCollection.getPars(mReconstructionPassDefault));
          }
          parameters.setMomentumChargeShiftParameters(paramCollection.getPars(mReconstructionPassDefault));
        }
      }
    } else { // Pass is available, load non standard parameters
      if (metadataInfo.isRun3()) {
        parameters.setResolutionParametrization(paramCollection.getPars(mReconstructionPass));
      } else {
        parameters.setResolutionParametrizationRun2(paramCollection.getPars(mReconstructionPass));
      }
      parameters.setMomentumChargeShiftParameters(paramCollection.getPars(mReconstructionPass));
    }
  } else if (!mEnableTimeDependentResponse) { // Loading it from CCDB
    LOG(info) << "Loading initial exp. sigma parametrization from CCDB, using path: " << mParametrizationPath << " for timestamp " << mTimestamp;
    o2::tof::ParameterCollection* paramCollection = mCcdb->getSpecific<o2::tof::ParameterCollection>(mParametrizationPath, mTimestamp);
    if (!paramCollection->retrieveParameters(parameters, mReconstructionPass)) { // Attempt at loading the parameters with the pass defined
      if (mFatalOnPassNotAvailable) {
        LOG(fatal) << "Pass '" << mReconstructionPass << "' not available in the retrieved CCDB object";
      } else {
        LOG(warning) << "Pass '" << mReconstructionPass << "' not available in the retrieved CCDB object, fetching '" << mReconstructionPassDefault << "'";
        if (!paramCollection->retrieveParameters(parameters, mReconstructionPassDefault)) {
          paramCollection->print();
          LOG(fatal) << "Cannot get default pass for calibration " << mReconstructionPassDefault;
        } else {
          if (metadataInfo.isRun3()) {
            parameters.setResolutionParametrization(paramCollection->getPars(mReconstructionPassDefault));
          } else {
            parameters.setResolutionParametrizationRun2(paramCollection->getPars(mReconstructionPassDefault));
          }
          parameters.setMomentumChargeShiftParameters(paramCollection->getPars(mReconstructionPassDefault));
        }
      }
    } else { // Pass is available, load non standard parameters
      if (metadataInfo.isRun3()) {
        parameters.setResolutionParametrization(paramCollection->getPars(mReconstructionPass));
      } else {
        parameters.setResolutionParametrizationRun2(paramCollection->getPars(mReconstructionPass));
      }
      parameters.setMomentumChargeShiftParameters(paramCollection->getPars(mReconstructionPass));
    }
  }

  // Loading additional calibration objects
  std::map<std::string, std::string> metadata;
  if (!mReconstructionPass.empty()) {
    metadata["RecoPassName"] = mReconstructionPass;
  }

  auto updateTimeShift = [&](const std::string& nameShift, bool isPositive) {
    if (nameShift.empty()) {
      return;
    }
    const bool isFromFile = nameShift.find(".root") != std::string::npos;
    if (isFromFile) {
      LOG(info) << "Initializing the time shift for " << (isPositive ? "positive" : "negative") << " from file '" << nameShift << "'";
      parameters.setTimeShiftParameters(nameShift, "ccdb_object", isPositive);
    } else if (!mEnableTimeDependentResponse) { // If the response is fixed fetch it at the init time
      LOG(info) << "Initializing the time shift for " << (isPositive ? "positive" : "negative")
                << " from ccdb '" << nameShift << "' and timestamp " << mTimestamp
                << " and pass '" << mReconstructionPass << "'";
      mCcdb->setFatalWhenNull(false);
      parameters.setTimeShiftParameters(mCcdb->getSpecific<TGraph>(nameShift, mTimestamp, metadata), isPositive);
      mCcdb->setFatalWhenNull(true);
    }
    LOG(info) << " test getTimeShift at 0 " << (isPositive ? "pos" : "neg") << ": "
              << parameters.getTimeShift(0, isPositive);
  };

  const std::string nameShiftPos = metadataInfo.isMC() ? mTimeShiftCCDBPathPosMC : mTimeShiftCCDBPathPos;
  updateTimeShift(nameShiftPos, true);
  const std::string nameShiftNeg = metadataInfo.isMC() ? mTimeShiftCCDBPathNegMC : mTimeShiftCCDBPathNeg;
  updateTimeShift(nameShiftNeg, false);

  // Calibration object is defined
  LOG(info) << "Parametrization at init time:";
  parameters.printFullConfig();
}

void o2::pid::tof::TOFResponseImpl::processSetup(const int runNumber, const int64_t timeStamp)
{
  LOG(debug) << "Processing setup for run number " << runNumber << " from run " << mLastRunNumber;
  // First we check if this run number was already processed
  if (mLastRunNumber == runNumber) {
    return;
  }
  LOG(info) << "Updating the parametrization from last run " << mLastRunNumber << " to " << runNumber << " and timestamp from " << mTimestamp << " " << timeStamp;
  mLastRunNumber = runNumber;
  mTimestamp = timeStamp;

  // Check the beam type
  if (mCollisionSystem == o2::common::core::CollisionSystemType::kCollSysUndef) {
    o2::parameters::GRPLHCIFData* grpo = mCcdb->getSpecific<o2::parameters::GRPLHCIFData>(mPathGrpLhcIf,
                                                                                          mTimestamp);
    mCollisionSystem = CollisionSystemType::getCollisionTypeFromGrp(grpo);
  } else {
    LOG(debug) << "Not setting collisions system as already set to " << mCollisionSystem << " " << CollisionSystemType::getCollisionSystemName(mCollisionSystem);
  }

  if (!mEnableTimeDependentResponse) {
    return;
  }
  LOG(info) << "Updating parametrization from path '" << mParametrizationPath << "' and timestamp " << mTimestamp << " and reconstruction pass '" << mReconstructionPass << "' for run number " << runNumber;
  if (mParamFileName.empty()) { // Not loading if parametrization was taken from file
    LOG(info) << "Updating parametrization from ccdb";
    const o2::tof::ParameterCollection* paramCollection = mCcdb->getSpecific<o2::tof::ParameterCollection>(mParametrizationPath, mTimestamp);
    if (!paramCollection->retrieveParameters(parameters, mReconstructionPass)) {
      if (mFatalOnPassNotAvailable) {
        LOGF(fatal, "Pass '%s' not available in the retrieved CCDB object", mReconstructionPass.data());
      } else {
        LOGF(warning, "Pass '%s' not available in the retrieved CCDB object, fetching '%s'", mReconstructionPass.data(), mReconstructionPassDefault.data());
        if (!paramCollection->retrieveParameters(parameters, mReconstructionPassDefault)) {
          paramCollection->print();
          LOG(fatal) << "Cannot get default pass for calibration " << mReconstructionPassDefault;
        } else { // Found the default case
          if (metadataInfo.isRun3()) {
            parameters.setResolutionParametrization(paramCollection->getPars(mReconstructionPassDefault));
          } else {
            parameters.setResolutionParametrizationRun2(paramCollection->getPars(mReconstructionPassDefault));
          }
          parameters.setMomentumChargeShiftParameters(paramCollection->getPars(mReconstructionPassDefault));
        }
      }
    } else { // Found the non default case
      if (metadataInfo.isRun3()) {
        parameters.setResolutionParametrization(paramCollection->getPars(mReconstructionPass));
      } else {
        parameters.setResolutionParametrizationRun2(paramCollection->getPars(mReconstructionPass));
      }
      parameters.setMomentumChargeShiftParameters(paramCollection->getPars(mReconstructionPass));
    }
  }

  // Loading additional calibration objects
  std::map<std::string, std::string> metadata;
  if (!mReconstructionPass.empty()) {
    metadata["RecoPassName"] = mReconstructionPass;
  }

  auto updateTimeShift = [&](const std::string& nameShift, bool isPositive) {
    if (nameShift.empty()) {
      return;
    }
    const bool isFromFile = nameShift.find(".root") != std::string::npos;
    if (isFromFile) {
      return;
    }
    LOG(info) << "Updating the time shift for " << (isPositive ? "positive" : "negative")
              << " from ccdb '" << nameShift << "' and timestamp " << mTimestamp
              << " and pass '" << mReconstructionPass << "'";
    mCcdb->setFatalWhenNull(false);
    parameters.setTimeShiftParameters(mCcdb->getSpecific<TGraph>(nameShift, mTimestamp, metadata), isPositive);
    mCcdb->setFatalWhenNull(true);
    LOG(info) << " test getTimeShift at 0 " << (isPositive ? "pos" : "neg") << ": "
              << parameters.getTimeShift(0, isPositive);
  };

  updateTimeShift(metadataInfo.isMC() ? mTimeShiftCCDBPathPosMC : mTimeShiftCCDBPathPos, true);
  updateTimeShift(metadataInfo.isMC() ? mTimeShiftCCDBPathNegMC : mTimeShiftCCDBPathNeg, false);

  LOG(info) << "Parametrization at setup time:";
  parameters.printFullConfig();
}

struct TOFSupport : o2::framework::ServicePlugin {
  o2::framework::ServiceSpec* create() final
  {
    return new ServiceSpec{
      .name = "tof-response",
      .init = [](ServiceRegistryRef, DeviceState&, fair::mq::ProgOptions&) -> ServiceHandle {
        auto* wrapper = new o2::pid::tof::TOFResponse();
        auto* ptr = new o2::pid::tof::TOFResponseImpl();
        wrapper->setInstance(ptr);
        return ServiceHandle{TypeIdHelpers::uniqueId<o2::pid::tof::TOFResponse>(), wrapper, ServiceKind::Serial, "database-pdg"};
      },
      .configure = CommonServices::noConfiguration(),
      .exit = [](ServiceRegistryRef, void* service) {
        auto* resp = reinterpret_cast<o2::pid::tof::TOFResponse*>(service);
        delete resp; },
      .kind = ServiceKind::Serial};
  }
};

DEFINE_DPL_PLUGINS_BEGIN
DEFINE_DPL_PLUGIN_INSTANCE(TOFSupport, CustomService);
DEFINE_DPL_PLUGINS_END
