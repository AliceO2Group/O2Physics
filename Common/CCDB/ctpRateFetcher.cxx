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

#include "ctpRateFetcher.h"

#include <map>
#include <vector>

#include "CommonConstants/LHCConstants.h"
#include "DataFormatsCTP/Configuration.h"
#include "DataFormatsCTP/Scalers.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "CCDB/BasicCCDBManager.h"

namespace o2
{
double ctpRateFetcher::fetch(o2::ccdb::BasicCCDBManager* ccdb, uint64_t timeStamp, int runNumber, std::string sourceName)
{
  setupRun(runNumber, ccdb, timeStamp);
  if (sourceName.find("ZNC") != std::string::npos) {
    if (runNumber < 544448) {
      return fetchCTPratesInputs(ccdb, timeStamp, runNumber, 25) / (sourceName.find("hadronic") != std::string::npos ? 28. : 1.);
    } else {
      return fetchCTPratesClasses(ccdb, timeStamp, runNumber, "C1ZNC-B-NOPF-CRU", 6) / (sourceName.find("hadronic") != std::string::npos ? 28. : 1.);
    }
  } else if (sourceName == "T0CE") {
    return fetchCTPratesClasses(ccdb, timeStamp, runNumber, "CMTVXTCE-B-NOPF-CRU");
  } else if (sourceName == "T0SC") {
    return fetchCTPratesClasses(ccdb, timeStamp, runNumber, "CMTVXTSC-B-NOPF-CRU");
  } else if (sourceName == "T0VTX") {
    if (runNumber < 534202) {
      return fetchCTPratesClasses(ccdb, timeStamp, runNumber, "minbias_TVX_L0"); // 2022
    } else {
      return fetchCTPratesClasses(ccdb, timeStamp, runNumber, "CMTVX-B-NOPF");
    }
  }
  LOG(error) << "CTP rate for " << sourceName << " not available";
  return -1.;
}

double ctpRateFetcher::fetchCTPratesClasses(o2::ccdb::BasicCCDBManager* ccdb, uint64_t timeStamp, int runNumber, const std::string& className, int inputType)
{
  std::vector<ctp::CTPClass> ctpcls = mConfig->getCTPClasses();
  std::vector<int> clslist = mConfig->getTriggerClassList();
  int classIndex = -1;
  for (size_t i = 0; i < clslist.size(); i++) {
    if (ctpcls[i].name == className) {
      classIndex = i;
      break;
    }
  }
  if (classIndex == -1) {
    LOG(fatal) << "Trigger class " << className << " not found in CTPConfiguration";
  }

  auto rate{mScalers->getRateGivenT(timeStamp * 1.e-3, classIndex, inputType)};

  return pileUpCorrection(rate.second);
}

double ctpRateFetcher::fetchCTPratesInputs(o2::ccdb::BasicCCDBManager* ccdb, uint64_t timeStamp, int runNumber, int input)
{
  std::vector<ctp::CTPScalerRecordO2> recs = mScalers->getScalerRecordO2();
  if (recs[0].scalersInps.size() == 48) {
    return pileUpCorrection(mScalers->getRateGivenT(timeStamp * 1.e-3, input, 7).second);
  } else {
    LOG(error) << "Inputs not available";
    return -1.;
  }
}

double ctpRateFetcher::pileUpCorrection(double triggerRate)
{
  auto bfilling = mLHCIFdata->getBunchFilling();
  std::vector<int> bcs = bfilling.getFilledBCs();
  double nbc = bcs.size();
  double nTriggersPerFilledBC = triggerRate / nbc / constants::lhc::LHCRevFreq;
  double mu = -std::log(1 - nTriggersPerFilledBC);
  return mu * nbc * constants::lhc::LHCRevFreq;
}

void ctpRateFetcher::setupRun(int runNumber, o2::ccdb::BasicCCDBManager* ccdb, uint64_t timeStamp)
{
  if (runNumber == mRunNumber) {
    return;
  }
  mRunNumber = runNumber;
  LOG(debug) << "Setting up CTP scalers for run " << mRunNumber;
  if (mManualCleanup) {
    delete mConfig;
    delete mScalers;
    delete mLHCIFdata;
  }
  std::map<string, string> metadata;
  mLHCIFdata = ccdb->getSpecific<parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", timeStamp, metadata);
  if (mLHCIFdata == nullptr) {
    LOG(fatal) << "GRPLHCIFData not in database, timestamp:" << timeStamp;
  }
  metadata["runNumber"] = std::to_string(mRunNumber);
  mConfig = ccdb->getSpecific<ctp::CTPConfiguration>("CTP/Config/Config", timeStamp, metadata);
  if (mConfig == nullptr) {
    LOG(fatal) << "CTPRunConfig not in database, timestamp:" << timeStamp;
  }
  mScalers = ccdb->getSpecific<ctp::CTPRunScalers>("CTP/Calib/Scalers", timeStamp, metadata);
  if (mScalers == nullptr) {
    LOG(fatal) << "CTPRunScalers not in database, timestamp:" << timeStamp;
  }
  mScalers->convertRawToO2();
}

} // namespace o2
