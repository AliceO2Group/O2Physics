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

using framework::Service;

double ctpRateFetcher::fetch(Service<o2::ccdb::BasicCCDBManager>& ccdb, uint64_t timeStamp, int runNumber, std::string sourceName)
{
  if (sourceName.find("ZNC") != std::string::npos) {
    if (runNumber < 544448) {
      return fetchCTPratesInputs(ccdb, timeStamp, runNumber, 26) / (sourceName.find("hadronic") != std::string::npos ? 28. : 1.);
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
      return fetchCTPratesClasses(ccdb, timeStamp, runNumber, "CMTVX-B-NOPF-CRU");
    }
  }
  LOG(error) << "CTP rate for " << sourceName << " not available";
  return -1.;
}

double ctpRateFetcher::fetchCTPratesClasses(Service<o2::ccdb::BasicCCDBManager>& ccdb, uint64_t timeStamp, int runNumber, std::string className, int inputType)
{
  getCTPscalers(ccdb, timeStamp, runNumber);
  getCTPconfig(ccdb, timeStamp, runNumber);

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

  auto rate{mScalers->getRateGivenT(timeStamp, classIndex, inputType)};

  return pileUpCorrection(rate.second);
}

double ctpRateFetcher::fetchCTPratesInputs(Service<o2::ccdb::BasicCCDBManager>& ccdb, uint64_t timeStamp, int runNumber, int input)
{
  getCTPscalers(ccdb, timeStamp, runNumber);
  getLHCIFdata(ccdb, timeStamp, runNumber);

  std::vector<ctp::CTPScalerRecordO2> recs = mScalers->getScalerRecordO2();
  if (recs[0].scalersInps.size() == 48) {
    return pileUpCorrection(mScalers->getRateGivenT(timeStamp, input, 7).second);
  } else {
    LOG(error) << "Inputs not available";
    return -1.;
  }
}

void ctpRateFetcher::getCTPscalers(Service<o2::ccdb::BasicCCDBManager>& ccdb, uint64_t timeStamp, int runNumber)
{
  if (runNumber == mRunNumber && mScalers != nullptr) {
    return;
  }
  std::map<string, string> metadata;
  metadata["runNumber"] = std::to_string(runNumber);
  mScalers = ccdb->getSpecific<ctp::CTPRunScalers>("CTP/Calib/Scalers", timeStamp, metadata);
  if (mScalers == nullptr) {
    LOG(fatal) << "CTPRunScalers not in database, timestamp:" << timeStamp;
  }
  mScalers->convertRawToO2();
}

void ctpRateFetcher::getLHCIFdata(Service<o2::ccdb::BasicCCDBManager>& ccdb, uint64_t timeStamp, int runNumber)
{
  if (runNumber == mRunNumber && mLHCIFdata != nullptr) {
    return;
  }
  std::map<string, string> metadata;
  mLHCIFdata = ccdb->getSpecific<parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", timeStamp, metadata);
  if (mLHCIFdata == nullptr) {
    LOG(fatal) << "GRPLHCIFData not in database, timestamp:" << timeStamp;
  }
}

void ctpRateFetcher::getCTPconfig(Service<o2::ccdb::BasicCCDBManager>& ccdb, uint64_t timeStamp, int runNumber)
{
  if (runNumber == mRunNumber && mConfig != nullptr) {
    return;
  }
  std::map<string, string> metadata;
  metadata["runNumber"] = std::to_string(runNumber);
  mConfig = ccdb->getSpecific<ctp::CTPConfiguration>("CTP/Config/Config", timeStamp, metadata);
  if (mConfig == nullptr) {
    LOG(fatal) << "CTPRunConfig not in database, timestamp:" << timeStamp;
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

} // namespace o2
