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

#ifndef COMMON_CCDB_CTPRATEFETCHER_H_
#define COMMON_CCDB_CTPRATEFETCHER_H_

#include <string>

#include "CCDB/BasicCCDBManager.h"
#include "Framework/AnalysisHelpers.h"

namespace o2
{

namespace ctp
{
class CTPRunScalers;
class CTPConfiguration;
} // namespace ctp

namespace parameters
{
class GRPLHCIFData;
}

class ctpRateFetcher
{
 public:
  ctpRateFetcher() = default;
  double fetch(framework::Service<o2::ccdb::BasicCCDBManager>& ccdb, uint64_t timeStamp, int runNumber, std::string sourceName);

 private:
  void getCTPconfig(framework::Service<o2::ccdb::BasicCCDBManager>& ccdb, uint64_t timeStamp, int runNumber);
  void getCTPscalers(framework::Service<o2::ccdb::BasicCCDBManager>& ccdb, uint64_t timeStamp, int runNumber);
  void getLHCIFdata(framework::Service<o2::ccdb::BasicCCDBManager>& ccdb, uint64_t timeStamp, int runNumber);
  double fetchCTPratesInputs(framework::Service<o2::ccdb::BasicCCDBManager>& ccdb, uint64_t timeStamp, int runNumber, int input);
  double fetchCTPratesClasses(framework::Service<o2::ccdb::BasicCCDBManager>& ccdb, uint64_t timeStamp, int runNumber, std::string className, int inputType = 1);
  double pileUpCorrection(double rate);

  int mRunNumber = -1;
  ctp::CTPConfiguration* mConfig = nullptr;
  ctp::CTPRunScalers* mScalers = nullptr;
  parameters::GRPLHCIFData* mLHCIFdata = nullptr;
};
} // namespace o2

#endif // COMMON_CCDB_CTPRATEFETCHER_H_
