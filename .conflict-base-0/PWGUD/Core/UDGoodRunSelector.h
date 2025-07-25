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

#ifndef PWGUD_CORE_UDGOODRUNSELECTOR_H_
#define PWGUD_CORE_UDGOODRUNSELECTOR_H_

#include <map>
#include <string>
#include <vector>

// A class to select good runs
struct UDGoodRunSelector {

 public:
  // constructor
  UDGoodRunSelector() {}
  explicit UDGoodRunSelector(std::string const& goodRunsFile);
  ~UDGoodRunSelector();

  // setters
  bool init(std::string const& goodRunsFile);
  void clear();

  // getters
  void Print();
  bool isGoodRun(int runNumber);
  std::vector<int> goodRuns() { return mgoodRuns; }
  std::vector<int> goodRuns(std::string runPeriod);
  int rnumMin() { return mrnMin; }
  int rnumMax() { return mrnMax; }

 private:
  bool misActive;
  std::string mgoodRunsFile;
  int mrnMin = -1, mrnMax = -1;
  std::vector<int> mgoodRuns;
  std::map<std::string, std::vector<int>> mrunMap;
};

#endif // PWGUD_CORE_UDGOODRUNSELECTOR_H_
