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

/// \file dqFlowAccWeights.C
/// \brief A simple macro to read produced accptance weights for Q-vector and submit them to CCDB
///
/// \author Chi ZHANG, CEA-Saclay, chi.zhang@cern.ch

#include "PWGCF/GenericFramework/Core/GFWWeights.h"

#include "CCDB/BasicCCDBManager.h"
#include "Framework/Logger.h"

#include <TFile.h>
#include <TSystem.h>

#include <cmath>
#include <gsl/span>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace o2;
using namespace std;

void dqFlowAccWeights(int64 tmin = 1546300800000, int64 tmax = 1577833200000, std::string Period = "LHC23zzh_pass2", std::string SubDir = "d-q-event-qvector", std::string FileName = "AnalysisResults.root")
{
  if (tmax < tmin) {
    LOG(fatal) << "Wrong validity syntax!";
  }

  const std::string& ccdbHost = "http://alice-ccdb.cern.ch";
  // long tmin = 1546300800000; // Jan 2019
  // long tmax = 1577833200000; // Dec 2019
  // long tmin = 1641038770000; // Jan 2022
  // long tmax = 1663632000000; // 20 Sep 2022
  const std::string& objectPath = Form("Users/c/chizh/Acceptance/%s", Period.c_str());

  // Load weights
  TFile* file;
  try {
    file = TFile::Open(FileName.c_str());
  } catch (std::exception const& e) {
    LOG(fatal) << "Cannot open input file!";
  }

  GFWWeights* weights;
  try {
    file->GetObject(Form("%s/%s", SubDir.c_str(), "weights"), weights);
  } catch (std::exception const& e) {
    LOG(fatal) << "Cannot get GFWWeights from inpout file!";
  }

  // Send to CCDB
  if (!ccdbHost.empty()) {
    LOGP(info, "Storing alignment object on {}/{}", ccdbHost, objectPath);
    o2::ccdb::CcdbApi api;
    map<string, string> metadata; // can be empty
    metadata.insert(std::pair{"comment", Form("Acceptance weights for %s", Period.c_str())});
    api.init(ccdbHost.c_str());
    try {
      api.storeAsTFileAny(weights, objectPath, metadata, tmin, tmax);
    } catch (std::exception const& e) {
      LOG(fatal) << "Failed at CCDB submission!";
    }
  }
  file->Close();
}
