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
/// \file materialBudgetWeights.cxx
///
/// \brief This code produces a table to retrieve material budget weights. The table is to be join with V0PhotonKF
///
/// \author Youssef El Mard (youssef.el.mard.bouziani@cern.ch)
///

#ifndef PWGEM_PHOTONMESON_UTILS_MATERIALBUDGETWEIGHTS_H_
#define PWGEM_PHOTONMESON_UTILS_MATERIALBUDGETWEIGHTS_H_

#include <map>
#include <string>

#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Logger.h"
#include <CCDB/BasicCCDBManager.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>

using namespace o2;
using namespace o2::soa;
using namespace o2::framework;

using MyV0PhotonsMB = o2::soa::Join<o2::aod::V0PhotonsKF, o2::aod::V0KFEMEventIds>;
using MyV0PhotonMB = MyV0PhotonsMB::iterator;

struct MaterialBudgetWeights {
  Produces<aod::V0PhotonOmegaMBWeights> omegaMBWeight;

  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://ccdb-test.cern.ch:8080", "CCDB url"};
  Configurable<std::string> mbWeightsPath{"mbWeightsPath", "Users/y/yelmard/MaterialBudget/OmegaMBWeights", "Path of the mb weights"};
  Configurable<bool> runWithMBWeights{"runWithMBWeights", false, "Run task using material-budget weights for PCM photons"};

  o2::ccdb::CcdbApi ccdbApi;
  TH1F* hOmegaMBFromCCDB = nullptr;

  void init(InitContext&)
  {
    if (!runWithMBWeights) {
      LOG(info) << "MaterialBudgetWeights: runWithMBWeights=false -> no CCDB query, will write weight=1";
      return;
    }

    ccdbApi.init(ccdbUrl.value);
    std::map<std::string, std::string> metadata;
    LOG(info) << "MaterialBudgetWeights: loading Omega MB histogram from CCDB at path: " << mbWeightsPath.value;

    hOmegaMBFromCCDB = ccdbApi.retrieveFromTFileAny<TH1F>(mbWeightsPath, metadata, -1);

    if (!hOmegaMBFromCCDB) {
      LOG(fatal) << "MaterialBudgetWeights: runWithMBWeights=true but CCDB object is missing. Path=" << mbWeightsPath.value;
    }
  }

  float computeMBWeight(float v0Rxy)
  {
    if (!hOmegaMBFromCCDB) {
      LOG(fatal) << "MaterialBudgetWeights: internal error: histogram nullptr while runWithMBWeights=true";
    }

    int binMBWeight = hOmegaMBFromCCDB->FindBin(v0Rxy);
    if (binMBWeight < 1 || binMBWeight > hOmegaMBFromCCDB->GetNbinsX()) {
      LOG(debug) << "MaterialBudgetWeights: v0Rxy out of histogram range, returning 1";
      return 1.f;
    }

    return hOmegaMBFromCCDB->GetBinContent(binMBWeight);
  }

  void process(MyV0PhotonMB const& v0)
  {
    if (!runWithMBWeights) {
      omegaMBWeight(1.f);
      return;
    }

    omegaMBWeight(computeMBWeight(v0.v0radius()));
  }
};

#endif // PWGEM_PHOTONMESON_UTILS_MATERIALBUDGETWEIGHTS_H_
