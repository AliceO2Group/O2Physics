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

// This task finds minimum and maximum orbit among all processed bcs

#include "OrbitRange.h"

#include <CommonConstants/LHCConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/runDataProcessing.h>

#include <TMathBase.h>

#include <cstdint>
using namespace o2;
using namespace o2::framework;

struct OrbitRangeTask {
  OutputObj<OrbitRange> orbitRange{OrbitRange("orbitRange")};
  void process(aod::BC const& bc)
  {
    uint32_t orbit = bc.globalBC() / o2::constants::lhc::LHCMaxBunches;
    orbitRange->SetRunNumber(bc.runNumber());
    orbitRange->SetMinOrbit(TMath::Min(orbit, orbitRange->GetMinOrbit()));
    orbitRange->SetMaxOrbit(TMath::Max(orbit, orbitRange->GetMaxOrbit()));
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<OrbitRangeTask>(cfgc)};
}
