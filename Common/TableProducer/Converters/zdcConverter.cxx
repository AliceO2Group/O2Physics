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

// ZDC converter to new format
// to be used with Run 2 converted data and older AO2Ds

#include <CommonConstants/ZDCConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/runDataProcessing.h>

#include <cstdint>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::zdc;

struct zdcConverter {
  Produces<aod::Zdcs_001> Zdcs_001;

  void process(aod::Zdcs_000 const& zdcLegacy, aod::BCs const&)
  {
    for (auto& zdcData : zdcLegacy) {
      // Get legacy information, please
      auto bc = zdcData.bc();
      auto energyZEM1 = zdcData.energyZEM1();
      auto energyZEM2 = zdcData.energyZEM2();
      auto energyCommonZNA = zdcData.energyCommonZNA();
      auto energyCommonZNC = zdcData.energyCommonZNC();
      auto energyCommonZPA = zdcData.energyCommonZPA();
      auto energyCommonZPC = zdcData.energyCommonZPC();
      auto energySectorZNA = zdcData.energySectorZNA();
      auto energySectorZNC = zdcData.energySectorZNC();
      auto energySectorZPA = zdcData.energySectorZPA();
      auto energySectorZPC = zdcData.energySectorZPC();
      auto timeZEM1 = zdcData.timeZEM1();
      auto timeZEM2 = zdcData.timeZEM2();
      auto timeZNA = zdcData.timeZNA();
      auto timeZNC = zdcData.timeZNC();
      auto timeZPA = zdcData.timeZPA();
      auto timeZPC = zdcData.timeZPC();

      // Create variables to initialize Zdcs_001 table
      std::vector<float> zdcEnergy, zdcAmplitudes, zdcTime;
      std::vector<uint8_t> zdcChannelsE, zdcChannelsT;

      // Tie variables in such that they get read correctly later
      zdcEnergy.emplace_back(energyZEM1);
      zdcChannelsE.emplace_back(IdZEM1);
      zdcAmplitudes.emplace_back(energyZEM1); // WARNING: DUMMY VALUE
      zdcTime.emplace_back(timeZEM1);
      zdcChannelsT.emplace_back(IdZEM1);

      zdcEnergy.emplace_back(energyZEM2);
      zdcChannelsE.emplace_back(IdZEM2);
      zdcAmplitudes.emplace_back(energyZEM2); // WARNING: DUMMY VALUE
      zdcTime.emplace_back(timeZEM2);
      zdcChannelsT.emplace_back(IdZEM2);

      zdcEnergy.emplace_back(energyCommonZNA);
      zdcChannelsE.emplace_back(IdZNAC);
      zdcAmplitudes.emplace_back(energyCommonZNA); // WARNING: DUMMY VALUE
      zdcTime.emplace_back(timeZNA);
      zdcChannelsT.emplace_back(IdZNAC);

      zdcEnergy.emplace_back(energyCommonZNC);
      zdcChannelsE.emplace_back(IdZNCC);
      zdcAmplitudes.emplace_back(energyCommonZNC); // WARNING: DUMMY VALUE
      zdcTime.emplace_back(timeZNC);
      zdcChannelsT.emplace_back(IdZNCC);

      zdcEnergy.emplace_back(energyCommonZPA);
      zdcChannelsE.emplace_back(IdZPAC);
      zdcAmplitudes.emplace_back(energyCommonZPA); // WARNING: DUMMY VALUE
      zdcTime.emplace_back(timeZPA);
      zdcChannelsT.emplace_back(IdZPAC);

      zdcEnergy.emplace_back(energyCommonZPC);
      zdcChannelsE.emplace_back(IdZPCC);
      zdcAmplitudes.emplace_back(energyCommonZPC); // WARNING: DUMMY VALUE
      zdcTime.emplace_back(timeZPC);
      zdcChannelsT.emplace_back(IdZPCC);

      for (uint64_t ic = 0; ic < 4; ic++) {
        zdcEnergy.emplace_back(energySectorZNA[ic]);
        zdcChannelsE.emplace_back(IdZNA1 + ic);
        zdcEnergy.emplace_back(energySectorZNC[ic]);
        zdcChannelsE.emplace_back(IdZNC1 + ic);
        zdcEnergy.emplace_back(energySectorZPA[ic]);
        zdcChannelsE.emplace_back(IdZPA1 + ic);
        zdcEnergy.emplace_back(energySectorZPC[ic]);
        zdcChannelsE.emplace_back(IdZPC1 + ic);
      }

      Zdcs_001(bc,
               zdcEnergy,
               zdcChannelsE,
               zdcAmplitudes,
               zdcTime,
               zdcChannelsT);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<zdcConverter>(cfgc),
  };
}
