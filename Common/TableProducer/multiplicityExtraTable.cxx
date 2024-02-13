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
#include "Framework/ConfigParamSpec.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "DataFormatsFIT/Triggers.h"
#include "TableHelper.h"

#include "CCDB/CcdbApi.h"
#include "CommonDataFormat/BunchFilling.h"
#include <CCDB/BasicCCDBManager.h>
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include <bitset>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using BCPattern = std::bitset<o2::constants::lhc::LHCMaxBunches>;
const int nBCsPerOrbit = o2::constants::lhc::LHCMaxBunches;

struct MultiplicityExtraTable {
  Produces<aod::MultsBC> multBC;

  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  BCPattern CollidingBunch;

  int newRunNumber = -999;
  int oldRunNumber = -999;

  void init(InitContext& context)
  {
    ccdbApi.init("http://alice-ccdb.cern.ch");
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
  }

  using BCsWithRun3Matchings = soa::Join<aod::BCs, aod::Timestamps, aod::Run3MatchedToBCSparse>;

  void process(BCsWithRun3Matchings::iterator const& bc, aod::FV0As const&, aod::FT0s const& ft0s, aod::FDDs const&, aod::Zdcs const&)
  {
    bool Tvx = false;
    bool isFV0OrA = false;
    float multFT0C = 0.f;
    float multFT0A = 0.f;
    float multFV0A = 0.f;
    float multFDDA = 0.f;
    float multFDDC = 0.f;

    // ZDC amplitudes
    float multZEM1 = -1.f;
    float multZEM2 = -1.f;
    float multZNA = -1.f;
    float multZNC = -1.f;
    float multZPA = -1.f;
    float multZPC = -1.f;

    uint8_t multFT0TriggerBits = 0;
    uint8_t multFV0TriggerBits = 0;
    uint8_t multFDDTriggerBits = 0;
    uint64_t multBCTriggerMask = bc.triggerMask();

    // initialize - from Arvind
    newRunNumber = bc.runNumber();
    int localBC = bc.globalBC() % nBCsPerOrbit;

    if (newRunNumber != oldRunNumber) {
      uint64_t ts{};
      std::map<string, string> metadataRCT, headers;
      headers = ccdbApi.retrieveHeaders(Form("RCT/Info/RunInformation/%i", newRunNumber), metadataRCT, -1);
      ts = atol(headers["SOR"].c_str());

      LOG(info) << " newRunNumber  " << newRunNumber << " time stamp " << ts;
      oldRunNumber = newRunNumber;
      std::map<std::string, std::string> mapMetadata;
      std::map<std::string, std::string> mapHeader;
      auto grplhcif = ccdb->getForTimeStamp<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", ts);
      CollidingBunch = grplhcif->getBunchFilling().getBCPattern();
    } // new run number

    bool collidingBC = CollidingBunch.test(localBC);

    if (bc.has_ft0()) {
      auto ft0 = bc.ft0();
      std::bitset<8> triggers = ft0.triggerMask();
      Tvx = triggers[o2::fit::Triggers::bitVertex];
      multFT0TriggerBits = static_cast<uint8_t>(triggers.to_ulong());

      // calculate T0 charge
      for (auto amplitude : ft0.amplitudeA()) {
        multFT0A += amplitude;
      }
      for (auto amplitude : ft0.amplitudeC()) {
        multFT0C += amplitude;
      }
    } else {
      multFT0A = -999.0f;
      multFT0C = -999.0f;
    }
    if (bc.has_fv0a()) {
      auto fv0 = bc.fv0a();
      std::bitset<8> fV0Triggers = fv0.triggerMask();
      multFV0TriggerBits = static_cast<uint8_t>(fV0Triggers.to_ulong());

      for (auto amplitude : fv0.amplitude()) {
        multFV0A += amplitude;
      }
      isFV0OrA = fV0Triggers[o2::fit::Triggers::bitA];
    } else {
      multFV0A = -999.0f;
    }

    if (bc.has_fdd()) {
      auto fdd = bc.fdd();
      std::bitset<8> fFDDTriggers = fdd.triggerMask();
      multFDDTriggerBits = static_cast<uint8_t>(fFDDTriggers.to_ulong());

      for (auto amplitude : fdd.chargeA()) {
        multFDDA += amplitude;
      }
      for (auto amplitude : fdd.chargeC()) {
        multFDDC += amplitude;
      }
    } else {
      multFDDA = -999.0f;
      multFDDC = -999.0f;
    }

    if (bc.has_zdc()) {
      multZNA = bc.zdc().amplitudeZNA();
      multZNC = bc.zdc().amplitudeZNC();
      multZEM1 = bc.zdc().amplitudeZEM1();
      multZEM2 = bc.zdc().amplitudeZEM2();
      multZPA = bc.zdc().amplitudeZPA();
      multZPC = bc.zdc().amplitudeZPC();
    } else {
      multZNA = -999.f;
      multZNC = -999.f;
      multZEM1 = -999.f;
      multZEM2 = -999.f;
      multZPA = -999.f;
      multZPC = -999.f;
    }

    multBC(multFT0A, multFT0C, multFV0A, multFDDA, multFDDC, multZNA, multZNC, multZEM1, multZEM2, multZPA, multZPC, Tvx, isFV0OrA, multFV0TriggerBits, multFT0TriggerBits, multFDDTriggerBits, multBCTriggerMask, collidingBC);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<MultiplicityExtraTable>(cfgc, TaskName{"multiplicity-extra-table"})};
}
