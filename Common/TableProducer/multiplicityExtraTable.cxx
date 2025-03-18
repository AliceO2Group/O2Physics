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
#include <CCDB/BasicCCDBManager.h> // megalinter thinks this is a C header...
#include <bitset>
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
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPLHCIFData.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using BCPattern = std::bitset<o2::constants::lhc::LHCMaxBunches>;
const int nBCsPerOrbit = o2::constants::lhc::LHCMaxBunches;

struct MultiplicityExtraTable {
  Produces<aod::MultBCs> multBC;
  Produces<aod::MultNeighs> multNeigh;

  Produces<aod::Mults2BC> mult2bc;
  Produces<aod::BC2Mults> bc2mult;

  // Allow for downscaling of BC table for less space use in derived data
  Configurable<float> bcDownscaleFactor{"bcDownscaleFactor", 2, "Downscale factor for BC table (0: save nothing, 1: save all)"};
  Configurable<float> minFT0CforBCTable{"minFT0CforBCTable", 25.0f, "Minimum FT0C amplitude to fill BC table to reduce data"};
  Configurable<bool> saveOnlyBCsWithCollisions{"saveOnlyBCsWithCollisions", true, "save only BCs with collisions in them"};

  Configurable<float> bcTableFloatPrecision{"bcTableFloatPrecision", 0.1, "float precision in bc table for data reduction"};

  float tru(float value)
  {
    if (bcTableFloatPrecision < 1e-4)
      return value; // make sure nothing bad happens in case zero (best precision)
    return bcTableFloatPrecision * std::round(value / bcTableFloatPrecision) + 0.5f * bcTableFloatPrecision;
  };

  // needed for downscale
  unsigned int randomSeed = 0;

  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  BCPattern CollidingBunch;

  int newRunNumber = -999;
  int oldRunNumber = -999;

  void init(InitContext&)
  {
    randomSeed = static_cast<unsigned int>(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());

    ccdbApi.init("http://alice-ccdb.cern.ch");
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
  }

  using BCsWithRun3Matchings = soa::Join<aod::BCs, aod::Timestamps, aod::Run3MatchedToBCSparse>;

  void processBCs(soa::Join<BCsWithRun3Matchings, aod::BCFlags> const& bcs, aod::FV0As const&, aod::FT0s const&, aod::FDDs const&, aod::Zdcs const&, soa::Join<aod::Collisions, aod::EvSels> const& collisions)
  {
    //+-<*>-+-<*>-+-<*>-+-<*>-+-<*>-+-<*>-+-<*>-+-<*>-+-<*>-+-<*>-+
    // determine saved BCs and corresponding new BC table index
    std::vector<int> bcHasCollision(bcs.size());
    std::vector<int> newBCindex(bcs.size());
    std::vector<int> bc2multArray(bcs.size());
    int atIndex = 0;
    for (const auto& bc : bcs) {
      bcHasCollision[bc.globalIndex()] = false;
      newBCindex[bc.globalIndex()] = -1;
      bc2multArray[bc.globalIndex()] = -1;
    }

    //+-<*>-+-<*>-+-<*>-+-<*>-+-<*>-+-<*>-+-<*>-+-<*>-+-<*>-+-<*>-+
    // tag BCs that have a collision (from evsel foundBC)
    for (const auto& collision : collisions) {
      bcHasCollision[collision.foundBCId()] = true;
    }
    //+-<*>-+-<*>-+-<*>-+-<*>-+-<*>-+-<*>-+-<*>-+-<*>-+-<*>-+-<*>-+

    for (const auto& bc : bcs) {
      // downscale if requested to do so
      if (bcDownscaleFactor < 1.f && (static_cast<float>(rand_r(&randomSeed)) / static_cast<float>(RAND_MAX)) > bcDownscaleFactor) {
        continue;
      }

      float multFT0C = 0.f;
      if (bc.has_ft0()) {
        auto ft0 = bc.ft0();
        for (auto amplitude : ft0.amplitudeC()) {
          multFT0C += amplitude;
        }
      } else {
        multFT0C = -999.0f;
      }

      if (multFT0C < minFT0CforBCTable) {
        continue; // skip this event
      }

      if (saveOnlyBCsWithCollisions && !bcHasCollision[bc.globalIndex()]) {
        continue; // skip if no collision is assigned to this BC (from evSel assignment)
      }

      newBCindex[bc.globalIndex()] = atIndex++;
    }
    //+-<*>-+-<*>-+-<*>-+-<*>-+-<*>-+-<*>-+-<*>-+-<*>-+-<*>-+-<*>-+

    //+-<*>-+-<*>-+-<*>-+-<*>-+-<*>-+-<*>-+-<*>-+-<*>-+-<*>-+-<*>-+
    // interlink: collision -> valid BC, BC -> collision
    for (const auto& collision : collisions) {
      mult2bc(newBCindex[collision.foundBCId()]);
      bc2multArray[collision.foundBCId()] = collision.globalIndex();
    }
    //+-<*>-+-<*>-+-<*>-+-<*>-+-<*>-+-<*>-+-<*>-+-<*>-+-<*>-+-<*>-+

    for (const auto& bc : bcs) {
      if (newBCindex[bc.globalIndex()] < 0) {
        continue; // don't keep if low mult or downsampled out
      }

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

      float posZFT0 = -1e+3;
      bool posZFT0valid = false;

      uint8_t multFT0TriggerBits = 0;
      uint8_t multFV0TriggerBits = 0;
      uint8_t multFDDTriggerBits = 0;
      uint64_t multBCTriggerMask = bc.triggerMask();

      // initialize - from Arvind
      newRunNumber = bc.runNumber();
      int localBC = bc.globalBC() % nBCsPerOrbit;

      if (newRunNumber != oldRunNumber) {
        auto soreor = o2::ccdb::BasicCCDBManager::getRunDuration(ccdbApi, newRunNumber);
        auto ts = soreor.first;

        LOG(info) << " newRunNumber  " << newRunNumber << " time stamp " << ts;
        oldRunNumber = newRunNumber;
        auto grplhcif = ccdb->getForTimeStamp<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", ts);
        CollidingBunch = grplhcif->getBunchFilling().getBCPattern();
      } // new run number

      bool collidingBC = CollidingBunch.test(localBC);

      if (bc.has_ft0()) {
        const auto& ft0 = bc.ft0();
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
        posZFT0 = ft0.posZ();
        posZFT0valid = ft0.isValidTime();
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

      bc2mult(bc2multArray[bc.globalIndex()]);
      multBC(
        tru(multFT0A), tru(multFT0C),
        tru(posZFT0), posZFT0valid, tru(multFV0A),
        tru(multFDDA), tru(multFDDC), tru(multZNA), tru(multZNC), tru(multZEM1),
        tru(multZEM2), tru(multZPA), tru(multZPC), Tvx, isFV0OrA,
        multFV0TriggerBits, multFT0TriggerBits, multFDDTriggerBits, multBCTriggerMask, collidingBC,
        bc.timestamp(),
        bc.flags());
    }
  }

  void processCollisionNeighbors(aod::Collisions const& collisions)
  {
    std::vector<float> timeArray;
    timeArray.resize(collisions.size(), 1e+3);

    for (const auto& collision : collisions) {
      timeArray[collision.globalIndex()] = collision.collisionTime();
    }

    float deltaPrevious = 1e+6, deltaPrePrevious = 1e+6;
    float deltaNext = 1e+6, deltaNeNext = 1e+6;
    for (const auto& collision : collisions) {
      int ii = collision.globalIndex();

      if (ii - 1 >= 0)
        deltaPrevious = timeArray[ii] - timeArray[ii - 1];
      if (ii - 2 >= 0)
        deltaPrePrevious = timeArray[ii] - timeArray[ii - 2];
      if (ii + 1 < collisions.size())
        deltaNext = timeArray[ii + 1] - timeArray[ii];
      if (ii + 2 < collisions.size())
        deltaNeNext = timeArray[ii + 2] - timeArray[ii];

      multNeigh(deltaPrePrevious, deltaPrevious, deltaNext, deltaNeNext);
    }
  }

  // Process switches
  PROCESS_SWITCH(MultiplicityExtraTable, processBCs, "Produce BC tables", true);
  PROCESS_SWITCH(MultiplicityExtraTable, processCollisionNeighbors, "Produce neighbor timing tables", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<MultiplicityExtraTable>(cfgc, TaskName{"multiplicity-extra-table"})};
}
