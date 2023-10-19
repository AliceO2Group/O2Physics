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

using namespace o2;
using namespace o2::framework;

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include <CCDB/BasicCCDBManager.h>
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "DataFormatsFIT/Triggers.h"
#include "TableHelper.h"
#include "iostream"

struct MultiplicityExtraTable {
  Produces<aod::MultsBC> multBC;

  unsigned int randomSeed = 0;
  void init(InitContext& context)
  {
    // empty for now
  }

  using BCsWithRun3Matchings = soa::Join<aod::BCs, aod::Timestamps, aod::Run3MatchedToBCSparse>;

  void process(BCsWithRun3Matchings::iterator const& bc, aod::FV0As const&, aod::FT0s const& ft0s, aod::FDDs const&)
  {
    bool Tvx = false;
    bool isFV0OrA = false;
    float multFT0C = 0.f;
    float multFT0A = 0.f;
    float multFV0A = 0.f;

    if (bc.has_ft0()) {
      auto ft0 = bc.ft0();
      std::bitset<8> triggers = ft0.triggerMask();
      Tvx = triggers[o2::fit::Triggers::bitVertex];

      // calculate T0 charge
      for (auto amplitude : ft0.amplitudeA()) {
        multFT0A += amplitude;
      }
      for (auto amplitude : ft0.amplitudeC()) {
        multFT0C += amplitude;
      }

      if (bc.has_fv0a()) {
        auto fv0 = bc.fv0a();
        std::bitset<8> fV0Triggers = fv0.triggerMask();

        for (auto amplitude : fv0.amplitude()) {
          multFV0A += amplitude;
        }
        isFV0OrA = fV0Triggers[o2::fit::Triggers::bitA];
      } // fv0
    }

    multBC(multFT0A, multFT0C, multFV0A, Tvx, isFV0OrA);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<MultiplicityExtraTable>(cfgc, TaskName{"multiplicity-extra-table"})};
}
