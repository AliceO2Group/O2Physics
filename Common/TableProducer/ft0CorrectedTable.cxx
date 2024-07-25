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

#include <bitset>
#include "Common/DataModel/FT0Corrected.h"
#include "Framework/ConfigParamSpec.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/EventSelection.h"
#include "Framework/AnalysisDataModel.h"
#include "CommonConstants/LHCConstants.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsFT0/Digit.h"

using namespace o2;
using namespace o2::framework;

using namespace o2::aod;
struct FT0CorrectedTable {
  Produces<o2::aod::FT0sCorrected> table;
  using BCsWithMatchings = soa::Join<aod::BCs, aod::Run3MatchedToBCSparse>;
  using CollisionEvSel = soa::Join<aod::Collisions, aod::EvSels>::iterator;

  void process(BCsWithMatchings const&, soa::Join<aod::Collisions, aod::EvSels> const& collisions, aod::FT0s const&)
  {
    for (auto& collision : collisions) {
      float vertexPV = collision.posZ();
      float vertex_corr = vertexPV / o2::constants::physics::LightSpeedCm2NS;
      float t0A = 1e10;
      float t0C = 1e10;
      constexpr float dummyTime = 30.; // Due to HW limitations time can be only within range (-25,25) ns, dummy time is around 32 ns
      if (collision.has_foundFT0()) {
        auto ft0 = collision.foundFT0();
        std::bitset<8> triggers = ft0.triggerMask();
        bool ora = triggers[o2::ft0::Triggers::bitA];
        bool orc = triggers[o2::ft0::Triggers::bitC];
        LOGF(debug, "triggers OrA %i OrC %i ", ora, orc);
        LOGF(debug, " T0A = %f, T0C %f, vertex_corr %f", ft0.timeA(), ft0.timeC(), vertex_corr);
        if (ora && ft0.timeA() < dummyTime) {
          t0A = ft0.timeA() + vertex_corr;
        }
        if (orc && ft0.timeC() < dummyTime) {
          t0C = ft0.timeC() - vertex_corr;
        }
      }
      LOGF(debug, " T0 collision time T0A = %f, T0C = %f", t0A, t0C);
      table(t0A, t0C);
    }
  }
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<FT0CorrectedTable>(cfgc, TaskName{"ft0-corrected-table"})};
}
