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
#include "Common/DataModel/EventSelection.h"
#include "Framework/AnalysisDataModel.h"
#include "CommonConstants/LHCConstants.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/DataModel/FT0Corrected.h"

using namespace o2::aod;
struct FT0CorrectedTable {
  Produces<o2::aod::FT0sCorrected> table;
  using BCsWithMatchings = soa::Join<aod::BCs, aod::Run3MatchedToBCSparse>;
  using CollisionEvSel = soa::Join<aod::Collisions, aod::EvSels>::iterator;

  void process(BCsWithMatchings const& bcs, soa::Join<aod::Collisions, aod::EvSels> const& collisions, aod::FT0s const& ft0s)
  {
    for (auto& collision : collisions) {
      float vertexPV = collision.posZ();
      float vertex_corr = vertexPV / o2::constants::physics::LightSpeedCm2NS;
      float t0A = 1e10;
      float t0C = 1e10;
      if (collision.has_foundFT0()) {
        auto ft0 = collision.foundFT0();
        int triggersignals = ft0.triggerMask();
        bool ora = (triggersignals & (1 << 0)) != 0;
        bool orc = (triggersignals & (1 << 1)) != 0;
        LOGF(debug, "triggers OrA %i OrC %i ", ora, orc);
        LOGF(debug, " T0A = %f, T0C %f, vertex_corr %f, triggersignals %i", ft0.timeA(), ft0.timeC(), vertex_corr, triggersignals);
        if (ora) {
          t0A = ft0.timeA() + vertex_corr;
        }
        if (orc) {
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
