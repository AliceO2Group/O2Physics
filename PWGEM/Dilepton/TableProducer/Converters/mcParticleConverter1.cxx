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
//
// ========================
//
// This code produces emmctable table 001 from 000.
//    Please write to: daiki.sekihata@cern.ch

#include "PWGEM/Dilepton/DataModel/dileptonTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct mcParticleConverter1 {
  Produces<aod::EMMCParticles_001> mcParticle_001;

  void process(aod::EMMCParticles_000 const& mcParticles)
  {
    for (const auto& mcParticle : mcParticles) {
      // LOGF(info, "mcParticles.emmceventId() = %d, mcParticle.mothersIds().size() = %d, mcParticle.daughtersIds().size() = %d", mcParticle.emmceventId(), mcParticle.mothersIds().size(), mcParticle.daughtersIds().size());

      std::vector<int> mothersIds;
      for (const auto& id : mcParticle.mothersIds()) {
        mothersIds.emplace_back(id);
      }

      std::vector<int> daughtersIds;
      for (const auto& id : mcParticle.daughtersIds()) {
        daughtersIds.emplace_back(id);
      }

      mcParticle_001(
        mcParticle.emmceventId(), mcParticle.pdgCode(), mcParticle.flags(), 0,
        mothersIds, daughtersIds,
        mcParticle.px(), mcParticle.py(), mcParticle.pz(), mcParticle.e(),
        mcParticle.vx(), mcParticle.vy(), mcParticle.vz());
    } // end of mc particle loop
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<mcParticleConverter1>(cfgc, TaskName{"mcparticle-converter1"})};
}
