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
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/runDataProcessing.h>

#include <vector>

using namespace o2;
using namespace o2::framework;

// Converts MCParticle table from version 000 to 001

struct McConverter {
  Produces<aod::StoredMcParticles_001> mcParticles_001;

  void process(aod::StoredMcParticles_000 const& mcParticles_000)
  {
    for (auto& p : mcParticles_000) {

      std::vector<int> mothers;
      if (p.mother0Id() >= 0) {
        mothers.push_back(p.mother0Id());
      }
      if (p.mother1Id() >= 0) {
        mothers.push_back(p.mother1Id());
      }

      int daughters[2] = {-1, -1};
      if (p.daughter0Id() >= 0 && p.daughter1Id() >= 0) {
        daughters[0] = p.daughter0Id();
        daughters[1] = p.daughter1Id();
      } else if (p.daughter0Id() >= 0) {
        daughters[0] = p.daughter0Id();
        daughters[1] = p.daughter0Id();
      }

      mcParticles_001(p.mcCollisionId(), p.pdgCode(), p.statusCode(), p.flags(),
                      mothers, daughters, p.weight(), p.px(), p.py(), p.pz(), p.e(),
                      p.vx(), p.vy(), p.vz(), p.vt());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<McConverter>(cfgc),
  };
}
