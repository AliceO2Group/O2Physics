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

// task to produce a self contained data format for jet analyses from the full AO2D
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include <vector>
#include <string>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "TDatabasePDG.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetHFUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetV0DataProducerTask {

  Produces<aod::JV0McParticles> jV0McParticlesTable;
  Produces<aod::JV0McPIs> jV0McParentIndexTable;

  Service<o2::framework::O2DatabasePDG> pdgDatabase;

  void init(InitContext const&)
  {
  }

  void processV0MC(aod::McParticle const& particle)
  { // can loop over McV0Labels tables if we want to only store matched V0Particles
    if (jethfutilities::isV0Particle(particle)) {
      std::vector<int> mothersId;
      if (particle.has_mothers()) {
        auto mothersIdTemps = particle.mothersIds();
        for (auto mothersIdTemp : mothersIdTemps) {
          mothersId.push_back(mothersIdTemp);
        }
      }
      int daughtersId[2] = {-1, -1};
      auto i = 0;
      if (particle.has_daughters()) {
        for (auto daughterId : particle.daughtersIds()) {
          if (i > 1) {
            break;
          }
          daughtersId[i] = daughterId;
          i++;
        }
      }
      auto pdgParticle = pdgDatabase->GetParticle(particle.pdgCode());
      jV0McParticlesTable(particle.pt(), particle.eta(), particle.phi(), particle.y(), particle.e(), pdgParticle->Mass(), particle.pdgCode(), particle.getGenStatusCode(), particle.getHepMCStatusCode(), particle.isPhysicalPrimary(), mothersId, daughtersId, jetv0utilities::setV0ParticleDecayBit(particle));
      jV0McParentIndexTable(particle.mcCollisionId(), particle.globalIndex());
    }
  }
  PROCESS_SWITCH(JetV0DataProducerTask, processV0MC, "produces V0 particles", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetV0DataProducerTask>(cfgc, TaskName{"jet-data-producer-v0"})};
}
