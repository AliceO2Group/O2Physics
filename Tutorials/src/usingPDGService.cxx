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

#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;

struct UsePdgDatabase {
  Service<o2::framework::O2DatabasePDG> pdg;
  OutputObj<TH1F> particleCharges{TH1F("charges", ";charge;entries", 201, -10.1, 10.1)};
  OutputObj<TH1F> particleMasses{TH1F("masses", ";mass (GeV/#it{c}^{2});entries", 1000, 0., 10.)};

  void process(aod::McCollision const&, aod::McParticles const& particles)
  {
    for (auto const& particle : particles) {
      auto pdgCode = particle.pdgCode();
      // Using the o2::framework::O2DatabasePDGImpl features
      auto mass = pdg->Mass(pdgCode);
      if (mass >= 0.) {
        particleMasses->Fill(mass);
      }
      // Using the TDatabasePDG features
      auto pdgInfo = pdg->GetParticle(pdgCode);
      if (pdgInfo != nullptr) {
        particleCharges->Fill(pdgInfo->Charge());
      } else {
        LOGF(warn, "[%d] unknown particle with PDG code %d", particle.globalIndex(), pdgCode);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UsePdgDatabase>(cfgc)};
}
