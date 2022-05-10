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

/// \brief write relevant information for photon conversion analysis to an AO2D.root file. This file is then the only necessary input to perform
/// pcm analysis.
/// dependencies: o2-analysis-lf-lambdakzerobuilder
/// \author stephan.friedrich.stiefelmaier@cern.ch

// runme like: o2-analysis-trackselection -b --aod-file ${sourceFile} --aod-writer-json ${writerFile} | o2-analysis-timestamp -b | o2-analysis-trackextension -b | o2-analysis-lf-lambdakzerobuilder -b | o2-analysis-pid-tpc -b | o2-analysis-em-skimmermc -b

// todo: remove reduantant information in GammaConversionsInfoTrue
#include "gammaTables.h"

#include "TVector3.h"

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct skimmerGammaConversionsTruthOnlyMc {

  Produces<aod::McGammasTrue> fFuncTableMcGammas;
  Produces<aod::McGammaDaughtersTrue> fFuncTableMcGammaDaughters;

  HistogramRegistry registry{
    "registry",
    {
      {"hCollisionZ_MCRec", "hCollisionZ_MCRec", {HistType::kTH1F, {{800, -10.f, 10.f}}}},
      {"hCollisionZ_all_MCTrue", "hCollisionZ_all_MCTrue", {HistType::kTH1F, {{800, -10.f, 10.f}}}},
      {"hCollisionZ_MCTrue", "hCollisionZ_MCTrue", {HistType::kTH1F, {{800, -10.f, 10.f}}}},
      {"hMcParticlesSize", "hMcParticlesSize", {HistType::kTH1F, {{100, 0.f, 1000000.f}}}},
      {"hEtaDiff", "hEtaDiff", {HistType::kTH1F, {{400, -2.f, 2.f}}}},
    },
  };

  void process(aod::McCollision const& theMcCollision,
               soa::SmallGroups<soa::Join<aod::McCollisionLabels,
                                          aod::Collisions>> const& theCollisions,
               aod::McParticles const& theMcParticles)
  {
    registry.fill(HIST("hCollisionZ_all_MCTrue"), theMcCollision.posZ());
    if (theCollisions.size() == 0) {
      return;
    }
    registry.fill(HIST("hCollisionZ_MCTrue"), theMcCollision.posZ());
    registry.fill(HIST("hMcParticlesSize"), theMcParticles.size());

    for (auto& lCollision : theCollisions) {
      registry.fill(HIST("hCollisionZ_MCRec"), lCollision.posZ());
    }

    for (auto& lMcParticle : theMcParticles) {
      if (lMcParticle.pdgCode() == 22) {

        size_t lNDaughters = 0;
        float lDaughter0Vx = -1.;
        float lDaughter0Vy = -1.;
        float lDaughter0Vz = -1.;
        float lV0Radius = -1.;

        if (lMcParticle.has_daughters()) {
          auto lDaughters = lMcParticle.daughters_as<aod::McParticles>();
          lNDaughters = lDaughters.size();
          auto lDaughter0 = lDaughters.begin();
          lDaughter0Vx = lDaughter0.vx();
          lDaughter0Vy = lDaughter0.vy();
          lDaughter0Vz = lDaughter0.vz();
          lV0Radius = sqrt(pow(lDaughter0Vx, 2) + pow(lDaughter0Vy, 2));

          for (auto& lDaughter : lMcParticle.daughters_as<aod::McParticles>()) {

            // SFS this can be removed once the eta integrity is checked
            TVector3 lDaughterVtx(lDaughter.vx(), lDaughter.vy(), lDaughter.vz());
            if (lMcParticle.isPhysicalPrimary()) {
              float_t lEtaDiff = lDaughterVtx.Eta() - lMcParticle.eta();
              registry.fill(HIST("hEtaDiff"), lEtaDiff);
            }
            fFuncTableMcGammaDaughters(lDaughter.mcCollisionId(),
                                       lMcParticle.globalIndex(),
                                       lDaughter.mothersIds().size(),
                                       lDaughter.pdgCode(), lDaughter.statusCode(), lDaughter.flags(),
                                       lDaughter.px(), lDaughter.py(), lDaughter.pz(), lDaughter.e(),
                                       lDaughter.vx(), lDaughter.vy(), lDaughter.vz(), lDaughter.vt());
          }
        }
        fFuncTableMcGammas(
          lMcParticle.mcCollisionId(),
          lMcParticle.globalIndex(),
          -1, // V0Id when running in reconstructed task
          lMcParticle.statusCode(),
          lMcParticle.flags(),
          lMcParticle.px(), lMcParticle.py(), lMcParticle.pz(),
          lMcParticle.vx(), lMcParticle.vy(), lMcParticle.vz(), lMcParticle.vt(),
          lNDaughters,
          lMcParticle.eta(), lMcParticle.phi(), lMcParticle.p(), lMcParticle.pt(), lMcParticle.y(),
          lDaughter0Vx, lDaughter0Vy, lDaughter0Vz,
          lV0Radius);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<skimmerGammaConversionsTruthOnlyMc>(cfgc)};
}
