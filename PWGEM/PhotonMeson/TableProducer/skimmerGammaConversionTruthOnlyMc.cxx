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
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/gammaConvDefinitions.h"

#include "TVector3.h"

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct skimmerGammaConversionTruthOnlyMc {

  Produces<aod::McGammasTrue> fFuncTableMcGammas;
  Produces<aod::McDaughterTrue> fFuncTableMcDaughter;

  HistogramRegistry registry{
    "registry",
    {
      gHistoSpec_hCollisionZ_all_MCTrue,
      gHistoSpec_hCollisionZ_MCTrue,
      gHistoSpec_hCollisionZ_MCRec,
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
        int lastIndex0 = -1;
        int lastIndex1 = -1;

        if (lMcParticle.has_daughters()) {
          auto lDaughters = lMcParticle.daughters_as<aod::McParticles>();
          lNDaughters = lDaughters.size();
          auto lDaughter0 = lDaughters.begin();
          lDaughter0Vx = lDaughter0.vx();
          lDaughter0Vy = lDaughter0.vy();
          lDaughter0Vz = lDaughter0.vz();
          lV0Radius = sqrt(pow(lDaughter0Vx, 2) + pow(lDaughter0Vy, 2));

          if (lNDaughters == 2) {
            auto lDaughter1 = lDaughters.iteratorAt(1);
            fFuncTableMcDaughter(lDaughter0.p());
            lastIndex0 = fFuncTableMcDaughter.lastIndex();
            fFuncTableMcDaughter(lDaughter1.p());
            lastIndex1 = fFuncTableMcDaughter.lastIndex();
          }
        }
        fFuncTableMcGammas(
          lMcParticle.mcCollisionId(),
          lMcParticle.globalIndex(),
          -1, // V0Id when running in reconstructed task
          lMcParticle.pdgCode(), lMcParticle.statusCode(), lMcParticle.flags(),
          lMcParticle.px(), lMcParticle.py(), lMcParticle.pz(),
          lMcParticle.vx(), lMcParticle.vy(), lMcParticle.vz(), lMcParticle.vt(),
          lNDaughters,
          lMcParticle.eta(), lMcParticle.phi(), lMcParticle.p(), lMcParticle.pt(), lMcParticle.y(),
          lDaughter0Vx, lDaughter0Vy, lDaughter0Vz,
          lV0Radius,
          lastIndex0, lastIndex1);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<skimmerGammaConversionTruthOnlyMc>(cfgc, TaskName{"skimmer-gamma-conversion-truthonlymc"})};
}
