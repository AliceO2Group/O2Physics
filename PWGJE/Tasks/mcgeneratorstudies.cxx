// Copyright 2019-2024 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

// Task that produces the generated pT spectrum of a given particle for MC studies
//
/// \author Nicolas Strangmann <nicolas.strangmann@cern.ch>, Goethe University Frankfurt / Oak Ridge National Laoratory

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/HistogramRegistry.h"

#include "PWGJE/DataModel/EMCALMatchedCollisions.h"

#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using MyMCCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::EMCALMatchedCollisions>;

struct MCGeneratorStudies {
  HistogramRegistry mHistManager{"MCGeneratorStudyHistograms"};

  Configurable<float> mVertexCut{"vertexCut", 10.f, "apply z-vertex cut with value in cm"};
  Configurable<float> mRapidityCut{"rapidityCut", 0.9f, "Maximum absolute rapidity of counted generated particles"};
  Configurable<double> mSelectedParticleCode{"particlePDGCode", 111, "PDG code of the particle to be investigated"};

  Filter collisionFilter = (aod::collision::posZ > -mVertexCut) && (aod::collision::posZ < mVertexCut);
  Filter mcParticleFilter = (aod::mcparticle::pdgCode == mSelectedParticleCode) && (aod::mcparticle::y > -mRapidityCut) && (aod::mcparticle::y < mRapidityCut);

  void init(InitContext const&)
  {
    AxisSpec pTAxis{250, 0., 25., "#it{p}_{T} (GeV/#it{c})"};

    mHistManager.add("hCollisionCounter", "Number of collisions after event cuts", HistType::kTH1F, {{5, 0.5, 5.5}});
    mHistManager.add("hpT_all", "All collisions", HistType::kTH1F, {pTAxis});
    mHistManager.add("hpT_T0Triggered", "T0 triggered collisions", HistType::kTH1F, {pTAxis});
    mHistManager.add("hpT_T0Triggered_EMCReadout", "T0 triggered and EMC readout collisions", HistType::kTH1F, {pTAxis});
    mHistManager.add("hpT_T0Triggered_EMCReadout_Unique", "Unique T0 triggered and EMC readout collisions", HistType::kTH1F, {pTAxis});
  }

  PresliceUnsorted<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;

  void process(soa::Filtered<MyMCCollisions>::iterator const& collision, soa::Filtered<aod::McParticles> const& mcParticles)
  {
    bool isT0Triggered = collision.sel8();
    bool isEMCReadout = collision.isemcreadout();
    bool isUniqueCollision = !collision.ambiguous();

    mHistManager.fill(HIST("hCollisionCounter"), 1);
    if (isT0Triggered) {
      mHistManager.fill(HIST("hCollisionCounter"), 2);
      // mHistManager.fill<aod::mcparticle::Pt>(HIST("hpT_T0Triggered"), mcParticles, aod::evsel::sel8 == true);
      if (isEMCReadout) {
        mHistManager.fill(HIST("hCollisionCounter"), 3);
        if (isUniqueCollision) {
          mHistManager.fill(HIST("hCollisionCounter"), 4);
        }
      }
    }
    auto mcParticles_inColl = mcParticles.sliceBy(perMcCollision, collision.globalIndex());
    for (auto& mcParticle : mcParticles_inColl) {
      mHistManager.fill(HIST("hpT_all"), mcParticle.pt());
      if (isT0Triggered) {
        mHistManager.fill(HIST("hpT_T0Triggered"), mcParticle.pt());
        if (isEMCReadout) {
          mHistManager.fill(HIST("hpT_T0Triggered_EMCReadout"), mcParticle.pt());
          if (isUniqueCollision) {
            mHistManager.fill(HIST("hpT_T0Triggered_EMCReadout_Unique"), mcParticle.pt());
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<MCGeneratorStudies>(cfgc, TaskName{"mc-generator-studies"})};
}
