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

/// \file femtoDreamDebugV0.cxx
/// \brief Tasks that reads the particle tables and fills QA histograms for V0s
/// \author Luca Barioglio, TU München, luca.barioglio@cern.ch

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "DataFormatsParameters/GRPObject.h"

#include "PWGCF/DataModel/FemtoDerived.h"
#include "FemtoDreamParticleHisto.h"
#include "FemtoDreamEventHisto.h"
#include "FemtoUtils.h"

using namespace o2;
using namespace o2::analysis::femtoDream;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct femtoDreamDebugV0 {

  Configurable<int> ConfPDGCodePartOne{"ConfPDGCodePartOne", 3122, "Particle 1 - PDG code"};
  Configurable<uint32_t> ConfCutPartOne{"ConfCutPartOne", 338, "Particle 1 - Selection bit from cutCulator"};

  using FemtoFullParticles = soa::Join<aod::FemtoDreamParticles, aod::FemtoDreamDebugParticles>;

  Partition<FemtoFullParticles> partsOne = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kV0)) &&
                                           // (aod::femtodreamparticle::pt < cfgCutTable->get("MaxPt")) &&
                                           ((aod::femtodreamparticle::cut & ConfCutPartOne) == ConfCutPartOne);

  /// Histogramming for Event
  FemtoDreamEventHisto eventHisto;

  /// The configurables need to be passed to an std::vector
  std::vector<int> vPIDPartOne;

  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry FullQaRegistry{"FullV0QA", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    eventHisto.init(&qaRegistry);

    FullQaRegistry.add("FullV0QA/hPt", "; #it{p}_{T} (GeV/#it{c}); Entries", kTH1F, {{240, 0, 6}});
    FullQaRegistry.add("FullV0QA/hEta", "; #eta; Entries", kTH1F, {{200, -1.5, 1.5}});
    FullQaRegistry.add("FullV0QA/hPhi", "; #phi; Entries", kTH1F, {{200, 0, 2. * M_PI}});
    FullQaRegistry.add("FullV0QA/hDaughDCA", "; DCA^{daugh} (cm); Entries", kTH1F, {{1000, 0, 10}});
    FullQaRegistry.add("FullV0QA/hTransRadius", "; #it{r}_{xy} (cm); Entries", kTH1F, {{1500, 0, 150}});
    FullQaRegistry.add("FullV0QA/hDecayVtxX", "; #it{Vtx}_{x} (cm); Entries", kTH1F, {{2000, 0, 200}});
    FullQaRegistry.add("FullV0QA/hDecayVtxY", "; #it{Vtx}_{y} (cm)); Entries", kTH1F, {{2000, 0, 200}});
    FullQaRegistry.add("FullV0QA/hDecayVtxZ", "; #it{Vtx}_{z} (cm); Entries", kTH1F, {{2000, 0, 200}});
    FullQaRegistry.add("FullV0QA/hCPA", "; #it{cos #theta_{p}}; Entries", kTH1F, {{1000, 0.9, 1.}});
    FullQaRegistry.add("FullV0QA/hCPAvsPt", "; #it{p}_{T} (GeV/#it{c}); #it{cos #theta_{p}}", kTH2F, {{8, 0.3, 4.3}, {1000, 0.9, 1.}});
    FullQaRegistry.add("FullV0QA/hInvMass", "; #it{m}_{#Lambda} (GeV/#it{c}^{2}); #it{m}_{#bar{#Lambda}} (GeV/#it{c}^{2})", kTH2F, {{150, 0.5, 2}, {150, 0.5, 2}});
  }

  /// Porduce QA plots for V0 selection in FemtoDream framework
  void process(o2::aod::FemtoDreamCollision& col, FemtoFullParticles& parts)
  {
    auto groupPartsOne = partsOne->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col.globalIndex());

    eventHisto.fillQA(col);

    for (auto& part : groupPartsOne) {

      FullQaRegistry.fill(HIST("FullV0QA/hPt"), part.pt());
      FullQaRegistry.fill(HIST("FullV0QA/hEta"), part.eta());
      FullQaRegistry.fill(HIST("FullV0QA/hPhi"), part.phi());
      FullQaRegistry.fill(HIST("FullV0QA/hDaughDCA"), part.daughDCA());
      FullQaRegistry.fill(HIST("FullV0QA/hTransRadius"), part.transRadius());
      FullQaRegistry.fill(HIST("FullV0QA/hDecayVtxX"), part.decayVtxX());
      FullQaRegistry.fill(HIST("FullV0QA/hDecayVtxY"), part.decayVtxY());
      FullQaRegistry.fill(HIST("FullV0QA/hDecayVtxZ"), part.decayVtxZ());
      FullQaRegistry.fill(HIST("FullV0QA/hCPA"), part.tempFitVar());
      FullQaRegistry.fill(HIST("FullV0QA/hCPAvsPt"), part.pt(), part.tempFitVar());
      FullQaRegistry.fill(HIST("FullV0QA/hInvMass"), part.mLambda(), part.mAntiLambda());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<femtoDreamDebugV0>(cfgc),
  };
  return workflow;
}
