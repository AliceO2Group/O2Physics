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

/// \file femtoWorldDebugV0.cxx
/// \brief Tasks that reads the particle tables and fills QA histograms for V0s
/// \author Luca Barioglio, TU München, luca.barioglio@cern.ch

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "DataFormatsParameters/GRPObject.h"

#include "PWGCF/FemtoWorld/DataModel/FemtoWorldDerived.h"
#include "PWGCF/FemtoWorld/Core/FemtoWorldParticleHisto.h"
#include "PWGCF/FemtoWorld/Core/FemtoWorldEventHisto.h"
#include "PWGCF/FemtoWorld/Core/FemtoWorldUtils.h"

using namespace o2;
using namespace o2::analysis::femtoWorld;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct femtoWorldDebugV0 {
  SliceCache cache;
  Preslice<aod::FemtoWorldParticles> perCol = aod::femtoworldparticle::femtoWorldCollisionId;

  Configurable<int> ConfPDGCodePartOne{"ConfPDGCodePartOne", 3122, "Particle 1 - PDG code"};
  Configurable<uint32_t> ConfCutPartOne{"ConfCutPartOne", 338, "Particle 1 - Selection bit from cutCulator"};

  using FemtoFullParticles = soa::Join<aod::FemtoWorldParticles, aod::FemtoWorldDebugParticles>;

  Partition<FemtoFullParticles> partsOne = (aod::femtoworldparticle::partType == uint8_t(aod::femtoworldparticle::ParticleType::kV0)) &&
                                           // (aod::femtoworldparticle::pt < cfgCutTable->get("MaxPt")) &&
                                           ((aod::femtoworldparticle::cut & ConfCutPartOne) == ConfCutPartOne);

  /// Histogramming for Event
  FemtoWorldEventHisto eventHisto;

  /// The configurables need to be passed to an std::vector
  std::vector<int> vPIDPartOne;

  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry FullQaRegistry{"FullV0QA", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    eventHisto.init(&qaRegistry);

    AxisSpec massAxisLambda = {600, 0.0f, 3.0f, "m_{#Lambda} (GeV/#it{c}^{2})"};
    AxisSpec massAxisAntiLambda = {600, 0.0f, 3.0f, "m_{#bar{#Lambda}} (GeV/#it{c}^{2})"};

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
    FullQaRegistry.add("FullV0QA/hInvMassLambda", "", kTH1F, {massAxisLambda});
    FullQaRegistry.add("FullV0QA/hInvMassAntiLambda", "", kTH1F, {massAxisAntiLambda});
    FullQaRegistry.add("FullV0QA/hInvMassLambdaAntiLambda", "", kTH2F, {massAxisLambda, massAxisAntiLambda});
  }

  /// Porduce QA plots for V0 selection in FemtoWorld framework
  void process(o2::aod::FemtoWorldCollision& col, FemtoFullParticles& parts)
  {
    auto groupPartsOne = partsOne->sliceByCached(aod::femtoworldparticle::femtoWorldCollisionId, col.globalIndex(), cache);

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
      FullQaRegistry.fill(HIST("FullV0QA/hInvMassLambda"), part.mLambda());
      FullQaRegistry.fill(HIST("FullV0QA/hInvMassAntiLambda"), part.mAntiLambda());
      FullQaRegistry.fill(HIST("FullV0QA/hInvMassLambdaAntiLambda"), part.mLambda(), part.mAntiLambda());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<femtoWorldDebugV0>(cfgc),
  };
  return workflow;
}
