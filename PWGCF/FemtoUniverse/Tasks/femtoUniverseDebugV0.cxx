// Copyright 2019-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file femtoUniverseDebugV0.cxx
/// \brief Tasks that reads the particle tables and fills QA histograms for V0s
/// \author Luca Barioglio, TU München, luca.barioglio@cern.ch
/// \author Zuzanna Chochulska, WUT Warsaw & CTU Prague, zchochul@cern.ch

#include <fairlogger/Logger.h>
#include <cstdint>
#include <vector>
#include <TDatabasePDG.h>
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "DataFormatsParameters/GRPObject.h"

#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseParticleHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseEventHisto.h"
#include "PWGCF/FemtoUniverse/Core/femtoUtils.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseMath.h"

using namespace o2;
using namespace o2::analysis::femto_universe;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct FemtoUniverseDebugV0 {

  Service<o2::framework::O2DatabasePDG> pdg;

  SliceCache cache;

  Configurable<int> confPDGCodeV0{"confPDGCodeV0", 3122, "V0 -- PDG code"};
  Configurable<int> confPDGCodePositiveChild{"confPDGCodePositiveChild", 2212, "Positive Child -- PDG code"};
  Configurable<int> confPDGCodeNegativeChild{"confPDGCodeNegativeChild", 211, "Negative Child -- PDG code"};
  Configurable<uint32_t> confCutV0{"confCutV0", 338, "V0 -- Selection bit from cutCulator"};
  ConfigurableAxis confV0TempFitVarBins{"confV0TempFitVarBins", {300, 0.95, 1.}, "V0: binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis confV0TempFitVarpTBins{"confV0TempFitVarpTBins", {20, 0.5, 4.05}, "V0: pT binning of the pT vs. TempFitVar plot"};

  Configurable<uint32_t> confCutPositiveChild{"confCutPositiveChild", 150, "Positive Child of V0 -- Selection bit from cutCulator"};
  Configurable<uint32_t> confCutNegativeChild{"confCutNegativeChild", 149, "Negative Child of V0 -- Selection bit from cutCulator"};
  Configurable<float> confPositiveChildPIDnSigmaMax{"confPositiveChildPIDnSigmaMax", 3.f, "Positive Child of V0 -- Selection bit from cutCulator"};
  Configurable<float> confNegativeChildPIDnSigmaMax{"confNegativeChildPIDnSigmaMax", 3.f, "Negative Child of V0 -- Selection bit from cutCulator"};
  Configurable<int> confPositiveChildIndex{"confPositiveChildIndex", 1, "Positive Child of V0 -- Index from cutCulator"};
  Configurable<int> confNegativeChildIndex{"confNegativeChildIndex", 0, "Negative Child of V0 -- Index from cutCulator"};
  Configurable<std::vector<float>> confChildPIDnSigmaMax{"confChildPIDnSigmaMax", std::vector<float>{4.f, 3.f}, "V0 child selection: max. PID nSigma TPC"};
  Configurable<int> confChildnSpecies{"confChildnSpecies", 2, "Number of particle spieces (for V0 children) with PID info"};
  ConfigurableAxis confChildTempFitVarBins{"confChildTempFitVarBins", {300, -0.15, 0.15}, "V0 child: binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis confChildTempFitVarpTBins{"confChildTempFitVarpTBins", {20, 0.5, 4.05}, "V0 child: pT binning of the pT vs. TempFitVar plot"};

  using FemtoFullParticles = soa::Join<aod::FDParticles, aod::FDExtParticles>;
  Partition<FemtoFullParticles> partsOne = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kV0)) && ((aod::femtouniverseparticle::cut & confCutV0) == confCutV0);
  Preslice<FemtoFullParticles> perCol = aod::femtouniverseparticle::fdCollisionId;

  /// Histogramming
  FemtoUniverseEventHisto eventHisto;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0Child, 3> posChildHistos;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0Child, 4> negChildHistos;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0> V0Histos;

  /// Histogram output
  HistogramRegistry EventRegistry{"Event", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry V0Registry{"FullV0QA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry thetaRegistry{"ThetaQA", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    eventHisto.init(&EventRegistry);
    posChildHistos.init(&V0Registry, confChildTempFitVarpTBins, confChildTempFitVarBins, false, confPDGCodePositiveChild.value, true);
    negChildHistos.init(&V0Registry, confChildTempFitVarpTBins, confChildTempFitVarBins, false, confPDGCodeNegativeChild, true);
    V0Histos.init(&V0Registry, confV0TempFitVarpTBins, confV0TempFitVarBins, false, confPDGCodeV0.value, true);

    thetaRegistry.add("Theta/hTheta", " ; p (GeV/#it{c}); cos(#theta)", kTH2F, {{100, 0, 10}, {50, -5, 5}});
  }

  /// Produce QA plots for V0 selection in FemtoUniverse framework
  void process(o2::aod::FdCollision const& col, FemtoFullParticles const& parts)
  {
    auto groupPartsOne = partsOne->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    eventHisto.fillQA(col);
    for (const auto& part : groupPartsOne) {
      if (!part.has_children()) {
        continue;
      }
      const auto& posChild = parts.iteratorAt(part.index() - 2);
      const auto& negChild = parts.iteratorAt(part.index() - 1);
      if (posChild.globalIndex() != part.childrenIds()[0] || negChild.globalIndex() != part.childrenIds()[1]) {
        LOG(warn) << "Indices of V0 children do not match";
        continue;
      }

      // Check cuts on V0 children
      if (posChild.partType() == uint8_t(aod::femtouniverseparticle::ParticleType::kV0Child) &&
          negChild.partType() == uint8_t(aod::femtouniverseparticle::ParticleType::kV0Child) &&
          isFullPIDSelected(posChild.pidCut(), posChild.p(), 999.f, confPositiveChildIndex.value, confChildnSpecies.value, confChildPIDnSigmaMax.value, confPositiveChildPIDnSigmaMax.value, 1.f) &&
          isFullPIDSelected(negChild.pidCut(), negChild.p(), 999.f, confNegativeChildIndex.value, confChildnSpecies.value, confChildPIDnSigmaMax.value, confNegativeChildPIDnSigmaMax.value, 1.f)) {
        auto pdgDB = TDatabasePDG::Instance();
        auto protonMass = pdgDB->GetParticle(confPDGCodePositiveChild)->Mass();
        auto pionMass = pdgDB->GetParticle(confPDGCodeNegativeChild)->Mass();
        auto protonBoosted = FemtoUniverseMath::boostPRF<decltype(posChild)>(posChild, protonMass, negChild, pionMass);
        auto cosineTheta = (protonBoosted.Px() * part.px() + protonBoosted.Py() * part.py() + protonBoosted.Pz() * part.pz()) / (protonBoosted.P() * part.p());

        V0Histos.fillQA<false, true>(part);
        posChildHistos.fillQA<false, true>(posChild);
        negChildHistos.fillQA<false, true>(negChild);
        thetaRegistry.fill(HIST("Theta/hTheta"), part.p(), cosineTheta);
      }
    }
  }
};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtoUniverseDebugV0>(cfgc),
  };
  return workflow;
}
