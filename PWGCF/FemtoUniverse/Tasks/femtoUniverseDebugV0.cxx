// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
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
/// \author Anna-Mariia Andrushko, WUT Warsaw, anna-mariia.andrushko@cern.ch

#include "PWGCF/FemtoUniverse/Core/FemtoUniverseEventHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseMath.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseParticleHisto.h"
#include "PWGCF/FemtoUniverse/Core/femtoUtils.h"
#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"

#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/PID.h"

#include <fairlogger/Logger.h>

#include <cstdint>
#include <vector>

using namespace o2;
using namespace o2::analysis::femto_universe;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct FemtoUniverseDebugV0 {

  Service<o2::framework::O2DatabasePDG> pdg;

  SliceCache cache;

  /// V0 configurables
  struct : o2::framework::ConfigurableGroup {
    Configurable<int> confPDGCodeV0{"confPDGCodeV0", 3122, "V0 -- PDG code"};
    Configurable<int> confPDGCodePositiveChild{"confPDGCodePositiveChild", 2212, "Positive Child -- PDG code"};
    Configurable<int> confPDGCodeNegativeChild{"confPDGCodeNegativeChild", 211, "Negative Child -- PDG code"};
    Configurable<uint32_t> confCutV0{"confCutV0", 338, "V0 -- Selection bit from cutCulator"};
    ConfigurableAxis confV0TempFitVarBins{"confV0TempFitVarBins", {300, 0.95, 1.}, "V0: binning of the TempFitVar in the pT vs. TempFitVar plot"};
    ConfigurableAxis confV0TempFitVarpTBins{"confV0TempFitVarpTBins", {20, 0.5, 4.05}, "V0: pT binning of the pT vs. TempFitVar plot"};
  } V0configs; // o2-linter: disable=name/function-variable

  /// Children configurables
  struct : o2::framework::ConfigurableGroup {
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
  } childconfigs;

  Configurable<bool> confIsMC{"confIsMC", false, "Enable additional histograms in the case of a Monte Carlo run"};

  /// Partitioning
  using FemtoFullParticles = soa::Join<aod::FDParticles, aod::FDExtParticles>;

  Partition<FemtoFullParticles> partsOne = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kV0)) && ((aod::femtouniverseparticle::cut & V0configs.confCutV0) == V0configs.confCutV0);

  Partition<soa::Join<FemtoFullParticles, aod::FDMCLabels>> partsOneMC = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kV0)) && ((aod::femtouniverseparticle::cut & V0configs.confCutV0) == V0configs.confCutV0);

  Preslice<FemtoFullParticles> perCol = aod::femtouniverseparticle::fdCollisionId;

  /// Histogramming
  FemtoUniverseEventHisto eventHisto;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0Child, 3> positiveChildHistos;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0Child, 4> negativeChildHistos;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0> V0Histos; // o2-linter: disable=name/function-variable

  /// Histogram output
  HistogramRegistry eventRegistry{"Event", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry V0Registry{"FullV0QA", {}, OutputObjHandlingPolicy::AnalysisObject}; // o2-linter: disable=name/function-variable
  HistogramRegistry thetaRegistry{"ThetaQA", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    eventHisto.init(&eventRegistry);
    positiveChildHistos.init(&V0Registry, childconfigs.confChildTempFitVarpTBins, childconfigs.confChildTempFitVarBins, confIsMC, V0configs.confPDGCodePositiveChild.value, true);
    negativeChildHistos.init(&V0Registry, childconfigs.confChildTempFitVarpTBins, childconfigs.confChildTempFitVarBins, confIsMC, V0configs.confPDGCodeNegativeChild.value, true);
    V0Histos.init(&V0Registry, V0configs.confV0TempFitVarpTBins, V0configs.confV0TempFitVarBins, confIsMC, V0configs.confPDGCodeV0.value, true);

    thetaRegistry.add("Theta/hTheta", " ; p (GeV/#it{c}); cos(#theta)", kTH2F, {{100, 0, 10}, {110, -1.1, 1.1}});
    thetaRegistry.add("Theta/PositiveChild/hThetaPt", " ; p_{T} (GeV/#it{c}); cos(#theta)", kTH2F, {{100, 0, 10}, {110, -1.1, 1.1}});
    thetaRegistry.add("Theta/PositiveChild/hThetaEta", " ; #eta; cos(#theta)", kTH2F, {{100, -1, 1}, {110, -1.1, 1.1}});
    thetaRegistry.add("Theta/PositiveChild/hThetaPhi", " ; #phi; cos(#theta)", kTH2F, {{100, -1, 7}, {110, -1.1, 1.1}});
    thetaRegistry.add("Theta/NegativeChild/hThetaPt", " ; p_{T} (GeV/#it{c}); cos(#theta)", kTH2F, {{100, 0, 10}, {110, -1.1, 1.1}});
    thetaRegistry.add("Theta/NegativeChild/hThetaEta", " ; #eta; cos(#theta)", kTH2F, {{100, -1, 1}, {110, -1.1, 1.1}});
    thetaRegistry.add("Theta/NegativeChild/hThetaPhi", " ; #phi; cos(#theta)", kTH2F, {{100, -1, 7}, {110, -1.1, 1.1}});
  }

  /// Produce QA plots for V0 and its children on real data
  void processData(o2::aod::FdCollision const& col, FemtoFullParticles const& parts)
  {
    auto groupPartsOne = partsOne->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    eventHisto.fillQA(col);
    for (const auto& part : groupPartsOne) {
      if (!part.has_children()) {
        continue;
      }
      const auto& positiveChild = parts.iteratorAt(part.index() - 2);
      const auto& negativeChild = parts.iteratorAt(part.index() - 1);
      if (positiveChild.globalIndex() != part.childrenIds()[0] || negativeChild.globalIndex() != part.childrenIds()[1]) {
        LOG(warn) << "Indices of V0 children do not match";
        continue;
      }

      // Check cuts on V0 children
      if (positiveChild.partType() == uint8_t(aod::femtouniverseparticle::ParticleType::kV0Child) &&
          negativeChild.partType() == uint8_t(aod::femtouniverseparticle::ParticleType::kV0Child) &&
          isFullPIDSelected(positiveChild.pidCut(), positiveChild.p(), 999.f, childconfigs.confPositiveChildIndex.value, childconfigs.confChildnSpecies.value, childconfigs.confChildPIDnSigmaMax.value, childconfigs.confPositiveChildPIDnSigmaMax.value, 1.f) &&
          isFullPIDSelected(negativeChild.pidCut(), negativeChild.p(), 999.f, childconfigs.confNegativeChildIndex.value, childconfigs.confChildnSpecies.value, childconfigs.confChildPIDnSigmaMax.value, childconfigs.confNegativeChildPIDnSigmaMax.value, 1.f)) {
        auto positiveChildMass = pdg->Mass(V0configs.confPDGCodePositiveChild);
        auto negativeChildMass = pdg->Mass(V0configs.confPDGCodeNegativeChild);
        auto positiveChildBoosted = FemtoUniverseMath::boostPRF<decltype(positiveChild)>(positiveChild, positiveChildMass, negativeChild, negativeChildMass);
        auto cosineTheta = (positiveChildBoosted.Px() * part.px() + positiveChildBoosted.Py() * part.py() + positiveChildBoosted.Pz() * part.pz()) / (positiveChildBoosted.P() * part.p());

        V0Histos.fillQA<false, true>(part);
        positiveChildHistos.fillQA<false, true>(positiveChild);
        negativeChildHistos.fillQA<false, true>(negativeChild);

        thetaRegistry.fill(HIST("Theta/hTheta"), part.p(), cosineTheta);
        thetaRegistry.fill(HIST("Theta/PositiveChild/hThetaPt"), positiveChild.pt(), cosineTheta);
        thetaRegistry.fill(HIST("Theta/PositiveChild/hThetaEta"), positiveChild.eta(), cosineTheta);
        thetaRegistry.fill(HIST("Theta/PositiveChild/hThetaPhi"), positiveChild.phi(), cosineTheta);
        thetaRegistry.fill(HIST("Theta/NegativeChild/hThetaPt"), negativeChild.pt(), cosineTheta);
        thetaRegistry.fill(HIST("Theta/NegativeChild/hThetaEta"), negativeChild.eta(), cosineTheta);
        thetaRegistry.fill(HIST("Theta/NegativeChild/hThetaPhi"), negativeChild.phi(), cosineTheta);
      }
    }
  }
  PROCESS_SWITCH(FemtoUniverseDebugV0, processData, "Enable processing on real data", true);

  /// Produce QA plots for V0 and its children on Monte Carlo
  void processMC(o2::aod::FdCollision const& col, soa::Join<FemtoFullParticles, aod::FDMCLabels> const& parts, o2::aod::FdMCParticles const&)
  {
    auto groupPartsOne = partsOneMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    eventHisto.fillQA(col);
    for (const auto& part : groupPartsOne) {
      if (!part.has_children()) {
        continue;
      }
      const auto& positiveChild = parts.iteratorAt(part.index() - 2);
      const auto& negativeChild = parts.iteratorAt(part.index() - 1);
      if (positiveChild.globalIndex() != part.childrenIds()[0] || negativeChild.globalIndex() != part.childrenIds()[1]) {
        LOG(warn) << "Indices of V0 children do not match";
        continue;
      }

      // Check cuts on V0 children
      if (positiveChild.partType() == uint8_t(aod::femtouniverseparticle::ParticleType::kV0Child) &&
          negativeChild.partType() == uint8_t(aod::femtouniverseparticle::ParticleType::kV0Child) &&
          isFullPIDSelected(positiveChild.pidCut(), positiveChild.p(), 999.f, childconfigs.confPositiveChildIndex.value, childconfigs.confChildnSpecies.value, childconfigs.confChildPIDnSigmaMax.value, childconfigs.confPositiveChildPIDnSigmaMax.value, 1.f) &&
          isFullPIDSelected(negativeChild.pidCut(), negativeChild.p(), 999.f, childconfigs.confNegativeChildIndex.value, childconfigs.confChildnSpecies.value, childconfigs.confChildPIDnSigmaMax.value, childconfigs.confNegativeChildPIDnSigmaMax.value, 1.f)) {
        auto positiveChildMass = pdg->Mass(V0configs.confPDGCodePositiveChild);
        auto negativeChildMass = pdg->Mass(V0configs.confPDGCodeNegativeChild);
        auto positiveChildBoosted = FemtoUniverseMath::boostPRF<decltype(positiveChild)>(positiveChild, positiveChildMass, negativeChild, negativeChildMass);
        auto cosineTheta = (positiveChildBoosted.Px() * part.px() + positiveChildBoosted.Py() * part.py() + positiveChildBoosted.Pz() * part.pz()) / (positiveChildBoosted.P() * part.p());

        V0Histos.fillQA<false, true>(part);
        positiveChildHistos.fillQA<false, true>(positiveChild);
        negativeChildHistos.fillQA<false, true>(negativeChild);

        thetaRegistry.fill(HIST("Theta/hTheta"), part.p(), cosineTheta);
        thetaRegistry.fill(HIST("Theta/PositiveChild/hThetaPt"), positiveChild.pt(), cosineTheta);
        thetaRegistry.fill(HIST("Theta/PositiveChild/hThetaEta"), positiveChild.eta(), cosineTheta);
        thetaRegistry.fill(HIST("Theta/PositiveChild/hThetaPhi"), positiveChild.phi(), cosineTheta);
        thetaRegistry.fill(HIST("Theta/NegativeChild/hThetaPt"), negativeChild.pt(), cosineTheta);
        thetaRegistry.fill(HIST("Theta/NegativeChild/hThetaEta"), negativeChild.eta(), cosineTheta);
        thetaRegistry.fill(HIST("Theta/NegativeChild/hThetaPhi"), negativeChild.phi(), cosineTheta);
      }
    }
  }
  PROCESS_SWITCH(FemtoUniverseDebugV0, processMC, "Enable processing on Monte Carlo", false);
};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtoUniverseDebugV0>(cfgc),
  };
  return workflow;
}
