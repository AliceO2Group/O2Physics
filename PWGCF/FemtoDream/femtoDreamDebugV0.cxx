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

/// \file femtoDreamDebugV0.cxx
/// \brief Tasks that reads the particle tables and fills QA histograms for V0s
/// \author Luca Barioglio, TU München, luca.barioglio@cern.ch

#include <fairlogger/Logger.h>
#include <cstdint>
#include <iostream>
#include <vector>
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
  SliceCache cache;

  Configurable<int> ConfPDGCodeV0{"ConfPDGCodePartOne", 3122, "V0 - PDG code"};
  Configurable<int> ConfPDGCodeChildPos{"ConfPDGCodeChildPos", 2212, "Positive Child - PDG code"};
  Configurable<int> ConfPDGCodeChildNeg{"ConfPDGCodeChildNeg", 211, "Negative Child- PDG code"};
  Configurable<uint32_t> ConfCutV0{"ConfCutV0", 338, "V0 - Selection bit from cutCulator"};
  ConfigurableAxis ConfV0TempFitVarBins{"ConfV0TempFitVarBins", {300, 0.95, 1.}, "V0: binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis ConfV0TempFitVarMomentumBins{"ConfV0TempFitVarMomentumBins", {20, 0.5, 4.05}, "V0: pT binning of the pT vs. TempFitVar plot"};

  Configurable<int> ConfTempFitVarMomentum{"ConfTempFitVarMomentum", 0, "Momentum used for binning: 0 -> pt; 1 -> preco; 2 -> ptpc"};

  ConfigurableAxis ConfV0InvMassBins{"ConfV0InvMassBins", {200, 1, 1.2}, "V0: InvMass binning"};

  ConfigurableAxis ConfChildTempFitVarMomentumBins{"ConfChildTempFitVarMomentumBins", {600, 0, 6}, "p binning for the p vs Nsigma TPC/TOF plot"};
  ConfigurableAxis ConfChildNsigmaTPCBins{"ConfChildNsigmaTPCBins", {1600, -8, 8}, "binning of Nsigma TPC plot"};
  ConfigurableAxis ConfChildNsigmaTOFBins{"ConfChildNsigmaTOFBins", {3000, -15, 15}, "binning of the Nsigma TOF plot"};
  ConfigurableAxis ConfChildNsigmaTPCTOFBins{"ConfChildNsigmaTPCTOFBins", {1000, 0, 10}, "binning of the Nsigma TPC+TOF plot"};

  Configurable<uint32_t> ConfCutChildPos{"ConfCutChildPos", 150, "Positive Child of V0 - Selection bit from cutCulator"};
  Configurable<uint32_t> ConfCutChildNeg{"ConfCutChildNeg", 149, "Negative Child of V0 - Selection bit from cutCulator"};
  Configurable<float> ConfChildPosPidnSigmaMax{"ConfChildPosPidnSigmaMax", 3.f, "Positive Child of V0 - Selection bit from cutCulator"};
  Configurable<float> ConfChildNegPidnSigmaMax{"ConfChildNegPidnSigmaMax", 3.f, "Negative Child of V0 - Selection bit from cutCulator"};
  Configurable<int> ConfChildPosIndex{"ConfChildPosIndex", 1, "Positive Child of V0 - Index from cutCulator"};
  Configurable<int> ConfChildNegIndex{"ConfChildNegIndex", 0, "Negative Child of V0 - Index from cutCulator"};
  Configurable<std::vector<float>> ConfChildPIDnSigmaMax{"ConfChildPIDnSigmaMax", std::vector<float>{4.f, 3.f}, "V0 child sel: Max. PID nSigma TPC"};
  Configurable<int> ConfChildnSpecies{"ConfChildnSpecies", 2, "Number of particle spieces (for V0 children) with PID info"};
  ConfigurableAxis ConfChildTempFitVarBins{"ConfChildTempFitVarBins", {300, -0.15, 0.15}, "V0 child: binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis ConfChildTempFitVarpTBins{"ConfChildTempFitVarpTBins", {20, 0.5, 4.05}, "V0 child: pT binning of the pT vs. TempFitVar plot"};

  using FemtoFullParticles = soa::Join<aod::FDParticles, aod::FDExtParticles>;
  Partition<FemtoFullParticles> partsOne = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kV0)) && ((aod::femtodreamparticle::cut & ConfCutV0) == ConfCutV0);
  Preslice<FemtoFullParticles> perCol = aod::femtodreamparticle::fdCollisionId;

  /// Histogramming
  FemtoDreamEventHisto eventHisto;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kV0Child, 3> posChildHistos;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kV0Child, 4> negChildHistos;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kV0> V0Histos;

  /// Histogram output
  HistogramRegistry EventRegistry{"Event", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry V0Registry{"FullV0QA", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    eventHisto.init(&EventRegistry);
    posChildHistos.init(&V0Registry, ConfChildTempFitVarMomentumBins, ConfChildTempFitVarBins, ConfChildNsigmaTPCBins, ConfChildNsigmaTOFBins, ConfChildNsigmaTPCTOFBins, ConfV0InvMassBins, false, ConfPDGCodeChildPos.value, true);
    negChildHistos.init(&V0Registry, ConfChildTempFitVarMomentumBins, ConfChildTempFitVarBins, ConfChildNsigmaTPCBins, ConfChildNsigmaTOFBins, ConfChildNsigmaTPCTOFBins, ConfV0InvMassBins, false, ConfPDGCodeChildNeg, true);
    V0Histos.init(&V0Registry, ConfV0TempFitVarMomentumBins, ConfV0TempFitVarBins, ConfChildNsigmaTPCBins, ConfChildNsigmaTOFBins, ConfChildNsigmaTPCTOFBins, ConfV0InvMassBins, false, ConfPDGCodeV0.value, true);
  }

  /// Porduce QA plots for V0 selection in FemtoDream framework
  void process(o2::aod::FDCollision const& col, FemtoFullParticles const& parts)
  {
    auto groupPartsOne = partsOne->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    eventHisto.fillQA(col);
    for (auto& part : groupPartsOne) {
      if (!part.has_children()) {
        continue;
      }
      const auto& posChild = parts.iteratorAt(part.index() - 2);
      const auto& negChild = parts.iteratorAt(part.index() - 1);
      if (posChild.globalIndex() != part.childrenIds()[0] || negChild.globalIndex() != part.childrenIds()[1]) {
        LOG(warn) << "Indices of V0 children do not match";
        continue;
      }
      // check cuts on V0 children
      if ((posChild.partType() == uint8_t(aod::femtodreamparticle::ParticleType::kV0Child) && (posChild.cut() & ConfCutChildPos) == ConfCutChildPos) &&
          (negChild.partType() == uint8_t(aod::femtodreamparticle::ParticleType::kV0Child) && (negChild.cut() & ConfCutChildNeg) == ConfCutChildNeg) &&
          isFullPIDSelected(posChild.pidcut(), posChild.p(), 999.f, ConfChildPosIndex.value, ConfChildnSpecies.value, ConfChildPIDnSigmaMax.value, ConfChildPosPidnSigmaMax.value, 1.f) && isFullPIDSelected(negChild.pidcut(), negChild.p(), 999.f, ConfChildNegIndex.value, ConfChildnSpecies.value, ConfChildPIDnSigmaMax.value, ConfChildNegPidnSigmaMax.value, 1.f)) {
        V0Histos.fillQA<false, true>(part, static_cast<aod::femtodreamparticle::MomentumType>(ConfTempFitVarMomentum.value));
        posChildHistos.fillQA<false, true>(posChild, static_cast<aod::femtodreamparticle::MomentumType>(ConfTempFitVarMomentum.value));
        negChildHistos.fillQA<false, true>(negChild, static_cast<aod::femtodreamparticle::MomentumType>(ConfTempFitVarMomentum.value));
      }
    }
  }
};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<femtoDreamDebugV0>(cfgc),
  };
  return workflow;
}
