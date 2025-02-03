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
#include "fairlogger/Logger.h"

#include "PWGCF/DataModel/FemtoDerived.h"
#include "PWGCF/FemtoDream/Core/femtoDreamParticleHisto.h"
#include "PWGCF/FemtoDream/Core/femtoDreamEventHisto.h"

using namespace o2;
using namespace o2::analysis::femtoDream;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct femtoDreamDebugV0 {
  SliceCache cache;

  Configurable<int> ConfV01_PDGCode{"ConfV01_PDGCode", 3122, "V0 - PDG code"};
  Configurable<int> ConfV01_ChildPos_PDGCode{"ConfV01_PosChild_PDGCode", 2212, "Positive Child - PDG code"};
  Configurable<int> ConfV01_ChildNeg_PDGCode{"ConfV01_NegChild_PDGCode", 211, "Negative Child- PDG code"};
  Configurable<aod::femtodreamparticle::cutContainerType> ConfV01_CutBit{"ConfV01_CutBit", 338, "V0 - Selection bit from cutCulator"};
  ConfigurableAxis ConfV0TempFitVarBins{"ConfV0TempFitVarBins", {300, 0.95, 1.}, "V0: binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis ConfV0TempFitVarMomentumBins{"ConfV0TempFitVarMomentumBins", {20, 0.5, 4.05}, "V0: pT binning of the pT vs. TempFitVar plot"};
  ConfigurableAxis ConfBinmult{"ConfBinmult", {1, 0, 1}, "multiplicity Binning"};
  ConfigurableAxis ConfDummy{"ConfDummy", {1, 0, 1}, "Dummy axis for inv mass"};

  Configurable<int> ConfV0TempFitVarMomentum{"ConfV0TempFitVarMomentum", 0, "Momentum used for binning: 0 -> pt; 1 -> preco; 2 -> ptpc"};

  ConfigurableAxis ConfV0InvMassBins{"ConfV0InvMassBins", {200, 1, 1.2}, "V0: InvMass binning"};

  ConfigurableAxis ConfV0ChildTempFitVarMomentumBins{"ConfV0ChildTempFitVarMomentumBins", {600, 0, 6}, "p binning for the p vs Nsigma TPC/TOF plot"};
  ConfigurableAxis ConfV0ChildNsigmaTPCBins{"ConfV0ChildNsigmaTPCBins", {1600, -8, 8}, "binning of Nsigma TPC plot"};
  ConfigurableAxis ConfV0ChildNsigmaTOFBins{"ConfV0ChildNsigmaTOFBins", {3000, -15, 15}, "binning of the Nsigma TOF plot"};
  ConfigurableAxis ConfV0ChildNsigmaTPCTOFBins{"ConfV0ChildNsigmaTPCTOFBins", {1000, 0, 10}, "binning of the Nsigma TPC+TOF plot"};
  ConfigurableAxis ConfV0ChildNsigmaITSBins{"ConfV0ChildNsigmaITSBins", {600, -3, 3}, "binning of the Nsigma ITS plot"};

  Configurable<aod::femtodreamparticle::cutContainerType> ConfV01_ChildPos_CutBit{"ConfV01_ChildPos_CutBit", 150, "Positive Child of V0 - Selection bit from cutCulator"};
  Configurable<aod::femtodreamparticle::cutContainerType> ConfV01_ChildPos_TPCBit{"ConfV01_ChildPos_TPCBit", 4, "Positive Child of V0 - PID bit from cutCulator"};
  Configurable<aod::femtodreamparticle::cutContainerType> ConfV01_ChildNeg_CutBit{"ConfV01_ChildNeg_CutBit", 149, "Negative Child of V0 - PID bit from cutCulator"};
  Configurable<aod::femtodreamparticle::cutContainerType> ConfV01_ChildNeg_TPCBit{"ConfV01_ChildNeg_TPCBit", 8, "Negative Child of V0 - PID bit from cutCulator"};
  ConfigurableAxis ConfChildTempFitVarBins{"ConfChildTempFitVarBins", {300, -0.15, 0.15}, "V0 child: binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis ConfChildTempFitVarpTBins{"ConfChildTempFitVarpTBins", {20, 0.5, 4.05}, "V0 child: pT binning of the pT vs. TempFitVar plot"};

  using FemtoFullParticles = soa::Join<aod::FDParticles, aod::FDExtParticles>;
  Partition<FemtoFullParticles> partsOne = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kV0)) && (ncheckbit(aod::femtodreamparticle::cut, ConfV01_CutBit));
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
    eventHisto.init(&EventRegistry, false);
    posChildHistos.init(&V0Registry, ConfBinmult, ConfDummy, ConfV0ChildTempFitVarMomentumBins, ConfDummy, ConfDummy, ConfChildTempFitVarBins, ConfV0ChildNsigmaTPCBins, ConfV0ChildNsigmaTOFBins, ConfV0ChildNsigmaTPCTOFBins, ConfV0ChildNsigmaITSBins, ConfV0InvMassBins, false, ConfV01_ChildPos_PDGCode.value, true);
    negChildHistos.init(&V0Registry, ConfBinmult, ConfDummy, ConfV0ChildTempFitVarMomentumBins, ConfDummy, ConfDummy, ConfChildTempFitVarBins, ConfV0ChildNsigmaTPCBins, ConfV0ChildNsigmaTOFBins, ConfV0ChildNsigmaTPCTOFBins, ConfV0ChildNsigmaITSBins, ConfV0InvMassBins, false, ConfV01_ChildNeg_PDGCode, true);
    V0Histos.init(&V0Registry, ConfBinmult, ConfDummy, ConfV0TempFitVarMomentumBins, ConfDummy, ConfDummy, ConfV0TempFitVarBins, ConfV0ChildNsigmaTPCBins, ConfV0ChildNsigmaTOFBins, ConfV0ChildNsigmaTPCTOFBins, ConfV0ChildNsigmaITSBins, ConfV0InvMassBins, false, ConfV01_PDGCode.value, true);
  }

  /// Porduce QA plots for V0 selection in FemtoDream framework
  void process(o2::aod::FDCollision const& col, FemtoFullParticles const& parts)
  {
    auto groupPartsOne = partsOne->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    eventHisto.fillQA<false>(col);
    for (auto& part : groupPartsOne) {
      if (!part.has_children()) {
        continue;
      }
      // check cut on v0 children
      // TODO: check if this should be possible
      // auto posChild = part.template children_as<FemtoFullParticles>().front();
      // auto negChild = part.template children_as<FemtoFullParticles>().back();
      const auto& posChild = parts.iteratorAt(part.index() - 2);
      const auto& negChild = parts.iteratorAt(part.index() - 1);
      if (posChild.globalIndex() != part.childrenIds()[0] || negChild.globalIndex() != part.childrenIds()[1]) {
        LOG(warn) << "Indices of V0 children do not match";
        continue;
      }
      // check cuts on V0 children
      if (posChild.partType() == uint8_t(aod::femtodreamparticle::ParticleType::kV0Child) &&
          (posChild.cut() & ConfV01_ChildPos_CutBit) == ConfV01_ChildPos_CutBit &&
          (posChild.pidcut() & ConfV01_ChildPos_TPCBit) == ConfV01_ChildPos_TPCBit &&
          negChild.partType() == uint8_t(aod::femtodreamparticle::ParticleType::kV0Child) &&
          (negChild.cut() & ConfV01_ChildNeg_CutBit) == ConfV01_ChildNeg_CutBit &&
          (negChild.pidcut() & ConfV01_ChildNeg_TPCBit) == ConfV01_ChildNeg_TPCBit) {
        V0Histos.fillQA<false, true>(part, static_cast<aod::femtodreamparticle::MomentumType>(ConfV0TempFitVarMomentum.value), col.multNtr(), col.multV0M());
        posChildHistos.fillQA<false, true>(posChild, static_cast<aod::femtodreamparticle::MomentumType>(ConfV0TempFitVarMomentum.value), col.multNtr(), col.multV0M());
        negChildHistos.fillQA<false, true>(negChild, static_cast<aod::femtodreamparticle::MomentumType>(ConfV0TempFitVarMomentum.value), col.multNtr(), col.multV0M());
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
