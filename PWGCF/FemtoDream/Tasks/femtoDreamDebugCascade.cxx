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
/// \author Georgios Mantzaridis, TU München, luca.barioglio@cern.ch

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
#include "PWGCF/FemtoDream/Core/femtoDreamMath.h"
#include "PWGCF/FemtoDream/Core/femtoDreamUtils.h"

using namespace o2;
using namespace o2::analysis::femtoDream;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct femtoDreamDebugCascade {
  SliceCache cache;

  Configurable<int> ConfCascade_PDGCode{"ConfCascade_PDGCode", 3312, "Cascade - PDG code"};
  Configurable<int> ConfCascade_ChildPos_PDGCode{"ConfCascade_PosChild_PDGCode", 2212, "Positive Child - PDG code"};
  Configurable<int> ConfCascade_ChildNeg_PDGCode{"ConfCascade_NegChild_PDGCode", 211, "Negative Child- PDG code"};
  Configurable<int> ConfCascade_Bach_PDGCode{"ConfCascade_Bach_PDGCode", 211, "Bachelor Child- PDG code"};

  Configurable<aod::femtodreamparticle::cutContainerType> ConfCascade_CutBit{"ConfCascade_CutBit", 338, "Cascade - Selection bit from cutCulator"};
  ConfigurableAxis ConfCascadeTempFitVarBins{"ConfCascadeTempFitVarBins", {300, 0.95, 1.}, "Cascade: binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis ConfCascadeTempFitVarMomentumBins{"ConfCascadeTempFitVarMomentumBins", {20, 0.5, 4.05}, "Cascade: pT binning of the pT vs. TempFitVar plot"};
  ConfigurableAxis ConfBinmult{"ConfBinmult", {1, 0, 1}, "multiplicity Binning"};
  ConfigurableAxis ConfDummy{"ConfDummy", {1, 0, 1}, "Dummy axis for inv mass"};

  Configurable<int> ConfCascadeTempFitVarMomentum{"ConfCascadeTempFitVarMomentum", 0, "Momentum used for binning: 0 -> pt; 1 -> preco; 2 -> ptpc"};

  ConfigurableAxis ConfCascadeInvMassBins{"ConfCascadeInvMassBins", {200, 1.25, 1.45}, "Cascade: InvMass binning"};
  ConfigurableAxis ConfCascadeInvMassCompetingBins{"ConfCascadeInvMassCompetingBins", {200, 1.57, 1.77}, "Cascade: InvMass binning of the competing candidate"};

  ConfigurableAxis ConfCascadeChildTempFitVarMomentumBins{"ConfCascadeChildTempFitVarMomentumBins", {600, 0, 6}, "p binning for the p vs Nsigma TPC/TOF plot"};
  ConfigurableAxis ConfCascadeChildNsigmaTPCBins{"ConfCascadeChildNsigmaTPCBins", {1600, -8, 8}, "binning of Nsigma TPC plot"};
  ConfigurableAxis ConfCascadeChildNsigmaTOFBins{"ConfCascadeChildNsigmaTOFBins", {3000, -15, 15}, "binning of the Nsigma TOF plot"};
  ConfigurableAxis ConfCascadeChildNsigmaTPCTOFBins{"ConfCascadeChildNsigmaTPCTOFBins", {1000, 0, 10}, "binning of the Nsigma TPC+TOF plot"};

  Configurable<aod::femtodreamparticle::cutContainerType> ConfCascade_ChildPos_CutBit{"ConfCascade_ChildPos_CutBit", 150, "Positive Child of Cascade - Selection bit from cutCulator"};
  Configurable<aod::femtodreamparticle::cutContainerType> ConfCascade_ChildPos_TPCBit{"ConfCascade_ChildPos_TPCBit", 4, "Positive Child of Cascade - PID bit from cutCulator"};
  Configurable<aod::femtodreamparticle::cutContainerType> ConfCascade_ChildNeg_CutBit{"ConfCascade_ChildNeg_CutBit", 149, "Negative Child of Cascade - PID bit from cutCulator"};
  Configurable<aod::femtodreamparticle::cutContainerType> ConfCascade_ChildNeg_TPCBit{"ConfCascade_ChildNeg_TPCBit", 8, "Negative Child of Cascade - PID bit from cutCulator"};
  Configurable<aod::femtodreamparticle::cutContainerType> ConfCascade_ChildBach_CutBit{"ConfCascade_ChildBach_CutBit", 149, "Bachelor Child of Cascade - PID bit from cutCulator"};
  Configurable<aod::femtodreamparticle::cutContainerType> ConfCascade_ChildBach_TPCBit{"ConfCascade_ChildBach_TPCBit", 8, "Bachelor Child of Cascade - PID bit from cutCulator"};
  Configurable<bool> ConfUseChildCuts{"ConfUseChildCuts", true, "Use cuts on the children of the Cascades additional to those of the selection of the cascade builder"};
  Configurable<bool> ConfUseChildPIDCuts{"ConfUseChildPIDCuts", true, "Use cuts on the children of the Cascades additional to those of the selection of the cascade builder"};

  Configurable<bool> ConfIsOmega{"ConfIsOmega", false, "Switch between Xi and Omaga Cascades: If true: Omega; else: Xi"};
  Configurable<bool> ConfRejectCompetingMass{"ConfRejectCompetingMass", false, "Reject the competing Cascade Mass (use only for debugging. More efficient to exclude it already at the producer level)"};
  Configurable<float> ConfCompetingCascadeMassLowLimit{"ConfCompetingCascadeMassLowLimit", 0., "Lower Limit of the invariant mass window within which to reject the cascade"};
  Configurable<float> ConfCompetingCascadeMassUpLimit{"ConfCompetingCascadeMassUpLimit", 0., "Upper Limit of the invariant mass window within which to reject the cascade"};

  ConfigurableAxis ConfCascadeChildTempFitVarBins{"ConfCascadeChildTempFitVarBins", {300, -0.15, 0.15}, "Cascade child: binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis ConfCascadeChildTempFitVarpTBins{"ConfCascadeChildTempFitVarpTBins", {20, 0.5, 4.05}, "Cascade child: pT binning of the pT vs. TempFitVar plot"};

  using FemtoFullParticles = soa::Join<aod::FDParticles, aod::FDExtParticles>;
  Partition<FemtoFullParticles> partsOne = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kCascade)) && (ncheckbit(aod::femtodreamparticle::cut, ConfCascade_CutBit));
  Preslice<FemtoFullParticles> perCol = aod::femtodreamparticle::fdCollisionId;

  /// Histogramming
  FemtoDreamEventHisto eventHisto;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kCascadeV0Child, 3> posChildHistos;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kCascadeV0Child, 4> negChildHistos;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kCascadeBachelor, 8> bachelorHistos;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kCascade> CascadeHistos;

  /// Histogram output
  HistogramRegistry EventRegistry{"Event", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry CascadeRegistry{"FullCascQA", {}, OutputObjHandlingPolicy::AnalysisObject};

  float massProton;
  float massPion;
  float massKaon;
  float massLambda;
  float massCompetingBach;

  void init(InitContext&)
  {
    eventHisto.init(&EventRegistry, false);
    posChildHistos.init(&CascadeRegistry, ConfBinmult, ConfDummy, ConfCascadeChildTempFitVarMomentumBins, ConfDummy, ConfDummy, ConfCascadeChildTempFitVarBins, ConfCascadeChildNsigmaTPCBins, ConfCascadeChildNsigmaTOFBins, ConfCascadeChildNsigmaTPCTOFBins, ConfDummy, ConfCascadeInvMassBins, ConfCascadeInvMassCompetingBins, false, ConfCascade_ChildPos_PDGCode.value, true);
    negChildHistos.init(&CascadeRegistry, ConfBinmult, ConfDummy, ConfCascadeChildTempFitVarMomentumBins, ConfDummy, ConfDummy, ConfCascadeChildTempFitVarBins, ConfCascadeChildNsigmaTPCBins, ConfCascadeChildNsigmaTOFBins, ConfCascadeChildNsigmaTPCTOFBins, ConfDummy, ConfCascadeInvMassBins, ConfCascadeInvMassCompetingBins, false, ConfCascade_ChildNeg_PDGCode.value, true);
    bachelorHistos.init(&CascadeRegistry, ConfBinmult, ConfDummy, ConfCascadeChildTempFitVarMomentumBins, ConfDummy, ConfDummy, ConfCascadeChildTempFitVarBins, ConfCascadeChildNsigmaTPCBins, ConfCascadeChildNsigmaTOFBins, ConfCascadeChildNsigmaTPCTOFBins, ConfDummy, ConfCascadeInvMassBins, ConfCascadeInvMassCompetingBins, false, ConfCascade_Bach_PDGCode.value, true);
    CascadeHistos.init(&CascadeRegistry, ConfBinmult, ConfDummy, ConfCascadeTempFitVarMomentumBins, ConfDummy, ConfDummy, ConfCascadeTempFitVarBins, ConfCascadeChildNsigmaTPCBins, ConfCascadeChildNsigmaTOFBins, ConfCascadeChildNsigmaTPCTOFBins, ConfDummy, ConfCascadeInvMassBins, ConfCascadeInvMassCompetingBins, false, ConfCascade_PDGCode.value, true);

    massProton = o2::analysis::femtoDream::getMass(2212);
    massPion = o2::analysis::femtoDream::getMass(211);
    massKaon = o2::analysis::femtoDream::getMass(321);
    massLambda = o2::analysis::femtoDream::getMass(3122);
    if (ConfIsOmega) { // if the Cascade is an Omega, then the bachelor is a Kaon
      massCompetingBach = o2::analysis::femtoDream::getMass(211);
    } else { // if the Cascade is a Xi, then the bachelor is a Pion
      massCompetingBach = o2::analysis::femtoDream::getMass(321);
    }
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
      const auto& posChild = parts.iteratorAt(part.index() - 3);
      const auto& negChild = parts.iteratorAt(part.index() - 2);
      const auto& bachChild = parts.iteratorAt(part.index() - 1);
      if (posChild.globalIndex() != part.childrenIds()[0] || negChild.globalIndex() != part.childrenIds()[1] || bachChild.globalIndex() != part.childrenIds()[2]) {
        LOG(warn) << "Indices of V0 children do not match";
        continue;
      }
      // check cuts on V0 children
      if (posChild.partType() == uint8_t(aod::femtodreamparticle::ParticleType::kCascadeV0Child) &&
          negChild.partType() == uint8_t(aod::femtodreamparticle::ParticleType::kCascadeV0Child) &&
          bachChild.partType() == uint8_t(aod::femtodreamparticle::ParticleType::kCascadeBachelor)) {

        if (ConfUseChildCuts) {
          if (!((posChild.cut() & ConfCascade_ChildPos_CutBit) == ConfCascade_ChildPos_CutBit &&
                (negChild.cut() & ConfCascade_ChildNeg_CutBit) == ConfCascade_ChildNeg_CutBit &&
                (bachChild.cut() & ConfCascade_ChildBach_CutBit) == ConfCascade_ChildBach_CutBit)) {
            continue;
          }
        }
        if (ConfUseChildPIDCuts) {
          if (!((posChild.pidcut() & ConfCascade_ChildPos_TPCBit) == ConfCascade_ChildPos_TPCBit &&
                (negChild.pidcut() & ConfCascade_ChildNeg_TPCBit) == ConfCascade_ChildNeg_TPCBit &&
                (bachChild.pidcut() & ConfCascade_ChildBach_TPCBit) == ConfCascade_ChildBach_TPCBit)) {
            continue;
          }
        }

        // Competing mass rejection
        if (ConfRejectCompetingMass) {
          float invMassCompetingCasc;
          if (part.sign() < 0) {
            invMassCompetingCasc = FemtoDreamMath::getInvMassCascade(posChild, massProton, negChild, massPion, bachChild, massCompetingBach, massLambda);
          } else {
            invMassCompetingCasc = FemtoDreamMath::getInvMassCascade(posChild, massPion, negChild, massProton, bachChild, massCompetingBach, massLambda);
          }
          if (invMassCompetingCasc > ConfCompetingCascadeMassLowLimit.value &&
              invMassCompetingCasc < ConfCompetingCascadeMassUpLimit.value) {
            continue;
          }
        }
        CascadeHistos.fillQA<false, true>(part, static_cast<aod::femtodreamparticle::MomentumType>(ConfCascadeTempFitVarMomentum.value), col.multNtr(), col.multV0M());
        posChildHistos.fillQA<false, true>(posChild, static_cast<aod::femtodreamparticle::MomentumType>(ConfCascadeTempFitVarMomentum.value), col.multNtr(), col.multV0M());
        negChildHistos.fillQA<false, true>(negChild, static_cast<aod::femtodreamparticle::MomentumType>(ConfCascadeTempFitVarMomentum.value), col.multNtr(), col.multV0M());
        bachelorHistos.fillQA<false, true>(bachChild, static_cast<aod::femtodreamparticle::MomentumType>(ConfCascadeTempFitVarMomentum.value), col.multNtr(), col.multV0M());
      }
    }
  }
};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<femtoDreamDebugCascade>(cfgc),
  };
  return workflow;
}
