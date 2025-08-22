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

#include "PWGCF/DataModel/FemtoDerived.h"
#include "PWGCF/FemtoDream/Core/femtoDreamEventHisto.h"
#include "PWGCF/FemtoDream/Core/femtoDreamParticleHisto.h"

#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"

#include "TVector3.h"

#include "fairlogger/Logger.h"

#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::analysis::femtoDream;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct femtoDreamDebugReso {
  SliceCache cache;

  struct : ConfigurableGroup {
    std::string prefix = std::string("resonance");

    Configurable<int> ConfReso_PDGCode{"ConfReso_PDGCode", 333, "Reso - PDG code"};
    Configurable<int> ConfReso_ChildPos_PDGCode{"ConfReso_PosChild_PDGCode", 321, "Positive Child - PDG code"};
    Configurable<int> ConfReso_ChildNeg_PDGCode{"ConfReso_NegChild_PDGCode", 321, "Negative Child- PDG code"};

    ConfigurableAxis ConfResoTempFitVarBins{"ConfResoTempFitVarBins", {300, 0.95, 1.}, "Reso: binning of the TempFitVar in the pT vs. TempFitVar plot"};
    ConfigurableAxis ConfResoTempFitVarMomentumBins{"ConfResoTempFitVarMomentumBins", {20, 0.5, 4.05}, "Reso: pT binning of the pT vs. TempFitVar plot"};
    ConfigurableAxis ConfBinmult{"ConfBinmult", {1, 0, 1}, "multiplicity Binning"};
    ConfigurableAxis ConfDummy{"ConfDummy", {1, 0, 1}, "Dummy axis for inv mass"};

    Configurable<int> ConfResoTempFitVarMomentum{"ConfResoTempFitVarMomentum", 0, "Momentum used for binning: 0 -> pt; 1 -> preco; 2 -> ptpc"};
    ConfigurableAxis ConfResoInvMassBins{"ConfResoInvMassBins", {200, 1, 1.2}, "Reso: InvMass binning"};

    ConfigurableAxis ConfResoChildTempFitVarMomentumBins{"ConfResoChildTempFitVarMomentumBins", {600, 0, 6}, "p binning for the p vs Nsigma TPC/TOF plot"};
    ConfigurableAxis ConfResoChildNsigmaTPCBins{"ConfResoChildNsigmaTPCBins", {1600, -8, 8}, "binning of Nsigma TPC plot"}; // TPC and TOf seperate doen't make sense really right??
    ConfigurableAxis ConfResoChildNsigmaTOFBins{"ConfResoChildNsigmaTOFBins", {3000, -15, 15}, "binning of the Nsigma TOF plot"};
    ConfigurableAxis ConfResoChildNsigmaTPCTOFBins{"ConfResoChildNsigmaTPCTOFBins", {1000, 0, 10}, "binning of the Nsigma TPC+TOF plot"};
    ConfigurableAxis ConfResoChildNsigmaITSBins{"ConfResoChildNsigmaITSBins", {600, -3, 3}, "binning of the Nsigma ITS plot"};

    Configurable<aod::femtodreamparticle::cutContainerType> ConfReso_ChildPos_CutBit{"ConfReso_ChildPos_CutBit", 4860458, "Positive Child of Reso - Selection bit from cutCulator"};
    Configurable<aod::femtodreamparticle::cutContainerType> ConfReso_ChildPos_TPCBit{"ConfReso_ChildPos_TPCBit", 16, "Positive Child of Reso - PID bit from cutCulator"};
    Configurable<aod::femtodreamparticle::cutContainerType> ConfReso_ChildPos_TPCTOFBit{"ConfReso_ChildPos_TPCTOFBit", 8, "Positive Child of Reso - PID bit from cutCulator"};
    Configurable<aod::femtodreamparticle::cutContainerType> ConfReso_ChildNeg_CutBit{"ConfReso_ChildNeg_CutBit", 4860457, "Negative Child of Reso - PID bit from cutCulator"};
    Configurable<aod::femtodreamparticle::cutContainerType> ConfReso_ChildNegMerged_TPCBit{"ConfReso_ChildNegMerged_TPCBit", 64, "Negative Child of Reso - PID bit from cutCulator"};       // change
    Configurable<aod::femtodreamparticle::cutContainerType> ConfReso_ChildNegMerged_TPCTOFBit{"ConfReso_ChildNegMerged_TPCTOFBit", 32, "Negative Child of Reso - PID bit from cutCulator"}; // change
    ConfigurableAxis ConfChildTempFitVarBins{"ConfChildTempFitVarBins", {300, -0.15, 0.15}, "V0 child: binning of the TempFitVar in the pT vs. TempFitVar plot"};
    ConfigurableAxis ConfChildTempFitVarpTBins{"ConfChildTempFitVarpTBins", {20, 0.5, 4.05}, "V0 child: pT binning of the pT vs. TempFitVar plot"};
  } resonance;

  using FemtoFullParticles = soa::Join<aod::FDParticles, aod::FDExtParticles>;

  Partition<FemtoFullParticles> partsTwo = (ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kPhiPosdaughTPC_NegdaughTPC), ncheckbit(aod::femtodreamparticle::pidcut, resonance.ConfReso_ChildPos_TPCBit) && ncheckbit(aod::femtodreamparticle::cut, resonance.ConfReso_ChildNegMerged_TPCBit), false) ||
                                            ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kPhiPosdaughTOF_NegdaughTOF), ncheckbit(aod::femtodreamparticle::pidcut, resonance.ConfReso_ChildPos_TPCTOFBit) && ncheckbit(aod::femtodreamparticle::cut, resonance.ConfReso_ChildNegMerged_TPCTOFBit), false) ||
                                            ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kPhiPosdaughTOF_NegdaughTPC), ncheckbit(aod::femtodreamparticle::pidcut, resonance.ConfReso_ChildPos_TPCTOFBit) && ncheckbit(aod::femtodreamparticle::cut, resonance.ConfReso_ChildNegMerged_TPCBit), false) ||
                                            ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kPhiPosdaughTPC_NegdaughTOF), ncheckbit(aod::femtodreamparticle::pidcut, resonance.ConfReso_ChildPos_TPCBit) && ncheckbit(aod::femtodreamparticle::cut, resonance.ConfReso_ChildNegMerged_TPCTOFBit), false));

  Preslice<FemtoFullParticles> perCol = aod::femtodreamparticle::fdCollisionId;

  /// Histogramming
  FemtoDreamEventHisto eventHisto;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kResoChild, 3> posResoChildHistos;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kResoChild, 4> negResoChildHistos;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kReso> ResoHistos;

  /// Histogram output
  HistogramRegistry EventRegistry{"Event", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry ResoRegistry{"FullResoQA", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    posResoChildHistos.init(&ResoRegistry, resonance.ConfBinmult, resonance.ConfDummy, resonance.ConfResoChildTempFitVarMomentumBins, resonance.ConfDummy, resonance.ConfDummy, resonance.ConfChildTempFitVarBins, resonance.ConfResoChildNsigmaTPCBins, resonance.ConfResoChildNsigmaTOFBins, resonance.ConfResoChildNsigmaTPCTOFBins, resonance.ConfResoChildNsigmaITSBins, resonance.ConfResoInvMassBins, resonance.ConfDummy, false, resonance.ConfReso_ChildPos_PDGCode.value, true); // isDebug == TRUE
    negResoChildHistos.init(&ResoRegistry, resonance.ConfBinmult, resonance.ConfDummy, resonance.ConfResoChildTempFitVarMomentumBins, resonance.ConfDummy, resonance.ConfDummy, resonance.ConfChildTempFitVarBins, resonance.ConfResoChildNsigmaTPCBins, resonance.ConfResoChildNsigmaTOFBins, resonance.ConfResoChildNsigmaTPCTOFBins, resonance.ConfResoChildNsigmaITSBins, resonance.ConfResoInvMassBins, resonance.ConfDummy, false, resonance.ConfReso_ChildNeg_PDGCode, true);       // isDebug == TRUE
    ResoHistos.init(&ResoRegistry, resonance.ConfBinmult, resonance.ConfDummy, resonance.ConfResoTempFitVarMomentumBins, resonance.ConfDummy, resonance.ConfDummy, resonance.ConfResoTempFitVarBins, resonance.ConfResoChildNsigmaTPCBins, resonance.ConfResoChildNsigmaTOFBins, resonance.ConfResoChildNsigmaTPCTOFBins, resonance.ConfResoChildNsigmaITSBins, resonance.ConfResoInvMassBins, resonance.ConfDummy, false, resonance.ConfReso_PDGCode.value, true);                        // isDebug == TRUE, isMc ==FALSE for all
    ResoRegistry.add("hArmenterosPodolanski/hArmenterosPodolanskiPlot", "; #alpha; p_{T} (MeV/#it{c})", kTH2F, {{100, -1, 1}, {500, -0.3, 2}});
  }

  /// Porduce QA plots for V0 & Reso selection in FemtoDream framework
  void process(o2::aod::FDCollision const& col, FemtoFullParticles const& parts)
  {

    auto groupPartsTwo = partsTwo->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache); // maybe . instead of -> ??
    for (auto& part : groupPartsTwo) {
      if (!part.has_children()) {
        LOG(warn) << " Reso has no children";
        continue;
      }

      const auto& posresoChild = parts.iteratorAt(part.index() - 2);
      const auto& negresoChild = parts.iteratorAt(part.index() - 1);
      if (posresoChild.globalIndex() != part.childrenIds()[0] || negresoChild.globalIndex() != part.childrenIds()[1]) {
        continue;
      }
      if (posresoChild.partType() == uint8_t(aod::femtodreamparticle::ParticleType::kResoChild) &&
          (posresoChild.cut() & resonance.ConfReso_ChildPos_CutBit) == resonance.ConfReso_ChildPos_CutBit &&
          negresoChild.partType() == uint8_t(aod::femtodreamparticle::ParticleType::kResoChild) &&
          (negresoChild.cut() & resonance.ConfReso_ChildNeg_CutBit) == resonance.ConfReso_ChildNeg_CutBit) {

        TVector3 p_parent(part.px(), part.py(), part.pz());                        // Parent momentum (px, py, pz)
        TVector3 p_plus(posresoChild.px(), posresoChild.py(), posresoChild.pz());  // Daughter 1 momentum (px, py, pz)
        TVector3 p_minus(negresoChild.px(), negresoChild.py(), negresoChild.pz()); // Daughter 2 momentum (px, py, pz)

        double pL_plus = p_plus.Dot(p_parent) / p_parent.Mag();
        double pL_minus = p_minus.Dot(p_parent) / p_parent.Mag();
        float alpha = (pL_plus - pL_minus) / (pL_plus + pL_minus);

        TVector3 p_perp = p_plus - (p_parent * (pL_plus / p_parent.Mag()));
        double qtarm = p_perp.Mag();

        ResoRegistry.fill(HIST("hArmenterosPodolanski/hArmenterosPodolanskiPlot"), alpha, qtarm);

        ResoHistos.fillQA<false, true>(part, static_cast<aod::femtodreamparticle::MomentumType>(resonance.ConfResoTempFitVarMomentum.value), col.multNtr(), col.multV0M());
        posResoChildHistos.fillQA<false, true>(posresoChild, static_cast<aod::femtodreamparticle::MomentumType>(resonance.ConfResoTempFitVarMomentum.value), col.multNtr(), col.multV0M());
        negResoChildHistos.fillQA<false, true>(negresoChild, static_cast<aod::femtodreamparticle::MomentumType>(resonance.ConfResoTempFitVarMomentum.value), col.multNtr(), col.multV0M());
      }
    }
  }
};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<femtoDreamDebugReso>(cfgc),
  };
  return workflow;
}
