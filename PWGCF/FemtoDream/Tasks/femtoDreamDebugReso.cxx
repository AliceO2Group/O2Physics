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

/// \file femtoDreamDebugReso.cxx
/// \brief Tasks that reads the particle tables and fills QA histograms for V0s
/// \author Christopher Klumm, TU München, christopher.klumm@cern.ch

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
#include <string>
#include <vector>

using namespace o2;
using namespace o2::analysis::femtoDream;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct FemtoDreamDebugReso {
  SliceCache cache;

  struct : ConfigurableGroup {
    std::string prefix = std::string("resonance");

    Configurable<int> confResoPDGCode{"confResoPDGCode", 333, "Reso - PDG code"};
    Configurable<int> confResoChildPosPDGCode{"confResoChildPosPDGCode", 321, "Positive Child - PDG code"};
    Configurable<int> confResoChildNegPDGCode{"confResoChildNegPDGCode", 321, "Negative Child- PDG code"};

    ConfigurableAxis confResoTempFitVarBins{"confResoTempFitVarBins", {300, 0.95, 1.}, "Reso: binning of the TempFitVar in the pT vs. TempFitVar plot"};
    ConfigurableAxis confResoTempFitVarMomentumBins{"confResoTempFitVarMomentumBins", {20, 0.5, 4.05}, "Reso: pT binning of the pT vs. TempFitVar plot"};
    ConfigurableAxis confBinmult{"confBinmult", {1, 0, 1}, "multiplicity Binning"};
    ConfigurableAxis confDummy{"confDummy", {1, 0, 1}, "Dummy axis for inv mass"};

    Configurable<int> confResoTempFitVarMomentum{"confResoTempFitVarMomentum", 0, "Momentum used for binning: 0 -> pt; 1 -> preco; 2 -> ptpc"};
    ConfigurableAxis confResoInvMassBins{"confResoInvMassBins", {200, 1, 1.2}, "Reso: InvMass binning"};

    ConfigurableAxis confResoChildTempFitVarMomentumBins{"confResoChildTempFitVarMomentumBins", {600, 0, 6}, "p binning for the p vs Nsigma TPC/TOF plot"};
    ConfigurableAxis confResoChildNsigmaTPCBins{"confResoChildNsigmaTPCBins", {1600, -8, 8}, "binning of Nsigma TPC plot"}; // TPC and TOf seperate doen't make sense really right??
    ConfigurableAxis confResoChildNsigmaTOFBins{"confResoChildNsigmaTOFBins", {3000, -15, 15}, "binning of the Nsigma TOF plot"};
    ConfigurableAxis confResoChildNsigmaTPCTOFBins{"confResoChildNsigmaTPCTOFBins", {1000, 0, 10}, "binning of the Nsigma TPC+TOF plot"};
    ConfigurableAxis confResoChildNsigmaITSBins{"confResoChildNsigmaITSBins", {600, -3, 3}, "binning of the Nsigma ITS plot"};

    Configurable<aod::femtodreamparticle::cutContainerType> confResoChildPosCutBit{"confResoChildPosCutBit", 4860458, "Positive Child of Reso - Selection bit from cutCulator"};
    Configurable<aod::femtodreamparticle::cutContainerType> confResoChildPosTPCBit{"confResoChildPosTPCBit", 64, "Positive Child of Reso - PID bit from cutCulator"};
    Configurable<aod::femtodreamparticle::cutContainerType> confResoChildPosTPCTOFBit{"confResoChildPosTPCTOFBit", 32, "Positive Child of Reso - PID bit from cutCulator"};
    Configurable<aod::femtodreamparticle::cutContainerType> confResoChildNegCutBit{"confResoChildNegCutBit", 4860457, "Negative Child of Reso - PID bit from cutCulator"};
    Configurable<aod::femtodreamparticle::cutContainerType> confResoChildNegMergedTPCBit{"confResoChildNegMergedTPCBit", 258, "Negative Child of Reso - PID bit from cutCulator"};       // change
    Configurable<aod::femtodreamparticle::cutContainerType> confResoChildNegMergedTPCTOFBit{"confResoChildNegMergedTPCTOFBit", 130, "Negative Child of Reso - PID bit from cutCulator"}; // change
    ConfigurableAxis confChildTempFitVarBins{"confChildTempFitVarBins", {300, -0.15, 0.15}, "V0 child: binning of the TempFitVar in the pT vs. TempFitVar plot"};
    ConfigurableAxis confChildTempFitVarpTBins{"confChildTempFitVarpTBins", {20, 0.5, 4.05}, "V0 child: pT binning of the pT vs. TempFitVar plot"};
  } resonance;

  using FemtoFullParticles = soa::Join<aod::FDParticles, aod::FDExtParticles>;

  Partition<FemtoFullParticles> partsTwo = (ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kResoPosdaughTPC_NegdaughTPC), ncheckbit(aod::femtodreamparticle::pidcut, resonance.confResoChildPosTPCBit) && ncheckbit(aod::femtodreamparticle::cut, resonance.confResoChildNegMergedTPCBit), false) ||
                                            ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kResoPosdaughTOF_NegdaughTOF), ncheckbit(aod::femtodreamparticle::pidcut, resonance.confResoChildPosTPCTOFBit) && ncheckbit(aod::femtodreamparticle::cut, resonance.confResoChildNegMergedTPCTOFBit), false) ||
                                            ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kResoPosdaughTOF_NegdaughTPC), ncheckbit(aod::femtodreamparticle::pidcut, resonance.confResoChildPosTPCTOFBit) && ncheckbit(aod::femtodreamparticle::cut, resonance.confResoChildNegMergedTPCBit), false) ||
                                            ifnode(aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kResoPosdaughTPC_NegdaughTOF), ncheckbit(aod::femtodreamparticle::pidcut, resonance.confResoChildPosTPCBit) && ncheckbit(aod::femtodreamparticle::cut, resonance.confResoChildNegMergedTPCTOFBit), false));

  Preslice<FemtoFullParticles> perCol = aod::femtodreamparticle::fdCollisionId;

  /// Histogramming
  FemtoDreamEventHisto eventHisto;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kResoChild, 3> posResoChildHistos;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kResoChild, 4> negResoChildHistos;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kReso> resoHistos;

  /// Histogram output
  HistogramRegistry eventRegistry{"Event", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resoRegistry{"FullResoQA", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    posResoChildHistos.init(&resoRegistry, resonance.confBinmult, resonance.confDummy, resonance.confResoChildTempFitVarMomentumBins, resonance.confDummy, resonance.confDummy, resonance.confChildTempFitVarBins, resonance.confResoChildNsigmaTPCBins, resonance.confResoChildNsigmaTOFBins, resonance.confResoChildNsigmaTPCTOFBins, resonance.confResoChildNsigmaITSBins, resonance.confResoInvMassBins, resonance.confDummy, false, resonance.confResoChildPosPDGCode.value, true); // isDebug == TRUE
    negResoChildHistos.init(&resoRegistry, resonance.confBinmult, resonance.confDummy, resonance.confResoChildTempFitVarMomentumBins, resonance.confDummy, resonance.confDummy, resonance.confChildTempFitVarBins, resonance.confResoChildNsigmaTPCBins, resonance.confResoChildNsigmaTOFBins, resonance.confResoChildNsigmaTPCTOFBins, resonance.confResoChildNsigmaITSBins, resonance.confResoInvMassBins, resonance.confDummy, false, resonance.confResoChildNegPDGCode, true);       // isDebug == TRUE
    resoHistos.init(&resoRegistry, resonance.confBinmult, resonance.confDummy, resonance.confResoTempFitVarMomentumBins, resonance.confDummy, resonance.confDummy, resonance.confResoTempFitVarBins, resonance.confResoChildNsigmaTPCBins, resonance.confResoChildNsigmaTOFBins, resonance.confResoChildNsigmaTPCTOFBins, resonance.confResoChildNsigmaITSBins, resonance.confResoInvMassBins, resonance.confDummy, false, resonance.confResoPDGCode.value, true);                       // isDebug == TRUE, isMc ==FALSE for all
    resoRegistry.add("hArmenterosPodolanski/hArmenterosPodolanskiPlot", "; #alpha; p_{T} (MeV/#it{c})", kTH2F, {{100, -1, 1}, {500, -0.3, 2}});
  }

  /// Porduce QA plots for V0 & Reso selection in FemtoDream framework
  void process(o2::aod::FDCollision const& col, FemtoFullParticles const& parts)
  {

    auto groupPartsTwo = partsTwo->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache); // maybe . instead of -> ??
    for (const auto& part : groupPartsTwo) {
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
          (posresoChild.cut() & resonance.confResoChildPosCutBit) == resonance.confResoChildPosCutBit &&
          negresoChild.partType() == uint8_t(aod::femtodreamparticle::ParticleType::kResoChild) &&
          (negresoChild.cut() & resonance.confResoChildNegCutBit) == resonance.confResoChildNegCutBit) {

        TVector3 pparent(part.px(), part.py(), part.pz());                        // Parent momentum (px, py, pz)
        TVector3 pplus(posresoChild.px(), posresoChild.py(), posresoChild.pz());  // Daughter 1 momentum (px, py, pz)
        TVector3 pminus(negresoChild.px(), negresoChild.py(), negresoChild.pz()); // Daughter 2 momentum (px, py, pz)

        double pLplus = pplus.Dot(pparent) / pparent.Mag();
        double pLminus = pminus.Dot(pparent) / pparent.Mag();
        float alpha = (pLplus - pLminus) / (pLplus + pLminus);

        TVector3 pperp = pplus - (pparent * (pLplus / pparent.Mag()));
        double qtarm = pperp.Mag();

        resoRegistry.fill(HIST("hArmenterosPodolanski/hArmenterosPodolanskiPlot"), alpha, qtarm);

        resoHistos.fillQA<false, true>(part, static_cast<aod::femtodreamparticle::MomentumType>(resonance.confResoTempFitVarMomentum.value), col.multNtr(), col.multV0M());
        posResoChildHistos.fillQA<false, true>(posresoChild, static_cast<aod::femtodreamparticle::MomentumType>(resonance.confResoTempFitVarMomentum.value), col.multNtr(), col.multV0M());
        negResoChildHistos.fillQA<false, true>(negresoChild, static_cast<aod::femtodreamparticle::MomentumType>(resonance.confResoTempFitVarMomentum.value), col.multNtr(), col.multV0M());
      }
    }
  }
};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtoDreamDebugReso>(cfgc),
  };
  return workflow;
}
