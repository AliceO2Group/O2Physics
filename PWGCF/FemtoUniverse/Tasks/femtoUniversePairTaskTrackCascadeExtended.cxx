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

/// \brief Task for cascade QA; in the future: for cascade correlations
/// \author Barbara Chytla, WUT Warsaw, barbara.chytla@cern.ch

#include <vector>
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseParticleHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseEventHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniversePairCleaner.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseContainer.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseDetaDphiStar.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUtils.h"

using namespace o2;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femtoUniverse;
using namespace o2::aod::pidutils;

struct femtoUniversePairTaskTrackCascadeExtended {

  SliceCache cache;
  using FemtoFullParticles = soa::Join<aod::FDParticles, aod::FDExtParticles>;
  Preslice<FemtoFullParticles> perCol = aod::femtouniverseparticle::fdCollisionId;

  Configurable<float> ConfZVertexCut{"ConfZVertexCut", 10.f, "Event sel: Maximum z-Vertex (cm)"};

  Filter collisionFilter = (nabs(aod::collision::posZ) < ConfZVertexCut);
  using FilteredFDCollisions = soa::Filtered<o2::aod::FDCollisions>;
  using FilteredFDCollision = FilteredFDCollisions::iterator;

  ConfigurableAxis ConfChildTempFitVarpTBins{"ConfChildTempFitVarpTBins", {20, 0.5, 4.05}, "V0 child: pT binning of the pT vs. TempFitVar plot"};
  ConfigurableAxis ConfChildTempFitVarBins{"ConfChildTempFitVarBins", {300, -0.15, 0.15}, "V0 child: binning of the TempFitVar in the pT vs. TempFitVar plot"};
  Configurable<float> ConfCascInvMassLowLimit{"ConfCascInvMassLowLimit", 1.315, "Lower limit of the Casc invariant mass"};
  Configurable<float> ConfCascInvMassUpLimit{"ConfCascInvMassUpLimit", 1.325, "Upper limit of the Casc invariant mass"};
  Configurable<float> ConfCascTranRad{"ConfCascTranRad", 0.5, "Cascade transverse radius"};

  Configurable<float> NSigmaTPCPion{"NSigmaTPCPion", 4, "NSigmaTPCPion"};
  Configurable<float> NSigmaTPCProton{"NSigmaTPCProton", 4, "NSigmaTPCProton"};

  /// Partition for cascades
  Partition<FemtoFullParticles> cascs = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kCascade));

  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0Child, 3> posChildHistos;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0Child, 4> negChildHistos;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kCascadeBachelor, 2> bachHistos;

  HistogramRegistry rXiQA{"xi", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};

  bool invMCascade(float invMassCascade, float invMassAntiCascade)
  {
    if ((invMassCascade < ConfCascInvMassLowLimit || invMassCascade > ConfCascInvMassUpLimit) && (invMassAntiCascade < ConfCascInvMassLowLimit || invMassAntiCascade > ConfCascInvMassUpLimit)) {
      return false;
    }
    return true;
  }

  void init(InitContext const&)
  {
    // Axes
    AxisSpec XiMassAxis = {200, 1.28f, 1.36f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec ptAxis = {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec etaAxis = {100, -2.0f, 2.0f, "#it{#eta}"};
    AxisSpec phiAxis = {100, 0.0f, 6.0f, "#it{#phi}"};

    // Histograms
    rXiQA.add("hMassXiMinus", "hMassXiMinus", {HistType::kTH1F, {XiMassAxis}});
    rXiQA.add("hMassXiPlus", "hMassXiPlus", {HistType::kTH1F, {XiMassAxis}});
    rXiQA.add("hMassXiMinusSelected", "hMassXiSelected", {HistType::kTH1F, {XiMassAxis}});
    rXiQA.add("hMassXiPlusSelected", "hMassXiSelected", {HistType::kTH1F, {XiMassAxis}});
    rXiQA.add("hPtXi", "hPtXi", {HistType::kTH1F, {{ptAxis}}});
    rXiQA.add("hEtaXi", "hEtaXi", {HistType::kTH1F, {{etaAxis}}});
    rXiQA.add("hPhiXi", "hPhiXi", {HistType::kTH1F, {{phiAxis}}});
    rXiQA.add("hCascCosPA", "hCascCosPA", {HistType::kTH1F, {{100, 0.9f, 1.f}}});
    // rXiQA.add("hCascDCAV0Daughters", "hCascDCAV0Daughters", {HistType::kTH1F, {{55, 0.0f, 2.2f}}});

    posChildHistos.init(&qaRegistry, ConfChildTempFitVarpTBins, ConfChildTempFitVarBins, false, 0, true);
    negChildHistos.init(&qaRegistry, ConfChildTempFitVarpTBins, ConfChildTempFitVarBins, false, 0, true);
    bachHistos.init(&qaRegistry, ConfChildTempFitVarpTBins, ConfChildTempFitVarBins, false, 0, true, "hBachelor");
  }

  void processCascades(FilteredFDCollision& col, FemtoFullParticles& parts)
  {
    auto groupCascs = cascs->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    // const int multCol = col.multNtr();

    for (auto& casc : groupCascs) {
      rXiQA.fill(HIST("hMassXiMinus"), casc.mLambda());
      rXiQA.fill(HIST("hMassXiPlus"), casc.mAntiLambda());

      // if (!invMCascade(casc.mLambda(), casc.mAntiLambda()))
      //   continue;

      const auto& posChild = parts.iteratorAt(casc.index() - 3);
      const auto& negChild = parts.iteratorAt(casc.index() - 2);
      const auto& bachelor = parts.iteratorAt(casc.index() - 1);

      // if (casc.transRadius() < ConfCascTranRad)
      //   continue;
      // std::cout<<std::endl;
      // std::cout<<"TYPE:"<<std::endl;
      // std::cout<<casc.partType()<<std::endl;
      //  nSigma selection for daughter and bachelor tracks
      if (casc.sign() < 0) {
        if (TMath::Abs(posChild.tpcNSigmaPr()) > NSigmaTPCProton) {
          continue;
        }
        if (TMath::Abs(negChild.tpcNSigmaPi()) > NSigmaTPCPion) {
          continue;
        }
      } else {
        if (TMath::Abs(negChild.tpcNSigmaPr()) > NSigmaTPCProton) {
          continue;
        }
        if (TMath::Abs(posChild.tpcNSigmaPi()) > NSigmaTPCPion) {
          continue;
        }
      }
      if (TMath::Abs(bachelor.tpcNSigmaPi()) > NSigmaTPCPion) {
        continue;
      }

      rXiQA.fill(HIST("hPtXi"), casc.pt());
      rXiQA.fill(HIST("hEtaXi"), casc.eta());
      rXiQA.fill(HIST("hPhiXi"), casc.phi());
      rXiQA.fill(HIST("hMassXiMinusSelected"), casc.mLambda());
      rXiQA.fill(HIST("hMassXiPlusSelected"), casc.mAntiLambda());
      rXiQA.fill(HIST("hCascCosPA"), casc.tempFitVar());
      // rXiQA.fill(HIST("hCascDCAV0Daughters"), casc.dcaV0daughters()); // nie ma miejsca na to w FemtoDerived

      posChildHistos.fillQA<false, true>(posChild);
      negChildHistos.fillQA<false, true>(negChild);
      bachHistos.fillQABase<false, true>(bachelor, HIST("hBachelor"));
    }
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackCascadeExtended, processCascades, "Enable processing cascades", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<femtoUniversePairTaskTrackCascadeExtended>(cfgc),
  };
  return workflow;
}
