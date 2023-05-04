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

  Configurable<int> ConfPDGCodePartOne{"ConfPDGCodePartOne", 3122, "Particle 1 - PDG code"};
  Configurable<uint32_t> ConfCutV0{"ConfCutV0", 338, "V0 - Selection bit from cutCulator"};

  Configurable<uint32_t> ConfCutChildPos{"ConfCutChildPos", 150, "Positive Child of V0 - Selection bit from cutCulator"};
  Configurable<uint32_t> ConfCutChildNeg{"ConfCutChildNeg", 149, "Negative Child of V0 - Selection bit from cutCulator"};

  Configurable<float> ConfChildPosPidnSigmaMax{"ConfChildPosPidnSigmaMax", 3.f, "Positive Child of V0 - Selection bit from cutCulator"};
  Configurable<float> ConfChildNegPidnSigmaMax{"ConfChildNegPidnSigmaMax", 3.f, "Negative Child of V0 - Selection bit from cutCulator"};
  Configurable<int> ConfChildPosPos{"ConfChildPosPos", 1, "Positive Child of V0 - Selection bit from cutCulator"};
  Configurable<int> ConfChildNegPos{"ConfChildNegPos", 0, "Negative Child of V0 - Selection bit from cutCulator"};

  Configurable<std::vector<float>> ConfChildPIDnSigmaMax{
    "ConfChildPIDnSigmaMax", std::vector<float>{4.f, 3.f},
    "V0 Child sel: Max. PID nSigma TPC"};
  Configurable<int> cfgNspecies{"ccfgNspecies", 2, "Number of particle spieces (for V0 children) with PID info"};

  using FemtoFullParticles = soa::Join<aod::FemtoDreamParticles, aod::FemtoDreamDebugParticles>;

  Partition<FemtoFullParticles> partsOne = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kV0)) && ((aod::femtodreamparticle::cut & ConfCutV0) == ConfCutV0);

  Preslice<FemtoFullParticles> perCol = aod::femtodreamparticle::femtoDreamCollisionId;

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

    FullQaRegistry.add("FullDaughterPosQA/hCharge", "; Q (e); Entries",
                       kTH1F, {{5, -2.5, 2.5}});
    FullQaRegistry.add("FullDaughterPosQA/hPt", "; #it{p}_{T} (GeV/#it{c}); Entries",
                       kTH1F, {{240, 0, 6}});
    FullQaRegistry.add("FullDaughterPosQA/hEta", "; #eta; Entries", kTH1F,
                       {{200, -1.5, 1.5}});
    FullQaRegistry.add("FullDaughterPosQA/hPhi", "; #phi; Entries", kTH1F,
                       {{200, 0, 2. * M_PI}});
    FullQaRegistry.add("FullDaughterPosQA/hTPCfindable",
                       "; TPC findable clusters; Entries", kTH1F,
                       {{163, -0.5, 162.5}});
    FullQaRegistry.add("FullDaughterPosQA/hTPCfound", "; TPC found clusters; Entries",
                       kTH1F, {{163, -0.5, 162.5}});
    FullQaRegistry.add("FullDaughterPosQA/hTPCcrossedOverFindalbe",
                       "; TPC ratio findable; Entries", kTH1F,
                       {{100, 0.5, 1.5}});
    FullQaRegistry.add("FullDaughterPosQA/hTPCcrossedRows",
                       "; TPC crossed rows; Entries", kTH1F,
                       {{163, -0.5, 162.5}});
    FullQaRegistry.add("FullDaughterPosQA/hTPCfindableVsCrossed",
                       ";TPC findable clusters ; TPC crossed rows;", kTH2F,
                       {{163, -0.5, 162.5}, {163, -0.5, 162.5}});
    FullQaRegistry.add("FullDaughterPosQA/hTPCshared",
                       "; TPC shared clusters; Entries", kTH1F,
                       {{163, -0.5, 162.5}});
    FullQaRegistry.add("FullDaughterPosQA/hITSclusters", "; ITS clusters; Entries",
                       kTH1F, {{10, -0.5, 9.5}});
    FullQaRegistry.add("FullDaughterPosQA/hITSclustersIB",
                       "; ITS clusters in IB; Entries", kTH1F,
                       {{10, -0.5, 9.5}});
    FullQaRegistry.add("FullDaughterPosQA/hDCAxy",
                       "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F,
                       {{20, 0.5, 4.05}, {500, -5, 5}});
    FullQaRegistry.add("FullDaughterPosQA/hDCAz",
                       "; #it{p}_{T} (GeV/#it{c}); DCA_{z} (cm)", kTH2F,
                       {{100, 0, 10}, {500, -5, 5}});
    FullQaRegistry.add("FullDaughterPosQA/hDCA",
                       "; #it{p}_{T} (GeV/#it{c}); DCA (cm)", kTH2F,
                       {{100, 0, 10}, {301, 0., 1.5}});
    FullQaRegistry.add("FullDaughterPosQA/hTPCdEdX",
                       "; #it{p} (GeV/#it{c}); TPC Signal", kTH2F,
                       {{100, 0, 10}, {1000, 0, 1000}});
    FullQaRegistry.add("FullDaughterPosQA/nSigmaTPC_el",
                       "; #it{p} (GeV/#it{c}); n#sigma_{TPC}^{e}", kTH2F,
                       {{100, 0, 10}, {100, -5, 5}});
    FullQaRegistry.add("FullDaughterPosQA/nSigmaTPC_pi",
                       "; #it{p} (GeV/#it{c}); n#sigma_{TPC}^{#pi}", kTH2F,
                       {{100, 0, 10}, {100, -5, 5}});
    FullQaRegistry.add("FullDaughterPosQA/nSigmaTPC_K",
                       "; #it{p} (GeV/#it{c}); n#sigma_{TPC}^{K}", kTH2F,
                       {{100, 0, 10}, {100, -5, 5}});
    FullQaRegistry.add("FullDaughterPosQA/nSigmaTPC_p",
                       "; #it{p} (GeV/#it{c}); n#sigma_{TPC}^{p}", kTH2F,
                       {{100, 0, 10}, {100, -5, 5}});
    FullQaRegistry.add("FullDaughterPosQA/nSigmaTPC_d",
                       "; #it{p} (GeV/#it{c}); n#sigma_{TPC}^{d}", kTH2F,
                       {{100, 0, 10}, {100, -5, 5}});
    FullQaRegistry.add("FullDaughterPosQA/nSigmaTOF_el",
                       "; #it{p} (GeV/#it{c}); n#sigma_{TOF}^{e}", kTH2F,
                       {{100, 0, 10}, {100, -5, 5}});
    FullQaRegistry.add("FullDaughterPosQA/nSigmaTOF_pi",
                       "; #it{p} (GeV/#it{c}); n#sigma_{TOF}^{#pi}", kTH2F,
                       {{100, 0, 10}, {100, -5, 5}});
    FullQaRegistry.add("FullDaughterPosQA/nSigmaTOF_K",
                       "; #it{p} (GeV/#it{c}); n#sigma_{TOF}^{K}", kTH2F,
                       {{100, 0, 10}, {100, -5, 5}});
    FullQaRegistry.add("FullDaughterPosQA/nSigmaTOF_p",
                       "; #it{p} (GeV/#it{c}); n#sigma_{TOF}^{p}", kTH2F,
                       {{100, 0, 10}, {100, -5, 5}});
    FullQaRegistry.add("FullDaughterPosQA/nSigmaTOF_d",
                       "; #it{p} (GeV/#it{c}); n#sigma_{TOF}^{d}", kTH2F,
                       {{100, 0, 10}, {100, -5, 5}});
    FullQaRegistry.add("FullDaughterPosQA/nSigmaComb_el",
                       "; #it{p} (GeV/#it{c}); n#sigma_{comb}^{e}", kTH2F,
                       {{100, 0, 10}, {100, -5, 5}});
    FullQaRegistry.add("FullDaughterPosQA/nSigmaComb_pi",
                       "; #it{p} (GeV/#it{c}); n#sigma_{comb}^{#pi}", kTH2F,
                       {{100, 0, 10}, {100, -5, 5}});
    FullQaRegistry.add("FullDaughterPosQA/nSigmaComb_K",
                       "; #it{p} (GeV/#it{c}); n#sigma_{comb}^{K}", kTH2F,
                       {{100, 0, 10}, {100, -5, 5}});
    FullQaRegistry.add("FullDaughterPosQA/nSigmaComb_p",
                       "; #it{p} (GeV/#it{c}); n#sigma_{comb}^{p}", kTH2F,
                       {{100, 0, 10}, {100, -5, 5}});
    FullQaRegistry.add("FullDaughterPosQA/nSigmaComb_d",
                       "; #it{p} (GeV/#it{c}); n#sigma_{comb}^{d}", kTH2F,
                       {{100, 0, 10}, {100, -5, 5}});

    FullQaRegistry.add("FullDaughterNegQA/hCharge", "; Q (e); Entries",
                       kTH1F, {{5, -2.5, 2.5}});
    FullQaRegistry.add("FullDaughterNegQA/hPt", "; #it{p}_{T} (GeV/#it{c}); Entries",
                       kTH1F, {{240, 0, 6}});
    FullQaRegistry.add("FullDaughterNegQA/hEta", "; #eta; Entries", kTH1F,
                       {{200, -1.5, 1.5}});
    FullQaRegistry.add("FullDaughterNegQA/hPhi", "; #phi; Entries", kTH1F,
                       {{200, 0, 2. * M_PI}});
    FullQaRegistry.add("FullDaughterNegQA/hTPCfindable",
                       "; TPC findable clusters; Entries", kTH1F,
                       {{163, -0.5, 162.5}});
    FullQaRegistry.add("FullDaughterNegQA/hTPCfound", "; TPC found clusters; Entries",
                       kTH1F, {{163, -0.5, 162.5}});
    FullQaRegistry.add("FullDaughterNegQA/hTPCcrossedOverFindalbe",
                       "; TPC ratio findable; Entries", kTH1F,
                       {{100, 0.5, 1.5}});
    FullQaRegistry.add("FullDaughterNegQA/hTPCcrossedRows",
                       "; TPC crossed rows; Entries", kTH1F,
                       {{163, -0.5, 162.5}});
    FullQaRegistry.add("FullDaughterNegQA/hTPCfindableVsCrossed",
                       ";TPC findable clusters ; TPC crossed rows;", kTH2F,
                       {{163, -0.5, 162.5}, {163, -0.5, 162.5}});
    FullQaRegistry.add("FullDaughterNegQA/hTPCshared",
                       "; TPC shared clusters; Entries", kTH1F,
                       {{163, -0.5, 162.5}});
    FullQaRegistry.add("FullDaughterNegQA/hITSclusters", "; ITS clusters; Entries",
                       kTH1F, {{10, -0.5, 9.5}});
    FullQaRegistry.add("FullDaughterNegQA/hITSclustersIB",
                       "; ITS clusters in IB; Entries", kTH1F,
                       {{10, -0.5, 9.5}});
    FullQaRegistry.add("FullDaughterNegQA/hDCAxy",
                       "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F,
                       {{20, 0.5, 4.05}, {500, -5, 5}});
    FullQaRegistry.add("FullDaughterNegQA/hDCAz",
                       "; #it{p}_{T} (GeV/#it{c}); DCA_{z} (cm)", kTH2F,
                       {{100, 0, 10}, {500, -5, 5}});
    FullQaRegistry.add("FullDaughterNegQA/hDCA",
                       "; #it{p}_{T} (GeV/#it{c}); DCA (cm)", kTH2F,
                       {{100, 0, 10}, {301, 0., 1.5}});
    FullQaRegistry.add("FullDaughterNegQA/hTPCdEdX",
                       "; #it{p} (GeV/#it{c}); TPC Signal", kTH2F,
                       {{100, 0, 10}, {1000, 0, 1000}});
    FullQaRegistry.add("FullDaughterNegQA/nSigmaTPC_el",
                       "; #it{p} (GeV/#it{c}); n#sigma_{TPC}^{e}", kTH2F,
                       {{100, 0, 10}, {100, -5, 5}});
    FullQaRegistry.add("FullDaughterNegQA/nSigmaTPC_pi",
                       "; #it{p} (GeV/#it{c}); n#sigma_{TPC}^{#pi}", kTH2F,
                       {{100, 0, 10}, {100, -5, 5}});
    FullQaRegistry.add("FullDaughterNegQA/nSigmaTPC_K",
                       "; #it{p} (GeV/#it{c}); n#sigma_{TPC}^{K}", kTH2F,
                       {{100, 0, 10}, {100, -5, 5}});
    FullQaRegistry.add("FullDaughterNegQA/nSigmaTPC_p",
                       "; #it{p} (GeV/#it{c}); n#sigma_{TPC}^{p}", kTH2F,
                       {{100, 0, 10}, {100, -5, 5}});
    FullQaRegistry.add("FullDaughterNegQA/nSigmaTPC_d",
                       "; #it{p} (GeV/#it{c}); n#sigma_{TPC}^{d}", kTH2F,
                       {{100, 0, 10}, {100, -5, 5}});
    FullQaRegistry.add("FullDaughterNegQA/nSigmaTOF_el",
                       "; #it{p} (GeV/#it{c}); n#sigma_{TOF}^{e}", kTH2F,
                       {{100, 0, 10}, {100, -5, 5}});
    FullQaRegistry.add("FullDaughterNegQA/nSigmaTOF_pi",
                       "; #it{p} (GeV/#it{c}); n#sigma_{TOF}^{#pi}", kTH2F,
                       {{100, 0, 10}, {100, -5, 5}});
    FullQaRegistry.add("FullDaughterNegQA/nSigmaTOF_K",
                       "; #it{p} (GeV/#it{c}); n#sigma_{TOF}^{K}", kTH2F,
                       {{100, 0, 10}, {100, -5, 5}});
    FullQaRegistry.add("FullDaughterNegQA/nSigmaTOF_p",
                       "; #it{p} (GeV/#it{c}); n#sigma_{TOF}^{p}", kTH2F,
                       {{100, 0, 10}, {100, -5, 5}});
    FullQaRegistry.add("FullDaughterNegQA/nSigmaTOF_d",
                       "; #it{p} (GeV/#it{c}); n#sigma_{TOF}^{d}", kTH2F,
                       {{100, 0, 10}, {100, -5, 5}});
    FullQaRegistry.add("FullDaughterNegQA/nSigmaComb_el",
                       "; #it{p} (GeV/#it{c}); n#sigma_{comb}^{e}", kTH2F,
                       {{100, 0, 10}, {100, -5, 5}});
    FullQaRegistry.add("FullDaughterNegQA/nSigmaComb_pi",
                       "; #it{p} (GeV/#it{c}); n#sigma_{comb}^{#pi}", kTH2F,
                       {{100, 0, 10}, {100, -5, 5}});
    FullQaRegistry.add("FullDaughterNegQA/nSigmaComb_K",
                       "; #it{p} (GeV/#it{c}); n#sigma_{comb}^{K}", kTH2F,
                       {{100, 0, 10}, {100, -5, 5}});
    FullQaRegistry.add("FullDaughterNegQA/nSigmaComb_p",
                       "; #it{p} (GeV/#it{c}); n#sigma_{comb}^{p}", kTH2F,
                       {{100, 0, 10}, {100, -5, 5}});
    FullQaRegistry.add("FullDaughterNegQA/nSigmaComb_d",
                       "; #it{p} (GeV/#it{c}); n#sigma_{comb}^{d}", kTH2F,
                       {{100, 0, 10}, {100, -5, 5}});
  }

  /// Porduce QA plots for V0 selection in FemtoDream framework
  void process(o2::aod::FemtoDreamCollision const& col, FemtoFullParticles const& parts)
  {
    auto groupPartsOne = partsOne->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col.globalIndex(), cache);

    eventHisto.fillQA(col);

    for (auto& part : groupPartsOne) {

      if (!part.has_children()) {
        continue;
      }

      const auto& posChild = parts.iteratorAt(part.index() - 2);
      const auto& negChild = parts.iteratorAt(part.index() - 1);

      if (posChild.globalIndex() != part.childrenIds()[0] || negChild.globalIndex() != part.childrenIds()[1]) {
        LOG(warn) << "Indices of V0 children do not match";
        LOG(info) << posChild.globalIndex();
        LOG(info) << part.childrenIds()[0];
        continue;
      }

      // check cuts on V0 children
      if ((posChild.partType() == uint8_t(aod::femtodreamparticle::ParticleType::kV0Child) && (posChild.cut() & ConfCutChildPos) == ConfCutChildPos) &&
          (negChild.partType() == uint8_t(aod::femtodreamparticle::ParticleType::kV0Child) && (negChild.cut() & ConfCutChildNeg) == ConfCutChildNeg)) {

        if (!isPIDSelected(posChild.pidcut(),
                           std::vector<int>(ConfChildPosPos.value),
                           cfgNspecies.value,
                           ConfChildPosPidnSigmaMax.value,
                           ConfChildPIDnSigmaMax.value,
                           o2::analysis::femtoDream::kDetector::kTPC) ||
            !isPIDSelected(negChild.pidcut(),
                           std::vector<int>(ConfChildNegPos.value),
                           cfgNspecies.value,
                           ConfChildNegPidnSigmaMax.value,
                           ConfChildPIDnSigmaMax.value,
                           o2::analysis::femtoDream::kDetector::kTPC)) {
          continue;
        }

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

        // asdf
        FullQaRegistry.fill(HIST("FullDaughterPosQA/hCharge"), posChild.sign());
        FullQaRegistry.fill(HIST("FullDaughterPosQA/hPt"), posChild.pt());
        FullQaRegistry.fill(HIST("FullDaughterPosQA/hEta"), posChild.eta());
        FullQaRegistry.fill(HIST("FullDaughterPosQA/hPhi"), posChild.phi());
        FullQaRegistry.fill(HIST("FullDaughterPosQA/hDCAxy"), posChild.pt(),
                            posChild.tempFitVar());
        FullQaRegistry.fill(HIST("FullDaughterPosQA/hTPCfindable"),
                            posChild.tpcNClsFindable());
        FullQaRegistry.fill(HIST("FullDaughterPosQA/hTPCfound"), posChild.tpcNClsFound());
        FullQaRegistry.fill(HIST("FullDaughterPosQA/hTPCcrossedOverFindalbe"),
                            posChild.tpcCrossedRowsOverFindableCls());
        FullQaRegistry.fill(HIST("FullDaughterPosQA/hTPCcrossedRows"),
                            posChild.tpcNClsCrossedRows());
        FullQaRegistry.fill(HIST("FullDaughterPosQA/hTPCfindableVsCrossed"),
                            posChild.tpcNClsFindable(), posChild.tpcNClsCrossedRows());
        FullQaRegistry.fill(HIST("FullDaughterPosQA/hTPCshared"), posChild.tpcNClsShared());
        FullQaRegistry.fill(HIST("FullDaughterPosQA/hITSclusters"), posChild.itsNCls());
        FullQaRegistry.fill(HIST("FullDaughterPosQA/hITSclustersIB"),
                            posChild.itsNClsInnerBarrel());
        FullQaRegistry.fill(HIST("FullDaughterPosQA/hDCAz"), posChild.pt(), posChild.dcaZ());
        FullQaRegistry.fill(
          HIST("FullDaughterPosQA/hDCA"), posChild.pt(),
          std::sqrt(pow(posChild.dcaXY(), 2.) + pow(posChild.dcaZ(), 2.)));
        FullQaRegistry.fill(HIST("FullDaughterPosQA/hTPCdEdX"), posChild.p(),
                            posChild.tpcSignal());
        FullQaRegistry.fill(HIST("FullDaughterPosQA/nSigmaTPC_el"), posChild.p(),
                            posChild.tpcNSigmaEl());
        FullQaRegistry.fill(HIST("FullDaughterPosQA/nSigmaTPC_pi"), posChild.p(),
                            posChild.tpcNSigmaPi());
        FullQaRegistry.fill(HIST("FullDaughterPosQA/nSigmaTPC_K"), posChild.p(),
                            posChild.tpcNSigmaKa());
        FullQaRegistry.fill(HIST("FullDaughterPosQA/nSigmaTPC_p"), posChild.p(),
                            posChild.tpcNSigmaPr());
        FullQaRegistry.fill(HIST("FullDaughterPosQA/nSigmaTPC_d"), posChild.p(),
                            posChild.tpcNSigmaDe());
        FullQaRegistry.fill(HIST("FullDaughterPosQA/nSigmaTOF_el"), posChild.p(),
                            posChild.tofNSigmaEl());
        FullQaRegistry.fill(HIST("FullDaughterPosQA/nSigmaTOF_pi"), posChild.p(),
                            posChild.tofNSigmaPi());
        FullQaRegistry.fill(HIST("FullDaughterPosQA/nSigmaTOF_K"), posChild.p(),
                            posChild.tofNSigmaKa());
        FullQaRegistry.fill(HIST("FullDaughterPosQA/nSigmaTOF_p"), posChild.p(),
                            posChild.tofNSigmaPr());
        FullQaRegistry.fill(HIST("FullDaughterPosQA/nSigmaTOF_d"), posChild.p(),
                            posChild.tofNSigmaDe());
        FullQaRegistry.fill(HIST("FullDaughterPosQA/nSigmaComb_el"), posChild.p(),
                            std::sqrt(posChild.tpcNSigmaEl() * posChild.tpcNSigmaEl() +
                                      posChild.tofNSigmaEl() * posChild.tofNSigmaEl()));
        FullQaRegistry.fill(HIST("FullDaughterPosQA/nSigmaComb_pi"), posChild.p(),
                            std::sqrt(posChild.tpcNSigmaPi() * posChild.tpcNSigmaPi() +
                                      posChild.tofNSigmaPi() * posChild.tofNSigmaPi()));
        FullQaRegistry.fill(HIST("FullDaughterPosQA/nSigmaComb_K"), posChild.p(),
                            std::sqrt(posChild.tpcNSigmaKa() * posChild.tpcNSigmaKa() +
                                      posChild.tofNSigmaKa() * posChild.tofNSigmaKa()));
        FullQaRegistry.fill(HIST("FullDaughterPosQA/nSigmaComb_p"), posChild.p(),
                            std::sqrt(posChild.tpcNSigmaPr() * posChild.tpcNSigmaPr() +
                                      posChild.tofNSigmaPr() * posChild.tofNSigmaPr()));
        FullQaRegistry.fill(HIST("FullDaughterPosQA/nSigmaComb_d"), posChild.p(),
                            std::sqrt(posChild.tpcNSigmaDe() * posChild.tpcNSigmaDe() +
                                      posChild.tofNSigmaDe() * posChild.tofNSigmaDe()));

        FullQaRegistry.fill(HIST("FullDaughterNegQA/hCharge"), negChild.sign());
        FullQaRegistry.fill(HIST("FullDaughterNegQA/hPt"), negChild.pt());
        FullQaRegistry.fill(HIST("FullDaughterNegQA/hEta"), negChild.eta());
        FullQaRegistry.fill(HIST("FullDaughterNegQA/hPhi"), negChild.phi());
        FullQaRegistry.fill(HIST("FullDaughterNegQA/hDCAxy"), negChild.pt(),
                            negChild.tempFitVar());
        FullQaRegistry.fill(HIST("FullDaughterNegQA/hTPCfindable"),
                            negChild.tpcNClsFindable());
        FullQaRegistry.fill(HIST("FullDaughterNegQA/hTPCfound"), negChild.tpcNClsFound());
        FullQaRegistry.fill(HIST("FullDaughterNegQA/hTPCcrossedOverFindalbe"),
                            negChild.tpcCrossedRowsOverFindableCls());
        FullQaRegistry.fill(HIST("FullDaughterNegQA/hTPCcrossedRows"),
                            negChild.tpcNClsCrossedRows());
        FullQaRegistry.fill(HIST("FullDaughterNegQA/hTPCfindableVsCrossed"),
                            negChild.tpcNClsFindable(), negChild.tpcNClsCrossedRows());
        FullQaRegistry.fill(HIST("FullDaughterNegQA/hTPCshared"), negChild.tpcNClsShared());
        FullQaRegistry.fill(HIST("FullDaughterNegQA/hITSclusters"), negChild.itsNCls());
        FullQaRegistry.fill(HIST("FullDaughterNegQA/hITSclustersIB"),
                            negChild.itsNClsInnerBarrel());
        FullQaRegistry.fill(HIST("FullDaughterNegQA/hDCAz"), negChild.pt(), negChild.dcaZ());
        FullQaRegistry.fill(
          HIST("FullDaughterNegQA/hDCA"), negChild.pt(),
          std::sqrt(pow(negChild.dcaXY(), 2.) + pow(negChild.dcaZ(), 2.)));
        FullQaRegistry.fill(HIST("FullDaughterNegQA/hTPCdEdX"), negChild.p(),
                            negChild.tpcSignal());
        FullQaRegistry.fill(HIST("FullDaughterNegQA/nSigmaTPC_el"), negChild.p(),
                            negChild.tpcNSigmaEl());
        FullQaRegistry.fill(HIST("FullDaughterNegQA/nSigmaTPC_pi"), negChild.p(),
                            negChild.tpcNSigmaPi());
        FullQaRegistry.fill(HIST("FullDaughterNegQA/nSigmaTPC_K"), negChild.p(),
                            negChild.tpcNSigmaKa());
        FullQaRegistry.fill(HIST("FullDaughterNegQA/nSigmaTPC_p"), negChild.p(),
                            negChild.tpcNSigmaPr());
        FullQaRegistry.fill(HIST("FullDaughterNegQA/nSigmaTPC_d"), negChild.p(),
                            negChild.tpcNSigmaDe());
        FullQaRegistry.fill(HIST("FullDaughterNegQA/nSigmaTOF_el"), negChild.p(),
                            negChild.tofNSigmaEl());
        FullQaRegistry.fill(HIST("FullDaughterNegQA/nSigmaTOF_pi"), negChild.p(),
                            negChild.tofNSigmaPi());
        FullQaRegistry.fill(HIST("FullDaughterNegQA/nSigmaTOF_K"), negChild.p(),
                            negChild.tofNSigmaKa());
        FullQaRegistry.fill(HIST("FullDaughterNegQA/nSigmaTOF_p"), negChild.p(),
                            negChild.tofNSigmaPr());
        FullQaRegistry.fill(HIST("FullDaughterNegQA/nSigmaTOF_d"), negChild.p(),
                            negChild.tofNSigmaDe());
        FullQaRegistry.fill(HIST("FullDaughterNegQA/nSigmaComb_el"), negChild.p(),
                            std::sqrt(negChild.tpcNSigmaEl() * negChild.tpcNSigmaEl() +
                                      negChild.tofNSigmaEl() * negChild.tofNSigmaEl()));
        FullQaRegistry.fill(HIST("FullDaughterNegQA/nSigmaComb_pi"), negChild.p(),
                            std::sqrt(negChild.tpcNSigmaPi() * negChild.tpcNSigmaPi() +
                                      negChild.tofNSigmaPi() * negChild.tofNSigmaPi()));
        FullQaRegistry.fill(HIST("FullDaughterNegQA/nSigmaComb_K"), negChild.p(),
                            std::sqrt(negChild.tpcNSigmaKa() * negChild.tpcNSigmaKa() +
                                      negChild.tofNSigmaKa() * negChild.tofNSigmaKa()));
        FullQaRegistry.fill(HIST("FullDaughterNegQA/nSigmaComb_p"), negChild.p(),
                            std::sqrt(negChild.tpcNSigmaPr() * negChild.tpcNSigmaPr() +
                                      negChild.tofNSigmaPr() * negChild.tofNSigmaPr()));
        FullQaRegistry.fill(HIST("FullDaughterNegQA/nSigmaComb_d"), negChild.p(),
                            std::sqrt(negChild.tpcNSigmaDe() * negChild.tpcNSigmaDe() +
                                      negChild.tofNSigmaDe() * negChild.tofNSigmaDe()));
      }
    }

    // for(auto& p: parts){}
    // LOG(info) << "MotherIndex: " << part.index();
    // LOG(info) << "MotherGlobalIndex: " << part.globalIndex();
    // uint64_t id1 = static_cast<uint64_t>(part.childrenIds()[0]);
    // uint64_t id2 = static_cast<uint64_t>(part.childrenIds()[1]);
    //
    // LOG(info) << "Index from mother: " << id1;
    // LOG(info) << "Index from mother: " << id2;
    //
    // uint64_t begin = static_cast<uint64_t>(parts.begin().index());
    // LOG(info) << "StartIndex: " << begin;
    //
    // LOG(info) << "New Index from mother: " << id1 - begin;
    // LOG(info) << "New Index from mother: " << id2 - begin;
    //
    // for (const auto& child : part.children_as<o2::aod::FemtoDreamParticles>()) {
    //   LOG(info) << "ChildIndex: " << static_cast<uint64_t>(child.index());
    //   LOG(info) << "New ChildIndex from mother: " << child.index() - begin;
    //   LOG(info) << "ChildGlobalIndex: " << static_cast<uint64_t>(child.globalIndex());
    //   LOG(info) << "New ChildGlobalIndex: " << child.globalIndex() - begin;
    // LOG(info) << allParts.pt();
    // }
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
