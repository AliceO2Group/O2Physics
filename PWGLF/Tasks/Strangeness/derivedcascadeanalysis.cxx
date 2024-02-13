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
//
/// \ post processing for Cascade analysis runing on derived data
/// \author Lucia Anna Tarasovicova (lucia.anna.husova@cern.ch)

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"

#include <TFile.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <Math/Vector4D.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <cmath>
#include <array>
#include <cstdlib>
#include "Framework/ASoAHelpers.h"

// constants
const float ctauxiPDG = 4.91;     // from PDG
const float ctauomegaPDG = 2.461; // from PDG

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

struct derivedCascadeAnalysis {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<float> zVertexCut{"zVertexCut", 10, "Cut on PV position"};

  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for QA histograms"};

  ConfigurableAxis vertexZ{"vertexZ", {30, -15.0f, 15.0f}, ""};
  ConfigurableAxis axisXiMass{"axisXiMass", {200, 1.222f, 1.422f}, ""};
  ConfigurableAxis axisOmegaMass{"axisOmegaMass", {200, 1.572f, 1.772f}, ""};

  Configurable<bool> isXi{"isXi", 1, "Apply cuts for Xi identification"};
  Configurable<bool> isMC{"isMC", false, "MC data are processed"};
  Configurable<bool> doPtDepCutStudy{"doPtDepCutStudy", false, "Fill histogram with a cutting paramer"};

  Configurable<float> minPt{"minPt", 0.0f, "minPt"};
  Configurable<float> masswin{"masswin", 0.05, "Mass window limit"};
  Configurable<float> lambdaMassWin{"lambdaMassWin", 0.005, "V0 Mass window limit"};
  Configurable<float> rapCut{"rapCut", 0.5, "Rapidity acceptance"};
  Configurable<float> etaDauCut{"etaDauCut", 0.8, "Pseudorapidity acceptance of the cascade daughters"};
  Configurable<float> dcaBaryonToPV{"dcaBaryonToPV", 0.05, "DCA of baryon doughter track To PV"};
  Configurable<float> dcaMesonToPV{"dcaMesonToPV", 0.1, "DCA of meson doughter track To PV"};
  Configurable<float> dcaBachToPV{"dcaBachToPV", 0.04, "DCA Bach To PV"};
  Configurable<double> casccospa{"casccospa", 0.97, "Casc CosPA"};
  Configurable<double> v0cospa{"v0cospa", 0.97, "V0 CosPA"};
  Configurable<float> dcacascdau{"dcacascdau", 1., "DCA Casc Daughters"};
  Configurable<float> dcav0dau{"dcav0dau", 1.5, "DCA V0 Daughters"};
  Configurable<float> dcaV0ToPV{"dcaV0ToPV", 0.06, "DCA V0 To PV"};
  Configurable<float> minRadius{"minRadius", 1.4f, "minRadius"};
  Configurable<float> maxRadius{"maxRadius", 100.0f, "maxRadius"};
  Configurable<float> minV0Radius{"minV0Radius", 1.2f, "V0 transverse decay radius, minimum"};
  Configurable<float> maxV0Radius{"maxV0Radius", 100.0f, "V0 transverse decay radius, maximum"};
  Configurable<float> nsigmatpcPi{"nsigmatpcPi", 5, "N sigma TPC Pion"};
  Configurable<float> nsigmatpcPr{"nsigmatpcPr", 5, "N sigma TPC Proton"};
  Configurable<float> nsigmatpcKa{"nsigmatpcKa", 5, "N sigma TPC Kaon"};
  Configurable<float> bachBaryonCosPA{"bachBaryonCosPA", 0.9999, "Bachelor baryon CosPA"};
  Configurable<float> bachBaryonDCAxyToPV{"bachBaryonDCAxyToPV", 0.05, "DCA bachelor baryon to PV"};
  Configurable<int> mintpccrrows{"mintpccrrows", 50, "min N TPC crossed rows"};
  Configurable<int> dooobrej{"dooobrej", 0, "OOB rejection: 0 no selection, 1 = ITS||TOF, 2 = TOF only for pT > ptthrtof"};
  Configurable<float> ptthrtof{"ptthrtof", 2, "Pt threshold for applying only tof oob rejection"};
  Configurable<float> proplifetime{"proplifetime", 3, "ctau/<ctau>"};
  Configurable<float> rejcomp{"rejcomp", 0.008, "Competing Cascade rejection"};

  Configurable<bool> doPtDepCosPaCut{"doPtDepCosPaCut", false, "Enable pt dependent cos PA cut"};
  Configurable<bool> doPtDepCascRadiusCut{"doPtDepCascRadiusCut", false, "Enable pt dependent cascade radius cut"};
  Configurable<bool> doPtDepV0RadiusCut{"doPtDepV0RadiusCut", false, "Enable pt dependent V0 radius cut"};
  Configurable<bool> doPtDepV0CosPaCut{"doPtDepV0CosPaCut", false, "Enable pt dependent cos PA cut of the V0 daughter"};
  Configurable<bool> doPtDepDCAcascDauCut{"doPtDepDCAcascDauCut", false, "Enable pt dependent DCA cascade daughter cut"};
  Configurable<bool> doDCAdauToPVCut{"doDCAdauToPVCut", true, "Enable cut DCA daughter track to PV"};
  Configurable<bool> doCascadeCosPaCut{"doCascadeCosPaCut", true, "Enable cos PA cut"};
  Configurable<bool> doV0CosPaCut{"doV0CosPaCut", true, "Enable cos PA cut for the V0 daughter"};
  Configurable<bool> doDCACascadeDauCut{"doDCACascadeDauCut", true, "Enable cut DCA betweenn daughter tracks"};
  Configurable<bool> doDCAV0DauCut{"doDCAV0DauCut", true, "Enable cut DCA betweenn V0 daughter tracks"};
  Configurable<bool> doCascadeRadiusCut{"doCascadeRadiusCut", true, "Enable cut on the cascade radius"};
  Configurable<bool> doV0RadiusCut{"doV0RadiusCut", true, "Enable cut on the V0 radius"};
  Configurable<bool> doDCAV0ToPVCut{"doDCAV0ToPVCut", true, "Enable cut DCA of V0 to PV"};
  Configurable<bool> doNTPCSigmaCut{"doNTPCSigmaCut", false, "Enable cut N sigma TPC"};
  Configurable<bool> doBachelorBaryonCut{"doBachelorBaryonCut", true, "Enable Bachelor-Baryon cut "};
  Configurable<bool> doProperLifeTimeCut{"doProperLifeTimeCut", true, "Enable proper life-time cut "};

  Configurable<float> cosPApar0{"cosPApar0", 0.247403, "const par for pt dep cosPA cut"};
  Configurable<float> cosPApar1{"cosPApar1", -0.068957, "linear par for pt dep cosPA cut"};
  Configurable<float> cosPApar2{"cosPApar2", 0.004934, "quadratic par for pt dep cosPA cut"};

  Configurable<float> cosPAV0par{"cosPAV0par", 0.162740, "par for pt dep V0cosPA cut"};

  Configurable<float> parCascRadius0{"parCascRadius0", 1.216159, "const par for pt dep radius cut"};
  Configurable<float> parCascRadius1{"parCascRadius1", 0.064462, "linear par for pt dep radius cut"};

  Configurable<float> parV0Radius0{"parV0Radius0", 2.136381, "const par for pt dep V0 radius cut"};
  Configurable<float> parV0Radius1{"parV0Radius1", 0.437074, "linear par for pt dep V0 radius cut"};

  Configurable<float> dcaCacsDauPar{"dcaCacsDauPar", 0.404424, " par for pt dep DCA cascade daughter cut"};

  Service<o2::framework::O2DatabasePDG> pdgDB;

  void init(InitContext const&)
  {
    histos.add("hEventVertexZ", "hEventVertexZ", kTH1F, {vertexZ});
    histos.add("hEventCentrality", "hEventCentrality", kTH1F, {{101, 0, 101}});
    histos.add("hEventSelection", "hEventSelection", kTH1F, {{4, 0, 4}});

    histos.add("hCandidate", "hCandidate", HistType::kTH1D, {{22, -0.5, 21.5}});

    TString CutLabel[22] = {"All", "MassWin", "y", "DCACascDau", "DCAV0Dau", "rCasc", "rCascMax", "rV0", "rV0Max", "LambdaMass", "Bach-baryon", "V0CosPA", "CompDecayMass", "DCADauToPV", "EtaDau", "CascCosPA", "DCAV0ToPV", "nSigmaTPCV0Dau", "NTPCrows", "OOBRej", "nSigmaTPCbachelor", "ctau"};
    for (Int_t i = 1; i <= histos.get<TH1>(HIST("hCandidate"))->GetNbinsX(); i++) {
      histos.get<TH1>(HIST("hCandidate"))->GetXaxis()->SetBinLabel(i, CutLabel[i - 1]);
    }

    histos.add("InvMassBefSel/hNegativeCascade", "hNegativeCascade", HistType::kTH3F, {axisPt, axisXiMass, {101, 0, 101}});
    histos.add("InvMassBefSel/hPositiveCascade", "hPositiveCascade", {HistType::kTH3F, {axisPt, axisXiMass, {101, 0, 101}}});

    if (!isXi) {
      histos.get<TH3>(HIST("InvMassBefSel/hNegativeCascade"))->GetYaxis()->Set(200, 1.572f, 1.772f);
      histos.get<TH3>(HIST("InvMassBefSel/hPositiveCascade"))->GetYaxis()->Set(200, 1.572f, 1.772f);
    }

    histos.addClone("InvMassBefSel/", "InvMassAfterSel/");
    if (isMC)
      histos.addClone("InvMassBefSel/", "InvMassAfterSelMCrecTruth/");

    if (doPtDepCutStudy && !doProperLifeTimeCut) {
      histos.add("PtDepCutStudy/hNegativeCascadeProperLifeTime", "hNegativeCascadeProperLifeTime", HistType::kTH3F, {axisPt, axisXiMass, {100, 0, 10}});
      histos.add("PtDepCutStudy/hPositiveCascadeProperLifeTime", "hPositiveCascadeProperLifeTime", {HistType::kTH3F, {axisPt, axisXiMass, {100, 0, 10}}});
      if (!isXi) {
        histos.get<TH3>(HIST("PtDepCutStudy/hNegativeCascadeProperLifeTime"))->GetYaxis()->Set(200, 1.572f, 1.772f);
        histos.get<TH3>(HIST("PtDepCutStudy/hPositiveCascadeProperLifeTime"))->GetYaxis()->Set(200, 1.572f, 1.772f);
      }
    }
    if (doPtDepCutStudy && !doBachelorBaryonCut) {
      histos.add("PtDepCutStudy/hNegativeBachelorBaryonDCA", "hNegativeBachelorBaryonDCA", HistType::kTH3F, {axisPt, axisXiMass, {40, 0, 1}});
      histos.add("PtDepCutStudy/hPositiveBachelorBaryonDCA", "hPositiveBachelorBaryonDCA", {HistType::kTH3F, {axisPt, axisXiMass, {40, 0, 1}}});
      if (!isXi) {
        histos.get<TH3>(HIST("PtDepCutStudy/hNegativeBachelorBaryonDCA"))->GetYaxis()->Set(200, 1.572f, 1.772f);
        histos.get<TH3>(HIST("PtDepCutStudy/hPositiveBachelorBaryonDCA"))->GetYaxis()->Set(200, 1.572f, 1.772f);
      }
    }
    if (doPtDepCutStudy && !doDCAV0ToPVCut) {
      histos.add("PtDepCutStudy/hNegativeDCAV0ToPV", "hNegativeDCAV0ToPV", HistType::kTH3F, {axisPt, axisXiMass, {40, 0, 1}});
      histos.add("PtDepCutStudy/hPositiveDCAV0ToPV", "hPositiveDCAV0ToPV", {HistType::kTH3F, {axisPt, axisXiMass, {40, 0, 1}}});
      if (!isXi) {
        histos.get<TH3>(HIST("PtDepCutStudy/hNegativeDCAV0ToPV"))->GetYaxis()->Set(200, 1.572f, 1.772f);
        histos.get<TH3>(HIST("PtDepCutStudy/hPositiveDCAV0ToPV"))->GetYaxis()->Set(200, 1.572f, 1.772f);
      }
    }
    if (doPtDepCutStudy && !doV0RadiusCut) {
      histos.add("PtDepCutStudy/hNegativeV0Radius", "hNegativeV0Radius", HistType::kTH3F, {axisPt, axisXiMass, {20, 0, 10}});
      histos.add("PtDepCutStudy/hPositiveV0Radius", "hPositiveV0Radius", {HistType::kTH3F, {axisPt, axisXiMass, {20, 0, 10}}});
      if (!isXi) {
        histos.get<TH3>(HIST("PtDepCutStudy/hNegativeV0Radius"))->GetYaxis()->Set(200, 1.572f, 1.772f);
        histos.get<TH3>(HIST("PtDepCutStudy/hPositiveV0Radius"))->GetYaxis()->Set(200, 1.572f, 1.772f);
      }
    }
    if (doPtDepCutStudy && !doCascadeRadiusCut) {
      histos.add("PtDepCutStudy/hNegativeCascadeRadius", "hNegativeCascadeRadius", HistType::kTH3F, {axisPt, axisXiMass, {50, 0, 5}});
      histos.add("PtDepCutStudy/hPositiveCascadeRadius", "hPositiveCascadeRadius", {HistType::kTH3F, {axisPt, axisXiMass, {50, 0, 5}}});
      if (!isXi) {
        histos.get<TH3>(HIST("PtDepCutStudy/hNegativeCascadeRadius"))->GetYaxis()->Set(200, 1.572f, 1.772f);
        histos.get<TH3>(HIST("PtDepCutStudy/hPositiveCascadeRadius"))->GetYaxis()->Set(200, 1.572f, 1.772f);
      }
    }

    if (doPtDepCutStudy && !doDCAV0DauCut) {
      histos.add("PtDepCutStudy/hNegativeDCAV0Daughters", "hNegativeDCAV0Daughters", HistType::kTH3F, {axisPt, axisXiMass, {50, 0, 5}});
      histos.add("PtDepCutStudy/hPositiveDCAV0Daughters", "hPositiveDCAV0Daughters", {HistType::kTH3F, {axisPt, axisXiMass, {50, 0, 5}}});
      if (!isXi) {
        histos.get<TH3>(HIST("PtDepCutStudy/hNegativeDCAV0Daughters"))->GetYaxis()->Set(200, 1.572f, 1.772f);
        histos.get<TH3>(HIST("PtDepCutStudy/hPositiveDCAV0Daughters"))->GetYaxis()->Set(200, 1.572f, 1.772f);
      }
    }

    if (doPtDepCutStudy && !doDCACascadeDauCut) {
      histos.add("PtDepCutStudy/hNegativeDCACascDaughters", "hNegativeDCACascDaughters", {HistType::kTH3F, {axisPt, axisXiMass, {20, 0, 1}}});
      histos.add("PtDepCutStudy/hPositiveDCACascDaughters", "hPositiveDCACascDaughters", {HistType::kTH3F, {axisPt, axisXiMass, {20, 0, 1}}});
      if (!isXi) {
        histos.get<TH3>(HIST("PtDepCutStudy/hPositiveDCACascDaughters"))->GetYaxis()->Set(200, 1.572f, 1.772f);
        histos.get<TH3>(HIST("PtDepCutStudy/hNegativeDCACascDaughters"))->GetYaxis()->Set(200, 1.572f, 1.772f);
      }
    }

    if (doPtDepCutStudy && !doV0CosPaCut) {
      histos.add("PtDepCutStudy/hNegativeV0pa", "hNegativeV0pa", HistType::kTH3F, {axisPt, axisXiMass, {40, 0, 0.4}});
      histos.add("PtDepCutStudy/hPositiveV0pa", "hPositiveV0pa", {HistType::kTH3F, {axisPt, axisXiMass, {40, 0, 0.4}}});
      if (!isXi) {
        histos.get<TH3>(HIST("PtDepCutStudy/hNegativeV0pa"))->GetYaxis()->Set(200, 1.572f, 1.772f);
        histos.get<TH3>(HIST("PtDepCutStudy/hPositiveV0pa"))->GetYaxis()->Set(200, 1.572f, 1.772f);
      }
    }
    if (doPtDepCutStudy && !doDCAdauToPVCut) {
      histos.add("PtDepCutStudy/hNegativeDCABachelorToPV", "hNegativeDCABachelorToPV", HistType::kTH3F, {axisPt, axisXiMass, {50, 0, 0.5}});
      histos.add("PtDepCutStudy/hNegativeDCABaryonToPV", "hNegativeDCABaryonToPV", HistType::kTH3F, {axisPt, axisXiMass, {50, 0, 0.5}});
      histos.add("PtDepCutStudy/hNegativeDCAMesonToPV", "hNegativeDCAMesonToPV", HistType::kTH3F, {axisPt, axisXiMass, {50, 0, 0.5}});
      histos.add("PtDepCutStudy/hPositiveDCABachelorToPV", "hPositiveDCABachelorToPV", {HistType::kTH3F, {axisPt, axisXiMass, {50, 0, 0.5}}});
      histos.add("PtDepCutStudy/hPositiveDCABaryonToPV", "hPositiveDCABaryonToPV", HistType::kTH3F, {axisPt, axisXiMass, {50, 0, 0.5}});
      histos.add("PtDepCutStudy/hPositiveDCAMesonToPV", "hPositiveDCAMesonToPV", HistType::kTH3F, {axisPt, axisXiMass, {50, 0, 0.5}});
      if (!isXi) {
        histos.get<TH3>(HIST("PtDepCutStudy/hNegativeDCABachelorToPV"))->GetYaxis()->Set(200, 1.572f, 1.772f);
        histos.get<TH3>(HIST("PtDepCutStudy/hNegativeDCABaryonToPV"))->GetYaxis()->Set(200, 1.572f, 1.772f);
        histos.get<TH3>(HIST("PtDepCutStudy/hNegativeDCAMesonToPV"))->GetYaxis()->Set(200, 1.572f, 1.772f);
        histos.get<TH3>(HIST("PtDepCutStudy/hPositiveDCABachelorToPV"))->GetYaxis()->Set(200, 1.572f, 1.772f);
        histos.get<TH3>(HIST("PtDepCutStudy/hPositiveDCABaryonToPV"))->GetYaxis()->Set(200, 1.572f, 1.772f);
        histos.get<TH3>(HIST("PtDepCutStudy/hPositiveDCAMesonToPV"))->GetYaxis()->Set(200, 1.572f, 1.772f);
      }
    }

    if (doPtDepCutStudy && !doCascadeCosPaCut) {
      histos.add("PtDepCutStudy/hNegativeCascPA", "hNegativeCascPA", HistType::kTH3F, {axisPt, axisXiMass, {40, 0, 0.4}});
      histos.add("PtDepCutStudy/hPositiveCascPA", "hPositiveCascPA", {HistType::kTH3F, {axisPt, axisXiMass, {40, 0, 0.4}}});
      if (!isXi) {
        histos.get<TH3>(HIST("PtDepCutStudy/hNegativeCascPA"))->GetYaxis()->Set(200, 1.572f, 1.772f);
        histos.get<TH3>(HIST("PtDepCutStudy/hPositiveCascPA"))->GetYaxis()->Set(200, 1.572f, 1.772f);
      }
    }

    if (isMC) {
      histos.addClone("PtDepCutStudy/", "PtDepCutStudyMCTruth/");
      histos.add("hNegativeCascadePtForEfficiency", "hNegativeCascadePtForEfficiency", HistType::kTH3F, {axisPt, axisXiMass, {101, 0, 101}});
      histos.add("hPositiveCascadePtForEfficiency", "hPositiveCascadePtForEfficiency", {HistType::kTH3F, {axisPt, axisXiMass, {101, 0, 101}}});
      if (!isXi) {
        histos.get<TH3>(HIST("hNegativeCascadePtForEfficiency"))->GetYaxis()->Set(200, 1.572f, 1.772f);
        histos.get<TH3>(HIST("hPositiveCascadePtForEfficiency"))->GetYaxis()->Set(200, 1.572f, 1.772f);
      }

      histos.add("h2dNVerticesVsCentrality", "h2dNVerticesVsCentrality", HistType::kTH2F, {{101, 0, 101}, {20, -0.5, 19.5}});
      histos.add("hGenNegativeXi", "hGenNegativeXi", HistType::kTH2F, {{101, 0, 101}, axisPt});
      histos.add("hGenPositiveXi", "hGenPositiveXi", HistType::kTH2F, {{101, 0, 101}, axisPt});
      histos.add("hGenNegativeOmega", "hGenNegativeOmega", HistType::kTH2F, {{101, 0, 101}, axisPt});
      histos.add("hGenPositiveOmega", "hGenPositiveOmega", HistType::kTH2F, {{101, 0, 101}, axisPt});
    }
  }
  template <typename TCascade>
  bool IsCosPAAccepted(TCascade casc, float x, float y, float z, bool ptdepcut, bool isCascPa)
  {

    if (ptdepcut) {
      double ptdepCut;
      if (isCascPa)
        ptdepCut = cosPApar0 + cosPApar1 * casc.pt() + cosPApar2 * TMath::Power(casc.pt(), 2);
      else
        ptdepCut = cosPAV0par / casc.pt();
      if (ptdepCut > 0.3 && casc.pt() < 0.5)
        ptdepCut = 0.3;
      if (ptdepCut < 0.012 || casc.pt() > 7)
        ptdepCut = 0.012;
      if (isCascPa && casc.casccosPA(x, y, z) < TMath::Cos(ptdepCut))
        return false;
      if (!isCascPa && casc.v0cosPA(x, y, z) < TMath::Cos(ptdepCut))
        return false;
    } else if (isCascPa && casc.casccosPA(x, y, z) < casccospa)
      return false;
    else if (!isCascPa && casc.v0cosPA(x, y, z) < v0cospa)
      return false;

    return true;
  }
  template <typename TCollision>
  bool IsEventAccepted(TCollision coll, bool sel)
  {
    histos.fill(HIST("hEventSelection"), 0.5 /* all collisions */);
    if (!sel) {
      return false;
    }
    histos.fill(HIST("hEventSelection"), 1.5 /* collisions  after sel*/);
    if (TMath::Abs(coll.posZ()) > zVertexCut) {
      return false;
    }
    histos.fill(HIST("hEventVertexZ"), coll.posZ());
    histos.fill(HIST("hEventSelection"), 2.5 /* collisions  after sel pvz sel*/);

    if (coll.centFT0C() > 100) {
      return false;
    }
    histos.fill(HIST("hEventCentrality"), coll.centFT0C());
    histos.fill(HIST("hEventSelection"), 3.5 /* collisions  after sel centrality sel*/);
    return true;
  }

  template <typename TCascade>
  bool IsCascadeCandidateAccepted(TCascade casc, int counter, float centrality)
  {

    if (isXi) {
      if (TMath::Abs(casc.mXi() - pdgDB->Mass(3312)) > masswin) {
        return false;
      }
      histos.fill(HIST("hCandidate"), ++counter);
      if (TMath::Abs(casc.yXi()) > rapCut)
        return false;
      histos.fill(HIST("hCandidate"), ++counter);
    } else {
      if (TMath::Abs(casc.mOmega() - pdgDB->Mass(3334)) > masswin) {
        return false;
      }
      histos.fill(HIST("hCandidate"), ++counter);
      if (TMath::Abs(casc.yOmega()) > rapCut)
        return false;
      histos.fill(HIST("hCandidate"), ++counter);
    }

    if (doDCACascadeDauCut) {
      if (doPtDepDCAcascDauCut) {
        float ptDepCut = dcaCacsDauPar / casc.pt();
        if (casc.dcacascdaughters() > ptDepCut)
          return false;
      } else if (casc.dcacascdaughters() > dcacascdau)
        return false;
      histos.fill(HIST("hCandidate"), ++counter);
    } else
      ++counter;

    if (doDCAV0DauCut) {
      if (casc.dcaV0daughters() > dcav0dau)
        return false;
      histos.fill(HIST("hCandidate"), ++counter);
    } else
      ++counter;

    if (doCascadeRadiusCut) {
      if (doPtDepCascRadiusCut) {
        double ptdepminRadius = parCascRadius0 + parCascRadius1 * casc.pt();
        if (casc.cascradius() < ptdepminRadius)
          return false;
      } else if (casc.cascradius() < minRadius)
        return false;
      histos.fill(HIST("hCandidate"), ++counter);

      if (casc.cascradius() > maxRadius)
        return false;
      histos.fill(HIST("hCandidate"), ++counter);
    } else
      counter += 2;

    if (doV0RadiusCut) {
      if (doPtDepV0RadiusCut) {
        float cut = parV0Radius0 + casc.pt() * parV0Radius1;
        if (casc.v0radius() < cut)
          return false;
      } else if (casc.v0radius() < minV0Radius)
        return false;
      histos.fill(HIST("hCandidate"), ++counter);
      if (casc.v0radius() > maxV0Radius)
        return false;
      histos.fill(HIST("hCandidate"), ++counter);
    } else
      counter += 2;

    if (TMath::Abs(casc.mLambda() - pdgDB->Mass(3122)) > lambdaMassWin)
      return false;
    histos.fill(HIST("hCandidate"), ++counter);

    if (doBachelorBaryonCut) {
      if ((casc.bachBaryonCosPA() > bachBaryonCosPA || TMath::Abs(casc.bachBaryonDCAxyToPV()) < bachBaryonDCAxyToPV)) { // Bach-baryon selection if required
        return false;
      }
      histos.fill(HIST("hCandidate"), ++counter);
    } else
      ++counter;

    if (doV0CosPaCut) {
      if (!IsCosPAAccepted(casc, casc.x(), casc.y(), casc.z(), doPtDepV0CosPaCut, false))
        return false;
      histos.fill(HIST("hCandidate"), ++counter);
    } else
      ++counter;

    if (isXi) {
      if (TMath::Abs(casc.mOmega() - pdgDB->Mass(3334)) < rejcomp)
        return false;
      histos.fill(HIST("hCandidate"), ++counter);
    } else {
      if (TMath::Abs(casc.mXi() - pdgDB->Mass(3312)) < rejcomp)
        return false;
      histos.fill(HIST("hCandidate"), ++counter);
    }

    if (doDCAdauToPVCut) {
      if (TMath::Abs(casc.dcabachtopv()) < dcaBachToPV)
        return false;
      if (casc.sign() > 0 && (TMath::Abs(casc.dcanegtopv()) < dcaBaryonToPV || TMath::Abs(casc.dcapostopv()) < dcaMesonToPV))
        return false;
      if (casc.sign() < 0 && (TMath::Abs(casc.dcapostopv()) < dcaBaryonToPV || TMath::Abs(casc.dcanegtopv()) < dcaMesonToPV))
        return false;
      histos.fill(HIST("hCandidate"), ++counter);
    } else
      ++counter;

    return true;
  }

  void processCascades(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels>::iterator const& coll, soa::Join<aod::CascCollRefs, aod::CascCores, aod::CascExtras, aod::CascBBs> const& Cascades, soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs> const&)
  {

    if (!IsEventAccepted(coll, coll.sel8()))
      return;

    for (auto& casc : Cascades) {

      int counter = -1;
      histos.fill(HIST("hCandidate"), ++counter);

      double invmass;
      if (isXi)
        invmass = casc.mXi();
      else
        invmass = casc.mOmega();
      // To have trace of how it was before selections
      if (casc.sign() < 0) {
        histos.fill(HIST("InvMassBefSel/hNegativeCascade"), casc.pt(), invmass, coll.centFT0C());
      }
      if (casc.sign() > 0) {
        histos.fill(HIST("InvMassBefSel/hPositiveCascade"), casc.pt(), invmass, coll.centFT0C());
      }

      if (!IsCascadeCandidateAccepted(casc, counter, coll.centFT0C()))
        continue;
      counter += 13;

      auto negExtra = casc.negTrackExtra_as<soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>>();
      auto posExtra = casc.posTrackExtra_as<soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>>();
      auto bachExtra = casc.bachTrackExtra_as<soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>>();

      auto poseta = RecoDecay::eta(std::array{casc.pxpos(), casc.pypos(), casc.pzpos()});
      auto negeta = RecoDecay::eta(std::array{casc.pxneg(), casc.pyneg(), casc.pzneg()});
      auto bacheta = RecoDecay::eta(std::array{casc.pxbach(), casc.pybach(), casc.pzbach()});
      if (TMath::Abs(poseta) > etaDauCut || TMath::Abs(negeta) > etaDauCut || TMath::Abs(bacheta) > etaDauCut)
        continue;
      histos.fill(HIST("hCandidate"), ++counter);

      if (doCascadeCosPaCut) {
        if (!IsCosPAAccepted(casc, coll.posX(), coll.posY(), coll.posZ(), doPtDepCosPaCut, true))
          continue;
        histos.fill(HIST("hCandidate"), ++counter);
      } else
        ++counter;

      if (doDCAV0ToPVCut) {
        if (TMath::Abs(casc.dcav0topv(coll.posX(), coll.posY(), coll.posZ())) < dcaV0ToPV)
          continue;
        histos.fill(HIST("hCandidate"), ++counter);
      } else
        ++counter;

      if (doNTPCSigmaCut) {
        if (casc.sign() < 0) {
          if (TMath::Abs(posExtra.tpcNSigmaPr()) > nsigmatpcPr || TMath::Abs(negExtra.tpcNSigmaPi()) > nsigmatpcPi)
            continue;
          histos.fill(HIST("hCandidate"), ++counter);
        } else if (casc.sign() > 0) {
          if (TMath::Abs(posExtra.tpcNSigmaPi()) > nsigmatpcPi || TMath::Abs(negExtra.tpcNSigmaPr()) > nsigmatpcPr)
            continue;
          histos.fill(HIST("hCandidate"), ++counter);
        }
      } else
        ++counter;

      if (posExtra.tpcCrossedRows() < mintpccrrows || negExtra.tpcCrossedRows() < mintpccrrows || bachExtra.tpcCrossedRows() < mintpccrrows)
        continue;
      histos.fill(HIST("hCandidate"), ++counter);

      bool kHasTOF = (posExtra.hasTOF() || negExtra.hasTOF() || bachExtra.hasTOF());
      bool kHasITS = (posExtra.hasITS() || negExtra.hasITS() || bachExtra.hasITS());
      if (dooobrej == 1) {
        if (!kHasTOF && !kHasITS)
          continue;
        histos.fill(HIST("hCandidate"), ++counter);
      } else if (dooobrej == 2) {
        if (!kHasTOF && (casc.pt() > ptthrtof))
          continue;
        histos.fill(HIST("hCandidate"), ++counter);
      } else {
        ++counter;
      }

      float cascpos = std::hypot(casc.x() - coll.posX(), casc.y() - coll.posY(), casc.z() - coll.posZ());
      float cascptotmom = std::hypot(casc.px(), casc.py(), casc.pz());
      float ctau = -10;

      if (isXi) {
        if (doNTPCSigmaCut) {
          if (TMath::Abs(bachExtra.tpcNSigmaPi()) > nsigmatpcPi)
            continue;
          histos.fill(HIST("hCandidate"), ++counter);
        } else
          ++counter;

        ctau = pdgDB->Mass(3312) * cascpos / ((cascptotmom + 1e-13) * ctauxiPDG);
        if (doProperLifeTimeCut) {
          if (ctau > proplifetime)
            continue;
          histos.fill(HIST("hCandidate"), ++counter);
        } else {
          ++counter;
        }
      } else {
        if (doNTPCSigmaCut) {
          if (TMath::Abs(bachExtra.tpcNSigmaKa()) > nsigmatpcKa)
            continue;
          histos.fill(HIST("hCandidate"), ++counter);
        } else
          ++counter;

        ctau = pdgDB->Mass(3334) * cascpos / ((cascptotmom + 1e-13) * ctauomegaPDG);
        if (doProperLifeTimeCut) {
          if (ctau > proplifetime)
            continue;
          histos.fill(HIST("hCandidate"), ++counter);
        } else
          ++counter;
      }

      if (casc.sign() < 0) {
        histos.fill(HIST("InvMassAfterSel/hNegativeCascade"), casc.pt(), invmass, coll.centFT0C());
        if (!doBachelorBaryonCut && doPtDepCutStudy)
          histos.fill(HIST("PtDepCutStudy/hNegativeBachelorBaryonDCA"), casc.pt(), invmass, casc.bachBaryonDCAxyToPV());
        if (!doDCAV0ToPVCut && doPtDepCutStudy)
          histos.fill(HIST("PtDepCutStudy/hNegativeDCAV0ToPV"), casc.pt(), invmass, TMath::Abs(casc.dcav0topv(casc.x(), casc.y(), casc.z())));
        if (!doV0RadiusCut && doPtDepCutStudy)
          histos.fill(HIST("PtDepCutStudy/hNegativeV0Radius"), casc.pt(), invmass, casc.v0radius());
        if (!doCascadeRadiusCut && doPtDepCutStudy)
          histos.fill(HIST("PtDepCutStudy/hNegativeCascadeRadius"), casc.pt(), invmass, casc.cascradius());
        if (!doDCAV0DauCut && doPtDepCutStudy)
          histos.fill(HIST("PtDepCutStudy/hNegativeDCAV0Daughters"), casc.pt(), invmass, casc.dcaV0daughters());
        if (!doDCACascadeDauCut && doPtDepCutStudy)
          histos.fill(HIST("PtDepCutStudy/hNegativeDCACascDaughters"), casc.pt(), invmass, casc.dcacascdaughters());
        if (!doV0CosPaCut && doPtDepCutStudy)
          histos.fill(HIST("PtDepCutStudy/hNegativeV0pa"), casc.pt(), invmass, TMath::ACos(casc.v0cosPA(casc.x(), casc.y(), casc.z())));
        if (!doCascadeCosPaCut && doPtDepCutStudy)
          histos.fill(HIST("PtDepCutStudy/hNegativeCascPA"), casc.pt(), invmass, TMath::ACos(casc.casccosPA(coll.posX(), coll.posY(), coll.posZ())));
        if (!doDCAdauToPVCut && doPtDepCutStudy) {
          histos.fill(HIST("PtDepCutStudy/hNegativeDCABachelorToPV"), casc.pt(), invmass, casc.dcabachtopv());
          histos.fill(HIST("PtDepCutStudy/hNegativeDCAMesonToPV"), casc.pt(), invmass, casc.dcanegtopv());
          histos.fill(HIST("PtDepCutStudy/hNegativeDCABaryonToPV"), casc.pt(), invmass, casc.dcapostopv());
        }
        if (!doProperLifeTimeCut && doPtDepCutStudy)
          histos.fill(HIST("PtDepCutStudy/hNegativeCascadeProperLifeTime"), casc.pt(), invmass, ctau);
      } else {
        histos.fill(HIST("InvMassAfterSel/hPositiveCascade"), casc.pt(), invmass, coll.centFT0C());
        if (!doBachelorBaryonCut && doPtDepCutStudy)
          histos.fill(HIST("PtDepCutStudy/hPositiveBachelorBaryonDCA"), casc.pt(), invmass, casc.bachBaryonDCAxyToPV());
        if (!doDCAV0ToPVCut && doPtDepCutStudy)
          histos.fill(HIST("PtDepCutStudy/hPositiveDCAV0ToPV"), casc.pt(), invmass, TMath::Abs(casc.dcav0topv(casc.x(), casc.y(), casc.z())));
        if (!doV0RadiusCut && doPtDepCutStudy)
          histos.fill(HIST("PtDepCutStudy/hPositiveV0Radius"), casc.pt(), invmass, casc.v0radius());
        if (!doCascadeRadiusCut && doPtDepCutStudy)
          histos.fill(HIST("PtDepCutStudy/hPositiveCascadeRadius"), casc.pt(), invmass, casc.cascradius());
        if (!doDCAV0DauCut && doPtDepCutStudy)
          histos.fill(HIST("PtDepCutStudy/hPositiveDCAV0Daughters"), casc.pt(), invmass, casc.dcaV0daughters());
        if (!doDCACascadeDauCut && doPtDepCutStudy)
          histos.fill(HIST("PtDepCutStudy/hPositiveDCACascDaughters"), casc.pt(), invmass, casc.dcacascdaughters());
        if (!doV0CosPaCut && doPtDepCutStudy)
          histos.fill(HIST("PtDepCutStudy/hPositiveV0pa"), casc.pt(), invmass, TMath::ACos(casc.v0cosPA(casc.x(), casc.y(), casc.z())));
        if (!doCascadeCosPaCut && doPtDepCutStudy)
          histos.fill(HIST("PtDepCutStudy/hPositiveCascPA"), casc.pt(), invmass, TMath::ACos(casc.casccosPA(coll.posX(), coll.posY(), coll.posZ())));
        if (!doDCAdauToPVCut && doPtDepCutStudy) {
          histos.fill(HIST("PtDepCutStudy/hPositiveDCABachelorToPV"), casc.pt(), invmass, casc.dcabachtopv());
          histos.fill(HIST("PtDepCutStudy/hPositiveDCAMesonToPV"), casc.pt(), invmass, casc.dcapostopv());
          histos.fill(HIST("PtDepCutStudy/hPositiveDCABaryonToPV"), casc.pt(), invmass, casc.dcanegtopv());
        }
        if (!doProperLifeTimeCut && doPtDepCutStudy)
          histos.fill(HIST("PtDepCutStudy/hPositiveCascadeProperLifeTime"), casc.pt(), invmass, ctau);
      }
    }
  }
  void processCascadesMCrec(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels>::iterator const& coll, soa::Join<aod::CascCollRefs, aod::CascCores, aod::CascMCCores, aod::CascExtras, aod::CascBBs> const& Cascades, soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs> const&)
  {
    if (!IsEventAccepted(coll, coll.sel8()))
      return;

    for (auto& casc : Cascades) {

      int counter = -1;
      histos.fill(HIST("hCandidate"), ++counter);

      // To have trace of how it was before selections
      if (casc.sign() < 0) {
        if (isXi)
          histos.fill(HIST("InvMassBefSel/hNegativeCascade"), casc.pt(), casc.mXi(), coll.centFT0C());
        else
          histos.fill(HIST("InvMassBefSel/hNegativeCascade"), casc.pt(), casc.mOmega(), coll.centFT0C());
      }
      if (casc.sign() > 0) {
        if (isXi)
          histos.fill(HIST("InvMassBefSel/hPositiveCascade"), casc.pt(), casc.mXi(), coll.centFT0C());
        else
          histos.fill(HIST("InvMassBefSel/hPositiveCascade"), casc.pt(), casc.mOmega(), coll.centFT0C());
      }

      if (!IsCascadeCandidateAccepted(casc, counter, coll.centFT0C()))
        continue;
      counter += 13;

      auto negExtra = casc.negTrackExtra_as<soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>>();
      auto posExtra = casc.posTrackExtra_as<soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>>();
      auto bachExtra = casc.bachTrackExtra_as<soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>>();

      auto poseta = RecoDecay::eta(std::array{casc.pxpos(), casc.pypos(), casc.pzpos()});
      auto negeta = RecoDecay::eta(std::array{casc.pxneg(), casc.pyneg(), casc.pzneg()});
      auto bacheta = RecoDecay::eta(std::array{casc.pxbach(), casc.pybach(), casc.pzbach()});
      if (TMath::Abs(poseta) > etaDauCut || TMath::Abs(negeta) > etaDauCut || TMath::Abs(bacheta) > etaDauCut)
        continue;
      histos.fill(HIST("hCandidate"), ++counter);

      if (doCascadeCosPaCut) {
        if (!IsCosPAAccepted(casc, coll.posX(), coll.posY(), coll.posZ(), doPtDepCosPaCut, true))
          continue;
        histos.fill(HIST("hCandidate"), ++counter);
      } else
        ++counter;

      if (doDCAV0ToPVCut) {
        if (TMath::Abs(casc.dcav0topv(coll.posX(), coll.posY(), coll.posZ())) < dcaV0ToPV)
          continue;
        histos.fill(HIST("hCandidate"), ++counter);
      } else
        ++counter;

      if (doNTPCSigmaCut) {
        if (casc.sign() < 0) {
          if (TMath::Abs(posExtra.tpcNSigmaPr()) > nsigmatpcPr || TMath::Abs(negExtra.tpcNSigmaPi()) > nsigmatpcPi)
            continue;
          histos.fill(HIST("hCandidate"), ++counter);
        } else if (casc.sign() > 0) {
          if (TMath::Abs(posExtra.tpcNSigmaPi()) > nsigmatpcPi || TMath::Abs(negExtra.tpcNSigmaPr()) > nsigmatpcPr)
            continue;
          histos.fill(HIST("hCandidate"), ++counter);
        }
      } else
        ++counter;

      if (posExtra.tpcCrossedRows() < mintpccrrows || negExtra.tpcCrossedRows() < mintpccrrows || bachExtra.tpcCrossedRows() < mintpccrrows)
        continue;
      histos.fill(HIST("hCandidate"), ++counter);

      bool kHasTOF = (posExtra.hasTOF() || negExtra.hasTOF() || bachExtra.hasTOF());
      bool kHasITS = (posExtra.hasITS() || negExtra.hasITS() || bachExtra.hasITS());
      if (dooobrej == 1) {
        if (!kHasTOF && !kHasITS)
          continue;
        histos.fill(HIST("hCandidate"), ++counter);
      } else if (dooobrej == 2) {
        if (!kHasTOF && (casc.pt() > ptthrtof))
          continue;
        histos.fill(HIST("hCandidate"), ++counter);
      } else {
        ++counter;
      }

      double invmass;
      float cascpos = std::hypot(casc.x() - coll.posX(), casc.y() - coll.posY(), casc.z() - coll.posZ());
      float cascptotmom = std::hypot(casc.px(), casc.py(), casc.pz());
      float ctau = -10;

      if (isXi) {
        if (doNTPCSigmaCut) {
          if (TMath::Abs(bachExtra.tpcNSigmaPi()) > nsigmatpcPi)
            continue;
          histos.fill(HIST("hCandidate"), ++counter);
        } else
          ++counter;

        ctau = pdgDB->Mass(3312) * cascpos / ((cascptotmom + 1e-13) * ctauxiPDG);
        if (doProperLifeTimeCut) {
          if (ctau > proplifetime)
            continue;
          histos.fill(HIST("hCandidate"), ++counter);
        } else {
          ++counter;
        }

        invmass = casc.mXi();
      } else {
        if (doNTPCSigmaCut) {
          if (TMath::Abs(bachExtra.tpcNSigmaKa()) > nsigmatpcKa)
            continue;
          histos.fill(HIST("hCandidate"), ++counter);
        } else
          ++counter;
        ctau = pdgDB->Mass(3334) * cascpos / ((cascptotmom + 1e-13) * ctauomegaPDG);
        if (doProperLifeTimeCut) {
          if (ctau > proplifetime)
            continue;
          histos.fill(HIST("hCandidate"), ++counter);
        } else
          ++counter;
        invmass = casc.mOmega();
      }

      if (casc.sign() < 0) {
        histos.fill(HIST("InvMassAfterSel/hNegativeCascade"), casc.pt(), invmass, coll.centFT0C());
        if (!doBachelorBaryonCut && doPtDepCutStudy)
          histos.fill(HIST("PtDepCutStudy/hNegativeBachelorBaryonDCA"), casc.pt(), invmass, casc.bachBaryonDCAxyToPV());
        if (!doDCAV0ToPVCut && doPtDepCutStudy)
          histos.fill(HIST("PtDepCutStudy/hNegativeDCAV0ToPV"), casc.pt(), invmass, TMath::Abs(casc.dcav0topv(casc.x(), casc.y(), casc.z())));
        if (!doV0RadiusCut && doPtDepCutStudy)
          histos.fill(HIST("PtDepCutStudy/hNegativeV0Radius"), casc.pt(), invmass, casc.v0radius());
        if (!doCascadeRadiusCut && doPtDepCutStudy)
          histos.fill(HIST("PtDepCutStudy/hNegativeCascadeRadius"), casc.pt(), invmass, casc.cascradius());
        if (!doDCAV0DauCut && doPtDepCutStudy)
          histos.fill(HIST("PtDepCutStudy/hNegativeDCAV0Daughters"), casc.pt(), invmass, casc.dcaV0daughters());
        if (!doDCACascadeDauCut && doPtDepCutStudy)
          histos.fill(HIST("PtDepCutStudy/hNegativeDCACascDaughters"), casc.pt(), invmass, casc.dcacascdaughters());
        if (!doV0CosPaCut && doPtDepCutStudy)
          histos.fill(HIST("PtDepCutStudy/hNegativeV0pa"), casc.pt(), invmass, TMath::ACos(casc.v0cosPA(casc.x(), casc.y(), casc.z())));
        if (!doCascadeCosPaCut && doPtDepCutStudy)
          histos.fill(HIST("PtDepCutStudy/hNegativeCascPA"), casc.pt(), invmass, TMath::ACos(casc.casccosPA(coll.posX(), coll.posY(), coll.posZ())));
        if (!doDCAdauToPVCut && doPtDepCutStudy) {
          histos.fill(HIST("PtDepCutStudy/hNegativeDCABachelorToPV"), casc.pt(), invmass, casc.dcabachtopv());
          histos.fill(HIST("PtDepCutStudy/hNegativeDCAMesonToPV"), casc.pt(), invmass, casc.dcanegtopv());
          histos.fill(HIST("PtDepCutStudy/hNegativeDCABaryonToPV"), casc.pt(), invmass, casc.dcapostopv());
        }
        if (!doProperLifeTimeCut && doPtDepCutStudy)
          histos.fill(HIST("PtDepCutStudy/hNegativeCascadeProperLifeTime"), casc.pt(), invmass, ctau);
        if (casc.isPhysicalPrimary()) {
          if ((isXi && casc.pdgCode() == 3312) || (!isXi && casc.pdgCode() == 3334)) {
            histos.fill(HIST("InvMassAfterSelMCrecTruth/hNegativeCascade"), casc.pt(), invmass, coll.centFT0C());
            if (!doBachelorBaryonCut && doPtDepCutStudy)
              histos.fill(HIST("PtDepCutStudyMCTruth/hNegativeBachelorBaryonDCA"), casc.pt(), invmass, casc.bachBaryonDCAxyToPV());
            if (!doDCAV0ToPVCut && doPtDepCutStudy)
              histos.fill(HIST("PtDepCutStudyMCTruth/hNegativeDCAV0ToPV"), casc.pt(), invmass, TMath::Abs(casc.dcav0topv(casc.x(), casc.y(), casc.z())));
            if (!doV0RadiusCut && doPtDepCutStudy)
              histos.fill(HIST("PtDepCutStudyMCTruth/hNegativeV0Radius"), casc.pt(), invmass, casc.v0radius());
            if (!doCascadeRadiusCut && doPtDepCutStudy)
              histos.fill(HIST("PtDepCutStudyMCTruth/hNegativeCascadeRadius"), casc.pt(), invmass, casc.cascradius());
            if (!doDCAV0DauCut && doPtDepCutStudy)
              histos.fill(HIST("PtDepCutStudyMCTruth/hNegativeDCAV0Daughters"), casc.pt(), invmass, casc.dcaV0daughters());
            if (!doDCACascadeDauCut && doPtDepCutStudy)
              histos.fill(HIST("PtDepCutStudyMCTruth/hNegativeDCACascDaughters"), casc.pt(), invmass, casc.dcacascdaughters());
            if (!doV0CosPaCut && doPtDepCutStudy)
              histos.fill(HIST("PtDepCutStudyMCTruth/hNegativeV0pa"), casc.pt(), invmass, TMath::ACos(casc.v0cosPA(casc.x(), casc.y(), casc.z())));
            if (!doCascadeCosPaCut && doPtDepCutStudy)
              histos.fill(HIST("PtDepCutStudyMCTruth/hNegativeCascPA"), casc.pt(), invmass, TMath::ACos(casc.casccosPA(coll.posX(), coll.posY(), coll.posZ())));
            if (!doDCAdauToPVCut && doPtDepCutStudy) {
              histos.fill(HIST("PtDepCutStudyMCTruth/hNegativeDCABachelorToPV"), casc.pt(), invmass, casc.dcabachtopv());
              histos.fill(HIST("PtDepCutStudyMCTruth/hNegativeDCAMesonToPV"), casc.pt(), invmass, casc.dcanegtopv());
              histos.fill(HIST("PtDepCutStudyMCTruth/hNegativeDCABaryonToPV"), casc.pt(), invmass, casc.dcapostopv());
            }
            if (!doProperLifeTimeCut && doPtDepCutStudy)
              histos.fill(HIST("PtDepCutStudyMCTruth/hNegativeCascadeProperLifeTime"), casc.pt(), invmass, ctau);
            if ((isXi && casc.pdgCodeV0() == 3122 && casc.pdgCodePositive() == 2212 && casc.pdgCodeNegative() == -211 && casc.pdgCodeBachelor() == -211) || (!isXi && casc.pdgCodeV0() == 3122 && casc.pdgCodePositive() == 2212 && casc.pdgCodeNegative() == -211 && casc.pdgCodeBachelor() == -321))
              histos.fill(HIST("hNegativeCascadePtForEfficiency"), casc.pt(), invmass, coll.centFT0C());
          }
        }
      } else {
        histos.fill(HIST("InvMassAfterSel/hPositiveCascade"), casc.pt(), invmass, coll.centFT0C());
        if (!doBachelorBaryonCut && doPtDepCutStudy)
          histos.fill(HIST("PtDepCutStudy/hPositiveBachelorBaryonDCA"), casc.pt(), invmass, casc.bachBaryonDCAxyToPV());
        if (!doDCAV0ToPVCut && doPtDepCutStudy)
          histos.fill(HIST("PtDepCutStudy/hPositiveDCAV0ToPV"), casc.pt(), invmass, TMath::Abs(casc.dcav0topv(casc.x(), casc.y(), casc.z())));
        if (!doV0RadiusCut && doPtDepCutStudy)
          histos.fill(HIST("PtDepCutStudy/hPositiveV0Radius"), casc.pt(), invmass, casc.v0radius());
        if (!doCascadeRadiusCut && doPtDepCutStudy)
          histos.fill(HIST("PtDepCutStudy/hPositiveCascadeRadius"), casc.pt(), invmass, casc.cascradius());
        if (!doDCAV0DauCut && doPtDepCutStudy)
          histos.fill(HIST("PtDepCutStudy/hPositiveDCAV0Daughters"), casc.pt(), invmass, casc.dcaV0daughters());
        if (!doDCACascadeDauCut && doPtDepCutStudy)
          histos.fill(HIST("PtDepCutStudy/hPositiveDCACascDaughters"), casc.pt(), invmass, casc.dcacascdaughters());
        if (!doV0CosPaCut && doPtDepCutStudy)
          histos.fill(HIST("PtDepCutStudy/hPositiveV0pa"), casc.pt(), invmass, TMath::ACos(casc.v0cosPA(casc.x(), casc.y(), casc.z())));
        if (!doCascadeCosPaCut && doPtDepCutStudy)
          histos.fill(HIST("PtDepCutStudy/hPositiveCascPA"), casc.pt(), invmass, TMath::ACos(casc.casccosPA(coll.posX(), coll.posY(), coll.posZ())));
        if (!doDCAdauToPVCut && doPtDepCutStudy) {
          histos.fill(HIST("PtDepCutStudy/hPositiveDCABachelorToPV"), casc.pt(), invmass, casc.dcabachtopv());
          histos.fill(HIST("PtDepCutStudy/hPositiveDCAMesonToPV"), casc.pt(), invmass, casc.dcapostopv());
          histos.fill(HIST("PtDepCutStudy/hPositiveDCABaryonToPV"), casc.pt(), invmass, casc.dcanegtopv());
        }
        if (!doProperLifeTimeCut && doPtDepCutStudy)
          histos.fill(HIST("PtDepCutStudy/hPositiveCascadeProperLifeTime"), casc.pt(), invmass, ctau);
        if (casc.isPhysicalPrimary()) {
          if ((isXi && casc.pdgCode() == -3312) || (!isXi && casc.pdgCode() == -3334)) {
            histos.fill(HIST("InvMassAfterSelMCrecTruth/hPositiveCascade"), casc.pt(), invmass, coll.centFT0C());
            if (!doBachelorBaryonCut && doPtDepCutStudy)
              histos.fill(HIST("PtDepCutStudyMCTruth/hPositiveBachelorBaryonDCA"), casc.pt(), invmass, casc.bachBaryonDCAxyToPV());
            if (!doDCAV0ToPVCut && doPtDepCutStudy)
              histos.fill(HIST("PtDepCutStudyMCTruth/hPositiveDCAV0ToPV"), casc.pt(), invmass, TMath::Abs(casc.dcav0topv(casc.x(), casc.y(), casc.z())));
            if (!doV0RadiusCut && doPtDepCutStudy)
              histos.fill(HIST("PtDepCutStudyMCTruth/hPositiveV0Radius"), casc.pt(), invmass, casc.v0radius());
            if (!doCascadeRadiusCut && doPtDepCutStudy)
              histos.fill(HIST("PtDepCutStudyMCTruth/hPositiveCascadeRadius"), casc.pt(), invmass, casc.cascradius());
            if (!doDCAV0DauCut && doPtDepCutStudy)
              histos.fill(HIST("PtDepCutStudyMCTruth/hPositiveDCAV0Daughters"), casc.pt(), invmass, casc.dcaV0daughters());
            if (!doDCACascadeDauCut && doPtDepCutStudy)
              histos.fill(HIST("PtDepCutStudyMCTruth/hPositiveDCACascDaughters"), casc.pt(), invmass, casc.dcacascdaughters());
            if (!doV0CosPaCut && doPtDepCutStudy)
              histos.fill(HIST("PtDepCutStudyMCTruth/hPositiveV0pa"), casc.pt(), invmass, TMath::ACos(casc.v0cosPA(casc.x(), casc.y(), casc.z())));
            if (!doCascadeCosPaCut && doPtDepCutStudy)
              histos.fill(HIST("PtDepCutStudyMCTruth/hPositiveCascPA"), casc.pt(), invmass, TMath::ACos(casc.casccosPA(coll.posX(), coll.posY(), coll.posZ())));
            if (!doDCAdauToPVCut && doPtDepCutStudy) {
              histos.fill(HIST("PtDepCutStudyMCTruth/hPositiveDCABachelorToPV"), casc.pt(), invmass, casc.dcabachtopv());
              histos.fill(HIST("PtDepCutStudyMCTruth/hPositiveDCAMesonToPV"), casc.pt(), invmass, casc.dcapostopv());
              histos.fill(HIST("PtDepCutStudyMCTruth/hPositiveDCABaryonToPV"), casc.pt(), invmass, casc.dcanegtopv());
            }
            if (!doProperLifeTimeCut && doPtDepCutStudy)
              histos.fill(HIST("PtDepCutStudyMCTruth/hPositiveCascadeProperLifeTime"), casc.pt(), invmass, ctau);
            if ((isXi && casc.pdgCodeV0() == -3122 && casc.pdgCodePositive() == 211 && casc.pdgCodeNegative() == -2212 && casc.pdgCodeBachelor() == 211) || (!isXi && casc.pdgCodeV0() == -3122 && casc.pdgCodePositive() == 211 && casc.pdgCodeNegative() == -2212 && casc.pdgCodeBachelor() == 321))
              histos.fill(HIST("hPositiveCascadePtForEfficiency"), casc.pt(), invmass, coll.centFT0C());
          }
        }
      }
    }
  }

  void processReconstructedSimulation(aod::McCollision const& mcCollision, soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFV0As>> const& collisions, aod::McParticles const& mcParticles)
  {
    // this process function also checks if a given collision was reconstructed and checks explicitly for splitting, etc

    // identify best-of collision
    int biggestNContribs = -1;
    float bestCentrality = 100.5;
    for (auto& collision : collisions) {
      if (biggestNContribs < collision.numContrib()) {
        biggestNContribs = collision.numContrib();
        bestCentrality = collision.centFT0C();
      }
    }
    histos.fill(HIST("h2dNVerticesVsCentrality"), bestCentrality, collisions.size());

    for (auto& mcp : mcParticles) {
      if (TMath::Abs(mcp.y()) < 0.5 && mcp.isPhysicalPrimary()) {
        if (mcp.pdgCode() == 3312)
          histos.fill(HIST("hGenNegativeXi"), bestCentrality, mcp.pt());
        if (mcp.pdgCode() == -3312)
          histos.fill(HIST("hGenPositiveXi"), bestCentrality, mcp.pt());
        if (mcp.pdgCode() == 3334)
          histos.fill(HIST("hGenNegativeOmega"), bestCentrality, mcp.pt());
        if (mcp.pdgCode() == -3334)
          histos.fill(HIST("hGenPositiveOmega"), bestCentrality, mcp.pt());
      }
    }
  }

  PROCESS_SWITCH(derivedCascadeAnalysis, processCascades, "cascade analysis, run3 data ", true);
  PROCESS_SWITCH(derivedCascadeAnalysis, processCascadesMCrec, "cascade analysis, run3 rec MC", false);
  PROCESS_SWITCH(derivedCascadeAnalysis, processReconstructedSimulation, "cascade analysis, run3 gen MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<derivedCascadeAnalysis>(cfgc)};
}
