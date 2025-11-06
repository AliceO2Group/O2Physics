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
/// \brief strangeness in pb-pb QC task
///
/// In case of questions please write to:
/// \author Roman Nepeivoda (roman.nepeivoda@cern.ch)

#include "PWGLF/DataModel/QC/strangenessTablesQC.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct strangenessQCPP {
  // Histogram registries
  HistogramRegistry rGeneral{"generalInfo", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  HistogramRegistry rVzero{"vzero", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rCascade{"cascade", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rK0S{"k0S", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rLambda{"lambda", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rAntiLambda{"antiLambda", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rOmega{"omega", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rAntiOmega{"antiomega", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rXi{"xi", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rAntiXi{"antixi", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Configurable parameters for V0 selection
  Configurable<float> v0setting_dcav0dau{"v0setting_dcav0dau", 1, "DCA V0 Daughters"};
  Configurable<float> v0setting_dcapostopv{"v0setting_dcapostopv", 0.1, "DCA Pos To PV"};
  Configurable<float> v0setting_dcanegtopv{"v0setting_dcanegtopv", 0.1, "DCA Neg To PV"};
  Configurable<float> v0setting_cospa{"v0setting_cospa", 0.97, "V0 CosPA"}; // should be double!
  Configurable<float> v0setting_radius{"v0setting_radius", 1, "v0radius"};
  Configurable<float> v0setting_rapidity{"v0setting_rapidity", 0.5, "rapidity"};

  static constexpr float defaultLifetimeCuts[1][2] = {{25., 20.}};
  Configurable<LabeledArray<float>> lifetimecut{"lifetimecut", {defaultLifetimeCuts[0], 2, {"lifetimecutLambda", "lifetimecutK0S"}}, "lifetimecut"};

  // V0 PID configurables
  Configurable<float> NSigmaV0Pion{"NSigmaV0Pion", 6, "NSigmaV0Pion"};
  Configurable<float> NSigmaV0Proton{"NSigmaV0Proton", 6, "NSigmaV0Proton"};

  // Configurable parameters for cascade selection
  Configurable<float> cascadesetting_cospa{"cascadesetting_cospa", 0.98, "Casc CosPA"};        // should be double!
  Configurable<float> cascadesetting_v0cospa{"cascadesetting_v0cospa", 0.98, "Casc V0 CosPA"}; // should be double!
  Configurable<float> cascadesetting_dcacascdau{"cascadesetting_dcacascdau", 1.0, "DCA cascade daughters"};
  Configurable<float> cascadesetting_dcav0dau{"cascadesetting_dcav0dau", 1.0, "DCA Cascade's V0 Daughters"};
  Configurable<float> cascadesetting_dcabachtopv{"cascadesetting_dcabachtopv", 0.1, "DCA bachelor to PV"};
  Configurable<float> cascadesetting_dcapostopv{"cascadesetting_dcapostopv", 0.2, "DCA Casc. V0's pos to PV"};
  Configurable<float> cascadesetting_dcanegtopv{"cascadesetting_dcanegtopv", 0.2, "DCA Casc V0's neg to PV"};
  Configurable<float> cascadesetting_mindcav0topv{"cascadesetting_mindcav0topv", 0.01, "minimum V0 DCA to PV"};
  Configurable<float> cascadesetting_cascradius{"cascadesetting_cascradius", 1.0, "cascradius"};
  Configurable<float> cascadesetting_v0masswindow{"cascadesetting_v0masswindow", 0.01, "v0 mass window"};
  Configurable<float> cascadesetting_v0radius{"cascadesetting_v0radius", 0.9, "v0 radius"};
  Configurable<float> cascadesetting_rapidity{"cascadesetting_rapidity", 0.5, "rapidity"};

  // Cascade PID configurables
  Configurable<float> NSigmaCascPion{"NSigmaCascPion", 6, "NSigmaCascPion"};
  Configurable<float> NSigmaCascProton{"NSigmaCascProton", 6, "NSigmaCascProton"};
  Configurable<float> NSigmaCascKaon{"NSigmaCascKaon", 6, "NSigmaCascKaon"};

  // General axes configurables
  ConfigurableAxis binPt{"binPt", {100, 0.0f, 10.0f}, ""};
  ConfigurableAxis binPtsmall{"binPtsmall", {50, 0.0f, 10.0f}, ""};
  ConfigurableAxis binCosPA{"binCosPA", {200, 0.8f, 1.0f}, ""};
  ConfigurableAxis binEtaSmall{"binEtaSmall", {2, -1.0f, 1.0f}, ""};
  ConfigurableAxis binEta{"binEta", {100, -1.0f, 1.0f}, ""};
  ConfigurableAxis binPhi{"binPhi", {(int)TMath::Pi() * 10 / 2, 0.0f, 2. * TMath::Pi()}, ""};
  ConfigurableAxis binRadius{"binRadius", {100, 0.0f, 50.0f}, ""};
  ConfigurableAxis binRadiussmall{"binRadiussmall", {30, 0.0f, 30.0f}, ""};
  ConfigurableAxis binITSMapDaughters{"binITSMapDaughters", {8, -0.5f, 7.5f}, ""};

  // V0 axes configurables
  ConfigurableAxis binV0PA{"binV0PA", {1000, 0.f, 1.0f}, ""};
  ConfigurableAxis binV0Radius{"binV0Radius", {100, 0.0f, 10.0f}, ""};
  ConfigurableAxis binV0DecayLength{"binV0DecayLength", {100, 0.0f, 10.0f}, ""};
  ConfigurableAxis binV0DCANegToPV{"binV0DCANegToPV", {100, -1.0f, 1.0f}, ""};
  ConfigurableAxis binV0DCAPosToPV{"binV0DCAPosToPV", {100, -1.0f, 1.0f}, ""};
  ConfigurableAxis binV0DCAV0Dau{"binV0DCAV0Dau", {55, 0.0f, 2.20f}, ""};
  ConfigurableAxis binCtauK0s{"binCtauK0s", {65, 0.0f, 13.0f}, ""};
  ConfigurableAxis binCtauLambda{"binCtauLambda", {100, 0.0f, 40.0f}, ""};
  ConfigurableAxis binCtauAntiLambda{"binCtauAntiLambda", {100, 0.0f, 40.0f}, ""};
  ConfigurableAxis binDecayLengthK0s{"binDecayLengthK0s", {100, 0.0f, 40.0f}, ""};
  ConfigurableAxis binDecayLengthLambda{"binDecayLengthLambda", {100, 0.0f, 80.0f}, ""};
  ConfigurableAxis binDecayLengthAntiLambda{"binDecayLengthAntiLambda", {100, 0.0f, 80.0f}, ""};
  ConfigurableAxis binV0DCAV0ToPVK0S{"binV0DCAV0ToPVK0S", {250, 0.0f, 0.25f}, ""};
  ConfigurableAxis binV0DCAV0ToPVLambda{"binV0DCAV0ToPVLambda", {250, 0.0f, 0.25f}, ""};
  ConfigurableAxis binV0DCAV0ToPVAntiLambda{"binV0DCAV0ToPVAntiLambda", {250, 0.0f, 0.25f}, ""};
  ConfigurableAxis binInvMassK0S{"binInvMassK0S", {200, 0.4f, 0.6f}, ""};
  ConfigurableAxis binInvMassLambda{"binInvMassLambda", {200, 1.07f, 1.17f}, ""};
  ConfigurableAxis binInvMassAntiLambda{"binInvMassAntiLambda", {200, 1.07f, 1.17f}, ""};
  ConfigurableAxis binResponsePionFromLambda{"binResponsePionFromLambda", {200, -20.f, 20.f}, ""};
  ConfigurableAxis binResponseProtonFromLambda{"binResponseProtonFromLambda", {200, -20.f, 20.f}, ""};

  // Cascade axes configurables
  ConfigurableAxis binCascSign{"binCascSign", {2, -1.5f, 1.5f}, ""};
  ConfigurableAxis binInvMassCasc{"binInvMassCasc", {1000, 0.f, 1.0f}, ""};
  ConfigurableAxis binCascDecayLength{"binCascDecayLength", {100, 0.0f, 10.0f}, ""};
  ConfigurableAxis binCascRadius{"binCascRadius", {100, 0.0f, 10.0f}, ""};
  ConfigurableAxis binCascRapidity{"binCascRapidity", {200, -2.0f, 2.0f}, ""};
  ConfigurableAxis binCascCtau{"binCascCtau", {100, 0.0f, 100.0f}, ""};
  ConfigurableAxis binDcaCascDaughters{"binDcaCascDaughters", {110, 0.0f, 2.2f}, ""};
  ConfigurableAxis binDcaV0ToPV{"binDcaV0ToPV", {200, 0.0f, 2.f}, ""};
  ConfigurableAxis binDcaBachToPV{"binDcaBachToPV", {80, -0.2f, 0.2f}, ""};
  ConfigurableAxis binInvMassXi{"binInvMassXi", {80, 1.28f, 1.36f}, ""};
  ConfigurableAxis binInvMassOmega{"binInvMassOmega", {80, 1.63f, 1.71f}, ""};

  // PDG data base
  Service<o2::framework::O2DatabasePDG> pdgDB;

  void init(InitContext const&)
  {
    // General axes
    const AxisSpec axisPt{binPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisPtsmall{binPtsmall, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisCosPA{binCosPA, "Cos(PA)"};
    const AxisSpec axisEta{binEta, "Eta"};
    const AxisSpec axisPhi{binPhi, "Phi"};
    const AxisSpec axisEtaSmall{binEtaSmall, "Eta"};
    const AxisSpec axisRadius{binRadius, "Radius"};
    const AxisSpec axisRadiussmall{binRadiussmall, "Radius"};
    const AxisSpec axisITSMapDaughters{binITSMapDaughters, "ITS Map Daughters"};

    // V0 axes
    const AxisSpec axisV0PA{binV0PA, "Pointing Angle"};
    const AxisSpec axisV0Radius{binV0Radius, "V0 Radius (cm)"};
    const AxisSpec axisV0DecayLength{binV0DecayLength, "V0 Decay Length (cm)"};
    const AxisSpec axisV0DCANegToPV{binV0DCANegToPV, "V0 DCA Neg To PV (cm)"};
    const AxisSpec axisV0DCAPosToPV{binV0DCAPosToPV, "V0 DCA Pos To PV (cm)"};
    const AxisSpec axisV0DCAV0Dau{binV0DCAV0Dau, "V0 DCA V0 Daughters (cm)"};
    const AxisSpec axisCtauK0s{binCtauK0s, "K0s c#tau (cm)"};
    const AxisSpec axisCtauLambda{binCtauLambda, "Lambda c#tau (cm)"};
    const AxisSpec axisCtauAntiLambda{binCtauAntiLambda, "AntiLambda c#tau (cm)"};
    const AxisSpec axisDecayLengthK0s{binDecayLengthK0s, "Decay length K0s (cm)"};
    const AxisSpec axisDecayLengthLambda{binDecayLengthLambda, "Decay length Lambda (cm)"};
    const AxisSpec axisDecayLengthAntiLambda{binDecayLengthAntiLambda, "Decay length AntiLambda (cm)"};
    const AxisSpec axisV0DCAV0ToPVK0S{binV0DCAV0ToPVK0S, "DCAV0ToPV K0s"};
    const AxisSpec axisV0DCAV0ToPVLambda{binV0DCAV0ToPVLambda, "DCAV0ToPV Lambda"};
    const AxisSpec axisV0DCAV0ToPVAntiLambda{binV0DCAV0ToPVAntiLambda, "DCAV0ToPV AntiLambda"};
    const AxisSpec axisInvMassK0S{binInvMassK0S, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    const AxisSpec axisInvMassLambda{binInvMassLambda, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    const AxisSpec axisInvMassAntiLambda{binInvMassAntiLambda, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    const AxisSpec axisInvMassCasc{binInvMassCasc, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    const AxisSpec axisResponsePionFromLambda{binResponsePionFromLambda, "Response Pion From Lambda"};
    const AxisSpec axisResponseProtonFromLambda{binResponseProtonFromLambda, "Response Proton From Lambda"};

    // Cascade axes
    const AxisSpec axisSign{binCascSign, "Casc sign"};
    const AxisSpec axisCascDecayLength{binCascDecayLength, "Casc Decay Length (cm)"};
    const AxisSpec axisCascRadius{binCascRadius, "Casc Radius (cm)"};
    const AxisSpec axisCascRapidity{binCascRapidity, "y"};
    const AxisSpec axisCascCTau{binCascCtau, "y"};
    const AxisSpec axisDcaCascDaughters{binDcaCascDaughters, "DCA Cascade Daughters (cm)"};
    const AxisSpec axisDcaV0ToPV{binDcaV0ToPV, "DCA V0 to PV (cm)"};
    const AxisSpec axisDcaBachToPV{binDcaBachToPV, "DCA Bach to PV (cm)"};
    const AxisSpec axisInvMassXi{binInvMassXi, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    const AxisSpec axisInvMassOmega{binInvMassOmega, "#it{M}_{inv} [GeV/#it{c}^{2}]"};

    // Histograms
    // Candidate counter
    rGeneral.add("nCandidates", "nCandidates", {HistType::kTH1F, {{2, -0.5f, 1.5f}}});

    // V0 general histograms
    rVzero.add("CosPA", "CosPA", kTH1F, {axisCosPA});
    rVzero.add("V0Radius", "V0Radius", kTH1D, {axisV0Radius});
    rVzero.add("DecayLength", "DecayLength", kTH1F, {axisV0DecayLength});
    rVzero.add("V0DCANegToPV", "V0DCANegToPV", kTH1F, {axisV0DCANegToPV});
    rVzero.add("V0DCAPosToPV", "V0DCAPosToPV", kTH1F, {axisV0DCAPosToPV});
    rVzero.add("V0DCAV0Daughters", "V0DCAV0Daughters", kTH1F, {axisV0DCAV0Dau});

    // K0s histograms
    rK0S.add("CtauK0s", "CtauK0s", kTH1F, {axisCtauK0s});
    rK0S.add("DecayLengthK0s", "DecayLengthK0s", kTH1F, {axisDecayLengthK0s});
    rK0S.add("V0DCAV0ToPVK0S", "V0DCAV0ToPVK0S", kTH1F, {axisV0DCAV0ToPVK0S});
    rK0S.add("InvMassK0S", "InvMassK0S", kTH3F, {axisPt, axisInvMassK0S, axisEtaSmall});
    rK0S.add("InvMassK0SVsPtVsPA", "InvMassK0SVsPtVsPA", kTH3F, {axisPt, axisV0PA, axisInvMassK0S});
    rK0S.add("InvMassK0S_Radius", "InvMassK0S_Radius", kTH2F, {axisRadius, axisInvMassK0S});
    rK0S.add("InvMassK0S_EtaDaughters", "InvMassK0S_EtaDaughters", kTH3F, {axisEta, axisEta, axisInvMassK0S});
    rK0S.add("InvMassK0S_PhiDaughters", "InvMassK0S_PhiDaughters", kTH3F, {axisPhi, axisPhi, axisInvMassK0S});
    rK0S.add("InvMassK0S_ITSMapDaughters", "InvMassK0S_ITSMapDaughters", kTH3F, {axisITSMapDaughters, axisITSMapDaughters, axisInvMassK0S});
    rK0S.add("InvMassK0S_PtRadius", "InvMassK0S_PtRadius", kTH3F, {axisPtsmall, axisRadiussmall, axisInvMassK0S});

    // Lambda histograms
    rLambda.add("CtauLambda", "CtauLambda", kTH1F, {axisCtauLambda});
    rLambda.add("DecayLengthLambda", "DecayLengthLambda", kTH1F, {axisDecayLengthLambda});
    rLambda.add("V0DCAV0ToPVLambda", "V0DCAV0ToPVLambda", kTH1F, {axisV0DCAV0ToPVLambda});
    rLambda.add("InvMassLambda", "InvMassLambda", kTH3F, {axisPt, axisInvMassLambda, axisEtaSmall});
    rLambda.add("InvMassLambdaVsPtVsPA", "InvMassLambdaVsPtVsPA", kTH3F, {axisPt, axisV0PA, axisInvMassLambda});
    rLambda.add("InvMassLambda_Radius", "InvMassLambda_Radius", kTH2F, {axisRadius, axisInvMassLambda});
    rLambda.add("InvMassLambda_EtaDaughters", "InvMassLambda_EtaDaughters", kTH3F, {axisEta, axisEta, axisInvMassLambda});
    rLambda.add("InvMassLambda_Ctau", "InvMassLambda_Ctau", kTH2F, {axisCtauLambda, axisInvMassLambda});
    rLambda.add("InvMassLambda_PhiDaughters", "InvMassLambda_PhiDaughters", kTH3F, {axisPhi, axisPhi, axisInvMassLambda});
    rLambda.add("InvMassLambda_ITSMapDaughters", "InvMassLambda_ITSMapDaughters", kTH3F, {axisITSMapDaughters, axisITSMapDaughters, axisInvMassLambda});
    rLambda.add("InvMassLambda_PtRadius", "InvMassLambda_PtRadius", kTH3F, {axisPtsmall, axisRadiussmall, axisInvMassLambda});
    rLambda.add("ResponsePionFromLambda", "ResponsePionFromLambda", kTH2F, {axisPt, axisResponsePionFromLambda});
    rLambda.add("ResponseProtonFromLambda", "ResponseProtonFromLambda", kTH2F, {axisPt, axisResponseProtonFromLambda});

    // Anti-Lambda histograms
    rAntiLambda.add("CtauAntiLambda", "CtauAntiLambda", kTH1F, {axisCtauAntiLambda});
    rAntiLambda.add("DecayLengthAntiLambda", "DecayLengthAntiLambda", kTH1F, {axisDecayLengthAntiLambda});
    rAntiLambda.add("V0DCAV0ToPVAntiLambda", "V0DCAV0ToPVAntiLambda", kTH1F, {axisV0DCAV0ToPVAntiLambda});
    rAntiLambda.add("InvMassAntiLambda", "InvMassAntiLambda", kTH3F, {axisPt, axisInvMassAntiLambda, axisEtaSmall});
    rAntiLambda.add("InvMassAntiLambdaVsPtVsPA", "InvMassAntiLambdaVsPtVsPA", kTH3F, {axisPt, axisV0PA, axisInvMassAntiLambda});
    rAntiLambda.add("InvMassAntiLambda_Radius", "InvMassAntiLambda_Radius", kTH2F, {axisRadius, axisInvMassAntiLambda});
    rAntiLambda.add("InvMassAntiLambda_EtaDaughters", "InvMassAntiLambda_EtaDaughters", kTH3F, {axisEta, axisEta, axisInvMassAntiLambda});
    rAntiLambda.add("InvMassAntiLambda_Ctau", "InvMassAntiLambda_Ctau", kTH2F, {axisCtauAntiLambda, axisInvMassAntiLambda});
    rAntiLambda.add("InvMassAntiLambda_PhiDaughters", "InvMassAntiLambda_PhiDaughters", kTH3F, {axisPhi, axisPhi, axisInvMassAntiLambda});
    rAntiLambda.add("InvMassAntiLambda_ITSMapDaughters", "InvMassAntiLambda_ITSMapDaughters", kTH3F, {axisITSMapDaughters, axisITSMapDaughters, axisInvMassAntiLambda});
    rAntiLambda.add("InvMassAntiLambda_PtRadius", "InvMassAntiLambda_PtRadius", kTH3F, {axisPtsmall, axisRadiussmall, axisInvMassAntiLambda});

    // Cascade general histograms
    rCascade.add("V0Radius", "V0Radius", {HistType::kTH2D, {axisV0Radius, axisSign}});
    rCascade.add("CascCosPA", "CascCosPA", {HistType::kTH2F, {axisCosPA, axisSign}});
    rCascade.add("V0CosPA", "V0CosPA", {HistType::kTH2F, {axisCosPA, axisSign}});
    rCascade.add("CascDecayLength", "CascDecayLength", {HistType::kTH2F, {axisCascDecayLength, axisSign}});
    rCascade.add("CascRadius", "CascRadius", {HistType::kTH2F, {axisCascRadius, axisSign}});
    rCascade.add("CascyXi", "CascyXi", {HistType::kTH2F, {axisCascRapidity, axisSign}});
    rCascade.add("CascyOmega", "CascyOmega", {HistType::kTH2F, {axisCascRapidity, axisSign}});
    rCascade.add("CascCtauXi", "CascCtauXi", {HistType::kTH2F, {axisCascCTau, axisSign}});
    rCascade.add("CascCtauOmega", "CascCtauOmega", {HistType::kTH2F, {axisCascCTau, axisSign}});
    rCascade.add("V0Ctau", "V0Ctau", {HistType::kTH2F, {axisCascCTau, axisSign}});
    rCascade.add("CascPt", "CascPt", {HistType::kTH2F, {binPt, axisSign}});
    rCascade.add("DcaV0Daughters", "DcaV0Daughters", {HistType::kTH2F, {axisV0DCAV0Dau, axisSign}});
    rCascade.add("DcaCascDaughters", "DcaCascDaughters", {HistType::kTH2F, {axisDcaCascDaughters, axisSign}});
    rCascade.add("DcaV0ToPV", "DcaV0ToPV", {HistType::kTH2F, {axisDcaV0ToPV, axisSign}});
    rCascade.add("DcaBachToPV", "DcaBachToPV", {HistType::kTH2F, {axisDcaBachToPV, axisSign}});
    rCascade.add("DcaPosToPV", "DcaPosToPV", {HistType::kTH2F, {axisV0DCAPosToPV, axisSign}});
    rCascade.add("DcaNegToPV", "DcaNegToPV", {HistType::kTH2F, {axisV0DCANegToPV, axisSign}});
    rCascade.add("InvMassLambdaDaughter", "InvMassLambdaDaughter", {HistType::kTH2F, {axisInvMassLambda, axisSign}});
    // rCascade.add("V0CosPAToXi", "V0CosPAToXi", {HistType::kTH2F, {{100, 0.9f, 1.0f}, axisSign}});

    // Xi histograms
    rXi.add("InvMassXiMinus", "InvMassXiMinus", {HistType::kTH3F, {axisPt, axisInvMassXi, axisEta}});
    rXi.add("InvMassXiMinus_Radius", "InvMassXiMinus_Radius", {HistType::kTH2F, {axisPt, axisInvMassXi}});

    // Anti-Xi histograms
    rAntiXi.add("InvMassXiPlus", "InvMassXiPlus", {HistType::kTH3F, {axisPt, axisInvMassXi, axisEta}});
    rAntiXi.add("InvMassXiPlus_Radius", "InvMassXiPlus_Radius", {HistType::kTH2F, {axisPt, axisInvMassXi}});

    // Omega histograms
    rOmega.add("InvMassOmegaMinus", "InvMassOmegaMinus", {HistType::kTH3F, {axisPt, axisInvMassOmega, axisEta}});

    // Anti-Omega histograms
    rAntiOmega.add("InvMassOmegaPlus", "InvMassOmegaPlus", {HistType::kTH3F, {axisPt, axisInvMassOmega, axisEta}});

    // Cut summary
    rGeneral.add("selectionSummary", "selectionSummary", HistType::kTH1F, {{18, -0.5, 17.5}});
    TString CutLabelSummary[18] = {"v0_dcav0dau", "v0_dcapostopv", "v0_dcanegtopv", "v0_cospa", "v0_radius", "v0_rapidity",
                                   "casc_cospa", "casc_v0cospa", "casc_dcacascdau", "casc_dcav0dau", "casc_dcabachtopv", "casc_dcapostopv", "casc_dcanegtopv", "casc_mindcav0topv", "casc_cascradius", "casc_v0masswindow", "casc_v0radius", "casc_rapidity"};
    for (Int_t n = 1; n <= rGeneral.get<TH1>(HIST("selectionSummary"))->GetNbinsX(); n++) {
      rGeneral.get<TH1>(HIST("selectionSummary"))->GetXaxis()->SetBinLabel(n, CutLabelSummary[n - 1]);
    }
    rGeneral.get<TH1>(HIST("selectionSummary"))->SetBinContent(1, v0setting_dcav0dau);
    rGeneral.get<TH1>(HIST("selectionSummary"))->SetBinContent(2, v0setting_dcapostopv);
    rGeneral.get<TH1>(HIST("selectionSummary"))->SetBinContent(3, v0setting_dcanegtopv);
    rGeneral.get<TH1>(HIST("selectionSummary"))->SetBinContent(4, v0setting_cospa);
    rGeneral.get<TH1>(HIST("selectionSummary"))->SetBinContent(5, v0setting_radius);
    rGeneral.get<TH1>(HIST("selectionSummary"))->SetBinContent(6, v0setting_rapidity);
    rGeneral.get<TH1>(HIST("selectionSummary"))->SetBinContent(7, cascadesetting_cospa);
    rGeneral.get<TH1>(HIST("selectionSummary"))->SetBinContent(8, cascadesetting_v0cospa);
    rGeneral.get<TH1>(HIST("selectionSummary"))->SetBinContent(9, cascadesetting_dcacascdau);
    rGeneral.get<TH1>(HIST("selectionSummary"))->SetBinContent(10, cascadesetting_dcav0dau);
    rGeneral.get<TH1>(HIST("selectionSummary"))->SetBinContent(11, cascadesetting_dcabachtopv);
    rGeneral.get<TH1>(HIST("selectionSummary"))->SetBinContent(12, cascadesetting_dcapostopv);
    rGeneral.get<TH1>(HIST("selectionSummary"))->SetBinContent(13, cascadesetting_dcanegtopv);
    rGeneral.get<TH1>(HIST("selectionSummary"))->SetBinContent(14, cascadesetting_mindcav0topv);
    rGeneral.get<TH1>(HIST("selectionSummary"))->SetBinContent(15, cascadesetting_cascradius);
    rGeneral.get<TH1>(HIST("selectionSummary"))->SetBinContent(16, cascadesetting_v0masswindow);
    rGeneral.get<TH1>(HIST("selectionSummary"))->SetBinContent(17, cascadesetting_v0radius);
    rGeneral.get<TH1>(HIST("selectionSummary"))->SetBinContent(18, cascadesetting_rapidity);
  }

  template <typename TV0>
  void fillV0Histograms(TV0 const& v0)
  {
    rVzero.fill(HIST("CosPA"), v0.v0cosPA());
    rVzero.fill(HIST("V0Radius"), v0.v0radius());
    rVzero.fill(HIST("V0DCANegToPV"), v0.dcanegtopv());
    rVzero.fill(HIST("V0DCAPosToPV"), v0.dcapostopv());
    rVzero.fill(HIST("V0DCAV0Daughters"), v0.dcaV0daughters());
    rVzero.fill(HIST("DecayLength"), v0.decayLength());
  }

  template <typename TCascade>
  void fillCascadeHistograms(TCascade const& casc)
  {
    rCascade.fill(HIST("V0Radius"), casc.v0radius(), casc.sign());
    rCascade.fill(HIST("CascCosPA"), casc.casccosPA(), casc.sign());
    rCascade.fill(HIST("V0CosPA"), casc.v0cosPA(), casc.sign());
    rCascade.fill(HIST("CascRadius"), casc.cascradius(), casc.sign());
    rCascade.fill(HIST("CascDecayLength"), casc.decayLength(), casc.sign());
    rCascade.fill(HIST("CascPt"), casc.pt(), casc.sign());
    rCascade.fill(HIST("DcaV0Daughters"), casc.dcaV0daughters(), casc.sign());
    rCascade.fill(HIST("DcaCascDaughters"), casc.dcacascdaughters(), casc.sign());
    rCascade.fill(HIST("DcaV0ToPV"), casc.dcav0topv(), casc.sign());
    rCascade.fill(HIST("DcaBachToPV"), casc.dcabachtopv(), casc.sign());
    rCascade.fill(HIST("DcaPosToPV"), casc.dcapostopv(), casc.sign());
    rCascade.fill(HIST("DcaNegToPV"), casc.dcanegtopv(), casc.sign());
    rCascade.fill(HIST("InvMassLambdaDaughter"), casc.mLambda(), casc.sign());
    rCascade.fill(HIST("V0Ctau"), casc.lifetimeV0(), casc.sign());
  }

  // Filters on V0s
  Filter preFilterV0 = (nabs(aod::vZerosQC::dcapostopv) > v0setting_dcapostopv &&
                        nabs(aod::vZerosQC::dcanegtopv) > v0setting_dcanegtopv &&
                        aod::vZerosQC::dcaV0daughters < v0setting_dcav0dau &&
                        aod::vZerosQC::v0cosPA > v0setting_cospa &&
                        aod::vZerosQC::v0radius > v0setting_radius);

  // Filters on Cascades
  Filter preFilterCascades = (nabs(aod::cascadesQC::dcabachtopv) > cascadesetting_dcabachtopv &&
                              aod::cascadesQC::dcaV0daughters < cascadesetting_dcav0dau &&
                              nabs(aod::cascadesQC::dcapostopv) > cascadesetting_dcapostopv &&
                              nabs(aod::cascadesQC::dcanegtopv) > cascadesetting_dcanegtopv &&
                              aod::cascadesQC::dcacascdaughters < cascadesetting_dcacascdau &&
                              aod::cascadesQC::dcaV0daughters < cascadesetting_dcav0dau &&
                              aod::cascadesQC::v0radius > cascadesetting_v0radius &&
                              aod::cascadesQC::cascradius > cascadesetting_cascradius &&
                              nabs(aod::cascadesQC::dcav0topv) > cascadesetting_mindcav0topv &&
                              aod::cascadesQC::v0cosPA > cascadesetting_v0cospa &&
                              aod::cascadesQC::casccosPA > cascadesetting_cospa);

  void processV0(soa::Filtered<aod::VZerosQC> const& v0s)
  {
    for (const auto& v0 : v0s) {
      // Fill the candidate counter
      rGeneral.fill(HIST("nCandidates"), 0);

      // K0Short
      if (TMath::Abs(v0.yK0Short()) < v0setting_rapidity &&
          v0.lifetimeK0s() < lifetimecut->get("lifetimecutK0S") &&
          TMath::Abs(v0.posNSigmaV0Pion()) < NSigmaV0Pion && TMath::Abs(v0.negNSigmaV0Pion()) < NSigmaV0Pion) {
        fillV0Histograms(v0);
        rK0S.fill(HIST("DecayLengthK0s"), v0.decayLength());
        rK0S.fill(HIST("CtauK0s"), v0.lifetimeK0s());
        rK0S.fill(HIST("InvMassK0S"), v0.pt(), v0.mK0Short(), v0.eta());
        rK0S.fill(HIST("InvMassK0SVsPtVsPA"), v0.pt(), TMath::ACos(v0.v0cosPA()), v0.mK0Short());
        rK0S.fill(HIST("V0DCAV0ToPVK0S"), v0.dcav0topv());
        rK0S.fill(HIST("InvMassK0S_Radius"), v0.v0radius(), v0.mK0Short());
        rK0S.fill(HIST("InvMassK0S_PtRadius"), v0.pt(), v0.v0radius(), v0.mK0Short());
        rK0S.fill(HIST("InvMassK0S_EtaDaughters"), v0.poseta(), v0.negeta(), v0.mK0Short());
        rK0S.fill(HIST("InvMassK0S_PhiDaughters"), v0.posphi(), v0.negphi(), v0.mK0Short());
        rK0S.fill(HIST("InvMassK0S_ITSMapDaughters"), v0.posITSNhits(), v0.negITSNhits(), v0.mK0Short());
      }

      // Lambda
      if (TMath::Abs(v0.yLambda()) < v0setting_rapidity &&
          v0.lifetimeLambda() < lifetimecut->get("lifetimecutLambda") &&
          TMath::Abs(v0.posNSigmaV0Proton()) < NSigmaV0Proton && TMath::Abs(v0.negNSigmaV0Pion()) < NSigmaV0Pion) {
        fillV0Histograms(v0);
        rLambda.fill(HIST("DecayLengthLambda"), v0.decayLength());
        rLambda.fill(HIST("CtauLambda"), v0.lifetimeLambda());
        rLambda.fill(HIST("InvMassLambda"), v0.pt(), v0.mLambda(), v0.eta());
        rLambda.fill(HIST("InvMassLambdaVsPtVsPA"), v0.pt(), TMath::ACos(v0.v0cosPA()), v0.mLambda());
        rLambda.fill(HIST("V0DCAV0ToPVLambda"), v0.dcav0topv());
        rLambda.fill(HIST("InvMassLambda_Radius"), v0.v0radius(), v0.mLambda());
        rLambda.fill(HIST("InvMassLambda_PtRadius"), v0.pt(), v0.v0radius(), v0.mLambda());
        rLambda.fill(HIST("InvMassLambda_EtaDaughters"), v0.poseta(), v0.negeta(), v0.mLambda());
        rLambda.fill(HIST("InvMassLambda_PhiDaughters"), v0.posphi(), v0.negphi(), v0.mLambda());
        rLambda.fill(HIST("InvMassLambda_Ctau"), v0.lifetimeLambda(), v0.mLambda());
        rLambda.fill(HIST("InvMassLambda_ITSMapDaughters"), v0.posITSNhits(), v0.negITSNhits(), v0.mLambda());
        if (v0.v0cosPA() > 0.999 && v0.dcaV0daughters() < 1 && TMath::Abs(v0.mK0Short() - pdgDB->Mass(310)) > 0.012 && TMath::Abs(v0.mAntiLambda() - o2::constants::physics::MassLambda0) > 0.08 && TMath::Abs(v0.mLambda() - o2::constants::physics::MassLambda0) < 0.002) {
          rLambda.fill(HIST("ResponsePionFromLambda"), v0.pt(), v0.negNSigmaV0Pion());
          rLambda.fill(HIST("ResponseProtonFromLambda"), v0.pt(), v0.posNSigmaV0Proton());
        }
      }

      // AntiLambda
      if (TMath::Abs(v0.yLambda()) < v0setting_rapidity &&
          v0.lifetimeLambda() < lifetimecut->get("lifetimecutLambda") &&
          TMath::Abs(v0.posNSigmaV0Pion()) < NSigmaV0Pion && TMath::Abs(v0.negNSigmaV0Proton()) < NSigmaV0Proton) {
        fillV0Histograms(v0);
        rAntiLambda.fill(HIST("DecayLengthAntiLambda"), v0.decayLength());
        rAntiLambda.fill(HIST("CtauAntiLambda"), v0.lifetimeLambda());
        rAntiLambda.fill(HIST("InvMassAntiLambda"), v0.pt(), v0.mAntiLambda(), v0.eta());
        rAntiLambda.fill(HIST("InvMassAntiLambdaVsPtVsPA"), v0.pt(), TMath::ACos(v0.v0cosPA()), v0.mAntiLambda());
        rAntiLambda.fill(HIST("V0DCAV0ToPVAntiLambda"), v0.dcav0topv());
        rAntiLambda.fill(HIST("InvMassAntiLambda_Radius"), v0.v0radius(), v0.mAntiLambda());
        rAntiLambda.fill(HIST("InvMassAntiLambda_PtRadius"), v0.pt(), v0.v0radius(), v0.mAntiLambda());
        rAntiLambda.fill(HIST("InvMassAntiLambda_EtaDaughters"), v0.poseta(), v0.negeta(), v0.mAntiLambda());
        rAntiLambda.fill(HIST("InvMassAntiLambda_PhiDaughters"), v0.posphi(), v0.negphi(), v0.mAntiLambda());
        rAntiLambda.fill(HIST("InvMassAntiLambda_Ctau"), v0.lifetimeLambda(), v0.mAntiLambda());
        rAntiLambda.fill(HIST("InvMassAntiLambda_ITSMapDaughters"), v0.posITSNhits(), v0.negITSNhits(), v0.mAntiLambda());
      }
    }
  }
  PROCESS_SWITCH(strangenessQCPP, processV0, "Process V0 candidates", true);

  void processCascade(soa::Filtered<aod::CascadesQC> const& cascades)
  {
    for (const auto& casc : cascades) {
      // Fill the candidate counter
      rGeneral.fill(HIST("nCandidates"), 1);

      if (casc.sign() < 0) {
        // Check lambda daughters` PID
        if (TMath::Abs(casc.posNSigmaV0Proton()) < NSigmaCascProton && TMath::Abs(casc.negNSigmaV0Pion()) < NSigmaCascPion) {
          // Xi
          if (TMath::Abs(casc.yXi()) < cascadesetting_rapidity && TMath::Abs(casc.bachNSigmaV0Pion()) < NSigmaCascPion) {
            fillCascadeHistograms(casc);
            rXi.fill(HIST("InvMassXiMinus"), casc.pt(), casc.mXi(), casc.eta());
            rXi.fill(HIST("InvMassXiMinus_Radius"), casc.cascradius(), casc.mXi());
            rCascade.fill(HIST("CascyXi"), casc.yXi(), casc.sign());
            rCascade.fill(HIST("CascCtauXi"), casc.lifetimeXi(), casc.sign());
          }
          // Omega
          if (TMath::Abs(casc.yOmega()) < cascadesetting_rapidity && TMath::Abs(casc.bachNSigmaV0Kaon()) < NSigmaCascKaon) {
            fillCascadeHistograms(casc);
            rOmega.fill(HIST("InvMassOmegaMinus"), casc.pt(), casc.mOmega(), casc.eta());
            rCascade.fill(HIST("CascCtauOmega"), casc.lifetimeOmega(), casc.sign());
            rCascade.fill(HIST("CascyOmega"), casc.yOmega(), casc.sign());
          }
        }
      } else {
        // Check anti-lambda daughters` PID
        if (TMath::Abs(casc.posNSigmaV0Pion()) < NSigmaCascPion && TMath::Abs(casc.negNSigmaV0Proton()) < NSigmaCascProton) {
          // Anti-Xi
          if (TMath::Abs(casc.yXi()) < cascadesetting_rapidity && TMath::Abs(casc.bachNSigmaV0Pion()) < NSigmaCascPion) {
            fillCascadeHistograms(casc);
            rAntiXi.fill(HIST("InvMassXiPlus"), casc.pt(), casc.mXi(), casc.eta());
            rAntiXi.fill(HIST("InvMassXiPlus_Radius"), casc.cascradius(), casc.mXi());
            rCascade.fill(HIST("CascyXi"), casc.yXi(), casc.sign());
            rCascade.fill(HIST("CascCtauXi"), casc.lifetimeXi(), casc.sign());
          }
          // Anti-Omega
          if (TMath::Abs(casc.yOmega()) < cascadesetting_rapidity && TMath::Abs(casc.bachNSigmaV0Kaon()) < NSigmaCascKaon) {
            fillCascadeHistograms(casc);
            rAntiOmega.fill(HIST("InvMassOmegaPlus"), casc.pt(), casc.mOmega(), casc.eta());
            rCascade.fill(HIST("CascCtauOmega"), casc.lifetimeOmega(), casc.sign());
            rCascade.fill(HIST("CascyOmega"), casc.yOmega(), casc.sign());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(strangenessQCPP, processCascade, "Process cascade candidates", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<strangenessQCPP>(cfgc, TaskName{"lf-strangenessqcpp"})};
}
