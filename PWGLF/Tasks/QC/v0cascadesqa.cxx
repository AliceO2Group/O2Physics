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
/// \brief QA task for V0s and Cascades
///
/// In case of questions please write to:
/// \author Aimeric Landou (aimeric.landou@cern.ch)
/// \author Chiara De Martin (chiara.de.martin@cern.ch)
/// \author Francesca Ercolessi (francesca.ercolessi@cern.ch)

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <cmath>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra>;
using MyTracksMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::McTrackLabels>;
using DaughterTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCPr, aod::pidTPCKa>;

struct v0cascadesQA {

  static constexpr float defaultLifetimeCuts[1][2] = {{25., 20.}};

  // configurable event properties
  Configurable<bool> isMC{"isMC", false, "does the data have MC info"};
  Configurable<bool> sel8{"sel8", 0, "Apply sel8 event selection"};
  Configurable<bool> doextraanalysis{"doextraanalysis", 0, "Add extra histograms"};

  // configurable track properties
  Configurable<bool> checkDauTPC{"checkDauTPC", false, "check if daughter tracks have TPC match"};

  // configurable binning of histograms
  ConfigurableAxis binPt{"binPt", {100, 0.0f, 10.0f}, ""};
  ConfigurableAxis binPtsmall{"binPtsmall", {50, 0.0f, 10.0f}, ""};
  ConfigurableAxis binV0CosPA{"binV0CosPA", {200, 0.8f, 1.0f}, ""};
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
  ConfigurableAxis binEtaFlag{"binEtaFlag", {3, -1.5f, 1.5f}, ""};
  ConfigurableAxis binEta{"binEta", {100, -1.0f, 1.0f}, ""};
  ConfigurableAxis binPhi{"binPhi", {static_cast<int>(TMath::Pi()) * 10 / 2, 0.0f, 2. * static_cast<int>(TMath::Pi())}, ""};
  ConfigurableAxis binRadius{"binRadius", {100, 0.0f, 50.0f}, ""};
  ConfigurableAxis binRadiussmall{"binRadiussmall", {30, 0.0f, 30.0f}, ""};
  ConfigurableAxis binITSMapDaughters{"binITSMapDaughters", {8, -0.5f, 7.5f}, ""};
  ConfigurableAxis binInvMassCasc{"binInvMassCasc", {1000, 0.f, 1.0f}, ""};

  // configurables for V0s
  Configurable<float> V0_rapidity{"V0_rapidity", 0.5, "rapidity"};
  Configurable<double> V0_cosPA{"V0_cosPA", 0.995, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<float> V0_dcav0dau{"V0_dcav0dau", 1, "DCA V0 Daughters"};
  Configurable<float> V0_dcapostopv{"V0_dcapostopv", 0.1, "DCA Pos To PV"};
  Configurable<float> V0_dcanegtopv{"V0_dcanegtopv", 0.1, "DCA Neg To PV"};
  Configurable<float> V0_radius{"V0_radius", 5, "v0radius"};
  Configurable<float> NSigmaV0Pion{"NSigmaV0Pion", 6, "NSigmaV0Pion"};
  Configurable<float> NSigmaV0Proton{"NSigmaV0Proton", 6, "NSigmaV0Proton"};
  Configurable<LabeledArray<float>> lifetimecut{"lifetimecut", {defaultLifetimeCuts[0], 2, {"lifetimecutLambda", "lifetimecutK0S"}}, "lifetimecut"};

  // configurables for Cascades
  Configurable<float> Casc_rapidity{"Casc_rapidity", 0.5, "rapidity"};
  Configurable<double> Casc_v0cospa{"Casc_V0cospa", 0.98, "V0 CosPA"};               // double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<double> Casc_casccospa{"Casc_casccospa", 0.98, "Cascade CosPA"};      // AliAnalysisTaskStrAODqa: 0.9992      //double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<float> Casc_dcav0dau{"Casc_dcav0dau", 1.0, "DCA V0 Daughters"};       // AliAnalysisTaskStrAODqa: 1. different scale
  Configurable<float> Casc_dcacascdau{"Casc_dcacascdau", 0.6, "DCA Casc Daughters"}; // AliAnalysisTaskStrAODqa: 0.3 different scale
  Configurable<float> Casc_dcav0topv{"Casc_dcav0topv", 0.1, "DCA Pos To PV"};        // AliAnalysisTaskStrAODqa: 0.15 different scale
  Configurable<float> Casc_dcabachtopv{"Casc_dcabachtopv", .1, "DCA Bach To PV"};    // AliAnalysisTaskStrAODqa: 0.17 different scale
  Configurable<float> Casc_dcapostopv{"Casc_dcapostopv", 0.1, "DCA V0 To PV"};       // AliAnalysisTaskStrAODqa:    if( fCasc_charge>0 &&(fCasc_DcaPosToPV < 0.3 || fCasc_DcaNegToPV < 0.11)) return kFALSE;  different scale
  Configurable<float> Casc_dcanegtopv{"Casc_dcanegtopv", 0.1, "DCA Neg To PV"};      // AliAnalysisTaskStrAODqa:    if( fCasc_charge<0 &&(fCasc_DcaPosToPV < 0.11 || fCasc_DcaNegToPV < 0.3)) return kFALSE;  different scale
  Configurable<float> Casc_v0radius{"Casc_v0radius", 0.9, "v0 radius"};              // AliAnalysisTaskStrAODqa: 5.
  Configurable<float> Casc_cascradius{"Casc_cascradius", 1.0, "cascade radius"};     // AliAnalysisTaskStrAODqa: 1.
  Configurable<float> NSigmaCascPion{"NSigmaCascPion", 6, "NSigmaCascPion"};
  Configurable<float> NSigmaCascProton{"NSigmaCascProton", 6, "NSigmaCascProton"};
  Configurable<float> NSigmaCascKaon{"NSigmaCascKaon", 6, "NSigmaCascKaon"};

  OutputObj<TH1F> V0SelectionSummary{TH1F("V0SelectionSummary", "V0SelectionSummary; Selections; Cut", 10, 0., 10.)};
  HistogramRegistry histos_eve{"histos-eve", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry histos_V0{"histos-V0", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry histos_Casc{"histos-Casc", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};

  void init(InitContext const&)
  {
    const AxisSpec axisPt{binPt, "p_{T} (GeV/c)"};
    const AxisSpec axisPtsmall{binPtsmall, "p_{T} (GeV/c)"};
    const AxisSpec axisV0CosPA{binV0CosPA, "V0 Cos(PA)"};
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
    const AxisSpec axisInvMassK0S{binInvMassK0S, "InvMass K0s"};
    const AxisSpec axisInvMassLambda{binInvMassLambda, "InvMass Lambda"};
    const AxisSpec axisInvMassAntiLambda{binInvMassAntiLambda, "InvMass AntiLambda"};
    const AxisSpec axisInvMassCasc{binInvMassCasc, "InvMass Cascades"};
    const AxisSpec axisResponsePionFromLambda{binResponsePionFromLambda, "Response Pion From Lambda"};
    const AxisSpec axisResponseProtonFromLambda{binResponseProtonFromLambda, "Response Proton From Lambda"};
    const AxisSpec axisEta{binEta, "Eta"};
    const AxisSpec axisPhi{binPhi, "Phi"};
    const AxisSpec axisEtaFlag{binEtaFlag, "Eta"};
    const AxisSpec axisRadius{binRadius, "Radius"};
    const AxisSpec axisRadiussmall{binRadiussmall, "Radius"};
    const AxisSpec axisITSMapDaughters{binITSMapDaughters, "ITS Map Daughters"};

    histos_eve.add("GeneratedParticles", "GeneratedParticles", {HistType::kTH3F, {{14, 0.0f, 14.0f}, {100, 0, 10}, {100, 0.f, 50.f}}});
    histos_eve.add("hEventCounter", "hEventCounter", {HistType::kTH1F, {{2, 0.0f, 2.0f}}});
    histos_eve.add("hEventCounterMC", "hEventCounterMC", {HistType::kTH1F, {{2, 0.0f, 2.0f}}});
    histos_eve.add("GenK0sPtVsEta", "GenK0sPtVsEta", {HistType::kTH2F, {axisPt, axisEta}});
    histos_eve.add("GenLambdaPtVsEta", "GenLambdaPtVsEta", {HistType::kTH2F, {axisPt, axisEta}});
    histos_eve.add("GenAntiLambdaPtVsEta", "GenAntiLambdaPtVsEta", {HistType::kTH2F, {axisPt, axisEta}});
    histos_eve.add("GenXiMinusPtVsEta", "GenXiMinusPtVsEta", {HistType::kTH2F, {axisPt, axisEta}});
    histos_eve.add("GenXiPlusPtVsEta", "GenXiPlusPtVsEta", {HistType::kTH2F, {axisPt, axisEta}});
    histos_eve.add("GenOmegaMinusPtVsEta", "GenOmegaMinusPtVsEta", {HistType::kTH2F, {axisPt, axisEta}});
    histos_eve.add("GenOmegaPlusPtVsEta", "GenOmegaPlusPtVsEta", {HistType::kTH2F, {axisPt, axisEta}});

    histos_V0.add("CosPA", "CosPA", HistType::kTH1F, {axisV0CosPA});
    histos_V0.add("V0Radius", "V0Radius", HistType::kTH1F, {axisV0Radius});
    histos_V0.add("DecayLength", "DecayLength", HistType::kTH1F, {axisV0DecayLength});
    histos_V0.add("V0DCANegToPV", "V0DCANegToPV", HistType::kTH1F, {axisV0DCANegToPV});
    histos_V0.add("V0DCAPosToPV", "V0DCAPosToPV", HistType::kTH1F, {axisV0DCAPosToPV});
    histos_V0.add("V0DCAV0Daughters", "V0DCAV0Daughters", HistType::kTH1F, {axisV0DCAV0Dau});
    histos_V0.add("CtauK0s", "CtauK0s", HistType::kTH1F, {axisCtauK0s});
    histos_V0.add("CtauLambda", "CtauLambda", HistType::kTH1F, {axisCtauLambda});
    histos_V0.add("CtauAntiLambda", "CtauAntiLambda", HistType::kTH1F, {axisCtauAntiLambda});
    histos_V0.add("DecayLengthK0s", "DecayLengthK0s", HistType::kTH1F, {axisDecayLengthK0s});
    histos_V0.add("DecayLengthLambda", "DecayLengthLambda", HistType::kTH1F, {axisDecayLengthLambda});
    histos_V0.add("DecayLengthAntiLambda", "DecayLengthAntiLambda", HistType::kTH1F, {axisDecayLengthAntiLambda});
    histos_V0.add("V0DCAV0ToPVK0S", "V0DCAV0ToPVK0S", HistType::kTH1F, {axisV0DCAV0ToPVK0S});
    histos_V0.add("V0DCAV0ToPVLambda", "V0DCAV0ToPVLambda", HistType::kTH1F, {axisV0DCAV0ToPVLambda});
    histos_V0.add("V0DCAV0ToPVAntiLambda", "V0DCAV0ToPVAntiLambda", HistType::kTH1F, {axisV0DCAV0ToPVAntiLambda});
    histos_V0.add("InvMassK0S", "InvMassK0S", HistType::kTH3F, {axisPt, axisInvMassK0S, axisEtaFlag});
    histos_V0.add("InvMassLambda", "InvMassLambda", HistType::kTH3F, {axisPt, axisInvMassLambda, axisEtaFlag});
    histos_V0.add("InvMassAntiLambda", "InvMassAntiLambda", HistType::kTH3F, {axisPt, axisInvMassAntiLambda, axisEtaFlag});
    histos_V0.add("ResponsePionFromLambda", "ResponsePionFromLambda", HistType::kTH2F, {axisPt, axisResponsePionFromLambda});
    histos_V0.add("ResponseProtonFromLambda", "ResponseProtonFromLambda", HistType::kTH2F, {axisPt, axisResponseProtonFromLambda});
    histos_V0.add("InvMassK0SVsPtVsPA", "InvMassK0SVsPtVsPA", HistType::kTH3F, {axisPt, axisV0PA, axisInvMassK0S});
    histos_V0.add("InvMassK0STrue", "InvMassK0STrue", {HistType::kTH3F, {{100, 0.0f, 10.0f}, {100, 0.f, 50.f}, {200, 0.4f, 0.6f}}});
    histos_V0.add("InvMassLambdaTrue", "InvMassLambdaTrue", {HistType::kTH3F, {{100, 0.0f, 10.0f}, {100, 0.f, 50.f}, {200, 1.07f, 1.17f}}});
    histos_V0.add("InvMassAntiLambdaTrue", "InvMassAntiLambdaTrue", {HistType::kTH3F, {{100, 0.0f, 10.0f}, {100, 0.f, 50.f}, {200, 1.07f, 1.17f}}});
    if (doextraanalysis) {
      histos_V0.add("InvMassK0S_Radius", "InvMassK0S_Radius", HistType::kTH2F, {axisRadius, axisInvMassK0S});
      histos_V0.add("InvMassLambda_Radius", "InvMassLambda_Radius", HistType::kTH2F, {axisRadius, axisInvMassLambda});
      histos_V0.add("InvMassAntiLambda_Radius", "InvMassAntiLambda_Radius", HistType::kTH2F, {axisRadius, axisInvMassAntiLambda});
      histos_V0.add("InvMassK0S_EtaDaughters", "InvMassK0S_EtaDaughters", HistType::kTH3F, {axisEta, axisEta, axisInvMassK0S});
      histos_V0.add("InvMassLambda_EtaDaughters", "InvMassLambda_EtaDaughters", HistType::kTH3F, {axisEta, axisEta, axisInvMassLambda});
      histos_V0.add("InvMassAntiLambda_EtaDaughters", "InvMassAntiLambda_EtaDaughters", HistType::kTH3F, {axisEta, axisEta, axisInvMassAntiLambda});
      histos_V0.add("InvMassK0S_Ctau", "InvMassK0S_Ctau", HistType::kTH2F, {axisCtauK0s, axisInvMassK0S});
      histos_V0.add("InvMassLambda_Ctau", "InvMassLambda_Ctau", HistType::kTH2F, {axisCtauLambda, axisInvMassLambda});
      histos_V0.add("InvMassAntiLambda_Ctau", "InvMassAntiLambda_Ctau", HistType::kTH2F, {axisCtauAntiLambda, axisInvMassAntiLambda});
      histos_V0.add("InvMassK0S_PhiDaughters", "InvMassK0S_PhiDaughters", HistType::kTH3F, {axisPhi, axisPhi, axisInvMassK0S});
      histos_V0.add("InvMassLambda_PhiDaughters", "InvMassLambda_PhiDaughters", HistType::kTH3F, {axisPhi, axisPhi, axisInvMassLambda});
      histos_V0.add("InvMassAntiLambda_PhiDaughters", "InvMassAntiLambda_PhiDaughters", HistType::kTH3F, {axisPhi, axisPhi, axisInvMassAntiLambda});
      histos_V0.add("InvMassK0S_ITSMapDaughters", "InvMassK0S_ITSMapDaughters", HistType::kTH3F, {axisITSMapDaughters, axisITSMapDaughters, axisInvMassK0S});
      histos_V0.add("InvMassLambda_ITSMapDaughters", "InvMassLambda_ITSMapDaughters", HistType::kTH3F, {axisITSMapDaughters, axisITSMapDaughters, axisInvMassLambda});
      histos_V0.add("InvMassAntiLambda_ITSMapDaughters", "InvMassAntiLambda_ITSMapDaughters", HistType::kTH3F, {axisITSMapDaughters, axisITSMapDaughters, axisInvMassAntiLambda});
      histos_V0.add("InvMassK0S_PtRadius", "InvMassK0S_PtRadius", HistType::kTH3F, {axisPtsmall, axisRadiussmall, axisInvMassK0S});
      histos_V0.add("InvMassLambda_PtRadius", "InvMassLambda_PtRadius", HistType::kTH3F, {axisPtsmall, axisRadiussmall, axisInvMassLambda});
      histos_V0.add("InvMassAntiLambda_PtRadius", "InvMassAntiLambda_PtRadius", HistType::kTH3F, {axisPtsmall, axisRadiussmall, axisInvMassAntiLambda});
      histos_V0.add("InvMassLambdaVsPtVsPA", "InvMassLambdaVsPtVsPA", HistType::kTH3F, {axisPt, axisV0PA, axisInvMassLambda});
      histos_V0.add("InvMassAntiLambdaVsPtVsPA", "InvMassAntiLambdaVsPtVsPA", HistType::kTH3F, {axisPt, axisV0PA, axisInvMassAntiLambda});
    }

    histos_Casc.add("CascSelectionSummary", "CascSelectionSummary", HistType::kTH1F, {{10, 0.f, 10.f}});
    histos_Casc.add("QA_XiMinusCandidates", "QA_XiMinusCandidates", HistType::kTH1F, {{10, 0.f, 10.f}});
    histos_Casc.add("XiProgSelections", "XiProgSelections", HistType::kTH2F, {{30, 0.5f, 30.5f}, {2, -2, 2}});
    histos_Casc.add("OmegaProgSelections", "OmegaProgSelections", HistType::kTH2F, {{30, 0.5f, 30.5f}, {2, -2, 2}});
    histos_Casc.add("CascCosPA", "CascCosPA", HistType::kTH2D, {{200, 0.6f, 1.0f}, {2, -2, 2}});
    histos_Casc.add("V0CosPA", "V0CosPA", HistType::kTH2D, {{300, 0.7f, 1.0f}, {2, -2, 2}});
    histos_Casc.add("V0CosPAToXi", "V0CosPAToXi", HistType::kTH2D, {{100, 0.9f, 1.0f}, {2, -2, 2}});
    histos_Casc.add("CascDecayLength", "CascDecayLength", HistType::kTH2F, {{100, 0.0f, 10.0f}, {2, -2, 2}});
    histos_Casc.add("CascDecayLengthXi", "CascDecayLengthXi", HistType::kTH2F, {{200, 0.0f, 20.0f}, {2, -2, 2}});
    histos_Casc.add("CascDecayLengthOmega", "CascDecayLengthOmega", HistType::kTH2F, {{200, 0.0f, 20.0f}, {2, -2, 2}});
    histos_Casc.add("CascRadius", "CascRadius", HistType::kTH2F, {{100, 0.0f, 10.0f}, {2, -2, 2}});
    histos_Casc.add("CascV0Radius", "CascV0Radius", HistType::kTH2D, {{100, 0.0f, 10.0f}, {2, -2, 2}});
    histos_Casc.add("CascyXi", "CascyXi", HistType::kTH2F, {{200, -2.0f, 2.0f}, {2, -2, 2}});
    histos_Casc.add("CascyOmega", "CascyOmega", HistType::kTH2F, {{200, -2.0f, 2.0f}, {2, -2, 2}});
    histos_Casc.add("CascCtauXi", "CascCtauXi", HistType::kTH2F, {{100, 0.0f, 100.0f}, {2, -2, 2}});
    histos_Casc.add("CascCtauOmega", "CascCtauOmega", HistType::kTH2F, {{100, 0.0f, 100.0f}, {2, -2, 2}});
    histos_Casc.add("V0Ctau", "V0Ctau", HistType::kTH2F, {{100, 0.0f, 100.0f}, {2, -2, 2}});
    histos_Casc.add("CascPt", "CascPt", HistType::kTH2F, {{100, 0.0f, 25.0f}, {2, -2, 2}});
    histos_Casc.add("DcaV0Daughters", "DcaV0Daughters", HistType::kTH2F, {{110, 0.0f, 2.2f}, {2, -2, 2}});
    histos_Casc.add("DcaCascDaughters", "DcaCascDaughters", HistType::kTH2F, {{110, 0.0f, 2.2f}, {2, -2, 2}});
    histos_Casc.add("DcaV0ToPV", "DcaV0ToPV", HistType::kTH2F, {{200, 0.0f, 2.f}, {2, -2, 2}});
    histos_Casc.add("DcaBachToPV", "DcaBachToPV", HistType::kTH2F, {{80, -0.2f, 0.2f}, {2, -2, 2}});
    histos_Casc.add("DcaPosToPV", "DcaPosToPV", HistType::kTH2F, {{80, -0.2f, 0.2f}, {2, -2, 2}});
    histos_Casc.add("DcaNegToPV", "DcaNegToPV", HistType::kTH2F, {{80, -0.2f, 0.2f}, {2, -2, 2}});
    histos_Casc.add("InvMassLambdaDaughter", "InvMassLambdaDaughter", HistType::kTH2F, {{100, 1.1f, 1.13f}, {2, -2, 2}});
    histos_Casc.add("InvMassXiPlus", "InvMassXiPlus", HistType::kTH3F, {{100, 0.f, 10.f}, {80, 1.28f, 1.36f}, {2, -1.0f, 1.0f}});
    histos_Casc.add("InvMassXiMinus", "InvMassXiMinus", HistType::kTH3F, {{100, 0.f, 10.f}, {80, 1.28f, 1.36f}, {2, -1.0f, 1.0f}});
    histos_Casc.add("InvMassXiPlus_Radius", "InvMassXiPlus_Radius", HistType::kTH2F, {{100, 0.f, 50.f}, {80, 1.28f, 1.36f}});
    histos_Casc.add("InvMassXiMinus_Radius", "InvMassXiMinus_Radius", HistType::kTH2F, {{100, 0.f, 50.f}, {80, 1.28f, 1.36f}});
    histos_Casc.add("InvMassOmegaPlus", "InvMassOmegaPlus", HistType::kTH3F, {{100, 0.f, 10.f}, {80, 1.63f, 1.71f}, {2, -1.0f, 1.0f}});
    histos_Casc.add("InvMassOmegaMinus", "InvMassOmegaMinus", HistType::kTH3F, {{100, 0.f, 10.f}, {80, 1.63f, 1.71f}, {2, -1.0f, 1.0f}});
    histos_Casc.add("InvMassXiPlusTrue", "InvMassXiPlusTrue", {HistType::kTH3F, {{100, 0.f, 10.f}, {100, 0.f, 50.f}, {80, 1.28f, 1.36f}}});
    histos_Casc.add("InvMassXiMinusTrue", "InvMassXiMinusTrue", {HistType::kTH3F, {{100, 0.f, 10.f}, {100, 0.f, 50.f}, {80, 1.28f, 1.36f}}});
    histos_Casc.add("InvMassOmegaPlusTrue", "InvMassOmegaPlusTrue", {HistType::kTH3F, {{100, 0.f, 10.f}, {100, 0.f, 50.f}, {80, 1.63f, 1.71f}}});
    histos_Casc.add("InvMassOmegaMinusTrue", "InvMassOmegaMinusTrue", {HistType::kTH3F, {{100, 0.f, 10.f}, {100, 0.f, 50.f}, {80, 1.63f, 1.71f}}});

    V0SelectionSummary->SetBinContent(1, V0_rapidity);
    V0SelectionSummary->SetBinContent(2, V0_cosPA);
    V0SelectionSummary->SetBinContent(3, V0_dcav0dau);
    V0SelectionSummary->SetBinContent(4, V0_dcapostopv);
    V0SelectionSummary->SetBinContent(5, V0_dcanegtopv);
    V0SelectionSummary->SetBinContent(6, V0_radius);
    V0SelectionSummary->SetBinContent(7, NSigmaV0Pion);
    V0SelectionSummary->SetBinContent(8, NSigmaV0Proton);
    V0SelectionSummary->SetBinContent(9, lifetimecut->get("lifetimecutLambda"));
    V0SelectionSummary->SetBinContent(10, lifetimecut->get("lifetimecutK0S"));

    V0SelectionSummary->GetXaxis()->SetBinLabel(1, "rapidity");
    V0SelectionSummary->GetXaxis()->SetBinLabel(2, "cosPA");
    V0SelectionSummary->GetXaxis()->SetBinLabel(3, "dcav0dau");
    V0SelectionSummary->GetXaxis()->SetBinLabel(4, "dcapostopv");
    V0SelectionSummary->GetXaxis()->SetBinLabel(5, "dcanegtopv");
    V0SelectionSummary->GetXaxis()->SetBinLabel(6, "radius");
    V0SelectionSummary->GetXaxis()->SetBinLabel(7, "NSigmaV0Pion");
    V0SelectionSummary->GetXaxis()->SetBinLabel(8, "NSigmaV0Proton");
    V0SelectionSummary->GetXaxis()->SetBinLabel(9, "lifetimecutLambda");
    V0SelectionSummary->GetXaxis()->SetBinLabel(10, "lifetimecutK0S");
  }

  ///////////////////////////////////////////////////
  ////////// Collisions QA - reconstructed //////////
  ///////////////////////////////////////////////////

  void processReconstructedEvent(soa::Join<aod::Collisions, aod::EvSels>::iterator const& Collision)
  {
    histos_eve.fill(HIST("hEventCounter"), 0.5);
    if (sel8 && !Collision.sel8()) {
      return;
    }
    histos_eve.fill(HIST("hEventCounter"), 1.5);
  }
  PROCESS_SWITCH(v0cascadesQA, processReconstructedEvent, "Process reconstructed level Event", true);

  ///////////////////////////////////////
  ////////// Collision QA - MC //////////
  ///////////////////////////////////////

  void processMcEvent(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles, const soa::SmallGroups<o2::soa::Join<o2::aod::Collisions, o2::aod::McCollisionLabels, o2::aod::EvSels>>& collisions)
  {
    histos_eve.fill(HIST("hEventCounterMC"), 0.5);

    std::vector<int64_t> SelectedEvents(collisions.size());
    int nevts = 0;
    for (const auto& collision : collisions) {
      if (sel8 && !collision.sel8()) {
        continue;
      }
      SelectedEvents[nevts++] = collision.mcCollision_as<aod::McCollisions>().globalIndex();
    }
    SelectedEvents.resize(nevts);

    const auto evtReconstructedAndSelected = std::find(SelectedEvents.begin(), SelectedEvents.end(), mcCollision.globalIndex()) != SelectedEvents.end();

    if (!evtReconstructedAndSelected) { // Check that the event is reconstructed and that the reconstructed events pass the selection
      return;
    }

    histos_eve.fill(HIST("hEventCounterMC"), 1.5);

    double posx = mcCollision.posX();
    double posy = mcCollision.posY();

    for (auto& mcparticle : mcParticles) {

      if (!mcparticle.has_daughters()) {
        continue;
      }

      double vx = 0;
      double vy = 0;
      for (auto& mcparticleDaughter0 : mcparticle.daughters_as<aod::McParticles>()) {
        vx = mcparticleDaughter0.vx() - posx;
        vy = mcparticleDaughter0.vy() - posy;
        if (vx != 0 && vy != 0)
          break;
      }
      double R_Decay = TMath::Sqrt(vx * vx + vy * vy);

      if (mcparticle.pdgCode() == PDG_t::kK0Short)
        histos_eve.fill(HIST("GenK0sPtVsEta"), mcparticle.pt(), mcparticle.eta());
      if (mcparticle.pdgCode() == PDG_t::kLambda0)
        histos_eve.fill(HIST("GenLambdaPtVsEta"), mcparticle.pt(), mcparticle.eta());
      if (mcparticle.pdgCode() == PDG_t::kLambda0Bar)
        histos_eve.fill(HIST("GenAntiLambdaPtVsEta"), mcparticle.pt(), mcparticle.eta());
      if (mcparticle.pdgCode() == PDG_t::kXiMinus)
        histos_eve.fill(HIST("GenXiMinusPtVsEta"), mcparticle.pt(), mcparticle.eta());
      if (mcparticle.pdgCode() == PDG_t::kXiPlusBar)
        histos_eve.fill(HIST("GenXiPlusPtVsEta"), mcparticle.pt(), mcparticle.eta());
      if (mcparticle.pdgCode() == PDG_t::kOmegaMinus)
        histos_eve.fill(HIST("GenOmegaMinusPtVsEta"), mcparticle.pt(), mcparticle.eta());
      if (mcparticle.pdgCode() == PDG_t::kOmegaPlusBar)
        histos_eve.fill(HIST("GenOmegaPlusPtVsEta"), mcparticle.pt(), mcparticle.eta());

      if (mcparticle.isPhysicalPrimary() && TMath::Abs(mcparticle.y()) < V0_rapidity) {
        if (mcparticle.pdgCode() == PDG_t::kK0Short)
          histos_eve.fill(HIST("GeneratedParticles"), 0.5, mcparticle.pt(), R_Decay);
        if (mcparticle.pdgCode() == PDG_t::kLambda0)
          histos_eve.fill(HIST("GeneratedParticles"), 2.5, mcparticle.pt(), R_Decay);
        if (mcparticle.pdgCode() == PDG_t::kLambda0Bar)
          histos_eve.fill(HIST("GeneratedParticles"), 4.5, mcparticle.pt(), R_Decay);
      }
      if (mcparticle.isPhysicalPrimary() && TMath::Abs(mcparticle.y()) < Casc_rapidity) {
        if (mcparticle.pdgCode() == PDG_t::kXiMinus)
          histos_eve.fill(HIST("GeneratedParticles"), 6.5, mcparticle.pt(), R_Decay);
        if (mcparticle.pdgCode() == PDG_t::kXiPlusBar)
          histos_eve.fill(HIST("GeneratedParticles"), 8.5, mcparticle.pt(), R_Decay);
        if (mcparticle.pdgCode() == PDG_t::kOmegaMinus)
          histos_eve.fill(HIST("GeneratedParticles"), 10.5, mcparticle.pt(), R_Decay);
        if (mcparticle.pdgCode() == PDG_t::kOmegaPlusBar)
          histos_eve.fill(HIST("GeneratedParticles"), 12.5, mcparticle.pt(), R_Decay);
      }
    }
  }
  PROCESS_SWITCH(v0cascadesQA, processMcEvent, "Process MC level Event", true);

  ////////////////////////////////////////////
  ////////// V0 QA - Reconstructed ///////////
  ////////////////////////////////////////////

  void processReconstructedV0(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::V0Datas const& fullV0s, DaughterTracks&)
  {
    if (sel8 && !collision.sel8()) {
      return;
    }

    int dauEtaFlag = 0;
    for (auto& v0 : fullV0s) {
      auto posdau = v0.posTrack_as<DaughterTracks>();
      auto negdau = v0.negTrack_as<DaughterTracks>();

      if (posdau.eta() < 0. && negdau.eta() < 0.) {
        dauEtaFlag = -1;
      } else if (posdau.eta() >= 0. && negdau.eta() >= 0.) {
        dauEtaFlag = 1;
      } else {
        dauEtaFlag = 0;
      }

      // check TPC
      if (checkDauTPC && (!posdau.hasTPC() || !negdau.hasTPC())) {
        continue;
      }

      Int_t posITSNhits = 0, negITSNhits = 0;
      for (unsigned int i = 0; i < 7; i++) {
        if (posdau.itsClusterMap() & (1 << i)) {
          posITSNhits++;
        }
        if (negdau.itsClusterMap() & (1 << i)) {
          negITSNhits++;
        }
      }

      histos_V0.fill(HIST("CosPA"), v0.v0cosPA());
      histos_V0.fill(HIST("V0Radius"), v0.v0radius());
      histos_V0.fill(HIST("V0DCANegToPV"), v0.dcanegtopv());
      histos_V0.fill(HIST("V0DCAPosToPV"), v0.dcapostopv());
      histos_V0.fill(HIST("V0DCAV0Daughters"), v0.dcaV0daughters());

      float decayLength = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * RecoDecay::sqrtSumOfSquares(v0.px(), v0.py(), v0.pz());
      histos_V0.fill(HIST("DecayLength"), decayLength);

      float CtauLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0;
      float CtauK0s = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short;

      if (v0.v0cosPA() > V0_cosPA &&
          v0.v0radius() > V0_radius &&
          v0.dcaV0daughters() < V0_dcav0dau &&
          TMath::Abs(v0.dcapostopv()) > V0_dcapostopv && TMath::Abs(v0.dcanegtopv()) > V0_dcanegtopv) {

        // K0Short
        if (TMath::Abs(v0.yK0Short()) < V0_rapidity &&
            CtauK0s < lifetimecut->get("lifetimecutK0S") &&
            TMath::Abs(posdau.tpcNSigmaPi()) < NSigmaV0Pion && TMath::Abs(negdau.tpcNSigmaPi()) < NSigmaV0Pion) {

          histos_V0.fill(HIST("CtauK0s"), CtauK0s);
          histos_V0.fill(HIST("DecayLengthK0s"), decayLength);
          histos_V0.fill(HIST("InvMassK0S"), v0.pt(), v0.mK0Short(), dauEtaFlag);
          histos_V0.fill(HIST("InvMassK0SVsPtVsPA"), v0.pt(), TMath::ACos(v0.v0cosPA()), v0.mK0Short());
          histos_V0.fill(HIST("V0DCAV0ToPVK0S"), v0.dcav0topv());
          if (doextraanalysis) {
            histos_V0.fill(HIST("InvMassK0S_Radius"), v0.v0radius(), v0.mK0Short());
            histos_V0.fill(HIST("InvMassK0S_PtRadius"), v0.pt(), v0.v0radius(), v0.mK0Short());
            histos_V0.fill(HIST("InvMassK0S_EtaDaughters"), posdau.eta(), negdau.eta(), v0.mK0Short());
            histos_V0.fill(HIST("InvMassK0S_PhiDaughters"), posdau.phi(), negdau.phi(), v0.mK0Short());
            histos_V0.fill(HIST("InvMassK0S_ITSMapDaughters"), posITSNhits, negITSNhits, v0.mK0Short());
          }
        }

        // Lambda
        if (TMath::Abs(v0.yLambda()) < V0_rapidity &&
            CtauLambda < lifetimecut->get("lifetimecutLambda") &&
            TMath::Abs(posdau.tpcNSigmaPr()) < NSigmaV0Proton && TMath::Abs(negdau.tpcNSigmaPi()) < NSigmaV0Pion) {

          histos_V0.fill(HIST("DecayLengthLambda"), decayLength);
          histos_V0.fill(HIST("CtauLambda"), CtauLambda);
          histos_V0.fill(HIST("InvMassLambda"), v0.pt(), v0.mLambda(), dauEtaFlag);
          if (doextraanalysis)
            histos_V0.fill(HIST("InvMassLambdaVsPtVsPA"), v0.pt(), TMath::ACos(v0.v0cosPA()), v0.mLambda());
          histos_V0.fill(HIST("V0DCAV0ToPVLambda"), v0.dcav0topv());
          if (v0.v0cosPA() > 0.999 && v0.dcaV0daughters() < 1 && TMath::Abs(v0.mK0Short() - o2::constants::physics::MassK0Short) > 0.012 && TMath::Abs(v0.mAntiLambda() - o2::constants::physics::MassLambda0) > 0.08 && TMath::Abs(v0.mLambda() - o2::constants::physics::MassLambda0) < 0.002) {
            histos_V0.fill(HIST("ResponsePionFromLambda"), v0.pt(), negdau.tpcNSigmaPi());
            histos_V0.fill(HIST("ResponseProtonFromLambda"), v0.pt(), posdau.tpcNSigmaPr());
          }
          if (doextraanalysis) {
            histos_V0.fill(HIST("InvMassLambda_Radius"), v0.v0radius(), v0.mLambda());
            histos_V0.fill(HIST("InvMassLambda_PtRadius"), v0.pt(), v0.v0radius(), v0.mLambda());
            histos_V0.fill(HIST("InvMassLambda_Ctau"), CtauLambda, v0.mLambda());
            histos_V0.fill(HIST("InvMassLambda_EtaDaughters"), posdau.eta(), negdau.eta(), v0.mLambda());
            histos_V0.fill(HIST("InvMassLambda_PhiDaughters"), posdau.phi(), negdau.phi(), v0.mLambda());
            histos_V0.fill(HIST("InvMassLambda_ITSMapDaughters"), posITSNhits, negITSNhits, v0.mLambda());
          }
        }

        // AntiLambda
        if (TMath::Abs(v0.yLambda()) < V0_rapidity &&
            CtauLambda < lifetimecut->get("lifetimecutLambda") &&
            TMath::Abs(posdau.tpcNSigmaPi()) < NSigmaV0Pion && TMath::Abs(negdau.tpcNSigmaPr()) < NSigmaV0Proton) {

          histos_V0.fill(HIST("DecayLengthAntiLambda"), decayLength);
          histos_V0.fill(HIST("CtauAntiLambda"), CtauLambda);
          histos_V0.fill(HIST("InvMassAntiLambda"), v0.pt(), v0.mAntiLambda(), dauEtaFlag);
          if (doextraanalysis)
            histos_V0.fill(HIST("InvMassAntiLambdaVsPtVsPA"), v0.pt(), TMath::ACos(v0.v0cosPA()), v0.mAntiLambda());
          histos_V0.fill(HIST("V0DCAV0ToPVAntiLambda"), v0.dcav0topv());
          if (doextraanalysis) {
            histos_V0.fill(HIST("InvMassAntiLambda_Radius"), v0.v0radius(), v0.mAntiLambda());
            histos_V0.fill(HIST("InvMassAntiLambda_PtRadius"), v0.pt(), v0.v0radius(), v0.mAntiLambda());
            histos_V0.fill(HIST("InvMassAntiLambda_Ctau"), CtauLambda, v0.mAntiLambda());
            histos_V0.fill(HIST("InvMassAntiLambda_EtaDaughters"), posdau.eta(), negdau.eta(), v0.mAntiLambda());
            histos_V0.fill(HIST("InvMassAntiLambda_PhiDaughters"), posdau.phi(), negdau.phi(), v0.mAntiLambda());
            histos_V0.fill(HIST("InvMassAntiLambda_ITSMapDaughters"), posITSNhits, negITSNhits, v0.mAntiLambda());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(v0cascadesQA, processReconstructedV0, "Process reconstructed level V0s", true);

  ////////////////////////////////
  ////////// V0 QA - MC //////////
  ////////////////////////////////

  void processMcV0(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Join<aod::V0Datas, aod::McV0Labels> const& fullV0s, aod::McParticles const&, MyTracksMC const&)
  {
    if (sel8 && !collision.sel8()) {
      return;
    }

    for (auto& v0 : fullV0s) {

      if (!v0.has_mcParticle()) {
        continue;
      }
      auto v0mcparticle = v0.mcParticle();
      Int_t lPDG = 0;
      if (TMath::Abs(v0mcparticle.pdgCode()) == 310 || TMath::Abs(v0mcparticle.pdgCode()) == PDG_t::kLambda0) {
        lPDG = v0mcparticle.pdgCode();
      }

      float CtauLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0;
      float CtauK0s = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short;

      if (v0.v0cosPA() > V0_cosPA &&
          v0.v0radius() > V0_radius &&
          v0.dcaV0daughters() < V0_dcav0dau &&
          TMath::Abs(v0.dcapostopv()) > V0_dcapostopv &&
          TMath::Abs(v0.dcanegtopv()) > V0_dcanegtopv) {

        // K0Short
        if (lPDG == 310) {
          if (TMath::Abs(v0.yK0Short()) < V0_rapidity && CtauK0s < lifetimecut->get("lifetimecutK0S")) {
            histos_V0.fill(HIST("InvMassK0STrue"), v0.pt(), v0.v0radius(), v0.mK0Short());
          }
        }
        if (lPDG == PDG_t::kLambda0) {
          if (TMath::Abs(v0.yLambda()) < V0_rapidity && CtauLambda < lifetimecut->get("lifetimecutLambda")) {
            histos_V0.fill(HIST("InvMassLambdaTrue"), v0.pt(), v0.v0radius(), v0.mLambda());
          }
        }
        if (lPDG == PDG_t::kLambda0Bar) {
          if (TMath::Abs(v0.yLambda()) < V0_rapidity && CtauLambda < lifetimecut->get("lifetimecutLambda")) {
            histos_V0.fill(HIST("InvMassAntiLambdaTrue"), v0.pt(), v0.v0radius(), v0.mAntiLambda());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(v0cascadesQA, processMcV0, "Process MC level V0s", false);

  //////////////////////////////////////
  ///// Cascade QA - Reconstructed /////
  //////////////////////////////////////

  void processReconstructedCascade(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::CascDataExt const& Cascades, aod::V0sLinked const&, aod::V0Datas const&, DaughterTracks&)
  {
    if (sel8 && !collision.sel8()) {
      return;
    }

    for (auto& casc : Cascades) {
      auto bachelor = casc.bachelor_as<DaughterTracks>();
      auto posdau = casc.posTrack_as<DaughterTracks>();
      auto negdau = casc.negTrack_as<DaughterTracks>();

      // check TPC
      if (checkDauTPC && (!posdau.hasTPC() || !negdau.hasTPC() || !bachelor.hasTPC())) {
        continue;
      }

      histos_Casc.fill(HIST("CascCosPA"), casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()), casc.sign());
      histos_Casc.fill(HIST("V0CosPA"), casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()), casc.sign());

      double v0cospatoxi = RecoDecay::cpa(array{casc.x(), casc.y(), casc.z()}, array{casc.xlambda(), casc.ylambda(), casc.zlambda()}, array{casc.pxpos() + casc.pxneg(), casc.pypos() + casc.pyneg(), casc.pzpos() + casc.pzneg()});

      histos_Casc.fill(HIST("V0CosPAToXi"), v0cospatoxi, casc.sign());
      histos_Casc.fill(HIST("CascRadius"), casc.cascradius(), casc.sign());
      histos_Casc.fill(HIST("CascV0Radius"), casc.v0radius(), casc.sign());
      histos_Casc.fill(HIST("CascyXi"), casc.yXi(), casc.sign());
      histos_Casc.fill(HIST("CascyOmega"), casc.yOmega(), casc.sign());

      float cascDecayLength = std::sqrt(std::pow(casc.x() - collision.posX(), 2) + std::pow(casc.y() - collision.posY(), 2) + std::pow(casc.z() - collision.posZ(), 2));
      histos_Casc.fill(HIST("CascDecayLength"), cascDecayLength, casc.sign());
      histos_Casc.fill(HIST("CascDecayLengthXi"), cascDecayLength, casc.sign());
      histos_Casc.fill(HIST("CascDecayLengthOmega"), cascDecayLength, casc.sign());

      float cascTotalMomentum = RecoDecay::sqrtSumOfSquares(casc.px(), casc.py(), casc.pz());
      float CtauXi = cascDecayLength / (cascTotalMomentum + 1E-10) * o2::constants::physics::MassXi0; // see O2Physics/Common/Core/MC.h for codes and names accepted
      float CtauOmega = cascDecayLength / (cascTotalMomentum + 1E-10) * o2::constants::physics::MassOmegaMinus;

      float v0TotalMomentum = RecoDecay::sqrtSumOfSquares(casc.pxpos() + casc.pxneg(), casc.pypos() + casc.pyneg(), casc.pzpos() + casc.pzneg());
      float v0DecayLength = std::sqrt(std::pow(casc.xlambda() - casc.x(), 2) + std::pow(casc.ylambda() - casc.y(), 2) + std::pow(casc.zlambda() - casc.z(), 2));
      float CtauV0 = v0DecayLength / (v0TotalMomentum + 1E-10) * o2::constants::physics::MassLambda0;

      histos_Casc.fill(HIST("CascCtauXi"), CtauXi, casc.sign());
      histos_Casc.fill(HIST("CascCtauOmega"), CtauOmega, casc.sign());
      histos_Casc.fill(HIST("V0Ctau"), CtauV0, casc.sign());
      histos_Casc.fill(HIST("CascPt"), casc.pt(), casc.sign());
      histos_Casc.fill(HIST("DcaV0Daughters"), casc.dcaV0daughters(), casc.sign());
      histos_Casc.fill(HIST("DcaCascDaughters"), casc.dcacascdaughters(), casc.sign());
      histos_Casc.fill(HIST("DcaV0ToPV"), casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()), casc.sign());
      histos_Casc.fill(HIST("DcaBachToPV"), casc.dcabachtopv(), casc.sign());
      histos_Casc.fill(HIST("DcaPosToPV"), casc.dcapostopv(), casc.sign());
      histos_Casc.fill(HIST("DcaNegToPV"), casc.dcanegtopv(), casc.sign());
      histos_Casc.fill(HIST("InvMassLambdaDaughter"), casc.mLambda(), casc.sign());

      if (casc.v0radius() > Casc_v0radius &&
          casc.cascradius() > Casc_cascradius &&
          casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > Casc_v0cospa &&
          casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) > Casc_casccospa &&
          TMath::Abs(casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ())) > Casc_dcav0topv &&
          TMath::Abs(casc.dcapostopv()) > Casc_dcapostopv && TMath::Abs(casc.dcanegtopv()) > Casc_dcanegtopv && TMath::Abs(casc.dcabachtopv()) > Casc_dcabachtopv &&
          casc.dcaV0daughters() < Casc_dcav0dau && casc.dcacascdaughters() < Casc_dcacascdau) {
        if (casc.sign() < 0) {
          if (TMath::Abs(posdau.tpcNSigmaPr()) < NSigmaCascProton && TMath::Abs(negdau.tpcNSigmaPi()) < NSigmaCascPion) {
            if (TMath::Abs(casc.yXi()) < Casc_rapidity && TMath::Abs(bachelor.tpcNSigmaPi()) < NSigmaCascPion) {
              histos_Casc.fill(HIST("InvMassXiMinus"), casc.pt(), casc.mXi(), casc.eta());
              histos_Casc.fill(HIST("InvMassXiMinus_Radius"), casc.cascradius(), casc.mXi());
            }
            if (TMath::Abs(casc.yOmega()) < Casc_rapidity && TMath::Abs(bachelor.tpcNSigmaKa()) < NSigmaCascKaon) {
              histos_Casc.fill(HIST("InvMassOmegaMinus"), casc.pt(), casc.mOmega(), casc.eta());
            }
          }
        } else {
          if (TMath::Abs(posdau.tpcNSigmaPi()) < NSigmaCascPion && TMath::Abs(negdau.tpcNSigmaPr()) < NSigmaCascProton) {
            if (TMath::Abs(casc.yXi()) < Casc_rapidity && TMath::Abs(bachelor.tpcNSigmaPi()) < NSigmaCascPion) {
              histos_Casc.fill(HIST("InvMassXiPlus"), casc.pt(), casc.mXi(), casc.eta());
              histos_Casc.fill(HIST("InvMassXiPlus_Radius"), casc.cascradius(), casc.mXi());
            }
            if (TMath::Abs(casc.yOmega()) < Casc_rapidity && TMath::Abs(bachelor.tpcNSigmaKa()) < NSigmaCascKaon) {
              histos_Casc.fill(HIST("InvMassOmegaPlus"), casc.pt(), casc.mOmega(), casc.eta());
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(v0cascadesQA, processReconstructedCascade, "Process reconstructed level Cascades", true);

  //////////////////////////////////////
  ////////// Cascade QA - MC ///////////
  //////////////////////////////////////

  void processMcCascade(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::CascDataExt const& Cascades, aod::V0sLinked const&, aod::V0Datas const&, MyTracksMC const&, aod::McParticles const&)
  {
    if (sel8 && !collision.sel8()) {
      return;
    }

    for (auto& casc : Cascades) {

      histos_Casc.fill(HIST("QA_XiMinusCandidates"), 0.5);

      if (casc.v0radius() > Casc_v0radius &&
          casc.cascradius() > Casc_cascradius &&
          casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > Casc_v0cospa &&
          casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) > Casc_casccospa &&
          TMath::Abs(casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ())) > Casc_dcav0topv &&
          TMath::Abs(casc.dcapostopv()) > Casc_dcapostopv && TMath::Abs(casc.dcanegtopv()) > Casc_dcanegtopv && TMath::Abs(casc.dcabachtopv()) > Casc_dcabachtopv &&
          casc.dcaV0daughters() < Casc_dcav0dau && casc.dcacascdaughters() < Casc_dcacascdau) {

        histos_Casc.fill(HIST("QA_XiMinusCandidates"), 1.5);

        auto reconegtrack = casc.negTrack_as<MyTracksMC>();
        auto recopostrack = casc.posTrack_as<MyTracksMC>();
        auto recobachelor = casc.bachelor_as<MyTracksMC>();
        if (!reconegtrack.has_mcParticle() || !recopostrack.has_mcParticle() || !recobachelor.has_mcParticle()) {
          continue;
        }
        histos_Casc.fill(HIST("QA_XiMinusCandidates"), 2.5);

        auto mcnegtrack = reconegtrack.mcParticle_as<aod::McParticles>();
        auto mcpostrack = recopostrack.mcParticle_as<aod::McParticles>();
        auto bachelor = recobachelor.mcParticle_as<aod::McParticles>();
        if (!mcnegtrack.has_mothers() || !mcpostrack.has_mothers() || !bachelor.has_mothers()) {
          continue;
        }
        histos_Casc.fill(HIST("QA_XiMinusCandidates"), 3.5);

        for (auto& particleMotherOfBach : bachelor.mothers_as<aod::McParticles>()) {
          for (auto& particleMotherOfNeg : mcnegtrack.mothers_as<aod::McParticles>()) {
            for (auto& particleMotherOfPos : mcpostrack.mothers_as<aod::McParticles>()) {
              for (auto& particleMotherOfV0 : particleMotherOfNeg.mothers_as<aod::McParticles>()) {

                bool MomOfBachIsPrimary = particleMotherOfBach.isPhysicalPrimary();
                bool MomOfNegIsPrimary = particleMotherOfNeg.isPhysicalPrimary();
                bool MomOfPosIsPrimary = particleMotherOfPos.isPhysicalPrimary();

                bool isXiMinusCascade = MomOfBachIsPrimary && !(MomOfNegIsPrimary) && !(MomOfPosIsPrimary) &&
                                        particleMotherOfNeg == particleMotherOfPos &&
                                        particleMotherOfV0 == particleMotherOfBach &&
                                        particleMotherOfBach.pdgCode() == PDG_t::kXiMinus &&
                                        bachelor.pdgCode() == PDG_t::kPiMinus &&
                                        particleMotherOfNeg.pdgCode() == PDG_t::kLambda0 &&
                                        mcnegtrack.pdgCode() == PDG_t::kPiMinus &&
                                        mcpostrack.pdgCode() == PDG_t::kProton;

                bool isOmegaMinusCascade = MomOfBachIsPrimary && !(MomOfNegIsPrimary) && !(MomOfPosIsPrimary) &&
                                           particleMotherOfNeg == particleMotherOfPos &&
                                           particleMotherOfV0 == particleMotherOfBach &&
                                           particleMotherOfBach.pdgCode() == PDG_t::kOmegaMinus &&
                                           bachelor.pdgCode() == PDG_t::kKMinus &&
                                           particleMotherOfNeg.pdgCode() == PDG_t::kLambda0 &&
                                           mcnegtrack.pdgCode() == PDG_t::kPiMinus &&
                                           mcpostrack.pdgCode() == PDG_t::kProton;

                bool isXiPlusCascade = MomOfBachIsPrimary && !(MomOfNegIsPrimary) && !(MomOfPosIsPrimary) &&
                                       particleMotherOfNeg == particleMotherOfPos &&
                                       particleMotherOfV0 == particleMotherOfBach &&
                                       particleMotherOfBach.pdgCode() == PDG_t::kXiPlusBar &&
                                       bachelor.pdgCode() == PDG_t::kPiPlus &&
                                       particleMotherOfNeg.pdgCode() == PDG_t::kLambda0Bar &&
                                       mcnegtrack.pdgCode() == PDG_t::kProtonBar &&
                                       mcpostrack.pdgCode() == PDG_t::kPiPlus;

                bool isOmegaPlusCascade = MomOfBachIsPrimary && !(MomOfNegIsPrimary) && !(MomOfPosIsPrimary) &&
                                          particleMotherOfNeg == particleMotherOfPos &&
                                          particleMotherOfV0 == particleMotherOfBach &&
                                          particleMotherOfBach.pdgCode() == PDG_t::kOmegaPlusBar &&
                                          bachelor.pdgCode() == PDG_t::kKPlus &&
                                          particleMotherOfNeg.pdgCode() == PDG_t::kLambda0Bar &&
                                          mcnegtrack.pdgCode() == PDG_t::kProtonBar &&
                                          mcpostrack.pdgCode() == PDG_t::kPiPlus;

                if (isXiMinusCascade) {
                  histos_Casc.fill(HIST("QA_XiMinusCandidates"), 8.5);
                  if (TMath::Abs(casc.yXi()) < Casc_rapidity) {
                    histos_Casc.fill(HIST("InvMassXiMinusTrue"), casc.pt(), casc.cascradius(), casc.mXi());
                  }
                }

                if (isOmegaMinusCascade) {
                  if (TMath::Abs(casc.yOmega()) < Casc_rapidity) {
                    histos_Casc.fill(HIST("InvMassOmegaMinusTrue"), casc.pt(), casc.cascradius(), casc.mOmega());
                  }
                }

                if (isXiPlusCascade) {
                  if (TMath::Abs(casc.yXi()) < Casc_rapidity) {
                    histos_Casc.fill(HIST("InvMassXiPlusTrue"), casc.pt(), casc.cascradius(), casc.mXi());
                  }
                }

                if (isOmegaPlusCascade) {
                  if (TMath::Abs(casc.yOmega()) < Casc_rapidity) {
                    histos_Casc.fill(HIST("InvMassOmegaPlusTrue"), casc.pt(), casc.cascradius(), casc.mOmega());
                  }
                }

                // QA section - XiMinus
                if (MomOfBachIsPrimary && !(MomOfNegIsPrimary) && !(MomOfPosIsPrimary)) {
                  histos_Casc.fill(HIST("QA_XiMinusCandidates"), 4.5);
                }
                if ((particleMotherOfNeg.pdgCode() == particleMotherOfPos.pdgCode())) {
                  histos_Casc.fill(HIST("QA_XiMinusCandidates"), 5.5);
                }
                if ((particleMotherOfV0 == particleMotherOfBach)) {
                  histos_Casc.fill(HIST("QA_XiMinusCandidates"), 6.5);
                }
                if (particleMotherOfBach.pdgCode() == PDG_t::kXiMinus &&
                    bachelor.pdgCode() == PDG_t::kPiMinus &&
                    particleMotherOfNeg.pdgCode() == PDG_t::kLambda0 &&
                    mcnegtrack.pdgCode() == PDG_t::kPiMinus &&
                    mcpostrack.pdgCode() == PDG_t::kProton) {
                  histos_Casc.fill(HIST("QA_XiMinusCandidates"), 7.5);
                }
              }
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(v0cascadesQA, processMcCascade, "Process MC level Cascades", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<v0cascadesQA>(cfgc)};
}
