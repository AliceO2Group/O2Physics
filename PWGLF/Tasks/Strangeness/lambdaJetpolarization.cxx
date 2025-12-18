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
///

// o2-linter: disable=name/workflow-file

/// \author Youpeng Su (yousu@cern.ch)
#include "PWGLF/DataModel/lambdaJetpolarization.h"

#include "PWGJE/Core/JetBkgSubUtils.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

#include "Math/GenVector/Boost.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TProfile2D.h"
#include <TFile.h>
// #include <TLorentzVector.h>
#include <TMatrixD.h>
#include <TTree.h>

#include <fastjet/AreaDefinition.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/Selector.hh>
#include <fastjet/tools/JetMedianBackgroundEstimator.hh>
#include <fastjet/tools/Subtractor.hh>

#include <cmath>
// #include <iostream>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct LambdaJetpolarization {

  HistogramRegistry registryData{"registryData", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  Configurable<double> zVtx{"zVtx", 10.0, "Maximum zVertex"};
  Configurable<double> rJet{"rJet", 0.4, "Jet resolution parameter R"};
  Configurable<float> etaMin{"etaMin", -0.9f, "eta min"};
  Configurable<float> etaMax{"etaMax", +0.9f, "eta max"};
  Configurable<double> deltaEtaEdge{"deltaEtaEdge", 0.00, "eta gap from the edge"};
  // track parameters
  Configurable<float> minITSnCls{"minITSnCls", 2.0f, "min number of ITS clusters"};
  Configurable<float> minTPCnClsFound{"minTPCnClsFound", 80.0f, "min number of found TPC clusters"};
  Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 80.0f, "min number of TPC crossed rows"};
  Configurable<float> minTpcNcrossedRowsOverFindable{"minTpcNcrossedRowsOverFindable", 0.8, "crossed rows/findable"};
  Configurable<bool> requireTOF{"requireTOF", false, "require TOF hit"};
  Configurable<bool> requireITS{"requireITS", false, "require ITS hit"};
  Configurable<bool> requireMaxTPCSharedCls{"requireMaxTPCSharedCls", false, "require max TPC shared clusters"};
  Configurable<float> maxTPCSharedCls{"maxTPCSharedCls", 100, "maxTPCSharedCls"};
  Configurable<float> maxChi2TPC{"maxChi2TPC", 4, "maxChi2TPC"};
  Configurable<float> maxChi2ITS{"maxChi2ITS", 36, "maxChi2ITS"};

  Configurable<float> ptMinV0Proton{"ptMinV0Proton", 0.3f, "pt min of proton from V0"};
  Configurable<float> ptMaxV0Proton{"ptMaxV0Proton", 10.0f, "pt max of proton from V0"};
  Configurable<float> ptMinV0Pion{"ptMinV0Pion", 0.1f, "pt min of pion from V0"};
  Configurable<float> ptMaxV0Pion{"ptMaxV0Pion", 1.5f, "pt max of pion from V0"};

  Configurable<float> nsigmaTPCmin{"nsigmaTPCmin", -5.0f, "Minimum nsigma TPC"};
  Configurable<float> nsigmaTPCmax{"nsigmaTPCmax", +5.0f, "Maximum nsigma TPC"};
  Configurable<float> nsigmaTOFmin{"nsigmaTOFmin", -5.0f, "Minimum nsigma TOF"};
  Configurable<float> nsigmaTOFmax{"nsigmaTOFmax", +5.0f, "Maximum nsigma TOF"};
  Configurable<double> cfgtrkMinPt{"cfgtrkMinPt", 0.10, "set track min pT"};

  // v0 parameters
  Configurable<float> v0Ptmin{"v0Ptmin", 0.6f, "Minimum V0 pT"};
  Configurable<float> v0cospaMin{"v0cospaMin", 0.995f, "Minimum V0 CosPA"};
  Configurable<float> v0cospainit{"v0cospainit", 0.97f, "Minimum V0 CosPA"};
  Configurable<float> minimumV0Radius{"minimumV0Radius", 0.5f, "Minimum V0 Radius"};
  Configurable<float> maximumV0Radius{"maximumV0Radius", 100000.0f, "Maximum V0 Radius"};
  Configurable<float> dcaV0DaughtersMax{"dcaV0DaughtersMax", 1.0f, "Maximum DCA Daughters"};
  Configurable<float> dcanegtoPVmin{"dcanegtoPVmin", 0.1f, "Minimum DCA Neg To PV"};
  Configurable<float> dcapostoPVmin{"dcapostoPVmin", 0.1f, "Minimum DCA Pos To PV"};
  Configurable<float> v0radius{"v0radius", 0.0, "Radius"};
  Configurable<float> dcav0dau{"dcav0dau", 1, "DCA V0 Daughters"};
  Configurable<float> dcanegtopv{"dcanegtopv", 0.05, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", 0.05, "DCA Pos To PV"};

  // jet selection
  Configurable<float> cfgjetPtMin{"cfgjetPtMin", 8.0, "minimum jet pT cut"};
  Configurable<bool> ispassdTrackSelectionForJetReconstruction{"ispassdTrackSelectionForJetReconstruction", 1, "do track selection"};

  // v0Event selection
  Configurable<bool> sel8{"sel8", 1, "Apply sel8 event selection"};
  Configurable<bool> isTriggerTVX{"isTriggerTVX", 1, "TVX trigger"};
  Configurable<bool> iscutzvertex{"iscutzvertex", 1, "Accepted z-vertex range (cm)"};
  Configurable<bool> isNoTimeFrameBorder{"isNoTimeFrameBorder", 1, "TF border cut"};
  Configurable<bool> isNoITSROFrameBorder{"isNoITSROFrameBorder", 1, "ITS ROF border cut"};
  Configurable<bool> isVertexTOFmatched{"isVertexTOFmatched", 1, "Is Vertex TOF matched"};
  Configurable<bool> isNoSameBunchPileup{"isNoSameBunchPileup", 0, "isNoSameBunchPileup"};
  Configurable<bool> isGoodZvtxFT0vsPV{"isGoodZvtxFT0vsPV", 1, "isGoodZvtxFT0vsPV"};
  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range (cm)"};
  Configurable<float> cTau{"cTau", 30, "C tau (cm)"};
  Configurable<bool> requirepassedSingleTrackSelection{"requirepassedSingleTrackSelection", false, "requirepassedSingleTrackSelection"};
  Configurable<float> v0TracketaMin{"v0TracketaMin", -0.8f, "eta min track"};
  Configurable<float> v0TracketaMax{"v0TracketaMax", +0.8f, "eta max track"};
  Configurable<bool> requireTPC{"requireTPC", true, "require TPC hit"};
  Configurable<float> yMin{"yMin", -0.5f, "minimum y"};
  Configurable<float> yMax{"yMax", +0.5f, "maximum y"};
  Configurable<float> v0rejLambda{"v0rejLambda", 0.01, "V0 rej Lambda"};
  Configurable<float> v0accLambda{"v0accLambda", 0.075, "V0 acc Lambda"};
  Configurable<bool> ifinitpasslambda{"ifinitpasslambda", 0, "ifinitpasslambda"};
  Configurable<bool> ifpasslambda{"ifpasslambda", 1, "ifpasslambda"};
  Configurable<float> paramArmenterosCut{"paramArmenterosCut", 0.2, "parameter Armenteros Cut"};
  Configurable<bool> doArmenterosCut{"doArmenterosCut", 1, "do Armenteros Cut"};
  Configurable<bool> noSameBunchPileUp{"noSameBunchPileUp", true, "reject SameBunchPileUp"};
  Configurable<int> v0TypeSelection{"v0TypeSelection", 1, "select on a certain V0 type (leave negative if no selection desired)"};
  Configurable<bool> notITSAfterburner{"notITSAfterburner", 0, "notITSAfterburner"};
  Configurable<bool> doQA{"doQA", 1, "fill QA histograms"};
  Configurable<bool> evSel{"evSel", 1, "evSel"};
  Configurable<bool> hasTOF2Leg{"hasTOF2Leg", 0, "hasTOF2Leg"};
  Configurable<bool> hasTOF1Leg{"hasTOF1Leg", 0, "hasTOF1Leg"};

  // Jet background subtraction
  JetBkgSubUtils backgroundSub;
  void init(InitContext const&)
  {
    const AxisSpec axisPx{100, -10, 10, "#px (GeV/c)"};
    const AxisSpec axisPy{100, -10, 10, "#py (GeV/c)"};
    const AxisSpec axisPz{100, -10, 10, "#pz (GeV/c)"};
    const AxisSpec axisPT{200, 0, 50, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisPhi{200, -TMath::Pi(), TMath::Pi(), "#Phi"};
    const AxisSpec axisDeltaPhi{200, 0, TMath::Pi(), "#Phi"};
    const AxisSpec axisTheta{100, -TMath::Pi(), TMath::Pi(), "#Theta"};
    const AxisSpec axisMass{100, 0.9, 1.0, "Mass(GeV/c^{2})"};
    const AxisSpec axisCostheta{100, -1, 1, "Cos(#theta^{*}_{p})"};
    const AxisSpec axisSinPhi{100, -1, 1, "Sin(#phi^{*}_{p})"};

    const AxisSpec JetaxisEta{30, -1.5, +1.5, "#eta"};
    const AxisSpec JetaxisPhi{200, -1, +7, "#phi"};
    const AxisSpec JetaxisPt{200, 0, +200, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec ptAxis{100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec invMassLambdaAxis{200, 1.016, 1.216, "m_{p#pi} (GeV/#it{c}^{2})"};

    ConfigurableAxis tprofile2DaxisPt{"tprofile2DaxisPt", {VARIABLE_WIDTH, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.2, 3.7, 4.2, 5, 6, 8, 10, 12}, "pt axis for histograms"};
    ConfigurableAxis tprofile2DaxisMass{"tprofile2DaxisMass", {VARIABLE_WIDTH, 1.10068, 1.10668, 1.11068, 1.11268, 1.11368, 1.11468, o2::constants::physics::MassLambda, 1.11668, 1.11768, 1.11868, 1.12068, 1.12468, 1.13068}, "Mass axis for histograms"};

    registryData.add("number_of_events_vsmultiplicity", "number of events in data vs multiplicity", HistType::kTH1D, {{101, 0, 101, "Multiplicity percentile"}});
    registryData.add("h_track_pt", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", kTH1F, {{200, 0., 200.}});
    registryData.add("h_track_eta", "track #eta;#eta_{track};entries", kTH1F, {{100, -1.f, 1.f}});
    registryData.add("h_track_phi", "track #varphi;#varphi_{track};entries", kTH1F, {{80, -1.f, 7.f}});
    registryData.add("h_track_pt_sel", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", kTH1F, {{200, 0., 200.}});
    registryData.add("h_track_eta_sel", "track #eta;#eta_{track};entries", kTH1F, {{100, -1.f, 1.f}});
    registryData.add("h_track_phi_sel", "track #varphi;#varphi_{track};entries", kTH1F, {{80, -1.f, 7.f}});

    registryData.add("FJetaHistogram", "FJetaHistogram", kTH1F, {JetaxisEta});
    registryData.add("FJphiHistogram", "FJphiHistogram", kTH1F, {JetaxisPhi});
    registryData.add("FJptHistogram", "FJptHistogram", kTH1F, {JetaxisPt});
    registryData.add("nJetsPerEvent", "nJetsPerEvent", kTH1F, {{10, 0.0, 10.0}});
    registryData.add("nJetsPerEventsel", "nJetsPerEventsel", kTH1F, {{10, 0.0, 10.0}});
    registryData.add("nV0sPerEvent", "nV0sPerEvent", kTH1F, {{10, 0.0, 10.0}});
    registryData.add("FJetaHistogramsel", "FJetaHistogramsel", kTH1F, {JetaxisEta});
    registryData.add("FJphiHistogramsel", "FJphiHistogramsel", kTH1F, {JetaxisPhi});
    registryData.add("FJptHistogramsel", "FJptHistogramsel", kTH1F, {JetaxisPt});

    registryData.add("FLeadingJetaHistogramsel", "FLeadingJetaHistogramsel", kTH1F, {JetaxisEta});
    registryData.add("FLeadingJphiHistogramsel", "FLeadingJphiHistogramsel", kTH1F, {JetaxisPhi});
    registryData.add("FLeadingJptHistogramsel", "FLeadingJptHistogramsel", kTH1F, {JetaxisPt});

    registryData.add("LambdaPtMass", "LambdaPtMass", HistType::kTH2F, {ptAxis, invMassLambdaAxis});
    registryData.add("AntiLambdaPtMass", "AntiLambdaPtMass", HistType::kTH2F, {ptAxis, invMassLambdaAxis});

    registryData.add("hMassLambda", "hMassLambda", {HistType::kTH1F, {invMassLambdaAxis}});
    registryData.add("hMassAntiLambda", "hMassAntiLambda", {HistType::kTH1F, {invMassLambdaAxis}});
    registryData.add("V0pTInLab", "V0pTInLab", kTH1F, {axisPT});

    registryData.add("V0pxInLab", "V0pxInLab", kTH1F, {axisPx});
    registryData.add("V0pyInLab", "V0pyInLab", kTH1F, {axisPy});
    registryData.add("V0pzInLab", "V0pzInLab", kTH1F, {axisPz});

    registryData.add("V0pxInRest_frame", "V0pxInRest_frame", kTH1F, {axisPx});
    registryData.add("V0pyInRest_frame", "V0pyInRest_frame", kTH1F, {axisPy});
    registryData.add("V0pzInRest_frame", "V0pzInRest_frame", kTH1F, {axisPz});

    registryData.add("V0pxInJetframe", "V0pxInJetframe", kTH1F, {axisPx});
    registryData.add("V0pyInJetframe", "V0pyInJetframe", kTH1F, {axisPy});
    registryData.add("V0pzInJetframe", "V0pzInJetframe", kTH1F, {axisPz});

    registryData.add("protonQA/V0protonpxInLab", "V0protonpxInLab", kTH1F, {axisPx});
    registryData.add("protonQA/V0protonpyInLab", "V0protonpyInLab", kTH1F, {axisPy});
    registryData.add("protonQA/V0protonpzInLab", "V0protonpzInLab", kTH1F, {axisPz});
    registryData.add("protonQA/V0protonMassInLab", "V0protonMassInLab", kTH1F, {axisMass});
    registryData.add("protonQA/V0protonphiInLab", "V0protonphiInLab", kTH1F, {axisPhi});
    registryData.add("protonQA/V0protonthetaInLab", "V0protonthetaInLab", kTH1F, {axisTheta});
    registryData.add("protonQA/V0protoncosthetaInLab", "V0protoncosthetaInLab", kTH1F, {axisCostheta});
    registryData.add("protonQA/profileprotonsinthetaInLab", "Invariant Mass vs sin(theta)", {HistType::kTProfile, {{200, 0.9, 1.2}}});
    registryData.add("protonQA/profileprotonsinphiInLab", "Invariant Mass vs sin(phi)", {HistType::kTProfile, {{200, 0.9, 1.2}}});
    registryData.add("protonQA/profileprotoncosSquarethetaInLab", "Invariant Mass vs cos^2(theta)", {HistType::kTProfile, {{200, 0.9, 1.2}}});
    registryData.add("JetQA/JetthetaInLab", "JetthetaInLab", kTH1F, {axisTheta});
    registryData.add("JetQA/JetphiInLab", "JetphiInLab", kTH1F, {axisPhi});
    registryData.add("JetQA/JetpxInLab", "JetpxInLab", kTH1F, {axisPx});
    registryData.add("JetQA/JetpyInLab", "JetpyInLab", kTH1F, {axisPy});
    registryData.add("JetQA/JetpzInLab", "JetpzInLab", kTH1F, {axisPz});
    registryData.add("JetQA/JetptInLab", "JetptInLab", kTH1F, {axisPT});

    registryData.add("protonQA/V0protonpxInRest_frame", "V0protonpxInRest_frame", kTH1F, {axisPx});
    registryData.add("protonQA/V0protonpyInRest_frame", "V0protonpyInRest_frame", kTH1F, {axisPy});
    registryData.add("protonQA/V0protonpzInRest_frame", "V0protonpzInRest_frame", kTH1F, {axisPz});
    registryData.add("protonQA/V0protonMassInRest_frame", "V0protonMassInRest_frame", kTH1F, {axisMass});
    registryData.add("protonQA/V0protonphiInRest_frame", "V0protonphiInRest_frame", kTH1F, {axisPhi});
    registryData.add("protonQA/V0protonthetaInRest_frame", "V0protonthetaInRest_frame", kTH1F, {axisTheta});
    registryData.add("protonQA/V0protoncosthetaInV0frame", "V0protoncosthetaInV0frame", kTH1F, {axisCostheta});
    registryData.add("protonQA/profileprotonsinthetaInV0frame", "Invariant Mass vs sin(theta)", {HistType::kTProfile, {{200, 0.9, 1.2}}});
    registryData.add("protonQA/profileprotonsinphiInV0frame", "Invariant Mass vs sin(phi)", {HistType::kTProfile, {{200, 0.9, 1.2}}});
    registryData.add("protonQA/profileprotoncosSquarethetaInV0frame", "Invariant Mass vs cos^2(theta)", {HistType::kTProfile, {{200, 0.9, 1.2}}});
    registryData.add("JetQA/JetthetaInV0", "JetthetaInV0", kTH1F, {axisTheta});
    registryData.add("JetQA/JetphiInV0", "JetphiInV0", kTH1F, {axisPhi});
    registryData.add("JetQA/JetpxInV0", "JetpxInV0", kTH1F, {axisPx});
    registryData.add("JetQA/JetpyInV0", "JetpyInV0", kTH1F, {axisPy});
    registryData.add("JetQA/JetpzInV0", "JetpzInV0", kTH1F, {axisPz});
    registryData.add("JetQA/JetptInV0", "JetptInV0", kTH1F, {axisPT});

    registryData.add("protonQA/V0protonpxInJetframe", "V0protonpxInJetframe", kTH1F, {axisPx});
    registryData.add("protonQA/V0protonpyInJetframe", "V0protonpyInJetframe", kTH1F, {axisPy});
    registryData.add("protonQA/V0protonpzInJetframe", "V0protonpzInJetframe", kTH1F, {axisPz});
    registryData.add("protonQA/V0protonphiInJetframe", "V0protonphiInJetframe", kTH1F, {axisPhi});
    registryData.add("protonQA/V0protonthetaInJetframe", "V0protonthetaInJetframe", kTH1F, {axisTheta});
    registryData.add("protonQA/V0protoncosthetaInJetframe", "V0protoncosthetaInJetframe", kTH1F, {axisCostheta});
    registryData.add("protonQA/profileprotonsinthetaInJetframe", "Invariant Mass vs sin(theta)", {HistType::kTProfile, {{200, 0.9, 1.2}}});
    registryData.add("protonQA/profileprotonsinphiInJetframe", "Invariant Mass vs sin(phi)", {HistType::kTProfile, {{200, 0.9, 1.2}}});
    registryData.add("protonQA/profileprotoncosSquarethetaInJetframe", "Invariant Mass vs cos^2(theta)", {HistType::kTProfile, {{200, 0.9, 1.2}}});
    registryData.add("protonQA/V0protonMassInJetframe", "V0protonMassInJetframe", kTH1F, {axisMass});
    registryData.add("JetQA/JetthetaInJetframe", "JetthetaInJetframe", kTH1F, {axisTheta});
    registryData.add("JetQA/JetphiInJetframe", "JetphiInJetframe", kTH1F, {axisPhi});
    registryData.add("JetQA/JetpxInJetframe", "JetpxInJetframe", kTH1F, {axisPx});
    registryData.add("JetQA/JetpyInJetframe", "JetpyInJetframe", kTH1F, {axisPy});
    registryData.add("JetQA/JetpzInJetframe", "JetpzInJetframe", kTH1F, {axisPz});
    registryData.add("JetQA/JetptInJetframe", "JetptInJetframe", kTH1F, {axisPT});

    registryData.add("protonQA/V0protonpxInJetV0frame", "V0protonpxInJetV0frame", kTH1F, {axisPx});
    registryData.add("protonQA/V0protonpyInJetV0frame", "V0protonpyInJetV0frame", kTH1F, {axisPy});
    registryData.add("protonQA/V0protonpzInJetV0frame", "V0protonpzInJetV0frame", kTH1F, {axisPz});
    registryData.add("protonQA/V0protonphiInJetV0frame", "V0protonphiInJetV0frame", kTH1F, {axisPhi});
    registryData.add("protonQA/V0protonthetaInJetV0frame", "V0protonthetaInJetV0frame", kTH1F, {axisTheta});
    registryData.add("protonQA/V0protoncosthetaInJetV0", "V0protoncosthetaInJetV0", kTH1F, {axisCostheta});
    registryData.add("protonQA/V0protonMassInJetV0frame", "V0protonMassInJetV0frame", kTH1F, {axisMass});
    registryData.add("protonQA/profileprotonsinthetaInJetV0frame", "Invariant Mass vs sin(theta)", {HistType::kTProfile, {{200, 0.9, 1.2}}});
    registryData.add("protonQA/profileprotonsinphiInJetV0frame", "Invariant Mass vs sin(phi)", {HistType::kTProfile, {{200, 0.9, 1.2}}});
    registryData.add("protonQA/profileprotoncosSquarethetaInJetV0frame", "Invariant Mass vs cos^2(theta)", {HistType::kTProfile, {{200, 0.9, 1.2}}});
    registryData.add("JetQA/JetthetaInJetV0frame", "JetthetaInJetV0frame", kTH1F, {axisTheta});
    registryData.add("JetQA/JetphiInJetV0frame", "JetphiInJetV0frame", kTH1F, {axisPhi});
    registryData.add("JetQA/JetpxInJetV0frame", "JetpxInJetV0frame", kTH1F, {axisPx});
    registryData.add("JetQA/JetpyInJetV0frame", "JetpyInJetV0frame", kTH1F, {axisPy});
    registryData.add("JetQA/JetpzInJetV0frame", "JetpzInJetV0frame", kTH1F, {axisPz});
    registryData.add("JetQA/JetptInJetV0frame", "JetptInJetV0frame", kTH1F, {axisPT});

    registryData.add("V0LambdapxInJetV0frame", "V0LambdapxInJetV0frame", kTH1F, {axisPx});
    registryData.add("V0LambdapyInJetV0frame", "V0LambdapyInJetV0frame", kTH1F, {axisPy});
    registryData.add("V0LambdapzInJetV0frame", "V0LambdapzInJetV0frame", kTH1F, {axisPz});

    registryData.add("hprotonPhi", "hprotonPhi", kTH1F, {axisPhi});
    registryData.add("hantiprotonPhi", "hantiprotonPhi", kTH1F, {axisPhi});

    registryData.add("hLambdamassandSinPhi", "hLambdamassandSinPhi", kTH2F, {{200, 0.9, 1.2}, {200, -1, 1}});
    registryData.add("profileLambda", "Invariant Mass vs sin(phi)", {HistType::kTProfile, {{200, 0.9, 1.2}}});
    registryData.add("hLambdaPhiandSinPhi", "hLambdaPhiandSinPhi", kTH2F, {{200, -TMath::Pi() / 2, TMath::Pi() / 2}, {200, -1, 1}});
    registryData.add("V0LambdaprotonPhi", "V0LambdaprotonPhi", {HistType::kTH1F, {{200, -TMath::Pi() / 2, TMath::Pi() / 2}}});

    registryData.add("profileAntiLambda", "Invariant Mass vs sin(phi)", {HistType::kTProfile, {{200, 0.9, 1.2}}});
    registryData.add("TProfile1DLambdasinphiInJet", "#Delta #theta vs sin(phi)", {HistType::kTProfile, {{200, 0.0, TMath::Pi()}}});
    registryData.add("hAntiLambdamassandSinPhi", "hAntiLambdaPhiandSinPhi", kTH2F, {{200, -TMath::Pi() / 2, TMath::Pi() / 2}, {200, -1, 1}});
    registryData.add("hprotonsinphiInJetV0frame", "hprotonsinphiInJetV0frame", kTH1F, {axisSinPhi});
    registryData.add("TProfile2DLambdaPtMassSinPhi", "", kTProfile2D, {tprofile2DaxisMass, tprofile2DaxisPt});
    registryData.add("TProfile2DAntiLambdaPtMassSinPhi", "", kTProfile2D, {tprofile2DaxisMass, tprofile2DaxisPt});
    registryData.add("TProfile2DLambdaPtMassSintheta", "", kTProfile2D, {tprofile2DaxisMass, tprofile2DaxisPt});
    registryData.add("TProfile2DAntiLambdaPtMassSintheta", "", kTProfile2D, {tprofile2DaxisMass, tprofile2DaxisPt});

    registryData.add("TProfile2DLambdaPtMassCosSquareTheta", "", kTProfile2D, {tprofile2DaxisMass, tprofile2DaxisPt});
    registryData.add("TProfile2DAntiLambdaPtMassCosSquareTheta", "", kTProfile2D, {tprofile2DaxisMass, tprofile2DaxisPt});
    registryData.add("TProfile2DLambdaMassDeltaPhi", "", kTProfile2D, {{200, -TMath::Pi(), TMath::Pi(), "#Delta#varphi"}, tprofile2DaxisMass});
    registryData.add("TProfile2DLambdaMassDeltaTheta", "", kTProfile2D, {{200, 0, TMath::Pi(), "#Delta#theta"}, tprofile2DaxisMass});
    registryData.add("TProfile2DAntiLambdaMassDeltaPhi", "", kTProfile2D, {{200, -TMath::Pi(), TMath::Pi(), "#Delta#varphi"}, tprofile2DaxisMass});
    registryData.add("hprotonThetaInLab", "hprotonThetaInLab", kTH1F, {axisTheta});
    registryData.add("hprotonThetaInV0", "hprotonThetaInV0", kTH1F, {axisTheta});
    registryData.add("hprotonThetaInJetV0", "hprotonThetaInJetV0", kTH1F, {axisTheta});

    registryData.add("LambdaQA/TH2FLambdaMassPhiInJet", "TH2FLambdaMassPhiInJet", kTH2F, {{200, 0, TMath::Pi()}, {200, 0.9, 1.2}});
    registryData.add("LambdaQA/hArmenterosPreAnalyserCuts", "hArmenterosPreAnalyserCuts", kTH2F, {{1000, -1.0f, 1.0f, "#alpha"}, {1000, 0.0f, 0.30f, "#it{Q}_{T}"}});
    registryData.add("AntiLambdaQA/hArmenterosPreAnalyserCuts", "hArmenterosPreAnalyserCuts", kTH2F, {{1000, -1.0f, 1.0f, "#alpha"}, {1000, 0.0f, 0.30f, "#it{Q}_{T}"}});

    // Lab frame measures
    registryData.add("LambdaQA/TH2FprotonCosThetaInLab", "TH2FprotonCosThetaInLab", kTH2F, {{200, 0.9, 1.2}, {200, -1.0, 1.0}});
    registryData.add("LambdaQA/TProfile1DprotonCosThetaInLab", "TProfile1DprotonCosThetaInLab", {HistType::kTProfile, {{200, 0.9, 1.2}}});
    registryData.add("LambdaQA/TProfile1DprotonCos2ThetaInLab", "TProfile1DprotonCos2ThetaInLab", {HistType::kTProfile, {{200, 0.9, 1.2}}});

    // V0(Lambda) frame messages
    registryData.add("LambdaQA/TH2FprotonCosThetaInV0", "TH2FprotonCosThetaInV0", kTH2F, {{200, 0.9, 1.2}, {200, -1.0, 1.0}});
    registryData.add("LambdaQA/TProfile1DprotonCosThetaInV0", "TProfile1DprotonCosThetaInV0", {HistType::kTProfile, {{200, 0.9, 1.2}}});
    registryData.add("LambdaQA/TProfile1DprotonCos2ThetaInV0", "TProfile1DprotonCos2ThetaInV0", {HistType::kTProfile, {{200, 0.9, 1.2}}});

    // jet frame messages
    registryData.add("LambdaQA/TH2FprotonCosThetaInJet", "TH2FprotonCosThetaInJet", kTH2F, {{200, 0.9, 1.2}, {200, -1.0, 1.0}});
    registryData.add("LambdaQA/TProfile1DprotonCosThetaInJet", "TProfile1DprotonCosThetaInJet", {HistType::kTProfile, {{200, 0.9, 1.2}}});
    registryData.add("LambdaQA/TProfile1DprotonCos2ThetaInJet", "TProfile1DprotonCos2ThetaInJet", {HistType::kTProfile, {{200, 0.9, 1.2}}});

    // Jet-V0 frame messages
    registryData.add("LambdaQA/TH2FprotonCosThetaInJetV0", "TH2FprotonCosThetaInJetV0", kTH2F, {{200, 0.9, 1.2}, {200, -1.0, 1.0}});
    registryData.add("LambdaQA/TProfile1DprotonCosThetaInJetV0", "TProfile1DprotonCosThetaInJetV0", {HistType::kTProfile, {{200, 0.9, 1.2}}});
    registryData.add("LambdaQA/TProfile1DprotonCos2ThetaInJetV0", "TProfile1DprotonCos2ThetaInJetV0", {HistType::kTProfile, {{200, 0.9, 1.2}}});
    registryData.add("LambdaQA/TProfile2DprotonCosThetaInJetV0", "TProfile2DprotonCosThetaInJetV0", kTProfile2D, {tprofile2DaxisMass, axisDeltaPhi});
    registryData.add("LambdaQA/TProfile2DprotonCos2ThetaInJetV0", "TProfile2DprotonCos2ThetaInJetV0", kTProfile2D, {tprofile2DaxisMass, axisDeltaPhi});

    registryData.add("hNEvents", "hNEvents", {HistType::kTH1D, {{10, 0.f, 10.f}}});
    registryData.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(1, "all");
    registryData.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(2, "sel8");
    registryData.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(3, "TVX");
    registryData.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(4, "zvertex");
    registryData.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(5, "TFBorder");
    registryData.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(6, "ITSROFBorder");
    registryData.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(7, "isTOFVertexMatched");
    registryData.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(8, "isGoodZvtxFT0vsPV");
    registryData.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(9, "Applied selected");

    if (doQA) {
      registryData.add("QA/hv0sSelection", ";Sel", {HistType::kTH1D, {{22, 0., 22.}}});
      registryData.get<TH1>(HIST("QA/hv0sSelection"))->GetXaxis()->SetBinLabel(1, "all");
      registryData.get<TH1>(HIST("QA/hv0sSelection"))->GetXaxis()->SetBinLabel(2, "Event selection");
      registryData.get<TH1>(HIST("QA/hv0sSelection"))->GetXaxis()->SetBinLabel(3, "Radius");
      registryData.get<TH1>(HIST("QA/hv0sSelection"))->GetXaxis()->SetBinLabel(4, "Eta Daughters");
      registryData.get<TH1>(HIST("QA/hv0sSelection"))->GetXaxis()->SetBinLabel(5, "Dau DCA to PV");
      registryData.get<TH1>(HIST("QA/hv0sSelection"))->GetXaxis()->SetBinLabel(6, "DCA Daughters");
      registryData.get<TH1>(HIST("QA/hv0sSelection"))->GetXaxis()->SetBinLabel(7, "min ITS hits");
      registryData.get<TH1>(HIST("QA/hv0sSelection"))->GetXaxis()->SetBinLabel(8, "has TOF 1 Leg");
      registryData.get<TH1>(HIST("QA/hv0sSelection"))->GetXaxis()->SetBinLabel(9, "has TOF 2 Legs");
      registryData.get<TH1>(HIST("QA/hv0sSelection"))->GetXaxis()->SetBinLabel(10, "TPC NCl");
      registryData.get<TH1>(HIST("QA/hv0sSelection"))->GetXaxis()->SetBinLabel(11, "TPC Cls Shared");
      registryData.get<TH1>(HIST("QA/hv0sSelection"))->GetXaxis()->SetBinLabel(12, "ITS Chi2");
      registryData.get<TH1>(HIST("QA/hv0sSelection"))->GetXaxis()->SetBinLabel(13, "TPC Chi2");
      registryData.get<TH1>(HIST("QA/hv0sSelection"))->GetXaxis()->SetBinLabel(14, "cosPA");
      registryData.get<TH1>(HIST("QA/hv0sSelection"))->GetXaxis()->SetBinLabel(15, "rapidity");
      registryData.get<TH1>(HIST("QA/hv0sSelection"))->GetXaxis()->SetBinLabel(16, "ctau");
      registryData.get<TH1>(HIST("QA/hv0sSelection"))->GetXaxis()->SetBinLabel(17, "v0 rej");
      registryData.get<TH1>(HIST("QA/hv0sSelection"))->GetXaxis()->SetBinLabel(18, "TPC nsigma Neg Dau");
      registryData.get<TH1>(HIST("QA/hv0sSelection"))->GetXaxis()->SetBinLabel(19, "TPC nsigma Pos Dau");
      registryData.get<TH1>(HIST("QA/hv0sSelection"))->GetXaxis()->SetBinLabel(20, "Armenteros-Podolansky");
    }
  }
  double massPr = o2::constants::physics::MassProton;
  double massLambda = o2::constants::physics::MassLambda;
  double massPi = o2::constants::physics::MassPionCharged;
  ROOT::Math::PxPyPzMVector ProtonVec, PionVec, LambdaVec, ProtonBoostedVec, LambdaBoostedVec;

  TMatrixD LorentzTransInV0frame(double ELambda, double Lambdapx, double Lambdapy, double Lambdapz)
  {
    double PLambda = std::sqrt(Lambdapx * Lambdapx + Lambdapy * Lambdapy + Lambdapz * Lambdapz);
    double LambdaMass = std::sqrt(ELambda * ELambda - PLambda * PLambda);
    double Alpha = 1 / (LambdaMass * (ELambda + LambdaMass));
    TMatrixD matrixLabToLambda(4, 4);
    matrixLabToLambda(0, 0) = ELambda / LambdaMass;
    matrixLabToLambda(0, 1) = -Lambdapx / LambdaMass;
    matrixLabToLambda(0, 2) = -Lambdapy / LambdaMass;
    matrixLabToLambda(0, 3) = -Lambdapz / LambdaMass;
    matrixLabToLambda(1, 0) = -Lambdapx / LambdaMass;
    matrixLabToLambda(1, 1) = 1 + Alpha * Lambdapx * Lambdapx;
    matrixLabToLambda(1, 2) = Alpha * Lambdapx * Lambdapy;
    matrixLabToLambda(1, 3) = Alpha * Lambdapx * Lambdapz;
    matrixLabToLambda(2, 0) = -Lambdapy / LambdaMass;
    matrixLabToLambda(2, 1) = Alpha * Lambdapy * Lambdapx;
    matrixLabToLambda(2, 2) = 1 + Alpha * Lambdapy * Lambdapy;
    matrixLabToLambda(2, 3) = Alpha * Lambdapy * Lambdapz;
    matrixLabToLambda(3, 0) = -Lambdapz / LambdaMass;
    matrixLabToLambda(3, 1) = Alpha * Lambdapz * Lambdapx;
    matrixLabToLambda(3, 2) = Alpha * Lambdapz * Lambdapy;
    matrixLabToLambda(3, 3) = 1 + Alpha * Lambdapz * Lambdapz;
    return matrixLabToLambda;
  }
  // The direction of jet is z axis, y is perpendicular to jet and lambda momentum
  TMatrixD MyTMatrixTranslationToJet(double Jetpx, double Jetpy, double Jetpz, double Lambdapx, double Lambdapy, double Lambdapz)
  {
    TVector3 UnitX(1.0, 0.0, 0.0);
    TVector3 UnitY(0.0, 1.0, 0.0);
    TVector3 UnitZ(0.0, 0.0, 1.0);
    TVector3 JetP(Jetpx, Jetpy, Jetpz);
    TVector3 V0LambdaP(Lambdapx, Lambdapy, Lambdapz);
    TVector3 vortex_y = (JetP.Cross(V0LambdaP));

    TVector3 z_hat = JetP.Unit();
    TVector3 y_hat = vortex_y.Unit();
    TVector3 x_hat1 = y_hat.Cross(z_hat);
    TVector3 x_hat = x_hat1.Unit();

    TMatrixD matrixLabToJet(4, 4);
    matrixLabToJet(0, 0) = 1;
    matrixLabToJet(0, 1) = 0.0;
    matrixLabToJet(0, 2) = 0.0;
    matrixLabToJet(0, 3) = 0.0;
    matrixLabToJet(1, 0) = 0.0;
    matrixLabToJet(1, 1) = x_hat.X();
    matrixLabToJet(1, 2) = x_hat.Y();
    matrixLabToJet(1, 3) = x_hat.Z();
    matrixLabToJet(2, 0) = 0.0;
    matrixLabToJet(2, 1) = y_hat.X();
    matrixLabToJet(2, 2) = y_hat.Y();
    matrixLabToJet(2, 3) = y_hat.Z();
    matrixLabToJet(3, 0) = 0.0;
    matrixLabToJet(3, 1) = z_hat.X();
    matrixLabToJet(3, 2) = z_hat.Y();
    matrixLabToJet(3, 3) = z_hat.Z();
    return matrixLabToJet;
  }
  // New transformation: The direction of jet is x axis, z is perpendicular to jet and lambda momentum
  TMatrixD TMatrixTranslationToJet(double Jetpx, double Jetpy, double Jetpz, double Lambdapx, double Lambdapy, double Lambdapz)
  {
    TVector3 UnitX(1.0, 0.0, 0.0);
    TVector3 UnitY(0.0, 1.0, 0.0);
    TVector3 UnitZ(0.0, 0.0, 1.0);
    TVector3 JetP(Jetpx, Jetpy, Jetpz);
    TVector3 V0LambdaP(Lambdapx, Lambdapy, Lambdapz);
    TVector3 vortex_z = (JetP.Cross(V0LambdaP));

    TVector3 x_hat = JetP.Unit();
    TVector3 z_hat = vortex_z.Unit();
    TVector3 y_hat = z_hat.Cross(x_hat);

    TMatrixD matrixLabToJet(4, 4);
    matrixLabToJet(0, 0) = 1;
    matrixLabToJet(0, 1) = 0.0;
    matrixLabToJet(0, 2) = 0.0;
    matrixLabToJet(0, 3) = 0.0;
    matrixLabToJet(1, 0) = 0.0;
    matrixLabToJet(1, 1) = x_hat.X();
    matrixLabToJet(1, 2) = x_hat.Y();
    matrixLabToJet(1, 3) = x_hat.Z();
    matrixLabToJet(2, 0) = 0.0;
    matrixLabToJet(2, 1) = y_hat.X();
    matrixLabToJet(2, 2) = y_hat.Y();
    matrixLabToJet(2, 3) = y_hat.Z();
    matrixLabToJet(3, 0) = 0.0;
    matrixLabToJet(3, 1) = z_hat.X();
    matrixLabToJet(3, 2) = z_hat.Y();
    matrixLabToJet(3, 3) = z_hat.Z();
    return matrixLabToJet;
  }
  // aod::MyCollision const& collision

  // ITS hit
  template <typename TrackIts>
  bool hasITSHit(const TrackIts& track, int layer)
  {
    int ibit = layer - 1;
    return (track.itsClusterMap() & (1 << ibit));
  }

  // Single-Track Selection for Particles inside Jets
  template <typename JetTrack>
  bool passedTrackSelectionForJetReconstruction(const JetTrack& track)
  {
    const int minTpcCr = 70;
    const double minCrFindable = 0.8;
    const double maxChi2Tpc = 4.0;
    const double maxChi2Its = 36.0;
    const double maxPseudorapidity = 0.9;
    const double minPtTrack = 0.1;
    const double dcaxyMaxTrackPar0 = 0.0105;
    const double dcaxyMaxTrackPar1 = 0.035;
    const double dcaxyMaxTrackPar2 = 1.1;
    const double dcazMaxTrack = 2.0;

    if (!track.hasITS())
      return false;
    if ((!hasITSHit(track, 1)) && (!hasITSHit(track, 2)) && (!hasITSHit(track, 3)))
      return false;
    if (!track.hasTPC())
      return false;
    if (track.tpcNClsCrossedRows() < minTpcCr)
      return false;
    if ((static_cast<double>(track.tpcNClsCrossedRows()) / static_cast<double>(track.tpcNClsFindable())) < minCrFindable)
      return false;
    if (track.tpcChi2NCl() > maxChi2Tpc)
      return false;
    if (track.itsChi2NCl() > maxChi2Its)
      return false;
    if (track.eta() < -maxPseudorapidity || track.eta() > maxPseudorapidity)
      return false;
    if (track.pt() < minPtTrack)
      return false;
    if (std::fabs(track.dcaXY()) > (dcaxyMaxTrackPar0 + dcaxyMaxTrackPar1 / std::pow(track.pt(), dcaxyMaxTrackPar2)))
      return false;
    if (std::fabs(track.dcaZ()) > dcazMaxTrack)
      return false;
    return true;
  }

  // init Selection
  template <typename Lambda, typename TrackPos, typename TrackNeg>
  bool passedInitLambdaSelection(const Lambda& v0, const TrackPos& ptrack, const TrackNeg& ntrack)
  {
    if (v0.v0radius() < v0radius || v0.v0cosPA() < v0cospainit ||
        TMath::Abs(ptrack.eta()) > v0TracketaMax ||
        TMath::Abs(ntrack.eta()) > v0TracketaMax) {
      return false;
    }
    if (v0.dcaV0daughters() > dcav0dau) {
      return false;
    }
    if (TMath::Abs(v0.dcanegtopv()) < dcanegtopv) {
      return false;
    }
    if (TMath::Abs(v0.dcapostopv()) < dcapostopv) {
      return false;
    }
    return true;
  }

  template <typename Lambda, typename TrackPos, typename TrackNeg, typename TCollision>
  bool AcceptV0Lambda(const Lambda& v0, const TrackPos& ptrack, const TrackNeg& ntrack, const TCollision& collision)
  {
    // Single-Track Selections
    if (requirepassedSingleTrackSelection && !passedSingleTrackSelection(ptrack))
      return false;
    if (requirepassedSingleTrackSelection && !passedSingleTrackSelection(ntrack))
      return false;

    int evFlag = 0;
    if (collision.isInelGt0()) {
      evFlag = 1;
    }

    if (evSel && evFlag < 1)
      return false;

    if (v0.v0radius() < minimumV0Radius || v0.v0radius() > maximumV0Radius)
      return false;

    if (TMath::Abs(ptrack.eta()) > v0TracketaMax || TMath::Abs(ntrack.eta()) > v0TracketaMax) {
      return false;
    }

    if (std::fabs(v0.dcapostopv()) < dcapostoPVmin)
      return false;
    if (std::fabs(v0.dcanegtopv()) < dcanegtoPVmin)
      return false;

    if (std::fabs(v0.dcaV0daughters()) > dcaV0DaughtersMax)
      return false;

    if (requireITS && ptrack.itsNCls() < minITSnCls)
      return false;
    if (requireITS && ntrack.itsNCls() < minITSnCls)
      return false;

    if (hasTOF1Leg && !ptrack.hasTOF() && !ntrack.hasTOF())
      return false;

    if (hasTOF2Leg && (!ptrack.hasTOF() || !ntrack.hasTOF()))
      return false;

    if (ptrack.tpcNClsCrossedRows() < minNCrossedRowsTPC)
      return false;
    if (ntrack.tpcNClsCrossedRows() < minNCrossedRowsTPC)
      return false;

    if (ptrack.tpcNClsShared() > maxTPCSharedCls)
      return false;
    if (ntrack.tpcNClsShared() > maxTPCSharedCls)
      return false;

    if (ptrack.itsChi2NCl() > maxChi2ITS)
      return false;
    if (ntrack.itsChi2NCl() > maxChi2ITS)
      return false;

    if (ptrack.tpcChi2NCl() > maxChi2TPC)
      return false;
    if (ntrack.tpcChi2NCl() > maxChi2TPC)
      return false;

    if (v0.v0cosPA() < v0cospaMin)
      return false;

    if (v0.yLambda() < yMin || v0.yLambda() > yMax) {
      return false;
    }

    float ctauLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0;
    if (ctauLambda >= cTau)
      return false;

    if (TMath::Abs(v0.mK0Short() - o2::constants::physics::MassK0Short) < v0rejLambda) {
      return false;
    }
    if (std::abs(v0.mLambda() - o2::constants::physics::MassLambda0) > v0accLambda) {
      return false;
    }

    if (ntrack.tpcNSigmaPi() < nsigmaTPCmin || ntrack.tpcNSigmaPi() > nsigmaTPCmax)
      return false;

    if (ptrack.tpcNSigmaPr() < nsigmaTPCmin || ptrack.tpcNSigmaPr() > nsigmaTPCmax)
      return false;

    if (doArmenterosCut && v0.qtarm() > (paramArmenterosCut * std::abs(v0.alpha())))
      return false;

    return true;
  }

  template <typename Lambda, typename TrackPos, typename TrackNeg, typename TCollision>
  bool AcceptV0AntiLambda(const Lambda& v0, const TrackPos& ptrack, const TrackNeg& ntrack, const TCollision& collision)
  {
    // Single-Track Selections
    if (requirepassedSingleTrackSelection && !passedSingleTrackSelection(ptrack))
      return false;
    if (requirepassedSingleTrackSelection && !passedSingleTrackSelection(ntrack))
      return false;

    int evFlag = 0;
    if (collision.isInelGt0()) {
      evFlag = 1;
    }

    if (evSel && evFlag < 1)
      return false;

    if (v0.v0radius() < minimumV0Radius || v0.v0radius() > maximumV0Radius)
      return false;

    if (TMath::Abs(ptrack.eta()) > v0TracketaMax || TMath::Abs(ntrack.eta()) > v0TracketaMax) {
      return false;
    }

    if (std::fabs(v0.dcapostopv()) < dcapostoPVmin)
      return false;
    if (std::fabs(v0.dcanegtopv()) < dcanegtoPVmin)
      return false;

    if (std::fabs(v0.dcaV0daughters()) > dcaV0DaughtersMax)
      return false;

    if (requireITS && ptrack.itsNCls() < minITSnCls)
      return false;
    if (requireITS && ntrack.itsNCls() < minITSnCls)
      return false;

    if (hasTOF1Leg && !ptrack.hasTOF() && !ntrack.hasTOF())
      return false;

    if (hasTOF2Leg && (!ptrack.hasTOF() || !ntrack.hasTOF()))
      return false;

    if (ptrack.tpcNClsCrossedRows() < minNCrossedRowsTPC)
      return false;
    if (ntrack.tpcNClsCrossedRows() < minNCrossedRowsTPC)
      return false;

    if (ptrack.tpcNClsShared() > maxTPCSharedCls)
      return false;
    if (ntrack.tpcNClsShared() > maxTPCSharedCls)
      return false;

    if (ptrack.itsChi2NCl() > maxChi2ITS)
      return false;
    if (ntrack.itsChi2NCl() > maxChi2ITS)
      return false;

    if (ptrack.tpcChi2NCl() > maxChi2TPC)
      return false;
    if (ntrack.tpcChi2NCl() > maxChi2TPC)
      return false;

    if (v0.v0cosPA() < v0cospaMin)
      return false;

    if (v0.yLambda() < yMin || v0.yLambda() > yMax) {
      return false;
    }

    float ctauAntiLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0Bar;
    if (ctauAntiLambda >= cTau)
      return false;

    if (TMath::Abs(v0.mK0Short() - o2::constants::physics::MassK0Short) < v0rejLambda) {
      return false;
    }
    if (std::abs(v0.mAntiLambda() - o2::constants::physics::MassLambda0) > v0accLambda) {
      return false;
    }
    if (ntrack.tpcNSigmaPr() < nsigmaTPCmin || ntrack.tpcNSigmaPr() > nsigmaTPCmax)
      return false;

    if (ptrack.tpcNSigmaPi() < nsigmaTPCmin || ptrack.tpcNSigmaPi() > nsigmaTPCmax)
      return false;

    if (doArmenterosCut && v0.qtarm() > (paramArmenterosCut * std::abs(v0.alpha())))
      return false;

    return true;
  }

  template <typename Lambda, typename TrackPos, typename TrackNeg, typename TCollision>
  bool registryDataAcceptV0Lambda(const Lambda& v0, const TrackPos& ptrack, const TrackNeg& ntrack, const TCollision& collision)
  {

    int evFlag = 0;
    if (collision.isInelGt0()) {
      evFlag = 1;
    }

    if (v0.pt() < v0Ptmin) {
      return false;
    }

    registryData.fill(HIST("QA/hv0sSelection"), 0.5);

    if (evSel && evFlag < 1)
      return false;

    registryData.fill(HIST("QA/hv0sSelection"), 1.5);

    if (v0.v0radius() < minimumV0Radius || v0.v0radius() > maximumV0Radius)
      return false;

    registryData.fill(HIST("QA/hv0sSelection"), 2.5);

    if (TMath::Abs(ptrack.eta()) > v0TracketaMax || TMath::Abs(ntrack.eta()) > v0TracketaMax) {
      return false;
    }
    registryData.fill(HIST("QA/hv0sSelection"), 3.5);

    if (std::fabs(v0.dcapostopv()) < dcapostoPVmin)
      return false;
    if (std::fabs(v0.dcanegtopv()) < dcanegtoPVmin)
      return false;
    registryData.fill(HIST("QA/hv0sSelection"), 4.5);

    if (std::fabs(v0.dcaV0daughters()) > dcaV0DaughtersMax)
      return false;
    registryData.fill(HIST("QA/hv0sSelection"), 5.5);

    if (requireITS && ptrack.itsNCls() < minITSnCls)
      return false;
    if (requireITS && ntrack.itsNCls() < minITSnCls)
      return false;
    registryData.fill(HIST("QA/hv0sSelection"), 6.5);

    if (hasTOF1Leg && !ptrack.hasTOF() && !ntrack.hasTOF())
      return false;
    registryData.fill(HIST("QA/hv0sSelection"), 7.5);

    if (hasTOF2Leg && (!ptrack.hasTOF() || !ntrack.hasTOF()))
      return false;
    registryData.fill(HIST("QA/hv0sSelection"), 8.5);

    if (ptrack.tpcNClsCrossedRows() < minNCrossedRowsTPC)
      return false;
    if (ntrack.tpcNClsCrossedRows() < minNCrossedRowsTPC)
      return false;
    registryData.fill(HIST("QA/hv0sSelection"), 9.5);

    if (requireMaxTPCSharedCls && ptrack.tpcNClsShared() > maxTPCSharedCls)
      return false;
    if (requireMaxTPCSharedCls && ntrack.tpcNClsShared() > maxTPCSharedCls)
      return false;
    registryData.fill(HIST("QA/hv0sSelection"), 10.5);

    if (ptrack.itsChi2NCl() > maxChi2ITS)
      return false;
    if (ntrack.itsChi2NCl() > maxChi2ITS)
      return false;
    registryData.fill(HIST("QA/hv0sSelection"), 11.5);

    if (ptrack.tpcChi2NCl() > maxChi2TPC)
      return false;
    if (ntrack.tpcChi2NCl() > maxChi2TPC)
      return false;
    registryData.fill(HIST("QA/hv0sSelection"), 12.5);

    if (v0.v0cosPA() < v0cospaMin)
      return false;
    registryData.fill(HIST("QA/hv0sSelection"), 13.5);

    if (v0.yLambda() < yMin || v0.yLambda() > yMax) {
      return false;
    }
    registryData.fill(HIST("QA/hv0sSelection"), 14.5);

    float ctauLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0;
    if (ctauLambda >= cTau)
      return false;
    registryData.fill(HIST("QA/hv0sSelection"), 15.5);

    if (TMath::Abs(v0.mK0Short() - o2::constants::physics::MassK0Short) < v0rejLambda) {
      return false;
    }
    if (std::abs(v0.mLambda() - o2::constants::physics::MassLambda0) > v0accLambda) {
      return false;
    }

    registryData.fill(HIST("QA/hv0sSelection"), 16.5);

    if (ntrack.tpcNSigmaPi() < nsigmaTPCmin || ntrack.tpcNSigmaPi() > nsigmaTPCmax)
      return false;
    registryData.fill(HIST("QA/hv0sSelection"), 17.5);

    if (ptrack.tpcNSigmaPr() < nsigmaTPCmin || ptrack.tpcNSigmaPr() > nsigmaTPCmax)
      return false;
    registryData.fill(HIST("QA/hv0sSelection"), 18.5);

    if (doArmenterosCut && v0.qtarm() > (paramArmenterosCut * std::abs(v0.alpha())))
      return false;
    registryData.fill(HIST("QA/hv0sSelection"), 19.5);

    return true;
  }

  template <typename Lambda, typename TrackPos, typename TrackNeg, typename TCollision>
  bool registryDataAcceptV0AntiLambda(const Lambda& v0, const TrackPos& ptrack, const TrackNeg& ntrack, const TCollision& collision)
  {

    int evFlag = 0;
    if (collision.isInelGt0()) {
      evFlag = 1;
    }
    if (v0.pt() < v0Ptmin) {
      return false;
    }
    if (evSel && evFlag < 1)
      return false;

    if (v0.v0radius() < minimumV0Radius || v0.v0radius() > maximumV0Radius)
      return false;

    if (TMath::Abs(ptrack.eta()) > v0TracketaMax || TMath::Abs(ntrack.eta()) > v0TracketaMax) {
      return false;
    }

    if (std::fabs(v0.dcapostopv()) < dcapostoPVmin)
      return false;
    if (std::fabs(v0.dcanegtopv()) < dcanegtoPVmin)
      return false;

    if (std::fabs(v0.dcaV0daughters()) > dcaV0DaughtersMax)
      return false;

    if (requireITS && ptrack.itsNCls() < minITSnCls)
      return false;
    if (requireITS && ntrack.itsNCls() < minITSnCls)
      return false;

    if (hasTOF1Leg && !ptrack.hasTOF() && !ntrack.hasTOF())
      return false;

    if (hasTOF2Leg && (!ptrack.hasTOF() || !ntrack.hasTOF()))
      return false;

    if (ptrack.tpcNClsCrossedRows() < minNCrossedRowsTPC)
      return false;
    if (ntrack.tpcNClsCrossedRows() < minNCrossedRowsTPC)
      return false;

    if (ptrack.tpcNClsShared() > maxTPCSharedCls)
      return false;
    if (ntrack.tpcNClsShared() > maxTPCSharedCls)
      return false;

    if (ptrack.itsChi2NCl() > maxChi2ITS)
      return false;
    if (ntrack.itsChi2NCl() > maxChi2ITS)
      return false;

    if (ptrack.tpcChi2NCl() > maxChi2TPC)
      return false;
    if (ntrack.tpcChi2NCl() > maxChi2TPC)
      return false;

    if (v0.v0cosPA() < v0cospaMin)
      return false;

    if (v0.yLambda() < yMin || v0.yLambda() > yMax) {
      return false;
    }

    float ctauAntiLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0Bar;
    if (ctauAntiLambda >= cTau)
      return false;

    if (TMath::Abs(v0.mK0Short() - o2::constants::physics::MassK0Short) < v0rejLambda) {
      return false;
    }
    if (std::abs(v0.mAntiLambda() - o2::constants::physics::MassLambda0) > v0accLambda) {
      return false;
    }
    if (ntrack.tpcNSigmaPr() < nsigmaTPCmin || ntrack.tpcNSigmaPr() > nsigmaTPCmax)
      return false;

    if (ptrack.tpcNSigmaPi() < nsigmaTPCmin || ptrack.tpcNSigmaPi() > nsigmaTPCmax)
      return false;

    if (doArmenterosCut && v0.qtarm() > (paramArmenterosCut * std::abs(v0.alpha())))
      return false;

    return true;
  }

  // Lambda Selections
  template <typename Lambda, typename TrackPos, typename TrackNeg>
  bool passedLambdaSelection(const Lambda& v0, const TrackPos& ptrack, const TrackNeg& ntrack)
  {
    // Single-Track Selections
    if (requirepassedSingleTrackSelection && !passedSingleTrackSelection(ptrack))
      return false;
    if (requirepassedSingleTrackSelection && !passedSingleTrackSelection(ntrack))
      return false;

    if (v0.v0radius() < minimumV0Radius || v0.v0radius() > maximumV0Radius)
      return false;

    if (TMath::Abs(ptrack.eta()) > v0TracketaMax || TMath::Abs(ntrack.eta()) > v0TracketaMax) {
      return false;
    }

    if (std::fabs(v0.dcapostopv()) < dcapostoPVmin)
      return false;
    if (std::fabs(v0.dcanegtopv()) < dcanegtoPVmin)
      return false;

    if (std::fabs(v0.dcaV0daughters()) > dcaV0DaughtersMax)
      return false;

    if (requireITS && ptrack.itsNCls() < minITSnCls)
      return false;
    if (requireITS && ntrack.itsNCls() < minITSnCls)
      return false;

    if (hasTOF1Leg && !ptrack.hasTOF() && !ntrack.hasTOF())
      return false;

    if (hasTOF2Leg && (!ptrack.hasTOF() || !ntrack.hasTOF()))
      return false;

    if (ptrack.tpcNClsCrossedRows() < minNCrossedRowsTPC)
      return false;
    if (ntrack.tpcNClsCrossedRows() < minNCrossedRowsTPC)
      return false;

    if (ptrack.tpcNClsShared() > maxTPCSharedCls)
      return false;
    if (ntrack.tpcNClsShared() > maxTPCSharedCls)
      return false;

    if (ptrack.itsChi2NCl() > maxChi2ITS)
      return false;
    if (ntrack.itsChi2NCl() > maxChi2ITS)
      return false;

    if (ptrack.tpcChi2NCl() > maxChi2TPC)
      return false;
    if (ntrack.tpcChi2NCl() > maxChi2TPC)
      return false;

    if (v0.v0cosPA() < v0cospaMin)
      return false;

    if (v0.yLambda() < yMin || v0.yLambda() > yMax) {
      return false;
    }

    // PID Selections (TPC)
    if (requireTPC) {
      if (ptrack.tpcNSigmaPr() < nsigmaTPCmin || ptrack.tpcNSigmaPr() > nsigmaTPCmax)
        return false;
      if (ntrack.tpcNSigmaPi() < nsigmaTPCmin || ntrack.tpcNSigmaPi() > nsigmaTPCmax)
        return false;
    }
    // PID Selections (TOF)
    if (requireTOF) {
      if (ptrack.tofNSigmaPr() < nsigmaTOFmin || ptrack.tofNSigmaPr() > nsigmaTOFmax)
        return false;
      if (ntrack.tofNSigmaPi() < nsigmaTOFmin || ntrack.tofNSigmaPi() > nsigmaTOFmax)
        return false;
    }
    // TLorentzVector lorentzVect;
    // lorentzVect.SetXYZM(v0.px(), v0.py(), v0.pz(), 1.115683);
    ROOT::Math::PxPyPzMVector lorentzVect(v0.px(), v0.py(), v0.pz(), o2::constants::physics::MassLambda0);
    if (lorentzVect.Rapidity() < yMin || lorentzVect.Rapidity() > yMax) {
      return false;
    }

    if (TMath::Abs(v0.mK0Short() - o2::constants::physics::MassK0Short) < v0rejLambda) {
      return false;
    }
    if (TMath::Abs(v0.mLambda() - o2::constants::physics::MassLambda0) > v0accLambda) {
      return false;
    }
    if (doArmenterosCut && v0.qtarm() > (paramArmenterosCut * std::abs(v0.alpha())))
      return false;

    return true;
  }

  // AntiLambda Selections
  template <typename AntiLambda, typename TrackPos, typename TrackNeg>
  bool passedAntiLambdaSelection(const AntiLambda& v0, const TrackPos& ptrack, const TrackNeg& ntrack)
  {
    // Single-Track Selections
    if (requirepassedSingleTrackSelection && !passedSingleTrackSelection(ptrack))
      return false;
    if (requirepassedSingleTrackSelection && !passedSingleTrackSelection(ntrack))
      return false;

    // Momentum AntiLambda Daughters
    TVector3 pion(v0.pxpos(), v0.pypos(), v0.pzpos());
    TVector3 proton(v0.pxneg(), v0.pyneg(), v0.pzneg());

    if (proton.Pt() < ptMinV0Proton)
      return false;
    if (proton.Pt() > ptMaxV0Proton)
      return false;
    if (pion.Pt() < ptMinV0Pion)
      return false;
    if (pion.Pt() > ptMaxV0Pion)
      return false;

    // V0 Selections
    if (v0.v0cosPA() < v0cospaMin)
      return false;
    if (v0.v0radius() < minimumV0Radius || v0.v0radius() > maximumV0Radius)
      return false;
    if (std::fabs(v0.dcaV0daughters()) > dcaV0DaughtersMax)
      return false;
    if (std::fabs(v0.dcapostopv()) < dcapostoPVmin)
      return false;
    if (std::fabs(v0.dcanegtopv()) < dcanegtoPVmin)
      return false;

    if (TMath::Abs(ptrack.eta()) > v0TracketaMax || TMath::Abs(ntrack.eta()) > v0TracketaMax) {
      return false;
    }

    // PID Selections (TPC)
    if (requireTPC) {
      if (ptrack.tpcNSigmaPi() < nsigmaTPCmin || ptrack.tpcNSigmaPi() > nsigmaTPCmax)
        return false;
      if (ntrack.tpcNSigmaPr() < nsigmaTPCmin || ntrack.tpcNSigmaPr() > nsigmaTPCmax)
        return false;
    }
    // PID Selections (TOF)
    if (requireTOF) {
      if (ptrack.tofNSigmaPi() < nsigmaTOFmin || ptrack.tofNSigmaPi() > nsigmaTOFmax)
        return false;
      if (ntrack.tofNSigmaPr() < nsigmaTOFmin || ntrack.tofNSigmaPr() > nsigmaTOFmax)
        return false;
    }
    // TLorentzVector lorentzVect;
    // lorentzVect.SetXYZM(v0.px(), v0.py(), v0.pz(), 1.115683);
    ROOT::Math::PxPyPzMVector lorentzVect(v0.px(), v0.py(), v0.pz(), o2::constants::physics::MassLambda0);

    if (lorentzVect.Rapidity() < yMin || lorentzVect.Rapidity() > yMax) {
      return false;
    }

    if (TMath::Abs(v0.mK0Short() - o2::constants::physics::MassK0Short) < v0rejLambda) {
      return false;
    }
    if (TMath::Abs(v0.mAntiLambda() - o2::constants::physics::MassLambda0) > v0accLambda) {
      return false;
    }
    if (doArmenterosCut && v0.qtarm() > (paramArmenterosCut * std::abs(v0.alpha())))
      return false;
    return true;
  }

  // Single-Track Selection
  template <typename Track>
  bool passedSingleTrackSelection(const Track& track)
  {
    if (requireITS && (!track.hasITS()))
      return false;
    if (requireITS && track.itsNCls() < minITSnCls)
      return false;
    if (!track.hasTPC())
      return false;
    if (track.tpcNClsFound() < minTPCnClsFound)
      return false;
    if (track.tpcNClsCrossedRows() < minNCrossedRowsTPC)
      return false;
    if ((static_cast<float>(track.tpcNClsCrossedRows()) / static_cast<float>(track.tpcNClsFindable())) < minTpcNcrossedRowsOverFindable)
      return false;
    if (track.tpcChi2NCl() > maxChi2TPC)
      return false;
    if (track.eta() < etaMin || track.eta() > etaMax)
      return false;
    if (requireTOF && (!track.hasTOF()))
      return false;
    return true;
  }

  ///////Event selection
  template <typename TCollision>
  bool AcceptEvent(TCollision const& collision)
  {
    if (sel8 && !collision.sel8()) {
      return false;
    }
    registryData.fill(HIST("hNEvents"), 1.5);

    if (isTriggerTVX && !collision.selection_bit(aod::evsel::kIsTriggerTVX)) {
      return false;
    }
    registryData.fill(HIST("hNEvents"), 2.5);

    if (iscutzvertex && TMath::Abs(collision.posZ()) > cutzvertex) {
      return false;
    }
    registryData.fill(HIST("hNEvents"), 3.5);

    if (isNoTimeFrameBorder && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }

    registryData.fill(HIST("hNEvents"), 4.5);

    if (isNoITSROFrameBorder && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    registryData.fill(HIST("hNEvents"), 5.5);
    if (isVertexTOFmatched && !collision.selection_bit(aod::evsel::kIsVertexTOFmatched)) {
      return false;
    }
    registryData.fill(HIST("hNEvents"), 6.5);
    if (isGoodZvtxFT0vsPV && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    registryData.fill(HIST("hNEvents"), 7.5);

    return true;
  }

  template <typename TCollision>
  bool AcceptEventForLongitudinalPolarization(TCollision const& collision)
  {

    if (sel8 && !collision.sel8()) {
      return false;
    }

    if (isTriggerTVX && !collision.selection_bit(aod::evsel::kIsTriggerTVX)) {
      return false;
    }

    if (std::abs(collision.posZ()) > cutzvertex) {
      return false;
    }

    if (isNoTimeFrameBorder && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }

    if (isNoITSROFrameBorder && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }

    if (isVertexTOFmatched && !collision.selection_bit(aod::evsel::kIsVertexTOFmatched)) {
      return false;
    }

    if (isNoSameBunchPileup && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return false;
    }

    return true;
  }

  using SelCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;
  using SelV0Collisions = soa::Join<aod::Collisions, aod::EvSels, aod::PVMults, aod::CentFT0Ms, aod::CentNGlobals>;
  using StrHadronDaughterTracks = soa::Join<aod::Tracks, aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
  void processData(SelV0Collisions::iterator const& collision, aod::V0Datas const& fullV0s, StrHadronDaughterTracks const& tracks)
  {
    registryData.fill(HIST("hNEvents"), 0.5);
    if (!AcceptEvent(collision)) {
      return;
    }
    registryData.fill(HIST("hNEvents"), 8.5);
    // event selection
    // loop over reconstructed tracks
    std::vector<fastjet::PseudoJet> fjParticles;
    for (auto const& track : tracks) {
      registryData.fill(HIST("h_track_pt"), track.pt());
      registryData.fill(HIST("h_track_eta"), track.eta());
      registryData.fill(HIST("h_track_phi"), track.phi());
      if (ispassdTrackSelectionForJetReconstruction && !passedTrackSelectionForJetReconstruction(track)) {
        continue;
      }
      registryData.fill(HIST("h_track_pt_sel"), track.pt());
      registryData.fill(HIST("h_track_eta_sel"), track.eta());
      registryData.fill(HIST("h_track_phi_sel"), track.phi());

      // 4-momentum representation of a particle
      fastjet::PseudoJet fourMomentum(track.px(), track.py(), track.pz(), track.energy(o2::constants::physics::MassPionCharged));
      fjParticles.emplace_back(fourMomentum);
    }
    // reject empty events
    if (fjParticles.size() < 1)
      return;
    // cluster particles using the anti-kt algorithm
    fastjet::RecombinationScheme recombScheme = fastjet::E_scheme;
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, rJet, recombScheme);
    fastjet::AreaDefinition areaDef(fastjet::active_area, fastjet::GhostedAreaSpec(1.0));
    fastjet::ClusterSequenceArea cs(fjParticles, jetDef, areaDef);
    std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());

    // jet selection
    bool isAtLeastOneJetSelected = false;
    int nJets = 0;
    int nJetssel = 0;
    // select most large momentum jet
    float maxJetpx = 0;
    float maxJetpy = 0;
    float maxJetpz = 0;
    float maxJeteta = 0;
    float maxJetphi = 0;
    float maxJetE = 0;
    float maxJetpT = 0;
    float maxJetPt = -999;
    for (const auto& jet : jets) {
      nJets++;
      registryData.fill(HIST("FJetaHistogram"), jet.eta());
      registryData.fill(HIST("FJphiHistogram"), jet.phi());
      registryData.fill(HIST("FJptHistogram"), jet.pt());
      // jet must be fully contained in the acceptance
      if ((std::fabs(jet.eta()) + rJet) > (etaMax - deltaEtaEdge)) {
        continue;
      }

      if (jet.pt() < cfgjetPtMin)
        continue;
      nJetssel++;
      registryData.fill(HIST("FJetaHistogramsel"), jet.eta());
      registryData.fill(HIST("FJphiHistogramsel"), jet.phi());
      registryData.fill(HIST("FJptHistogramsel"), jet.pt());

      if (jet.pt() > maxJetPt) {
        maxJetpx = jet.px();
        maxJetpy = jet.py();
        maxJetpz = jet.pz();
        maxJeteta = jet.eta();
        maxJetE = jet.E();
        maxJetphi = jet.phi();
        maxJetpT = jet.pt();
        maxJetPt = maxJetpT;
      }
    }
    if (maxJetpT > 0) {
      registryData.fill(HIST("FLeadingJetaHistogramsel"), maxJeteta);
      registryData.fill(HIST("FLeadingJphiHistogramsel"), maxJetphi);
      registryData.fill(HIST("FLeadingJptHistogramsel"), maxJetpT);
    }
    registryData.fill(HIST("nJetsPerEvent"), nJets);
    registryData.fill(HIST("nJetsPerEventsel"), nJetssel);
    isAtLeastOneJetSelected = true;
    if (!isAtLeastOneJetSelected) {
      return;
    }

    // Event multiplicity
    const float multiplicity = collision.centFT0M();
    registryData.fill(HIST("number_of_events_vsmultiplicity"), multiplicity);
    // v0 loop
    int V0Numbers = 0;
    int AntiV0Numbers = 0;
    for (const auto& v0 : fullV0s) {
      const auto& pos = v0.posTrack_as<StrHadronDaughterTracks>();
      const auto& neg = v0.negTrack_as<StrHadronDaughterTracks>();
      TVector3 v0dir(v0.px(), v0.py(), v0.pz());
      if (registryDataAcceptV0Lambda(v0, pos, neg, collision)) {
        V0Numbers = V0Numbers + 1;
        registryData.fill(HIST("LambdaPtMass"), v0.pt(), v0.mLambda());
      }
      if (registryDataAcceptV0AntiLambda(v0, pos, neg, collision)) {
        AntiV0Numbers = AntiV0Numbers + 1;
        registryData.fill(HIST("AntiLambdaPtMass"), v0.pt(), v0.mAntiLambda());
      }
    }
    registryData.fill(HIST("nV0sPerEvent"), V0Numbers);

    // calculate lambda polarization induced by jet

    if (V0Numbers == 0 && AntiV0Numbers == 0) {
      return;
    }
    if (maxJetpx == 0) {
      return;
    }
    double protonsinPhiInJetV0frame = 0;
    double AntiprotonsinPhiInJetV0frame = 0;
    for (const auto& candidate : fullV0s) {
      const auto& pos = candidate.posTrack_as<StrHadronDaughterTracks>();
      const auto& neg = candidate.negTrack_as<StrHadronDaughterTracks>();
      TVector3 v0dir(candidate.px(), candidate.py(), candidate.pz());

      if (registryDataAcceptV0Lambda(candidate, pos, neg, collision)) {
        registryData.fill(HIST("hMassLambda"), candidate.mLambda());
        registryData.fill(HIST("V0pTInLab"), candidate.pt());
        registryData.fill(HIST("V0pxInLab"), candidate.px());
        registryData.fill(HIST("V0pyInLab"), candidate.py());
        registryData.fill(HIST("V0pzInLab"), candidate.pz());
        registryData.fill(HIST("protonQA/V0protonpxInLab"), pos.px());
        registryData.fill(HIST("protonQA/V0protonpyInLab"), pos.py());
        registryData.fill(HIST("protonQA/V0protonpzInLab"), pos.pz());

        double PLambda = std::sqrt(candidate.px() * candidate.px() + candidate.py() * candidate.py() + candidate.pz() * candidate.pz());
        double ELambda = std::sqrt(candidate.mLambda() * candidate.mLambda() + PLambda * PLambda);
        double protonE = std::sqrt(massPr * massPr + pos.px() * pos.px() + pos.py() * pos.py() + pos.pz() * pos.pz());

        TMatrixD pLabJet(4, 1);
        pLabJet(0, 0) = maxJetE;
        pLabJet(1, 0) = maxJetpx;
        pLabJet(2, 0) = maxJetpy;
        pLabJet(3, 0) = maxJetpz;

        TMatrixD pLabV0(4, 1);
        pLabV0(0, 0) = ELambda;
        pLabV0(1, 0) = candidate.px();
        pLabV0(2, 0) = candidate.py();
        pLabV0(3, 0) = candidate.pz();

        TMatrixD V0InV0(4, 1);
        V0InV0 = LorentzTransInV0frame(ELambda, candidate.px(), candidate.py(), candidate.pz()) * pLabV0;
        registryData.fill(HIST("V0pxInRest_frame"), V0InV0(1, 0));
        registryData.fill(HIST("V0pyInRest_frame"), V0InV0(2, 0));
        registryData.fill(HIST("V0pzInRest_frame"), V0InV0(3, 0));

        TMatrixD lambdaInJet(4, 1);
        lambdaInJet = MyTMatrixTranslationToJet(maxJetpx, maxJetpy, maxJetpz, candidate.px(), candidate.py(), candidate.pz()) * pLabV0;
        double cosThetaLambdaInJet = lambdaInJet(3, 0) / std::sqrt(lambdaInJet(1, 0) * lambdaInJet(1, 0) + lambdaInJet(2, 0) * lambdaInJet(2, 0) + lambdaInJet(3, 0) * lambdaInJet(3, 0));
        double lambdasinphiInJet = lambdaInJet(2, 0) / std::sqrt(lambdaInJet(1, 0) * lambdaInJet(1, 0) + lambdaInJet(2, 0) * lambdaInJet(2, 0));
        registryData.fill(HIST("TProfile2DLambdaMassDeltaTheta"), TMath::ACos(cosThetaLambdaInJet), candidate.mLambda(), lambdasinphiInJet);
        registryData.fill(HIST("TProfile1DLambdasinphiInJet"), TMath::ACos(cosThetaLambdaInJet), lambdasinphiInJet);

        registryData.fill(HIST("V0pxInJetframe"), lambdaInJet(1, 0));
        registryData.fill(HIST("V0pyInJetframe"), lambdaInJet(2, 0));
        registryData.fill(HIST("V0pzInJetframe"), lambdaInJet(3, 0));

        TMatrixD lambdaInJetV0(4, 1);
        lambdaInJetV0 = LorentzTransInV0frame(ELambda, lambdaInJet(1, 0), lambdaInJet(2, 0), lambdaInJet(3, 0)) * MyTMatrixTranslationToJet(maxJetpx, maxJetpy, maxJetpz, candidate.px(), candidate.py(), candidate.pz()) * pLabV0;
        registryData.fill(HIST("V0LambdapxInJetV0frame"), lambdaInJetV0(1, 0));
        registryData.fill(HIST("V0LambdapyInJetV0frame"), lambdaInJetV0(2, 0));
        registryData.fill(HIST("V0LambdapzInJetV0frame"), lambdaInJetV0(3, 0));

        TMatrixD pLabproton(4, 1);
        pLabproton(0, 0) = protonE;
        pLabproton(1, 0) = pos.px();
        pLabproton(2, 0) = pos.py();
        pLabproton(3, 0) = pos.pz();
        double protonsinPhiInLab = pLabproton(2, 0) / std::sqrt(pLabproton(1, 0) * pLabproton(1, 0) + pLabproton(2, 0) * pLabproton(2, 0));
        double protoncosthetaInLab = pLabproton(3, 0) / std::sqrt(pLabproton(1, 0) * pLabproton(1, 0) + pLabproton(2, 0) * pLabproton(2, 0) + pLabproton(3, 0) * pLabproton(3, 0));
        double protonPtInLab = std::sqrt(pLabproton(1, 0) * pLabproton(1, 0) + pLabproton(2, 0) * pLabproton(2, 0));
        double protonPInLab = std::sqrt(pLabproton(1, 0) * pLabproton(1, 0) + pLabproton(2, 0) * pLabproton(2, 0) + pLabproton(3, 0) * pLabproton(3, 0));
        double protonsinThetaInLab = protonPtInLab / protonPInLab;
        double protonMassInLab = std::sqrt(pLabproton(0, 0) * pLabproton(0, 0) - pLabproton(1, 0) * pLabproton(1, 0) - pLabproton(2, 0) * pLabproton(2, 0) - pLabproton(3, 0) * pLabproton(3, 0));
        double jettheta = maxJetpz / std::sqrt(pLabJet(1, 0) * pLabJet(1, 0) + pLabJet(2, 0) * pLabJet(2, 0) + pLabJet(3, 0) * pLabJet(3, 0));
        double jetphi = maxJetpy / std::sqrt(pLabJet(1, 0) * pLabJet(1, 0) + pLabJet(2, 0) * pLabJet(2, 0));
        double jetptInLab = std::sqrt(pLabJet(1, 0) * pLabJet(1, 0) + pLabJet(2, 0) * pLabJet(2, 0));
        registryData.fill(HIST("JetQA/JetthetaInLab"), TMath::ASin(jettheta));
        registryData.fill(HIST("JetQA/JetphiInLab"), TMath::ASin(jetphi));
        registryData.fill(HIST("JetQA/JetpxInLab"), pLabJet(1, 0));
        registryData.fill(HIST("JetQA/JetpyInLab"), pLabJet(2, 0));
        registryData.fill(HIST("JetQA/JetpzInLab"), pLabJet(3, 0));
        registryData.fill(HIST("JetQA/JetptInLab"), jetptInLab);

        registryData.fill(HIST("protonQA/V0protonphiInLab"), TMath::ASin(protonsinPhiInLab));
        registryData.fill(HIST("protonQA/V0protonthetaInLab"), TMath::ACos(protoncosthetaInLab));
        registryData.fill(HIST("protonQA/V0protoncosthetaInLab"), protoncosthetaInLab);
        registryData.fill(HIST("protonQA/profileprotonsinthetaInLab"), candidate.mLambda(), protonsinThetaInLab);
        registryData.fill(HIST("protonQA/profileprotonsinphiInLab"), candidate.mLambda(), protonsinPhiInLab);
        registryData.fill(HIST("protonQA/profileprotoncosSquarethetaInLab"), candidate.mLambda(), protoncosthetaInLab * protoncosthetaInLab);
        registryData.fill(HIST("protonQA/V0protonMassInLab"), protonMassInLab);

        TMatrixD protonInV0(4, 1);
        protonInV0 = LorentzTransInV0frame(ELambda, candidate.px(), candidate.py(), candidate.pz()) * pLabproton;
        double protonMassInV0 = std::sqrt(protonInV0(0, 0) * protonInV0(0, 0) - protonInV0(1, 0) * protonInV0(1, 0) - protonInV0(2, 0) * protonInV0(2, 0) - protonInV0(3, 0) * protonInV0(3, 0));
        double protonPInV0 = std::sqrt(protonInV0(1, 0) * protonInV0(1, 0) + protonInV0(2, 0) * protonInV0(2, 0) + protonInV0(3, 0) * protonInV0(3, 0));
        double protonPtInV0 = std::sqrt(protonInV0(1, 0) * protonInV0(1, 0) + protonInV0(2, 0) * protonInV0(2, 0));
        double protonsinThetaInV0 = protonPtInV0 / protonPInV0;

        TMatrixD JetInV0(4, 1);
        JetInV0 = LorentzTransInV0frame(ELambda, candidate.px(), candidate.py(), candidate.pz()) * pLabJet;
        double jetthetaInV0 = JetInV0(3, 0) / std::sqrt(JetInV0(1, 0) * JetInV0(1, 0) + JetInV0(2, 0) * JetInV0(2, 0) + JetInV0(3, 0) * JetInV0(3, 0));
        double jetphiInV0 = JetInV0(2, 0) / std::sqrt(JetInV0(1, 0) * JetInV0(1, 0) + JetInV0(2, 0) * JetInV0(2, 0));
        double jetptInV0 = std::sqrt(JetInV0(1, 0) * JetInV0(1, 0) + JetInV0(2, 0) * JetInV0(2, 0));
        registryData.fill(HIST("JetQA/JetthetaInV0"), TMath::ASin(jetthetaInV0));
        registryData.fill(HIST("JetQA/JetphiInV0"), TMath::ASin(jetphiInV0));
        registryData.fill(HIST("JetQA/JetpxInV0"), JetInV0(1, 0));
        registryData.fill(HIST("JetQA/JetpyInV0"), JetInV0(2, 0));
        registryData.fill(HIST("JetQA/JetpzInV0"), JetInV0(3, 0));
        registryData.fill(HIST("JetQA/JetptInV0"), jetptInV0);

        registryData.fill(HIST("protonQA/V0protonMassInRest_frame"), protonMassInV0);
        registryData.fill(HIST("protonQA/V0protonpxInRest_frame"), protonInV0(1, 0));
        registryData.fill(HIST("protonQA/V0protonpyInRest_frame"), protonInV0(2, 0));
        registryData.fill(HIST("protonQA/V0protonpzInRest_frame"), protonInV0(3, 0));
        double protonsinPhiInV0frame = protonInV0(2, 0) / std::sqrt(protonInV0(1, 0) * protonInV0(1, 0) + protonInV0(2, 0) * protonInV0(2, 0));
        double protoncosthetaInV0frame = protonInV0(3, 0) / std::sqrt(protonInV0(1, 0) * protonInV0(1, 0) + protonInV0(2, 0) * protonInV0(2, 0) + protonInV0(3, 0) * protonInV0(3, 0));
        registryData.fill(HIST("protonQA/V0protonphiInRest_frame"), TMath::ASin(protonsinPhiInV0frame));
        registryData.fill(HIST("protonQA/V0protonthetaInRest_frame"), TMath::ACos(protoncosthetaInV0frame));
        registryData.fill(HIST("protonQA/V0protoncosthetaInV0frame"), protoncosthetaInV0frame);
        registryData.fill(HIST("protonQA/profileprotonsinthetaInV0frame"), candidate.mLambda(), protonsinThetaInV0);
        registryData.fill(HIST("protonQA/profileprotonsinphiInV0frame"), candidate.mLambda(), protonsinPhiInV0frame);
        registryData.fill(HIST("protonQA/profileprotoncosSquarethetaInV0frame"), candidate.mLambda(), protoncosthetaInV0frame * protoncosthetaInV0frame);

        TMatrixD protonInJet(4, 1);
        protonInJet = MyTMatrixTranslationToJet(maxJetpx, maxJetpy, maxJetpz, candidate.px(), candidate.py(), candidate.pz()) * pLabproton;
        double protoncosthetaInJet = protonInJet(3, 0) / std::sqrt(protonInJet(1, 0) * protonInJet(1, 0) + protonInJet(2, 0) * protonInJet(2, 0) + protonInJet(3, 0) * protonInJet(3, 0));
        double protonsinPhiInJet = protonInJet(2, 0) / std::sqrt(protonInJet(1, 0) * protonInJet(1, 0) + protonInJet(2, 0) * protonInJet(2, 0));
        double protonPtinJet = std::sqrt(protonInJet(1, 0) * protonInJet(1, 0) + protonInJet(2, 0) * protonInJet(2, 0));
        double protonPinJet = std::sqrt(protonInJet(1, 0) * protonInJet(1, 0) + protonInJet(2, 0) * protonInJet(2, 0) + protonInJet(3, 0) * protonInJet(3, 0));
        double protonSinThetainJet = protonPtinJet / protonPinJet;
        double protonMassInJetframe = std::sqrt(protonInJet(0, 0) * protonInJet(0, 0) - protonInJet(1, 0) * protonInJet(1, 0) - protonInJet(2, 0) * protonInJet(2, 0) - protonInJet(3, 0) * protonInJet(3, 0));

        TMatrixD pInJet(4, 1);
        pInJet = MyTMatrixTranslationToJet(maxJetpx, maxJetpy, maxJetpz, candidate.px(), candidate.py(), candidate.pz()) * pLabJet;
        double jetthetaInJet = pInJet(3, 0) / std::sqrt(pInJet(1, 0) * pInJet(1, 0) + pInJet(2, 0) * pInJet(2, 0) + pInJet(3, 0) * pInJet(3, 0));
        double jetphiInJet = pInJet(2, 0) / std::sqrt(pInJet(1, 0) * pInJet(1, 0) + pInJet(2, 0) * pInJet(2, 0));
        double jetptInJet = std::sqrt(pInJet(1, 0) * pInJet(1, 0) + pInJet(2, 0) * pInJet(2, 0));
        registryData.fill(HIST("JetQA/JetthetaInJetframe"), TMath::ASin(jetthetaInJet));
        registryData.fill(HIST("JetQA/JetphiInJetframe"), TMath::ASin(jetphiInJet));
        registryData.fill(HIST("JetQA/JetpxInJetframe"), pInJet(1, 0));
        registryData.fill(HIST("JetQA/JetpyInJetframe"), pInJet(2, 0));
        registryData.fill(HIST("JetQA/JetpzInJetframe"), pInJet(3, 0));
        registryData.fill(HIST("JetQA/JetptInJetframe"), jetptInJet);

        registryData.fill(HIST("protonQA/V0protonpxInJetframe"), protonInJet(1, 0));
        registryData.fill(HIST("protonQA/V0protonpyInJetframe"), protonInJet(2, 0));
        registryData.fill(HIST("protonQA/V0protonpzInJetframe"), protonInJet(3, 0));
        registryData.fill(HIST("protonQA/V0protonphiInJetframe"), TMath::ASin(protonsinPhiInJet));
        registryData.fill(HIST("protonQA/V0protonthetaInJetframe"), TMath::ACos(protoncosthetaInJet));
        registryData.fill(HIST("protonQA/V0protoncosthetaInJetframe"), protoncosthetaInJet);
        registryData.fill(HIST("protonQA/profileprotonsinthetaInJetframe"), candidate.mLambda(), protonSinThetainJet);
        registryData.fill(HIST("protonQA/profileprotonsinphiInJetframe"), candidate.mLambda(), protonsinPhiInJet);
        registryData.fill(HIST("protonQA/profileprotoncosSquarethetaInJetframe"), candidate.mLambda(), protoncosthetaInJet * protoncosthetaInJet);
        registryData.fill(HIST("protonQA/V0protonMassInJetframe"), protonMassInJetframe);

        TMatrixD protonInJetV0(4, 1);
        protonInJetV0 = LorentzTransInV0frame(ELambda, lambdaInJet(1, 0), lambdaInJet(2, 0), lambdaInJet(3, 0)) * MyTMatrixTranslationToJet(maxJetpx, maxJetpy, maxJetpz, candidate.px(), candidate.py(), candidate.pz()) * pLabproton;
        double protoncosthetaInJetV0 = protonInJetV0(3, 0) / std::sqrt(protonInJetV0(1, 0) * protonInJetV0(1, 0) + protonInJetV0(2, 0) * protonInJetV0(2, 0) + protonInJetV0(3, 0) * protonInJetV0(3, 0));
        double protonsinphiInJetV0 = protonInJetV0(2, 0) / std::sqrt(protonInJetV0(1, 0) * protonInJetV0(1, 0) + protonInJetV0(2, 0) * protonInJetV0(2, 0));
        double protonPtinJetV0 = std::sqrt(protonInJetV0(1, 0) * protonInJetV0(1, 0) + protonInJetV0(2, 0) * protonInJetV0(2, 0));
        double protonPinJetV0 = std::sqrt(protonInJetV0(1, 0) * protonInJetV0(1, 0) + protonInJetV0(2, 0) * protonInJetV0(2, 0) + protonInJetV0(3, 0) * protonInJetV0(3, 0));
        double protonSinThetainJetV0 = protonPtinJetV0 / protonPinJetV0;
        double protonMassInJetV0frame = std::sqrt(protonInJetV0(0, 0) * protonInJetV0(0, 0) - protonInJetV0(1, 0) * protonInJetV0(1, 0) - protonInJetV0(2, 0) * protonInJetV0(2, 0) - protonInJetV0(3, 0) * protonInJetV0(3, 0));

        TMatrixD JetInJetV0(4, 1);
        JetInJetV0 = LorentzTransInV0frame(ELambda, lambdaInJet(1, 0), lambdaInJet(2, 0), lambdaInJet(3, 0)) * MyTMatrixTranslationToJet(maxJetpx, maxJetpy, maxJetpz, candidate.px(), candidate.py(), candidate.pz()) * pLabJet;
        double jetthetaInJetV0 = JetInJetV0(3, 0) / std::sqrt(JetInJetV0(1, 0) * JetInJetV0(1, 0) + JetInJetV0(2, 0) * JetInJetV0(2, 0) + JetInJetV0(3, 0) * JetInJetV0(3, 0));
        double jetphiInJetV0 = JetInJetV0(2, 0) / std::sqrt(JetInJetV0(1, 0) * JetInJetV0(1, 0) + JetInJetV0(2, 0) * JetInJetV0(2, 0));
        double jetptInJetV0 = std::sqrt(JetInJetV0(1, 0) * JetInJetV0(1, 0) + JetInJetV0(2, 0) * JetInJetV0(2, 0));
        registryData.fill(HIST("JetQA/JetthetaInJetV0frame"), TMath::ASin(jetthetaInJetV0));
        registryData.fill(HIST("JetQA/JetphiInJetV0frame"), TMath::ASin(jetphiInJetV0));
        registryData.fill(HIST("JetQA/JetpxInJetV0frame"), JetInJetV0(1, 0));
        registryData.fill(HIST("JetQA/JetpyInJetV0frame"), JetInJetV0(2, 0));
        registryData.fill(HIST("JetQA/JetpzInJetV0frame"), JetInJetV0(3, 0));
        registryData.fill(HIST("JetQA/JetptInJetV0frame"), jetptInJetV0);

        registryData.fill(HIST("protonQA/V0protonpxInJetV0frame"), protonInJetV0(1, 0));
        registryData.fill(HIST("protonQA/V0protonpyInJetV0frame"), protonInJetV0(2, 0));
        registryData.fill(HIST("protonQA/V0protonpzInJetV0frame"), protonInJetV0(3, 0));
        registryData.fill(HIST("protonQA/V0protonphiInJetV0frame"), TMath::ASin(protonsinphiInJetV0));
        registryData.fill(HIST("protonQA/V0protonthetaInJetV0frame"), TMath::ACos(protoncosthetaInJetV0));
        registryData.fill(HIST("protonQA/V0protoncosthetaInJetV0"), protoncosthetaInJetV0);
        registryData.fill(HIST("protonQA/V0protonMassInJetV0frame"), protonMassInJetV0frame);
        registryData.fill(HIST("protonQA/profileprotonsinthetaInJetV0frame"), candidate.mLambda(), protonSinThetainJetV0);
        registryData.fill(HIST("protonQA/profileprotonsinphiInJetV0frame"), candidate.mLambda(), protonsinphiInJetV0);
        registryData.fill(HIST("protonQA/profileprotoncosSquarethetaInJetV0frame"), candidate.mLambda(), protoncosthetaInJetV0 * protoncosthetaInJetV0);

        double protonCosThetainJetV0 = protonInJetV0(3, 0) / protonPinJetV0;

        protonsinPhiInJetV0frame = protonsinPhiInJetV0frame + protonInJetV0(2, 0) / std::sqrt(protonInJetV0(1, 0) * protonInJetV0(1, 0) + protonInJetV0(2, 0) * protonInJetV0(2, 0));

        registryData.fill(HIST("hprotonsinphiInJetV0frame"), protonsinPhiInJetV0frame);

        registryData.fill(HIST("TProfile2DLambdaPtMassSinPhi"), candidate.mLambda(), candidate.pt(), protonInJetV0(2, 0) / std::sqrt(protonInJetV0(1, 0) * protonInJetV0(1, 0) + protonInJetV0(2, 0) * protonInJetV0(2, 0)));
        registryData.fill(HIST("TProfile2DLambdaPtMassSintheta"), candidate.mLambda(), candidate.pt(), (4.0 / TMath::Pi()) * protonSinThetainJetV0);
        registryData.fill(HIST("TProfile2DLambdaPtMassCosSquareTheta"), candidate.mLambda(), candidate.pt(), 3.0 * protonCosThetainJetV0 * protonCosThetainJetV0);
        registryData.fill(HIST("TProfile2DLambdaMassDeltaPhi"), TMath::ASin(protonsinPhiInJetV0frame), candidate.mLambda(), protonsinPhiInJetV0frame);
        registryData.fill(HIST("hprotonPhi"), TMath::ASin(protonsinPhiInJetV0frame));

        double protonCosThetaInLab = pLabproton(3, 0) / std::sqrt(pLabproton(1, 0) * pLabproton(1, 0) + pLabproton(2, 0) * pLabproton(2, 0) + pLabproton(3, 0) * pLabproton(3, 0));     // cos(theta) of lambda in lab frame
        double protonCosThetaInV0frame = protonInV0(3, 0) / std::sqrt(protonInV0(1, 0) * protonInV0(1, 0) + protonInV0(2, 0) * protonInV0(2, 0) + protonInV0(3, 0) * protonInV0(3, 0)); // cos(theta) of lambda in V0 frame
        double protonCosThetaInJetV0frame = protonCosThetainJetV0;                                                                                                                      // cos(theta) of lambda in jet V0 frame
        registryData.fill(HIST("hprotonThetaInLab"), TMath::ACos(protonCosThetaInLab));
        registryData.fill(HIST("hprotonThetaInV0"), TMath::ACos(protonCosThetaInV0frame));
        registryData.fill(HIST("hprotonThetaInJetV0"), TMath::ACos(protonCosThetaInJetV0frame));
      }
      if (registryDataAcceptV0AntiLambda(candidate, pos, neg, collision)) {
        registryData.fill(HIST("hMassAntiLambda"), candidate.mAntiLambda());
        double PAntiLambda = std::sqrt(candidate.px() * candidate.px() + candidate.py() * candidate.py() + candidate.pz() * candidate.pz());
        double EAntiLambda = std::sqrt(candidate.mAntiLambda() * candidate.mAntiLambda() + PAntiLambda * PAntiLambda);
        double AntiprotonE = std::sqrt(massPr * massPr + neg.px() * neg.px() + neg.py() * neg.py() + neg.pz() * neg.pz());
        TMatrixD pLabAntiV0(4, 1);
        pLabAntiV0(0, 0) = EAntiLambda;
        pLabAntiV0(1, 0) = candidate.px();
        pLabAntiV0(2, 0) = candidate.py();
        pLabAntiV0(3, 0) = candidate.pz();

        TMatrixD AntilambdaInJet(4, 1);
        AntilambdaInJet = MyTMatrixTranslationToJet(maxJetpx, maxJetpy, maxJetpz, candidate.px(), candidate.py(), candidate.pz()) * pLabAntiV0;

        TMatrixD pLabAntiproton(4, 1);
        pLabAntiproton(0, 0) = AntiprotonE;
        pLabAntiproton(1, 0) = neg.px();
        pLabAntiproton(2, 0) = neg.py();
        pLabAntiproton(3, 0) = neg.pz();
        TMatrixD AntiprotonInJetV0(4, 1);
        AntiprotonInJetV0 = LorentzTransInV0frame(EAntiLambda, AntilambdaInJet(1, 0), AntilambdaInJet(2, 0), AntilambdaInJet(3, 0)) * MyTMatrixTranslationToJet(maxJetpx, maxJetpy, maxJetpz, candidate.px(), candidate.py(), candidate.pz()) * pLabAntiproton;
        AntiprotonsinPhiInJetV0frame = AntiprotonsinPhiInJetV0frame + AntiprotonInJetV0(2, 0) / std::sqrt(AntiprotonInJetV0(1, 0) * AntiprotonInJetV0(1, 0) + AntiprotonInJetV0(2, 0) * AntiprotonInJetV0(2, 0));
        TMatrixD AntiprotonInV0(4, 1);
        AntiprotonInV0 = LorentzTransInV0frame(EAntiLambda, candidate.px(), candidate.py(), candidate.pz()) * pLabAntiproton;
        double AntiprotonPinJetV0 = std::sqrt(AntiprotonInJetV0(1, 0) * AntiprotonInJetV0(1, 0) + AntiprotonInJetV0(2, 0) * AntiprotonInJetV0(2, 0) + AntiprotonInJetV0(3, 0) * AntiprotonInJetV0(3, 0));
        double AntiprotonPtinJetV0 = std::sqrt(AntiprotonInJetV0(1, 0) * AntiprotonInJetV0(1, 0) + AntiprotonInJetV0(2, 0) * AntiprotonInJetV0(2, 0));
        double AntiprotonCosThetainJetV0 = AntiprotonInJetV0(3, 0) / AntiprotonPinJetV0;
        double AntiprotonSinThetainJetV0 = AntiprotonPtinJetV0 / AntiprotonPinJetV0;
        registryData.fill(HIST("TProfile2DAntiLambdaPtMassSinPhi"), candidate.mAntiLambda(), candidate.pt(), AntiprotonInJetV0(2, 0) / std::sqrt(AntiprotonInJetV0(1, 0) * AntiprotonInJetV0(1, 0) + AntiprotonInJetV0(2, 0) * AntiprotonInJetV0(2, 0)));
        registryData.fill(HIST("TProfile2DAntiLambdaPtMassSintheta"), candidate.mAntiLambda(), candidate.pt(), (4.0 / TMath::Pi()) * AntiprotonSinThetainJetV0);
        registryData.fill(HIST("TProfile2DAntiLambdaPtMassCosSquareTheta"), candidate.mAntiLambda(), candidate.pt(), 3.0 * AntiprotonCosThetainJetV0 * AntiprotonCosThetainJetV0);
        registryData.fill(HIST("TProfile2DAntiLambdaMassDeltaPhi"), TMath::ASin(AntiprotonsinPhiInJetV0frame), candidate.mAntiLambda(), AntiprotonsinPhiInJetV0frame);
        registryData.fill(HIST("hantiprotonPhi"), TMath::ASin(AntiprotonsinPhiInJetV0frame));
      }
    }

    for (const auto& candidate : fullV0s) {
      const auto& pos = candidate.posTrack_as<StrHadronDaughterTracks>();
      const auto& neg = candidate.negTrack_as<StrHadronDaughterTracks>();
      if (passedLambdaSelection(candidate, pos, neg)) {
        registryData.fill(HIST("hLambdamassandSinPhi"), candidate.mLambda(), protonsinPhiInJetV0frame / V0Numbers);
        registryData.fill(HIST("hLambdaPhiandSinPhi"), TMath::ASin(protonsinPhiInJetV0frame / V0Numbers), protonsinPhiInJetV0frame / V0Numbers);
        registryData.fill(HIST("V0LambdaprotonPhi"), TMath::ASin(protonsinPhiInJetV0frame / V0Numbers));
        registryData.fill(HIST("profileLambda"), candidate.mLambda(), protonsinPhiInJetV0frame / V0Numbers);
      }
      if (passedAntiLambdaSelection(candidate, pos, neg)) {
        registryData.fill(HIST("hAntiLambdamassandSinPhi"), candidate.mAntiLambda(), AntiprotonsinPhiInJetV0frame / AntiV0Numbers);
        registryData.fill(HIST("profileAntiLambda"), candidate.mAntiLambda(), AntiprotonsinPhiInJetV0frame / AntiV0Numbers);
      }
    }
  }
  PROCESS_SWITCH(LambdaJetpolarization, processData, "processData", false);

  // V0Collisions
  // SelCollisions
  using V0Collisions = soa::Join<aod::Collisions, aod::EvSels, aod::PVMults, aod::CentFT0Ms, aod::CentNGlobals>;
  void processLongitudinalPolarization(V0Collisions::iterator const& collision, aod::V0Datas const& fullV0s, StrHadronDaughterTracks const&)
  {

    if (!AcceptEventForLongitudinalPolarization(collision)) {
      return;
    }

    for (const auto& v0 : fullV0s) { // loop over V0s

      if (v0.v0Type() != v0TypeSelection) {
        continue;
      }

      const auto& pos = v0.posTrack_as<StrHadronDaughterTracks>();
      const auto& neg = v0.negTrack_as<StrHadronDaughterTracks>();

      if (notITSAfterburner && (v0.negTrack_as<StrHadronDaughterTracks>().isITSAfterburner() || v0.posTrack_as<StrHadronDaughterTracks>().isITSAfterburner())) {
        continue;
      }

      if (AcceptV0Lambda(v0, pos, neg, collision) && ifpasslambda) {

        ProtonVec = ROOT::Math::PxPyPzMVector(v0.pxpos(), v0.pypos(), v0.pzpos(), massPr);
        PionVec = ROOT::Math::PxPyPzMVector(v0.pxneg(), v0.pyneg(), v0.pzneg(), massPi);
        LambdaVec = ProtonVec + PionVec;
        LambdaVec.SetM(massLambda);
        ROOT::Math::Boost boost{LambdaVec.BoostToCM()};
        ProtonBoostedVec = boost(ProtonVec);
        LambdaBoostedVec = boost(LambdaVec);
      }
      if (AcceptV0AntiLambda(v0, pos, neg, collision) && ifpasslambda) {

        ProtonVec = ROOT::Math::PxPyPzMVector(v0.pxneg(), v0.pyneg(), v0.pzneg(), massPr);
        PionVec = ROOT::Math::PxPyPzMVector(v0.pxpos(), v0.pypos(), v0.pzpos(), massPi);
        LambdaVec = ProtonVec + PionVec;
        LambdaVec.SetM(massLambda);
        ROOT::Math::Boost boost{LambdaVec.BoostToCM()};
        ProtonBoostedVec = boost(ProtonVec);
        LambdaBoostedVec = boost(LambdaVec);
      }
    }
  }
  PROCESS_SWITCH(LambdaJetpolarization, processLongitudinalPolarization, "processLongitudinalPolarization", false);

  void processLambdaJetPolarization(SelV0Collisions::iterator const& collision, aod::V0Datas const& fullV0s, StrHadronDaughterTracks const& tracks)
  {
    registryData.fill(HIST("hNEvents"), 0.5);
    if (!AcceptEvent(collision)) {
      return;
    }
    registryData.fill(HIST("hNEvents"), 8.5);
    // event selection
    // loop over reconstructed tracks
    std::vector<fastjet::PseudoJet> fjParticles;
    for (auto const& track : tracks) {
      registryData.fill(HIST("h_track_pt"), track.pt());
      registryData.fill(HIST("h_track_eta"), track.eta());
      registryData.fill(HIST("h_track_phi"), track.phi());
      if (ispassdTrackSelectionForJetReconstruction && !passedTrackSelectionForJetReconstruction(track)) {
        continue;
      }
      registryData.fill(HIST("h_track_pt_sel"), track.pt());
      registryData.fill(HIST("h_track_eta_sel"), track.eta());
      registryData.fill(HIST("h_track_phi_sel"), track.phi());

      // 4-momentum representation of a particle
      fastjet::PseudoJet fourMomentum(track.px(), track.py(), track.pz(), track.energy(o2::constants::physics::MassPionCharged));
      fjParticles.emplace_back(fourMomentum);
    }
    // reject empty events
    if (fjParticles.size() < 1)
      return;
    // cluster particles using the anti-kt algorithm
    fastjet::RecombinationScheme recombScheme = fastjet::E_scheme;
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, rJet, recombScheme);
    fastjet::AreaDefinition areaDef(fastjet::active_area, fastjet::GhostedAreaSpec(1.0));
    fastjet::ClusterSequenceArea cs(fjParticles, jetDef, areaDef);
    std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());
    // jet selection
    bool isAtLeastOneJetSelected = false;
    int nJets = 0;
    int nJetssel = 0;
    // select most large momentum jet
    float maxJetpx = 0;
    float maxJetpy = 0;
    float maxJetpz = 0;
    float maxJeteta = 0;
    float maxJetphi = 0;
    float maxJetE = 0;
    float maxJetpT = 0;
    float maxJetPt = -999;
    for (const auto& jet : jets) {
      nJets++;
      registryData.fill(HIST("FJetaHistogram"), jet.eta());
      registryData.fill(HIST("FJphiHistogram"), jet.phi());
      registryData.fill(HIST("FJptHistogram"), jet.pt());
      // jet must be fully contained in the acceptance
      if ((std::fabs(jet.eta()) + rJet) > (etaMax - deltaEtaEdge)) {
        continue;
      }

      if (jet.pt() < cfgjetPtMin)
        continue;
      nJetssel++;
      registryData.fill(HIST("FJetaHistogramsel"), jet.eta());
      registryData.fill(HIST("FJphiHistogramsel"), jet.phi());
      registryData.fill(HIST("FJptHistogramsel"), jet.pt());

      if (jet.pt() > maxJetPt) {
        maxJetpx = jet.px();
        maxJetpy = jet.py();
        maxJetpz = jet.pz();
        maxJeteta = jet.eta();
        maxJetE = jet.E();
        maxJetphi = jet.phi();
        maxJetpT = jet.pt();
        maxJetPt = maxJetpT;
      }
    }
    if (maxJetpT > 0) {
      registryData.fill(HIST("FLeadingJetaHistogramsel"), maxJeteta);
      registryData.fill(HIST("FLeadingJphiHistogramsel"), maxJetphi);
      registryData.fill(HIST("FLeadingJptHistogramsel"), maxJetpT);
    }
    registryData.fill(HIST("nJetsPerEvent"), nJets);
    registryData.fill(HIST("nJetsPerEventsel"), nJetssel);
    isAtLeastOneJetSelected = true;
    if (!isAtLeastOneJetSelected) {
      return;
    }

    // Event multiplicity
    const float multiplicity = collision.centFT0M();
    registryData.fill(HIST("number_of_events_vsmultiplicity"), multiplicity);
    // v0 loop
    int V0Numbers = 0;
    int AntiV0Numbers = 0;
    for (const auto& v0 : fullV0s) {
      const auto& pos = v0.posTrack_as<StrHadronDaughterTracks>();
      const auto& neg = v0.negTrack_as<StrHadronDaughterTracks>();
      TVector3 v0dir(v0.px(), v0.py(), v0.pz());
      if (registryDataAcceptV0Lambda(v0, pos, neg, collision)) {
        V0Numbers = V0Numbers + 1;
        registryData.fill(HIST("LambdaPtMass"), v0.pt(), v0.mLambda());
        registryData.fill(HIST("LambdaQA/hArmenterosPreAnalyserCuts"), v0.alpha(), v0.qtarm());
      }
      if (registryDataAcceptV0AntiLambda(v0, pos, neg, collision)) {
        AntiV0Numbers = AntiV0Numbers + 1;
        registryData.fill(HIST("AntiLambdaPtMass"), v0.pt(), v0.mAntiLambda());
        registryData.fill(HIST("AntiLambdaQA/hArmenterosPreAnalyserCuts"), v0.alpha(), v0.qtarm());
      }
    }
    registryData.fill(HIST("nV0sPerEvent"), V0Numbers);

    // calculate lambda polarization induced by jet

    if (V0Numbers == 0 && AntiV0Numbers == 0) {
      return;
    }
    if (maxJetpx == 0) {
      return;
    }
    double protonsinPhiInJetV0frame = 0;
    double AntiprotonsinPhiInJetV0frame = 0;
    for (const auto& candidate : fullV0s) {
      const auto& pos = candidate.posTrack_as<StrHadronDaughterTracks>();
      const auto& neg = candidate.negTrack_as<StrHadronDaughterTracks>();
      TVector3 v0dir(candidate.px(), candidate.py(), candidate.pz());

      if (registryDataAcceptV0Lambda(candidate, pos, neg, collision)) {
        registryData.fill(HIST("hMassLambda"), candidate.mLambda());
        registryData.fill(HIST("V0pTInLab"), candidate.pt());
        registryData.fill(HIST("V0pxInLab"), candidate.px());
        registryData.fill(HIST("V0pyInLab"), candidate.py());
        registryData.fill(HIST("V0pzInLab"), candidate.pz());
        registryData.fill(HIST("protonQA/V0protonpxInLab"), pos.px());
        registryData.fill(HIST("protonQA/V0protonpyInLab"), pos.py());
        registryData.fill(HIST("protonQA/V0protonpzInLab"), pos.pz());

        double PLambda = std::sqrt(candidate.px() * candidate.px() + candidate.py() * candidate.py() + candidate.pz() * candidate.pz());
        double ELambda = std::sqrt(candidate.mLambda() * candidate.mLambda() + PLambda * PLambda);
        double protonE = std::sqrt(massPr * massPr + pos.px() * pos.px() + pos.py() * pos.py() + pos.pz() * pos.pz());

        TMatrixD pLabJet(4, 1);
        pLabJet(0, 0) = maxJetE;
        pLabJet(1, 0) = maxJetpx;
        pLabJet(2, 0) = maxJetpy;
        pLabJet(3, 0) = maxJetpz;

        TMatrixD pLabV0(4, 1);
        pLabV0(0, 0) = ELambda;
        pLabV0(1, 0) = candidate.px();
        pLabV0(2, 0) = candidate.py();
        pLabV0(3, 0) = candidate.pz();

        TMatrixD V0InV0(4, 1);
        V0InV0 = LorentzTransInV0frame(ELambda, candidate.px(), candidate.py(), candidate.pz()) * pLabV0;
        registryData.fill(HIST("V0pxInRest_frame"), V0InV0(1, 0));
        registryData.fill(HIST("V0pyInRest_frame"), V0InV0(2, 0));
        registryData.fill(HIST("V0pzInRest_frame"), V0InV0(3, 0));

        TMatrixD lambdaInJet(4, 1);
        lambdaInJet = TMatrixTranslationToJet(maxJetpx, maxJetpy, maxJetpz, candidate.px(), candidate.py(), candidate.pz()) * pLabV0;

        registryData.fill(HIST("LambdaQA/TH2FLambdaMassPhiInJet"), TMath::ATan2(lambdaInJet(2, 0), lambdaInJet(1, 0)), candidate.mLambda());

        registryData.fill(HIST("V0pxInJetframe"), lambdaInJet(1, 0));
        registryData.fill(HIST("V0pyInJetframe"), lambdaInJet(2, 0));
        registryData.fill(HIST("V0pzInJetframe"), lambdaInJet(3, 0));

        TMatrixD lambdaInJetV0(4, 1);
        lambdaInJetV0 = LorentzTransInV0frame(ELambda, lambdaInJet(1, 0), lambdaInJet(2, 0), lambdaInJet(3, 0)) * MyTMatrixTranslationToJet(maxJetpx, maxJetpy, maxJetpz, candidate.px(), candidate.py(), candidate.pz()) * pLabV0;
        registryData.fill(HIST("V0LambdapxInJetV0frame"), lambdaInJetV0(1, 0));
        registryData.fill(HIST("V0LambdapyInJetV0frame"), lambdaInJetV0(2, 0));
        registryData.fill(HIST("V0LambdapzInJetV0frame"), lambdaInJetV0(3, 0));

        TMatrixD pLabproton(4, 1);
        pLabproton(0, 0) = protonE;
        pLabproton(1, 0) = pos.px();
        pLabproton(2, 0) = pos.py();
        pLabproton(3, 0) = pos.pz();
        double protonsinPhiInLab = pLabproton(2, 0) / std::sqrt(pLabproton(1, 0) * pLabproton(1, 0) + pLabproton(2, 0) * pLabproton(2, 0));
        double protoncosthetaInLab = pLabproton(3, 0) / std::sqrt(pLabproton(1, 0) * pLabproton(1, 0) + pLabproton(2, 0) * pLabproton(2, 0) + pLabproton(3, 0) * pLabproton(3, 0));
        double protonPtInLab = std::sqrt(pLabproton(1, 0) * pLabproton(1, 0) + pLabproton(2, 0) * pLabproton(2, 0));
        double protonPInLab = std::sqrt(pLabproton(1, 0) * pLabproton(1, 0) + pLabproton(2, 0) * pLabproton(2, 0) + pLabproton(3, 0) * pLabproton(3, 0));
        double protonsinThetaInLab = protonPtInLab / protonPInLab;
        double protonMassInLab = std::sqrt(pLabproton(0, 0) * pLabproton(0, 0) - pLabproton(1, 0) * pLabproton(1, 0) - pLabproton(2, 0) * pLabproton(2, 0) - pLabproton(3, 0) * pLabproton(3, 0));
        double jettheta = maxJetpz / std::sqrt(pLabJet(1, 0) * pLabJet(1, 0) + pLabJet(2, 0) * pLabJet(2, 0) + pLabJet(3, 0) * pLabJet(3, 0));
        double jetphi = maxJetpy / std::sqrt(pLabJet(1, 0) * pLabJet(1, 0) + pLabJet(2, 0) * pLabJet(2, 0));
        double jetptInLab = std::sqrt(pLabJet(1, 0) * pLabJet(1, 0) + pLabJet(2, 0) * pLabJet(2, 0));
        registryData.fill(HIST("JetQA/JetthetaInLab"), TMath::ASin(jettheta));
        registryData.fill(HIST("JetQA/JetphiInLab"), TMath::ASin(jetphi));
        registryData.fill(HIST("JetQA/JetpxInLab"), pLabJet(1, 0));
        registryData.fill(HIST("JetQA/JetpyInLab"), pLabJet(2, 0));
        registryData.fill(HIST("JetQA/JetpzInLab"), pLabJet(3, 0));
        registryData.fill(HIST("JetQA/JetptInLab"), jetptInLab);

        registryData.fill(HIST("protonQA/V0protonphiInLab"), TMath::ASin(protonsinPhiInLab));
        registryData.fill(HIST("protonQA/V0protonthetaInLab"), TMath::ACos(protoncosthetaInLab));
        registryData.fill(HIST("protonQA/V0protoncosthetaInLab"), protoncosthetaInLab);
        registryData.fill(HIST("protonQA/profileprotonsinthetaInLab"), candidate.mLambda(), protonsinThetaInLab);
        registryData.fill(HIST("protonQA/profileprotonsinphiInLab"), candidate.mLambda(), protonsinPhiInLab);
        registryData.fill(HIST("protonQA/profileprotoncosSquarethetaInLab"), candidate.mLambda(), protoncosthetaInLab * protoncosthetaInLab);
        registryData.fill(HIST("protonQA/V0protonMassInLab"), protonMassInLab);

        TMatrixD protonInV0(4, 1);
        protonInV0 = LorentzTransInV0frame(ELambda, candidate.px(), candidate.py(), candidate.pz()) * pLabproton;
        double protonMassInV0 = std::sqrt(protonInV0(0, 0) * protonInV0(0, 0) - protonInV0(1, 0) * protonInV0(1, 0) - protonInV0(2, 0) * protonInV0(2, 0) - protonInV0(3, 0) * protonInV0(3, 0));
        double protonPInV0 = std::sqrt(protonInV0(1, 0) * protonInV0(1, 0) + protonInV0(2, 0) * protonInV0(2, 0) + protonInV0(3, 0) * protonInV0(3, 0));
        double protonPtInV0 = std::sqrt(protonInV0(1, 0) * protonInV0(1, 0) + protonInV0(2, 0) * protonInV0(2, 0));
        double protonsinThetaInV0 = protonPtInV0 / protonPInV0;

        TMatrixD JetInV0(4, 1);
        JetInV0 = LorentzTransInV0frame(ELambda, candidate.px(), candidate.py(), candidate.pz()) * pLabJet;
        double jetthetaInV0 = JetInV0(3, 0) / std::sqrt(JetInV0(1, 0) * JetInV0(1, 0) + JetInV0(2, 0) * JetInV0(2, 0) + JetInV0(3, 0) * JetInV0(3, 0));
        double jetphiInV0 = JetInV0(2, 0) / std::sqrt(JetInV0(1, 0) * JetInV0(1, 0) + JetInV0(2, 0) * JetInV0(2, 0));
        double jetptInV0 = std::sqrt(JetInV0(1, 0) * JetInV0(1, 0) + JetInV0(2, 0) * JetInV0(2, 0));
        registryData.fill(HIST("JetQA/JetthetaInV0"), TMath::ASin(jetthetaInV0));
        registryData.fill(HIST("JetQA/JetphiInV0"), TMath::ASin(jetphiInV0));
        registryData.fill(HIST("JetQA/JetpxInV0"), JetInV0(1, 0));
        registryData.fill(HIST("JetQA/JetpyInV0"), JetInV0(2, 0));
        registryData.fill(HIST("JetQA/JetpzInV0"), JetInV0(3, 0));
        registryData.fill(HIST("JetQA/JetptInV0"), jetptInV0);

        registryData.fill(HIST("protonQA/V0protonMassInRest_frame"), protonMassInV0);
        registryData.fill(HIST("protonQA/V0protonpxInRest_frame"), protonInV0(1, 0));
        registryData.fill(HIST("protonQA/V0protonpyInRest_frame"), protonInV0(2, 0));
        registryData.fill(HIST("protonQA/V0protonpzInRest_frame"), protonInV0(3, 0));
        double protonsinPhiInV0frame = protonInV0(2, 0) / std::sqrt(protonInV0(1, 0) * protonInV0(1, 0) + protonInV0(2, 0) * protonInV0(2, 0));
        double protoncosthetaInV0frame = protonInV0(3, 0) / std::sqrt(protonInV0(1, 0) * protonInV0(1, 0) + protonInV0(2, 0) * protonInV0(2, 0) + protonInV0(3, 0) * protonInV0(3, 0));
        registryData.fill(HIST("protonQA/V0protonphiInRest_frame"), TMath::ASin(protonsinPhiInV0frame));
        registryData.fill(HIST("protonQA/V0protonthetaInRest_frame"), TMath::ACos(protoncosthetaInV0frame));
        registryData.fill(HIST("protonQA/V0protoncosthetaInV0frame"), protoncosthetaInV0frame);
        registryData.fill(HIST("protonQA/profileprotonsinthetaInV0frame"), candidate.mLambda(), protonsinThetaInV0);
        registryData.fill(HIST("protonQA/profileprotonsinphiInV0frame"), candidate.mLambda(), protonsinPhiInV0frame);
        registryData.fill(HIST("protonQA/profileprotoncosSquarethetaInV0frame"), candidate.mLambda(), protoncosthetaInV0frame * protoncosthetaInV0frame);

        TMatrixD protonInJet(4, 1);
        protonInJet = TMatrixTranslationToJet(maxJetpx, maxJetpy, maxJetpz, candidate.px(), candidate.py(), candidate.pz()) * pLabproton;
        double protoncosthetaInJet = protonInJet(3, 0) / std::sqrt(protonInJet(1, 0) * protonInJet(1, 0) + protonInJet(2, 0) * protonInJet(2, 0) + protonInJet(3, 0) * protonInJet(3, 0));
        double protonsinPhiInJet = protonInJet(2, 0) / std::sqrt(protonInJet(1, 0) * protonInJet(1, 0) + protonInJet(2, 0) * protonInJet(2, 0));
        double protonPtinJet = std::sqrt(protonInJet(1, 0) * protonInJet(1, 0) + protonInJet(2, 0) * protonInJet(2, 0));
        double protonPinJet = std::sqrt(protonInJet(1, 0) * protonInJet(1, 0) + protonInJet(2, 0) * protonInJet(2, 0) + protonInJet(3, 0) * protonInJet(3, 0));
        double protonSinThetainJet = protonPtinJet / protonPinJet;
        double protonMassInJetframe = std::sqrt(protonInJet(0, 0) * protonInJet(0, 0) - protonInJet(1, 0) * protonInJet(1, 0) - protonInJet(2, 0) * protonInJet(2, 0) - protonInJet(3, 0) * protonInJet(3, 0));

        TMatrixD pInJet(4, 1);
        pInJet = TMatrixTranslationToJet(maxJetpx, maxJetpy, maxJetpz, candidate.px(), candidate.py(), candidate.pz()) * pLabJet;
        double jetthetaInJet = pInJet(3, 0) / std::sqrt(pInJet(1, 0) * pInJet(1, 0) + pInJet(2, 0) * pInJet(2, 0) + pInJet(3, 0) * pInJet(3, 0));
        double jetphiInJet = pInJet(2, 0) / std::sqrt(pInJet(1, 0) * pInJet(1, 0) + pInJet(2, 0) * pInJet(2, 0));
        double jetptInJet = std::sqrt(pInJet(1, 0) * pInJet(1, 0) + pInJet(2, 0) * pInJet(2, 0));
        registryData.fill(HIST("JetQA/JetthetaInJetframe"), TMath::ASin(jetthetaInJet));
        registryData.fill(HIST("JetQA/JetphiInJetframe"), TMath::ASin(jetphiInJet));
        registryData.fill(HIST("JetQA/JetpxInJetframe"), pInJet(1, 0));
        registryData.fill(HIST("JetQA/JetpyInJetframe"), pInJet(2, 0));
        registryData.fill(HIST("JetQA/JetpzInJetframe"), pInJet(3, 0));
        registryData.fill(HIST("JetQA/JetptInJetframe"), jetptInJet);

        registryData.fill(HIST("protonQA/V0protonpxInJetframe"), protonInJet(1, 0));
        registryData.fill(HIST("protonQA/V0protonpyInJetframe"), protonInJet(2, 0));
        registryData.fill(HIST("protonQA/V0protonpzInJetframe"), protonInJet(3, 0));
        registryData.fill(HIST("protonQA/V0protonphiInJetframe"), TMath::ASin(protonsinPhiInJet));
        registryData.fill(HIST("protonQA/V0protonthetaInJetframe"), TMath::ACos(protoncosthetaInJet));
        registryData.fill(HIST("protonQA/V0protoncosthetaInJetframe"), protoncosthetaInJet);
        registryData.fill(HIST("protonQA/profileprotonsinthetaInJetframe"), candidate.mLambda(), protonSinThetainJet);
        registryData.fill(HIST("protonQA/profileprotonsinphiInJetframe"), candidate.mLambda(), protonsinPhiInJet);
        registryData.fill(HIST("protonQA/profileprotoncosSquarethetaInJetframe"), candidate.mLambda(), protoncosthetaInJet * protoncosthetaInJet);
        registryData.fill(HIST("protonQA/V0protonMassInJetframe"), protonMassInJetframe);

        TMatrixD protonInJetV0(4, 1);
        protonInJetV0 = LorentzTransInV0frame(ELambda, lambdaInJet(1, 0), lambdaInJet(2, 0), lambdaInJet(3, 0)) * TMatrixTranslationToJet(maxJetpx, maxJetpy, maxJetpz, candidate.px(), candidate.py(), candidate.pz()) * pLabproton;
        double protoncosthetaInJetV0 = protonInJetV0(3, 0) / std::sqrt(protonInJetV0(1, 0) * protonInJetV0(1, 0) + protonInJetV0(2, 0) * protonInJetV0(2, 0) + protonInJetV0(3, 0) * protonInJetV0(3, 0));
        double protonsinphiInJetV0 = protonInJetV0(2, 0) / std::sqrt(protonInJetV0(1, 0) * protonInJetV0(1, 0) + protonInJetV0(2, 0) * protonInJetV0(2, 0));
        double protonPtinJetV0 = std::sqrt(protonInJetV0(1, 0) * protonInJetV0(1, 0) + protonInJetV0(2, 0) * protonInJetV0(2, 0));
        double protonPinJetV0 = std::sqrt(protonInJetV0(1, 0) * protonInJetV0(1, 0) + protonInJetV0(2, 0) * protonInJetV0(2, 0) + protonInJetV0(3, 0) * protonInJetV0(3, 0));
        double protonSinThetainJetV0 = protonPtinJetV0 / protonPinJetV0;
        double protonMassInJetV0frame = std::sqrt(protonInJetV0(0, 0) * protonInJetV0(0, 0) - protonInJetV0(1, 0) * protonInJetV0(1, 0) - protonInJetV0(2, 0) * protonInJetV0(2, 0) - protonInJetV0(3, 0) * protonInJetV0(3, 0));

        TMatrixD JetInJetV0(4, 1);
        JetInJetV0 = LorentzTransInV0frame(ELambda, lambdaInJet(1, 0), lambdaInJet(2, 0), lambdaInJet(3, 0)) * TMatrixTranslationToJet(maxJetpx, maxJetpy, maxJetpz, candidate.px(), candidate.py(), candidate.pz()) * pLabJet;
        double jetthetaInJetV0 = JetInJetV0(3, 0) / std::sqrt(JetInJetV0(1, 0) * JetInJetV0(1, 0) + JetInJetV0(2, 0) * JetInJetV0(2, 0) + JetInJetV0(3, 0) * JetInJetV0(3, 0));
        double jetphiInJetV0 = JetInJetV0(2, 0) / std::sqrt(JetInJetV0(1, 0) * JetInJetV0(1, 0) + JetInJetV0(2, 0) * JetInJetV0(2, 0));
        double jetptInJetV0 = std::sqrt(JetInJetV0(1, 0) * JetInJetV0(1, 0) + JetInJetV0(2, 0) * JetInJetV0(2, 0));
        registryData.fill(HIST("JetQA/JetthetaInJetV0frame"), TMath::ASin(jetthetaInJetV0));
        registryData.fill(HIST("JetQA/JetphiInJetV0frame"), TMath::ASin(jetphiInJetV0));
        registryData.fill(HIST("JetQA/JetpxInJetV0frame"), JetInJetV0(1, 0));
        registryData.fill(HIST("JetQA/JetpyInJetV0frame"), JetInJetV0(2, 0));
        registryData.fill(HIST("JetQA/JetpzInJetV0frame"), JetInJetV0(3, 0));
        registryData.fill(HIST("JetQA/JetptInJetV0frame"), jetptInJetV0);

        registryData.fill(HIST("protonQA/V0protonpxInJetV0frame"), protonInJetV0(1, 0));
        registryData.fill(HIST("protonQA/V0protonpyInJetV0frame"), protonInJetV0(2, 0));
        registryData.fill(HIST("protonQA/V0protonpzInJetV0frame"), protonInJetV0(3, 0));
        registryData.fill(HIST("protonQA/V0protonphiInJetV0frame"), TMath::ASin(protonsinphiInJetV0));
        registryData.fill(HIST("protonQA/V0protonthetaInJetV0frame"), TMath::ACos(protoncosthetaInJetV0));
        registryData.fill(HIST("protonQA/V0protoncosthetaInJetV0"), protoncosthetaInJetV0);
        registryData.fill(HIST("protonQA/V0protonMassInJetV0frame"), protonMassInJetV0frame);
        registryData.fill(HIST("protonQA/profileprotonsinthetaInJetV0frame"), candidate.mLambda(), protonSinThetainJetV0);
        registryData.fill(HIST("protonQA/profileprotonsinphiInJetV0frame"), candidate.mLambda(), protonsinphiInJetV0);
        registryData.fill(HIST("protonQA/profileprotoncosSquarethetaInJetV0frame"), candidate.mLambda(), protoncosthetaInJetV0 * protoncosthetaInJetV0);

        double protonCosThetainJetV0 = protonInJetV0(3, 0) / protonPinJetV0;

        protonsinPhiInJetV0frame = protonsinPhiInJetV0frame + protonInJetV0(2, 0) / std::sqrt(protonInJetV0(1, 0) * protonInJetV0(1, 0) + protonInJetV0(2, 0) * protonInJetV0(2, 0));

        registryData.fill(HIST("hprotonsinphiInJetV0frame"), protonsinPhiInJetV0frame);

        registryData.fill(HIST("TProfile2DLambdaPtMassSinPhi"), candidate.mLambda(), candidate.pt(), protonInJetV0(2, 0) / std::sqrt(protonInJetV0(1, 0) * protonInJetV0(1, 0) + protonInJetV0(2, 0) * protonInJetV0(2, 0)));
        registryData.fill(HIST("TProfile2DLambdaPtMassSintheta"), candidate.mLambda(), candidate.pt(), (4.0 / TMath::Pi()) * protonSinThetainJetV0);
        registryData.fill(HIST("TProfile2DLambdaPtMassCosSquareTheta"), candidate.mLambda(), candidate.pt(), 3.0 * protonCosThetainJetV0 * protonCosThetainJetV0);
        registryData.fill(HIST("TProfile2DLambdaMassDeltaPhi"), TMath::ASin(protonsinPhiInJetV0frame), candidate.mLambda(), protonsinPhiInJetV0frame);
        registryData.fill(HIST("hprotonPhi"), TMath::ASin(protonsinPhiInJetV0frame));

        double protonCosThetaInLab = pLabproton(3, 0) / std::sqrt(pLabproton(1, 0) * pLabproton(1, 0) + pLabproton(2, 0) * pLabproton(2, 0) + pLabproton(3, 0) * pLabproton(3, 0));     // cos(theta) of lambda in lab frame
        double protonCosThetaInV0frame = protonInV0(3, 0) / std::sqrt(protonInV0(1, 0) * protonInV0(1, 0) + protonInV0(2, 0) * protonInV0(2, 0) + protonInV0(3, 0) * protonInV0(3, 0)); // cos(theta) of lambda in V0 frame
        double protonCosThetaInJetV0frame = protonCosThetainJetV0;
        double protonCosThetaInJet = protonInJet(3, 0) / std::sqrt(protonInJet(1, 0) * protonInJet(1, 0) + protonInJet(2, 0) * protonInJet(2, 0) + protonInJet(3, 0) * protonInJet(3, 0)); // cos(theta) of lambda in Jet frame

        registryData.fill(HIST("hprotonThetaInLab"), TMath::ACos(protonCosThetaInLab));
        registryData.fill(HIST("hprotonThetaInV0"), TMath::ACos(protonCosThetaInV0frame));
        registryData.fill(HIST("hprotonThetaInJetV0"), TMath::ACos(protonCosThetaInJetV0frame));

        registryData.fill(HIST("LambdaQA/TH2FprotonCosThetaInLab"), candidate.mLambda(), protonCosThetaInLab);
        registryData.fill(HIST("LambdaQA/TProfile1DprotonCosThetaInLab"), candidate.mLambda(), protonCosThetaInLab);
        registryData.fill(HIST("LambdaQA/TProfile1DprotonCos2ThetaInLab"), candidate.mLambda(), protonCosThetaInLab * protonCosThetaInLab);

        registryData.fill(HIST("LambdaQA/TH2FprotonCosThetaInV0"), candidate.mLambda(), protonCosThetaInV0frame);
        registryData.fill(HIST("LambdaQA/TProfile1DprotonCosThetaInV0"), candidate.mLambda(), protonCosThetaInV0frame);
        registryData.fill(HIST("LambdaQA/TProfile1DprotonCos2ThetaInV0"), candidate.mLambda(), protonCosThetaInV0frame * protonCosThetaInV0frame);

        registryData.fill(HIST("LambdaQA/TH2FprotonCosThetaInJet"), candidate.mLambda(), protonCosThetaInJet);
        registryData.fill(HIST("LambdaQA/TProfile1DprotonCosThetaInJet"), candidate.mLambda(), protonCosThetaInJet);
        registryData.fill(HIST("LambdaQA/TProfile1DprotonCos2ThetaInJet"), candidate.mLambda(), protonCosThetaInJet * protonCosThetaInJet);

        registryData.fill(HIST("LambdaQA/TH2FprotonCosThetaInJetV0"), candidate.mLambda(), protonCosThetaInJetV0frame);
        registryData.fill(HIST("LambdaQA/TProfile1DprotonCosThetaInJetV0"), candidate.mLambda(), protonCosThetaInJetV0frame);
        registryData.fill(HIST("LambdaQA/TProfile1DprotonCos2ThetaInJetV0"), candidate.mLambda(), protonCosThetaInJetV0frame * protonCosThetaInJetV0frame);
        registryData.fill(HIST("LambdaQA/TProfile2DprotonCosThetaInJetV0"), candidate.mLambda(), TMath::ATan2(lambdaInJet(2, 0), lambdaInJet(1, 0)), protonCosThetaInJetV0frame);
        registryData.fill(HIST("LambdaQA/TProfile2DprotonCos2ThetaInJetV0"), candidate.mLambda(), TMath::ATan2(lambdaInJet(2, 0), lambdaInJet(1, 0)), protonCosThetaInJetV0frame * protonCosThetaInJetV0frame);
      }
      if (registryDataAcceptV0AntiLambda(candidate, pos, neg, collision)) {
        registryData.fill(HIST("hMassAntiLambda"), candidate.mAntiLambda());
        double PAntiLambda = std::sqrt(candidate.px() * candidate.px() + candidate.py() * candidate.py() + candidate.pz() * candidate.pz());
        double EAntiLambda = std::sqrt(candidate.mAntiLambda() * candidate.mAntiLambda() + PAntiLambda * PAntiLambda);
        double AntiprotonE = std::sqrt(massPr * massPr + neg.px() * neg.px() + neg.py() * neg.py() + neg.pz() * neg.pz());
        TMatrixD pLabAntiV0(4, 1);
        pLabAntiV0(0, 0) = EAntiLambda;
        pLabAntiV0(1, 0) = candidate.px();
        pLabAntiV0(2, 0) = candidate.py();
        pLabAntiV0(3, 0) = candidate.pz();

        TMatrixD AntilambdaInJet(4, 1);
        AntilambdaInJet = MyTMatrixTranslationToJet(maxJetpx, maxJetpy, maxJetpz, candidate.px(), candidate.py(), candidate.pz()) * pLabAntiV0;

        TMatrixD pLabAntiproton(4, 1);
        pLabAntiproton(0, 0) = AntiprotonE;
        pLabAntiproton(1, 0) = neg.px();
        pLabAntiproton(2, 0) = neg.py();
        pLabAntiproton(3, 0) = neg.pz();
        TMatrixD AntiprotonInJetV0(4, 1);
        AntiprotonInJetV0 = LorentzTransInV0frame(EAntiLambda, AntilambdaInJet(1, 0), AntilambdaInJet(2, 0), AntilambdaInJet(3, 0)) * MyTMatrixTranslationToJet(maxJetpx, maxJetpy, maxJetpz, candidate.px(), candidate.py(), candidate.pz()) * pLabAntiproton;
        AntiprotonsinPhiInJetV0frame = AntiprotonsinPhiInJetV0frame + AntiprotonInJetV0(2, 0) / std::sqrt(AntiprotonInJetV0(1, 0) * AntiprotonInJetV0(1, 0) + AntiprotonInJetV0(2, 0) * AntiprotonInJetV0(2, 0));
        TMatrixD AntiprotonInV0(4, 1);
        AntiprotonInV0 = LorentzTransInV0frame(EAntiLambda, candidate.px(), candidate.py(), candidate.pz()) * pLabAntiproton;
        double AntiprotonPinJetV0 = std::sqrt(AntiprotonInJetV0(1, 0) * AntiprotonInJetV0(1, 0) + AntiprotonInJetV0(2, 0) * AntiprotonInJetV0(2, 0) + AntiprotonInJetV0(3, 0) * AntiprotonInJetV0(3, 0));
        double AntiprotonPtinJetV0 = std::sqrt(AntiprotonInJetV0(1, 0) * AntiprotonInJetV0(1, 0) + AntiprotonInJetV0(2, 0) * AntiprotonInJetV0(2, 0));
        double AntiprotonCosThetainJetV0 = AntiprotonInJetV0(3, 0) / AntiprotonPinJetV0;
        double AntiprotonSinThetainJetV0 = AntiprotonPtinJetV0 / AntiprotonPinJetV0;
        registryData.fill(HIST("TProfile2DAntiLambdaPtMassSinPhi"), candidate.mAntiLambda(), candidate.pt(), AntiprotonInJetV0(2, 0) / std::sqrt(AntiprotonInJetV0(1, 0) * AntiprotonInJetV0(1, 0) + AntiprotonInJetV0(2, 0) * AntiprotonInJetV0(2, 0)));
        registryData.fill(HIST("TProfile2DAntiLambdaPtMassSintheta"), candidate.mAntiLambda(), candidate.pt(), (4.0 / TMath::Pi()) * AntiprotonSinThetainJetV0);
        registryData.fill(HIST("TProfile2DAntiLambdaPtMassCosSquareTheta"), candidate.mAntiLambda(), candidate.pt(), 3.0 * AntiprotonCosThetainJetV0 * AntiprotonCosThetainJetV0);
        registryData.fill(HIST("TProfile2DAntiLambdaMassDeltaPhi"), TMath::ASin(AntiprotonsinPhiInJetV0frame), candidate.mAntiLambda(), AntiprotonsinPhiInJetV0frame);
        registryData.fill(HIST("hantiprotonPhi"), TMath::ASin(AntiprotonsinPhiInJetV0frame));
      }
    }

    for (const auto& candidate : fullV0s) {
      const auto& pos = candidate.posTrack_as<StrHadronDaughterTracks>();
      const auto& neg = candidate.negTrack_as<StrHadronDaughterTracks>();
      if (passedLambdaSelection(candidate, pos, neg)) {
        registryData.fill(HIST("hLambdamassandSinPhi"), candidate.mLambda(), protonsinPhiInJetV0frame / V0Numbers);
        registryData.fill(HIST("hLambdaPhiandSinPhi"), TMath::ASin(protonsinPhiInJetV0frame / V0Numbers), protonsinPhiInJetV0frame / V0Numbers);
        registryData.fill(HIST("V0LambdaprotonPhi"), TMath::ASin(protonsinPhiInJetV0frame / V0Numbers));
        registryData.fill(HIST("profileLambda"), candidate.mLambda(), protonsinPhiInJetV0frame / V0Numbers);
      }
      if (passedAntiLambdaSelection(candidate, pos, neg)) {
        registryData.fill(HIST("hAntiLambdamassandSinPhi"), candidate.mAntiLambda(), AntiprotonsinPhiInJetV0frame / AntiV0Numbers);
        registryData.fill(HIST("profileAntiLambda"), candidate.mAntiLambda(), AntiprotonsinPhiInJetV0frame / AntiV0Numbers);
      }
    }
  }
  PROCESS_SWITCH(LambdaJetpolarization, processLambdaJetPolarization, "processLambdaJetPolarization", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<LambdaJetpolarization>(cfgc)}; // TaskName{"lambdaJetpolarization"}
}
