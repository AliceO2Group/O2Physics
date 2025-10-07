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
#include <TLorentzVector.h>
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
#include <iostream>
#include <string>
#include <vector>

using std::cout;
using std::endl;
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct LfMyV0s {
  HistogramRegistry registry{"registry"};
  HistogramRegistry registryData{"registryData", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryV0Data{"registryV0Data", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryLongitudinalPolarization{"registryLongitudinalPolarization", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  Configurable<double> zVtx{"zVtx", 10.0, "Maximum zVertex"};
  Configurable<double> rJet{"rJet", 0.4, "Jet resolution parameter R"};
  Configurable<float> etaMin{"etaMin", -0.9f, "eta min"};
  Configurable<float> etaMax{"etaMax", +0.9f, "eta max"};
  Configurable<double> deltaEtaEdge{"deltaEtaEdge", 0.00, "eta gap from the edge"};
  // track parameters
  Configurable<float> minITSnCls{"minITSnCls", 4.0f, "min number of ITS clusters"};
  Configurable<float> minTPCnClsFound{"minTPCnClsFound", 80.0f, "min number of found TPC clusters"};
  Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 80.0f, "min number of TPC crossed rows"};
  Configurable<float> minTpcNcrossedRowsOverFindable{"minTpcNcrossedRowsOverFindable", 0.8, "crossed rows/findable"};
  Configurable<float> maxChi2TPC{"maxChi2TPC", 4.0f, "max chi2 per cluster TPC"};
  Configurable<float> maxChi2ITS{"maxChi2ITS", 36.0f, "max chi2 per cluster ITS"};
  Configurable<bool> requireTOF{"requireTOF", false, "require TOF hit"};
  Configurable<bool> requireITS{"requireITS", false, "require ITS hit"};
  Configurable<float> max_tpcSharedCls{"max_tpcSharedCls", 100, "max_tpcSharedCls"};
  Configurable<float> max_chi2_TPC{"max_chi2_TPC", 4, "max_chi2_TPC"};
  Configurable<float> max_chi2_ITS{"max_chi2_ITS", 36, "max_chi2_ITS"};

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
  Configurable<float> v0cospaMin{"v0cospaMin", 0.995f, "Minimum V0 CosPA"};
  Configurable<float> v0cospainit{"v0cospainit", 0.97f, "Minimum V0 CosPA"};
  Configurable<float> minimumV0Radius{"minimumV0Radius", 0.2f, "Minimum V0 Radius"};
  Configurable<float> maximumV0Radius{"maximumV0Radius", 100000.0f, "Maximum V0 Radius"};
  Configurable<float> dcaV0DaughtersMax{"dcaV0DaughtersMax", 1.0f, "Maximum DCA Daughters"};
  Configurable<float> dcanegtoPVmin{"dcanegtoPVmin", 0.1f, "Minimum DCA Neg To PV"};
  Configurable<float> dcapostoPVmin{"dcapostoPVmin", 0.1f, "Minimum DCA Pos To PV"};
  Configurable<float> v0radius{"v0radius", 0.0, "Radius"};
  Configurable<float> dcav0dau{"dcav0dau", 10, "DCA V0 Daughters"};
  Configurable<float> dcanegtopv{"dcanegtopv", 0.0, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", 0.0, "DCA Pos To PV"};

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
  Configurable<float> CtauLambda{"ctauLambda", 30, "C tau Lambda (cm)"};
  Configurable<bool> requirepassedSingleTrackSelection{"requirepassedSingleTrackSelection", false, "requirepassedSingleTrackSelection"};
  Configurable<float> V0tracketaMin{"V0tracketaMin", -0.8f, "eta min track"};
  Configurable<float> V0tracketaMax{"V0tracketaMax", +0.8f, "eta max track"};
  Configurable<bool> requireTPC{"requireTPC", true, "require TPC hit"};
  Configurable<float> yMin{"V0yMin", -0.5f, "minimum y"};
  Configurable<float> yMax{"V0yMax", +0.5f, "maximum y"};
  Configurable<float> v0rejLambda{"v0rejLambda", 0.01, "V0 rej Lambda"};
  Configurable<float> v0accLambda{"v0accLambda", 0.075, "V0 acc Lambda"};
  Configurable<bool> ifinitpasslambda{"ifinitpasslambda", 0, "ifinitpasslambda"};
  Configurable<bool> ifpasslambda{"passedLambdaSelection", 1, "passedLambdaSelection"};
  Configurable<float> paramArmenterosCut{"paramArmenterosCut", 0.2, "parameter Armenteros Cut"};
  Configurable<bool> doArmenterosCut{"doArmenterosCut", 0, "do Armenteros Cut"};
  Configurable<bool> noSameBunchPileUp{"noSameBunchPileUp", true, "reject SameBunchPileUp"};
  Configurable<int> v0TypeSelection{"v0TypeSelection", 1, "select on a certain V0 type (leave negative if no selection desired)"};
  Configurable<bool> NotITSAfterburner{"NotITSAfterburner", 0, "NotITSAfterburner"};
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
    const AxisSpec axisPhi{100, -3.14, 3.14, "#Phi"};
    const AxisSpec axisMass{100, 0, 2, "Mass(GeV/c^{2})"};

    const AxisSpec JetaxisEta{30, -1.5, +1.5, "#eta"};
    const AxisSpec JetaxisPhi{200, -1, +7, "#phi"};
    const AxisSpec JetaxisPt{200, 0, +200, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec ptAxis{100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec invMassLambdaAxis{200, 1.016, 1.216, "m_{p#pi} (GeV/#it{c}^{2})"};

    ConfigurableAxis TProfile2DaxisPt{"#it{p}_{T} (GeV/#it{c})", {VARIABLE_WIDTH, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.2, 3.7, 4.2, 5, 6, 8, 10, 12}, "pt axis for histograms"};
    ConfigurableAxis TProfile2DaxisMass{"Mass p#pi (GeV/#it{c^{2}})", {VARIABLE_WIDTH, 1.10068, 1.10668, 1.11068, 1.11268, 1.11368, 1.11468, 1.11568, 1.11668, 1.11768, 1.11868, 1.12068, 1.12468, 1.13068}, "Mass axis for histograms"};

    registry.add("hMassLambda", "hMassLambda", {HistType::kTH1F, {invMassLambdaAxis}});
    registry.add("V0pTInLab", "V0pTInLab", kTH1F, {axisPT});
    registry.add("hMassVsPtLambda", "hMassVsPtLambda", {HistType::kTH2F, {ptAxis, invMassLambdaAxis}});
    registry.add("hMassVsPtAntiLambda", "hMassVsPtAntiLambda", {HistType::kTH2F, {ptAxis, invMassLambdaAxis}});

    registry.add("V0pxInLab", "V0pxInLab", kTH1F, {axisPx});
    registry.add("V0pyInLab", "V0pyInLab", kTH1F, {axisPy});
    registry.add("V0pzInLab", "V0pzInLab", kTH1F, {axisPz});

    registry.add("V0pxInRest_frame", "V0pxInRest_frame", kTH1F, {axisPx});
    registry.add("V0pyInRest_frame", "V0pyInRest_frame", kTH1F, {axisPy});
    registry.add("V0pzInRest_frame", "V0pzInRest_frame", kTH1F, {axisPz});

    registry.add("JetpxInLab", "JetpxInLab", kTH1F, {axisPx});
    registry.add("JetpyInLab", "JetpyInLab", kTH1F, {axisPy});
    registry.add("JetpzInLab", "JetpzInLab", kTH1F, {axisPz});
    registry.add("JetpTInLab", "JetpTInLab", kTH1F, {axisPT});

    registry.add("LeadingJetpx", "LeadingJetpx", kTH1F, {axisPx});
    registry.add("LeadingJetpy", "LeadingJetpy", kTH1F, {axisPy});
    registry.add("LeadingJetpz", "LeadingJetpz", kTH1F, {axisPz});
    registry.add("LeadingJetpT", "LeadingJetpT", kTH1F, {axisPT});

    registry.add("V0protonpxInLab", "V0protonpxInLab", kTH1F, {axisPx});
    registry.add("V0protonpyInLab", "V0protonpyInLab", kTH1F, {axisPy});
    registry.add("V0protonpzInLab", "V0protonpzInLab", kTH1F, {axisPz});
    registry.add("V0protonphiInLab", "V0protonphiInLab", kTH1F, {axisPhi});

    registry.add("V0protonpxInRest_frame", "V0protonpxInRest_frame", kTH1F, {axisPx});
    registry.add("V0protonpyInRest_frame", "V0protonpyInRest_frame", kTH1F, {axisPy});
    registry.add("V0protonpzInRest_frame", "V0protonpzInRest_frame", kTH1F, {axisPz});
    registry.add("V0protonMassInRest_frame", "V0protonMassInRest_frame", kTH1F, {axisMass});
    registry.add("V0protonphiInRest_frame", "V0protonphiInRest_frame", kTH1F, {axisPhi});

    registry.add("V0protonpxInJetV0frame", "V0protonpxInJetV0frame", kTH1F, {axisPx});
    registry.add("V0protonpyInJetV0frame", "V0protonpyInJetV0frame", kTH1F, {axisPy});
    registry.add("V0protonpzInJetV0frame", "V0protonpzInJetV0frame", kTH1F, {axisPz});
    registry.add("V0protonphiInJetV0frame", "V0protonphiInJetV0frame", kTH1F, {axisPhi});
    registry.add("V0antiprotonphiInJetV0frame", "V0antiprotonphiInJetV0frame", kTH1F, {axisPhi});

    registry.add("V0LambdapxInJetV0frame", "V0LambdapxInJetV0frame", kTH1F, {axisPx});
    registry.add("V0LambdapyInJetV0frame", "V0LambdapyInJetV0frame", kTH1F, {axisPy});
    registry.add("V0LambdapzInJetV0frame", "V0LambdapzInJetV0frame", kTH1F, {axisPz});

    registry.add("hLambdamassandSinPhi", "hLambdamassandSinPhi", kTH2F, {{200, 0.9, 1.2}, {200, -1, 1}});
    registry.add("hAntiLambdamassandSinPhi", "hAntiLambdamassandSinPhi", kTH2F, {{200, 0.9, 1.2}, {200, -1, 1}});
    registry.add("profile", "Invariant Mass vs sin(phi)", {HistType::kTProfile, {{200, 0.9, 1.2}}});
    registry.add("profileAntiV0", "Invariant Mass vs sin(phi)", {HistType::kTProfile, {{200, 0.9, 1.2}}});
    registry.add("hLambdaPhiandSinPhi", "hLambdaPhiandSinPhi", kTH2F, {{200, -TMath::Pi() / 2, TMath::Pi() / 2}, {200, -1, 1}});
    registry.add("hAntiLambdaPhiandSinPhi", "hAntiLambdaPhiandSinPhi", kTH2F, {{200, -TMath::Pi() / 2, TMath::Pi() / 2}, {200, -1, 1}});

    registry.add("V0LambdaprotonPhi", "V0LambdaprotonPhi", {HistType::kTH1F, {{200, -TMath::Pi() / 2, TMath::Pi() / 2}}});
    registry.add("V0AntiLambdaprotonPhi", "V0AntiLambdaprotonPhi", {HistType::kTH1F, {{200, -TMath::Pi() / 2, TMath::Pi() / 2}}});

    registryData.add("number_of_events_data", "number of events in data", HistType::kTH1D, {{20, 0, 20, "Event Cuts"}});
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

    registryData.add("hMassLambda", "hMassLambda", {HistType::kTH1F, {{200, 0.9f, 1.2f}}});
    registryData.add("hMassAntiLambda", "hMassAntiLambda", {HistType::kTH1F, {{200, 0.9f, 1.2f}}});
    registryData.add("V0pTInLab", "V0pTInLab", kTH1F, {axisPT});

    registryData.add("V0pxInLab", "V0pxInLab", kTH1F, {axisPx});
    registryData.add("V0pyInLab", "V0pyInLab", kTH1F, {axisPy});
    registryData.add("V0pzInLab", "V0pzInLab", kTH1F, {axisPz});

    registryData.add("V0pxInRest_frame", "V0pxInRest_frame", kTH1F, {axisPx});
    registryData.add("V0pyInRest_frame", "V0pyInRest_frame", kTH1F, {axisPy});
    registryData.add("V0pzInRest_frame", "V0pzInRest_frame", kTH1F, {axisPz});

    registryData.add("V0protonpxInLab", "V0protonpxInLab", kTH1F, {axisPx});
    registryData.add("V0protonpyInLab", "V0protonpyInLab", kTH1F, {axisPy});
    registryData.add("V0protonpzInLab", "V0protonpzInLab", kTH1F, {axisPz});
    registryData.add("V0protonphiInLab", "V0protonphiInLab", kTH1F, {axisPhi});

    registryData.add("V0protonpxInRest_frame", "V0protonpxInRest_frame", kTH1F, {axisPx});
    registryData.add("V0protonpyInRest_frame", "V0protonpyInRest_frame", kTH1F, {axisPy});
    registryData.add("V0protonpzInRest_frame", "V0protonpzInRest_frame", kTH1F, {axisPz});
    registryData.add("V0protonMassInRest_frame", "V0protonMassInRest_frame", kTH1F, {axisMass});
    registryData.add("V0protonphiInRest_frame", "V0protonphiInRest_frame", kTH1F, {axisPhi});

    registryData.add("V0protonpxInJetV0frame", "V0protonpxInJetV0frame", kTH1F, {axisPx});
    registryData.add("V0protonpyInJetV0frame", "V0protonpyInJetV0frame", kTH1F, {axisPy});
    registryData.add("V0protonpzInJetV0frame", "V0protonpzInJetV0frame", kTH1F, {axisPz});

    registryData.add("V0LambdapxInJetV0frame", "V0LambdapxInJetV0frame", kTH1F, {axisPx});
    registryData.add("V0LambdapyInJetV0frame", "V0LambdapyInJetV0frame", kTH1F, {axisPy});
    registryData.add("V0LambdapzInJetV0frame", "V0LambdapzInJetV0frame", kTH1F, {axisPz});
    registryData.add("hLambdamassandSinPhi", "hLambdamassandSinPhi", kTH2F, {{200, 0.9, 1.2}, {200, -1, 1}});
    registryData.add("profileLambda", "Invariant Mass vs sin(phi)", {HistType::kTProfile, {{200, 0.9, 1.2}}});
    registryData.add("hLambdaPhiandSinPhi", "hLambdaPhiandSinPhi", kTH2F, {{200, -TMath::Pi() / 2, TMath::Pi() / 2}, {200, -1, 1}});
    registryData.add("V0LambdaprotonPhi", "V0LambdaprotonPhi", {HistType::kTH1F, {{200, -TMath::Pi() / 2, TMath::Pi() / 2}}});

    registryData.add("profileAntiLambda", "Invariant Mass vs sin(phi)", {HistType::kTProfile, {{200, 0.9, 1.2}}});
    registryData.add("hAntiLambdamassandSinPhi", "hAntiLambdaPhiandSinPhi", kTH2F, {{200, -TMath::Pi() / 2, TMath::Pi() / 2}, {200, -1, 1}});

    registryData.add("TProfile2DLambdaPtMassSinPhi", "", kTProfile2D, {TProfile2DaxisMass, TProfile2DaxisPt});
    registryData.add("TProfile2DAntiLambdaPtMassSinPhi", "", kTProfile2D, {TProfile2DaxisMass, TProfile2DaxisPt});
    registryData.add("TProfile2DLambdaPtMassSintheta", "", kTProfile2D, {TProfile2DaxisMass, TProfile2DaxisPt});
    registryData.add("TProfile2DAntiLambdaPtMassSintheta", "", kTProfile2D, {TProfile2DaxisMass, TProfile2DaxisPt});

    registryData.add("TProfile2DLambdaPtMassCosSquareTheta", "", kTProfile2D, {TProfile2DaxisMass, TProfile2DaxisPt});
    registryData.add("TProfile2DAntiLambdaPtMassCosSquareTheta", "", kTProfile2D, {TProfile2DaxisMass, TProfile2DaxisPt});

    registryData.add("hNEvents", "hNEvents", {HistType::kTH1I, {{10, 0.f, 10.f}}});
    registryData.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(1, "all");
    registryData.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(2, "sel8");
    registryData.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(3, "TVX");
    registryData.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(4, "zvertex");
    registryData.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(5, "TFBorder");
    registryData.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(6, "ITSROFBorder");
    registryData.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(7, "isTOFVertexMatched");
    registryData.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(8, "isGoodZvtxFT0vsPV");
    registryData.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(9, "Applied selected");

    registryV0Data.add("hLambdaPt", "hLambdaPt", {HistType::kTH1F, {ptAxis}});
    registryV0Data.add("hAntiLambdaPt", "hAntiLambdaPt", {HistType::kTH1F, {ptAxis}});

    registryV0Data.add("hMassVsPtLambda", "hMassVsPtLambda", {HistType::kTH2F, {ptAxis, invMassLambdaAxis}});
    registryV0Data.add("hMassVsPtAntiLambda", "hMassVsPtAntiLambda", {HistType::kTH2F, {ptAxis, invMassLambdaAxis}});
    registryV0Data.add("hMassLambda", "hMassLambda", {HistType::kTH1F, {invMassLambdaAxis}});
    registryV0Data.add("hMassAntiLambda", "hMassAntiLambda", {HistType::kTH1F, {invMassLambdaAxis}});
    registryV0Data.add("nV0sPerEvent", "nV0sPerEvent", kTH1F, {{10, 0.0, 10.0}});
    registryV0Data.add("nV0sPerEventsel", "nV0sPerEventsel", kTH1F, {{10, 0.0, 10.0}});

    registryV0Data.add("hprotoncosthetainLab", "hprotoncosthetainLab", kTH1F, {{200, -1.f, 1.f}});
    registryV0Data.add("hprotonsinthetainLab", "hprotonsinthetainLab", kTH1F, {{200, -1.f, 1.f}});
    registryV0Data.add("hprotonthetainLab", "hprotonthetainLab", kTH1F, {{200, 0.f, TMath::Pi()}});

    registryV0Data.add("hprotoncosthetainV0", "hprotoncosthetainV0", kTH1F, {{200, -1.f, 1.f}});
    registryV0Data.add("hprotonsinthetainV0", "hprotonsinthetainV0", kTH1F, {{200, -1.f, 1.f}});
    registryV0Data.add("hprotonthetainV0", "hprotonthetainV0", kTH1F, {{200, 0.f, TMath::Pi()}});

    registryV0Data.add("hprotoncosthetainJetV0", "hprotoncosthetainJetV0", kTH1F, {{200, -1.f, 1.f}});
    registryV0Data.add("hprotonsinthetainJetV0", "hprotonsinthetainJetV0", kTH1F, {{200, -1.f, 1.f}});
    registryV0Data.add("hprotonthetainJetV0", "hprotonthetainJetV0", kTH1F, {{200, 0.f, TMath::Pi()}});

    registryV0Data.add("hprotoncosSquarethetainLab", "hprotoncosSquarethetainLab", kTH1F, {{200, -1.f, 1.f}});
    registryV0Data.add("hprotoncosSquarethetainV0", "hprotoncosSquarethetainV0", kTH1F, {{200, -1.f, 1.f}});
    registryV0Data.add("hprotoncosSquarethetainJetV0", "hprotoncosSquarethetainJetV0", kTH1F, {{200, -1.f, 1.f}});

    registryV0Data.add("hLambdamassandSinthetainV0", "hLambdamassandSinthetainV0", kTH2F, {{200, 0.9, 1.2}, {200, -1, 1}});
    registryV0Data.add("hLambdamassandCosthetainV0", "hLambdamassandCosthetainV0", kTH2F, {{200, 0.9, 1.2}, {200, -1, 1}});
    registryV0Data.add("hLambdamassandCosSquarethetainV0", "hLambdamassandCosSquarethetainV0", kTH2F, {{200, 0.9, 1.2}, {200, -1, 1}});

    registryV0Data.add("hLambdamassandSinthetainJetV0", "hLambdamassandSinthetainJetV0", kTH2F, {{200, 0.9, 1.2}, {200, -1, 1}});
    registryV0Data.add("hLambdamassandCosthetainJetV0", "hLambdamassandCosthetainJetV0", kTH2F, {{200, 0.9, 1.2}, {200, -1, 1}});
    registryV0Data.add("hLambdamassandCosSquarethetainJetV0", "hLambdamassandCosSquarethetainJetV0", kTH2F, {{200, 0.9, 1.2}, {200, -1, 1}});

    registryV0Data.add("AverageSinthetainV0", "AverageSinthetainV0", {HistType::kTProfile, {{200, 0.9, 1.2}}});
    registryV0Data.add("AverageCosSquarethetainV0", "AverageCosSquarethetainV0", {HistType::kTProfile, {{200, 0.9, 1.2}}});

    registryV0Data.add("AverageSinthetainJetV0", "AverageSinthetainJetV0", {HistType::kTProfile, {{200, 0.9, 1.2}}});
    registryV0Data.add("AverageCosSquarethetainJetV0", "AverageCosSquarethetainJetV0", {HistType::kTProfile, {{200, 0.9, 1.2}}});

    // LongitudinalPolarization event selection
    registryLongitudinalPolarization.add("hNEvents", "hNEvents", {HistType::kTH1D, {{10, 0.f, 10.f}}});
    registryLongitudinalPolarization.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(1, "all");
    registryLongitudinalPolarization.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(2, "sel8");
    registryLongitudinalPolarization.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(3, "TVX");
    registryLongitudinalPolarization.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(4, "zvertex");
    registryLongitudinalPolarization.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(5, "TFBorder");
    registryLongitudinalPolarization.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(6, "ITSROFBorder");
    registryLongitudinalPolarization.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(7, "isTOFVertexMatched");
    registryLongitudinalPolarization.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(8, "isNoSameBunchPileup");
    registryLongitudinalPolarization.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(9, "Applied selection");

    registryLongitudinalPolarization.add("hMassVsPtLambda", "hMassVsPtLambda", {HistType::kTH2F, {ptAxis, invMassLambdaAxis}});
    registryLongitudinalPolarization.add("hMassLambda", "hMassLambda", {HistType::kTH1F, {invMassLambdaAxis}});
    registryLongitudinalPolarization.add("hMassVsPtAntiLambda", "hMassVsPtAntiLambda", {HistType::kTH2F, {ptAxis, invMassLambdaAxis}});

    registryLongitudinalPolarization.add("hMassLambdasel", "hMassLambdasel", {HistType::kTH1F, {invMassLambdaAxis}});
    registryLongitudinalPolarization.add("hMassAntiLambdasel", "hMassAntiLambdasel", {HistType::kTH1F, {invMassLambdaAxis}});
    registryLongitudinalPolarization.add("hMassVsPtLambdasel", "hMassVsPtLambdasel", {HistType::kTH2F, {ptAxis, invMassLambdaAxis}});
    registryLongitudinalPolarization.add("hMassVsPtAntiLambdasel", "hMassVsPtAntiLambdasel", {HistType::kTH2F, {ptAxis, invMassLambdaAxis}});

    registryLongitudinalPolarization.add("V0pxInRest_frame", "V0pxInRest_frame", kTH1F, {axisPx});
    registryLongitudinalPolarization.add("V0pyInRest_frame", "V0pyInRest_frame", kTH1F, {axisPy});
    registryLongitudinalPolarization.add("V0pzInRest_frame", "V0pzInRest_frame", kTH1F, {axisPz});

    registryLongitudinalPolarization.add("nV0sPerEvent", "nV0sPerEvent", kTH1F, {{10, 0.0, 10.0}});
    registryLongitudinalPolarization.add("nV0sPerEventsel", "nV0sPerEventsel", kTH1F, {{10, 0.0, 10.0}});

    registryLongitudinalPolarization.add("hprotoncosthetainV0", "hprotoncosthetainV0", kTH1F, {{200, -1.f, 1.f}});
    registryLongitudinalPolarization.add("hprotoncosSquarethetainV0", "hprotoncosSquarethetainV0", kTH1F, {{200, -1.f, 1.f}});
    registryLongitudinalPolarization.add("hLambdamassandCosthetaInV0", "hLambdamassandCosthetaInV0", kTH2F, {{200, 0.9, 1.2}, {200, -1, 1}});
    registryLongitudinalPolarization.add("TProfile2DLambdaPtMassCostheta", "", kTProfile2D, {TProfile2DaxisMass, TProfile2DaxisPt});
    registryLongitudinalPolarization.add("TProfile2DLambdaPtMassCosSquareTheta", "", kTProfile2D, {TProfile2DaxisMass, TProfile2DaxisPt});

    registryLongitudinalPolarization.add("hantiprotoncosthetainV0", "hantiprotoncosthetainV0", kTH1F, {{200, -1.f, 1.f}});
    registryLongitudinalPolarization.add("hantiprotoncosSquarethetainV0", "hantiprotoncosSquarethetainV0", kTH1F, {{200, -1.f, 1.f}});
    registryLongitudinalPolarization.add("hAntiLambdamassandCosthetaInV0", "hAntiLambdamassandCosthetaInV0", kTH2F, {{200, 0.9, 1.2}, {200, -1, 1}});
    registryLongitudinalPolarization.add("TProfile2DAntiLambdaPtMassCostheta", "TProfile2DAntiLambdaPtMassCostheta", kTProfile2D, {TProfile2DaxisMass, TProfile2DaxisPt});
    registryLongitudinalPolarization.add("TProfile2DAntiLambdaPtMassCosSquareTheta", "TProfile2DAntiLambdaPtMassCosSquareTheta", kTProfile2D, {TProfile2DaxisMass, TProfile2DaxisPt});
    registryLongitudinalPolarization.add("TProfile1DLambdaPtMassCostheta", "Invariant Mass vs cos(#theta)", {HistType::kTProfile, {{200, 0.9, 1.2}}});
    registryLongitudinalPolarization.add("TProfile1DAntiLambdaPtMassCostheta", "Invariant Mass vs cos(#theta)", {HistType::kTProfile, {{200, 0.9, 1.2}}});

    if (doQA) {
      registryLongitudinalPolarization.add("QA/hv0sSelection", ";Sel", {HistType::kTH1D, {{22, 0., 22.}}});
      registryLongitudinalPolarization.get<TH1>(HIST("QA/hv0sSelection"))->GetXaxis()->SetBinLabel(1, "all");
      registryLongitudinalPolarization.get<TH1>(HIST("QA/hv0sSelection"))->GetXaxis()->SetBinLabel(2, "Event selection");
      registryLongitudinalPolarization.get<TH1>(HIST("QA/hv0sSelection"))->GetXaxis()->SetBinLabel(3, "Radius");
      registryLongitudinalPolarization.get<TH1>(HIST("QA/hv0sSelection"))->GetXaxis()->SetBinLabel(4, "Eta Daughters");
      registryLongitudinalPolarization.get<TH1>(HIST("QA/hv0sSelection"))->GetXaxis()->SetBinLabel(5, "Dau DCA to PV");
      registryLongitudinalPolarization.get<TH1>(HIST("QA/hv0sSelection"))->GetXaxis()->SetBinLabel(6, "DCA Daughters");
      registryLongitudinalPolarization.get<TH1>(HIST("QA/hv0sSelection"))->GetXaxis()->SetBinLabel(7, "min ITS hits");
      registryLongitudinalPolarization.get<TH1>(HIST("QA/hv0sSelection"))->GetXaxis()->SetBinLabel(8, "has TOF 1 Leg");
      registryLongitudinalPolarization.get<TH1>(HIST("QA/hv0sSelection"))->GetXaxis()->SetBinLabel(9, "has TOF 2 Legs");
      registryLongitudinalPolarization.get<TH1>(HIST("QA/hv0sSelection"))->GetXaxis()->SetBinLabel(10, "TPC NCl");
      registryLongitudinalPolarization.get<TH1>(HIST("QA/hv0sSelection"))->GetXaxis()->SetBinLabel(11, "TPC Cls Shared");
      registryLongitudinalPolarization.get<TH1>(HIST("QA/hv0sSelection"))->GetXaxis()->SetBinLabel(12, "ITS Chi2");
      registryLongitudinalPolarization.get<TH1>(HIST("QA/hv0sSelection"))->GetXaxis()->SetBinLabel(13, "TPC Chi2");
      registryLongitudinalPolarization.get<TH1>(HIST("QA/hv0sSelection"))->GetXaxis()->SetBinLabel(14, "cosPA");
      registryLongitudinalPolarization.get<TH1>(HIST("QA/hv0sSelection"))->GetXaxis()->SetBinLabel(15, "rapidity");
      registryLongitudinalPolarization.get<TH1>(HIST("QA/hv0sSelection"))->GetXaxis()->SetBinLabel(16, "ctau");
      registryLongitudinalPolarization.get<TH1>(HIST("QA/hv0sSelection"))->GetXaxis()->SetBinLabel(17, "v0 rej");
      registryLongitudinalPolarization.get<TH1>(HIST("QA/hv0sSelection"))->GetXaxis()->SetBinLabel(18, "TPC nsigma Neg Dau");
      registryLongitudinalPolarization.get<TH1>(HIST("QA/hv0sSelection"))->GetXaxis()->SetBinLabel(19, "TPC nsigma Pos Dau");
      registryLongitudinalPolarization.get<TH1>(HIST("QA/hv0sSelection"))->GetXaxis()->SetBinLabel(20, "Armenteros-Podolansky");
    }
  }
  double massPr = o2::constants::physics::MassProton;
  double massLambda = o2::constants::physics::MassLambda;
  double massPi = o2::constants::physics::MassPionCharged;
  ROOT::Math::PxPyPzMVector ProtonVec, PionVec, LambdaVec, ProtonBoostedVec, LambdaBoostedVec;

  TMatrixD LorentzTransInV0frame(double ELambda, double Lambdapx, double Lambdapy, double Lambdapz)
  {
    double PLambda = sqrt(Lambdapx * Lambdapx + Lambdapy * Lambdapy + Lambdapz * Lambdapz);
    double LambdaMass = sqrt(ELambda * ELambda - PLambda * PLambda);
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
  // aod::MyCollision const& collision

  void processJetV0Analysis(aod::MyTable const& myv0s, aod::MyTableJet const& myJets)
  {

    for (auto& candidate : myv0s) {
      registry.fill(HIST("hMassLambda"), candidate.v0Lambdamass());
      registry.fill(HIST("V0pTInLab"), candidate.v0pt());
      registry.fill(HIST("hMassVsPtLambda"), candidate.v0pt(), candidate.v0Lambdamass());
      registry.fill(HIST("V0pxInLab"), candidate.v0px());
      registry.fill(HIST("V0pyInLab"), candidate.v0py());
      registry.fill(HIST("V0pzInLab"), candidate.v0pz());
      registry.fill(HIST("V0protonpxInLab"), candidate.v0protonpx());
      registry.fill(HIST("V0protonpyInLab"), candidate.v0protonpy());
      registry.fill(HIST("V0protonpzInLab"), candidate.v0protonpz());
      double protonsinPhiInLab = candidate.v0protonpy() / sqrt(candidate.v0protonpx() * candidate.v0protonpx() + candidate.v0protonpy() * candidate.v0protonpy());
      registry.fill(HIST("V0protonphiInLab"), protonsinPhiInLab);
      double PLambda = sqrt(candidate.v0px() * candidate.v0px() + candidate.v0py() * candidate.v0py() + candidate.v0pz() * candidate.v0pz());
      double ELambda = sqrt(candidate.v0Lambdamass() * candidate.v0Lambdamass() + PLambda * PLambda);
      TMatrixD pLabV0(4, 1);
      pLabV0(0, 0) = ELambda;
      pLabV0(1, 0) = candidate.v0px();
      pLabV0(2, 0) = candidate.v0py();
      pLabV0(3, 0) = candidate.v0pz();
      TMatrixD V0InV0(4, 1);
      V0InV0 = LorentzTransInV0frame(ELambda, candidate.v0px(), candidate.v0py(), candidate.v0pz()) * pLabV0;
      registry.fill(HIST("V0pxInRest_frame"), V0InV0(1, 0));
      registry.fill(HIST("V0pyInRest_frame"), V0InV0(2, 0));
      registry.fill(HIST("V0pzInRest_frame"), V0InV0(3, 0));
    }
    for (auto& candidate : myv0s) {
      double PLambda = sqrt(candidate.v0px() * candidate.v0px() + candidate.v0py() * candidate.v0py() + candidate.v0pz() * candidate.v0pz());
      double ELambda = sqrt(candidate.v0Lambdamass() * candidate.v0Lambdamass() + PLambda * PLambda);
      TMatrixD pLabproton(4, 1);
      double protonE = sqrt(massPr * massPr + candidate.v0protonpx() * candidate.v0protonpx() + candidate.v0protonpy() * candidate.v0protonpy() + candidate.v0protonpz() * candidate.v0protonpz());
      pLabproton(0, 0) = protonE;
      pLabproton(1, 0) = candidate.v0protonpx();
      pLabproton(2, 0) = candidate.v0protonpy();
      pLabproton(3, 0) = candidate.v0protonpz();
      TMatrixD protonInV0(4, 1);
      protonInV0 = LorentzTransInV0frame(ELambda, candidate.v0px(), candidate.v0py(), candidate.v0pz()) * pLabproton;
      double protonMassInV0 = sqrt(protonInV0(0, 0) * protonInV0(0, 0) - protonInV0(1, 0) * protonInV0(1, 0) - protonInV0(2, 0) * protonInV0(2, 0) - protonInV0(3, 0) * protonInV0(3, 0));
      registry.fill(HIST("V0protonMassInRest_frame"), protonMassInV0);
      registry.fill(HIST("V0protonpxInRest_frame"), protonInV0(1, 0));
      registry.fill(HIST("V0protonpyInRest_frame"), protonInV0(2, 0));
      registry.fill(HIST("V0protonpzInRest_frame"), protonInV0(3, 0));
      double protonsinPhiInV0frame = protonInV0(2, 0) / sqrt(protonInV0(1, 0) * protonInV0(1, 0) + protonInV0(2, 0) * protonInV0(2, 0));
      registry.fill(HIST("V0protonphiInRest_frame"), protonsinPhiInV0frame);
    }

    for (auto& Jet : myJets) {
      registry.fill(HIST("JetpxInLab"), Jet.jetpx());
      registry.fill(HIST("JetpyInLab"), Jet.jetpy());
      registry.fill(HIST("JetpzInLab"), Jet.jetpz());
      registry.fill(HIST("JetpTInLab"), Jet.jetpt());
    }
  }
  PROCESS_SWITCH(LfMyV0s, processJetV0Analysis, "processJetV0Analysis", true);
  void processLeadingJetV0Analysis(aod::MyTable const& myv0s, aod::MyTableLeadingJet const& myleadingJets)
  {
    for (auto& LeadingJet : myleadingJets) {
      int V0Numbers = 0;
      double protonsinPhiInJetV0frame = 0;
      for (auto& candidate : myv0s) {
        if (candidate.mycollisionv0() == LeadingJet.mycollisionleadingjet()) {
          V0Numbers = V0Numbers + 1;
          double PLambda = sqrt(candidate.v0px() * candidate.v0px() + candidate.v0py() * candidate.v0py() + candidate.v0pz() * candidate.v0pz());
          double ELambda = sqrt(candidate.v0Lambdamass() * candidate.v0Lambdamass() + PLambda * PLambda);
          double protonE = sqrt(massPr * massPr + candidate.v0protonpx() * candidate.v0protonpx() + candidate.v0protonpy() * candidate.v0protonpy() + candidate.v0protonpz() * candidate.v0protonpz());

          TMatrixD pLabV0(4, 1);
          pLabV0(0, 0) = ELambda;
          pLabV0(1, 0) = candidate.v0px();
          pLabV0(2, 0) = candidate.v0py();
          pLabV0(3, 0) = candidate.v0pz();

          TMatrixD lambdaInJet(4, 1);
          lambdaInJet = MyTMatrixTranslationToJet(LeadingJet.leadingjetpx(), LeadingJet.leadingjetpy(), LeadingJet.leadingjetpz(), candidate.v0px(), candidate.v0py(), candidate.v0pz()) * pLabV0;

          TMatrixD lambdaInJetV0(4, 1);
          lambdaInJetV0 = LorentzTransInV0frame(ELambda, lambdaInJet(1, 0), lambdaInJet(2, 0), lambdaInJet(3, 0)) * MyTMatrixTranslationToJet(LeadingJet.leadingjetpx(), LeadingJet.leadingjetpy(), LeadingJet.leadingjetpz(), candidate.v0px(), candidate.v0py(), candidate.v0pz()) * pLabV0;
          registry.fill(HIST("V0LambdapxInJetV0frame"), lambdaInJetV0(1, 0));
          registry.fill(HIST("V0LambdapyInJetV0frame"), lambdaInJetV0(2, 0));
          registry.fill(HIST("V0LambdapzInJetV0frame"), lambdaInJetV0(3, 0));

          TMatrixD pLabproton(4, 1);
          pLabproton(0, 0) = protonE;
          pLabproton(1, 0) = candidate.v0protonpx();
          pLabproton(2, 0) = candidate.v0protonpy();
          pLabproton(3, 0) = candidate.v0protonpz();
          TMatrixD protonInJetV0(4, 1);
          protonInJetV0 = LorentzTransInV0frame(ELambda, lambdaInJet(1, 0), lambdaInJet(2, 0), lambdaInJet(3, 0)) * MyTMatrixTranslationToJet(LeadingJet.leadingjetpx(), LeadingJet.leadingjetpy(), LeadingJet.leadingjetpz(), candidate.v0px(), candidate.v0py(), candidate.v0pz()) * pLabproton;
          registry.fill(HIST("V0protonpxInJetV0frame"), protonInJetV0(1, 0));
          registry.fill(HIST("V0protonpyInJetV0frame"), protonInJetV0(2, 0));
          registry.fill(HIST("V0protonpzInJetV0frame"), protonInJetV0(3, 0));
          protonsinPhiInJetV0frame = protonsinPhiInJetV0frame + protonInJetV0(2, 0) / sqrt(protonInJetV0(1, 0) * protonInJetV0(1, 0) + protonInJetV0(2, 0) * protonInJetV0(2, 0));
        }
      }
      for (auto& candidate : myv0s) {
        if (candidate.mycollisionv0() == LeadingJet.mycollisionleadingjet()) {
          registry.fill(HIST("V0protonphiInJetV0frame"), protonsinPhiInJetV0frame / V0Numbers);
          registry.fill(HIST("hLambdamassandSinPhi"), candidate.v0Lambdamass(), protonsinPhiInJetV0frame / V0Numbers);
          registry.fill(HIST("hLambdaPhiandSinPhi"), TMath::ASin(protonsinPhiInJetV0frame / V0Numbers), protonsinPhiInJetV0frame / V0Numbers);
          registry.fill(HIST("V0LambdaprotonPhi"), TMath::ASin(protonsinPhiInJetV0frame / V0Numbers));
          registry.fill(HIST("profile"), candidate.v0Lambdamass(), protonsinPhiInJetV0frame / V0Numbers);
        }
      }
    }
    for (auto& LeadingJet : myleadingJets) {
      registry.fill(HIST("LeadingJetpx"), LeadingJet.leadingjetpx());
      registry.fill(HIST("LeadingJetpy"), LeadingJet.leadingjetpy());
      registry.fill(HIST("LeadingJetpz"), LeadingJet.leadingjetpz());
      registry.fill(HIST("LeadingJetpT"), LeadingJet.leadingjetpt());
    }
  }
  PROCESS_SWITCH(LfMyV0s, processLeadingJetV0Analysis, "processLeadingJetV0Analysis", true);

  void processLeadingJetAntiV0Analysis(aod::MyTableAnti const& myv0s, aod::MyTableLeadingJet const& myleadingJets)
  {
    for (auto& LeadingJet : myleadingJets) {
      int V0Numbers = 0;
      double protonsinPhiInJetV0frame = 0;
      for (auto& candidate : myv0s) {
        if (candidate.mycollisionv0() == LeadingJet.mycollisionleadingjet()) {
          V0Numbers = V0Numbers + 1;
          double PLambda = sqrt(candidate.v0px() * candidate.v0px() + candidate.v0py() * candidate.v0py() + candidate.v0pz() * candidate.v0pz());
          double ELambda = sqrt(candidate.v0Lambdamass() * candidate.v0Lambdamass() + PLambda * PLambda);
          double protonE = sqrt(massPr * massPr + candidate.v0protonpx() * candidate.v0protonpx() + candidate.v0protonpy() * candidate.v0protonpy() + candidate.v0protonpz() * candidate.v0protonpz());

          TMatrixD pLabV0(4, 1);
          pLabV0(0, 0) = ELambda;
          pLabV0(1, 0) = candidate.v0px();
          pLabV0(2, 0) = candidate.v0py();
          pLabV0(3, 0) = candidate.v0pz();

          TMatrixD lambdaInJet(4, 1);
          lambdaInJet = MyTMatrixTranslationToJet(LeadingJet.leadingjetpx(), LeadingJet.leadingjetpy(), LeadingJet.leadingjetpz(), candidate.v0px(), candidate.v0py(), candidate.v0pz()) * pLabV0;

          TMatrixD lambdaInJetV0(4, 1);
          lambdaInJetV0 = LorentzTransInV0frame(ELambda, lambdaInJet(1, 0), lambdaInJet(2, 0), lambdaInJet(3, 0)) * MyTMatrixTranslationToJet(LeadingJet.leadingjetpx(), LeadingJet.leadingjetpy(), LeadingJet.leadingjetpz(), candidate.v0px(), candidate.v0py(), candidate.v0pz()) * pLabV0;

          TMatrixD pLabproton(4, 1);
          pLabproton(0, 0) = protonE;
          pLabproton(1, 0) = candidate.v0protonpx();
          pLabproton(2, 0) = candidate.v0protonpy();
          pLabproton(3, 0) = candidate.v0protonpz();
          TMatrixD protonInJetV0(4, 1);
          protonInJetV0 = LorentzTransInV0frame(ELambda, lambdaInJet(1, 0), lambdaInJet(2, 0), lambdaInJet(3, 0)) * MyTMatrixTranslationToJet(LeadingJet.leadingjetpx(), LeadingJet.leadingjetpy(), LeadingJet.leadingjetpz(), candidate.v0px(), candidate.v0py(), candidate.v0pz()) * pLabproton;
          protonsinPhiInJetV0frame = protonsinPhiInJetV0frame + protonInJetV0(2, 0) / sqrt(protonInJetV0(1, 0) * protonInJetV0(1, 0) + protonInJetV0(2, 0) * protonInJetV0(2, 0));
        }
      }
      for (auto& candidate : myv0s) {
        if (candidate.mycollisionv0() == LeadingJet.mycollisionleadingjet()) {
          registry.fill(HIST("V0antiprotonphiInJetV0frame"), protonsinPhiInJetV0frame / V0Numbers);
          registry.fill(HIST("hAntiLambdamassandSinPhi"), candidate.v0Lambdamass(), protonsinPhiInJetV0frame / V0Numbers);
          registry.fill(HIST("hAntiLambdaPhiandSinPhi"), TMath::ASin(protonsinPhiInJetV0frame / V0Numbers), protonsinPhiInJetV0frame / V0Numbers);
          registry.fill(HIST("V0AntiLambdaprotonPhi"), TMath::ASin(protonsinPhiInJetV0frame / V0Numbers));
          registry.fill(HIST("profileAntiV0"), candidate.v0Lambdamass(), protonsinPhiInJetV0frame / V0Numbers);
        }
      }
    }
    for (auto& candidate : myv0s) {
      registry.fill(HIST("hMassVsPtAntiLambda"), candidate.v0pt(), candidate.v0Lambdamass());
    }
  }
  PROCESS_SWITCH(LfMyV0s, processLeadingJetAntiV0Analysis, "processLeadingJetAntiV0Analysis", true);

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
        TMath::Abs(ptrack.eta()) > V0tracketaMax ||
        TMath::Abs(ntrack.eta()) > V0tracketaMax) {
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

    registryLongitudinalPolarization.fill(HIST("QA/hv0sSelection"), 0.5);

    if (evSel && evFlag < 1)
      return false;

    registryLongitudinalPolarization.fill(HIST("QA/hv0sSelection"), 1.5);

    if (v0.v0radius() < minimumV0Radius || v0.v0radius() > maximumV0Radius)
      return false;

    registryLongitudinalPolarization.fill(HIST("QA/hv0sSelection"), 2.5);

    if (TMath::Abs(ptrack.eta()) > V0tracketaMax || TMath::Abs(ntrack.eta()) > V0tracketaMax) {
      return false;
    }
    registryLongitudinalPolarization.fill(HIST("QA/hv0sSelection"), 3.5);

    if (std::fabs(v0.dcapostopv()) < dcapostoPVmin)
      return false;
    if (std::fabs(v0.dcanegtopv()) < dcanegtoPVmin)
      return false;
    registryLongitudinalPolarization.fill(HIST("QA/hv0sSelection"), 4.5);

    if (std::fabs(v0.dcaV0daughters()) > dcaV0DaughtersMax)
      return false;
    registryLongitudinalPolarization.fill(HIST("QA/hv0sSelection"), 5.5);

    if (requireITS && ptrack.itsNCls() < minITSnCls)
      return false;
    if (requireITS && ntrack.itsNCls() < minITSnCls)
      return false;
    registryLongitudinalPolarization.fill(HIST("QA/hv0sSelection"), 6.5);

    if (hasTOF1Leg && !ptrack.hasTOF() && !ntrack.hasTOF())
      return false;
    registryLongitudinalPolarization.fill(HIST("QA/hv0sSelection"), 7.5);

    if (hasTOF2Leg && (!ptrack.hasTOF() || !ntrack.hasTOF()))
      return false;
    registryLongitudinalPolarization.fill(HIST("QA/hv0sSelection"), 8.5);

    if (ptrack.tpcNClsCrossedRows() < minNCrossedRowsTPC)
      return false;
    if (ntrack.tpcNClsCrossedRows() < minNCrossedRowsTPC)
      return false;
    registryLongitudinalPolarization.fill(HIST("QA/hv0sSelection"), 9.5);

    if (ptrack.tpcNClsShared() > max_tpcSharedCls)
      return false;
    if (ntrack.tpcNClsShared() > max_tpcSharedCls)
      return false;
    registryLongitudinalPolarization.fill(HIST("QA/hv0sSelection"), 10.5);

    if (ptrack.itsChi2NCl() > max_chi2_ITS)
      return false;
    if (ntrack.itsChi2NCl() > max_chi2_ITS)
      return false;
    registryLongitudinalPolarization.fill(HIST("QA/hv0sSelection"), 11.5);

    if (ptrack.tpcChi2NCl() > max_chi2_TPC)
      return false;
    if (ntrack.tpcChi2NCl() > max_chi2_TPC)
      return false;
    registryLongitudinalPolarization.fill(HIST("QA/hv0sSelection"), 12.5);

    if (v0.v0cosPA() < v0cospaMin)
      return false;
    registryLongitudinalPolarization.fill(HIST("QA/hv0sSelection"), 13.5);

    if (v0.yLambda() < yMin || v0.yLambda() > yMax) {
      return false;
    }
    registryLongitudinalPolarization.fill(HIST("QA/hv0sSelection"), 14.5);

    float ctauLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0;
    if (ctauLambda >= CtauLambda)
      return false;
    registryLongitudinalPolarization.fill(HIST("QA/hv0sSelection"), 15.5);

    if (TMath::Abs(v0.mK0Short() - o2::constants::physics::MassK0Short) < v0rejLambda) {
      return false;
    }
    if (std::abs(v0.mLambda() - o2::constants::physics::MassLambda0) > v0accLambda) {
      return false;
    }

    registryLongitudinalPolarization.fill(HIST("QA/hv0sSelection"), 16.5);

    if (ntrack.tpcNSigmaPi() < nsigmaTPCmin || ntrack.tpcNSigmaPi() > nsigmaTPCmax)
      return false;
    registryLongitudinalPolarization.fill(HIST("QA/hv0sSelection"), 17.5);

    if (ptrack.tpcNSigmaPr() < nsigmaTPCmin || ptrack.tpcNSigmaPr() > nsigmaTPCmax)
      return false;
    registryLongitudinalPolarization.fill(HIST("QA/hv0sSelection"), 18.5);

    if (doArmenterosCut && v0.qtarm() > (paramArmenterosCut * std::abs(v0.alpha())))
      return false;
    registryLongitudinalPolarization.fill(HIST("QA/hv0sSelection"), 19.5);

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

    if (TMath::Abs(ptrack.eta()) > V0tracketaMax || TMath::Abs(ntrack.eta()) > V0tracketaMax) {
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

    if (ptrack.tpcNClsShared() > max_tpcSharedCls)
      return false;
    if (ntrack.tpcNClsShared() > max_tpcSharedCls)
      return false;

    if (ptrack.itsChi2NCl() > max_chi2_ITS)
      return false;
    if (ntrack.itsChi2NCl() > max_chi2_ITS)
      return false;

    if (ptrack.tpcChi2NCl() > max_chi2_TPC)
      return false;
    if (ntrack.tpcChi2NCl() > max_chi2_TPC)
      return false;

    if (v0.v0cosPA() < v0cospaMin)
      return false;

    if (v0.yLambda() < yMin || v0.yLambda() > yMax) {
      return false;
    }

    float ctauAntiLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0Bar;
    if (ctauAntiLambda >= CtauLambda)
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

    registryLongitudinalPolarization.fill(HIST("QA/hv0sSelection"), 1.5);

    if (v0.v0radius() < minimumV0Radius || v0.v0radius() > maximumV0Radius)
      return false;

    registryLongitudinalPolarization.fill(HIST("QA/hv0sSelection"), 2.5);

    if (TMath::Abs(ptrack.eta()) > V0tracketaMax || TMath::Abs(ntrack.eta()) > V0tracketaMax) {
      return false;
    }
    registryLongitudinalPolarization.fill(HIST("QA/hv0sSelection"), 3.5);

    if (std::fabs(v0.dcapostopv()) < dcapostoPVmin)
      return false;
    if (std::fabs(v0.dcanegtopv()) < dcanegtoPVmin)
      return false;
    registryLongitudinalPolarization.fill(HIST("QA/hv0sSelection"), 4.5);

    if (std::fabs(v0.dcaV0daughters()) > dcaV0DaughtersMax)
      return false;
    registryLongitudinalPolarization.fill(HIST("QA/hv0sSelection"), 5.5);

    if (requireITS && ptrack.itsNCls() < minITSnCls)
      return false;
    if (requireITS && ntrack.itsNCls() < minITSnCls)
      return false;
    registryLongitudinalPolarization.fill(HIST("QA/hv0sSelection"), 6.5);

    if (hasTOF1Leg && !ptrack.hasTOF() && !ntrack.hasTOF())
      return false;
    registryLongitudinalPolarization.fill(HIST("QA/hv0sSelection"), 7.5);

    if (hasTOF2Leg && (!ptrack.hasTOF() || !ntrack.hasTOF()))
      return false;
    registryLongitudinalPolarization.fill(HIST("QA/hv0sSelection"), 8.5);

    if (ptrack.tpcNClsCrossedRows() < minNCrossedRowsTPC)
      return false;
    if (ntrack.tpcNClsCrossedRows() < minNCrossedRowsTPC)
      return false;
    registryLongitudinalPolarization.fill(HIST("QA/hv0sSelection"), 9.5);

    if (ptrack.tpcNClsShared() > max_tpcSharedCls)
      return false;
    if (ntrack.tpcNClsShared() > max_tpcSharedCls)
      return false;
    registryLongitudinalPolarization.fill(HIST("QA/hv0sSelection"), 10.5);

    if (ptrack.itsChi2NCl() > max_chi2_ITS)
      return false;
    if (ntrack.itsChi2NCl() > max_chi2_ITS)
      return false;
    registryLongitudinalPolarization.fill(HIST("QA/hv0sSelection"), 11.5);

    if (ptrack.tpcChi2NCl() > max_chi2_TPC)
      return false;
    if (ntrack.tpcChi2NCl() > max_chi2_TPC)
      return false;
    registryLongitudinalPolarization.fill(HIST("QA/hv0sSelection"), 12.5);

    if (v0.v0cosPA() < v0cospaMin)
      return false;
    registryLongitudinalPolarization.fill(HIST("QA/hv0sSelection"), 13.5);

    if (v0.yLambda() < yMin || v0.yLambda() > yMax) {
      return false;
    }
    registryLongitudinalPolarization.fill(HIST("QA/hv0sSelection"), 14.5);

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
    TLorentzVector lorentzVect;
    lorentzVect.SetXYZM(v0.px(), v0.py(), v0.pz(), 1.115683);
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

    if (TMath::Abs(ptrack.eta()) > V0tracketaMax || TMath::Abs(ntrack.eta()) > V0tracketaMax) {
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
    TLorentzVector lorentzVect;
    lorentzVect.SetXYZM(v0.px(), v0.py(), v0.pz(), 1.115683);
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
    registryLongitudinalPolarization.fill(HIST("hNEvents"), 0.5);
    if (sel8 && !collision.sel8()) {
      return false;
    }
    registryLongitudinalPolarization.fill(HIST("hNEvents"), 1.5);
    if (isTriggerTVX && !collision.selection_bit(aod::evsel::kIsTriggerTVX)) {
      return false;
    }
    registryLongitudinalPolarization.fill(HIST("hNEvents"), 2.5);
    if (std::abs(collision.posZ()) > cutzvertex) {
      return false;
    }
    registryLongitudinalPolarization.fill(HIST("hNEvents"), 3.5);
    if (isNoTimeFrameBorder && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    registryLongitudinalPolarization.fill(HIST("hNEvents"), 4.5);
    if (isNoITSROFrameBorder && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    registryLongitudinalPolarization.fill(HIST("hNEvents"), 5.5);
    if (isVertexTOFmatched && !collision.selection_bit(aod::evsel::kIsVertexTOFmatched)) {
      return false;
    }
    registryLongitudinalPolarization.fill(HIST("hNEvents"), 6.5);
    if (isNoSameBunchPileup && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    registryLongitudinalPolarization.fill(HIST("hNEvents"), 7.5);

    return true;
  }

  using SelCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;
  using StrHadronDaughterTracks = soa::Join<aod::Tracks, aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
  void processData(SelCollisions::iterator const& collision, aod::V0Datas const& fullV0s, StrHadronDaughterTracks const& tracks)
  {
    registryData.fill(HIST("number_of_events_data"), 0.5);
    // event selection
    if (!collision.sel8() || std::fabs(collision.posZ()) > zVtx) {
      return;
    }
    // event counter: after event selection
    registryData.fill(HIST("number_of_events_data"), 1.5);
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
    registryData.fill(HIST("number_of_events_data"), 2.5);

    // cluster particles using the anti-kt algorithm
    fastjet::RecombinationScheme recombScheme = fastjet::E_scheme;
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, rJet, recombScheme);
    fastjet::AreaDefinition areaDef(fastjet::active_area, fastjet::GhostedAreaSpec(1.0));
    fastjet::ClusterSequenceArea cs(fjParticles, jetDef, areaDef);
    std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());
    auto [rhoPerp, rhoMPerp] = backgroundSub.estimateRhoPerpCone(fjParticles, jets);

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
    float maxJetpT = 0;
    float maxJetPt = -999;
    for (auto& jet : jets) {
      nJets++;
      registryData.fill(HIST("FJetaHistogram"), jet.eta());
      registryData.fill(HIST("FJphiHistogram"), jet.phi());
      registryData.fill(HIST("FJptHistogram"), jet.pt());
      // jet must be fully contained in the acceptance
      fastjet::PseudoJet jetMinusBkg = backgroundSub.doRhoAreaSub(jet, rhoPerp, rhoMPerp);
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
    registryData.fill(HIST("number_of_events_data"), 3.5);
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
      if (passedLambdaSelection(v0, pos, neg)) {
        V0Numbers = V0Numbers + 1;
        registryData.fill(HIST("LambdaPtMass"), v0.pt(), v0.mLambda());
      }
      if (passedAntiLambdaSelection(v0, pos, neg)) {
        AntiV0Numbers = AntiV0Numbers + 1;
        registryData.fill(HIST("AntiLambdaPtMass"), v0.pt(), v0.mAntiLambda());
      }
    }
    registryData.fill(HIST("nV0sPerEvent"), V0Numbers);

    // calculate lambda polarization induced by jet

    if (V0Numbers == 0) {
      return;
    }
    if (maxJetpx == 0) {
      return;
    }
    double protonsinPhiInJetV0frame = 0;
    double AntiprotonsinPhiInJetV0frame = 0;
    cout << maxJetpx << endl;
    for (const auto& candidate : fullV0s) {
      const auto& pos = candidate.posTrack_as<StrHadronDaughterTracks>();
      const auto& neg = candidate.negTrack_as<StrHadronDaughterTracks>();
      TVector3 v0dir(candidate.px(), candidate.py(), candidate.pz());

      if (passedLambdaSelection(candidate, pos, neg)) {
        registryData.fill(HIST("hMassLambda"), candidate.mLambda());
        registryData.fill(HIST("V0pTInLab"), candidate.pt());
        registryData.fill(HIST("V0pxInLab"), candidate.px());
        registryData.fill(HIST("V0pyInLab"), candidate.py());
        registryData.fill(HIST("V0pzInLab"), candidate.pz());
        registryData.fill(HIST("V0protonpxInLab"), pos.px());
        registryData.fill(HIST("V0protonpyInLab"), pos.py());
        registryData.fill(HIST("V0protonpzInLab"), pos.pz());

        double PLambda = sqrt(candidate.px() * candidate.px() + candidate.py() * candidate.py() + candidate.pz() * candidate.pz());
        double ELambda = sqrt(candidate.mLambda() * candidate.mLambda() + PLambda * PLambda);
        double protonE = sqrt(massPr * massPr + pos.px() * pos.px() + pos.py() * pos.py() + pos.pz() * pos.pz());

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

        double protonsinPhiInLab = candidate.py() / sqrt(candidate.px() * candidate.px() + candidate.py() * candidate.py());
        registryData.fill(HIST("V0protonphiInLab"), protonsinPhiInLab);

        TMatrixD lambdaInJet(4, 1);
        lambdaInJet = MyTMatrixTranslationToJet(maxJetpx, maxJetpy, maxJetpz, candidate.px(), candidate.py(), candidate.pz()) * pLabV0;

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

        TMatrixD protonInV0(4, 1);
        protonInV0 = LorentzTransInV0frame(ELambda, candidate.px(), candidate.py(), candidate.pz()) * pLabproton;
        double protonMassInV0 = sqrt(protonInV0(0, 0) * protonInV0(0, 0) - protonInV0(1, 0) * protonInV0(1, 0) - protonInV0(2, 0) * protonInV0(2, 0) - protonInV0(3, 0) * protonInV0(3, 0));
        registryData.fill(HIST("V0protonMassInRest_frame"), protonMassInV0);
        registryData.fill(HIST("V0protonpxInRest_frame"), protonInV0(1, 0));
        registryData.fill(HIST("V0protonpyInRest_frame"), protonInV0(2, 0));
        registryData.fill(HIST("V0protonpzInRest_frame"), protonInV0(3, 0));
        double protonsinPhiInV0frame = protonInV0(2, 0) / sqrt(protonInV0(1, 0) * protonInV0(1, 0) + protonInV0(2, 0) * protonInV0(2, 0));
        registryData.fill(HIST("V0protonphiInRest_frame"), protonsinPhiInV0frame);

        TMatrixD protonInJetV0(4, 1);
        protonInJetV0 = LorentzTransInV0frame(ELambda, lambdaInJet(1, 0), lambdaInJet(2, 0), lambdaInJet(3, 0)) * MyTMatrixTranslationToJet(maxJetpx, maxJetpy, maxJetpz, candidate.px(), candidate.py(), candidate.pz()) * pLabproton;
        registryData.fill(HIST("V0protonpxInJetV0frame"), protonInJetV0(1, 0));
        registryData.fill(HIST("V0protonpyInJetV0frame"), protonInJetV0(2, 0));
        registryData.fill(HIST("V0protonpzInJetV0frame"), protonInJetV0(3, 0));

        double protonPinJetV0 = sqrt(protonInJetV0(1, 0) * protonInJetV0(1, 0) + protonInJetV0(2, 0) * protonInJetV0(2, 0) + protonInJetV0(3, 0) * protonInJetV0(3, 0));
        double protonPtinJetV0 = sqrt(protonInJetV0(1, 0) * protonInJetV0(1, 0) + protonInJetV0(2, 0) * protonInJetV0(2, 0));

        double protonCosThetainJetV0 = protonInJetV0(3, 0) / protonPinJetV0;
        double protonSinThetainJetV0 = protonPtinJetV0 / protonPinJetV0;
        double protonthetainJetV0 = TMath::ASin(protonSinThetainJetV0);
        registryV0Data.fill(HIST("hprotoncosthetainJetV0"), protonCosThetainJetV0);
        registryV0Data.fill(HIST("hprotonsinthetainJetV0"), protonSinThetainJetV0);
        registryV0Data.fill(HIST("hprotonthetainJetV0"), protonthetainJetV0);
        registryV0Data.fill(HIST("hprotoncosSquarethetainJetV0"), protonCosThetainJetV0 * protonCosThetainJetV0);

        registryV0Data.fill(HIST("hLambdamassandSinthetainJetV0"), candidate.mLambda(), protonSinThetainJetV0);
        registryV0Data.fill(HIST("hLambdamassandCosthetainJetV0"), candidate.mLambda(), protonCosThetainJetV0);
        registryV0Data.fill(HIST("hLambdamassandCosSquarethetainJetV0"), candidate.mLambda(), protonCosThetainJetV0 * protonCosThetainJetV0);

        registryV0Data.fill(HIST("AverageSinthetainJetV0"), candidate.mLambda(), protonSinThetainJetV0);
        registryV0Data.fill(HIST("AverageCosSquarethetainJetV0"), candidate.mLambda(), protonCosThetainJetV0 * protonCosThetainJetV0);
        protonsinPhiInJetV0frame = protonsinPhiInJetV0frame + protonInJetV0(2, 0) / sqrt(protonInJetV0(1, 0) * protonInJetV0(1, 0) + protonInJetV0(2, 0) * protonInJetV0(2, 0));

        registryData.fill(HIST("TProfile2DLambdaPtMassSinPhi"), candidate.mLambda(), candidate.pt(), protonInJetV0(2, 0) / sqrt(protonInJetV0(1, 0) * protonInJetV0(1, 0) + protonInJetV0(2, 0) * protonInJetV0(2, 0)));
        registryData.fill(HIST("TProfile2DLambdaPtMassSintheta"), candidate.mLambda(), candidate.pt(), (4.0 / TMath::Pi()) * protonSinThetainJetV0);
        registryData.fill(HIST("TProfile2DLambdaPtMassCosSquareTheta"), candidate.mLambda(), candidate.pt(), 3.0 * protonCosThetainJetV0 * protonCosThetainJetV0);
      }
      if (passedAntiLambdaSelection(candidate, pos, neg)) {
        registryData.fill(HIST("hMassAntiLambda"), candidate.mAntiLambda());
        double PAntiLambda = sqrt(candidate.px() * candidate.px() + candidate.py() * candidate.py() + candidate.pz() * candidate.pz());
        double EAntiLambda = sqrt(candidate.mAntiLambda() * candidate.mAntiLambda() + PAntiLambda * PAntiLambda);
        double AntiprotonE = sqrt(massPr * massPr + neg.px() * neg.px() + neg.py() * neg.py() + neg.pz() * neg.pz());
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
        AntiprotonsinPhiInJetV0frame = AntiprotonsinPhiInJetV0frame + AntiprotonInJetV0(2, 0) / sqrt(AntiprotonInJetV0(1, 0) * AntiprotonInJetV0(1, 0) + AntiprotonInJetV0(2, 0) * AntiprotonInJetV0(2, 0));
        TMatrixD AntiprotonInV0(4, 1);
        AntiprotonInV0 = LorentzTransInV0frame(EAntiLambda, candidate.px(), candidate.py(), candidate.pz()) * pLabAntiproton;
        double AntiprotonPinJetV0 = sqrt(AntiprotonInJetV0(1, 0) * AntiprotonInJetV0(1, 0) + AntiprotonInJetV0(2, 0) * AntiprotonInJetV0(2, 0) + AntiprotonInJetV0(3, 0) * AntiprotonInJetV0(3, 0));
        double AntiprotonPtinJetV0 = sqrt(AntiprotonInJetV0(1, 0) * AntiprotonInJetV0(1, 0) + AntiprotonInJetV0(2, 0) * AntiprotonInJetV0(2, 0));
        double AntiprotonCosThetainJetV0 = AntiprotonInJetV0(3, 0) / AntiprotonPinJetV0;
        double AntiprotonSinThetainJetV0 = AntiprotonPtinJetV0 / AntiprotonPinJetV0;
        registryData.fill(HIST("TProfile2DAntiLambdaPtMassSinPhi"), candidate.mAntiLambda(), candidate.pt(), AntiprotonInJetV0(2, 0) / sqrt(AntiprotonInJetV0(1, 0) * AntiprotonInJetV0(1, 0) + AntiprotonInJetV0(2, 0) * AntiprotonInJetV0(2, 0)));
        registryData.fill(HIST("TProfile2DAntiLambdaPtMassSintheta"), candidate.mAntiLambda(), candidate.pt(), (4.0 / TMath::Pi()) * AntiprotonSinThetainJetV0);
        registryData.fill(HIST("TProfile2DAntiLambdaPtMassCosSquareTheta"), candidate.mAntiLambda(), candidate.pt(), 3.0 * AntiprotonCosThetainJetV0 * AntiprotonCosThetainJetV0);
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
  PROCESS_SWITCH(LfMyV0s, processData, "processData", true);

  void processDataV0(SelCollisions::iterator const& collision, aod::V0Datas const& fullV0s, StrHadronDaughterTracks const&)
  {
    registryData.fill(HIST("hNEvents"), 0.5);
    if (!AcceptEvent(collision)) {
      return;
    }
    registryData.fill(HIST("hNEvents"), 8.5);
    int V0NumbersPerEvent = 0;
    int V0NumbersPerEventsel = 0;
    for (const auto& v0 : fullV0s) {
      V0NumbersPerEvent++;
      float ctauLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0;
      float ctauAntiLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0Bar;
      const auto& pos = v0.posTrack_as<StrHadronDaughterTracks>();
      const auto& neg = v0.negTrack_as<StrHadronDaughterTracks>();
      if (passedLambdaSelection(v0, pos, neg) && ctauLambda < CtauLambda && ifpasslambda) {
        V0NumbersPerEventsel++;
        registryV0Data.fill(HIST("hLambdaPt"), v0.pt());
        registryV0Data.fill(HIST("hMassVsPtLambda"), v0.pt(), v0.mLambda());
        registryV0Data.fill(HIST("hMassLambda"), v0.mLambda());
      } else if (passedInitLambdaSelection(v0, pos, neg) && ifinitpasslambda) {
        registryV0Data.fill(HIST("hLambdaPt"), v0.pt());
        registryV0Data.fill(HIST("hMassVsPtLambda"), v0.pt(), v0.mLambda());
        registryV0Data.fill(HIST("hMassLambda"), v0.mLambda());
        double PLambda = sqrt(v0.px() * v0.px() + v0.py() * v0.py() + v0.pz() * v0.pz());
        double ELambda = sqrt(v0.mLambda() * v0.mLambda() + PLambda * PLambda);
        double protonE = sqrt(massPr * massPr + pos.px() * pos.px() + pos.py() * pos.py() + pos.pz() * pos.pz());
        TMatrixD pLabproton(4, 1);
        pLabproton(0, 0) = protonE;
        pLabproton(1, 0) = pos.px();
        pLabproton(2, 0) = pos.py();
        pLabproton(3, 0) = pos.pz();
        double protonCosThetainLab = pLabproton(3, 0) / pos.p();
        double protonSinThetainLab = pos.pt() / pos.p();
        double protonthetainLab = TMath::ASin(protonSinThetainLab);
        registryV0Data.fill(HIST("hprotoncosthetainLab"), protonCosThetainLab);
        registryV0Data.fill(HIST("hprotonsinthetainLab"), protonSinThetainLab);
        registryV0Data.fill(HIST("hprotonthetainLab"), protonthetainLab);
        registryV0Data.fill(HIST("hprotoncosSquarethetainLab"), protonCosThetainLab * protonCosThetainLab);

        TMatrixD protonInV0(4, 1);
        protonInV0 = LorentzTransInV0frame(ELambda, v0.px(), v0.py(), v0.pz()) * pLabproton;
        double protonPinV0 = sqrt(protonInV0(1, 0) * protonInV0(1, 0) + protonInV0(2, 0) * protonInV0(2, 0) + protonInV0(3, 0) * protonInV0(3, 0));
        double protonPtinV0 = sqrt(protonInV0(1, 0) * protonInV0(1, 0) + protonInV0(2, 0) * protonInV0(2, 0));
        double protonCosThetainV0 = protonInV0(3, 0) / protonPinV0;
        double protonSinThetainV0 = protonPtinV0 / protonPinV0;
        double protonthetainV0 = TMath::ASin(protonSinThetainV0);
        registryV0Data.fill(HIST("hprotoncosthetainV0"), protonCosThetainV0);
        registryV0Data.fill(HIST("hprotonsinthetainV0"), protonSinThetainV0);
        registryV0Data.fill(HIST("hprotonthetainV0"), protonthetainV0);
        registryV0Data.fill(HIST("hprotoncosSquarethetainV0"), protonCosThetainV0 * protonCosThetainV0);

        registryV0Data.fill(HIST("hLambdamassandSinthetainV0"), v0.mLambda(), protonSinThetainV0);
        registryV0Data.fill(HIST("hLambdamassandCosthetainV0"), v0.mLambda(), protonCosThetainV0);
        registryV0Data.fill(HIST("hLambdamassandCosSquarethetainV0"), v0.mLambda(), protonCosThetainV0 * protonCosThetainV0);

        registryV0Data.fill(HIST("AverageSinthetainV0"), v0.mLambda(), protonSinThetainV0);
        registryV0Data.fill(HIST("AverageCosSquarethetainV0"), v0.mLambda(), protonCosThetainV0 * protonCosThetainV0);
      }
      if (passedAntiLambdaSelection(v0, pos, neg) && ctauAntiLambda < CtauLambda && ifpasslambda) {
        registryV0Data.fill(HIST("hAntiLambdaPt"), v0.pt());
        registryV0Data.fill(HIST("hMassVsPtAntiLambda"), v0.pt(), v0.mAntiLambda());
        registryV0Data.fill(HIST("hMassAntiLambda"), v0.mAntiLambda());
      } else if (passedInitLambdaSelection(v0, pos, neg) && ifinitpasslambda) {
        registryV0Data.fill(HIST("hAntiLambdaPt"), v0.pt());
        registryV0Data.fill(HIST("hMassVsPtAntiLambda"), v0.pt(), v0.mAntiLambda());
        registryV0Data.fill(HIST("hMassAntiLambda"), v0.mAntiLambda());
      }
    }
    registryV0Data.fill(HIST("nV0sPerEvent"), V0NumbersPerEvent);
    registryV0Data.fill(HIST("nV0sPerEventsel"), V0NumbersPerEventsel);
  }
  PROCESS_SWITCH(LfMyV0s, processDataV0, "processDataV0", true);

  // V0Collisions
  // SelCollisions
  using V0Collisions = soa::Join<aod::Collisions, aod::EvSels, aod::PVMults, aod::CentFT0Ms, aod::CentNGlobals>;
  void processLongitudinalPolarization(V0Collisions::iterator const& collision, aod::V0Datas const& fullV0s, StrHadronDaughterTracks const&)
  {

    if (!AcceptEventForLongitudinalPolarization(collision)) {
      return;
    }
    registryLongitudinalPolarization.fill(HIST("hNEvents"), 8.5);

    int V0NumbersPerEvent = 0;
    int V0NumbersPerEventsel = 0;
    for (const auto& v0 : fullV0s) { // loop over V0s

      if (v0.v0Type() != v0TypeSelection) {
        continue;
      }

      V0NumbersPerEvent++;
      const auto& pos = v0.posTrack_as<StrHadronDaughterTracks>();
      const auto& neg = v0.negTrack_as<StrHadronDaughterTracks>();

      if (NotITSAfterburner && (v0.negTrack_as<StrHadronDaughterTracks>().isITSAfterburner() || v0.posTrack_as<StrHadronDaughterTracks>().isITSAfterburner())) {
        continue;
      }
      if (passedInitLambdaSelection(v0, pos, neg) && ifinitpasslambda) {
        registryLongitudinalPolarization.fill(HIST("hMassVsPtLambda"), v0.pt(), v0.mLambda());
        registryLongitudinalPolarization.fill(HIST("hMassVsPtAntiLambda"), v0.pt(), v0.mAntiLambda());
        registryLongitudinalPolarization.fill(HIST("hMassLambda"), v0.mLambda());
      }

      if (AcceptV0Lambda(v0, pos, neg, collision) && ifpasslambda) {
        V0NumbersPerEventsel++;
        registryLongitudinalPolarization.fill(HIST("hMassLambdasel"), v0.mLambda());
        registryLongitudinalPolarization.fill(HIST("hMassVsPtLambdasel"), v0.pt(), v0.mLambda());

        ProtonVec = ROOT::Math::PxPyPzMVector(v0.pxpos(), v0.pypos(), v0.pzpos(), massPr);
        PionVec = ROOT::Math::PxPyPzMVector(v0.pxneg(), v0.pyneg(), v0.pzneg(), massPi);
        LambdaVec = ProtonVec + PionVec;
        LambdaVec.SetM(massLambda);
        ROOT::Math::Boost boost{LambdaVec.BoostToCM()};
        ProtonBoostedVec = boost(ProtonVec);
        LambdaBoostedVec = boost(LambdaVec);

        registryLongitudinalPolarization.fill(HIST("V0pxInRest_frame"), LambdaBoostedVec.Px());
        registryLongitudinalPolarization.fill(HIST("V0pyInRest_frame"), LambdaBoostedVec.Py());
        registryLongitudinalPolarization.fill(HIST("V0pzInRest_frame"), LambdaBoostedVec.Pz());

        double protonCosThetainV0 = ProtonBoostedVec.Pz() / ProtonBoostedVec.P();

        registryLongitudinalPolarization.fill(HIST("hprotoncosthetainV0"), protonCosThetainV0);
        registryLongitudinalPolarization.fill(HIST("hprotoncosSquarethetainV0"), protonCosThetainV0 * protonCosThetainV0);
        registryLongitudinalPolarization.fill(HIST("hLambdamassandCosthetaInV0"), v0.mLambda(), protonCosThetainV0);

        registryLongitudinalPolarization.fill(HIST("TProfile2DLambdaPtMassCostheta"), v0.mLambda(), v0.pt(), protonCosThetainV0);
        registryLongitudinalPolarization.fill(HIST("TProfile2DLambdaPtMassCosSquareTheta"), v0.mLambda(), v0.pt(), protonCosThetainV0 * protonCosThetainV0);

        registryLongitudinalPolarization.fill(HIST("TProfile1DLambdaPtMassCostheta"), v0.mLambda(), protonCosThetainV0);
      }
      if (AcceptV0AntiLambda(v0, pos, neg, collision) && ifpasslambda) {
        registryLongitudinalPolarization.fill(HIST("hMassAntiLambdasel"), v0.mAntiLambda());
        registryLongitudinalPolarization.fill(HIST("hMassVsPtAntiLambdasel"), v0.pt(), v0.mAntiLambda());

        ProtonVec = ROOT::Math::PxPyPzMVector(v0.pxneg(), v0.pyneg(), v0.pzneg(), massPr);
        PionVec = ROOT::Math::PxPyPzMVector(v0.pxpos(), v0.pypos(), v0.pzpos(), massPi);
        LambdaVec = ProtonVec + PionVec;
        LambdaVec.SetM(massLambda);
        ROOT::Math::Boost boost{LambdaVec.BoostToCM()};
        ProtonBoostedVec = boost(ProtonVec);
        LambdaBoostedVec = boost(LambdaVec);

        double protonCosThetainV0 = ProtonBoostedVec.Pz() / ProtonBoostedVec.P();

        registryLongitudinalPolarization.fill(HIST("hantiprotoncosthetainV0"), protonCosThetainV0);
        registryLongitudinalPolarization.fill(HIST("hantiprotoncosSquarethetainV0"), protonCosThetainV0 * protonCosThetainV0);
        registryLongitudinalPolarization.fill(HIST("hAntiLambdamassandCosthetaInV0"), v0.mAntiLambda(), protonCosThetainV0);

        registryLongitudinalPolarization.fill(HIST("TProfile2DAntiLambdaPtMassCostheta"), v0.mAntiLambda(), v0.pt(), protonCosThetainV0);
        registryLongitudinalPolarization.fill(HIST("TProfile2DAntiLambdaPtMassCosSquareTheta"), v0.mAntiLambda(), v0.pt(), protonCosThetainV0 * protonCosThetainV0);
        registryLongitudinalPolarization.fill(HIST("TProfile1DAntiLambdaPtMassCostheta"), v0.mAntiLambda(), protonCosThetainV0);
      }
    }
    registryLongitudinalPolarization.fill(HIST("nV0sPerEvent"), V0NumbersPerEvent);
    registryLongitudinalPolarization.fill(HIST("nV0sPerEventsel"), V0NumbersPerEventsel);
  }
  PROCESS_SWITCH(LfMyV0s, processLongitudinalPolarization, "processLongitudinalPolarization", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<LfMyV0s>(cfgc, TaskName{"lf-my-v0s"}),
  };
}
