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

/// \file higherMassResonances.cxx
/// \brief glueball resonance
/// \author Sawan <sawan.sawan@cern.ch>

// #include <TDatabasePDG.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THn.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TPDGCode.h>
#include <array>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include "TF1.h"
#include "TRandom3.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Boost.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/StepTHn.h"
#include "ReconstructionDataFormats/Track.h"
#include "Framework/O2DatabasePDGPlugin.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h" //
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h" //
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"              //
#include "Framework/runDataProcessing.h"         //
#include "PWGLF/DataModel/LFStrangenessTables.h" //

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
// using namespace o2::constants::physics;
using std::array;

struct HigherMassResonances {
  SliceCache cache;
  HistogramRegistry rEventSelection{"eventSelection", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rKzeroShort{"kzeroShort", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry hglue{"hglueball", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry hMChists{"hMChists", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  struct : ConfigurableGroup {
    // PID and QA
    Configurable<bool> qAv0{"qAv0", false, "qAv0"};
    Configurable<bool> qAPID{"qAPID", true, "qAPID"};
    // Configurable<bool> qAv0Daughters{"qAv0Daughters", false, "QA of v0 daughters"};
    Configurable<bool> qAevents{"qAevents", false, "QA of events"};
    // Configurable<bool> invMass1D{"invMass1D", false, "1D invariant mass histograms"};
    Configurable<bool> correlation2Dhist{"correlation2Dhist", true, "Lamda K0 mass correlation"};
    Configurable<bool> cDCAv0topv{"cDCAv0topv", false, "DCA V0 to PV"};
    Configurable<bool> armcut{"armcut", false, "arm cut"};
    Configurable<bool> globalTracks{"globalTracks", false, "Global tracks"};
    Configurable<bool> hasTPC{"hasTPC", false, "TPC"};
    Configurable<bool> selectTWOKsOnly{"selectTWOKsOnly", true, "Select only events with two K0s"};

    // Configurables for event selection
    Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range (cm)"};
    Configurable<float> cfgETAcut{"cfgETAcut", 0.8f, "Track ETA cut"};
    Configurable<bool> timFrameEvsel{"timFrameEvsel", false, "TPC Time frame boundary cut"};
    Configurable<bool> piluprejection{"piluprejection", false, "Pileup rejection"};
    Configurable<bool> goodzvertex{"goodzvertex", false, "removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference."};
    Configurable<bool> itstpctracks{"itstpctracks", false, "selects collisions with at least one ITS-TPC track,"};
    Configurable<bool> additionalEvsel{"additionalEvsel", false, "Additional event selcection"};
    // Configurable<bool> applyOccupancyCut{"applyOccupancyCut", false, "Apply occupancy cut"};
    // Configurable<int> occupancyCut{"occupancyCut", 1000, "Mimimum Occupancy cut"};

    // Configurable parameters for V0 selection
    Configurable<float> confV0DCADaughMax{"confV0DCADaughMax", 1.0f, "DCA b/w V0 daughters"};
    Configurable<float> v0settingDcapostopv{"v0settingDcapostopv", 0.06, "DCA Pos To PV"};
    Configurable<float> v0settingDcanegtopv{"v0settingDcanegtopv", 0.06, "DCA Neg To PV"};
    Configurable<double> cMaxV0DCA{"cMaxV0DCA", 1, "DCA V0 to PV"};
    // Configurable<bool> isStandarv0{"isStandarv0", false, "Standard V0"};
    // Configurable<float> ConfDaughDCAMin{"ConfDaughDCAMin", 0.06f, "V0 Daugh sel:  Max. DCA Daugh to PV (cm)"}; // same as DCA pos to pv and neg to pv
    Configurable<float> confV0PtMin{"confV0PtMin", 0.f, "Minimum transverse momentum of V0"};
    Configurable<float> confV0CPAMin{"confV0CPAMin", 0.97f, "Minimum CPA of V0"};
    Configurable<float> confV0TranRadV0Min{"confV0TranRadV0Min", 0.5f, "Minimum transverse radius"};
    Configurable<float> confV0TranRadV0Max{"confV0TranRadV0Max", 200.f, "Maximum transverse radius"};
    Configurable<double> cMaxV0LifeTime{"cMaxV0LifeTime", 15, "Maximum V0 life time"};
    Configurable<double> cSigmaMassKs0{"cSigmaMassKs0", 4, "n Sigma cut on Ks0 mass (Mass (Ks) - cSigmaMassKs0*cWidthKs0)"};
    Configurable<double> cWidthKs0{"cWidthKs0", 0.005, "Width of KS0"};
    Configurable<float> confDaughEta{"confDaughEta", 0.8f, "V0 Daugh sel: max eta"};
    Configurable<float> confDaughTPCnclsMin{"confDaughTPCnclsMin", 70.f, "V0 Daugh sel: Min. nCls TPC"};
    Configurable<float> confDaughPIDCuts{"confDaughPIDCuts", 5, "PID selections for KS0 daughters"};
    Configurable<float> confarmcut{"confarmcut", 0.2f, "Armenteros cut"};
    Configurable<float> confKsrapidity{"confKsrapidity", 0.5f, "Rapidity cut on K0s"};
    // Configurable<float> lowmasscutks0{"lowmasscutks0", 0.497 - 4 * 0.005, "Low mass cut on K0s"};
    // Configurable<float> highmasscutks0{"highmasscutks0", 0.497 + 4 * 0.005, "High mass cut on K0s"};
    Configurable<bool> applyAngSepCut{"applyAngSepCut", false, "Apply angular separation cut"};
    Configurable<float> angSepCut{"angSepCut", 0.01f, "Angular separation cut"};

    // Configurable for track selection and multiplicity
    Configurable<float> cfgPTcut{"cfgPTcut", 0.2f, "Track PT cut"};
    Configurable<int> cfgNmixedEvents{"cfgNmixedEvents", 5, "Number of mixed events"};
    Configurable<bool> cfgMultFOTM{"cfgMultFOTM", true, "Use FOTM multiplicity if pp else use 0 here for PbPb (FT0C)"};
    ConfigurableAxis binsCent{"binsCent", {VARIABLE_WIDTH, 0., 5., 10., 30., 50., 70., 100., 110., 150.}, "Binning of the centrality axis"};

    // Configurable for MC
    Configurable<bool> allGenCollisions{"allGenCollisions", true, "To fill all generated collisions for the signal loss calculations"};
    Configurable<bool> cTVXEvsel{"cTVXEvsel", true, "Triggger selection"};
    Configurable<bool> avoidsplitrackMC{"avoidsplitrackMC", false, "avoid split track in MC"};
    Configurable<int> selectMCparticles{"selectMCparticles", 1, "0: f0(1710), 1: f2(1525), 2: a2(1320), 3: f0(1370), 4: f0(1500)"};
    std::vector<int> pdgCodes = {10331, 335, 115, 10221, 9030221};

    // output THnSparses
    Configurable<bool> activateTHnSparseCosThStarHelicity{"activateTHnSparseCosThStarHelicity", false, "Activate the THnSparse with cosThStar w.r.t. helicity axis"};
    Configurable<bool> activateTHnSparseCosThStarProduction{"activateTHnSparseCosThStarProduction", false, "Activate the THnSparse with cosThStar w.r.t. production axis"};
    Configurable<bool> activateTHnSparseCosThStarBeam{"activateTHnSparseCosThStarBeam", true, "Activate the THnSparse with cosThStar w.r.t. beam axis (Gottified jackson frame)"};
    Configurable<bool> activateTHnSparseCosThStarRandom{"activateTHnSparseCosThStarRandom", false, "Activate the THnSparse with cosThStar w.r.t. random axis"};
    Configurable<int> cRotations{"cRotations", 3, "Number of random rotations in the rotational background"};

    // Other cuts on Ks and glueball
    // Configurable<bool> rapidityks{"rapidityks", true, "rapidity cut on K0s"};
    Configurable<bool> applyCompetingcut{"applyCompetingcut", false, "Competing cascade rejection cut"};
    Configurable<float> competingcascrejlambda{"competingcascrejlambda", 0.005, "rejecting competing cascade lambda"};
    Configurable<float> competingcascrejlambdaanti{"competingcascrejlambdaanti", 0.005, "rejecting competing cascade anti-lambda"}; // If one of the pions is misidentified as a proton, then instead of Ks we reconstruct lambda, therefore the competing cascade rejection cut is applied in which if the reconstrcted mass of a pion and proton (which we are assuming to be misidentified as proton) is close to lambda or anti-lambda, then the track is rejected
    Configurable<int> tpcCrossedrows{"tpcCrossedrows", 70, "TPC crossed rows"};
    Configurable<float> tpcCrossedrowsOverfcls{"tpcCrossedrowsOverfcls", 0.8, "TPC crossed rows over findable clusters"};

    // Mass and pT axis as configurables
    Configurable<float> cPtMin{"cPtMin", 0.0f, "Minimum pT"};
    Configurable<float> cPtMax{"cPtMax", 30.0f, "Maximum pT"};
    Configurable<int> cPtBins{"cPtBins", 300, "Number of pT bins"};
    Configurable<float> cMassMin{"cMassMin", 0.9f, "Minimum mass of glueball"};
    Configurable<float> cMassMax{"cMassMax", 3.0f, "Maximum mass of glueball"};
    Configurable<int> cMassBins{"cMassBins", 210, "Number of mass bins for glueball"};
    Configurable<float> ksMassMin{"ksMassMin", 0.45f, "Minimum mass of K0s"};
    Configurable<float> ksMassMax{"ksMassMax", 0.55f, "Maximum mass of K0s"};
    Configurable<int> ksMassBins{"ksMassBins", 200, "Number of mass bins for K0s"};
    Configurable<int> rotationalCut{"rotationalCut", 10, "Cut value (Rotation angle pi - pi/cut and pi + pi/cut)"};
    ConfigurableAxis configThnAxisPOL{"configThnAxisPOL", {20, -1.0, 1.0}, "Costheta axis"};
    // ConfigurableAxis axisdEdx{"axisdEdx", {20000, 0.0f, 200.0f}, "dE/dx (a.u.)"};
    // ConfigurableAxis axisPtfordEbydx{"axisPtfordEbydx", {2000, 0, 20}, "pT (GeV/c)"};
    // ConfigurableAxis axisMultdist{"axisMultdist", {3500, 0, 70000}, "Multiplicity distribution"};
    // ConfigurableAxis occupancyBins{"occupancyBins", {VARIABLE_WIDTH, 0.0, 100, 500, 600, 1000, 1100, 1500, 1600, 2000, 2100, 2500, 2600, 3000, 3100, 3500, 3600, 4000, 4100, 4500, 4600, 5000, 5100, 9999}, "Binning: occupancy axis"};
  } config;

  // Service<o2::framework::O2DatabasePDG> PDGdatabase;
  TRandom* rn = new TRandom();

  // variables declaration
  TLorentzVector lv1, lv2, lv3, lv4, lv5;
  float multiplicity = 0.0f;
  float theta2;
  ROOT::Math::PxPyPzMVector daughter1, daughter2, fourVecDau1, fourVecMother, fourVecDauCM;
  ROOT::Math::XYZVector threeVecDauCM, helicityVec, randomVec, beamVec, normalVec;
  // const double massK0s = o2::constants::physics::MassK0Short;
  bool isMix = false;

  void init(InitContext const&)
  {
    // Axes
    AxisSpec k0ShortMassAxis = {config.ksMassBins, config.ksMassMin, config.ksMassMax, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec glueballMassAxis = {config.cMassBins, config.cMassMin, config.cMassMax, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec vertexZAxis = {60, -15.f, 15.f, "vrtx_{Z} [cm]"}; // for histogram
    AxisSpec ptAxis = {config.cPtBins, config.cPtMin, config.cPtMax, "#it{p}_{T} (GeV/#it{c})"};
    // AxisSpec multiplicityAxis = {110, 0.0f, 150.0f, "Multiplicity Axis"};
    AxisSpec multiplicityAxis = {config.binsCent, "Multiplicity Axis"};
    AxisSpec thnAxisPOL{config.configThnAxisPOL, "Configurabel theta axis"};
    // AxisSpec occupancyAxis = {occupancyBins, "Occupancy [-40,100]"};

    //  THnSparses
    std::array<int, 4> sparses = {config.activateTHnSparseCosThStarHelicity, config.activateTHnSparseCosThStarProduction, config.activateTHnSparseCosThStarBeam, config.activateTHnSparseCosThStarRandom};

    // std::array<int, 1> sparses = {config.activateTHnSparseCosThStarHelicity};

    if (std::accumulate(sparses.begin(), sparses.end(), 0) == 0) {
      LOGP(fatal, "No output THnSparses enabled");
    } else {
      if (config.activateTHnSparseCosThStarHelicity) {
        LOGP(info, "THnSparse with cosThStar w.r.t. helicity axis active.");
      }
      if (config.activateTHnSparseCosThStarProduction) {
        LOGP(info, "THnSparse with cosThStar w.r.t. production axis active.");
      }
      if (config.activateTHnSparseCosThStarBeam) {
        LOGP(info, "THnSparse with cosThStar w.r.t. beam axis active. (Gottified jackson frame)");
      }
      if (config.activateTHnSparseCosThStarRandom) {
        LOGP(info, "THnSparse with cosThStar w.r.t. random axis active.");
      }
    }

    // Event selection
    if (config.qAevents) {
      rEventSelection.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});
      rEventSelection.add("hmultiplicity", "multiplicity percentile distribution", {HistType::kTH1F, {{150, 0.0f, 150.0f}}});
      // rEventSelection.add("multdist_FT0M", "FT0M Multiplicity distribution", kTH1F, {config.axisMultdist});
      // rEventSelection.add("multdist_FT0A", "FT0A Multiplicity distribution", kTH1F, {config.axisMultdist});
      // rEventSelection.add("multdist_FT0C", "FT0C Multiplicity distribution", kTH1F, {config.axisMultdist});
      // rEventSelection.add("hNcontributor", "Number of primary vertex contributor", kTH1F, {{2000, 0.0f, 10000.0f}});
    }

    hglue.add("h3glueInvMassDS", "h3glueInvMassDS", kTHnSparseF, {multiplicityAxis, ptAxis, glueballMassAxis, thnAxisPOL}, true);
    hglue.add("h3glueInvMassME", "h3glueInvMassME", kTHnSparseF, {multiplicityAxis, ptAxis, glueballMassAxis, thnAxisPOL}, true);
    hglue.add("h3glueInvMassRot", "h3glueInvMassRot", kTHnSparseF, {multiplicityAxis, ptAxis, glueballMassAxis, thnAxisPOL}, true);
    hglue.add("heventscheck", "heventscheck", kTH1I, {{10, 0, 10}});
    hglue.add("htrackscheck_v0", "htrackscheck_v0", kTH1I, {{15, 0, 15}});
    hglue.add("htrackscheck_v0_daughters", "htrackscheck_v0_daughters", kTH1I, {{15, 0, 15}});

    // K0s topological/PID cuts
    if (config.correlation2Dhist) {
      rKzeroShort.add("mass_lambda_kshort_before", "mass under lambda hypotheses and Kshort mass", kTH2F, {{100, 0.2, 0.8}, {100, 0.9, 1.5}});
      // rKzeroShort.add("mass_lambda_kshort_after1", "mass under lambda hypotheses and Kshort mass", kTH2F, {{100, 0.2, 0.8}, {100, 0.9, 1.5}});
      // rKzeroShort.add("mass_lambda_kshort_after2", "mass under lambda hypotheses and Kshort mass", kTH2F, {{100, 0.2, 0.8}, {100, 0.9, 1.5}});
      // rKzeroShort.add("mass_lambda_kshort_after3", "mass under lambda hypotheses and Kshort mass", kTH2F, {{100, 0.2, 0.8}, {100, 0.9, 1.5}});
      // rKzeroShort.add("mass_lambda_kshort_after4", "mass under lambda hypotheses and Kshort mass", kTH2F, {{100, 0.2, 0.8}, {100, 0.9, 1.5}});
      // rKzeroShort.add("mass_lambda_kshort_after5", "mass under lambda hypotheses and Kshort mass", kTH2F, {{100, 0.2, 0.8}, {100, 0.9, 1.5}});
      // rKzeroShort.add("mass_lambda_kshort_after6", "mass under lambda hypotheses and Kshort mass", kTH2F, {{100, 0.2, 0.8}, {100, 0.9, 1.5}});
      // rKzeroShort.add("mass_lambda_kshort_after7", "mass under lambda hypotheses and Kshort mass", kTH2F, {{100, 0.2, 0.8}, {100, 0.9, 1.5}});
      // rKzeroShort.add("mass_lambda_kshort_after8", "mass under lambda hypotheses and Kshort mass", kTH2F, {{100, 0.2, 0.8}, {100, 0.9, 1.5}});
      // rKzeroShort.add("mass_lambda_kshort_after9", "mass under lambda hypotheses and Kshort mass", kTH2F, {{100, 0.2, 0.8}, {100, 0.9, 1.5}});
      rKzeroShort.add("mass_lambda_kshort_after10", "mass under lambda hypotheses and Kshort mass", kTH2F, {{100, 0.2, 0.8}, {100, 0.9, 1.5}});
    }
    if (config.qAv0) {
      // Invariant Mass
      rKzeroShort.add("hMassK0Shortbefore", "hMassK0Shortbefore", kTHnSparseF, {k0ShortMassAxis, ptAxis});
      rKzeroShort.add("hMasscorrelationbefore", "hMasscorrelationbefore", kTH2F, {k0ShortMassAxis, k0ShortMassAxis});
      rKzeroShort.add("hMassK0ShortSelected", "hMassK0ShortSelected", kTHnSparseF, {k0ShortMassAxis, ptAxis});
      // Topological histograms (after the selection)
      rKzeroShort.add("hDCAV0Daughters", "DCA between v0 daughters", {HistType::kTH1F, {{60, -3.0f, 3.0f}}});
      rKzeroShort.add("hV0CosPA", "hV0CosPA", {HistType::kTH1F, {{100, 0.96f, 1.1f}}});
      rKzeroShort.add("hLT", "hLT", {HistType::kTH1F, {{100, 0.0f, 50.0f}}});
      rKzeroShort.add("angularSeparation", "Angular distribution between two K0s vs pT", {HistType::kTH1F, {{200, 0.0f, 4.0f}}});
      // rKzeroShort.add("Mass_lambda", "Mass under lambda hypothesis", kTH1F, {glueballMassAxis});
      // rKzeroShort.add("mass_AntiLambda", "Mass under anti-lambda hypothesis", kTH1F, {glueballMassAxis});
      // rKzeroShort.add("mass_Gamma", "Mass under Gamma hypothesis", kTH1F, {glueballMassAxis});

      // rKzeroShort.add("mass_Hypertriton", "Mass under hypertriton hypothesis", kTH1F, {glueballMassAxis});
      // rKzeroShort.add("mass_AnitHypertriton", "Mass under anti-hypertriton hypothesis", kTH1F, {glueballMassAxis});
      // rKzeroShort.add("rapidity", "Rapidity distribution", kTH1F, {{100, -1.0f, 1.0f}});
      // rKzeroShort.add("hv0radius", "hv0radius", kTH1F, {{100, 0.0f, 200.0f}});
      // rKzeroShort.add("hDCApostopv", "DCA positive daughter to PV", kTH1F, {{1000, -10.0f, 10.0f}});
      // rKzeroShort.add("hDCAnegtopv", "DCA negative daughter to PV", kTH1F, {{1000, -10.0f, 10.0f}});
      // rKzeroShort.add("hcDCAv0topv", "DCA V0 to PV", kTH1F, {{60, -3.0f, 3.0f}});
      // rKzeroShort.add("halpha", "Armenteros alpha", kTH1F, {{100, -5.0f, 5.0f}});
      // rKzeroShort.add("hqtarmbyalpha", "qtarm/alpha", kTH1F, {{100, 0.0f, 1.0f}});
      // rKzeroShort.add("hpsipair", "psi pair angle", kTH1F, {{100, -5.0f, 5.0f}});

      // // Topological histograms (before the selection)
      // rKzeroShort.add("hDCAV0Daughters_before", "DCA between v0 daughters before the selection", {HistType::kTH1F, {{60, -3.0f, 3.0f}}});
      // rKzeroShort.add("hV0CosPA_before", "hV0CosPA_before", {HistType::kTH1F, {{200, 0.91f, 1.1f}}});
      // rKzeroShort.add("hLT_before", "hLT_before", {HistType::kTH1F, {{100, 0.0f, 50.0f}}});
    }
    rKzeroShort.add("NksProduced", "Number of K0s produced", kTH1I, {{15, -0.5, 14.5}});

    if (config.qAPID) {
      rKzeroShort.add("hNSigmaPosPionK0s_before", "hNSigmaPosPionK0s_before", {HistType::kTH2F, {{ptAxis}, {100, -5.f, 5.f}}});
      rKzeroShort.add("hNSigmaPosPionK0s_after", "hNSigmaPosPionK0s_after", {HistType::kTH2F, {{ptAxis}, {100, -5.f, 5.f}}});
      rKzeroShort.add("hNSigmaNegPionK0s_before", "hNSigmaNegPionK0s_before", {HistType::kTH2F, {{ptAxis}, {100, -5.f, 5.f}}});
      rKzeroShort.add("hNSigmaNegPionK0s_after", "hNSigmaNegPionK0s_after", {HistType::kTH2F, {{ptAxis}, {100, -5.f, 5.f}}});
      // rKzeroShort.add("dE_by_dx_TPC", "dE/dx signal in the TPC as a function of pT", kTH2F, {config.axisPtfordEbydx, config.axisdEdx});
    }
    // if (config.qAv0Daughters) {
    //   rKzeroShort.add("negative_pt", "Negative daughter pT", kTH1F, {ptAxis});
    //   rKzeroShort.add("positive_pt", "Positive daughter pT", kTH1F, {ptAxis});
    //   rKzeroShort.add("negative_eta", "Negative daughter eta", kTH1F, {{100, -1.0f, 1.0f}});
    //   rKzeroShort.add("positive_eta", "Positive daughter eta", kTH1F, {{100, -1.0f, 1.0f}});
    //   rKzeroShort.add("negative_phi", "Negative daughter phi", kTH1F, {{70, 0.0f, 7.0f}});
    //   rKzeroShort.add("positive_phi", "Positive daughter phi", kTH1F, {{70, 0.0f, 7.0f}});
    // }

    // For MC
    hMChists.add("events_check", "No. of events in the generated MC", kTH1I, {{20, 0, 20}});
    hMChists.add("events_checkrec", "No. of events in the reconstructed MC", kTH1I, {{20, 0, 20}});
    hMChists.add("Genf1710", "Gen f_{0}(1710)", kTHnSparseF, {multiplicityAxis, ptAxis});
    hMChists.add("Recf1710_pt1", "Rec f_{0}(1710) p_{T}", kTHnSparseF, {multiplicityAxis, ptAxis, glueballMassAxis});
    hMChists.add("Recf1710_pt2", "Rec f_{0}(1710) p_{T}", kTHnSparseF, {multiplicityAxis, ptAxis, glueballMassAxis});
    hMChists.add("Recf1710_p", "Rec f_{0}(1710) p", kTH1F, {ptAxis});
    hMChists.add("h1Recsplit", "Rec p_{T}2", kTH1F, {ptAxis});

    hMChists.add("Recf1710_mass", "Rec f_{0}(1710) mass", kTH1F, {glueballMassAxis});
    hMChists.add("Genf1710_mass", "Gen f_{0}(1710) mass", kTH1F, {glueballMassAxis});
    hMChists.add("GenEta", "Gen Eta", kTHnSparseF, {ptAxis, {100, -1.0f, 1.0f}});
    hMChists.add("GenPhi", "Gen Phi", kTH1F, {{70, -3.5f, 3.5f}});
    hMChists.add("GenRapidity", "Gen Rapidity", kTHnSparseF, {ptAxis, {100, -1.0f, 1.0f}});
    hMChists.add("RecEta", "Rec Eta", kTH1F, {{100, -1.0f, 1.0f}});
    hMChists.add("RecPhi", "Rec Phi", kTH1F, {{70, 0.0f, 7.0f}});
    hMChists.add("RecRapidity", "Rec Rapidity", kTH1F, {{100, -1.0f, 1.0f}});
    hMChists.add("MC_mult", "Multiplicity in MC", kTH1F, {multiplicityAxis});
    hMChists.add("MC_mult_after_event_sel", "Multiplicity in MC", kTH1F, {multiplicityAxis});
    // hMChists.add("GenPx", "Gen Px", kTH1F, {{100, -10.0f, 10.0f}});
    // hMChists.add("GenPy", "Gen Py", kTH1F, {{100, -10.0f, 10.0f}});
    // hMChists.add("GenPz", "Gen Pz", kTH1F, {{100, -10.0f, 10.0f}});
  }

  template <typename Collision>
  bool eventselection(Collision const& collision)
  {
    hglue.fill(HIST("heventscheck"), 1.5);

    if (config.timFrameEvsel && (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
      return false;
    }
    hglue.fill(HIST("heventscheck"), 2.5);

    if (!collision.sel8()) {
      return false;
    }
    hglue.fill(HIST("heventscheck"), 3.5);

    if (config.piluprejection && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    hglue.fill(HIST("heventscheck"), 4.5);

    if (config.goodzvertex && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    hglue.fill(HIST("heventscheck"), 5.5);

    if (config.itstpctracks && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return false;
    }
    hglue.fill(HIST("heventscheck"), 6.5);

    return true;
  }

  template <typename Collision, typename V0>
  bool selectionV0(Collision const& collision, V0 const& candidate, float /*multiplicity*/)
  {
    const float qtarm = candidate.qtarm();
    const float alph = candidate.alpha();
    float arm = qtarm / alph;
    const float pT = candidate.pt();
    const float tranRad = candidate.v0radius();
    const float dcaDaughv0 = candidate.dcaV0daughters();
    const float cpav0 = candidate.v0cosPA();

    float ctauK0s = candidate.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short;
    float lowmasscutks0 = 0.497 - config.cWidthKs0 * config.cSigmaMassKs0;
    float highmasscutks0 = 0.497 + config.cWidthKs0 * config.cSigmaMassKs0;
    // float decayLength = candidate.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * RecoDecay::sqrtSumOfSquares(candidate.px(), candidate.py(), candidate.pz());

    if (config.qAv0) {
      rKzeroShort.fill(HIST("hMassK0Shortbefore"), candidate.mK0Short(), candidate.pt());
      rKzeroShort.fill(HIST("hLT"), ctauK0s);
      rKzeroShort.fill(HIST("hDCAV0Daughters"), candidate.dcaV0daughters());
      rKzeroShort.fill(HIST("hV0CosPA"), candidate.v0cosPA());
      // rKzeroShort.fill(HIST("Mass_lambda"), candidate.mLambda());
      // rKzeroShort.fill(HIST("mass_AntiLambda"), candidate.mAntiLambda());
      // rKzeroShort.fill(HIST("mass_Gamma"), candidate.mGamma());
      // rKzeroShort.fill(HIST("mass_Hypertriton"), candidate.mHypertriton());
      // rKzeroShort.fill(HIST("mass_AnitHypertriton"), candidate.mAntiHypertriton());
      // rKzeroShort.fill(HIST("rapidity"), candidate.yK0Short());
      // rKzeroShort.fill(HIST("hv0radius"), candidate.v0radius());
      // rKzeroShort.fill(HIST("hDCApostopv"), candidate.dcapostopv());
      // rKzeroShort.fill(HIST("hDCAnegtopv"), candidate.dcanegtopv());
      // rKzeroShort.fill(HIST("hcDCAv0topv"), candidate.dcav0topv());
      // rKzeroShort.fill(HIST("halpha"), candidate.alpha());
      // rKzeroShort.fill(HIST("hqtarmbyalpha"), arm);
      // rKzeroShort.fill(HIST("hpsipair"), candidate.psipair());
    }
    if (config.correlation2Dhist)
      rKzeroShort.fill(HIST("mass_lambda_kshort_before"), candidate.mK0Short(), candidate.mLambda());

    hglue.fill(HIST("htrackscheck_v0"), 0.5);

    if (config.cDCAv0topv && std::fabs(candidate.dcav0topv()) > config.cMaxV0DCA) {
      return false;
    }
    hglue.fill(HIST("htrackscheck_v0"), 1.5);
    // if (config.correlation2Dhist)
    //   rKzeroShort.fill(HIST("mass_lambda_kshort_after1"), candidate.mK0Short(), candidate.mLambda());

    // if (config.rapidityks && std::abs(candidate.yK0Short()) >= config.confKsrapidity) {
    //   return false;
    // }
    if (std::abs(candidate.yK0Short()) >= config.confKsrapidity) {
      return false;
    }
    hglue.fill(HIST("htrackscheck_v0"), 2.5);
    // if (config.correlation2Dhist)
    //   rKzeroShort.fill(HIST("mass_lambda_kshort_after2"), candidate.mK0Short(), candidate.mLambda());

    if (pT < config.confV0PtMin) {
      return false;
    }
    hglue.fill(HIST("htrackscheck_v0"), 3.5);
    // if (config.correlation2Dhist)
    //   rKzeroShort.fill(HIST("mass_lambda_kshort_after3"), candidate.mK0Short(), candidate.mLambda());

    if (dcaDaughv0 > config.confV0DCADaughMax) {
      return false;
    }
    hglue.fill(HIST("htrackscheck_v0"), 4.5);
    // if (config.correlation2Dhist)
    //   rKzeroShort.fill(HIST("mass_lambda_kshort_after4"), candidate.mK0Short(), candidate.mLambda());

    if (cpav0 < config.confV0CPAMin) {
      return false;
    }
    hglue.fill(HIST("htrackscheck_v0"), 5.5);
    // if (config.correlation2Dhist)
    //   rKzeroShort.fill(HIST("mass_lambda_kshort_after5"), candidate.mK0Short(), candidate.mLambda());

    if (tranRad < config.confV0TranRadV0Min) {
      return false;
    }
    hglue.fill(HIST("htrackscheck_v0"), 6.5);
    // if (config.correlation2Dhist)
    //   rKzeroShort.fill(HIST("mass_lambda_kshort_after6"), candidate.mK0Short(), candidate.mLambda());

    if (tranRad > config.confV0TranRadV0Max) {
      return false;
    }
    hglue.fill(HIST("htrackscheck_v0"), 7.5);
    // if (config.correlation2Dhist)
    //   rKzeroShort.fill(HIST("mass_lambda_kshort_after7"), candidate.mK0Short(), candidate.mLambda());

    if (std::fabs(ctauK0s) > config.cMaxV0LifeTime) {
      return false;
    }
    hglue.fill(HIST("htrackscheck_v0"), 8.5);
    // if (config.correlation2Dhist)
    //   rKzeroShort.fill(HIST("mass_lambda_kshort_after8"), candidate.mK0Short(), candidate.mLambda());

    if (config.armcut && arm < config.confarmcut) {
      return false;
    }
    hglue.fill(HIST("htrackscheck_v0"), 9.5);
    // if (config.correlation2Dhist)
    //   rKzeroShort.fill(HIST("mass_lambda_kshort_after9"), candidate.mK0Short(), candidate.mLambda());

    // if (config.applyCompetingcut && (std::abs(candidate.mLambda() - PDGdatabase->Mass(3122)) <= config.competingcascrejlambda || std::abs(candidate.mAntiLambda() - PDGdatabase->Mass(-3122)) <= config.competingcascrejlambdaanti))
    if (config.applyCompetingcut && (std::abs(candidate.mLambda() - o2::constants::physics::MassLambda0) <= config.competingcascrejlambda || std::abs(candidate.mAntiLambda() - o2::constants::physics::MassLambda0) <= config.competingcascrejlambdaanti)) {
      return false;
    }
    hglue.fill(HIST("htrackscheck_v0"), 10.5);
    if (config.correlation2Dhist)
      rKzeroShort.fill(HIST("mass_lambda_kshort_after10"), candidate.mK0Short(), candidate.mLambda());

    if (config.qAv0) {
      rKzeroShort.fill(HIST("hMassK0ShortSelected"), candidate.mK0Short(), candidate.pt());
      // rKzeroShort.fill(HIST("mass_lambda_kshort_after"), candidate.mK0Short(), candidate.mLambda());
    }

    if (candidate.mK0Short() < lowmasscutks0 || candidate.mK0Short() > highmasscutks0) {
      return false;
    }
    return true;
  }

  template <typename T, typename V0s>
  bool isSelectedV0Daughter(T const& track, float charge, double nsigmaV0Daughter, V0s const& /*candidate*/)
  {
    if (config.qAPID) {
      // Filling the PID of the V0 daughters in the region of the K0 peak.
      (charge == 1) ? rKzeroShort.fill(HIST("hNSigmaPosPionK0s_before"), track.tpcInnerParam(), track.tpcNSigmaPi()) : rKzeroShort.fill(HIST("hNSigmaNegPionK0s_before"), track.tpcInnerParam(), track.tpcNSigmaPi());
      // rKzeroShort.fill(HIST("dE_by_dx_TPC"), track.p(), track.tpcSignal());
    }
    const auto eta = track.eta();
    const auto tpcNClsF = track.tpcNClsFound();
    const auto sign = track.sign();

    hglue.fill(HIST("htrackscheck_v0_daughters"), 0.5);

    if (config.hasTPC && !track.hasTPC())
      return false;
    hglue.fill(HIST("htrackscheck_v0_daughters"), 1.5);

    if (!config.globalTracks) {
      if (track.tpcNClsCrossedRows() < config.tpcCrossedrows)
        return false;
      hglue.fill(HIST("htrackscheck_v0_daughters"), 2.5);

      if (track.tpcCrossedRowsOverFindableCls() < config.tpcCrossedrowsOverfcls)
        return false;
      hglue.fill(HIST("htrackscheck_v0_daughters"), 3.5);

      if (tpcNClsF < config.confDaughTPCnclsMin) {
        return false;
      }
      hglue.fill(HIST("htrackscheck_v0_daughters"), 4.5);
    } else {
      if (!track.isGlobalTrack())
        return false;
      hglue.fill(HIST("htrackscheck_v0_daughters"), 4.5);
    }

    if (charge < 0 && sign > 0) {
      return false;
    }
    hglue.fill(HIST("htrackscheck_v0_daughters"), 5.5);

    if (charge > 0 && sign < 0) {
      return false;
    }
    hglue.fill(HIST("htrackscheck_v0_daughters"), 6.5);

    if (std::abs(eta) > config.confDaughEta) {
      return false;
    }
    hglue.fill(HIST("htrackscheck_v0_daughters"), 7.5);

    if (std::abs(nsigmaV0Daughter) > config.confDaughPIDCuts) {
      return false;
    }
    hglue.fill(HIST("htrackscheck_v0_daughters"), 8.5);

    if (config.qAPID) {
      (charge == 1) ? rKzeroShort.fill(HIST("hNSigmaPosPionK0s_after"), track.tpcInnerParam(), track.tpcNSigmaPi()) : rKzeroShort.fill(HIST("hNSigmaNegPionK0s_after"), track.tpcInnerParam(), track.tpcNSigmaPi());
    }

    return true;
  }

  // Angular separation cut on KsKs pairs
  template <typename T1, typename T2>
  bool applyAngSep(const T1& candidate1, const T2& candidate2)
  {
    double eta1, eta2, phi1, phi2;
    eta1 = candidate1.eta();
    eta2 = candidate2.eta();
    phi1 = candidate1.phi();
    phi2 = candidate2.phi();

    double angle = std::sqrt(std::pow(eta1 - eta2, 2) + std::pow(phi1 - phi2, 2));
    rKzeroShort.fill(HIST("angularSeparation"), angle);
    if (config.applyAngSepCut && angle > config.angSepCut) {
      return false;
    }
    return true;
  }

  // Defining filters for events (event selection)
  // Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < config.cutzvertex);
  Filter acceptenceFilter = (nabs(aod::track::eta) < config.cfgETAcut && nabs(aod::track::pt) > config.cfgPTcut);

  // Filters on V0s
  Filter preFilterV0 = (nabs(aod::v0data::dcapostopv) > config.v0settingDcapostopv && nabs(aod::v0data::dcanegtopv) > config.v0settingDcanegtopv);

  // Defining the type of the daughter tracks
  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::MultZeqs, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::Mults>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTOFFullPi>>;
  using V0TrackCandidate = aod::V0Datas;
  // For Monte Carlo
  using EventCandidatesMC = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::CentFT0Cs, aod::CentFT0Ms>;
  using TrackCandidatesMC = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::McTrackLabels>>;
  using V0TrackCandidatesMC = soa::Join<aod::V0Datas, aod::McV0Labels>;

  template <typename T2, typename T3, typename T4>
  void fillInvMass(const T2& lv3, const T3& lv5, float multiplicity, const T4& daughter1, bool isMix)
  {

    // polarization calculations

    fourVecDau1 = ROOT::Math::PxPyPzMVector(daughter1.Px(), daughter1.Py(), daughter1.Pz(), o2::constants::physics::MassK0Short); // Kaon or Pion

    fourVecMother = ROOT::Math::PxPyPzMVector(lv3.Px(), lv3.Py(), lv3.Pz(), lv3.M()); // mass of KshortKshort pair
    ROOT::Math::Boost boost{fourVecMother.BoostToCM()};                               // boost mother to center of mass frame
    fourVecDauCM = boost(fourVecDau1);                                                // boost the frame of daughter same as mother
    threeVecDauCM = fourVecDauCM.Vect();                                              // get the 3 vector of daughter in the frame of mother

    if (std::abs(lv3.Rapidity()) < 0.5) {
      if (config.activateTHnSparseCosThStarHelicity) {
        helicityVec = fourVecMother.Vect(); // 3 vector of mother in COM frame
        auto cosThetaStarHelicity = helicityVec.Dot(threeVecDauCM) / (std::sqrt(threeVecDauCM.Mag2()) * std::sqrt(helicityVec.Mag2()));
        if (!isMix) {
          hglue.fill(HIST("h3glueInvMassDS"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarHelicity);

          for (int i = 0; i < config.cRotations; i++) {
            theta2 = rn->Uniform(o2::constants::math::PI - o2::constants::math::PI / config.rotationalCut, o2::constants::math::PI + o2::constants::math::PI / config.rotationalCut);
            hglue.fill(HIST("h3glueInvMassRot"), multiplicity, lv5.Pt(), lv5.M(), cosThetaStarHelicity);
          }
        } else {
          hglue.fill(HIST("h3glueInvMassME"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarHelicity);
        }

      } else if (config.activateTHnSparseCosThStarProduction) {
        normalVec = ROOT::Math::XYZVector(lv3.Py(), -lv3.Px(), 0.f);
        auto cosThetaStarProduction = normalVec.Dot(threeVecDauCM) / (std::sqrt(threeVecDauCM.Mag2()) * std::sqrt(normalVec.Mag2()));
        if (!isMix) {
          hglue.fill(HIST("h3glueInvMassDS"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarProduction);
          for (int i = 0; i < config.cRotations; i++) {
            theta2 = rn->Uniform(0, o2::constants::math::PI);
            hglue.fill(HIST("h3glueInvMassRot"), multiplicity, lv5.Pt(), lv5.M(), cosThetaStarProduction);
          }
        } else {
          hglue.fill(HIST("h3glueInvMassME"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarProduction);
        }

      } else if (config.activateTHnSparseCosThStarBeam) {
        beamVec = ROOT::Math::XYZVector(0.f, 0.f, 1.f);
        auto cosThetaStarBeam = beamVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());
        if (!isMix) {
          hglue.fill(HIST("h3glueInvMassDS"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarBeam);
          for (int i = 0; i < config.cRotations; i++) {
            theta2 = rn->Uniform(0, o2::constants::math::PI);
            hglue.fill(HIST("h3glueInvMassRot"), multiplicity, lv5.Pt(), lv5.M(), cosThetaStarBeam);
          }
        } else {
          hglue.fill(HIST("h3glueInvMassME"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarBeam);
        }

      } else if (config.activateTHnSparseCosThStarRandom) {
        auto phiRandom = gRandom->Uniform(0.f, constants::math::TwoPI);
        auto thetaRandom = gRandom->Uniform(0.f, constants::math::PI);

        randomVec = ROOT::Math::XYZVector(std::sin(thetaRandom) * std::cos(phiRandom), std::sin(thetaRandom) * std::sin(phiRandom), std::cos(thetaRandom));
        auto cosThetaStarRandom = randomVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());
        if (!isMix) {
          hglue.fill(HIST("h3glueInvMassDS"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarRandom);
          for (int i = 0; i < config.cRotations; i++) {
            theta2 = rn->Uniform(0, o2::constants::math::PI);
            hglue.fill(HIST("h3glueInvMassRot"), multiplicity, lv5.Pt(), lv5.M(), cosThetaStarRandom);
          }
        } else {
          hglue.fill(HIST("h3glueInvMassME"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarRandom);
        }
      }
    }
  }

  void processSE(EventCandidates::iterator const& collision, TrackCandidates const& /*tracks*/, aod::V0Datas const& V0s)
  {
    hglue.fill(HIST("heventscheck"), 0.5);
    multiplicity = 0.0;
    if (config.cfgMultFOTM) {
      multiplicity = collision.centFT0M();
    } else {
      multiplicity = collision.centFT0C();
    }
    if (!eventselection(collision)) {
      return;
    }

    // auto occupancyNumber = collision.trackOccupancyInTimeRange();
    // if (applyOccupancyCut && occupancyNumber < occupancyCut) {
    //   return;
    // }

    if (config.qAevents) {
      rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());
      rEventSelection.fill(HIST("hmultiplicity"), multiplicity);
      // rEventSelection.fill(HIST("multdist_FT0M"), collision.multFT0M());
      // rEventSelection.fill(HIST("multdist_FT0A"), collision.multFT0A());
      // rEventSelection.fill(HIST("multdist_FT0C"), collision.multFT0C());
      // rEventSelection.fill(HIST("hNcontributor"), collision.numContrib());
    }

    std::vector<int> v0indexes;
    bool allConditionsMet = 0;

    for (const auto& [v1, v2] : combinations(CombinationsFullIndexPolicy(V0s, V0s))) {

      if (v1.size() == 0 || v2.size() == 0) {
        continue;
      }

      if (!selectionV0(collision, v1, multiplicity)) {
        continue;
      }
      if (!selectionV0(collision, v2, multiplicity)) {
        continue;
      }

      auto postrack1 = v1.template posTrack_as<TrackCandidates>();
      auto negtrack1 = v1.template negTrack_as<TrackCandidates>();
      auto postrack2 = v2.template posTrack_as<TrackCandidates>();
      auto negtrack2 = v2.template negTrack_as<TrackCandidates>();

      double nTPCSigmaPos1{postrack1.tpcNSigmaPi()};
      double nTPCSigmaNeg1{negtrack1.tpcNSigmaPi()};
      double nTPCSigmaPos2{postrack2.tpcNSigmaPi()};
      double nTPCSigmaNeg2{negtrack2.tpcNSigmaPi()};

      if (!(isSelectedV0Daughter(negtrack1, -1, nTPCSigmaNeg1, v1) && isSelectedV0Daughter(postrack1, 1, nTPCSigmaPos1, v1))) {
        continue;
      }
      if (!(isSelectedV0Daughter(postrack2, 1, nTPCSigmaPos2, v2) && isSelectedV0Daughter(negtrack2, -1, nTPCSigmaNeg2, v2))) {
        continue;
      }

      if (std::find(v0indexes.begin(), v0indexes.end(), v1.globalIndex()) == v0indexes.end()) {
        v0indexes.push_back(v1.globalIndex());
      }
      // if (!(std::find(v0indexes.begin(), v0indexes.end(), v2.globalIndex()) != v0indexes.end())) {
      //   v0indexes.push_back(v2.globalIndex());
      // }

      if (v2.globalIndex() <= v1.globalIndex()) {
        continue;
      }

      // if (config.qAv0Daughters) {
      //   rKzeroShort.fill(HIST("negative_pt"), negtrack1.pt());
      //   rKzeroShort.fill(HIST("positive_pt"), postrack1.pt());
      //   rKzeroShort.fill(HIST("negative_eta"), negtrack1.eta());
      //   rKzeroShort.fill(HIST("positive_eta"), postrack1.eta());
      //   rKzeroShort.fill(HIST("negative_phi"), negtrack1.phi());
      //   rKzeroShort.fill(HIST("positive_phi"), postrack1.phi());
      // }

      if (postrack1.globalIndex() == postrack2.globalIndex()) {
        continue;
      }
      if (negtrack1.globalIndex() == negtrack2.globalIndex()) {
        continue;
      }

      if (!applyAngSep(v1, v2)) {
        continue;
      }

      if (config.qAv0) {
        rKzeroShort.fill(HIST("hMasscorrelationbefore"), v1.mK0Short(), v2.mK0Short());
      }
      allConditionsMet = 1;
      daughter1 = ROOT::Math::PxPyPzMVector(v1.px(), v1.py(), v1.pz(), o2::constants::physics::MassK0Short); // Kshort

      lv3.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
      lv5.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
      lv1.SetPtEtaPhiM(v1.pt(), v1.eta(), v1.phi(), o2::constants::physics::MassK0Short);
      lv2.SetPtEtaPhiM(v2.pt(), v2.eta(), v2.phi(), o2::constants::physics::MassK0Short);
      lv4.SetPtEtaPhiM(v1.pt(), v1.eta(), v1.phi() + theta2, o2::constants::physics::MassK0Short); // for rotated background
      lv3 = lv1 + lv2;
      lv5 = lv2 + lv4;
      isMix = false;

      if (!config.selectTWOKsOnly)
        fillInvMass(lv3, lv5, multiplicity, daughter1, isMix);
    }
    int sizeofv0indexes = v0indexes.size();
    rKzeroShort.fill(HIST("NksProduced"), sizeofv0indexes);
    if (config.selectTWOKsOnly && sizeofv0indexes == 2 && allConditionsMet) {
      fillInvMass(lv3, lv5, multiplicity, daughter1, false);
    }
    v0indexes.clear();
  }

  PROCESS_SWITCH(HigherMassResonances, processSE, "same event process", true);

  array<float, 3> pvec0;
  array<float, 3> pvec1;
  // use any one of 3 alias depending on the dataset. If pp then FT0M and if pbpb then FTOC
  using BinningTypeTPCMultiplicity = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultTPC>;
  using BinningTypeCentralityM = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;
  using BinningTypeVertexContributor = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
  ConfigurableAxis mevz = {"mevz", {10, -10., 10.}, "mixed event vertex z binning"};
  ConfigurableAxis memult = {"memult", {2000, 0, 10000}, "mixed event multiplicity binning"};

  void processME(EventCandidates const& collisions, TrackCandidates const& /*tracks*/, V0TrackCandidate const& v0s)
  {
    auto tracksTuple = std::make_tuple(v0s);
    BinningTypeVertexContributor binningOnPositions1{{mevz, memult}, true};
    BinningTypeCentralityM binningOnPositions2{{mevz, memult}, true};

    SameKindPair<EventCandidates, V0TrackCandidate, BinningTypeVertexContributor> pair1{binningOnPositions1, config.cfgNmixedEvents, -1, collisions, tracksTuple, &cache}; // for PbPb
    SameKindPair<EventCandidates, V0TrackCandidate, BinningTypeCentralityM> pair2{binningOnPositions2, config.cfgNmixedEvents, -1, collisions, tracksTuple, &cache};       // for pp

    if (config.cfgMultFOTM) {
      for (const auto& [c1, tracks1, c2, tracks2] : pair2) // two different centrality c1 and c2 and tracks corresponding to them
      {

        multiplicity = 0.0;
        multiplicity = c1.centFT0M();

        if (!eventselection(c1) || !eventselection(c2)) {
          continue;
        }
        // auto occupancyNumber = c1.trackOccupancyInTimeRange();
        // auto occupancyNumber2 = c2.trackOccupancyInTimeRange();
        // if (applyOccupancyCut && (occupancyNumber < occupancyCut || occupancyNumber2 < occupancyCut)) {
        //   return;
        // }

        for (const auto& [t1, t2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {

          if (t1.size() == 0 || t2.size() == 0) {
            continue;
          }

          if (!selectionV0(c1, t1, multiplicity))
            continue;
          if (!selectionV0(c2, t2, multiplicity))
            continue;

          auto postrack1 = t1.template posTrack_as<TrackCandidates>();
          auto negtrack1 = t1.template negTrack_as<TrackCandidates>();
          auto postrack2 = t2.template posTrack_as<TrackCandidates>();
          auto negtrack2 = t2.template negTrack_as<TrackCandidates>();
          if (postrack1.globalIndex() == postrack2.globalIndex()) {
            continue;
          }
          if (negtrack1.globalIndex() == negtrack2.globalIndex()) {
            continue;
          }
          double nTPCSigmaPos1{postrack1.tpcNSigmaPi()};
          double nTPCSigmaNeg1{negtrack1.tpcNSigmaPi()};
          double nTPCSigmaPos2{postrack2.tpcNSigmaPi()};
          double nTPCSigmaNeg2{negtrack2.tpcNSigmaPi()};

          if (!isSelectedV0Daughter(postrack1, 1, nTPCSigmaPos1, t1)) {
            continue;
          }
          if (!isSelectedV0Daughter(postrack2, 1, nTPCSigmaPos2, t2)) {
            continue;
          }
          if (!isSelectedV0Daughter(negtrack1, -1, nTPCSigmaNeg1, t1)) {
            continue;
          }
          if (!isSelectedV0Daughter(negtrack2, -1, nTPCSigmaNeg2, t2)) {
            continue;
          }

          daughter1 = ROOT::Math::PxPyPzMVector(t1.px(), t1.py(), t1.pz(), o2::constants::physics::MassK0Short); // Kshort

          lv3.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
          lv1.SetPtEtaPhiM(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassK0Short);
          lv2.SetPtEtaPhiM(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassK0Short);
          lv4.SetPtEtaPhiM(t1.pt(), t1.eta(), t1.phi() + theta2, o2::constants::physics::MassK0Short); // for rotated background
          lv3 = lv1 + lv2;
          lv5 = lv2 + lv4;
          isMix = true;
          fillInvMass(lv3, lv5, multiplicity, daughter1, isMix);
        }
      }
    } else {
      for (const auto& [c1, tracks1, c2, tracks2] : pair1) // two different centrality c1 and c2 and tracks corresponding to them
      {
        multiplicity = 0.0f;
        multiplicity = c1.centFT0C();

        if (!eventselection(c1) || !eventselection(c2)) {
          continue;
        }
        // auto occupancyNumber = c1.trackOccupancyInTimeRange();
        // auto occupancyNumber2 = c2.trackOccupancyInTimeRange();
        // if (applyOccupancyCut && (occupancyNumber < occupancyCut || occupancyNumber2 < occupancyCut)) {
        //   return;
        // }

        for (const auto& [t1, t2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
          if (t1.size() == 0 || t2.size() == 0) {
            continue;
          }

          if (!selectionV0(c1, t1, multiplicity))
            continue;
          if (!selectionV0(c2, t2, multiplicity))
            continue;

          auto postrack1 = t1.template posTrack_as<TrackCandidates>();
          auto negtrack1 = t1.template negTrack_as<TrackCandidates>();
          auto postrack2 = t2.template posTrack_as<TrackCandidates>();
          auto negtrack2 = t2.template negTrack_as<TrackCandidates>();
          if (postrack1.globalIndex() == postrack2.globalIndex()) {
            continue;
          }
          if (negtrack1.globalIndex() == negtrack2.globalIndex()) {
            continue;
          }
          double nTPCSigmaPos1{postrack1.tpcNSigmaPi()};
          double nTPCSigmaNeg1{negtrack1.tpcNSigmaPi()};
          double nTPCSigmaPos2{postrack2.tpcNSigmaPi()};
          double nTPCSigmaNeg2{negtrack2.tpcNSigmaPi()};

          if (!isSelectedV0Daughter(postrack1, 1, nTPCSigmaPos1, t1)) {
            continue;
          }
          if (!isSelectedV0Daughter(postrack2, 1, nTPCSigmaPos2, t2)) {
            continue;
          }
          if (!isSelectedV0Daughter(negtrack1, -1, nTPCSigmaNeg1, t1)) {
            continue;
          }
          if (!isSelectedV0Daughter(negtrack2, -1, nTPCSigmaNeg2, t2)) {
            continue;
          }

          daughter1 = ROOT::Math::PxPyPzMVector(t1.px(), t1.py(), t1.pz(), o2::constants::physics::MassK0Short); // Kshort

          lv3.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
          lv1.SetPtEtaPhiM(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassK0Short);
          lv2.SetPtEtaPhiM(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassK0Short);
          lv4.SetPtEtaPhiM(t1.pt(), t1.eta(), t1.phi() + theta2, o2::constants::physics::MassK0Short); // for rotated background
          lv3 = lv1 + lv2;
          lv5 = lv2 + lv4;
          isMix = true;
          fillInvMass(lv3, lv5, multiplicity, daughter1, isMix);
        }
      }
    }
  }
  PROCESS_SWITCH(HigherMassResonances, processME, "mixed event process", true);

  int counter = 0;
  float multiplicityGen = 0.0;
  void processGen(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles, const soa::SmallGroups<EventCandidatesMC>& collisions)
  {
    TLorentzVector genvec;
    hMChists.fill(HIST("events_check"), 0.5);
    if (std::abs(mcCollision.posZ()) < config.cutzvertex) {
      hMChists.fill(HIST("events_check"), 1.5);
    }
    // int Nchinel = 0;
    // for (const auto& mcParticle : mcParticles) {
    //   auto pdgcode = std::abs(mcParticle.pdgCode());
    //   if (mcParticle.isPhysicalPrimary() && (pdgcode == 211 || pdgcode == 321 || pdgcode == 2212 || pdgcode == 11 || pdgcode == 13)) {
    //     if (std::abs(mcParticle.eta()) < 1.0) {
    //       Nchinel = Nchinel + 1;
    //     }
    //   }
    // }
    // if (Nchinel > 0 && std::abs(mcCollision.posZ()) < config.cutzvertex)
    hMChists.fill(HIST("events_check"), 2.5);

    std::vector<int64_t> selectedEvents(collisions.size());
    int nevts = 0;
    multiplicityGen = 0.0;
    for (const auto& collision : collisions) {
      if (std::abs(collision.mcCollision().posZ()) > config.cutzvertex) {
        continue;
      }

      if (config.timFrameEvsel && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
        continue;
      }
      if (config.cTVXEvsel && (!collision.selection_bit(aod::evsel::kIsTriggerTVX))) {
        continue;
      }

      multiplicityGen = collision.centFT0M();

      selectedEvents[nevts++] = collision.mcCollision_as<aod::McCollisions>().globalIndex();
    }
    selectedEvents.resize(nevts);
    hMChists.fill(HIST("events_check"), 3.5);
    const auto evtReconstructedAndSelected = std::find(selectedEvents.begin(), selectedEvents.end(), mcCollision.globalIndex()) != selectedEvents.end();

    if (!config.allGenCollisions && !evtReconstructedAndSelected) { // Check that the event is reconstructed and that the reconstructed events pass the selection
      return;
    }
    hMChists.fill(HIST("events_check"), 4.5);
    for (const auto& mcParticle : mcParticles) {

      if (std::abs(mcParticle.pdgCode()) != config.pdgCodes[config.selectMCparticles]) // f2(1525), f0(1710)
      {
        continue;
      }
      hMChists.fill(HIST("events_check"), 5.5);

      // if (counter < 1e3)
      //   std::cout << "px " << mcParticle.px() << " py " << mcParticle.py() << " pz " << mcParticle.pz() << " y " << mcParticle.y() << std::endl;
      // counter++;

      hMChists.fill(HIST("GenRapidity"), mcParticle.pt(), mcParticle.y());
      hMChists.fill(HIST("GenPhi"), mcParticle.phi());
      hMChists.fill(HIST("GenEta"), mcParticle.pt(), mcParticle.eta());
      // hMChists.fill(HIST("GenPx"), mcParticle.px());
      // hMChists.fill(HIST("GenPy"), mcParticle.py());
      // hMChists.fill(HIST("GenPz"), mcParticle.pz());

      if (std::abs(mcParticle.y()) >= 0.5) {
        continue;
      }
      hMChists.fill(HIST("events_check"), 6.5);

      auto kDaughters = mcParticle.daughters_as<aod::McParticles>();
      if (kDaughters.size() != 2) {
        continue;
      }
      hMChists.fill(HIST("events_check"), 7.5);

      auto passKs = false;
      for (const auto& kCurrentDaughter : kDaughters) {
        // int daupdg = std::abs(kCurrentDaughter.pdgCode());

        if (!kCurrentDaughter.isPhysicalPrimary()) {
          continue;
        }
        hMChists.fill(HIST("events_check"), 8.5);

        if (std::abs(kCurrentDaughter.pdgCode()) == 310) {
          passKs = true;
          hMChists.fill(HIST("events_check"), 9.5);
        }
      }
      if (passKs) {
        genvec.SetPtEtaPhiE(mcParticle.pt(), mcParticle.eta(), mcParticle.phi(), mcParticle.e());
        hMChists.fill(HIST("Genf1710_mass"), genvec.M());
        hMChists.fill(HIST("Genf1710"), multiplicityGen, mcParticle.pt());
      }
    }
  }
  PROCESS_SWITCH(HigherMassResonances, processGen, "Process Generated", false);

  int counter2 = 0;
  int eventCounter = 0;
  std::vector<int> gindex1, gindex2;
  void processRec(EventCandidatesMC::iterator const& collision, TrackCandidatesMC const&, V0TrackCandidatesMC const& V0s, aod::McParticles const&, aod::McCollisions const& /*mcCollisions*/)
  {

    TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonance;
    // lDecayDaughter1.SetXYZM(0.0,0.0,0.0,0.0);
    // lDecayDaughter2.SetXYZM(0.0,0.0,0.0,0.0);
    // lDecayDaughter1.SetXYZM(0.0,0.0,0.0,0.0);
    auto multiplicity = collision.centFT0C();
    hMChists.fill(HIST("MC_mult"), multiplicity);

    hMChists.fill(HIST("events_checkrec"), 0.5);
    if (!collision.has_mcCollision()) {
      return;
    }
    hMChists.fill(HIST("events_checkrec"), 1.5);
    // // if (std::abs(collision.mcCollision().posZ()) > config.cutzvertex || !collision.sel8()) {
    if (std::abs(collision.mcCollision().posZ()) > config.cutzvertex) {
      return;
    }
    hMChists.fill(HIST("events_checkrec"), 2.5);

    if (config.timFrameEvsel && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return;
    }
    hMChists.fill(HIST("events_checkrec"), 3.5);
    if (config.cTVXEvsel && (!collision.selection_bit(aod::evsel::kIsTriggerTVX))) {
      return;
    }
    hMChists.fill(HIST("events_checkrec"), 4.5);
    hMChists.fill(HIST("MC_mult_after_event_sel"), multiplicity);
    eventCounter++;
    auto oldindex = -999;

    for (const auto& v01 : V0s) {

      for (const auto& v02 : V0s) {

        hMChists.fill(HIST("events_checkrec"), 5.5);

        if (v02.index() <= v01.index()) {
          continue;
        }

        if (!v01.has_mcParticle() || !v02.has_mcParticle()) {
          continue;
        }
        hMChists.fill(HIST("events_checkrec"), 6.5);

        auto postrack1 = v01.template posTrack_as<TrackCandidatesMC>();
        auto negtrack1 = v01.template negTrack_as<TrackCandidatesMC>();

        auto postrack2 = v02.template posTrack_as<TrackCandidatesMC>();
        auto negtrack2 = v02.template negTrack_as<TrackCandidatesMC>();

        if (!postrack1.has_mcParticle() || !postrack2.has_mcParticle())
          continue; // Checking that the daughter tracks come from particles and are not fake
        hMChists.fill(HIST("events_checkrec"), 7.5);

        if (!negtrack1.has_mcParticle() || !negtrack2.has_mcParticle())
          continue;
        hMChists.fill(HIST("events_checkrec"), 8.5);

        double nTPCSigmaPos1[1]{postrack1.tpcNSigmaPi()};
        double nTPCSigmaNeg1[1]{negtrack1.tpcNSigmaPi()};

        double nTPCSigmaPos2[1]{postrack2.tpcNSigmaPi()};
        double nTPCSigmaNeg2[1]{negtrack2.tpcNSigmaPi()};

        if (!isSelectedV0Daughter(postrack1, 1, nTPCSigmaPos1[0], v01) || !isSelectedV0Daughter(postrack2, 1, nTPCSigmaPos2[0], v02)) {
          continue;
        }
        hMChists.fill(HIST("events_checkrec"), 9.5);

        if (!isSelectedV0Daughter(negtrack1, -1, nTPCSigmaNeg1[0], v01) || !isSelectedV0Daughter(negtrack2, -1, nTPCSigmaNeg2[0], v02)) {
          continue;
        }
        hMChists.fill(HIST("events_checkrec"), 10.5);

        if (!selectionV0(collision, v01, multiplicity) || !selectionV0(collision, v02, multiplicity)) {
          continue;
        }
        hMChists.fill(HIST("events_checkrec"), 11.5);

        auto mctrackv01 = v01.mcParticle();
        auto mctrackv02 = v02.mcParticle();

        int trackv0PDG1 = std::abs(mctrackv01.pdgCode());
        int trackv0PDG2 = std::abs(mctrackv02.pdgCode());

        if (std::abs(trackv0PDG1) != 310 || std::abs(trackv0PDG2) != 310) {
          continue;
        }
        hMChists.fill(HIST("events_checkrec"), 12.5);

        for (const auto& mothertrack1 : mctrackv01.mothers_as<aod::McParticles>()) {

          // int motpdgs = std::abs(mothertrack1.pdgCode());
          gindex1.push_back(mothertrack1.globalIndex());
          if (gindex1.size() > 1) {
            if (std::find(gindex1.begin(), gindex1.end(), mothertrack1.globalIndex()) != gindex1.end()) {
              continue;
            }
          }
          // if (counter2 < 1e4)
          //   std::cout << "Mother1 pdg code: " << motpdgs << " p_{T} " << mothertrack1.pt() << "Global index " << mothertrack1.globalIndex() << " event " << eventCounter << std::endl;
          // counter2++;

          // int counter_check = 0;

          for (const auto& mothertrack2 : mctrackv02.mothers_as<aod::McParticles>()) {

            hMChists.fill(HIST("events_checkrec"), 13.5);

            if (mothertrack1.pdgCode() != mothertrack2.pdgCode()) {
              continue;
            }
            hMChists.fill(HIST("events_checkrec"), 14.5);

            // int motpdgs2 = std::abs(mothertrack2.pdgCode());
            gindex2.push_back(mothertrack2.globalIndex());
            if (gindex2.size() > 1) {
              if (std::find(gindex2.begin(), gindex2.end(), mothertrack2.globalIndex()) != gindex2.end()) {
                continue;
              }
            }
            // if (counter2 < 1e4)
            //   std::cout << "Mother2 pdg code: " << motpdgs2 << " p_{T} " << mothertrack2.pt() << "Global index " << mothertrack1.globalIndex() << " event " << eventCounter << std::endl;

            if (mothertrack1.pdgCode() != config.pdgCodes[config.selectMCparticles]) {
              continue;
            }
            hMChists.fill(HIST("events_checkrec"), 15.5);

            if (mothertrack1.globalIndex() != mothertrack2.globalIndex()) {
              continue;
            }
            hMChists.fill(HIST("events_checkrec"), 16.5);

            if (!mothertrack1.producedByGenerator()) {
              continue;
            }
            hMChists.fill(HIST("events_checkrec"), 17.5);

            if (std::abs(mothertrack1.y()) >= 0.5) {
              continue;
            }
            hMChists.fill(HIST("events_checkrec"), 18.5);

            if (config.avoidsplitrackMC && oldindex == mothertrack1.globalIndex()) {
              hMChists.fill(HIST("h1Recsplit"), mothertrack1.pt());
              continue;
            }
            oldindex = mothertrack1.globalIndex();

            // counter_check++;
            // if (counter_check > 1) {
            //   std::cout << "Total mothers is " << counter_check << std::endl;
            // }
            // std::cout << "After selection " << " p_{T} " << mothertrack2.pt() << " event " << eventCounter << std::endl;

            pvec0 = std::array{v01.px(), v01.py(), v01.pz()};
            pvec1 = std::array{v02.px(), v02.py(), v02.pz()};
            auto arrMomrec = std::array{pvec0, pvec1};
            auto motherP = mothertrack1.p();
            // auto motherE = mothertrack1.e();
            // auto genMass = std::sqrt(motherE * motherE - motherP * motherP);
            auto recMass = RecoDecay::m(arrMomrec, std::array{o2::constants::physics::MassK0Short, o2::constants::physics::MassK0Short});
            // auto recpt = TMath::Sqrt((track1.px() + track2.px()) * (track1.px() + track2.px()) + (track1.py() + track2.py()) * (track1.py() + track2.py()));
            //// Resonance reconstruction
            lDecayDaughter1.SetXYZM(v01.px(), v01.py(), v01.pz(), o2::constants::physics::MassK0Short);
            lDecayDaughter2.SetXYZM(v02.px(), v02.py(), v02.pz(), o2::constants::physics::MassK0Short);
            lResonance = lDecayDaughter1 + lDecayDaughter2;

            hMChists.fill(HIST("Recf1710_p"), motherP);
            hMChists.fill(HIST("Recf1710_mass"), recMass);
            hMChists.fill(HIST("Recf1710_pt1"), multiplicity, mothertrack1.pt(), recMass);
            // hMChists.fill(HIST("Genf1710_mass"), genMass);
            hMChists.fill(HIST("Recf1710_pt2"), multiplicity, lResonance.Pt(), recMass);

            hMChists.fill(HIST("RecRapidity"), mothertrack1.y());
            hMChists.fill(HIST("RecPhi"), mothertrack1.phi());
            hMChists.fill(HIST("RecEta"), mothertrack1.eta());
          }
          gindex2.clear();
        }
        gindex1.clear();
      }
    }
  }
  PROCESS_SWITCH(HigherMassResonances, processRec, "Process Reconstructed", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HigherMassResonances>(cfgc)};
}
