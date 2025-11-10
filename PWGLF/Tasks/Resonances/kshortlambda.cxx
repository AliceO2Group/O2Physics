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

/// \file kshortlambda.cxx
/// \brief higher mass resonance search in non-identical V0 pairs (K0s-L)
/// \author  dukhishyam Mallick (dukhishyam.mallick@cern.ch)

#include "PWGLF/DataModel/LFStrangenessTables.h" //

#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h" //
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h" //
#include "Common/DataModel/PIDResponseTPC.h" //
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h" //
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h" //
#include "ReconstructionDataFormats/Track.h"

#include "Math/GenVector/Boost.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TF1.h"
#include "TRandom3.h"
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

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
// using namespace o2::constants::physics;
using std::array;

struct Kshortlambda {
  SliceCache cache;
  HistogramRegistry rEventSelection{"eventSelection", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rKzeroShort{"kzeroShort", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry hvzero{"hvzeroball", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  Configurable<bool> qAv0{"qAv0", false, "qAv0"};
  Configurable<bool> qAPID{"qAPID", true, "qAPID"};
  Configurable<bool> qAv0daughters{"qAv0daughters", false, "qA of v0 daughters"};
  Configurable<bool> qAevents{"qAevents", false, "QA of events"};
  Configurable<bool> invMass1D{"invMass1D", false, "1D invariant mass histograms"};
  Configurable<bool> correlation2Dhist{"correlation2Dhist", true, "Lamda K0 mass correlation"};
  Configurable<bool> cDCAv0topv{"cDCAv0topv", false, "DCA V0 to PV"};
  Configurable<bool> armcut{"armcut", true, "arm cut"};
  Configurable<bool> globalTracks{"globalTracks", false, "Global tracks"};
  Configurable<bool> hasTPC{"hasTPC", false, "TPC"};

  // Configurable for event selection
  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range (cm)"};
  Configurable<float> cfgETAcut{"cfgETAcut", 0.8f, "Track ETA cut"};
  Configurable<bool> timFrameEvsel{"timFrameEvsel", false, "TPC Time frame boundary cut"};
  Configurable<bool> piluprejection{"piluprejection", false, "Pileup rejection"};
  Configurable<bool> goodzvertex{"goodzvertex", false, "removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference."};
  Configurable<bool> itstpctracks{"itstpctracks", false, "selects collisions with at least one ITS-TPC track,"};
  Configurable<bool> additionalEvsel{"additionalEvsel", false, "Additional event selcection"};
  Configurable<bool> applyOccupancyCut{"applyOccupancyCut", false, "Apply occupancy cut"};
  Configurable<int> occupancyCut{"occupancyCut", 1000, "Mimimum Occupancy cut"};

  // Configurable parameters for V0 selection
  Configurable<float> confV0DCADaughMax{"confV0DCADaughMax", 1.0f, "DCA b/w V0 daughters"};
  Configurable<float> v0settingdcapostopv{"v0settingdcapostopv", 0.06, "DCA Pos To PV"};
  Configurable<float> v0settingdcanegtopv{"v0settingdcanegtopv", 0.06, "DCA Neg To PV"};
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

  // Configurable for track selection and multiplicity
  Configurable<float> cfgPTcut{"cfgPTcut", 0.2f, "Track PT cut"};
  Configurable<int> cfgNmixedEvents{"cfgNmixedEvents", 5, "Number of mixed events"};
  Configurable<bool> cfgMultFOTM{"cfgMultFOTM", true, "Use FOTM multiplicity if pp else use 0 here for PbPb (FT0C)"};
  ConfigurableAxis binsCent{"binsCent", {VARIABLE_WIDTH, 0., 5., 10., 30., 50., 70., 100., 110., 150.}, "Binning of the centrality axis"};

  // output THnSparses
  Configurable<bool> activateTHnSparseCosThStarHelicity{"activateTHnSparseCosThStarHelicity", false, "Activate the THnSparse with cosThStar w.r.t. helicity axis"};
  Configurable<bool> activateTHnSparseCosThStarProduction{"activateTHnSparseCosThStarProduction", false, "Activate the THnSparse with cosThStar w.r.t. production axis"};
  Configurable<bool> activateTHnSparseCosThStarBeam{"activateTHnSparseCosThStarBeam", true, "Activate the THnSparse with cosThStar w.r.t. beam axis (Gottified jackson frame)"};
  Configurable<bool> activateTHnSparseCosThStarRandom{"activateTHnSparseCosThStarRandom", false, "Activate the THnSparse with cosThStar w.r.t. random axis"};
  Configurable<int> cnofrotations{"cnofrotations", 3, "Number of random rotations in the rotational background"};

  // Other cuts on Ks and glueball
  Configurable<bool> rapidityks{"rapidityks", true, "rapidity cut on K0s"};
  Configurable<bool> applyCompetingcut{"applyCompetingcut", false, "Competing cascade rejection cut"};
  Configurable<float> competingCascrejlambda{"competingCascrejlambda", 0.005, "rejecting competing cascade lambda"};
  Configurable<float> competingCascrejlambdaanti{"competingCascrejlambdaanti", 0.005, "rejecting competing cascade anti-lambda"}; // If one of the pions is misidentified as a proton, then instead of Ks we reconstruct lambda, therefore the competing cascade rejection cut is applied in which if the reconstrcted mass of a pion and proton (which we are assuming to be misidentified as proton) is close to lambda or anti-lambda, then the track is rejected.
  Configurable<int> tpcCrossedrows{"tpcCrossedrows", 70, "TPC crossed rows"};
  Configurable<float> tpcCrossedrowsOverfcls{"tpcCrossedrowsOverfcls", 0.8, "TPC crossed rows over findable clusters"};

  // Mass and pT axis as configurables
  Configurable<float> cPtMin{"cPtMin", 0.0f, "Minimum pT"};
  Configurable<float> cPtMax{"cPtMax", 50.0f, "Maximum pT"};
  Configurable<int> cPtBins{"cPtBins", 500, "Number of pT bins"};
  Configurable<float> cMassMin{"cMassMin", 0.9f, "Minimum mass bin"};
  Configurable<float> cMassMax{"cMassMax", 3.0f, "Maximum mass bin"};
  Configurable<int> cMassBins{"cMassBins", 210, "Number of mass binsl"};
  Configurable<float> ksMassMin{"ksMassMin", 0.45f, "Minimum mass of K0s"};
  Configurable<float> ksMassMax{"ksMassMax", 0.55f, "Maximum mass of K0s"};
  Configurable<int> ksMassBins{"ksMassBins", 200, "Number of mass bins for K0s"};
  Configurable<int> rotationalCut{"rotationalCut", 10, "Cut value (Rotation angle pi - pi/cut and pi + pi/cut)"};
  ConfigurableAxis configThnAxisPOL{"configThnAxisPOL", {20, -1.0, 1.0}, "Costheta axis"};
  ConfigurableAxis axisdEdx{"axisdEdx", {20000, 0.0f, 200.0f}, "dE/dx (a.u.)"};
  ConfigurableAxis axisPtfordEbydx{"axisPtfordEbydx", {2000, 0, 20}, "pT (GeV/c)"};
  ConfigurableAxis axisMultdist{"axisMultdist", {3500, 0, 70000}, "Multiplicity distribution"};
  ConfigurableAxis occupancyBins{"occupancyBins", {VARIABLE_WIDTH, 0.0, 100, 500, 600, 1000, 1100, 1500, 1600, 2000, 2100, 2500, 2600, 3000, 3100, 3500, 3600, 4000, 4100, 4500, 4600, 5000, 5100, 9999}, "Binning of the occupancy axis"};

  // Event selection cuts - Alex (Temporary, need to fix!)
  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;
  TF1* fMultCutLow = nullptr;
  TF1* fMultCutHigh = nullptr;
  TF1* fMultMultPVCut = nullptr;

  TRandom* rn = new TRandom();

  void init(InitContext const&)
  {
    // Axes
    AxisSpec k0ShortMassAxis = {ksMassBins, ksMassMin, ksMassMax, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec vzeroMassAxis = {cMassBins, cMassMin, cMassMax, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec vertexZAxis = {60, -15.f, 15.f, "vrtx_{Z} [cm]"}; // for histogram
    AxisSpec ptAxis = {cPtBins, cPtMin, cPtMax, "#it{p}_{T} (GeV/#it{c})"};
    // AxisSpec multiplicityAxis = {110, 0.0f, 150.0f, "Multiplicity Axis"};
    AxisSpec multiplicityAxis = {binsCent, "Multiplicity Axis"};
    AxisSpec thnAxisPOL{configThnAxisPOL, "Configurabel theta axis"};
    AxisSpec occupancyAxis = {occupancyBins, "Occupancy [-40,100]"};

    //  THnSparses
    std::array<int, 4> sparses = {activateTHnSparseCosThStarHelicity, activateTHnSparseCosThStarProduction, activateTHnSparseCosThStarBeam, activateTHnSparseCosThStarRandom};
    // std::array<int, 1> sparses = {activateTHnSparseCosThStarHelicity};
    if (std::accumulate(sparses.begin(), sparses.end(), 0) == 0) {
      LOGP(fatal, "No output THnSparses enabled");
    } else {
      if (activateTHnSparseCosThStarHelicity) {
        LOGP(info, "THnSparse with cosThStar w.r.t. helicity axis active.");
      }
      if (activateTHnSparseCosThStarProduction) {
        LOGP(info, "THnSparse with cosThStar w.r.t. production axis active.");
      }
      if (activateTHnSparseCosThStarBeam) {
        LOGP(info, "THnSparse with cosThStar w.r.t. beam axis active. (Gottified jackson frame)");
      }
      if (activateTHnSparseCosThStarRandom) {
        LOGP(info, "THnSparse with cosThStar w.r.t. random axis active.");
      }
    }

    // Event selection
    if (qAevents) {
      rEventSelection.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});
      rEventSelection.add("hmultiplicity", "multiplicity percentile distribution", {HistType::kTH1F, {{150, 0.0f, 150.0f}}});
      rEventSelection.add("multdist_FT0M", "FT0M Multiplicity distribution", kTH1F, {axisMultdist});
      rEventSelection.add("multdist_FT0A", "FT0A Multiplicity distribution", kTH1F, {axisMultdist});
      rEventSelection.add("multdist_FT0C", "FT0C Multiplicity distribution", kTH1F, {axisMultdist});
      rEventSelection.add("hNcontributor", "Number of primary vertex contributor", kTH1F, {{2000, 0.0f, 10000.0f}});
    }

    if (invMass1D) {
      hvzero.add("h3vzeropairInvMassDS", "h1vzeropairInvMassDS", kTH1F, {vzeroMassAxis});
      hvzero.add("h3vzeropairInvMassME", "h1vzeropairInvMassME", kTH1F, {vzeroMassAxis});
      hvzero.add("h3vzeropairInvMassRot", "h1vzeropairInvMassRot", kTH1F, {vzeroMassAxis});
    }
    hvzero.add("h3vzeropairInvMassDS", "h3vzeropairInvMassDS", kTHnSparseF, {multiplicityAxis, ptAxis, vzeroMassAxis, thnAxisPOL, occupancyAxis}, true);
    hvzero.add("h3vzeropairInvMassME", "h3vzeropairInvMassME", kTHnSparseF, {multiplicityAxis, ptAxis, vzeroMassAxis, thnAxisPOL, occupancyAxis}, true);
    hvzero.add("h3vzeropairInvMassRot", "h3vzeropairInvMassRot", kTHnSparseF, {multiplicityAxis, ptAxis, vzeroMassAxis, thnAxisPOL, occupancyAxis}, true);
    hvzero.add("heventscheck", "heventscheck", kTH1I, {{10, 0, 10}});
    hvzero.add("htrackscheckv0", "htrackscheckv0", kTH1I, {{15, 0, 15}});
    hvzero.add("htrackscheckv0daughters", "htrackscheckv0daughters", kTH1I, {{15, 0, 15}});

    // K0s topological/PID cuts
    if (correlation2Dhist) {
      rKzeroShort.add("mass_lambda_kshort_before", "mass under lambda hypotheses and Kshort mass", kTH2F, {{100, 0.2, 0.8}, {100, 0.9, 1.5}});
      rKzeroShort.add("mass_lambda_kshort_after1", "mass under lambda hypotheses and Kshort mass", kTH2F, {{100, 0.2, 0.8}, {100, 0.9, 1.5}});
      rKzeroShort.add("mass_lambda_kshort_after2", "mass under lambda hypotheses and Kshort mass", kTH2F, {{100, 0.2, 0.8}, {100, 0.9, 1.5}});
      rKzeroShort.add("mass_lambda_kshort_after3", "mass under lambda hypotheses and Kshort mass", kTH2F, {{100, 0.2, 0.8}, {100, 0.9, 1.5}});
      rKzeroShort.add("mass_lambda_kshort_after4", "mass under lambda hypotheses and Kshort mass", kTH2F, {{100, 0.2, 0.8}, {100, 0.9, 1.5}});
      rKzeroShort.add("mass_lambda_kshort_after5", "mass under lambda hypotheses and Kshort mass", kTH2F, {{100, 0.2, 0.8}, {100, 0.9, 1.5}});
      rKzeroShort.add("mass_lambda_kshort_after6", "mass under lambda hypotheses and Kshort mass", kTH2F, {{100, 0.2, 0.8}, {100, 0.9, 1.5}});
      rKzeroShort.add("mass_lambda_kshort_after7", "mass under lambda hypotheses and Kshort mass", kTH2F, {{100, 0.2, 0.8}, {100, 0.9, 1.5}});
      rKzeroShort.add("mass_lambda_kshort_after8", "mass under lambda hypotheses and Kshort mass", kTH2F, {{100, 0.2, 0.8}, {100, 0.9, 1.5}});
      rKzeroShort.add("mass_lambda_kshort_after9", "mass under lambda hypotheses and Kshort mass", kTH2F, {{100, 0.2, 0.8}, {100, 0.9, 1.5}});
      rKzeroShort.add("mass_lambda_kshort_after10", "mass under lambda hypotheses and Kshort mass", kTH2F, {{100, 0.2, 0.8}, {100, 0.9, 1.5}});
    }
    if (qAv0) {
      // Invariant Mass
      rKzeroShort.add("hMassK0Shortbefore", "hMassK0Shortbefore", kTHnSparseF, {k0ShortMassAxis, ptAxis});
      rKzeroShort.add("hMasscorrelationbefore", "hMasscorrelationbefore", kTH2F, {k0ShortMassAxis, k0ShortMassAxis});
      rKzeroShort.add("hMassK0ShortSelected", "hMassK0ShortSelected", kTHnSparseF, {k0ShortMassAxis, ptAxis});

      // Topological histograms (after the selection)
      rKzeroShort.add("hDCAV0Daughters", "DCA between v0 daughters", {HistType::kTH1F, {{60, -3.0f, 3.0f}}});
      rKzeroShort.add("hV0CosPA", "hV0CosPA", {HistType::kTH1F, {{100, 0.96f, 1.1f}}});
      rKzeroShort.add("hLT", "hLT", {HistType::kTH1F, {{100, 0.0f, 50.0f}}});
      rKzeroShort.add("Mass_lambda", "Mass under lambda hypothesis", kTH1F, {vzeroMassAxis});
      rKzeroShort.add("mass_AntiLambda", "Mass under anti-lambda hypothesis", kTH1F, {vzeroMassAxis});
      rKzeroShort.add("mass_Gamma", "Mass under Gamma hypothesis", kTH1F, {vzeroMassAxis});

      rKzeroShort.add("rapidity", "Rapidity distribution", kTH1F, {{100, -1.0f, 1.0f}});
      rKzeroShort.add("hv0radius", "hv0radius", kTH1F, {{100, 0.0f, 200.0f}});
      rKzeroShort.add("hDCApostopv", "DCA positive daughter to PV", kTH1F, {{1000, -10.0f, 10.0f}});
      rKzeroShort.add("hDCAnegtopv", "DCA negative daughter to PV", kTH1F, {{1000, -10.0f, 10.0f}});
      rKzeroShort.add("hDCAv0topv", "DCA V0 to PV", kTH1F, {{60, -3.0f, 3.0f}});
      rKzeroShort.add("halpha", "Armenteros alpha", kTH1F, {{100, -5.0f, 5.0f}});
      rKzeroShort.add("hqtarmbyalpha", "qtarm/alpha", kTH1F, {{100, 0.0f, 1.0f}});
      rKzeroShort.add("hpsipair", "psi pair angle", kTH1F, {{100, -5.0f, 5.0f}});
      rKzeroShort.add("NksProduced", "Number of K0s produced", kTH1I, {{15, 0, 15}});
    }
    if (qAPID) {
      rKzeroShort.add("hNSigmaPosPionK0s_before", "hNSigmaPosPionK0s_before", {HistType::kTH2F, {{ptAxis}, {100, -5.f, 5.f}}});
      rKzeroShort.add("hNSigmaNegPionK0s_before", "hNSigmaNegPionK0s_before", {HistType::kTH2F, {{ptAxis}, {100, -5.f, 5.f}}});
      rKzeroShort.add("dE_by_dx_TPC", "dE/dx signal in the TPC as a function of pT", kTH2F, {axisPtfordEbydx, axisdEdx});
    }
    if (qAv0daughters) {
      rKzeroShort.add("negative_pt", "Negative daughter pT", kTH1F, {ptAxis});
      rKzeroShort.add("positive_pt", "Positive daughter pT", kTH1F, {ptAxis});
      rKzeroShort.add("negative_eta", "Negative daughter eta", kTH1F, {{100, -1.0f, 1.0f}});
      rKzeroShort.add("positive_eta", "Positive daughter eta", kTH1F, {{100, -1.0f, 1.0f}});
      rKzeroShort.add("negative_phi", "Negative daughter phi", kTH1F, {{70, 0.0f, 7.0f}});
      rKzeroShort.add("positive_phi", "Positive daughter phi", kTH1F, {{70, 0.0f, 7.0f}});
    }
    if (additionalEvsel) {
      fMultPVCutLow = new TF1("fMultPVCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 2.5*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultPVCutLow->SetParameters(2834.66, -87.0127, 0.915126, -0.00330136, 332.513, -12.3476, 0.251663, -0.00272819, 1.12242e-05);
      fMultPVCutHigh = new TF1("fMultPVCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 2.5*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultPVCutHigh->SetParameters(2834.66, -87.0127, 0.915126, -0.00330136, 332.513, -12.3476, 0.251663, -0.00272819, 1.12242e-05);
      // fMultCutLow = new TF1("fMultCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 2.5*([4]+[5]*x)", 0, 100);
      // fMultCutLow->SetParameters(1893.94, -53.86, 0.502913, -0.0015122, 109.625, -1.19253);
      // fMultCutHigh = new TF1("fMultCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 3.*([4]+[5]*x)", 0, 100);
      // fMultCutHigh->SetParameters(1893.94, -53.86, 0.502913, -0.0015122, 109.625, -1.19253);
      // fMultMultPVCut = new TF1("fMultMultPVCut", "[0]+[1]*x+[2]*x*x", 0, 5000);
      // fMultMultPVCut->SetParameters(-0.1, 0.785, -4.7e-05);
    }
  }

  template <typename Collision>
  bool eventselection(Collision const& collision, float& multiplicity)
  {
    hvzero.fill(HIST("heventscheck"), 1.5);

    if (timFrameEvsel && (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
      return false;
    }
    hvzero.fill(HIST("heventscheck"), 2.5);

    if (!collision.sel8()) {
      return false;
    }
    hvzero.fill(HIST("heventscheck"), 3.5);

    if (piluprejection && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    hvzero.fill(HIST("heventscheck"), 4.5);

    if (goodzvertex && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    hvzero.fill(HIST("heventscheck"), 5.5);

    if (itstpctracks && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return false;
    }
    hvzero.fill(HIST("heventscheck"), 6.5);

    auto multNTracksPV = collision.multNTracksPV();
    if (additionalEvsel && multNTracksPV < fMultPVCutLow->Eval(multiplicity)) {
      return false;
    }
    hvzero.fill(HIST("heventscheck"), 7.5);
    if (additionalEvsel && multNTracksPV > fMultPVCutHigh->Eval(multiplicity)) {
      return false;
    }
    hvzero.fill(HIST("heventscheck"), 8.5);

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
    float lowmasscutks0 = 0.497 - cWidthKs0 * cSigmaMassKs0;
    float highmasscutks0 = 0.497 + cWidthKs0 * cSigmaMassKs0;
    // float decayLength = candidate.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * RecoDecay::sqrtSumOfSquares(candidate.px(), candidate.py(), candidate.pz());

    if (qAv0) {
      rKzeroShort.fill(HIST("hMassK0Shortbefore"), candidate.mK0Short(), candidate.pt());
      rKzeroShort.fill(HIST("hMasscorrelationbefore"), candidate.mK0Short(), candidate.mK0Short());
      rKzeroShort.fill(HIST("hLT"), ctauK0s);
      rKzeroShort.fill(HIST("hDCAV0Daughters"), candidate.dcaV0daughters());
      rKzeroShort.fill(HIST("hV0CosPA"), candidate.v0cosPA());
      rKzeroShort.fill(HIST("Mass_lambda"), candidate.mLambda());
      rKzeroShort.fill(HIST("mass_AntiLambda"), candidate.mAntiLambda());
      rKzeroShort.fill(HIST("mass_Gamma"), candidate.mGamma());
      rKzeroShort.fill(HIST("rapidity"), candidate.yK0Short());
      rKzeroShort.fill(HIST("hv0radius"), candidate.v0radius());
      rKzeroShort.fill(HIST("hDCApostopv"), candidate.dcapostopv());
      rKzeroShort.fill(HIST("hDCAnegtopv"), candidate.dcanegtopv());
      rKzeroShort.fill(HIST("hDCAv0topv"), candidate.dcav0topv());
      rKzeroShort.fill(HIST("halpha"), candidate.alpha());
      rKzeroShort.fill(HIST("hqtarmbyalpha"), arm);
      rKzeroShort.fill(HIST("hpsipair"), candidate.psipair());
    }
    if (correlation2Dhist)
      rKzeroShort.fill(HIST("mass_lambda_kshort_before"), candidate.mK0Short(), candidate.mLambda());

    hvzero.fill(HIST("htrackscheckv0"), 0.5);

    if (cDCAv0topv && std::fabs(candidate.dcav0topv()) > cMaxV0DCA) {
      return false;
    }
    hvzero.fill(HIST("htrackscheckv0"), 1.5);
    if (correlation2Dhist)
      rKzeroShort.fill(HIST("mass_lambda_kshort_after1"), candidate.mK0Short(), candidate.mLambda());

    if (rapidityks && std::abs(candidate.yK0Short()) >= confKsrapidity) {
      return false;
    }
    hvzero.fill(HIST("htrackscheckv0"), 2.5);
    if (correlation2Dhist)
      rKzeroShort.fill(HIST("mass_lambda_kshort_after2"), candidate.mK0Short(), candidate.mLambda());

    if (pT < confV0PtMin) {
      return false;
    }
    hvzero.fill(HIST("htrackscheckv0"), 3.5);
    if (correlation2Dhist)
      rKzeroShort.fill(HIST("mass_lambda_kshort_after3"), candidate.mK0Short(), candidate.mLambda());

    if (dcaDaughv0 > confV0DCADaughMax) {
      return false;
    }
    hvzero.fill(HIST("htrackscheckv0"), 4.5);
    if (correlation2Dhist)
      rKzeroShort.fill(HIST("mass_lambda_kshort_after4"), candidate.mK0Short(), candidate.mLambda());

    if (cpav0 < confV0CPAMin) {
      return false;
    }
    hvzero.fill(HIST("htrackscheckv0"), 5.5);
    if (correlation2Dhist)
      rKzeroShort.fill(HIST("mass_lambda_kshort_after5"), candidate.mK0Short(), candidate.mLambda());

    if (tranRad < confV0TranRadV0Min) {
      return false;
    }
    hvzero.fill(HIST("htrackscheckv0"), 6.5);
    if (correlation2Dhist)
      rKzeroShort.fill(HIST("mass_lambda_kshort_after6"), candidate.mK0Short(), candidate.mLambda());

    if (tranRad > confV0TranRadV0Max) {
      return false;
    }
    hvzero.fill(HIST("htrackscheckv0"), 7.5);
    if (correlation2Dhist)
      rKzeroShort.fill(HIST("mass_lambda_kshort_after7"), candidate.mK0Short(), candidate.mLambda());

    if (std::fabs(ctauK0s) > cMaxV0LifeTime) {
      return false;
    }
    hvzero.fill(HIST("htrackscheckv0"), 8.5);
    if (correlation2Dhist)
      rKzeroShort.fill(HIST("mass_lambda_kshort_after8"), candidate.mK0Short(), candidate.mLambda());

    if (armcut && arm < confarmcut) {
      return false;
    }
    hvzero.fill(HIST("htrackscheckv0"), 9.5);
    if (correlation2Dhist)
      rKzeroShort.fill(HIST("mass_lambda_kshort_after9"), candidate.mK0Short(), candidate.mLambda());

    if (applyCompetingcut && (std::abs(candidate.mLambda() - o2::constants::physics::MassLambda0) <= competingCascrejlambda || std::abs(candidate.mAntiLambda() - o2::constants::physics::MassLambda0) <= competingCascrejlambdaanti)) {
      return false;
    }

    hvzero.fill(HIST("htrackscheckv0"), 10.5);
    if (correlation2Dhist)
      rKzeroShort.fill(HIST("mass_lambda_kshort_after10"), candidate.mK0Short(), candidate.mLambda());

    if (qAv0) {
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
    if (qAPID) {
      // Filling the PID of the V0 daughters in the region of the K0 peak.
      (charge == 1) ? rKzeroShort.fill(HIST("hNSigmaPosPionK0s_before"), track.tpcInnerParam(), track.tpcNSigmaPi()) : rKzeroShort.fill(HIST("hNSigmaNegPionK0s_before"), track.tpcInnerParam(), track.tpcNSigmaPi());
      rKzeroShort.fill(HIST("dE_by_dx_TPC"), track.p(), track.tpcSignal());
    }
    const auto eta = track.eta();
    const auto tpcNClsF = track.tpcNClsFound();
    const auto sign = track.sign();

    hvzero.fill(HIST("htrackscheckv0daughters"), 0.5);

    if (hasTPC && !track.hasTPC())
      return false;
    hvzero.fill(HIST("htrackscheckv0daughters"), 1.5);

    if (!globalTracks) {
      if (track.tpcNClsCrossedRows() < tpcCrossedrows)
        return false;
      hvzero.fill(HIST("htrackscheckv0daughters"), 2.5);

      if (track.tpcCrossedRowsOverFindableCls() < tpcCrossedrowsOverfcls)
        return false;
      hvzero.fill(HIST("htrackscheckv0daughters"), 3.5);

      if (tpcNClsF < confDaughTPCnclsMin) {
        return false;
      }
      hvzero.fill(HIST("htrackscheckv0daughters"), 4.5);
    } else {
      if (!track.isGlobalTrack())
        return false;
      hvzero.fill(HIST("htrackscheckv0daughters"), 4.5);
    }

    if (charge < 0 && sign > 0) {
      return false;
    }
    hvzero.fill(HIST("htrackscheckv0daughters"), 5.5);

    if (charge > 0 && sign < 0) {
      return false;
    }
    hvzero.fill(HIST("htrackscheckv0daughters"), 6.5);

    if (std::abs(eta) > confDaughEta) {
      return false;
    }
    hvzero.fill(HIST("htrackscheckv0daughters"), 7.5);

    if (std::abs(nsigmaV0Daughter) > confDaughPIDCuts) {
      return false;
    }
    hvzero.fill(HIST("htrackscheckv0daughters"), 8.5);

    return true;
  }

  double massK0s = o2::constants::physics::MassK0Short;
  double massLambda = o2::constants::physics::MassLambda;
  double massPr = o2::constants::physics::MassProton;
  double massPi = o2::constants::physics::MassPionCharged;

  // Defining filters for events (event selection)
  // Filter eventFilter = (o2::aod::evsel::sel8 == true);

  Filter posZFilter = (nabs(o2::aod::collision::posZ) < cutzvertex);
  Filter acceptenceFilter = (nabs(aod::track::eta) < cfgETAcut && nabs(aod::track::pt) > cfgPTcut);

  // Filters on V0s
  Filter preFilterV0 = (nabs(aod::v0data::dcapostopv) > v0settingdcapostopv &&
                        nabs(aod::v0data::dcanegtopv) > v0settingdcanegtopv);

  // Defining the type of the daughter tracks
  using DaughterTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTOFFullPi>;
  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::MultZeqs, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::Mults>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullPr, aod::pidTOFFullPr>>;
  using V0TrackCandidate = aod::V0Datas;

  ROOT::Math::PxPyPzMVector daughter1, daughter2;
  ROOT::Math::PxPyPzMVector protonVec, pionVec, lambdaVec;
  ROOT::Math::PxPyPzMVector fourVecDau, fourVecMother, fourVecDauCM;
  ROOT::Math::XYZVector threeVecDauCM;
  ROOT::Math::XYZVector helicityVec, normalVec, randomVec, beamVec;

  void processSE(EventCandidates::iterator const& collision, TrackCandidates const& /*tracks*/, aod::V0Datas const& V0s)
  {
    hvzero.fill(HIST("heventscheck"), 0.5);

    float multiplicity = 0.0f;
    if (cfgMultFOTM) {
      multiplicity = collision.centFT0M();
    } else {
      multiplicity = collision.centFT0C();
    }
    if (!eventselection(collision, multiplicity)) {
      return;
    }

    auto occupancyno = collision.trackOccupancyInTimeRange();
    if (applyOccupancyCut && occupancyno < occupancyCut) {
      return;
    }

    if (qAevents) {
      rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());
      rEventSelection.fill(HIST("hmultiplicity"), multiplicity);
      rEventSelection.fill(HIST("multdist_FT0M"), collision.multFT0M());
      rEventSelection.fill(HIST("multdist_FT0A"), collision.multFT0A());
      rEventSelection.fill(HIST("multdist_FT0C"), collision.multFT0C());
      rEventSelection.fill(HIST("hNcontributor"), collision.numContrib());
    }

    std::vector<int> v0indexes;
    TLorentzVector lv1, lv3, lvLambda, lv4, lv5;
    for (const auto& [v1, v2] : combinations(CombinationsStrictlyUpperIndexPolicy(V0s, V0s))) {
      hvzero.fill(HIST("heventscheck"), 2.5);
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

      // selection for kshort

      double nTPCSigmaPosPi = postrack1.tpcNSigmaPi();
      double nTPCSigmaNegPi = negtrack1.tpcNSigmaPi();

      if (!(isSelectedV0Daughter(negtrack1, -1, nTPCSigmaNegPi, v1) && isSelectedV0Daughter(postrack1, 1, nTPCSigmaPosPi, v1))) {
        continue;
      }

      // for lamba baryon selection

      auto postrack2 = v2.template posTrack_as<TrackCandidates>();
      auto negtrack2 = v2.template negTrack_as<TrackCandidates>();

      double nTPCSigmaPosPr = postrack2.tpcNSigmaPr();
      double nTPCSigmaNegPiTrk2 = negtrack2.tpcNSigmaPi();

      double nTPCSigmaNegPr = negtrack2.tpcNSigmaPr();
      double nTPCSigmaPosPiTrk2 = postrack2.tpcNSigmaPi();

      int lambdaTag = 0;
      int alambdaTag = 0;

      if (isSelectedV0Daughter(postrack2, 1, nTPCSigmaPosPr, v2) && isSelectedV0Daughter(negtrack2, -1, nTPCSigmaNegPiTrk2, v2)) {
        lambdaTag = 1;
      }
      if (isSelectedV0Daughter(postrack2, 1, nTPCSigmaPosPiTrk2, v2) && isSelectedV0Daughter(negtrack2, -1, nTPCSigmaNegPr, v2)) {
        alambdaTag = 1;
      }

      if (!lambdaTag && !alambdaTag)
        continue;

      if (lambdaTag) {
        protonVec = ROOT::Math::PxPyPzMVector(postrack2.px(), postrack2.py(), postrack2.pz(), massPr);
        pionVec = ROOT::Math::PxPyPzMVector(negtrack2.px(), negtrack2.py(), negtrack2.pz(), massPi);
      }
      if (alambdaTag) {
        protonVec = ROOT::Math::PxPyPzMVector(negtrack2.px(), negtrack2.py(), negtrack2.pz(), massPr);
        pionVec = ROOT::Math::PxPyPzMVector(postrack2.px(), postrack2.py(), postrack2.pz(), massPi);
      }
      lambdaVec = protonVec + pionVec;

      lv1.SetPtEtaPhiM(v1.pt(), v1.eta(), v1.phi(), massK0s);
      lvLambda.SetPtEtaPhiM(lambdaVec.pt(), lambdaVec.eta(), lambdaVec.phi(), massLambda);

      lv3 = lv1 + lvLambda;

      if (qAv0daughters) {
        rKzeroShort.fill(HIST("negative_pt"), negtrack1.pt());
        rKzeroShort.fill(HIST("positive_pt"), postrack1.pt());
        rKzeroShort.fill(HIST("negative_eta"), negtrack1.eta());
        rKzeroShort.fill(HIST("positive_eta"), postrack1.eta());
        rKzeroShort.fill(HIST("negative_phi"), negtrack1.phi());
        rKzeroShort.fill(HIST("positive_phi"), postrack1.phi());
      }

      daughter1 = ROOT::Math::PxPyPzMVector(v1.px(), v1.py(), v1.pz(), massK0s);    // V01
      daughter2 = ROOT::Math::PxPyPzMVector(v2.px(), v2.py(), v2.pz(), massLambda); // V02

      // polarization calculations
      fourVecDau = ROOT::Math::PxPyPzMVector(daughter1.Px(), daughter1.Py(), daughter1.Pz(), massK0s); // Kshort
      fourVecMother = ROOT::Math::PxPyPzMVector(lv3.Px(), lv3.Py(), lv3.Pz(), lv3.M());                // mass of KshortKshort pair
      ROOT::Math::Boost boost{fourVecMother.BoostToCM()};                                              // boost mother to center of mass frame
      fourVecDauCM = boost(fourVecDau);                                                                // boost the frame of daughter same as mother
      threeVecDauCM = fourVecDauCM.Vect();                                                             // get the 3 vector of daughter in the frame of mother

      if (std::abs(lv3.Rapidity()) < 0.5) {
        if (invMass1D) {
          hvzero.fill(HIST("h1vzeropairInvMassRot"), lv3.M());
        }

        if (activateTHnSparseCosThStarHelicity) {
          helicityVec = fourVecMother.Vect(); // 3 vector of mother in COM frame
          auto cosThetaStarHelicity = helicityVec.Dot(threeVecDauCM) / (std::sqrt(threeVecDauCM.Mag2()) * std::sqrt(helicityVec.Mag2()));
          hvzero.fill(HIST("h3vzeropairInvMassDS"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarHelicity, occupancyno);

          for (int i = 0; i < cnofrotations; i++) {
            float theta2 = rn->Uniform(o2::constants::math::PI - o2::constants::math::PI / rotationalCut, o2::constants::math::PI + o2::constants::math::PI / rotationalCut);
            lv4.SetPtEtaPhiM(lambdaVec.pt(), lambdaVec.eta(), lambdaVec.phi() + theta2, massLambda); // for rotated background
            lv5 = lv1 + lv4;
            hvzero.fill(HIST("h3vzeropairInvMassRot"), multiplicity, lv5.Pt(), lv5.M(), cosThetaStarHelicity, occupancyno);
          }

        } else if (activateTHnSparseCosThStarProduction) {
          normalVec = ROOT::Math::XYZVector(lv3.Py(), -lv3.Px(), 0.f);
          auto cosThetaStarProduction = normalVec.Dot(threeVecDauCM) / (std::sqrt(threeVecDauCM.Mag2()) * std::sqrt(normalVec.Mag2()));
          hvzero.fill(HIST("h3vzeropairInvMassDS"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarProduction, occupancyno);
          for (int i = 0; i < cnofrotations; i++) {

            auto theta2 = rn->Uniform(o2::constants::math::PI - o2::constants::math::PI / rotationalCut, o2::constants::math::PI + o2::constants::math::PI / rotationalCut);

            lv4.SetPtEtaPhiM(lambdaVec.pt(), lambdaVec.eta(), lambdaVec.phi() + theta2, massLambda); // for rotated background
            lv5 = lv1 + lv4;
            hvzero.fill(HIST("h3vzeropairInvMassRot"), multiplicity, lv5.Pt(), lv5.M(), cosThetaStarProduction, occupancyno);
          }
        } else if (activateTHnSparseCosThStarBeam) {
          beamVec = ROOT::Math::XYZVector(0.f, 0.f, 1.f);
          auto cosThetaStarBeam = beamVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());
          hvzero.fill(HIST("h3vzeropairInvMassDS"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarBeam, occupancyno);
          for (int i = 0; i < cnofrotations; i++) {
            auto theta2 = rn->Uniform(o2::constants::math::PI - o2::constants::math::PI / rotationalCut, o2::constants::math::PI + o2::constants::math::PI / rotationalCut);
            lv4.SetPtEtaPhiM(lambdaVec.pt(), lambdaVec.eta(), lambdaVec.phi() + theta2, massLambda); // for rotated background
            lv5 = lv1 + lv4;
            hvzero.fill(HIST("h3vzeropairInvMassRot"), multiplicity, lv5.Pt(), lv5.M(), cosThetaStarBeam, occupancyno);
          }
        } else if (activateTHnSparseCosThStarRandom) {
          auto phiRandom = gRandom->Uniform(0.f, constants::math::TwoPI);
          auto thetaRandom = gRandom->Uniform(0.f, constants::math::PI);
          randomVec = ROOT::Math::XYZVector(std::sin(thetaRandom) * std::cos(phiRandom), std::sin(thetaRandom) * std::sin(phiRandom), std::cos(thetaRandom));
          auto cosThetaStarRandom = randomVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());
          hvzero.fill(HIST("h3vzeropairInvMassDS"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarRandom, occupancyno);

          for (int i = 0; i < cnofrotations; i++) {
            auto theta2 = rn->Uniform(o2::constants::math::PI - o2::constants::math::PI / rotationalCut, o2::constants::math::PI + o2::constants::math::PI / rotationalCut);
            lv4.SetPtEtaPhiM(lambdaVec.pt(), lambdaVec.eta(), lambdaVec.phi() + theta2, massLambda); // for rotated background
            lv5 = lv1 + lv4;
            hvzero.fill(HIST("h3vzeropairInvMassRot"), multiplicity, lv5.Pt(), lv5.M(), cosThetaStarRandom, occupancyno);
          }
        }
      }
    }
    if (qAv0) {
      int sizeofv0indexes = v0indexes.size();
      rKzeroShort.fill(HIST("NksProduced"), sizeofv0indexes);
      // std::cout << "Size of v0indexes: " << sizeofv0indexes << std::endl;
    }
  }
  PROCESS_SWITCH(Kshortlambda, processSE, "same event process", true);

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

    SameKindPair<EventCandidates, V0TrackCandidate, BinningTypeVertexContributor> pair1{binningOnPositions1, cfgNmixedEvents, -1, collisions, tracksTuple, &cache}; // for PbPb
    SameKindPair<EventCandidates, V0TrackCandidate, BinningTypeCentralityM> pair2{binningOnPositions2, cfgNmixedEvents, -1, collisions, tracksTuple, &cache};       // for pp

    if (cfgMultFOTM) {
      TLorentzVector lv1, lv3, lvLambda, lv4, lv5;

      for (const auto& [c1, tracks1, c2, tracks2] : pair2) // two different centrality c1 and c2 and tracks corresponding to them
      {

        float multiplicity = 0.0f;

        multiplicity = c1.centFT0M();

        if (!eventselection(c1, multiplicity) || !eventselection(c2, multiplicity)) {
          continue;
        }
        auto occupancyno = c1.trackOccupancyInTimeRange();
        auto occupancyno2 = c2.trackOccupancyInTimeRange();
        if (applyOccupancyCut && (occupancyno < occupancyCut || occupancyno2 < occupancyCut)) {
          continue;
        }

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

          // selection for kshort

          double nTPCSigmaPosPi = postrack1.tpcNSigmaPi();
          double nTPCSigmaNegPi = negtrack1.tpcNSigmaPi();

          if (!(isSelectedV0Daughter(negtrack1, -1, nTPCSigmaNegPi, t1) && isSelectedV0Daughter(postrack1, 1, nTPCSigmaPosPi, t1))) {
            continue;
          }

          // for lamba baryon selection

          auto postrack2 = t2.template posTrack_as<TrackCandidates>();
          auto negtrack2 = t2.template negTrack_as<TrackCandidates>();

          double nTPCSigmaPosPr = postrack2.tpcNSigmaPr();
          double nTPCSigmaNegPiTrk2 = negtrack2.tpcNSigmaPi();

          double nTPCSigmaNegPr = negtrack2.tpcNSigmaPr();
          double nTPCSigmaPosPiTrk2 = postrack2.tpcNSigmaPi();

          int lambdaTag = 0;
          int alambdaTag = 0;

          if (isSelectedV0Daughter(postrack2, 1, nTPCSigmaPosPr, t2) && isSelectedV0Daughter(negtrack2, -1, nTPCSigmaNegPiTrk2, t2)) {
            lambdaTag = 1;
          }
          if (isSelectedV0Daughter(postrack2, 1, nTPCSigmaPosPiTrk2, t2) && isSelectedV0Daughter(negtrack2, -1, nTPCSigmaNegPr, t2)) {
            alambdaTag = 1;
          }

          if (!lambdaTag && !alambdaTag)
            continue;

          if (lambdaTag) {
            protonVec = ROOT::Math::PxPyPzMVector(postrack2.px(), postrack2.py(), postrack2.pz(), massPr);
            pionVec = ROOT::Math::PxPyPzMVector(negtrack2.px(), negtrack2.py(), negtrack2.pz(), massPi);
          }
          if (alambdaTag) {
            protonVec = ROOT::Math::PxPyPzMVector(negtrack2.px(), negtrack2.py(), negtrack2.pz(), massPr);
            pionVec = ROOT::Math::PxPyPzMVector(postrack2.px(), postrack2.py(), postrack2.pz(), massPi);
          }
          lambdaVec = protonVec + pionVec;
          lv1.SetPtEtaPhiM(t1.pt(), t1.eta(), t1.phi(), massK0s);
          lvLambda.SetPtEtaPhiM(lambdaVec.pt(), lambdaVec.eta(), lambdaVec.phi(), massLambda);
          lv3 = lv1 + lvLambda;
          fourVecDau = ROOT::Math::PxPyPzMVector(daughter1.Px(), daughter1.Py(), daughter1.Pz(), massK0s); // Kshort
          fourVecMother = ROOT::Math::PxPyPzMVector(lv3.Px(), lv3.Py(), lv3.Pz(), lv3.M());                // mass of KshortKshort pair
          ROOT::Math::Boost boost{fourVecMother.BoostToCM()};                                              // boost mother to center of mass frame
          fourVecDauCM = boost(fourVecDau);                                                                // boost the frame of daughter same as mother
          threeVecDauCM = fourVecDauCM.Vect();                                                             // get the 3 vector of daughter in the frame of mother

          if (std::abs(lv3.Rapidity()) < 0.5) {
            if (activateTHnSparseCosThStarHelicity) {
              helicityVec = fourVecMother.Vect(); // 3 vector of mother in COM frame
              auto cosThetaStarHelicity = helicityVec.Dot(threeVecDauCM) / (std::sqrt(threeVecDauCM.Mag2()) * std::sqrt(helicityVec.Mag2()));
              hvzero.fill(HIST("h3vzeropairInvMassME"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarHelicity, occupancyno);
            } else if (activateTHnSparseCosThStarProduction) {
              normalVec = ROOT::Math::XYZVector(lv3.Py(), -lv3.Px(), 0.f);
              auto cosThetaStarProduction = normalVec.Dot(threeVecDauCM) / (std::sqrt(threeVecDauCM.Mag2()) * std::sqrt(normalVec.Mag2()));
              hvzero.fill(HIST("h3vzeropairInvMassME"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarProduction, occupancyno);
            } else if (activateTHnSparseCosThStarBeam) {
              beamVec = ROOT::Math::XYZVector(0.f, 0.f, 1.f);
              auto cosThetaStarBeam = beamVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());
              hvzero.fill(HIST("h3vzeropairInvMassME"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarBeam, occupancyno);
            } else if (activateTHnSparseCosThStarRandom) {
              auto phiRandom = gRandom->Uniform(0.f, constants::math::TwoPI);
              auto thetaRandom = gRandom->Uniform(0.f, constants::math::PI);
              randomVec = ROOT::Math::XYZVector(std::sin(thetaRandom) * std::cos(phiRandom), std::sin(thetaRandom) * std::sin(phiRandom), std::cos(thetaRandom));
              auto cosThetaStarRandom = randomVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());
              hvzero.fill(HIST("h3vzeropairInvMassME"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarRandom, occupancyno);
            }
          }
        }
      }
    } else {
      for (const auto& [c1, tracks1, c2, tracks2] : pair1) // two different centrality c1 and c2 and tracks corresponding to them
      {
        float multiplicity = 0.0f;
        multiplicity = c1.centFT0C();

        if (!eventselection(c1, multiplicity) || !eventselection(c2, multiplicity)) {
          continue;
        }
        auto occupancyno = c1.trackOccupancyInTimeRange();
        auto occupancyno2 = c2.trackOccupancyInTimeRange();

        if (applyOccupancyCut && (occupancyno < occupancyCut || occupancyno2 < occupancyCut)) {
          continue;
        }
        TLorentzVector lv1, lv3, lvLambda, lv4, lv5;
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

          // selection for kshort

          double nTPCSigmaPosPi = postrack1.tpcNSigmaPi();
          double nTPCSigmaNegPi = negtrack1.tpcNSigmaPi();

          if (!(isSelectedV0Daughter(negtrack1, -1, nTPCSigmaNegPi, t1) && isSelectedV0Daughter(postrack1, 1, nTPCSigmaPosPi, t1))) {
            continue;
          }

          // for lamba baryon selection

          auto postrack2 = t2.template posTrack_as<TrackCandidates>();
          auto negtrack2 = t2.template negTrack_as<TrackCandidates>();

          double nTPCSigmaPosPr = postrack2.tpcNSigmaPr();
          double nTPCSigmaNegPiTrk2 = negtrack2.tpcNSigmaPi();

          double nTPCSigmaNegPr = negtrack2.tpcNSigmaPr();
          double nTPCSigmaPosPiTrk2 = postrack2.tpcNSigmaPi();

          int lambdaTag = 0;
          int alambdaTag = 0;

          if (isSelectedV0Daughter(postrack2, 1, nTPCSigmaPosPr, t2) && isSelectedV0Daughter(negtrack2, -1, nTPCSigmaNegPiTrk2, t2)) {
            lambdaTag = 1;
          }
          if (isSelectedV0Daughter(postrack2, 1, nTPCSigmaPosPiTrk2, t2) && isSelectedV0Daughter(negtrack2, -1, nTPCSigmaNegPr, t2)) {
            alambdaTag = 1;
          }

          if (!lambdaTag && !alambdaTag)
            continue;

          if (lambdaTag) {
            protonVec = ROOT::Math::PxPyPzMVector(postrack2.px(), postrack2.py(), postrack2.pz(), massPr);
            pionVec = ROOT::Math::PxPyPzMVector(negtrack2.px(), negtrack2.py(), negtrack2.pz(), massPi);
          }
          if (alambdaTag) {
            protonVec = ROOT::Math::PxPyPzMVector(negtrack2.px(), negtrack2.py(), negtrack2.pz(), massPr);
            pionVec = ROOT::Math::PxPyPzMVector(postrack2.px(), postrack2.py(), postrack2.pz(), massPi);
          }
          lambdaVec = protonVec + pionVec;

          lv1.SetPtEtaPhiM(t1.pt(), t1.eta(), t1.phi(), massK0s);
          lvLambda.SetPtEtaPhiM(lambdaVec.pt(), lambdaVec.eta(), lambdaVec.phi(), massLambda);
          lv3 = lv1 + lvLambda;

          fourVecDau = ROOT::Math::PxPyPzMVector(daughter1.Px(), daughter1.Py(), daughter1.Pz(), massK0s); // Kshort
          fourVecMother = ROOT::Math::PxPyPzMVector(lv3.Px(), lv3.Py(), lv3.Pz(), lv3.M());                // mass of KshortKshort pair
          ROOT::Math::Boost boost{fourVecMother.BoostToCM()};                                              // boost mother to center of mass frame
          fourVecDauCM = boost(fourVecDau);                                                                // boost the frame of daughter same as mother
          threeVecDauCM = fourVecDauCM.Vect();                                                             // get the 3 vector of daughter in the frame of mother

          if (std::abs(lv3.Rapidity()) < 0.5) {

            if (activateTHnSparseCosThStarHelicity) {
              helicityVec = fourVecMother.Vect(); // 3 vector of mother in COM frame
              auto cosThetaStarHelicity = helicityVec.Dot(threeVecDauCM) / (std::sqrt(threeVecDauCM.Mag2()) * std::sqrt(helicityVec.Mag2()));
              hvzero.fill(HIST("h3vzeropairInvMassME"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarHelicity, occupancyno);
            } else if (activateTHnSparseCosThStarProduction) {
              normalVec = ROOT::Math::XYZVector(lv3.Py(), -lv3.Px(), 0.f);
              auto cosThetaStarProduction = normalVec.Dot(threeVecDauCM) / (std::sqrt(threeVecDauCM.Mag2()) * std::sqrt(normalVec.Mag2()));
              hvzero.fill(HIST("h3vzeropairInvMassME"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarProduction, occupancyno);
            } else if (activateTHnSparseCosThStarBeam) {
              beamVec = ROOT::Math::XYZVector(0.f, 0.f, 1.f);
              auto cosThetaStarBeam = beamVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());
              hvzero.fill(HIST("h3vzeropairInvMassME"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarBeam, occupancyno);
            } else if (activateTHnSparseCosThStarRandom) {
              auto phiRandom = gRandom->Uniform(0.f, constants::math::TwoPI);
              auto thetaRandom = gRandom->Uniform(0.f, constants::math::PI);
              randomVec = ROOT::Math::XYZVector(std::sin(thetaRandom) * std::cos(phiRandom), std::sin(thetaRandom) * std::sin(phiRandom), std::cos(thetaRandom));
              auto cosThetaStarRandom = randomVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());
              hvzero.fill(HIST("h3vzeropairInvMassME"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarRandom, occupancyno);
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(Kshortlambda, processME, "mixed event process", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<Kshortlambda>(cfgc)};
}
