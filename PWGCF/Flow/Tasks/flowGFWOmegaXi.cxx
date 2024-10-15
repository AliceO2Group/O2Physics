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

/// In case of questions please write to:
/// \author Fuchun Cui(fcui@cern.ch)

#include <CCDB/BasicCCDBManager.h>
#include <cmath>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/HistogramRegistry.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include "GFWPowerArray.h"
#include "GFW.h"
#include "GFWCumulant.h"
#include "GFWWeights.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/Core/EventPlaneHelper.h"
#include "ReconstructionDataFormats/Track.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "TList.h"
#include <TProfile.h>
#include <TRandom3.h>
#include <TF1.h>
#include <TF2.h>
#include <TPDGCode.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct FlowGFWOmegaXi {

  O2_DEFINE_CONFIGURABLE(cfgCutVertex, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMin, float, 0.2f, "Minimal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMax, float, 10.0f, "Maximal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMin, float, 0.2f, "Minimal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMax, float, 3.0f, "Maximal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutChi2prTPCcls, float, 2.5, "Chi2 per TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgUseNch, bool, false, "Use Nch for flow observables")
  O2_DEFINE_CONFIGURABLE(cfgNbootstrap, int, 10, "Number of subsamples")
  O2_DEFINE_CONFIGURABLE(cfgEfficiency, std::string, "", "CCDB path to efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgAcceptance, std::string, "", "CCDB path to acceptance object")
  O2_DEFINE_CONFIGURABLE(cfgOmegaMassbins, int, 16, "Number of Omega mass axis bins for c22")
  O2_DEFINE_CONFIGURABLE(cfgXiMassbins, int, 14, "Number of Xi mass axis bins for c22")
  O2_DEFINE_CONFIGURABLE(cfgK0sMassbins, int, 80, "Number of K0s mass axis bins for c22")
  O2_DEFINE_CONFIGURABLE(cfgLambdaMassbins, int, 32, "Number of Lambda mass axis bins for c22")
  // topological cut for V0
  O2_DEFINE_CONFIGURABLE(cfgv0_radius, float, 5.0f, "minimum decay radius")
  O2_DEFINE_CONFIGURABLE(cfgv0_v0cospa, float, 0.995f, "minimum cosine of pointing angle")
  O2_DEFINE_CONFIGURABLE(cfgv0_dcav0topv, float, 0.1f, "minimum daughter DCA to PV")
  O2_DEFINE_CONFIGURABLE(cfgv0_dcav0dau, float, 0.5f, "maximum DCA among V0 daughters")
  O2_DEFINE_CONFIGURABLE(cfgv0_mk0swindow, float, 0.1f, "Invariant mass window of K0s")
  O2_DEFINE_CONFIGURABLE(cfgv0_mlambdawindow, float, 0.04f, "Invariant mass window of lambda")
  O2_DEFINE_CONFIGURABLE(cfgv0_ArmPodocut, float, 0.2f, "Armenteros Podolski cut for K0")
  // topological cut for cascade
  O2_DEFINE_CONFIGURABLE(cfgcasc_radius, float, 0.5f, "minimum decay radius")
  O2_DEFINE_CONFIGURABLE(cfgcasc_casccospa, float, 0.999f, "minimum cosine of pointing angle")
  O2_DEFINE_CONFIGURABLE(cfgcasc_v0cospa, float, 0.998f, "minimum cosine of pointing angle")
  O2_DEFINE_CONFIGURABLE(cfgcasc_dcav0topv, float, 0.01f, "minimum daughter DCA to PV")
  O2_DEFINE_CONFIGURABLE(cfgcasc_dcabachtopv, float, 0.01f, "minimum bachelor DCA to PV")
  O2_DEFINE_CONFIGURABLE(cfgcasc_dcacascdau, float, 0.3f, "maximum DCA among cascade daughters")
  O2_DEFINE_CONFIGURABLE(cfgcasc_dcav0dau, float, 1.0f, "maximum DCA among V0 daughters")
  O2_DEFINE_CONFIGURABLE(cfgcasc_mlambdawindow, float, 0.04f, "Invariant mass window of lambda")
  // track quality and type selections
  O2_DEFINE_CONFIGURABLE(cfgtpcclusters, int, 70, "minimum number of TPC clusters requirement")
  O2_DEFINE_CONFIGURABLE(cfgitsclusters, int, 1, "minimum number of ITS clusters requirement")

  ConfigurableAxis cfgaxisVertex{"axisVertex", {20, -10, 10}, "vertex axis for histograms"};
  ConfigurableAxis cfgaxisPhi{"axisPhi", {60, 0.0, constants::math::TwoPI}, "phi axis for histograms"};
  ConfigurableAxis cfgaxisEta{"axisEta", {40, -1., 1.}, "eta axis for histograms"};
  ConfigurableAxis cfgaxisPtXi{"axisPtXi", {VARIABLE_WIDTH, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.9, 4.9, 5.9, 9.9}, "pt (GeV)"};
  ConfigurableAxis cfgaxisPtV0{"axisPtV0", {VARIABLE_WIDTH, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.9, 4.9, 5.9, 9.9}, "pt (GeV)"};

  Configurable<std::vector<float>> cfgvecwacc{"vecwacc", std::vector<float>{0.879543, 0.893808, 0.993375, 1.09663, 0.983883, 0.984094, 1.11362, 0.963896, 0.911212, 1.02934, 1.00295, 0.950711, 0.996856, 1.11934, 0.993665, 0.99087, 1.11915, 1.0198, 0.966849, 1.03237, 0.989367, 0.948312, 0.970883, 0.984305, 0.920335, 0.929722, 1.07467, 1.00862, 0.977185, 0.870868, 1.06552, 0.962393, 1.01025, 1.09959, 0.984226, 0.986361, 1.0931, 0.994377, 0.976051, 1.05249, 0.995538, 0.886452, 0.936763, 0.993613, 0.94491, 0.966559, 1.10829, 1.01998, 0.991503, 1.07918, 1.05655, 0.973784, 1.00914, 1.11678, 1.00092, 0.95232, 1.09814, 1.02322, 0.958543, 0.947231}, "wacc in phi bins"};
  AxisSpec axisPt{{0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.50, 4.00, 4.50, 5.00, 5.50, 6.00, 7.00, 8.00, 9.00, 10.0}, "pt(GeV)"};
  AxisSpec axisMultiplicity{{0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90}, "Centrality (%)"};
  AxisSpec axisOmegaminusMass = {80, 1.63f, 1.71f, "Inv. Mass (GeV)"};
  AxisSpec axisXiminusMass = {70, 1.3f, 1.37f, "Inv. Mass (GeV)"};
  AxisSpec axisK0sMass = {400, 0.4f, 0.6f, "Inv. Mass (GeV)"};
  AxisSpec axisLambdaMass = {160, 1.08f, 1.16f, "Inv. Mass (GeV)"};
  AxisSpec axisOmegaminusMassforflow = {cfgOmegaMassbins, 1.63f, 1.71f, "Inv. Mass (GeV)"};
  AxisSpec axisXiminusMassforflow = {cfgXiMassbins, 1.3f, 1.37f, "Inv. Mass (GeV)"};
  AxisSpec axisK0sMassforflow = {cfgK0sMassbins, 0.4f, 0.6f, "Inv. Mass (GeV)"};
  AxisSpec axisLambdaMassforflow = {cfgLambdaMassbins, 1.08f, 1.16f, "Inv. Mass (GeV)"};

  Configurable<bool> cfgcheckDauTPC{"checkDauTPC", false, "check if daughter tracks have TPC match"};
  Configurable<float> cfgCasc_rapidity{"Casc_rapidity", 0.5, "rapidity"};
  Configurable<float> cfgNSigmaCascPion{"NSigmaCascPion", 3, "NSigmaCascPion"};
  Configurable<float> cfgNSigmaCascProton{"NSigmaCascProton", 3, "NSigmaCascProton"};
  Configurable<float> cfgNSigmaCascKaon{"NSigmaCascKaon", 3, "NSigmaCascKaon"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtPOIMin) && (aod::track::pt < cfgCutPtPOIMax) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && (aod::track::tpcChi2NCl < cfgCutChi2prTPCcls);

  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<int64_t> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> url{"ccdb-url", "http://ccdb-test.cern.ch:8080", "url of the ccdb repository"};

  // Define output
  HistogramRegistry registry{"registry"};

  // define global variables
  GFW* fGFW = new GFW(); // GFW class used from main src
  std::vector<GFW::CorrConfig> corrconfigs;

  TH1D* mEfficiency = nullptr;
  GFWWeights* mAcceptance = nullptr;
  bool correctionsLoaded = false;

  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;
  TF1* fMultCutLow = nullptr;
  TF1* fMultCutHigh = nullptr;
  TF1* fT0AV0AMean = nullptr;
  TF1* fT0AV0ASigma = nullptr;

  EventPlaneHelper helperEP;
  std::vector<float> vecwa = cfgvecwacc;

  using TracksPID = soa::Join<aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr>;
  using aodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, TracksPID>>; // tracks filter
  using aodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::Mults>>;  // collisions filter
  using DaughterTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCPr, aod::pidTPCKa>;

  // Set the pt, mult and phi Axis;
  int nPtBins = axisPt.binEdges.size() - 1;
  double* PtBins = &(axisPt.binEdges)[0];
  TAxis* fPtAxis = new TAxis(nPtBins, PtBins);

  o2::framework::AxisSpec axisXiPt = cfgaxisPtXi;
  int nXiPtBins = axisXiPt.binEdges.size() - 1;
  double* XiPtBins = &(axisXiPt.binEdges)[0];
  TAxis* fXiPtAxis = new TAxis(nXiPtBins, XiPtBins);

  o2::framework::AxisSpec axisV0Pt = cfgaxisPtV0;
  int nV0PtBins = axisV0Pt.binEdges.size() - 1;
  double* V0PtBins = &(axisV0Pt.binEdges)[0];
  TAxis* fV0PtAxis = new TAxis(nV0PtBins, V0PtBins);

  o2::framework::AxisSpec axisMult = axisMultiplicity;
  int nMultBins = axisMult.binEdges.size() - 1;
  double* MultBins = &(axisMult.binEdges)[0];
  TAxis* fMultAxis = new TAxis(nMultBins, MultBins);

  int nPhiBins = 60;
  TAxis* fPhiAxis = new TAxis(nPhiBins, 0, constants::math::TwoPI);

  int nOmegaMassBins = cfgOmegaMassbins;
  TAxis* fOmegaMass = new TAxis(nOmegaMassBins, 1.63, 1.71);

  int nXiMassBins = cfgXiMassbins;
  TAxis* fXiMass = new TAxis(nXiMassBins, 1.3, 1.37);

  int nK0sMassBins = cfgK0sMassbins;
  TAxis* fK0sMass = new TAxis(nK0sMassBins, 0.4, 0.6);

  int nLambdaMassBins = cfgLambdaMassbins;
  TAxis* fLambdaMass = new TAxis(nLambdaMassBins, 1.08, 1.16);

  void init(InitContext const&) // Initialization
  {
    ccdb->setURL(url.value);
    ccdb->setCaching(true);
    ccdb->setCreatedNotAfter(nolaterthan.value);

    // Add some output objects to the histogram registry
    registry.add("hPhi", "", {HistType::kTH1D, {cfgaxisPhi}});
    registry.add("hEta", "", {HistType::kTH1D, {cfgaxisEta}});
    registry.add("hVtxZ", "", {HistType::kTH1D, {cfgaxisVertex}});
    registry.add("hMult", "", {HistType::kTH1D, {{3000, 0.5, 3000.5}}});
    registry.add("hCent", "", {HistType::kTH1D, {{90, 0, 90}}});
    registry.add("hPt", "", {HistType::kTH1D, {axisPt}});
    registry.add("hEtaPhiVtxzREF", "", {HistType::kTH3D, {cfgaxisEta, cfgaxisPhi, {20, -10, 10}}});
    registry.add("hEtaPhiVtxzPOIXi", "", {HistType::kTH3D, {cfgaxisEta, cfgaxisPhi, {20, -10, 10}}});
    registry.add("hEtaPhiVtxzPOIOmega", "", {HistType::kTH3D, {cfgaxisEta, cfgaxisPhi, {20, -10, 10}}});
    registry.add("hEtaPhiVtxzPOIK0s", "", {HistType::kTH3D, {cfgaxisEta, cfgaxisPhi, {20, -10, 10}}});
    registry.add("hEtaPhiVtxzPOILambda", "", {HistType::kTH3D, {cfgaxisEta, cfgaxisPhi, {20, -10, 10}}});
    // cumulant of flow
    registry.add("c22", ";Centrality  (%) ; C_{2}{2}", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("c24", ";Centrality  (%) ; C_{2}{4}", {HistType::kTProfile, {axisMultiplicity}});
    // pt-diff cumulant of flow
    registry.add("Xic22dpt", ";pt ; C_{2}{2} ", {HistType::kTProfile3D, {cfgaxisPtXi, axisXiminusMassforflow, axisMultiplicity}});
    registry.add("Omegac22dpt", ";pt ; C_{2}{2} ", {HistType::kTProfile3D, {cfgaxisPtXi, axisOmegaminusMassforflow, axisMultiplicity}});
    registry.add("K0sc22dpt", ";pt ; C_{2}{2} ", {HistType::kTProfile3D, {cfgaxisPtV0, axisK0sMassforflow, axisMultiplicity}});
    registry.add("Lambdac22dpt", ";pt ; C_{2}{2} ", {HistType::kTProfile3D, {cfgaxisPtV0, axisLambdaMassforflow, axisMultiplicity}});
    // InvMass(GeV) of casc and v0
    registry.add("InvMassXiMinus_all", "", {HistType::kTHnSparseF, {cfgaxisPtXi, axisXiminusMass, cfgaxisEta, axisMultiplicity}});
    registry.add("InvMassOmegaMinus_all", "", {HistType::kTHnSparseF, {cfgaxisPtXi, axisOmegaminusMass, cfgaxisEta, axisMultiplicity}});
    registry.add("InvMassOmegaMinus", "", {HistType::kTHnSparseF, {cfgaxisPtXi, axisOmegaminusMass, cfgaxisEta, axisMultiplicity}});
    registry.add("InvMassXiMinus", "", {HistType::kTHnSparseF, {cfgaxisPtXi, axisXiminusMass, cfgaxisEta, axisMultiplicity}});
    registry.add("InvMassK0s_all", "", {HistType::kTHnSparseF, {cfgaxisPtV0, axisK0sMass, cfgaxisEta, axisMultiplicity}});
    registry.add("InvMassLambda_all", "", {HistType::kTHnSparseF, {cfgaxisPtV0, axisLambdaMass, cfgaxisEta, axisMultiplicity}});
    registry.add("InvMassK0s", "", {HistType::kTHnSparseF, {cfgaxisPtV0, axisK0sMass, cfgaxisEta, axisMultiplicity}});
    registry.add("InvMassLambda", "", {HistType::kTHnSparseF, {cfgaxisPtV0, axisLambdaMass, cfgaxisEta, axisMultiplicity}});

    fGFW->AddRegion("full", -0.8, 0.8, 1, 1); // ("name", etamin, etamax, ptbinnum, bitmask)eta region -0.8 to 0.8
    // with (-0.5, 0.5) eta gap
    fGFW->AddRegion("refN10", -0.8, -0.5, 1, 1);
    fGFW->AddRegion("refP10", 0.5, 0.8, 1, 1);
    int nXiptMassBins = nXiPtBins * nXiMassBins;
    fGFW->AddRegion("poiXiP", 0.5, 0.8, nXiptMassBins, 2);
    int nOmegaptMassBins = nXiPtBins * nOmegaMassBins;
    fGFW->AddRegion("poiOmegaP", 0.5, 0.8, nOmegaptMassBins, 4);
    int nK0sptMassBins = nV0PtBins * nK0sMassBins;
    fGFW->AddRegion("poiK0sP", 0.5, 0.8, nK0sptMassBins, 8);
    int nLambdaptMassBins = nV0PtBins * nLambdaMassBins;
    fGFW->AddRegion("poiLambdaP", 0.5, 0.8, nLambdaptMassBins, 16);
    // pushback
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refP10 {2} refN10 {-2}", "ChFull220", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refP10 {2 2} refN10 {-2 -2}", "ChFull240", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiXiP {2} refN10 {-2}", "Ch10Gap22", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiOmegaP {2} refN10 {-2}", "Ch10Gap22", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiK0sP {2} refN10 {-2}", "Ch10Gap22", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiLambdaP {2} refN10 {-2}", "Ch10Gap22", kTRUE));
    //
    fGFW->CreateRegions(); // finalize the initialization

    // used for event selection
    fMultPVCutLow = new TF1("fMultPVCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x - 3.5*([5]+[6must ]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)", 0, 100);
    fMultPVCutLow->SetParameters(3257.29, -121.848, 1.98492, -0.0172128, 6.47528e-05, 154.756, -1.86072, -0.0274713, 0.000633499, -3.37757e-06);
    fMultPVCutHigh = new TF1("fMultPVCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x + 3.5*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)", 0, 100);
    fMultPVCutHigh->SetParameters(3257.29, -121.848, 1.98492, -0.0172128, 6.47528e-05, 154.756, -1.86072, -0.0274713, 0.000633499, -3.37757e-06);
    fMultCutLow = new TF1("fMultCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 2.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
    fMultCutLow->SetParameters(1654.46, -47.2379, 0.449833, -0.0014125, 150.773, -3.67334, 0.0530503, -0.000614061, 3.15956e-06);
    fMultCutHigh = new TF1("fMultCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 3.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
    fMultCutHigh->SetParameters(1654.46, -47.2379, 0.449833, -0.0014125, 150.773, -3.67334, 0.0530503, -0.000614061, 3.15956e-06);
    fT0AV0AMean = new TF1("fT0AV0AMean", "[0]+[1]*x", 0, 200000);
    fT0AV0AMean->SetParameters(-1601.0581, 9.417652e-01);
    fT0AV0ASigma = new TF1("fT0AV0ASigma", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 200000);
    fT0AV0ASigma->SetParameters(463.4144, 6.796509e-02, -9.097136e-07, 7.971088e-12, -2.600581e-17);
  }

  template <char... chars>
  void FillProfile(const GFW::CorrConfig& corrconf, const ConstStr<chars...>& tarName, const double& cent)
  {
    float dnx = 0;
    float val = 0;
    dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
    if (dnx == 0)
      return;
    if (!corrconf.pTDif) {
      val = fGFW->Calculate(corrconf, 0, kFALSE).real() / dnx;
      if (TMath::Abs(val) < 1) {
        registry.fill(tarName, cent, val, dnx);
      }
      return;
    }
    return;
  }

  template <char... chars>
  void FillProfilepT(const GFW::CorrConfig& corrconf, const ConstStr<chars...>& tarName, const int& ptbin, const int& PDGCode, const float& cent)
  {
    int nMassBins = 0;
    int nptbins = nPtBins;
    TAxis* fMass = nullptr;
    TAxis* fpt = nullptr;
    if (PDGCode == kXiMinus) {
      nMassBins = nXiMassBins;
      nptbins = nXiPtBins;
      fpt = fXiPtAxis;
      fMass = fXiMass;
    } else if (PDGCode == kOmegaMinus) {
      nMassBins = nOmegaMassBins;
      nptbins = nXiPtBins;
      fpt = fXiPtAxis;
      fMass = fOmegaMass;
    } else if (PDGCode == kK0Short) {
      nMassBins = nK0sMassBins;
      nptbins = nV0PtBins;
      fpt = fV0PtAxis;
      fMass = fK0sMass;
    } else if (PDGCode == kLambda0) {
      nMassBins = nLambdaMassBins;
      nptbins = nV0PtBins;
      fpt = fV0PtAxis;
      fMass = fLambdaMass;
    } else {
      LOGF(error, "Error, please put in correct PDGCode of K0s, Lambda, Xi or Omega");
      return;
    }
    for (int massbin = 1; massbin <= nMassBins; massbin++) {
      float dnx = 0;
      float val = 0;
      dnx = fGFW->Calculate(corrconf, (ptbin - 1) + ((massbin - 1) * nptbins), kTRUE).real();
      if (dnx == 0)
        continue;
      val = fGFW->Calculate(corrconf, (ptbin - 1) + ((massbin - 1) * nptbins), kFALSE).real() / dnx;
      if (TMath::Abs(val) < 1) {
        registry.fill(tarName, fpt->GetBinCenter(ptbin), fMass->GetBinCenter(massbin), cent, val, dnx);
      }
    }
    return;
  }

  void loadCorrections(uint64_t timestamp)
  {
    if (correctionsLoaded)
      return;
    if (cfgAcceptance.value.empty() == false) {
      mAcceptance = ccdb->getForTimeStamp<GFWWeights>(cfgAcceptance, timestamp);
      if (mAcceptance)
        LOGF(info, "Loaded acceptance weights from %s (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance);
      else
        LOGF(warning, "Could not load acceptance weights from %s (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance);
    }
    if (cfgEfficiency.value.empty() == false) {
      mEfficiency = ccdb->getForTimeStamp<TH1D>(cfgEfficiency, timestamp);
      if (mEfficiency == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", cfgEfficiency.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram from %s (%p)", cfgEfficiency.value.c_str(), (void*)mEfficiency);
    }
    correctionsLoaded = true;
  }

  template <typename TrackObject>
  bool setCurrentParticleWeights(float& weight_nue, float& weight_nua, TrackObject track, float vtxz)
  {
    float eff = 1.;
    if (mEfficiency)
      eff = mEfficiency->GetBinContent(mEfficiency->FindBin(track.pt()));
    else
      eff = 1.0;
    if (eff == 0)
      return false;
    weight_nue = 1. / eff;
    if (mAcceptance)
      weight_nua = mAcceptance->GetNUA(track.phi(), track.eta(), vtxz);
    else
      weight_nua = 1;
    return true;
  }
  // event selection
  template <typename TCollision>
  bool eventSelected(TCollision collision, const int multTrk, const float centrality)
  {
    if (collision.alias_bit(kTVXinTRD)) {
      // TRD triggered
      return false;
    }
    if (!collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      // reject collisions close to Time Frame borders
      // https://its.cern.ch/jira/browse/O2-4623
      return false;
    }
    if (!collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      // reject events affected by the ITS ROF border
      // https://its.cern.ch/jira/browse/O2-4309
      return false;
    }
    if (!collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      // rejects collisions which are associated with the same "found-by-T0" bunch crossing
      // https://indico.cern.ch/event/1396220/#1-event-selection-with-its-rof
      return false;
    }
    if (!collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      // removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference
      // use this cut at low multiplicities with caution
      return false;
    }
    float vtxz = -999;
    if (collision.numContrib() > 1) {
      vtxz = collision.posZ();
      float zRes = TMath::Sqrt(collision.covZZ());
      if (zRes > 0.25 && collision.numContrib() < 20)
        vtxz = -999;
    }
    auto multNTracksPV = collision.multNTracksPV();

    if (abs(vtxz) > cfgCutVertex)
      return false;
    if (multNTracksPV < fMultPVCutLow->Eval(centrality))
      return false;
    if (multNTracksPV > fMultPVCutHigh->Eval(centrality))
      return false;
    if (multTrk < fMultCutLow->Eval(centrality))
      return false;
    if (multTrk > fMultCutHigh->Eval(centrality))
      return false;

    // V0A T0A 5 sigma cut
    if (abs(collision.multFV0A() - fT0AV0AMean->Eval(collision.multFT0A())) > 5 * fT0AV0ASigma->Eval(collision.multFT0A()))
      return 0;

    return true;
  }

  void process(aodCollisions::iterator const& collision, aod::BCsWithTimestamps const&, aodTracks const& tracks, aod::CascDataExt const& Cascades, aod::V0Datas const& V0s, DaughterTracks&)
  {
    int Ntot = tracks.size();
    if (Ntot < 1)
      return;
    fGFW->Clear();
    const auto cent = collision.centFT0C();

    if (eventSelected(collision, tracks.size(), cent))
      return;
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    loadCorrections(bc.timestamp());
    float vtxz = collision.posZ();
    registry.fill(HIST("hVtxZ"), vtxz);
    registry.fill(HIST("hMult"), Ntot);
    registry.fill(HIST("hCent"), collision.centFT0C());

    float weff = 1;
    float wacc = 1;
    // fill GFW ref flow
    for (auto& track : tracks) {
      if (!setCurrentParticleWeights(weff, wacc, track, vtxz))
        continue;
      int phibin = fPhiAxis->FindBin(track.phi()) - 1;
      wacc = 1 / vecwa[phibin];
      registry.fill(HIST("hPhi"), track.phi());
      registry.fill(HIST("hEta"), track.eta());
      registry.fill(HIST("hEtaPhiVtxzREF"), track.eta(), track.phi(), vtxz);
      registry.fill(HIST("hPt"), track.pt());
      int ptbin = fPtAxis->FindBin(track.pt()) - 1;
      if ((track.pt() > cfgCutPtMin) && (track.pt() < cfgCutPtMax)) {
        fGFW->Fill(track.eta(), ptbin, track.phi(), wacc * weff, 1); //(eta, ptbin, phi, wacc*weff, bitmask)
      }
    }
    // fill GFW of V0 flow
    for (auto& v0 : V0s) {
      auto v0posdau = v0.posTrack_as<DaughterTracks>();
      auto v0negdau = v0.negTrack_as<DaughterTracks>();
      // check tpc
      int PDGCode = 0;
      if (v0.qtarm() / TMath::Abs(v0.alpha()) > cfgv0_ArmPodocut && TMath::Abs(v0posdau.tpcNSigmaPi()) < cfgNSigmaCascPion && TMath::Abs(v0negdau.tpcNSigmaPi()) < cfgNSigmaCascPion) {
        registry.fill(HIST("InvMassK0s_all"), v0.pt(), v0.mK0Short(), v0.eta(), cent);
        PDGCode = kK0Short;
      } else if (v0.qtarm() / TMath::Abs(v0.alpha()) < cfgv0_ArmPodocut && TMath::Abs(v0posdau.tpcNSigmaPr()) < cfgNSigmaCascProton && TMath::Abs(v0negdau.tpcNSigmaPi()) < cfgNSigmaCascPion) {
        registry.fill(HIST("InvMassLambda_all"), v0.pt(), v0.mLambda(), v0.eta(), cent);
        PDGCode = kLambda0;
      }
      // track quality check
      if (v0posdau.tpcNClsFound() < cfgtpcclusters)
        continue;
      if (v0negdau.tpcNClsFound() < cfgtpcclusters)
        continue;
      if (v0posdau.itsNCls() < cfgitsclusters)
        continue;
      if (v0negdau.itsNCls() < cfgitsclusters)
        continue;
      // topological cut
      if (v0.v0radius() < cfgv0_radius)
        continue;
      if (v0.v0cosPA() < cfgv0_v0cospa)
        continue;
      if (v0.dcaV0daughters() > cfgv0_dcav0dau)
        continue;
      if (PDGCode == kK0Short) {
        if (TMath::Abs(v0.mK0Short() - 0.49761) < cfgv0_mk0swindow) {
          registry.fill(HIST("InvMassK0s"), v0.pt(), v0.mK0Short(), v0.eta(), cent);
          registry.fill(HIST("hEtaPhiVtxzPOIK0s"), v0.eta(), v0.phi(), vtxz);
          fGFW->Fill(v0.eta(), fV0PtAxis->FindBin(v0.pt()) - 1 + ((fK0sMass->FindBin(v0.mK0Short()) - 1) * nV0PtBins), v0.phi(), wacc * weff, 8);
        }
      } else if (PDGCode == kLambda0) {
        if (TMath::Abs(v0.mLambda() - 1.115683) < cfgv0_mlambdawindow) {
          registry.fill(HIST("InvMassLambda"), v0.pt(), v0.mLambda(), v0.eta(), cent);
          registry.fill(HIST("hEtaPhiVtxzPOILambda"), v0.eta(), v0.phi(), vtxz);
          fGFW->Fill(v0.eta(), fV0PtAxis->FindBin(v0.pt()) - 1 + ((fLambdaMass->FindBin(v0.mLambda()) - 1) * nV0PtBins), v0.phi(), wacc * weff, 16);
        }
      }
    }
    // fill GFW of casc flow
    for (auto& casc : Cascades) {
      auto bachelor = casc.bachelor_as<DaughterTracks>();
      auto posdau = casc.posTrack_as<DaughterTracks>();
      auto negdau = casc.negTrack_as<DaughterTracks>();
      // check TPC
      if (cfgcheckDauTPC && (!posdau.hasTPC() || !negdau.hasTPC() || !bachelor.hasTPC())) {
        continue;
      }
      int PDGCode = 0;
      if (casc.sign() < 0 && TMath::Abs(casc.yOmega()) < cfgCasc_rapidity && TMath::Abs(bachelor.tpcNSigmaKa()) < cfgNSigmaCascKaon && TMath::Abs(posdau.tpcNSigmaPr()) < cfgNSigmaCascProton && TMath::Abs(negdau.tpcNSigmaPi()) < cfgNSigmaCascPion) {
        registry.fill(HIST("InvMassOmegaMinus_all"), casc.pt(), casc.mOmega(), casc.eta(), cent);
        PDGCode = kOmegaMinus;
      } else if (casc.sign() < 0 && TMath::Abs(casc.yXi()) < cfgCasc_rapidity && TMath::Abs(bachelor.tpcNSigmaPi()) < cfgNSigmaCascPion && TMath::Abs(posdau.tpcNSigmaPr()) < cfgNSigmaCascProton && TMath::Abs(negdau.tpcNSigmaPi()) < cfgNSigmaCascPion) {
        registry.fill(HIST("InvMassXiMinus_all"), casc.pt(), casc.mXi(), casc.eta(), cent);
        PDGCode = kXiMinus;
      } else {
        continue;
      }
      // topological cut
      if (casc.cascradius() < cfgcasc_radius)
        continue;
      if (casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) < cfgcasc_casccospa)
        continue;
      if (casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < cfgcasc_v0cospa)
        continue;
      if (casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()) < cfgcasc_dcav0topv)
        continue;
      if (casc.dcabachtopv() < cfgcasc_dcabachtopv)
        continue;
      if (casc.dcacascdaughters() > cfgcasc_dcacascdau)
        continue;
      if (casc.dcaV0daughters() > cfgcasc_dcav0dau)
        continue;
      if (TMath::Abs(casc.mLambda() - 1.115683) > cfgcasc_mlambdawindow)
        continue;
      // track quality check
      if (bachelor.tpcNClsFound() < cfgtpcclusters)
        continue;
      if (posdau.tpcNClsFound() < cfgtpcclusters)
        continue;
      if (negdau.tpcNClsFound() < cfgtpcclusters)
        continue;
      if (bachelor.itsNCls() < cfgitsclusters)
        continue;
      if (posdau.itsNCls() < cfgitsclusters)
        continue;
      if (negdau.itsNCls() < cfgitsclusters)
        continue;
      if (PDGCode == kOmegaMinus) {
        registry.fill(HIST("hEtaPhiVtxzPOIOmega"), casc.eta(), casc.phi(), vtxz);
        registry.fill(HIST("InvMassOmegaMinus"), casc.pt(), casc.mOmega(), casc.eta(), cent);
        if ((casc.pt() < cfgCutPtPOIMax) && (casc.pt() > cfgCutPtPOIMin) && (casc.mOmega() > 1.63) && (casc.mOmega() < 1.71)) {
          fGFW->Fill(casc.eta(), fXiPtAxis->FindBin(casc.pt()) - 1 + ((fOmegaMass->FindBin(casc.mOmega()) - 1) * nXiPtBins), casc.phi(), wacc * weff, 4);
        }
      }
      if (PDGCode == kXiMinus) {
        registry.fill(HIST("hEtaPhiVtxzPOIXi"), casc.eta(), casc.phi(), vtxz);
        registry.fill(HIST("InvMassXiMinus"), casc.pt(), casc.mXi(), casc.eta(), cent);
        if ((casc.pt() < cfgCutPtPOIMax) && (casc.pt() > cfgCutPtPOIMin) && (casc.mXi() > 1.30) && (casc.mXi() < 1.37)) {
          fGFW->Fill(casc.eta(), fXiPtAxis->FindBin(casc.pt()) - 1 + ((fXiMass->FindBin(casc.mXi()) - 1) * nXiPtBins), casc.phi(), wacc * weff, 2);
        }
      }
    }
    // Filling cumulant with ROOT TProfile
    FillProfile(corrconfigs.at(0), HIST("c22"), cent);
    FillProfile(corrconfigs.at(1), HIST("c24"), cent);
    for (int i = 1; i <= nV0PtBins; i++) // loop for all ptBins
    {
      FillProfilepT(corrconfigs.at(4), HIST("K0sc22dpt"), i, kK0Short, cent);
      FillProfilepT(corrconfigs.at(5), HIST("Lambdac22dpt"), i, kLambda0, cent);
    }
    for (int i = 1; i <= nXiPtBins; i++) {
      FillProfilepT(corrconfigs.at(2), HIST("Xic22dpt"), i, kXiMinus, cent);
      FillProfilepT(corrconfigs.at(3), HIST("Omegac22dpt"), i, kOmegaMinus, cent);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowGFWOmegaXi>(cfgc)};
}
