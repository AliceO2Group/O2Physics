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
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"

#include "TList.h"
#include <TProfile.h>
#include <TRandom3.h>
#include <TF1.h>
#include <TF2.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{ #NAME, DEFAULT, HELP };

struct FlowGFWPbPb {

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

  ConfigurableAxis axisVertex{ "axisVertex", { 20, -10, 10 }, "vertex axis for histograms" };
  ConfigurableAxis axisPhi{ "axisPhi", { 60, 0.0, constants::math::TwoPI }, "phi axis for histograms" };
  ConfigurableAxis axisEta{ "axisEta", { 40, -1., 1. }, "eta axis for histograms" };
  // ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.2, 0.25, 0.30, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00}, "pt axis for histograms"};
  // ConfigurableAxis axisPt{"axisPt", {0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00}, "pt axis for histograms"};
  // ConfigurableAxis axisPt{"axisPt", {0.01, 0.1, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.80, 2.90, 3.00}, "pt axis for histograms"};
  // ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90}, "centrality axis for histograms"};

  Configurable<std::vector<float>> vecCnMult{ "vecCnMult", std::vector<float>{ 0.00171171, 0.00520243, 0.00438017, 0.00351376, -0.000248092, 0.00133333, 0.00340237, -0.000133215, 0.0020018, -0.00188645, -0.042973, -0.0441498, -0.0446074, -0.0426147, -0.0485496, -0.0411961, -0.0380375, -0.0361012, -0.0277648, -0.0305861 }, "Cn of MultBin" };
  Configurable<std::vector<float>> vecSnMult{ "vecSnMult", std::vector<float>{ 0.00301802, -0.000668016, -0.00239669, -0.00277982, 0.000534351, -0.00280392, -0.00967456, -0.0032238, 0.00713645, 0.00344322, -0.000765766, -0.00346154, -0.00485537, -0.00180734, 0.000629771, 0.000235294, -0.00910256, -0.000275311, -0.00166068, 0.000421245 }, "Sn of MultBin" };
  Configurable<std::vector<float>> vecL2nspMult{ "vecL2nspMult", std::vector<float>{ 0.00463219, 0.00517462, 0.00562754, 0.00710129, 0.00762803, 0.00569557, 0.00730412, 0.00750648, 0.018422, 0.0134009, -0.0211638, -0.0207159, -0.0190476, -0.0220175, -0.0209654, -0.0184701, -0.0268248, -0.0284841, -0.0234963, -0.0357846 }, "lamda_2n^s+ of MultBin" };
  Configurable<std::vector<float>> vecL2nsmMult{ "vecL2nsmMult", std::vector<float>{ 0.00374913, 0.00426522, 0.00470312, 0.00600039, 0.00647403, 0.00481039, 0.00621549, 0.00641559, 0.0160814, 0.0114422, -0.0195323, -0.0190689, -0.0175432, -0.0202407, -0.0193653, -0.0170386, -0.0248667, -0.0268234, -0.0219932, -0.0328736 }, "lamda_2n^s- of MultBin" };
  Configurable<std::vector<float>> vecA2npMult{ "vecA2npMult", std::vector<float>{ 0.89464, 0.903664, 0.910517, 0.915972, 0.918168, 0.915745, 0.919477, 0.921643, 0.932163, 0.921154, 0.95991, 0.958603, 0.958884, 0.957954, 0.960324, 0.959686, 0.96212, 0.969973, 0.966957, 0.957601 }, "a_2n^+ of MultBin" };
  Configurable<std::vector<float>> vecA2nmMult{ "vecA2nmMult", std::vector<float>{ 1.10536, 1.09634, 1.08948, 1.08403, 1.08183, 1.08425, 1.08052, 1.07836, 1.06784, 1.07885, 1.04009, 1.0414, 1.04112, 1.04205, 1.03968, 1.04031, 1.03788, 1.03003, 1.03304, 1.0424 }, "a_2n^- of MultBin" };
  Configurable<std::vector<float>> vecwacc{ "vecwacc", std::vector<float>{ 0.879543, 0.893808, 0.993375, 1.09663, 0.983883, 0.984094, 1.11362, 0.963896, 0.911212, 1.02934, 1.00295, 0.950711, 0.996856, 1.11934, 0.993665, 0.99087, 1.11915, 1.0198, 0.966849, 1.03237, 0.989367, 0.948312, 0.970883, 0.984305, 0.920335, 0.929722, 1.07467, 1.00862, 0.977185, 0.870868, 1.06552, 0.962393, 1.01025, 1.09959, 0.984226, 0.986361, 1.0931, 0.994377, 0.976051, 1.05249, 0.995538, 0.886452, 0.936763, 0.993613, 0.94491, 0.966559, 1.10829, 1.01998, 0.991503, 1.07918, 1.05655, 0.973784, 1.00914, 1.11678, 1.00092, 0.95232, 1.09814, 1.02322, 0.958543, 0.947231 }, "wacc in phi bins" };
  AxisSpec axisPt{ { 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.50, 4.00, 4.50, 5.00, 5.50, 6.00, 7.00, 8.00, 9.00, 10.0 }, "pt(GeV)" };
  AxisSpec axisMultiplicity{ { 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90 }, "Centrality (%)" };
  AxisSpec axisOmegaminusMass{ { 1.63, 1.635, 1.64, 1.645, 1.65, 1.655, 1.66, 1.665, 1.67, 1.675, 1.68, 1.685, 1.69, 1.695, 1.7, 1.705, 1.71 }, "Mass(GeV)" };
  AxisSpec axisXiminusMass{ { 1.3, 1.305, 1.31, 1.315, 1.32, 1.325, 1.33, 1.335, 1.34, 1.345, 1.35, 1.355, 1.36, 1.365, 1.37 }, "Mass(GeV)" };

  Configurable<bool> checkDauTPC{ "checkDauTPC", false, "check if daughter tracks have TPC match" };
  Configurable<float> Casc_rapidity{ "Casc_rapidity", 0.5, "rapidity" };
  Configurable<float> NSigmaCascPion{ "NSigmaCascPion", 6, "NSigmaCascPion" };
  Configurable<float> NSigmaCascProton{ "NSigmaCascProton", 6, "NSigmaCascProton" };
  Configurable<float> NSigmaCascKaon{ "NSigmaCascKaon", 6, "NSigmaCascKaon" };

  // AxisSpec axisC22{{0,0.02}, "C_2{2}"};

  //  std::vector<double> centBinning = {0, 5., 10., 20., 30., 40., 50., 60., 70.};
  //  AxisSpec axisCentBins{centBinning, "centrality percentile"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtPOIMin) && (aod::track::pt < cfgCutPtPOIMax) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && (aod::track::tpcChi2NCl < cfgCutChi2prTPCcls);

  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<int64_t> nolaterthan{ "ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object" };
  Configurable<std::string> url{ "ccdb-url", "http://ccdb-test.cern.ch:8080", "url of the ccdb repository" };

  // Define output
  HistogramRegistry registry{ "registry" };

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
  std::vector<float> veccn = vecCnMult;
  std::vector<float> vecsn = vecSnMult;
  std::vector<float> veclp = vecL2nspMult;
  std::vector<float> veclm = vecL2nsmMult;
  std::vector<float> vecap = vecA2npMult;
  std::vector<float> vecam = vecA2nmMult;
  std::vector<float> vecwa = vecwacc;

  using TracksPID = soa::Join<aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr>;
  using aodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, TracksPID>>; // tracks filter
  using aodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::Mults>>;  // collisions filter
  using DaughterTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCPr, aod::pidTPCKa>;

  // Set the pt, mult and phi Axis;
  o2::framework::AxisSpec axis = axisPt;
  int nPtBins = axis.binEdges.size() - 1;
  double* PtBins = &(axis.binEdges)[0];
  TAxis* fPtAxis = new TAxis(nPtBins, PtBins);

  o2::framework::AxisSpec axisMult = axisMultiplicity;
  int nMultBins = axisMult.binEdges.size() - 1;
  double* MultBins = &(axisMult.binEdges)[0];
  TAxis* fMultAxis = new TAxis(nMultBins, MultBins);

  o2::framework::AxisSpec axisphi = axisPhi;
  int nPhiBins = axisphi.binEdges.size() - 1;
  double* PhiBins = &(axisphi.binEdges)[0];
  TAxis* fPhiAxis = new TAxis(nPhiBins, PhiBins);

  o2::framework::AxisSpec axisOmegamass = axisOmegaminusMass;
  int nOmegaMassBins = axisOmegamass.binEdges.size() - 1;
  double* OmegaMassBins = &(axisOmegamass.binEdges)[0];
  TAxis* fOmegaMass = new TAxis(nOmegaMassBins, OmegaMassBins);

  o2::framework::AxisSpec axisXimass = axisXiminusMass;
  int nXiMassBins = axisXimass.binEdges.size() - 1;
  double* XiMassBins = &(axisXimass.binEdges)[0];
  TAxis* fXiMass = new TAxis(nXiMassBins, XiMassBins);

  void init(InitContext const&) // Initialization
  {
    ccdb->setURL(url.value);
    ccdb->setCaching(true);
    ccdb->setCreatedNotAfter(nolaterthan.value);

    // Add some output objects to the histogram registry
    registry.add("hPhi", "", { HistType::kTH1D, { axisPhi } });
    registry.add("hEta", "", { HistType::kTH1D, { axisEta } });
    registry.add("hEtaPhi", "", { HistType::kTH2D, { axisEta, axisPhi } });
    registry.add("hVtxZ", "", { HistType::kTH1D, { axisVertex } });
    registry.add("hMult", "", { HistType::kTH1D, { { 3000, 0.5, 3000.5 } } });
    registry.add("hCent", "", { HistType::kTH1D, { { 90, 0, 90 } } });
    registry.add("hPt", "", { HistType::kTH1D, { axisPt } });
    // cumulant of flow
    registry.add("c22", ";Centrality  (%) ; C_{2}{2}", { HistType::kTProfile, { axisMultiplicity } });
    registry.add("c24", ";Centrality  (%) ; C_{2}{4}", { HistType::kTProfile, { axisMultiplicity } });
    // pt-diff cumulant of flow
    registry.add("Xic22dpt", ";pt ; C_{2}{2} ", { HistType::kTProfile2D, { axisPt, axisXiminusMass } });
    registry.add("Omegac22dpt", ";pt ; C_{2}{2} ", { HistType::kTProfile2D, { axisPt, axisOmegaminusMass } });
    // InvMass(GeV) of casc
    registry.add("InvMassOmegaMinus", "", { HistType::kTH3D, { axisPt, axisOmegaminusMass, axisEta } });
    registry.add("InvMassXiMinus", "", { HistType::kTH3D, { axisPt, axisXiminusMass, axisEta } });

    fGFW->AddRegion("full", -0.8, 0.8, 1, 1); // ("name", etamin, etamax, ptbinnum, bitmask)eta region -0.8 to 0.8
    // with (-0.5, 0.5) eta gap
    fGFW->AddRegion("refN10", -0.8, -0.5, 1, 1);
    fGFW->AddRegion("refP10", 0.5, 0.8, 1, 1);
    int nXiptMassBins = nPtBins * nXiMassBins;
    fGFW->AddRegion("poiXiP", 0.5, 0.8, nXiptMassBins, 2);
    int nOmegaptMassBins = nPtBins * nOmegaMassBins;
    fGFW->AddRegion("poiOmegaP", 0.5, 0.8, nOmegaptMassBins, 4);
    // pushback
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refP10 {2} refN10 {-2}", "ChFull220", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refP10 {2 2} refN10 {-2 -2}", "ChFull240", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiXiP {2} refN10 {-2}", "Ch10Gap22", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiOmegaP {2} refN10 {-2}", "Ch10Gap22", kTRUE));
    //
    fGFW->CreateRegions(); // finalize the initialization

    // used for event selection
    fMultPVCutLow = new TF1("fMultPVCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x - 3.5*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)", 0, 100);
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
  void FillProfilepTOmega(const GFW::CorrConfig& corrconf, const ConstStr<chars...>& tarName, const double& cent, const int& ptbin)
  {
    float dnx = 0;
    float val = 0;
    for (int massbin = 1; massbin <= nOmegaMassBins; massbin++) {
      dnx = fGFW->Calculate(corrconf, (ptbin - 1) + ((massbin - 1) * nPtBins), kTRUE).real();
      if (dnx == 0)
        continue;
      val = fGFW->Calculate(corrconf, (ptbin - 1) + ((massbin - 1) * nPtBins), kFALSE).real() / dnx;
      if (TMath::Abs(val) < 1) {
        registry.fill(tarName, fPtAxis->GetBinCenter(ptbin), fOmegaMass->GetBinCenter(massbin), val, dnx);
      }
    }
    return;
  }
  template <char... chars>
  void FillProfilepTXi(const GFW::CorrConfig& corrconf, const ConstStr<chars...>& tarName, const double& cent, const int& ptbin)
  {
    float dnx = 0;
    float val = 0;
    for (int massbin = 1; massbin <= nXiMassBins; massbin++) {
      dnx = fGFW->Calculate(corrconf, (ptbin - 1) + ((massbin - 1) * nPtBins), kTRUE).real();
      if (dnx == 0)
        continue;
      val = fGFW->Calculate(corrconf, (ptbin - 1) + ((massbin - 1) * nPtBins), kFALSE).real() / dnx;
      if (TMath::Abs(val) < 1) {
        registry.fill(tarName, fPtAxis->GetBinCenter(ptbin), fXiMass->GetBinCenter(massbin), val, dnx);
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

  bool setCurrentParticleWeights(float& weight_nue, float& weight_nua, float phi, float eta, float pt, float vtxz)
  {
    float eff = 1.;
    if (mEfficiency)
      eff = mEfficiency->GetBinContent(mEfficiency->FindBin(pt));
    else
      eff = 1.0;
    if (eff == 0)
      return false;
    weight_nue = 1. / eff;
    if (mAcceptance)
      weight_nua = mAcceptance->GetNUA(phi, eta, vtxz);
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
      return 0;
    }
    if (!collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      // reject collisions close to Time Frame borders
      // https://its.cern.ch/jira/browse/O2-4623
      return 0;
    }
    if (!collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      // reject events affected by the ITS ROF border
      // https://its.cern.ch/jira/browse/O2-4309
      return 0;
    }
    if (!collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      // rejects collisions which are associated with the same "found-by-T0" bunch crossing
      // https://indico.cern.ch/event/1396220/#1-event-selection-with-its-rof
      return 0;
    }
    if (!collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      // removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference
      // use this cut at low multiplicities with caution
      return 0;
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
      return 0;
    if (multNTracksPV < fMultPVCutLow->Eval(centrality))
      return 0;
    if (multNTracksPV > fMultPVCutHigh->Eval(centrality))
      return 0;
    if (multTrk < fMultCutLow->Eval(centrality))
      return 0;
    if (multTrk > fMultCutHigh->Eval(centrality))
      return 0;

    // V0A T0A 5 sigma cut
    if (abs(collision.multFV0A() - fT0AV0AMean->Eval(collision.multFT0A())) > 5 * fT0AV0ASigma->Eval(collision.multFT0A()))
      return 0;

    return 1;
  }

  void process(aodCollisions::iterator const& collision, aod::BCsWithTimestamps const&, aodTracks const& tracks, aod::CascDataExt const& Cascades, DaughterTracks& dtracks)
  {
    int Ntot = tracks.size();
    if (Ntot < 1)
      return;
    fGFW->Clear();
    const auto cent = collision.centFT0C();

    // auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    if (eventSelected(collision, tracks.size(), cent))
      return;
    // loadCorrections(bc.timestamp());
    float vtxz = collision.posZ();
    registry.fill(HIST("hVtxZ"), vtxz);
    registry.fill(HIST("hMult"), Ntot);
    registry.fill(HIST("hCent"), collision.centFT0C());

    float weff = 1;
    float wacc = 1;
    // fill GFW ref flow
    for (auto& track : tracks) {
      //添加内容
      if (!setCurrentParticleWeights(weff, wacc, track.phi(), track.eta(), track.pt(), vtxz))
        continue;
      //
      int phibin = fPhiAxis->FindBin(track.phi()) - 1;
      wacc = 1 / vecwa[phibin];
      registry.fill(HIST("hPhi"), track.phi());
      registry.fill(HIST("hEta"), track.eta());
      registry.fill(HIST("hEtaPhi"), track.eta(), track.phi());
      registry.fill(HIST("hPt"), track.pt());
      int ptbin = fPtAxis->FindBin(track.pt()) - 1;
      if ((track.pt() > cfgCutPtMin) && (track.pt() < cfgCutPtMax)) {
        fGFW->Fill(track.eta(), ptbin, track.phi(), wacc * weff, 1); //(eta, ptbin, phi, wacc*weff, bitmask)
      }
    }
    // fill GFW of casc flow
    for (auto& casc : Cascades) {
      auto bachelor = casc.bachelor_as<DaughterTracks>();
      auto posdau = casc.posTrack_as<DaughterTracks>();
      auto negdau = casc.negTrack_as<DaughterTracks>();
      // check TPC
      if (checkDauTPC && (!posdau.hasTPC() || !negdau.hasTPC() || !bachelor.hasTPC())) {
        continue;
      }
      if (TMath::Abs(casc.yOmega()) < Casc_rapidity && TMath::Abs(bachelor.tpcNSigmaKa()) < NSigmaCascKaon && TMath::Abs(posdau.tpcNSigmaPr()) < NSigmaCascProton && TMath::Abs(negdau.tpcNSigmaPi()) < NSigmaCascPion) {
        registry.fill(HIST("InvMassOmegaMinus"), casc.pt(), casc.mOmega(), casc.eta());
        if ((casc.pt() < cfgCutPtPOIMax) && (casc.pt() > cfgCutPtPOIMin) && (casc.mOmega() > 1.63) && (casc.mOmega() < 1.71)) {
          fGFW->Fill(casc.eta(), fPtAxis->FindBin(casc.pt()) - 1 + ((fOmegaMass->FindBin(casc.mOmega()) - 1) * nPtBins), casc.phi(), wacc * weff, 4);
        }
      }
      if (TMath::Abs(casc.yXi()) < Casc_rapidity && TMath::Abs(bachelor.tpcNSigmaKa()) < NSigmaCascKaon && TMath::Abs(posdau.tpcNSigmaPr()) < NSigmaCascProton && TMath::Abs(negdau.tpcNSigmaPi()) < NSigmaCascPion) {
        registry.fill(HIST("InvMassXiMinus"), casc.pt(), casc.mXi(), casc.eta());
        if ((casc.pt() < cfgCutPtPOIMax) && (casc.pt() > cfgCutPtPOIMin) && (casc.mXi() > 1.30) && (casc.mXi() < 1.37)) {
          fGFW->Fill(casc.eta(), fPtAxis->FindBin(casc.pt()) - 1 + ((fXiMass->FindBin(casc.mXi()) - 1) * nPtBins), casc.phi(), wacc * weff, 2);
        }
      }
    }
    // Filling cumulant with ROOT TProfile
    FillProfile(corrconfigs.at(0), HIST("c22"), cent);
    FillProfile(corrconfigs.at(1), HIST("c24"), cent);
    if (cent > 40)
      return;
    if (cent < 30)
      return;
    for (int i = 1; i <= nPtBins; i++) // loop for all ptBins
    {
      FillProfilepTXi(corrconfigs.at(2), HIST("Xic22dpt"), cent, i);
      FillProfilepTOmega(corrconfigs.at(3), HIST("Omegac22dpt"), cent, i);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowGFWPbPb>(cfgc)
  };
}