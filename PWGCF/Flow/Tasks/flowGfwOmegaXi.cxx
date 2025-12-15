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

/// \file   flowGfwOmegaXi.cxx
/// \author Fuchun Cui(fcui@cern.ch)
/// \since  Sep/13/2024
/// \brief  This task is to caculate V0s and cascades flow by GenericFramework

#include "GFW.h"
#include "GFWCumulant.h"
#include "GFWPowerArray.h"
#include "GFWWeights.h"

#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGMM/Mult/DataModel/Index.h"

#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/EventPlaneHelper.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"
#include <CCDB/BasicCCDBManager.h>

#include "TList.h"
#include <TF1.h>
#include <TF2.h>
#include <TPDGCode.h>
#include <TProfile.h>
#include <TRandom3.h>

#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
namespace
{
std::shared_ptr<TProfile> refc22[10];
std::shared_ptr<TProfile> refc24[10];
std::shared_ptr<TProfile> refc22Full[10];
std::shared_ptr<TProfile> refc32[10];
std::shared_ptr<TProfile3D> k0sc22[10];
std::shared_ptr<TProfile3D> k0sc24[10];
std::shared_ptr<TProfile3D> k0sc22Full[10];
std::shared_ptr<TProfile3D> k0sc32[10];
std::shared_ptr<TProfile3D> lambdac22[10];
std::shared_ptr<TProfile3D> lambdac24[10];
std::shared_ptr<TProfile3D> lambdac22Full[10];
std::shared_ptr<TProfile3D> lambdac32[10];
std::shared_ptr<TProfile3D> xic22[10];
std::shared_ptr<TProfile3D> xic24[10];
std::shared_ptr<TProfile3D> xic22Full[10];
std::shared_ptr<TProfile3D> xic32[10];
std::shared_ptr<TProfile3D> omegac22[10];
std::shared_ptr<TProfile3D> omegac24[10];
std::shared_ptr<TProfile3D> omegac22Full[10];
std::shared_ptr<TProfile3D> omegac32[10];

} // namespace

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct FlowGfwOmegaXi {

  O2_DEFINE_CONFIGURABLE(cfgCutVertex, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutChi2prTPCcls, float, 2.5, "Chi2 per TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgMassBins, std::vector<int>, (std::vector<int>{80, 32, 14, 16}), "Number of K0s, Lambda, Xi, Omega mass axis bins for c22")
  O2_DEFINE_CONFIGURABLE(cfgDeltaPhiLocDen, int, 3, "Number of delta phi for local density, 200 bins in 2 pi")

  struct : ConfigurableGroup {
    std::string prefix = "v0BuilderOpts";
    // topological cut for V0
    O2_DEFINE_CONFIGURABLE(cfgv0_radius, float, 5.0f, "minimum decay radius")
    O2_DEFINE_CONFIGURABLE(cfgv0_v0cospa, float, 0.995f, "minimum cosine of pointing angle")
    O2_DEFINE_CONFIGURABLE(cfgv0_dcadautopv, float, 0.1f, "minimum daughter DCA to PV")
    O2_DEFINE_CONFIGURABLE(cfgv0_dcav0dau, float, 0.5f, "maximum DCA among V0 daughters")
    O2_DEFINE_CONFIGURABLE(cfgv0_mk0swindow, float, 0.1f, "Invariant mass window of K0s")
    O2_DEFINE_CONFIGURABLE(cfgv0_mlambdawindow, float, 0.04f, "Invariant mass window of lambda")
    O2_DEFINE_CONFIGURABLE(cfgv0_ArmPodocut, float, 0.2f, "Armenteros Podolski cut for K0")
    O2_DEFINE_CONFIGURABLE(cfgv0_compmassrejLambda, float, 0.01f, "competing mass rejection of lambda")
    O2_DEFINE_CONFIGURABLE(cfgv0_compmassrejK0s, float, 0.005f, "competing mass rejection of K0s")
  } v0BuilderOpts;

  struct : ConfigurableGroup {
    std::string prefix = "cascBuilderOpts";
    // topological cut for cascade
    O2_DEFINE_CONFIGURABLE(cfgcasc_radius, float, 0.5f, "minimum decay radius")
    O2_DEFINE_CONFIGURABLE(cfgcasc_casccospa, float, 0.999f, "minimum cosine of pointing angle")
    O2_DEFINE_CONFIGURABLE(cfgcasc_v0cospa, float, 0.998f, "minimum cosine of pointing angle")
    O2_DEFINE_CONFIGURABLE(cfgcasc_dcav0topv, float, 0.01f, "minimum daughter DCA to PV")
    O2_DEFINE_CONFIGURABLE(cfgcasc_dcabachtopv, float, 0.01f, "minimum bachelor DCA to PV")
    O2_DEFINE_CONFIGURABLE(cfgcasc_dcacascdau, float, 0.3f, "maximum DCA among cascade daughters")
    O2_DEFINE_CONFIGURABLE(cfgcasc_dcav0dau, float, 1.0f, "maximum DCA among V0 daughters")
    O2_DEFINE_CONFIGURABLE(cfgcasc_mlambdawindow, float, 0.04f, "Invariant mass window of lambda")
    O2_DEFINE_CONFIGURABLE(cfgcasc_compmassrej, float, 0.008f, "competing mass rejection of cascades")
  } cascBuilderOpts;

  struct : ConfigurableGroup {
    std::string prefix = "trkQualityOpts";
    // track selections
    O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta range for tracks")
    O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMin, float, 0.2f, "Minimal pT for poi tracks")
    O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMax, float, 10.0f, "Maximal pT for poi tracks")
    O2_DEFINE_CONFIGURABLE(cfgCutPtMin, float, 0.2f, "Minimal pT for ref tracks")
    O2_DEFINE_CONFIGURABLE(cfgCutPtMax, float, 10.0f, "Maximal pT for ref tracks")
    O2_DEFINE_CONFIGURABLE(cfgCutPtDauMin, float, 0.2f, "Minimal pT for daughter tracks")
    O2_DEFINE_CONFIGURABLE(cfgCutPtDauMax, float, 10.0f, "Maximal pT for daughter tracks")
    O2_DEFINE_CONFIGURABLE(cfgCutPtK0sMin, float, 0.2f, "Minimal pT for K0s")
    O2_DEFINE_CONFIGURABLE(cfgCutPtK0sMax, float, 10.0f, "Maximal pT for K0s")
    O2_DEFINE_CONFIGURABLE(cfgCutPtLambdaMin, float, 0.2f, "Minimal pT for Lambda")
    O2_DEFINE_CONFIGURABLE(cfgCutPtLambdaMax, float, 10.0f, "Maximal pT for Lambda")
    O2_DEFINE_CONFIGURABLE(cfgCutPtXiMin, float, 0.2f, "Minimal pT for Xi")
    O2_DEFINE_CONFIGURABLE(cfgCutPtXiMax, float, 10.0f, "Maximal pT for Xi")
    O2_DEFINE_CONFIGURABLE(cfgCutPtOmegaMin, float, 0.2f, "Minimal pT for Omega")
    O2_DEFINE_CONFIGURABLE(cfgCutPtOmegaMax, float, 10.0f, "Maximal pT for Omega")
    O2_DEFINE_CONFIGURABLE(cfgCutPtPIDDauMin, float, 0.15f, "Minimal pT for daughter PID")
    O2_DEFINE_CONFIGURABLE(cfgCutPtPIDbachMin, float, 0.15f, "Minimal pT for daughter PID")
    O2_DEFINE_CONFIGURABLE(cfgCutPtPIDdauLaPrMin, float, 0.15f, "Minimal pT for daughter PID")
    O2_DEFINE_CONFIGURABLE(cfgCutPtPIDdauLaPiMin, float, 0.15f, "Minimal pT for daughter PID")
    // track quality selections for daughter track
    O2_DEFINE_CONFIGURABLE(cfgITSNCls, int, 5, "check minimum number of ITS clusters")
    O2_DEFINE_CONFIGURABLE(cfgTPCNCls, int, 50, "check minimum number of TPC hits")
    O2_DEFINE_CONFIGURABLE(cfgTPCCrossedRows, int, 70, "check minimum number of TPC crossed rows")
    O2_DEFINE_CONFIGURABLE(cfgCheckGlobalTrack, bool, false, "check global track")
  } trkQualityOpts;

  struct : ConfigurableGroup {
    std::string prefix = "evtSelOpts";
    O2_DEFINE_CONFIGURABLE(cfgDoTVXinTRD, bool, true, "check kTVXinTRD")
    O2_DEFINE_CONFIGURABLE(cfgDoNoTimeFrameBorder, bool, true, "check kNoTimeFrameBorder")
    O2_DEFINE_CONFIGURABLE(cfgDoNoITSROFrameBorder, bool, true, "check kNoITSROFrameBorder")
    O2_DEFINE_CONFIGURABLE(cfgDoNoSameBunchPileup, bool, true, "check kNoITSROFrameBorder")
    O2_DEFINE_CONFIGURABLE(cfgDoIsGoodZvtxFT0vsPV, bool, true, "check kIsGoodZvtxFT0vsPV")
    O2_DEFINE_CONFIGURABLE(cfgDoNoCollInTimeRangeStandard, bool, true, "check kNoCollInTimeRangeStandard")
    O2_DEFINE_CONFIGURABLE(cfgDoIsGoodITSLayersAll, bool, true, "check kIsGoodITSLayersAll")
    O2_DEFINE_CONFIGURABLE(cfgCutOccupancyHigh, int, 500, "High cut on TPC occupancy")
    O2_DEFINE_CONFIGURABLE(cfgDoMultPVCut, bool, true, "do multNTracksPV vs cent cut")
    O2_DEFINE_CONFIGURABLE(cfgMultPVCut, std::vector<float>, (std::vector<float>{3074.43, -106.192, 1.46176, -0.00968364, 2.61923e-05, 182.128, -7.43492, 0.193901, -0.00256715, 1.22594e-05}), "Used MultPVCut function parameter")
    O2_DEFINE_CONFIGURABLE(cfgDoV0AT0Acut, bool, true, "do V0A-T0A cut")
    O2_DEFINE_CONFIGURABLE(cfgCutminIR, float, -1, "cut min IR")
    O2_DEFINE_CONFIGURABLE(cfgCutmaxIR, float, -1, "cut max IR")
  } evtSeleOpts;

  O2_DEFINE_CONFIGURABLE(cfgCasc_rapidity, float, 0.5, "rapidity")
  O2_DEFINE_CONFIGURABLE(cfgNSigmapid, std::vector<float>, (std::vector<float>{9, 9, 9, 3, 3, 3, 3, 3, 3}), "tpc, tof and its NSigma for Pion Proton Kaon")
  O2_DEFINE_CONFIGURABLE(cfgAcceptancePath, std::vector<std::string>, (std::vector<std::string>{"Users/f/fcui/NUA/NUAREFPartical", "Users/f/fcui/NUA/NUAK0s", "Users/f/fcui/NUA/NUALambda", "Users/f/fcui/NUA/NUAXi", "Users/f/fcui/NUA/NUAOmega"}), "CCDB path to acceptance object")
  O2_DEFINE_CONFIGURABLE(cfgEfficiencyPath, std::vector<std::string>, (std::vector<std::string>{"PathtoRef"}), "CCDB path to efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgLocDenParaXi, std::vector<double>, (std::vector<double>{-0.000986187, -3.86861, -0.000912481, -3.29206, -0.000859271, -2.89389, -0.000817039, -2.61201, -0.000788792, -2.39079, -0.000780182, -2.19276, -0.000750457, -2.07205, -0.000720279, -1.96865, -0.00073247, -1.85642, -0.000695091, -1.82625, -0.000693332, -1.72679, -0.000681225, -1.74305, -0.000652818, -1.92608, -0.000618892, -2.31985}), "Local density efficiency function parameter for Xi, exp(Ax + B)")
  O2_DEFINE_CONFIGURABLE(cfgLocDenParaOmega, std::vector<double>, (std::vector<double>{-0.000444324, -6.0424, -0.000566208, -5.42168, -0.000580338, -4.96967, -0.000721054, -4.41994, -0.000626394, -4.27934, -0.000652167, -3.9543, -0.000592327, -3.79053, -0.000544721, -3.73292, -0.000613419, -3.43849, -0.000402506, -3.47687, -0.000602687, -3.24491, -0.000460848, -3.056, -0.00039428, -2.35188, -0.00041908, -2.03642}), "Local density efficiency function parameter for Omega, exp(Ax + B)")
  O2_DEFINE_CONFIGURABLE(cfgLocDenParaK0s, std::vector<double>, (std::vector<double>{-0.00043057, -3.2435, -0.000385085, -2.97687, -0.000350298, -2.81502, -0.000326159, -2.71091, -0.000299563, -2.65448, -0.000294284, -2.60865, -0.000277938, -2.589, -0.000277091, -2.56983, -0.000272783, -2.56825, -0.000252706, -2.58996, -0.000247834, -2.63158, -0.00024379, -2.76976, -0.000286468, -2.92484, -0.000310149, -3.27746}), "Local density efficiency function parameter for K0s, exp(Ax + B)")
  O2_DEFINE_CONFIGURABLE(cfgLocDenParaLambda, std::vector<double>, (std::vector<double>{-0.000510948, -4.4846, -0.000460629, -4.14465, -0.000433729, -3.94173, -0.000412751, -3.81839, -0.000411211, -3.72502, -0.000401511, -3.68426, -0.000407461, -3.67005, -0.000379371, -3.71153, -0.000392828, -3.73214, -0.000403996, -3.80717, -0.000403376, -3.90917, -0.000354624, -4.34629, -0.000477606, -4.66307, -0.000541139, -4.61364}), "Local density efficiency function parameter for Lambda, exp(Ax + B)")
  O2_DEFINE_CONFIGURABLE(cfgRunNumbers, std::vector<int>, (std::vector<int>{544095, 544098, 544116, 544121, 544122, 544123, 544124}), "Preconfigured run numbers")
  // switch
  O2_DEFINE_CONFIGURABLE(cfgDoAccEffCorr, bool, false, "do acc and eff corr")
  O2_DEFINE_CONFIGURABLE(cfgDoLocDenCorr, bool, false, "do local density corr")
  O2_DEFINE_CONFIGURABLE(cfgDoJackknife, bool, false, "do jackknife")
  O2_DEFINE_CONFIGURABLE(cfgOutputV0, bool, true, "Fill and output V0s flow")
  O2_DEFINE_CONFIGURABLE(cfgOutputCasc, bool, true, "Fill and output cascades flow")
  O2_DEFINE_CONFIGURABLE(cfgOutputNUAWeights, bool, false, "Fill and output NUA weights")
  O2_DEFINE_CONFIGURABLE(cfgOutputrunbyrun, bool, false, "Fill and output NUA weights run by run")
  O2_DEFINE_CONFIGURABLE(cfgOutputLocDenWeights, bool, false, "Fill and output local density weights")
  O2_DEFINE_CONFIGURABLE(cfgOutputQA, bool, false, "do QA")

  ConfigurableAxis cfgaxisVertex{"cfgaxisVertex", {20, -10, 10}, "vertex axis for histograms"};
  ConfigurableAxis cfgaxisPhi{"cfgaxisPhi", {60, 0.0, constants::math::TwoPI}, "phi axis for histograms"};
  ConfigurableAxis cfgaxisEta{"cfgaxisEta", {40, -1., 1.}, "eta axis for histograms"};
  ConfigurableAxis cfgaxisPt{"cfgaxisPt", {VARIABLE_WIDTH, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.50, 4.00, 4.50, 5.00, 5.50, 6.00, 10.0}, "pt (GeV)"};
  ConfigurableAxis cfgaxisPtXi{"cfgaxisPtXi", {VARIABLE_WIDTH, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.9, 4.9, 5.9, 9.9}, "pt (GeV)"};
  ConfigurableAxis cfgaxisPtOmega{"cfgaxisPtOmega", {VARIABLE_WIDTH, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.9, 4.9, 5.9, 9.9}, "pt (GeV)"};
  ConfigurableAxis cfgaxisPtK0s{"cfgaxisPtK0s", {VARIABLE_WIDTH, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.9, 4.9, 5.9, 9.9}, "pt (GeV)"};
  ConfigurableAxis cfgaxisPtLambda{"cfgaxisPtLambda", {VARIABLE_WIDTH, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.9, 4.9, 5.9, 9.9}, "pt (GeV)"};
  ConfigurableAxis cfgaxisOmegaMassforflow{"cfgaxisOmegaMassforflow", {16, 1.63f, 1.71f}, "Inv. Mass (GeV)"};
  ConfigurableAxis cfgaxisXiMassforflow{"cfgaxisXiMassforflow", {14, 1.3f, 1.37f}, "Inv. Mass (GeV)"};
  ConfigurableAxis cfgaxisK0sMassforflow{"cfgaxisK0sMassforflow", {40, 0.4f, 0.6f}, "Inv. Mass (GeV)"};
  ConfigurableAxis cfgaxisLambdaMassforflow{"cfgaxisLambdaMassforflow", {32, 1.08f, 1.16f}, "Inv. Mass (GeV)"};
  ConfigurableAxis cfgaxisNch{"cfgaxisNch", {3000, 0.5, 3000.5}, "Nch"};
  ConfigurableAxis cfgaxisLocalDensity{"cfgaxisLocalDensity", {200, 0, 600}, "local density"};

  AxisSpec axisMultiplicity{{0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90}, "Centrality (%)"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = (nabs(aod::track::eta) < trkQualityOpts.cfgCutEta.value) && (aod::track::pt > trkQualityOpts.cfgCutPtPOIMin.value) && (aod::track::pt < trkQualityOpts.cfgCutPtPOIMax.value) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && (aod::track::tpcChi2NCl < cfgCutChi2prTPCcls);

  using TracksPID = soa::Join<aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr>;
  using AodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, o2::aod::TrackSelectionExtension, aod::TracksExtra, TracksPID, aod::TracksIU>>; // tracks filter
  using AodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::MultsRun3>>;                                               // collisions filter
  using DaughterTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, TracksPID, aod::TrackSelection, o2::aod::TrackSelectionExtension>;

  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;
  ctpRateFetcher rateFetcher;
  O2_DEFINE_CONFIGURABLE(cfgnolaterthan, int64_t, std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object")
  O2_DEFINE_CONFIGURABLE(cfgurl, std::string, "http://alice-ccdb.cern.ch", "url of the ccdb repository")

  // Define output
  HistogramRegistry registry{"registry"};
  OutputObj<GFWWeights> fWeightsREF{GFWWeights("weightsREF")};
  OutputObj<GFWWeights> fWeightsK0s{GFWWeights("weightsK0s")};
  OutputObj<GFWWeights> fWeightsLambda{GFWWeights("weightsLambda")};
  OutputObj<GFWWeights> fWeightsXi{GFWWeights("weightsXi")};
  OutputObj<GFWWeights> fWeightsOmega{GFWWeights("weightsOmega")};

  // define global variables
  GFW* fGFW = new GFW(); // GFW class used from main src
  std::vector<GFW::CorrConfig> corrconfigs;
  std::vector<std::string> cfgAcceptance;
  std::vector<std::string> cfgEfficiency;
  std::vector<float> cfgNSigma;
  std::vector<float> cfgMultPVCutPara;
  std::vector<int> cfgmassbins;
  std::vector<int> runNumbers;
  std::map<int, std::vector<std::shared_ptr<TH1>>> th1sList;
  std::map<int, std::vector<std::shared_ptr<TH3>>> th3sList;
  enum OutputTH1Names {
    // here are TProfiles for vn-pt correlations that are not implemented in GFW
    hPhi = 0,
    hPhicorr,
    hPhiK0s,
    hPhiLambda,
    hPhiXi,
    hPhiOmega,
    hPhiK0scorr,
    hPhiLambdacorr,
    hPhiXicorr,
    hPhiOmegacorr,
    kCount_TH1Names
  };

  enum OutputTH3Names {
    hPhiEtaVtxz = 0,
    hPhiEtaVtxzK0s,
    hPhiEtaVtxzLambda,
    hPhiEtaVtxzXi,
    hPhiEtaVtxzOmega,
    kCount_TH3Names
  };

  std::vector<GFWWeights*> mAcceptance;
  std::vector<TH1D*> mEfficiency;
  bool correctionsLoaded = false;

  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;
  TF1* fT0AV0AMean = nullptr;
  TF1* fT0AV0ASigma = nullptr;

  // Declare the pt, mult and phi Axis;
  int nPtBins = 0;
  TAxis* fPtAxis = nullptr;

  int nXiPtBins = 0;
  TAxis* fXiPtAxis = nullptr;

  int nOmegaPtBins = 0;
  TAxis* fOmegaPtAxis = nullptr;

  int nK0sPtBins = 0;
  TAxis* fK0sPtAxis = nullptr;

  int nLambdaPtBins = 0;
  TAxis* fLambdaPtAxis = nullptr;

  TAxis* fMultAxis = nullptr;

  TAxis* fOmegaMass = nullptr;

  TAxis* fXiMass = nullptr;

  TAxis* fK0sMass = nullptr;

  TAxis* fLambdaMass = nullptr;

  void init(InitContext const&) // Initialization
  {
    ccdb->setURL(cfgurl.value);
    ccdb->setCaching(true);
    ccdb->setCreatedNotAfter(cfgnolaterthan.value);

    cfgAcceptance = cfgAcceptancePath;
    cfgEfficiency = cfgEfficiencyPath;
    cfgNSigma = cfgNSigmapid;
    cfgmassbins = cfgMassBins;
    cfgMultPVCutPara = evtSeleOpts.cfgMultPVCut;

    // Set the pt, mult and phi Axis;
    o2::framework::AxisSpec axisPt = cfgaxisPt;
    nPtBins = axisPt.binEdges.size() - 1;
    fPtAxis = new TAxis(nPtBins, &(axisPt.binEdges)[0]);

    o2::framework::AxisSpec axisXiPt = cfgaxisPtXi;
    nXiPtBins = axisXiPt.binEdges.size() - 1;
    fXiPtAxis = new TAxis(nXiPtBins, &(axisXiPt.binEdges)[0]);

    o2::framework::AxisSpec axisOmegaPt = cfgaxisPtOmega;
    nOmegaPtBins = axisOmegaPt.binEdges.size() - 1;
    fOmegaPtAxis = new TAxis(nOmegaPtBins, &(axisOmegaPt.binEdges)[0]);

    o2::framework::AxisSpec axisK0sPt = cfgaxisPtK0s;
    nK0sPtBins = axisK0sPt.binEdges.size() - 1;
    fK0sPtAxis = new TAxis(nK0sPtBins, &(axisK0sPt.binEdges)[0]);

    o2::framework::AxisSpec axisLambdaPt = cfgaxisPtLambda;
    nLambdaPtBins = axisLambdaPt.binEdges.size() - 1;
    fLambdaPtAxis = new TAxis(nLambdaPtBins, &(axisLambdaPt.binEdges)[0]);

    o2::framework::AxisSpec axisMult = axisMultiplicity;
    int nMultBins = axisMult.binEdges.size() - 1;
    fMultAxis = new TAxis(nMultBins, &(axisMult.binEdges)[0]);

    fOmegaMass = new TAxis(cfgmassbins[3], 1.63, 1.71);

    fXiMass = new TAxis(cfgmassbins[2], 1.29, 1.36);

    fK0sMass = new TAxis(cfgmassbins[0], 0.4, 0.6);

    fLambdaMass = new TAxis(cfgmassbins[1], 1.08, 1.16);

    // Add some output objects to the histogram registry
    registry.add("hPhi", "", {HistType::kTH1D, {cfgaxisPhi}});
    registry.add("hPhicorr", "", {HistType::kTH1D, {cfgaxisPhi}});
    registry.add("hPhiK0s", "", {HistType::kTH1D, {cfgaxisPhi}});
    registry.add("hPhiLambda", "", {HistType::kTH1D, {cfgaxisPhi}});
    registry.add("hPhiXi", "", {HistType::kTH1D, {cfgaxisPhi}});
    registry.add("hPhiOmega", "", {HistType::kTH1D, {cfgaxisPhi}});
    registry.add("hPhiK0scorr", "", {HistType::kTH1D, {cfgaxisPhi}});
    registry.add("hPhiLambdacorr", "", {HistType::kTH1D, {cfgaxisPhi}});
    registry.add("hPhiXicorr", "", {HistType::kTH1D, {cfgaxisPhi}});
    registry.add("hPhiOmegacorr", "", {HistType::kTH1D, {cfgaxisPhi}});
    registry.add("hEta", "", {HistType::kTH1D, {cfgaxisEta}});
    registry.add("hVtxZ", "", {HistType::kTH1D, {cfgaxisVertex}});
    registry.add("hMult", "", {HistType::kTH1D, {cfgaxisNch}});
    registry.add("hMultTPC", "", {HistType::kTH1D, {cfgaxisNch}});
    registry.add("hCent", "", {HistType::kTH1D, {{90, 0, 90}}});
    registry.add("hCentvsNch", "", {HistType::kTH2D, {{18, 0, 90}, cfgaxisNch}});
    registry.add("MC/hCentvsNchMC", "", {HistType::kTH2D, {{18, 0, 90}, cfgaxisNch}});
    registry.add("hCentvsMultTPC", "", {HistType::kTH2D, {{18, 0, 90}, cfgaxisNch}});
    registry.add("MC/hCentvsMultTPCMC", "", {HistType::kTH2D, {{18, 0, 90}, cfgaxisNch}});
    registry.add("hNTracksPVvsCentrality", "", {HistType::kTH2D, {{5000, 0, 5000}, axisMultiplicity}});
    registry.add("hPt", "", {HistType::kTH1D, {cfgaxisPt}});

    if (cfgOutputrunbyrun) {
      runNumbers = cfgRunNumbers;
      for (const auto& runNumber : runNumbers) {
        if (cfgOutputQA) {
          std::vector<std::shared_ptr<TH1>> histosPhi(kCount_TH1Names);
          histosPhi[hPhi] = registry.add<TH1>(Form("%d/hPhi", runNumber), "", {HistType::kTH1D, {cfgaxisPhi}});
          histosPhi[hPhicorr] = registry.add<TH1>(Form("%d/hPhicorr", runNumber), "", {HistType::kTH1D, {cfgaxisPhi}});
          histosPhi[hPhiK0s] = registry.add<TH1>(Form("%d/hPhiK0s", runNumber), "", {HistType::kTH1D, {cfgaxisPhi}});
          histosPhi[hPhiLambda] = registry.add<TH1>(Form("%d/hPhiLambda", runNumber), "", {HistType::kTH1D, {cfgaxisPhi}});
          histosPhi[hPhiXi] = registry.add<TH1>(Form("%d/hPhiXi", runNumber), "", {HistType::kTH1D, {cfgaxisPhi}});
          histosPhi[hPhiOmega] = registry.add<TH1>(Form("%d/hPhiOmega", runNumber), "", {HistType::kTH1D, {cfgaxisPhi}});
          histosPhi[hPhiK0scorr] = registry.add<TH1>(Form("%d/hPhiK0scorr", runNumber), "", {HistType::kTH1D, {cfgaxisPhi}});
          histosPhi[hPhiLambdacorr] = registry.add<TH1>(Form("%d/hPhiLambdacorr", runNumber), "", {HistType::kTH1D, {cfgaxisPhi}});
          histosPhi[hPhiXicorr] = registry.add<TH1>(Form("%d/hPhiXicorr", runNumber), "", {HistType::kTH1D, {cfgaxisPhi}});
          histosPhi[hPhiOmegacorr] = registry.add<TH1>(Form("%d/hPhiOmegacorr", runNumber), "", {HistType::kTH1D, {cfgaxisPhi}});
          th1sList.insert(std::make_pair(runNumber, histosPhi));
        }

        std::vector<std::shared_ptr<TH3>> nuaTH3(kCount_TH3Names);
        nuaTH3[hPhiEtaVtxz] = registry.add<TH3>(Form("%d/hPhiEtaVtxz", runNumber), ";#varphi;#eta;v_{z}", {HistType::kTH3D, {cfgaxisPhi, {64, -1.6, 1.6}, cfgaxisVertex}});
        nuaTH3[hPhiEtaVtxzK0s] = registry.add<TH3>(Form("%d/hPhiEtaVtxzK0s", runNumber), ";#varphi;#eta;v_{z}", {HistType::kTH3D, {cfgaxisPhi, {64, -1.6, 1.6}, cfgaxisVertex}});
        nuaTH3[hPhiEtaVtxzLambda] = registry.add<TH3>(Form("%d/hPhiEtaVtxzLambda", runNumber), ";#varphi;#eta;v_{z}", {HistType::kTH3D, {cfgaxisPhi, {64, -1.6, 1.6}, cfgaxisVertex}});
        nuaTH3[hPhiEtaVtxzXi] = registry.add<TH3>(Form("%d/hPhiEtaVtxzXi", runNumber), ";#varphi;#eta;v_{z}", {HistType::kTH3D, {cfgaxisPhi, {64, -1.6, 1.6}, cfgaxisVertex}});
        nuaTH3[hPhiEtaVtxzOmega] = registry.add<TH3>(Form("%d/hPhiEtaVtxzOmega", runNumber), ";#varphi;#eta;v_{z}", {HistType::kTH3D, {cfgaxisPhi, {64, -1.6, 1.6}, cfgaxisVertex}});
        th3sList.insert(std::make_pair(runNumber, nuaTH3));
      }
    }

    registry.add("hEventCount", "", {HistType::kTH1D, {{14, 0, 14}}});
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(1, "Filtered event");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(2, "after sel8");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(3, "after kTVXinTRD");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(4, "after kNoTimeFrameBorder");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(5, "after kNoITSROFrameBorder");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(6, "after kDoNoSameBunchPileup");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(7, "after kIsGoodZvtxFT0vsPV");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(8, "after kNoCollInTimeRangeStandard");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(9, "after kIsGoodITSLayersAll");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(10, "after MultPVCut");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(11, "after TPC occupancy cut");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(12, "after V0AT0Acut");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(13, "after IRmincut");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(14, "after IRmaxcut");
    registry.add("hInteractionRate", "", {HistType::kTH1D, {{1000, 0, 1000}}});

    // QA
    if (cfgOutputQA) {
      // V0 QA
      registry.add("QAhisto/V0/hqaV0radiusbefore", "", {HistType::kTH1D, {{200, 0, 200}}});
      registry.add("QAhisto/V0/hqaV0radiusafter", "", {HistType::kTH1D, {{200, 0, 200}}});
      registry.add("QAhisto/V0/hqaV0cosPAbefore", "", {HistType::kTH1D, {{1000, 0.95, 1}}});
      registry.add("QAhisto/V0/hqaV0cosPAafter", "", {HistType::kTH1D, {{1000, 0.95, 1}}});
      registry.add("QAhisto/V0/hqadcaV0daubefore", "", {HistType::kTH1D, {{100, 0, 1}}});
      registry.add("QAhisto/V0/hqadcaV0dauafter", "", {HistType::kTH1D, {{100, 0, 1}}});
      registry.add("QAhisto/V0/hqaarm_podobefore", "", {HistType::kTH2D, {{100, -1, 1}, {50, 0, 0.3}}});
      registry.add("QAhisto/V0/hqaarm_podoafter", "", {HistType::kTH2D, {{100, -1, 1}, {50, 0, 0.3}}});
      registry.add("QAhisto/V0/hqadcapostoPVbefore", "", {HistType::kTH1D, {{1000, -10, 10}}});
      registry.add("QAhisto/V0/hqadcapostoPVafter", "", {HistType::kTH1D, {{1000, -10, 10}}});
      registry.add("QAhisto/V0/hqadcanegtoPVbefore", "", {HistType::kTH1D, {{1000, -10, 10}}});
      registry.add("QAhisto/V0/hqadcanegtoPVafter", "", {HistType::kTH1D, {{1000, -10, 10}}});
      // Cascade QA
      registry.add("QAhisto/Xi/hqaCascRadiusbefore", "", {HistType::kTH1D, {{200, -10, 10}}});
      registry.add("QAhisto/Xi/hqaCascRadiusafter", "", {HistType::kTH1D, {{200, -10, 10}}});
      registry.add("QAhisto/Xi/hqaCasccosPAbefore", "", {HistType::kTH1D, {{1000, 0.95, 1}}});
      registry.add("QAhisto/Xi/hqaCasccosPAafter", "", {HistType::kTH1D, {{1000, 0.95, 1}}});
      registry.add("QAhisto/Xi/hqaCascV0cosPAbefore", "", {HistType::kTH1D, {{1000, 0.95, 1}}});
      registry.add("QAhisto/Xi/hqaCascV0cosPAafter", "", {HistType::kTH1D, {{1000, 0.95, 1}}});
      registry.add("QAhisto/Xi/hqadcaCascV0toPVbefore", "", {HistType::kTH1D, {{1000, -10, 10}}});
      registry.add("QAhisto/Xi/hqadcaCascV0toPVafter", "", {HistType::kTH1D, {{1000, -10, 10}}});
      registry.add("QAhisto/Xi/hqadcaCascBachtoPVbefore", "", {HistType::kTH1D, {{1000, -10, 10}}});
      registry.add("QAhisto/Xi/hqadcaCascBachtoPVafter", "", {HistType::kTH1D, {{1000, -10, 10}}});
      registry.add("QAhisto/Xi/hqadcaCascdaubefore", "", {HistType::kTH1D, {{100, 0, 1}}});
      registry.add("QAhisto/Xi/hqadcaCascdauafter", "", {HistType::kTH1D, {{100, 0, 1}}});
      registry.add("QAhisto/Xi/hqadcaCascV0daubefore", "", {HistType::kTH1D, {{100, 0, 1}}});
      registry.add("QAhisto/Xi/hqadcaCascV0dauafter", "", {HistType::kTH1D, {{100, 0, 1}}});

      registry.add("QAhisto/Omega/hqaCascRadiusbefore", "", {HistType::kTH1D, {{200, -10, 10}}});
      registry.add("QAhisto/Omega/hqaCascRadiusafter", "", {HistType::kTH1D, {{200, -10, 10}}});
      registry.add("QAhisto/Omega/hqaCasccosPAbefore", "", {HistType::kTH1D, {{1000, 0.95, 1}}});
      registry.add("QAhisto/Omega/hqaCasccosPAafter", "", {HistType::kTH1D, {{1000, 0.95, 1}}});
      registry.add("QAhisto/Omega/hqaCascV0cosPAbefore", "", {HistType::kTH1D, {{1000, 0.95, 1}}});
      registry.add("QAhisto/Omega/hqaCascV0cosPAafter", "", {HistType::kTH1D, {{1000, 0.95, 1}}});
      registry.add("QAhisto/Omega/hqadcaCascV0toPVbefore", "", {HistType::kTH1D, {{1000, -10, 10}}});
      registry.add("QAhisto/Omega/hqadcaCascV0toPVafter", "", {HistType::kTH1D, {{1000, -10, 10}}});
      registry.add("QAhisto/Omega/hqadcaCascBachtoPVbefore", "", {HistType::kTH1D, {{1000, -10, 10}}});
      registry.add("QAhisto/Omega/hqadcaCascBachtoPVafter", "", {HistType::kTH1D, {{1000, -10, 10}}});
      registry.add("QAhisto/Omega/hqadcaCascdaubefore", "", {HistType::kTH1D, {{100, 0, 1}}});
      registry.add("QAhisto/Omega/hqadcaCascdauafter", "", {HistType::kTH1D, {{100, 0, 1}}});
      registry.add("QAhisto/Omega/hqadcaCascV0daubefore", "", {HistType::kTH1D, {{100, 0, 1}}});
      registry.add("QAhisto/Omega/hqadcaCascV0dauafter", "", {HistType::kTH1D, {{100, 0, 1}}});
    }

    // cumulant of flow
    registry.add("c22", ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("c32", ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("c24", ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("c22Full", ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("c22dpt", ";Centrality  (%) ; C_{2}{2}", {HistType::kTProfile2D, {cfgaxisPt, axisMultiplicity}});
    registry.add("c24dpt", ";Centrality  (%) ; C_{2}{4}", {HistType::kTProfile2D, {cfgaxisPt, axisMultiplicity}});
    registry.add("c22Fulldpt", ";Centrality  (%) ; C_{2}{2}", {HistType::kTProfile2D, {cfgaxisPt, axisMultiplicity}});
    // pt-diff cumulant of flow
    if (cfgOutputCasc) {
      // v2
      registry.add("Xic22dpt", ";pt ; C_{2}{2} ", {HistType::kTProfile3D, {cfgaxisPtXi, cfgaxisXiMassforflow, axisMultiplicity}});
      registry.add("Omegac22dpt", ";pt ; C_{2}{2} ", {HistType::kTProfile3D, {cfgaxisPtOmega, cfgaxisOmegaMassforflow, axisMultiplicity}});
      registry.add("Xic24dpt", ";pt ; C_{2}{4} ", {HistType::kTProfile3D, {cfgaxisPtXi, cfgaxisXiMassforflow, axisMultiplicity}});
      registry.add("Omegac24dpt", ";pt ; C_{2}{4} ", {HistType::kTProfile3D, {cfgaxisPtOmega, cfgaxisOmegaMassforflow, axisMultiplicity}});
      registry.add("Xic22Fulldpt", ";pt ; C_{2}{2} ", {HistType::kTProfile3D, {cfgaxisPtXi, cfgaxisXiMassforflow, axisMultiplicity}});
      registry.add("Omegac22Fulldpt", ";pt ; C_{2}{2} ", {HistType::kTProfile3D, {cfgaxisPtOmega, cfgaxisOmegaMassforflow, axisMultiplicity}});
      registry.add("Xic24_gapdpt", ";pt ; C_{2}{4} ", {HistType::kTProfile3D, {cfgaxisPtXi, cfgaxisXiMassforflow, axisMultiplicity}});
      registry.add("Omegac24_gapdpt", ";pt ; C_{2}{4} ", {HistType::kTProfile3D, {cfgaxisPtOmega, cfgaxisOmegaMassforflow, axisMultiplicity}});
      // v3
      registry.add("Xic32dpt", ";pt ; C_{2}{2} ", {HistType::kTProfile3D, {cfgaxisPtXi, cfgaxisXiMassforflow, axisMultiplicity}});
      registry.add("Omegac32dpt", ";pt ; C_{2}{2} ", {HistType::kTProfile3D, {cfgaxisPtOmega, cfgaxisOmegaMassforflow, axisMultiplicity}});
    }
    if (cfgOutputV0) {
      // v2
      registry.add("K0sc22dpt", ";pt ; C_{2}{2} ", {HistType::kTProfile3D, {cfgaxisPtK0s, cfgaxisK0sMassforflow, axisMultiplicity}});
      registry.add("Lambdac22dpt", ";pt ; C_{2}{2} ", {HistType::kTProfile3D, {cfgaxisPtLambda, cfgaxisLambdaMassforflow, axisMultiplicity}});
      registry.add("K0sc24dpt", ";pt ; C_{2}{4} ", {HistType::kTProfile3D, {cfgaxisPtK0s, cfgaxisK0sMassforflow, axisMultiplicity}});
      registry.add("Lambdac24dpt", ";pt ; C_{2}{4} ", {HistType::kTProfile3D, {cfgaxisPtK0s, cfgaxisLambdaMassforflow, axisMultiplicity}});
      registry.add("K0sc22Fulldpt", ";pt ; C_{2}{2} ", {HistType::kTProfile3D, {cfgaxisPtK0s, cfgaxisK0sMassforflow, axisMultiplicity}});
      registry.add("Lambdac22Fulldpt", ";pt ; C_{2}{2} ", {HistType::kTProfile3D, {cfgaxisPtLambda, cfgaxisLambdaMassforflow, axisMultiplicity}});
      // v3
      registry.add("K0sc32dpt", ";pt ; C_{2}{2} ", {HistType::kTProfile3D, {cfgaxisPtK0s, cfgaxisK0sMassforflow, axisMultiplicity}});
      registry.add("Lambdac32dpt", ";pt ; C_{2}{2} ", {HistType::kTProfile3D, {cfgaxisPtLambda, cfgaxisLambdaMassforflow, axisMultiplicity}});
    }
    // for Jackknife
    if (cfgDoJackknife) {
      int nsubevent = 10;
      for (int i = 1; i <= nsubevent; i++) {
        refc22[i - 1] = registry.add<TProfile>(Form("Jackknife/REF/c22_%d", i), ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile, {axisMultiplicity}});
        refc24[i - 1] = registry.add<TProfile>(Form("Jackknife/REF/c24_%d", i), ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile, {axisMultiplicity}});
        refc22Full[i - 1] = registry.add<TProfile>(Form("Jackknife/REF/c22Full_%d", i), ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile, {axisMultiplicity}});
        xic22[i - 1] = registry.add<TProfile3D>(Form("Jackknife/Xi/Xic22dpt_%d", i), ";pt ; C_{2}{2} ", {HistType::kTProfile3D, {cfgaxisPtXi, cfgaxisXiMassforflow, axisMultiplicity}});
        omegac22[i - 1] = registry.add<TProfile3D>(Form("Jackknife/Omega/Omegac22dpt_%d", i), ";pt ; C_{2}{2} ", {HistType::kTProfile3D, {cfgaxisPtOmega, cfgaxisOmegaMassforflow, axisMultiplicity}});
        k0sc22[i - 1] = registry.add<TProfile3D>(Form("Jackknife/K0s/K0sc22dpt_%d", i), ";pt ; C_{2}{2} ", {HistType::kTProfile3D, {cfgaxisPtK0s, cfgaxisK0sMassforflow, axisMultiplicity}});
        lambdac22[i - 1] = registry.add<TProfile3D>(Form("Jackknife/Lambda/Lambdac22dpt_%d", i), ";pt ; C_{2}{2} ", {HistType::kTProfile3D, {cfgaxisPtLambda, cfgaxisLambdaMassforflow, axisMultiplicity}});
        xic24[i - 1] = registry.add<TProfile3D>(Form("Jackknife/Xi/Xic24dpt_%d", i), ";pt ; C_{2}{4} ", {HistType::kTProfile3D, {cfgaxisPtXi, cfgaxisXiMassforflow, axisMultiplicity}});
        omegac24[i - 1] = registry.add<TProfile3D>(Form("Jackknife/Omega/Omegac24dpt_%d", i), ";pt ; C_{2}{4} ", {HistType::kTProfile3D, {cfgaxisPtOmega, cfgaxisOmegaMassforflow, axisMultiplicity}});
        k0sc24[i - 1] = registry.add<TProfile3D>(Form("Jackknife/K0s/K0sc24dpt_%d", i), ";pt ; C_{2}{4} ", {HistType::kTProfile3D, {cfgaxisPtK0s, cfgaxisK0sMassforflow, axisMultiplicity}});
        lambdac24[i - 1] = registry.add<TProfile3D>(Form("Jackknife/Lambda/Lambdac24dpt_%d", i), ";pt ; C_{2}{4} ", {HistType::kTProfile3D, {cfgaxisPtLambda, cfgaxisLambdaMassforflow, axisMultiplicity}});
        refc32[i - 1] = registry.add<TProfile>(Form("Jackknife/REF/c32_%d", i), ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile, {axisMultiplicity}});
        xic22Full[i - 1] = registry.add<TProfile3D>(Form("Jackknife/Xi/Xic22Fulldpt_%d", i), ";pt ; C_{2}{2} ", {HistType::kTProfile3D, {cfgaxisPtXi, cfgaxisXiMassforflow, axisMultiplicity}});
        omegac22Full[i - 1] = registry.add<TProfile3D>(Form("Jackknife/Omega/Omegac22Fulldpt_%d", i), ";pt ; C_{2}{2} ", {HistType::kTProfile3D, {cfgaxisPtOmega, cfgaxisOmegaMassforflow, axisMultiplicity}});
        k0sc22Full[i - 1] = registry.add<TProfile3D>(Form("Jackknife/K0s/K0sc22Fulldpt_%d", i), ";pt ; C_{2}{2} ", {HistType::kTProfile3D, {cfgaxisPtK0s, cfgaxisK0sMassforflow, axisMultiplicity}});
        lambdac22Full[i - 1] = registry.add<TProfile3D>(Form("Jackknife/Lambda/Lambdac22Fulldpt_%d", i), ";pt ; C_{2}{2} ", {HistType::kTProfile3D, {cfgaxisPtLambda, cfgaxisLambdaMassforflow, axisMultiplicity}});
        xic32[i - 1] = registry.add<TProfile3D>(Form("Jackknife/Xi/Xic32dpt_%d", i), ";pt ; C_{2}{2} ", {HistType::kTProfile3D, {cfgaxisPtXi, cfgaxisXiMassforflow, axisMultiplicity}});
        omegac32[i - 1] = registry.add<TProfile3D>(Form("Jackknife/Omega/Omegac32dpt_%d", i), ";pt ; C_{2}{2} ", {HistType::kTProfile3D, {cfgaxisPtOmega, cfgaxisOmegaMassforflow, axisMultiplicity}});
        k0sc32[i - 1] = registry.add<TProfile3D>(Form("Jackknife/K0s/K0sc32dpt_%d", i), ";pt ; C_{2}{2} ", {HistType::kTProfile3D, {cfgaxisPtK0s, cfgaxisK0sMassforflow, axisMultiplicity}});
        lambdac32[i - 1] = registry.add<TProfile3D>(Form("Jackknife/Lambda/Lambdac32dpt_%d", i), ";pt ; C_{2}{2} ", {HistType::kTProfile3D, {cfgaxisPtLambda, cfgaxisLambdaMassforflow, axisMultiplicity}});
      }
    }
    // MC True flow
    registry.add("MC/c22MC", ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile, {axisMultiplicity}});
    if (cfgOutputCasc) {
      registry.add("MC/Xic22dptMC", ";pt ; C_{2}{2} ", {HistType::kTProfile2D, {cfgaxisPtXi, axisMultiplicity}});
      registry.add("MC/Omegac22dptMC", ";pt ; C_{2}{2} ", {HistType::kTProfile2D, {cfgaxisPtOmega, axisMultiplicity}});
    }
    if (cfgOutputV0) {
      registry.add("MC/K0sc22dptMC", ";pt ; C_{2}{2} ", {HistType::kTProfile2D, {cfgaxisPtK0s, axisMultiplicity}});
      registry.add("MC/Lambdac22dptMC", ";pt ; C_{2}{2} ", {HistType::kTProfile2D, {cfgaxisPtLambda, axisMultiplicity}});
    }
    // InvMass(GeV) of casc and v0
    AxisSpec axisOmegaMass = {80, 1.63f, 1.71f, "Inv. Mass (GeV)"};
    AxisSpec axisXiMass = {80, 1.29f, 1.37f, "Inv. Mass (GeV)"};
    AxisSpec axisK0sMass = {400, 0.4f, 0.6f, "Inv. Mass (GeV)"};
    AxisSpec axisLambdaMass = {160, 1.08f, 1.16f, "Inv. Mass (GeV)"};
    if (cfgOutputCasc) {
      registry.add("InvMassXi_all", "", {HistType::kTHnSparseF, {cfgaxisPtXi, axisXiMass, cfgaxisEta, axisMultiplicity}});
      registry.add("InvMassOmega_all", "", {HistType::kTHnSparseF, {cfgaxisPtOmega, axisOmegaMass, cfgaxisEta, axisMultiplicity}});
      registry.add("InvMassXi", "", {HistType::kTHnSparseF, {cfgaxisPtXi, axisXiMass, cfgaxisEta, axisMultiplicity}});
      registry.add("InvMassOmega", "", {HistType::kTHnSparseF, {cfgaxisPtOmega, axisOmegaMass, cfgaxisEta, axisMultiplicity}});
    }
    if (cfgOutputV0) {
      registry.add("InvMassK0s_all", "", {HistType::kTHnSparseF, {cfgaxisPtK0s, axisK0sMass, cfgaxisEta, axisMultiplicity}});
      registry.add("InvMassLambda_all", "", {HistType::kTHnSparseF, {cfgaxisPtLambda, axisLambdaMass, cfgaxisEta, axisMultiplicity}});
      registry.add("InvMassK0s", "", {HistType::kTHnSparseF, {cfgaxisPtK0s, axisK0sMass, cfgaxisEta, axisMultiplicity}});
      registry.add("InvMassLambda", "", {HistType::kTHnSparseF, {cfgaxisPtLambda, axisLambdaMass, cfgaxisEta, axisMultiplicity}});
    }
    // for local density correction
    if (cfgOutputLocDenWeights) {
      registry.add("MC/densityMCGenK0s", "", {HistType::kTH3D, {cfgaxisPtK0s, cfgaxisNch, cfgaxisLocalDensity}});
      registry.add("MC/densityMCGenLambda", "", {HistType::kTH3D, {cfgaxisPtLambda, cfgaxisNch, cfgaxisLocalDensity}});
      registry.add("MC/densityMCGenXi", "", {HistType::kTH3D, {cfgaxisPtXi, cfgaxisNch, cfgaxisLocalDensity}});
      registry.add("MC/densityMCGenOmega", "", {HistType::kTH3D, {cfgaxisPtOmega, cfgaxisNch, cfgaxisLocalDensity}});
      registry.add("MC/densityMCRecK0s", "", {HistType::kTHnSparseF, {cfgaxisPtK0s, cfgaxisNch, cfgaxisLocalDensity, axisK0sMass}});
      registry.add("MC/densityMCRecLambda", "", {HistType::kTHnSparseF, {cfgaxisPtLambda, cfgaxisNch, cfgaxisLocalDensity, axisLambdaMass}});
      registry.add("MC/densityMCRecXi", "", {HistType::kTHnSparseF, {cfgaxisPtXi, cfgaxisNch, cfgaxisLocalDensity, axisXiMass}});
      registry.add("MC/densityMCRecOmega", "", {HistType::kTHnSparseF, {cfgaxisPtOmega, cfgaxisNch, cfgaxisLocalDensity, axisOmegaMass}});

      registry.add("MC/densityMCGenK0sMultTPC", "", {HistType::kTH3D, {cfgaxisPtK0s, cfgaxisNch, cfgaxisLocalDensity}});
      registry.add("MC/densityMCGenLambdaMultTPC", "", {HistType::kTH3D, {cfgaxisPtLambda, cfgaxisNch, cfgaxisLocalDensity}});
      registry.add("MC/densityMCGenXiMultTPC", "", {HistType::kTH3D, {cfgaxisPtXi, cfgaxisNch, cfgaxisLocalDensity}});
      registry.add("MC/densityMCGenOmegaMultTPC", "", {HistType::kTH3D, {cfgaxisPtOmega, cfgaxisNch, cfgaxisLocalDensity}});
      registry.add("MC/densityMCRecK0sMultTPC", "", {HistType::kTHnSparseF, {cfgaxisPtK0s, cfgaxisNch, cfgaxisLocalDensity, axisK0sMass}});
      registry.add("MC/densityMCRecLambdaMultTPC", "", {HistType::kTHnSparseF, {cfgaxisPtLambda, cfgaxisNch, cfgaxisLocalDensity, axisLambdaMass}});
      registry.add("MC/densityMCRecXiMultTPC", "", {HistType::kTHnSparseF, {cfgaxisPtXi, cfgaxisNch, cfgaxisLocalDensity, axisXiMass}});
      registry.add("MC/densityMCRecOmegaMultTPC", "", {HistType::kTHnSparseF, {cfgaxisPtOmega, cfgaxisNch, cfgaxisLocalDensity, axisOmegaMass}});
    }

    // Data
    fGFW->AddRegion("reffull", -0.8, 0.8, 1, 1); // ("name", etamin, etamax, ptbinnum, bitmask)eta region -0.8 to 0.8
    fGFW->AddRegion("refN10", -0.8, -0.4, 1, 1);
    fGFW->AddRegion("refP10", 0.4, 0.8, 1, 1);
    // POI
    fGFW->AddRegion("poiN10dpt", -0.8, -0.4, nPtBins, 32);
    fGFW->AddRegion("poiP10dpt", 0.4, 0.8, nPtBins, 32);
    fGFW->AddRegion("poifulldpt", -0.8, 0.8, nPtBins, 32);
    fGFW->AddRegion("poioldpt", -0.8, 0.8, nPtBins, 1);

    int nXiptMassBins = nXiPtBins * cfgmassbins[2];
    fGFW->AddRegion("poiXiPdpt", 0.4, 0.8, nXiptMassBins, 2);
    fGFW->AddRegion("poiXiNdpt", -0.8, -0.4, nXiptMassBins, 2);
    fGFW->AddRegion("poiXifulldpt", -0.8, 0.8, nXiptMassBins, 2);
    int nOmegaptMassBins = nXiPtBins * cfgmassbins[3];
    fGFW->AddRegion("poiOmegaPdpt", 0.4, 0.8, nOmegaptMassBins, 4);
    fGFW->AddRegion("poiOmegaNdpt", -0.8, -0.4, nOmegaptMassBins, 4);
    fGFW->AddRegion("poiOmegafulldpt", -0.8, 0.8, nOmegaptMassBins, 4);
    int nK0sptMassBins = nK0sPtBins * cfgmassbins[0];
    fGFW->AddRegion("poiK0sPdpt", 0.4, 0.8, nK0sptMassBins, 8);
    fGFW->AddRegion("poiK0sNdpt", -0.8, -0.4, nK0sptMassBins, 8);
    fGFW->AddRegion("poiK0sfulldpt", -0.8, 0.8, nK0sptMassBins, 8);
    int nLambdaptMassBins = nLambdaPtBins * cfgmassbins[1];
    fGFW->AddRegion("poiLambdaPdpt", 0.4, 0.8, nLambdaptMassBins, 16);
    fGFW->AddRegion("poiLambdaNdpt", -0.8, -0.4, nLambdaptMassBins, 16);
    fGFW->AddRegion("poiLambdafulldpt", -0.8, 0.8, nLambdaptMassBins, 16);
    // MC
    fGFW->AddRegion("refN10MC", -0.8, -0.4, 1, 64);
    fGFW->AddRegion("refP10MC", 0.4, 0.8, 1, 64);
    fGFW->AddRegion("poiXiPdptMC", 0.4, 0.8, nXiptMassBins, 128);
    fGFW->AddRegion("poiXiNdptMC", -0.8, -0.4, nXiptMassBins, 128);
    fGFW->AddRegion("poiOmegaPdptMC", 0.4, 0.8, nOmegaptMassBins, 256);
    fGFW->AddRegion("poiOmegaNdptMC", -0.8, -0.4, nOmegaptMassBins, 256);
    fGFW->AddRegion("poiK0sPdptMC", 0.4, 0.8, nK0sptMassBins, 512);
    fGFW->AddRegion("poiK0sNdptMC", -0.8, -0.4, nK0sptMassBins, 512);
    fGFW->AddRegion("poiLambdaPdptMC", 0.4, 0.8, nLambdaptMassBins, 1024);
    fGFW->AddRegion("poiLambdaNdptMC", -0.8, -0.4, nLambdaptMassBins, 1024);
    // pushback
    // Data
    // v2
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiP10dpt {2} refN10 {-2}", "Poi10Gap22dpta", kTRUE)); // 0
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiN10dpt {2} refP10 {-2}", "Poi10Gap22dptb", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poifulldpt reffull | poioldpt {2 2 -2 -2}", "Poi10Gap24dpt", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poifulldpt reffull | poioldpt {2 -2}", "PoiFull22dpt", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiXiPdpt {2} refN10 {-2}", "Xi10Gap22a", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiXiNdpt {2} refP10 {-2}", "Xi10Gap22b", kTRUE)); // 5
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiXifulldpt reffull {2 2 -2 -2}", "Xi10Gap24", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiXifulldpt {2} reffull {-2}", "XiFull22", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiOmegaPdpt {2} refN10 {-2}", "Omega10Gap22a", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiOmegaNdpt {2} refP10 {-2}", "Omega10Gap22b", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiOmegafulldpt reffull {2 2 -2 -2}", "Omega10Gap24", kTRUE)); // 10
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiOmegafulldpt {2} reffull {-2}", "OmegaFull22", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiK0sPdpt {2} refN10 {-2}", "K0short10Gap22a", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiK0sNdpt {2} refP10 {-2}", "K0short10Gap22b", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiK0sfulldpt reffull {2 2 -2 -2}", "K0short10Gap24", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiK0sfulldpt {2} reffull {-2}", "K0shortFull22", kTRUE)); // 15
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiLambdaPdpt {2} refN10 {-2}", "Lambda10Gap22a", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiLambdaNdpt {2} refP10 {-2}", "Lambda10Gap22b", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiLambdafulldpt reffull {2 2 -2 -2}", "LambdaFull24", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiLambdafulldpt {2} reffull {-2}", "LambdaFull22", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refP10 {2} refN10 {-2}", "Ref10Gap22a", kFALSE)); // 20
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("reffull reffull {2 2 -2 -2}", "Ref10Gap24", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("reffull reffull {2 -2}", "RefFull22", kFALSE));
    // v3
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiXiPdpt {3} refN10 {-3}", "Xi10Gap32a", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiXiNdpt {3} refP10 {-3}", "Xi10Gap32b", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiOmegaPdpt {3} refN10 {-3}", "Omega10Gap32a", kTRUE)); // 25
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiOmegaNdpt {3} refP10 {-3}", "Omega10Gap32b", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiK0sPdpt {3} refN10 {-3}", "K0short10Gap32a", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiK0sNdpt {3} refP10 {-3}", "K0short10Gap32b", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiLambdaPdpt {3} refN10 {-3}", "Lambda10Gap32a", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiLambdaNdpt {3} refP10 {-3}", "Lambda10Gap32b", kTRUE)); // 30
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refP10 {3} refN10 {-3}", "Ref10Gap32a", kFALSE));
    // MC
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiXiPdptMC {2} refN10MC {-2}", "MCXi10Gap22a", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiXiNdptMC {2} refP10MC {-2}", "MCXi10Gap22b", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiOmegaPdptMC {2} refN10MC {-2}", "MCOmega10Gap22a", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiOmegaNdptMC {2} refP10MC {-2}", "MCOmega10Gap22b", kTRUE)); // 35
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiK0sPdptMC {2} refN10MC {-2}", "MCK0s10Gap22a", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiK0sNdptMC {2} refP10MC {-2}", "MCK0s10Gap22b", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiLambdaPdptMC {2} refN10MC {-2}", "MCLambda10Gap22a", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiLambdaNdptMC {2} refP10MC {-2}", "MCLambda10Gap22b", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refP10MC {2} refN10MC {-2}", "MCRef10Gap22a", kFALSE)); // 40

    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiXiPdpt refP10 {2, 2} refN10 {-2 -2}", "Xi10Gap24a", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiXiNdpt refN10 {2, 2} refP10 {-2 -2}", "Xi10Gap24b", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiOmegaPdpt refP10 {2, 2} refN10 {-2 -2}", "Omega10Gap24a", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiOmegaNdpt refN10 {2, 2} refP10 {-2 -2}", "Omega10Gap24b", kTRUE));
    fGFW->CreateRegions(); // finalize the initialization

    // used for event selection
    fMultPVCutLow = new TF1("fMultPVCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x - 3.5*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)", 0, 100);
    fMultPVCutLow->SetParameters(cfgMultPVCutPara[0], cfgMultPVCutPara[1], cfgMultPVCutPara[2], cfgMultPVCutPara[3], cfgMultPVCutPara[4], cfgMultPVCutPara[5], cfgMultPVCutPara[6], cfgMultPVCutPara[7], cfgMultPVCutPara[8], cfgMultPVCutPara[9]);
    fMultPVCutHigh = new TF1("fMultPVCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x + 3.5*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)", 0, 100);
    fMultPVCutHigh->SetParameters(cfgMultPVCutPara[0], cfgMultPVCutPara[1], cfgMultPVCutPara[2], cfgMultPVCutPara[3], cfgMultPVCutPara[4], cfgMultPVCutPara[5], cfgMultPVCutPara[6], cfgMultPVCutPara[7], cfgMultPVCutPara[8], cfgMultPVCutPara[9]);

    fT0AV0AMean = new TF1("fT0AV0AMean", "[0]+[1]*x", 0, 200000);
    fT0AV0AMean->SetParameters(-1601.0581, 9.417652e-01);
    fT0AV0ASigma = new TF1("fT0AV0ASigma", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 200000);
    fT0AV0ASigma->SetParameters(463.4144, 6.796509e-02, -9.097136e-07, 7.971088e-12, -2.600581e-17);

    // fWeight output
    if (cfgOutputNUAWeights) {
      fWeightsREF->setPtBins(nPtBins, &(axisPt.binEdges)[0]);
      fWeightsREF->init(true, false);
      fWeightsK0s->setPtBins(nPtBins, &(axisPt.binEdges)[0]);
      fWeightsK0s->init(true, false);
      fWeightsLambda->setPtBins(nPtBins, &(axisPt.binEdges)[0]);
      fWeightsLambda->init(true, false);
      fWeightsXi->setPtBins(nPtBins, &(axisPt.binEdges)[0]);
      fWeightsXi->init(true, false);
      fWeightsOmega->setPtBins(nPtBins, &(axisPt.binEdges)[0]);
      fWeightsOmega->init(true, false);
    }
  }

  // input HIST("name")
  template <char... chars>
  void fillProfile(const GFW::CorrConfig& corrconf, const ConstStr<chars...>& tarName, const double& cent)
  {
    double dnx, val;
    dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
    if (dnx == 0)
      return;
    if (!corrconf.pTDif) {
      val = fGFW->Calculate(corrconf, 0, kFALSE).real() / dnx;
      if (std::fabs(val) < 1)
        registry.fill(tarName, cent, val, dnx);
      return;
    }
    return;
  }

  // input shared_ptr<TProfile>
  void fillProfile(const GFW::CorrConfig& corrconf, std::shared_ptr<TProfile> TProfile, const double& cent)
  {
    double dnx, val;
    dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
    if (dnx == 0)
      return;
    if (!corrconf.pTDif) {
      val = fGFW->Calculate(corrconf, 0, kFALSE).real() / dnx;
      if (std::fabs(val) < 1)
        TProfile->Fill(cent, val, dnx);
      return;
    }
    return;
  }

  template <char... chars>
  void fillProfilepT(const GFW::CorrConfig& corrconf, const ConstStr<chars...>& tarName, const int& ptbin, const double& cent)
  {
    float dnx = 0;
    float val = 0;
    dnx = fGFW->Calculate(corrconf, ptbin - 1, kTRUE).real();
    if (dnx == 0)
      return;
    val = fGFW->Calculate(corrconf, ptbin - 1, kFALSE).real() / dnx;
    if (std::fabs(val) < 1) {
      registry.fill(tarName, fPtAxis->GetBinCenter(ptbin), cent, val, dnx);
    }
    return;
  }

  template <char... chars>
  void fillProfilepTMC(const GFW::CorrConfig& corrconf, const ConstStr<chars...>& tarName, const int& ptbin, const int& PDGCode, const double& cent)
  {
    TAxis* fpt = nullptr;
    if (PDGCode == kXiMinus) {
      fpt = fXiPtAxis;
    } else if (PDGCode == kOmegaMinus) {
      fpt = fOmegaPtAxis;
    } else if (PDGCode == kK0Short) {
      fpt = fK0sPtAxis;
    } else if (PDGCode == kLambda0) {
      fpt = fLambdaPtAxis;
    } else {
      LOGF(error, "Error, please put in correct PDGCode of K0s, Lambda, Xi or Omega");
      return;
    }
    float dnx = 0;
    float val = 0;
    dnx = fGFW->Calculate(corrconf, ptbin - 1, kTRUE).real();
    if (dnx == 0)
      return;
    val = fGFW->Calculate(corrconf, ptbin - 1, kFALSE).real() / dnx;
    if (std::fabs(val) < 1) {
      registry.fill(tarName, fpt->GetBinCenter(ptbin), cent, val, dnx);
    }
    return;
  }

  // input HIST("name")
  template <char... chars>
  void fillProfilepTMass(const GFW::CorrConfig& corrconf, const ConstStr<chars...>& tarName, const int& ptbin, const int& PDGCode, const float& cent)
  {
    int nMassBins = 0;
    int nptbins = 0;
    TAxis* fMass = nullptr;
    TAxis* fpt = nullptr;
    if (PDGCode == kXiMinus) {
      nMassBins = cfgmassbins[2];
      nptbins = nXiPtBins;
      fpt = fXiPtAxis;
      fMass = fXiMass;
    } else if (PDGCode == kOmegaMinus) {
      nMassBins = cfgmassbins[3];
      nptbins = nOmegaPtBins;
      fpt = fOmegaPtAxis;
      fMass = fOmegaMass;
    } else if (PDGCode == kK0Short) {
      nMassBins = cfgmassbins[0];
      nptbins = nK0sPtBins;
      fpt = fK0sPtAxis;
      fMass = fK0sMass;
    } else if (PDGCode == kLambda0) {
      nMassBins = cfgmassbins[1];
      nptbins = nLambdaPtBins;
      fpt = fLambdaPtAxis;
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
      if (std::fabs(val) < 1) {
        registry.fill(tarName, fpt->GetBinCenter(ptbin), fMass->GetBinCenter(massbin), cent, val, dnx);
      }
    }
    return;
  }

  // input shared_ptr<TProfile3D>
  void fillProfilepTMass(const GFW::CorrConfig& corrconf, std::shared_ptr<TProfile3D> TProfile3D, const int& ptbin, const int& PDGCode, const float& cent)
  {
    int nMassBins = 0;
    int nptbins = 0;
    TAxis* fMass = nullptr;
    TAxis* fpt = nullptr;
    if (PDGCode == kXiMinus) {
      nMassBins = cfgmassbins[2];
      nptbins = nXiPtBins;
      fpt = fXiPtAxis;
      fMass = fXiMass;
    } else if (PDGCode == kOmegaMinus) {
      nMassBins = cfgmassbins[3];
      nptbins = nOmegaPtBins;
      fpt = fOmegaPtAxis;
      fMass = fOmegaMass;
    } else if (PDGCode == kK0Short) {
      nMassBins = cfgmassbins[0];
      nptbins = nK0sPtBins;
      fpt = fK0sPtAxis;
      fMass = fK0sMass;
    } else if (PDGCode == kLambda0) {
      nMassBins = cfgmassbins[1];
      nptbins = nLambdaPtBins;
      fpt = fLambdaPtAxis;
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
      if (std::fabs(val) < 1) {
        TProfile3D->Fill(fpt->GetBinCenter(ptbin), fMass->GetBinCenter(massbin), cent, val, dnx);
      }
    }
    return;
  }

  void loadCorrections(uint64_t timestamp)
  {
    if (correctionsLoaded)
      return;
    int nspecies = 5;
    if (cfgAcceptance.size() == static_cast<uint64_t>(nspecies)) {
      for (int i = 0; i <= nspecies - 1; i++) {
        mAcceptance.push_back(ccdb->getForTimeStamp<GFWWeights>(cfgAcceptance[i], timestamp));
      }
      if (mAcceptance.size() == static_cast<uint64_t>(nspecies))
        LOGF(info, "Loaded acceptance weights");
      else
        LOGF(warning, "Could not load acceptance weights");
    }

    if (cfgEfficiency.size() == 1)
      mEfficiency.push_back(ccdb->getForTimeStamp<TH1D>(cfgEfficiency[0], timestamp));

    correctionsLoaded = true;
  }

  template <typename TrackObject>
  bool setCurrentParticleWeights(float& weight_nue, float& weight_nua, TrackObject track, float vtxz, int ispecies)
  {
    float eff = 1.;
    int nspecies = 5;
    if (ispecies == 0 && cfgEfficiency.size() == 1)
      eff = mEfficiency[ispecies]->GetBinContent(mEfficiency[ispecies]->FindBin(track.pt()));
    if (eff <= 0)
      return false;
    weight_nue = 1. / eff;
    if (mAcceptance.size() == static_cast<uint64_t>(nspecies))
      weight_nua = mAcceptance[ispecies]->getNUA(track.phi(), track.eta(), vtxz);
    else
      weight_nua = 1;
    return true;
  }

  template <typename TrackObject>
  bool setCurrentLocalDensityWeights(float& weight_loc, TrackObject track, double locDensity, int ispecies)
  {
    auto cfgLocDenPara = (std::vector<std::vector<double>>){cfgLocDenParaK0s, cfgLocDenParaLambda, cfgLocDenParaXi, cfgLocDenParaOmega};
    int ptbin = fXiPtAxis->FindBin(track.pt());
    if (ptbin == 0 || ptbin == (fXiPtAxis->GetNbins() + 1)) {
      weight_loc = 1.0;
      return true;
    }
    double paraA = cfgLocDenPara[ispecies - 1][2 * ptbin - 2];
    double paraB = cfgLocDenPara[ispecies - 1][2 * ptbin - 1];
    double density = locDensity * 200 / (2 * cfgDeltaPhiLocDen + 1);
    double eff = std::exp(paraA * density + paraB);
    weight_loc = 1 / eff;
    return true;
  }

  // event selection
  template <typename TCollision>
  bool eventSelected(TCollision collision, const float centrality, float interactionRate = -1)
  {
    if (evtSeleOpts.cfgDoTVXinTRD.value && collision.alias_bit(kTVXinTRD)) {
      // TRD triggered
      return false;
    }
    registry.fill(HIST("hEventCount"), 2.5);
    if (evtSeleOpts.cfgDoNoTimeFrameBorder.value && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      // reject collisions close to Time Frame borders
      // https://its.cern.ch/jira/browse/O2-4623
      return false;
    }
    registry.fill(HIST("hEventCount"), 3.5);
    if (evtSeleOpts.cfgDoNoITSROFrameBorder.value && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      // reject events affected by the ITS ROF border
      // https://its.cern.ch/jira/browse/O2-4309
      return false;
    }
    registry.fill(HIST("hEventCount"), 4.5);
    if (evtSeleOpts.cfgDoNoSameBunchPileup.value && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      // rejects collisions which are associated with the same "found-by-T0" bunch crossing
      // https://indico.cern.ch/event/1396220/#1-event-selection-with-its-rof
      return false;
    }
    registry.fill(HIST("hEventCount"), 5.5);
    if (evtSeleOpts.cfgDoIsGoodZvtxFT0vsPV.value && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      // removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference
      // use this cut at low multiplicities with caution
      return false;
    }
    registry.fill(HIST("hEventCount"), 6.5);
    if (evtSeleOpts.cfgDoNoCollInTimeRangeStandard.value && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      // no collisions in specified time range
      return 0;
    }
    registry.fill(HIST("hEventCount"), 7.5);
    if (evtSeleOpts.cfgDoIsGoodITSLayersAll.value && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      // cut time intervals with dead ITS staves
      return 0;
    }
    registry.fill(HIST("hEventCount"), 8.5);
    float vtxz = -999;
    if (collision.numContrib() > 1) {
      vtxz = collision.posZ();
      float zRes = std::sqrt(collision.covZZ());
      double zResMin = 0.25;
      int numContMax = 20;
      if (zRes > zResMin && collision.numContrib() < numContMax)
        vtxz = -999;
    }
    auto multNTracksPV = collision.multNTracksPV();
    auto occupancy = collision.trackOccupancyInTimeRange();

    if (std::fabs(vtxz) > cfgCutVertex)
      return false;

    registry.fill(HIST("hNTracksPVvsCentrality"), multNTracksPV, centrality);
    if (evtSeleOpts.cfgDoMultPVCut.value) {
      if (multNTracksPV < fMultPVCutLow->Eval(centrality))
        return false;
      if (multNTracksPV > fMultPVCutHigh->Eval(centrality))
        return false;
    }

    registry.fill(HIST("hEventCount"), 9.5);

    if (occupancy > evtSeleOpts.cfgCutOccupancyHigh.value)
      return 0;
    registry.fill(HIST("hEventCount"), 10.5);

    // V0A T0A 5 sigma cut
    if (evtSeleOpts.cfgDoV0AT0Acut.value) {
      int nsigma = 5;
      if (std::fabs(collision.multFV0A() - fT0AV0AMean->Eval(collision.multFT0A())) > nsigma * fT0AV0ASigma->Eval(collision.multFT0A()))
        return 0;
    }
    registry.fill(HIST("hEventCount"), 11.5);

    registry.fill(HIST("hInteractionRate"), interactionRate);
    if (interactionRate > 0 && interactionRate < evtSeleOpts.cfgCutminIR.value)
      return false;
    registry.fill(HIST("hEventCount"), 12.5);
    if (interactionRate > evtSeleOpts.cfgCutmaxIR.value)
      return false;
    registry.fill(HIST("hEventCount"), 13.5);

    return true;
  }

  void processData(AodCollisions::iterator const& collision, aod::BCsWithTimestamps const&, AodTracks const& tracks, soa::Join<aod::CascDataExt, aod::CascTOFNSigmas> const& Cascades, aod::V0Datas const& V0s, DaughterTracks const&)
  {
    o2::aod::ITSResponse itsResponse;
    int nTot = tracks.size();
    float nMultTPC = collision.multTPC();
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int runNumber = bc.runNumber();
    double interactionRate = rateFetcher.fetch(ccdb.service, bc.timestamp(), runNumber, "ZNC hadronic") * 1.e-3;

    registry.fill(HIST("hEventCount"), 0.5);
    if (nTot < 1)
      return;
    fGFW->Clear();
    const auto cent = collision.centFT0C();
    if (!collision.sel8())
      return;
    registry.fill(HIST("hEventCount"), 1.5);

    if (!eventSelected(collision, cent, interactionRate))
      return;
    TH1D* hLocalDensity = new TH1D("hphi", "hphi", 400, -constants::math::TwoPI, constants::math::TwoPI);
    loadCorrections(bc.timestamp());
    float vtxz = collision.posZ();
    registry.fill(HIST("hVtxZ"), vtxz);
    registry.fill(HIST("hMult"), nTot);
    registry.fill(HIST("hMultTPC"), nMultTPC);
    registry.fill(HIST("hCent"), cent);

    float weff = 1;
    float wacc = 1;
    float wloc = 1;
    double nch = 0;
    // fill GFW ref flow
    for (const auto& track : tracks) {
      if (cfgDoAccEffCorr) {
        if (!setCurrentParticleWeights(weff, wacc, track, vtxz, 0))
          continue;
      }
      registry.fill(HIST("hPhi"), track.phi());
      registry.fill(HIST("hPhicorr"), track.phi(), wacc);
      registry.fill(HIST("hEta"), track.eta());
      registry.fill(HIST("hPt"), track.pt());
      int ptbin = fPtAxis->FindBin(track.pt()) - 1;
      if ((track.pt() > trkQualityOpts.cfgCutPtMin.value) && (track.pt() < trkQualityOpts.cfgCutPtMax.value)) {
        fGFW->Fill(track.eta(), ptbin, track.phi(), wacc * weff, 1); //(eta, ptbin, phi, wacc*weff, bitmask)
      }
      if ((track.pt() > trkQualityOpts.cfgCutPtPOIMin.value) && (track.pt() < trkQualityOpts.cfgCutPtPOIMax.value)) {
        fGFW->Fill(track.eta(), ptbin, track.phi(), wacc * weff, 32);
        if (cfgDoLocDenCorr) {
          hLocalDensity->Fill(track.phi(), wacc * weff);
          hLocalDensity->Fill(RecoDecay::constrainAngle(track.phi(), -constants::math::TwoPI), wacc * weff);
          nch += wacc * weff;
        }
      }
      if (cfgOutputNUAWeights)
        fWeightsREF->fill(track.phi(), track.eta(), vtxz, track.pt(), cent, 0);

      if (cfgOutputrunbyrun) {
        if (cfgOutputQA) {
          th1sList[runNumber][hPhi]->Fill(track.phi());
          th1sList[runNumber][hPhicorr]->Fill(track.phi(), wacc);
        }
        th3sList[runNumber][hPhiEtaVtxz]->Fill(track.phi(), track.eta(), vtxz);
      }
    }
    if (cfgDoLocDenCorr) {
      registry.fill(HIST("hCentvsNch"), cent, nch);
      registry.fill(HIST("hCentvsMultTPC"), cent, nMultTPC);
    }
    // fill GFW of V0 flow
    double lowpt = trkQualityOpts.cfgCutPtPIDDauMin.value;
    double bachPtcut = trkQualityOpts.cfgCutPtPIDbachMin.value;
    double dauLaPrPtcut = trkQualityOpts.cfgCutPtPIDdauLaPrMin.value;
    double dauLaPiPtcut = trkQualityOpts.cfgCutPtPIDdauLaPiMin.value;

    if (cfgOutputV0) {
      for (const auto& v0 : V0s) {
        auto v0posdau = v0.posTrack_as<DaughterTracks>();
        auto v0negdau = v0.negTrack_as<DaughterTracks>();
        // check tpc
        bool isK0s = false;
        bool isLambda = false;
        if (v0posdau.pt() < trkQualityOpts.cfgCutPtDauMin.value || v0posdau.pt() > trkQualityOpts.cfgCutPtDauMax.value)
          continue;
        if (v0negdau.pt() < trkQualityOpts.cfgCutPtDauMin.value || v0negdau.pt() > trkQualityOpts.cfgCutPtDauMax.value)
          continue;

        // fill QA
        if (cfgOutputQA) {
          registry.fill(HIST("QAhisto/V0/hqaarm_podobefore"), v0.alpha(), v0.qtarm());
        }
        // check daughter ITS, TPC and TOF
        // K0short
        if (v0.pt() > trkQualityOpts.cfgCutPtK0sMin.value && v0.pt() < trkQualityOpts.cfgCutPtK0sMax.value) {
          if (v0.qtarm() / std::fabs(v0.alpha()) > v0BuilderOpts.cfgv0_ArmPodocut.value &&
              std::fabs(v0.mK0Short() - o2::constants::physics::MassK0Short) < v0BuilderOpts.cfgv0_mk0swindow.value &&
              (std::fabs(v0posdau.tpcNSigmaPi()) < cfgNSigma[0] && std::fabs(v0negdau.tpcNSigmaPi()) < cfgNSigma[0]) &&
              ((std::fabs(v0posdau.tofNSigmaPi()) < cfgNSigma[3] || v0posdau.pt() < lowpt) && (std::fabs(v0negdau.tofNSigmaPi()) < cfgNSigma[3] || v0negdau.pt() < lowpt)) &&
              ((std::fabs(itsResponse.nSigmaITS<o2::track::PID::Pion>(v0posdau)) < cfgNSigma[6]) && v0posdau.pt() < lowpt) && ((std::fabs(itsResponse.nSigmaITS<o2::track::PID::Pion>(v0negdau)) < cfgNSigma[6]) && v0negdau.pt() < lowpt)) {
            registry.fill(HIST("InvMassK0s_all"), v0.pt(), v0.mK0Short(), v0.eta(), cent);
            isK0s = true;
            if (cfgOutputQA) {
              registry.fill(HIST("QAhisto/V0/hqaarm_podoafter"), v0.alpha(), v0.qtarm());
            }
          }
        }
        // Lambda and antiLambda
        if (v0.pt() > trkQualityOpts.cfgCutPtLambdaMin.value && v0.pt() < trkQualityOpts.cfgCutPtLambdaMax.value) {
          if (std::fabs(v0.mLambda() - o2::constants::physics::MassLambda) < v0BuilderOpts.cfgv0_mlambdawindow.value &&
              (std::fabs(v0posdau.tpcNSigmaPr()) < cfgNSigma[1] && std::fabs(v0negdau.tpcNSigmaPi()) < cfgNSigma[0]) &&
              ((std::fabs(v0posdau.tofNSigmaPr()) < cfgNSigma[4] || v0posdau.pt() < lowpt) && (std::fabs(v0negdau.tofNSigmaPi()) < cfgNSigma[3] || v0negdau.pt() < lowpt)) &&
              ((std::fabs(itsResponse.nSigmaITS<o2::track::PID::Proton>(v0posdau)) < cfgNSigma[7]) && v0posdau.pt() < lowpt) && ((std::fabs(itsResponse.nSigmaITS<o2::track::PID::Pion>(v0negdau)) < cfgNSigma[6]) && v0negdau.pt() < lowpt)) {
            registry.fill(HIST("InvMassLambda_all"), v0.pt(), v0.mLambda(), v0.eta(), cent);
            isLambda = true;
          } else if (std::fabs(v0.mLambda() - o2::constants::physics::MassLambda) < v0BuilderOpts.cfgv0_mlambdawindow.value &&
                     (std::fabs(v0negdau.tpcNSigmaPr()) < cfgNSigma[1] && std::fabs(v0posdau.tpcNSigmaPi()) < cfgNSigma[0]) &&
                     ((std::fabs(v0negdau.tofNSigmaPr()) < cfgNSigma[4] || v0negdau.pt() < lowpt) && (std::fabs(v0posdau.tofNSigmaPi()) < cfgNSigma[3] || v0posdau.pt() < lowpt)) &&
                     ((std::fabs(itsResponse.nSigmaITS<o2::track::PID::Proton>(v0posdau)) < cfgNSigma[7]) && v0posdau.pt() < lowpt) && ((std::fabs(itsResponse.nSigmaITS<o2::track::PID::Pion>(v0negdau)) < cfgNSigma[6]) && v0negdau.pt() < lowpt)) {
            registry.fill(HIST("InvMassLambda_all"), v0.pt(), v0.mLambda(), v0.eta(), cent);
            isLambda = true;
          }
        }
        // fill QA before cut
        if (cfgOutputQA) {
          registry.fill(HIST("QAhisto/V0/hqaV0radiusbefore"), v0.v0radius());
          registry.fill(HIST("QAhisto/V0/hqaV0cosPAbefore"), v0.v0cosPA());
          registry.fill(HIST("QAhisto/V0/hqadcaV0daubefore"), v0.dcaV0daughters());
          registry.fill(HIST("QAhisto/V0/hqadcapostoPVbefore"), v0.dcapostopv());
          registry.fill(HIST("QAhisto/V0/hqadcanegtoPVbefore"), v0.dcanegtopv());
        }
        if (!isK0s && !isLambda)
          continue;
        // track quality check
        if (v0posdau.itsNCls() <= trkQualityOpts.cfgITSNCls.value)
          continue;
        if (v0negdau.itsNCls() <= trkQualityOpts.cfgITSNCls.value)
          continue;
        if (v0posdau.tpcNClsFound() <= trkQualityOpts.cfgTPCNCls.value)
          continue;
        if (v0negdau.tpcNClsFound() <= trkQualityOpts.cfgTPCNCls.value)
          continue;
        if (v0posdau.tpcNClsCrossedRows() <= trkQualityOpts.cfgTPCCrossedRows.value)
          continue;
        if (v0negdau.tpcNClsCrossedRows() <= trkQualityOpts.cfgTPCCrossedRows.value)
          continue;
        if (trkQualityOpts.cfgCheckGlobalTrack.value) {
          if (!v0posdau.hasTPC() || !v0posdau.hasITS())
            continue;
          if (!v0negdau.hasTPC() || !v0negdau.hasITS())
            continue;
        }
        // // topological cut
        if (v0.v0radius() < v0BuilderOpts.cfgv0_radius.value)
          continue;
        if (v0.v0cosPA() < v0BuilderOpts.cfgv0_v0cospa.value)
          continue;
        if (v0.dcaV0daughters() > v0BuilderOpts.cfgv0_dcav0dau.value)
          continue;
        if (std::fabs(v0.dcapostopv()) < v0BuilderOpts.cfgv0_dcadautopv.value)
          continue;
        if (std::fabs(v0.dcanegtopv()) < v0BuilderOpts.cfgv0_dcadautopv.value)
          continue;
        if (isK0s && std::fabs(v0.mLambda() - o2::constants::physics::MassLambda0) < v0BuilderOpts.cfgv0_compmassrejLambda.value)
          isK0s = false;
        if (isLambda && std::fabs(v0.mK0Short() - o2::constants::physics::MassK0Short) < v0BuilderOpts.cfgv0_compmassrejK0s.value)
          isLambda = false;

        // fill QA after cut
        if (cfgOutputQA) {
          registry.fill(HIST("QAhisto/V0/hqaV0radiusafter"), v0.v0radius());
          registry.fill(HIST("QAhisto/V0/hqaV0cosPAafter"), v0.v0cosPA());
          registry.fill(HIST("QAhisto/V0/hqadcaV0dauafter"), v0.dcaV0daughters());
          registry.fill(HIST("QAhisto/V0/hqadcapostoPVafter"), v0.dcapostopv());
          registry.fill(HIST("QAhisto/V0/hqadcanegtoPVafter"), v0.dcanegtopv());
        }
        if (isK0s) {
          if (cfgDoAccEffCorr)
            setCurrentParticleWeights(weff, wacc, v0, vtxz, 1);
          if (cfgDoLocDenCorr) {
            int phibin = -999;
            phibin = hLocalDensity->FindBin(RecoDecay::constrainAngle(v0.phi(), -constants::math::PI));
            double density = hLocalDensity->Integral(phibin - cfgDeltaPhiLocDen, phibin + cfgDeltaPhiLocDen);
            setCurrentLocalDensityWeights(wloc, v0, density, 1);
            if (cfgOutputLocDenWeights) {
              registry.fill(HIST("MC/densityMCRecK0s"), v0.pt(), nch, density, v0.mK0Short());
              registry.fill(HIST("MC/densityMCRecK0sMultTPC"), v0.pt(), nMultTPC, density, v0.mK0Short());
            }
          }
          registry.fill(HIST("InvMassK0s"), v0.pt(), v0.mK0Short(), v0.eta(), cent);
          registry.fill(HIST("hPhiK0s"), v0.phi());
          registry.fill(HIST("hPhiK0scorr"), v0.phi(), wacc);
          fGFW->Fill(v0.eta(), fK0sPtAxis->FindBin(v0.pt()) - 1 + ((fK0sMass->FindBin(v0.mK0Short()) - 1) * nK0sPtBins), v0.phi(), wacc * weff * wloc, 8);
          if (cfgOutputNUAWeights)
            fWeightsK0s->fill(v0.phi(), v0.eta(), vtxz, v0.pt(), cent, 0);
          if (cfgOutputrunbyrun) {
            if (cfgOutputQA) {
              th1sList[runNumber][hPhiK0s]->Fill(v0.phi());
              th1sList[runNumber][hPhiK0scorr]->Fill(v0.phi(), wacc);
            }
            th3sList[runNumber][hPhiEtaVtxzK0s]->Fill(v0.phi(), v0.eta(), vtxz);
          }
        }
        if (isLambda) {
          if (cfgDoAccEffCorr)
            setCurrentParticleWeights(weff, wacc, v0, vtxz, 2);
          if (cfgDoLocDenCorr) {
            int phibin = -999;
            phibin = hLocalDensity->FindBin(RecoDecay::constrainAngle(v0.phi(), -constants::math::PI));
            double density = hLocalDensity->Integral(phibin - cfgDeltaPhiLocDen, phibin + cfgDeltaPhiLocDen);
            setCurrentLocalDensityWeights(wloc, v0, density, 2);
            if (cfgOutputLocDenWeights) {
              registry.fill(HIST("MC/densityMCRecLambda"), v0.pt(), nch, density, v0.mLambda());
              registry.fill(HIST("MC/densityMCRecLambdaMultTPC"), v0.pt(), nMultTPC, density, v0.mLambda());
            }
          }
          registry.fill(HIST("InvMassLambda"), v0.pt(), v0.mLambda(), v0.eta(), cent);
          registry.fill(HIST("hPhiLambda"), v0.phi());
          registry.fill(HIST("hPhiLambdacorr"), v0.phi(), wacc);
          fGFW->Fill(v0.eta(), fLambdaPtAxis->FindBin(v0.pt()) - 1 + ((fLambdaMass->FindBin(v0.mLambda()) - 1) * nLambdaPtBins), v0.phi(), wacc * weff * wloc, 16);
          if (cfgOutputNUAWeights)
            fWeightsLambda->fill(v0.phi(), v0.eta(), vtxz, v0.pt(), cent, 0);
          if (cfgOutputrunbyrun) {
            if (cfgOutputQA) {
              th1sList[runNumber][hPhiLambda]->Fill(v0.phi());
              th1sList[runNumber][hPhiLambdacorr]->Fill(v0.phi(), wacc);
            }
            th3sList[runNumber][hPhiEtaVtxzLambda]->Fill(v0.phi(), v0.eta(), vtxz);
          }
        }
      }
    }

    // fill GFW of casc flow
    if (cfgOutputCasc) {
      for (const auto& casc : Cascades) {
        auto bachelor = casc.bachelor_as<DaughterTracks>();
        auto posdau = casc.posTrack_as<DaughterTracks>();
        auto negdau = casc.negTrack_as<DaughterTracks>();
        // check TPC
        bool isOmega = false;
        bool isXi = false;

        if (bachelor.pt() < trkQualityOpts.cfgCutPtDauMin.value || bachelor.pt() > trkQualityOpts.cfgCutPtDauMax.value)
          continue;
        if (posdau.pt() < trkQualityOpts.cfgCutPtDauMin.value || posdau.pt() > trkQualityOpts.cfgCutPtDauMax.value)
          continue;
        if (negdau.pt() < trkQualityOpts.cfgCutPtDauMin.value || negdau.pt() > trkQualityOpts.cfgCutPtDauMax.value)
          continue;

        // Omega and antiOmega
        if (casc.pt() > trkQualityOpts.cfgCutPtOmegaMin.value && casc.pt() < trkQualityOpts.cfgCutPtOmegaMax.value) {
          if (casc.sign() < 0 && std::fabs(casc.yOmega()) < cfgCasc_rapidity &&
              (std::fabs(bachelor.tpcNSigmaKa()) < cfgNSigma[2] && std::fabs(posdau.tpcNSigmaPr()) < cfgNSigma[1] && std::fabs(negdau.tpcNSigmaPi()) < cfgNSigma[0]) &&
              ((std::fabs(casc.tofNSigmaOmKa()) < cfgNSigma[5] || bachelor.pt() < bachPtcut) && (std::fabs(casc.tofNSigmaOmLaPr()) < cfgNSigma[4] || posdau.pt() < dauLaPrPtcut) && (std::fabs(casc.tofNSigmaOmLaPi()) < cfgNSigma[3] || negdau.pt() < dauLaPiPtcut)) &&
              ((std::fabs(itsResponse.nSigmaITS<o2::track::PID::Kaon>(bachelor)) < cfgNSigma[8]) || bachelor.pt() < bachPtcut) && ((std::fabs(itsResponse.nSigmaITS<o2::track::PID::Proton>(posdau)) < cfgNSigma[7]) || posdau.pt() < dauLaPrPtcut) && ((std::fabs(itsResponse.nSigmaITS<o2::track::PID::Pion>(negdau)) < cfgNSigma[6]) || negdau.pt() < dauLaPiPtcut)) {
            registry.fill(HIST("InvMassOmega_all"), casc.pt(), casc.mOmega(), casc.eta(), cent);
            isOmega = true;
          } else if (casc.sign() > 0 && std::fabs(casc.yOmega()) < cfgCasc_rapidity &&
                     (std::fabs(bachelor.tpcNSigmaKa()) < cfgNSigma[2] && std::fabs(negdau.tpcNSigmaPr()) < cfgNSigma[1] && std::fabs(posdau.tpcNSigmaPi()) < cfgNSigma[0]) &&
                     ((std::fabs(casc.tofNSigmaOmKa()) < cfgNSigma[5] || bachelor.pt() < bachPtcut) && (std::fabs(casc.tofNSigmaOmLaPr()) < cfgNSigma[4] || negdau.pt() < dauLaPrPtcut) && (std::fabs(casc.tofNSigmaOmLaPi()) < cfgNSigma[3] || posdau.pt() < dauLaPiPtcut)) &&
                     ((std::fabs(itsResponse.nSigmaITS<o2::track::PID::Kaon>(bachelor)) < cfgNSigma[8]) || bachelor.pt() < bachPtcut) && ((std::fabs(itsResponse.nSigmaITS<o2::track::PID::Proton>(negdau)) < cfgNSigma[7]) || negdau.pt() < dauLaPrPtcut) && ((std::fabs(itsResponse.nSigmaITS<o2::track::PID::Pion>(posdau)) < cfgNSigma[6]) || posdau.pt() < dauLaPiPtcut)) {
            registry.fill(HIST("InvMassOmega_all"), casc.pt(), casc.mOmega(), casc.eta(), cent);
            isOmega = true;
          }
        }
        // Xi and antiXi
        if (casc.pt() > trkQualityOpts.cfgCutPtXiMin.value && casc.pt() < trkQualityOpts.cfgCutPtXiMax.value) {
          if (casc.sign() < 0 && std::fabs(casc.yXi()) < cfgCasc_rapidity &&
              (std::fabs(bachelor.tpcNSigmaPi()) < cfgNSigma[0] && std::fabs(posdau.tpcNSigmaPr()) < cfgNSigma[1] && std::fabs(negdau.tpcNSigmaPi()) < cfgNSigma[0]) &&
              ((std::fabs(casc.tofNSigmaXiPi()) < cfgNSigma[3] || bachelor.pt() < bachPtcut) && (std::fabs(casc.tofNSigmaXiLaPr()) < cfgNSigma[4] || posdau.pt() < dauLaPrPtcut) && (std::fabs(casc.tofNSigmaXiLaPi()) < cfgNSigma[3] || negdau.pt() < dauLaPiPtcut)) &&
              ((std::fabs(itsResponse.nSigmaITS<o2::track::PID::Pion>(bachelor)) < cfgNSigma[6]) || bachelor.pt() < bachPtcut) && ((std::fabs(itsResponse.nSigmaITS<o2::track::PID::Proton>(posdau)) < cfgNSigma[7]) || posdau.pt() < dauLaPrPtcut) && ((std::fabs(itsResponse.nSigmaITS<o2::track::PID::Pion>(negdau)) < cfgNSigma[6]) || negdau.pt() < dauLaPiPtcut)) {
            registry.fill(HIST("InvMassXi_all"), casc.pt(), casc.mXi(), casc.eta(), cent);
            isXi = true;
          } else if (casc.sign() > 0 && std::fabs(casc.yXi()) < cfgCasc_rapidity &&
                     (std::fabs(bachelor.tpcNSigmaPi()) < cfgNSigma[0] && std::fabs(negdau.tpcNSigmaPr()) < cfgNSigma[1] && std::fabs(posdau.tpcNSigmaPi()) < cfgNSigma[0]) &&
                     ((std::fabs(casc.tofNSigmaXiPi()) < cfgNSigma[3] || bachelor.pt() < bachPtcut) && (std::fabs(casc.tofNSigmaXiLaPr()) < cfgNSigma[4] || negdau.pt() < dauLaPrPtcut) && (std::fabs(casc.tofNSigmaXiLaPi()) < cfgNSigma[3] || posdau.pt() < dauLaPiPtcut)) &&
                     ((std::fabs(itsResponse.nSigmaITS<o2::track::PID::Pion>(bachelor)) < cfgNSigma[6]) || bachelor.pt() < bachPtcut) && ((std::fabs(itsResponse.nSigmaITS<o2::track::PID::Proton>(negdau)) < cfgNSigma[7]) || negdau.pt() < dauLaPrPtcut) && ((std::fabs(itsResponse.nSigmaITS<o2::track::PID::Pion>(posdau)) < cfgNSigma[6]) || posdau.pt() < dauLaPiPtcut)) {
            registry.fill(HIST("InvMassXi_all"), casc.pt(), casc.mXi(), casc.eta(), cent);
            isXi = true;
          }
        }
        // fill QA
        if (cfgOutputQA) {
          if (isXi) {
            registry.fill(HIST("QAhisto/Xi/hqaCascRadiusbefore"), casc.cascradius());
            registry.fill(HIST("QAhisto/Xi/hqaCasccosPAbefore"), casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
            registry.fill(HIST("QAhisto/Xi/hqaCascV0cosPAbefore"), casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
            registry.fill(HIST("QAhisto/Xi/hqadcaCascV0toPVbefore"), casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()));
            registry.fill(HIST("QAhisto/Xi/hqadcaCascBachtoPVbefore"), casc.dcabachtopv());
            registry.fill(HIST("QAhisto/Xi/hqadcaCascdaubefore"), casc.dcacascdaughters());
            registry.fill(HIST("QAhisto/Xi/hqadcaCascV0daubefore"), casc.dcaV0daughters());
          }
          if (isOmega) {
            registry.fill(HIST("QAhisto/Omega/hqaCascRadiusbefore"), casc.cascradius());
            registry.fill(HIST("QAhisto/Omega/hqaCasccosPAbefore"), casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
            registry.fill(HIST("QAhisto/Omega/hqaCascV0cosPAbefore"), casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
            registry.fill(HIST("QAhisto/Omega/hqadcaCascV0toPVbefore"), casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()));
            registry.fill(HIST("QAhisto/Omega/hqadcaCascBachtoPVbefore"), casc.dcabachtopv());
            registry.fill(HIST("QAhisto/Omega/hqadcaCascdaubefore"), casc.dcacascdaughters());
            registry.fill(HIST("QAhisto/Omega/hqadcaCascV0daubefore"), casc.dcaV0daughters());
          }
        }

        if (!isXi && !isOmega)
          continue;
        // // topological cut
        if (casc.cascradius() < cascBuilderOpts.cfgcasc_radius.value)
          continue;
        if (casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) < cascBuilderOpts.cfgcasc_casccospa.value)
          continue;
        if (casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < cascBuilderOpts.cfgcasc_v0cospa.value)
          continue;
        if (std::fabs(casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ())) < cascBuilderOpts.cfgcasc_dcav0topv.value)
          continue;
        if (std::fabs(casc.dcabachtopv()) < cascBuilderOpts.cfgcasc_dcabachtopv.value)
          continue;
        if (casc.dcacascdaughters() > cascBuilderOpts.cfgcasc_dcacascdau.value)
          continue;
        if (casc.dcaV0daughters() > cascBuilderOpts.cfgcasc_dcav0dau.value)
          continue;
        if (std::fabs(casc.mLambda() - o2::constants::physics::MassLambda0) > cascBuilderOpts.cfgcasc_mlambdawindow.value)
          continue;
        // track quality check
        if (bachelor.itsNCls() <= trkQualityOpts.cfgITSNCls.value)
          continue;
        if (posdau.itsNCls() <= trkQualityOpts.cfgITSNCls.value)
          continue;
        if (negdau.itsNCls() <= trkQualityOpts.cfgITSNCls.value)
          continue;
        if (bachelor.tpcNClsFound() <= trkQualityOpts.cfgTPCNCls.value)
          continue;
        if (posdau.tpcNClsFound() <= trkQualityOpts.cfgTPCNCls.value)
          continue;
        if (negdau.tpcNClsFound() <= trkQualityOpts.cfgTPCNCls.value)
          continue;
        if (bachelor.tpcNClsCrossedRows() <= trkQualityOpts.cfgTPCCrossedRows.value)
          continue;
        if (posdau.tpcNClsCrossedRows() <= trkQualityOpts.cfgTPCCrossedRows.value)
          continue;
        if (negdau.tpcNClsCrossedRows() <= trkQualityOpts.cfgTPCCrossedRows.value)
          continue;
        if (trkQualityOpts.cfgCheckGlobalTrack.value) {
          if (!bachelor.hasTPC() || !bachelor.hasITS())
            continue;
          if (!posdau.hasTPC() || !posdau.hasITS())
            continue;
          if (!negdau.hasTPC() || !negdau.hasITS())
            continue;
        }
        if (isXi && std::fabs(casc.mOmega() - o2::constants::physics::MassOmegaMinus) < cascBuilderOpts.cfgcasc_compmassrej.value) {
          isXi = false;
        }
        if (isOmega && std::fabs(casc.mXi() - o2::constants::physics::MassXiMinus) < cascBuilderOpts.cfgcasc_compmassrej.value) {
          isOmega = false;
        }
        // fill QA
        if (cfgOutputQA) {
          if (isXi) {
            registry.fill(HIST("QAhisto/Xi/hqaCascRadiusafter"), casc.cascradius());
            registry.fill(HIST("QAhisto/Xi/hqaCasccosPAafter"), casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
            registry.fill(HIST("QAhisto/Xi/hqaCascV0cosPAafter"), casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
            registry.fill(HIST("QAhisto/Xi/hqadcaCascV0toPVafter"), casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()));
            registry.fill(HIST("QAhisto/Xi/hqadcaCascBachtoPVafter"), casc.dcabachtopv());
            registry.fill(HIST("QAhisto/Xi/hqadcaCascdauafter"), casc.dcacascdaughters());
            registry.fill(HIST("QAhisto/Xi/hqadcaCascV0dauafter"), casc.dcaV0daughters());
          }
          if (isOmega) {
            registry.fill(HIST("QAhisto/Omega/hqaCascRadiusafter"), casc.cascradius());
            registry.fill(HIST("QAhisto/Omega/hqaCasccosPAafter"), casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
            registry.fill(HIST("QAhisto/Omega/hqaCascV0cosPAafter"), casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
            registry.fill(HIST("QAhisto/Omega/hqadcaCascV0toPVafter"), casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()));
            registry.fill(HIST("QAhisto/Omega/hqadcaCascBachtoPVafter"), casc.dcabachtopv());
            registry.fill(HIST("QAhisto/Omega/hqadcaCascdauafter"), casc.dcacascdaughters());
            registry.fill(HIST("QAhisto/Omega/hqadcaCascV0dauafter"), casc.dcaV0daughters());
          }
        }

        if (isOmega) {
          if (cfgDoAccEffCorr) {
            setCurrentParticleWeights(weff, wacc, casc, vtxz, 4);
          }
          if (cfgDoLocDenCorr) {
            int phibin = -999;
            phibin = hLocalDensity->FindBin(RecoDecay::constrainAngle(casc.phi(), -constants::math::PI));
            double density = hLocalDensity->Integral(phibin - cfgDeltaPhiLocDen, phibin + cfgDeltaPhiLocDen);
            setCurrentLocalDensityWeights(wloc, casc, density, 4);
            if (cfgOutputLocDenWeights) {
              registry.fill(HIST("MC/densityMCRecOmega"), casc.pt(), nch, density, casc.mOmega());
              registry.fill(HIST("MC/densityMCRecOmegaMultTPC"), casc.pt(), nMultTPC, density, casc.mOmega());
            }
          }
          registry.fill(HIST("hPhiOmega"), casc.phi());
          registry.fill(HIST("hPhiOmegacorr"), casc.phi(), wacc);
          registry.fill(HIST("InvMassOmega"), casc.pt(), casc.mOmega(), casc.eta(), cent);
          fGFW->Fill(casc.eta(), fOmegaPtAxis->FindBin(casc.pt()) - 1 + ((fOmegaMass->FindBin(casc.mOmega()) - 1) * nOmegaPtBins), casc.phi(), wacc * weff * wloc, 4);

          if (cfgOutputNUAWeights)
            fWeightsOmega->fill(casc.phi(), casc.eta(), vtxz, casc.pt(), cent, 0);
          if (cfgOutputrunbyrun) {
            if (cfgOutputQA) {
              th1sList[runNumber][hPhiOmega]->Fill(casc.phi());
              th1sList[runNumber][hPhiOmegacorr]->Fill(casc.phi(), wacc);
            }
            th3sList[runNumber][hPhiEtaVtxzOmega]->Fill(casc.phi(), casc.eta(), vtxz);
          }
        }
        if (isXi) {
          if (cfgDoAccEffCorr) {
            setCurrentParticleWeights(weff, wacc, casc, vtxz, 3);
          }
          if (cfgDoLocDenCorr) {
            int phibin = -999;
            phibin = hLocalDensity->FindBin(RecoDecay::constrainAngle(casc.phi(), -constants::math::PI));
            double density = hLocalDensity->Integral(phibin - cfgDeltaPhiLocDen, phibin + cfgDeltaPhiLocDen);
            setCurrentLocalDensityWeights(wloc, casc, density, 3);
            if (cfgOutputLocDenWeights) {
              registry.fill(HIST("MC/densityMCRecXi"), casc.pt(), nch, density, casc.mXi());
              registry.fill(HIST("MC/densityMCRecXiMultTPC"), casc.pt(), nMultTPC, density, casc.mXi());
            }
          }
          registry.fill(HIST("hPhiXi"), casc.phi());
          registry.fill(HIST("hPhiXicorr"), casc.phi(), wacc);
          registry.fill(HIST("InvMassXi"), casc.pt(), casc.mXi(), casc.eta(), cent);
          fGFW->Fill(casc.eta(), fXiPtAxis->FindBin(casc.pt()) - 1 + ((fXiMass->FindBin(casc.mXi()) - 1) * nXiPtBins), casc.phi(), wacc * weff * wloc, 2);

          if (cfgOutputNUAWeights)
            fWeightsXi->fill(casc.phi(), casc.eta(), vtxz, casc.pt(), cent, 0);
          if (cfgOutputrunbyrun) {
            if (cfgOutputQA) {
              th1sList[runNumber][hPhiXi]->Fill(casc.phi());
              th1sList[runNumber][hPhiXicorr]->Fill(casc.phi(), wacc);
            }
            th3sList[runNumber][hPhiEtaVtxzXi]->Fill(casc.phi(), casc.eta(), vtxz);
          }
        }
      }
    }
    delete hLocalDensity;
    // Filling cumulant with ROOT TProfile and loop for all ptBins
    fillProfile(corrconfigs.at(20), HIST("c22"), cent);
    fillProfile(corrconfigs.at(21), HIST("c24"), cent);
    fillProfile(corrconfigs.at(22), HIST("c22Full"), cent);
    fillProfile(corrconfigs.at(31), HIST("c32"), cent);
    for (int i = 1; i <= nPtBins; i++) {
      fillProfilepT(corrconfigs.at(0), HIST("c22dpt"), i, cent);
      fillProfilepT(corrconfigs.at(1), HIST("c22dpt"), i, cent);
      fillProfilepT(corrconfigs.at(2), HIST("c24dpt"), i, cent);
      fillProfilepT(corrconfigs.at(3), HIST("c22Fulldpt"), i, cent);
    }
    if (cfgOutputV0) {
      for (int i = 1; i <= nK0sPtBins; i++) {
        fillProfilepTMass(corrconfigs.at(12), HIST("K0sc22dpt"), i, kK0Short, cent);
        fillProfilepTMass(corrconfigs.at(13), HIST("K0sc22dpt"), i, kK0Short, cent);
        fillProfilepTMass(corrconfigs.at(14), HIST("K0sc24dpt"), i, kK0Short, cent);
        fillProfilepTMass(corrconfigs.at(15), HIST("K0sc22Fulldpt"), i, kK0Short, cent);
        fillProfilepTMass(corrconfigs.at(27), HIST("K0sc32dpt"), i, kK0Short, cent);
        fillProfilepTMass(corrconfigs.at(28), HIST("K0sc32dpt"), i, kK0Short, cent);
      }
      for (int i = 1; i <= nLambdaPtBins; i++) {
        fillProfilepTMass(corrconfigs.at(16), HIST("Lambdac22dpt"), i, kLambda0, cent);
        fillProfilepTMass(corrconfigs.at(17), HIST("Lambdac22dpt"), i, kLambda0, cent);
        fillProfilepTMass(corrconfigs.at(18), HIST("Lambdac24dpt"), i, kLambda0, cent);
        fillProfilepTMass(corrconfigs.at(19), HIST("Lambdac22Fulldpt"), i, kLambda0, cent);
        fillProfilepTMass(corrconfigs.at(29), HIST("Lambdac32dpt"), i, kLambda0, cent);
        fillProfilepTMass(corrconfigs.at(30), HIST("Lambdac32dpt"), i, kLambda0, cent);
      }
    }
    if (cfgOutputCasc) {
      for (int i = 1; i <= nXiPtBins; i++) {
        fillProfilepTMass(corrconfigs.at(4), HIST("Xic22dpt"), i, kXiMinus, cent);
        fillProfilepTMass(corrconfigs.at(5), HIST("Xic22dpt"), i, kXiMinus, cent);
        fillProfilepTMass(corrconfigs.at(6), HIST("Xic24dpt"), i, kXiMinus, cent);
        fillProfilepTMass(corrconfigs.at(7), HIST("Xic22Fulldpt"), i, kXiMinus, cent);
        fillProfilepTMass(corrconfigs.at(23), HIST("Xic32dpt"), i, kXiMinus, cent);
        fillProfilepTMass(corrconfigs.at(24), HIST("Xic32dpt"), i, kXiMinus, cent);

        fillProfilepTMass(corrconfigs.at(41), HIST("Xic24_gapdpt"), i, kXiMinus, cent);
        fillProfilepTMass(corrconfigs.at(42), HIST("Xic24_gapdpt"), i, kXiMinus, cent);
      }
      for (int i = 1; i <= nOmegaPtBins; i++) {
        fillProfilepTMass(corrconfigs.at(8), HIST("Omegac22dpt"), i, kOmegaMinus, cent);
        fillProfilepTMass(corrconfigs.at(9), HIST("Omegac22dpt"), i, kOmegaMinus, cent);
        fillProfilepTMass(corrconfigs.at(10), HIST("Omegac24dpt"), i, kOmegaMinus, cent);
        fillProfilepTMass(corrconfigs.at(11), HIST("Omegac22Fulldpt"), i, kOmegaMinus, cent);
        fillProfilepTMass(corrconfigs.at(25), HIST("Omegac32dpt"), i, kOmegaMinus, cent);
        fillProfilepTMass(corrconfigs.at(26), HIST("Omegac32dpt"), i, kOmegaMinus, cent);

        fillProfilepTMass(corrconfigs.at(43), HIST("Omegac24_gapdpt"), i, kOmegaMinus, cent);
        fillProfilepTMass(corrconfigs.at(44), HIST("Omegac24_gapdpt"), i, kOmegaMinus, cent);
      }
    }
    // Fill subevents flow
    if (cfgDoJackknife) {
      TRandom3* fRdm = new TRandom3(0);
      int nsubevent = 10;
      double eventrdm = nsubevent * fRdm->Rndm();
      for (int j = 1; j <= nsubevent; j++) {
        if (eventrdm > (j - 1) && eventrdm < j)
          continue;
        fillProfile(corrconfigs.at(20), refc22[j - 1], cent);
        fillProfile(corrconfigs.at(21), refc24[j - 1], cent);
        fillProfile(corrconfigs.at(22), refc22Full[j - 1], cent);
        fillProfile(corrconfigs.at(31), refc32[j - 1], cent);
        if (cfgOutputV0) {
          for (int i = 1; i <= nK0sPtBins; i++) {
            fillProfilepTMass(corrconfigs.at(12), k0sc22[j - 1], i, kK0Short, cent);
            fillProfilepTMass(corrconfigs.at(13), k0sc22[j - 1], i, kK0Short, cent);
            fillProfilepTMass(corrconfigs.at(14), k0sc24[j - 1], i, kK0Short, cent);
            fillProfilepTMass(corrconfigs.at(15), k0sc22Full[j - 1], i, kK0Short, cent);
            fillProfilepTMass(corrconfigs.at(27), k0sc32[j - 1], i, kK0Short, cent);
            fillProfilepTMass(corrconfigs.at(28), k0sc32[j - 1], i, kK0Short, cent);
          }
          for (int i = 1; i <= nLambdaPtBins; i++) {
            fillProfilepTMass(corrconfigs.at(16), lambdac22[j - 1], i, kLambda0, cent);
            fillProfilepTMass(corrconfigs.at(17), lambdac22[j - 1], i, kLambda0, cent);
            fillProfilepTMass(corrconfigs.at(18), lambdac24[j - 1], i, kLambda0, cent);
            fillProfilepTMass(corrconfigs.at(19), lambdac22Full[j - 1], i, kLambda0, cent);
            fillProfilepTMass(corrconfigs.at(29), lambdac32[j - 1], i, kLambda0, cent);
            fillProfilepTMass(corrconfigs.at(30), lambdac32[j - 1], i, kLambda0, cent);
          }
        }
        if (cfgOutputCasc) {
          for (int i = 1; i <= nXiPtBins; i++) {
            fillProfilepTMass(corrconfigs.at(4), xic22[j - 1], i, kXiMinus, cent);
            fillProfilepTMass(corrconfigs.at(5), xic22[j - 1], i, kXiMinus, cent);
            fillProfilepTMass(corrconfigs.at(6), xic24[j - 1], i, kXiMinus, cent);
            fillProfilepTMass(corrconfigs.at(7), xic22Full[j - 1], i, kXiMinus, cent);
            fillProfilepTMass(corrconfigs.at(23), xic32[j - 1], i, kXiMinus, cent);
            fillProfilepTMass(corrconfigs.at(24), xic32[j - 1], i, kXiMinus, cent);
          }
          for (int i = 1; i <= nOmegaPtBins; i++) {
            fillProfilepTMass(corrconfigs.at(8), omegac22[j - 1], i, kOmegaMinus, cent);
            fillProfilepTMass(corrconfigs.at(9), omegac22[j - 1], i, kOmegaMinus, cent);
            fillProfilepTMass(corrconfigs.at(10), omegac24[j - 1], i, kOmegaMinus, cent);
            fillProfilepTMass(corrconfigs.at(11), omegac22Full[j - 1], i, kOmegaMinus, cent);
            fillProfilepTMass(corrconfigs.at(25), omegac32[j - 1], i, kOmegaMinus, cent);
            fillProfilepTMass(corrconfigs.at(26), omegac32[j - 1], i, kOmegaMinus, cent);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(FlowGfwOmegaXi, processData, "", true);

  void processMCGen(aod::McCollisions::iterator const&, soa::Join<aod::McParticles, aod::ParticlesToTracks> const& tracksGen, soa::SmallGroups<soa::Join<aod::McCollisionLabels, AodCollisions>> const& collisionsRec, AodTracks const&)
  {
    fGFW->Clear();
    int nch = 0;
    float nMultTPC = 0;
    double cent = -1;
    TH1D* hLocalDensity = new TH1D("hphi", "hphi", 400, -constants::math::TwoPI, constants::math::TwoPI);
    for (const auto& collision : collisionsRec) {
      if (!collision.sel8())
        return;
      if (!eventSelected(collision, cent))
        return;
      cent = collision.centFT0C();
      nMultTPC = collision.multTPC();
      registry.fill(HIST("MC/hCentvsMultTPCMC"), cent, nMultTPC);
    }
    if (cent < 0)
      return;

    for (auto const& mcParticle : tracksGen) {
      if (!mcParticle.isPhysicalPrimary())
        continue;

      if (mcParticle.has_tracks()) {
        auto const& tracks = mcParticle.tracks_as<AodTracks>();
        for (const auto& track : tracks) {
          if (std::fabs(track.eta()) > trkQualityOpts.cfgCutEta.value) {
            continue;
          }
          if (!(track.isGlobalTrack())) {
            continue;
          }
          if (track.tpcChi2NCl() > cfgCutChi2prTPCcls) {
            continue;
          }
          int ptbin = fPtAxis->FindBin(mcParticle.pt()) - 1;
          if ((mcParticle.pt() > trkQualityOpts.cfgCutPtMin.value) && (mcParticle.pt() < trkQualityOpts.cfgCutPtMax.value)) {
            fGFW->Fill(mcParticle.eta(), ptbin, mcParticle.phi(), 1, 64); //(eta, ptbin, phi, wacc*weff, bitmask)
          }
          if ((mcParticle.pt() > trkQualityOpts.cfgCutPtPOIMin.value) && (mcParticle.pt() < trkQualityOpts.cfgCutPtPOIMax.value)) {
            hLocalDensity->Fill(mcParticle.phi(), 1);
            hLocalDensity->Fill(RecoDecay::constrainAngle(mcParticle.phi(), -constants::math::TwoPI), 1);
            nch++;
          }
        }
      }
    }
    registry.fill(HIST("MC/hCentvsNchMC"), cent, nch);

    for (const auto& straGen : tracksGen) {
      if (!straGen.isPhysicalPrimary())
        continue;
      int pdgCode = std::abs(straGen.pdgCode());
      if (pdgCode != PDG_t::kXiMinus && pdgCode != PDG_t::kOmegaMinus && pdgCode != PDG_t::kK0Short && pdgCode != PDG_t::kLambda0)
        continue;
      if (std::fabs(straGen.eta()) > trkQualityOpts.cfgCutEta.value)
        continue;

      if (pdgCode == PDG_t::kXiMinus) {
        int phibin = hLocalDensity->FindBin(RecoDecay::constrainAngle(straGen.phi(), -constants::math::PI));
        double density = hLocalDensity->Integral(phibin - cfgDeltaPhiLocDen, phibin + cfgDeltaPhiLocDen);
        if (cfgOutputLocDenWeights) {
          registry.fill(HIST("MC/densityMCGenXi"), straGen.pt(), nch, density);
          registry.fill(HIST("MC/densityMCGenXiMultTPC"), straGen.pt(), nMultTPC, density);
        }
        fGFW->Fill(straGen.eta(), fXiPtAxis->FindBin(straGen.pt()) - 1, straGen.phi(), 1, 128);
      }
      if (pdgCode == PDG_t::kOmegaMinus) {
        int phibin = hLocalDensity->FindBin(RecoDecay::constrainAngle(straGen.phi(), -constants::math::PI));
        double density = hLocalDensity->Integral(phibin - cfgDeltaPhiLocDen, phibin + cfgDeltaPhiLocDen);
        if (cfgOutputLocDenWeights) {
          registry.fill(HIST("MC/densityMCGenOmega"), straGen.pt(), nch, density);
          registry.fill(HIST("MC/densityMCGenOmegaMultTPC"), straGen.pt(), nMultTPC, density);
        }
        fGFW->Fill(straGen.eta(), fOmegaPtAxis->FindBin(straGen.pt()) - 1, straGen.phi(), 1, 256);
      }

      if (pdgCode == PDG_t::kK0Short) {
        int phibin = hLocalDensity->FindBin(RecoDecay::constrainAngle(straGen.phi(), -constants::math::PI));
        double density = hLocalDensity->Integral(phibin - cfgDeltaPhiLocDen, phibin + cfgDeltaPhiLocDen);
        if (cfgOutputLocDenWeights) {
          registry.fill(HIST("MC/densityMCGenK0s"), straGen.pt(), nch, density);
          registry.fill(HIST("MC/densityMCGenK0sMultTPC"), straGen.pt(), nMultTPC, density);
        }
        fGFW->Fill(straGen.eta(), fK0sPtAxis->FindBin(straGen.pt()) - 1, straGen.phi(), 1, 512);
      }
      if (pdgCode == PDG_t::kLambda0) {
        int phibin = hLocalDensity->FindBin(RecoDecay::constrainAngle(straGen.phi(), -constants::math::PI));
        double density = hLocalDensity->Integral(phibin - cfgDeltaPhiLocDen, phibin + cfgDeltaPhiLocDen);
        if (cfgOutputLocDenWeights) {
          registry.fill(HIST("MC/densityMCGenLambda"), straGen.pt(), nch, density);
          registry.fill(HIST("MC/densityMCGenLambdaMultTPC"), straGen.pt(), nMultTPC, density);
        }
        fGFW->Fill(straGen.eta(), fLambdaPtAxis->FindBin(straGen.pt()) - 1, straGen.phi(), 1, 1024);
      }
    }
    fillProfile(corrconfigs.at(40), HIST("MC/c22MC"), cent);
    for (int i = 1; i <= nK0sPtBins; i++) {
      fillProfilepTMC(corrconfigs.at(36), HIST("MC/K0sc22dptMC"), i, kK0Short, cent);
      fillProfilepTMC(corrconfigs.at(37), HIST("MC/K0sc22dptMC"), i, kK0Short, cent);
    }
    for (int i = 1; i <= nLambdaPtBins; i++) {
      fillProfilepTMC(corrconfigs.at(38), HIST("MC/Lambdac22dptMC"), i, kLambda0, cent);
      fillProfilepTMC(corrconfigs.at(39), HIST("MC/Lambdac22dptMC"), i, kLambda0, cent);
    }
    for (int i = 1; i <= nXiPtBins; i++) {
      fillProfilepTMC(corrconfigs.at(32), HIST("MC/Xic22dptMC"), i, kXiMinus, cent);
      fillProfilepTMC(corrconfigs.at(33), HIST("MC/Xic22dptMC"), i, kXiMinus, cent);
    }
    for (int i = 1; i <= nOmegaPtBins; i++) {
      fillProfilepTMC(corrconfigs.at(34), HIST("MC/Omegac22dptMC"), i, kOmegaMinus, cent);
      fillProfilepTMC(corrconfigs.at(35), HIST("MC/Omegac22dptMC"), i, kOmegaMinus, cent);
    }

    delete hLocalDensity;
  }
  PROCESS_SWITCH(FlowGfwOmegaXi, processMCGen, "", true);

  void processMCRec(AodCollisions::iterator const& collision, soa::Join<AodTracks, aod::McTrackLabels> const& tracks, aod::BCsWithTimestamps const&, soa::Join<aod::V0Datas, aod::McV0Labels> const& V0s, soa::Join<aod::CascDataExt, aod::McCascLabels> const& Cascades, DaughterTracks const&, aod::McParticles const&)
  {
    fGFW->Clear();
    const auto cent = collision.centFT0C();
    if (!collision.sel8())
      return;
    if (eventSelected(collision, cent))
      return;
    TH1D* hLocalDensity = new TH1D("hphi", "hphi", 400, -constants::math::TwoPI, constants::math::TwoPI);
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    loadCorrections(bc.timestamp());
    float vtxz = collision.posZ();
    registry.fill(HIST("hVtxZ"), vtxz);
    registry.fill(HIST("hCent"), cent);

    float weff = 1;
    float wacc = 1;
    float wloc = 1;
    double nch = 0;

    for (const auto& track : tracks) {
      if (track.pt() < trkQualityOpts.cfgCutPtPOIMin.value || track.pt() > trkQualityOpts.cfgCutPtPOIMax.value)
        continue;
      if (std::fabs(track.eta()) > trkQualityOpts.cfgCutEta.value)
        continue;
      if (!(track.isGlobalTrack()))
        continue;
      if (track.tpcChi2NCl() > cfgCutChi2prTPCcls)
        continue;
      if (cfgDoAccEffCorr) {
        if (!setCurrentParticleWeights(weff, wacc, track, vtxz, 0))
          continue;
      }
      nch += wacc * weff;
      if (!track.has_mcParticle())
        continue;
      auto mcParticle = track.mcParticle_as<aod::McParticles>();
      registry.fill(HIST("hPhi"), track.phi());
      registry.fill(HIST("hPhicorr"), track.phi(), wacc);
      registry.fill(HIST("hEta"), track.eta());
      registry.fill(HIST("hPt"), track.pt());
      int ptbin = fPtAxis->FindBin(track.pt()) - 1;
      if ((track.pt() > trkQualityOpts.cfgCutPtMin.value) && (track.pt() < trkQualityOpts.cfgCutPtMax.value)) {
        fGFW->Fill(track.eta(), ptbin, track.phi(), wacc * weff, 1); //(eta, ptbin, phi, wacc*weff, bitmask)
      }
      if ((track.pt() > trkQualityOpts.cfgCutPtPOIMin.value) && (track.pt() < trkQualityOpts.cfgCutPtPOIMax.value)) {
        fGFW->Fill(track.eta(), ptbin, track.phi(), wacc * weff, 32);
        if (cfgDoLocDenCorr) {
          hLocalDensity->Fill(mcParticle.phi(), wacc * weff);
          hLocalDensity->Fill(RecoDecay::constrainAngle(mcParticle.phi(), -constants::math::TwoPI), wacc * weff);
        }
      }
    }

    if (cfgDoLocDenCorr) {
      registry.fill(HIST("hCentvsNch"), cent, nch);
    }

    for (const auto& casc : Cascades) {
      if (!casc.has_mcParticle())
        continue;
      auto cascMC = casc.mcParticle_as<aod::McParticles>();
      auto negdau = casc.negTrack_as<DaughterTracks>();
      auto posdau = casc.posTrack_as<DaughterTracks>();
      auto bachelor = casc.bachelor_as<DaughterTracks>();
      int pdgCode{cascMC.pdgCode()};
      if (bachelor.pt() < trkQualityOpts.cfgCutPtDauMin.value || bachelor.pt() > trkQualityOpts.cfgCutPtDauMax.value)
        continue;
      if (posdau.pt() < trkQualityOpts.cfgCutPtDauMin.value || posdau.pt() > trkQualityOpts.cfgCutPtDauMax.value)
        continue;
      if (negdau.pt() < trkQualityOpts.cfgCutPtDauMin.value || negdau.pt() > trkQualityOpts.cfgCutPtDauMax.value)
        continue;
      // fill QA
      if (cfgOutputQA) {
        if (std::abs(pdgCode) == kXiMinus) {
          registry.fill(HIST("QAhisto/Xi/hqaCascRadiusbefore"), casc.dcabachtopv());
          registry.fill(HIST("QAhisto/Xi/hqaCasccosPAbefore"), casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
          registry.fill(HIST("QAhisto/Xi/hqaCascV0cosPAbefore"), casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
          registry.fill(HIST("QAhisto/Xi/hqadcaCascV0toPVbefore"), casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()));
          registry.fill(HIST("QAhisto/Xi/hqadcaCascBachtoPVbefore"), casc.dcabachtopv());
          registry.fill(HIST("QAhisto/Xi/hqadcaCascdaubefore"), casc.dcacascdaughters());
          registry.fill(HIST("QAhisto/Xi/hqadcaCascV0daubefore"), casc.dcaV0daughters());
        }
        if (std::abs(pdgCode) == kOmegaMinus) {
          registry.fill(HIST("QAhisto/Omega/hqaCascRadiusbefore"), casc.dcabachtopv());
          registry.fill(HIST("QAhisto/Omega/hqaCasccosPAbefore"), casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
          registry.fill(HIST("QAhisto/Omega/hqaCascV0cosPAbefore"), casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
          registry.fill(HIST("QAhisto/Omega/hqadcaCascV0toPVbefore"), casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()));
          registry.fill(HIST("QAhisto/Omega/hqadcaCascBachtoPVbefore"), casc.dcabachtopv());
          registry.fill(HIST("QAhisto/Omega/hqadcaCascdaubefore"), casc.dcacascdaughters());
          registry.fill(HIST("QAhisto/Omega/hqadcaCascV0daubefore"), casc.dcaV0daughters());
        }
      }
      // track quality check
      if (bachelor.itsNCls() <= trkQualityOpts.cfgITSNCls.value)
        continue;
      if (posdau.itsNCls() <= trkQualityOpts.cfgITSNCls.value)
        continue;
      if (negdau.itsNCls() <= trkQualityOpts.cfgITSNCls.value)
        continue;
      if (bachelor.tpcNClsFound() <= trkQualityOpts.cfgTPCNCls.value)
        continue;
      if (posdau.tpcNClsFound() <= trkQualityOpts.cfgTPCNCls.value)
        continue;
      if (negdau.tpcNClsFound() <= trkQualityOpts.cfgTPCNCls.value)
        continue;
      if (bachelor.tpcNClsCrossedRows() <= trkQualityOpts.cfgTPCCrossedRows.value)
        continue;
      if (posdau.tpcNClsCrossedRows() <= trkQualityOpts.cfgTPCCrossedRows.value)
        continue;
      if (negdau.tpcNClsCrossedRows() <= trkQualityOpts.cfgTPCCrossedRows.value)
        continue;
      if (trkQualityOpts.cfgCheckGlobalTrack.value) {
        if (!bachelor.hasTPC() || !bachelor.hasITS())
          continue;
        if (!posdau.hasTPC() || !posdau.hasITS())
          continue;
        if (!negdau.hasTPC() || !negdau.hasITS())
          continue;
      }
      // // topological cut
      if (casc.cascradius() < cascBuilderOpts.cfgcasc_radius.value)
        continue;
      if (casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) < cascBuilderOpts.cfgcasc_casccospa.value)
        continue;
      if (casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < cascBuilderOpts.cfgcasc_v0cospa.value)
        continue;
      if (std::fabs(casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ())) < cascBuilderOpts.cfgcasc_dcav0topv.value)
        continue;
      if (std::fabs(casc.dcabachtopv()) < cascBuilderOpts.cfgcasc_dcabachtopv.value)
        continue;
      if (casc.dcacascdaughters() > cascBuilderOpts.cfgcasc_dcacascdau.value)
        continue;
      if (casc.dcaV0daughters() > cascBuilderOpts.cfgcasc_dcav0dau.value)
        continue;
      if (std::fabs(casc.mLambda() - o2::constants::physics::MassLambda0) > cascBuilderOpts.cfgcasc_mlambdawindow.value)
        continue;
      // Omega and antiOmega
      double cascPt{cascMC.pt()};
      double cascPhi{cascMC.phi()};
      double cascEta{cascMC.eta()};
      if (std::abs(pdgCode) == kOmegaMinus) {
        if (casc.sign() < 0 && std::fabs(casc.yOmega()) < cfgCasc_rapidity &&
            (std::fabs(bachelor.tpcNSigmaKa()) < cfgNSigma[2] && std::fabs(posdau.tpcNSigmaPr()) < cfgNSigma[1] && std::fabs(negdau.tpcNSigmaPi()) < cfgNSigma[0])) {
          if (cfgDoAccEffCorr)
            setCurrentParticleWeights(weff, wacc, casc, vtxz, 4);
          if (cfgDoLocDenCorr) {
            int phibin = -999;
            phibin = hLocalDensity->FindBin(RecoDecay::constrainAngle(cascPhi, -constants::math::PI));
            double density = hLocalDensity->Integral(phibin - cfgDeltaPhiLocDen, phibin + cfgDeltaPhiLocDen);
            setCurrentLocalDensityWeights(wloc, casc, density, 4);
            if (cfgOutputLocDenWeights)
              registry.fill(HIST("MC/densityMCRecOmega"), cascPt, nch, density, casc.mOmega());
          }
          // fill QA
          if (cfgOutputQA) {
            registry.fill(HIST("QAhisto/Omega/hqaCascRadiusafter"), casc.cascradius());
            registry.fill(HIST("QAhisto/Omega/hqaCasccosPAafter"), casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
            registry.fill(HIST("QAhisto/Omega/hqaCascV0cosPAafter"), casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
            registry.fill(HIST("QAhisto/Omega/hqadcaCascV0toPVafter"), casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()));
            registry.fill(HIST("QAhisto/Omega/hqadcaCascBachtoPVafter"), casc.dcabachtopv());
            registry.fill(HIST("QAhisto/Omega/hqadcaCascdauafter"), casc.dcacascdaughters());
            registry.fill(HIST("QAhisto/Omega/hqadcaCascV0dauafter"), casc.dcaV0daughters());
          }
          fGFW->Fill(cascEta, fOmegaPtAxis->FindBin(cascPt) - 1, cascPhi, wacc * weff * wloc, 4);
        } else if (casc.sign() > 0 && std::fabs(casc.yOmega()) < cfgCasc_rapidity &&
                   (std::fabs(bachelor.tpcNSigmaKa()) < cfgNSigma[2] && std::fabs(negdau.tpcNSigmaPr()) < cfgNSigma[1] && std::fabs(posdau.tpcNSigmaPi()) < cfgNSigma[0])) {
          if (cfgDoAccEffCorr)
            setCurrentParticleWeights(weff, wacc, casc, vtxz, 4);
          if (cfgDoLocDenCorr) {
            int phibin = -999;
            phibin = hLocalDensity->FindBin(RecoDecay::constrainAngle(cascPhi, -constants::math::PI));
            double density = hLocalDensity->Integral(phibin - cfgDeltaPhiLocDen, phibin + cfgDeltaPhiLocDen);
            setCurrentLocalDensityWeights(wloc, casc, density, 4);
            if (cfgOutputLocDenWeights)
              registry.fill(HIST("MC/densityMCRecOmega"), cascPt, nch, density, casc.mOmega());
          }
          // fill QA
          if (cfgOutputQA) {
            registry.fill(HIST("QAhisto/Omega/hqaCascRadiusafter"), casc.cascradius());
            registry.fill(HIST("QAhisto/Omega/hqaCasccosPAafter"), casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
            registry.fill(HIST("QAhisto/Omega/hqaCascV0cosPAafter"), casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
            registry.fill(HIST("QAhisto/Omega/hqadcaCascV0toPVafter"), casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()));
            registry.fill(HIST("QAhisto/Omega/hqadcaCascBachtoPVafter"), casc.dcabachtopv());
            registry.fill(HIST("QAhisto/Omega/hqadcaCascdauafter"), casc.dcacascdaughters());
            registry.fill(HIST("QAhisto/Omega/hqadcaCascV0dauafter"), casc.dcaV0daughters());
          }
          fGFW->Fill(cascEta, fOmegaPtAxis->FindBin(cascPt) - 1, cascPhi, wacc * weff * wloc, 4);
        }
      }
      // Xi and antiXi
      if (std::abs(pdgCode) == kXiMinus) {
        if (casc.sign() < 0 && std::fabs(casc.yXi()) < cfgCasc_rapidity &&
            (std::fabs(bachelor.tpcNSigmaPi()) < cfgNSigma[0] && std::fabs(posdau.tpcNSigmaPr()) < cfgNSigma[1] && std::fabs(negdau.tpcNSigmaPi()) < cfgNSigma[0])) {
          if (cfgDoAccEffCorr)
            setCurrentParticleWeights(weff, wacc, casc, vtxz, 3);
          if (cfgDoLocDenCorr) {
            int phibin = -999;
            phibin = hLocalDensity->FindBin(RecoDecay::constrainAngle(cascPhi, -constants::math::PI));
            double density = hLocalDensity->Integral(phibin - cfgDeltaPhiLocDen, phibin + cfgDeltaPhiLocDen);
            setCurrentLocalDensityWeights(wloc, casc, density, 3);
            if (cfgOutputLocDenWeights)
              registry.fill(HIST("MC/densityMCRecXi"), cascPt, nch, density, casc.mXi());
          }
          // fill QA
          if (cfgOutputQA) {
            registry.fill(HIST("QAhisto/Xi/hqaCascRadiusafter"), casc.cascradius());
            registry.fill(HIST("QAhisto/Xi/hqaCasccosPAafter"), casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
            registry.fill(HIST("QAhisto/Xi/hqaCascV0cosPAafter"), casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
            registry.fill(HIST("QAhisto/Xi/hqadcaCascV0toPVafter"), casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()));
            registry.fill(HIST("QAhisto/Xi/hqadcaCascBachtoPVafter"), casc.dcabachtopv());
            registry.fill(HIST("QAhisto/Xi/hqadcaCascdauafter"), casc.dcacascdaughters());
            registry.fill(HIST("QAhisto/Xi/hqadcaCascV0dauafter"), casc.dcaV0daughters());
          }
          fGFW->Fill(cascEta, fXiPtAxis->FindBin(cascPt) - 1, cascPhi, wacc * weff * wloc, 2);
        } else if (casc.sign() > 0 && std::fabs(casc.yXi()) < cfgCasc_rapidity &&
                   (std::fabs(bachelor.tpcNSigmaPi()) < cfgNSigma[0] && std::fabs(negdau.tpcNSigmaPr()) < cfgNSigma[1] && std::fabs(posdau.tpcNSigmaPi()) < cfgNSigma[0])) {
          if (cfgDoAccEffCorr)
            setCurrentParticleWeights(weff, wacc, casc, vtxz, 3);
          if (cfgDoLocDenCorr) {
            int phibin = -999;
            phibin = hLocalDensity->FindBin(RecoDecay::constrainAngle(cascPhi, -constants::math::PI));
            double density = hLocalDensity->Integral(phibin - cfgDeltaPhiLocDen, phibin + cfgDeltaPhiLocDen);
            setCurrentLocalDensityWeights(wloc, casc, density, 3);
            if (cfgOutputLocDenWeights)
              registry.fill(HIST("MC/densityMCRecXi"), cascPt, nch, density, casc.mXi());
          }
          if (cfgOutputQA) {
            registry.fill(HIST("QAhisto/Xi/hqaCascRadiusafter"), casc.cascradius());
            registry.fill(HIST("QAhisto/Xi/hqaCasccosPAafter"), casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
            registry.fill(HIST("QAhisto/Xi/hqaCascV0cosPAafter"), casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
            registry.fill(HIST("QAhisto/Xi/hqadcaCascV0toPVafter"), casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()));
            registry.fill(HIST("QAhisto/Xi/hqadcaCascBachtoPVafter"), casc.dcabachtopv());
            registry.fill(HIST("QAhisto/Xi/hqadcaCascdauafter"), casc.dcacascdaughters());
            registry.fill(HIST("QAhisto/Xi/hqadcaCascV0dauafter"), casc.dcaV0daughters());
          }
          fGFW->Fill(cascEta, fXiPtAxis->FindBin(cascPt) - 1, cascPhi, wacc * weff * wloc, 2);
        }
      }
    }

    for (const auto& v0 : V0s) {
      if (!v0.has_mcParticle())
        continue;
      auto v0MC = v0.mcParticle_as<aod::McParticles>();
      auto v0negdau = v0.negTrack_as<DaughterTracks>();
      auto v0posdau = v0.posTrack_as<DaughterTracks>();

      if (v0posdau.pt() < trkQualityOpts.cfgCutPtDauMin.value || v0posdau.pt() > trkQualityOpts.cfgCutPtDauMax.value)
        continue;
      if (v0negdau.pt() < trkQualityOpts.cfgCutPtDauMin.value || v0negdau.pt() > trkQualityOpts.cfgCutPtDauMax.value)
        continue;

      // fill QA before cut
      if (cfgOutputQA) {
        registry.fill(HIST("QAhisto/V0/hqaV0radiusbefore"), v0.v0radius());
        registry.fill(HIST("QAhisto/V0/hqaV0cosPAbefore"), v0.v0cosPA());
        registry.fill(HIST("QAhisto/V0/hqadcaV0daubefore"), v0.dcaV0daughters());
        registry.fill(HIST("QAhisto/V0/hqadcapostoPVbefore"), v0.dcapostopv());
        registry.fill(HIST("QAhisto/V0/hqadcanegtoPVbefore"), v0.dcanegtopv());
        registry.fill(HIST("QAhisto/V0/hqaarm_podobefore"), v0.alpha(), v0.qtarm());
      }
      // // track quality check
      if (v0posdau.itsNCls() <= trkQualityOpts.cfgITSNCls.value)
        continue;
      if (v0negdau.itsNCls() <= trkQualityOpts.cfgITSNCls.value)
        continue;
      if (v0posdau.tpcNClsFound() <= trkQualityOpts.cfgTPCNCls.value)
        continue;
      if (v0negdau.tpcNClsFound() <= trkQualityOpts.cfgTPCNCls.value)
        continue;
      if (v0posdau.tpcNClsCrossedRows() <= trkQualityOpts.cfgTPCCrossedRows.value)
        continue;
      if (v0negdau.tpcNClsCrossedRows() <= trkQualityOpts.cfgTPCCrossedRows.value)
        continue;
      if (trkQualityOpts.cfgCheckGlobalTrack.value) {
        if (!v0posdau.hasTPC() || !v0posdau.hasITS())
          continue;
        if (!v0negdau.hasTPC() || !v0negdau.hasITS())
          continue;
      }
      // topological cut
      if (v0.v0radius() < v0BuilderOpts.cfgv0_radius.value)
        continue;
      if (v0.v0cosPA() < v0BuilderOpts.cfgv0_v0cospa.value)
        continue;
      if (v0.dcaV0daughters() > v0BuilderOpts.cfgv0_dcav0dau.value)
        continue;
      if (std::fabs(v0.dcapostopv()) < v0BuilderOpts.cfgv0_dcadautopv.value)
        continue;
      if (std::fabs(v0.dcanegtopv()) < v0BuilderOpts.cfgv0_dcadautopv.value)
        continue;
      // fill QA after cut
      if (cfgOutputQA) {
        registry.fill(HIST("QAhisto/V0/hqaV0radiusafter"), v0.v0radius());
        registry.fill(HIST("QAhisto/V0/hqaV0cosPAafter"), v0.v0cosPA());
        registry.fill(HIST("QAhisto/V0/hqadcaV0dauafter"), v0.dcaV0daughters());
        registry.fill(HIST("QAhisto/V0/hqadcapostoPVafter"), v0.dcapostopv());
        registry.fill(HIST("QAhisto/V0/hqadcanegtoPVafter"), v0.dcanegtopv());
      }

      int pdgCode{v0MC.pdgCode()};
      double v0Pt{v0MC.pt()};
      double v0Phi{v0MC.phi()};
      double v0Eta{v0MC.eta()};
      // K0short
      if (std::abs(pdgCode) == kK0Short) {
        if (v0.qtarm() / std::fabs(v0.alpha()) > v0BuilderOpts.cfgv0_ArmPodocut.value &&
            std::fabs(v0.mK0Short() - o2::constants::physics::MassK0Short) < v0BuilderOpts.cfgv0_mk0swindow.value &&
            (std::fabs(v0posdau.tpcNSigmaPi()) < cfgNSigma[0] && std::fabs(v0negdau.tpcNSigmaPi()) < cfgNSigma[0])) {
          if (cfgDoAccEffCorr)
            setCurrentParticleWeights(weff, wacc, v0, vtxz, 1);
          if (cfgDoLocDenCorr) {
            int phibin = -999;
            phibin = hLocalDensity->FindBin(RecoDecay::constrainAngle(v0Phi, -constants::math::PI));
            double density = hLocalDensity->Integral(phibin - cfgDeltaPhiLocDen, phibin + cfgDeltaPhiLocDen);
            setCurrentLocalDensityWeights(wloc, v0, density, 1);
            if (cfgOutputLocDenWeights)
              registry.fill(HIST("MC/densityMCRecK0s"), v0Pt, nch, density, v0.mK0Short());
          }
          fGFW->Fill(v0Eta, fK0sPtAxis->FindBin(v0Pt) - 1, v0Phi, wacc * weff * wloc, 8);
        }
      }
      // Lambda and antiLambda
      if (std::fabs(v0.mLambda() - o2::constants::physics::MassLambda) < v0BuilderOpts.cfgv0_mlambdawindow.value &&
          (std::fabs(v0posdau.tpcNSigmaPr()) < cfgNSigma[1] && std::fabs(v0negdau.tpcNSigmaPi()) < cfgNSigma[0])) {
        if (std::abs(pdgCode) == kLambda0) {
          if (cfgDoAccEffCorr)
            setCurrentParticleWeights(weff, wacc, v0, vtxz, 2);
          if (cfgDoLocDenCorr) {
            int phibin = -999;
            phibin = hLocalDensity->FindBin(RecoDecay::constrainAngle(v0Phi, -constants::math::PI));
            double density = hLocalDensity->Integral(phibin - cfgDeltaPhiLocDen, phibin + cfgDeltaPhiLocDen);
            setCurrentLocalDensityWeights(wloc, v0, density, 2);
            if (cfgOutputLocDenWeights)
              registry.fill(HIST("MC/densityMCRecLambda"), v0Pt, nch, density, v0.mLambda());
          }
          fGFW->Fill(v0Eta, fLambdaPtAxis->FindBin(v0Pt) - 1, v0Phi, wacc * weff * wloc, 16);
        }
      } else if (std::fabs(v0.mLambda() - o2::constants::physics::MassLambda) < v0BuilderOpts.cfgv0_mlambdawindow.value &&
                 (std::fabs(v0negdau.tpcNSigmaPr()) < cfgNSigma[1] && std::fabs(v0posdau.tpcNSigmaPi()) < cfgNSigma[0])) {
        if (std::abs(pdgCode) == kLambda0) {
          if (cfgDoAccEffCorr)
            setCurrentParticleWeights(weff, wacc, v0, vtxz, 2);
          if (cfgDoLocDenCorr) {
            int phibin = -999;
            phibin = hLocalDensity->FindBin(RecoDecay::constrainAngle(v0Phi, -constants::math::PI));
            double density = hLocalDensity->Integral(phibin - cfgDeltaPhiLocDen, phibin + cfgDeltaPhiLocDen);
            setCurrentLocalDensityWeights(wloc, v0, density, 2);
            if (cfgOutputLocDenWeights)
              registry.fill(HIST("MC/densityMCRecLambda"), v0Pt, nch, density, v0.mLambda());
          }
          fGFW->Fill(v0Eta, fLambdaPtAxis->FindBin(v0Pt) - 1, v0Phi, wacc * weff * wloc, 16);
        }
      }
    }
    delete hLocalDensity;
    fillProfile(corrconfigs.at(20), HIST("c22"), cent);
    fillProfile(corrconfigs.at(31), HIST("c32"), cent);
    for (int i = 1; i <= nK0sPtBins; i++) {
      fillProfilepTMass(corrconfigs.at(12), HIST("K0sc22dpt"), i, kK0Short, cent);
      fillProfilepTMass(corrconfigs.at(13), HIST("K0sc22dpt"), i, kK0Short, cent);
      fillProfilepTMass(corrconfigs.at(27), HIST("K0sc32dpt"), i, kK0Short, cent);
      fillProfilepTMass(corrconfigs.at(28), HIST("K0sc32dpt"), i, kK0Short, cent);
    }
    for (int i = 1; i <= nLambdaPtBins; i++) {
      fillProfilepTMass(corrconfigs.at(16), HIST("Lambdac22dpt"), i, kLambda0, cent);
      fillProfilepTMass(corrconfigs.at(17), HIST("Lambdac22dpt"), i, kLambda0, cent);
      fillProfilepTMass(corrconfigs.at(29), HIST("Lambdac32dpt"), i, kLambda0, cent);
      fillProfilepTMass(corrconfigs.at(30), HIST("Lambdac32dpt"), i, kLambda0, cent);
    }
    for (int i = 1; i <= nXiPtBins; i++) {
      fillProfilepTMass(corrconfigs.at(4), HIST("Xic22dpt"), i, kXiMinus, cent);
      fillProfilepTMass(corrconfigs.at(5), HIST("Xic22dpt"), i, kXiMinus, cent);
      fillProfilepTMass(corrconfigs.at(23), HIST("Xic32dpt"), i, kXiMinus, cent);
      fillProfilepTMass(corrconfigs.at(24), HIST("Xic32dpt"), i, kXiMinus, cent);
    }
    for (int i = 1; i <= nOmegaPtBins; i++) {
      fillProfilepTMass(corrconfigs.at(8), HIST("Omegac22dpt"), i, kOmegaMinus, cent);
      fillProfilepTMass(corrconfigs.at(9), HIST("Omegac22dpt"), i, kOmegaMinus, cent);
      fillProfilepTMass(corrconfigs.at(25), HIST("Omegac32dpt"), i, kOmegaMinus, cent);
      fillProfilepTMass(corrconfigs.at(26), HIST("Omegac32dpt"), i, kOmegaMinus, cent);
    }
  }
  PROCESS_SWITCH(FlowGfwOmegaXi, processMCRec, "", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowGfwOmegaXi>(cfgc)};
}
