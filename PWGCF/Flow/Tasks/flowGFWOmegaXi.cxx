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

/// \file   flowGFWOmegaXi.cxx
/// \author Fuchun Cui(fcui@cern.ch)
/// \since  Sep/13/2024
/// \brief  This task is to caculate V0s and cascades flow by GenericFramework

#include <CCDB/BasicCCDBManager.h>
#include <vector>
#include <string>
#include <cmath>
#include <memory>
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
namespace
{
std::shared_ptr<TProfile> refc22[10];
std::shared_ptr<TProfile> refc24[10];
std::shared_ptr<TProfile3D> k0sc22[10];
std::shared_ptr<TProfile3D> k0sc24[10];
std::shared_ptr<TProfile3D> lambdac22[10];
std::shared_ptr<TProfile3D> lambdac24[10];
std::shared_ptr<TProfile3D> xic22[10];
std::shared_ptr<TProfile3D> xic24[10];
std::shared_ptr<TProfile3D> omegac22[10];
std::shared_ptr<TProfile3D> omegac24[10];
} // namespace

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct FlowGFWOmegaXi {

  O2_DEFINE_CONFIGURABLE(cfgCutVertex, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMin, float, 0.2f, "Minimal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMax, float, 10.0f, "Maximal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMin, float, 0.2f, "Minimal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMax, float, 10.0f, "Maximal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutChi2prTPCcls, float, 2.5, "Chi2 per TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgCutOccupancyHigh, int, 500, "High cut on TPC occupancy")
  O2_DEFINE_CONFIGURABLE(cfgCutOccupancyLow, int, 0, "Low cut on TPC occupancy")
  O2_DEFINE_CONFIGURABLE(cfgOmegaMassbins, int, 16, "Number of Omega mass axis bins for c22")
  O2_DEFINE_CONFIGURABLE(cfgXiMassbins, int, 14, "Number of Xi mass axis bins for c22")
  O2_DEFINE_CONFIGURABLE(cfgK0sMassbins, int, 80, "Number of K0s mass axis bins for c22")
  O2_DEFINE_CONFIGURABLE(cfgLambdaMassbins, int, 32, "Number of Lambda mass axis bins for c22")
  // topological cut for V0
  O2_DEFINE_CONFIGURABLE(cfgv0_radius, float, 5.0f, "minimum decay radius")
  O2_DEFINE_CONFIGURABLE(cfgv0_v0cospa, float, 0.995f, "minimum cosine of pointing angle")
  O2_DEFINE_CONFIGURABLE(cfgv0_dcadautopv, float, 0.1f, "minimum daughter DCA to PV")
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
  O2_DEFINE_CONFIGURABLE(cfgtpcclufindable, int, 1, "minimum number of findable TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgtpccrossoverfindable, int, 1, "minimum number of Ratio crossed rows over findable clusters")
  O2_DEFINE_CONFIGURABLE(cfgcheckDauTPC, bool, true, "check daughter tracks TPC or not")
  O2_DEFINE_CONFIGURABLE(cfgcheckDauTOF, bool, false, "check daughter tracks TOF or not")
  O2_DEFINE_CONFIGURABLE(cfgCasc_rapidity, float, 0.5, "rapidity")
  O2_DEFINE_CONFIGURABLE(cfgtpcNSigmaCascPion, float, 3, "NSigmaCascPion")
  O2_DEFINE_CONFIGURABLE(cfgtpcNSigmaCascProton, float, 3, "NSigmaCascProton")
  O2_DEFINE_CONFIGURABLE(cfgtpcNSigmaCascKaon, float, 3, "NSigmaCascKaon")
  O2_DEFINE_CONFIGURABLE(cfgtofNSigmaCascPion, float, 3, "NSigmaCascPion")
  O2_DEFINE_CONFIGURABLE(cfgtofNSigmaCascProton, float, 3, "NSigmaCascProton")
  O2_DEFINE_CONFIGURABLE(cfgtofNSigmaCascKaon, float, 3, "NSigmaCascKaon")
  O2_DEFINE_CONFIGURABLE(cfgOutputNUAWeights, bool, true, "Fill and output NUA weights")
  O2_DEFINE_CONFIGURABLE(cfgAcceptancePath, std::vector<std::string>, (std::vector<std::string>{"Users/f/fcui/NUA/NUAREFPartical", "Users/f/fcui/NUA/NUAK0s", "Users/f/fcui/NUA/NUALambda", "Users/f/fcui/NUA/NUAXi", "Users/f/fcui/NUA/NUAOmega"}), "CCDB path to acceptance object")
  O2_DEFINE_CONFIGURABLE(cfgEfficiencyPath, std::vector<std::string>, (std::vector<std::string>{"PathtoRef"}), "CCDB path to efficiency object")

  ConfigurableAxis cfgaxisVertex{"cfgaxisVertex", {20, -10, 10}, "vertex axis for histograms"};
  ConfigurableAxis cfgaxisPhi{"cfgaxisPhi", {60, 0.0, constants::math::TwoPI}, "phi axis for histograms"};
  ConfigurableAxis cfgaxisEta{"cfgaxisEta", {40, -1., 1.}, "eta axis for histograms"};
  ConfigurableAxis cfgaxisPt{"cfgaxisPt", {VARIABLE_WIDTH, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.50, 4.00, 4.50, 5.00, 5.50, 6.00, 10.0}, "pt (GeV)"};
  ConfigurableAxis cfgaxisPtXi{"cfgaxisPtXi", {VARIABLE_WIDTH, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.9, 4.9, 5.9, 9.9}, "pt (GeV)"};
  ConfigurableAxis cfgaxisPtOmega{"cfgaxisPtOmega", {VARIABLE_WIDTH, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.9, 4.9, 5.9, 9.9}, "pt (GeV)"};
  ConfigurableAxis cfgaxisPtV0{"cfgaxisPtV0", {VARIABLE_WIDTH, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.9, 4.9, 5.9, 9.9}, "pt (GeV)"};
  ConfigurableAxis cfgaxisOmegaMassforflow{"cfgaxisOmegaMassforflow", {16, 1.63f, 1.71f}, "Inv. Mass (GeV)"};
  ConfigurableAxis cfgaxisXiMassforflow{"cfgaxisXiMassforflow", {14, 1.3f, 1.37f}, "Inv. Mass (GeV)"};
  ConfigurableAxis cfgaxisK0sMassforflow{"cfgaxisK0sMassforflow", {40, 0.4f, 0.6f}, "Inv. Mass (GeV)"};
  ConfigurableAxis cfgaxisLambdaMassforflow{"cfgaxisLambdaMassforflow", {32, 1.08f, 1.16f}, "Inv. Mass (GeV)"};

  AxisSpec axisMultiplicity{{0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90}, "Centrality (%)"};
  AxisSpec axisOmegaMass = {80, 1.63f, 1.71f, "Inv. Mass (GeV)"};
  AxisSpec axisXiMass = {70, 1.3f, 1.37f, "Inv. Mass (GeV)"};
  AxisSpec axisK0sMass = {400, 0.4f, 0.6f, "Inv. Mass (GeV)"};
  AxisSpec axisLambdaMass = {160, 1.08f, 1.16f, "Inv. Mass (GeV)"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtPOIMin) && (aod::track::pt < cfgCutPtPOIMax) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && (aod::track::tpcChi2NCl < cfgCutChi2prTPCcls);

  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;
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
  std::vector<std::string> cfgAcceptance = cfgAcceptancePath;
  std::vector<std::string> cfgEfficiency = cfgEfficiencyPath;

  std::vector<TH1D*> mEfficiency;
  std::vector<GFWWeights*> mAcceptance;
  bool correctionsLoaded = false;

  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;
  TF1* fMultCutLow = nullptr;
  TF1* fMultCutHigh = nullptr;
  TF1* fT0AV0AMean = nullptr;
  TF1* fT0AV0ASigma = nullptr;

  using TracksPID = soa::Join<aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr>;
  using AodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, TracksPID>>; // tracks filter
  using AodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::Mults>>;  // collisions filter
  using DaughterTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, TracksPID>;

  // Declare the pt, mult and phi Axis;
  int nPtBins = 0;
  TAxis* fPtAxis = nullptr;

  int nXiPtBins = 0;
  TAxis* fXiPtAxis = nullptr;

  int nV0PtBins = 0;
  TAxis* fV0PtAxis = nullptr;

  int nMultBins = 0;
  TAxis* fMultAxis = nullptr;

  int nPhiBins = 60;
  TAxis* fPhiAxis = new TAxis(nPhiBins, 0, constants::math::TwoPI);

  TAxis* fOmegaMass = nullptr;

  TAxis* fXiMass = nullptr;

  TAxis* fK0sMass = nullptr;

  TAxis* fLambdaMass = nullptr;

  void init(InitContext const&) // Initialization
  {
    ccdb->setURL(cfgurl.value);
    ccdb->setCaching(true);
    ccdb->setCreatedNotAfter(cfgnolaterthan.value);

    // Add some output objects to the histogram registry
    registry.add("hPhi", "", {HistType::kTH1D, {cfgaxisPhi}});
    registry.add("hPhicorr", "", {HistType::kTH1D, {cfgaxisPhi}});
    registry.add("hEta", "", {HistType::kTH1D, {cfgaxisEta}});
    registry.add("hVtxZ", "", {HistType::kTH1D, {cfgaxisVertex}});
    registry.add("hMult", "", {HistType::kTH1D, {{3000, 0.5, 3000.5}}});
    registry.add("hCent", "", {HistType::kTH1D, {{90, 0, 90}}});
    registry.add("hPt", "", {HistType::kTH1D, {cfgaxisPt}});
    registry.add("hEtaPhiVtxzREF", "", {HistType::kTH3D, {cfgaxisPhi, cfgaxisEta, {20, -10, 10}}});
    registry.add("hEtaPhiVtxzPOIXi", "", {HistType::kTH3D, {cfgaxisPhi, cfgaxisEta, {20, -10, 10}}});
    registry.add("hEtaPhiVtxzPOIOmega", "", {HistType::kTH3D, {cfgaxisPhi, cfgaxisEta, {20, -10, 10}}});
    registry.add("hEtaPhiVtxzPOIK0s", "", {HistType::kTH3D, {cfgaxisPhi, cfgaxisEta, {20, -10, 10}}});
    registry.add("hEtaPhiVtxzPOILambda", "", {HistType::kTH3D, {cfgaxisPhi, cfgaxisEta, {20, -10, 10}}});
    registry.add("hEventCount", "", {HistType::kTH2D, {{4, 0, 4}, {4, 0, 4}}});
    registry.get<TH2>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(1, "Filtered event");
    registry.get<TH2>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(2, "after sel8");
    registry.get<TH2>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(3, "before topological cut");
    registry.get<TH2>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(4, "after topological cut");
    registry.get<TH2>(HIST("hEventCount"))->GetYaxis()->SetBinLabel(1, "K0s");
    registry.get<TH2>(HIST("hEventCount"))->GetYaxis()->SetBinLabel(2, "Lambda");
    registry.get<TH2>(HIST("hEventCount"))->GetYaxis()->SetBinLabel(3, "XiMinus");
    registry.get<TH2>(HIST("hEventCount"))->GetYaxis()->SetBinLabel(4, "Omega");

    // QA
    registry.add("hqaV0radiusbefore", "", {HistType::kTH1D, {{200, 0, 200}}});
    registry.add("hqaV0radiusafter", "", {HistType::kTH1D, {{200, 0, 200}}});
    registry.add("hqaV0cosPAbefore", "", {HistType::kTH1D, {{1000, 0.95, 1}}});
    registry.add("hqaV0cosPAafter", "", {HistType::kTH1D, {{1000, 0.95, 1}}});
    registry.add("hqadcaV0daubefore", "", {HistType::kTH1D, {{100, 0, 1}}});
    registry.add("hqadcaV0dauafter", "", {HistType::kTH1D, {{100, 0, 1}}});
    registry.add("hqaarm_podobefore", "", {HistType::kTH2D, {{100, -1, 1}, {50, 0, 0.3}}});
    registry.add("hqaarm_podoafter", "", {HistType::kTH2D, {{100, -1, 1}, {50, 0, 0.3}}});
    registry.add("hqadcapostoPVbefore", "", {HistType::kTH1D, {{1000, -10, 10}}});
    registry.add("hqadcapostoPVafter", "", {HistType::kTH1D, {{1000, -10, 10}}});
    registry.add("hqadcanegtoPVbefore", "", {HistType::kTH1D, {{1000, -10, 10}}});
    registry.add("hqadcanegtoPVafter", "", {HistType::kTH1D, {{1000, -10, 10}}});

    // cumulant of flow
    registry.add("c22", ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("c24", ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("K0sc22", ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("Lambdac22", ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("c22dpt", ";Centrality  (%) ; C_{2}{2}", {HistType::kTProfile2D, {cfgaxisPt, axisMultiplicity}});
    registry.add("c24dpt", ";Centrality  (%) ; C_{2}{4}", {HistType::kTProfile2D, {cfgaxisPt, axisMultiplicity}});
    // pt-diff cumulant of flow
    registry.add("Xic22dpt", ";pt ; C_{2}{2} ", {HistType::kTProfile3D, {cfgaxisPtXi, cfgaxisXiMassforflow, axisMultiplicity}});
    registry.add("Omegac22dpt", ";pt ; C_{2}{2} ", {HistType::kTProfile3D, {cfgaxisPtXi, cfgaxisOmegaMassforflow, axisMultiplicity}});
    registry.add("K0sc22dpt", ";pt ; C_{2}{2} ", {HistType::kTProfile3D, {cfgaxisPtV0, cfgaxisK0sMassforflow, axisMultiplicity}});
    registry.add("Lambdac22dpt", ";pt ; C_{2}{2} ", {HistType::kTProfile3D, {cfgaxisPtV0, cfgaxisLambdaMassforflow, axisMultiplicity}});
    registry.add("Xic24dpt", ";pt ; C_{2}{4} ", {HistType::kTProfile3D, {cfgaxisPtXi, cfgaxisXiMassforflow, axisMultiplicity}});
    registry.add("Omegac24dpt", ";pt ; C_{2}{4} ", {HistType::kTProfile3D, {cfgaxisPtXi, cfgaxisOmegaMassforflow, axisMultiplicity}});
    registry.add("K0sc24dpt", ";pt ; C_{2}{4} ", {HistType::kTProfile3D, {cfgaxisPtV0, cfgaxisK0sMassforflow, axisMultiplicity}});
    registry.add("Lambdac24dpt", ";pt ; C_{2}{4} ", {HistType::kTProfile3D, {cfgaxisPtV0, cfgaxisLambdaMassforflow, axisMultiplicity}});
    // for Jackknife
    for (int i = 1; i <= 10; i++) {
      refc22[i - 1] = registry.add<TProfile>(Form("Jackknife/REF/c22_%d", i), ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile, {axisMultiplicity}});
      refc24[i - 1] = registry.add<TProfile>(Form("Jackknife/REF/c24_%d", i), ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile, {axisMultiplicity}});
      xic22[i - 1] = registry.add<TProfile3D>(Form("Jackknife/Xi/Xic22dpt_%d", i), ";pt ; C_{2}{2} ", {HistType::kTProfile3D, {cfgaxisPtXi, cfgaxisXiMassforflow, axisMultiplicity}});
      omegac22[i - 1] = registry.add<TProfile3D>(Form("Jackknife/Omega/Omegac22dpt_%d", i), ";pt ; C_{2}{2} ", {HistType::kTProfile3D, {cfgaxisPtXi, cfgaxisOmegaMassforflow, axisMultiplicity}});
      k0sc22[i - 1] = registry.add<TProfile3D>(Form("Jackknife/K0s/K0sc22dpt_%d", i), ";pt ; C_{2}{2} ", {HistType::kTProfile3D, {cfgaxisPtV0, cfgaxisK0sMassforflow, axisMultiplicity}});
      lambdac22[i - 1] = registry.add<TProfile3D>(Form("Jackknife/Lambda/Lambdac22dpt_%d", i), ";pt ; C_{2}{2} ", {HistType::kTProfile3D, {cfgaxisPtV0, cfgaxisLambdaMassforflow, axisMultiplicity}});
      xic24[i - 1] = registry.add<TProfile3D>(Form("Jackknife/Xi/Xic24dpt_%d", i), ";pt ; C_{2}{4} ", {HistType::kTProfile3D, {cfgaxisPtXi, cfgaxisXiMassforflow, axisMultiplicity}});
      omegac24[i - 1] = registry.add<TProfile3D>(Form("Jackknife/Omega/Omegac24dpt_%d", i), ";pt ; C_{2}{4} ", {HistType::kTProfile3D, {cfgaxisPtXi, cfgaxisOmegaMassforflow, axisMultiplicity}});
      k0sc24[i - 1] = registry.add<TProfile3D>(Form("Jackknife/K0s/K0sc24dpt_%d", i), ";pt ; C_{2}{4} ", {HistType::kTProfile3D, {cfgaxisPtV0, cfgaxisK0sMassforflow, axisMultiplicity}});
      lambdac24[i - 1] = registry.add<TProfile3D>(Form("Jackknife/Lambda/Lambdac24dpt_%d", i), ";pt ; C_{2}{4} ", {HistType::kTProfile3D, {cfgaxisPtV0, cfgaxisLambdaMassforflow, axisMultiplicity}});
    }
    // InvMass(GeV) of casc and v0
    registry.add("InvMassXi_all", "", {HistType::kTHnSparseF, {cfgaxisPtXi, axisXiMass, cfgaxisEta, axisMultiplicity}});
    registry.add("InvMassOmega_all", "", {HistType::kTHnSparseF, {cfgaxisPtXi, axisOmegaMass, cfgaxisEta, axisMultiplicity}});
    registry.add("InvMassOmega", "", {HistType::kTHnSparseF, {cfgaxisPtXi, axisOmegaMass, cfgaxisEta, axisMultiplicity}});
    registry.add("InvMassXi", "", {HistType::kTHnSparseF, {cfgaxisPtXi, axisXiMass, cfgaxisEta, axisMultiplicity}});
    registry.add("InvMassK0s_all", "", {HistType::kTHnSparseF, {cfgaxisPtV0, axisK0sMass, cfgaxisEta, axisMultiplicity}});
    registry.add("InvMassLambda_all", "", {HistType::kTHnSparseF, {cfgaxisPtV0, axisLambdaMass, cfgaxisEta, axisMultiplicity}});
    registry.add("InvMassK0s", "", {HistType::kTHnSparseF, {cfgaxisPtV0, axisK0sMass, cfgaxisEta, axisMultiplicity}});
    registry.add("InvMassLambda", "", {HistType::kTHnSparseF, {cfgaxisPtV0, axisLambdaMass, cfgaxisEta, axisMultiplicity}});

    // Set the pt, mult and phi Axis;
    o2::framework::AxisSpec axisPt = cfgaxisPt;
    nPtBins = axisPt.binEdges.size() - 1;
    fPtAxis = new TAxis(nPtBins, &(axisPt.binEdges)[0]);

    o2::framework::AxisSpec axisXiPt = cfgaxisPtXi;
    nXiPtBins = axisXiPt.binEdges.size() - 1;
    fXiPtAxis = new TAxis(nXiPtBins, &(axisXiPt.binEdges)[0]);

    o2::framework::AxisSpec axisV0Pt = cfgaxisPtV0;
    nV0PtBins = axisV0Pt.binEdges.size() - 1;
    fV0PtAxis = new TAxis(nV0PtBins, &(axisV0Pt.binEdges)[0]);

    o2::framework::AxisSpec axisMult = axisMultiplicity;
    nMultBins = axisMult.binEdges.size() - 1;
    fMultAxis = new TAxis(nMultBins, &(axisMult.binEdges)[0]);

    fOmegaMass = new TAxis(cfgOmegaMassbins, 1.63, 1.71);

    fXiMass = new TAxis(cfgXiMassbins, 1.3, 1.37);

    fK0sMass = new TAxis(cfgK0sMassbins, 0.4, 0.6);

    fLambdaMass = new TAxis(cfgLambdaMassbins, 1.08, 1.16);

    fGFW->AddRegion("reffull", -0.8, 0.8, 1, 1); // ("name", etamin, etamax, ptbinnum, bitmask)eta region -0.8 to 0.8
    fGFW->AddRegion("refN10", -0.8, -0.4, 1, 1);
    fGFW->AddRegion("refP10", 0.4, 0.8, 1, 1);
    // POI
    fGFW->AddRegion("poiN10dpt", -0.8, -0.4, nPtBins, 32);
    fGFW->AddRegion("poiP10dpt", 0.4, 0.8, nPtBins, 32);
    fGFW->AddRegion("poifulldpt", -0.8, 0.8, nPtBins, 32);
    fGFW->AddRegion("poioldpt", -0.8, 0.8, nPtBins, 1);

    int nXiptMassBins = nXiPtBins * cfgXiMassbins;
    fGFW->AddRegion("poiXiPdpt", 0.4, 0.8, nXiptMassBins, 2);
    fGFW->AddRegion("poiXiNdpt", -0.8, -0.4, nXiptMassBins, 2);
    fGFW->AddRegion("poiXifulldpt", -0.8, 0.8, nXiptMassBins, 2);
    fGFW->AddRegion("poiXiP", 0.4, 0.8, 1, 2);
    fGFW->AddRegion("poiXiN", -0.8, -0.4, 1, 2);
    int nOmegaptMassBins = nXiPtBins * cfgOmegaMassbins;
    fGFW->AddRegion("poiOmegaPdpt", 0.4, 0.8, nOmegaptMassBins, 4);
    fGFW->AddRegion("poiOmegaNdpt", -0.8, -0.4, nOmegaptMassBins, 4);
    fGFW->AddRegion("poiOmegafulldpt", -0.8, 0.8, nOmegaptMassBins, 4);
    fGFW->AddRegion("poiOmegaP", 0.4, 0.8, 1, 4);
    fGFW->AddRegion("poiOmegaN", -0.8, -0.4, 1, 4);
    int nK0sptMassBins = nV0PtBins * cfgK0sMassbins;
    fGFW->AddRegion("poiK0sPdpt", 0.4, 0.8, nK0sptMassBins, 8);
    fGFW->AddRegion("poiK0sNdpt", -0.8, -0.4, nK0sptMassBins, 8);
    fGFW->AddRegion("poiK0sfulldpt", -0.8, 0.8, nK0sptMassBins, 8);
    fGFW->AddRegion("poiK0sP", 0.4, 0.8, 1, 8);
    fGFW->AddRegion("poiK0sN", -0.8, 0.4, 1, 8);
    int nLambdaptMassBins = nV0PtBins * cfgLambdaMassbins;
    fGFW->AddRegion("poiLambdaPdpt", 0.4, 0.8, nLambdaptMassBins, 16);
    fGFW->AddRegion("poiLambdaNdpt", -0.8, -0.4, nLambdaptMassBins, 16);
    fGFW->AddRegion("poiLambdafulldpt", -0.8, 0.8, nLambdaptMassBins, 16);
    fGFW->AddRegion("poiLambdaP", 0.4, 0.8, 1, 16);
    fGFW->AddRegion("poiLambdaN", -0.8, -0.4, 1, 16);
    // pushback
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiP10dpt {2} refN10 {-2}", "Poi10Gap22dpta", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiN10dpt {2} refP10 {-2}", "Poi10Gap22dptb", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poifulldpt reffull | poioldpt {2 2 -2 -2}", "Poi10Gap24dpt", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiXiPdpt {2} refN10 {-2}", "Xi10Gap22a", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiXiNdpt {2} refP10 {-2}", "Xi10Gap22b", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiXifulldpt reffull {2 2 -2 -2}", "Xi10Gap24", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiOmegaPdpt {2} refN10 {-2}", "Omega10Gap22a", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiOmegaNdpt {2} refP10 {-2}", "Omega10Gap22b", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiOmegafulldpt reffull {2 2 -2 -2}", "Xi10Gap24", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiK0sPdpt {2} refN10 {-2}", "K0short10Gap22a", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiK0sNdpt {2} refP10 {-2}", "K0short10Gap22b", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiK0sfulldpt reffull {2 2 -2 -2}", "Xi10Gap24", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiLambdaPdpt {2} refN10 {-2}", "Lambda10Gap22a", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiLambdaNdpt {2} refP10 {-2}", "Lambda10Gap22b", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiLambdafulldpt reffull {2 2 -2 -2}", "Xi10Gap24a", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refP10 {2} refN10 {-2}", "Ref10Gap22a", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("reffull reffull {2 2 -2 -2}", "Ref10Gap24", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiK0sP {2} refN10 {-2}", "K0s10Gap22inta", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiK0sN {2} refP10 {-2}", "K0s10Gap22intb", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiLambdaP {2} refN10 {-2}", "Lambda10Gap22inta", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiLambdaN {2} refP10 {-2}", "Lambda10Gap22intb", kFALSE));
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

  // input HIST("name")
  template <char... chars>
  void fillProfilepTMass(const GFW::CorrConfig& corrconf, const ConstStr<chars...>& tarName, const int& ptbin, const int& PDGCode, const float& cent)
  {
    int nMassBins = 0;
    int nptbins = 0;
    TAxis* fMass = nullptr;
    TAxis* fpt = nullptr;
    if (PDGCode == kXiMinus) {
      nMassBins = cfgXiMassbins;
      nptbins = nXiPtBins;
      fpt = fXiPtAxis;
      fMass = fXiMass;
    } else if (PDGCode == kOmegaMinus) {
      nMassBins = cfgOmegaMassbins;
      nptbins = nXiPtBins;
      fpt = fXiPtAxis;
      fMass = fOmegaMass;
    } else if (PDGCode == kK0Short) {
      nMassBins = cfgK0sMassbins;
      nptbins = nV0PtBins;
      fpt = fV0PtAxis;
      fMass = fK0sMass;
    } else if (PDGCode == kLambda0) {
      nMassBins = cfgLambdaMassbins;
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
      nMassBins = cfgXiMassbins;
      nptbins = nXiPtBins;
      fpt = fXiPtAxis;
      fMass = fXiMass;
    } else if (PDGCode == kOmegaMinus) {
      nMassBins = cfgOmegaMassbins;
      nptbins = nXiPtBins;
      fpt = fXiPtAxis;
      fMass = fOmegaMass;
    } else if (PDGCode == kK0Short) {
      nMassBins = cfgK0sMassbins;
      nptbins = nV0PtBins;
      fpt = fV0PtAxis;
      fMass = fK0sMass;
    } else if (PDGCode == kLambda0) {
      nMassBins = cfgLambdaMassbins;
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
    if (cfgAcceptance.size() == 5) {
      for (int i = 0; i <= 4; i++) {
        mAcceptance.push_back(ccdb->getForTimeStamp<GFWWeights>(cfgAcceptance[i], timestamp));
      }
      if (mAcceptance.size() == 5)
        LOGF(info, "Loaded acceptance weights");
      else
        LOGF(warning, "Could not load acceptance weights");
    }
    if (cfgEfficiency.size() == 5) {
      for (int i = 0; i <= 4; i++) {
        mAcceptance.push_back(ccdb->getForTimeStamp<GFWWeights>(cfgAcceptance[i], timestamp));
      }
      if (mEfficiency.size() == 5)
        LOGF(info, "Loaded efficiency histogram");
      else
        LOGF(fatal, "Could not load efficiency histogram");
    }
    correctionsLoaded = true;
  }

  template <typename TrackObject>
  bool setCurrentParticleWeights(float& weight_nue, float& weight_nua, TrackObject track, float vtxz, int ispecies)
  {
    float eff = 1.;
    if (mEfficiency.size() == 5)
      eff = mEfficiency[ispecies]->GetBinContent(mEfficiency[ispecies]->FindBin(track.pt()));
    else
      eff = 1.0;
    if (eff == 0)
      return false;
    weight_nue = 1. / eff;
    if (mAcceptance.size() == 5)
      weight_nua = mAcceptance[ispecies]->getNUA(track.phi(), track.eta(), vtxz);
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
    if (!collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      // no collisions in specified time range
      return 0;
    }
    if (!collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      // cut time intervals with dead ITS staves
      return 0;
    }
    float vtxz = -999;
    if (collision.numContrib() > 1) {
      vtxz = collision.posZ();
      float zRes = std::sqrt(collision.covZZ());
      if (zRes > 0.25 && collision.numContrib() < 20)
        vtxz = -999;
    }
    auto multNTracksPV = collision.multNTracksPV();
    auto occupancy = collision.trackOccupancyInTimeRange();

    if (std::fabs(vtxz) > cfgCutVertex)
      return false;
    if (multNTracksPV < fMultPVCutLow->Eval(centrality))
      return false;
    if (multNTracksPV > fMultPVCutHigh->Eval(centrality))
      return false;
    if (multTrk < fMultCutLow->Eval(centrality))
      return false;
    if (multTrk > fMultCutHigh->Eval(centrality))
      return false;
    if (occupancy < cfgCutOccupancyLow || occupancy > cfgCutOccupancyHigh)
      return 0;

    // V0A T0A 5 sigma cut
    if (std::fabs(collision.multFV0A() - fT0AV0AMean->Eval(collision.multFT0A())) > 5 * fT0AV0ASigma->Eval(collision.multFT0A()))
      return 0;

    return true;
  }

  void process(AodCollisions::iterator const& collision, aod::BCsWithTimestamps const&, AodTracks const& tracks, aod::CascDataExt const& Cascades, aod::V0Datas const& V0s, DaughterTracks const&)
  {
    int nTot = tracks.size();
    int candNumAll[4] = {0, 0, 0, 0};
    int candNum[4] = {0, 0, 0, 0};
    for (int i = 0; i < 4; i++) {
      registry.fill(HIST("hEventCount"), 0.5, i + 0.5);
    }
    if (nTot < 1)
      return;
    fGFW->Clear();
    const auto cent = collision.centFT0C();
    if (!collision.sel8())
      return;
    if (eventSelected(collision, tracks.size(), cent))
      return;
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    loadCorrections(bc.timestamp());
    float vtxz = collision.posZ();
    registry.fill(HIST("hVtxZ"), vtxz);
    registry.fill(HIST("hMult"), nTot);
    registry.fill(HIST("hCent"), collision.centFT0C());
    for (int i = 0; i < 4; i++) {
      registry.fill(HIST("hEventCount"), 1.5, i + 0.5);
    }

    float weff = 1;
    float wacc = 1;
    // fill GFW ref flow
    for (const auto& track : tracks) {
      if (!setCurrentParticleWeights(weff, wacc, track, vtxz, 0))
        continue;
      registry.fill(HIST("hPhi"), track.phi());
      registry.fill(HIST("hPhicorr"), track.phi(), wacc);
      registry.fill(HIST("hEta"), track.eta());
      registry.fill(HIST("hEtaPhiVtxzREF"), track.phi(), track.eta(), vtxz, wacc);
      registry.fill(HIST("hPt"), track.pt());
      int ptbin = fPtAxis->FindBin(track.pt()) - 1;
      if ((track.pt() > cfgCutPtMin) && (track.pt() < cfgCutPtMax)) {
        fGFW->Fill(track.eta(), ptbin, track.phi(), wacc * weff, 1); //(eta, ptbin, phi, wacc*weff, bitmask)
      }
      if ((track.pt() > cfgCutPtPOIMin) && (track.pt() < cfgCutPtPOIMax)) {
        fGFW->Fill(track.eta(), ptbin, track.phi(), wacc * weff, 32);
      }
      if (cfgOutputNUAWeights)
        fWeightsREF->fill(track.phi(), track.eta(), vtxz, track.pt(), cent, 0);
    }
    // fill GFW of V0 flow
    for (const auto& v0 : V0s) {
      auto v0posdau = v0.posTrack_as<DaughterTracks>();
      auto v0negdau = v0.negTrack_as<DaughterTracks>();
      // check tpc
      bool isK0s = false;
      bool isLambda = false;
      // fill QA
      registry.fill(HIST("hqaarm_podobefore"), v0.alpha(), v0.qtarm());
      // check daughter TPC and TOF
      // K0short
      if (v0.qtarm() / std::fabs(v0.alpha()) > cfgv0_ArmPodocut && std::fabs(v0.y()) < 0.5 && std::fabs(v0.mK0Short() - o2::constants::physics::MassK0Short) < cfgv0_mk0swindow &&
          (!cfgcheckDauTPC || (std::fabs(v0posdau.tpcNSigmaPi()) < cfgtpcNSigmaCascPion && std::fabs(v0negdau.tpcNSigmaPi()) < cfgtpcNSigmaCascPion)) &&
          (!cfgcheckDauTOF || ((std::fabs(v0posdau.tofNSigmaPi()) < cfgtofNSigmaCascPion || v0posdau.pt() < 0.4) && (std::fabs(v0negdau.tofNSigmaPi()) < cfgtofNSigmaCascPion || v0negdau.pt() < 0.4)))) {
        registry.fill(HIST("InvMassK0s_all"), v0.pt(), v0.mK0Short(), v0.eta(), cent);
        if (!setCurrentParticleWeights(weff, wacc, v0, vtxz, 1))
          continue;
        isK0s = true;
        candNumAll[0] = candNumAll[0] + 1;
        registry.fill(HIST("hqaarm_podoafter"), v0.alpha(), v0.qtarm());
      }
      // Lambda and antiLambda
      if (std::fabs(v0.y()) < 0.5 && std::fabs(v0.mLambda() - o2::constants::physics::MassLambda) < cfgv0_mlambdawindow &&
          (!cfgcheckDauTPC || (std::fabs(v0posdau.tpcNSigmaPr()) < cfgtpcNSigmaCascProton && std::fabs(v0negdau.tpcNSigmaPi()) < cfgtpcNSigmaCascPion)) &&
          (!cfgcheckDauTOF || ((std::fabs(v0posdau.tofNSigmaPr()) < cfgtofNSigmaCascProton || v0posdau.pt() < 0.4) && (std::fabs(v0negdau.tofNSigmaPi()) < cfgtofNSigmaCascPion || v0negdau.pt() < 0.4)))) {
        registry.fill(HIST("InvMassLambda_all"), v0.pt(), v0.mLambda(), v0.eta(), cent);
        if (!setCurrentParticleWeights(weff, wacc, v0, vtxz, 2))
          continue;
        isLambda = true;
        candNumAll[1] = candNumAll[1] + 1;
      } else if (std::fabs(v0.y()) < 0.5 && std::fabs(v0.mLambda() - o2::constants::physics::MassLambda) < cfgv0_mlambdawindow &&
                 (!cfgcheckDauTPC || (std::fabs(v0negdau.tpcNSigmaPr()) < cfgtpcNSigmaCascProton && std::fabs(v0posdau.tpcNSigmaPi()) < cfgtpcNSigmaCascPion)) &&
                 (!cfgcheckDauTOF || ((std::fabs(v0negdau.tofNSigmaPr()) < cfgtofNSigmaCascProton || v0negdau.pt() < 0.4) && (std::fabs(v0posdau.tofNSigmaPi()) < cfgtofNSigmaCascPion || v0posdau.pt() < 0.4)))) {
        registry.fill(HIST("InvMassLambda_all"), v0.pt(), v0.mLambda(), v0.eta(), cent);
        if (!setCurrentParticleWeights(weff, wacc, v0, vtxz, 2))
          continue;
        isLambda = true;
        candNumAll[1] = candNumAll[1] + 1;
      }
      // fill QA before cut
      registry.fill(HIST("hqaV0radiusbefore"), v0.v0radius());
      registry.fill(HIST("hqaV0cosPAbefore"), v0.v0cosPA());
      registry.fill(HIST("hqadcaV0daubefore"), v0.dcaV0daughters());
      registry.fill(HIST("hqadcapostoPVbefore"), v0.dcapostopv());
      registry.fill(HIST("hqadcanegtoPVbefore"), v0.dcanegtopv());

      // track quality check
      if (v0posdau.tpcNClsFound() < cfgtpcclusters)
        continue;
      if (v0negdau.tpcNClsFound() < cfgtpcclusters)
        continue;
      if (v0posdau.tpcNClsFindable() < cfgtpcclufindable)
        continue;
      if (v0negdau.tpcNClsFindable() < cfgtpcclufindable)
        continue;
      if (v0posdau.tpcCrossedRowsOverFindableCls() < cfgtpccrossoverfindable)
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
      if (std::fabs(v0.dcapostopv()) < cfgv0_dcadautopv)
        continue;
      if (std::fabs(v0.dcanegtopv()) < cfgv0_dcadautopv)
        continue;
      // fill QA after cut
      registry.fill(HIST("hqaV0radiusafter"), v0.v0radius());
      registry.fill(HIST("hqaV0cosPAafter"), v0.v0cosPA());
      registry.fill(HIST("hqadcaV0dauafter"), v0.dcaV0daughters());
      registry.fill(HIST("hqadcapostoPVafter"), v0.dcapostopv());
      registry.fill(HIST("hqadcanegtoPVafter"), v0.dcanegtopv());
      if (isK0s) {
        candNum[0] = candNum[0] + 1;
        registry.fill(HIST("InvMassK0s"), v0.pt(), v0.mK0Short(), v0.eta(), cent);
        registry.fill(HIST("hEtaPhiVtxzPOIK0s"), v0.phi(), v0.eta(), vtxz, wacc);
        fGFW->Fill(v0.eta(), fV0PtAxis->FindBin(v0.pt()) - 1 + ((fK0sMass->FindBin(v0.mK0Short()) - 1) * nV0PtBins), v0.phi(), wacc * weff, 8);
        if (cfgOutputNUAWeights)
          fWeightsK0s->fill(v0.phi(), v0.eta(), vtxz, v0.pt(), cent, 0);
      }
      if (isLambda) {
        candNum[1] = candNum[1] + 1;
        registry.fill(HIST("InvMassLambda"), v0.pt(), v0.mLambda(), v0.eta(), cent);
        registry.fill(HIST("hEtaPhiVtxzPOILambda"), v0.phi(), v0.eta(), vtxz, wacc);
        fGFW->Fill(v0.eta(), fV0PtAxis->FindBin(v0.pt()) - 1 + ((fLambdaMass->FindBin(v0.mLambda()) - 1) * nV0PtBins), v0.phi(), wacc * weff, 16);
        if (cfgOutputNUAWeights)
          fWeightsLambda->fill(v0.phi(), v0.eta(), vtxz, v0.pt(), cent, 0);
      }
    }
    // fill GFW of casc flow
    for (const auto& casc : Cascades) {
      auto bachelor = casc.bachelor_as<DaughterTracks>();
      auto posdau = casc.posTrack_as<DaughterTracks>();
      auto negdau = casc.negTrack_as<DaughterTracks>();
      // check TPC
      if (cfgcheckDauTPC && (!posdau.hasTPC() || !negdau.hasTPC() || !bachelor.hasTPC())) {
        continue;
      }
      bool isOmega = false;
      bool isXi = false;
      // Omega and antiOmega
      if (casc.sign() < 0 && (casc.mOmega() > 1.63) && (casc.mOmega() < 1.71) && std::fabs(casc.yOmega()) < cfgCasc_rapidity &&
          (!cfgcheckDauTPC || (std::fabs(bachelor.tpcNSigmaKa()) < cfgtpcNSigmaCascKaon && std::fabs(posdau.tpcNSigmaPr()) < cfgtpcNSigmaCascProton && std::fabs(negdau.tpcNSigmaPi()) < cfgtpcNSigmaCascPion)) &&
          (!cfgcheckDauTOF || ((std::fabs(bachelor.tofNSigmaKa()) < cfgtofNSigmaCascKaon || bachelor.pt() < 0.4) && (std::fabs(posdau.tofNSigmaPr()) < cfgtofNSigmaCascProton || posdau.pt() < 0.4) && (std::fabs(negdau.tofNSigmaPi()) < cfgtofNSigmaCascPion || negdau.pt() < 0.4)))) {
        registry.fill(HIST("InvMassOmega_all"), casc.pt(), casc.mOmega(), casc.eta(), cent);
        if (!setCurrentParticleWeights(weff, wacc, casc, vtxz, 4))
          continue;
        isOmega = true;
        candNumAll[3] = candNumAll[3] + 1;
      } else if (casc.sign() < 0 && (casc.mOmega() > 1.63) && (casc.mOmega() < 1.71) && std::fabs(casc.yOmega()) < cfgCasc_rapidity &&
                 (!cfgcheckDauTPC || (std::fabs(bachelor.tpcNSigmaKa()) < cfgtpcNSigmaCascKaon && std::fabs(negdau.tpcNSigmaPr()) < cfgtpcNSigmaCascProton && std::fabs(posdau.tpcNSigmaPi()) < cfgtpcNSigmaCascPion)) &&
                 (!cfgcheckDauTOF || ((std::fabs(bachelor.tofNSigmaKa()) < cfgtofNSigmaCascKaon || bachelor.pt() < 0.4) && (std::fabs(negdau.tofNSigmaPr()) < cfgtofNSigmaCascProton || negdau.pt() < 0.4) && (std::fabs(posdau.tofNSigmaPi()) < cfgtofNSigmaCascPion || posdau.pt() < 0.4)))) {
        registry.fill(HIST("InvMassOmega_all"), casc.pt(), casc.mOmega(), casc.eta(), cent);
        if (!setCurrentParticleWeights(weff, wacc, casc, vtxz, 4))
          continue;
        isOmega = true;
        candNumAll[3] = candNumAll[3] + 1;
      }
      // Xi and antiXi
      if (casc.sign() < 0 && (casc.mXi() > 1.30) && (casc.mXi() < 1.37) && std::fabs(casc.yXi()) < cfgCasc_rapidity &&
          (!cfgcheckDauTPC || (std::fabs(bachelor.tpcNSigmaPi()) < cfgtpcNSigmaCascPion && std::fabs(posdau.tpcNSigmaPr()) < cfgtpcNSigmaCascProton && std::fabs(negdau.tpcNSigmaPi()) < cfgtpcNSigmaCascPion)) &&
          (!cfgcheckDauTOF || ((std::fabs(bachelor.tofNSigmaPi()) < cfgtofNSigmaCascPion || bachelor.pt() < 0.4) && (std::fabs(posdau.tofNSigmaPr()) < cfgtofNSigmaCascProton || posdau.pt() < 0.4) && (std::fabs(negdau.tofNSigmaPi()) < cfgtofNSigmaCascPion || negdau.pt() < 0.4)))) {
        registry.fill(HIST("InvMassXi_all"), casc.pt(), casc.mXi(), casc.eta(), cent);
        if (!setCurrentParticleWeights(weff, wacc, casc, vtxz, 3))
          continue;
        isXi = true;
        candNumAll[2] = candNumAll[2] + 1;
      } else if (casc.sign() < 0 && (casc.mXi() > 1.30) && (casc.mXi() < 1.37) && std::fabs(casc.yXi()) < cfgCasc_rapidity &&
                 (!cfgcheckDauTPC || (std::fabs(bachelor.tpcNSigmaPi()) < cfgtpcNSigmaCascPion && std::fabs(negdau.tpcNSigmaPr()) < cfgtpcNSigmaCascProton && std::fabs(posdau.tpcNSigmaPi()) < cfgtpcNSigmaCascPion)) &&
                 (!cfgcheckDauTOF || ((std::fabs(bachelor.tofNSigmaPi()) < cfgtofNSigmaCascPion || bachelor.pt() < 0.4) && (std::fabs(negdau.tofNSigmaPr()) < cfgtofNSigmaCascProton || negdau.pt() < 0.4) && (std::fabs(posdau.tofNSigmaPi()) < cfgtofNSigmaCascPion || posdau.pt() < 0.4)))) {
        registry.fill(HIST("InvMassXi_all"), casc.pt(), casc.mXi(), casc.eta(), cent);
        if (!setCurrentParticleWeights(weff, wacc, casc, vtxz, 3))
          continue;
        isXi = true;
        candNumAll[2] = candNumAll[2] + 1;
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
      if (std::fabs(casc.mLambda() - o2::constants::physics::MassLambda0) > cfgcasc_mlambdawindow)
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
      if (isOmega) {
        candNum[3] = candNum[3] + 1;
        registry.fill(HIST("hEtaPhiVtxzPOIOmega"), casc.phi(), casc.eta(), vtxz, wacc);
        registry.fill(HIST("InvMassOmega"), casc.pt(), casc.mOmega(), casc.eta(), cent);
        fGFW->Fill(casc.eta(), fXiPtAxis->FindBin(casc.pt()) - 1 + ((fOmegaMass->FindBin(casc.mOmega()) - 1) * nXiPtBins), casc.phi(), wacc * weff, 4);
        if (cfgOutputNUAWeights)
          fWeightsOmega->fill(casc.phi(), casc.eta(), vtxz, casc.pt(), cent, 0);
      }
      if (isXi) {
        candNum[2] = candNum[2] + 1;
        registry.fill(HIST("hEtaPhiVtxzPOIXi"), casc.phi(), casc.eta(), vtxz, wacc);
        registry.fill(HIST("InvMassXi"), casc.pt(), casc.mXi(), casc.eta(), cent);
        fGFW->Fill(casc.eta(), fXiPtAxis->FindBin(casc.pt()) - 1 + ((fXiMass->FindBin(casc.mXi()) - 1) * nXiPtBins), casc.phi(), wacc * weff, 2);
        if (cfgOutputNUAWeights)
          fWeightsXi->fill(casc.phi(), casc.eta(), vtxz, casc.pt(), cent, 0);
      }
    }
    for (int i = 0; i < 4; i++) {
      if (candNumAll[i] > 0) {
        registry.fill(HIST("hEventCount"), 2.5, i + 0.5);
      }
      if (candNum[i] > 0) {
        registry.fill(HIST("hEventCount"), 3.5, i + 0.5);
      }
    }
    // Filling cumulant with ROOT TProfile and loop for all ptBins
    fillProfile(corrconfigs.at(15), HIST("c22"), cent);
    fillProfile(corrconfigs.at(16), HIST("c24"), cent);
    fillProfile(corrconfigs.at(17), HIST("K0sc22"), cent);
    fillProfile(corrconfigs.at(18), HIST("K0sc22"), cent);
    fillProfile(corrconfigs.at(19), HIST("Lambdac22"), cent);
    fillProfile(corrconfigs.at(20), HIST("Lambdac22"), cent);
    for (int i = 1; i <= nPtBins; i++) {
      fillProfilepT(corrconfigs.at(0), HIST("c22dpt"), i, cent);
      fillProfilepT(corrconfigs.at(1), HIST("c22dpt"), i, cent);
      fillProfilepT(corrconfigs.at(2), HIST("c24dpt"), i, cent);
    }
    for (int i = 1; i <= nV0PtBins; i++) {
      fillProfilepTMass(corrconfigs.at(9), HIST("K0sc22dpt"), i, kK0Short, cent);
      fillProfilepTMass(corrconfigs.at(10), HIST("K0sc22dpt"), i, kK0Short, cent);
      fillProfilepTMass(corrconfigs.at(11), HIST("K0sc24dpt"), i, kK0Short, cent);
      fillProfilepTMass(corrconfigs.at(12), HIST("Lambdac22dpt"), i, kLambda0, cent);
      fillProfilepTMass(corrconfigs.at(13), HIST("Lambdac22dpt"), i, kLambda0, cent);
      fillProfilepTMass(corrconfigs.at(14), HIST("Lambdac24dpt"), i, kLambda0, cent);
    }
    for (int i = 1; i <= nXiPtBins; i++) {
      fillProfilepTMass(corrconfigs.at(3), HIST("Xic22dpt"), i, kXiMinus, cent);
      fillProfilepTMass(corrconfigs.at(4), HIST("Xic22dpt"), i, kXiMinus, cent);
      fillProfilepTMass(corrconfigs.at(5), HIST("Xic24dpt"), i, kXiMinus, cent);
      fillProfilepTMass(corrconfigs.at(6), HIST("Omegac22dpt"), i, kOmegaMinus, cent);
      fillProfilepTMass(corrconfigs.at(7), HIST("Omegac22dpt"), i, kOmegaMinus, cent);
      fillProfilepTMass(corrconfigs.at(8), HIST("Omegac24dpt"), i, kOmegaMinus, cent);
    }
    // Fill subevents flow
    TRandom3* fRdm = new TRandom3(0);
    double eventrdm = 10 * fRdm->Rndm();
    for (int j = 1; j <= 10; j++) {
      if (eventrdm > (j - 1) && eventrdm < j)
        continue;
      fillProfile(corrconfigs.at(15), refc22[j - 1], cent);
      fillProfile(corrconfigs.at(16), refc24[j - 1], cent);
      for (int i = 1; i <= nV0PtBins; i++) {
        fillProfilepTMass(corrconfigs.at(9), k0sc22[j - 1], i, kK0Short, cent);
        fillProfilepTMass(corrconfigs.at(10), k0sc22[j - 1], i, kK0Short, cent);
        fillProfilepTMass(corrconfigs.at(11), k0sc24[j - 1], i, kK0Short, cent);
        fillProfilepTMass(corrconfigs.at(12), lambdac22[j - 1], i, kLambda0, cent);
        fillProfilepTMass(corrconfigs.at(13), lambdac22[j - 1], i, kLambda0, cent);
        fillProfilepTMass(corrconfigs.at(14), lambdac24[j - 1], i, kLambda0, cent);
      }
      for (int i = 1; i <= nXiPtBins; i++) {
        fillProfilepTMass(corrconfigs.at(3), xic22[j - 1], i, kXiMinus, cent);
        fillProfilepTMass(corrconfigs.at(4), xic22[j - 1], i, kXiMinus, cent);
        fillProfilepTMass(corrconfigs.at(5), xic24[j - 1], i, kXiMinus, cent);
        fillProfilepTMass(corrconfigs.at(6), omegac22[j - 1], i, kOmegaMinus, cent);
        fillProfilepTMass(corrconfigs.at(7), omegac22[j - 1], i, kOmegaMinus, cent);
        fillProfilepTMass(corrconfigs.at(8), omegac24[j - 1], i, kOmegaMinus, cent);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowGFWOmegaXi>(cfgc)};
}
