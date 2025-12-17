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

/// \file   pidFlowPtCorr.cxx
/// \author Fuchun Cui(fcui@cern.ch) Qiuyu Xia(qiuyu.xia@cern.ch)
/// \since  Nov/24/2025
/// \brief  This task is to caculate vn-[pt] correlation of PID particles

#include "GFW.h"
#include "GFWCumulant.h"
#include "GFWPowerArray.h"
#include "GFWWeights.h"

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

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct PidFlowPtCorr {

  O2_DEFINE_CONFIGURABLE(cfgCutVertex, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutChi2prTPCcls, float, 2.5, "Chi2 per TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgDeltaPhiLocDen, int, 3, "Number of delta phi for local density, 200 bins in 2 pi")

  struct : ConfigurableGroup {
    std::string prefix = "trkQualityOpts";
    // track selections
    O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta range for tracks")
    O2_DEFINE_CONFIGURABLE(cfgRangeEta, float, 0.4f, "Eta range for mean Pt")
    O2_DEFINE_CONFIGURABLE(cfgCutPtMin, float, 0.2f, "Minimal pT for ref tracks")
    O2_DEFINE_CONFIGURABLE(cfgCutPtMax, float, 10.0f, "Maximal pT for ref tracks")
    // track quality selections for daughter track
    O2_DEFINE_CONFIGURABLE(cfgITSNCls, int, 5, "check minimum number of ITS clusters")
    O2_DEFINE_CONFIGURABLE(cfgTPCNCls, int, 50, "check minimum number of TPC hits")
    O2_DEFINE_CONFIGURABLE(cfgTPCCrossedRows, int, 70, "check minimum number of TPC crossed rows")
    O2_DEFINE_CONFIGURABLE(cfgITSChi2NDF, double, 2.5, "check ITS Chi2NDF")
    O2_DEFINE_CONFIGURABLE(cfgCheckGlobalTrack, bool, false, "check global track")
  } trkQualityOpts;

  struct : ConfigurableGroup {
    std::string prefix = "evtSelOpts";
    O2_DEFINE_CONFIGURABLE(cfgDoTVXinTRD, bool, false, "check kTVXinTRD")
    O2_DEFINE_CONFIGURABLE(cfgDoNoTimeFrameBorder, bool, true, "check kNoTimeFrameBorder")
    O2_DEFINE_CONFIGURABLE(cfgDoNoITSROFrameBorder, bool, true, "check kNoITSROFrameBorder")
    O2_DEFINE_CONFIGURABLE(cfgDoNoSameBunchPileup, bool, true, "check kNoITSROFrameBorder")
    O2_DEFINE_CONFIGURABLE(cfgDoIsGoodZvtxFT0vsPV, bool, true, "check kIsGoodZvtxFT0vsPV")
    O2_DEFINE_CONFIGURABLE(cfgDoNoCollInTimeRangeStandard, bool, true, "check kNoCollInTimeRangeStandard")
    O2_DEFINE_CONFIGURABLE(cfgDoIsGoodITSLayersAll, bool, true, "check kIsGoodITSLayersAll")
    O2_DEFINE_CONFIGURABLE(cfgCutOccupancyHigh, int, 3000, "High cut on TPC occupancy")
    O2_DEFINE_CONFIGURABLE(cfgDoMultPVCut, bool, true, "do multNTracksPV vs cent cut")
    O2_DEFINE_CONFIGURABLE(cfgMultPVCut, std::vector<float>, (std::vector<float>{3074.43, -106.192, 1.46176, -0.00968364, 2.61923e-05, 182.128, -7.43492, 0.193901, -0.00256715, 1.22594e-05}), "Used MultPVCut function parameter")
    O2_DEFINE_CONFIGURABLE(cfgDoV0AT0Acut, bool, true, "do V0A-T0A cut")
    O2_DEFINE_CONFIGURABLE(cfgCutminIR, float, -1, "cut min IR")
    O2_DEFINE_CONFIGURABLE(cfgCutmaxIR, float, 3000, "cut max IR")
  } evtSeleOpts;

  O2_DEFINE_CONFIGURABLE(cfgCasc_rapidity, float, 0.5, "rapidity")
  O2_DEFINE_CONFIGURABLE(cfgNSigmapid, std::vector<float>, (std::vector<float>{3, 3, 3, 9, 9, 9, 9, 9, 9}), "tpc, tof and its NSigma for Pion Proton Kaon")
  O2_DEFINE_CONFIGURABLE(cfgMeanPtcent, std::vector<float>, (std::vector<float>{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}), "mean Pt in different cent bin")
  O2_DEFINE_CONFIGURABLE(cfgAcceptancePath, std::vector<std::string>, (std::vector<std::string>{"Users/f/fcui/NUA/NUAREFPartical", "Users/f/fcui/NUA/NUAK0s", "Users/f/fcui/NUA/NUALambda", "Users/f/fcui/NUA/NUAXi", "Users/f/fcui/NUA/NUAOmega"}), "CCDB path to acceptance object")
  O2_DEFINE_CONFIGURABLE(cfgEfficiencyPath, std::vector<std::string>, (std::vector<std::string>{"PathtoRef"}), "CCDB path to efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgRunNumbers, std::vector<int>, (std::vector<int>{544095, 544098, 544116, 544121, 544122, 544123, 544124}), "Preconfigured run numbers")

  // switch
  O2_DEFINE_CONFIGURABLE(cfgDoAccEffCorr, bool, false, "do acc and eff corr")
  O2_DEFINE_CONFIGURABLE(cfgDoLocDenCorr, bool, false, "do local density corr")
  O2_DEFINE_CONFIGURABLE(cfgOutputNUAWeights, bool, false, "Fill and output NUA weights")
  O2_DEFINE_CONFIGURABLE(cfgOutputrunbyrun, bool, false, "Fill and output NUA weights run by run")
  O2_DEFINE_CONFIGURABLE(cfgOutPutMC, bool, false, "Fill MC graphs, note that if the processMCgen is open,this MUST be open")
  O2_DEFINE_CONFIGURABLE(cfgOutputLocDenWeights, bool, false, "Fill and output local density weights")
  O2_DEFINE_CONFIGURABLE(cfgOutputQA, bool, false, "do QA")

  /**
   * @brief cfg for PID pt range
   * @details default datas are from run2, note that separate pi-k and k-p needs to use difference pt range
   */
  // separate pi and k
  O2_DEFINE_CONFIGURABLE(cfgPtMin4ITS_pi_k, float, 0.2, "pt min for ITS to separate pi and k");
  O2_DEFINE_CONFIGURABLE(cfgPtMax4ITS_pi_k, float, 0.8, "pt max for ITS to separate pi and k");
  O2_DEFINE_CONFIGURABLE(cfgPtMin4TOF_pi_k, float, 0.5, "pt min for TOF to separate pi and k");
  O2_DEFINE_CONFIGURABLE(cfgPtMax4TOF_pi_k, float, 3.0, "pt max for TOF to separate pi and k");
  O2_DEFINE_CONFIGURABLE(cfgPtMin4TPC_pi_k, float, 0.2, "pt min for TPC to separate pi and k");
  O2_DEFINE_CONFIGURABLE(cfgPtMax4TPC_pi_k, float, 0.6, "pt max for TPC to separate pi and k");
  // end separate pi and k

  // separate k-p
  O2_DEFINE_CONFIGURABLE(cfgPtMin4ITS_k_p, float, 0.4, "pt min for ITS to separate k and p");
  O2_DEFINE_CONFIGURABLE(cfgPtMax4ITS_k_p, float, 1.4, "pt max for ITS to separate k and p");
  O2_DEFINE_CONFIGURABLE(cfgPtMin4TOF_k_p, float, 0.6, "pt min for TOF to separate k and p");
  O2_DEFINE_CONFIGURABLE(cfgPtMax4TOF_k_p, float, 3.0, "pt max for TOF to separate k and p");
  O2_DEFINE_CONFIGURABLE(cfgPtMin4TPC_k_p, float, 0.2, "pt min for TPC to separate k and p");
  O2_DEFINE_CONFIGURABLE(cfgPtMax4TPC_k_p, float, 1.0, "pt max for TPC to separate k and p");
  // end separate k-p
  // end cfg for PID pt range

  ConfigurableAxis cfgaxisVertex{"cfgaxisVertex", {20, -10, 10}, "vertex axis for histograms"};
  ConfigurableAxis cfgaxisPhi{"cfgaxisPhi", {60, 0.0, constants::math::TwoPI}, "phi axis for histograms"};
  ConfigurableAxis cfgaxisEta{"cfgaxisEta", {40, -1., 1.}, "eta axis for histograms"};
  ConfigurableAxis cfgaxisPt{"cfgaxisPt", {VARIABLE_WIDTH, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.50, 4.00, 4.50, 5.00, 5.50, 6.00, 10.0}, "pt (GeV)"};
  ConfigurableAxis cfgaxisMeanPt{"cfgaxisMeanPt", {300, 0, 3}, "pt (GeV)"};
  ConfigurableAxis cfgaxisNch{"cfgaxisNch", {3000, 0.5, 3000.5}, "Nch"};
  ConfigurableAxis cfgaxisLocalDensity{"cfgaxisLocalDensity", {200, 0, 600}, "local density"};
  ConfigurableAxis cfgaxisRun{"cfgaxisRun", {7, 0, 7}, "axis of runs in the data"};
  Configurable<std::vector<double>> cfgTrackDensityP0{"cfgTrackDensityP0", std::vector<double>{0.7217476707, 0.7384792571, 0.7542625668, 0.7640680200, 0.7701951667, 0.7755299053, 0.7805901710, 0.7849446786, 0.7957356586, 0.8113039262, 0.8211968966, 0.8280558878, 0.8329342135}, "parameter 0 for track density efficiency correction"};
  Configurable<std::vector<double>> cfgTrackDensityP1{"cfgTrackDensityP1", std::vector<double>{-2.169488e-05, -2.191913e-05, -2.295484e-05, -2.556538e-05, -2.754463e-05, -2.816832e-05, -2.846502e-05, -2.843857e-05, -2.705974e-05, -2.477018e-05, -2.321730e-05, -2.203315e-05, -2.109474e-05}, "parameter 1 for track density efficiency correction"};

  AxisSpec axisMultiplicity{{0, 5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90}, "Centrality (%)"};

  // filter and using
  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = (nabs(aod::track::eta) < trkQualityOpts.cfgCutEta.value) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && (aod::track::tpcChi2NCl < cfgCutChi2prTPCcls);

  using TracksPID = soa::Join<aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr>;
  using AodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, o2::aod::TrackSelectionExtension, aod::TracksExtra, TracksPID, aod::TracksIU>>; // tracks filter
  using AodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::MultsRun3>>;                                               // collisions filter

  Filter mccollisionFilter = nabs(aod::mccollision::posZ) < cfgCutVertex;
  using FilteredMcCollisions = soa::Filtered<aod::McCollisions>;

  Filter particleFilter = nabs(aod::mcparticle::eta) < trkQualityOpts.cfgCutEta.value;
  using FilteredMcParticles = soa::Filtered<soa::Join<aod::McParticles, aod::ParticlesToTracks>>;
  // end using and filter

  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;
  ctpRateFetcher rateFetcher;
  O2_DEFINE_CONFIGURABLE(cfgnolaterthan, int64_t, std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object")
  O2_DEFINE_CONFIGURABLE(cfgurl, std::string, "http://alice-ccdb.cern.ch", "url of the ccdb repository")

  // Define output
  HistogramRegistry registry{"registry"};
  OutputObj<GFWWeights> fWeightsREF{GFWWeights("weightsREF")};

  // define global variables
  GFW* fGFW = new GFW(); // GFW class used from main src
  std::vector<GFW::CorrConfig> corrconfigs;
  std::vector<std::string> cfgAcceptance;
  std::vector<std::string> cfgEfficiency;
  std::vector<float> cfgMultPVCutPara;
  std::vector<float> cfgNSigma;
  std::vector<float> cfgMeanPt;
  std::vector<int> runNumbers;
  std::map<int, std::vector<std::shared_ptr<TH1>>> th1sList;
  std::map<int, std::vector<std::shared_ptr<TH3>>> th3sList;
  enum OutputTH1Names {
    // here are TProfiles for vn-pt correlations that are not implemented in GFW
    hPhi = 0,
    hPhicorr,
    kCount_TH1Names
  };

  enum OutputTH3Names {
    hPhiEtaVtxz = 0,
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

  TAxis* fMultAxis = nullptr;

  std::vector<TF1*> funcEff;
  TH1D* hFindPtBin;
  TF1* funcV2;
  TF1* funcV3;
  TF1* funcV4;

  void init(InitContext const&) // Initialization
  {
    ccdb->setURL(cfgurl.value);
    ccdb->setCaching(true);
    ccdb->setCreatedNotAfter(cfgnolaterthan.value);

    cfgAcceptance = cfgAcceptancePath;
    cfgEfficiency = cfgEfficiencyPath;
    cfgNSigma = cfgNSigmapid;
    cfgMeanPt = cfgMeanPtcent;
    cfgMultPVCutPara = evtSeleOpts.cfgMultPVCut;

    // Set the pt, mult and phi Axis;
    o2::framework::AxisSpec axisPt = cfgaxisPt;
    nPtBins = axisPt.binEdges.size() - 1;
    fPtAxis = new TAxis(nPtBins, &(axisPt.binEdges)[0]);

    o2::framework::AxisSpec axisMult = axisMultiplicity;
    int nMultBins = axisMult.binEdges.size() - 1;
    fMultAxis = new TAxis(nMultBins, &(axisMult.binEdges)[0]);

    // Add some output objects to the histogram registry
    registry.add("hPhi", "", {HistType::kTH1D, {cfgaxisPhi}});
    registry.add("hPhicorr", "", {HistType::kTH1D, {cfgaxisPhi}});
    registry.add("hEta", "", {HistType::kTH1D, {cfgaxisEta}});
    registry.add("hVtxZ", "", {HistType::kTH1D, {cfgaxisVertex}});
    registry.add("hMult", "", {HistType::kTH1D, {cfgaxisNch}});
    registry.add("hMultTPC", "", {HistType::kTH1D, {cfgaxisNch}});
    registry.add("hCent", "", {HistType::kTH1D, {{90, 0, 90}}});
    registry.add("hCentvsNch", "", {HistType::kTH2D, {{18, 0, 90}, cfgaxisNch}});
    registry.add("MC/hCentvsNchMC", "", {HistType::kTH2D, {{18, 0, 90}, cfgaxisNch}});
    registry.add("hCentvsMultTPC", "", {HistType::kTH2D, {{18, 0, 90}, cfgaxisNch}});
    registry.add("MC/hCentvsMultTPCMC", "", {HistType::kTH2D, {{18, 0, 90}, cfgaxisNch}});
    registry.add("hPt", "", {HistType::kTH1D, {cfgaxisPt}});
    registry.add("hEtaPhiVtxzREF", "", {HistType::kTH3D, {cfgaxisPhi, cfgaxisEta, {20, -10, 10}}});
    registry.add("hNTracksPVvsCentrality", "", {HistType::kTH2D, {{5000, 0, 5000}, axisMultiplicity}});

    // TPC vs TOF vs its, comparation graphs, check the PID performance in difference pt
    if (cfgOutputQA) {
      registry.add("DetectorPidPerformace/TPCvsTOF/Pi", "", {HistType::kTH3D, {{600, -30, 30}, {600, -30, 30}, cfgaxisPt}});
      registry.add("DetectorPidPerformace/TPCvsTOF/Pr", "", {HistType::kTH3D, {{600, -30, 30}, {600, -30, 30}, cfgaxisPt}});
      registry.add("DetectorPidPerformace/TPCvsTOF/Ka", "", {HistType::kTH3D, {{600, -30, 30}, {600, -30, 30}, cfgaxisPt}});

      registry.add("DetectorPidPerformace/TPCvsITS/Pi", "", {HistType::kTH3D, {{600, -30, 30}, {600, -30, 30}, cfgaxisPt}});
      registry.add("DetectorPidPerformace/TPCvsITS/Pr", "", {HistType::kTH3D, {{600, -30, 30}, {600, -30, 30}, cfgaxisPt}});
      registry.add("DetectorPidPerformace/TPCvsITS/Ka", "", {HistType::kTH3D, {{600, -30, 30}, {600, -30, 30}, cfgaxisPt}});

      registry.add("DetectorPidPerformace/ITSvsTOF/Pi", "", {HistType::kTH3D, {{600, -30, 30}, {600, -30, 30}, cfgaxisPt}});
      registry.add("DetectorPidPerformace/ITSvsTOF/Pr", "", {HistType::kTH3D, {{600, -30, 30}, {600, -30, 30}, cfgaxisPt}});
      registry.add("DetectorPidPerformace/ITSvsTOF/Ka", "", {HistType::kTH3D, {{600, -30, 30}, {600, -30, 30}, cfgaxisPt}});
    } // end TPC vs TOF vs its graph add

    if (cfgOutPutMC) {
      // hist for NUE correction
      registry.add("correction/hCentPtMC", "", {HistType::kTH2D, {axisMultiplicity, cfgaxisPt}});
      registry.add("correction/hCentPtData", "", {HistType::kTH2D, {axisMultiplicity, cfgaxisPt}});
    } // cfgoutputMC

    runNumbers = cfgRunNumbers;
    if (cfgOutputrunbyrun) {
      // hist for NUA
      registry.add("correction/hRunNumberPhiEtaVertex", "", {HistType::kTHnSparseF, {cfgaxisRun, cfgaxisPhi, cfgaxisEta, cfgaxisVertex}});
      // set "correction/hRunNumberPhiEtaVertex" axis0 label
      for (uint64_t idx = 1; idx <= runNumbers.size(); idx++) {
        registry.get<THnSparse>(HIST("correction/hRunNumberPhiEtaVertex"))->GetAxis(0)->SetBinLabel(idx, std::to_string(runNumbers[idx - 1]).c_str());
      }
      // end set "correction/hRunNumberPhiEtaVertex" axis0 label
    } // cfgooutputrunbyrun

    // set bin label for hEventCount
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
    // end set bin label for eventcount

    // cumulant of flow
    registry.add("c22", ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("c32", ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("c24", ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("c34", ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("c22Full", ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("c22TrackWeight", ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("c32TrackWeight", ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("c24TrackWeight", ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("c34TrackWeight", ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("c22FullTrackWeight", ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile, {axisMultiplicity}});

    registry.add("pi/c22", ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("ka/c22", ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("pr/c22", ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("pi/c24", ";Centrality  (%) ; C_{2}{4} ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("ka/c24", ";Centrality  (%) ; C_{2}{4} ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("pr/c24", ";Centrality  (%) ; C_{2}{4} ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("pi/c32", ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("ka/c32", ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("pr/c32", ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("pi/c34", ";Centrality  (%) ; C_{2}{4} ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("ka/c34", ";Centrality  (%) ; C_{2}{4} ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("pr/c34", ";Centrality  (%) ; C_{2}{4} ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("pi/c22PiPi", ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("ka/c22KaKa", ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("pr/c22PrPr", ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("pi/c32PiPi", ";Centrality  (%) ; C_{2}{4} ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("ka/c32KaKa", ";Centrality  (%) ; C_{2}{4} ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("pr/c32PrPr", ";Centrality  (%) ; C_{2}{4} ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("pi/c22TrackWeight", ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("ka/c22TrackWeight", ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("pr/c22TrackWeight", ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("pi/c32TrackWeight", ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("ka/c32TrackWeight", ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("pr/c32TrackWeight", ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile, {axisMultiplicity}});

    // vn-pt corr
    registry.add("covV2Pt", ";Centrality  (%) ; cov(v_{2}^{2}{2}, P_{T}) ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("pi/covV2Pt", ";Centrality  (%) ; cov(v_{2}^{2}{2}, P_{T}) ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("ka/covV2Pt", ";Centrality  (%) ; cov(v_{2}^{2}{2}, P_{T}) ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("pr/covV2Pt", ";Centrality  (%) ; cov(v_{2}^{2}{2}, P_{T}) ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("covV3Pt", ";Centrality  (%) ; cov(v_{2}^{2}{2}, P_{T}) ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("pi/covV3Pt", ";Centrality  (%) ; cov(v_{2}^{2}{2}, P_{T}) ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("ka/covV3Pt", ";Centrality  (%) ; cov(v_{2}^{2}{2}, P_{T}) ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("pr/covV3Pt", ";Centrality  (%) ; cov(v_{2}^{2}{2}, P_{T}) ", {HistType::kTProfile, {axisMultiplicity}});

    registry.add("covV2Pt_diffpt", ";Centrality  (%) ; cov(v_{2}^{2}{2}, P_{T}) ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("pi/covV2Pt_diffpt", ";Centrality  (%) ; cov(v_{2}^{2}{2}, P_{T}) ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("ka/covV2Pt_diffpt", ";Centrality  (%) ; cov(v_{2}^{2}{2}, P_{T}) ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("pr/covV2Pt_diffpt", ";Centrality  (%) ; cov(v_{2}^{2}{2}, P_{T}) ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("covV3Pt_diffpt", ";Centrality  (%) ; cov(v_{2}^{2}{2}, P_{T}) ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("pi/covV3Pt_diffpt", ";Centrality  (%) ; cov(v_{2}^{2}{2}, P_{T}) ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("ka/covV3Pt_diffpt", ";Centrality  (%) ; cov(v_{2}^{2}{2}, P_{T}) ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("pr/covV3Pt_diffpt", ";Centrality  (%) ; cov(v_{2}^{2}{2}, P_{T}) ", {HistType::kTProfile, {axisMultiplicity}});

    registry.add("hMeanPt", ";Centrality  (%) ; [P_{T}]} ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("ptSquareAve", ";Centrality  (%) ; [P_{T}^{2}] ", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("ptAve", ";Centrality  (%) ; [P_{T}] ", {HistType::kTProfile, {axisMultiplicity}});

    registry.add("c22dmeanpt", ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile2D, {axisMultiplicity, cfgaxisMeanPt}});
    registry.add("pi/c22dmeanpt", ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile2D, {axisMultiplicity, cfgaxisMeanPt}});
    registry.add("ka/c22dmeanpt", ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile2D, {axisMultiplicity, cfgaxisMeanPt}});
    registry.add("pr/c22dmeanpt", ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile2D, {axisMultiplicity, cfgaxisMeanPt}});

    // Data
    fGFW->AddRegion("reffull", -0.8, 0.8, 1, 1); // ("name", etamin, etamax, ptbinnum, bitmask)eta region -0.8 to 0.8
    fGFW->AddRegion("refN08", -0.8, -0.4, 1, 1);
    fGFW->AddRegion("refP08", 0.4, 0.8, 1, 1);
    fGFW->AddRegion("refN", -0.8, 0, 1, 1);
    fGFW->AddRegion("refP", 0, 0.8, 1, 1);

    fGFW->AddRegion("poiPiN08", -0.8, -0.4, 1, 2);
    fGFW->AddRegion("poiPiP08", 0.4, 0.8, 1, 2);
    fGFW->AddRegion("poiPiN", -0.8, 0, 1, 2);
    fGFW->AddRegion("poiPiP", 0, 0.8, 1, 2);
    fGFW->AddRegion("olPiN", -0.8, 0, 1, 16);
    fGFW->AddRegion("olPiP", 0, 0.8, 1, 16);

    fGFW->AddRegion("poiKaN08", -0.8, -0.4, 1, 4);
    fGFW->AddRegion("poiKaP08", 0.4, 0.8, 1, 4);
    fGFW->AddRegion("poiKaN", -0.8, 0, 1, 4);
    fGFW->AddRegion("poiKaP", 0, 0.8, 1, 4);
    fGFW->AddRegion("olKaN", -0.8, 0, 1, 32);
    fGFW->AddRegion("olKaP", 0, 0.8, 1, 32);

    fGFW->AddRegion("poiPrN08", -0.8, -0.4, 1, 8);
    fGFW->AddRegion("poiPrP08", 0.4, 0.8, 1, 8);
    fGFW->AddRegion("poiPrN", -0.8, 0, 1, 8);
    fGFW->AddRegion("poiPrP", 0, 0.8, 1, 8);
    fGFW->AddRegion("olPrN", -0.8, 0, 1, 64);
    fGFW->AddRegion("olPrP", 0, 0.8, 1, 64);

    // pushback
    // Data
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refP08 {2} refN08 {-2}", "Ref08Gap22", kFALSE)); // 0
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN {2 2} refP {-2 -2}", "Ref0Gap24", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN {2} refP {-2}", "Ref0Gap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refP08 {3} refN08 {-3}", "Ref08Gap32", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refP08 {3 3} refN08 {-3 -3}", "Ref08Gap34", kFALSE));

    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPiN08 {2} refP08 {-2}", "Pion08gap22a", kFALSE)); // 5
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPiP08 {2} refN08 {-2}", "Pion08gap22b", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiKaN08 {2} refP08 {-2}", "Kaon08gap22a", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiKaP08 {2} refN08 {-2}", "Kaon08gap22b", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPrN08 {2} refP08 {-2}", "Prot08gap22a", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPrP08 {2} refN08 {-2}", "Prot08gap22b", kFALSE)); // 10
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPiN refN | olPiN {2 2} refP {-2 -2}", "Pion0gap24a", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPiP refP | olPiP {2 2} refN {-2 -2}", "Pion0gap24b", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiKaN refN | olKaN {2 2} refP {-2 -2}", "Kaon0gap24a", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiKaP refP | olKaP {2 2} refN {-2 -2}", "Kaon0gap24b", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPrN refN | olPrN {2 2} refP {-2 -2}", "Prot0gap24a", kFALSE)); // 15
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPrP refP | olPaP {2 2} refN {-2 -2}", "Prot0gap24b", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPiN08 {3} refP08 {-3}", "Pion08gap32a", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPiP08 {3} refN08 {-3}", "Pion08gap32b", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiKaN08 {3} refP08 {-3}", "Kaon08gap32a", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiKaP08 {3} refN08 {-3}", "Kaon08gap32b", kFALSE)); // 20
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPrN08 {3} refP08 {-3}", "Prot08gap32a", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPrP08 {3} refN08 {-3}", "Prot08gap32b", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPiN refN | olPiN {3 3} refP {-3 -3}", "Pion0gap34a", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPiP refP | olPiP {3 3} refN {-3 -3}", "Pion0gap34b", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiKaN refN | olKaN {3 3} refP {-3 -3}", "Kaon0gap34a", kFALSE)); // 25
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiKaP refP | olKaP {3 3} refN {-3 -3}", "Kaon0gap34b", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPrN refN | olPrN {3 3} refP {-3 -3}", "Prot0gap34a", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPrP refP | olPaP {3 3} refN {-3 -3}", "Prot0gap34b", kFALSE));

    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPiN08 {2} poiPiP08 {-2}", "PiPi08gap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiKaN08 {2} poiKaP08 {-2}", "KaKa08gap22", kFALSE)); // 30
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPrN08 {2} poiPrP08 {-2}", "PrPr08gap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPiN08 {3} poiPiP08 {-3}", "PiPi08gap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiKaN08 {3} poiKaP08 {-3}", "KaKa08gap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPrN08 {3} poiPrP08 {-3}", "PrPr08gap22", kFALSE));

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
    }

    std::vector<double> pTEffBins = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.4, 1.8, 2.2, 2.6, 3.0};
    hFindPtBin = new TH1D("hFindPtBin", "hFindPtBin", pTEffBins.size() - 1, &pTEffBins[0]);
    funcEff.resize(pTEffBins.size() - 1);
    // LHC24g3 Eff
    std::vector<double> f1p0 = cfgTrackDensityP0;
    std::vector<double> f1p1 = cfgTrackDensityP1;
    for (uint ifunc = 0; ifunc < pTEffBins.size() - 1; ifunc++) {
      funcEff[ifunc] = new TF1(Form("funcEff%i", ifunc), "[0]+[1]*x", 0, 3000);
      funcEff[ifunc]->SetParameters(f1p0[ifunc], f1p1[ifunc]);
    }
    funcV2 = new TF1("funcV2", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 100);
    funcV2->SetParameters(0.0186111, 0.00351907, -4.38264e-05, 1.35383e-07, -3.96266e-10);
    funcV3 = new TF1("funcV3", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 100);
    funcV3->SetParameters(0.0174056, 0.000703329, -1.45044e-05, 1.91991e-07, -1.62137e-09);
    funcV4 = new TF1("funcV4", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 100);
    funcV4->SetParameters(0.008845, 0.000259668, -3.24435e-06, 4.54837e-08, -6.01825e-10);
  }

  /**
   * @brief Identify whether the input track is a Pion
   *
   * @param track Input track object to be identified
   * @return true  The track is identified as Pion
   * @return false The track is NOT identified as Pion
   * @note The result is the logical AND of valid detector Sigma checks.
   *       If the track pt is out of a detector's valid range, that detector's check is skipped.
   */
  template <typename TrackObject>
  bool isPion(TrackObject track)
  {
    bool resultPion = true;

    // Declare ITSResponse object internally to get ITS Sigma
    o2::aod::ITSResponse itsResponse;
    // Extract sigma values and pt from track
    const float itsSigma = std::fabs(itsResponse.nSigmaITS<o2::track::PID::Pion>(track));
    const float tofSigma = std::fabs(track.tofNSigmaPi());
    const float tpcSigma = std::fabs(track.tpcNSigmaPi());
    const float pt = track.pt();

    // ITS detector check (pi-k separation pt range)
    if (pt > cfgPtMin4ITS_pi_k && pt < cfgPtMax4ITS_pi_k) {
      resultPion &= (itsSigma < cfgNSigma[6]);
    }
    // end ITS

    // TOF detector check (pi-k separation pt range)
    if (pt > cfgPtMin4TOF_pi_k && pt < cfgPtMax4TOF_pi_k) {
      resultPion &= (tofSigma < cfgNSigma[3]);
    }
    // end TOF

    // TPC detector check (pi-k separation pt range)
    if (pt > cfgPtMin4TPC_pi_k && pt < cfgPtMax4TPC_pi_k) {
      resultPion &= (tpcSigma < cfgNSigma[0]);
    }
    // end TPC

    return resultPion;
  }

  /**
   * @brief Identify whether the input track is a Proton
   *
   * @param track Input track object to be identified
   * @return true  The track is identified as Proton
   * @return false The track is NOT identified as Proton
   * @note The result is the logical AND of valid detector Sigma checks.
   *       If the track pt is out of a detector's valid range, that detector's check is skipped.
   */
  template <typename TrackObject>
  bool isProton(TrackObject track)
  {
    bool resultProton = true;

    // Declare ITSResponse object internally to get ITS Sigma
    o2::aod::ITSResponse itsResponse;
    // Extract sigma values and pt from track
    const float itsSigma = std::fabs(itsResponse.nSigmaITS<o2::track::PID::Proton>(track));
    const float tofSigma = std::fabs(track.tofNSigmaPr());
    const float tpcSigma = std::fabs(track.tpcNSigmaPr());
    const float pt = track.pt();

    // ITS detector check (k-p separation pt range)
    if (pt > cfgPtMin4ITS_k_p && pt < cfgPtMax4ITS_k_p) {
      resultProton &= (itsSigma < cfgNSigma[7]);
    }
    // end ITS

    // TOF detector check (k-p separation pt range)
    if (pt > cfgPtMin4TOF_k_p && pt < cfgPtMax4TOF_k_p) {
      resultProton &= (tofSigma < cfgNSigma[4]);
    }
    // end TOF

    // TPC detector check (k-p separation pt range)
    if (pt > cfgPtMin4TPC_k_p && pt < cfgPtMax4TPC_k_p) {
      resultProton &= (tpcSigma < cfgNSigma[1]);
    }
    // end TPC

    return resultProton;
  }

  /**
   * @brief Identify whether the input track is a Kaon (separate from Pion and Proton)
   *
   * @param track Input track object (aod::Track) to be identified
   * @return true  The track is identified as Kaon
   * @return false The track is NOT identified as Kaon
   * @note The result is the logical AND of valid detector Sigma checks.
   *       Only pt range that overlaps with both pi-k and k-p separation ranges is checked.
   *       If the track pt is out of the overlapping range, that detector's check is skipped.
   */
  template <typename TrackObject>
  bool isKaon(TrackObject track)
  {
    bool resultKaon = true;

    // Declare ITSResponse object internally to get ITS Sigma
    o2::aod::ITSResponse itsResponse;
    // Extract sigma values and pt from track
    const float itsSigma = std::fabs(itsResponse.nSigmaITS<o2::track::PID::Kaon>(track));
    const float tofSigma = std::fabs(track.tofNSigmaKa());
    const float tpcSigma = std::fabs(track.tpcNSigmaKa());
    const float pt = track.pt();

    // ITS detector check (overlap of pi-k and k-p separation pt ranges)
    if (pt > cfgPtMin4ITS_k_p && pt > cfgPtMin4ITS_pi_k && pt < cfgPtMax4ITS_k_p && pt < cfgPtMax4ITS_pi_k) {
      resultKaon &= (itsSigma < cfgNSigma[8]);
    }
    // end ITS

    // TOF detector check (overlap of pi-k and k-p separation pt ranges)
    if (pt > cfgPtMin4TOF_k_p && pt > cfgPtMin4TOF_pi_k && pt < cfgPtMax4TOF_k_p && pt < cfgPtMax4TOF_pi_k) {
      resultKaon &= (tofSigma < cfgNSigma[5]);
    }
    // end TOF

    // TPC detector check (overlap of pi-k and k-p separation pt ranges)
    if (pt > cfgPtMin4TPC_k_p && pt > cfgPtMin4TPC_pi_k && pt < cfgPtMax4TPC_k_p && pt < cfgPtMax4TPC_pi_k) {
      resultKaon &= (tpcSigma < cfgNSigma[2]);
    }
    // end TPC

    return resultKaon;
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
      if (std::fabs(val) < 1) {
        registry.fill(tarName, cent, val, dnx);
        return;
      }
    }
    return;
  }

  template <char... chars>
  void fillProfilevnpt(const GFW::CorrConfig& corrconf, const ConstStr<chars...>& tarName, const double& cent, const double& ptSum, const double& nch, const double& meanPt = 0)
  {
    double dnx, val;
    dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
    if (dnx == 0)
      return;
    val = fGFW->Calculate(corrconf, 0, kFALSE).real() / dnx;
    if (std::fabs(val) < 1)
      registry.fill(tarName, cent, val * (ptSum / nch - meanPt), dnx * nch);
    return;
  }

  template <char... chars>
  void fillProfilePOIvnpt(const GFW::CorrConfig& corrconf, const ConstStr<chars...>& tarName, const double& cent, const double& ptSum, const double& nch)
  {
    double dnx, val;
    dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
    if (dnx == 0)
      return;
    val = fGFW->Calculate(corrconf, 0, kFALSE).real() / dnx;

    if (std::fabs(val) < 1)
      registry.fill(tarName, cent, ptSum / nch, val, dnx);
    return;
  }

  void loadCorrections(uint64_t timestamp)
  {
    if (correctionsLoaded)
      return;
    int nspecies = 1;
    if (cfgAcceptance.size() == static_cast<uint64_t>(nspecies)) {
      for (int i = 0; i <= nspecies - 1; i++) {
        mAcceptance.push_back(ccdb->getForTimeStamp<GFWWeights>(cfgAcceptance[i], timestamp));
      }
      if (mAcceptance.size() == static_cast<uint64_t>(nspecies))
        LOGF(info, "Loaded acceptance weights");
      else
        LOGF(warning, "Could not load acceptance weights");
    }
    if (cfgEfficiency.size() == static_cast<uint64_t>(nspecies)) {
      for (int i = 0; i <= nspecies - 1; i++) {
        mEfficiency.push_back(ccdb->getForTimeStamp<TH1D>(cfgEfficiency[i], timestamp));
      }
      if (mEfficiency.size() == static_cast<uint64_t>(nspecies))
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
    int nspecies = 1;
    if (mEfficiency.size() == static_cast<uint64_t>(nspecies))
      eff = mEfficiency[ispecies]->GetBinContent(mEfficiency[ispecies]->FindBin(track.pt()));
    else
      eff = 1.0;
    if (eff == 0)
      return false;
    weight_nue = 1. / eff;
    if (mAcceptance.size() == static_cast<uint64_t>(nspecies))
      weight_nua = mAcceptance[ispecies]->getNUA(track.phi(), track.eta(), vtxz);
    else
      weight_nua = 1;
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

  void processData(AodCollisions::iterator const& collision, aod::BCsWithTimestamps const&, AodTracks const& tracks)
  {
    // init
    int nTot = tracks.size();
    float nMultTPC = collision.multTPC();
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int runNumber = bc.runNumber();
    double interactionRate = rateFetcher.fetch(ccdb.service, bc.timestamp(), runNumber, "ZNC hadronic") * 1.e-3;
    // end init

    // collision cut
    // include : 1.track.size 2.collision.sel8 3.this->evenSelected
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
    // end collision cut

    // correction
    loadCorrections(bc.timestamp());
    float vtxz = collision.posZ();
    registry.fill(HIST("hVtxZ"), vtxz);
    registry.fill(HIST("hMult"), nTot);
    registry.fill(HIST("hMultTPC"), nMultTPC);
    registry.fill(HIST("hCent"), cent);
    // end correction

    double psi2Est = 0, psi3Est = 0, psi4Est = 0;
    float wEPeff = 1;
    double v2 = 0, v3 = 0, v4 = 0;
    if (cfgDoLocDenCorr) {
      double q2x = 0, q2y = 0;
      double q3x = 0, q3y = 0;
      double q4x = 0, q4y = 0;
      for (const auto& track : tracks) {
        // pt cut
        bool withinPtRef = (trkQualityOpts.cfgCutPtMin.value < track.pt()) && (track.pt() < trkQualityOpts.cfgCutPtMax.value); // within RF pT rang
        if (withinPtRef) {
          q2x += std::cos(2 * track.phi());
          q2y += std::sin(2 * track.phi());
          q3x += std::cos(3 * track.phi());
          q3y += std::sin(3 * track.phi());
          q4x += std::cos(4 * track.phi());
          q4y += std::sin(4 * track.phi());
        }
      }
      psi2Est = std::atan2(q2y, q2x) / 2.;
      psi3Est = std::atan2(q3y, q3x) / 3.;
      psi4Est = std::atan2(q4y, q4x) / 4.;
      v2 = funcV2->Eval(cent);
      v3 = funcV3->Eval(cent);
      v4 = funcV4->Eval(cent);
    } // cfgDoLocDenCorr

    float weff = 1;
    float wacc = 1;
    double ptSum = 0;
    double ptSumw2 = 0;
    double nch = 0;
    double nchSquare = 0;
    double ptSquareSum = 0;
    // fill GFW ref flow
    for (const auto& track : tracks) {
      if (cfgDoAccEffCorr) {
        if (!setCurrentParticleWeights(weff, wacc, track, vtxz, 0))
          continue;
      } // cfgDoAccEffCorr
      if (cfgDoLocDenCorr) {
        bool withinPtRef = (trkQualityOpts.cfgCutPtMin.value < track.pt()) && (track.pt() < trkQualityOpts.cfgCutPtMax.value);
        if (withinPtRef) {
          double fphi = v2 * std::cos(2 * (track.phi() - psi2Est)) + v3 * std::cos(3 * (track.phi() - psi3Est)) + v4 * std::cos(4 * (track.phi() - psi4Est));
          fphi = (1 + 2 * fphi);
          int pTBinForEff = hFindPtBin->FindBin(track.pt());
          if (pTBinForEff >= 1 && pTBinForEff <= hFindPtBin->GetNbinsX()) {
            wEPeff = funcEff[pTBinForEff - 1]->Eval(fphi * tracks.size());
            if (wEPeff > 0.) {
              wEPeff = 1. / wEPeff;
              weff *= wEPeff;
            }
          }
        }
      } // cfgDoLocDenCorr

      // track cut
      if (track.itsNCls() <= trkQualityOpts.cfgITSNCls.value)
        continue;
      if (track.tpcNClsCrossedRows() <= trkQualityOpts.cfgTPCCrossedRows.value)
        continue;
      if (track.tpcNClsFound() <= trkQualityOpts.cfgTPCNCls.value)
        continue;
      // end track cut

      registry.fill(HIST("hPhi"), track.phi());
      registry.fill(HIST("hPhicorr"), track.phi(), wacc);
      registry.fill(HIST("hEta"), track.eta());
      registry.fill(HIST("hEtaPhiVtxzREF"), track.phi(), track.eta(), vtxz, wacc);
      registry.fill(HIST("hPt"), track.pt());
      if (!((track.pt() > trkQualityOpts.cfgCutPtMin.value) && (track.pt() < trkQualityOpts.cfgCutPtMax.value)))
        continue;
      fGFW->Fill(track.eta(), 0, track.phi(), wacc * weff, 1); //(eta, ptbin, phi, wacc*weff, bitmask)

      if (this->isPion(track)) {
        fGFW->Fill(track.eta(), 0, track.phi(), wacc * weff, 18);
      }

      if (this->isKaon(track)) {
        fGFW->Fill(track.eta(), 0, track.phi(), wacc * weff, 36);
      }

      if (this->isProton(track)) {
        fGFW->Fill(track.eta(), 0, track.phi(), wacc * weff, 72);
      }

      if (cfgOutputNUAWeights)
        fWeightsREF->fill(track.phi(), track.eta(), vtxz, track.pt(), cent, 0);

      if (std::fabs(track.eta()) < trkQualityOpts.cfgRangeEta.value) {
        nch += weff;
        nchSquare += weff * weff;
        ptSum += weff * track.pt();
        ptSumw2 += weff * weff * track.pt();
        ptSquareSum += weff * weff * track.pt() * track.pt();
      }
    }

    if (nch > 0) {
      int centbin = 0;
      centbin = fMultAxis->FindBin(cent) - 1;

      fillProfile(corrconfigs.at(0), HIST("c22"), cent);
      fillProfile(corrconfigs.at(1), HIST("c24"), cent);
      fillProfile(corrconfigs.at(2), HIST("c22Full"), cent);
      fillProfile(corrconfigs.at(3), HIST("c32"), cent);
      fillProfile(corrconfigs.at(4), HIST("c34"), cent);

      fillProfile(corrconfigs.at(5), HIST("pi/c22"), cent);
      fillProfile(corrconfigs.at(6), HIST("pi/c22"), cent);
      fillProfile(corrconfigs.at(7), HIST("ka/c22"), cent);
      fillProfile(corrconfigs.at(8), HIST("ka/c22"), cent);
      fillProfile(corrconfigs.at(9), HIST("pr/c22"), cent);
      fillProfile(corrconfigs.at(10), HIST("pr/c22"), cent);

      fillProfile(corrconfigs.at(11), HIST("pi/c24"), cent);
      fillProfile(corrconfigs.at(12), HIST("pi/c24"), cent);
      fillProfile(corrconfigs.at(13), HIST("ka/c24"), cent);
      fillProfile(corrconfigs.at(14), HIST("ka/c24"), cent);
      fillProfile(corrconfigs.at(15), HIST("pr/c24"), cent);
      fillProfile(corrconfigs.at(16), HIST("pr/c24"), cent);

      fillProfile(corrconfigs.at(17), HIST("pi/c32"), cent);
      fillProfile(corrconfigs.at(18), HIST("pi/c32"), cent);
      fillProfile(corrconfigs.at(19), HIST("ka/c32"), cent);
      fillProfile(corrconfigs.at(20), HIST("ka/c32"), cent);
      fillProfile(corrconfigs.at(21), HIST("pr/c32"), cent);
      fillProfile(corrconfigs.at(22), HIST("pr/c32"), cent);

      fillProfile(corrconfigs.at(23), HIST("pi/c34"), cent);
      fillProfile(corrconfigs.at(24), HIST("pi/c34"), cent);
      fillProfile(corrconfigs.at(25), HIST("ka/c34"), cent);
      fillProfile(corrconfigs.at(26), HIST("ka/c34"), cent);
      fillProfile(corrconfigs.at(27), HIST("pr/c34"), cent);
      fillProfile(corrconfigs.at(28), HIST("pr/c34"), cent);

      fillProfile(corrconfigs.at(29), HIST("pi/c22PiPi"), cent);
      fillProfile(corrconfigs.at(30), HIST("ka/c22KaKa"), cent);
      fillProfile(corrconfigs.at(31), HIST("pr/c22PrPr"), cent);
      fillProfile(corrconfigs.at(32), HIST("pi/c32PiPi"), cent);
      fillProfile(corrconfigs.at(33), HIST("ka/c32KaKa"), cent);
      fillProfile(corrconfigs.at(34), HIST("pr/c32PrPr"), cent);

      fillProfilevnpt(corrconfigs.at(0), HIST("c22TrackWeight"), cent, nch, nch, 0);
      fillProfilevnpt(corrconfigs.at(1), HIST("c24TrackWeight"), cent, nch, nch, 0);
      fillProfilevnpt(corrconfigs.at(2), HIST("c22FullTrackWeight"), cent, nch, nch, 0);
      fillProfilevnpt(corrconfigs.at(3), HIST("c32TrackWeight"), cent, nch, nch, 0);
      fillProfilevnpt(corrconfigs.at(4), HIST("c34TrackWeight"), cent, nch, nch, 0);

      fillProfilevnpt(corrconfigs.at(29), HIST("pi/c22TrackWeight"), cent, nch, nch, 0);
      fillProfilevnpt(corrconfigs.at(30), HIST("ka/c22TrackWeight"), cent, nch, nch, 0);
      fillProfilevnpt(corrconfigs.at(31), HIST("pr/c22TrackWeight"), cent, nch, nch, 0);

      fillProfilevnpt(corrconfigs.at(32), HIST("pi/c32TrackWeight"), cent, nch, nch, 0);
      fillProfilevnpt(corrconfigs.at(33), HIST("ka/c32TrackWeight"), cent, nch, nch, 0);
      fillProfilevnpt(corrconfigs.at(34), HIST("pr/c32TrackWeight"), cent, nch, nch, 0);

      fillProfilevnpt(corrconfigs.at(0), HIST("covV2Pt"), cent, ptSum, nch, 0);
      fillProfilevnpt(corrconfigs.at(0), HIST("covV2Pt_diffpt"), cent, ptSum, nch, cfgMeanPt[centbin]);
      fillProfilevnpt(corrconfigs.at(29), HIST("pi/covV2Pt"), cent, ptSum, nch, 0);
      fillProfilevnpt(corrconfigs.at(29), HIST("pi/covV2Pt_diffpt"), cent, ptSum, nch, cfgMeanPt[centbin]);
      fillProfilevnpt(corrconfigs.at(30), HIST("ka/covV2Pt"), cent, ptSum, nch, 0);
      fillProfilevnpt(corrconfigs.at(30), HIST("ka/covV2Pt_diffpt"), cent, ptSum, nch, cfgMeanPt[centbin]);
      fillProfilevnpt(corrconfigs.at(31), HIST("pr/covV2Pt"), cent, ptSum, nch, 0);
      fillProfilevnpt(corrconfigs.at(31), HIST("pr/covV2Pt_diffpt"), cent, ptSum, nch, cfgMeanPt[centbin]);

      fillProfilevnpt(corrconfigs.at(3), HIST("covV3Pt"), cent, ptSum, nch, 0);
      fillProfilevnpt(corrconfigs.at(3), HIST("covV3Pt_diffpt"), cent, ptSum, nch, cfgMeanPt[centbin]);
      fillProfilevnpt(corrconfigs.at(32), HIST("pi/covV3Pt"), cent, ptSum, nch, 0);
      fillProfilevnpt(corrconfigs.at(32), HIST("pi/covV3Pt_diffpt"), cent, ptSum, nch, cfgMeanPt[centbin]);
      fillProfilevnpt(corrconfigs.at(33), HIST("ka/covV3Pt"), cent, ptSum, nch, 0);
      fillProfilevnpt(corrconfigs.at(33), HIST("ka/covV3Pt_diffpt"), cent, ptSum, nch, cfgMeanPt[centbin]);
      fillProfilevnpt(corrconfigs.at(34), HIST("pr/covV3Pt"), cent, ptSum, nch, 0);
      fillProfilevnpt(corrconfigs.at(34), HIST("pr/covV3Pt_diffpt"), cent, ptSum, nch, cfgMeanPt[centbin]);

      fillProfilePOIvnpt(corrconfigs.at(0), HIST("c22dmeanpt"), cent, ptSum, nch);
      fillProfilePOIvnpt(corrconfigs.at(5), HIST("pi/c22dmeanpt"), cent, ptSum, nch);
      fillProfilePOIvnpt(corrconfigs.at(6), HIST("pi/c22dmeanpt"), cent, ptSum, nch);
      fillProfilePOIvnpt(corrconfigs.at(7), HIST("ka/c22dmeanpt"), cent, ptSum, nch);
      fillProfilePOIvnpt(corrconfigs.at(8), HIST("ka/c22dmeanpt"), cent, ptSum, nch);
      fillProfilePOIvnpt(corrconfigs.at(9), HIST("pr/c22dmeanpt"), cent, ptSum, nch);
      fillProfilePOIvnpt(corrconfigs.at(10), HIST("pr/c22dmeanpt"), cent, ptSum, nch);
      registry.fill(HIST("hMeanPt"), cent, (ptSum / nch), nch);

      double nchDiff = nch * nch - nchSquare;
      if (nchDiff) {
        registry.fill(HIST("ptSquareAve"), cent,
                      (ptSum * ptSum - ptSquareSum) / nchDiff,
                      nchDiff);
        registry.fill(HIST("ptAve"), cent,
                      (nch * ptSum - ptSumw2) / nchDiff,
                      nchDiff);
      }
    }
  }
  PROCESS_SWITCH(PidFlowPtCorr, processData, "", true);

  /**
   * @brief this function is used to fill THn hist for NUA correction and for NUE correction
   * @details hist THn: (runNumberIDX, phi, eta, Vz), note that different runNumber will be put in the same hist
   *
   * @param collision
   * @param tracks
   */
  void fillCorrectionGraph(AodCollisions::iterator const& collision, aod::BCsWithTimestamps const&, AodTracks const& tracks)
  {
    // init
    int nTot = tracks.size();
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int runNumber = bc.runNumber();
    double interactionRate = rateFetcher.fetch(ccdb.service, bc.timestamp(), runNumber, "ZNC hadronic") * 1.e-3;
    // end init

    // collision cut
    // include : 1.track.size 2.collision.sel8 3.this->evenSelected
    if (nTot < 1)
      return;
    fGFW->Clear();
    const auto cent = collision.centFT0C();
    if (!collision.sel8())
      return;

    if (!eventSelected(collision, cent, interactionRate))
      return;
    // end collision cut

    // loop the vector, find the place to put (phi eta Vz)
    // if the run number is new, create one
    int matchedPosition = -1;
    for (uint64_t idxPosition = 0; idxPosition < this->runNumbers.size(); idxPosition++) {
      if (this->runNumbers[idxPosition] == runNumber) {
        matchedPosition = idxPosition;
        break;
      }
    }
    if (matchedPosition == -1) {
      return;
    }
    // end find place to put run data

    // loop all the track
    for (const auto& track : tracks) {
      // track cut
      if (!((track.pt() > trkQualityOpts.cfgCutPtMin.value) && (track.pt() < trkQualityOpts.cfgCutPtMax.value)))
        continue;
      if (track.itsNCls() <= trkQualityOpts.cfgITSNCls.value)
        continue;
      if (track.tpcNClsCrossedRows() <= trkQualityOpts.cfgTPCCrossedRows.value)
        continue;
      if (track.tpcNClsFound() <= trkQualityOpts.cfgTPCNCls.value)
        continue;
      // end track cut

      // fill the THn
      registry.fill(HIST("correction/hRunNumberPhiEtaVertex"), matchedPosition, track.phi(), track.eta(), collision.posZ());
      // end fill the THn

    } // end loop all the track
  }
  PROCESS_SWITCH(PidFlowPtCorr, fillCorrectionGraph, "", true);

  /**
   * @brief this main function is used to check the PID performance of ITS TOC TPC
   *
   * @param collision
   * @param tracks
   */
  void detectorPidQA(AodCollisions::iterator const& collision, aod::BCsWithTimestamps const&, AodTracks const& tracks)
  {
    // cut and correction
    o2::aod::ITSResponse itsResponse;
    int nTot = tracks.size();
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int runNumber = bc.runNumber();
    double interactionRate = rateFetcher.fetch(ccdb.service, bc.timestamp(), runNumber, "ZNC hadronic") * 1.e-3;

    if (nTot < 1)
      return;
    fGFW->Clear();
    const auto cent = collision.centFT0C();
    if (!collision.sel8())
      return;

    if (!eventSelected(collision, cent, interactionRate))
      return;
    // loadCorrections(bc.timestamp());
    // end cut and correction

    // start filling graphs
    for (const auto& track : tracks) {
      // track cut
      if (!((track.pt() > trkQualityOpts.cfgCutPtMin.value) && (track.pt() < trkQualityOpts.cfgCutPtMax.value)))
        continue;
      if (track.itsNCls() <= trkQualityOpts.cfgITSNCls.value)
        continue;
      if (track.tpcNClsCrossedRows() <= trkQualityOpts.cfgTPCCrossedRows.value)
        continue;
      if (track.tpcNClsFound() <= trkQualityOpts.cfgTPCNCls.value)
        continue;
      // end track cut

      // TPC TOF
      registry.fill(HIST("DetectorPidPerformace/TPCvsTOF/Pi"), track.tpcNSigmaPi(), track.tofNSigmaPi(), track.pt());
      registry.fill(HIST("DetectorPidPerformace/TPCvsTOF/Pr"), track.tpcNSigmaPr(), track.tofNSigmaPr(), track.pt());
      registry.fill(HIST("DetectorPidPerformace/TPCvsTOF/Ka"), track.tpcNSigmaKa(), track.tofNSigmaKa(), track.pt());

      // TPC ITS
      registry.fill(HIST("DetectorPidPerformace/TPCvsITS/Pi"), track.tpcNSigmaPi(), itsResponse.nSigmaITS<o2::track::PID::Pion>(track), track.pt());
      registry.fill(HIST("DetectorPidPerformace/TPCvsITS/Pr"), track.tpcNSigmaPr(), itsResponse.nSigmaITS<o2::track::PID::Proton>(track), track.pt());
      registry.fill(HIST("DetectorPidPerformace/TPCvsITS/Ka"), track.tpcNSigmaKa(), itsResponse.nSigmaITS<o2::track::PID::Kaon>(track), track.pt());

      // ITS vs TOF
      registry.fill(HIST("DetectorPidPerformace/ITSvsTOF/Pi"), itsResponse.nSigmaITS<o2::track::PID::Pion>(track), track.tofNSigmaPi(), track.pt());
      registry.fill(HIST("DetectorPidPerformace/ITSvsTOF/Pr"), itsResponse.nSigmaITS<o2::track::PID::Proton>(track), track.tofNSigmaPr(), track.pt());
      registry.fill(HIST("DetectorPidPerformace/ITSvsTOF/Ka"), itsResponse.nSigmaITS<o2::track::PID::Kaon>(track), track.tofNSigmaKa(), track.pt());

    } // end filling graphs
  }
  PROCESS_SWITCH(PidFlowPtCorr, detectorPidQA, "", true);

  /**
   * @brief this function is used to loop mc events
   * @note implemented:
   *       1.fill mc pt hist to calculate efficiency(correction/hCentPtMC, correction/hCentPtData)
   *
   */
  void processMCGen(FilteredMcCollisions::iterator const&,
                    aod::BCsWithTimestamps const&,
                    FilteredMcParticles const& mcParticles,
                    soa::SmallGroups<soa::Join<aod::McCollisionLabels, AodCollisions>> const& collisions,
                    AodTracks const&)
  {
    // cut && init
    // loop all the collisions reco matched to MC
    double cent = -1;
    for (const auto& oneColl : collisions) {
      auto bc = oneColl.bc_as<aod::BCsWithTimestamps>();
      int runNumber = bc.runNumber();
      double interactionRate = rateFetcher.fetch(ccdb.service, bc.timestamp(), runNumber, "ZNC hadronic") * 1.e-3;
      cent = oneColl.centFT0C();

      // collision cut
      if (!oneColl.sel8())
        return;
      if (!eventSelected(oneColl, cent, interactionRate))
        return;
      // end collision cut
    }
    // end loop all the collisions reco matched to MC
    // end cut && init

    // loop all the particle
    for (const auto& mcPart : mcParticles) {
      // track cut
      if (!mcPart.isPhysicalPrimary())
        continue;
      // end track cut

      registry.fill(HIST("correction/hCentPtMC"), cent, mcPart.pt());

      // loop real track
      if (mcPart.has_tracks()) {
        auto const& tracks = mcPart.tracks_as<AodTracks>();
        for (const auto& track : tracks) {
          // track select
          // end track select
          registry.fill(HIST("correction/hCentPtData"), cent, track.pt());
        }
      }
      // end loop real track
    }
    // end loop all the particle
  }
  PROCESS_SWITCH(PidFlowPtCorr, processMCGen, "", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<PidFlowPtCorr>(cfgc)};
}
