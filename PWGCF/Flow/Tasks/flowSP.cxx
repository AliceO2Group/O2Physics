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

/// \file   flowSP.cxx
/// \author Noor Koster
/// \since  01/12/2024
/// \brief  task to evaluate flow with respect to spectator plane.

#include <CCDB/BasicCCDBManager.h>
#include <DataFormatsParameters/GRPObject.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <algorithm>
#include <numeric>
#include <vector>
#include <string>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/HistogramRegistry.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/Core/RecoDecay.h"

#include "PWGCF/DataModel/SPTableZDC.h"
#include "GFWWeights.h"
#include "TF1.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
// using namespace o2::analysis;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct FlowSP {

  O2_DEFINE_CONFIGURABLE(cfgDCAxy, float, 0.2, "Cut on DCA in the transverse direction (cm)");
  O2_DEFINE_CONFIGURABLE(cfgDCAz, float, 2, "Cut on DCA in the longitudinal direction (cm)");
  O2_DEFINE_CONFIGURABLE(cfgNcls, float, 70, "Cut on number of TPC clusters found");
  O2_DEFINE_CONFIGURABLE(cfgPtmin, float, 0.2, "minimum pt (GeV/c)");
  O2_DEFINE_CONFIGURABLE(cfgPtmax, float, 10, "maximum pt (GeV/c)");
  O2_DEFINE_CONFIGURABLE(cfgEta, float, 0.8, "eta cut");
  O2_DEFINE_CONFIGURABLE(cfgVtxZ, float, 10, "vertex cut (cm)");
  O2_DEFINE_CONFIGURABLE(cfgMagField, float, 99999, "Configurable magnetic field;default CCDB will be queried");
  O2_DEFINE_CONFIGURABLE(cfgUseAdditionalEventCut, bool, true, "Bool to enable Additional Event Cut");
  O2_DEFINE_CONFIGURABLE(cfgUseAdditionalTrackCut, bool, true, "Bool to enable Additional Track Cut");
  O2_DEFINE_CONFIGURABLE(cfgCentMin, float, 0, "Minimum cenrality for selected events");
  O2_DEFINE_CONFIGURABLE(cfgCentMax, float, 90, "Maximum cenrality for selected events");
  O2_DEFINE_CONFIGURABLE(cfgFillWeights, bool, true, "Fill NUA weights");
  O2_DEFINE_CONFIGURABLE(cfgFillWeightsPOS, bool, false, "Fill NUA weights only for positive charges");
  O2_DEFINE_CONFIGURABLE(cfgFillWeightsNEG, bool, false, "Fill NUA weights only for negative charges");
  O2_DEFINE_CONFIGURABLE(cfgAcceptance, std::string, "", "ccdb dir for NUA corrections");
  O2_DEFINE_CONFIGURABLE(cfgEfficiency, std::string, "", "ccdb dir for NUE corrections");
  O2_DEFINE_CONFIGURABLE(cfgDoubleTrackFunction, bool, true, "Include track cut at low pt");
  O2_DEFINE_CONFIGURABLE(cfgTrackCutSize, float, 0.06, "Spread of track cut");
  O2_DEFINE_CONFIGURABLE(cfgMaxOccupancy, int, 10000, "Maximum occupancy of selected events");
  O2_DEFINE_CONFIGURABLE(cfgNoSameBunchPileupCut, bool, true, "kNoSameBunchPileupCut");
  O2_DEFINE_CONFIGURABLE(cfgIsGoodZvtxFT0vsPV, bool, true, "kIsGoodZvtxFT0vsPV");
  O2_DEFINE_CONFIGURABLE(cfgNoCollInTimeRangeStandard, bool, true, "kNoCollInTimeRangeStandard");
  O2_DEFINE_CONFIGURABLE(cfgDoOccupancySel, bool, true, "Bool for event selection on detector occupancy");
  O2_DEFINE_CONFIGURABLE(cfgMultCut, bool, true, "Use additional evenr cut on mult correlations");
  O2_DEFINE_CONFIGURABLE(cfgTVXinTRD, bool, false, "Use kTVXinTRD (reject TRD triggered events)");
  O2_DEFINE_CONFIGURABLE(cfgIsVertexITSTPC, bool, true, "Selects collisions with at least one ITS-TPC track");
  O2_DEFINE_CONFIGURABLE(cfgIsGoodITSLayersAll, bool, true, "Cut time intervals with dead ITS staves");
  O2_DEFINE_CONFIGURABLE(cfgCCDBdir, std::string, "Users/c/ckoster/ZDC/LHC23_zzh_pass4_small/meanQQ", "ccdb dir for average QQ values in 1% centrality bins");
  O2_DEFINE_CONFIGURABLE(cfgLoadAverageQQ, bool, true, "Load average values for QQ (in centrality bins)");
  O2_DEFINE_CONFIGURABLE(cfgHarm, int, 1, "Flow harmonic n for ux and uy: (Cos(n*phi), Sin(n*phi))");
  O2_DEFINE_CONFIGURABLE(cfgHarmMixed, int, 2, "Flow harmonic n for ux and uy in mixed harmonics (MH): (Cos(n*phi), Sin(n*phi))");
  O2_DEFINE_CONFIGURABLE(cfgLoadSPPlaneRes, bool, false, "Load ZDC spectator plane resolution");
  O2_DEFINE_CONFIGURABLE(cfgCCDBdir_SP, std::string, "Users/c/ckoster/ZDC/LHC23_zzh_pass4_small/SPPlaneRes", "ccdb dir for average event plane resolution in 1% centrality bins");

  ConfigurableAxis axisDCAz{"axisDCAz", {200, -.5, .5}, "DCA_{z} (cm)"};
  ConfigurableAxis axisDCAxy{"axisDCAxy", {200, -.5, .5}, "DCA_{xy} (cm)"};
  ConfigurableAxis axisPhiMod = {"axisPhiMod", {100, 0, constants::math::PI / 9}, "fmod(#varphi,#pi/9)"};
  ConfigurableAxis axisPhi = {"axisPhi", {60, 0, constants::math::TwoPI}, "#varphi"};
  ConfigurableAxis axisEta = {"axisEta", {64, -1.8, 1.8}, "#eta"};
  ConfigurableAxis axisEtaVn = {"axisEtaVn", {8, -.8, .8}, "#eta"};
  ConfigurableAxis axisVx = {"axisVx", {40, -0.01, 0.01}, "v_{x}"};
  ConfigurableAxis axisVy = {"axisVy", {40, -0.01, 0.01}, "v_{y}"};
  ConfigurableAxis axisVz = {"axisVz", {40, -10, 10}, "v_{z}"};
  ConfigurableAxis axisCent = {"axisCent", {90, 0, 90}, "Centrality(%)"};
  ConfigurableAxis axisPhiPlane = {"axisPhiPlane", {100, -constants::math::PI, constants::math::PI}, "#Psi"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgVtxZ;
  Filter trackFilter = nabs(aod::track::eta) < cfgEta && aod::track::pt > cfgPtmin&& aod::track::pt < cfgPtmax && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && nabs(aod::track::dcaXY) < cfgDCAxy&& nabs(aod::track::dcaZ) < cfgDCAz;
  using UsedCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs, aod::SPTableZDC>>;
  using UsedTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>>;

  //  Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;

  // from Generic Framework
  struct Config {
    std::vector<TH1D*> mEfficiency = {};
    std::vector<GFWWeights*> mAcceptance = {};
    bool correctionsLoaded = false;
    int lastRunNumber = 0;
  } cfg;

  OutputObj<GFWWeights> fWeights{GFWWeights("weights")};
  OutputObj<GFWWeights> fWeightsPOS{GFWWeights("weights_positive")};
  OutputObj<GFWWeights> fWeightsNEG{GFWWeights("weights_negative")};
  HistogramRegistry registry{"registry"};

  // Event selection cuts - Alex
  TF1* fPhiCutLow = nullptr;
  TF1* fPhiCutHigh = nullptr;
  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;
  TF1* fMultCutLow = nullptr;
  TF1* fMultCutHigh = nullptr;
  TF1* fMultMultPVCut = nullptr;

  enum SelectionCriteria {
    evSel_FilteredEvent,
    evSel_sel8,
    evSel_occupancy,
    evSel_kTVXinTRD,
    evSel_kNoSameBunchPileup,
    evSel_kIsGoodZvtxFT0vsPV,
    evSel_kNoCollInTimeRangeStandard,
    evSel_kIsVertexITSTPC,
    evSel_MultCuts,
    evSel_CentCuts,
    evSel_kIsGoodITSLayersAll,
    evSel_isSelectedZDC,
    nEventSelections
  };

  enum ChargeType {
    kInclusive,
    kPositive,
    kNegative
  };

  void init(InitContext const&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    std::vector<double> ptbinning = {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.5, 4, 5, 6, 8, 10};
    AxisSpec axisPt = {ptbinning, "#it{p}_{T} GeV/#it{c}"};
    AxisSpec nchAxis = {4000, 0, 4000, "N_{ch}"};
    AxisSpec t0cAxis = {70, 0, 70000, "N_{ch} (T0C)"};
    AxisSpec t0aAxis = {200, 0, 200, "N_{ch}"};
    AxisSpec multpvAxis = {4000, 0, 4000, "N_{ch} (PV)"};

    int ptbins = ptbinning.size() - 1;

    if (cfgFillWeights) {
      fWeights->SetPtBins(ptbins, &ptbinning[0]);
      fWeights->Init(true, false);

      fWeightsPOS->SetPtBins(ptbins, &ptbinning[0]);
      fWeightsPOS->Init(true, false);

      fWeightsNEG->SetPtBins(ptbins, &ptbinning[0]);
      fWeightsNEG->Init(true, false);
    }

    registry.add<TH1>("hSPplaneA", "hSPplaneA", kTH1D, {axisPhiPlane});
    registry.add<TH1>("hSPplaneC", "hSPplaneC", kTH1D, {axisPhiPlane});
    registry.add<TH1>("hSPplaneFull", "hSPplaneFull", kTH1D, {axisPhiPlane});

    registry.add<TProfile>("hCosPhiACosPhiC", "hCosPhiACosPhiC; Centrality(%); #LT Cos(#Psi^{A})Cos(#Psi^{C})#GT", kTProfile, {axisCent});
    registry.add<TProfile>("hSinPhiASinPhiC", "hSinPhiASinPhiC; Centrality(%); #LT Sin(#Psi^{A})Sin(#Psi^{C})#GT", kTProfile, {axisCent});
    registry.add<TProfile>("hSinPhiACosPhiC", "hSinPhiACosPhiC; Centrality(%); #LT Sin(#Psi^{A})Cos(#Psi^{C})#GT", kTProfile, {axisCent});
    registry.add<TProfile>("hCosPhiASinsPhiC", "hCosPhiASinsPhiC; Centrality(%); #LT Cos(#Psi^{A})Sin(#Psi^{C})#GT", kTProfile, {axisCent});
    registry.add<TProfile>("hFullEvPlaneRes", "hFullEvPlaneRes; Centrality(%); -#LT Cos(#Psi^{A} - #Psi^{C})#GT ", kTProfile, {axisCent});

    registry.add("QA/after/hCent", "", {HistType::kTH1D, {axisCent}});
    registry.add("QA/after/pt_phi", "", {HistType::kTH2D, {axisPt, axisPhiMod}});
    registry.add("QA/after/hPt_inclusive", "", {HistType::kTH1D, {axisPt}});
    registry.add("QA/after/globalTracks_centT0C", "", {HistType::kTH2D, {axisCent, nchAxis}});
    registry.add("QA/after/PVTracks_centT0C", "", {HistType::kTH2D, {axisCent, multpvAxis}});
    registry.add("QA/after/globalTracks_PVTracks", "", {HistType::kTH2D, {multpvAxis, nchAxis}});
    registry.add("QA/after/globalTracks_multT0A", "", {HistType::kTH2D, {t0aAxis, nchAxis}});
    registry.add("QA/after/globalTracks_multV0A", "", {HistType::kTH2D, {t0aAxis, nchAxis}});
    registry.add("QA/after/multV0A_multT0A", "", {HistType::kTH2D, {t0aAxis, t0aAxis}});
    registry.add("QA/after/multT0C_centT0C", "", {HistType::kTH2D, {axisCent, t0cAxis}});

    registry.add("QA/after/PsiA_vs_Cent", "", {HistType::kTH2D, {axisPhiPlane, axisCent}});
    registry.add("QA/after/PsiC_vs_Cent", "", {HistType::kTH2D, {axisPhiPlane, axisCent}});
    registry.add("QA/after/PsiFull_vs_Cent", "", {HistType::kTH2D, {axisPhiPlane, axisCent}});

    registry.add("QA/after/PsiA_vs_Vx", "", {HistType::kTH2D, {axisPhiPlane, axisVx}});
    registry.add("QA/after/PsiC_vs_Vx", "", {HistType::kTH2D, {axisPhiPlane, axisVx}});
    registry.add("QA/after/PsiFull_vs_Vx", "", {HistType::kTH2D, {axisPhiPlane, axisVx}});

    registry.add("QA/after/PsiA_vs_Vy", "", {HistType::kTH2D, {axisPhiPlane, axisVy}});
    registry.add("QA/after/PsiC_vs_Vy", "", {HistType::kTH2D, {axisPhiPlane, axisVy}});
    registry.add("QA/after/PsiFull_vs_Vy", "", {HistType::kTH2D, {axisPhiPlane, axisVy}});

    registry.add("QA/after/PsiA_vs_Vz", "", {HistType::kTH2D, {axisPhiPlane, axisVz}});
    registry.add("QA/after/PsiC_vs_Vz", "", {HistType::kTH2D, {axisPhiPlane, axisVz}});
    registry.add("QA/after/PsiFull_vs_Vz", "", {HistType::kTH2D, {axisPhiPlane, axisVz}});

    registry.addClone("QA/after/", "QA/before/");

    // track properties per centrality and per eta, pt bin
    registry.add<TProfile>("incl/vnAx_eta", "", kTProfile, {axisEtaVn});
    registry.add<TProfile>("incl/vnAy_eta", "", kTProfile, {axisEtaVn});
    registry.add<TProfile>("incl/vnCx_eta", "", kTProfile, {axisEtaVn});
    registry.add<TProfile>("incl/vnCy_eta", "", kTProfile, {axisEtaVn});
    registry.add<TProfile>("incl/vnC_eta", "", kTProfile, {axisEtaVn});
    registry.add<TProfile>("incl/vnA_eta", "", kTProfile, {axisEtaVn});
    registry.add<TProfile>("incl/vnA_eta_EP", "", kTProfile, {axisEtaVn});
    registry.add<TProfile>("incl/vnC_eta_EP", "", kTProfile, {axisEtaVn});
    registry.add<TProfile>("incl/vnFull_eta_EP", "", kTProfile, {axisEtaVn});
    registry.add<TProfile>("incl/vnAxCxUx_eta_MH", "", kTProfile, {axisEtaVn});
    registry.add<TProfile>("incl/vnAxCyUx_eta_MH", "", kTProfile, {axisEtaVn});
    registry.add<TProfile>("incl/vnAxCyUy_eta_MH", "", kTProfile, {axisEtaVn});
    registry.add<TProfile>("incl/vnAyCxUy_eta_MH", "", kTProfile, {axisEtaVn});

    registry.add<TProfile>("incl/vnAx_pt", "", kTProfile, {axisPt});
    registry.add<TProfile>("incl/vnAy_pt", "", kTProfile, {axisPt});
    registry.add<TProfile>("incl/vnCx_pt", "", kTProfile, {axisPt});
    registry.add<TProfile>("incl/vnCy_pt", "", kTProfile, {axisPt});
    registry.add<TProfile>("incl/vnC_pt", "", kTProfile, {axisPt});
    registry.add<TProfile>("incl/vnA_pt", "", kTProfile, {axisPt});
    registry.add<TProfile>("incl/vnC_pt_odd", "", kTProfile, {axisPt});
    registry.add<TProfile>("incl/vnA_pt_odd", "", kTProfile, {axisPt});
    registry.add<TProfile>("incl/vnA_pt_EP", "", kTProfile, {axisPt});
    registry.add<TProfile>("incl/vnC_pt_EP", "", kTProfile, {axisPt});
    registry.add<TProfile>("incl/vnFull_pt_EP", "", kTProfile, {axisPt});
    registry.add<TProfile>("incl/vnAxCxUx_pt_MH", "", kTProfile, {axisPt});
    registry.add<TProfile>("incl/vnAxCyUx_pt_MH", "", kTProfile, {axisPt});
    registry.add<TProfile>("incl/vnAxCyUy_pt_MH", "", kTProfile, {axisPt});
    registry.add<TProfile>("incl/vnAyCxUy_pt_MH", "", kTProfile, {axisPt});

    registry.add<TProfile>("incl/vnC_cent_minEta", "", kTProfile, {axisCent});
    registry.add<TProfile>("incl/vnA_cent_minEta", "", kTProfile, {axisCent});
    registry.add<TProfile>("incl/vnC_cent_plusEta", "", kTProfile, {axisCent});
    registry.add<TProfile>("incl/vnA_cent_plusEta", "", kTProfile, {axisCent});
    registry.add<TProfile>("incl/vnA_cent_EP", "", kTProfile, {axisCent});
    registry.add<TProfile>("incl/vnC_cent_EP", "", kTProfile, {axisCent});
    registry.add<TProfile>("incl/vnFull_cent_EP", "", kTProfile, {axisCent});
    registry.add<TProfile>("incl/vnAxCxUx_cent_MH", "", kTProfile, {axisCent});
    registry.add<TProfile>("incl/vnAxCyUx_cent_MH", "", kTProfile, {axisCent});
    registry.add<TProfile>("incl/vnAxCyUy_cent_MH", "", kTProfile, {axisCent});
    registry.add<TProfile>("incl/vnAyCxUy_cent_MH", "", kTProfile, {axisCent});

    // track QA for pos, neg, incl
    registry.add<TH1>("incl/QA/hPt", "", kTH1D, {axisPt});
    registry.add<TH1>("incl/QA/hPhi", "", kTH1D, {axisPhi});
    registry.add<TH1>("incl/QA/hEta", "", kTH1D, {axisEta});
    registry.add<TH3>("incl/QA/hPhi_Eta_vz", "", kTH3D, {axisPhi, axisEta, axisVz});
    registry.add<TH1>("incl/QA/hDCAxy", "", kTH1D, {axisDCAxy});
    registry.add<TH1>("incl/QA/hDCAz", "", kTH1D, {axisDCAz});

    registry.addClone("incl/", "pos/");
    registry.addClone("incl/", "neg/");

    registry.add<TProfile>("qAqCX", "", kTProfile, {axisCent});
    registry.add<TProfile>("qAqCY", "", kTProfile, {axisCent});
    registry.add<TProfile>("qAqCXY", "", kTProfile, {axisCent});

    registry.add("hEventCount", "Number of Event; Cut; #Events Passed Cut", {HistType::kTH1D, {{nEventSelections, 0, nEventSelections}}});
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_FilteredEvent + 1, "Filtered event");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_sel8 + 1, "Sel8");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_occupancy + 1, "kOccupancy");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_kTVXinTRD + 1, "kTVXinTRD");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_kNoSameBunchPileup + 1, "kNoSameBunchPileup");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_kIsGoodZvtxFT0vsPV + 1, "kIsGoodZvtxFT0vsPV");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_kNoCollInTimeRangeStandard + 1, "kNoCollInTimeRangeStandard");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_kIsVertexITSTPC + 1, "kIsVertexITSTPC");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_MultCuts + 1, "Mult cuts (Alex)");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_CentCuts + 1, "Cenrality range");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_kIsGoodITSLayersAll + 1, "kkIsGoodITSLayersAll");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_isSelectedZDC + 1, "isSelected");

    if (cfgUseAdditionalEventCut) {
      fMultPVCutLow = new TF1("fMultPVCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x - 3.5*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)", 0, 100);
      fMultPVCutLow->SetParameters(3257.29, -121.848, 1.98492, -0.0172128, 6.47528e-05, 154.756, -1.86072, -0.0274713, 0.000633499, -3.37757e-06);
      fMultPVCutHigh = new TF1("fMultPVCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x + 3.5*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)", 0, 100);
      fMultPVCutHigh->SetParameters(3257.29, -121.848, 1.98492, -0.0172128, 6.47528e-05, 154.756, -1.86072, -0.0274713, 0.000633499, -3.37757e-06);

      fMultCutLow = new TF1("fMultCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 2.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultCutLow->SetParameters(1654.46, -47.2379, 0.449833, -0.0014125, 150.773, -3.67334, 0.0530503, -0.000614061, 3.15956e-06);
      fMultCutHigh = new TF1("fMultCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 3.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultCutHigh->SetParameters(1654.46, -47.2379, 0.449833, -0.0014125, 150.773, -3.67334, 0.0530503, -0.000614061, 3.15956e-06);
    }

    if (cfgUseAdditionalTrackCut) {
      fPhiCutLow = new TF1("fPhiCutLow", "0.06/x+pi/18.0-0.06", 0, 100);
      fPhiCutHigh = new TF1("fPhiCutHigh", "0.1/x+pi/18.0+0.06", 0, 100);
    }
  }

  int getMagneticField(uint64_t timestamp)
  {
    // TODO done only once (and not per run). Will be replaced by CCDBConfigurable
    // static o2::parameters::GRPObject* grpo = nullptr;
    static o2::parameters::GRPMagField* grpo = nullptr;
    if (grpo == nullptr) {
      // grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>("GLO/GRP/GRP", timestamp);
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 0;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
    }
    return grpo->getNominalL3Field();
  }

  // From Generic Framework
  void loadCorrections(uint64_t timestamp)
  {
    // corrections saved on CCDB as TList {incl, pos, neg} of GFWWeights (acc) TH1D (eff) objects!
    if (cfg.correctionsLoaded)
      return;

    if (cfgAcceptance.value.empty() == false) {
      TList* listCorrections = ccdb->getForTimeStamp<TList>(cfgAcceptance, timestamp);
      cfg.mAcceptance.push_back(reinterpret_cast<GFWWeights*>(listCorrections->FindObject("weights")));
      cfg.mAcceptance.push_back(reinterpret_cast<GFWWeights*>(listCorrections->FindObject("weights_positive")));
      cfg.mAcceptance.push_back(reinterpret_cast<GFWWeights*>(listCorrections->FindObject("weights_negative")));
      int sizeAcc = cfg.mAcceptance.size();
      if (sizeAcc < 3)
        LOGF(warning, "Could not load acceptance weights from %s", cfgAcceptance.value.c_str());
      else
        LOGF(info, "Loaded acceptance weights from %s", cfgAcceptance.value.c_str());
    } else {
      LOGF(info, "cfgAcceptance empty! No corrections loaded");
    }
    if (cfgEfficiency.value.empty() == false) {
      TList* listCorrections = ccdb->getForTimeStamp<TList>(cfgEfficiency, timestamp);
      cfg.mEfficiency.push_back(reinterpret_cast<TH1D*>(listCorrections->FindObject("Efficiency")));
      cfg.mEfficiency.push_back(reinterpret_cast<TH1D*>(listCorrections->FindObject("Efficiency_pos")));
      cfg.mEfficiency.push_back(reinterpret_cast<TH1D*>(listCorrections->FindObject("Efficiency_neg")));
      int sizeEff = cfg.mEfficiency.size();
      if (sizeEff < 3) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", cfgEfficiency.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram from %s", cfgEfficiency.value.c_str());
    } else {
      LOGF(info, "cfgEfficiency empty! No corrections loaded");
    }
    cfg.correctionsLoaded = true;
  }

  // From Generic Framework
  bool setCurrentParticleWeights(int pID, float& weight_nue, float& weight_nua, const float& phi, const float& eta, const float& pt, const float& vtxz)
  {
    float eff = 1.;
    int sizeEff = cfg.mEfficiency.size();
    if (sizeEff > pID)
      eff = cfg.mEfficiency[pID]->GetBinContent(cfg.mEfficiency[pID]->FindBin(pt));
    else
      eff = 1.0;
    if (eff == 0)
      return false;
    weight_nue = 1. / eff;
    int sizeAcc = cfg.mAcceptance.size();
    if (sizeAcc > pID)
      weight_nua = cfg.mAcceptance[pID]->GetNUA(phi, eta, vtxz);
    else
      weight_nua = 1;
    return true;
  }

  template <typename TCollision>
  bool eventSelected(TCollision collision, const int& multTrk, const float& centrality)
  {
    if (!collision.sel8())
      return 0;
    registry.fill(HIST("hEventCount"), evSel_sel8);

    // Occupancy
    if (cfgDoOccupancySel) {
      auto occupancy = collision.trackOccupancyInTimeRange();
      if (occupancy > cfgMaxOccupancy) {
        return 0;
      }
      registry.fill(HIST("hEventCount"), evSel_occupancy);
    }

    if (cfgTVXinTRD) {
      if (collision.alias_bit(kTVXinTRD)) {
        // TRD triggered
        // "CMTVX-B-NOPF-TRD,minbias_TVX"
        return 0;
      }
      registry.fill(HIST("hEventCount"), evSel_kTVXinTRD);
    }

    if (cfgNoSameBunchPileupCut) {
      if (!collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
        // rejects collisions which are associated with the same "found-by-T0" bunch crossing
        // https://indico.cern.ch/event/1396220/#1-event-selection-with-its-rof
        return 0;
      }
      registry.fill(HIST("hEventCount"), evSel_kNoSameBunchPileup);
    }
    if (cfgIsGoodZvtxFT0vsPV) {
      if (!collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
        // removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference
        // use this cut at low multiplicities with caution
        return 0;
      }
      registry.fill(HIST("hEventCount"), evSel_kIsGoodZvtxFT0vsPV);
    }
    if (cfgNoCollInTimeRangeStandard) {
      if (!collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
        //  Rejection of the collisions which have other events nearby
        return 0;
      }
      registry.fill(HIST("hEventCount"), evSel_kNoCollInTimeRangeStandard);
    }

    if (cfgIsVertexITSTPC) {
      if (!collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
        // selects collisions with at least one ITS-TPC track, and thus rejects vertices built from ITS-only tracks
        return 0;
      }
      registry.fill(HIST("hEventCount"), evSel_kIsVertexITSTPC);
    }

    if (cfgUseAdditionalEventCut) {
      float vtxz = -999;
      if (collision.numContrib() > 1) {
        vtxz = collision.posZ();
        float zRes = std::sqrt(collision.covZZ());
        if (zRes > 0.25 && collision.numContrib() < 20)
          vtxz = -999;
      }

      auto multNTracksPV = collision.multNTracksPV();

      if (vtxz > 10 || vtxz < -10)
        return 0;
      if (multNTracksPV < fMultPVCutLow->Eval(centrality))
        return 0;
      if (multNTracksPV > fMultPVCutHigh->Eval(centrality))
        return 0;
      if (multTrk < fMultCutLow->Eval(centrality))
        return 0;
      if (multTrk > fMultCutHigh->Eval(centrality))
        return 0;

      registry.fill(HIST("hEventCount"), evSel_MultCuts);
    }

    if (centrality > cfgCentMax || centrality < cfgCentMin)
      return 0;
    registry.fill(HIST("hEventCount"), evSel_CentCuts);

    if (cfgIsGoodITSLayersAll) {
      if (!collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
        // New event selection bits to cut time intervals with dead ITS staves
        // https://indico.cern.ch/event/1493023/ (09-01-2025)
        return 0;
      }
      registry.fill(HIST("hEventCount"), evSel_kIsGoodITSLayersAll);
    }

    return 1;
  }

  template <typename TTrack>
  bool trackSelected(TTrack track, const int& field)
  {

    if (track.tpcNClsFound() < cfgNcls)
      return false;

    double phimodn = track.phi();
    if (field < 0) // for negative polarity field
      phimodn = o2::constants::math::TwoPI - phimodn;
    if (track.sign() < 0) // for negative charge
      phimodn = o2::constants::math::TwoPI - phimodn;
    if (phimodn < 0)
      LOGF(warning, "phi < 0: %g", phimodn);

    phimodn += o2::constants::math::PI / 18.0; // to center gap in the middle
    phimodn = fmod(phimodn, o2::constants::math::PI / 9.0);
    registry.fill(HIST("QA/before/pt_phi"), track.pt(), phimodn);

    if (cfgUseAdditionalTrackCut) {
      if (phimodn < fPhiCutHigh->Eval(track.pt()) && phimodn > fPhiCutLow->Eval(track.pt()))
        return false; // reject track
    }
    registry.fill(HIST("QA/after/pt_phi"), track.pt(), phimodn);
    return true;
  }

  template <typename CollisionObject, typename TracksObject>
  inline void fillEventQA(CollisionObject collision, TracksObject tracks, bool before)
  {
    if (before) {
      registry.fill(HIST("QA/before/hCent"), collision.centFT0C());
      registry.fill(HIST("QA/before/globalTracks_centT0C"), collision.centFT0C(), tracks.size());
      registry.fill(HIST("QA/before/PVTracks_centT0C"), collision.centFT0C(), collision.multNTracksPV());
      registry.fill(HIST("QA/before/globalTracks_PVTracks"), collision.multNTracksPV(), tracks.size());
      registry.fill(HIST("QA/before/globalTracks_multT0A"), collision.multFT0A(), tracks.size());
      registry.fill(HIST("QA/before/globalTracks_multV0A"), collision.multFV0A(), tracks.size());
      registry.fill(HIST("QA/before/multV0A_multT0A"), collision.multFT0A(), collision.multFV0A());
      registry.fill(HIST("QA/before/multT0C_centT0C"), collision.centFT0C(), collision.multFT0C());
    } else {
      registry.fill(HIST("QA/after/hCent"), collision.centFT0C());
      registry.fill(HIST("QA/after/globalTracks_centT0C"), collision.centFT0C(), tracks.size());
      registry.fill(HIST("QA/after/PVTracks_centT0C"), collision.centFT0C(), collision.multNTracksPV());
      registry.fill(HIST("QA/after/globalTracks_PVTracks"), collision.multNTracksPV(), tracks.size());
      registry.fill(HIST("QA/after/globalTracks_multT0A"), collision.multFT0A(), tracks.size());
      registry.fill(HIST("QA/after/globalTracks_multV0A"), collision.multFV0A(), tracks.size());
      registry.fill(HIST("QA/after/multV0A_multT0A"), collision.multFT0A(), collision.multFV0A());
      registry.fill(HIST("QA/after/multT0C_centT0C"), collision.centFT0C(), collision.multFT0C());

      double psiA = 1.0 * std::atan2(collision.qyA(), collision.qxA());
      double psiC = 1.0 * std::atan2(collision.qyC(), collision.qxC());
      double psiFull = 1.0 * std::atan2(collision.qyA() + collision.qyC(), collision.qxA() + collision.qxC());

      registry.fill(HIST("QA/after/PsiA_vs_Cent"), psiA, collision.centFT0C());
      registry.fill(HIST("QA/after/PsiC_vs_Cent"), psiC, collision.centFT0C());
      registry.fill(HIST("QA/after/PsiFull_vs_Cent"), psiFull, collision.centFT0C());
      registry.fill(HIST("QA/after/PsiA_vs_Vx"), psiA, collision.vx());
      registry.fill(HIST("QA/after/PsiC_vs_Vx"), psiC, collision.vx());
      registry.fill(HIST("QA/after/PsiFull_vs_Vx"), psiFull, collision.vx());
      registry.fill(HIST("QA/after/PsiA_vs_Vy"), psiA, collision.vy());
      registry.fill(HIST("QA/after/PsiC_vs_Vy"), psiC, collision.vy());
      registry.fill(HIST("QA/after/PsiFull_vs_Vy"), psiFull, collision.vy());
      registry.fill(HIST("QA/after/PsiA_vs_Vz"), psiA, collision.posZ());
      registry.fill(HIST("QA/after/PsiC_vs_Vz"), psiC, collision.posZ());
      registry.fill(HIST("QA/after/PsiFull_vs_Vz"), psiFull, collision.posZ());
    }
    return;
  }

  template <ChargeType ct, typename TrackObject>
  inline void fillHistograms(TrackObject track, float wacc, float weff, double ux, double uy, double uxMH, double uyMH, double qxA, double qyA, double qxC, double qyC, double corrQQx, double corrQQy, double corrQQ, double vnA, double vnC, double vnFull, double centrality, double vz)
  {
    static constexpr std::string_view Charge[] = {"incl/", "pos/", "neg/"};

    registry.fill(HIST(Charge[ct]) + HIST("QA/hPt"), track.pt());
    registry.fill(HIST(Charge[ct]) + HIST("QA/hPhi"), track.phi());
    registry.fill(HIST(Charge[ct]) + HIST("QA/hEta"), track.eta());
    registry.fill(HIST(Charge[ct]) + HIST("QA/hPhi_Eta_vz"), track.phi(), track.eta(), vz);
    registry.fill(HIST(Charge[ct]) + HIST("QA/hDCAxy"), track.dcaXY());
    registry.fill(HIST(Charge[ct]) + HIST("QA/hDCAz"), track.dcaZ());

    registry.fill(HIST(Charge[ct]) + HIST("vnAx_eta"), track.eta(), (ux * qxA) / std::sqrt(std::fabs(corrQQx)), wacc * weff);
    registry.fill(HIST(Charge[ct]) + HIST("vnAy_eta"), track.eta(), (uy * qyA) / std::sqrt(std::fabs(corrQQy)), wacc * weff);
    registry.fill(HIST(Charge[ct]) + HIST("vnCx_eta"), track.eta(), (ux * qxC) / std::sqrt(std::fabs(corrQQx)), wacc * weff);
    registry.fill(HIST(Charge[ct]) + HIST("vnCy_eta"), track.eta(), (uy * qyC) / std::sqrt(std::fabs(corrQQy)), wacc * weff);
    registry.fill(HIST(Charge[ct]) + HIST("vnA_eta"), track.eta(), (uy * qyA + ux * qxA) / std::sqrt(std::fabs(corrQQ)), wacc * weff);
    registry.fill(HIST(Charge[ct]) + HIST("vnC_eta"), track.eta(), (uy * qyC + ux * qxC) / std::sqrt(std::fabs(corrQQ)), wacc * weff);

    registry.fill(HIST(Charge[ct]) + HIST("vnAxCxUx_eta_MH"), track.eta(), (uxMH * qxA * qxC) / corrQQx, wacc * weff);
    registry.fill(HIST(Charge[ct]) + HIST("vnAxCyUx_eta_MH"), track.eta(), (uxMH * qyA * qyC) / corrQQy, wacc * weff);
    registry.fill(HIST(Charge[ct]) + HIST("vnAxCyUy_eta_MH"), track.eta(), (uyMH * qxA * qyC) / corrQQx, wacc * weff);
    registry.fill(HIST(Charge[ct]) + HIST("vnAyCxUy_eta_MH"), track.eta(), (uyMH * qyA * qxC) / corrQQy, wacc * weff);

    registry.fill(HIST(Charge[ct]) + HIST("vnAx_pt"), track.pt(), (ux * qxA) / std::sqrt(std::fabs(corrQQx)), wacc * weff);
    registry.fill(HIST(Charge[ct]) + HIST("vnAy_pt"), track.pt(), (uy * qyA) / std::sqrt(std::fabs(corrQQy)), wacc * weff);
    registry.fill(HIST(Charge[ct]) + HIST("vnCx_pt"), track.pt(), (ux * qxC) / std::sqrt(std::fabs(corrQQx)), wacc * weff);
    registry.fill(HIST(Charge[ct]) + HIST("vnCy_pt"), track.pt(), (uy * qyC) / std::sqrt(std::fabs(corrQQy)), wacc * weff);
    registry.fill(HIST(Charge[ct]) + HIST("vnA_pt"), track.pt(), (uy * qyA + ux * qxA) / std::sqrt(std::fabs(corrQQ)), wacc * weff);
    registry.fill(HIST(Charge[ct]) + HIST("vnC_pt"), track.pt(), (uy * qyC + ux * qxC) / std::sqrt(std::fabs(corrQQ)), wacc * weff);

    registry.fill(HIST(Charge[ct]) + HIST("vnAxCxUx_pt_MH"), track.pt(), (uxMH * qxA * qxC) / corrQQx, wacc * weff);
    registry.fill(HIST(Charge[ct]) + HIST("vnAxCyUx_pt_MH"), track.pt(), (uxMH * qyA * qyC) / corrQQy, wacc * weff);
    registry.fill(HIST(Charge[ct]) + HIST("vnAxCyUy_pt_MH"), track.pt(), (uyMH * qxA * qyC) / corrQQx, wacc * weff);
    registry.fill(HIST(Charge[ct]) + HIST("vnAyCxUy_pt_MH"), track.pt(), (uyMH * qyA * qxC) / corrQQy, wacc * weff);

    registry.fill(HIST(Charge[ct]) + HIST("vnA_eta_EP"), track.eta(), vnA, wacc * weff);
    registry.fill(HIST(Charge[ct]) + HIST("vnC_eta_EP"), track.eta(), vnC, wacc * weff);
    registry.fill(HIST(Charge[ct]) + HIST("vnFull_eta_EP"), track.eta(), vnFull, wacc * weff);

    registry.fill(HIST(Charge[ct]) + HIST("vnA_pt_EP"), track.pt(), vnA, wacc * weff);
    registry.fill(HIST(Charge[ct]) + HIST("vnC_pt_EP"), track.pt(), vnC, wacc * weff);
    registry.fill(HIST(Charge[ct]) + HIST("vnFull_pt_EP"), track.pt(), vnFull, wacc * weff);

    // For integrated v1 take only tracks from eta>0.
    // Following https://arxiv.org/pdf/1306.4145
    if (track.eta() < 0 && cfgHarm == 1) {
      registry.fill(HIST(Charge[ct]) + HIST("vnA_cent_minEta"), centrality, -1.0 * (uy * qyA + ux * qxA) / std::sqrt(std::fabs(corrQQ)), wacc * weff);
      registry.fill(HIST(Charge[ct]) + HIST("vnC_cent_minEta"), centrality, -1.0 * (uy * qyC + ux * qxC) / std::sqrt(std::fabs(corrQQ)), wacc * weff);

      registry.fill(HIST(Charge[ct]) + HIST("vnA_pt_odd"), track.pt(), -1.0 * (uy * qyA + ux * qxA) / std::sqrt(std::fabs(corrQQ)), wacc * weff);
      registry.fill(HIST(Charge[ct]) + HIST("vnC_pt_odd"), track.pt(), -1.0 * (uy * qyC + ux * qxC) / std::sqrt(std::fabs(corrQQ)), wacc * weff);
    } else {
      registry.fill(HIST(Charge[ct]) + HIST("vnA_cent_plusEta"), centrality, (uy * qyA + ux * qxA) / std::sqrt(std::fabs(corrQQ)), wacc * weff);
      registry.fill(HIST(Charge[ct]) + HIST("vnC_cent_plusEta"), centrality, (uy * qyC + ux * qxC) / std::sqrt(std::fabs(corrQQ)), wacc * weff);

      registry.fill(HIST(Charge[ct]) + HIST("vnA_pt_odd"), track.pt(), (uy * qyA + ux * qxA) / std::sqrt(std::fabs(corrQQ)), wacc * weff);
      registry.fill(HIST(Charge[ct]) + HIST("vnC_pt_odd"), track.pt(), (uy * qyC + ux * qxC) / std::sqrt(std::fabs(corrQQ)), wacc * weff);
    }

    registry.fill(HIST(Charge[ct]) + HIST("vnAxCxUx_cent_MH"), centrality, (uxMH * qxA * qxC) / corrQQx, wacc * weff);
    registry.fill(HIST(Charge[ct]) + HIST("vnAxCyUx_cent_MH"), centrality, (uxMH * qyA * qyC) / corrQQy, wacc * weff);
    registry.fill(HIST(Charge[ct]) + HIST("vnAxCyUy_cent_MH"), centrality, (uyMH * qxA * qyC) / corrQQx, wacc * weff);
    registry.fill(HIST(Charge[ct]) + HIST("vnAyCxUy_cent_MH"), centrality, (uyMH * qyA * qxC) / corrQQy, wacc * weff);

    registry.fill(HIST(Charge[ct]) + HIST("vnA_cent_EP"), centrality, vnA, wacc * weff);
    registry.fill(HIST(Charge[ct]) + HIST("vnC_cent_EP"), centrality, vnC, wacc * weff);
    registry.fill(HIST(Charge[ct]) + HIST("vnFull_cent_EP"), centrality, vnFull, wacc * weff);
  }

  void process(UsedCollisions::iterator const& collision, aod::BCsWithTimestamps const&, UsedTracks const& tracks)
  {
    registry.fill(HIST("hEventCount"), evSel_FilteredEvent);

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    auto field = (cfgMagField == 99999) ? getMagneticField(bc.timestamp()) : cfgMagField;

    if (bc.runNumber() != cfg.lastRunNumber) {
      // load corrections again for new run!
      cfg.correctionsLoaded = false;
      cfg.lastRunNumber = bc.runNumber();
    }

    loadCorrections(bc.timestamp());

    auto centrality = collision.centFT0C();

    if (!eventSelected(collision, tracks.size(), centrality))
      return;

    fillEventQA(collision, tracks, true);

    if (collision.isSelected()) {

      registry.fill(HIST("hEventCount"), evSel_isSelectedZDC);

      double qxA = collision.qxA();
      double qyA = collision.qyA();
      double qxC = collision.qxC();
      double qyC = collision.qyC();

      double vtxz = collision.posZ();

      double psiA = 1.0 * std::atan2(qyA, qxA);
      registry.fill(HIST("hSPplaneA"), psiA, 1);

      double psiC = 1.0 * std::atan2(qyC, qxC);
      registry.fill(HIST("hSPplaneC"), psiC, 1);

      // https://twiki.cern.ch/twiki/pub/ALICE/DirectedFlowAnalysisNote/vn_ZDC_ALICE_INT_NOTE_version02.pdf
      double psiFull = 1.0 * std::atan2(qyA + qyC, qxA + qxC);
      registry.fill(HIST("hSPplaneFull"), psiFull, 1);

      fillEventQA(collision, tracks, false);

      registry.fill(HIST("hCosPhiACosPhiC"), centrality, std::cos(psiA) * std::cos(psiC));
      registry.fill(HIST("hSinPhiASinPhiC"), centrality, std::sin(psiA) * std::sin(psiC));
      registry.fill(HIST("hSinPhiACosPhiC"), centrality, std::sin(psiA) * std::cos(psiC));
      registry.fill(HIST("hCosPhiASinsPhiC"), centrality, std::cos(psiA) * std::sin(psiC));

      registry.fill(HIST("hFullEvPlaneRes"), centrality, -1 * std::cos(psiA - psiC));

      registry.fill(HIST("qAqCXY"), centrality, qxA * qxC + qyA * qyC);
      registry.fill(HIST("qAqCX"), centrality, qxA * qxC);
      registry.fill(HIST("qAqCY"), centrality, qyA * qyC);

      double corrQQ = 1.;
      double corrQQx = 1.;
      double corrQQy = 1.;

      if (cfgLoadAverageQQ) {
        TList* hcorrList = ccdb->getForTimeStamp<TList>(cfgCCDBdir.value, bc.timestamp());
        TProfile* hcorrQQ = reinterpret_cast<TProfile*>(hcorrList->FindObject("qAqCXY"));
        TProfile* hcorrQQx = reinterpret_cast<TProfile*>(hcorrList->FindObject("qAqCX"));
        TProfile* hcorrQQy = reinterpret_cast<TProfile*>(hcorrList->FindObject("qAqCY"));
        corrQQ = hcorrQQ->GetBinContent(hcorrQQ->FindBin(centrality));
        corrQQx = hcorrQQx->GetBinContent(hcorrQQx->FindBin(centrality));
        corrQQy = hcorrQQy->GetBinContent(hcorrQQy->FindBin(centrality));
      }

      double evPlaneRes = 1.;
      if (cfgLoadSPPlaneRes) {
        TProfile* hEvPlaneRes = ccdb->getForTimeStamp<TProfile>(cfgCCDBdir_SP.value, bc.timestamp());
        evPlaneRes = hEvPlaneRes->GetBinContent(hEvPlaneRes->FindBin(centrality));
        if (evPlaneRes < 0)
          LOGF(fatal, "<Cos(PsiA-PsiC)> > 0 for centrality %.2f! Cannot determine resolution.. Change centrality ranges!!!", centrality);
        evPlaneRes = std::sqrt(evPlaneRes);
      }

      for (const auto& track : tracks) {
        registry.fill(HIST("QA/before/hPt_inclusive"), track.pt());

        float weff = 1, wacc = 1;
        float weffP = 1, waccP = 1;
        float weffN = 1, waccN = 1;

        if (!trackSelected(track, field))
          continue;

        if (track.sign() == 0.0)
          continue;
        bool pos = (track.sign() > 0) ? true : false;

        if (cfgFillWeights) {
          fWeights->Fill(track.phi(), track.eta(), vtxz, track.pt(), centrality, 0);
        } else if (cfgFillWeightsPOS) {
          if (pos)
            fWeightsPOS->Fill(track.phi(), track.eta(), vtxz, track.pt(), centrality, 0);
        } else if (cfgFillWeightsNEG) {
          if (!pos)
            fWeightsNEG->Fill(track.phi(), track.eta(), vtxz, track.pt(), centrality, 0);
        }

        if (!setCurrentParticleWeights(kInclusive, weff, wacc, track.phi(), track.eta(), track.pt(), vtxz))
          continue;

        if (pos && !setCurrentParticleWeights(kPositive, weffP, waccP, track.phi(), track.eta(), track.pt(), vtxz))
          continue;

        if (!pos && !setCurrentParticleWeights(kNegative, weffN, waccN, track.phi(), track.eta(), track.pt(), vtxz))
          continue;

        registry.fill(HIST("QA/after/hPt_inclusive"), track.pt());
        // // constrain angle to 0 -> [0,0+2pi]
        auto phi = RecoDecay::constrainAngle(track.phi(), 0);

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        auto ux = std::cos(cfgHarm * phi);
        auto uy = std::sin(cfgHarm * phi);

        auto uxMH = std::cos(cfgHarmMixed * phi);
        auto uyMH = std::sin(cfgHarmMixed * phi);
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        double vnA = std::cos(cfgHarm * (phi - psiA)) / evPlaneRes;
        double vnC = std::cos(cfgHarm * (phi - psiC)) / evPlaneRes;
        double vnFull = std::cos(cfgHarm * (phi - psiFull)) / evPlaneRes;

        fillHistograms<kInclusive>(track, wacc, weff, ux, uy, uxMH, uyMH, qxA, qyA, qxC, qyC, corrQQx, corrQQy, corrQQ, vnA, vnC, vnFull, centrality, vtxz);
        if (pos) {
          fillHistograms<kPositive>(track, waccP, weffP, ux, uy, uxMH, uyMH, qxA, qyA, qxC, qyC, corrQQx, corrQQy, corrQQ, vnA, vnC, vnFull, centrality, vtxz);
        } else {
          fillHistograms<kNegative>(track, waccN, weffN, ux, uy, uxMH, uyMH, qxA, qyA, qxC, qyC, corrQQx, corrQQy, corrQQ, vnA, vnC, vnFull, centrality, vtxz);
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowSP>(cfgc)};
}
