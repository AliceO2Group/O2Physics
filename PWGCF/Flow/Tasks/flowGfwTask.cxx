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

/// \file   flowGfwTask.cxx
/// \author Iris Likmeta (iris.likmeta@cern.ch)
/// \since  Mar 28, 2024
/// \brief  Multiparticle flow measurements with FT0 and ZDC

#include <CCDB/BasicCCDBManager.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <cmath>
#include <vector>
#include <string>
#include <memory>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/AnalysisDataModel.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"

#include "GFWPowerArray.h"
#include "GFW.h"
#include "GFWCumulant.h"
#include "GFWWeights.h"
#include "FlowContainer.h"
#include "TList.h"
#include <TProfile.h>
#include <TRandom3.h>
#include <TF1.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::track;
using namespace o2::aod::evsel;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct FlowGfwTask {

  O2_DEFINE_CONFIGURABLE(cfgCutVertex, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMin, float, 0.2f, "Minimal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMax, float, 3.0f, "Maximal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutChi2prTPCcls, float, 2.5, "Chi2 per TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgCutTPCclu, float, 70.0f, "minimum TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgCutITSclu, float, 5.0f, "minimum ITS clusters")
  O2_DEFINE_CONFIGURABLE(cfgUseAdditionalEventCut, bool, false, "Use additional event cut on mult correlations")
  O2_DEFINE_CONFIGURABLE(cfgUseAdditionalTrackCut, bool, false, "Use additional track cut on phi")
  O2_DEFINE_CONFIGURABLE(cfgUseNch, bool, false, "Use Nch for flow observables")
  O2_DEFINE_CONFIGURABLE(cfgNbootstrap, int, 10, "Number of subsamples")
  O2_DEFINE_CONFIGURABLE(cfgOutputNUAWeights, bool, false, "Fill and output NUA weights")
  O2_DEFINE_CONFIGURABLE(cfgEfficiency, std::string, "", "CCDB path to efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgAcceptance, std::string, "", "CCDB path to acceptance object")
  O2_DEFINE_CONFIGURABLE(cfgMagnetField, std::string, "GLO/Config/GRPMagField", "CCDB path to Magnet field object")
  O2_DEFINE_CONFIGURABLE(cfgCutOccupancyHigh, int, 500, "High cut on TPC occupancy")
  O2_DEFINE_CONFIGURABLE(cfgCutOccupancyLow, int, 0, "Low cut on TPC occupancy")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAz, float, 2, "Custom DCA Z cut")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAxy, float, 0.2f, "Custom DCA XY cut")
  O2_DEFINE_CONFIGURABLE(cfgTVXinTRD, bool, false, "Use kTVXinTRD (reject TRD triggered events)");
  O2_DEFINE_CONFIGURABLE(cfgNoTimeFrameBorder, bool, false, "kNoTimeFrameBorder");
  O2_DEFINE_CONFIGURABLE(cfgNoITSROFrameBorder, bool, false, "kNoITSROFrameBorder");
  O2_DEFINE_CONFIGURABLE(cfgNoSameBunchPileup, bool, false, "kNoSameBunchPileup");
  O2_DEFINE_CONFIGURABLE(cfgIsGoodZvtxFT0vsPV, bool, false, "kIsGoodZvtxFT0vsPV");
  O2_DEFINE_CONFIGURABLE(cfgNoCollInTimeRangeStandard, bool, false, "kNoCollInTimeRangeStandard");
  O2_DEFINE_CONFIGURABLE(cfgOccupancy, bool, false, "Bool for event selection on detector occupancy");
  O2_DEFINE_CONFIGURABLE(cfgMultCut, bool, false, "Use additional event cut on mult correlations");
  O2_DEFINE_CONFIGURABLE(FineBinning, bool, false, "Manually change to fine binning")
  O2_DEFINE_CONFIGURABLE(cfgTrackSelRun3ITSMatch, bool, false, "System check: Run3ITSMatch")
  O2_DEFINE_CONFIGURABLE(cfgTrackSel, bool, false, "System check: track selection")

  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for histograms"};
  ConfigurableAxis axisPhi{"axisPhi", {60, 0.0, constants::math::TwoPI}, "phi axis for histograms"};
  ConfigurableAxis axisPhiMod{"axisPhiMod", {100, 0, constants::math::PI / 9}, "fmod(#varphi,#pi/9)"};
  ConfigurableAxis axisEta{"axisEta", {40, -1., 1.}, "eta axis for histograms"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.2, 0.25, 0.30, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00}, "pt axis for histograms"};
  ConfigurableAxis axisPtHist{"axisPtHist", {100, 0., 10.}, "pt axis for histograms"};
  ConfigurableAxis axisCentrality{"axisCentrality", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100}, "centrality axis for histograms"};
  ConfigurableAxis axisNch{"axisNch", {4000, 0, 4000}, "N_{ch}"};
  ConfigurableAxis axisCentForQA{"axisCentForQA", {100, 0, 100}, "centrality for QA"};
  ConfigurableAxis axisT0C{"axisT0C", {70, 0, 70000}, "N_{ch} (T0C)"};
  ConfigurableAxis axisT0A{"axisT0A", {200, 0, 200000}, "N_{ch} (T0A)"};
  ConfigurableAxis axisFT0CAmp{"axisFT0CAmp", {5000, 0, 5000}, "axisFT0CAmp"};
  ConfigurableAxis axisFT0AAmp{"axisFT0AAmp", {5000, 0, 5000}, "axisFT0AAmp"};
  ConfigurableAxis axisFT0MAmp{"axisFT0MAmp", {10000, 0, 10000}, "axisFT0MAmp"};
  ConfigurableAxis axisNchPV{"axisNchPV", {4000, 0, 4000}, "N_{ch} (PV)"};
  ConfigurableAxis axisDCAz{"axisDCAz", {200, -2, 2}, "DCA_{z} (cm)"};
  ConfigurableAxis axisDCAxy{"axisDCAxy", {200, -1, 1}, "DCA_{xy} (cm)"};

  // Corrections
  TH1D* mEfficiency = nullptr;
  GFWWeights* mAcceptance = nullptr;
  bool correctionsLoaded = false;

  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<int64_t> ccdbNoLaterThan{"ccdbNoLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  // Define output
  OutputObj<FlowContainer> fFC{FlowContainer("FlowContainer")};
  OutputObj<GFWWeights> fWeights{GFWWeights("weights")};
  HistogramRegistry registry{"registry"};

  // define global variables
  GFW* fGFW = new GFW(); // GFW class used from main src
  std::vector<GFW::CorrConfig> corrconfigs;
  TRandom3* fRndm = new TRandom3(0);
  TAxis* fPtAxis;
  std::vector<std::vector<std::shared_ptr<TProfile>>> bootstrapArray; // TProfile is a shared pointer

  enum ExtraProfile {

    // here are TProfiles for vn-ft0 correlations that are not implemented in GFW
    kc22,
    kc24,
    kc26,
    kc28,
    kc22etagap,
    kc32,
    kc32etagap,
    kc34,
    kc22Nch,
    kc24Nch,
    kc26Nch,
    kc28Nch,
    kc22Nchetagap,
    kc32Nch,
    kc32Nchetagap,
    kc34Nch,
    kc22Nch05,
    kc24Nch05,
    kc26Nch05,
    kc28Nch05,
    kc22Nch05etagap,
    kc32Nch05,
    kc32Nch05etagap,
    kc34Nch05,

    // Count the total number of enum
    kCount_ExtraProfile
  };

  enum EventProgress {
    kFILTERED,
    kSEL8,
    kOCCUPANCY,
    kTVXINTRD,
    kNOTIMEFRAMEBORDER,
    kNOITSROFRAMEBORDER,
    kNOPSAMEBUNCHPILEUP,
    kISGOODZVTXFT0VSPV,
    kNOCOLLINTIMERANGESTANDART,
    kAFTERMULTCUTS,
    kCENTRALITY,
    kNOOFEVENTSTEPS
  };

  // Additional Event selection cuts - Copy from flowGenericFramework.cxx
  TrackSelection myTrackSel;
  TF1* fPhiCutLow = nullptr;
  TF1* fPhiCutHigh = nullptr;
  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;
  TF1* fMultCutLow = nullptr;
  TF1* fMultCutHigh = nullptr;
  TF1* fMultMultPVCut = nullptr;
  TF1* fT0AV0AMean = nullptr;
  TF1* fT0AV0ASigma = nullptr;

  void init(InitContext const&) // Initialization
  {
    ccdb->setURL(ccdbUrl.value);
    ccdb->setCaching(true);
    ccdb->setCreatedNotAfter(ccdbNoLaterThan.value);

    // Add some output objects to the histogram registry
    registry.add("hEventCount", "Number of Events;; No. of Events", {HistType::kTH1D, {{kNOOFEVENTSTEPS, -0.5, static_cast<int>(kNOOFEVENTSTEPS) - 0.5}}});
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kFILTERED + 1, "Filtered events");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kSEL8 + 1, "Sel8");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kOCCUPANCY + 1, "Occupancy");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kTVXINTRD + 1, "kTVXinTRD");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kNOTIMEFRAMEBORDER + 1, "kNoTimeFrameBorder");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kNOITSROFRAMEBORDER + 1, "kNoITSROFrameBorder");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kNOPSAMEBUNCHPILEUP + 1, "kNoSameBunchPileup");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kISGOODZVTXFT0VSPV + 1, "kIsGoodZvtxFT0vsPV");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kNOCOLLINTIMERANGESTANDART + 1, "kNoCollInTimeRangeStandard");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kAFTERMULTCUTS + 1, "After Mult cuts");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kCENTRALITY + 1, "Centrality");
    registry.add("hPhi", "#phi distribution", {HistType::kTH1D, {axisPhi}});
    registry.add("hPhiWeighted", "corrected #phi distribution", {HistType::kTH1D, {axisPhi}});
    registry.add("hEta", "", {HistType::kTH1D, {axisEta}});
    registry.add("hVtxZ", "Vexter Z distribution", {HistType::kTH1D, {axisVertex}});
    registry.add("hMult", "Multiplicity distribution", {HistType::kTH1D, {{3000, 0.5, 3000.5}}});
    registry.add("hCent", "Centrality distribution", {HistType::kTH1D, {{90, 0, 90}}});
    registry.add("cent_vs_Nch", ";Centrality (%); M (|#eta| < 0.8);", {HistType::kTH2D, {axisCentrality, axisNch}});

    // Before cuts
    registry.add("BeforeCut_globalTracks_centT0C", "before cut;Centrality T0C;mulplicity global tracks", {HistType::kTH2D, {axisCentForQA, axisNch}});
    registry.add("BeforeCut_PVTracks_centT0C", "before cut;Centrality T0C;mulplicity PV tracks", {HistType::kTH2D, {axisCentForQA, axisNchPV}});
    registry.add("BeforeCut_globalTracks_PVTracks", "before cut;mulplicity PV tracks;mulplicity global tracks", {HistType::kTH2D, {axisNchPV, axisNch}});
    registry.add("BeforeCut_globalTracks_multT0A", "before cut;mulplicity T0A;mulplicity global tracks", {HistType::kTH2D, {axisT0A, axisNch}});
    registry.add("BeforeCut_globalTracks_multV0A", "before cut;mulplicity V0A;mulplicity global tracks", {HistType::kTH2D, {axisT0A, axisNch}});
    registry.add("BeforeCut_multV0A_multT0A", "before cut;mulplicity T0A;mulplicity V0A", {HistType::kTH2D, {axisT0A, axisT0A}});
    registry.add("BeforeCut_multT0C_centT0C", "before cut;Centrality T0C;mulplicity T0C", {HistType::kTH2D, {axisCentForQA, axisT0C}});

    // After cuts
    registry.add("globalTracks_centT0C_Aft", "after cut;Centrality T0C;mulplicity global tracks", {HistType::kTH2D, {axisCentForQA, axisNch}});
    registry.add("PVTracks_centT0C_Aft", "after cut;Centrality T0C;mulplicity PV tracks", {HistType::kTH2D, {axisCentForQA, axisNchPV}});
    registry.add("globalTracks_PVTracks_Aft", "after cut;mulplicity PV tracks;mulplicity global tracks", {HistType::kTH2D, {axisNchPV, axisNch}});
    registry.add("globalTracks_multT0A_Aft", "after cut;mulplicity T0A;mulplicity global tracks", {HistType::kTH2D, {axisT0A, axisNch}});
    registry.add("globalTracks_multV0A_Aft", "after cut;mulplicity V0A;mulplicity global tracks", {HistType::kTH2D, {axisT0A, axisNch}});
    registry.add("multV0A_multT0A_Aft", "after cut;mulplicity T0A;mulplicity V0A", {HistType::kTH2D, {axisT0A, axisT0A}});
    registry.add("multT0C_centT0C_Aft", "after cut;Centrality T0C;mulplicity T0C", {HistType::kTH2D, {axisCentForQA, axisT0C}});

    // FT0 plots
    registry.add("FT0CAmp", ";FT0C amplitude;Events", kTH1F, {axisFT0CAmp});
    registry.add("FT0AAmp", ";FT0A amplitude;Events", kTH1F, {axisFT0AAmp});
    registry.add("FT0MAmp", ";FT0M amplitude;Events", kTH1F, {axisFT0MAmp});

    // Track plots
    registry.add("Events_per_Centrality_Bin", "Events_per_Centrality_Bin;Centrality FT0C;No. of Events", kTH1F, {axisCentrality});
    registry.add("Global_Tracks_Nch_vs_Cent", "Global Tracks;Centrality (%); M (|#eta| < 0.8);", {HistType::kTH2D, {axisCentrality, axisNch}});

    // Track QA
    registry.add("hPt", "p_{T} distribution before cut", {HistType::kTH1D, {axisPtHist}});
    registry.add("hPtRef", "p_{T} distribution after cut", {HistType::kTH1D, {axisPtHist}});
    registry.add("pt_phi_bef", "before cut;p_{T};#phi_{modn}", {HistType::kTH2D, {axisPt, axisPhiMod}});
    registry.add("pt_phi_aft", "after cut;p_{T};#phi_{modn}", {HistType::kTH2D, {axisPt, axisPhiMod}});
    registry.add("hChi2prTPCcls", "#chi^{2}/cluster for the TPC track segment", {HistType::kTH1D, {{100, 0., 5.}}});
    registry.add("hnTPCClu", "Number of found TPC clusters", {HistType::kTH1D, {{100, 40, 180}}});
    registry.add("hnTPCCrossedRow", "Number of crossed TPC Rows", {HistType::kTH1D, {{100, 40, 180}}});
    registry.add("hDCAz", "DCAz after cuts", {HistType::kTH1D, {{100, -3, 3}}});
    registry.add("hDCAxy", "DCAxy after cuts; DCAxy (cm); Pt", {HistType::kTH2D, {{50, -1, 1}, {50, 0, 10}}});

    // Additional Output histograms
    registry.add("c22", ";Centrality  (%) ; C_{2}{2} ", {HistType::kTProfile, {axisCentrality}});
    registry.add("c24", ";Centrality  (%) ; C_{2}{4}", {HistType::kTProfile, {axisCentrality}});
    registry.add("c26", ";Centrality  (%) ; C_{2}{6}", {HistType::kTProfile, {axisCentrality}});
    registry.add("c28", ";Centrality  (%) ; C_{2}{8}", {HistType::kTProfile, {axisCentrality}});
    registry.add("c22etagap", ";Centrality  (%) ; C_{2}{2} (|#eta| < 0.8) ", {HistType::kTProfile, {axisCentrality}});
    registry.add("c32", ";Centrality  (%) ; C_{3}{2} ", {HistType::kTProfile, {axisCentrality}});
    registry.add("c32etagap", ";Centrality  (%) ; C_{3}{2} (|#eta| < 0.8) ", {HistType::kTProfile, {axisCentrality}});
    registry.add("c34", ";Centrality  (%) ; C_{3}{4} ", {HistType::kTProfile, {axisCentrality}});

    registry.add("c22Nch", ";N_{ch}(|#eta| < 0.8) ; C_{2}{2} ", {HistType::kTProfile, {axisNch}});
    registry.add("c24Nch", ";N_{ch}(|#eta| < 0.8) ; C_{2}{4}", {HistType::kTProfile, {axisNch}});
    registry.add("c26Nch", ";N_{ch}(|#eta| < 0.8) ; C_{2}{6}", {HistType::kTProfile, {axisNch}});
    registry.add("c28Nch", ";N_{ch}(|#eta| < 0.8) ; C_{2}{8}", {HistType::kTProfile, {axisNch}});
    registry.add("c22Nchetagap", ";N_ch(|#eta| < 0.8) ; C_{2}{2} (|#eta| < 0.8) ", {HistType::kTProfile, {axisNch}});
    registry.add("c32Nch", ";N_{ch}(|#eta| < 0.8) ; C_{3}{2} ", {HistType::kTProfile, {axisNch}});
    registry.add("c32Nchetagap", ";N_ch(|#eta| < 0.8) ; C_{3}{2} (|#eta| < 0.8) ", {HistType::kTProfile, {axisNch}});
    registry.add("c34Nch", ";N_{ch}(|#eta| < 0.8) ; C_{3}{4} ", {HistType::kTProfile, {axisNch}});

    registry.add("c22Nch05", ";N_{ch 0-5%}(|#eta| < 0.8) ; C_{2}{2} ", {HistType::kTProfile, {axisNch}});
    registry.add("c24Nch05", ";N_{ch 0-5%}(|#eta| < 0.8) ; C_{2}{4}", {HistType::kTProfile, {axisNch}});
    registry.add("c26Nch05", ";N_{ch 0-5%}(|#eta| < 0.8) ; C_{2}{6}", {HistType::kTProfile, {axisNch}});
    registry.add("c28Nch05", ";N_{ch 0-5%}(|#eta| < 0.8) ; C_{2}{8}", {HistType::kTProfile, {axisNch}});
    registry.add("c22Nch05etagap", ";N_{ch 0-5%}(|#eta| < 0.8) ; C_{2}{2} (|#eta| < 0.8) ", {HistType::kTProfile, {axisNch}});
    registry.add("c32Nch05", ";N_{ch 0-5%}(|#eta| < 0.8) ; C_{3}{2} ", {HistType::kTProfile, {axisNch}});
    registry.add("c32Nch05etagap", ";N_{ch 0-5%}(|#eta| < 0.8) ; C_{3}{2} (|#eta| < 0.8) ", {HistType::kTProfile, {axisNch}});
    registry.add("c34Nch05", ";N_{ch 0-5%}(|#eta| < 0.8) ; C_{3}{4} ", {HistType::kTProfile, {axisNch}});

    // Create histograms for Reco and MC
    const AxisSpec axisCounter{1, 0, +1, ""};
    registry.add("eventCounter", "eventCounter", kTH1F, {axisCounter});
    registry.add("hPtMCRec", "Monte Carlo Reco", {HistType::kTH1D, {axisPt}});
    registry.add("hCenMCRec", "Monte Carlo Reco", {HistType::kTH1D, {axisCentrality}});
    registry.add("hPtMCRec05", "Monte Carlo Reco", {HistType::kTH1D, {axisPt}});
    registry.add("hPtMCRec5060", "Monte Carlo Reco", {HistType::kTH1D, {axisPt}});

    registry.add("mcEventCounter", "Monte Carlo Truth EventCounter", kTH1F, {axisCounter});
    registry.add("hPtMCGen", "Monte Carlo Truth", {HistType::kTH1D, {axisPt}});
    registry.add("hCenMCGen", "Monte Carlo Truth", {HistType::kTH1D, {axisCentrality}});
    registry.add("hPtMCGen05", "Monte Carlo Truth", {HistType::kTH1D, {axisPt}});
    registry.add("hPtMCGen5060", "Monte Carlo Truth", {HistType::kTH1D, {axisPt}});

    // initial array
    bootstrapArray.resize(cfgNbootstrap);
    for (int i = 0; i < cfgNbootstrap; i++) {
      bootstrapArray[i].resize(kCount_ExtraProfile);
    }

    for (int i = 0; i < cfgNbootstrap; i++) {
      bootstrapArray[i][kc22] = registry.add<TProfile>(Form("BootstrapContainer_%d/c22", i), ";Centrality  (%) ; C_{2}{2}", {HistType::kTProfile, {axisCentrality}});
      bootstrapArray[i][kc24] = registry.add<TProfile>(Form("BootstrapContainer_%d/c24", i), ";Centrality  (%) ; C_{2}{4}", {HistType::kTProfile, {axisCentrality}});
      bootstrapArray[i][kc26] = registry.add<TProfile>(Form("BootstrapContainer_%d/c26", i), ";Centrality  (%) ; C_{2}{6}", {HistType::kTProfile, {axisCentrality}});
      bootstrapArray[i][kc28] = registry.add<TProfile>(Form("BootstrapContainer_%d/c28", i), ";Centrality  (%) ; C_{2}{8}", {HistType::kTProfile, {axisCentrality}});
      bootstrapArray[i][kc22etagap] = registry.add<TProfile>(Form("BootstrapContainer_%d/c22etagap", i), ";Centrality (%) ; C_{2}{2} (|#eta| < 0.8)", {HistType::kTProfile, {axisCentrality}});
      bootstrapArray[i][kc32] = registry.add<TProfile>(Form("BootstrapContainer_%d/c32", i), ";Centrality  (%) ; C_{3}{2}", {HistType::kTProfile, {axisCentrality}});
      bootstrapArray[i][kc32etagap] = registry.add<TProfile>(Form("BootstrapContainer_%d/c32etagap", i), ";Centrality (%) ; C_{3}{2} (|#eta| < 0.8)", {HistType::kTProfile, {axisCentrality}});
      bootstrapArray[i][kc34] = registry.add<TProfile>(Form("BootstrapContainer_%d/c34", i), ";Centrality (%) ; C_{3}{4}", {HistType::kTProfile, {axisCentrality}});

      bootstrapArray[i][kc22Nch] = registry.add<TProfile>(Form("BootstrapContainer_%d/c22Nch", i), ";N_ch(|#eta| < 0.8) ; C_{2}{2}", {HistType::kTProfile, {axisNch}});
      bootstrapArray[i][kc24Nch] = registry.add<TProfile>(Form("BootstrapContainer_%d/c24Nch", i), ";N_ch(|#eta| < 0.8) ; C_{2}{4}", {HistType::kTProfile, {axisNch}});
      bootstrapArray[i][kc26Nch] = registry.add<TProfile>(Form("BootstrapContainer_%d/c26Nch", i), ";N_ch(|#eta| < 0.8) ; C_{2}{6}", {HistType::kTProfile, {axisNch}});
      bootstrapArray[i][kc28Nch] = registry.add<TProfile>(Form("BootstrapContainer_%d/c28Nch", i), ";N_ch(|#eta| < 0.8) ; C_{2}{8}", {HistType::kTProfile, {axisNch}});
      bootstrapArray[i][kc22Nchetagap] = registry.add<TProfile>(Form("BootstrapContainer_%d/c22Nchetagap", i), ";N_ch(|#eta| < 0.8) ; C_{2}{2} (|#eta| < 0.8)", {HistType::kTProfile, {axisNch}});
      bootstrapArray[i][kc32Nch] = registry.add<TProfile>(Form("BootstrapContainer_%d/c32Nch", i), ";N_ch(|#eta| < 0.8) ; C_{3}{2}", {HistType::kTProfile, {axisNch}});
      bootstrapArray[i][kc32Nchetagap] = registry.add<TProfile>(Form("BootstrapContainer_%d/c32Nchetagap", i), ";N_ch(|#eta| < 0.8) ; C_{3}{2} (|#eta| < 0.8)", {HistType::kTProfile, {axisNch}});
      bootstrapArray[i][kc34Nch] = registry.add<TProfile>(Form("BootstrapContainer_%d/c34Nch", i), ";N_ch(|#eta| < 0.8) ; C_{3}{4}", {HistType::kTProfile, {axisNch}});

      bootstrapArray[i][kc22Nch05] = registry.add<TProfile>(Form("BootstrapContainer_%d/c22Nch05", i), ";N_ch05(|#eta| < 0.8) ; C_{2}{2}", {HistType::kTProfile, {axisNch}});
      bootstrapArray[i][kc24Nch05] = registry.add<TProfile>(Form("BootstrapContainer_%d/c24Nch05", i), ";N_ch05(|#eta| < 0.8) ; C_{2}{4}", {HistType::kTProfile, {axisNch}});
      bootstrapArray[i][kc26Nch05] = registry.add<TProfile>(Form("BootstrapContainer_%d/c26Nch05", i), ";N_ch05(|#eta| < 0.8) ; C_{2}{6}", {HistType::kTProfile, {axisNch}});
      bootstrapArray[i][kc28Nch05] = registry.add<TProfile>(Form("BootstrapContainer_%d/c28Nch05", i), ";N_ch05(|#eta| < 0.8) ; C_{2}{8}", {HistType::kTProfile, {axisNch}});
      bootstrapArray[i][kc22Nch05etagap] = registry.add<TProfile>(Form("BootstrapContainer_%d/c22Nch05etagap", i), ";N_ch05(|#eta| < 0.8) ; C_{2}{2} (|#eta| < 0.8)", {HistType::kTProfile, {axisNch}});
      bootstrapArray[i][kc32Nch05] = registry.add<TProfile>(Form("BootstrapContainer_%d/c32Nch05", i), ";N_ch05(|#eta| < 0.8) ; C_{3}{2}", {HistType::kTProfile, {axisNch}});
      bootstrapArray[i][kc32Nch05etagap] = registry.add<TProfile>(Form("BootstrapContainer_%d/c32Nch05etagap", i), ";N_ch05(|#eta| < 0.8) ; C_{3}{2} (|#eta| < 0.8)", {HistType::kTProfile, {axisNch}});
      bootstrapArray[i][kc34Nch05] = registry.add<TProfile>(Form("BootstrapContainer_%d/c34Nch05", i), ";N_ch05(|#eta| < 0.8) ; C_{3}{4}", {HistType::kTProfile, {axisNch}});
    }

    o2::framework::AxisSpec axis = axisPt;
    int nPtBins = axis.binEdges.size() - 1;
    double* ptBins = &(axis.binEdges)[0];
    fPtAxis = new TAxis(nPtBins, ptBins);

    if (cfgOutputNUAWeights) {
      fWeights->SetPtBins(nPtBins, ptBins);
      fWeights->Init(true, false);
    }

    // add in FlowContainer to Get boostrap sample automatically
    TObjArray* oba = new TObjArray();
    fFC->SetXAxis(fPtAxis);
    fFC->SetName("FlowContainer");
    fFC->Initialize(oba, axisCentrality, cfgNbootstrap);
    delete oba;

    fGFW->AddRegion("full", -0.8, 0.8, 1, 1); // eta region -0.8 to 0.8
    fGFW->AddRegion("refN10", -0.8, -0.5, 1, 1);
    fGFW->AddRegion("refP10", 0.5, 0.8, 1, 1);

    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {2 -2}", "ChFull22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {2 2 -2 -2}", "ChFull24", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {2 2 2 -2 -2 -2}", "ChFull26", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {2 2 2 2  -2 -2 -2 -2}", "ChFull28", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN10 {2} refP10 {-2}", "Ch10Gap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {3 -3}", "ChFull32", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN10 {3} refP10 {-3}", "Ch10Gap32", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {3 3 -3 -3}", "ChFull34", kFALSE));
    fGFW->CreateRegions(); // finalize the initialization

    if (cfgUseAdditionalEventCut) {
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

    if (cfgUseAdditionalTrackCut) {
      fPhiCutLow = new TF1("fPhiCutLow", "0.06/x+pi/18.0-0.06", 0, 100);
      fPhiCutHigh = new TF1("fPhiCutHigh", "0.1/x+pi/18.0+0.06", 0, 100);
    }

    if (cfgTrackSelRun3ITSMatch) {
      myTrackSel = getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSall7Layers, TrackSelection::GlobalTrackRun3DCAxyCut::Default);
    } else {
      myTrackSel = getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSibAny, TrackSelection::GlobalTrackRun3DCAxyCut::Default);
    }

    myTrackSel.SetMinNClustersTPC(cfgCutTPCclu);
    myTrackSel.SetMinNClustersITS(cfgCutITSclu);

  } // end of Initialization

  template <char... chars>
  void fillProfile(const GFW::CorrConfig& corrconf, const ConstStr<chars...>& tarName, const double& cent)
  {
    double dnx, val;
    dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
    if (dnx == 0)
      return;
    if (!corrconf.pTDif) {
      val = fGFW->Calculate(corrconf, 0, kFALSE).real() / dnx;
      if (std::abs(val) < 1)
        registry.fill(tarName, cent, val, dnx);
      return;
    }
    return;
  }

  void fillProfile(const GFW::CorrConfig& corrconf, std::shared_ptr<TProfile> tarName, const double& cent)
  {
    double dnx, val;
    dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
    if (dnx == 0)
      return;
    if (!corrconf.pTDif) {
      val = fGFW->Calculate(corrconf, 0, kFALSE).real() / dnx;
      if (std::abs(val) < 1) {
        tarName->Fill(cent, val, dnx);
      }
      return;
    }
    return;
  }

  void fillFC(const GFW::CorrConfig& corrconf, const double& cent, const double& rndm)
  {
    double dnx, val;
    dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
    if (dnx == 0)
      return;
    if (!corrconf.pTDif) {
      val = fGFW->Calculate(corrconf, 0, kFALSE).real() / dnx;
      if (std::abs(val) < 1)
        fFC->FillProfile(corrconf.Head.c_str(), cent, val, dnx, rndm);
      return;
    }
    for (int i = 1; i <= fPtAxis->GetNbins(); i++) {
      dnx = fGFW->Calculate(corrconf, i - 1, kTRUE).real();
      if (dnx == 0)
        continue;
      val = fGFW->Calculate(corrconf, i - 1, kFALSE).real() / dnx;
      if (std::abs(val) < 1)
        fFC->FillProfile(Form("%s_pt_%i", corrconf.Head.c_str(), i), cent, val, dnx, rndm);
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

  template <typename TCollision>
  bool eventSelected(o2::aod::mult::MultNTracksPV, TCollision collision, const int multTrk, const float centrality)
  {
    if (cfgTVXinTRD) {
      if (collision.alias_bit(kTVXinTRD)) {
        // TRD triggered
        return false;
      }
      registry.fill(HIST("hEventCount"), kTVXINTRD);
    }
    if (cfgNoTimeFrameBorder) {
      if (!collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
        // reject collisions close to Time Frame borders
        // https://its.cern.ch/jira/browse/O2-4623
        return false;
      }
      registry.fill(HIST("hEventCount"), kNOTIMEFRAMEBORDER);
    }
    if (cfgNoITSROFrameBorder) {
      if (!collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
        // reject events affected by the ITS ROF border
        // https://its.cern.ch/jira/browse/O2-4309
        return false;
      }
      registry.fill(HIST("hEventCount"), kNOITSROFRAMEBORDER);
    }
    if (cfgNoSameBunchPileup) {
      if (!collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
        // rejects collisions which are associated with the same "found-by-T0" bunch crossing
        // https://indico.cern.ch/event/1396220/#1-event-selection-with-its-rof
        return false;
      }
      registry.fill(HIST("hEventCount"), kNOPSAMEBUNCHPILEUP);
    }
    if (cfgIsGoodZvtxFT0vsPV) {
      if (!collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
        // removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference
        // use this cut at low multiplicities with caution
        return false;
      }
      registry.fill(HIST("hEventCount"), kISGOODZVTXFT0VSPV);
    }
    if (cfgNoCollInTimeRangeStandard) {
      if (!collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
        // no collisions in specified time range
        return false;
      }
      registry.fill(HIST("hEventCount"), kNOCOLLINTIMERANGESTANDART);
    }

    float vtxz = -999;
    if (collision.numContrib() > 1) {
      vtxz = collision.posZ();
      float zRes = std::sqrt(collision.covZZ());
      if (zRes > 0.25 && collision.numContrib() < 20)
        vtxz = -999;
    }

    auto multNTracksPV = collision.multNTracksPV();

    if (std::abs(vtxz) > cfgCutVertex)
      return false;

    if (cfgMultCut) {
      if (multNTracksPV < fMultPVCutLow->Eval(centrality))
        return false;
      if (multNTracksPV > fMultPVCutHigh->Eval(centrality))
        return false;
      if (multTrk < fMultCutLow->Eval(centrality))
        return false;
      if (multTrk > fMultCutHigh->Eval(centrality))
        return false;
      registry.fill(HIST("hEventCount"), kAFTERMULTCUTS);
    }

    // V0A T0A 5 sigma cut
    if (std::abs(collision.multFV0A() - fT0AV0AMean->Eval(collision.multFT0A())) > 5 * fT0AV0ASigma->Eval(collision.multFT0A()))
      return false;

    return true;
  }

  int getMagneticField(uint64_t timestamp)
  {
    static o2::parameters::GRPMagField* grpo = nullptr;
    if (grpo == nullptr) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(cfgMagnetField, timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found in %s for timestamp %llu", cfgMagnetField.value.c_str(), timestamp);
        return 0;
      }
      LOGF(info, "Retrieved GRP from %s for timestamp %llu with magnetic field of %d kG", cfgMagnetField.value.c_str(), timestamp, grpo->getNominalL3Field());
    }
    return grpo->getNominalL3Field();
  }

  template <typename TTrack>
  bool trackSelected(TTrack track, const int field)
  {
    double phimodn = track.phi();
    if (field < 0) // for negative polarity field
      phimodn = o2::constants::math::TwoPI - phimodn;
    if (track.sign() < 0) // for negative charge
      phimodn = o2::constants::math::TwoPI - phimodn;
    if (phimodn < 0)
      LOGF(warning, "phi < 0: %g", phimodn);

    phimodn += o2::constants::math::PI / 18.0; // to center gap in the middle
    phimodn = fmod(phimodn, o2::constants::math::PI / 9.0);
    registry.fill(HIST("pt_phi_bef"), track.pt(), phimodn);
    if (phimodn < fPhiCutHigh->Eval(track.pt()) && phimodn > fPhiCutLow->Eval(track.pt()))
      return false; // reject track
    registry.fill(HIST("pt_phi_aft"), track.pt(), phimodn);
    return true;
  }

  template <typename TTrack>
  bool trackSelected(TTrack track)
  {

    if (cfgTrackSel) {
      return myTrackSel.IsSelected(track);
    } else {
      return (track.tpcNClsFound() >= cfgCutTPCclu);
    }
  }

  // Apply process filters
  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtMin) && (aod::track::pt < cfgCutPtMax) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && (aod::track::tpcChi2NCl < cfgCutChi2prTPCcls) && (nabs(aod::track::dcaZ) < cfgCutDCAz) && (nabs(aod::track::dcaXY) < cfgCutDCAxy);

  using Colls = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs>>;               // collisions filter
  using AodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksDCA, aod::TracksExtra>>; // tracks filter

  void processData(Colls::iterator const& collision, aod::BCsWithTimestamps const&, AodTracks const& tracks, aod::FT0s const&)
  {
    registry.fill(HIST("hEventCount"), kFILTERED);
    if (!collision.sel8())
      return;

    int nTotal = tracks.size();
    if (nTotal < 1)
      return;

    // fill event QA before cuts
    registry.fill(HIST("BeforeCut_globalTracks_centT0C"), collision.centFT0C(), tracks.size());
    registry.fill(HIST("BeforeCut_PVTracks_centT0C"), collision.centFT0C(), collision.multNTracksPV());
    registry.fill(HIST("BeforeCut_globalTracks_PVTracks"), collision.multNTracksPV(), tracks.size());
    registry.fill(HIST("BeforeCut_globalTracks_multT0A"), collision.multFT0A(), tracks.size());
    registry.fill(HIST("BeforeCut_globalTracks_multV0A"), collision.multFV0A(), tracks.size());
    registry.fill(HIST("BeforeCut_multV0A_multT0A"), collision.multFT0A(), collision.multFV0A());
    registry.fill(HIST("BeforeCut_multT0C_centT0C"), collision.centFT0C(), collision.multFT0C());
    registry.fill(HIST("hEventCount"), kSEL8);

    const auto centrality = collision.centFT0C();

    if (cfgOccupancy) {
      int occupancy = collision.trackOccupancyInTimeRange();
      if (occupancy < cfgCutOccupancyLow || occupancy > cfgCutOccupancyHigh)
        return;
      registry.fill(HIST("hEventCount"), kOCCUPANCY);
    }

    if (cfgUseAdditionalEventCut && !eventSelected(o2::aod::mult::MultNTracksPV(), collision, tracks.size(), centrality)) {
      return;
    }

    if (centrality < 0 || centrality >= 70.)
      return;

    float vtxz = collision.posZ();
    float lRandom = fRndm->Rndm();
    registry.fill(HIST("hVtxZ"), vtxz);
    registry.fill(HIST("hMult"), nTotal);
    registry.fill(HIST("hCent"), centrality);
    registry.fill(HIST("cent_vs_Nch"), centrality, nTotal);

    fGFW->Clear();

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    loadCorrections(bc.timestamp());
    registry.fill(HIST("hEventCount"), kCENTRALITY);

    // fill event QA after cuts
    registry.fill(HIST("globalTracks_centT0C_Aft"), collision.centFT0C(), tracks.size());
    registry.fill(HIST("PVTracks_centT0C_Aft"), collision.centFT0C(), collision.multNTracksPV());
    registry.fill(HIST("globalTracks_PVTracks_Aft"), collision.multNTracksPV(), tracks.size());
    registry.fill(HIST("globalTracks_multT0A_Aft"), collision.multFT0A(), tracks.size());
    registry.fill(HIST("globalTracks_multV0A_Aft"), collision.multFV0A(), tracks.size());
    registry.fill(HIST("multV0A_multT0A_Aft"), collision.multFT0A(), collision.multFV0A());
    registry.fill(HIST("multT0C_centT0C_Aft"), collision.centFT0C(), collision.multFT0C());

    // FT0 amplitude to use in fine binning
    double ft0aAmp = 0;
    double ft0cAmp = 0;

    if (collision.has_foundFT0()) {
      auto ft0 = collision.foundFT0();
      for (const auto& amplitude : ft0.amplitudeA()) {
        ft0aAmp += amplitude;
      }
      for (const auto& amplitude : ft0.amplitudeC()) {
        ft0cAmp += amplitude;
      }
    }

    registry.fill(HIST("FT0AAmp"), ft0aAmp);
    registry.fill(HIST("FT0CAmp"), ft0cAmp);

    double ft0mAmp = ft0aAmp + ft0cAmp;
    registry.fill(HIST("FT0MAmp"), ft0mAmp);

    // track weights
    float weff = 1, wacc = 1;
    int magnetfield = 0;

    if (cfgUseAdditionalTrackCut) {
      // magnet field dependence cut
      magnetfield = getMagneticField(bc.timestamp());
    }

    // track loop
    int globalTracksNch = 0;

    for (const auto& track : tracks) {

      if (!trackSelected(track))
        continue;

      if (cfgUseAdditionalTrackCut && !trackSelected(track, magnetfield))
        continue;

      if (cfgOutputNUAWeights)
        fWeights->Fill(track.phi(), track.eta(), vtxz, track.pt(), centrality, 0);

      if (!setCurrentParticleWeights(weff, wacc, track.phi(), track.eta(), track.pt(), vtxz))
        continue;

      bool withinPtRef = (cfgCutPtMin < track.pt()) && (track.pt() < cfgCutPtMax); // within RF pT range
      registry.fill(HIST("hPt"), track.pt());

      if (withinPtRef) {
        registry.fill(HIST("hPhi"), track.phi());
        registry.fill(HIST("hPhiWeighted"), track.phi(), wacc);
        registry.fill(HIST("hEta"), track.eta());
        registry.fill(HIST("hPtRef"), track.pt());
        registry.fill(HIST("hChi2prTPCcls"), track.tpcChi2NCl());
        registry.fill(HIST("hnTPCClu"), track.tpcNClsFound());
        registry.fill(HIST("hnTPCCrossedRow"), track.tpcNClsCrossedRows());
        registry.fill(HIST("hDCAz"), track.dcaZ());
        registry.fill(HIST("hDCAxy"), track.dcaXY(), track.pt());
      }

      globalTracksNch++;

      if (withinPtRef)
        fGFW->Fill(track.eta(), 1, track.phi(), wacc * weff, 1);

      if (FineBinning == true)
        fGFW->Fill(track.eta(), 1, track.phi(), wacc * weff, 1);

    } // End of track loop

    registry.fill(HIST("Events_per_Centrality_Bin"), centrality);
    registry.fill(HIST("Global_Tracks_Nch_vs_Cent"), centrality, globalTracksNch);

    // Filling c22 with ROOT TProfile
    fillProfile(corrconfigs.at(0), HIST("c22"), centrality);
    fillProfile(corrconfigs.at(1), HIST("c24"), centrality);
    fillProfile(corrconfigs.at(2), HIST("c26"), centrality);
    fillProfile(corrconfigs.at(3), HIST("c28"), centrality);
    fillProfile(corrconfigs.at(4), HIST("c22etagap"), centrality);
    fillProfile(corrconfigs.at(5), HIST("c32"), centrality);
    fillProfile(corrconfigs.at(6), HIST("c32etagap"), centrality);
    fillProfile(corrconfigs.at(7), HIST("c34"), centrality);

    fillProfile(corrconfigs.at(0), HIST("c22Nch"), globalTracksNch);
    fillProfile(corrconfigs.at(1), HIST("c24Nch"), globalTracksNch);
    fillProfile(corrconfigs.at(2), HIST("c26Nch"), globalTracksNch);
    fillProfile(corrconfigs.at(3), HIST("c28Nch"), globalTracksNch);
    fillProfile(corrconfigs.at(4), HIST("c22Nchetagap"), globalTracksNch);
    fillProfile(corrconfigs.at(5), HIST("c32Nch"), globalTracksNch);
    fillProfile(corrconfigs.at(6), HIST("c32Nchetagap"), globalTracksNch);
    fillProfile(corrconfigs.at(7), HIST("c34Nch"), globalTracksNch);

    // 0-5% centrality Nch
    if (centrality >= 0 && centrality <= 5) {
      fillProfile(corrconfigs.at(0), HIST("c22Nch05"), globalTracksNch);
      fillProfile(corrconfigs.at(1), HIST("c24Nch05"), globalTracksNch);
      fillProfile(corrconfigs.at(2), HIST("c26Nch05"), globalTracksNch);
      fillProfile(corrconfigs.at(3), HIST("c28Nch05"), globalTracksNch);
      fillProfile(corrconfigs.at(4), HIST("c22Nch05etagap"), globalTracksNch);
      fillProfile(corrconfigs.at(5), HIST("c32Nch05"), globalTracksNch);
      fillProfile(corrconfigs.at(6), HIST("c32Nch05etagap"), globalTracksNch);
      fillProfile(corrconfigs.at(7), HIST("c34Nch05"), globalTracksNch);
    }

    // Filling Bootstrap Samples
    int sampleIndex = static_cast<int>(cfgNbootstrap * lRandom);
    fillProfile(corrconfigs.at(0), bootstrapArray[sampleIndex][kc22], centrality);
    fillProfile(corrconfigs.at(1), bootstrapArray[sampleIndex][kc24], centrality);
    fillProfile(corrconfigs.at(2), bootstrapArray[sampleIndex][kc26], centrality);
    fillProfile(corrconfigs.at(3), bootstrapArray[sampleIndex][kc28], centrality);
    fillProfile(corrconfigs.at(4), bootstrapArray[sampleIndex][kc22etagap], centrality);
    fillProfile(corrconfigs.at(5), bootstrapArray[sampleIndex][kc32], centrality);
    fillProfile(corrconfigs.at(6), bootstrapArray[sampleIndex][kc32etagap], centrality);
    fillProfile(corrconfigs.at(7), bootstrapArray[sampleIndex][kc34], centrality);

    fillProfile(corrconfigs.at(0), bootstrapArray[sampleIndex][kc22Nch], globalTracksNch);
    fillProfile(corrconfigs.at(1), bootstrapArray[sampleIndex][kc24Nch], globalTracksNch);
    fillProfile(corrconfigs.at(2), bootstrapArray[sampleIndex][kc26Nch], globalTracksNch);
    fillProfile(corrconfigs.at(3), bootstrapArray[sampleIndex][kc28Nch], globalTracksNch);
    fillProfile(corrconfigs.at(4), bootstrapArray[sampleIndex][kc22Nchetagap], globalTracksNch);
    fillProfile(corrconfigs.at(5), bootstrapArray[sampleIndex][kc32Nch], globalTracksNch);
    fillProfile(corrconfigs.at(6), bootstrapArray[sampleIndex][kc32Nchetagap], globalTracksNch);
    fillProfile(corrconfigs.at(7), bootstrapArray[sampleIndex][kc34Nch], globalTracksNch);

    if (centrality >= 0 && centrality <= 5) {
      fillProfile(corrconfigs.at(0), bootstrapArray[sampleIndex][kc22Nch05], globalTracksNch);
      fillProfile(corrconfigs.at(1), bootstrapArray[sampleIndex][kc24Nch05], globalTracksNch);
      fillProfile(corrconfigs.at(2), bootstrapArray[sampleIndex][kc26Nch05], globalTracksNch);
      fillProfile(corrconfigs.at(3), bootstrapArray[sampleIndex][kc28Nch05], globalTracksNch);
      fillProfile(corrconfigs.at(4), bootstrapArray[sampleIndex][kc22Nch05etagap], globalTracksNch);
      fillProfile(corrconfigs.at(5), bootstrapArray[sampleIndex][kc32Nch05], globalTracksNch);
      fillProfile(corrconfigs.at(6), bootstrapArray[sampleIndex][kc32Nch05etagap], globalTracksNch);
      fillProfile(corrconfigs.at(7), bootstrapArray[sampleIndex][kc34Nch05], globalTracksNch);
    }

    // Filling Flow Container
    for (uint l_ind = 0; l_ind < corrconfigs.size(); l_ind++) {
      fillFC(corrconfigs.at(l_ind), centrality, lRandom);
    }

  } // End of process
  PROCESS_SWITCH(FlowGfwTask, processData, "Process analysis for Run 3 data", false);

  // Filter the Reconstructed tracks
  Filter mytrackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtMin) && (aod::track::pt < cfgCutPtMax) && (nabs(aod::track::dcaXY) < cfgCutDCAxy);
  using MyTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::McTrackLabels>>;
  using MyCollisions = soa::Join<aod::Collisions, aod::CentFT0Cs>;

  void processMCReco(MyCollisions::iterator const& collision, MyTracks const& tracks, aod::McParticles const&)
  {
    registry.fill(HIST("eventCounter"), 0.5);
    const auto centrality = collision.centFT0C();
    registry.fill(HIST("hCenMCRec"), centrality);
    for (const auto& track : tracks) {
      if (track.tpcNClsCrossedRows() < 70)
        continue;

      if (track.has_mcParticle()) {
        registry.fill(HIST("hPtMCRec"), track.pt());
        if (centrality > 0 && centrality <= 5) {
          registry.fill(HIST("hPtMCRec05"), track.pt());
        }
        if (centrality >= 50 && centrality <= 60) {
          registry.fill(HIST("hPtMCRec5060"), track.pt());
        }
      }
    }
  }
  PROCESS_SWITCH(FlowGfwTask, processMCReco, "process reconstructed information", false);

  // Filter for MCParticle simulation
  Filter particleFilter = (nabs(aod::mcparticle::eta) < cfgCutEta) && (aod::mcparticle::pt > cfgCutPtMin) && (aod::mcparticle::pt < cfgCutPtMax);
  using MyMcParticles = soa::Filtered<aod::McParticles>;
  using MyMcCollisionsFT0Cs = soa::Join<o2::aod::Collisions, o2::aod::CentFT0Cs>;

  void processMCGEN(aod::McCollision const&, soa::SmallGroups<soa::Join<o2::aod::Collisions, o2::aod::McCollisionLabels>> const& collisions, MyMcParticles const& mcParticles, MyMcCollisionsFT0Cs const& mcCollisionsFT0Cs)
  {
    if (collisions.size() > -1) {
      registry.fill(HIST("mcEventCounter"), 0.5);
      for (const auto& mcCollisionsFT0C : mcCollisionsFT0Cs) {
        registry.fill(HIST("hCenMCGen"), mcCollisionsFT0C.centFT0C());
      }

      for (const auto& mcCollisionsFT0C : mcCollisionsFT0Cs) {
        const auto centrality = mcCollisionsFT0C.centFT0C();
        for (const auto& mcParticle : mcParticles) {
          registry.fill(HIST("hPtMCGen"), mcParticle.pt());
          if (centrality > 0 && centrality <= 5) {
            registry.fill(HIST("hPtMCGen05"), mcParticle.pt());
          }
          if (centrality >= 50 && centrality <= 60) {
            registry.fill(HIST("hPtMCGen5060"), mcParticle.pt());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(FlowGfwTask, processMCGEN, "process pure simulation information", false);

}; // End of struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowGfwTask>(cfgc)};
}
