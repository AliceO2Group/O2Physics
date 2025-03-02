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
  O2_DEFINE_CONFIGURABLE(cfgTrackSel, bool, false, "ITS and TPC cluster selection")
  O2_DEFINE_CONFIGURABLE(cfgMinCentFT0C, float, 0.0f, "Minimum FT0C Centrality")
  O2_DEFINE_CONFIGURABLE(cfgMaxCentFT0C, float, 100.0f, "Maximum FT0C Centrality")
  O2_DEFINE_CONFIGURABLE(cfgcentEstFt0c, bool, false, "Centrality estimator based on FT0C signal")
  O2_DEFINE_CONFIGURABLE(cfgcentEstFt0a, bool, false, "Centrality estimator based on FT0A signal")
  O2_DEFINE_CONFIGURABLE(cfgcentEstFt0m, bool, false, " A centrality estimator based on FT0A+FT0C signals.")
  O2_DEFINE_CONFIGURABLE(cfgcentEstFv0a, bool, false, "Centrality estimator based on FV0A signal")
  O2_DEFINE_CONFIGURABLE(cfgcentEstFt0cVariant1, bool, false, "A variant of FT0C")
  O2_DEFINE_CONFIGURABLE(cfgUseAdditionalEventCut, bool, false, "Use additional event cut on mult correlations")
  O2_DEFINE_CONFIGURABLE(cfgUseAdditionalTrackCut, bool, false, "Use additional track cut on phi")
  O2_DEFINE_CONFIGURABLE(cfgTrackSelRun3ITSMatch, bool, false, "Track selection for ITS matches")
  O2_DEFINE_CONFIGURABLE(cfgUseNch, bool, false, "Use Nch for flow observables")
  O2_DEFINE_CONFIGURABLE(cfgNbootstrap, int, 10, "Number of subsamples")
  O2_DEFINE_CONFIGURABLE(cfgOutputNUAWeights, bool, false, "Fill and output NUA weights")
  O2_DEFINE_CONFIGURABLE(cfgEfficiency, std::string, "", "CCDB path to efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgAcceptance, std::string, "", "CCDB path to acceptance object")
  O2_DEFINE_CONFIGURABLE(cfgMagnetField, std::string, "GLO/Config/GRPMagField", "CCDB path to Magnet field object")
  O2_DEFINE_CONFIGURABLE(cfgCutOccupancyHigh, int, 3000, "High cut on TPC occupancy")
  O2_DEFINE_CONFIGURABLE(cfgCutOccupancyLow, int, 0, "Low cut on TPC occupancy")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAz, float, 2, "Custom DCA Z cut")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAxy, float, 0.2f, "Custom DCA XY cut")
  O2_DEFINE_CONFIGURABLE(cfgDCAzPt, bool, false, "switch for DCAz pt dependent")
  O2_DEFINE_CONFIGURABLE(cfgNoTimeFrameBorder, bool, false, "kNoTimeFrameBorder");
  O2_DEFINE_CONFIGURABLE(cfgNoITSROFrameBorder, bool, false, "kNoITSROFrameBorder");
  O2_DEFINE_CONFIGURABLE(cfgNoSameBunchPileup, bool, false, "kNoSameBunchPileup");
  O2_DEFINE_CONFIGURABLE(cfgIsGoodZvtxFT0vsPV, bool, false, "kIsGoodZvtxFT0vsPV");
  O2_DEFINE_CONFIGURABLE(cfgIsVertexITSTPC, bool, false, "kIsVertexITSTPC");
  O2_DEFINE_CONFIGURABLE(cfgNoCollInTimeRangeStandard, bool, false, "kNoCollInTimeRangeStandard");
  O2_DEFINE_CONFIGURABLE(cfgEvSelkIsGoodITSLayersAll, bool, false, "kIsGoodITSLayersAll")
  O2_DEFINE_CONFIGURABLE(cfgOccupancy, bool, false, "Bool for event selection on detector occupancy");
  O2_DEFINE_CONFIGURABLE(cfgMultCut, bool, false, "Use additional event cut on mult correlations");
  O2_DEFINE_CONFIGURABLE(cfgV0AT0A5Sigma, bool, true, "V0A T0A 5 sigma cut")
  O2_DEFINE_CONFIGURABLE(cfgGlobalplusITS, bool, false, "Global and ITS tracks")
  O2_DEFINE_CONFIGURABLE(cfgGlobalonly, bool, false, "Global only tracks")
  O2_DEFINE_CONFIGURABLE(cfgITSonly, bool, false, "ITS only tracks")
  O2_DEFINE_CONFIGURABLE(cfgFineBinning, bool, false, "Manually change to fine binning")

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
  ConfigurableAxis axisT0M{"axisT0M", {70, 0, 70000}, "N_{ch} (T0M)"};
  ConfigurableAxis axisFT0CAmp{"axisFT0CAmp", {50000, 0, 50000}, "axisFT0CAmp"};
  ConfigurableAxis axisFT0AAmp{"axisFT0AAmp", {50000, 0, 50000}, "axisFT0AAmp"};
  ConfigurableAxis axisFT0MAmp{"axisFT0MAmp", {50000, 0, 50000}, "axisFT0MAmp"};
  ConfigurableAxis axisNchPV{"axisNchPV", {4000, 0, 4000}, "N_{ch} (PV)"};
  ConfigurableAxis axisDCAz{"axisDCAz", {200, -2, 2}, "DCA_{z} (cm)"};
  ConfigurableAxis axisDCAxy{"axisDCAxy", {200, -1, 1}, "DCA_{xy} (cm)"};

  // Configurables for ZDC
  Configurable<int> nBinsAmp{"nBinsAmp", 1025, "nbinsAmp"};
  Configurable<float> maxZN{"maxZN", 4099.5, "Max ZN signal"};
  Configurable<float> maxZP{"maxZP", 3099.5, "Max ZP signal"};
  Configurable<float> maxZEM{"maxZEM", 3099.5, "Max ZEM signal"};
  Configurable<int> nBinsFit{"nBinsFit", 1000, "nbinsFit"};
  Configurable<float> maxMultFT0{"maxMultFT0", 5000, "Max FT0 signal"};

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
    kNOTIMEFRAMEBORDER,
    kNOITSROFRAMEBORDER,
    kNOPSAMEBUNCHPILEUP,
    kISGOODZVTXFT0VSPV,
    kISVERTEXITSTPC,
    kNOCOLLINTIMERANGESTANDART,
    kISGOODITSLAYERSALL,
    kAFTERMULTCUTS,
    kCENTRALITY,
    kNOOFEVENTSTEPS
  };

  enum CentEstimators {
    kCentFT0C,
    kCentFT0A,
    kCentFT0M,
    kCentFV0A,
    kCentFT0CVariant1,
    kNoCentEstimators
  };

  // Contruct Global+ITS sample
  static constexpr TrackSelectionFlags::flagtype TrackSelectionITS =
    TrackSelectionFlags::kITSNCls | TrackSelectionFlags::kITSChi2NDF |
    TrackSelectionFlags::kITSHits;
  static constexpr TrackSelectionFlags::flagtype TrackSelectionTPC =
    TrackSelectionFlags::kTPCNCls |
    TrackSelectionFlags::kTPCCrossedRowsOverNCls |
    TrackSelectionFlags::kTPCChi2NDF;
  static constexpr TrackSelectionFlags::flagtype TrackSelectionDCA =
    TrackSelectionFlags::kDCAz | TrackSelectionFlags::kDCAxy;
  static constexpr TrackSelectionFlags::flagtype TrackSelectionDCAXYonly =
    TrackSelectionFlags::kDCAxy;

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
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kNOTIMEFRAMEBORDER + 1, "kNoTimeFrameBorder");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kNOITSROFRAMEBORDER + 1, "kNoITSROFrameBorder");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kNOPSAMEBUNCHPILEUP + 1, "kNoSameBunchPileup");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kISGOODZVTXFT0VSPV + 1, "kIsGoodZvtxFT0vsPV");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kISVERTEXITSTPC + 1, "kIsVertexITSTPC");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kNOCOLLINTIMERANGESTANDART + 1, "kNoCollInTimeRangeStandard");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kISGOODITSLAYERSALL + 1, "kIsGoodITSLayersAll");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kAFTERMULTCUTS + 1, "After Mult cuts");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kCENTRALITY + 1, "Centrality");
    registry.add("hPhi", "#phi distribution", {HistType::kTH1D, {axisPhi}});
    registry.add("hPhiWeighted", "corrected #phi distribution", {HistType::kTH1D, {axisPhi}});
    registry.add("hEta", "", {HistType::kTH1D, {axisEta}});
    registry.add("hVtxZ", "Vexter Z distribution", {HistType::kTH1D, {axisVertex}});
    registry.add("hMult", "Multiplicity distribution", {HistType::kTH1D, {{3000, 0.5, 3000.5}}});
    registry.add("hCent", "Centrality distribution", {HistType::kTH1D, {{90, 0, 90}}});
    registry.add("cent_vs_Nch", ";Centrality (%); M (|#eta| < 0.8);", {HistType::kTH2D, {axisCentrality, axisNch}});

    // Centrality estimators
    registry.add("hCentEstimators", "Number of Unfiltered Events;; No. of Events", {HistType::kTH1D, {{kNoCentEstimators, -0.5, static_cast<int>(kNoCentEstimators) - 0.5}}});
    registry.get<TH1>(HIST("hCentEstimators"))->GetXaxis()->SetBinLabel(kCentFT0C + 1, "FT0C");
    registry.get<TH1>(HIST("hCentEstimators"))->GetXaxis()->SetBinLabel(kCentFT0A + 1, "FT0A");
    registry.get<TH1>(HIST("hCentEstimators"))->GetXaxis()->SetBinLabel(kCentFT0M + 1, "FT0M");
    registry.get<TH1>(HIST("hCentEstimators"))->GetXaxis()->SetBinLabel(kCentFV0A + 1, "FV0A");
    registry.get<TH1>(HIST("hCentEstimators"))->GetXaxis()->SetBinLabel(kCentFT0CVariant1 + 1, "FT0CVar1");
    registry.add("hCentFT0C", "Uncorrected FT0C;Centrality FT0C ;Events", kTH1F, {axisCentrality});
    registry.add("hCentFT0A", "Uncorrected FT0A;Centrality FT0A ;Events", kTH1F, {axisCentrality});
    registry.add("hCentFT0M", "Uncorrected FT0M;Centrality FT0M ;Events", kTH1F, {axisCentrality});
    registry.add("hCentFV0A", "Uncorrected FV0A;Centrality FV0A ;Events", kTH1F, {axisCentrality});
    registry.add("hCentFT0CVariant1", "Uncorrected FT0CVariant1;Centrality FT0CVariant1 ;Events", kTH1F, {axisCentrality});

    // Before cuts
    registry.add("BeforeCut_globalTracks_centT0C", "before cut;Centrality T0C;mulplicity global tracks", {HistType::kTH2D, {axisCentForQA, axisNch}});
    registry.add("BeforeCut_PVTracks_centT0C", "before cut;Centrality T0C;mulplicity PV tracks", {HistType::kTH2D, {axisCentForQA, axisNchPV}});
    registry.add("BeforeCut_globalTracks_PVTracks", "before cut;mulplicity PV tracks;mulplicity global tracks", {HistType::kTH2D, {axisNchPV, axisNch}});
    registry.add("BeforeCut_globalTracks_multT0A", "before cut;mulplicity T0A;mulplicity global tracks", {HistType::kTH2D, {axisT0A, axisNch}});
    registry.add("BeforeCut_globalTracks_multV0A", "before cut;mulplicity V0A;mulplicity global tracks", {HistType::kTH2D, {axisT0A, axisNch}});
    registry.add("BeforeCut_multV0A_multT0A", "before cut;mulplicity T0A;mulplicity V0A", {HistType::kTH2D, {axisT0A, axisT0A}});
    registry.add("BeforeCut_multT0C_centT0C", "before cut;Centrality T0C;mulplicity T0C", {HistType::kTH2D, {axisCentForQA, axisT0C}});
    registry.add("BeforeCut_multT0A_centT0A", "before cut;Centrality T0C;mulplicity T0A", {HistType::kTH2D, {axisCentForQA, axisT0A}});
    registry.add("BeforeCut_multFT0M_centFT0M", "before cut;Centrality FT0M;mulplicity FT0M", {HistType::kTH2D, {axisCentForQA, axisT0M}});

    // After cuts
    registry.add("globalTracks_centT0C_Aft", "after cut;Centrality T0C;mulplicity global tracks", {HistType::kTH2D, {axisCentForQA, axisNch}});
    registry.add("PVTracks_centT0C_Aft", "after cut;Centrality T0C;mulplicity PV tracks", {HistType::kTH2D, {axisCentForQA, axisNchPV}});
    registry.add("globalTracks_PVTracks_Aft", "after cut;mulplicity PV tracks;mulplicity global tracks", {HistType::kTH2D, {axisNchPV, axisNch}});
    registry.add("globalTracks_multT0A_Aft", "after cut;mulplicity T0A;mulplicity global tracks", {HistType::kTH2D, {axisT0A, axisNch}});
    registry.add("globalTracks_multV0A_Aft", "after cut;mulplicity V0A;mulplicity global tracks", {HistType::kTH2D, {axisT0A, axisNch}});
    registry.add("multV0A_multT0A_Aft", "after cut;mulplicity T0A;mulplicity V0A", {HistType::kTH2D, {axisT0A, axisT0A}});
    registry.add("multT0C_centT0C_Aft", "after cut;Centrality T0C;mulplicity T0C", {HistType::kTH2D, {axisCentForQA, axisT0C}});
    registry.add("multT0A_centT0A_Aft", "after cut;Centrality T0A;mulplicity T0A", {HistType::kTH2D, {axisCentForQA, axisT0A}});
    registry.add("multFT0M_centFT0M_Aft", "after cut;Centrality FT0M;mulplicity FT0M", {HistType::kTH2D, {axisCentForQA, axisT0M}});

    // FT0 plots
    registry.add("FT0CAmp", ";FT0C amplitude;Events", kTH1F, {axisFT0CAmp});
    registry.add("FT0AAmp", ";FT0A amplitude;Events", kTH1F, {axisFT0AAmp});
    registry.add("FT0MAmp", ";FT0M amplitude;Events", kTH1F, {axisFT0MAmp});

    // ZDC plots
    const AxisSpec axisEvent{3, 0., +3.0, ""};
    registry.add("hEventCounterForZDC", "Event counter", kTH1F, {axisEvent});
    registry.add("ZNAcoll", "ZNAcoll; ZNA amplitude; Entries", {HistType::kTH1F, {{nBinsAmp, -0.5, maxZN}}});
    registry.add("ZPAcoll", "ZPAcoll; ZPA amplitude; Entries", {HistType::kTH1F, {{nBinsAmp, -0.5, maxZP}}});
    registry.add("ZNCcoll", "ZNCcoll; ZNC amplitude; Entries", {HistType::kTH1F, {{nBinsAmp, -0.5, maxZN}}});
    registry.add("ZPCcoll", "ZPCcoll; ZPC amplitude; Entries", {HistType::kTH1F, {{nBinsAmp, -0.5, maxZP}}});
    registry.add("ZNvsFT0correl", "ZNvsFT0correl; FT0 amplitude; ZN", {HistType::kTH2F, {{{nBinsFit, 0., maxMultFT0}, {nBinsAmp, -0.5, 2. * maxZN}}}});
    registry.add("ZDCAmp", "ZDC Amplitude; ZDC Amplitude; Events", {HistType::kTH1F, {{nBinsAmp, -0.5, maxZP}}});
    registry.add("ZNAmp", "ZNA+ZNC Amplitude; ZN Amplitude; Events", {HistType::kTH1F, {{nBinsAmp, -0.5, maxZN}}});
    registry.add("ZPAmp", "ZPA+ZPC Amplitude; ZP Amplitude; Events", {HistType::kTH1F, {{nBinsAmp, -0.5, maxZP}}});
    registry.add("ZNvsZEMcoll", "ZNvsZEMcoll; ZEM; ZDC energy (GeV)", {HistType::kTH2F, {{{nBinsAmp, -0.5, maxZEM}, {nBinsAmp, -0.5, 2. * maxZN}}}});
    registry.add("ZNvsZEMcoll05", "ZNvsZEMcoll; ZEM; ZDC energy (GeV)", {HistType::kTH2F, {{{nBinsAmp, -0.5, maxZEM}, {nBinsAmp, -0.5, 2. * maxZN}}}});
    registry.add("ZNvsZEMcoll510", "ZNvsZEMcoll; ZEM; ZDC energy (GeV)", {HistType::kTH2F, {{{nBinsAmp, -0.5, maxZEM}, {nBinsAmp, -0.5, 2. * maxZN}}}});
    registry.add("ZNvsZEMcoll1020", "ZNvsZEMcoll; ZEM; ZDC energy (GeV)", {HistType::kTH2F, {{{nBinsAmp, -0.5, maxZEM}, {nBinsAmp, -0.5, 2. * maxZN}}}});
    registry.add("ZNvsZEMcoll2030", "ZNvsZEMcoll; ZEM; ZDC energy (GeV)", {HistType::kTH2F, {{{nBinsAmp, -0.5, maxZEM}, {nBinsAmp, -0.5, 2. * maxZN}}}});
    registry.add("ZNvsZEMcollrest", "ZNvsZEMcoll; ZEM; ZDC energy (GeV)", {HistType::kTH2F, {{{nBinsAmp, -0.5, maxZEM}, {nBinsAmp, -0.5, 2. * maxZN}}}});

    // Track plots
    registry.add("Nch", "N_{ch 0-5%} vs #Events;N_{ch 0-5%};No. of Events", {HistType::kTH1D, {axisNch}});
    registry.add("Events_per_Centrality_Bin", "Events_per_Centrality_Bin;Centrality FT0C;No. of Events", kTH1F, {axisCentrality});
    registry.add("Tracks_per_Centrality_Bin", "Tracks_per_Centrality_Bin;Centrality FT0C;No. of Tracks", kTH1F, {axisCentrality});
    registry.add("pt_Cen_GlobalOnly", "pt_Cen_Global;Centrality (%); p_{T} (GeV/c);", {HistType::kTH2D, {axisCentrality, axisPt}});
    registry.add("phi_Cen_GlobalOnly", "phi_Cen_Global;Centrality (%); #phi;", {HistType::kTH2D, {axisCentrality, axisPhi}});
    registry.add("pt_Cen_ITSOnly", "pt_Cen_ITS;Centrality (%); p_{T} (GeV/c);", {HistType::kTH2D, {axisCentrality, axisPt}});
    registry.add("phi_Cen_ITSOnly", "phi_Cen_ITS;Centrality (%); #phi;", {HistType::kTH2D, {axisCentrality, axisPhi}});

    // Track types
    registry.add("GlobalplusITS", "Global plus ITS;Centrality FT0C;Nch", {HistType::kTH2D, {axisCentrality, axisNch}});
    registry.add("Globalonly", "Global only;Centrality FT0C;Nch", {HistType::kTH2D, {axisCentrality, axisNch}});
    registry.add("ITSonly", "ITS only;Centrality FT0C;Nch", {HistType::kTH2D, {axisCentrality, axisNch}});

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
      fWeights->setPtBins(nPtBins, ptBins);
      fWeights->init(true, false);
    }

    // add in FlowContainer to Get boostrap sample automatically -- Use post process flow task
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
      weight_nua = mAcceptance->getNUA(phi, eta, vtxz);
    else
      weight_nua = 1;
    return true;
  }

  template <typename TCollision>
  bool eventSelected(o2::aod::mult::MultNTracksPV, TCollision collision, const int multTrk, const float centrality)
  {
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
    if (cfgIsVertexITSTPC) {
      if (!collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
        // removes collisions without vertex match between ITS-TPC
        return false;
      }
      registry.fill(HIST("hEventCount"), kISVERTEXITSTPC);
      if (cfgNoCollInTimeRangeStandard) {
        if (!collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
          // no collisions in specified time range
          return false;
        }
        registry.fill(HIST("hEventCount"), kNOCOLLINTIMERANGESTANDART);
      }
      if (cfgEvSelkIsGoodITSLayersAll) {
        if (cfgEvSelkIsGoodITSLayersAll && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
          // removes dead staves of ITS
          return false;
        }
        registry.fill(HIST("hEventCount"), kISGOODITSLAYERSALL);
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
      if (cfgV0AT0A5Sigma)
        if (std::abs(collision.multFV0A() - fT0AV0AMean->Eval(collision.multFT0A())) > 5 * fT0AV0ASigma->Eval(collision.multFT0A()))
          return false;
    }

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
    if (cfgDCAzPt && (std::fabs(track.dcaZ()) > (0.004f + 0.013f / track.pt())))
      return false;

    if (cfgTrackSel) {
      return myTrackSel.IsSelected(track);
    } else if (cfgGlobalplusITS) {
      return ((track.tpcNClsFound() >= cfgCutTPCclu) && (track.itsNCls() >= cfgCutITSclu));
    } else if (cfgGlobalonly) {
      return ((track.tpcNClsFound() >= cfgCutTPCclu));
    } else if (cfgITSonly) {
      return ((track.itsNCls() >= cfgCutITSclu));
    } else {
      return false;
    }
  }

  // Apply process filters
  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex && (aod::cent::centFT0C > cfgMinCentFT0C) && (aod::cent::centFT0C < cfgMaxCentFT0C);
  Filter trackFilter = ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::ITS) &&
                       ncheckbit(aod::track::trackCutFlag, TrackSelectionITS) &&
                       ifnode(ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::TPC),
                              ncheckbit(aod::track::trackCutFlag, TrackSelectionTPC), true) &&
                       ifnode(dcaZ > 0.f, nabs(aod::track::dcaZ) <= dcaZ && ncheckbit(aod::track::trackCutFlag, TrackSelectionDCAXYonly),
                              ncheckbit(aod::track::trackCutFlag, TrackSelectionDCA)) &&
                       (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtMin) && (aod::track::pt < cfgCutPtMax) && (aod::track::tpcChi2NCl < cfgCutChi2prTPCcls);

  using Colls = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs, aod::CentFT0As, aod::CentFT0Ms, aod::CentFV0As, aod::CentFT0CVariant1s>>; // collisions filter
  using AodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksDCA, aod::TracksExtra>>; // tracks filter

  using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;

  void processData(Colls::iterator const& collision, aod::BCsWithTimestamps const&, AodTracks const& tracks, aod::FT0s const&, aod::Zdcs const&, BCsRun3 const&)
  {
    registry.fill(HIST("hEventCount"), kFILTERED);
    if (!collision.sel8())
      return;

    int nTotal = tracks.size();
    if (nTotal < 1)
      return;

    // Choose centrality estimator -- Only one can be true
    auto centrality = -1;
    if (cfgcentEstFt0c) {
      centrality = collision.centFT0C();
      registry.fill(HIST("hCentEstimators"), kCentFT0C);
      registry.fill(HIST("hCentFT0C"), centrality);
    }
    if (cfgcentEstFt0a) {
      centrality = collision.centFT0A();
      registry.fill(HIST("hCentEstimators"), kCentFT0A);
      registry.fill(HIST("hCentFT0A"), centrality);
    }
    if (cfgcentEstFt0m) {
      centrality = collision.centFT0M();
      registry.fill(HIST("hCentEstimators"), kCentFT0M);
      registry.fill(HIST("hCentFT0M"), centrality);
    }
    if (cfgcentEstFv0a) {
      centrality = collision.centFV0A();
      registry.fill(HIST("hCentEstimators"), kCentFV0A);
      registry.fill(HIST("hCentFV0A"), centrality);
    }
    if (cfgcentEstFt0cVariant1) {
      centrality = collision.centFT0CVariant1();
      registry.fill(HIST("hCentEstimators"), kCentFT0CVariant1);
      registry.fill(HIST("hCentFT0CVariant1"), centrality);
    }

    // fill event QA before cuts
    registry.fill(HIST("BeforeCut_globalTracks_centT0C"), collision.centFT0C(), tracks.size());
    registry.fill(HIST("BeforeCut_PVTracks_centT0C"), collision.centFT0C(), collision.multNTracksPV());
    registry.fill(HIST("BeforeCut_globalTracks_PVTracks"), collision.multNTracksPV(), tracks.size());
    registry.fill(HIST("BeforeCut_globalTracks_multT0A"), collision.multFT0A(), tracks.size());
    registry.fill(HIST("BeforeCut_globalTracks_multV0A"), collision.multFV0A(), tracks.size());
    registry.fill(HIST("BeforeCut_multV0A_multT0A"), collision.multFT0A(), collision.multFV0A());
    registry.fill(HIST("BeforeCut_multT0C_centT0C"), collision.centFT0C(), collision.multFT0C());
    registry.fill(HIST("BeforeCut_multT0A_centT0A"), collision.centFT0A(), collision.multFT0A());
    registry.fill(HIST("BeforeCut_multFT0M_centFT0M"), collision.centFT0M(), collision.multFT0M());
    registry.fill(HIST("hEventCount"), kSEL8);

    if (cfgOccupancy) {
      int occupancy = collision.trackOccupancyInTimeRange();
      if (occupancy < cfgCutOccupancyLow || occupancy > cfgCutOccupancyHigh)
        return;
      registry.fill(HIST("hEventCount"), kOCCUPANCY);
    }

    if (cfgUseAdditionalEventCut && !eventSelected(o2::aod::mult::MultNTracksPV(), collision, tracks.size(), centrality)) {
      return;
    }

    const auto& foundBC = collision.foundBC_as<BCsRun3>();
    if (foundBC.has_zdc()) {
      registry.fill(HIST("hEventCounterForZDC"), 1);

      // FT0 amplitude to use in fine binning
      double ft0aAmp = 0;
      double ft0cAmp = 0;
      double ft0mAmp = 0;

      if (foundBC.has_ft0()) {
        for (const auto& amplitude : foundBC.ft0().amplitudeA()) {
          ft0aAmp += amplitude;
        }
        for (const auto& amplitude : foundBC.ft0().amplitudeC()) {
          ft0cAmp += amplitude;
        }
      } else {
        ft0aAmp = ft0cAmp = -999;
      }

      registry.fill(HIST("FT0AAmp"), ft0aAmp);
      registry.fill(HIST("FT0CAmp"), ft0cAmp);

      ft0mAmp = ft0aAmp + ft0cAmp;
      registry.fill(HIST("FT0MAmp"), ft0mAmp);

      // ZDC amplitude to use in fine binning
      const auto& zdcread = foundBC.zdc();
      auto aZNA = zdcread.amplitudeZNA();
      auto aZNC = zdcread.amplitudeZNC();
      auto aZPA = zdcread.amplitudeZPA();
      auto aZPC = zdcread.amplitudeZPC();
      auto aZEM1 = zdcread.amplitudeZEM1();
      auto aZEM2 = zdcread.amplitudeZEM2();

      registry.fill(HIST("ZNAcoll"), aZNA);
      registry.fill(HIST("ZNCcoll"), aZNC);
      registry.fill(HIST("ZPAcoll"), aZPA);
      registry.fill(HIST("ZPCcoll"), aZPC);

      registry.fill(HIST("ZNvsFT0correl"), (ft0aAmp + ft0cAmp) / 100., aZNC + aZNA);

      double aZDC = aZNC + aZNA + aZPA + aZPC;
      registry.fill(HIST("ZDCAmp"), aZDC);
      registry.fill(HIST("ZNAmp"), aZNC + aZNA);
      registry.fill(HIST("ZPAmp"), aZPA + aZPC);

      registry.fill(HIST("ZNvsZEMcoll"), aZEM1 + aZEM2, aZNA + aZNC);

      if (centrality >= 0 && centrality <= 5) {
        registry.fill(HIST("ZNvsZEMcoll05"), aZEM1 + aZEM2, aZNA + aZNC);
      } else if (centrality > 5 && centrality <= 10) {
        registry.fill(HIST("ZNvsZEMcoll510"), aZEM1 + aZEM2, aZNA + aZNC);
      } else if (centrality > 10 && centrality <= 20) {
        registry.fill(HIST("ZNvsZEMcoll1020"), aZEM1 + aZEM2, aZNA + aZNC);
      } else if (centrality > 20 && centrality <= 30) {
        registry.fill(HIST("ZNvsZEMcoll2030"), aZEM1 + aZEM2, aZNA + aZNC);
      } else {
        registry.fill(HIST("ZNvsZEMcollrest"), aZEM1 + aZEM2, aZNA + aZNC);
      }
    } // End of ZDC

    float vtxz = collision.posZ();
    float lRandom = fRndm->Rndm();
    registry.fill(HIST("hVtxZ"), vtxz);
    registry.fill(HIST("hMult"), nTotal);
    registry.fill(HIST("hCent"), centrality);
    registry.fill(HIST("cent_vs_Nch"), centrality, nTotal);

    fGFW->Clear();

    auto bc = collision.bc_as<BCsRun3>();
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
    registry.fill(HIST("multT0A_centT0A_Aft"), collision.centFT0A(), collision.multFT0A());
    registry.fill(HIST("multFT0M_centFT0M_Aft"), collision.centFT0M(), collision.multFT0M());

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
        fWeights->fill(track.phi(), track.eta(), vtxz, track.pt(), centrality, 0);

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

      if (cfgGlobalplusITS) {
        if (withinPtRef) {
          globalTracksNch++;
          registry.fill(HIST("GlobalplusITS"), centrality, globalTracksNch);
          fGFW->Fill(track.eta(), fPtAxis->FindBin(track.pt()) - 1, track.phi(), wacc * weff, 1);
        }
      }

      if (track.hasTPC()) {
        if (cfgGlobalonly) {
          if (withinPtRef) {
            globalTracksNch++;
            registry.fill(HIST("Globalonly"), centrality, globalTracksNch);
            registry.fill(HIST("pt_Cen_GlobalOnly"), centrality, track.pt());
            registry.fill(HIST("phi_Cen_GlobalOnly"), centrality, track.phi());
            fGFW->Fill(track.eta(), fPtAxis->FindBin(track.pt()) - 1, track.phi(), wacc * weff, 1);
          }
        }
      } else {
        if (cfgITSonly) {
          if (withinPtRef) {
            globalTracksNch++;
            registry.fill(HIST("ITSonly"), centrality, globalTracksNch);
            registry.fill(HIST("pt_Cen_ITSOnly"), centrality, track.pt());
            registry.fill(HIST("phi_Cen_ITSOnly"), centrality, track.phi());
            fGFW->Fill(track.eta(), fPtAxis->FindBin(track.pt()) - 1, track.phi(), wacc * weff, 1);
          }
        }
      }

      if (cfgFineBinning)
        fGFW->Fill(track.eta(), 1, track.phi(), wacc * weff, 1);

    } // End of track loop

    // Only one type of track will be plotted
    registry.fill(HIST("Events_per_Centrality_Bin"), centrality);
    registry.fill(HIST("Tracks_per_Centrality_Bin"), centrality, globalTracksNch);

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

      registry.fill(HIST("Nch"), globalTracksNch);
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
