// Copyright 2020-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file flattenictyPikp.cxx
/// \author Gyula Bencedi, gyula.bencedi@cern.ch
/// \brief Task to produce pion, kaon, proton high-pT
///        distributions as a function of charged-particle flattenicity
/// \since 26 June 2025

#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/mcCentrality.h"
#include "PWGLF/Utils/inelGt.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DataFormatsFIT/Triggers.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/StaticFor.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/PID.h>

#include <TEfficiency.h>
#include <TF1.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH3.h>
#include <THashList.h>
#include <THnSparse.h>
#include <TProfile2D.h>
#include <TString.h>

#include <fmt/format.h>

#include <algorithm>
#include <array>
#include <bitset>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iterator>
#include <map>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

using std::string;
using namespace o2;
using namespace o2::framework;
using namespace o2::constants::math;
using namespace o2::aod::rctsel;

auto static constexpr CminCharge = 3.f;
static constexpr float CnullInt = 0;
static constexpr float Cnull = 0.0f;
static constexpr float ConeInt = 1;
static constexpr float Cone = 1.0f;

// FV0 specific constants
static constexpr int CmaxRingsFV0 = 5;
static constexpr int CnCellsFV0 = 48;
static constexpr int CinnerFV0 = 32;
static constexpr float Cfv0IndexPhi[5] = {0., 8., 16., 24., 32.};
static constexpr float CmaxEtaFV0 = 5.1;
static constexpr float CminEtaFV0 = 2.2;
static constexpr float CdEtaFV0 = (CmaxEtaFV0 - CminEtaFV0) / CmaxRingsFV0;
auto static constexpr CminAccFT0A = 3.5f;
auto static constexpr CmaxAccFT0A = 4.9f;
auto static constexpr CminAccFT0C = -3.3f;
auto static constexpr CmaxAccFT0C = -2.1f;
// PID names
static constexpr int CprocessIdWeak = 4;
static constexpr int Ncharges = 2;
static constexpr o2::track::PID::ID Npart = 5;
static constexpr o2::track::PID::ID NpartChrg = Npart * Ncharges;
static constexpr int PDGs[] = {11, 13, 211, 321, 2212};
static constexpr int PidSgn[NpartChrg] = {11, 13, 211, 321, 2212, -11, -13, -211, -321, -2212};
static constexpr const char* Pid[Npart] = {"el", "mu", "pi", "ka", "pr"};
static constexpr const char* PidChrg[NpartChrg] = {"e^{-}", "#mu^{-}", "#pi^{+}", "K^{+}", "p", "e^{+}", "#mu^{+}", "#pi^{-}", "K^{-}", "#bar{p}"};
static constexpr std::string_view CspeciesAll[Npart] = {"El", "Mu", "Pi", "Ka", "Pr"};
// histogram naming
static constexpr std::string_view PidDir[] = {"el/", "pi/", "pr/"};
static constexpr std::string_view V0Dir[] = {"Ga/", "K0s/", "La/", "ALa/"};
static constexpr std::string_view Ccharge[] = {"all/", "pos/", "neg/"};
static constexpr std::string_view Cprefix = "Tracks/";
static constexpr std::string_view CprefixCleanTof = "Tracks/CleanTof/";
static constexpr std::string_view CprefixCleanV0 = "Tracks/CleanV0/";
static constexpr std::string_view CprefixV0qa = "Tracks/V0qa/";
static constexpr std::string_view Cstatus[] = {"preSel/", "postSel/"};
static constexpr std::string_view CstatCalib[] = {"preCalib/", "postCalib/"};
static constexpr std::string_view CdEdxMcRecPrim = "/hdEdxMcRecPrim";
static constexpr std::string_view CdEdxMcRecPrimF = "Tracks/{}/hdEdxMcRecPrim";
static constexpr std::string_view CdEdxMcRecPrimSel = "/hdEdxMcRecPrimSel";
static constexpr std::string_view CdEdxMcRecPrimSelF = "Tracks/{}/hdEdxMcRecPrimSel";
static constexpr std::string_view CpTvsDCAxy = "/hPtVsDCAxy";
static constexpr std::string_view CpTvsDCAxyF = "Tracks/{}/hPtVsDCAxy";
static constexpr std::string_view CpTvsDCAxyAll = "/hPtVsDCAxyAll";
static constexpr std::string_view CpTvsDCAxyAllF = "Tracks/{}/hPtVsDCAxyAll";
static constexpr std::string_view CpTvsDCAxyPrimAll = "/hPtVsDCAxyRecPrimAll";
static constexpr std::string_view CpTvsDCAxyPrimAllF = "Tracks/{}/hPtVsDCAxyRecPrimAll";
static constexpr std::string_view CpTvsDCAxyWeakAll = "/hPtVsDCAxyRecWeakAll";
static constexpr std::string_view CpTvsDCAxyWeakAllF = "Tracks/{}/hPtVsDCAxyRecWeakAll";
static constexpr std::string_view CpTvsDCAxyMatAll = "/hPtVsDCAxyRecMatAll";
static constexpr std::string_view CpTvsDCAxyMatAllF = "Tracks/{}/hPtVsDCAxyRecMatAll";
static constexpr std::string_view CpTgenPrimSgn = "/hPtGenPrimSgn";
static constexpr std::string_view CpTgenPrimSgnF = "Tracks/{}/hPtGenPrimSgn";
static constexpr std::string_view CpTrecCollPrimSgn = "/hPtRecCollPrimSgn";
static constexpr std::string_view CpTrecCollPrimSgnF = "Tracks/{}/hPtRecCollPrimSgn";
static constexpr std::string_view CpTmcClosureGenPrim = "/hPtMCclosureGenPrim";
static constexpr std::string_view CpTmcClosureGenPrimF = "Tracks/{}/hPtMCclosureGenPrim";
static constexpr std::string_view CpTmcClosureRec = "/hPtMCclosureRec";
static constexpr std::string_view CpTmcClosureRecF = "Tracks/{}/hPtMCclosureRec";
static constexpr std::string_view CpTeffPrimRecEvt = "/hPtEffPrimRecEvt";
static constexpr std::string_view CpTeffPrimRecEvtF = "Tracks/{}/hPtEffPrimRecEvt";
static constexpr std::string_view CpTeffGenPrimRecEvt = "/hPtEffGenPrimRecEvt";
static constexpr std::string_view CpTeffGenPrimRecEvtF = "Tracks/{}/hPtEffGenPrimRecEvt";

enum PidType {
  kEl = 0,
  kPi,
  kPr
};

enum V0sSel {
  kNaN = -1,
  kGa = 0,
  kKz = 1,
  kLam = 2,
  kaLam = 3
};

enum FillType {
  kBefore,
  kAfter
};

enum ChargeType {
  kAll,
  kPos,
  kNeg
};

enum TrkSel {
  trkSelAll,
  trkSelEta,
  trkSelPt,
  trkSelDCA,
  trkNRowsTPC,
  trkSelNClsFound,
  trkSelNClsPID,
  trkSelTPCBndr,
  nTrkSel
};

enum V0Sel {
  v0SelAll,
  v0SelRejectSameSign,
  v0SelRejectV0sAtTPCSector,
  v0SelCosPA,
  v0SelV0radius,
  v0SelDCAposToPV,
  v0SelDaughters,
  v0SelDCAv0daughter,
  nV0Sel
};

enum EvtSel {
  evtSelAll,
  evtSelSel8,
  evtSelNoITSROFrameBorder,
  evtSelkNoTimeFrameBorder,
  evtSelkNoSameBunchPileup,
  evtSelkIsGoodZvtxFT0vsPV,
  evtSelkIsVertexITSTPC,
  evtSelkIsVertexTOFmatched,
  evtSelVtxZ,
  evtSelINELgt0,
  evtSelRCTFlagChecker,
  nEvtSel
};

struct MultE {
  static constexpr int CnoMult = 0;
  static constexpr int CmultFT0M = 1;
  static constexpr int CmultTPC = 2;
};

std::array<float, CnCellsFV0> rhoLatticeFV0{0};
std::array<float, CnCellsFV0> fv0AmplitudeWoCalib{0};
std::array<std::shared_ptr<THnSparse>, NpartChrg> hPtGenRecEvt{};
std::array<std::shared_ptr<THnSparse>, NpartChrg> hPtGenPrimRecEvt{};
std::array<std::shared_ptr<THnSparse>, NpartChrg> hPtGenRecEvtGtZero{};
std::array<std::shared_ptr<THnSparse>, NpartChrg> hPtGenPrimRecEvtGtZero{};
std::array<std::shared_ptr<TH1>, NpartChrg> hPtEffRec{};
std::array<std::shared_ptr<TH1>, NpartChrg> hPtEffGen{};
std::array<std::shared_ptr<THnSparse>, NpartChrg> hPtEffRecGoodCollPrim{};
std::array<std::shared_ptr<THnSparse>, NpartChrg> hPtEffRecGoodCollWeak{};
std::array<std::shared_ptr<THnSparse>, NpartChrg> hPtEffRecGoodCollMat{};
std::array<std::shared_ptr<THnSparse>, NpartChrg> hPtEffGenPrim{};
std::array<std::shared_ptr<THnSparse>, NpartChrg> hPtEffGenWeak{};
std::array<std::shared_ptr<THnSparse>, NpartChrg> hPtEffGenMat{};
std::array<std::shared_ptr<THnSparse>, NpartChrg> hPtEffGenPrimEvtSelGen{};
std::array<std::shared_ptr<THnSparse>, NpartChrg> hPtEffGenWeakEvtSelGen{};
std::array<std::shared_ptr<THnSparse>, NpartChrg> hPtEffGenMatEvtSelGen{};

std::array<std::shared_ptr<TH2>, NpartChrg> hDCAxyRecBadCollPrim{};
std::array<std::shared_ptr<TH2>, NpartChrg> hDCAxyRecBadCollWeak{};
std::array<std::shared_ptr<TH2>, NpartChrg> hDCAxyRecBadCollMat{};
std::array<std::shared_ptr<TH2>, NpartChrg> hPtVsDCAxyRecGoodCollPrim{};
std::array<std::shared_ptr<TH2>, NpartChrg> hPtVsDCAxyRecGoodCollWeak{};
std::array<std::shared_ptr<TH2>, NpartChrg> hPtVsDCAxyRecGoodCollMat{};

struct FlattenictyPikp {

  HistogramRegistry registryData{"registryData", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryMC{"registryMC", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryQC{"registryQC", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  OutputObj<THashList> listEfficiency{"Efficiency"};
  Service<o2::framework::O2DatabasePDG> pdg;

  std::vector<float> fv0AmplCorr{};
  TProfile2D* zVtxMap = nullptr;
  float magField;
  int runNumber{-1};
  o2::parameters::GRPMagField* grpmag = nullptr;

  Configurable<int> multEst{"multEst", 1, "0: without multiplicity; 1: MultFT0M; 2: MultTPC"};
  Configurable<bool> applyCalibGain{"applyCalibGain", false, "equalize detector amplitudes"};
  Configurable<bool> applyCalibVtx{"applyCalibVtx", false, "equalize Amp vs vtx"};
  Configurable<bool> applyCalibDeDx{"applyCalibDeDx", false, "calibration of dedx signal"};
  Configurable<bool> applyCalibDeDxFromCCDB{"applyCalibDeDxFromCCDB", false, "use CCDB-based calibration of dedx signal"};
  Configurable<bool> cfgFillTrackQaHist{"cfgFillTrackQaHist", false, "fill track QA histograms"};
  Configurable<bool> cfgFillNclVsPhiCutQaHist{"cfgFillNclVsPhiCutQaHist", false, "fill TPC cluster vs geometrical cut QA histograms"};
  Configurable<bool> cfgFilldEdxCalibHist{"cfgFilldEdxCalibHist", false, "fill dEdx calibration histograms"};
  Configurable<bool> cfgFilldEdxQaHist{"cfgFilldEdxQaHist", false, "fill dEdx QA histograms"};
  Configurable<bool> cfgFillDCAxyHist{"cfgFillDCAxyHist", false, "fill nsigma QA histograms"};
  Configurable<bool> cfgFillV0Hist{"cfgFillV0Hist", false, "fill V0 histograms"};
  Configurable<bool> cfgFillChrgType{"cfgFillChrgType", false, "fill histograms per charge types"};
  Configurable<bool> cfgFillChrgTypeV0s{"cfgFillChrgTypeV0s", false, "fill V0s histograms per charge types"};
  Configurable<std::string> cfgCalibDeDxFunction{"cfgCalibDeDxFunction", "pol8", "Functional form for dEdx calibration"};
  Configurable<std::vector<float>> paramsFuncMIPposEtaP{"paramsFuncMIPposEtaP", std::vector<float>{-1.f}, "function parameters"};
  Configurable<std::vector<float>> paramsFuncMIPnegEtaP{"paramsFuncMIPnegEtaP", std::vector<float>{-1.f}, "function parameters"};
  Configurable<std::vector<float>> paramsFuncMIPallEtaP{"paramsFuncMIPallEtaP", std::vector<float>{-1.f}, "function parameters"};
  Configurable<std::vector<float>> paramsFuncMIPposEtaN{"paramsFuncMIPposEtaN", std::vector<float>{-1.f}, "function parameters"};
  Configurable<std::vector<float>> paramsFuncMIPnegEtaN{"paramsFuncMIPnegEtaN", std::vector<float>{-1.f}, "function parameters"};
  Configurable<std::vector<float>> paramsFuncMIPallEtaN{"paramsFuncMIPallEtaN", std::vector<float>{-1.f}, "function parameters"};
  Configurable<std::vector<float>> paramsFuncPlateaUposEtaP{"paramsFuncPlateaUposEtaP", std::vector<float>{-1.f}, "function parameters"};
  Configurable<std::vector<float>> paramsFuncPlateaUnegEtaP{"paramsFuncPlateaUnegEtaP", std::vector<float>{-1.f}, "function parameters"};
  Configurable<std::vector<float>> paramsFuncPlateaUallEtaP{"paramsFuncPlateaUallEtaP", std::vector<float>{-1.f}, "function parameters"};
  Configurable<std::vector<float>> paramsFuncPlateaUposEtaN{"paramsFuncPlateaUposEtaN", std::vector<float>{-1.f}, "function parameters"};
  Configurable<std::vector<float>> paramsFuncPlateaUnegEtaN{"paramsFuncPlateaUnegEtaN", std::vector<float>{-1.f}, "function parameters"};
  Configurable<std::vector<float>> paramsFuncPlateaUallEtaN{"paramsFuncPlateaUallEtaN", std::vector<float>{-1.f}, "function parameters"};
  Configurable<std::string> cfgGainEqCcdbPath{"cfgGainEqCcdbPath", "Users/g/gbencedi/flattenicity/GainEq", "CCDB path for gain equalization constants"};
  Configurable<std::string> cfgVtxEqCcdbPath{"cfgVtxEqCcdbPath", "Users/g/gbencedi/flattenicity/ZvtxEq", "CCDB path for z-vertex equalization constants"};
  Configurable<std::string> cfgDeDxCalibCcdbPath{"cfgDeDxCalibCcdbPath", "Users/g/gbencedi/flattenicity/dEdxCalib", "CCDB path for dEdx calibration"};
  Configurable<bool> cfgUseCcdbForRun{"cfgUseCcdbForRun", true, "Get ccdb object based on run number instead of timestamp"};
  Configurable<bool> cfgStoreThnSparse{"cfgStoreThnSparse", false, "Store histograms as THnSparse"};

  struct : ConfigurableGroup {
    Configurable<bool> cfgCustomTVX{"cfgCustomTVX", false, "Ask for custom TVX instead of sel8"};
    Configurable<bool> cfgRemoveNoTimeFrameBorder{"cfgRemoveNoTimeFrameBorder", false, "Bunch crossing is far from Time Frame borders"};
    Configurable<bool> cfgRemoveITSROFrameBorder{"cfgRemoveITSROFrameBorder", false, "Bunch crossing is far from ITS RO Frame border"};
    Configurable<float> cfgCutVtxZ{"cfgCutVtxZ", 10.0f, "Accepted z-vertex range"};
    Configurable<bool> useZVtxCutMC{"useZVtxCutMC", true, "use Zvtx cut in MC"};
    Configurable<bool> useINELCutMC{"useINELCutMC", true, "use INEL>0 cut in MC"};
    Configurable<bool> cfgRemoveNoSameBunchPileup{"cfgRemoveNoSameBunchPileup", true, "Reject collisions in case of pileup with another collision in the same foundBC"};
    Configurable<bool> cfgRequireIsGoodZvtxFT0vsPV{"cfgRequireIsGoodZvtxFT0vsPV", true, "Small difference between z-vertex from PV and from FT0"};
    Configurable<bool> cfgRequireIsVertexITSTPC{"cfgRequireIsVertexITSTPC", false, "At least one ITS-TPC track (reject vertices built from ITS-only tracks)"};
    Configurable<bool> cfgRequirekIsVertexTOFmatched{"cfgRequirekIsVertexTOFmatched", false, "Require kIsVertexTOFmatched: at least one of vertex contributors is matched to TOF"};
    Configurable<bool> useMultMCmidrap{"useMultMCmidrap", true, "use generated Nch in ∣eta∣ < 0.8"};
    Configurable<bool> cfgUseInelgt0wTVX{"cfgUseInelgt0wTVX", true, "Use INEL > 0 condition with TVX trigger, i.e. FT0A and FT0C acceptance"};
    Configurable<bool> cfgRemoveSplitVertex{"cfgRemoveSplitVertex", true, "Remove split vertices"};
  } evtSelOpt;

  struct : ConfigurableGroup {
    ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.4, 4.8, 5.2, 5.6, 6, 6.5, 7, 7.5, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30}, "pT binning"};
    ConfigurableAxis axisPtV0s{"axisPtV0s", {VARIABLE_WIDTH, 0, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.4, 4.8, 5.2, 5.6, 6, 6.4, 6.8, 7.2, 7.6, 8, 8.4, 8.8, 9.2, 9.6, 10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 14.5, 15, 20}, "pT V0s binning"};
    ConfigurableAxis axisFlatPerc{"axisFlatPerc", {102, -0.01, 1.01}, "Flattenicity percentiles binning"};
    ConfigurableAxis axisMultPerc{"axisMultPerc", {VARIABLE_WIDTH, 0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.}, "T0 percentiles binning"};
    ConfigurableAxis axisVertexZ{"axisVertexZ", {80, -20., 20.}, "Vertex z binning"};
    ConfigurableAxis axisMult{"axisMult", {301, -0.5, 300.5}, "Multiplicity binning"};
    ConfigurableAxis axisDCAxy{"axisDCAxy", {200, -5., 5.}, "DCAxy binning"};
    ConfigurableAxis axisDCAz{"axisDCAz", {200, -5., 5.}, "DCAz binning"};
    ConfigurableAxis axisPhi = {"axisPhi", {60, 0, constants::math::TwoPI}, "#varphi binning"};
    ConfigurableAxis axisPhiMod = {"axisPhiMod", {100, 0, constants::math::PI / 9}, "fmod(#varphi,#pi/9)"};
    ConfigurableAxis axisEta = {"axisEta", {50, -1.0, 1.0}, "#eta binning"};
    ConfigurableAxis axisRapidity = {"axisRapidity", {50, -1.0, 1.0}, "#it{y} binning"};
    ConfigurableAxis axisDedx{"axisDedx", {100, 0, 100}, "dE/dx binning"};
    ConfigurableAxis axisNsigmaTPC{"axisNsigmaTPC", {200, -10, 10}, "nsigmaTPC binning"};
    ConfigurableAxis axisNsigmaTOF{"axisNsigmaTOF", {200, -10, 10}, "nsigmaTOF binning"};
    ConfigurableAxis axisAmplFV0{"axisAmplFV0", {4096, 0, 4096}, "FV0 amplitude (ADC) binning"};
    ConfigurableAxis axisAmplFV0Sum{"axisAmplFV0Sum", {4096, 0, 4096 * 49}, "FV0 amplitude sum (ADC) binning"};
    ConfigurableAxis axisChannelFV0{"axisChannelFV0", {49, 0., 49.}, "FV0 channel ID binning"};
  } binOpt;

  struct : ConfigurableGroup {
    Configurable<float> cfgTrkEtaMax{"cfgTrkEtaMax", 0.8f, "Eta range for tracks"};
    Configurable<float> cfgRapMax{"cfgRapMax", 0.5f, "Maximum range of rapidity for tracks"};
    Configurable<float> cfgTrkPtMin{"cfgTrkPtMin", 0.1f, "Minimum pT of tracks"};
    Configurable<bool> cfgApplyNcl{"cfgApplyNcl", false, "Apply cut on TPC clusters"};
    Configurable<float> cfgNclTPCMin{"cfgNclTPCMin", 135.0f, "Minimum of number of TPC found clusters"};
    Configurable<bool> cfgApplyNclPID{"cfgApplyNclPID", true, "Apply cut on TPC PID clusters"};
    Configurable<float> cfgNclPidTPCMin{"cfgNclPidTPCMin", 135.0f, "Minimum of number of TPC PID clusters"};
    Configurable<float> cfgPhiCutPtMin{"cfgPhiCutPtMin", 2.0f, "Minimum pT for phi cut"};
    Configurable<float> cfgTOFBetaPion{"cfgTOFBetaPion", 1.0f, "Minimum beta for TOF pions"};
    Configurable<float> cfgTofBetaPiMax{"cfgTofBetaPiMax", 5E-5, "Maximum beta for TOF pion selection"};
    Configurable<bool> cfgRejectTrkAtTPCSector{"cfgRejectTrkAtTPCSector", true, "Reject tracks close to the TPC sector boundaries"};
    Configurable<std::string> cfgGeoTrkCutMin{"cfgGeoTrkCutMin", "0.06/x+pi/18.0-0.06", "ROOT TF1 formula for minimum phi cut in TPC"};
    Configurable<std::string> cfgGeoTrkCutMax{"cfgGeoTrkCutMax", "0.1/x+pi/18.0+0.06", "ROOT TF1 formula for maximum phi cut in TPC"};
    Configurable<float> cfgMomMIPMax{"cfgMomMIPMax", 0.6f, "Maximum momentum of MIP pions"};
    Configurable<float> cfgMomMIPMin{"cfgMomMIPMin", 0.4f, "Minimum momentum of MIP pions"};
    Configurable<float> cfgDeDxMIPMax{"cfgDeDxMIPMax", 60.0f, "Maximum range of MIP dedx"};
    Configurable<float> cfgDeDxMIPMin{"cfgDeDxMIPMin", 40.0f, "Maximum range of MIP dedx"};
    Configurable<float> cfgNsigmaMax{"cfgNsigmaMax", 100.0f, "Maximum range of nsgima for tracks"};
    Configurable<float> cfgDcaNsigmaCombinedMax{"cfgDcaNsigmaCombinedMax", 3.0f, "Maximum range of combined nsgima of tracks for DCA"};
    Configurable<float> cfgMomSelPiTOF{"cfgMomSelPiTOF", 0.4f, "Minimum momentum cut for TOF pions"};
    Configurable<float> cfgNsigSelKaTOF{"cfgNsigSelKaTOF", 3.0f, "Nsigma cut for TOF kaons"};
    Configurable<float> cfgBetaPlateuMax{"cfgBetaPlateuMax", 0.1f, "Beta max for Plateau electrons"};
  } trkSelOpt;

  struct : ConfigurableGroup {
    // common selection
    Configurable<int> cfgV0TypeSel{"cfgV0TypeSel", 1, "select on a certain V0 type (leave negative if no selection desired)"};
    Configurable<float> cfgV0Ymax{"cfgV0Ymax", 0.5f, "Maximum rapidity of V0s"};
    Configurable<bool> cfgRejectV0sAtTPCSector{"cfgRejectV0sAtTPCSector", true, "Reject V0s close to the TPC sector boundaries"};
    Configurable<bool> cfgRequireITS{"cfgRequireITS", true, "Additional cut on the ITS requirement"};
    Configurable<float> cfgNsigmaElTPC{"cfgNsigmaElTPC", 5.0, "max nsigma of TPC for electorn"};
    Configurable<float> cfgNsigmaPiTPC{"cfgNsigmaPiTPC", 5.0, "max nsigma of TPC for pion"};
    Configurable<float> cfgNsigmaPrTPC{"cfgNsigmaPrTPC", 5.0, "max nsigma of TPC for proton"};
    Configurable<float> cfgNsigmaElTOF{"cfgNsigmaElTOF", 3.0, "max nsigma of TOF for electorn"};
    Configurable<float> cfgNsigmaPiTOF{"cfgNsigmaPiTOF", 3.0, "max nsigma of TOF for pion"};
    Configurable<float> cfgNsigmaPrTOF{"cfgNsigmaPrTOF", 3.0, "max nsigma of TOF for proton"};
    ConfigurableAxis axisArmPodAlpha{"axisArmPodAlpha", {200, -1.0, 1.0}, "Armenteros-Podolanski alpha"};
    ConfigurableAxis axisArmPodqT{"axisArmPodqT", {600, 0.0f, 0.3f}, "Armenteros-Podolanski qT"};
    // standad parameters for V0 selection
    Configurable<float> cfgV0etamin{"cfgV0etamin", -0.8f, "min eta of V0s"};
    Configurable<float> cfgV0etamax{"cfgV0etamax", +0.8f, "max eta of V0s"};
    Configurable<float> cfgminNCrossedRowsTPC{"cfgminNCrossedRowsTPC", 70.f, "Additional cut on the minimum number of crossed rows in the TPC"};
    Configurable<bool> cfgApplyV0sNclFound{"cfgApplyV0sNclFound", false, "Apply cut on TPC Found clusters"};
    Configurable<float> cfgV0NclTPCMin{"cfgV0NclTPCMin", 135.0f, "Minimum of number of TPC found clusters"};
    Configurable<bool> cfgApplyV0sNclPID{"cfgApplyV0sNclPID", true, "Apply cut on TPC PID clusters"};
    Configurable<float> cfgV0NclPidTPCMin{"cfgV0NclPidTPCMin", 135.0f, "Minimum of number of TPC PID clusters"};
    Configurable<float> cfgmaxChi2PerClusterTPC{"cfgmaxChi2PerClusterTPC", 4.f, "Additional cut on the maximum value of the chi2 per cluster in the TPC"};
    Configurable<float> cfgmaxChi2PerClusterITS{"cfgmaxChi2PerClusterITS", 36.f, "Additional cut on the maximum value of the chi2 per cluster in the ITS"};
    Configurable<int> cfgminITSnClusters{"cfgminITSnClusters", 4, "minimum number of found ITS clusters"};
    Configurable<float> cfgminNCrossedRowsOverFindableClustersTPC{"cfgminNCrossedRowsOverFindableClustersTPC", 0.8f, "Additional cut on the minimum value of the ratio between crossed rows and findable clusters in the TPC"};
    Configurable<float> cfgDCAv0daughter{"cfgDCAv0daughter", 1.0, "max DCA of V0 daughter tracks (cm)"};
    Configurable<float> cfgv0cospa{"cfgv0cospa", 0.995, "min V0 CosPA"};
    Configurable<float> cfgDCAposToPV{"cfgDCAposToPV", 0.05f, "min DCA Pos To PV (cm)"};
    Configurable<float> cfgDCAnegToPV{"cfgDCAnegToPV", 0.05f, "min DCA Neg To PV (cm)"};
    Configurable<float> cfgv0Rmin{"cfgv0Rmin", 1.2, "min V0 radius (cm)"};
    Configurable<float> cfgv0Rmax{"cfgv0Rmax", 1E5, "max V0 radius (cm)"};
    // parameters for selection KOs
    Configurable<float> cfgcTauK0s{"cfgcTauK0s", 20, "v0ctau for K0s"};
    Configurable<float> cfgCosPAK0s{"cfgCosPAK0s", 0.995, "V0 CosPA for K0s"};
    Configurable<float> cfgV0radiusK0s{"cfgV0radiusK0s", 0.5, "v0radius for K0s"};
    Configurable<float> cfgdmassK{"cfgdmassK", 0.01f, "Competing Mass Rejection cut for K0s"};
    Configurable<float> cfgArmPodK0s{"cfgArmPodK0s", 5.0f, "pT * (cut) > |alpha|, Armenteros-Podolanski cut for K0s"};
    ConfigurableAxis axisK0sMass{"axisK0sMass", {200, 0.4f, 0.6f}, "K0Short mass binning"};
    // parameters for selection Lambda / antiLambda
    Configurable<float> cfgcTauLambda{"cfgcTauLambda", 30, "v0ctau for Lambda"};
    Configurable<float> cfgCosPALambda{"cfgCosPALambda", 0.995, "V0 CosPA for Lambda"};
    Configurable<float> cfgV0radiusLambda{"cfgV0radiusLambda", 0.5, "v0radius for Lambda"};
    Configurable<float> cfgdmassL{"cfgdmassL", 0.01f, "Competing Mass Rejection cut for Lambda"};
    ConfigurableAxis axisLambdaMass{"axisLambdaMass", {200, 1.101f, 1.131f}, "Lambda mass binning"};
    // parameters for selection Gamma
    Configurable<float> cfgdmassG{"cfgdmassG", 0.1f, "max mass for Gammas"};
    Configurable<float> cfgArmPodGammasalpha{"cfgArmPodGammasalpha", 0.45f, "Armenteros-Podolanski alpha cut for Gammas"};
    Configurable<float> cfgArmPodGammasqT{"cfgArmPodGammasqT", 0.01f, "Armenteros-Podolanski qT cut for Gammas"};
    ConfigurableAxis axisGammaMass{"axisGammaMass", {200, 0.0f, 0.5f}, "Gamma mass binning"};
    Configurable<float> cfgdEdxPlateauSel{"cfgdEdxPlateauSel", 50, "dEdx selection sensitivity for electrons"};
  } v0SelOpt;

  Service<ccdb::BasicCCDBManager> ccdb;
  struct : ConfigurableGroup {
    Configurable<float> cfgMagField{"cfgMagField", 99999, "Configurable magnetic field;default CCDB will be queried"};
    Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
    Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  } ccdbConf;

  struct DeDxCalib {
    TList* lCalibObjects = nullptr;
    TH1F* hMIPcalibPos = nullptr;
    TH1F* hMIPcalibNeg = nullptr;
    TH1F* hMIPcalibAll = nullptr;
    TH1F* hPlateauCalibPos = nullptr;
    TH1F* hPlateauCalibNeg = nullptr;
    TH1F* hPlateauCalibAll = nullptr;
    bool lCalibLoaded = false;
  } dedxcalib;

  struct : ConfigurableGroup {
    Configurable<bool> requireRCTFlagChecker{"requireRCTFlagChecker", false, "Check event quality in run condition table"};
    Configurable<std::string> cfgEvtRCTFlagCheckerLabel{"cfgEvtRCTFlagCheckerLabel", "CBT_hadronPID", "Evt sel: RCT flag checker label"};
    Configurable<bool> cfgEvtRCTFlagCheckerZDCCheck{"cfgEvtRCTFlagCheckerZDCCheck", false, "Evt sel: RCT flag checker ZDC check"};
    Configurable<bool> cfgEvtRCTFlagCheckerLimitAcceptAsBad{"cfgEvtRCTFlagCheckerLimitAcceptAsBad", true, "Evt sel: RCT flag checker treat Limited Acceptance As Bad"};
  } rctCuts;

  RCTFlagsChecker rctChecker;

  TrackSelection selTrkGlobal;
  Configurable<bool> isCustomTracks{"isCustomTracks", true, "Use custom track cuts"};
  Configurable<float> minPt{"minPt", 0.1f, "Set minimum pT of tracks"};
  Configurable<float> maxPt{"maxPt", 1e10f, "Set maximum pT of tracks"};
  Configurable<float> requireEta{"requireEta", 0.8f, "Set eta range of tracks"};
  Configurable<int> setITSreq{"setITSreq", 0, "0 = Run3ITSibAny, 1 = Run3ITSallAny, 2 = Run3ITSall7Layers, 3 = Run3ITSibTwo"};
  Configurable<bool> requireITSminCl{"requireITSminCl", false, "Require additional cut on ITS clusters"};
  Configurable<int> setITSminCl{"setITSminCl", 7, "Additional cut on ITS clusters"};
  Configurable<bool> requireITS{"requireITS", true, "Additional cut on the ITS requirement"};
  Configurable<bool> requireTPC{"requireTPC", true, "Additional cut on the TPC requirement"};
  Configurable<bool> requireGoldenChi2{"requireGoldenChi2", true, "Additional cut on the GoldenChi2"};
  Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 70.f, "Additional cut on the minimum number of crossed rows in the TPC"};
  Configurable<float> minNCrossedRowsOverFindableClustersTPC{"minNCrossedRowsOverFindableClustersTPC", 0.8f, "Additional cut on the minimum value of the ratio between crossed rows and findable clusters in the TPC"};
  Configurable<float> maxChi2PerClusterTPC{"maxChi2PerClusterTPC", 4.f, "Additional cut on the maximum value of the chi2 per cluster in the TPC"};
  Configurable<float> maxChi2PerClusterITS{"maxChi2PerClusterITS", 36.f, "Additional cut on the maximum value of the chi2 per cluster in the ITS"};
  Configurable<int> minITSnClusters{"minITSnClusters", 5, "minimum number of found ITS clusters"};
  Configurable<float> maxDcaXYFactor{"maxDcaXYFactor", 1.f, "Multiplicative factor on the maximum value of the DCA xy"};
  Configurable<float> maxDcaZ{"maxDcaZ", 2.f, "Additional cut on the maximum value of the DCA z"};

  TF1* fPhiCutLow = nullptr;
  TF1* fPhiCutHigh = nullptr;

  std::vector<std::unique_ptr<TF1>> fDeDxVsEta;
  std::vector<std::vector<float>> vecParamsMIP;
  std::vector<std::unique_ptr<TF1>> fEDeDxVsEta;
  std::vector<std::vector<float>> vecParamsPLA;

  using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::FT0sCorrected, aod::CentFT0As, aod::CentFT0Cs>;
  using Colls = soa::Join<aod::Collisions, aod::EvSels, aod::TPCMults, aod::PVMults, aod::MultZeqs, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>;
  using CollsGen = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::TPCMults, aod::PVMults, aod::MultZeqs, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>;
  using MCColls = soa::Join<aod::McCollisions, aod::McCentFT0Ms, aod::MultsExtraMC>;
  using CollsMCExtraMult = soa::Join<aod::McCollisions, aod::McCentFT0Ms, aod::MultMCExtras, aod::McCollsExtra>;
  using CollsGenSgn = soa::SmallGroups<soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::TPCMults, aod::PVMults, aod::MultZeqs, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>>;
  using MyPIDTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::TrackSelectionExtension, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFbeta, aod::TOFSignal, aod::pidTOFFlags>;
  using MyLabeledTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksDCA, aod::McTrackLabels>;
  using MyFiltLabeledTracks = soa::Filtered<MyLabeledTracks>;
  using MyLabeledPIDTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::TrackSelectionExtension, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFbeta, aod::McTrackLabels>;

  void init(InitContext&)
  {
    auto vecParamsMIPposEtaP = (std::vector<float>)paramsFuncMIPposEtaP;
    auto vecParamsMIPposEtaN = (std::vector<float>)paramsFuncMIPposEtaN;
    auto vecParamsMIPnegEtaP = (std::vector<float>)paramsFuncMIPnegEtaP;
    auto vecParamsMIPnegEtaN = (std::vector<float>)paramsFuncMIPnegEtaN;
    auto vecParamsMIPallEtaP = (std::vector<float>)paramsFuncMIPallEtaP;
    auto vecParamsMIPallEtaN = (std::vector<float>)paramsFuncMIPallEtaN;

    auto vecParamsPLAposEtaP = (std::vector<float>)paramsFuncPlateaUposEtaP;
    auto vecParamsPLAposEtaN = (std::vector<float>)paramsFuncPlateaUposEtaN;
    auto vecParamsPLAnegEtaP = (std::vector<float>)paramsFuncPlateaUnegEtaP;
    auto vecParamsPLAnegEtaN = (std::vector<float>)paramsFuncPlateaUnegEtaN;
    auto vecParamsPLAallEtaP = (std::vector<float>)paramsFuncPlateaUallEtaP;
    auto vecParamsPLAallEtaN = (std::vector<float>)paramsFuncPlateaUallEtaN;

    auto addVec = [&](std::vector<std::vector<float>>& targetVec, const std::string& name, bool isMIP) {
      if (isMIP) {
        targetVec.emplace_back(vecParamsMIPposEtaP);
        targetVec.emplace_back(vecParamsMIPnegEtaP);
        targetVec.emplace_back(vecParamsMIPallEtaP);
        targetVec.emplace_back(vecParamsMIPposEtaN);
        targetVec.emplace_back(vecParamsMIPnegEtaN);
        targetVec.emplace_back(vecParamsMIPallEtaN);
        if (!vecParamsMIP.size()) {
          LOG(info) << "size of " << name << "is zero.";
        }
      } else {
        targetVec.emplace_back(vecParamsPLAposEtaP);
        targetVec.emplace_back(vecParamsPLAnegEtaP);
        targetVec.emplace_back(vecParamsPLAallEtaP);
        targetVec.emplace_back(vecParamsPLAposEtaN);
        targetVec.emplace_back(vecParamsPLAnegEtaN);
        targetVec.emplace_back(vecParamsPLAallEtaN);
        if (!vecParamsPLA.size()) {
          LOG(info) << "size of " << name << "is zero.";
        }
      }
    };
    addVec(vecParamsMIP, "vecParamsMIP", true);
    std::transform(std::begin(vecParamsMIP), std::end(vecParamsMIP), std::back_inserter(fDeDxVsEta), [&](auto const& params) {
      return setFuncPars(params);
    });
    addVec(vecParamsPLA, "vecParamsPLA", false);
    std::transform(std::begin(vecParamsPLA), std::end(vecParamsPLA), std::back_inserter(fEDeDxVsEta), [&](auto const& params) {
      return setFuncPars(params);
    });

    ccdb->setURL(ccdbConf.ccdbUrl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    ccdb->setFatalWhenNull(false);

    rctChecker.init(rctCuts.cfgEvtRCTFlagCheckerLabel, rctCuts.cfgEvtRCTFlagCheckerZDCCheck, rctCuts.cfgEvtRCTFlagCheckerLimitAcceptAsBad);

    if (isCustomTracks.value) {
      selTrkGlobal = getGlobalTrackSelectionRun3ITSMatch(setITSreq.value);
      selTrkGlobal.SetPtRange(minPt.value, maxPt.value);
      selTrkGlobal.SetEtaRange(-requireEta.value, requireEta.value);
      selTrkGlobal.SetRequireITSRefit(requireITS.value);
      selTrkGlobal.SetRequireTPCRefit(requireTPC.value);
      selTrkGlobal.SetRequireGoldenChi2(requireGoldenChi2.value);
      selTrkGlobal.SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC.value);
      selTrkGlobal.SetMaxChi2PerClusterITS(maxChi2PerClusterITS.value);
      selTrkGlobal.SetMinNClustersITS(minITSnClusters.value);
      if (requireITSminCl.value) {
        selTrkGlobal.SetRequireHitsInITSLayers(setITSminCl.value, {0, 1, 2, 3, 4, 5, 6});
      }
      selTrkGlobal.SetMinNCrossedRowsTPC(minNCrossedRowsTPC.value);
      selTrkGlobal.SetMinNCrossedRowsOverFindableClustersTPC(minNCrossedRowsOverFindableClustersTPC.value);
      // //     selTrkGlobal.SetMaxDcaXYPtDep([](float pt) { return 0.0105f + 0.0350f / pow(pt, 1.1f); });
      selTrkGlobal.SetMaxDcaXYPtDep([](float /*pt*/) { return 10000.f; });
      selTrkGlobal.SetMaxDcaZ(maxDcaZ.value);
      selTrkGlobal.print();
    }

    const AxisSpec chargeAxis{2, -2.f, 2.f, "Charge"};
    const AxisSpec dEdxAxis{binOpt.axisDedx, "TPC dEdx (a.u.)"};
    const AxisSpec vtxzAxis{binOpt.axisVertexZ, "Z_{vtx} (cm)"};
    const AxisSpec flatAxis{binOpt.axisFlatPerc, "Flat FV0"};
    const AxisSpec nChAxis{binOpt.axisMult, "Nch, |eta|<0.8"};
    const AxisSpec etaAxis{binOpt.axisEta, "#eta"};
    const AxisSpec rapidityAxis{binOpt.axisRapidity, "#it{y}"};
    const AxisSpec phiAxis{binOpt.axisPhi, "#varphi"};
    const AxisSpec phiAxisMod{binOpt.axisPhiMod, "fmod(#varphi,#pi/9)"};
    const AxisSpec pAxis{binOpt.axisPt, "#it{p} (GeV/#it{c})"};
    const AxisSpec ptAxis{binOpt.axisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec ptAxisV0s{binOpt.axisPtV0s, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec dcaXYAxis{binOpt.axisDCAxy, "DCA_{xy} (cm)"};
    const AxisSpec dcaZAxis{binOpt.axisDCAz, "DCA_{z} (cm)"};
    const AxisSpec shCluserAxis{200, 0, 1, "Fraction of shared TPC clusters"};
    const AxisSpec clTpcAxis{160, 0, 160, "Number of clusters in TPC"};
    const AxisSpec nSigmaTPCAxis{binOpt.axisNsigmaTPC, "n#sigma_{TPC}"};
    const AxisSpec nSigmaTOFAxis{binOpt.axisNsigmaTOF, "n#sigma_{TOF}"};
    const AxisSpec amplitudeFV0{binOpt.axisAmplFV0, "FV0 amplitude (ADC)"};
    const AxisSpec amplitudeFV0Sum{binOpt.axisAmplFV0Sum, "FV0 amplitude sm (ADC)"};
    const AxisSpec channelFV0Axis{binOpt.axisChannelFV0, "FV0 channel ID"};

    AxisSpec multAxis{binOpt.axisMultPerc, "multiplicity estimator"};

    switch (multEst) {
      case MultE::CnoMult:
        break;
      case MultE::CmultFT0M:
        multAxis.name = "multFT0M";
        break;
      case MultE::CmultTPC:
        multAxis.name = "multTPC";
        break;
      default:
        LOG(fatal) << "No valid option for mult estimator " << multEst;
    }

    if (trkSelOpt.cfgRejectTrkAtTPCSector || v0SelOpt.cfgRejectV0sAtTPCSector) {
      fPhiCutLow = new TF1("fPhiCutLow", trkSelOpt.cfgGeoTrkCutMin.value.c_str(), 0, 100);
      fPhiCutHigh = new TF1("fPhiCutHigh", trkSelOpt.cfgGeoTrkCutMax.value.c_str(), 0, 100);
    }

    registryQC.add("Events/hVtxZ", "Measured vertex z position", kTH1F, {vtxzAxis});

    // Event counter
    registryQC.add("Events/hEvtSel", "Number of events; Cut; #Events Passed Cut", {kTH1F, {{nEvtSel, 0, nEvtSel}}});
    registryQC.get<TH1>(HIST("Events/hEvtSel"))->GetXaxis()->SetBinLabel(evtSelAll + 1, "Events read");
    registryQC.get<TH1>(HIST("Events/hEvtSel"))->GetXaxis()->SetBinLabel(evtSelSel8 + 1, "Evt. sel8");
    registryQC.get<TH1>(HIST("Events/hEvtSel"))->GetXaxis()->SetBinLabel(evtSelNoITSROFrameBorder + 1, "NoITSROFrameBorder");
    registryQC.get<TH1>(HIST("Events/hEvtSel"))->GetXaxis()->SetBinLabel(evtSelkNoTimeFrameBorder + 1, "NoTimeFrameBorder");
    registryQC.get<TH1>(HIST("Events/hEvtSel"))->GetXaxis()->SetBinLabel(evtSelkNoSameBunchPileup + 1, "NoSameBunchPileup");
    registryQC.get<TH1>(HIST("Events/hEvtSel"))->GetXaxis()->SetBinLabel(evtSelkIsGoodZvtxFT0vsPV + 1, "IsGoodZvtxFT0vsPV");
    registryQC.get<TH1>(HIST("Events/hEvtSel"))->GetXaxis()->SetBinLabel(evtSelkIsVertexITSTPC + 1, "IsVertexITSTPC");
    registryQC.get<TH1>(HIST("Events/hEvtSel"))->GetXaxis()->SetBinLabel(evtSelkIsVertexTOFmatched + 1, "IsVertexTOFmatched");
    registryQC.get<TH1>(HIST("Events/hEvtSel"))->GetXaxis()->SetBinLabel(evtSelVtxZ + 1, "Vtx-z pos");
    registryQC.get<TH1>(HIST("Events/hEvtSel"))->GetXaxis()->SetBinLabel(evtSelINELgt0 + 1, "INEL>0");
    registryQC.get<TH1>(HIST("Events/hEvtSel"))->GetXaxis()->SetBinLabel(evtSelRCTFlagChecker + 1, "RCT Flag Checker");
    // Number of tracks vs centrality
    registryQC.add("Events/hNchVsCent", "Measured Nch vs Cent; centrality; Nch (|#eta|<0.8)", {kTH2F, {nChAxis, multAxis}});
    // FV0 QA
    registryQC.add("FV0/hFV0AmplWCalib", "", {kTH2F, {channelFV0Axis, amplitudeFV0}});
    registryQC.add("FV0/hFV0AmplvsVtxzWoCalib", "", {kTH2F, {vtxzAxis, amplitudeFV0Sum}});
    registryQC.add("FV0/hFV0AmplvsVtxzCalib", "", {kTH2F, {vtxzAxis, amplitudeFV0Sum}});
    registryQC.add("FV0/hFV0amp", "", {kTH2F, {channelFV0Axis, amplitudeFV0}});
    registryQC.add("FV0/pFV0amp", "", kTProfile, {channelFV0Axis});
    registryQC.add("FV0/hFV0ampCorr", "", {kTH2F, {channelFV0Axis, amplitudeFV0}});

    LOG(info) << "Size of the QC histograms:";
    registryQC.print();

    if (doprocessFlat) {
      // Track counter
      registryData.add("Tracks/hTrkSel", "Number of tracks; Cut; #Tracks Passed Cut", {kTH1F, {{nTrkSel, 0, nTrkSel}}});
      registryData.get<TH1>(HIST("Tracks/hTrkSel"))->GetXaxis()->SetBinLabel(trkSelAll + 1, "All");
      registryData.get<TH1>(HIST("Tracks/hTrkSel"))->GetXaxis()->SetBinLabel(trkSelEta + 1, "Eta");
      registryData.get<TH1>(HIST("Tracks/hTrkSel"))->GetXaxis()->SetBinLabel(trkSelPt + 1, "Pt");
      registryData.get<TH1>(HIST("Tracks/hTrkSel"))->GetXaxis()->SetBinLabel(trkSelDCA + 1, "DCA");
      registryData.get<TH1>(HIST("Tracks/hTrkSel"))->GetXaxis()->SetBinLabel(trkNRowsTPC + 1, "trkNRowsTPC");
      registryData.get<TH1>(HIST("Tracks/hTrkSel"))->GetXaxis()->SetBinLabel(trkSelNClsFound + 1, "NClsTPCFound");
      registryData.get<TH1>(HIST("Tracks/hTrkSel"))->GetXaxis()->SetBinLabel(trkSelNClsPID + 1, "NClsTPCPid");
      registryData.get<TH1>(HIST("Tracks/hTrkSel"))->GetXaxis()->SetBinLabel(trkSelTPCBndr + 1, "TPC Boundary");
      // V0 counter
      registryData.add("Tracks/V0qa/hV0Sel", "Number of V0s; Cut; #Tracks Passed Cut", {kTH1F, {{nV0Sel, 0, nV0Sel}}});
      registryData.get<TH1>(HIST("Tracks/V0qa/hV0Sel"))->GetXaxis()->SetBinLabel(v0SelAll + 1, "All");
      registryData.get<TH1>(HIST("Tracks/V0qa/hV0Sel"))->GetXaxis()->SetBinLabel(v0SelRejectSameSign + 1, "Reject same sign");
      registryData.get<TH1>(HIST("Tracks/V0qa/hV0Sel"))->GetXaxis()->SetBinLabel(v0SelRejectV0sAtTPCSector + 1, "Reject V0s at TPC sector");
      registryData.get<TH1>(HIST("Tracks/V0qa/hV0Sel"))->GetXaxis()->SetBinLabel(v0SelCosPA + 1, "Cos PA");
      registryData.get<TH1>(HIST("Tracks/V0qa/hV0Sel"))->GetXaxis()->SetBinLabel(v0SelV0radius + 1, "V0 radius");
      registryData.get<TH1>(HIST("Tracks/V0qa/hV0Sel"))->GetXaxis()->SetBinLabel(v0SelDCAposToPV + 1, "DCA pos to PV");
      registryData.get<TH1>(HIST("Tracks/V0qa/hV0Sel"))->GetXaxis()->SetBinLabel(v0SelDaughters + 1, "V0 daughters' sel.");
      registryData.get<TH1>(HIST("Tracks/V0qa/hV0Sel"))->GetXaxis()->SetBinLabel(v0SelDCAv0daughter + 1, "DCA v0 daughter");

      registryData.add("Events/hFlatVsMultEst", "hFlatVsMultEst", kTH2F, {flatAxis, multAxis});
      registryData.add("Tracks/postSel/hPVsPtEta", "; #it{p} (GeV/#it{c}); #it{p}_{T} (GeV/#it{c}); #eta;", {kTH3F, {pAxis, ptAxis, etaAxis}});
      if (cfgFillNclVsPhiCutQaHist || cfgFillTrackQaHist || cfgFilldEdxQaHist || cfgFillDCAxyHist) {
        if (cfgFillNclVsPhiCutQaHist) {
          registryData.add("Tracks/postSel/hPtPhi", "; #it{p}_{T} (GeV/#it{c}); fmod(#varphi,#pi/9)", {kTH2F, {ptAxis, phiAxisMod}});
          registryData.add("Tracks/postSel/hPtPhiNclTPC", "; #{eta}; #it{p}_{T} (GeV/#it{c}); fmod(#varphi,#pi/9); N_{cluster}", {kTHnSparseF, {etaAxis, ptAxis, phiAxisMod, clTpcAxis}});
          registryData.add("Tracks/postSel/hPtPhiNclPIDTPC", "; #{eta}; #it{p}_{T} (GeV/#it{c}); fmod(#varphi,#pi/9); N_{PID cluster}", {kTHnSparseF, {etaAxis, ptAxis, phiAxisMod, clTpcAxis}});
          registryData.add("Tracks/postSel/hPtNclTPC", "; #it{p}_{T} (GeV/#it{c}); N_{cluster}", {kTH2F, {ptAxis, clTpcAxis}});
          registryData.add("Tracks/postSel/pPtNclTPC", "; #it{p}_{T} (GeV/#it{c}); N_{cluster}", {kTProfile, {ptAxis}});
          registryData.add("Tracks/postSel/hPtNclPIDTPC", "; #it{p}_{T} (GeV/#it{c}); N_{PID cluster}", {kTH2F, {ptAxis, clTpcAxis}});
          registryData.add("Tracks/postSel/pPtNclPIDTPC", "; #it{p}_{T} (GeV/#it{c}); N_{PID cluster}", {kTProfile, {ptAxis}});
        }
        if (cfgFillTrackQaHist) {
          registryData.add("Tracks/postSel/hShTpcClvsPt", "", {kTH2F, {ptAxis, shCluserAxis}});
          registryData.add("Tracks/postSel/hNclTPCFoundvsPt", "", {kTH2F, {ptAxis, clTpcAxis}});
          registryData.add("Tracks/postSel/hNClTPCPidvsPt", "", {kTH2F, {ptAxis, clTpcAxis}});
          registryData.add("Tracks/postSel/hNclTPCFoundvsEta", "", {kTH2F, {etaAxis, clTpcAxis}});
          registryData.add("Tracks/postSel/hNClTPCPidvsEta", "", {kTH2F, {etaAxis, clTpcAxis}});
          registryData.add("Tracks/postSel/hTPCCluster", "N_{cluster}", kTH1F, {clTpcAxis});
          registryData.add("Tracks/postSel/hPtVsWOcutDCA", "hPtVsWOcutDCA", kTH2F, {ptAxis, dcaXYAxis});
          registryData.add("Tracks/postSel/hPt", "", kTH1F, {ptAxis});
          registryData.add("Tracks/postSel/hPhi", "", kTH1F, {phiAxis});
          registryData.add("Tracks/postSel/hEta", "", kTH1F, {etaAxis});
          registryData.add("Tracks/postSel/hDCAXYvsPt", "", kTH2F, {ptAxis, dcaXYAxis});
          registryData.add("Tracks/postSel/hDCAZvsPt", "", kTH2F, {ptAxis, dcaZAxis});
          // tpc
          registryData.add("Tracks/postSel/hTPCnClsShared", " ; # shared TPC clusters TPC", kTH1F, {{165, -0.5, 164.5}});
          registryData.add("Tracks/postSel/hTPCcrossedRows", " ; # crossed TPC rows", kTH1F, {{165, -0.5, 164.5}});
          registryData.add("Tracks/postSel/hTPCcrossedRowsOverFindableCls", " ; crossed rows / findable TPC clusters", kTH1F, {{60, 0.7, 1.3}});
          // its
          registryData.add("Tracks/postSel/hITSnCls", " ; # ITS clusters", kTH1F, {{8, -0.5, 7.5}});
          registryData.add("Tracks/postSel/hChi2ITSTrkSegment", "chi2ITS", kTH1F, {{100, -0.5, 99.5}});
          // tof
          registryData.add("Tracks/postSel/hTOFPvsBeta", "Beta from TOF; #it{p} (GeV/#it{c}); #beta", {kTH2F, {pAxis, {120, 0.0, 1.2}}});
          registryData.add("Tracks/postSel/hTOFpi", "Primary Pions from TOF; #eta; #it{p} (GeV/#it{c}); dEdx", {kTHnSparseF, {etaAxis, pAxis, dEdxAxis}});
        }
        if (cfgFilldEdxQaHist) {
          if (cfgStoreThnSparse) {
            registryData.add("Tracks/postCalib/all/hMIP", "; mult; flat; #eta; #LT dE/dx #GT_{MIP, primary tracks};", {kTHnSparseF, {multAxis, flatAxis, etaAxis, dEdxAxis}});
            registryData.add("Tracks/postCalib/all/hPlateau", "; mult; flat; #eta; #LT dE/dx #GT_{Plateau, primary tracks};", {kTHnSparseF, {multAxis, flatAxis, etaAxis, dEdxAxis}});
          } else {
            registryData.add("Tracks/postCalib/all/hMIP", "; #eta; #LT dE/dx #GT_{MIP, primary tracks};", {kTH2F, {etaAxis, dEdxAxis}});
            registryData.add("Tracks/postCalib/all/hPlateau", "; #eta; #LT dE/dx #GT_{Plateau, primary tracks};", {kTH2F, {etaAxis, dEdxAxis}});
          }
          registryData.add("Tracks/postCalib/all/hMIPVsPhi", "; #varphi; #LT dE/dx #GT_{MIP, primary tracks};", {kTH2F, {phiAxis, dEdxAxis}});
          registryData.add("Tracks/postCalib/all/pMIPVsPhi", "; #varphi; #LT dE/dx #GT_{MIP, primary tracks};", {kTProfile, {phiAxis}});
          registryData.add("Tracks/postCalib/all/hMIPVsPhiVsEta", "; #varphi; #LT dE/dx #GT_{MIP, primary tracks}; #eta;", {kTH3F, {phiAxis, dEdxAxis, etaAxis}});
          registryData.add("Tracks/postCalib/all/hPlateauVsPhi", "; #varphi; #LT dE/dx #GT_{Plateau, primary tracks};", {kTH2F, {phiAxis, dEdxAxis}});
          registryData.add("Tracks/postCalib/all/pPlateauVsPhi", "; #varphi; #LT dE/dx #GT_{Plateau, primary tracks};", {kTProfile, {phiAxis}});
          registryData.add("Tracks/postCalib/all/hPlateauVsPhiVsEta", "; #varphi; #LT dE/dx #GT_{Plateau, primary tracks}; #eta;", {kTH3F, {phiAxis, dEdxAxis, etaAxis}});
          registryData.addClone("Tracks/postCalib/all/", "Tracks/preCalib/all/");
          if (cfgFillChrgType) {
            registryData.addClone("Tracks/postCalib/all/", "Tracks/postCalib/pos/");
            registryData.addClone("Tracks/postCalib/all/", "Tracks/postCalib/neg/");
            registryData.addClone("Tracks/preCalib/all/", "Tracks/preCalib/pos/");
            registryData.addClone("Tracks/preCalib/all/", "Tracks/preCalib/neg/");
          }
        }
        if (cfgFillDCAxyHist) {
          for (int i = 0; i < Npart; i++) {
            registryData.add({fmt::format(CpTvsDCAxyF.data(), CspeciesAll[i]).c_str(), "; mult; flat; #it{p} (GeV/#it{c}); DCA_{xy} (cm)", {kTHnSparseF, {multAxis, flatAxis, ptAxis, dcaXYAxis}}});
          }
        }
      }
      registryData.addClone("Tracks/postSel/", "Tracks/preSel/");
      // V0's QA
      registryData.add("Tracks/V0qa/hV0Pt", "pT", kTH1F, {ptAxisV0s});
      registryData.add("Tracks/V0qa/hV0ArmPod", ";#alpha; #it{q}_T (GeV/c)", kTH2F, {v0SelOpt.axisArmPodAlpha, v0SelOpt.axisArmPodqT});
      // daughters' QA
      registryData.add("Tracks/V0qa/el/Ga/hArmPod", ";#alpha; #it{q}_T (GeV/c)", kTH2F, {v0SelOpt.axisArmPodAlpha, v0SelOpt.axisArmPodqT});
      registryData.add("Tracks/V0qa/pi/K0s/hArmPod", ";#alpha; #it{q}_T (GeV/c)", kTH2F, {v0SelOpt.axisArmPodAlpha, v0SelOpt.axisArmPodqT});
      registryData.add("Tracks/V0qa/el/Ga/hNclVsEta", ";#eta; #it{N}^{TPC}_cl", kTH2F, {etaAxis, clTpcAxis});
      registryData.add("Tracks/V0qa/pi/K0s/hNclVsEta", ";#eta; #it{N}^{TPC}_cl", kTH2F, {etaAxis, clTpcAxis});
      registryData.add("Tracks/V0qa/el/Ga/hdEdxMIPVsEta", ";#eta; dE/dx", kTH2F, {etaAxis, dEdxAxis});
      registryData.add("Tracks/V0qa/pi/K0s/hdEdxMIPVsEta", ";#eta; dE/dx", kTH2F, {etaAxis, dEdxAxis});
      registryData.addClone("Tracks/V0qa/pi/K0s/", "Tracks/V0qa/pi/La/");
      registryData.addClone("Tracks/V0qa/pi/K0s/", "Tracks/V0qa/pi/ALa/");
      registryData.addClone("Tracks/V0qa/pi/La/", "Tracks/V0qa/pr/La/");
      registryData.addClone("Tracks/V0qa/pi/ALa/", "Tracks/V0qa/pr/ALa/");

      // dEdx PID
      registryData.add({"Tracks/all/hdEdx", "; #eta; mult; flat; #it{p} (GeV/#it{c}); dEdx", {kTHnSparseF, {etaAxis, multAxis, flatAxis, pAxis, dEdxAxis}}});
      // Clean samples
      if (cfgFillV0Hist) {
        if (cfgStoreThnSparse) {
          registryData.add({"Tracks/CleanTof/all/hPiTof", "; #eta; mult; flat; #it{p} (GeV/#it{c}); dEdx", {kTHnSparseF, {etaAxis, multAxis, flatAxis, pAxis, dEdxAxis}}});
          registryData.add({"Tracks/CleanV0/all/hEV0", "; #eta; mult; flat; #it{p} (GeV/#it{c}); dEdx", {kTHnSparseF, {etaAxis, multAxis, flatAxis, pAxis, dEdxAxis}}});
          registryData.add({"Tracks/CleanV0/all/hPiV0", "; #eta; mult; flat; #it{p} (GeV/#it{c}); dEdx", {kTHnSparseF, {etaAxis, multAxis, flatAxis, pAxis, dEdxAxis}}});
          registryData.add({"Tracks/CleanV0/all/hPV0", "; #eta; mult; flat; #it{p} (GeV/#it{c}); dEdx", {kTHnSparseF, {etaAxis, multAxis, flatAxis, pAxis, dEdxAxis}}});
        } else {
          registryData.add({"Tracks/CleanTof/all/hPiTof", "; #eta; mult; flat; #it{p} (GeV/#it{c}); dEdx", {kTH3F, {etaAxis, pAxis, dEdxAxis}}});
          registryData.add({"Tracks/CleanV0/all/hEV0", "; #eta; mult; flat; #it{p} (GeV/#it{c}); dEdx", {kTH3F, {etaAxis, pAxis, dEdxAxis}}});
          registryData.add({"Tracks/CleanV0/all/hPiV0", "; #eta; mult; flat; #it{p} (GeV/#it{c}); dEdx", {kTH3F, {etaAxis, pAxis, dEdxAxis}}});
          registryData.add({"Tracks/CleanV0/all/hPV0", "; #eta; mult; flat; #it{p} (GeV/#it{c}); dEdx", {kTH3F, {etaAxis, pAxis, dEdxAxis}}});
        }
        registryData.add("Tracks/CleanTof/all/hBetaVsP", ";Momentum (GeV/#it{c}); #beta", kTH2F, {{{ptAxisV0s}, {120, 0., 1.2}}});
        registryData.add("Tracks/CleanTof/all/hTofExpPi", ";Momentum (GeV/#it{c});#it{t}^{#pi}_{Exp}/#it{t}_{TOF}", kTH2F, {{{ptAxisV0s}, {100, 0.2, 1.2}}});
        if (cfgFillChrgType) {
          registryData.addClone("Tracks/CleanTof/all/", "Tracks/CleanTof/pos/");
          registryData.addClone("Tracks/CleanTof/all/", "Tracks/CleanTof/neg/");
          registryData.addClone("Tracks/CleanV0/all/", "Tracks/CleanV0/pos/");
          registryData.addClone("Tracks/CleanV0/all/", "Tracks/CleanV0/neg/");
        }
      }
      if (cfgFillChrgType) {
        registryData.addClone("Tracks/all/", "Tracks/pos/");
        registryData.addClone("Tracks/all/", "Tracks/neg/");
      }
      LOG(info) << "Size of the Data histograms:";
      registryData.print();
    }

    if (doprocessMC) {
      registryMC.add({"Events/ResponseGen", ";N_{ch,FV0};1-#rho_{FV0};", {kTHnSparseF, {multAxis, flatAxis}}});
      registryMC.add("Events/h1flatencityFV0MCGen", "", {kTH1F, {flatAxis}});
      registryMC.add("Events/hFlatMCGen", "Events/hFlatMCGen", {kTH1F, {flatAxis}});
      registryMC.add("Events/hFlatMCRec", "Events/hFlatMCRec", {kTH1F, {flatAxis}});
      // Event counter
      auto h = registryMC.add<TH1>("Events/hEvtGenRec", "Generated and Reconstructed MC Collisions", kTH1F, {{3, 0.5, 3.5}});
      h->GetXaxis()->SetBinLabel(1, "Gen coll");
      h->GetXaxis()->SetBinLabel(2, "Rec coll");
      h->GetXaxis()->SetBinLabel(3, "INEL>0");
      registryMC.add("Events/hEvtMcGen", "Events/hEvtMcGen", {kTH1F, {{4, 0.f, 4.f}}});
      registryMC.get<TH1>(HIST("Events/hEvtMcGen"))->GetXaxis()->SetBinLabel(1, "all");
      registryMC.get<TH1>(HIST("Events/hEvtMcGen"))->GetXaxis()->SetBinLabel(2, "z-vtx");
      registryMC.get<TH1>(HIST("Events/hEvtMcGen"))->GetXaxis()->SetBinLabel(3, "INELgt0");
      registryMC.get<TH1>(HIST("Events/hEvtMcGen"))->GetXaxis()->SetBinLabel(4, "INELgt0TVX");
      registryMC.add("Events/hEvtMCRec", "Events/hEvtMCRec", {kTH1F, {{3, 0.f, 3.f}}});
      registryMC.get<TH1>(HIST("Events/hEvtMCRec"))->GetXaxis()->SetBinLabel(1, "all");
      registryMC.get<TH1>(HIST("Events/hEvtMCRec"))->GetXaxis()->SetBinLabel(2, "evt sel");
      registryMC.get<TH1>(HIST("Events/hEvtMCRec"))->GetXaxis()->SetBinLabel(3, "INELgt0");
      registryMC.add("Events/hEvtMcGenColls", "Number of events; Cut; #Events Passed Cut", {kTH1F, {{4, 0.5, 4.5}}});
      registryMC.get<TH1>(HIST("Events/hEvtMcGenColls"))->GetXaxis()->SetBinLabel(1, "Gen. coll");
      registryMC.get<TH1>(HIST("Events/hEvtMcGenColls"))->GetXaxis()->SetBinLabel(2, "At least 1 reco");
      registryMC.get<TH1>(HIST("Events/hEvtMcGenColls"))->GetXaxis()->SetBinLabel(3, "Reco. coll.");
      registryMC.get<TH1>(HIST("Events/hEvtMcGenColls"))->GetXaxis()->SetBinLabel(4, "Reco. good coll.");
      //
      registryMC.add("Events/hNchGenVsCent", "Gen Nch vs Cent; mult; Gen Nch (|#eta|<0.8)", {kTH2F, {nChAxis, multAxis}});
      registryMC.add("Events/hVtxZRec", "MC Rec vertex z position", kTH1F, {vtxzAxis});
      registryMC.add("Events/hVtxZGen", "Generated vertex z position", kTH1F, {vtxzAxis});
      registryMC.add("Events/hNchTVX", "Nch in FT0A+FT0C; Nch; status", {kTH2F, {nChAxis, {2, 0, 2}}});
      registryMC.add({"Tracks/hPtFakes", "Fake tracks; Gen Nch (|#eta|<0.8); flat; #it{p}_{T} (GeV/#it{c})", {kTHnSparseF, {nChAxis, flatAxis, ptAxis}}});
      // Event loss
      registryMC.add("Events/hNchVsFlatGenINELgt0", "Gen Nch w/o Evt sel; Gen Nch (|#eta|<0.8); flat", {kTH2F, {nChAxis, flatAxis}});
      registryMC.add("Events/hNchVsFlatGenINELgt0wRecEvtSel", "Gen Nch w/ Nrec > 0 + Evt sel; Gen Nch (|#eta|<0.8); flat", {kTH2F, {nChAxis, flatAxis}});
      // Event split
      registryMC.add("Events/hCentVsFlatRecINELgt0", "Gen evt w/o Evt sel; mult; flat", {kTH2F, {multAxis, flatAxis}});
      registryMC.add("Events/hCentVsFlatRecINELgt0wRecEvt", "Gen evt w/ Nrec > 0; mult; flat", {kTH2F, {multAxis, flatAxis}});
      registryMC.add("Events/hCentVsFlatRecINELgt0wRecEvtSel", "Gen evt w/ Nrec > 0 + Evt sel; mult; flat", {kTH2F, {multAxis, flatAxis}});
      for (int i = 0; i < Npart; ++i) {
        // Signal loss
        registryMC.add({fmt::format(CpTgenPrimSgnF.data(), CspeciesAll[i]).c_str(), "Gen evt w/o Evt sel; Gen Nch (|#eta|<0.8); flat; #it{p}_{T} (GeV/#it{c})", {kTHnSparseF, {nChAxis, flatAxis, ptAxis}}});
        registryMC.add({fmt::format(CpTrecCollPrimSgnF.data(), CspeciesAll[i]).c_str(), "Gen Nch w/ Nrec > 0; Gen Nch (|#eta|<0.8); flat; #it{p}_{T} (GeV/#it{c})", {kTHnSparseF, {nChAxis, flatAxis, ptAxis}}});
        // Closure test
        registryMC.add({fmt::format(CpTmcClosureGenPrimF.data(), CspeciesAll[i]).c_str(), "Gen evt w/o Evt sel; Gen Nch (|#eta|<0.8); flat; #it{p}_{T} (GeV/#it{c})", {kTHnSparseF, {nChAxis, flatAxis, ptAxis}}});
        registryMC.add({fmt::format(CpTmcClosureRecF.data(), CspeciesAll[i]).c_str(), "Gen Nch w/ Nrec > 0 + Evt. sel; Gen Nch (|#eta|<0.8); flat; #it{p}_{T} (GeV/#it{c})", {kTHnSparseF, {nChAxis, flatAxis, ptAxis}}});
      }
      if (cfgFillNclVsPhiCutQaHist) {
        registryMC.add("Tracks/postSel/hPtPhi", "; #it{p}_{T} (GeV/#it{c}); fmod(#varphi,#pi/9)", {kTH2F, {ptAxis, phiAxisMod}});
        registryMC.add("Tracks/postSel/hPtPhiNclTPC", "; #eta; #it{p}_{T} (GeV/#it{c}); fmod(#varphi,#pi/9); N_{cluster}", {kTHnSparseF, {etaAxis, ptAxis, phiAxisMod, clTpcAxis}});
        registryMC.add("Tracks/postSel/hPtPhiNclPIDTPC", "; #eta; #it{p}_{T} (GeV/#it{c}); fmod(#varphi,#pi/9); N_{PID cluster}", {kTHnSparseF, {etaAxis, ptAxis, phiAxisMod, clTpcAxis}});
        registryMC.add("Tracks/postSel/hPtNclTPC", "; #it{p}_{T} (GeV/#it{c}); N_{cluster}", {kTH2F, {ptAxis, clTpcAxis}});
        registryMC.add("Tracks/postSel/pPtNclTPC", "; #it{p}_{T} (GeV/#it{c}); N_{cluster}", {kTProfile, {ptAxis}});
        registryMC.add("Tracks/postSel/hPtNclPIDTPC", "; #it{p}_{T} (GeV/#it{c}); N_{PID cluster}", {kTH2F, {ptAxis, clTpcAxis}});
        registryMC.add("Tracks/postSel/pPtNclPIDTPC", "; #it{p}_{T} (GeV/#it{c}); N_{PID cluster}", {kTProfile, {ptAxis}});
      }
      registryMC.addClone("Tracks/postSel/", "Tracks/preSel/");

      for (int i = 0; i < NpartChrg; i++) {
        const std::string strID = Form("/%s/%s", (i < Npart) ? "pos" : "neg", Pid[i % Npart]);
        hPtGenRecEvt[i] = registryMC.add<THnSparse>("Tracks/hPtGenRecEvt" + strID, "Gen evt w/ Nrec > 0; mult; flat; #it{p}_{T} (GeV/#it{c})", kTHnSparseF, {multAxis, flatAxis, ptAxis});
        hPtGenPrimRecEvt[i] = registryMC.add<THnSparse>("Tracks/hPtGenPrimRecEvt" + strID, "Gen evt w/ Nrec > 0 (primary); mult; flat; #it{p}_{T} (GeV/#it{c})", kTHnSparseF, {multAxis, flatAxis, ptAxis});
        hPtGenRecEvtGtZero[i] = registryMC.add<THnSparse>("Tracks/hPtGenRecEvtGtZero" + strID, "Gen evt w/ Nrec > 0; mult; flat; #it{p}_{T} (GeV/#it{c})", kTHnSparseF, {multAxis, flatAxis, ptAxis});
        hPtGenPrimRecEvtGtZero[i] = registryMC.add<THnSparse>("Tracks/hPtGenPrimRecEvtGtZero" + strID, "Gen evt w/ Nrec > 0 (primary); mult; flat; #it{p}_{T} (GeV/#it{c})", kTHnSparseF, {multAxis, flatAxis, ptAxis});
        hPtEffGenPrim[i] = registryMC.add<THnSparse>("Tracks/hPtEffGenPrim" + strID, "Gen evt w/o rec Evt; mult; flat; #it{p}_{T} (GeV/#it{c})", kTHnSparseF, {multAxis, flatAxis, ptAxis});
        hPtEffGenWeak[i] = registryMC.add<THnSparse>("Tracks/hPtEffGenWeak" + strID, "Gen evt w/o rec Evt; mult; flat; #it{p}_{T} (GeV/#it{c})", kTHnSparseF, {multAxis, flatAxis, ptAxis});
        hPtEffGenMat[i] = registryMC.add<THnSparse>("Tracks/hPtEffGenMat" + strID, "Gen evt w/o rec Evt; mult; flat; #it{p}_{T} (GeV/#it{c})", kTHnSparseF, {multAxis, flatAxis, ptAxis});
        hPtEffGenPrimEvtSelGen[i] = registryMC.add<THnSparse>("Tracks/hPtEffGenPrimEvtSelGen" + strID, "Gen evt w/o rec Evt; mult; flat; #it{p}_{T} (GeV/#it{c})", kTHnSparseF, {multAxis, flatAxis, ptAxis});
        hPtEffGenWeakEvtSelGen[i] = registryMC.add<THnSparse>("Tracks/hPtEffGenWeakEvtSelGen" + strID, "Gen evt w/o rec Evt; mult; flat; #it{p}_{T} (GeV/#it{c})", kTHnSparseF, {multAxis, flatAxis, ptAxis});
        hPtEffGenMatEvtSelGen[i] = registryMC.add<THnSparse>("Tracks/hPtEffGenMatEvtSelGen" + strID, "Gen evt w/o rec Evt; mult; flat; #it{p}_{T} (GeV/#it{c})", kTHnSparseF, {multAxis, flatAxis, ptAxis});
        hPtEffRecGoodCollPrim[i] = registryMC.add<THnSparse>("Tracks/hPtEffRecGoodCollPrim" + strID, "; mult; flat; #it{p}_{T} (GeV/#it{c})", kTHnSparseF, {multAxis, flatAxis, ptAxis});
        hPtEffRecGoodCollWeak[i] = registryMC.add<THnSparse>("Tracks/hPtEffRecGoodCollWeak" + strID, "; mult; flat; #it{p}_{T} (GeV/#it{c})", kTHnSparseF, {multAxis, flatAxis, ptAxis});
        hPtEffRecGoodCollMat[i] = registryMC.add<THnSparse>("Tracks/hPtEffRecGoodCollMat" + strID, "; mult; flat; #it{p}_{T} (GeV/#it{c})", kTHnSparseF, {multAxis, flatAxis, ptAxis});
        hDCAxyRecBadCollPrim[i] = registryMC.add<TH2>("Tracks/hDCAxyRecBadCollPrim" + strID, "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {ptAxis, dcaXYAxis});
        hDCAxyRecBadCollWeak[i] = registryMC.add<TH2>("Tracks/hDCAxyRecBadCollWeak" + strID, "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {ptAxis, dcaXYAxis});
        hDCAxyRecBadCollMat[i] = registryMC.add<TH2>("Tracks/hDCAxyRecBadCollMat" + strID, "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {ptAxis, dcaXYAxis});
        hPtVsDCAxyRecGoodCollPrim[i] = registryMC.add<TH2>("Tracks/hPtVsDCAxyRecGoodCollPrim" + strID, "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {ptAxis, dcaXYAxis});
        hPtVsDCAxyRecGoodCollWeak[i] = registryMC.add<TH2>("Tracks/hPtVsDCAxyRecGoodCollWeak" + strID, "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {ptAxis, dcaXYAxis});
        hPtVsDCAxyRecGoodCollMat[i] = registryMC.add<TH2>("Tracks/hPtVsDCAxyRecGoodCollMat" + strID, "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {ptAxis, dcaXYAxis});
      }

      for (int i = 0; i < Npart; i++) {
        registryMC.add({fmt::format(CpTeffGenPrimRecEvtF.data(), CspeciesAll[i]).c_str(), "Gen evt w/ Nrec > 0; mult; flat; #it{p}_{T} (GeV/#it{c})", {kTHnSparseF, {multAxis, flatAxis, ptAxis}}});
        registryMC.add({fmt::format(CpTeffPrimRecEvtF.data(), CspeciesAll[i]).c_str(), "Gen evt w/ Nrec > 0 + Evt sel; mult; flat; #it{p}_{T} (GeV/#it{c})", {kTHnSparseF, {multAxis, flatAxis, ptAxis}}});
        registryMC.add({fmt::format(CpTvsDCAxyAllF.data(), CspeciesAll[i]).c_str(), "; mult; flat; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", {kTHnSparseF, {multAxis, flatAxis, ptAxis, dcaXYAxis}}});
        registryMC.add({fmt::format(CpTvsDCAxyPrimAllF.data(), CspeciesAll[i]).c_str(), "; mult; flat; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", {kTHnSparseF, {multAxis, flatAxis, ptAxis, dcaXYAxis}}});
        registryMC.add({fmt::format(CpTvsDCAxyWeakAllF.data(), CspeciesAll[i]).c_str(), "; mult; flat; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", {kTHnSparseF, {multAxis, flatAxis, ptAxis, dcaXYAxis}}});
        registryMC.add({fmt::format(CpTvsDCAxyMatAllF.data(), CspeciesAll[i]).c_str(), "; mult; flat; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", {kTHnSparseF, {multAxis, flatAxis, ptAxis, dcaXYAxis}}});
        registryMC.add({fmt::format(CdEdxMcRecPrimF.data(), CspeciesAll[i]).c_str(), "; #eta; mult; flat; #it{p} (GeV/#it{c}); dEdx", {kTHnSparseF, {etaAxis, multAxis, flatAxis, pAxis, dEdxAxis}}});
        registryMC.add({fmt::format(CdEdxMcRecPrimSelF.data(), CspeciesAll[i]).c_str(), "; #eta; mult; flat; #it{p} (GeV/#it{c}); dEdx", {kTHnSparseF, {etaAxis, multAxis, flatAxis, pAxis, dEdxAxis}}});
      }

      // Hash list for efficiency
      listEfficiency.setObject(new THashList);
      static_for<0, 1>([&](auto pidSgn) {
        bookMcHist<pidSgn, o2::track::PID::Pion>();
        bookMcHist<pidSgn, o2::track::PID::Kaon>();
        bookMcHist<pidSgn, o2::track::PID::Proton>();
        initEfficiency<pidSgn, o2::track::PID::Pion>();
        initEfficiency<pidSgn, o2::track::PID::Kaon>();
        initEfficiency<pidSgn, o2::track::PID::Proton>();
      });

      LOG(info) << "Size of the MC histograms:";
      registryMC.print();
    }
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    fv0AmplCorr.clear();
    fv0AmplCorr = {};
    std::string fullPathCalibGain;
    std::string fullPathCalibVtx;
    std::string fullPathCalibDeDxMip;
    std::string fullPathCalibDeDxPlateau;
    auto timestamp = bc.timestamp();
    auto runnumber = bc.runNumber();

    if (trkSelOpt.cfgRejectTrkAtTPCSector || v0SelOpt.cfgRejectV0sAtTPCSector) {
      grpmag = ccdb->getForRun<o2::parameters::GRPMagField>(ccdbConf.grpmagPath, runnumber);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << ccdbConf.grpmagPath << " of object GRPMagField and " << ccdbConf.grpPath << " of object GRPObject for run " << runnumber;
      }
      magField = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for run " << runnumber << " with magnetic field of " << magField << " kZG";
    }
    if (applyCalibGain) {
      fullPathCalibGain = cfgGainEqCcdbPath;
      fullPathCalibGain += "/FV0";
      const auto* objfv0Gain = getForTsOrRun<std::vector<float>>(fullPathCalibGain, timestamp, runnumber);
      if (!objfv0Gain) {
        for (auto i{0u}; i < CnCellsFV0; i++) {
          fv0AmplCorr.push_back(1.);
          LOGF(warning, "Setting FV0 calibration object values to 1");
        }
      } else {
        fv0AmplCorr = *(objfv0Gain);
      }
    }
    if (applyCalibVtx) {
      fullPathCalibVtx = cfgVtxEqCcdbPath;
      fullPathCalibVtx += "/FV0";
      zVtxMap = getForTsOrRun<TProfile2D>(fullPathCalibVtx, timestamp, runnumber);
    }

    if (applyCalibDeDxFromCCDB) {
      fullPathCalibDeDxMip = cfgDeDxCalibCcdbPath;
      fullPathCalibDeDxMip += "/MIP";
      fullPathCalibDeDxPlateau = cfgDeDxCalibCcdbPath;
      fullPathCalibDeDxPlateau += "/Plateau";
      if (!fullPathCalibDeDxMip.empty()) {
        dedxcalib.lCalibObjects = getForTsOrRun<TList>(fullPathCalibDeDxMip, timestamp, runnumber);
        if (dedxcalib.lCalibObjects) {
          LOG(info) << "CCDB objects loaded successfully";
          dedxcalib.hMIPcalibPos = static_cast<TH1F*>(dedxcalib.lCalibObjects->FindObject("hMIPcalibPos"));
          dedxcalib.hMIPcalibNeg = static_cast<TH1F*>(dedxcalib.lCalibObjects->FindObject("hMIPcalibNeg"));
          dedxcalib.hMIPcalibAll = static_cast<TH1F*>(dedxcalib.lCalibObjects->FindObject("hMIPcalibAll"));
          dedxcalib.lCalibLoaded = true;
          if (!dedxcalib.hMIPcalibPos || !dedxcalib.hMIPcalibNeg || !dedxcalib.hMIPcalibAll) {
            LOGF(error, "Problem loading CCDB objects! Please check");
            dedxcalib.lCalibLoaded = false;
          }
        } else {
          LOGF(fatal, "Could not load hMIPcalib from %s", fullPathCalibDeDxMip.c_str());
          dedxcalib.lCalibLoaded = false;
        }
      }
      if (!fullPathCalibDeDxPlateau.empty()) {
        dedxcalib.lCalibObjects = getForTsOrRun<TList>(fullPathCalibDeDxPlateau, timestamp, runnumber);
        if (dedxcalib.lCalibObjects) {
          LOG(info) << "CCDB objects loaded successfully";
          dedxcalib.hPlateauCalibPos = static_cast<TH1F*>(dedxcalib.lCalibObjects->FindObject("hPlateauCalibPos"));
          dedxcalib.hPlateauCalibNeg = static_cast<TH1F*>(dedxcalib.lCalibObjects->FindObject("hPlateauCalibNeg"));
          dedxcalib.hPlateauCalibAll = static_cast<TH1F*>(dedxcalib.lCalibObjects->FindObject("hPlateauCalibAll"));
          dedxcalib.lCalibLoaded = true;
          if (!dedxcalib.hPlateauCalibPos || !dedxcalib.hPlateauCalibNeg || !dedxcalib.hPlateauCalibAll) {
            LOGF(error, "Problem loading CCDB objects! Please check");
            dedxcalib.lCalibLoaded = false;
          }
        } else {
          LOGF(fatal, "Could not load hPlateauCalib from %s", fullPathCalibDeDxPlateau.c_str());
        }
      }
    }
  }

  template <typename T>
  std::unique_ptr<TF1> setFuncPars(T const& vecPars)
  {
    std::unique_ptr<TF1> fCalibDeDxFunc(new TF1("fCalibDeDxFunc", cfgCalibDeDxFunction.value.c_str(), -1., 1.));
    if (vecPars.size() >= 1) {
      for (typename T::size_type i = 0; i < vecPars.size(); i++) {
        fCalibDeDxFunc->SetParameter(i, vecPars[i]);
      }
    }
    return fCalibDeDxFunc;
  }

  template <int pidSgn, o2::track::PID::ID id, typename P>
  bool isPID(const P& mcParticle)
  {
    static_assert(pidSgn == CnullInt || pidSgn == ConeInt);
    static_assert(id > CnullInt || id < Npart);
    constexpr int Cidx = id + pidSgn * Npart;
    return mcParticle.pdgCode() == PidSgn[Cidx];
  }

  template <ChargeType chrg, typename T>
  bool selTOFPi(T const& track)
  {
    if (track.hasTOF() && track.goodTOFMatch()) {
      const float tTOF = track.tofSignal();
      const float trkLength = track.length();
      const float tExpPiTOF = track.tofExpSignalPi(tTOF);
      if (track.p() >= trkSelOpt.cfgMomSelPiTOF && trkLength > Cnull && tTOF > Cnull) {
        registryData.fill(HIST(CprefixCleanTof) + HIST(Ccharge[chrg]) + HIST("hTofExpPi"), track.p(), tExpPiTOF / tTOF);
        if (std::abs((tExpPiTOF / tTOF) - Cone) < trkSelOpt.cfgTofBetaPiMax) {
          registryData.fill(HIST(CprefixCleanTof) + HIST(Ccharge[chrg]) + HIST("hBetaVsP"), track.p(), track.beta());
          // if (std::abs(track.tpcNSigmaPi()) < v0SelOpt.cfgNsigmaPiTPC && std::abs(track.tofNSigmaPi()) < v0SelOpt.cfgNsigmaPiTOF) {
          return true;
          // }
        }
      }
    }
    return false;
  }

  template <int id, typename T, typename C>
  void fillDCA(T const& tracks, C const& collision, aod::BCsWithTimestamps const& /*bcs*/)
  {
    if (trkSelOpt.cfgRejectTrkAtTPCSector) {
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      int currentRun = bc.runNumber();
      if (runNumber != currentRun) {
        initCCDB(bc);
        runNumber = currentRun;
      }
    }
    const float mult = getMult(collision);
    const float flat = fillFlat<false>(collision);
    for (const auto& track : tracks) {
      if (std::abs(track.eta()) > trkSelOpt.cfgTrkEtaMax) {
        continue;
      }
      if (track.pt() < trkSelOpt.cfgTrkPtMin) {
        continue;
      }
      if (trkSelOpt.cfgApplyNcl && track.tpcNClsFound() < trkSelOpt.cfgNclTPCMin) {
        continue;
      }
      if (trkSelOpt.cfgApplyNclPID && track.tpcNClsPID() < trkSelOpt.cfgNclPidTPCMin) {
        continue;
      }
      float phiModn = track.phi();
      phiMod(phiModn, magField, track.sign());
      if (trkSelOpt.cfgRejectTrkAtTPCSector && (track.pt() >= trkSelOpt.cfgPhiCutPtMin && phiModn < fPhiCutHigh->Eval(track.pt()) && phiModn > fPhiCutLow->Eval(track.pt()))) {
        continue;
      }
      if (!isDCAxyWoCut(track)) {
        continue;
      }
      if (track.hasTOF() && (std::sqrt(std::pow(std::fabs(o2::aod::pidutils::tpcNSigma<id>(track)), 2) + std::pow(std::fabs(o2::aod::pidutils::tofNSigma<id>(track)), 2) < trkSelOpt.cfgDcaNsigmaCombinedMax))) {
        registryData.fill(HIST(Cprefix) + HIST(CspeciesAll[id]) + HIST(CpTvsDCAxy), mult, flat, track.pt(), track.dcaXY());
      }
    }
  }

  template <typename T, typename V, typename C>
  void filldEdx(T const& tracks, V const& v0s, C const& collision, aod::BCsWithTimestamps const& bcs)
  {
    if (trkSelOpt.cfgRejectTrkAtTPCSector || v0SelOpt.cfgRejectV0sAtTPCSector || applyCalibGain || applyCalibVtx) {
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      int currentRun = bc.runNumber();
      if (runNumber != currentRun) {
        initCCDB(bc);
        runNumber = currentRun;
      }
    }

    const float mult = getMult(collision);
    const float flat = fillFlat<true>(collision);
    registryData.fill(HIST("Events/hFlatVsMultEst"), flat, mult);

    countTracks<true>(tracks, collision, bcs, mult);

    for (const auto& track : tracks) {
      float dEdx = track.tpcSignal();
      if (cfgFillTrackQaHist) {
        fillTrackQA<kBefore, true>(track);
      }
      if (!isGoodTrack(track, magField)) {
        continue;
      }
      if (cfgFillTrackQaHist) {
        fillTrackQA<kAfter, true>(track);
      }
      if (cfgFilldEdxCalibHist && cfgFilldEdxQaHist) {
        if (cfgFillChrgType) {
          if (track.sign() * track.p() > Cnull) {
            filldEdxQA<kPos, kBefore, true>(track, collision, dEdx);
          } else {
            filldEdxQA<kNeg, kBefore, true>(track, collision, dEdx);
          }
        } else {
          filldEdxQA<kAll, kBefore, true>(track, collision, dEdx);
        }
      }
      if (applyCalibDeDx) {
        if (cfgFillChrgType) {
          if (track.sign() * track.p() > Cnull) {
            if (applyCalibDeDxFromCCDB) {
              dEdx *= (50.0 / dedxcalib.hMIPcalibPos->GetBinContent(dedxcalib.hMIPcalibPos->FindBin(track.eta())));
            } else {
              dEdx *= (50.0 / getCalibration(fDeDxVsEta, track));
            }
            if (cfgFilldEdxQaHist) {
              filldEdxQA<kPos, kAfter, true>(track, collision, dEdx);
            }
          } else {
            if (applyCalibDeDxFromCCDB) {
              dEdx *= (50.0 / dedxcalib.hMIPcalibNeg->GetBinContent(dedxcalib.hMIPcalibNeg->FindBin(track.eta())));
            } else {
              dEdx *= (50.0 / getCalibration(fDeDxVsEta, track));
            }
            if (cfgFilldEdxQaHist) {
              filldEdxQA<kNeg, kAfter, true>(track, collision, dEdx);
            }
          }
        } else {
          if (applyCalibDeDxFromCCDB) {
            dEdx *= (50.0 / dedxcalib.hMIPcalibAll->GetBinContent(dedxcalib.hMIPcalibAll->FindBin(track.eta())));
          } else {
            dEdx *= (50.0 / getCalibration<false>(fDeDxVsEta, track));
          }
          if (cfgFilldEdxQaHist) {
            filldEdxQA<kAll, kAfter, true>(track, collision, dEdx);
          }
        }
      }

      // PID TPC dEdx
      if (cfgFillChrgType) {
        if (track.sign() * track.p() > Cnull) {
          registryData.fill(HIST(Cprefix) + HIST(Ccharge[kPos]) + HIST("hdEdx"), track.eta(), mult, flat, track.p(), dEdx);
        } else {
          registryData.fill(HIST(Cprefix) + HIST(Ccharge[kNeg]) + HIST("hdEdx"), track.eta(), mult, flat, track.p(), dEdx);
        }
      } else {
        registryData.fill(HIST(Cprefix) + HIST(Ccharge[kAll]) + HIST("hdEdx"), track.eta(), mult, flat, track.p(), dEdx);
      }

      // TOF pions
      if (cfgFillV0Hist) {
        if (selTOFPi<kAll>(track)) {
          if (cfgFillChrgType) {
            if (track.sign() * track.p() > Cnull) {
              if (cfgStoreThnSparse) {
                registryData.fill(HIST(CprefixCleanTof) + HIST(Ccharge[kPos]) + HIST("hPiTof"), track.eta(), mult, flat, track.p(), dEdx);
              } else {
                registryData.fill(HIST(CprefixCleanTof) + HIST(Ccharge[kPos]) + HIST("hPiTof"), track.eta(), track.p(), dEdx);
              }
            } else {
              if (cfgStoreThnSparse) {
                registryData.fill(HIST(CprefixCleanTof) + HIST(Ccharge[kNeg]) + HIST("hPiTof"), track.eta(), mult, flat, track.p(), dEdx);
              } else {
                registryData.fill(HIST(CprefixCleanTof) + HIST(Ccharge[kNeg]) + HIST("hPiTof"), track.eta(), track.p(), dEdx);
              }
            }
          } else {
            if (cfgStoreThnSparse) {
              registryData.fill(HIST(CprefixCleanTof) + HIST(Ccharge[kAll]) + HIST("hPiTof"), track.eta(), mult, flat, track.p(), dEdx);
            } else {
              registryData.fill(HIST(CprefixCleanTof) + HIST(Ccharge[kAll]) + HIST("hPiTof"), track.eta(), track.p(), dEdx);
            }
          }
        }
      }
    }

    // V0s
    if (cfgFillV0Hist) {
      for (const auto& v0 : v0s) {
        if (v0.v0Type() != v0SelOpt.cfgV0TypeSel && v0SelOpt.cfgV0TypeSel > -1) {
          continue;
        }
        if (!isGoodV0Track(v0, tracks, magField)) {
          continue;
        }

        const auto& posTrack = v0.template posTrack_as<T>();
        const auto& negTrack = v0.template negTrack_as<T>();
        float dEdxPos = posTrack.tpcSignal();
        float dEdxNeg = negTrack.tpcSignal();

        if (applyCalibDeDx) {
          if (cfgFillChrgTypeV0s) {
            if (applyCalibDeDxFromCCDB) {
              dEdxPos *= (50.0 / dedxcalib.hMIPcalibPos->GetBinContent(dedxcalib.hMIPcalibPos->FindBin(posTrack.eta())));
              dEdxNeg *= (50.0 / dedxcalib.hMIPcalibNeg->GetBinContent(dedxcalib.hMIPcalibNeg->FindBin(negTrack.eta())));
            } else {
              dEdxPos *= (50.0 / getCalibration(fDeDxVsEta, posTrack));
              dEdxNeg *= (50.0 / getCalibration(fDeDxVsEta, negTrack));
            }
          } else {
            if (applyCalibDeDxFromCCDB) {
              dEdxPos *= (50.0 / dedxcalib.hMIPcalibAll->GetBinContent(dedxcalib.hMIPcalibAll->FindBin(posTrack.eta())));
              dEdxNeg *= (50.0 / dedxcalib.hMIPcalibAll->GetBinContent(dedxcalib.hMIPcalibAll->FindBin(negTrack.eta())));
            } else {
              dEdxPos *= (50.0 / getCalibration<false>(fDeDxVsEta, posTrack));
              dEdxNeg *= (50.0 / getCalibration<false>(fDeDxVsEta, negTrack));
            }
          }
        }

        if (selectTypeV0s(collision, v0, posTrack, negTrack) == kGa) { // Gamma selection
          fillV0QA<kEl, kGa>(v0, posTrack);
          fillV0QA<kEl, kGa>(v0, negTrack);
          if (applyCalibDeDx) {
            if (cfgFillChrgTypeV0s) {
              if (applyCalibDeDxFromCCDB) {
                const float dEdxPosGa = dedxcalib.hMIPcalibPos->GetBinContent(dedxcalib.hMIPcalibPos->FindBin(posTrack.eta()));
                const float dEdxNegGa = dedxcalib.hMIPcalibNeg->GetBinContent(dedxcalib.hMIPcalibNeg->FindBin(negTrack.eta()));
                if (std::abs(dEdxPos - dEdxPosGa) >= v0SelOpt.cfgdEdxPlateauSel || std::abs(dEdxNeg - dEdxNegGa) >= v0SelOpt.cfgdEdxPlateauSel) {
                  continue;
                }
              } else {
                const float dEdxPosGa = getCalibration(fEDeDxVsEta, posTrack);
                const float dEdxNegGa = getCalibration(fEDeDxVsEta, negTrack);
                if (std::abs(dEdxPos - dEdxPosGa) >= v0SelOpt.cfgdEdxPlateauSel || std::abs(dEdxNeg - dEdxNegGa) >= v0SelOpt.cfgdEdxPlateauSel) {
                  continue;
                }
              }
            } else {
              if (applyCalibDeDxFromCCDB) {
                const float dEdxPosGa = dedxcalib.hPlateauCalibAll->GetBinContent(dedxcalib.hPlateauCalibAll->FindBin(posTrack.eta()));
                const float dEdxNegGa = dedxcalib.hPlateauCalibAll->GetBinContent(dedxcalib.hPlateauCalibAll->FindBin(negTrack.eta()));
                if (std::abs(dEdxPos - dEdxPosGa) >= v0SelOpt.cfgdEdxPlateauSel || std::abs(dEdxNeg - dEdxNegGa) >= v0SelOpt.cfgdEdxPlateauSel) {
                  continue;
                }
              } else {
                const float dEdxPosGa = getCalibration<false>(fEDeDxVsEta, posTrack);
                const float dEdxNegGa = getCalibration<false>(fEDeDxVsEta, negTrack);
                if (std::abs(dEdxPos - dEdxPosGa) >= v0SelOpt.cfgdEdxPlateauSel || std::abs(dEdxNeg - dEdxNegGa) >= v0SelOpt.cfgdEdxPlateauSel) {
                  continue;
                }
              }
            }
          }
          if (cfgStoreThnSparse) {
            if (cfgFillChrgType) {
              registryData.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kPos]) + HIST("hEV0"), posTrack.eta(), mult, flat, posTrack.sign() * posTrack.p(), dEdxPos);
              registryData.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kNeg]) + HIST("hEV0"), negTrack.eta(), mult, flat, negTrack.sign() * negTrack.p(), dEdxNeg);
            } else {
              registryData.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kAll]) + HIST("hEV0"), posTrack.eta(), mult, flat, posTrack.sign() * posTrack.p(), dEdxPos);
              registryData.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kAll]) + HIST("hEV0"), negTrack.eta(), mult, flat, negTrack.sign() * negTrack.p(), dEdxNeg);
            }
          } else {
            if (cfgFillChrgType) {
              registryData.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kPos]) + HIST("hEV0"), posTrack.eta(), posTrack.sign() * posTrack.p(), dEdxPos);
              registryData.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kNeg]) + HIST("hEV0"), negTrack.eta(), negTrack.sign() * negTrack.p(), dEdxNeg);
            } else {
              registryData.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kAll]) + HIST("hEV0"), posTrack.eta(), posTrack.sign() * posTrack.p(), dEdxPos);
              registryData.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kAll]) + HIST("hEV0"), negTrack.eta(), negTrack.sign() * negTrack.p(), dEdxNeg);
            }
          }
        }
        if (selectTypeV0s(collision, v0, posTrack, negTrack) == kKz) { // K0S -> pi + pi
          fillV0QA<kPi, kKz>(v0, posTrack);
          fillV0QA<kPi, kKz>(v0, negTrack);
          if (cfgStoreThnSparse) {
            if (cfgFillChrgType) {
              registryData.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kPos]) + HIST("hPiV0"), posTrack.eta(), mult, flat, posTrack.sign() * posTrack.p(), dEdxPos);
              registryData.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kNeg]) + HIST("hPiV0"), negTrack.eta(), mult, flat, negTrack.sign() * negTrack.p(), dEdxNeg);
            } else {
              registryData.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kAll]) + HIST("hPiV0"), posTrack.eta(), mult, flat, posTrack.sign() * posTrack.p(), dEdxPos);
              registryData.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kAll]) + HIST("hPiV0"), negTrack.eta(), mult, flat, negTrack.sign() * negTrack.p(), dEdxNeg);
            }
          } else {
            if (cfgFillChrgType) {
              registryData.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kPos]) + HIST("hPiV0"), posTrack.eta(), posTrack.sign() * posTrack.p(), dEdxPos);
              registryData.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kNeg]) + HIST("hPiV0"), negTrack.eta(), negTrack.sign() * negTrack.p(), dEdxNeg);
            } else {
              registryData.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kAll]) + HIST("hPiV0"), posTrack.eta(), posTrack.sign() * posTrack.p(), dEdxPos);
              registryData.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kAll]) + HIST("hPiV0"), negTrack.eta(), negTrack.sign() * negTrack.p(), dEdxNeg);
            }
          }
        }
        if (selectTypeV0s(collision, v0, posTrack, negTrack) == kLam) { // L -> p + pi-
          fillV0QA<kPi, kLam>(v0, negTrack);
          fillV0QA<kPr, kLam>(v0, posTrack);
          if (cfgStoreThnSparse) {
            if (cfgFillChrgType) {
              registryData.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kPos]) + HIST("hPV0"), posTrack.eta(), mult, flat, posTrack.sign() * posTrack.p(), dEdxPos);
              registryData.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kNeg]) + HIST("hPiV0"), negTrack.eta(), mult, flat, negTrack.sign() * negTrack.p(), dEdxNeg);
            } else {
              registryData.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kAll]) + HIST("hPV0"), posTrack.eta(), mult, flat, posTrack.sign() * posTrack.p(), dEdxPos);
              registryData.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kAll]) + HIST("hPiV0"), negTrack.eta(), mult, flat, negTrack.sign() * negTrack.p(), dEdxNeg);
            }
          } else {
            if (cfgFillChrgType) {
              registryData.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kPos]) + HIST("hPV0"), posTrack.eta(), posTrack.sign() * posTrack.p(), dEdxPos);
              registryData.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kNeg]) + HIST("hPiV0"), negTrack.eta(), negTrack.sign() * negTrack.p(), dEdxNeg);
            } else {
              registryData.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kAll]) + HIST("hPV0"), posTrack.eta(), posTrack.sign() * posTrack.p(), dEdxPos);
              registryData.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kAll]) + HIST("hPiV0"), negTrack.eta(), negTrack.sign() * negTrack.p(), dEdxNeg);
            }
          }
        }
        if (selectTypeV0s(collision, v0, posTrack, negTrack) == kaLam) { // antiLambda -> pbar + pi+
          fillV0QA<kPi, kaLam>(v0, posTrack);
          fillV0QA<kPr, kaLam>(v0, negTrack);
          if (cfgStoreThnSparse) {
            if (cfgFillChrgType) {
              registryData.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kPos]) + HIST("hPiV0"), posTrack.eta(), mult, flat, posTrack.sign() * posTrack.p(), dEdxPos);
              registryData.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kNeg]) + HIST("hPV0"), negTrack.eta(), mult, flat, negTrack.sign() * negTrack.p(), dEdxNeg);
            } else {
              registryData.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kAll]) + HIST("hPiV0"), posTrack.eta(), mult, flat, posTrack.sign() * posTrack.p(), dEdxPos);
              registryData.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kAll]) + HIST("hPV0"), negTrack.eta(), mult, flat, negTrack.sign() * negTrack.p(), dEdxNeg);
            }
          } else {
            if (cfgFillChrgType) {
              registryData.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kPos]) + HIST("hPiV0"), posTrack.eta(), posTrack.sign() * posTrack.p(), dEdxPos);
              registryData.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kNeg]) + HIST("hPV0"), negTrack.eta(), negTrack.sign() * negTrack.p(), dEdxNeg);
            } else {
              registryData.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kAll]) + HIST("hPiV0"), posTrack.eta(), posTrack.sign() * posTrack.p(), dEdxPos);
              registryData.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kAll]) + HIST("hPV0"), negTrack.eta(), negTrack.sign() * negTrack.p(), dEdxNeg);
            }
          }
        }
      }
    }
  }

  template <bool isChrg = true, typename T, typename V>
  float getCalibration(std::vector<std::unique_ptr<T>> const& fCalib, V const& track)
  {
    float valCalib = -1.;
    if constexpr (isChrg) {
      if (track.sign() * track.p() > Cnull) {
        if (track.eta() >= 0) {
          valCalib = fCalib.at(0)->Eval(track.eta());
        } else {
          valCalib = fCalib.at(3)->Eval(track.eta());
        }
      } else {
        if (track.eta() >= 0) {
          valCalib = fCalib.at(1)->Eval(track.eta());
        } else {
          valCalib = fCalib.at(4)->Eval(track.eta());
        }
      }
    } else {
      if (track.eta() >= 0) {
        valCalib = fCalib.at(2)->Eval(track.eta());
      } else {
        valCalib = fCalib.at(5)->Eval(track.eta());
      }
    }
    return valCalib;
  }

  bool isChrgParticle(int code)
  {
    auto p = pdg->GetParticle(code);
    auto charge = 0.;
    if (p != nullptr) {
      charge = p->Charge();
    }
    return std::abs(charge) >= CminCharge;
  }

  template <typename P>
  int countPart(P const& particles)
  {
    auto nCharged = 0;
    for (auto const& particle : particles) {
      if (!isChrgParticle(particle.pdgCode())) {
        continue;
      }
      if (!particle.isPhysicalPrimary()) {
        continue;
      }
      if (std::abs(particle.eta()) > trkSelOpt.cfgTrkEtaMax) {
        continue;
      }
      nCharged++;
    }
    return nCharged;
  }

  template <typename P>
  bool isInelGt0wTVX(P const& particles)
  {
    int nChrgMc = 0;
    int nChrgFT0A = 0;
    int nChrgFT0C = 0;
    for (auto const& particle : particles) {
      if (!isChrgParticle(particle.pdgCode())) {
        continue;
      }
      if (!particle.isPhysicalPrimary()) {
        continue;
      }
      // trigger TVX
      if (particle.eta() > CminAccFT0A && particle.eta() < CmaxAccFT0A) {
        nChrgFT0A++;
      }
      if (particle.eta() > CminAccFT0C && particle.eta() < CmaxAccFT0C) {
        nChrgFT0C++;
      }
      nChrgMc++;
    }
    if (nChrgFT0A == CnullInt || nChrgFT0C == CnullInt) {
      registryMC.fill(HIST("Events/hNchTVX"), nChrgMc, 0.5);
      return false;
    }
    registryMC.fill(HIST("Events/hNchTVX"), nChrgMc, 1.5);
    if (nChrgMc == CnullInt) {
      return false;
    }
    return true;
  }

  template <bool fillHis = false, typename T, typename C>
  int countTracks(T const& tracks, C const& collision, aod::BCsWithTimestamps const& /*bcs*/, float mult)
  {
    if (trkSelOpt.cfgRejectTrkAtTPCSector || v0SelOpt.cfgRejectV0sAtTPCSector || applyCalibGain || applyCalibVtx) {
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      int currentRun = bc.runNumber();
      if (runNumber != currentRun) {
        initCCDB(bc);
        runNumber = currentRun;
      }
    }

    auto nTrk = 0;
    for (auto const& track : tracks) {
      if (!isGoodTrack(track, magField)) {
        continue;
      }
      nTrk++;
    }
    if (fillHis) {
      registryQC.fill(HIST("Events/hNchVsCent"), nTrk, mult);
    }
    return nTrk;
  }

  void phiMod(float& phimodn, const int& mag, const int& charge)
  {
    if (mag < Cnull) // for negative polarity field
      phimodn = o2::constants::math::TwoPI - phimodn;
    if (charge < Cnull) // for negative charge
      phimodn = o2::constants::math::TwoPI - phimodn;
    if (phimodn < Cnull)
      LOGF(warning, "phi < Cnull: %g", phimodn);

    phimodn += o2::constants::math::PI / 18.0f; // to center gap in the middle
    phimodn = std::fmod(phimodn, o2::constants::math::PI / 9.0f);
  }

  template <FillType ft, typename T>
  inline void fillNclVsPhiCutQaHist(T const& track, const float phimodn)
  {
    if (cfgFillNclVsPhiCutQaHist) {
      registryData.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hPtPhi"), track.pt(), phimodn);
      registryData.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hPtPhiNclTPC"), track.eta(), track.pt(), phimodn, track.tpcNClsFound());
      registryData.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hPtPhiNclPIDTPC"), track.eta(), track.pt(), phimodn, track.tpcNClsPID());
      registryData.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hPtNclTPC"), track.pt(), track.tpcNClsFound());
      registryData.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("pPtNclTPC"), track.pt(), track.tpcNClsFound());
      registryData.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hPtNclPIDTPC"), track.pt(), track.tpcNClsPID());
      registryData.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("pPtNclPIDTPC"), track.pt(), track.tpcNClsPID());
    }
  }

  template <typename T>
  bool isGoodTrack(T const& track, const int magfield)
  {
    registryData.fill(HIST("Tracks/hTrkSel"), trkSelAll);
    if (std::abs(track.eta()) > trkSelOpt.cfgTrkEtaMax) {
      return false;
    }
    registryData.fill(HIST("Tracks/hTrkSel"), trkSelEta);
    if (track.pt() < trkSelOpt.cfgTrkPtMin) {
      return false;
    }
    registryData.fill(HIST("Tracks/hTrkSel"), trkSelPt);
    if (!isDCAxyCut(track)) {
      return false;
    }
    registryData.fill(HIST("Tracks/hTrkSel"), trkSelDCA);
    if (track.tpcNClsCrossedRows() < minNCrossedRowsTPC) {
      return false;
    }
    registryData.fill(HIST("Tracks/hTrkSel"), trkNRowsTPC);
    if (trkSelOpt.cfgApplyNcl && track.tpcNClsFound() < trkSelOpt.cfgNclTPCMin) {
      return false;
    }
    registryData.fill(HIST("Tracks/hTrkSel"), trkSelNClsFound);

    if (trkSelOpt.cfgApplyNclPID && track.tpcNClsPID() < trkSelOpt.cfgNclPidTPCMin) {
      return false;
    }
    registryData.fill(HIST("Tracks/hTrkSel"), trkSelNClsPID);

    float phimodn = track.phi();
    phiMod(phimodn, magfield, track.sign());
    if (cfgFillNclVsPhiCutQaHist) {
      fillNclVsPhiCutQaHist<kBefore>(track, phimodn);
    }
    if (trkSelOpt.cfgRejectTrkAtTPCSector && (track.pt() >= trkSelOpt.cfgPhiCutPtMin && phimodn < fPhiCutHigh->Eval(track.pt()) && phimodn > fPhiCutLow->Eval(track.pt()))) {
      return false;
    }
    if (cfgFillNclVsPhiCutQaHist) {
      fillNclVsPhiCutQaHist<kAfter>(track, phimodn);
    }
    registryData.fill(HIST("Tracks/hTrkSel"), trkSelTPCBndr);
    return true;
  }

  template <int id, int typeMother, typename V, typename U>
  void fillV0QA(V const& v0, U const& track)
  {
    registryData.fill(HIST(CprefixV0qa) + HIST(PidDir[id]) + HIST(V0Dir[typeMother]) + HIST("hArmPod"), v0.alpha(), v0.qtarm());
    registryData.fill(HIST(CprefixV0qa) + HIST(PidDir[id]) + HIST(V0Dir[typeMother]) + HIST("hNclVsEta"), track.eta(), track.tpcNClsPID());
    registryData.fill(HIST(CprefixV0qa) + HIST(PidDir[id]) + HIST(V0Dir[typeMother]) + HIST("hdEdxMIPVsEta"), track.eta(), track.tpcSignal());
  }

  template <typename C, typename T1, typename T2>
  int selectTypeV0s(C const& collision, T1 const& v0, T2 const& postrk, T2 const& negtrk)
  {
    const float dMassK0s = std::abs(v0.mK0Short() - o2::constants::physics::MassK0Short);
    const float dMassL = std::abs(v0.mLambda() - o2::constants::physics::MassLambda0);
    const float dMassAL = std::abs(v0.mAntiLambda() - o2::constants::physics::MassLambda0);
    const float dMassG = std::abs(v0.mGamma() - o2::constants::physics::MassGamma);

    bool isMassG = false;
    bool isMassK0s = false;
    bool isMassL = false;
    bool isMassAL = false;

    if (dMassK0s > v0SelOpt.cfgdmassK && dMassL > v0SelOpt.cfgdmassL && dMassAL > v0SelOpt.cfgdmassL && dMassG < v0SelOpt.cfgdmassG) {
      isMassG = true;
    }
    if (dMassK0s < v0SelOpt.cfgdmassK && dMassL > v0SelOpt.cfgdmassL && dMassAL > v0SelOpt.cfgdmassL && dMassG > v0SelOpt.cfgdmassG) {
      isMassK0s = true;
    }
    if (dMassK0s > v0SelOpt.cfgdmassK && dMassL < v0SelOpt.cfgdmassL && dMassG > v0SelOpt.cfgdmassG) {
      isMassL = true;
    }
    if (dMassK0s > v0SelOpt.cfgdmassK && dMassAL < v0SelOpt.cfgdmassL && dMassG > v0SelOpt.cfgdmassG) {
      isMassAL = true;
    }

    // Gamma selection
    const float yGamma = RecoDecay::y(std::array{v0.px(), v0.py(), v0.pz()}, o2::constants::physics::MassGamma);
    if (isMassG) {
      if (std::abs(yGamma) < v0SelOpt.cfgV0Ymax) {                                                             // rapidity cut
        if (std::abs(v0.alpha()) < v0SelOpt.cfgArmPodGammasalpha && v0.qtarm() < v0SelOpt.cfgArmPodGammasqT) { //
          if (postrk.hasTPC() && std::abs(postrk.tpcNSigmaEl()) < v0SelOpt.cfgNsigmaElTPC) {
            if (postrk.hasTOF() && std::abs(postrk.tofNSigmaEl()) < v0SelOpt.cfgNsigmaElTOF) {
              return kGa;
            }
          }
          if (negtrk.hasTPC() && std::abs(negtrk.tpcNSigmaEl()) < v0SelOpt.cfgNsigmaElTPC) {
            if (negtrk.hasTOF() && std::abs(negtrk.tofNSigmaEl()) < v0SelOpt.cfgNsigmaElTOF) {
              return kGa;
            }
          }
        }
      }
    }
    // K0S selection, K0S -> pi + pi
    if (isMassK0s) {
      if (std::abs(v0.yK0Short()) < v0SelOpt.cfgV0Ymax) {                                                                                                      // rapidity cut
        if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short < v0SelOpt.cfgcTauK0s) {             // ctau cut
          if (v0.v0cosPA() >= v0SelOpt.cfgCosPAK0s && v0.v0radius() >= v0SelOpt.cfgV0radiusK0s && v0.qtarm() * v0SelOpt.cfgArmPodK0s > std::abs(v0.alpha())) { //
            if (postrk.hasTPC() && std::abs(postrk.tpcNSigmaPi()) < v0SelOpt.cfgNsigmaPiTPC) {
              if (postrk.hasTOF() && std::abs(postrk.tofNSigmaPi()) < v0SelOpt.cfgNsigmaPiTOF) {
                return kKz;
              }
            }
            if (negtrk.hasTPC() && std::abs(negtrk.tpcNSigmaPi()) < v0SelOpt.cfgNsigmaPiTPC) {
              if (negtrk.hasTOF() && std::abs(negtrk.tofNSigmaPi()) < v0SelOpt.cfgNsigmaPiTOF) {
                return kKz;
              }
            }
          }
        }
      }
    }
    // Lambda selection, L -> p + pi-
    if (isMassL) {
      if (std::abs(v0.yLambda()) < v0SelOpt.cfgV0Ymax) {                                                                                              // rapidity cut
        if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0 < v0SelOpt.cfgcTauLambda) { // ctau cut
          if (v0.v0cosPA() >= v0SelOpt.cfgCosPALambda && v0.v0radius() >= v0SelOpt.cfgV0radiusLambda) {                                               //
            if (postrk.hasTPC() && std::abs(postrk.tpcNSigmaPr()) < v0SelOpt.cfgNsigmaPrTPC && negtrk.hasTPC() && std::abs(negtrk.tpcNSigmaPi()) < v0SelOpt.cfgNsigmaPiTPC) {
              if (postrk.hasTOF() && std::abs(postrk.tofNSigmaPr()) < v0SelOpt.cfgNsigmaPrTOF && negtrk.hasTOF() && std::abs(negtrk.tofNSigmaPi()) < v0SelOpt.cfgNsigmaPiTOF) {
                return kLam;
              }
            }
          }
        }
      }
    }
    // antiLambda -> pbar + pi+
    if (isMassAL) {
      if (std::abs(v0.yLambda()) < v0SelOpt.cfgV0Ymax) {                                                                                              // rapidity cut
        if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0 < v0SelOpt.cfgcTauLambda) { // ctau cut
          if (v0.v0cosPA() >= v0SelOpt.cfgCosPALambda && v0.v0radius() >= v0SelOpt.cfgV0radiusLambda) {                                               //
            if (postrk.hasTPC() && std::abs(postrk.tpcNSigmaPi()) < v0SelOpt.cfgNsigmaPiTPC && negtrk.hasTPC() && std::abs(negtrk.tpcNSigmaPr()) < v0SelOpt.cfgNsigmaPrTPC) {
              if (postrk.hasTOF() && std::abs(postrk.tofNSigmaPi()) < v0SelOpt.cfgNsigmaPiTOF && negtrk.hasTOF() && std::abs(negtrk.tofNSigmaPr()) < v0SelOpt.cfgNsigmaPrTOF) {
                return kaLam;
              }
            }
          }
        }
      }
    }
    return kNaN;
  }

  template <bool fillHist = true, typename T1, typename T2>
  bool isGoodV0Track(T1 const& v0, T2 const& /*track*/, const int magfield)
  {
    const auto& posTrack = v0.template posTrack_as<T2>();
    const auto& negTrack = v0.template negTrack_as<T2>();

    registryData.fill(HIST("Tracks/V0qa/hV0Sel"), v0SelAll);
    if (posTrack.sign() * negTrack.sign() > Cnull) {
      return false;
    }
    registryData.fill(HIST("Tracks/V0qa/hV0Sel"), v0SelRejectSameSign);

    float posTrackPhiModn = posTrack.phi();
    float negTrackPhiModn = negTrack.phi();
    phiMod(posTrackPhiModn, magfield, posTrack.sign());
    phiMod(negTrackPhiModn, magfield, negTrack.sign());
    if (v0SelOpt.cfgRejectV0sAtTPCSector) {
      if ((posTrack.pt() >= trkSelOpt.cfgPhiCutPtMin && posTrackPhiModn < fPhiCutHigh->Eval(posTrack.pt()) && posTrackPhiModn > fPhiCutLow->Eval(posTrack.pt())) && (negTrack.pt() >= trkSelOpt.cfgPhiCutPtMin && negTrackPhiModn < fPhiCutHigh->Eval(negTrack.pt()) && negTrackPhiModn > fPhiCutLow->Eval(negTrack.pt()))) {
        return false;
      }
    }
    registryData.fill(HIST("Tracks/V0qa/hV0Sel"), v0SelRejectV0sAtTPCSector);
    // V0 topological selections
    if (v0.v0cosPA() < v0SelOpt.cfgv0cospa) {
      return false;
    }
    registryData.fill(HIST("Tracks/V0qa/hV0Sel"), v0SelCosPA);
    if (v0.v0radius() < v0SelOpt.cfgv0Rmin || v0.v0radius() > v0SelOpt.cfgv0Rmax) {
      return false;
    }
    registryData.fill(HIST("Tracks/V0qa/hV0Sel"), v0SelV0radius);
    if (std::abs(v0.dcapostopv()) < v0SelOpt.cfgDCAposToPV || std::abs(v0.dcanegtopv()) < v0SelOpt.cfgDCAnegToPV) {
      return false;
    }
    registryData.fill(HIST("Tracks/V0qa/hV0Sel"), v0SelDCAposToPV);
    // selection of V0 daughters
    if (!(isGoodV0DaughterTrack(posTrack) && isGoodV0DaughterTrack(negTrack))) {
      return false;
    }
    registryData.fill(HIST("Tracks/V0qa/hV0Sel"), v0SelDaughters);
    if (v0.dcaV0daughters() > v0SelOpt.cfgDCAv0daughter) {
      return false;
    }
    registryData.fill(HIST("Tracks/V0qa/hV0Sel"), v0SelDCAv0daughter);
    if constexpr (fillHist) {
      registryData.fill(HIST("Tracks/V0qa/hV0Pt"), v0.pt());
      registryData.fill(HIST("Tracks/V0qa/hV0ArmPod"), v0.alpha(), v0.qtarm());
    }
    return true;
  }

  template <typename T>
  bool isGoodV0DaughterTrack(const T& track)
  {
    if (track.eta() < v0SelOpt.cfgV0etamin || track.eta() > v0SelOpt.cfgV0etamax) {
      return false;
    }
    if (!track.hasTPC()) {
      return false;
    }
    if (track.tpcNClsCrossedRows() < v0SelOpt.cfgminNCrossedRowsTPC) {
      return false;
    }
    if (track.tpcCrossedRowsOverFindableCls() < v0SelOpt.cfgminNCrossedRowsOverFindableClustersTPC) {
      return false;
    }
    if (v0SelOpt.cfgApplyV0sNclFound) {
      if (track.tpcNClsFound() < v0SelOpt.cfgV0NclTPCMin) {
        return false;
      }
    }
    if (v0SelOpt.cfgApplyV0sNclPID) {
      if (track.tpcNClsPID() < v0SelOpt.cfgV0NclPidTPCMin) {
        return false;
      }
    }
    if (track.tpcChi2NCl() > v0SelOpt.cfgmaxChi2PerClusterTPC) {
      return false;
    }
    if (v0SelOpt.cfgRequireITS && (!track.hasITS())) {
      return false;
    }
    if (v0SelOpt.cfgRequireITS && track.itsNCls() < v0SelOpt.cfgminITSnClusters) {
      return false;
    }
    if (track.itsChi2NCl() > v0SelOpt.cfgmaxChi2PerClusterITS) {
      return false;
    }
    return true;
  }

  template <FillType ft, bool fillHist = false, typename T>
  inline void fillTrackQA(T const& track)
  {
    if constexpr (fillHist) {
      registryData.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hPt"), track.pt());
      registryData.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hPhi"), track.phi());
      registryData.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hEta"), track.eta());
      registryData.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hDCAXYvsPt"), track.pt(), track.dcaXY());
      registryData.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hDCAZvsPt"), track.pt(), track.dcaZ());

      if (track.hasTPC() && track.hasITS()) {
        registryData.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hTPCCluster"), track.tpcNClsFindable() - track.tpcNClsFindableMinusFound());
        registryData.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hShTpcClvsPt"), track.pt(), track.tpcFractionSharedCls());
        registryData.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hNclTPCFoundvsPt"), track.pt(), track.tpcNClsFound());
        registryData.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hNClTPCPidvsPt"), track.pt(), track.tpcNClsPID());
        registryData.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hNclTPCFoundvsEta"), track.eta(), track.tpcNClsFound());
        registryData.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hNClTPCPidvsEta"), track.eta(), track.tpcNClsPID());
        registryData.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hTPCnClsShared"), track.tpcNClsShared());
        registryData.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hTPCcrossedRows"), track.tpcNClsCrossedRows());
        registryData.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hTPCcrossedRowsOverFindableCls"), track.tpcCrossedRowsOverFindableCls());
        registryData.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hChi2ITSTrkSegment"), track.itsChi2NCl());
        registryData.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hITSnCls"), track.itsNCls());
      }

      registryData.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hTOFPvsBeta"), track.p(), track.beta());
      if (track.beta() > trkSelOpt.cfgTOFBetaPion && track.beta() < trkSelOpt.cfgTOFBetaPion + 0.05) { // TOF pions
        registryData.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hTOFpi"), track.eta(), track.p(), track.tpcSignal());
      }

      if (std::abs(track.eta()) < trkSelOpt.cfgTrkEtaMax) {
        if (isDCAxyWoCut(track)) {
          registryData.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hPtVsWOcutDCA"), track.pt(), track.dcaXY());
        }
      }
    }
    registryData.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hPVsPtEta"), track.p(), track.pt(), track.eta());
  }

  template <ChargeType chrg, FillType ft, bool fillHist = false, typename T, typename C>
  inline void filldEdxQA(T const& track, C const& collision, const float dEdx)
  {
    const float mult = getMult(collision);
    const float flat = fillFlat<false>(collision);
    if constexpr (fillHist) {
      if (track.p() >= trkSelOpt.cfgMomMIPMin && track.p() <= trkSelOpt.cfgMomMIPMax) {
        if (dEdx > trkSelOpt.cfgDeDxMIPMin && dEdx < trkSelOpt.cfgDeDxMIPMax) { // MIP pions
          if (cfgStoreThnSparse) {
            registryData.fill(HIST(Cprefix) + HIST(CstatCalib[ft]) + HIST(Ccharge[chrg]) + HIST("hMIP"), mult, flat, track.eta(), dEdx);
          } else {
            registryData.fill(HIST(Cprefix) + HIST(CstatCalib[ft]) + HIST(Ccharge[chrg]) + HIST("hMIP"), track.eta(), dEdx);
          }
          registryData.fill(HIST(Cprefix) + HIST(CstatCalib[ft]) + HIST(Ccharge[chrg]) + HIST("hMIPVsPhi"), track.phi(), dEdx);
          registryData.fill(HIST(Cprefix) + HIST(CstatCalib[ft]) + HIST(Ccharge[chrg]) + HIST("hMIPVsPhiVsEta"), track.phi(), dEdx, track.eta());
          registryData.fill(HIST(Cprefix) + HIST(CstatCalib[ft]) + HIST(Ccharge[chrg]) + HIST("pMIPVsPhi"), track.phi(), dEdx);
        }
        if (dEdx > trkSelOpt.cfgDeDxMIPMax + 10. && dEdx < trkSelOpt.cfgDeDxMIPMax + 30.) { // Plateau electrons
          if (std::abs(track.beta() - 1) < trkSelOpt.cfgBetaPlateuMax) {
            if (cfgStoreThnSparse) {
              registryData.fill(HIST(Cprefix) + HIST(CstatCalib[ft]) + HIST(Ccharge[chrg]) + HIST("hPlateau"), mult, flat, track.eta(), dEdx);
            } else {
              registryData.fill(HIST(Cprefix) + HIST(CstatCalib[ft]) + HIST(Ccharge[chrg]) + HIST("hPlateau"), track.eta(), dEdx);
            }
            registryData.fill(HIST(Cprefix) + HIST(CstatCalib[ft]) + HIST(Ccharge[chrg]) + HIST("hPlateauVsPhi"), track.phi(), dEdx);
            registryData.fill(HIST(Cprefix) + HIST(CstatCalib[ft]) + HIST(Ccharge[chrg]) + HIST("hPlateauVsPhiVsEta"), track.phi(), dEdx, track.eta());
            registryData.fill(HIST(Cprefix) + HIST(CstatCalib[ft]) + HIST(Ccharge[chrg]) + HIST("pPlateauVsPhi"), track.phi(), dEdx);
          }
        }
      }
    }
  }

  template <bool fillHist = false, typename C>
  bool isGoodEvent(C const& collision)
  {
    if constexpr (fillHist) {
      registryQC.fill(HIST("Events/hEvtSel"), evtSelAll);
    }
    if (evtSelOpt.cfgCustomTVX) {
      if (!collision.selection_bit(aod::evsel::kIsTriggerTVX)) {
        return false;
      }
    } else {
      if (!collision.sel8()) {
        return false;
      }
    }
    if constexpr (fillHist) {
      registryQC.fill(HIST("Events/hEvtSel"), evtSelSel8);
    }
    if (evtSelOpt.cfgRemoveITSROFrameBorder && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    if constexpr (fillHist) {
      registryQC.fill(HIST("Events/hEvtSel"), evtSelNoITSROFrameBorder);
    }
    if (evtSelOpt.cfgRemoveNoTimeFrameBorder && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    if constexpr (fillHist) {
      registryQC.fill(HIST("Events/hEvtSel"), evtSelkNoTimeFrameBorder);
    }
    if (evtSelOpt.cfgRemoveNoSameBunchPileup && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if constexpr (fillHist) {
      registryQC.fill(HIST("Events/hEvtSel"), evtSelkNoSameBunchPileup);
    }
    if (evtSelOpt.cfgRequireIsGoodZvtxFT0vsPV && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if constexpr (fillHist) {
      registryQC.fill(HIST("Events/hEvtSel"), evtSelkIsGoodZvtxFT0vsPV);
    }
    if (evtSelOpt.cfgRequireIsVertexITSTPC && !collision.selection_bit(aod::evsel::kIsVertexITSTPC)) {
      return false;
    }
    if constexpr (fillHist) {
      registryQC.fill(HIST("Events/hEvtSel"), evtSelkIsVertexITSTPC);
    }
    if (evtSelOpt.cfgRequirekIsVertexTOFmatched && !collision.selection_bit(aod::evsel::kIsVertexTOFmatched)) {
      return false;
    }
    if constexpr (fillHist) {
      registryQC.fill(HIST("Events/hEvtSel"), evtSelkIsVertexTOFmatched);
    }
    if (std::abs(collision.posZ()) > evtSelOpt.cfgCutVtxZ) {
      return false;
    }
    if constexpr (fillHist) {
      registryQC.fill(HIST("Events/hEvtSel"), evtSelVtxZ);
      registryQC.fill(HIST("Events/hVtxZ"), collision.posZ());
    }
    if (!collision.isInelGt0()) {
      return false;
    }
    if constexpr (fillHist) {
      registryQC.fill(HIST("Events/hEvtSel"), evtSelINELgt0);
    }
    if (rctCuts.requireRCTFlagChecker && !rctChecker(collision)) {
      return false;
    }
    if constexpr (fillHist) {
      registryQC.fill(HIST("Events/hEvtSel"), evtSelRCTFlagChecker);
    }
    return true;
  }

  int getFV0IndexPhi(int i_ch)
  {
    int iRing = -1;

    if (i_ch >= Cfv0IndexPhi[0] && i_ch < Cfv0IndexPhi[1]) {
      if (i_ch < Cfv0IndexPhi[0] + 4) {
        iRing = i_ch;
      } else {
        if (i_ch == Cfv0IndexPhi[1] - 1) {
          iRing = i_ch - 3; // 4;
        } else if (i_ch == Cfv0IndexPhi[1] - 2) {
          iRing = i_ch - 1; // 5;
        } else if (i_ch == Cfv0IndexPhi[1] - 3) {
          iRing = i_ch + 1; // 6;
        } else if (i_ch == Cfv0IndexPhi[1] - 4) {
          iRing = i_ch + 3; // 7;
        }
      }
    } else if (i_ch >= Cfv0IndexPhi[1] && i_ch < Cfv0IndexPhi[2]) {
      if (i_ch < Cfv0IndexPhi[2] - 4) {
        iRing = i_ch;
      } else {
        if (i_ch == Cfv0IndexPhi[2] - 1) {
          iRing = i_ch - 3; // 12;
        } else if (i_ch == Cfv0IndexPhi[2] - 2) {
          iRing = i_ch - 1; // 13;
        } else if (i_ch == Cfv0IndexPhi[2] - 3) {
          iRing = i_ch + 1; // 14;
        } else if (i_ch == Cfv0IndexPhi[2] - 4) {
          iRing = i_ch + 3; // 15;
        }
      }
    } else if (i_ch >= Cfv0IndexPhi[2] && i_ch < Cfv0IndexPhi[3]) {
      if (i_ch < Cfv0IndexPhi[3] - 4) {
        iRing = i_ch;
      } else {
        if (i_ch == Cfv0IndexPhi[3] - 1) {
          iRing = i_ch - 3; // 20;
        } else if (i_ch == Cfv0IndexPhi[3] - 2) {
          iRing = i_ch - 1; // 21;
        } else if (i_ch == Cfv0IndexPhi[3] - 3) {
          iRing = i_ch + 1; // 22;
        } else if (i_ch == Cfv0IndexPhi[3] - 4) {
          iRing = i_ch + 3; // 23;
        }
      }
    } else if (i_ch >= Cfv0IndexPhi[3] && i_ch < Cfv0IndexPhi[4]) {
      if (i_ch < Cfv0IndexPhi[3] + 4) {
        iRing = i_ch;
      } else {
        if (i_ch == Cfv0IndexPhi[4] - 5) {
          iRing = i_ch - 3; // 28;
        } else if (i_ch == Cfv0IndexPhi[4] - 6) {
          iRing = i_ch - 1; // 29;
        } else if (i_ch == Cfv0IndexPhi[4] - 7) {
          iRing = i_ch + 1; // 30;
        } else if (i_ch == Cfv0IndexPhi[4] - 8) {
          iRing = i_ch + 3; // 31;
        }
      }
    } else if (i_ch == Cfv0IndexPhi[4]) {
      iRing = Cfv0IndexPhi[4];
    } else if (i_ch == Cfv0IndexPhi[4] + 8) {
      iRing = i_ch - 7; // 33;
    } else if (i_ch == Cfv0IndexPhi[4] + 1) {
      iRing = i_ch + 1; // 34;
    } else if (i_ch == Cfv0IndexPhi[4] + 9) {
      iRing = i_ch - 6; // 35;
    } else if (i_ch == Cfv0IndexPhi[4] + 2) {
      iRing = i_ch + 2; // 36;
    } else if (i_ch == Cfv0IndexPhi[4] + 10) {
      iRing = i_ch - 5; // 37;
    } else if (i_ch == Cfv0IndexPhi[4] + 3) {
      iRing = i_ch + 3; // 38;
    } else if (i_ch == Cfv0IndexPhi[4] + 11) {
      iRing = i_ch - 4; // 39;
    } else if (i_ch == Cfv0IndexPhi[4] + 15) {
      iRing = i_ch - 7; // 40;
    } else if (i_ch == Cfv0IndexPhi[4] + 7) {
      iRing = i_ch + 2; // 41;
    } else if (i_ch == Cfv0IndexPhi[4] + 14) {
      iRing = i_ch - 4; // 42;
    } else if (i_ch == Cfv0IndexPhi[4] + 6) {
      iRing = i_ch + 5; // 43;
    } else if (i_ch == Cfv0IndexPhi[4] + 13) {
      iRing = i_ch - 1; // 44;
    } else if (i_ch == Cfv0IndexPhi[4] + 5) {
      iRing = i_ch + 8; // 45;
    } else if (i_ch == Cfv0IndexPhi[4] + 12) {
      iRing = i_ch + 2; // 46;
    } else if (i_ch == Cfv0IndexPhi[4] + 4) {
      iRing = i_ch + 11; // 47;
    }
    return iRing;
  }

  template <typename C, bool isMC = false>
  float getMult(C const& collision)
  {
    float val = -999.0;
    switch (multEst) {
      case MultE::CnoMult:
        return val;
        break;
      case MultE::CmultFT0M:
        if constexpr (!isMC) {
          return collision.centFT0M();
        } else {
          return collision.multMCFT0C() + collision.multMCFT0A();
        }
        break;
      case MultE::CmultTPC:
        if constexpr (!isMC) {
          return collision.multTPC();
        } else {
          LOG(fatal) << "No valid multiplicity estimator: " << multEst;
          return val;
        }
        break;
      default:
        return collision.centFT0M();
        break;
    }
  }

  float getMultMC(CollsMCExtraMult::iterator const& collision)
  {
    return getMult<CollsMCExtraMult::iterator, true>(collision);
  }

  template <typename T, std::size_t S>
  float calcFlatenicity(std::array<T, S> const& signals)
  {
    static_assert(S != 0);

    int entries = signals.size();
    float flat{-1};
    float mRho{0};
    for (int iCell = 0; iCell < entries; ++iCell) {
      if (signals[iCell] > Cnull) {
        mRho += 1.0 * signals[iCell];
      }
    }
    // average activity per cell
    mRho /= (1.0 * entries);
    if (mRho <= 0) {
      return -1;
    }
    // get sigma
    float sRhoTmp{0};
    float sRho{0};
    for (int iCell = 0; iCell < entries; ++iCell) {
      if (signals[iCell] > Cnull) {
        sRhoTmp += std::pow(1.0 * signals[iCell] - mRho, 2);
      }
    }
    sRhoTmp /= (1.0 * entries * entries);
    sRho = std::sqrt(sRhoTmp);
    if (mRho > Cnull) {
      flat = sRho / mRho;
    } else {
      flat = -1;
    }
    return flat;
  }

  template <bool fillHist = true, typename C>
  float fillFlat(C const& collision)
  {
    rhoLatticeFV0.fill(0);
    fv0AmplitudeWoCalib.fill(0);
    if (collision.has_foundFV0()) {
      auto fv0 = collision.foundFV0();
      std::bitset<8> fV0Triggers = fv0.triggerMask();
      bool isOkFV0OrA = fV0Triggers[o2::fit::Triggers::bitA];
      if (isOkFV0OrA) {
        for (std::size_t ich = 0; ich < fv0.channel().size(); ich++) {
          float amplCh = fv0.amplitude()[ich];
          int chv0 = fv0.channel()[ich];
          int chv0phi = getFV0IndexPhi(chv0);
          if constexpr (fillHist) {
            registryQC.fill(HIST("FV0/hFV0amp"), chv0, amplCh);
            registryQC.fill(HIST("FV0/pFV0amp"), chv0, amplCh);
            if (applyCalibGain) {
              registryQC.fill(HIST("FV0/hFV0ampCorr"), chv0, amplCh / fv0AmplCorr[chv0]);
            }
          }
          if (amplCh > Cnull) {
            if (applyCalibGain) { // equalize gain channel-by-channel
              amplCh /= fv0AmplCorr[chv0];
            }
            if (chv0phi > Cnull) {
              fv0AmplitudeWoCalib[chv0phi] = amplCh;
              if constexpr (fillHist) {
                registryQC.fill(HIST("FV0/hFV0AmplWCalib"), ich, fv0AmplitudeWoCalib[ich]);
              }
              if (chv0 < CinnerFV0) {
                rhoLatticeFV0[chv0phi] += amplCh;
              } else { // two channels per bin
                rhoLatticeFV0[chv0phi] += amplCh / 2.;
              }
              if constexpr (fillHist) {
                registryQC.fill(HIST("FV0/hFV0AmplvsVtxzWoCalib"), collision.posZ(), rhoLatticeFV0[chv0phi]);
              }
              if (applyCalibVtx) {
                rhoLatticeFV0[chv0phi] *= zVtxMap->GetBinContent(zVtxMap->GetXaxis()->FindBin(chv0phi), zVtxMap->GetYaxis()->FindBin(collision.posZ()));
                if constexpr (fillHist) {
                  registryQC.fill(HIST("FV0/hFV0AmplvsVtxzCalib"), collision.posZ(), rhoLatticeFV0[chv0phi]);
                }
              }
            }
          }
        }
        float flattenicityFV0 = calcFlatenicity(rhoLatticeFV0);
        return 1. - flattenicityFV0;
      } else {
        return 9999;
      }
    } else {
      return 9999;
    }
  }

  template <typename T>
  bool isDCAxyCut(T const& track) const
  {
    if (isCustomTracks.value) {
      for (int i = 0; i < static_cast<int>(TrackSelection::TrackCuts::kNCuts); i++) {
        if (i == static_cast<int>(TrackSelection::TrackCuts::kDCAxy)) {
          continue;
        }
        if (!selTrkGlobal.IsSelected(track, static_cast<TrackSelection::TrackCuts>(i))) {
          return false;
        }
      }
      return (std::abs(track.dcaXY()) <= (maxDcaXYFactor.value * (0.0105f + 0.0350f / std::pow(track.pt(), 1.1f))));
    }
    return track.isGlobalTrack();
  }

  template <typename T>
  bool isDCAxyWoCut(T const& track) const
  {
    if (isCustomTracks.value) {
      for (int i = 0; i < static_cast<int>(TrackSelection::TrackCuts::kNCuts); i++) {
        if (i == static_cast<int>(TrackSelection::TrackCuts::kDCAxy)) {
          continue;
        }
        if (i == static_cast<int>(TrackSelection::TrackCuts::kDCAz)) {
          continue;
        }
        if (!selTrkGlobal.IsSelected(track, static_cast<TrackSelection::TrackCuts>(i))) {
          return false;
        }
      }
      return true;
    }
    return track.isGlobalTrackWoDCA();
  }

  Preslice<MyPIDTracks> perCol = aod::track::collisionId;
  Preslice<aod::V0Datas> perColV0s = aod::v0::collisionId;

  template <typename C>
  void processData(C const& collisions,
                   MyPIDTracks const& tracks,
                   aod::V0Datas const& v0s,
                   aod::BCsWithTimestamps const& bcs)
  {
    for (const auto& collision : collisions) {
      if (!isGoodEvent<true>(collision)) {
        continue;
      }
      auto tracksPerCollision = tracks.sliceBy(perCol, collision.globalIndex());
      auto v0sPerCollision = v0s.sliceBy(perColV0s, collision.globalIndex());
      v0sPerCollision.bindExternalIndices(&tracks);
      filldEdx(tracksPerCollision, v0sPerCollision, collision, bcs);
      if (cfgFillDCAxyHist) {
        static_for<0, 4>([&](auto i) {
          fillDCA<i>(tracksPerCollision, collision, bcs);
        });
      }
    }
  }

  void processFlat(
    Colls const& collisions,
    MyPIDTracks const& tracks,
    aod::V0Datas const& v0s,
    aod::BCsWithTimestamps const& bcs,
    aod::FT0s const&, aod::FV0As const&)
  {
    processData<Colls>(collisions, tracks, v0s, bcs);
  }
  PROCESS_SWITCH(FlattenictyPikp, processFlat, "process Flat data inclusive", true);

  template <bool fillHist = true, typename McPart>
  float fillFlatMC(McPart const& mcparts)
  {
    int nFV0sectors = 0;
    float minPhi = 0;
    float maxPhi = 0;
    float dPhi = 0;

    float etaMinFV0bins[CmaxRingsFV0] = {0.0};
    float etaMaxFV0bins[CmaxRingsFV0] = {0.0};
    for (int i = 0; i < CmaxRingsFV0; ++i) {
      etaMaxFV0bins[i] = CmaxEtaFV0 - i * CdEtaFV0;
      if (i < CmaxRingsFV0 - 1) {
        etaMinFV0bins[i] = CmaxEtaFV0 - (i + 1) * CdEtaFV0;
      } else {
        etaMinFV0bins[i] = CminEtaFV0;
      }
    }

    rhoLatticeFV0.fill(0);
    std::vector<float> vNch;
    float nChFV0{0};
    for (const auto& mcPart : mcparts) {
      if (!isChrgParticle(mcPart.pdgCode())) {
        continue;
      }
      if (!mcPart.isPhysicalPrimary()) {
        continue;
      }
      if (mcPart.pt() <= 0.) {
        continue;
      }

      auto etaMc = mcPart.eta();
      auto phiMc = mcPart.phi();

      int isegment = 0;
      for (int ieta = 0; ieta < CmaxRingsFV0; ieta++) {

        nFV0sectors = CnCellsFV0 / 6.;
        if (ieta == CmaxRingsFV0 - 1) {
          nFV0sectors = CnCellsFV0 / 3.;
        }

        for (int iphi = 0; iphi < nFV0sectors; iphi++) {

          minPhi = iphi * TwoPI / nFV0sectors;
          maxPhi = (iphi + 1) * TwoPI / nFV0sectors;
          dPhi = std::abs(maxPhi - minPhi);

          if (etaMc >= etaMinFV0bins[ieta] && etaMc < etaMaxFV0bins[ieta] && phiMc >= minPhi && phiMc < maxPhi) {
            rhoLatticeFV0[isegment] += 1. / std::abs(CdEtaFV0 * dPhi);
            nChFV0++;
          }
          isegment++;
        }
      }
    }

    vNch.push_back(nChFV0);
    auto flatFV0 = 1. - calcFlatenicity(rhoLatticeFV0);

    if constexpr (fillHist) {
      registryMC.fill(HIST("Events/ResponseGen"), vNch[0], flatFV0);
      registryMC.fill(HIST("Events/h1flatencityFV0MCGen"), flatFV0);
    }
    vNch.clear();
    return flatFV0;
  }

  template <int pidSgn, o2::track::PID::ID id>
  void bookMcHist()
  {
    AxisSpec ptAxis{binOpt.axisPt, "#it{p}_{T} (GeV/#it{c})"};
    constexpr int ChistIdx = id + pidSgn * Npart;
    auto idx = static_cast<int>(id);
    const std::string strID = Form("/%s/%s", (pidSgn == CnullInt && id < Npart) ? "pos" : "neg", Pid[idx]);
    hPtEffRec[ChistIdx] = registryMC.add<TH1>("Tracks/hPtEffRec" + strID, " ; #it{p}_{T} (GeV/#it{c})", kTH1F, {ptAxis});
    hPtEffGen[ChistIdx] = registryMC.add<TH1>("Tracks/hPtEffGen" + strID, " ; #it{p}_{T} (GeV/#it{c})", kTH1F, {ptAxis});
  }

  template <int pidSgn, o2::track::PID::ID id>
  void initEfficiency()
  {
    static_assert(pidSgn == CnullInt || pidSgn == ConeInt);
    static_assert(id > CnullInt || id < Npart);
    constexpr int Cidx = id + pidSgn * Npart;
    const TString partName = PidChrg[Cidx];
    THashList* lhash = new THashList();
    lhash->SetName(partName);
    listEfficiency->Add(lhash);

    auto bookEff = [&](const TString eName, auto h) {
      const TAxis* axis = h->GetXaxis();
      TString eTitle = h->GetTitle();
      eTitle.ReplaceAll("Numerator", "").Strip(TString::kBoth);
      eTitle = Form("%s;%s;Efficiency", eTitle.Data(), axis->GetTitle());
      lhash->Add(new TEfficiency(eName, eTitle, axis->GetNbins(), axis->GetXbins()->GetArray()));
    };

    const int idx = id + pidSgn * Npart;
    bookEff("hEffvsPt", hPtEffRec[idx]);
  }

  template <int pidSgn, o2::track::PID::ID id>
  void fillEfficiency()
  {
    static_assert(pidSgn == CnullInt || pidSgn == ConeInt);
    constexpr int ChistIdx = id + pidSgn * Npart;
    const char* partName = PidChrg[ChistIdx];
    THashList* lhash = static_cast<THashList*>(listEfficiency->FindObject(partName));
    if (!lhash) {
      LOG(warning) << "No efficiency object found for particle " << partName;
      return;
    }

    auto fillEff = [&](const TString eName, auto num, auto den) {
      TEfficiency* eff = static_cast<TEfficiency*>(lhash->FindObject(eName));
      if (!eff) {
        LOG(warning) << "Cannot find TEfficiency " << eName;
        return;
      }
      eff->SetTotalHistogram(*den, "f");
      eff->SetPassedHistogram(*num, "f");
    };
    fillEff("hEffvsPt", hPtEffRec[ChistIdx], hPtEffGen[ChistIdx]);
  }

  template <int pidSgn, o2::track::PID::ID id>
  void fillMCRecTrack(MyLabeledPIDTracks::iterator const& track, const float mult, const float flat)
  {
    static_assert(pidSgn == CnullInt || pidSgn == ConeInt);
    constexpr int ChistIdx = id + pidSgn * Npart;
    // LOG(debug) << "fillMCRecTrack for pidSgn '" << pidSgn << "' and id '" << static_cast<int>(id) << " with index " << ChistIdx;
    const aod::McParticles::iterator& mcParticle = track.mcParticle();
    const CollsGen::iterator& collision = track.collision_as<soa::SmallGroups<CollsGen>>();
    // const CollsGen::iterator& collision = track.collision_as<CollsGen>();

    if (!isChrgParticle(mcParticle.pdgCode())) {
      return;
    }
    if (std::abs(mcParticle.eta()) > trkSelOpt.cfgTrkEtaMax) {
      return;
    }
    if (mcParticle.pt() < trkSelOpt.cfgTrkPtMin) {
      return;
    }
    if (!isPID<pidSgn, id>(mcParticle)) {
      return;
    }

    if ((collision.has_mcCollision() && (mcParticle.mcCollisionId() != collision.mcCollisionId())) || !collision.has_mcCollision()) {
      if (!mcParticle.isPhysicalPrimary()) {
        if (mcParticle.getProcess() == CprocessIdWeak) {
          hDCAxyRecBadCollWeak[ChistIdx]->Fill(track.pt(), track.dcaXY());
        } else {
          hDCAxyRecBadCollMat[ChistIdx]->Fill(track.pt(), track.dcaXY());
        }
      } else {
        hDCAxyRecBadCollPrim[ChistIdx]->Fill(track.pt(), track.dcaXY());
      }
    }

    if (collision.has_mcCollision() && (mcParticle.mcCollisionId() == collision.mcCollisionId())) {
      if (!mcParticle.isPhysicalPrimary()) {
        if (mcParticle.getProcess() == CprocessIdWeak) {
          hPtEffRecGoodCollWeak[ChistIdx]->Fill(mult, flat, track.pt());
          hPtVsDCAxyRecGoodCollWeak[ChistIdx]->Fill(track.pt(), track.dcaXY());
        } else {
          hPtEffRecGoodCollMat[ChistIdx]->Fill(mult, flat, track.pt());
          hPtVsDCAxyRecGoodCollMat[ChistIdx]->Fill(track.pt(), track.dcaXY());
        }
      } else {
        hPtEffRecGoodCollPrim[ChistIdx]->Fill(mult, flat, track.pt());
        hPtVsDCAxyRecGoodCollPrim[ChistIdx]->Fill(track.pt(), track.dcaXY());
        if (isDCAxyCut(track)) {
          hPtEffRec[ChistIdx]->Fill(mcParticle.pt());
        }
      }
    }
  }

  template <int pidSgn, o2::track::PID::ID id, bool isGtZeroColl = false>
  void fillMCGenRecEvt(aod::McParticles::iterator const& mcParticle, const float mult, const float flat)
  {
    static_assert(pidSgn == CnullInt || pidSgn == ConeInt);
    constexpr int ChistIdx = id + pidSgn * Npart;

    if (!isPID<pidSgn, id>(mcParticle)) {
      return;
    }
    if constexpr (isGtZeroColl) {
      hPtGenRecEvtGtZero[ChistIdx]->Fill(mult, flat, mcParticle.pt());
      if (mcParticle.isPhysicalPrimary()) {
        hPtGenPrimRecEvtGtZero[ChistIdx]->Fill(mult, flat, mcParticle.pt());
        hPtEffGen[ChistIdx]->Fill(mcParticle.pt());
      }
    } else {
      hPtGenRecEvt[ChistIdx]->Fill(mult, flat, mcParticle.pt());
      if (mcParticle.isPhysicalPrimary()) {
        hPtGenPrimRecEvt[ChistIdx]->Fill(mult, flat, mcParticle.pt());
      }
    }
  }

  template <int pidSgn, o2::track::PID::ID id, bool evtSel = false>
  void fillMCGen(aod::McParticles::iterator const& mcParticle, const float mult, const float flat)
  {
    static_assert(pidSgn == CnullInt || pidSgn == ConeInt);
    constexpr int ChistIdx = id + pidSgn * Npart;

    if (!isPID<pidSgn, id>(mcParticle)) {
      return;
    }

    if constexpr (evtSel) {
      if (!mcParticle.isPhysicalPrimary()) {
        if (mcParticle.getProcess() == CprocessIdWeak) {
          hPtEffGenWeakEvtSelGen[ChistIdx]->Fill(mult, flat, mcParticle.pt());
        } else {
          hPtEffGenMatEvtSelGen[ChistIdx]->Fill(mult, flat, mcParticle.pt());
        }
      } else {
        hPtEffGenPrimEvtSelGen[ChistIdx]->Fill(mult, flat, mcParticle.pt());
      }
    } else {
      if (!mcParticle.isPhysicalPrimary()) {
        if (mcParticle.getProcess() == CprocessIdWeak) {
          hPtEffGenWeak[ChistIdx]->Fill(mult, flat, mcParticle.pt());
        } else {
          hPtEffGenMat[ChistIdx]->Fill(mult, flat, mcParticle.pt());
        }
      } else {
        hPtEffGenPrim[ChistIdx]->Fill(mult, flat, mcParticle.pt());
        // hPtEffGen[ChistIdx]->Fill(mcParticle.pt());
      }
    }
  }

  Preslice<MyLabeledPIDTracks> perCollTrk = aod::track::collisionId;

  void processMC(CollsMCExtraMult::iterator const& mcCollision,
                 soa::SmallGroups<CollsGen> const& collisions,
                 aod::BCsWithTimestamps const& /*bcs*/,
                 aod::FV0As const& /*fv0s*/,
                 aod::McParticles const& particles,
                 MyLabeledPIDTracks const& tracks)
  {
    LOGP(debug, "MC col {} has {} reco cols", mcCollision.globalIndex(), collisions.size());
    auto multMC = -1.;
    if (evtSelOpt.useMultMCmidrap || multEst == 2) { // use generated Nch in ∣eta∣ < 0.8
      multMC = countPart(particles);
    } else {
      multMC = getMultMC(mcCollision); // using McCentFT0Ms
    }
    const float flatMC = fillFlatMC<true>(particles);
    registryMC.fill(HIST("Events/hFlatMCGen"), flatMC);

    // Loop on rec collisions
    // Obtain here: Denominator of tracking efficiency; Numerator event and signal loss
    //
    bool gtZeroColl = false;
    auto multRecGt1 = -999;
    auto flatRec = -999;
    for (const auto& collision : collisions) {
      if (!isGoodEvent<false>(collision)) {
        continue;
      }
      registryMC.fill(HIST("Events/hCentVsFlatRecINELgt0"), multRecGt1, flatRec); // Evt split den
      if (evtSelOpt.cfgRemoveSplitVertex && collision.globalIndex() != mcCollision.bestCollisionIndex()) {
        continue;
      }
      registryMC.fill(HIST("Events/hCentVsFlatRecINELgt0wRecEvt"), multRecGt1, flatRec); // Evt split num,  w/ Nrec > 0
      gtZeroColl = true;
      multRecGt1 = getMult(collision);
      flatRec = fillFlat<true>(collision);
    }

    if (gtZeroColl) {
      registryMC.fill(HIST("Events/hCentVsFlatRecINELgt0wRecEvtSel"), multRecGt1, flatRec); // Evt split num,  w/ Nrec > 0 + Evt. sel
      registryMC.fill(HIST("Events/hNchGenVsCent"), multMC, multRecGt1);
      registryMC.fill(HIST("Events/hNchVsFlatGenINELgt0wRecEvtSel"), multMC, flatMC); // Evt loss num,   w/ Nrec > 0 + Evt. sel
    }

    for (const auto& particle : particles) {
      if (!isChrgParticle(particle.pdgCode())) {
        continue;
      }
      if (!particle.isPhysicalPrimary()) {
        continue;
      }
      if (std::abs(particle.eta()) > trkSelOpt.cfgTrkEtaMax) {
        continue;
      }
      if (particle.pt() < trkSelOpt.cfgTrkPtMin) {
        continue;
      }
      if (gtZeroColl) {
        static_for<0, 1>([&](auto pidSgn) {
          fillMCGenRecEvt<pidSgn, o2::track::PID::Pion, true>(particle, multMC, flatMC);
          fillMCGenRecEvt<pidSgn, o2::track::PID::Kaon, true>(particle, multMC, flatMC);
          fillMCGenRecEvt<pidSgn, o2::track::PID::Proton, true>(particle, multMC, flatMC);
        });
        static_for<0, 4>([&](auto i) {
          constexpr int Cidx = i.value;
          if (std::fabs(particle.pdgCode()) == PDGs[Cidx]) {
            registryMC.fill(HIST(Cprefix) + HIST(CspeciesAll[Cidx]) + HIST(CpTrecCollPrimSgn), multMC, flatMC, particle.pt());        // Sgn loss num
            registryMC.fill(HIST(Cprefix) + HIST(CspeciesAll[Cidx]) + HIST(CpTeffGenPrimRecEvt), multRecGt1, flatRec, particle.pt()); // Tracking eff. den
          }
        });
      } else {
        static_for<0, 1>([&](auto pidSgn) {
          fillMCGenRecEvt<pidSgn, o2::track::PID::Pion>(particle, multMC, flatMC);
          fillMCGenRecEvt<pidSgn, o2::track::PID::Kaon>(particle, multMC, flatMC);
          fillMCGenRecEvt<pidSgn, o2::track::PID::Proton>(particle, multMC, flatMC);
        });
      }
    }

    // Loop on rec collisions
    // Obtain here: Numerator of tracking efficiency; Secondary contamination correction
    for (const auto& collision : collisions) {
      if (trkSelOpt.cfgRejectTrkAtTPCSector || applyCalibGain || applyCalibVtx) {
        auto bc = collision.bc_as<aod::BCsWithTimestamps>();
        int currentRun = bc.runNumber();
        if (runNumber != currentRun) {
          initCCDB(bc);
          runNumber = currentRun;
        }
      }
      if (!isGoodEvent<false>(collision)) {
        continue;
      }
      registryMC.fill(HIST("Events/hVtxZRec"), collision.posZ());
      if (evtSelOpt.cfgRemoveSplitVertex && collision.globalIndex() != mcCollision.bestCollisionIndex()) {
        continue;
      }
      // Rec tracks; track selection w/ DCA open (for secondaries), w/ DCA close (for efficiency)
      // Obtain here: DCAxy for sec contamination, MC closure
      const auto& groupedTrks = tracks.sliceBy(perCollTrk, collision.globalIndex());
      int nTrk = 0;
      for (const auto& track : groupedTrks) {
        if (!track.has_collision()) {
          continue;
        }
        if (std::abs(track.eta()) > trkSelOpt.cfgTrkEtaMax) {
          continue;
        }
        if (track.pt() < trkSelOpt.cfgTrkPtMin) {
          continue;
        }
        if (trkSelOpt.cfgApplyNcl && track.tpcNClsFound() < trkSelOpt.cfgNclTPCMin) {
          continue;
        }
        if (trkSelOpt.cfgApplyNclPID && track.tpcNClsPID() < trkSelOpt.cfgNclPidTPCMin) {
          continue;
        }
        float phiModn = track.phi();
        phiMod(phiModn, magField, track.sign());
        if (trkSelOpt.cfgRejectTrkAtTPCSector && (track.pt() >= trkSelOpt.cfgPhiCutPtMin && phiModn < fPhiCutHigh->Eval(track.pt()) && phiModn > fPhiCutLow->Eval(track.pt()))) {
          continue;
        }
        if (!isDCAxyWoCut(track)) {
          continue;
        }
        if (!track.has_mcParticle()) {
          registryMC.fill(HIST("Tracks/hPtFakes"), multMC, flatMC, track.pt());
          continue;
        }
        auto particle = track.mcParticle_as<aod::McParticles>();
        if (collision.mcCollisionId() != particle.mcCollisionId()) {
          continue;
        }
        if (!isChrgParticle(particle.pdgCode())) {
          continue;
        }
        if (std::abs(particle.eta()) > trkSelOpt.cfgTrkEtaMax) {
          continue;
        }
        if (particle.pt() < trkSelOpt.cfgTrkPtMin) {
          continue;
        }
        static_for<0, 1>([&](auto pidSgn) { // for checking purposes only: use gen Nch, gen Flat
          fillMCRecTrack<pidSgn, o2::track::PID::Pion>(track, multMC, flatMC);
          fillMCRecTrack<pidSgn, o2::track::PID::Kaon>(track, multMC, flatMC);
          fillMCRecTrack<pidSgn, o2::track::PID::Proton>(track, multMC, flatMC);
        });
        static_for<0, 4>([&](auto i) {
          constexpr int Cidx = i.value;
          if (std::sqrt(std::pow(std::fabs(o2::aod::pidutils::tpcNSigma<Cidx>(track)), 2) + std::pow(std::fabs(o2::aod::pidutils::tofNSigma<Cidx>(track)), 2) < trkSelOpt.cfgDcaNsigmaCombinedMax)) {
            if (std::fabs(particle.pdgCode()) == PDGs[Cidx]) {
              if (!particle.isPhysicalPrimary()) {
                if (particle.getProcess() == CprocessIdWeak) {
                  registryMC.fill(HIST(Cprefix) + HIST(CspeciesAll[Cidx]) + HIST(CpTvsDCAxyWeakAll), multRecGt1, flatRec, track.pt(), track.dcaXY());
                } else {
                  registryMC.fill(HIST(Cprefix) + HIST(CspeciesAll[Cidx]) + HIST(CpTvsDCAxyMatAll), multRecGt1, flatRec, track.pt(), track.dcaXY());
                }
              } else {
                registryMC.fill(HIST(Cprefix) + HIST(CspeciesAll[Cidx]) + HIST(CpTvsDCAxyPrimAll), multRecGt1, flatRec, track.pt(), track.dcaXY());
                registryMC.fill(HIST(Cprefix) + HIST(CspeciesAll[Cidx]) + HIST(CdEdxMcRecPrim), track.eta(), multRecGt1, flatRec, track.p(), track.tpcSignal());
              }
              registryMC.fill(HIST(Cprefix) + HIST(CspeciesAll[Cidx]) + HIST(CpTvsDCAxyAll), multRecGt1, flatRec, track.pt(), track.dcaXY());
            }
          }
        });
        if (isDCAxyCut(track)) {
          static_for<0, 4>([&](auto i) {
            constexpr int Cidx = i.value;
            if (std::fabs(particle.pdgCode()) == PDGs[Cidx]) {
              if (particle.isPhysicalPrimary()) {
                registryMC.fill(HIST(Cprefix) + HIST(CspeciesAll[Cidx]) + HIST(CdEdxMcRecPrimSel), track.eta(), multRecGt1, flatRec, track.p(), track.tpcSignal());
                registryMC.fill(HIST(Cprefix) + HIST(CspeciesAll[Cidx]) + HIST(CpTeffPrimRecEvt), multRecGt1, flatRec, track.pt()); // Tracking eff. num
                registryMC.fill(HIST(Cprefix) + HIST(CspeciesAll[Cidx]) + HIST(CpTmcClosureRec), multMC, flatMC, track.pt());       // closure
              }
            }
          });
          nTrk++;
        }
      }
      registryQC.fill(HIST("Events/hNchVsCent"), nTrk, multRecGt1);
    }
    static_for<0, 1>([&](auto pidSgn) {
      fillEfficiency<pidSgn, o2::track::PID::Pion>();
      fillEfficiency<pidSgn, o2::track::PID::Kaon>();
      fillEfficiency<pidSgn, o2::track::PID::Proton>();
    });

    // Loop on generated particles (no requirement on availaability of reconstructed collision; no event selection)
    //
    for (const auto& particle : particles) {
      if (!isChrgParticle(particle.pdgCode())) {
        continue;
      }
      if (std::abs(particle.eta()) > trkSelOpt.cfgTrkEtaMax) {
        continue;
      }
      if (particle.pt() < trkSelOpt.cfgTrkPtMin) {
        continue;
      }
      static_for<0, 1>([&](auto pidSgn) {
        fillMCGen<pidSgn, o2::track::PID::Pion>(particle, multMC, flatMC);
        fillMCGen<pidSgn, o2::track::PID::Kaon>(particle, multMC, flatMC);
        fillMCGen<pidSgn, o2::track::PID::Proton>(particle, multMC, flatMC);
      });
    }

    // Obtain here: Denominator of signal loss and event loss; MC closure
    //
    registryMC.fill(HIST("Events/hEvtMcGen"), 0.5);
    if (evtSelOpt.useZVtxCutMC && std::abs(mcCollision.posZ()) > evtSelOpt.cfgCutVtxZ) {
      return;
    }
    registryMC.fill(HIST("Events/hVtxZGen"), mcCollision.posZ());
    registryMC.fill(HIST("Events/hEvtMcGen"), 1.5);
    if (evtSelOpt.useINELCutMC) {
      if (!o2::pwglf::isINELgt0mc(particles, pdg)) {
        return;
      }
    }
    registryMC.fill(HIST("Events/hEvtMcGen"), 2.5);
    if (evtSelOpt.cfgUseInelgt0wTVX && !isInelGt0wTVX(particles)) { // TVX trigger: FT0A + FT0C acceptance
      return;
    }
    registryMC.fill(HIST("Events/hEvtMcGen"), 3.5);
    registryMC.fill(HIST("Events/hNchVsFlatGenINELgt0"), multMC, flatMC); // Evt loss den

    for (const auto& particle : particles) {
      if (!isChrgParticle(particle.pdgCode())) {
        continue;
      }
      if (std::abs(particle.eta()) > trkSelOpt.cfgTrkEtaMax) {
        continue;
      }
      if (particle.pt() < trkSelOpt.cfgTrkPtMin) {
        continue;
      }
      static_for<0, 1>([&](auto pidSgn) {
        fillMCGen<pidSgn, o2::track::PID::Pion, true>(particle, multMC, flatMC);
        fillMCGen<pidSgn, o2::track::PID::Kaon, true>(particle, multMC, flatMC);
        fillMCGen<pidSgn, o2::track::PID::Proton, true>(particle, multMC, flatMC);
      });
      static_for<0, 4>([&](auto i) {
        constexpr int Cidx = i.value;
        // LOG(debug) << "fillMCGen for pidSgn '" << pidSgn << "' and id '" << static_cast<int>(id) << " with index " << ChistIdx;
        if (particle.isPhysicalPrimary()) {
          if (std::fabs(particle.pdgCode()) == PDGs[Cidx]) {
            registryMC.fill(HIST(Cprefix) + HIST(CspeciesAll[Cidx]) + HIST(CpTgenPrimSgn), multMC, flatMC, particle.pt());       // Sgn loss den
            registryMC.fill(HIST(Cprefix) + HIST(CspeciesAll[Cidx]) + HIST(CpTmcClosureGenPrim), multMC, flatMC, particle.pt()); // closure
          }
        }
      });
    }
  }
  PROCESS_SWITCH(FlattenictyPikp, processMC, "process MC", false);

  template <typename ObjType>
  ObjType* getForTsOrRun(std::string const& fullPath, int64_t timestamp, int runNumber)
  {
    if (cfgUseCcdbForRun) {
      return ccdb->getForRun<ObjType>(fullPath, runNumber);
    } else {
      return ccdb->getForTimeStamp<ObjType>(fullPath, timestamp);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<FlattenictyPikp>(cfgc)};
}
