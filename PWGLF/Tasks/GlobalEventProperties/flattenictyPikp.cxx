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

#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/mcCentrality.h"
#include "PWGLF/Utils/inelGt.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/McCollisionExtra.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/MathConstants.h"
#include "DataFormatsFIT/Triggers.h"
#include "DetectorsCommonDataFormats/AlignParam.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/StaticFor.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/Track.h"
#include <CCDB/BasicCCDBManager.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DataFormatsParameters/GRPObject.h>

#include "TEfficiency.h"
#include "THashList.h"
#include "TPDGCode.h"
#include <TF1.h>
#include <TGraph.h>
#include <TH1.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <map>
#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

using std::string;
using namespace o2;
using namespace o2::framework;
using namespace o2::constants::math;
using namespace o2::aod::rctsel;

auto static constexpr CminCharge = 3.f;
static constexpr float Cnull = 0.0f;
static constexpr float Cone = 1.0f;

// FV0 specific constants
static constexpr int CmaxRingsFV0 = 5;
static constexpr int CnCellsFV0 = 48;
static constexpr int CinnerFV0 = 32;
static constexpr float Cfv0IndexPhi[5] = {0., 8., 16., 24., 32.};
static constexpr float CmaxEtaFV0 = 5.1;
static constexpr float CminEtaFV0 = 2.2;
static constexpr float CdEtaFV0 = (CmaxEtaFV0 - CminEtaFV0) / CmaxRingsFV0;
// PID names
static constexpr int CprocessIdWeak = 4;
static constexpr int Ncharges = 2;
static constexpr o2::track::PID::ID Npart = 5;
static constexpr o2::track::PID::ID NpartChrg = Npart * Ncharges;
static constexpr int PDGs[] = {11, 13, 211, 321, 2212};
static constexpr int PidSgn[NpartChrg] = {11, 13, 211, 321, 2212, -11, -13, -211, -321, -2212};
static constexpr const char* Pid[Npart] = {"el", "mu", "pi", "ka", "pr"};
static constexpr const char* PidChrg[NpartChrg] = {"e^{-}", "#mu^{-}", "#pi^{+}", "K^{+}", "p", "e^{+}", "#mu^{+}", "#pi^{-}", "K^{-}", "#bar{p}"};
static constexpr std::string_view Cspecies[NpartChrg] = {"Elminus", "Muplus", "PiPlus", "KaPlus", "Pr", "ElPlus", "MuMinus", "PiMinus", "KaMinus", "PrBar"};
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
static constexpr std::string_view CpTgenPrimSgn = "/hPtGenPrimSgn";
static constexpr std::string_view CpTgenPrimSgnF = "Tracks/{}/hPtGenPrimSgn";
static constexpr std::string_view CpTgenPrimSgnINEL = "/hPtGenPrimSgnINEL";
static constexpr std::string_view CpTgenPrimSgnINELF = "Tracks/{}/hPtGenPrimSgnINEL";
static constexpr std::string_view CpTrecCollPrimSgn = "/hPtRecCollPrimSgn";
static constexpr std::string_view CpTrecCollPrimSgnF = "Tracks/{}/hPtRecCollPrimSgn";
static constexpr std::string_view CpTrecCollPrimSgnINEL = "/hPtRecCollPrimSgnINEL";
static constexpr std::string_view CpTrecCollPrimSgnINELF = "Tracks/{}/hPtRecCollPrimSgnINEL";
static constexpr std::string_view CpTGenRecCollPrimSgn = "/hPtGenRecCollPrimSgn";
static constexpr std::string_view CpTGenRecCollPrimSgnF = "Tracks/{}/hPtGenRecCollPrimSgn";
static constexpr std::string_view CpTGenRecCollPrimSgnINEL = "/hPtGenRecCollPrimSgnINEL";
static constexpr std::string_view CpTGenRecCollPrimSgnINELF = "Tracks/{}/hPtGenRecCollPrimSgnINEL";
static constexpr std::string_view CpTmcClosurePrim = "/hPtMCclosurePrim";
static constexpr std::string_view CpTmcClosurePrimF = "Tracks/{}/hPtMCclosurePrim";

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

std::array<std::shared_ptr<TH3>, NpartChrg> hPtEffRecPrim{};
std::array<std::shared_ptr<TH3>, NpartChrg> hPtEffRecWeak{};
std::array<std::shared_ptr<TH3>, NpartChrg> hPtEffRecMat{};
std::array<std::shared_ptr<TH1>, NpartChrg> hPtEffRec{};
std::array<std::shared_ptr<TH1>, NpartChrg> hPtEffGen{};
std::array<std::shared_ptr<TH1>, NpartChrg> hPtGenRecEvt{};
std::array<std::shared_ptr<TH1>, NpartChrg> hPtGenPrimRecEvt{};
std::array<std::shared_ptr<TH3>, NpartChrg> hPtEffGenPrim{};
std::array<std::shared_ptr<TH3>, NpartChrg> hPtEffGenWeak{};
std::array<std::shared_ptr<TH3>, NpartChrg> hPtEffGenMat{};
std::array<std::shared_ptr<TH2>, NpartChrg> hPtVsDCAxyPrim{};
std::array<std::shared_ptr<TH2>, NpartChrg> hPtVsDCAxyWeak{};
std::array<std::shared_ptr<TH2>, NpartChrg> hPtVsDCAxyMat{};
std::array<std::shared_ptr<TH2>, NpartChrg> hDCAxyBadCollPrim{};
std::array<std::shared_ptr<TH2>, NpartChrg> hDCAxyBadCollWeak{};
std::array<std::shared_ptr<TH2>, NpartChrg> hDCAxyBadCollMat{};
std::array<std::shared_ptr<TH2>, NpartChrg> hPtNsigmaTPC{};
std::array<std::shared_ptr<THnSparse>, NpartChrg> hThPtNsigmaTPC{};
std::array<std::shared_ptr<TH2>, NpartChrg> hPtNsigmaTOF{};
std::array<std::shared_ptr<TH2>, NpartChrg> hPtNsigmaTPCTOF{};

struct FlattenictyPikp {

  HistogramRegistry flatchrg{"flatchrg", {}, OutputObjHandlingPolicy::AnalysisObject, true, false};
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
  Configurable<bool> cfgFilldEdxCalibHist{"cfgFilldEdxCalibHist", false, "fill dEdx calibration histograms"};
  Configurable<bool> cfgFilldEdxQaHist{"cfgFilldEdxQaHist", false, "fill dEdx QA histograms"};
  Configurable<bool> cfgFillNsigmaQAHist{"cfgFillNsigmaQAHist", false, "fill nsigma QA histograms"};
  Configurable<bool> cfgFillV0Hist{"cfgFillV0Hist", false, "fill V0 histograms"};
  Configurable<bool> cfgFillChrgType{"cfgFillChrgType", false, "fill histograms per charge types"};
  Configurable<bool> cfgFillChrgTypeV0s{"cfgFillChrgTypeV0s", false, "fill V0s histograms per charge types"};
  Configurable<std::vector<float>> paramsFuncMIPposEtaP{"paramsFuncMIPposEtaP", std::vector<float>{-1.f}, "parameters of pol2"};
  Configurable<std::vector<float>> paramsFuncMIPnegEtaP{"paramsFuncMIPnegEtaP", std::vector<float>{-1.f}, "parameters of pol2"};
  Configurable<std::vector<float>> paramsFuncMIPallEtaP{"paramsFuncMIPallEtaP", std::vector<float>{-1.f}, "parameters of pol2"};
  Configurable<std::vector<float>> paramsFuncMIPposEtaN{"paramsFuncMIPposEtaN", std::vector<float>{-1.f}, "parameters of pol2"};
  Configurable<std::vector<float>> paramsFuncMIPnegEtaN{"paramsFuncMIPnegEtaN", std::vector<float>{-1.f}, "parameters of pol2"};
  Configurable<std::vector<float>> paramsFuncMIPallEtaN{"paramsFuncMIPallEtaN", std::vector<float>{-1.f}, "parameters of pol2"};
  Configurable<std::vector<float>> paramsFuncPlateaUposEtaP{"paramsFuncPlateaUposEtaP", std::vector<float>{-1.f}, "parameters of pol2"};
  Configurable<std::vector<float>> paramsFuncPlateaUnegEtaP{"paramsFuncPlateaUnegEtaP", std::vector<float>{-1.f}, "parameters of pol2"};
  Configurable<std::vector<float>> paramsFuncPlateaUallEtaP{"paramsFuncPlateaUallEtaP", std::vector<float>{-1.f}, "parameters of pol2"};
  Configurable<std::vector<float>> paramsFuncPlateaUposEtaN{"paramsFuncPlateaUposEtaN", std::vector<float>{-1.f}, "parameters of pol2"};
  Configurable<std::vector<float>> paramsFuncPlateaUnegEtaN{"paramsFuncPlateaUnegEtaN", std::vector<float>{-1.f}, "parameters of pol2"};
  Configurable<std::vector<float>> paramsFuncPlateaUallEtaN{"paramsFuncPlateaUallEtaN", std::vector<float>{-1.f}, "parameters of pol2"};
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
    Configurable<bool> cfgINELCut{"cfgINELCut", true, "INEL event selection"};
    Configurable<bool> cfgRemoveNoSameBunchPileup{"cfgRemoveNoSameBunchPileup", true, "Reject collisions in case of pileup with another collision in the same foundBC"};
    Configurable<bool> cfgRequireIsGoodZvtxFT0vsPV{"cfgRequireIsGoodZvtxFT0vsPV", true, "Small difference between z-vertex from PV and from FT0"};
    Configurable<bool> cfgRequireIsVertexITSTPC{"cfgRequireIsVertexITSTPC", false, "At least one ITS-TPC track (reject vertices built from ITS-only tracks)"};
    Configurable<bool> cfgRequirekIsVertexTOFmatched{"cfgRequirekIsVertexTOFmatched", false, "Require kIsVertexTOFmatched: at least one of vertex contributors is matched to TOF"};
  } evtSelOpt;

  struct : ConfigurableGroup {
    Configurable<bool> useFlatData{"useFlatData", true, "use flattenicity from rec collisions"};
  } flatSelOpt;

  struct : ConfigurableGroup {
    ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0}, "pT binning"};
    ConfigurableAxis axisPtV0s{"axisPtV0s", {VARIABLE_WIDTH, 0, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.5, 3.0, 3.5, 4, 5, 7, 9, 12, 15, 20}, "pT V0s binning"};
    ConfigurableAxis axisFlatPerc{"axisFlatPerc", {102, -0.01, 1.01}, "Flattenicity percentiles binning"};
    ConfigurableAxis axisMultPerc{"axisMultPerc", {20, 0, 100}, "Multiplicity percentiles binning"};
    ConfigurableAxis axisVertexZ{"axisVertexZ", {80, -20., 20.}, "Vertex z binning"};
    ConfigurableAxis axisMult{"axisMult", {301, -0.5, 300.5}, "Multiplicity binning"};
    ConfigurableAxis axisDCAxy{"axisDCAxy", {200, -5, 5}, "DCAxy binning"};
    ConfigurableAxis axisDCAz{"axisDCAz", {200, -5, 5}, "DCAz binning"};
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
    Configurable<float> cfgMomSelPiTOF{"cfgMomSelPiTOF", 0.4f, "Minimum momentum cut for TOF pions"};
    Configurable<float> cfgNsigSelKaTOF{"cfgNsigSelKaTOF", 3.0f, "Nsigma cut for TOF kaons"};
    Configurable<float> cfgBetaPlateuMax{"cfgBetaPlateuMax", 0.1f, "Beta max for Plateau electrons"};
  } trkSelOpt;

  struct : ConfigurableGroup {
    // common selection
    Configurable<int> cfgV0TypeSel{"cfgV0TypeSel", 1, "select on a certain V0 type (leave negative if no selection desired)"};
    Configurable<float> cfgV0Ymax{"cfgV0Ymax", 0.5f, "Maximum rapidity of V0s"};
    Configurable<float> cfgPtDaughterMin{"cfgPtDaughterMin", 0.1f, "minimum pT of the V0 daughter tracks"};
    Configurable<float> cfgPtDaughterMax{"cfgPtDaughterMax", 20.0f, "maximum pT of the V0 daughter tracks"};
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
  Configurable<float> minPt{"minPt", 0.15f, "Set minimum pT of tracks"};
  Configurable<float> maxPt{"maxPt", 20.0f, "Set maximum pT of tracks"};
  Configurable<float> requireEta{"requireEta", 0.8f, "Set eta range of tracks"};
  Configurable<int> setITSreq{"setITSreq", 0, "0 = Run3ITSibAny, 1 = Run3ITSallAny, 2 = Run3ITSall7Layers, 3 = Run3ITSibTwo"};
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

    // Event counter
    flatchrg.add("Events/hEvtSel", "Number of events; Cut; #Events Passed Cut", {kTH1F, {{nEvtSel, 0, nEvtSel}}});
    flatchrg.get<TH1>(HIST("Events/hEvtSel"))->GetXaxis()->SetBinLabel(evtSelAll + 1, "Events read");
    flatchrg.get<TH1>(HIST("Events/hEvtSel"))->GetXaxis()->SetBinLabel(evtSelSel8 + 1, "Evt. sel8");
    flatchrg.get<TH1>(HIST("Events/hEvtSel"))->GetXaxis()->SetBinLabel(evtSelNoITSROFrameBorder + 1, "NoITSROFrameBorder");
    flatchrg.get<TH1>(HIST("Events/hEvtSel"))->GetXaxis()->SetBinLabel(evtSelkNoTimeFrameBorder + 1, "NoTimeFrameBorder");
    flatchrg.get<TH1>(HIST("Events/hEvtSel"))->GetXaxis()->SetBinLabel(evtSelkNoSameBunchPileup + 1, "NoSameBunchPileup");
    flatchrg.get<TH1>(HIST("Events/hEvtSel"))->GetXaxis()->SetBinLabel(evtSelkIsGoodZvtxFT0vsPV + 1, "IsGoodZvtxFT0vsPV");
    flatchrg.get<TH1>(HIST("Events/hEvtSel"))->GetXaxis()->SetBinLabel(evtSelkIsVertexITSTPC + 1, "IsVertexITSTPC");
    flatchrg.get<TH1>(HIST("Events/hEvtSel"))->GetXaxis()->SetBinLabel(evtSelkIsVertexTOFmatched + 1, "IsVertexTOFmatched");
    flatchrg.get<TH1>(HIST("Events/hEvtSel"))->GetXaxis()->SetBinLabel(evtSelVtxZ + 1, "Vtx-z pos");
    flatchrg.get<TH1>(HIST("Events/hEvtSel"))->GetXaxis()->SetBinLabel(evtSelINELgt0 + 1, "INEL>0");
    flatchrg.get<TH1>(HIST("Events/hEvtSel"))->GetXaxis()->SetBinLabel(evtSelRCTFlagChecker + 1, "RCT Flag Checker");
    // Track counter
    flatchrg.add("Tracks/hTrkSel", "Number of tracks; Cut; #Tracks Passed Cut", {kTH1F, {{nTrkSel, 0, nTrkSel}}});
    flatchrg.get<TH1>(HIST("Tracks/hTrkSel"))->GetXaxis()->SetBinLabel(trkSelAll + 1, "All");
    flatchrg.get<TH1>(HIST("Tracks/hTrkSel"))->GetXaxis()->SetBinLabel(trkSelEta + 1, "Eta");
    flatchrg.get<TH1>(HIST("Tracks/hTrkSel"))->GetXaxis()->SetBinLabel(trkSelPt + 1, "Pt");
    flatchrg.get<TH1>(HIST("Tracks/hTrkSel"))->GetXaxis()->SetBinLabel(trkSelDCA + 1, "DCA");
    flatchrg.get<TH1>(HIST("Tracks/hTrkSel"))->GetXaxis()->SetBinLabel(trkNRowsTPC + 1, "trkNRowsTPC");
    flatchrg.get<TH1>(HIST("Tracks/hTrkSel"))->GetXaxis()->SetBinLabel(trkSelNClsFound + 1, "NClsTPCFound");
    flatchrg.get<TH1>(HIST("Tracks/hTrkSel"))->GetXaxis()->SetBinLabel(trkSelNClsPID + 1, "NClsTPCPid");
    flatchrg.get<TH1>(HIST("Tracks/hTrkSel"))->GetXaxis()->SetBinLabel(trkSelTPCBndr + 1, "TPC Boundary");
    // V0 counter
    flatchrg.add("Tracks/V0qa/hV0Sel", "Number of V0s; Cut; #Tracks Passed Cut", {kTH1F, {{nV0Sel, 0, nV0Sel}}});
    flatchrg.get<TH1>(HIST("Tracks/V0qa/hV0Sel"))->GetXaxis()->SetBinLabel(v0SelAll + 1, "All");
    flatchrg.get<TH1>(HIST("Tracks/V0qa/hV0Sel"))->GetXaxis()->SetBinLabel(v0SelRejectSameSign + 1, "Reject same sign");
    flatchrg.get<TH1>(HIST("Tracks/V0qa/hV0Sel"))->GetXaxis()->SetBinLabel(v0SelRejectV0sAtTPCSector + 1, "Reject V0s at TPC sector");
    flatchrg.get<TH1>(HIST("Tracks/V0qa/hV0Sel"))->GetXaxis()->SetBinLabel(v0SelCosPA + 1, "Cos PA");
    flatchrg.get<TH1>(HIST("Tracks/V0qa/hV0Sel"))->GetXaxis()->SetBinLabel(v0SelV0radius + 1, "V0 radius");
    flatchrg.get<TH1>(HIST("Tracks/V0qa/hV0Sel"))->GetXaxis()->SetBinLabel(v0SelDCAposToPV + 1, "DCA pos to PV");
    flatchrg.get<TH1>(HIST("Tracks/V0qa/hV0Sel"))->GetXaxis()->SetBinLabel(v0SelDaughters + 1, "V0 daughters' sel.");
    flatchrg.get<TH1>(HIST("Tracks/V0qa/hV0Sel"))->GetXaxis()->SetBinLabel(v0SelDCAv0daughter + 1, "DCA v0 daughter");

    if (trkSelOpt.cfgRejectTrkAtTPCSector || v0SelOpt.cfgRejectV0sAtTPCSector) {
      fPhiCutLow = new TF1("fPhiCutLow", trkSelOpt.cfgGeoTrkCutMin.value.c_str(), 0, 100);
      fPhiCutHigh = new TF1("fPhiCutHigh", trkSelOpt.cfgGeoTrkCutMax.value.c_str(), 0, 100);
    }

    if (doprocessFlat) {
      flatchrg.add("Events/hVtxZ", "Measured vertex z position", kTH1F, {vtxzAxis});
      flatchrg.add("Events/hFlatVsMultEst", "hFlatVsMultEst", kTH2F, {flatAxis, multAxis});
      flatchrg.add("Tracks/postSel/hPVsPtEta", "; #it{p} (GeV/#it{c}); #it{p}_{T} (GeV/#it{c}); #eta;", {kTH3F, {pAxis, ptAxis, etaAxis}});
      if (cfgFillTrackQaHist || cfgFilldEdxQaHist || cfgFillNsigmaQAHist) {
        if (cfgFillTrackQaHist) {
          flatchrg.add("Tracks/postSel/hPtPhi", "; #it{p}_{T} (GeV/#it{c}); fmod(#varphi,#pi/9)", {kTH2F, {ptAxis, phiAxisMod}});
          flatchrg.add("Tracks/postSel/hPtVsWOcutDCA", "hPtVsWOcutDCA", kTH2F, {ptAxis, dcaXYAxis});
          flatchrg.add("Tracks/postSel/hPt", "", kTH1F, {ptAxis});
          flatchrg.add("Tracks/postSel/hPhi", "", kTH1F, {phiAxis});
          flatchrg.add("Tracks/postSel/hEta", "", kTH1F, {etaAxis});
          flatchrg.add("Tracks/postSel/hDCAXYvsPt", "", kTH2F, {ptAxis, dcaXYAxis});
          flatchrg.add("Tracks/postSel/hDCAZvsPt", "", kTH2F, {ptAxis, dcaZAxis});
          // tpc
          if (cfgStoreThnSparse) {
            flatchrg.add("Tracks/postSel/hPtPhiNclTPC", "; #it{p}_{T} (GeV/#it{c}); fmod(#varphi,#pi/9); N_{cluster}", {kTHnSparseF, {ptAxis, phiAxisMod, clTpcAxis}});
            flatchrg.add("Tracks/postSel/hPtPhiNclPIDTPC", "; #it{p}_{T} (GeV/#it{c}); fmod(#varphi,#pi/9); N_{PID cluster}", {kTHnSparseF, {ptAxis, phiAxisMod, clTpcAxis}});
          } else {
            flatchrg.add("Tracks/postSel/hPtPhiNclTPC", "; #it{p}_{T} (GeV/#it{c}); fmod(#varphi,#pi/9); N_{cluster}", {kTH3F, {ptAxis, phiAxisMod, clTpcAxis}});
            flatchrg.add("Tracks/postSel/hPtPhiNclPIDTPC", "; #it{p}_{T} (GeV/#it{c}); fmod(#varphi,#pi/9); N_{PID cluster}", {kTH3F, {ptAxis, phiAxisMod, clTpcAxis}});
          }
          flatchrg.add("Tracks/postSel/hPtNclTPC", "; #it{p}_{T} (GeV/#it{c}); N_{cluster}", {kTH2F, {ptAxis, clTpcAxis}});
          flatchrg.add("Tracks/postSel/pPtNclTPC", "; #it{p}_{T} (GeV/#it{c}); N_{cluster}", {kTProfile, {ptAxis}});
          flatchrg.add("Tracks/postSel/hPtNclPIDTPC", "; #it{p}_{T} (GeV/#it{c}); N_{PID cluster}", {kTH2F, {ptAxis, clTpcAxis}});
          flatchrg.add("Tracks/postSel/pPtNclPIDTPC", "; #it{p}_{T} (GeV/#it{c}); N_{PID cluster}", {kTProfile, {ptAxis}});
          flatchrg.add("Tracks/postSel/hShTpcClvsPt", "", {kTH2F, {ptAxis, shCluserAxis}});
          flatchrg.add("Tracks/postSel/hNclTPCFoundvsPt", "", {kTH2F, {ptAxis, clTpcAxis}});
          flatchrg.add("Tracks/postSel/hNClTPCPidvsPt", "", {kTH2F, {ptAxis, clTpcAxis}});
          flatchrg.add("Tracks/postSel/hTPCCluster", "N_{cluster}", kTH1F, {clTpcAxis});
          flatchrg.add("Tracks/postSel/hTPCnClsShared", " ; # shared TPC clusters TPC", kTH1F, {{165, -0.5, 164.5}});
          flatchrg.add("Tracks/postSel/hTPCcrossedRows", " ; # crossed TPC rows", kTH1F, {{165, -0.5, 164.5}});
          flatchrg.add("Tracks/postSel/hTPCcrossedRowsOverFindableCls", " ; crossed rows / findable TPC clusters", kTH1F, {{60, 0.7, 1.3}});
          // its
          flatchrg.add("Tracks/postSel/hITSnCls", " ; # ITS clusters", kTH1F, {{8, -0.5, 7.5}});
          flatchrg.add("Tracks/postSel/hChi2ITSTrkSegment", "chi2ITS", kTH1F, {{100, -0.5, 99.5}});
          // tof
          flatchrg.add("Tracks/postSel/hTOFPvsBeta", "Beta from TOF; #it{p} (GeV/#it{c}); #beta", {kTH2F, {pAxis, {120, 0.0, 1.2}}});
          if (cfgStoreThnSparse) {
            flatchrg.add("Tracks/postSel/hTOFpi", "Primary Pions from TOF; #eta; #it{p} (GeV/#it{c}); dEdx", {kTHnSparseF, {etaAxis, pAxis, dEdxAxis}});
          } else {
            flatchrg.add("Tracks/postSel/hTOFpi", "Primary Pions from TOF; #eta; #it{p} (GeV/#it{c}); dEdx", {kTH3F, {etaAxis, pAxis, dEdxAxis}});
          }
        }
        if (cfgFilldEdxQaHist) {
          if (cfgStoreThnSparse) {
            flatchrg.add("Tracks/postCalib/all/hMIP", "; mult; flat; #eta; #LT dE/dx #GT_{MIP, primary tracks};", {kTHnSparseF, {multAxis, flatAxis, etaAxis, dEdxAxis}});
            flatchrg.add("Tracks/postCalib/all/hPlateau", "; mult; flat; #eta; #LT dE/dx #GT_{Plateau, primary tracks};", {kTHnSparseF, {multAxis, flatAxis, etaAxis, dEdxAxis}});
          } else {
            flatchrg.add("Tracks/postCalib/all/hMIP", "; #eta; #LT dE/dx #GT_{MIP, primary tracks};", {kTH2F, {etaAxis, dEdxAxis}});
            flatchrg.add("Tracks/postCalib/all/hPlateau", "; #eta; #LT dE/dx #GT_{Plateau, primary tracks};", {kTH2F, {etaAxis, dEdxAxis}});
          }
          flatchrg.add("Tracks/postCalib/all/hMIPVsPhi", "; #varphi; #LT dE/dx #GT_{MIP, primary tracks};", {kTH2F, {phiAxis, dEdxAxis}});
          flatchrg.add("Tracks/postCalib/all/pMIPVsPhi", "; #varphi; #LT dE/dx #GT_{MIP, primary tracks};", {kTProfile, {phiAxis}});
          flatchrg.add("Tracks/postCalib/all/hMIPVsPhiVsEta", "; #varphi; #LT dE/dx #GT_{MIP, primary tracks}; #eta;", {kTH3F, {phiAxis, dEdxAxis, etaAxis}});
          flatchrg.add("Tracks/postCalib/all/hPlateauVsPhi", "; #varphi; #LT dE/dx #GT_{Plateau, primary tracks};", {kTH2F, {phiAxis, dEdxAxis}});
          flatchrg.add("Tracks/postCalib/all/pPlateauVsPhi", "; #varphi; #LT dE/dx #GT_{Plateau, primary tracks};", {kTProfile, {phiAxis}});
          flatchrg.add("Tracks/postCalib/all/hPlateauVsPhiVsEta", "; #varphi; #LT dE/dx #GT_{Plateau, primary tracks}; #eta;", {kTH3F, {phiAxis, dEdxAxis, etaAxis}});
          flatchrg.addClone("Tracks/postCalib/all/", "Tracks/preCalib/all/");
          if (cfgFillChrgType) {
            flatchrg.addClone("Tracks/postCalib/all/", "Tracks/postCalib/pos/");
            flatchrg.addClone("Tracks/postCalib/all/", "Tracks/postCalib/neg/");
            flatchrg.addClone("Tracks/preCalib/all/", "Tracks/preCalib/pos/");
            flatchrg.addClone("Tracks/preCalib/all/", "Tracks/preCalib/neg/");
          }
        }
        if (cfgFillNsigmaQAHist) {
          for (int i = 0; i < NpartChrg; i++) {
            const std::string strID = Form("/%s/%s", (i < Npart) ? "pos" : "neg", Pid[i % Npart]);
            if (cfgStoreThnSparse) {
              hThPtNsigmaTPC[i] = flatchrg.add<THnSparse>("Tracks/hThPtNsigmaTPC" + strID, " ; p_{T} (GeV/c)", kTHnSparseF, {ptAxis, nSigmaTPCAxis, multAxis, flatAxis});
            } else {
              hPtNsigmaTPC[i] = flatchrg.add<TH2>("Tracks/hPtNsigmaTPC" + strID, " ; p_{T} (GeV/c)", kTH2F, {ptAxis, nSigmaTPCAxis});
            }
            hPtNsigmaTOF[i] = flatchrg.add<TH2>("Tracks/hPtNsigmaTOF" + strID, " ; p_{T} (GeV/c)", kTH2F, {ptAxis, nSigmaTOFAxis});
            hPtNsigmaTPCTOF[i] = flatchrg.add<TH2>("Tracks/hPtNsigmaTPCTOF" + strID, PidChrg[i], kTH2F, {nSigmaTPCAxis, nSigmaTOFAxis});
          }
        }
      }
      flatchrg.addClone("Tracks/postSel/", "Tracks/preSel/");
      // FV0 QA
      flatchrg.add("FV0/hFV0AmplWCalib", "", {kTH2F, {channelFV0Axis, amplitudeFV0}});
      flatchrg.add("FV0/hFV0AmplvsVtxzWoCalib", "", {kTH2F, {vtxzAxis, amplitudeFV0Sum}});
      flatchrg.add("FV0/hFV0AmplvsVtxzCalib", "", {kTH2F, {vtxzAxis, amplitudeFV0Sum}});
      flatchrg.add("FV0/hFV0amp", "", {kTH2F, {channelFV0Axis, amplitudeFV0}});
      flatchrg.add("FV0/pFV0amp", "", kTProfile, {channelFV0Axis});
      flatchrg.add("FV0/hFV0ampCorr", "", {kTH2F, {channelFV0Axis, amplitudeFV0}});
      // V0's QA
      flatchrg.add("Tracks/V0qa/hV0Pt", "pT", kTH1F, {ptAxisV0s});
      flatchrg.add("Tracks/V0qa/hV0ArmPod", ";#alpha; #it{q}_T (GeV/c)", kTH2F, {v0SelOpt.axisArmPodAlpha, v0SelOpt.axisArmPodqT});
      // daughters' QA
      flatchrg.add("Tracks/V0qa/el/Ga/hArmPod", ";#alpha; #it{q}_T (GeV/c)", kTH2F, {v0SelOpt.axisArmPodAlpha, v0SelOpt.axisArmPodqT});
      flatchrg.add("Tracks/V0qa/pi/K0s/hArmPod", ";#alpha; #it{q}_T (GeV/c)", kTH2F, {v0SelOpt.axisArmPodAlpha, v0SelOpt.axisArmPodqT});
      flatchrg.add("Tracks/V0qa/el/Ga/hNclVsEta", ";#eta; #it{N}^{TPC}_cl", kTH2F, {etaAxis, clTpcAxis});
      flatchrg.add("Tracks/V0qa/pi/K0s/hNclVsEta", ";#eta; #it{N}^{TPC}_cl", kTH2F, {etaAxis, clTpcAxis});
      flatchrg.add("Tracks/V0qa/el/Ga/hdEdxMIPVsEta", ";#eta; dE/dx", kTH2F, {etaAxis, dEdxAxis});
      flatchrg.add("Tracks/V0qa/pi/K0s/hdEdxMIPVsEta", ";#eta; dE/dx", kTH2F, {etaAxis, dEdxAxis});
      flatchrg.addClone("Tracks/V0qa/pi/K0s/", "Tracks/V0qa/pi/La/");
      flatchrg.addClone("Tracks/V0qa/pi/K0s/", "Tracks/V0qa/pi/ALa/");
      flatchrg.addClone("Tracks/V0qa/pi/La/", "Tracks/V0qa/pr/La/");
      flatchrg.addClone("Tracks/V0qa/pi/ALa/", "Tracks/V0qa/pr/ALa/");

      // dEdx PID
      flatchrg.add({"Tracks/all/hdEdx", "; #eta; mult; flat; #it{p} (GeV/#it{c}); dEdx", {kTHnSparseF, {etaAxis, multAxis, flatAxis, pAxis, dEdxAxis}}});
      // Clean samples
      if (cfgFillV0Hist) {
        if (cfgStoreThnSparse) {
          flatchrg.add({"Tracks/CleanTof/all/hPiTof", "; #eta; mult; flat; #it{p} (GeV/#it{c}); dEdx", {kTHnSparseF, {etaAxis, multAxis, flatAxis, pAxis, dEdxAxis}}});
          flatchrg.add({"Tracks/CleanV0/all/hEV0", "; #eta; mult; flat; #it{p} (GeV/#it{c}); dEdx", {kTHnSparseF, {etaAxis, multAxis, flatAxis, pAxis, dEdxAxis}}});
          flatchrg.add({"Tracks/CleanV0/all/hPiV0", "; #eta; mult; flat; #it{p} (GeV/#it{c}); dEdx", {kTHnSparseF, {etaAxis, multAxis, flatAxis, pAxis, dEdxAxis}}});
          flatchrg.add({"Tracks/CleanV0/all/hPV0", "; #eta; mult; flat; #it{p} (GeV/#it{c}); dEdx", {kTHnSparseF, {etaAxis, multAxis, flatAxis, pAxis, dEdxAxis}}});
        } else {
          flatchrg.add({"Tracks/CleanTof/all/hPiTof", "; #eta; mult; flat; #it{p} (GeV/#it{c}); dEdx", {kTH3F, {etaAxis, pAxis, dEdxAxis}}});
          flatchrg.add({"Tracks/CleanV0/all/hEV0", "; #eta; mult; flat; #it{p} (GeV/#it{c}); dEdx", {kTH3F, {etaAxis, pAxis, dEdxAxis}}});
          flatchrg.add({"Tracks/CleanV0/all/hPiV0", "; #eta; mult; flat; #it{p} (GeV/#it{c}); dEdx", {kTH3F, {etaAxis, pAxis, dEdxAxis}}});
          flatchrg.add({"Tracks/CleanV0/all/hPV0", "; #eta; mult; flat; #it{p} (GeV/#it{c}); dEdx", {kTH3F, {etaAxis, pAxis, dEdxAxis}}});
        }
        flatchrg.add("Tracks/CleanTof/all/hBetaVsP", ";Momentum (GeV/#it{c}); #beta", kTH2F, {{{ptAxisV0s}, {120, 0., 1.2}}});
        flatchrg.add("Tracks/CleanTof/all/hTofExpPi", ";Momentum (GeV/#it{c});#it{t}^{#pi}_{Exp}/#it{t}_{TOF}", kTH2F, {{{ptAxisV0s}, {100, 0.2, 1.2}}});
        if (cfgFillChrgType) {
          flatchrg.addClone("Tracks/CleanTof/all/", "Tracks/CleanTof/pos/");
          flatchrg.addClone("Tracks/CleanTof/all/", "Tracks/CleanTof/neg/");
          flatchrg.addClone("Tracks/CleanV0/all/", "Tracks/CleanV0/pos/");
          flatchrg.addClone("Tracks/CleanV0/all/", "Tracks/CleanV0/neg/");
        }
      }
      if (cfgFillChrgType) {
        flatchrg.addClone("Tracks/all/", "Tracks/pos/");
        flatchrg.addClone("Tracks/all/", "Tracks/neg/");
      }
    }

    if (doprocessMC) {
      auto h = flatchrg.add<TH1>("hEvtGenRec", "Generated and Reconstructed MC Collisions", kTH1F, {{3, 0.5, 3.5}});
      h->GetXaxis()->SetBinLabel(1, "Gen coll");
      h->GetXaxis()->SetBinLabel(2, "Rec coll");
      h->GetXaxis()->SetBinLabel(3, "INEL>0");

      flatchrg.add("hEvtMcGenColls", "Number of events; Cut; #Events Passed Cut", {kTH1F, {{5, 0.5, 5.5}}});
      flatchrg.get<TH1>(HIST("hEvtMcGenColls"))->GetXaxis()->SetBinLabel(1, "Gen. coll");
      flatchrg.get<TH1>(HIST("hEvtMcGenColls"))->GetXaxis()->SetBinLabel(2, "At least 1 reco");
      flatchrg.get<TH1>(HIST("hEvtMcGenColls"))->GetXaxis()->SetBinLabel(3, "Reco. coll.");
      flatchrg.get<TH1>(HIST("hEvtMcGenColls"))->GetXaxis()->SetBinLabel(4, "Reco. good coll.");

      flatchrg.add("Events/hVtxZRec", "MC Rec vertex z position", kTH1F, {vtxzAxis});
      flatchrg.add("Events/hVtxZGen", "Generated vertex z position", kTH1F, {vtxzAxis});

      for (int i = 0; i < NpartChrg; i++) {
        const std::string strID = Form("/%s/%s", (i < Npart) ? "pos" : "neg", Pid[i % Npart]);
        hPtGenRecEvt[i] = flatchrg.add<TH1>("Tracks/hPtGenRecEvt" + strID, " ; p_{T} (GeV/c)", kTH1F, {ptAxis});
        hPtGenPrimRecEvt[i] = flatchrg.add<TH1>("Tracks/hPtGenPrimRecEvt" + strID, " ; p_{T} (GeV/c)", kTH1F, {ptAxis});
        hPtEffGenPrim[i] = flatchrg.add<TH3>("Tracks/hPtEffGenPrim" + strID, " ; p_{T} (GeV/c)", kTH3F, {multAxis, flatAxis, ptAxis});
        hPtEffGenWeak[i] = flatchrg.add<TH3>("Tracks/hPtEffGenWeak" + strID, " ; p_{T} (GeV/c)", kTH3F, {multAxis, flatAxis, ptAxis});
        hPtEffGenMat[i] = flatchrg.add<TH3>("Tracks/hPtEffGenMat" + strID, " ; p_{T} (GeV/c)", kTH3F, {multAxis, flatAxis, ptAxis});
        hPtEffRecPrim[i] = flatchrg.add<TH3>("Tracks/hPtEffRecPrim" + strID, " ; p_{T} (GeV/c)", kTH3F, {multAxis, flatAxis, ptAxis});
        hPtEffRecWeak[i] = flatchrg.add<TH3>("Tracks/hPtEffRecWeak" + strID, " ; p_{T} (GeV/c)", kTH3F, {multAxis, flatAxis, ptAxis});
        hPtEffRecMat[i] = flatchrg.add<TH3>("Tracks/hPtEffRecMat" + strID, " ; p_{T} (GeV/c)", kTH3F, {multAxis, flatAxis, ptAxis});
        hDCAxyBadCollPrim[i] = flatchrg.add<TH2>("Tracks/hDCAxyBadCollPrim" + strID, " ; p_{T} (GeV/c)", kTH2F, {ptAxis, dcaXYAxis});
        hDCAxyBadCollWeak[i] = flatchrg.add<TH2>("Tracks/hDCAxyBadCollWeak" + strID, " ; p_{T} (GeV/c)", kTH2F, {ptAxis, dcaXYAxis});
        hDCAxyBadCollMat[i] = flatchrg.add<TH2>("Tracks/hDCAxyBadCollMat" + strID, " ; p_{T} (GeV/c)", kTH2F, {ptAxis, dcaXYAxis});
        hPtVsDCAxyPrim[i] = flatchrg.add<TH2>("Tracks/hPtVsDCAxyPrim" + strID, " ; p_{T} (GeV/c)", kTH2F, {ptAxis, dcaXYAxis});
        hPtVsDCAxyWeak[i] = flatchrg.add<TH2>("Tracks/hPtVsDCAxyWeak" + strID, " ; p_{T} (GeV/c)", kTH2F, {ptAxis, dcaXYAxis});
        hPtVsDCAxyMat[i] = flatchrg.add<TH2>("Tracks/hPtVsDCAxyMat" + strID, " ; p_{T} (GeV/c)", kTH2F, {ptAxis, dcaXYAxis});
      }

      flatchrg.add({"hPtOut", " ; p_{T} (GeV/c)", {kTH1F, {ptAxis}}});
      flatchrg.add({"hPtOutPrim", " ; p_{T} (GeV/c)", {kTH1F, {ptAxis}}});
      flatchrg.add({"hPtOutNoEtaCut", " ; p_{T} (GeV/c)", {kTH1F, {ptAxis}}});
      flatchrg.add({"PtOutFakes", " ; p_{T} (GeV/c)", {kTH1F, {ptAxis}}});
      flatchrg.add("hPtVsDCAxyPrimAll", "hPtVsDCAxyPrimAll", kTH2F, {ptAxis, dcaXYAxis});
      flatchrg.add("hPtVsDCAxyWeakAll", "hPtVsDCAxyWeakAll", kTH2F, {ptAxis, dcaXYAxis});
      flatchrg.add("hPtVsDCAxyMatAll", "hPtVsDCAxyMatAll", kTH2F, {ptAxis, dcaXYAxis});
      flatchrg.add("hPtVsDCAxyAll", "hPtVsDCAxyAll", kTH2F, {ptAxis, dcaXYAxis});
      flatchrg.add({"ResponseGen", " ; N_{part}; F_{FV0};", {kTHnSparseF, {multAxis, flatAxis}}});
      flatchrg.add("h1flatencityFV0MCGen", "", kTH1F, {{102, -0.01, 1.01, "1-flatencityFV0"}});

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
    }

    if (doprocessMCclosure) {
      for (int i = 0; i < Npart; i++) {
        flatchrg.add({fmt::format(CpTmcClosurePrimF.data(), CspeciesAll[i]).c_str(), " ; p_{T} (GeV/c)", {kTH3F, {multAxis, flatAxis, ptAxis}}});
      }
    }

    if (doprocessSgnLoss) {
      flatchrg.add("hFlatMCGenRecColl", "hFlatMCGenRecColl", {kTH1F, {flatAxis}});
      flatchrg.add("hFlatMCGen", "hFlatMCGen", {kTH1F, {flatAxis}});
      // Event counter
      flatchrg.add("hEvtMcGen", "hEvtMcGen", {kTH1F, {{4, 0.f, 4.f}}});
      flatchrg.get<TH1>(HIST("hEvtMcGen"))->GetXaxis()->SetBinLabel(1, "all");
      flatchrg.get<TH1>(HIST("hEvtMcGen"))->GetXaxis()->SetBinLabel(2, "z-vtx");
      flatchrg.get<TH1>(HIST("hEvtMcGen"))->GetXaxis()->SetBinLabel(3, "INELgt0");
      flatchrg.add("hEvtMCRec", "hEvtMCRec", {kTH1F, {{4, 0.f, 4.f}}});
      flatchrg.get<TH1>(HIST("hEvtMCRec"))->GetXaxis()->SetBinLabel(1, "all");
      flatchrg.get<TH1>(HIST("hEvtMCRec"))->GetXaxis()->SetBinLabel(2, "evt sel");
      flatchrg.get<TH1>(HIST("hEvtMCRec"))->GetXaxis()->SetBinLabel(3, "INELgt0");
      flatchrg.add("hEvtMcGenRecColl", "hEvtMcGenRecColl", {kTH1F, {{2, 0.f, 2.f}}});
      flatchrg.get<TH1>(HIST("hEvtMcGenRecColl"))->GetXaxis()->SetBinLabel(1, "INEL");
      flatchrg.get<TH1>(HIST("hEvtMcGenRecColl"))->GetXaxis()->SetBinLabel(2, "INELgt0");

      flatchrg.add("hFlatGenINELgt0", "hFlatGenINELgt0", {kTH1F, {flatAxis}});

      for (int i = 0; i < NpartChrg; ++i) {
        flatchrg.add({fmt::format(CpTgenPrimSgnF.data(), Cspecies[i]).c_str(), " ; p_{T} (GeV/c)", {kTH3F, {multAxis, flatAxis, ptAxis}}});
        flatchrg.add({fmt::format(CpTgenPrimSgnINELF.data(), Cspecies[i]).c_str(), " ; p_{T} (GeV/c)", {kTH3F, {multAxis, flatAxis, ptAxis}}});
        flatchrg.add({fmt::format(CpTrecCollPrimSgnF.data(), Cspecies[i]).c_str(), " ; p_{T} (GeV/c)", {kTH3F, {multAxis, flatAxis, ptAxis}}});
        flatchrg.add({fmt::format(CpTrecCollPrimSgnINELF.data(), Cspecies[i]).c_str(), " ; p_{T} (GeV/c)", {kTH3F, {multAxis, flatAxis, ptAxis}}});
        flatchrg.add({fmt::format(CpTGenRecCollPrimSgnF.data(), Cspecies[i]).c_str(), " ; p_{T} (GeV/c)", {kTH3F, {multAxis, flatAxis, ptAxis}}});
        flatchrg.add({fmt::format(CpTGenRecCollPrimSgnINELF.data(), Cspecies[i]).c_str(), " ; p_{T} (GeV/c)", {kTH3F, {multAxis, flatAxis, ptAxis}}});
      }
    }

    LOG(info) << "Size of the histograms:";
    flatchrg.print();
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
    std::unique_ptr<TF1> fCalibFunc(new TF1("fCalibFunc", "pol2", -1., 1.));
    if (vecPars.size() >= 1) {
      for (typename T::size_type i = 0; i < vecPars.size(); i++) {
        fCalibFunc->SetParameter(i, vecPars[i]);
      }
    }
    return fCalibFunc;
  }

  template <int pidSgn, o2::track::PID::ID id, typename P>
  bool isPID(const P& mcParticle)
  {
    static_assert(pidSgn == Cnull || pidSgn == 1);
    static_assert(id > Cnull || id < Npart);
    constexpr int Cidx = id + pidSgn * Npart;
    return mcParticle.pdgCode() == PidSgn[Cidx];
  }

  template <typename T>
  bool selTPCtrack(T const& track)
  {
    return (track.hasTPC() &&
            track.passedTPCCrossedRowsOverNCls() &&
            track.passedTPCChi2NDF() &&
            track.passedTPCNCls() &&
            track.passedTPCCrossedRows() &&
            track.passedTPCRefit());
  }

  template <ChargeType chrg, typename T>
  bool selTOFPi(T const& track)
  {
    if (track.hasTOF() && track.goodTOFMatch()) {
      const float tTOF = track.tofSignal();
      const float trkLength = track.length();
      const float tExpPiTOF = track.tofExpSignalPi(tTOF);
      if (track.p() >= trkSelOpt.cfgMomSelPiTOF && trkLength > Cnull && tTOF > Cnull) {
        flatchrg.fill(HIST(CprefixCleanTof) + HIST(Ccharge[chrg]) + HIST("hTofExpPi"), track.p(), tExpPiTOF / tTOF);
        if (std::abs((tExpPiTOF / tTOF) - Cone) < trkSelOpt.cfgTofBetaPiMax) {
          flatchrg.fill(HIST(CprefixCleanTof) + HIST(Ccharge[chrg]) + HIST("hBetaVsP"), track.p(), track.beta());
          // if (std::abs(track.tpcNSigmaPi()) < v0SelOpt.cfgNsigmaPiTPC && std::abs(track.tofNSigmaPi()) < v0SelOpt.cfgNsigmaPiTOF) {
          return true;
          // }
        }
      }
    }
    return false;
  }

  template <o2::track::PID::ID id, typename T, typename C>
  void fillNsigma(T const& tracks, C const& collision)
  {
    const float mult = getMult(collision);
    const float flat = fillFlat<false>(collision);
    for (const auto& track : tracks) {
      checkNsigma<id>(track, mult, flat);
    }
  }

  template <typename T, typename V, typename C>
  void filldEdx(T const& tracks, V const& v0s, C const& collision, aod::BCsWithTimestamps const& /*bcs*/)
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
    flatchrg.fill(HIST("Events/hFlatVsMultEst"), flat, mult);

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
          flatchrg.fill(HIST(Cprefix) + HIST(Ccharge[kPos]) + HIST("hdEdx"), track.eta(), mult, flat, track.p(), dEdx);
        } else {
          flatchrg.fill(HIST(Cprefix) + HIST(Ccharge[kNeg]) + HIST("hdEdx"), track.eta(), mult, flat, track.p(), dEdx);
        }
      } else {
        flatchrg.fill(HIST(Cprefix) + HIST(Ccharge[kAll]) + HIST("hdEdx"), track.eta(), mult, flat, track.p(), dEdx);
      }

      // TOF pions
      if (cfgFillV0Hist) {
        if (selTOFPi<kAll>(track)) {
          if (cfgFillChrgType) {
            if (track.sign() * track.p() > Cnull) {
              if (cfgStoreThnSparse) {
                flatchrg.fill(HIST(CprefixCleanTof) + HIST(Ccharge[kPos]) + HIST("hPiTof"), track.eta(), mult, flat, track.p(), dEdx);
              } else {
                flatchrg.fill(HIST(CprefixCleanTof) + HIST(Ccharge[kPos]) + HIST("hPiTof"), track.eta(), track.p(), dEdx);
              }
            } else {
              if (cfgStoreThnSparse) {
                flatchrg.fill(HIST(CprefixCleanTof) + HIST(Ccharge[kNeg]) + HIST("hPiTof"), track.eta(), mult, flat, track.p(), dEdx);
              } else {
                flatchrg.fill(HIST(CprefixCleanTof) + HIST(Ccharge[kNeg]) + HIST("hPiTof"), track.eta(), track.p(), dEdx);
              }
            }
          } else {
            if (cfgStoreThnSparse) {
              flatchrg.fill(HIST(CprefixCleanTof) + HIST(Ccharge[kAll]) + HIST("hPiTof"), track.eta(), mult, flat, track.p(), dEdx);
            } else {
              flatchrg.fill(HIST(CprefixCleanTof) + HIST(Ccharge[kAll]) + HIST("hPiTof"), track.eta(), track.p(), dEdx);
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
              flatchrg.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kPos]) + HIST("hEV0"), posTrack.eta(), mult, flat, posTrack.sign() * posTrack.p(), dEdxPos);
              flatchrg.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kNeg]) + HIST("hEV0"), negTrack.eta(), mult, flat, negTrack.sign() * negTrack.p(), dEdxNeg);
            } else {
              flatchrg.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kAll]) + HIST("hEV0"), posTrack.eta(), mult, flat, posTrack.sign() * posTrack.p(), dEdxPos);
              flatchrg.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kAll]) + HIST("hEV0"), negTrack.eta(), mult, flat, negTrack.sign() * negTrack.p(), dEdxNeg);
            }
          } else {
            if (cfgFillChrgType) {
              flatchrg.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kPos]) + HIST("hEV0"), posTrack.eta(), posTrack.sign() * posTrack.p(), dEdxPos);
              flatchrg.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kNeg]) + HIST("hEV0"), negTrack.eta(), negTrack.sign() * negTrack.p(), dEdxNeg);
            } else {
              flatchrg.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kAll]) + HIST("hEV0"), posTrack.eta(), posTrack.sign() * posTrack.p(), dEdxPos);
              flatchrg.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kAll]) + HIST("hEV0"), negTrack.eta(), negTrack.sign() * negTrack.p(), dEdxNeg);
            }
          }
        }
        if (selectTypeV0s(collision, v0, posTrack, negTrack) == kKz) { // K0S -> pi + pi
          fillV0QA<kPi, kKz>(v0, posTrack);
          fillV0QA<kPi, kKz>(v0, negTrack);
          if (cfgStoreThnSparse) {
            if (cfgFillChrgType) {
              flatchrg.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kPos]) + HIST("hPiV0"), posTrack.eta(), mult, flat, posTrack.sign() * posTrack.p(), dEdxPos);
              flatchrg.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kNeg]) + HIST("hPiV0"), negTrack.eta(), mult, flat, negTrack.sign() * negTrack.p(), dEdxNeg);
            } else {
              flatchrg.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kAll]) + HIST("hPiV0"), posTrack.eta(), mult, flat, posTrack.sign() * posTrack.p(), dEdxPos);
              flatchrg.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kAll]) + HIST("hPiV0"), negTrack.eta(), mult, flat, negTrack.sign() * negTrack.p(), dEdxNeg);
            }
          } else {
            if (cfgFillChrgType) {
              flatchrg.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kPos]) + HIST("hPiV0"), posTrack.eta(), posTrack.sign() * posTrack.p(), dEdxPos);
              flatchrg.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kNeg]) + HIST("hPiV0"), negTrack.eta(), negTrack.sign() * negTrack.p(), dEdxNeg);
            } else {
              flatchrg.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kAll]) + HIST("hPiV0"), posTrack.eta(), posTrack.sign() * posTrack.p(), dEdxPos);
              flatchrg.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kAll]) + HIST("hPiV0"), negTrack.eta(), negTrack.sign() * negTrack.p(), dEdxNeg);
            }
          }
        }
        if (selectTypeV0s(collision, v0, posTrack, negTrack) == kLam) { // L -> p + pi-
          fillV0QA<kPi, kLam>(v0, negTrack);
          fillV0QA<kPr, kLam>(v0, posTrack);
          if (cfgStoreThnSparse) {
            if (cfgFillChrgType) {
              flatchrg.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kPos]) + HIST("hPV0"), posTrack.eta(), mult, flat, posTrack.sign() * posTrack.p(), dEdxPos);
              flatchrg.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kNeg]) + HIST("hPiV0"), negTrack.eta(), mult, flat, negTrack.sign() * negTrack.p(), dEdxNeg);
            } else {
              flatchrg.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kAll]) + HIST("hPV0"), posTrack.eta(), mult, flat, posTrack.sign() * posTrack.p(), dEdxPos);
              flatchrg.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kAll]) + HIST("hPiV0"), negTrack.eta(), mult, flat, negTrack.sign() * negTrack.p(), dEdxNeg);
            }
          } else {
            if (cfgFillChrgType) {
              flatchrg.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kPos]) + HIST("hPV0"), posTrack.eta(), posTrack.sign() * posTrack.p(), dEdxPos);
              flatchrg.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kNeg]) + HIST("hPiV0"), negTrack.eta(), negTrack.sign() * negTrack.p(), dEdxNeg);
            } else {
              flatchrg.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kAll]) + HIST("hPV0"), posTrack.eta(), posTrack.sign() * posTrack.p(), dEdxPos);
              flatchrg.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kAll]) + HIST("hPiV0"), negTrack.eta(), negTrack.sign() * negTrack.p(), dEdxNeg);
            }
          }
        }
        if (selectTypeV0s(collision, v0, posTrack, negTrack) == kaLam) { // antiLambda -> pbar + pi+
          fillV0QA<kPi, kaLam>(v0, posTrack);
          fillV0QA<kPr, kaLam>(v0, negTrack);
          if (cfgStoreThnSparse) {
            if (cfgFillChrgType) {
              flatchrg.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kPos]) + HIST("hPiV0"), posTrack.eta(), mult, flat, posTrack.sign() * posTrack.p(), dEdxPos);
              flatchrg.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kNeg]) + HIST("hPV0"), negTrack.eta(), mult, flat, negTrack.sign() * negTrack.p(), dEdxNeg);
            } else {
              flatchrg.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kAll]) + HIST("hPiV0"), posTrack.eta(), mult, flat, posTrack.sign() * posTrack.p(), dEdxPos);
              flatchrg.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kAll]) + HIST("hPV0"), negTrack.eta(), mult, flat, negTrack.sign() * negTrack.p(), dEdxNeg);
            }
          } else {
            if (cfgFillChrgType) {
              flatchrg.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kPos]) + HIST("hPiV0"), posTrack.eta(), posTrack.sign() * posTrack.p(), dEdxPos);
              flatchrg.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kNeg]) + HIST("hPV0"), negTrack.eta(), negTrack.sign() * negTrack.p(), dEdxNeg);
            } else {
              flatchrg.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kAll]) + HIST("hPiV0"), posTrack.eta(), posTrack.sign() * posTrack.p(), dEdxPos);
              flatchrg.fill(HIST(CprefixCleanV0) + HIST(Ccharge[kAll]) + HIST("hPV0"), negTrack.eta(), negTrack.sign() * negTrack.p(), dEdxNeg);
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

  template <typename T>
  bool phiCut(T const& track, float mag, TF1* fphiCutLow, TF1* fphiCutHigh)
  {
    if (track.pt() < trkSelOpt.cfgPhiCutPtMin)
      return true;
    // cut to remove tracks at TPC boundaries
    double phimodn = track.phi();
    if (mag < Cnull) // for negative polarity field
      phimodn = o2::constants::math::TwoPI - phimodn;
    if (track.sign() < Cnull) // for negative charge
      phimodn = o2::constants::math::TwoPI - phimodn;
    if (phimodn < Cnull)
      LOGF(warning, "phi < Cnull: %g", phimodn);

    phimodn += o2::constants::math::PI / 18.0f; // to center gap in the middle
    phimodn = std::fmod(phimodn, o2::constants::math::PI / 9.0f);

    if (cfgFillTrackQaHist) {
      flatchrg.fill(HIST("Tracks/preSel/hPtPhi"), track.pt(), phimodn);
      if (track.hasTPC() && track.hasITS()) {
        if (cfgStoreThnSparse) {
          flatchrg.fill(HIST("Tracks/preSel/hPtPhiNclTPC"), track.pt(), phimodn, track.tpcNClsFound());
          flatchrg.fill(HIST("Tracks/preSel/hPtPhiNclPIDTPC"), track.pt(), phimodn, track.tpcNClsPID());
        } else {
          flatchrg.fill(HIST("Tracks/preSel/hPtPhiNclTPC"), track.pt(), phimodn, track.tpcNClsFound());
          flatchrg.fill(HIST("Tracks/preSel/hPtPhiNclPIDTPC"), track.pt(), phimodn, track.tpcNClsPID());
          flatchrg.fill(HIST("Tracks/preSel/hPtNclTPC"), track.pt(), track.tpcNClsFound());
          flatchrg.fill(HIST("Tracks/preSel/pPtNclTPC"), track.pt(), track.tpcNClsFound());
          flatchrg.fill(HIST("Tracks/preSel/hPtNclPIDTPC"), track.pt(), track.tpcNClsPID());
          flatchrg.fill(HIST("Tracks/preSel/pPtNclPIDTPC"), track.pt(), track.tpcNClsPID());
        }
      }
    }
    if (phimodn < fphiCutHigh->Eval(track.pt()) && phimodn > fphiCutLow->Eval(track.pt())) {
      return false;
    }
    if (cfgFillTrackQaHist) {
      flatchrg.fill(HIST("Tracks/postSel/hPtPhi"), track.pt(), phimodn);
      if (track.hasTPC() && track.hasITS()) {
        if (cfgStoreThnSparse) {
          flatchrg.fill(HIST("Tracks/postSel/hPtPhiNclTPC"), track.pt(), phimodn, track.tpcNClsFound());
          flatchrg.fill(HIST("Tracks/postSel/hPtPhiNclPIDTPC"), track.pt(), phimodn, track.tpcNClsPID());
        } else {
          flatchrg.fill(HIST("Tracks/postSel/hPtPhiNclTPC"), track.pt(), phimodn, track.tpcNClsFound());
          flatchrg.fill(HIST("Tracks/postSel/hPtPhiNclPIDTPC"), track.pt(), phimodn, track.tpcNClsPID());
          flatchrg.fill(HIST("Tracks/postSel/hPtNclTPC"), track.pt(), track.tpcNClsFound());
          flatchrg.fill(HIST("Tracks/postSel/pPtNclTPC"), track.pt(), track.tpcNClsFound());
          flatchrg.fill(HIST("Tracks/postSel/hPtNclPIDTPC"), track.pt(), track.tpcNClsPID());
          flatchrg.fill(HIST("Tracks/postSel/pPtNclPIDTPC"), track.pt(), track.tpcNClsPID());
        }
      }
    }
    return true;
  }

  template <typename T>
  bool isGoodTrack(T const& track, int const magfield)
  {
    flatchrg.fill(HIST("Tracks/hTrkSel"), trkSelAll);
    if (std::abs(track.eta()) > trkSelOpt.cfgTrkEtaMax) {
      return false;
    }
    flatchrg.fill(HIST("Tracks/hTrkSel"), trkSelEta);
    if (track.pt() < trkSelOpt.cfgTrkPtMin) {
      return false;
    }
    flatchrg.fill(HIST("Tracks/hTrkSel"), trkSelPt);
    if (!isDCAxyCut(track)) {
      return false;
    }
    flatchrg.fill(HIST("Tracks/hTrkSel"), trkSelDCA);
    if (track.tpcNClsCrossedRows() < minNCrossedRowsTPC) {
      return false;
    }
    flatchrg.fill(HIST("Tracks/hTrkSel"), trkNRowsTPC);
    if (trkSelOpt.cfgApplyNcl && track.tpcNClsFound() < trkSelOpt.cfgNclTPCMin) {
      return false;
    }
    flatchrg.fill(HIST("Tracks/hTrkSel"), trkSelNClsFound);

    if (trkSelOpt.cfgApplyNclPID && track.tpcNClsPID() < trkSelOpt.cfgNclPidTPCMin) {
      return false;
    }
    flatchrg.fill(HIST("Tracks/hTrkSel"), trkSelNClsPID);
    if (trkSelOpt.cfgRejectTrkAtTPCSector && !phiCut(track, magfield, fPhiCutLow, fPhiCutHigh)) {
      return false;
    }
    flatchrg.fill(HIST("Tracks/hTrkSel"), trkSelTPCBndr);
    return true;
  }

  template <int id, int typeMother, typename V, typename U>
  void fillV0QA(V const& v0, U const& track)
  {
    flatchrg.fill(HIST(CprefixV0qa) + HIST(PidDir[id]) + HIST(V0Dir[typeMother]) + HIST("hArmPod"), v0.alpha(), v0.qtarm());
    flatchrg.fill(HIST(CprefixV0qa) + HIST(PidDir[id]) + HIST(V0Dir[typeMother]) + HIST("hNclVsEta"), track.eta(), track.tpcNClsPID());
    flatchrg.fill(HIST(CprefixV0qa) + HIST(PidDir[id]) + HIST(V0Dir[typeMother]) + HIST("hdEdxMIPVsEta"), track.eta(), track.tpcSignal());
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
  bool isGoodV0Track(T1 const& v0, T2 const& /*track*/, int const magfield)
  {
    const auto& posTrack = v0.template posTrack_as<T2>();
    const auto& negTrack = v0.template negTrack_as<T2>();

    flatchrg.fill(HIST("Tracks/V0qa/hV0Sel"), v0SelAll);
    if (posTrack.sign() * negTrack.sign() > Cnull) {
      return false;
    }
    flatchrg.fill(HIST("Tracks/V0qa/hV0Sel"), v0SelRejectSameSign);
    if (v0SelOpt.cfgRejectV0sAtTPCSector) {
      if (!(phiCut(posTrack, magfield, fPhiCutLow, fPhiCutHigh) && phiCut(negTrack, magfield, fPhiCutLow, fPhiCutHigh))) {
        return false;
      }
    }
    flatchrg.fill(HIST("Tracks/V0qa/hV0Sel"), v0SelRejectV0sAtTPCSector);
    // V0 topological selections
    if (v0.v0cosPA() < v0SelOpt.cfgv0cospa) {
      return false;
    }
    flatchrg.fill(HIST("Tracks/V0qa/hV0Sel"), v0SelCosPA);
    if (v0.v0radius() < v0SelOpt.cfgv0Rmin || v0.v0radius() > v0SelOpt.cfgv0Rmax) {
      return false;
    }
    flatchrg.fill(HIST("Tracks/V0qa/hV0Sel"), v0SelV0radius);
    if (std::abs(v0.dcapostopv()) < v0SelOpt.cfgDCAposToPV || std::abs(v0.dcanegtopv()) < v0SelOpt.cfgDCAnegToPV) {
      return false;
    }
    flatchrg.fill(HIST("Tracks/V0qa/hV0Sel"), v0SelDCAposToPV);
    // selection of V0 daughters
    if (!(isGoodV0DaughterTrack(posTrack) && isGoodV0DaughterTrack(negTrack))) {
      return false;
    }
    flatchrg.fill(HIST("Tracks/V0qa/hV0Sel"), v0SelDaughters);
    if (v0.dcaV0daughters() > v0SelOpt.cfgDCAv0daughter) {
      return false;
    }
    flatchrg.fill(HIST("Tracks/V0qa/hV0Sel"), v0SelDCAv0daughter);
    if constexpr (fillHist) {
      flatchrg.fill(HIST("Tracks/V0qa/hV0Pt"), v0.pt());
      flatchrg.fill(HIST("Tracks/V0qa/hV0ArmPod"), v0.alpha(), v0.qtarm());
    }
    return true;
  }

  template <typename T>
  bool isGoodV0DaughterTrack(const T& track)
  {
    if (track.eta() < v0SelOpt.cfgV0etamin || track.eta() > v0SelOpt.cfgV0etamax) {
      return false;
    }
    if (track.pt() < v0SelOpt.cfgPtDaughterMin || track.pt() > v0SelOpt.cfgPtDaughterMax) {
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

  template <o2::track::PID::ID pid, typename T>
  void checkNsigma(const T& track, const float mult, const float flat)
  {
    if (std::abs(track.rapidity(o2::track::PID::getMass(pid))) > trkSelOpt.cfgRapMax) {
      return;
    }

    float valTPCnsigma = -999, valTOFnsigma = -999;
    switch (pid) {
      case o2::track::PID::Pion:
        valTPCnsigma = track.tpcNSigmaPi();
        valTOFnsigma = track.tofNSigmaPi();
        break;
      case o2::track::PID::Kaon:
        valTPCnsigma = track.tpcNSigmaKa();
        valTOFnsigma = track.tofNSigmaKa();
        break;
      case o2::track::PID::Proton:
        valTPCnsigma = track.tpcNSigmaPr();
        valTOFnsigma = track.tofNSigmaPr();
        break;
      case o2::track::PID::Electron:
        valTPCnsigma = track.tpcNSigmaEl();
        valTOFnsigma = track.tofNSigmaEl();
        break;
      case o2::track::PID::Muon:
        valTPCnsigma = track.tpcNSigmaMu();
        valTOFnsigma = track.tofNSigmaMu();
        break;
      default:
        valTPCnsigma = -999, valTOFnsigma = -999;
        break;
    }

    if (track.sign() > Cnull) {
      if (cfgStoreThnSparse) {
        hThPtNsigmaTPC[pid]->Fill(track.pt(), valTPCnsigma, mult, flat);
      } else {
        hPtNsigmaTPC[pid]->Fill(track.pt(), valTPCnsigma);
      }
    } else {
      if (cfgStoreThnSparse) {
        hThPtNsigmaTPC[pid + Npart]->Fill(track.pt(), valTPCnsigma, mult, flat);
      } else {
        hPtNsigmaTPC[pid + Npart]->Fill(track.pt(), valTPCnsigma);
      }
    }
    if (!track.hasTOF()) {
      return;
    }
    if (track.sign() > Cnull) {
      hPtNsigmaTOF[pid]->Fill(track.pt(), valTOFnsigma);
      hPtNsigmaTPCTOF[pid]->Fill(valTPCnsigma, valTOFnsigma);
    } else {
      hPtNsigmaTOF[pid + Npart]->Fill(track.pt(), valTOFnsigma);
      hPtNsigmaTPCTOF[pid + Npart]->Fill(valTPCnsigma, valTOFnsigma);
    }
  }

  template <FillType ft, bool fillHist = false, typename T>
  inline void fillTrackQA(T const& track)
  {
    if constexpr (fillHist) {
      flatchrg.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hPt"), track.pt());
      flatchrg.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hPhi"), track.phi());
      flatchrg.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hEta"), track.eta());
      flatchrg.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hDCAXYvsPt"), track.pt(), track.dcaXY());
      flatchrg.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hDCAZvsPt"), track.pt(), track.dcaZ());

      if (track.hasTPC() && track.hasITS()) {
        flatchrg.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hTPCCluster"), track.tpcNClsFindable() - track.tpcNClsFindableMinusFound());
        flatchrg.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hShTpcClvsPt"), track.pt(), track.tpcFractionSharedCls());
        flatchrg.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hNclTPCFoundvsPt"), track.pt(), track.tpcNClsFound());
        flatchrg.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hNClTPCPidvsPt"), track.pt(), track.tpcNClsPID());
        flatchrg.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hTPCnClsShared"), track.tpcNClsShared());
        flatchrg.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hTPCcrossedRows"), track.tpcNClsCrossedRows());
        flatchrg.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hTPCcrossedRowsOverFindableCls"), track.tpcCrossedRowsOverFindableCls());
        flatchrg.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hChi2ITSTrkSegment"), track.itsChi2NCl());
        flatchrg.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hITSnCls"), track.itsNCls());
      }

      flatchrg.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hTOFPvsBeta"), track.p(), track.beta());
      if (track.beta() > trkSelOpt.cfgTOFBetaPion && track.beta() < trkSelOpt.cfgTOFBetaPion + 0.05) { // TOF pions
        flatchrg.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hTOFpi"), track.eta(), track.p(), track.tpcSignal());
      }

      if (std::abs(track.eta()) < trkSelOpt.cfgTrkEtaMax) {
        if (isDCAxyWoCut(track)) {
          flatchrg.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hPtVsWOcutDCA"), track.pt(), track.dcaXY());
        }
      }
    }
    flatchrg.fill(HIST(Cprefix) + HIST(Cstatus[ft]) + HIST("hPVsPtEta"), track.p(), track.pt(), track.eta());
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
            flatchrg.fill(HIST(Cprefix) + HIST(CstatCalib[ft]) + HIST(Ccharge[chrg]) + HIST("hMIP"), mult, flat, track.eta(), dEdx);
          } else {
            flatchrg.fill(HIST(Cprefix) + HIST(CstatCalib[ft]) + HIST(Ccharge[chrg]) + HIST("hMIP"), track.eta(), dEdx);
          }
          flatchrg.fill(HIST(Cprefix) + HIST(CstatCalib[ft]) + HIST(Ccharge[chrg]) + HIST("hMIPVsPhi"), track.phi(), dEdx);
          flatchrg.fill(HIST(Cprefix) + HIST(CstatCalib[ft]) + HIST(Ccharge[chrg]) + HIST("hMIPVsPhiVsEta"), track.phi(), dEdx, track.eta());
          flatchrg.fill(HIST(Cprefix) + HIST(CstatCalib[ft]) + HIST(Ccharge[chrg]) + HIST("pMIPVsPhi"), track.phi(), dEdx);
        }
        if (dEdx > trkSelOpt.cfgDeDxMIPMax + 10. && dEdx < trkSelOpt.cfgDeDxMIPMax + 30.) { // Plateau electrons
          if (std::abs(track.beta() - 1) < trkSelOpt.cfgBetaPlateuMax) {
            if (cfgStoreThnSparse) {
              flatchrg.fill(HIST(Cprefix) + HIST(CstatCalib[ft]) + HIST(Ccharge[chrg]) + HIST("hPlateau"), mult, flat, track.eta(), dEdx);
            } else {
              flatchrg.fill(HIST(Cprefix) + HIST(CstatCalib[ft]) + HIST(Ccharge[chrg]) + HIST("hPlateau"), track.eta(), dEdx);
            }
            flatchrg.fill(HIST(Cprefix) + HIST(CstatCalib[ft]) + HIST(Ccharge[chrg]) + HIST("hPlateauVsPhi"), track.phi(), dEdx);
            flatchrg.fill(HIST(Cprefix) + HIST(CstatCalib[ft]) + HIST(Ccharge[chrg]) + HIST("hPlateauVsPhiVsEta"), track.phi(), dEdx, track.eta());
            flatchrg.fill(HIST(Cprefix) + HIST(CstatCalib[ft]) + HIST(Ccharge[chrg]) + HIST("pPlateauVsPhi"), track.phi(), dEdx);
          }
        }
      }
    }
  }

  template <bool fillHist = false, typename C>
  bool isGoodEvent(C const& collision)
  {
    if constexpr (fillHist) {
      flatchrg.fill(HIST("Events/hEvtSel"), evtSelAll);
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
      flatchrg.fill(HIST("Events/hEvtSel"), evtSelSel8);
    }
    if (evtSelOpt.cfgRemoveITSROFrameBorder && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    if constexpr (fillHist) {
      flatchrg.fill(HIST("Events/hEvtSel"), evtSelNoITSROFrameBorder);
    }
    if (evtSelOpt.cfgRemoveNoTimeFrameBorder && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    if constexpr (fillHist) {
      flatchrg.fill(HIST("Events/hEvtSel"), evtSelkNoTimeFrameBorder);
    }
    if (evtSelOpt.cfgRemoveNoSameBunchPileup && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if constexpr (fillHist) {
      flatchrg.fill(HIST("Events/hEvtSel"), evtSelkNoSameBunchPileup);
    }
    if (evtSelOpt.cfgRequireIsGoodZvtxFT0vsPV && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if constexpr (fillHist) {
      flatchrg.fill(HIST("Events/hEvtSel"), evtSelkIsGoodZvtxFT0vsPV);
    }
    if (evtSelOpt.cfgRequireIsVertexITSTPC && !collision.selection_bit(aod::evsel::kIsVertexITSTPC)) {
      return false;
    }
    if constexpr (fillHist) {
      flatchrg.fill(HIST("Events/hEvtSel"), evtSelkIsVertexITSTPC);
    }
    if (evtSelOpt.cfgRequirekIsVertexTOFmatched && !collision.selection_bit(aod::evsel::kIsVertexTOFmatched)) {
      return false;
    }
    if constexpr (fillHist) {
      flatchrg.fill(HIST("Events/hEvtSel"), evtSelkIsVertexTOFmatched);
    }
    if (std::abs(collision.posZ()) > evtSelOpt.cfgCutVtxZ) {
      return false;
    }
    if constexpr (fillHist) {
      flatchrg.fill(HIST("Events/hEvtSel"), evtSelVtxZ);
      flatchrg.fill(HIST("Events/hVtxZ"), collision.posZ());
    }
    if (evtSelOpt.cfgINELCut && !collision.isInelGt0()) {
      return false;
    }
    if constexpr (fillHist) {
      flatchrg.fill(HIST("Events/hEvtSel"), evtSelINELgt0);
    }
    if (rctCuts.requireRCTFlagChecker && !rctChecker(collision)) {
      return false;
    }
    if constexpr (fillHist) {
      flatchrg.fill(HIST("Events/hEvtSel"), evtSelRCTFlagChecker);
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
        return collision.centFT0M();
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

  float getMultMC(MCColls::iterator const& collision)
  {
    return getMult<MCColls::iterator, true>(collision);
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
            flatchrg.fill(HIST("FV0/hFV0amp"), chv0, amplCh);
            flatchrg.fill(HIST("FV0/pFV0amp"), chv0, amplCh);
            if (applyCalibGain) {
              flatchrg.fill(HIST("FV0/hFV0ampCorr"), chv0, amplCh / fv0AmplCorr[chv0]);
            }
          }
          if (amplCh > Cnull) {
            if (applyCalibGain) { // equalize gain channel-by-channel
              amplCh /= fv0AmplCorr[chv0];
            }
            if (chv0phi > Cnull) {
              fv0AmplitudeWoCalib[chv0phi] = amplCh;
              if constexpr (fillHist) {
                flatchrg.fill(HIST("FV0/hFV0AmplWCalib"), ich, fv0AmplitudeWoCalib[ich]);
              }
              if (chv0 < CinnerFV0) {
                rhoLatticeFV0[chv0phi] += amplCh;
              } else { // two channels per bin
                rhoLatticeFV0[chv0phi] += amplCh / 2.;
              }
              if constexpr (fillHist) {
                flatchrg.fill(HIST("FV0/hFV0AmplvsVtxzWoCalib"), collision.posZ(), rhoLatticeFV0[chv0phi]);
              }
              if (applyCalibVtx) {
                rhoLatticeFV0[chv0phi] *= zVtxMap->GetBinContent(zVtxMap->GetXaxis()->FindBin(chv0phi), zVtxMap->GetYaxis()->FindBin(collision.posZ()));
                if constexpr (fillHist) {
                  flatchrg.fill(HIST("FV0/hFV0AmplvsVtxzCalib"), collision.posZ(), rhoLatticeFV0[chv0phi]);
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
      if (cfgFillNsigmaQAHist) {
        fillNsigma<o2::track::PID::Pion>(tracksPerCollision, collision);
        fillNsigma<o2::track::PID::Kaon>(tracksPerCollision, collision);
        fillNsigma<o2::track::PID::Proton>(tracksPerCollision, collision);
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

    double etaMinFV0bins[CmaxRingsFV0] = {0.0};
    double etaMaxFV0bins[CmaxRingsFV0] = {0.0};
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
    float nCharged{0};
    for (const auto& mcPart : mcparts) {
      if (!isChrgParticle(mcPart.pdgCode())) {
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
          }
          isegment++;
        }
      }
      nCharged++;
    }

    vNch.push_back(nCharged);
    auto flatFV0 = calcFlatenicity(rhoLatticeFV0);

    if constexpr (fillHist) {
      flatchrg.fill(HIST("ResponseGen"), vNch[0], 1. - flatFV0);
      flatchrg.fill(HIST("h1flatencityFV0MCGen"), 1. - flatFV0);
    }
    vNch.clear();
    return 1. - flatFV0;
  }

  template <int pidSgn, o2::track::PID::ID id>
  void bookMcHist()
  {
    AxisSpec ptAxis{binOpt.axisPt, "#it{p}_{T} (GeV/#it{c})"};
    constexpr int ChistIdx = id + pidSgn * Npart;
    auto idx = static_cast<int>(id);
    const std::string strID = Form("/%s/%s", (pidSgn == Cnull && id < Npart) ? "pos" : "neg", Pid[idx]);
    hPtEffRec[ChistIdx] = flatchrg.add<TH1>("Tracks/hPtEffRec" + strID, " ; p_{T} (GeV/c)", kTH1F, {ptAxis});
    hPtEffGen[ChistIdx] = flatchrg.add<TH1>("Tracks/hPtEffGen" + strID, " ; p_{T} (GeV/c)", kTH1F, {ptAxis});
  }

  template <int pidSgn, o2::track::PID::ID id>
  void initEfficiency()
  {
    static_assert(pidSgn == Cnull || pidSgn == 1);
    static_assert(id > Cnull || id < Npart);
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
    static_assert(pidSgn == Cnull || pidSgn == 1);
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
    static_assert(pidSgn == Cnull || pidSgn == 1);
    constexpr int ChistIdx = id + pidSgn * Npart;
    const aod::McParticles::iterator& mcParticle = track.mcParticle();
    const CollsGen::iterator& collision = track.collision_as<CollsGen>();

    if (!isPID<pidSgn, id>(mcParticle)) {
      return;
    }
    flatchrg.fill(HIST("hPtOutNoEtaCut"), track.pt());
    if (std::abs(track.eta()) > trkSelOpt.cfgTrkEtaMax) {
      return;
    }
    if (std::abs(mcParticle.y()) > trkSelOpt.cfgRapMax) {
      return;
    }

    if ((collision.has_mcCollision() && (mcParticle.mcCollisionId() != collision.mcCollisionId())) || !collision.has_mcCollision()) {
      if (!mcParticle.isPhysicalPrimary()) {
        if (mcParticle.getProcess() == CprocessIdWeak) {
          hDCAxyBadCollWeak[ChistIdx]->Fill(track.pt(), track.dcaXY());
        } else {
          hDCAxyBadCollMat[ChistIdx]->Fill(track.pt(), track.dcaXY());
        }
      } else {
        hDCAxyBadCollPrim[ChistIdx]->Fill(track.pt(), track.dcaXY());
      }
    }

    if (!isDCAxyCut(track)) {
      return;
    }
    flatchrg.fill(HIST("hPtVsDCAxyAll"), track.pt(), track.dcaXY());

    if (selTPCtrack(track)) {
      hPtEffRec[ChistIdx]->Fill(mcParticle.pt());
    }

    if (!mcParticle.isPhysicalPrimary()) {
      if (mcParticle.getProcess() == CprocessIdWeak) {
        hPtEffRecWeak[ChistIdx]->Fill(mult, flat, track.pt());
        hPtVsDCAxyWeak[ChistIdx]->Fill(track.pt(), track.dcaXY());
        flatchrg.fill(HIST("hPtVsDCAxyWeakAll"), track.pt(), track.dcaXY());
      } else {
        hPtEffRecMat[ChistIdx]->Fill(mult, flat, track.pt());
        hPtVsDCAxyMat[ChistIdx]->Fill(track.pt(), track.dcaXY());
        flatchrg.fill(HIST("hPtVsDCAxyMatAll"), track.pt(), track.dcaXY());
      }
    } else {
      hPtEffRecPrim[ChistIdx]->Fill(mult, flat, track.pt());
      hPtVsDCAxyPrim[ChistIdx]->Fill(track.pt(), track.dcaXY());
      flatchrg.fill(HIST("hPtVsDCAxyPrimAll"), track.pt(), track.dcaXY());
    }
  }

  template <int pidSgn, o2::track::PID::ID id, bool recoEvt = false>
  void fillMCGen(aod::McParticles::iterator const& mcParticle, const float mult, const float flat)
  {
    static_assert(pidSgn == Cnull || pidSgn == 1);
    constexpr int ChistIdx = id + pidSgn * Npart;

    if (!isPID<pidSgn, id>(mcParticle)) {
      return;
    }

    if constexpr (recoEvt) {
      hPtGenRecEvt[ChistIdx]->Fill(mcParticle.pt());
      if (mcParticle.isPhysicalPrimary()) {
        hPtGenPrimRecEvt[ChistIdx]->Fill(mcParticle.pt());
      }
      return;
    }

    if (!mcParticle.isPhysicalPrimary()) {
      if (mcParticle.getProcess() == CprocessIdWeak) {
        hPtEffGenWeak[ChistIdx]->Fill(mult, flat, mcParticle.pt());
      } else {
        hPtEffGenMat[ChistIdx]->Fill(mult, flat, mcParticle.pt());
      }
    } else {
      hPtEffGenPrim[ChistIdx]->Fill(mult, flat, mcParticle.pt());
      hPtEffGen[ChistIdx]->Fill(mcParticle.pt());
    }
  }

  void processSgnLoss(MCColls::iterator const& mcCollision,
                      CollsGenSgn const& collisions,
                      aod::FV0As const& /*fv0s*/,
                      aod::McParticles const& particles)
  {
    float flat;
    float mult;
    if (flatSelOpt.useFlatData) {
      float flatRec = 999.0;
      float multRec = 999.0;
      for (const auto& collision : collisions) {
        multRec = getMult(collision);
        flatRec = fillFlat<false>(collision);
      }
      flat = flatRec;
      mult = multRec;
      flatchrg.fill(HIST("hFlatMCGenRecColl"), flatRec);
    } else {
      float flatGen = fillFlatMC<false>(particles);
      flat = flatGen;
      flatchrg.fill(HIST("hFlatMCGen"), flatGen);
      float multGen = getMultMC(mcCollision);
      mult = multGen;
    }

    // Evt loss den
    flatchrg.fill(HIST("hEvtMcGen"), 0.5);
    if (std::abs(mcCollision.posZ()) > evtSelOpt.cfgCutVtxZ) {
      return;
    }
    flatchrg.fill(HIST("hEvtMcGen"), 1.5);

    bool isINELgt0mc = false;
    if (pwglf::isINELgtNmc(particles, 0, pdg)) {
      isINELgt0mc = true;
      flatchrg.fill(HIST("hEvtMcGen"), 2.5);
      flatchrg.fill(HIST("hFlatGenINELgt0"), flat);
    }

    // Sgn loss den
    for (const auto& particle : particles) {
      if (!particle.isPhysicalPrimary()) {
        continue;
      }
      if (std::abs(particle.y()) > trkSelOpt.cfgRapMax) {
        continue;
      }
      static_for<0, 5>([&](auto i) {
        constexpr int Cidx = i.value;
        if (particle.pdgCode() == PidSgn[Cidx]) {
          flatchrg.fill(HIST(Cprefix) + HIST(Cspecies[Cidx]) + HIST(CpTgenPrimSgn), mult, flat, particle.pt());
          if (isINELgt0mc) {
            flatchrg.fill(HIST(Cprefix) + HIST(Cspecies[Cidx]) + HIST(CpTgenPrimSgnINEL), mult, flat, particle.pt());
          }
        }
      });
    }

    int nRecCollINEL = 0;
    int nRecCollINELgt0 = 0;
    for (const auto& collision : collisions) {
      // Evt split num
      flatchrg.fill(HIST("hEvtMCRec"), 0.5);
      if (!isGoodEvent<false>(collision)) {
        continue;
      }
      flatchrg.fill(HIST("hEvtMCRec"), 1.5);

      nRecCollINEL++;

      if (collision.isInelGt0() && isINELgt0mc) {
        flatchrg.fill(HIST("hEvtMCRec"), 2.5);
        nRecCollINELgt0++;
      }
      // Sgn split num
      for (const auto& particle : particles) {
        if (!particle.isPhysicalPrimary()) {
          continue;
        }
        if (std::abs(particle.y()) > trkSelOpt.cfgRapMax) {
          continue;
        }
        static_for<0, 5>([&](auto i) {
          constexpr int Cidx = i.value;
          if (particle.pdgCode() == PidSgn[Cidx]) {
            flatchrg.fill(HIST(Cprefix) + HIST(Cspecies[Cidx]) + HIST(CpTrecCollPrimSgn), mult, flat, particle.pt());
            if (nRecCollINELgt0) {
              flatchrg.fill(HIST(Cprefix) + HIST(Cspecies[Cidx]) + HIST(CpTrecCollPrimSgnINEL), mult, flat, particle.pt());
            }
          }
        });
      }
    }

    if (nRecCollINEL < 1) {
      return;
    }
    // Evt loss num
    flatchrg.fill(HIST("hEvtMcGenRecColl"), 0.5);
    if (nRecCollINELgt0 > Cnull) {
      flatchrg.fill(HIST("hEvtMcGenRecColl"), 1.5);
    }

    // Sgn loss num
    for (const auto& particle : particles) {
      if (!particle.isPhysicalPrimary()) {
        continue;
      }
      if (std::abs(particle.y()) > trkSelOpt.cfgRapMax) {
        continue;
      }
      static_for<0, 5>([&](auto i) {
        constexpr int Cidx = i.value;
        if (particle.pdgCode() == PidSgn[Cidx]) {
          flatchrg.fill(HIST(Cprefix) + HIST(Cspecies[Cidx]) + HIST(CpTGenRecCollPrimSgn), mult, flat, particle.pt());
          if (nRecCollINELgt0) {
            flatchrg.fill(HIST(Cprefix) + HIST(Cspecies[Cidx]) + HIST(CpTGenRecCollPrimSgnINEL), mult, flat, particle.pt());
          }
        }
      });
    }
  }
  PROCESS_SWITCH(FlattenictyPikp, processSgnLoss, "process to calcuate signal/event lossses", false);

  // using Particles = soa::Filtered<aod::McParticles>;
  // expressions::Filter primaries = (aod::mcparticle::flags & (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary) == (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary;
  void processMCclosure(Colls::iterator const& collision,
                        MyPIDTracks const& tracks,
                        MyLabeledTracks const& mcTrackLabels,
                        aod::McParticles const& particles,
                        // Particles const& particles,
                        aod::FV0As const& /*fv0s*/,
                        aod::BCsWithTimestamps const& /*bcs*/)
  {
    const float multRec = getMult(collision);
    const float flatRec = fillFlat<false>(collision);
    for (const auto& track : tracks) {
      if (!track.has_collision()) {
        continue;
      }
      const auto& coll = track.collision_as<Colls>();
      if (trkSelOpt.cfgRejectTrkAtTPCSector) {
        auto bc = coll.template bc_as<aod::BCsWithTimestamps>();
        int currentRun = bc.runNumber();
        if (runNumber != currentRun) {
          initCCDB(bc);
          runNumber = currentRun;
        }
      }
      if (!isGoodEvent<false>(coll)) {
        continue;
      }
      if (!isGoodTrack(track, magField)) {
        continue;
      }
      if (!isDCAxyCut(track)) {
        continue;
      }
      const auto& mcLabel = mcTrackLabels.iteratorAt(track.globalIndex());
      const auto& mcParticle = particles.iteratorAt(mcLabel.mcParticleId());

      static_for<0, 4>([&](auto i) {
        constexpr int Cidx = i.value;
        if ((std::abs(o2::aod::pidutils::tpcNSigma<Cidx>(track)) < trkSelOpt.cfgNsigmaMax) && std::abs(track.rapidity(o2::track::PID::getMass(Cidx))) <= trkSelOpt.cfgRapMax) {
          if (std::fabs(mcParticle.pdgCode()) == PDGs[Cidx]) {
            flatchrg.fill(HIST(Cprefix) + HIST(CspeciesAll[Cidx]) + HIST(CpTmcClosurePrim), multRec, flatRec, track.pt());
          }
        }
      });
    }
  }
  PROCESS_SWITCH(FlattenictyPikp, processMCclosure, "process MC closure test", false);

  Preslice<MyLabeledPIDTracks> perCollTrk = aod::track::collisionId;
  PresliceUnsorted<CollsGen> perCollMcLabel = aod::mccollisionlabel::mcCollisionId;
  Preslice<aod::McParticles> perCollMcPart = aod::mcparticle::mcCollisionId;

  void processMC(MCColls const& mcCollisions,
                 CollsGen const& collisions,
                 MyLabeledPIDTracks const& tracks,
                 aod::McParticles const& mcparticles)
  {
    flatchrg.fill(HIST("hEvtGenRec"), 1.f, mcCollisions.size());
    flatchrg.fill(HIST("hEvtGenRec"), 2.f, collisions.size());

    for (const auto& mcCollision : mcCollisions) {
      if (mcCollision.isInelGt0()) {
        flatchrg.fill(HIST("hEvtGenRec"), 3.f);
      }
      flatchrg.fill(HIST("hEvtMcGenColls"), 1);
      const auto groupedColls = collisions.sliceBy(perCollMcLabel, mcCollision.globalIndex());
      const auto groupedParts = mcparticles.sliceBy(perCollMcPart, mcCollision.globalIndex());
      const float flatMC = fillFlatMC<true>(groupedParts);
      const float multMC = getMultMC(mcCollision);
      if (groupedColls.size() < 1) { // if MC events have no rec collisions
        continue;
      }
      flatchrg.fill(HIST("hEvtMcGenColls"), 2);
      for (const auto& collision : groupedColls) {
        flatchrg.fill(HIST("hEvtMcGenColls"), 3);
        if (!isGoodEvent<false>(collision)) {
          continue;
        }
        flatchrg.fill(HIST("hEvtMcGenColls"), 4);
        const auto groupedTrks = tracks.sliceBy(perCollTrk, collision.globalIndex());
        for (const auto& track : groupedTrks) {
          if (!isDCAxyWoCut(track)) {
            continue;
          }
          if (!track.has_mcParticle()) {
            flatchrg.fill(HIST("PtOutFakes"), track.pt());
            continue;
          }
          const auto& mcParticle = track.mcParticle();
          if (std::abs(mcParticle.y()) > trkSelOpt.cfgRapMax) {
            continue;
          }
          if (!track.has_collision()) {
            continue;
          }
          static_for<0, 1>([&](auto pidSgn) {
            fillMCRecTrack<pidSgn, o2::track::PID::Pion>(track, multMC, flatMC);
            fillMCRecTrack<pidSgn, o2::track::PID::Kaon>(track, multMC, flatMC);
            fillMCRecTrack<pidSgn, o2::track::PID::Proton>(track, multMC, flatMC);
          });
        }

        if (std::abs(mcCollision.posZ()) > evtSelOpt.cfgCutVtxZ) {
          continue;
        }

        flatchrg.fill(HIST("Events/hVtxZRec"), collision.posZ());
        flatchrg.fill(HIST("Events/hVtxZGen"), mcCollision.posZ());

        if (evtSelOpt.cfgINELCut.value) {
          if (!o2::pwglf::isINELgt0mc(groupedParts, pdg)) {
            continue;
          }
        }

        for (const auto& particle : groupedParts) {
          if (std::abs(particle.y()) > trkSelOpt.cfgRapMax) {
            continue;
          }
          static_for<0, 1>([&](auto pidSgn) {
            fillMCGen<pidSgn, o2::track::PID::Pion, true>(particle, multMC, flatMC);
            fillMCGen<pidSgn, o2::track::PID::Kaon, true>(particle, multMC, flatMC);
            fillMCGen<pidSgn, o2::track::PID::Proton, true>(particle, multMC, flatMC);
          });
        }
      } // reco collisions

      if (evtSelOpt.cfgINELCut.value) {
        if (!o2::pwglf::isINELgt0mc(groupedParts, pdg)) {
          continue;
        }
      }

      for (const auto& mcParticle : groupedParts) {
        if (std::abs(mcParticle.y()) > trkSelOpt.cfgRapMax) {
          continue;
        }
        static_for<0, 1>([&](auto pidSgn) {
          fillMCGen<pidSgn, o2::track::PID::Pion>(mcParticle, multMC, flatMC);
          fillMCGen<pidSgn, o2::track::PID::Kaon>(mcParticle, multMC, flatMC);
          fillMCGen<pidSgn, o2::track::PID::Proton>(mcParticle, multMC, flatMC);
        });
      }

      static_for<0, 1>([&](auto pidSgn) {
        fillEfficiency<pidSgn, o2::track::PID::Pion>();
        fillEfficiency<pidSgn, o2::track::PID::Kaon>();
        fillEfficiency<pidSgn, o2::track::PID::Proton>();
      });

    } // gen collisions
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
