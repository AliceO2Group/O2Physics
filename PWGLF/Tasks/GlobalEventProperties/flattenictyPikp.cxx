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
#include "Common/DataModel/PIDResponse.h"
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

#include <cmath>
#include <cstdlib>
#include <map>
#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::constants::math;

auto static constexpr kMinCharge = 3.f;
// FV0 specific constants
static constexpr int kMaxRingsFV0 = 5;
static constexpr int kNCellsFV0 = 48;
static constexpr int kInnerFV0 = 32;
static constexpr float kFV0IndexPhi[5] = {0., 8., 16., 24., 32.};
static constexpr float kEtaMaxFV0 = 5.1;
static constexpr float kEtaMinFV0 = 2.2;
static constexpr float kDEtaFV0 = (kEtaMaxFV0 - kEtaMinFV0) / kMaxRingsFV0;
// PID names
static constexpr int kProcessIdWeak = 4;
static constexpr int Ncharges = 2;
static constexpr o2::track::PID::ID Npart = 5;
static constexpr o2::track::PID::ID NpartChrg = Npart * Ncharges;
static constexpr int PDGs[] = {11, 13, 211, 321, 2212};
static constexpr int PidSgn[NpartChrg] = {11, 13, 211, 321, 2212, -11, -13, -211, -321, -2212};
static constexpr const char* Pid[Npart] = {"el", "mu", "pi", "ka", "pr"};
static constexpr const char* PidChrg[NpartChrg] = {"e^{-}", "#mu^{-}", "#pi^{+}", "K^{+}", "p", "e^{+}", "#mu^{+}", "#pi^{-}", "K^{-}", "#bar{p}"};
static constexpr std::string_view kSpecies[NpartChrg] = {"Elminus", "Muplus", "PiPlus", "KaPlus", "Pr", "ElPlus", "MuMinus", "PiMinus", "KaMinus", "PrBar"};
static constexpr std::string_view kSpeciesAll[Npart] = {"El", "Mu", "Pi", "Ka", "Pr"};
// histogram naming
static constexpr std::string_view kCharge[] = {"all/", "pos/", "neg/"};
static constexpr std::string_view kPrefix = "Tracks/";
static constexpr std::string_view kPrefixCleanTof = "Tracks/CleanTof/";
static constexpr std::string_view kPrefixCleanV0 = "Tracks/CleanV0/";
static constexpr std::string_view kStatus[] = {"preSel/", "postSel/"};
static constexpr std::string_view kStatCalib[] = {"preCalib/", "postCalib/"};
static constexpr std::string_view kPtGenPrimSgn = "/hPtGenPrimSgn";
static constexpr std::string_view kPtGenPrimSgnF = "Tracks/{}/hPtGenPrimSgn";
static constexpr std::string_view kPtGenPrimSgnINEL = "/hPtGenPrimSgnINEL";
static constexpr std::string_view kPtGenPrimSgnINELF = "Tracks/{}/hPtGenPrimSgnINEL";
static constexpr std::string_view kPtRecCollPrimSgn = "/hPtRecCollPrimSgn";
static constexpr std::string_view kPtRecCollPrimSgnF = "Tracks/{}/hPtRecCollPrimSgn";
static constexpr std::string_view kPtRecCollPrimSgnINEL = "/hPtRecCollPrimSgnINEL";
static constexpr std::string_view kPtRecCollPrimSgnINELF = "Tracks/{}/hPtRecCollPrimSgnINEL";
static constexpr std::string_view kPtGenRecCollPrimSgn = "/hPtGenRecCollPrimSgn";
static constexpr std::string_view kPtGenRecCollPrimSgnF = "Tracks/{}/hPtGenRecCollPrimSgn";
static constexpr std::string_view kPtGenRecCollPrimSgnINEL = "/hPtGenRecCollPrimSgnINEL";
static constexpr std::string_view kPtGenRecCollPrimSgnINELF = "Tracks/{}/hPtGenRecCollPrimSgnINEL";
static constexpr std::string_view kPtMCclosurePrim = "/hPtMCclosurePrim";
static constexpr std::string_view kPtMCclosurePrimF = "Tracks/{}/hPtMCclosurePrim";

enum V0sSel {
  kNaN = -1,
  kGa = 0,
  kKz = 1,
  kLam = 2,
  kaLam = 3
};

enum TrkSelNoFilt {
  trkSelEta,
  trkSelPt,
  trkSelDCA,
  trkNRowsTPC,
  trkSelNCls,
  trkSelTPCBndr,
  nTrkSel
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
  nEvtSel
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

struct MultE {
  static constexpr int kNoMult = 0;
  static constexpr int kMultFT0M = 1;
  static constexpr int kMultTPC = 2;
};

std::array<float, kNCellsFV0> rhoLatticeFV0{0};
std::array<float, kNCellsFV0> fv0AmplitudeWoCalib{0};

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
  int runNumber{-1};

  Configurable<int> multEst{"multEst", 1, "0: without multiplicity; 1: MultFT0M; 2: MultTPC"};
  Configurable<bool> applyCalibGain{"applyCalibGain", false, "equalize detector amplitudes"};
  Configurable<bool> applyCalibVtx{"applyCalibVtx", false, "equalize Amp vs vtx"};
  Configurable<bool> applyCalibDeDx{"applyCalibDeDx", false, "calibration of dedx signal"};
  Configurable<bool> cfgFillTrackQaHist{"cfgFillTrackQaHist", false, "fill track QA histograms"};
  Configurable<bool> cfgFilldEdxCalibHist{"cfgFilldEdxCalibHist", false, "fill dEdx calibration histograms"};
  Configurable<bool> cfgFilldEdxQaHist{"cfgFilldEdxQaHist", false, "fill dEdx QA histograms"};
  Configurable<bool> cfgFillNsigmaQAHist{"cfgFillNsigmaQAHist", false, "fill nsigma QA histograms"};
  Configurable<bool> cfgFillV0Hist{"cfgFillV0Hist", false, "fill V0 histograms"};
  Configurable<bool> cfgFillChrgType{"cfgFillChrgType", false, "fill histograms per charge types"};
  Configurable<std::vector<float>> paramsFuncMIPpos{"paramsFuncMIPpos", std::vector<float>{-1.f}, "parameters of pol2"};
  Configurable<std::vector<float>> paramsFuncMIPneg{"paramsFuncMIPneg", std::vector<float>{-1.f}, "parameters of pol2"};
  Configurable<std::vector<float>> paramsFuncMIPall{"paramsFuncMIPall", std::vector<float>{-1.f}, "parameters of pol2"};
  Configurable<std::string> cfgGainEqCcdbPath{"cfgGainEqCcdbPath", "Users/g/gbencedi/flattenicity/GainEq", "CCDB path for gain equalization constants"};
  Configurable<std::string> cfgVtxEqCcdbPath{"cfgVtxEqCcdbPath", "Users/g/gbencedi/flattenicity/ZvtxEq", "CCDB path for z-vertex equalization constants"};
  Configurable<bool> cfgUseCcdbForRun{"cfgUseCcdbForRun", true, "Get ccdb object based on run number instead of timestamp"};
  Configurable<bool> cfgStoreThnSparse{"cfgStoreThnSparse", false, "Store histograms as THnSparse"};

  struct : ConfigurableGroup {
    Configurable<bool> cfgCustomTVX{"cfgCustomTVX", false, "Ask for custom TVX instead of sel8"};
    Configurable<bool> cfgRemoveNoTimeFrameBorder{"cfgRemoveNoTimeFrameBorder", true, "Bunch crossing is far from Time Frame borders"};
    Configurable<bool> cfgRemoveITSROFrameBorder{"cfgRemoveITSROFrameBorder", true, "Bunch crossing is far from ITS RO Frame border"};
    Configurable<float> cfgCutVtxZ{"cfgCutVtxZ", 10.0f, "Accepted z-vertex range"};
    Configurable<bool> cfgINELCut{"cfgINELCut", true, "INEL event selection"};
    Configurable<bool> cfgRemoveNoSameBunchPileup{"cfgRemoveNoSameBunchPileup", false, "Reject collisions in case of pileup with another collision in the same foundBC"};
    Configurable<bool> cfgRequireIsGoodZvtxFT0vsPV{"cfgRequireIsGoodZvtxFT0vsPV", false, "Small difference between z-vertex from PV and from FT0"};
    Configurable<bool> cfgRequireIsVertexITSTPC{"cfgRequireIsVertexITSTPC", false, "At least one ITS-TPC track (reject vertices built from ITS-only tracks)"};
    Configurable<bool> cfgRequirekIsVertexTOFmatched{"cfgRequirekIsVertexTOFmatched", false, "Require kIsVertexTOFmatched: at least one of vertex contributors is matched to TOF"};
  } evtSelOpt;

  struct : ConfigurableGroup {
    Configurable<bool> useFlatData{"useFlatData", true, "use flattenicity from rec collisions"};
  } flatSelOpt;

  struct : ConfigurableGroup {
    ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0}, "pT binning"};
    ConfigurableAxis axisFlatPerc{"axisFlatPerc", {102, -0.01, 1.01}, "Flattenicity percentiles binning"};
    ConfigurableAxis axisMultPerc{"axisMultPerc", {20, 0, 100}, "Multiplicity percentiles binning"};
    ConfigurableAxis axisVertexZ{"axisVertexZ", {80, -20., 20.}, "Vertex z binning"};
    ConfigurableAxis axisMult{"axisMult", {301, -0.5, 300.5}, "Multiplicity binning"};
    ConfigurableAxis axisDCAxy{"axisDCAxy", {200, -5, 5}, "DCAxy binning"};
    ConfigurableAxis axisDCAz{"axisDCAz", {200, -5, 5}, "DCAz binning"};
    ConfigurableAxis axisPhi = {"axisPhi", {60, 0, constants::math::TwoPI}, "#varphi binning"};
    ConfigurableAxis axisPhiMod = {"axisPhiMod", {100, 0, constants::math::PI / 9}, "fmod(#varphi,#pi/9)"};
    ConfigurableAxis axisEta = {"axisEta", {8, -0.8, 0.8}, "#eta binning"};
    ConfigurableAxis axisDedx{"axisDedx", {100, 0, 100}, "dE/dx binning"};
    ConfigurableAxis axisNsigmaTPC{"axisNsigmaTPC", {200, -10, 10}, "nsigmaTPC binning"};
    ConfigurableAxis axisNsigmaTOF{"axisNsigmaTOF", {200, -10, 10}, "nsigmaTOF binning"};
    ConfigurableAxis axisAmplFV0{"axsAmplFV0", {4096, 0, 4096}, "FV0 amplitude (ADC) binning"};
    ConfigurableAxis axisAmplFV0Sum{"axisAmplFV0Sum", {4096, 0, 4096 * 49}, "FV0 amplitude sum (ADC) binning"};
    ConfigurableAxis axisChannelFV0{"axisChannelFV0", {49, 0., 49.}, "FV0 channel ID binning"};
  } binOpt;

  struct : ConfigurableGroup {
    Configurable<float> cfgTrkEtaMax{"cfgTrkEtaMax", 0.8f, "Eta range for tracks"};
    Configurable<float> cfgRapMax{"cfgRapMax", 0.5f, "Maximum range of rapidity for tracks"};
    Configurable<float> cfgTrkPtMin{"cfgTrkPtMin", 0.15f, "Minimum pT of tracks"};
    Configurable<float> cfgNclTPCMin{"cfgNclTPCMin", 100.0f, "Minimum of number of TPC clusters"};
    Configurable<float> cfgPhiCutPtMin{"cfgPhiCutPtMin", 2.0f, "Minimum pT for phi cut"};
    Configurable<float> cfgTOFBetaPion{"cfgTOFBetaPion", 1.0f, "Minimum beta for TOF pions"};
    Configurable<bool> cfgUseExtraTrkCut{"cfgUseExtraTrkCut", true, "Use extra track cut"};
    Configurable<std::string> cfgGeoTrkCutMin{"cfgGeoTrkCutMin", "0.06/x+pi/18.0-0.06", "ROOT TF1 formula for minimum phi cut in TPC"};
    Configurable<std::string> cfgGeoTrkCutMax{"cfgGeoTrkCutMax", "0.1/x+pi/18.0+0.06", "ROOT TF1 formula for maximum phi cut in TPC"};
    Configurable<float> cfgMomMIPMax{"cfgMomMIPMax", 0.6f, "Maximum momentum of MIP pions"};
    Configurable<float> cfgMomMIPMin{"cfgMomMIPMin", 0.4f, "Minimum momentum of MIP pions"};
    Configurable<float> cfgDeDxMIPMax{"cfgDeDxMIPMax", 60.0f, "Maximum range of MIP dedx"};
    Configurable<float> cfgDeDxMIPMin{"cfgDeDxMIPMin", 40.0f, "Maximum range of MIP dedx"};
    Configurable<float> cfgNsigmaMax{"cfgNsigmaMax", 100.0f, "Maximum range of nsgima for tracks"};
    Configurable<float> cfgMomSelPiTOF{"cfgMomSelPiTOF", 0.7f, "Momentum cut for TOF pions"};
    Configurable<float> cfgNsigSelKaTOF{"cfgNsigSelKaTOF", 3.0f, "Nsigma cut for TOF kaons"};
    Configurable<float> cfgBetaPlateuMax{"cfgBetaPlateuMax", 0.1f, "Beta max for Plateau electrons"};
  } trkSelOpt;

  struct : ConfigurableGroup {
    Configurable<float> cfgDCAv0daughter{"cfgDCAv0daughter", 0.3, "DCA of V0 daughter tracks"};
    Configurable<float> cfgV0etamax{"cfgV0etamax", 0.9f, "max eta of V0s"};
    Configurable<float> cfgV0DaughterTpcMomMax{"cfgV0DaughterTpcMomMax", 0.6f, "Maximum momentum of V0 daughter tracks in TPC"};
    Configurable<float> cfgTPCnClsmin{"cfgTPCnClsmin", 50.0f, "cfgTPCnClsmin"};
    Configurable<float> cfgDCAposToPV{"cfgDCAposToPV", 0.05f, "cfgDCAposToPV"};
    Configurable<float> cfgv0cospa{"cfgv0cospa", 0.998, "V0 CosPA"};
    Configurable<float> cfgv0Rmin{"cfgv0Rmin", 0.0, "cfgv0Rmin"};
    Configurable<float> cfgv0Rmax{"cfgv0Rmax", 90.0, "cfgv0Rmax"};
    Configurable<float> cfgdmassG{"cfgdmassG", 0.1, "max mass for photon"};
    Configurable<float> cfgdmassK{"cfgdmassK", 0.1, "max mass for K0s"};
    Configurable<float> cfgdmassL{"cfgdmassL", 0.1, "max mass for Lambda"};
    Configurable<float> cfgdmassAL{"cfgdmassAL", 0.1, "max mass for anti Lambda"};
    Configurable<float> cfgNsigmaElTPC{"cfgNsigmaElTPC", 5.0, "max nsigma of TPC for electorn"};
    Configurable<float> cfgNsigmaPiTPC{"cfgNsigmaPiTPC", 5.0, "max nsigma of TPC for pion"};
    Configurable<float> cfgNsigmaPrTPC{"cfgNsigmaPrTPC", 5.0, "max nsigma of TPC for proton"};
    Configurable<float> cfgNsigmaElTOF{"cfgNsigmaElTOF", 3.0, "max nsigma of TOF for electorn"};
    Configurable<float> cfgNsigmaPiTOF{"cfgNsigmaPiTOF", 3.0, "max nsigma of TOF for pion"};
    Configurable<float> cfgNsigmaPrTOF{"cfgNsigmaPrTOF", 3.0, "max nsigma of TOF for proton"};
  } v0SelOpt;

  Service<ccdb::BasicCCDBManager> ccdb;
  struct : ConfigurableGroup {
    Configurable<float> cfgMagField{"cfgMagField", 99999, "Configurable magnetic field;default CCDB will be queried"};
    Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
    Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  } ccdbConf;

  TrackSelection mTrackSelector;
  Configurable<bool> isCustomTracks{"isCustomTracks", false, "Use custom track cuts"};
  Configurable<float> requirePt{"requirePt", 0.15f, "Set minimum pT of tracks"};
  Configurable<float> requireEta{"requireEta", 0.8f, "Set eta range of tracks"};
  Configurable<int> setITSreq{"setITSreq", 1, "0 = Run3ITSibAny, 1 = Run3ITSallAny, 2 = Run3ITSall7Layers, 3 = Run3ITSibTwo"};
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
  Configurable<float> minTPCNClsFound{"minTPCNClsFound", 70.0f, "Additional cut on the minimum value of the number of found clusters in the TPC"};

  TF1* fPhiCutLow = nullptr;
  TF1* fPhiCutHigh = nullptr;

  std::vector<std::unique_ptr<TF1>> fDeDxVsEta;
  std::vector<std::vector<float>> vecParams;

  void init(InitContext&)
  {
    auto vecParamsMIPpos = (std::vector<float>)paramsFuncMIPpos;
    auto vecParamsMIPneg = (std::vector<float>)paramsFuncMIPneg;
    auto vecParamsMIPall = (std::vector<float>)paramsFuncMIPall;

    auto addVec = [&](std::vector<std::vector<float>>& targetVec, const std::string& name) {
      targetVec.emplace_back(vecParamsMIPpos);
      targetVec.emplace_back(vecParamsMIPneg);
      targetVec.emplace_back(vecParamsMIPall);
      if (!vecParams.size()) {
        LOG(info) << "size of " << name << "is zero.";
      }
    };
    addVec(vecParams, "vecParams");
    for (const auto& params : vecParams) {
      fDeDxVsEta.emplace_back(setFuncPars(params));
    }

    ccdb->setURL(ccdbConf.ccdbUrl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    ccdb->setFatalWhenNull(false);

    if (isCustomTracks.value) {
      mTrackSelector = getGlobalTrackSelectionRun3ITSMatch(setITSreq.value);
      mTrackSelector.SetPtRange(requirePt.value, 1e10f);
      mTrackSelector.SetEtaRange(-requireEta.value, requireEta.value);
      mTrackSelector.SetRequireITSRefit(requireITS.value);
      mTrackSelector.SetRequireTPCRefit(requireTPC.value);
      mTrackSelector.SetRequireGoldenChi2(requireGoldenChi2.value);
      mTrackSelector.SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC.value);
      mTrackSelector.SetMaxChi2PerClusterITS(maxChi2PerClusterITS.value);
      mTrackSelector.SetMinNClustersITS(minITSnClusters.value);
      mTrackSelector.SetMinNCrossedRowsTPC(minNCrossedRowsTPC.value);
      mTrackSelector.SetMinNClustersTPC(minTPCNClsFound.value);
      mTrackSelector.SetMinNCrossedRowsOverFindableClustersTPC(minNCrossedRowsOverFindableClustersTPC.value);
      // //     mTrackSelector.SetMaxDcaXYPtDep([](float pt) { return 0.0105f + 0.0350f / pow(pt, 1.1f); });
      mTrackSelector.SetMaxDcaXYPtDep([](float /*pt*/) { return 10000.f; });
      mTrackSelector.SetMaxDcaZ(maxDcaZ.value);
      mTrackSelector.print();
    }

    const AxisSpec chargeAxis{2, -2.f, 2.f, "Charge"};
    const AxisSpec dEdxAxis{binOpt.axisDedx, "TPC dEdx (a.u.)"};
    const AxisSpec vtxzAxis{binOpt.axisVertexZ, "Z_{vtx} (cm)"};
    const AxisSpec flatAxis{binOpt.axisFlatPerc, "Flat FV0"};
    const AxisSpec etaAxis{binOpt.axisEta, "#eta"};
    const AxisSpec phiAxis{binOpt.axisPhi, "#varphi"};
    const AxisSpec phiAxisMod{binOpt.axisPhiMod, "fmod(#varphi,#pi/9)"};
    const AxisSpec pAxis{binOpt.axisPt, "#it{p} (GeV/#it{c})"};
    const AxisSpec ptAxis{binOpt.axisPt, "#it{p}_{T} (GeV/#it{c})"};
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
      case MultE::kNoMult:
        break;
      case MultE::kMultFT0M:
        multAxis.name = "multFT0M";
        break;
      case MultE::kMultTPC:
        multAxis.name = "multTPC";
        break;
      default:
        LOG(fatal) << "No valid option for mult estimator " << multEst;
    }

    // Event counter
    flatchrg.add("Events/hEvtSel", "Number of events; Cut; #Events Passed Cut", {HistType::kTH1F, {{nEvtSel, 0, nEvtSel}}});
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
    // Track counter
    flatchrg.add("Tracks/hTrkSel", "Number of tracks; Cut; #Tracks Passed Cut", {HistType::kTH1F, {{nTrkSel, 0, nTrkSel}}});
    flatchrg.get<TH1>(HIST("Tracks/hTrkSel"))->GetXaxis()->SetBinLabel(trkSelEta + 1, "Eta");
    flatchrg.get<TH1>(HIST("Tracks/hTrkSel"))->GetXaxis()->SetBinLabel(trkSelPt + 1, "Pt");
    flatchrg.get<TH1>(HIST("Tracks/hTrkSel"))->GetXaxis()->SetBinLabel(trkSelDCA + 1, "DCA");
    flatchrg.get<TH1>(HIST("Tracks/hTrkSel"))->GetXaxis()->SetBinLabel(trkNRowsTPC + 1, "trkNRowsTPC");
    flatchrg.get<TH1>(HIST("Tracks/hTrkSel"))->GetXaxis()->SetBinLabel(trkSelNCls + 1, "NClsTPC");
    flatchrg.get<TH1>(HIST("Tracks/hTrkSel"))->GetXaxis()->SetBinLabel(trkSelTPCBndr + 1, "TPC Boundary");

    if (trkSelOpt.cfgUseExtraTrkCut) {
      fPhiCutLow = new TF1("fPhiCutLow", trkSelOpt.cfgGeoTrkCutMin.value.c_str(), 0, 100);
      fPhiCutHigh = new TF1("fPhiCutHigh", trkSelOpt.cfgGeoTrkCutMax.value.c_str(), 0, 100);
    }

    if (doprocessFlat) {
      flatchrg.add("Events/hVtxZ", "Measured vertex z position", HistType::kTH1F, {vtxzAxis});
      flatchrg.add("Events/hFlatVsMultEst", "hFlatVsMultEst", HistType::kTH2F, {flatAxis, multAxis});
      flatchrg.add("Tracks/postSel/hPVsPtEta", "; #it{p} (GeV/#it{c}); #it{p}_{T} (GeV/#it{c}); #eta;", {HistType::kTH3F, {pAxis, ptAxis, etaAxis}});
      if (cfgFillTrackQaHist || cfgFilldEdxQaHist || cfgFillNsigmaQAHist) {
        if (cfgFillTrackQaHist) {
          flatchrg.add("Tracks/postSel/hPtPhi", "; #it{p}_{T} (GeV/#it{c}); fmod(#varphi,#pi/9)", {HistType::kTH2F, {ptAxis, phiAxisMod}});
          flatchrg.add("Tracks/postSel/hPtVsWOcutDCA", "hPtVsWOcutDCA", HistType::kTH2F, {ptAxis, dcaXYAxis});
          flatchrg.add("Tracks/postSel/hPt", "", HistType::kTH1F, {ptAxis});
          flatchrg.add("Tracks/postSel/hPhi", "", HistType::kTH1F, {phiAxis});
          flatchrg.add("Tracks/postSel/hEta", "", HistType::kTH1F, {etaAxis});
          flatchrg.add("Tracks/postSel/hDCAXYvsPt", "", HistType::kTH2F, {ptAxis, dcaXYAxis});
          flatchrg.add("Tracks/postSel/hDCAZvsPt", "", HistType::kTH2F, {ptAxis, dcaZAxis});
          // tpc
          if (cfgStoreThnSparse) {
            flatchrg.add("Tracks/postSel/hPtPhiTPCCluster", "; #it{p}_{T} (GeV/#it{c}); fmod(#varphi,#pi/9); N_{cluster}", {HistType::kTHnSparseF, {ptAxis, phiAxisMod, clTpcAxis}});
          } else {
            flatchrg.add("Tracks/postSel/hPtPhiTPCCluster", "; #it{p}_{T} (GeV/#it{c}); fmod(#varphi,#pi/9); N_{cluster}", {HistType::kTH2F, {ptAxis, phiAxisMod}});
          }
          flatchrg.add("Tracks/postSel/hShTpcClvsPt", "", {HistType::kTH2F, {ptAxis, shCluserAxis}});
          flatchrg.add("Tracks/postSel/hCrossTPCvsPt", "", {HistType::kTH2F, {ptAxis, clTpcAxis}});
          flatchrg.add("Tracks/postSel/hTPCCluster", "N_{cluster}", HistType::kTH1F, {clTpcAxis});
          flatchrg.add("Tracks/postSel/tpcNClsShared", " ; # shared TPC clusters TPC", HistType::kTH1F, {{165, -0.5, 164.5}});
          flatchrg.add("Tracks/postSel/tpcCrossedRows", " ; # crossed TPC rows", HistType::kTH1F, {{165, -0.5, 164.5}});
          flatchrg.add("Tracks/postSel/tpcCrossedRowsOverFindableCls", " ; crossed rows / findable TPC clusters", HistType::kTH1F, {{60, 0.7, 1.3}});
          // its
          flatchrg.add("Tracks/postSel/itsNCls", " ; # ITS clusters", HistType::kTH1F, {{8, -0.5, 7.5}});
          flatchrg.add("Tracks/postSel/hChi2ITSTrkSegment", "chi2ITS", HistType::kTH1F, {{100, -0.5, 99.5}});
          // tof
          flatchrg.add("Tracks/postSel/hTOFPvsBeta", "Beta from TOF; #it{p} (GeV/#it{c}); #beta", {HistType::kTH2F, {pAxis, {120, 0.0, 1.2}}});
          if (cfgStoreThnSparse) {
            flatchrg.add("Tracks/postSel/hTOFpi", "Primary Pions from TOF; #eta; #it{p} (GeV/#it{c}); dEdx", {HistType::kTHnSparseF, {etaAxis, pAxis, dEdxAxis}});
          } else {
            flatchrg.add("Tracks/postSel/hTOFpi", "Primary Pions from TOF; #eta; #it{p} (GeV/#it{c}); dEdx", {HistType::kTH3F, {etaAxis, pAxis, dEdxAxis}});
          }
        }
        if (cfgFilldEdxQaHist) {
          if (cfgStoreThnSparse) {
            flatchrg.add("Tracks/postCalib/all/hMIP", "; mult; flat; #eta; #LT dE/dx #GT_{MIP, primary tracks};", {HistType::kTHnSparseF, {multAxis, flatAxis, etaAxis, dEdxAxis}});
            flatchrg.add("Tracks/postCalib/all/hPlateau", "; mult; flat; #eta; #LT dE/dx #GT_{Plateau, primary tracks};", {HistType::kTHnSparseF, {multAxis, flatAxis, etaAxis, dEdxAxis}});
          } else {
            flatchrg.add("Tracks/postCalib/all/hMIP", "; #eta; #LT dE/dx #GT_{MIP, primary tracks};", {HistType::kTH2F, {etaAxis, dEdxAxis}});
            flatchrg.add("Tracks/postCalib/all/hPlateau", "; #eta; #LT dE/dx #GT_{Plateau, primary tracks};", {HistType::kTH2F, {etaAxis, dEdxAxis}});
          }
          flatchrg.add("Tracks/postCalib/all/hMIPVsPhi", "; #varphi; #LT dE/dx #GT_{MIP, primary tracks};", {HistType::kTH2F, {phiAxis, dEdxAxis}});
          flatchrg.add("Tracks/postCalib/all/pMIPVsPhi", "; #varphi; #LT dE/dx #GT_{MIP, primary tracks};", {HistType::kTProfile, {phiAxis}});
          flatchrg.add("Tracks/postCalib/all/hMIPVsPhiVsEta", "; #varphi; #LT dE/dx #GT_{MIP, primary tracks}; #eta;", {HistType::kTH3F, {phiAxis, dEdxAxis, etaAxis}});
          flatchrg.add("Tracks/postCalib/all/hPlateauVsPhi", "; #varphi; #LT dE/dx #GT_{Plateau, primary tracks};", {HistType::kTH2F, {phiAxis, dEdxAxis}});
          flatchrg.add("Tracks/postCalib/all/pPlateauVsPhi", "; #varphi; #LT dE/dx #GT_{Plateau, primary tracks};", {HistType::kTProfile, {phiAxis}});
          flatchrg.add("Tracks/postCalib/all/hPlateauVsPhiVsEta", "; #varphi; #LT dE/dx #GT_{Plateau, primary tracks}; #eta;", {HistType::kTH3F, {phiAxis, dEdxAxis, etaAxis}});
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
              hThPtNsigmaTPC[i] = flatchrg.add<THnSparse>("Tracks/hThPtNsigmaTPC" + strID, " ; p_{T} (GeV/c)", HistType::kTHnSparseF, {ptAxis, nSigmaTPCAxis, multAxis, flatAxis});
            } else {
              hPtNsigmaTPC[i] = flatchrg.add<TH2>("Tracks/hPtNsigmaTPC" + strID, " ; p_{T} (GeV/c)", HistType::kTH2F, {ptAxis, nSigmaTPCAxis});
            }
            hPtNsigmaTOF[i] = flatchrg.add<TH2>("Tracks/hPtNsigmaTOF" + strID, " ; p_{T} (GeV/c)", HistType::kTH2F, {ptAxis, nSigmaTOFAxis});
            hPtNsigmaTPCTOF[i] = flatchrg.add<TH2>("Tracks/hPtNsigmaTPCTOF" + strID, PidChrg[i], HistType::kTH2F, {nSigmaTPCAxis, nSigmaTOFAxis});
          }
        }
      }
      flatchrg.addClone("Tracks/postSel/", "Tracks/preSel/");
      // FV0 QA
      flatchrg.add("FV0/hFV0AmplWCalib", "", {HistType::kTH2F, {channelFV0Axis, amplitudeFV0}});
      flatchrg.add("FV0/hFV0AmplvsVtxzWoCalib", "", {HistType::kTH2F, {vtxzAxis, amplitudeFV0Sum}});
      flatchrg.add("FV0/hFV0AmplvsVtxzCalib", "", {HistType::kTH2F, {vtxzAxis, amplitudeFV0Sum}});
      flatchrg.add("FV0/hFV0amp", "", {HistType::kTH2F, {channelFV0Axis, amplitudeFV0}});
      flatchrg.add("FV0/pFV0amp", "", HistType::kTProfile, {channelFV0Axis});
      flatchrg.add("FV0/hFV0ampCorr", "", {HistType::kTH2F, {channelFV0Axis, amplitudeFV0}});
      // V0's QA
      flatchrg.add("Tracks/V0qa/hV0Pt", "pT", HistType::kTH1F, {{100, 0.0f, 10}});
      flatchrg.add("Tracks/V0qa/hV0ArmPod", ";#alpha; #it{q}_T", HistType::kTH2F, {{200, -1.0f, +1.0f}, {250, 0.0f, 0.25f}});
      // dEdx PID
      flatchrg.add({"Tracks/all/hdEdx", "; #eta; mult; flat; #it{p} (GeV/#it{c}); dEdx", {HistType::kTHnSparseF, {etaAxis, multAxis, flatAxis, pAxis, dEdxAxis}}});
      // Clean samples
      if (cfgFillV0Hist) {
        if (cfgStoreThnSparse) {
          flatchrg.add({"Tracks/CleanTof/all/hPiTof", "; #eta; mult; flat; #it{p} (GeV/#it{c}); dEdx", {HistType::kTHnSparseF, {etaAxis, multAxis, flatAxis, pAxis, dEdxAxis}}});
          flatchrg.add({"Tracks/CleanV0/pos/hEV0", "; #eta; mult; flat; #it{p} (GeV/#it{c}); dEdx", {HistType::kTHnSparseF, {etaAxis, multAxis, flatAxis, pAxis, dEdxAxis}}});
          flatchrg.add({"Tracks/CleanV0/pos/hPiV0", "; #eta; mult; flat; #it{p} (GeV/#it{c}); dEdx", {HistType::kTHnSparseF, {etaAxis, multAxis, flatAxis, pAxis, dEdxAxis}}});
          flatchrg.add({"Tracks/CleanV0/pos/hPV0", "; #eta; mult; flat; #it{p} (GeV/#it{c}); dEdx", {HistType::kTHnSparseF, {etaAxis, multAxis, flatAxis, pAxis, dEdxAxis}}});
        } else {
          flatchrg.add({"Tracks/CleanTof/all/hPiTof", "; #eta; mult; flat; #it{p} (GeV/#it{c}); dEdx", {HistType::kTH3F, {etaAxis, pAxis, dEdxAxis}}});
          flatchrg.add({"Tracks/CleanV0/pos/hEV0", "; #eta; mult; flat; #it{p} (GeV/#it{c}); dEdx", {HistType::kTH3F, {etaAxis, pAxis, dEdxAxis}}});
          flatchrg.add({"Tracks/CleanV0/pos/hPiV0", "; #eta; mult; flat; #it{p} (GeV/#it{c}); dEdx", {HistType::kTH3F, {etaAxis, pAxis, dEdxAxis}}});
          flatchrg.add({"Tracks/CleanV0/pos/hPV0", "; #eta; mult; flat; #it{p} (GeV/#it{c}); dEdx", {HistType::kTH3F, {etaAxis, pAxis, dEdxAxis}}});
        }
        flatchrg.addClone("Tracks/CleanV0/pos/", "Tracks/CleanV0/neg/");
        if (cfgFillChrgType) {
          flatchrg.addClone("Tracks/CleanTof/all/", "Tracks/CleanTof/pos/");
          flatchrg.addClone("Tracks/CleanTof/all/", "Tracks/CleanTof/neg/");
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

      flatchrg.add("hEvtMcGenColls", "Number of events; Cut; #Events Passed Cut", {HistType::kTH1F, {{5, 0.5, 5.5}}});
      flatchrg.get<TH1>(HIST("hEvtMcGenColls"))->GetXaxis()->SetBinLabel(1, "Gen. coll");
      flatchrg.get<TH1>(HIST("hEvtMcGenColls"))->GetXaxis()->SetBinLabel(2, "At least 1 reco");
      flatchrg.get<TH1>(HIST("hEvtMcGenColls"))->GetXaxis()->SetBinLabel(3, "Reco. coll.");
      flatchrg.get<TH1>(HIST("hEvtMcGenColls"))->GetXaxis()->SetBinLabel(4, "Reco. good coll.");

      for (int i = 0; i < NpartChrg; i++) {
        const std::string strID = Form("/%s/%s", (i < Npart) ? "pos" : "neg", Pid[i % Npart]);
        hPtGenRecEvt[i] = flatchrg.add<TH1>("Tracks/hPtGenRecEvt" + strID, " ; p_{T} (GeV/c)", HistType::kTH1F, {ptAxis});
        hPtGenPrimRecEvt[i] = flatchrg.add<TH1>("Tracks/hPtGenPrimRecEvt" + strID, " ; p_{T} (GeV/c)", HistType::kTH1F, {ptAxis});
        hPtEffGenPrim[i] = flatchrg.add<TH3>("Tracks/hPtEffGenPrim" + strID, " ; p_{T} (GeV/c)", HistType::kTH3F, {multAxis, flatAxis, ptAxis});
        hPtEffGenWeak[i] = flatchrg.add<TH3>("Tracks/hPtEffGenWeak" + strID, " ; p_{T} (GeV/c)", HistType::kTH3F, {multAxis, flatAxis, ptAxis});
        hPtEffGenMat[i] = flatchrg.add<TH3>("Tracks/hPtEffGenMat" + strID, " ; p_{T} (GeV/c)", HistType::kTH3F, {multAxis, flatAxis, ptAxis});
        hPtEffRecPrim[i] = flatchrg.add<TH3>("Tracks/hPtEffRecPrim" + strID, " ; p_{T} (GeV/c)", HistType::kTH3F, {multAxis, flatAxis, ptAxis});
        hPtEffRecWeak[i] = flatchrg.add<TH3>("Tracks/hPtEffRecWeak" + strID, " ; p_{T} (GeV/c)", HistType::kTH3F, {multAxis, flatAxis, ptAxis});
        hPtEffRecMat[i] = flatchrg.add<TH3>("Tracks/hPtEffRecMat" + strID, " ; p_{T} (GeV/c)", HistType::kTH3F, {multAxis, flatAxis, ptAxis});
        hDCAxyBadCollPrim[i] = flatchrg.add<TH2>("Tracks/hDCAxyBadCollPrim" + strID, " ; p_{T} (GeV/c)", HistType::kTH2F, {ptAxis, dcaXYAxis});
        hDCAxyBadCollWeak[i] = flatchrg.add<TH2>("Tracks/hDCAxyBadCollWeak" + strID, " ; p_{T} (GeV/c)", HistType::kTH2F, {ptAxis, dcaXYAxis});
        hDCAxyBadCollMat[i] = flatchrg.add<TH2>("Tracks/hDCAxyBadCollMat" + strID, " ; p_{T} (GeV/c)", HistType::kTH2F, {ptAxis, dcaXYAxis});
        hPtVsDCAxyPrim[i] = flatchrg.add<TH2>("Tracks/hPtVsDCAxyPrim" + strID, " ; p_{T} (GeV/c)", HistType::kTH2F, {ptAxis, dcaXYAxis});
        hPtVsDCAxyWeak[i] = flatchrg.add<TH2>("Tracks/hPtVsDCAxyWeak" + strID, " ; p_{T} (GeV/c)", HistType::kTH2F, {ptAxis, dcaXYAxis});
        hPtVsDCAxyMat[i] = flatchrg.add<TH2>("Tracks/hPtVsDCAxyMat" + strID, " ; p_{T} (GeV/c)", HistType::kTH2F, {ptAxis, dcaXYAxis});
      }

      flatchrg.add({"hPtOut", " ; p_{T} (GeV/c)", {HistType::kTH1F, {ptAxis}}});
      flatchrg.add({"hPtOutPrim", " ; p_{T} (GeV/c)", {HistType::kTH1F, {ptAxis}}});
      flatchrg.add({"hPtOutNoEtaCut", " ; p_{T} (GeV/c)", {HistType::kTH1F, {ptAxis}}});
      flatchrg.add({"PtOutFakes", " ; p_{T} (GeV/c)", {HistType::kTH1F, {ptAxis}}});
      flatchrg.add("hPtVsDCAxyPrimAll", "hPtVsDCAxyPrimAll", HistType::kTH2F, {ptAxis, dcaXYAxis});
      flatchrg.add("hPtVsDCAxyWeakAll", "hPtVsDCAxyWeakAll", HistType::kTH2F, {ptAxis, dcaXYAxis});
      flatchrg.add("hPtVsDCAxyMatAll", "hPtVsDCAxyMatAll", HistType::kTH2F, {ptAxis, dcaXYAxis});
      flatchrg.add("hPtVsDCAxyAll", "hPtVsDCAxyAll", HistType::kTH2F, {ptAxis, dcaXYAxis});
      flatchrg.add({"ResponseGen", " ; N_{part}; F_{FV0};", {HistType::kTHnSparseF, {multAxis, flatAxis}}});
      flatchrg.add("h1flatencityFV0MCGen", "", HistType::kTH1F, {{102, -0.01, 1.01, "1-flatencityFV0"}});

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
        flatchrg.add({fmt::format(kPtMCclosurePrimF.data(), kSpeciesAll[i]).c_str(), " ; p_{T} (GeV/c)", {HistType::kTH3F, {multAxis, flatAxis, ptAxis}}});
      }
    }

    if (doprocessSgnLoss) {
      flatchrg.add("hFlatMCGenRecColl", "hFlatMCGenRecColl", {HistType::kTH1F, {flatAxis}});
      flatchrg.add("hFlatMCGen", "hFlatMCGen", {HistType::kTH1F, {flatAxis}});
      // Event counter
      flatchrg.add("hEvtMcGen", "hEvtMcGen", {HistType::kTH1F, {{4, 0.f, 4.f}}});
      flatchrg.get<TH1>(HIST("hEvtMcGen"))->GetXaxis()->SetBinLabel(1, "all");
      flatchrg.get<TH1>(HIST("hEvtMcGen"))->GetXaxis()->SetBinLabel(2, "z-vtx");
      flatchrg.get<TH1>(HIST("hEvtMcGen"))->GetXaxis()->SetBinLabel(3, "INELgt0");
      flatchrg.add("hEvtMCRec", "hEvtMCRec", {HistType::kTH1F, {{4, 0.f, 4.f}}});
      flatchrg.get<TH1>(HIST("hEvtMCRec"))->GetXaxis()->SetBinLabel(1, "all");
      flatchrg.get<TH1>(HIST("hEvtMCRec"))->GetXaxis()->SetBinLabel(2, "evt sel");
      flatchrg.get<TH1>(HIST("hEvtMCRec"))->GetXaxis()->SetBinLabel(3, "INELgt0");
      flatchrg.add("hEvtMcGenRecColl", "hEvtMcGenRecColl", {HistType::kTH1F, {{2, 0.f, 2.f}}});
      flatchrg.get<TH1>(HIST("hEvtMcGenRecColl"))->GetXaxis()->SetBinLabel(1, "INEL");
      flatchrg.get<TH1>(HIST("hEvtMcGenRecColl"))->GetXaxis()->SetBinLabel(2, "INELgt0");

      flatchrg.add("hFlatGenINELgt0", "hFlatGenINELgt0", {HistType::kTH1F, {flatAxis}});

      for (int i = 0; i < NpartChrg; ++i) {
        flatchrg.add({fmt::format(kPtGenPrimSgnF.data(), kSpecies[i]).c_str(), " ; p_{T} (GeV/c)", {HistType::kTH3F, {multAxis, flatAxis, ptAxis}}});
        flatchrg.add({fmt::format(kPtGenPrimSgnINELF.data(), kSpecies[i]).c_str(), " ; p_{T} (GeV/c)", {HistType::kTH3F, {multAxis, flatAxis, ptAxis}}});
        flatchrg.add({fmt::format(kPtRecCollPrimSgnF.data(), kSpecies[i]).c_str(), " ; p_{T} (GeV/c)", {HistType::kTH3F, {multAxis, flatAxis, ptAxis}}});
        flatchrg.add({fmt::format(kPtRecCollPrimSgnINELF.data(), kSpecies[i]).c_str(), " ; p_{T} (GeV/c)", {HistType::kTH3F, {multAxis, flatAxis, ptAxis}}});
        flatchrg.add({fmt::format(kPtGenRecCollPrimSgnF.data(), kSpecies[i]).c_str(), " ; p_{T} (GeV/c)", {HistType::kTH3F, {multAxis, flatAxis, ptAxis}}});
        flatchrg.add({fmt::format(kPtGenRecCollPrimSgnINELF.data(), kSpecies[i]).c_str(), " ; p_{T} (GeV/c)", {HistType::kTH3F, {multAxis, flatAxis, ptAxis}}});
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
    auto timestamp = bc.timestamp();
    auto runnumber = bc.runNumber();

    if (applyCalibGain) {
      fullPathCalibGain = cfgGainEqCcdbPath;
      fullPathCalibGain += "/FV0";
      auto objfv0Gain = getForTsOrRun<std::vector<float>>(fullPathCalibGain, timestamp, runnumber);
      if (!objfv0Gain) {
        for (auto i{0u}; i < kNCellsFV0; i++) {
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
  }

  using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::FT0sCorrected, aod::CentFT0As, aod::CentFT0Cs>;
  using Colls = soa::Join<aod::Collisions, aod::EvSels, aod::TPCMults, aod::PVMults, aod::MultZeqs, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>;
  using CollsGen = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::TPCMults, aod::PVMults, aod::MultZeqs, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>;
  using MCColls = soa::Join<aod::McCollisions, aod::McCentFT0Ms, aod::MultsExtraMC>;
  using CollsGenSgn = soa::SmallGroups<soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::TPCMults, aod::PVMults, aod::MultZeqs, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>>;
  using MyPIDTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::TrackSelectionExtension, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFbeta>;
  using MyLabeledTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksDCA, aod::McTrackLabels>;
  using MyFiltLabeledTracks = soa::Filtered<MyLabeledTracks>;
  using MyLabeledPIDTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::TrackSelectionExtension, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFbeta, aod::McTrackLabels>;

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
    static_assert(pidSgn == 0 || pidSgn == 1);
    static_assert(id > 0 || id < Npart);
    constexpr int kIdx = id + pidSgn * Npart;
    return mcParticle.pdgCode() == PidSgn[kIdx];
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

  template <typename T>
  bool selTOFPi(T const& track)
  {
    if (track.p() < trkSelOpt.cfgMomSelPiTOF) {
      if (track.hasTOF()) {
        if (std::abs(track.tpcNSigmaPi()) < v0SelOpt.cfgNsigmaPiTPC && std::abs(track.tofNSigmaPi()) < v0SelOpt.cfgNsigmaPiTOF) {
          return true;
        }
      }
    }
    return false;
  }

  template <o2::track::PID::ID id, typename T, typename C>
  void fillNsigma(T const& tracks, const C& collision)
  {
    const float mult = getMult(collision);
    const float flat = fillFlat<false>(collision);
    for (const auto& track : tracks) {
      checkNsigma<id>(track, mult, flat);
    }
  }

  template <typename T, typename V, typename C>
  void filldEdx(T const& tracks, V const& v0s, const C& collision, aod::BCsWithTimestamps const& /*bcs*/)
  {
    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    auto magField = (ccdbConf.cfgMagField == 0) ? getMagneticField(bc.timestamp()) : ccdbConf.cfgMagField;

    if (applyCalibGain || applyCalibVtx) {
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
          if (track.sign() * track.tpcInnerParam() > 0) {
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
          dEdx *= (50.0 / getCalibration(fDeDxVsEta, track));
          if (cfgFilldEdxQaHist) {
            if (track.sign() * track.tpcInnerParam() > 0) {
              filldEdxQA<kPos, kAfter, true>(track, collision, dEdx);
            } else {
              filldEdxQA<kNeg, kAfter, true>(track, collision, dEdx);
            }
          }
        } else {
          dEdx *= (50.0 / getCalibration<false>(fDeDxVsEta, track));
          if (cfgFilldEdxQaHist) {
            filldEdxQA<kAll, kAfter, true>(track, collision, dEdx);
          }
        }
      }

      // PID TPC dEdx
      if (cfgFillChrgType) {
        if (track.sign() * track.tpcInnerParam() > 0) {
          flatchrg.fill(HIST(kPrefix) + HIST(kCharge[kPos]) + HIST("hdEdx"), track.eta(), mult, flat, track.tpcInnerParam(), dEdx);
        } else {
          flatchrg.fill(HIST(kPrefix) + HIST(kCharge[kNeg]) + HIST("hdEdx"), track.eta(), mult, flat, track.tpcInnerParam(), dEdx);
        }
      } else {
        flatchrg.fill(HIST(kPrefix) + HIST(kCharge[kAll]) + HIST("hdEdx"), track.eta(), mult, flat, track.tpcInnerParam(), dEdx);
      }

      // TOF pions
      if (cfgFillV0Hist) {
        if (track.hasTOF() && track.beta() > 1.) {
          if (selTOFPi(track)) {
            if (cfgFillChrgType) {
              if (track.sign() * track.tpcInnerParam() > 0) {
                if (cfgStoreThnSparse) {
                  flatchrg.fill(HIST(kPrefixCleanTof) + HIST(kCharge[kPos]) + HIST("hPiTof"), track.eta(), mult, flat, track.tpcInnerParam(), dEdx);
                } else {
                  flatchrg.fill(HIST(kPrefixCleanTof) + HIST(kCharge[kPos]) + HIST("hPiTof"), track.eta(), track.tpcInnerParam(), dEdx);
                }
              } else {
                if (cfgStoreThnSparse) {
                  flatchrg.fill(HIST(kPrefixCleanTof) + HIST(kCharge[kNeg]) + HIST("hPiTof"), track.eta(), mult, flat, track.tpcInnerParam(), dEdx);
                } else {
                  flatchrg.fill(HIST(kPrefixCleanTof) + HIST(kCharge[kNeg]) + HIST("hPiTof"), track.eta(), track.tpcInnerParam(), dEdx);
                }
              }
            } else {
              if (cfgStoreThnSparse) {
                flatchrg.fill(HIST(kPrefixCleanTof) + HIST(kCharge[kAll]) + HIST("hPiTof"), track.eta(), mult, flat, track.tpcInnerParam(), dEdx);
              } else {
                flatchrg.fill(HIST(kPrefixCleanTof) + HIST(kCharge[kAll]) + HIST("hPiTof"), track.eta(), track.tpcInnerParam(), dEdx);
              }
            }
          }
        }
      }
    }

    // V0s
    if (cfgFillV0Hist) {
      for (const auto& v0 : v0s) {
        if (!isGoodV0Track(v0, tracks)) {
          continue;
        }

        const auto& posTrack = v0.template posTrack_as<T>();
        const auto& negTrack = v0.template negTrack_as<T>();
        float dEdxPos = posTrack.tpcSignal();
        float dEdxNeg = negTrack.tpcSignal();

        if (applyCalibDeDx) {
          dEdxPos *= (50.0 / getCalibration(fDeDxVsEta, posTrack));
          dEdxNeg *= (50.0 / getCalibration(fDeDxVsEta, negTrack));
        }

        if (selectTypeV0s(v0, posTrack, negTrack) == kGa) { // Gamma selection
          if (cfgStoreThnSparse) {
            flatchrg.fill(HIST(kPrefixCleanV0) + HIST(kCharge[kPos]) + HIST("hEV0"), posTrack.eta(), mult, flat, posTrack.sign() * posTrack.tpcInnerParam(), dEdxPos);
            flatchrg.fill(HIST(kPrefixCleanV0) + HIST(kCharge[kNeg]) + HIST("hEV0"), negTrack.eta(), mult, flat, negTrack.sign() * negTrack.tpcInnerParam(), dEdxNeg);
          } else {
            flatchrg.fill(HIST(kPrefixCleanV0) + HIST(kCharge[kPos]) + HIST("hEV0"), posTrack.eta(), posTrack.sign() * posTrack.tpcInnerParam(), dEdxPos);
            flatchrg.fill(HIST(kPrefixCleanV0) + HIST(kCharge[kNeg]) + HIST("hEV0"), negTrack.eta(), negTrack.sign() * negTrack.tpcInnerParam(), dEdxNeg);
          }
        }
        if (selectTypeV0s(v0, posTrack, negTrack) == kKz) { // K0S -> pi + pi
          if (cfgStoreThnSparse) {
            flatchrg.fill(HIST(kPrefixCleanV0) + HIST(kCharge[kPos]) + HIST("hPiV0"), posTrack.eta(), mult, flat, posTrack.sign() * posTrack.tpcInnerParam(), dEdxPos);
            flatchrg.fill(HIST(kPrefixCleanV0) + HIST(kCharge[kNeg]) + HIST("hPiV0"), negTrack.eta(), mult, flat, negTrack.sign() * negTrack.tpcInnerParam(), dEdxNeg);
          } else {
            flatchrg.fill(HIST(kPrefixCleanV0) + HIST(kCharge[kPos]) + HIST("hPiV0"), posTrack.eta(), posTrack.sign() * posTrack.tpcInnerParam(), dEdxPos);
            flatchrg.fill(HIST(kPrefixCleanV0) + HIST(kCharge[kNeg]) + HIST("hPiV0"), negTrack.eta(), negTrack.sign() * negTrack.tpcInnerParam(), dEdxNeg);
          }
        }
        if (selectTypeV0s(v0, posTrack, negTrack) == kLam) { // L -> p + pi-
          if (cfgStoreThnSparse) {
            flatchrg.fill(HIST(kPrefixCleanV0) + HIST(kCharge[kPos]) + HIST("hPV0"), posTrack.eta(), mult, flat, posTrack.sign() * posTrack.tpcInnerParam(), dEdxPos);
            flatchrg.fill(HIST(kPrefixCleanV0) + HIST(kCharge[kNeg]) + HIST("hPV0"), negTrack.eta(), mult, flat, negTrack.sign() * negTrack.tpcInnerParam(), dEdxNeg);
          } else {
            flatchrg.fill(HIST(kPrefixCleanV0) + HIST(kCharge[kPos]) + HIST("hPV0"), posTrack.eta(), posTrack.sign() * posTrack.tpcInnerParam(), dEdxPos);
            flatchrg.fill(HIST(kPrefixCleanV0) + HIST(kCharge[kNeg]) + HIST("hPV0"), negTrack.eta(), negTrack.sign() * negTrack.tpcInnerParam(), dEdxNeg);
          }
        }
        if (selectTypeV0s(v0, posTrack, negTrack) == kaLam) { // L -> p + pi-
          if (cfgStoreThnSparse) {
            flatchrg.fill(HIST(kPrefixCleanV0) + HIST(kCharge[kPos]) + HIST("hPiV0"), posTrack.eta(), mult, flat, posTrack.sign() * posTrack.tpcInnerParam(), dEdxPos);
            flatchrg.fill(HIST(kPrefixCleanV0) + HIST(kCharge[kNeg]) + HIST("hPiV0"), negTrack.eta(), mult, flat, negTrack.sign() * negTrack.tpcInnerParam(), dEdxNeg);
          } else {
            flatchrg.fill(HIST(kPrefixCleanV0) + HIST(kCharge[kPos]) + HIST("hPiV0"), posTrack.eta(), posTrack.sign() * posTrack.tpcInnerParam(), dEdxPos);
            flatchrg.fill(HIST(kPrefixCleanV0) + HIST(kCharge[kNeg]) + HIST("hPiV0"), negTrack.eta(), negTrack.sign() * negTrack.tpcInnerParam(), dEdxNeg);
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
      if (track.sign() * track.tpcInnerParam() > 0) {
        valCalib = fCalib.at(0)->Eval(track.eta());
      } else {
        valCalib = fCalib.at(1)->Eval(track.eta());
      }
    } else {
      valCalib = fCalib.at(2)->Eval(track.eta());
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
    return std::abs(charge) >= kMinCharge;
  }

  template <typename T>
  bool phiCut(T const& track, float mag, TF1* fphiCutLow, TF1* fphiCutHigh)
  {
    if (track.pt() < trkSelOpt.cfgPhiCutPtMin)
      return true;
    // cut to remove tracks at TPC boundaries
    double phimodn = track.phi();
    if (mag < 0) // for negative polarity field
      phimodn = o2::constants::math::TwoPI - phimodn;
    if (track.sign() < 0) // for negative charge
      phimodn = o2::constants::math::TwoPI - phimodn;
    if (phimodn < 0)
      LOGF(warning, "phi < 0: %g", phimodn);

    phimodn += o2::constants::math::PI / 18.0; // to center gap in the middle
    phimodn = std::fmod(phimodn, o2::constants::math::PI / 9.0);

    if (cfgFillTrackQaHist) {
      flatchrg.fill(HIST("Tracks/preSel/hPtPhi"), track.pt(), phimodn);
      if (track.hasTPC() && track.hasITS()) {
        if (cfgStoreThnSparse) {
          flatchrg.fill(HIST("Tracks/preSel/hPtPhiTPCCluster"), track.pt(), phimodn, track.tpcNClsFindable() - track.tpcNClsFindableMinusFound());
        } else {
          flatchrg.fill(HIST("Tracks/preSel/hPtPhiTPCCluster"), track.pt(), phimodn);
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
          flatchrg.fill(HIST("Tracks/postSel/hPtPhiTPCCluster"), track.pt(), phimodn, track.tpcNClsFindable() - track.tpcNClsFindableMinusFound());
        } else {
          flatchrg.fill(HIST("Tracks/postSel/hPtPhiTPCCluster"), track.pt(), phimodn);
        }
      }
    }
    return true;
  }

  template <typename T>
  bool isGoodTrack(T const& track, int const magfield)
  {
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
    auto nClusterTPC = track.tpcNClsFindable() - track.tpcNClsFindableMinusFound();
    if (nClusterTPC < trkSelOpt.cfgNclTPCMin) {
      return false;
    }
    flatchrg.fill(HIST("Tracks/hTrkSel"), trkSelNCls);
    if (trkSelOpt.cfgUseExtraTrkCut && !phiCut(track, magfield, fPhiCutLow, fPhiCutHigh)) {
      return false;
    }
    flatchrg.fill(HIST("Tracks/hTrkSel"), trkSelTPCBndr);
    return true;
  }

  template <typename T1, typename T2>
  int selectTypeV0s(T1 const& v0, T2 const& postrk, T2 const& negtrk)
  {
    // Gamma selection
    if (v0.mGamma() < v0SelOpt.cfgdmassG) {
      if (postrk.tpcInnerParam() < v0SelOpt.cfgV0DaughterTpcMomMax && postrk.hasTPC() && std::abs(postrk.tpcNSigmaEl()) < v0SelOpt.cfgNsigmaElTPC) {
        if (postrk.hasTOF() && std::abs(postrk.tofNSigmaEl()) < v0SelOpt.cfgNsigmaElTOF) {
          return kGa;
        }
      }
      if (negtrk.tpcInnerParam() < v0SelOpt.cfgV0DaughterTpcMomMax && negtrk.hasTPC() && std::abs(negtrk.tpcNSigmaEl()) < v0SelOpt.cfgNsigmaElTPC) {
        if (negtrk.hasTOF() && std::abs(negtrk.tofNSigmaEl()) < v0SelOpt.cfgNsigmaElTOF) {
          return kGa;
        }
      }
    }
    // K0S selection, K0S -> pi + pi
    if (std::abs(v0.mK0Short() - o2::constants::physics::MassK0Short) < v0SelOpt.cfgdmassK) {
      if (postrk.tpcInnerParam() < v0SelOpt.cfgV0DaughterTpcMomMax && postrk.hasTPC() && std::abs(postrk.tpcNSigmaPi()) < v0SelOpt.cfgNsigmaPiTPC) {
        if (postrk.hasTOF() && std::abs(postrk.tofNSigmaPi()) < v0SelOpt.cfgNsigmaElTOF) {
          return kKz;
        }
      }
      if (negtrk.tpcInnerParam() < v0SelOpt.cfgV0DaughterTpcMomMax && negtrk.hasTPC() && std::abs(negtrk.tpcNSigmaPi()) < v0SelOpt.cfgNsigmaPiTPC) {
        if (negtrk.hasTOF() && std::abs(negtrk.tofNSigmaPi()) < v0SelOpt.cfgNsigmaPiTOF) {
          return kKz;
        }
      }
    }
    // Lambda selection, L -> p + pi-
    if (std::abs(v0.mLambda() - o2::constants::physics::MassLambda0) < v0SelOpt.cfgdmassL) {
      if (postrk.tpcInnerParam() < v0SelOpt.cfgV0DaughterTpcMomMax && postrk.hasTPC() && std::abs(postrk.tpcNSigmaPr()) < v0SelOpt.cfgNsigmaPrTPC && negtrk.hasTPC() && std::abs(negtrk.tpcNSigmaPi()) < v0SelOpt.cfgNsigmaPiTPC) {
        if (postrk.hasTOF() && std::abs(postrk.tofNSigmaPr()) < v0SelOpt.cfgNsigmaPrTOF && negtrk.hasTOF() && std::abs(negtrk.tofNSigmaPi()) < v0SelOpt.cfgNsigmaPiTOF) {
          return kLam;
        }
      }
    }
    // antiLambda -> pbar + pi+
    if (std::abs(v0.mAntiLambda() - o2::constants::physics::MassLambda0) < v0SelOpt.cfgdmassL) {
      if (postrk.tpcInnerParam() < v0SelOpt.cfgV0DaughterTpcMomMax && postrk.hasTPC() && std::abs(postrk.tpcNSigmaPi()) < v0SelOpt.cfgNsigmaPiTPC && negtrk.hasTPC() && std::abs(negtrk.tpcNSigmaPr()) < v0SelOpt.cfgNsigmaPrTPC) {
        if (postrk.hasTOF() && std::abs(postrk.tofNSigmaPi()) < v0SelOpt.cfgNsigmaPiTOF && negtrk.hasTOF() && std::abs(negtrk.tofNSigmaPr()) < v0SelOpt.cfgNsigmaPrTOF) {
          return kaLam;
        }
      }
    }
    return kNaN;
  }

  template <bool fillHist = true, typename T1, typename T2>
  bool isGoodV0Track(T1 const& v0, T2 const& /*track*/)
  {
    const auto& posTrack = v0.template posTrack_as<T2>();
    const auto& negTrack = v0.template negTrack_as<T2>();

    if (std::abs(posTrack.eta()) > v0SelOpt.cfgV0etamax || std::abs(negTrack.eta()) > v0SelOpt.cfgV0etamax) {
      return false;
    }
    if (posTrack.tpcNClsFound() < v0SelOpt.cfgTPCnClsmin || negTrack.tpcNClsFound() < v0SelOpt.cfgTPCnClsmin) {
      return false;
    }
    if (posTrack.sign() * negTrack.sign() > 0) { // reject same sign pair
      return false;
    }
    if (v0.dcaV0daughters() > v0SelOpt.cfgDCAv0daughter) {
      return false;
    }
    if (v0.v0cosPA() < v0SelOpt.cfgv0cospa) {
      return false;
    }
    if (v0.v0radius() < v0SelOpt.cfgv0Rmin || v0.v0radius() > v0SelOpt.cfgv0Rmax) {
      return false;
    }
    if (std::abs(v0.dcapostopv()) < v0SelOpt.cfgDCAposToPV || std::abs(v0.dcanegtopv()) < v0SelOpt.cfgDCAposToPV) {
      return false;
    }
    if constexpr (fillHist) {
      flatchrg.fill(HIST("Tracks/V0qa/hV0Pt"), v0.pt());
      flatchrg.fill(HIST("Tracks/V0qa/hV0ArmPod"), v0.alpha(), v0.qtarm());
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

    if (track.sign() > 0) {
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
    if (track.sign() > 0) {
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
      flatchrg.fill(HIST(kPrefix) + HIST(kStatus[ft]) + HIST("hPt"), track.pt());
      flatchrg.fill(HIST(kPrefix) + HIST(kStatus[ft]) + HIST("hPhi"), track.phi());
      flatchrg.fill(HIST(kPrefix) + HIST(kStatus[ft]) + HIST("hEta"), track.eta());
      flatchrg.fill(HIST(kPrefix) + HIST(kStatus[ft]) + HIST("hDCAXYvsPt"), track.pt(), track.dcaXY());
      flatchrg.fill(HIST(kPrefix) + HIST(kStatus[ft]) + HIST("hDCAZvsPt"), track.pt(), track.dcaZ());

      if (track.hasTPC() && track.hasITS()) {
        int nFindable = track.tpcNClsFindable();
        int nMinusFound = track.tpcNClsFindableMinusFound();
        int nCluster = nFindable - nMinusFound;
        flatchrg.fill(HIST(kPrefix) + HIST(kStatus[ft]) + HIST("hTPCCluster"), nCluster);
        flatchrg.fill(HIST(kPrefix) + HIST(kStatus[ft]) + HIST("hShTpcClvsPt"), track.pt(), track.tpcFractionSharedCls());
        flatchrg.fill(HIST(kPrefix) + HIST(kStatus[ft]) + HIST("hCrossTPCvsPt"), track.pt(), track.tpcNClsFound());
        flatchrg.fill(HIST(kPrefix) + HIST(kStatus[ft]) + HIST("tpcNClsShared"), track.tpcNClsShared());
        flatchrg.fill(HIST(kPrefix) + HIST(kStatus[ft]) + HIST("tpcCrossedRows"), track.tpcNClsCrossedRows());
        flatchrg.fill(HIST(kPrefix) + HIST(kStatus[ft]) + HIST("tpcCrossedRowsOverFindableCls"), track.tpcCrossedRowsOverFindableCls());
        flatchrg.fill(HIST(kPrefix) + HIST(kStatus[ft]) + HIST("hChi2ITSTrkSegment"), track.itsChi2NCl());
        flatchrg.fill(HIST(kPrefix) + HIST(kStatus[ft]) + HIST("itsNCls"), track.itsNCls());
      }

      flatchrg.fill(HIST(kPrefix) + HIST(kStatus[ft]) + HIST("hTOFPvsBeta"), track.tpcInnerParam(), track.beta());
      if (track.beta() > trkSelOpt.cfgTOFBetaPion && track.beta() < trkSelOpt.cfgTOFBetaPion + 0.05) { // TOF pions
        flatchrg.fill(HIST(kPrefix) + HIST(kStatus[ft]) + HIST("hTOFpi"), track.eta(), track.tpcInnerParam(), track.tpcSignal());
      }

      if (std::abs(track.eta()) < trkSelOpt.cfgTrkEtaMax) {
        if (isDCAxyWoCut(track)) {
          flatchrg.fill(HIST(kPrefix) + HIST(kStatus[ft]) + HIST("hPtVsWOcutDCA"), track.pt(), track.dcaXY());
        }
      }
    }
    flatchrg.fill(HIST(kPrefix) + HIST(kStatus[ft]) + HIST("hPVsPtEta"), track.tpcInnerParam(), track.pt(), track.eta());
  }

  template <ChargeType chrg, FillType ft, bool fillHist = false, typename T, typename C>
  inline void filldEdxQA(T const& track, C const& collision, const float dEdx)
  {
    const float mult = getMult(collision);
    const float flat = fillFlat<false>(collision);
    // float dEdx = track.tpcSignal();
    if constexpr (fillHist) {
      if (track.tpcInnerParam() >= trkSelOpt.cfgMomMIPMin && track.tpcInnerParam() <= trkSelOpt.cfgMomMIPMax) {
        if (dEdx > trkSelOpt.cfgDeDxMIPMin && dEdx < trkSelOpt.cfgDeDxMIPMax) { // MIP pions
          if (cfgStoreThnSparse) {
            flatchrg.fill(HIST(kPrefix) + HIST(kStatCalib[ft]) + HIST(kCharge[chrg]) + HIST("hMIP"), mult, flat, track.eta(), dEdx);
          } else {
            flatchrg.fill(HIST(kPrefix) + HIST(kStatCalib[ft]) + HIST(kCharge[chrg]) + HIST("hMIP"), track.eta(), dEdx);
          }
          flatchrg.fill(HIST(kPrefix) + HIST(kStatCalib[ft]) + HIST(kCharge[chrg]) + HIST("hMIPVsPhi"), track.phi(), dEdx);
          flatchrg.fill(HIST(kPrefix) + HIST(kStatCalib[ft]) + HIST(kCharge[chrg]) + HIST("hMIPVsPhiVsEta"), track.phi(), dEdx, track.eta());
          flatchrg.fill(HIST(kPrefix) + HIST(kStatCalib[ft]) + HIST(kCharge[chrg]) + HIST("pMIPVsPhi"), track.phi(), dEdx);
        }
        if (dEdx > trkSelOpt.cfgDeDxMIPMax + 10. && dEdx < trkSelOpt.cfgDeDxMIPMax + 30.) { // Plateau electrons
          if (std::abs(track.beta() - 1) < trkSelOpt.cfgBetaPlateuMax) {
            if (cfgStoreThnSparse) {
              flatchrg.fill(HIST(kPrefix) + HIST(kStatCalib[ft]) + HIST(kCharge[chrg]) + HIST("hPlateau"), mult, flat, track.eta(), dEdx);
            } else {
              flatchrg.fill(HIST(kPrefix) + HIST(kStatCalib[ft]) + HIST(kCharge[chrg]) + HIST("hPlateau"), track.eta(), dEdx);
            }
            flatchrg.fill(HIST(kPrefix) + HIST(kStatCalib[ft]) + HIST(kCharge[chrg]) + HIST("hPlateauVsPhi"), track.phi(), dEdx);
            flatchrg.fill(HIST(kPrefix) + HIST(kStatCalib[ft]) + HIST(kCharge[chrg]) + HIST("hPlateauVsPhiVsEta"), track.phi(), dEdx, track.eta());
            flatchrg.fill(HIST(kPrefix) + HIST(kStatCalib[ft]) + HIST(kCharge[chrg]) + HIST("pPlateauVsPhi"), track.phi(), dEdx);
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
    return true;
  }

  int getFV0IndexPhi(int i_ch)
  {
    int iRing = -1;

    if (i_ch >= kFV0IndexPhi[0] && i_ch < kFV0IndexPhi[1]) {
      if (i_ch < kFV0IndexPhi[0] + 4) {
        iRing = i_ch;
      } else {
        if (i_ch == kFV0IndexPhi[1] - 1) {
          iRing = i_ch - 3; // 4;
        } else if (i_ch == kFV0IndexPhi[1] - 2) {
          iRing = i_ch - 1; // 5;
        } else if (i_ch == kFV0IndexPhi[1] - 3) {
          iRing = i_ch + 1; // 6;
        } else if (i_ch == kFV0IndexPhi[1] - 4) {
          iRing = i_ch + 3; // 7;
        }
      }
    } else if (i_ch >= kFV0IndexPhi[1] && i_ch < kFV0IndexPhi[2]) {
      if (i_ch < kFV0IndexPhi[2] - 4) {
        iRing = i_ch;
      } else {
        if (i_ch == kFV0IndexPhi[2] - 1) {
          iRing = i_ch - 3; // 12;
        } else if (i_ch == kFV0IndexPhi[2] - 2) {
          iRing = i_ch - 1; // 13;
        } else if (i_ch == kFV0IndexPhi[2] - 3) {
          iRing = i_ch + 1; // 14;
        } else if (i_ch == kFV0IndexPhi[2] - 4) {
          iRing = i_ch + 3; // 15;
        }
      }
    } else if (i_ch >= kFV0IndexPhi[2] && i_ch < kFV0IndexPhi[3]) {
      if (i_ch < kFV0IndexPhi[3] - 4) {
        iRing = i_ch;
      } else {
        if (i_ch == kFV0IndexPhi[3] - 1) {
          iRing = i_ch - 3; // 20;
        } else if (i_ch == kFV0IndexPhi[3] - 2) {
          iRing = i_ch - 1; // 21;
        } else if (i_ch == kFV0IndexPhi[3] - 3) {
          iRing = i_ch + 1; // 22;
        } else if (i_ch == kFV0IndexPhi[3] - 4) {
          iRing = i_ch + 3; // 23;
        }
      }
    } else if (i_ch >= kFV0IndexPhi[3] && i_ch < kFV0IndexPhi[4]) {
      if (i_ch < kFV0IndexPhi[3] + 4) {
        iRing = i_ch;
      } else {
        if (i_ch == kFV0IndexPhi[4] - 5) {
          iRing = i_ch - 3; // 28;
        } else if (i_ch == kFV0IndexPhi[4] - 6) {
          iRing = i_ch - 1; // 29;
        } else if (i_ch == kFV0IndexPhi[4] - 7) {
          iRing = i_ch + 1; // 30;
        } else if (i_ch == kFV0IndexPhi[4] - 8) {
          iRing = i_ch + 3; // 31;
        }
      }
    } else if (i_ch == kFV0IndexPhi[4]) {
      iRing = kFV0IndexPhi[4];
    } else if (i_ch == kFV0IndexPhi[4] + 8) {
      iRing = i_ch - 7; // 33;
    } else if (i_ch == kFV0IndexPhi[4] - 3) {
      iRing = i_ch + 1; // 34;
    } else if (i_ch == kFV0IndexPhi[4] + 5) {
      iRing = i_ch - 6; // 35;
    } else if (i_ch == kFV0IndexPhi[4] - 2) {
      iRing = i_ch + 2; // 36;
    } else if (i_ch == kFV0IndexPhi[4] + 6) {
      iRing = i_ch - 5; // 37;
    } else if (i_ch == kFV0IndexPhi[4] - 1) {
      iRing = i_ch + 3; // 38;
    } else if (i_ch == kFV0IndexPhi[4] + 7) {
      iRing = i_ch - 4; // 39;
    } else if (i_ch == kFV0IndexPhi[4] + 11) {
      iRing = i_ch + 7; // 40;
    } else if (i_ch == kFV0IndexPhi[4] + 3) {
      iRing = i_ch + 2; // 41;
    } else if (i_ch == kFV0IndexPhi[4] + 10) {
      iRing = i_ch - 4; // 42;
    } else if (i_ch == kFV0IndexPhi[4] + 1) {
      iRing = i_ch + 5; // 43;
    } else if (i_ch == kFV0IndexPhi[4] + 9) {
      iRing = i_ch - 1; // 44;
    } else if (i_ch == kFV0IndexPhi[4] + 1) {
      iRing = i_ch + 8; // 45;
    } else if (i_ch == kFV0IndexPhi[4] + 8) {
      iRing = i_ch + 2; // 46;
    } else if (i_ch == kFV0IndexPhi[4]) {
      iRing = i_ch + 11; // 47;
    }
    return iRing;
  }

  int getMagneticField(uint64_t timestamp)
  {
    o2::parameters::GRPMagField* grpmag = nullptr;
    grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(ccdbConf.grpmagPath, timestamp);
    if (!grpmag) {
      return 0;
    }
    return grpmag->getNominalL3Field();
  }

  template <typename C, bool isMC = false>
  float getMult(C const& collision)
  {
    float val = -999.0;
    switch (multEst) {
      case MultE::kNoMult:
        return val;
        break;
      case MultE::kMultFT0M:
        return collision.centFT0M();
        break;
      case MultE::kMultTPC:
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
      if (signals[iCell] > 0.) {
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
      if (signals[iCell] > 0.) {
        sRhoTmp += std::pow(1.0 * signals[iCell] - mRho, 2);
      }
    }
    sRhoTmp /= (1.0 * entries * entries);
    sRho = std::sqrt(sRhoTmp);
    if (mRho > 0.) {
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
    bool isOkFV0OrA = false;
    if (collision.has_foundFV0()) {
      auto fv0 = collision.foundFV0();
      std::bitset<8> fV0Triggers = fv0.triggerMask();
      isOkFV0OrA = fV0Triggers[o2::fit::Triggers::bitA];
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
          if (amplCh > 0.) {
            if (applyCalibGain) { // equalize gain channel-by-channel
              amplCh /= fv0AmplCorr[chv0];
            }
            if (chv0phi > 0) {
              fv0AmplitudeWoCalib[chv0phi] = amplCh;
              if constexpr (fillHist) {
                flatchrg.fill(HIST("FV0/hFV0AmplWCalib"), ich, fv0AmplitudeWoCalib[ich]);
              }
              if (chv0 < kInnerFV0) {
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
        if (!mTrackSelector.IsSelected(track, static_cast<TrackSelection::TrackCuts>(i))) {
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
        if (!mTrackSelector.IsSelected(track, static_cast<TrackSelection::TrackCuts>(i))) {
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

    double etaMinFV0bins[kMaxRingsFV0] = {0.0};
    double etaMaxFV0bins[kMaxRingsFV0] = {0.0};
    for (int i = 0; i < kMaxRingsFV0; ++i) {
      etaMaxFV0bins[i] = kEtaMaxFV0 - i * kDEtaFV0;
      if (i < kMaxRingsFV0 - 1) {
        etaMinFV0bins[i] = kEtaMaxFV0 - (i + 1) * kDEtaFV0;
      } else {
        etaMinFV0bins[i] = kEtaMinFV0;
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
      for (int ieta = 0; ieta < kMaxRingsFV0; ieta++) {

        nFV0sectors = kNCellsFV0 / 6.;
        if (ieta == kMaxRingsFV0 - 1) {
          nFV0sectors = kNCellsFV0 / 3.;
        }

        for (int iphi = 0; iphi < nFV0sectors; iphi++) {

          minPhi = iphi * TwoPI / nFV0sectors;
          maxPhi = (iphi + 1) * TwoPI / nFV0sectors;
          dPhi = std::abs(maxPhi - minPhi);

          if (etaMc >= etaMinFV0bins[ieta] && etaMc < etaMaxFV0bins[ieta] && phiMc >= minPhi && phiMc < maxPhi) {
            rhoLatticeFV0[isegment] += 1. / std::abs(kDEtaFV0 * dPhi);
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
    constexpr int kHistIdx = id + pidSgn * Npart;
    auto kIdx = static_cast<int>(id);
    const std::string strID = Form("/%s/%s", (pidSgn == 0 && id < Npart) ? "pos" : "neg", Pid[kIdx]);
    hPtEffRec[kHistIdx] = flatchrg.add<TH1>("Tracks/hPtEffRec" + strID, " ; p_{T} (GeV/c)", HistType::kTH1F, {ptAxis});
    hPtEffGen[kHistIdx] = flatchrg.add<TH1>("Tracks/hPtEffGen" + strID, " ; p_{T} (GeV/c)", HistType::kTH1F, {ptAxis});
  }

  template <int pidSgn, o2::track::PID::ID id>
  void initEfficiency()
  {
    static_assert(pidSgn == 0 || pidSgn == 1);
    static_assert(id > 0 || id < Npart);
    constexpr int kIdx = id + pidSgn * Npart;
    const TString partName = PidChrg[kIdx];
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

    const int kHistIdx = id + pidSgn * Npart;
    bookEff("hEffvsPt", hPtEffRec[kHistIdx]);
  }

  template <int pidSgn, o2::track::PID::ID id>
  void fillEfficiency()
  {
    static_assert(pidSgn == 0 || pidSgn == 1);
    constexpr int kHistIdx = id + pidSgn * Npart;
    const char* partName = PidChrg[kHistIdx];
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
    fillEff("hEffvsPt", hPtEffRec[kHistIdx], hPtEffGen[kHistIdx]);
  }

  template <int pidSgn, o2::track::PID::ID id>
  void fillMCRecTrack(MyLabeledPIDTracks::iterator const& track, const float mult, const float flat)
  {
    static_assert(pidSgn == 0 || pidSgn == 1);
    constexpr int kHistIdx = id + pidSgn * Npart;
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
        if (mcParticle.getProcess() == kProcessIdWeak) {
          hDCAxyBadCollWeak[kHistIdx]->Fill(track.pt(), track.dcaXY());
        } else {
          hDCAxyBadCollMat[kHistIdx]->Fill(track.pt(), track.dcaXY());
        }
      } else {
        hDCAxyBadCollPrim[kHistIdx]->Fill(track.pt(), track.dcaXY());
      }
    }

    if (!isDCAxyCut(track)) {
      return;
    }
    flatchrg.fill(HIST("hPtVsDCAxyAll"), track.pt(), track.dcaXY());

    if (selTPCtrack(track)) {
      hPtEffRec[kHistIdx]->Fill(mcParticle.pt());
    }

    if (!mcParticle.isPhysicalPrimary()) {
      if (mcParticle.getProcess() == kProcessIdWeak) {
        hPtEffRecWeak[kHistIdx]->Fill(mult, flat, track.pt());
        hPtVsDCAxyWeak[kHistIdx]->Fill(track.pt(), track.dcaXY());
        flatchrg.fill(HIST("hPtVsDCAxyWeakAll"), track.pt(), track.dcaXY());
      } else {
        hPtEffRecMat[kHistIdx]->Fill(mult, flat, track.pt());
        hPtVsDCAxyMat[kHistIdx]->Fill(track.pt(), track.dcaXY());
        flatchrg.fill(HIST("hPtVsDCAxyMatAll"), track.pt(), track.dcaXY());
      }
    } else {
      hPtEffRecPrim[kHistIdx]->Fill(mult, flat, track.pt());
      hPtVsDCAxyPrim[kHistIdx]->Fill(track.pt(), track.dcaXY());
      flatchrg.fill(HIST("hPtVsDCAxyPrimAll"), track.pt(), track.dcaXY());
    }
  }

  template <int pidSgn, o2::track::PID::ID id, bool recoEvt = false>
  void fillMCGen(aod::McParticles::iterator const& mcParticle, const float mult, const float flat)
  {
    static_assert(pidSgn == 0 || pidSgn == 1);
    constexpr int kHistIdx = id + pidSgn * Npart;

    if (!isPID<pidSgn, id>(mcParticle)) {
      return;
    }

    if constexpr (recoEvt) {
      hPtGenRecEvt[kHistIdx]->Fill(mcParticle.pt());
      if (mcParticle.isPhysicalPrimary()) {
        hPtGenPrimRecEvt[kHistIdx]->Fill(mcParticle.pt());
      }
      return;
    }

    if (!mcParticle.isPhysicalPrimary()) {
      if (mcParticle.getProcess() == kProcessIdWeak) {
        hPtEffGenWeak[kHistIdx]->Fill(mult, flat, mcParticle.pt());
      } else {
        hPtEffGenMat[kHistIdx]->Fill(mult, flat, mcParticle.pt());
      }
    } else {
      hPtEffGenPrim[kHistIdx]->Fill(mult, flat, mcParticle.pt());
      hPtEffGen[kHistIdx]->Fill(mcParticle.pt());
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
        constexpr int kIdx = i.value;
        if (particle.pdgCode() == PidSgn[kIdx]) {
          flatchrg.fill(HIST(kPrefix) + HIST(kSpecies[kIdx]) + HIST(kPtGenPrimSgn), mult, flat, particle.pt());
          if (isINELgt0mc) {
            flatchrg.fill(HIST(kPrefix) + HIST(kSpecies[kIdx]) + HIST(kPtGenPrimSgnINEL), mult, flat, particle.pt());
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
          constexpr int kIdx = i.value;
          if (particle.pdgCode() == PidSgn[kIdx]) {
            flatchrg.fill(HIST(kPrefix) + HIST(kSpecies[kIdx]) + HIST(kPtRecCollPrimSgn), mult, flat, particle.pt());
            if (nRecCollINELgt0) {
              flatchrg.fill(HIST(kPrefix) + HIST(kSpecies[kIdx]) + HIST(kPtRecCollPrimSgnINEL), mult, flat, particle.pt());
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
    if (nRecCollINELgt0 > 0) {
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
        constexpr int kIdx = i.value;
        if (particle.pdgCode() == PidSgn[kIdx]) {
          flatchrg.fill(HIST(kPrefix) + HIST(kSpecies[kIdx]) + HIST(kPtGenRecCollPrimSgn), mult, flat, particle.pt());
          if (nRecCollINELgt0) {
            flatchrg.fill(HIST(kPrefix) + HIST(kSpecies[kIdx]) + HIST(kPtGenRecCollPrimSgnINEL), mult, flat, particle.pt());
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
      auto bc = coll.template bc_as<aod::BCsWithTimestamps>();
      auto magField = (ccdbConf.cfgMagField == 0) ? getMagneticField(bc.timestamp()) : ccdbConf.cfgMagField;

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
        constexpr int kIdx = i.value;
        if ((std::abs(o2::aod::pidutils::tpcNSigma<kIdx>(track)) < trkSelOpt.cfgNsigmaMax) && std::abs(track.rapidity(o2::track::PID::getMass(kIdx))) <= trkSelOpt.cfgRapMax) {
          if (std::fabs(mcParticle.pdgCode()) == PDGs[kIdx]) {
            flatchrg.fill(HIST(kPrefix) + HIST(kSpeciesAll[kIdx]) + HIST(kPtMCclosurePrim), multRec, flatRec, track.pt());
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
