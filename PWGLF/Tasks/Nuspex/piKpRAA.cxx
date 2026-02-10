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

/// \file piKpRAA.cxx
///
/// \brief task for analysis of piKp RAA
/// \author Omar Vazquez (omar.vazquez.rueda@cern.ch)
/// \since August 10, 2025

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/MathConstants.h"
#include "CommonConstants/ZDCConstants.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h" // required for Filter op.
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "ReconstructionDataFormats/Track.h"
#include <CCDB/BasicCCDBManager.h>

#include "TMCProcess.h"
#include "TPDGCode.h"
#include "TVector3.h"
#include <TString.h>

#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <string>
#include <string_view>
#include <typeinfo>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::evsel;
using namespace o2::constants::math;
using namespace o2::framework::expressions;
using namespace o2::aod::rctsel;

using ColEvSels = soa::Join<aod::Collisions, aod::EvSels, aod::FT0MultZeqs, o2::aod::CentFT0Cs, o2::aod::CentFT0Ms, aod::TPCMults, o2::aod::BarrelMults>;
using BCsRun3 = soa::Join<aod::BCsWithTimestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;

using ColEvSelsMC = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, o2::aod::CentFT0Cs, o2::aod::CentFT0Ms, o2::aod::CentFT0Ms, o2::aod::BarrelMults>;

using TracksFull = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelectionExtension, aod::TracksDCA, aod::TrackSelection, aod::TracksCovIU, aod::pidTPCPi, aod::pidTPCPr, aod::pidTOFPr, aod::pidTPCEl, aod::pidTOFFlags, aod::pidTOFbeta, aod::TOFSignal, aod::pidTOFFullPi, aod::pidTOFFullEl>;

using TracksMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelectionExtension, aod::TracksDCA, aod::TrackSelection, aod::TracksCovIU, aod::pidTPCPi, aod::pidTPCPr, aod::pidTOFPr, aod::pidTPCEl, aod::pidTOFFlags, aod::pidTOFbeta, aod::TOFSignal, aod::pidTOFFullPi, aod::pidTOFFullEl, aod::McTrackLabels>;
// using SimTracks = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels>;

static constexpr int kNEtaHists{8};

std::array<std::shared_ptr<TH3>, kNEtaHists> dEdxPiV0{};
std::array<std::shared_ptr<TH3>, kNEtaHists> dEdxPrV0{};
std::array<std::shared_ptr<TH3>, kNEtaHists> dEdxElV0{};
std::array<std::shared_ptr<TH3>, kNEtaHists> dEdxPiTOF{};
std::array<std::shared_ptr<TH3>, kNEtaHists> dEdx{};
std::array<std::shared_ptr<TH2>, kNEtaHists> nClVsdEdxPiV0{};
std::array<std::shared_ptr<TH2>, kNEtaHists> nClVsdEdxElV0{};
std::array<std::shared_ptr<TH2>, kNEtaHists> nClVsdEdxPrV0{};
std::array<std::shared_ptr<TH2>, kNEtaHists> pTVsP{};
std::array<std::shared_ptr<TH2>, kNEtaHists> nClVsP{};
std::array<std::shared_ptr<TH2>, kNEtaHists> nClVsPElV0{};
std::array<std::shared_ptr<TH2>, kNEtaHists> nClVsPPiV0{};
std::array<std::shared_ptr<TH2>, kNEtaHists> nClVsPPrV0{};
std::array<std::shared_ptr<TProfile>, kNEtaHists> nClVsPp{};
std::array<std::shared_ptr<TProfile>, kNEtaHists> nClVsPpElV0{};
std::array<std::shared_ptr<TProfile>, kNEtaHists> nClVsPpPiV0{};
std::array<std::shared_ptr<TProfile>, kNEtaHists> nClVsPpPrV0{};
std::array<std::shared_ptr<TProfile>, kNEtaHists> nClVsdEdxpPiV0{};
std::array<std::shared_ptr<TProfile>, kNEtaHists> nClVsdEdxpElV0{};
std::array<std::shared_ptr<TProfile>, kNEtaHists> nClVsdEdxpPrV0{};

struct PiKpRAA {

  static constexpr int kZeroInt{0};
  static constexpr int kSevenInt{7};

  static constexpr float kZero{0.0f};
  static constexpr float kOne{1.0f};
  static constexpr float kTwoPtGeVSel{2.0f};
  static constexpr float kThree{3.0f};
  static constexpr float kTenToMinusNine{1e-9};
  static constexpr float kMinPtNchSel{0.1f};
  static constexpr float kMaxPtNchSel{3.0f};
  static constexpr float kMinCharge{3.f};
  static constexpr float kMinPElMIP{0.3f};
  static constexpr float kMaxPElMIP{0.45f};
  static constexpr float kMinPMIP{0.4f};
  static constexpr float kMaxPMIP{0.6f};
  static constexpr float kMindEdxMIP{40.0f};
  static constexpr float kMaxdEdxMIP{60.0f};
  static constexpr float kMindEdxMIPPlateau{70.0f};
  static constexpr float kMaxdEdxMIPPlateau{90.0f};
  static constexpr float kMinFT0A{3.5f};
  static constexpr float kMaxFT0A{4.9f};
  static constexpr float kMinFT0C{-3.3f};
  static constexpr float kMaxFT0C{-2.1f};

  static constexpr float kLowEta[kNEtaHists] = {-0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6};
  static constexpr float kHighEta[kNEtaHists] = {-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8};
  // static constexpr float kLowEta[kNEtaHists] = {0.0, 0.2, 0.4, 0.6};
  // static constexpr float kHighEta[kNEtaHists] = {0.2, 0.4, 0.6, 0.8};

  static constexpr float DefaultLifetimeCuts[1][2] = {{30., 20.}};
  Configurable<LabeledArray<float>> lifetimecut{"lifetimecut", {DefaultLifetimeCuts[0], 2, {"lifetimecutLambda", "lifetimecutK0S"}}, "lifetimecut"};

  struct : ConfigurableGroup {
    Configurable<uint8_t> v0TypeSelection{"v0TypeSelection", 1, "select on a certain V0 type (leave negative if no selection desired)"};
    Configurable<bool> useOfficialV0sSelOfDaughters{"useOfficialV0sSelOfDaughters", true, "Use the same track selection for daughters as the V0s analysis in OO"};

    // Selection criteria: acceptance
    Configurable<float> rapidityCut{"rapidityCut", 0.5, "rapidity"};
    Configurable<float> minEtaDaughter{"minEtaDaughter", -0.8, "Daughter minimum-eta selection"};
    Configurable<float> maxEtaDaughter{"maxEtaDaughter", +0.8, "Daughter maximum-eta selection"};
    Configurable<float> minPt{"minPt", 0.15, "minimum pt of the tracks"};
    Configurable<float> maxPt{"maxPt", 20.0, "maximum pt of the tracks"};
    Configurable<float> minPtDaughter{"minPtDaughter", 0.15, "minimum pt of the tracks"};
    Configurable<float> maxPtDaughter{"maxPtDaughter", 20.0, "maximum pt of the tracks"};
    Configurable<bool> useNclsPID{"useNclsPID", true, "Use Ncl for PID?"};
    Configurable<int16_t> minNcl{"minNcl", 135, "minimum found Ncl in TPC"};
    Configurable<int16_t> minNCrossedRows{"minNCrossedRows", 70, "minimum number of crossed rows"};
    Configurable<float> minNCrossedRowsOverFindableCls{"minNCrossedRowsOverFindableCls", 0.8, "min N crossed rows over findable Cls"};
    Configurable<float> maxChi2ClsTPC{"maxChi2ClsTPC", 4.0, "Max chi2 per Cls TPC"};
    Configurable<float> maxChi2ClsITS{"maxChi2ClsITS", 36.0, "chi2 per Cls ITS"};
    Configurable<float> maxDCAZ{"maxDCAZ", 2.0, "Max DCA Z"};
    Configurable<bool> itsRefit{"itsRefit", true, "Require ITS refit"};
    Configurable<bool> tpcRefit{"tpcRefit", true, "Require TPC refit"};
    Configurable<bool> chi2Golden{"chi2Golden", true, "Require Chi2 golde selection"};
    Configurable<bool> its1HitIB{"its1HitIB", true, "Require one hit in the ITS IB"};
    Configurable<bool> requireITShit{"requireITShit", true, "Apply requirement of one hit in the ITS IB?"};

    // Standard 5 topological criteria
    Configurable<float> v0cospa{"v0cospa", 0.995, "min V0 CosPA"};
    Configurable<float> dcav0dau{"dcav0dau", 1.0, "max DCA V0 Daughters (cm)"};
    Configurable<float> dcanegtopv{"dcanegtopv", .1, "min DCA Neg To PV (cm)"};
    Configurable<float> dcapostopv{"dcapostopv", .1, "min DCA Pos To PV (cm)"};
    Configurable<float> dcaProtonFromLambda{"dcaProtonFromLambda", .05, "min DCA Proton (from Lambda) To PV (cm)"};
    Configurable<float> dcaPionFromLambda{"dcaPionFromLambda", .2, "min DCA Pion (from Lambda) To PV (cm)"};
    Configurable<float> dcaElectronFromGamma{"dcaElectronFromGamma", .1, "min DCA Electron (from Gamma conversion) To PV (cm)"};

    Configurable<float> v0radius{"v0radius", 1.2, "minimum V0 radius (cm)"};
    Configurable<float> v0radiusMax{"v0radiusMax", 1E5, "maximum V0 radius (cm)"};
    Configurable<float> lifeTimeCutK0s{"lifeTimeCutK0s", 20.0, "lifetime cut K0s (cm)"};
    Configurable<float> lifeTimeCutLambda{"lifeTimeCutLambda", 30.0, "lifetime cut Lambda and AntiLambda (cm)"};

    // Additional selection on the AP plot (exclusive for K0Short)
    // original equation: lArmPt*5>TMath::Abs(lArmAlpha)
    Configurable<float> armPodCut{"armPodCut", 5.0f, "pT * (cut) > |alpha|, AP cut. Negative: no cut"};
    Configurable<float> armAlphaSel{"armAlphaSel", 0.45f, "Armenteros alpha selection (Gammas)"};
    Configurable<float> qTSel{"qTSel", 0.01f, "Armenteros qT select (Gammas)"};

    // Selection
    Configurable<bool> applyInvMassSel{"applyInvMassSel", true, "Select V0s close to the Inv. mass value"};
    Configurable<bool> selElecFromGammas{"selElecFromGammas", true, "track selection for electrons"};
    Configurable<float> dMassSel{"dMassSel", 0.01f, "Invariant mass selection"};
    Configurable<float> dMassSelG{"dMassSelG", 0.1f, "Inv mass selection gammas"};

    // PID (TPC/TOF)
    Configurable<float> dEdxPlateauSel{"dEdxPlateauSel", 50, "dEdx selection for electrons"};
    Configurable<float> tpcPidNsigmaCut{"tpcPidNsigmaCut", 5, "tpcPidNsigmaCut"};
    Configurable<float> maxExpTOFPi{"maxExpTOFPi", 0.00005, "Maximum beta TOF selection"};

    // Phi cut
    Configurable<bool> applyPhiCut{"applyPhiCut", false, "Apply geometrical cut?"};
    Configurable<bool> applyEtaCal{"applyEtaCal", false, "Apply eta calibration?"};
    Configurable<bool> applyPlateauSel{"applyPlateauSel", false, "Apply eta calibration?"};
    Configurable<bool> usePinPhiSelection{"usePinPhiSelection", true, "Uses Phi selection as a function of P or Pt?"};
    Configurable<bool> applyNclSel{"applyNclSel", false, "Apply Min. found Ncl in TPC?"};
  } v0Selections;

  // Configurables Event Selection
  Configurable<bool> isNoCollInTimeRangeStrict{"isNoCollInTimeRangeStrict", true, "use isNoCollInTimeRangeStrict?"};
  Configurable<bool> selNoSameBunchPileup{"selNoSameBunchPileup", true, "selNoSameBunchPileup?"};
  Configurable<bool> selIsGoodZvtxFT0vsPV{"selIsGoodZvtxFT0vsPV", true, "selIsGoodZvtxFT0vsPV?"};
  Configurable<bool> isNoCollInTimeRangeStandard{"isNoCollInTimeRangeStandard", false, "use isNoCollInTimeRangeStandard?"};
  Configurable<bool> isNoCollInRofStrict{"isNoCollInRofStrict", true, "use isNoCollInRofStrict?"};
  Configurable<bool> isNoCollInRofStandard{"isNoCollInRofStandard", false, "use isNoCollInRofStandard?"};
  Configurable<bool> isNoHighMultCollInPrevRof{"isNoHighMultCollInPrevRof", true, "use isNoHighMultCollInPrevRof?"};
  Configurable<bool> isNoCollInTimeRangeNarrow{"isNoCollInTimeRangeNarrow", false, "use isNoCollInTimeRangeNarrow?"};
  Configurable<bool> isOccupancyCut{"isOccupancyCut", true, "Occupancy cut?"};
  Configurable<bool> isCentSel{"isCentSel", true, "Centrality selection?"};
  Configurable<bool> selHasFT0{"selHasFT0", true, "Has FT0?"};
  Configurable<bool> isT0Ccent{"isT0Ccent", true, "Use T0C-based centrality?"};
  Configurable<bool> isZvtxPosSel{"isZvtxPosSel", true, "Zvtx position selection?"};
  Configurable<bool> isZvtxPosSelMC{"isZvtxPosSelMC", true, "Zvtx position selection for MC events?"};
  Configurable<bool> selTVXMC{"selTVXMC", true, "apply TVX selection in MC?"};
  Configurable<bool> selINELgt0{"selINELgt0", true, "Select INEL > 0?"};
  Configurable<bool> isApplyFT0CbasedOccupancy{"isApplyFT0CbasedOccupancy", false, "T0C Occu cut"};
  Configurable<bool> applyNchSel{"applyNchSel", false, "Use mid-rapidity-based Nch selection"};
  Configurable<bool> skipRecoColGTOne{"skipRecoColGTOne", true, "Remove collisions if reconstructed more than once"};

  // Event selection
  Configurable<float> posZcut{"posZcut", +10.0, "z-vertex position cut"};
  Configurable<float> minT0CcentCut{"minT0CcentCut", 0.0, "Min T0C Cent. cut"};
  Configurable<float> maxT0CcentCut{"maxT0CcentCut", 100.0, "Max T0C Cent. cut"};
  Configurable<float> minOccCut{"minOccCut", 0., "min Occu cut"};
  Configurable<float> maxOccCut{"maxOccCut", 500., "max Occu cut"};
  Configurable<float> nSigmaNchCut{"nSigmaNchCut", 3., "nSigma Nch selection"};

  ConfigurableAxis binsPtPhiCut{"binsPtPhiCut", {VARIABLE_WIDTH, 0.0, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.6, 3.0, 3.5, 4.0, 5.0, 7.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 45.0, 50.0}, "pT"};
  ConfigurableAxis binsPtV0s{"binsPtV0s", {VARIABLE_WIDTH, 0, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.5, 3.0, 3.5, 4, 5, 7, 9, 12, 15, 20}, "pT"};
  ConfigurableAxis binsPtNcl{"binsPtNcl", {VARIABLE_WIDTH, 0.0, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.0, 9.0, 12.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0}, "pT"};
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.12}, "pT binning"};
  ConfigurableAxis binsCent{"binsCent", {VARIABLE_WIDTH, 0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.}, "T0C binning"};
  ConfigurableAxis binsZpos{"binsZpos", {60, -30.0, 30.0}, "Z pos axis"};
  ConfigurableAxis axisEta{"axisEta", {50, -1.0, 1.0}, "Eta axis"};
  ConfigurableAxis axisY{"axisY", {50, -1.0, 1.0}, "rapidity axis"};
  ConfigurableAxis axisArmAlpha{"axisArmAlpha", {200, -1.0, 1.0}, "Armenteros alpha"};
  ConfigurableAxis axisArmqT{"axisArmqT", {600, 0.0f, 0.3f}, "Armenteros qT"};
  ConfigurableAxis axisK0Mass{"axisK0Mass", {200, 0.4f, 0.6f}, "Mass K0Short"};
  ConfigurableAxis axisLambdaMass{"axisLambdaMass", {200, 1.101f, 1.131f}, "Mass Lambda"};
  ConfigurableAxis axisGammaMass{"axisGammaMass", {150, 0.0f, 0.15f}, "Mass Gamma"};
  ConfigurableAxis axisNsigmaTPC{"axisNsigmaTPC", {200, -10.0f, 10.0f}, "N sigma TPC"};
  ConfigurableAxis axisdEdx{"axisdEdx", {140, 20.0, 160.0}, "dEdx binning"};
  ConfigurableAxis axisDCAxy{"axisDCAxy", {105, -1.05f, 1.05f}, "DCAxy axis"};
  Configurable<int> nBinsNch{"nBinsNch", 400, "N bins Nch (|eta|<0.8)"};
  Configurable<int> nBinsNPV{"nBinsNPV", 600, "N bins ITS tracks"};
  Configurable<float> minNch{"minNch", 0, "Min Nch (|eta|<0.8)"};
  Configurable<float> maxNch{"maxNch", 400, "Max Nch (|eta|<0.8)"};
  Configurable<float> minNpv{"minNpv", 0, "Min NPV"};
  Configurable<float> maxNpv{"maxNpv", 600, "Max NPV"};

  // CCDB paths
  Configurable<std::string> pathMeanNch{"pathMeanNch", "Users/o/omvazque/MeanNch/OO/Pass2/PerTimeStamp/Aug20", "base path to the ccdb object"};
  Configurable<std::string> pathSigmaNch{"pathSigmaNch", "Users/o/omvazque/SigmaNch/OO/Pass2/PerTimeStamp/Aug20", "base path to the ccdb object"};
  Configurable<std::string> pathEtaCal{"pathEtaCal", "Users/o/omvazque/EtaCal/OO/Global", "base path to the ccdb object"};
  Configurable<std::string> pathEtaCalPlateau{"pathEtaCalPlateau", "Users/o/omvazque/EtaCal/OO/Global", "base path to the ccdb object"};
  Configurable<std::string> pathPhiCutHigh{"pathPhiCutHigh", "Users/o/omvazque/PhiCut/OO/Global/High", "base path to the ccdb object"};
  Configurable<std::string> pathPhiCutLow{"pathPhiCutLow", "Users/o/omvazque/PhiCut/OO/Global/Low", "base path to the ccdb object"};
  Configurable<int64_t> ccdbNoLaterThan{"ccdbNoLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> rctLabel{"rctLabel", "CBT_hadronPID", "RCT selection flag (CBT, CBT_hadronPID, CBT_electronPID, CBT_calo, CBT_muon, CBT_muon_glo)"};
  Configurable<bool> rctCheckZDC{"rctCheckZDC", false, "RCT flag to check whether the ZDC is present or not"};
  Configurable<bool> rctTreatLimitedAcceptanceAsBad{"rctTreatLimitedAcceptanceAsBad", false, "RCT flag to reject events with limited acceptance for selected detectors"};
  Configurable<bool> requireGoodRct{"requireGoodRct", true, "RCT flag to reject events with limited acceptance for selected detectors"};
  Configurable<bool> requireGoodPIDRct{"requireGoodPIDRct", true, "RCT flag to reject events with limited acceptance for selected detectors"};

  // RCT Checker instance
  RCTFlagsChecker rctChecker;

  enum EvCutLabel {
    All = 1,
    SelEigth,
    NoSameBunchPileup,
    IsGoodZvtxFT0vsPV,
    NoCollInTimeRangeStrict,
    NoCollInTimeRangeStandard,
    NoCollInRofStrict,
    NoCollInRofStandard,
    NoHighMultCollInPrevRof,
    NoCollInTimeRangeNarrow,
    OccuCut,
    HasFT0,
    Centrality,
    VtxZ,
    NchSel,
    INELgt0
  };

  enum TrkSelLabel {
    AllTrks = 1,
    Eta,
    Pt,
    XRows,
    XRowsOverFindableCls,
    Chi2TPC,
    Chi2ITS,
    Itsrefit,
    Tpcrefit,
    Golden,
    Itshit,
    PassedAll
  };

  enum V0sCounter {
    K0s = 1,
    Lambda,
    AntiLambda,
    Gamma
  };

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  Service<ccdb::BasicCCDBManager> ccdb;

  struct ConfigNch {
    TH1F* hMeanNch = nullptr;
    TH1F* hSigmaNch = nullptr;
    bool calibrationsLoaded = false;
  } cfgNch;

  struct ConfigPhiCut {
    TH1F* hPhiCutHigh = nullptr;
    TH1F* hPhiCutLow = nullptr;
    bool isPhiCutLoaded = false;
  } phiCut;

  struct ConfigEtaCalib {
    TProfile* pEtaCal = nullptr;
    TProfile* pEtaCalPlateau = nullptr;
    bool isMIPCalLoaded = false;
    bool isCalPlateauLoaded = false;
  } etaCal;

  TrackSelection trkSelGlobalOpenDCAxy;
  // TrackSelection trkSelDaugthers;
  TrackSelection trkSelGlobal;
  // TrackSelection trkSelDaugthersV0s() {
  //     TrackSelection selectedTracks;
  //     selectedTracks.SetEtaRange(-0.8f, 0.8f);
  //     selectedTracks.SetMinNCrossedRowsTPC(70);
  //     return selectedTracks;
  // }

  TrackSelection trkSelOpenDCAxy()
  {
    TrackSelection selectedTracks;
    selectedTracks.SetTrackType(o2::aod::track::TrackTypeEnum::Track); // Run 2 track asked by default
    selectedTracks.SetPtRange(0.1f, 1e10f);
    selectedTracks.SetEtaRange(-0.8f, 0.8f);
    selectedTracks.SetRequireITSRefit(true);
    selectedTracks.SetRequireTPCRefit(true);
    selectedTracks.SetRequireGoldenChi2(true);
    selectedTracks.SetMinNCrossedRowsTPC(70);
    selectedTracks.SetMinNCrossedRowsOverFindableClustersTPC(0.8f);
    selectedTracks.SetMaxChi2PerClusterTPC(4.f);
    selectedTracks.SetRequireHitsInITSLayers(1, {0, 1}); // one hit in any SPD layer
    selectedTracks.SetMaxChi2PerClusterITS(36.f);
    // selectedTracks.SetMaxDcaXYPtDep([](float pt) { return 0.0105f + 0.0350f / pow(pt, 1.1f); });
    selectedTracks.SetMaxDcaZ(2.f);
    return selectedTracks;
  }

  int currentRunNumberNchSel;
  int currentRunNumberPhiSel;
  void init(InitContext const&)
  {

    // Initialize the rct checker
    if (requireGoodRct) {
      rctChecker.init(rctLabel.value, rctCheckZDC.value, rctTreatLimitedAcceptanceAsBad.value);
    }

    currentRunNumberNchSel = -1;
    currentRunNumberPhiSel = -1;
    trkSelGlobalOpenDCAxy = trkSelOpenDCAxy();
    // trkSelDaugthers = trkSelDaugthersV0s();
    trkSelGlobal = getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSibAny, TrackSelection::GlobalTrackRun3DCAxyCut::Default);

    // define axes you want to use
    const std::string titlePorPt{v0Selections.usePinPhiSelection ? "#it{p} (GeV/#it{c})" : "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisZpos{binsZpos, "Vtx_{z} (cm)"};
    const AxisSpec axisEvent{17, 0.5, 17.5, ""};
    const AxisSpec axisNcl{161, -0.5, 160.5, "#it{N}_{cl} TPC"};
    const AxisSpec axisPt{binsPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisPtV0s{binsPtV0s, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisPtNcl{binsPtNcl, Form("%s", titlePorPt.data())};
    const AxisSpec axisXPhiCut{binsPtPhiCut, Form("%s", titlePorPt.data())};
    const AxisSpec axisCent{binsCent, "Centrality Perc."};
    const char* endingEta[kNEtaHists] = {"86", "64", "42", "20", "02", "24", "46", "68"};
    const char* latexEta[kNEtaHists] = {"-0.8<#eta<-0.6", "-0.6<#eta<-0.4", "-0.4<#eta<-0.2", "-0.2<#eta<0", "0<#eta<0.2", "0.2<#eta<0.4", "0.4<#eta<0.6", "0.6<#eta<0.8"};

    registry.add("EventCounter", ";;Events", kTH1F, {axisEvent});
    registry.add("zPos", "With Event Selection;;Entries;", kTH1F, {axisZpos});
    registry.add("T0Ccent", ";;Entries", kTH1F, {axisCent});
    registry.add("T0CcentVsFoundFT0", "Found(=1.5) NOT Found(=0.5);;Status;", kTH2F, {{{axisCent}, {2, 0, 2}}});
    registry.add("NchVsCent", "Measured Nch v.s. Centrality (At least Once Rec. Coll. + Sel. criteria);;Nch", kTH2F, {{axisCent, {nBinsNch, minNch, maxNch}}});
    registry.add("NclVsEtaPID", ";#eta;Ncl used for PID", kTH2F, {{{axisEta}, {161, -0.5, 160.5}}});
    registry.add("NclVsEtaPIDp", ";#eta;#LTNcl#GT used for PID", kTProfile, {axisEta});
    registry.add("dcaVsPtPi", "Primary pions;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm);Centrality Perc.;", kTH3F, {axisPt, axisDCAxy, axisCent});
    registry.add("dcaVsPtPr", "Primary protons;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm);Centrality Perc.;", kTH3F, {axisPt, axisDCAxy, axisCent});

    auto hstat = registry.get<TH1>(HIST("EventCounter"));
    auto* x = hstat->GetXaxis();
    x->SetBinLabel(1, "All");
    x->SetBinLabel(2, "SelEigth");
    x->SetBinLabel(3, "NoSameBunchPileup");
    x->SetBinLabel(4, "GoodZvtxFT0vsPV");
    x->SetBinLabel(5, "NoCollInTimeRangeStrict");
    x->SetBinLabel(6, "NoCollInTimeRangeStandard");
    x->SetBinLabel(7, "NoCollInRofStrict");
    x->SetBinLabel(8, "NoCollInRofStandard");
    x->SetBinLabel(9, "NoHighMultCollInPrevRof");
    x->SetBinLabel(10, "NoCollInTimeRangeNarrow");
    x->SetBinLabel(11, "Occupancy Cut");
    x->SetBinLabel(12, "Has FT0?");
    x->SetBinLabel(13, "Cent. Sel.");
    x->SetBinLabel(14, "VtxZ Sel.");
    x->SetBinLabel(15, "Nch Sel.");
    x->SetBinLabel(16, "INEL > 0");

    if (doprocessCalibrationAndV0s) {
      registry.add("T0CcentVsRCTSel", "Bad RCT(=0.5) Good RCT(=1.5) Good RCT & Good PID RCT(=2.5);;RCT Status;", kTH2F, {{{axisCent}, {3, 0, 3}}});
      registry.add("NchVsNPV", ";Nch; NPV;", kTH2F, {{{nBinsNPV, minNpv, maxNpv}, {nBinsNch, minNch, maxNch}}});
      registry.add("ExcludedEvtVsNch", ";Nch;Entries;", kTH1F, {{nBinsNch, minNch, maxNch}});
      registry.add("ExcludedEvtVsNPV", ";NPV;Entries;", kTH1F, {{nBinsNPV, minNpv, maxNpv}});
      registry.add("TrackDaughterCounter", "itsrefit, and itshit NOT appplied for electrons sel.;Trk Sel.; Entries;", kTH1F, {{14, 0.5, 14.5}});
      registry.add("V0sCounter", ";V0 type; Entries;", kTH1F, {{4, 0.5, 4.5}});
      registry.add("dcaDauVsPt", ";V0 #it{p}_{T} (GeV/#it{c});DCA_{xy} (cm) daughters;", kTH2F, {axisPt, axisDCAxy});
      registry.add("nSigPiFromK0s", ";#it{n#sigma};;", kTH2F, {axisPtV0s, axisNsigmaTPC});
      registry.add("nSigPiFromL", ";#it{n#sigma};;", kTH2F, {axisPtV0s, axisNsigmaTPC});
      registry.add("nSigPrFromL", ";#it{n#sigma};;", kTH2F, {axisPtV0s, axisNsigmaTPC});
      registry.add("nSigPiFromAL", ";#it{n#sigma};;", kTH2F, {axisPtV0s, axisNsigmaTPC});
      registry.add("nSigPrFromAL", ";#it{n#sigma};;", kTH2F, {axisPtV0s, axisNsigmaTPC});
      registry.add("nSigElFromG", ";#it{n#sigma};;", kTH2F, {axisPtV0s, axisNsigmaTPC});
      registry.add("ArmAll", "Armenteros-Podolanski;#alpha;q_{T} (GeV/c)", kTH2F, {axisArmAlpha, axisArmqT});
      registry.add("ArmAfterTopoSel", "Armenteros-Podolanski anfter topological selection;#alpha;q_{T} (GeV/c)", kTH2F, {axisArmAlpha, axisArmqT});
      registry.add("ArmK0NOSel", "Armenteros-Podolanski WITH OUT 5 #times q_{T} > #alpha selection;#alpha;q_{T} (GeV/c)", kTH2F, {axisArmAlpha, axisArmqT});
      registry.add("ArmK0", "Armenteros-Podolanski WITH 5 #times q_{T} > #alpha selection;#alpha;q_{T} (GeV/c)", kTH2F, {axisArmAlpha, axisArmqT});
      registry.add("ArmL", "Armenteros-Podolanski;#alpha;q_{T} (GeV/c)", kTH2F, {axisArmAlpha, axisArmqT});
      registry.add("ArmAL", "Armenteros-Podolanski;#alpha;q_{T} (GeV/c)", kTH2F, {axisArmAlpha, axisArmqT});
      registry.add("ArmG", "Armenteros-Podolanski;#alpha;q_{T} (GeV/c)", kTH2F, {axisArmAlpha, axisArmqT});
      registry.add("MassK0sVsPt", ";;Inv. Mass (GeV/#it{c}^{2});", kTH2F, {axisPtV0s, axisK0Mass});
      registry.add("MassLVsPt", ";;Inv. Mass (GeV/#it{c}^{2});", kTH2F, {axisPtV0s, axisLambdaMass});
      registry.add("MassALVsPt", ";;Inv. Mass (GeV/#it{c}^{2});", kTH2F, {axisPtV0s, axisLambdaMass});
      registry.add("MassGVsPt", ";;Inv. Mass (GeV/#it{c}^{2});", kTH2F, {axisPtV0s, axisGammaMass});

      registry.add("NclVsPhipBeforeCut", Form("Found #LTNcl#GT TPC;%s (GeV/#it{c});#varphi", titlePorPt.data()), kTProfile2D, {{{axisXPhiCut}, {350, 0.0, 0.35}}});
      registry.add("NclVsPhipBeforeCutPID", Form("#LTNcl#GT used for PID;%s (GeV/#it{c});#varphi", titlePorPt.data()), kTProfile2D, {{{axisXPhiCut}, {350, 0.0, 0.35}}});
      registry.add("NclVsPhipAfterCut", Form("Found #LTNcl#GT TPC;%s (GeV/#it{c});#varphi", titlePorPt.data()), kTProfile2D, {{{axisXPhiCut}, {350, 0.0, 0.35}}});
      registry.add("NclVsPhipAfterCutPID", Form("#LTNcl#GT used for PID;%s (GeV/#it{c});#varphi", titlePorPt.data()), kTProfile2D, {{{axisXPhiCut}, {350, 0.0, 0.35}}});
      registry.add("NclVsEta", ";#eta;Found Ncl TPC", kTH2F, {{{axisEta}, {161, -0.5, 160.5}}});
      registry.add("NclVsEtap", ";#eta;Found #LTNcl#GT TPC", kTProfile, {axisEta});

      registry.add("NclVsEtaPiMIP", "MIP #pi^{+} + #pi^{-} (0.4 < #it{p} < 0.6 GeV/#it{c}, 40 < dE/dx < 60);#eta;Ncl TPC", kTH2F, {{{axisEta}, {161, -0.5, 160.5}}});
      registry.add("NclVsEtaPiMIPp", "MIP #pi^{+} + #pi^{-} (0.4 < #it{p} < 0.6 GeV/#it{c}, 40 < dE/dx < 60);#eta;#LTNcl#GT TPC", kTProfile, {axisEta});
      registry.add("NclVsEtaPiV0", ";#eta;Ncl TPC", kTH2F, {axisEta, axisNcl});
      registry.add("NclVsEtaPiV0p", ";#eta;#LTNcl#GT TPC", kTProfile, {axisEta});
      registry.add("NclVsEtaPrV0", ";#eta;Ncl TPC", kTH2F, {axisEta, axisNcl});
      registry.add("NclVsEtaPrV0p", ";#eta;#LTNcl#GT TPC", kTProfile, {axisEta});
      registry.add("NclVsEtaElV0", ";#eta;Ncl TPC", kTH2F, {axisEta, axisNcl});
      registry.add("NclVsEtaElV0p", ";#eta;#LTNcl#GT TPC", kTProfile, {axisEta});

      registry.add("EtaVsPhi", ";#eta;#varphi;", kTH2F, {{axisEta}, {100, 0, o2::constants::math::TwoPI}});
      registry.add("EtaVsYK0s", ";#eta;#it{y};", kTH2F, {axisEta, axisY});
      registry.add("EtaVsYPiL", ";#eta;#it{y};", kTH2F, {axisEta, axisY});
      registry.add("EtaVsYPrL", ";#eta;#it{y};", kTH2F, {axisEta, axisY});
      registry.add("EtaVsYPiAL", ";#eta;#it{y};", kTH2F, {axisEta, axisY});
      registry.add("EtaVsYPrAL", ";#eta;#it{y};", kTH2F, {axisEta, axisY});
      registry.add("EtaVsYG", ";#eta;#it{y};", kTH2F, {axisEta, axisY});

      // registry.add("TOFExpPi2TOF", ";Momentum (GeV/#it{c});t^{#pi}_{Exp}/t_{TOF}", kTH2F, {{{axisPtV0s}, {100, 0.2, 1.2}}});
      registry.add("DCAxyPtPiK0s", ";DCA_{xy} (cm); p_{T} (GeV/c)", kTH2F, {axisDCAxy, axisPt});
      registry.add("DCAxyPtPrL", ";DCA_{xy} (cm); p_{T} (GeV/c)", kTH2F, {axisDCAxy, axisPt});
      registry.add("DCAxyPtPrAL", ";DCA_{xy} (cm); p_{T} (GeV/c)", kTH2F, {axisDCAxy, axisPt});

      // registry.add("betaVsMomentum", ";Momentum (GeV/#it{c}); #beta", kTH2F, {{{axisPtV0s}, {500, 0, 1.2}}});
      registry.add("dEdxVsEtaPiMIP", "MIP #pi^{+} + #pi^{-} (0.4 < #it{p} < 0.6 GeV/#it{c});#eta; dE/dx;", kTH2F, {{{axisEta}, {100, 0, 100}}});
      registry.add("dEdxVsEtaPiMIPp", "MIP #pi^{+} + #pi^{-} (0.4 < #it{p} < 0.6 GeV/#it{c});#eta; #LTdE/dx#GT", kTProfile, {axisEta});
      registry.add("dEdxVsEtaElMIP", "MIP e^{+} + e^{-} (0.3 < #it{p} < 0.45 GeV/#it{c});#eta; dE/dx;", kTH2F, {{{axisEta}, {100, 0, 100}}});
      registry.add("dEdxVsEtaElMIPp", "MIP e^{+} + e^{-} (0.3 < #it{p} < 0.45 GeV/#it{c});#eta; #LTdE/dx#GT", kTProfile, {axisEta});
      registry.add("dEdxVsEtaPiMIPV0", "MIP #pi^{+} + #pi^{-} (0.4 < #it{p} < 0.6 GeV/#it{c});#eta; dE/dx", kTH2F, {{{axisEta}, {100, 0, 100}}});
      registry.add("dEdxVsEtaPiMIPV0p", "MIP #pi^{+} + #pi^{-} (0.4 < #it{p} < 0.6 GeV/#it{c});#eta; #LTdE/dx#GT", kTProfile, {axisEta});
      registry.add("dEdxVsEtaElMIPV0", "e^{+} + e^{-} (0.15 <#it{p}_{T} < 50 GeV/#it{c});#eta; dE/dx", kTH2F, {{{axisEta}, {100, 0, 100}}});
      registry.add("dEdxVsEtaElMIPV0p", "e^{+} + e^{-} (0.15 <#it{p}_{T} < 50 GeV/#it{c});#eta; #LTdE/dx#GT", kTProfile, {axisEta});

      registry.add("pTVsCent", "", kTH2F, {axisPt, axisCent});

      for (int i = 0; i < kNEtaHists; ++i) {
        dEdx[i] = registry.add<TH3>(Form("dEdx_%s", endingEta[i]), Form("%s;Momentum (GeV/#it{c});dE/dx;", latexEta[i]), kTH3F, {axisPt, axisdEdx, axisCent});
        pTVsP[i] = registry.add<TH2>(Form("pTVsP_%s", endingEta[i]), Form("%s;Momentum (GeV/#it{c});#it{p}_{T} (GeV/#it{c});", latexEta[i]), kTH2F, {axisPt, axisPt});
        dEdxPiV0[i] = registry.add<TH3>(Form("dEdxPiV0_%s", endingEta[i]), Form("#pi^{+} + #pi^{-}, %s;Momentum (GeV/#it{c});dE/dx;", latexEta[i]), kTH3F, {axisPtV0s, axisdEdx, axisCent});
        dEdxPrV0[i] = registry.add<TH3>(Form("dEdxPrV0_%s", endingEta[i]), Form("p + #bar{p}, %s;Momentum (GeV/#it{c});dE/dx;", latexEta[i]), kTH3F, {axisPtV0s, axisdEdx, axisCent});
        dEdxElV0[i] = registry.add<TH3>(Form("dEdxElV0_%s", endingEta[i]), Form("e^{+} + e^{-}, %s;Momentum (GeV/#it{c});dE/dx;", latexEta[i]), kTH3F, {axisPtV0s, axisdEdx, axisCent});
        dEdxPiTOF[i] = registry.add<TH3>(Form("dEdxPiTOF_%s", endingEta[i]), Form("#pi^{+} + #pi^{-}, %s;Momentum (GeV/#it{c});dE/dx;", latexEta[i]), kTH3F, {axisPtV0s, axisdEdx, axisCent});
        nClVsdEdxPiV0[i] = registry.add<TH2>(Form("NclVsdEdxPiV0_%s", endingEta[i]), Form("%s;#it{N}_{cl} used for PID;dE/dx;", latexEta[i]), kTH2F, {axisNcl, axisdEdx});
        nClVsdEdxElV0[i] = registry.add<TH2>(Form("NclVsdEdxElV0_%s", endingEta[i]), Form("%s;#it{N}_{cl} used for PID;dE/dx;", latexEta[i]), kTH2F, {axisNcl, axisdEdx});
        nClVsdEdxPrV0[i] = registry.add<TH2>(Form("NclVsdEdxPrV0_%s", endingEta[i]), Form("%s;#it{N}_{cl} used for PID;dE/dx;", latexEta[i]), kTH2F, {axisNcl, axisdEdx});
        nClVsP[i] = registry.add<TH2>(Form("NclVsPPrimaries_%s", endingEta[i]), Form("%s;;Ncl TPC", latexEta[i]), kTH2F, {axisPtNcl, axisNcl});
        nClVsPElV0[i] = registry.add<TH2>(Form("NclVsPElV0_%s", endingEta[i]), Form("%s;;Ncl TPC", latexEta[i]), kTH2F, {axisPtNcl, axisNcl});
        nClVsPPiV0[i] = registry.add<TH2>(Form("NclVsPPiV0_%s", endingEta[i]), Form("%s;;Ncl TPC", latexEta[i]), kTH2F, {axisPtNcl, axisNcl});
        nClVsPPrV0[i] = registry.add<TH2>(Form("NclVsPPrV0_%s", endingEta[i]), Form("%s;;Ncl TPC", latexEta[i]), kTH2F, {axisPtNcl, axisNcl});
        nClVsPp[i] = registry.add<TProfile>(Form("NclVsPrimariesp_%s", endingEta[i]), Form("%s;;#LT#it{N}_{cl}#GT TPC", latexEta[i]), kTProfile, {axisPtNcl});
        nClVsPpElV0[i] = registry.add<TProfile>(Form("NclVsPElV0p_%s", endingEta[i]), Form("%s;;#LT#it{N}_{cl}#GT TPC", latexEta[i]), kTProfile, {axisPtNcl});
        nClVsPpPiV0[i] = registry.add<TProfile>(Form("NclVsPPiV0p_%s", endingEta[i]), Form("%s;;#LT#it{N}_{cl}#GT TPC", latexEta[i]), kTProfile, {axisPtNcl});
        nClVsPpPrV0[i] = registry.add<TProfile>(Form("NclVsPPrV0p_%s", endingEta[i]), Form("%s;;#LT#it{N}_{cl}#GT TPC", latexEta[i]), kTProfile, {axisPtNcl});
        nClVsdEdxpElV0[i] = registry.add<TProfile>(Form("NclVsdEdxElV0p_%s", endingEta[i]), Form("%s;;#LTd#it{E}/d#it{x}#GT", latexEta[i]), kTProfile, {axisNcl});
        nClVsdEdxpPiV0[i] = registry.add<TProfile>(Form("NclVsdEdxPiV0p_%s", endingEta[i]), Form("%s;;#LTd#it{E}/d#it{x}#GT", latexEta[i]), kTProfile, {axisNcl});
        nClVsdEdxpPrV0[i] = registry.add<TProfile>(Form("NclVsdEdxPrV0p_%s", endingEta[i]), Form("%s;;#LTd#it{E}/d#it{x}#GT", latexEta[i]), kTProfile, {axisNcl});
      }

      auto htrkSel = registry.get<TH1>(HIST("TrackDaughterCounter"));
      auto* xtrkSel = htrkSel->GetXaxis();
      xtrkSel->SetBinLabel(1, "All");
      xtrkSel->SetBinLabel(2, "Eta");
      xtrkSel->SetBinLabel(3, "Pt");
      xtrkSel->SetBinLabel(4, "XRows");
      xtrkSel->SetBinLabel(5, "XRowsOverFindableCls");
      xtrkSel->SetBinLabel(6, "Chi2TPC");
      xtrkSel->SetBinLabel(7, "Chi2ITS");
      xtrkSel->SetBinLabel(8, "Itsrefit");
      xtrkSel->SetBinLabel(9, "Tpcrefit");
      xtrkSel->SetBinLabel(10, "Golden");
      xtrkSel->SetBinLabel(11, "Itshit");
      xtrkSel->SetBinLabel(12, "Passed all");
    }

    if (doprocessMC || doprocessSim) {
      registry.add("zPosMC", "Generated Events With at least One Rec. Collision + Sel. criteria;;Entries;", kTH1F, {axisZpos});
      registry.add("dcaVsPtPiDec", "Secondary pions from decays;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm);Centrality Perc.;", kTH3F, {axisPt, axisDCAxy, axisCent});
      registry.add("dcaVsPtPrDec", "Secondary protons from decays;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm);Centrality Perc.;", kTH3F, {axisPt, axisDCAxy, axisCent});
      registry.add("dcaVsPtPiMat", "Secondary pions from material interactions;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm);Centrality Perc.;", kTH3F, {axisPt, axisDCAxy, axisCent});
      registry.add("dcaVsPtPrMat", "Secondary protons from material interactions;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm);Centrality Perc.;", kTH3F, {axisPt, axisDCAxy, axisCent});
      registry.add("NclVsPhip", Form("#LTNcl#GT used for PID;%s (GeV/#it{c});#varphi", titlePorPt.data()), kTProfile2D, {{{axisXPhiCut}, {350, 0.0, 0.35}}});
    }

    if (doprocessMC) {

      registry.add("EventCounterMC", "Event counter", kTH1F, {axisEvent});

      registry.add("PtPiVsCent", "", kTH2F, {axisPt, axisCent});
      registry.add("PtKaVsCent", "", kTH2F, {axisPt, axisCent});
      registry.add("PtPrVsCent", "", kTH2F, {axisPt, axisCent});

      registry.add("PtPiVsCentMC", "", kTH2F, {axisPt, axisCent});
      registry.add("PtKaVsCentMC", "", kTH2F, {axisPt, axisCent});
      registry.add("PtPrVsCentMC", "", kTH2F, {axisPt, axisCent});
    }

    if (doprocessSim) {

      // MC events passing the TVX requirement
      registry.add("NchMCcentVsTVX", ";Passed(=1.5) NOT Passed(=0.5);", kTH2F, {{{nBinsNch, minNch, maxNch}, {2, 0, 2}}});

      registry.add("NumberOfRecoCollisions", "Number of times Gen. Coll.are reconstructed;N;Entries", kTH1F, {{10, -0.5, 9.5}});

      // Pt resolution
      registry.add("PtResolution", "p_{T} resolution;;(pt_{rec} - pt_{gen})/pt_{gen};", kTH2F, {axisPt, {100, -1.0, 1.0}});

      // Needed to calculate the numerator of the Acceptance X Efficiency
      registry.add("PtPiVsCent_WithRecoEvt", "Generated Events With at least One Rec. Collision + Sel. criteria;;;", kTH2F, {axisPt, axisCent});
      registry.add("PtKaVsCent_WithRecoEvt", "Generated Events With at least One Rec. Collision + Sel. criteria;;;", kTH2F, {axisPt, axisCent});
      registry.add("PtPrVsCent_WithRecoEvt", "Generated Events With at least One Rec. Collision + Sel. criteria;;;", kTH2F, {axisPt, axisCent});

      // Needed to calculate the denominator of the Acceptance X Efficiency
      registry.add("PtPiVsCentMC_WithRecoEvt", "Generated Events With at least One Rec. Collision;;;", kTH2F, {axisPt, axisCent});
      registry.add("PtKaVsCentMC_WithRecoEvt", "Generated Events With at least One Rec. Collision;;;", kTH2F, {axisPt, axisCent});
      registry.add("PtPrVsCentMC_WithRecoEvt", "Generated Events With at least One Rec. Collision;;;", kTH2F, {axisPt, axisCent});

      // Needed for the Gen. Nch to Centrality conversion
      registry.add("NchMCVsCent", "Generated Nch v.s. Centrality (At least Once Rec. Coll. + Sel. criteria);;Gen. Nch", kTH2F, {{axisCent, {nBinsNch, minNch, maxNch}}});

      // Needed to measure Event Loss
      registry.add("NchMC_WithRecoEvt", "Generated Nch of Evts With at least one Rec. Coll. + Sel. criteria;Gen. Nch MC;Entries", kTH1F, {{nBinsNch, minNch, maxNch}});
      registry.add("NchMC_AllGen", "Generated Nch of All Gen. Evts.;Gen. Nch;Entries", kTH1F, {{nBinsNch, minNch, maxNch}});

      // Needed to measure Event Splitting
      registry.add("Centrality_WRecoEvt", "Generated Events With at least One Rec. Collision And NO Sel. criteria;;Entries", kTH1F, {axisCent});
      registry.add("Centrality_WRecoEvtWSelCri", "Generated Events With at least One Rec. Collision + Sel. criteria;;Entries", kTH1F, {axisCent});
      registry.add("Centrality_AllRecoEvt", "Generated Events Irrespective of the number of times it was reconstructed + Evt. Selections;;Entries", kTH1F, {axisCent});

      // Needed to calculate the numerator of the Signal Loss correction
      registry.add("PtPiVsNchMC_WithRecoEvt", "Generated Events With at least One Rec. Collision;;Gen. Nch;", kTH2F, {{axisPt, {nBinsNch, minNch, maxNch}}});
      registry.add("PtKaVsNchMC_WithRecoEvt", "Generated Events With at least One Rec. Collision;;Gen. Nch;", kTH2F, {{axisPt, {nBinsNch, minNch, maxNch}}});
      registry.add("PtPrVsNchMC_WithRecoEvt", "Generated Events With at least One Rec. Collision;;Gen. Nch;", kTH2F, {{axisPt, {nBinsNch, minNch, maxNch}}});

      // Needed to calculate the denominator of the Signal Loss correction
      registry.add("PtPiVsNchMC_AllGen", "All Generated Events;;Gen. Nch;", kTH2F, {{axisPt, {nBinsNch, minNch, maxNch}}});
      registry.add("PtKaVsNchMC_AllGen", "All Generated Events;;Gen. Nch;", kTH2F, {{axisPt, {nBinsNch, minNch, maxNch}}});
      registry.add("PtPrVsNchMC_AllGen", "All Generated Events;;Gen. Nch;", kTH2F, {{axisPt, {nBinsNch, minNch, maxNch}}});

      registry.add("MCclosure_PtMCPiVsNchMC", "All Generated Events 4 MC closure;;Gen. Nch;", kTH2F, {{axisPt, {nBinsNch, minNch, maxNch}}});
      registry.add("MCclosure_PtMCKaVsNchMC", "All Generated Events 4 MC closure;;Gen. Nch;", kTH2F, {{axisPt, {nBinsNch, minNch, maxNch}}});
      registry.add("MCclosure_PtMCPrVsNchMC", "All Generated Events 4 MC closure;;Gen. Nch;", kTH2F, {{axisPt, {nBinsNch, minNch, maxNch}}});

      registry.add("MCclosure_PtPiVsNchMC", "Gen Evts With at least one Rec. Coll. + Sel. criteria 4 MC closure;;Gen. Nch;", kTH2F, {{axisPt, {nBinsNch, minNch, maxNch}}});
      registry.add("MCclosure_PtKaVsNchMC", "Gen Evts With at least one Rec. Coll. + Sel. criteria 4 MC closure;;Gen. Nch;", kTH2F, {{axisPt, {nBinsNch, minNch, maxNch}}});
      registry.add("MCclosure_PtPrVsNchMC", "Gen Evts With at least one Rec. Coll. + Sel. criteria 4 MC closure;;Gen. Nch;", kTH2F, {{axisPt, {nBinsNch, minNch, maxNch}}});
    }

    LOG(info) << "\trequireGoodRct=" << requireGoodRct.value;
    LOG(info) << "\tccdbNoLaterThan=" << ccdbNoLaterThan.value;
    LOG(info) << "\tapplyNchSel=" << applyNchSel.value;
    LOG(info) << "\tselINELgt0=" << selINELgt0.value;
    LOG(info) << "\tv0TypeSelection=" << static_cast<int>(v0Selections.v0TypeSelection);
    LOG(info) << "\tselElecFromGammas=" << v0Selections.selElecFromGammas;
    LOG(info) << "\tapplyInvMassSel=" << v0Selections.applyInvMassSel;
    LOG(info) << "\trequireITShit=" << v0Selections.requireITShit;
    LOG(info) << "\tminPt=" << v0Selections.minPt;
    LOG(info) << "\tmaxPt=" << v0Selections.maxPt;
    LOG(info) << "\tminPtDaughter=" << v0Selections.minPtDaughter;
    LOG(info) << "\tmaxPtDaughter=" << v0Selections.maxPtDaughter;
    LOG(info) << "\tuseNclsPID=" << v0Selections.useNclsPID;
    LOG(info) << "\tqTSel=" << v0Selections.qTSel;
    LOG(info) << "\tarmAlphaSel=" << v0Selections.armAlphaSel;
    LOG(info) << "\tapplyNclSel=" << v0Selections.applyNclSel;
    LOG(info) << "\tapplyPhiCut=" << v0Selections.applyPhiCut;
    LOG(info) << "\tusePinPhiSelection=" << v0Selections.usePinPhiSelection;
    LOG(info) << "\ttitlePorPt=" << titlePorPt;
    LOG(info) << "\tcurrentRunNumberNchSel=" << currentRunNumberNchSel;
    LOG(info) << "\tcurrentRunNumberPhiSel=" << currentRunNumberPhiSel;

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    ccdb->setCreatedNotAfter(ccdbNoLaterThan.value);

    if (applyNchSel.value) {
      LOG(info) << "\tLoading Nch-based selections!";
      LOG(info) << "\tpathMeanNch=" << pathMeanNch.value;
      LOG(info) << "\tpathSigmaNch=" << pathSigmaNch.value;
    }

    if (v0Selections.applyPhiCut) {
      LOG(info) << "\tLoading Phi cut!";
      LOG(info) << "\tpathPhiCutLow=" << pathPhiCutLow.value;
      LOG(info) << "\tpathPhiCutHigh=" << pathPhiCutHigh.value;
    }

    if (v0Selections.applyEtaCal) {
      LOG(info) << "\tLoading Eta Cal!";
      LOG(info) << "\tpathEtaCal=" << pathEtaCal.value;
      loadEtaCalibration();
      LOG(info) << "\tisMIPCalLoaded=" << etaCal.isMIPCalLoaded;
    }

    if (v0Selections.applyPlateauSel) {
      LOG(info) << "\tLoading Eta Plateau Cal!";
      LOG(info) << "\tpathEtaCalPlateau=" << pathEtaCalPlateau.value;
      loadEtaPlateauCalibration();
      LOG(info) << "\tisCalPlateauLoaded=" << etaCal.isCalPlateauLoaded;
    }

    if (v0Selections.applyNclSel)
      LOG(info) << "\t minNcl=" << v0Selections.minNcl;
  }

  void processCalibrationAndV0s(ColEvSels::iterator const& collision, BCsRun3 const& /**/, aod::V0Datas const& v0s, aod::FV0As const& /**/, aod::FT0s const& /**/, TracksFull const& tracks)
  {
    // LOG(info) << " Collisions size: " << collisions.size() << "
    // Table's size: " << collisions.tableSize() << "\n";
    // LOG(info) << "Run number: " << foundBC.runNumber() << "\n";

    if (!isEventSelected(collision)) {
      return;
    }

    const auto& foundBC = collision.foundBC_as<BCsRun3>();
    const uint64_t timeStamp{foundBC.timestamp()};
    const int magField{getMagneticField(timeStamp)};
    const double nPV{collision.multNTracksPVeta1() / 1.};
    const float centrality{isT0Ccent ? collision.centFT0C() : collision.centFT0M()};

    //---------------------------
    // Control histogram
    //---------------------------
    if (selHasFT0 && !collision.has_foundFT0()) {
      registry.fill(HIST("T0CcentVsFoundFT0"), centrality, 0.5);
    }
    registry.fill(HIST("T0CcentVsFoundFT0"), centrality, 1.5);

    // Apply RCT selection?
    if (requireGoodRct) {

      // Checks if collisions passes RCT selection
      if (!rctChecker(*collision)) {
        registry.fill(HIST("T0CcentVsRCTSel"), centrality, 0.5);
        return;
      }

      registry.fill(HIST("T0CcentVsRCTSel"), centrality, 1.5);
      // Checks if collisions passes good PID RCT status
      if (requireGoodPIDRct && collision.rct_bit(kTPCBadPID)) {
        return;
      }
      registry.fill(HIST("T0CcentVsRCTSel"), centrality, 2.5);
    }

    if (applyNchSel) {
      const int nextRunNumber{foundBC.runNumber()};
      if (currentRunNumberNchSel != nextRunNumber) {
        loadNchCalibrations(timeStamp);
        currentRunNumberNchSel = nextRunNumber;
        LOG(info) << "\tcurrentRunNumberNchSel= " << currentRunNumberNchSel << " timeStamp = " << timeStamp;
      }

      // return if Nch selection objects are nullptr
      if (!(cfgNch.hMeanNch && cfgNch.hSigmaNch))
        return;
    }

    int nch{0};
    for (const auto& track : tracks) {
      // Track Selection
      if (!trkSelGlobal.IsSelected(track)) {
        continue;
      }
      if (track.pt() < kMinPtNchSel || track.pt() > kMaxPtNchSel) {
        continue;
      }
      nch++;
    }

    bool skipEvent{false};
    if (applyNchSel) {
      if (!cfgNch.calibrationsLoaded)
        return;

      const double xEval{nPV};
      const int bin4Calibration{cfgNch.hMeanNch->FindBin(xEval)};
      const double meanNch{cfgNch.hMeanNch->GetBinContent(bin4Calibration)};
      const double sigmaNch{cfgNch.hSigmaNch->GetBinContent(bin4Calibration)};
      const double nSigmaSelection{nSigmaNchCut * sigmaNch};
      const double diffMeanNch{meanNch - nch};

      if (std::abs(diffMeanNch) > nSigmaSelection) {
        registry.fill(HIST("ExcludedEvtVsNch"), nch);
        registry.fill(HIST("ExcludedEvtVsNPV"), nPV);
        skipEvent = true;
      }
    }

    if (applyNchSel && skipEvent) {
      return;
    }

    if (applyNchSel)
      registry.fill(HIST("EventCounter"), EvCutLabel::NchSel);

    registry.fill(HIST("NchVsNPV"), nPV, nch);
    registry.fill(HIST("zPos"), collision.posZ());
    registry.fill(HIST("NchVsCent"), centrality, nch);

    if (v0Selections.applyPhiCut) {
      const int nextRunNumber{foundBC.runNumber()};
      if (currentRunNumberPhiSel != nextRunNumber) {
        loadPhiCutSelections(timeStamp);
        currentRunNumberPhiSel = nextRunNumber;
        LOG(info) << "\tcurrentRunNumberPhiSel= " << currentRunNumberPhiSel << " timeStamp = " << timeStamp;
      }

      // return if phi cut objects are nullptr
      if (!(phiCut.hPhiCutHigh && phiCut.hPhiCutLow))
        return;
    }

    registry.fill(HIST("T0Ccent"), centrality);
    // Fill DCAxy vs pT for secondary-particle contamination correction
    for (const auto& track : tracks) {
      // Track Selection
      if (!trkSelGlobalOpenDCAxy.IsSelected(track))
        continue;

      if (track.pt() < v0Selections.minPt || track.pt() > v0Selections.maxPt)
        continue;

      if (track.eta() < v0Selections.minEtaDaughter || track.eta() > v0Selections.maxEtaDaughter)
        continue;

      const float momentum{track.p()};
      const float pt{track.pt()};
      const float phi{track.phi()};
      const float pOrPt{v0Selections.usePinPhiSelection ? momentum : pt};
      const int16_t nclFound{track.tpcNClsFound()};
      const int16_t nclPID{track.tpcNClsPID()};
      const int16_t ncl{v0Selections.useNclsPID ? nclPID : nclFound};

      if (v0Selections.applyNclSel && ncl < v0Selections.minNcl)
        continue;

      float phiPrime{phi};
      const int charge{track.sign()};
      phiPrimeFunc(phiPrime, magField, charge);

      if (v0Selections.applyPhiCut) {
        if (!passesPhiSelection(pOrPt, phiPrime))
          continue;
      }

      const float piTPCNsigma{std::fabs(track.tpcNSigmaPi())};
      const float prTPCNsigma{std::fabs(track.tpcNSigmaPr())};
      const float piTOFNsigma{std::fabs(track.tofNSigmaPi())};
      const float prTOFNsigma{std::fabs(track.tofNSigmaPr())};
      const double piRadiusNsigma{std::sqrt(std::pow(piTPCNsigma, 2.) + std::pow(piTOFNsigma, 2.))};
      const double prRadiusNsigma{std::sqrt(std::pow(prTPCNsigma, 2.) + std::pow(prTOFNsigma, 2.))};

      if (piRadiusNsigma < kThree)
        registry.fill(HIST("dcaVsPtPi"), track.pt(), track.dcaXY(), centrality);
      if (prRadiusNsigma < kThree)
        registry.fill(HIST("dcaVsPtPr"), track.pt(), track.dcaXY(), centrality);
    }

    for (const auto& track : tracks) {

      if (!trkSelGlobal.IsSelected(track))
        continue;

      if (track.pt() < v0Selections.minPt || track.pt() > v0Selections.maxPt)
        continue;

      if (track.eta() < v0Selections.minEtaDaughter || track.eta() > v0Selections.maxEtaDaughter)
        continue;

      const float momentum{track.p()};
      const float pt{track.pt()};
      const float phi{track.phi()};
      const float eta{track.eta()};
      float dedx{track.tpcSignal()};
      const int charge{track.sign()};
      const float pOrPt{v0Selections.usePinPhiSelection ? momentum : pt};
      const int16_t nclFound{track.tpcNClsFound()};
      const int16_t nclPID{track.tpcNClsPID()};
      const int16_t ncl{v0Selections.useNclsPID ? nclPID : nclFound};

      if (v0Selections.applyNclSel && ncl < v0Selections.minNcl)
        continue;

      float phiPrime{phi};
      phiPrimeFunc(phiPrime, magField, charge);
      registry.fill(HIST("NclVsPhipBeforeCut"), pOrPt, phiPrime, nclFound);
      registry.fill(HIST("NclVsPhipBeforeCutPID"), pOrPt, phiPrime, nclPID);

      if (v0Selections.applyPhiCut) {
        if (!passesPhiSelection(pOrPt, phiPrime))
          continue;
      }

      if (v0Selections.applyEtaCal && etaCal.isMIPCalLoaded) {
        const double dedxCal{etaCal.pEtaCal->GetBinContent(etaCal.pEtaCal->FindBin(eta))};
        if (dedxCal > kMindEdxMIP && dedxCal < kMaxdEdxMIP)
          dedx *= (50.0 / dedxCal);
        else
          continue;
      }

      int indexEta{-999};
      for (int i = 0; i < kNEtaHists; ++i) {
        if (eta >= kLowEta[i] && eta < kHighEta[i]) {
          indexEta = i;
          break;
        }
      }

      if (indexEta < kZeroInt || indexEta > kSevenInt)
        continue;

      if (momentum > kMinPMIP && momentum < kMaxPMIP && dedx > kMindEdxMIP && dedx < kMaxdEdxMIP) {
        registry.fill(HIST("dEdxVsEtaPiMIP"), eta, dedx);
        registry.fill(HIST("dEdxVsEtaPiMIPp"), eta, dedx);
        registry.fill(HIST("NclVsEtaPiMIP"), eta, ncl);
        registry.fill(HIST("NclVsEtaPiMIPp"), eta, ncl);
      }

      if (momentum > kMinPElMIP && momentum < kMaxPElMIP && dedx > kMindEdxMIPPlateau && dedx < kMaxdEdxMIPPlateau) {
        registry.fill(HIST("dEdxVsEtaElMIP"), eta, dedx);
        registry.fill(HIST("dEdxVsEtaElMIPp"), eta, dedx);
      }

      dEdx[indexEta]->Fill(momentum, dedx, centrality);
      pTVsP[indexEta]->Fill(momentum, pt);
      nClVsP[indexEta]->Fill(pOrPt, ncl);
      nClVsPp[indexEta]->Fill(pOrPt, ncl);
      registry.fill(HIST("pTVsCent"), pt, centrality);
      registry.fill(HIST("EtaVsPhi"), eta, track.phi());
      registry.fill(HIST("NclVsEta"), eta, nclFound);
      registry.fill(HIST("NclVsEtap"), eta, nclFound);
      registry.fill(HIST("NclVsEtaPID"), eta, nclPID);
      registry.fill(HIST("NclVsEtaPIDp"), eta, nclPID);
      registry.fill(HIST("NclVsPhipAfterCut"), pOrPt, phiPrime, nclFound);
      registry.fill(HIST("NclVsPhipAfterCutPID"), pOrPt, phiPrime, nclPID);

      if (track.hasTOF() && track.goodTOFMatch()) {
        const float tTOF{track.tofSignal()};
        const float trkLength{track.length()};
        const float tExpPiTOF{track.tofExpSignalPi(tTOF)};
        // const float tExpElTOF{track.tofExpSignalEl(tTOF)};

        if (trkLength > kZero && tTOF > kZero) {
          // registry.fill(HIST("betaVsMomentum"), momentum, track.beta());
          // registry.fill(HIST("TOFExpPi2TOF"), momentum, tExpPiTOF / tTOF);
          // registry.fill(HIST("TOFExpEl2TOF"), momentum, tExpElTOF / tTOF);

          if (std::abs((tExpPiTOF / tTOF) - kOne) < v0Selections.maxExpTOFPi) {
            dEdxPiTOF[indexEta]->Fill(momentum, dedx, centrality);
          }
        }
      }
    }

    for (const auto& v0 : v0s) {

      // Select V0 type
      if (v0.v0Type() != v0Selections.v0TypeSelection)
        continue;

      // Positive-(negative-)charged tracks (daughters)
      const auto& posTrack = v0.posTrack_as<TracksFull>();
      const auto& negTrack = v0.negTrack_as<TracksFull>();
      const int posTrackCharge{posTrack.sign()};
      const int negTrackCharge{negTrack.sign()};
      const float posTrkP{posTrack.p()};
      const float negTrkP{negTrack.p()};
      const float posTrkPt{posTrack.pt()};
      const float negTrkPt{negTrack.pt()};
      const float posTrkEta{posTrack.eta()};
      const float negTrkEta{negTrack.eta()};
      float posTrkdEdx{posTrack.tpcSignal()};
      float negTrkdEdx{negTrack.tpcSignal()};
      float posTrackPhiPrime{posTrack.phi()};
      float negTrackPhiPrime{negTrack.phi()};
      const int16_t posNclFound{posTrack.tpcNClsFound()};
      const int16_t negNclFound{negTrack.tpcNClsFound()};
      const int16_t posNclPID{posTrack.tpcNClsPID()};
      const int16_t negNclPID{negTrack.tpcNClsPID()};
      const int16_t posNcl = v0Selections.useNclsPID ? posNclPID : posNclFound;
      const int16_t negNcl = v0Selections.useNclsPID ? negNclPID : negNclFound;

      phiPrimeFunc(posTrackPhiPrime, magField, posTrackCharge);
      phiPrimeFunc(negTrackPhiPrime, magField, negTrackCharge);
      const float posPorPt{v0Selections.usePinPhiSelection ? posTrkP : posTrkPt};
      const float negPorPt{v0Selections.usePinPhiSelection ? negTrkP : negTrkPt};

      // Skip v0s with like-sig daughters
      if (posTrack.sign() == negTrack.sign())
        continue;

      // Passes Geometrical (Phi) cut?
      if (v0Selections.applyPhiCut) {
        if (!(passesPhiSelection(posPorPt, posTrackPhiPrime) && passesPhiSelection(negPorPt, negTrackPhiPrime)))
          continue;
      }

      // Passes daughters track-selection?
      if (!(passesTrackSelectionDaughters(posTrack) && passesTrackSelectionDaughters(negTrack)))
        continue;

      // Selection based on Ncl for PID
      if (v0Selections.applyNclSel && !(posNcl >= v0Selections.minNcl && negNcl >= v0Selections.minNcl))
        continue;

      // Eta calibration positive-charge track
      if (v0Selections.applyEtaCal && etaCal.isMIPCalLoaded) {
        const double dedxCal{etaCal.pEtaCal->GetBinContent(etaCal.pEtaCal->FindBin(posTrkEta))};
        if (dedxCal > kMindEdxMIP && dedxCal < kMaxdEdxMIP)
          posTrkdEdx *= (50.0 / dedxCal);
        else
          continue;
      }

      // Eta calibration negative-charge track
      if (v0Selections.applyEtaCal && etaCal.isMIPCalLoaded) {
        const double dedxCal{etaCal.pEtaCal->GetBinContent(etaCal.pEtaCal->FindBin(negTrkEta))};
        if (dedxCal > kMindEdxMIP && dedxCal < kMaxdEdxMIP)
          negTrkdEdx *= (50.0 / dedxCal);
        else
          continue;
      }

      const TVector3 ppos(posTrack.px(), posTrack.py(), posTrack.pz());
      const TVector3 pneg(negTrack.px(), negTrack.py(), negTrack.pz());
      double alpha, qT;

      getArmeterosVariables(ppos, pneg, alpha, qT);
      registry.fill(HIST("ArmAll"), alpha, qT);

      bool passesTopoSel{false};
      // Passes V0 topological cuts?
      if (passesV0TopologicalSelection(v0))
        passesTopoSel = true;

      if (!passesTopoSel)
        continue;

      int posIndexEta{-999};
      int negIndexEta{-999};
      for (int i = 0; i < kNEtaHists; ++i) {
        if (posTrkEta >= kLowEta[i] && posTrkEta < kHighEta[i]) {
          posIndexEta = i;
          break;
        }
      }

      for (int i = 0; i < kNEtaHists; ++i) {
        if (negTrkEta >= kLowEta[i] && negTrkEta < kHighEta[i]) {
          negIndexEta = i;
          break;
        }
      }

      if (posIndexEta < kZeroInt || posIndexEta > kSevenInt)
        continue;

      if (negIndexEta < kZeroInt || negIndexEta > kSevenInt)
        continue;

      registry.fill(HIST("ArmAfterTopoSel"), alpha, qT);
      registry.fill(HIST("dcaDauVsPt"), v0.pt(), v0.dcapostopv());
      registry.fill(HIST("dcaDauVsPt"), v0.pt(), v0.dcanegtopv());

      if (passesK0Selection(collision, v0)) {
        registry.fill(HIST("ArmK0NOSel"), alpha, qT);
        if (v0Selections.armPodCut * qT > std::abs(alpha)) { // Armenters selection
          registry.fill(HIST("V0sCounter"), V0sCounter::K0s);
          registry.fill(HIST("ArmK0"), alpha, qT);
          registry.fill(HIST("MassK0sVsPt"), v0.pt(), v0.mK0Short());
          registry.fill(HIST("nSigPiFromK0s"), posTrkPt, posTrack.tpcNSigmaPi());
          registry.fill(HIST("nSigPiFromK0s"), negTrkPt, negTrack.tpcNSigmaPi());
          registry.fill(HIST("NclVsEtaPiV0"), posTrkEta, posNcl);
          registry.fill(HIST("NclVsEtaPiV0p"), posTrkEta, posNcl);
          registry.fill(HIST("NclVsEtaPiV0"), negTrkEta, negNcl);
          registry.fill(HIST("NclVsEtaPiV0p"), negTrkEta, negNcl);
          nClVsPPiV0[posIndexEta]->Fill(posPorPt, posNcl);
          nClVsPpPiV0[posIndexEta]->Fill(posPorPt, posNcl);
          nClVsPPiV0[negIndexEta]->Fill(negPorPt, negNcl);
          nClVsdEdxPiV0[negIndexEta]->Fill(negNcl, negTrkdEdx);
          nClVsdEdxpPiV0[negIndexEta]->Fill(negNcl, negTrkdEdx);
          nClVsdEdxPiV0[posIndexEta]->Fill(posNcl, posTrkdEdx);
          nClVsdEdxpPiV0[posIndexEta]->Fill(posNcl, posTrkdEdx);
          nClVsPpPiV0[negIndexEta]->Fill(negPorPt, negNcl);
          dEdxPiV0[posIndexEta]->Fill(posTrkP, posTrkdEdx, centrality);
          dEdxPiV0[negIndexEta]->Fill(negTrkP, negTrkdEdx, centrality);

          if (posTrkP > kMinPMIP && posTrkP < kMaxPMIP && posTrkdEdx > kMindEdxMIP && posTrkdEdx < kMaxdEdxMIP) {
            registry.fill(HIST("dEdxVsEtaPiMIPV0"), posTrkEta, posTrkdEdx);
            registry.fill(HIST("dEdxVsEtaPiMIPV0p"), posTrkEta, posTrkdEdx);
          }
          if (negTrkP > kMinPMIP && negTrkP < kMaxPMIP && negTrkdEdx > kMindEdxMIP && negTrkdEdx < kMaxdEdxMIP) {
            registry.fill(HIST("dEdxVsEtaPiMIPV0"), negTrkEta, negTrkdEdx);
            registry.fill(HIST("dEdxVsEtaPiMIPV0p"), negTrkEta, negTrkdEdx);
          }
        }
      }

      if (passesLambdaSelection(collision, v0)) {
        registry.fill(HIST("V0sCounter"), V0sCounter::Lambda);
        registry.fill(HIST("ArmL"), alpha, qT);
        registry.fill(HIST("MassLVsPt"), v0.pt(), v0.mLambda());
        if (posTrackCharge > kZero) {
          registry.fill(HIST("nSigPrFromL"), posTrkPt, posTrack.tpcNSigmaPr());
          registry.fill(HIST("NclVsEtaPrV0"), posTrkEta, posNcl);
          registry.fill(HIST("NclVsEtaPrV0p"), posTrkEta, posNcl);
          nClVsPPrV0[posIndexEta]->Fill(posPorPt, posNcl);
          nClVsPpPrV0[posIndexEta]->Fill(posPorPt, posNcl);
          nClVsdEdxPrV0[posIndexEta]->Fill(posNcl, posTrkdEdx);
          nClVsdEdxpPrV0[posIndexEta]->Fill(posNcl, posTrkdEdx);
          dEdxPrV0[posIndexEta]->Fill(posTrkP, posTrkdEdx, centrality);
        }
        if (negTrackCharge < kZero) {
          registry.fill(HIST("nSigPiFromL"), negTrkPt, negTrack.tpcNSigmaPi());
          // registry.fill(HIST("NclVsEtaPiV0"), negTrkEta, negNcl);
          // registry.fill(HIST("NclVsEtaPiV0p"), negTrkEta, negNcl);
          // nClVsPPiV0[negIndexEta]->Fill(negPorPt, negNcl);
          // nClVsPpPiV0[negIndexEta]->Fill(negPorPt, negNcl);
          // nClVsdEdxPiV0[negIndexEta]->Fill(negNcl, negTrkdEdx);
          // nClVsdEdxpPiV0[negIndexEta]->Fill(negNcl, negTrkdEdx);
          // dEdxPiV0[negIndexEta]->Fill(negTrkP, negTrkdEdx, centrality);
        }
      }

      if (passesAntiLambdaSelection(collision, v0)) {
        registry.fill(HIST("V0sCounter"), V0sCounter::AntiLambda);
        registry.fill(HIST("ArmAL"), alpha, qT);
        registry.fill(HIST("MassALVsPt"), v0.pt(), v0.mAntiLambda());
        if (posTrackCharge > kZero) {
          registry.fill(HIST("nSigPiFromAL"), posTrkPt, posTrack.tpcNSigmaPi());
          // registry.fill(HIST("NclVsEtaPiV0"), posTrkEta, posNcl);
          // registry.fill(HIST("NclVsEtaPiV0p"), posTrkEta, posNcl);
          // nClVsPPiV0[posIndexEta]->Fill(posPorPt, posNcl);
          // nClVsPpPiV0[posIndexEta]->Fill(posPorPt, posNcl);
          // nClVsdEdxPiV0[posIndexEta]->Fill(posNcl, posTrkdEdx);
          // nClVsdEdxpPiV0[posIndexEta]->Fill(posNcl, posTrkdEdx);
          // dEdxPiV0[posIndexEta]->Fill(posTrkP, posTrkdEdx, centrality);
        }
        if (negTrackCharge < kZero) {
          registry.fill(HIST("nSigPrFromAL"), negTrkPt, negTrack.tpcNSigmaPr());
          registry.fill(HIST("NclVsEtaPrV0"), negTrkEta, negNcl);
          registry.fill(HIST("NclVsEtaPrV0p"), negTrkEta, negNcl);
          nClVsPPrV0[negIndexEta]->Fill(negPorPt, negNcl);
          nClVsPpPrV0[negIndexEta]->Fill(negPorPt, negNcl);
          nClVsdEdxPrV0[negIndexEta]->Fill(negNcl, negTrkdEdx);
          nClVsdEdxpPrV0[negIndexEta]->Fill(negNcl, negTrkdEdx);
          dEdxPrV0[negIndexEta]->Fill(negTrkP, negTrkdEdx, centrality);
        }
      }

      if (passesGammaSelection(collision, v0)) {
        if (std::abs(alpha) < v0Selections.armAlphaSel && qT < v0Selections.qTSel) {

          if (v0Selections.applyPlateauSel && etaCal.isCalPlateauLoaded) {
            const double posDedxCal{etaCal.pEtaCalPlateau->GetBinContent(etaCal.pEtaCalPlateau->FindBin(posTrkEta))};
            const double negDedxCal{etaCal.pEtaCalPlateau->GetBinContent(etaCal.pEtaCalPlateau->FindBin(negTrkEta))};
            if (!(std::abs(posTrkdEdx - posDedxCal) < v0Selections.dEdxPlateauSel && std::abs(negTrkdEdx - negDedxCal) < v0Selections.dEdxPlateauSel))
              continue;
          }

          registry.fill(HIST("V0sCounter"), V0sCounter::Gamma);
          registry.fill(HIST("ArmG"), alpha, qT);
          registry.fill(HIST("MassGVsPt"), v0.pt(), v0.mGamma());
          registry.fill(HIST("nSigElFromG"), negTrkPt, negTrack.tpcNSigmaEl());
          registry.fill(HIST("nSigElFromG"), posTrkPt, posTrack.tpcNSigmaEl());
          registry.fill(HIST("NclVsEtaElV0"), posTrkEta, posNcl);
          registry.fill(HIST("NclVsEtaElV0p"), posTrkEta, posNcl);
          registry.fill(HIST("NclVsEtaElV0"), negTrkEta, negNcl);
          registry.fill(HIST("NclVsEtaElV0p"), negTrkEta, negNcl);
          nClVsPElV0[negIndexEta]->Fill(negPorPt, negNcl);
          nClVsPpElV0[negIndexEta]->Fill(negPorPt, negNcl);
          nClVsPElV0[posIndexEta]->Fill(posPorPt, posNcl);
          nClVsPpElV0[posIndexEta]->Fill(posPorPt, posNcl);
          nClVsdEdxElV0[negIndexEta]->Fill(negNcl, negTrkdEdx);
          nClVsdEdxpElV0[negIndexEta]->Fill(negNcl, negTrkdEdx);
          nClVsdEdxElV0[posIndexEta]->Fill(posNcl, posTrkdEdx);
          nClVsdEdxpElV0[posIndexEta]->Fill(posNcl, posTrkdEdx);
          registry.fill(HIST("dEdxVsEtaElMIPV0"), posTrkEta, posTrkdEdx);
          registry.fill(HIST("dEdxVsEtaElMIPV0p"), posTrkEta, posTrkdEdx);
          registry.fill(HIST("dEdxVsEtaElMIPV0"), negTrkEta, negTrkdEdx);
          registry.fill(HIST("dEdxVsEtaElMIPV0p"), negTrkEta, negTrkdEdx);
          dEdxElV0[posIndexEta]->Fill(posTrkP, posTrkdEdx, centrality);
          dEdxElV0[negIndexEta]->Fill(negTrkP, negTrkdEdx, centrality);
        }
      }
    } // v0s
  }
  PROCESS_SWITCH(PiKpRAA, processCalibrationAndV0s, "Process QA", true);

  Preslice<TracksMC> perCollision = aod::track::collisionId;
  Service<o2::framework::O2DatabasePDG> pdg;
  void processMC(aod::McCollisions::iterator const& mccollision, soa::SmallGroups<ColEvSelsMC> const& collisions, BCsRun3 const& /*bcs*/, aod::FT0s const& /*ft0s*/, aod::McParticles const& mcParticles, TracksMC const& tracksMC)
  {
    for (const auto& collision : collisions) {
      // Event selection
      if (!isEventSelected(collision)) {
        continue;
      }
      // MC collision?
      if (!collision.has_mcCollision()) {
        continue;
      }

      registry.fill(HIST("EventCounterMC"), EvCutLabel::All);

      if (std::fabs(mccollision.posZ()) > posZcut)
        continue;

      registry.fill(HIST("zPos"), collision.posZ());
      registry.fill(HIST("zPosMC"), mccollision.posZ());
      registry.fill(HIST("EventCounterMC"), EvCutLabel::VtxZ);

      const auto& foundBC = collision.foundBC_as<BCsRun3>();
      uint64_t timeStamp{foundBC.timestamp()};
      const int magField{getMagneticField(timeStamp)};

      if (v0Selections.applyPhiCut) {
        const int nextRunNumber{foundBC.runNumber()};
        if (currentRunNumberPhiSel != nextRunNumber) {
          loadPhiCutSelections(timeStamp);
          currentRunNumberPhiSel = nextRunNumber;
          LOG(info) << "\tcurrentRunNumberPhiSel= " << currentRunNumberPhiSel << " timeStamp = " << timeStamp;
        }

        // return if phi cut objects are nullptr
        if (!(phiCut.hPhiCutHigh && phiCut.hPhiCutLow))
          return;
      }

      // Remove collisions if reconstructed more than once
      if (skipRecoColGTOne && (collisions.size() > kOne))
        continue;

      const float centrality{isT0Ccent ? collision.centFT0C() : collision.centFT0M()};
      registry.fill(HIST("T0Ccent"), centrality);

      const auto& groupedTracks{tracksMC.sliceBy(perCollision, collision.globalIndex())};

      // Track selection with Open DCAxy
      for (const auto& track : groupedTracks) {
        // Track Selection
        if (track.eta() < v0Selections.minEtaDaughter || track.eta() > v0Selections.maxEtaDaughter)
          continue;

        if (track.pt() < v0Selections.minPt || track.pt() > v0Selections.maxPt)
          continue;

        if (!trkSelGlobalOpenDCAxy.IsSelected(track))
          continue;

        if (!track.has_mcParticle())
          continue;

        // Get the MC particle
        const auto& particle{track.mcParticle()};
        auto charge{0.};
        auto* pdgParticle = pdg->GetParticle(particle.pdgCode());
        if (pdgParticle != nullptr)
          charge = pdgParticle->Charge();
        else
          continue;

        // Is it a charged particle?
        if (std::abs(charge) < kMinCharge)
          continue;

        float phiPrime{track.phi()};
        phiPrimeFunc(phiPrime, magField, charge);

        const float pOrPt{v0Selections.usePinPhiSelection ? track.p() : track.pt()};
        if (v0Selections.applyPhiCut) {
          if (!passesPhiSelection(pOrPt, phiPrime))
            continue;
        }

        const int16_t nclFound{track.tpcNClsFound()};
        const int16_t nclPID{track.tpcNClsPID()};
        const int16_t ncl = v0Selections.useNclsPID ? nclPID : nclFound;
        if (v0Selections.applyNclSel && ncl < v0Selections.minNcl)
          continue;

        bool isPrimary{false};
        bool isDecay{false};
        bool isMaterial{false};
        if (particle.isPhysicalPrimary())
          isPrimary = true;
        else if (particle.getProcess() == TMCProcess::kPDecay)
          isDecay = true;
        else
          isMaterial = true;

        bool isPi{false};
        bool isPr{false};
        if (particle.pdgCode() == PDG_t::kPiPlus || particle.pdgCode() == PDG_t::kPiMinus)
          isPi = true;
        else if (particle.pdgCode() == PDG_t::kProton || particle.pdgCode() == PDG_t::kProtonBar)
          isPr = true;
        else
          continue;

        if (isPrimary && !isDecay && !isMaterial) {
          if (isPi && !isPr)
            registry.fill(HIST("dcaVsPtPi"), track.pt(), track.dcaXY(), centrality);
          if (isPr && !isPi)
            registry.fill(HIST("dcaVsPtPr"), track.pt(), track.dcaXY(), centrality);
        }

        if (isDecay && !isPrimary && !isMaterial) {
          if (isPi && !isPr)
            registry.fill(HIST("dcaVsPtPiDec"), track.pt(), track.dcaXY(), centrality);
          if (isPr && !isPi)
            registry.fill(HIST("dcaVsPtPrDec"), track.pt(), track.dcaXY(), centrality);
        }

        if (isMaterial && !isPrimary && !isDecay) {
          if (isPi && !isPr)
            registry.fill(HIST("dcaVsPtPiMat"), track.pt(), track.dcaXY(), centrality);
          if (isPr && !isPi)
            registry.fill(HIST("dcaVsPtPrMat"), track.pt(), track.dcaXY(), centrality);
        }
      }

      // Global track + DCAxy selections
      for (const auto& track : groupedTracks) {
        // Track Selection
        if (track.eta() < v0Selections.minEtaDaughter || track.eta() > v0Selections.maxEtaDaughter)
          continue;

        if (track.pt() < v0Selections.minPt || track.pt() > v0Selections.maxPt)
          continue;

        if (!trkSelGlobal.IsSelected(track))
          continue;

        // Has MC particle?
        if (!track.has_mcParticle())
          continue;

        // Get the MC particle
        const auto& particle{track.mcParticle()};
        auto charge{0.};
        auto* pdgParticle = pdg->GetParticle(particle.pdgCode());
        if (pdgParticle != nullptr)
          charge = pdgParticle->Charge();
        else
          continue;

        // Is it a charged particle?
        if (std::abs(charge) < kMinCharge)
          continue;

        float phiPrime{track.phi()};
        phiPrimeFunc(phiPrime, magField, charge);

        const float pOrPt{v0Selections.usePinPhiSelection ? track.p() : track.pt()};
        if (v0Selections.applyPhiCut) {
          if (!passesPhiSelection(pOrPt, phiPrime))
            continue;
        }

        const int16_t nclFound{track.tpcNClsFound()};
        const int16_t nclPID{track.tpcNClsPID()};
        const int16_t ncl = v0Selections.useNclsPID ? nclPID : nclFound;
        if (v0Selections.applyNclSel && ncl < v0Selections.minNcl)
          continue;

        int indexEta{-999};
        const float eta{track.eta()};
        for (int i = 0; i < kNEtaHists; ++i) {
          if (eta >= kLowEta[i] && eta < kHighEta[i]) {
            indexEta = i;
            break;
          }
        }

        if (indexEta < kZeroInt || indexEta > kSevenInt)
          continue;

        registry.fill(HIST("NclVsPhip"), pOrPt, phiPrime, ncl);
        registry.fill(HIST("NclVsEtaPID"), eta, ncl);
        registry.fill(HIST("NclVsEtaPIDp"), eta, ncl);

        bool isPrimary{false};
        if (particle.isPhysicalPrimary())
          isPrimary = true;

        if (!isPrimary)
          continue;

        bool isPi{false};
        bool isKa{false};
        bool isPr{false};

        if (particle.pdgCode() == PDG_t::kPiPlus || particle.pdgCode() == PDG_t::kPiMinus)
          isPi = true;
        else if (particle.pdgCode() == PDG_t::kKPlus || particle.pdgCode() == PDG_t::kKMinus)
          isKa = true;
        else if (particle.pdgCode() == PDG_t::kProton || particle.pdgCode() == PDG_t::kProtonBar)
          isPr = true;
        else
          continue;

        if (isPi && !isKa && !isPr)
          registry.fill(HIST("PtPiVsCent"), track.pt(), centrality);
        if (isKa && !isPi && !isPr)
          registry.fill(HIST("PtKaVsCent"), track.pt(), centrality);
        if (isPr && !isPi && !isKa)
          registry.fill(HIST("PtPrVsCent"), track.pt(), centrality);
      }

      // Generated MC
      for (const auto& particle : mcParticles) {
        if (particle.eta() < v0Selections.minEtaDaughter || particle.eta() > v0Selections.maxEtaDaughter)
          continue;

        if (particle.pt() < v0Selections.minPt || particle.pt() > v0Selections.maxPt)
          continue;

        auto charge{0.};
        // Get the MC particle
        auto* pdgParticle = pdg->GetParticle(particle.pdgCode());
        if (pdgParticle != nullptr)
          charge = pdgParticle->Charge();
        else
          continue;

        // Is it a charged particle?
        if (std::abs(charge) < kMinCharge)
          continue;

        // Is it a primary particle?
        bool isPrimary{true};
        if (!particle.isPhysicalPrimary())
          isPrimary = false;

        if (isPrimary) {
          if (particle.pdgCode() == PDG_t::kPiPlus || particle.pdgCode() == PDG_t::kPiMinus) // pion
            registry.fill(HIST("PtPiVsCentMC"), particle.pt(), centrality);
          else if (particle.pdgCode() == PDG_t::kKPlus || particle.pdgCode() == PDG_t::kKMinus) // kaon
            registry.fill(HIST("PtKaVsCentMC"), particle.pt(), centrality);
          else if (particle.pdgCode() == PDG_t::kProton || particle.pdgCode() == PDG_t::kProtonBar) // proton
            registry.fill(HIST("PtPrVsCentMC"), particle.pt(), centrality);
          else
            continue;
        }
      }
    } // Collisions
  }
  PROCESS_SWITCH(PiKpRAA, processMC, "Process MC closure", false);

  void processSim(aod::McCollisions::iterator const& mccollision, soa::SmallGroups<ColEvSelsMC> const& collisions, BCsRun3 const& /*bcs*/, aod::FT0s const& /*ft0s*/, aod::McParticles const& mcParticles, TracksMC const& tracksMC)
  {

    //---------------------------
    // Only INEL > 0 generated collisions
    // By counting number of primary charged particles in |eta| < 1
    //---------------------------
    int nChMC{0};
    int nChFT0A{0};
    int nChFT0C{0};
    for (const auto& particle : mcParticles) {

      auto charge{0.};
      // Get the MC particle
      const auto* pdgParticle = pdg->GetParticle(particle.pdgCode());
      if (pdgParticle != nullptr) {
        charge = pdgParticle->Charge();
      } else {
        continue;
      }

      // Is it a charged particle?
      if (std::abs(charge) < kMinCharge)
        continue;

      // Is it a primary particle?
      if (!particle.isPhysicalPrimary())
        continue;

      const float eta{particle.eta()};

      // TVX requirement
      if (eta > kMinFT0A && eta < kMaxFT0A) {
        nChFT0A++;
      }

      if (eta > kMinFT0C && eta < kMaxFT0C) {
        nChFT0C++;
      }

      // INEL > 0
      if (std::abs(eta) > kOne)
        continue;

      nChMC++;
    }

    //---------------------------
    // Only events with at least one charged particle in the FT0A and FT0C acceptances
    //---------------------------
    if (selTVXMC) {
      if (!(nChFT0A > kZeroInt && nChFT0C > kZeroInt)) {
        registry.fill(HIST("NchMCcentVsTVX"), nChMC, 0.5);
        return;
      }
      registry.fill(HIST("NchMCcentVsTVX"), nChMC, 1.5);
    }

    //---------------------------
    // Only MC events with |Vtx Z| < 10 cm
    //---------------------------
    if (isZvtxPosSelMC && (std::fabs(mccollision.posZ()) > posZcut)) {
      return;
    }

    //---------------------------
    // Only INEL > 0 generated events
    //---------------------------
    if (selINELgt0) {
      if (!(nChMC > kZeroInt)) {
        return;
      }
    }

    const auto& nRecColls{collisions.size()};
    registry.fill(HIST("NumberOfRecoCollisions"), nRecColls);

    //---------------------------
    // Only Generated evets with at least one reconstrued collision
    //---------------------------
    if (nRecColls > kZeroInt) {

      // Finds the collisions with the largest number of contributors
      // in case nRecColls is larger than One
      int biggestNContribs{-1};
      int bestCollisionIndex{-1};
      for (const auto& collision : collisions) {

        const float centrality{isT0Ccent ? collision.centFT0C() : collision.centFT0M()};
        if (selHasFT0 && !collision.has_foundFT0()) {
          continue;
        }

        if (biggestNContribs < collision.numContrib()) {
          biggestNContribs = collision.numContrib();
          bestCollisionIndex = collision.globalIndex();
        }

        // Needed to calculate denominator of the Event Splitting correction
        if (isEventSelected(collision)) {
          registry.fill(HIST("Centrality_AllRecoEvt"), centrality);
        }
      }

      //---------------------------
      // Loop over the reconstructed collisions
      // Only that one with the largest number of contributors is considered
      //---------------------------
      for (const auto& collision : collisions) {

        const float centrality{isT0Ccent ? collision.centFT0C() : collision.centFT0M()};

        //---------------------------
        // Reject collisions if has_foundFT0() returns false
        //---------------------------
        if (selHasFT0 && !collision.has_foundFT0()) {
          registry.fill(HIST("T0CcentVsFoundFT0"), centrality, 0.5);
          continue;
        }
        registry.fill(HIST("T0CcentVsFoundFT0"), centrality, 1.5);

        //---------------------------
        // Pick the collisions with the largest number of contributors
        //---------------------------
        if (bestCollisionIndex != collision.globalIndex()) {
          continue;
        }

        // Needed to load the Phi selection from the CCDB
        const auto& foundBC = collision.foundBC_as<BCsRun3>();
        uint64_t timeStamp{foundBC.timestamp()};
        const int magField{getMagneticField(timeStamp)};

        if (v0Selections.applyPhiCut) {
          const int nextRunNumber{foundBC.runNumber()};
          if (currentRunNumberPhiSel != nextRunNumber) {
            loadPhiCutSelections(timeStamp);
            currentRunNumberPhiSel = nextRunNumber;
            LOG(info) << "\tcurrentRunNumberPhiSel= " << currentRunNumberPhiSel << " timeStamp = " << timeStamp;
          }

          // return if phi cut objects are nullptr
          if (!(phiCut.hPhiCutHigh && phiCut.hPhiCutLow))
            return;
        }

        //---------------------------
        // Needed to construct the correlation between MC Nch v.s. centrality
        //---------------------------

        registry.fill(HIST("Centrality_WRecoEvt"), centrality);
        registry.fill(HIST("zPosMC"), mccollision.posZ());

        //---------------------------
        // All Generated events with at least one associated reconstructed collision
        // The Generated events are not subjected to any selection criteria
        // This histograms are used for the denominator of the tracking efficiency
        //---------------------------
        for (const auto& particle : mcParticles) {
          if (particle.eta() < v0Selections.minEtaDaughter || particle.eta() > v0Selections.maxEtaDaughter)
            continue;

          if (particle.pt() < v0Selections.minPt || particle.pt() > v0Selections.maxPt)
            continue;

          auto charge{0.};
          // Get the MC particle
          auto* pdgParticle = pdg->GetParticle(particle.pdgCode());
          if (pdgParticle != nullptr) {
            charge = pdgParticle->Charge();
          } else {
            continue;
          }

          // Is it a charged particle?
          if (std::abs(charge) < kMinCharge)
            continue;

          // Is it a primary particle?
          bool isPrimary{true};
          if (!particle.isPhysicalPrimary())
            isPrimary = false;

          if (isPrimary) {
            if (particle.pdgCode() == PDG_t::kPiPlus || particle.pdgCode() == PDG_t::kPiMinus) {
              registry.fill(HIST("PtPiVsCentMC_WithRecoEvt"), particle.pt(), centrality); // Denominator of tracking efficiency
              registry.fill(HIST("PtPiVsNchMC_WithRecoEvt"), particle.pt(), nChMC);       // Numerator of signal loss
            } else if (particle.pdgCode() == PDG_t::kKPlus || particle.pdgCode() == PDG_t::kKMinus) {
              registry.fill(HIST("PtKaVsCentMC_WithRecoEvt"), particle.pt(), centrality);
              registry.fill(HIST("PtKaVsNchMC_WithRecoEvt"), particle.pt(), nChMC);
            } else if (particle.pdgCode() == PDG_t::kProton || particle.pdgCode() == PDG_t::kProtonBar) {
              registry.fill(HIST("PtPrVsCentMC_WithRecoEvt"), particle.pt(), centrality);
              registry.fill(HIST("PtPrVsNchMC_WithRecoEvt"), particle.pt(), nChMC);
            } else {
              continue;
            }
          }
        } // Loop over generated particles per generated collision

        //---------------------------
        // Reconstructed collisions subjected to selection criteria
        //---------------------------

        // Event selection
        if (!isEventSelected(collision)) {
          continue;
        }

        registry.fill(HIST("Centrality_WRecoEvtWSelCri"), centrality);
        registry.fill(HIST("NchMCVsCent"), centrality, nChMC);
        registry.fill(HIST("NchMC_WithRecoEvt"), nChMC);
        registry.fill(HIST("T0Ccent"), centrality);
        registry.fill(HIST("zPos"), collision.posZ());

        const auto& groupedTracks{tracksMC.sliceBy(perCollision, collision.globalIndex())};

        //---------------------------
        // Track selection with Open DCAxy
        //  This is needed for the Secondary Particle Correction
        //---------------------------
        for (const auto& track : groupedTracks) {
          // Track Selection
          if (track.eta() < v0Selections.minEtaDaughter || track.eta() > v0Selections.maxEtaDaughter)
            continue;

          if (track.pt() < v0Selections.minPt || track.pt() > v0Selections.maxPt)
            continue;

          if (!trkSelGlobalOpenDCAxy.IsSelected(track))
            continue;

          if (!track.has_mcParticle())
            continue;

          // Get the MC particle
          const auto& particle{track.mcParticle()};
          auto charge{0.};
          auto* pdgParticle = pdg->GetParticle(particle.pdgCode());
          if (pdgParticle != nullptr) {
            charge = pdgParticle->Charge();
          } else {
            continue;
          }

          // Is it a charged particle?
          if (std::abs(charge) < kMinCharge)
            continue;

          float phiPrime{track.phi()};
          phiPrimeFunc(phiPrime, magField, charge);

          const float pOrPt{v0Selections.usePinPhiSelection ? track.p() : track.pt()};
          if (v0Selections.applyPhiCut) {
            if (!passesPhiSelection(pOrPt, phiPrime))
              continue;
          }

          const int16_t nclFound{track.tpcNClsFound()};
          const int16_t nclPID{track.tpcNClsPID()};
          const int16_t ncl = v0Selections.useNclsPID ? nclPID : nclFound;
          if (v0Selections.applyNclSel && ncl < v0Selections.minNcl)
            continue;

          bool isPrimary{false};
          bool isDecay{false};
          bool isMaterial{false};
          if (particle.isPhysicalPrimary()) {
            isPrimary = true;
          } else if (particle.getProcess() == TMCProcess::kPDecay) {
            isDecay = true;
          } else {
            isMaterial = true;
          }

          bool isPi{false};
          bool isPr{false};
          if (particle.pdgCode() == PDG_t::kPiPlus || particle.pdgCode() == PDG_t::kPiMinus) {
            isPi = true;
          } else if (particle.pdgCode() == PDG_t::kProton || particle.pdgCode() == PDG_t::kProtonBar) {
            isPr = true;
          } else {
            continue;
          }

          if (isPrimary && !isDecay && !isMaterial) {
            if (isPi && !isPr)
              registry.fill(HIST("dcaVsPtPi"), track.pt(), track.dcaXY(), centrality);
            if (isPr && !isPi)
              registry.fill(HIST("dcaVsPtPr"), track.pt(), track.dcaXY(), centrality);
          }

          if (isDecay && !isPrimary && !isMaterial) {
            if (isPi && !isPr)
              registry.fill(HIST("dcaVsPtPiDec"), track.pt(), track.dcaXY(), centrality);
            if (isPr && !isPi)
              registry.fill(HIST("dcaVsPtPrDec"), track.pt(), track.dcaXY(), centrality);
          }

          if (isMaterial && !isPrimary && !isDecay) {
            if (isPi && !isPr)
              registry.fill(HIST("dcaVsPtPiMat"), track.pt(), track.dcaXY(), centrality);
            if (isPr && !isPi)
              registry.fill(HIST("dcaVsPtPrMat"), track.pt(), track.dcaXY(), centrality);
          }
        }

        //---------------------------
        // Global track + DCAxy selections
        // This is needed for the number of the Tracking Efficiency
        // and the spectra to be corrected
        //---------------------------
        int nCh{0};
        for (const auto& track : groupedTracks) {
          // Track Selection
          if (track.eta() < v0Selections.minEtaDaughter || track.eta() > v0Selections.maxEtaDaughter)
            continue;

          if (track.pt() < v0Selections.minPt || track.pt() > v0Selections.maxPt)
            continue;

          if (!trkSelGlobal.IsSelected(track))
            continue;

          // Has MC particle?
          if (!track.has_mcParticle())
            continue;

          // Get the MC particle
          const auto& particle{track.mcParticle()};
          auto charge{0.};
          auto* pdgParticle = pdg->GetParticle(particle.pdgCode());
          if (pdgParticle != nullptr) {
            charge = pdgParticle->Charge();
          } else {
            continue;
          }

          // Is it a charged particle?
          if (std::abs(charge) < kMinCharge)
            continue;

          float phiPrime{track.phi()};
          phiPrimeFunc(phiPrime, magField, charge);

          const float pOrPt{v0Selections.usePinPhiSelection ? track.p() : track.pt()};
          if (v0Selections.applyPhiCut) {
            if (!passesPhiSelection(pOrPt, phiPrime))
              continue;
          }

          const int16_t nclFound{track.tpcNClsFound()};
          const int16_t nclPID{track.tpcNClsPID()};
          const int16_t ncl = v0Selections.useNclsPID ? nclPID : nclFound;
          if (v0Selections.applyNclSel && ncl < v0Selections.minNcl)
            continue;

          int indexEta{-999};
          const float eta{track.eta()};
          for (int i = 0; i < kNEtaHists; ++i) {
            if (eta >= kLowEta[i] && eta < kHighEta[i]) {
              indexEta = i;
              break;
            }
          }

          if (indexEta < kZeroInt || indexEta > kSevenInt)
            continue;

          registry.fill(HIST("NclVsPhip"), pOrPt, phiPrime, ncl);
          registry.fill(HIST("NclVsEtaPID"), eta, ncl);
          registry.fill(HIST("NclVsEtaPIDp"), eta, ncl);
          nCh++;

          bool isPrimary{false};
          if (particle.isPhysicalPrimary())
            isPrimary = true;

          if (!isPrimary)
            continue;

          bool isPi{false};
          bool isKa{false};
          bool isPr{false};

          if (particle.pdgCode() == PDG_t::kPiPlus || particle.pdgCode() == PDG_t::kPiMinus) {
            isPi = true;
          } else if (particle.pdgCode() == PDG_t::kKPlus || particle.pdgCode() == PDG_t::kKMinus) {
            isKa = true;
          } else if (particle.pdgCode() == PDG_t::kProton || particle.pdgCode() == PDG_t::kProtonBar) {
            isPr = true;
          } else {
            continue;
          }

          if (isPi && !isKa && !isPr) {
            registry.fill(HIST("PtPiVsCent_WithRecoEvt"), track.pt(), centrality); // Numerator of tracking efficiency
            registry.fill(HIST("MCclosure_PtPiVsNchMC"), track.pt(), nChMC);
          }
          if (isKa && !isPi && !isPr) {
            registry.fill(HIST("PtKaVsCent_WithRecoEvt"), track.pt(), centrality);
            registry.fill(HIST("MCclosure_PtKaVsNchMC"), track.pt(), nChMC);
          }
          if (isPr && !isPi && !isKa) {
            registry.fill(HIST("PtPrVsCent_WithRecoEvt"), track.pt(), centrality);
            registry.fill(HIST("MCclosure_PtPrVsNchMC"), track.pt(), nChMC);
          }
          registry.fill(HIST("PtResolution"), particle.pt(), (track.pt() - particle.pt()) / particle.pt());
        } // Loop over reconstructed tracks
        registry.fill(HIST("NchVsCent"), centrality, nCh);
      } // Loop over Reco. Collisions: Only the collisions with the largest number of contributors
    } // If condition: Only simulated evets with at least one reconstrued collision

    //---------------------------
    // All Generated events irrespective of whether there is an associated reconstructed collision
    // Consequently, the centrality being a reconstructed quantity, might not always be available
    // Therefore it is expressed as a function of the generated pT and the generated Nch in ∣eta∣ < 0.8
    //---------------------------

    //---------------------------
    // Generated Pt spectra of all INEL > 0 Generated evets
    // irrespective of whether there is a reconstructed collision
    // This is used for the denominator of the signal loss correction
    // Also for MC closure: True Pt vs Generated Nch
    //---------------------------
    for (const auto& particle : mcParticles) {
      if (particle.eta() < v0Selections.minEtaDaughter || particle.eta() > v0Selections.maxEtaDaughter)
        continue;

      if (particle.pt() < v0Selections.minPt || particle.pt() > v0Selections.maxPt)
        continue;

      auto charge{0.};
      // Get the MC particle
      auto* pdgParticle = pdg->GetParticle(particle.pdgCode());
      if (pdgParticle != nullptr) {
        charge = pdgParticle->Charge();
      } else {
        continue;
      }

      // Is it a charged particle?
      if (std::abs(charge) < kMinCharge)
        continue;

      // Is it a primary particle?
      bool isPrimary{true};
      if (!particle.isPhysicalPrimary())
        isPrimary = false;

      if (isPrimary) {
        if (particle.pdgCode() == PDG_t::kPiPlus || particle.pdgCode() == PDG_t::kPiMinus) {
          registry.fill(HIST("PtPiVsNchMC_AllGen"), particle.pt(), nChMC);
          registry.fill(HIST("MCclosure_PtMCPiVsNchMC"), particle.pt(), nChMC);
        } else if (particle.pdgCode() == PDG_t::kKPlus || particle.pdgCode() == PDG_t::kKMinus) {
          registry.fill(HIST("PtKaVsNchMC_AllGen"), particle.pt(), nChMC);
          registry.fill(HIST("MCclosure_PtMCKaVsNchMC"), particle.pt(), nChMC);
        } else if (particle.pdgCode() == PDG_t::kProton || particle.pdgCode() == PDG_t::kProtonBar) {
          registry.fill(HIST("PtPrVsNchMC_AllGen"), particle.pt(), nChMC);
          registry.fill(HIST("MCclosure_PtMCPrVsNchMC"), particle.pt(), nChMC);
        } else {
          continue;
        }
      }
    } // Loop over Generated Particles

    //---------------------------
    //  This is used for the denominator of the event loss correction
    //---------------------------
    registry.fill(HIST("NchMC_AllGen"), nChMC);
  }
  PROCESS_SWITCH(PiKpRAA, processSim, "Process Sim", false);

  template <typename T, typename U>
  void getArmeterosVariables(const T& ppos, const T& pneg, U& alpha, U& qT)
  {

    alpha = 0., qT = 0.;
    TVector3 pV0 = ppos + pneg;
    double pV0mag = pV0.Mag();
    if (pV0mag < kTenToMinusNine)
      return; // protect against zero momentum

    const TVector3 u = pV0 * (1.0 / pV0mag);

    double pLpos = ppos.Dot(u);
    double pLneg = pneg.Dot(u);

    // qT: transverse momentum of the + track w.r.t. V0 direction
    TVector3 pTpos = ppos - pLpos * u;
    qT = pTpos.Mag();

    // α: longitudinal asymmetry (uses + and − labels by charge)
    double denom = pLpos + pLneg;
    if (std::abs(denom) < kTenToMinusNine)
      return; // avoid 0 division (unphysical for V0s)

    alpha = (pLpos - pLneg) / denom; // equivalently / pV0mag
  }

  template <typename T>
  bool passesTrackSelectionDaughters(const T& track)
  {

    // Decay daughters used Global tracks excluding the DCAxy selection; this approach was used in the Run 2 analyses
    // https://github.com/AliceO2Group/O2Physics/blob/b178c96d03ede873ee769ef8a4d7c1e21bf78332/Common/Core/TrackSelectionDefaults.cxx

    // In the V0s analysis only the selection on eta (|eta|<0.8) and the selection on number of crossed rows was used: nCrossedRows >= 70
    const float pt{track.pt()};
    const float eta{track.eta()};
    const int16_t nCrossedRows{track.tpcNClsCrossedRows()};
    const float nCrossedRowsOverFindableCls{track.tpcCrossedRowsOverFindableCls()};
    const float chi2PerClsTPC{track.tpcChi2NCl()};
    const float chi2PerClsITS{track.itsChi2NCl()};
    const bool refitITS{track.passedITSRefit()};
    const bool refitTPC{track.passedTPCRefit()};
    const bool goldeChi2{track.passedGoldenChi2()};
    const bool oneHitITSib{track.passedITSHitsFB1()};

    bool etaSel{false};
    bool ptSel{false};
    bool xRows{false};
    bool xRowsToFindCls{false};
    bool chi2TPC{false};
    bool chi2ITS{false};
    bool itsrefit{false};
    bool tpcrefit{false};
    bool golden{false};
    bool itshit{false};

    registry.fill(HIST("TrackDaughterCounter"), TrkSelLabel::AllTrks);
    if (v0Selections.useOfficialV0sSelOfDaughters) {
      if (eta > v0Selections.minEtaDaughter && eta < v0Selections.maxEtaDaughter) {
        registry.fill(HIST("TrackDaughterCounter"), TrkSelLabel::Eta);
        etaSel = true;
      }
      if (nCrossedRows >= v0Selections.minNCrossedRows) {
        registry.fill(HIST("TrackDaughterCounter"), TrkSelLabel::XRows);
        xRows = true;
      }
    } else {
      if (eta > v0Selections.minEtaDaughter && eta < v0Selections.maxEtaDaughter) {
        etaSel = true;
        registry.fill(HIST("TrackDaughterCounter"), TrkSelLabel::Eta);
      }
      if (pt > v0Selections.minPtDaughter && pt < v0Selections.maxPtDaughter) {
        ptSel = true;
        registry.fill(HIST("TrackDaughterCounter"), TrkSelLabel::Pt);
      }
      if (nCrossedRows >= v0Selections.minNCrossedRows) {
        xRows = true;
        registry.fill(HIST("TrackDaughterCounter"), TrkSelLabel::XRows);
      }
      if (nCrossedRowsOverFindableCls >= v0Selections.minNCrossedRowsOverFindableCls) {
        xRowsToFindCls = true;
        registry.fill(HIST("TrackDaughterCounter"), TrkSelLabel::XRowsOverFindableCls);
      }
      if (chi2PerClsTPC < v0Selections.maxChi2ClsTPC) {
        chi2TPC = true;
        registry.fill(HIST("TrackDaughterCounter"), TrkSelLabel::Chi2TPC);
      }
      if (chi2PerClsITS < v0Selections.maxChi2ClsITS) {
        chi2ITS = true;
        registry.fill(HIST("TrackDaughterCounter"), TrkSelLabel::Chi2ITS);
      }
      if (refitITS == v0Selections.itsRefit) {
        itsrefit = true;
        registry.fill(HIST("TrackDaughterCounter"), TrkSelLabel::Itsrefit);
      }
      if (refitTPC == v0Selections.tpcRefit) {
        tpcrefit = true;
        registry.fill(HIST("TrackDaughterCounter"), TrkSelLabel::Tpcrefit);
      }
      if (goldeChi2 == v0Selections.chi2Golden) {
        golden = true;
        registry.fill(HIST("TrackDaughterCounter"), TrkSelLabel::Golden);
      }
      if (oneHitITSib == v0Selections.its1HitIB) {
        itshit = true;
        registry.fill(HIST("TrackDaughterCounter"), TrkSelLabel::Itshit);
      }
    }

    bool isSelected{false};
    if (v0Selections.useOfficialV0sSelOfDaughters) {
      isSelected = etaSel && xRows ? true : false;
    } else {
      if (!v0Selections.selElecFromGammas && v0Selections.requireITShit)
        isSelected = etaSel && ptSel && xRows && xRowsToFindCls && chi2TPC && chi2ITS && itsrefit && tpcrefit && golden && itshit ? true : false;
      if (!v0Selections.selElecFromGammas && !v0Selections.requireITShit)
        isSelected = etaSel && ptSel && xRows && xRowsToFindCls && chi2TPC && chi2ITS && itsrefit && tpcrefit && golden ? true : false;
      if (v0Selections.selElecFromGammas)
        isSelected = etaSel && ptSel && xRows && xRowsToFindCls && chi2TPC && chi2ITS && tpcrefit && golden ? true : false;
    }

    if (isSelected == true)
      registry.fill(HIST("TrackDaughterCounter"), TrkSelLabel::PassedAll);

    return isSelected;
  }

  // V0 topological selection:
  // *) V0 decay radius (cm)
  // *) V0 cosine of pointing angle
  // *) DCA between V0 daughters (cm)
  template <typename T>
  bool passesV0TopologicalSelection(const T& v0)
  {
    // bool isSelected = v0.v0radius() > v0Selections.v0radius && v0.v0radius() < v0Selections.v0radiusMax && v0.v0cosPA() > v0Selections.v0cospa && v0.dcaV0daughters() < v0Selections.dcav0dau && passesDCASelectionDaughters(v0) ? true : false;
    bool isSelected = v0.v0radius() > v0Selections.v0radius && v0.v0radius() < v0Selections.v0radiusMax && v0.v0cosPA() > v0Selections.v0cospa && v0.dcaV0daughters() < v0Selections.dcav0dau ? true : false;
    return isSelected;
  }

  // K0s topological selection:
  // *) V0 Lifetime (cm)
  // *) Rapidity selection
  // *) DCA of pions to prim. vtx (cm)
  template <typename C, typename T>
  bool passesK0Selection(const C& collision, const T& v0)
  {

    const auto& posTrack = v0.template posTrack_as<TracksFull>();
    const auto& negTrack = v0.template negTrack_as<TracksFull>();

    const double dcaPos{std::fabs(v0.dcapostopv())};
    const double dcaNeg{std::fabs(v0.dcanegtopv())};
    const float posTPCNsigma{std::fabs(posTrack.tpcNSigmaPi())};
    const float negTPCNsigma{std::fabs(negTrack.tpcNSigmaPi())};
    const double dMassK0s{std::abs(v0.mK0Short() - o2::constants::physics::MassK0Short)};
    const double dMassL{std::abs(v0.mLambda() - o2::constants::physics::MassLambda0)};
    const double dMassAL{std::abs(v0.mAntiLambda() - o2::constants::physics::MassLambda0)};
    const double dMassG{std::abs(v0.mGamma() - o2::constants::physics::MassGamma)};
    const double rapidity{std::abs(v0.yK0Short())};
    const double lifeTime{v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short};

    // Checks if DCA of daughters to PV passes selection
    const bool dcaDaugToPV{dcaPos > v0Selections.dcapostopv && dcaNeg > v0Selections.dcanegtopv ? true : false};

    // Rejects V0 if its invariant mass is not compatible with the K0s proper mass
    if (v0Selections.applyInvMassSel) {
      if (!(dMassK0s < v0Selections.dMassSel && dMassL > v0Selections.dMassSel && dMassAL > v0Selections.dMassSel && dMassG > v0Selections.dMassSelG))
        return false;
    }

    bool isSelected{false};
    if (v0Selections.useOfficialV0sSelOfDaughters) {
      isSelected = lifeTime < v0Selections.lifeTimeCutK0s && rapidity < v0Selections.rapidityCut && dcaDaugToPV ? true : false;
    } else {
      isSelected = lifeTime < v0Selections.lifeTimeCutK0s && rapidity < v0Selections.rapidityCut && posTPCNsigma < v0Selections.tpcPidNsigmaCut && negTPCNsigma < v0Selections.tpcPidNsigmaCut && dcaDaugToPV ? true : false;
    }

    if (isSelected) {
      registry.fill(HIST("EtaVsYK0s"), negTrack.eta(), v0.yK0Short());
      registry.fill(HIST("EtaVsYK0s"), posTrack.eta(), v0.yK0Short());
      registry.fill(HIST("DCAxyPtPiK0s"), v0.dcapostopv(), posTrack.pt());
      registry.fill(HIST("DCAxyPtPiK0s"), v0.dcanegtopv(), negTrack.pt());
    }
    return isSelected;
  }

  // Lambda topological selection:
  // *) V0 Lifetime (cm)
  // *) Rapidity selection
  // *) DCA of pions and protons to prim. vtx (cm)
  template <typename C, typename T>
  bool passesLambdaSelection(const C& collision, const T& v0)
  {

    const auto& posTrack = v0.template posTrack_as<TracksFull>();
    const auto& negTrack = v0.template negTrack_as<TracksFull>();

    const double dcaPos{std::fabs(v0.dcapostopv())};
    const double dcaNeg{std::fabs(v0.dcanegtopv())};
    const float posTPCNsigma{std::fabs(posTrack.tpcNSigmaPr())};
    const float negTPCNsigma{std::fabs(negTrack.tpcNSigmaPi())};
    const double dMassK0s{std::abs(v0.mK0Short() - o2::constants::physics::MassK0Short)};
    const double dMassL{std::abs(v0.mLambda() - o2::constants::physics::MassLambda0)};
    const double dMassG{std::abs(v0.mGamma() - o2::constants::physics::MassGamma)};
    const double rapidity{std::abs(v0.yLambda())};
    const double lifeTime{v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0};

    // Checks if DCA of daughters to PV passes selection
    const bool dcaDaugToPV{dcaPos > v0Selections.dcaProtonFromLambda && dcaNeg > v0Selections.dcaPionFromLambda ? true : false};

    // Rejects V0 if the invariant mass is not compatible with the Lambda proper mass
    if (v0Selections.applyInvMassSel) {
      if (!(dMassL < v0Selections.dMassSel && dMassK0s > v0Selections.dMassSel && dMassG > v0Selections.dMassSelG))
        return false;
    }

    bool isSelected{false};
    if (v0Selections.useOfficialV0sSelOfDaughters) {
      isSelected = lifeTime < v0Selections.lifeTimeCutLambda && rapidity < v0Selections.rapidityCut && dcaDaugToPV ? true : false;
    } else {
      isSelected = lifeTime < v0Selections.lifeTimeCutLambda && rapidity < v0Selections.rapidityCut && posTPCNsigma < v0Selections.tpcPidNsigmaCut && negTPCNsigma < v0Selections.tpcPidNsigmaCut && dcaDaugToPV ? true : false;
    }

    if (isSelected) {
      registry.fill(HIST("DCAxyPtPrL"), v0.dcapostopv(), posTrack.pt());
      registry.fill(HIST("EtaVsYPiL"), negTrack.eta(), v0.yLambda());
      registry.fill(HIST("EtaVsYPrL"), posTrack.eta(), v0.yLambda());
    }

    return isSelected;
  }

  // Anti Lambda topological selection:
  // *) V0 Lifetime (cm)
  // *) Rapidity selection
  // *) DCA of pions and protons to prim. vtx (cm)
  template <typename C, typename T>
  bool passesAntiLambdaSelection(const C& collision, const T& v0)
  {

    const auto& posTrack = v0.template posTrack_as<TracksFull>();
    const auto& negTrack = v0.template negTrack_as<TracksFull>();
    const double dcaPos{std::fabs(v0.dcapostopv())};
    const double dcaNeg{std::fabs(v0.dcanegtopv())};
    const float posTPCNsigma{std::fabs(posTrack.tpcNSigmaPi())};
    const float negTPCNsigma{std::fabs(negTrack.tpcNSigmaPr())};
    const double dMassK0s{std::abs(v0.mK0Short() - o2::constants::physics::MassK0Short)};
    const double dMassAL{std::abs(v0.mAntiLambda() - o2::constants::physics::MassLambda0)};
    const double dMassG{std::abs(v0.mGamma() - o2::constants::physics::MassGamma)};
    const double rapidity{std::abs(v0.yLambda())};
    const double lifeTime{v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0};

    // Checks if DCA of daughters to PV passes selection
    const bool dcaDaugToPV{dcaPos > v0Selections.dcaPionFromLambda && dcaNeg > v0Selections.dcaProtonFromLambda ? true : false};

    // Rejects V0 if the invariant mass is not compatible with the Lambda proper mass
    if (v0Selections.applyInvMassSel) {
      if (!(dMassAL < v0Selections.dMassSel && dMassK0s > v0Selections.dMassSel && dMassG > v0Selections.dMassSelG))
        return false;
    }

    bool isSelected{false};
    if (v0Selections.useOfficialV0sSelOfDaughters) {
      isSelected = lifeTime < v0Selections.lifeTimeCutLambda && rapidity < v0Selections.rapidityCut && dcaDaugToPV ? true : false;
    } else {
      isSelected = lifeTime < v0Selections.lifeTimeCutLambda && rapidity < v0Selections.rapidityCut && posTPCNsigma < v0Selections.tpcPidNsigmaCut && negTPCNsigma < v0Selections.tpcPidNsigmaCut && dcaDaugToPV ? true : false;
    }

    if (isSelected) {
      registry.fill(HIST("DCAxyPtPrAL"), v0.dcanegtopv(), negTrack.pt());
      registry.fill(HIST("EtaVsYPiAL"), posTrack.eta(), v0.yLambda());
      registry.fill(HIST("EtaVsYPrAL"), negTrack.eta(), v0.yLambda());
    }

    return isSelected;
  }

  template <typename C, typename T>
  bool passesGammaSelection(const C& /*collision*/, const T& v0)
  {
    const auto& posTrack = v0.template posTrack_as<TracksFull>();
    const auto& negTrack = v0.template negTrack_as<TracksFull>();

    const double dcaPos{std::fabs(v0.dcapostopv())};
    const double dcaNeg{std::fabs(v0.dcanegtopv())};
    const float posTPCNsigma{std::fabs(posTrack.tpcNSigmaEl())};
    const float negTPCNsigma{std::fabs(negTrack.tpcNSigmaEl())};
    const double dMassK0s{std::abs(v0.mK0Short() - o2::constants::physics::MassK0Short)};
    const double dMassL{std::abs(v0.mLambda() - o2::constants::physics::MassLambda0)};
    const double dMassAL{std::abs(v0.mAntiLambda() - o2::constants::physics::MassLambda0)};
    const double dMassG{std::abs(v0.mGamma() - o2::constants::physics::MassGamma)};
    const float yGamma = RecoDecay::y(std::array{v0.px(), v0.py(), v0.pz()}, o2::constants::physics::MassGamma);

    // Checks if DCA of daughters to PV passes selection
    const bool dcaDaugToPV{dcaPos > v0Selections.dcaElectronFromGamma && dcaNeg > v0Selections.dcaElectronFromGamma ? true : false};

    if (v0Selections.applyInvMassSel) {
      if (!(dMassK0s > v0Selections.dMassSel && dMassL > v0Selections.dMassSel && dMassAL > v0Selections.dMassSel && dMassG < v0Selections.dMassSel))
        return false;
    }

    if (!(std::abs(yGamma) < v0Selections.rapidityCut))
      return false;

    bool isSelected{false};
    if (v0Selections.useOfficialV0sSelOfDaughters)
      isSelected = dcaDaugToPV ? true : false;
    else
      isSelected = posTPCNsigma < v0Selections.tpcPidNsigmaCut && negTPCNsigma < v0Selections.tpcPidNsigmaCut ? true : false;

    if (isSelected) {
      registry.fill(HIST("EtaVsYG"), negTrack.eta(), yGamma);
      registry.fill(HIST("EtaVsYG"), posTrack.eta(), yGamma);
    }

    return isSelected;
  }

  int getMagneticField(uint64_t timestamp)
  {
    // TODO done only once (and not per run). Will be replaced by CCDBConfigurable
    static o2::parameters::GRPMagField* grpo = nullptr;
    if (grpo == nullptr) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 0;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
    }
    return grpo->getNominalL3Field();
  }

  void phiPrimeFunc(float& phi, const int& magField, const int& charge)
  {

    if (magField < 0)
      phi = o2::constants::math::TwoPI - phi;
    if (charge < 0)
      phi = o2::constants::math::TwoPI - phi;

    phi += o2::constants::math::PI / 18.0f;
    phi = std::fmod(phi, o2::constants::math::PI / 9.0f);
  }

  bool passesPhiSelection(const float& pt, const float& phi)
  {

    // Do not apply Phi Sel if pt < 2 GeV/c
    if (pt < kTwoPtGeVSel)
      return true;

    bool isSelected{true};
    if (phiCut.isPhiCutLoaded) {
      const int binLow{phiCut.hPhiCutLow->FindBin(pt)};
      const int binHigh{phiCut.hPhiCutHigh->FindBin(pt)};
      const double phiCutLow{phiCut.hPhiCutLow->GetBinContent(binLow)};
      const double phiCutHigh{phiCut.hPhiCutHigh->GetBinContent(binHigh)};
      if (phi >= phiCutLow && phi <= phiCutHigh)
        isSelected = false;
    }
    return isSelected;
  }

  template <typename CheckCol>
  bool isEventSelected(CheckCol const& col)
  {
    registry.fill(HIST("EventCounter"), EvCutLabel::All);
    if (!col.sel8()) {
      return false;
    }
    registry.fill(HIST("EventCounter"), EvCutLabel::SelEigth);

    if (selNoSameBunchPileup) {
      if (!col.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
        return false;
      }
      registry.fill(HIST("EventCounter"), EvCutLabel::NoSameBunchPileup);
    }

    if (selIsGoodZvtxFT0vsPV) {
      if (!col.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
        return false;
      }
      registry.fill(HIST("EventCounter"), EvCutLabel::IsGoodZvtxFT0vsPV);
    }

    if (isNoCollInTimeRangeStrict) {
      if (!col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict)) {
        return false;
      }
      registry.fill(HIST("EventCounter"), EvCutLabel::NoCollInTimeRangeStrict);
    }

    if (isNoCollInTimeRangeStandard) {
      if (!col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
        return false;
      }
      registry.fill(HIST("EventCounter"), EvCutLabel::NoCollInTimeRangeStandard);
    }

    if (isNoCollInRofStrict) {
      if (!col.selection_bit(o2::aod::evsel::kNoCollInRofStrict)) {
        return false;
      }
      registry.fill(HIST("EventCounter"), EvCutLabel::NoCollInRofStrict);
    }

    if (isNoCollInRofStandard) {
      if (!col.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
        return false;
      }
      registry.fill(HIST("EventCounter"), EvCutLabel::NoCollInRofStandard);
    }

    if (isNoHighMultCollInPrevRof) {
      if (!col.selection_bit(
            o2::aod::evsel::kNoHighMultCollInPrevRof)) {
        return false;
      }
      registry.fill(HIST("EventCounter"), EvCutLabel::NoHighMultCollInPrevRof);
    }

    // To be used in combination with FT0C-based occupancy
    if (isNoCollInTimeRangeNarrow) {
      if (!col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow)) {
        return false;
      }
      registry.fill(HIST("EventCounter"), EvCutLabel::NoCollInTimeRangeNarrow);
    }

    if (isOccupancyCut) {
      auto occuValue{isApplyFT0CbasedOccupancy ? col.ft0cOccupancyInTimeRange() : col.trackOccupancyInTimeRange()};
      if (occuValue < minOccCut || occuValue > maxOccCut) {
        return false;
      }
      registry.fill(HIST("EventCounter"), EvCutLabel::OccuCut);
    }

    if (selHasFT0) {
      if (!col.has_foundFT0()) {
        return false;
      }
      registry.fill(HIST("EventCounter"), EvCutLabel::HasFT0);
    }

    if (isCentSel) {
      if (col.centFT0C() < minT0CcentCut || col.centFT0C() > maxT0CcentCut) {
        return false;
      }
      registry.fill(HIST("EventCounter"), EvCutLabel::Centrality);
    }

    if (isZvtxPosSel) {
      if (std::fabs(col.posZ()) > posZcut) {
        return false;
      }
      registry.fill(HIST("EventCounter"), EvCutLabel::VtxZ);
    }

    if (selINELgt0) {
      if (!col.isInelGt0()) {
        return false;
      }
      registry.fill(HIST("EventCounter"), EvCutLabel::INELgt0);
    }

    return true;
  }

  template <typename T, typename U>
  void getPTpowers(const T& pTs, const T& vecEff, const T& vecFD, U& pOne,
                   U& wOne, U& pTwo, U& wTwo, U& pThree, U& wThree,
                   U& pFour, U& wFour)
  {
    pOne = wOne = pTwo = wTwo = pThree = wThree = pFour = wFour = 0.;
    for (std::size_t i = 0; i < pTs.size(); ++i) {
      const double pTi{pTs.at(i)};
      const double eFFi{vecEff.at(i)};
      const double fDi{vecFD.at(i)};
      const double wEighti{std::pow(eFFi, -1.) * fDi};
      pOne += wEighti * pTi;
      wOne += wEighti;
      pTwo += std::pow(wEighti * pTi, 2.);
      wTwo += std::pow(wEighti, 2.);
      pThree += std::pow(wEighti * pTi, 3.);
      wThree += std::pow(wEighti, 3.);
      pFour += std::pow(wEighti * pTi, 4.);
      wFour += std::pow(wEighti, 4.);
    }
  }

  void loadNchCalibrations(uint64_t timeStamp)
  {
    if (pathMeanNch.value.empty() == false) {
      cfgNch.hMeanNch = ccdb->getForTimeStamp<TH1F>(pathMeanNch, timeStamp);
      if (cfgNch.hMeanNch == nullptr) {
        LOGF(fatal, "Could not load hMeanNch histogram from %s", pathMeanNch.value.c_str());
      }
    }

    if (pathSigmaNch.value.empty() == false) {
      cfgNch.hSigmaNch = ccdb->getForTimeStamp<TH1F>(pathSigmaNch, timeStamp);
      if (cfgNch.hSigmaNch == nullptr) {
        LOGF(fatal, "Could not load hSigmaNch histogram from %s", pathSigmaNch.value.c_str());
      }
    }
    if (cfgNch.hMeanNch && cfgNch.hSigmaNch)
      cfgNch.calibrationsLoaded = true;
  }

  void loadPhiCutSelections(const uint64_t& timeStamp)
  {

    if (pathPhiCutHigh.value.empty() == false) {
      phiCut.hPhiCutHigh = ccdb->getForTimeStamp<TH1F>(pathPhiCutHigh, timeStamp);
      if (phiCut.hPhiCutHigh == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram from %s", pathPhiCutHigh.value.c_str());
      }
    }

    if (pathPhiCutLow.value.empty() == false) {
      phiCut.hPhiCutLow = ccdb->getForTimeStamp<TH1F>(pathPhiCutLow, timeStamp);
      if (phiCut.hPhiCutLow == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram from %s", pathPhiCutLow.value.c_str());
      }
    }

    if (phiCut.hPhiCutHigh && phiCut.hPhiCutLow)
      phiCut.isPhiCutLoaded = true;
  }

  void loadEtaCalibration()
  {
    if (pathEtaCal.value.empty() == false) {
      etaCal.pEtaCal = ccdb->getForTimeStamp<TProfile>(pathEtaCal, ccdbNoLaterThan.value);
      if (etaCal.pEtaCal == nullptr)
        LOGF(fatal, "Could not load pEtaCal from %s", pathEtaCal.value.c_str());
    }

    if (etaCal.pEtaCal)
      etaCal.isMIPCalLoaded = true;
  }

  void loadEtaPlateauCalibration()
  {
    if (pathEtaCalPlateau.value.empty() == false) {
      etaCal.pEtaCalPlateau = ccdb->getForTimeStamp<TProfile>(pathEtaCalPlateau, ccdbNoLaterThan.value);

      if (etaCal.pEtaCalPlateau == nullptr)
        LOGF(fatal, "Could not load pEtaCalPlateau from %s", pathEtaCalPlateau.value.c_str());
    }

    if (etaCal.pEtaCalPlateau)
      etaCal.isCalPlateauLoaded = true;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PiKpRAA>(cfgc)};
}
