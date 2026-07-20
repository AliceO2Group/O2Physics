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
#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TMCProcess.h>
#include <TPDGCode.h>
#include <TProfile.h>
#include <TString.h>
#include <TVector3.h>

#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::evsel;
using namespace o2::constants::math;
using namespace o2::framework::expressions;
using namespace o2::aod::rctsel;

using ColEvSels = soa::Join<aod::Collisions, aod::EvSels, aod::FT0MultZeqs, o2::aod::CentFT0Cs, o2::aod::CentFT0Ms, o2::aod::CentFV0As, aod::TPCMults, o2::aod::BarrelMults>;
using BCsRun3 = soa::Join<aod::BCsWithTimestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;

using ColEvSelsMC = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, o2::aod::CentFT0Cs, o2::aod::CentFT0Ms, o2::aod::CentFV0As, o2::aod::CentFT0Ms, o2::aod::BarrelMults>;

using TracksFull = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelectionExtension, aod::TracksDCA, aod::TrackSelection, aod::TracksCovIU, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTPCEl, aod::pidTOFFlags, aod::pidTOFbeta, aod::TOFSignal, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullEl>;

using TracksMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelectionExtension, aod::TracksDCA, aod::TrackSelection, aod::TracksCovIU, aod::pidTPCPi, aod::pidTPCPr, aod::pidTOFPr, aod::pidTPCEl, aod::pidTOFFlags, aod::pidTOFbeta, aod::TOFSignal, aod::pidTOFFullPi, aod::pidTOFFullEl, aod::McTrackLabels>;

static constexpr int KnEtaHists{8};

std::array<std::shared_ptr<TH2>, KnEtaHists> dEdxPiTOF{};
std::array<std::shared_ptr<TH2>, KnEtaHists> dEdxPiTOF2{};
std::array<std::shared_ptr<TH2>, KnEtaHists> dEdxPiV0{};
std::array<std::shared_ptr<TH2>, KnEtaHists> dEdxPrV0{};
std::array<std::shared_ptr<TH2>, KnEtaHists> dEdxElV0{};
std::array<std::shared_ptr<TH3>, KnEtaHists> dEdx{};
std::array<std::shared_ptr<TH2>, KnEtaHists> pTVsP{};
std::array<std::shared_ptr<TH2>, KnEtaHists> nClVsP{};
std::array<std::shared_ptr<TH2>, KnEtaHists> nClVsPElV0{};
std::array<std::shared_ptr<TH2>, KnEtaHists> nClVsPPiV0{};
std::array<std::shared_ptr<TH2>, KnEtaHists> nClVsPPrV0{};
std::array<std::shared_ptr<TProfile>, KnEtaHists> nClVsPp{};
std::array<std::shared_ptr<TProfile>, KnEtaHists> nClVsPpElV0{};
std::array<std::shared_ptr<TProfile>, KnEtaHists> nClVsPpPiV0{};
std::array<std::shared_ptr<TProfile>, KnEtaHists> nClVsPpPrV0{};
std::array<std::shared_ptr<TH1>, KnEtaHists> etaTest{};

template <typename T>
auto printArray(std::vector<T> const& vec)
{
  std::stringstream ss;
  ss << "[" << vec[0];
  for (auto i = 1u; i < vec.size(); ++i) {
    ss << ", " << vec[i];
  }
  ss << "]";
  return ss.str();
}

struct PiKpRAA {

  static constexpr int KzeroInt{0};
  static constexpr int KsevenInt{7};

  static constexpr float Kzero{0.0f};
  static constexpr float Kone{1.0f};
  static constexpr float KtwoPtGeVSel{2.0f};
  static constexpr float Kthree{3.0f};
  static constexpr float KtEnToMinusNine{1e-9};
  static constexpr float KminCharge{3.f};
  static constexpr float KminPMIP{0.4f};
  static constexpr float KmaxPMIP{0.6f};
  static constexpr float KmindEdxMIP{40.0f};
  static constexpr float KmaxdEdxMIP{60.0f};
  static constexpr float KmindEdxMIPPlateau{70.0f};
  static constexpr float KmaxdEdxMIPPlateau{90.0f};
  static constexpr float KminFT0A{3.5f};
  static constexpr float KmaxFT0A{4.9f};
  static constexpr float KminFT0C{-3.3f};
  static constexpr float KmaxFT0C{-2.1f};

  static constexpr float DefaultLifetimeCuts[1][2] = {{30., 20.}};
  Configurable<LabeledArray<float>> lifetimecut{"lifetimecut", {DefaultLifetimeCuts[0], 2, {"lifetimecutLambda", "lifetimecutK0S"}}, "lifetimecut"};

  struct : ConfigurableGroup {
    Configurable<float> minEta{"minEta", -0.8, "minimum eta selection"};
    Configurable<float> maxEta{"maxEta", +0.8, "maximum eta selection"};
    Configurable<float> minPt{"minPt", 0.1, "minimum pt of the tracks"};
    Configurable<float> maxPt{"maxPt", 10000000000.0, "maximum pt of the tracks"};
    Configurable<float> nSigmaDCAxy{"nSigmaDCAxy", 1.0, "nSigma DCAxy selection"};
    Configurable<float> nSigmaDCAz{"nSigmaDCAz", 1.0, "nSigma DCAz selection"};

    Configurable<int16_t> minNCrossedRows{"minNCrossedRows", 70, "minimum number of crossed rows"};
    Configurable<int16_t> minNcl{"minNcl", 135, "minimum number of clusters for primary-particle selection"};
    Configurable<uint8_t> minNClusITS{"minNClusITS", 7, "minimum number of ITS clusters"};
    Configurable<float> minChi2ClsTPC{"minChi2ClsTPC", 0.5, "Max chi2 per Cls TPC"};
    Configurable<float> maxChi2ClsTPC{"maxChi2ClsTPC", 4.0, "Max chi2 per Cls TPC"};
    Configurable<float> chi2ClsITS{"chi2ClsITS", 36.0, "chi2 per Cls ITS selection"};
    Configurable<float> maxElTOFBeta{"maxElTOFBeta", 0.1, "Maximum beta TOF selection"};
    Configurable<float> maxPiTOFBeta{"maxPiTOFBeta", 0.005, "Maximum beta TOF selection for Pions"};
    Configurable<float> nSigmaPIDselPrim{"nSigmaPIDselPrim", 1.0, "N sigma selection for primary pions with TOF"};

    // Phi cut
    Configurable<bool> applyPhiCut{"applyPhiCut", false, "Apply geometrical cut?"};
    Configurable<bool> applyEtaCal{"applyEtaCal", false, "Apply eta calibration?"};
    Configurable<bool> applyPlateauSel{"applyPlateauSel", false, "Apply eta calibration?"};
    Configurable<bool> usePinPhiSelection{"usePinPhiSelection", true, "Uses Phi selection as a function of p or p at inner wall of TPC?"};
    Configurable<bool> applyNclSel{"applyNclSel", false, "Apply Min. N clusters selection?"};
    Configurable<bool> useNclsPID{"useNclsPID", true, "Use Ncl for PID?"};
    Configurable<std::string> signCharge{"signCharge", "Positive", "Perform analysis as a function of the particle charge"};
  } trackSelections;

  struct : ConfigurableGroup {
    Configurable<uint8_t> v0TypeSelection{"v0TypeSelection", 1, "select on a certain V0 type (leave negative if no selection desired)"};
    Configurable<bool> useOfficialV0sSelOfDaughters{"useOfficialV0sSelOfDaughters", true, "Use the same track selection for daughters as the V0s analysis in OO"};
    Configurable<bool> useTPCNsigma{"useTPCNsigma", false, "Include TPC Nsigma selection on decay products"};

    // Selection criteria: acceptance
    Configurable<float> minY{"minY", -0.1, "minimum rapidity selection"};
    Configurable<float> maxY{"maxY", 0.8, "maximum rapidity selection"};
    Configurable<int16_t> minNclV0Daugther{"minNclV0Daugther", 135, "minimum number of clusters for V0 daughters"};
    Configurable<uint8_t> nSharedClusTpc{"nSharedClusTpc", 5, "maximum number of shared clusters"};

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
    Configurable<float> dMassSel{"dMassSel", 0.01f, "Invariant mass selection"};
    Configurable<float> dMassSelG{"dMassSelG", 0.1f, "Inv mass selection gammas"};

    // PID (TPC/TOF)
    Configurable<float> dEdxPlateauSel{"dEdxPlateauSel", 50, "dEdx selection for electrons"};
    Configurable<float> pidNsigmaCut{"pidNsigmaCut", 3.0, "N sigma selection for V0 daughters"};
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
  Configurable<bool> selHasBC{"selHasBC", true, "Has BC?"};
  Configurable<bool> selHasFT0{"selHasFT0", true, "Has FT0?"};
  Configurable<std::string> centralitySelector{"centralitySelector", "FT0C", "which centrality selector?"};

  Configurable<bool> useSel8{"useSel8", false, "Use sel8?"};
  Configurable<bool> selTriggerTVX{"selTriggerTVX", true, "selTriggerTVX?"};
  Configurable<bool> selNoITSROFrameBorder{"selNoITSROFrameBorder", true, "selNoITSROFrameBorder?"};
  Configurable<bool> selNoTimeFrameBorder{"selNoTimeFrameBorder", true, "selNoTimeFrameBorder?"};
  Configurable<bool> isZvtxPosSel{"isZvtxPosSel", true, "Zvtx position selection?"};
  Configurable<bool> isZvtxPosSelMC{"isZvtxPosSelMC", true, "Zvtx position selection for MC events?"};
  Configurable<bool> selTVXMC{"selTVXMC", true, "apply TVX selection in MC?"};
  Configurable<bool> selINELgt0{"selINELgt0", true, "Select INEL > 0?"};
  Configurable<bool> isApplyFT0CbasedOccupancy{"isApplyFT0CbasedOccupancy", false, "T0C Occu cut"};
  Configurable<bool> loadHisWithDCASel{"loadHisWithDCASel", false, "load histograms with DCA selections"};

  // Event selection
  Configurable<float> posZcut{"posZcut", +10.0, "z-vertex position cut"};
  Configurable<float> minCentCut{"minCentCut", 0.0, "Min Cent. cut"};
  Configurable<float> maxCentCut{"maxCentCut", 100.0, "Max Cent. cut"};
  Configurable<float> minOccCut{"minOccCut", 0., "min Occu cut"};
  Configurable<float> maxOccCut{"maxOccCut", 500., "max Occu cut"};
  Configurable<float> tpcNchAcceptance{"tpcNchAcceptance", 0.5, "Eta window to measure Nch MC for Nch vs Cent distribution"};

  ConfigurableAxis binsPtV0s{"binsPtV0s", {VARIABLE_WIDTH, 0, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.5, 3.0, 3.5, 4, 5, 7, 9, 12, 15, 20}, "pT"};
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.12}, "pT binning"};
  ConfigurableAxis binsCent{"binsCent", {VARIABLE_WIDTH, 0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.}, "Centrality bins"};
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
  ConfigurableAxis axisPtFineFixedWidth{"axisPtFineFixedWidth", {250, 0.1f, 20.1f}, "DCAxy axis"};
  ConfigurableAxis axisNPV{"axisNPV", {600, 0.0f, 600.0f}, "Number of Primary Vertex Contributors"};
  ConfigurableAxis axisNch{"axisNch", {400, 0.0f, 400.0f}, "N_{ch} Multiplicity"};

  Configurable<std::vector<double>> vecLowEta{"vecLowEta", {-0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6}, "lower eta cuts"};
  Configurable<std::vector<double>> vecUpEta{"vecUpEta", {-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8}, "upper eta cuts"};
  Configurable<std::vector<std::string>> vecEndingEta{"vecEndingEta", {"86", "64", "42", "20", "02", "24", "46", "68"}, "string representing names for histogram distiction"};
  Configurable<std::vector<std::string>> vecLatexEta{"vecLatexEta", {"-0.8<#eta<-0.6", "-0.6<#eta<-0.4", "-0.4<#eta<-0.2", "-0.2<#eta< 0", "0<#eta<0.2", "0.2<#eta<0.4", "0.4<#eta< 0.6", "0.6<#eta< 0.8"}, "string representing an eta interval"};

  // CCDB paths
  Configurable<std::string> pathDCAxy{"pathDCAxy", "Users/o/omvazque/DCAxy/Test", "base path to the ccdb object"};
  Configurable<std::string> pathDCAz{"pathDCAz", "Users/o/omvazque/DCAz/Test", "base path to the ccdb object"};
  Configurable<std::string> pathEtaCal{"pathEtaCal", "Users/o/omvazque/EtaCal/OO/Global", "base path to the ccdb object"};
  Configurable<std::string> pathEtaCalPlateau{"pathEtaCalPlateau", "Users/o/omvazque/EtaCal/OO/Global", "base path to the ccdb object"};
  Configurable<std::string> pathPhiCutHigh{"pathPhiCutHigh", "Users/o/omvazque/PhiCut/OO/Global/High", "base path to the ccdb object"};
  Configurable<std::string> pathPhiCutLow{"pathPhiCutLow", "Users/o/omvazque/PhiCut/OO/Global/Low", "base path to the ccdb object"};
  Configurable<int64_t> ccdbNoLaterThan{"ccdbNoLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> rctLabel{"rctLabel", "CBT_hadronPID", "RCT selection flag (CBT, CBT_hadronPID, CBT_electronPID, CBT_calo, CBT_muon, CBT_muon_glo)"};
  Configurable<bool> rctCheckZDC{"rctCheckZDC", false, "RCT flag to check whether the ZDC is present or not"};
  Configurable<bool> rctTreatLimitedAcceptanceAsBad{"rctTreatLimitedAcceptanceAsBad", false, "RCT flag to reject events with limited acceptance for selected detectors"};
  Configurable<bool> requireGoodRct{"requireGoodRct", true, "RCT flag to reject events with limited acceptance for selected detectors"};
  Configurable<bool> requireBCRct{"requireBCRct", true, "RCT flag to reject events with limited acceptance for selected detectors"};

  // RCT Checker instance
  RCTFlagsChecker rctChecker;

  enum EvCutLabel {
    All = 1,
    HasBC,
    HasFT0,
    SelEigth,
    SelTriggerTVX,
    SelNoITSROFrameBorder,
    SelNoTimeFrameBorder,
    VtxZ,
    IsGoodZvtxFT0vsPV,
    NoSameBunchPileup,
    NoCollInTimeRangeStrict,
    NoCollInTimeRangeStandard,
    NoCollInRofStrict,
    NoCollInRofStandard,
    NoHighMultCollInPrevRof,
    NoCollInTimeRangeNarrow,
    OccuCut,
    Centrality,
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
  Service<ccdb::BasicCCDBManager> ccdb{};

  struct ConfigDCA {
    TH1F* hDCAxy = nullptr;
    TH1F* hDCAz = nullptr;
    bool dcaSelectionsLoaded = false;
  } cfgDCA;

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

  int currentRunNumberPhiSel = -1;
  void init(InitContext const&)
  {

    // Initialize the rct checker
    if (requireGoodRct) {
      rctChecker.init(rctLabel.value, rctCheckZDC.value, rctTreatLimitedAcceptanceAsBad.value);
    }

    currentRunNumberPhiSel = -1;

    // define axes you want to use
    const std::string titlePorPTPC{trackSelections.usePinPhiSelection ? "#it{p} (GeV/#it{c})" : "#it{p} inner TPC wall (GeV/#it{c})"};
    const AxisSpec axisZpos{binsZpos, "Vtx_{z} (cm)"};
    const AxisSpec axisEvent{22, 0.5, 22.5, ""};
    const AxisSpec axisNcl{161, -0.5, 160.5, "#it{N}_{cl}"};
    const AxisSpec axisPt{binsPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisPtV0s{binsPtV0s, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisCent{binsCent, "Centrality Perc."};
    const auto latexEta = (std::vector<std::string>)vecLatexEta;
    const auto endingEta = (std::vector<std::string>)vecEndingEta;

    registry.add("EventCounter", ";;Events", kTH1F, {axisEvent});
    registry.add("HasBCVsFT0VsTVXVsEvSel", "Alls=1 | BC=2 | FT0=3 | TVX=4 | EvSel=5;;", kTH1F, {{5, 0.5, 5.5}});
    registry.add("zPos", "With Event Selection;;Entries;", kTH1F, {axisZpos});
    registry.add("Centrality", ";;Entries", kTH1F, {axisCent});
    registry.add("RCTSel", "Event accepted if flag=false: All=1 | RTC sel=2;;;", kTH1F, {{2, 0.5, 2.5}});
    registry.add("CentralityVsRCTSel", "Event accepted if flag=false;;RCT Status;", kTH2F, {{{axisCent}, {9, 0.5, 9.5}}});
    registry.add("CentralityVsFoundFT0", "Found(=1.5) NOT Found(=0.5);;Status;", kTH2F, {{{axisCent}, {2, 0, 2}}});
    registry.add("NchVsCent", "Measured Nch v.s. Centrality (At least Once Rec. Coll. + Sel. criteria);;Nch", kTH2F, {axisCent, axisNch});
    registry.add("NclVsEtaPID", "", kTH2F, {axisEta, axisNcl});
    registry.add("NclVsEtaPIDp", "", kTProfile, {axisEta});
    registry.add("NclVsEta", "", kTH2F, {axisEta, axisNcl});
    registry.add("NclVsEtap", "", kTProfile, {axisEta});
    registry.add("NclVsPhiHist", "", kTH2F, {{100, 0, o2::constants::math::TwoPI}, {axisNcl}});
    registry.add("NclVsPhiProfile", "", kTProfile, {{100, 0, o2::constants::math::TwoPI}});
    registry.add("EtaVsPhi", ";#eta;#varphi;", kTH2F, {{axisEta}, {100, 0, o2::constants::math::TwoPI}});
    registry.add("MomentumTPCVsP", ";Global track momentum (GeV/#it{c});Momentum at inner wall of the TPC (GeV/#it{c});", kTH2F, {axisPtV0s, axisPtV0s});
    registry.add("DCAxyVsPt", ";#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm);", kTH2F, {{axisPtFineFixedWidth}, {axisDCAxy}});
    registry.add("DCAzVsPt", ";#it{p}_{T} (GeV/#it{c});DCA_{z} (cm);", kTH2F, {{axisPtFineFixedWidth}, {axisDCAxy}});
    registry.add("dcaVsPtPi", "Primary pions;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm);Centrality Perc.;", kTH3F, {axisPt, axisDCAxy, axisCent});
    registry.add("dcaVsPtPr", "Primary protons;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm);Centrality Perc.;", kTH3F, {axisPt, axisDCAxy, axisCent});

    auto hPvsPTPC = registry.get<TH2>(HIST("MomentumTPCVsP"));
    hPvsPTPC->GetXaxis()->SetTitle("p (GeV/#it{c})");
    hPvsPTPC->GetYaxis()->SetTitle("p at TPC inner wall (GeV/#it{c})");

    auto hstat = registry.get<TH1>(HIST("EventCounter"));
    auto* x = hstat->GetXaxis();
    x->SetBinLabel(1, "All");
    x->SetBinLabel(2, "Has BC?");
    x->SetBinLabel(3, "Has FT0?");
    x->SetBinLabel(4, "SelEigth");
    x->SetBinLabel(5, "SelTriggerTVX");
    x->SetBinLabel(6, "SelNoITSROFrameBorder");
    x->SetBinLabel(7, "SelNoTimeFrameBorder");
    x->SetBinLabel(8, "VtxZ Sel.");
    x->SetBinLabel(9, "GoodZvtxFT0vsPV");
    x->SetBinLabel(10, "NoSameBunchPileup");
    x->SetBinLabel(11, "NoCollInTimeRangeStrict");
    x->SetBinLabel(12, "NoCollInTimeRangeStandard");
    x->SetBinLabel(13, "NoCollInRofStrict");
    x->SetBinLabel(14, "NoCollInRofStandard");
    x->SetBinLabel(15, "NoHighMultCollInPrevRof");
    x->SetBinLabel(16, "NoCollInTimeRangeNarrow");
    x->SetBinLabel(17, "Occupancy Cut");
    x->SetBinLabel(18, "Cent. Sel.");
    x->SetBinLabel(19, "Nch Sel.");
    x->SetBinLabel(20, "INEL > 0");

    auto hrct = registry.get<TH2>(HIST("CentralityVsRCTSel"));
    auto* y = hrct->GetYaxis();
    y->SetBinLabel(1, "All");
    y->SetBinLabel(2, "kFT0Bad");
    y->SetBinLabel(3, "kITSBad");
    y->SetBinLabel(4, "kITSLimAccMCRepr");
    y->SetBinLabel(5, "kTOFBad");
    y->SetBinLabel(6, "kTOFLimAccMCRepr");
    y->SetBinLabel(7, "kTPCBadTracking");
    y->SetBinLabel(8, "kTPCBadPID");
    y->SetBinLabel(9, "kTPCLimAccMCRepr");

    if (doprocessCalibrationAndV0s) {
      registry.add("NchVsNPV", ";Nch; NPV;", kTH2F, {axisNPV, axisNch});
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
      registry.add("NclVsPhip", Form("Found #LTNcl#GT TPC;%s (GeV/#it{c});#varphi", titlePorPTPC.data()), kTProfile2D, {{{axisPt}, {350, 0.0, 0.35}}});
      registry.add("NclPIDVsPhip", Form("#LTNcl#GT used for PID;%s (GeV/#it{c});#varphi", titlePorPTPC.data()), kTProfile2D, {{{axisPt}, {350, 0.0, 0.35}}});
      registry.add("NclVsEtaPiMIP", "MIP #pi^{+} + #pi^{-}: 0.4 < #it{p} < 0.6 [GeV/#it{c}], 40 < dE/dx < 60;#eta;Ncl TPC", kTH2F, {{{axisEta}, {161, -0.5, 160.5}}});
      registry.add("NclVsEtaPiMIPp", "MIP #pi^{+} + #pi^{-}: 0.4 < #it{p} < 0.6 [GeV/#it{c}], 40 < dE/dx < 60;#eta;#LTNcl#GT TPC", kTProfile, {axisEta});
      registry.add("NclVsEtaPiV0", ";#eta;Ncl TPC", kTH2F, {axisEta, axisNcl});
      registry.add("NclVsEtaPiV0p", ";#eta;#LTNcl#GT TPC", kTProfile, {axisEta});
      registry.add("NclVsEtaPrV0", ";#eta;Ncl TPC", kTH2F, {axisEta, axisNcl});
      registry.add("NclVsEtaPrV0p", ";#eta;#LTNcl#GT TPC", kTProfile, {axisEta});
      registry.add("NclVsEtaElV0", ";#eta;Ncl TPC", kTH2F, {axisEta, axisNcl});
      registry.add("NclVsEtaElV0p", ";#eta;#LTNcl#GT TPC", kTProfile, {axisEta});
      registry.add("EtaVsYK0s", ";#eta;#it{y};", kTH2F, {axisEta, axisY});
      registry.add("EtaVsYPiL", ";#eta;#it{y};", kTH2F, {axisEta, axisY});
      registry.add("EtaVsYPrL", ";#eta;#it{y};", kTH2F, {axisEta, axisY});
      registry.add("EtaVsYPiAL", ";#eta;#it{y};", kTH2F, {axisEta, axisY});
      registry.add("EtaVsYPrAL", ";#eta;#it{y};", kTH2F, {axisEta, axisY});
      registry.add("EtaVsYG", ";#eta;#it{y};", kTH2F, {axisEta, axisY});
      registry.add("DCAxyPtPiK0s", ";DCA_{xy} (cm); p_{T} [GeV/#it{c}]", kTH2F, {axisDCAxy, axisPt});
      registry.add("DCAxyPtPrL", ";DCA_{xy} (cm); p_{T} [GeV/#it{c}]", kTH2F, {axisDCAxy, axisPt});
      registry.add("DCAxyPtPrAL", ";DCA_{xy} (cm); p_{T} [GeV/#it{c}]", kTH2F, {axisDCAxy, axisPt});
      registry.add("dEdxVsEtaPiMIP", "MIP #pi^{+} + #pi^{-}: 0.4 < #it{p} < 0.6 [GeV/#it{c}];#eta; dE/dx;", kTH2F, {{{axisEta}, {100, 0, 100}}});
      registry.add("dEdxVsEtaPiMIPp", "MIP #pi^{+} + #pi^{-}: 0.4 < #it{p} < 0.6 [GeV/#it{c}];#eta; #LTdE/dx#GT", kTProfile, {axisEta});
      registry.add("dEdxVsEtaElMIP", "MIP e^{+} + e^{-}: 0.4 < #it{p} < 0.6 [GeV/#it{c}];#eta; dE/dx;", kTH2F, {{{axisEta}, {100, 0, 100}}});
      registry.add("dEdxVsEtaElMIPp", "MIP e^{+} + e^{-}: 0.4 < #it{p} < 0.6 [GeV/#it{c}];#eta; #LTdE/dx#GT", kTProfile, {axisEta});
      registry.add("dEdxVsEtaPiMIPV0", "MIP #pi^{+} + #pi^{-}: 0.4 < #it{p} < 0.6 [GeV/#it{c}];#eta; dE/dx", kTH2F, {{{axisEta}, {100, 0, 100}}});
      registry.add("dEdxVsEtaPiMIPV0p", "MIP #pi^{+} + #pi^{-}: 0.4 < #it{p} < 0.6 [GeV/#it{c}];#eta; #LTdE/dx#GT", kTProfile, {axisEta});
      registry.add("dEdxVsEtaElMIPV0", "e^{+} + e^{-}: 0.4 <#it{p}_{T} < 0.6 [GeV/#it{c}];#eta; dE/dx", kTH2F, {{{axisEta}, {100, 0, 100}}});
      registry.add("dEdxVsEtaElMIPV0p", "e^{+} + e^{-}: 0.4 <#it{p}_{T} < 0.6 [GeV/#it{c}];#eta; #LTdE/dx#GT", kTProfile, {axisEta});
      registry.add("pTVsCent", "", kTH2F, {axisPt, axisCent});

      for (int i = 0; i < KnEtaHists; ++i) {
        dEdx[i] = registry.add<TH3>(Form("dEdx_%s", endingEta[i].c_str()), Form("%s;Momentum;dE/dx;", latexEta[i].c_str()), kTH3F, {axisPt, axisdEdx, axisCent});
        pTVsP[i] = registry.add<TH2>(Form("pTVsP_%s", endingEta[i].c_str()), Form("%s;Momentum;#it{p}_{T} (GeV/#it{c});", latexEta[i].c_str()), kTH2F, {axisPt, axisPt});
        dEdxPiTOF[i] = registry.add<TH2>(Form("dEdxPiTOF_%s", endingEta[i].c_str()), Form("#pi^{+} + #pi^{-}, %s;Momentum;dE/dx;", latexEta[i].c_str()), kTH2F, {axisPtV0s, axisdEdx});
        dEdxPiTOF2[i] = registry.add<TH2>(Form("dEdxPiTOF2_%s", endingEta[i].c_str()), Form("#pi^{+} + #pi^{-}, %s;Momentum;dE/dx;", latexEta[i].c_str()), kTH2F, {axisPtV0s, axisdEdx});
        dEdxPiV0[i] = registry.add<TH2>(Form("dEdxPiV0_%s", endingEta[i].c_str()), Form("#pi^{+} + #pi^{-}, %s;Momentum;dE/dx;", latexEta[i].c_str()), kTH2F, {axisPtV0s, axisdEdx});
        dEdxPrV0[i] = registry.add<TH2>(Form("dEdxPrV0_%s", endingEta[i].c_str()), Form("p + #bar{p}, %s;Momentum;dE/dx;", latexEta[i].c_str()), kTH2F, {axisPtV0s, axisdEdx});
        dEdxElV0[i] = registry.add<TH2>(Form("dEdxElV0_%s", endingEta[i].c_str()), Form("e^{+} + e^{-}, %s;Momentum;dE/dx;", latexEta[i].c_str()), kTH2F, {axisPtV0s, axisdEdx});
        nClVsP[i] = registry.add<TH2>(Form("NclVsPPrimaries_%s", endingEta[i].c_str()), Form("%s;;Ncl TPC", latexEta[i].c_str()), kTH2F, {axisPtV0s, axisNcl});
        nClVsPElV0[i] = registry.add<TH2>(Form("NclVsPElV0_%s", endingEta[i].c_str()), Form("%s;;Ncl TPC", latexEta[i].c_str()), kTH2F, {axisPtV0s, axisNcl});
        nClVsPPiV0[i] = registry.add<TH2>(Form("NclVsPPiV0_%s", endingEta[i].c_str()), Form("%s;;Ncl TPC", latexEta[i].c_str()), kTH2F, {axisPtV0s, axisNcl});
        nClVsPPrV0[i] = registry.add<TH2>(Form("NclVsPPrV0_%s", endingEta[i].c_str()), Form("%s;;Ncl TPC", latexEta[i].c_str()), kTH2F, {axisPtV0s, axisNcl});
        nClVsPp[i] = registry.add<TProfile>(Form("NclVsPrimariesp_%s", endingEta[i].c_str()), Form("%s;;#LT#it{N}_{cl}#GT TPC", latexEta[i].c_str()), kTProfile, {axisPtV0s});
        nClVsPpElV0[i] = registry.add<TProfile>(Form("NclVsPElV0p_%s", endingEta[i].c_str()), Form("%s;;#LT#it{N}_{cl}#GT TPC", latexEta[i].c_str()), kTProfile, {axisPtV0s});
        nClVsPpPiV0[i] = registry.add<TProfile>(Form("NclVsPPiV0p_%s", endingEta[i].c_str()), Form("%s;;#LT#it{N}_{cl}#GT TPC", latexEta[i].c_str()), kTProfile, {axisPtV0s});
        nClVsPpPrV0[i] = registry.add<TProfile>(Form("NclVsPPrV0p_%s", endingEta[i].c_str()), Form("%s;;#LT#it{N}_{cl}#GT TPC", latexEta[i].c_str()), kTProfile, {axisPtV0s});
        etaTest[i] = registry.add<TH1>(Form("etaTest_%s", endingEta[i].c_str()), ";Eta;Entries;", kTH1F, {axisEta});
      }
    }

    if (doprocessSim) {
      registry.add("zPosMC", "Generated Events With at least One Rec. Collision + Sel. criteria;;Entries;", kTH1F, {axisZpos});
      registry.add("dcaVsPtPiDec", "Secondary pions from decays;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm);Centrality Perc.;", kTH3F, {axisPt, axisDCAxy, axisCent});
      registry.add("dcaVsPtPrDec", "Secondary protons from decays;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm);Centrality Perc.;", kTH3F, {axisPt, axisDCAxy, axisCent});
      registry.add("dcaVsPtPiMat", "Secondary pions from material interactions;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm);Centrality Perc.;", kTH3F, {axisPt, axisDCAxy, axisCent});
      registry.add("dcaVsPtPrMat", "Secondary protons from material interactions;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm);Centrality Perc.;", kTH3F, {axisPt, axisDCAxy, axisCent});
      registry.add("DCAxyVsPtWithSelection", ";#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm);", kTH2F, {{axisPtFineFixedWidth}, {axisDCAxy}});
      registry.add("DCAzVsPtWithSelection", ";#it{p}_{T} (GeV/#it{c});DCA_{z} (cm);", kTH2F, {{axisPtFineFixedWidth}, {axisDCAxy}});

      registry.add("CentralityVsBCVsFT0VsTVXVsEvSel", "All=1 | BC=2 | FT0=3 | TVX=4 | EvSel=5;;Status;", kTH2F, {{axisCent}, {5, 0.5, 5.5}});

      // MC events passing the TVX requirement
      registry.add("NchMCcentVsTVX", ";Passed(=1.5) NOT Passed(=0.5);", kTH2F, {axisNch, {2, 0, 2}});

      registry.add("NumberOfRecoCollisions", "Number of times Gen. Coll.are reconstructed;N;Entries", kTH1F, {{10, -0.5, 9.5}});

      // Pt resolution
      registry.add("PtResolution", "p_{T} resolution;;(pt_{rec} - pt_{gen})/pt_{gen};", kTH2F, {axisPt, {100, -1.0, 1.0}});

      // Needed to calculate the numerator of the Acceptance X Efficiency
      registry.add("PtPiVsCent_WithRecoEvt", "Generated Events With at least One Rec. Collision + Sel. criteria;;;", kTH2F, {axisPt, axisCent});
      registry.add("PtKaVsCent_WithRecoEvt", "Generated Events With at least One Rec. Collision + Sel. criteria;;;", kTH2F, {axisPt, axisCent});
      registry.add("PtPrVsCent_WithRecoEvt", "Generated Events With at least One Rec. Collision + Sel. criteria;;;", kTH2F, {axisPt, axisCent});

      registry.add("PtGenPiVsCent_WithRecoEvt", "Generated Events With at least One Rec. Collision + Sel. criteria;;;", kTH2F, {axisPt, axisCent});
      registry.add("PtGenKaVsCent_WithRecoEvt", "Generated Events With at least One Rec. Collision + Sel. criteria;;;", kTH2F, {axisPt, axisCent});
      registry.add("PtGenPrVsCent_WithRecoEvt", "Generated Events With at least One Rec. Collision + Sel. criteria;;;", kTH2F, {axisPt, axisCent});

      // Needed to calculate the denominator of the Acceptance X Efficiency
      registry.add("PtPiVsCentMC_WithRecoEvt", "Generated Events With at least One Rec. Collision;;;", kTH2F, {axisPt, axisCent});
      registry.add("PtKaVsCentMC_WithRecoEvt", "Generated Events With at least One Rec. Collision;;;", kTH2F, {axisPt, axisCent});
      registry.add("PtPrVsCentMC_WithRecoEvt", "Generated Events With at least One Rec. Collision;;;", kTH2F, {axisPt, axisCent});

      // Needed for the Gen. Nch to Centrality conversion
      registry.add("NchMCVsCent", "Generated Nch v.s. Centrality (At least Once Rec. Coll. + Sel. criteria);;Gen. Nch MC", kTH2F, {axisCent, axisNch});

      // Needed to measure Event Loss
      registry.add("NchMC_WithRecoEvt", "Generated Nch of Evts With at least one Rec. Coll. + Sel. criteria;Gen. Nch MC;Entries", kTH1F, {axisNch});
      registry.add("NchMC_AllGen", "Generated Nch of All Gen. Evts.;Gen. Nch;Entries", kTH1F, {axisNch});

      // Needed to measure Event Splitting
      registry.add("Centrality_WRecoEvt", "Generated Events With at least One Rec. Collision And NO Sel. criteria;;Entries", kTH1F, {axisCent});
      registry.add("Centrality_WRecoEvtWSelCri", "Generated Events With at least One Rec. Collision + Sel. criteria;;Entries", kTH1F, {axisCent});
      registry.add("Centrality_AllRecoEvt", "Generated Events Irrespective of the number of times it was reconstructed + Evt. Selections;;Entries", kTH1F, {axisCent});

      // Needed to calculate the numerator of the Signal Loss correction
      registry.add("PtPiVsNchMC_WithRecoEvt", "Generated Events With at least One Rec. Collision;;Gen. Nch;", kTH2F, {axisPt, axisNch});
      registry.add("PtKaVsNchMC_WithRecoEvt", "Generated Events With at least One Rec. Collision;;Gen. Nch;", kTH2F, {axisPt, axisNch});
      registry.add("PtPrVsNchMC_WithRecoEvt", "Generated Events With at least One Rec. Collision;;Gen. Nch;", kTH2F, {axisPt, axisNch});

      // Needed to calculate the denominator of the Signal Loss correction
      registry.add("PtPiVsNchMC_AllGen", "All Generated Events;;Gen. Nch;", kTH2F, {axisPt, axisNch});
      registry.add("PtKaVsNchMC_AllGen", "All Generated Events;;Gen. Nch;", kTH2F, {axisPt, axisNch});
      registry.add("PtPrVsNchMC_AllGen", "All Generated Events;;Gen. Nch;", kTH2F, {axisPt, axisNch});

      registry.add("MCclosure_PtMCPiVsNchMC", "All Generated Events 4 MC closure;;Gen. Nch;", kTH2F, {axisPt, axisNch});
      registry.add("MCclosure_PtMCKaVsNchMC", "All Generated Events 4 MC closure;;Gen. Nch;", kTH2F, {axisPt, axisNch});
      registry.add("MCclosure_PtMCPrVsNchMC", "All Generated Events 4 MC closure;;Gen. Nch;", kTH2F, {axisPt, axisNch});

      registry.add("MCclosure_PtPiVsNchMC", "Gen Evts With at least one Rec. Coll. + Sel. criteria 4 MC closure;;Gen. Nch (|#eta|<0.8);", kTH2F, {axisPt, axisNch});
      registry.add("MCclosure_PtKaVsNchMC", "Gen Evts With at least one Rec. Coll. + Sel. criteria 4 MC closure;;Gen. Nch (|#eta|<0.8);", kTH2F, {axisPt, axisNch});
      registry.add("MCclosure_PtPrVsNchMC", "Gen Evts With at least one Rec. Coll. + Sel. criteria 4 MC closure;;Gen. Nch (|#eta|<0.8);", kTH2F, {axisPt, axisNch});
    }

    LOG(info) << "\trequireGoodRct=" << requireGoodRct.value;
    LOG(info) << "\tccdbNoLaterThan=" << ccdbNoLaterThan.value;
    LOG(info) << "\tselINELgt0=" << selINELgt0.value;
    LOG(info) << "\tcentralitySelector=" << centralitySelector.value;
    LOG(info) << "\tv0TypeSelection=" << static_cast<int>(v0Selections.v0TypeSelection);
    LOG(info) << "\tapplyInvMassSel=" << v0Selections.applyInvMassSel;
    LOG(info) << "\tsignCharge =" << trackSelections.signCharge.value;
    LOG(info) << "\tminPt=" << trackSelections.minPt;
    LOG(info) << "\tmaxPt=" << trackSelections.maxPt;
    LOG(info) << "\tuseNclsPID=" << trackSelections.useNclsPID;
    LOG(info) << "\tqTSel=" << v0Selections.qTSel;
    LOG(info) << "\tarmAlphaSel=" << v0Selections.armAlphaSel;
    LOG(info) << "\tapplyNclSel=" << trackSelections.applyNclSel;
    LOG(info) << "\tapplyPhiCut=" << trackSelections.applyPhiCut;
    LOG(info) << "\tusePinPhiSelection=" << trackSelections.usePinPhiSelection;
    LOG(info) << "\tcurrentRunNumberPhiSel=" << currentRunNumberPhiSel;
    LOG(info) << "\tpathDCAxy=" << pathDCAxy.value;
    LOG(info) << "\tpathDCAz=" << pathDCAz.value;

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    ccdb->setCreatedNotAfter(ccdbNoLaterThan.value);

    if (trackSelections.applyPhiCut) {
      LOG(info) << "\tLoading Phi cut!";
      LOG(info) << "\tpathPhiCutLow=" << pathPhiCutLow.value;
      LOG(info) << "\tpathPhiCutHigh=" << pathPhiCutHigh.value;
    }

    if (trackSelections.applyEtaCal) {
      LOG(info) << "\tLoading Eta Cal!";
      LOG(info) << "\tpathEtaCal=" << pathEtaCal.value;
      loadEtaCalibration();
      LOG(info) << "\tisMIPCalLoaded=" << etaCal.isMIPCalLoaded;
    }

    if (trackSelections.applyPlateauSel) {
      LOG(info) << "\tLoading Eta Plateau Cal!";
      LOG(info) << "\tpathEtaCalPlateau=" << pathEtaCalPlateau.value;
      loadEtaPlateauCalibration();
      LOG(info) << "\tisCalPlateauLoaded=" << etaCal.isCalPlateauLoaded;
    }

    if (trackSelections.applyNclSel) {
      LOG(info) << "\ttrackSelections.minNcl=" << trackSelections.minNcl;
      LOG(info) << "\tv0Selections.minNclV0Daugther=" << v0Selections.minNclV0Daugther;
    }

    if (loadHisWithDCASel) {
      LOG(info) << "\tUsing DCAxy and DCAz selections";
      LOG(info) << "\tpathDCAxy.value=" << pathDCAxy.value;
      LOG(info) << "\tpathDCAz.value=" << pathDCAz.value;
      loadDCAselections();
      LOG(info) << "\tdcaSelectionsLoaded=" << cfgDCA.dcaSelectionsLoaded;
    }

    auto vecLowEtaSel = (std::vector<double>)vecLowEta;
    auto vecUpEtaSel = (std::vector<double>)vecUpEta;
    LOGF(info, "vecLowEta: %s", printArray(vecLowEtaSel).c_str());
    LOGF(info, "vecUpEta: %s", printArray(vecUpEtaSel).c_str());
    LOGF(info, "vecEndingEta: %s", printArray(endingEta).c_str());
    LOGF(info, "vecLatexEta: %s", printArray(latexEta).c_str());
  }

  void processCalibrationAndV0s(ColEvSels::iterator const& collision, BCsRun3 const& /**/, aod::V0Datas const& v0s, aod::FV0As const& /**/, aod::FT0s const& /**/, TracksFull const& tracks)
  {
    // LOG(info) << " Collisions size: " << collisions.size() << "
    // Table's size: " << collisions.tableSize() << "\n";
    // LOG(info) << "Run number: " << foundBC.runNumber() << "\n";

    const auto& kLowEta = (std::vector<double>)vecLowEta;
    const auto& kHighEta = (std::vector<double>)vecUpEta;

    const auto& foundBC = collision.foundBC_as<BCsRun3>();
    const uint64_t timeStamp{foundBC.timestamp()};
    const int magField{getMagneticField(timeStamp)};
    const double nPV{collision.multNTracksPVeta1() / 1.};
    float centrality{-999.0};
    if (centralitySelector.value == "FT0C") {
      centrality = collision.centFT0C();
    } else if (centralitySelector.value == "FT0M") {
      centrality = collision.centFT0M();
    } else if (centralitySelector.value == "FV0A") {
      centrality = collision.centFV0A();
    } else {
      centrality = -999.0;
    }

    // Apply RCT selection?
    if (requireGoodRct) {
      // Checks if collisions passes RCT selection
      const bool isFT0Bad{requireBCRct ? foundBC.rct_bit(kFT0Bad) : collision.rct_bit(kFT0Bad)};
      const bool isITSBad{requireBCRct ? foundBC.rct_bit(kITSBad) : collision.rct_bit(kITSBad)};
      const bool isITSLimAcc{requireBCRct ? foundBC.rct_bit(kITSLimAccMCRepr) : collision.rct_bit(kITSLimAccMCRepr)};
      const bool isTOFBad{requireBCRct ? foundBC.rct_bit(kTOFBad) : collision.rct_bit(kTOFBad)};
      const bool isTOFLimAcc{requireBCRct ? foundBC.rct_bit(kTOFLimAccMCRepr) : collision.rct_bit(kTOFLimAccMCRepr)};
      const bool isTPCTrackingBad{requireBCRct ? foundBC.rct_bit(kTPCBadTracking) : collision.rct_bit(kTPCBadTracking)};
      const bool isTPCPIDBad{requireBCRct ? foundBC.rct_bit(kTPCBadPID) : collision.rct_bit(kTPCBadPID)};
      const bool isTPCLimAcc{requireBCRct ? foundBC.rct_bit(kTPCLimAccMCRepr) : collision.rct_bit(kTPCLimAccMCRepr)};

      registry.fill(HIST("CentralityVsRCTSel"), centrality, 1.0);
      if (!isFT0Bad) {
        registry.fill(HIST("CentralityVsRCTSel"), centrality, 2.0);
      }
      if (!isITSBad) {
        registry.fill(HIST("CentralityVsRCTSel"), centrality, 3.0);
      }
      if (!isITSLimAcc) {
        registry.fill(HIST("CentralityVsRCTSel"), centrality, 4.0);
      }
      if (!isTOFBad) {
        registry.fill(HIST("CentralityVsRCTSel"), centrality, 5.0);
      }
      if (!isTOFLimAcc) {
        registry.fill(HIST("CentralityVsRCTSel"), centrality, 6.0);
      }
      if (!isTPCTrackingBad) {
        registry.fill(HIST("CentralityVsRCTSel"), centrality, 7.0);
      }
      if (!isTPCPIDBad) {
        registry.fill(HIST("CentralityVsRCTSel"), centrality, 8.0);
      }
      if (!isTPCLimAcc) {
        registry.fill(HIST("CentralityVsRCTSel"), centrality, 9.0);
      }

      registry.fill(HIST("RCTSel"), 1.0);
      if (!rctChecker(collision)) {
        return;
      }

      registry.fill(HIST("RCTSel"), 2.0);
    }

    if (!isEventSelected(collision)) {
      return;
    }

    //---------------------------
    // Control histogram
    //---------------------------
    if (selHasFT0 && !collision.has_foundFT0()) {
      registry.fill(HIST("CentralityVsFoundFT0"), centrality, 0.5);
    }

    registry.fill(HIST("CentralityVsFoundFT0"), centrality, 1.5);

    //---------------------------
    // Load Phi selection
    //---------------------------
    if (trackSelections.applyPhiCut) {
      const int nextRunNumber{foundBC.runNumber()};
      if (currentRunNumberPhiSel != nextRunNumber) {
        loadPhiCutSelections(timeStamp);
        currentRunNumberPhiSel = nextRunNumber;
        LOG(info) << "\tcurrentRunNumberPhiSel= " << currentRunNumberPhiSel << " timeStamp = " << timeStamp;
      }

      // return if phi cut objects are nullptr
      if (!(phiCut.hPhiCutHigh && phiCut.hPhiCutLow)) {
        return;
      }
    }

    registry.fill(HIST("zPos"), collision.posZ());
    registry.fill(HIST("Centrality"), centrality);

    // ================
    // Track selection WITHOUT DCAxy and DCAz selections
    // ================
    int nch{0};
    for (const auto& track : tracks) {

      const bool applyDca{false};
      if (!selectPrimary(track, applyDca)) {
        continue;
      }

      const int charge{track.sign()};
      const float momentum{track.p()};
      const float pTPC{track.tpcInnerParam()};
      const float pt{track.pt()};
      const float phi{track.phi()};
      const float dcaXy{track.dcaXY()};
      const float pOrpTPC{trackSelections.usePinPhiSelection ? momentum : pTPC};

      // Reject track on the basis of its charge
      if (trackSelections.signCharge.value == "Positive" && charge < KzeroInt) {
        continue;
      }
      if (trackSelections.signCharge.value == "Negative" && charge > KzeroInt) {
        continue;
      }

      float phiPrime{phi};
      if (trackSelections.applyPhiCut) {
        phiPrimeFunc(phiPrime, magField, charge);
        if (!passesPhiSelection(pOrpTPC, phiPrime)) {
          continue;
        }
      }

      const float piTPCNsigma{std::fabs(track.tpcNSigmaPi())};
      const float prTPCNsigma{std::fabs(track.tpcNSigmaPr())};
      const float piTOFNsigma{std::fabs(track.tofNSigmaPi())};
      const float prTOFNsigma{std::fabs(track.tofNSigmaPr())};
      const double piRadiusNsigma{std::sqrt(std::pow(piTPCNsigma, 2.) + std::pow(piTOFNsigma, 2.))};
      const double prRadiusNsigma{std::sqrt(std::pow(prTPCNsigma, 2.) + std::pow(prTOFNsigma, 2.))};

      if (piRadiusNsigma < Kthree) {
        registry.fill(HIST("dcaVsPtPi"), pt, dcaXy, centrality);
      }
      if (prRadiusNsigma < Kthree) {
        registry.fill(HIST("dcaVsPtPr"), pt, dcaXy, centrality);
      }
    }

    // ================
    // Track selection with DCAxy and DCAz selections
    // ================
    for (const auto& track : tracks) {

      const float momentum{track.p()};
      const float pTPC{track.tpcInnerParam()};
      const float pt{track.pt()};
      const float phi{track.phi()};
      const float eta{track.eta()};
      float dedx{track.tpcSignal()};
      const int charge{track.sign()};
      const float pOrpTPC{trackSelections.usePinPhiSelection ? momentum : pTPC};
      const int16_t nclFound{track.tpcNClsFound()};
      const int16_t nclPID{track.tpcNClsPID()};
      const int16_t ncl{trackSelections.useNclsPID ? nclPID : nclFound};

      // Reject track on the basis of its charge
      if (trackSelections.signCharge.value == "Positive" && charge < KzeroInt) {
        continue;
      }
      if (trackSelections.signCharge.value == "Negative" && charge > KzeroInt) {
        continue;
      }

      // ================
      // DCAxy & DCAz Selections
      // ================
      const bool applyDca{true};
      if (!selectPrimary(track, applyDca)) {
        continue;
      }

      float phiPrime{phi};
      if (trackSelections.applyPhiCut) {
        phiPrimeFunc(phiPrime, magField, charge);
        if (!passesPhiSelection(pOrpTPC, phiPrime)) {
          continue;
        }
      }

      if (trackSelections.applyEtaCal && etaCal.isMIPCalLoaded) {
        const double dedxCal{etaCal.pEtaCal->GetBinContent(etaCal.pEtaCal->FindBin(eta))};
        if (dedxCal > KmindEdxMIP && dedxCal < KmaxdEdxMIP) {
          dedx *= (50.0 / dedxCal);
        } else {
          continue;
        }
      }

      int indexEta{-999};
      for (int i = 0; i < KnEtaHists; ++i) {
        if (eta >= kLowEta[i] && eta < kHighEta[i]) {
          indexEta = i;
          break;
        }
      }

      if (indexEta < KzeroInt || indexEta > KsevenInt) {
        continue;
      }

      if (momentum > KminPMIP && momentum < KmaxPMIP && dedx > KmindEdxMIP && dedx < KmaxdEdxMIP) {
        registry.fill(HIST("dEdxVsEtaPiMIP"), eta, dedx);
        registry.fill(HIST("dEdxVsEtaPiMIPp"), eta, dedx);
        registry.fill(HIST("NclVsEtaPiMIP"), eta, ncl);
        registry.fill(HIST("NclVsEtaPiMIPp"), eta, ncl);
      }

      if (momentum > KminPMIP && momentum < KmaxPMIP && dedx > KmindEdxMIPPlateau && dedx < KmaxdEdxMIPPlateau) {
        if (track.hasTOF() && track.goodTOFMatch()) {
          const float tTOF{track.tofSignal()};
          const float trkLength{track.length()};
          const float tExpElTOF{track.tofExpSignalEl(tTOF)};

          if ((std::abs((tExpElTOF / tTOF) - Kone) < trackSelections.maxElTOFBeta) && tTOF > Kzero && trkLength > Kzero) {
            registry.fill(HIST("dEdxVsEtaElMIP"), eta, dedx);
            registry.fill(HIST("dEdxVsEtaElMIPp"), eta, dedx);
          }
        }
      }

      if (track.hasTOF() && track.goodTOFMatch()) {
        const float tTOF{track.tofSignal()};
        const float trkLength{track.length()};
        const float tExpPiTOF{track.tofExpSignalPi(tTOF)};
        const float tOFNsigmaPi{std::fabs(track.tofNSigmaPi())};
        const float tPCNsigmaPi{std::fabs(track.tpcNSigmaPi())};

        if (tTOF > Kzero && trkLength > Kzero) {
          if (std::abs((tExpPiTOF / tTOF) - Kone) < trackSelections.maxPiTOFBeta) {
            dEdxPiTOF2[indexEta]->Fill(momentum, dedx);
          }

          if (std::sqrt(std::pow(tOFNsigmaPi, 2.0) + std::pow(tPCNsigmaPi, 2.0)) < trackSelections.nSigmaPIDselPrim) {
            dEdxPiTOF[indexEta]->Fill(momentum, dedx);
          }
        }
      }

      etaTest[indexEta]->Fill(eta);
      dEdx[indexEta]->Fill(momentum, dedx, centrality);
      pTVsP[indexEta]->Fill(momentum, pt);
      nClVsP[indexEta]->Fill(pOrpTPC, ncl);
      nClVsPp[indexEta]->Fill(pOrpTPC, ncl);
      registry.fill(HIST("DCAxyVsPt"), pt, track.dcaXY());
      registry.fill(HIST("DCAzVsPt"), pt, track.dcaZ());
      registry.fill(HIST("MomentumTPCVsP"), momentum, pTPC);
      registry.fill(HIST("pTVsCent"), pt, centrality);
      registry.fill(HIST("EtaVsPhi"), eta, track.phi());
      registry.fill(HIST("NclVsPhiHist"), track.phi(), nclFound);
      registry.fill(HIST("NclVsPhiProfile"), track.phi(), nclFound);
      registry.fill(HIST("NclVsEta"), eta, nclFound);
      registry.fill(HIST("NclVsEtap"), eta, nclFound);
      registry.fill(HIST("NclVsEtaPID"), eta, nclPID);
      registry.fill(HIST("NclVsEtaPIDp"), eta, nclPID);
      registry.fill(HIST("NclVsPhip"), pOrpTPC, phiPrime, nclFound);
      registry.fill(HIST("NclPIDVsPhip"), pOrpTPC, phiPrime, nclPID);
      nch++;
    }

    registry.fill(HIST("NchVsCent"), centrality, nch);
    registry.fill(HIST("NchVsNPV"), nPV, nch);

    for (const auto& v0 : v0s) {

      // Select V0 type
      if (v0.v0Type() != v0Selections.v0TypeSelection) {
        continue;
      }

      // Positive-(negative-)charged tracks (daughters)
      const auto& posTrack = v0.posTrack_as<TracksFull>();
      const auto& negTrack = v0.negTrack_as<TracksFull>();
      const int posTrackCharge{posTrack.sign()};
      const int negTrackCharge{negTrack.sign()};
      const float posTrkP{posTrack.p()};
      const float negTrkP{negTrack.p()};
      const float posTrkPt{posTrack.pt()};
      const float negTrkPt{negTrack.pt()};
      const float posTrkPTPC{posTrack.tpcInnerParam()};
      const float negTrkPTPC{negTrack.tpcInnerParam()};
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
      const int16_t posNcl{trackSelections.useNclsPID ? posNclPID : posNclFound};
      const int16_t negNcl{trackSelections.useNclsPID ? negNclPID : negNclFound};
      const float posPorPTPC{trackSelections.usePinPhiSelection ? posTrkP : posTrkPTPC};
      const float negPorPTPC{trackSelections.usePinPhiSelection ? negTrkP : negTrkPTPC};

      // Skip v0s with like-sig daughters
      if (posTrack.sign() == negTrack.sign()) {
        continue;
      }

      // Passes Geometrical (Phi) cut?
      if (trackSelections.applyPhiCut) {
        phiPrimeFunc(posTrackPhiPrime, magField, posTrackCharge);
        phiPrimeFunc(negTrackPhiPrime, magField, negTrackCharge);
        if (!(passesPhiSelection(posPorPTPC, posTrackPhiPrime) && passesPhiSelection(negPorPTPC, negTrackPhiPrime))) {
          continue;
        }
      }

      // Passes daughters track-selection?
      if (!(selectV0Daughter(posTrack) && selectV0Daughter(negTrack))) {
        continue;
      }

      // Eta calibration positive-charge track
      if (trackSelections.applyEtaCal && etaCal.isMIPCalLoaded) {
        const double dedxCal{etaCal.pEtaCal->GetBinContent(etaCal.pEtaCal->FindBin(posTrkEta))};
        if (dedxCal > KmindEdxMIP && dedxCal < KmaxdEdxMIP) {
          posTrkdEdx *= (50.0 / dedxCal);
        } else {
          continue;
        }
      }

      // Eta calibration negative-charge track
      if (trackSelections.applyEtaCal && etaCal.isMIPCalLoaded) {
        const double dedxCal{etaCal.pEtaCal->GetBinContent(etaCal.pEtaCal->FindBin(negTrkEta))};
        if (dedxCal > KmindEdxMIP && dedxCal < KmaxdEdxMIP) {
          negTrkdEdx *= (50.0 / dedxCal);
        } else {
          continue;
        }
      }

      const TVector3 ppos(posTrack.px(), posTrack.py(), posTrack.pz());
      const TVector3 pneg(negTrack.px(), negTrack.py(), negTrack.pz());
      double alpha{-999.0};
      double qT{-999.0};

      getArmeterosVariables(ppos, pneg, alpha, qT);
      registry.fill(HIST("ArmAll"), alpha, qT);

      bool passesTopoSel{false};
      // Passes V0 topological cuts?
      if (passesV0TopologicalSelection(v0)) {
        passesTopoSel = true;
      }

      if (!passesTopoSel) {
        continue;
      }

      int posIndexEta{-999};
      int negIndexEta{-999};
      for (int i = 0; i < KnEtaHists; ++i) {
        if (posTrkEta >= kLowEta[i] && posTrkEta < kHighEta[i]) {
          posIndexEta = i;
          break;
        }
      }

      for (int i = 0; i < KnEtaHists; ++i) {
        if (negTrkEta >= kLowEta[i] && negTrkEta < kHighEta[i]) {
          negIndexEta = i;
          break;
        }
      }

      if (posIndexEta < KzeroInt || posIndexEta > KsevenInt) {
        continue;
      }

      if (negIndexEta < KzeroInt || negIndexEta > KsevenInt) {
        continue;
      }

      registry.fill(HIST("ArmAfterTopoSel"), alpha, qT);
      registry.fill(HIST("dcaDauVsPt"), v0.pt(), v0.dcapostopv());
      registry.fill(HIST("dcaDauVsPt"), v0.pt(), v0.dcanegtopv());

      if (passesK0Selection(collision, v0)) {
        registry.fill(HIST("ArmK0NOSel"), alpha, qT);
        if (v0Selections.armPodCut * qT > std::abs(alpha)) { // Armenters selection
          registry.fill(HIST("V0sCounter"), V0sCounter::K0s);
          registry.fill(HIST("ArmK0"), alpha, qT);
          registry.fill(HIST("MassK0sVsPt"), v0.pt(), v0.mK0Short());

          if (trackSelections.signCharge.value == "Positive" && posTrackCharge > KzeroInt) {
            registry.fill(HIST("nSigPiFromK0s"), posTrkPt, posTrack.tpcNSigmaPi());
            registry.fill(HIST("NclVsEtaPiV0"), posTrkEta, posNcl);
            registry.fill(HIST("NclVsEtaPiV0p"), posTrkEta, posNcl);
            nClVsPPiV0[posIndexEta]->Fill(posPorPTPC, posNcl);
            nClVsPpPiV0[posIndexEta]->Fill(posPorPTPC, posNcl);
            dEdxPiV0[posIndexEta]->Fill(posPorPTPC, posTrkdEdx);
            if (posPorPTPC > KminPMIP && posPorPTPC < KmaxPMIP && posTrkdEdx > KmindEdxMIP && posTrkdEdx < KmaxdEdxMIP) {
              registry.fill(HIST("dEdxVsEtaPiMIPV0"), posTrkEta, posTrkdEdx);
              registry.fill(HIST("dEdxVsEtaPiMIPV0p"), posTrkEta, posTrkdEdx);
            }
          } else if (trackSelections.signCharge.value == "Negative" && negTrackCharge < KzeroInt) {
            registry.fill(HIST("nSigPiFromK0s"), negTrkPt, negTrack.tpcNSigmaPi());
            registry.fill(HIST("NclVsEtaPiV0"), negTrkEta, negNcl);
            registry.fill(HIST("NclVsEtaPiV0p"), negTrkEta, negNcl);
            nClVsPPiV0[negIndexEta]->Fill(negPorPTPC, negNcl);
            nClVsPpPiV0[negIndexEta]->Fill(negPorPTPC, negNcl);
            dEdxPiV0[negIndexEta]->Fill(negPorPTPC, negTrkdEdx);
            if (negPorPTPC > KminPMIP && negPorPTPC < KmaxPMIP && negTrkdEdx > KmindEdxMIP && negTrkdEdx < KmaxdEdxMIP) {
              registry.fill(HIST("dEdxVsEtaPiMIPV0"), negTrkEta, negTrkdEdx);
              registry.fill(HIST("dEdxVsEtaPiMIPV0p"), negTrkEta, negTrkdEdx);
            }
          } else {
            registry.fill(HIST("nSigPiFromK0s"), posTrkPt, posTrack.tpcNSigmaPi());
            registry.fill(HIST("nSigPiFromK0s"), negTrkPt, negTrack.tpcNSigmaPi());
            registry.fill(HIST("NclVsEtaPiV0"), posTrkEta, posNcl);
            registry.fill(HIST("NclVsEtaPiV0p"), posTrkEta, posNcl);
            registry.fill(HIST("NclVsEtaPiV0"), negTrkEta, negNcl);
            registry.fill(HIST("NclVsEtaPiV0p"), negTrkEta, negNcl);
            nClVsPPiV0[posIndexEta]->Fill(posPorPTPC, posNcl);
            nClVsPpPiV0[posIndexEta]->Fill(posPorPTPC, posNcl);
            nClVsPPiV0[negIndexEta]->Fill(negPorPTPC, negNcl);
            nClVsPpPiV0[negIndexEta]->Fill(negPorPTPC, negNcl);
            dEdxPiV0[posIndexEta]->Fill(posPorPTPC, posTrkdEdx);
            dEdxPiV0[negIndexEta]->Fill(negPorPTPC, negTrkdEdx);
            if (posPorPTPC > KminPMIP && posPorPTPC < KmaxPMIP && posTrkdEdx > KmindEdxMIP && posTrkdEdx < KmaxdEdxMIP) {
              registry.fill(HIST("dEdxVsEtaPiMIPV0"), posTrkEta, posTrkdEdx);
              registry.fill(HIST("dEdxVsEtaPiMIPV0p"), posTrkEta, posTrkdEdx);
            }
            if (negPorPTPC > KminPMIP && negPorPTPC < KmaxPMIP && negTrkdEdx > KmindEdxMIP && negTrkdEdx < KmaxdEdxMIP) {
              registry.fill(HIST("dEdxVsEtaPiMIPV0"), negTrkEta, negTrkdEdx);
              registry.fill(HIST("dEdxVsEtaPiMIPV0p"), negTrkEta, negTrkdEdx);
            }
          }
        }
      }

      if (passesLambdaSelection(collision, v0)) {
        registry.fill(HIST("V0sCounter"), V0sCounter::Lambda);
        registry.fill(HIST("ArmL"), alpha, qT);
        registry.fill(HIST("MassLVsPt"), v0.pt(), v0.mLambda());
        // Distributions of positive-charge particles
        if (trackSelections.signCharge.value == "Positive" && posTrackCharge > KzeroInt) {
          registry.fill(HIST("nSigPrFromL"), posTrkPt, posTrack.tpcNSigmaPr());
          registry.fill(HIST("NclVsEtaPrV0"), posTrkEta, posNcl);
          registry.fill(HIST("NclVsEtaPrV0p"), posTrkEta, posNcl);
          nClVsPPrV0[posIndexEta]->Fill(posPorPTPC, posNcl);
          nClVsPpPrV0[posIndexEta]->Fill(posPorPTPC, posNcl);
          dEdxPrV0[posIndexEta]->Fill(posPorPTPC, posTrkdEdx);
          // Distributions of negative-charge particles
        } else if (trackSelections.signCharge.value == "Negative" && negTrackCharge < KzeroInt) {
          registry.fill(HIST("nSigPiFromL"), negTrkPt, negTrack.tpcNSigmaPi());
          registry.fill(HIST("NclVsEtaPiV0"), negTrkEta, negNcl);
          registry.fill(HIST("NclVsEtaPiV0p"), negTrkEta, negNcl);
          nClVsPpPiV0[negIndexEta]->Fill(negPorPTPC, negNcl);
          dEdxPiV0[negIndexEta]->Fill(negPorPTPC, negTrkdEdx);
          // Distributions of positive- and negative-charge particles
        } else {
          if (posTrackCharge > KzeroInt) {
            registry.fill(HIST("nSigPrFromL"), posTrkPt, posTrack.tpcNSigmaPr());
            registry.fill(HIST("NclVsEtaPrV0"), posTrkEta, posNcl);
            registry.fill(HIST("NclVsEtaPrV0p"), posTrkEta, posNcl);
            nClVsPPrV0[posIndexEta]->Fill(posPorPTPC, posNcl);
            nClVsPpPrV0[posIndexEta]->Fill(posPorPTPC, posNcl);
            dEdxPrV0[posIndexEta]->Fill(posPorPTPC, posTrkdEdx);
          }
          if (negTrackCharge < KzeroInt) {
            registry.fill(HIST("nSigPiFromL"), negTrkPt, negTrack.tpcNSigmaPi());
            registry.fill(HIST("NclVsEtaPiV0"), negTrkEta, negNcl);
            registry.fill(HIST("NclVsEtaPiV0p"), negTrkEta, negNcl);
            nClVsPpPiV0[negIndexEta]->Fill(negPorPTPC, negNcl);
            dEdxPiV0[negIndexEta]->Fill(negPorPTPC, negTrkdEdx);
          }
        }
      }

      if (passesAntiLambdaSelection(collision, v0)) {
        registry.fill(HIST("V0sCounter"), V0sCounter::AntiLambda);
        registry.fill(HIST("ArmAL"), alpha, qT);
        registry.fill(HIST("MassALVsPt"), v0.pt(), v0.mAntiLambda());
        // Distributions of positive-charge particles
        if (trackSelections.signCharge.value == "Positive" && posTrackCharge > KzeroInt) {
          registry.fill(HIST("nSigPiFromAL"), posTrkPt, posTrack.tpcNSigmaPi());
          registry.fill(HIST("NclVsEtaPiV0"), posTrkEta, posNcl);
          registry.fill(HIST("NclVsEtaPiV0p"), posTrkEta, posNcl);
          nClVsPpPiV0[posIndexEta]->Fill(posPorPTPC, posNcl);
          dEdxPiV0[posIndexEta]->Fill(posPorPTPC, posTrkdEdx);
          // Distributions of positive-charge particles
        } else if (trackSelections.signCharge.value == "Negative" && negTrackCharge < KzeroInt) {
          registry.fill(HIST("nSigPrFromAL"), negTrkPt, negTrack.tpcNSigmaPr());
          registry.fill(HIST("NclVsEtaPrV0"), negTrkEta, negNcl);
          registry.fill(HIST("NclVsEtaPrV0p"), negTrkEta, negNcl);
          nClVsPPrV0[negIndexEta]->Fill(negPorPTPC, negNcl);
          nClVsPpPrV0[negIndexEta]->Fill(negPorPTPC, negNcl);
          dEdxPrV0[negIndexEta]->Fill(negPorPTPC, negTrkdEdx);
          // Distributions of positive- and negative-charge particles
        } else {
          if (posTrackCharge > KzeroInt) {
            registry.fill(HIST("nSigPiFromAL"), posTrkPt, posTrack.tpcNSigmaPi());
            registry.fill(HIST("NclVsEtaPiV0"), posTrkEta, posNcl);
            registry.fill(HIST("NclVsEtaPiV0p"), posTrkEta, posNcl);
            nClVsPpPiV0[posIndexEta]->Fill(posPorPTPC, posNcl);
            dEdxPiV0[posIndexEta]->Fill(posPorPTPC, posTrkdEdx);
          }
          if (negTrackCharge < KzeroInt) {
            registry.fill(HIST("nSigPrFromAL"), negTrkPt, negTrack.tpcNSigmaPr());
            registry.fill(HIST("NclVsEtaPrV0"), negTrkEta, negNcl);
            registry.fill(HIST("NclVsEtaPrV0p"), negTrkEta, negNcl);
            nClVsPPrV0[negIndexEta]->Fill(negPorPTPC, negNcl);
            nClVsPpPrV0[negIndexEta]->Fill(negPorPTPC, negNcl);
            dEdxPrV0[negIndexEta]->Fill(negPorPTPC, negTrkdEdx);
          }
        }
      }

      if (passesGammaSelection(collision, v0)) {
        if (std::abs(alpha) < v0Selections.armAlphaSel && qT < v0Selections.qTSel) {

          if (trackSelections.applyPlateauSel && etaCal.isCalPlateauLoaded) {
            const double posDedxCal{etaCal.pEtaCalPlateau->GetBinContent(etaCal.pEtaCalPlateau->FindBin(posTrkEta))};
            const double negDedxCal{etaCal.pEtaCalPlateau->GetBinContent(etaCal.pEtaCalPlateau->FindBin(negTrkEta))};
            if (!(std::abs(posTrkdEdx - posDedxCal) < v0Selections.dEdxPlateauSel && std::abs(negTrkdEdx - negDedxCal) < v0Selections.dEdxPlateauSel)) {
              continue;
            }
          }

          registry.fill(HIST("V0sCounter"), V0sCounter::Gamma);
          registry.fill(HIST("ArmG"), alpha, qT);
          registry.fill(HIST("MassGVsPt"), v0.pt(), v0.mGamma());

          if (trackSelections.signCharge.value == "Positive" && posTrackCharge > KzeroInt) {
            registry.fill(HIST("nSigElFromG"), posTrkPt, posTrack.tpcNSigmaEl());
            registry.fill(HIST("NclVsEtaElV0"), posTrkEta, posNcl);
            registry.fill(HIST("NclVsEtaElV0p"), posTrkEta, posNcl);
            nClVsPElV0[posIndexEta]->Fill(posPorPTPC, posNcl);
            nClVsPpElV0[posIndexEta]->Fill(posPorPTPC, posNcl);
            if (posPorPTPC > KminPMIP && posPorPTPC < KmaxPMIP && posTrkdEdx > KmindEdxMIPPlateau && posTrkdEdx < KmaxdEdxMIPPlateau) {
              registry.fill(HIST("dEdxVsEtaElMIPV0"), posTrkEta, posTrkdEdx);
              registry.fill(HIST("dEdxVsEtaElMIPV0p"), posTrkEta, posTrkdEdx);
            }
            dEdxElV0[posIndexEta]->Fill(posPorPTPC, posTrkdEdx);
          } else if (trackSelections.signCharge.value == "Negative" && negTrackCharge < KzeroInt) {
            registry.fill(HIST("nSigElFromG"), negTrkPt, negTrack.tpcNSigmaEl());
            registry.fill(HIST("NclVsEtaElV0"), negTrkEta, negNcl);
            registry.fill(HIST("NclVsEtaElV0p"), negTrkEta, negNcl);
            nClVsPElV0[negIndexEta]->Fill(negPorPTPC, negNcl);
            nClVsPpElV0[negIndexEta]->Fill(negPorPTPC, negNcl);
            if (negPorPTPC > KminPMIP && negPorPTPC < KmaxPMIP && negTrkdEdx > KmindEdxMIPPlateau && negTrkdEdx < KmaxdEdxMIPPlateau) {
              registry.fill(HIST("dEdxVsEtaElMIPV0"), negTrkEta, negTrkdEdx);
              registry.fill(HIST("dEdxVsEtaElMIPV0p"), negTrkEta, negTrkdEdx);
            }
            dEdxElV0[negIndexEta]->Fill(negPorPTPC, negTrkdEdx);
          } else {
            registry.fill(HIST("nSigElFromG"), negTrkPt, negTrack.tpcNSigmaEl());
            registry.fill(HIST("nSigElFromG"), posTrkPt, posTrack.tpcNSigmaEl());
            registry.fill(HIST("NclVsEtaElV0"), posTrkEta, posNcl);
            registry.fill(HIST("NclVsEtaElV0p"), posTrkEta, posNcl);
            registry.fill(HIST("NclVsEtaElV0"), negTrkEta, negNcl);
            registry.fill(HIST("NclVsEtaElV0p"), negTrkEta, negNcl);
            nClVsPElV0[negIndexEta]->Fill(negPorPTPC, negNcl);
            nClVsPpElV0[negIndexEta]->Fill(negPorPTPC, negNcl);
            nClVsPElV0[posIndexEta]->Fill(posPorPTPC, posNcl);
            nClVsPpElV0[posIndexEta]->Fill(posPorPTPC, posNcl);
            if (posPorPTPC > KminPMIP && posPorPTPC < KmaxPMIP && posTrkdEdx > KmindEdxMIPPlateau && posTrkdEdx < KmaxdEdxMIPPlateau) {
              registry.fill(HIST("dEdxVsEtaElMIPV0"), posTrkEta, posTrkdEdx);
              registry.fill(HIST("dEdxVsEtaElMIPV0p"), posTrkEta, posTrkdEdx);
            }
            if (negPorPTPC > KminPMIP && negPorPTPC < KmaxPMIP && negTrkdEdx > KmindEdxMIPPlateau && negTrkdEdx < KmaxdEdxMIPPlateau) {
              registry.fill(HIST("dEdxVsEtaElMIPV0"), negTrkEta, negTrkdEdx);
              registry.fill(HIST("dEdxVsEtaElMIPV0p"), negTrkEta, negTrkdEdx);
            }
            dEdxElV0[posIndexEta]->Fill(posPorPTPC, posTrkdEdx);
            dEdxElV0[negIndexEta]->Fill(negPorPTPC, negTrkdEdx);
          }
        }
      }
    } // v0s
  }
  PROCESS_SWITCH(PiKpRAA, processCalibrationAndV0s, "Process QA", true);

  Preslice<TracksMC> perCollision = aod::track::collisionId;
  Service<o2::framework::O2DatabasePDG> pdg;

  void processSim(aod::McCollisions::iterator const& mccollision, soa::SmallGroups<ColEvSelsMC> const& collisions, BCsRun3 const& /*bcs*/, aod::FT0s const& /*ft0s*/, aod::McParticles const& mcParticles, TracksMC const& tracksMC)
  {

    const auto& kLowEta = (std::vector<double>)vecLowEta;
    const auto& kHighEta = (std::vector<double>)vecUpEta;

    //---------------------------
    // Only INEL > 0 generated collisions
    // By counting number of primary charged particles in |eta| < 1
    //---------------------------
    int nChMC{0};
    int nChMCTPCAcc{0};
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
      if (std::abs(charge) < KminCharge) {
        continue;
      }

      // Select particle based on its charge
      if (trackSelections.signCharge.value == "Positive" && charge < Kzero) {
        continue;
      }
      if (trackSelections.signCharge.value == "Negative" && charge > Kzero) {
        continue;
      }

      // Is it a primary particle?
      if (!particle.isPhysicalPrimary()) {
        continue;
      }

      const float eta{particle.eta()};

      // TVX requirement
      if (eta > KminFT0A && eta < KmaxFT0A) {
        nChFT0A++;
      }

      if (eta > KminFT0C && eta < KmaxFT0C) {
        nChFT0C++;
      }

      if (std::abs(eta) < tpcNchAcceptance) {
        nChMCTPCAcc++;
      }

      // INEL > 0
      if (std::abs(eta) > Kone) {
        continue;
      }

      nChMC++;
    }

    //---------------------------
    // Select only events with at least one charged particle in the FT0A and FT0C acceptances?
    //---------------------------
    if (selTVXMC) {
      if (!(nChFT0A > KzeroInt && nChFT0C > KzeroInt)) {
        registry.fill(HIST("NchMCcentVsTVX"), nChMC, 0.5);
        return;
      }
      registry.fill(HIST("NchMCcentVsTVX"), nChMC, 1.5);
    }

    //---------------------------
    // Select only MC events with |Vtx Z| < 10 cm?
    //---------------------------
    if (isZvtxPosSelMC && (std::fabs(mccollision.posZ()) > posZcut)) {
      return;
    }

    //---------------------------
    // Select only INEL > 0 generated events?
    //---------------------------
    if (selINELgt0) {
      if (!(nChMC > KzeroInt)) {
        return;
      }
    }

    //---------------------------
    // All Generated events irrespective of whether there is an associated reconstructed collision
    // Consequently, the centrality being a reconstructed quantity, might not always be available
    // Therefore it is expressed as a function of the generated pT and the generated Nch in ∣eta∣ < 0.8
    // This is used for the denominator of the signal loss correction
    // Also for MC closure: True Pt vs Generated Nch
    //---------------------------
    for (const auto& particle : mcParticles) {
      if (particle.eta() < trackSelections.minEta || particle.eta() > trackSelections.maxEta) {
        continue;
      }

      if (particle.pt() < trackSelections.minPt || particle.pt() > trackSelections.maxPt) {
        continue;
      }

      auto charge{0.};
      // Get the MC particle
      auto* pdgParticle = pdg->GetParticle(particle.pdgCode());
      if (pdgParticle != nullptr) {
        charge = pdgParticle->Charge();
      } else {
        continue;
      }

      // Is it a charged particle?
      if (std::abs(charge) < KminCharge) {
        continue;
      }

      // Select particle based on its charge
      if (trackSelections.signCharge.value == "Positive" && charge < Kzero) {
        continue;
      }
      if (trackSelections.signCharge.value == "Negative" && charge > Kzero) {
        continue;
      }

      // Is it a primary particle?
      bool isPrimary{true};
      if (!particle.isPhysicalPrimary()) {
        isPrimary = false;
      }

      if (isPrimary) {
        if (particle.pdgCode() == PDG_t::kPiPlus || particle.pdgCode() == PDG_t::kPiMinus) {
          registry.fill(HIST("PtPiVsNchMC_AllGen"), particle.pt(), nChMCTPCAcc);
          registry.fill(HIST("MCclosure_PtMCPiVsNchMC"), particle.pt(), nChMCTPCAcc);
        } else if (particle.pdgCode() == PDG_t::kKPlus || particle.pdgCode() == PDG_t::kKMinus) {
          registry.fill(HIST("PtKaVsNchMC_AllGen"), particle.pt(), nChMCTPCAcc);
          registry.fill(HIST("MCclosure_PtMCKaVsNchMC"), particle.pt(), nChMCTPCAcc);
        } else if (particle.pdgCode() == PDG_t::kProton || particle.pdgCode() == PDG_t::kProtonBar) {
          registry.fill(HIST("PtPrVsNchMC_AllGen"), particle.pt(), nChMCTPCAcc);
          registry.fill(HIST("MCclosure_PtMCPrVsNchMC"), particle.pt(), nChMCTPCAcc);
        } else {
          continue;
        }
      }
    } // Loop over Generated Particles

    //---------------------------
    //  This is used for the denominator of the event loss correction
    //  Only charge particles within tpcNchAcceptance and without pT selection
    //---------------------------
    registry.fill(HIST("NchMC_AllGen"), nChMCTPCAcc);

    //---------------------------
    // How many times the Generated evet was reconstrued?
    //---------------------------
    const auto& nRecColls{collisions.size()};
    registry.fill(HIST("NumberOfRecoCollisions"), nRecColls);

    //---------------------------
    // Only Generated evets with at least one reconstrued collision
    //---------------------------
    if (nRecColls > KzeroInt) {

      //---------------------------
      // Looks for the collision with the largest number of contributors
      // The selected collisions is identified with its index (bestCollisionIndex)
      //---------------------------
      int biggestNContribs{-1};
      int bestCollisionIndex{-1};
      for (const auto& collision : collisions) {

        float centrality{-999.0};
        if (centralitySelector.value == "FT0C") {
          centrality = collision.centFT0C();
        } else if (centralitySelector.value == "FT0M") {
          centrality = collision.centFT0M();
        } else if (centralitySelector.value == "FV0A") {
          centrality = collision.centFV0A();
        } else {
          centrality = -999.0;
        }

        if (selHasFT0 && !collision.has_foundFT0()) {
          continue;
        }

        if (biggestNContribs < collision.numContrib()) {
          biggestNContribs = collision.numContrib();
          bestCollisionIndex = collision.globalIndex();
        }

        if (selHasBC && !collision.has_foundBC()) {
          continue;
        }

        if (useSel8 && !collision.sel8()) {
          continue;
        }

        // kIsTriggerTVX
        if (selTriggerTVX && !collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
          continue;
        }

        // kNoITSROFrameBorder
        if (selNoITSROFrameBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
          continue;
        }

        // kNoTimeFrameBorder
        if (selNoTimeFrameBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
          continue;
        }

        // Zvtx
        if (isZvtxPosSel && std::fabs(collision.posZ()) > posZcut) {
          continue;
        }

        if (selIsGoodZvtxFT0vsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
          continue;
        }

        if (selNoSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
          continue;
        }

        //---------------------------
        // Needed to calculate denominator of the Event Splitting correction
        //---------------------------
        registry.fill(HIST("Centrality_AllRecoEvt"), centrality);
      }

      //---------------------------
      // Loop over the reconstructed collisions
      // Only that one with the largest number of contributors is considered
      //---------------------------
      for (const auto& collision : collisions) {

        float centrality{-999.0};
        if (centralitySelector.value == "FT0C") {
          centrality = collision.centFT0C();
        } else if (centralitySelector.value == "FT0M") {
          centrality = collision.centFT0M();
        } else if (centralitySelector.value == "FV0A") {
          centrality = collision.centFV0A();
        } else {
          centrality = -999.0;
        }

        //---------------------------
        // Pick the collisions with the largest number of contributors
        //---------------------------
        if (bestCollisionIndex != collision.globalIndex()) {
          continue;
        }

        // Needed to load the Phi selection from the CCDB
        const auto& foundBC = collision.foundBC_as<BCsRun3>();
        uint64_t timeStamp{foundBC.timestamp()};
        // const int magField{getMagneticField(timeStamp)};

        if (trackSelections.applyPhiCut) {
          const int nextRunNumber{foundBC.runNumber()};
          if (currentRunNumberPhiSel != nextRunNumber) {
            loadPhiCutSelections(timeStamp);
            currentRunNumberPhiSel = nextRunNumber;
            LOG(info) << "\tcurrentRunNumberPhiSel= " << currentRunNumberPhiSel << " timeStamp = " << timeStamp;
          }

          // return if phi cut objects are nullptr
          if (!(phiCut.hPhiCutHigh && phiCut.hPhiCutLow)) {
            return;
          }
        }

        //---------------------------
        // Needed to construct the correlation between MC Nch v.s. centrality
        //---------------------------
        registry.fill(HIST("Centrality_WRecoEvt"), centrality);
        registry.fill(HIST("zPosMC"), mccollision.posZ());

        registry.fill(HIST("CentralityVsBCVsFT0VsTVXVsEvSel"), centrality, 1.0);
        registry.fill(HIST("HasBCVsFT0VsTVXVsEvSel"), 1.0);

        if (collision.has_foundBC()) {
          registry.fill(HIST("CentralityVsBCVsFT0VsTVXVsEvSel"), centrality, 2.0);
          registry.fill(HIST("HasBCVsFT0VsTVXVsEvSel"), 2.0);
        }

        if (collision.has_foundBC() && collision.has_foundFT0()) {
          registry.fill(HIST("CentralityVsBCVsFT0VsTVXVsEvSel"), centrality, 3.0);
          registry.fill(HIST("HasBCVsFT0VsTVXVsEvSel"), 3.0);
        }

        if (collision.has_foundBC() && collision.has_foundFT0() && collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
          registry.fill(HIST("CentralityVsBCVsFT0VsTVXVsEvSel"), centrality, 4.0);
          registry.fill(HIST("HasBCVsFT0VsTVXVsEvSel"), 4.0);
        }

        //---------------------------
        // RCT Selection
        //---------------------------
        if (requireGoodRct) {
          // Checks if collisions passes RCT selection
          const bool isFT0Bad{requireBCRct ? foundBC.rct_bit(kFT0Bad) : collision.rct_bit(kFT0Bad)};
          const bool isITSBad{requireBCRct ? foundBC.rct_bit(kITSBad) : collision.rct_bit(kITSBad)};
          const bool isITSLimAcc{requireBCRct ? foundBC.rct_bit(kITSLimAccMCRepr) : collision.rct_bit(kITSLimAccMCRepr)};
          const bool isTOFBad{requireBCRct ? foundBC.rct_bit(kTOFBad) : collision.rct_bit(kTOFBad)};
          const bool isTOFLimAcc{requireBCRct ? foundBC.rct_bit(kTOFLimAccMCRepr) : collision.rct_bit(kTOFLimAccMCRepr)};
          const bool isTPCTrackingBad{requireBCRct ? foundBC.rct_bit(kTPCBadTracking) : collision.rct_bit(kTPCBadTracking)};
          const bool isTPCPIDBad{requireBCRct ? foundBC.rct_bit(kTPCBadPID) : collision.rct_bit(kTPCBadPID)};
          const bool isTPCLimAcc{requireBCRct ? foundBC.rct_bit(kTPCLimAccMCRepr) : collision.rct_bit(kTPCLimAccMCRepr)};

          registry.fill(HIST("CentralityVsRCTSel"), centrality, 1.0);
          if (!isFT0Bad) {
            registry.fill(HIST("CentralityVsRCTSel"), centrality, 2.0);
          }
          if (!isITSBad) {
            registry.fill(HIST("CentralityVsRCTSel"), centrality, 3.0);
          }
          if (!isITSLimAcc) {
            registry.fill(HIST("CentralityVsRCTSel"), centrality, 4.0);
          }
          if (!isTOFBad) {
            registry.fill(HIST("CentralityVsRCTSel"), centrality, 5.0);
          }
          if (!isTOFLimAcc) {
            registry.fill(HIST("CentralityVsRCTSel"), centrality, 6.0);
          }
          if (!isTPCTrackingBad) {
            registry.fill(HIST("CentralityVsRCTSel"), centrality, 7.0);
          }
          if (!isTPCPIDBad) {
            registry.fill(HIST("CentralityVsRCTSel"), centrality, 8.0);
          }
          if (!isTPCLimAcc) {
            registry.fill(HIST("CentralityVsRCTSel"), centrality, 9.0);
          }

          registry.fill(HIST("RCTSel"), 1.0);
          if (!rctChecker(collision)) {
            return;
          }

          registry.fill(HIST("RCTSel"), 2.0);
        }

        //---------------------------
        // Event Selection
        //---------------------------
        if (!isEventSelected(collision)) {
          return;
        }

        registry.fill(HIST("Centrality_WRecoEvtWSelCri"), centrality);
        registry.fill(HIST("NchMCVsCent"), centrality, nChMCTPCAcc);
        registry.fill(HIST("NchMC_WithRecoEvt"), nChMCTPCAcc); // Numerator of event loss correction
        registry.fill(HIST("zPos"), collision.posZ());
        registry.fill(HIST("Centrality"), centrality);

        //---------------------------
        // has_foundFT0() ?
        //---------------------------

        if (collision.has_foundBC() && collision.has_foundFT0() && collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
          registry.fill(HIST("CentralityVsBCVsFT0VsTVXVsEvSel"), centrality, 5.0);
          registry.fill(HIST("HasBCVsFT0VsTVXVsEvSel"), 5.0);
        }

        if (collision.has_foundFT0()) {
          registry.fill(HIST("CentralityVsFoundFT0"), centrality, 1.5);
        } else {
          registry.fill(HIST("CentralityVsFoundFT0"), centrality, 0.5);
        }

        //---------------------------
        // All Generated events with at least one associated reconstructed collision
        // The Generated events are not subjected to any selection criteria
        // However, the associated reconstructed collisions pass the selection criteria
        // This histograms are used for the denominator of the tracking efficiency
        //---------------------------
        for (const auto& particle : mcParticles) {
          if (particle.eta() < trackSelections.minEta || particle.eta() > trackSelections.maxEta) {
            continue;
          }

          if (particle.pt() < trackSelections.minPt || particle.pt() > trackSelections.maxPt) {
            continue;
          }

          auto charge{0.};
          // Get the MC particle
          auto* pdgParticle = pdg->GetParticle(particle.pdgCode());
          if (pdgParticle != nullptr) {
            charge = pdgParticle->Charge();
          } else {
            continue;
          }

          // Is it a charged particle?
          if (std::abs(charge) < KminCharge) {
            continue;
          }

          // Select particle based on its charge
          if (trackSelections.signCharge.value == "Positive" && charge < Kzero) {
            continue;
          }
          if (trackSelections.signCharge.value == "Negative" && charge > Kzero) {
            continue;
          }

          // Is it a primary particle?
          bool isPrimary{true};
          if (!particle.isPhysicalPrimary()) {
            isPrimary = false;
          }

          if (isPrimary) {
            if (particle.pdgCode() == PDG_t::kPiPlus || particle.pdgCode() == PDG_t::kPiMinus) {
              registry.fill(HIST("PtPiVsCentMC_WithRecoEvt"), particle.pt(), centrality); // Denominator of tracking efficiency
              registry.fill(HIST("PtPiVsNchMC_WithRecoEvt"), particle.pt(), nChMCTPCAcc); // Numerator of signal loss
            } else if (particle.pdgCode() == PDG_t::kKPlus || particle.pdgCode() == PDG_t::kKMinus) {
              registry.fill(HIST("PtKaVsCentMC_WithRecoEvt"), particle.pt(), centrality);
              registry.fill(HIST("PtKaVsNchMC_WithRecoEvt"), particle.pt(), nChMCTPCAcc);
            } else if (particle.pdgCode() == PDG_t::kProton || particle.pdgCode() == PDG_t::kProtonBar) {
              registry.fill(HIST("PtPrVsCentMC_WithRecoEvt"), particle.pt(), centrality);
              registry.fill(HIST("PtPrVsNchMC_WithRecoEvt"), particle.pt(), nChMCTPCAcc);
            } else {
              continue;
            }
          }
        } // Loop over generated particles per generated collision

        const auto& groupedTracks{tracksMC.sliceBy(perCollision, collision.globalIndex())};

        // ================
        // Track selection WITHOUT DCAxy and DCAz selections
        // Needed for the Secondary Particle Correction
        // ================
        for (const auto& track : groupedTracks) {

          const bool applyDca{false};
          if (!selectPrimary(track, applyDca)) {
            continue;
          }

          if (!track.has_mcParticle()) {
            continue;
          }

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
          if (std::abs(charge) < KminCharge) {
            continue;
          }

          // Select particle based on its charge
          if (trackSelections.signCharge.value == "Positive" && charge < Kzero) {
            continue;
          }
          if (trackSelections.signCharge.value == "Negative" && charge > Kzero) {
            continue;
          }

          registry.fill(HIST("DCAxyVsPt"), track.pt(), track.dcaXY());
          registry.fill(HIST("DCAzVsPt"), track.pt(), track.dcaZ());

          // float phiPrime{track.phi()};
          // phiPrimeFunc(phiPrime, magField, charge);

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
            if (isPi && !isPr) {
              registry.fill(HIST("dcaVsPtPi"), track.pt(), track.dcaXY(), centrality);
            }
            if (isPr && !isPi) {
              registry.fill(HIST("dcaVsPtPr"), track.pt(), track.dcaXY(), centrality);
            }
          }

          if (isDecay && !isPrimary && !isMaterial) {
            if (isPi && !isPr) {
              registry.fill(HIST("dcaVsPtPiDec"), track.pt(), track.dcaXY(), centrality);
            }
            if (isPr && !isPi) {
              registry.fill(HIST("dcaVsPtPrDec"), track.pt(), track.dcaXY(), centrality);
            }
          }

          if (isMaterial && !isPrimary && !isDecay) {
            if (isPi && !isPr) {
              registry.fill(HIST("dcaVsPtPiMat"), track.pt(), track.dcaXY(), centrality);
            }
            if (isPr && !isPi) {
              registry.fill(HIST("dcaVsPtPrMat"), track.pt(), track.dcaXY(), centrality);
            }
          }
        }

        // ================
        // Track selection WITH DCAxy and DCAz selections
        // Needed for the numerator in the tracking efficiency
        // ================
        int nCh{0};
        for (const auto& track : groupedTracks) {

          const bool applyDca{true};
          if (!selectPrimary(track, applyDca)) {
            continue;
          }

          // Has MC particle?
          if (!track.has_mcParticle()) {
            continue;
          }

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
          if (std::abs(charge) < KminCharge) {
            continue;
          }

          // Select particle based on its charge
          if (trackSelections.signCharge.value == "Positive" && charge < Kzero) {
            continue;
          }
          if (trackSelections.signCharge.value == "Negative" && charge > Kzero) {
            continue;
          }

          // float phiPrime{track.phi()};
          // phiPrimeFunc(phiPrime, magField, charge);

          int indexEta{-999};
          const float eta{track.eta()};
          for (int i = 0; i < KnEtaHists; ++i) {
            if (eta >= kLowEta[i] && eta < kHighEta[i]) {
              indexEta = i;
              break;
            }
          }

          if (indexEta < KzeroInt || indexEta > KsevenInt) {
            continue;
          }

          nCh++;

          registry.fill(HIST("MomentumTPCVsP"), track.p(), track.tpcInnerParam());
          registry.fill(HIST("DCAxyVsPtWithSelection"), track.pt(), track.dcaXY());
          registry.fill(HIST("DCAzVsPtWithSelection"), track.pt(), track.dcaZ());
          registry.fill(HIST("EtaVsPhi"), track.eta(), track.phi());
          registry.fill(HIST("NclVsPhiHist"), track.phi(), track.tpcNClsFound());
          registry.fill(HIST("NclVsPhiProfile"), track.phi(), track.tpcNClsFound());
          registry.fill(HIST("NclVsEta"), track.eta(), track.tpcNClsFound());
          registry.fill(HIST("NclVsEtap"), track.eta(), track.tpcNClsFound());
          registry.fill(HIST("NclVsEtaPID"), track.eta(), track.tpcNClsPID());
          registry.fill(HIST("NclVsEtaPIDp"), track.eta(), track.tpcNClsPID());

          bool isPrimary{false};
          if (particle.isPhysicalPrimary()) {
            isPrimary = true;
          }

          if (!isPrimary) {
            continue;
          }

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
            registry.fill(HIST("PtPiVsCent_WithRecoEvt"), track.pt(), centrality);       // Numerator of tracking efficiency
            registry.fill(HIST("PtGenPiVsCent_WithRecoEvt"), particle.pt(), centrality); // Numerator of tracking efficiency
            registry.fill(HIST("MCclosure_PtPiVsNchMC"), track.pt(), nChMCTPCAcc);
          }
          if (isKa && !isPi && !isPr) {
            registry.fill(HIST("PtKaVsCent_WithRecoEvt"), track.pt(), centrality);
            registry.fill(HIST("PtGenKaVsCent_WithRecoEvt"), particle.pt(), centrality);
            registry.fill(HIST("MCclosure_PtKaVsNchMC"), track.pt(), nChMCTPCAcc);
          }
          if (isPr && !isPi && !isKa) {
            registry.fill(HIST("PtPrVsCent_WithRecoEvt"), track.pt(), centrality);
            registry.fill(HIST("PtGenPrVsCent_WithRecoEvt"), particle.pt(), centrality);
            registry.fill(HIST("MCclosure_PtPrVsNchMC"), track.pt(), nChMCTPCAcc);
          }
          registry.fill(HIST("PtResolution"), particle.pt(), (track.pt() - particle.pt()) / particle.pt());
        } // Loop over reconstructed tracks
        registry.fill(HIST("NchVsCent"), centrality, nCh);
      } // Loop over Reco. Collisions: Only the collisions with the largest number of contributors
    } // If condition: Only simulated evets with at least one reconstrued collision
  }
  PROCESS_SWITCH(PiKpRAA, processSim, "Process Sim", false);

  template <typename T, typename U>
  void getArmeterosVariables(const T& ppos, const T& pneg, U& alpha, U& qT)
  {
    alpha = 0., qT = 0.;
    TVector3 pV0 = ppos + pneg;
    double pV0mag = pV0.Mag();
    if (pV0mag < KtEnToMinusNine)
      return; // protect against zero momentum

    const TVector3 u = pV0 * (1.0 / pV0mag);

    double pLpos = ppos.Dot(u);
    double pLneg = pneg.Dot(u);

    // qT: transverse momentum of the + track w.r.t. V0 direction
    TVector3 pTpos = ppos - pLpos * u;
    qT = pTpos.Mag();

    // α: longitudinal asymmetry (uses + and − labels by charge)
    double denom = pLpos + pLneg;
    if (std::abs(denom) < KtEnToMinusNine) {
      return;
    } // avoid 0 division (unphysical for V0s)

    alpha = (pLpos - pLneg) / denom; // equivalently / pV0mag
  }

  template <typename T>
  bool selectPrimary(T const& track, const bool& applyDca)
  {
    const float dcaXY{track.dcaXY()};
    const float dcaZ{track.dcaZ()};
    const float pt{track.pt()};
    const int16_t nclFound{track.tpcNClsFound()};
    const int16_t nclPID{track.tpcNClsPID()};
    const int16_t ncl{trackSelections.useNclsPID ? nclPID : nclFound};
    double dcaXYcut{0.0};
    double dcaZcut{0.0};

    if (cfgDCA.dcaSelectionsLoaded) {
      dcaXYcut = cfgDCA.hDCAxy->GetBinContent(1) + cfgDCA.hDCAxy->GetBinContent(2) / std::pow(std::abs(pt), cfgDCA.hDCAxy->GetBinContent(3));
      dcaZcut = cfgDCA.hDCAz->GetBinContent(1) + cfgDCA.hDCAz->GetBinContent(2) / std::pow(std::abs(pt), cfgDCA.hDCAz->GetBinContent(3));
      dcaXYcut *= trackSelections.nSigmaDCAxy;
      dcaZcut *= trackSelections.nSigmaDCAz;
    }

    if (track.eta() < trackSelections.minEta || track.eta() > trackSelections.maxEta) {
      return false;
    }

    if (pt < trackSelections.minPt || pt > trackSelections.maxPt) {
      return false;
    }

    // If the track produces a hit in neither the first nor the second ITS layer, it is excluded
    // The track must produce a hit in either the 1st or the 2nd ITS layer to be accepted
    if (!(track.itsClusterMap() & 0x01) && !(track.itsClusterMap() & 0x02)) {
      return false;
    }

    if (track.itsNCls() < trackSelections.minNClusITS ||
        track.tpcNClsCrossedRows() < trackSelections.minNCrossedRows ||
        track.tpcChi2NCl() > trackSelections.maxChi2ClsTPC ||
        track.tpcChi2NCl() < trackSelections.minChi2ClsTPC ||
        track.itsChi2NCl() > trackSelections.chi2ClsITS) {
      return false;
    }

    // ==== Ncl selection ==== //
    if (trackSelections.applyNclSel && ncl < trackSelections.minNcl) {
      return false;
    }

    // ==== DCAxy & DCAz selections ==== //
    if (applyDca && cfgDCA.dcaSelectionsLoaded) {
      if (std::abs(dcaZ) > dcaZcut || std::abs(dcaXY) > dcaXYcut) {
        return false;
      }
    }

    // Flag to check that the DCA selections are loaded
    // when asking to apply DCA
    // ==== DCAxy & DCAz selections ==== //
    if (loadHisWithDCASel && applyDca && !cfgDCA.dcaSelectionsLoaded) {
      return false;
    }

    return true;
  }

  template <typename T>
  bool selectV0Daughter(const T& track)
  {

    // const float eta{track.eta()};
    const uint8_t nclShared{track.tpcNClsShared()};
    const int16_t nCrossedRows{track.tpcNClsCrossedRows()};
    const int16_t nclFound{track.tpcNClsFound()};
    const int16_t nclPID{track.tpcNClsPID()};
    const int16_t ncl{trackSelections.useNclsPID ? nclPID : nclFound};

    if (track.eta() < trackSelections.minEta || track.eta() > trackSelections.maxEta) {
      return false;
    }

    if (track.pt() < trackSelections.minPt || track.pt() > trackSelections.maxPt) {
      return false;
    }

    // ==== Crossed rows TPC ==== //
    if (nCrossedRows < trackSelections.minNCrossedRows) {
      return false;
    }

    // ==== Ncl selection ==== //
    if (trackSelections.applyNclSel && ncl < v0Selections.minNclV0Daugther) {
      return false;
    }

    // ==== Ncl shared ==== //
    if (nclShared > v0Selections.nSharedClusTpc) {
      return false;
    }

    return true;
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

    const bool hasToF{posTrack.hasTOF() && negTrack.hasTOF() ? true : false};
    const bool goodToFmatch{posTrack.goodTOFMatch() && negTrack.goodTOFMatch() ? true : false};
    const double dcaPos{std::fabs(v0.dcapostopv())};
    const double dcaNeg{std::fabs(v0.dcanegtopv())};
    const float posTPCNsigma{std::fabs(posTrack.tpcNSigmaPi())};
    const float negTPCNsigma{std::fabs(negTrack.tpcNSigmaPi())};
    const float posTOFNsigma{std::fabs(posTrack.tofNSigmaPi())};
    const float negTOFNsigma{std::fabs(negTrack.tofNSigmaPi())};
    const double dMassK0s{std::abs(v0.mK0Short() - o2::constants::physics::MassK0Short)};
    const double dMassL{std::abs(v0.mLambda() - o2::constants::physics::MassLambda0)};
    const double dMassAL{std::abs(v0.mAntiLambda() - o2::constants::physics::MassLambda0)};
    const double dMassG{std::abs(v0.mGamma() - o2::constants::physics::MassGamma)};
    const double rapidity{v0.yK0Short()};
    const double lifeTime{v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short};

    // Checks if DCA of daughters to PV passes selection
    const bool dcaDaugToPV{dcaPos > v0Selections.dcapostopv && dcaNeg > v0Selections.dcanegtopv ? true : false};

    // Rejects V0 if its invariant mass is not compatible with the K0s proper mass
    if (v0Selections.applyInvMassSel) {
      if (!(dMassK0s < v0Selections.dMassSel && dMassL > v0Selections.dMassSel && dMassAL > v0Selections.dMassSel && dMassG > v0Selections.dMassSelG)) {
        return false;
      }
    }

    bool isSelected{false};
    if (v0Selections.useOfficialV0sSelOfDaughters) {
      isSelected = lifeTime < v0Selections.lifeTimeCutK0s && (rapidity > v0Selections.minY && rapidity < v0Selections.maxY) && dcaDaugToPV ? true : false;
    }
    if (!v0Selections.useOfficialV0sSelOfDaughters) {
      if (v0Selections.useTPCNsigma) {
        isSelected = lifeTime < v0Selections.lifeTimeCutK0s && (rapidity > v0Selections.minY && rapidity < v0Selections.maxY) && posTPCNsigma < v0Selections.pidNsigmaCut && negTPCNsigma < v0Selections.pidNsigmaCut && dcaDaugToPV ? true : false;
      }
      if (!v0Selections.useTPCNsigma) {
        isSelected = lifeTime < v0Selections.lifeTimeCutK0s && (rapidity > v0Selections.minY && rapidity < v0Selections.maxY) && posTOFNsigma < v0Selections.pidNsigmaCut && negTOFNsigma < v0Selections.pidNsigmaCut && hasToF && goodToFmatch && dcaDaugToPV ? true : false;
      }
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

    const bool hasToF{posTrack.hasTOF() && negTrack.hasTOF() ? true : false};
    const bool goodToFmatch{posTrack.goodTOFMatch() && negTrack.goodTOFMatch() ? true : false};
    const double dcaPos{std::fabs(v0.dcapostopv())};
    const double dcaNeg{std::fabs(v0.dcanegtopv())};
    const float posTPCNsigma{std::fabs(posTrack.tpcNSigmaPr())};
    const float negTPCNsigma{std::fabs(negTrack.tpcNSigmaPi())};
    const float posTOFNsigma{std::fabs(posTrack.tofNSigmaPr())};
    const float negTOFNsigma{std::fabs(negTrack.tofNSigmaPi())};
    const double dMassK0s{std::abs(v0.mK0Short() - o2::constants::physics::MassK0Short)};
    const double dMassL{std::abs(v0.mLambda() - o2::constants::physics::MassLambda0)};
    const double dMassG{std::abs(v0.mGamma() - o2::constants::physics::MassGamma)};
    const double rapidity{v0.yLambda()};
    const double lifeTime{v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0};

    // Checks if DCA of daughters to PV passes selection
    const bool dcaDaugToPV{dcaPos > v0Selections.dcaProtonFromLambda && dcaNeg > v0Selections.dcaPionFromLambda ? true : false};

    // Rejects V0 if the invariant mass is not compatible with the Lambda proper mass
    if (v0Selections.applyInvMassSel) {
      if (!(dMassL < v0Selections.dMassSel && dMassK0s > v0Selections.dMassSel && dMassG > v0Selections.dMassSelG)) {
        return false;
      }
    }

    bool isSelected{false};
    if (v0Selections.useOfficialV0sSelOfDaughters) {
      isSelected = lifeTime < v0Selections.lifeTimeCutLambda && (rapidity > v0Selections.minY && rapidity < v0Selections.maxY) && dcaDaugToPV ? true : false;
    }
    if (!v0Selections.useOfficialV0sSelOfDaughters) {
      if (v0Selections.useTPCNsigma) {
        isSelected = lifeTime < v0Selections.lifeTimeCutLambda && (rapidity > v0Selections.minY && rapidity < v0Selections.maxY) && posTPCNsigma < v0Selections.pidNsigmaCut && negTPCNsigma < v0Selections.pidNsigmaCut && dcaDaugToPV ? true : false;
      }
      if (!v0Selections.useTPCNsigma) {
        isSelected = lifeTime < v0Selections.lifeTimeCutLambda && (rapidity > v0Selections.minY && rapidity < v0Selections.maxY) && posTOFNsigma < v0Selections.pidNsigmaCut && negTOFNsigma < v0Selections.pidNsigmaCut && hasToF && goodToFmatch && dcaDaugToPV ? true : false;
      }
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

    const bool hasToF{posTrack.hasTOF() && negTrack.hasTOF() ? true : false};
    const bool goodToFmatch{posTrack.goodTOFMatch() && negTrack.goodTOFMatch() ? true : false};
    const double dcaPos{std::fabs(v0.dcapostopv())};
    const double dcaNeg{std::fabs(v0.dcanegtopv())};
    const float posTPCNsigma{std::fabs(posTrack.tpcNSigmaPi())};
    const float negTPCNsigma{std::fabs(negTrack.tpcNSigmaPr())};
    const float posTOFNsigma{std::fabs(posTrack.tofNSigmaPi())};
    const float negTOFNsigma{std::fabs(negTrack.tofNSigmaPr())};
    const double dMassK0s{std::abs(v0.mK0Short() - o2::constants::physics::MassK0Short)};
    const double dMassAL{std::abs(v0.mAntiLambda() - o2::constants::physics::MassLambda0)};
    const double dMassG{std::abs(v0.mGamma() - o2::constants::physics::MassGamma)};
    const double rapidity{v0.yLambda()};
    const double lifeTime{v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0};

    // Checks if DCA of daughters to PV passes selection
    const bool dcaDaugToPV{dcaPos > v0Selections.dcaPionFromLambda && dcaNeg > v0Selections.dcaProtonFromLambda ? true : false};

    // Rejects V0 if the invariant mass is not compatible with the Lambda proper mass
    if (v0Selections.applyInvMassSel) {
      if (!(dMassAL < v0Selections.dMassSel && dMassK0s > v0Selections.dMassSel && dMassG > v0Selections.dMassSelG)) {
        return false;
      }
    }

    bool isSelected{false};
    if (v0Selections.useOfficialV0sSelOfDaughters) {
      isSelected = lifeTime < v0Selections.lifeTimeCutLambda && (rapidity > v0Selections.minY && rapidity < v0Selections.maxY) && dcaDaugToPV ? true : false;
    }
    if (!v0Selections.useOfficialV0sSelOfDaughters) {
      if (v0Selections.useTPCNsigma) {
        isSelected = lifeTime < v0Selections.lifeTimeCutLambda && (rapidity > v0Selections.minY && rapidity < v0Selections.maxY) && posTPCNsigma < v0Selections.pidNsigmaCut && negTPCNsigma < v0Selections.pidNsigmaCut && dcaDaugToPV ? true : false;
      }
      if (!v0Selections.useTPCNsigma) {
        isSelected = lifeTime < v0Selections.lifeTimeCutLambda && (rapidity > v0Selections.minY && rapidity < v0Selections.maxY) && posTOFNsigma < v0Selections.pidNsigmaCut && negTOFNsigma < v0Selections.pidNsigmaCut && hasToF && goodToFmatch && dcaDaugToPV ? true : false;
      }
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
    const float rapidity = RecoDecay::y(std::array{v0.px(), v0.py(), v0.pz()}, o2::constants::physics::MassGamma);

    // Checks if DCA of daughters to PV passes selection
    const bool dcaDaugToPV{dcaPos > v0Selections.dcaElectronFromGamma && dcaNeg > v0Selections.dcaElectronFromGamma ? true : false};

    if (v0Selections.applyInvMassSel) {
      if (!(dMassK0s > v0Selections.dMassSel && dMassL > v0Selections.dMassSel && dMassAL > v0Selections.dMassSel && dMassG < v0Selections.dMassSel)) {
        return false;
      }
    }

    if (!(rapidity > v0Selections.minY && rapidity < v0Selections.maxY)) {
      return false;
    }

    bool isSelected{false};
    if (v0Selections.useOfficialV0sSelOfDaughters) {
      isSelected = dcaDaugToPV ? true : false;
    }
    if (!v0Selections.useOfficialV0sSelOfDaughters) {
      isSelected = dcaDaugToPV && posTPCNsigma < v0Selections.pidNsigmaCut && negTPCNsigma < v0Selections.pidNsigmaCut ? true : false;
    }

    if (isSelected) {
      registry.fill(HIST("EtaVsYG"), negTrack.eta(), rapidity);
      registry.fill(HIST("EtaVsYG"), posTrack.eta(), rapidity);
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
    if (magField < 0) {
      phi = o2::constants::math::TwoPI - phi;
    }
    if (charge < 0) {
      phi = o2::constants::math::TwoPI - phi;
    }

    phi += o2::constants::math::PI / 18.0f;
    phi = std::fmod(phi, o2::constants::math::PI / 9.0f);
  }

  bool passesPhiSelection(const float& pt, const float& phi)
  {
    // Do not apply Phi Sel if pt < 2 GeV/c
    if (pt < KtwoPtGeVSel) {
      return true;
    }

    bool isSelected{true};
    if (phiCut.isPhiCutLoaded) {
      const int binLow{phiCut.hPhiCutLow->FindBin(pt)};
      const int binHigh{phiCut.hPhiCutHigh->FindBin(pt)};
      const double phiCutLow{phiCut.hPhiCutLow->GetBinContent(binLow)};
      const double phiCutHigh{phiCut.hPhiCutHigh->GetBinContent(binHigh)};
      if (phi >= phiCutLow && phi <= phiCutHigh) {
        isSelected = false;
      }
    }
    return isSelected;
  }

  template <typename CheckCol>
  bool isEventSelected(CheckCol const& col)
  {
    registry.fill(HIST("EventCounter"), EvCutLabel::All);

    // Has BC?
    if (selHasBC) {
      if (!col.has_foundBC()) {
        return false;
      }
      registry.fill(HIST("EventCounter"), EvCutLabel::HasBC);
    }

    // Has FT0 information?
    if (selHasFT0) {
      if (!col.has_foundFT0()) {
        return false;
      }
      registry.fill(HIST("EventCounter"), EvCutLabel::HasFT0);
    }

    if (useSel8) {
      if (!col.sel8()) {
        return false;
      }
      registry.fill(HIST("EventCounter"), EvCutLabel::SelEigth);
    }

    // kIsTriggerTVX
    if (selTriggerTVX) {
      if (!col.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
        return false;
      }
      registry.fill(HIST("EventCounter"), EvCutLabel::SelTriggerTVX);
    }

    // kNoITSROFrameBorder
    if (selNoITSROFrameBorder) {
      if (!col.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
        return false;
      }
      registry.fill(HIST("EventCounter"), EvCutLabel::SelNoITSROFrameBorder);
    }

    // kNoTimeFrameBorder
    if (selNoTimeFrameBorder) {
      if (!col.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
        return false;
      }
      registry.fill(HIST("EventCounter"), EvCutLabel::SelNoTimeFrameBorder);
    }

    // Zvtx
    if (isZvtxPosSel) {
      if (std::fabs(col.posZ()) > posZcut) {
        return false;
      }
      registry.fill(HIST("EventCounter"), EvCutLabel::VtxZ);
    }

    if (selIsGoodZvtxFT0vsPV) {
      if (!col.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
        return false;
      }
      registry.fill(HIST("EventCounter"), EvCutLabel::IsGoodZvtxFT0vsPV);
    }

    if (selNoSameBunchPileup) {
      if (!col.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
        return false;
      }
      registry.fill(HIST("EventCounter"), EvCutLabel::NoSameBunchPileup);
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

    if (isCentSel) {
      if (col.centFT0C() < minCentCut || col.centFT0C() > maxCentCut) {
        return false;
      }
      registry.fill(HIST("EventCounter"), EvCutLabel::Centrality);
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

  void loadDCAselections()
  {
    if (pathDCAxy.value.empty() == false) {
      cfgDCA.hDCAxy = ccdb->getForTimeStamp<TH1F>(pathDCAxy, ccdbNoLaterThan.value);
      if (cfgDCA.hDCAxy == nullptr) {
        LOGF(fatal, "Could not load hDCAxy histogram from %s", pathDCAxy.value.c_str());
      }
    }

    if (pathDCAz.value.empty() == false) {
      cfgDCA.hDCAz = ccdb->getForTimeStamp<TH1F>(pathDCAz, ccdbNoLaterThan.value);
      if (cfgDCA.hDCAz == nullptr) {
        LOGF(fatal, "Could not load hDCAz histogram from %s", pathDCAz.value.c_str());
      }
    }
    if (cfgDCA.hDCAxy && cfgDCA.hDCAz) {
      cfgDCA.dcaSelectionsLoaded = true;
    }
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

    if (phiCut.hPhiCutHigh && phiCut.hPhiCutLow) {
      phiCut.isPhiCutLoaded = true;
    }
  }

  void loadEtaCalibration()
  {
    if (pathEtaCal.value.empty() == false) {
      etaCal.pEtaCal = ccdb->getForTimeStamp<TProfile>(pathEtaCal, ccdbNoLaterThan.value);
      if (etaCal.pEtaCal == nullptr) {
        LOGF(fatal, "Could not load pEtaCal from %s", pathEtaCal.value.c_str());
      }
    }

    if (etaCal.pEtaCal) {
      etaCal.isMIPCalLoaded = true;
    }
  }

  void loadEtaPlateauCalibration()
  {
    if (pathEtaCalPlateau.value.empty() == false) {
      etaCal.pEtaCalPlateau = ccdb->getForTimeStamp<TProfile>(pathEtaCalPlateau, ccdbNoLaterThan.value);

      if (etaCal.pEtaCalPlateau == nullptr) {
        LOGF(fatal, "Could not load pEtaCalPlateau from %s", pathEtaCalPlateau.value.c_str());
      }
    }

    if (etaCal.pEtaCalPlateau) {
      etaCal.isCalPlateauLoaded = true;
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PiKpRAA>(cfgc)};
}
