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

#include "TPDGCode.h"
#include "TVector3.h"
#include <TString.h>

#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
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

using ColEvSels = soa::Join<aod::Collisions, aod::EvSels, aod::FT0MultZeqs, o2::aod::CentFT0Cs, aod::TPCMults, o2::aod::BarrelMults>;
using BCsRun3 = soa::Join<aod::BCsWithTimestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;

using TracksFull = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelectionExtension, aod::TracksDCA, aod::TrackSelection, aod::TracksCovIU, aod::pidTPCPi, aod::pidTPCPr, aod::pidTOFPr, aod::pidTPCEl, aod::pidTOFFlags, aod::pidTOFbeta, aod::TOFSignal, aod::pidTOFFullPi, aod::pidTOFFullEl, aod::pidTOFFullMu>;

static constexpr int kNEtaHists{8};

std::array<std::shared_ptr<TH3>, kNEtaHists> dEdxPiV0{};
std::array<std::shared_ptr<TH3>, kNEtaHists> dEdxPrV0{};
std::array<std::shared_ptr<TH3>, kNEtaHists> dEdxElV0{};
std::array<std::shared_ptr<TH3>, kNEtaHists> dEdxPiTOF{};
std::array<std::shared_ptr<TH3>, kNEtaHists> dEdx{};
std::array<std::shared_ptr<TH2>, kNEtaHists> pTVsP{};
std::array<std::shared_ptr<TH2>, kNEtaHists> nClVsP{};
std::array<std::shared_ptr<TH2>, kNEtaHists> nClVsPElV0{};
std::array<std::shared_ptr<TH2>, kNEtaHists> nClVsPPiV0{};
std::array<std::shared_ptr<TH2>, kNEtaHists> nClVsPPrV0{};
std::array<std::shared_ptr<TProfile>, kNEtaHists> nClVsPp{};
std::array<std::shared_ptr<TProfile>, kNEtaHists> nClVsPpElV0{};
std::array<std::shared_ptr<TProfile>, kNEtaHists> nClVsPpPiV0{};
std::array<std::shared_ptr<TProfile>, kNEtaHists> nClVsPpPrV0{};

struct PiKpRAA {

  static constexpr int kZeroInt{0};
  static constexpr int kSevenInt{7};

  static constexpr float kZero{0.0f};
  static constexpr float kOne{1.0f};
  static constexpr float kTenToMinusNine{1e-9};
  static constexpr float kMinPtNchSel{0.1f};
  static constexpr float kMaxPtNchSel{3.0f};
  static constexpr float kMinCharge{3.f};
  static constexpr float kMinPElMIP{0.3f};
  static constexpr float kMaxPElMIP{0.4f};
  static constexpr float kMinPMIP{0.4f};
  static constexpr float kMaxPMIP{0.6f};
  static constexpr float kMindEdxMIP{40.0f};
  static constexpr float kMaxdEdxMIP{60.0f};
  static constexpr float kMindEdxMIPPlateau{70.0f};
  static constexpr float kMaxdEdxMIPPlateau{90.0f};

  static constexpr float kLowEta[kNEtaHists] = {-0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6};
  static constexpr float kHighEta[kNEtaHists] = {-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8};
  // static constexpr float kLowEta[kNEtaHists] = {0.0, 0.2, 0.4, 0.6};
  // static constexpr float kHighEta[kNEtaHists] = {0.2, 0.4, 0.6, 0.8};

  static constexpr float DefaultLifetimeCuts[1][2] = {{30., 20.}};
  Configurable<LabeledArray<float>> lifetimecut{"lifetimecut", {DefaultLifetimeCuts[0], 2, {"lifetimecutLambda", "lifetimecutK0S"}}, "lifetimecut"};

  struct : ConfigurableGroup {
    Configurable<uint8_t> v0TypeSelection{"v0TypeSelection", 1, "select on a certain V0 type (leave negative if no selection desired)"};

    // Selection criteria: acceptance
    Configurable<float> rapidityCut{"rapidityCut", 0.5, "rapidity"};
    Configurable<float> minEtaDaughter{"minEtaDaughter", -0.8, "Daughter minimum-eta selection"};
    Configurable<float> maxEtaDaughter{"maxEtaDaughter", +0.8, "Daughter maximum-eta selection"};
    Configurable<float> minPt{"minPt", 0.15, "minimum pt of the tracks"};
    Configurable<float> maxPt{"maxPt", 20.0, "maximum pt of the tracks"};
    Configurable<int16_t> minNclFound{"minNclFound", 135, "minimum found Ncl in TPC"};

    // Standard 5 topological criteria
    Configurable<float> v0cospa{"v0cospa", 0.995, "min V0 CosPA"};
    Configurable<float> dcav0dau{"dcav0dau", 1.0, "max DCA V0 Daughters (cm)"};
    Configurable<float> dcanegtopv{"dcanegtopv", .1, "min DCA Neg To PV (cm)"};
    Configurable<float> dcapostopv{"dcapostopv", .1, "min DCA Pos To PV (cm)"};
    Configurable<float> v0radius{"v0radius", 1.2, "minimum V0 radius (cm)"};
    Configurable<float> v0radiusMax{"v0radiusMax", 1E5, "maximum V0 radius (cm)"};

    // Additional selection on the AP plot (exclusive for K0Short)
    // original equation: lArmPt*5>TMath::Abs(lArmAlpha)
    Configurable<float> armPodCut{"armPodCut", 5.0f, "pT * (cut) > |alpha|, AP cut. Negative: no cut"};
    Configurable<float> armAlphaSel{"armAlphaSel", 0.45f, "Armenteros alpha selection (Gammas)"};
    Configurable<float> qTSel{"qTSel", 0.01f, "Armenteros qT select (Gammas)"};

    // Selection
    Configurable<bool> applyInvMassSel{"applyInvMassSel", false, "Select V0s close to the Inv. mass value"};
    Configurable<float> dMassSel{"dMassSel", 0.01f, "Invariant mass selection"};
    Configurable<float> dMassSelG{"dMassSelG", 0.1f, "Inv mass selection gammas"};
    Configurable<float> dMassGcut{"dMassGcut", 0.1f, "Inv mass selection gammas"};

    // PID (TPC/TOF)
    Configurable<float> tpcPidNsigmaCut{"tpcPidNsigmaCut", 5, "tpcPidNsigmaCut"};
    Configurable<float> maxPiTOFBeta{"maxPiTOFBeta", 0.00005, "Maximum beta TOF selection"};
    Configurable<float> maxElTOFBeta{"maxElTOFBeta", 0.1, "Maximum beta TOF selection"};
    Configurable<bool> applyTPCTOFCombinedCut{"applyTPCTOFCombinedCut", false, " Apply geometrical cut ? "};

    // Phi cut
    Configurable<bool> applyPhiCut{"applyPhiCut", false, "Apply geometrical cut?"};
    Configurable<bool> applyEtaCal{"applyEtaCal", false, "Apply eta calibration?"};
    Configurable<bool> usePinPhiSelection{"usePinPhiSelection", true, "Uses Phi selection as a function of P or Pt?"};
    Configurable<bool> applyNclSel{"applyNclSel", false, "Apply Min. found Ncl in TPC?"};
  } v0Selections;

  // Configurables Event Selection
  Configurable<bool> isNoCollInTimeRangeStrict{"isNoCollInTimeRangeStrict", true, "use isNoCollInTimeRangeStrict?"};
  Configurable<bool> isNoCollInTimeRangeStandard{"isNoCollInTimeRangeStandard", false, "use isNoCollInTimeRangeStandard?"};
  Configurable<bool> isNoCollInRofStrict{"isNoCollInRofStrict", true, "use isNoCollInRofStrict?"};
  Configurable<bool> isNoCollInRofStandard{"isNoCollInRofStandard", false, "use isNoCollInRofStandard?"};
  Configurable<bool> isNoHighMultCollInPrevRof{"isNoHighMultCollInPrevRof", true, "use isNoHighMultCollInPrevRof?"};
  Configurable<bool> isNoCollInTimeRangeNarrow{"isNoCollInTimeRangeNarrow", false, "use isNoCollInTimeRangeNarrow?"};
  Configurable<bool> isOccupancyCut{"isOccupancyCut", true, "Occupancy cut?"};
  Configurable<bool> isApplyFT0CbasedOccupancy{"isApplyFT0CbasedOccupancy", false, "T0C Occu cut"};
  Configurable<bool> applyNchSel{"applyNchSel", false, "Use mid-rapidity-based Nch selection"};
  Configurable<bool> skipRecoColGTOne{"skipRecoColGTOne", true, "Remove collisions if reconstructed more than once"};
  Configurable<std::string> detector4Calibration{"detector4Calibration", "T0M", "Detector for nSigma-Nch rejection"};

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
  ConfigurableAxis axisEta{"axisEta", {50, -1.0, 1.0}, "Eta axix"};
  ConfigurableAxis axisY{"axisY", {50, -1.0, 1.0}, "rapidity axix"};
  ConfigurableAxis axisArmAlpha{"axisArmAlpha", {200, -1.0, 1.0}, "Armenteros alpha"};
  ConfigurableAxis axisArmqT{"axisArmqT", {600, 0.0f, 0.3f}, "Armenteros qT"};
  ConfigurableAxis axisK0Mass{"axisK0Mass", {200, 0.4f, 0.6f}, "Mass K0Short"};
  ConfigurableAxis axisLambdaMass{"axisLambdaMass", {200, 1.101f, 1.131f}, "Mass Lambda"};
  ConfigurableAxis axisGammaMass{"axisGammaMass", {200, 0.0f, 0.5f}, "Mass Gamma"};
  ConfigurableAxis axisNsigmaTPC{"axisNsigmaTPC", {200, -10.0f, 10.0f}, "N sigma TPC"};
  ConfigurableAxis axisdEdx{"axisdEdx", {140, 20.0, 160.0}, "dEdx binning"};
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
  Configurable<std::string> pathPhiCutHigh{"pathPhiCutHigh", "Users/o/omvazque/PhiCut/OO/Global/High", "base path to the ccdb object"};
  Configurable<std::string> pathPhiCutLow{"pathPhiCutLow", "Users/o/omvazque/PhiCut/OO/Global/Low", "base path to the ccdb object"};
  Configurable<int64_t> ccdbNoLaterThan{"ccdbNoLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};

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
    Centrality,
    VtxZ,
    NchSel
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
    bool isCalLoaded = false;
  } etaCal;

  TrackSelection trkSelDaugthers;
  TrackSelection trkSelGlobal;
  TrackSelection trkSelDaugthersV0s()
  {
    TrackSelection selectedTracks;
    selectedTracks.SetEtaRange(-0.8f, 0.8f);
    selectedTracks.SetMinNCrossedRowsTPC(70);
    return selectedTracks;
  }

  int currentRunNumberNchSel;
  int currentRunNumberPhiSel;
  void init(InitContext const&)
  {
    currentRunNumberNchSel = -1;
    currentRunNumberPhiSel = -1;
    trkSelDaugthers = trkSelDaugthersV0s();
    trkSelGlobal = getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSibAny, TrackSelection::GlobalTrackRun3DCAxyCut::Default);

    // define axes you want to use
    const std::string titlePorPt{v0Selections.usePinPhiSelection ? "#it{p} (GeV/#it{c})" : "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisZpos{48, -12., 12., "Vtx_{z} (cm)"};
    const AxisSpec axisEvent{15, 0.5, 15.5, ""};
    const AxisSpec axisNcl{161, -0.5, 160.5, "#it{N}_{cl} TPC"};
    const AxisSpec axisPt{binsPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisPtV0s{binsPtV0s, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisPtNcl{binsPtNcl, Form("%s", titlePorPt.data())};
    const AxisSpec axisXPhiCut{binsPtPhiCut, Form("%s", titlePorPt.data())};
    const AxisSpec axisCent{binsCent, "T0C centrality"};
    // const char* endingEta[kNEtaHists] = {"02", "24", "46", "68"};
    // const char* latexEta[kNEtaHists] = {"0<|#eta|<0.2", "0.2<#eta<0.4", "0.4<#eta<0.6", "0.6<#eta<0.8"};
    const char* endingEta[kNEtaHists] = {"86", "64", "42", "20", "02", "24", "46", "68"};
    const char* latexEta[kNEtaHists] = {"-0.8<#eta<-0.6", "-0.6<#eta<-0.4", "-0.4<#eta<-0.2", "-0.2<#eta<0", "0<#eta<0.2", "0.2<#eta<0.4", "0.4<#eta<0.6", "0.6<#eta<0.8"};

    registry.add("EventCounter", ";;Events", kTH1F, {axisEvent});

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
    x->SetBinLabel(12, "Cent. Sel.");
    x->SetBinLabel(13, "VtxZ Sel.");
    x->SetBinLabel(14, "Nch Sel.");

    if (doprocessCalibrationAndV0s) {
      registry.add("zPos", ";;Entries;", kTH1F, {axisZpos});
      registry.add("T0Ccent", ";;Entries", kTH1F, {axisCent});
      registry.add("NchVsNPV", ";Nch; NPV;", kTH2F, {{{nBinsNPV, minNpv, maxNpv}, {nBinsNch, minNch, maxNch}}});
      registry.add("ExcludedEvtVsNch", ";Nch;Entries;", kTH1F, {{nBinsNch, minNch, maxNch}});
      registry.add("ExcludedEvtVsNPV", ";NPV;Entries;", kTH1F, {{nBinsNPV, minNpv, maxNpv}});

      registry.add("V0sCounter", ";V0 type; Entries;", kTH1F, {{4, 0.5, 4.5}});
      registry.add("dcaVsPt", "Primary particles;#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm);", kTH2F, {{{axisPt}, {40, -1.0, 1.0}}});
      registry.add("dcaDauVsPt", ";V0 #it{p}_{T} (GeV/#it{c});DCA_{xy} (cm) daughters;", kTH2F, {{{axisPtV0s}, {200, -10., 10.}}});
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

      registry.add("NclFindable", ";;Findable Ncl TPC", kTH2F, {axisPtNcl, axisNcl});
      registry.add("NclFindablep", ";;Findable #LTNcl#GT TPC", kTProfile, {axisPtNcl});
      registry.add("NclFoundVsPhipBeforeCut", Form("Found #LTNcl#GT TPC;%s (GeV/#it{c});#varphi", titlePorPt.data()), kTProfile2D, {{{axisXPhiCut}, {350, 0.0, 0.35}}});
      registry.add("NclFoundVsPhipAfterCut", Form("Found #LTNcl#GT TPC;%s (GeV/#it{c});#varphi", titlePorPt.data()), kTProfile2D, {{{axisXPhiCut}, {350, 0.0, 0.35}}});
      registry.add("NclVsEta", ";#eta;Found Ncl TPC", kTH2F, {{{axisEta}, {161, -0.5, 160.5}}});
      registry.add("NclVsEtap", ";#eta;Found #LTNcl#GT TPC", kTProfile, {axisEta});

      registry.add("NclVsEtaPiMIP", "MIP #pi^{+} + #pi^{-} (0.4 < #it{p} < 0.6 GeV/#it{c}, 40 < dE/dx < 60);#eta; Found Ncl TPC", kTH2F, {{{axisEta}, {161, -0.5, 160.5}}});
      registry.add("NclVsEtaPiMIPp", "MIP #pi^{+} + #pi^{-} (0.4 < #it{p} < 0.6 GeV/#it{c}, 40 < dE/dx < 60);#eta; Found #LTNcl#GT TPC", kTProfile, {axisEta});
      registry.add("NclVsEtaPiV0", ";#eta; Found Ncl TPC", kTH2F, {axisEta, axisNcl});
      registry.add("NclVsEtaPiV0p", ";#eta; Found #LTNcl#GT TPC", kTProfile, {axisEta});
      registry.add("NclVsEtaPrV0", ";#eta; Found Ncl TPC", kTH2F, {axisEta, axisNcl});
      registry.add("NclVsEtaPrV0p", ";#eta; Found #LTNcl#GT TPC", kTProfile, {axisEta});
      registry.add("NclVsEtaElV0", ";#eta; Found Ncl TPC", kTH2F, {axisEta, axisNcl});
      registry.add("NclVsEtaElV0p", ";#eta; Found #LTNcl#GT TPC", kTProfile, {axisEta});

      registry.add("EtaVsYK0s", ";#eta;#it{y};", kTH2F, {axisEta, axisY});
      registry.add("EtaVsYPiL", ";#eta;#it{y};", kTH2F, {axisEta, axisY});
      registry.add("EtaVsYPrL", ";#eta;#it{y};", kTH2F, {axisEta, axisY});
      registry.add("EtaVsYPiAL", ";#eta;#it{y};", kTH2F, {axisEta, axisY});
      registry.add("EtaVsYPrAL", ";#eta;#it{y};", kTH2F, {axisEta, axisY});
      registry.add("EtaVsYG", ";#eta;#it{y};", kTH2F, {axisEta, axisY});

      registry.add("TOFExpPi2TOF", ";Momentum (GeV/#it{c});t^{e}_{Exp}/t_{TOF}", kTH2F, {{{axisPtV0s}, {100, 0.2, 1.2}}});
      registry.add("TOFExpEl2TOF", ";Momentum (GeV/#it{c});t^{#pi}_{Exp}/t_{TOF}", kTH2F, {{{axisPtV0s}, {100, 0.2, 1.2}}});

      registry.add("betaVsMomentum", ";Momentum (GeV/#it{c}); #beta", kTH2F, {{{axisPtV0s}, {500, 0, 1.2}}});
      registry.add("dEdxVsEtaPiMIP", "MIP #pi^{+} + #pi^{-} (0.4 < #it{p} < 0.6 GeV/#it{c});#eta; dE/dx;", kTH2F, {{{axisEta}, {100, 0, 100}}});
      registry.add("dEdxVsEtaPiMIPp", "MIP #pi^{+} + #pi^{-} (0.4 < #it{p} < 0.6 GeV/#it{c});#eta; #LTdE/dx#GT", kTProfile, {axisEta});
      registry.add("dEdxVsEtaElMIP", "MIP e^{+} + e^{-} (0.3 < #it{p} < 0.4 GeV/#it{c});#eta; dE/dx;", kTH2F, {{{axisEta}, {100, 0, 100}}});
      registry.add("dEdxVsEtaElMIPp", "MIP e^{+} + e^{-} (0.3 < #it{p} < 0.4 GeV/#it{c});#eta; #LTdE/dx#GT", kTProfile, {axisEta});
      registry.add("dEdxVsEtaPiMIPV0", "MIP #pi^{+} + #pi^{-} (0.4 < #it{p} < 0.6 GeV/#it{c});#eta; dE/dx", kTH2F, {{{axisEta}, {100, 0, 100}}});
      registry.add("dEdxVsEtaPiMIPV0p", "MIP #pi^{+} + #pi^{-} (0.4 < #it{p} < 0.6 GeV/#it{c});#eta; #LTdE/dx#GT", kTProfile, {axisEta});
      registry.add("dEdxVsEtaElMIPV0", "e^{+} + e^{-} (0.4 < #it{p} < 0.6 GeV/#it{c});#eta; dE/dx", kTH2F, {{{axisEta}, {100, 0, 100}}});
      registry.add("dEdxVsEtaElMIPV0p", "e^{+} + e^{-} (0.4 < #it{p} < 0.6 GeV/#it{c});#eta; #LTdE/dx#GT", kTProfile, {axisEta});

      for (int i = 0; i < kNEtaHists; ++i) {
        dEdx[i] = registry.add<TH3>(Form("dEdx_%s", endingEta[i]), Form("%s;Momentum (GeV/#it{c});dE/dx;", latexEta[i]), kTH3F, {axisPt, axisdEdx, axisCent});
        pTVsP[i] = registry.add<TH2>(Form("pTVsP_%s", endingEta[i]), Form("%s;Momentum (GeV/#it{c});#it{p}_{T} (GeV/#it{c});", latexEta[i]), kTH2F, {axisPt, axisPt});
        dEdxPiV0[i] = registry.add<TH3>(Form("dEdxPiV0_%s", endingEta[i]), Form("#pi^{+} + #pi^{-}, %s;Momentum (GeV/#it{c});dE/dx;", latexEta[i]), kTH3F, {axisPtV0s, axisdEdx, axisCent});
        dEdxPrV0[i] = registry.add<TH3>(Form("dEdxPrV0_%s", endingEta[i]), Form("p + #bar{p}, %s;Momentum (GeV/#it{c});dE/dx;", latexEta[i]), kTH3F, {axisPtV0s, axisdEdx, axisCent});
        dEdxElV0[i] = registry.add<TH3>(Form("dEdxElV0_%s", endingEta[i]), Form("e^{+} + e^{-}, %s;Momentum (GeV/#it{c});dE/dx;", latexEta[i]), kTH3F, {axisPtV0s, axisdEdx, axisCent});
        dEdxPiTOF[i] = registry.add<TH3>(Form("dEdxPiTOF_%s", endingEta[i]), Form("#pi^{+} + #pi^{-}, %s;Momentum (GeV/#it{c});dE/dx;", latexEta[i]), kTH3F, {axisPtV0s, axisdEdx, axisCent});
        nClVsP[i] = registry.add<TH2>(Form("NclFound_%s", endingEta[i]), Form("%s;;Found Ncl TPC", latexEta[i]), kTH2F, {axisPtNcl, axisNcl});
        nClVsPElV0[i] = registry.add<TH2>(Form("NclFoundElV0_%s", endingEta[i]), Form("%s;;Found Ncl TPC", latexEta[i]), kTH2F, {axisPtNcl, axisNcl});
        nClVsPPiV0[i] = registry.add<TH2>(Form("NclFoundPiV0_%s", endingEta[i]), Form("%s;;Found Ncl TPC", latexEta[i]), kTH2F, {axisPtNcl, axisNcl});
        nClVsPPrV0[i] = registry.add<TH2>(Form("NclFoundPrV0_%s", endingEta[i]), Form("%s;;Found Ncl TPC", latexEta[i]), kTH2F, {axisPtNcl, axisNcl});
        nClVsPp[i] = registry.add<TProfile>(Form("NclFoundp_%s", endingEta[i]), Form("%s;;Found #LT#it{N}_{cl}#GT TPC", latexEta[i]), kTProfile, {axisPtNcl});
        nClVsPpElV0[i] = registry.add<TProfile>(Form("NclFoundpElV0_%s", endingEta[i]), Form("%s;;Found #LT#it{N}_{cl}#GT TPC", latexEta[i]), kTProfile, {axisPtNcl});
        nClVsPpPiV0[i] = registry.add<TProfile>(Form("NclFoundpPiV0_%s", endingEta[i]), Form("%s;;Found #LT#it{N}_{cl}#GT TPC", latexEta[i]), kTProfile, {axisPtNcl});
        nClVsPpPrV0[i] = registry.add<TProfile>(Form("NclFoundpPrV0_%s", endingEta[i]), Form("%s;;Found #LT#it{N}_{cl}#GT TPC", latexEta[i]), kTProfile, {axisPtNcl});
      }
    }

    LOG(info) << "\tccdbNoLaterThan=" << ccdbNoLaterThan.value;
    LOG(info) << "\tapplyNchSel=" << applyNchSel.value;
    LOG(info) << "\tdetector4Calibration=" << detector4Calibration.value;
    LOG(info) << "\tv0TypeSelection=" << static_cast<int>(v0Selections.v0TypeSelection);
    LOG(info) << "\tminPt=" << v0Selections.minPt;
    LOG(info) << "\tmaxPt=" << v0Selections.maxPt;
    LOG(info) << "\tqTSel=" << v0Selections.qTSel;
    LOG(info) << "\tarmAlphaSel=" << v0Selections.armAlphaSel;
    LOG(info) << "\tapplyTPCTOFCombinedCut=" << v0Selections.applyTPCTOFCombinedCut;
    LOG(info) << "\tapplyInvMassSel=" << v0Selections.applyInvMassSel;
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
      LOG(info) << "\t pathMeanNch=" << pathMeanNch.value;
      LOG(info) << "\t pathSigmaNch=" << pathSigmaNch.value;
    }

    if (v0Selections.applyPhiCut) {
      LOG(info) << "\tLoading Phi cut!";
      LOG(info) << "\t pathPhiCutLow=" << pathPhiCutLow.value;
      LOG(info) << "\t pathPhiCutHigh=" << pathPhiCutHigh.value;
    }

    if (v0Selections.applyEtaCal) {
      LOG(info) << "\tLoading Eta Cal!";
      LOG(info) << "\t pathEtaCal=" << pathEtaCal.value;
      loadEtaCalibration();
    }

    if (v0Selections.applyNclSel)
      LOG(info) << "\t minNclFound=" << v0Selections.minNclFound;
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
    registry.fill(HIST("T0Ccent"), collision.centFT0C());
    const float centrality{collision.centFT0C()};

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

    for (const auto& track : tracks) {

      if (!trkSelGlobal.IsSelected(track))
        continue;

      if (track.pt() < v0Selections.minPt || track.pt() > v0Selections.maxPt)
        continue;

      const float momentum{track.p()};
      const float pt{track.pt()};
      const float phi{track.phi()};
      const float eta{track.eta()};
      float dedx{track.tpcSignal()};
      const int charge{track.sign()};
      const float pOrPt{v0Selections.usePinPhiSelection ? momentum : pt};
      const int16_t nclFound{track.tpcNClsFound()};

      if (v0Selections.applyNclSel && nclFound < v0Selections.minNclFound)
        continue;

      float phiPrime{phi};
      phiPrimeFunc(phiPrime, magField, charge);
      registry.fill(HIST("NclFoundVsPhipBeforeCut"), pOrPt, phiPrime, track.tpcNClsFound());

      if (v0Selections.applyPhiCut) {
        if (!passesPhiSelection(pOrPt, phiPrime))
          continue;
      }

      if (v0Selections.applyEtaCal) {
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
        registry.fill(HIST("NclVsEtaPiMIP"), eta, track.tpcNClsFound());
        registry.fill(HIST("NclVsEtaPiMIPp"), eta, track.tpcNClsFound());
      }

      if (momentum > kMinPElMIP && momentum < kMaxPElMIP && dedx > kMindEdxMIPPlateau && dedx < kMaxdEdxMIPPlateau) {
        registry.fill(HIST("dEdxVsEtaElMIP"), eta, dedx);
        registry.fill(HIST("dEdxVsEtaElMIPp"), eta, dedx);
      }

      dEdx[indexEta]->Fill(momentum, dedx, centrality);
      pTVsP[indexEta]->Fill(momentum, pt);
      nClVsP[indexEta]->Fill(pOrPt, track.tpcNClsFound());
      nClVsPp[indexEta]->Fill(pOrPt, track.tpcNClsFound());
      registry.fill(HIST("dcaVsPt"), pt, track.dcaXY());
      registry.fill(HIST("NclVsEta"), eta, track.tpcNClsFound());
      registry.fill(HIST("NclVsEtap"), eta, track.tpcNClsFound());
      // registry.fill(HIST("NclFoundp"), pOrPt, track.tpcNClsFound());
      registry.fill(HIST("NclFindable"), pOrPt, track.tpcNClsFindable());
      registry.fill(HIST("NclFindablep"), pOrPt, track.tpcNClsFindable());
      registry.fill(HIST("NclFoundVsPhipAfterCut"), pOrPt, phiPrime, track.tpcNClsFound());

      if (track.hasTOF() && track.goodTOFMatch()) {
        const float tTOF{track.tofSignal()};
        const float trkLength{track.length()};
        const float tExpPiTOF{track.tofExpSignalPi(tTOF)};
        const float tExpElTOF{track.tofExpSignalEl(tTOF)};
        // const float dTOFPi{tTOF - tExpPiTOF - colTime};
        // const float dTOFEl{tTOF - tExpElTOF - colTime};

        if (trkLength > kZero && tTOF > kZero) {
          registry.fill(HIST("betaVsMomentum"), momentum, track.beta());
          registry.fill(HIST("TOFExpPi2TOF"), momentum, tExpPiTOF / tTOF);
          registry.fill(HIST("TOFExpEl2TOF"), momentum, tExpElTOF / tTOF);

          if (std::abs((tExpPiTOF / tTOF) - kOne) < v0Selections.maxPiTOFBeta) {
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

      if (v0Selections.applyNclSel && (posNclFound < v0Selections.minNclFound || negNclFound < v0Selections.minNclFound))
        continue;

      if (v0Selections.applyEtaCal) {
        const double dedxCal{etaCal.pEtaCal->GetBinContent(etaCal.pEtaCal->FindBin(posTrkEta))};
        if (dedxCal > kMindEdxMIP && dedxCal < kMaxdEdxMIP)
          posTrkdEdx *= (50.0 / dedxCal);
        else
          continue;
      }

      if (v0Selections.applyEtaCal) {
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

      const double dMassK0s{std::abs(v0.mK0Short() - o2::constants::physics::MassK0Short)};
      const double dMassL{std::abs(v0.mLambda() - o2::constants::physics::MassLambda0)};
      const double dMassAL{std::abs(v0.mAntiLambda() - o2::constants::physics::MassLambda0)};
      const double dMassG{std::abs(v0.mGamma() - o2::constants::physics::MassGamma)};

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

      if (passesTopoSel) {
        registry.fill(HIST("ArmAfterTopoSel"), alpha, qT);
        registry.fill(HIST("dcaDauVsPt"), v0.pt(), v0.dcapostopv());
        registry.fill(HIST("dcaDauVsPt"), v0.pt(), v0.dcanegtopv());

        if (v0Selections.applyInvMassSel) {                                                                                                               // apply Inv. Mass selection?
          if (dMassK0s < v0Selections.dMassSel && dMassL > v0Selections.dMassSel && dMassAL > v0Selections.dMassSel && dMassG > v0Selections.dMassSelG) { // Mass cut
            if (passesK0Selection(collision, v0)) {                                                                                                       // nSigma TPC and y cuts
              registry.fill(HIST("ArmK0NOSel"), alpha, qT);
              if (v0Selections.armPodCut * qT > std::abs(alpha)) { // Armenters selection
                registry.fill(HIST("V0sCounter"), V0sCounter::K0s);
                registry.fill(HIST("ArmK0"), alpha, qT);
                registry.fill(HIST("MassK0sVsPt"), v0.pt(), v0.mK0Short());
                registry.fill(HIST("nSigPiFromK0s"), posTrkPt, posTrack.tpcNSigmaPi());
                registry.fill(HIST("nSigPiFromK0s"), negTrkPt, negTrack.tpcNSigmaPi());
                registry.fill(HIST("NclVsEtaPiV0"), posTrkEta, posTrack.tpcNClsFound());
                registry.fill(HIST("NclVsEtaPiV0p"), posTrkEta, posTrack.tpcNClsFound());
                registry.fill(HIST("NclVsEtaPiV0"), negTrkEta, negTrack.tpcNClsFound());
                registry.fill(HIST("NclVsEtaPiV0p"), negTrkEta, negTrack.tpcNClsFound());
                nClVsPPiV0[posIndexEta]->Fill(posPorPt, posTrack.tpcNClsFound());
                nClVsPpPiV0[posIndexEta]->Fill(posPorPt, posTrack.tpcNClsFound());
                nClVsPPiV0[negIndexEta]->Fill(negPorPt, negTrack.tpcNClsFound());
                nClVsPpPiV0[negIndexEta]->Fill(negPorPt, negTrack.tpcNClsFound());
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
          }
        }

        if (v0Selections.applyInvMassSel) {
          if (dMassL < v0Selections.dMassSel && dMassK0s > v0Selections.dMassSel && dMassG > v0Selections.dMassSelG) {
            if (passesLambdaSelection(collision, v0)) {
              registry.fill(HIST("V0sCounter"), V0sCounter::Lambda);
              registry.fill(HIST("ArmL"), alpha, qT);
              registry.fill(HIST("MassLVsPt"), v0.pt(), v0.mLambda());
              registry.fill(HIST("nSigPrFromL"), posTrkPt, posTrack.tpcNSigmaPr());
              registry.fill(HIST("nSigPiFromL"), negTrkPt, negTrack.tpcNSigmaPi());
              registry.fill(HIST("NclVsEtaPrV0"), posTrkEta, posTrack.tpcNClsFound());
              registry.fill(HIST("NclVsEtaPrV0p"), posTrkEta, posTrack.tpcNClsFound());
              registry.fill(HIST("NclVsEtaPiV0"), negTrkEta, negTrack.tpcNClsFound());
              registry.fill(HIST("NclVsEtaPiV0p"), negTrkEta, negTrack.tpcNClsFound());
              nClVsPPrV0[posIndexEta]->Fill(posPorPt, posTrack.tpcNClsFound());
              nClVsPpPrV0[posIndexEta]->Fill(posPorPt, posTrack.tpcNClsFound());
              nClVsPPiV0[negIndexEta]->Fill(negPorPt, negTrack.tpcNClsFound());
              nClVsPpPiV0[negIndexEta]->Fill(negPorPt, negTrack.tpcNClsFound());
              dEdxPrV0[posIndexEta]->Fill(posTrkP, posTrkdEdx, centrality);
              dEdxPiV0[negIndexEta]->Fill(negTrkP, negTrkdEdx, centrality);
            }
          }
        }

        if (v0Selections.applyInvMassSel && dMassAL < v0Selections.dMassSel && dMassK0s > v0Selections.dMassSel && dMassG > v0Selections.dMassSelG) {
          if (passesAntiLambdaSelection(collision, v0)) {
            registry.fill(HIST("V0sCounter"), V0sCounter::AntiLambda);
            registry.fill(HIST("ArmAL"), alpha, qT);
            registry.fill(HIST("MassALVsPt"), v0.pt(), v0.mAntiLambda());
            registry.fill(HIST("nSigPrFromAL"), negTrkPt, negTrack.tpcNSigmaPr());
            registry.fill(HIST("nSigPiFromAL"), posTrkPt, posTrack.tpcNSigmaPi());
            registry.fill(HIST("NclVsEtaPiV0"), posTrkEta, posTrack.tpcNClsFound());
            registry.fill(HIST("NclVsEtaPiV0p"), posTrkEta, posTrack.tpcNClsFound());
            registry.fill(HIST("NclVsEtaPrV0"), negTrkEta, negTrack.tpcNClsFound());
            registry.fill(HIST("NclVsEtaPrV0p"), negTrkEta, negTrack.tpcNClsFound());
            nClVsPPrV0[negIndexEta]->Fill(negPorPt, negTrack.tpcNClsFound());
            nClVsPpPrV0[negIndexEta]->Fill(negPorPt, negTrack.tpcNClsFound());
            nClVsPPiV0[posIndexEta]->Fill(posPorPt, posTrack.tpcNClsFound());
            nClVsPpPiV0[posIndexEta]->Fill(posPorPt, posTrack.tpcNClsFound());
            dEdxPrV0[negIndexEta]->Fill(negTrkP, negTrkdEdx, centrality);
            dEdxPiV0[posIndexEta]->Fill(posTrkP, posTrkdEdx, centrality);
          }
        }
      }

      if (v0Selections.applyInvMassSel && dMassK0s > v0Selections.dMassSel && dMassL > v0Selections.dMassSel && dMassAL > v0Selections.dMassSel && dMassG < v0Selections.dMassGcut) {
        if (passesGammaSelection(collision, v0)) {
          if (std::abs(alpha) < v0Selections.armAlphaSel && qT < v0Selections.qTSel) {
            registry.fill(HIST("V0sCounter"), V0sCounter::Gamma);
            registry.fill(HIST("ArmG"), alpha, qT);
            registry.fill(HIST("MassGVsPt"), v0.pt(), v0.mGamma());
            registry.fill(HIST("nSigElFromG"), negTrkPt, negTrack.tpcNSigmaEl());
            registry.fill(HIST("nSigElFromG"), posTrkPt, posTrack.tpcNSigmaEl());
            registry.fill(HIST("NclVsEtaElV0"), posTrkEta, posTrack.tpcNClsFound());
            registry.fill(HIST("NclVsEtaElV0p"), posTrkEta, posTrack.tpcNClsFound());
            registry.fill(HIST("NclVsEtaElV0"), negTrkEta, negTrack.tpcNClsFound());
            registry.fill(HIST("NclVsEtaElV0p"), negTrkEta, negTrack.tpcNClsFound());
            nClVsPElV0[negIndexEta]->Fill(negPorPt, negTrack.tpcNClsFound());
            nClVsPpElV0[negIndexEta]->Fill(negPorPt, negTrack.tpcNClsFound());
            nClVsPElV0[posIndexEta]->Fill(posPorPt, posTrack.tpcNClsFound());
            nClVsPpElV0[posIndexEta]->Fill(posPorPt, posTrack.tpcNClsFound());
            registry.fill(HIST("dEdxVsEtaElMIPV0"), posTrkEta, posTrkdEdx);
            registry.fill(HIST("dEdxVsEtaElMIPV0p"), posTrkEta, posTrkdEdx);
            registry.fill(HIST("dEdxVsEtaElMIPV0"), negTrkEta, negTrkdEdx);
            registry.fill(HIST("dEdxVsEtaElMIPV0p"), negTrkEta, negTrkdEdx);
            dEdxElV0[posIndexEta]->Fill(posTrkP, posTrkdEdx, centrality);
            dEdxElV0[negIndexEta]->Fill(negTrkP, negTrkdEdx, centrality);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(PiKpRAA, processCalibrationAndV0s, "Process QA", true);

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

  // Daughters DCA selection
  template <typename T>
  bool passesDCASelectionDaughters(const T& v0)
  {

    bool isSelected{false};
    // const double ptPos{std::sqrt(std::pow(v0.pxpos(), 2.0) + std::pow(v0.pypos(), 2.0))};
    // const double ptNeg{std::sqrt(std::pow(v0.pxneg(), 2.0) + std::pow(v0.pyneg(), 2.0))};
    // const double dcaPtDepPos{0.0105 + 0.035 * std::pow(ptPos, -1.1)};
    // const double dcaPtDepNeg{0.0105 + 0.035 * std::pow(ptNeg, -1.1)};
    // const double dcaSelPos = std::max(0.1, dcaPtDepPos);
    // const double dcaSelNeg = std::max(0.1, dcaPtDepNeg);

    const double dcaPos{std::fabs(v0.dcapostopv())};
    const double dcaNeg{std::fabs(v0.dcanegtopv())};

    isSelected = dcaPos > v0Selections.dcapostopv && dcaNeg > v0Selections.dcanegtopv ? true : false;
    return isSelected;
  }

  template <typename T>
  bool passesTrackSelectionDaughters(const T& track)
  {

    bool isSelected = trkSelDaugthers.IsSelected(track) ? true : false;

    return isSelected;
  }

  // V0 topological selection
  template <typename T>
  bool passesV0TopologicalSelection(const T& v0)
  {

    bool isSelected = v0.v0radius() > v0Selections.v0radius && v0.v0radius() < v0Selections.v0radiusMax && passesDCASelectionDaughters(v0) && v0.v0cosPA() > v0Selections.v0cospa && v0.dcaV0daughters() < v0Selections.dcav0dau ? true : false;

    return isSelected;
  }

  template <typename C, typename T>
  bool passesK0Selection(const C& collision, const T& v0)
  {
    // Selection on rapiditty, proper lifetime, and Nsigma Pion

    const auto& posTrack = v0.template posTrack_as<TracksFull>();
    const auto& negTrack = v0.template negTrack_as<TracksFull>();

    const float posTPCNsigma{std::fabs(posTrack.tpcNSigmaPi())};
    const float negTPCNsigma{std::fabs(negTrack.tpcNSigmaPi())};
    const float posTOFNsigma{std::fabs(posTrack.tofNSigmaPi())};
    const float negTOFNsigma{std::fabs(negTrack.tofNSigmaPi())};
    const double posRadiusNsigma{std::sqrt(std::pow(posTPCNsigma, 2.) + std::pow(posTOFNsigma, 2.))};
    const double negRadiusNsigma{std::sqrt(std::pow(negTPCNsigma, 2.) + std::pow(negTOFNsigma, 2.))};

    registry.fill(HIST("EtaVsYK0s"), negTrack.eta(), v0.yK0Short());
    registry.fill(HIST("EtaVsYK0s"), posTrack.eta(), v0.yK0Short());

    bool isSelected{false};
    if (v0Selections.applyTPCTOFCombinedCut)
      isSelected = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short < lifetimecut->get("lifetimecutK0S") && std::abs(v0.yK0Short()) < v0Selections.rapidityCut && posTrack.hasTOF() && negTrack.hasTOF() && posRadiusNsigma < v0Selections.tpcPidNsigmaCut && negRadiusNsigma < v0Selections.tpcPidNsigmaCut ? true : false;
    else
      isSelected = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short < lifetimecut->get("lifetimecutK0S") && std::abs(v0.yK0Short()) < v0Selections.rapidityCut && posTPCNsigma < v0Selections.tpcPidNsigmaCut && negTPCNsigma < v0Selections.tpcPidNsigmaCut ? true : false;

    return isSelected;
  }

  template <typename C, typename T>
  bool passesLambdaSelection(const C& collision, const T& v0)
  {
    // Selection on rapiditty, proper lifetime, and Nsigma Pion

    const auto& posTrack = v0.template posTrack_as<TracksFull>();
    const auto& negTrack = v0.template negTrack_as<TracksFull>();

    const float posTPCNsigma{std::fabs(posTrack.tpcNSigmaPr())};
    const float negTPCNsigma{std::fabs(negTrack.tpcNSigmaPi())};
    const float posTOFNsigma{std::fabs(posTrack.tofNSigmaPr())};
    const float negTOFNsigma{std::fabs(negTrack.tofNSigmaPi())};
    const double posRadiusNsigma{std::sqrt(std::pow(posTPCNsigma, 2.) + std::pow(posTOFNsigma, 2.))};
    const double negRadiusNsigma{std::sqrt(std::pow(negTPCNsigma, 2.) + std::pow(negTOFNsigma, 2.))};

    registry.fill(HIST("EtaVsYPiL"), negTrack.eta(), v0.yLambda());
    registry.fill(HIST("EtaVsYPrL"), posTrack.eta(), v0.yLambda());

    bool isSelected{false};
    if (v0Selections.applyTPCTOFCombinedCut)
      isSelected = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0 < lifetimecut->get("lifetimecutLambda") && std::abs(v0.yLambda()) < v0Selections.rapidityCut && posTrack.hasTOF() && negTrack.hasTOF() && posRadiusNsigma < v0Selections.tpcPidNsigmaCut && negRadiusNsigma < v0Selections.tpcPidNsigmaCut ? true : false;
    else
      isSelected = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0 < lifetimecut->get("lifetimecutLambda") && std::abs(v0.yLambda()) < v0Selections.rapidityCut && posTPCNsigma < v0Selections.tpcPidNsigmaCut && negTPCNsigma < v0Selections.tpcPidNsigmaCut ? true : false;

    return isSelected;
  }

  template <typename C, typename T>
  bool passesAntiLambdaSelection(const C& collision, const T& v0)
  {
    // Selection on rapiditty, proper lifetime, and Nsigma Pion

    const auto& posTrack = v0.template posTrack_as<TracksFull>();
    const auto& negTrack = v0.template negTrack_as<TracksFull>();

    const float posTPCNsigma{std::fabs(posTrack.tpcNSigmaPi())};
    const float negTPCNsigma{std::fabs(negTrack.tpcNSigmaPr())};
    const float posTOFNsigma{std::fabs(posTrack.tofNSigmaPi())};
    const float negTOFNsigma{std::fabs(negTrack.tofNSigmaPr())};
    const double posRadiusNsigma{std::sqrt(std::pow(posTPCNsigma, 2.) + std::pow(posTOFNsigma, 2.))};
    const double negRadiusNsigma{std::sqrt(std::pow(negTPCNsigma, 2.) + std::pow(negTOFNsigma, 2.))};

    registry.fill(HIST("EtaVsYPiAL"), posTrack.eta(), v0.yLambda());
    registry.fill(HIST("EtaVsYPrAL"), negTrack.eta(), v0.yLambda());

    bool isSelected{false};
    if (v0Selections.applyTPCTOFCombinedCut)
      isSelected = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0 < lifetimecut->get("lifetimecutLambda") && std::abs(v0.yLambda()) < v0Selections.rapidityCut && posTrack.hasTOF() && negTrack.hasTOF() && posRadiusNsigma < v0Selections.tpcPidNsigmaCut && negRadiusNsigma < v0Selections.tpcPidNsigmaCut ? true : false;
    else
      isSelected = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0 < lifetimecut->get("lifetimecutLambda") && std::abs(v0.yLambda()) < v0Selections.rapidityCut && posTPCNsigma < v0Selections.tpcPidNsigmaCut && negTPCNsigma < v0Selections.tpcPidNsigmaCut ? true : false;

    return isSelected;
  }

  template <typename C, typename T>
  bool passesGammaSelection(const C& /*collision*/, const T& v0)
  {
    const auto& posTrack = v0.template posTrack_as<TracksFull>();
    const auto& negTrack = v0.template negTrack_as<TracksFull>();

    const float posTPCNsigma{std::fabs(posTrack.tpcNSigmaEl())};
    const float negTPCNsigma{std::fabs(negTrack.tpcNSigmaEl())};
    const float posTOFNsigma{std::fabs(posTrack.tofNSigmaEl())};
    const float negTOFNsigma{std::fabs(negTrack.tofNSigmaEl())};
    const double posRadiusNsigma{std::sqrt(std::pow(posTPCNsigma, 2.) + std::pow(posTOFNsigma, 2.))};
    const double negRadiusNsigma{std::sqrt(std::pow(negTPCNsigma, 2.) + std::pow(negTOFNsigma, 2.))};
    const float yGamma = RecoDecay::y(std::array{v0.px(), v0.py(), v0.pz()}, o2::constants::physics::MassGamma);

    registry.fill(HIST("EtaVsYG"), negTrack.eta(), yGamma);
    registry.fill(HIST("EtaVsYG"), posTrack.eta(), yGamma);

    if (!(std::abs(yGamma) < v0Selections.rapidityCut))
      return false;

    bool isSelected{false};
    if (v0Selections.applyTPCTOFCombinedCut)
      isSelected = posTrack.hasTOF() && negTrack.hasTOF() && posRadiusNsigma < v0Selections.tpcPidNsigmaCut && negRadiusNsigma < v0Selections.tpcPidNsigmaCut ? true : false;
    else
      isSelected = posTPCNsigma < v0Selections.tpcPidNsigmaCut && negTPCNsigma < v0Selections.tpcPidNsigmaCut ? true : false;

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

    bool isSelected{false};
    if (phiCut.isPhiCutLoaded) {
      const int binLow{phiCut.hPhiCutLow->FindBin(pt)};
      const int binHigh{phiCut.hPhiCutHigh->FindBin(pt)};
      const double phiCutLow{phiCut.hPhiCutLow->GetBinContent(binLow)};
      const double phiCutHigh{phiCut.hPhiCutHigh->GetBinContent(binHigh)};
      if (phi < phiCutLow || phi > phiCutHigh)
        isSelected = true;
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

    if (!col.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    registry.fill(HIST("EventCounter"), EvCutLabel::NoSameBunchPileup);

    if (!col.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    registry.fill(HIST("EventCounter"), EvCutLabel::IsGoodZvtxFT0vsPV);

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
    }
    registry.fill(HIST("EventCounter"), EvCutLabel::OccuCut);

    if (col.centFT0C() < minT0CcentCut || col.centFT0C() > maxT0CcentCut) {
      return false;
    }
    registry.fill(HIST("EventCounter"), EvCutLabel::Centrality);

    // Z-vertex position cut
    if (std::fabs(col.posZ()) > posZcut) {
      return false;
    }
    registry.fill(HIST("EventCounter"), EvCutLabel::VtxZ);

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
      if (etaCal.pEtaCal == nullptr) {
        LOGF(fatal, "Could not load pEtaCal from %s", pathEtaCal.value.c_str());
      }
    }

    if (etaCal.pEtaCal)
      etaCal.isCalLoaded = true;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PiKpRAA>(cfgc)};
}
