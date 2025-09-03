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

/// \file uccZdc.cxx
///
/// \brief task for analysis of UCC with the ZDC
/// \author Omar Vazquez (omar.vazquez.rueda@cern.ch)
/// \since January 29, 2025

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/MathConstants.h"
#include "CommonConstants/ZDCConstants.h"
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
#include <TRandom3.h>
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

using namespace std;
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::evsel;
using namespace o2::constants::physics;
using namespace o2::constants::math;

namespace o2::aod
{
using ColEvSels = soa::Join<aod::Collisions, aod::EvSels, aod::FT0MultZeqs, o2::aod::CentFT0Cs, aod::TPCMults, o2::aod::BarrelMults>;
using BCsRun3 = soa::Join<aod::BCsWithTimestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
using TracksSel = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCovIU, aod::TrackSelection, aod::TracksDCA>;
using SimCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, o2::aod::CentFT0Cs>;
using SimTracks = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels>;
} // namespace o2::aod

static constexpr int kSizeBootStrapEnsemble{8};

std::array<std::shared_ptr<TH1>, kSizeBootStrapEnsemble> hPoisson{};
std::array<std::shared_ptr<TH2>, kSizeBootStrapEnsemble> hNchVsT0M{};
std::array<std::shared_ptr<TH2>, kSizeBootStrapEnsemble> hNchVsV0A{};
std::array<std::shared_ptr<TH2>, kSizeBootStrapEnsemble> hNchVsZN{};

std::array<std::shared_ptr<TProfile2D>, kSizeBootStrapEnsemble> pNchVsOneParCorrVsZN{};
std::array<std::shared_ptr<TProfile2D>, kSizeBootStrapEnsemble> pNchVsTwoParCorrVsZN{};
std::array<std::shared_ptr<TProfile2D>, kSizeBootStrapEnsemble> pNchVsThreeParCorrVsZN{};

std::array<std::shared_ptr<TProfile2D>, kSizeBootStrapEnsemble> pNchVsOneParCorrVsT0M{};
std::array<std::shared_ptr<TProfile2D>, kSizeBootStrapEnsemble> pNchVsTwoParCorrVsT0M{};
std::array<std::shared_ptr<TProfile2D>, kSizeBootStrapEnsemble> pNchVsThreeParCorrVsT0M{};

std::array<std::shared_ptr<TProfile2D>, kSizeBootStrapEnsemble> pNchVsOneParCorrVsV0A{};
std::array<std::shared_ptr<TProfile2D>, kSizeBootStrapEnsemble> pNchVsTwoParCorrVsV0A{};
std::array<std::shared_ptr<TProfile2D>, kSizeBootStrapEnsemble> pNchVsThreeParCorrVsV0A{};

std::array<std::shared_ptr<TProfile>, kSizeBootStrapEnsemble> pOneParCorrVsNch{};
std::array<std::shared_ptr<TProfile>, kSizeBootStrapEnsemble> pTwoParCorrVsNch{};
std::array<std::shared_ptr<TProfile>, kSizeBootStrapEnsemble> pThreeParCorrVsNch{};

std::array<std::shared_ptr<TH1>, kSizeBootStrapEnsemble> hPoissonMC{};
std::array<std::shared_ptr<TH1>, kSizeBootStrapEnsemble> hNchGen{};
std::array<std::shared_ptr<TH1>, kSizeBootStrapEnsemble> hNch{};

std::array<std::shared_ptr<TProfile>, kSizeBootStrapEnsemble> pOneParCorrVsNchGen{};
std::array<std::shared_ptr<TProfile>, kSizeBootStrapEnsemble> pTwoParCorrVsNchGen{};
std::array<std::shared_ptr<TProfile>, kSizeBootStrapEnsemble> pThreeParCorrVsNchGen{};

struct UccZdc {

  static constexpr float kCollEnergy{2.68};
  static constexpr float kZero{0.};
  static constexpr float kOne{1.};
  static constexpr float kMinCharge{3.f};

  // Configurables Event Selection
  Configurable<bool> isNoCollInTimeRangeStrict{"isNoCollInTimeRangeStrict", true, "use isNoCollInTimeRangeStrict?"};
  Configurable<bool> isNoCollInTimeRangeStandard{"isNoCollInTimeRangeStandard", false, "use isNoCollInTimeRangeStandard?"};
  Configurable<bool> isNoCollInRofStrict{"isNoCollInRofStrict", true, "use isNoCollInRofStrict?"};
  Configurable<bool> isNoCollInRofStandard{"isNoCollInRofStandard", false, "use isNoCollInRofStandard?"};
  Configurable<bool> isNoHighMultCollInPrevRof{"isNoHighMultCollInPrevRof", true, "use isNoHighMultCollInPrevRof?"};
  Configurable<bool> isNoCollInTimeRangeNarrow{"isNoCollInTimeRangeNarrow", false, "use isNoCollInTimeRangeNarrow?"};
  Configurable<bool> isOccupancyCut{"isOccupancyCut", true, "Occupancy cut?"};
  Configurable<bool> isApplyFT0CbasedOccupancy{"isApplyFT0CbasedOccupancy", false, "T0C Occu cut"};
  Configurable<bool> isTDCcut{"isTDCcut", false, "Use TDC cut"};
  Configurable<bool> useMidRapNchSel{"useMidRapNchSel", true, "Use mid-rapidit Nch selection"};
  Configurable<bool> applyEff{"applyEff", true, "Apply track-by-track efficiency correction"};
  Configurable<bool> applyFD{"applyFD", false, "Apply track-by-track feed down correction"};
  Configurable<bool> correctNch{"correctNch", true, "Correct also Nch"};
  Configurable<bool> skipRecoColGTOne{"skipRecoColGTOne", true, "Remove collisions if reconstructed more than once"};
  Configurable<std::string> detector4Calibration{"detector4Calibration", "T0M", "Detector for nSigma-Nch rejection"};
  Configurable<std::string> detectorZDC{"detectorZDC", "ZN", "Detector for Cent. Selec. based on spectator neutrons"};

  // Event selection
  Configurable<float> posZcut{"posZcut", +10.0, "z-vertex position cut"};
  Configurable<float> minT0CcentCut{"minT0CcentCut", 0.0, "Min T0C Cent. cut"};
  Configurable<float> maxT0CcentCut{"maxT0CcentCut", 90.0, "Max T0C Cent. cut"};
  Configurable<float> nSigmaNchCut{"nSigmaNchCut", 1., "nSigma Nch selection"};
  Configurable<float> zemCut{"zemCut", 1000., "ZEM cut"};
  Configurable<float> tdcCut{"tdcCut", 1., "TDC cut"};
  Configurable<float> minOccCut{"minOccCut", 0., "min Occu cut"};
  Configurable<float> maxOccCut{"maxOccCut", 500., "max Occu cut"};
  Configurable<float> minNchSel{"minNchSel", 5., "min Nch Selection"};
  Configurable<float> evtFracMCcl{"evtFracMCcl", 0.5, "fraction of events for MC closure"};

  // Track-kinematics selection
  Configurable<float> minPt{"minPt", 0.1, "minimum pt of the tracks"};
  Configurable<float> maxPt{"maxPt", 3., "maximum pt of the tracks"};
  Configurable<float> maxPtSpectra{"maxPtSpectra", 50., "maximum pt of the tracks"};
  Configurable<float> minEta{"minEta", -0.8, "minimum eta"};
  Configurable<float> maxEta{"maxEta", +0.8, "maximum eta"};

  // Configurables, binning
  Configurable<float> minNch{"minNch", 0, "Min Nch (|eta|<0.8)"};
  Configurable<float> minZN{"minZN", -0.5, "Min ZN signal"};
  Configurable<float> minTdc{"minTdc", -15.0, "minimum TDC"};
  Configurable<float> arbScale{"arbScale", 100.0, "Scale factor for forward multiplicity"};

  Configurable<float> maxNch{"maxNch", 3000, "Max Nch (|eta|<0.8)"};
  Configurable<float> maxITSTrack{"maxITSTrack", 6000., "Min ITS tracks"};
  Configurable<float> maxAmpV0A{"maxAmpV0A", 2000, "Max FV0 amp"};
  Configurable<float> maxAmpFT0{"maxAmpFT0", 2500, "Max FT0 amp"};
  Configurable<float> maxAmpFT0A{"maxAmpFT0A", 200, "Max FT0 amp"};
  Configurable<float> maxAmpFT0C{"maxAmpFT0C", 60, "Max FT0 amp"};
  Configurable<float> maxZN{"maxZN", 150, "Max ZN signal"};
  Configurable<float> maxZP{"maxZP", 60, "Max ZP signal"};
  Configurable<float> maxZEM{"maxZEM", 2200, "Max ZEM signal"};
  Configurable<float> maxTdc{"maxTdc", 15.0, "maximum TDC"};

  Configurable<int> nBinsNch{"nBinsNch", 2501, "N bins Nch (|eta|<0.8)"};
  Configurable<int> nBinsITSTrack{"nBinsITSTrack", 2000, "N bins ITS tracks"};
  Configurable<int> nBinsAmpV0A{"nBinsAmpV0A", 100, "N bins V0A amp"};
  Configurable<int> nBinsAmpFT0{"nBinsAmpFT0", 100, "N bins FT0 amp"};
  Configurable<int> nBinsAmpFT0A{"nBinsAmpFT0A", 100, "N bins FT0A amp"};
  Configurable<int> nBinsAmpFT0C{"nBinsAmpFT0C", 100, "N bins FT0C amp"};
  Configurable<int> nBinsZN{"nBinsZN", 400, "N bins ZN"};
  Configurable<int> nBinsZP{"nBinsZP", 160, "N bins ZP"};
  Configurable<int> nBinsTDC{"nBinsTDC", 150, "nbinsTDC"};

  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.12}, "pT binning"};
  ConfigurableAxis binsCent{"binsCent", {VARIABLE_WIDTH, 0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.}, "T0C binning"};

  // CCDB paths
  Configurable<std::string> paTHEff{"paTHEff", "Users/o/omvazque/MCcorrection/perTimeStamp/TrackingEff", "base path to the ccdb object"};
  Configurable<std::string> paTHFD{"paTHFD", "Users/o/omvazque/MCcorrection/perTimeStamp/FeedDown", "base path to the ccdb object"};
  Configurable<std::string> paTHmeanNch{"paTHmeanNch", "Users/o/omvazque/FitMeanNch_9May2025", "base path to the ccdb object"};
  Configurable<std::string> paTHsigmaNch{"paTHsigmaNch", "Users/o/omvazque/FitSigmaNch_9May2025", "base path to the ccdb object"};
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
    Zdc,
    TZero,
    Tdc,
    Zem
  };

  Filter trackFilter = ((aod::track::eta > minEta) && (aod::track::eta < maxEta));

  // Apply Filters
  using TheFilteredTracks = soa::Filtered<o2::aod::TracksSel>;
  using TheFilteredSimTracks = soa::Filtered<o2::aod::SimTracks>;

  // Histograms: Data
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  Service<ccdb::BasicCCDBManager> ccdb;

  struct Config {
    TH2F* hEfficiency = nullptr;
    TH2F* hFeedDown = nullptr;
    bool correctionsLoaded = false;
  } cfg;

  struct NchConfig {
    TH1F* hMeanNch = nullptr;
    TH1F* hSigmaNch = nullptr;
    bool calibrationsLoaded = false;
  } cfgNch;

  int currentRunNumber;

  void init(InitContext const&)
  {
    currentRunNumber = -1;
    const char* tiT0A{"T0A (#times 1/100, 3.5 < #eta < 4.9)"};
    const char* tiT0C{"T0C (#times 1/100, -3.3 < #eta < -2.1)"};
    const char* tiT0M{"T0A+T0C (#times 1/100, -3.3 < #eta < -2.1 and 3.5 < #eta < 4.9)"};
    const char* tiNch{"#it{N}_{ch} (|#eta| < 0.8)"};
    const char* tiNPV{"#it{N}_{PV} (|#eta|<1)"};
    const char* tiV0A{"V0A (#times 1/100, 2.2 < #eta < 5)"};
    const char* tiZNs{"ZNA + ZNC"};
    const char* tiZPs{"ZPA + ZPC"};
    const char* tiPt{"#it{p}_{T} (GeV/#it{c})"};
    const char* tiOneParCorr{"#LT[#it{p}_{T}^{(1)}]#GT (GeV/#it{c})"};
    const char* tiTwoParCorr{"Two-Particle #it{p}_{T} correlation"};
    const char* tiThreeParCorr{"Three-Particle #it{p}_{T} correlation"};

    // define axes you want to use
    const AxisSpec axisZpos{48, -12., 12., "Vtx_{z} (cm)"};
    const AxisSpec axisEvent{18, 0.5, 18.5, ""};
    const AxisSpec axisEta{40, -1., +1., "#eta"};
    const AxisSpec axisPt{binsPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisDeltaPt{100, -1.0, +1.0, "#Delta(p_{T})"};
    const AxisSpec axisCent{binsCent, "T0C centrality"};
    const AxisSpec axisAmpCh{250, 0., 2500., "Amplitude of non-zero channels"};
    const AxisSpec axisEneCh{300, 0., 300., "Energy of non-zero channels"};

    registry.add("NchUncorrected", Form(";%s;Entries;", tiNch), kTH1F, {{nBinsNch, minNch, maxNch}});
    registry.add("hEventCounter", ";;Events", kTH1F, {axisEvent});
    registry.add("ExcludedEvtVsFT0M", Form(";%s;Entries;", tiT0M), kTH1F, {{nBinsAmpFT0, 0., maxAmpFT0}});
    registry.add("ExcludedEvtVsFV0A", Form(";%s;Entries;", tiT0M), kTH1F, {{nBinsAmpV0A, 0., maxAmpV0A}});
    registry.add("ExcludedEvtVsNch", Form(";%s;Entries;", tiNch), kTH1F, {{nBinsNch, minNch, maxNch}});
    registry.add("ExcludedEvtVsNPV", Form(";%s;Entries;", tiNPV), kTH1F, {{nBinsITSTrack, 0, maxITSTrack}});
    registry.add("Nch", Form(";%s;Entries;", tiNch), kTH1F, {{nBinsNch, minNch, maxNch}});
    registry.add("NchVsOneParCorr", Form(";%s;%s;", tiNch, tiOneParCorr), kTProfile, {{nBinsNch, minNch, maxNch}});
    registry.add("NchVsTwoParCorr", Form(";%s;#LT[#it{p}_{T}^{(2)}]#GT;", tiNch), kTProfile, {{nBinsNch, minNch, maxNch}});
    registry.add("NchVsThreeParCorr", Form(";%s;#LT[#it{p}_{T}^{(3)}]#GT;", tiNch), kTProfile, {{nBinsNch, minNch, maxNch}});

    auto hstat = registry.get<TH1>(HIST("hEventCounter"));
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
    x->SetBinLabel(13, "VtxZ cut");
    x->SetBinLabel(14, "has ZDC?");
    x->SetBinLabel(15, "has T0?");
    x->SetBinLabel(16, "Within TDC cut?");
    x->SetBinLabel(17, "Within ZEM cut?");

    if (doprocessZdcCollAss) {
      registry.add("NchVsT0M", Form(";%s;%s;", tiNch, tiT0M), kTH2F, {{{nBinsNch, minNch, maxNch}, {nBinsAmpFT0, 0., maxAmpFT0}}});
      registry.add("NchVsV0A", Form(";%s;%s;", tiNch, tiV0A), kTH2F, {{{nBinsNch, minNch, maxNch}, {nBinsAmpV0A, 0., maxAmpV0A}}});
      registry.add("NchVsZN", Form(";%s;%s;", tiNch, tiZNs), kTH2F, {{{nBinsNch, minNch, maxNch}, {nBinsZN, minZN, maxZN}}});
      registry.add("NchVsZP", Form(";%s;%s;", tiNch, tiZPs), kTH2F, {{{nBinsNch, minNch, maxNch}, {nBinsZP, minZN, maxZP}}});
      registry.add("NchVsZNVsPt", Form(";%s;%s;%s", tiNch, tiZNs, tiPt), kTH3F, {{{nBinsNch, minNch, maxNch}, {nBinsZN, minZN, maxZN}, {axisPt}}});

      registry.add("NchVsOneParCorrVsZN", Form(";%s;%s;%s", tiNch, tiZNs, tiOneParCorr), kTProfile2D, {{{nBinsNch, minNch, maxNch}, {nBinsZN, minZN, maxZN}}});
      registry.add("NchVsTwoParCorrVsZN", Form(";%s;%s;%s", tiNch, tiZNs, tiTwoParCorr), kTProfile2D, {{{nBinsNch, minNch, maxNch}, {nBinsZN, minZN, maxZN}}});
      registry.add("NchVsThreeParCorrVsZN", Form(";%s;%s;%s", tiNch, tiZNs, tiThreeParCorr), kTProfile2D, {{{nBinsNch, minNch, maxNch}, {nBinsZN, minZN, maxZN}}});

      registry.add("NchVsOneParCorrVsT0M", Form(";%s;%s;%s", tiNch, tiT0M, tiOneParCorr), kTProfile2D, {{{nBinsNch, minNch, maxNch}, {nBinsAmpFT0, 0., maxAmpFT0}}});
      registry.add("NchVsTwoParCorrVsT0M", Form(";%s;%s;%s", tiNch, tiT0M, tiTwoParCorr), kTProfile2D, {{{nBinsNch, minNch, maxNch}, {nBinsAmpFT0, 0., maxAmpFT0}}});
      registry.add("NchVsThreeParCorrVsT0M", Form(";%s;%s;%s", tiNch, tiT0M, tiThreeParCorr), kTProfile2D, {{{nBinsNch, minNch, maxNch}, {nBinsAmpFT0, 0., maxAmpFT0}}});

      registry.add("NchVsOneParCorrVsV0A", Form(";%s;%s;%s", tiNch, tiV0A, tiOneParCorr), kTProfile2D, {{{nBinsNch, minNch, maxNch}, {nBinsAmpV0A, 0., maxAmpV0A}}});
      registry.add("NchVsTwoParCorrVsV0A", Form(";%s;%s;%s", tiNch, tiV0A, tiTwoParCorr), kTProfile2D, {{{nBinsNch, minNch, maxNch}, {nBinsAmpV0A, 0., maxAmpV0A}}});
      registry.add("NchVsThreeParCorrVsV0A", Form(";%s;%s;%s", tiNch, tiV0A, tiThreeParCorr), kTProfile2D, {{{nBinsNch, minNch, maxNch}, {nBinsAmpV0A, 0., maxAmpV0A}}});

      for (int i = 0; i < kSizeBootStrapEnsemble; i++) {
        hNchVsZN[i] = registry.add<TH2>(Form("NchVsZN_Rep%d", i), Form(";%s;%s", tiNch, tiZNs), kTH2F, {{{nBinsNch, minNch, maxNch}, {nBinsZN, minZN, maxZN}}});
        hNchVsV0A[i] = registry.add<TH2>(Form("NchVsV0A_Rep%d", i), Form(";%s;%s", tiNch, tiV0A), kTH2F, {{{nBinsNch, minNch, maxNch}, {nBinsAmpV0A, 0., maxAmpV0A}}});
        hNchVsT0M[i] = registry.add<TH2>(Form("NchVsT0M_Rep%d", i), Form(";%s;%s", tiNch, tiT0M), kTH2F, {{{nBinsNch, minNch, maxNch}, {nBinsAmpFT0, 0., maxAmpFT0}}});
        hPoisson[i] = registry.add<TH1>(Form("Poisson_Rep%d", i), ";#it{k};Entries", kTH1F, {{11, -0.5, 10.5}});

        pNchVsOneParCorrVsZN[i] = registry.add<TProfile2D>(Form("NchVsOneParCorrVsZN_Rep%d", i), Form(";%s;%s;%s", tiNch, tiZNs, tiOneParCorr), kTProfile2D, {{{nBinsNch, minNch, maxNch}, {nBinsZN, minZN, maxZN}}});
        pNchVsTwoParCorrVsZN[i] = registry.add<TProfile2D>(Form("NchVsTwoParCorrVsZN_Rep%d", i), Form(";%s;%s;%s", tiNch, tiZNs, tiTwoParCorr), kTProfile2D, {{{nBinsNch, minNch, maxNch}, {nBinsZN, minZN, maxZN}}});
        pNchVsThreeParCorrVsZN[i] = registry.add<TProfile2D>(Form("NchVsThreeParCorrVsZN_Rep%d", i), Form(";%s;%s;%s", tiNch, tiZNs, tiThreeParCorr), kTProfile2D, {{{nBinsNch, minNch, maxNch}, {nBinsZN, minZN, maxZN}}});

        pNchVsOneParCorrVsT0M[i] = registry.add<TProfile2D>(Form("NchVsOneParCorrVsT0M_Rep%d", i), Form(";%s;%s;%s", tiNch, tiT0M, tiOneParCorr), kTProfile2D, {{{nBinsNch, minNch, maxNch}, {nBinsAmpFT0, 0., maxAmpFT0}}});
        pNchVsTwoParCorrVsT0M[i] = registry.add<TProfile2D>(Form("NchVsTwoParCorrVsT0M_Rep%d", i), Form(";%s;%s;%s", tiNch, tiT0M, tiTwoParCorr), kTProfile2D, {{{nBinsNch, minNch, maxNch}, {nBinsAmpFT0, 0., maxAmpFT0}}});
        pNchVsThreeParCorrVsT0M[i] = registry.add<TProfile2D>(Form("NchVsThreeParCorrVsT0M_Rep%d", i), Form(";%s;%s;%s", tiNch, tiT0M, tiThreeParCorr), kTProfile2D, {{{nBinsNch, minNch, maxNch}, {nBinsAmpFT0, 0., maxAmpFT0}}});

        pNchVsOneParCorrVsV0A[i] = registry.add<TProfile2D>(Form("NchVsOneParCorrVsV0A_Rep%d", i), Form(";%s;%s;%s", tiNch, tiV0A, tiOneParCorr), kTProfile2D, {{{nBinsNch, minNch, maxNch}, {nBinsAmpV0A, 0., maxAmpV0A}}});
        pNchVsTwoParCorrVsV0A[i] = registry.add<TProfile2D>(Form("NchVsTwoParCorrVsV0A_Rep%d", i), Form(";%s;%s;%s", tiNch, tiV0A, tiTwoParCorr), kTProfile2D, {{{nBinsNch, minNch, maxNch}, {nBinsAmpV0A, 0., maxAmpV0A}}});
        pNchVsThreeParCorrVsV0A[i] = registry.add<TProfile2D>(Form("NchVsThreeParCorrVsV0A_Rep%d", i), Form(";%s;%s;%s", tiNch, tiV0A, tiThreeParCorr), kTProfile2D, {{{nBinsNch, minNch, maxNch}, {nBinsAmpV0A, 0., maxAmpV0A}}});
      }
    }

    if (doprocessMCclosure) {
      registry.add("zPos", ";;Entries;", kTH1F, {axisZpos});
      registry.add("T0Ccent", ";;Entries", kTH1F, {axisCent});
      registry.add("EtaVsPhi", ";#eta;#varphi", kTH2F, {{{axisEta}, {100, -0.1 * PI, +2.1 * PI}}});
      registry.add("ZposVsEta", "", kTProfile, {axisZpos});
      registry.add("sigma1Pt", ";;#sigma(p_{T})/p_{T};", kTProfile, {axisPt});
      registry.add("dcaXYvspT", ";DCA_{xy} (cm);;", kTH2F, {{{50, -1., 1.}, {axisPt}}});
      registry.add("RandomNumber", "", kTH1F, {{50, 0., 1.}});
      registry.add("EvtsDivided", ";Event type;Entries;", kTH1F, {{2, -0.5, 1.5}});
      auto hEvtsDiv = registry.get<TH1>(HIST("EvtsDivided"));
      auto* xEvtsDiv = hEvtsDiv->GetXaxis();
      xEvtsDiv->SetBinLabel(1, "MC closure");
      xEvtsDiv->SetBinLabel(2, "Corrections");
      // MC closure
      registry.add("NchGen", Form("MC Closure;%s;Entries", tiNch), kTH1F, {{nBinsNch, minNch, maxNch}});
      registry.add("NchvsOneParCorrGen", Form("MC Closure;%s;%s", tiNch, tiOneParCorr), kTProfile, {{nBinsNch, minNch, maxNch}});
      registry.add("NchvsTwoParCorrGen", Form("MC Closure;%s;%s", tiNch, tiTwoParCorr), kTProfile, {{nBinsNch, minNch, maxNch}});
      registry.add("NchvsThreeParCorrGen", Form("MC Closure;%s;%s", tiNch, tiThreeParCorr), kTProfile, {{nBinsNch, minNch, maxNch}});
      // Corrections
      registry.add("zPosMC", "Filled at MC closure + Corrections;;Entries;", kTH1F, {axisZpos});
      registry.add("hEventCounterMC", "Event counter", kTH1F, {axisEvent});
      registry.add("nRecColvsCent", "", kTH2F, {{6, -0.5, 5.5}, {{axisCent}}});
      registry.add("Pt_all_ch", Form("Corrections;%s;%s", tiNch, tiPt), kTH2F, {{nBinsNch, minNch, maxNch}, {axisPt}});
      registry.add("Pt_ch", Form("Corrections;%s;%s", tiNch, tiPt), kTH2F, {{nBinsNch, minNch, maxNch}, {axisPt}});
      registry.add("Pt_pi", Form("Corrections;%s;%s", tiNch, tiPt), kTH2F, {{nBinsNch, minNch, maxNch}, {axisPt}});
      registry.add("Pt_ka", Form("Corrections;%s;%s", tiNch, tiPt), kTH2F, {{nBinsNch, minNch, maxNch}, {axisPt}});
      registry.add("Pt_pr", Form("Corrections;%s;%s", tiNch, tiPt), kTH2F, {{nBinsNch, minNch, maxNch}, {axisPt}});
      registry.add("Pt_sigpos", Form("Corrections;%s;%s", tiNch, tiPt), kTH2F, {{nBinsNch, minNch, maxNch}, {axisPt}});
      registry.add("Pt_signeg", Form("Corrections;%s;%s", tiNch, tiPt), kTH2F, {{nBinsNch, minNch, maxNch}, {axisPt}});
      registry.add("Pt_re", Form("Corrections;%s;%s", tiNch, tiPt), kTH2F, {{nBinsNch, minNch, maxNch}, {axisPt}});
      registry.add("PtMC_ch", Form("Corrections;%s;%s", tiNch, tiPt), kTH2F, {{nBinsNch, minNch, maxNch}, {axisPt}});
      registry.add("PtMC_pi", Form("Corrections;%s;%s", tiNch, tiPt), kTH2F, {{nBinsNch, minNch, maxNch}, {axisPt}});
      registry.add("PtMC_ka", Form("Corrections;%s;%s", tiNch, tiPt), kTH2F, {{nBinsNch, minNch, maxNch}, {axisPt}});
      registry.add("PtMC_pr", Form("Corrections;%s;%s", tiNch, tiPt), kTH2F, {{nBinsNch, minNch, maxNch}, {axisPt}});
      registry.add("PtMC_sigpos", Form("Corrections;%s;%s", tiNch, tiPt), kTH2F, {{nBinsNch, minNch, maxNch}, {axisPt}});
      registry.add("PtMC_signeg", Form("Corrections;%s;%s", tiNch, tiPt), kTH2F, {{nBinsNch, minNch, maxNch}, {axisPt}});
      registry.add("PtMC_re", Form("Corrections;%s;%s", tiNch, tiPt), kTH2F, {{nBinsNch, minNch, maxNch}, {axisPt}});
      registry.add("McNchVsFT0M", Form("Corrections;%s;%s", tiT0M, tiNch), kTH2F, {{{nBinsAmpFT0, 0., maxAmpFT0}, {nBinsNch, minNch, maxNch}}});

      auto hECMC = registry.get<TH1>(HIST("hEventCounterMC"));
      auto* x = hECMC->GetXaxis();
      x->SetBinLabel(1, "All");
      x->SetBinLabel(13, "VtxZ cut");

      for (int i = 0; i < kSizeBootStrapEnsemble; i++) {

        hPoissonMC[i] = registry.add<TH1>(Form("PoissonMC_Rep%d", i), ";#it{k};Entries", kTH1F, {{11, -0.5, 10.5}});
        hNchGen[i] = registry.add<TH1>(Form("NchGen_Rep%d", i), Form(";%s;Entries", tiNch), kTH1F, {{nBinsNch, minNch, maxNch}});
        pOneParCorrVsNchGen[i] = registry.add<TProfile>(Form("OneParCorrVsNchGen_Rep%d", i), Form(";%s;%s;", tiNch, tiOneParCorr), kTProfile, {{nBinsNch, minNch, maxNch}});
        pTwoParCorrVsNchGen[i] = registry.add<TProfile>(Form("TwoParCorrVsNchGen_Rep%d", i), Form(";%s;%s;", tiNch, tiTwoParCorr), kTProfile, {{nBinsNch, minNch, maxNch}});
        pThreeParCorrVsNchGen[i] = registry.add<TProfile>(Form("ThreeParCorrVsNchGen_Rep%d", i), Form(";%s;%s;", tiNch, tiThreeParCorr), kTProfile, {{nBinsNch, minNch, maxNch}});

        hNch[i] = registry.add<TH1>(Form("Nch_Rep%d", i), Form(";%s;Entries", tiNch), kTH1F, {{nBinsNch, minNch, maxNch}});
        hPoisson[i] = registry.add<TH1>(Form("Poisson_Rep%d", i), ";#it{k};Entries", kTH1F, {{11, -0.5, 10.5}});

        pOneParCorrVsNch[i] = registry.add<TProfile>(Form("OneParCorrVsNch_Rep%d", i), Form(";%s;%s;", tiNch, tiOneParCorr), kTProfile, {{nBinsNch, minNch, maxNch}});
        pTwoParCorrVsNch[i] = registry.add<TProfile>(Form("TwoParCorrVsNch_Rep%d", i), Form(";%s;%s;", tiNch, tiTwoParCorr), kTProfile, {{nBinsNch, minNch, maxNch}});
        pThreeParCorrVsNch[i] = registry.add<TProfile>(Form("ThreeParCorrVsNch_Rep%d", i), Form(";%s;%s;", tiNch, tiTwoParCorr), kTProfile, {{nBinsNch, minNch, maxNch}});
      }
    }
    if (doprocessQA) {
      registry.add("zPos", ";;Entries;", kTH1F, {axisZpos});
      registry.add("T0Ccent", ";;Entries", kTH1F, {axisCent});
      registry.add("EtaVsPhi", ";#eta;#varphi", kTH2F, {{{axisEta}, {100, -0.1 * PI, +2.1 * PI}}});
      registry.add("ZposVsEta", "", kTProfile, {axisZpos});
      registry.add("sigma1Pt", ";;#sigma(p_{T})/p_{T};", kTProfile, {axisPt});
      registry.add("dcaXYvspT", ";DCA_{xy} (cm);;", kTH2F, {{{50, -1., 1.}, {axisPt}}});
      registry.add("Debunch", ";t_{ZDC}-t_{ZDA};t_{ZDC}+t_{ZDA}", kTH2F, {{{nBinsTDC, minTdc, maxTdc}, {nBinsTDC, minTdc, maxTdc}}});
      registry.add("NchVsFT0M", Form(";%s;%s;", tiT0M, tiNch), kTH2F, {{{nBinsAmpFT0, 0., maxAmpFT0}, {nBinsNch, minNch, maxNch}}});
      registry.add("NchVsFT0A", Form(";%s;%s;", tiT0A, tiNch), kTH2F, {{{nBinsAmpFT0A, 0., maxAmpFT0A}, {nBinsNch, minNch, maxNch}}});
      registry.add("NchVsFT0C", Form(";%s;%s;", tiT0C, tiNch), kTH2F, {{{nBinsAmpFT0C, 0., maxAmpFT0C}, {nBinsNch, minNch, maxNch}}});
      registry.add("NchVsFV0A", Form(";%s;%s;", tiV0A, tiNch), kTH2F, {{{nBinsAmpV0A, 0., maxAmpV0A}, {nBinsNch, minNch, maxNch}}});
      registry.add("NchVsNPV", ";#it{N}_{PV} (|#eta|<1);ITS+TPC tracks (|#eta|<0.8);", kTH2F, {{{nBinsITSTrack, 0, maxITSTrack}, {nBinsNch, minNch, maxNch}}});
      registry.add("NchVsITStracks", ";ITS trks (|#eta|<0.8);TITS+TPC tracks (|#eta|<0.8);", kTH2F, {{{nBinsITSTrack, 0, maxITSTrack}, {nBinsNch, minNch, maxNch}}});
      registry.add("ZNVsFT0A", Form(";%s;%s;", tiT0A, tiZNs), kTH2F, {{{nBinsAmpFT0A, 0., maxAmpFT0A}, {nBinsZN, minZN, maxZN}}});
      registry.add("ZNVsFT0C", Form(";%s;%s;", tiT0C, tiZNs), kTH2F, {{{nBinsAmpFT0C, 0., maxAmpFT0C}, {nBinsZN, minZN, maxZN}}});
      registry.add("ZNVsFT0M", Form(";%s;%s;", tiT0M, tiZNs), kTH2F, {{{nBinsAmpFT0, 0., maxAmpFT0}, {nBinsZN, minZN, maxZN}}});
      registry.add("ZNAamp", ";ZNA;Entries;", kTH1F, {{nBinsZN, minZN, maxZN}});
      registry.add("ZPAamp", ";ZPA;Entries;", kTH1F, {{nBinsZP, minZN, maxZP}});
      registry.add("ZNCamp", ";ZNC;Entries;", kTH1F, {{nBinsZN, minZN, maxZN}});
      registry.add("ZPCamp", ";ZPC;Entries;", kTH1F, {{nBinsZP, minZN, maxZP}});
      registry.add("ZNAVsZNC", ";ZNC;ZNA", kTH2F, {{{nBinsZN, minZN, maxZN}, {nBinsZN, minZN, maxZN}}});
      registry.add("ZPAVsZPC", ";ZPC;ZPA;", kTH2F, {{{nBinsZP, minZN, maxZP}, {nBinsZP, minZN, maxZP}}});
      registry.add("ZNAVsZPA", ";ZPA;ZNA;", kTH2F, {{{nBinsZP, minZN, maxZP}, {nBinsZN, minZN, maxZN}}});
      registry.add("ZNCVsZPC", ";ZPC;ZNC;", kTH2F, {{{nBinsZP, minZN, maxZP}, {nBinsZN, minZN, maxZN}}});
      registry.add("ZNVsZEM", ";ZEM;ZNA+ZNC;", kTH2F, {{{nBinsZN, minZN, maxZEM}, {nBinsZN, minZN, maxZN}}});
      registry.add("ZNCVsNch", Form(";%s;ZNC;", tiNch), kTH2F, {{{nBinsNch, minNch, maxNch}, {nBinsZN, minZN, maxZN}}});
      registry.add("ZNAVsNch", Form(";%s;ZNA;", tiNch), kTH2F, {{{nBinsNch, minNch, maxNch}, {nBinsZN, minZN, maxZN}}});
      registry.add("ZNVsNch", Form(";%s;%s;", tiNch, tiZNs), kTH2F, {{{nBinsNch, minNch, maxNch}, {nBinsZN, minZN, maxZN}}});
      registry.add("ZNDifVsNch", Form(";%s;ZNA-ZNC;", tiNch), kTH2F, {{{nBinsNch, minNch, maxNch}, {100, -50., 50.}}});
    }

    LOG(info) << "\tccdbNoLaterThan=" << ccdbNoLaterThan.value;
    LOG(info) << "\tapplyEff=" << applyEff.value;
    LOG(info) << "\tapplyFD=" << applyFD.value;
    LOG(info) << "\tcorrectNch=" << correctNch.value;
    LOG(info) << "\tpaTHEff=" << paTHEff.value;
    LOG(info) << "\tpaTHFD=" << paTHFD.value;
    LOG(info) << "\tdetectorZDC=" << detectorZDC.value;
    LOG(info) << "\tuseMidRapNchSel=" << useMidRapNchSel.value;
    LOG(info) << "\tdetector4Calibration=" << detector4Calibration.value;
    LOG(info) << "\tnSigmaNchCut=" << nSigmaNchCut.value;
    LOG(info) << "\tpaTHmeanNch=" << paTHmeanNch.value;
    LOG(info) << "\tpaTHsigmaNch=" << paTHsigmaNch.value;
    LOG(info) << "\tminPt=" << minPt.value;
    LOG(info) << "\tmaxPt=" << maxPt.value;
    LOG(info) << "\tmaxPtSpectra=" << maxPtSpectra.value;
    LOG(info) << "\tcurrentRunNumber= " << currentRunNumber;

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    ccdb->setCreatedNotAfter(ccdbNoLaterThan.value);
  }

  template <typename CheckCol>
  bool isEventSelected(CheckCol const& col)
  {
    registry.fill(HIST("hEventCounter"), EvCutLabel::All);
    if (!col.sel8()) {
      return false;
    }
    registry.fill(HIST("hEventCounter"), EvCutLabel::SelEigth);

    if (!col.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    registry.fill(HIST("hEventCounter"), EvCutLabel::NoSameBunchPileup);

    if (!col.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    registry.fill(HIST("hEventCounter"), EvCutLabel::IsGoodZvtxFT0vsPV);

    if (isNoCollInTimeRangeStrict) {
      if (!col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict)) {
        return false;
      }
      registry.fill(HIST("hEventCounter"), EvCutLabel::NoCollInTimeRangeStrict);
    }

    if (isNoCollInTimeRangeStandard) {
      if (!col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
        return false;
      }
      registry.fill(HIST("hEventCounter"), EvCutLabel::NoCollInTimeRangeStandard);
    }

    if (isNoCollInRofStrict) {
      if (!col.selection_bit(o2::aod::evsel::kNoCollInRofStrict)) {
        return false;
      }
      registry.fill(HIST("hEventCounter"), EvCutLabel::NoCollInRofStrict);
    }

    if (isNoCollInRofStandard) {
      if (!col.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
        return false;
      }
      registry.fill(HIST("hEventCounter"), EvCutLabel::NoCollInRofStandard);
    }

    if (isNoHighMultCollInPrevRof) {
      if (!col.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
        return false;
      }
      registry.fill(HIST("hEventCounter"), EvCutLabel::NoHighMultCollInPrevRof);
    }

    // To be used in combination with FT0C-based occupancy
    if (isNoCollInTimeRangeNarrow) {
      if (!col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow)) {
        return false;
      }
      registry.fill(HIST("hEventCounter"), EvCutLabel::NoCollInTimeRangeNarrow);
    }

    if (isOccupancyCut) {
      auto occuValue{isApplyFT0CbasedOccupancy ? col.ft0cOccupancyInTimeRange() : col.trackOccupancyInTimeRange()};
      if (occuValue < minOccCut || occuValue > maxOccCut) {
        return false;
      }
      registry.fill(HIST("hEventCounter"), EvCutLabel::OccuCut);
    }

    if (col.centFT0C() < minT0CcentCut || col.centFT0C() > maxT0CcentCut) {
      return false;
    }
    registry.fill(HIST("hEventCounter"), EvCutLabel::Centrality);

    // Z-vertex position cut
    if (std::fabs(col.posZ()) > posZcut) {
      return false;
    }
    registry.fill(HIST("hEventCounter"), EvCutLabel::VtxZ);

    return true;
  }
  void processQA(o2::aod::ColEvSels::iterator const& collision, o2::aod::BCsRun3 const& /**/, aod::Zdcs const& /**/, aod::FV0As const& /**/, aod::FT0s const& /**/, TheFilteredTracks const& tracks)
  {
    // LOG(info) << " Collisions size: " << collisions.size() << " Table's size: " << collisions.tableSize() << "\n";
    const auto& foundBC = collision.foundBC_as<o2::aod::BCsRun3>();
    // LOG(info) << "Run number: " << foundBC.runNumber() << "\n";
    if (!isEventSelected(collision)) {
      return;
    }

    // has ZDC?
    if (!foundBC.has_zdc()) {
      return;
    }
    registry.fill(HIST("hEventCounter"), EvCutLabel::Zdc);
    auto zdc = foundBC.zdc();

    double aT0A{0.};
    double aT0C{0.};
    if (foundBC.has_ft0()) {
      for (const auto& amplitude : foundBC.ft0().amplitudeA()) {
        aT0A += amplitude;
      }
      for (const auto& amplitude : foundBC.ft0().amplitudeC()) {
        aT0C += amplitude;
      }
    } else {
      return;
    }
    registry.fill(HIST("hEventCounter"), EvCutLabel::TZero);

    double aV0A{0.};
    if (foundBC.has_fv0a()) {
      for (const auto& amplitude : foundBC.fv0a().amplitude()) {
        aV0A += amplitude;
      }
    }

    const double nPV{collision.multNTracksPVeta1() / 1.};
    const double normT0M{(aT0A + aT0C) / arbScale};
    const double normV0A{aV0A / arbScale};
    const double normT0A{aT0A / arbScale};
    const double normT0C{aT0C / arbScale};
    float znA{zdc.amplitudeZNA()};
    float znC{zdc.amplitudeZNC()};
    float zpA{zdc.amplitudeZPA()};
    float zpC{zdc.amplitudeZPC()};
    float aZEM1{zdc.amplitudeZEM1()};
    float aZEM2{zdc.amplitudeZEM2()};
    float tZNA{zdc.timeZNA()};
    float tZNC{zdc.timeZNC()};
    float tZPA{zdc.timeZPA()};
    float tZPC{zdc.timeZPC()};
    float tZDCdif{tZNC + tZPC - tZNA - tZPA};
    float tZDCsum{tZNC + tZPC + tZNA + tZPA};
    znA /= kCollEnergy;
    znC /= kCollEnergy;
    zpA /= kCollEnergy;
    zpC /= kCollEnergy;
    float sumZNs{znA + znC};
    float sumZEMs{aZEM1 + aZEM2};
    // TDC cut
    if (isTDCcut) {
      if (std::sqrt(std::pow(tZDCdif, 2.) + std::pow(tZDCsum, 2.)) > tdcCut) {
        return;
      }
      registry.fill(HIST("hEventCounter"), EvCutLabel::Tdc);
    }

    int itsTracks = 0, glbTracks = 0;
    for (const auto& track : tracks) {
      if (track.hasITS() && ((track.eta() > minEta) && (track.eta() < maxEta))) {
        itsTracks++;
      }
      // Track Selection
      if (!track.isGlobalTrack()) {
        continue;
      }
      if ((track.pt() < minPt) || (track.pt() > maxPt)) {
        continue;
      }
      if ((track.eta() < minEta) || (track.eta() > maxEta)) {
        continue;
      }
      glbTracks++;
    }

    bool skipEvent{false};
    if (useMidRapNchSel) {

      const int nextRunNumber{foundBC.runNumber()};
      if (currentRunNumber != nextRunNumber) {
        loadNchCalibrations(foundBC.timestamp());
        currentRunNumber = nextRunNumber;
        LOG(info) << "\tcurrentRunNumber= " << currentRunNumber << " timeStamp = " << foundBC.timestamp();
      }

      if (!(cfgNch.hMeanNch && cfgNch.hSigmaNch))
        return;

      TString s1 = TString(detector4Calibration.value);
      double xEval{1.};
      if (s1 == "T0M") {
        xEval = normT0M;
      }
      if (s1 == "V0A") {
        xEval = normV0A;
      }
      if (s1 == "NPV") {
        xEval = nPV;
      }

      const int bin4Calibration{cfgNch.hMeanNch->FindBin(xEval)};
      const double meanNch{cfgNch.hMeanNch->GetBinContent(bin4Calibration)};
      const double sigmaNch{cfgNch.hSigmaNch->GetBinContent(bin4Calibration)};
      const double nSigmaSelection{nSigmaNchCut * sigmaNch};
      const double diffMeanNch{meanNch - glbTracks};

      if (std::abs(diffMeanNch) > nSigmaSelection) {
        registry.fill(HIST("ExcludedEvtVsFT0M"), normT0M);
        registry.fill(HIST("ExcludedEvtVsFV0A"), normV0A);
        registry.fill(HIST("ExcludedEvtVsNch"), glbTracks);
        registry.fill(HIST("ExcludedEvtVsNPV"), nPV);
        skipEvent = true;
      }
    }

    if (useMidRapNchSel && skipEvent) {
      return;
    }

    double sumpt{0.};
    for (const auto& track : tracks) {
      // Track Selection
      if (!track.isGlobalTrack()) {
        continue;
      }
      if ((track.pt() < minPt) || (track.pt() > maxPtSpectra)) {
        continue;
      }
      if ((track.eta() < minEta) || (track.eta() > maxEta)) {
        continue;
      }

      registry.fill(HIST("ZposVsEta"), collision.posZ(), track.eta());
      registry.fill(HIST("EtaVsPhi"), track.eta(), track.phi());
      registry.fill(HIST("sigma1Pt"), track.pt(), track.sigma1Pt());
      registry.fill(HIST("dcaXYvspT"), track.dcaXY(), track.pt());
      sumpt += track.pt();
    }

    registry.fill(HIST("zPos"), collision.posZ());
    registry.fill(HIST("T0Ccent"), collision.centFT0C());
    registry.fill(HIST("Debunch"), tZDCdif, tZDCsum);
    registry.fill(HIST("NchVsFV0A"), normV0A, glbTracks);
    registry.fill(HIST("NchVsFT0A"), normT0A, glbTracks);
    registry.fill(HIST("NchVsFT0C"), normT0C, glbTracks);
    registry.fill(HIST("NchVsFT0M"), normT0M, glbTracks);
    registry.fill(HIST("NchUncorrected"), glbTracks);
    registry.fill(HIST("Nch"), glbTracks);
    registry.fill(HIST("NchVsNPV"), collision.multNTracksPVeta1(), glbTracks);
    registry.fill(HIST("NchVsITStracks"), itsTracks, glbTracks);
    if (glbTracks >= minNchSel)
      registry.fill(HIST("NchVsOneParCorr"), glbTracks, sumpt / glbTracks);

    // ZEM cut
    if (sumZEMs > zemCut) {
      registry.fill(HIST("hEventCounter"), EvCutLabel::Zem);
      registry.fill(HIST("ZNAamp"), znA);
      registry.fill(HIST("ZNCamp"), znC);
      registry.fill(HIST("ZPAamp"), zpA);
      registry.fill(HIST("ZPCamp"), zpC);
      registry.fill(HIST("ZNAVsZNC"), znC, znA);
      registry.fill(HIST("ZNAVsZPA"), zpA, znA);
      registry.fill(HIST("ZNCVsZPC"), zpC, znC);
      registry.fill(HIST("ZPAVsZPC"), zpC, zpA);
      registry.fill(HIST("ZNVsZEM"), sumZEMs, sumZNs);
      registry.fill(HIST("ZNVsFT0A"), normT0A, sumZNs);
      registry.fill(HIST("ZNVsFT0C"), normT0C, sumZNs);
      registry.fill(HIST("ZNVsFT0M"), normT0M, sumZNs);
      registry.fill(HIST("ZNAVsNch"), glbTracks, znA);
      registry.fill(HIST("ZNCVsNch"), glbTracks, znC);
      registry.fill(HIST("ZNVsNch"), glbTracks, sumZNs);
      registry.fill(HIST("ZNDifVsNch"), glbTracks, znA - znC);
    }
  }
  PROCESS_SWITCH(UccZdc, processQA, "Process QA", true);
  void processZdcCollAss(o2::aod::ColEvSels::iterator const& collision, o2::aod::BCsRun3 const& /*bcs*/, aod::Zdcs const& /*zdcs*/, aod::FV0As const& /*fv0as*/, aod::FT0s const& /*ft0s*/, TheFilteredTracks const& tracks)
  {
    if (!isEventSelected(collision)) {
      return;
    }

    const auto& foundBC = collision.foundBC_as<o2::aod::BCsRun3>();
    // LOGF(info, "Getting object %s for run number %i from timestamp=%llu", paTH.value.data(), foundBC.runNumber(), foundBC.timestamp());

    // has ZDC?
    if (!foundBC.has_zdc()) {
      return;
    }
    registry.fill(HIST("hEventCounter"), EvCutLabel::Zdc);

    double aT0A{0.};
    double aT0C{0.};
    if (foundBC.has_ft0()) {
      for (const auto& amplitude : foundBC.ft0().amplitudeA()) {
        aT0A += amplitude;
      }
      for (const auto& amplitude : foundBC.ft0().amplitudeC()) {
        aT0C += amplitude;
      }
    } else {
      return;
    }
    registry.fill(HIST("hEventCounter"), EvCutLabel::TZero);

    double aV0A{0.};
    if (foundBC.has_fv0a()) {
      for (const auto& amplitude : foundBC.fv0a().amplitude()) {
        aV0A += amplitude;
      }
    }

    const double nPV{collision.multNTracksPVeta1() / 1.};
    const double normT0M{(aT0A + aT0C) / arbScale};
    const double normV0A{aV0A / arbScale};
    float znA{foundBC.zdc().amplitudeZNA()};
    float znC{foundBC.zdc().amplitudeZNC()};
    float zpA{foundBC.zdc().amplitudeZPA()};
    float zpC{foundBC.zdc().amplitudeZPC()};
    float aZEM1{foundBC.zdc().amplitudeZEM1()};
    float aZEM2{foundBC.zdc().amplitudeZEM2()};
    float tZNA{foundBC.zdc().timeZNA()};
    float tZNC{foundBC.zdc().timeZNC()};
    float tZPA{foundBC.zdc().timeZPA()};
    float tZPC{foundBC.zdc().timeZPC()};
    float tZDCdif{tZNC + tZPC - tZNA - tZPA};
    float tZDCsum{tZNC + tZPC + tZNA + tZPA};
    znA /= kCollEnergy;
    znC /= kCollEnergy;
    zpA /= kCollEnergy;
    zpC /= kCollEnergy;
    double sumZNs{-999.};
    const double sumZPs{zpA + zpC};
    const double sumZEMs{aZEM1 + aZEM2};

    TString sZDC = TString(detectorZDC.value);
    if (sZDC == "ZNA") {
      sumZNs = znA;
    } else if (sZDC == "ZNC") {
      sumZNs = znC;
    } else {
      sumZNs = (znA + znC);
    }

    // TDC cut
    if (isTDCcut) {
      if (std::sqrt(std::pow(tZDCdif, 2.) + std::pow(tZDCsum, 2.)) > tdcCut) {
        return;
      }
      registry.fill(HIST("hEventCounter"), EvCutLabel::Tdc);
    }

    // Nch-based selection
    double glbTracks{0.0};
    for (const auto& track : tracks) {
      // Track Selection
      if (!track.isGlobalTrack()) {
        continue;
      }
      if ((track.pt() < minPt) || (track.pt() > maxPt)) {
        continue;
      }
      if ((track.eta() < minEta) || (track.eta() > maxEta)) {
        continue;
      }
      glbTracks += 1.0;
    }

    bool skipEvent{false};
    if (useMidRapNchSel) {

      const int nextRunNumber{foundBC.runNumber()};
      if (currentRunNumber != nextRunNumber) {
        loadNchCalibrations(foundBC.timestamp());
        currentRunNumber = nextRunNumber;
        LOG(info) << "\tcurrentRunNumber= " << currentRunNumber << " timeStamp = " << foundBC.timestamp();
      }

      if (!(cfgNch.hMeanNch && cfgNch.hSigmaNch))
        return;

      TString s1 = TString(detector4Calibration.value);
      double xEval{1.};
      if (s1 == "T0M") {
        xEval = normT0M;
      }
      if (s1 == "V0A") {
        xEval = normV0A;
      }
      if (s1 == "NPV") {
        xEval = nPV;
      }

      const int bin4Calibration{cfgNch.hMeanNch->FindBin(xEval)};
      const double meanNch{cfgNch.hMeanNch->GetBinContent(bin4Calibration)};
      const double sigmaNch{cfgNch.hSigmaNch->GetBinContent(bin4Calibration)};
      const double nSigmaSelection{nSigmaNchCut * sigmaNch};
      const double diffMeanNch{meanNch - glbTracks};

      if (std::abs(diffMeanNch) > nSigmaSelection) {
        registry.fill(HIST("ExcludedEvtVsFT0M"), normT0M);
        registry.fill(HIST("ExcludedEvtVsFV0A"), normV0A);
        registry.fill(HIST("ExcludedEvtVsNch"), glbTracks);
        registry.fill(HIST("ExcludedEvtVsNPV"), nPV);
        skipEvent = true;
      }
    }

    // Skip event based on number of Nch sigmas
    if (useMidRapNchSel && skipEvent) {
      return;
    }

    // Reject low-multiplcicity events
    if (glbTracks < minNchSel) {
      return;
    }

    double nchMult{glbTracks};
    std::vector<double> pTs;
    std::vector<double> vecFD;
    std::vector<double> vecEff;

    // apply corrections
    if (applyEff || applyFD) {
      nchMult = 0.;
      loadCorrections(foundBC.timestamp());
      if (!(cfg.hEfficiency && cfg.hFeedDown))
        return;

      const int foundNchBin{cfg.hEfficiency->GetXaxis()->FindBin(glbTracks)};

      // Calculates the Corrected Nch
      for (const auto& track : tracks) {
        // Track Selection
        if (!track.isGlobalTrack()) {
          continue;
        }
        if ((track.pt() < minPt) || (track.pt() > maxPt)) {
          continue;
        }
        if ((track.eta() < minEta) || (track.eta() > maxEta)) {
          continue;
        }

        const float pt{track.pt()};
        const int foundPtBin{cfg.hEfficiency->GetYaxis()->FindBin(pt)};
        const double effValue{cfg.hEfficiency->GetBinContent(foundNchBin, foundPtBin)};
        double fdValue{1.};

        if (applyFD)
          fdValue = cfg.hFeedDown->GetBinContent(foundNchBin, foundPtBin);

        if ((effValue > 0.) && (fdValue > 0.)) {
          nchMult += (std::pow(effValue, -1.) * fdValue);
        }
      }

      if (applyEff && !correctNch) {
        nchMult = glbTracks;
      }

      // Reject low-multiplcicity events
      if (nchMult < minNchSel) {
        return;
      }

      // Fill vectors for [pT] measurement
      for (const auto& track : tracks) {
        // Track Selection
        if (!track.isGlobalTrack()) {
          continue;
        }
        if ((track.pt() < minPt) || (track.pt() > maxPtSpectra)) {
          continue;
        }
        if ((track.eta() < minEta) || (track.eta() > maxEta)) {
          continue;
        }

        const float pt{track.pt()};
        const int foundPtBin{cfg.hEfficiency->GetYaxis()->FindBin(pt)};
        const double effValue{cfg.hEfficiency->GetBinContent(foundNchBin, foundPtBin)};
        double fdValue{1.};

        if (applyFD)
          fdValue = cfg.hFeedDown->GetBinContent(foundNchBin, foundPtBin);

        if ((effValue > 0.) && (fdValue > 0.)) {
          pTs.emplace_back(pt);
          vecEff.emplace_back(effValue);
          vecFD.emplace_back(fdValue);
          // To calculate event-averaged <pt>
          registry.fill(HIST("NchVsZNVsPt"), nchMult, sumZNs, pt * (fdValue / effValue));
        }
      }
    } else {
      // Fill vectors for [pT] measurement
      for (const auto& track : tracks) {
        // Track Selection
        if (!track.isGlobalTrack()) {
          continue;
        }
        if ((track.pt() < minPt) || (track.pt() > maxPtSpectra)) {
          continue;
        }
        if ((track.eta() < minEta) || (track.eta() > maxEta)) {
          continue;
        }

        pTs.emplace_back(track.pt());
        vecEff.emplace_back(1.);
        vecFD.emplace_back(1.);
        // To calculate event-averaged <pt>
        registry.fill(HIST("NchVsZNVsPt"), nchMult, sumZNs, track.pt());
      }
    }

    double p1, p2, p3, p4, w1, w2, w3, w4;
    p1 = p2 = p3 = p4 = w1 = w2 = w3 = w4 = 0.0;
    getPTpowers(pTs, vecEff, vecFD, p1, w1, p2, w2, p3, w3, p4, w4);

    // EbE one-particle pT correlation
    const double oneParCorr{p1 / w1};

    // EbE two-particle pT correlation
    const double denTwoParCorr{std::pow(w1, 2.) - w2};
    const double numTwoParCorr{std::pow(p1, 2.) - p2};
    const double twoParCorr{numTwoParCorr / denTwoParCorr};

    // EbE three-particle pT correlation
    const double denThreeParCorr{std::pow(w1, 3.) - (3. * w2 * w1) + (2. * w3)};
    const double numThreeParCorr{std::pow(p1, 3.) - (3. * p2 * p1) + (2. * p3)};
    const double threeParCorr{numThreeParCorr / denThreeParCorr};

    // One-dimensional distributions
    registry.fill(HIST("Nch"), nchMult);
    registry.fill(HIST("NchUncorrected"), glbTracks);
    registry.fill(HIST("NchVsV0A"), nchMult, normV0A);
    registry.fill(HIST("NchVsT0M"), nchMult, normT0M);

    registry.fill(HIST("NchVsOneParCorr"), nchMult, oneParCorr);
    registry.fill(HIST("NchVsTwoParCorr"), nchMult, twoParCorr);
    registry.fill(HIST("NchVsThreeParCorr"), nchMult, threeParCorr);

    registry.fill(HIST("NchVsOneParCorrVsT0M"), nchMult, normT0M, oneParCorr);
    registry.fill(HIST("NchVsTwoParCorrVsT0M"), nchMult, normT0M, twoParCorr);
    registry.fill(HIST("NchVsThreeParCorrVsT0M"), nchMult, normT0M, threeParCorr);

    registry.fill(HIST("NchVsOneParCorrVsV0A"), nchMult, normV0A, oneParCorr);
    registry.fill(HIST("NchVsTwoParCorrVsV0A"), nchMult, normV0A, twoParCorr);
    registry.fill(HIST("NchVsThreeParCorrVsV0A"), nchMult, normV0A, threeParCorr);

    if (sumZEMs > zemCut) {
      registry.fill(HIST("hEventCounter"), EvCutLabel::Zem);
      registry.fill(HIST("NchVsZN"), nchMult, sumZNs);
      registry.fill(HIST("NchVsZP"), nchMult, sumZPs);
      registry.fill(HIST("NchVsOneParCorrVsZN"), nchMult, sumZNs, oneParCorr);
      registry.fill(HIST("NchVsTwoParCorrVsZN"), nchMult, sumZNs, twoParCorr);
      registry.fill(HIST("NchVsThreeParCorrVsZN"), nchMult, sumZNs, threeParCorr);
    }

    const uint64_t timeStamp{foundBC.timestamp()};
    eventSampling(tracks, normV0A, normT0M, sumZNs, sumZEMs, timeStamp);
  }
  PROCESS_SWITCH(UccZdc, processZdcCollAss, "Process ZDC W/Coll Ass.", true);

  Preslice<TheFilteredSimTracks> perCollision = aod::track::collisionId;
  Service<o2::framework::O2DatabasePDG> pdg;
  void processMCclosure(aod::McCollisions::iterator const& mccollision, soa::SmallGroups<o2::aod::SimCollisions> const& collisions, o2::aod::BCsRun3 const& /*bcs*/, aod::FT0s const& /*ft0s*/, aod::FV0As const& /*fv0as*/, aod::McParticles const& mcParticles, TheFilteredSimTracks const& simTracks)
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

      registry.fill(HIST("hEventCounterMC"), EvCutLabel::All);
      // Vtx_z selection MC
      if (std::fabs(mccollision.posZ()) > posZcut) {
        continue;
      }

      const auto& foundBC = collision.foundBC_as<o2::aod::BCsRun3>();

      uint64_t timeStamp{foundBC.timestamp()};
      TRandom3 rndGen(timeStamp);
      const double rndNum{rndGen.Uniform(0.0, 1.0)};
      registry.fill(HIST("RandomNumber"), rndNum);

      double aT0A = 0., aT0C = 0.;
      if (foundBC.has_ft0()) {
        for (const auto& amplitude : foundBC.ft0().amplitudeA()) {
          aT0A += amplitude;
        }
        for (const auto& amplitude : foundBC.ft0().amplitudeC()) {
          aT0C += amplitude;
        }
      } else {
        continue;
      }

      const double normT0M{(aT0A + aT0C) / arbScale};
      double nchRaw{0.};
      double nchMult{0.};
      double nchMC{0.};

      registry.fill(HIST("zPos"), collision.posZ());
      registry.fill(HIST("zPosMC"), mccollision.posZ());
      registry.fill(HIST("hEventCounterMC"), EvCutLabel::VtxZ);

      if (skipRecoColGTOne && (collisions.size() > kOne)) {
        continue;
      }

      registry.fill(HIST("nRecColvsCent"), collisions.size(), collision.centFT0C());

      const auto& cent{collision.centFT0C()};
      registry.fill(HIST("T0Ccent"), cent);

      const auto& groupedTracks{simTracks.sliceBy(perCollision, collision.globalIndex())};

      // Half of the statistics for MC closure
      if (rndNum >= kZero && rndNum < evtFracMCcl) {

        registry.fill(HIST("EvtsDivided"), 0);

        // Run-by-run efficiency
        loadCorrections(foundBC.timestamp());
        if (!(cfg.hEfficiency && cfg.hFeedDown)) {
          continue;
        }

        std::vector<double> pTs;
        std::vector<double> vecFD;
        std::vector<double> vecEff;

        // Calculates the event's Nch to evaluate the efficiency
        for (const auto& track : groupedTracks) {
          // Track Selection
          if (track.eta() < minEta || track.eta() > maxEta) {
            continue;
          }
          if (track.pt() < minPt || track.pt() > maxPt) {
            continue;
          }
          if (!track.isGlobalTrack()) {
            continue;
          }
          nchRaw++;
        }

        // Reject event if nchRaw less than a lower cutoff
        if (nchRaw < minNchSel) {
          continue;
        }

        const int foundNchBin{cfg.hEfficiency->GetXaxis()->FindBin(nchRaw)};

        // Calculates the event weight, W_k
        for (const auto& track : groupedTracks) {
          // Track Selection
          if (track.eta() < minEta || track.eta() > maxEta) {
            continue;
          }
          if (track.pt() < minPt || track.pt() > maxPtSpectra) {
            continue;
          }
          if (!track.isGlobalTrack()) {
            continue;
          }
          if (!track.has_mcParticle()) {
            continue;
          }
          const auto& particle{track.mcParticle()};

          auto charge{0.};
          // Get the MC particle
          auto* pdgParticle = pdg->GetParticle(particle.pdgCode());
          if (pdgParticle != nullptr) {
            charge = pdgParticle->Charge();
          } else {
            continue;
          }

          // Is it a charged particle?
          if (std::abs(charge) < kMinCharge) {
            continue;
          }

          const double pt{static_cast<double>(track.pt())};
          const int foundPtBin{cfg.hEfficiency->GetYaxis()->FindBin(pt)};
          const double effValue{cfg.hEfficiency->GetBinContent(foundNchBin, foundPtBin)};
          double fdValue{1.};

          if (applyFD)
            fdValue = cfg.hFeedDown->GetBinContent(foundNchBin, foundPtBin);
          if ((effValue > 0.) && (fdValue > 0.)) {
            pTs.emplace_back(pt);
            vecEff.emplace_back(effValue);
            vecFD.emplace_back(fdValue);
            nchMult += (std::pow(effValue, -1.0) * fdValue);
          }
        }

        double p1, p2, p3, p4, w1, w2, w3, w4;
        p1 = p2 = p3 = p4 = w1 = w2 = w3 = w4 = 0.0;
        getPTpowers(pTs, vecEff, vecFD, p1, w1, p2, w2, p3, w3, p4, w4);

        const double denTwoParCorr{std::pow(w1, 2.) - w2};
        const double numTwoParCorr{std::pow(p1, 2.) - p2};
        const double denThreeParCorr{std::pow(w1, 3.) - (3. * w2 * w1) + (2. * w3)};
        const double numThreeParCorr{std::pow(p1, 3.) - (3. * p2 * p1) + (2. * p3)};

        const double oneParCorr{p1 / w1};
        const double twoParCorr{numTwoParCorr / denTwoParCorr};
        const double threeParCorr{numThreeParCorr / denThreeParCorr};

        registry.fill(HIST("Nch"), nchMult);
        registry.fill(HIST("NchUncorrected"), nchRaw);
        registry.fill(HIST("NchVsOneParCorr"), nchMult, oneParCorr);
        registry.fill(HIST("NchVsTwoParCorr"), nchMult, twoParCorr);
        registry.fill(HIST("NchVsThreeParCorr"), nchMult, threeParCorr);

        //--------------------------- Generated MC ---------------------------
        std::vector<float> pTsMC;
        std::vector<float> vecFullEff;
        std::vector<float> vecFDEqualOne;

        // calculates the  true Nch
        for (const auto& particle : mcParticles) {
          if (particle.eta() < minEta || particle.eta() > maxEta) {
            continue;
          }
          if (particle.pt() < minPt || particle.pt() > maxPt) {
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
          if (std::abs(charge) < kMinCharge) {
            continue;
          }
          // Is it a primary particle?
          if (!particle.isPhysicalPrimary()) {
            continue;
          }
          nchMC++;
        }

        if (nchMC < minNchSel) {
          continue;
        }

        // Calculates the event weight, W_k
        for (const auto& particle : mcParticles) {
          if (particle.eta() < minEta || particle.eta() > maxEta) {
            continue;
          }
          if (particle.pt() < minPt || particle.pt() > maxPtSpectra) {
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
          if (std::abs(charge) < kMinCharge) {
            continue;
          }
          // Is it a primary particle?
          if (!particle.isPhysicalPrimary()) {
            continue;
          }

          const float pt{particle.pt()};
          pTsMC.emplace_back(pt);
          vecFullEff.emplace_back(1.);
          vecFDEqualOne.emplace_back(1.);
        }

        double p1MC, p2MC, p3MC, p4MC, w1MC, w2MC, w3MC, w4MC;
        p1MC = p2MC = p3MC = p4MC = w1MC = w2MC = w3MC = w4MC = 0.0;
        getPTpowers(pTsMC, vecFullEff, vecFDEqualOne, p1MC, w1MC, p2MC, w2MC, p3MC, w3MC, p4MC, w4MC);

        const double denTwoParCorrMC{std::pow(w1MC, 2.) - w2MC};
        const double numTwoParCorrMC{std::pow(p1MC, 2.) - p2MC};
        const double denThreeParCorrMC{std::pow(w1MC, 3.) - (3. * w2MC * w1MC) + (2. * w3MC)};
        const double numThreeParCorrMC{std::pow(p1MC, 3.) - (3. * p2MC * p1MC) + (2. * p3MC)};

        const double oneParCorrMC{p1MC / w1MC};
        const double twoParCorrMC{numTwoParCorrMC / denTwoParCorrMC};
        const double threeParCorrMC{numThreeParCorrMC / denThreeParCorrMC};

        registry.fill(HIST("NchGen"), nchMC);
        registry.fill(HIST("NchvsOneParCorrGen"), nchMC, oneParCorrMC);
        registry.fill(HIST("NchvsTwoParCorrGen"), nchMC, twoParCorrMC);
        registry.fill(HIST("NchvsThreeParCorrGen"), nchMC, threeParCorrMC);

        //------------------ Poisson sampling
        eventSamplingMC(mcParticles, timeStamp);
        eventSamplingMCRec(groupedTracks, timeStamp);
      } else { // Correction with the remaining half of the sample
        registry.fill(HIST("EvtsDivided"), 1);
        //----- MC reconstructed -----//
        for (const auto& track : groupedTracks) {
          // Track Selection
          if (track.eta() < minEta || track.eta() > maxEta) {
            continue;
          }
          if (track.pt() < minPt || track.pt() > maxPt) {
            continue;
          }
          if (!track.isGlobalTrack()) {
            continue;
          }
          nchRaw++;
        }

        for (const auto& track : groupedTracks) {
          // Track Selection
          if (track.eta() < minEta || track.eta() > maxEta) {
            continue;
          }
          if (track.pt() < minPt || track.pt() > maxPtSpectra) {
            continue;
          }
          if (!track.isGlobalTrack()) {
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
          if (std::abs(charge) < kMinCharge) {
            continue;
          }

          // All charged particles
          registry.fill(HIST("Pt_all_ch"), nchRaw, track.pt());
          registry.fill(HIST("ZposVsEta"), collision.posZ(), track.eta());
          registry.fill(HIST("EtaVsPhi"), track.eta(), track.phi());
          registry.fill(HIST("dcaXYvspT"), track.dcaXY(), track.pt());

          // Is it a primary particle?
          if (!particle.isPhysicalPrimary()) {
            continue;
          }

          registry.fill(HIST("Pt_ch"), nchRaw, track.pt());
          if (particle.pdgCode() == PDG_t::kPiPlus || particle.pdgCode() == PDG_t::kPiMinus) {
            registry.fill(HIST("Pt_pi"), nchRaw, track.pt());
          } else if (particle.pdgCode() == PDG_t::kKPlus || particle.pdgCode() == PDG_t::kKMinus) {
            registry.fill(HIST("Pt_ka"), nchRaw, track.pt());
          } else if (particle.pdgCode() == PDG_t::kProton || particle.pdgCode() == PDG_t::kProtonBar) {
            registry.fill(HIST("Pt_pr"), nchRaw, track.pt());
          } else if (particle.pdgCode() == PDG_t::kSigmaPlus || particle.pdgCode() == PDG_t::kSigmaBarMinus) {
            registry.fill(HIST("Pt_sigpos"), nchRaw, track.pt());
          } else if (particle.pdgCode() == PDG_t::kSigmaMinus || particle.pdgCode() == PDG_t::kSigmaBarPlus) {
            registry.fill(HIST("Pt_signeg"), nchRaw, track.pt());
          } else {
            registry.fill(HIST("Pt_re"), nchRaw, track.pt());
          }
        }

        // Generated MC
        for (const auto& particle : mcParticles) {
          if (particle.eta() < minEta || particle.eta() > maxEta) {
            continue;
          }
          if (particle.pt() < minPt || particle.pt() > maxPtSpectra) {
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
          if (std::abs(charge) < kMinCharge) {
            continue;
          }
          // Is it a primary particle?
          if (!particle.isPhysicalPrimary()) {
            continue;
          }

          registry.fill(HIST("PtMC_ch"), nchRaw, particle.pt());
          if (particle.pdgCode() == PDG_t::kPiPlus || particle.pdgCode() == PDG_t::kPiMinus) { // pion
            registry.fill(HIST("PtMC_pi"), nchRaw, particle.pt());
          } else if (particle.pdgCode() == PDG_t::kKPlus || particle.pdgCode() == PDG_t::kKMinus) { // kaon
            registry.fill(HIST("PtMC_ka"), nchRaw, particle.pt());
          } else if (particle.pdgCode() == PDG_t::kProton || particle.pdgCode() == PDG_t::kProtonBar) { // proton
            registry.fill(HIST("PtMC_pr"), nchRaw, particle.pt());
          } else if (particle.pdgCode() == PDG_t::kSigmaPlus || particle.pdgCode() == PDG_t::kSigmaBarMinus) { // positive sigma
            registry.fill(HIST("PtMC_sigpos"), nchRaw, particle.pt());
          } else if (particle.pdgCode() == PDG_t::kSigmaMinus || particle.pdgCode() == PDG_t::kSigmaBarPlus) { // negative sigma
            registry.fill(HIST("PtMC_signeg"), nchRaw, particle.pt());
          } else { // rest
            registry.fill(HIST("PtMC_re"), nchRaw, particle.pt());
          }
        }
      } // Half of statistics for corrections
      registry.fill(HIST("McNchVsFT0M"), normT0M, nchRaw);
    } // Collisions
  }
  PROCESS_SWITCH(UccZdc, processMCclosure, "Process MC closure", false);

  template <typename T, typename U>
  void getPTpowers(const T& pTs, const T& vecEff, const T& vecFD, U& pOne, U& wOne, U& pTwo, U& wTwo, U& pThree, U& wThree, U& pFour, U& wFour)
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

  template <typename T, typename V>
  void eventSamplingMC(const T& mcParticles, const V& timeStamp)
  {
    TRandom3 rndGen(timeStamp);
    std::vector<uint64_t> vPoisson;
    for (int replica = 0; replica < kSizeBootStrapEnsemble; ++replica)
      vPoisson.emplace_back(rndGen.Poisson(1.));

    for (int replica = 0; replica < kSizeBootStrapEnsemble; ++replica) {

      hPoissonMC[replica]->Fill(vPoisson.at(replica));

      for (uint64_t evtRep = 0; evtRep < vPoisson.at(replica); ++evtRep) {

        double nchMult{0.};
        std::vector<float> pTs;
        std::vector<float> vecFD;
        std::vector<float> vecEff;

        // Calculates the true Nch
        for (const auto& particle : mcParticles) {
          if (particle.eta() < minEta || particle.eta() > maxEta) {
            continue;
          }
          if (particle.pt() < minPt || particle.pt() > maxPt) {
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
          if (std::abs(charge) < kMinCharge) {
            continue;
          }
          // Is it a primary particle?
          if (!particle.isPhysicalPrimary()) {
            continue;
          }
          nchMult++;
        }

        if (nchMult < minNchSel) {
          continue;
        }

        // Calculates the event weight, W_k
        for (const auto& particle : mcParticles) {
          if (particle.eta() < minEta || particle.eta() > maxEta) {
            continue;
          }
          if (particle.pt() < minPt || particle.pt() > maxPtSpectra) {
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
          if (std::abs(charge) < kMinCharge) {
            continue;
          }
          // Is it a primary particle?
          if (!particle.isPhysicalPrimary()) {
            continue;
          }

          const float pt{particle.pt()};
          pTs.emplace_back(pt);
          vecEff.emplace_back(1.);
          vecFD.emplace_back(1.);
        }

        double p1, p2, p3, p4, w1, w2, w3, w4;
        p1 = p2 = p3 = p4 = w1 = w2 = w3 = w4 = 0.0;
        getPTpowers(pTs, vecEff, vecFD, p1, w1, p2, w2, p3, w3, p4, w4);

        // EbE one-particle pT correlation
        const double oneParCorr{p1 / w1};

        // EbE two-particle pT correlation
        const double denTwoParCorr{std::pow(w1, 2.) - w2};
        const double numTwoParCorr{std::pow(p1, 2.) - p2};
        const double twoParCorr{numTwoParCorr / denTwoParCorr};

        // EbE three-particle pT correlation
        const double denThreeParCorr{std::pow(w1, 3.) - (3. * w2 * w1) + (2. * w3)};
        const double numThreeParCorr{std::pow(p1, 3.) - (3. * p2 * p1) + (2. * p3)};
        const double threeParCorr{numThreeParCorr / denThreeParCorr};

        hNchGen[replica]->Fill(nchMult);
        pOneParCorrVsNchGen[replica]->Fill(nchMult, oneParCorr);
        pTwoParCorrVsNchGen[replica]->Fill(nchMult, twoParCorr);
        pThreeParCorrVsNchGen[replica]->Fill(nchMult, threeParCorr);
      } // event per replica
    } // replica's loop
  }

  template <typename T, typename V>
  void eventSamplingMCRec(const T& tracks, const V& timeStamp)
  {
    TRandom3 rndGen(timeStamp);
    std::vector<uint64_t> vPoisson;
    for (int replica = 0; replica < kSizeBootStrapEnsemble; ++replica)
      vPoisson.emplace_back(rndGen.Poisson(1.));

    for (int replica = 0; replica < kSizeBootStrapEnsemble; ++replica) {

      hPoisson[replica]->Fill(vPoisson.at(replica));

      for (uint64_t evtRep = 0; evtRep < vPoisson.at(replica); ++evtRep) {

        double nchMult{0.};
        int glbTracks{0};
        std::vector<float> pTs;
        std::vector<float> vecFD;
        std::vector<float> vecEff;

        // Calculates the uncorrected Nch
        for (const auto& track : tracks) {
          // Track Selection
          if (!track.isGlobalTrack()) {
            continue;
          }
          if ((track.pt() < minPt) || (track.pt() > maxPt)) {
            continue;
          }
          if ((track.eta() < minEta) || (track.eta() > maxEta)) {
            continue;
          }
          glbTracks++;
        }

        if (glbTracks < minNchSel) {
          continue;
        }

        // Calculates the Nch multiplicity if corrections are loaded
        if (cfg.correctionsLoaded) {
          const int foundNchBin{cfg.hEfficiency->GetXaxis()->FindBin(glbTracks)};
          for (const auto& track : tracks) {
            // Track Selection
            if (!track.isGlobalTrack()) {
              continue;
            }
            if ((track.pt() < minPt) || (track.pt() > maxPt)) {
              continue;
            }
            if ((track.eta() < minEta) || (track.eta() > maxEta)) {
              continue;
            }

            const float pt{track.pt()};
            double fdValue{1.};
            const int foundPtBin{cfg.hEfficiency->GetYaxis()->FindBin(pt)};
            const double effValue{cfg.hEfficiency->GetBinContent(foundNchBin, foundPtBin)};

            if (applyFD)
              fdValue = cfg.hFeedDown->GetBinContent(foundNchBin, foundPtBin);
            if ((effValue > 0.) && (fdValue > 0.)) {
              nchMult += (std::pow(effValue, -1.) * fdValue);
            }
          }
        }

        if (!applyEff)
          nchMult = static_cast<double>(glbTracks);
        if (applyEff && !correctNch)
          nchMult = static_cast<double>(glbTracks);

        // Fill vectors for [pT] measurement
        if (cfg.correctionsLoaded) {
          const int foundNchBin{cfg.hEfficiency->GetXaxis()->FindBin(glbTracks)};
          // Fill vectors for [pT] measurement
          for (const auto& track : tracks) {
            // Track Selection
            if (!track.isGlobalTrack()) {
              continue;
            }
            if ((track.pt() < minPt) || (track.pt() > maxPtSpectra)) {
              continue;
            }
            if ((track.eta() < minEta) || (track.eta() > maxEta)) {
              continue;
            }

            const float pt{track.pt()};
            double fdValue{1.};
            const int foundPtBin{cfg.hEfficiency->GetYaxis()->FindBin(pt)};
            const double effValue{cfg.hEfficiency->GetBinContent(foundNchBin, foundPtBin)};

            if (applyFD)
              fdValue = cfg.hFeedDown->GetBinContent(foundNchBin, foundPtBin);

            if ((effValue > 0.) && (fdValue > 0.)) {
              pTs.emplace_back(pt);
              vecEff.emplace_back(effValue);
              vecFD.emplace_back(fdValue);
            }
          }
        } else {
          for (const auto& track : tracks) {
            // Track Selection
            if (!track.isGlobalTrack()) {
              continue;
            }
            if ((track.pt() < minPt) || (track.pt() > maxPtSpectra)) {
              continue;
            }

            pTs.emplace_back(track.pt());
            vecEff.emplace_back(1.);
            vecFD.emplace_back(1.);
          }
        }

        double p1, p2, p3, p4, w1, w2, w3, w4;
        p1 = p2 = p3 = p4 = w1 = w2 = w3 = w4 = 0.0;
        getPTpowers(pTs, vecEff, vecFD, p1, w1, p2, w2, p3, w3, p4, w4);

        // EbE one-particle pT correlation
        const double oneParCorr{p1 / w1};

        // EbE two-particle pT correlation
        const double denTwoParCorr{std::pow(w1, 2.) - w2};
        const double numTwoParCorr{std::pow(p1, 2.) - p2};
        const double twoParCorr{numTwoParCorr / denTwoParCorr};

        // EbE three-particle pT correlation
        const double denThreeParCorr{std::pow(w1, 3.) - (3. * w2 * w1) + (2. * w3)};
        const double numThreeParCorr{std::pow(p1, 3.) - (3. * p2 * p1) + (2. * p3)};
        const double threeParCorr{numThreeParCorr / denThreeParCorr};

        hNch[replica]->Fill(nchMult);
        pOneParCorrVsNch[replica]->Fill(nchMult, oneParCorr);
        pTwoParCorrVsNch[replica]->Fill(nchMult, twoParCorr);
        pThreeParCorrVsNch[replica]->Fill(nchMult, threeParCorr);
      } // event per replica
    } // replica's loop
  }

  template <typename T, typename U, typename V>
  void eventSampling(const T& tracks, const U& normV0A, const U& normT0M, const U& sumZNs, const U& sumZEMs, const V& timeStamp)
  {
    TRandom3 rndGen(timeStamp);
    std::vector<uint64_t> vPoisson;
    for (int replica = 0; replica < kSizeBootStrapEnsemble; ++replica)
      vPoisson.emplace_back(rndGen.Poisson(1.));

    for (int replica = 0; replica < kSizeBootStrapEnsemble; ++replica) {

      hPoisson[replica]->Fill(vPoisson.at(replica));

      for (uint64_t evtRep = 0; evtRep < vPoisson.at(replica); ++evtRep) {

        double nchMult{0.};
        int glbTracks{0};
        std::vector<float> pTs;
        std::vector<float> vecFD;
        std::vector<float> vecEff;

        // Calculates the uncorrected Nch multiplicity
        for (const auto& track : tracks) {
          // Track Selection
          if (!track.isGlobalTrack()) {
            continue;
          }
          if ((track.pt() < minPt) || (track.pt() > maxPt)) {
            continue;
          }
          if ((track.eta() < minEta) || (track.eta() > maxEta)) {
            continue;
          }
          glbTracks++;
        }

        if (glbTracks < minNchSel) {
          continue;
        }

        // Calculates the Nch multiplicity if corrections are loaded
        if (cfg.correctionsLoaded) {
          const int foundNchBin{cfg.hEfficiency->GetXaxis()->FindBin(glbTracks)};
          for (const auto& track : tracks) {
            // Track Selection
            if (!track.isGlobalTrack()) {
              continue;
            }
            if ((track.pt() < minPt) || (track.pt() > maxPt)) {
              continue;
            }
            if ((track.eta() < minEta) || (track.eta() > maxEta)) {
              continue;
            }

            float pt{track.pt()};
            double fdValue{1.};
            int foundPtBin{cfg.hEfficiency->GetYaxis()->FindBin(pt)};
            double effValue{cfg.hEfficiency->GetBinContent(foundNchBin, foundPtBin)};

            if (applyFD)
              fdValue = cfg.hFeedDown->GetBinContent(foundNchBin, foundPtBin);
            if ((effValue > 0.) && (fdValue > 0.)) {
              nchMult += (std::pow(effValue, -1.) * fdValue);
            }
          }
        }

        if (!applyEff)
          nchMult = static_cast<double>(glbTracks);
        if (applyEff && !correctNch)
          nchMult = static_cast<double>(glbTracks);

        // Fill vectors for [pT] measurement
        if (cfg.correctionsLoaded) {
          const int foundNchBin{cfg.hEfficiency->GetXaxis()->FindBin(glbTracks)};
          // Fill vectors for [pT] measurement
          for (const auto& track : tracks) {
            // Track Selection
            if (!track.isGlobalTrack()) {
              continue;
            }
            if ((track.pt() < minPt) || (track.pt() > maxPtSpectra)) {
              continue;
            }
            if ((track.eta() < minEta) || (track.eta() > maxEta)) {
              continue;
            }

            float pt{track.pt()};
            double fdValue{1.};
            int foundPtBin{cfg.hEfficiency->GetYaxis()->FindBin(pt)};
            double effValue{cfg.hEfficiency->GetBinContent(foundNchBin, foundPtBin)};

            if (applyFD)
              fdValue = cfg.hFeedDown->GetBinContent(foundNchBin, foundPtBin);

            if ((effValue > 0.) && (fdValue > 0.)) {
              pTs.emplace_back(pt);
              vecEff.emplace_back(effValue);
              vecFD.emplace_back(fdValue);
            }
          }
        } else {
          for (const auto& track : tracks) {
            // Track Selection
            if (!track.isGlobalTrack()) {
              continue;
            }
            if ((track.pt() < minPt) || (track.pt() > maxPtSpectra)) {
              continue;
            }

            pTs.emplace_back(track.pt());
            vecEff.emplace_back(1.);
            vecFD.emplace_back(1.);
          }
        }

        double p1, p2, p3, p4, w1, w2, w3, w4;
        p1 = p2 = p3 = p4 = w1 = w2 = w3 = w4 = 0.0;
        getPTpowers(pTs, vecEff, vecFD, p1, w1, p2, w2, p3, w3, p4, w4);

        // EbE one-particle pT correlation
        const double oneParCorr{p1 / w1};

        // EbE two-particle pT correlation
        const double denTwoParCorr{std::pow(w1, 2.) - w2};
        const double numTwoParCorr{std::pow(p1, 2.) - p2};
        const double twoParCorr{numTwoParCorr / denTwoParCorr};

        // EbE three-particle pT correlation
        const double denThreeParCorr{std::pow(w1, 3.) - (3. * w2 * w1) + (2. * w3)};
        const double numThreeParCorr{std::pow(p1, 3.) - (3. * p2 * p1) + (2. * p3)};
        const double threeParCorr{numThreeParCorr / denThreeParCorr};

        hNchVsV0A[replica]->Fill(nchMult, normV0A);
        hNchVsT0M[replica]->Fill(nchMult, normT0M);

        pNchVsOneParCorrVsT0M[replica]->Fill(nchMult, normT0M, oneParCorr);
        pNchVsTwoParCorrVsT0M[replica]->Fill(nchMult, normT0M, twoParCorr);
        pNchVsThreeParCorrVsT0M[replica]->Fill(nchMult, normT0M, threeParCorr);

        pNchVsOneParCorrVsV0A[replica]->Fill(nchMult, normV0A, oneParCorr);
        pNchVsTwoParCorrVsV0A[replica]->Fill(nchMult, normV0A, twoParCorr);
        pNchVsThreeParCorrVsV0A[replica]->Fill(nchMult, normV0A, threeParCorr);

        if (sumZEMs > zemCut) {
          hNchVsZN[replica]->Fill(nchMult, sumZNs);
          pNchVsOneParCorrVsZN[replica]->Fill(nchMult, sumZNs, oneParCorr);
          pNchVsTwoParCorrVsZN[replica]->Fill(nchMult, sumZNs, twoParCorr);
          pNchVsThreeParCorrVsZN[replica]->Fill(nchMult, sumZNs, threeParCorr);
        }
      } // event per replica
    } // replica's loop
  }

  void loadCorrections(uint64_t timeStamp)
  {
    if (paTHEff.value.empty() == false) {
      cfg.hEfficiency = ccdb->getForTimeStamp<TH2F>(paTHEff, timeStamp);
      if (cfg.hEfficiency == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram from %s", paTHEff.value.c_str());
      }
    }

    if (paTHFD.value.empty() == false) {
      cfg.hFeedDown = ccdb->getForTimeStamp<TH2F>(paTHFD, timeStamp);
      if (cfg.hFeedDown == nullptr) {
        LOGF(fatal, "Could not load feed down histogram from %s", paTHFD.value.c_str());
      }
    }
    cfg.correctionsLoaded = true;
  }

  void loadNchCalibrations(uint64_t timeStamp)
  {
    if (paTHmeanNch.value.empty() == false) {
      cfgNch.hMeanNch = ccdb->getForTimeStamp<TH1F>(paTHmeanNch, timeStamp);
      if (cfgNch.hMeanNch == nullptr) {
        LOGF(fatal, "Could not load hMeanNch histogram from %s", paTHmeanNch.value.c_str());
      }
    }

    if (paTHsigmaNch.value.empty() == false) {
      cfgNch.hSigmaNch = ccdb->getForTimeStamp<TH1F>(paTHsigmaNch, timeStamp);
      if (cfgNch.hSigmaNch == nullptr) {
        LOGF(fatal, "Could not load hSigmaNch histogram from %s", paTHsigmaNch.value.c_str());
      }
    }
    cfgNch.calibrationsLoaded = true;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<UccZdc>(cfgc)};
}
