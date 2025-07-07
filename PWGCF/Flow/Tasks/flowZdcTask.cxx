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

/// \file   flowZdcTask.cxx
/// \author Sabrina Hernandez
/// \since  10/01/2024
/// \brief  task to evaluate flow and neutron skin with information from ZDC

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/MathConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"
#include <CCDB/BasicCCDBManager.h>

#include "TList.h"
#include <TComplex.h>
#include <TF1.h>
#include <TMath.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TRandom3.h>

#include <cmath>
#include <complex>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::mult;
using namespace o2::constants::math;
using namespace o2::aod::evsel;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct FlowZdcTask {
  SliceCache cache;

  O2_DEFINE_CONFIGURABLE(cfgCutVertex, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMin, float, 0.2f, "Minimal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMax, float, 10.0f, "Maximal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutChi2prTPCcls, float, 2.5, "Chi2 per TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAz, float, 2, "DCA Z cut")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAxy, float, 0.2f, "DCA XY cut")

  Configurable<int> eventSelection{"eventSelection", 1, "event selection"};
  Configurable<float> maxZp{"maxZp", 125.5, "Max ZP signal"};
  Configurable<float> maxZem{"maxZem", 3099.5, "Max ZEM signal"};
  // for ZDC info and analysis
  Configurable<int> nBinsAmp{"nBinsAmp", 1025, "nbinsAmp"};
  Configurable<float> maxZn{"maxZn", 125.5, "Max ZN signal"};
  Configurable<float> vtxRange{"vtxRange", 10.0f, "Vertex Z range to consider"};
  Configurable<float> etaRange{"etaRange", 1.0f, "Eta range to consider"};
  // configs for process QA
  Configurable<int> nBinsNch{"nBinsNch", 2501, "N bins Nch (|eta|<0.8)"};
  Configurable<int> nBinsAmpFT0{"nBinsAmpFT0", 100, "N bins FT0 amp"};
  Configurable<float> maxAmpFT0{"maxAmpFT0", 2500, "Max FT0 amp"};
  Configurable<int> nBinsAmpFV0{"nBinsAmpFV0", 100, "N bins FV0 amp"};
  Configurable<float> maxAmpFV0{"maxAmpFV0", 2000, "Max FV0 amp"};
  Configurable<int> nBinsZDC{"nBinsZDC", 400, "nBinsZDC"};
  Configurable<int> nBinsZN{"nBinsZN", 400, "N bins ZN"};
  Configurable<int> nBinsZP{"nBinsZP", 160, "N bins ZP"};
  Configurable<float> minNch{"minNch", 0, "Min Nch (|eta|<0.8)"};
  Configurable<float> maxNch{"maxNch", 2500, "Max Nch (|eta|<0.8)"};
  Configurable<int> nBinsTDC{"nBinsTDC", 150, "nbinsTDC"};
  Configurable<float> minTdc{"minTdc", -15.0, "minimum TDC"};
  Configurable<float> maxTdc{"maxTdc", 15.0, "maximum TDC"};
  Configurable<float> cfgCollisionEnergy{"cfgCollisionEnergy", 2.68, "cfgCollisionEnergy"};
  // event selection
  Configurable<bool> isNoCollInTimeRangeStrict{"isNoCollInTimeRangeStrict", true, "isNoCollInTimeRangeStrict?"};
  Configurable<bool> isNoCollInTimeRangeStandard{"isNoCollInTimeRangeStandard", false, "isNoCollInTimeRangeStandard?"};
  Configurable<bool> isNoCollInRofStrict{"isNoCollInRofStrict", true, "isNoCollInRofStrict?"};
  Configurable<bool> isNoCollInRofStandard{"isNoCollInRofStandard", false, "isNoCollInRofStandard?"};
  Configurable<bool> isNoHighMultCollInPrevRof{"isNoHighMultCollInPrevRof", true, "isNoHighMultCollInPrevRof?"};
  Configurable<bool> isNoCollInTimeRangeNarrow{"isNoCollInTimeRangeNarrow", false, "isNoCollInTimeRangeNarrow?"};
  Configurable<bool> isOccupancyCut{"isOccupancyCut", true, "Occupancy cut?"};
  Configurable<bool> isApplyFT0CbasedOccupancy{"isApplyFT0CbasedOccupancy", false, "T0C Occu cut?"};
  Configurable<bool> isTDCcut{"isTDCcut", false, "Use TDC cut?"};
  Configurable<bool> isZEMcut{"isZEMcut", false, "Use ZEM cut?"};
  Configurable<bool> useMidRapNchSel{"useMidRapNchSel", false, "Use mid-rapidity Nch selection"};
  Configurable<bool> applyEff{"applyEff", true, "Apply track-by-track efficiency correction"};
  Configurable<bool> applyFD{"applyFD", false, "Apply track-by-track feed down correction"};
  Configurable<bool> correctNch{"correctNch", true, "Correct also Nch"};

  Configurable<float> nSigmaNchCut{"nSigmaNchCut", 1., "nSigma Nch selection"};
  Configurable<double> minNchSel{"minNchSel", 5., "min Nch Selection"};
  Configurable<float> znBasedCut{"znBasedCut", 100, "ZN-based cut"};
  Configurable<float> zemCut{"zemCut", 1000., "ZEM cut"};
  Configurable<float> tdcCut{"tdcCut", 1., "TDC cut"};
  Configurable<float> minOccCut{"minOccCut", 0, "min Occu cut"};
  Configurable<float> maxOccCut{"maxOccCut", 500, "max Occu cut"};
  Configurable<float> minPt{"minPt", 0.1, "minimum pt of the tracks"};
  Configurable<float> maxPt{"maxPt", 3., "maximum pt of the tracks"};
  Configurable<float> maxPtSpectra{"maxPtSpectra", 50., "maximum pt of the tracks"};
  // axis configs
  ConfigurableAxis axisPhi{"axisPhi", {60, 0.0, constants::math::TwoPI}, "phi axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {3500, 0, 3500}, "centrality axis for histograms"};
  ConfigurableAxis axisZN{"axisZN", {5000, 0, 500}, "axisZN"};
  ConfigurableAxis axisZP{"axisZP", {5000, 0, 500}, "axisZP"};
  ConfigurableAxis axisFT0CAmp{"axisFT0CAmp", {5000, 0, 5000}, "axisFT0CAmp"};
  ConfigurableAxis axisFT0AAmp{"axisFT0AAmp", {5000, 0, 5000}, "axisFT0AAmp"};
  ConfigurableAxis axisFT0MAmp{"axisFT0MAmp", {10000, 0, 10000}, "axisFT0MAmp"};
  ConfigurableAxis multHistBin{"multHistBin", {501, -0.5, 500.5}, ""};
  ConfigurableAxis axisCent{"axisCent", {10, 0, 100}, "axisCent"};
  ConfigurableAxis ft0cMultHistBin{"ft0cMultHistBin", {501, -0.5, 500.5}, ""};
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.12}, "pT binning"};
  Configurable<float> posZcut{"posZcut", +10.0, "z-vertex position cut"};
  Configurable<float> minT0CcentCut{"minT0CcentCut", 0.0, "Min T0C Cent. cut"};
  Configurable<float> maxT0CcentCut{"maxT0CcentCut", 90.0, "Max T0C Cent. cut"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtMin) && (aod::track::pt < cfgCutPtMax) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && (aod::track::tpcChi2NCl < cfgCutChi2prTPCcls) && (nabs(aod::track::dcaZ) < cfgCutDCAz) && (nabs(aod::track::dcaXY) < cfgCutDCAxy);
  using ColEvSels = soa::Join<aod::Collisions, aod::EvSels, o2::aod::CentFT0Cs>;
  using AodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::Mults>>;
  using AodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA>>;
  Partition<AodTracks> tracksIUWithTPC = (aod::track::tpcNClsFindable > (uint8_t)0);
  using TracksSel = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCovIU, aod::TrackSelection, aod::TracksDCA>;
  using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
  using AodZDCs = soa::Join<aod::ZDCMults, aod::Zdcs>;
  using CollisionDataTable = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::CentFT0Ms>;
  using TrackDataTable = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;
  using FilTrackDataTable = soa::Filtered<TrackDataTable>;

  // CCDB paths
  Configurable<std::string> paTH{"paTH", "Users/s/sahernan/test", "base path to the ccdb object"};
  Configurable<std::string> paTHmeanNch{"paTHmeanNch", "Users/s/shernan/test", "base path to the ccdb object"};
  Configurable<std::string> paTHsigmaNch{"paTHsigmaNch", "Users/s/shernan/testSigma", "base path to the ccdb object"};
  Configurable<std::string> paTHEff{"paTHEff", "Users/s/shernan/TrackingEff", "base path to the ccdb object"};
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
  // Begin Histogram Registry

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Service<ccdb::BasicCCDBManager> ccdb;
  OutputObj<TProfile> pZNvsFT0Ccent{TProfile("pZNvsFT0Ccent", "ZN Energy vs FT0C Centrality", 100, 0, 100, 0, 500)};
  OutputObj<TProfile> pZPvsFT0Ccent{TProfile("pZPvsFT0Ccent", "ZP Energy vs FT0C Centrality", 100, 0, 100, 0, 500)};
  OutputObj<TProfile> pZNratiovscent{TProfile("pZNratiovscent", "Ratio ZNC/ZNA vs FT0C Centrality", 100, 0, 100, 0, 5)};
  OutputObj<TProfile> pZPratiovscent{TProfile("pZPratiovscent", "Ratio ZPC/ZPA vs FT0C Centrality", 100, 0, 100, 0, 5)};

  void init(InitContext const&)
  {
    // define axes
    const AxisSpec axisCounter{1, 0, +1, ""};
    const AxisSpec axisEvent{18, 0.5, 18.5, ""};
    const AxisSpec axisZpos{48, -12., 12., "Vtx_{z} (cm)"};
    const AxisSpec axisEta{40, -1., +1., "#eta"};
    const AxisSpec axisPt{binsPt, "#it{p}_{T} (GeV/#it{c})"};

    AxisSpec axisVtxZ{40, -20, 20, "Vertex Z", "VzAxis"};
    AxisSpec axisMult = {multHistBin, "Mult", "MultAxis"};
    AxisSpec axisFT0CMult = {ft0cMultHistBin, "ft0c", "FT0CMultAxis"};

    // create histograms
    histos.add("hEventCounter", "Event counter", kTH1F, {axisEvent});
    histos.add("zPos", ";;Entries;", kTH1F, {axisZpos});

    histos.add("eventCounter", "eventCounter", kTH1F, {axisCounter});
    histos.add("hZNvsFT0Ccent",
               "ZN Energy vs FT0C Centrality",
               kTH2F,
               {axisCent, axisZN});
    histos.add("hZPvsFT0Ccent",
               "ZP Energy vs FT0C Centrality;Centrality [%];ZP Energy",
               kTH2F,
               {axisCent, axisZP});
    histos.add("hNchvsNPV", ";NPVTracks (|#eta|<1);N_{ch} (|#eta|<0.8);",
               kTH2F,
               {{{nBinsNch, -0.5, maxNch}, {nBinsNch, -0.5, maxNch}}});
    histos.add("T0Ccent", ";;Entries", kTH1F, {axisCent});
    histos.add("NchUncorrected", ";#it{N}_{ch} (|#eta| < 0.8);Entries;", kTH1F, {{300, 0., 3000.}});
    histos.add("ZNamp", ";ZNA+ZNC;Entries;", kTH1F, {{nBinsZN, -0.5, maxZn}});
    histos.add("ExcludedEvtVsFT0M", ";T0A+T0C (#times 1/100, -3.3 < #eta < -2.1 and 3.5 < #eta < 4.9);Entries;", kTH1F, {{nBinsAmpFT0, 0., maxAmpFT0}});
    histos.add("ExcludedEvtVsNch", ";#it{N}_{ch} (|#eta|<0.8);Entries;", kTH1F, {{300, 0, 3000}});
    histos.add("Nch", ";#it{N}_{ch} (|#eta| < 0.8, Corrected);", kTH1F, {{nBinsNch, minNch, maxNch}});
    histos.add("NchVsOneParCorr", ";#it{N}_{ch} (|#eta| < 0.8, Corrected);#LT[#it{p}_{T}^{(1)}]#GT (GeV/#it{c})", kTProfile, {{nBinsNch, minNch, maxNch}});
    histos.add("EtaVsPhi", ";#eta;#varphi", kTH2F, {{{axisEta}, {100, -0.1 * PI, +2.1 * PI}}});
    histos.add("ZposVsEta", "", kTProfile, {axisZpos});
    histos.add("sigma1Pt", ";;#sigma(p_{T})/p_{T};", kTProfile, {axisPt});
    histos.add("dcaXYvspT", ";DCA_{xy} (cm);;", kTH2F, {{{50, -1., 1.}, {axisPt}}});
    histos.add("GlobalMult_vs_FT0C", "GlobalMult_vs_FT0C", kTH2F, {axisMult, axisFT0CMult});
    histos.add("VtxZHist", "VtxZHist", kTH1D, {axisVtxZ});

    // event selection steps
    histos.add("eventSelectionSteps", "eventSelectionSteps", kTH1D, {axisEvent});
    auto hstat = histos.get<TH1>(HIST("eventSelectionSteps"));
    auto* xAxis = hstat->GetXaxis();
    xAxis->SetBinLabel(1, "All events");
    xAxis->SetBinLabel(2, "SelEigth");
    xAxis->SetBinLabel(3, "NoSameBunchPileup"); // reject collisions in case of pileup with another collision in the same foundBC
    xAxis->SetBinLabel(4, "GoodZvtxFT0vsPV");   // small difference between z-vertex from PV and from FT0
    xAxis->SetBinLabel(5, "NoCollInTimeRangeStrict");
    xAxis->SetBinLabel(6, "NoCollInTimeRangeStandard");
    xAxis->SetBinLabel(7, "NoCollInRofStrict");
    xAxis->SetBinLabel(8, "NoCollInRofStandard");
    xAxis->SetBinLabel(9, "NoHighMultCollInPrevRof");
    xAxis->SetBinLabel(10, "NoCollInTimeRangeNarrow");
    xAxis->SetBinLabel(11, "Occupancy Cut");
    xAxis->SetBinLabel(12, "Cent. Sel.");
    xAxis->SetBinLabel(13, "VtxZ cut");
    xAxis->SetBinLabel(14, "has ZDC?");
    xAxis->SetBinLabel(15, "has T0?");
    xAxis->SetBinLabel(16, "Within TDC cut?");
    xAxis->SetBinLabel(17, "Within ZEM cut?");

    if (doprocessZdcCollAssoc) { // Check if the process function for ZDCCollAssoc is enabled
      histos.add("ZNAcoll", "ZNAcoll; ZNA amplitude; Entries", {HistType::kTH1F, {{nBinsAmp, -0.5, maxZn}}});
      histos.add("ZNCcoll", "ZNCcoll; ZNC amplitude; Entries", {HistType::kTH1F, {{nBinsAmp, -0.5, maxZn}}});
      histos.add("ZPCcoll", "ZPCcoll; ZPC amplitude; Entries", {HistType::kTH1F, {{nBinsAmp, -0.5, maxZn}}});
      histos.add("ZPAcoll", "ZPAcoll; ZPA amplitude; Entries", {HistType::kTH1F, {{nBinsAmp, -0.5, maxZn}}});
      histos.add("ZEM1coll", "ZEM1coll; ZEM1 amplitude; Entries", {HistType::kTH1F, {{nBinsAmp, -0.5, maxZem}}});
      histos.add("ZEM2coll", "ZEM2coll; ZEM2 amplitude; Entries", {HistType::kTH1F, {{nBinsAmp, -0.5, maxZem}}});
      histos.add("ZNvsZEMcoll", "ZNvsZEMcoll; ZEM; ZNA+ZNC", {HistType::kTH2F, {{{nBinsAmp, -0.5, maxZem}, {nBinsAmp, -0.5, 2. * maxZn}}}});
      histos.add("ZNAvsZNCcoll", "ZNAvsZNCcoll; ZNC; ZNA", {HistType::kTH2F, {{{nBinsAmp, -0.5, maxZn}, {nBinsAmp, -0.5, maxZn}}}});
      histos.add("ZDC_energy_vs_ZEM", "ZDCvsZEM; ZEM; ZNA+ZNC+ZPA+ZPC", {HistType::kTH2F, {{{nBinsAmp, -0.5, maxZem}, {nBinsAmp, -0.5, 2. * maxZn}}}});
      // common energies information for ZDC
      histos.add("ZNCenergy", "common sum ZN energy side c", kTH1F, {axisZN});
      histos.add("ZNAenergy", "common sum ZN energy side a", kTH1F, {axisZN});
      histos.add("ZPCenergy", "common sum ZP energy side c", kTH1F, {axisZP});
      histos.add("ZPAenergy", "common sum ZP energy side a", kTH1F, {axisZP});
      histos.add("ZNenergy", "common sum zn (a + c sides) energy", kTH1F, {axisZN});
      histos.add("ZPenergy", "common sum zp energy (a + c sides)", kTH1F, {axisZP});
      histos.add("hZNvsFT0CAmp", "ZN Energy vs FT0C Amplitude", kTH2F, {axisFT0CAmp, axisZN});
      histos.add("hZPvsFT0CAmp", "ZP Energy vs FT0C Amplitude", kTH2F, {axisFT0CAmp, axisZP});
      histos.add("hZNvsMult", "ZN Energy vs Multiplicity", kTH2F, {axisMultiplicity, axisZN});
      histos.add("hZPvsMult", "ZP Energy vs Multiplicity", kTH2F, {axisMultiplicity, axisZP});
    }

    if (doprocessQA) {
      histos.add("ZNVsFT0A", ";T0A (#times 1/100);ZNA+ZNC;", kTH2F, {{{nBinsAmpFT0, 0., maxAmpFT0}, {nBinsZDC, -0.5, maxZn}}});
      histos.add("ZNVsFT0C", ";T0C (#times 1/100);ZNA+ZNC;", kTH2F, {{{nBinsAmpFT0, 0., maxAmpFT0}, {nBinsZDC, -0.5, maxZn}}});
      histos.add("ZNVsFT0M", ";T0A+T0C (#times 1/100);ZNA+ZNC;", kTH2F, {{{nBinsAmpFT0, 0., 3000.}, {nBinsZDC, -0.5, maxZn}}});
      histos.add("ZN", ";ZNA+ZNC;Entries;", kTH1F, {{nBinsZDC, -0.5, maxZn}});
      histos.add("ZNA", ";ZNA;Entries;", kTH1F, {{nBinsZDC, -0.5, maxZn}});
      histos.add("ZPA", ";ZPA;Entries;", kTH1F, {{nBinsZDC, -0.5, maxZp}});
      histos.add("ZNC", ";ZNC;Entries;", kTH1F, {{nBinsZDC, -0.5, maxZn}});
      histos.add("ZPC", ";ZPC;Entries;", kTH1F, {{nBinsZDC, -0.5, maxZp}});
      histos.add("ZNAVsZNC", ";ZNC;ZNA", kTH2F, {{{30, -0.5, maxZn}, {30, -0.5, maxZn}}});
      histos.add("ZPAVsZPC", ";ZPC;ZPA;", kTH2F, {{{100, -0.5, maxZp}, {100, -0.5, maxZp}}});
      histos.add("ZNAVsZPA", ";ZPA;ZNA;", kTH2F, {{{20, -0.5, maxZp}, {30, -0.5, maxZn}}});
      histos.add("ZNCVsZPC", ";ZPC;ZNC;", kTH2F, {{{20, -0.5, maxZp}, {30, -0.5, maxZn}}});
      histos.add("ZNASector", ";ZNA;Entries;", kTH1F, {{nBinsZDC, -0.5, maxZn}});
      histos.add("ZPASector", ";ZPA;Entries;", kTH1F, {{nBinsZDC, -0.5, maxZp}});
      histos.add("ZNCSector", ";ZNC;Entries;", kTH1F, {{nBinsZDC, -0.5, maxZn}});
      histos.add("ZPCSector", ";ZPC;Entries;", kTH1F, {{nBinsZDC, -0.5, maxZp}});
      histos.add("ZNCcvsZNCsum", ";ZNC common;ZNC sum towers;", kTH2F, {{{30, -0.5, maxZn}, {30, -0.5, maxZn}}});
      histos.add("ZNAcvsZNAsum", ";ZNA common;ZNA sum towers;", kTH2F, {{{30, -0.5, maxZn}, {30, -0.5, maxZn}}});
      histos.add("ZPCcvsZPCsum", ";ZPC common;ZPC sum towers;", kTH2F, {{{30, -0.5, maxZp}, {30, -0.5, maxZp}}});
      histos.add("ZPAcvsZPAsum", ";ZPA common;ZPA sum towers;", kTH2F, {{{30, -0.5, maxZp}, {30, -0.5, maxZp}}});
      histos.add("ZNVsZEM", ";ZEM;ZNA+ZNC;", kTH2F, {{{60, -0.5, maxZem}, {60, -0.5, maxZn}}});
      histos.add("ZNCVstdc", ";t_{ZNC};ZNC;", kTH2F, {{{30, -15., 15.}, {nBinsZDC, -0.5, maxZn}}});
      histos.add("ZNAVstdc", ";t_{ZNA};ZNA;", kTH2F, {{{30, -15., 15.}, {30, -0.5, maxZn}}});
      histos.add("ZPCVstdc", ";t_{ZPC};ZPC;", kTH2F, {{{30, -15., 15}, {20, -0.5, maxZp}}});
      histos.add("ZPAVstdc", ";t_{ZPA};ZPA;", kTH2F, {{{30, -15., 15.}, {20, -0.5, maxZp}}});
      histos.add("ZEM1Vstdc", ";t_{ZEM1};ZEM1;", kTH2F, {{{30, -15., 15.}, {30, -0.5, 2000.5}}});
      histos.add("ZEM2Vstdc", ";t_{ZEM2};ZEM2;", kTH2F, {{{30, -15., 15.}, {30, -0.5, 2000.5}}});
      histos.add("debunch", ";t_{ZDC}-t_{ZDA};t_{ZDC}+t_{ZDA}", kTH2F, {{{nBinsTDC, minTdc, maxTdc}, {nBinsTDC, minTdc, maxTdc}}});

      histos.add("GlbTracks", "Nch", kTH1F, {{nBinsNch, minNch, maxNch}});
      histos.add("NchVsFT0C", ";T0C (#times 1/100, -3.3 < #eta < -2.1);#it{N}_{ch} (|#eta|<0.8);", kTH2F, {{{nBinsAmpFT0, 0., 950.}, {nBinsNch, minNch, maxNch}}});
      histos.add("NchVsFT0M", ";T0A+T0C (#times 1/100, -3.3 < #eta < -2.1 and 3.5 < #eta < 4.9);#it{N}_{ch} (|#eta|<0.8);", kTH2F, {{{nBinsAmpFT0, 0., 3000.}, {nBinsNch, minNch, maxNch}}});
      histos.add("NchVsFT0A", ";T0A (#times 1/100, 3.5 < #eta < 4.9);#it{N}_{ch} (|#eta|<0.8);", kTH2F, {{{nBinsAmpFT0, 0., maxAmpFT0}, {nBinsNch, minNch, maxNch}}});
      histos.add("NchVsFV0A", ";V0A (#times 1/100, 2.2 < #eta < 5);#it{N}_{ch} (|#eta|<0.8);", kTH2F, {{{nBinsAmpFV0, 0., maxAmpFV0}, {nBinsNch, minNch, maxNch}}});

      histos.add("NchVsEt", ";#it{E}_{T} (|#eta|<0.8);#LTITS+TPC tracks#GT (|#eta|<0.8);", kTH2F, {{{nBinsNch, minNch, maxNch}, {nBinsNch, minNch, maxNch}}});
      histos.add("NchVsMeanPt", ";#it{N}_{ch} (|#eta|<0.8);#LT[#it{p}_{T}]#GT (|#eta|<0.8);", kTProfile, {{nBinsNch, minNch, maxNch}});
      histos.add("NchVsNPV", ";#it{N}_{PV} (|#eta|<1);ITS+TPC tracks (|#eta|<0.8);", kTH2F, {{{300, -0.5, 5999.5}, {nBinsNch, minNch, maxNch}}});
      histos.add("NchVsITStracks", ";ITS tracks nCls >= 5;TITS+TPC tracks (|#eta|<0.8);", kTH2F, {{{300, -0.5, 5999.5}, {nBinsNch, minNch, maxNch}}});
      histos.add("ZNCVsNch", ";#it{N}_{ch} (|#eta|<0.8);ZNC;", kTH2F, {{{nBinsNch, minNch, maxNch}, {nBinsZDC, minNch, maxZn}}});
      histos.add("ZNAVsNch", ";#it{N}_{ch} (|#eta|<0.8);ZNA;", kTH2F, {{{nBinsNch, minNch, maxNch}, {nBinsZDC, minNch, maxZn}}});
      histos.add("ZNVsNch", ";#it{N}_{ch} (|#eta|<0.8);ZNA+ZNC;", kTH2F, {{{nBinsNch, minNch, maxNch}, {nBinsZDC, minNch, maxZn}}});
      histos.add("ZNDifVsNch", ";#it{N}_{ch} (|#eta|<0.8);ZNA-ZNC;", kTH2F, {{{nBinsNch, minNch, maxNch}, {100, -50., 50.}}});
    }

    ccdb->setURL("http://alice-ccdb.cern.ch");
    // Enabling object caching, otherwise each call goes to the CCDB server
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    // Not later than now, will be replaced by the value of the train creation
    // This avoids that users can replace objects **while** a train is running
    ccdb->setCreatedNotAfter(ccdbNoLaterThan.value);
  }
  template <typename EventCuts>
  bool isEventSelected(EventCuts const& col)
  {
    histos.fill(HIST("hEventCounter"), EvCutLabel::All);
    if (!col.sel8()) {
      return false;
    }
    histos.fill(HIST("hEventCounter"), EvCutLabel::SelEigth);

    if (!col.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    histos.fill(HIST("hEventCounter"), EvCutLabel::NoSameBunchPileup);

    if (!col.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    histos.fill(HIST("hEventCounter"), EvCutLabel::IsGoodZvtxFT0vsPV);

    if (isNoCollInTimeRangeStrict) {
      if (!col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict)) {
        return false;
      }
      histos.fill(HIST("hEventCounter"), EvCutLabel::NoCollInTimeRangeStrict);
    }

    if (isNoCollInTimeRangeStandard) {
      if (!col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
        return false;
      }
      histos.fill(HIST("hEventCounter"), EvCutLabel::NoCollInTimeRangeStandard);
    }

    if (isNoCollInRofStrict) {
      if (!col.selection_bit(o2::aod::evsel::kNoCollInRofStrict)) {
        return false;
      }
      histos.fill(HIST("hEventCounter"), EvCutLabel::NoCollInRofStrict);
    }

    if (isNoCollInRofStandard) {
      if (!col.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
        return false;
      }
      histos.fill(HIST("hEventCounter"), EvCutLabel::NoCollInRofStandard);
    }

    if (isNoHighMultCollInPrevRof) {
      if (!col.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
        return false;
      }
      histos.fill(HIST("hEventCounter"), EvCutLabel::NoHighMultCollInPrevRof);
    }

    // To be used in combination with FT0C-based occupancy
    if (isNoCollInTimeRangeNarrow) {
      if (!col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow)) {
        return false;
      }
      histos.fill(HIST("hEventCounter"), EvCutLabel::NoCollInTimeRangeNarrow);
    }

    if (isOccupancyCut) {
      auto occuValue{isApplyFT0CbasedOccupancy ? col.ft0cOccupancyInTimeRange() : col.trackOccupancyInTimeRange()};
      if (occuValue < minOccCut || occuValue > maxOccCut) {
        return false;
      }
    }
    histos.fill(HIST("hEventCounter"), EvCutLabel::OccuCut);

    if (col.centFT0C() < minT0CcentCut || col.centFT0C() > maxT0CcentCut) {
      return false;
    }
    histos.fill(HIST("hEventCounter"), EvCutLabel::Centrality);

    // Z-vertex position cut
    if (std::fabs(col.posZ()) > posZcut) {
      return false;
    }
    histos.fill(HIST("hEventCounter"), EvCutLabel::VtxZ);

    return true;
  }

  void processQA(ColEvSels::iterator const& collision, BCsRun3 const& /*bcs*/, aod::Zdcs const& /*zdcsData*/, aod::FV0As const& /*fv0as*/, aod::FT0s const& /*ft0s*/, AodTracks const& tracks)
  {
    const auto& foundBC = collision.foundBC_as<BCsRun3>();
    if (!isEventSelected(collision)) {
      return;
    }
    if (!foundBC.has_zdc()) {
      return;
    }
    histos.fill(HIST("hEventCounter"), EvCutLabel::Zdc);
    auto zdc = foundBC.zdc();

    float aT0A = 0., aT0C = 0., aV0A = 0.;
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
    histos.fill(HIST("hEventCounter"), EvCutLabel::TZero);
    if (foundBC.has_fv0a()) {
      for (const auto& amplitude : foundBC.fv0a().amplitude()) {
        aV0A += amplitude;
      }
    } else {
      aV0A = -999.;
    }
    float tZNA{zdc.timeZNA()};
    float tZNC{zdc.timeZNC()};
    float tZPA{zdc.timeZPA()};
    float tZPC{zdc.timeZPC()};
    float tZDCdif{tZNC + tZPC - tZNA - tZPA};
    float tZDCsum{tZNC + tZPC + tZNA + tZPA};
    const double normT0M{(aT0A + aT0C) / 100.};
    float znA = zdc.amplitudeZNA() / cfgCollisionEnergy;
    float znC = zdc.amplitudeZNC() / cfgCollisionEnergy;
    float zpA = zdc.amplitudeZPA() / cfgCollisionEnergy;
    float zpC = zdc.amplitudeZPC() / cfgCollisionEnergy;
    float aZEM1{zdc.amplitudeZEM1()};
    float aZEM2{zdc.amplitudeZEM2()};
    float sumZEMs{aZEM1 + aZEM2};
    float tZEM1{zdc.timeZEM1()};
    float tZEM2{zdc.timeZEM2()};
    float sumZNs{znA + znC};

    // TDC cut
    if (isTDCcut) {
      if (std::sqrt(std::pow(tZDCdif, 2.) + std::pow(tZDCsum, 2.)) > tdcCut) {
        return;
      }
      histos.fill(HIST("hEventCounter"), EvCutLabel::Tdc);
    }

    // ZEM cut
    if (isZEMcut) {
      if (sumZEMs < zemCut) {
        return;
      }
      histos.fill(HIST("hEventCounter"), EvCutLabel::Zem);
    }

    float sumZNC = (zdc.energySectorZNC())[0] + (zdc.energySectorZNC())[1] + (zdc.energySectorZNC())[2] + (zdc.energySectorZNC())[3];
    float sumZNA = (zdc.energySectorZNA())[0] + (zdc.energySectorZNA())[1] + (zdc.energySectorZNA())[2] + (zdc.energySectorZNA())[3];
    float sumZPC = (zdc.energySectorZPC())[0] + (zdc.energySectorZPC())[1] + (zdc.energySectorZPC())[2] + (zdc.energySectorZPC())[3];
    float sumZPA = (zdc.energySectorZPA())[0] + (zdc.energySectorZPA())[1] + (zdc.energySectorZPA())[2] + (zdc.energySectorZPA())[3];

    float et = 0., meanpt = 0.;
    int itsTracks = 0, glbTracks = 0;
    for (const auto& track : tracks) {
      if (track.hasITS()) {
        itsTracks++;
      }
      // Track Selection
      if (!track.isGlobalTrack()) {
        continue;
      }
      if ((track.pt() < minPt) || (track.pt() > maxPt)) {
        continue;
      }
      glbTracks++;
    }
    bool skipEvent{false};
    if (useMidRapNchSel) {
      auto hMeanNch = ccdb->getForTimeStamp<TH1F>(paTHmeanNch.value, foundBC.timestamp());
      auto hSigmaNch = ccdb->getForTimeStamp<TH1F>(paTHsigmaNch.value, foundBC.timestamp());
      if (!hMeanNch) {
        LOGF(info, "hMeanNch NOT LOADED!");
        return;
      }
      if (!hSigmaNch) {
        LOGF(info, "hSigmaNch NOT LOADED!");
        return;
      }

      const int binT0M{hMeanNch->FindBin(normT0M)};
      const double meanNch{hMeanNch->GetBinContent(binT0M)};
      const double sigmaNch{hSigmaNch->GetBinContent(binT0M)};
      const double nSigmaSelection{nSigmaNchCut * sigmaNch};
      const double diffMeanNch{meanNch - glbTracks};

      if (!(std::abs(diffMeanNch) < nSigmaSelection)) {
        histos.fill(HIST("ExcludedEvtVsNch"), glbTracks);
      } else {
        skipEvent = true;
      }
    } else {
      skipEvent = true;
    }
    if (!skipEvent) {
      return;
    }

    for (const auto& track : tracks) {
      // Track Selection
      if (!track.isGlobalTrack()) {
        continue;
      }
      if ((track.pt() < minPt) || (track.pt() > maxPtSpectra)) {
        continue;
      }

      histos.fill(HIST("ZposVsEta"), collision.posZ(), track.eta());
      histos.fill(HIST("EtaVsPhi"), track.eta(), track.phi());
      histos.fill(HIST("dcaXYvspT"), track.dcaXY(), track.pt());
      et += std::sqrt(std::pow(track.pt(), 2.) + std::pow(o2::constants::physics::MassPionCharged, 2.));
      meanpt += track.pt();
    }
    histos.fill(HIST("zPos"), collision.posZ());
    histos.fill(HIST("T0Ccent"), collision.centFT0C());

    histos.fill(HIST("ZNCcvsZNCsum"), sumZNC / cfgCollisionEnergy, zdc.energyCommonZNC() / cfgCollisionEnergy);
    histos.fill(HIST("ZNAcvsZNAsum"), sumZNA / cfgCollisionEnergy, zdc.energyCommonZNA() / cfgCollisionEnergy);
    histos.fill(HIST("ZPCcvsZPCsum"), sumZPC / cfgCollisionEnergy, zdc.energyCommonZPC() / cfgCollisionEnergy);
    histos.fill(HIST("ZPAcvsZPAsum"), sumZPA / cfgCollisionEnergy, zdc.energyCommonZPA() / cfgCollisionEnergy);

    histos.fill(HIST("GlbTracks"), glbTracks);
    histos.fill(HIST("ZNA"), znA);
    histos.fill(HIST("ZNC"), znC);
    histos.fill(HIST("ZPA"), zpA);
    histos.fill(HIST("ZPC"), zpC);
    histos.fill(HIST("ZNASector"), sumZNA / cfgCollisionEnergy);
    histos.fill(HIST("ZNCSector"), sumZNC / cfgCollisionEnergy);
    histos.fill(HIST("ZPASector"), sumZPA / cfgCollisionEnergy);
    histos.fill(HIST("ZPCSector"), sumZPC / cfgCollisionEnergy);
    histos.fill(HIST("ZN"), znA + znC);
    histos.fill(HIST("ZNAVsZNC"), znC, znA);
    histos.fill(HIST("ZNAVsZPA"), zpA, znA);
    histos.fill(HIST("ZNCVsZPC"), zpC, znC);
    histos.fill(HIST("ZPAVsZPC"), zpC, zpA);
    histos.fill(HIST("ZNVsZEM"), sumZEMs, sumZNs);
    histos.fill(HIST("ZNCVstdc"), tZNC, znC);
    histos.fill(HIST("ZNAVstdc"), tZNA, znA);
    histos.fill(HIST("ZPCVstdc"), tZPC, zpC);
    histos.fill(HIST("ZPAVstdc"), tZPA, zpA);
    histos.fill(HIST("ZEM1Vstdc"), tZEM1, aZEM1);
    histos.fill(HIST("ZEM2Vstdc"), tZEM2, aZEM2);
    histos.fill(HIST("debunch"), tZDCdif, tZDCsum);

    histos.fill(HIST("ZNVsFT0A"), aT0A / 100., sumZNs);
    histos.fill(HIST("ZNVsFT0C"), aT0C / 100., sumZNs);
    histos.fill(HIST("ZNVsFT0M"), (aT0A + aT0C) / 100., sumZNs);

    if (sumZNs > znBasedCut) {
      return;
    }
    histos.fill(HIST("NchVsFV0A"), aV0A / 100., glbTracks);
    histos.fill(HIST("NchVsFT0A"), aT0A / 100., glbTracks);
    histos.fill(HIST("NchVsFT0C"), aT0C / 100., glbTracks);
    histos.fill(HIST("NchVsFT0M"), (aT0A + aT0C) / 100., glbTracks);

    histos.fill(HIST("NchVsEt"), et, glbTracks);
    histos.fill(HIST("NchVsITStracks"), itsTracks, glbTracks);
    histos.fill(HIST("ZNAVsNch"), glbTracks, znA);
    histos.fill(HIST("ZNCVsNch"), glbTracks, znC);
    histos.fill(HIST("ZNVsNch"), glbTracks, sumZNs);
    histos.fill(HIST("ZNDifVsNch"), glbTracks, znA - znC);
    if (glbTracks >= minNchSel) {
      histos.fill(HIST("NchVsMeanPt"), glbTracks, meanpt / glbTracks);
    }
  }

  void processZdcCollAssoc(
    AodCollisions::iterator const& collision,
    AodTracks const& tracks,
    BCsRun3 const& /*bcs*/,
    aod::Zdcs const& /*zdcs*/,
    aod::FT0s const& /*ft0s*/)
  {
    if (!isEventSelected(collision)) {
      return;
    }
    const auto& foundBC = collision.foundBC_as<BCsRun3>();
    if (!foundBC.has_zdc()) {
      return;
    }
    int nTot = tracks.size();
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
    const double normT0M{(ft0aAmp + ft0aAmp) / 100.};

    const auto& zdcread = foundBC.zdc();
    const auto cent = collision.centFT0C();

    // ZDC data and histogram filling
    float znA = zdcread.amplitudeZNA();
    float znC = zdcread.amplitudeZNC();
    float zpA = zdcread.amplitudeZPA();
    float zpC = zdcread.amplitudeZPC();
    float tZNA{zdcread.timeZNA()};
    float tZNC{zdcread.timeZNC()};
    float tZPA{zdcread.timeZPA()};
    float tZPC{zdcread.timeZPC()};
    float tZDCdif{tZNC + tZPC - tZNA - tZPA};
    float tZDCsum{tZNC + tZPC + tZNA + tZPA};
    float sumZNC = (zdcread.energySectorZNC())[0] + (zdcread.energySectorZNC())[1] + (zdcread.energySectorZNC())[2] + (zdcread.energySectorZNC())[3];
    float sumZNA = (zdcread.energySectorZNA())[0] + (zdcread.energySectorZNA())[1] + (zdcread.energySectorZNA())[2] + (zdcread.energySectorZNA())[3];
    float sumZPC = (zdcread.energySectorZPC())[0] + (zdcread.energySectorZPC())[1] + (zdcread.energySectorZPC())[2] + (zdcread.energySectorZPC())[3];
    float sumZPA = (zdcread.energySectorZPA())[0] + (zdcread.energySectorZPA())[1] + (zdcread.energySectorZPA())[2] + (zdcread.energySectorZPA())[3];
    float sumZDC = sumZPA + sumZPC + sumZNA + sumZNC;
    float sumZEM = zdcread.amplitudeZEM1() + zdcread.amplitudeZEM2();
    znA /= cfgCollisionEnergy;
    znC /= cfgCollisionEnergy;
    zpA /= cfgCollisionEnergy;
    zpC /= cfgCollisionEnergy;
    float sumZNs{znA + znC};
    float sumZPs{zpA + zpC};
    // TDC cut
    if (isTDCcut) {
      if (std::sqrt(std::pow(tZDCdif, 2.) + std::pow(tZDCsum, 2.)) > tdcCut) {
        return;
      }
      histos.fill(HIST("hEventCounter"), EvCutLabel::Tdc);
    }
    // ZEM cut
    if (isZEMcut) {
      if (sumZEM < zemCut) {
        return;
      }
    }
    // common energies
    float commonSumZnc = (zdcread.energyCommonZNC());
    float commonSumZna = (zdcread.energyCommonZNA());
    float commonSumZpc = (zdcread.energyCommonZPC());
    float commonSumZpa = (zdcread.energyCommonZPA());
    float sumZN = (sumZNC) + (sumZNA);
    float sumZP = (sumZPC) + (sumZPA);

    int itsTracks = 0, glbTracks = 0;
    for (const auto& track : tracks) {
      // Track Selection
      if (track.hasITS()) {
        itsTracks++;
      }
      if (!track.isGlobalTrack()) {
        continue;
      }
      if ((track.pt() < minPt) || (track.pt() > maxPt)) {
        continue;
      }
      histos.fill(HIST("ZposVsEta"), collision.posZ(), track.eta());
      histos.fill(HIST("EtaVsPhi"), track.eta(), track.phi());
      histos.fill(HIST("dcaXYvspT"), track.dcaXY(), track.pt());
      glbTracks++;
    }
    bool skipEvent{false};
    if (useMidRapNchSel) {
      auto hMeanNch = ccdb->getForTimeStamp<TH1F>(paTHmeanNch.value, foundBC.timestamp());
      auto hSigmaNch = ccdb->getForTimeStamp<TH1F>(paTHsigmaNch.value, foundBC.timestamp());
      if (!hMeanNch) {
        LOGF(info, "hMeanNch NOT LOADED!");
        return;
      }
      if (!hSigmaNch) {
        LOGF(info, "hSigmaNch NOT LOADED!");
        return;
      }
      const int binT0M{hMeanNch->FindBin(normT0M)};
      const double meanNch{hMeanNch->GetBinContent(binT0M)};
      const double sigmaNch{hSigmaNch->GetBinContent(binT0M)};
      const double nSigmaSelection{nSigmaNchCut * sigmaNch};
      const double diffMeanNch{meanNch - glbTracks};
      if (!(std::abs(diffMeanNch) < nSigmaSelection)) {
        histos.fill(HIST("ExcludedEvtVsFT0M"), normT0M);
        histos.fill(HIST("ExcludedEvtVsNch"), glbTracks);
      } else {
        skipEvent = true;
      }
    }
    // Skip event based on number of Nch sigmas
    if (!skipEvent) {
      return;
    }
    std::vector<float> vecOneOverEff;
    auto efficiency = ccdb->getForTimeStamp<TH1F>(paTHEff.value, foundBC.timestamp());
    if (!efficiency) {
      return;
    }
    // Calculates the Nch multiplicity
    for (const auto& track : tracks) {
      // Track Selection
      if (!track.isGlobalTrack()) {
        continue;
      }
      if ((track.pt() < minPt) || (track.pt() > maxPt)) {
        continue;
      }

      float pt{track.pt()};
      float effValue{1.0};
      if (applyEff) {
        effValue = efficiency->GetBinContent(efficiency->FindBin(pt));
      }
      if (effValue > 0.) {
        vecOneOverEff.emplace_back(1. / effValue);
      }
    }

    double nchMult{0.};
    nchMult = std::accumulate(vecOneOverEff.begin(), vecOneOverEff.end(), 0);
    if (!applyEff)
      nchMult = static_cast<double>(glbTracks);
    if (applyEff && !correctNch)
      nchMult = static_cast<double>(glbTracks);
    if (nchMult < minNchSel) {
      return;
    }
    histos.get<TH2>(HIST("ZNvsZEMcoll"))->Fill(zdcread.amplitudeZEM1() + zdcread.amplitudeZEM2(), zdcread.amplitudeZNA() + zdcread.amplitudeZNC());
    histos.get<TH2>(HIST("ZNAvsZNCcoll"))->Fill(zdcread.amplitudeZNC(), zdcread.amplitudeZNA());
    histos.get<TH1>(HIST("ZEM1coll"))->Fill(zdcread.amplitudeZEM1());
    histos.get<TH1>(HIST("ZEM2coll"))->Fill(zdcread.amplitudeZEM2());
    histos.fill(HIST("ZNenergy"), sumZN);
    histos.fill(HIST("ZPenergy"), sumZP);
    histos.fill(HIST("ZNCenergy"), commonSumZnc);
    histos.fill(HIST("ZNAenergy"), commonSumZna);
    histos.fill(HIST("ZPAenergy"), commonSumZpa);
    histos.fill(HIST("ZPCenergy"), commonSumZpc);
    histos.fill(HIST("hZNvsFT0Ccent"), cent, sumZN);
    histos.fill(HIST("hZPvsFT0Ccent"), cent, sumZP);
    histos.fill(HIST("hZNvsFT0CAmp"), ft0cAmp, sumZN);
    histos.fill(HIST("hZPvsFT0CAmp"), ft0cAmp, sumZP);
    histos.fill(HIST("hZNvsMult"), nTot, sumZN);
    histos.fill(HIST("hZPvsMult"), nTot, sumZP);
    histos.fill(HIST("hNchvsNPV"), collision.multNTracksPVeta1(), nTot);
    histos.fill(HIST("Nch"), nchMult);
    histos.fill(HIST("ZNamp"), sumZNs);
    histos.fill(HIST("NchVsZN"), nchMult, sumZNs);
    histos.fill(HIST("NchVsZP"), nchMult, sumZPs);
    histos.fill(HIST("NITSTacksVsZN"), itsTracks, sumZNs);
    histos.fill(HIST("NITSTacksVsZP"), itsTracks, sumZPs);
    histos.fill(HIST("T0MVsZN"), normT0M, sumZNs);
    histos.fill(HIST("T0MVsZP"), normT0M, sumZPs);
    histos.fill(HIST("NchUncorrected"), glbTracks);

    float ratioZN = sumZNC / sumZNA;
    float ratioZP = sumZPC / sumZPA;
    pZNratiovscent->Fill(cent, ratioZN);
    pZPratiovscent->Fill(cent, ratioZP);
    pZNvsFT0Ccent->Fill(cent, sumZN);
    pZPvsFT0Ccent->Fill(cent, sumZP);
    histos.get<TH2>(HIST("ZDC_energy_vs_ZEM"))->Fill(sumZEM, sumZDC);
  }

  void processCorrelation(CollisionDataTable::iterator const& collision, FilTrackDataTable const& tracks)
  {
    if (!isEventSelected(collision)) {
      return;
    }
    if (std::abs(collision.posZ()) >= vtxRange) {
      return;
    }
    histos.fill(HIST("VtxZHist"), collision.posZ());
    auto nchTracks = 0;
    for (const auto& track : tracks) {
      if (std::abs(track.eta()) >= etaRange) {
        continue;
      }
      nchTracks++;
    }
    histos.fill(HIST("GlobalMult_vs_FT0C"), nchTracks, collision.multFT0C());
  }

  PROCESS_SWITCH(FlowZdcTask, processZdcCollAssoc, "Processing ZDC w. collision association", false);
  PROCESS_SWITCH(FlowZdcTask, processQA, "Process QA", true);
  PROCESS_SWITCH(FlowZdcTask, processCorrelation, "Process correlations", true);

}; // end of struct function

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowZdcTask>(cfgc)};
}
