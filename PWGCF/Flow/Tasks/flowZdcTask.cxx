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
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
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
  Configurable<int> eventSelection{"eventSelection", 1, "event selection"};
  Configurable<float> maxZp{"maxZp", 125.5, "Max ZP signal"};
  Configurable<float> maxZem{"maxZem", 3099.5, "Max ZEM signal"};
  // for ZDC info and analysis
  Configurable<int> nBinsAmp{"nBinsAmp", 1025, "nbinsAmp"};
  Configurable<int> nBinsADC{"nBinsADC", 1000, "nbinsADC"};
  Configurable<int> nBinsCent{"nBinsCent", 90, "nBinsCent"};
  Configurable<float> maxZn{"maxZn", 125.5, "Max ZN signal"};
  // configs for process QA
  Configurable<int> nBinsNch{"nBinsNch", 2501, "N bins Nch (|eta|<0.8)"};
  Configurable<int> nBinsAmpFT0{"nBinsAmpFT0", 100, "N bins FT0 amp"};
  Configurable<float> maxAmpFT0{"maxAmpFT0", 2500, "Max FT0 amp"};
  Configurable<float> maxAmpFT0M{"maxAmpFT0M", 2500, "Max FT0M amp"};
  Configurable<int> nBinsAmpFV0{"nBinsAmpFV0", 100, "N bins FV0 amp"};
  Configurable<float> maxAmpFV0{"maxAmpFV0", 2000, "Max FV0 amp"};
  Configurable<int> nBinsZDC{"nBinsZDC", 400, "nBinsZDC"};
  Configurable<float> minNch{"minNch", 0, "Min Nch (|eta|<0.8)"};
  Configurable<float> maxNch{"maxNch", 2500, "Max Nch (|eta|<0.8)"};
  Configurable<int> nBinsTDC{"nBinsTDC", 150, "nbinsTDC"};
  Configurable<float> minTdcZn{"minTdcZn", -4.0, "minimum TDC for ZN"};
  Configurable<float> maxTdcZn{"maxTdcZn", -4.0, "maximum TDC for ZN"};
  Configurable<float> minTdcZp{"minTdcZp", -4.0, "minimum TDC for ZP"};
  Configurable<float> maxTdcZp{"maxTdcZp", -4.0, "maximum TDC for ZP"};
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
  Configurable<bool> useMidRapNchSel{"useMidRapNchSel", false, "Use mid-rapidity Nch selection"};
  Configurable<bool> applyEff{"applyEff", true, "Apply track-by-track efficiency correction"};
  Configurable<bool> correctNch{"correctNch", true, "Correct also Nch"};

  Configurable<float> nSigmaNchCut{"nSigmaNchCut", 1., "nSigma Nch selection"};
  Configurable<double> minNchSel{"minNchSel", 5., "min Nch Selection"};
  Configurable<float> tdcCut{"tdcCut", 1., "TDC cut"};
  Configurable<float> minOccCut{"minOccCut", 0, "min Occu cut"};
  Configurable<float> maxOccCut{"maxOccCut", 500, "max Occu cut"};
  Configurable<float> minPt{"minPt", 0.1, "minimum pt of the tracks"};
  Configurable<float> maxPt{"maxPt", 3., "maximum pt of the tracks"};
  Configurable<float> maxPtSpectra{"maxPtSpectra", 50., "maximum pt of the tracks"};
  Configurable<float> zemCut{"zemCut", 100., "ZEM cut"};
  // axis configs
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {3500, 0, 3500}, "centrality axis for histograms"};
  ConfigurableAxis axisZN{"axisZN", {5000, 0, 500}, "axisZN"};
  ConfigurableAxis axisZP{"axisZP", {5000, 0, 500}, "axisZP"};
  ConfigurableAxis axisCent{"axisCent", {10, 0, 100}, "axisCent"};
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.12}, "pT binning"};
  Configurable<float> posZcut{"posZcut", +10.0, "z-vertex position cut"};
  Configurable<float> minEta{"minEta", -0.8, "minimum eta"};
  Configurable<float> maxEta{"maxEta", +0.8, "maximum eta"};
  Configurable<float> minT0CcentCut{"minT0CcentCut", 0.0, "Min T0C Cent. cut"};
  Configurable<float> maxT0CcentCut{"maxT0CcentCut", 90.0, "Max T0C Cent. cut"};

  using ColEvents = soa::Join<aod::Collisions, aod::EvSels>;
  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = ((aod::track::eta > minEta) && (aod::track::eta < maxEta));
  using ColEvSels = soa::Join<aod::Collisions, aod::EvSels, o2::aod::CentFT0Cs, aod::TPCMults, o2::aod::BarrelMults, aod::FT0MultZeqs>;
  using AodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::Mults>>;
  using AodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA>>;
  Partition<AodTracks> tracksIUWithTPC = (aod::track::tpcNClsFindable > (uint8_t)0);
  using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
  using TracksSel = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCovIU, aod::TrackSelection, aod::TracksDCA>;
  using CollisionDataTable = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::CentFT0Ms>;
  using TrackDataTable = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;
  using FilTrackDataTable = soa::Filtered<TrackDataTable>;
  using TheFilteredTracks = soa::Filtered<TracksSel>;

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
    const AxisSpec axisEvent{18, 0.5, 18.5, ""};
    const AxisSpec axisZpos{48, -12., 12., "Vtx_{z} (cm)"};
    const AxisSpec axisEta{40, -1., +1., "#eta"};
    const AxisSpec axisPt{binsPt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec axisVtxZ{40, -20, 20, "Vertex Z", "VzAxis"};

    // create histograms
    histos.add("hEventCounter", "Event counter", kTH1F, {axisEvent});
    histos.add("zPos", ";;Entries;", kTH1F, {axisZpos});
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
    histos.add("ZNamp", ";ZNA+ZNC;Entries;", kTH1F, {{nBinsZDC, -0.5, maxZn}});
    histos.add("ExcludedEvtVsFT0M", ";T0A+T0C (#times 1/100, -3.3 < #eta < -2.1 and 3.5 < #eta < 4.9);Entries;", kTH1F, {{nBinsAmpFT0, 0., maxAmpFT0}});
    histos.add("ExcludedEvtVsNch", ";#it{N}_{ch} (|#eta|<0.8);Entries;", kTH1F, {{300, 0, 3000}});
    histos.add("Nch", ";#it{N}_{ch} (|#eta| < 0.8, Corrected);", kTH1F, {{nBinsNch, minNch, maxNch}});
    histos.add("EtaVsPhi", ";#eta;#varphi", kTH2F, {{{axisEta}, {100, -0.1 * PI, +2.1 * PI}}});
    histos.add("ZposVsEta", "", kTProfile, {axisZpos});
    histos.add("dcaXYvspT", ";DCA_{xy} (cm);;", kTH2F, {{{50, -1., 1.}, {axisPt}}});

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
      histos.add("hZNvsFT0CAmp", "ZN Energy vs FT0C Amplitude", kTH2F, {{nBinsAmpFT0, 0., maxAmpFT0}, axisZN});
      histos.add("hZPvsFT0CAmp", "ZP Energy vs FT0C Amplitude", kTH2F, {{nBinsAmpFT0, 0., maxAmpFT0}, axisZP});
      histos.add("hZNvsMult", "ZN Energy vs Multiplicity", kTH2F, {axisMultiplicity, axisZN});
      histos.add("hZPvsMult", "ZP Energy vs Multiplicity", kTH2F, {axisMultiplicity, axisZP});
    }

    if (doprocessQA) {
      histos.add("ZNVsFT0A", ";T0A (#times 1/100);ZNA+ZNC Amplitude;", kTH2F, {{{nBinsAmpFT0, 0., maxAmpFT0}, {nBinsZDC, -0.5, maxZn}}});
      histos.add("ZNVsFT0C", ";T0C (#times 1/100);ZNA+ZNC Amplitude;", kTH2F, {{{nBinsAmpFT0, 0., maxAmpFT0}, {nBinsZDC, -0.5, maxZn}}});
      histos.add("ZNVsFT0M", ";T0A+T0C (#times 1/100);ZNA+ZNC Amplitude;", kTH2F, {{{nBinsAmpFT0, 0., maxAmpFT0M}, {nBinsZDC, -0.5, maxZn}}});
      histos.add("ZPVsFT0A", ";T0A (#times 1/100);ZPA+ZPC Amplitude;", kTH2F, {{{nBinsAmpFT0, 0., maxAmpFT0}, {nBinsZDC, -0.5, maxZp}}});
      histos.add("ZPVsFT0C", ";T0C (#times 1/100);ZPA+ZPC Amplitude;", kTH2F, {{{nBinsAmpFT0, 0., maxAmpFT0}, {nBinsZDC, -0.5, maxZp}}});
      histos.add("ZPVsFT0M", ";T0A+T0C (#times 1/100);ZPA+ZPC Amplitude;", kTH2F, {{{nBinsAmpFT0, 0., maxAmpFT0M}, {nBinsZDC, -0.5, maxZp}}});
      histos.add("ZNAVsFT0A", ";T0A (#times 1/100);ZNA Amplitude;", kTH2F, {{{nBinsAmpFT0, 0., maxAmpFT0}, {nBinsZDC, -0.5, maxZn}}});
      histos.add("ZNAVsFT0C", ";T0C (#times 1/100);ZNA Amplitude;", kTH2F, {{{nBinsAmpFT0, 0., maxAmpFT0}, {nBinsZDC, -0.5, maxZn}}});
      histos.add("ZNAVsFT0M", ";T0A+T0C (#times 1/100);ZNA Amplitude;", kTH2F, {{{nBinsAmpFT0, 0., maxAmpFT0M}, {nBinsZDC, -0.5, maxZn}}});
      histos.add("ZNCVsFT0A", ";T0A (#times 1/100);ZNC Amplitude;", kTH2F, {{{nBinsAmpFT0, 0., maxAmpFT0}, {nBinsZDC, -0.5, maxZn}}});
      histos.add("ZNCVsFT0C", ";T0C (#times 1/100);ZNC Amplitude;", kTH2F, {{{nBinsAmpFT0, 0., maxAmpFT0}, {nBinsZDC, -0.5, maxZn}}});
      histos.add("ZNCVsFT0M", ";T0A+T0C (#times 1/100);ZNC Amplitude;", kTH2F, {{{nBinsAmpFT0, 0., maxAmpFT0M}, {nBinsZDC, -0.5, maxZn}}});
      histos.add("ZPAVsFT0A", ";T0A (#times 1/100);ZPA Amplitude;", kTH2F, {{{nBinsAmpFT0, 0., maxAmpFT0}, {nBinsZDC, -0.5, maxZp}}});
      histos.add("ZPAVsFT0C", ";T0C (#times 1/100);ZPA Amplitude;", kTH2F, {{{nBinsAmpFT0, 0., maxAmpFT0}, {nBinsZDC, -0.5, maxZp}}});
      histos.add("ZPAVsFT0M", ";T0A+T0C (#times 1/100);ZPA Amplitude;", kTH2F, {{{nBinsAmpFT0, 0., maxAmpFT0M}, {nBinsZDC, -0.5, maxZp}}});
      histos.add("ZPCVsFT0A", ";T0A (#times 1/100);ZPC Amplitude;", kTH2F, {{{nBinsAmpFT0, 0., maxAmpFT0}, {nBinsZDC, -0.5, maxZp}}});
      histos.add("ZPCVsFT0C", ";T0C (#times 1/100);ZPC Amplitude;", kTH2F, {{{nBinsAmpFT0, 0., maxAmpFT0}, {nBinsZDC, -0.5, maxZp}}});
      histos.add("ZPCVsFT0M", ";T0A+T0C (#times 1/100);ZPC Amplitude;", kTH2F, {{{nBinsAmpFT0, 0., maxAmpFT0M}, {nBinsZDC, -0.5, maxZp}}});
      histos.add("ZN", ";ZNA+ZNC;Entries;", kTH1F, {{nBinsZDC, -0.5, maxZn}});
      histos.add("ZNA", ";ZNA Amplitude;Entries;", kTH1F, {{nBinsZDC, -0.5, maxZn}});
      histos.add("ZPA", ";ZPA Amplitude;Entries;", kTH1F, {{nBinsZDC, -0.5, maxZp}});
      histos.add("ZNC", ";ZNC Amplitude;Entries;", kTH1F, {{nBinsZDC, -0.5, maxZn}});
      histos.add("ZPC", ";ZPC Amplitude;Entries;", kTH1F, {{nBinsZDC, -0.5, maxZp}});
      histos.add("ZNACommon", ";ZNA Common Energy;Entries;", kTH1F, {{nBinsZDC, -0.5, maxZn}});
      histos.add("ZPACommon", ";ZPA Common Energy;Entries;", kTH1F, {{nBinsZDC, -0.5, maxZp}});
      histos.add("ZNCCommon", ";ZNC Common Energy;Entries;", kTH1F, {{nBinsZDC, -0.5, maxZn}});
      histos.add("ZPCCommon", ";ZPC Common Energy;Entries;", kTH1F, {{nBinsZDC, -0.5, maxZp}});
      histos.add("ZNAVsZNC", ";ZNC;ZNA", kTH2F, {{{nBinsZDC, -0.5, maxZn}, {nBinsZDC, -0.5, maxZn}}});
      histos.add("ZPAVsZPC", ";ZPC;ZPA;", kTH2F, {{{nBinsZDC, -0.5, maxZp}, {nBinsZDC, -0.5, maxZp}}});
      histos.add("ZNAVsZPA", ";ZPA;ZNA;", kTH2F, {{{nBinsZDC, -0.5, maxZp}, {nBinsZDC, -0.5, maxZn}}});
      histos.add("ZNCVsZPC", ";ZPC;ZNC;", kTH2F, {{{nBinsZDC, -0.5, maxZp}, {nBinsZDC, -0.5, maxZn}}});
      histos.add("ZNASector", ";ZNA;Entries;", kTH1F, {{nBinsZDC, -0.5, maxZn}});
      histos.add("ZPASector", ";ZPA;Entries;", kTH1F, {{nBinsZDC, -0.5, maxZp}});
      histos.add("ZNCSector", ";ZNC;Entries;", kTH1F, {{nBinsZDC, -0.5, maxZn}});
      histos.add("ZPCSector", ";ZPC;Entries;", kTH1F, {{nBinsZDC, -0.5, maxZp}});
      histos.add("ZNCcvsZNCsum", ";ZNC common;ZNC sum towers;", kTH2F, {{{30, -0.5, maxZn}, {30, -0.5, maxZn}}});
      histos.add("ZNAcvsZNAsum", ";ZNA common;ZNA sum towers;", kTH2F, {{{30, -0.5, maxZn}, {30, -0.5, maxZn}}});
      histos.add("ZPCcvsZPCsum", ";ZPC common;ZPC sum towers;", kTH2F, {{{30, -0.5, maxZp}, {30, -0.5, maxZp}}});
      histos.add("ZPAcvsZPAsum", ";ZPA common;ZPA sum towers;", kTH2F, {{{30, -0.5, maxZp}, {30, -0.5, maxZp}}});
      histos.add("ZNVsZEM", ";ZEM;ZNA+ZNC;", kTH2F, {{{60, -0.5, maxZem}, {60, -0.5, maxZn}}});
      histos.add("ZNCVstdccoll", ";t_{ZNC};ZNC;", kTH2F, {{{nBinsTDC, -13.5, 11.45}, {nBinsZDC, -0.5, maxZn}}});
      histos.add("ZNAVstdccoll", ";t_{ZNA};ZNA;", kTH2F, {{{nBinsTDC, -13.5, 11.45}, {nBinsZDC, -0.5, maxZn}}});
      histos.add("ZPCVstdccoll", ";t_{ZPC};ZPC;", kTH2F, {{{nBinsTDC, -13.5, 11.45}, {nBinsZDC, -0.5, maxZp}}});
      histos.add("ZPAVstdccoll", ";t_{ZPA};ZPA;", kTH2F, {{{nBinsTDC, -13.5, 11.45}, {nBinsZDC, -0.5, maxZp}}});
      histos.add("ZEM1", ";ZEM1 Amplitude;Entries;", kTH1F, {{nBinsZDC, -0.5, maxZem}});
      histos.add("ZEM2", ";ZEM2 Amplitude;Entries;", kTH1F, {{nBinsZDC, -0.5, maxZem}});
      histos.add("ZEM1Vstdc", ";t_{ZEM1};ZEM1;", kTH2F, {{{480, -13.5, 11.45}, {30, -0.5, 2000.5}}});
      histos.add("ZEM2Vstdc", ";t_{ZEM2};ZEM2;", kTH2F, {{{480, -13.5, 11.45}, {30, -0.5, 2000.5}}});
      histos.add("debunch", ";t_{ZDC}-t_{ZDA};t_{ZDC}+t_{ZDA}", kTH2F, {{{nBinsTDC, minTdcZn, maxTdcZn}, {nBinsTDC, minTdcZp, maxTdcZp}}});
      histos.add("GlbTracks", "Nch", kTH1F, {{nBinsNch, minNch, maxNch}});
      histos.add("ampFT0C", ";T0C (#times 1/100);", kTH1F, {{nBinsAmpFT0, 0., maxAmpFT0}});
      histos.add("ampFT0A", ";T0A (#times 1/100);", kTH1F, {{nBinsAmpFT0, 0., maxAmpFT0}});
      histos.add("ampFT0M", ";T0A+T0C (#times 1/100);", kTH1F, {{nBinsAmpFT0, 0., maxAmpFT0M}});
      histos.add("ampFV0A", ";V0A (#times 1/100);", kTH1F, {{nBinsAmpFV0, 0., maxAmpFV0}});
      histos.add("NchVsFT0C", ";T0C (#times 1/100, -3.3 < #eta < -2.1);#it{N}_{ch} (|#eta|<0.8);", kTH2F, {{{nBinsAmpFT0, 0., 950.}, {nBinsNch, minNch, maxNch}}});
      histos.add("NchVsFT0M", ";T0A+T0C (#times 1/100, -3.3 < #eta < -2.1 and 3.5 < #eta < 4.9);#it{N}_{ch} (|#eta|<0.8);", kTH2F, {{{nBinsAmpFT0, 0., maxAmpFT0M}, {nBinsNch, minNch, maxNch}}});
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
      histos.add("ZPAvsCent", ";centFT0C;ZPA", kTH2F, {{{axisCent}, {nBinsZDC, -0.5, maxZp}}});
      histos.add("ZPCvsCent", ";centFT0C;ZPC", kTH2F, {{{axisCent}, {nBinsZDC, -0.5, maxZp}}});
      histos.add("pZPAvsFT0Ccent", ";FT0C centrality;ZPA Amplitude", kTProfile, {{nBinsCent, minT0CcentCut, maxT0CcentCut}});
      histos.add("pZPCvsFT0Ccent", ";FT0C centrality;ZPC Amplitude", kTProfile, {{nBinsCent, minT0CcentCut, maxT0CcentCut}});
      histos.add("pZPAvsGlbTrack", ";Global Tracks (ITS + TPC);ZPA Amplitude", kTProfile, {{nBinsNch, minNch, maxNch}});
      histos.add("pZPCvsGlbTrack", ";Global Tracks (ITS + TPC);ZPC Amplitude", kTProfile, {{nBinsNch, minNch, maxNch}});
      histos.add("hZPASectorvsGlbTrack", ";Global Tracks (ITS + TPC);ZPA Sector Energy", kTH2F, {{{nBinsNch, minNch, maxNch}, {nBinsZDC, minNch, maxZp}}});
      histos.add("hZPCSectorvsGlbTrack", ";Global Tracks (ITS + TPC);ZPC Sector Energy", kTH2F, {{{nBinsNch, minNch, maxNch}, {nBinsZDC, minNch, maxZp}}});
      histos.add("hZNASectorvsGlbTrack", ";Global Tracks (ITS + TPC);ZNA Sector Energy", kTH2F, {{{nBinsNch, minNch, maxNch}, {nBinsZDC, minNch, maxZn}}});
      histos.add("hZNCSectorvsGlbTrack", ";Global Tracks (ITS + TPC);ZNC Sector Energy", kTH2F, {{{nBinsNch, minNch, maxNch}, {nBinsZDC, minNch, maxZn}}});
      histos.add("hZPSectorvsGlbTrack", ";Global Tracks (ITS + TPC);(ZPA + ZPC) Sector Energy", kTH2F, {{{nBinsNch, minNch, maxNch}, {nBinsZDC, minNch, maxZp}}});
      histos.add("hZNSectorvsGlbTrack", ";Global Tracks (ITS + TPC);(ZNA + ZNC) Sector Energy", kTH2F, {{{nBinsNch, minNch, maxNch}, {nBinsZDC, minNch, maxZn}}});
    }
    if (doprocessZdc) {
      histos.add("ampZna", ";ZNA Amplitude;Entries;", kTH1F, {{nBinsZDC, -0.5, maxZn}});
      histos.add("ampZpa", ";ZPA Amplitude;Entries;", kTH1F, {{nBinsZDC, -0.5, maxZp}});
      histos.add("ampZnc", ";ZNC Amplitude;Entries;", kTH1F, {{nBinsZDC, -0.5, maxZn}});
      histos.add("ampZpc", ";ZPC Amplitude;Entries;", kTH1F, {{nBinsZDC, -0.5, maxZp}});
      histos.add("ampZEM1", ";ZEM1 Amplitude;Entries;", kTH1F, {{nBinsZDC, -0.5, maxZem}});
      histos.add("ampZEM2", ";ZEM2 Amplitude;Entries;", kTH1F, {{nBinsZDC, -0.5, maxZem}});
      histos.add("ZnVsZem", "ZnVsZEM; ZEM; ZNA + ZNC", kTH2F, {{{nBinsZDC, -0.5, maxZem}, {nBinsZDC, -0.5, maxZn}}});
      histos.add("ZnaVsZnc", "ZNAvsZNC; ZNC; ZNA;", kTH2F, {{{nBinsZDC, -0.5, maxZn}, {nBinsZDC, -0.5, maxZn}}});
      histos.add("ZpaVsZpc", "ZPAvsZPC; ZPC; ZPA;", kTH2F, {{{nBinsZDC, -0.5, maxZp}, {nBinsZDC, -0.5, maxZp}}});
      histos.add("ZnaVsZpa", "ZNAvsZPA; ZNA; ZPA;", kTH2F, {{{nBinsZDC, -0.5, maxZn}, {nBinsZDC, -0.5, maxZp}}});
      histos.add("ZncVsZpc", "ZNCvsZPC; ZNC; ZPC;", kTH2F, {{{nBinsZDC, -0.5, maxZn}, {nBinsZDC, -0.5, maxZp}}});
      histos.add("ZnccVsZncSum", "ZNCcVsZNCsum; ZNCC ADC; ZNCsum", kTH2F, {{{nBinsADC, -0.5, 3. * maxZn}, {nBinsADC, -0.5, 3. * maxZn}}});
      histos.add("ZnacVsZnaSum", "ZNAcVsZNAsum; ZNAC ADC; ZNAsum", kTH2F, {{{nBinsADC, -0.5, 3. * maxZn}, {nBinsADC, -0.5, 3. * maxZn}}});
      histos.add("ZpacVsZpaSum", "ZPAcVsZPAsum; ZPAC ADC; ZPAsum", kTH2F, {{{nBinsADC, -0.5, 3. * maxZp}, {nBinsADC, -0.5, 3. * maxZp}}});
      histos.add("ZpccVsZpcSum", "ZPCcVsZPCsum; ZPCC ADC; ZPCsum", kTH2F, {{{nBinsADC, -0.5, 3. * maxZp}, {nBinsADC, -0.5, 3. * maxZp}}});
      histos.add("ZncVsTdc", "ZNCvsTDC; ZNC Amp; ZNC TDC", kTH2F, {{{480, -13.5, 11.45}, {nBinsZDC, -0.5, maxZn}}});
      histos.add("ZnaVsTdc", "ZNAvsTDC; ZNA Amp; ZNA TDC", kTH2F, {{{480, -13.5, 11.45}, {nBinsZDC, -0.5, maxZn}}});
      histos.add("ZpcVsTdc", "ZPCvsTDC; ZPC Amp; ZPC TDC", kTH2F, {{{480, -13.5, 11.45}, {nBinsZDC, -0.5, maxZp}}});
      histos.add("ZpaVsTdc", "ZPAvsTDC; ZPA Amp; ZPA TDC", kTH2F, {{{480, -13.5, 11.45}, {nBinsZDC, -0.5, maxZp}}});
      histos.add("Zem1VsTdc", "ZEM1vsTDC; ZEM1 Amp; ZEM1 TDC", kTH2F, {{{480, -13.5, 11.45}, {nBinsZDC, -0.5, maxZem}}});
      histos.add("Zem2VsTdc", "ZEM2vsTDC; ZEM2 Amp; ZEM2 TDC", kTH2F, {{{480, -13.5, 11.45}, {nBinsZDC, -0.5, maxZem}}});
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

  void processQA(ColEvSels::iterator const& collision, BCsRun3 const& /*bcs*/, aod::Zdcs const& /*zdcsData*/, aod::FV0As const& /*fv0as*/, aod::FT0s const& /*ft0s*/, TheFilteredTracks const& tracks)
  {
    const auto& foundBC = collision.foundBC_as<BCsRun3>();
    const auto cent = collision.centFT0C();
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
    }
    histos.fill(HIST("hEventCounter"), EvCutLabel::TZero);
    if (foundBC.has_fv0a()) {
      for (const auto& amplitude : foundBC.fv0a().amplitude()) {
        aV0A += amplitude;
      }
    }
    float tZNA{zdc.timeZNA()};
    float tZNC{zdc.timeZNC()};
    float tZPA{zdc.timeZPA()};
    float tZPC{zdc.timeZPC()};
    const double normT0M{(aT0A + aT0C) / 100.};
    float znA = zdc.amplitudeZNA() / cfgCollisionEnergy;
    float znC = zdc.amplitudeZNC() / cfgCollisionEnergy;
    float zpA = zdc.amplitudeZPA() / cfgCollisionEnergy;
    float zpC = zdc.amplitudeZPC() / cfgCollisionEnergy;
    float commonSumZnc = zdc.energyCommonZNC() / cfgCollisionEnergy;
    float commonSumZna = zdc.energyCommonZNA() / cfgCollisionEnergy;
    float commonSumZpc = zdc.energyCommonZPC() / cfgCollisionEnergy;
    float commonSumZpa = zdc.energyCommonZPA() / cfgCollisionEnergy;
    float aZEM1{zdc.amplitudeZEM1()};
    float aZEM2{zdc.amplitudeZEM2()};
    float sumZEMs{aZEM1 + aZEM2};
    float tZEM1{zdc.timeZEM1()};
    float tZEM2{zdc.timeZEM2()};
    float sumZNs{znA + znC};
    float sumZNC = ((zdc.energySectorZNC())[0] + (zdc.energySectorZNC())[1] + (zdc.energySectorZNC())[2] + (zdc.energySectorZNC())[3]) / cfgCollisionEnergy;
    float sumZNA = ((zdc.energySectorZNA())[0] + (zdc.energySectorZNA())[1] + (zdc.energySectorZNA())[2] + (zdc.energySectorZNA())[3]) / cfgCollisionEnergy;
    float sumZPC = ((zdc.energySectorZPC())[0] + (zdc.energySectorZPC())[1] + (zdc.energySectorZPC())[2] + (zdc.energySectorZPC())[3]) / cfgCollisionEnergy;
    float sumZPA = ((zdc.energySectorZPA())[0] + (zdc.energySectorZPA())[1] + (zdc.energySectorZPA())[2] + (zdc.energySectorZPA())[3]) / cfgCollisionEnergy;
    float sumSectZN = (sumZNC + sumZNA);
    float sumSectZP = (sumZPC + sumZPA);

    if (sumZEMs > zemCut) {
      if (isTDCcut) {
        if ((tZNA >= minTdcZn) && (tZNA <= maxTdcZn)) {
          histos.fill(HIST("ZNA"), znA);
          histos.fill(HIST("ZNACommon"), commonSumZna);
          histos.fill(HIST("ZNASector"), sumZNA);
        }
        if ((tZNC >= minTdcZn) && (tZNC <= maxTdcZn)) {
          histos.fill(HIST("ZNC"), znC);
          histos.fill(HIST("ZNCCommon"), commonSumZnc);
          histos.fill(HIST("ZNCSector"), sumZNC);
        }
        if ((tZPA >= minTdcZp) && (tZPA <= maxTdcZp)) {
          histos.fill(HIST("ZPA"), zpA);
          histos.fill(HIST("ZPACommon"), commonSumZpa);
          histos.fill(HIST("ZPASector"), sumZPA);
        }
        if ((tZPC >= minTdcZp) && (tZPC <= maxTdcZp)) {
          histos.fill(HIST("ZPC"), zpC);
          histos.fill(HIST("ZPCCommon"), commonSumZpc);
          histos.fill(HIST("ZPCSector"), sumZPC);
        }
        if (((tZNA >= minTdcZn) && (tZNA <= maxTdcZn)) && ((tZNC >= minTdcZn) && (tZNC <= maxTdcZn)))
          histos.fill(HIST("ZNVsZEM"), sumZEMs, sumZNs);
        if (((tZNA >= minTdcZn) && (tZNA <= maxTdcZn)) && ((tZNC >= minTdcZn) && (tZNC <= maxTdcZn))) {
          histos.fill(HIST("ZNAVsZNC"), znC, znA);
          histos.fill(HIST("ZN"), znA + znC);
        }
        if ((tZNA >= minTdcZn) && (tZNA <= maxTdcZn))
          histos.fill(HIST("ZNAVsZPA"), zpA, znA);
        if ((tZNC >= minTdcZn) && (tZNC <= maxTdcZn))
          histos.fill(HIST("ZNCVsZPC"), zpC, znC);
        if (((tZPA >= minTdcZp) && (tZPA <= maxTdcZp)) && ((tZPC >= minTdcZp) && (tZPC <= maxTdcZp)))
          histos.fill(HIST("ZPAVsZPC"), zpC, zpA);
      } else {
        histos.fill(HIST("ZNA"), znA);
        histos.fill(HIST("ZNC"), znC);
        histos.fill(HIST("ZPA"), zpA);
        histos.fill(HIST("ZPC"), zpC);
        histos.fill(HIST("ZNVsZEM"), sumZEMs, sumZNs);
        histos.fill(HIST("ZNAVsZNC"), znC, znA);
        histos.fill(HIST("ZNAVsZPA"), zpA, znA);
        histos.fill(HIST("ZNCVsZPC"), zpC, znC);
        histos.fill(HIST("ZPAVsZPC"), zpC, zpA);
        histos.fill(HIST("ZNACommon"), commonSumZna);
        histos.fill(HIST("ZNASector"), sumZNA);
        histos.fill(HIST("ZNCCommon"), commonSumZnc);
        histos.fill(HIST("ZNCSector"), sumZNC);
        histos.fill(HIST("ZPACommon"), commonSumZpa);
        histos.fill(HIST("ZPASector"), sumZPA);
        histos.fill(HIST("ZPCCommon"), commonSumZpc);
        histos.fill(HIST("ZPCSector"), sumZPC);
        histos.fill(HIST("ZN"), znA + znC);
      }
      histos.fill(HIST("ZEM1"), aZEM1);
      histos.fill(HIST("ZEM2"), aZEM2);
      histos.fill(HIST("ZNCVstdccoll"), tZNC, znC);
      histos.fill(HIST("ZNAVstdccoll"), tZNA, znA);
      histos.fill(HIST("ZPCVstdccoll"), tZPC, zpC);
      histos.fill(HIST("ZPAVstdccoll"), tZPA, zpA);
      histos.fill(HIST("ZEM1Vstdc"), tZEM1, aZEM1);
      histos.fill(HIST("ZEM2Vstdc"), tZEM2, aZEM2);
      histos.fill(HIST("debunch"), tZNA - tZNC, tZNA + tZNC);
    }
    float et = 0., meanpt = 0.;
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
      if ((track.eta() < minEta) || (track.eta() > maxEta)) {
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
    histos.fill(HIST("GlbTracks"), glbTracks);

    if (sumZEMs > zemCut) {
      histos.fill(HIST("ZNVsFT0C"), aT0C / 100., znA + znC);
      histos.fill(HIST("ZNVsFT0M"), (aT0A + aT0C) / 100., znA + znC);
      histos.fill(HIST("ZPVsFT0A"), aT0A / 100., zpA + zpC);
      histos.fill(HIST("ZPVsFT0C"), aT0C / 100., zpA + zpC);
      histos.fill(HIST("ZPVsFT0M"), (aT0A + aT0C) / 100., zpA + zpC);
      histos.fill(HIST("ZPAVsFT0A"), aT0A / 100., zpA);
      histos.fill(HIST("ZPAVsFT0C"), aT0C / 100., zpA);
      histos.fill(HIST("ZPAVsFT0M"), (aT0A + aT0C) / 100., zpA);
      histos.fill(HIST("ZPCVsFT0A"), aT0A / 100., zpC);
      histos.fill(HIST("ZPCVsFT0C"), aT0C / 100., zpC);
      histos.fill(HIST("ZPCVsFT0M"), (aT0A + aT0C) / 100., zpC);
      histos.fill(HIST("ZNCVsFT0A"), aT0A / 100., znC);
      histos.fill(HIST("ZNCVsFT0C"), aT0C / 100., znC);
      histos.fill(HIST("ZNCVsFT0M"), (aT0A + aT0C) / 100., znC);
      histos.fill(HIST("ZNAVsFT0A"), aT0A / 100., znA);
      histos.fill(HIST("ZNAVsFT0C"), aT0C / 100., znA);
      histos.fill(HIST("ZNAVsFT0M"), (aT0A + aT0C) / 100., znA);
      histos.fill(HIST("ZPAvsCent"), cent, zpA);
      histos.fill(HIST("ZPCvsCent"), cent, zpC);
      if (std::isfinite(zpA) && !std::isnan(zpA) && cent >= minT0CcentCut && cent < maxT0CcentCut && glbTracks >= minNch && glbTracks < maxNch) {
        histos.fill(HIST("pZPAvsFT0Ccent"), cent, zpA);
        histos.fill(HIST("pZPAvsGlbTrack"), glbTracks, zpA);
        histos.fill(HIST("hZPASectorvsGlbTrack"), glbTracks, sumZPA);
      }
      if (std::isfinite(zpC) && !std::isnan(zpC) && cent >= minT0CcentCut && cent < maxT0CcentCut && glbTracks >= minNch && glbTracks < maxNch) {
        histos.fill(HIST("pZPCvsFT0Ccent"), cent, zpC);
        histos.fill(HIST("pZPCvsGlbTrack"), glbTracks, zpC);
        histos.fill(HIST("hZPCSectorvsGlbTrack"), glbTracks, sumZPC);
      }
      histos.fill(HIST("hZNASectorvsGlbTrack"), glbTracks, sumZNA);
      histos.fill(HIST("hZNCSectorvsGlbTrack"), glbTracks, sumZNC);
      histos.fill(HIST("hZPSectorvsGlbTrack"), glbTracks, sumSectZP);
      histos.fill(HIST("hZNSectorvsGlbTrack"), glbTracks, sumSectZN);
      // ZDC Correlations
      histos.fill(HIST("ZNAVsNch"), glbTracks, znA);
      histos.fill(HIST("ZNCVsNch"), glbTracks, znC);
      histos.fill(HIST("ZNVsNch"), glbTracks, sumZNs);
      histos.fill(HIST("ZNDifVsNch"), glbTracks, znA - znC);
      histos.fill(HIST("ZNCcvsZNCsum"), sumZNC, zdc.energyCommonZNC());
      histos.fill(HIST("ZNAcvsZNAsum"), sumZNA, zdc.energyCommonZNA());
      histos.fill(HIST("ZPCcvsZPCsum"), sumZPC, zdc.energyCommonZPC());
      histos.fill(HIST("ZPAcvsZPAsum"), sumZPA, zdc.energyCommonZPA());
    }

    histos.fill(HIST("ampFT0C"), aT0C / 100.);
    histos.fill(HIST("ampFT0A"), aT0A / 100.);
    histos.fill(HIST("ampFT0M"), (aT0A + aT0C) / 100.);
    histos.fill(HIST("ampFV0A"), aV0A / 100.);
    // charged particle correlations
    histos.fill(HIST("NchVsFV0A"), aV0A / 100., glbTracks);
    histos.fill(HIST("NchVsFT0A"), aT0A / 100., glbTracks);
    histos.fill(HIST("NchVsFT0C"), aT0C / 100., glbTracks);
    histos.fill(HIST("NchVsFT0M"), (aT0A + aT0C) / 100., glbTracks);
    histos.fill(HIST("hNchvsNPV"), collision.multNTracksPVeta1(), tracks.size());
    histos.fill(HIST("NchVsEt"), et, glbTracks);
    histos.fill(HIST("NchVsITStracks"), itsTracks, glbTracks);
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
    float znA = zdcread.amplitudeZNA() / cfgCollisionEnergy;
    float znC = zdcread.amplitudeZNC() / cfgCollisionEnergy;
    float zpA = zdcread.amplitudeZPA() / cfgCollisionEnergy;
    float zpC = zdcread.amplitudeZPC() / cfgCollisionEnergy;
    float tZNA{zdcread.timeZNA()};
    float tZNC{zdcread.timeZNC()};
    float tZPA{zdcread.timeZPA()};
    float tZPC{zdcread.timeZPC()};
    float tZDCdif{tZNC + tZPC - tZNA - tZPA};
    float tZDCsum{tZNC + tZPC + tZNA + tZPA};
    float sumZNC = ((zdcread.energySectorZNC())[0] + (zdcread.energySectorZNC())[1] + (zdcread.energySectorZNC())[2] + (zdcread.energySectorZNC())[3]) / cfgCollisionEnergy;
    float sumZNA = ((zdcread.energySectorZNA())[0] + (zdcread.energySectorZNA())[1] + (zdcread.energySectorZNA())[2] + (zdcread.energySectorZNA())[3]) / cfgCollisionEnergy;
    float sumZPC = ((zdcread.energySectorZPC())[0] + (zdcread.energySectorZPC())[1] + (zdcread.energySectorZPC())[2] + (zdcread.energySectorZPC())[3]) / cfgCollisionEnergy;
    float sumZPA = ((zdcread.energySectorZPA())[0] + (zdcread.energySectorZPA())[1] + (zdcread.energySectorZPA())[2] + (zdcread.energySectorZPA())[3]) / cfgCollisionEnergy;
    float sumZDC = sumZPA + sumZPC + sumZNA + sumZNC;
    float sumZEM = zdcread.amplitudeZEM1() + zdcread.amplitudeZEM2();
    float sumZNs{znA + znC};
    float sumZPs{zpA + zpC};
    // TDC cut
    if (isTDCcut) {
      if (std::sqrt(std::pow(tZDCdif, 2.) + std::pow(tZDCsum, 2.)) > tdcCut) {
        return;
      }
      histos.fill(HIST("hEventCounter"), EvCutLabel::Tdc);
    }
    // common energies
    float commonSumZnc = zdcread.energyCommonZNC() / cfgCollisionEnergy;
    float commonSumZna = zdcread.energyCommonZNA() / cfgCollisionEnergy;
    float commonSumZpc = zdcread.energyCommonZPC() / cfgCollisionEnergy;
    float commonSumZpa = zdcread.energyCommonZPA() / cfgCollisionEnergy;
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

  void processZdc(
    ColEvents const& cols,
    BCsRun3 const& /*bcs*/,
    aod::Zdcs const& /*zdcs*/)
  {
    for (const auto& collision : cols) {
      const auto& foundBC = collision.foundBC_as<BCsRun3>();
      if (foundBC.has_zdc()) {
        const auto& zdc = foundBC.zdc();
        auto znA = zdc.amplitudeZNA() / cfgCollisionEnergy;
        auto znC = zdc.amplitudeZNC() / cfgCollisionEnergy;
        auto zpA = zdc.amplitudeZPA() / cfgCollisionEnergy;
        auto zpC = zdc.amplitudeZPC() / cfgCollisionEnergy;
        float sumZNC = ((zdc.energySectorZNC())[0] + (zdc.energySectorZNC())[1] + (zdc.energySectorZNC())[2] + (zdc.energySectorZNC())[3]) / cfgCollisionEnergy;
        float sumZNA = ((zdc.energySectorZNA())[0] + (zdc.energySectorZNA())[1] + (zdc.energySectorZNA())[2] + (zdc.energySectorZNA())[3]) / cfgCollisionEnergy;
        float sumZPC = ((zdc.energySectorZPC())[0] + (zdc.energySectorZPC())[1] + (zdc.energySectorZPC())[2] + (zdc.energySectorZPC())[3]) / cfgCollisionEnergy;
        float sumZPA = ((zdc.energySectorZPA())[0] + (zdc.energySectorZPA())[1] + (zdc.energySectorZPA())[2] + (zdc.energySectorZPA())[3]) / cfgCollisionEnergy;
        float commonSumZnc = zdc.energyCommonZNC() / cfgCollisionEnergy;
        float commonSumZna = zdc.energyCommonZNA() / cfgCollisionEnergy;
        float commonSumZpc = zdc.energyCommonZPC() / cfgCollisionEnergy;
        float commonSumZpa = zdc.energyCommonZPA() / cfgCollisionEnergy;
        float aZEM1 = zdc.amplitudeZEM1();
        float aZEM2 = zdc.amplitudeZEM2();
        float sumZEMs = aZEM1 + aZEM2;
        auto tZNA = zdc.timeZNA();
        auto tZNC = zdc.timeZNC();
        auto tZPA = zdc.timeZPA();
        auto tZPC = zdc.timeZPC();
        if (isTDCcut) {
          if ((tZNA >= minTdcZn) && (tZNA <= maxTdcZn))
            histos.fill(HIST("ampZna"), znA);
          if ((tZNC >= minTdcZn) && (tZNC <= minTdcZn))
            histos.fill(HIST("ampZnc"), znC);
          if ((tZPA >= minTdcZp) && (tZPA <= maxTdcZp))
            histos.fill(HIST("ampZpa"), zpA);
          if ((tZPC >= minTdcZp) && (tZPC <= maxTdcZp))
            histos.fill(HIST("ampZpc"), zpC);
          if (((tZNC >= minTdcZn) && (tZNC <= maxTdcZn)) && ((tZNA >= minTdcZn) && (tZNA <= maxTdcZn)))
            histos.fill(HIST("ZnVsZem"), sumZEMs, znC + znA);
          if (((tZNC >= minTdcZn) && (tZNC <= maxTdcZn)) && ((tZNA >= minTdcZn) && (tZNA <= maxTdcZn)))
            histos.fill(HIST("ZnaVsZnc"), znA, znC);
          if (((tZPC >= minTdcZp) && (tZPC <= maxTdcZp)) && ((tZPA >= minTdcZp) && (tZPA <= maxTdcZp)))
            histos.fill(HIST("ZpaVsZpc"), zpA, zpC);
          if ((tZNA >= minTdcZn) && (tZNA <= maxTdcZn))
            histos.fill(HIST("ZnaVsZpa"), znA, zpA);
          if ((tZNC >= minTdcZn) && (tZNC <= maxTdcZn))
            histos.fill(HIST("ZncVsZpc"), znC, zpC);
        } else {
          histos.fill(HIST("ampZna"), znA);
          histos.fill(HIST("ampZnc"), znC);
          histos.fill(HIST("ampZpa"), zpA);
          histos.fill(HIST("ampZpc"), zpC);
          histos.fill(HIST("ZnVsZem"), sumZEMs, znC + znA);
          histos.fill(HIST("ZnaVsZnc"), znA, znC);
          histos.fill(HIST("ZpaVsZpc"), zpA, zpC);
          histos.fill(HIST("ZnaVsZpa"), znA, zpA);
          histos.fill(HIST("ZncVsZpc"), znC, zpC);
        }
        histos.fill(HIST("ampZEM1"), aZEM1);
        histos.fill(HIST("ampZEM2"), aZEM2);
        histos.fill(HIST("ZnccVsZncSum"), sumZNC, commonSumZnc);
        histos.fill(HIST("ZnacVsZnaSum"), sumZNA, commonSumZna);
        histos.fill(HIST("ZpccVsZpcSum"), sumZPC, commonSumZpc);
        histos.fill(HIST("ZpacVsZpaSum"), sumZPA, commonSumZpa);
        histos.fill(HIST("ZncVsTdc"), zdc.timeZNC(), znC);
        histos.fill(HIST("ZnaVsTdc"), zdc.timeZNA(), znA);
        histos.fill(HIST("ZpcVsTdc"), zdc.timeZPC(), zpC);
        histos.fill(HIST("ZpaVsTdc"), zdc.timeZPA(), zpA);
        histos.fill(HIST("Zem1VsTdc"), zdc.timeZEM1(), aZEM1);
        histos.fill(HIST("Zem2VsTdc"), zdc.timeZEM2(), aZEM2);
      }
    }
  }

  PROCESS_SWITCH(FlowZdcTask, processQA, "Process QA", true);
  PROCESS_SWITCH(FlowZdcTask, processZdcCollAssoc, "Processing ZDC w. collision association", false);
  PROCESS_SWITCH(FlowZdcTask, processZdc, "Process ZDC without corrections or associations", true);

}; // end of struct function

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowZdcTask>(cfgc)};
}
