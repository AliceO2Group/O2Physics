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

/// \file   bjetTaggingGnn.cxx
/// \brief  b-jet tagging using GNN
///
/// \author Changhwan Choi <changhwan.choi@cern.ch>, Pusan National University

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetTaggingUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetSubtraction.h"
#include "PWGJE/DataModel/JetTagging.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/DataModel/EventSelection.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/runDataProcessing.h>

#include <TAxis.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace bjet_tagging_gnn_evtsel
{
enum class EvtSelFlag : uint8_t {
  kNone = 0,
  kINEL = 1 << 0,
  kColl = 1 << 1,
  kTVX = 1 << 2,
  kNoTFB = 1 << 3,
  kNoITSROFB = 1 << 4,
  kNoPileup = 1 << 5,
  kIsGoodZvtx = 1 << 6,
  kZvtx = 1 << 7,
  // kINELgt0 = 1 << 8,
  // kINELgt0rec = 1 << 9,

  INEL = kINEL,
  INELZvtx = kINEL | kZvtx,
  Coll = kINEL | kColl,
  CollZvtx = kINEL | kColl | kZvtx,
  TVX = kINEL | kColl | kTVX,
  TVXZvtx = kINEL | kColl | kTVX | kZvtx,
  SelMC = kINEL | kColl | kTVX | kNoTFB,
  SelMCZvtx = kINEL | kColl | kTVX | kNoTFB | kZvtx,
  Sel8 = kINEL | kColl | kTVX | kNoTFB | kNoITSROFB,
  Sel8Zvtx = kINEL | kColl | kTVX | kNoTFB | kNoITSROFB | kZvtx,
  Sel8Full = kINEL | kColl | kTVX | kNoTFB | kNoITSROFB | kNoPileup,
  Sel8FullZvtx = kINEL | kColl | kTVX | kNoTFB | kNoITSROFB | kNoPileup | kZvtx,
  Sel8FullGood = kINEL | kColl | kTVX | kNoTFB | kNoITSROFB | kNoPileup | kIsGoodZvtx,
  Sel8FullGoodZvtx = kINEL | kColl | kTVX | kNoTFB | kNoITSROFB | kNoPileup | kIsGoodZvtx | kZvtx,
  // INELgt0 = kINEL | kZvtx | kINELgt0
  // INELgt0rec = kINEL | kZvtx | kColl | kTVX | kNoTFB | kNoITSROFB | kINELgt0rec
};
constexpr EvtSelFlag operator|(EvtSelFlag a, EvtSelFlag b)
{
  return static_cast<EvtSelFlag>(
    static_cast<uint8_t>(a) | static_cast<uint8_t>(b));
}
constexpr EvtSelFlag operator|=(EvtSelFlag& a, EvtSelFlag b)
{
  return a = a | b;
}
constexpr EvtSelFlag operator&(EvtSelFlag a, EvtSelFlag b)
{
  return static_cast<EvtSelFlag>(
    static_cast<uint8_t>(a) & static_cast<uint8_t>(b));
}
constexpr bool hasAll(EvtSelFlag value, EvtSelFlag required)
{
  return (value & required) == required;
}
constexpr bool hasAny(EvtSelFlag value, EvtSelFlag mask)
{
  return (value & mask) != EvtSelFlag::kNone;
}
enum class EvtSel {
  None,
  INEL,
  INELZvtx,
  Coll,
  CollZvtx,
  TVX,
  TVXZvtx,
  SelMC,
  SelMCZvtx,
  Sel8,
  Sel8Zvtx,
  Sel8Full,
  Sel8FullZvtx,
  Sel8FullGood,
  Sel8FullGoodZvtx
  // INELgt0,
  // INELgt0rec
};
// Ordered list of event-selection stages shared by all "evtsel"-differential histograms below;
// bin index in the histogram == index in this table (all "evtsel" histograms span the full table,
// so a given event-selection stage always lands on the same bin at particle- and reco-level, matching
// hCollCounter/hMcCollCounter). Reco-level histograms only ever get bins filled for stages
// [kEvtSelStageRecoFirst, kEvtSelStageRecoLast]; bins before kEvtSelStageRecoFirst are unfillable
// placeholders (see setEvtSelStageAxisLabels()), mirroring hCollCounter's "_1"/"_2" bins.
// MC-truth (MCP) ones fill the full range.
struct EvtSelStage {
  EvtSelFlag flag;
  const char* label;
};
constexpr std::array<EvtSelStage, 14> kEvtSelStages{{{.flag = EvtSelFlag::INEL, .label = "INEL"},
                                                     {.flag = EvtSelFlag::INELZvtx, .label = "INEL+Zvtx"},
                                                     {.flag = EvtSelFlag::Coll, .label = "Coll"},
                                                     {.flag = EvtSelFlag::CollZvtx, .label = "Coll+Zvtx"},
                                                     {.flag = EvtSelFlag::TVX, .label = "TVX"},
                                                     {.flag = EvtSelFlag::TVXZvtx, .label = "TVX+Zvtx"},
                                                     {.flag = EvtSelFlag::SelMC, .label = "SelMC"},
                                                     {.flag = EvtSelFlag::SelMCZvtx, .label = "SelMC+Zvtx"},
                                                     {.flag = EvtSelFlag::Sel8, .label = "Sel8"},
                                                     {.flag = EvtSelFlag::Sel8Zvtx, .label = "Sel8+Zvtx"},
                                                     {.flag = EvtSelFlag::Sel8Full, .label = "Sel8Full"},
                                                     {.flag = EvtSelFlag::Sel8FullZvtx, .label = "Sel8Full+Zvtx"},
                                                     {.flag = EvtSelFlag::Sel8FullGood, .label = "Sel8FullGood"},
                                                     {.flag = EvtSelFlag::Sel8FullGoodZvtx, .label = "Sel8FullGood+Zvtx"}}};
constexpr int kEvtSelStageRecoFirst = 2; // Coll
constexpr int kEvtSelStageRecoLast = 13; // Sel8FullGood+Zvtx

// Labels a full-range [0, kEvtSelStages.size()-1] evtsel axis. Bins below `first` (unreachable at
// reco level) get numbered placeholder labels ("_1", "_2", ...) instead of a stage label, matching
// hCollCounter's placeholder bins for the same reason: keeps every "evtsel" histogram's bin index
// aligned to the same event-selection stage across particle- and reco-level histograms.
void setEvtSelStageAxisLabels(TAxis* axis, int first = 0, int last = kEvtSelStages.size() - 1)
{
  for (int i = 0; i < first; ++i) {
    axis->SetBinLabel(i + 1, ("_" + std::to_string(i + 1)).c_str());
  }
  for (int i = first; i <= last; ++i) {
    axis->SetBinLabel(i + 1, kEvtSelStages[i].label);
  }
}
}; // namespace bjet_tagging_gnn_evtsel
using namespace bjet_tagging_gnn_evtsel;

struct BjetTaggingGnn {

  HistogramRegistry registry;

  // event level configurables
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<bool> useEventWeight{"useEventWeight", true, "Flag whether to scale histograms with the event weight"};

  // track level configurables
  Configurable<std::string> trackSelections{"trackSelections", "QualityTracks", "set track selections"};
  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum track pT"};
  Configurable<float> trackPtMax{"trackPtMax", 1000.0, "maximum track pT"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum track eta"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum track eta"};
  Configurable<float> trackPtMinGnn{"trackPtMinGnn", 0.15, "minimum track pT for GNN inputs"};

  Configurable<float> maxIPxy{"maxIPxy", 10, "maximum track DCA in xy plane"};
  Configurable<float> maxIPz{"maxIPz", 10, "maximum track DCA in z direction"};

  Configurable<float> trackNppCrit{"trackNppCrit", 0.95, "track not physical primary ratio"};

  // jet level configurables
  Configurable<float> jetPtMin{"jetPtMin", 0.0, "minimum jet pT"};
  Configurable<float> jetPtMax{"jetPtMax", 1000.0, "maximum jet pT"};
  Configurable<float> jetEtaMin{"jetEtaMin", -0.9, "minimum jet pseudorapidity"};
  Configurable<float> jetEtaMax{"jetEtaMax", 0.9, "maximum jet pseudorapidity"};
  Configurable<float> leadingConstituentPtMin{"leadingConstituentPtMin", -99.0, "minimum pT selection on jet constituent"};
  Configurable<float> leadingConstituentPtMax{"leadingConstituentPtMax", 9999.0, "maximum pT selection on jet constituent"};
  Configurable<float> jetAreaFractionMin{"jetAreaFractionMin", -99.0, "used to make a cut on the jet areas"};

  Configurable<std::vector<double>> jetRadii{"jetRadii", std::vector<double>{0.4}, "jet resolution parameters"};

  Configurable<double> dbMin{"dbMin", -15., "minimum GNN Db"};
  Configurable<double> dbMax{"dbMax", 15., "maximum GNN Db"};
  Configurable<int> dbNbins{"dbNbins", 3000, "number of bins in axisDbFine"};

  Configurable<bool> doDataDriven{"doDataDriven", false, "Flag whether to use fill THnSpase for data driven methods"};
  Configurable<bool> doDataDrivenExtra{"doDataDrivenExtra", false, "Flag whether to add extra axes to THnSparses"};
  Configurable<bool> doDataDrivenSV{"doDataDrivenSV", false, "Flag whether to subscribe SV tables and to use fill THnSparse for data driven methods for SV"};
  Configurable<bool> callSumw2{"callSumw2", false, "Flag whether to call THnSparse::Sumw2() for error calculation"};

  Configurable<int> trainingDatasetRatioParam{"trainingDatasetRatioParam", 0, "Parameter for splitting training/evaluation datasets by collisionId"};

  // Software Trigger configurables
  Configurable<bool> doSoftwareTriggerSelection{"doSoftwareTriggerSelection", false, "set to true when using triggered datasets"};
  Configurable<std::string> triggerMasks{"triggerMasks", "fJetChLowPt", "possible JE Trigger masks: fJetChLowPt,fJetChHighPt,fEMCALReadout,fJetFullHighPt,fJetFullLowPt"};

  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  // Service
  Service<o2::framework::O2DatabasePDG> pdg{};
  Service<o2::ccdb::BasicCCDBManager> ccdb{};

  // Event selection bits
  std::vector<int> eventSelectionBits;
  std::vector<int> eventSelectionBitsTVX;
  std::vector<int> eventSelectionBitsSelMC;
  std::vector<int> eventSelectionBitsSel8;
  std::vector<int> eventSelectionBitsSel8Full;
  std::vector<int> eventSelectionBitsSel8FullGood;

  int trackSelectionBits{};

  std::vector<double> jetRadiiValues;

  void init(InitContext const&)
  {
    jetRadiiValues = (std::vector<double>)jetRadii;

    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    eventSelectionBitsTVX = jetderiveddatautilities::initialiseEventSelectionBits("TVX");
    eventSelectionBitsSel8 = jetderiveddatautilities::initialiseEventSelectionBits("sel8");
    eventSelectionBitsSelMC = jetderiveddatautilities::initialiseEventSelectionBits("selMC");
    eventSelectionBitsSel8Full = jetderiveddatautilities::initialiseEventSelectionBits("sel8Full");
    eventSelectionBitsSel8FullGood = jetderiveddatautilities::initialiseEventSelectionBits("sel8Full+IsGoodZvtxFT0vsPV");

    trackSelectionBits = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

    if (doprocessDataJetsTrig) {
      zorroSummary.setObject(zorro.getZorroSummary());
    }

    registry.add("h_event_counter", ";analysis collision counter", {HistType::kTH1F, {{1, 0.0, 1.0}}}, callSumw2);
    registry.add("h_event_counter_mcp", ";analysis collision matched MC collision counter", {HistType::kTH1F, {{1, 0.0, 1.0}}}, callSumw2);
    registry.add("h_vertexZ", "Vertex Z;#it{Z} (cm)", {HistType::kTH1F, {{100, -20.0, 20.0}}}, callSumw2);
    registry.add("hCollCounter", ";collision counter", {HistType::kTH1F, {{14, 1.0, 15.0}}}, callSumw2);
    auto hCollCounter = registry.get<TH1>(HIST("hCollCounter"));
    hCollCounter->GetXaxis()->SetBinLabel(1, "_1");
    hCollCounter->GetXaxis()->SetBinLabel(2, "_2");
    hCollCounter->GetXaxis()->SetBinLabel(3, "Coll");
    hCollCounter->GetXaxis()->SetBinLabel(4, "Coll+Zvtx");
    hCollCounter->GetXaxis()->SetBinLabel(5, "Coll+TVX");
    hCollCounter->GetXaxis()->SetBinLabel(6, "Coll+TVX+Zvtx");
    hCollCounter->GetXaxis()->SetBinLabel(7, "Coll+TVX+NoTFB"); // selMC
    hCollCounter->GetXaxis()->SetBinLabel(8, "Coll+TVX+NoTFB+Zvtx");
    hCollCounter->GetXaxis()->SetBinLabel(9, "Coll+TVX+NoTFB+NoITSROFB"); // sel8
    hCollCounter->GetXaxis()->SetBinLabel(10, "Coll+TVX+NoTFB+NoITSROFB+Zvtx");
    hCollCounter->GetXaxis()->SetBinLabel(11, "Coll+TVX+NoTFB+NoITSROFB+NoPileup"); // sel8Full
    hCollCounter->GetXaxis()->SetBinLabel(12, "Coll+TVX+NoTFB+NoITSROFB+NoPileup+Zvtx");
    hCollCounter->GetXaxis()->SetBinLabel(13, "Coll+TVX+NoTFB+NoITSROFB+NoPileup+GoodZvtx"); // sel8FullGood
    hCollCounter->GetXaxis()->SetBinLabel(14, "Coll+TVX+NoTFB+NoITSROFB+NoPileup+GoodZvtx+Zvtx");
    // hCollCounter->GetXaxis()->SetBinLabel(11, "_11");
    // hCollCounter->GetXaxis()->SetBinLabel(12, "INELgt0+Zvtx(rec)"); // sel8
    registry.add("hMcCollCounter", ";MC collision counter", {HistType::kTH1F, {{14, 1.0, 15.0}}}, callSumw2);
    auto hMcCollCounter = registry.get<TH1>(HIST("hMcCollCounter"));
    hMcCollCounter->GetXaxis()->SetBinLabel(1, "McColl(INEL)");
    hMcCollCounter->GetXaxis()->SetBinLabel(2, "McColl+Zvtx");
    hMcCollCounter->GetXaxis()->SetBinLabel(3, "McColl(-> Coll)");
    hMcCollCounter->GetXaxis()->SetBinLabel(4, "McColl+Zvtx(-> Coll+Zvtx)");
    hMcCollCounter->GetXaxis()->SetBinLabel(5, "McColl(-> Coll+TVX)");
    hMcCollCounter->GetXaxis()->SetBinLabel(6, "McColl+Zvtx(-> Coll+TVX+Zvtx)");
    hMcCollCounter->GetXaxis()->SetBinLabel(7, "McColl(-> Coll+TVX+NoTFB)"); // selMC
    hMcCollCounter->GetXaxis()->SetBinLabel(8, "McColl+Zvtx(-> Coll+TVX+NoTFB+Zvtx)");
    hMcCollCounter->GetXaxis()->SetBinLabel(9, "McColl(-> Coll+TVX+NoTFB+NoITSROFB)"); // sel8
    hMcCollCounter->GetXaxis()->SetBinLabel(10, "McColl+Zvtx(-> Coll+TVX+NoTFB+NoITSROFB+Zvtx)");
    hMcCollCounter->GetXaxis()->SetBinLabel(11, "McColl(-> Coll+TVX+NoTFB+NoITSROFB+NoPileup)"); // sel8Full
    hMcCollCounter->GetXaxis()->SetBinLabel(12, "McColl+Zvtx(-> Coll+TVX+NoTFB+NoITSROFB+NoPileup+Zvtx)");
    hMcCollCounter->GetXaxis()->SetBinLabel(13, "McColl(-> Coll+TVX+NoTFB+NoITSROFB+NoPileup+GoodZvtx)"); // sel8FullGood
    hMcCollCounter->GetXaxis()->SetBinLabel(14, "McColl+Zvtx(-> Coll+TVX+NoTFB+NoITSROFB+NoPileup+GoodZvtx+Zvtx)");
    // hMcCollCounter->GetXaxis()->SetBinLabel(11, "INELgt0+Zvtx");
    // hMcCollCounter->GetXaxis()->SetBinLabel(12, "INELgt0+Zvtx(-> INELgt0+Zvtx(rec))"); // sel8
    // Same 14-bin layout as hCollCounter/hMcCollCounter (see fillMcCollCounterBC()), but the ladder is
    // evaluated on the MC-truth BC's raw evsel bits instead of a matched reco collision - bins 1-2
    // (INEL/INEL+Zvtx) are placeholders here too since that's already covered by hMcCollCounter.
    registry.add("hMcCollCounterBC", ";MC collision counter (BC-based)", {HistType::kTH1F, {{14, 1.0, 15.0}}}, callSumw2);
    auto hMcCollCounterBC = registry.get<TH1>(HIST("hMcCollCounterBC"));
    hMcCollCounterBC->GetXaxis()->SetBinLabel(1, "_1");
    hMcCollCounterBC->GetXaxis()->SetBinLabel(2, "_2");
    hMcCollCounterBC->GetXaxis()->SetBinLabel(3, "McColl(-> BC)");
    hMcCollCounterBC->GetXaxis()->SetBinLabel(4, "McColl+Zvtx(-> BC)");
    hMcCollCounterBC->GetXaxis()->SetBinLabel(5, "McColl(-> BC+TVX)");
    hMcCollCounterBC->GetXaxis()->SetBinLabel(6, "McColl+Zvtx(-> BC+TVX)");
    hMcCollCounterBC->GetXaxis()->SetBinLabel(7, "McColl(-> BC+TVX+NoTFB)"); // selMC
    hMcCollCounterBC->GetXaxis()->SetBinLabel(8, "McColl+Zvtx(-> BC+TVX+NoTFB)");
    hMcCollCounterBC->GetXaxis()->SetBinLabel(9, "McColl(-> BC+TVX+NoTFB+NoITSROFB)"); // sel8
    hMcCollCounterBC->GetXaxis()->SetBinLabel(10, "McColl+Zvtx(-> BC+TVX+NoTFB+NoITSROFB)");
    hMcCollCounterBC->GetXaxis()->SetBinLabel(11, "McColl(-> BC+TVX+NoTFB+NoITSROFB+NoPileup)"); // sel8Full
    hMcCollCounterBC->GetXaxis()->SetBinLabel(12, "McColl+Zvtx(-> BC+TVX+NoTFB+NoITSROFB+NoPileup)");
    hMcCollCounterBC->GetXaxis()->SetBinLabel(13, "McColl(-> BC+TVX+NoTFB+NoITSROFB+NoPileup+GoodZvtx)"); // sel8FullGood
    hMcCollCounterBC->GetXaxis()->SetBinLabel(14, "McColl+Zvtx(-> BC+TVX+NoTFB+NoITSROFB+NoPileup+GoodZvtx)");
    registry.add("hBCCounter", "", {HistType::kTH1F, {{10, 0.0, 10.}}}, callSumw2);
    auto hBCCounter = registry.get<TH1>(HIST("hBCCounter"));
    hBCCounter->GetXaxis()->SetBinLabel(1, "BC");
    hBCCounter->GetXaxis()->SetBinLabel(2, "BC+TVX");
    hBCCounter->GetXaxis()->SetBinLabel(3, "BC+TVX+NoTFB");
    hBCCounter->GetXaxis()->SetBinLabel(4, "BC+TVX+NoTFB+NoITSROFB");
    hBCCounter->GetXaxis()->SetBinLabel(5, "CollinBC");
    hBCCounter->GetXaxis()->SetBinLabel(6, "CollinBC+Sel8");
    hBCCounter->GetXaxis()->SetBinLabel(7, "CollinBC+Sel8+VtxZ");
    hBCCounter->GetXaxis()->SetBinLabel(8, "CollinBC+Sel8Full");
    hBCCounter->GetXaxis()->SetBinLabel(9, "CollinBC+Sel8Full+GoodZvtx");
    hBCCounter->GetXaxis()->SetBinLabel(10, "CollinBC+Sel8Full+VtxZ+GoodZvtx");

    const AxisSpec axisTrackpT{200, 0., 200., "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisTrackpTFine{1000, 0., 10., "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisJetpT{250, 0., 250., "#it{p}_{T, ch jet} (GeV/#it{c})"};
    // Used in place of axisJetpT for every "_sub"-suffixed histogram (UE-subtracted jet pT can go negative).
    const AxisSpec axisJetpTSub{300, -50., 250., "#it{p}_{T, ch jet}^{sub} (GeV/#it{c})"};
    const AxisSpec axisJetEta{200, -0.8, 0.8, "#it{#eta}_{jet}"};
    const AxisSpec axisDb{200, dbMin, dbMax, "#it{D}_{b}"};
    const AxisSpec axisDbFine{dbNbins, dbMin, dbMax, "#it{D}_{b}"};
    const AxisSpec axisJetMass{200, 0., 50., "#it{m}_{jet} (GeV/#it{c}^{2})"};
    const AxisSpec axisJetProb{200, 0., 40., "-ln(JP)"};
    const AxisSpec axisNTracks{42, 0, 42, "#it{N}_{tracks}"};
    // Jet-flavour category axis shared by every "_flavor"-suffixed histogram below (see getJetFlavorCat()/setJetFlavorCatAxisLabels()).
    const AxisSpec axisJetFlavorCat{3, -0.5, 2.5, "flavour"};
    const AxisSpec axisLfMatchStatus{2, -0.5, 1.5, "lf match status"};
    // Event-selection-stage axis shared by every "_evtsel"-suffixed histogram below (reco and particle-level
    // alike span the full table so a given stage always lands on the same bin, see kEvtSelStages/fillEvtSelStages()).
    const AxisSpec axisEvtSelStage{static_cast<int>(kEvtSelStages.size()), -0.5, static_cast<double>(kEvtSelStages.size()) - 0.5, "evt. sel. stage"};
    const AxisSpec axisSVMass{200, 0., 10., "#it{m}_{SV} (GeV/#it{c}^{2})"};
    const AxisSpec axisSVLxyS{200, 0., 100., "L_{xy}/#sigma"};
    const AxisSpec axisSVDispersion{100, 0., 0.5, "SV dispersion"};

    // `suffix` ("" or "_sub") is appended to jet-pT-dependent histogram names below. Histograms whose value
    // doesn't depend on jet pT (h_SVMass, h_jetPhi, ...) are booked once, unsuffixed, outside these lambdas.
    if (doprocessDataJetsSV || doprocessMCDJetsSV || doprocessDataJetsSVSub || doprocessMCDJetsSVSub) {
      registry.add("h_SVMass", "", {HistType::kTH1F, {axisSVMass}}, callSumw2);
    }
    auto addSVHistograms = [&](const AxisSpec& axisJetpT, const std::string& suffix) {
      // Best SV per jet (highest decay length significance); see fillDataJetHistogramsSV()/fillMCDJetHistogramsSV().
      registry.add("h2_SVMass_jetpT" + suffix, "", HistType::kTH2F, {axisJetpT, axisSVMass}, callSumw2);
      registry.add("h2_SVLxyS_jetpT" + suffix, "", HistType::kTH2F, {axisJetpT, axisSVLxyS}, callSumw2);
      registry.add("h2_SVDispersion_jetpT" + suffix, "", HistType::kTH2F, {axisJetpT, axisSVDispersion}, callSumw2);
    };
    if (doprocessMCDJetsSV || doprocessMCDJetsSVSub) {
      registry.add("h2_SVMass_flavor", "", {HistType::kTH2F, {axisSVMass, axisJetFlavorCat}}, callSumw2);
    }
    // Flavour split into three TH2 (b/c/lf) instead of one TH3 with a flavour axis: GRID-scale hadd
    // merging of TH3 has been observed to dump all content into the flavour-axis overflow bin, while the
    // equivalent TH2 merges correctly (see fillSVHistograms()/fillMCDJetHistograms() fill sites).
    auto addSVFlavorHistograms = [&](const AxisSpec& axisJetpT, const std::string& suffix) {
      registry.add("h2_SVMass_jetpT_b" + suffix, "", HistType::kTH2F, {axisJetpT, axisSVMass}, callSumw2);
      registry.add("h2_SVMass_jetpT_c" + suffix, "", HistType::kTH2F, {axisJetpT, axisSVMass}, callSumw2);
      registry.add("h2_SVMass_jetpT_lf" + suffix, "", HistType::kTH2F, {axisJetpT, axisSVMass}, callSumw2);
      registry.add("h2_SVLxyS_jetpT_b" + suffix, "", HistType::kTH2F, {axisJetpT, axisSVLxyS}, callSumw2);
      registry.add("h2_SVLxyS_jetpT_c" + suffix, "", HistType::kTH2F, {axisJetpT, axisSVLxyS}, callSumw2);
      registry.add("h2_SVLxyS_jetpT_lf" + suffix, "", HistType::kTH2F, {axisJetpT, axisSVLxyS}, callSumw2);
      registry.add("h2_SVDispersion_jetpT_b" + suffix, "", HistType::kTH2F, {axisJetpT, axisSVDispersion}, callSumw2);
      registry.add("h2_SVDispersion_jetpT_c" + suffix, "", HistType::kTH2F, {axisJetpT, axisSVDispersion}, callSumw2);
      registry.add("h2_SVDispersion_jetpT_lf" + suffix, "", HistType::kTH2F, {axisJetpT, axisSVDispersion}, callSumw2);
    };
    auto addCoreJetHistograms = [&](const AxisSpec& axisJetpT, const std::string& suffix) {
      registry.add("h_jetpT" + suffix, "", HistType::kTH1F, {axisJetpT}, callSumw2);
      registry.add("h2_jetpT_Db" + suffix, "", HistType::kTH2F, {axisJetpT, axisDb});
    };

    if (doprocessDataJetsSV || doprocessMCDJetsSV) {
      addSVHistograms(axisJetpT, "");
    }
    if (doprocessDataJetsSVSub || doprocessMCDJetsSVSub) {
      addSVHistograms(axisJetpTSub, "_sub");
    }
    if (doprocessMCDJetsSV) {
      addSVFlavorHistograms(axisJetpT, "");
    }
    if (doprocessMCDJetsSVSub) {
      addSVFlavorHistograms(axisJetpTSub, "_sub");
    }

    // h_jetEta/h_jetPhi/h_jetMass/h_Db/h2_nTracks_Db don't depend on jet pT - one shared copy regardless of
    // `withSub` (see addCoreJetHistograms()'s comment above).
    registry.add("h_jetEta", "", {HistType::kTH1F, {axisJetEta}}, callSumw2);
    registry.add("h_jetPhi", "", {HistType::kTH1F, {{200, 0., o2::constants::math::TwoPI, "#it{phi}_{jet}"}}});
    registry.add("h_jetMass", "", {HistType::kTH1F, {axisJetMass}});
    registry.add("h_Db", "", {HistType::kTH1F, {axisDbFine}});
    registry.add("h2_nTracks_Db", "", {HistType::kTH2F, {axisNTracks, axisDb}});

    addCoreJetHistograms(axisJetpT, "");
    if (doprocessDataJetsSub || doprocessDataJetsSVSub || doprocessMCDJetsSub || doprocessMCDJetsSVSub || doprocessDataJetsTrigSub) {
      addCoreJetHistograms(axisJetpTSub, "_sub");
    }

    registry.add("h_gnnfeat_trackpT", "", {HistType::kTH1F, {{200, 0., 100., "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("h_gnnfeat_trackPhi", "", {HistType::kTH1F, {{200, 0., o2::constants::math::TwoPI, "#it{#phi}"}}});
    registry.add("h_gnnfeat_trackEta", "", {HistType::kTH1F, {{200, -0.9, 0.9, "#it{#eta}"}}});
    registry.add("h_gnnfeat_trackCharge", "", {HistType::kTH1F, {{3, -1., 2., "#it{q}"}}});
    registry.add("h_gnnfeat_trackDCAxy", "", {HistType::kTH1F, {{200, -5., 5., "DCA_{#it{xy}} (cm)"}}});
    registry.add("h_gnnfeat_trackSigmaDCAxy", "", {HistType::kTH1F, {{200, 0., 5., "#it{#sigma}_{{DCA_{#it{xy}}} (cm)"}}});
    registry.add("h_gnnfeat_trackDCAz", "", {HistType::kTH1F, {{200, -5., 5., "DCA_{#it{z}} (cm)"}}});
    registry.add("h_gnnfeat_trackSigmaDCAz", "", {HistType::kTH1F, {{200, 0., 5., "#it{#sigma}_{{DCA_{#it{z}}} (cm)"}}});
    registry.add("h_gnnfeat_trackITSChi2NCl", "", {HistType::kTH1F, {{200, 0., 40., "ITS #it{#chi}^{2}/ndf"}}});
    registry.add("h_gnnfeat_trackTPCChi2NCl", "", {HistType::kTH1F, {{200, 0., 5., "TPC #it{#chi}^{2}/ndf"}}});
    registry.add("h_gnnfeat_trackITSNCls", "", {HistType::kTH1F, {{8, 0., 8., "ITS NCls"}}});
    registry.add("h_gnnfeat_trackTPCNCls", "", {HistType::kTH1F, {{153, 0., 153., "TPC NCls"}}});
    registry.add("h_gnnfeat_trackTPCNCrossedRows", "", {HistType::kTH1F, {{153, 0., 153., "TPC NCrossedRows"}}});
    registry.add("h_gnnfeat_tracksIPxy", "", {HistType::kTH1F, {{200, -5., 5., "{sIP}_{#it{xy}}"}}});
    registry.add("h_gnnfeat_tracksIPz", "", {HistType::kTH1F, {{200, -5., 5., "{sIP}_{#it{z}}"}}});

    if (doprocessDataTracks || doprocessMCDTracks) {
      registry.add("h_trackpT", "", {HistType::kTH1F, {axisTrackpT}}, callSumw2);
      registry.add("h_tracketa", "", {HistType::kTH1F, {{100, trackEtaMin, trackEtaMax, "#it{#eta}"}}}, callSumw2);
      registry.add("h_trackphi", "", {HistType::kTH1F, {{100, 0.0, o2::constants::math::TwoPI, "#it{#phi}"}}}, callSumw2);
      registry.add("h_dcaXY", "", {HistType::kTH1F, {{200, 0., 4., "|DCA_{#it{xy}}| (cm)"}}}, callSumw2);
      registry.add("h_dcaZ", "", {HistType::kTH1F, {{200, 0., 4., "|DCA_{#it{z}}| (cm)"}}}, callSumw2);
      registry.add("hSparse_dca_pt", "", {HistType::kTHnSparseF, {{1000, 0., 100., "#it{p}_{T} (GeV/#it{c})"}, {200, 0., 4., "|DCA_{#it{z}}| (cm)"}, {200, 0., 4., "|DCA_{#it{xy}}| (cm)"}}}, callSumw2);
    }

    if (doprocessMCDTracks) {
      registry.add("h2_trackpT_partpT", "", {HistType::kTH2F, {axisTrackpT, axisTrackpT}}, callSumw2);
      registry.add("h_partpT_matched_fine", "", {HistType::kTH1F, {axisTrackpTFine}}, callSumw2);
      registry.add("h_partpT", "", {HistType::kTH1F, {axisTrackpT}}, callSumw2);
      registry.add("h_partpT_fine", "", {HistType::kTH1F, {axisTrackpTFine}}, callSumw2);
      // DCA cut study histograms (pT > pTMin)
      registry.add("h_dcaXY_coll_fake", "", {HistType::kTH1F, {{200, 0., 4., "|DCA_{#it{xy}}| (cm)"}}}, callSumw2);       // tracks from fake collisions
      registry.add("h_dcaXY_fake", "", {HistType::kTH1F, {{200, 0., 4., "|DCA_{#it{xy}}| (cm)"}}}, callSumw2);            // fake tracks (no matched particle)
      registry.add("h_dcaXY_coll_matched", "", {HistType::kTH1F, {{200, 0., 4., "|DCA_{#it{xy}}| (cm)"}}}, callSumw2);    // tracks from matched collisions
      registry.add("h_dcaXY_coll_matched_b", "", {HistType::kTH1F, {{200, 0., 4., "|DCA_{#it{xy}}| (cm)"}}}, callSumw2);  // tracks from matched collisions, b hadron decay
      registry.add("h_dcaXY_coll_matched_c", "", {HistType::kTH1F, {{200, 0., 4., "|DCA_{#it{xy}}| (cm)"}}}, callSumw2);  // tracks from matched collisions, c hadron decay
      registry.add("h_dcaXY_coll_matched_lf", "", {HistType::kTH1F, {{200, 0., 4., "|DCA_{#it{xy}}| (cm)"}}}, callSumw2); // tracks from matched collisions, others
      registry.add("h_dcaXY_coll_mismatched", "", {HistType::kTH1F, {{200, 0., 4., "|DCA_{#it{xy}}| (cm)"}}}, callSumw2); // tracks from mismatched collisions
      registry.add("h_dcaXY_npp", "", {HistType::kTH1F, {{200, 0., 4., "|DCA_{#it{xy}}| (cm)"}}}, callSumw2);             // non-physical primary tracks (GenStatusCode=-1) from matched collisions
      registry.add("h_dcaXY_npp_mismatched", "", {HistType::kTH1F, {{200, 0., 4., "|DCA_{#it{xy}}| (cm)"}}}, callSumw2);  // non-physical primary tracks from mismatched collisions
      registry.add("h_dcaZ_coll_fake", "", {HistType::kTH1F, {{200, 0., 4., "|DCA_{#it{z}}| (cm)"}}}, callSumw2);
      registry.add("h_dcaZ_fake", "", {HistType::kTH1F, {{200, 0., 4., "|DCA_{#it{z}}| (cm)"}}}, callSumw2);
      registry.add("h_dcaZ_coll_matched", "", {HistType::kTH1F, {{200, 0., 4., "|DCA_{#it{z}}| (cm)"}}}, callSumw2);
      registry.add("h_dcaZ_coll_matched_b", "", {HistType::kTH1F, {{200, 0., 4., "|DCA_{#it{z}}| (cm)"}}}, callSumw2);
      registry.add("h_dcaZ_coll_matched_c", "", {HistType::kTH1F, {{200, 0., 4., "|DCA_{#it{z}}| (cm)"}}}, callSumw2);
      registry.add("h_dcaZ_coll_matched_lf", "", {HistType::kTH1F, {{200, 0., 4., "|DCA_{#it{z}}| (cm)"}}}, callSumw2);
      registry.add("h_dcaZ_coll_mismatched", "", {HistType::kTH1F, {{200, 0., 4., "|DCA_{#it{z}}| (cm)"}}}, callSumw2);
      registry.add("h_dcaZ_npp", "", {HistType::kTH1F, {{200, 0., 4., "|DCA_{#it{z}}| (cm)"}}}, callSumw2);
      registry.add("h_dcaZ_npp_mismatched", "", {HistType::kTH1F, {{200, 0., 4., "|DCA_{#it{z}}| (cm)"}}}, callSumw2);
      registry.add("hSparse_dca_pt_coll_fake", "", {HistType::kTHnSparseF, {{1000, 0., 100., "#it{p}_{T} (GeV/#it{c})"}, {200, 0., 4., "|DCA_{#it{z}}| (cm)"}, {200, 0., 4., "|DCA_{#it{xy}}| (cm)"}}}, callSumw2);
      registry.add("hSparse_dca_pt_fake", "", {HistType::kTHnSparseF, {{1000, 0., 100., "#it{p}_{T} (GeV/#it{c})"}, {200, 0., 4., "|DCA_{#it{z}}| (cm)"}, {200, 0., 4., "|DCA_{#it{xy}}| (cm)"}}}, callSumw2);
      registry.add("hSparse_dca_pt_coll_matched", "", {HistType::kTHnSparseF, {{1000, 0., 100., "#it{p}_{T} (GeV/#it{c})"}, {200, 0., 4., "|DCA_{#it{z}}| (cm)"}, {200, 0., 4., "|DCA_{#it{xy}}| (cm)"}}}, callSumw2);
      registry.add("hSparse_dca_pt_coll_matched_b", "", {HistType::kTHnSparseF, {{1000, 0., 100., "#it{p}_{T} (GeV/#it{c})"}, {200, 0., 4., "|DCA_{#it{z}}| (cm)"}, {200, 0., 4., "|DCA_{#it{xy}}| (cm)"}}}, callSumw2);
      registry.add("hSparse_dca_pt_coll_matched_c", "", {HistType::kTHnSparseF, {{1000, 0., 100., "#it{p}_{T} (GeV/#it{c})"}, {200, 0., 4., "|DCA_{#it{z}}| (cm)"}, {200, 0., 4., "|DCA_{#it{xy}}| (cm)"}}}, callSumw2);
      registry.add("hSparse_dca_pt_coll_matched_lf", "", {HistType::kTHnSparseF, {{1000, 0., 100., "#it{p}_{T} (GeV/#it{c})"}, {200, 0., 4., "|DCA_{#it{z}}| (cm)"}, {200, 0., 4., "|DCA_{#it{xy}}| (cm)"}}}, callSumw2);
      registry.add("hSparse_dca_pt_coll_mismatched", "", {HistType::kTHnSparseF, {{1000, 0., 100., "#it{p}_{T} (GeV/#it{c})"}, {200, 0., 4., "|DCA_{#it{z}}| (cm)"}, {200, 0., 4., "|DCA_{#it{xy}}| (cm)"}}}, callSumw2);
      registry.add("hSparse_dca_pt_npp", "", {HistType::kTHnSparseF, {{1000, 0., 100., "#it{p}_{T} (GeV/#it{c})"}, {200, 0., 4., "|DCA_{#it{z}}| (cm)"}, {200, 0., 4., "|DCA_{#it{xy}}| (cm)"}}}, callSumw2);
      registry.add("hSparse_dca_pt_npp_mismatched", "", {HistType::kTHnSparseF, {{1000, 0., 100., "#it{p}_{T} (GeV/#it{c})"}, {200, 0., 4., "|DCA_{#it{z}}| (cm)"}, {200, 0., 4., "|DCA_{#it{xy}}| (cm)"}}}, callSumw2);
    }

    if (doprocessDataJetsSel || doprocessMCDJetsSel) {
      auto hEvtsel = registry.add<TH2>("h2_jetpT_evtsel", "", HistType::kTH2F, {axisJetpT, axisEvtSelStage}, callSumw2);
      setEvtSelStageAxisLabels(hEvtsel->GetYaxis(), kEvtSelStageRecoFirst, kEvtSelStageRecoLast);
    }
    if (doprocessDataJetsSelSub || doprocessMCDJetsSelSub) {
      auto hEvtselSub = registry.add<TH2>("h2_jetpT_evtsel_sub", "", HistType::kTH2F, {axisJetpTSub, axisEvtSelStage}, callSumw2);
      setEvtSelStageAxisLabels(hEvtselSub->GetYaxis(), kEvtSelStageRecoFirst, kEvtSelStageRecoLast);
    }

    // b/c/lf categories collapsed onto axisJetFlavorCat instead of one histogram per category. Also includes
    // the geometric matched-response and pTHat-sparse histograms, which get UE-subtracted variants too, not
    // just the core Db/flavor spectra. Histograms whose value doesn't depend on jetpT are booked once, above.
    if (doprocessMCDJets || doprocessMCDJetsSV || doprocessMCDJetsSub || doprocessMCDJetsSVSub) {
      registry.add("h2_Db_flavor", "", {HistType::kTH2F, {axisDbFine, axisJetFlavorCat}});
      registry.add("h2_nTracks_Db_b", "", {HistType::kTH2F, {axisNTracks, axisDb}});
      registry.add("h2_nTracks_Db_c", "", {HistType::kTH2F, {axisNTracks, axisDb}});
      registry.add("h2_nTracks_Db_lf", "", {HistType::kTH2F, {axisNTracks, axisDb}});
      registry.add("h2_Db_lfmatch", "lf-jet", {HistType::kTH2F, {axisDbFine, axisLfMatchStatus}});
      registry.add("h_Db_npp", "NotPhysPrim", {HistType::kTH1F, {axisDbFine}});
      registry.add("h2_Db_npp_flavor", "NotPhysPrim", {HistType::kTH2F, {axisDbFine, axisJetFlavorCat}});
      registry.add("h2_Db_npp_lfmatch", "NotPhysPrim lf-jet", {HistType::kTH2F, {axisDbFine, axisLfMatchStatus}});
      setLfMatchStatusAxisLabels(registry.get<TH2>(HIST("h2_Db_lfmatch"))->GetYaxis());
      setLfMatchStatusAxisLabels(registry.get<TH2>(HIST("h2_Db_npp_lfmatch"))->GetYaxis());
      setJetFlavorCatAxisLabels(registry.get<TH2>(HIST("h2_Db_flavor"))->GetYaxis());
      setJetFlavorCatAxisLabels(registry.get<TH2>(HIST("h2_Db_npp_flavor"))->GetYaxis());
    }
    auto addMCDFlavorAndMatchedHistograms = [&](const AxisSpec& axisJetpT, const std::string& suffix) {
      auto hJetpTFlavor = registry.add<TH2>("h2_jetpT_flavor" + suffix, "", HistType::kTH2F, {axisJetpT, axisJetFlavorCat}, callSumw2);
      registry.add("h2_jetpT_Db_b" + suffix, "", HistType::kTH2F, {axisJetpT, axisDb});
      registry.add("h2_jetpT_Db_c" + suffix, "", HistType::kTH2F, {axisJetpT, axisDb});
      registry.add("h2_jetpT_Db_lf" + suffix, "", HistType::kTH2F, {axisJetpT, axisDb});
      auto hJetpTLfmatch = registry.add<TH2>("h2_jetpT_lfmatch" + suffix, "lf-jet", HistType::kTH2F, {axisJetpT, axisLfMatchStatus}, callSumw2);
      setLfMatchStatusAxisLabels(hJetpTLfmatch->GetYaxis());
      registry.add("h2_jetpT_Db_lfmatched" + suffix, "lf-jet", HistType::kTH2F, {axisJetpT, axisDb}, callSumw2);
      registry.add("h2_jetpT_Db_lfnone" + suffix, "lf-jet", HistType::kTH2F, {axisJetpT, axisDb}, callSumw2);
      auto hResponse = registry.add<TH2>("h2_Response_DetjetpT_PartjetpT" + suffix, "", HistType::kTH2F, {axisJetpT, axisJetpT}, callSumw2);
      registry.add("h2_Response_DetjetpT_PartjetpT_b" + suffix, "", HistType::kTH2F, {axisJetpT, axisJetpT}, callSumw2);
      registry.add("h2_Response_DetjetpT_PartjetpT_c" + suffix, "", HistType::kTH2F, {axisJetpT, axisJetpT}, callSumw2);
      registry.add("h2_Response_DetjetpT_PartjetpT_lf" + suffix, "", HistType::kTH2F, {axisJetpT, axisJetpT}, callSumw2);
      registry.add("h2_jetpT_Db_npp" + suffix, "NotPhysPrim", HistType::kTH2F, {axisJetpT, axisDb});
      registry.add("h2_jetpT_Db_npp_b" + suffix, "NotPhysPrim", HistType::kTH2F, {axisJetpT, axisDb});
      registry.add("h2_jetpT_Db_npp_c" + suffix, "NotPhysPrim", HistType::kTH2F, {axisJetpT, axisDb});
      registry.add("h2_jetpT_Db_npp_lf" + suffix, "NotPhysPrim", HistType::kTH2F, {axisJetpT, axisDb});
      setJetFlavorCatAxisLabels(hJetpTFlavor->GetYaxis());
      (void)hResponse;

      registry.add("h_jetpT_matched" + suffix, "", HistType::kTH1F, {axisJetpT}, callSumw2);
      registry.add("h_jetpT_particle_matched" + suffix, "", HistType::kTH1F, {axisJetpT}, callSumw2);
      auto hJetpTMatchedFlavor = registry.add<TH2>("h2_jetpT_matched_flavor" + suffix, "", HistType::kTH2F, {axisJetpT, axisJetFlavorCat}, callSumw2);
      auto hJetpTParticleMatchedFlavor = registry.add<TH2>("h2_jetpT_particle_matched_flavor" + suffix, "", HistType::kTH2F, {axisJetpT, axisJetFlavorCat}, callSumw2);
      setJetFlavorCatAxisLabels(hJetpTMatchedFlavor->GetYaxis());
      setJetFlavorCatAxisLabels(hJetpTParticleMatchedFlavor->GetYaxis());
      // MCP jets
      registry.add("h_jetpT_particle" + suffix, "", HistType::kTH1F, {axisJetpT}, callSumw2);
      auto hJetpTParticleFlavor = registry.add<TH2>("h2_jetpT_particle_flavor" + suffix, "", HistType::kTH2F, {axisJetpT, axisJetFlavorCat}, callSumw2);
      setJetFlavorCatAxisLabels(hJetpTParticleFlavor->GetYaxis());
      // pTHat study
      registry.add("hSparse_pthat_jetpT" + suffix, "", HistType::kTHnSparseF, {{300, 0., 300., "#hat{#it{p}}_{T} (GeV/#it{c})"}, axisJetpT, axisJetpT}, callSumw2);
      registry.add("hSparse_pthat_jetpT_b" + suffix, "b-jet", HistType::kTHnSparseF, {{300, 0., 300., "#hat{#it{p}}_{T} (GeV/#it{c})"}, axisJetpT, axisJetpT}, callSumw2);
      registry.add("hSparse_pthat_jetpT_c" + suffix, "c-jet", HistType::kTHnSparseF, {{300, 0., 300., "#hat{#it{p}}_{T} (GeV/#it{c})"}, axisJetpT, axisJetpT}, callSumw2);
    };

    if (doprocessMCDJets || doprocessMCDJetsSV) {
      addMCDFlavorAndMatchedHistograms(axisJetpT, "");
    }
    if (doprocessMCDJetsSub || doprocessMCDJetsSVSub) {
      addMCDFlavorAndMatchedHistograms(axisJetpTSub, "_sub");
    }

    if (doprocessMCDJetsSel) {
      auto hEvtselB = registry.add<TH2>("h2_jetpT_evtsel_b", "b-jet", HistType::kTH2F, {axisJetpT, axisEvtSelStage}, callSumw2);
      auto hEvtselC = registry.add<TH2>("h2_jetpT_evtsel_c", "c-jet", HistType::kTH2F, {axisJetpT, axisEvtSelStage}, callSumw2);
      setEvtSelStageAxisLabels(hEvtselB->GetYaxis(), kEvtSelStageRecoFirst, kEvtSelStageRecoLast);
      setEvtSelStageAxisLabels(hEvtselC->GetYaxis(), kEvtSelStageRecoFirst, kEvtSelStageRecoLast);
    }
    if (doprocessMCDJetsSelSub) {
      auto hEvtselBSub = registry.add<TH2>("h2_jetpT_evtsel_b_sub", "b-jet", HistType::kTH2F, {axisJetpTSub, axisEvtSelStage}, callSumw2);
      auto hEvtselCSub = registry.add<TH2>("h2_jetpT_evtsel_c_sub", "c-jet", HistType::kTH2F, {axisJetpTSub, axisEvtSelStage}, callSumw2);
      setEvtSelStageAxisLabels(hEvtselBSub->GetYaxis(), kEvtSelStageRecoFirst, kEvtSelStageRecoLast);
      setEvtSelStageAxisLabels(hEvtselCSub->GetYaxis(), kEvtSelStageRecoFirst, kEvtSelStageRecoLast);
    }

    // Used by both processMCPJets ("", axisJetpT) and processMCPJetsSub ("_sub", axisJetpTSub).
    auto addMCPEvtSelHistograms = [&](const AxisSpec& axisJetpT, const std::string& suffix) {
      // All INEL/Coll/TVX/SelMC/Sel8(+Zvtx) variants collapsed onto axisEvtSelStage.
      auto hEvtsel = registry.add<TH2>("h2_jetpT_particle_evtsel" + suffix, "", HistType::kTH2F, {axisJetpT, axisEvtSelStage}, callSumw2);
      auto hEvtselB = registry.add<TH2>("h2_jetpT_particle_evtsel_b" + suffix, "particle b-jet", HistType::kTH2F, {axisJetpT, axisEvtSelStage}, callSumw2);
      auto hEvtselC = registry.add<TH2>("h2_jetpT_particle_evtsel_c" + suffix, "particle c-jet", HistType::kTH2F, {axisJetpT, axisEvtSelStage}, callSumw2);
      setEvtSelStageAxisLabels(hEvtsel->GetYaxis());
      setEvtSelStageAxisLabels(hEvtselB->GetYaxis());
      setEvtSelStageAxisLabels(hEvtselC->GetYaxis());

      // Same as above, but the event-selection stage is evaluated on the MC-truth BC directly
      // (see fillMcCollCounterBC()) instead of via a matched reco collision.
      auto hEvtselbc = registry.add<TH2>("h2_jetpT_particle_evtselbc" + suffix, "", HistType::kTH2F, {axisJetpT, axisEvtSelStage}, callSumw2);
      auto hEvtselbcB = registry.add<TH2>("h2_jetpT_particle_evtselbc_b" + suffix, "particle b-jet", HistType::kTH2F, {axisJetpT, axisEvtSelStage}, callSumw2);
      auto hEvtselbcC = registry.add<TH2>("h2_jetpT_particle_evtselbc_c" + suffix, "particle c-jet", HistType::kTH2F, {axisJetpT, axisEvtSelStage}, callSumw2);
      setEvtSelStageAxisLabels(hEvtselbc->GetYaxis());
      setEvtSelStageAxisLabels(hEvtselbcB->GetYaxis());
      setEvtSelStageAxisLabels(hEvtselbcC->GetYaxis());
    };

    // processMCPJetsCommon() fills these unconditionally regardless of `withSub` (see fillMcCollCounterBC()/
    // hMcCollCounter's convention) - must be booked whenever either the raw or Sub process is on, not just
    // processMCPJets (processMCPJetsSub alone would otherwise fill an unbooked histogram and throw at runtime).
    if (doprocessMCPJets || doprocessMCPJetsSub) {
      registry.add("h_vertexZ_truth", "Vertex Z truth;#it{Z} (cm)", {HistType::kTH1F, {{100, -20.0, 20.0}}}, callSumw2);
      registry.add("h_vertexZ_truth_coll", "Vertex Z truth (Coll);#it{Z} (cm)", {HistType::kTH1F, {{100, -20.0, 20.0}}}, callSumw2);
      registry.add("h_vertexZ_truth_tvx", "Vertex Z truth (TVX);#it{Z} (cm)", {HistType::kTH1F, {{100, -20.0, 20.0}}}, callSumw2);
      registry.add("h_vertexZ_truth_bctvx", "Vertex Z truth (BC+TVX);#it{Z} (cm)", {HistType::kTH1F, {{100, -20.0, 20.0}}}, callSumw2);
    }
    if (doprocessMCPJets) {
      addMCPEvtSelHistograms(axisJetpT, "");
    }
    if (doprocessMCPJetsSub) {
      addMCPEvtSelHistograms(axisJetpTSub, "_sub");
    }

    auto addDataDrivenHistograms = [&](const AxisSpec& axisJetpT, const std::string& suffix, bool mcdJetsEnabled) {
      // hSparse_Incljets' last axis is the best-SV mass instead of the jet mass when doDataDrivenSV is on
      // (only meaningful together with processDataJetsSV/processMCDJetsSV, see fill*JetHistogramsSV()).
      if (doDataDrivenSV) {
        registry.add("hSparse_Incljets" + suffix, "", HistType::kTHnSparseF, {axisJetpT, axisDbFine, axisNTracks, axisSVMass}, callSumw2);
      } else if (doDataDrivenExtra) {
        registry.add("hSparse_Incljets" + suffix, "", HistType::kTHnSparseF, {axisJetpT, axisDbFine, axisNTracks, axisJetMass}, callSumw2);
      } else {
        registry.add("hSparse_Incljets" + suffix, "", HistType::kTHnSparseF, {axisJetpT, axisDbFine, axisNTracks}, callSumw2);
      }
      if (doDataDrivenSV) {
        if (mcdJetsEnabled) {
          registry.add("hSparse_bjets" + suffix, "", HistType::kTHnSparseF, {axisJetpT, axisDbFine, axisNTracks, axisSVMass}, callSumw2);
          registry.add("hSparse_cjets" + suffix, "", HistType::kTHnSparseF, {axisJetpT, axisDbFine, axisNTracks, axisSVMass}, callSumw2);
          registry.add("hSparse_lfjets" + suffix, "", HistType::kTHnSparseF, {axisJetpT, axisDbFine, axisNTracks, axisSVMass}, callSumw2);
          registry.add("hSparse_lfjets_none" + suffix, "", HistType::kTHnSparseF, {axisJetpT, axisDbFine, axisNTracks, axisSVMass}, callSumw2);
          registry.add("hSparse_lfjets_matched" + suffix, "", HistType::kTHnSparseF, {axisJetpT, axisDbFine, axisNTracks, axisSVMass}, callSumw2);
        }
      } else if (doDataDrivenExtra) {
        if (mcdJetsEnabled) {
          registry.add("hSparse_bjets" + suffix, "", HistType::kTHnSparseF, {axisJetpT, axisDbFine, axisNTracks, axisJetMass}, callSumw2);
          registry.add("hSparse_cjets" + suffix, "", HistType::kTHnSparseF, {axisJetpT, axisDbFine, axisNTracks, axisJetMass}, callSumw2);
          registry.add("hSparse_lfjets" + suffix, "", HistType::kTHnSparseF, {axisJetpT, axisDbFine, axisNTracks, axisJetMass}, callSumw2);
          registry.add("hSparse_lfjets_none" + suffix, "", HistType::kTHnSparseF, {axisJetpT, axisDbFine, axisNTracks, axisJetMass}, callSumw2);
          registry.add("hSparse_lfjets_matched" + suffix, "", HistType::kTHnSparseF, {axisJetpT, axisDbFine, axisNTracks, axisJetMass}, callSumw2);
        }
      } else {
        if (mcdJetsEnabled) {
          registry.add("hSparse_bjets" + suffix, "", HistType::kTHnSparseF, {axisJetpT, axisDbFine, axisNTracks}, callSumw2);
          registry.add("hSparse_cjets" + suffix, "", HistType::kTHnSparseF, {axisJetpT, axisDbFine, axisNTracks}, callSumw2);
          registry.add("hSparse_lfjets" + suffix, "", HistType::kTHnSparseF, {axisJetpT, axisDbFine, axisNTracks}, callSumw2);
          registry.add("hSparse_lfjets_none" + suffix, "", HistType::kTHnSparseF, {axisJetpT, axisDbFine, axisNTracks}, callSumw2);
          registry.add("hSparse_lfjets_matched" + suffix, "", HistType::kTHnSparseF, {axisJetpT, axisDbFine, axisNTracks}, callSumw2);
        }
      }
    };

    if (doDataDriven) {
      addDataDrivenHistograms(axisJetpT, "", doprocessMCDJets || doprocessMCDJetsSV);
      if (doprocessDataJetsSub || doprocessDataJetsSVSub || doprocessMCDJetsSub || doprocessMCDJetsSVSub || doprocessDataJetsTrigSub) {
        addDataDrivenHistograms(axisJetpTSub, "_sub", doprocessMCDJetsSub || doprocessMCDJetsSVSub);
      }
    }
  }

  // Filters
  Filter collisionFilter = nabs(aod::jcollision::posZ) < vertexZCut;
  Filter trackFilter = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax);
  Filter jetFilter = (aod::jet::pt >= jetPtMin && aod::jet::pt < jetPtMax && aod::jet::eta < jetEtaMax - aod::jet::r / 100.f && aod::jet::eta > jetEtaMin + aod::jet::r / 100.f);

  // Data
  using AnalysisCollisions = soa::Join<aod::JetCollisions, aod::JCollisionPIs, aod::BkgChargedRhos>;
  using FilteredCollisions = soa::Filtered<AnalysisCollisions>;
  using AnalysisCollisionsTriggered = soa::Join<aod::JetCollisions, aod::JCollisionPIs, aod::BkgChargedRhos>;
  using FilteredCollisionsTriggered = soa::Filtered<AnalysisCollisionsTriggered>;
  // using OrigCollisions = soa::Join<aod::Collisions, aod::PVMults>;

  using DataJets = soa::Join<aod::ChargedJets, aod::ChargedJetConstituents, aod::ChargedJetTags>;
  using FilteredDataJets = soa::Filtered<DataJets>;
  using DataJetsSV = soa::Join<DataJets, aod::DataSecondaryVertex3ProngIndices>;
  using FilteredDataJetsSV = soa::Filtered<DataJetsSV>;

  using AnalysisTracks = soa::Join<aod::JetTracks, aod::JTrackExtras, aod::JTrackPIs>;
  using FilteredTracks = soa::Filtered<AnalysisTracks>;
  // using OriginalTracks = soa::Join<aod::Tracks, aod::TracksCov, aod::TrackSelection, aod::TracksDCA, aod::TracksDCACov, aod::TracksExtra>;

  // MCD
  using AnalysisCollisionsMCD = soa::Join<aod::JetCollisionsMCD, aod::JCollisionPIs, aod::BkgChargedRhos, aod::JCollisionOutliers>;
  using FilteredCollisionsMCD = soa::Filtered<AnalysisCollisionsMCD>;

  using MCDJets = soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets, aod::ChargedMCDetectorLevelJetFlavourDef, aod::ChargedMCDetectorLevelJetTags>;
  using FilteredMCDJets = soa::Filtered<MCDJets>;
  using MCDJetsSV = soa::Join<MCDJets, aod::MCDSecondaryVertex3ProngIndices>;
  using FilteredMCDJetsSV = soa::Filtered<MCDJetsSV>;

  using AnalysisTracksMCD = soa::Join<aod::JetTracksMCD, aod::JTrackExtras, aod::JTrackPIs>;
  using FilteredTracksMCD = soa::Filtered<AnalysisTracksMCD>;

  // MCP
  using AnalysisCollisionsMCP = soa::Join<aod::JetMcCollisions, aod::JMcCollisionPIs, aod::BkgChargedMcRhos, aod::JMcCollisionOutliers>;

  using MCPJets = soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets, aod::ChargedMCParticleLevelJetFlavourDef>;
  using FilteredMCPJets = soa::Filtered<MCPJets>;

  template <typename AnalysisJet, typename AnyTracks>
  int analyzeJetTrackInfo(AnalysisJet const& analysisJet, AnyTracks const& /*allTracks*/ /*, int8_t jetFlavor = 0, double weight = 1.0*/)
  {
    int nTracks = 0;
    for (const auto& constituent : analysisJet.template tracks_as<AnyTracks>()) {

      if (constituent.pt() < trackPtMin || !jettaggingutilities::trackAcceptanceWithDca(constituent, maxIPxy, maxIPz)) {
        continue;
      }

      int sign = jettaggingutilities::getGeoSign(analysisJet, constituent);
      // auto origConstit = constituent.template track_as<AnyOriginalTracks>();

      registry.fill(HIST("h_gnnfeat_trackpT"), constituent.pt());
      registry.fill(HIST("h_gnnfeat_trackPhi"), constituent.phi());
      registry.fill(HIST("h_gnnfeat_trackEta"), constituent.eta());
      registry.fill(HIST("h_gnnfeat_trackCharge"), static_cast<float>(constituent.sign()));
      registry.fill(HIST("h_gnnfeat_trackDCAxy"), std::abs(constituent.dcaXY()) * sign);
      registry.fill(HIST("h_gnnfeat_trackSigmaDCAxy"), constituent.sigmadcaXY());
      registry.fill(HIST("h_gnnfeat_trackDCAz"), std::abs(constituent.dcaZ()) * sign);
      registry.fill(HIST("h_gnnfeat_trackSigmaDCAz"), constituent.sigmadcaZ());
      // registry.fill(HIST("h_gnnfeat_trackITSNCls"), static_cast<float>(origConstit.itsNCls()));
      // registry.fill(HIST("h_gnnfeat_trackTPCNCls"), static_cast<float>(origConstit.tpcNClsFound()));
      // registry.fill(HIST("h_gnnfeat_trackTPCNCrossedRows"), static_cast<float>(origConstit.tpcNClsCrossedRows()));
      // registry.fill(HIST("h_gnnfeat_trackITSChi2NCl"), origConstit.itsChi2NCl());
      // registry.fill(HIST("h_gnnfeat_trackTPCChi2NCl"), origConstit.tpcChi2NCl());

      registry.fill(HIST("h_gnnfeat_tracksIPxy"), std::abs(constituent.dcaXY()) * sign / constituent.sigmadcaXY());
      registry.fill(HIST("h_gnnfeat_tracksIPz"), std::abs(constituent.dcaZ()) * sign / constituent.sigmadcaZ());

      ++nTracks;
    }
    return nTracks;
  }

  static constexpr float largeNegativeNumber = -98.0f;
  static constexpr float largePositiveNumber = 9999.0f;

  template <typename AnyTracks, typename AnalysisJet>
  bool isAcceptedJet(AnalysisJet const& jet)
  {
    if (jetAreaFractionMin > largeNegativeNumber) {
      if (jet.area() < jetAreaFractionMin * o2::constants::math::PI * (jet.r() / 100.0) * (jet.r() / 100.0)) {
        return false;
      }
    }
    bool checkConstituentPt = true;
    bool checkConstituentMinPt = (leadingConstituentPtMin > largeNegativeNumber);
    bool checkConstituentMaxPt = (leadingConstituentPtMax < largePositiveNumber);
    if (!checkConstituentMinPt && !checkConstituentMaxPt) {
      checkConstituentPt = false;
    }

    if (checkConstituentPt) {
      bool isMinLeadingConstituent = !checkConstituentMinPt;
      bool isMaxLeadingConstituent = true;

      for (const auto& constituent : jet.template tracks_as<AnyTracks>()) {
        double pt = constituent.pt();

        if (checkConstituentMinPt && pt >= leadingConstituentPtMin) {
          isMinLeadingConstituent = true;
        }
        if (checkConstituentMaxPt && pt > leadingConstituentPtMax) {
          isMaxLeadingConstituent = false;
        }
      }
      return isMinLeadingConstituent && isMaxLeadingConstituent;
    }

    return true;
  }

  // Jet-flavour categories collapsed onto a single axis (instead of one histogram per category): b, c, lf.
  static double getJetFlavorCat(int8_t jetFlavor)
  {
    if (jetFlavor == JetTaggingSpecies::beauty) {
      return 0.;
    }
    if (jetFlavor == JetTaggingSpecies::charm) {
      return 1.;
    }
    return 2.; // lf
  }

  static void setJetFlavorCatAxisLabels(TAxis* axis)
  {
    axis->SetBinLabel(1, "b");
    axis->SetBinLabel(2, "c");
    axis->SetBinLabel(3, "lf");
  }

  // Fills one bin per satisfied event-selection stage in [first, last] (see kEvtSelStages). Bin index
  // == stage index i (not offset by `first`), since the histogram's axis always spans the full table.
  // Callers pass the "_sub"-suffixed HIST name directly for a UE-subtracted or MC-particle-context variant
  // of the same histogram (single `registry`, see its declaration).
  void fillEvtSelStages(auto const& histName, double xVal, EvtSelFlag evtselCode, double weight, int first, int last)
  {
    for (int i = first; i <= last; ++i) {
      if (hasAll(evtselCode, kEvtSelStages[i].flag)) {
        registry.fill(histName, xVal, static_cast<double>(i), weight);
      }
    }
  }

  // lf-jets further split by MC-particle-jet match status (see h2_jetpT_Db_lfmatched/h2_jetpT_Db_lfnone).
  static double getLfMatchStatus(int8_t jetFlavor)
  {
    return (jetFlavor == JetTaggingSpecies::none) ? 1. : 0.;
  }

  static void setLfMatchStatusAxisLabels(TAxis* axis)
  {
    axis->SetBinLabel(1, "matched");
    axis->SetBinLabel(2, "none");
  }

  // Selects the SV with the highest decay length significance out of all SVs matched to the jet,
  // shared by fillDataJetHistogramsSV()/fillMCDJetHistogramsSV(). Fills the inclusive SV histograms
  // (plus the flavour-split variants when jetFlavor >= 0, i.e. for MCD jets) and returns the
  // selected SV's mass (-1 if the jet has no SV). `jetpT` is passed in already raw-or-UE-subtracted by
  // the caller (per `withSub`) rather than recomputed here, since this function doesn't otherwise need rho.
  template <bool withSub, typename AnalysisJet, typename SecondaryVertices>
  float fillSVHistograms(AnalysisJet const& analysisJet, SecondaryVertices const& /*allSVs*/, float jetpT, double weightEvt = 1.0, int8_t jetFlavor = -1)
  {
    auto svs = analysisJet.template secondaryVertices_as<SecondaryVertices>();
    if (svs.size() == 0) {
      // Not jet-pT-dependent (no-SV sentinel) - one shared histogram regardless of `withSub`.
      registry.fill(HIST("h_SVMass"), -1.f, weightEvt);
      return -1.f;
    }
    const auto& sv = *std::max_element(svs.begin(), svs.end(), [](auto const& svA, auto const& svB) {
      return (svA.decayLengthXY() / svA.errorDecayLengthXY()) < (svB.decayLengthXY() / svB.errorDecayLengthXY());
    });
    float massSV = sv.m();
    float lxyS = sv.decayLengthXY() / sv.errorDecayLengthXY();
    // h_SVMass/h2_SVMass_flavor don't use jetpT as a fill value - shared regardless of `withSub`. The
    // jetpT-axis histograms below get the "_sub" name instead of a second registry (see `registry`'s
    // declaration).
    registry.fill(HIST("h_SVMass"), massSV, weightEvt);
    if (jetFlavor >= 0) {
      registry.fill(HIST("h2_SVMass_flavor"), massSV, getJetFlavorCat(jetFlavor), weightEvt);
    }
    if constexpr (withSub) {
      registry.fill(HIST("h2_SVMass_jetpT_sub"), jetpT, massSV, weightEvt);
      registry.fill(HIST("h2_SVLxyS_jetpT_sub"), jetpT, lxyS, weightEvt);
      registry.fill(HIST("h2_SVDispersion_jetpT_sub"), jetpT, sv.dispersion(), weightEvt);
      if (jetFlavor == JetTaggingSpecies::beauty) {
        registry.fill(HIST("h2_SVMass_jetpT_b_sub"), jetpT, massSV, weightEvt);
        registry.fill(HIST("h2_SVLxyS_jetpT_b_sub"), jetpT, lxyS, weightEvt);
        registry.fill(HIST("h2_SVDispersion_jetpT_b_sub"), jetpT, sv.dispersion(), weightEvt);
      } else if (jetFlavor == JetTaggingSpecies::charm) {
        registry.fill(HIST("h2_SVMass_jetpT_c_sub"), jetpT, massSV, weightEvt);
        registry.fill(HIST("h2_SVLxyS_jetpT_c_sub"), jetpT, lxyS, weightEvt);
        registry.fill(HIST("h2_SVDispersion_jetpT_c_sub"), jetpT, sv.dispersion(), weightEvt);
      } else if (jetFlavor >= 0) {
        registry.fill(HIST("h2_SVMass_jetpT_lf_sub"), jetpT, massSV, weightEvt);
        registry.fill(HIST("h2_SVLxyS_jetpT_lf_sub"), jetpT, lxyS, weightEvt);
        registry.fill(HIST("h2_SVDispersion_jetpT_lf_sub"), jetpT, sv.dispersion(), weightEvt);
      }
    } else {
      registry.fill(HIST("h2_SVMass_jetpT"), jetpT, massSV, weightEvt);
      registry.fill(HIST("h2_SVLxyS_jetpT"), jetpT, lxyS, weightEvt);
      registry.fill(HIST("h2_SVDispersion_jetpT"), jetpT, sv.dispersion(), weightEvt);
      if (jetFlavor == JetTaggingSpecies::beauty) {
        registry.fill(HIST("h2_SVMass_jetpT_b"), jetpT, massSV, weightEvt);
        registry.fill(HIST("h2_SVLxyS_jetpT_b"), jetpT, lxyS, weightEvt);
        registry.fill(HIST("h2_SVDispersion_jetpT_b"), jetpT, sv.dispersion(), weightEvt);
      } else if (jetFlavor == JetTaggingSpecies::charm) {
        registry.fill(HIST("h2_SVMass_jetpT_c"), jetpT, massSV, weightEvt);
        registry.fill(HIST("h2_SVLxyS_jetpT_c"), jetpT, lxyS, weightEvt);
        registry.fill(HIST("h2_SVDispersion_jetpT_c"), jetpT, sv.dispersion(), weightEvt);
      } else if (jetFlavor >= 0) {
        registry.fill(HIST("h2_SVMass_jetpT_lf"), jetpT, massSV, weightEvt);
        registry.fill(HIST("h2_SVLxyS_jetpT_lf"), jetpT, lxyS, weightEvt);
        registry.fill(HIST("h2_SVDispersion_jetpT_lf"), jetpT, sv.dispersion(), weightEvt);
      }
    }
    return massSV;
  }

  // `withSub` selects UE(rho*area)-subtracted jet pT + the "_sub"-suffixed histogram set in the single
  // `registry` instead of raw pT + the unsuffixed names - used by processDataJets(SV)/processDataJets(SV)Sub.
  template <bool withSub, typename AnalysisJet, typename AnyTracks>
  int fillDataJetHistograms(AnalysisJet const& analysisJet, AnyTracks const& allTracks, float rho = 0.0)
  {
    int nTracks = analyzeJetTrackInfo(analysisJet, allTracks);

    float jetpT = withSub ? (analysisJet.pt() - rho * analysisJet.area()) : analysisJet.pt();
    // h_jetPhi/h_jetEta/h_jetMass/h_Db/h2_nTracks_Db don't use jetpT as a fill value - one shared histogram
    // regardless of `withSub` (see addCoreJetHistograms()'s init() comment). Only h_jetpT/h2_jetpT_Db/
    // hSparse_Incljets actually differ between the raw and Sub call, so only those get the "_sub" name.
    registry.fill(HIST("h_jetPhi"), analysisJet.phi());
    registry.fill(HIST("h_jetEta"), analysisJet.eta());
    registry.fill(HIST("h_jetMass"), analysisJet.mass());
    registry.fill(HIST("h_Db"), analysisJet.scoreML());
    registry.fill(HIST("h2_nTracks_Db"), nTracks, analysisJet.scoreML());

    if constexpr (withSub) {
      registry.fill(HIST("h_jetpT_sub"), jetpT);
      registry.fill(HIST("h2_jetpT_Db_sub"), jetpT, analysisJet.scoreML());

      if (doDataDriven && !doDataDrivenSV) {
        if (doDataDrivenExtra) {
          registry.fill(HIST("hSparse_Incljets_sub"), jetpT, analysisJet.scoreML(), nTracks, analysisJet.mass());
        } else {
          registry.fill(HIST("hSparse_Incljets_sub"), jetpT, analysisJet.scoreML(), nTracks);
        }
      }
    } else {
      registry.fill(HIST("h_jetpT"), jetpT);
      registry.fill(HIST("h2_jetpT_Db"), jetpT, analysisJet.scoreML());

      if (doDataDriven && !doDataDrivenSV) {
        if (doDataDrivenExtra) {
          registry.fill(HIST("hSparse_Incljets"), jetpT, analysisJet.scoreML(), nTracks, analysisJet.mass());
        } else {
          registry.fill(HIST("hSparse_Incljets"), jetpT, analysisJet.scoreML(), nTracks);
        }
      }
    }
    return nTracks;
  }

  template <bool withSub, typename AnalysisJet, typename AnyTracks, typename AnySVs>
  void fillDataJetHistogramsSV(AnalysisJet const& analysisJet, AnyTracks const& allTracks, AnySVs const& allSVs, float rho = 0.0)
  {
    int nTracks = fillDataJetHistograms<withSub>(analysisJet, allTracks, rho);
    float jetpT = withSub ? (analysisJet.pt() - rho * analysisJet.area()) : analysisJet.pt();
    float massSV = fillSVHistograms<withSub>(analysisJet, allSVs, jetpT);

    if (doDataDriven && doDataDrivenSV) {
      if constexpr (withSub) {
        registry.fill(HIST("hSparse_Incljets_sub"), jetpT, analysisJet.scoreML(), nTracks, massSV);
      } else {
        registry.fill(HIST("hSparse_Incljets"), jetpT, analysisJet.scoreML(), nTracks, massSV);
      }
    }
  }

  // jetFlavor: JetTaggingSpecies of the jet; nTracks: number of GNN-input constituents (reused by fillMCDJetHistogramsSV()).
  struct JetHistFillResult {
    int8_t jetFlavor;
    // cppcheck-suppress unusedStructMember
    int nTracks;
  };

  // Constituent loop factored out of fillMCDJetHistograms(): the h_gnnfeat_* track-QA histograms it fills
  // don't depend on the jet's own (possibly UE-subtracted) pT, so they always use the unsuffixed names
  // regardless of `withSub` - no "_sub" duplicates needed.
  template <typename AnalysisJet, typename AnyTracks>
  std::pair<int, int> analyzeMCDJetTrackInfo(AnalysisJet const& analysisJet, AnyTracks const& /*allTracks*/, double weightEvt = 1.0)
  {
    int nTracks = 0;
    int nNppTracks = 0;
    for (const auto& constituent : analysisJet.template tracks_as<AnyTracks>()) {
      if (constituent.pt() < trackPtMinGnn) {
        continue;
      }
      if (!constituent.has_mcParticle() || !constituent.template mcParticle_as<aod::JetParticles>().isPhysicalPrimary()) {
        ++nNppTracks;
      }

      int sign = jettaggingutilities::getGeoSign(analysisJet, constituent);
      // auto origConstit = constituent.template track_as<AnyOriginalTracks>();

      registry.fill(HIST("h_gnnfeat_trackpT"), constituent.pt(), weightEvt);
      registry.fill(HIST("h_gnnfeat_trackPhi"), constituent.phi(), weightEvt);
      registry.fill(HIST("h_gnnfeat_trackEta"), constituent.eta(), weightEvt);
      registry.fill(HIST("h_gnnfeat_trackCharge"), static_cast<float>(constituent.sign()), weightEvt);
      registry.fill(HIST("h_gnnfeat_trackDCAxy"), std::abs(constituent.dcaXY()) * sign, weightEvt);
      registry.fill(HIST("h_gnnfeat_trackSigmaDCAxy"), constituent.sigmadcaXY(), weightEvt);
      registry.fill(HIST("h_gnnfeat_trackDCAz"), std::abs(constituent.dcaZ()) * sign, weightEvt);
      registry.fill(HIST("h_gnnfeat_trackSigmaDCAz"), constituent.sigmadcaZ(), weightEvt);
      // registry.fill(HIST("h_gnnfeat_trackITSNCls"), static_cast<float>(origConstit.itsNCls()), weightEvt);
      // registry.fill(HIST("h_gnnfeat_trackTPCNCls"), static_cast<float>(origConstit.tpcNClsFound()), weightEvt);
      // registry.fill(HIST("h_gnnfeat_trackTPCNCrossedRows"), static_cast<float>(origConstit.tpcNClsCrossedRows()), weightEvt);
      // registry.fill(HIST("h_gnnfeat_trackITSChi2NCl"), origConstit.itsChi2NCl(), weightEvt);
      // registry.fill(HIST("h_gnnfeat_trackTPCChi2NCl"), origConstit.tpcChi2NCl(), weightEvt);

      registry.fill(HIST("h_gnnfeat_tracksIPxy"), std::abs(constituent.dcaXY()) * sign / constituent.sigmadcaXY(), weightEvt);
      registry.fill(HIST("h_gnnfeat_tracksIPz"), std::abs(constituent.dcaZ()) * sign / constituent.sigmadcaZ(), weightEvt);

      ++nTracks;
    }
    return {nTracks, nNppTracks};
  }

  // `withSub` selects UE(rho*area)-subtracted jet pT + the "_sub"-suffixed histogram set, same convention as
  // fillDataJetHistograms(). Histograms whose value doesn't depend on jetpT are shared, unsuffixed.
  template <bool withSub, typename AnalysisJet, typename AnyTracks>
  JetHistFillResult fillMCDJetHistograms(AnalysisJet const& analysisJet, AnyTracks const& allTracks, float rho = 0.f, double weightEvt = 1.0)
  {
    int8_t jetFlavor = analysisJet.origin();
    auto [nTracks, nNppTracks] = analyzeMCDJetTrackInfo(analysisJet, allTracks, weightEvt);

    float jetpT = withSub ? (analysisJet.pt() - rho * analysisJet.area()) : analysisJet.pt();
    double flavorCat = getJetFlavorCat(jetFlavor);

    registry.fill(HIST("h_jetPhi"), analysisJet.phi(), weightEvt);
    registry.fill(HIST("h_jetEta"), analysisJet.eta(), weightEvt);
    registry.fill(HIST("h_jetMass"), analysisJet.mass(), weightEvt);
    registry.fill(HIST("h_Db"), analysisJet.scoreML(), weightEvt);
    registry.fill(HIST("h2_nTracks_Db"), nTracks, analysisJet.scoreML(), weightEvt);
    registry.fill(HIST("h2_Db_flavor"), analysisJet.scoreML(), flavorCat, weightEvt);
    if (jetFlavor == JetTaggingSpecies::beauty) {
      registry.fill(HIST("h2_nTracks_Db_b"), nTracks, analysisJet.scoreML(), weightEvt);
    } else if (jetFlavor == JetTaggingSpecies::charm) {
      registry.fill(HIST("h2_nTracks_Db_c"), nTracks, analysisJet.scoreML(), weightEvt);
    } else {
      registry.fill(HIST("h2_nTracks_Db_lf"), nTracks, analysisJet.scoreML(), weightEvt);
    }
    bool isLf = jetFlavor != JetTaggingSpecies::beauty && jetFlavor != JetTaggingSpecies::charm;
    if (isLf) {
      registry.fill(HIST("h2_Db_lfmatch"), analysisJet.scoreML(), getLfMatchStatus(jetFlavor), weightEvt);
    }
    bool isNpp = static_cast<float>(nNppTracks) / nTracks > trackNppCrit;
    if (isNpp) {
      registry.fill(HIST("h_Db_npp"), analysisJet.scoreML(), weightEvt);
      registry.fill(HIST("h2_Db_npp_flavor"), analysisJet.scoreML(), flavorCat, weightEvt);
      if (isLf) {
        registry.fill(HIST("h2_Db_npp_lfmatch"), analysisJet.scoreML(), getLfMatchStatus(jetFlavor), weightEvt);
      }
    }

    if constexpr (withSub) {
      registry.fill(HIST("h_jetpT_sub"), jetpT, weightEvt);
      registry.fill(HIST("h2_jetpT_Db_sub"), jetpT, analysisJet.scoreML(), weightEvt);
      registry.fill(HIST("h2_jetpT_flavor_sub"), jetpT, flavorCat, weightEvt);
      if (jetFlavor == JetTaggingSpecies::beauty) {
        registry.fill(HIST("h2_jetpT_Db_b_sub"), jetpT, analysisJet.scoreML(), weightEvt);
      } else if (jetFlavor == JetTaggingSpecies::charm) {
        registry.fill(HIST("h2_jetpT_Db_c_sub"), jetpT, analysisJet.scoreML(), weightEvt);
      } else {
        registry.fill(HIST("h2_jetpT_Db_lf_sub"), jetpT, analysisJet.scoreML(), weightEvt);
      }
      if (isLf) {
        registry.fill(HIST("h2_jetpT_lfmatch_sub"), jetpT, getLfMatchStatus(jetFlavor), weightEvt);
        if (jetFlavor == JetTaggingSpecies::none) {
          registry.fill(HIST("h2_jetpT_Db_lfnone_sub"), jetpT, analysisJet.scoreML(), weightEvt);
        } else {
          registry.fill(HIST("h2_jetpT_Db_lfmatched_sub"), jetpT, analysisJet.scoreML(), weightEvt);
        }
      }
      if (isNpp) {
        registry.fill(HIST("h2_jetpT_Db_npp_sub"), jetpT, analysisJet.scoreML(), weightEvt);
        if (jetFlavor == JetTaggingSpecies::beauty) {
          registry.fill(HIST("h2_jetpT_Db_npp_b_sub"), jetpT, analysisJet.scoreML(), weightEvt);
        } else if (jetFlavor == JetTaggingSpecies::charm) {
          registry.fill(HIST("h2_jetpT_Db_npp_c_sub"), jetpT, analysisJet.scoreML(), weightEvt);
        } else {
          registry.fill(HIST("h2_jetpT_Db_npp_lf_sub"), jetpT, analysisJet.scoreML(), weightEvt);
        }
      }

      if (doDataDriven && !doDataDrivenSV) {
        if (doDataDrivenExtra) {
          registry.fill(HIST("hSparse_Incljets_sub"), jetpT, analysisJet.scoreML(), nTracks, analysisJet.mass(), weightEvt);
          if (jetFlavor == JetTaggingSpecies::beauty) {
            registry.fill(HIST("hSparse_bjets_sub"), jetpT, analysisJet.scoreML(), nTracks, analysisJet.mass(), weightEvt);
          } else if (jetFlavor == JetTaggingSpecies::charm) {
            registry.fill(HIST("hSparse_cjets_sub"), jetpT, analysisJet.scoreML(), nTracks, analysisJet.mass(), weightEvt);
          } else {
            registry.fill(HIST("hSparse_lfjets_sub"), jetpT, analysisJet.scoreML(), nTracks, analysisJet.mass(), weightEvt);
            if (jetFlavor == JetTaggingSpecies::none) {
              registry.fill(HIST("hSparse_lfjets_none_sub"), jetpT, analysisJet.scoreML(), nTracks, analysisJet.mass(), weightEvt);
            } else {
              registry.fill(HIST("hSparse_lfjets_matched_sub"), jetpT, analysisJet.scoreML(), nTracks, analysisJet.mass(), weightEvt);
            }
          }
        } else {
          registry.fill(HIST("hSparse_Incljets_sub"), jetpT, analysisJet.scoreML(), nTracks, weightEvt);
          if (jetFlavor == JetTaggingSpecies::beauty) {
            registry.fill(HIST("hSparse_bjets_sub"), jetpT, analysisJet.scoreML(), nTracks, weightEvt);
          } else if (jetFlavor == JetTaggingSpecies::charm) {
            registry.fill(HIST("hSparse_cjets_sub"), jetpT, analysisJet.scoreML(), nTracks, weightEvt);
          } else {
            registry.fill(HIST("hSparse_lfjets_sub"), jetpT, analysisJet.scoreML(), nTracks, weightEvt);
            if (jetFlavor == JetTaggingSpecies::none) {
              registry.fill(HIST("hSparse_lfjets_none_sub"), jetpT, analysisJet.scoreML(), nTracks, weightEvt);
            } else {
              registry.fill(HIST("hSparse_lfjets_matched_sub"), jetpT, analysisJet.scoreML(), nTracks, weightEvt);
            }
          }
        }
      }
    } else {
      registry.fill(HIST("h_jetpT"), jetpT, weightEvt);
      registry.fill(HIST("h2_jetpT_Db"), jetpT, analysisJet.scoreML(), weightEvt);
      registry.fill(HIST("h2_jetpT_flavor"), jetpT, flavorCat, weightEvt);
      if (jetFlavor == JetTaggingSpecies::beauty) {
        registry.fill(HIST("h2_jetpT_Db_b"), jetpT, analysisJet.scoreML(), weightEvt);
      } else if (jetFlavor == JetTaggingSpecies::charm) {
        registry.fill(HIST("h2_jetpT_Db_c"), jetpT, analysisJet.scoreML(), weightEvt);
      } else {
        registry.fill(HIST("h2_jetpT_Db_lf"), jetpT, analysisJet.scoreML(), weightEvt);
      }
      if (isLf) {
        registry.fill(HIST("h2_jetpT_lfmatch"), jetpT, getLfMatchStatus(jetFlavor), weightEvt);
        if (jetFlavor == JetTaggingSpecies::none) {
          registry.fill(HIST("h2_jetpT_Db_lfnone"), jetpT, analysisJet.scoreML(), weightEvt);
        } else {
          registry.fill(HIST("h2_jetpT_Db_lfmatched"), jetpT, analysisJet.scoreML(), weightEvt);
        }
      }
      if (isNpp) {
        registry.fill(HIST("h2_jetpT_Db_npp"), jetpT, analysisJet.scoreML(), weightEvt);
        if (jetFlavor == JetTaggingSpecies::beauty) {
          registry.fill(HIST("h2_jetpT_Db_npp_b"), jetpT, analysisJet.scoreML(), weightEvt);
        } else if (jetFlavor == JetTaggingSpecies::charm) {
          registry.fill(HIST("h2_jetpT_Db_npp_c"), jetpT, analysisJet.scoreML(), weightEvt);
        } else {
          registry.fill(HIST("h2_jetpT_Db_npp_lf"), jetpT, analysisJet.scoreML(), weightEvt);
        }
      }

      if (doDataDriven && !doDataDrivenSV) {
        if (doDataDrivenExtra) {
          registry.fill(HIST("hSparse_Incljets"), jetpT, analysisJet.scoreML(), nTracks, analysisJet.mass(), weightEvt);
          if (jetFlavor == JetTaggingSpecies::beauty) {
            registry.fill(HIST("hSparse_bjets"), jetpT, analysisJet.scoreML(), nTracks, analysisJet.mass(), weightEvt);
          } else if (jetFlavor == JetTaggingSpecies::charm) {
            registry.fill(HIST("hSparse_cjets"), jetpT, analysisJet.scoreML(), nTracks, analysisJet.mass(), weightEvt);
          } else {
            registry.fill(HIST("hSparse_lfjets"), jetpT, analysisJet.scoreML(), nTracks, analysisJet.mass(), weightEvt);
            if (jetFlavor == JetTaggingSpecies::none) {
              registry.fill(HIST("hSparse_lfjets_none"), jetpT, analysisJet.scoreML(), nTracks, analysisJet.mass(), weightEvt);
            } else {
              registry.fill(HIST("hSparse_lfjets_matched"), jetpT, analysisJet.scoreML(), nTracks, analysisJet.mass(), weightEvt);
            }
          }
        } else {
          registry.fill(HIST("hSparse_Incljets"), jetpT, analysisJet.scoreML(), nTracks, weightEvt);
          if (jetFlavor == JetTaggingSpecies::beauty) {
            registry.fill(HIST("hSparse_bjets"), jetpT, analysisJet.scoreML(), nTracks, weightEvt);
          } else if (jetFlavor == JetTaggingSpecies::charm) {
            registry.fill(HIST("hSparse_cjets"), jetpT, analysisJet.scoreML(), nTracks, weightEvt);
          } else {
            registry.fill(HIST("hSparse_lfjets"), jetpT, analysisJet.scoreML(), nTracks, weightEvt);
            if (jetFlavor == JetTaggingSpecies::none) {
              registry.fill(HIST("hSparse_lfjets_none"), jetpT, analysisJet.scoreML(), nTracks, weightEvt);
            } else {
              registry.fill(HIST("hSparse_lfjets_matched"), jetpT, analysisJet.scoreML(), nTracks, weightEvt);
            }
          }
        }
      }
    }

    return {jetFlavor, nTracks};
  }

  // Same as fillMCDJetHistograms(), plus SV histograms (see fillSVHistograms()) and the SV-mass variant
  // of hSparse_Incljets (only when doDataDrivenSV is on); used by processMCDJetsSV.
  template <bool withSub, typename AnalysisJet, typename AnyTracks, typename AnySVs>
  int8_t fillMCDJetHistogramsSV(AnalysisJet const& analysisJet, AnyTracks const& allTracks, AnySVs const& allSVs, float rho = 0.f, double weightEvt = 1.0)
  {
    auto [jetFlavor, nTracks] = fillMCDJetHistograms<withSub>(analysisJet, allTracks, rho, weightEvt);
    float jetpT = withSub ? (analysisJet.pt() - rho * analysisJet.area()) : analysisJet.pt();
    float massSV = fillSVHistograms<withSub>(analysisJet, allSVs, jetpT, weightEvt, jetFlavor);

    if (doDataDriven && doDataDrivenSV) {
      if constexpr (withSub) {
        registry.fill(HIST("hSparse_Incljets_sub"), jetpT, analysisJet.scoreML(), nTracks, massSV, weightEvt);
      } else {
        registry.fill(HIST("hSparse_Incljets"), jetpT, analysisJet.scoreML(), nTracks, massSV, weightEvt);
      }
    }

    return jetFlavor;
  }

  // Check if the collision is INEL>0
  static constexpr int nPartInel0 = 3;
  template <typename MCColl, typename MCPart>
  bool isTrueINEL0(MCColl const& /*mccoll*/, MCPart const& mcparts)
  {
    for (const auto& mcparticle : mcparts) {
      if (!mcparticle.isPhysicalPrimary()) {
        continue;
      }
      const auto p = pdg->GetParticle(mcparticle.pdgCode());
      if (p != nullptr) {
        if (std::abs(p->Charge()) >= nPartInel0) {
          if (std::abs(mcparticle.eta()) < 1) {
            return true;
          }
        }
      }
    }
    return false;
  }

  // Shared entry checks for MCD collisions, used by processMCDJets(SV) and processMCDTracks:
  // event selection, MC-outlier rejection (own and matched mcCollision's), and eval/train split.
  template <typename AnalysisCollision>
  bool isRejectedMCDCollision(AnalysisCollision const& collision)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return true;
    }
    // Reject outlier MC collisions
    if (collision.isOutlier()) {
      return true;
    }
    if (collision.has_mcCollision() && collision.template mcCollision_as<AnalysisCollisionsMCP>().isOutlier()) {
      return true;
    }
    // Uses only collisionId % trainingDatasetRaioParam != 0 for evaluation dataset
    if (trainingDatasetRatioParam && collision.collisionId() % trainingDatasetRatioParam == 0) {
      return true;
    }
    return false;
  }

  template <typename AnalysisCollision>
  EvtSelFlag fillCollCounter(AnalysisCollision const& collision, float weightEvt = 1.f)
  {
    EvtSelFlag evtselCode = EvtSelFlag::Coll;
    registry.fill(HIST("hCollCounter"), static_cast<int>(EvtSel::Coll), weightEvt); // Coll

    bool zvtx = std::fabs(collision.posZ()) < vertexZCut;
    if (zvtx) {
      evtselCode |= EvtSelFlag::kZvtx;
      registry.fill(HIST("hCollCounter"), static_cast<int>(EvtSel::CollZvtx), weightEvt); // Coll+Zvtx
    }

    if (jetderiveddatautilities::selectCollision(collision, eventSelectionBitsTVX)) {
      evtselCode |= EvtSelFlag::kTVX;
      registry.fill(HIST("hCollCounter"), static_cast<int>(EvtSel::TVX), weightEvt); // Coll+TVX
      if (zvtx) {
        registry.fill(HIST("hCollCounter"), static_cast<int>(EvtSel::TVXZvtx), weightEvt); // Coll+TVX+Zvtx
      }
      if (jetderiveddatautilities::selectCollision(collision, eventSelectionBitsSelMC)) {
        evtselCode |= EvtSelFlag::kNoTFB;
        registry.fill(HIST("hCollCounter"), static_cast<int>(EvtSel::SelMC), weightEvt); // Coll+TVX+NoTFB
        if (zvtx) {
          registry.fill(HIST("hCollCounter"), static_cast<int>(EvtSel::SelMCZvtx), weightEvt); // Coll+TVX+NoTFB+Zvtx
        }
        if (jetderiveddatautilities::selectCollision(collision, eventSelectionBitsSel8)) {
          evtselCode |= EvtSelFlag::kNoITSROFB;
          registry.fill(HIST("hCollCounter"), static_cast<int>(EvtSel::Sel8), weightEvt); // Coll+TVX+NoTFB+NoITSROFB
          if (zvtx) {
            registry.fill(HIST("hCollCounter"), static_cast<int>(EvtSel::Sel8Zvtx), weightEvt); // Coll+TVX+NoTFB+NoITSROFB+Zvtx
            // if (collision.template collision_as<OrigCollisions>().isInelGt0()) {
            //   evtselCode |= EvtSelFlag::kINELgt0rec;
            //   registry.fill(HIST("hCollCounter"), static_cast<int>(EvtSel::INELgt0rec), weightEvt); // INELgt0+Zvtx(rec)
            // }
          }
          if (jetderiveddatautilities::selectCollision(collision, eventSelectionBitsSel8Full)) {
            evtselCode |= EvtSelFlag::kNoPileup;
            registry.fill(HIST("hCollCounter"), static_cast<int>(EvtSel::Sel8Full), weightEvt); // Coll+TVX+NoTFB+NoITSROFB+NoPileup
            if (zvtx) {
              registry.fill(HIST("hCollCounter"), static_cast<int>(EvtSel::Sel8FullZvtx), weightEvt); // Coll+TVX+NoTFB+NoITSROFB+NoPileup+Zvtx
            }
            if (jetderiveddatautilities::selectCollision(collision, eventSelectionBitsSel8FullGood)) {
              evtselCode |= EvtSelFlag::kIsGoodZvtx;
              registry.fill(HIST("hCollCounter"), static_cast<int>(EvtSel::Sel8FullGood), weightEvt); // Coll+TVX+NoTFB+NoITSROFB+NoPileup+GoodZvtx
              if (zvtx) {
                registry.fill(HIST("hCollCounter"), static_cast<int>(EvtSel::Sel8FullGoodZvtx), weightEvt); // Coll+TVX+NoTFB+NoITSROFB+NoPileup+GoodZvtx+Zvtx
              }
            }
          }
        }
      }
    }

    return evtselCode;
  }

  // Same stage ladder as fillCollCounter(), but evaluated directly on the MC-truth BC's raw evsel bits
  // (mcCollision.bc_as<aod::JBCs>(), a mandatory 1:1 index -> always resolves, no "has BC" gate needed)
  // instead of a matched reco collision's packed EventSel word. Reuses EvtSelFlag::kColl/EvtSel::Coll(+Zvtx)
  // to mean "BC" here (a BC always exists for an McCollision, so that stage is trivially satisfied) purely
  // so this cascade lines up bin-for-bin with hCollCounter/hMcCollCounter.
  template <typename AnalysisBC>
  EvtSelFlag fillMcCollCounterBC(AnalysisBC const& bc, float posZ, float weightEvt = 1.f)
  {
    EvtSelFlag evtselCode = EvtSelFlag::Coll;                                           // "Coll" == BC here, see above
    registry.fill(HIST("hMcCollCounterBC"), static_cast<int>(EvtSel::Coll), weightEvt); // BC

    bool zvtx = std::fabs(posZ) < vertexZCut;
    if (zvtx) {
      evtselCode |= EvtSelFlag::kZvtx;
      registry.fill(HIST("hMcCollCounterBC"), static_cast<int>(EvtSel::CollZvtx), weightEvt); // BC+Zvtx
    }

    if (bc.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
      evtselCode |= EvtSelFlag::kTVX;
      registry.fill(HIST("h_vertexZ_truth_bctvx"), posZ, weightEvt);
      registry.fill(HIST("hMcCollCounterBC"), static_cast<int>(EvtSel::TVX), weightEvt); // BC+TVX
      if (zvtx) {
        registry.fill(HIST("hMcCollCounterBC"), static_cast<int>(EvtSel::TVXZvtx), weightEvt); // BC+TVX, +Zvtx
      }
      if (bc.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
        evtselCode |= EvtSelFlag::kNoTFB;
        registry.fill(HIST("hMcCollCounterBC"), static_cast<int>(EvtSel::SelMC), weightEvt); // BC+TVX+NoTFB
        if (zvtx) {
          registry.fill(HIST("hMcCollCounterBC"), static_cast<int>(EvtSel::SelMCZvtx), weightEvt); // BC+TVX+NoTFB, +Zvtx
        }
        if (bc.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
          evtselCode |= EvtSelFlag::kNoITSROFB;
          registry.fill(HIST("hMcCollCounterBC"), static_cast<int>(EvtSel::Sel8), weightEvt); // BC+TVX+NoTFB+NoITSROFB
          if (zvtx) {
            registry.fill(HIST("hMcCollCounterBC"), static_cast<int>(EvtSel::Sel8Zvtx), weightEvt); // BC+TVX+NoTFB+NoITSROFB, +Zvtx
          }
          if (bc.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
            evtselCode |= EvtSelFlag::kNoPileup;
            registry.fill(HIST("hMcCollCounterBC"), static_cast<int>(EvtSel::Sel8Full), weightEvt); // BC+TVX+NoTFB+NoITSROFB+NoPileup
            if (zvtx) {
              registry.fill(HIST("hMcCollCounterBC"), static_cast<int>(EvtSel::Sel8FullZvtx), weightEvt); // BC+TVX+NoTFB+NoITSROFB+NoPileup, +Zvtx
            }
            if (bc.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
              evtselCode |= EvtSelFlag::kIsGoodZvtx;
              registry.fill(HIST("hMcCollCounterBC"), static_cast<int>(EvtSel::Sel8FullGood), weightEvt); // BC+TVX+NoTFB+NoITSROFB+NoPileup+GoodZvtx
              if (zvtx) {
                registry.fill(HIST("hMcCollCounterBC"), static_cast<int>(EvtSel::Sel8FullGoodZvtx), weightEvt); // BC+TVX+NoTFB+NoITSROFB+NoPileup+GoodZvtx, +Zvtx
              }
            }
          }
        }
      }
    }

    return evtselCode;
  }

  // Initialize CCDB access and histogram registry for Zorro processing
  template <typename BCType>
  void initCCDB(const BCType& bc)
  {
    if (doSoftwareTriggerSelection) {
      zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), triggerMasks.value);
      zorro.populateHistRegistry(registry, bc.runNumber());
    }
  }

  void processDummy(FilteredCollisions::iterator const& /*collision*/)
  {
  }
  PROCESS_SWITCH(BjetTaggingGnn, processDummy, "Dummy process function turned on by default", true);

  // Shared body of processDataJets(Sub)()/processDataJetsSV(Sub)(); see processMCDJetsCommon() for the
  // withSV/withSub convention.
  template <bool withSV, bool withSub, typename AnalysisJets, typename SecondaryVertices = void>
  void processDataJetsCommon(FilteredCollisions::iterator const& collision, AnalysisJets const& alljets, FilteredTracks const& allTracks, SecondaryVertices const* allSVs = nullptr)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    // Keeps hCollCounter populated consistently across every process[Data,MCD]Jets(SV)(Sub) variant, not
    // just the Sel ones - these should all behave identically except for which jet histograms get filled.
    fillCollCounter(collision);
    registry.fill(HIST("h_event_counter"), 0.0); // Coll+TVX+NoTFB+NoITSROFB+Zvtx

    registry.fill(HIST("h_vertexZ"), collision.posZ());

    float rho = collision.rho();

    for (const auto& analysisJet : alljets) {
      if (!isAcceptedJet<FilteredTracks>(analysisJet)) {
        continue;
      }

      if constexpr (withSV) {
        fillDataJetHistogramsSV<withSub>(analysisJet, allTracks, *allSVs, rho);
      } else {
        fillDataJetHistograms<withSub>(analysisJet, allTracks, rho);
      }
    }
  }

  void processDataJets(FilteredCollisions::iterator const& collision, FilteredDataJets const& alljets, FilteredTracks const& allTracks)
  {
    processDataJetsCommon<false, false>(collision, alljets, allTracks);
  }
  PROCESS_SWITCH(BjetTaggingGnn, processDataJets, "jet information in Data", false);

  void processDataJetsSub(FilteredCollisions::iterator const& collision, FilteredDataJets const& alljets, FilteredTracks const& allTracks)
  {
    processDataJetsCommon<false, true>(collision, alljets, allTracks);
  }
  PROCESS_SWITCH(BjetTaggingGnn, processDataJetsSub, "jet information in Data, UE(rho*area)-subtracted jet pT", false);

  void processDataJetsSV(FilteredCollisions::iterator const& collision, FilteredDataJetsSV const& alljets, aod::DataSecondaryVertex3Prongs const& allSVs, FilteredTracks const& allTracks)
  {
    processDataJetsCommon<true, false>(collision, alljets, allTracks, &allSVs);
  }
  PROCESS_SWITCH(BjetTaggingGnn, processDataJetsSV, "jet information in Data with secondary vertex info", false);

  void processDataJetsSVSub(FilteredCollisions::iterator const& collision, FilteredDataJetsSV const& alljets, aod::DataSecondaryVertex3Prongs const& allSVs, FilteredTracks const& allTracks)
  {
    processDataJetsCommon<true, true>(collision, alljets, allTracks, &allSVs);
  }
  PROCESS_SWITCH(BjetTaggingGnn, processDataJetsSVSub, "jet information in Data with secondary vertex info, UE(rho*area)-subtracted jet pT", false);

  void processDataJetsTrig(FilteredCollisionsTriggered::iterator const& collision, FilteredDataJets const& alljets, FilteredTracks const& allTracks, aod::JBCs const& /*bcInfo*/)
  {
    // Get BC info associated with the collision before applying any event selections
    const auto& bc = collision.bc_as<aod::JBCs>();
    // Initialize CCDB objects using the BC info
    initCCDB(bc);
    // If SoftwareTriggerSelection (i.e. skimming) is enabled, skip this event unless it passes Zorro selection
    if (doSoftwareTriggerSelection && !zorro.isSelected(bc.globalBC())) {
      return;
    }

    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    registry.fill(HIST("h_event_counter"), 0.0); // Coll+TVX+NoTFB+NoITSROFB+Zvtx

    registry.fill(HIST("h_vertexZ"), collision.posZ());

    float rho = collision.rho();

    for (const auto& analysisJet : alljets) {
      if (!isAcceptedJet<FilteredTracks>(analysisJet)) {
        continue;
      }

      fillDataJetHistograms<false>(analysisJet, allTracks, rho);
    }
  }
  PROCESS_SWITCH(BjetTaggingGnn, processDataJetsTrig, "jet information in software triggered Data", false);

  void processDataJetsTrigSub(FilteredCollisionsTriggered::iterator const& collision, FilteredDataJets const& alljets, FilteredTracks const& allTracks, aod::JBCs const& /*bcInfo*/)
  {
    // Get BC info associated with the collision before applying any event selections
    const auto& bc = collision.bc_as<aod::JBCs>();
    // Initialize CCDB objects using the BC info
    initCCDB(bc);
    // If SoftwareTriggerSelection (i.e. skimming) is enabled, skip this event unless it passes Zorro selection
    if (doSoftwareTriggerSelection && !zorro.isSelected(bc.globalBC())) {
      return;
    }

    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    registry.fill(HIST("h_event_counter"), 0.0); // Coll+TVX+NoTFB+NoITSROFB+Zvtx

    registry.fill(HIST("h_vertexZ"), collision.posZ());

    float rho = collision.rho();

    for (const auto& analysisJet : alljets) {
      if (!isAcceptedJet<FilteredTracks>(analysisJet)) {
        continue;
      }

      fillDataJetHistograms<true>(analysisJet, allTracks, rho);
    }
  }
  PROCESS_SWITCH(BjetTaggingGnn, processDataJetsTrigSub, "jet information in software triggered Data, UE(rho*area)-subtracted jet pT", false);

  void processDataJetsSel(AnalysisCollisions::iterator const& collision, FilteredDataJets const& alljets, FilteredTracks const& /*allTracks*/)
  {
    EvtSelFlag evtselCode = fillCollCounter(collision);

    for (const auto& analysisJet : alljets) {
      if (!isAcceptedJet<FilteredTracks>(analysisJet)) {
        continue;
      }

      fillEvtSelStages(HIST("h2_jetpT_evtsel"), analysisJet.pt(), evtselCode, 1.0, kEvtSelStageRecoFirst, kEvtSelStageRecoLast);
    }
  }
  PROCESS_SWITCH(BjetTaggingGnn, processDataJetsSel, "jet information in Data (event selection)", false);

  void processDataJetsSelSub(AnalysisCollisions::iterator const& collision, FilteredDataJets const& alljets, FilteredTracks const& /*allTracks*/)
  {
    EvtSelFlag evtselCode = fillCollCounter(collision);

    float rho = collision.rho();

    for (const auto& analysisJet : alljets) {
      if (!isAcceptedJet<FilteredTracks>(analysisJet)) {
        continue;
      }

      float jetptSub = analysisJet.pt() - rho * analysisJet.area();
      fillEvtSelStages(HIST("h2_jetpT_evtsel_sub"), jetptSub, evtselCode, 1.0, kEvtSelStageRecoFirst, kEvtSelStageRecoLast);
    }
  }
  PROCESS_SWITCH(BjetTaggingGnn, processDataJetsSelSub, "jet information in Data (event selection), UE(rho*area)-subtracted jet pT", false);

  void processDataTracks(FilteredCollisions::iterator const& collision, AnalysisTracks const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }

    for (const auto& track : tracks) {
      if (track.eta() <= trackEtaMin || track.eta() >= trackEtaMax) {
        continue;
      }
      registry.fill(HIST("h_trackpT"), track.pt());
      registry.fill(HIST("h_tracketa"), track.eta());
      registry.fill(HIST("h_trackphi"), track.phi());
      if (track.pt() >= trackPtMin) {
        registry.fill(HIST("h_dcaXY"), std::fabs(track.dcaXY()));
        registry.fill(HIST("h_dcaZ"), std::fabs(track.dcaZ()));
        registry.fill(HIST("hSparse_dca_pt"), track.pt(), std::fabs(track.dcaZ()), std::fabs(track.dcaXY()));
      }
    }
  }
  PROCESS_SWITCH(BjetTaggingGnn, processDataTracks, "track information in Data", false);

  Preslice<FilteredMCPJets> mcpjetsPerMCPCollision = aod::jet::mcCollisionId;

  // Shared body of processMCDJets(Sub)()/processMCDJetsSV(Sub)(): identical except for whether each jet's
  // histograms are filled via fillMCDJetHistograms() or fillMCDJetHistogramsSV(), and whether raw or
  // UE(rho*area)-subtracted jet pT (`withSub`) + the "_sub"-suffixed names are used throughout, including the
  // matched geometric-response and pTHat-sparse histograms below, not just the core Db/flavor spectra.
  // `withSV` selects the branch at compile time via `if constexpr` (bool-template-flag + defaulted-pointer
  // pattern, matches PWGHF/D2H/Tasks/taskCharmPolarisation.cxx's `WithEp`/`QVecs*`).
  template <bool withSV, bool withSub, typename AnalysisJets, typename SecondaryVertices = void>
  void processMCDJetsCommon(FilteredCollisionsMCD::iterator const& collision, AnalysisJets const& MCDjets, FilteredTracksMCD const& allTracks, FilteredMCPJets const& MCPjets, SecondaryVertices const* allSVs = nullptr)
  {
    if (isRejectedMCDCollision(collision)) {
      return;
    }
    bool matchedMcColl = collision.has_mcCollision();

    float rho = collision.rho();
    float weightEvt = useEventWeight ? collision.weight() : 1.f;

    // Keeps hCollCounter populated consistently across every process[Data,MCD]Jets(SV)(Sub) variant, not
    // just the Sel ones - these should all behave identically except for which jet histograms get filled.
    fillCollCounter(collision, weightEvt);
    registry.fill(HIST("h_event_counter"), 0.0, weightEvt);

    registry.fill(HIST("h_vertexZ"), collision.posZ(), weightEvt);

    // Only valid once matchedMcColl is confirmed true below; computed eagerly here (guarded by the ternary)
    // since it's needed inside the MCDjets loop's matched-mcpjet block, not just after it.
    float rhoMc = matchedMcColl ? collision.template mcCollision_as<AnalysisCollisionsMCP>().rho() : 0.f;

    // Store matched particle jet indices to avoid double-counting in mcpjets loop
    std::unordered_set<uint32_t> matchedMcpJetIndices;

    for (const auto& analysisJet : MCDjets) {
      if (!isAcceptedJet<FilteredTracksMCD>(analysisJet)) {
        continue;
      }

      int8_t jetFlavor = 0;
      if constexpr (withSV) {
        jetFlavor = fillMCDJetHistogramsSV<withSub>(analysisJet, allTracks, *allSVs, rho, weightEvt);
      } else {
        jetFlavor = fillMCDJetHistograms<withSub>(analysisJet, allTracks, rho, weightEvt).jetFlavor;
      }

      if (!matchedMcColl) {
        continue;
      }
      float jetpT = withSub ? (analysisJet.pt() - rho * analysisJet.area()) : analysisJet.pt();
      bool matchedJet = false;
      // Everything below is jetpT/mcpJetpT-dependent, so the whole block is duplicated once per `withSub`
      // with the "_sub"-suffixed HIST names.
      for (const auto& mcpjet : analysisJet.template matchedJetGeo_as<FilteredMCPJets>()) {
        // matchedJetGeo_as is not Filtered.
        if (mcpjet.pt() < jetPtMin || mcpjet.pt() >= jetPtMax || mcpjet.eta() >= jetEtaMax - mcpjet.r() / 100.f || mcpjet.eta() <= jetEtaMin + mcpjet.r() / 100.f) {
          continue;
        }
        matchedJet = true;
        matchedMcpJetIndices.insert(mcpjet.globalIndex());
        float mcpJetpT = withSub ? (mcpjet.pt() - rhoMc * mcpjet.area()) : mcpjet.pt();
        if constexpr (withSub) {
          registry.fill(HIST("h2_Response_DetjetpT_PartjetpT_sub"), jetpT, mcpJetpT, weightEvt);
          registry.fill(HIST("h_jetpT_matched_sub"), jetpT, weightEvt);
          registry.fill(HIST("h_jetpT_particle_matched_sub"), mcpJetpT, weightEvt);
          registry.fill(HIST("hSparse_pthat_jetpT_sub"), collision.template mcCollision_as<AnalysisCollisionsMCP>().ptHard(), jetpT, mcpJetpT, weightEvt); // Matched jets
          registry.fill(HIST("h2_jetpT_matched_flavor_sub"), jetpT, getJetFlavorCat(jetFlavor), weightEvt);
          registry.fill(HIST("h2_jetpT_particle_matched_flavor_sub"), mcpJetpT, getJetFlavorCat(jetFlavor), weightEvt);
          if (jetFlavor == JetTaggingSpecies::beauty) {
            registry.fill(HIST("h2_Response_DetjetpT_PartjetpT_b_sub"), jetpT, mcpJetpT, weightEvt);
            registry.fill(HIST("hSparse_pthat_jetpT_b_sub"), collision.template mcCollision_as<AnalysisCollisionsMCP>().ptHard(), jetpT, mcpJetpT, weightEvt); // Matched b-jets
          } else if (jetFlavor == JetTaggingSpecies::charm) {
            registry.fill(HIST("h2_Response_DetjetpT_PartjetpT_c_sub"), jetpT, mcpJetpT, weightEvt);
            registry.fill(HIST("hSparse_pthat_jetpT_c_sub"), collision.template mcCollision_as<AnalysisCollisionsMCP>().ptHard(), jetpT, mcpJetpT, weightEvt); // Matched c-jets
          } else {
            registry.fill(HIST("h2_Response_DetjetpT_PartjetpT_lf_sub"), jetpT, mcpJetpT, weightEvt);
          }
        } else {
          registry.fill(HIST("h2_Response_DetjetpT_PartjetpT"), jetpT, mcpJetpT, weightEvt);
          registry.fill(HIST("h_jetpT_matched"), jetpT, weightEvt);
          registry.fill(HIST("h_jetpT_particle_matched"), mcpJetpT, weightEvt);
          registry.fill(HIST("hSparse_pthat_jetpT"), collision.template mcCollision_as<AnalysisCollisionsMCP>().ptHard(), jetpT, mcpJetpT, weightEvt); // Matched jets
          registry.fill(HIST("h2_jetpT_matched_flavor"), jetpT, getJetFlavorCat(jetFlavor), weightEvt);
          registry.fill(HIST("h2_jetpT_particle_matched_flavor"), mcpJetpT, getJetFlavorCat(jetFlavor), weightEvt);
          if (jetFlavor == JetTaggingSpecies::beauty) {
            registry.fill(HIST("h2_Response_DetjetpT_PartjetpT_b"), jetpT, mcpJetpT, weightEvt);
            registry.fill(HIST("hSparse_pthat_jetpT_b"), collision.template mcCollision_as<AnalysisCollisionsMCP>().ptHard(), jetpT, mcpJetpT, weightEvt); // Matched b-jets
          } else if (jetFlavor == JetTaggingSpecies::charm) {
            registry.fill(HIST("h2_Response_DetjetpT_PartjetpT_c"), jetpT, mcpJetpT, weightEvt);
            registry.fill(HIST("hSparse_pthat_jetpT_c"), collision.template mcCollision_as<AnalysisCollisionsMCP>().ptHard(), jetpT, mcpJetpT, weightEvt); // Matched c-jets
          } else {
            registry.fill(HIST("h2_Response_DetjetpT_PartjetpT_lf"), jetpT, mcpJetpT, weightEvt);
          }
        }
      }
      if (!matchedJet) {
        if constexpr (withSub) {
          registry.fill(HIST("hSparse_pthat_jetpT_sub"), collision.template mcCollision_as<AnalysisCollisionsMCP>().ptHard(), jetpT, -1.f, weightEvt); // Fake jets, overflow-pTpart jets
          if (jetFlavor == JetTaggingSpecies::beauty) {
            registry.fill(HIST("hSparse_pthat_jetpT_b_sub"), collision.template mcCollision_as<AnalysisCollisionsMCP>().ptHard(), jetpT, -1.f, weightEvt); // Overflow-pTpart b-jets
          }
        } else {
          registry.fill(HIST("hSparse_pthat_jetpT"), collision.template mcCollision_as<AnalysisCollisionsMCP>().ptHard(), jetpT, -1.f, weightEvt); // Fake jets, overflow-pTpart jets
          if (jetFlavor == JetTaggingSpecies::beauty) {
            registry.fill(HIST("hSparse_pthat_jetpT_b"), collision.template mcCollision_as<AnalysisCollisionsMCP>().ptHard(), jetpT, -1.f, weightEvt); // Overflow-pTpart b-jets
          }
        }
      }
    }

    if (!matchedMcColl) {
      return;
    }

    // Fill histograms for jets matched to the analysis event selection
    const auto& mcpjetspermcpcollision = MCPjets.sliceBy(mcpjetsPerMCPCollision, collision.mcCollisionId());
    for (const auto& mcpjet : mcpjetspermcpcollision) {
      float mcpJetpT = withSub ? (mcpjet.pt() - rhoMc * mcpjet.area()) : mcpjet.pt();
      int8_t jetFlavor = mcpjet.origin();
      if constexpr (withSub) {
        registry.fill(HIST("h_jetpT_particle_sub"), mcpJetpT, weightEvt);
        // Fill hSparse_pthat_jetpT only for unmatched particle jets (reco pT = -1)
        if (!matchedMcpJetIndices.contains(mcpjet.globalIndex())) {
          registry.fill(HIST("hSparse_pthat_jetpT_sub"), collision.template mcCollision_as<AnalysisCollisionsMCP>().ptHard(), -1.f, mcpJetpT, weightEvt); // Missing jets, overflow-pTreco jets
          if (jetFlavor == JetTaggingSpecies::beauty) {
            registry.fill(HIST("hSparse_pthat_jetpT_b_sub"), collision.template mcCollision_as<AnalysisCollisionsMCP>().ptHard(), -1.f, mcpJetpT, weightEvt); // Missing b-jets, overflow-pTpart b-jets
          } else if (jetFlavor == JetTaggingSpecies::charm) {
            registry.fill(HIST("hSparse_pthat_jetpT_c_sub"), collision.template mcCollision_as<AnalysisCollisionsMCP>().ptHard(), -1.f, mcpJetpT, weightEvt); // Missing c-jets, overflow-pTpart c-jets
          }
        }
        registry.fill(HIST("h2_jetpT_particle_flavor_sub"), mcpJetpT, getJetFlavorCat(jetFlavor), weightEvt);
      } else {
        registry.fill(HIST("h_jetpT_particle"), mcpJetpT, weightEvt);
        // Fill hSparse_pthat_jetpT only for unmatched particle jets (reco pT = -1)
        if (!matchedMcpJetIndices.contains(mcpjet.globalIndex())) {
          registry.fill(HIST("hSparse_pthat_jetpT"), collision.template mcCollision_as<AnalysisCollisionsMCP>().ptHard(), -1.f, mcpJetpT, weightEvt); // Missing jets, overflow-pTreco jets
          if (jetFlavor == JetTaggingSpecies::beauty) {
            registry.fill(HIST("hSparse_pthat_jetpT_b"), collision.template mcCollision_as<AnalysisCollisionsMCP>().ptHard(), -1.f, mcpJetpT, weightEvt); // Missing b-jets, overflow-pTpart b-jets
          } else if (jetFlavor == JetTaggingSpecies::charm) {
            registry.fill(HIST("hSparse_pthat_jetpT_c"), collision.template mcCollision_as<AnalysisCollisionsMCP>().ptHard(), -1.f, mcpJetpT, weightEvt); // Missing c-jets, overflow-pTpart c-jets
          }
        }
        registry.fill(HIST("h2_jetpT_particle_flavor"), mcpJetpT, getJetFlavorCat(jetFlavor), weightEvt);
      }
    }
  }

  void processMCDJets(FilteredCollisionsMCD::iterator const& collision, FilteredMCDJets const& MCDjets, FilteredTracksMCD const& allTracks, FilteredMCPJets const& MCPjets, aod::JetParticles const& /*mcParticles*/, AnalysisCollisionsMCP const& /*mcCollisions*/)
  {
    processMCDJetsCommon<false, false>(collision, MCDjets, allTracks, MCPjets);
  }
  PROCESS_SWITCH(BjetTaggingGnn, processMCDJets, "jet information in MC", false);

  void processMCDJetsSub(FilteredCollisionsMCD::iterator const& collision, FilteredMCDJets const& MCDjets, FilteredTracksMCD const& allTracks, FilteredMCPJets const& MCPjets, aod::JetParticles const& /*mcParticles*/, AnalysisCollisionsMCP const& /*mcCollisions*/)
  {
    processMCDJetsCommon<false, true>(collision, MCDjets, allTracks, MCPjets);
  }
  PROCESS_SWITCH(BjetTaggingGnn, processMCDJetsSub, "jet information in MC, UE(rho*area)-subtracted jet pT", false);

  void processMCDJetsSV(FilteredCollisionsMCD::iterator const& collision, FilteredMCDJetsSV const& MCDjets, aod::MCDSecondaryVertex3Prongs const& allSVs, FilteredTracksMCD const& allTracks, FilteredMCPJets const& MCPjets, aod::JetParticles const& /*mcParticles*/, AnalysisCollisionsMCP const& /*mcCollisions*/)
  {
    processMCDJetsCommon<true, false>(collision, MCDjets, allTracks, MCPjets, &allSVs);
  }
  PROCESS_SWITCH(BjetTaggingGnn, processMCDJetsSV, "jet information in MC with secondary vertex info", false);

  void processMCDJetsSVSub(FilteredCollisionsMCD::iterator const& collision, FilteredMCDJetsSV const& MCDjets, aod::MCDSecondaryVertex3Prongs const& allSVs, FilteredTracksMCD const& allTracks, FilteredMCPJets const& MCPjets, aod::JetParticles const& /*mcParticles*/, AnalysisCollisionsMCP const& /*mcCollisions*/)
  {
    processMCDJetsCommon<true, true>(collision, MCDjets, allTracks, MCPjets, &allSVs);
  }
  PROCESS_SWITCH(BjetTaggingGnn, processMCDJetsSVSub, "jet information in MC with secondary vertex info, UE(rho*area)-subtracted jet pT", false);

  void processMCDJetsSel(AnalysisCollisionsMCD::iterator const& collision, FilteredMCDJets const& MCDjets, FilteredTracksMCD const& /*allTracks*/, FilteredMCPJets const& /*MCPjets*/, AnalysisCollisionsMCP const& /*mcCollisions*/)
  {
    // Reject outlier MC collisions
    if (collision.isOutlier()) {
      return;
    }
    if (collision.has_mcCollision() && collision.template mcCollision_as<AnalysisCollisionsMCP>().isOutlier()) {
      return;
    }

    float weightEvt = useEventWeight ? collision.weight() : 1.f;
    EvtSelFlag evtselCode = fillCollCounter(collision, weightEvt);

    for (const auto& analysisJet : MCDjets) {
      if (!isAcceptedJet<FilteredTracksMCD>(analysisJet)) {
        continue;
      }

      int8_t jetFlavor = analysisJet.origin();

      fillEvtSelStages(HIST("h2_jetpT_evtsel"), analysisJet.pt(), evtselCode, weightEvt, kEvtSelStageRecoFirst, kEvtSelStageRecoLast);
      if (jetFlavor == JetTaggingSpecies::beauty) {
        fillEvtSelStages(HIST("h2_jetpT_evtsel_b"), analysisJet.pt(), evtselCode, weightEvt, kEvtSelStageRecoFirst, kEvtSelStageRecoLast);
      } else if (jetFlavor == JetTaggingSpecies::charm) {
        fillEvtSelStages(HIST("h2_jetpT_evtsel_c"), analysisJet.pt(), evtselCode, weightEvt, kEvtSelStageRecoFirst, kEvtSelStageRecoLast);
      }
    }
  }
  PROCESS_SWITCH(BjetTaggingGnn, processMCDJetsSel, "jet information in MC (event selection)", false);

  void processMCDJetsSelSub(AnalysisCollisionsMCD::iterator const& collision, FilteredMCDJets const& MCDjets, FilteredTracksMCD const& /*allTracks*/, FilteredMCPJets const& /*MCPjets*/, AnalysisCollisionsMCP const& /*mcCollisions*/)
  {
    // Reject outlier MC collisions
    if (collision.isOutlier()) {
      return;
    }
    if (collision.has_mcCollision() && collision.template mcCollision_as<AnalysisCollisionsMCP>().isOutlier()) {
      return;
    }

    float weightEvt = useEventWeight ? collision.weight() : 1.f;
    EvtSelFlag evtselCode = fillCollCounter(collision, weightEvt);

    float rho = collision.rho();

    for (const auto& analysisJet : MCDjets) {
      if (!isAcceptedJet<FilteredTracksMCD>(analysisJet)) {
        continue;
      }

      int8_t jetFlavor = analysisJet.origin();
      float jetptSub = analysisJet.pt() - rho * analysisJet.area();

      fillEvtSelStages(HIST("h2_jetpT_evtsel_sub"), jetptSub, evtselCode, weightEvt, kEvtSelStageRecoFirst, kEvtSelStageRecoLast);
      if (jetFlavor == JetTaggingSpecies::beauty) {
        fillEvtSelStages(HIST("h2_jetpT_evtsel_b_sub"), jetptSub, evtselCode, weightEvt, kEvtSelStageRecoFirst, kEvtSelStageRecoLast);
      } else if (jetFlavor == JetTaggingSpecies::charm) {
        fillEvtSelStages(HIST("h2_jetpT_evtsel_c_sub"), jetptSub, evtselCode, weightEvt, kEvtSelStageRecoFirst, kEvtSelStageRecoLast);
      }
    }
  }
  PROCESS_SWITCH(BjetTaggingGnn, processMCDJetsSelSub, "jet information in MC (event selection), UE(rho*area)-subtracted jet pT", false);

  PresliceUnsorted<AnalysisCollisionsMCD> collisionsPerMCPCollision = aod::jmccollisionlb::mcCollisionId;

  // `withSub` selects UE(rho*area)-subtracted mcp-jet pT + the "_sub"-suffixed evtsel histogram names below
  // instead of raw pT + the unsuffixed names. hMcCollCounter/h_vertexZ_truth*/h_event_counter_mcp aren't
  // jet-pT-dependent and always use the plain name.
  template <bool withSub>
  void processMCPJetsCommon(AnalysisCollisionsMCP const& mcCollisions, FilteredMCPJets const& mcpjets, AnalysisCollisionsMCD const& collisions, aod::JBCs const&)
  {
    // Subscribing AnalysisCollisionsMCP::iterator causes an issue related to unsorted JMcCollisionLbs index.
    for (const auto& mcCollision : mcCollisions) {
      if (mcCollision.isOutlier()) {
        continue;
      }
      const auto& matchedCollisions = collisions.sliceBy(collisionsPerMCPCollision, mcCollision.mcCollisionId());
      if (matchedCollisions.size() >= 1) {
        if (matchedCollisions.begin().isOutlier()) {
          continue;
        }
      }
      // mcCollision -> BC is a mandatory 1:1 index (jmccollision::JBCId), unlike the reco-collision match
      // above (which can be 0, 1, or split into several matchedCollisions) - no slicing/size check needed.
      const auto& matchedBC = mcCollision.bc_as<aod::JBCs>();

      float rho = mcCollision.rho();
      float weightEvt = useEventWeight ? mcCollision.weight() : 1.f;

      EvtSelFlag evtselCode = EvtSelFlag::INEL;
      EvtSelFlag evtselBCCode = EvtSelFlag::INEL;
      registry.fill(HIST("hMcCollCounter"), static_cast<int>(EvtSel::INEL), weightEvt); // INEL
      registry.fill(HIST("h_vertexZ_truth"), mcCollision.posZ(), weightEvt);

      bool zvtx = std::fabs(mcCollision.posZ()) < vertexZCut;
      bool zvtxMatched = false;

      // bool isTrueINELgt0 = isTrueINEL0(mcCollision, mcParticles);
      if (zvtx) {
        evtselCode |= EvtSelFlag::kZvtx;
        registry.fill(HIST("hMcCollCounter"), static_cast<int>(EvtSel::INELZvtx), weightEvt); // INEL+Zvtx
        // if (isTrueINELgt0) {
        //   evtselCode |= EvtSelFlag::kINELgt0;
        //   registry.fill(HIST("hMcCollCounter"), static_cast<int>(EvtSel::INELgt0), weightEvt); // INELgt0
        // }
      }

      // Coll
      if (matchedCollisions.size() >= 1) {
        zvtxMatched = std::fabs(matchedCollisions.begin().posZ()) < vertexZCut;
        evtselCode |= EvtSelFlag::kColl;
        registry.fill(HIST("h_vertexZ_truth_coll"), mcCollision.posZ(), weightEvt);
        registry.fill(HIST("hMcCollCounter"), static_cast<int>(EvtSel::Coll), weightEvt); // McColl(-> Coll)
        if (zvtxMatched) {
          registry.fill(HIST("hMcCollCounter"), static_cast<int>(EvtSel::CollZvtx), weightEvt); // McColl(-> Coll+Zvtx)
        }
        if (jetderiveddatautilities::selectCollision(matchedCollisions.begin(), eventSelectionBitsTVX)) {
          evtselCode |= EvtSelFlag::kTVX;
          registry.fill(HIST("h_vertexZ_truth_tvx"), mcCollision.posZ(), weightEvt);
          registry.fill(HIST("hMcCollCounter"), static_cast<int>(EvtSel::TVX), weightEvt); // McColl(-> Coll+TVX)
          if (zvtxMatched) {
            registry.fill(HIST("hMcCollCounter"), static_cast<int>(EvtSel::TVXZvtx), weightEvt); // McColl(-> Coll+TVX+Zvtx)
          }
          if (jetderiveddatautilities::selectCollision(matchedCollisions.begin(), eventSelectionBitsSelMC)) {
            evtselCode |= EvtSelFlag::kNoTFB;
            registry.fill(HIST("hMcCollCounter"), static_cast<int>(EvtSel::SelMC), weightEvt); // McColl(-> Coll+TVX+NoTFB)
            if (zvtxMatched) {
              registry.fill(HIST("hMcCollCounter"), static_cast<int>(EvtSel::SelMCZvtx), weightEvt); // McColl(-> Coll+TVX+NoTFB+Zvtx)
            }
            if (jetderiveddatautilities::selectCollision(matchedCollisions.begin(), eventSelectionBitsSel8)) {
              evtselCode |= EvtSelFlag::kNoITSROFB;
              registry.fill(HIST("hMcCollCounter"), static_cast<int>(EvtSel::Sel8), weightEvt); // McColl(-> Coll+TVX+NoTFB+NoITSROFB)
              if (zvtxMatched) {
                registry.fill(HIST("hMcCollCounter"), static_cast<int>(EvtSel::Sel8Zvtx), weightEvt); // McColl(-> Coll+TVX+NoTFB+NoITSROFB+Zvtx)
                // if (matchedCollisions.begin().template collision_as<OrigCollisions>().isInelGt0()) {
                //   evtselCode |= EvtSelFlag::kINELgt0rec;
                //   registry.fill(HIST("hMcCollCounter"), static_cast<int>(EvtSel::INELgt0rec), weightEvt); // INELgt0+Zvtx(rec)
                // }
              }
              if (jetderiveddatautilities::selectCollision(matchedCollisions.begin(), eventSelectionBitsSel8Full)) {
                evtselCode |= EvtSelFlag::kNoPileup;
                registry.fill(HIST("hMcCollCounter"), static_cast<int>(EvtSel::Sel8Full), weightEvt); // McColl(-> Coll+TVX+NoTFB+NoITSROFB+NoPileup)
                if (zvtxMatched) {
                  registry.fill(HIST("hMcCollCounter"), static_cast<int>(EvtSel::Sel8FullZvtx), weightEvt); // McColl(-> Coll+TVX+NoTFB+NoITSROFB+NoPileup+Zvtx)
                }
                if (jetderiveddatautilities::selectCollision(matchedCollisions.begin(), eventSelectionBitsSel8FullGood)) {
                  evtselCode |= EvtSelFlag::kIsGoodZvtx;
                  registry.fill(HIST("hMcCollCounter"), static_cast<int>(EvtSel::Sel8FullGood), weightEvt); // McColl(-> Coll+TVX+NoTFB+NoITSROFB+NoPileup+GoodZvtx)
                  if (zvtxMatched) {
                    registry.fill(HIST("hMcCollCounter"), static_cast<int>(EvtSel::Sel8FullGoodZvtx), weightEvt); // McColl(-> Coll+TVX+NoTFB+NoITSROFB+NoPileup+GoodZvtx+Zvtx)
                  }
                }
              }
            }
          }
        }
        if (jetderiveddatautilities::selectCollision(matchedCollisions.begin(), eventSelectionBits) && zvtxMatched) {
          registry.fill(HIST("h_event_counter_mcp"), 0.0, weightEvt);
        }
      }

      // BC (mcCollision's own truth-level BC, always resolved - see fillMcCollCounterBC()).
      evtselBCCode |= fillMcCollCounterBC(matchedBC, mcCollision.posZ(), weightEvt);

      const auto& mcpjetspermcpcollision = mcpjets.sliceBy(mcpjetsPerMCPCollision, mcCollision.mcCollisionId());
      for (const auto& mcpjet : mcpjetspermcpcollision) {

        int8_t jetFlavor = mcpjet.origin();
        float mcpJetpT = withSub ? (mcpjet.pt() - rho * mcpjet.area()) : mcpjet.pt();

        if constexpr (withSub) {
          fillEvtSelStages(HIST("h2_jetpT_particle_evtsel_sub"), mcpJetpT, evtselCode, weightEvt, 0, kEvtSelStages.size() - 1);
          if (jetFlavor == JetTaggingSpecies::beauty) {
            fillEvtSelStages(HIST("h2_jetpT_particle_evtsel_b_sub"), mcpJetpT, evtselCode, weightEvt, 0, kEvtSelStages.size() - 1);
          } else if (jetFlavor == JetTaggingSpecies::charm) {
            fillEvtSelStages(HIST("h2_jetpT_particle_evtsel_c_sub"), mcpJetpT, evtselCode, weightEvt, 0, kEvtSelStages.size() - 1);
          }

          fillEvtSelStages(HIST("h2_jetpT_particle_evtselbc_sub"), mcpJetpT, evtselBCCode, weightEvt, 0, kEvtSelStages.size() - 1);
          if (jetFlavor == JetTaggingSpecies::beauty) {
            fillEvtSelStages(HIST("h2_jetpT_particle_evtselbc_b_sub"), mcpJetpT, evtselBCCode, weightEvt, 0, kEvtSelStages.size() - 1);
          } else if (jetFlavor == JetTaggingSpecies::charm) {
            fillEvtSelStages(HIST("h2_jetpT_particle_evtselbc_c_sub"), mcpJetpT, evtselBCCode, weightEvt, 0, kEvtSelStages.size() - 1);
          }
        } else {
          fillEvtSelStages(HIST("h2_jetpT_particle_evtsel"), mcpJetpT, evtselCode, weightEvt, 0, kEvtSelStages.size() - 1);
          if (jetFlavor == JetTaggingSpecies::beauty) {
            fillEvtSelStages(HIST("h2_jetpT_particle_evtsel_b"), mcpJetpT, evtselCode, weightEvt, 0, kEvtSelStages.size() - 1);
          } else if (jetFlavor == JetTaggingSpecies::charm) {
            fillEvtSelStages(HIST("h2_jetpT_particle_evtsel_c"), mcpJetpT, evtselCode, weightEvt, 0, kEvtSelStages.size() - 1);
          }

          fillEvtSelStages(HIST("h2_jetpT_particle_evtselbc"), mcpJetpT, evtselBCCode, weightEvt, 0, kEvtSelStages.size() - 1);
          if (jetFlavor == JetTaggingSpecies::beauty) {
            fillEvtSelStages(HIST("h2_jetpT_particle_evtselbc_b"), mcpJetpT, evtselBCCode, weightEvt, 0, kEvtSelStages.size() - 1);
          } else if (jetFlavor == JetTaggingSpecies::charm) {
            fillEvtSelStages(HIST("h2_jetpT_particle_evtselbc_c"), mcpJetpT, evtselBCCode, weightEvt, 0, kEvtSelStages.size() - 1);
          }
        }
      }
    }
  }

  void processMCPJets(AnalysisCollisionsMCP const& mcCollisions, FilteredMCPJets const& mcpjets, AnalysisCollisionsMCD const& collisions, aod::JBCs const& bcs)
  {
    processMCPJetsCommon<false>(mcCollisions, mcpjets, collisions, bcs);
  }
  PROCESS_SWITCH(BjetTaggingGnn, processMCPJets, "mc collision information", false);

  void processMCPJetsSub(AnalysisCollisionsMCP const& mcCollisions, FilteredMCPJets const& mcpjets, AnalysisCollisionsMCD const& collisions, aod::JBCs const& bcs)
  {
    processMCPJetsCommon<true>(mcCollisions, mcpjets, collisions, bcs);
  }
  PROCESS_SWITCH(BjetTaggingGnn, processMCPJetsSub, "mc collision information, UE(rho*area)-subtracted jet pT", false);

  Preslice<aod::JetParticles> mcparticlesPerMCPCollision = aod::jmcparticle::mcCollisionId;

  void processMCDTracks(FilteredCollisionsMCD::iterator const& collision, AnalysisTracksMCD const& tracks, AnalysisCollisionsMCP const& /*mcCollisions*/, aod::JetParticles const& allParticles)
  {
    if (isRejectedMCDCollision(collision)) {
      return;
    }
    bool matchedMcColl = collision.has_mcCollision();

    float weightEvt = useEventWeight ? collision.weight() : 1.f;

    for (const auto& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelectionBits) || track.eta() <= trackEtaMin || track.eta() >= trackEtaMax) {
        continue;
      }
      registry.fill(HIST("h_trackpT"), track.pt(), weightEvt);
      registry.fill(HIST("h_tracketa"), track.eta(), weightEvt);
      registry.fill(HIST("h_trackphi"), track.phi(), weightEvt);
      if (track.pt() >= trackPtMin) {
        registry.fill(HIST("h_dcaXY"), std::fabs(track.dcaXY()), weightEvt);
        registry.fill(HIST("h_dcaZ"), std::fabs(track.dcaZ()), weightEvt);
        registry.fill(HIST("hSparse_dca_pt"), track.pt(), std::fabs(track.dcaZ()), std::fabs(track.dcaXY()), weightEvt);
      }

      if (!matchedMcColl) {
        if (track.pt() >= trackPtMin) {
          registry.fill(HIST("h_dcaXY_coll_fake"), std::fabs(track.dcaXY()), weightEvt);
          registry.fill(HIST("h_dcaZ_coll_fake"), std::fabs(track.dcaZ()), weightEvt);
          registry.fill(HIST("hSparse_dca_pt_coll_fake"), track.pt(), std::fabs(track.dcaZ()), std::fabs(track.dcaXY()), weightEvt);
        }
        continue;
      }
      if (!track.has_mcParticle()) {
        if (track.pt() >= trackPtMin) {
          registry.fill(HIST("h_dcaXY_fake"), std::fabs(track.dcaXY()), weightEvt);
          registry.fill(HIST("h_dcaZ_fake"), std::fabs(track.dcaZ()), weightEvt);
          registry.fill(HIST("hSparse_dca_pt_fake"), track.pt(), std::fabs(track.dcaZ()), std::fabs(track.dcaXY()), weightEvt);
        }
        continue;
      }
      const auto& particle = track.template mcParticle_as<aod::JetParticles>();
      if (particle.eta() > trackEtaMin && particle.eta() < trackEtaMax) {
        if (particle.isPhysicalPrimary()) {
          registry.fill(HIST("h2_trackpT_partpT"), track.pt(), particle.pt(), weightEvt);
          registry.fill(HIST("h_partpT_matched_fine"), particle.pt(), weightEvt);
        }
        if (track.pt() >= trackPtMin) {
          // Track association accuracy as a function of DCA
          if (particle.isPhysicalPrimary()) {
            if (particle.mcCollisionId() == collision.mcCollisionId()) {
              registry.fill(HIST("h_dcaXY_coll_matched"), std::fabs(track.dcaXY()), weightEvt); // Matched to particle from the same MC collision
              registry.fill(HIST("h_dcaZ_coll_matched"), std::fabs(track.dcaZ()), weightEvt);
              registry.fill(HIST("hSparse_dca_pt_coll_matched"), track.pt(), std::fabs(track.dcaZ()), std::fabs(track.dcaXY()), weightEvt);
              int origin = RecoDecay::getParticleOrigin(allParticles, particle, false);
              if (origin == RecoDecay::OriginType::NonPrompt) {
                registry.fill(HIST("h_dcaXY_coll_matched_b"), std::fabs(track.dcaXY()), weightEvt);
                registry.fill(HIST("h_dcaZ_coll_matched_b"), std::fabs(track.dcaZ()), weightEvt);
                registry.fill(HIST("hSparse_dca_pt_coll_matched_b"), track.pt(), std::fabs(track.dcaZ()), std::fabs(track.dcaXY()), weightEvt);
              } else if (origin == RecoDecay::OriginType::Prompt) {
                registry.fill(HIST("h_dcaXY_coll_matched_c"), std::fabs(track.dcaXY()), weightEvt);
                registry.fill(HIST("h_dcaZ_coll_matched_c"), std::fabs(track.dcaZ()), weightEvt);
                registry.fill(HIST("hSparse_dca_pt_coll_matched_c"), track.pt(), std::fabs(track.dcaZ()), std::fabs(track.dcaXY()), weightEvt);
              } else {
                registry.fill(HIST("h_dcaXY_coll_matched_lf"), std::fabs(track.dcaXY()), weightEvt);
                registry.fill(HIST("h_dcaZ_coll_matched_lf"), std::fabs(track.dcaZ()), weightEvt);
                registry.fill(HIST("hSparse_dca_pt_coll_matched_lf"), track.pt(), std::fabs(track.dcaZ()), std::fabs(track.dcaXY()), weightEvt);
              }
            } else {
              registry.fill(HIST("h_dcaXY_coll_mismatched"), std::fabs(track.dcaXY()), weightEvt); // Matched to particle from a different MC collision
              registry.fill(HIST("h_dcaZ_coll_mismatched"), std::fabs(track.dcaZ()), weightEvt);
              registry.fill(HIST("hSparse_dca_pt_coll_mismatched"), track.pt(), std::fabs(track.dcaZ()), std::fabs(track.dcaXY()), weightEvt);
            }
          } else {
            if (particle.mcCollisionId() == collision.mcCollisionId()) {
              registry.fill(HIST("h_dcaXY_npp"), std::fabs(track.dcaXY()), weightEvt);
              registry.fill(HIST("h_dcaZ_npp"), std::fabs(track.dcaZ()), weightEvt);
              registry.fill(HIST("hSparse_dca_pt_npp"), track.pt(), std::fabs(track.dcaZ()), std::fabs(track.dcaXY()), weightEvt);
            } else {
              registry.fill(HIST("h_dcaXY_npp_mismatched"), std::fabs(track.dcaXY()), weightEvt);
              registry.fill(HIST("h_dcaZ_npp_mismatched"), std::fabs(track.dcaZ()), weightEvt);
              registry.fill(HIST("hSparse_dca_pt_npp_mismatched"), track.pt(), std::fabs(track.dcaZ()), std::fabs(track.dcaXY()), weightEvt);
            }
          }
        }
      }
    }

    if (!matchedMcColl) {
      return;
    }

    const auto& particles = allParticles.sliceBy(mcparticlesPerMCPCollision, collision.mcCollisionId());

    for (const auto& particle : particles) {
      const auto pdgParticle = pdg->GetParticle(particle.pdgCode());
      if (!pdgParticle || pdgParticle->Charge() == 0.0) {
        continue;
      }
      if (particle.isPhysicalPrimary() && particle.eta() > trackEtaMin && particle.eta() < trackEtaMax) {
        registry.fill(HIST("h_partpT"), particle.pt(), weightEvt);
        registry.fill(HIST("h_partpT_fine"), particle.pt(), weightEvt);
      }
    }
  }
  PROCESS_SWITCH(BjetTaggingGnn, processMCDTracks, "track information in MCD", false);

  PresliceUnsorted<o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>> perFoundBC = aod::evsel::foundBCId;

  void processBCs(soa::Join<aod::BCs, aod::BcSels> const& bcs, soa::Join<aod::Collisions, aod::EvSels> const& collisions)
  {
    if (bcs.size() == 0) {
      return;
    }
    for (const auto& bc : bcs) {
      registry.fill(HIST("hBCCounter"), 0.5); // All BC
      if (bc.selection_bit(aod::evsel::kIsTriggerTVX)) {
        registry.fill(HIST("hBCCounter"), 1.5); // BC+TVX
        if (bc.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
          registry.fill(HIST("hBCCounter"), 2.5); // BC+TVX+NoTFB
          if (bc.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
            registry.fill(HIST("hBCCounter"), 3.5); // BC+TVX+NoTFB+NoITSROFB ----> this goes to Lumi i.e. hLumiAfterBCcuts in eventSelection task
          }
        }
      }
      const auto& collisionsInBC = collisions.sliceBy(perFoundBC, bc.globalIndex());
      for (const auto& collision : collisionsInBC) {
        registry.fill(HIST("hBCCounter"), 4.5); // CollinBC
        if (collision.sel8()) {
          registry.fill(HIST("hBCCounter"), 5.5); // CollinBC+sel8
          if (std::fabs(collision.posZ()) < vertexZCut) {
            registry.fill(HIST("hBCCounter"), 6.5); // CollinBC+sel8+VtxZ
          }
          if (collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
            registry.fill(HIST("hBCCounter"), 7.5); // CollinBC+sel8Full
            if (collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
              registry.fill(HIST("hBCCounter"), 8.5); // CollinBC+sel8Full+GoodZvtx
              if (std::fabs(collision.posZ()) < vertexZCut) {
                registry.fill(HIST("hBCCounter"), 9.5); // CollinBC+sel8Full+VtxZ+GoodZvtx ----> this goes to my analysis task for jet events selection
              }
            }
          }
        }
      } // collision loop
    } // bc loop
  }
  PROCESS_SWITCH(BjetTaggingGnn, processBCs, "BCs for 0 vertex QA", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<BjetTaggingGnn>(cfgc)};
}
