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
#include "PWGJE/DataModel/JetTagging.h"

#include "Common/CCDB/TriggerAliases.h"
#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/DataModel/EventSelection.h"

#include "CCDB/BasicCCDBManager.h"
#include "Framework/ASoA.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include <Framework/Configurable.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <type_traits>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace BjetTaggingGnnEvtSel
{
enum class EvtSelFlag : uint8_t {
  kNone = 0,
  kINEL = 1 << 0,
  kColl = 1 << 1,
  kTVX = 1 << 2,
  kNoTFB = 1 << 3,
  kNoITSROFB = 1 << 4,
  kZvtx = 1 << 5,
  kINELgt0 = 1 << 6,
  kINELgt0rec = 1 << 7,

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
  INELgt0 = kINEL | kZvtx | kINELgt0,
  INELgt0rec = kINEL | kZvtx | kColl | kTVX | kNoTFB | kNoITSROFB | kINELgt0rec
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
  INELgt0,
  INELgt0rec
};
}; // namespace BjetTaggingGnnEvtSel
using namespace BjetTaggingGnnEvtSel;

struct BjetTaggingGnn {

  HistogramRegistry registry;

  // event level configurables
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<bool> useEventWeight{"useEventWeight", true, "Flag whether to scale histograms with the event weight"};

  Configurable<float> pTHatMaxMCD{"pTHatMaxMCD", 999.0, "maximum fraction of hard scattering for jet acceptance in detector MC"};
  Configurable<float> pTHatMaxMCP{"pTHatMaxMCP", 999.0, "maximum fraction of hard scattering for jet acceptance in particle MC"};
  Configurable<float> pTHatExponent{"pTHatExponent", 6.0, "exponent of the event weight for the calculation of pTHat"};

  // track level configurables
  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum track pT"};
  Configurable<float> trackPtMax{"trackPtMax", 1000.0, "maximum track pT"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum track eta"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum track eta"};
  Configurable<float> trackPtMinGnn{"trackPtMinGnn", 0.5, "minimum track pT for GNN inputs"};

  Configurable<float> maxIPxy{"maxIPxy", 10, "maximum track DCA in xy plane"};
  Configurable<float> maxIPz{"maxIPz", 10, "maximum track DCA in z direction"};

  Configurable<float> trackNppCrit{"trackNppCrit", 0.95, "track not physical primary ratio"};

  // jet level configurables
  Configurable<float> jetPtMin{"jetPtMin", 5.0, "minimum jet pT"};
  Configurable<float> jetPtMax{"jetPtMax", 1000.0, "maximum jet pT"};
  Configurable<float> jetEtaMin{"jetEtaMin", -99.0, "minimum jet pseudorapidity"};
  Configurable<float> jetEtaMax{"jetEtaMax", 99.0, "maximum jet pseudorapidity"};
  Configurable<float> leadingConstituentPtMin{"leadingConstituentPtMin", -99.0, "minimum pT selection on jet constituent"};
  Configurable<float> leadingConstituentPtMax{"leadingConstituentPtMax", 9999.0, "maximum pT selection on jet constituent"};
  Configurable<float> jetAreaFractionMin{"jetAreaFractionMin", -99.0, "used to make a cut on the jet areas"};

  Configurable<std::vector<double>> jetRadii{"jetRadii", std::vector<double>{0.4}, "jet resolution parameters"};

  Configurable<double> dbMin{"dbMin", -10., "minimum GNN Db"};
  Configurable<double> dbMax{"dbMax", 20., "maximum GNN Db"};
  Configurable<int> dbNbins{"dbNbins", 3000, "number of bins in axisDbFine"};

  Configurable<bool> doDataDriven{"doDataDriven", false, "Flag whether to use fill THnSpase for data driven methods"};
  Configurable<bool> doDataDrivenExtra{"doDataDrivenExtra", false, "Flag whether to add extra axes to THnSparses"};
  Configurable<bool> callSumw2{"callSumw2", false, "Flag whether to call THnSparse::Sumw2() for error calculation"};

  Configurable<int> trainingDatasetRatioParam{"trainingDatasetRatioParam", 0, "Parameter for splitting training/evaluation datasets by collisionId"};

  // Software Trigger configurables
  Configurable<bool> doSoftwareTriggerSelection{"doSoftwareTriggerSelection", false, "set to true when using triggered datasets"};
  Configurable<std::string> triggerMasks{"triggerMasks", "fJetChLowPt", "possible JE Trigger masks: fJetChLowPt,fJetChHighPt,fEMCALReadout,fJetFullHighPt,fJetFullLowPt"};

  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  // Service
  Service<o2::framework::O2DatabasePDG> pdg;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // Event selection bits
  std::vector<int> eventSelectionBits;
  std::vector<int> eventSelectionBitsTVX;
  std::vector<int> eventSelectionBitsSelMC;
  std::vector<int> eventSelectionBitsSel8;

  std::vector<double> jetRadiiValues;

  void init(InitContext const&)
  {
    jetRadiiValues = (std::vector<double>)jetRadii;

    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    eventSelectionBitsTVX = jetderiveddatautilities::initialiseEventSelectionBits("TVX");
    eventSelectionBitsSel8 = jetderiveddatautilities::initialiseEventSelectionBits("sel8");
    eventSelectionBitsSelMC = jetderiveddatautilities::initialiseEventSelectionBits("selMC");

    if (doprocessDataJetsTrig) {
      zorroSummary.setObject(zorro.getZorroSummary());
    }

    registry.add("h_event_counter", ";analysis collision counter", {HistType::kTH1F, {{1, 0.0, 1.0}}}, callSumw2);
    registry.add("h_event_counter_mcp", ";analysis collision matched MC collision counter", {HistType::kTH1F, {{1, 0.0, 1.0}}}, callSumw2);
    registry.add("h_vertexZ", "Vertex Z;#it{Z} (cm)", {HistType::kTH1F, {{100, -20.0, 20.0}}}, callSumw2);
    registry.add("hCollCounter", ";collision counter", {HistType::kTH1F, {{12, 1.0, 13.0}}}, callSumw2);
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
    hCollCounter->GetXaxis()->SetBinLabel(11, "_11");
    hCollCounter->GetXaxis()->SetBinLabel(12, "INELgt0+Zvtx(rec)"); // sel8
    registry.add("hMcCollCounter", ";MC collision counter", {HistType::kTH1F, {{12, 1.0, 13.0}}}, callSumw2);
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
    hMcCollCounter->GetXaxis()->SetBinLabel(11, "INELgt0+Zvtx");
    hMcCollCounter->GetXaxis()->SetBinLabel(12, "INELgt0+Zvtx(-> INELgt0+Zvtx(rec))"); // sel8
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
    const AxisSpec axisJetpT{200, 0., 200., "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisJetEta{200, -0.8, 0.8, "#it{#eta}_{jet}"};
    const AxisSpec axisDb{200, dbMin, dbMax, "#it{D}_{b}"};
    const AxisSpec axisDbFine{dbNbins, dbMin, dbMax, "#it{D}_{b}"};
    const AxisSpec axisJetMass{200, 0., 50., "#it{m}_{jet} (GeV/#it{c}^{2})"};
    const AxisSpec axisJetProb{200, 0., 40., "-ln(JP)"};
    const AxisSpec axisNTracks{42, 0, 42, "#it{n}_{tracks}"};

    registry.add("h_jetpT", "", {HistType::kTH1F, {axisJetpT}}, callSumw2);
    registry.add("h_jetEta", "", {HistType::kTH1F, {axisJetEta}}, callSumw2);
    registry.add("h_jetPhi", "", {HistType::kTH1F, {{200, 0., 2. * M_PI, "#it{phi}_{jet}"}}});
    registry.add("h_jetMass", "", {HistType::kTH1F, {axisJetMass}});
    registry.add("h_Db", "", {HistType::kTH1F, {axisDbFine}});
    registry.add("h2_jetpT_Db", "", {HistType::kTH2F, {axisJetpT, axisDb}});

    registry.add("h_gnnfeat_trackpT", "", {HistType::kTH1F, {{200, 0., 100., "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("h_gnnfeat_trackPhi", "", {HistType::kTH1F, {{200, 0., 2. * M_PI, "#it{#phi}"}}});
    registry.add("h_gnnfeat_trackEta", "", {HistType::kTH1F, {{200, -0.9, 0.9, "#it{#eta}"}}});
    registry.add("h_gnnfeat_trackCharge", "", {HistType::kTH1F, {{3, -1., 2., "#it{q}"}}});
    registry.add("h_gnnfeat_trackDCAxy", "", {HistType::kTH1F, {{200, -5., 5., "DCA_#it{xy} (cm)"}}});
    registry.add("h_gnnfeat_trackSigmaDCAxy", "", {HistType::kTH1F, {{200, 0., 5., "#it{#sigma}_{{DCA_#it{xy}} (cm)"}}});
    registry.add("h_gnnfeat_trackDCAz", "", {HistType::kTH1F, {{200, -5., 5., "DCA_#it{z} (cm)"}}});
    registry.add("h_gnnfeat_trackSigmaDCAz", "", {HistType::kTH1F, {{200, 0., 5., "#it{#sigma}_{{DCA_#it{z}} (cm)"}}});
    registry.add("h_gnnfeat_trackITSChi2NCl", "", {HistType::kTH1F, {{200, 0., 40., "ITS #it{#chi}^{2}/ndf"}}});
    registry.add("h_gnnfeat_trackTPCChi2NCl", "", {HistType::kTH1F, {{200, 0., 5., "TPC #it{#chi}^{2}/ndf"}}});
    registry.add("h_gnnfeat_trackITSNCls", "", {HistType::kTH1F, {{8, 0., 8., "ITS NCls"}}});
    registry.add("h_gnnfeat_trackTPCNCls", "", {HistType::kTH1F, {{153, 0., 153., "TPC NCls"}}});
    registry.add("h_gnnfeat_trackTPCNCrossedRows", "", {HistType::kTH1F, {{153, 0., 153., "TPC NCrossedRows"}}});
    registry.add("h_gnnfeat_tracksIPxy", "", {HistType::kTH1F, {{200, -5., 5., "{sIP}_#it{xy}"}}});
    registry.add("h_gnnfeat_tracksIPz", "", {HistType::kTH1F, {{200, -5., 5., "{sIP}_#it{z}"}}});

    if (doprocessDataTracks || doprocessMCDTracks) {
      registry.add("h_trackpT", "", {HistType::kTH1F, {axisTrackpT}}, callSumw2);
      registry.add("h_tracketa", "", {HistType::kTH1F, {{100, trackEtaMin, trackEtaMax, "#it{#eta}"}}}, callSumw2);
      registry.add("h_trackphi", "", {HistType::kTH1F, {{100, 0.0, 2.0 * M_PI, "#it{#phi}"}}}, callSumw2);
    }

    if (doprocessMCDTracks) {
      registry.add("h2_trackpT_partpT", "", {HistType::kTH2F, {axisTrackpT, axisTrackpT}}, callSumw2);
      registry.add("h_partpT_matched_fine", "", {HistType::kTH1F, {axisTrackpTFine}}, callSumw2);
      registry.add("h_partpT", "", {HistType::kTH1F, {axisTrackpT}}, callSumw2);
      registry.add("h_partpT_fine", "", {HistType::kTH1F, {axisTrackpTFine}}, callSumw2);
    }

    if (doprocessDataJetsSel || doprocessMCDJetsSel) {
      registry.add("h_jetpT_sel8", "", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_selmc", "", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_tvx", "", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_coll", "", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_inelgt0rec", "", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_sel8_zvtx", "", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_selmc_zvtx", "", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_tvx_zvtx", "", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_coll_zvtx", "", {HistType::kTH1F, {axisJetpT}}, callSumw2);
    }

    if (doprocessMCDJets) {
      registry.add("h_jetpT_b", "b-jet", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_c", "c-jet", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_lf", "lf-jet", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_Db_b", "b-jet", {HistType::kTH1F, {axisDbFine}});
      registry.add("h_Db_c", "c-jet", {HistType::kTH1F, {axisDbFine}});
      registry.add("h_Db_lf", "lf-jet", {HistType::kTH1F, {axisDbFine}});
      registry.add("h2_jetpT_Db_b", "b-jet", {HistType::kTH2F, {axisJetpT, axisDb}});
      registry.add("h2_jetpT_Db_c", "c-jet", {HistType::kTH2F, {axisJetpT, axisDb}});
      registry.add("h2_jetpT_Db_lf", "lf-jet", {HistType::kTH2F, {axisJetpT, axisDb}});
      registry.add("h2_Response_DetjetpT_PartjetpT", "", {HistType::kTH2F, {axisJetpT, axisJetpT}}, callSumw2);
      registry.add("h2_Response_DetjetpT_PartjetpT_b", "b-jet", {HistType::kTH2F, {axisJetpT, axisJetpT}}, callSumw2);
      registry.add("h2_Response_DetjetpT_PartjetpT_c", "c-jet", {HistType::kTH2F, {axisJetpT, axisJetpT}}, callSumw2);
      registry.add("h2_Response_DetjetpT_PartjetpT_lf", "lf-jet", {HistType::kTH2F, {axisJetpT, axisJetpT}, callSumw2});
      registry.add("h2_jetpT_Db_lf_none", "lf-jet (none)", {HistType::kTH2F, {axisJetpT, axisDb}}, callSumw2);
      registry.add("h2_jetpT_Db_lf_matched", "lf-jet (matched)", {HistType::kTH2F, {axisJetpT, axisDb}}, callSumw2);
      registry.add("h2_jetpT_Db_npp", "NotPhysPrim", {HistType::kTH2F, {axisJetpT, axisDb}});
      registry.add("h2_jetpT_Db_npp_b", "NotPhysPrim b-jet", {HistType::kTH2F, {axisJetpT, axisDb}});
      registry.add("h2_jetpT_Db_npp_c", "NotPhysPrim c-jet", {HistType::kTH2F, {axisJetpT, axisDb}});
      registry.add("h2_jetpT_Db_npp_lf", "NotPhysPrim lf-jet", {HistType::kTH2F, {axisJetpT, axisDb}});
      registry.add("h_Db_npp", "NotPhysPrim", {HistType::kTH1F, {axisDbFine}});
      registry.add("h_Db_npp_b", "NotPhysPrim b-jet", {HistType::kTH1F, {axisDbFine}});
      registry.add("h_Db_npp_c", "NotPhysPrim c-jet", {HistType::kTH1F, {axisDbFine}});
      registry.add("h_Db_npp_lf", "NotPhysPrim lf-jet", {HistType::kTH1F, {axisDbFine}});
      registry.add("h_jetpT_matched", "", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_particle_matched", "", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_b_matched", "b-jet", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_particle_b_matched", "b-jet", {HistType::kTH1F, {axisJetpT}}, callSumw2);
    }

    if (doprocessMCDJetsSel) {
      registry.add("h_jetpT_b_sel8", "b-jet", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_b_sel8_zvtx", "b-jet", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h2_Response_DetjetpT_PartjetpT_sel8", "", {HistType::kTH2F, {axisJetpT, axisJetpT}}, callSumw2);
      registry.add("h2_Response_DetjetpT_PartjetpT_b_sel8", "b-jet", {HistType::kTH2F, {axisJetpT, axisJetpT}}, callSumw2);
      registry.add("h_jetpT_b_selmc", "b-jet", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_b_selmc_zvtx", "b-jet", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h2_Response_DetjetpT_PartjetpT_selmc", "", {HistType::kTH2F, {axisJetpT, axisJetpT}}, callSumw2);
      registry.add("h2_Response_DetjetpT_PartjetpT_b_selmc", "b-jet", {HistType::kTH2F, {axisJetpT, axisJetpT}}, callSumw2);
      registry.add("h_jetpT_b_tvx", "b-jet", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_b_tvx_zvtx", "b-jet", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h2_Response_DetjetpT_PartjetpT_tvx", "", {HistType::kTH2F, {axisJetpT, axisJetpT}}, callSumw2);
      registry.add("h2_Response_DetjetpT_PartjetpT_b_tvx", "b-jet", {HistType::kTH2F, {axisJetpT, axisJetpT}}, callSumw2);
      registry.add("h_jetpT_b_coll", "b-jet", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_b_coll_zvtx", "b-jet", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_b_inelgt0rec", "b-jet", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h2_Response_DetjetpT_PartjetpT_inelgt0", "", {HistType::kTH2F, {axisJetpT, axisJetpT}}, callSumw2);
      registry.add("h2_Response_DetjetpT_PartjetpT_b_inelgt0", "b-jet", {HistType::kTH2F, {axisJetpT, axisJetpT}}, callSumw2);
    }

    if (doprocessMCPJets) {
      registry.add("h_vertexZ_truth", "Vertex Z truth;#it{Z} (cm)", {HistType::kTH1F, {{100, -20.0, 20.0}}}, callSumw2);
      registry.add("h_vertexZ_truth_coll", "Vertex Z truth (Coll);#it{Z} (cm)", {HistType::kTH1F, {{100, -20.0, 20.0}}}, callSumw2);
      registry.add("h_vertexZ_truth_tvx", "Vertex Z truth (TVX);#it{Z} (cm)", {HistType::kTH1F, {{100, -20.0, 20.0}}}, callSumw2);

      registry.add("h_jetpT_particle", "", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_particle_b", "particle b-jet", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_particle_c", "particle c-jet", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_particle_lf", "particle lf-jet", {HistType::kTH1F, {axisJetpT}}, callSumw2);

      registry.add("h_jetpT_particle_sel8", "", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_particle_b_sel8", "particle b-jet", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_particle_selmc", "", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_particle_b_selmc", "particle b-jet", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_particle_tvx", "", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_particle_b_tvx", "particle b-jet", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_particle_coll", "", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_particle_b_coll", "particle b-jet", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_particle_inel", "", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_particle_b_inel", "particle b-jet", {HistType::kTH1F, {axisJetpT}}, callSumw2);

      registry.add("h_jetpT_particle_sel8_zvtx", "", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_particle_b_sel8_zvtx", "particle b-jet", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_particle_selmc_zvtx", "", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_particle_b_selmc_zvtx", "particle b-jet", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_particle_tvx_zvtx", "", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_particle_b_tvx_zvtx", "particle b-jet", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_particle_coll_zvtx", "", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_particle_b_coll_zvtx", "particle b-jet", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_particle_inel_zvtx", "", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_particle_b_inel_zvtx", "particle b-jet", {HistType::kTH1F, {axisJetpT}}, callSumw2);

      registry.add("h_jetpT_particle_inelgt0", "", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_particle_b_inelgt0", "particle b-jet", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_particle_inelgt0rec", "", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_particle_b_inelgt0rec", "particle b-jet", {HistType::kTH1F, {axisJetpT}}, callSumw2);
    }

    if (doDataDriven) {
      if (doDataDrivenExtra) {
        registry.add("hSparse_Incljets", "", {HistType::kTHnSparseF, {axisJetpT, axisDbFine, axisNTracks, axisJetMass}}, callSumw2);
        if (doprocessMCDJets) {
          registry.add("hSparse_bjets", "", {HistType::kTHnSparseF, {axisJetpT, axisDbFine, axisNTracks, axisJetMass}}, callSumw2);
          registry.add("hSparse_cjets", "", {HistType::kTHnSparseF, {axisJetpT, axisDbFine, axisNTracks, axisJetMass}}, callSumw2);
          registry.add("hSparse_lfjets", "", {HistType::kTHnSparseF, {axisJetpT, axisDbFine, axisNTracks, axisJetMass}}, callSumw2);
          registry.add("hSparse_lfjets_none", "", {HistType::kTHnSparseF, {axisJetpT, axisDbFine, axisNTracks, axisJetMass}}, callSumw2);
          registry.add("hSparse_lfjets_matched", "", {HistType::kTHnSparseF, {axisJetpT, axisDbFine, axisNTracks, axisJetMass}}, callSumw2);
        }
      } else {
        registry.add("hSparse_Incljets", "", {HistType::kTHnSparseF, {axisJetpT, axisDbFine, axisNTracks}}, callSumw2);
        if (doprocessMCDJets) {
          registry.add("hSparse_bjets", "", {HistType::kTHnSparseF, {axisJetpT, axisDbFine, axisNTracks}}, callSumw2);
          registry.add("hSparse_cjets", "", {HistType::kTHnSparseF, {axisJetpT, axisDbFine, axisNTracks}}, callSumw2);
          registry.add("hSparse_lfjets", "", {HistType::kTHnSparseF, {axisJetpT, axisDbFine, axisNTracks}}, callSumw2);
          registry.add("hSparse_lfjets_none", "", {HistType::kTHnSparseF, {axisJetpT, axisDbFine, axisNTracks}}, callSumw2);
          registry.add("hSparse_lfjets_matched", "", {HistType::kTHnSparseF, {axisJetpT, axisDbFine, axisNTracks}}, callSumw2);
        }
      }
    }
  }

  Filter collisionFilter = nabs(aod::jcollision::posZ) < vertexZCut;
  Filter trackFilter = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax);
  Filter jetFilter = (aod::jet::pt >= jetPtMin && aod::jet::pt < jetPtMax && aod::jet::eta < jetEtaMax - aod::jet::r / 100.f && aod::jet::eta > jetEtaMin + aod::jet::r / 100.f);

  using AnalysisCollisions = soa::Join<aod::JetCollisions, aod::JCollisionPIs>;
  using FilteredCollisions = soa::Filtered<AnalysisCollisions>;
  using AnalysisCollisionsTriggered = soa::Join<aod::JetCollisions, aod::JCollisionPIs, aod::JCollisionBCs>;
  using FilteredCollisionsTriggered = soa::Filtered<AnalysisCollisionsTriggered>;
  using OrigCollisions = soa::Join<aod::Collisions, aod::PVMults>;
  using DataJets = soa::Join<aod::ChargedJets, aod::ChargedJetConstituents, aod::ChargedJetTags>;
  using FilteredDataJets = soa::Filtered<DataJets>;
  using AnalysisTracks = soa::Join<aod::JetTracks, aod::JTrackExtras, aod::JTrackPIs>;
  using FilteredTracks = soa::Filtered<AnalysisTracks>;
  using OriginalTracks = soa::Join<aod::Tracks, aod::TracksCov, aod::TrackSelection, aod::TracksDCA, aod::TracksDCACov, aod::TracksExtra>;

  using MCDJets = soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets, aod::ChargedMCDetectorLevelJetFlavourDef, aod::ChargedMCDetectorLevelJetTags>;
  using FilteredMCDJets = soa::Filtered<MCDJets>;
  using AnalysisTracksMCD = soa::Join<aod::JetTracksMCD, aod::JTrackExtras, aod::JTrackPIs>;
  using FilteredTracksMCD = soa::Filtered<AnalysisTracksMCD>;

  using AnalysisCollisionsMCD = soa::Join<aod::JetCollisionsMCD, aod::JCollisionPIs>;
  using FilteredCollisionsMCD = soa::Filtered<AnalysisCollisionsMCD>;

  Filter mccollisionFilter = nabs(aod::jmccollision::posZ) < vertexZCut;
  using FilteredCollisionsMCP = soa::Filtered<aod::JetMcCollisions>;
  using MCPJets = soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets, aod::ChargedMCParticleLevelJetFlavourDef>;
  using FilteredMCPJets = soa::Filtered<MCPJets>;

  template <typename AnalysisJet, typename AnyTracks, typename AnyOriginalTracks>
  int analyzeJetTrackInfo(AnalysisJet const& analysisJet, AnyTracks const& /*allTracks*/, AnyOriginalTracks const /*origTracks*/ /*, int8_t jetFlavor = 0, double weight = 1.0*/)
  {
    int nTracks = 0;
    for (const auto& constituent : analysisJet.template tracks_as<AnyTracks>()) {

      if (constituent.pt() < trackPtMin || !jettaggingutilities::trackAcceptanceWithDca(constituent, maxIPxy, maxIPz)) {
        continue;
      }

      int sign = jettaggingutilities::getGeoSign(analysisJet, constituent);
      auto origConstit = constituent.template track_as<AnyOriginalTracks>();

      registry.fill(HIST("h_gnnfeat_trackpT"), constituent.pt());
      registry.fill(HIST("h_gnnfeat_trackPhi"), origConstit.phi());
      registry.fill(HIST("h_gnnfeat_trackEta"), constituent.eta());
      registry.fill(HIST("h_gnnfeat_trackCharge"), static_cast<float>(constituent.sign()));
      registry.fill(HIST("h_gnnfeat_trackDCAxy"), std::abs(constituent.dcaXY()) * sign);
      registry.fill(HIST("h_gnnfeat_trackSigmaDCAxy"), constituent.sigmadcaXY());
      registry.fill(HIST("h_gnnfeat_trackDCAz"), std::abs(constituent.dcaZ()) * sign);
      registry.fill(HIST("h_gnnfeat_trackSigmaDCAz"), constituent.sigmadcaZ());
      registry.fill(HIST("h_gnnfeat_trackITSNCls"), static_cast<float>(origConstit.itsNCls()));
      registry.fill(HIST("h_gnnfeat_trackTPCNCls"), static_cast<float>(origConstit.tpcNClsFound()));
      registry.fill(HIST("h_gnnfeat_trackTPCNCrossedRows"), static_cast<float>(origConstit.tpcNClsCrossedRows()));
      registry.fill(HIST("h_gnnfeat_trackITSChi2NCl"), origConstit.itsChi2NCl());
      registry.fill(HIST("h_gnnfeat_trackTPCChi2NCl"), origConstit.tpcChi2NCl());

      registry.fill(HIST("h_gnnfeat_tracksIPxy"), std::abs(constituent.dcaXY()) * sign / constituent.sigmadcaXY());
      registry.fill(HIST("h_gnnfeat_tracksIPz"), std::abs(constituent.dcaZ()) * sign / constituent.sigmadcaZ());

      ++nTracks;
    }
    return nTracks;
  }

  template <typename AnyTracks, typename AnalysisJet>
  bool isAcceptedJet(AnalysisJet const& jet)
  {

    if (jetAreaFractionMin > -98.0) {
      if (jet.area() < jetAreaFractionMin * M_PI * (jet.r() / 100.0) * (jet.r() / 100.0)) {
        return false;
      }
    }
    bool checkConstituentPt = true;
    bool checkConstituentMinPt = (leadingConstituentPtMin > -98.0);
    bool checkConstituentMaxPt = (leadingConstituentPtMax < 9998.0);
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

  template <typename AnalysisJet, typename AnyTracks, typename AnyOriginalTracks>
  void fillDataJetHistograms(AnalysisJet const& analysisJet, AnyTracks const& allTracks, AnyOriginalTracks const& origTracks)
  {
    int nTracks = analyzeJetTrackInfo(analysisJet, allTracks, origTracks);

    registry.fill(HIST("h_jetpT"), analysisJet.pt());
    registry.fill(HIST("h_jetPhi"), analysisJet.phi());
    registry.fill(HIST("h_jetEta"), analysisJet.eta());
    registry.fill(HIST("h_jetMass"), analysisJet.mass());
    registry.fill(HIST("h_Db"), analysisJet.scoreML());
    registry.fill(HIST("h2_jetpT_Db"), analysisJet.pt(), analysisJet.scoreML());

    if (doDataDriven) {
      if (doDataDrivenExtra) {
        registry.fill(HIST("hSparse_Incljets"), analysisJet.pt(), analysisJet.scoreML(), nTracks, analysisJet.mass());
      } else {
        registry.fill(HIST("hSparse_Incljets"), analysisJet.pt(), analysisJet.scoreML(), nTracks);
      }
    }
  }

  template <typename AnalysisJet, typename AnyTracks, typename AnyOriginalTracks>
  int8_t fillMCDJetHistograms(AnalysisJet const& analysisJet, AnyTracks const& /*allTracks*/, AnyOriginalTracks const& /*origTracks*/, double weightEvt = 1.0)
  {
    int8_t jetFlavor = analysisJet.origin();

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
      auto origConstit = constituent.template track_as<AnyOriginalTracks>();

      registry.fill(HIST("h_gnnfeat_trackpT"), constituent.pt());
      registry.fill(HIST("h_gnnfeat_trackPhi"), origConstit.phi());
      registry.fill(HIST("h_gnnfeat_trackEta"), constituent.eta());
      registry.fill(HIST("h_gnnfeat_trackCharge"), static_cast<float>(constituent.sign()));
      registry.fill(HIST("h_gnnfeat_trackDCAxy"), std::abs(constituent.dcaXY()) * sign);
      registry.fill(HIST("h_gnnfeat_trackSigmaDCAxy"), constituent.sigmadcaXY());
      registry.fill(HIST("h_gnnfeat_trackDCAz"), std::abs(constituent.dcaZ()) * sign);
      registry.fill(HIST("h_gnnfeat_trackSigmaDCAz"), constituent.sigmadcaZ());
      registry.fill(HIST("h_gnnfeat_trackITSNCls"), static_cast<float>(origConstit.itsNCls()));
      registry.fill(HIST("h_gnnfeat_trackTPCNCls"), static_cast<float>(origConstit.tpcNClsFound()));
      registry.fill(HIST("h_gnnfeat_trackTPCNCrossedRows"), static_cast<float>(origConstit.tpcNClsCrossedRows()));
      registry.fill(HIST("h_gnnfeat_trackITSChi2NCl"), origConstit.itsChi2NCl());
      registry.fill(HIST("h_gnnfeat_trackTPCChi2NCl"), origConstit.tpcChi2NCl());

      registry.fill(HIST("h_gnnfeat_tracksIPxy"), std::abs(constituent.dcaXY()) * sign / constituent.sigmadcaXY());
      registry.fill(HIST("h_gnnfeat_tracksIPz"), std::abs(constituent.dcaZ()) * sign / constituent.sigmadcaZ());

      ++nTracks;
    }

    registry.fill(HIST("h_jetpT"), analysisJet.pt(), weightEvt);
    registry.fill(HIST("h_jetPhi"), analysisJet.phi(), weightEvt);
    registry.fill(HIST("h_jetEta"), analysisJet.eta(), weightEvt);
    registry.fill(HIST("h_jetMass"), analysisJet.mass(), weightEvt);
    registry.fill(HIST("h_Db"), analysisJet.scoreML(), weightEvt);
    registry.fill(HIST("h2_jetpT_Db"), analysisJet.pt(), analysisJet.scoreML(), weightEvt);

    if (jetFlavor == JetTaggingSpecies::beauty) {
      registry.fill(HIST("h_jetpT_b"), analysisJet.pt(), weightEvt);
      registry.fill(HIST("h_Db_b"), analysisJet.scoreML(), weightEvt);
      registry.fill(HIST("h2_jetpT_Db_b"), analysisJet.pt(), analysisJet.scoreML(), weightEvt);
    } else if (jetFlavor == JetTaggingSpecies::charm) {
      registry.fill(HIST("h_jetpT_c"), analysisJet.pt(), weightEvt);
      registry.fill(HIST("h_Db_c"), analysisJet.scoreML(), weightEvt);
      registry.fill(HIST("h2_jetpT_Db_c"), analysisJet.pt(), analysisJet.scoreML(), weightEvt);
    } else {
      registry.fill(HIST("h_jetpT_lf"), analysisJet.pt(), weightEvt);
      registry.fill(HIST("h_Db_lf"), analysisJet.scoreML(), weightEvt);
      registry.fill(HIST("h2_jetpT_Db_lf"), analysisJet.pt(), analysisJet.scoreML(), weightEvt);
      if (jetFlavor == JetTaggingSpecies::none) {
        registry.fill(HIST("h2_jetpT_Db_lf_none"), analysisJet.pt(), analysisJet.scoreML(), weightEvt);
      } else {
        registry.fill(HIST("h2_jetpT_Db_lf_matched"), analysisJet.pt(), analysisJet.scoreML(), weightEvt);
      }
    }

    // Inspection for jets with predominant non-physical primary tracks
    if (static_cast<float>(nNppTracks) / nTracks > trackNppCrit) {
      registry.fill(HIST("h_Db_npp"), analysisJet.scoreML(), weightEvt);
      registry.fill(HIST("h2_jetpT_Db_npp"), analysisJet.pt(), analysisJet.scoreML(), weightEvt);
      if (jetFlavor == JetTaggingSpecies::beauty) {
        registry.fill(HIST("h_Db_npp_b"), analysisJet.scoreML(), weightEvt);
        registry.fill(HIST("h2_jetpT_Db_npp_b"), analysisJet.pt(), analysisJet.scoreML(), weightEvt);
      } else if (jetFlavor == JetTaggingSpecies::charm) {
        registry.fill(HIST("h_Db_npp_c"), analysisJet.scoreML(), weightEvt);
        registry.fill(HIST("h2_jetpT_Db_npp_c"), analysisJet.pt(), analysisJet.scoreML(), weightEvt);
      } else {
        registry.fill(HIST("h_Db_npp_lf"), analysisJet.scoreML(), weightEvt);
        registry.fill(HIST("h2_jetpT_Db_npp_lf"), analysisJet.pt(), analysisJet.scoreML(), weightEvt);
      }
    }

    if (doDataDriven) {
      if (doDataDrivenExtra) {
        registry.fill(HIST("hSparse_Incljets"), analysisJet.pt(), analysisJet.scoreML(), nTracks, analysisJet.mass(), weightEvt);
        if (jetFlavor == JetTaggingSpecies::beauty) {
          registry.fill(HIST("hSparse_bjets"), analysisJet.pt(), analysisJet.scoreML(), nTracks, analysisJet.mass(), weightEvt);
        } else if (jetFlavor == JetTaggingSpecies::charm) {
          registry.fill(HIST("hSparse_cjets"), analysisJet.pt(), analysisJet.scoreML(), nTracks, analysisJet.mass(), weightEvt);
        } else {
          registry.fill(HIST("hSparse_lfjets"), analysisJet.pt(), analysisJet.scoreML(), nTracks, analysisJet.mass(), weightEvt);
          if (jetFlavor == JetTaggingSpecies::none) {
            registry.fill(HIST("hSparse_lfjets_none"), analysisJet.pt(), analysisJet.scoreML(), nTracks, analysisJet.mass(), weightEvt);
          } else {
            registry.fill(HIST("hSparse_lfjets_matched"), analysisJet.pt(), analysisJet.scoreML(), nTracks, analysisJet.mass(), weightEvt);
          }
        }
      } else {
        registry.fill(HIST("hSparse_Incljets"), analysisJet.pt(), analysisJet.scoreML(), nTracks, weightEvt);
        if (jetFlavor == JetTaggingSpecies::beauty) {
          registry.fill(HIST("hSparse_bjets"), analysisJet.pt(), analysisJet.scoreML(), nTracks, weightEvt);
        } else if (jetFlavor == JetTaggingSpecies::charm) {
          registry.fill(HIST("hSparse_cjets"), analysisJet.pt(), analysisJet.scoreML(), nTracks, weightEvt);
        } else {
          registry.fill(HIST("hSparse_lfjets"), analysisJet.pt(), analysisJet.scoreML(), nTracks, weightEvt);
          if (jetFlavor == JetTaggingSpecies::none) {
            registry.fill(HIST("hSparse_lfjets_none"), analysisJet.pt(), analysisJet.scoreML(), nTracks, weightEvt);
          } else {
            registry.fill(HIST("hSparse_lfjets_matched"), analysisJet.pt(), analysisJet.scoreML(), nTracks, weightEvt);
          }
        }
      }
    }

    return jetFlavor;
  }

  // Check if the collision is INEL>0
  template <typename MCColl, typename MCPart>
  bool isTrueINEL0(MCColl const& /*mccoll*/, MCPart const& mcparts)
  {
    for (auto const& mcparticle : mcparts) {
      if (!mcparticle.isPhysicalPrimary())
        continue;
      auto p = pdg->GetParticle(mcparticle.pdgCode());
      if (p != nullptr) {
        if (std::abs(p->Charge()) >= 3) {
          if (std::abs(mcparticle.eta()) < 1)
            return true;
        }
      }
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
            if (collision.template collision_as<OrigCollisions>().isInelGt0()) {
              evtselCode |= EvtSelFlag::kINELgt0rec;
              registry.fill(HIST("hCollCounter"), static_cast<int>(EvtSel::INELgt0rec), weightEvt); // INELgt0+Zvtx(rec)
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

  void processDataJets(FilteredCollisions::iterator const& collision, FilteredDataJets const& alljets, FilteredTracks const& allTracks, OriginalTracks const& origTracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    registry.fill(HIST("h_event_counter"), 0.0); // Coll+TVX+NoTFB+NoITSROFB+Zvtx

    registry.fill(HIST("h_vertexZ"), collision.posZ());

    for (const auto& analysisJet : alljets) {
      if (!isAcceptedJet<FilteredTracks>(analysisJet)) {
        continue;
      }

      fillDataJetHistograms(analysisJet, allTracks, origTracks);
    }
  }
  PROCESS_SWITCH(BjetTaggingGnn, processDataJets, "jet information in Data", false);

  void processDataJetsTrig(FilteredCollisionsTriggered::iterator const& collision, FilteredDataJets const& alljets, FilteredTracks const& allTracks, OriginalTracks const& origTracks, aod::JBCs const& /*bcInfo*/)
  {
    // Get BC info associated with the collision before applying any event selections
    auto bc = collision.bc_as<aod::JBCs>();
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

    for (const auto& analysisJet : alljets) {
      if (!isAcceptedJet<FilteredTracks>(analysisJet)) {
        continue;
      }

      fillDataJetHistograms(analysisJet, allTracks, origTracks);
    }
  }
  PROCESS_SWITCH(BjetTaggingGnn, processDataJetsTrig, "jet information in software triggered Data", false);

  void processDataJetsSel(AnalysisCollisions::iterator const& collision, FilteredDataJets const& alljets, FilteredTracks const& /*allTracks*/, OrigCollisions const& /*origCollisions*/)
  {
    EvtSelFlag evtselCode = fillCollCounter(collision);

    for (const auto& analysisJet : alljets) {
      if (!isAcceptedJet<FilteredTracks>(analysisJet)) {
        continue;
      }

      registry.fill(HIST("h_jetpT_coll"), analysisJet.pt(), hasAll(evtselCode, EvtSelFlag::Coll));
      registry.fill(HIST("h_jetpT_coll_zvtx"), analysisJet.pt(), hasAll(evtselCode, EvtSelFlag::CollZvtx));
      registry.fill(HIST("h_jetpT_tvx"), analysisJet.pt(), hasAll(evtselCode, EvtSelFlag::TVX));
      registry.fill(HIST("h_jetpT_tvx_zvtx"), analysisJet.pt(), hasAll(evtselCode, EvtSelFlag::TVXZvtx));
      registry.fill(HIST("h_jetpT_selmc"), analysisJet.pt(), hasAll(evtselCode, EvtSelFlag::SelMC));
      registry.fill(HIST("h_jetpT_selmc_zvtx"), analysisJet.pt(), hasAll(evtselCode, EvtSelFlag::SelMCZvtx));
      registry.fill(HIST("h_jetpT_sel8"), analysisJet.pt(), hasAll(evtselCode, EvtSelFlag::Sel8));
      registry.fill(HIST("h_jetpT_sel8_zvtx"), analysisJet.pt(), hasAll(evtselCode, EvtSelFlag::Sel8Zvtx));
      registry.fill(HIST("h_jetpT_inelgt0rec"), analysisJet.pt(), hasAll(evtselCode, EvtSelFlag::INELgt0rec));
    }
  }
  PROCESS_SWITCH(BjetTaggingGnn, processDataJetsSel, "jet information in Data (event selection)", false);

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
    }
  }
  PROCESS_SWITCH(BjetTaggingGnn, processDataTracks, "track information in Data", false);

  void processMCDJets(FilteredCollisionsMCD::iterator const& collision, FilteredMCDJets const& MCDjets, FilteredTracksMCD const& allTracks, OriginalTracks const& origTracks, FilteredMCPJets const& /*MCPjets*/, aod::JetParticles const& /*mcParticles*/, FilteredCollisionsMCP const& /*mcCollisions*/)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    // Uses only collisionId % trainingDatasetRaioParam != 0 for evaluation dataset
    if (trainingDatasetRatioParam && collision.collisionId() % trainingDatasetRatioParam == 0) {
      return;
    }

    float weightEvt = useEventWeight ? collision.weight() : 1.f;
    bool matchedMcColl = collision.has_mcCollision() && std::fabs(collision.template mcCollision_as<FilteredCollisionsMCP>().posZ()) < vertexZCut;

    registry.fill(HIST("h_event_counter"), 0.0, weightEvt);

    registry.fill(HIST("h_vertexZ"), collision.posZ(), weightEvt);

    for (const auto& analysisJet : MCDjets) {
      if (!isAcceptedJet<FilteredTracksMCD>(analysisJet)) {
        continue;
      }

      float pTHat = 10. / (std::pow(weightEvt, 1.0 / pTHatExponent));
      if (useEventWeight && analysisJet.pt() > pTHatMaxMCD * pTHat) {
        continue;
      }

      int8_t jetFlavor = fillMCDJetHistograms(analysisJet, allTracks, origTracks, weightEvt);

      if (!matchedMcColl) {
        continue;
      }
      for (const auto& mcpjet : analysisJet.template matchedJetGeo_as<FilteredMCPJets>()) {
        if (mcpjet.pt() > pTHatMaxMCP * pTHat) {
          continue;
        }

        registry.fill(HIST("h2_Response_DetjetpT_PartjetpT"), analysisJet.pt(), mcpjet.pt(), weightEvt);
        registry.fill(HIST("h_jetpT_matched"), analysisJet.pt(), weightEvt);
        registry.fill(HIST("h_jetpT_particle_matched"), mcpjet.pt(), weightEvt);
        if (jetFlavor == JetTaggingSpecies::beauty) {
          registry.fill(HIST("h2_Response_DetjetpT_PartjetpT_b"), analysisJet.pt(), mcpjet.pt(), weightEvt);
          registry.fill(HIST("h_jetpT_b_matched"), analysisJet.pt(), weightEvt);
          registry.fill(HIST("h_jetpT_particle_b_matched"), mcpjet.pt(), weightEvt);
        } else if (jetFlavor == JetTaggingSpecies::charm) {
          registry.fill(HIST("h2_Response_DetjetpT_PartjetpT_c"), analysisJet.pt(), mcpjet.pt(), weightEvt);
        } else {
          registry.fill(HIST("h2_Response_DetjetpT_PartjetpT_lf"), analysisJet.pt(), mcpjet.pt(), weightEvt);
        }
      }
    }
  }
  PROCESS_SWITCH(BjetTaggingGnn, processMCDJets, "jet information in MC", false);

  void processMCDJetsSel(AnalysisCollisionsMCD::iterator const& collision, FilteredMCDJets const& MCDjets, FilteredTracksMCD const& /*allTracks*/, FilteredMCPJets const& /*MCPjets*/, OrigCollisions const& /*origCollisions*/, FilteredCollisionsMCP const& /*mcCollisions*/, aod::JetParticles const& mcParticles)
  {
    float weightEvt = useEventWeight ? collision.weight() : 1.f;
    EvtSelFlag evtselCode = fillCollCounter(collision, weightEvt);
    bool isTrueINELgt0 = collision.has_mcCollision() && isTrueINEL0(collision.template mcCollision_as<FilteredCollisionsMCP>(), mcParticles);

    for (const auto& analysisJet : MCDjets) {
      if (!isAcceptedJet<FilteredTracksMCD>(analysisJet)) {
        continue;
      }

      float pTHat = 10. / (std::pow(weightEvt, 1.0 / pTHatExponent));
      if (useEventWeight && analysisJet.pt() > pTHatMaxMCD * pTHat) {
        continue;
      }

      int8_t jetFlavor = analysisJet.origin();

      // Get matched particle-level jet pT
      double mcpjetpT = -1.0;
      for (const auto& mcpjet : analysisJet.template matchedJetGeo_as<FilteredMCPJets>()) {
        if (useEventWeight && mcpjet.pt() > pTHatMaxMCP * pTHat) {
          continue;
        }
        mcpjetpT = mcpjet.pt();
      }
      bool isMatched = mcpjetpT > 0.0;

      registry.fill(HIST("h_jetpT_coll"), analysisJet.pt(), hasAll(evtselCode, EvtSelFlag::Coll) ? weightEvt : 0.0);
      registry.fill(HIST("h_jetpT_coll_zvtx"), analysisJet.pt(), hasAll(evtselCode, EvtSelFlag::CollZvtx) ? weightEvt : 0.0);
      registry.fill(HIST("h_jetpT_tvx"), analysisJet.pt(), hasAll(evtselCode, EvtSelFlag::TVX) ? weightEvt : 0.0);
      registry.fill(HIST("h_jetpT_tvx_zvtx"), analysisJet.pt(), hasAll(evtselCode, EvtSelFlag::TVXZvtx) ? weightEvt : 0.0);
      registry.fill(HIST("h_jetpT_selmc"), analysisJet.pt(), hasAll(evtselCode, EvtSelFlag::SelMC) ? weightEvt : 0.0);
      registry.fill(HIST("h_jetpT_selmc_zvtx"), analysisJet.pt(), hasAll(evtselCode, EvtSelFlag::SelMCZvtx) ? weightEvt : 0.0);
      registry.fill(HIST("h_jetpT_sel8"), analysisJet.pt(), hasAll(evtselCode, EvtSelFlag::Sel8) ? weightEvt : 0.0);
      registry.fill(HIST("h_jetpT_sel8_zvtx"), analysisJet.pt(), hasAll(evtselCode, EvtSelFlag::Sel8Zvtx) ? weightEvt : 0.0);
      registry.fill(HIST("h_jetpT_inelgt0rec"), analysisJet.pt(), hasAll(evtselCode, EvtSelFlag::INELgt0rec) ? weightEvt : 0.0);
      if (isMatched) {
        registry.fill(HIST("h2_Response_DetjetpT_PartjetpT_selmc"), analysisJet.pt(), mcpjetpT, hasAll(evtselCode, EvtSelFlag::SelMCZvtx) ? weightEvt : 0.0);
        registry.fill(HIST("h2_Response_DetjetpT_PartjetpT_sel8"), analysisJet.pt(), mcpjetpT, hasAll(evtselCode, EvtSelFlag::Sel8Zvtx) ? weightEvt : 0.0);
        registry.fill(HIST("h2_Response_DetjetpT_PartjetpT_inelgt0"), analysisJet.pt(), mcpjetpT, isTrueINELgt0 && (hasAll(evtselCode, EvtSelFlag::INELgt0rec)) ? weightEvt : 0.0);
      }
      if (jetFlavor == JetTaggingSpecies::beauty) {
        registry.fill(HIST("h_jetpT_b_coll"), analysisJet.pt(), hasAll(evtselCode, EvtSelFlag::Coll) ? weightEvt : 0.0);
        registry.fill(HIST("h_jetpT_b_coll_zvtx"), analysisJet.pt(), hasAll(evtselCode, EvtSelFlag::CollZvtx) ? weightEvt : 0.0);
        registry.fill(HIST("h_jetpT_b_tvx"), analysisJet.pt(), hasAll(evtselCode, EvtSelFlag::TVX) ? weightEvt : 0.0);
        registry.fill(HIST("h_jetpT_b_tvx_zvtx"), analysisJet.pt(), hasAll(evtselCode, EvtSelFlag::TVXZvtx) ? weightEvt : 0.0);
        registry.fill(HIST("h_jetpT_b_selmc"), analysisJet.pt(), hasAll(evtselCode, EvtSelFlag::SelMC) ? weightEvt : 0.0);
        registry.fill(HIST("h_jetpT_b_selmc_zvtx"), analysisJet.pt(), hasAll(evtselCode, EvtSelFlag::SelMCZvtx) ? weightEvt : 0.0);
        registry.fill(HIST("h_jetpT_b_sel8"), analysisJet.pt(), hasAll(evtselCode, EvtSelFlag::Sel8) ? weightEvt : 0.0);
        registry.fill(HIST("h_jetpT_b_sel8_zvtx"), analysisJet.pt(), hasAll(evtselCode, EvtSelFlag::Sel8Zvtx) ? weightEvt : 0.0);
        registry.fill(HIST("h_jetpT_b_inelgt0rec"), analysisJet.pt(), hasAll(evtselCode, EvtSelFlag::INELgt0rec) ? weightEvt : 0.0);
        if (isMatched) {
          registry.fill(HIST("h2_Response_DetjetpT_PartjetpT_b_selmc"), analysisJet.pt(), mcpjetpT, hasAll(evtselCode, EvtSelFlag::SelMCZvtx) ? weightEvt : 0.0);
          registry.fill(HIST("h2_Response_DetjetpT_PartjetpT_b_sel8"), analysisJet.pt(), mcpjetpT, hasAll(evtselCode, EvtSelFlag::Sel8Zvtx) ? weightEvt : 0.0);
          registry.fill(HIST("h2_Response_DetjetpT_PartjetpT_b_inelgt0"), analysisJet.pt(), mcpjetpT, isTrueINELgt0 && (hasAll(evtselCode, EvtSelFlag::INELgt0rec)) ? weightEvt : 0.0);
        }
      }
    }
  }
  PROCESS_SWITCH(BjetTaggingGnn, processMCDJetsSel, "jet information in MC (event selection)", false);

  PresliceUnsorted<AnalysisCollisionsMCD> collisionsPerMCPCollision = aod::jmccollisionlb::mcCollisionId;
  Preslice<FilteredMCPJets> mcpjetsPerMCPCollision = aod::jmccollisionlb::mcCollisionId;

  void processMCPJets(aod::McCollisions::iterator const& mcCollision, FilteredMCPJets const& mcpjets, AnalysisCollisionsMCD const& collisions, OrigCollisions const& /*origCollisions*/, aod::JetParticles const& mcParticles)
  {
    float weightEvt = useEventWeight ? mcCollision.weight() : 1.f;
    auto matchedCollisions = collisions.sliceBy(collisionsPerMCPCollision, mcCollision.globalIndex());

    EvtSelFlag evtselCode = EvtSelFlag::INEL;
    registry.fill(HIST("hMcCollCounter"), static_cast<int>(EvtSel::INEL), weightEvt); // INEL
    registry.fill(HIST("h_vertexZ_truth"), mcCollision.posZ(), weightEvt);

    bool zvtx = std::fabs(mcCollision.posZ()) < vertexZCut;
    bool zvtxMatched = false;

    bool isTrueINELgt0 = isTrueINEL0(mcCollision, mcParticles);
    if (zvtx) {
      evtselCode |= EvtSelFlag::kZvtx;
      registry.fill(HIST("hMcCollCounter"), static_cast<int>(EvtSel::INELZvtx), weightEvt); // INEL+Zvtx
      if (isTrueINELgt0) {
        evtselCode |= EvtSelFlag::kINELgt0;
        registry.fill(HIST("hMcCollCounter"), static_cast<int>(EvtSel::INELgt0), weightEvt); // INELgt0
      }
    }

    bool isMatchedToAnalysisSelection = false;

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
              if (matchedCollisions.begin().template collision_as<OrigCollisions>().isInelGt0()) {
                evtselCode |= EvtSelFlag::kINELgt0rec;
                registry.fill(HIST("hMcCollCounter"), static_cast<int>(EvtSel::INELgt0rec), weightEvt); // INELgt0+Zvtx(rec)
              }
            }
          }
        }
      }
      if (jetderiveddatautilities::selectCollision(matchedCollisions.begin(), eventSelectionBits) && zvtxMatched) {
        isMatchedToAnalysisSelection = true;
        registry.fill(HIST("h_event_counter_mcp"), 0.0, weightEvt);
      }
    }

    auto mcpjetspermcpcollision = mcpjets.sliceBy(mcpjetsPerMCPCollision, mcCollision.globalIndex());
    for (const auto& mcpjet : mcpjetspermcpcollision) {
      bool jetIncluded = false;
      for (const auto& jetR : jetRadiiValues) {
        if (mcpjet.r() == static_cast<int>(jetR * 100)) {
          jetIncluded = true;
          break;
        }
      }

      if (!jetIncluded) {
        continue;
      }

      float pTHat = 10. / (std::pow(weightEvt, 1.0 / pTHatExponent));
      if (useEventWeight && mcpjet.pt() > pTHatMaxMCP * pTHat) {
        continue;
      }

      int8_t jetFlavor = mcpjet.origin();

      registry.fill(HIST("h_jetpT_particle_inel"), mcpjet.pt(), hasAll(evtselCode, EvtSelFlag::INEL) ? weightEvt : 0.0);
      registry.fill(HIST("h_jetpT_particle_inel_zvtx"), mcpjet.pt(), hasAll(evtselCode, EvtSelFlag::INELZvtx) ? weightEvt : 0.0);
      registry.fill(HIST("h_jetpT_particle_coll"), mcpjet.pt(), hasAll(evtselCode, EvtSelFlag::Coll) ? weightEvt : 0.0);
      registry.fill(HIST("h_jetpT_particle_coll_zvtx"), mcpjet.pt(), hasAll(evtselCode, EvtSelFlag::CollZvtx) ? weightEvt : 0.0);
      registry.fill(HIST("h_jetpT_particle_tvx"), mcpjet.pt(), hasAll(evtselCode, EvtSelFlag::TVX) ? weightEvt : 0.0);
      registry.fill(HIST("h_jetpT_particle_tvx_zvtx"), mcpjet.pt(), hasAll(evtselCode, EvtSelFlag::TVXZvtx) ? weightEvt : 0.0);
      registry.fill(HIST("h_jetpT_particle_selmc"), mcpjet.pt(), hasAll(evtselCode, EvtSelFlag::SelMC) ? weightEvt : 0.0);
      registry.fill(HIST("h_jetpT_particle_selmc_zvtx"), mcpjet.pt(), hasAll(evtselCode, EvtSelFlag::SelMCZvtx) ? weightEvt : 0.0);
      registry.fill(HIST("h_jetpT_particle_sel8"), mcpjet.pt(), hasAll(evtselCode, EvtSelFlag::Sel8) ? weightEvt : 0.0);
      registry.fill(HIST("h_jetpT_particle_sel8_zvtx"), mcpjet.pt(), hasAll(evtselCode, EvtSelFlag::Sel8Zvtx) ? weightEvt : 0.0);
      registry.fill(HIST("h_jetpT_particle_inelgt0"), mcpjet.pt(), hasAll(evtselCode, EvtSelFlag::INELgt0) ? weightEvt : 0.0);
      registry.fill(HIST("h_jetpT_particle_inelgt0rec"), mcpjet.pt(), hasAll(evtselCode, EvtSelFlag::INELgt0rec) ? weightEvt : 0.0);
      if (jetFlavor == JetTaggingSpecies::beauty) {
        registry.fill(HIST("h_jetpT_particle_b_inel"), mcpjet.pt(), hasAll(evtselCode, EvtSelFlag::INEL) ? weightEvt : 0.0);
        registry.fill(HIST("h_jetpT_particle_b_inel_zvtx"), mcpjet.pt(), hasAll(evtselCode, EvtSelFlag::INELZvtx) ? weightEvt : 0.0);
        registry.fill(HIST("h_jetpT_particle_b_coll"), mcpjet.pt(), hasAll(evtselCode, EvtSelFlag::Coll) ? weightEvt : 0.0);
        registry.fill(HIST("h_jetpT_particle_b_coll_zvtx"), mcpjet.pt(), hasAll(evtselCode, EvtSelFlag::CollZvtx) ? weightEvt : 0.0);
        registry.fill(HIST("h_jetpT_particle_b_tvx"), mcpjet.pt(), hasAll(evtselCode, EvtSelFlag::TVX) ? weightEvt : 0.0);
        registry.fill(HIST("h_jetpT_particle_b_tvx_zvtx"), mcpjet.pt(), hasAll(evtselCode, EvtSelFlag::TVXZvtx) ? weightEvt : 0.0);
        registry.fill(HIST("h_jetpT_particle_b_selmc"), mcpjet.pt(), hasAll(evtselCode, EvtSelFlag::SelMC) ? weightEvt : 0.0);
        registry.fill(HIST("h_jetpT_particle_b_selmc_zvtx"), mcpjet.pt(), hasAll(evtselCode, EvtSelFlag::SelMCZvtx) ? weightEvt : 0.0);
        registry.fill(HIST("h_jetpT_particle_b_sel8"), mcpjet.pt(), hasAll(evtselCode, EvtSelFlag::Sel8) ? weightEvt : 0.0);
        registry.fill(HIST("h_jetpT_particle_b_sel8_zvtx"), mcpjet.pt(), hasAll(evtselCode, EvtSelFlag::Sel8Zvtx) ? weightEvt : 0.0);
        registry.fill(HIST("h_jetpT_particle_b_inelgt0"), mcpjet.pt(), hasAll(evtselCode, EvtSelFlag::INELgt0) ? weightEvt : 0.0);
        registry.fill(HIST("h_jetpT_particle_b_inelgt0rec"), mcpjet.pt(), hasAll(evtselCode, EvtSelFlag::INELgt0rec) ? weightEvt : 0.0);
      }

      // Fill histograms for jets matched to the analysis event selection
      if (isMatchedToAnalysisSelection) {
        registry.fill(HIST("h_jetpT_particle"), mcpjet.pt(), weightEvt);

        if (jetFlavor == JetTaggingSpecies::beauty) {
          registry.fill(HIST("h_jetpT_particle_b"), mcpjet.pt(), weightEvt);
        } else if (jetFlavor == JetTaggingSpecies::charm) {
          registry.fill(HIST("h_jetpT_particle_c"), mcpjet.pt(), weightEvt);
        } else {
          registry.fill(HIST("h_jetpT_particle_lf"), mcpjet.pt(), weightEvt);
        }
      }
    }
  }
  PROCESS_SWITCH(BjetTaggingGnn, processMCPJets, "mc collision information", false);

  Preslice<aod::JetParticles> mcparticlesPerMCPCollision = aod::jmcparticle::mcCollisionId;

  void processMCDTracks(FilteredCollisionsMCD::iterator const& collision, AnalysisTracksMCD const& tracks, FilteredCollisionsMCP const& /*mcCollisions*/, aod::JetParticles const& allParticles)
  {
    float weightEvt = useEventWeight ? collision.weight() : 1.f;
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    // Uses only collisionId % trainingDatasetRaioParam != 0 for evaluation dataset
    if (trainingDatasetRatioParam && collision.collisionId() % trainingDatasetRatioParam == 0) {
      return;
    }

    bool matchedMcColl = collision.has_mcCollision() && std::fabs(collision.template mcCollision_as<FilteredCollisionsMCP>().posZ()) < vertexZCut;

    for (const auto& track : tracks) {
      if (track.eta() <= trackEtaMin || track.eta() >= trackEtaMax) {
        continue;
      }
      registry.fill(HIST("h_trackpT"), track.pt(), weightEvt);
      registry.fill(HIST("h_tracketa"), track.eta(), weightEvt);
      registry.fill(HIST("h_trackphi"), track.phi(), weightEvt);

      if (!matchedMcColl || !track.has_mcParticle()) {
        continue;
      }
      auto particle = track.template mcParticle_as<aod::JetParticles>();
      if (particle.isPhysicalPrimary() && particle.eta() > trackEtaMin && particle.eta() < trackEtaMax) {
        registry.fill(HIST("h2_trackpT_partpT"), track.pt(), particle.pt(), weightEvt);
        registry.fill(HIST("h_partpT_matched_fine"), particle.pt(), weightEvt);
      }
    }

    if (!matchedMcColl) {
      return;
    }

    auto const particles = allParticles.sliceBy(mcparticlesPerMCPCollision, collision.mcCollisionId());

    for (const auto& particle : particles) {
      auto pdgParticle = pdg->GetParticle(particle.pdgCode());
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
      auto collisionsInBC = collisions.sliceBy(perFoundBC, bc.globalIndex());
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
