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
///
/// \file longrangecorrDerived.cxx
///
/// \brief task for long range correlation analysis based on derived table
/// \author Abhi Modak (abhi.modak@cern.ch)
/// \since November 05, 2025

#include "PWGCF/Core/CorrelationContainer.h"
#include "PWGCF/TwoParticleCorrelations/DataModel/LongRangeDerived.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGUD/Core/SGSelector.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/McCollisionExtra.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH3.h>
#include <TPDGCode.h>
#include <TRandom.h>

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::track;
using namespace o2::aod::fwdtrack;
using namespace o2::aod::evsel;
using namespace o2::constants::math;

auto static constexpr KminCharge = 3.0f;
auto static constexpr KPidMaskPion = 2;
auto static constexpr KPidMaskKaon = 4;
auto static constexpr KPidMaskProton = 8;

struct LongrangecorrDerived {

  SliceCache cache;
  SGSelector sgSelector;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  TrackSelection myTrackFilter;
  Service<o2::framework::O2DatabasePDG> pdg{};
  Service<o2::ccdb::BasicCCDBManager> ccdb{};

  struct : ConfigurableGroup {
    Configurable<int> cfgNmixedevent{"cfgNmixedevent", 5, "how many events are mixed"};
    Configurable<double> cfgSampleSize{"cfgSampleSize", 10, "Sample size for bootstrapping"};
    Configurable<int> cfgPidMask{"cfgPidMask", 0, "Selection bitmask for the TPC particle"};
    Configurable<int> cfgV0Mask{"cfgV0Mask", 0, "Selection bitmask for the V0 particle"};
    Configurable<float> cfgVtxCut{"cfgVtxCut", 10.0f, "Vertex Z range to consider"};
    Configurable<bool> isUseCentEst{"isUseCentEst", false, "Centrality based classification"};
    Configurable<int> isUseDataLikeMult{"isUseDataLikeMult", 0, "Data like mult/cent classification"};
    Configurable<int> cfgVerbosity{"cfgVerbosity", 0, "print statement"};

    Configurable<float> cfgEtaCut{"cfgEtaCut", 0.8f, "Eta range to consider"};
    Configurable<float> cfgPtCutMin{"cfgPtCutMin", 0.2f, "minimum accepted track pT"};
    Configurable<float> cfgPtCutMax{"cfgPtCutMax", 10.0f, "maximum accepted track pT"};
    Configurable<float> cfgPtCutMinForMult{"cfgPtCutMinForMult", 0.2f, "minimum track pT for fallback multiplicity classification"};
    Configurable<float> cfgPtCutMaxForMult{"cfgPtCutMaxForMult", 3.0f, "maximum track pT for fallback multiplicity classification"};
    Configurable<float> cfgTpcMinNclsFound{"cfgTpcMinNclsFound", 50.0f, ""};
    Configurable<float> cfgTpcMinNCrossedRows{"cfgTpcMinNCrossedRows", 70.0f, ""};
    Configurable<float> cfgTpcMaxChi2PerCluster{"cfgTpcMaxChi2PerCluster", 4.0f, ""};
    Configurable<float> cfgTpcMaxDcaZ{"cfgTpcMaxDcaZ", 1.0f, ""};
    Configurable<bool> applyEffCorr{"applyEffCorr", true, "Enable efficiency correction"};
    Configurable<bool> applyAccCorr{"applyAccCorr", false, "Enable NUA correction"};
    Configurable<std::string> cfgEffccdbPath{"cfgEffccdbPath", "Users/a/abmodak/Efficiency/OO/default", "Browse track eff object from CCDB"};
    Configurable<std::string> cfgAccccdbPath{"cfgAccccdbPath", "Users/a/abmodak/Acceptance/OO/default", "Browse track eff object from CCDB"};

    Configurable<bool> cfgRequireFt0aOuterRing{"cfgRequireFt0aOuterRing", false, "Consider FT0A Outer Ring"};
    Configurable<bool> cfgRequireFt0aInnerRing{"cfgRequireFt0aInnerRing", false, "Consider FT0A Inner Ring"};
    Configurable<bool> cfgRequireFt0cOuterRing{"cfgRequireFt0cOuterRing", false, "Consider FT0C Outer Ring"};
    Configurable<bool> cfgRequireFt0cInnerRing{"cfgRequireFt0cInnerRing", false, "Consider FT0C Inner Ring"};

    Configurable<float> cfgTofPidPtCut{"cfgTofPidPtCut", 0.4f, "Minimum pt to use TOF N-sigma"};
    Configurable<float> cfgPidNsigmaMax{"cfgPidNsigmaMax", 1.5f, "Maximum n-sigma for PID"};
    Configurable<float> cfgPidNsigmaMin{"cfgPidNsigmaMin", -1.5f, "Minimum n-sigma for PID"};
    Configurable<float> cfgCleanNsigmaTpcMax{"cfgCleanNsigmaTpcMax", 10.0f, "Max TPC n-sigma for clean templates (wide to avoid dE/dx truncation)"};
    Configurable<float> cfgCleanNsigmaTpcMin{"cfgCleanNsigmaTpcMin", -10.0f, "Min TPC n-sigma for clean templates (wide to avoid dE/dx truncation)"};
    Configurable<float> cfgCleanNsigmaTofMax{"cfgCleanNsigmaTofMax", 1.5f, "Max TOF n-sigma for clean templates"};
    Configurable<float> cfgCleanNsigmaTofMin{"cfgCleanNsigmaTofMin", -1.5f, "Min TOF n-sigma for clean templates"};
    Configurable<float> cfgCleanTofVetoNsigma{"cfgCleanTofVetoNsigma", 3.0f, "Veto other species within this |nSigmaTOF| for track-level TOF-clean templates"};
    Configurable<bool> cfgGetNsigmaQA{"cfgGetNsigmaQA", true, "Get QA histograms for PID selection"};
    Configurable<bool> cfgGetdEdx{"cfgGetdEdx", true, "Get dEdx histograms for TPC signal"};
    Configurable<bool> cfgGetAllSpeciesQA{"cfgGetAllSpeciesQA", true, "Fill unified all-species dE/dx and expected sigma histograms"};
    Configurable<bool> cfgGetV0CleanSample{"cfgGetV0CleanSample", true, "Fill clean V0 daughter (K0s -> pi pi, Lambda -> p pi) dE/dx templates"};
    Configurable<bool> cfgFillAfterCutWithV0{"cfgFillAfterCutWithV0", true, "Fill TpcdEdx_ptwise_afterCut and TofTpcNsigma_after for pi and p using V0 clean daughters"};
    Configurable<bool> cfgRequireTofForV0AfterCut{"cfgRequireTofForV0AfterCut", true, "Require hasTOF and |nSigmaTOF| < cfgCleanNsigmaTofMax when filling afterCut from V0s"};
    Configurable<bool> cfgIncludeLambdaPiInCleanPi{"cfgIncludeLambdaPiInCleanPi", true, "Also include clean pions from Lambda -> p pi in TpcdEdx_ptwise_afterCut for pion wagon"};
    Configurable<bool> cfgAllowFallbackMultInPIDQA{"cfgAllowFallbackMultInPIDQA", false, "If true and LRCollisions is empty, use fallback multiplicity in processPIDQA"};
    Configurable<bool> cfgRequireTofInMCdEdx{"cfgRequireTofInMCdEdx", true, "Require hasTOF and |nSigmaTOF| < 1.5 in MC_TpcdEdx_ptwise"};

    // V0 topological cuts for clean daughter template extraction (tag-and-probe without TPC cut on probe daughter)
    Configurable<float> cfgCleanV0MinCosPA{"cfgCleanV0MinCosPA", 0.995f, "Min cos(PA) for clean V0 daughters"};
    Configurable<float> cfgCleanV0MaxDcaDau{"cfgCleanV0MaxDcaDau", 0.8f, "Max DCA between V0 daughters (cm)"};
    Configurable<float> cfgCleanV0MinRadius{"cfgCleanV0MinRadius", 1.2f, "Min V0 2D decay radius (cm)"};
    Configurable<float> cfgCleanV0MaxRadius{"cfgCleanV0MaxRadius", 30.0f, "Max V0 2D decay radius (cm)"};
    Configurable<float> cfgCleanK0sMassWin{"cfgCleanK0sMassWin", 0.008f, "Tight K0s mass window around PDG mass (GeV/c^2)"};
    Configurable<float> cfgCleanK0sMinArmQt{"cfgCleanK0sMinArmQt", 0.2f, "Min Armenteros qT/|alpha| for clean K0s"};
    Configurable<float> cfgCleanK0sMaxCtau{"cfgCleanK0sMaxCtau", 20.0f, "Max proper lifetime c*tau for clean K0s (cm)"};
    Configurable<float> cfgCleanK0sMinDcaToPv{"cfgCleanK0sMinDcaToPv", 0.1f, "Min pion DCA to PV for clean K0s (cm)"};
    Configurable<float> cfgCleanLambdaMassWin{"cfgCleanLambdaMassWin", 0.005f, "Tight Lambda mass window around PDG mass (GeV/c^2)"};
    Configurable<float> cfgCleanLambdaRejK0sWin{"cfgCleanLambdaRejK0sWin", 0.012f, "Competing K0s mass rejection window for clean Lambda (GeV/c^2)"};
    Configurable<float> cfgCleanLambdaMinArmQt{"cfgCleanLambdaMinArmQt", 0.02f, "Min Armenteros qT for clean Lambda to reject photon conversions (GeV/c)"};
    Configurable<float> cfgCleanLambdaMaxCtau{"cfgCleanLambdaMaxCtau", 30.0f, "Max proper lifetime c*tau for clean Lambda (cm)"};
    Configurable<float> cfgCleanLambdaMinDcaPr{"cfgCleanLambdaMinDcaPr", 0.07f, "Min proton DCA to PV for clean Lambda (cm)"};
    Configurable<float> cfgCleanLambdaMinDcaPi{"cfgCleanLambdaMinDcaPi", 0.2f, "Min pion DCA to PV for clean Lambda (cm)"};
    Configurable<float> cfgCleanV0TagTpcNsigma{"cfgCleanV0TagTpcNsigma", 3.0f, "Max |nSigmaTPC| on companion tag daughter"};
  } cfgSel;

  struct : ConfigurableGroup {
    ConfigurableAxis axisMultQA{"axisMultQA", {500, -0.5, 499.5}, "multiplicity QA axis"};
    ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 26, 46, 133, 1000}, "multiplicity axis"};
    ConfigurableAxis axisCentrality{"axisCentrality", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100}, "Centrality (%)"};

    ConfigurableAxis axisPhi{"axisPhi", {96, 0, TwoPI}, "#phi axis"};
    ConfigurableAxis axisEtaTrig{"axisEtaTrig", {40, -1., 1.}, "#eta trig axis"};
    ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt trigger axis for histograms"};
    ConfigurableAxis axisPtQA{"axisPtQA", {80, 0.0, 8.0}, "pt axis for PID QA histograms"};

    ConfigurableAxis axisVtxZ{"axisVtxZ", {40, -20, 20}, "vertex axis"};
    ConfigurableAxis axisEtaAssoc{"axisEtaAssoc", {96, 3.5, 4.9}, "#eta assoc axis"};
    ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {72, -PIHalf, PIHalf * 3}, "delta phi axis for histograms"};
    ConfigurableAxis axisDeltaEta{"axisDeltaEta", {40, -6, -2}, "delta eta axis for histograms"};
    ConfigurableAxis axisInvMass{"axisInvMass", {VARIABLE_WIDTH, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2.0}, "invariant mass axis"};
    ConfigurableAxis axisInvMassQA{"axisInvMassQA", {20, 0.45, 0.55}, "invariant mass axis for QA"};
    ConfigurableAxis axisAmplitude{"axisAmplitude", {5000, 0, 10000}, "FT0 amplitude"};
    ConfigurableAxis axisChannel{"axisChannel", {208, 0, 208}, "FT0 channel"};
    ConfigurableAxis axisMultME{"axisMultME", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 1000}, "Mixing bins - multiplicity"};
    ConfigurableAxis axisVtxZME{"axisVtxZME", {VARIABLE_WIDTH, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10}, "Mixing bins - z-vertex"};
    ConfigurableAxis axisSample{"axisSample", {10, 0, 10}, "sample axis for histograms"};

    ConfigurableAxis axisTPCNClsFound{"axisTPCNClsFound", {200, -0.5, 199.5}, "TPC Cluster axis"};
    ConfigurableAxis axisTPCNClsCrossedRows{"axisTPCNClsCrossedRows", {200, -0.5, 199.5}, "TPC NCrossedRow axis"};
    ConfigurableAxis axisTPCChi2NCl{"axisTPCChi2NCl", {20, 0.0, 20.0}, "TPC Chi2/NCl axis"};
    ConfigurableAxis axisTPCdcaZ{"axisTPCdcaZ", {200, -10.0, 10.0}, "TPC dcaZ axis"};

    ConfigurableAxis axisVertexEfficiency{"axisVertexEfficiency", {10, -10, 10}, "vertex axis for efficiency histograms"};
    ConfigurableAxis axisEtaEfficiency{"axisEtaEfficiency", {20, -1.0, 1.0}, "eta axis for efficiency histograms"};
    ConfigurableAxis axisPtEfficiency{"axisPtEfficiency", {1, 0.5, 8.0}, "pt axis for efficiency histograms"};

    ConfigurableAxis axisNsigmaTPC{"axisNsigmaTPC", {80, -5, 5}, "nsigmaTPC axis"};
    ConfigurableAxis axisNsigmaTOF{"axisNsigmaTOF", {80, -5, 5}, "nsigmaTOF axis"};
    ConfigurableAxis axisTpcSignal{"axisTpcSignal", {250, 0, 250}, "dEdx axis for TPC"};
    ConfigurableAxis axisTpcExpSigma{"axisTpcExpSigma", {200, 0, 20}, "Expected sigma axis for TPC dE/dx"};

    ConfigurableAxis axisMcPt{"axisMcPt", {100, 0.2, 10.0}, "p_{T} (GeV/c)"};
    ConfigurableAxis axisMcEta{"axisMcEta", {100, -0.8, 0.8}, "#eta"};
    ConfigurableAxis axisMcPhi{"axisMcPhi", {100, 0., TwoPI}, "#phi"};
  } cfgAxis;

  Configurable<float> cfgFv0Cut{"cfgFv0Cut", 50.0f, "FV0A threshold"};
  Configurable<float> cfgFt0aCut{"cfgFt0aCut", 100.0f, "FT0A threshold"};
  Configurable<float> cfgFt0cCut{"cfgFt0cCut", 50.0f, "FT0C threshold"};
  Configurable<float> cfgZdcCut{"cfgZdcCut", 0.1f, "ZDC threshold"};
  Configurable<int> cfgGapSideCut{"cfgGapSideCut", 0, "Gap-side A=0, C=1, AC = 2, No Gap = -1, All events = 3"};

  OutputObj<CorrelationContainer> same{"sameEvent"};
  OutputObj<CorrelationContainer> mixed{"mixedEvent"};

  // corrections
  TH3D* hTrkEff = nullptr;
  TH3D* hTrkAcc = nullptr;
  bool fLoadTrkEffCorr = false;

  using CollsTable = aod::LRCollisions;
  using TrksTable = aod::LRMidTracks;
  using Ft0aTrksTable = aod::LRFt0aTracks;
  using Ft0cTrksTable = aod::LRFt0cTracks;
  using V0TrksTable = aod::LRV0Tracks;

  using McCollsTable = aod::LRMcCollisions;
  using McTrksTable = aod::LRMidMcTracks;
  using McFt0aTrksTable = aod::LRFt0aMcTracks;
  using McFt0cTrksTable = aod::LRFt0cMcTracks;

  using UpcCollsTable = soa::Join<aod::UpcLRCollisions, aod::UpcSgLRCollisions, aod::LRZdcs>;
  using TrksUpcTable = aod::UpcLRMidTracks;
  using Ft0aTrksUpcTable = aod::UpcLRFt0aTracks;
  using Ft0cTrksUpcTable = aod::UpcLRFt0cTracks;
  using V0TrksUpcTable = aod::UpcLRV0Tracks;

  Preslice<TrksTable> perColTpc = aod::lrcorrtrktable::lrCollisionId;
  Preslice<Ft0aTrksTable> perColFt0a = aod::lrcorrtrktable::lrCollisionId;
  Preslice<Ft0cTrksTable> perColFt0c = aod::lrcorrtrktable::lrCollisionId;
  Preslice<V0TrksTable> perColV0 = aod::lrcorrtrktable::lrCollisionId;

  Preslice<TrksUpcTable> perUpcColTpc = aod::lrcorrtrktable::upcLRCollisionId;
  Preslice<Ft0aTrksUpcTable> perUpcColFt0a = aod::lrcorrtrktable::upcLRCollisionId;
  Preslice<Ft0cTrksUpcTable> perUpcColFt0c = aod::lrcorrtrktable::upcLRCollisionId;
  Preslice<V0TrksUpcTable> perUpcColV0 = aod::lrcorrtrktable::upcLRCollisionId;

  Preslice<McTrksTable> perMcColTpc = aod::lrcorrmctrktable::lrMcCollisionId;
  Preslice<McFt0aTrksTable> perMcColFt0a = aod::lrcorrmctrktable::lrMcCollisionId;
  Preslice<McFt0cTrksTable> perMcColFt0c = aod::lrcorrmctrktable::lrMcCollisionId;

  void init(InitContext const&)
  {
    std::vector<AxisSpec> corrAxis = {{cfgAxis.axisSample, "Sample"},
                                      {cfgAxis.axisVtxZ, "z-vtx (cm)"},
                                      {cfgAxis.axisMultiplicity, "multiplicity"},
                                      {cfgAxis.axisPtTrigger, "p_{T} (GeV/c)"},
                                      {cfgAxis.axisDeltaPhi, "#Delta#varphi (rad)"},
                                      {cfgAxis.axisDeltaEta, "#Delta#eta"}};
    std::vector<AxisSpec> effAxis = {{cfgAxis.axisVertexEfficiency, "z-vtx (cm)"},
                                     {cfgAxis.axisPtEfficiency, "p_{T} (GeV/c)"},
                                     {cfgAxis.axisEtaEfficiency, "#eta"}};
    std::vector<AxisSpec> userAxis = {{cfgAxis.axisInvMass, "m (GeV/c^2)"}};

    if (!(doprocessTPCtrackEff)) {

      same.setObject(new CorrelationContainer("sameEvent", "sameEvent", corrAxis, effAxis, userAxis));
      mixed.setObject(new CorrelationContainer("mixedEvent", "mixedEvent", corrAxis, effAxis, userAxis));

      histos.add("hMultiplicity", "hMultiplicity", kTH1D, {cfgAxis.axisMultQA});
      histos.add("hCentrality", "hCentrality", kTH1D, {cfgAxis.axisMultQA});
      histos.add("hVertexZ", "hVertexZ", kTH1D, {cfgAxis.axisVtxZ});

      histos.add("hGapSide", "hGapSide", kTH1I, {{5, -0.5, 4.5}});
      histos.add("hTrueGapSide", "hTrueGapSide", kTH1I, {{6, -1.5, 4.5}});
      histos.add("hTrueGapSide_AfterSel", "hTrueGapSide_AfterSel", kTH1I, {{6, -1.5, 4.5}});

      histos.add("Trig_eta", "Trig_eta", kTH1D, {cfgAxis.axisEtaTrig});
      histos.add("Trig_eta_corrected", "Trig_eta_corrected", kTH1D, {cfgAxis.axisEtaTrig});
      histos.add("Trig_phi", "Trig_phi", kTH1D, {cfgAxis.axisPhi});
      histos.add("Trig_phi_corrected", "Trig_phi_corrected", kTH1D, {cfgAxis.axisPhi});
      histos.add("Trig_etavsphi", "Trig_etavsphi", kTH2D, {cfgAxis.axisPhi, cfgAxis.axisEtaTrig});
      histos.add("Trig_pt", "Trig_pt", kTH1D, {cfgAxis.axisPtTrigger});
      histos.add("Trig_pt_corrected", "Trig_pt_corrected", kTH1D, {cfgAxis.axisPtTrigger});
      histos.add("Trig_invMass", "Trig_invMass", kTH1D, {cfgAxis.axisInvMassQA});
      histos.add("Trig_hist", "Trigger Track Properties;Sample;Vertex Z (cm);Multiplicity;p_{T} (GeV/c);Invariant Mass (GeV/c^{2})", kTHnSparseF, {cfgAxis.axisSample, cfgAxis.axisVtxZ, cfgAxis.axisMultiplicity, cfgAxis.axisPtTrigger, cfgAxis.axisInvMass});
      histos.add("Trig_amp", "Trig_amp", kTH1D, {cfgAxis.axisAmplitude});
      histos.add("Channel_vs_Trig_amp", "Channel_vs_Trig_amp", kTH2D, {cfgAxis.axisChannel, cfgAxis.axisAmplitude});

      histos.add("Assoc_eta", "Assoc_eta", kTH1D, {cfgAxis.axisEtaAssoc});
      histos.add("Assoc_phi", "Assoc_phi", kTH1D, {cfgAxis.axisPhi});
      histos.add("Assoc_etavsphi", "Assoc_etavsphi", kTH2D, {cfgAxis.axisPhi, cfgAxis.axisEtaAssoc});
      histos.add("Assoc_amp", "Assoc_amp", kTH1D, {cfgAxis.axisAmplitude});
      histos.add("Channel_vs_Assoc_amp", "Channel_vs_Assoc_amp", kTH2D, {cfgAxis.axisChannel, cfgAxis.axisAmplitude});

      histos.add("deltaEta_deltaPhi_same", "deltaEta_deltaPhi_same", kTH2D, {cfgAxis.axisDeltaPhi, cfgAxis.axisDeltaEta});
      histos.add("deltaEta_deltaPhi_mixed", "deltaEta_deltaPhi_mixed", kTH2D, {cfgAxis.axisDeltaPhi, cfgAxis.axisDeltaEta});

      histos.add("TPCNClsFound", "TPCNClsFound", kTH1D, {cfgAxis.axisTPCNClsFound});
      histos.add("TPCNClsCrossedRows", "TPCNClsCrossedRows", kTH1D, {cfgAxis.axisTPCNClsCrossedRows});
      histos.add("TPCChi2NCl", "TPCChi2NCl", kTH1D, {cfgAxis.axisTPCChi2NCl});
      histos.add("TPCdcaZ", "TPCdcaZ", kTH1D, {cfgAxis.axisTPCdcaZ});
    }

    if (doprocessTPCtrackEff) {
      histos.add("hGenMCdndpt", "hGenMCdndpt", kTH3D, {cfgAxis.axisVtxZ, cfgAxis.axisEtaEfficiency, cfgAxis.axisPtEfficiency});
      histos.add("hRecMCdndpt", "hRecMCdndpt", kTH3D, {cfgAxis.axisVtxZ, cfgAxis.axisEtaEfficiency, cfgAxis.axisPtEfficiency});

      if (cfgSel.cfgPidMask == KPidMaskPion || cfgSel.cfgPidMask == KPidMaskKaon || cfgSel.cfgPidMask == KPidMaskProton) {
        histos.add("hMCGen_PidPtEtaPhi", "MC Gen Target", kTH3D, {cfgAxis.axisMcPt, cfgAxis.axisMcEta, cfgAxis.axisMcPhi});
        histos.add("hMCRec_PidPtEtaPhi", "MC Rec Target", kTH3D, {cfgAxis.axisMcPt, cfgAxis.axisMcEta, cfgAxis.axisMcPhi});
      }

      if (cfgSel.cfgGetdEdx) {
        histos.add("MC_TpcdEdx_ptwise", "MC True TPC dE/dx;Particle Species;p_{T} (GeV/c);TPC dE/dx",
                   kTHnSparseD, {{3, 0.5, 3.5}, cfgAxis.axisPtQA, cfgAxis.axisTpcSignal});
        histos.add("MC_TpcdEdx_TofNsigma_ptwise", "MC True TPC dE/dx vs TOF n#sigma;Particle Species;p_{T} (GeV/c);TPC dE/dx;n#sigma_{TOF}",
                   kTHnSparseD, {{3, 0.5, 3.5}, cfgAxis.axisPtQA, cfgAxis.axisTpcSignal, cfgAxis.axisNsigmaTOF});
        histos.add("MC_ExpSigma_ptwise", "MC Expected TPC #sigma;Particle Species;p_{T} (GeV/c);#sigma_{exp}^{TPC};n#sigma_{TOF}",
                   kTHnSparseD, {{3, 0.5, 3.5}, cfgAxis.axisPtQA, cfgAxis.axisTpcExpSigma, cfgAxis.axisNsigmaTOF});
      }
    }

    myTrackFilter = getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSibAny,
                                                        TrackSelection::GlobalTrackRun3DCAxyCut::Default);
    myTrackFilter.SetPtRange(cfgSel.cfgPtCutMin, cfgSel.cfgPtCutMax);
    myTrackFilter.SetEtaRange(-cfgSel.cfgEtaCut, cfgSel.cfgEtaCut);
    myTrackFilter.SetMinNCrossedRowsTPC(static_cast<int>(cfgSel.cfgTpcMinNCrossedRows));
    myTrackFilter.SetMinNClustersTPC(static_cast<int>(cfgSel.cfgTpcMinNclsFound));
    myTrackFilter.SetMaxChi2PerClusterTPC(cfgSel.cfgTpcMaxChi2PerCluster);
    myTrackFilter.SetMaxDcaZ(cfgSel.cfgTpcMaxDcaZ);
    myTrackFilter.print();

    if (doprocessPIDQA) {
      const AxisSpec axisEvtClass = cfgSel.isUseCentEst ? AxisSpec{cfgAxis.axisCentrality, "Centrality (%)"}
                                                        : AxisSpec{cfgAxis.axisMultiplicity, "Multiplicity (N_{ch})"};
      const AxisSpec axisSpecies{3, 0.5, 3.5, "Particle Species (1=#pi, 2=K, 3=p)"};
      const AxisSpec axisV0CleanType{3, 0.5, 3.5, "V0 Clean Daughter (1=#pi from K_{S}^{0}, 2=p from #Lambda, 3=#pi from #Lambda)"};

      histos.add("hPIDQA_EventCount", "PIDQA Event Matching Counter;;Events", kTH1D, {{4, 0.5, 4.5}});
      histos.get<TH1>(HIST("hPIDQA_EventCount"))->GetXaxis()->SetBinLabel(1, "Raw collisions");
      histos.get<TH1>(HIST("hPIDQA_EventCount"))->GetXaxis()->SetBinLabel(2, "Pass basic EvSel");
      histos.get<TH1>(HIST("hPIDQA_EventCount"))->GetXaxis()->SetBinLabel(3, "Matched to LRCollisions");
      histos.get<TH1>(HIST("hPIDQA_EventCount"))->GetXaxis()->SetBinLabel(4, "Unmatched (skipped)");

      if (cfgSel.cfgGetNsigmaQA && (cfgSel.cfgPidMask == KPidMaskPion || cfgSel.cfgPidMask == KPidMaskKaon || cfgSel.cfgPidMask == KPidMaskProton)) {
        histos.add("TofTpcNsigma_before", "TPC vs TOF n#sigma Before Cuts;Event Class;p_{T} (GeV/c);n#sigma_{TPC};n#sigma_{TOF}",
                   kTHnSparseD, {axisEvtClass, cfgAxis.axisPtQA, cfgAxis.axisNsigmaTPC, cfgAxis.axisNsigmaTOF});

        histos.add("TofTpcNsigma_after", "TPC vs TOF n#sigma After Cuts;Event Class;p_{T} (GeV/c);n#sigma_{TPC};n#sigma_{TOF}",
                   kTHnSparseD, {axisEvtClass, cfgAxis.axisPtQA, cfgAxis.axisNsigmaTPC, cfgAxis.axisNsigmaTOF});

        if (cfgSel.cfgGetdEdx) {
          histos.add("TpcdEdx_ptwise", "TPC dE/dx Before Cuts;Event Class;p_{T} (GeV/c);TPC dE/dx;n#sigma_{TOF}",
                     kTHnSparseD, {axisEvtClass, cfgAxis.axisPtQA, cfgAxis.axisTpcSignal, cfgAxis.axisNsigmaTOF});

          histos.add("TpcdEdx_ptwise_afterCut", "TPC dE/dx After Clean TOF Cuts;Event Class;p_{T} (GeV/c);TPC dE/dx;n#sigma_{TOF}",
                     kTHnSparseD, {axisEvtClass, cfgAxis.axisPtQA, cfgAxis.axisTpcSignal, cfgAxis.axisNsigmaTOF});

          histos.add("ExpTpcdEdx_ptwise", "Expected TPC dE/dx;Event Class;p_{T} (GeV/c);Expected TPC dE/dx;n#sigma_{TOF}",
                     kTHnSparseD, {axisEvtClass, cfgAxis.axisPtQA, cfgAxis.axisTpcSignal, cfgAxis.axisNsigmaTOF});

          histos.add("ExpSigma_ptwise", "Expected TPC #sigma;Event Class;p_{T} (GeV/c);#sigma_{exp}^{TPC};n#sigma_{TOF}",
                     kTHnSparseD, {axisEvtClass, cfgAxis.axisPtQA, cfgAxis.axisTpcExpSigma, cfgAxis.axisNsigmaTOF});
        }
      }

      if (cfgSel.cfgGetAllSpeciesQA && cfgSel.cfgGetdEdx) {
        histos.add("TpcdEdx_allSpecies_ptwise", "All-Species TPC dE/dx;Species (1=#pi,2=K,3=p);Event Class;p_{T} (GeV/c);TPC dE/dx;n#sigma_{TOF}",
                   kTHnSparseD, {axisSpecies, axisEvtClass, cfgAxis.axisPtQA, cfgAxis.axisTpcSignal, cfgAxis.axisNsigmaTOF});

        histos.add("TpcdEdx_allSpecies_afterCut", "All-Species Clean TOF TPC dE/dx;Species (1=#pi,2=K,3=p);Event Class;p_{T} (GeV/c);TPC dE/dx;n#sigma_{TOF}",
                   kTHnSparseD, {axisSpecies, axisEvtClass, cfgAxis.axisPtQA, cfgAxis.axisTpcSignal, cfgAxis.axisNsigmaTOF});

        histos.add("ExpSigma_allSpecies_ptwise", "All-Species Expected TPC #sigma;Species (1=#pi,2=K,3=p);Event Class;p_{T} (GeV/c);#sigma_{exp}^{TPC};n#sigma_{TOF}",
                   kTHnSparseD, {axisSpecies, axisEvtClass, cfgAxis.axisPtQA, cfgAxis.axisTpcExpSigma, cfgAxis.axisNsigmaTOF});
      }

      if (cfgSel.cfgGetV0CleanSample && cfgSel.cfgGetdEdx) {
        histos.add("V0_TpcdEdx_ptwise", "Clean V0 Daughter TPC dE/dx;Daughter Type (1=#pi(K_{S}^{0}), 2=p(#Lambda), 3=#pi(#Lambda));Event Class;p_{T} (GeV/c);TPC dE/dx;n#sigma_{TOF}",
                   kTHnSparseD, {axisV0CleanType, axisEvtClass, cfgAxis.axisPtQA, cfgAxis.axisTpcSignal, cfgAxis.axisNsigmaTOF});

        histos.add("V0_ExpSigma_ptwise", "Clean V0 Daughter Expected TPC #sigma;Daughter Type (1=#pi(K_{S}^{0}), 2=p(#Lambda), 3=#pi(#Lambda));Event Class;p_{T} (GeV/c);#sigma_{exp}^{TPC};n#sigma_{TOF}",
                   kTHnSparseD, {axisV0CleanType, axisEvtClass, cfgAxis.axisPtQA, cfgAxis.axisTpcExpSigma, cfgAxis.axisNsigmaTOF});

        histos.add("hCleanK0s_MassVsPt", "Selected Clean K_{S}^{0};p_{T} (GeV/c);M(#pi^{+}#pi^{-}) (GeV/c^{2})",
                   kTH2D, {cfgAxis.axisPtQA, {100, 0.45, 0.55}});
        histos.add("hCleanLambda_MassVsPt", "Selected Clean #Lambda/#bar{#Lambda};p_{T} (GeV/c);M(p#pi) (GeV/c^{2})",
                   kTH2D, {cfgAxis.axisPtQA, {100, 1.09, 1.14}});
      }
    }
  }

  void loadEffCorrection(uint64_t timestamp)
  {
    if (fLoadTrkEffCorr) {
      return;
    }
    if (!cfgSel.cfgEffccdbPath.value.empty()) {
      hTrkEff = ccdb->getForTimeStamp<TH3D>(cfgSel.cfgEffccdbPath, timestamp);
      if (hTrkEff == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", cfgSel.cfgEffccdbPath.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram from %s (%p)", cfgSel.cfgEffccdbPath.value.c_str(), static_cast<void*>(hTrkEff));
    }

    if (!cfgSel.cfgAccccdbPath.value.empty()) {
      hTrkAcc = ccdb->getForTimeStamp<TH3D>(cfgSel.cfgAccccdbPath, timestamp);
      if (hTrkAcc == nullptr) {
        LOGF(fatal, "Could not load NUA histogram for trigger particles from %s", cfgSel.cfgAccccdbPath.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram from %s (%p)", cfgSel.cfgAccccdbPath.value.c_str(), static_cast<void*>(hTrkAcc));
    }

    fLoadTrkEffCorr = true;
  }

  float getTrkEffCorr(float posZ, float eta, float pt)
  {
    if (!cfgSel.applyEffCorr || !hTrkEff) {
      return 1.0;
    }
    int zBin = hTrkEff->GetXaxis()->FindBin(posZ);
    int etaBin = hTrkEff->GetYaxis()->FindBin(eta);
    int ptBin = hTrkEff->GetZaxis()->FindBin(pt);
    float effweight = 1.0 / hTrkEff->GetBinContent(zBin, etaBin, ptBin);
    if (!std::isfinite(effweight) || effweight <= 0) {
      return 1.0;
    }
    return effweight;
  }

  float getTrkAccCorr(float posZ, float eta, float phi)
  {
    if (!cfgSel.applyAccCorr || !hTrkAcc) {
      return 1.0;
    }
    int zBin = hTrkAcc->GetXaxis()->FindBin(posZ);
    int etaBin = hTrkAcc->GetYaxis()->FindBin(eta);
    int phiBin = hTrkAcc->GetZaxis()->FindBin(phi);
    float nua = hTrkAcc->GetBinContent(zBin, etaBin, phiBin);
    if (!std::isfinite(nua) || nua <= 0) {
      return 1.0;
    }
    return nua;
  }

  template <typename TTrack>
  bool isTrackSelected(TTrack const& track)
  {
    if constexpr (requires { track.tpcNClsFound(); }) {
      if (track.tpcNClsFound() < cfgSel.cfgTpcMinNclsFound) {
        return false;
      }
      if (track.tpcNClsCrossedRows() < cfgSel.cfgTpcMinNCrossedRows) {
        return false;
      }
      if (track.tpcChi2NCl() > cfgSel.cfgTpcMaxChi2PerCluster) {
        return false;
      }
      if (std::abs(track.dcaZ()) > cfgSel.cfgTpcMaxDcaZ) {
        return false;
      }
      return true;
    } else {
      return true;
    }
  }

  template <typename TCollision>
  void fillCollQA(TCollision const& col)
  {
    histos.fill(HIST("hMultiplicity"), col.multiplicity());
    if constexpr (requires { col.centrality(); }) {
      histos.fill(HIST("hCentrality"), col.centrality());
    }
    histos.fill(HIST("hVertexZ"), col.posZ());
  }

  template <typename TTrack>
  void fillTrigTrackQA(TTrack const& track, float trigAmpl, float trkeff, float trkAcc)
  {
    histos.fill(HIST("Trig_etavsphi"), track.phi(), track.eta());
    histos.fill(HIST("Trig_eta"), track.eta());
    histos.fill(HIST("Trig_phi"), track.phi());
    histos.fill(HIST("Trig_eta_corrected"), track.eta(), trkeff);
    histos.fill(HIST("Trig_phi_corrected"), track.phi(), trkAcc);
    if constexpr (requires { track.channelID(); }) {
      histos.fill(HIST("Trig_amp"), trigAmpl);
      histos.fill(HIST("Channel_vs_Trig_amp"), track.channelID(), trigAmpl);
    } else {
      histos.fill(HIST("Trig_pt"), track.pt());
      histos.fill(HIST("Trig_pt_corrected"), track.pt(), trkeff);
    }
    if constexpr (requires { track.invMass(); }) {
      histos.fill(HIST("Trig_invMass"), track.invMass());
    }
    if constexpr (requires { track.tpcNClsFound(); }) {
      histos.fill(HIST("TPCNClsFound"), track.tpcNClsFound());
      histos.fill(HIST("TPCNClsCrossedRows"), track.tpcNClsCrossedRows());
      histos.fill(HIST("TPCChi2NCl"), track.tpcChi2NCl());
      histos.fill(HIST("TPCdcaZ"), track.dcaZ());
    }
  }

  template <typename TTrack>
  void fillAssocTrackQA(TTrack const& track, float assoAmpl)
  {
    histos.fill(HIST("Assoc_etavsphi"), track.phi(), track.eta());
    histos.fill(HIST("Assoc_eta"), track.eta());
    histos.fill(HIST("Assoc_phi"), track.phi());
    if constexpr (requires { track.channelID(); }) {
      histos.fill(HIST("Assoc_amp"), assoAmpl);
      histos.fill(HIST("Channel_vs_Assoc_amp"), track.channelID(), assoAmpl);
    }
  }

  template <bool fillHist = true, typename CheckCol>
  bool isUpcEventSelected(CheckCol const& col)
  {
    if constexpr (fillHist) {
      histos.fill(HIST("hGapSide"), col.gapSide());
    }
    int truegapSide = sgSelector.trueGap(col, cfgFv0Cut, cfgFt0aCut, cfgFt0cCut, cfgZdcCut);
    if constexpr (fillHist) {
      histos.fill(HIST("hTrueGapSide"), truegapSide);
    }
    if (truegapSide != cfgGapSideCut) {
      return false;
    }
    if constexpr (fillHist) {
      histos.fill(HIST("hTrueGapSide_AfterSel"), truegapSide);
    }
    return true;
  }

  template <CorrelationContainer::CFStep step, typename TTarget, typename TTriggers, typename TAssocs>
  void fillCorrHist(TTarget target, TTriggers const& triggers, TAssocs const& assocs, bool mixing, float vz, float multiplicity, float eventWeight)
  {
    int fSampleIndex = static_cast<int>(gRandom->Uniform(0, cfgSel.cfgSampleSize));
    for (auto const& triggerTrack : triggers) {
      auto trigAmpl = 1.0f;
      auto trkeff = 1.0f;
      auto trkAcc = 1.0f;

      if (!isTrackSelected(triggerTrack)) {
        continue;
      }

      if constexpr (requires { triggerTrack.trackType(); }) {
        if (cfgSel.cfgPidMask != 0 && (cfgSel.cfgPidMask & (1u << static_cast<uint32_t>(triggerTrack.trackType()))) == 0u) {
          continue;
        }
      } else if constexpr (requires { triggerTrack.v0Type(); }) {
        if (cfgSel.cfgV0Mask != 0 && (cfgSel.cfgV0Mask & (1u << static_cast<uint32_t>(triggerTrack.v0Type()))) == 0u) {
          continue;
        }
      }

      if constexpr (requires { triggerTrack.channelID(); }) {
        if (cfgSel.cfgRequireFt0aOuterRing && !triggerTrack.isTrackFT0Outer()) {
          continue;
        }
        if (cfgSel.cfgRequireFt0aInnerRing && triggerTrack.isTrackFT0Outer()) {
          continue;
        }
        trigAmpl = triggerTrack.amplitude();
      } else {
        trigAmpl = 1.0;
      }

      if constexpr (step == CorrelationContainer::kCFStepCorrected) {
        if constexpr (requires { triggerTrack.trackType(); }) {
          trkeff = getTrkEffCorr(vz, triggerTrack.eta(), triggerTrack.pt());
          trkAcc = getTrkAccCorr(vz, triggerTrack.eta(), triggerTrack.phi());
        } else {
          trkeff = 1.0;
          trkAcc = 1.0;
        }
      }

      if (cfgSel.cfgVerbosity > 0) {
        LOGF(info, "NUE correction factor: %f | NUA correction factor: %f", trkeff, trkAcc);
      }

      if (!mixing) {
        fillTrigTrackQA(triggerTrack, trigAmpl, trkeff, trkAcc);
        if constexpr (requires { triggerTrack.channelID(); }) {
          histos.fill(HIST("Trig_hist"), fSampleIndex, vz, multiplicity, 1.0, 1.0, eventWeight * trigAmpl * trkeff * trkAcc);
        } else if constexpr (requires { triggerTrack.invMass(); }) {
          histos.fill(HIST("Trig_hist"), fSampleIndex, vz, multiplicity, triggerTrack.pt(), triggerTrack.invMass(), eventWeight * trigAmpl * trkeff * trkAcc);
        } else {
          histos.fill(HIST("Trig_hist"), fSampleIndex, vz, multiplicity, triggerTrack.pt(), 1.0, eventWeight * trigAmpl * trkeff * trkAcc);
        }
      }
      for (auto const& assoTrack : assocs) {
        auto assoAmpl = 1.0f;

        if (!isTrackSelected(assoTrack)) {
          continue;
        }

        if constexpr (requires { assoTrack.channelID(); }) {
          if (cfgSel.cfgRequireFt0cOuterRing && !assoTrack.isTrackFT0Outer()) {
            continue;
          }
          if (cfgSel.cfgRequireFt0cInnerRing && assoTrack.isTrackFT0Outer()) {
            continue;
          }
          assoAmpl = assoTrack.amplitude();
        } else {
          assoAmpl = 1.0f;
        }

        float deltaPhi = RecoDecay::constrainAngle(triggerTrack.phi() - assoTrack.phi(), -PIHalf);
        float deltaEta = triggerTrack.eta() - assoTrack.eta();
        if (!mixing) {
          fillAssocTrackQA(assoTrack, assoAmpl);
          histos.fill(HIST("deltaEta_deltaPhi_same"), deltaPhi, deltaEta, eventWeight * trigAmpl * assoAmpl * trkeff * trkAcc);
        } else {
          histos.fill(HIST("deltaEta_deltaPhi_mixed"), deltaPhi, deltaEta, eventWeight * trigAmpl * assoAmpl * trkeff * trkAcc);
        }
        if constexpr (requires { triggerTrack.channelID(); }) {
          target->getPairHist()->Fill(step, fSampleIndex, vz, multiplicity, 1.0, deltaPhi, deltaEta, 1.0, eventWeight * trigAmpl * assoAmpl * trkeff * trkAcc);
        } else if constexpr (requires { triggerTrack.invMass(); }) {
          target->getPairHist()->Fill(step, fSampleIndex, vz, multiplicity, triggerTrack.pt(), deltaPhi, deltaEta, triggerTrack.invMass(), eventWeight * trigAmpl * assoAmpl * trkeff * trkAcc);
        } else {
          target->getPairHist()->Fill(step, fSampleIndex, vz, multiplicity, triggerTrack.pt(), deltaPhi, deltaEta, 1.0, eventWeight * trigAmpl * assoAmpl * trkeff * trkAcc);
        }
      } // associated tracks
    } // trigger tracks
  } // fill correlation

  template <typename TCollision, typename TTriggers, typename TAssocs>
  void processSame(TCollision const& col, TTriggers const& triggers, TAssocs const& assocs)
  {
    if (std::abs(col.posZ()) >= cfgSel.cfgVtxCut) {
      return;
    }
    loadEffCorrection(col.timestamp());
    fillCollQA(col);
    auto multiplicity = 1.0f;
    if constexpr (requires { col.centrality(); }) {
      if (cfgSel.isUseCentEst) {
        multiplicity = col.centrality();
      } else {
        multiplicity = col.multiplicity();
      }
    } else {
      multiplicity = col.multiplicity();
    }
    if (cfgSel.applyEffCorr) {
      fillCorrHist<CorrelationContainer::kCFStepCorrected>(same, triggers, assocs, false, col.posZ(), multiplicity, 1.0);
    } else {
      fillCorrHist<CorrelationContainer::kCFStepReconstructed>(same, triggers, assocs, false, col.posZ(), multiplicity, 1.0);
    }
  } // process same

  template <typename TCollision, typename... TrackTypes>
  void processMixed(TCollision const& cols, TrackTypes&&... tracks)
  {
    auto getMultiplicity = [this](auto& col) {
      if constexpr (requires { col.gapSide(); }) {
        if (!isUpcEventSelected<false>(col)) {
          return -1.0f;
        }
      } else {
        (void)this;
      }
      auto multiplicity = 1.0f;
      if constexpr (requires { col.centrality(); }) {
        if (cfgSel.isUseCentEst) {
          multiplicity = col.centrality();
        } else {
          multiplicity = col.multiplicity();
        }
      } else {
        multiplicity = col.multiplicity();
      }
      return multiplicity;
    };
    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getMultiplicity)>, aod::collision::PosZ, decltype(getMultiplicity)>;
    MixedBinning binningOnVtxAndMult{{getMultiplicity}, {cfgAxis.axisVtxZME, cfgAxis.axisMultME}, true};
    auto tracksTuple = std::make_tuple(std::forward<TrackTypes>(tracks)...);
    using TupleAtrack = std::tuple_element<0, decltype(tracksTuple)>::type;
    using TupleBtrack = std::tuple_element<std::tuple_size_v<decltype(tracksTuple)> - 1, decltype(tracksTuple)>::type;
    Pair<TCollision, TupleAtrack, TupleBtrack, MixedBinning> pairs{binningOnVtxAndMult, cfgSel.cfgNmixedevent, -1, cols, tracksTuple, &cache};
    for (auto it = pairs.begin(); it != pairs.end(); it++) {
      auto& [col1, tracks1, col2, tracks2] = *it;
      if constexpr (requires { col1.gapSide(); } || requires { col2.gapSide(); }) {
        if (!isUpcEventSelected<false>(col1) || !isUpcEventSelected<false>(col2)) {
          continue;
        }
      }
      if (std::abs(col1.posZ()) >= cfgSel.cfgVtxCut || std::abs(col2.posZ()) >= cfgSel.cfgVtxCut) {
        continue;
      }
      float eventweight = 1.0f / it.currentWindowNeighbours();
      auto multiplicity = getMultiplicity(col1);
      int bin = binningOnVtxAndMult.getBin(std::tuple(col1.posZ(), multiplicity));
      loadEffCorrection(col1.timestamp());

      if (cfgSel.cfgVerbosity > 0) {
        LOGF(info, "processMixed: Mixed collisions bin: %d pair: [%d, %d] %d (%.3f, %.3f), %d (%.3f, %.3f)", bin, it.isNewWindow(), it.currentWindowNeighbours(), col1.globalIndex(), col1.posZ(), col1.multiplicity(), col2.globalIndex(), col2.posZ(), col2.multiplicity());
      }

      if (cfgSel.applyEffCorr) {
        fillCorrHist<CorrelationContainer::kCFStepCorrected>(mixed, tracks1, tracks2, true, col1.posZ(), multiplicity, eventweight);
      } else {
        fillCorrHist<CorrelationContainer::kCFStepReconstructed>(mixed, tracks1, tracks2, true, col1.posZ(), multiplicity, eventweight);
      }
    } // pair loop
  } // process mixed

  template <typename TTriggers, typename TAssocs>
  void processMcSame(McCollsTable::iterator const& mccollision, soa::SmallGroups<aod::LRCollisionsWithLabel> const& collisions, TTriggers const& triggers, TAssocs const& assocs)
  {
    if (std::abs(mccollision.posZ()) >= cfgSel.cfgVtxCut) {
      return;
    }
    fillCollQA(mccollision);
    auto multiplicity = mccollision.multiplicity();
    if (cfgSel.isUseDataLikeMult > 0) {
      for (const auto& collision : collisions) {
        if (cfgSel.isUseCentEst) {
          multiplicity = collision.centrality();
        } else {
          multiplicity = collision.multiplicity();
        }
      }
    }
    fillCorrHist<CorrelationContainer::kCFStepAll>(same, triggers, assocs, false, mccollision.posZ(), multiplicity, 1.0);
  } // process MC same

  template <typename... TrackTypes>
  void processMcMixed(McCollsTable const& mccollisions, aod::LRCollisionsWithLabel const& collisions, TrackTypes&&... tracks)
  {
    bool useMCMultiplicity = (cfgSel.isUseDataLikeMult == 0);
    auto getMultiplicity =
      [&collisions, &useMCMultiplicity, this](auto& col) {
        if (useMCMultiplicity) {
          return col.multiplicity();
        }
        auto groupedCollisions = collisions.sliceByCached(aod::lrcorrcolltable::lrMcCollisionId, col.globalIndex(), this->cache);
        if (groupedCollisions.size() == 0) {
          return -1.0f;
        }
        if (cfgSel.isUseCentEst) {
          return groupedCollisions.begin().centrality();
        }
        return groupedCollisions.begin().multiplicity();
      };

    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getMultiplicity)>, aod::mccollision::PosZ, decltype(getMultiplicity)>;
    MixedBinning binningOnVtxAndMult{{getMultiplicity}, {cfgAxis.axisVtxZME, cfgAxis.axisMultME}, true};
    auto tracksTuple = std::make_tuple(std::forward<TrackTypes>(tracks)...);
    using TupleAtrack = std::tuple_element<0, decltype(tracksTuple)>::type;
    using TupleBtrack = std::tuple_element<std::tuple_size_v<decltype(tracksTuple)> - 1, decltype(tracksTuple)>::type;
    Pair<McCollsTable, TupleAtrack, TupleBtrack, MixedBinning> pairs{binningOnVtxAndMult, cfgSel.cfgNmixedevent, -1, mccollisions, tracksTuple, &cache};
    for (auto it = pairs.begin(); it != pairs.end(); it++) {
      auto& [col1, tracks1, col2, tracks2] = *it;
      float eventweight = 1.0f / it.currentWindowNeighbours();
      auto multiplicity = getMultiplicity(col1);
      fillCorrHist<CorrelationContainer::kCFStepAll>(mixed, tracks1, tracks2, true, col1.posZ(), multiplicity, eventweight);
    } // pair loop
  } // process MC mixed

  template <typename TTriggers, typename TAssocs>
  void processMcGenSame(McCollsTable::iterator const& mccollision, TTriggers const& triggers, TAssocs const& assocs)
  {
    if (std::abs(mccollision.posZ()) >= cfgSel.cfgVtxCut) {
      return;
    }
    fillCollQA(mccollision);
    auto multiplicity = mccollision.multiplicity();
    fillCorrHist<CorrelationContainer::kCFStepAll>(same, triggers, assocs, false, mccollision.posZ(), multiplicity, 1.0);
  } // process MC gen same

  template <typename... TrackTypes>
  void processMcGenMixed(McCollsTable const& mccollisions, TrackTypes&&... tracks)
  {
    auto getMultiplicity = [this](auto& collision) {
      (void)this;
      return collision.multiplicity();
    };
    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getMultiplicity)>, aod::mccollision::PosZ, decltype(getMultiplicity)>;
    MixedBinning binningOnVtxAndMult{{getMultiplicity}, {cfgAxis.axisVtxZME, cfgAxis.axisMultME}, true};
    auto tracksTuple = std::make_tuple(std::forward<TrackTypes>(tracks)...);
    using TupleAtrack = std::tuple_element<0, decltype(tracksTuple)>::type;
    using TupleBtrack = std::tuple_element<std::tuple_size_v<decltype(tracksTuple)> - 1, decltype(tracksTuple)>::type;
    Pair<McCollsTable, TupleAtrack, TupleBtrack, MixedBinning> pairs{binningOnVtxAndMult, cfgSel.cfgNmixedevent, -1, mccollisions, tracksTuple, &cache};
    for (auto it = pairs.begin(); it != pairs.end(); it++) {
      auto& [col1, tracks1, col2, tracks2] = *it;
      float eventweight = 1.0f / it.currentWindowNeighbours();
      auto multiplicity = getMultiplicity(col1);
      fillCorrHist<CorrelationContainer::kCFStepAll>(mixed, tracks1, tracks2, true, col1.posZ(), multiplicity, eventweight);
    } // pair loop
  } // process MC gen mixed

  void processTpcft0aSE(CollsTable::iterator const& col, TrksTable const& tracks, Ft0aTrksTable const& ft0as)
  {
    processSame(col, tracks, ft0as);
  }

  void processTpcft0cSE(CollsTable::iterator const& col, TrksTable const& tracks, Ft0cTrksTable const& ft0cs)
  {
    processSame(col, tracks, ft0cs);
  }

  void processV0ft0aSE(CollsTable::iterator const& col, V0TrksTable const& tracks, Ft0aTrksTable const& ft0as)
  {
    processSame(col, tracks, ft0as);
  }

  void processV0ft0cSE(CollsTable::iterator const& col, V0TrksTable const& tracks, Ft0cTrksTable const& ft0cs)
  {
    processSame(col, tracks, ft0cs);
  }

  void processV0tpcSE(CollsTable::iterator const& col, V0TrksTable const& v0s, TrksTable const& tracks)
  {
    processSame(col, v0s, tracks);
  }

  void processFt0aft0cSE(CollsTable::iterator const& col, Ft0aTrksTable const& ft0as, Ft0cTrksTable const& ft0cs)
  {
    processSame(col, ft0as, ft0cs);
  }

  void processTpcft0aME(CollsTable const& cols, TrksTable const& tracks, Ft0aTrksTable const& ft0as)
  {
    processMixed(cols, tracks, ft0as);
  }

  void processTpcft0cME(CollsTable const& cols, TrksTable const& tracks, Ft0cTrksTable const& ft0cs)
  {
    processMixed(cols, tracks, ft0cs);
  }

  void processV0ft0aME(CollsTable const& cols, V0TrksTable const& tracks, Ft0aTrksTable const& ft0as)
  {
    processMixed(cols, tracks, ft0as);
  }

  void processV0ft0cME(CollsTable const& cols, V0TrksTable const& tracks, Ft0cTrksTable const& ft0cs)
  {
    processMixed(cols, tracks, ft0cs);
  }

  void processV0tpcME(CollsTable const& cols, V0TrksTable const& v0s, TrksTable const& tracks)
  {
    processMixed(cols, v0s, tracks);
  }

  void processFt0aft0cME(CollsTable const& cols, Ft0aTrksTable const& ft0as, Ft0cTrksTable const& ft0cs)
  {
    processMixed(cols, ft0as, ft0cs);
  }

  void processUpcTpcft0aSE(UpcCollsTable::iterator const& col, TrksUpcTable const& tracks, Ft0aTrksUpcTable const& ft0as)
  {
    if (!isUpcEventSelected<true>(col)) {
      return;
    }
    processSame(col, tracks, ft0as);
  }

  void processUpcTpcft0cSE(UpcCollsTable::iterator const& col, TrksUpcTable const& tracks, Ft0cTrksUpcTable const& ft0cs)
  {
    if (!isUpcEventSelected<true>(col)) {
      return;
    }
    processSame(col, tracks, ft0cs);
  }

  void processUpcV0ft0aSE(UpcCollsTable::iterator const& col, V0TrksUpcTable const& tracks, Ft0aTrksUpcTable const& ft0as)
  {
    if (!isUpcEventSelected<true>(col)) {
      return;
    }
    processSame(col, tracks, ft0as);
  }

  void processUpcV0ft0cSE(UpcCollsTable::iterator const& col, V0TrksUpcTable const& tracks, Ft0cTrksUpcTable const& ft0cs)
  {
    processSame(col, tracks, ft0cs);
  }

  void processUpcTpcft0aME(UpcCollsTable const& cols, TrksUpcTable const& tracks, Ft0aTrksUpcTable const& ft0as)
  {
    processMixed(cols, tracks, ft0as);
  }

  void processUpcTpcft0cME(UpcCollsTable const& cols, TrksUpcTable const& tracks, Ft0cTrksUpcTable const& ft0cs)
  {
    processMixed(cols, tracks, ft0cs);
  }

  void processUpcV0ft0aME(UpcCollsTable const& cols, V0TrksUpcTable const& tracks, Ft0aTrksUpcTable const& ft0as)
  {
    processMixed(cols, tracks, ft0as);
  }

  void processUpcV0ft0cME(UpcCollsTable const& cols, V0TrksUpcTable const& tracks, Ft0cTrksUpcTable const& ft0cs)
  {
    processMixed(cols, tracks, ft0cs);
  }

  void processMcTpcft0aSE(McCollsTable::iterator const& mccollision, soa::SmallGroups<aod::LRCollisionsWithLabel> const& collisions, McTrksTable const& tracks, McFt0aTrksTable const& ft0as)
  {
    processMcSame(mccollision, collisions, tracks, ft0as);
  }

  void processMcTpcft0cSE(McCollsTable::iterator const& mccollision, soa::SmallGroups<aod::LRCollisionsWithLabel> const& collisions, McTrksTable const& tracks, McFt0cTrksTable const& ft0cs)
  {
    processMcSame(mccollision, collisions, tracks, ft0cs);
  }

  void processMcFt0aft0cSE(McCollsTable::iterator const& mccollision, soa::SmallGroups<aod::LRCollisionsWithLabel> const& collisions, McFt0aTrksTable const& ft0as, McFt0cTrksTable const& ft0cs)
  {
    processMcSame(mccollision, collisions, ft0as, ft0cs);
  }

  void processMcTpcft0aME(McCollsTable const& mccollisions, aod::LRCollisionsWithLabel const& collisions, McTrksTable const& tracks, McFt0aTrksTable const& ft0as)
  {
    processMcMixed(mccollisions, collisions, tracks, ft0as);
  }

  void processMcTpcft0cME(McCollsTable const& mccollisions, aod::LRCollisionsWithLabel const& collisions, McTrksTable const& tracks, McFt0cTrksTable const& ft0cs)
  {
    processMcMixed(mccollisions, collisions, tracks, ft0cs);
  }

  void processMcFt0aft0cME(McCollsTable const& mccollisions, aod::LRCollisionsWithLabel const& collisions, McFt0aTrksTable const& ft0as, McFt0cTrksTable const& ft0cs)
  {
    processMcMixed(mccollisions, collisions, ft0as, ft0cs);
  }

  void processMcGenTpcft0aSE(McCollsTable::iterator const& mccollision, McTrksTable const& tracks, McFt0aTrksTable const& ft0as)
  {
    processMcGenSame(mccollision, tracks, ft0as);
  }

  void processMcGenTpcft0cSE(McCollsTable::iterator const& mccollision, McTrksTable const& tracks, McFt0cTrksTable const& ft0cs)
  {
    processMcGenSame(mccollision, tracks, ft0cs);
  }

  void processMcGenV0ft0aSE(McCollsTable::iterator const& mccollision, McTrksTable const& tracks, McFt0aTrksTable const& ft0as)
  {
    processMcGenSame(mccollision, tracks, ft0as);
  }

  void processMcGenV0ft0cSE(McCollsTable::iterator const& mccollision, McTrksTable const& tracks, McFt0cTrksTable const& ft0cs)
  {
    processMcGenSame(mccollision, tracks, ft0cs);
  }

  void processMcGenFt0aft0cSE(McCollsTable::iterator const& mccollision, McFt0aTrksTable const& ft0as, McFt0cTrksTable const& ft0cs)
  {
    processMcGenSame(mccollision, ft0as, ft0cs);
  }

  void processMcGenTpcft0aME(McCollsTable const& mccollisions, McTrksTable const& tracks, McFt0aTrksTable const& ft0as)
  {
    processMcGenMixed(mccollisions, tracks, ft0as);
  }

  void processMcGenTpcft0cME(McCollsTable const& mccollisions, McTrksTable const& tracks, McFt0cTrksTable const& ft0cs)
  {
    processMcGenMixed(mccollisions, tracks, ft0cs);
  }

  void processMcGenV0ft0aME(McCollsTable const& mccollisions, McTrksTable const& tracks, McFt0aTrksTable const& ft0as)
  {
    processMcGenMixed(mccollisions, tracks, ft0as);
  }

  void processMcGenV0ft0cME(McCollsTable const& mccollisions, McTrksTable const& tracks, McFt0cTrksTable const& ft0cs)
  {
    processMcGenMixed(mccollisions, tracks, ft0cs);
  }

  void processMcGenFt0aft0cME(McCollsTable const& mccollisions, McFt0aTrksTable const& ft0as, McFt0cTrksTable const& ft0cs)
  {
    processMcGenMixed(mccollisions, ft0as, ft0cs);
  }

  using ColMCTrueTable = soa::Join<aod::McCollisions, aod::McCollsExtra>;
  using ColMCRecTable = soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::EvSels>>;
  using TrksMCRecTable = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFbeta, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
  Preslice<TrksMCRecTable> perColMidtrack = aod::track::collisionId;

  template <typename CheckCol>
  bool isEventSelected(CheckCol const& col)
  {
    if (!col.sel8()) {
      return false;
    }
    if (!col.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if (!col.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if (!col.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      return false;
    }
    if (std::abs(col.posZ()) >= cfgSel.cfgVtxCut) {
      return false;
    }
    return true;
  }

  template <typename CheckGenPart>
  bool isGenPartSelected(CheckGenPart const& particle)
  {
    if (!particle.isPhysicalPrimary() || !particle.producedByGenerator()) {
      return false;
    }
    auto pdgParticle = pdg->GetParticle(particle.pdgCode());
    if (pdgParticle == nullptr) {
      return false;
    }
    if (std::abs(pdgParticle->Charge()) < KminCharge) {
      return false;
    }
    return true;
  }

  void processTPCtrackEff(ColMCTrueTable::iterator const& mcCollision, ColMCRecTable const& RecCols,
                          TrksMCRecTable const& RecTracks, aod::McParticles const& mcparticles)
  {
    if (std::abs(mcCollision.posZ()) >= cfgSel.cfgVtxCut) {
      return;
    }
    bool atLeastOne = false;
    for (const auto& RecCol : RecCols) {
      if (!isEventSelected(RecCol)) {
        continue;
      }
      if (RecCol.globalIndex() != mcCollision.bestCollisionIndex()) {
        continue;
      }
      atLeastOne = true;
    }
    for (const auto& particle : mcparticles) {
      if (!isGenPartSelected(particle) ||
          std::abs(particle.eta()) > cfgSel.cfgEtaCut ||
          particle.pt() < cfgSel.cfgPtCutMin ||
          particle.pt() > cfgSel.cfgPtCutMax) {
        continue;
      }
      if (atLeastOne) {
        histos.fill(HIST("hGenMCdndpt"), mcCollision.posZ(), particle.eta(), particle.pt());

        if (cfgSel.cfgPidMask == KPidMaskPion || cfgSel.cfgPidMask == KPidMaskKaon || cfgSel.cfgPidMask == KPidMaskProton) {
          auto pdgcode = std::abs(particle.pdgCode());

          bool isTargetGen = (cfgSel.cfgPidMask == KPidMaskPion && pdgcode == PDG_t::kPiPlus) ||
                             (cfgSel.cfgPidMask == KPidMaskKaon && pdgcode == PDG_t::kKPlus) ||
                             (cfgSel.cfgPidMask == KPidMaskProton && pdgcode == PDG_t::kProton);

          if (isTargetGen) {
            histos.fill(HIST("hMCGen_PidPtEtaPhi"), particle.pt(), particle.eta(), particle.phi());
          }
        }
      }
    }

    for (const auto& RecCol : RecCols) {
      if (!isEventSelected(RecCol)) {
        continue;
      }
      if (RecCol.globalIndex() != mcCollision.bestCollisionIndex()) {
        continue;
      }
      auto recTracksPart = RecTracks.sliceBy(perColMidtrack, RecCol.globalIndex());
      for (const auto& track : recTracksPart) {
        if (!track.isGlobalTrack()) {
          continue;
        }
        if (!myTrackFilter.IsSelected(track)) {
          continue;
        }
        if (!track.has_mcParticle()) {
          continue;
        }
        auto particle = track.mcParticle();
        if (RecCol.mcCollisionId() != particle.mcCollisionId()) {
          continue;
        }
        if (particle.isPhysicalPrimary()) {
          histos.fill(HIST("hRecMCdndpt"), mcCollision.posZ(), particle.eta(), particle.pt());

          if (cfgSel.cfgGetdEdx) {
            int pdgIdx = 0;
            float tofNsigMC = -999.0f;
            float expSigmaMC = 0.0f;
            auto absPdg = std::abs(particle.pdgCode());

            if (absPdg == PDG_t::kPiPlus) {
              pdgIdx = 1; // o2-linter: disable=pdg/explicit-code (histogram species index)
              tofNsigMC = track.tofNSigmaPi();
              expSigmaMC = track.tpcExpSigmaPi();
            } else if (absPdg == PDG_t::kKPlus) {
              pdgIdx = 2; // o2-linter: disable=pdg/explicit-code (histogram species index)
              tofNsigMC = track.tofNSigmaKa();
              expSigmaMC = track.tpcExpSigmaKa();
            } else if (absPdg == PDG_t::kProton) {
              pdgIdx = 3; // o2-linter: disable=pdg/explicit-code (histogram species index)
              tofNsigMC = track.tofNSigmaPr();
              expSigmaMC = track.tpcExpSigmaPr();
            }

            if (pdgIdx > 0) {
              if (track.hasTOF()) {
                histos.fill(HIST("MC_TpcdEdx_TofNsigma_ptwise"), static_cast<double>(pdgIdx), track.pt(), track.tpcSignal(), tofNsigMC);
                histos.fill(HIST("MC_ExpSigma_ptwise"), static_cast<double>(pdgIdx), track.pt(), expSigmaMC, tofNsigMC);
              }
              bool passTofForMC = !cfgSel.cfgRequireTofInMCdEdx ||
                                  (track.hasTOF() && tofNsigMC > cfgSel.cfgCleanNsigmaTofMin && tofNsigMC < cfgSel.cfgCleanNsigmaTofMax);
              if (passTofForMC) {
                histos.fill(HIST("MC_TpcdEdx_ptwise"), static_cast<double>(pdgIdx), track.pt(), track.tpcSignal());
              }
            }
          }

          bool isTpcPion = (track.tpcNSigmaPi() > cfgSel.cfgPidNsigmaMin && track.tpcNSigmaPi() < cfgSel.cfgPidNsigmaMax);
          bool isTpcKaon = (track.tpcNSigmaKa() > cfgSel.cfgPidNsigmaMin && track.tpcNSigmaKa() < cfgSel.cfgPidNsigmaMax);
          bool isTpcProton = (track.tpcNSigmaPr() > cfgSel.cfgPidNsigmaMin && track.tpcNSigmaPr() < cfgSel.cfgPidNsigmaMax);

          bool isTofPion = (track.tofNSigmaPi() > cfgSel.cfgPidNsigmaMin && track.tofNSigmaPi() < cfgSel.cfgPidNsigmaMax);
          bool isTofKaon = (track.tofNSigmaKa() > cfgSel.cfgPidNsigmaMin && track.tofNSigmaKa() < cfgSel.cfgPidNsigmaMax);
          bool isTofProton = (track.tofNSigmaPr() > cfgSel.cfgPidNsigmaMin && track.tofNSigmaPr() < cfgSel.cfgPidNsigmaMax);

          bool isPion = false, isKaon = false, isProton = false;

          if (track.pt() > cfgSel.cfgTofPidPtCut && track.hasTOF()) {
            isPion = isTofPion && isTpcPion;
            isKaon = isTofKaon && isTpcKaon;
            isProton = isTofProton && isTpcProton;
          } else if (track.pt() <= cfgSel.cfgTofPidPtCut || track.hasTOF()) {
            isPion = isTpcPion;
            isKaon = isTpcKaon;
            isProton = isTpcProton;
          }

          if ((isPion && isKaon) || (isPion && isProton) || (isKaon && isProton)) {
            isPion = isKaon = isProton = false;
          }

          if (cfgSel.cfgPidMask == KPidMaskPion || cfgSel.cfgPidMask == KPidMaskKaon || cfgSel.cfgPidMask == KPidMaskProton) {

            auto pdgcodeRec = std::abs(particle.pdgCode());
            bool isTargetRec = (cfgSel.cfgPidMask == KPidMaskPion && isPion && pdgcodeRec == PDG_t::kPiPlus) ||
                               (cfgSel.cfgPidMask == KPidMaskKaon && isKaon && pdgcodeRec == PDG_t::kKPlus) ||
                               (cfgSel.cfgPidMask == KPidMaskProton && isProton && pdgcodeRec == PDG_t::kProton);

            if (isTargetRec) {
              histos.fill(HIST("hMCRec_PidPtEtaPhi"), particle.pt(), particle.eta(), particle.phi());
            }
          }
        }
      }
    }
  }

  using CollStandardTable = soa::Join<aod::Collisions, aod::EvSels>;
  using TrksStandardTable = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFbeta, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;

  template <typename TTracks>
  float countFallbackMultiplicity(TTracks const& tracks, float vz)
  {
    float nTrk = 0.0f;
    for (auto const& trk : tracks) {
      if (!trk.isGlobalTrack()) {
        continue;
      }
      if (!myTrackFilter.IsSelected(trk)) {
        continue;
      }
      if (trk.pt() < cfgSel.cfgPtCutMinForMult || trk.pt() > cfgSel.cfgPtCutMaxForMult) {
        continue;
      }
      nTrk += getTrkEffCorr(vz, trk.eta(), trk.pt());
    }
    return nTrk;
  }

  template <typename TCollision, typename TV0, typename TTrack>
  void fillCleanV0Daughters(TCollision const& rawCol, TV0 const& v0, TTrack const& posTrack, TTrack const& negTrack, float multiplicity)
  {
    if (!posTrack.hasTPC() || !negTrack.hasTPC()) {
      return;
    }
    if (posTrack.tpcNClsCrossedRows() < cfgSel.cfgTpcMinNCrossedRows || negTrack.tpcNClsCrossedRows() < cfgSel.cfgTpcMinNCrossedRows) {
      return;
    }
    if (std::abs(v0.positiveeta()) > cfgSel.cfgEtaCut || std::abs(v0.negativeeta()) > cfgSel.cfgEtaCut) {
      return;
    }
    if (v0.v0cosPA() < cfgSel.cfgCleanV0MinCosPA) {
      return;
    }
    if (v0.dcaV0daughters() > cfgSel.cfgCleanV0MaxDcaDau) {
      return;
    }
    if (v0.v0radius() < cfgSel.cfgCleanV0MinRadius || v0.v0radius() > cfgSel.cfgCleanV0MaxRadius) {
      return;
    }

    float distOverP = v0.distovertotmom(rawCol.posX(), rawCol.posY(), rawCol.posZ());

    auto fillV0Probe = [&](double v0CleanType, TTrack const& probeTrack, float tpcNsig, float tofNsigReal, float expSigma, bool isPionProbe, bool isProtonProbe) {
      if (!probeTrack.isGlobalTrack() || !myTrackFilter.IsSelected(probeTrack)) {
        return;
      }

      const float tofVal = probeTrack.hasTOF() ? tofNsigReal : -999.0f;

      if (cfgSel.cfgGetV0CleanSample && cfgSel.cfgGetdEdx) {
        histos.fill(HIST("V0_TpcdEdx_ptwise"), v0CleanType, multiplicity, probeTrack.pt(), probeTrack.tpcSignal(), tofVal);
        histos.fill(HIST("V0_ExpSigma_ptwise"), v0CleanType, multiplicity, probeTrack.pt(), expSigma, tofVal);
      }

      if (cfgSel.cfgFillAfterCutWithV0 && cfgSel.cfgGetNsigmaQA) {
        const bool passTofForAfterCut = !cfgSel.cfgRequireTofForV0AfterCut ||
                                        (probeTrack.hasTOF() && tofNsigReal > cfgSel.cfgCleanNsigmaTofMin && tofNsigReal < cfgSel.cfgCleanNsigmaTofMax);
        const bool matchWagon = (((cfgSel.cfgPidMask & KPidMaskPion) != 0) && isPionProbe) ||
                                (((cfgSel.cfgPidMask & KPidMaskProton) != 0) && isProtonProbe);
        if (matchWagon && passTofForAfterCut) {
          histos.fill(HIST("TofTpcNsigma_after"), multiplicity, probeTrack.pt(), tpcNsig, tofVal);
          if (cfgSel.cfgGetdEdx) {
            histos.fill(HIST("TpcdEdx_ptwise_afterCut"), multiplicity, probeTrack.pt(), probeTrack.tpcSignal(), tofVal);
          }
        }
      }
    };

    float ctauK0s = distOverP * o2::constants::physics::MassK0;
    bool isK0sTopo = (std::abs(v0.mK0Short() - o2::constants::physics::MassK0) < cfgSel.cfgCleanK0sMassWin) &&
                     (v0.qtarm() >= cfgSel.cfgCleanK0sMinArmQt * std::abs(v0.alpha())) &&
                     (std::abs(ctauK0s) <= cfgSel.cfgCleanK0sMaxCtau) &&
                     (std::abs(v0.dcapostopv()) >= cfgSel.cfgCleanK0sMinDcaToPv) &&
                     (std::abs(v0.dcanegtopv()) >= cfgSel.cfgCleanK0sMinDcaToPv);

    if (isK0sTopo) {
      if (cfgSel.cfgGetV0CleanSample && cfgSel.cfgGetdEdx) {
        histos.fill(HIST("hCleanK0s_MassVsPt"), v0.pt(), v0.mK0Short());
      }
      if (std::abs(negTrack.tpcNSigmaPi()) < cfgSel.cfgCleanV0TagTpcNsigma) {
        fillV0Probe(1.0, posTrack, posTrack.tpcNSigmaPi(), posTrack.tofNSigmaPi(), posTrack.tpcExpSigmaPi(), true, false);
      }
      if (std::abs(posTrack.tpcNSigmaPi()) < cfgSel.cfgCleanV0TagTpcNsigma) {
        fillV0Probe(1.0, negTrack, negTrack.tpcNSigmaPi(), negTrack.tofNSigmaPi(), negTrack.tpcExpSigmaPi(), true, false);
      }
    }

    float ctauLambda = distOverP * o2::constants::physics::MassLambda;
    bool passK0sRej = std::abs(v0.mK0Short() - o2::constants::physics::MassK0) >= cfgSel.cfgCleanLambdaRejK0sWin;
    bool passLambdaArm = v0.qtarm() >= cfgSel.cfgCleanLambdaMinArmQt;
    bool passLambdaCtau = std::abs(ctauLambda) <= cfgSel.cfgCleanLambdaMaxCtau;

    if (passK0sRej && passLambdaArm && passLambdaCtau) {
      bool isLambdaTopo = (std::abs(v0.mLambda() - o2::constants::physics::MassLambda) < cfgSel.cfgCleanLambdaMassWin) &&
                          (std::abs(v0.dcapostopv()) >= cfgSel.cfgCleanLambdaMinDcaPr) &&
                          (std::abs(v0.dcanegtopv()) >= cfgSel.cfgCleanLambdaMinDcaPi);
      if (isLambdaTopo) {
        if (cfgSel.cfgGetV0CleanSample && cfgSel.cfgGetdEdx) {
          histos.fill(HIST("hCleanLambda_MassVsPt"), v0.pt(), v0.mLambda());
        }
        if (std::abs(negTrack.tpcNSigmaPi()) < cfgSel.cfgCleanV0TagTpcNsigma) {
          fillV0Probe(2.0, posTrack, posTrack.tpcNSigmaPr(), posTrack.tofNSigmaPr(), posTrack.tpcExpSigmaPr(), false, true);
        }
        if (std::abs(posTrack.tpcNSigmaPr()) < cfgSel.cfgCleanV0TagTpcNsigma) {
          fillV0Probe(3.0, negTrack, negTrack.tpcNSigmaPi(), negTrack.tofNSigmaPi(), negTrack.tpcExpSigmaPi(), cfgSel.cfgIncludeLambdaPiInCleanPi, false);
        }
      }

      bool isAntiLambdaTopo = (std::abs(v0.mAntiLambda() - o2::constants::physics::MassLambda) < cfgSel.cfgCleanLambdaMassWin) &&
                              (std::abs(v0.dcapostopv()) >= cfgSel.cfgCleanLambdaMinDcaPi) &&
                              (std::abs(v0.dcanegtopv()) >= cfgSel.cfgCleanLambdaMinDcaPr);
      if (isAntiLambdaTopo) {
        if (cfgSel.cfgGetV0CleanSample && cfgSel.cfgGetdEdx) {
          histos.fill(HIST("hCleanLambda_MassVsPt"), v0.pt(), v0.mAntiLambda());
        }
        if (std::abs(posTrack.tpcNSigmaPi()) < cfgSel.cfgCleanV0TagTpcNsigma) {
          fillV0Probe(2.0, negTrack, negTrack.tpcNSigmaPr(), negTrack.tofNSigmaPr(), negTrack.tpcExpSigmaPr(), false, true);
        }
        if (std::abs(negTrack.tpcNSigmaPr()) < cfgSel.cfgCleanV0TagTpcNsigma) {
          fillV0Probe(3.0, posTrack, posTrack.tpcNSigmaPi(), posTrack.tofNSigmaPi(), posTrack.tpcExpSigmaPi(), cfgSel.cfgIncludeLambdaPiInCleanPi, false);
        }
      }
    }
  }

  void processPIDQA(CollStandardTable const& rawCollisions, CollsTable const& lrCollisions,
                    TrksStandardTable const& standardTracks, aod::V0Datas const& V0s)
  {
    int64_t lrIdx = 0;
    const int64_t nLrCols = lrCollisions.size();

    for (auto const& rawCol : rawCollisions) {
      histos.fill(HIST("hPIDQA_EventCount"), 1.0);

      bool matched = false;
      float multiplicity = -1.0f;
      if (nLrCols > 0 && lrIdx < nLrCols) {
        if (std::abs(lrCollisions.iteratorAt(lrIdx).posZ() - rawCol.posZ()) <= 1e-5f) {
          auto const& lrCol = lrCollisions.iteratorAt(lrIdx);
          multiplicity = cfgSel.isUseCentEst ? lrCol.centrality() : lrCol.multiplicity();
          lrIdx++;
          matched = true;
        }
      }

      if (!isEventSelected(rawCol)) {
        continue;
      }
      histos.fill(HIST("hPIDQA_EventCount"), 2.0);

      if (matched) {
        histos.fill(HIST("hPIDQA_EventCount"), 3.0);
      } else if (nLrCols == 0 && cfgSel.cfgAllowFallbackMultInPIDQA) {
        auto tracksForMult = standardTracks.sliceByCached(aod::track::collisionId, rawCol.globalIndex(), this->cache);
        multiplicity = countFallbackMultiplicity(tracksForMult, rawCol.posZ());
        histos.fill(HIST("hPIDQA_EventCount"), 3.0);
      } else {
        histos.fill(HIST("hPIDQA_EventCount"), 4.0);
        continue;
      }

      auto tracks = standardTracks.sliceByCached(aod::track::collisionId, rawCol.globalIndex(), this->cache);

      for (auto const& track : tracks) {
        if (!track.isGlobalTrack()) {
          continue;
        }
        if (!myTrackFilter.IsSelected(track)) {
          continue;
        }
        if (!track.hasTOF()) {
          continue;
        }

        const float tpcNsigPi = track.tpcNSigmaPi(), tofNsigPi = track.tofNSigmaPi(), expSigPi = track.tpcExpSigmaPi();
        const float tpcNsigKa = track.tpcNSigmaKa(), tofNsigKa = track.tofNSigmaKa(), expSigKa = track.tpcExpSigmaKa();
        const float tpcNsigPr = track.tpcNSigmaPr(), tofNsigPr = track.tofNSigmaPr(), expSigPr = track.tpcExpSigmaPr();

        const bool isTpcCleanPi = (tpcNsigPi > cfgSel.cfgCleanNsigmaTpcMin && tpcNsigPi < cfgSel.cfgCleanNsigmaTpcMax);
        const bool isTpcCleanKa = (tpcNsigKa > cfgSel.cfgCleanNsigmaTpcMin && tpcNsigKa < cfgSel.cfgCleanNsigmaTpcMax);
        const bool isTpcCleanPr = (tpcNsigPr > cfgSel.cfgCleanNsigmaTpcMin && tpcNsigPr < cfgSel.cfgCleanNsigmaTpcMax);

        const bool isTofCleanPi = (tofNsigPi > cfgSel.cfgCleanNsigmaTofMin && tofNsigPi < cfgSel.cfgCleanNsigmaTofMax) &&
                                  (std::abs(tofNsigKa) > cfgSel.cfgCleanTofVetoNsigma) &&
                                  (std::abs(tofNsigPr) > cfgSel.cfgCleanTofVetoNsigma);
        const bool isTofCleanKa = (tofNsigKa > cfgSel.cfgCleanNsigmaTofMin && tofNsigKa < cfgSel.cfgCleanNsigmaTofMax) &&
                                  (std::abs(tofNsigPi) > cfgSel.cfgCleanTofVetoNsigma) &&
                                  (std::abs(tofNsigPr) > cfgSel.cfgCleanTofVetoNsigma);
        const bool isTofCleanPr = (tofNsigPr > cfgSel.cfgCleanNsigmaTofMin && tofNsigPr < cfgSel.cfgCleanNsigmaTofMax) &&
                                  (std::abs(tofNsigPi) > cfgSel.cfgCleanTofVetoNsigma) &&
                                  (std::abs(tofNsigKa) > cfgSel.cfgCleanTofVetoNsigma);

        const bool isCleanPi = isTofCleanPi && isTpcCleanPi;
        const bool isCleanKa = isTofCleanKa && isTpcCleanKa;
        const bool isCleanPr = isTofCleanPr && isTpcCleanPr;

        if (cfgSel.cfgGetAllSpeciesQA && cfgSel.cfgGetdEdx) {
          histos.fill(HIST("TpcdEdx_allSpecies_ptwise"), 1.0, multiplicity, track.pt(), track.tpcSignal(), tofNsigPi);
          histos.fill(HIST("TpcdEdx_allSpecies_ptwise"), 2.0, multiplicity, track.pt(), track.tpcSignal(), tofNsigKa);
          histos.fill(HIST("TpcdEdx_allSpecies_ptwise"), 3.0, multiplicity, track.pt(), track.tpcSignal(), tofNsigPr);

          histos.fill(HIST("ExpSigma_allSpecies_ptwise"), 1.0, multiplicity, track.pt(), expSigPi, tofNsigPi);
          histos.fill(HIST("ExpSigma_allSpecies_ptwise"), 2.0, multiplicity, track.pt(), expSigKa, tofNsigKa);
          histos.fill(HIST("ExpSigma_allSpecies_ptwise"), 3.0, multiplicity, track.pt(), expSigPr, tofNsigPr);

          if (isCleanPi) {
            histos.fill(HIST("TpcdEdx_allSpecies_afterCut"), 1.0, multiplicity, track.pt(), track.tpcSignal(), tofNsigPi);
          }
          if (isCleanKa) {
            histos.fill(HIST("TpcdEdx_allSpecies_afterCut"), 2.0, multiplicity, track.pt(), track.tpcSignal(), tofNsigKa);
          }
          if (isCleanPr) {
            histos.fill(HIST("TpcdEdx_allSpecies_afterCut"), 3.0, multiplicity, track.pt(), track.tpcSignal(), tofNsigPr);
          }
        }

        if (cfgSel.cfgGetNsigmaQA && (cfgSel.cfgPidMask == KPidMaskPion || cfgSel.cfgPidMask == KPidMaskKaon || cfgSel.cfgPidMask == KPidMaskProton)) {
          float tpcNsig = 0.0f, tofNsig = 0.0f, expSig = 0.0f;
          bool isTargetClean = false;

          if ((cfgSel.cfgPidMask & KPidMaskPion) != 0) {
            tpcNsig = tpcNsigPi;
            tofNsig = tofNsigPi;
            expSig = expSigPi;
            isTargetClean = isCleanPi;
          } else if ((cfgSel.cfgPidMask & KPidMaskKaon) != 0) {
            tpcNsig = tpcNsigKa;
            tofNsig = tofNsigKa;
            expSig = expSigKa;
            isTargetClean = isCleanKa;
          } else if ((cfgSel.cfgPidMask & KPidMaskProton) != 0) {
            tpcNsig = tpcNsigPr;
            tofNsig = tofNsigPr;
            expSig = expSigPr;
            isTargetClean = isCleanPr;
          }

          const float expSignal = track.tpcSignal() - (tpcNsig * expSig);

          histos.fill(HIST("TofTpcNsigma_before"), multiplicity, track.pt(), tpcNsig, tofNsig);
          if (cfgSel.cfgGetdEdx) {
            histos.fill(HIST("TpcdEdx_ptwise"), multiplicity, track.pt(), track.tpcSignal(), tofNsig);
            histos.fill(HIST("ExpTpcdEdx_ptwise"), multiplicity, track.pt(), expSignal, tofNsig);
            histos.fill(HIST("ExpSigma_ptwise"), multiplicity, track.pt(), expSig, tofNsig);
          }

          const bool useTrackCleanHere = !cfgSel.cfgFillAfterCutWithV0 || (cfgSel.cfgPidMask == KPidMaskKaon);
          if (useTrackCleanHere && isTargetClean) {
            histos.fill(HIST("TofTpcNsigma_after"), multiplicity, track.pt(), tpcNsig, tofNsig);
            if (cfgSel.cfgGetdEdx) {
              histos.fill(HIST("TpcdEdx_ptwise_afterCut"), multiplicity, track.pt(), track.tpcSignal(), tofNsig);
            }
          }
        }
      }

      if ((cfgSel.cfgGetV0CleanSample || cfgSel.cfgFillAfterCutWithV0) && cfgSel.cfgGetdEdx) {
        auto v0sThisCol = V0s.sliceByCached(aod::v0data::collisionId, rawCol.globalIndex(), this->cache);
        for (auto const& v0 : v0sThisCol) {
          auto const& posTrack = v0.template posTrack_as<TrksStandardTable>();
          auto const& negTrack = v0.template negTrack_as<TrksStandardTable>();
          fillCleanV0Daughters(rawCol, v0, posTrack, negTrack, multiplicity);
        }
      }
    }
  }

  PROCESS_SWITCH(LongrangecorrDerived, processTpcft0aSE, "same event TPC vs FT0A", false);
  PROCESS_SWITCH(LongrangecorrDerived, processTpcft0aME, "mixed event TPC vs FT0A", false);
  PROCESS_SWITCH(LongrangecorrDerived, processTpcft0cSE, "same event TPC vs FT0C", false);
  PROCESS_SWITCH(LongrangecorrDerived, processTpcft0cME, "mixed event TPC vs FT0C", false);
  PROCESS_SWITCH(LongrangecorrDerived, processV0ft0aSE, "same event V0 vs FT0A", false);
  PROCESS_SWITCH(LongrangecorrDerived, processV0ft0aME, "mixed event V0 vs FT0A", false);
  PROCESS_SWITCH(LongrangecorrDerived, processV0ft0cSE, "same event V0 vs FT0C", false);
  PROCESS_SWITCH(LongrangecorrDerived, processV0ft0cME, "mixed event V0 vs FT0C", false);
  PROCESS_SWITCH(LongrangecorrDerived, processV0tpcSE, "same event V0 vs TPC", false);
  PROCESS_SWITCH(LongrangecorrDerived, processV0tpcME, "mixed event V0 vs TPC", false);
  PROCESS_SWITCH(LongrangecorrDerived, processFt0aft0cSE, "same event FT0A vs FT0C", false);
  PROCESS_SWITCH(LongrangecorrDerived, processFt0aft0cME, "mixed event FT0A vs FT0C", false);
  PROCESS_SWITCH(LongrangecorrDerived, processUpcTpcft0aSE, "same UPC event TPC vs FT0A", false);
  PROCESS_SWITCH(LongrangecorrDerived, processUpcTpcft0aME, "mixed UPC event TPC vs FT0A", false);
  PROCESS_SWITCH(LongrangecorrDerived, processUpcTpcft0cSE, "same UPC event TPC vs FT0C", false);
  PROCESS_SWITCH(LongrangecorrDerived, processUpcTpcft0cME, "mixed UPC event TPC vs FT0C", false);
  PROCESS_SWITCH(LongrangecorrDerived, processUpcV0ft0aSE, "same UPC event V0 vs FT0A", false);
  PROCESS_SWITCH(LongrangecorrDerived, processUpcV0ft0aME, "mixed UPC event V0 vs FT0A", false);
  PROCESS_SWITCH(LongrangecorrDerived, processUpcV0ft0cSE, "same UPC event V0 vs FT0C", false);
  PROCESS_SWITCH(LongrangecorrDerived, processUpcV0ft0cME, "mixed UPC event V0 vs FT0C", false);
  PROCESS_SWITCH(LongrangecorrDerived, processMcTpcft0aSE, "same MC event TPC vs FT0A", false);
  PROCESS_SWITCH(LongrangecorrDerived, processMcTpcft0aME, "mixed MC event TPC vs FT0A", false);
  PROCESS_SWITCH(LongrangecorrDerived, processMcTpcft0cSE, "same MC event TPC vs FT0C", false);
  PROCESS_SWITCH(LongrangecorrDerived, processMcTpcft0cME, "mixed MC event TPC vs FT0C", false);
  PROCESS_SWITCH(LongrangecorrDerived, processMcFt0aft0cSE, "same MC event FT0A vs FT0C", false);
  PROCESS_SWITCH(LongrangecorrDerived, processMcFt0aft0cME, "mixed MC event FT0A vs FT0C", false);
  PROCESS_SWITCH(LongrangecorrDerived, processMcGenTpcft0aSE, "same MC gen event TPC vs FT0A", false);
  PROCESS_SWITCH(LongrangecorrDerived, processMcGenTpcft0aME, "mixed MC gen event TPC vs FT0A", false);
  PROCESS_SWITCH(LongrangecorrDerived, processMcGenTpcft0cSE, "same MC gen event TPC vs FT0C", false);
  PROCESS_SWITCH(LongrangecorrDerived, processMcGenTpcft0cME, "mixed MC gen event TPC vs FT0C", false);
  PROCESS_SWITCH(LongrangecorrDerived, processMcGenV0ft0aSE, "same MC gen event V0 vs FT0A", false);
  PROCESS_SWITCH(LongrangecorrDerived, processMcGenV0ft0aME, "mixed MC gen event V0 vs FT0A", false);
  PROCESS_SWITCH(LongrangecorrDerived, processMcGenV0ft0cSE, "same MC gen event V0 vs FT0C", false);
  PROCESS_SWITCH(LongrangecorrDerived, processMcGenV0ft0cME, "mixed MC gen event V0 vs FT0C", false);
  PROCESS_SWITCH(LongrangecorrDerived, processMcGenFt0aft0cSE, "same MC gen event FT0A vs FT0C", false);
  PROCESS_SWITCH(LongrangecorrDerived, processMcGenFt0aft0cME, "mixed MC gen event FT0A vs FT0C", false);
  PROCESS_SWITCH(LongrangecorrDerived, processTPCtrackEff, "process TPC track efficiency", false);
  PROCESS_SWITCH(LongrangecorrDerived, processPIDQA, "process PID QA on standard AODs", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& context)
{
  return WorkflowSpec{
    adaptAnalysisTask<LongrangecorrDerived>(context)};
}
