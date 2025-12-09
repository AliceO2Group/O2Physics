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

/// \file kaonIsospinFluctuations.cxx
/// \brief Kaon Isospin fluctuations
///
/// \author Rahul Verma (rahul.verma@iitb.ac.in) :: Sadhana Dash (sadhana@phy.iitb.ac.in)

#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/mcCentrality.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

#include <algorithm>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics; // for constants

#define ID_BIT_PI 0 // Identificationi bits for PID checks
#define ID_BIT_KA 1
#define ID_BIT_PR 2
#define ID_BIT_EL 3
#define ID_BIT_DE 4

#define BIT_IS_K0S 0
#define BIT_IS_LAMBDA 1
#define BIT_IS_ANTILAMBDA 2

// #define kPAIRBIT_ISLAMBDA

#define BIT_POS_DAU_HAS_SAME_COLL 0
#define BIT_NEG_DAU_HAS_SAME_COLL 1
#define BIT_BOTH_DAU_HAS_SAME_COLL 2

#define BITSET(mask, ithBit) ((mask) |= (1 << (ithBit)))  // avoid name bitset as std::bitset is already there
#define BITCHECK(mask, ithBit) ((mask) & (1 << (ithBit))) // bit check will return int value, not bool, use BITCHECK != 0 in Analysi

struct KaonIsospinFluctuations {
  // Hisogram registry:
  HistogramRegistry recoV0s{"recoV0s", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry recoEvent{"recoEvent", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry recoK0s{"recoK0s", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry recoTracks{"recoTracks", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry recoAnalysis{"recoAnalysis", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry genAnalysis{"genAnalysis", {}, OutputObjHandlingPolicy::AnalysisObject};

  // PDG data base
  Service<o2::framework::O2DatabasePDG> pdgDB;

  // Configurables
  // Event Selection
  Configurable<float> cutZvertex{"cutZvertex", 8.0f, "Accepted z-vertex range (cm)"};

  // Configurable parameters for V0 selection
  Configurable<float> v0settingDcaPosToPV{"v0settingDcaPosToPV", 0.06, "DCA Pos to PV"};
  Configurable<float> v0settingDcaNegToPV{"v0settingDcaNegToPV", 0.06, "DCA Neg to PV"};
  Configurable<float> v0settingDcaV0Dau{"v0settingDcaV0Dau", 1, "DCA V0 Daughters"};
  Configurable<double> v0settingCosPA{"v0settingCosPA", 0.98, "V0 CosPA"};
  Configurable<float> v0settingRadius{"v0settingRadius", 0.5, "v0radius"};

  // Configurable K0s
  struct : ConfigurableGroup {
    Configurable<float> cfgK0sMLow{"cfgK0sMLow", 0.48, "cfgK0sMLow"};
    Configurable<float> cfgK0sMHigh{"cfgK0sMHigh", 0.515, "cfgK0sMHigh"};
    Configurable<float> cfgK0sLowPt{"cfgK0sLowPt", 0.1, "cfgK0sLowPt"};
    Configurable<double> cfgK0sHighPt{"cfgK0sHighPt", 1.5, "cfgK0sHighPt"};
    Configurable<double> cfgK0sRapitidy{"cfgK0sRapitidy", 0.5, "cfgK0sRapitidy"};
    Configurable<double> cfgK0sARMcut{"cfgK0sARMcut", 0.2, "cfgK0sARMcut"};
  } k0sSelCut;

  // Histogram Configurables
  struct : ConfigurableGroup {
    Configurable<int> centBins{"centBins", 1020, "No of bins in centrality axis"};
    Configurable<double> centBinsxLow{"centBinsxLow", -1.0, "centBinsxLow"};
    Configurable<double> centBinsxUp{"centBinsxUp", 101.0, "centBinsxUp"};
    Configurable<int> centAxisType{"centAxisType", 0, "centAxisType"};
  } cfgCentAxis;
  // Track Configurables
  struct : ConfigurableGroup {
    Configurable<int> cfgTrkTpcNClsCrossedRows{"cfgTrkTpcNClsCrossedRows", 70, "cfgTrkTpcNClsCrossedRows"};
    Configurable<float> cfgTrkdcaXY{"cfgTrkdcaXY", 0.2, "cfgTrkdcaXY"};
    Configurable<bool> cfgDoVGselTrackCheck{"cfgDoVGselTrackCheck", false, "cfgDoVGselTrackCheck"};
    Configurable<float> cfgTrackEta{"cfgTrackEta", 0.8, "cfgTrackEta"};
    Configurable<float> cfgTrackPtLow{"cfgTrackPtLow", 0.15, "cfgTrackPtLow"};
    Configurable<float> cfgTrackPtHigh{"cfgTrackPtHigh", 2.0, "cfgTrackPtHigh"};
  } cfgTrackCuts;

  // Configurables for particle Identification
  Configurable<bool> cfgCheckVetoCut{"cfgCheckVetoCut", false, "cfgCheckVetoCut"};
  Configurable<bool> cfgDoElRejection{"cfgDoElRejection", true, "cfgDoElRejection"};
  Configurable<bool> cfgDoDeRejection{"cfgDoDeRejection", false, "cfgDoDeRejection"};
  Configurable<bool> cfgDoPdependentId{"cfgDoPdependentId", true, "cfgDoPdependentId"};
  Configurable<bool> cfgDoTpcInnerParamId{"cfgDoTpcInnerParamId", false, "cfgDoTpcInnerParamId"};

  Configurable<float> cfgPiThrPforTOF{"cfgPiThrPforTOF", 0.7, "cfgPiThrPforTOF"};
  Configurable<int> cfgPiIdCutTypeLowP{"cfgPiIdCutTypeLowP", 0, "cfgPiIdCutTypeLowP"};
  Configurable<float> cfgPiNSigmaTPCLowP{"cfgPiNSigmaTPCLowP", 3.0, "cfgPiNSigmaTPCLowP"};
  Configurable<float> cfgPiNSigmaTOFLowP{"cfgPiNSigmaTOFLowP", 3.0, "cfgPiNSigmaTOFLowP"};
  Configurable<float> cfgPiNSigmaRadLowP{"cfgPiNSigmaRadLowP", 9.0, "cfgPiNSigmaRadLowP"};
  Configurable<int> cfgPiIdCutTypeHighP{"cfgPiIdCutTypeHighP", 0, "cfgPiIdCutTypeHighP"};
  Configurable<float> cfgPiNSigmaTPCHighP{"cfgPiNSigmaTPCHighP", 3.0, "cfgPiNSigmaTPCHighP"};
  Configurable<float> cfgPiNSigmaTOFHighP{"cfgPiNSigmaTOFHighP", 3.0, "cfgPiNSigmaTOFHighP"};
  Configurable<float> cfgPiNSigmaRadHighP{"cfgPiNSigmaRadHighP", 9.0, "cfgPiNSigmaRadHighP"};

  Configurable<float> cfgKaThrPforTOF{"cfgKaThrPforTOF", 0.8, "cfgKaThrPforTOF"};
  Configurable<int> cfgKaIdCutTypeLowP{"cfgKaIdCutTypeLowP", 0, "cfgKaIdCutTypeLowP"};
  Configurable<float> cfgKaNSigmaTPCLowP{"cfgKaNSigmaTPCLowP", 3.0, "cfgKaNSigmaTPCLowP"};
  Configurable<float> cfgKaNSigmaTOFLowP{"cfgKaNSigmaTOFLowP", 3.0, "cfgKaNSigmaTOFLowP"};
  Configurable<float> cfgKaNSigmaRadLowP{"cfgKaNSigmaRadLowP", 9.0, "cfgKaNSigmaRadLowP"};
  Configurable<int> cfgKaIdCutTypeHighP{"cfgKaIdCutTypeHighP", 0, "cfgKaIdCutTypeHighP"};
  Configurable<float> cfgKaNSigmaTPCHighP{"cfgKaNSigmaTPCHighP", 3.0, "cfgKaNSigmaTPCHighP"};
  Configurable<float> cfgKaNSigmaTOFHighP{"cfgKaNSigmaTOFHighP", 3.0, "cfgKaNSigmaTOFHighP"};
  Configurable<float> cfgKaNSigmaRadHighP{"cfgKaNSigmaRadHighP", 9.0, "cfgKaNSigmaRadHighP"};

  Configurable<float> cfgPrThrPforTOF{"cfgPrThrPforTOF", 0.8, "cfgPrThrPforTOF"};
  Configurable<int> cfgPrIdCutTypeLowP{"cfgPrIdCutTypeLowP", 0, "cfgPrIdCutTypeLowP"};
  Configurable<float> cfgPrNSigmaTPCLowP{"cfgPrNSigmaTPCLowP", 3.0, "cfgPrNSigmaTPCLowP"};
  Configurable<float> cfgPrNSigmaTOFLowP{"cfgPrNSigmaTOFLowP", 3.0, "cfgPrNSigmaTOFLowP"};
  Configurable<float> cfgPrNSigmaRadLowP{"cfgPrNSigmaRadLowP", 9.0, "cfgPrNSigmaRadLowP"};
  Configurable<int> cfgPrIdCutTypeHighP{"cfgPrIdCutTypeHighP", 0, "cfgPrIdCutTypeHighP"};
  Configurable<float> cfgPrNSigmaTPCHighP{"cfgPrNSigmaTPCHighP", 3.0, "cfgPrNSigmaTPCHighP"};
  Configurable<float> cfgPrNSigmaTOFHighP{"cfgPrNSigmaTOFHighP", 3.0, "cfgPrNSigmaTOFHighP"};
  Configurable<float> cfgPrNSigmaRadHighP{"cfgPrNSigmaRadHighP", 9.0, "cfgPrNSigmaRadHighP"};

  // configurable for process functions to  reduce memory usage
  Configurable<bool> cfgFillV0TableFull{"cfgFillV0TableFull", true, "cfgFillV0TableFull"};
  Configurable<bool> cfgFillV0TablePostK0sCheck{"cfgFillV0TablePostK0sCheck", false, "cfgFillV0TablePostK0sCheck"};
  Configurable<bool> cfgFillV0TablePostMassCut{"cfgFillV0TablePostMassCut", false, "cfgFillV0TablePostMassCut"};
  Configurable<bool> cfgFillV0TablePostSelectionCut{"cfgFillV0TablePostSelectionCut", true, "cfgFillV0TablePostSelectionCut"};

  Configurable<bool> cfgFillRecoK0sPreSel{"cfgFillRecoK0sPreSel", false, "cfgFillRecoK0sPreSel"};
  Configurable<bool> cfgFillRecoK0sPostSel{"cfgFillRecoK0sPostSel", true, "cfgFillRecoK0sPostSel"};

  Configurable<bool> cfgFillRecoTrackPreSel{"cfgFillRecoTrackPreSel", false, "cfgFillRecoTrackPreSel"};
  Configurable<bool> cfgFillRecoTrackPostSel{"cfgFillRecoTrackPostSel", true, "cfgFillRecoTrackPostSel"};

  Configurable<bool> cfgFillPiQA{"cfgFillPiQA", true, "cfgFillPiQA"};
  Configurable<bool> cfgFillKaQA{"cfgFillKaQA", true, "cfgFillKaQA"};
  Configurable<bool> cfgFillPrQA{"cfgFillPrQA", true, "cfgFillPrQA"};
  Configurable<bool> cfgFillElQA{"cfgFillElQA", true, "cfgFillElQA"};
  Configurable<bool> cfgFillDeQA{"cfgFillDeQA", true, "cfgFillDeQA"};

  Configurable<bool> cfgFillSparseFullK0sPiKa{"cfgFillSparseFullK0sPiKa", true, "cfgFillSparseFullK0sPiKa"};
  Configurable<bool> cfgFillSparseFullK0sPrDe{"cfgFillSparseFullK0sPrDe", true, "cfgFillSparseFullK0sPrDe"};
  Configurable<bool> cfgFillSparseFullK0sKaEl{"cfgFillSparseFullK0sKaEl", false, "cfgFillSparseFullK0sKaEl"};
  Configurable<bool> cfgFillSparseFullPiKaPr{"cfgFillSparseFullPiKaPr", false, "cfgFillSparseFullPiKaPr"};
  Configurable<bool> cfgFillSparseFullPiElDe{"cfgFillSparseFullPiElDe", false, "cfgFillSparseFullPiElDe"};
  Configurable<bool> cfgFillSparseFullKaPrDe{"cfgFillSparseFullKaPrDe", false, "cfgFillSparseFullKaPrDe"};
  Configurable<bool> cfgFillSparseFullPrElDe{"cfgFillSparseFullPrElDe", false, "cfgFillSparseFullPrElDe"};
  Configurable<bool> cfgFillSparsenewDynmK0sKa{"cfgFillSparsenewDynmK0sKa", true, "cfgFillSparsenewDynmK0sKa"};
  Configurable<bool> cfgFillSparsenewDynmKpKm{"cfgFillSparsenewDynmKpKm", true, "cfgFillSparsenewDynmKpKm"};

  Configurable<bool> cfgVtxZCheck{"cfgVtxZCheck", 0, "cfgVtxZCheck"};
  Configurable<bool> cfgCountFinalParticles{"cfgCountFinalParticles", 1, "cfgCountFinalParticles"};
  Configurable<bool> cfgCountNonFinalParticles{"cfgCountNonFinalParticles", 0, "cfgCountNonFinalParticles"};
  Configurable<bool> cfgCountPhysicalPrimAndFinalParticles{"cfgCountPhysicalPrimAndFinalParticles", 0, "cfgCountPhysicalPrimAndFinalParticles"};
  Configurable<bool> doFWDPtDependentCheck{"doFWDPtDependentCheck", 1, "doFWDPtDependentCheck"};
  Configurable<float> cfgFWDPtCut{"cfgFWDPtCut", 0.2, "cfgFWDPtCut"};
  Configurable<std::vector<int>> cfgFinalParticleIdList{"cfgFinalParticleIdList", {11, -11, 13, -13, 15, -15, 211, -211, 321, -321, 2212, -2212}, "cfgFinalParticleIdList"};
  Configurable<std::vector<int>> cfgNonFinalParticleIdList{"cfgNonFinalParticleIdList", {11, -11, 13, -13, 15, -15, 211, -211, 321, -321, 2212, -2212}, "cfgNonFinalParticleIdList"};

  void init(InitContext const&)
  {
    // Axes
    const AxisSpec axisK0sMass = {200, 0.40f, 0.60f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    const AxisSpec axisLambdaMass = {200, 1.f, 1.2f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};

    const AxisSpec axisVertexZ = {30, -15., 15., "vrtx_{Z} [cm]"};
    AxisSpec axisCent = {cfgCentAxis.centBins, cfgCentAxis.centBinsxLow, cfgCentAxis.centBinsxUp, "centFT0C(percentile)"};
    if (cfgCentAxis.centAxisType == 1) {
      axisCent = {cfgCentAxis.centBins, cfgCentAxis.centBinsxLow, cfgCentAxis.centBinsxUp, "centFT0M(percentile)"};
    }
    if (cfgCentAxis.centAxisType == 2) {
      axisCent = {cfgCentAxis.centBins, cfgCentAxis.centBinsxLow, cfgCentAxis.centBinsxUp, "multFT0M"};
    }
    if (cfgCentAxis.centAxisType == 3) {
      axisCent = {cfgCentAxis.centBins, cfgCentAxis.centBinsxLow, cfgCentAxis.centBinsxUp, "multFT0C"};
    }

    const AxisSpec axisP = {200, 0.0f, 10.0f, "#it{p} (GeV/#it{c})"};
    const AxisSpec axisPt = {200, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisTPCInnerParam = {200, 0.0f, 10.0f, "#it{p}_{tpcInnerParam} (GeV/#it{c})"};
    const AxisSpec axisTOFExpMom = {200, 0.0f, 10.0f, "#it{p}_{tofExpMom} (GeV/#it{c})"};

    const AxisSpec axisEta = {100, -5, 5, "#eta"};
    const AxisSpec axisPhi = {90, -1, 8, "#phi (radians)"};
    const AxisSpec axisRapidity = {200, -5, 5, "Rapidity (y)"};
    const AxisSpec axisDcaXY = {100, -5, 5, "dcaXY"};
    const AxisSpec axisDcaZ = {100, -5, 5, "dcaZ"};
    const AxisSpec axisDcaXYwide = {2000, -100, 100, "dcaXY"};
    const AxisSpec axisDcaZwide = {2000, -100, 100, "dcaZ"};
    const AxisSpec axisSign = {10, -5, 5, "track.sign"};

    const AxisSpec axisTPCSignal = {100, -1, 1000, "tpcSignal"};
    const AxisSpec axisTOFBeta = {40, -2.0, 2.0, "tofBeta"};

    const AxisSpec axisTPCSignalFine = {10010, -1, 1000, "tpcSignal"};
    const AxisSpec axisTOFBetaFine = {10010, -1, 1000, "tpcSignal"};

    const AxisSpec axisTPCNSigma = {200, -10.0, 10.0, "n#sigma_{TPC}"};
    const AxisSpec axisTOFNSigma = {200, -10.0, 10.0, "n#sigma_{TOF}"};

    const AxisSpec axisTPCNSigmaPi = {200, -10.0, 10.0, "n#sigma_{TPC}^{Pi}"};
    const AxisSpec axisTOFNSigmaPi = {200, -10.0, 10.0, "n#sigma_{TOF}^{Pi}"};

    const AxisSpec axisTPCNClsCrossedRows = {200, -1.5, 198.5, "tpcNClsCrossedRows"};
    const AxisSpec axisIsPVContributor = {4, -1, 3, "isPVContributor"};
    const AxisSpec axisIsGlobalTrack = {4, -1, 3, "isGobalTrack"};
    const AxisSpec axisIsK0sDau = {4, -1, 3, "isK0sDau"};

    const AxisSpec axisDcapostopv = {100000, -50, 50, "dcapostopv"};
    const AxisSpec axisDcanegtopv = {100000, -50, 50, "dcanegtopv"};
    const AxisSpec axisDcaV0daughters = {2000, -10.0, 10.0, "dcaV0daughters"};
    const AxisSpec axisV0cosPA = {3000, -1.5, 1.5, "v0cosPA"};
    const AxisSpec axisV0radius = {100000, -50, 50, "v0radius"};

    const AxisSpec axisParticleCount1 = {60, -10, 50, "particleCount"};
    const AxisSpec axisParticleCount2 = {260, -10, 250, "particleCount"};
    const AxisSpec axisParticleCount3 = {1060, -10, 1050, "particleCount"};

    const AxisSpec axisArmenterosAlpha = {100, -1.0, 1.0, "ArmenterosAlpha"};
    const AxisSpec axisArmenterosQt = {150, 0, 0.3, "ArmenterosQt"};

    HistogramConfigSpec histPDcaXY({HistType::kTH2F, {axisP, axisDcaXY}});
    HistogramConfigSpec histPtDcaXY({HistType::kTH2F, {axisPt, axisDcaXY}});
    HistogramConfigSpec histTpcInnerParamDcaXY({HistType::kTH2F, {axisTPCInnerParam, axisDcaXY}});
    HistogramConfigSpec histTofExpMomDcaXY({HistType::kTH2F, {axisTOFExpMom, axisDcaXY}});

    HistogramConfigSpec histPDcaZ({HistType::kTH2F, {axisP, axisDcaZ}});
    HistogramConfigSpec histPtDcaZ({HistType::kTH2F, {axisPt, axisDcaZ}});
    HistogramConfigSpec histTpcInnerParamDcaZ({HistType::kTH2F, {axisTPCInnerParam, axisDcaZ}});
    HistogramConfigSpec histTofExpMomDcaZ({HistType::kTH2F, {axisTOFExpMom, axisDcaZ}});

    HistogramConfigSpec histPPt({HistType::kTH2F, {axisP, axisPt}});
    HistogramConfigSpec histPTpcInnerParam({HistType::kTH2F, {axisP, axisTPCInnerParam}});
    HistogramConfigSpec histPTofExpMom({HistType::kTH2F, {axisP, axisTOFExpMom}});

    HistogramConfigSpec histPTpcSignal({HistType::kTH2F, {axisP, axisTPCSignal}});
    HistogramConfigSpec histTpcInnerParamTpcSignal({HistType::kTH2F, {axisTPCInnerParam, axisTPCSignal}});
    HistogramConfigSpec histTofExpMomTpcSignal({HistType::kTH2F, {axisTOFExpMom, axisTPCSignal}});

    HistogramConfigSpec histPBeta({HistType::kTH2F, {axisP, axisTOFBeta}});
    HistogramConfigSpec histTpcInnerParamBeta({HistType::kTH2F, {axisTPCInnerParam, axisTOFBeta}});
    HistogramConfigSpec histTofExpMomBeta({HistType::kTH2F, {axisTOFExpMom, axisTOFBeta}});

    HistogramConfigSpec histPTpcNSigma({HistType::kTH2F, {axisP, axisTPCNSigma}});
    HistogramConfigSpec histPtTpcNSigma({HistType::kTH2F, {axisPt, axisTPCNSigma}});
    HistogramConfigSpec histTpcInnerParamTpcNSigma({HistType::kTH2F, {axisTPCInnerParam, axisTPCNSigma}});
    HistogramConfigSpec histTofExpMomTpcNSigma({HistType::kTH2F, {axisTOFExpMom, axisTPCNSigma}});
    HistogramConfigSpec histPTofNSigma({HistType::kTH2F, {axisP, axisTOFNSigma}});
    HistogramConfigSpec histPtTofNSigma({HistType::kTH2F, {axisPt, axisTOFNSigma}});
    HistogramConfigSpec histTpcInnerParamTofNSigma({HistType::kTH2F, {axisTPCInnerParam, axisTOFNSigma}});
    HistogramConfigSpec histTofExpMomTofNSigma({HistType::kTH2F, {axisTOFExpMom, axisTOFNSigma}});
    HistogramConfigSpec histTpcNSigmaTofNSigma({HistType::kTH2F, {axisTPCNSigma, axisTOFNSigma}});

    recoV0s.add("v0Table/Full/h01_K0s_Mass", "K0s_Mass", {HistType::kTH1F, {axisK0sMass}});
    recoV0s.add("v0Table/Full/h02_Lambda_Mass", "Lambda_Mass", {HistType::kTH1F, {axisLambdaMass}});
    recoV0s.add("v0Table/Full/h03_AntiLambda_Mass", "AntiLambda_Mass", {HistType::kTH1F, {axisLambdaMass}});
    recoV0s.add("v0Table/Full/h04_v0DaughterCollisionIndexTag", "hV0s_K0s_v0DaughterCollisionIndexTag", {HistType::kTH1D, {{22, -1.0, 10.0}}});
    recoV0s.add("v0Table/Full/h05_V0Tag", "V0Tag", {HistType::kTH1F, {{12, -2, 10}}}); // 001 = Kaon, 010 = Lambda, 100 = AnitLambda

    // Topological Cuts
    recoV0s.add("v0Table/Full/h06_dcapostopv", "dcapostopv", kTH1F, {axisDcapostopv});
    recoV0s.add("v0Table/Full/h07_dcanegtopv", "dcanegtopv", kTH1F, {axisDcanegtopv});
    recoV0s.add("v0Table/Full/h08_dcaV0daughters", "dcaV0daughters", kTH1F, {axisDcaV0daughters});
    recoV0s.add("v0Table/Full/h09_v0cosPA", "v0cosPA", kTH1F, {axisV0cosPA});
    recoV0s.add("v0Table/Full/h10_v0radius", "v0radius", kTH1F, {axisV0radius});

    // K0s-FullInformation
    recoV0s.add("v0Table/Full/h11_mass", "mass", kTH1F, {axisK0sMass});
    recoV0s.add("v0Table/Full/h12_p", "p", kTH1F, {axisP});
    recoV0s.add("v0Table/Full/h13_pt", "pt", kTH1F, {axisPt});
    recoV0s.add("v0Table/Full/h14_eta", "eta", kTH1F, {axisEta});
    recoV0s.add("v0Table/Full/h15_phi", "phi", kTH1F, {axisPhi});
    recoV0s.add("v0Table/Full/h16_rapidity", "rapidity", kTH1F, {axisRapidity});
    recoV0s.add("v0Table/Full/h17_alpha", "alpha", kTH1F, {axisArmenterosAlpha});
    recoV0s.add("v0Table/Full/h18_qtarm", "qtarm", kTH1F, {axisArmenterosQt});
    recoV0s.add("v0Table/Full/h19_alpha_qtarm", "alpha_qtarm", kTH2F, {axisArmenterosAlpha, axisArmenterosQt});
    recoV0s.add("v0Table/Full/h20_pt_eta", "pt_eta", kTH2F, {axisPt, axisEta});

    // K0s-Daughter Info
    recoV0s.add("v0Table/Full/Pi/tpcId/h01_p", "p", kTH1F, {axisP});
    recoV0s.add("v0Table/Full/Pi/tpcId/h02_pt", "pt", kTH1F, {axisPt});
    recoV0s.add("v0Table/Full/Pi/tpcId/h03_tpcInnerParam", "tpcInnerParam", kTH1F, {axisTPCInnerParam});
    recoV0s.add("v0Table/Full/Pi/tpcId/h04_tofExpMom", "tofExpMom", kTH1F, {axisTOFExpMom});
    recoV0s.add("v0Table/Full/Pi/tpcId/h05_eta", "eta", kTH1F, {axisEta});
    recoV0s.add("v0Table/Full/Pi/tpcId/h06_phi", "phi", kTH1F, {axisPhi});
    recoV0s.add("v0Table/Full/Pi/tpcId/h07_rapidity", "rapidity", kTH1F, {axisRapidity});
    recoV0s.add("v0Table/Full/Pi/tpcId/h08_isPVContributor", "isPVContributor", kTH1F, {axisIsPVContributor});
    recoV0s.add("v0Table/Full/Pi/tpcId/h09_isGlobalTrack", "isGlobalTrack", kTH1F, {axisIsGlobalTrack});
    recoV0s.add("v0Table/Full/Pi/tpcId/h10_dcaXY", "dcaXY", kTH1F, {axisDcaXY});
    recoV0s.add("v0Table/Full/Pi/tpcId/h11_dcaZ", "dcaZ", kTH1F, {axisDcaZ});

    recoV0s.add("v0Table/Full/Pi/tpcId/h12_p_dcaXY", "p_dcaXY", kTH2F, {axisP, axisDcaXY});
    recoV0s.add("v0Table/Full/Pi/tpcId/h13_p_dcaZ", "p_dcaZ", kTH2F, {axisP, axisDcaZ});
    recoV0s.add("v0Table/Full/Pi/tpcId/h14_pt_dcaXY", "pt_dcaXY", kTH2F, {axisP, axisDcaXY});
    recoV0s.add("v0Table/Full/Pi/tpcId/h15_pt_dcaZ", "pt_dcaZ", kTH2F, {axisP, axisDcaZ});
    recoV0s.add("v0Table/Full/Pi/tpcId/h16_dcaXYwide", "dcaXYwide", kTH1F, {axisDcaXYwide});
    recoV0s.add("v0Table/Full/Pi/tpcId/h17_dcaZwide", "dcaZwide", kTH1F, {axisDcaZwide});
    recoV0s.add("v0Table/Full/Pi/tpcId/h20_pt_eta", "pt_eta", kTH2F, {axisPt, axisEta});
    // K0s-Daughter identification
    // momemtum
    recoV0s.add("v0Table/Full/Pi/tpcId/h20_p_pt", "p_pt", histPPt);
    recoV0s.add("v0Table/Full/Pi/tpcId/h21_p_tpcInnerParam", "p_tpcInnerParam", histPTpcInnerParam);
    recoV0s.add("v0Table/Full/Pi/tpcId/h22_p_tofExpMom", "p_tofExpMom", histPTofExpMom);
    // tpcSignal
    recoV0s.add("v0Table/Full/Pi/tpcId/h23_p_tpcSignal", "p_tpcSignal", histPTpcSignal);
    recoV0s.add("v0Table/Full/Pi/tpcId/h24_tpcInnerParam_tpcSignal", "tpcInnerParam_tpcSignal", histTpcInnerParamTpcSignal);
    recoV0s.add("v0Table/Full/Pi/tpcId/h25_tofExpMom_tpcSignal", "tofExpMom_tpcSignal", histTofExpMomTpcSignal);
    // tofBeta
    recoV0s.add("v0Table/Full/Pi/tpcId/h26_p_beta", "p_beta", histPBeta);
    recoV0s.add("v0Table/Full/Pi/tpcId/h27_tpcInnerParam_beta", "tpcInnerParam_beta", histTpcInnerParamBeta);
    recoV0s.add("v0Table/Full/Pi/tpcId/h28_tofExpMom_beta", "tofExpMom_beta", histTofExpMomBeta);
    // Look at Pion
    recoV0s.add("v0Table/Full/Pi/tpcId/h29_p_tpcNSigma", "p_tpcNSigma", histPTpcNSigma);
    recoV0s.add("v0Table/Full/Pi/tpcId/h30_pt_tpcNSigma", "pt_tpcNSigma", histPtTpcNSigma);
    recoV0s.add("v0Table/Full/Pi/tpcId/h31_tpcInnerParam_tpcNSigma", "tpcInnerParam_tpcNSigma", histTpcInnerParamTpcNSigma);
    recoV0s.add("v0Table/Full/Pi/tpcId/h32_tofExpMom_tpcNSigma", "tofExpMom_tpcNSigma", histTofExpMomTpcNSigma);
    recoV0s.add("v0Table/Full/Pi/tpcId/h33_p_tofNSigma", "p_tofNSigma", histPTofNSigma);
    recoV0s.add("v0Table/Full/Pi/tpcId/h34_pt_tofNSigma", "pt_tofNSigma", histPtTofNSigma);
    recoV0s.add("v0Table/Full/Pi/tpcId/h35_tpcInnerParam_tofNSigma", "tpcInnerParam_tofNSigma", histTpcInnerParamTofNSigma);
    recoV0s.add("v0Table/Full/Pi/tpcId/h36_tofExpMom_tofNSigma", "tofExpMom_tofNSigma", histTofExpMomTofNSigma);
    recoV0s.add("v0Table/Full/Pi/tpcId/h37_tpcNSigma_tofNSigma", "tpcNSigma_tofNSigma", histTpcNSigmaTofNSigma);

    recoV0s.addClone("v0Table/Full/Pi/tpcId/", "v0Table/Full/Pi/tpctofId/"); // for identification using tof+tpc
    recoV0s.addClone("v0Table/Full/Pi/tpcId/", "v0Table/Full/Pi/NoId/");     // for unidentified case // to observe and debug

    if (cfgFillV0TablePostK0sCheck) {
      recoV0s.addClone("v0Table/Full/", "v0Table/postK0sCheck/");
    }
    if (cfgFillV0TablePostMassCut) {
      recoV0s.addClone("v0Table/Full/", "v0Table/postMassCut/");
    }
    if (cfgFillV0TablePostSelectionCut) {
      recoV0s.addClone("v0Table/Full/", "v0Table/postSelectionCut/");
    }

    recoV0s.add("v0Table/postSelectionCut/hTrueV0TagCount", "hTrueV0TagCount", {HistType::kTH1F, {{12, -2, 10}}}); // 001 = Kaon, 010 = Lambda, 100 = AnitLambda
    recoV0s.add("v0Table/postSelectionCut/nCommonPionOfDifferentK0s", "nCommonPionOfDifferentK0s", {HistType::kTH1D, {{44, -2, 20}}});

    // Event Selection
    recoEvent.add("recoEvent/ProcessType", "ProcessType", {HistType::kTH1D, {{20, -1, 9}}});
    recoEvent.add("recoEvent/h01_CollisionCount", "CollisionCount", {HistType::kTH1D, {{1, 0, 1}}});
    recoEvent.add("recoEvent/h02_VertexXRec", "VertexXRec", {HistType::kTH1D, {{1000, -0.2, 0.2}}});
    recoEvent.add("recoEvent/h03_VertexYRec", "VertexYRec", {HistType::kTH1D, {{1000, -0.2, 0.2}}});
    recoEvent.add("recoEvent/h04_VertexZRec", "VertexZRec", {HistType::kTH1F, {axisVertexZ}});
    recoEvent.add("recoEvent/h05_Centrality", "Centrality", {HistType::kTH1F, {axisCent}});
    recoEvent.add("recoEvent/h06_V0Size", "V0Size", {HistType::kTH1F, {{60, -10, 50}}});
    recoEvent.add("recoEvent/h07_TracksSize", "TracksSize", {HistType::kTH1F, {axisParticleCount2}});
    recoEvent.add("recoEvent/h08_nTrack", "nTrack", {HistType::kTH1F, {axisParticleCount2}});
    recoEvent.add("recoEvent/h09_nK0s", "nK0s", {HistType::kTH1F, {axisParticleCount1}});
    recoEvent.add("recoEvent/h10_nPiPlus", "nPiPlus", {HistType::kTH1F, {axisParticleCount2}});
    recoEvent.add("recoEvent/h11_nPiMinus", "nPiMinus", {HistType::kTH1F, {axisParticleCount2}});
    recoEvent.add("recoEvent/h12_nKaPlus", "nKaPlus", {HistType::kTH1F, {axisParticleCount1}});
    recoEvent.add("recoEvent/h13_nKaMinus", "nKaMinus", {HistType::kTH1F, {axisParticleCount1}});
    recoEvent.add("recoEvent/h14_nProton", "nProton", {HistType::kTH1F, {axisParticleCount1}});
    recoEvent.add("recoEvent/h15_nPBar", "nPBar", {HistType::kTH1F, {axisParticleCount1}});
    recoEvent.add("recoEvent/h16_nElPlus", "nElPlus", {HistType::kTH1F, {axisParticleCount1}});
    recoEvent.add("recoEvent/h17_nElMinus", "nElMinus", {HistType::kTH1F, {axisParticleCount1}});
    recoEvent.add("recoEvent/h18_nDePlus", "nDePlus", {HistType::kTH1F, {axisParticleCount1}});
    recoEvent.add("recoEvent/h19_nDeMinus", "nDeMinus", {HistType::kTH1F, {axisParticleCount1}});

    //
    // K0s reconstruction
    recoK0s.add("recoK0s/PreSel/h01_K0s_Mass", "K0s_Mass", {HistType::kTH1F, {axisK0sMass}});
    recoK0s.add("recoK0s/PreSel/h02_Lambda_Mass", "Lambda_Mass", {HistType::kTH1F, {axisLambdaMass}});
    recoK0s.add("recoK0s/PreSel/h03_AntiLambda_Mass", "AntiLambda_Mass", {HistType::kTH1F, {axisLambdaMass}});
    recoK0s.add("recoK0s/PreSel/h04_v0DaughterCollisionIndexTag", "hV0s_K0s_v0DaughterCollisionIndexTag", {HistType::kTH1D, {{22, -1.0, 10.0}}});
    recoK0s.add("recoK0s/PreSel/h05_V0Tag", "V0Tag", {HistType::kTH1F, {{12, -2, 10}}}); // 001 = Kaon, 010 = Lambda, 100 = AnitLambda

    // Topological Cuts
    recoK0s.add("recoK0s/PreSel/h06_dcapostopv", "dcapostopv", kTH1F, {axisDcapostopv});
    recoK0s.add("recoK0s/PreSel/h07_dcanegtopv", "dcanegtopv", kTH1F, {axisDcanegtopv});
    recoK0s.add("recoK0s/PreSel/h08_dcaV0daughters", "dcaV0daughters", kTH1F, {axisDcaV0daughters});
    recoK0s.add("recoK0s/PreSel/h09_v0cosPA", "v0cosPA", kTH1F, {axisV0cosPA});
    recoK0s.add("recoK0s/PreSel/h10_v0radius", "v0radius", kTH1F, {axisV0radius});

    // K0s-FullInformation
    recoK0s.add("recoK0s/PreSel/h11_mass", "mass", kTH1F, {axisK0sMass});
    recoK0s.add("recoK0s/PreSel/h12_p", "p", kTH1F, {axisP});
    recoK0s.add("recoK0s/PreSel/h13_pt", "pt", kTH1F, {axisPt});
    recoK0s.add("recoK0s/PreSel/h14_eta", "eta", kTH1F, {axisEta});
    recoK0s.add("recoK0s/PreSel/h15_phi", "phi", kTH1F, {axisPhi});
    recoK0s.add("recoK0s/PreSel/h16_rapidity", "rapidity", kTH1F, {axisRapidity});
    recoK0s.add("recoK0s/PreSel/h17_alpha", "alpha", kTH1F, {axisArmenterosAlpha});
    recoK0s.add("recoK0s/PreSel/h18_qtarm", "qtarm", kTH1F, {axisArmenterosQt});
    recoK0s.add("recoK0s/PreSel/h19_alpha_qtarm", "alpha_qtarm", kTH2F, {axisArmenterosAlpha, axisArmenterosQt});
    recoK0s.add("recoK0s/PreSel/h20_pt_eta", "pt_eta", kTH2F, {axisPt, axisEta});

    // K0s-Daughter Info
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h01_p", "p", kTH1F, {axisP});
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h02_pt", "pt", kTH1F, {axisPt});
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h03_tpcInnerParam", "tpcInnerParam", kTH1F, {axisTPCInnerParam});
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h04_tofExpMom", "tofExpMom", kTH1F, {axisTOFExpMom});
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h05_eta", "eta", kTH1F, {axisEta});
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h06_phi", "phi", kTH1F, {axisPhi});
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h07_rapidity", "rapidity", kTH1F, {axisRapidity});
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h08_isPVContributor", "isPVContributor", kTH1F, {axisIsPVContributor});
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h09_isGlobalTrack", "isGlobalTrack", kTH1F, {axisIsGlobalTrack});
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h10_dcaXY", "dcaXY", kTH1F, {axisDcaXY});
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h11_dcaZ", "dcaZ", kTH1F, {axisDcaZ});
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h12_p_dcaXY", "p_dcaXY", kTH2F, {axisP, axisDcaXY});
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h13_p_dcaZ", "p_dcaZ", kTH2F, {axisP, axisDcaZ});
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h14_pt_dcaXY", "pt_dcaXY", kTH2F, {axisP, axisDcaXY});
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h15_pt_dcaZ", "pt_dcaZ", kTH2F, {axisP, axisDcaZ});
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h16_dcaXYwide", "dcaXYwide", kTH1F, {axisDcaXYwide});
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h17_dcaZwide", "dcaZwide", kTH1F, {axisDcaZwide});

    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h20_pt_eta", "pt_eta", kTH2F, {axisPt, axisEta});
    // K0s-Daughter identification
    // momemtum
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h20_p_pt", "p_pt", histPPt);
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h21_p_tpcInnerParam", "p_tpcInnerParam", histPTpcInnerParam);
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h22_p_tofExpMom", "p_tofExpMom", histPTofExpMom);
    // tpcSignal
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h23_p_tpcSignal", "p_tpcSignal", histPTpcSignal);
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h24_tpcInnerParam_tpcSignal", "tpcInnerParam_tpcSignal", histTpcInnerParamTpcSignal);
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h25_tofExpMom_tpcSignal", "tofExpMom_tpcSignal", histTofExpMomTpcSignal);
    // tofBeta
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h26_p_beta", "p_beta", histPBeta);
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h27_tpcInnerParam_beta", "tpcInnerParam_beta", histTpcInnerParamBeta);
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h28_tofExpMom_beta", "tofExpMom_beta", histTofExpMomBeta);
    // Look at Pion
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h29_p_tpcNSigma", "p_tpcNSigma", histPTpcNSigma);
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h30_pt_tpcNSigma", "pt_tpcNSigma", histPtTpcNSigma);
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h31_tpcInnerParam_tpcNSigma", "tpcInnerParam_tpcNSigma", histTpcInnerParamTpcNSigma);
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h32_tofExpMom_tpcNSigma", "tofExpMom_tpcNSigma", histTofExpMomTpcNSigma);
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h33_p_tofNSigma", "p_tofNSigma", histPTofNSigma);
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h34_pt_tofNSigma", "pt_tofNSigma", histPtTofNSigma);
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h35_tpcInnerParam_tofNSigma", "tpcInnerParam_tofNSigma", histTpcInnerParamTofNSigma);
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h36_tofExpMom_tofNSigma", "tofExpMom_tofNSigma", histTofExpMomTofNSigma);
    recoK0s.add("recoK0s/PreSel/Pi/tpcId/h37_tpcNSigma_tofNSigma", "tpcNSigma_tofNSigma", histTpcNSigmaTofNSigma);

    recoK0s.addClone("recoK0s/PreSel/Pi/tpcId/", "recoK0s/PreSel/Pi/tpctofId/"); // for identification using tof+tpc
    recoK0s.addClone("recoK0s/PreSel/Pi/tpcId/", "recoK0s/PreSel/Pi/NoId/");     // for unidentified case // to observe and debug

    recoK0s.addClone("recoK0s/PreSel/", "recoK0s/PostSel/"); // for unidentified case // to observe and debug

    recoK0s.add("recoK0s/PostSel/mK0s_vs_cent", "mK0s_vs_cent", kTH2F, {axisCent, axisK0sMass});

    // Tracks reconstruction
    // FullTrack
    recoTracks.add("recoTracks/PreSel/h01_p", "p", {HistType::kTH1F, {axisP}});
    recoTracks.add("recoTracks/PreSel/h02_pt", "pt", {HistType::kTH1F, {axisPt}});
    recoTracks.add("recoTracks/PreSel/h03_tpcInnerParam", "tpcInnerParam", {HistType::kTH1F, {axisTPCInnerParam}});
    recoTracks.add("recoTracks/PreSel/h04_tofExpMom", "tofExpMom", {HistType::kTH1F, {axisTOFExpMom}});
    recoTracks.add("recoTracks/PreSel/h05_eta", "eta", {HistType::kTH1F, {axisEta}});
    recoTracks.add("recoTracks/PreSel/h06_phi", "phi", {HistType::kTH1F, {axisPhi}});
    recoTracks.add("recoTracks/PreSel/h07_dcaXY", "dcaXY", {HistType::kTH1F, {axisDcaXY}});
    recoTracks.add("recoTracks/PreSel/h08_dcaZ", "dcaZ", {HistType::kTH1F, {axisDcaZ}});
    recoTracks.add("recoTracks/PreSel/h09_sign", "sign", {HistType::kTH1D, {axisSign}});

    // DcaXY
    recoTracks.add("recoTracks/PreSel/h10_p_dcaXY", "p_dcaXY", histPDcaXY);
    recoTracks.add("recoTracks/PreSel/h11_pt_dcaXY", "pt_dcaXY", histPtDcaXY);
    recoTracks.add("recoTracks/PreSel/h12_tpcInnerParam_dcaXY", "tpcInnerParam_dcaXY", histTpcInnerParamDcaXY);
    recoTracks.add("recoTracks/PreSel/h13_tofExpMom_dcaXY", "tofExpMom_dcaXY", histTofExpMomDcaXY);

    // DcaZ
    recoTracks.add("recoTracks/PreSel/h14_p_dcaZ", "p_dcaZ", histPDcaZ);
    recoTracks.add("recoTracks/PreSel/h15_pt_dcaZ", "pt_dcaZ", histPtDcaZ);
    recoTracks.add("recoTracks/PreSel/h16_tpcInnerParam_dcaZ", "tpcInnerParam_dcaZ", histTpcInnerParamDcaZ);
    recoTracks.add("recoTracks/PreSel/h17_tofExpMom_dcaZ", "tofExpMom_dcaZ", histTofExpMomDcaZ);

    recoTracks.add("recoTracks/PreSel/h20_pt_eta", "pt_eta", kTH2F, {axisPt, axisEta});
    // momemtum
    recoTracks.add("recoTracks/PreSel/h20_p_pt", "p_pt", histPPt);
    recoTracks.add("recoTracks/PreSel/h21_p_tpcInnerParam", "p_tpcInnerParam", histPTpcInnerParam);
    recoTracks.add("recoTracks/PreSel/h22_p_tofExpMom", "p_tofExpMom", histPTofExpMom);

    // tpcSignal
    recoTracks.add("recoTracks/PreSel/h23_p_tpcSignal", "p_tpcSignal", histPTpcSignal);
    recoTracks.add("recoTracks/PreSel/h24_tpcInnerParam_tpcSignal", "tpcInnerParam_tpcSignal", histTpcInnerParamTpcSignal);
    recoTracks.add("recoTracks/PreSel/h25_tofExpMom_tpcSignal", "tofExpMom_tpcSignal", histTofExpMomTpcSignal);

    // tofBeta
    recoTracks.add("recoTracks/PreSel/h26_p_beta", "p_beta", histPBeta);
    recoTracks.add("recoTracks/PreSel/h27_tpcInnerParam_beta", "tpcInnerParam_beta", histTpcInnerParamBeta);
    recoTracks.add("recoTracks/PreSel/h28_tofExpMom_beta", "tofExpMom_beta", histTofExpMomBeta);

    // Look at Pion
    recoTracks.add("recoTracks/PreSel/Pi/NoId/h29_p_tpcNSigma", "p_tpcNSigma", histPTpcNSigma);
    recoTracks.add("recoTracks/PreSel/Pi/NoId/h30_pt_tpcNSigma", "pt_tpcNSigma", histPtTpcNSigma);
    recoTracks.add("recoTracks/PreSel/Pi/NoId/h31_tpcInnerParam_tpcNSigma", "tpcInnerParam_tpcNSigma", histTpcInnerParamTpcNSigma);
    recoTracks.add("recoTracks/PreSel/Pi/NoId/h32_tofExpMom_tpcNSigma", "tofExpMom_tpcNSigma", histTofExpMomTpcNSigma);
    recoTracks.add("recoTracks/PreSel/Pi/NoId/h33_p_tofNSigma", "p_tofNSigma", histPTofNSigma);
    recoTracks.add("recoTracks/PreSel/Pi/NoId/h34_pt_tofNSigma", "pt_tofNSigma", histPtTofNSigma);
    recoTracks.add("recoTracks/PreSel/Pi/NoId/h35_tpcInnerParam_tofNSigma", "tpcInnerParam_tofNSigma", histTpcInnerParamTofNSigma);
    recoTracks.add("recoTracks/PreSel/Pi/NoId/h36_tofExpMom_tofNSigma", "tofExpMom_tofNSigma", histTofExpMomTofNSigma);
    recoTracks.add("recoTracks/PreSel/Pi/NoId/h37_tpcNSigma_tofNSigma", "tpcNSigma_tofNSigma", histTpcNSigmaTofNSigma);
    // Pion

    recoTracks.addClone("recoTracks/PreSel/Pi/", "recoTracks/PreSel/Ka/"); // Kaon
    recoTracks.addClone("recoTracks/PreSel/Pi/", "recoTracks/PreSel/Pr/"); // Proton
    recoTracks.addClone("recoTracks/PreSel/Pi/", "recoTracks/PreSel/El/"); // Electron
    recoTracks.addClone("recoTracks/PreSel/Pi/", "recoTracks/PreSel/De/"); // Deuteron

    // Write Code for naming the axis for Identified Particles

    // Selection
    recoTracks.addClone("recoTracks/PreSel/", "recoTracks/PostSel/");
    //

    // Analysis
    recoAnalysis.add("recoAnalysis/Pi/tpcId/h20_pt_eta", "pt_eta", kTH2F, {axisPt, axisEta});
    // momemtum
    recoAnalysis.add("recoAnalysis/Pi/tpcId/h20_p_pt", "p_pt", histPPt);
    recoAnalysis.add("recoAnalysis/Pi/tpcId/h21_p_tpcInnerParam", "p_tpcInnerParam", histPTpcInnerParam);
    recoAnalysis.add("recoAnalysis/Pi/tpcId/h22_p_tofExpMom", "p_tofExpMom", histPTofExpMom);
    // tpcSignal
    recoAnalysis.add("recoAnalysis/Pi/tpcId/h23_p_tpcSignal", "p_tpcSignal", histPTpcSignal);
    recoAnalysis.add("recoAnalysis/Pi/tpcId/h24_tpcInnerParam_tpcSignal", "tpcInnerParam_tpcSignal", histTpcInnerParamTpcSignal);
    recoAnalysis.add("recoAnalysis/Pi/tpcId/h25_tofExpMom_tpcSignal", "tofExpMom_tpcSignal", histTofExpMomTpcSignal);
    // tofBeta
    recoAnalysis.add("recoAnalysis/Pi/tpcId/h26_p_beta", "p_beta", histPBeta);
    recoAnalysis.add("recoAnalysis/Pi/tpcId/h27_tpcInnerParam_beta", "tpcInnerParam_beta", histTpcInnerParamBeta);
    recoAnalysis.add("recoAnalysis/Pi/tpcId/h28_tofExpMom_beta", "tofExpMom_beta", histTofExpMomBeta);
    // Pion
    recoAnalysis.add("recoAnalysis/Pi/tpcId/h29_p_tpcNSigma", "p_tpcNSigma", histPTpcNSigma);
    recoAnalysis.add("recoAnalysis/Pi/tpcId/h30_pt_tpcNSigma", "pt_tpcNSigma", histPtTpcNSigma);
    recoAnalysis.add("recoAnalysis/Pi/tpcId/h31_tpcInnerParam_tpcNSigma", "tpcInnerParam_tpcNSigma", histTpcInnerParamTpcNSigma);
    recoAnalysis.add("recoAnalysis/Pi/tpcId/h32_tofExpMom_tpcNSigma", "tofExpMom_tpcNSigma", histTofExpMomTpcNSigma);
    recoAnalysis.add("recoAnalysis/Pi/tpcId/h33_p_tofNSigma", "p_tofNSigma", histPTofNSigma);
    recoAnalysis.add("recoAnalysis/Pi/tpcId/h34_pt_tofNSigma", "pt_tofNSigma", histPtTofNSigma);
    recoAnalysis.add("recoAnalysis/Pi/tpcId/h35_tpcInnerParam_tofNSigma", "tpcInnerParam_tofNSigma", histTpcInnerParamTofNSigma);
    recoAnalysis.add("recoAnalysis/Pi/tpcId/h36_tofExpMom_tofNSigma", "tofExpMom_tofNSigma", histTofExpMomTofNSigma);
    recoAnalysis.add("recoAnalysis/Pi/tpcId/h37_tpcNSigma_tofNSigma", "tpcNSigma_tofNSigma", histTpcNSigmaTofNSigma);

    recoAnalysis.addClone("recoAnalysis/Pi/tpcId/", "recoAnalysis/Pi/tpctofId/");
    recoAnalysis.addClone("recoAnalysis/Pi/tpcId/", "recoAnalysis/Pi/NoId/");
    recoAnalysis.addClone("recoAnalysis/Pi/", "recoAnalysis/Ka/"); // Kaon
    recoAnalysis.addClone("recoAnalysis/Pi/", "recoAnalysis/Pr/"); // Proton
    recoAnalysis.addClone("recoAnalysis/Pi/", "recoAnalysis/El/"); // Electron
    recoAnalysis.addClone("recoAnalysis/Pi/", "recoAnalysis/De/"); // Deuteron

    recoAnalysis.add("recoAnalysis/SelectedTrack_IdentificationTag", "SelectedTrack_IdentificationTag", kTH1D, {{34, -1.5, 32.5, "trackTAG"}});
    recoAnalysis.add("recoAnalysis/RejectedTrack_RejectionTag", "RejectedTrack_RejectionTag", kTH1D, {{16, -1.5, 6.5, "rejectionTAG"}});

    recoAnalysis.add("recoAnalysis/Sparse_Full_K0sPiKa", "Sparse_Full_K0sPiKa", kTHnSparseD, {axisCent, {2000, -1.5, 1998.5, "nTrack"}, {100, -1.5, 98.5, "nK0s"}, {100, -1.5, 98.5, "nRejectedPiPlus"}, {100, -1.5, 98.5, "nRejectedPiMinus"}, {500, -1.5, 498.5, "nPiPlus"}, {500, -1.5, 498.5, "nPiMinus"}, {500, -1.5, 498.5, "nKaPlus"}, {500, -1.5, 498.5, "nKaMinus"}});
    recoAnalysis.add("recoAnalysis/Sparse_Full_K0sPrDe", "Sparse_Full_K0sPrDe", kTHnSparseD, {axisCent, {2000, -1.5, 1998.5, "nTrack"}, {100, -1.5, 98.5, "nK0s"}, {100, -1.5, 98.5, "nRejectedPiPlus"}, {100, -1.5, 98.5, "nRejectedPiMinus"}, {500, -1.5, 498.5, "nProton"}, {500, -1.5, 498.5, "nPBar"}, {500, -1.5, 498.5, "nDePlus"}, {500, -1.5, 498.5, "nDeMinus"}});
    recoAnalysis.add("recoAnalysis/Sparse_Full_K0sKaEl", "Sparse_Full_K0sKaEl", kTHnSparseD, {axisCent, {2000, -1.5, 1998.5, "nTrack"}, {100, -1.5, 98.5, "nK0s"}, {100, -1.5, 98.5, "nRejectedPiPlus"}, {100, -1.5, 98.5, "nRejectedPiMinus"}, {500, -1.5, 498.5, "nKaPlus"}, {500, -1.5, 498.5, "nKaMinus"}, {500, -1.5, 498.5, "nElPlus"}, {500, -1.5, 498.5, "nElMinus"}});
    recoAnalysis.add("recoAnalysis/Sparse_Full_PiKaPr", "Sparse_Full_PiKaPr", kTHnSparseD, {axisCent, {2000, -1.5, 1998.5, "nTrack"}, {100, -1.5, 98.5, "nRejectedPiPlus"}, {100, -1.5, 98.5, "nRejectedPiMinus"}, {500, -1.5, 498.5, "nPiPlus"}, {500, -1.5, 498.5, "nPiMinus"}, {500, -1.5, 498.5, "nKaPlus"}, {500, -1.5, 498.5, "nKaMinus"}, {500, -1.5, 498.5, "nProton"}, {500, -1.5, 498.5, "nPBar"}});
    recoAnalysis.add("recoAnalysis/Sparse_Full_PiElDe", "Sparse_Full_PiElDe", kTHnSparseD, {axisCent, {2000, -1.5, 1998.5, "nTrack"}, {100, -1.5, 98.5, "nRejectedPiPlus"}, {100, -1.5, 98.5, "nRejectedPiMinus"}, {500, -1.5, 498.5, "nPiPlus"}, {500, -1.5, 498.5, "nPiMinus"}, {500, -1.5, 498.5, "nElPlus"}, {500, -1.5, 498.5, "nElMinus"}, {500, -1.5, 498.5, "nDePlus"}, {500, -1.5, 498.5, "nDeMinus"}});
    recoAnalysis.add("recoAnalysis/Sparse_Full_KaPrDe", "Sparse_Full_KaPrDe", kTHnSparseD, {axisCent, {2000, -1.5, 1998.5, "nTrack"}, {100, -1.5, 98.5, "nRejectedPiPlus"}, {100, -1.5, 98.5, "nRejectedPiMinus"}, {500, -1.5, 498.5, "nKaPlus"}, {500, -1.5, 498.5, "nKaMinus"}, {500, -1.5, 498.5, "nProton"}, {500, -1.5, 498.5, "nPBar"}, {500, -1.5, 498.5, "nDePlus"}, {500, -1.5, 498.5, "nDeMinus"}});
    recoAnalysis.add("recoAnalysis/Sparse_Full_PrElDe", "Sparse_Full_PrElDe", kTHnSparseD, {axisCent, {2000, -1.5, 1998.5, "nTrack"}, {100, -1.5, 98.5, "nRejectedPiPlus"}, {100, -1.5, 98.5, "nRejectedPiMinus"}, {500, -1.5, 498.5, "nProton"}, {500, -1.5, 498.5, "nPBar"}, {500, -1.5, 498.5, "nElPlus"}, {500, -1.5, 498.5, "nElMinus"}, {500, -1.5, 498.5, "nDePlus"}, {500, -1.5, 498.5, "nDeMinus"}});

    recoAnalysis.add("recoAnalysis/Sparse_newDynm_K0s_Ka", "Sparse_newDynm_K0s_Ka", kTHnSparseD, {axisCent, {2000, -1.5, 1998.5, "nTrack"}, {100, -1.5, 98.5, "nK0s"}, {500, -1.5, 498.5, "nKaon"}, {10000, -1.5, 9998.5, "(nK0s)^{2}"}, {250000, -1.5, 249998.5, "(nKaon)^{2}"}, {500, -1.5, 498.5, "(nK0s*nKaon)"}});
    recoAnalysis.add("recoAnalysis/Sparse_newDynm_Kp_Km", "Sparse_newDynm_Kp_Km", kTHnSparseD, {axisCent, {2000, -1.5, 1998.5, "nTrack"}, {500, -1.5, 498.5, "nKaPlus"}, {500, -1.5, 498.5, "nKaMinus"}, {250000, -1.5, 249998.5, "(nKaPlus)^{2}"}, {250000, -1.5, 249998.5, "(nKaMinus)^{2}"}, {250000, -1.5, 249998.5, "(nKaPlus*nKaMinus)"}});
    //

    genAnalysis.add("genAnalysis/K0s/h12_p", "p", kTH1F, {axisP});
    genAnalysis.add("genAnalysis/K0s/h13_pt", "pt", kTH1F, {axisPt});
    genAnalysis.add("genAnalysis/K0s/h14_eta", "eta", kTH1F, {axisEta});
    genAnalysis.add("genAnalysis/K0s/h15_phi", "phi", kTH1F, {axisPhi});
    genAnalysis.add("genAnalysis/K0s/h16_rapidity", "rapidity", kTH1F, {axisRapidity});
    genAnalysis.add("genAnalysis/K0s/h20_pt_eta", "pt_eta", kTH2F, {axisPt, axisEta});
    genAnalysis.addClone("genAnalysis/K0s/", "genAnalysis/Pi/");
    genAnalysis.addClone("genAnalysis/K0s/", "genAnalysis/Ka/");
    genAnalysis.addClone("genAnalysis/K0s/", "genAnalysis/Pr/");
    genAnalysis.addClone("genAnalysis/K0s/", "genAnalysis/El/");
    genAnalysis.addClone("genAnalysis/K0s/", "genAnalysis/De/");

    // Printing the Stored Registry information
    LOG(info) << "Printing Stored Registry Information";
    LOG(info) << " DEBUG :: 01- recoV0s.print()";
    recoV0s.print();
    LOG(info) << " DEBUG :: 02- recoEvent.print()";
    recoEvent.print();
    LOG(info) << " DEBUG :: 03- recoK0s.print()";
    recoK0s.print();
    LOG(info) << " DEBUG :: 04- recoTracks.print()";
    recoTracks.print();
    LOG(info) << " DEBUG :: 05- recoAnalysis.print()";
    recoAnalysis.print();
    LOG(info) << " DEBUG :: 06- genAnalysis.print()";
    genAnalysis.print();
  }

  enum RejectionTagEnum {
    kPassed = 0,
    kFailTpcNClsCrossedRows,
    kFailTrkdcaXY,
    kFailGlobalTrack,
    kFailVGSelCheck,
    kFailK0ShortDaughter,
    kFailPhiDaughter,
  };

  enum IdentificationType {
    kTPCidentified = 0,
    kTOFidentified,
    kTPCTOFidentified,
    kUnidentified
  };

  enum TpcTofCutType {
    kRectangularCut = 0,
    kCircularCut,
    kEllipsoidalCut
  };

  enum ProcessTypeEnum {
    doDataProcessing = 0,
    doRecoProcessing,
    doPurityProcessing,
    doGenProcessing,
    doSimProcessing
  };

  enum HistRegEnum {
    v0TableFull = 0,
    v0TablePostK0sCheck,
    v0TablePostMassCut,
    v0TablePostSelectionCut,
    recoK0sPreSel,
    recoK0sPostSel,
    recoTrackPreSel,
    recoTrackPostSel,
    recoAnalysisDir,
    genAnalysisDir
  };

  static constexpr std::string_view HistRegDire[] = {
    "v0Table/Full/",
    "v0Table/postK0sCheck/",
    "v0Table/postMassCut/",
    "v0Table/postSelectionCut/",
    "recoK0s/PreSel/",
    "recoK0s/PostSel/",
    "recoTracks/PreSel/",
    "recoTracks/PostSel/",
    "recoAnalysis/",
    "genAnalysis/"};

  enum PidEnum {
    kPi = 0, // dont use kPion, kKaon, as these enumeration
    kKa,     // are already defined in $ROOTSYS/root/include/TPDGCode.h
    kPr,
    kEl,
    kDe,
    kK0s
  };

  static constexpr std::string_view PidDire[] = {
    "Pi/",
    "Ka/",
    "Pr/",
    "El/",
    "De/",
    "K0s/"};

  enum DetEnum {
    tpcId = 0,
    tofId,
    tpctofId,
    NoId
  };

  static constexpr std::string_view DetDire[] = {
    "tpcId/",
    "tofId/",
    "tpctofId/",
    "NoId/"};

  // vetoRejection for particles //From Victor Luis Gonzalez Sebastian's analysis note for balance functions
  template <typename T>
  bool selTrackForId(const T& track)
  {
    if (-3.0 < track.tpcNSigmaEl() && track.tpcNSigmaEl() < 5.0 &&
        std::fabs(track.tpcNSigmaPi()) > 3.0 &&
        std::fabs(track.tpcNSigmaKa()) > 3.0 &&
        std::fabs(track.tpcNSigmaPr()) > 3.0) {
      return false;
    } else {
      return true;
    }
  }

  template <int pidMode, typename T>
  bool vetoIdOthersTPC(const T& track)
  {
    if (pidMode != kPi) {
      if (std::fabs(track.tpcNSigmaPi()) < 3.0)
        return false;
    }
    if (pidMode != kKa) {
      if (std::fabs(track.tpcNSigmaKa()) < 3.0)
        return false;
    }
    if (pidMode != kPr) {
      if (std::fabs(track.tpcNSigmaPr()) < 3.0)
        return false;
    }
    if (cfgDoElRejection) {
      if (pidMode != kEl) {
        if (std::fabs(track.tpcNSigmaEl()) < 3.0)
          return false;
      }
    }
    if (cfgDoDeRejection) {
      if (pidMode != kDe) {
        if (std::fabs(track.tpcNSigmaDe()) < 3.0)
          return false;
      }
    }
    return true;
  }

  template <int pidMode, typename T>
  bool vetoIdOthersTOF(const T& track)
  {
    if (pidMode != kPi) {
      if (std::fabs(track.tofNSigmaPi()) < 3.0)
        return false;
    }
    if (pidMode != kKa) {
      if (std::fabs(track.tofNSigmaKa()) < 3.0)
        return false;
    }
    if (pidMode != kPr) {
      if (std::fabs(track.tofNSigmaPr()) < 3.0)
        return false;
    }
    if (cfgDoElRejection) {
      if (pidMode != kEl) {
        if (std::fabs(track.tofNSigmaEl()) < 3.0)
          return false;
      }
    }
    if (cfgDoDeRejection) {
      if (pidMode != kDe) {
        if (std::fabs(track.tofNSigmaDe()) < 3.0)
          return false;
      }
    }
    return true;
  }

  template <int pidMode, typename T>
  bool vetoIdOthersTPCTOF(const T& track)
  {
    if (pidMode != kPi) {
      if (std::fabs(track.tpcNSigmaPi()) < 3.0 && std::fabs(track.tofNSigmaPi()) < 3.0)
        return false;
    }
    if (pidMode != kKa) {
      if (std::fabs(track.tpcNSigmaKa()) < 3.0 && std::fabs(track.tofNSigmaKa()) < 3.0)
        return false;
    }
    if (pidMode != kPr) {
      if (std::fabs(track.tpcNSigmaPr()) < 3.0 && std::fabs(track.tofNSigmaPr()) < 3.0)
        return false;
    }
    if (cfgDoElRejection) {
      if (pidMode != kEl) {
        if (std::fabs(track.tpcNSigmaEl()) < 3.0 && std::fabs(track.tofNSigmaEl()) < 3.0)
          return false;
      }
    }
    if (cfgDoDeRejection) {
      if (pidMode != kDe) {
        if (std::fabs(track.tpcNSigmaDe()) < 3.0 && std::fabs(track.tofNSigmaDe()) < 3.0)
          return false;
      }
    }
    return true;
  }

  template <int pidMode, typename T>
  bool selIdRectangularCut(const T& track, const float& nSigmaTPC, const float& nSigmaTOF)
  {
    switch (pidMode) {
      case kPi:
        if (std::fabs(track.tpcNSigmaPi()) < nSigmaTPC &&
            std::fabs(track.tofNSigmaPi()) < nSigmaTOF) {
          return true;
        }
        break;
      case kKa:
        if (std::fabs(track.tpcNSigmaKa()) < nSigmaTPC &&
            std::fabs(track.tofNSigmaKa()) < nSigmaTOF) {
          return true;
        }
        break;
      case kPr:
        if (std::fabs(track.tpcNSigmaPr()) < nSigmaTPC &&
            std::fabs(track.tofNSigmaPr()) < nSigmaTOF) {
          return true;
        }
        break;
      default:
        return false;
        break;
    }
    return false;
  }

  template <int pidMode, typename T>
  bool selIdEllipsoidalCut(const T& track, const float& nSigmaTPC, const float& nSigmaTOF)
  {
    switch (pidMode) {
      case kPi:
        if (std::pow(track.tpcNSigmaPi() / nSigmaTPC, 2) + std::pow(track.tofNSigmaPi() / nSigmaTOF, 2) < 1.0)
          return true;
        break;
      case kKa:
        if (std::pow(track.tpcNSigmaKa() / nSigmaTPC, 2) + std::pow(track.tofNSigmaKa() / nSigmaTOF, 2) < 1.0)
          return true;
        break;
      case kPr:
        if (std::pow(track.tpcNSigmaPr() / nSigmaTPC, 2) + std::pow(track.tofNSigmaPr() / nSigmaTOF, 2) < 1.0)
          return true;
        break;
      default:
        return false;
        break;
    }
    return false;
  }

  template <int pidMode, typename T>
  bool selIdCircularCut(const T& track, const float& nSigmaSquaredRad)
  {
    switch (pidMode) {
      case kPi:
        if (std::pow(track.tpcNSigmaPi(), 2) + std::pow(track.tofNSigmaPi(), 2) < nSigmaSquaredRad)
          return true;
        break;
      case kKa:
        if (std::pow(track.tpcNSigmaKa(), 2) + std::pow(track.tofNSigmaKa(), 2) < nSigmaSquaredRad)
          return true;
        break;
      case kPr:
        if (std::pow(track.tpcNSigmaPr(), 2) + std::pow(track.tofNSigmaPr(), 2) < nSigmaSquaredRad)
          return true;
        break;
      default:
        return false;
        break;
    }
    return false;
  }

  template <typename T>
  bool checkReliableTOF(const T& track)
  {
    if (track.hasTOF())
      return true; // which check makes the information of TOF relaiable? should track.beta() be checked?
    else
      return false;
  }

  template <int pidMode, typename T>
  bool idTPC(const T& track, const float& nSigmaTPC)
  {
    if (cfgCheckVetoCut && !vetoIdOthersTPC<pidMode>(track))
      return false;
    switch (pidMode) {
      case kPi:
        if (std::fabs(track.tpcNSigmaPi()) < nSigmaTPC)
          return true;
        break;
      case kKa:
        if (std::fabs(track.tpcNSigmaKa()) < nSigmaTPC)
          return true;
        break;
      case kPr:
        if (std::fabs(track.tpcNSigmaPr()) < nSigmaTPC)
          return true;
        break;
      default:
        return false;
        break;
    }
    return false;
  }

  template <int pidMode, typename T>
  bool idTPCTOF(const T& track, const int& pidCutType, const float& nSigmaTPC, const float& nSigmaTOF, const float& nSigmaSquaredRad)
  {
    if (cfgCheckVetoCut && !vetoIdOthersTPCTOF<pidMode>(track))
      return false;
    if (pidCutType == kRectangularCut) {
      return selIdRectangularCut<pidMode>(track, nSigmaTPC, nSigmaTOF);
    } else if (pidCutType == kCircularCut) {
      return selIdCircularCut<pidMode>(track, nSigmaSquaredRad);
    } else if (pidCutType == kEllipsoidalCut) {
      return selIdEllipsoidalCut<pidMode>(track, nSigmaTPC, nSigmaTOF);
    }
    return false;
  }

  template <typename T>
  bool selPiPdependent(const T& track, int& IdMethod)
  {
    if (track.p() < cfgPiThrPforTOF) {
      if (checkReliableTOF(track)) {
        if (idTPCTOF<kPi>(track, cfgPiIdCutTypeLowP, cfgPiNSigmaTPCLowP, cfgPiNSigmaTOFLowP, cfgPiNSigmaRadLowP)) {
          IdMethod = kTPCTOFidentified;
          return true;
        }
        return false;
      } else {
        if (idTPC<kPi>(track, cfgPiNSigmaTPCLowP)) {
          IdMethod = kTPCidentified;
          return true;
        }
        return false;
      }
    } else {
      if (checkReliableTOF(track)) {
        if (idTPCTOF<kPi>(track, cfgPiIdCutTypeHighP, cfgPiNSigmaTPCHighP, cfgPiNSigmaTOFHighP, cfgPiNSigmaRadHighP)) {
          IdMethod = kTPCTOFidentified;
          return true;
        }
        return false;
      }
      return false;
    }
  }

  template <typename T>
  bool selKaPdependent(const T& track, int& IdMethod)
  {
    if (track.p() < cfgKaThrPforTOF) {
      if (checkReliableTOF(track)) {
        if (idTPCTOF<kKa>(track, cfgKaIdCutTypeLowP, cfgKaNSigmaTPCLowP, cfgKaNSigmaTOFLowP, cfgKaNSigmaRadLowP)) {
          IdMethod = kTPCTOFidentified;
          return true;
        }
        return false;
      } else {
        if (idTPC<kKa>(track, cfgKaNSigmaTPCLowP)) {
          IdMethod = kTPCidentified;
          return true;
        }
        return false;
      }
    } else {
      if (checkReliableTOF(track)) {
        if (idTPCTOF<kKa>(track, cfgKaIdCutTypeHighP, cfgKaNSigmaTPCHighP, cfgKaNSigmaTOFHighP, cfgKaNSigmaRadHighP)) {
          IdMethod = kTPCTOFidentified;
          return true;
        }
        return false;
      }
      return false;
    }
  }

  template <typename T>
  bool selPrPdependent(const T& track, int& IdMethod)
  {
    if (track.p() < cfgPrThrPforTOF) {
      if (checkReliableTOF(track)) {
        if (idTPCTOF<kPr>(track, cfgPrIdCutTypeLowP, cfgPrNSigmaTPCLowP, cfgPrNSigmaTOFLowP, cfgPrNSigmaRadLowP)) {
          IdMethod = kTPCTOFidentified;
          return true;
        }
        return false;
      } else {
        if (idTPC<kPr>(track, cfgPrNSigmaTPCLowP)) {
          IdMethod = kTPCidentified;
          return true;
        }
        return false;
      }
    } else {
      if (checkReliableTOF(track)) {
        if (idTPCTOF<kPr>(track, cfgPrIdCutTypeHighP, cfgPrNSigmaTPCHighP, cfgPrNSigmaTOFHighP, cfgPrNSigmaRadHighP)) {
          IdMethod = kTPCTOFidentified;
          return true;
        }
        return false;
      }
      return false;
    }
  }

  //_______________________________Identification Funtions Depending on the tpcInnerParam _______________________________
  // tpc Selections
  template <typename T>
  bool selPionTPCInnerParam(T track)
  {
    if (vetoIdOthersTPC<kPi>(track)) {
      if (0.05 <= track.tpcInnerParam() && track.tpcInnerParam() < 0.70 && std::abs(track.tpcNSigmaPi()) < cfgPiNSigmaTPCLowP) {
        return true;
      }
      if (0.70 <= track.tpcInnerParam() && std::abs(track.tpcNSigmaPi()) < cfgPiNSigmaTPCHighP) {
        return true;
      }
    }
    return false;
  }

  template <typename T>
  bool selKaonTPCInnerParam(T track)
  {
    if (vetoIdOthersTPC<kKa>(track)) {
      if (0.05 <= track.tpcInnerParam() && track.tpcInnerParam() < 0.70 && std::abs(track.tpcNSigmaKa()) < cfgKaNSigmaTPCLowP) {
        return true;
      }
      if (0.70 <= track.tpcInnerParam() && std::abs(track.tpcNSigmaKa()) < cfgKaNSigmaTPCHighP) {
        return true;
      }
    }
    return false;
  }

  template <typename T>
  bool selProtonTPCInnerParam(T track)
  {
    if (vetoIdOthersTPC<kPr>(track)) {
      if (0.05 <= track.tpcInnerParam() && track.tpcInnerParam() < 1.60 && std::abs(track.tpcNSigmaPr()) < cfgPrNSigmaTPCLowP) {
        return true;
      }
      if (1.60 <= track.tpcInnerParam() && std::abs(track.tpcNSigmaPr()) < cfgPrNSigmaTPCHighP) {
        return true;
      }
    }
    return false;
  }

  template <typename T>
  bool selDeuteronTPCInnerParam(T track)
  {
    if (vetoIdOthersTPC<kDe>(track)) {
      if (0.05 <= track.tpcInnerParam() && track.tpcInnerParam() < 1.80 && std::abs(track.tpcNSigmaDe()) < 3.0) {
        return true;
      }
      if (1.80 <= track.tpcInnerParam() && std::abs(track.tpcNSigmaDe()) < 2.0) {
        return true;
      }
    }
    return false;
  }

  template <typename T>
  bool selElectronTPCInnerParam(T track)
  {
    if (track.tpcNSigmaEl() < 3.0 && track.tpcNSigmaPi() > 3.0 && track.tpcNSigmaKa() > 3.0 && track.tpcNSigmaPr() > 3.0 && track.tpcNSigmaDe() > 3.0) {
      return true;
    }
    return false;
  }
  //
  //_____________________________________________________TOF selection Functions _______________________________________________________________________
  // TOF Selections
  // Pion
  template <typename T>
  bool selPionTOF(T track)
  {
    if (vetoIdOthersTOF<kPi>(track)) {
      if (track.p() <= 0.75 && std::abs(track.tpcNSigmaPi()) < cfgPiNSigmaTPCLowP && std::abs(track.tofNSigmaPi()) < cfgPiNSigmaTOFLowP) {
        return true;
      } else if (0.75 < track.p() // after p = 0.75, Pi and Ka lines of nSigma 3.0 will start intersecting
                 && std::abs(track.tpcNSigmaPi()) < cfgPiNSigmaTPCHighP && std::abs(track.tofNSigmaPi()) < cfgPiNSigmaTOFHighP) {
        return true;
      }
    }
    return false;
  }

  // Kaon
  template <typename T>
  bool selKaonTOF(T track)
  {
    if (vetoIdOthersTOF<kKa>(track)) {
      if (track.p() <= 0.75 && std::abs(track.tpcNSigmaKa()) < cfgKaNSigmaTPCLowP && std::abs(track.tofNSigmaKa()) < cfgKaNSigmaTOFLowP) {
        return true;
      }
      if (0.75 < track.p() && track.p() <= 1.30 // after 0.75 Pi and Ka lines of nSigma 3.0 will start intersecting
          && std::abs(track.tpcNSigmaKa()) < cfgKaNSigmaTPCLowP && std::abs(track.tofNSigmaKa()) < cfgKaNSigmaTOFLowP) {
        return true;
      }
      if (1.30 < track.p() // after 1.30 Pr and Ka lines of nSigma 3.0 will start intersecting
          && std::abs(track.tpcNSigmaKa()) < cfgKaNSigmaTPCHighP && std::abs(track.tofNSigmaKa()) < cfgKaNSigmaTOFHighP) {
        return true;
      }
    }
    return false;
  }

  // Proton
  template <typename T>
  bool selProtonTOF(T track)
  {
    if (vetoIdOthersTOF<kPr>(track)) {
      if (track.p() <= 1.30 && std::abs(track.tpcNSigmaPr()) < cfgPrNSigmaTPCLowP && std::abs(track.tofNSigmaPr()) < cfgPrNSigmaTOFLowP) {
        return true;
      }
      if (1.30 < track.p() && track.p() <= 3.10                                                                       // after 1.30 Pr and Ka lines of nSigma 3.0 will start intersecting
          && std::abs(track.tpcNSigmaPr()) < cfgPrNSigmaTPCLowP && std::abs(track.tofNSigmaPr()) < cfgPrNSigmaTOFLowP // Some Deuteron contamination is still coming in p dependent cuts
      ) {
        return true;
      }
      if (3.10 < track.p() // after 3.10 Pr and De lines of nSigma 3.0 will start intersecting
          && std::abs(track.tpcNSigmaPr()) < cfgPrNSigmaTPCHighP && std::abs(track.tofNSigmaPr()) < cfgPrNSigmaTOFHighP) {
        return true;
      }
    }
    return false;
  }

  // Deuteron
  template <typename T>
  bool selDeuteronTOF(T track)
  {
    if (vetoIdOthersTOF<kDe>(track)) {
      if (track.p() <= 3.10 && std::abs(track.tpcNSigmaDe()) < 3.0 && std::abs(track.tofNSigmaDe()) < 3.0) {
        return true;
      }
      if (3.10 < track.p() // after 3.10 De and Pr lines of nSigma 3.0 will start intersecting
          && std::abs(track.tpcNSigmaDe()) < 2.0 && std::abs(track.tofNSigmaDe()) < 2.0) {
        return true;
      }
    }
    return false;
  }

  // Electron
  template <typename T>
  bool selElectronTOF(T track)
  {
    if ((std::pow(track.tpcNSigmaEl(), 2) + std::pow(track.tofNSigmaEl(), 2)) < 9.00 && vetoIdOthersTOF<kEl>(track)) {
      return true;
    }
    return false;
  }
  //

  //______________________________Identification Functions________________________________________________________________
  // Pion
  template <typename T>
  bool selPion(T track, int& IdMethod)
  {
    if (cfgDoPdependentId) {
      return selPiPdependent(track, IdMethod);
    } else if (cfgDoTpcInnerParamId) {
      if (selPionTPCInnerParam(track)) {
        IdMethod = kTPCidentified;
        return true;
      } else if (track.hasTOF() && track.beta() > 0.0 && selPionTOF(track)) {
        IdMethod = kTOFidentified;
        return true;
      }
      return false;
    }
    return false;
  }

  // Kaon
  template <typename T>
  bool selKaon(T track, int& IdMethod)
  {
    if (cfgDoPdependentId) {
      return selKaPdependent(track, IdMethod);
    } else if (cfgDoTpcInnerParamId) {
      if (selKaonTPCInnerParam(track)) {
        IdMethod = kTPCidentified;
        return true;
      } else if (track.hasTOF() && track.beta() > 0.0 && selKaonTOF(track)) {
        IdMethod = kTOFidentified;
        return true;
      }
      return false;
    }
    return false;
  }

  // Proton
  template <typename T>
  bool selProton(T track, int& IdMethod)
  {
    if (cfgDoPdependentId) {
      return selPrPdependent(track, IdMethod);
    } else if (cfgDoTpcInnerParamId) {
      if (selProtonTPCInnerParam(track)) {
        IdMethod = kTPCidentified;
        return true;
      } else if (track.hasTOF() && track.beta() > 0.0 && selProtonTOF(track)) {
        IdMethod = kTOFidentified;
        return true;
      }
      return false;
    }
    return false;
  }

  // Deuteron
  template <typename T>
  bool selDeuteron(T track, int& IdMethod)
  {
    if (cfgDoPdependentId) {
      return false;
    } else if (cfgDoTpcInnerParamId) {
      if (selDeuteronTPCInnerParam(track)) {
        IdMethod = kTPCidentified;
        return true;
      } else if (track.hasTOF() && track.beta() > 0.0 && selDeuteronTOF(track)) {
        IdMethod = kTOFidentified;
        return true;
      }
      return false;
    }
    return false;
  }

  // Electron
  template <typename T>
  bool selElectron(T track, int& IdMethod)
  {
    if (cfgDoPdependentId) {
      return false;
    } else if (cfgDoTpcInnerParamId) {
      if (selElectronTPCInnerParam(track)) {
        IdMethod = kTPCidentified;
        return true;
      } else if (track.hasTOF() && track.beta() > 0.0 && selElectronTOF(track)) {
        IdMethod = kTOFidentified;
        return true;
      }
      return false;
    }
    return false;
  }

  template <typename T>
  int findV0Tag(const T& posDaughterTrack, const T& negDaughterTrack, int& posPiIdMethod, int& posPrIdMethod, int& negPiIdMethod, int& negPrIdMethod)
  {
    bool posIsPion = false;
    bool posIsProton = false;
    bool negIsPion = false;
    bool negIsProton = false;

    int v0TagValue = 0;

    // Check if positive track is pion or proton
    if (selPion(posDaughterTrack, posPiIdMethod))
      posIsPion = true; // Coming From K0s    -> PiPlus + PiMinus and AntiLambda -> PiPlus + AntiProton
    if (selProton(posDaughterTrack, posPrIdMethod))
      posIsProton = true; // Coming From Lambda -> proton + PiMinus
    if (selPion(negDaughterTrack, negPiIdMethod))
      negIsPion = true; // Coming From K0s       -> PiPlus + PiMinus and Lambda -> proton + PiMinus
    if (selProton(negDaughterTrack, negPrIdMethod))
      negIsProton = true; // Coming From AntiLambda -> PiPlus + AntiProton

    if (posIsPion && negIsPion) {
      BITSET(v0TagValue, BIT_IS_K0S);
    }
    if (posIsProton && negIsPion) {
      BITSET(v0TagValue, BIT_IS_LAMBDA);
    }
    if (posIsPion && negIsProton) {
      BITSET(v0TagValue, BIT_IS_ANTILAMBDA);
    }
    return v0TagValue;
  }

  template <typename T, typename U>
  int findCollisionIndexTag(const T& v0, const U& posDaughterTrack, const U& negDaughterTrack)
  {
    int v0daughterCollisionIndexTag = 0;
    if (v0.collisionId() == posDaughterTrack.collisionId()) {
      BITSET(v0daughterCollisionIndexTag, BIT_POS_DAU_HAS_SAME_COLL);
    }
    if (v0.collisionId() == negDaughterTrack.collisionId()) {
      BITSET(v0daughterCollisionIndexTag, BIT_NEG_DAU_HAS_SAME_COLL);
    }
    if (posDaughterTrack.collisionId() == negDaughterTrack.collisionId()) {
      BITSET(v0daughterCollisionIndexTag, BIT_BOTH_DAU_HAS_SAME_COLL);
    }
    return v0daughterCollisionIndexTag;
  }

  template <typename T>
  bool selK0s(T v0)
  {
    if (k0sSelCut.cfgK0sMLow < v0.mK0Short() && v0.mK0Short() < k0sSelCut.cfgK0sMHigh &&
        k0sSelCut.cfgK0sLowPt < v0.pt() && v0.pt() < k0sSelCut.cfgK0sHighPt &&
        std::abs(v0.rapidity(MassK0Short)) < k0sSelCut.cfgK0sRapitidy &&
        v0.qtarm() > (k0sSelCut.cfgK0sARMcut * std::abs(v0.alpha()))) {
      return true;
    } else {
      return false;
    }
  }

  template <typename T>
  void findRepeatedEntries(std::vector<int64_t> ParticleList, T hist)
  {
    for (uint ii = 0; ii < ParticleList.size(); ii++) {
      int nCommonCount = 0; // checking the repeat number of track
      for (uint jj = 0; jj < ParticleList.size(); jj++) {
        if (ParticleList[jj] == ParticleList[ii]) {
          if (jj < ii) {
            break;
          } // break if it was already counted
          nCommonCount++; // To Calculate no of times the entry was repeated
        }
      }
      hist->Fill(nCommonCount);
    }
  }

  template <typename T>
  bool checkTrackSelection(const T& track, int& rejectionTag)
  {
    if (track.tpcNClsCrossedRows() < cfgTrackCuts.cfgTrkTpcNClsCrossedRows) {
      rejectionTag = kFailTpcNClsCrossedRows;
      return false;
    }
    if (std::fabs(track.dcaXY()) > cfgTrackCuts.cfgTrkdcaXY) {
      rejectionTag = kFailTrkdcaXY;
      return false;
    }
    if (!track.isGlobalTrack()) {
      rejectionTag = kFailGlobalTrack;
      return false;
    }
    if (cfgTrackCuts.cfgDoVGselTrackCheck) {
      if (!selTrackForId(track)) {
        rejectionTag = kFailVGSelCheck;
        return false;
      }
    }
    return true;
  }

  template <typename T>
  bool checkTrackInList(const T& track, const std::vector<int64_t>& vecList, int& rejectionTag, const int& listTagValue)
  {
    if (std::binary_search(vecList.begin(), vecList.end(), track.globalIndex())) {
      rejectionTag = listTagValue;
      return true; // Binary Search is fastest search in a sorted array.
    }
    return false;
  }

  template <int Mode, int pidMode, int detMode, bool fillSignal, typename H, typename T>
  void fillIdentificationQA(H histReg, const T& track)
  {

    float tpcNSigmaVal = -999, tofNSigmaVal = -999;
    switch (pidMode) {
      case kPi:
        if (!cfgFillPiQA)
          return;
        tpcNSigmaVal = track.tpcNSigmaPi();
        tofNSigmaVal = track.tofNSigmaPi();
        break;
      case kKa:
        if (!cfgFillKaQA)
          return;
        tpcNSigmaVal = track.tpcNSigmaKa();
        tofNSigmaVal = track.tofNSigmaKa();
        break;
      case kPr:
        if (!cfgFillPrQA)
          return;
        tpcNSigmaVal = track.tpcNSigmaPr();
        tofNSigmaVal = track.tofNSigmaPr();
        break;
      case kEl:
        if (!cfgFillElQA)
          return;
        tpcNSigmaVal = track.tpcNSigmaEl();
        tofNSigmaVal = track.tofNSigmaEl();
        break;
      case kDe:
        if (!cfgFillDeQA)
          return;
        tpcNSigmaVal = track.tpcNSigmaDe();
        tofNSigmaVal = track.tofNSigmaDe();
        break;
      default:
        tpcNSigmaVal = -999, tofNSigmaVal = -999;
        break;
    }

    if (fillSignal) {
      // momemtum
      histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h20_p_pt"), track.p(), track.pt());
      histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h21_p_tpcInnerParam"), track.p(), track.tpcInnerParam());
      histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h22_p_tofExpMom"), track.p(), track.tofExpMom());
      // tpcSignal
      histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h23_p_tpcSignal"), track.p(), track.tpcSignal());
      histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h24_tpcInnerParam_tpcSignal"), track.tpcInnerParam(), track.tpcSignal());
      histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h25_tofExpMom_tpcSignal"), track.tofExpMom(), track.tpcSignal());
      // tofBeta
      histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h26_p_beta"), track.p(), track.beta());
      histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h27_tpcInnerParam_beta"), track.tpcInnerParam(), track.beta());
      histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h28_tofExpMom_beta"), track.tofExpMom(), track.beta());
    }
    // NSigma
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h29_p_tpcNSigma"), track.p(), tpcNSigmaVal);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h30_pt_tpcNSigma"), track.pt(), tpcNSigmaVal);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h31_tpcInnerParam_tpcNSigma"), track.tpcInnerParam(), tpcNSigmaVal);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h32_tofExpMom_tpcNSigma"), track.tofExpMom(), tpcNSigmaVal);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h33_p_tofNSigma"), track.p(), tofNSigmaVal);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h34_pt_tofNSigma"), track.pt(), tofNSigmaVal);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h35_tpcInnerParam_tofNSigma"), track.tpcInnerParam(), tofNSigmaVal);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h36_tofExpMom_tofNSigma"), track.tofExpMom(), tofNSigmaVal);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h37_tpcNSigma_tofNSigma"), tpcNSigmaVal, tofNSigmaVal);
  }

  template <int Mode, typename T>
  void fillTrackQA(T track)
  {
    // FullTrack
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h01_p"), track.p());
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h02_pt"), track.pt());
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h03_tpcInnerParam"), track.tpcInnerParam());
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h04_tofExpMom"), track.tofExpMom());
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h05_eta"), track.eta());
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h06_phi"), track.phi());
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h07_dcaXY"), track.dcaXY());
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h08_dcaZ"), track.dcaZ());
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h09_sign"), track.sign());
    // DcaXY
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h10_p_dcaXY"), track.p(), track.dcaXY());
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h11_pt_dcaXY"), track.pt(), track.dcaXY());
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h12_tpcInnerParam_dcaXY"), track.tpcInnerParam(), track.dcaXY());
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h13_tofExpMom_dcaXY"), track.tofExpMom(), track.dcaXY());

    // DcaZ
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h14_p_dcaZ"), track.p(), track.dcaZ());
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h15_pt_dcaZ"), track.pt(), track.dcaZ());
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h16_tpcInnerParam_dcaZ"), track.tpcInnerParam(), track.dcaZ());
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h17_tofExpMom_dcaZ"), track.tofExpMom(), track.dcaZ());

    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h20_pt_eta"), track.pt(), track.eta());
    // momemtum
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h20_p_pt"), track.p(), track.pt());
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h21_p_tpcInnerParam"), track.p(), track.tpcInnerParam());
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h22_p_tofExpMom"), track.p(), track.tofExpMom());

    // tpcSignal
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h23_p_tpcSignal"), track.p(), track.tpcSignal());
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h24_tpcInnerParam_tpcSignal"), track.tpcInnerParam(), track.tpcSignal());
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h25_tofExpMom_tpcSignal"), track.tofExpMom(), track.tpcSignal());

    // tofBeta
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h26_p_beta"), track.p(), track.beta());
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h27_tpcInnerParam_beta"), track.tpcInnerParam(), track.beta());
    recoTracks.fill(HIST(HistRegDire[Mode]) + HIST("h28_tofExpMom_beta"), track.tofExpMom(), track.beta());

    fillIdentificationQA<Mode, kPi, NoId, false>(recoTracks, track); // Look at Pion
    fillIdentificationQA<Mode, kKa, NoId, false>(recoTracks, track); // Look at Kaon
    fillIdentificationQA<Mode, kPr, NoId, false>(recoTracks, track); // Look at Proton
    fillIdentificationQA<Mode, kEl, NoId, false>(recoTracks, track); // Look at Electron
    fillIdentificationQA<Mode, kDe, NoId, false>(recoTracks, track); // Look at Deuteron
  }

  template <int Mode, int pidMode, int detMode, typename H, typename T>
  void fillV0DaughterQA(H histReg, const T& track, double particleMass)
  {
    // K0s-Daughter Info
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h01_p"), track.p());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h02_pt"), track.pt());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h03_tpcInnerParam"), track.tpcInnerParam());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h04_tofExpMom"), track.tofExpMom());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h05_eta"), track.eta());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h06_phi"), track.phi());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h07_rapidity"), track.rapidity(particleMass));
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h08_isPVContributor"), track.isPVContributor());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h09_isGlobalTrack"), track.isGlobalTrack());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h10_dcaXY"), track.dcaXY());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h11_dcaZ"), track.dcaZ());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h12_p_dcaXY"), track.p(), track.dcaXY());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h13_p_dcaZ"), track.p(), track.dcaZ());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h14_pt_dcaXY"), track.pt(), track.dcaXY());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h15_pt_dcaZ"), track.pt(), track.dcaZ());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h16_dcaXYwide"), track.dcaXY());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h17_dcaZwide"), track.dcaZ());

    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("h20_pt_eta"), track.pt(), track.eta());

    fillIdentificationQA<Mode, kPi, detMode, true>(histReg, track);
  }

  template <int Mode, typename H, typename T, typename U>
  void fillV0QA(H histReg, const T& v0, const U& posDaughterTrack, const U& negDaughterTrack, const int& v0Tag, const int& v0DauCollisionIndexTag, const int& posPiIdMethod, const int& negPiIdMethod)
  {
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h01_K0s_Mass"), v0.mK0Short());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h02_Lambda_Mass"), v0.mLambda());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h03_AntiLambda_Mass"), v0.mAntiLambda());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h04_v0DaughterCollisionIndexTag"), v0DauCollisionIndexTag);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h05_V0Tag"), v0Tag);

    // Topological Cuts
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h06_dcapostopv"), v0.dcapostopv());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h07_dcanegtopv"), v0.dcanegtopv());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h08_dcaV0daughters"), v0.dcaV0daughters());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h09_v0cosPA"), v0.v0cosPA());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h10_v0radius"), v0.v0radius());

    // K0s-FullInformation
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h11_mass"), v0.mK0Short());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h12_p"), v0.p());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h13_pt"), v0.pt());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h14_eta"), v0.eta());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h15_phi"), v0.phi());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h16_rapidity"), v0.rapidity(MassK0Short));
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h17_alpha"), v0.alpha());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h18_qtarm"), v0.qtarm());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h19_alpha_qtarm"), v0.alpha(), v0.qtarm());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h20_pt_eta"), v0.pt(), v0.eta());

    if (posPiIdMethod == kTPCidentified) {
      fillV0DaughterQA<Mode, kPi, tpcId>(histReg, posDaughterTrack, MassProton);
    } else if (posPiIdMethod == kTPCTOFidentified) {
      fillV0DaughterQA<Mode, kPi, tpctofId>(histReg, posDaughterTrack, MassProton);
    } else if (posPiIdMethod == kUnidentified) {
      fillV0DaughterQA<Mode, kPi, NoId>(histReg, posDaughterTrack, MassProton);
    }

    if (negPiIdMethod == kTPCidentified) {
      fillV0DaughterQA<Mode, kPi, tpcId>(histReg, negDaughterTrack, MassProton);
    } else if (negPiIdMethod == kTPCTOFidentified) {
      fillV0DaughterQA<Mode, kPi, tpctofId>(histReg, negDaughterTrack, MassProton);
    } else if (negPiIdMethod == kUnidentified) {
      fillV0DaughterQA<Mode, kPi, NoId>(histReg, negDaughterTrack, MassProton);
    }
  }

  template <int Mode, int pidMode, typename H, typename T>
  void fillGenTrackQA(H& histReg, const T& mcTrack)
  {
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST("h12_p"), mcTrack.p());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST("h13_pt"), mcTrack.pt());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST("h14_eta"), mcTrack.eta());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST("h15_phi"), mcTrack.phi());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST("h16_rapidity"), mcTrack.y());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST("h20_pt_eta"), mcTrack.pt(), mcTrack.eta());
  }

  template <typename T, typename U, typename H>
  void executeV0loop(const T& posDaughterTrack, const T& negDaughterTrack, const U& v0, H& recoV0s,
                     int& posPiIdMethod, int& posPrIdMethod, int& negPiIdMethod, int& negPrIdMethod,
                     int& v0Tag, int& trueV0TagValue, bool& isK0s, int& v0DauCollisionIndexTag,
                     auto& k0sPosDauList, auto& k0sNegDauList)
  {
    posPiIdMethod = kUnidentified;
    posPrIdMethod = kUnidentified;
    negPiIdMethod = kUnidentified;
    negPrIdMethod = kUnidentified;
    v0Tag = findV0Tag(posDaughterTrack, negDaughterTrack, posPiIdMethod, posPrIdMethod, negPiIdMethod, negPrIdMethod);
    v0DauCollisionIndexTag = findCollisionIndexTag(v0, posDaughterTrack, negDaughterTrack);

    isK0s = false;
    if (BITCHECK(v0Tag, BIT_IS_K0S) != 0)
      isK0s = true;
    trueV0TagValue = 0;
    if (cfgFillV0TableFull) {
      fillV0QA<v0TableFull>(recoV0s, v0, posDaughterTrack, negDaughterTrack, v0Tag, v0DauCollisionIndexTag, posPiIdMethod, negPiIdMethod);
    }
    // cut on dynamic columns for v0 particles
    if (v0.v0cosPA() < v0settingCosPA)
      return; // in place of continue;
    if (v0.v0radius() < v0settingRadius)
      return; // in place of continue;

    // K0s Analysis
    if (isK0s) {
      if (cfgFillV0TablePostK0sCheck) {
        fillV0QA<v0TablePostK0sCheck>(recoV0s, v0, posDaughterTrack, negDaughterTrack, v0Tag, v0DauCollisionIndexTag, posPiIdMethod, negPiIdMethod);
      }
      // K0s mass cut
      if (cfgFillV0TablePostMassCut) {
        if (k0sSelCut.cfgK0sMLow < v0.mK0Short() && v0.mK0Short() < k0sSelCut.cfgK0sMHigh) {
          fillV0QA<v0TablePostMassCut>(recoV0s, v0, posDaughterTrack, negDaughterTrack, v0Tag, v0DauCollisionIndexTag, posPiIdMethod, negPiIdMethod);
        }
      }
      // Final K0s Selection.
      if (selK0s(v0)) {
        if (cfgFillV0TablePostSelectionCut) {
          fillV0QA<v0TablePostSelectionCut>(recoV0s, v0, posDaughterTrack, negDaughterTrack, v0Tag, v0DauCollisionIndexTag, posPiIdMethod, negPiIdMethod);
        }
        trueV0TagValue += 1;
        k0sPosDauList.push_back(posDaughterTrack.globalIndex());
        k0sNegDauList.push_back(negDaughterTrack.globalIndex());
      }
      recoV0s.fill(HIST(HistRegDire[v0TablePostSelectionCut]) + HIST("hTrueV0TagCount"), trueV0TagValue); // 001 = Kaon, 010 = Lambda, 100 = AnitLambda
    } // End of K0s block
  }

  template <typename T>
  void executeSortPairDaughters(const auto& posDauList, const auto& negDauList, auto& fullDauList, const T& hist)
  {
    findRepeatedEntries(posDauList, hist);
    findRepeatedEntries(negDauList, hist);

    // Obtain one single new daughter vector to remove double counting
    fullDauList.insert(fullDauList.end(), posDauList.begin(), posDauList.end());
    fullDauList.insert(fullDauList.end(), negDauList.begin(), negDauList.end());

    // Sort and Remove repeated entries
    std::sort(fullDauList.begin(), fullDauList.end());
    auto last = std::unique(fullDauList.begin(), fullDauList.end()); // std::unique only moves duplicates to end of the vector
    fullDauList.erase(last, fullDauList.end());                      // last is the iterator position from where duplicate entries start

    // Check sorting
    if (!std::is_sorted(fullDauList.begin(), fullDauList.end())) {
      LOG(error) << "fullDauList is unsorted, will give wrong results when v0 and collisions will be checked";
    }
  }

  template <typename T, typename U> //, typename H>
  void executeV0InCollisionloop(const T& posDaughterTrack, const T& negDaughterTrack, const U& v0,
                                int& posPiIdMethod, int& posPrIdMethod, int& negPiIdMethod, int& negPrIdMethod,
                                int& v0Tag, bool& isK0s, int& v0DauCollisionIndexTag, float& nK0s, const float& centrality)
  {
    if (v0.v0cosPA() < v0settingCosPA)
      return; // for continue; // cut on dynamic columns for v0 particles
    if (v0.v0radius() < v0settingRadius)
      return; // for continue;
    isK0s = false;

    posPiIdMethod = kUnidentified;
    posPrIdMethod = kUnidentified;
    negPiIdMethod = kUnidentified;
    negPrIdMethod = kUnidentified;
    v0Tag = findV0Tag(posDaughterTrack, negDaughterTrack, posPiIdMethod, posPrIdMethod, negPiIdMethod, negPrIdMethod);
    v0DauCollisionIndexTag = findCollisionIndexTag(v0, posDaughterTrack, negDaughterTrack);

    if (BITCHECK(v0Tag, BIT_IS_K0S) != 0)
      isK0s = true;
    if (cfgFillRecoK0sPreSel) {
      fillV0QA<recoK0sPreSel>(recoK0s, v0, posDaughterTrack, negDaughterTrack, v0Tag, v0DauCollisionIndexTag, posPiIdMethod, negPiIdMethod);
    }
    // K0s Analysis
    if (isK0s && selK0s(v0)) {
      if (cfgFillRecoK0sPostSel) {
        fillV0QA<recoK0sPostSel>(recoK0s, v0, posDaughterTrack, negDaughterTrack, v0Tag, v0DauCollisionIndexTag, posPiIdMethod, negPiIdMethod);
      }
      recoK0s.fill(HIST(HistRegDire[recoK0sPostSel]) + HIST("mK0s_vs_cent"), centrality, v0.mK0Short()); // centrality dependent mass
      nK0s++;
    } // End of K0s block
  }

  template <typename T>
  void executeTrackQAPart(const T& track, const auto& fullDauList, int& rejectionTag, float& nRejectedPiMinus, float& nRejectedPiPlus, int& nTrack, bool& isAcceptedTrack)
  {
    if (cfgFillRecoTrackPreSel) {
      fillTrackQA<recoTrackPreSel>(track);
    }
    rejectionTag = 0;
    if (!checkTrackSelection(track, rejectionTag)) {
      recoAnalysis.fill(HIST("recoAnalysis/RejectedTrack_RejectionTag"), rejectionTag);
      isAcceptedTrack = false;
      return; // for continue;
    } else if (checkTrackInList(track, fullDauList, rejectionTag, kFailK0ShortDaughter)) {
      recoAnalysis.fill(HIST("recoAnalysis/RejectedTrack_RejectionTag"), rejectionTag);
      if (track.signed1Pt() > 0) {
        nRejectedPiPlus++; // DOEFFCORR
      } else if (track.signed1Pt() < 0) {
        nRejectedPiMinus++; // DOEFFCORR
      }
      isAcceptedTrack = false;
      return; // for continue;
    }
    isAcceptedTrack = true;
    if (cfgFillRecoTrackPostSel) {
      fillTrackQA<recoTrackPostSel>(track);
    }
    nTrack++;
  }

  template <typename T>
  void executeTrackAnalysisPart(const T& track, const int& trackIdTag,
                                const int& idMethodPi, const bool& trackIsPion, float& nPiMinus, float& nPiPlus,
                                const int& idMethodKa, const bool& trackIsKaon, float& nKaMinus, float& nKaPlus,
                                const int& idMethodPr, const bool& trackIsProton, float& nProton, float& nPBar,
                                const int& idMethodEl, const bool& trackIsElectron, float& nElPlus, float& nElMinus,
                                const int& idMethodDe, const bool& trackIsDeuteron, float& nDePlus, float& nDeMinus)
  {
    if (trackIsPion) {
      if (idMethodPi == kTPCidentified) {
        fillIdentificationQA<recoAnalysisDir, kPi, tpcId, true>(recoAnalysis, track);
      } else if (idMethodPi == kTPCTOFidentified) {
        fillIdentificationQA<recoAnalysisDir, kPi, tpctofId, true>(recoAnalysis, track);
      } else if (idMethodPi == kUnidentified) {
        fillIdentificationQA<recoAnalysisDir, kPi, NoId, true>(recoAnalysis, track);
      }
      if (track.sign() > 0) {
        nPiPlus++;
      } else if (track.sign() < 0) {
        nPiMinus++;
      }
    }
    if (trackIsKaon) {
      if (idMethodKa == kTPCidentified) {
        fillIdentificationQA<recoAnalysisDir, kKa, tpcId, true>(recoAnalysis, track);
      } else if (idMethodKa == kTPCTOFidentified) {
        fillIdentificationQA<recoAnalysisDir, kKa, tpctofId, true>(recoAnalysis, track);
      } else if (idMethodKa == kUnidentified) {
        fillIdentificationQA<recoAnalysisDir, kKa, NoId, true>(recoAnalysis, track);
      }
      if (track.sign() > 0) {
        nKaPlus++;
      } else if (track.sign() < 0) {
        nKaMinus++;
      }
    }
    if (trackIsProton) {
      if (idMethodPr == kTPCidentified) {
        fillIdentificationQA<recoAnalysisDir, kPr, tpcId, true>(recoAnalysis, track);
      } else if (idMethodPr == kTPCTOFidentified) {
        fillIdentificationQA<recoAnalysisDir, kPr, tpctofId, true>(recoAnalysis, track);
      } else if (idMethodPr == kUnidentified) {
        fillIdentificationQA<recoAnalysisDir, kPr, NoId, true>(recoAnalysis, track);
      }
      if (track.sign() > 0) {
        nProton++;
      } else if (track.sign() < 0) {
        nPBar++;
      }
    }
    if (trackIsElectron) {
      if (idMethodEl == kTPCidentified) {
        fillIdentificationQA<recoAnalysisDir, kEl, tpcId, true>(recoAnalysis, track);
      } else if (idMethodEl == kTPCTOFidentified) {
        fillIdentificationQA<recoAnalysisDir, kEl, tpctofId, true>(recoAnalysis, track);
      } else if (idMethodEl == kUnidentified) {
        fillIdentificationQA<recoAnalysisDir, kEl, NoId, true>(recoAnalysis, track);
      }
      if (track.sign() > 0) {
        nElPlus++;
      } else if (track.sign() < 0) {
        nElMinus++;
      }
    }
    if (trackIsDeuteron) {
      if (idMethodDe == kTPCidentified) {
        fillIdentificationQA<recoAnalysisDir, kDe, tpcId, true>(recoAnalysis, track);
      } else if (idMethodDe == kTPCTOFidentified) {
        fillIdentificationQA<recoAnalysisDir, kDe, tpctofId, true>(recoAnalysis, track);
      } else if (idMethodDe == kUnidentified) {
        fillIdentificationQA<recoAnalysisDir, kDe, NoId, true>(recoAnalysis, track);
      }
      if (track.sign() > 0) {
        nDePlus++;
      }
      if (track.sign() < 0) {
        nDeMinus++;
      }
    }
    recoAnalysis.fill(HIST("recoAnalysis/SelectedTrack_IdentificationTag"), trackIdTag);
  }

  void executeSparseAnalysisPart(const float& centFT0C, const float& nTrack, const float& nK0s,
                                 const float& nRejectedPiPlus, const float& nRejectedPiMinus, float& nKaon,
                                 const float& nPiPlus, const float& nKaPlus, const float& nProton, const float& nElPlus, const float& nDePlus,
                                 const float& nPiMinus, const float& nKaMinus, const float& nPBar, const float& nElMinus, const float& nDeMinus)
  {
    nKaon = nKaPlus + nKaMinus;
    if (cfgFillSparseFullK0sPiKa) {
      recoAnalysis.fill(HIST("recoAnalysis/Sparse_Full_K0sPiKa"),
                        centFT0C, nTrack, nK0s,
                        nRejectedPiPlus, nRejectedPiMinus,
                        nPiPlus, nPiMinus, nKaPlus, nKaMinus);
    }
    if (cfgFillSparseFullK0sPrDe) {
      recoAnalysis.fill(HIST("recoAnalysis/Sparse_Full_K0sPrDe"),
                        centFT0C, nTrack, nK0s,
                        nRejectedPiPlus, nRejectedPiMinus,
                        nProton, nPBar, nDePlus, nDeMinus);
    }
    if (cfgFillSparseFullK0sKaEl) {
      recoAnalysis.fill(HIST("recoAnalysis/Sparse_Full_K0sKaEl"),
                        centFT0C, nTrack, nK0s, nRejectedPiPlus, nRejectedPiMinus,
                        nKaPlus, nKaMinus, nElPlus, nElMinus);
    }
    if (cfgFillSparseFullPiKaPr) {
      recoAnalysis.fill(HIST("recoAnalysis/Sparse_Full_PiKaPr"),
                        centFT0C, nTrack,
                        nRejectedPiPlus, nRejectedPiMinus,
                        nPiPlus, nPiMinus, nKaPlus, nKaMinus, nProton, nPBar);
    }
    if (cfgFillSparseFullPiElDe) {
      recoAnalysis.fill(HIST("recoAnalysis/Sparse_Full_PiElDe"),
                        centFT0C, nTrack,
                        nRejectedPiPlus, nRejectedPiMinus,
                        nPiPlus, nPiMinus, nElPlus, nElMinus, nDePlus, nDeMinus);
    }
    if (cfgFillSparseFullKaPrDe) {
      recoAnalysis.fill(HIST("recoAnalysis/Sparse_Full_KaPrDe"),
                        centFT0C, nTrack,
                        nRejectedPiPlus, nRejectedPiMinus,
                        nKaPlus, nKaMinus, nProton, nPBar, nDePlus, nDeMinus);
    }
    if (cfgFillSparseFullPrElDe) {
      recoAnalysis.fill(HIST("recoAnalysis/Sparse_Full_PrElDe"),
                        centFT0C, nTrack,
                        nRejectedPiPlus, nRejectedPiMinus,
                        nProton, nPBar, nElPlus, nElMinus, nDePlus, nDeMinus);
    }
    if (cfgFillSparsenewDynmK0sKa) {
      recoAnalysis.fill(HIST("recoAnalysis/Sparse_newDynm_K0s_Ka"),
                        centFT0C, nTrack, nK0s, nKaon,
                        nK0s * nK0s, nKaon * nKaon, nK0s * nKaon);
    }
    if (cfgFillSparsenewDynmKpKm) {
      recoAnalysis.fill(HIST("recoAnalysis/Sparse_newDynm_Kp_Km"),
                        centFT0C, nTrack, nKaPlus, nKaMinus,
                        nKaPlus * nKaPlus, nKaMinus * nKaMinus, nKaPlus * nKaMinus);
    }
  }

  template <typename C, typename T>
  void executeEventInfoPart(const C& collision, const float& centrality, const int v0TableSize, const T& tracksTablePerColl,
                            const int& nTrack, const int& nK0s,
                            const int& nPiPlus, const int& nKaPlus, const int& nProton, const int& nElPlus, const int& nDePlus,
                            const int& nPiMinus, const int& nKaMinus, const int& nPBar, const int& nElMinus, const int& nDeMinus)
  {
    // Collisions QA
    recoEvent.fill(HIST("recoEvent/h01_CollisionCount"), 0.5);
    recoEvent.fill(HIST("recoEvent/h02_VertexXRec"), collision.posX());
    recoEvent.fill(HIST("recoEvent/h03_VertexYRec"), collision.posY());
    recoEvent.fill(HIST("recoEvent/h04_VertexZRec"), collision.posZ());
    recoEvent.fill(HIST("recoEvent/h05_Centrality"), centrality);
    recoEvent.fill(HIST("recoEvent/h06_V0Size"), v0TableSize);
    recoEvent.fill(HIST("recoEvent/h07_TracksSize"), tracksTablePerColl.size());
    recoEvent.fill(HIST("recoEvent/h08_nTrack"), nTrack);
    recoEvent.fill(HIST("recoEvent/h09_nK0s"), nK0s);
    recoEvent.fill(HIST("recoEvent/h10_nPiPlus"), nPiPlus);
    recoEvent.fill(HIST("recoEvent/h11_nPiMinus"), nPiMinus);
    recoEvent.fill(HIST("recoEvent/h12_nKaPlus"), nKaPlus);
    recoEvent.fill(HIST("recoEvent/h13_nKaMinus"), nKaMinus);
    recoEvent.fill(HIST("recoEvent/h14_nProton"), nProton);
    recoEvent.fill(HIST("recoEvent/h15_nPBar"), nPBar);
    recoEvent.fill(HIST("recoEvent/h16_nElPlus"), nElPlus);
    recoEvent.fill(HIST("recoEvent/h17_nElMinus"), nElMinus);
    recoEvent.fill(HIST("recoEvent/h18_nDePlus"), nDePlus);
    recoEvent.fill(HIST("recoEvent/h19_nDeMinus"), nDeMinus);
  }

  // Event Filter
  Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < cutZvertex);

  // Track Filter
  Filter ptFilter = (o2::aod::track::pt) > cfgTrackCuts.cfgTrackPtLow && (o2::aod::track::pt) < cfgTrackCuts.cfgTrackPtHigh;
  Filter etaFilter = (nabs(o2::aod::track::eta) < cfgTrackCuts.cfgTrackEta);

  // Filters on V0s
  Filter preFilterv0 = (nabs(aod::v0data::dcapostopv) > v0settingDcaPosToPV &&
                        nabs(aod::v0data::dcanegtopv) > v0settingDcaNegToPV &&
                        aod::v0data::dcaV0daughters < v0settingDcaV0Dau);

  using MyCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels,
                                               aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::Mults>>;

  using MyTracks = soa::Filtered<soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                           aod::TOFSignal, aod::pidTOFbeta, aod::pidTOFmass,
                                           aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullEl, aod::pidTPCFullDe,
                                           aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullEl, aod::pidTOFFullDe>>;

  using MyV0s = soa::Filtered<aod::V0Datas>;

  // For manual sliceBy
  Preslice<MyTracks> tracksPerCollisionPreslice = o2::aod::track::collisionId;
  Preslice<MyV0s> v0sPerCollisionPreslice = o2::aod::track::collisionId;

  using MyCollisionsWithMcLabels = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels,
                                                           aod::CentFT0Ms, aod::CentFT0Cs, aod::Mults, aod::McCollisionLabels>>;

  using MyTracksWithMcLabels = soa::Filtered<soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                                       aod::TOFSignal, aod::pidTOFbeta, // aod::pidTOFmass,
                                                       aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullEl, aod::pidTPCFullDe,
                                                       aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullEl, aod::pidTOFFullDe, aod::McTrackLabels>>;

  using MyV0sWithMcLabels = soa::Filtered<soa::Join<aod::V0Datas, aod::McV0Labels>>;

  // Declaring vectors outside the process to avoid slight overhead for stack allocation and deallocation during each iteration.
  std::vector<int64_t> k0sPosDauList;
  std::vector<int64_t> k0sNegDauList;
  std::vector<int64_t> k0sFullDauList;

  template <int analysisType, typename C, typename V, typename T>
  void executeAnalysis(const C& collisions, const V& V0s, const T& tracks)
  {
    k0sPosDauList.clear();
    k0sNegDauList.clear();
    k0sFullDauList.clear();

    int posPiIdMethod = kUnidentified;
    int posPrIdMethod = kUnidentified;
    int negPiIdMethod = kUnidentified;
    int negPrIdMethod = kUnidentified;

    bool isK0s = false;

    int v0Tag = 0;
    int trueV0TagValue = 0;
    int v0DauCollisionIndexTag = 0;

    // Declaring variables outside the loop to avoid slight overhead for stack allocation and deallocation during each iteration.

    float nK0s = 0;
    float nPiPlus = 0;
    float nPiMinus = 0;
    float nKaPlus = 0;
    float nKaMinus = 0;
    float nProton = 0;
    float nPBar = 0;
    float nElPlus = 0;
    float nElMinus = 0;
    float nDePlus = 0;
    float nDeMinus = 0;

    float nKaon = 0;
    int nTrack = 0;

    float centrality = 0;

    float nRejectedPiPlus = 0;
    float nRejectedPiMinus = 0;
    int rejectionTag = 0;

    bool trackIsPion = false;
    bool trackIsKaon = false;
    bool trackIsProton = false;
    bool trackIsElectron = false;
    bool trackIsDeuteron = false;

    int trackIdTag = 0;
    int idMethodPi = kUnidentified;
    int idMethodKa = kUnidentified;
    int idMethodPr = kUnidentified;
    int idMethodEl = kUnidentified;
    int idMethodDe = kUnidentified;

    if constexpr (analysisType == doDataProcessing) {
      for (const auto& v0 : V0s) {
        const auto& posDaughterTrack = v0.template posTrack_as<MyTracks>();
        const auto& negDaughterTrack = v0.template negTrack_as<MyTracks>();

        executeV0loop(posDaughterTrack, negDaughterTrack, v0, recoV0s,
                      posPiIdMethod, posPrIdMethod, negPiIdMethod, negPrIdMethod,
                      v0Tag, trueV0TagValue, isK0s, v0DauCollisionIndexTag,
                      k0sPosDauList, k0sNegDauList);
      } // End of V0s Loop

      executeSortPairDaughters(k0sPosDauList, k0sNegDauList, k0sFullDauList, recoV0s.get<TH1>(HIST(HistRegDire[v0TablePostSelectionCut]) + HIST("nCommonPionOfDifferentK0s")));

      for (const auto& collision : collisions) {

        nK0s = 0;
        nPiPlus = 0;
        nPiMinus = 0;
        nKaPlus = 0;
        nKaMinus = 0;
        nProton = 0;
        nPBar = 0;
        nElPlus = 0;
        nElMinus = 0;
        nDePlus = 0;
        nDeMinus = 0;
        nTrack = 0;
        nKaon = 0;

        centrality = collision.centFT0C();
        if (cfgCentAxis.centAxisType == 1) {
          centrality = collision.centFT0M();
        } else if (cfgCentAxis.centAxisType == 2) {
          centrality = collision.multFT0M();
        } else if (cfgCentAxis.centAxisType == 3) {
          centrality = collision.multFT0C();
        }

        // group tracks, v0s manually
        const uint64_t collIdx = collision.globalIndex();
        const auto tracksTablePerColl = tracks.sliceBy(tracksPerCollisionPreslice, collIdx);
        const auto v0sTablePerColl = V0s.sliceBy(v0sPerCollisionPreslice, collIdx);

        for (const auto& v0 : v0sTablePerColl) {
          const auto& posDaughterTrack = v0.template posTrack_as<MyTracks>();
          const auto& negDaughterTrack = v0.template negTrack_as<MyTracks>();

          executeV0InCollisionloop(posDaughterTrack, negDaughterTrack, v0,
                                   posPiIdMethod, posPrIdMethod, negPiIdMethod, negPrIdMethod,
                                   v0Tag, isK0s, v0DauCollisionIndexTag, nK0s, centrality);

        } // End of V0s Loop

        nTrack = 0;
        nRejectedPiPlus = 0;
        nRejectedPiMinus = 0;
        for (const auto& track : tracksTablePerColl) {
          bool isAcceptedTrack = true;
          executeTrackQAPart(track, k0sFullDauList, rejectionTag, nRejectedPiMinus, nRejectedPiPlus, nTrack, isAcceptedTrack);
          if (!isAcceptedTrack) {
            continue;
          }

          // Do Proper Track Identification
          trackIsPion = false;
          trackIsKaon = false;
          trackIsProton = false;
          trackIsElectron = false;
          trackIsDeuteron = false;

          trackIdTag = 0;
          idMethodPi = kUnidentified;
          idMethodKa = kUnidentified;
          idMethodPr = kUnidentified;
          idMethodEl = kUnidentified;
          idMethodDe = kUnidentified;

          if (selPion(track, idMethodPi)) {
            trackIsPion = true;
            BITSET(trackIdTag, ID_BIT_PI);
          }
          if (selKaon(track, idMethodKa)) {
            trackIsKaon = true;
            BITSET(trackIdTag, ID_BIT_KA);
          }
          if (selProton(track, idMethodPr)) {
            trackIsProton = true;
            BITSET(trackIdTag, ID_BIT_PR);
          }
          if (selElectron(track, idMethodEl)) {
            trackIsElectron = true;
            BITSET(trackIdTag, ID_BIT_EL);
          }
          if (selDeuteron(track, idMethodDe)) {
            trackIsDeuteron = true;
            BITSET(trackIdTag, ID_BIT_DE);
          }

          executeTrackAnalysisPart(track, trackIdTag,
                                   idMethodPi, trackIsPion, nPiMinus, nPiPlus,
                                   idMethodKa, trackIsKaon, nKaMinus, nKaPlus,
                                   idMethodPr, trackIsProton, nProton, nPBar,
                                   idMethodEl, trackIsElectron, nElPlus, nElMinus,
                                   idMethodDe, trackIsDeuteron, nDePlus, nDeMinus);
        } // track loop ends

        executeSparseAnalysisPart(centrality, nTrack, nK0s,
                                  nRejectedPiPlus, nRejectedPiMinus, nKaon,
                                  nPiPlus, nKaPlus, nProton, nElPlus, nDePlus,
                                  nPiMinus, nKaMinus, nPBar, nElMinus, nDeMinus);

        executeEventInfoPart(collision, centrality, v0sTablePerColl.size(), tracksTablePerColl,
                             nTrack, nK0s,
                             nPiPlus, nKaPlus, nProton, nElPlus, nDePlus,
                             nPiMinus, nKaMinus, nPBar, nElMinus, nDeMinus);
      } // collision loop ends
    } else if constexpr (analysisType == doRecoProcessing || analysisType == doPurityProcessing) {
      for (const auto& v0 : V0s) {
        const auto& posDaughterTrack = v0.template posTrack_as<MyTracksWithMcLabels>();
        const auto& negDaughterTrack = v0.template negTrack_as<MyTracksWithMcLabels>();

        //__________________________Reco Level ____________________________________________________
        if (!v0.has_mcParticle() || !posDaughterTrack.has_mcParticle() || !negDaughterTrack.has_mcParticle()) {
          continue;
        }

        auto v0mcparticle = v0.mcParticle(); // if (v0mcparticle.pdgCode() != 310 || !v0mcparticle.isPhysicalPrimary())
        if (!v0mcparticle.isPhysicalPrimary())
          continue;

        if constexpr (analysisType == doPurityProcessing) {
          auto posDauMcPart = posDaughterTrack.mcParticle();
          auto negDauMcPart = negDaughterTrack.mcParticle();
          if (v0mcparticle.pdgCode() != kK0Short || posDauMcPart.pdgCode() != kPiPlus || negDauMcPart.pdgCode() != kPiMinus)
            continue;
        }

        executeV0loop(posDaughterTrack, negDaughterTrack, v0, recoV0s,
                      posPiIdMethod, posPrIdMethod, negPiIdMethod, negPrIdMethod,
                      v0Tag, trueV0TagValue, isK0s, v0DauCollisionIndexTag,
                      k0sPosDauList, k0sNegDauList);
      } // End of V0s Loop
      executeSortPairDaughters(k0sPosDauList, k0sNegDauList, k0sFullDauList, recoV0s.get<TH1>(HIST(HistRegDire[v0TablePostSelectionCut]) + HIST("nCommonPionOfDifferentK0s")));

      for (const auto& collision : collisions) {
        if (!collision.has_mcCollision()) {
          LOG(warning) << "No MC collision for this collision, skip...";
          continue;
        }

        nK0s = 0;
        nPiPlus = 0;
        nPiMinus = 0;
        nKaPlus = 0;
        nKaMinus = 0;
        nProton = 0;
        nPBar = 0;
        nElPlus = 0;
        nElMinus = 0;
        nDePlus = 0;
        nDeMinus = 0;
        nTrack = 0;
        nKaon = 0;

        centrality = collision.centFT0C();
        if (cfgCentAxis.centAxisType == 1) {
          centrality = collision.centFT0M();
        } else if (cfgCentAxis.centAxisType == 2) {
          centrality = collision.multFT0M();
        } else if (cfgCentAxis.centAxisType == 3) {
          centrality = collision.multFT0C();
        }

        // group tracks, v0s manually
        const uint64_t collIdx = collision.globalIndex();
        const auto tracksTablePerColl = tracks.sliceBy(tracksPerCollisionPreslice, collIdx);
        const auto v0sTablePerColl = V0s.sliceBy(v0sPerCollisionPreslice, collIdx);

        for (const auto& v0 : v0sTablePerColl) {
          const auto& posDaughterTrack = v0.template posTrack_as<MyTracksWithMcLabels>();
          const auto& negDaughterTrack = v0.template negTrack_as<MyTracksWithMcLabels>();

          //__________________________Reco Level ____________________________________________________
          if (!v0.has_mcParticle() || !posDaughterTrack.has_mcParticle() || !negDaughterTrack.has_mcParticle()) {
            continue;
          }

          auto v0mcparticle = v0.mcParticle();
          if (!v0mcparticle.isPhysicalPrimary())
            continue;

          if constexpr (analysisType == doPurityProcessing) {
            auto posDauMcPart = posDaughterTrack.mcParticle();
            auto negDauMcPart = negDaughterTrack.mcParticle();
            if (v0mcparticle.pdgCode() != kK0Short || posDauMcPart.pdgCode() != kPiPlus || negDauMcPart.pdgCode() != kPiMinus)
              continue;
          }

          executeV0InCollisionloop(posDaughterTrack, negDaughterTrack, v0,
                                   posPiIdMethod, posPrIdMethod, negPiIdMethod, negPrIdMethod,
                                   v0Tag, isK0s, v0DauCollisionIndexTag, nK0s, centrality);
        } // End of V0s Loop

        nTrack = 0;
        nRejectedPiPlus = 0;
        nRejectedPiMinus = 0;
        for (const auto& track : tracksTablePerColl) {
          if (!track.has_mcParticle()) {
            LOG(warning) << "No MC Particle for this track, skip...";
            continue;
          }

          auto mcPart = track.mcParticle();
          if (!mcPart.isPhysicalPrimary()) {
            continue;
          }

          bool isAcceptedTrack = true;
          executeTrackQAPart(track, k0sFullDauList, rejectionTag, nRejectedPiMinus, nRejectedPiPlus, nTrack, isAcceptedTrack);
          if (!isAcceptedTrack) {
            continue;
          }

          // Do Proper Track Identification
          trackIsPion = false;
          trackIsKaon = false;
          trackIsProton = false;
          trackIsElectron = false;
          trackIsDeuteron = false;

          trackIdTag = 0;
          idMethodPi = kUnidentified;
          idMethodKa = kUnidentified;
          idMethodPr = kUnidentified;
          idMethodEl = kUnidentified;
          idMethodDe = kUnidentified;

          if (selPion(track, idMethodPi)) {
            trackIsPion = true;
            BITSET(trackIdTag, ID_BIT_PI);
          }
          if (selKaon(track, idMethodKa)) {
            trackIsKaon = true;
            BITSET(trackIdTag, ID_BIT_KA);
          }
          if (selProton(track, idMethodPr)) {
            trackIsProton = true;
            BITSET(trackIdTag, ID_BIT_PR);
          }
          if (selElectron(track, idMethodEl)) {
            trackIsElectron = true;
            BITSET(trackIdTag, ID_BIT_EL);
          }
          if (selDeuteron(track, idMethodDe)) {
            trackIsDeuteron = true;
            BITSET(trackIdTag, ID_BIT_DE);
          }

          if constexpr (analysisType == doPurityProcessing) {
            if (trackIsPion) {
              if (track.sign() > 0 && mcPart.pdgCode() != kPiPlus) {
                trackIsPion = false;
              }
              if (track.sign() < 0 && mcPart.pdgCode() != kPiMinus) {
                trackIsPion = false;
              }
            }
            if (trackIsKaon) {
              if (track.sign() > 0 && mcPart.pdgCode() != kKPlus) {
                trackIsKaon = false;
              }
              if (track.sign() < 0 && mcPart.pdgCode() != kKMinus) {
                trackIsKaon = false;
              }
            }
            if (trackIsProton) {
              if (track.sign() > 0 && mcPart.pdgCode() != kProton) {
                trackIsProton = false;
              }
              if (track.sign() < 0 && mcPart.pdgCode() != kProtonBar) {
                trackIsProton = false;
              }
            }
            if (trackIsElectron) {
              if (track.sign() > 0 && mcPart.pdgCode() != kPositron) {
                trackIsElectron = false;
              }
              if (track.sign() < 0 && mcPart.pdgCode() != kElectron) {
                trackIsElectron = false;
              }
            }
            if (trackIsDeuteron) {
              if (track.sign() > 0 && mcPart.pdgCode() != kDeuteron) {
                trackIsDeuteron = false;
              }
              if (track.sign() < 0 && mcPart.pdgCode() != -kDeuteron) {
                trackIsDeuteron = false;
              }
            }
          }

          executeTrackAnalysisPart(track, trackIdTag,
                                   idMethodPi, trackIsPion, nPiMinus, nPiPlus,
                                   idMethodKa, trackIsKaon, nKaMinus, nKaPlus,
                                   idMethodPr, trackIsProton, nProton, nPBar,
                                   idMethodEl, trackIsElectron, nElPlus, nElMinus,
                                   idMethodDe, trackIsDeuteron, nDePlus, nDeMinus);
        } // track loop ends

        executeSparseAnalysisPart(centrality, nTrack, nK0s,
                                  nRejectedPiPlus, nRejectedPiMinus, nKaon,
                                  nPiPlus, nKaPlus, nProton, nElPlus, nDePlus,
                                  nPiMinus, nKaMinus, nPBar, nElMinus, nDeMinus);

        executeEventInfoPart(collision, centrality, v0sTablePerColl.size(), tracksTablePerColl,
                             nTrack, nK0s,
                             nPiPlus, nKaPlus, nProton, nElPlus, nDePlus,
                             nPiMinus, nKaMinus, nPBar, nElMinus, nDeMinus);
      } // collision loop ends
    }
  } //

  //____________________________________Process Funtion For Analysis Starts Here____________________________________//

  void processData(MyCollisions const& collisions, MyV0s const& V0s, MyTracks const& tracks)
  {
    recoEvent.fill(HIST("recoEvent/ProcessType"), doDataProcessing);
    executeAnalysis<doDataProcessing>(collisions, V0s, tracks);

  } // Process Function Ends
  PROCESS_SWITCH(KaonIsospinFluctuations, processData, "Process for Data", true);

  void processReco(MyCollisionsWithMcLabels const& collisions, MyV0sWithMcLabels const& V0s, MyTracksWithMcLabels const& tracks, aod::McParticles const&)
  {
    recoEvent.fill(HIST("recoEvent/ProcessType"), doRecoProcessing);
    executeAnalysis<doRecoProcessing>(collisions, V0s, tracks);

  } // Process function is over
  PROCESS_SWITCH(KaonIsospinFluctuations, processReco, "Process for Reco", false);

  void processPurity(MyCollisionsWithMcLabels const& collisions, MyV0sWithMcLabels const& V0s, MyTracksWithMcLabels const& tracks, aod::McParticles const&)
  {
    recoEvent.fill(HIST("recoEvent/ProcessType"), doPurityProcessing);
    executeAnalysis<doPurityProcessing>(collisions, V0s, tracks);

  } // Process function is over
  PROCESS_SWITCH(KaonIsospinFluctuations, processPurity, "Process for Purity", false);

  Preslice<aod::McParticles> mcTracksPerMcCollisionPreslice = o2::aod::mcparticle::mcCollisionId;

  using MyMcCollisions = aod::McCollisions;
  void processGen(MyMcCollisions const&, MyCollisionsWithMcLabels const& collisions, aod::McParticles const& mcParticles)
  {
    recoEvent.fill(HIST("recoEvent/ProcessType"), doGenProcessing);
    float centrality = -1;
    for (const auto& collision : collisions) {
      if (!collision.has_mcCollision()) {
        continue;
      }
      centrality = -1;
      const auto& mcColl = collision.mcCollision();

      centrality = collision.centFT0C();
      if (cfgCentAxis.centAxisType == 1) {
        centrality = collision.centFT0M();
      } else if (cfgCentAxis.centAxisType == 2) {
        centrality = collision.multFT0M();
      } else if (cfgCentAxis.centAxisType == 3) {
        centrality = collision.multFT0C();
      }

      // group over mcParticles
      const auto mcTracksTablePerMcColl = mcParticles.sliceBy(mcTracksPerMcCollisionPreslice, mcColl.globalIndex());

      float nRejectedPiPlus = 0;
      float nRejectedPiMinus = 0;

      float nK0s = 0;
      float nPiPlus = 0;
      float nPiMinus = 0;
      float nKaPlus = 0;
      float nKaMinus = 0;
      float nProton = 0;
      float nPBar = 0;
      float nElPlus = 0;
      float nElMinus = 0;
      float nDePlus = 0;
      float nDeMinus = 0;
      float nTrack = 0;
      float nKaon = 0;

      for (const auto& mcTrack : mcTracksTablePerMcColl) {
        if (!mcTrack.isPhysicalPrimary()) {
          continue;
        }

        if (mcTrack.pdgCode() == kK0Short &&
            k0sSelCut.cfgK0sLowPt < mcTrack.pt() && mcTrack.pt() < k0sSelCut.cfgK0sHighPt &&
            std::abs(mcTrack.y()) < k0sSelCut.cfgK0sRapitidy) {
          nK0s++;
          fillGenTrackQA<genAnalysisDir, kK0s>(genAnalysis, mcTrack);
        }

        if (mcTrack.pt() <= cfgTrackCuts.cfgTrackPtLow || mcTrack.pt() >= cfgTrackCuts.cfgTrackPtHigh || std::abs(mcTrack.eta()) >= cfgTrackCuts.cfgTrackEta) {
          continue;
        }

        if (mcTrack.pdgCode() == kPiPlus) {
          fillGenTrackQA<genAnalysisDir, kPi>(genAnalysis, mcTrack);
          nPiPlus++;
        } else if (mcTrack.pdgCode() == kPiMinus) {
          fillGenTrackQA<genAnalysisDir, kPi>(genAnalysis, mcTrack);
          nPiMinus++;
        } else if (mcTrack.pdgCode() == kKPlus) {
          fillGenTrackQA<genAnalysisDir, kKa>(genAnalysis, mcTrack);
          nKaPlus++;
        } else if (mcTrack.pdgCode() == kKMinus) {
          fillGenTrackQA<genAnalysisDir, kKa>(genAnalysis, mcTrack);
          nKaMinus++;
        } else if (mcTrack.pdgCode() == kProton) {
          fillGenTrackQA<genAnalysisDir, kPr>(genAnalysis, mcTrack);
          nProton++;
        } else if (mcTrack.pdgCode() == kProtonBar) {
          fillGenTrackQA<genAnalysisDir, kPr>(genAnalysis, mcTrack);
          nPBar++;
        } else if (mcTrack.pdgCode() == kElectron) {
          fillGenTrackQA<genAnalysisDir, kEl>(genAnalysis, mcTrack);
          nElPlus++;
        } else if (mcTrack.pdgCode() == kPositron) {
          fillGenTrackQA<genAnalysisDir, kEl>(genAnalysis, mcTrack);
          nElMinus++;
        } else if (mcTrack.pdgCode() == kDeuteron) {
          fillGenTrackQA<genAnalysisDir, kDe>(genAnalysis, mcTrack);
          nDePlus++;
        } else if (mcTrack.pdgCode() == -kDeuteron) {
          fillGenTrackQA<genAnalysisDir, kDe>(genAnalysis, mcTrack);
          nDeMinus++;
        }

        nTrack++;
      } // mcTrack loop is over
      nKaon = nKaPlus + nKaMinus;
      executeSparseAnalysisPart(centrality, nTrack, nK0s,
                                nRejectedPiPlus, nRejectedPiMinus, nKaon,
                                nPiPlus, nKaPlus, nProton, nElPlus, nDePlus,
                                nPiMinus, nKaMinus, nPBar, nElMinus, nDeMinus);

      executeEventInfoPart(mcColl, centrality, 0, mcTracksTablePerMcColl,
                           nTrack, nK0s,
                           nPiPlus, nKaPlus, nProton, nElPlus, nDePlus,
                           nPiMinus, nKaMinus, nPBar, nElMinus, nDeMinus);
    } // collision loop is over
  }
  PROCESS_SWITCH(KaonIsospinFluctuations, processGen, "Process for Gen", false);

  template <typename T>
  void getV0MCount(const T& mcTrack, float& multV0M)
  {
    if ((-3.7 < mcTrack.eta() && mcTrack.eta() < -1.7) || (2.8 < mcTrack.eta() && mcTrack.eta() < 5.1)) {
      if (doFWDPtDependentCheck) {
        if (mcTrack.pt() > cfgFWDPtCut) {
          multV0M++; // V0C: at -3.7 < η < -1.7 (backward direction).
                     // V0A: at 2.8 < η < 5.1 (forward direction).
        }
      } else {
        multV0M++;
      }
    }
  }

  void processSim(MyMcCollisions const& mcCollisions, aod::McParticles const& mcParticles)
  {
    auto finalParticleIdList = (std::vector<int>)cfgFinalParticleIdList;
    auto nonFinalParticleIdList = (std::vector<int>)cfgNonFinalParticleIdList;
    std::sort(finalParticleIdList.begin(), finalParticleIdList.end());
    std::sort(nonFinalParticleIdList.begin(), nonFinalParticleIdList.end());

    recoEvent.fill(HIST("recoEvent/ProcessType"), doSimProcessing);
    float centrality = -1;
    for (const auto& mcColl : mcCollisions) {
      centrality = -1;

      if (cfgVtxZCheck) {
        if (std::abs(mcColl.posZ()) >= cutZvertex) {
          continue;
        }
      }
      // group over mcParticles
      const auto mcTracksTablePerMcColl = mcParticles.sliceBy(mcTracksPerMcCollisionPreslice, mcColl.globalIndex());

      float nRejectedPiPlus = 0;
      float nRejectedPiMinus = 0;
      float nK0s = 0;
      float nPiPlus = 0;
      float nPiMinus = 0;
      float nKaPlus = 0;
      float nKaMinus = 0;
      float nProton = 0;
      float nPBar = 0;
      float nElPlus = 0;
      float nElMinus = 0;
      float nDePlus = 0;
      float nDeMinus = 0;
      float nTrack = 0;
      float nKaon = 0;

      float multV0M = 0;

      for (const auto& mcTrack : mcTracksTablePerMcColl) {

        if (cfgCountFinalParticles) {
          if (!mcTrack.has_daughters() && std::binary_search(finalParticleIdList.begin(), finalParticleIdList.end(), mcTrack.pdgCode())) {
            getV0MCount(mcTrack, multV0M);
          }
        }
        if (cfgCountNonFinalParticles) {
          if (mcTrack.has_daughters() && std::binary_search(nonFinalParticleIdList.begin(), nonFinalParticleIdList.end(), mcTrack.pdgCode())) {
            getV0MCount(mcTrack, multV0M);
          }
        }

        if (cfgCountPhysicalPrimAndFinalParticles) {
          if (!mcTrack.has_daughters() && mcTrack.isPhysicalPrimary()) {
            if (!(std::abs(mcTrack.pdgCode()) == kNuE || std::abs(mcTrack.pdgCode()) == kNuMu || std::abs(mcTrack.pdgCode()) == kNuTau)) {
              // Removed invisible neutrinos;
              getV0MCount(mcTrack, multV0M);
            }
          }
        }

        if (!mcTrack.isPhysicalPrimary()) {
          continue;
        }

        if (mcTrack.pdgCode() == kK0Short &&
            k0sSelCut.cfgK0sLowPt < mcTrack.pt() && mcTrack.pt() < k0sSelCut.cfgK0sHighPt &&
            std::abs(mcTrack.y()) < k0sSelCut.cfgK0sRapitidy) {
          nK0s++;
          fillGenTrackQA<genAnalysisDir, kK0s>(genAnalysis, mcTrack);
        }

        if (mcTrack.pt() <= cfgTrackCuts.cfgTrackPtLow || mcTrack.pt() >= cfgTrackCuts.cfgTrackPtHigh || std::abs(mcTrack.eta()) >= cfgTrackCuts.cfgTrackEta) {
          continue;
        }

        if (mcTrack.pdgCode() == kPiPlus) {
          fillGenTrackQA<genAnalysisDir, kPi>(genAnalysis, mcTrack);
          nPiPlus++;
        } else if (mcTrack.pdgCode() == kPiMinus) {
          fillGenTrackQA<genAnalysisDir, kPi>(genAnalysis, mcTrack);
          nPiMinus++;
        } else if (mcTrack.pdgCode() == kKPlus) {
          fillGenTrackQA<genAnalysisDir, kKa>(genAnalysis, mcTrack);
          nKaPlus++;
        } else if (mcTrack.pdgCode() == kKMinus) {
          fillGenTrackQA<genAnalysisDir, kKa>(genAnalysis, mcTrack);
          nKaMinus++;
        } else if (mcTrack.pdgCode() == kProton) {
          fillGenTrackQA<genAnalysisDir, kPr>(genAnalysis, mcTrack);
          nProton++;
        } else if (mcTrack.pdgCode() == kProtonBar) {
          fillGenTrackQA<genAnalysisDir, kPr>(genAnalysis, mcTrack);
          nPBar++;
        } else if (mcTrack.pdgCode() == kElectron) {
          fillGenTrackQA<genAnalysisDir, kEl>(genAnalysis, mcTrack);
          nElPlus++;
        } else if (mcTrack.pdgCode() == kPositron) {
          fillGenTrackQA<genAnalysisDir, kEl>(genAnalysis, mcTrack);
          nElMinus++;
        } else if (mcTrack.pdgCode() == kDeuteron) {
          fillGenTrackQA<genAnalysisDir, kDe>(genAnalysis, mcTrack);
          nDePlus++;
        } else if (mcTrack.pdgCode() == -kDeuteron) {
          fillGenTrackQA<genAnalysisDir, kDe>(genAnalysis, mcTrack);
          nDeMinus++;
        }

        nTrack++;
      } // mcTrack loop is over
      nKaon = nKaPlus + nKaMinus;
      centrality = multV0M;
      executeSparseAnalysisPart(centrality, nTrack, nK0s,
                                nRejectedPiPlus, nRejectedPiMinus, nKaon,
                                nPiPlus, nKaPlus, nProton, nElPlus, nDePlus,
                                nPiMinus, nKaMinus, nPBar, nElMinus, nDeMinus);

      executeEventInfoPart(mcColl, centrality, 0, mcTracksTablePerMcColl,
                           nTrack, nK0s,
                           nPiPlus, nKaPlus, nProton, nElPlus, nDePlus,
                           nPiMinus, nKaMinus, nPBar, nElMinus, nDeMinus);
    } // collision loop is over
  }
  PROCESS_SWITCH(KaonIsospinFluctuations, processSim, "Process for Sim", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<KaonIsospinFluctuations>(cfgc)};
}
