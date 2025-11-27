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

/// \file   upcPhotonuclearAnalysisJMG.cxx
/// \brief  Task for photonuclear UPC analysis for azimuthal correlation: selection, histograms and observables.
/// \author Josué Martínez García <josuem@cern.ch>

#include "PWGCF/Core/CorrelationContainer.h"
#include "PWGUD/Core/UPCPairCuts.h"
#include "PWGUD/Core/UPCTauCentralBarrelHelperRL.h"
#include "PWGUD/DataModel/UDTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/MathConstants.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"

#include <TTree.h>

#include <algorithm>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace constants::math;
namespace o2::aod
{
namespace tree
{
DECLARE_SOA_COLUMN(PtSideA, ptSideA, std::vector<float>);
DECLARE_SOA_COLUMN(RapSideA, rapSideA, std::vector<float>);
DECLARE_SOA_COLUMN(PhiSideA, phiSideA, std::vector<float>);
DECLARE_SOA_COLUMN(TpcSignalSideA, tpcSignalSideA, std::vector<float>);
DECLARE_SOA_COLUMN(TofSignalSideA, tofSignalSideA, std::vector<float>);
DECLARE_SOA_COLUMN(TpcNSigmaPiSideA, tpcNSigmaPiSideA, std::vector<float>);
DECLARE_SOA_COLUMN(TofNSigmaPiSideA, tofNSigmaPiSideA, std::vector<float>);
DECLARE_SOA_COLUMN(TpcNSigmaKaSideA, tpcNSigmaKaSideA, std::vector<float>);
DECLARE_SOA_COLUMN(TofNSigmaKaSideA, tofNSigmaKaSideA, std::vector<float>);
DECLARE_SOA_COLUMN(PtSideC, ptSideC, std::vector<float>);
DECLARE_SOA_COLUMN(RapSideC, rapSideC, std::vector<float>);
DECLARE_SOA_COLUMN(PhiSideC, phiSideC, std::vector<float>);
DECLARE_SOA_COLUMN(TpcSignalSideC, tpcSignalSideC, std::vector<float>);
DECLARE_SOA_COLUMN(TofSignalSideC, tofSignalSideC, std::vector<float>);
DECLARE_SOA_COLUMN(TpcNSigmaPiSideC, tpcNSigmaPiSideC, std::vector<float>);
DECLARE_SOA_COLUMN(TofNSigmaPiSideC, tofNSigmaPiSideC, std::vector<float>);
DECLARE_SOA_COLUMN(TpcNSigmaKaSideC, tpcNSigmaKaSideC, std::vector<float>);
DECLARE_SOA_COLUMN(TofNSigmaKaSideC, tofNSigmaKaSideC, std::vector<float>);
DECLARE_SOA_COLUMN(NchSideA, nchSideA, int);
DECLARE_SOA_COLUMN(MultiplicitySideA, multiplicitySideA, int);
DECLARE_SOA_COLUMN(NchSideC, nchSideC, int);
DECLARE_SOA_COLUMN(MultiplicitySideC, multiplicitySideC, int);
} // namespace tree
DECLARE_SOA_TABLE(TREE, "AOD", "Tree",
                  tree::PtSideA,
                  tree::RapSideA,
                  tree::PhiSideA,
                  tree::TpcSignalSideA,
                  tree::TofSignalSideA,
                  tree::TpcNSigmaPiSideA,
                  tree::TofNSigmaPiSideA,
                  tree::TpcNSigmaKaSideA,
                  tree::TofNSigmaKaSideA,
                  tree::PtSideC,
                  tree::RapSideC,
                  tree::PhiSideC,
                  tree::TpcSignalSideC,
                  tree::TofSignalSideC,
                  tree::TpcNSigmaPiSideC,
                  tree::TofNSigmaPiSideC,
                  tree::TpcNSigmaKaSideC,
                  tree::TofNSigmaKaSideC,
                  tree::NchSideA,
                  tree::MultiplicitySideA,
                  tree::NchSideC,
                  tree::MultiplicitySideC);
} // namespace o2::aod

static constexpr float CFGPairCutDefaults[1][5] = {{-1, -1, -1, -1, -1}};
constexpr float kThreeHalfPi = 1.5f * PI;

struct UpcPhotonuclearAnalysisJMG {

  Produces<aod::TREE> tree;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Declare configurables on events/collisions
  Configurable<int> minMultiplicity{"minMultiplicity", 2, {"Range on multiplicity"}};
  Configurable<int> range1Max{"range1Max", 10, {"Range on multiplicity"}};
  Configurable<int> range2Min{"range2Min", 11, {"Range on multiplicity"}};
  Configurable<int> range2Max{"range2Max", 20, {"Range on multiplicity"}};
  Configurable<int> range3Min{"range3Min", 21, {"Range on multiplicity"}};
  Configurable<int> range3Max{"range3Max", 30, {"Range on multiplicity"}};
  Configurable<int> range4Min{"range4Min", 31, {"Range on multiplicity"}};
  Configurable<int> range4Max{"range4Max", 40, {"Range on multiplicity"}};
  Configurable<int> range5Min{"range5Min", 41, {"Range on multiplicity"}};
  Configurable<int> range5Max{"range5Max", 50, {"Range on multiplicity"}};
  Configurable<int> nEventsMixed{"nEventsMixed", 3, {"Events to be Mixed"}};
  Configurable<int> factorEventsMixed{"factorEventsMixed", 100, {"factorEventsMixed to events mixed"}};
  Configurable<float> myZVtxCut{"myZVtxCut", 10., {"My collision cut"}};
  Configurable<float> myTimeZNACut{"myTimeZNACut", 2., {"My collision cut"}};
  Configurable<float> myTimeZNCCut{"myTimeZNCCut", 2., {"My collision cut"}};
  // Declare configurables on side A gap
  Configurable<float> cutGapAMyEnergyZNA{"cutGapAMyEnergyZNA", 0., {"My collision cut. A Gap"}};
  // Configurable<float> cutAGapMyAmplitudeFT0AMax{"cutAGapMyAmplitudeFT0AMax", 200., {"My collision cut. A Gap"}};
  Configurable<float> cutGapAMyEnergyZNC{"cutGapAMyEnergyZNC", 1., {"My collision cut. A Gap"}};
  // Configurable<float> cutAGapMyAmplitudeFT0CMin{"cutAGapMyAmplitudeFT0CMin", 0., {"My collision cut. A Gap"}};
  // Declare configurables on side C gap
  Configurable<float> cutGapCMyEnergyZNA{"cutGapCMyEnergyZNA", 1., {"My collision cut. C Gap"}};
  // Configurable<float> cutCGapMyAmplitudeFT0AMin{"cutCGapMyAmplitudeFT0AMin", 0., {"My collision cut. A Gap"}};
  Configurable<float> cutGapCMyEnergyZNC{"cutGapCMyEnergyZNC", 0., {"My collision cut. C Gap"}};
  // Configurable<float> cutCGapMyAmplitudeFT0CMax{"cutCGapMyAmplitudeFT0CMax", 200., {"My collision cut. A Gap"}};
  // Declare configurables on tracks
  Configurable<float> cutMyptMin{"cutMyptMin", 0.2, {"My Track cut"}};
  Configurable<float> cutMyptMax{"cutMyptMax", 3., {"My Track cut"}};
  Configurable<float> cutMyetaMin{"cutMyetaMin", -0.8, {"My Track cut"}};
  Configurable<float> cutMyetaMax{"cutMyetaMax", 0.8, {"My Track cut"}};
  Configurable<float> cutMydcaZmax{"cutMydcaZmax", 2.f, {"My Track cut"}};
  Configurable<float> cutMydcaXYmax{"cutMydcaXYmax", 1e0f, {"My Track cut"}};
  Configurable<bool> cutMydcaXYusePt{"cutMydcaXYusePt", false, {"My Track cut"}};
  Configurable<bool> cutMyHasITS{"cutMyHasITS", true, {"My Track cut"}};
  Configurable<int> cutMyITSNClsMin{"cutMyITSNClsMin", 1, {"My Track cut"}};
  Configurable<float> cutMyITSChi2NClMax{"cutMyITSChi2NClMax", 36.f, {"My Track cut"}};
  Configurable<bool> cutMyHasTPC{"cutMyHasTPC", true, {"MyGlobalTrack cut"}};
  Configurable<int> cutMyTPCNClsCrossedRowsMin{"cutMyTPCNClsCrossedRowsMin", 70, {"My Track cut"}};
  Configurable<int> cutMyTPCNClsFindableMin{"cutMyTPCNClsFindableMin", 50, {"My Track cut"}};
  Configurable<int> cutMyTPCNClsMin{"cutMyTPCNClsMin", 1, {"My Track cut"}};
  Configurable<float> cutMyTPCNClsCrossedRowsOverNClsFindableMin{"cutMyTPCNClsCrossedRowsOverNClsFindableMin", 0.8f, {"My Track cut"}};
  Configurable<float> cutMyTPCNClsOverFindableNClsMin{"cutMyTPCNClsOverFindableNClsMin", 0.5f, {"My Track cut"}};
  Configurable<float> cutMyTPCChi2NclMax{"cutMyTPCChi2NclMax", 4.f, {"My Track cut"}};
  Configurable<float> myWeightMin{"myWeightMin", 0.2f, {"My Track cut"}};
  Configurable<float> myWeightMax{"myWeightMax", 5.f, {"My Track cut"}};
  Configurable<float> myEpsilonToWeight{"myEpsilonToWeight", 1e-6f, {"NUA correction"}};
  Configurable<bool> useEpsilon{"useEpsilon", false, {"NUA correction"}};
  Configurable<bool> useNMax{"useNMax", true, {"NUA correction"}};
  Configurable<LabeledArray<float>> cfgPairCut{"cfgPairCut",
                                               {CFGPairCutDefaults[0],
                                                5,
                                                {"Photon", "K0", "Lambda", "Phi", "Rho"}},
                                               "Pair cuts on various particles"};
  Configurable<float> cfgTwoTrackCut{"cfgTwoTrackCut", -1, {"Two track cut"}};
  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for histograms"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {32, -PIHalf, kThreeHalfPi}, "delta phi axis for histograms"};
  ConfigurableAxis axisDeltaEta{"axisDeltaEta", {32, -1.6, 1.6}, "delta eta axis for histograms"};
  ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt trigger axis for histograms"};
  ConfigurableAxis axisPtAssoc{"axisPtAssoc", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0}, "pt associated axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 110.1}, "multiplicity / multiplicity axis for histograms"};
  ConfigurableAxis axisVertexEfficiency{"axisVertexEfficiency", {10, -10, 10}, "vertex axis for efficiency histograms"};
  ConfigurableAxis axisEtaEfficiency{"axisEtaEfficiency", {20, -1.0, 1.0}, "eta axis for efficiency histograms"};
  ConfigurableAxis axisPtEfficiency{"axisPtEfficiency", {VARIABLE_WIDTH, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0}, "pt axis for efficiency histograms"};

  Filter collisionZVtxFilter = nabs(aod::collision::posZ) < myZVtxCut;
  Filter collisionZNTimeFilter = nabs(aod::udzdc::timeZNA) < myTimeZNACut && nabs(aod::udzdc::timeZNC) < myTimeZNCCut;
  Filter collisionZNeEnergyFilter = (aod::udzdc::energyCommonZNA < cutGapAMyEnergyZNA && aod::udzdc::energyCommonZNC >= cutGapAMyEnergyZNC) || (aod::udzdc::energyCommonZNA >= cutGapCMyEnergyZNA && aod::udzdc::energyCommonZNC < cutGapCMyEnergyZNC);
  Filter collisioSGFilter = aod::udcollision::gapSide == uint8_t(0) || aod::udcollision::gapSide == uint8_t(1);

  using FullSGUDCollision = soa::Filtered<soa::Join<aod::UDCollisions, aod::UDCollisionsSels, aod::SGCollisions, aod::UDZdcsReduced>>;
  using FullUDTracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksPID, aod::UDTracksDCA, aod::UDTracksFlags>;

  // Output definitions
  OutputObj<CorrelationContainer> same{"sameEvent"};
  OutputObj<CorrelationContainer> mixed{"mixedEvent"};
  // OutputObj<CorrelationContainer> sameGapSideA{"sameEventGapSideA"};
  // OutputObj<CorrelationContainer> mixedGapSideA{"mixedEventGapSideA"};
  // OutputObj<CorrelationContainer> sameGapSideC{"sameEventGapSideC"};
  // OutputObj<CorrelationContainer> mixedGapSideC{"mixedEventGapSideC"};

  UPCPairCuts mPairCuts;
  bool doPairCuts = false;

  void init(InitContext const&)
  {
    const AxisSpec axisCollision{4, -0.5, 3.5};
    const AxisSpec axisZvtx{20, -10., 10.};
    const AxisSpec axisPt{402, -0.05, 20.05};
    const AxisSpec axisP{402, -10.05, 10.05};
    const AxisSpec axisTPCSignal{802, -0.05, 400.05};
    const AxisSpec axisPhi{64, 0., TwoPI};
    const AxisSpec axisEta{32, -0.8, 0.8};
    const AxisSpec axisNch{601, -0.5, 600.5};
    const AxisSpec axisZNEnergy{1002, -0.5, 500.5};
    const AxisSpec axisZNTime{21, -10.5, 10.5};
    const AxisSpec axisFT0Amplitud{201, -0.5, 200.5};
    const AxisSpec axisNCls{201, -0.5, 200.5};
    const AxisSpec axisChi2NCls{100, 0, 50};
    const AxisSpec axisTPCNClsCrossedRowsMin{100, -0.05, 2.05};
    const AxisSpec axisCountTracks{17, -0.5, 16.5};

    histos.add("yields", "multiplicity vs pT vs eta", {HistType::kTH3F, {{100, 0, 100, "multiplicity"}, {40, 0, 20, "p_{T}"}, {100, -2, 2, "#eta"}}});
    histos.add("etaphi", "multiplicity vs eta vs phi", {HistType::kTH3F, {{100, 0, 100, "multiplicity"}, {100, -2, 2, "#eta"}, {64, 0., TwoPI, "#varphi"}}});
    histos.add("etaphiVtx", "vertex Z vs eta vs phi", {HistType::kTH3F, {{20, -10., 10., "vertex Z"}, {32, -0.8, 0.8, "#eta"}, {64, 0., TwoPI, "#varphi"}}});
    histos.add("sameEvent2D", "#Delta #eta vs #Delta #phi", {HistType::kTH2F, {axisDeltaEta, axisDeltaPhi}});
    histos.add("sameEvent_2_10", "#Delta #eta vs #Delta #phi", {HistType::kTH2F, {axisDeltaEta, axisDeltaPhi}});
    histos.add("sameEvent_11_20", "#Delta #eta vs #Delta #phi", {HistType::kTH2F, {axisDeltaEta, axisDeltaPhi}});
    histos.add("sameEvent_21_30", "#Delta #eta vs #Delta #phi", {HistType::kTH2F, {axisDeltaEta, axisDeltaPhi}});
    histos.add("sameEvent_31_40", "#Delta #eta vs #Delta #phi", {HistType::kTH2F, {axisDeltaEta, axisDeltaPhi}});
    histos.add("sameEvent_41_50", "#Delta #eta vs #Delta #phi", {HistType::kTH2F, {axisDeltaEta, axisDeltaPhi}});
    histos.add("mixedEvent2D", "#Delta #eta vs #Delta #phi", {HistType::kTH2F, {axisDeltaEta, axisDeltaPhi}});
    histos.add("mixedEvent_2_10", "#Delta #eta vs #Delta #phi", {HistType::kTH2F, {axisDeltaEta, axisDeltaPhi}});
    histos.add("mixedEvent_11_20", "#Delta #eta vs #Delta #phi", {HistType::kTH2F, {axisDeltaEta, axisDeltaPhi}});
    histos.add("mixedEvent_21_30", "#Delta #eta vs #Delta #phi", {HistType::kTH2F, {axisDeltaEta, axisDeltaPhi}});
    histos.add("mixedEvent_31_40", "#Delta #eta vs #Delta #phi", {HistType::kTH2F, {axisDeltaEta, axisDeltaPhi}});
    histos.add("mixedEvent_41_50", "#Delta #eta vs #Delta #phi", {HistType::kTH2F, {axisDeltaEta, axisDeltaPhi}});

    const int maxMixBin = axisMultiplicity->size() * axisVertex->size();
    histos.add("eventcount", "bin", {HistType::kTH1F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}}});
    mPairCuts.setHistogramRegistry(&histos);
    if (cfgPairCut->get("Photon") > 0 || cfgPairCut->get("K0") > 0 || cfgPairCut->get("Lambda") > 0 ||
        cfgPairCut->get("Phi") > 0 || cfgPairCut->get("Rho") > 0) {
      mPairCuts.setPairCut(UPCPairCuts::Photon, cfgPairCut->get("Photon"));
      mPairCuts.setPairCut(UPCPairCuts::K0, cfgPairCut->get("K0"));
      mPairCuts.setPairCut(UPCPairCuts::Lambda, cfgPairCut->get("Lambda"));
      mPairCuts.setPairCut(UPCPairCuts::Phi, cfgPairCut->get("Phi"));
      mPairCuts.setPairCut(UPCPairCuts::Rho, cfgPairCut->get("Rho"));
      doPairCuts = true;
    }
    histos.add("Events/hCountCollisions", "0 total - 1 side A - 2 side C - 3 both side; Number of analysed collision; counts", kTH1F, {axisCollision});
    histos.add("Events/hCountCollisionsMixed", "0 total - 1 side A - 2 side C - 3 both side; Number of analysed collision; counts", kTH1F, {axisCollision});
    histos.add("Tracks/hTracksAfterCuts", " ; ; counts", kTH1F, {axisCountTracks});

    // histos to selection gap in side A
    histos.add("Tracks/SGsideA/hTrackPt", "#it{p_{T}} distribution; #it{p_{T}}; counts", kTH1F, {axisPt});
    histos.add("Tracks/SGsideA/hTrackPhi", "#it{#phi} distribution; #it{#phi}; counts", kTH1F, {axisPhi});
    histos.add("Tracks/SGsideA/hTrackEta", "#it{#eta} distribution; #it{#eta}; counts", kTH1F, {axisEta});
    histos.add("Tracks/SGsideA/hTrackTPCSignnalP", "#it{TPC dE/dx vs p}; #it{p*charge}; #it{TPC dE/dx}", kTH2F, {axisP, axisTPCSignal});
    histos.add("Tracks/SGsideA/hTrackTOFSignnalP", "#it{TOF signal vs p}; #it{p*charge}; #it{TOF signal}", kTH2F, {axisP, axisTPCSignal});
    histos.add("Tracks/SGsideA/hTrackITSNCls", "#it{N Clusters ITS} distribution; #it{N Clusters ITS}; counts", kTH1F, {axisNCls});
    histos.add("Tracks/SGsideA/hTrackITSChi2NCls", "#it{N Clusters Chi2 ITS} distribution; #it{N Clusters Chi2 ITS}; counts", kTH1F, {axisChi2NCls});
    histos.add("Tracks/SGsideA/hTrackNClsCrossedRowsOverNClsFindable", "#it{NClsCrossedRows/FindableNCls} distribution in TPC; #it{NClsCrossedRows/FindableNCls}; counts", kTH1F, {axisTPCNClsCrossedRowsMin});
    histos.add("Tracks/SGsideA/hTrackNClsCrossedRowsOverNCls", "#it{NClsCrossedRows/NCls} distribution in TPC; #it{NClsCrossedRows/NCls}; counts", kTH1F, {axisTPCNClsCrossedRowsMin});
    histos.add("Tracks/SGsideA/hTrackTPCNClsCrossedRows", "#it{Number of crossed TPC Rows} distribution; #it{Number of crossed TPC Rows}; counts", kTH1F, {axisNCls});
    histos.add("Tracks/SGsideA/hTrackTPCNClsFindable", "#it{Findable TPC clusters per track} distribution; #it{Findable TPC clusters per track}; counts", kTH1F, {axisNCls});
    histos.add("Tracks/SGsideA/hTrackTPCNClsFound", "#it{Found TPC clusters per track} distribution; #it{Found TPC clusters per track}; counts", kTH1F, {axisNCls});
    histos.add("Tracks/SGsideA/hTrackTPCNClsFindableMinusFound", "#it{TPCNCls: Findable - Found per track} distribution; #it{TPCNCls: Findable - Found per track}; counts", kTH1F, {axisNCls});
    histos.add("Tracks/SGsideA/hTrackTPCNClsFindableMinusCrossedRows", "#it{TPCNCls: Findable - CrossedRows per track} distribution; #it{TPCNCls: Findable - CrossedRows per track}; counts", kTH1F, {axisNCls});
    histos.add("Tracks/SGsideA/hTrackTPCChi2NCls", "#it{N Clusters Chi2 TPC} distribution; #it{N Clusters Chi2 TPC}; counts", kTH1F, {axisChi2NCls});
    histos.add("Tracks/SGsideA/hTrackITSNClsTPCCls", "#it{ITS Clusters vs TPC Clusters}; #it{TPC Clusters}; #it{ITS Clusters}", kTH2F, {axisNCls, axisNCls});

    histos.add("Events/SGsideA/hZVtx", "vertex in z; z (cm); counts", kTH1F, {axisZvtx});
    histos.add("Events/SGsideA/hNch", "#it{Charged Tracks Multiplicity} distribution; #it{Charged Tracks Multiplicity}; counts", kTH1F, {axisNch});
    histos.add("Events/SGsideA/hMultiplicity", "#it{Multiplicity} distribution; #it{Multiplicity}; counts", kTH1F, {axisNch});
    histos.add("Events/SGsideA/hPtVSNch", "#it{ #LT p_{T} #GT } vs #it{Charged Tracks Multiplicity}; #it{Charged Tracks Multiplicity}; #it{ #LT p_{T} #GT }", kTH2F, {axisNch, axisPt});
    histos.add("Events/SGsideA/hEnergyZNA", "Energy in side A distribution; Energy in side A; counts", kTH1F, {axisZNEnergy});
    histos.add("Events/SGsideA/hEnergyZNC", "Energy in side C distribution; Energy in side C; counts", kTH1F, {axisZNEnergy});
    histos.add("Events/SGsideA/hEnergyRelationSides", "Energy in side A vs energy in side C; Energy in side A; Energy in side C", kTH2F, {axisZNEnergy, axisZNEnergy});
    histos.add("Events/SGsideA/hTimeZNA", "Time in side A distribution; Time in side A; counts", kTH1F, {axisZNTime});
    histos.add("Events/SGsideA/hTimeZNC", "Time in side C distribution; Time in side C; counts", kTH1F, {axisZNTime});
    histos.add("Events/SGsideA/hTimeRelationSides", "Time in side A vs time in side C; Time in side A; Time in side C", kTH2F, {axisZNTime, axisZNTime});
    histos.add("Events/SGsideA/hAmplitudFT0A", "Amplitud in side A distribution; Amplitud in side A; counts", kTH1F, {axisFT0Amplitud});
    histos.add("Events/SGsideA/hAmplitudFT0C", "Amplitud in side C distribution; Amplitud in side C; counts", kTH1F, {axisFT0Amplitud});
    histos.add("Events/SGsideA/hTrackPV", "#it{Nch vs Nch(PV)}; #it{N_{ch}(PV)}; #it{N_{ch}}", kTH2F, {axisNch, axisNch});

    // histos to selection gap in side C
    histos.add("Tracks/SGsideC/hTrackPt", "#it{p_{T}} distribution; #it{p_{T}}; counts", kTH1F, {axisPt});
    histos.add("Tracks/SGsideC/hTrackPhi", "#it{#phi} distribution; #it{#phi}; counts", kTH1F, {axisPhi});
    histos.add("Tracks/SGsideC/hTrackEta", "#it{#eta} distribution; #it{#eta}; counts", kTH1F, {axisEta});
    histos.add("Tracks/SGsideC/hTrackTPCSignnalP", "#it{TPC dE/dx vs p}; #it{p*charge}; #it{TPC dE/dx}", kTH2F, {axisP, axisTPCSignal});
    histos.add("Tracks/SGsideC/hTrackTOFSignnalP", "#it{TOF signal vs p}; #it{p*charge}; #it{TOF signal}", kTH2F, {axisP, axisTPCSignal});
    histos.add("Tracks/SGsideC/hTrackITSNCls", "#it{N Clusters ITS} distribution; #it{N Clusters ITS}; counts", kTH1F, {axisNCls});
    histos.add("Tracks/SGsideC/hTrackITSChi2NCls", "#it{N Clusters Chi2 ITS} distribution; #it{N Clusters Chi2 ITS}; counts", kTH1F, {axisChi2NCls});
    histos.add("Tracks/SGsideC/hTrackNClsCrossedRowsOverNClsFindable", "#it{NClsCrossedRows/FindableNCls} distribution in TPC; #it{NClsCrossedRows/FindableNCls}; counts", kTH1F, {axisTPCNClsCrossedRowsMin});
    histos.add("Tracks/SGsideC/hTrackNClsCrossedRowsOverNCls", "#it{NClsCrossedRows/NCls} distribution in TPC; #it{NClsCrossedRows/NCls}; counts", kTH1F, {axisTPCNClsCrossedRowsMin});
    histos.add("Tracks/SGsideC/hTrackTPCNClsCrossedRows", "#it{Number of crossed TPC Rows} distribution; #it{Number of crossed TPC Rows}; counts", kTH1F, {axisNCls});
    histos.add("Tracks/SGsideC/hTrackTPCNClsFindable", "#it{Findable TPC clusters per track} distribution; #it{Findable TPC clusters per track}; counts", kTH1F, {axisNCls});
    histos.add("Tracks/SGsideC/hTrackTPCNClsFound", "#it{Found TPC clusters per track} distribution; #it{Found TPC clusters per track}; counts", kTH1F, {axisNCls});
    histos.add("Tracks/SGsideC/hTrackTPCNClsFindableMinusFound", "#it{TPCNCls: Findable - Found per track} distribution; #it{TPCNCls: Findable - Found per track}; counts", kTH1F, {axisNCls});
    histos.add("Tracks/SGsideC/hTrackTPCNClsFindableMinusCrossedRows", "#it{TPCNCls: Findable - CrossedRows per track} distribution; #it{TPCNCls: Findable - CrossedRows per track}; counts", kTH1F, {axisNCls});
    histos.add("Tracks/SGsideC/hTrackTPCChi2NCls", "#it{N Clusters Chi2 TPC} distribution; #it{N Clusters Chi2 TPC}; counts", kTH1F, {axisChi2NCls});
    histos.add("Tracks/SGsideC/hTrackITSNClsTPCCls", "#it{ITS Clusters vs TPC Clusters}; #it{TPC Clusters}; #it{ITS Clusters}", kTH2F, {axisNCls, axisNCls});

    histos.add("Events/SGsideC/hZVtx", "vertex in z; z (cm); counts", kTH1F, {axisZvtx});
    histos.add("Events/SGsideC/hNch", "#it{Charged Tracks Multiplicity} distribution; #it{Charged Tracks Multiplicity}; counts", kTH1F, {axisNch});
    histos.add("Events/SGsideC/hMultiplicity", "#it{Multiplicity} distribution; #it{Multiplicity}; counts", kTH1F, {axisNch});
    histos.add("Events/SGsideC/hPtVSNch", "#it{ #LT p_{T} #GT } vs #it{Charged Tracks Multiplicity}; #it{Charged Tracks Multiplicity}; #it{ #LT p_{T} #GT }", kTH2F, {axisNch, axisPt});
    histos.add("Events/SGsideC/hEnergyZNA", "Energy in side A distribution; Energy in side A; counts", kTH1F, {axisZNEnergy});
    histos.add("Events/SGsideC/hEnergyZNC", "Energy in side C distribution; Energy in side C; counts", kTH1F, {axisZNEnergy});
    histos.add("Events/SGsideC/hEnergyRelationSides", "Energy in side A vs energy in side C; Energy in side A; Energy in side C", kTH2F, {axisZNEnergy, axisZNEnergy});
    histos.add("Events/SGsideC/hTimeZNA", "Time in side A distribution; Time in side A; counts", kTH1F, {axisZNTime});
    histos.add("Events/SGsideC/hTimeZNC", "Time in side C distribution; Time in side C; counts", kTH1F, {axisZNTime});
    histos.add("Events/SGsideC/hTimeRelationSides", "Time in side A vs time in side C; Time in side A; Time in side C", kTH2F, {axisZNTime, axisZNTime});
    histos.add("Events/SGsideC/hAmplitudFT0A", "Amplitud in side A distribution; Amplitud in side A; counts", kTH1F, {axisFT0Amplitud});
    histos.add("Events/SGsideC/hAmplitudFT0C", "Amplitud in side C distribution; Amplitud in side C; counts", kTH1F, {axisFT0Amplitud});
    histos.add("Events/SGsideC/hTrackPV", "#it{Nch vs Nch(PV)}; #it{N_{ch}(PV)}; #it{N_{ch}}", kTH2F, {axisNch, axisNch});

    std::vector<AxisSpec> corrAxis = {{axisDeltaEta, "#Delta#eta"},
                                      {axisPtAssoc, "p_{T} (GeV/c)"},
                                      {axisPtTrigger, "p_{T} (GeV/c)"},
                                      {axisMultiplicity, "multiplicity / multiplicity"},
                                      {axisDeltaPhi, "#Delta#varphi (rad)"},
                                      {axisVertex, "z-vtx (cm)"}};
    std::vector<AxisSpec> effAxis = {{axisEtaEfficiency, "#eta"},
                                     {axisEtaEfficiency, "#eta"},
                                     {axisPtEfficiency, "p_{T} (GeV/c)"},
                                     {axisVertexEfficiency, "z-vtx (cm)"}};
    same.setObject(new CorrelationContainer("sameEvent", "sameEvent", corrAxis, effAxis, {}));
    mixed.setObject(new CorrelationContainer("mixedEvent", "mixedEvent", corrAxis, effAxis, {}));
    // sameGapSideA.setObject(new CorrelationContainer("sameEventGapSideA", "sameEventGapSideA", corrAxis, effAxis, {}));
    // mixedGapSideA.setObject(new CorrelationContainer("mixedEventGapSideA", "mixedEventGapSideA", corrAxis, effAxis, {}));
    // sameGapSideC.setObject(new CorrelationContainer("sameEventGapSideC", "sameEventGapSideC", corrAxis, effAxis, {}));
    // mixedGapSideC.setObject(new CorrelationContainer("mixedEventGapSideC", "mixedEventGapSideC", corrAxis, effAxis, {}));
  }

  std::vector<double> vtxBinsEdges{VARIABLE_WIDTH, -10.0f, -7.0f, -5.0f, -2.5f, 0.0f, 2.5f, 5.0f, 7.0f, 10.0f};
  std::vector<double> gapSideBinsEdges{VARIABLE_WIDTH, -0.5, 0.5, 1.5};

  struct SameEventTag {
  };
  struct MixedEventTag {
  };

  SliceCache cache;
  // int countEvents = 0;
  // int countGapA = 0;
  // int countGapC = 0;

  // Binning only on PosZ without multiplicity
  using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::udcollision::GapSide>;
  // BinningType bindingOnVtx{{vtxBinsEdges, gapSideBinsEdges}, true};
  // using BinningType = ColumnBinningPolicy<aod::collision::PosZ>;
  // BinningType bindingOnVtx{{vtxBinsEdges}, true};
  // SameKindPair<FullSGUDCollision, FullUDTracks, BinningType> pairs{bindingOnVtx, nEventsMixed, -1, &cache};

  template <typename CSG>
  bool isCollisionCutSG(CSG const& collision, int SideGap)
  {
    bool gapSideA = (collision.energyCommonZNA() < cutGapAMyEnergyZNA) && (collision.energyCommonZNC() >= cutGapAMyEnergyZNC);
    bool gapSideC = (collision.energyCommonZNA() >= cutGapCMyEnergyZNA) && (collision.energyCommonZNC() < cutGapCMyEnergyZNC);

    switch (SideGap) {
      case 0:            // Gap in A side
        return gapSideA; // 0n - A side && Xn - C Side
        // if ((collision.totalFT0AmplitudeA() < cutAGapMyAmplitudeFT0AMax && collision.totalFT0AmplitudeC() >= cutAGapMyAmplitudeFT0CMin) == false) {
        //   return false;
        // }
        break;
      case 1:            // Gap in C side
        return gapSideC; // Xn - A side && 0n - C Side
        // if ((collision.totalFT0AmplitudeA() >= cutCGapMyAmplitudeFT0AMin && collision.totalFT0AmplitudeC() < cutCGapMyAmplitudeFT0CMax) == false) {
        //   return false;
        // }
        break;
      default:
        return false;
        break;
    }
  }

  template <typename CSG>
  bool isCollisionCutSG(CSG const& collision)
  {
    return isCollisionCutSG(collision, 0) || isCollisionCutSG(collision, 1);
  }

  template <typename T>
  bool isTrackCut(T const& track)
  {
    if (track.sign() != 1 && track.sign() != -1) {
      return false;
    }
    if (track.pt() < cutMyptMin || track.pt() > cutMyptMax) {
      return false;
    }
    if (eta(track.px(), track.py(), track.pz()) < cutMyetaMin || eta(track.px(), track.py(), track.pz()) > cutMyetaMax) {
      return false;
    }
    if (std::abs(track.dcaZ()) > cutMydcaZmax) {
      return false;
    }
    if (cutMydcaXYusePt) {
      float maxDCA = 0.0105f + 0.0350f / std::pow(track.pt(), 1.1f);
      if (std::abs(track.dcaXY()) > maxDCA) {
        return false;
      }
    } else {
      if (std::abs(track.dcaXY()) > cutMydcaXYmax) {
        return false;
      }
    }
    if (track.isPVContributor() == false) {
      return false;
    }
    // Quality Track
    // ITS
    if (cutMyHasITS && !track.hasITS()) {
      return false; // ITS refit
    }
    if (track.itsNCls() < cutMyITSNClsMin) {
      return false;
    }
    if (track.itsChi2NCl() > cutMyITSChi2NClMax) {
      return false;
    }
    // TPC
    if (cutMyHasTPC && !track.hasTPC()) {
      return false; // TPC refit
    }
    if (track.tpcNClsCrossedRows() < cutMyTPCNClsCrossedRowsMin) {
      return false;
    }
    if ((track.tpcNClsFindable() - track.tpcNClsFindableMinusFound()) < cutMyTPCNClsMin) {
      return false; // tpcNClsFound()
    }
    if (track.tpcNClsFindable() < cutMyTPCNClsFindableMin) {
      return false;
    }
    if ((static_cast<float>(track.tpcNClsCrossedRows()) / static_cast<float>(track.tpcNClsFindable())) < cutMyTPCNClsCrossedRowsOverNClsFindableMin) {
      return false; //
    }
    if ((static_cast<float>(track.tpcNClsFindable() - track.tpcNClsFindableMinusFound()) / static_cast<float>(track.tpcNClsFindable())) < cutMyTPCNClsCrossedRowsOverNClsFindableMin) {
      return false; //
    }
    if (track.tpcChi2NCl() > cutMyTPCChi2NclMax) {
      return false; // TPC chi2
    }
    return true;
  }

  template <typename TTracks>
  void fillQAUD(const TTracks tracks, float multiplicity)
  {
    for (const auto& track : tracks) {
      if (isTrackCut(track) == false) {
        continue;
      }
      float phiVal = RecoDecay::constrainAngle(phi(track.px(), track.py()), 0.f);
      histos.fill(HIST("yields"), multiplicity, track.pt(), eta(track.px(), track.py(), track.pz()));
      histos.fill(HIST("etaphi"), multiplicity, eta(track.px(), track.py(), track.pz()), phiVal);
    }
  }

  template <typename TTarget>
  bool fillCollisionUD(TTarget target, float multiplicity)
  {
    target->fillEvent(multiplicity, CorrelationContainer::kCFStepAll);
    target->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
    return true;
  }

  template <typename TTarget, typename TTracks, typename TTag>
  void fillCorrelationsUD(TTarget target, const TTracks& tracks1, const TTracks& tracks2, float multiplicity, float posZ, TTag)
  {
    for (const auto& track1 : tracks1) {
      if (isTrackCut(track1) == false) {
        return;
      }
      float phi1 = phi(track1.px(), track1.py());
      phi1 = RecoDecay::constrainAngle(phi1, 0.f);
      float eta1 = eta(track1.px(), track1.py(), track1.pz());
      target->getTriggerHist()->Fill(CorrelationContainer::kCFStepReconstructed, track1.pt(), multiplicity, posZ, 1.0);
      for (const auto& track2 : tracks2) {
        if (track1 == track2) {
          continue;
        }
        if (isTrackCut(track2) == false) {
          return;
        }
        float phi2 = phi(track2.px(), track2.py());
        phi2 = RecoDecay::constrainAngle(phi2, 0.f);
        float eta2 = eta(track2.px(), track2.py(), track2.pz());
        /*if (doPairCuts && mPairCuts.conversionCuts(track1, track2)) {
          continue;
        }*/
        float deltaPhi = phi1 - phi2;
        float deltaEta = eta1 - eta2;
        deltaPhi = RecoDecay::constrainAngle(deltaPhi, -PIHalf);
        target->getPairHist()->Fill(CorrelationContainer::kCFStepReconstructed,
                                    deltaEta,
                                    track2.pt(), track1.pt(),
                                    multiplicity,
                                    deltaPhi,
                                    posZ);
        if constexpr (std::is_same_v<TTag, SameEventTag>) {
          if (minMultiplicity <= multiplicity) {
            histos.fill(HIST("sameEvent2D"), deltaEta, deltaPhi);
          }
          if (minMultiplicity <= multiplicity && multiplicity <= range1Max) {
            histos.fill(HIST("sameEvent_2_10"), deltaEta, deltaPhi);
          }
          if (range2Min <= multiplicity && multiplicity <= range2Max) {
            histos.fill(HIST("sameEvent_11_20"), deltaEta, deltaPhi);
          }
          if (range3Min <= multiplicity && multiplicity <= range3Max) {
            histos.fill(HIST("sameEvent_21_30"), deltaEta, deltaPhi);
          }
          if (range4Min <= multiplicity && multiplicity <= range4Max) {
            histos.fill(HIST("sameEvent_31_40"), deltaEta, deltaPhi);
          }
          if (range5Min <= multiplicity && multiplicity <= range5Max) {
            histos.fill(HIST("sameEvent_41_50"), deltaEta, deltaPhi);
          }
        } else if constexpr (std::is_same_v<TTag, MixedEventTag>) {
          if (minMultiplicity <= multiplicity) {
            histos.fill(HIST("mixedEvent2D"), deltaEta, deltaPhi);
          }
          if (minMultiplicity <= multiplicity && multiplicity <= range1Max) {
            histos.fill(HIST("mixedEvent_2_10"), deltaEta, deltaPhi);
          }
          if (range2Min <= multiplicity && multiplicity <= range2Max) {
            histos.fill(HIST("mixedEvent_11_20"), deltaEta, deltaPhi);
          }
          if (range3Min <= multiplicity && multiplicity <= range3Max) {
            histos.fill(HIST("mixedEvent_21_30"), deltaEta, deltaPhi);
          }
          if (range4Min <= multiplicity && multiplicity <= range4Max) {
            histos.fill(HIST("mixedEvent_31_40"), deltaEta, deltaPhi);
          }
          if (range5Min <= multiplicity && multiplicity <= range5Max) {
            histos.fill(HIST("mixedEvent_41_50"), deltaEta, deltaPhi);
          }
        }
      }
    }
  }

  void processSG(FullSGUDCollision::iterator const& reconstructedCollision, FullUDTracks const& reconstructedTracks)
  {
    histos.fill(HIST("Events/hCountCollisions"), 0);
    int sgSide = reconstructedCollision.gapSide();
    int nTracksCharged = 0;
    float sumPt = 0;
    int nchPVGapSideA = 0;
    int nchPVGapSideC = 0;
    int nchGapSideA = 0;
    int nchGapSideC = 0;
    std::vector<float> vTrackPtSideA, vTrackEtaSideA, vTrackPhiSideA, vTrackTPCSignalSideA, vTrackTOFSignalSideA, vTrackTPCNSigmaPiSideA, vTrackTOFNSigmaPiSideA, vTrackTPCNSigmaKaSideA, vTrackTOFNSigmaKaSideA;
    std::vector<float> vTrackPtSideC, vTrackEtaSideC, vTrackPhiSideC, vTrackTPCSignalSideC, vTrackTOFSignalSideC, vTrackTPCNSigmaPiSideC, vTrackTOFNSigmaPiSideC, vTrackTPCNSigmaKaSideC, vTrackTOFNSigmaKaSideC;

    int nTracksChargedSideA(-222), nTracksChargedSideC(-222);
    int multiplicitySideA(-222), multiplicitySideC(-222);

    for (const auto& track : reconstructedTracks) {
      if (isTrackCut(track) == false) {
        continue;
      }
      float phiVal = RecoDecay::constrainAngle(phi(track.px(), track.py()), 0.f);
      histos.fill(HIST("etaphiVtx"), reconstructedCollision.posZ(), eta(track.px(), track.py(), track.pz()), phiVal);
    }

    switch (sgSide) {
      case 0: // gap for side A
        if (isCollisionCutSG(reconstructedCollision, 0) == false) {
          return;
        }
        histos.fill(HIST("Events/hCountCollisions"), 1);
        histos.fill(HIST("Events/SGsideA/hEnergyZNA"), reconstructedCollision.energyCommonZNA());
        histos.fill(HIST("Events/SGsideA/hEnergyZNC"), reconstructedCollision.energyCommonZNC());
        histos.fill(HIST("Events/SGsideA/hEnergyRelationSides"), reconstructedCollision.energyCommonZNA(), reconstructedCollision.energyCommonZNC());
        histos.fill(HIST("Events/SGsideA/hTimeZNA"), reconstructedCollision.timeZNA());
        histos.fill(HIST("Events/SGsideA/hTimeZNC"), reconstructedCollision.timeZNC());
        histos.fill(HIST("Events/SGsideA/hTimeRelationSides"), reconstructedCollision.timeZNA(), reconstructedCollision.timeZNC());
        histos.fill(HIST("Events/SGsideA/hZVtx"), reconstructedCollision.posZ());
        histos.fill(HIST("Events/SGsideA/hAmplitudFT0A"), reconstructedCollision.totalFT0AmplitudeA());
        histos.fill(HIST("Events/SGsideA/hAmplitudFT0C"), reconstructedCollision.totalFT0AmplitudeC());
        for (const auto& track : reconstructedTracks) {
          // LOGF(debug, "Filling tracks. Gap Side A");
          ++nchGapSideA;
          if (track.isPVContributor() == true) {
            ++nchPVGapSideA;
          }
          if (isTrackCut(track) == false) {
            continue;
          }
          nTracksCharged++;
          sumPt += track.pt();
          float phiVal = RecoDecay::constrainAngle(phi(track.px(), track.py()), 0.f);
          histos.fill(HIST("Tracks/SGsideA/hTrackPt"), track.pt());
          histos.fill(HIST("Tracks/SGsideA/hTrackPhi"), phiVal);
          histos.fill(HIST("Tracks/SGsideA/hTrackEta"), eta(track.px(), track.py(), track.pz()));
          histos.fill(HIST("Tracks/SGsideA/hTrackTPCSignnalP"), momentum(track.px(), track.py(), track.pz()) * track.sign(), track.tpcSignal());
          histos.fill(HIST("Tracks/SGsideA/hTrackTOFSignnalP"), momentum(track.px(), track.py(), track.pz()) * track.sign(), track.tofSignal());
          vTrackPtSideA.push_back(track.pt());
          vTrackEtaSideA.push_back(eta(track.px(), track.py(), track.pz()));
          vTrackPhiSideA.push_back(phiVal);
          vTrackTPCSignalSideA.push_back(track.tpcSignal());
          vTrackTOFSignalSideA.push_back(track.tofSignal());
          vTrackTPCNSigmaPiSideA.push_back(track.tpcNSigmaPi());
          vTrackTOFNSigmaPiSideA.push_back(track.tofNSigmaPi());
          vTrackTPCNSigmaKaSideA.push_back(track.tpcNSigmaKa());
          vTrackTOFNSigmaKaSideA.push_back(track.tofNSigmaKa());

          histos.fill(HIST("Tracks/SGsideA/hTrackITSNCls"), track.itsNCls());
          histos.fill(HIST("Tracks/SGsideA/hTrackITSChi2NCls"), track.itsChi2NCl());
          histos.fill(HIST("Tracks/SGsideA/hTrackNClsCrossedRowsOverNClsFindable"), (static_cast<float>(track.tpcNClsCrossedRows()) / static_cast<float>(track.tpcNClsFindable())));
          histos.fill(HIST("Tracks/SGsideA/hTrackNClsCrossedRowsOverNCls"), (static_cast<float>(track.tpcNClsFindable() - track.tpcNClsFindableMinusFound()) / static_cast<float>(track.tpcNClsFindable())));
          histos.fill(HIST("Tracks/SGsideA/hTrackTPCNClsCrossedRows"), track.tpcNClsCrossedRows());
          histos.fill(HIST("Tracks/SGsideA/hTrackTPCNClsFindable"), track.tpcNClsFindable());
          histos.fill(HIST("Tracks/SGsideA/hTrackTPCNClsFound"), track.tpcNClsFindable() - track.tpcNClsFindableMinusFound());
          histos.fill(HIST("Tracks/SGsideA/hTrackTPCNClsFindableMinusFound"), track.tpcNClsFindableMinusFound());
          histos.fill(HIST("Tracks/SGsideA/hTrackTPCNClsFindableMinusCrossedRows"), track.tpcNClsFindableMinusCrossedRows());
          histos.fill(HIST("Tracks/SGsideA/hTrackTPCChi2NCls"), track.tpcChi2NCl());
          histos.fill(HIST("Tracks/SGsideA/hTrackITSNClsTPCCls"), track.tpcNClsFindable() - track.tpcNClsFindableMinusFound(), track.itsNCls());
        }
        histos.fill(HIST("Events/SGsideA/hNch"), nTracksCharged);
        histos.fill(HIST("Events/SGsideA/hMultiplicity"), reconstructedTracks.size());
        histos.fill(HIST("Events/SGsideA/hPtVSNch"), nTracksCharged, (sumPt / nTracksCharged));
        histos.fill(HIST("Events/SGsideA/hTrackPV"), nchPVGapSideA, nchGapSideA);
        nTracksChargedSideA = nTracksCharged;
        multiplicitySideA = reconstructedTracks.size();
        nTracksCharged = sumPt = 0;
        break;
      case 1: // gap for side C
        if (isCollisionCutSG(reconstructedCollision, 1) == false) {
          return;
        }
        histos.fill(HIST("Events/hCountCollisions"), 2);
        histos.fill(HIST("Events/SGsideC/hEnergyZNA"), reconstructedCollision.energyCommonZNA());
        histos.fill(HIST("Events/SGsideC/hEnergyZNC"), reconstructedCollision.energyCommonZNC());
        histos.fill(HIST("Events/SGsideC/hEnergyRelationSides"), reconstructedCollision.energyCommonZNA(), reconstructedCollision.energyCommonZNC());
        histos.fill(HIST("Events/SGsideC/hTimeZNA"), reconstructedCollision.timeZNA());
        histos.fill(HIST("Events/SGsideC/hTimeZNC"), reconstructedCollision.timeZNC());
        histos.fill(HIST("Events/SGsideC/hTimeRelationSides"), reconstructedCollision.timeZNA(), reconstructedCollision.timeZNC());
        histos.fill(HIST("Events/SGsideC/hZVtx"), reconstructedCollision.posZ());
        histos.fill(HIST("Events/SGsideC/hAmplitudFT0A"), reconstructedCollision.totalFT0AmplitudeA());
        histos.fill(HIST("Events/SGsideC/hAmplitudFT0C"), reconstructedCollision.totalFT0AmplitudeC());
        for (const auto& track : reconstructedTracks) {
          ++nchGapSideC;
          if (track.isPVContributor() == true) {
            ++nchPVGapSideC;
          }
          if (isTrackCut(track) == false) {
            continue;
          }
          nTracksCharged++;
          sumPt += track.pt();
          float phiVal = RecoDecay::constrainAngle(phi(track.px(), track.py()), 0.f);
          histos.fill(HIST("Tracks/SGsideC/hTrackPt"), track.pt());
          histos.fill(HIST("Tracks/SGsideC/hTrackPhi"), phiVal);
          histos.fill(HIST("Tracks/SGsideC/hTrackEta"), eta(track.px(), track.py(), track.pz()));
          histos.fill(HIST("Tracks/SGsideC/hTrackTPCSignnalP"), momentum(track.px(), track.py(), track.pz()) * track.sign(), track.tpcSignal());
          histos.fill(HIST("Tracks/SGsideC/hTrackTOFSignnalP"), momentum(track.px(), track.py(), track.pz()) * track.sign(), track.tofSignal());
          vTrackPtSideC.push_back(track.pt());
          vTrackEtaSideC.push_back(eta(track.px(), track.py(), track.pz()));
          vTrackPhiSideC.push_back(phiVal);
          vTrackTPCSignalSideC.push_back(track.tpcSignal());
          vTrackTOFSignalSideC.push_back(track.tofSignal());
          vTrackTPCNSigmaPiSideC.push_back(track.tpcNSigmaPi());
          vTrackTOFNSigmaPiSideC.push_back(track.tofNSigmaPi());
          vTrackTPCNSigmaKaSideC.push_back(track.tpcNSigmaKa());
          vTrackTOFNSigmaKaSideC.push_back(track.tofNSigmaKa());

          histos.fill(HIST("Tracks/SGsideC/hTrackITSNCls"), track.itsNCls());
          histos.fill(HIST("Tracks/SGsideC/hTrackITSChi2NCls"), track.itsChi2NCl());
          histos.fill(HIST("Tracks/SGsideC/hTrackNClsCrossedRowsOverNClsFindable"), (static_cast<float>(track.tpcNClsCrossedRows()) / static_cast<float>(track.tpcNClsFindable())));
          histos.fill(HIST("Tracks/SGsideC/hTrackNClsCrossedRowsOverNCls"), (static_cast<float>(track.tpcNClsFindable() - track.tpcNClsFindableMinusFound()) / static_cast<float>(track.tpcNClsFindable())));
          histos.fill(HIST("Tracks/SGsideC/hTrackTPCNClsCrossedRows"), track.tpcNClsCrossedRows());
          histos.fill(HIST("Tracks/SGsideC/hTrackTPCNClsFindable"), track.tpcNClsFindable());
          histos.fill(HIST("Tracks/SGsideC/hTrackTPCNClsFound"), track.tpcNClsFindable() - track.tpcNClsFindableMinusFound());
          histos.fill(HIST("Tracks/SGsideC/hTrackTPCNClsFindableMinusFound"), track.tpcNClsFindableMinusFound());
          histos.fill(HIST("Tracks/SGsideC/hTrackTPCNClsFindableMinusCrossedRows"), track.tpcNClsFindableMinusCrossedRows());
          histos.fill(HIST("Tracks/SGsideC/hTrackTPCChi2NCls"), track.tpcChi2NCl());
          histos.fill(HIST("Tracks/SGsideC/hTrackITSNClsTPCCls"), track.tpcNClsFindable() - track.tpcNClsFindableMinusFound(), track.itsNCls());
        }
        histos.fill(HIST("Events/SGsideC/hNch"), nTracksCharged);
        histos.fill(HIST("Events/SGsideC/hMultiplicity"), reconstructedTracks.size());
        histos.fill(HIST("Events/SGsideC/hPtVSNch"), nTracksCharged, (sumPt / nTracksCharged));
        histos.fill(HIST("Events/SGsideC/hTrackPV"), nchPVGapSideC, nchGapSideC);
        nTracksChargedSideC = nTracksCharged;
        multiplicitySideC = reconstructedTracks.size();
        nTracksCharged = sumPt = 0;
        break;
      default:
        return;
        break;
    }
    tree(vTrackPtSideA, vTrackEtaSideA, vTrackPhiSideA, vTrackTPCSignalSideA, vTrackTOFSignalSideA, vTrackTPCNSigmaPiSideA, vTrackTOFNSigmaPiSideA, vTrackTPCNSigmaKaSideA, vTrackTOFNSigmaKaSideA, vTrackPtSideC, vTrackEtaSideC, vTrackPhiSideC, vTrackTPCSignalSideA, vTrackTOFSignalSideA, vTrackTPCNSigmaPiSideA, vTrackTOFNSigmaPiSideA, vTrackTPCNSigmaKaSideA, vTrackTOFNSigmaKaSideA, nTracksChargedSideA, multiplicitySideA, nTracksChargedSideC, multiplicitySideC);
    // nTracksChargedSideA = nTracksChargedSideC = multiplicitySideA = multiplicitySideC = 0;
  }

  PROCESS_SWITCH(UpcPhotonuclearAnalysisJMG, processSG, "Process in UD tables", true);

  void processMixed(FullSGUDCollision const& reconstructedCollision, FullUDTracks const& reconstructedTracks)
  {
    // (void)reconstructedCollision;
    // int sgSide = reconstructedCollision.gapSide();
    // int sgSide = 0;

    // int maxCount = 0;
    // int maxCountGapA = 0;
    // int maxCountGapC = 0;

    // if (auto histEventCount = histos.get<TH1>(HIST("eventcount"))) {
    //   int binA = histEventCount->GetXaxis()->FindBin(-2);  Gap A
    //   int binC = histEventCount->GetXaxis()->FindBin(-1);  Gap C

    //   maxCount = histEventCount->GetBinContent(binA) * factorEventsMixed;
    //   maxCountGapA = histEventCount->GetBinContent(binA) * factorEventsMixed;
    //   maxCountGapC = histEventCount->GetBinContent(binC) * factorEventsMixed;
    // }

    BinningType bindingOnVtx{{vtxBinsEdges, gapSideBinsEdges}, true};
    // BinningType bindingOnVtx{{vtxBinsEdges}, true};
    auto tracksTuple = std::make_tuple(reconstructedTracks);
    SameKindPair<FullSGUDCollision, FullUDTracks, BinningType> pairs{bindingOnVtx, nEventsMixed, -1, reconstructedCollision, tracksTuple, &cache};

    for (const auto& [collision1, tracks1, collision2, tracks2] : pairs) {
      if (collision1.size() == 0 || collision2.size() == 0) {
        // LOGF(info, "One or both collisions are empty.");
        continue;
      }

      // if (countGapA >= maxCountGapA && countGapC >= maxCountGapC) {
      //   break;
      // }

      float multiplicity = 0;

      histos.fill(HIST("Events/hCountCollisionsMixed"), 0);

      if (isCollisionCutSG(collision1) == false || isCollisionCutSG(collision2) == false) {
        continue;
      }
      histos.fill(HIST("Events/hCountCollisionsMixed"), 1);
      // ++countEvents;
      // LOGF(info, "In the pairs loop");
      for (const auto& track : tracks1) {
        if (isTrackCut(track) == false) {
          continue;
        }
        ++multiplicity;
      }
      // multiplicity = tracks1.size();
      if (fillCollisionUD(mixed, multiplicity) == false) {
        return;
      }
      histos.fill(HIST("Events/hCountCollisionsMixed"), 2);
      // histos.fill(HIST("eventcount"), bindingOnVtx.getBin({collision1.posZ()}));
      // histos.fill(HIST("eventcount"), bindingOnVtx.getBin({collision1.posZ(), collision1.gapSide()}));
      fillCorrelationsUD(mixed, tracks1, tracks2, multiplicity, collision1.posZ(), MixedEventTag{});
      // LOGF(info, "Filling mixed events");

      // if (collision1.gapSide() == 0 && collision2.gapSide() == 0) { gap on side A
      // if (isCollisionCutSG(collision1, 0) == false && isCollisionCutSG(collision2, 0) == false) {
      // continue;
      // }
      // std::cout << "Counts for Gap A: " << countGapA << " Maximum Count for Gap A " << maxCountGapA << std::endl;
      // ++countGapA;
      // LOGF(info, "In the pairs loop, gap side A");
      // multiplicity = tracks1.size();
      // if (fillCollisionUD(mixedGapSideA, multiplicity) == false) {
      // return;
      // }
      // histos.fill(HIST("eventcount"), bindingOnVtx.getBin({collision1.posZ()}));
      // histos.fill(HIST("eventcount"), bindingOnVtx.getBin({collision1.posZ(), collision1.gapSide()}));
      // fillCorrelationsUD(mixedGapSideA, tracks1, tracks2, multiplicity, collision1.posZ());
      // LOGF(info, "Filling mixedGapSideA events, Gap for side A");
      // }

      // if (collision1.gapSide() == 1 && collision2.gapSide() == 1) { gap on side C
      // if (isCollisionCutSG(collision1, 1) == false && isCollisionCutSG(collision2, 1) == false) {
      // continue;
      // }
      // std::cout << "Counts for Gap C: " << countGapC << " Maximum Count for Gap C" << maxCountGapC << std::endl;
      // ++countGapC;
      // LOGF(info, "In the pairs loop, gap side C");
      // multiplicity = tracks1.size();
      // if (fillCollisionUD(mixedGapSideC, multiplicity) == false) {
      // return;
      // }
      // fillCorrelationsUD(mixedGapSideC, tracks1, tracks2, multiplicity, collision1.posZ());
      // LOGF(info, "Filling mixedGapSideC events, Gap for side C");
      // } else {
      // continue;
      // }
    }
  }

  PROCESS_SWITCH(UpcPhotonuclearAnalysisJMG, processMixed, "Process mixed events", true);

  void processSame(FullSGUDCollision::iterator const& reconstructedCollision, FullUDTracks const& reconstructedTracks)
  {
    // int sgSide = reconstructedCollision.gapSide();
    float multiplicity = 0;

    if (isCollisionCutSG(reconstructedCollision) == false) {
      return;
    }

    // Configure track flow histogram labels
    auto hFlow = histos.get<TH1>(HIST("Tracks/hTracksAfterCuts"));
    hFlow->GetXaxis()->SetBinLabel(1, "All tracks");
    hFlow->GetXaxis()->SetBinLabel(2, "Track sign");
    hFlow->GetXaxis()->SetBinLabel(3, "p_{T} range");
    hFlow->GetXaxis()->SetBinLabel(4, "#eta range");
    hFlow->GetXaxis()->SetBinLabel(5, "dcaZ");
    hFlow->GetXaxis()->SetBinLabel(6, "dcaXY");
    hFlow->GetXaxis()->SetBinLabel(7, "PV contrib cut");
    hFlow->GetXaxis()->SetBinLabel(8, "has ITS cut");
    hFlow->GetXaxis()->SetBinLabel(9, "N clusters ITS cut");
    hFlow->GetXaxis()->SetBinLabel(10, "#chi^{2} N cluster ITS cut");
    hFlow->GetXaxis()->SetBinLabel(11, "has TPC cut");
    hFlow->GetXaxis()->SetBinLabel(12, "N clusters crossed row TPC cut");
    hFlow->GetXaxis()->SetBinLabel(13, "(N cluster findable - N cluster minus findable) TPC cut");
    hFlow->GetXaxis()->SetBinLabel(14, "N cluster findable TPC cut");
    hFlow->GetXaxis()->SetBinLabel(15, "(N cluster crossed row / N cluster findable) TPC cut");
    hFlow->GetXaxis()->SetBinLabel(16, "(N cluster findable - N cluster minus findable) / N cluster findable cut");
    hFlow->GetXaxis()->SetBinLabel(17, "#chi^{2} N cluster TPC cut");

    for (const auto& track : reconstructedTracks) {
      histos.fill(HIST("Tracks/hTracksAfterCuts"), 0);
      if (track.sign() != 1 && track.sign() != -1) {
        continue;
      }
      histos.fill(HIST("Tracks/hTracksAfterCuts"), 1);
      if (track.pt() < cutMyptMin || track.pt() > cutMyptMax) {
        continue;
      }
      histos.fill(HIST("Tracks/hTracksAfterCuts"), 2);
      if (eta(track.px(), track.py(), track.pz()) < cutMyetaMin || eta(track.px(), track.py(), track.pz()) > cutMyetaMax) {
        continue;
      }
      histos.fill(HIST("Tracks/hTracksAfterCuts"), 3);
      if (std::abs(track.dcaZ()) > cutMydcaZmax) {
        continue;
      }
      histos.fill(HIST("Tracks/hTracksAfterCuts"), 4);
      if (cutMydcaXYusePt) {
        float maxDCA = 0.0105f + 0.0350f / std::pow(track.pt(), 1.1f);
        if (std::abs(track.dcaXY()) > maxDCA) {
          continue;
        }
      } else {
        if (std::abs(track.dcaXY()) > cutMydcaXYmax) {
          continue;
        }
      }
      histos.fill(HIST("Tracks/hTracksAfterCuts"), 5);
      if (track.isPVContributor() == false) {
        continue;
      }
      histos.fill(HIST("Tracks/hTracksAfterCuts"), 6);
      // Quality Track
      // ITS
      if (cutMyHasITS && !track.hasITS()) {
        continue; // ITS refit
      }
      histos.fill(HIST("Tracks/hTracksAfterCuts"), 7);
      if (track.itsNCls() < cutMyITSNClsMin) {
        continue;
      }
      histos.fill(HIST("Tracks/hTracksAfterCuts"), 8);
      if (track.itsChi2NCl() > cutMyITSChi2NClMax) {
        continue;
      }
      histos.fill(HIST("Tracks/hTracksAfterCuts"), 9);
      // TPC
      if (cutMyHasTPC && !track.hasTPC()) {
        continue; // TPC refit
      }
      histos.fill(HIST("Tracks/hTracksAfterCuts"), 10);
      if (track.tpcNClsCrossedRows() < cutMyTPCNClsCrossedRowsMin) {
        continue;
      }
      histos.fill(HIST("Tracks/hTracksAfterCuts"), 11);
      if ((track.tpcNClsFindable() - track.tpcNClsFindableMinusFound()) < cutMyTPCNClsMin) {
        continue; // tpcNClsFound()
      }
      histos.fill(HIST("Tracks/hTracksAfterCuts"), 12);
      if (track.tpcNClsFindable() < cutMyTPCNClsFindableMin) {
        continue;
      }
      histos.fill(HIST("Tracks/hTracksAfterCuts"), 13);
      if ((static_cast<float>(track.tpcNClsCrossedRows()) / static_cast<float>(track.tpcNClsFindable())) < cutMyTPCNClsCrossedRowsOverNClsFindableMin) {
        continue; //
      }
      histos.fill(HIST("Tracks/hTracksAfterCuts"), 14);
      if ((static_cast<float>(track.tpcNClsFindable() - track.tpcNClsFindableMinusFound()) / static_cast<float>(track.tpcNClsFindable())) < cutMyTPCNClsCrossedRowsOverNClsFindableMin) {
        continue; //
      }
      histos.fill(HIST("Tracks/hTracksAfterCuts"), 15);
      if (track.tpcChi2NCl() > cutMyTPCChi2NclMax) {
        continue; // TPC chi2
      }
      histos.fill(HIST("Tracks/hTracksAfterCuts"), 16);

      if (isTrackCut(track) == false) {
        continue;
      }
      ++multiplicity;
    }
    // multiplicity = reconstructedTracks.size();
    if (fillCollisionUD(same, multiplicity) == false) {
      return;
    }
    // LOGF(debug, "Filling same events");
    histos.fill(HIST("eventcount"), -2);
    if (minMultiplicity <= multiplicity && multiplicity <= range1Max) {
      histos.fill(HIST("eventcount"), 1);
    } else if (range2Min <= multiplicity && multiplicity <= range2Max) {
      histos.fill(HIST("eventcount"), 2);
    } else if (range3Min <= multiplicity && multiplicity <= range3Max) {
      histos.fill(HIST("eventcount"), 3);
    } else if (range4Min <= multiplicity && multiplicity <= range4Max) {
      histos.fill(HIST("eventcount"), 4);
    } else if (range5Min <= multiplicity && multiplicity <= range5Max) {
      histos.fill(HIST("eventcount"), 5);
    }
    fillQAUD(reconstructedTracks, multiplicity);
    fillCorrelationsUD(same, reconstructedTracks, reconstructedTracks, multiplicity, reconstructedCollision.posZ(), SameEventTag{});

    /*switch (sgSide) {
      case 0: // gap for side A
        if (isCollisionCutSG(reconstructedCollision, 0) == false) {
          return;
        }
        multiplicity = reconstructedTracks.size();
        if (fillCollisionUD(sameGapSideA, multiplicity) == false) {
          return;
        }
        LOGF(info, "Filling sameGapSideA events");
        histos.fill(HIST("eventcount"), -2);
        fillQAUD(reconstructedTracks);
        fillCorrelationsUD(sameGapSideA, reconstructedTracks, reconstructedTracks, multiplicity, reconstructedCollision.posZ());
        break;
      case 1: // gap for side C
        if (isCollisionCutSG(reconstructedCollision, 1) == false) {
          return;
        }
        multiplicity = reconstructedTracks.size();
        if (fillCollisionUD(sameGapSideC, multiplicity) == false) {
          return;
        }
        histos.fill(HIST("eventcount"), -1);
        // LOGF(info, "Filling sameGapSideC events");
        fillCorrelationsUD(sameGapSideC, reconstructedTracks, reconstructedTracks, multiplicity, reconstructedCollision.posZ());
        break;
      default:
        return;
        break;
    }*/
  }

  PROCESS_SWITCH(UpcPhotonuclearAnalysisJMG, processSame, "Process same event", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<UpcPhotonuclearAnalysisJMG>(cfgc)};
}
