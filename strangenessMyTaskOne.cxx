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
/// \brief Step2 of the Strangeness tutorial
/// \author Nepeivoda Roman (roman.nepeivoda@cern.ch)
/// \author Chiara De Martin (chiara.de.martin@cern.ch)

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/EventSelection.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/DataModel/PIDResponse.h"

#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/Centrality.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;


#define bitset(var, nbit) ((var) |= (static_cast<uint64_t>(1) << static_cast<uint64_t>(nbit)))
#define bitcheck(var, nbit) ((var) & (static_cast<uint64_t>(1) << static_cast<uint64_t>(nbit)))


//修改部分

#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
using std::array;

using dauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;
using dauMCTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackMCIds, aod::DauTrackTPCPIDs>;
using v0Candidates = soa::Join<aod::V0CollRefs, aod::V0Cores, aod::V0Extras, aod::V0TOFPIDs, aod::V0TOFNSigmas>;
using v0MCCandidates = soa::Join<aod::V0CollRefs, aod::V0Cores, aod::V0MCCores, aod::V0Extras, aod::V0TOFPIDs, aod::V0TOFNSigmas, aod::V0MCMothers, aod::V0MCCollRefs>;




// STEP 0
// Starting point: loop over all V0s and fill invariant mass histogram
// STEP 1
// Apply selections on topological variables of V0s
// STEP 2
// Apply PID selections on V0 daughter tracks

struct strangeness_tutorial {
  // Histograms are defined with HistogramRegistry
  HistogramRegistry rEventSelection{"eventSelection", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rKzeroShort{"kzeroShort", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rLambda{"Lambda", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  //master analysis switches
  Configurable<bool> analyseK0Short{"analyseK0Short", true, "process K0Short-like candidates"};
  Configurable<bool> analyseLambda{"analyseLambda", true, "process Lambda-like candidates"};
  Configurable<bool> analyseAntiLambda{"analyseAntiLambda", true, "process AntiLambda-like candidates"};



  //Selection criteria:acceptance
  Configurable<float> rapidityCut{"rapidityCut", 0.5, "rapidity"};
  Configurable<float> daughterEtaCut{"daughterEtaCut", 0.8, "max eta for daughters"};

  // Configurable for histograms
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};

  // Configurable for event selection
  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range (cm)"};

  // Standard 5 topological criteria
  Configurable<double> v0setting_cospa{"v0setting_cospa", 0.995, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0
  Configurable<float> v0setting_dcav0dau{"v0setting_dcav0dau", 1, "DCA V0 Daughters"};
  Configurable<float> dcanegtopv{"dcanegtopv", 0.05, "min DCA Neg To PV (cm)"};
  Configurable<float> dcapostopv{"dcapostopv", .2, "min DCA Pos To PV (cm)"};
  Configurable<float> v0setting_radius{"v0setting_radius", 1.2, "v0radius"};
  Configurable<float> v0radiusMax{"v0radiusMax", 1.2, "minimum V0 radius (cm)"};

  Configurable<float> armPodCut{"armPodCut", 5.0f, "pT * (cut) > |alpha|, AP cut. Negative: no cut"};


  //Track quality
  Configurable<int> minTPCrows{"minTPCrows", 70, "minimum TPC crossed rows"};
  Configurable<int> minITSclusters{"minITSclusters", -1, "minimum ITS clusters"};
  Configurable<bool> skipTPConly{"skipTPConly", false, "skip V0s comprised of at least one TPC only prong"};
  Configurable<bool> requirePosITSonly{"requirePosITSonly", false, "require that positive track is ITSonly (overrides TPC quality)"};
  Configurable<bool> requireNegITSonly{"requireNegITSonly", false, "require that negative track is ITSonly (overrides TPC quality)"};


  //PID
  Configurable<float> TpcPidNsigmaCut{"TpcPidNsigmaCut", 5, "TpcPidNsigmaCut"};

  
  //for MC
  static constexpr float defaultLifetimeCuts[1][2] = {{30., 20.}};
  Configurable<LabeledArray<float>> lifetimecut{"lifetimecut", {defaultLifetimeCuts[0], 2, {"lifetimecutLambda", "lifetimecutK0S"}}, "lifetimecut"};


  // Configurable parameters for PID selection
  Configurable<float> NSigmaTPCPion{"NSigmaTPCPion", 5, "NSigmaTPCPion"};

  ConfigurableAxis axisAPAlpha{"axisAPAlpha", {220, -1.1f, 1.1f}, "V0 AP alpha"};
  ConfigurableAxis axisAPQt{"axisAPQt", {220, 0.0f, 0.5f}, "V0 AP alpha"};

  enum selection : uint64_t { selCosPA = 0,
                              selRadius,
                              selRadiusMax,
                              selDCANegToPV,
                              selDCAPosToPV,
                              selDCAV0Dau,
                              selK0ShortRapidity,
                              selLambdaRapidity,
                              selTPCPIDPositivePion,
                              selTPCPIDNegativePion,
                              selTPCPIDPositiveProton,
                              selTPCPIDNegativeProton,
                              selLambdaCTau,
                              selK0ShortCTau,
                              selK0ShortArmenteros,
                              selLambdaArmenteros,
                              selAntiLambdaArmenteros,
                              selPosGoodTPCTrack, // at least min # TPC rows
                              selNegGoodTPCTrack,// at least min # TPC rows
                              selPosItsOnly,
                              selNegItsOnly,  
                              selPosGoodITSTrack,
                              selNegGoodITSTrack,
                              selConsiderK0Short,    // for mc tagging
                              selConsiderLambda,     // for mc tagging
                              selConsiderAntiLambda, // for mc tagging
                              selPhysPrimK0Short,    // for mc tagging
                              selPhysPrimLambda,     // for mc tagging
                              selPhysPrimAntiLambda, // for mc tagging
  };

  uint64_t maskTopological;
  uint64_t maskTopoNoV0Radius;
  uint64_t maskTopoNoDCANegToPV;
  uint64_t maskTopoNoDCAPosToPV;
  uint64_t maskTopoNoCosPA;
  uint64_t maskTopoNoDCAV0Dau;
  uint64_t maskTrackProperties;

  uint64_t maskK0ShortSpecific;
  uint64_t maskLambdaSpecific;
  uint64_t maskAntiLambdaSpecific;

  
  uint64_t maskSelectionK0Short;
  uint64_t maskSelectionLambda;
  uint64_t maskSelectionAntiLambda;

  uint64_t secondaryMaskSelectionLambda;



  void init(InitContext const&)
  {
    // initialise bit masks
    maskTopological = (uint64_t(1) << selCosPA) | (uint64_t(1) << selRadius) | (uint64_t(1) << selDCANegToPV) | (uint64_t(1) << selDCAPosToPV) | (uint64_t(1) << selDCAV0Dau) | (uint64_t(1) << selRadiusMax);
    maskTopoNoV0Radius = (uint64_t(1) << selCosPA) | (uint64_t(1) << selDCANegToPV) | (uint64_t(1) << selDCAPosToPV) | (uint64_t(1) << selDCAV0Dau) | (uint64_t(1) << selRadiusMax);
    maskTopoNoDCANegToPV = (uint64_t(1) << selCosPA) | (uint64_t(1) << selRadius) | (uint64_t(1) << selDCAPosToPV) | (uint64_t(1) << selDCAV0Dau) | (uint64_t(1) << selRadiusMax);
    maskTopoNoDCAPosToPV = (uint64_t(1) << selCosPA) | (uint64_t(1) << selRadius) | (uint64_t(1) << selDCANegToPV) | (uint64_t(1) << selDCAV0Dau) | (uint64_t(1) << selRadiusMax);
    maskTopoNoCosPA = (uint64_t(1) << selRadius) | (uint64_t(1) << selDCANegToPV) | (uint64_t(1) << selDCAPosToPV) | (uint64_t(1) << selDCAV0Dau) | (uint64_t(1) << selRadiusMax);
    maskTopoNoDCAV0Dau = (uint64_t(1) << selCosPA) | (uint64_t(1) << selRadius) | (uint64_t(1) << selDCANegToPV) | (uint64_t(1) << selDCAPosToPV) | (uint64_t(1) << selRadiusMax);
    
    maskK0ShortSpecific = (uint64_t(1) << selK0ShortRapidity) | (uint64_t(1) << selK0ShortCTau) | (uint64_t(1) << selK0ShortArmenteros) |(uint64_t(1) << selConsiderK0Short);
    maskLambdaSpecific = (uint64_t(1) << selLambdaRapidity) | (uint64_t(1) << selLambdaCTau)| (uint64_t(1) << selConsiderLambda) ;
    maskAntiLambdaSpecific = (uint64_t(1) << selLambdaRapidity) | (uint64_t(1) << selLambdaCTau) | (uint64_t(1) << selConsiderAntiLambda);

    // ask for specific TPC/TOF PID selections
    


    // Axes
    AxisSpec K0ShortMassAxis = {200, 0.45f, 0.55f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec LambdaMassAxis = {200, 0.9f, 1.3f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec vertexZAxis = {nBins, -15., 15., "vrtx_{Z} [cm]"};
    AxisSpec ptAxis = {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};

    ConfigurableAxis axisTPCrows{"axisTPCrows", {160, 0.0f, 160.0f}, "N TPC rows"};
    ConfigurableAxis axisITSclus{"axisITSclus", {7, 0.0f, 7.0f}, "N ITS Clusters"};

    ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for analysis"};
    ConfigurableAxis axisLambdaMass{"axisLambdaMass", {200, 1.101f, 1.131f}, ""};
    // Histograms
    // Event selection
    rEventSelection.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});

    // K0s reconstruction
    // Mass
    rKzeroShort.add("hMassK0Short", "hMassK0Short", {HistType::kTH1F, {K0ShortMassAxis}});
    rKzeroShort.add("hMassK0ShortSelected", "hMassK0ShortSelected", {HistType::kTH1F, {K0ShortMassAxis}});

    // K0s topological/PID cuts
    rKzeroShort.add("hDCAV0Daughters", "hDCAV0Daughters", {HistType::kTH1F, {{55, 0.0f, 2.2f}}});
    rKzeroShort.add("hV0CosPA", "hV0CosPA", {HistType::kTH1F, {{100, 0.95f, 1.f}}});
    rKzeroShort.add("hNSigmaPosPionFromK0s", "hNSigmaPosPionFromK0s", {HistType::kTH2F, {{100, -5.f, 5.f}, {ptAxis}}});
    rKzeroShort.add("hNSigmaNegPionFromK0s", "hNSigmaNegPionFromK0s", {HistType::kTH2F, {{100, -5.f, 5.f}, {ptAxis}}});

    
    rLambda.add("hEventSelection", "hEventSelection", {HistType::kTH1F, {{5,-0.5f,+4.5f}}});
    rLambda.add("hLambda", "hLambda", {HistType::kTH1F, {LambdaMassAxis}});
    rLambda.add("hMassLambdaSelected", "hMassLambdaSelected", {HistType::kTH1F, {LambdaMassAxis}});
    rLambda.add("hDCAV0Daughters", "hDCAV0Daughters", {HistType::kTH1F, {{55, 0.0f, 2.2f}}});
    rLambda.add("hV0CosPA", "hV0CosPA", {HistType::kTH1F, {{100, 0.95f, 1.f}}});
    rLambda.add("hNSigmaPosProtonFromLambdas", "hNSigmaPosProtonFromLambdas", {HistType::kTH2F, {{100, -5.f, 5.f}, {ptAxis}}});
    rLambda.add("hNSigmaNegPionFromLambdas", "hNSigmaNegPionFromLambdas", {HistType::kTH2F, {{100, -5.f, 5.f}, {ptAxis}}});
    
    rLambda.add("hPosDCAToPV", "hPosDCAToPV", {HistType::kTH1F, {{1000, 0.f, 100.f}}});
    rLambda.add("hNegDCAToPV", "hNegDCAToPV", {HistType::kTH1F, {{1000, 0.f, 100.f}}});
    rLambda.add("hPostpcCrossedRows", "hPostpcCrossedRows", {HistType::kTH1F, {{100, 0.f, 160.f}}});
    rLambda.add("hNegtpcCrossedRows", "hNegtpcCrossedRows", {HistType::kTH1F, {{100, 0.f, 160.f}}});
    rLambda.add("hPostpcCrossedRowsOverFindableCls", "hPostpcCrossedRowsOverFindableCls", {HistType::kTH1F, {{100, 0.f, 1.2f}}});
    rLambda.add("hNegtpcCrossedRowsOverFindableCls", "hNegtpcCrossedRowsOverFindableCls", {HistType::kTH1F, {{100, 0.f, 1.2f}}});

    rLambda.add("hCentrality", "hCentrality", {HistType::kTH1F, {{100, 0.f, 100.f}}});
    rLambda.add("hposativeeta", "hposativeeta", {HistType::kTH1F, {{100, -1.f, 1.f}}});
    rLambda.add("hnegativeeta", "hnegativeeta", {HistType::kTH1F, {{100, -1.f, 1.f}}});
     
    rLambda.add("h2dArmenterosAll", "h2dArmenterosAll", kTH2F, {axisAPAlpha, axisAPQt});

    rLambda.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(1, "All collisions");
    rLambda.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(2, "sel8 cut");
    rLambda.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(3, "posZ cut");
    rLambda.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(4, "kNoITSROFrameBorder");
    rLambda.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(5, "kNoTimeFrameBorder");
    //rLambda.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(6, "kIsVertexITSTPC");
    //rLambda.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(7, "kIsGoodZvtxFT0vsPV");
    //rLambda.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(8, "kIsVertexTOFmatched");
    //rLambda.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(9, "kIsVertexTRDmatched");
    //rLambda.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(10, "kNoSameBunchPileup");

    
    histos.add("hPosDCAToPV", "hPosDCAToPV", {HistType::kTH1F, {{1000, 0.f, 100.f}}});
    histos.add("hNegDCAToPV", "hNegDCAToPV", {HistType::kTH1F, {{1000, 0.f, 100.f}}});
    histos.add("hDCADaughters", "hDCADaughters", {HistType::kTH1F, {{1000, 0.0f, 2.2f}}});
    histos.add("hPointingAngle", "hPointingAngle", {HistType::kTH1F, {{100, 0.95f, 1.f}}});
    histos.add("hV0Radius", "hV0Radius", {HistType::kTH1F, {{2000, 0.f, 100.f}}});
    histos.add("h2dPositiveITSvsTPCpts", "h2dPositiveITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
    histos.add("h2dNegativeITSvsTPCpts", "h2dNegativeITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
    histos.add("hPostpcCrossedRows", "hPostpcCrossedRows", {HistType::kTH1F, {{100, 0.f, 160.f}}});
    histos.add("hNegtpcCrossedRows", "hNegtpcCrossedRows", {HistType::kTH1F, {{100, 0.f, 160.f}}});
    histos.add("hMassLambda", "hMassLambda", kTH1F, {axisLambdaMass});
    histos.add("hMassAntiLambda", "hMassAntiLambda", kTH1F, {axisLambdaMass});
    histos.add("hKs0", "hKs0", kTH1F, {{200,0.45f,0.55f}});
    histos.add("h3dMassLambda", "h3dMassLambda", kTH2F, {axisPt, axisLambdaMass});
    histos.add("hpositiveeta", "hpositiveeta", {HistType::kTH1F, {{100, -1.f, 1.f}}});
    histos.add("hnegativeeta", "hnegativeeta", {HistType::kTH1F, {{100, -1.f, 1.f}}});
    histos.add("hPostpcCrossedRowsOverFindableCls", "hPostpcCrossedRowsOverFindableCls", {HistType::kTH1F, {{100, 0.f, 1.5f}}});
    histos.add("hNegtpcCrossedRowsOverFindableCls", "hNegtpcCrossedRowsOverFindableCls", {HistType::kTH1F, {{100, 0.f, 1.5f}}});
    histos.add("h2dArmenterosAll", "h2dArmenterosAll", kTH2F, {axisAPAlpha, axisAPQt});
    histos.add("h2dArmenterosLambda", "h2dArmenterosLambda", kTH2F, {axisAPAlpha, axisAPQt});
    histos.add("h2dArmenterosKshort", "h2dArmenterosKshort", kTH2F, {axisAPAlpha, axisAPQt});
    histos.add("h2dArmenterosAntiLambda", "h2dArmenterosAntiLambda", kTH2F, {axisAPAlpha, axisAPQt});
    histos.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});

    // ask for specific TPC/TOF PID selections
    maskTrackProperties = 0;
    
    if (requirePosITSonly) {
      maskTrackProperties = maskTrackProperties | (uint64_t(1) << selPosItsOnly) | (uint64_t(1) << selPosGoodITSTrack);
    } else {
      maskTrackProperties = maskTrackProperties | (uint64_t(1) << selPosGoodTPCTrack) | (uint64_t(1) << selPosGoodITSTrack);
      // TPC signal is available: ask for positive track PID
      if (TpcPidNsigmaCut < 1e+5) { // safeguard for no cut
        maskK0ShortSpecific = maskK0ShortSpecific | (uint64_t(1) << selTPCPIDPositivePion);
        maskLambdaSpecific = maskLambdaSpecific | (uint64_t(1) << selTPCPIDPositiveProton);
        maskAntiLambdaSpecific = maskAntiLambdaSpecific | (uint64_t(1) << selTPCPIDPositivePion);

      }
      // TOF PID
     
    }
    if (requireNegITSonly) {
      maskTrackProperties = maskTrackProperties | (uint64_t(1) << selNegItsOnly) | (uint64_t(1) << selNegGoodITSTrack);
    } else {
      maskTrackProperties = maskTrackProperties | (uint64_t(1) << selNegGoodTPCTrack) | (uint64_t(1) << selNegGoodITSTrack);
      // TPC signal is available: ask for negative track PID
      if (TpcPidNsigmaCut < 1e+5) { // safeguard for no cut
       maskK0ShortSpecific = maskK0ShortSpecific | (uint64_t(1) << selTPCPIDNegativePion);
       maskLambdaSpecific = maskLambdaSpecific | (uint64_t(1) << selTPCPIDNegativePion);
       maskAntiLambdaSpecific = maskAntiLambdaSpecific | (uint64_t(1) << selTPCPIDNegativeProton);

      }
      // TOF PID
    }
    //Primary particle selection,central to analysis
    maskSelectionK0Short = maskTopological | maskTrackProperties | maskK0ShortSpecific | (uint64_t(1) << selPhysPrimK0Short);
    maskSelectionLambda = maskTopological | maskTrackProperties| maskLambdaSpecific | (uint64_t(1) << selLambdaArmenteros) | (uint64_t(1) << selPhysPrimLambda) ;
    maskSelectionAntiLambda = maskTopological | maskTrackProperties | maskAntiLambdaSpecific | (uint64_t(1) << selAntiLambdaArmenteros) | (uint64_t(1) << selPhysPrimAntiLambda);

  
  }

  template <typename TV0, typename TCollision>

  uint64_t computeReconstructionBitmap(TV0 v0, TCollision collision, float rapidityLambda, float rapidityK0Short,float /*pT*/){
    uint64_t bitMap = 0;
    // Base topological variables
    if (v0.v0radius() > v0setting_radius)
      bitset(bitMap, selRadius);
    if (v0.v0radius() > v0radiusMax)
      bitset(bitMap, selRadiusMax);
    if (TMath::Abs(v0.dcapostopv()) >= dcapostopv)
      bitset(bitMap, selDCAPosToPV);
    if (TMath::Abs(v0.dcanegtopv()) >= dcanegtopv)
      bitset(bitMap, selDCANegToPV);
    if (v0.v0cosPA() > v0setting_cospa)
      bitset(bitMap, selCosPA);
    if (v0.dcaV0daughters() < v0setting_dcav0dau)
      bitset(bitMap, selDCAV0Dau);

      //rapidity
    if (TMath::Abs(rapidityLambda) < rapidityCut)
      bitset(bitMap, selLambdaRapidity);
    if (TMath::Abs(rapidityK0Short) < rapidityCut)
      bitset(bitMap, selK0ShortRapidity);



    const auto& posDaughterTrack = v0.template posTrack_as<DaughterTracks>();
    const auto& negDaughterTrack = v0.template negTrack_as<DaughterTracks>();

    //const auto& posDaughterTrack = v0.template posTrackExtra_as<DaughterTracks>();
    //const auto& negDaughterTrack = v0.template negTrackExtra_as<DaughterTracks>();



    // TPC quality flags
    if ((posDaughterTrack.tpcNClsCrossedRows()) > minTPCrows)
      bitset(bitMap, selPosGoodTPCTrack);
    if (negDaughterTrack.tpcNClsCrossedRows() > minTPCrows)
      bitset(bitMap, selNegGoodTPCTrack);


    //TPC PID
    if (fabs(posDaughterTrack.tpcNSigmaPi()) < TpcPidNsigmaCut)
      bitset(bitMap, selTPCPIDPositivePion);
    if (fabs(posDaughterTrack.tpcNSigmaPr()) < TpcPidNsigmaCut)
      bitset(bitMap, selTPCPIDPositiveProton);
    if (fabs(negDaughterTrack.tpcNSigmaPi()) < TpcPidNsigmaCut)
      bitset(bitMap, selTPCPIDNegativePion);
    if (fabs(negDaughterTrack.tpcNSigmaPr()) < TpcPidNsigmaCut)
      bitset(bitMap, selTPCPIDNegativeProton);

    //ITS quality flags
    if (posDaughterTrack.itsNCls() >= minITSclusters)
      bitset(bitMap, selPosGoodITSTrack);
    if (negDaughterTrack.itsNCls() >= minITSclusters)
      bitset(bitMap, selNegGoodITSTrack);

    // proper lifetime
    if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0 < lifetimecut->get("lifetimecutLambda"))
      bitset(bitMap, selLambdaCTau);
    if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short < lifetimecut->get("lifetimecutK0S"))
      bitset(bitMap, selK0ShortCTau);

    // armenteros
    
    if (v0.qtarm() * armPodCut < TMath::Abs(v0.alpha()) || armPodCut < 1e-4)
      bitset(bitMap, selLambdaArmenteros);
    if (v0.qtarm() * armPodCut < TMath::Abs(v0.alpha()) || armPodCut < 1e-4)
      bitset(bitMap, selAntiLambdaArmenteros);
    if (v0.qtarm() * armPodCut > TMath::Abs(v0.alpha()) || armPodCut < 1e-4)
      bitset(bitMap, selK0ShortArmenteros);
      
    return bitMap;
  }

  bool verifyMask(uint64_t bitmap, uint64_t mask)
  {
    return (bitmap & mask) == mask;
  }

  template <typename TV0>
  void analyseCandidate(TV0 v0, float pt, uint64_t selMap){
    const auto& posDaughterTrack = v0.template posTrack_as<DaughterTracks>();
    const auto& negDaughterTrack = v0.template negTrack_as<DaughterTracks>();

    //const auto& posDaughterTrack = v0.template posTrackExtra_as<DaughterTracks>();
    //const auto& negDaughterTrack = v0.template negTrackExtra_as<DaughterTracks>();

    //bool Lambdaselection = verifyMask(selMap, maskSelectionLambda);

    if(verifyMask(selMap, maskSelectionLambda) && analyseLambda){
      histos.fill(HIST("hMassLambda"), v0.mLambda());
      //histos.fill(HIST("hKs0"), v0.mK0Short());
      histos.fill(HIST("h3dMassLambda"),pt, v0.mLambda());

      histos.fill(HIST("hPosDCAToPV"), v0.dcapostopv());
      histos.fill(HIST("hNegDCAToPV"), v0.dcanegtopv());
      histos.fill(HIST("hDCADaughters"), v0.dcaV0daughters());
      histos.fill(HIST("hPointingAngle"), v0.v0cosPA());
      histos.fill(HIST("hV0Radius"), v0.v0radius());
      histos.fill(HIST("h2dPositiveITSvsTPCpts"), posDaughterTrack.tpcNClsCrossedRows(), posDaughterTrack.itsNCls());
      histos.fill(HIST("h2dNegativeITSvsTPCpts"), negDaughterTrack.tpcNClsCrossedRows(), negDaughterTrack.itsNCls());
      
      histos.fill(HIST("hPostpcCrossedRows"), posDaughterTrack.tpcNClsCrossedRows());
      histos.fill(HIST("hNegtpcCrossedRows"), negDaughterTrack.tpcNClsCrossedRows());

      histos.fill(HIST("hpositiveeta"), v0.positiveeta());
      histos.fill(HIST("hnegativeeta"), v0.negativeeta());
      histos.fill(HIST("hPostpcCrossedRowsOverFindableCls"), posDaughterTrack.tpcCrossedRowsOverFindableCls());
      histos.fill(HIST("hNegtpcCrossedRowsOverFindableCls"), negDaughterTrack.tpcCrossedRowsOverFindableCls());
      histos.fill(HIST("h2dArmenterosAll"), v0.alpha(), v0.qtarm());
      histos.fill(HIST("h2dArmenterosLambda"), v0.alpha(), v0.qtarm());
    }
    if(verifyMask(selMap,maskSelectionAntiLambda) && analyseAntiLambda){
      histos.fill(HIST("hMassAntiLambda"), v0.mAntiLambda());
      histos.fill(HIST("h2dArmenterosAntiLambda"), v0.alpha(), v0.qtarm());
    }

    if(verifyMask(selMap,maskSelectionK0Short) && analyseK0Short){
      histos.fill(HIST("hKs0"), v0.mK0Short());
      histos.fill(HIST("h2dArmenterosKshort"), v0.alpha(), v0.qtarm());
    }
    

  }



  // Defining filters for events (event selection)
  // Processed events will be already fulfilling the event selection requirements
  //Filter eventFilter = (o2::aod::evsel::sel8 == true);
  //Filter posZFilter = (nabs(o2::aod::collision::posZ) < cutzvertex);

  // Filters on V0s
  // Cannot filter on dynamic columns, so we cut on DCA to PV and DCA between daughters only
  /*
  Filter preFilterV0 = (nabs(aod::v0data::dcapostopv) > v0setting_dcapostopv &&
                        nabs(aod::v0data::dcanegtopv) > v0setting_dcanegtopv &&
                        aod::v0data::dcaV0daughters < v0setting_dcav0dau);
 */
  // Defining the type of the daughter tracks
  using DaughterTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCPi,o2::aod::pidTPCPr,o2::aod::StoredTracksExtra_001>;

  using Collisionstable = soa::Join<aod::StraCollisions, aod::StraCents, aod::StraRawCents, aod::StraEvSels>::iterator;
  //soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator
  //soa::Filtered<soa::Join<aod::Collisions
  //aod::V0Datas const& V0s,
  //DaughterTracks
  //dauTracks
  void process( soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
               aod::V0Datas const& V0s,
               DaughterTracks const&)
  {
    rLambda.fill(HIST("hEventSelection"), 0. /* all collisions */);
    // Fill the event counter
    ////////////////////////////////Event selection///////////////////////////
    rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());
    
    if (!collision.sel8()) {
        return;
      }
    rLambda.fill(HIST("hEventSelection"), 1 /* sel8 collisions */);
    if (std::abs(collision.posZ()) >= 10.f) {
      return;
    }
    histos.fill(HIST("hVertexZRec"), collision.posZ());
    rLambda.fill(HIST("hEventSelection"), 2 /* vertex-Z selected */);
    if (!collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return;
    }
    rLambda.fill(HIST("hEventSelection"), 3 /* Not at ITS ROF border */);

    if (!collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return;
    }
    
    
    rLambda.fill(HIST("hEventSelection"), 4 );
    /*
    if (!collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return;
    }
    rLambda.fill(HIST("hEventSelection"), 5 );

    if (!collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return;
    }
    rLambda.fill(HIST("hEventSelection"), 6);

    if (!collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
      return;
    }
    rLambda.fill(HIST("hEventSelection"), 7 );

    if (!collision.selection_bit(o2::aod::evsel::kIsVertexTRDmatched)) {
      return;
    }
    rLambda.fill(HIST("hEventSelection"), 8 );

    if (!collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return;
    }
    rLambda.fill(HIST("hEventSelection"), 9 );
  */
 ////////////////////////////////Event selection///////////////////////////

    // perform main analysis
    for (const auto& v0 : V0s) {
      if (std::abs(v0.negativeeta()) > daughterEtaCut || std::abs(v0.positiveeta()) > daughterEtaCut){
        continue;
      }
      rLambda.fill(HIST("h2dArmenterosAll"), v0.alpha(), v0.qtarm());

      uint64_t selMap = computeReconstructionBitmap(v0, collision, v0.yLambda(), v0.yK0Short(), v0.pt());
      selMap = selMap | (uint64_t(1) << selConsiderK0Short) | (uint64_t(1) << selConsiderLambda) | (uint64_t(1) << selConsiderAntiLambda);
      selMap = selMap | (uint64_t(1) << selPhysPrimK0Short) | (uint64_t(1) << selPhysPrimLambda) | (uint64_t(1) << selPhysPrimAntiLambda);

      //selMap = selMap |   (uint64_t(1) << selPhysPrimLambda) ;

      analyseCandidate(v0, v0.pt(), selMap);


      //const auto& posDaughterTrack = v0.posTrack_as<DaughterTracks>();
      //const auto& negDaughterTrack = v0.negTrack_as<DaughterTracks>();

      rKzeroShort.fill(HIST("hMassK0Short"), v0.mK0Short());
      rLambda.fill(HIST("hLambda"), v0.mLambda());
      
      // Cut on dynamic columns
      /*
      if (v0.v0cosPA() < v0setting_cospa)
        continue;
      if (v0.v0radius() < v0setting_radius)
        continue;

      if (TMath::Abs(posDaughterTrack.tpcNSigmaPi()) > NSigmaTPCPion) {
        continue;
      }
      if (TMath::Abs(negDaughterTrack.tpcNSigmaPi()) > NSigmaTPCPion) {
        continue;
      }
      if (TMath::Abs(v0.dcapostopv()) < dcapostopv){
        continue;
      }
      if (TMath::Abs(v0.dcanegtopv()) < dcanegtopv){
        continue;
      }
      if (v0.dcaV0daughters() > v0setting_dcav0dau){
        continue;
      }
      

      rKzeroShort.fill(HIST("hMassK0ShortSelected"), v0.mK0Short());
      rKzeroShort.fill(HIST("hDCAV0Daughters"), v0.dcaV0daughters());
      rKzeroShort.fill(HIST("hV0CosPA"), v0.v0cosPA());

      rLambda.fill(HIST("hMassLambdaSelected"), v0.mLambda());
      rLambda.fill(HIST("hDCAV0Daughters"), v0.dcaV0daughters());
      rLambda.fill(HIST("hV0CosPA"), v0.v0cosPA());
      rLambda.fill(HIST("hPosDCAToPV"), v0.dcapostopv());
      rLambda.fill(HIST("hNegDCAToPV"), v0.dcanegtopv());

      rLambda.fill(HIST("hPostpcCrossedRows"), posDaughterTrack.tpcNClsCrossedRows());
      rLambda.fill(HIST("hNegtpcCrossedRows"), negDaughterTrack.tpcNClsCrossedRows());

      rLambda.fill(HIST("hPostpcCrossedRowsOverFindableCls"), posDaughterTrack.tpcCrossedRowsOverFindableCls());
      rLambda.fill(HIST("hNegtpcCrossedRowsOverFindableCls"), negDaughterTrack.tpcCrossedRowsOverFindableCls());

      rLambda.fill(HIST("hposativeeta"), v0.positiveeta());
      rLambda.fill(HIST("hnegativeeta"), v0.negativeeta());

      rLambda.fill(HIST("hNSigmaPosProtonFromLambdas"), posDaughterTrack.tpcNSigmaPr(), posDaughterTrack.tpcInnerParam());
      rLambda.fill(HIST("hNSigmaNegPionFromLambdas"), negDaughterTrack.tpcNSigmaPi(), negDaughterTrack.tpcInnerParam());
      */
      
    
    }
  }
  PROCESS_SWITCH(strangeness_tutorial , process, "process as if real data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<strangeness_tutorial>(cfgc)};
}
