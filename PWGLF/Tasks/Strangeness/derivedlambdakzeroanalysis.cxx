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
//
// V0 analysis task
// ================
//
// This code loops over a V0Cores table and produces some
// standard analysis output. It is meant to be run over
// derived data.
//
//    Comments, questions, complaints, suggestions?
//    Please write to:
//    romain.schotter@cern.ch
//    david.dobrigkeit.chinellato@cern.ch
//

#include <Math/Vector4D.h>
#include <cmath>
#include <array>
#include <cstdlib>

#include <TFile.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

using dauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;
using v0Candidates = soa::Join<aod::V0CollRefs, aod::V0Cores, aod::V0Extras, aod::V0TOFPIDs, aod::V0TOFNSigmas>;
using v0MCCandidates = soa::Join<aod::V0CollRefs, aod::V0Cores, aod::V0MCCores, aod::V0Extras, aod::V0TOFPIDs, aod::V0TOFNSigmas, aod::V0MCMothers, aod::V0MCCollRefs>;

// simple checkers, but ensure 64 bit integers
#define bitset(var, nbit) ((var) |= (static_cast<uint64_t>(1) << static_cast<uint64_t>(nbit)))
#define bitcheck(var, nbit) ((var) & (static_cast<uint64_t>(1) << static_cast<uint64_t>(nbit)))

struct derivedlambdakzeroanalysis {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // master analysis switches
  Configurable<bool> analyseK0Short{"analyseK0Short", true, "process K0Short-like candidates"};
  Configurable<bool> analyseLambda{"analyseLambda", true, "process Lambda-like candidates"};
  Configurable<bool> analyseAntiLambda{"analyseAntiLambda", true, "process AntiLambda-like candidates"};
  Configurable<bool> calculateFeeddownMatrix{"calculateFeeddownMatrix", true, "fill feeddown matrix if MC"};
  Configurable<bool> rejectITSROFBorder{"rejectITSROFBorder", true, "reject events at ITS ROF border"};
  Configurable<bool> rejectTFBorder{"rejectTFBorder", true, "reject events at TF border"};

  // Selection criteria: acceptance
  Configurable<float> rapidityCut{"rapidityCut", 0.5, "rapidity"};
  Configurable<float> daughterEtaCut{"daughterEtaCut", 0.8, "max eta for daughters"};

  // Standard 5 topological criteria
  Configurable<float> v0cospa{"v0cospa", 0.97, "min V0 CosPA"};
  Configurable<float> dcav0dau{"dcav0dau", 1.0, "max DCA V0 Daughters (cm)"};
  Configurable<float> dcanegtopv{"dcanegtopv", .05, "min DCA Neg To PV (cm)"};
  Configurable<float> dcapostopv{"dcapostopv", .05, "min DCA Pos To PV (cm)"};
  Configurable<float> v0radius{"v0radius", 1.2, "minimum V0 radius (cm)"};
  Configurable<float> v0radiusMax{"v0radiusMax", 1E5, "maximum V0 radius (cm)"};

  // Additional selection on the AP plot (exclusive for K0Short)
  // original equation: lArmPt*5>TMath::Abs(lArmAlpha)
  Configurable<float> armPodCut{"armPodCut", 5.0f, "pT * (cut) > |alpha|, AP cut. Negative: no cut"};

  // Track quality
  Configurable<int> minTPCrows{"minTPCrows", 70, "minimum TPC crossed rows"};
  Configurable<int> minITSclusters{"minITSclusters", -1, "minimum ITS clusters"};
  Configurable<bool> skipTPConly{"skipTPConly", false, "skip V0s comprised of at least one TPC only prong"};
  Configurable<bool> requirePosITSonly{"requirePosITSonly", false, "require that positive track is ITSonly (overrides TPC quality)"};
  Configurable<bool> requireNegITSonly{"requireNegITSonly", false, "require that negative track is ITSonly (overrides TPC quality)"};

  // PID (TPC/TOF)
  Configurable<float> TpcPidNsigmaCut{"TpcPidNsigmaCut", 5, "TpcPidNsigmaCut"};
  Configurable<float> TofPidNsigmaCutLaPr{"TofPidNsigmaCutLaPr", 1e+6, "TofPidNsigmaCutLaPr"};
  Configurable<float> TofPidNsigmaCutLaPi{"TofPidNsigmaCutLaPi", 1e+6, "TofPidNsigmaCutLaPi"};
  Configurable<float> TofPidNsigmaCutK0Pi{"TofPidNsigmaCutK0Pi", 1e+6, "TofPidNsigmaCutK0Pi"};

  Configurable<bool> doCompleteTopoQA{"doCompleteTopoQA", false, "do topological variable QA histograms"};
  Configurable<bool> doTPCQA{"doTPCQA", false, "do TPC QA histograms"};
  Configurable<bool> doTOFQA{"doTOFQA", false, "do TOF QA histograms"};
  Configurable<int> doDetectPropQA{"doDetectPropQA", 0, "do Detector/ITS map QA: 0: no, 1: 4D, 2: 5D with mass"};

  Configurable<bool> doPlainTopoQA{"doPlainTopoQA", true, "do simple 1D QA of candidates"};
  Configurable<float> qaMinPt{"qaMinPt", 0.0f, "minimum pT for QA plots"};
  Configurable<float> qaMaxPt{"qaMaxPt", 1000.0f, "maximum pT for QA plots"};
  Configurable<bool> qaCentrality{"qaCentrality", false, "qa centrality flag: check base raw values"};

  // PID (TOF)
  Configurable<float> maxDeltaTimeProton{"maxDeltaTimeProton", 1e+9, "check maximum allowed time"};
  Configurable<float> maxDeltaTimePion{"maxDeltaTimePion", 1e+9, "check maximum allowed time"};

  // for MC
  Configurable<bool> doMCAssociation{"doMCAssociation", true, "if MC, do MC association"};
  Configurable<bool> doCollisionAssociationQA{"doCollisionAssociationQA", true, "check collision association"};

  static constexpr float defaultLifetimeCuts[1][2] = {{30., 20.}};
  Configurable<LabeledArray<float>> lifetimecut{"lifetimecut", {defaultLifetimeCuts[0], 2, {"lifetimecutLambda", "lifetimecutK0S"}}, "lifetimecut"};

  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for analysis"};
  ConfigurableAxis axisPtXi{"axisPtXi", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for feeddown from Xi"};
  ConfigurableAxis axisPtCoarse{"axisPtCoarse", {VARIABLE_WIDTH, 0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 7.0f, 10.0f, 15.0f}, "pt axis for QA"};
  ConfigurableAxis axisK0Mass{"axisK0Mass", {200, 0.4f, 0.6f}, ""};
  ConfigurableAxis axisLambdaMass{"axisLambdaMass", {200, 1.101f, 1.131f}, ""};
  ConfigurableAxis axisCentrality{"axisCentrality", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f}, "Centrality"};

  ConfigurableAxis axisRawCentrality{"axisRawCentrality", {VARIABLE_WIDTH, 0.000f, 52.320f, 75.400f, 95.719f, 115.364f, 135.211f, 155.791f, 177.504f, 200.686f, 225.641f, 252.645f, 281.906f, 313.850f, 348.302f, 385.732f, 426.307f, 470.146f, 517.555f, 568.899f, 624.177f, 684.021f, 748.734f, 818.078f, 892.577f, 973.087f, 1058.789f, 1150.915f, 1249.319f, 1354.279f, 1465.979f, 1584.790f, 1710.778f, 1844.863f, 1985.746f, 2134.643f, 2291.610f, 2456.943f, 2630.653f, 2813.959f, 3006.631f, 3207.229f, 3417.641f, 3637.318f, 3865.785f, 4104.997f, 4354.938f, 4615.786f, 4885.335f, 5166.555f, 5458.021f, 5762.584f, 6077.881f, 6406.834f, 6746.435f, 7097.958f, 7462.579f, 7839.165f, 8231.629f, 8635.640f, 9052.000f, 9484.268f, 9929.111f, 10389.350f, 10862.059f, 11352.185f, 11856.823f, 12380.371f, 12920.401f, 13476.971f, 14053.087f, 14646.190f, 15258.426f, 15890.617f, 16544.433f, 17218.024f, 17913.465f, 18631.374f, 19374.983f, 20136.700f, 20927.783f, 21746.796f, 22590.880f, 23465.734f, 24372.274f, 25314.351f, 26290.488f, 27300.899f, 28347.512f, 29436.133f, 30567.840f, 31746.818f, 32982.664f, 34276.329f, 35624.859f, 37042.588f, 38546.609f, 40139.742f, 41837.980f, 43679.429f, 45892.130f, 400000.000f}, "raw centrality signal"}; // for QA

  // topological variable QA axes
  ConfigurableAxis axisDCAtoPV{"axisDCAtoPV", {20, 0.0f, 1.0f}, "DCA (cm)"};
  ConfigurableAxis axisDCAdau{"axisDCAdau", {20, 0.0f, 2.0f}, "DCA (cm)"};
  ConfigurableAxis axisPointingAngle{"axisPointingAngle", {20, 0.0f, 2.0f}, "pointing angle (rad)"};
  ConfigurableAxis axisV0Radius{"axisV0Radius", {20, 0.0f, 60.0f}, "V0 2D radius (cm)"};
  ConfigurableAxis axisNsigmaTPC{"axisNsigmaTPC", {200, -10.0f, 10.0f}, "N sigma TPC"};
  ConfigurableAxis axisTPCsignal{"axisTPCsignal", {200, 0.0f, 200.0f}, "TPC signal"};
  ConfigurableAxis axisTOFdeltaT{"axisTOFdeltaT", {200, -5000.0f, 5000.0f}, "TOF Delta T (ps)"};

  // AP plot axes
  ConfigurableAxis axisAPAlpha{"axisAPAlpha", {220, -1.1f, 1.1f}, "V0 AP alpha"};
  ConfigurableAxis axisAPQt{"axisAPQt", {220, 0.0f, 0.5f}, "V0 AP alpha"};

  // Track quality axes
  ConfigurableAxis axisTPCrows{"axisTPCrows", {160, 0.0f, 160.0f}, "N TPC rows"};
  ConfigurableAxis axisITSclus{"axisITSclus", {7, 0.0f, 7.0f}, "N ITS Clusters"};
  ConfigurableAxis axisITScluMap{"axisITSMap", {128, -0.5f, 127.5f}, "ITS Cluster map"};
  ConfigurableAxis axisDetMap{"axisDetMap", {16, -0.5f, 15.5f}, "Detector use map"};
  ConfigurableAxis axisITScluMapCoarse{"axisITScluMapCoarse", {13, -0.5f, 12.5f}, "ITS Coarse cluster map"};
  ConfigurableAxis axisDetMapCoarse{"axisDetMapCoarse", {4, -0.5f, 3.5f}, "Detector Coarse user map"};

  // MC coll assoc QA axis
  ConfigurableAxis axisMonteCarloNch{"axisMonteCarloNch", {300, 0.0f, 3000.0f}, "N_{ch} MC"};

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
                              selTOFDeltaTPositiveProtonLambda,
                              selTOFDeltaTPositivePionLambda,
                              selTOFDeltaTPositivePionK0Short,
                              selTOFDeltaTNegativeProtonLambda,
                              selTOFDeltaTNegativePionLambda,
                              selTOFDeltaTNegativePionK0Short,
                              selTOFNSigmaPositiveProtonLambda, // Nsigma
                              selTOFNSigmaPositivePionLambda,   // Nsigma
                              selTOFNSigmaPositivePionK0Short,  // Nsigma
                              selTOFNSigmaNegativeProtonLambda, // Nsigma
                              selTOFNSigmaNegativePionLambda,   // Nsigma
                              selTOFNSigmaNegativePionK0Short,  // Nsigma
                              selK0ShortCTau,
                              selLambdaCTau,
                              selK0ShortArmenteros,
                              selPosGoodTPCTrack, // at least min # TPC rows
                              selNegGoodTPCTrack, // at least min # TPC rows
                              selPosGoodITSTrack, // at least min # ITS clusters
                              selNegGoodITSTrack, // at least min # ITS clusters
                              selPosItsOnly,
                              selNegItsOnly,
                              selPosNotTPCOnly,
                              selNegNotTPCOnly,
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
  uint64_t secondaryMaskSelectionAntiLambda;

  void init(InitContext const&)
  {
    // initialise bit masks
    maskTopological = (uint64_t(1) << selCosPA) | (uint64_t(1) << selRadius) | (uint64_t(1) << selDCANegToPV) | (uint64_t(1) << selDCAPosToPV) | (uint64_t(1) << selDCAV0Dau) | (uint64_t(1) << selRadiusMax);
    maskTopoNoV0Radius = (uint64_t(1) << selCosPA) | (uint64_t(1) << selDCANegToPV) | (uint64_t(1) << selDCAPosToPV) | (uint64_t(1) << selDCAV0Dau) | (uint64_t(1) << selRadiusMax);
    maskTopoNoDCANegToPV = (uint64_t(1) << selCosPA) | (uint64_t(1) << selRadius) | (uint64_t(1) << selDCAPosToPV) | (uint64_t(1) << selDCAV0Dau) | (uint64_t(1) << selRadiusMax);
    maskTopoNoDCAPosToPV = (uint64_t(1) << selCosPA) | (uint64_t(1) << selRadius) | (uint64_t(1) << selDCANegToPV) | (uint64_t(1) << selDCAV0Dau) | (uint64_t(1) << selRadiusMax);
    maskTopoNoCosPA = (uint64_t(1) << selRadius) | (uint64_t(1) << selDCANegToPV) | (uint64_t(1) << selDCAPosToPV) | (uint64_t(1) << selDCAV0Dau) | (uint64_t(1) << selRadiusMax);
    maskTopoNoDCAV0Dau = (uint64_t(1) << selCosPA) | (uint64_t(1) << selRadius) | (uint64_t(1) << selDCANegToPV) | (uint64_t(1) << selDCAPosToPV) | (uint64_t(1) << selRadiusMax);

    maskK0ShortSpecific = (uint64_t(1) << selK0ShortRapidity) | (uint64_t(1) << selK0ShortCTau) | (uint64_t(1) << selK0ShortArmenteros) | (uint64_t(1) << selConsiderK0Short);
    maskLambdaSpecific = (uint64_t(1) << selLambdaRapidity) | (uint64_t(1) << selLambdaCTau) | (uint64_t(1) << selConsiderLambda);
    maskAntiLambdaSpecific = (uint64_t(1) << selLambdaRapidity) | (uint64_t(1) << selLambdaCTau) | (uint64_t(1) << selConsiderAntiLambda);

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
      if (TofPidNsigmaCutK0Pi < 1e+5) // safeguard for no cut
        maskK0ShortSpecific = maskK0ShortSpecific | (uint64_t(1) << selTOFNSigmaPositivePionK0Short) | (uint64_t(1) << selTOFDeltaTPositivePionK0Short);
      if (TofPidNsigmaCutLaPr < 1e+5) // safeguard for no cut
        maskLambdaSpecific = maskLambdaSpecific | (uint64_t(1) << selTOFNSigmaPositiveProtonLambda) | (uint64_t(1) << selTOFDeltaTPositiveProtonLambda);
      if (TofPidNsigmaCutLaPi < 1e+5) // safeguard for no cut
        maskAntiLambdaSpecific = maskAntiLambdaSpecific | (uint64_t(1) << selTOFNSigmaPositivePionLambda) | (uint64_t(1) << selTOFDeltaTPositivePionLambda);
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
      if (TofPidNsigmaCutK0Pi < 1e+5) // safeguard for no cut
        maskK0ShortSpecific = maskK0ShortSpecific | (uint64_t(1) << selTOFNSigmaNegativePionK0Short) | (uint64_t(1) << selTOFDeltaTNegativePionK0Short);
      if (TofPidNsigmaCutLaPr < 1e+5) // safeguard for no cut
        maskLambdaSpecific = maskLambdaSpecific | (uint64_t(1) << selTOFNSigmaNegativePionLambda) | (uint64_t(1) << selTOFDeltaTNegativePionLambda);
      if (TofPidNsigmaCutLaPi < 1e+5) // safeguard for no cut
        maskAntiLambdaSpecific = maskAntiLambdaSpecific | (uint64_t(1) << selTOFNSigmaNegativeProtonLambda) | (uint64_t(1) << selTOFDeltaTNegativeProtonLambda);
    }

    if (skipTPConly) {
      maskK0ShortSpecific = maskK0ShortSpecific | (uint64_t(1) << selPosNotTPCOnly) | (uint64_t(1) << selNegNotTPCOnly);
      maskLambdaSpecific = maskLambdaSpecific | (uint64_t(1) << selPosNotTPCOnly) | (uint64_t(1) << selNegNotTPCOnly);
      maskAntiLambdaSpecific = maskAntiLambdaSpecific | (uint64_t(1) << selPosNotTPCOnly) | (uint64_t(1) << selNegNotTPCOnly);
    }

    // Primary particle selection, central to analysis
    maskSelectionK0Short = maskTopological | maskTrackProperties | maskK0ShortSpecific | (uint64_t(1) << selPhysPrimK0Short);
    maskSelectionLambda = maskTopological | maskTrackProperties | maskLambdaSpecific | (uint64_t(1) << selPhysPrimLambda);
    maskSelectionAntiLambda = maskTopological | maskTrackProperties | maskAntiLambdaSpecific | (uint64_t(1) << selPhysPrimAntiLambda);

    // No primary requirement for feeddown matrix
    secondaryMaskSelectionLambda = maskTopological | maskTrackProperties | maskLambdaSpecific;
    secondaryMaskSelectionAntiLambda = maskTopological | maskTrackProperties | maskAntiLambdaSpecific;

    // Event Counters
    histos.add("hEventSelection", "hEventSelection", kTH1F, {{10, -0.5f, +9.5f}});
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(1, "All collisions");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(2, "sel8 cut");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(3, "posZ cut");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(4, "kNoITSROFrameBorder");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(5, "kNoTimeFrameBorder");

    histos.add("hEventCentrality", "hEventCentrality", kTH1F, {{100, 0.0f, +100.0f}});
    histos.add("hCentralityVsNch", "hCentralityVsNch", kTH2F, {axisCentrality, {500, 0.0f, +2000.0f}});

    // for QA and test purposes
    auto hRawCentrality = histos.add<TH1>("hRawCentrality", "hRawCentrality", kTH1F, {axisRawCentrality});

    for (int ii = 1; ii < 101; ii++) {
      float value = 100.5f - static_cast<float>(ii);
      hRawCentrality->SetBinContent(ii, value);
    }

    // histograms versus mass
    if (analyseK0Short) {
      histos.add("h3dMassK0Short", "h3dMassK0Short", kTH3F, {axisCentrality, axisPt, axisK0Mass});
      if (doTPCQA) {
        histos.add("K0Short/h3dPosNsigmaTPC", "h3dPosNsigmaTPC", kTH3F, {axisCentrality, axisPtCoarse, axisNsigmaTPC});
        histos.add("K0Short/h3dNegNsigmaTPC", "h3dNegNsigmaTPC", kTH3F, {axisCentrality, axisPtCoarse, axisNsigmaTPC});
        histos.add("K0Short/h3dPosTPCsignal", "h3dPosTPCsignal", kTH3F, {axisCentrality, axisPtCoarse, axisTPCsignal});
        histos.add("K0Short/h3dNegTPCsignal", "h3dNegTPCsignal", kTH3F, {axisCentrality, axisPtCoarse, axisTPCsignal});
        histos.add("K0Short/h3dPosNsigmaTPCvsTrackPtot", "h3dPosNsigmaTPCvsTrackPtot", kTH3F, {axisCentrality, axisPtCoarse, axisNsigmaTPC});
        histos.add("K0Short/h3dNegNsigmaTPCvsTrackPtot", "h3dNegNsigmaTPCvsTrackPtot", kTH3F, {axisCentrality, axisPtCoarse, axisNsigmaTPC});
        histos.add("K0Short/h3dPosTPCsignalVsTrackPtot", "h3dPosTPCsignalVsTrackPtot", kTH3F, {axisCentrality, axisPtCoarse, axisTPCsignal});
        histos.add("K0Short/h3dNegTPCsignalVsTrackPtot", "h3dNegTPCsignalVsTrackPtot", kTH3F, {axisCentrality, axisPtCoarse, axisTPCsignal});
        histos.add("K0Short/h3dPosNsigmaTPCvsTrackPt", "h3dPosNsigmaTPCvsTrackPt", kTH3F, {axisCentrality, axisPtCoarse, axisNsigmaTPC});
        histos.add("K0Short/h3dNegNsigmaTPCvsTrackPt", "h3dNegNsigmaTPCvsTrackPt", kTH3F, {axisCentrality, axisPtCoarse, axisNsigmaTPC});
        histos.add("K0Short/h3dPosTPCsignalVsTrackPt", "h3dPosTPCsignalVsTrackPt", kTH3F, {axisCentrality, axisPtCoarse, axisTPCsignal});
        histos.add("K0Short/h3dNegTPCsignalVsTrackPt", "h3dNegTPCsignalVsTrackPt", kTH3F, {axisCentrality, axisPtCoarse, axisTPCsignal});
      }
      if (doTOFQA) {
        histos.add("K0Short/h3dPosTOFdeltaT", "h3dPosTOFdeltaT", kTH3F, {axisCentrality, axisPtCoarse, axisTOFdeltaT});
        histos.add("K0Short/h3dNegTOFdeltaT", "h3dNegTOFdeltaT", kTH3F, {axisCentrality, axisPtCoarse, axisTOFdeltaT});
        histos.add("K0Short/h3dPosTOFdeltaTvsTrackPtot", "h3dPosTOFdeltaTvsTrackPtot", kTH3F, {axisCentrality, axisPtCoarse, axisTOFdeltaT});
        histos.add("K0Short/h3dNegTOFdeltaTvsTrackPtot", "h3dNegTOFdeltaTvsTrackPtot", kTH3F, {axisCentrality, axisPtCoarse, axisTOFdeltaT});
        histos.add("K0Short/h3dPosTOFdeltaTvsTrackPt", "h3dPosTOFdeltaTvsTrackPt", kTH3F, {axisCentrality, axisPtCoarse, axisTOFdeltaT});
        histos.add("K0Short/h3dNegTOFdeltaTvsTrackPt", "h3dNegTOFdeltaTvsTrackPt", kTH3F, {axisCentrality, axisPtCoarse, axisTOFdeltaT});
      }
      if (doCollisionAssociationQA) {
        histos.add("K0Short/h2dPtVsNch", "h2dPtVsNch", kTH2F, {axisMonteCarloNch, axisPt});
        histos.add("K0Short/h2dPtVsNch_BadCollAssig", "h2dPtVsNch_BadCollAssig", kTH2F, {axisMonteCarloNch, axisPt});
      }
      if (doDetectPropQA == 1) {
        histos.add("K0Short/h6dDetectPropVsCentrality", "h6dDetectPropVsCentrality", kTHnF, {axisCentrality, axisDetMapCoarse, axisITScluMapCoarse, axisDetMapCoarse, axisITScluMapCoarse, axisPtCoarse});
        histos.add("K0Short/h4dPosDetectPropVsCentrality", "h4dPosDetectPropVsCentrality", kTHnF, {axisCentrality, axisDetMap, axisITScluMap, axisPtCoarse});
        histos.add("K0Short/h4dNegDetectPropVsCentrality", "h4dNegDetectPropVsCentrality", kTHnF, {axisCentrality, axisDetMap, axisITScluMap, axisPtCoarse});
      }
      if (doDetectPropQA == 2) {
        histos.add("K0Short/h7dDetectPropVsCentrality", "h7dDetectPropVsCentrality", kTHnF, {axisCentrality, axisDetMapCoarse, axisITScluMapCoarse, axisDetMapCoarse, axisITScluMapCoarse, axisPtCoarse, axisK0Mass});
        histos.add("K0Short/h5dPosDetectPropVsCentrality", "h5dPosDetectPropVsCentrality", kTHnF, {axisCentrality, axisDetMap, axisITScluMap, axisPtCoarse, axisK0Mass});
        histos.add("K0Short/h5dNegDetectPropVsCentrality", "h5dNegDetectPropVsCentrality", kTHnF, {axisCentrality, axisDetMap, axisITScluMap, axisPtCoarse, axisK0Mass});
      }
    }
    if (analyseLambda) {
      histos.add("h3dMassLambda", "h3dMassLambda", kTH3F, {axisCentrality, axisPt, axisLambdaMass});
      if (doTPCQA) {
        histos.add("Lambda/h3dPosNsigmaTPC", "h3dPosNsigmaTPC", kTH3F, {axisCentrality, axisPtCoarse, axisNsigmaTPC});
        histos.add("Lambda/h3dNegNsigmaTPC", "h3dNegNsigmaTPC", kTH3F, {axisCentrality, axisPtCoarse, axisNsigmaTPC});
        histos.add("Lambda/h3dPosTPCsignal", "h3dPosTPCsignal", kTH3F, {axisCentrality, axisPtCoarse, axisTPCsignal});
        histos.add("Lambda/h3dNegTPCsignal", "h3dNegTPCsignal", kTH3F, {axisCentrality, axisPtCoarse, axisTPCsignal});
        histos.add("Lambda/h3dPosNsigmaTPCvsTrackPtot", "h3dPosNsigmaTPCvsTrackPtot", kTH3F, {axisCentrality, axisPtCoarse, axisNsigmaTPC});
        histos.add("Lambda/h3dNegNsigmaTPCvsTrackPtot", "h3dNegNsigmaTPCvsTrackPtot", kTH3F, {axisCentrality, axisPtCoarse, axisNsigmaTPC});
        histos.add("Lambda/h3dPosTPCsignalVsTrackPtot", "h3dPosTPCsignalVsTrackPtot", kTH3F, {axisCentrality, axisPtCoarse, axisTPCsignal});
        histos.add("Lambda/h3dNegTPCsignalVsTrackPtot", "h3dNegTPCsignalVsTrackPtot", kTH3F, {axisCentrality, axisPtCoarse, axisTPCsignal});
        histos.add("Lambda/h3dPosNsigmaTPCvsTrackPt", "h3dPosNsigmaTPCvsTrackPt", kTH3F, {axisCentrality, axisPtCoarse, axisNsigmaTPC});
        histos.add("Lambda/h3dNegNsigmaTPCvsTrackPt", "h3dNegNsigmaTPCvsTrackPt", kTH3F, {axisCentrality, axisPtCoarse, axisNsigmaTPC});
        histos.add("Lambda/h3dPosTPCsignalVsTrackPt", "h3dPosTPCsignalVsTrackPt", kTH3F, {axisCentrality, axisPtCoarse, axisTPCsignal});
        histos.add("Lambda/h3dNegTPCsignalVsTrackPt", "h3dNegTPCsignalVsTrackPt", kTH3F, {axisCentrality, axisPtCoarse, axisTPCsignal});
      }
      if (doTOFQA) {
        histos.add("Lambda/h3dPosTOFdeltaT", "h3dPosTOFdeltaT", kTH3F, {axisCentrality, axisPtCoarse, axisTOFdeltaT});
        histos.add("Lambda/h3dNegTOFdeltaT", "h3dNegTOFdeltaT", kTH3F, {axisCentrality, axisPtCoarse, axisTOFdeltaT});
        histos.add("Lambda/h3dPosTOFdeltaTvsTrackPtot", "h3dPosTOFdeltaTvsTrackPtot", kTH3F, {axisCentrality, axisPtCoarse, axisTOFdeltaT});
        histos.add("Lambda/h3dNegTOFdeltaTvsTrackPtot", "h3dNegTOFdeltaTvsTrackPtot", kTH3F, {axisCentrality, axisPtCoarse, axisTOFdeltaT});
        histos.add("Lambda/h3dPosTOFdeltaTvsTrackPt", "h3dPosTOFdeltaTvsTrackPt", kTH3F, {axisCentrality, axisPtCoarse, axisTOFdeltaT});
        histos.add("Lambda/h3dNegTOFdeltaTvsTrackPt", "h3dNegTOFdeltaTvsTrackPt", kTH3F, {axisCentrality, axisPtCoarse, axisTOFdeltaT});
      }
      if (doCollisionAssociationQA) {
        histos.add("Lambda/h2dPtVsNch", "h2dPtVsNch", kTH2F, {axisMonteCarloNch, axisPt});
        histos.add("Lambda/h2dPtVsNch_BadCollAssig", "h2dPtVsNch_BadCollAssig", kTH2F, {axisMonteCarloNch, axisPt});
      }
      if (doDetectPropQA == 1) {
        histos.add("Lambda/h6dDetectPropVsCentrality", "h6dDetectPropVsCentrality", kTHnF, {axisCentrality, axisDetMapCoarse, axisITScluMapCoarse, axisDetMapCoarse, axisITScluMapCoarse, axisPtCoarse});
        histos.add("Lambda/h4dPosDetectPropVsCentrality", "h4dPosDetectPropVsCentrality", kTHnF, {axisCentrality, axisDetMap, axisITScluMap, axisPtCoarse});
        histos.add("Lambda/h4dNegDetectPropVsCentrality", "h4dNegDetectPropVsCentrality", kTHnF, {axisCentrality, axisDetMap, axisITScluMap, axisPtCoarse});
      }
      if (doDetectPropQA == 2) {
        histos.add("Lambda/h7dDetectPropVsCentrality", "h7dDetectPropVsCentrality", kTHnF, {axisCentrality, axisDetMapCoarse, axisITScluMapCoarse, axisDetMapCoarse, axisITScluMapCoarse, axisPtCoarse, axisLambdaMass});
        histos.add("Lambda/h5dPosDetectPropVsCentrality", "h5dPosDetectPropVsCentrality", kTHnF, {axisCentrality, axisDetMap, axisITScluMap, axisPtCoarse, axisLambdaMass});
        histos.add("Lambda/h5dNegDetectPropVsCentrality", "h5dNegDetectPropVsCentrality", kTHnF, {axisCentrality, axisDetMap, axisITScluMap, axisPtCoarse, axisLambdaMass});
      }
    }
    if (analyseAntiLambda) {
      histos.add("h3dMassAntiLambda", "h3dMassAntiLambda", kTH3F, {axisCentrality, axisPt, axisLambdaMass});
      if (doTPCQA) {
        histos.add("AntiLambda/h3dPosNsigmaTPC", "h3dPosNsigmaTPC", kTH3F, {axisCentrality, axisPtCoarse, axisNsigmaTPC});
        histos.add("AntiLambda/h3dNegNsigmaTPC", "h3dNegNsigmaTPC", kTH3F, {axisCentrality, axisPtCoarse, axisNsigmaTPC});
        histos.add("AntiLambda/h3dPosTPCsignal", "h3dPosTPCsignal", kTH3F, {axisCentrality, axisPtCoarse, axisTPCsignal});
        histos.add("AntiLambda/h3dNegTPCsignal", "h3dNegTPCsignal", kTH3F, {axisCentrality, axisPtCoarse, axisTPCsignal});
        histos.add("AntiLambda/h3dPosNsigmaTPCvsTrackPtot", "h3dPosNsigmaTPCvsTrackPtot", kTH3F, {axisCentrality, axisPtCoarse, axisNsigmaTPC});
        histos.add("AntiLambda/h3dNegNsigmaTPCvsTrackPtot", "h3dNegNsigmaTPCvsTrackPtot", kTH3F, {axisCentrality, axisPtCoarse, axisNsigmaTPC});
        histos.add("AntiLambda/h3dPosTPCsignalVsTrackPtot", "h3dPosTPCsignalVsTrackPtot", kTH3F, {axisCentrality, axisPtCoarse, axisTPCsignal});
        histos.add("AntiLambda/h3dNegTPCsignalVsTrackPtot", "h3dNegTPCsignalVsTrackPtot", kTH3F, {axisCentrality, axisPtCoarse, axisTPCsignal});
        histos.add("AntiLambda/h3dPosNsigmaTPCvsTrackPt", "h3dPosNsigmaTPCvsTrackPt", kTH3F, {axisCentrality, axisPtCoarse, axisNsigmaTPC});
        histos.add("AntiLambda/h3dNegNsigmaTPCvsTrackPt", "h3dNegNsigmaTPCvsTrackPt", kTH3F, {axisCentrality, axisPtCoarse, axisNsigmaTPC});
        histos.add("AntiLambda/h3dPosTPCsignalVsTrackPt", "h3dPosTPCsignalVsTrackPt", kTH3F, {axisCentrality, axisPtCoarse, axisTPCsignal});
        histos.add("AntiLambda/h3dNegTPCsignalVsTrackPt", "h3dNegTPCsignalVsTrackPt", kTH3F, {axisCentrality, axisPtCoarse, axisTPCsignal});
      }
      if (doTOFQA) {
        histos.add("AntiLambda/h3dPosTOFdeltaT", "h3dPosTOFdeltaT", kTH3F, {axisCentrality, axisPtCoarse, axisTOFdeltaT});
        histos.add("AntiLambda/h3dNegTOFdeltaT", "h3dNegTOFdeltaT", kTH3F, {axisCentrality, axisPtCoarse, axisTOFdeltaT});
        histos.add("AntiLambda/h3dPosTOFdeltaTvsTrackPtot", "h3dPosTOFdeltaTvsTrackPtot", kTH3F, {axisCentrality, axisPtCoarse, axisTOFdeltaT});
        histos.add("AntiLambda/h3dNegTOFdeltaTvsTrackPtot", "h3dNegTOFdeltaTvsTrackPtot", kTH3F, {axisCentrality, axisPtCoarse, axisTOFdeltaT});
        histos.add("AntiLambda/h3dPosTOFdeltaTvsTrackPt", "h3dPosTOFdeltaTvsTrackPt", kTH3F, {axisCentrality, axisPtCoarse, axisTOFdeltaT});
        histos.add("AntiLambda/h3dNegTOFdeltaTvsTrackPt", "h3dNegTOFdeltaTvsTrackPt", kTH3F, {axisCentrality, axisPtCoarse, axisTOFdeltaT});
      }
      if (doCollisionAssociationQA) {
        histos.add("AntiLambda/h2dPtVsNch", "h2dPtVsNch", kTH2F, {axisMonteCarloNch, axisPt});
        histos.add("AntiLambda/h2dPtVsNch_BadCollAssig", "h2dPtVsNch_BadCollAssig", kTH2F, {axisMonteCarloNch, axisPt});
      }
      if (doDetectPropQA == 1) {
        histos.add("AntiLambda/h6dDetectPropVsCentrality", "h6dDetectPropVsCentrality", kTHnF, {axisCentrality, axisDetMapCoarse, axisITScluMapCoarse, axisDetMapCoarse, axisITScluMapCoarse, axisPtCoarse});
        histos.add("AntiLambda/h4dPosDetectPropVsCentrality", "h4dPosDetectPropVsCentrality", kTHnF, {axisCentrality, axisDetMap, axisITScluMap, axisPtCoarse});
        histos.add("AntiLambda/h4dNegDetectPropVsCentrality", "h4dNegDetectPropVsCentrality", kTHnF, {axisCentrality, axisDetMap, axisITScluMap, axisPtCoarse});
      }
      if (doDetectPropQA == 2) {
        histos.add("AntiLambda/h7dDetectPropVsCentrality", "h7dDetectPropVsCentrality", kTHnF, {axisCentrality, axisDetMapCoarse, axisITScluMapCoarse, axisDetMapCoarse, axisITScluMapCoarse, axisPtCoarse, axisLambdaMass});
        histos.add("AntiLambda/h5dPosDetectPropVsCentrality", "h5dPosDetectPropVsCentrality", kTHnF, {axisCentrality, axisDetMap, axisITScluMap, axisPtCoarse, axisLambdaMass});
        histos.add("AntiLambda/h5dNegDetectPropVsCentrality", "h5dNegDetectPropVsCentrality", kTHnF, {axisCentrality, axisDetMap, axisITScluMap, axisPtCoarse, axisLambdaMass});
      }
    }

    if (analyseLambda && calculateFeeddownMatrix && doprocessMonteCarlo)
      histos.add("h3dLambdaFeeddown", "h3dLambdaFeeddown", kTH3F, {axisCentrality, axisPt, axisPtXi});
    if (analyseAntiLambda && calculateFeeddownMatrix && doprocessMonteCarlo)
      histos.add("h3dAntiLambdaFeeddown", "h3dAntiLambdaFeeddown", kTH3F, {axisCentrality, axisPt, axisPtXi});

    // demo // fast
    histos.add("hMassK0Short", "hMassK0Short", kTH1F, {axisK0Mass});

    // QA histograms if requested
    if (doCompleteTopoQA) {
      // initialize for K0short...
      if (analyseK0Short) {
        histos.add("K0Short/h4dPosDCAToPV", "h4dPosDCAToPV", kTHnF, {axisCentrality, axisPtCoarse, axisK0Mass, axisDCAtoPV});
        histos.add("K0Short/h4dNegDCAToPV", "h4dNegDCAToPV", kTHnF, {axisCentrality, axisPtCoarse, axisK0Mass, axisDCAtoPV});
        histos.add("K0Short/h4dDCADaughters", "h4dDCADaughters", kTHnF, {axisCentrality, axisPtCoarse, axisK0Mass, axisDCAdau});
        histos.add("K0Short/h4dPointingAngle", "h4dPointingAngle", kTHnF, {axisCentrality, axisPtCoarse, axisK0Mass, axisPointingAngle});
        histos.add("K0Short/h4dV0Radius", "h4dV0Radius", kTHnF, {axisCentrality, axisPtCoarse, axisK0Mass, axisV0Radius});
      }
      if (analyseLambda) {
        histos.add("Lambda/h4dPosDCAToPV", "h4dPosDCAToPV", kTHnF, {axisCentrality, axisPtCoarse, axisLambdaMass, axisDCAtoPV});
        histos.add("Lambda/h4dNegDCAToPV", "h4dNegDCAToPV", kTHnF, {axisCentrality, axisPtCoarse, axisLambdaMass, axisDCAtoPV});
        histos.add("Lambda/h4dDCADaughters", "h4dDCADaughters", kTHnF, {axisCentrality, axisPtCoarse, axisLambdaMass, axisDCAdau});
        histos.add("Lambda/h4dPointingAngle", "h4dPointingAngle", kTHnF, {axisCentrality, axisPtCoarse, axisLambdaMass, axisPointingAngle});
        histos.add("Lambda/h4dV0Radius", "h4dV0Radius", kTHnF, {axisCentrality, axisPtCoarse, axisLambdaMass, axisV0Radius});
      }
      if (analyseAntiLambda) {
        histos.add("AntiLambda/h4dPosDCAToPV", "h4dPosDCAToPV", kTHnF, {axisCentrality, axisPtCoarse, axisLambdaMass, axisDCAtoPV});
        histos.add("AntiLambda/h4dNegDCAToPV", "h4dNegDCAToPV", kTHnF, {axisCentrality, axisPtCoarse, axisLambdaMass, axisDCAtoPV});
        histos.add("AntiLambda/h4dDCADaughters", "h4dDCADaughters", kTHnF, {axisCentrality, axisPtCoarse, axisLambdaMass, axisDCAdau});
        histos.add("AntiLambda/h4dPointingAngle", "h4dPointingAngle", kTHnF, {axisCentrality, axisPtCoarse, axisLambdaMass, axisPointingAngle});
        histos.add("AntiLambda/h4dV0Radius", "h4dV0Radius", kTHnF, {axisCentrality, axisPtCoarse, axisLambdaMass, axisV0Radius});
      }
    }

    if (doPlainTopoQA) {
      // All candidates received
      histos.add("hPosDCAToPV", "hPosDCAToPV", kTH1F, {axisDCAtoPV});
      histos.add("hNegDCAToPV", "hNegDCAToPV", kTH1F, {axisDCAtoPV});
      histos.add("hDCADaughters", "hDCADaughters", kTH1F, {axisDCAdau});
      histos.add("hPointingAngle", "hPointingAngle", kTH1F, {axisPointingAngle});
      histos.add("hV0Radius", "hV0Radius", kTH1F, {axisV0Radius});
      histos.add("h2dPositiveITSvsTPCpts", "h2dPositiveITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
      histos.add("h2dNegativeITSvsTPCpts", "h2dNegativeITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
      if (analyseK0Short) {
        histos.add("K0Short/hPosDCAToPV", "hPosDCAToPV", kTH1F, {axisDCAtoPV});
        histos.add("K0Short/hNegDCAToPV", "hNegDCAToPV", kTH1F, {axisDCAtoPV});
        histos.add("K0Short/hDCADaughters", "hDCADaughters", kTH1F, {axisDCAdau});
        histos.add("K0Short/hPointingAngle", "hPointingAngle", kTH1F, {axisPointingAngle});
        histos.add("K0Short/hV0Radius", "hV0Radius", kTH1F, {axisV0Radius});
        histos.add("K0Short/h2dPositiveITSvsTPCpts", "h2dPositiveITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
        histos.add("K0Short/h2dNegativeITSvsTPCpts", "h2dNegativeITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
      }
      if (analyseLambda) {
        histos.add("Lambda/hPosDCAToPV", "hPosDCAToPV", kTH1F, {axisDCAtoPV});
        histos.add("Lambda/hNegDCAToPV", "hNegDCAToPV", kTH1F, {axisDCAtoPV});
        histos.add("Lambda/hDCADaughters", "hDCADaughters", kTH1F, {axisDCAdau});
        histos.add("Lambda/hPointingAngle", "hPointingAngle", kTH1F, {axisPointingAngle});
        histos.add("Lambda/hV0Radius", "hV0Radius", kTH1F, {axisV0Radius});
        histos.add("Lambda/h2dPositiveITSvsTPCpts", "h2dPositiveITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
        histos.add("Lambda/h2dNegativeITSvsTPCpts", "h2dNegativeITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
      }
      if (analyseAntiLambda) {
        histos.add("AntiLambda/hPosDCAToPV", "hPosDCAToPV", kTH1F, {axisDCAtoPV});
        histos.add("AntiLambda/hNegDCAToPV", "hNegDCAToPV", kTH1F, {axisDCAtoPV});
        histos.add("AntiLambda/hDCADaughters", "hDCADaughters", kTH1F, {axisDCAdau});
        histos.add("AntiLambda/hPointingAngle", "hPointingAngle", kTH1F, {axisPointingAngle});
        histos.add("AntiLambda/hV0Radius", "hV0Radius", kTH1F, {axisV0Radius});
        histos.add("AntiLambda/h2dPositiveITSvsTPCpts", "h2dPositiveITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
        histos.add("AntiLambda/h2dNegativeITSvsTPCpts", "h2dNegativeITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
      }
    }

    // Check if doing the right thing in AP space please
    histos.add("GeneralQA/h2dArmenterosAll", "h2dArmenterosAll", kTH2F, {axisAPAlpha, axisAPQt});
    histos.add("GeneralQA/h2dArmenterosSelected", "h2dArmenterosSelected", kTH2F, {axisAPAlpha, axisAPQt});

    // Creation of histograms: MC generated
    if (doprocessBinnedGenerated) {
      histos.add("h2dGenK0Short", "h2dGenK0Short", kTH2D, {axisCentrality, axisPt});
      histos.add("h2dGenLambda", "h2dGenLambda", kTH2D, {axisCentrality, axisPt});
      histos.add("h2dGenAntiLambda", "h2dGenAntiLambda", kTH2D, {axisCentrality, axisPt});
      histos.add("h2dGenXiMinus", "h2dGenXiMinus", kTH2D, {axisCentrality, axisPt});
      histos.add("h2dGenXiPlus", "h2dGenXiPlus", kTH2D, {axisCentrality, axisPt});
      histos.add("h2dGenOmegaMinus", "h2dGenOmegaMinus", kTH2D, {axisCentrality, axisPt});
      histos.add("h2dGenOmegaPlus", "h2dGenOmegaPlus", kTH2D, {axisCentrality, axisPt});
    }

    // inspect histogram sizes, please
    histos.print();
  }

  template <typename TV0, typename TCollision>
  uint64_t computeReconstructionBitmap(TV0 v0, TCollision collision, float rapidityLambda, float rapidityK0Short, float /*pT*/)
  // precalculate this information so that a check is one mask operation, not many
  {
    uint64_t bitMap = 0;
    // Base topological variables
    if (v0.v0radius() > v0radius)
      bitset(bitMap, selRadius);
    if (v0.v0radius() < v0radiusMax)
      bitset(bitMap, selRadiusMax);
    if (TMath::Abs(v0.dcapostopv()) > dcapostopv)
      bitset(bitMap, selDCAPosToPV);
    if (TMath::Abs(v0.dcanegtopv()) > dcanegtopv)
      bitset(bitMap, selDCANegToPV);
    if (v0.v0cosPA() > v0cospa)
      bitset(bitMap, selCosPA);
    if (v0.dcaV0daughters() < dcav0dau)
      bitset(bitMap, selDCAV0Dau);

    // rapidity
    if (TMath::Abs(rapidityLambda) < rapidityCut)
      bitset(bitMap, selLambdaRapidity);
    if (TMath::Abs(rapidityK0Short) < rapidityCut)
      bitset(bitMap, selK0ShortRapidity);

    auto posTrackExtra = v0.template posTrackExtra_as<dauTracks>();
    auto negTrackExtra = v0.template negTrackExtra_as<dauTracks>();

    // ITS quality flags
    if (posTrackExtra.itsNCls() >= minITSclusters)
      bitset(bitMap, selPosGoodITSTrack);
    if (negTrackExtra.itsNCls() >= minITSclusters)
      bitset(bitMap, selNegGoodITSTrack);

    // TPC quality flags
    if (posTrackExtra.tpcCrossedRows() >= minTPCrows)
      bitset(bitMap, selPosGoodTPCTrack);
    if (negTrackExtra.tpcCrossedRows() >= minTPCrows)
      bitset(bitMap, selNegGoodTPCTrack);

    // TPC PID
    if (fabs(posTrackExtra.tpcNSigmaPi()) < TpcPidNsigmaCut)
      bitset(bitMap, selTPCPIDPositivePion);
    if (fabs(posTrackExtra.tpcNSigmaPr()) < TpcPidNsigmaCut)
      bitset(bitMap, selTPCPIDPositiveProton);
    if (fabs(negTrackExtra.tpcNSigmaPi()) < TpcPidNsigmaCut)
      bitset(bitMap, selTPCPIDNegativePion);
    if (fabs(negTrackExtra.tpcNSigmaPr()) < TpcPidNsigmaCut)
      bitset(bitMap, selTPCPIDNegativeProton);

    // TOF PID in DeltaT
    // Positive track
    if (fabs(v0.posTOFDeltaTLaPr()) < maxDeltaTimeProton)
      bitset(bitMap, selTOFDeltaTPositiveProtonLambda);
    if (fabs(v0.posTOFDeltaTLaPi()) < maxDeltaTimePion)
      bitset(bitMap, selTOFDeltaTPositivePionLambda);
    if (fabs(v0.posTOFDeltaTK0Pi()) < maxDeltaTimePion)
      bitset(bitMap, selTOFDeltaTPositivePionK0Short);
    // Negative track
    if (fabs(v0.negTOFDeltaTLaPr()) < maxDeltaTimeProton)
      bitset(bitMap, selTOFDeltaTNegativeProtonLambda);
    if (fabs(v0.negTOFDeltaTLaPi()) < maxDeltaTimePion)
      bitset(bitMap, selTOFDeltaTNegativePionLambda);
    if (fabs(v0.negTOFDeltaTK0Pi()) < maxDeltaTimePion)
      bitset(bitMap, selTOFDeltaTNegativePionK0Short);

    // TOF PID in NSigma
    // Positive track
    if (fabs(v0.tofNSigmaLaPr()) < TofPidNsigmaCutLaPr)
      bitset(bitMap, selTOFNSigmaPositiveProtonLambda);
    if (fabs(v0.tofNSigmaALaPi()) < TofPidNsigmaCutLaPi)
      bitset(bitMap, selTOFNSigmaPositivePionLambda);
    if (fabs(v0.tofNSigmaK0PiPlus()) < TofPidNsigmaCutK0Pi)
      bitset(bitMap, selTOFNSigmaPositivePionK0Short);
    // Negative track
    if (fabs(v0.tofNSigmaALaPr()) < TofPidNsigmaCutLaPr)
      bitset(bitMap, selTOFNSigmaNegativeProtonLambda);
    if (fabs(v0.tofNSigmaLaPi()) < TofPidNsigmaCutLaPi)
      bitset(bitMap, selTOFNSigmaNegativePionLambda);
    if (fabs(v0.tofNSigmaK0PiMinus()) < TofPidNsigmaCutK0Pi)
      bitset(bitMap, selTOFNSigmaNegativePionK0Short);

    // ITS only tag
    if (posTrackExtra.tpcCrossedRows() < 1)
      bitset(bitMap, selPosItsOnly);
    if (negTrackExtra.tpcCrossedRows() < 1)
      bitset(bitMap, selNegItsOnly);

    // TPC only tag
    if (posTrackExtra.detectorMap() != o2::aod::track::TPC)
      bitset(bitMap, selPosNotTPCOnly);
    if (negTrackExtra.detectorMap() != o2::aod::track::TPC)
      bitset(bitMap, selNegNotTPCOnly);

    // proper lifetime
    if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0 < lifetimecut->get("lifetimecutLambda"))
      bitset(bitMap, selLambdaCTau);
    if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short < lifetimecut->get("lifetimecutK0S"))
      bitset(bitMap, selK0ShortCTau);

    // armenteros
    if (v0.qtarm() * armPodCut > TMath::Abs(v0.alpha()) || armPodCut < 1e-4)
      bitset(bitMap, selK0ShortArmenteros);

    return bitMap;
  }

  template <typename TV0>
  uint64_t computeMCAssociation(TV0 v0)
  // precalculate this information so that a check is one mask operation, not many
  {
    uint64_t bitMap = 0;
    // check for specific particle species

    if (v0.pdgCode() == 310 && v0.pdgCodePositive() == 211 && v0.pdgCodeNegative() == -211) {
      bitset(bitMap, selConsiderK0Short);
      if (v0.isPhysicalPrimary())
        bitset(bitMap, selPhysPrimK0Short);
    }
    if (v0.pdgCode() == 3122 && v0.pdgCodePositive() == 2212 && v0.pdgCodeNegative() == -211) {
      bitset(bitMap, selConsiderLambda);
      if (v0.isPhysicalPrimary())
        bitset(bitMap, selPhysPrimLambda);
    }
    if (v0.pdgCode() == -3122 && v0.pdgCodePositive() == 211 && v0.pdgCodeNegative() == -2212) {
      bitset(bitMap, selConsiderAntiLambda);
      if (v0.isPhysicalPrimary())
        bitset(bitMap, selPhysPrimAntiLambda);
    }
    return bitMap;
  }

  bool verifyMask(uint64_t bitmap, uint64_t mask)
  {
    return (bitmap & mask) == mask;
  }

  uint computeITSclusBitmap(uint8_t itsClusMap)
  // Focus on the 12 dominant ITS cluster configurations
  {
    uint bitMap = 0;

    if (verifyMask(itsClusMap, ((uint8_t(1) << 0) | (uint8_t(1) << 1) | (uint8_t(1) << 2) | (uint8_t(1) << 3) | (uint8_t(1) << 4) | (uint8_t(1) << 5) | (uint8_t(1) << 6)))) {
      // ITS :    IB         OB
      // ITS : L0 L1 L2 L3 L4 L5 L6
      // ITS : x  x  x  x  x  x  x
      bitMap = 12;
    } else if (verifyMask(itsClusMap, ((uint8_t(1) << 1) | (uint8_t(1) << 2) | (uint8_t(1) << 3) | (uint8_t(1) << 4) | (uint8_t(1) << 5) | (uint8_t(1) << 6)))) {
      // ITS :    IB         OB
      // ITS : L0 L1 L2 L3 L4 L5 L6
      // ITS :    x  x  x  x  x  x
      bitMap = 11;
    } else if (verifyMask(itsClusMap, ((uint8_t(1) << 2) | (uint8_t(1) << 3) | (uint8_t(1) << 4) | (uint8_t(1) << 5) | (uint8_t(1) << 6)))) {
      // ITS :    IB         OB
      // ITS : L0 L1 L2 L3 L4 L5 L6
      // ITS :       x  x  x  x  x
      bitMap = 10;
    } else if (verifyMask(itsClusMap, ((uint8_t(1) << 3) | (uint8_t(1) << 4) | (uint8_t(1) << 5) | (uint8_t(1) << 6)))) {
      // ITS :    IB         OB
      // ITS : L0 L1 L2 L3 L4 L5 L6
      // ITS :          x  x  x  x
      bitMap = 9;
    } else if (verifyMask(itsClusMap, ((uint8_t(1) << 4) | (uint8_t(1) << 5) | (uint8_t(1) << 6)))) {
      // ITS :    IB         OB
      // ITS : L0 L1 L2 L3 L4 L5 L6
      // ITS :             x  x  x
      bitMap = 8;
    } else if (verifyMask(itsClusMap, ((uint8_t(1) << 5) | (uint8_t(1) << 6)))) {
      // ITS :    IB         OB
      // ITS : L0 L1 L2 L3 L4 L5 L6
      // ITS :                x  x
      bitMap = 7;
    } else if (verifyMask(itsClusMap, ((uint8_t(1) << 0) | (uint8_t(1) << 1) | (uint8_t(1) << 2) | (uint8_t(1) << 3) | (uint8_t(1) << 4) | (uint8_t(1) << 5)))) {
      // ITS :    IB         OB
      // ITS : L0 L1 L2 L3 L4 L5 L6
      // ITS : x  x  x  x  x  x
      bitMap = 6;
    } else if (verifyMask(itsClusMap, ((uint8_t(1) << 1) | (uint8_t(1) << 2) | (uint8_t(1) << 3) | (uint8_t(1) << 4) | (uint8_t(1) << 5)))) {
      // ITS :    IB         OB
      // ITS : L0 L1 L2 L3 L4 L5 L6
      // ITS :    x  x  x  x  x
      bitMap = 5;
    } else if (verifyMask(itsClusMap, ((uint8_t(1) << 2) | (uint8_t(1) << 3) | (uint8_t(1) << 4) | (uint8_t(1) << 5)))) {
      // ITS :    IB         OB
      // ITS : L0 L1 L2 L3 L4 L5 L6
      // ITS :       x  x  x  x
      bitMap = 4;
    } else if (verifyMask(itsClusMap, ((uint8_t(1) << 0) | (uint8_t(1) << 1) | (uint8_t(1) << 2) | (uint8_t(1) << 3) | (uint8_t(1) << 4)))) {
      // ITS :    IB         OB
      // ITS : L0 L1 L2 L3 L4 L5 L6
      // ITS : x  x  x  x  x
      bitMap = 3;
    } else if (verifyMask(itsClusMap, ((uint8_t(1) << 1) | (uint8_t(1) << 2) | (uint8_t(1) << 3) | (uint8_t(1) << 4)))) {
      // ITS :    IB         OB
      // ITS : L0 L1 L2 L3 L4 L5 L6
      // ITS :    x  x  x  x
      bitMap = 2;
    } else if (verifyMask(itsClusMap, ((uint8_t(1) << 0) | (uint8_t(1) << 1) | (uint8_t(1) << 2) | (uint8_t(1) << 3)))) {
      // ITS :    IB         OB
      // ITS : L0 L1 L2 L3 L4 L5 L6
      // ITS : x  x  x  x
      bitMap = 1;
    } else {
      // ITS : other configurations
      bitMap = 0;
    }

    return bitMap;
  }

  uint computeDetBitmap(uint8_t detMap)
  // Focus on the 4 dominant track configurations :
  //  ITS-TPC
  //  ITS-TPC-TRD
  //  ITS-TPC-TOF
  //  ITS-TPC-TRD-TOF
  {
    uint bitMap = 0;

    if (verifyMask(detMap, (o2::aod::track::ITS | o2::aod::track::TPC | o2::aod::track::TRD | o2::aod::track::TOF))) {
      // ITS-TPC-TRD-TOF
      bitMap = 3;
    } else if (verifyMask(detMap, (o2::aod::track::ITS | o2::aod::track::TPC | o2::aod::track::TOF))) {
      // ITS-TPC-TOF
      bitMap = 2;
    } else if (verifyMask(detMap, (o2::aod::track::ITS | o2::aod::track::TPC | o2::aod::track::TRD))) {
      // ITS-TPC-TRD
      bitMap = 1;
    } else if (verifyMask(detMap, (o2::aod::track::ITS | o2::aod::track::TPC))) {
      // ITS-TPC
      bitMap = 0;
    }

    return bitMap;
  }

  template <typename TV0>
  void analyseCandidate(TV0 v0, float pt, float centrality, uint64_t selMap)
  // precalculate this information so that a check is one mask operation, not many
  {
    auto posTrackExtra = v0.template posTrackExtra_as<dauTracks>();
    auto negTrackExtra = v0.template negTrackExtra_as<dauTracks>();

    uint posDetMap = computeDetBitmap(posTrackExtra.detectorMap());
    uint posITSclusMap = computeITSclusBitmap(posTrackExtra.itsClusterMap());
    uint negDetMap = computeDetBitmap(negTrackExtra.detectorMap());
    uint negITSclusMap = computeITSclusBitmap(negTrackExtra.itsClusterMap());

    // __________________________________________
    // fill with no selection if plain QA requested
    if (doPlainTopoQA) {
      histos.fill(HIST("hPosDCAToPV"), v0.dcapostopv());
      histos.fill(HIST("hNegDCAToPV"), v0.dcanegtopv());
      histos.fill(HIST("hDCADaughters"), v0.dcaV0daughters());
      histos.fill(HIST("hPointingAngle"), TMath::ACos(v0.v0cosPA()));
      histos.fill(HIST("hV0Radius"), v0.v0radius());
      histos.fill(HIST("h2dPositiveITSvsTPCpts"), posTrackExtra.tpcCrossedRows(), posTrackExtra.itsNCls());
      histos.fill(HIST("h2dNegativeITSvsTPCpts"), negTrackExtra.tpcCrossedRows(), negTrackExtra.itsNCls());
    }

    // __________________________________________
    // main analysis
    if (verifyMask(selMap, maskSelectionK0Short) && analyseK0Short) {
      histos.fill(HIST("GeneralQA/h2dArmenterosSelected"), v0.alpha(), v0.qtarm()); // cross-check
      histos.fill(HIST("h3dMassK0Short"), centrality, pt, v0.mK0Short());
      histos.fill(HIST("hMassK0Short"), v0.mK0Short());
      if (doPlainTopoQA) {
        histos.fill(HIST("K0Short/hPosDCAToPV"), v0.dcapostopv());
        histos.fill(HIST("K0Short/hNegDCAToPV"), v0.dcanegtopv());
        histos.fill(HIST("K0Short/hDCADaughters"), v0.dcaV0daughters());
        histos.fill(HIST("K0Short/hPointingAngle"), TMath::ACos(v0.v0cosPA()));
        histos.fill(HIST("K0Short/hV0Radius"), v0.v0radius());
        histos.fill(HIST("K0Short/h2dPositiveITSvsTPCpts"), posTrackExtra.tpcCrossedRows(), posTrackExtra.itsNCls());
        histos.fill(HIST("K0Short/h2dNegativeITSvsTPCpts"), negTrackExtra.tpcCrossedRows(), negTrackExtra.itsNCls());
      }
      if (doDetectPropQA == 1) {
        histos.fill(HIST("K0Short/h6dDetectPropVsCentrality"), centrality, posDetMap, posITSclusMap, negDetMap, negITSclusMap, pt);
        histos.fill(HIST("K0Short/h4dPosDetectPropVsCentrality"), centrality, posTrackExtra.detectorMap(), posTrackExtra.itsClusterMap(), pt);
        histos.fill(HIST("K0Short/h4dNegDetectPropVsCentrality"), centrality, negTrackExtra.detectorMap(), negTrackExtra.itsClusterMap(), pt);
      }
      if (doDetectPropQA == 2) {
        histos.fill(HIST("K0Short/h7dPosDetectPropVsCentrality"), centrality, posDetMap, posITSclusMap, negDetMap, negITSclusMap, pt, v0.mK0Short());
        histos.fill(HIST("K0Short/h5dPosDetectPropVsCentrality"), centrality, posTrackExtra.detectorMap(), posTrackExtra.itsClusterMap(), pt, v0.mK0Short());
        histos.fill(HIST("K0Short/h5dNegDetectPropVsCentrality"), centrality, negTrackExtra.detectorMap(), negTrackExtra.itsClusterMap(), pt, v0.mK0Short());
      }
      if (doTPCQA) {
        histos.fill(HIST("K0Short/h3dPosNsigmaTPC"), centrality, pt, posTrackExtra.tpcNSigmaPi());
        histos.fill(HIST("K0Short/h3dNegNsigmaTPC"), centrality, pt, negTrackExtra.tpcNSigmaPi());
        histos.fill(HIST("K0Short/h3dPosTPCsignal"), centrality, pt, posTrackExtra.tpcSignal());
        histos.fill(HIST("K0Short/h3dNegTPCsignal"), centrality, pt, negTrackExtra.tpcSignal());
        histos.fill(HIST("K0Short/h3dPosNsigmaTPCvsTrackPtot"), centrality, v0.positivept() * TMath::CosH(v0.positiveeta()), posTrackExtra.tpcNSigmaPi());
        histos.fill(HIST("K0Short/h3dNegNsigmaTPCvsTrackPtot"), centrality, v0.negativept() * TMath::CosH(v0.negativeeta()), negTrackExtra.tpcNSigmaPi());
        histos.fill(HIST("K0Short/h3dPosTPCsignalVsTrackPtot"), centrality, v0.positivept() * TMath::CosH(v0.positiveeta()), posTrackExtra.tpcSignal());
        histos.fill(HIST("K0Short/h3dNegTPCsignalVsTrackPtot"), centrality, v0.negativept() * TMath::CosH(v0.negativeeta()), negTrackExtra.tpcSignal());
        histos.fill(HIST("K0Short/h3dPosNsigmaTPCvsTrackPt"), centrality, v0.positivept(), posTrackExtra.tpcNSigmaPi());
        histos.fill(HIST("K0Short/h3dNegNsigmaTPCvsTrackPt"), centrality, v0.negativept(), negTrackExtra.tpcNSigmaPi());
        histos.fill(HIST("K0Short/h3dPosTPCsignalVsTrackPt"), centrality, v0.positivept(), posTrackExtra.tpcSignal());
        histos.fill(HIST("K0Short/h3dNegTPCsignalVsTrackPt"), centrality, v0.negativept(), negTrackExtra.tpcSignal());
      }
      if (doTOFQA) {
        histos.fill(HIST("K0Short/h3dPosTOFdeltaT"), centrality, pt, v0.posTOFDeltaTK0Pi());
        histos.fill(HIST("K0Short/h3dNegTOFdeltaT"), centrality, pt, v0.negTOFDeltaTK0Pi());
        histos.fill(HIST("K0Short/h3dPosTOFdeltaTvsTrackPtot"), centrality, v0.positivept() * TMath::CosH(v0.positiveeta()), v0.posTOFDeltaTK0Pi());
        histos.fill(HIST("K0Short/h3dNegTOFdeltaTvsTrackPtot"), centrality, v0.negativept() * TMath::CosH(v0.negativeeta()), v0.negTOFDeltaTK0Pi());
        histos.fill(HIST("K0Short/h3dPosTOFdeltaTvsTrackPt"), centrality, v0.positivept(), v0.posTOFDeltaTK0Pi());
        histos.fill(HIST("K0Short/h3dNegTOFdeltaTvsTrackPt"), centrality, v0.negativept(), v0.negTOFDeltaTK0Pi());
      }
    }
    if (verifyMask(selMap, maskSelectionLambda) && analyseLambda) {
      histos.fill(HIST("h3dMassLambda"), centrality, pt, v0.mLambda());
      if (doPlainTopoQA) {
        histos.fill(HIST("Lambda/hPosDCAToPV"), v0.dcapostopv());
        histos.fill(HIST("Lambda/hNegDCAToPV"), v0.dcanegtopv());
        histos.fill(HIST("Lambda/hDCADaughters"), v0.dcaV0daughters());
        histos.fill(HIST("Lambda/hPointingAngle"), TMath::ACos(v0.v0cosPA()));
        histos.fill(HIST("Lambda/hV0Radius"), v0.v0radius());
        histos.fill(HIST("Lambda/h2dPositiveITSvsTPCpts"), posTrackExtra.tpcCrossedRows(), posTrackExtra.itsNCls());
        histos.fill(HIST("Lambda/h2dNegativeITSvsTPCpts"), negTrackExtra.tpcCrossedRows(), negTrackExtra.itsNCls());
      }
      if (doDetectPropQA == 1) {
        histos.fill(HIST("Lambda/h6dDetectPropVsCentrality"), centrality, posDetMap, posITSclusMap, negDetMap, negITSclusMap, pt);
        histos.fill(HIST("Lambda/h4dPosDetectPropVsCentrality"), centrality, posTrackExtra.detectorMap(), posTrackExtra.itsClusterMap(), pt);
        histos.fill(HIST("Lambda/h4dNegDetectPropVsCentrality"), centrality, negTrackExtra.detectorMap(), negTrackExtra.itsClusterMap(), pt);
      }
      if (doDetectPropQA == 2) {
        histos.fill(HIST("Lambda/h7dDetectPropVsCentrality"), centrality, posDetMap, posITSclusMap, negDetMap, negITSclusMap, pt, v0.mLambda());
        histos.fill(HIST("Lambda/h5dPosDetectPropVsCentrality"), centrality, posTrackExtra.detectorMap(), posTrackExtra.itsClusterMap(), pt, v0.mLambda());
        histos.fill(HIST("Lambda/h5dNegDetectPropVsCentrality"), centrality, negTrackExtra.detectorMap(), negTrackExtra.itsClusterMap(), pt, v0.mLambda());
      }
      if (doTPCQA) {
        histos.fill(HIST("Lambda/h3dPosNsigmaTPC"), centrality, pt, posTrackExtra.tpcNSigmaPr());
        histos.fill(HIST("Lambda/h3dNegNsigmaTPC"), centrality, pt, negTrackExtra.tpcNSigmaPi());
        histos.fill(HIST("Lambda/h3dPosTPCsignal"), centrality, pt, posTrackExtra.tpcSignal());
        histos.fill(HIST("Lambda/h3dNegTPCsignal"), centrality, pt, negTrackExtra.tpcSignal());
        histos.fill(HIST("Lambda/h3dPosNsigmaTPCvsTrackPtot"), centrality, v0.positivept() * TMath::CosH(v0.positiveeta()), posTrackExtra.tpcNSigmaPr());
        histos.fill(HIST("Lambda/h3dNegNsigmaTPCvsTrackPtot"), centrality, v0.negativept() * TMath::CosH(v0.negativeeta()), negTrackExtra.tpcNSigmaPi());
        histos.fill(HIST("Lambda/h3dPosTPCsignalVsTrackPtot"), centrality, v0.positivept() * TMath::CosH(v0.positiveeta()), posTrackExtra.tpcSignal());
        histos.fill(HIST("Lambda/h3dNegTPCsignalVsTrackPtot"), centrality, v0.negativept() * TMath::CosH(v0.negativeeta()), negTrackExtra.tpcSignal());
        histos.fill(HIST("Lambda/h3dPosNsigmaTPCvsTrackPt"), centrality, v0.positivept(), posTrackExtra.tpcNSigmaPr());
        histos.fill(HIST("Lambda/h3dNegNsigmaTPCvsTrackPt"), centrality, v0.negativept(), negTrackExtra.tpcNSigmaPi());
        histos.fill(HIST("Lambda/h3dPosTPCsignalVsTrackPt"), centrality, v0.positivept(), posTrackExtra.tpcSignal());
        histos.fill(HIST("Lambda/h3dNegTPCsignalVsTrackPt"), centrality, v0.negativept(), negTrackExtra.tpcSignal());
      }
      if (doTOFQA) {
        histos.fill(HIST("Lambda/h3dPosTOFdeltaT"), centrality, pt, v0.posTOFDeltaTLaPr());
        histos.fill(HIST("Lambda/h3dNegTOFdeltaT"), centrality, pt, v0.negTOFDeltaTLaPi());
        histos.fill(HIST("Lambda/h3dPosTOFdeltaTvsTrackPtot"), centrality, v0.positivept() * TMath::CosH(v0.positiveeta()), v0.posTOFDeltaTLaPr());
        histos.fill(HIST("Lambda/h3dNegTOFdeltaTvsTrackPtot"), centrality, v0.negativept() * TMath::CosH(v0.negativeeta()), v0.negTOFDeltaTLaPi());
        histos.fill(HIST("Lambda/h3dPosTOFdeltaTvsTrackPt"), centrality, v0.positivept(), v0.posTOFDeltaTLaPr());
        histos.fill(HIST("Lambda/h3dNegTOFdeltaTvsTrackPt"), centrality, v0.negativept(), v0.negTOFDeltaTLaPi());
      }
    }
    if (verifyMask(selMap, maskSelectionAntiLambda) && analyseAntiLambda) {
      histos.fill(HIST("h3dMassAntiLambda"), centrality, pt, v0.mAntiLambda());
      if (doPlainTopoQA) {
        histos.fill(HIST("AntiLambda/hPosDCAToPV"), v0.dcapostopv());
        histos.fill(HIST("AntiLambda/hNegDCAToPV"), v0.dcanegtopv());
        histos.fill(HIST("AntiLambda/hDCADaughters"), v0.dcaV0daughters());
        histos.fill(HIST("AntiLambda/hPointingAngle"), TMath::ACos(v0.v0cosPA()));
        histos.fill(HIST("AntiLambda/hV0Radius"), v0.v0radius());
        histos.fill(HIST("AntiLambda/h2dPositiveITSvsTPCpts"), posTrackExtra.tpcCrossedRows(), posTrackExtra.itsNCls());
        histos.fill(HIST("AntiLambda/h2dNegativeITSvsTPCpts"), negTrackExtra.tpcCrossedRows(), negTrackExtra.itsNCls());
      }
      if (doDetectPropQA == 1) {
        histos.fill(HIST("AntiLambda/h6dDetectPropVsCentrality"), centrality, posDetMap, posITSclusMap, negDetMap, negITSclusMap, pt);
        histos.fill(HIST("AntiLambda/h4dPosDetectPropVsCentrality"), centrality, posTrackExtra.detectorMap(), posTrackExtra.itsClusterMap(), pt);
        histos.fill(HIST("AntiLambda/h4dNegDetectPropVsCentrality"), centrality, negTrackExtra.detectorMap(), negTrackExtra.itsClusterMap(), pt);
      }
      if (doDetectPropQA == 2) {
        histos.fill(HIST("AntiLambda/h7dDetectPropVsCentrality"), centrality, posDetMap, posITSclusMap, negDetMap, negITSclusMap, pt, v0.mAntiLambda());
        histos.fill(HIST("AntiLambda/h5dPosDetectPropVsCentrality"), centrality, posTrackExtra.detectorMap(), posTrackExtra.itsClusterMap(), pt, v0.mAntiLambda());
        histos.fill(HIST("AntiLambda/h5dNegDetectPropVsCentrality"), centrality, negTrackExtra.detectorMap(), negTrackExtra.itsClusterMap(), pt, v0.mAntiLambda());
      }
      if (doTPCQA) {
        histos.fill(HIST("AntiLambda/h3dPosNsigmaTPC"), centrality, pt, posTrackExtra.tpcNSigmaPi());
        histos.fill(HIST("AntiLambda/h3dNegNsigmaTPC"), centrality, pt, negTrackExtra.tpcNSigmaPr());
        histos.fill(HIST("AntiLambda/h3dPosTPCsignal"), centrality, pt, posTrackExtra.tpcSignal());
        histos.fill(HIST("AntiLambda/h3dNegTPCsignal"), centrality, pt, negTrackExtra.tpcSignal());
        histos.fill(HIST("AntiLambda/h3dPosNsigmaTPCvsTrackPtot"), centrality, v0.positivept() * TMath::CosH(v0.positiveeta()), posTrackExtra.tpcNSigmaPi());
        histos.fill(HIST("AntiLambda/h3dNegNsigmaTPCvsTrackPtot"), centrality, v0.negativept() * TMath::CosH(v0.negativeeta()), negTrackExtra.tpcNSigmaPr());
        histos.fill(HIST("AntiLambda/h3dPosTPCsignalVsTrackPtot"), centrality, v0.positivept() * TMath::CosH(v0.positiveeta()), posTrackExtra.tpcSignal());
        histos.fill(HIST("AntiLambda/h3dNegTPCsignalVsTrackPtot"), centrality, v0.negativept() * TMath::CosH(v0.negativeeta()), negTrackExtra.tpcSignal());
        histos.fill(HIST("AntiLambda/h3dPosNsigmaTPCvsTrackPt"), centrality, v0.positivept(), posTrackExtra.tpcNSigmaPi());
        histos.fill(HIST("AntiLambda/h3dNegNsigmaTPCvsTrackPt"), centrality, v0.negativept(), negTrackExtra.tpcNSigmaPr());
        histos.fill(HIST("AntiLambda/h3dPosTPCsignalVsTrackPt"), centrality, v0.positivept(), posTrackExtra.tpcSignal());
        histos.fill(HIST("AntiLambda/h3dNegTPCsignalVsTrackPt"), centrality, v0.negativept(), negTrackExtra.tpcSignal());
      }
      if (doTOFQA) {
        histos.fill(HIST("AntiLambda/h3dPosTOFdeltaT"), centrality, pt, v0.posTOFDeltaTLaPi());
        histos.fill(HIST("AntiLambda/h3dNegTOFdeltaT"), centrality, pt, v0.negTOFDeltaTLaPr());
        histos.fill(HIST("AntiLambda/h3dPosTOFdeltaTvsTrackPtot"), centrality, v0.positivept() * TMath::CosH(v0.positiveeta()), v0.posTOFDeltaTLaPi());
        histos.fill(HIST("AntiLambda/h3dNegTOFdeltaTvsTrackPtot"), centrality, v0.negativept() * TMath::CosH(v0.negativeeta()), v0.negTOFDeltaTLaPr());
        histos.fill(HIST("AntiLambda/h3dPosTOFdeltaTvsTrackPt"), centrality, v0.positivept(), v0.posTOFDeltaTLaPi());
        histos.fill(HIST("AntiLambda/h3dNegTOFdeltaTvsTrackPt"), centrality, v0.negativept(), v0.negTOFDeltaTLaPr());
      }
    }

    // __________________________________________
    // do systematics / qa plots
    if (doCompleteTopoQA) {
      if (analyseK0Short) {
        if (verifyMask(selMap, maskTopoNoV0Radius | maskK0ShortSpecific))
          histos.fill(HIST("K0Short/h4dV0Radius"), centrality, pt, v0.mK0Short(), v0.v0radius());
        if (verifyMask(selMap, maskTopoNoDCAPosToPV | maskK0ShortSpecific))
          histos.fill(HIST("K0Short/h4dPosDCAToPV"), centrality, pt, v0.mK0Short(), TMath::Abs(v0.dcapostopv()));
        if (verifyMask(selMap, maskTopoNoDCANegToPV | maskK0ShortSpecific))
          histos.fill(HIST("K0Short/h4dNegDCAToPV"), centrality, pt, v0.mK0Short(), TMath::Abs(v0.dcanegtopv()));
        if (verifyMask(selMap, maskTopoNoCosPA | maskK0ShortSpecific))
          histos.fill(HIST("K0Short/h4dPointingAngle"), centrality, pt, v0.mK0Short(), TMath::ACos(v0.v0cosPA()));
        if (verifyMask(selMap, maskTopoNoDCAV0Dau | maskK0ShortSpecific))
          histos.fill(HIST("K0Short/h4dDCADaughters"), centrality, pt, v0.mK0Short(), v0.dcaV0daughters());
      }

      if (analyseLambda) {
        if (verifyMask(selMap, maskTopoNoV0Radius | maskLambdaSpecific))
          histos.fill(HIST("Lambda/h4dV0Radius"), centrality, pt, v0.mLambda(), v0.v0radius());
        if (verifyMask(selMap, maskTopoNoDCAPosToPV | maskLambdaSpecific))
          histos.fill(HIST("Lambda/h4dPosDCAToPV"), centrality, pt, v0.mLambda(), TMath::Abs(v0.dcapostopv()));
        if (verifyMask(selMap, maskTopoNoDCANegToPV | maskLambdaSpecific))
          histos.fill(HIST("Lambda/h4dNegDCAToPV"), centrality, pt, v0.mLambda(), TMath::Abs(v0.dcanegtopv()));
        if (verifyMask(selMap, maskTopoNoCosPA | maskLambdaSpecific))
          histos.fill(HIST("Lambda/h4dPointingAngle"), centrality, pt, v0.mLambda(), TMath::ACos(v0.v0cosPA()));
        if (verifyMask(selMap, maskTopoNoDCAV0Dau | maskLambdaSpecific))
          histos.fill(HIST("Lambda/h4dDCADaughters"), centrality, pt, v0.mLambda(), v0.dcaV0daughters());
      }
      if (analyseAntiLambda) {
        if (verifyMask(selMap, maskTopoNoV0Radius | maskAntiLambdaSpecific))
          histos.fill(HIST("AntiLambda/h4dV0Radius"), centrality, pt, v0.mAntiLambda(), v0.v0radius());
        if (verifyMask(selMap, maskTopoNoDCAPosToPV | maskAntiLambdaSpecific))
          histos.fill(HIST("AntiLambda/h4dPosDCAToPV"), centrality, pt, v0.mAntiLambda(), TMath::Abs(v0.dcapostopv()));
        if (verifyMask(selMap, maskTopoNoDCANegToPV | maskAntiLambdaSpecific))
          histos.fill(HIST("AntiLambda/h4dNegDCAToPV"), centrality, pt, v0.mAntiLambda(), TMath::Abs(v0.dcanegtopv()));
        if (verifyMask(selMap, maskTopoNoCosPA | maskAntiLambdaSpecific))
          histos.fill(HIST("AntiLambda/h4dPointingAngle"), centrality, pt, v0.mAntiLambda(), TMath::ACos(v0.v0cosPA()));
        if (verifyMask(selMap, maskTopoNoDCAV0Dau | maskAntiLambdaSpecific))
          histos.fill(HIST("AntiLambda/h4dDCADaughters"), centrality, pt, v0.mAntiLambda(), v0.dcaV0daughters());
      }
    } // end systematics / qa
  }

  template <typename TV0>
  void analyseCollisionAssociation(TV0 /*v0*/, float pt, int mcNch, bool correctAssociation, uint64_t selMap)
  // analyse collision association
  {
    // __________________________________________
    // main analysis
    if (verifyMask(selMap, maskSelectionK0Short) && analyseK0Short) {
      histos.fill(HIST("K0Short/h2dPtVsNch"), mcNch, pt);
      if (!correctAssociation)
        histos.fill(HIST("K0Short/h2dPtVsNch_BadCollAssig"), mcNch, pt);
    }
    if (verifyMask(selMap, maskSelectionLambda) && analyseLambda) {
      histos.fill(HIST("Lambda/h2dPtVsNch"), mcNch, pt);
      if (!correctAssociation)
        histos.fill(HIST("Lambda/h2dPtVsNch_BadCollAssig"), mcNch, pt);
    }
    if (verifyMask(selMap, maskSelectionAntiLambda) && analyseAntiLambda) {
      histos.fill(HIST("AntiLambda/h2dPtVsNch"), mcNch, pt);
      if (!correctAssociation)
        histos.fill(HIST("AntiLambda/h2dPtVsNch_BadCollAssig"), mcNch, pt);
    }
  }

  template <typename TV0>
  void fillFeeddownMatrix(TV0 v0, float pt, float centrality, uint64_t selMap)
  // fill feeddown matrix for Lambdas or AntiLambdas
  // fixme: a potential improvement would be to consider mass windows for the l/al
  {
    if (!v0.has_motherMCPart())
      return; // does not have mother particle in record, skip

    auto v0mother = v0.motherMCPart();
    float rapidityXi = RecoDecay::y(std::array{v0mother.px(), v0mother.py(), v0mother.pz()}, o2::constants::physics::MassXiMinus);
    if (fabs(rapidityXi) > 0.5f)
      return; // not a valid mother rapidity (PDG selection is later)

    // __________________________________________
    if (verifyMask(selMap, secondaryMaskSelectionLambda) && analyseLambda) {
      if (v0mother.pdgCode() == 3312 && v0mother.isPhysicalPrimary())
        histos.fill(HIST("h3dLambdaFeeddown"), centrality, pt, std::hypot(v0mother.px(), v0mother.py()));
    }
    if (verifyMask(selMap, secondaryMaskSelectionAntiLambda) && analyseAntiLambda) {
      if (v0mother.pdgCode() == -3312 && v0mother.isPhysicalPrimary())
        histos.fill(HIST("h3dAntiLambdaFeeddown"), centrality, pt, std::hypot(v0mother.px(), v0mother.py()));
    }
  }

  // ______________________________________________________
  // Real data processing - no MC subscription
  void processRealData(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraRawCents, aod::StraEvSels>::iterator const& collision, v0Candidates const& fullV0s, dauTracks const&)
  {
    histos.fill(HIST("hEventSelection"), 0. /* all collisions */);
    if (!collision.sel8()) {
      return;
    }
    histos.fill(HIST("hEventSelection"), 1 /* sel8 collisions */);

    if (std::abs(collision.posZ()) > 10.f) {
      return;
    }
    histos.fill(HIST("hEventSelection"), 2 /* vertex-Z selected */);

    if (rejectITSROFBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return;
    }
    histos.fill(HIST("hEventSelection"), 3 /* Not at ITS ROF border */);

    if (rejectTFBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return;
    }
    histos.fill(HIST("hEventSelection"), 4 /* Not at TF border */);

    float centrality = collision.centFT0C();
    if (qaCentrality) {
      auto hRawCentrality = histos.get<TH1>(HIST("hRawCentrality"));
      centrality = hRawCentrality->GetBinContent(hRawCentrality->FindBin(collision.multFT0C()));
    }

    histos.fill(HIST("hEventCentrality"), centrality);

    histos.fill(HIST("hCentralityVsNch"), centrality, collision.multNTracksPVeta1());

    // __________________________________________
    // perform main analysis
    for (auto& v0 : fullV0s) {
      if (std::abs(v0.negativeeta()) > daughterEtaCut || std::abs(v0.positiveeta()) > daughterEtaCut)
        continue; // remove acceptance that's badly reproduced by MC / superfluous in future

      // fill AP plot for all V0s
      histos.fill(HIST("GeneralQA/h2dArmenterosAll"), v0.alpha(), v0.qtarm());

      uint64_t selMap = computeReconstructionBitmap(v0, collision, v0.yLambda(), v0.yK0Short(), v0.pt());

      // consider for histograms for all species
      selMap = selMap | (uint64_t(1) << selConsiderK0Short) | (uint64_t(1) << selConsiderLambda) | (uint64_t(1) << selConsiderAntiLambda);
      selMap = selMap | (uint64_t(1) << selPhysPrimK0Short) | (uint64_t(1) << selPhysPrimLambda) | (uint64_t(1) << selPhysPrimAntiLambda);

      analyseCandidate(v0, v0.pt(), centrality, selMap);
    } // end v0 loop
  }

  // ______________________________________________________
  // Simulated processing (subscribes to MC information too)
  void processMonteCarlo(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraRawCents, aod::StraEvSels, aod::StraCollLabels>::iterator const& collision, v0MCCandidates const& fullV0s, dauTracks const&, aod::MotherMCParts const&, soa::Join<aod::StraMCCollisions, aod::StraMCCollMults> const& /*mccollisions*/)
  {
    histos.fill(HIST("hEventSelection"), 0. /* all collisions */);
    if (!collision.sel8()) {
      return;
    }
    histos.fill(HIST("hEventSelection"), 1 /* sel8 collisions */);

    if (std::abs(collision.posZ()) > 10.f) {
      return;
    }
    histos.fill(HIST("hEventSelection"), 2 /* vertex-Z selected */);

    if (rejectITSROFBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return;
    }
    histos.fill(HIST("hEventSelection"), 3 /* Not at ITS ROF border */);

    if (rejectTFBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return;
    }
    histos.fill(HIST("hEventSelection"), 4 /* Not at TF border */);

    float centrality = collision.centFT0C();
    if (qaCentrality) {
      auto hRawCentrality = histos.get<TH1>(HIST("hRawCentrality"));
      centrality = hRawCentrality->GetBinContent(hRawCentrality->FindBin(collision.multFT0C()));
    }

    histos.fill(HIST("hEventCentrality"), centrality);

    histos.fill(HIST("hCentralityVsNch"), centrality, collision.multNTracksPVeta1());

    // __________________________________________
    // perform main analysis
    for (auto& v0 : fullV0s) {
      if (std::abs(v0.negativeeta()) > daughterEtaCut || std::abs(v0.positiveeta()) > daughterEtaCut)
        continue; // remove acceptance that's badly reproduced by MC / superfluous in future

      // fill AP plot for all V0s
      histos.fill(HIST("GeneralQA/h2dArmenterosAll"), v0.alpha(), v0.qtarm());

      float ptmc = RecoDecay::sqrtSumOfSquares(v0.pxPosMC() + v0.pxNegMC(), v0.pyPosMC() + v0.pyNegMC());
      float ymc = 1e-3;
      if (v0.pdgCode() == 310)
        ymc = RecoDecay::y(std::array{v0.pxPosMC() + v0.pxNegMC(), v0.pyPosMC() + v0.pyNegMC(), v0.pzPosMC() + v0.pzNegMC()}, o2::constants::physics::MassKaonNeutral);
      else if (TMath::Abs(v0.pdgCode()) == 3122)
        ymc = RecoDecay::y(std::array{v0.pxPosMC() + v0.pxNegMC(), v0.pyPosMC() + v0.pyNegMC(), v0.pzPosMC() + v0.pzNegMC()}, o2::constants::physics::MassLambda);

      uint64_t selMap = computeReconstructionBitmap(v0, collision, ymc, ymc, ptmc);
      selMap = selMap | computeMCAssociation(v0);

      // feeddown matrix always with association
      if (calculateFeeddownMatrix)
        fillFeeddownMatrix(v0, ptmc, centrality, selMap);

      // consider only associated candidates if asked to do so, disregard association
      if (!doMCAssociation) {
        selMap = selMap | (uint64_t(1) << selConsiderK0Short) | (uint64_t(1) << selConsiderLambda) | (uint64_t(1) << selConsiderAntiLambda);
        selMap = selMap | (uint64_t(1) << selPhysPrimK0Short) | (uint64_t(1) << selPhysPrimLambda) | (uint64_t(1) << selPhysPrimAntiLambda);
      }

      analyseCandidate(v0, ptmc, centrality, selMap);

      if (doCollisionAssociationQA) {
        // check collision association explicitly
        bool correctCollision = false;
        int mcNch = -1;
        if (collision.has_straMCCollision()) {
          auto mcCollision = collision.straMCCollision_as<soa::Join<aod::StraMCCollisions, aod::StraMCCollMults>>();
          mcNch = mcCollision.multMCNParticlesEta05();
          correctCollision = (v0.straMCCollisionId() == mcCollision.globalIndex());
        }
        analyseCollisionAssociation(v0, ptmc, mcNch, correctCollision, selMap);
      }

    } // end v0 loop
  }

  // ______________________________________________________
  // Simulated processing (subscribes to MC information too)
  void processBinnedGenerated(
    aod::GeK0Short const& geK0Short, aod::GeLambda const& geLambda, aod::GeAntiLambda const& geAntiLambda,
    aod::GeXiMinus const& geXiMinus, aod::GeXiPlus const& geXiPlus,
    aod::GeOmegaMinus const& geOmegaMinus, aod::GeOmegaPlus const& geOmegaPlus)
  {
    auto hK0Short = histos.get<TH2>(HIST("h2dGenK0Short"));
    auto hLambda = histos.get<TH2>(HIST("h2dGenLambda"));
    auto hAntiLambda = histos.get<TH2>(HIST("h2dGenAntiLambda"));
    auto hXiMinus = histos.get<TH2>(HIST("h2dGenXiMinus"));
    auto hXiPlus = histos.get<TH2>(HIST("h2dGenXiPlus"));
    auto hOmegaMinus = histos.get<TH2>(HIST("h2dGenOmegaMinus"));
    auto hOmegaPlus = histos.get<TH2>(HIST("h2dGenOmegaPlus"));
    for (auto& gVec : geK0Short) {
      if (gVec.generatedK0Short().size() != hK0Short->GetNcells())
        LOGF(fatal, "K0Short: Number of elements in generated array and number of cells in receiving histogram differ: %i vs %i!", gVec.generatedK0Short().size(), hK0Short->GetNcells());
      for (uint64_t iv = 0; iv < hK0Short->GetNcells(); iv++) {
        hK0Short->SetBinContent(iv, hK0Short->GetBinContent(iv) + gVec.generatedK0Short()[iv]);
      }
    }
    for (auto& gVec : geLambda) {
      if (gVec.generatedLambda().size() != hLambda->GetNcells())
        LOGF(fatal, "Lambda: Number of elements in generated array and number of cells in receiving histogram differ: %i vs %i!", gVec.generatedLambda().size(), hLambda->GetNcells());
      for (uint64_t iv = 0; iv < hLambda->GetNcells(); iv++) {
        hLambda->SetBinContent(iv, hLambda->GetBinContent(iv) + gVec.generatedLambda()[iv]);
      }
    }
    for (auto& gVec : geAntiLambda) {
      if (gVec.generatedAntiLambda().size() != hAntiLambda->GetNcells())
        LOGF(fatal, "AntiLambda: Number of elements in generated array and number of cells in receiving histogram differ: %i vs %i!", gVec.generatedAntiLambda().size(), hAntiLambda->GetNcells());
      for (uint64_t iv = 0; iv < hAntiLambda->GetNcells(); iv++) {
        hAntiLambda->SetBinContent(iv, hAntiLambda->GetBinContent(iv) + gVec.generatedAntiLambda()[iv]);
      }
    }
    for (auto& gVec : geXiMinus) {
      if (gVec.generatedXiMinus().size() != hXiMinus->GetNcells())
        LOGF(fatal, "XiMinus: Number of elements in generated array and number of cells in receiving histogram differ: %i vs %i!", gVec.generatedXiMinus().size(), hXiMinus->GetNcells());
      for (uint64_t iv = 0; iv < hXiMinus->GetNcells(); iv++) {
        hXiMinus->SetBinContent(iv, hXiMinus->GetBinContent(iv) + gVec.generatedXiMinus()[iv]);
      }
    }
    for (auto& gVec : geXiPlus) {
      if (gVec.generatedXiPlus().size() != hXiPlus->GetNcells())
        LOGF(fatal, "XiPlus: Number of elements in generated array and number of cells in receiving histogram differ: %i vs %i!", gVec.generatedXiPlus().size(), hXiPlus->GetNcells());
      for (uint64_t iv = 0; iv < hXiPlus->GetNcells(); iv++) {
        hXiPlus->SetBinContent(iv, hXiPlus->GetBinContent(iv) + gVec.generatedXiPlus()[iv]);
      }
    }
    for (auto& gVec : geOmegaMinus) {
      if (gVec.generatedOmegaMinus().size() != hOmegaMinus->GetNcells())
        LOGF(fatal, "OmegaMinus: Number of elements in generated array and number of cells in receiving histogram differ: %i vs %i!", gVec.generatedOmegaMinus().size(), hOmegaMinus->GetNcells());
      for (uint64_t iv = 0; iv < hOmegaMinus->GetNcells(); iv++) {
        hOmegaMinus->SetBinContent(iv, hOmegaMinus->GetBinContent(iv) + gVec.generatedOmegaMinus()[iv]);
      }
    }
    for (auto& gVec : geOmegaPlus) {
      if (gVec.generatedOmegaPlus().size() != hOmegaPlus->GetNcells())
        LOGF(fatal, "OmegaPlus: Number of elements in generated array and number of cells in receiving histogram differ: %i vs %i!", gVec.generatedOmegaPlus().size(), hOmegaPlus->GetNcells());
      for (uint64_t iv = 0; iv < hOmegaPlus->GetNcells(); iv++) {
        hOmegaPlus->SetBinContent(iv, hOmegaPlus->GetBinContent(iv) + gVec.generatedOmegaPlus()[iv]);
      }
    }
  }

  PROCESS_SWITCH(derivedlambdakzeroanalysis, processRealData, "process as if real data", true);
  PROCESS_SWITCH(derivedlambdakzeroanalysis, processMonteCarlo, "process as if MC", false);
  PROCESS_SWITCH(derivedlambdakzeroanalysis, processBinnedGenerated, "process MC generated", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<derivedlambdakzeroanalysis>(cfgc)};
}
