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
/// \file hStrangeCorrelation.cxx
/// \brief This task serves to do hadron-(strange hadron) correlation studies.
///  The yield will be calculated using the two-particle correlation method.
///  Trigger particle : Hadrons
///  Associated Particles : V0s or Cascades
///  this task requires the hStrangeCorrelationFilter to have been run before.
///
/// \author Kai Cui (kaicui@mails.ccnu.edu.cn)
/// \author Lucia Anna Tarasovicova (lucia.anna.husova@cern.ch)
/// \author David Dobrigkeit Chinellato (david.dobrigkeit.chinellato@cern.ch)
/// \author Zhongbao Yin (Zhong-Bao.Yin@cern.ch)

#include "PWGLF/DataModel/LFHStrangeCorrelationTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/mcCentrality.h"
#include "PWGLF/Utils/inelGt.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <Framework/ASoA.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/StaticFor.h>
#include <Framework/runDataProcessing.h>

#include <TAxis.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THn.h>
#include <THnSparse.h>
#include <TList.h>
#include <TPDGCode.h>
#include <TString.h>

#include <fmt/format.h>

#include <Rtypes.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <variant>
#include <vector>

using namespace o2;
using namespace o2::constants::math;
using namespace o2::framework;
using namespace o2::framework::expressions;

using TracksComplete = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA>;
using TracksCompleteMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels>;
using V0DatasWithoutTrackX = soa::Join<aod::V0Indices, aod::V0Cores>;
using V0DatasWithoutTrackXMC = soa::Join<aod::V0Indices, aod::V0Cores, aod::McV0Labels>;
using V0DatasWithoutTrackXWithMCCores = soa::Join<aod::V0Indices, aod::V0Cores, aod::V0MCCores>;

struct HStrangeCorrelation {
  // for efficiency corrections if requested
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  Service<o2::framework::O2DatabasePDG> pdgDB;
  o2::pwglf::ParticleCounter<o2::framework::O2DatabasePDG> mCounter;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // event filtering
  Configurable<std::string> zorroMask{"zorroMask", "", "zorro trigger class to select on (empty: none)"};

  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  struct : ConfigurableGroup {
    std::string prefix = "masterConfigurations";
    Configurable<bool> doPPAnalysis{"doPPAnalysis", true, "if in pp, set to true"};
    Configurable<int> collisionHasTriggOrAssoc{"collisionHasTriggOrAssoc", 0, "require the collisions containing (0:no requirement 1:trig 2:assoc 3:trig or assoc 4:trig and assoc"};
    Configurable<bool> doFullCorrelationStudy{"doFullCorrelationStudy", true, "if true, do full correlation study by creating all THnSparse histograms for the correlation function"};
    Configurable<bool> doCorrelationHadron{"doCorrelationHadron", false, "do Hadron correlation"};
    Configurable<bool> doCorrelationK0Short{"doCorrelationK0Short", true, "do K0Short correlation"};
    Configurable<bool> doCorrelationLambda{"doCorrelationLambda", false, "do Lambda correlation"};
    Configurable<bool> doCorrelationAntiLambda{"doCorrelationAntiLambda", false, "do AntiLambda correlation"};
    Configurable<bool> doCorrelationXiMinus{"doCorrelationXiMinus", false, "do XiMinus correlation"};
    Configurable<bool> doCorrelationXiPlus{"doCorrelationXiPlus", false, "do XiMinus correlation"};
    Configurable<bool> doCorrelationOmegaMinus{"doCorrelationOmegaMinus", false, "do OmegaMinus correlation"};
    Configurable<bool> doCorrelationOmegaPlus{"doCorrelationOmegaPlus", false, "do OmegaPlus correlation"};
    Configurable<bool> doCorrelationPion{"doCorrelationPion", false, "do Pion correlation"};
    Configurable<bool> doGenEventSelection{"doGenEventSelection", true, "use event selections when performing closure test for the gen events"};
    Configurable<bool> selectINELgtZERO{"selectINELgtZERO", true, "select INEL>0 events"};
    Configurable<bool> selectINELgtONE{"selectINELgtONE", false, "select INEL>1 events (at least 2 charged particles in |eta| < 1)"};
    Configurable<float> zVertexCut{"zVertexCut", 10, "Cut on PV position"};
    Configurable<bool> requireAllGoodITSLayers{"requireAllGoodITSLayers", false, " require that in the event all ITS are good"};
    Configurable<bool> rejectSameBunchPileup{"rejectSameBunchPileup", false, "reject collisions associated with the same found-by-T0 bunch crossing"};
    Configurable<bool> requireGoodTriggerTVX{"requireGoodTriggerTVX", false, " require acceptable FT0C-FT0A time difference"};
    Configurable<bool> requireGoodZvtxFT0vsPV{"requireGoodZvtxFT0vsPV", false, " require small difference between z-vertex from PV and from FT0"};
    Configurable<bool> skipUnderOverflowInTHn{"skipUnderOverflowInTHn", false, "skip under/overflow in THns"};
    Configurable<int> mixingParameter{"mixingParameter", 10, "how many events are mixed"};
    Configurable<bool> doMCassociation{"doMCassociation", false, "fill everything only for MC associated"};
    Configurable<bool> doTriggPhysicalPrimary{"doTriggPhysicalPrimary", false, "require physical primary for trigger particles"};
    Configurable<bool> doClosureTestTriggerWithRecoTrackMatch{"doClosureTestTriggerWithRecoTrackMatch", false, "add closure-test correlations requiring the truth trigger to have at least one reconstructed track with a matching MC label"};
    Configurable<bool> applyNewMCSelection{"applyNewMCSelection", false, "apply new MC Generated selection"};
    Configurable<bool> doSeparateFT0Prediction{"doSeparateFT0Prediction", false, "separate FT0M to FT0A and FT0C in prediction process"};
    Configurable<bool> useCentralityinPrediction{"useCentralityinPrediction", false, "if true, use centrality instead of multiplisity"};
    Configurable<bool> doMirroringInDelataEta{"doMirroringInDelataEta", false, "if true, fill only positive delta eta and mirror the negative side in post processing, Adjust the delta axis!"};
    Configurable<bool> doMassSpectrumCheck{"doMassSpectrumCheck", false, "if true, add and fill invariant-mass spectrum"};
    Configurable<bool> doCorrelationsHadronV0daughter{"doCorrelationsHadronV0daughter", false, "if true, do correlations of hadrons with V0 daughters"};
    Configurable<bool> doLocalDensityStudy{"doLocalDensityStudy", false, "if true, create and fill the pt vs eta vs local density spectra of triggers and V0s"};
    Configurable<float> localDensityConeRadius{"localDensityConeRadius", 0.4, "radius of the cone in which the local density is counted"};
  } masterConfigurations;

  // master analysis switches
  Configurable<bool> doAssocPhysicalPrimary{"doAssocPhysicalPrimary", false, "require physical primary for associated particles"};
  Configurable<bool> doAssocPhysicalPrimaryInGen{"doAssocPhysicalPrimaryInGen", false, "require physical primary for associated particles in Generated Partilces"};
  Configurable<bool> doLambdaPrimary{"doLambdaPrimary", false, "do primary selection for lambda"};
  Configurable<bool> doAutocorrelationRejection{"doAutocorrelationRejection", true, "reject pairs where trigger Id is the same as daughter particle Id"};
  Configurable<bool> doMixingQAandEventQA{"doMixingQAandEventQA", true, "if true, add EvnetQA and MixingQA hist to histos"};
  Configurable<bool> doITSClustersQA{"doITSClustersQA", true, "if true, add ITSCluster hist to histos"};
  Configurable<bool> doDeltaPhiStarCheck{"doDeltaPhiStarCheck", false, "if true, create and fill delta phi star histograms"};

  Configurable<int> triggerBinToSelect{"triggerBinToSelect", 0, "trigger bin to select on if processSelectEventWithTrigger enabled"};
  Configurable<int> triggerParticleCharge{"triggerParticleCharge", 0, "For checks, if 0 all charged tracks, if -1 only neg., if 1 only positive"};
  Configurable<float> etaSel{"etaSel", 0.8, "Selection in eta for trigger and associated particles"};
  Configurable<float> ySel{"ySel", 0.5, "Selection in rapidity for consistency checks"};

  Configurable<bool> useTheLeadingParticleAsTrigger{"useTheLeadingParticleAsTrigger", false, "if true, use the leading particle in the event as trigger particle"};
  // used for event selections in Pb-Pb
  Configurable<int> cfgCutOccupancyHigh{"cfgCutOccupancyHigh", 3000, "High cut on TPC occupancy"};
  Configurable<int> cfgCutOccupancyLow{"cfgCutOccupancyLow", 0, "Low cut on TPC occupancy"};

  // Axes - configurable for smaller sizes
  struct : ConfigurableGroup {
    std::string prefix = "axesConfigurations";
    ConfigurableAxis axisMult{"axisMult", {VARIABLE_WIDTH, 0.0f, 0.01f, 1.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 70.0f, 100.0f}, "Mixing bins - multiplicity"};
    ConfigurableAxis axisVtxZ{"axisVtxZ", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
    ConfigurableAxis axisPhi{"axisPhi", {72, 0, TwoPI}, "#phi"};
    ConfigurableAxis axisEta{"axisEta", {80, -0.8, +0.8}, "#eta"};
    ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {72, -PIHalf, PIHalf * 3}, "delta #varphi axis for histograms"};
    ConfigurableAxis axisDeltaEta{"axisDeltaEta", {50, -1.6, 1.6}, "delta eta axis for histograms"};
    ConfigurableAxis axisPtAssoc{"axisPtAssoc", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt associated axis for histograms"};
    ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.0, 1.0, 2.0, 3.0, 100}, "pt associated axis for histograms"};
    ConfigurableAxis axisPtQA{"axisPtQA", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for QA histograms"};
    ConfigurableAxis axisMassNSigma{"axisMassNSigma", {40, -2, 2}, "Axis for mass Nsigma"};
    ConfigurableAxis axisLocalDensity{"axisLocalDensity", {21, -0.5, 20.5}, "local density (number of tracks in the cone)"};
    ConfigurableAxis axisK0ShortMass{"axisK0ShortMass", {200, 0.400f, 0.600f}, "Inv. Mass (GeV/c^{2})"};
    ConfigurableAxis axisLambdaMass{"axisLambdaMass", {200, 1.01f, 1.21f}, "Inv. Mass (GeV/c^{2})"};
    ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 20, 40, 60, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300}, "Binning of the Multiplicity axis in model prediction process"};
    ConfigurableAxis axisMidrapidityMultiplicity{"axisMidrapidityMultiplicity", {VARIABLE_WIDTH, 0, 20, 40, 60, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300}, "Binning of the Midrapidity Multiplicity axis in model prediction process"};

  } axesConfigurations;

  // for topo var QA
  struct : ConfigurableGroup {
    std::string prefix = "massWindowConfigurations";
    Configurable<float> maxPeakNSigma{"maxPeakNSigma", 5, "Peak region edge definition (in sigma)"};
    Configurable<float> minBgNSigma{"minBgNSigma", 5, "Bg region edge closest to peak (in sigma)"};
    Configurable<float> maxBgNSigma{"maxBgNSigma", 10, "Bg region edge furthest to peak (in sigma)"};
    Configurable<float> nSigmaNearXiMassCenter{"nSigmaNearXiMassCenter", 1, "for Oemga analysis only, to check if candidate mass is around Xi"};
  } massWindowConfigurations; // allows for gap between peak and bg in case someone wants to

  // Implementation of on-the-spot efficiency correction
  struct : ConfigurableGroup {
    std::string prefix = "efficiencyFlags";
    Configurable<bool> applyEfficiencyCorrection{"applyEfficiencyCorrection", false, "apply efficiency correction"};
    Configurable<bool> applyEfficiencyForTrigger{"applyEfficiencyForTrigger", false, "apply efficiency correction for the trigger particle"};
    Configurable<bool> applyEfficiencyPropagation{"applyEfficiencyPropagation", false, "propagate also the efficiency uncertainty"};
    Configurable<bool> applyPurityHadron{"applyPurityHadron", false, "apply the purity correction for associated hadrons"};
    Configurable<bool> applyPurityTrigger{"applyPurityTrigger", false, "apply the purity correction for trigger particle"};
    Configurable<bool> applyEffAsFunctionOfMult{"applyEffAsFunctionOfMult", false, "apply efficiency as a function of multiplicity as well"};
    Configurable<bool> applyEffAsFunctionOfMultAndPhi{"applyEffAsFunctionOfMultAndPhi", false, "apply efficiency as a function of multiplicity and phi"};
  } efficiencyFlags;
  Configurable<std::string> ccdburl{"ccdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository to use"};
  Configurable<std::string> efficiencyCCDBPath{"efficiencyCCDBPath", "GLO/Config/GeometryAligned", "Path of the efficiency corrections"};

  // Configurables for doing subwagon systematics
  struct : ConfigurableGroup {
    std::string prefix = "trackSelection";
    // --- Track quality variations (single track, both trigger and assoc daughters)
    Configurable<int> minTPCNCrossedRowsTrigger{"minTPCNCrossedRowsTrigger", 70, "Minimum TPC crossed rows (trigger)"};
    Configurable<int> minTPCNCrossedRowsAssociated{"minTPCNCrossedRowsAssociated", 70, "Minimum TPC crossed rows (associated)"};
    Configurable<bool> triggerRequireITS{"triggerRequireITS", true, "require ITS signal in trigger tracks"};
    Configurable<bool> assocRequireITS{"assocRequireITS", true, "require ITS signal in associated primary tracks"};
    Configurable<int> triggerMaxTPCSharedClusters{"triggerMaxTPCSharedClusters", 200, "maximum number of shared TPC clusters (inclusive)"};
    Configurable<int> assocMaxTPCSharedClusters{"assocMaxTPCSharedClusters", 200, "maximum number of shared TPC clusters (inclusive) for assoc primary tracks"};
    Configurable<bool> triggerRequireL0{"triggerRequireL0", false, "require ITS L0 cluster for trigger"};
    Configurable<bool> assocRequireL0{"assocRequireL0", true, "require ITS L0 cluster for assoc primary track"};
    Configurable<float> minTPCChi2PerClusterAssociated{"minTPCChi2PerClusterAssociated", 4.0f, "Minimum TPC chi2 per cluster for associated primary tracks"};
    Configurable<bool> checksRequireTPCChi2{"checksRequireTPCChi2", false, "require TPC chi2 per cluster for trigger and associated primary tracks"};
    Configurable<bool> requireClusterInITS{"requireClusterInITS", false, "require cluster in ITS for V0 and cascade daughter tracks"};
    Configurable<int> minITSClustersForDaughterTracks{"minITSClustersForDaughterTracks", 1, "Minimum number of ITS clusters for V0 daughter tracks"};
    Configurable<bool> requireDCAzCut{"requireDCAzCut", false, "require DCAz cut for trigger and associated primary tracks"};
    Configurable<bool> doITSChi2Selection{"doITSChi2Selection", false, "require ITS chi2 cut for trigger and associated primary tracks"};
    Configurable<bool> checkForITSTPCMissmatchMC{"checkForITSTPCMissmatchMC", false, "if true, reject tracks with MC mask bit 13 (ITS-TPC mismatch) for trigger and associated primary tracks"};
    // --- Trigger: DCA variation from basic formula: |DCAxy| <  0.004f + (0.013f / pt)
    Configurable<float> dcaXYconstant{"dcaXYconstant", 0.004, "[0] in |DCAxy| < [0]+[1]/pT"};
    Configurable<float> dcaXYpTdep{"dcaXYpTdep", 0.013, "[1] in |DCAxy| < [0]+[1]/pT"};
    // --- Assoc track: DCA variation from basic formula: |DCAxy| <  0.004f + (0.013f / pt)
    Configurable<float> dcaXYconstantAssoc{"dcaXYconstantAssoc", 0.004, "[0] in |DCAxy| < [0]+[1]/pT"};
    Configurable<float> dcaXYpTdepAssoc{"dcaXYpTdepAssoc", 0.013, "[1] in |DCAxy| < [0]+[1]/pT"};

    Configurable<float> dcaZconstant{"dcaZconstant", 0.004, "[0] in |DCAz| < [0]+[1]/pT"};
    Configurable<float> dcaZpTdep{"dcaZpTdep", 0.013, "[1] in |DCAz| < [0]+[1]/pT"};
    Configurable<float> dcaZconstantAssoc{"dcaZconstantAssoc", 0.004, "[0] in |DCAz| < [0]+[1]/pT"};
    Configurable<float> dcaZpTdepAssoc{"dcaZpTdepAssoc", 0.013, "[1] in |DCAz| < [0]+[1]/pT"};
    Configurable<int> dEdxCompatibility{"dEdxCompatibility", 1, "0: loose, 1: normal, 2: tight. Defined in HStrangeCorrelationFilter"};
    Configurable<float> chi2ITSfit{"chi2ITSfit", 36, "chi2 ITS fit for trigger and associated primary tracks"};

  } trackSelection;
  struct : ConfigurableGroup {
    std::string prefix = "v0Selection";
    // --- V0 selections (both K0Short and Lambda, but some can be Lambda specific, see below)
    // --- Associated: topological variable variation (OK to vary all-at-once, at least for first study)
    Configurable<double> v0cospa{"v0cospa", 0.97, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0)
    Configurable<float> dcaV0dau{"dcaV0dau", 1.0, "DCA V0 Daughters"};
    Configurable<float> dcanegtopv{"dcanegtopv", 0.06, "DCA Neg To PV"};
    Configurable<float> dcapostopv{"dcapostopv", 0.06, "DCA Pos To PV"};
    Configurable<float> v0RadiusMin{"v0RadiusMin", 0.5, "v0radius"};
    Configurable<float> v0RadiusMax{"v0RadiusMax", 200, "v0radius"};
    Configurable<float> dcaBaryonToPV{"dcaBaryonToPV", 0.2, "DCA of baryon daughter track To PV"};
    Configurable<float> dcaMesonToPV{"dcaMesonToPV", 0.05, "DCA of meson daughter track To PV"};
    Configurable<float> dcaDaugToPVForK0s{"dcaDaugToPVForK0s", 0, "DCA of K0s daughter tracks To PV"};
    Configurable<float> lifetimecutK0S{"lifetimecutK0S", 20, "lifetimecutK0S"};
    Configurable<float> lifetimecutLambda{"lifetimecutLambda", 30, "lifetimecutLambda"};
    // original equation: lArmPt*2>TMath::Abs(lArmAlpha) only for K0S
    Configurable<float> armPodCut{"armPodCut", 5.0f, "pT * (cut) > |alpha|, AP cut. Negative: no cut"};
  } v0Selection;

  struct : ConfigurableGroup {
    std::string prefix = "cascadeSelections";
    Configurable<double> cascCospa{"cascCospa", 0.95, "cascCospa"};
    Configurable<float> dcaCascDaughters{"dcaCascDaughters", 1.0, "DCA between the V0 candidate and the bachelor track"};
    // Configurable<float> cascDcabachtopv{"cascDcabachtopv", 0.1, "cascDcabachtopv"};
    Configurable<float> cascRadius{"cascRadius", 0.5, "cascRadius"};
    Configurable<float> cascV0masswindow{"cascV0masswindow", 0.01, "cascV0masswindow"};
    Configurable<float> bachBaryonCosPA{"bachBaryonCosPA", 0.9999, "Bachelor baryon CosPA"};
    Configurable<float> bachBaryonDCAxyToPV{"bachBaryonDCAxyToPV", 0.08, "DCA bachelor baryon to PV"};
    Configurable<float> cascdcaV0dau{"cascdcaV0dau", 0.5, "DCA V0 Daughters"};
    Configurable<float> proplifetime{"proplifetime", 3, "ctau/<ctau>"};
    Configurable<float> rejcomp{"rejcomp", 0.008, "Competing Cascade rejection"};
    Configurable<float> rapCut{"rapCut", 0.8, "Rapidity acceptance"};
    Configurable<float> dcaCacsDauPar0{"dcaCacsDauPar0", 0.8, " par for pt dep DCA cascade daughter cut, p_T < 1 GeV/c"};
    Configurable<float> dcaCacsDauPar1{"dcaCacsDauPar1", 0.5, " par for pt dep DCA cascade daughter cut, 1< p_T < 4 GeV/c"};
    Configurable<float> dcaCacsDauPar2{"dcaCacsDauPar2", 0.2, " par for pt dep DCA cascade daughter cut, p_T > 4 GeV/c"};
    Configurable<float> cascdcaV0ToPV{"cascdcaV0ToPV", 0.06, "DCA V0 To PV"};
    Configurable<double> cascv0cospa{"cascv0cospa", 0.98, "V0 CosPA"};
    Configurable<float> cascv0RadiusMin{"cascv0RadiusMin", 2.5, "v0radius"};
    Configurable<float> dcaBaryonToPV{"dcaBaryonToPV", 0.05, "DCA of baryon doughter track To PV"};
    Configurable<float> dcaMesonToPV{"dcaMesonToPV", 0.1, "DCA of meson doughter track To PV"};
    Configurable<float> dcaBachToPV{"dcaBachToPV", 0.07, "DCA Bach To PV"};
    // pt Range for pt dep cuts
    Configurable<float> highPtForCascDaugPtDep{"highPtForCascDaugPtDep", 4.0, "high pt range for pt dep cuts"};
    Configurable<float> lowPtForCascDaugPtDep{"lowPtForCascDaugPtDep", 1.0, "low pt range for pt dep cuts"};

  } cascadeSelections;

  struct : ConfigurableGroup {
    std::string prefix = "checks";
    // cascade selections
    // more cascade selections in PbPb
    // Configurable<float> bachBaryonCosPA{"bachBaryonCosPA", 0.9999, "Bachelor baryon CosPA"};
    // Configurable<float> bachBaryonDCAxyToPV{"bachBaryonDCAxyToPV", 0.08, "DCA bachelor baryon to PV"};
    //
    // on the fly correction instead of mixingParameter
    Configurable<bool> doOnTheFlyFlattening{"doOnTheFlyFlattening", false, "enable an on-the-fly correction instead of using mixing"};

    // delta eta ranges for toward and transverse checks
    Configurable<float> transwerseDeltaEtaRangeMin{"transwerseDeltaEtaRangeMin", 1.0, "minimum delta eta for transverse region"};
    Configurable<float> transwerseDeltaEtaRangeMax{"transwerseDeltaEtaRangeMax", 2.0, "maximum delta eta for transverse region"};
    Configurable<float> towardDeltaEtaRange{"towardDeltaEtaRange", 0.8, "delta eta range for toward region"};
    // (N.B.: sources that can be investigated in post are not listed!)
  } checks;

  struct : ConfigurableGroup {
    std::string prefix = "pairLossK0Configurations";
    // processPairLossK0MC is split into three independent parts. Each has its own
    // switch, they write into disjoint folders, and any combination of them may
    // run in the same job -- including all three at once.
    //   doStageDiagnostics  PairLossK0/{Stage,State,Geometry,Matching,Response,
    //                       TrackQA,V0QA}: the truth-pair reconstruction ladder
    //                       and its close-pair diagnostics
    //   doRecComparison     PairLossK0/Comparison: the cumulative
    //                       Rec/Truth/Gen/Final variant ladder. Runs the exact
    //                       reconstructed correlation path internally, so
    //                       processSameEventHV0s must be off when it is on.
    //   doGenLevelStudy     PairLossK0/GenStudy: generator-level only, see the
    //                       comment on runGenLevelStudy in processPairLossK0MC
    Configurable<bool> doStageDiagnostics{"doStageDiagnostics", true, "part 1: fill the PairLossK0 truth-pair reconstruction ladder and its diagnostics"};
    Configurable<bool> doRecComparison{"doRecComparison", false, "part 2: fill the PairLossK0/Comparison cumulative Rec/Truth/Gen/Final ladder (runs the exact reconstructed path internally)"};
    Configurable<bool> doGenLevelStudy{"doGenLevelStudy", false, "part 3: fill the PairLossK0/GenStudy generated-vs-reconstructed split, using generator-level event selection only"};
    // Generated charged multiplicity of the MC collision, counted in |eta| < 0.8
    // by mCounter. Plain ConfigurableAxis: unlike the correlation axes it is NOT
    // trimmed by skipUnderOverflowInTHn, so what you configure is what you get.
    ConfigurableAxis axisGenStudyNch{"axisGenStudyNch", {VARIABLE_WIDTH, 0.0f, 2.0f, 5.0f, 10.0f, 15.0f, 20.0f, 25.0f, 30.0f, 40.0f, 60.0f, 100.0f}, "generated charged multiplicity in |#eta| < 0.8"};
    Configurable<bool> doClosureTestStages{"doClosureTestStages", true, "create and fill the whole ClosureTest/PairLossK0 folder: the truth and any-reconstructed-object stages of the truth h-K0 pair, mirroring the first processPairLossK0MC stages"};
    Configurable<bool> applyRecoEventSelection{"applyRecoEventSelection", true, "apply the standard reconstructed-event selection in the K0 pair-loss diagnostic"};
    Configurable<bool> fillFinalPairOnce{"fillFinalPairOnce", true, "part 1: fill each truth h-K0 pair at most once at the final stage. Set false to loop over all reconstructed trigger x K0 combinations, as the data path does"};
    Configurable<bool> finalPairUseBestCollisionOnly{"finalPairUseBestCollisionOnly", true, "part 1: build the final trigger and K0 objects only in the best collision. Set false to use every associated reconstructed collision that passes the event selection, pairing within one collision, as the data path does"};
    Configurable<float> daughterPtMin{"daughterPtMin", 0.05f, "minimum generated daughter pT for the findable K0 category"};
    Configurable<float> daughterEtaMax{"daughterEtaMax", 0.9f, "maximum absolute generated daughter eta for the findable K0 category"};
    Configurable<float> minRadius{"minRadius", 0.8f, "minimum radius in metres used for the minimum delta-phi-star search"};
    Configurable<float> maxRadius{"maxRadius", 2.5f, "maximum radius in metres used for the minimum delta-phi-star search"};
    Configurable<float> radiusStep{"radiusStep", 0.05f, "radius step in metres used for the minimum delta-phi-star search"};
    ConfigurableAxis axisMinDeltaPhiStar{"axisMinDeltaPhiStar", {120, 0.0f, 0.3f}, "minimum absolute delta phi star"};
    ConfigurableAxis axisSignedDeltaPhiStar{"axisSignedDeltaPhiStar", {120, -0.3f, 0.3f}, "signed delta phi star at the minimum"};
    ConfigurableAxis axisDaughterDeltaEta{"axisDaughterDeltaEta", {100, -0.2f, 0.2f}, "trigger-daughter delta eta"};
    ConfigurableAxis axisDaughterDeltaPhi{"axisDaughterDeltaPhi", {120, -0.3f, 0.3f}, "trigger-daughter delta phi"};
    ConfigurableAxis axisDecayRadius{"axisDecayRadius", {100, 0.0f, 200.0f}, "generated K0 decay radius (cm)"};
    // The following mirror the isValidTrigger() cuts applied by hStrangeCorrelationFilter.cxx
    // when building the TriggerTracks table. They must be kept in sync by hand with whatever
    // values are configured for that (separate) filter task, since this diagnostic cannot read
    // them directly: it only sees the already-filtered TriggerTracks table, not why a given
    // best-collision track failed to enter it.
    Configurable<float> triggerTracksEtaMin{"triggerTracksEtaMin", -0.8f, "must match triggerEtaMin in hStrangeCorrelationFilter.cxx"};
    Configurable<float> triggerTracksEtaMax{"triggerTracksEtaMax", 0.8f, "must match triggerEtaMax in hStrangeCorrelationFilter.cxx"};
    Configurable<float> triggerTracksPtMin{"triggerTracksPtMin", 3.0f, "must match triggerPtCutMin in hStrangeCorrelationFilter.cxx"};
    Configurable<float> triggerTracksPtMax{"triggerTracksPtMax", 20.0f, "must match triggerPtCutMax in hStrangeCorrelationFilter.cxx"};
    Configurable<int> triggerTracksMinCrossedRows{"triggerTracksMinCrossedRows", 70, "must match minTPCNCrossedRows in hStrangeCorrelationFilter.cxx"};
    Configurable<bool> triggerTracksRequireITS{"triggerTracksRequireITS", true, "must match triggerRequireITS in hStrangeCorrelationFilter.cxx"};
    Configurable<int> triggerTracksMaxSharedClusters{"triggerTracksMaxSharedClusters", 200, "must match triggerMaxTPCSharedClusters in hStrangeCorrelationFilter.cxx"};
    Configurable<bool> triggerTracksRequireLayer0{"triggerTracksRequireLayer0", false, "must match triggerRequireL0 in hStrangeCorrelationFilter.cxx"};
    // Separate from trackSelection.dcaXY*, which is this task's systematic-variation knob: the
    // replica must follow the value the filter was run with, not the varied analysis-level one.
    Configurable<float> triggerTracksDcaXYconstant{"triggerTracksDcaXYconstant", 0.004f, "must match dcaXYconstant in hStrangeCorrelationFilter.cxx"};
    Configurable<float> triggerTracksDcaXYpTdep{"triggerTracksDcaXYpTdep", 0.013f, "must match dcaXYpTdep in hStrangeCorrelationFilter.cxx"};
  } pairLossK0Configurations;

  struct ValidCollision {
    struct ValidParticle {
      float eta;
      float phi;
      float pt;
      int region;
      float efficiency;
      float efficiencyError;
      int type;
    };
    float pvz = 0.0f;
    float mult = 0.0f;
    std::vector<ValidParticle> trigParticles;
    std::vector<ValidParticle> assocParticles;
    void addValidParticle(float eta, float phi, float pt, int region, float efficiency, float efficiencyError, int type)
    {
      ValidParticle particle{.eta = eta, .phi = phi, .pt = pt, .region = region, .efficiency = efficiency, .efficiencyError = efficiencyError, .type = type};

      if (type == -1) {
        trigParticles.push_back(particle);
      } else {
        assocParticles.push_back(particle);
      }
    }
  };

  using ValidCollisions = std::vector<std::vector<ValidCollision>>;
  ValidCollisions validCollisions;

  // objects to use for efficiency corrections
  TH2F* hEfficiencyTrigger = nullptr;
  TH3F* hEfficiencyTriggerMult = nullptr;
  THnT<float>* hEfficiencyTriggerMultVsPhi = nullptr;
  TH2F* hEfficiencyPion = nullptr;
  TH2F* hEfficiencyK0Short = nullptr;
  THnT<float>* hEfficiencyK0ShortMultVsPhi = nullptr;
  TH2F* hEfficiencyLambda = nullptr;
  THnT<float>* hEfficiencyLambdaMultVsPhi = nullptr;
  TH2F* hEfficiencyAntiLambda = nullptr;
  THnT<float>* hEfficiencyAntiLambdaMultVsPhi = nullptr;
  TH2F* hEfficiencyXiMinus = nullptr;
  THnT<float>* hEfficiencyXiMinusMultVsPhi = nullptr;
  TH2F* hEfficiencyXiPlus = nullptr;
  THnT<float>* hEfficiencyXiPlusMultVsPhi = nullptr;
  TH2F* hEfficiencyOmegaMinus = nullptr;
  THnT<float>* hEfficiencyOmegaMinusMultVsPhi = nullptr;
  TH2F* hEfficiencyOmegaPlus = nullptr;
  THnT<float>* hEfficiencyOmegaPlusMultVsPhi = nullptr;
  TH2F* hEfficiencyHadron = nullptr;
  TH3F* hEfficiencyHadronMult = nullptr;
  TH1F* hPurityHadron = nullptr;
  TH2F* hPurityHadronMult = nullptr;
  // objects to propagate the efficiency uncertainty
  TH2F* hEfficiencyUncertaintyTrigger = nullptr;
  TH3F* hEfficiencyUncertaintyTriggerMult = nullptr;
  TH2F* hEfficiencyUncertaintyPion = nullptr;
  TH2F* hEfficiencyUncertaintyK0Short = nullptr;
  TH2F* hEfficiencyUncertaintyLambda = nullptr;
  TH2F* hEfficiencyUncertaintyAntiLambda = nullptr;
  TH2F* hEfficiencyUncertaintyXiMinus = nullptr;
  TH2F* hEfficiencyUncertaintyXiPlus = nullptr;
  TH2F* hEfficiencyUncertaintyOmegaMinus = nullptr;
  TH2F* hEfficiencyUncertaintyOmegaPlus = nullptr;
  TH2F* hEfficiencyUncertaintyHadron = nullptr;
  TH3F* hEfficiencyUncertaintyHadronMult = nullptr;
  TH1F* hPurityUncertaintyHadron = nullptr;
  TH2F* hPurityUncertaintyHadronMult = nullptr;

  using BinningTypePP = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;
  using BinningTypePbPb = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;

  // collision slicing for mixed events
  Preslice<aod::TriggerTracks> collisionSliceTracks = aod::triggerTracks::collisionId;
  Preslice<aod::AssocV0s> collisionSliceV0s = aod::assocV0s::collisionId;
  Preslice<aod::AssocCascades> collisionSliceCascades = aod::assocCascades::collisionId;
  // Preslice<aod::AssocPions> collisionSlicePions = aod::assocHadrons::collisionId;
  Preslice<aod::AssocHadrons> collisionSliceHadrons = aod::assocHadrons::collisionId;
  Preslice<aod::McParticles> perCollision = aod::mcparticle::mcCollisionId;
  Preslice<TracksCompleteMC> pairLossTracksPerCollision = aod::track::collisionId;
  Preslice<V0DatasWithoutTrackXWithMCCores> pairLossV0sPerCollision = aod::v0data::collisionId;

  static constexpr std::array<std::string_view, 3> V0names = {"K0Short", "Lambda", "AntiLambda"};
  static constexpr std::array<std::string_view, 4> Cascadenames = {"XiMinus", "XiPlus", "OmegaMinus", "OmegaPlus"};
  static constexpr std::array<std::string_view, 9> Particlenames = {"K0Short", "Lambda", "AntiLambda", "XiMinus", "XiPlus", "OmegaMinus", "OmegaPlus", "Pion", "Hadron"};
  static constexpr std::array<int, 8> PdgCodes = {310, 3122, -3122, 3312, -3312, 3334, -3334, 211};

  static constexpr int IndexPion = 7;
  static constexpr int IndexXiMinus = 0;
  static constexpr int IndexXiPlus = 1;
  static constexpr int IndexOmegaMinus = 2;
  static constexpr int IndexOmegaPlus = 3;
  static constexpr int IndexK0 = 0;
  static constexpr int IndexLambda = 1;
  static constexpr int IndexAntiLambda = 2;

  template <int Index, typename TV0>
  static float getV0InvariantMass(TV0 const& v0)
  {
    static_assert(Index >= IndexK0 && Index <= IndexAntiLambda);
    if constexpr (Index == IndexK0) {
      return v0.mK0Short();
    } else if constexpr (Index == IndexLambda) {
      return v0.mLambda();
    } else {
      return v0.mAntiLambda();
    }
  }

  uint16_t doCorrelation = 0;
  int mRunNumber = -1;
  int mRunNumberZorro = -1;

  std::vector<std::vector<float>> axisRanges;

  static constexpr float ctauxi = 4.91;     // from PDG
  static constexpr float ctauomega = 2.461; // from PDG

  static constexpr float MinRadiusTPC = 0.8;
  static constexpr float MaxRadiusTPC = 2.5;

  static constexpr float Neutral = 0.0;

  static constexpr int AssocParticleTypes = 9;         // K0S, Lambda, AntiLambda, Xi-, Xi+, Omega-, Omega+, Pion, Hadron
  static constexpr int AssocParticleTypesNoHadron = 8; // K0S, Lambda, AntiLambda, Xi-, Xi+, Omega-, Omega+, Pion
  static constexpr int AssocV0Types = 3;               // K0S, Lambda, AntiLambda,
  static constexpr int AssocCascadeTypes = 4;          // Xi-, Xi+, Omega-, Omega+

  // Cumulative reconstruction ladder for one truth h-K0 pair: every stage is a
  // strictly narrower requirement than the one before it, so the ratio of two
  // neighbouring stages is the efficiency of exactly the step between them.
  //
  // Two "pure reconstruction" levels anchor the chain. Both mean "the object is
  // present in the reconstruction with no selection applied whatsoever":
  //   PairLossTriggerPureReco  a track carrying the trigger's MC label exists
  //   PairLossV0PureReco       both K0 daughters have a reconstructed track
  // The V0 one sits deliberately at daughter-track level rather than at V0Datas
  // level: a row in V0Datas has already survived the V0 builder's own cuts
  // (cos(PA), daughter DCA, radius, crossed rows), so the step
  // PairLossV0PureReco -> PairLossV0Candidate isolates exactly what the builder
  // throws away, which no other stage can show.
  //
  // Everything is evaluated in the best collision. The "any reconstructed
  // collision" variants live in the ClosureTest/PairLossK0/AnyTrack* folders
  // instead, so keeping them here as well would only duplicate them.
  enum PairLossK0Stage : int {
    PairLossGenPair = 0,
    PairLossFindablePair,
    PairLossTriggerPureReco,
    PairLossTriggerInTable,
    PairLossTriggerFinal,
    PairLossV0PureReco,
    PairLossV0Candidate,
    PairLossV0InTable,
    PairLossV0Final,
    PairLossFinalPair,
    PairLossK0NStages
  };

  static constexpr std::array<std::string_view, PairLossK0NStages> PairLossK0StageNames = {
    "Gen pair",
    "Findable K0->pi+pi-",
    "Trigger pure reco (best collision)",
    "Trigger in TriggerTracks",
    "Trigger final selection",
    "V0 pure reco (both daughters)",
    "V0 candidate (best collision)",
    "V0 in AssocV0s",
    "V0 final selection",
    "Final reconstructed pair"};

  struct PairLossTrackInfo {
    int64_t globalIndex = -1;
    float pt = 0.0f;
    float eta = 0.0f;
    float phi = 0.0f;
    int sign = 0;
    int tpcCrossedRows = 0;
    int tpcSharedClusters = 0;
    int itsClusters = 0;
    bool hasLayer0 = false;
    // Needed to replicate the declarative preFilterTracks DCAxy cut of
    // hStrangeCorrelationFilter.cxx, which uses exactly these two columns.
    float dcaXY = 0.0f;
    float signed1Pt = 0.0f;
    int64_t collisionId = -1;
  };

  // First-failing-condition breakdown for the stage3->4 gate ("Trigger track, best collision"
  // -> "Trigger in TriggerTracks"). Order mirrors the order in which hStrangeCorrelationFilter.cxx
  // applies its conditions -- the declarative preFilterTracks DCAxy cut first, then the early
  // returns of isValidTrigger() -- so "first reason to fail" means the same thing here.
  //
  // Pass/fail itself is decided by the real TriggerTracks table, not by the replica below;
  // the replica only attributes a reason once the table has already said "not in". That keeps
  // the Passed bin numerically identical to stage 4 even when the replica's configurables have
  // drifted away from the filter's, and routes any such drift into the two dedicated bins.
  enum PairLossTriggerTracksFailureReason : int {
    PairLossTriggerTracksPassed = 0,
    // First, because preFilterTracks is applied to the track table before isValidTrigger() ever
    // sees the track, so DCAxy is the earliest condition a trigger candidate can fail.
    PairLossTriggerTracksFailDcaXY,
    PairLossTriggerTracksFailEta,
    PairLossTriggerTracksFailPt,
    PairLossTriggerTracksFailCrossedRows,
    PairLossTriggerTracksFailITS,
    PairLossTriggerTracksFailSharedClusters,
    PairLossTriggerTracksFailLayer0,
    // Really not in the TriggerTracks table, yet the replica below claims that at least one
    // candidate track should have passed. The only possible source is a mismatch between this
    // task's triggerTracks* configurables and hStrangeCorrelationFilter.cxx, so this bin acts
    // directly as a "are the configurations aligned?" monitor: it should be zero.
    PairLossTriggerTracksNotInTableUnexplained,
    // The opposite drift: really in the table, but the replica would reject it.
    PairLossTriggerTracksInTableButCutRejects,
    PairLossTriggerTracksNReasons
  };

  static constexpr std::array<std::string_view, PairLossTriggerTracksNReasons> PairLossTriggerTracksFailureNames = {
    "Passed (really in TriggerTracks)",
    "Failed DCAxy prefilter",
    "Failed eta window",
    "Failed pT window",
    "Failed min TPC crossed rows",
    "Failed hasITS requirement",
    "Failed max TPC shared clusters",
    "Failed ITS layer-0 requirement",
    "Not in table, no cut explains it",
    "In table, but replica would reject"};

  struct PairLossV0Info {
    int64_t globalIndex = -1;
    int64_t positiveTrackId = -1;
    int64_t negativeTrackId = -1;
    float pt = 0.0f;
    float eta = 0.0f;
    float phi = 0.0f;
    float radius = 0.0f;
    float cosPA = 0.0f;
    float dcaDaughters = 0.0f;
    float massNSigma = 0.0f;
    int64_t collisionId = -1;
  };

  struct PairLossTruthTrackInfo {
    int64_t globalIndex = -1;
    float pt = 0.0f;
    float eta = 0.0f;
    float phi = 0.0f;
    int sign = 0;
    bool physicalPrimary = false;
  };

  // One object of a GenStudy h-K0 pair: generated kinematics plus whether it has a
  // reconstructed counterpart, in exactly the sense the GenStudy single-particle
  // folders use.
  struct GenStudyPairObject {
    float pt = 0.0f;
    float eta = 0.0f;
    float phi = 0.0f;
    int64_t globalIndex = -1;
    int64_t motherIndex = -1;
    bool reconstructed = false;
  };

  struct PairLossTruthK0Info {
    int64_t globalIndex = -1;
    float pt = 0.0f;
    float eta = 0.0f;
    float phi = 0.0f;
    float decayRadius = -1.0f;
    bool findable = false;
    bool physicalPrimary = false;
    PairLossTruthTrackInfo positiveDaughter;
    PairLossTruthTrackInfo negativeDaughter;
  };

  struct PairLossDeltaPhiStarInfo {
    bool valid = false;
    float minAbs = std::numeric_limits<float>::max();
    float signedAtMin = 0.0f;
    float radiusAtMin = -1.0f;
  };

  using PairLossTrackMap = std::unordered_map<int64_t, std::vector<PairLossTrackInfo>>;
  using PairLossV0Map = std::unordered_map<int64_t, std::vector<PairLossV0Info>>;

  struct PairLossPairKey {
    int64_t triggerId = -1;
    int64_t k0Id = -1;
    bool operator==(PairLossPairKey const&) const = default;
  };
  struct PairLossPairKeyHash {
    size_t operator()(PairLossPairKey const& key) const
    {
      const auto h1 = std::hash<int64_t>{}(key.triggerId);
      const auto h2 = std::hash<int64_t>{}(key.k0Id);
      return h1 ^ (h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2));
    }
  };

  // Per-MC-collision context used only while the ordinary Rec implementation
  // is running inside the PairLoss Rec comparison part. It is thread-local
  // static because adding another task data member exceeds the number of
  // elements supported by the O2 task-reflection machinery.
  struct PairLossComparisonContext {
    bool active = false;
    int64_t bestCollisionId = -1;
    float vtxZ = 0.0f;
    float multiplicity = 0.0f;
    std::unordered_map<int64_t, int64_t> trackToMc;
    std::unordered_map<int64_t, int64_t> v0ToMc;
    std::unordered_map<int64_t, PairLossTruthTrackInfo> truthTriggers;
    std::unordered_map<int64_t, PairLossTruthK0Info> truthK0s;
    std::unordered_set<int64_t> primaryTruthTriggerIds;
    std::unordered_set<int64_t> primaryTruthK0Ids;
    std::unordered_set<int64_t> oldTruthTriggerIds;
    std::unordered_set<int64_t> oldTruthK0Ids;
    // Exact old-Final truth-pair keys for the current MC collision. Used at
    // Rec variant 14 to form the Rec-minus-Final set difference.
    std::unordered_set<PairLossPairKey, PairLossPairKeyHash> finalTruthPairs;
    std::unordered_set<PairLossPairKey, PairLossPairKeyHash> countedPairs;
    std::unordered_set<int64_t> countedTriggers;

    void clear()
    {
      active = false;
      bestCollisionId = -1;
      vtxZ = 0.0f;
      multiplicity = 0.0f;
      trackToMc.clear();
      v0ToMc.clear();
      truthTriggers.clear();
      truthK0s.clear();
      primaryTruthTriggerIds.clear();
      primaryTruthK0Ids.clear();
      oldTruthTriggerIds.clear();
      oldTruthK0Ids.clear();
      finalTruthPairs.clear();
      countedPairs.clear();
      countedTriggers.clear();
    }
  };
  static thread_local PairLossComparisonContext pairLossComparison;

  int mPairLossRunNumber = -1;
  double mPairLossMagneticField = 0.0;

  /// Function to aid in calculating delta-phi
  /// \param phi1 first phi value
  /// \param phi2 second phi value
  double computeDeltaPhi(double phi1, double phi2)
  {
    double deltaPhi = phi1 - phi2;
    double shiftedDeltaPhi = RecoDecay::constrainAngle(deltaPhi, -PIHalf);
    return shiftedDeltaPhi;
  }

  /// Counts how many associated-quality tracks lie inside a cone of radius localDensityConeRadius around (etaRef, phiRef), skipping the listed track indices
  template <typename TTracks>
  int computeLocalDensity(float etaRef, float phiRef, TTracks const& tracks, std::initializer_list<int64_t> skipIds)
  {
    int localDensity = 0;
    for (auto const& track : tracks) {
      if (std::find(skipIds.begin(), skipIds.end(), track.globalIndex()) != skipIds.end() || !isValidAssocHadron(track)) {
        continue;
      }
      double deltaEta = track.eta() - etaRef;
      double deltaPhi = RecoDecay::constrainAngle(track.phi() - phiRef, -PI);
      if (std::hypot(deltaEta, deltaPhi) < masterConfigurations.localDensityConeRadius) {
        localDensity++;
      }
    }
    return localDensity;
  }

  /// Collects the MC index of a particle and of its decay products, so that a particle never contributes to its own local density
  template <typename TMcParticle>
  void collectDescendantIds(TMcParticle const& mcParticle, std::vector<int64_t>& ids, int depth = 0)
  {
    ids.push_back(mcParticle.globalIndex());
    if (depth >= 3 || !mcParticle.has_daughters()) {
      return;
    }
    for (auto const& daughter : mcParticle.template daughters_as<aod::McParticles>()) {
      collectDescendantIds(daughter, ids, depth + 1);
    }
  }

  /// Generated-level counterpart: the density the reconstruction would have measured around the generated
  /// direction, i.e. the same associated-quality tracks of the same collision, so that the axis means the
  /// same thing here as in the reconstructed histograms and can be used to correct data binned in it
  template <typename TMcParticle, typename TTracks>
  int computeLocalDensityGen(TMcParticle const& mcParticle, TTracks const& tracks)
  {
    std::vector<int64_t> skipIds;
    collectDescendantIds(mcParticle, skipIds);
    int localDensity = 0;
    for (auto const& track : tracks) {
      if (!isValidAssocHadron(track)) {
        continue;
      }
      if (track.has_mcParticle() && std::find(skipIds.begin(), skipIds.end(), track.mcParticleId()) != skipIds.end()) {
        continue;
      }
      double deltaEta = track.eta() - mcParticle.eta();
      double deltaPhi = RecoDecay::constrainAngle(track.phi() - mcParticle.phi(), -PI);
      if (std::hypot(deltaEta, deltaPhi) < masterConfigurations.localDensityConeRadius) {
        localDensity++;
      }
    }
    return localDensity;
  }

  static bool isPairLossTriggerPdg(int pdgCode)
  {
    const int absPdg = std::abs(pdgCode);
    return absPdg == PDG_t::kPiPlus || absPdg == PDG_t::kKPlus || absPdg == PDG_t::kProton || absPdg == PDG_t::kElectron || absPdg == PDG_t::kMuonMinus;
  }

  template <typename TTrack>
  PairLossTrackInfo makePairLossTrackInfo(TTrack const& track) const
  {
    return PairLossTrackInfo{
      .globalIndex = static_cast<int64_t>(track.globalIndex()),
      .pt = track.pt(),
      .eta = track.eta(),
      .phi = track.phi(),
      .sign = track.sign(),
      .tpcCrossedRows = track.tpcNClsCrossedRows(),
      .tpcSharedClusters = track.tpcNClsShared(),
      .itsClusters = track.itsNCls(),
      .hasLayer0 = static_cast<bool>(TESTBIT(track.itsClusterMap(), 0)),
      .dcaXY = track.dcaXY(),
      .signed1Pt = track.signed1Pt(),
      .collisionId = static_cast<int64_t>(track.collisionId())};
  }

  // Replicates preFilterTracks and isValidTrigger() from hStrangeCorrelationFilter.cxx
  // condition-by-condition (same order) so that, for a trigger already known to exist in the best
  // collision (stage 3), we can tell which single cut is responsible for it not making it
  // into the TriggerTracks table (stage 4), instead of only knowing that it failed overall.
  int classifyTriggerTracksFailure(PairLossTrackInfo const& info)
  {
    // Negation of "Filter preFilterTracks = nabs(dcaXY) < dcaXYconstant + dcaXYpTdep * nabs(signed1Pt)":
    // same columns, same formula, and >= so that the boundary is rejected exactly as the filter does.
    if (std::abs(info.dcaXY) >= pairLossK0Configurations.triggerTracksDcaXYconstant + pairLossK0Configurations.triggerTracksDcaXYpTdep * std::abs(info.signed1Pt)) {
      return PairLossTriggerTracksFailDcaXY;
    }
    if (info.eta > pairLossK0Configurations.triggerTracksEtaMax || info.eta < pairLossK0Configurations.triggerTracksEtaMin) {
      return PairLossTriggerTracksFailEta;
    }
    if (info.pt > pairLossK0Configurations.triggerTracksPtMax || info.pt < pairLossK0Configurations.triggerTracksPtMin) {
      return PairLossTriggerTracksFailPt;
    }
    if (info.tpcCrossedRows < pairLossK0Configurations.triggerTracksMinCrossedRows) {
      return PairLossTriggerTracksFailCrossedRows;
    }
    if (info.itsClusters <= 0 && pairLossK0Configurations.triggerTracksRequireITS) {
      return PairLossTriggerTracksFailITS;
    }
    if (info.tpcSharedClusters > pairLossK0Configurations.triggerTracksMaxSharedClusters) {
      return PairLossTriggerTracksFailSharedClusters;
    }
    if (!info.hasLayer0 && pairLossK0Configurations.triggerTracksRequireLayer0) {
      return PairLossTriggerTracksFailLayer0;
    }
    return PairLossTriggerTracksPassed;
  }

  PairLossDeltaPhiStarInfo calculateMinimumDeltaPhiStar(PairLossTruthTrackInfo const& trigger, PairLossTruthTrackInfo const& daughter, double magneticField, float decayRadiusCm)
  {
    PairLossDeltaPhiStarInfo result;
    if (trigger.pt <= 0.0f || daughter.pt <= 0.0f || pairLossK0Configurations.radiusStep <= 0.0f) {
      return result;
    }

    const double triggerPhase = (-0.3 * magneticField * trigger.sign) / (2.0 * trigger.pt);
    const double daughterPhase = (-0.3 * magneticField * daughter.sign) / (2.0 * daughter.pt);
    const double deltaPhi = daughter.phi - trigger.phi;
    const double startRadius = std::max(static_cast<double>(pairLossK0Configurations.minRadius), static_cast<double>(std::max(decayRadiusCm, 0.0f)) / 100.0);
    const double maxRadius = pairLossK0Configurations.maxRadius;
    const double radiusStep = pairLossK0Configurations.radiusStep;

    for (double radius = startRadius; radius <= maxRadius + 0.5 * radiusStep; radius += radiusStep) {
      const double daughterArgument = daughterPhase * radius;
      const double triggerArgument = triggerPhase * radius;
      if (std::abs(daughterArgument) > 1.0 || std::abs(triggerArgument) > 1.0) {
        continue;
      }
      const double signedDeltaPhiStar = RecoDecay::constrainAngle(deltaPhi + std::asin(daughterArgument) - std::asin(triggerArgument), -PI);
      const float absDeltaPhiStar = std::abs(signedDeltaPhiStar);
      if (absDeltaPhiStar < result.minAbs) {
        result.valid = true;
        result.minAbs = absDeltaPhiStar;
        result.signedAtMin = signedDeltaPhiStar;
        result.radiusAtMin = radius;
      }
    }
    return result;
  }

  /// Function to load zorro
  /// \param bc provided such that the run number + timestamp can be used
  void initZorro(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumberZorro == bc.runNumber()) {
      return;
    }

    zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), zorroMask.value);
    zorro.populateHistRegistry(histos, bc.runNumber());

    mRunNumberZorro = bc.runNumber();
  }

  /// Function to load efficiencies to memory from CCDB
  /// \param bc provided such that the run number can be used
  void initEfficiencyFromCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    mRunNumber = bc.runNumber();
    LOG(info) << "Loading efficiencies from CCDB for run " << mRunNumber << " now...";
    auto timeStamp = bc.timestamp();

    auto* listEfficiencies = ccdb->getForTimeStamp<TList>(efficiencyCCDBPath, timeStamp);

    if (!listEfficiencies) {
      LOG(fatal) << "Problem getting TList object with efficiencies!";
    }

    hEfficiencyTrigger = dynamic_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyTrigger"));
    hEfficiencyTriggerMult = dynamic_cast<TH3F*>(listEfficiencies->FindObject("hEfficiencyTriggerMult"));
    hEfficiencyTriggerMultVsPhi = dynamic_cast<THnT<float>*>(listEfficiencies->FindObject("hEfficiencyTriggerMultVsPhi"));
    hEfficiencyK0Short = dynamic_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyK0Short"));
    hEfficiencyK0ShortMultVsPhi = dynamic_cast<THnT<float>*>(listEfficiencies->FindObject("hEfficiencyK0ShortMultVsPhi"));
    hEfficiencyLambda = dynamic_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyLambda"));
    hEfficiencyLambdaMultVsPhi = dynamic_cast<THnT<float>*>(listEfficiencies->FindObject("hEfficiencyLambdaMultVsPhi"));
    hEfficiencyAntiLambda = dynamic_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyAntiLambda"));
    hEfficiencyAntiLambdaMultVsPhi = dynamic_cast<THnT<float>*>(listEfficiencies->FindObject("hEfficiencyAntiLambdaMultVsPhi"));
    hEfficiencyXiMinus = dynamic_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyXiMinus"));
    hEfficiencyXiMinusMultVsPhi = dynamic_cast<THnT<float>*>(listEfficiencies->FindObject("hEfficiencyXiMinusMultVsPhi"));
    hEfficiencyXiPlus = dynamic_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyXiPlus"));
    hEfficiencyXiPlusMultVsPhi = dynamic_cast<THnT<float>*>(listEfficiencies->FindObject("hEfficiencyXiPlusMultVsPhi"));
    hEfficiencyOmegaMinus = dynamic_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyOmegaMinus"));
    hEfficiencyOmegaMinusMultVsPhi = dynamic_cast<THnT<float>*>(listEfficiencies->FindObject("hEfficiencyOmegaMinusMultVsPhi"));
    hEfficiencyOmegaPlus = dynamic_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyOmegaPlus"));
    hEfficiencyOmegaPlusMultVsPhi = dynamic_cast<THnT<float>*>(listEfficiencies->FindObject("hEfficiencyOmegaPlusMultVsPhi"));
    hEfficiencyHadron = dynamic_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyHadron"));
    hEfficiencyHadronMult = dynamic_cast<TH3F*>(listEfficiencies->FindObject("hEfficiencyHadronMult"));
    hEfficiencyPion = dynamic_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyPion"));
    hPurityHadron = dynamic_cast<TH1F*>(listEfficiencies->FindObject("hPurityHadron"));
    hPurityHadronMult = dynamic_cast<TH2F*>(listEfficiencies->FindObject("hPurityHadronMult"));
    hEfficiencyUncertaintyTrigger = dynamic_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyUncertaintyTrigger"));
    hEfficiencyUncertaintyTriggerMult = dynamic_cast<TH3F*>(listEfficiencies->FindObject("hEfficiencyUncertaintyTriggerMult"));
    hEfficiencyUncertaintyK0Short = dynamic_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyUncertaintyK0Short"));
    hEfficiencyUncertaintyLambda = dynamic_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyUncertaintyLambda"));
    hEfficiencyUncertaintyAntiLambda = dynamic_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyUncertaintyAntiLambda"));
    hEfficiencyUncertaintyXiMinus = dynamic_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyUncertaintyXiMinus"));
    hEfficiencyUncertaintyXiPlus = dynamic_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyUncertaintyXiPlus"));
    hEfficiencyUncertaintyOmegaMinus = dynamic_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyUncertaintyOmegaMinus"));
    hEfficiencyUncertaintyOmegaPlus = dynamic_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyUncertaintyOmegaPlus"));
    hEfficiencyUncertaintyPion = dynamic_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyUncertaintyPion"));
    hEfficiencyUncertaintyHadron = dynamic_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyUncertaintyHadron"));
    hEfficiencyUncertaintyHadronMult = dynamic_cast<TH3F*>(listEfficiencies->FindObject("hEfficiencyUncertaintyHadronMult"));
    hPurityUncertaintyHadron = dynamic_cast<TH1F*>(listEfficiencies->FindObject("hPurityUncertaintyHadron"));
    hPurityUncertaintyHadronMult = dynamic_cast<TH2F*>(listEfficiencies->FindObject("hPurityUncertaintyHadronMult"));
    if (efficiencyFlags.applyEfficiencyPropagation && !efficiencyFlags.applyEffAsFunctionOfMultAndPhi && !hEfficiencyUncertaintyTrigger) {
      LOG(fatal) << "Problem getting hEfficiencyUncertaintyTrigger!";
    }
    if (efficiencyFlags.applyEffAsFunctionOfMult && !hEfficiencyTriggerMult) {
      LOG(fatal) << "Problem getting hEfficiencyTriggerMult!";
    }
    LOG(info) << "Efficiencies now loaded for " << mRunNumber;
  }

  template <typename TV0>
  uint64_t v0selectionBitmap(TV0 const& v0, float pvx, float pvy, float pvz)
  // precalculate this information so that a check is one mask operation, not many
  {
    uint64_t bitMap = 0;
    // proper lifetime , DCA daughter to prim.vtx
    if (masterConfigurations.doCorrelationK0Short) {
      // proper lifetime
      if (v0.distovertotmom(pvx, pvy, pvz) * o2::constants::physics::MassK0Short < v0Selection.lifetimecutK0S) {
        SETBIT(bitMap, 0);
      }
      // DCA daughter to prim.vtx and armenteros
      float dcaCut = v0Selection.dcaDaugToPVForK0s == 0 ? v0Selection.dcaMesonToPV : v0Selection.dcaDaugToPVForK0s;
      if (std::abs(v0.dcapostopv()) > dcaCut && std::abs(v0.dcanegtopv()) > dcaCut && v0.qtarm() * v0Selection.armPodCut > std::abs(v0.alpha())) {
        SETBIT(bitMap, 3);
      }
    }
    if (masterConfigurations.doCorrelationLambda) {
      // proper lifetime
      if (v0.distovertotmom(pvx, pvy, pvz) * o2::constants::physics::MassLambda0 < v0Selection.lifetimecutLambda) {
        SETBIT(bitMap, 1);
      }
      // DCA daughter to prim.vtx
      if (std::abs(v0.dcapostopv()) > v0Selection.dcaBaryonToPV && std::abs(v0.dcanegtopv()) > v0Selection.dcaMesonToPV) {
        SETBIT(bitMap, 4);
      }
    }
    if (masterConfigurations.doCorrelationAntiLambda) {
      // proper lifetime
      if (v0.distovertotmom(pvx, pvy, pvz) * o2::constants::physics::MassLambda0 < v0Selection.lifetimecutLambda) {
        SETBIT(bitMap, 2);
      }
      // DCA daughter to prim.vtx
      if (std::abs(v0.dcapostopv()) > v0Selection.dcaMesonToPV && std::abs(v0.dcanegtopv()) > v0Selection.dcaBaryonToPV) {
        SETBIT(bitMap, 5);
      }
    }
    return bitMap;
  }

  template <typename TCascade>
  uint64_t cascadeselectionBitmap(TCascade const& casc, float pvx, float pvy, float pvz)
  {
    uint64_t bitMap = 0;
    float cascpos = std::hypot(casc.x() - pvx, casc.y() - pvy, casc.z() - pvz);
    float cascptotmom = std::hypot(casc.px(), casc.py(), casc.pz());
    float ctauXi = o2::constants::physics::MassXiMinus * cascpos / ((cascptotmom + 1e-13) * ctauxi);
    float ctauOmega = o2::constants::physics::MassOmegaMinus * cascpos / ((cascptotmom + 1e-13) * ctauomega);

    bool isGoodDCANegCasc = std::abs(casc.dcabachtopv()) > cascadeSelections.dcaBachToPV && std::abs(casc.dcapostopv()) > cascadeSelections.dcaBaryonToPV &&
                            std::abs(casc.dcanegtopv()) > cascadeSelections.dcaMesonToPV;
    bool isGoodDCAPosCasc = std::abs(casc.dcabachtopv()) > cascadeSelections.dcaBachToPV && std::abs(casc.dcapostopv()) > cascadeSelections.dcaMesonToPV &&
                            std::abs(casc.dcanegtopv()) > cascadeSelections.dcaBaryonToPV;
    // TPC PID and DCA daughter to prim.vtx and comopeting casc.rej and life time
    if (masterConfigurations.doCorrelationXiMinus) {
      // DCA daughter to prim.vtx
      if (isGoodDCANegCasc) {
        SETBIT(bitMap, 0);
      }
      // comopeting casc.rej
      if (std::abs(casc.mOmega() - o2::constants::physics::MassOmegaMinus) > cascadeSelections.rejcomp) {
        SETBIT(bitMap, 4);
      }
      if (ctauXi < cascadeSelections.proplifetime) {
        SETBIT(bitMap, 8);
      }
      // y cut
      if (std::abs(casc.yXi()) < cascadeSelections.rapCut) {
        SETBIT(bitMap, 12);
      }
    }
    if (masterConfigurations.doCorrelationXiPlus) {
      // DCA daughter to prim.vtx
      if (isGoodDCAPosCasc) {
        SETBIT(bitMap, 1);
      }
      // comopeting casc.rej
      if (std::abs(casc.mOmega() - o2::constants::physics::MassOmegaMinus) > cascadeSelections.rejcomp) {
        SETBIT(bitMap, 5);
      }
      // life time
      if (ctauXi < cascadeSelections.proplifetime) {
        SETBIT(bitMap, 9);
      }
      // y cut
      if (std::abs(casc.yXi()) > cascadeSelections.rapCut) {
        SETBIT(bitMap, 13);
      }
    }
    if (masterConfigurations.doCorrelationOmegaMinus) {
      // DCA daughter to prim.vtx
      if (isGoodDCANegCasc) {
        SETBIT(bitMap, 2);
      }
      // comopeting casc.rej
      if (std::abs(casc.mXi() - o2::constants::physics::MassXiMinus) > cascadeSelections.rejcomp) {
        SETBIT(bitMap, 6);
      }
      // life time
      if (ctauOmega < cascadeSelections.proplifetime) {
        SETBIT(bitMap, 10);
      }
      // y cut
      if (std::abs(casc.yOmega()) < cascadeSelections.rapCut) {
        SETBIT(bitMap, 14);
      }
    }
    if (masterConfigurations.doCorrelationOmegaPlus) {
      // DCA daughter to prim.vtx
      if (isGoodDCAPosCasc) {
        SETBIT(bitMap, 3);
      }
      // comopeting casc.rej
      if (std::abs(casc.mXi() - o2::constants::physics::MassXiMinus) > cascadeSelections.rejcomp) {
        SETBIT(bitMap, 7);
      }
      // life time
      if (ctauOmega < cascadeSelections.proplifetime) {
        SETBIT(bitMap, 11);
      }
      // y cut
      if (std::abs(casc.yOmega()) > cascadeSelections.rapCut) {
        SETBIT(bitMap, 15);
      }
    }
    return bitMap;
  }

  template <class TTrack>
  bool isValidTrigger(TTrack const& track, bool isLeading)
  {
    if (track.tpcNClsCrossedRows() < trackSelection.minTPCNCrossedRowsTrigger) {
      return false; // crossed rows
    }
    if (!track.hasITS() && trackSelection.triggerRequireITS) {
      return false; // skip, doesn't have ITS signal (skips lots of TPC-only!)
    }
    if (track.tpcNClsShared() > trackSelection.triggerMaxTPCSharedClusters) {
      return false; // skip, has shared clusters
    }
    if (!(TESTBIT(track.itsClusterMap(), 0)) && trackSelection.triggerRequireL0) {
      return false; // skip, doesn't have cluster in ITS L0
    }
    // systematic variations: trigger DCAxy
    if (std::abs(track.dcaXY()) > trackSelection.dcaXYconstant + trackSelection.dcaXYpTdep * std::abs(track.signed1Pt())) {
      return false;
    }
    // systematic variations: trigger DCAz
    if (trackSelection.requireDCAzCut && std::abs(track.dcaZ()) > trackSelection.dcaZconstant + trackSelection.dcaZpTdep * std::abs(track.signed1Pt())) {
      return false;
    }
    if (track.pt() > axisRanges[3][1] || track.pt() < axisRanges[3][0]) {
      return false;
    }
    if (triggerParticleCharge > 0 && track.sign() < 0) {
      return false;
    }
    if (triggerParticleCharge < 0 && track.sign() > 0) {
      return false;
    }
    if (useTheLeadingParticleAsTrigger && !isLeading) {
      return false;
    }
    if (trackSelection.doITSChi2Selection && track.itsChi2NCl() > trackSelection.chi2ITSfit) {
      return false;
    }
    return true;
  }
  template <class TTrack>
  bool isValidAssocHadron(TTrack const& track)
  {
    if (track.tpcNClsCrossedRows() < trackSelection.minTPCNCrossedRowsAssociated) {
      return false; // crossed rows
    }
    if (!track.hasITS() && trackSelection.assocRequireITS) {
      return false; // skip, doesn't have ITS signal (skips lots of TPC-only!)
    }
    if (track.tpcNClsShared() > trackSelection.assocMaxTPCSharedClusters) {
      return false; // skip, has shared clusters
    }
    if (!(TESTBIT(track.itsClusterMap(), 0)) && trackSelection.assocRequireL0) {
      return false; // skip, doesn't have cluster in ITS L0
    }
    // systematic variations: trigger DCAxy
    if (std::abs(track.dcaXY()) > trackSelection.dcaXYconstantAssoc + trackSelection.dcaXYpTdepAssoc * std::abs(track.signed1Pt())) {
      return false;
    }
    // systematic variations: trigger DCAz
    if (trackSelection.requireDCAzCut && std::abs(track.dcaZ()) > trackSelection.dcaZconstantAssoc + trackSelection.dcaZpTdepAssoc * std::abs(track.signed1Pt())) {
      return false;
    }
    if (track.pt() > axisRanges[2][1] || track.pt() < axisRanges[2][0]) {
      return false;
    }
    if (trackSelection.doITSChi2Selection && track.itsChi2NCl() > trackSelection.chi2ITSfit) {
      return false;
    }
    return true;
  }
  // V0selection in PbPb
  template <typename TV0>
  bool v0SelectedPbPb(TV0 const& v0)
  {
    // v0radius
    if (v0.v0radius() < v0Selection.v0RadiusMin) {
      return false;
    }
    if (v0.v0radius() > v0Selection.v0RadiusMax) {
      return false;
    }
    // v0cosPA
    if (v0.v0cosPA() < v0Selection.v0cospa) {
      return false;
    }
    // dcaV0daughters
    if (v0.dcaV0daughters() > v0Selection.dcaV0dau) {
      return false;
    }
    return true;
  }

  // cascadeselection in PbPb
  template <typename TCascade>
  bool cascadeSelectedPbPb(TCascade const& casc, float pvx, float pvy, float pvz)
  {
    // bachBaryonCosPA
    if (casc.bachBaryonCosPA() < cascadeSelections.bachBaryonCosPA) {
      return false;
    }
    // bachBaryonDCAxyToPV
    if (std::abs(casc.bachBaryonDCAxyToPV()) > cascadeSelections.bachBaryonDCAxyToPV) {
      return false;
    }
    // casccosPA
    if (casc.casccosPA(pvx, pvy, pvz) < cascadeSelections.cascCospa) {
      return false;
    }
    // dcacascdaughters
    float ptDepCut = cascadeSelections.dcaCacsDauPar0;
    if (casc.pt() > cascadeSelections.lowPtForCascDaugPtDep && casc.pt() < cascadeSelections.highPtForCascDaugPtDep) {
      ptDepCut = cascadeSelections.dcaCacsDauPar1;
    } else if (casc.pt() > cascadeSelections.highPtForCascDaugPtDep) {
      ptDepCut = cascadeSelections.dcaCacsDauPar2;
    }
    if (casc.dcacascdaughters() > ptDepCut) {
      return false;
    }
    // dcaV0daughters
    if (casc.dcaV0daughters() > cascadeSelections.cascdcaV0dau) {
      return false;
    }
    // dcav0topv
    if (std::abs(casc.dcav0topv(pvx, pvy, pvz)) < cascadeSelections.cascdcaV0ToPV) {
      return false;
    }
    // cascradius
    if (casc.cascradius() < cascadeSelections.cascRadius) {
      return false;
    }
    // v0radius
    if (casc.v0radius() < cascadeSelections.cascv0RadiusMin) {
      return false;
    }
    // v0cosPA
    if (casc.v0cosPA(pvx, pvy, pvz) < cascadeSelections.cascv0cospa) {
      return false;
    }
    // lambdaMassWin
    if (std::abs(casc.mLambda() - o2::constants::physics::MassLambda0) > cascadeSelections.cascV0masswindow) {
      return false;
    }
    return true;
  }
  double calculateAverageDeltaPhiStar(std::array<double, 3> const& trigg, std::array<double, 3> const& assoc, double B)
  {
    double dPhiStar = 0;
    double dPhiStarMean = 0;

    double dPhi = assoc[0] - trigg[0];
    double phaseProton = (-0.3 * B * assoc[2]) / (2 * assoc[1]);
    double phaseTrack = (-0.3 * B * trigg[2]) / (2 * trigg[1]);

    for (double r = MinRadiusTPC; r <= MaxRadiusTPC; r += 0.05) {
      dPhiStar = dPhi + std::asin(phaseProton * r) - std::asin(phaseTrack * r);
      dPhiStarMean += (dPhiStar / 34);
    }

    return dPhiStarMean;
  }
  void fillTriggerHistogram(std::shared_ptr<TH2> const& hist, double pt, double mult, float eff, float effUncert, float purity, float purityErr)
  {
    int binx = hist->GetXaxis()->FindBin(pt);
    int biny = hist->GetYaxis()->FindBin(mult);
    float previousContent = hist->GetBinContent(binx, biny);
    float previousUncert = hist->GetBinError(binx, biny);
    float newContent = previousContent + purity / eff;
    float newUncert = std::sqrt(previousUncert * previousUncert + std::pow(purity / eff, 2) + std::pow(purityErr / eff, 2) + std::pow(effUncert, 2) / std::pow(eff, 4));
    hist->SetBinContent(binx, biny, newContent);
    hist->SetBinError(binx, biny, newUncert);
  }
  void fillCorrelationHistogram(std::shared_ptr<THn> const& hist, std::array<double, 6> const& binFillThn, float etaWeight, float efficiency, float totalEffUncert, float purity, float totalPurityUncert)
  {
    float previousContent = 0.0f, previousError2 = 0.0f, currentContent = 0.0f, currentError2 = 0.0f;
    int bin = hist->GetBin(binFillThn.data());
    previousContent = hist->GetBinContent(bin);
    previousError2 = hist->GetBinError2(bin);
    currentContent = previousContent + etaWeight * purity / (efficiency);
    currentError2 = previousError2 + std::pow(etaWeight * purity / (efficiency), 2) + std::pow(etaWeight * totalPurityUncert / (efficiency), 2) + std::pow(totalEffUncert * purity * etaWeight, 2) / std::pow(efficiency, 4);
    hist->SetBinContent(bin, currentContent);
    hist->SetBinError2(bin, currentError2);
  }
  void fillCorrelationsV0(aod::TriggerTracks const& triggers, aod::AssocV0s const& assocs, bool mixing, bool mixingInBf, float pvx, float pvy, float pvz, float mult, double bField)
  {
    ValidCollision currentCollision;
    int binMult = 0;
    int nBinsMult = 0;
    int nBinsVtxZ = 0;
    int binVtxZ = 0;
    currentCollision.pvz = pvz;
    currentCollision.mult = mult;
    if (mixingInBf) {
      nBinsMult = histos.get<TH1>(HIST("axes/hMultAxis"))->GetNbinsX();
      binMult = histos.get<TH1>(HIST("axes/hMultAxis"))->GetXaxis()->FindBin(mult) - 1;
      nBinsVtxZ = histos.get<TH1>(HIST("axes/hVertexZAxis"))->GetNbinsX();
      binVtxZ = histos.get<TH1>(HIST("axes/hVertexZAxis"))->GetXaxis()->FindBin(pvz) - 1;
    }
    bool firstLoop = false;
    for (auto const& triggerTrack : triggers) {
      if (masterConfigurations.doTriggPhysicalPrimary && !triggerTrack.mcPhysicalPrimary()) {
        continue;
      }
      auto trigg = triggerTrack.track_as<TracksComplete>();
      if (!isValidTrigger(trigg, triggerTrack.isLeading())) {
        continue;
      }
      if (trackSelection.checkForITSTPCMissmatchMC && bitcheck(triggerTrack.mcMask(), 13)) {
        continue;
      }
      float efficiencyTrigg = 1.0f;
      float efficiencyTriggError = 0.0f;
      float purityTrigg = 1.0f;
      float purityTriggErr = 0.0;
      std::array<double, 4> bintrig = {trigg.pt(), trigg.eta(), trigg.phi(), mult};
      if (efficiencyFlags.applyEfficiencyForTrigger) {
        if (!efficiencyFlags.applyEffAsFunctionOfMultAndPhi) {
          efficiencyTrigg = hEfficiencyTrigger->Interpolate(trigg.pt(), trigg.eta());
          if (efficiencyFlags.applyPurityTrigger) {
            purityTrigg = hPurityHadron->Interpolate(trigg.pt());
          }
          if (efficiencyFlags.applyEfficiencyPropagation) {
            efficiencyTriggError = hEfficiencyUncertaintyTrigger->Interpolate(trigg.pt(), trigg.eta());
            if (efficiencyFlags.applyPurityTrigger) {
              purityTriggErr = hPurityHadron->Interpolate(trigg.pt());
            }
          }
        } else {
          efficiencyTrigg = hEfficiencyTriggerMultVsPhi->GetBinContent(hEfficiencyTriggerMultVsPhi->GetBin(bintrig.data()));
          if (efficiencyFlags.applyEfficiencyPropagation) {
            efficiencyTriggError = hEfficiencyTriggerMultVsPhi->GetBinError(hEfficiencyTriggerMultVsPhi->GetBin(bintrig.data()));
          }
        }
        if (efficiencyTrigg == 0) { // check for zero efficiency, do not apply if the case
          efficiencyTrigg = 1;
          efficiencyTriggError = 0;
        }
      }
      if (!mixing) {

        fillTriggerHistogram(histos.get<TH2>(HIST("sameEvent/TriggerParticlesV0")), trigg.pt(), mult, efficiencyTrigg, efficiencyTriggError, purityTrigg, purityTriggErr);

        // Trigger-side control chain, filled from the exact Rec trigger loop so
        // that it is the same denominator the ordinary Rec correlation uses.
        // Each pair variant needs the matching trigger variant, otherwise only
        // shapes and not per-trigger amplitudes can be compared.
        if (pairLossComparison.active) {
          const double triggerWeight = purityTrigg / efficiencyTrigg;
          auto fillTriggerVariant = [&](int variant, double fillPtTrigger, double fillVtxZ, double fillMult, double weight) {
            histos.fill(HIST("PairLossK0/Comparison/TriggerVariants"), variant, fillPtTrigger, fillVtxZ, fillMult, weight);
          };
          fillTriggerVariant(0, trigg.pt(), pvz, mult, triggerWeight);
          const auto triggerLadderMcIt = pairLossComparison.trackToMc.find(trigg.globalIndex());
          if (triggerLadderMcIt != pairLossComparison.trackToMc.end()) {
            const auto triggerLadderTruthIt = pairLossComparison.truthTriggers.find(triggerLadderMcIt->second);
            if (triggerLadderTruthIt != pairLossComparison.truthTriggers.end()) {
              auto const& ladderTruthTrigger = triggerLadderTruthIt->second;
              fillTriggerVariant(1, trigg.pt(), pvz, mult, triggerWeight);
              if (pairLossComparison.primaryTruthTriggerIds.contains(ladderTruthTrigger.globalIndex)) {
                fillTriggerVariant(2, trigg.pt(), pvz, mult, triggerWeight);
                if (pairLossComparison.oldTruthTriggerIds.contains(ladderTruthTrigger.globalIndex)) {
                  fillTriggerVariant(3, trigg.pt(), pvz, mult, triggerWeight);
                  if (static_cast<int64_t>(triggerTrack.collisionId()) == pairLossComparison.bestCollisionId) {
                    fillTriggerVariant(4, trigg.pt(), pvz, mult, triggerWeight);
                    fillTriggerVariant(5, ladderTruthTrigger.pt, pvz, mult, triggerWeight);
                    fillTriggerVariant(6, ladderTruthTrigger.pt, pairLossComparison.vtxZ, pairLossComparison.multiplicity, triggerWeight);
                    fillTriggerVariant(7, ladderTruthTrigger.pt, pairLossComparison.vtxZ, pairLossComparison.multiplicity, 1.0);
                    if (pairLossComparison.countedTriggers.insert(ladderTruthTrigger.globalIndex).second) {
                      fillTriggerVariant(8, ladderTruthTrigger.pt, pairLossComparison.vtxZ, pairLossComparison.multiplicity, 1.0);
                    }
                  }
                }
              }
            }
          }
        }
      }

      double triggSign = trigg.sign();
      std::array<double, 3> triggForDeltaPhiStar = {trigg.phi(), trigg.pt(), triggSign};

      if (mixingInBf) {
        currentCollision.addValidParticle(trigg.eta(), trigg.phi(), trigg.pt(), -1, efficiencyTrigg, efficiencyTriggError, -1);
        if (firstLoop) {
          continue;
        }
      }
      for (auto const& assocCandidate : assocs) {
        firstLoop = true;
        auto assoc = assocCandidate.v0Core_as<V0DatasWithoutTrackX>();

        //---] syst cuts [---
        if ((masterConfigurations.doPPAnalysis && (assoc.v0radius() < v0Selection.v0RadiusMin || assoc.v0radius() > v0Selection.v0RadiusMax ||
                                                   std::abs(assoc.dcapostopv()) < v0Selection.dcapostopv || std::abs(assoc.dcanegtopv()) < v0Selection.dcanegtopv ||
                                                   assoc.v0cosPA() < v0Selection.v0cospa || assoc.dcaV0daughters() > v0Selection.dcaV0dau))) {
          continue;
        }

        if (!masterConfigurations.doPPAnalysis && !v0SelectedPbPb(assoc)) {
          continue;
        }

        uint64_t selMap = v0selectionBitmap(assoc, pvx, pvy, pvz);

        //---] removing autocorrelations [---
        auto postrack = assoc.posTrack_as<TracksComplete>();
        auto negtrack = assoc.negTrack_as<TracksComplete>();
        if (doAutocorrelationRejection) {
          if (trigg.globalIndex() == postrack.globalIndex()) {
            histos.fill(HIST("hNumberOfRejectedPairsV0"), 0.5);
            continue;
          }
          if (trigg.globalIndex() == negtrack.globalIndex()) {
            histos.fill(HIST("hNumberOfRejectedPairsV0"), 0.5);
            continue;
          }
        }

        //---] track quality check [---
        if (postrack.tpcNClsCrossedRows() < trackSelection.minTPCNCrossedRowsAssociated || negtrack.tpcNClsCrossedRows() < trackSelection.minTPCNCrossedRowsAssociated) {
          continue;
        }
        if (trackSelection.checksRequireTPCChi2 && (postrack.tpcChi2NCl() < trackSelection.minTPCChi2PerClusterAssociated || negtrack.tpcChi2NCl() < trackSelection.minTPCChi2PerClusterAssociated)) {
          continue;
        }
        if (trackSelection.requireClusterInITS && (postrack.itsNCls() < trackSelection.minITSClustersForDaughterTracks || negtrack.itsNCls() < trackSelection.minITSClustersForDaughterTracks)) {
          continue;
        }

        float deltaphi = computeDeltaPhi(trigg.phi(), assoc.phi());
        float deltaeta = trigg.eta() - assoc.eta();
        if (masterConfigurations.doMirroringInDelataEta) {
          deltaeta = std::abs(deltaeta);
        }
        float ptassoc = assoc.pt();
        float pttrigger = trigg.pt();

        // skip if basic ranges not met
        if (deltaphi < axisRanges[0][0] || deltaphi > axisRanges[0][1]) {
          continue;
        }
        if (deltaeta < axisRanges[1][0] || deltaeta > axisRanges[1][1]) {
          continue;
        }
        if (ptassoc < axisRanges[2][0] || ptassoc > axisRanges[2][1]) {
          continue;
        }

        std::array<TH2F*, AssocV0Types> hEfficiencyV0{nullptr, nullptr, nullptr};
        std::array<TH2F*, AssocV0Types> hEfficiencyUncertaintyV0{nullptr, nullptr, nullptr};
        std::array<THnF*, AssocV0Types> hEfficiencyV0MultVsPhi{nullptr, nullptr, nullptr};
        if (efficiencyFlags.applyEffAsFunctionOfMultAndPhi) {
          hEfficiencyV0MultVsPhi[0] = hEfficiencyK0ShortMultVsPhi;
          hEfficiencyV0MultVsPhi[1] = hEfficiencyLambdaMultVsPhi;
          hEfficiencyV0MultVsPhi[2] = hEfficiencyAntiLambdaMultVsPhi;
        } else {
          hEfficiencyV0[0] = hEfficiencyK0Short;
          hEfficiencyV0[1] = hEfficiencyLambda;
          hEfficiencyV0[2] = hEfficiencyAntiLambda;

          hEfficiencyUncertaintyV0[0] = hEfficiencyUncertaintyK0Short;
          hEfficiencyUncertaintyV0[1] = hEfficiencyUncertaintyLambda;
          hEfficiencyUncertaintyV0[2] = hEfficiencyUncertaintyAntiLambda;
        }

        float etaWeight = 1;
        if (checks.doOnTheFlyFlattening) {
          float preWeight = 1 - std::abs(deltaeta) / 1.6;
          etaWeight = preWeight != 0 ? 1.0f / preWeight : 1.0f;
        }

        double phiProton = postrack.phi(); // in Case of K0, both are pions, but the one in proton tagged is the positive one
        double phiPion = negtrack.phi();
        double etaProton = postrack.eta();
        double etaPion = negtrack.eta();
        double ptProton = postrack.pt();
        double ptPion = negtrack.pt();
        double signProton = postrack.sign();
        double signPion = negtrack.sign();
        if (assocCandidate.compatible(2, trackSelection.dEdxCompatibility)) {
          phiProton = negtrack.phi();
          etaProton = negtrack.eta();
          ptProton = negtrack.pt();
          signProton = negtrack.sign();
          phiPion = postrack.phi();
          etaPion = postrack.eta();
          ptPion = postrack.pt();
          signPion = postrack.sign();
        }
        std::array<double, 3> assocForDeltaPhiStar = {phiProton, ptProton, signProton};
        std::array<double, 3> assocForDeltaPhiStarPion = {phiPion, ptPion, signPion};

        static_for<0, 2>([&](auto i) {
          constexpr int Index = i.value;
          float efficiency = 1.0f;
          float totalEffUncert = 0.0;
          float efficiencyError = 0.0f;
          if (efficiencyFlags.applyEfficiencyCorrection) {
            if (efficiencyFlags.applyEffAsFunctionOfMultAndPhi) {
              std::array<double, 4> bin = {ptassoc, assoc.eta(), assoc.phi(), mult};
              efficiency = hEfficiencyV0MultVsPhi[Index]->GetBinContent(hEfficiencyV0MultVsPhi[Index]->GetBin(bin.data()));
              if (efficiencyFlags.applyEfficiencyPropagation) {
                efficiencyError = hEfficiencyV0MultVsPhi[Index]->GetBinError(hEfficiencyV0MultVsPhi[Index]->GetBin(bin.data()));
              }
            } else {
              efficiency = hEfficiencyV0[Index]->Interpolate(ptassoc, assoc.eta());
              if (efficiencyFlags.applyEfficiencyPropagation) {
                efficiencyError = hEfficiencyUncertaintyV0[Index]->Interpolate(ptassoc, assoc.eta());
              }
            }
          }
          if (efficiency == 0) { // check for zero efficiency, do not apply if the case
            efficiency = 1;
            efficiencyError = 0;
          }
          if (efficiencyFlags.applyEfficiencyPropagation) {
            totalEffUncert = std::sqrt(std::pow(efficiencyTrigg * efficiencyError, 2) + std::pow(efficiencyTriggError * efficiency, 2));
          }

          std::array<double, 6> binFillThnOppositeSignDaugher = {computeDeltaPhi(trigg.phi(), phiPion), trigg.eta() - etaPion, ptPion, pttrigger, pvz, mult};
          std::array<double, 6> binFillThnSameSignDaugher = {computeDeltaPhi(trigg.phi(), phiProton), trigg.eta() - etaProton, ptProton, pttrigger, pvz, mult};
          if ((triggSign < Neutral && !assocCandidate.compatible(2, trackSelection.dEdxCompatibility)) || (triggSign > Neutral && assocCandidate.compatible(2, trackSelection.dEdxCompatibility))) {
            std::swap(binFillThnOppositeSignDaugher[1], binFillThnSameSignDaugher[1]);
            std::swap(binFillThnOppositeSignDaugher[0], binFillThnSameSignDaugher[0]);
            std::swap(binFillThnOppositeSignDaugher[2], binFillThnSameSignDaugher[2]);
          }
          std::array<double, 6> binFillThn = {deltaphi, deltaeta, ptassoc, pttrigger, pvz, mult};
          if (TESTBIT(doCorrelation, Index) && (!efficiencyFlags.applyEfficiencyCorrection || efficiency != 0) && (masterConfigurations.doPPAnalysis || (TESTBIT(selMap, Index) && TESTBIT(selMap, Index + 3)))) {
            if (assocCandidate.compatible(Index, trackSelection.dEdxCompatibility) && (!masterConfigurations.doMCassociation || assocCandidate.mcTrue(Index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && !mixing &&
                ((-massWindowConfigurations.maxBgNSigma < assocCandidate.invMassNSigma(Index) && assocCandidate.invMassNSigma(Index) < -massWindowConfigurations.minBgNSigma) ||
                 (-massWindowConfigurations.maxPeakNSigma < assocCandidate.invMassNSigma(Index) && assocCandidate.invMassNSigma(Index) < +massWindowConfigurations.maxPeakNSigma) ||
                 (+massWindowConfigurations.minBgNSigma < assocCandidate.invMassNSigma(Index) && assocCandidate.invMassNSigma(Index) < +massWindowConfigurations.maxBgNSigma))) {
              if (std::abs(deltaphi) < 0.5) {
                histos.fill(HIST("sameEvent/InvariantMass/") + HIST(V0names[Index]) + HIST("/hNearSide"), ptassoc, pttrigger, getV0InvariantMass<Index>(assoc));
              }
              if (std::abs(PI - deltaphi) < 0.5) {
                histos.fill(HIST("sameEvent/InvariantMass/") + HIST(V0names[Index]) + HIST("/hAwaySide"), ptassoc, pttrigger, getV0InvariantMass<Index>(assoc));
              }
              if (deltaphi > 1.0 && deltaphi < 1.5) {
                histos.fill(HIST("sameEvent/InvariantMass/") + HIST(V0names[Index]) + HIST("/hUnderlyingEvent"), ptassoc, pttrigger, getV0InvariantMass<Index>(assoc));
              }
            }
            if (assocCandidate.compatible(Index, trackSelection.dEdxCompatibility) && (!masterConfigurations.doMCassociation || assocCandidate.mcTrue(Index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && !mixing && -massWindowConfigurations.maxBgNSigma < assocCandidate.invMassNSigma(Index) && assocCandidate.invMassNSigma(Index) < -massWindowConfigurations.minBgNSigma) {
              fillCorrelationHistogram(histos.get<THn>(HIST("sameEvent/LeftBg/") + HIST(V0names[Index])), binFillThn, etaWeight, efficiency * efficiencyTrigg, totalEffUncert, purityTrigg, purityTriggErr);
              if (pairLossComparison.active && Index == IndexK0) {
                // Exact Rec left sideband, with the same reconstructed
                // selections, coordinates and weight as the ordinary Rec path.
                fillCorrelationHistogram(histos.get<THn>(HIST("PairLossK0/Comparison/RecLeftBg")), binFillThn, etaWeight, efficiency * efficiencyTrigg, totalEffUncert, purityTrigg, purityTriggErr);
              }
              if (doDeltaPhiStarCheck) {
                double deltaPhiStar = calculateAverageDeltaPhiStar(triggForDeltaPhiStar, assocForDeltaPhiStar, bField);
                double deltaPhiStarPion = calculateAverageDeltaPhiStar(triggForDeltaPhiStar, assocForDeltaPhiStarPion, bField);
                if ((Index == IndexK0 && triggSign > Neutral) || (Index == IndexLambda && triggSign > Neutral) || (Index == IndexAntiLambda && triggSign < Neutral)) {
                  histos.fill(HIST("sameEvent/LeftBg/") + HIST(V0names[Index]) + HIST("DeltaPhiStar"), deltaPhiStar, trigg.eta() - etaProton, 0.5);
                  if (Index == IndexK0) {
                    histos.fill(HIST("sameEvent/LeftBg/") + HIST(V0names[Index]) + HIST("DeltaPhiStar"), deltaPhiStarPion, trigg.eta() - etaPion, -0.5);
                  }
                } else {
                  histos.fill(HIST("sameEvent/LeftBg/") + HIST(V0names[Index]) + HIST("DeltaPhiStar"), deltaPhiStar, trigg.eta() - etaProton, -0.5);
                  if (Index == IndexK0) {
                    histos.fill(HIST("sameEvent/LeftBg/") + HIST(V0names[Index]) + HIST("DeltaPhiStar"), deltaPhiStarPion, trigg.eta() - etaPion, 0.5);
                  }
                }
              }
            }
            if (assocCandidate.compatible(Index, trackSelection.dEdxCompatibility) && (!masterConfigurations.doMCassociation || assocCandidate.mcTrue(Index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && !mixing && ((-massWindowConfigurations.maxPeakNSigma < assocCandidate.invMassNSigma(Index) && assocCandidate.invMassNSigma(Index) < +massWindowConfigurations.maxPeakNSigma))) {
              if (masterConfigurations.doCorrelationsHadronV0daughter) {
                fillCorrelationHistogram(histos.get<THn>(HIST("sameEvent/Signal/") + HIST(V0names[Index]) + HIST("_hSameSign")), binFillThnSameSignDaugher, 1, 1, 1, 1, 1);
                fillCorrelationHistogram(histos.get<THn>(HIST("sameEvent/Signal/") + HIST(V0names[Index]) + HIST("_hOppositeSign")), binFillThnOppositeSignDaugher, 1, 1, 1, 1, 1);
              }
              fillCorrelationHistogram(histos.get<THn>(HIST("sameEvent/Signal/") + HIST(V0names[Index])), binFillThn, etaWeight, efficiency * efficiencyTrigg, totalEffUncert, purityTrigg, purityTriggErr);

              // PairLoss Rec control study. Stage 0 is filled from this
              // exact Rec signal branch, so it is not a hand-written
              // approximation of the reconstructed pair selection. Each next
              // stage adds exactly one cumulative condition.
              if (pairLossComparison.active && Index == IndexK0) {
                fillCorrelationHistogram(histos.get<THn>(HIST("PairLossK0/Comparison/Rec")), binFillThn, etaWeight, efficiency * efficiencyTrigg, totalEffUncert, purityTrigg, purityTriggErr);
                // Rec is the complete peak-window sample. Split it into true
                // and fake K0 components without changing any reconstructed
                // selection, coordinate or weight. This makes the MC-truth
                // background removal independently observable.
                if (assocCandidate.mcTrue(IndexK0)) {
                  fillCorrelationHistogram(histos.get<THn>(HIST("PairLossK0/Comparison/RecPeakTrueK0")), binFillThn, etaWeight, efficiency * efficiencyTrigg, totalEffUncert, purityTrigg, purityTriggErr);
                } else {
                  fillCorrelationHistogram(histos.get<THn>(HIST("PairLossK0/Comparison/RecPeakFakeK0")), binFillThn, etaWeight, efficiency * efficiencyTrigg, totalEffUncert, purityTrigg, purityTriggErr);
                }
                const double pairWeight = etaWeight * purityTrigg / (efficiency * efficiencyTrigg);
                auto fillWeightedRecVariant = [&](int variant, double fillDeltaPhi, double fillDeltaEta, double fillPtAssoc, double fillPtTrigger, double fillVtxZ, double fillMult) {
                  histos.fill(HIST("PairLossK0/Comparison/RecVariants"), variant, fillDeltaPhi, fillDeltaEta, fillPtAssoc, fillPtTrigger, fillVtxZ, fillMult, pairWeight);
                };

                // 0: the exact Rec entry, including reconstructed coordinates
                // and the ordinary Rec weight.
                fillWeightedRecVariant(0, deltaphi, deltaeta, ptassoc, pttrigger, pvz, mult);

                const auto triggerMcIdIt = pairLossComparison.trackToMc.find(trigg.globalIndex());
                if (triggerMcIdIt != pairLossComparison.trackToMc.end()) {
                  const auto truthTriggerIt = pairLossComparison.truthTriggers.find(triggerMcIdIt->second);
                  if (truthTriggerIt != pairLossComparison.truthTriggers.end()) {
                    auto const& truthTrigger = truthTriggerIt->second;

                    // 1: add only a valid PairLoss truth-trigger match. The
                    // truth-trigger map already contains its species, charge,
                    // eta and pT acceptance. Primary and truth-leading choices
                    // are deliberately kept for later, separate stages.
                    fillWeightedRecVariant(1, deltaphi, deltaeta, ptassoc, pttrigger, pvz, mult);

                    if (assocCandidate.mcTrue(IndexK0)) {
                      // 2: remove only the fake/combinatorial peak component.
                      // No truth-K0 eta or pT acceptance is imposed yet.
                      fillWeightedRecVariant(2, deltaphi, deltaeta, ptassoc, pttrigger, pvz, mult);

                      const auto k0McIdIt = pairLossComparison.v0ToMc.find(assoc.globalIndex());
                      if (k0McIdIt != pairLossComparison.v0ToMc.end()) {
                        const auto truthK0It = pairLossComparison.truthK0s.find(k0McIdIt->second);
                        if (truthK0It != pairLossComparison.truthK0s.end()) {
                          auto const& truthK0 = truthK0It->second;

                          // 3: additionally require the matched truth K0 to
                          // satisfy the old PairLoss truth acceptance.
                          fillWeightedRecVariant(3, deltaphi, deltaeta, ptassoc, pttrigger, pvz, mult);

                          const bool triggerIsK0Daughter = truthTrigger.globalIndex == truthK0.positiveDaughter.globalIndex || truthTrigger.globalIndex == truthK0.negativeDaughter.globalIndex;
                          if (!triggerIsK0Daughter) {
                            // 4: additionally apply the PairLoss truth-level
                            // trigger-versus-K0-daughter exclusion.
                            fillWeightedRecVariant(4, deltaphi, deltaeta, ptassoc, pttrigger, pvz, mult);

                            if (pairLossComparison.primaryTruthTriggerIds.contains(truthTrigger.globalIndex)) {
                              // 5: additionally impose only the configured truth-
                              // trigger primary requirement.
                              fillWeightedRecVariant(5, deltaphi, deltaeta, ptassoc, pttrigger, pvz, mult);

                              if (pairLossComparison.primaryTruthK0Ids.contains(truthK0.globalIndex)) {
                                // 6: additionally impose only the configured truth-
                                // K0 primary requirement.
                                fillWeightedRecVariant(6, deltaphi, deltaeta, ptassoc, pttrigger, pvz, mult);

                                if (pairLossComparison.oldTruthTriggerIds.contains(truthTrigger.globalIndex) &&
                                    pairLossComparison.oldTruthK0Ids.contains(truthK0.globalIndex)) {
                                  // 7: additionally reproduce the old PairLoss
                                  // truth-leading choice, made after primary cuts.
                                  fillWeightedRecVariant(7, deltaphi, deltaeta, ptassoc, pttrigger, pvz, mult);

                                  if (static_cast<int64_t>(assocCandidate.collisionId()) == pairLossComparison.bestCollisionId) {
                                    // 8: additionally retain only pairs from the
                                    // largest-numContrib reconstructed collision.
                                    fillWeightedRecVariant(8, deltaphi, deltaeta, ptassoc, pttrigger, pvz, mult);

                                    // 9 changes ONLY the two pT coordinates. Its
                                    // ratio to stage 8 isolates pT-bin migration:
                                    // delta eta/phi and the event coordinates are
                                    // deliberately left reconstructed.
                                    fillWeightedRecVariant(9, deltaphi, deltaeta, truthK0.pt, truthTrigger.pt, pvz, mult);

                                    // 10 additionally moves the event coordinates
                                    // to the best reconstructed collision, so that
                                    // vertex-z and multiplicity binning is separable
                                    // from pT migration.
                                    fillWeightedRecVariant(10, deltaphi, deltaeta, truthK0.pt, truthTrigger.pt, pairLossComparison.vtxZ, pairLossComparison.multiplicity);

                                    float truthDeltaEta = truthTrigger.eta - truthK0.eta;
                                    if (masterConfigurations.doMirroringInDelataEta) {
                                      truthDeltaEta = std::abs(truthDeltaEta);
                                    }
                                    const float truthDeltaPhi = computeDeltaPhi(truthTrigger.phi, truthK0.phi);
                                    const bool truthDeltaEtaInRange = truthDeltaEta >= axisRanges[1][0] && truthDeltaEta <= axisRanges[1][1];
                                    const bool truthDeltaPhiInRange = truthDeltaPhi >= axisRanges[0][0] && truthDeltaPhi <= axisRanges[0][1];
                                    if (truthDeltaEtaInRange) {
                                      // 11 changes only delta eta and applies its
                                      // truth range. Delta phi stays reconstructed,
                                      // so a step here is delta-eta resolution alone.
                                      fillWeightedRecVariant(11, deltaphi, truthDeltaEta, truthK0.pt, truthTrigger.pt, pairLossComparison.vtxZ, pairLossComparison.multiplicity);

                                      if (truthDeltaPhiInRange) {
                                        // 12 additionally changes delta phi and
                                        // completes the truth angular range.
                                        fillWeightedRecVariant(12, truthDeltaPhi, truthDeltaEta, truthK0.pt, truthTrigger.pt, pairLossComparison.vtxZ, pairLossComparison.multiplicity);

                                        // 13 removes only the ordinary Rec pair
                                        // weight. Reconstructed matches remain separate.
                                        histos.fill(HIST("PairLossK0/Comparison/RecVariants"), 13, truthDeltaPhi, truthDeltaEta, truthK0.pt, truthTrigger.pt, pairLossComparison.vtxZ, pairLossComparison.multiplicity);

                                        const PairLossPairKey truthPairKey{truthTrigger.globalIndex, truthK0.globalIndex};
                                        if (pairLossComparison.countedPairs.insert(truthPairKey).second) {
                                          // 14 additionally keeps one reconstructed
                                          // match per truth pair. Stage 15 and Final
                                          // are filled by the literal old finalPair.
                                          histos.fill(HIST("PairLossK0/Comparison/RecVariants"), 14, truthDeltaPhi, truthDeltaEta, truthK0.pt, truthTrigger.pt, pairLossComparison.vtxZ, pairLossComparison.multiplicity);
                                          if (!pairLossComparison.finalTruthPairs.contains(truthPairKey)) {
                                            histos.fill(HIST("PairLossK0/Comparison/RecNotInFinal"), truthDeltaPhi, truthDeltaEta, truthK0.pt, truthTrigger.pt, pairLossComparison.vtxZ, pairLossComparison.multiplicity);
                                          }
                                        }
                                      }
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
              if (std::abs(deltaphi) < checks.towardDeltaEtaRange && doITSClustersQA) {
                histos.fill(HIST("hITSClusters") + HIST(V0names[Index]) + HIST("NegativeDaughterToward"), ptassoc, negtrack.itsNCls(), assoc.v0radius());
                histos.fill(HIST("hITSClusters") + HIST(V0names[Index]) + HIST("PositiveDaughterToward"), ptassoc, postrack.itsNCls(), assoc.v0radius());
              }
              if (std::abs(deltaphi) > checks.transwerseDeltaEtaRangeMin && std::abs(deltaphi) < checks.transwerseDeltaEtaRangeMax && doITSClustersQA) {
                histos.fill(HIST("hITSClusters") + HIST(V0names[Index]) + HIST("NegativeDaughterTransverse"), ptassoc, negtrack.itsNCls(), assoc.v0radius());
                histos.fill(HIST("hITSClusters") + HIST(V0names[Index]) + HIST("PositiveDaughterTransverse"), ptassoc, postrack.itsNCls(), assoc.v0radius());
              }
              if (doDeltaPhiStarCheck) {
                double deltaPhiStar = calculateAverageDeltaPhiStar(triggForDeltaPhiStar, assocForDeltaPhiStar, bField);
                double deltaPhiStarPion = calculateAverageDeltaPhiStar(triggForDeltaPhiStar, assocForDeltaPhiStarPion, bField);
                if ((Index == IndexK0 && triggSign > Neutral) || (Index == IndexLambda && triggSign > Neutral) || (Index == IndexAntiLambda && triggSign < Neutral)) {
                  histos.fill(HIST("sameEvent/Signal/") + HIST(V0names[Index]) + HIST("DeltaPhiStar"), deltaPhiStar, trigg.eta() - etaProton, 0.5);
                  if (Index == IndexK0) {
                    histos.fill(HIST("sameEvent/Signal/") + HIST(V0names[Index]) + HIST("DeltaPhiStar"), deltaPhiStarPion, trigg.eta() - etaPion, -0.5);
                  }
                } else {
                  histos.fill(HIST("sameEvent/Signal/") + HIST(V0names[Index]) + HIST("DeltaPhiStar"), deltaPhiStar, trigg.eta() - etaProton, -0.5);
                  if (Index == IndexK0) {
                    histos.fill(HIST("sameEvent/Signal/") + HIST(V0names[Index]) + HIST("DeltaPhiStar"), deltaPhiStarPion, trigg.eta() - etaPion, 0.5);
                  }
                }
              }
            }
            if (assocCandidate.compatible(Index, trackSelection.dEdxCompatibility) && (!masterConfigurations.doMCassociation || assocCandidate.mcTrue(Index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && !mixing && +massWindowConfigurations.minBgNSigma < assocCandidate.invMassNSigma(Index) && assocCandidate.invMassNSigma(Index) < +massWindowConfigurations.maxBgNSigma) {
              fillCorrelationHistogram(histos.get<THn>(HIST("sameEvent/RightBg/") + HIST(V0names[Index])), binFillThn, etaWeight, efficiency * efficiencyTrigg, totalEffUncert, purityTrigg, purityTriggErr);
              if (pairLossComparison.active && Index == IndexK0) {
                // Exact Rec right sideband, with the same reconstructed
                // selections, coordinates and weight as the ordinary Rec path.
                fillCorrelationHistogram(histos.get<THn>(HIST("PairLossK0/Comparison/RecRightBg")), binFillThn, etaWeight, efficiency * efficiencyTrigg, totalEffUncert, purityTrigg, purityTriggErr);
              }
              if (doDeltaPhiStarCheck) {
                double deltaPhiStar = calculateAverageDeltaPhiStar(triggForDeltaPhiStar, assocForDeltaPhiStar, bField);
                double deltaPhiStarPion = calculateAverageDeltaPhiStar(triggForDeltaPhiStar, assocForDeltaPhiStarPion, bField);
                if ((Index == IndexK0 && triggSign > Neutral) || (Index == IndexLambda && triggSign > Neutral) || (Index == IndexAntiLambda && triggSign < Neutral)) {
                  histos.fill(HIST("sameEvent/RightBg/") + HIST(V0names[Index]) + HIST("DeltaPhiStar"), deltaPhiStar, trigg.eta() - etaProton, 0.5);
                  if (Index == IndexK0) {
                    histos.fill(HIST("sameEvent/RightBg/") + HIST(V0names[Index]) + HIST("DeltaPhiStar"), deltaPhiStarPion, trigg.eta() - etaPion, -0.5);
                  }
                } else {
                  histos.fill(HIST("sameEvent/RightBg/") + HIST(V0names[Index]) + HIST("DeltaPhiStar"), deltaPhiStar, trigg.eta() - etaProton, -0.5);
                  if (Index == IndexK0) {
                    histos.fill(HIST("sameEvent/RightBg/") + HIST(V0names[Index]) + HIST("DeltaPhiStar"), deltaPhiStarPion, trigg.eta() - etaPion, 0.5);
                  }
                }
              }
            }
            if (assocCandidate.compatible(Index, trackSelection.dEdxCompatibility) && (!masterConfigurations.doMCassociation || assocCandidate.mcTrue(Index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && mixing && -massWindowConfigurations.maxBgNSigma < assocCandidate.invMassNSigma(Index) && assocCandidate.invMassNSigma(Index) < -massWindowConfigurations.minBgNSigma) {
              if (mixingInBf) {
                currentCollision.addValidParticle(assoc.eta(), assoc.phi(), assoc.pt(), 0, efficiency, efficiencyError, Index);
              } else {
                fillCorrelationHistogram(histos.get<THn>(HIST("mixedEvent/LeftBg/") + HIST(V0names[Index])), binFillThn, 1, efficiency * efficiencyTrigg, totalEffUncert, purityTrigg, purityTriggErr);
              }
            }
            if (assocCandidate.compatible(Index, trackSelection.dEdxCompatibility) && (!masterConfigurations.doMCassociation || assocCandidate.mcTrue(Index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && mixing && -massWindowConfigurations.maxPeakNSigma < assocCandidate.invMassNSigma(Index) && assocCandidate.invMassNSigma(Index) < +massWindowConfigurations.maxPeakNSigma) {
              if (mixingInBf) {
                currentCollision.addValidParticle(assoc.eta(), assoc.phi(), assoc.pt(), 1, efficiency, efficiencyError, Index);
              } else {
                fillCorrelationHistogram(histos.get<THn>(HIST("mixedEvent/Signal/") + HIST(V0names[Index])), binFillThn, 1, efficiency * efficiencyTrigg, totalEffUncert, purityTrigg, purityTriggErr);
              }
              if (masterConfigurations.doCorrelationsHadronV0daughter) {
                fillCorrelationHistogram(histos.get<THn>(HIST("mixedEvent/Signal/") + HIST(V0names[Index]) + HIST("_hSameSign")), binFillThnSameSignDaugher, 1, 1, 1, 1, 1);
                fillCorrelationHistogram(histos.get<THn>(HIST("mixedEvent/Signal/") + HIST(V0names[Index]) + HIST("_hOppositeSign")), binFillThnOppositeSignDaugher, 1, 1, 1, 1, 1);
              }
            }
            if (assocCandidate.compatible(Index, trackSelection.dEdxCompatibility) && (!masterConfigurations.doMCassociation || assocCandidate.mcTrue(Index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && mixing && +massWindowConfigurations.minBgNSigma < assocCandidate.invMassNSigma(Index) && assocCandidate.invMassNSigma(Index) < +massWindowConfigurations.maxBgNSigma) {
              if (mixingInBf) {
                currentCollision.addValidParticle(assoc.eta(), assoc.phi(), assoc.pt(), 2, efficiency, efficiencyError, Index);
              } else {
                fillCorrelationHistogram(histos.get<THn>(HIST("mixedEvent/RightBg/") + HIST(V0names[Index])), binFillThn, 1, efficiency * efficiencyTrigg, totalEffUncert, purityTrigg, purityTriggErr);
              }
            }
          }
        });
      }
    }
    if (!mixingInBf || binVtxZ < 0 || binVtxZ > nBinsVtxZ - 1 || binMult < 0 || binMult > nBinsMult - 1) {
      return;
    }
    int binnumb = binMult * nBinsVtxZ + binVtxZ;
    int hastirgorassoc = masterConfigurations.collisionHasTriggOrAssoc;
    if ((hastirgorassoc == 1 && currentCollision.trigParticles.empty()) ||
        (hastirgorassoc == 2 && currentCollision.assocParticles.empty()) ||
        (hastirgorassoc == 3 && currentCollision.trigParticles.empty() && currentCollision.assocParticles.empty()) ||
        (hastirgorassoc == 4 && (currentCollision.trigParticles.empty() || currentCollision.assocParticles.empty()))) {
      return;
    }
    for (const auto& collision : validCollisions[binnumb]) {
      BinningTypePP colBinning{{axesConfigurations.axisVtxZ, axesConfigurations.axisMult}, true};
      // When 'collisionHasTriggOrAssoc' = 0:
      // binContent(hMECollisionBins) = Σ(Submitted jobs(done))[binContent(the same bin of hSECollisionBins) * masterConfigurations.mixingParameter - Σ(k=0 to min(masterConfigurations.mixingParameter,binContent)) k]
      // When 'collisionHasTriggOrAssoc' = 3
      // More collision loss at higher peripheral centrality (fewer target particles); higher bincontent of HIST("mixedEvent/Signal/") + HIST(Cascadenames[Index])
      // due to avoiding vector occupancy by collisions with no trigger&associated particles
      histos.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision.pvz, collision.mult}));
      for (const auto& trigger : collision.trigParticles) {
        for (const auto& assoc : currentCollision.assocParticles) {
          float deltaeta = trigger.eta - assoc.eta;
          float deltaphi = computeDeltaPhi(trigger.phi, assoc.phi);
          float efficiencyTrigg = trigger.efficiency;
          float efficiencyAssoc = assoc.efficiency;
          float efficiencyTriggError = trigger.efficiencyError;
          float efficiencyAssocError = assoc.efficiencyError;
          float totalEffUncert = 0.0;
          float ptassoc = assoc.pt;
          float pttrigger = trigger.pt;
          if (efficiencyFlags.applyEfficiencyPropagation) {
            totalEffUncert = std::sqrt(std::pow(efficiencyTrigg * efficiencyAssocError, 2) + std::pow(efficiencyTriggError * efficiencyAssoc, 2));
          }
          if (masterConfigurations.doMirroringInDelataEta) {
            deltaeta = std::abs(deltaeta);
          }
          std::array<double, 6> binFillThn = {deltaphi, deltaeta, ptassoc, pttrigger, pvz, mult};
          static_for<0, 2>([&](auto i) {
            constexpr int Index = i.value;
            if (Index == assoc.type && assoc.region == 0) {
              fillCorrelationHistogram(histos.get<THn>(HIST("mixedEvent/LeftBg/") + HIST(V0names[Index])), binFillThn, 1, efficiencyTrigg * efficiencyAssoc, totalEffUncert, 1., 0.);
            }
            if (Index == assoc.type && assoc.region == 1) {
              fillCorrelationHistogram(histos.get<THn>(HIST("mixedEvent/Signal/") + HIST(V0names[Index])), binFillThn, 1, efficiencyTrigg * efficiencyAssoc, totalEffUncert, 1., 0.);
            }
            if (Index == assoc.type && assoc.region == 2) {
              fillCorrelationHistogram(histos.get<THn>(HIST("mixedEvent/RightBg/") + HIST(V0names[Index])), binFillThn, 1, efficiencyTrigg * efficiencyAssoc, totalEffUncert, 1., 0.);
            }
          });
        }
      }
    }
    if (validCollisions[binnumb].size() >= static_cast<size_t>(masterConfigurations.mixingParameter)) {
      validCollisions[binnumb].erase(validCollisions[binnumb].begin());
    }
    if (!currentCollision.trigParticles.empty()) {
      validCollisions[binnumb].push_back(currentCollision);
    }
  }

  void fillCorrelationsCascade(aod::TriggerTracks const& triggers, aod::AssocCascades const& assocs, bool mixing, bool mixingInBf, float pvx, float pvy, float pvz, float mult, double bField)
  {
    ValidCollision currentCollision;
    int nBinsMult = 0;
    int binMult = 0;
    int nBinsVtxZ = 0;
    int binVtxZ = 0;
    currentCollision.pvz = pvz;
    currentCollision.mult = mult;
    if (mixingInBf) {
      nBinsMult = histos.get<TH1>(HIST("axes/hMultAxis"))->GetNbinsX();
      binMult = histos.get<TH1>(HIST("axes/hMultAxis"))->GetXaxis()->FindBin(mult) - 1;
      nBinsVtxZ = histos.get<TH1>(HIST("axes/hVertexZAxis"))->GetNbinsX();
      binVtxZ = histos.get<TH1>(HIST("axes/hVertexZAxis"))->GetXaxis()->FindBin(pvz) - 1;
    }
    bool firstLoop = false;
    for (auto const& triggerTrack : triggers) {
      if (masterConfigurations.doTriggPhysicalPrimary && !triggerTrack.mcPhysicalPrimary()) {
        continue;
      }
      auto trigg = triggerTrack.track_as<TracksComplete>();
      if (!isValidTrigger(trigg, triggerTrack.isLeading())) {
        continue;
      }
      if (trackSelection.checkForITSTPCMissmatchMC && bitcheck(triggerTrack.mcMask(), 13)) {
        continue;
      }

      float efficiencyTrigg = 1.0f;
      float efficiencyTriggError = 0.0f;
      float purityTrigg = 1.0f;
      float purityTriggErr = 0.0f;
      if (efficiencyFlags.applyEfficiencyForTrigger) {
        std::array<double, 4> bintrig = {trigg.pt(), trigg.eta(), trigg.phi(), mult};
        if (efficiencyFlags.applyEffAsFunctionOfMult) {
          efficiencyTrigg = hEfficiencyTriggerMult->Interpolate(trigg.pt(), trigg.eta(), mult);
        } else if (efficiencyFlags.applyEffAsFunctionOfMultAndPhi) {
          efficiencyTrigg = hEfficiencyTriggerMultVsPhi->GetBinContent(hEfficiencyTriggerMultVsPhi->GetBin(bintrig.data()));
        } else {
          efficiencyTrigg = hEfficiencyTrigger->Interpolate(trigg.pt(), trigg.eta());
        }
        if (efficiencyFlags.applyPurityTrigger) {
          if (efficiencyFlags.applyEffAsFunctionOfMult) {
            purityTrigg = hPurityHadronMult->Interpolate(trigg.pt(), mult);
          } else {
            purityTrigg = hPurityHadron->Interpolate(trigg.pt());
          }
        }
        if (efficiencyFlags.applyEfficiencyPropagation) {
          if (efficiencyFlags.applyEffAsFunctionOfMult) {
            efficiencyTriggError = hEfficiencyUncertaintyTriggerMult->Interpolate(trigg.pt(), trigg.eta(), mult);
          } else if (efficiencyFlags.applyEffAsFunctionOfMultAndPhi) {
            efficiencyTriggError = hEfficiencyTriggerMultVsPhi->GetBinError(hEfficiencyTriggerMultVsPhi->GetBin(bintrig.data()));
          } else {
            efficiencyTriggError = hEfficiencyUncertaintyTrigger->Interpolate(trigg.pt(), trigg.eta());
          }
          if (efficiencyFlags.applyPurityTrigger) {
            if (efficiencyFlags.applyEffAsFunctionOfMult) {
              purityTriggErr = hPurityUncertaintyHadronMult->Interpolate(trigg.pt(), mult);
            } else {
              purityTriggErr = hPurityUncertaintyHadron->Interpolate(trigg.pt());
            }
          }
        }
        if (efficiencyTrigg == 0) { // check for zero efficiency, do not apply if the case
          efficiencyTrigg = 1;
          efficiencyTriggError = 0;
        }
      }
      if (!mixing) {
        fillTriggerHistogram(histos.get<TH2>(HIST("sameEvent/TriggerParticlesCascade")), trigg.pt(), mult, efficiencyTrigg, efficiencyTriggError, purityTrigg, purityTriggErr);
      }
      double triggSign = trigg.sign();
      std::array<double, 3> triggForDeltaPhiStar = {trigg.phi(), trigg.pt(), triggSign};

      if (mixingInBf) {
        currentCollision.addValidParticle(trigg.eta(), trigg.phi(), trigg.pt(), -1, efficiencyTrigg, efficiencyTriggError, -1);
        if (firstLoop) {
          continue;
        }
      }
      for (auto const& assocCandidate : assocs) {
        firstLoop = true;
        auto assoc = assocCandidate.cascData();

        //---] syst cuts [---
        if (masterConfigurations.doPPAnalysis && (std::abs(assoc.dcapostopv()) < v0Selection.dcapostopv ||
                                                  std::abs(assoc.dcanegtopv()) < v0Selection.dcanegtopv ||
                                                  std::abs(assoc.dcabachtopv()) < cascadeSelections.dcaBachToPV ||
                                                  assoc.dcaV0daughters() > cascadeSelections.cascdcaV0dau ||
                                                  assoc.dcacascdaughters() > cascadeSelections.dcaCascDaughters ||
                                                  assoc.v0cosPA(pvx, pvy, pvz) < cascadeSelections.cascv0cospa ||
                                                  assoc.casccosPA(pvx, pvy, pvz) < cascadeSelections.cascCospa ||
                                                  assoc.cascradius() < cascadeSelections.cascRadius ||
                                                  std::abs(assoc.dcav0topv(pvx, pvy, pvz)) < cascadeSelections.cascdcaV0ToPV ||
                                                  std::abs(assoc.mLambda() - o2::constants::physics::MassLambda0) > cascadeSelections.cascV0masswindow)) {
          continue;
        }
        if (!masterConfigurations.doPPAnalysis && !cascadeSelectedPbPb(assoc, pvx, pvy, pvz)) {
          continue;
        }
        uint64_t cascselMap = cascadeselectionBitmap(assoc, pvx, pvy, pvz);
        //---] removing autocorrelations [---
        auto postrack = assoc.posTrack_as<TracksComplete>();
        auto negtrack = assoc.negTrack_as<TracksComplete>();
        auto bachtrack = assoc.bachelor_as<TracksComplete>();
        if (doAutocorrelationRejection) {
          if (trigg.globalIndex() == postrack.globalIndex()) {
            histos.fill(HIST("hNumberOfRejectedPairsCascades"), 0.5);
            continue;
          }
          if (trigg.globalIndex() == negtrack.globalIndex()) {
            histos.fill(HIST("hNumberOfRejectedPairsCascades"), 0.5);
            continue;
          }
          if (trigg.globalIndex() == bachtrack.globalIndex()) {
            histos.fill(HIST("hNumberOfRejectedPairsCascades"), 0.5);
            continue;
          }
        }

        double phiProton = postrack.phi();
        double etaProton = postrack.eta();
        double ptProton = postrack.pt();
        double signProton = postrack.sign();
        if (assoc.sign() > 0) {
          phiProton = negtrack.phi();
          etaProton = negtrack.eta();
          ptProton = negtrack.pt();
          signProton = negtrack.sign();
        }
        std::array<double, 3> assocForDeltaPhiStar = {phiProton, ptProton, signProton};
        //---] track quality check [---
        if (postrack.tpcNClsCrossedRows() < trackSelection.minTPCNCrossedRowsAssociated || negtrack.tpcNClsCrossedRows() < trackSelection.minTPCNCrossedRowsAssociated || bachtrack.tpcNClsCrossedRows() < trackSelection.minTPCNCrossedRowsAssociated) {
          continue;
        }
        if (trackSelection.checksRequireTPCChi2 && (postrack.tpcChi2NCl() < trackSelection.minTPCChi2PerClusterAssociated || negtrack.tpcChi2NCl() < trackSelection.minTPCChi2PerClusterAssociated || bachtrack.tpcChi2NCl() < trackSelection.minTPCChi2PerClusterAssociated)) {
          continue;
        }
        if (trackSelection.requireClusterInITS && (postrack.itsNCls() < trackSelection.minITSClustersForDaughterTracks || negtrack.itsNCls() < trackSelection.minITSClustersForDaughterTracks || bachtrack.itsNCls() < trackSelection.minITSClustersForDaughterTracks)) {
          continue;
        }

        float deltaphi = computeDeltaPhi(trigg.phi(), assoc.phi());
        float deltaeta = trigg.eta() - assoc.eta();
        if (masterConfigurations.doMirroringInDelataEta) {
          deltaeta = std::abs(deltaeta);
        }
        float ptassoc = assoc.pt();
        float pttrigger = trigg.pt();

        // skip if basic ranges not met
        if (deltaphi < axisRanges[0][0] || deltaphi > axisRanges[0][1]) {
          continue;
        }
        if (deltaeta < axisRanges[1][0] || deltaeta > axisRanges[1][1]) {
          continue;
        }
        if (ptassoc < axisRanges[2][0] || ptassoc > axisRanges[2][1]) {
          continue;
        }

        std::array<TH2F*, AssocCascadeTypes> hEfficiencyCascade{nullptr, nullptr, nullptr, nullptr};
        std::array<THnF*, AssocCascadeTypes> hEfficiencyCascadeMultVsPhi{nullptr, nullptr, nullptr, nullptr};
        if (efficiencyFlags.applyEffAsFunctionOfMultAndPhi) {
          hEfficiencyCascadeMultVsPhi[0] = hEfficiencyXiMinusMultVsPhi;
          hEfficiencyCascadeMultVsPhi[1] = hEfficiencyXiPlusMultVsPhi;
          hEfficiencyCascadeMultVsPhi[2] = hEfficiencyOmegaMinusMultVsPhi;
          hEfficiencyCascadeMultVsPhi[3] = hEfficiencyOmegaPlusMultVsPhi;
        } else {
          hEfficiencyCascade[0] = hEfficiencyXiMinus;
          hEfficiencyCascade[1] = hEfficiencyXiPlus;
          hEfficiencyCascade[2] = hEfficiencyOmegaMinus;
          hEfficiencyCascade[3] = hEfficiencyOmegaPlus;
        }

        std::array<TH2F*, AssocCascadeTypes> hEfficiencyUncertaintyCascade{nullptr, nullptr, nullptr, nullptr};
        hEfficiencyUncertaintyCascade[0] = hEfficiencyUncertaintyXiMinus;
        hEfficiencyUncertaintyCascade[1] = hEfficiencyUncertaintyXiPlus;
        hEfficiencyUncertaintyCascade[2] = hEfficiencyUncertaintyOmegaMinus;
        hEfficiencyUncertaintyCascade[3] = hEfficiencyUncertaintyOmegaPlus;

        float etaWeight = 1;
        if (checks.doOnTheFlyFlattening) {
          float preWeight = 1 - std::abs(deltaeta) / 1.6;
          etaWeight = preWeight != 0 ? 1.0f / preWeight : 1.0f;
        }

        static_for<0, 3>([&](auto i) {
          constexpr int Index = i.value;
          if ((Index == IndexOmegaMinus || Index == IndexOmegaPlus) && assocCandidate.compatible(Index, trackSelection.dEdxCompatibility) && std::abs(assocCandidate.invMassNSigma(Index - 2)) < massWindowConfigurations.nSigmaNearXiMassCenter) {
            return;
          }
          float efficiency = 1.0f;
          float totalEffUncert = 0.0;
          float efficiencyError = 0.0f;
          if (efficiencyFlags.applyEfficiencyCorrection) {
            if (efficiencyFlags.applyEffAsFunctionOfMultAndPhi) {
              std::array<double, 4> bin = {ptassoc, assoc.eta(), assoc.phi(), mult};
              efficiency = hEfficiencyCascadeMultVsPhi[Index]->GetBinContent(hEfficiencyCascadeMultVsPhi[Index]->GetBin(bin.data()));
              if (efficiencyFlags.applyEfficiencyPropagation) {
                efficiencyError = hEfficiencyCascadeMultVsPhi[Index]->GetBinError(hEfficiencyCascadeMultVsPhi[Index]->GetBin(bin.data()));
              }
            } else {
              efficiency = hEfficiencyCascade[Index]->Interpolate(ptassoc, assoc.eta());
              if (efficiencyFlags.applyEfficiencyPropagation) {
                efficiencyError = hEfficiencyUncertaintyCascade[Index]->Interpolate(ptassoc, assoc.eta());
              }
            }
          }
          if (efficiency == 0) { // check for zero efficiency, do not apply if the case
            efficiency = 1;
            efficiencyError = 0;
          }
          if (efficiencyFlags.applyEfficiencyPropagation) {
            totalEffUncert = std::sqrt(std::pow(efficiencyTrigg * efficiencyError, 2) + std::pow(efficiencyTriggError * efficiency, 2));
          }

          std::array<double, 6> binFillThn = {deltaphi, deltaeta, ptassoc, pttrigger, pvz, mult};
          if (TESTBIT(doCorrelation, Index + 3) && (!efficiencyFlags.applyEfficiencyCorrection || efficiency != 0) && (masterConfigurations.doPPAnalysis || (TESTBIT(cascselMap, Index) && TESTBIT(cascselMap, Index + 4) && TESTBIT(cascselMap, Index + 8) && TESTBIT(cascselMap, Index + 12)))) {
            if (assocCandidate.compatible(Index, trackSelection.dEdxCompatibility) && (!masterConfigurations.doMCassociation || assocCandidate.mcTrue(Index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && !mixing && -massWindowConfigurations.maxBgNSigma < assocCandidate.invMassNSigma(Index) && assocCandidate.invMassNSigma(Index) < -massWindowConfigurations.minBgNSigma) {
              fillCorrelationHistogram(histos.get<THn>(HIST("sameEvent/LeftBg/") + HIST(Cascadenames[Index])), binFillThn, etaWeight, efficiency * efficiencyTrigg, totalEffUncert, purityTrigg, purityTriggErr);
              if (doDeltaPhiStarCheck) {
                double deltaPhiStar = calculateAverageDeltaPhiStar(triggForDeltaPhiStar, assocForDeltaPhiStar, bField);
                if ((Index == IndexXiMinus && triggSign > Neutral) || (Index == IndexXiPlus && triggSign < Neutral) || (Index == IndexOmegaMinus && triggSign > Neutral) || (Index == IndexOmegaPlus && triggSign < 0)) {
                  histos.fill(HIST("sameEvent/LeftBg/") + HIST(Cascadenames[Index]) + HIST("DeltaPhiStar"), deltaPhiStar, trigg.eta() - etaProton, 0.5);
                } else {
                  histos.fill(HIST("sameEvent/LeftBg/") + HIST(Cascadenames[Index]) + HIST("DeltaPhiStar"), deltaPhiStar, trigg.eta() - etaProton, -0.5);
                }
              }
            }
            if (assocCandidate.compatible(Index, trackSelection.dEdxCompatibility) && (!masterConfigurations.doMCassociation || assocCandidate.mcTrue(Index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && !mixing && -massWindowConfigurations.maxPeakNSigma < assocCandidate.invMassNSigma(Index) && assocCandidate.invMassNSigma(Index) < +massWindowConfigurations.maxPeakNSigma) {
              fillCorrelationHistogram(histos.get<THn>(HIST("sameEvent/Signal/") + HIST(Cascadenames[Index])), binFillThn, etaWeight, efficiency * efficiencyTrigg, totalEffUncert, purityTrigg, purityTriggErr);
              if (doDeltaPhiStarCheck) {
                double deltaPhiStar = calculateAverageDeltaPhiStar(triggForDeltaPhiStar, assocForDeltaPhiStar, bField);
                if ((Index == IndexXiMinus && triggSign > Neutral) || (Index == IndexXiPlus && triggSign < Neutral) || (Index == IndexOmegaMinus && triggSign > Neutral) || (Index == IndexOmegaPlus && triggSign < 0)) {
                  histos.fill(HIST("sameEvent/Signal/") + HIST(Cascadenames[Index]) + HIST("DeltaPhiStar"), deltaPhiStar, trigg.eta() - etaProton, 0.5);
                } else {
                  histos.fill(HIST("sameEvent/Signal/") + HIST(Cascadenames[Index]) + HIST("DeltaPhiStar"), deltaPhiStar, trigg.eta() - etaProton, -0.5);
                }
              }
            }
            if (assocCandidate.compatible(Index, trackSelection.dEdxCompatibility) && (!masterConfigurations.doMCassociation || assocCandidate.mcTrue(Index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && !mixing && +massWindowConfigurations.minBgNSigma < assocCandidate.invMassNSigma(Index) && assocCandidate.invMassNSigma(Index) < +massWindowConfigurations.maxBgNSigma) {
              fillCorrelationHistogram(histos.get<THn>(HIST("sameEvent/RightBg/") + HIST(Cascadenames[Index])), binFillThn, etaWeight, efficiency * efficiencyTrigg, totalEffUncert, purityTrigg, purityTriggErr);
              if (doDeltaPhiStarCheck) {
                double deltaPhiStar = calculateAverageDeltaPhiStar(triggForDeltaPhiStar, assocForDeltaPhiStar, bField);
                if ((Index == IndexXiMinus && triggSign > Neutral) || (Index == IndexXiPlus && triggSign < Neutral) || (Index == IndexOmegaMinus && triggSign > Neutral) || (Index == IndexOmegaPlus && triggSign < 0)) {
                  histos.fill(HIST("sameEvent/RightBg/") + HIST(Cascadenames[Index]) + HIST("DeltaPhiStar"), deltaPhiStar, trigg.eta() - etaProton, 0.5);
                } else {
                  histos.fill(HIST("sameEvent/RightBg/") + HIST(Cascadenames[Index]) + HIST("DeltaPhiStar"), deltaPhiStar, trigg.eta() - etaProton, -0.5);
                }
              }
            }
            if (assocCandidate.compatible(Index, trackSelection.dEdxCompatibility) && (!masterConfigurations.doMCassociation || assocCandidate.mcTrue(Index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && mixing && -massWindowConfigurations.maxBgNSigma < assocCandidate.invMassNSigma(Index) && assocCandidate.invMassNSigma(Index) < -massWindowConfigurations.minBgNSigma) {
              if (mixingInBf) {
                currentCollision.addValidParticle(assoc.eta(), assoc.phi(), assoc.pt(), 0, efficiency, efficiencyError, Index);
              } else {
                fillCorrelationHistogram(histos.get<THn>(HIST("mixedEvent/LeftBg/") + HIST(Cascadenames[Index])), binFillThn, 1, efficiency * efficiencyTrigg, totalEffUncert, purityTrigg, purityTriggErr);
              }
            }
            if (assocCandidate.compatible(Index, trackSelection.dEdxCompatibility) && (!masterConfigurations.doMCassociation || assocCandidate.mcTrue(Index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && mixing && -massWindowConfigurations.maxPeakNSigma < assocCandidate.invMassNSigma(Index) && assocCandidate.invMassNSigma(Index) < +massWindowConfigurations.maxPeakNSigma) {
              if (mixingInBf) {
                currentCollision.addValidParticle(assoc.eta(), assoc.phi(), assoc.pt(), 1, efficiency, efficiencyError, Index);
              } else {
                fillCorrelationHistogram(histos.get<THn>(HIST("mixedEvent/Signal/") + HIST(Cascadenames[Index])), binFillThn, 1, efficiency * efficiencyTrigg, totalEffUncert, purityTrigg, purityTriggErr);
              }
            }
            if (assocCandidate.compatible(Index, trackSelection.dEdxCompatibility) && (!masterConfigurations.doMCassociation || assocCandidate.mcTrue(Index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && mixing && +massWindowConfigurations.minBgNSigma < assocCandidate.invMassNSigma(Index) && assocCandidate.invMassNSigma(Index) < +massWindowConfigurations.maxBgNSigma) {
              if (mixingInBf) {
                currentCollision.addValidParticle(assoc.eta(), assoc.phi(), assoc.pt(), 2, efficiency, efficiencyError, Index);
              } else {
                fillCorrelationHistogram(histos.get<THn>(HIST("mixedEvent/RightBg/") + HIST(Cascadenames[Index])), binFillThn, 1, efficiency * efficiencyTrigg, totalEffUncert, purityTrigg, purityTriggErr);
              }
            }
          }
        });
      }
    }
    if (!mixingInBf || binVtxZ < 0 || binVtxZ > nBinsVtxZ - 1 || binMult < 0 || binMult > nBinsMult - 1) {
      return;
    }
    int binnumb = binMult * nBinsVtxZ + binVtxZ;
    int hastirgorassoc = masterConfigurations.collisionHasTriggOrAssoc;
    if ((hastirgorassoc == 1 && currentCollision.trigParticles.empty()) ||
        (hastirgorassoc == 2 && currentCollision.assocParticles.empty()) ||
        (hastirgorassoc == 3 && currentCollision.trigParticles.empty() && currentCollision.assocParticles.empty()) ||
        (hastirgorassoc == 4 && (currentCollision.trigParticles.empty() || currentCollision.assocParticles.empty()))) {
      return;
    }
    for (const auto& collision : validCollisions[binnumb]) {
      BinningTypePP colBinning{{axesConfigurations.axisVtxZ, axesConfigurations.axisMult}, true};
      histos.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision.pvz, collision.mult}));
      for (const auto& trigger : collision.trigParticles) {
        for (const auto& assoc : currentCollision.assocParticles) {
          float deltaeta = trigger.eta - assoc.eta;
          float deltaphi = computeDeltaPhi(trigger.phi, assoc.phi);
          float efficiencyTrigg = trigger.efficiency;
          float efficiencyAssoc = assoc.efficiency;
          float efficiencyTriggError = trigger.efficiencyError;
          float efficiencyAssocError = assoc.efficiencyError;
          float totalEffUncert = 0.0;
          float ptassoc = assoc.pt;
          float pttrigger = trigger.pt;
          if (efficiencyFlags.applyEfficiencyPropagation) {
            totalEffUncert = std::sqrt(std::pow(efficiencyTrigg * efficiencyAssocError, 2) + std::pow(efficiencyTriggError * efficiencyAssoc, 2));
          }
          if (masterConfigurations.doMirroringInDelataEta) {
            deltaeta = std::abs(deltaeta);
          }
          std::array<double, 6> binFillThn = {deltaphi, deltaeta, ptassoc, pttrigger, pvz, mult};
          static_for<0, 3>([&](auto i) {
            constexpr int Index = i.value;
            if (Index == assoc.type && assoc.region == 0) {
              fillCorrelationHistogram(histos.get<THn>(HIST("mixedEvent/LeftBg/") + HIST(Cascadenames[Index])), binFillThn, 1, efficiencyTrigg * efficiencyAssoc, totalEffUncert, 1., 0.);
            }
            if (Index == assoc.type && assoc.region == 1) {
              fillCorrelationHistogram(histos.get<THn>(HIST("mixedEvent/Signal/") + HIST(Cascadenames[Index])), binFillThn, 1, efficiencyTrigg * efficiencyAssoc, totalEffUncert, 1., 0.);
            }
            if (Index == assoc.type && assoc.region == 2) {
              fillCorrelationHistogram(histos.get<THn>(HIST("mixedEvent/RightBg/") + HIST(Cascadenames[Index])), binFillThn, 1, efficiencyTrigg * efficiencyAssoc, totalEffUncert, 1., 0.);
            }
          });
        }
      }
    }
    if (validCollisions[binnumb].size() >= static_cast<size_t>(masterConfigurations.mixingParameter)) {
      validCollisions[binnumb].erase(validCollisions[binnumb].begin());
    }
    if (!currentCollision.trigParticles.empty()) {
      validCollisions[binnumb].push_back(currentCollision);
    }
  }
  template <typename TTriggers, typename THadrons>
  void fillCorrelationsHadron(TTriggers const& triggers, THadrons const& assocs, bool mixing, float pvz, float mult, double bField)
  {

    for (auto const& triggerTrack : triggers) {
      if (masterConfigurations.doTriggPhysicalPrimary && !triggerTrack.mcPhysicalPrimary()) {
        continue;
      }
      auto trigg = triggerTrack.template track_as<TracksComplete>();
      if (!isValidTrigger(trigg, triggerTrack.isLeading())) {
        continue;
      }
      if (trackSelection.checkForITSTPCMissmatchMC && bitcheck(triggerTrack.mcMask(), 13)) {
        continue;
      }

      float efficiencyTrigger = 1.0f;
      float efficiencyTriggerError = 0.0f;
      float purityTrigger = 1.0f;
      float purityTriggerError = 0.0f;
      if (efficiencyFlags.applyEfficiencyForTrigger) {
        if (efficiencyFlags.applyEffAsFunctionOfMult) {
          efficiencyTrigger = hEfficiencyTriggerMult->Interpolate(trigg.pt(), trigg.eta(), mult);
        } else {
          efficiencyTrigger = hEfficiencyTrigger->Interpolate(trigg.pt(), trigg.eta());
        }
        if (efficiencyFlags.applyPurityTrigger) {
          if (efficiencyFlags.applyEffAsFunctionOfMult) {
            purityTrigger = hPurityHadronMult->Interpolate(trigg.pt(), mult);
          } else {
            purityTrigger = hPurityHadron->Interpolate(trigg.pt());
          }
        }
        if (efficiencyFlags.applyEfficiencyPropagation) {
          if (efficiencyFlags.applyEffAsFunctionOfMult) {
            efficiencyTriggerError = hEfficiencyUncertaintyTriggerMult->Interpolate(trigg.pt(), trigg.eta(), mult);
          } else {
            efficiencyTriggerError = hEfficiencyUncertaintyTrigger->Interpolate(trigg.pt(), trigg.eta());
          }
          if (efficiencyFlags.applyPurityTrigger) {
            if (efficiencyFlags.applyEffAsFunctionOfMult) {
              purityTriggerError = hPurityUncertaintyHadronMult->Interpolate(trigg.pt(), mult);
            } else {
              purityTriggerError = hPurityUncertaintyHadron->Interpolate(trigg.pt());
            }
          }
        }
        if (efficiencyTrigger == 0) { // check for zero efficiency, do not apply if the case
          efficiencyTrigger = 1;
          efficiencyTriggerError = 0;
        }
      }
      if (!mixing) {
        if constexpr (requires { triggerTrack.extra(); }) {
          fillTriggerHistogram(histos.get<TH2>(HIST("sameEvent/TriggerParticlesPion")), trigg.pt(), mult, efficiencyTrigger, efficiencyTriggerError, purityTrigger, purityTriggerError);
        } else {
          fillTriggerHistogram(histos.get<TH2>(HIST("sameEvent/TriggerParticlesHadron")), trigg.pt(), mult, efficiencyTrigger, efficiencyTriggerError, purityTrigger, purityTriggerError);
        }
      }
      double triggSign = trigg.sign();
      std::array<double, 3> triggForDeltaPhiStar = {trigg.phi(), trigg.pt(), triggSign};
      for (auto const& assocTrack : assocs) {
        auto assoc = assocTrack.template track_as<TracksComplete>();

        //---] removing autocorrelations [---
        if (doAutocorrelationRejection) {
          if (trigg.globalIndex() == assoc.globalIndex()) {
            if constexpr (requires { assocTrack.nSigmaTPCPi(); }) {
              histos.fill(HIST("hNumberOfRejectedPairsPion"), 0.5);
            } else {
              histos.fill(HIST("hNumberOfRejectedPairsHadron"), 0.5);
            }
            continue;
          }
        }
        //---] track quality check [---
        if (!isValidAssocHadron(assoc)) {
          continue;
        }
        if (trackSelection.checkForITSTPCMissmatchMC && bitcheck(assocTrack.mcMask(), 13)) {
          continue;
        }
        if (doAssocPhysicalPrimary && !assocTrack.mcPhysicalPrimary()) {
          continue;
        }
        float deltaphi = computeDeltaPhi(trigg.phi(), assoc.phi());
        float deltaeta = trigg.eta() - assoc.eta();
        if (masterConfigurations.doMirroringInDelataEta) {
          deltaeta = std::abs(deltaeta);
        }
        float ptassoc = assoc.pt();
        float pttrigger = trigg.pt();

        double assocSign = assoc.sign();
        std::array<double, 3> assocForDeltaPhiStar = {assoc.phi(), assoc.pt(), assocSign};

        float etaWeight = 1.;
        if (checks.doOnTheFlyFlattening) {
          float preWeight = 1 - std::abs(deltaeta) / 1.6;
          etaWeight = preWeight != 0 ? 1.0f / preWeight : 1.0f;
        }

        // skip if basic ranges not met
        if (deltaphi < axisRanges[0][0] || deltaphi > axisRanges[0][1]) {
          continue;
        }
        if (deltaeta < axisRanges[1][0] || deltaeta > axisRanges[1][1]) {
          continue;
        }
        if (ptassoc < axisRanges[2][0] || ptassoc > axisRanges[2][1]) {
          continue;
        }

        float efficiency = 1;
        float purity = 1.0f;
        float purityUncertainty = 0.0f;
        float totalEffUncert = 0.0;
        float efficiencyUncertainty = 0.0f;
        float totalPurityUncert = 0.0;
        if (efficiencyFlags.applyEfficiencyCorrection) {
          if constexpr (requires { assocTrack.nSigmaTPCPi(); }) {
            efficiency = hEfficiencyPion->Interpolate(ptassoc, assoc.eta());
            if (efficiencyFlags.applyEfficiencyPropagation) {
              efficiencyUncertainty = hEfficiencyUncertaintyPion->Interpolate(ptassoc, assoc.eta());
            }
          } else {
            if (efficiencyFlags.applyEffAsFunctionOfMult) {
              efficiency = hEfficiencyHadronMult->Interpolate(ptassoc, assoc.eta(), mult);
            } else {
              efficiency = hEfficiencyHadron->Interpolate(ptassoc, assoc.eta());
            }
            if (efficiencyFlags.applyPurityHadron) {
              if (efficiencyFlags.applyEffAsFunctionOfMult) {
                purity = hPurityHadronMult->Interpolate(ptassoc, mult);
              } else {
                purity = hPurityHadron->Interpolate(ptassoc);
              }
            }
            if (efficiencyFlags.applyEfficiencyPropagation) {
              if (efficiencyFlags.applyEffAsFunctionOfMult) {
                efficiencyUncertainty = hEfficiencyUncertaintyHadronMult->Interpolate(ptassoc, assoc.eta(), mult);
              } else {
                efficiencyUncertainty = hEfficiencyUncertaintyHadron->Interpolate(ptassoc, assoc.eta());
              }
              if (efficiencyFlags.applyPurityHadron) {
                if (efficiencyFlags.applyEffAsFunctionOfMult) {
                  purityUncertainty = hPurityUncertaintyHadronMult->Interpolate(ptassoc, mult);
                } else {
                  purityUncertainty = hPurityUncertaintyHadron->Interpolate(ptassoc);
                }
              }
            }
          }
        }
        if (efficiency == 0) { // check for zero efficiency, do not apply if the case
          efficiency = 1;
          efficiencyUncertainty = 0.0;
        }
        if (efficiencyFlags.applyEfficiencyPropagation) {
          totalEffUncert = std::sqrt(std::pow(efficiencyTrigger * efficiencyUncertainty, 2) + std::pow(efficiencyTriggerError * efficiency, 2));
          totalPurityUncert = std::sqrt(std::pow(purityTrigger * purityUncertainty, 2) + std::pow(purity * purityTriggerError, 2));
        }

        std::array<double, 6> binFillThn = {deltaphi, deltaeta, ptassoc, pttrigger, pvz, mult};
        double deltaPhiStar = calculateAverageDeltaPhiStar(triggForDeltaPhiStar, assocForDeltaPhiStar, bField);
        if (!mixing) {
          if constexpr (requires { assocTrack.nSigmaTPCPi(); }) {
            fillCorrelationHistogram(histos.get<THn>(HIST("sameEvent/Signal/Pion")), binFillThn, etaWeight, efficiency * efficiencyTrigger, totalEffUncert, purity * purityTrigger, totalPurityUncert);
            if (triggSign == assocSign && doDeltaPhiStarCheck) {
              histos.fill(HIST("sameEvent/Signal/Pion") + HIST("DeltaPhiStar"), deltaPhiStar, trigg.eta() - assoc.eta(), 0.5);
            } else if (doDeltaPhiStarCheck) {
              histos.fill(HIST("sameEvent/Signal/Pion") + HIST("DeltaPhiStar"), deltaPhiStar, trigg.eta() - assoc.eta(), -0.5);
            }
          } else {
            if (triggSign == assocSign && doDeltaPhiStarCheck) {
              histos.fill(HIST("sameEvent/Signal/Hadron") + HIST("DeltaPhiStar"), deltaPhiStar, trigg.eta() - assoc.eta(), 0.5);
            } else if (doDeltaPhiStarCheck) {
              histos.fill(HIST("sameEvent/Signal/Hadron") + HIST("DeltaPhiStar"), deltaPhiStar, trigg.eta() - assoc.eta(), -0.5);
            }
            fillCorrelationHistogram(histos.get<THn>(HIST("sameEvent/Signal/Hadron")), binFillThn, etaWeight, efficiency * efficiencyTrigger, totalEffUncert, purity * purityTrigger, totalPurityUncert);
          }
        } else {
          if constexpr (requires { assocTrack.nSigmaTPCPi(); }) {
            fillCorrelationHistogram(histos.get<THn>(HIST("mixedEvent/Signal/Pion")), binFillThn, 1, efficiency * efficiencyTrigger, totalEffUncert, purity * purityTrigger, totalPurityUncert);
          } else {
            fillCorrelationHistogram(histos.get<THn>(HIST("mixedEvent/Signal/Hadron")), binFillThn, 1, efficiency * efficiencyTrigger, totalEffUncert, purity * purityTrigger, totalPurityUncert);
          }
        }
      }
    }
  }

  void init(InitContext const&)
  {
    zorroSummary.setObject(zorro.getZorroSummary());
    hEfficiencyPion = nullptr;
    hEfficiencyK0Short = nullptr;
    hEfficiencyLambda = nullptr;
    hEfficiencyAntiLambda = nullptr;
    hEfficiencyXiMinus = nullptr;
    hEfficiencyXiPlus = nullptr;
    hEfficiencyOmegaMinus = nullptr;
    hEfficiencyOmegaPlus = nullptr;
    hEfficiencyUncertaintyTrigger = nullptr;
    hEfficiencyUncertaintyXiMinus = nullptr;
    hEfficiencyUncertaintyXiPlus = nullptr;
    hEfficiencyUncertaintyOmegaMinus = nullptr;
    hEfficiencyUncertaintyOmegaPlus = nullptr;
    hEfficiencyUncertaintyPion = nullptr;
    hEfficiencyUncertaintyK0Short = nullptr;
    hEfficiencyUncertaintyLambda = nullptr;
    hEfficiencyUncertaintyAntiLambda = nullptr;

    hEfficiencyHadron = nullptr;
    hPurityHadron = nullptr;
    hPurityUncertaintyHadron = nullptr;
    hEfficiencyUncertaintyHadron = nullptr;

    // set bitmap for convenience
    doCorrelation = 0;
    if (masterConfigurations.doCorrelationK0Short) {
      SETBIT(doCorrelation, 0);
    }
    if (masterConfigurations.doCorrelationLambda) {
      SETBIT(doCorrelation, 1);
    }
    if (masterConfigurations.doCorrelationAntiLambda) {
      SETBIT(doCorrelation, 2);
    }
    if (masterConfigurations.doCorrelationXiMinus) {
      SETBIT(doCorrelation, 3);
    }
    if (masterConfigurations.doCorrelationXiPlus) {
      SETBIT(doCorrelation, 4);
    }
    if (masterConfigurations.doCorrelationOmegaMinus) {
      SETBIT(doCorrelation, 5);
    }
    if (masterConfigurations.doCorrelationOmegaPlus) {
      SETBIT(doCorrelation, 6);
    }
    if (masterConfigurations.doCorrelationPion) {
      SETBIT(doCorrelation, 7);
    }
    if (masterConfigurations.doCorrelationHadron) {
      SETBIT(doCorrelation, 8);
    }

    // Store axis ranges to prevent spurious filling
    // axis status:
    // --- Delta-phi is safe -> math forbids insanity
    // --- Delta-eta depends on pre-filter -> check
    // --- pT assoc depends on binning -> check
    // --- vertex Z is safe -> skipped at evsel level
    // --- multiplicity -> check

    // grab axis edge from ConfigurableAxes
    const AxisSpec preAxisDeltaPhi{axesConfigurations.axisDeltaPhi, "#Delta#varphi"};
    const AxisSpec preAxisDeltaEta{axesConfigurations.axisDeltaEta, "#Delta#eta"};
    const AxisSpec preAxisPtAssoc{axesConfigurations.axisPtAssoc, "#it{p}_{T}^{assoc} (GeV/c)"};
    const AxisSpec preAxisPtTrigger{axesConfigurations.axisPtTrigger, "#it{p}_{T}^{trigger} (GeV/c)"};
    const AxisSpec preAxisVtxZ{axesConfigurations.axisVtxZ, "vertex Z (cm)"};
    const AxisSpec preAxisMult{axesConfigurations.axisMult, "mult percentile"};
    const AxisSpec axisPtLambda{axesConfigurations.axisPtAssoc, "#it{p}_{T}^{#Lambda} (GeV/c)"};
    const AxisSpec axisPtCascade{axesConfigurations.axisPtAssoc, "#it{p}_{T}^{Mother} (GeV/c)"};
    const AxisSpec preAxisMultiplicity{axesConfigurations.axisMultiplicity, "multiplicity"};

    // store the original axes in specific TH1Cs for completeness
    histos.add("axes/hDeltaPhiAxis", "", kTH1C, {preAxisDeltaPhi});
    histos.add("axes/hDeltaEtaAxis", "", kTH1C, {preAxisDeltaEta});
    histos.add("axes/hPtAssocAxis", "", kTH1C, {preAxisPtAssoc});
    histos.add("axes/hPtTriggerAxis", "", kTH1C, {preAxisPtTrigger});
    histos.add("axes/hVertexZAxis", "", kTH1C, {preAxisVtxZ});
    histos.add("axes/hMultAxis", "", kTH1C, {preAxisMult});
    histos.add("axes/hMultiplicityAxis", "", kTH1C, {preAxisMultiplicity});

    std::vector<double> edgesDeltaPhiOrig = preAxisDeltaPhi.binEdges;
    std::vector<double> edgesDeltaEtaOrig = preAxisDeltaEta.binEdges;
    std::vector<double> edgesPtAssocOrig = preAxisPtAssoc.binEdges;
    std::vector<double> edgesPtTriggerOrig = preAxisPtTrigger.binEdges;
    std::vector<double> edgesVtxZOrig = preAxisVtxZ.binEdges;
    std::vector<double> edgesMultOrig = preAxisMult.binEdges;
    std::vector<double> edgesMultiplicityOrig = preAxisMultiplicity.binEdges;

    std::vector<float> rangesDeltaPhi = {static_cast<float>(edgesDeltaPhiOrig[0]), static_cast<float>(edgesDeltaPhiOrig[edgesDeltaPhiOrig.size() - 1])};
    std::vector<float> rangesDeltaEta = {static_cast<float>(edgesDeltaEtaOrig[0]), static_cast<float>(edgesDeltaEtaOrig[edgesDeltaEtaOrig.size() - 1])};
    std::vector<float> rangesPtAssoc = {static_cast<float>(edgesPtAssocOrig[0]), static_cast<float>(edgesPtAssocOrig[edgesPtAssocOrig.size() - 1])};
    std::vector<float> rangesPtTrigger = {static_cast<float>(edgesPtTriggerOrig[0]), static_cast<float>(edgesPtTriggerOrig[edgesPtTriggerOrig.size() - 1])};
    std::vector<float> rangesVtxZ = {static_cast<float>(edgesVtxZOrig[0]), static_cast<float>(edgesVtxZOrig[edgesVtxZOrig.size() - 1])};
    std::vector<float> rangesMult = {static_cast<float>(edgesMultOrig[0]), static_cast<float>(edgesMultOrig[edgesMultOrig.size() - 1])};
    std::vector<float> rangesMultiplicity = {static_cast<float>(edgesMultiplicityOrig[0]), static_cast<float>(edgesMultiplicityOrig[edgesMultiplicityOrig.size() - 1])};

    axisRanges.emplace_back(rangesDeltaPhi);
    axisRanges.emplace_back(rangesDeltaEta);
    axisRanges.emplace_back(rangesPtAssoc);
    axisRanges.emplace_back(rangesPtTrigger);
    axisRanges.emplace_back(rangesVtxZ);
    axisRanges.emplace_back(rangesMult);
    axisRanges.emplace_back(rangesMultiplicity);

    std::vector<double> edgesDeltaPhi;
    std::vector<double> edgesDeltaEta;
    std::vector<double> edgesPtAssoc;
    std::vector<double> edgesPtTrigger;
    std::vector<double> edgesVtxZ;
    std::vector<double> edgesMult;
    std::vector<double> edgesMultiplicity;

    // v--- skipUnderOverflowInTHn ---v
    //
    // if enabled, this will change the axes such that they will solely cover the interval from
    // edge[1] to edge[n-1]; this will mean that the bin 1 and bin N will be stored in
    // under / overflow bins and will have to be manually unpacked. Do not forget to do the manual
    // unpacking a posteriori!
    //
    // this feature is meant to save memory conveniently.
    // it should actually be implemented centrally in ROOT but ok, this will do it for now.

    int offset = masterConfigurations.skipUnderOverflowInTHn ? 1 : 0;
    // ===] delta-phi [===
    if (!preAxisDeltaPhi.nBins.has_value()) {
      // variable binning, use bins provided
      for (int i = offset; i < static_cast<int>(edgesDeltaPhiOrig.size()) - offset; i++) {
        edgesDeltaPhi.emplace_back(edgesDeltaPhiOrig[i]);
      }
    } else {
      // fixed binning, generate the bin edges on-the-spot
      double min = edgesDeltaPhiOrig[0];
      double delta = (edgesDeltaPhiOrig[1] - edgesDeltaPhiOrig[0]) / preAxisDeltaPhi.nBins.value();
      for (int i = offset; i < preAxisDeltaPhi.nBins.value() + 1 - offset; i++) {
        edgesDeltaPhi.emplace_back(min + static_cast<double>(i) * delta);
      }
    }
    // ===] delta-eta [===
    if (!preAxisDeltaEta.nBins.has_value()) {
      // variable binning, use bins provided
      for (int i = offset; i < static_cast<int>(edgesDeltaEtaOrig.size()) - offset; i++) {
        edgesDeltaEta.emplace_back(edgesDeltaEtaOrig[i]);
      }
    } else {
      // fixed binning, generate the bin edges on-the-spot
      double min = edgesDeltaEtaOrig[0];
      double delta = (edgesDeltaEtaOrig[1] - edgesDeltaEtaOrig[0]) / preAxisDeltaEta.nBins.value();
      for (int i = offset; i < preAxisDeltaEta.nBins.value() + 1 - offset; i++) {
        edgesDeltaEta.emplace_back(min + static_cast<double>(i) * delta);
      }
    }
    // ===] pt assoc [===
    if (!preAxisPtAssoc.nBins.has_value()) {
      // variable binning, use bins provided
      for (int i = offset; i < static_cast<int>(edgesPtAssocOrig.size()) - offset; i++) {
        edgesPtAssoc.emplace_back(edgesPtAssocOrig[i]);
      }
    } else {
      // fixed binning, generate the bin edges on-the-spot
      double min = edgesPtAssocOrig[0];
      double delta = (edgesPtAssocOrig[1] - edgesPtAssocOrig[0]) / preAxisPtAssoc.nBins.value();
      for (int i = offset; i < preAxisPtAssoc.nBins.value() + 1 - offset; i++) {
        edgesPtAssoc.emplace_back(min + static_cast<double>(i) * delta);
      }
    }
    // ===] pt trigger [===
    if (!preAxisPtTrigger.nBins.has_value()) {
      // variable binning, use bins provided
      for (int i = offset; i < static_cast<int>(edgesPtTriggerOrig.size()) - offset; i++) {
        edgesPtTrigger.emplace_back(edgesPtTriggerOrig[i]);
      }
    } else {
      // fixed binning, generate the bin edges on-the-spot
      double min = edgesPtTriggerOrig[0];
      double delta = (edgesPtTriggerOrig[1] - edgesPtTriggerOrig[0]) / preAxisPtTrigger.nBins.value();
      for (int i = offset; i < preAxisPtTrigger.nBins.value() + 1 - offset; i++) {
        edgesPtTrigger.emplace_back(min + static_cast<double>(i) * delta);
      }
    }
    // ===] vtx Z [===
    if (!preAxisVtxZ.nBins.has_value()) {
      // variable binning, use bins provided
      for (int i = offset; i < static_cast<int>(edgesVtxZOrig.size()) - offset; i++) {
        edgesVtxZ.emplace_back(edgesVtxZOrig[i]);
      }
    } else {
      // fixed binning, generate the bin edges on-the-spot
      double min = edgesVtxZOrig[0];
      double delta = (edgesVtxZOrig[1] - edgesVtxZOrig[0]) / preAxisVtxZ.nBins.value();
      for (int i = offset; i < preAxisVtxZ.nBins.value() + 1 - offset; i++) {
        edgesVtxZ.emplace_back(min + static_cast<double>(i) * delta);
      }
    }
    // ===] mult percentile [===
    if (!preAxisMult.nBins.has_value()) {
      // variable binning, use bins provided
      for (int i = offset; i < static_cast<int>(edgesMultOrig.size()) - offset; i++) {
        edgesMult.emplace_back(edgesMultOrig[i]);
      }
    } else {
      // fixed binning, generate the bin edges on-the-spot
      double min = edgesMultOrig[0];
      double delta = (edgesMultOrig[1] - edgesMultOrig[0]) / preAxisMult.nBins.value();
      for (int i = offset; i < preAxisMult.nBins.value() + 1 - offset; i++) {
        edgesMult.emplace_back(min + static_cast<double>(i) * delta);
      }
    }
    // ===] multiplicity count [===
    if (!preAxisMultiplicity.nBins.has_value()) {
      // variable binning, use bins provided
      for (int i = offset; i < static_cast<int>(edgesMultiplicityOrig.size()) - offset; i++) {
        edgesMultiplicity.emplace_back(edgesMultiplicityOrig[i]);
      }
    } else {
      // fixed binning, generate the bin edges on-the-spot
      double min = edgesMultiplicityOrig[0];
      double delta = (edgesMultiplicityOrig[1] - edgesMultiplicityOrig[0]) / preAxisMultiplicity.nBins.value();
      for (int i = offset; i < preAxisMultiplicity.nBins.value() + 1 - offset; i++) {
        edgesMultiplicity.emplace_back(min + static_cast<double>(i) * delta);
      }
    }

    LOGF(info, "Initialized THnF axis delta-phi with %i bins.", edgesDeltaPhi.size() - 1);
    LOGF(info, "Initialized THnF axis delta-eta with %i bins.", edgesDeltaEta.size() - 1);
    LOGF(info, "Initialized THnF axis pTassoc with %i bins.", edgesPtAssoc.size() - 1);
    LOGF(info, "Initialized THnF axis pTtrigger with %i bins.", edgesPtTrigger.size() - 1);
    LOGF(info, "Initialized THnF axis vertex-Z with %i bins.", edgesVtxZ.size() - 1);
    LOGF(info, "Initialized THnF axis mult cent with %i bins.", edgesMult.size() - 1);
    LOGF(info, "Initialized THnF axis multiplicity with %i bins.", edgesMultiplicity.size() - 1);

    const AxisSpec axisDeltaPhiNDim{edgesDeltaPhi, "#Delta#varphi"};
    const AxisSpec axisDeltaEtaNDim{edgesDeltaEta, "#Delta#eta"};
    const AxisSpec axisPtAssocNDim{edgesPtAssoc, "#it{p}_{T}^{assoc} (GeV/c)"};
    const AxisSpec axisPtTriggerNDim{edgesPtTrigger, "#it{p}_{T}^{trigger} (GeV/c)"};
    const AxisSpec axisVtxZNDim{edgesVtxZ, "vertex Z (cm)"};
    const AxisSpec axisMultNDim{edgesMult, "mult percentile"};
    const AxisSpec axisMultiplicityNDim{edgesMultiplicity, "Multiplicity"};

    if (doprocessPairLossK0MC && pairLossK0Configurations.doStageDiagnostics) {
      const AxisSpec axisPairLossEventStage{6, -0.5, 5.5, "Event-selection stage"};
      const AxisSpec axisPairLossNRecoCollisions{11, -0.5, 10.5, "#it{N}_{reco collisions} per MC collision"};
      const AxisSpec axisPairLossStage{PairLossK0NStages, -0.5, static_cast<double>(PairLossK0NStages) - 0.5, "Reconstruction stage"};
      const AxisSpec axisPairLossFinalObjectState{4, -0.5, 3.5, "Final h-K^{0}_{S} object state"};
      const AxisSpec axisPairLossTrackLevelState{4, -0.5, 3.5, "Trigger and K^{0}_{S}-daughter track state"};
      const AxisSpec axisPairLossDaughterState{4, -0.5, 3.5, "K^{0}_{S}-daughter track state"};
      const AxisSpec axisPairLossFieldSign{3, -1.5, 1.5, "sgn(#it{B}_{z})"};
      const AxisSpec axisPairLossChargeProduct{3, -1.5, 1.5, "sgn(#it{q}_{h}#it{q}_{#pi})"};
      const AxisSpec axisPairLossMatchType{7, -0.5, 6.5, "Matched-object category"};
      const AxisSpec axisPairLossMatchCount{11, -0.5, 10.5, "#it{N}_{matched objects} per truth pair"};
      const AxisSpec axisPairLossTPCCrossedRows{161, -0.5, 160.5, "#it{N}_{TPC crossed rows}"};
      const AxisSpec axisPairLossTPCSharedClusters{161, -0.5, 160.5, "#it{N}_{TPC shared clusters}"};
      const AxisSpec axisPairLossITSClusters{8, -0.5, 7.5, "#it{N}_{ITS clusters}"};
      const AxisSpec axisPairLossLayer0{2, -0.5, 1.5, "Hit in ITS layer 0"};
      const AxisSpec axisPairLossMassNSigma{120, -12.0, 12.0, "#it{n}_{#sigma}(m_{#pi^{+}#pi^{-}}; K^{0}_{S})"};
      const AxisSpec axisPairLossCosPA{100, 0.9, 1.0, "cos(#theta_{PA})"};
      const AxisSpec axisPairLossDcaDaughters{100, 0.0, 2.0, "DCA(#pi^{+}, #pi^{-}) (cm)"};
      const AxisSpec axisPairLossDeltaPhiStarValidity{2, -0.5, 1.5, "Validity of min |#Delta#varphi^{*}_{h-#pi}|"};
      const AxisSpec axisPairLossTruthDeltaPhi{edgesDeltaPhi, "#Delta#varphi_{h-K^{0}_{S}}^{gen} (rad)"};
      const AxisSpec axisPairLossRecoDeltaPhi{edgesDeltaPhi, "#Delta#varphi_{h-K^{0}_{S}}^{reco} (rad)"};
      const AxisSpec axisPairLossTruthDeltaEta{edgesDeltaEta, "#Delta#eta_{h-K^{0}_{S}}^{gen}"};
      const AxisSpec axisPairLossRecoDeltaEta{edgesDeltaEta, "#Delta#eta_{h-K^{0}_{S}}^{reco}"};
      const AxisSpec axisPairLossTruthK0Pt{edgesPtAssoc, "#it{p}_{T,K^{0}_{S}}^{gen} (GeV/#it{c})"};
      const AxisSpec axisPairLossRecoK0Pt{edgesPtAssoc, "#it{p}_{T,K^{0}_{S}}^{reco} (GeV/#it{c})"};
      const AxisSpec axisPairLossTruthTriggerPt{edgesPtTrigger, "#it{p}_{T,h}^{gen} (GeV/#it{c})"};
      const AxisSpec axisPairLossRecoTriggerPt{edgesPtTrigger, "#it{p}_{T,h}^{reco} (GeV/#it{c})"};
      const AxisSpec axisPairLossMinDeltaPhiStar{pairLossK0Configurations.axisMinDeltaPhiStar, "min_{#it{r}} |#Delta#varphi^{*}_{h-#pi}(#it{r})| (rad)"};
      const AxisSpec axisPairLossSignedDeltaPhiStar{pairLossK0Configurations.axisSignedDeltaPhiStar, "#Delta#varphi^{*}_{h-#pi}(#it{r}_{min}) (rad)"};
      const AxisSpec axisPairLossDaughterDeltaEta{pairLossK0Configurations.axisDaughterDeltaEta, "#Delta#eta_{h-#pi}^{gen}"};
      const AxisSpec axisPairLossDaughterDeltaPhi{pairLossK0Configurations.axisDaughterDeltaPhi, "#Delta#varphi_{h-#pi}^{gen} (rad)"};
      const AxisSpec axisPairLossTruthDecayRadius{pairLossK0Configurations.axisDecayRadius, "#it{R}_{decay}^{gen}(K^{0}_{S}) (cm)"};
      const AxisSpec axisPairLossRecoV0Radius{pairLossK0Configurations.axisDecayRadius, "#it{R}_{V0}^{reco} (cm)"};
      const double pairLossRadiusMinimum = pairLossK0Configurations.minRadius;
      const double pairLossRadiusMaximum = std::max(static_cast<double>(pairLossK0Configurations.maxRadius), pairLossRadiusMinimum + 0.01);
      const double pairLossRadiusStep = std::max(static_cast<double>(pairLossK0Configurations.radiusStep), 0.001);
      const int pairLossRadiusBins = std::max(1, static_cast<int>(std::ceil((pairLossRadiusMaximum - pairLossRadiusMinimum) / pairLossRadiusStep)));
      const AxisSpec axisPairLossRadiusAtMinimum{pairLossRadiusBins, pairLossRadiusMinimum, pairLossRadiusMaximum, "#it{r}_{min} (m)"};
      const AxisSpec axisPairLossTriggerTracksFailureReason{PairLossTriggerTracksNReasons, -0.5, static_cast<double>(PairLossTriggerTracksNReasons) - 0.5, "TriggerTracks first-failing condition"};

      histos.add("PairLossK0/Event/hCounter", "K0 pair-loss event counter", kTH1F, {axisPairLossEventStage});
      histos.add("PairLossK0/Event/hNRecoCollisions", "reconstructed collisions per MC collision", kTH1F, {axisPairLossNRecoCollisions});
      histos.add("PairLossK0/Stage/hCounts", "pair-loss diagnostic stage counts", kTH1F, {axisPairLossStage});
      histos.add("PairLossK0/Stage/hCountsFindable", "pair-loss diagnostic stage counts for findable K0", kTH1F, {axisPairLossStage});
      histos.add("PairLossK0/Stage/hPhysics", "stages in h-K0 physics variables", kTHnF, {axisPairLossStage, axisPairLossTruthDeltaPhi, axisPairLossTruthDeltaEta, axisPairLossTruthK0Pt, axisPairLossTruthTriggerPt, axisMultNDim});
      histos.add("PairLossK0/Stage/hPhysicsFindable", "stages in h-K0 physics variables for findable K0", kTHnF, {axisPairLossStage, axisPairLossTruthDeltaPhi, axisPairLossTruthDeltaEta, axisPairLossTruthK0Pt, axisPairLossTruthTriggerPt, axisMultNDim});
      histos.add("PairLossK0/Stage/hClose", "stages in trigger-daughter close-pair variables", kTHnF, {axisPairLossStage, axisPairLossMinDeltaPhiStar, axisPairLossDaughterDeltaEta, axisPairLossTruthK0Pt, axisPairLossTruthTriggerPt, axisPairLossFieldSign, axisPairLossChargeProduct});
      histos.add("PairLossK0/Stage/hTriggerTracksFailureReason", "first-failing TriggerTracks condition for best-collision triggers, in h-K0 physics variables", kTHnF, {axisPairLossTriggerTracksFailureReason, axisPairLossTruthDeltaPhi, axisPairLossTruthDeltaEta, axisPairLossTruthK0Pt, axisPairLossTruthTriggerPt});
      // 1/N_trigger denominator of the stage ladder: same stage axis as hCounts/hPhysics, generated trigger pT, eta, phi.
      // Sparse because the eta-phi part would cost tens of MB dense for an occupancy of a few percent.
      histos.add("PairLossK0/Stage/hTriggers", "triggers per stage: the 1/#it{N}_{trigger} normalisation of the stage ladder", kTHnSparseF, {axisPairLossStage, axisPairLossTruthTriggerPt, axesConfigurations.axisEta, axesConfigurations.axisPhi, axisVtxZNDim, axisMultNDim});

      histos.add("PairLossK0/State/hFinalObjectStatePhysics", "00/01/10/11 final trigger-K0 object state", kTHnF, {axisPairLossFinalObjectState, axisPairLossTruthDeltaPhi, axisPairLossTruthDeltaEta, axisPairLossTruthK0Pt, axisPairLossTruthTriggerPt, axisPairLossFieldSign});
      histos.add("PairLossK0/State/hFinalObjectStateClose", "00/01/10/11 final trigger-K0 object state in close-pair variables", kTHnF, {axisPairLossFinalObjectState, axisPairLossMinDeltaPhiStar, axisPairLossDaughterDeltaEta, axisPairLossTruthK0Pt, axisPairLossTruthTriggerPt, axisPairLossFieldSign, axisPairLossChargeProduct});
      histos.add("PairLossK0/State/hTrackLevelStateClose", "00/01/10/11 trigger-track and both-daughters state", kTHnF, {axisPairLossTrackLevelState, axisPairLossMinDeltaPhiStar, axisPairLossDaughterDeltaEta, axisPairLossTruthK0Pt, axisPairLossTruthTriggerPt, axisPairLossFieldSign, axisPairLossChargeProduct});
      histos.add("PairLossK0/State/hDaughterTrackStateClose", "daughter state: none/negative/positive/both", kTHnF, {axisPairLossDaughterState, axisPairLossMinDeltaPhiStar, axisPairLossDaughterDeltaEta, axisPairLossTruthK0Pt, axisPairLossTruthTriggerPt, axisPairLossFieldSign, axisPairLossChargeProduct});

      histos.add("PairLossK0/Geometry/hGeneratedDaughters", "generated trigger-daughter geometry", kTHnF, {axisPairLossSignedDeltaPhiStar, axisPairLossDaughterDeltaEta, axisPairLossTruthK0Pt, axisPairLossTruthTriggerPt, axisPairLossFieldSign, axisPairLossChargeProduct});
      histos.add("PairLossK0/Geometry/hFinalDaughters", "geometry for final reconstructed h-K0 pairs", kTHnF, {axisPairLossSignedDeltaPhiStar, axisPairLossDaughterDeltaEta, axisPairLossTruthK0Pt, axisPairLossTruthTriggerPt, axisPairLossFieldSign, axisPairLossChargeProduct});
      histos.add("PairLossK0/Geometry/hRawDeltaPhiVsSignedDeltaPhiStar", "closest daughter raw delta phi versus signed delta phi star", kTH3F, {axisPairLossDaughterDeltaPhi, axisPairLossSignedDeltaPhiStar, axisPairLossFieldSign});
      histos.add("PairLossK0/Geometry/hGeneratedMinDeltaPhiStarVsDecayRadius", "generated minimum delta phi star versus K0 decay radius", kTH3F, {axisPairLossMinDeltaPhiStar, axisPairLossTruthDecayRadius, axisPairLossTruthK0Pt});
      histos.add("PairLossK0/Geometry/hFinalMinDeltaPhiStarVsDecayRadius", "final-pair minimum delta phi star versus K0 decay radius", kTH3F, {axisPairLossMinDeltaPhiStar, axisPairLossTruthDecayRadius, axisPairLossTruthK0Pt});
      histos.add("PairLossK0/Geometry/hGeneratedRadiusAtMinimum", "radius of the closest trigger-daughter approach", kTH3F, {axisPairLossRadiusAtMinimum, axisPairLossMinDeltaPhiStar, axisPairLossTruthK0Pt});
      histos.add("PairLossK0/Geometry/hFinalRadiusAtMinimum", "radius of the closest approach for final pairs", kTH3F, {axisPairLossRadiusAtMinimum, axisPairLossMinDeltaPhiStar, axisPairLossTruthK0Pt});
      histos.add("PairLossK0/Geometry/hDeltaPhiStarValidity", "minimum delta-phi-star validity", kTH2F, {axisPairLossDeltaPhiStarValidity, axisPairLossTruthK0Pt});
      histos.add("PairLossK0/Geometry/hAutocorrelationRejected", "pairs rejected only by reconstructed track-ID autocorrelation", kTHnF, {axisPairLossMinDeltaPhiStar, axisPairLossDaughterDeltaEta, axisPairLossTruthK0Pt, axisPairLossTruthTriggerPt, axisPairLossFieldSign});

      histos.add("PairLossK0/Matching/hNMatches", "number of matches for each truth pair", kTH2F, {axisPairLossMatchType, axisPairLossMatchCount});
      histos.add("PairLossK0/Response/hTriggerPt", "trigger pT response", kTH2F, {axisPairLossTruthTriggerPt, axisPairLossRecoTriggerPt});
      histos.add("PairLossK0/Response/hK0Pt", "K0 pT response", kTH2F, {axisPairLossTruthK0Pt, axisPairLossRecoK0Pt});
      histos.add("PairLossK0/Response/hDeltaPhi", "h-K0 delta-phi response", kTH2F, {axisPairLossTruthDeltaPhi, axisPairLossRecoDeltaPhi});
      histos.add("PairLossK0/Response/hDeltaEta", "h-K0 delta-eta response", kTH2F, {axisPairLossTruthDeltaEta, axisPairLossRecoDeltaEta});

      histos.add("PairLossK0/TrackQA/hTriggerTPCPairWeighted", "pair-weighted trigger TPC quality", kTH3F, {axisPairLossMinDeltaPhiStar, axisPairLossTPCCrossedRows, axisPairLossTPCSharedClusters});
      histos.add("PairLossK0/TrackQA/hTriggerITSPairWeighted", "pair-weighted trigger ITS quality", kTH3F, {axisPairLossMinDeltaPhiStar, axisPairLossITSClusters, axisPairLossLayer0});
      histos.add("PairLossK0/TrackQA/hPositiveDaughterTPCPairWeighted", "pair-weighted positive-daughter TPC quality", kTH3F, {axisPairLossMinDeltaPhiStar, axisPairLossTPCCrossedRows, axisPairLossTPCSharedClusters});
      histos.add("PairLossK0/TrackQA/hPositiveDaughterITSPairWeighted", "pair-weighted positive-daughter ITS quality", kTH3F, {axisPairLossMinDeltaPhiStar, axisPairLossITSClusters, axisPairLossLayer0});
      histos.add("PairLossK0/TrackQA/hNegativeDaughterTPCPairWeighted", "pair-weighted negative-daughter TPC quality", kTH3F, {axisPairLossMinDeltaPhiStar, axisPairLossTPCCrossedRows, axisPairLossTPCSharedClusters});
      histos.add("PairLossK0/TrackQA/hNegativeDaughterITSPairWeighted", "pair-weighted negative-daughter ITS quality", kTH3F, {axisPairLossMinDeltaPhiStar, axisPairLossITSClusters, axisPairLossLayer0});

      histos.add("PairLossK0/V0QA/hCosPAVsMinDeltaPhiStar", "matched V0 cos(PA) versus minimum delta phi star", kTH2F, {axisPairLossMinDeltaPhiStar, axisPairLossCosPA});
      histos.add("PairLossK0/V0QA/hDcaDaughtersVsMinDeltaPhiStar", "matched V0 daughter DCA versus minimum delta phi star", kTH2F, {axisPairLossMinDeltaPhiStar, axisPairLossDcaDaughters});
      histos.add("PairLossK0/V0QA/hRadiusVsMinDeltaPhiStar", "matched V0 radius versus minimum delta phi star", kTH2F, {axisPairLossMinDeltaPhiStar, axisPairLossRecoV0Radius});
      histos.add("PairLossK0/V0QA/hMassNSigmaVsMinDeltaPhiStar", "AssocV0 K0 mass n-sigma versus minimum delta phi star", kTH3F, {axisPairLossMinDeltaPhiStar, axisPairLossMassNSigma, axisPairLossTruthK0Pt});

      auto eventCounter = histos.get<TH1>(HIST("PairLossK0/Event/hCounter"));
      const std::array<std::string_view, 6> eventLabels = {"MC collision", "has rec collision", "best collision selected", "has truth trigger", "has truth K0", "has truth pair"};
      for (size_t i = 0; i < eventLabels.size(); ++i) {
        eventCounter->GetXaxis()->SetBinLabel(i + 1, eventLabels[i].data());
      }
      eventCounter->GetYaxis()->SetTitle("Entries");
      histos.get<TH1>(HIST("PairLossK0/Event/hNRecoCollisions"))->GetYaxis()->SetTitle("MC collisions");
      auto setStageAxisLabels = [&](TAxis* stageAxis) {
        for (int i = 0; i < PairLossK0NStages; ++i) {
          stageAxis->SetBinLabel(i + 1, PairLossK0StageNames[i].data());
        }
      };
      auto setStageLabels = [&](auto const& stageHistogram) {
        setStageAxisLabels(stageHistogram->GetXaxis());
      };
      auto stageCounts = histos.get<TH1>(HIST("PairLossK0/Stage/hCounts"));
      auto stageCountsFindable = histos.get<TH1>(HIST("PairLossK0/Stage/hCountsFindable"));
      setStageLabels(stageCounts);
      setStageLabels(stageCountsFindable);
      setStageAxisLabels(histos.get<THnSparse>(HIST("PairLossK0/Stage/hTriggers"))->GetAxis(0));
      stageCounts->GetYaxis()->SetTitle("Truth-pair entries");
      stageCountsFindable->GetYaxis()->SetTitle("Findable truth-pair entries");

      auto triggerTracksFailureHistogram = histos.get<THn>(HIST("PairLossK0/Stage/hTriggerTracksFailureReason"));
      for (int i = 0; i < PairLossTriggerTracksNReasons; ++i) {
        triggerTracksFailureHistogram->GetAxis(0)->SetBinLabel(i + 1, PairLossTriggerTracksFailureNames[i].data());
      }

      const std::array<std::string_view, 4> objectStateLabels = {"00 neither", "01 K0 only", "10 trigger only", "11 both"};
      for (auto const& histogram : {histos.get<THn>(HIST("PairLossK0/State/hFinalObjectStatePhysics")), histos.get<THn>(HIST("PairLossK0/State/hFinalObjectStateClose"))}) {
        for (size_t i = 0; i < objectStateLabels.size(); ++i) {
          histogram->GetAxis(0)->SetBinLabel(i + 1, objectStateLabels[i].data());
        }
      }
      const std::array<std::string_view, 4> trackLevelStateLabels = {"00 neither", "01 both daughters only", "10 trigger only", "11 trigger + both daughters"};
      auto trackLevelStateHistogram = histos.get<THn>(HIST("PairLossK0/State/hTrackLevelStateClose"));
      for (size_t i = 0; i < trackLevelStateLabels.size(); ++i) {
        trackLevelStateHistogram->GetAxis(0)->SetBinLabel(i + 1, trackLevelStateLabels[i].data());
      }
      const std::array<std::string_view, 4> daughterStateLabels = {"none", "negative only", "positive only", "both"};
      auto daughterStateHistogram = histos.get<THn>(HIST("PairLossK0/State/hDaughterTrackStateClose"));
      for (size_t i = 0; i < daughterStateLabels.size(); ++i) {
        daughterStateHistogram->GetAxis(0)->SetBinLabel(i + 1, daughterStateLabels[i].data());
      }

      auto deltaPhiStarValidityHistogram = histos.get<TH2>(HIST("PairLossK0/Geometry/hDeltaPhiStarValidity"));
      deltaPhiStarValidityHistogram->GetXaxis()->SetBinLabel(1, "invalid");
      deltaPhiStarValidityHistogram->GetXaxis()->SetBinLabel(2, "valid");
      const std::array<std::string_view, 7> matchLabels = {"trigger all collisions", "trigger best collision", "trigger final", "V0 all collisions", "V0 best collision", "V0 AssocV0s", "V0 final"};
      auto matchHistogram = histos.get<TH2>(HIST("PairLossK0/Matching/hNMatches"));
      for (size_t i = 0; i < matchLabels.size(); ++i) {
        matchHistogram->GetXaxis()->SetBinLabel(i + 1, matchLabels[i].data());
      }

      for (auto const& histogram : {histos.get<TH3>(HIST("PairLossK0/TrackQA/hTriggerITSPairWeighted")),
                                    histos.get<TH3>(HIST("PairLossK0/TrackQA/hPositiveDaughterITSPairWeighted")),
                                    histos.get<TH3>(HIST("PairLossK0/TrackQA/hNegativeDaughterITSPairWeighted"))}) {
        histogram->GetZaxis()->SetBinLabel(1, "no");
        histogram->GetZaxis()->SetBinLabel(2, "yes");
      }

      for (auto const& histogram : {deltaPhiStarValidityHistogram,
                                    matchHistogram,
                                    histos.get<TH2>(HIST("PairLossK0/Response/hTriggerPt")),
                                    histos.get<TH2>(HIST("PairLossK0/Response/hK0Pt")),
                                    histos.get<TH2>(HIST("PairLossK0/Response/hDeltaPhi")),
                                    histos.get<TH2>(HIST("PairLossK0/Response/hDeltaEta")),
                                    histos.get<TH2>(HIST("PairLossK0/V0QA/hCosPAVsMinDeltaPhiStar")),
                                    histos.get<TH2>(HIST("PairLossK0/V0QA/hDcaDaughtersVsMinDeltaPhiStar")),
                                    histos.get<TH2>(HIST("PairLossK0/V0QA/hRadiusVsMinDeltaPhiStar"))}) {
        histogram->GetZaxis()->SetTitle("Entries");
      }
    }

    if (doprocessPairLossK0MC && pairLossK0Configurations.doRecComparison) {
      if (doprocessSameEventHV0s) {
        LOGF(fatal, "pairLossK0Configurations.doRecComparison already runs the exact Rec path internally; set processSameEventHV0s=false to avoid double filling");
      }
      constexpr int PairLossComparisonNVariants = 16;
      const AxisSpec axisPairLossComparisonVariant{PairLossComparisonNVariants, -0.5, PairLossComparisonNVariants - 0.5, "cumulative Rec control variant"};
      constexpr int PairLossComparisonNTriggerVariants = 10;
      const AxisSpec axisPairLossComparisonTriggerVariant{PairLossComparisonNTriggerVariants, -0.5, PairLossComparisonNTriggerVariants - 0.5, "cumulative Rec trigger control variant"};
      histos.add("PairLossK0/Comparison/Rec", "exact Rec peak-window h-K0 pairs; all candidates", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
      histos.add("PairLossK0/Comparison/RecPeakTrueK0", "exact Rec peak-window pairs with a true K0", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
      histos.add("PairLossK0/Comparison/RecPeakFakeK0", "exact Rec peak-window pairs without a true K0", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
      histos.add("PairLossK0/Comparison/RecLeftBg", "exact Rec left invariant-mass sideband pairs", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
      histos.add("PairLossK0/Comparison/RecRightBg", "exact Rec right invariant-mass sideband pairs", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
      histos.add("PairLossK0/Comparison/Truth", "PairLoss truth pairs; same axes as Rec", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
      histos.add("PairLossK0/Comparison/Gen", "Closure-test Gen pairs; same axes as Rec", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
      histos.add("PairLossK0/Comparison/Final", "PairLoss final truth-matched pairs; same axes as Rec", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
      histos.add("PairLossK0/Comparison/RecVariants", "cumulative controlled variants of the Rec pair definition", kTHnF, {axisPairLossComparisonVariant, axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
      // The two directional set differences between Rec variant 14 and the
      // exact old Final definition. They make non-nested pair definitions
      // explicit instead of hiding the mismatch in a single yield ratio.
      histos.add("PairLossK0/Comparison/FinalNotInRec", "old-Final truth pairs absent from Rec variant 14", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
      histos.add("PairLossK0/Comparison/RecNotInFinal", "Rec variant-14 truth pairs absent from old Final", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
      // Trigger-side counterpart of the pair ladder. Pair variant i has to be
      // divided by the trigger variant listed in the label to obtain a
      // per-trigger yield that is narrowed by the same cumulative conditions.
      histos.add("PairLossK0/Comparison/TriggerVariants", "cumulative controlled variants of the Rec trigger normalisation", kTHnF, {axisPairLossComparisonTriggerVariant, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});

      const std::array<std::string_view, PairLossComparisonNVariants> variantLabels = {
        "0 exact Rec",
        "1 + truth-trigger match/acceptance",
        "2 + true K0 (remove peak background)",
        "3 + truth-K0 match/acceptance",
        "4 + truth daughter exclusion",
        "5 + truth-trigger primary",
        "6 + truth-K0 primary",
        "7 + old truth-leading choice",
        "8 + best collision",
        "9 + truth pT coordinates only",
        "10 + truth event coordinates (vtxZ, mult)",
        "11 + truth delta eta only",
        "12 + truth delta phi and truth angular range",
        "13 + unit pair weight",
        "14 + one entry per truth pair",
        "15 exact old Final reference"};
      auto recVariants = histos.get<THn>(HIST("PairLossK0/Comparison/RecVariants"));
      for (size_t i = 0; i < variantLabels.size(); ++i) {
        recVariants->GetAxis(0)->SetBinLabel(i + 1, variantLabels[i].data());
      }

      // Pair variant -> trigger variant pairing for per-trigger yields:
      // pair 0-3 with trigger 0-1, pair 4-6 with trigger 2, pair 7 with
      // trigger 3, pair 8-9 with trigger 4, pair 10-12 with trigger 6,
      // pair 13-14 with trigger 7-8, pair 15 with trigger 9.
      const std::array<std::string_view, PairLossComparisonNTriggerVariants> triggerVariantLabels = {
        "0 exact Rec trigger",
        "1 + truth-trigger match/acceptance",
        "2 + truth-trigger primary",
        "3 + old truth-leading choice",
        "4 + best collision",
        "5 + truth pT coordinate only",
        "6 + truth event coordinates (vtxZ, mult)",
        "7 + unit trigger weight",
        "8 + one entry per truth trigger",
        "9 exact old Final trigger reference"};
      auto triggerVariants = histos.get<THn>(HIST("PairLossK0/Comparison/TriggerVariants"));
      for (size_t i = 0; i < triggerVariantLabels.size(); ++i) {
        triggerVariants->GetAxis(0)->SetBinLabel(i + 1, triggerVariantLabels[i].data());
      }
    }

    if (doprocessPairLossK0MC && pairLossK0Configurations.doGenLevelStudy) {
      const AxisSpec axisGenStudyNch{pairLossK0Configurations.axisGenStudyNch, "#it{N}_{ch}^{gen} (|#eta| < 0.8)"};
      const AxisSpec axisGenStudyEventStage{4, -0.5, 3.5, "Generated-event selection stage"};
      // Findability of the K0, in exactly the sense the stage ladder uses: it
      // decayed to pi+ pi- and both charged daughters are inside the tracking
      // acceptance set by daughterPtMin / daughterEtaMax. Kept as an axis rather
      // than as a separate folder so that the inclusive and the findable-only
      // answer come out of one and the same object.
      const AxisSpec axisGenStudyFindable{2, -0.5, 1.5, "K^{0}_{S} findable"};

      histos.add("PairLossK0/GenStudy/hEventCounter", "generator-level event selection", kTH1F, {axisGenStudyEventStage});
      histos.add("PairLossK0/GenStudy/hNch", "generated charged multiplicity of selected MC collisions", kTH1F, {axisGenStudyNch});
      histos.add("PairLossK0/GenStudy/hNRecoCollisions", "reconstructed collisions per selected MC collision", kTH1F, {{11, -0.5, 10.5}});

      // Gen/ holds every generated object that passes the generator-level
      // selection; Reconstructed/ and NotReconstructed/ split that same set by
      // whether the object has a reconstructed counterpart. All three are filled
      // with generated coordinates, so Gen == Reconstructed + NotReconstructed
      // bin by bin and NotReconstructed/Gen reads directly as the loss.
      histos.add("PairLossK0/GenStudy/Gen/hTrigger", "generated triggers;#it{p}_{T}^{gen} (GeV/#it{c});#eta^{gen};#varphi^{gen};#it{N}_{ch}^{gen}", kTHnF, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi, axisGenStudyNch});
      histos.add("PairLossK0/GenStudy/Gen/hK0Short", "generated K0s;#it{p}_{T}^{gen} (GeV/#it{c});#eta^{gen};#varphi^{gen};#it{N}_{ch}^{gen};findable", kTHnF, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi, axisGenStudyNch, axisGenStudyFindable});
      // h-K0 correlations of the very same objects. Gen/ is every generated pair and the
      // four exclusive classes below split it by which of the two objects was
      // reconstructed, so Reconstructed + OnlyTriggerReconstructed + OnlyK0Reconstructed
      // + NotReconstructed equals Gen bin by bin.
      histos.add("PairLossK0/GenStudy/Gen/hCorrelation", "generated h-K0s pairs;#Delta#eta;#Delta#varphi;#it{p}_{T}^{trigger} (GeV/#it{c});#it{p}_{T}^{K^{0}_{S}} (GeV/#it{c});#it{N}_{ch}^{gen}", kTHnF, {axisDeltaEtaNDim, axisDeltaPhiNDim, axisPtTriggerNDim, axisPtAssocNDim, axisGenStudyNch});
      histos.addClone("PairLossK0/GenStudy/Gen/", "PairLossK0/GenStudy/Reconstructed/");
      histos.addClone("PairLossK0/GenStudy/Gen/", "PairLossK0/GenStudy/NotReconstructed/");
      // Only the correlation exists for the two mixed classes -- a single particle is
      // either reconstructed or not, so cloning the single-particle folders here would
      // only produce histograms with no meaning.
      histos.add("PairLossK0/GenStudy/OnlyTriggerReconstructed/hCorrelation", "h-K0s pairs with only the trigger reconstructed;#Delta#eta;#Delta#varphi;#it{p}_{T}^{trigger} (GeV/#it{c});#it{p}_{T}^{K^{0}_{S}} (GeV/#it{c});#it{N}_{ch}^{gen}", kTHnF, {axisDeltaEtaNDim, axisDeltaPhiNDim, axisPtTriggerNDim, axisPtAssocNDim, axisGenStudyNch});
      histos.add("PairLossK0/GenStudy/OnlyK0Reconstructed/hCorrelation", "h-K0s pairs with only the K0s reconstructed;#Delta#eta;#Delta#varphi;#it{p}_{T}^{trigger} (GeV/#it{c});#it{p}_{T}^{K^{0}_{S}} (GeV/#it{c});#it{N}_{ch}^{gen}", kTHnF, {axisDeltaEtaNDim, axisDeltaPhiNDim, axisPtTriggerNDim, axisPtAssocNDim, axisGenStudyNch});

      for (auto const& histogram : {histos.get<THn>(HIST("PairLossK0/GenStudy/Gen/hK0Short")),
                                    histos.get<THn>(HIST("PairLossK0/GenStudy/Reconstructed/hK0Short")),
                                    histos.get<THn>(HIST("PairLossK0/GenStudy/NotReconstructed/hK0Short"))}) {
        histogram->GetAxis(4)->SetBinLabel(1, "not findable");
        histogram->GetAxis(4)->SetBinLabel(2, "findable");
      }

      auto genStudyEventCounter = histos.get<TH1>(HIST("PairLossK0/GenStudy/hEventCounter"));
      const std::array<std::string_view, 4> genStudyEventLabels = {"MC collisions", "INEL>0 (generated)", "|vtx z| < cut (generated)", "has >= 1 rec collision"};
      for (size_t i = 0; i < genStudyEventLabels.size(); ++i) {
        genStudyEventCounter->GetXaxis()->SetBinLabel(i + 1, genStudyEventLabels[i].data());
      }
      genStudyEventCounter->GetYaxis()->SetTitle("MC collisions");
    }

    if (doprocessMixedEventHV0sInBuffer || doprocessMixedEventHCascadesInBuffer) {
      validCollisions.resize(histos.get<TH1>(HIST("axes/hMultAxis"))->GetNbinsX() * histos.get<TH1>(HIST("axes/hVertexZAxis"))->GetNbinsX());
      for (size_t i = 0; i < validCollisions.size(); ++i) {
        validCollisions[i].reserve(masterConfigurations.mixingParameter);
      }
    }
    if (!masterConfigurations.doPPAnalysis) {
      // event selections in Pb-Pb
      histos.add("hEventSelection", "hEventSelection", kTH1F, {{10, 0, 10}});
      std::array<TString, 10> eventSelLabel = {"all", "sel8", "kIsTriggerTVX", "PV_{z}", "kIsGoodITSLayersAll", "kIsGoodZvtxFT0vsPV", "OccupCut", "kNoTimeFrameBorder", "kNoITSROFrameBorder", "kNoSameBunchPileup "};
      for (int i = 1; i <= histos.get<TH1>(HIST("hEventSelection"))->GetNbinsX(); i++) {
        histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(i, eventSelLabel[i - 1].Data());
      }
    }
    // Some QA plots
    if (doprocessMCGenerated) {
      histos.add("hGeneratedQAPtTrigger", "hGeneratedQAPtTrigger", kTH2F, {axesConfigurations.axisPtQA, {5, -0.5f, 4.5f}});
      histos.add("hGeneratedQAPtAssociatedK0", "hGeneratedQAPtAssociatedK0", kTH2F, {axesConfigurations.axisPtQA, {5, -0.5f, 4.5f}});
    }
    if (doprocessClosureTest) {
      histos.add("hClosureQAPtTrigger", "hClosureQAPtTrigger", kTH2F, {axesConfigurations.axisPtQA, {5, -0.5f, 4.5f}});
      histos.add("hClosureQAPtAssociatedK0", "hClosureQAPtAssociatedK0", kTH2F, {axesConfigurations.axisPtQA, {5, -0.5f, 4.5f}});
    }
    if (doprocessMCGenerated || doprocessClosureTest) {
      histos.add("hClosureTestEventCounter", "hClosureTestEventCounter", kTH1F, {{10, 0, 10}});
    }
    if (doprocessSameEventHV0s || doprocessPairLossK0MC || doprocessSameEventHCascades || doprocessSameEventHPions || doprocessSameEventHHadrons) {
      histos.add("hTriggerAllSelectedEtaVsPt", "hTriggerAllSelectedEtaVsPt", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
      histos.add("hPositiveTriggerPrimaryEtaVsPt", "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
      histos.add("hNegativeTriggerPrimaryEtaVsPt", "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
      // QA and THn Histograms
      histos.add("hTriggerPtResolution", ";p_{T}^{reconstructed} (GeV/c); p_{T}^{generated} (GeV/c)", kTH2F, {axesConfigurations.axisPtQA, axesConfigurations.axisPtQA});
      if (!efficiencyFlags.applyEffAsFunctionOfMultAndPhi) {
        histos.add("hTriggerPrimaryEtaVsPt", "hTriggerPrimaryEtaVsPt", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
      } else {
        histos.add("hTriggerPrimaryEtaVsPt", "hTriggerPrimaryEtaVsPt", kTHnF, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi, axesConfigurations.axisMult});
      }
      histos.add("hTrackEtaVsPtVsPhi", "hTrackEtaVsPtVsPhi", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi});
      histos.add("hAssocTrackEtaVsPtVsPhi", "hAssocTrackEtaVsPtVsPhi", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi});
      if (masterConfigurations.doLocalDensityStudy) {
        histos.add("hTrackEtaVsPtVsLocalDensity", "hTrackEtaVsPtVsLocalDensity", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisLocalDensity});
      }
      // histos.add("hTrackAttempt", "Attempt", kTH3F, {axisPtQA, axisEta, axisPhi});
    }
    if (doprocessSameEventHPions || doprocessSameEventHHadrons || doprocessMixedEventHPions || doprocessMixedEventHHadrons) {
      histos.add("hNumberOfRejectedPairsHadron", "hNumberOfRejectedPairsHadron", kTH1F, {{1, 0, 1}});
      histos.add("hNumberOfRejectedPairsPion", "hNumberOfRejectedPairsPion", kTH1F, {{1, 0, 1}});
    }
    if (doprocessSameEventHHadrons) {
      histos.add("hDCAzTriggerHadron", "hDCAzTriggerHadron", kTH2F, {{200, -0.5, 0.5}, axesConfigurations.axisPtQA});
      histos.add("hDCAxyTriggerHadron", "hDCAxyTriggerHadron", kTH2F, {{200, -0.5, 0.5}, axesConfigurations.axisPtQA});
      histos.add("hDCAzAssociatedHadron", "hDCAzAssociatedHadron", kTH2F, {{200, -0.5, 0.5}, axesConfigurations.axisPtQA});
      histos.add("hDCAxyAssociatedHadron", "hDCAxyAssociatedHadron", kTH2F, {{200, -0.5, 0.5}, axesConfigurations.axisPtQA});
    }
    if (doprocessSameEventHV0s || doprocessPairLossK0MC || doprocessMixedEventHV0s) {
      histos.add("hNumberOfRejectedPairsV0", "hNumberOfRejectedPairsV0", kTH1F, {{1, 0, 1}});
    }
    if (doprocessSameEventHCascades || doprocessMixedEventHCascades) {
      histos.add("hNumberOfRejectedPairsCascades", "hNumberOfRejectedPairsCascades", kTH1F, {{1, 0, 1}});
    }

    if (doMixingQAandEventQA) {
      // mixing QA
      histos.add("MixingQA/hSECollisionBins", ";bin;Entries", kTH1F, {{140, -0.5, 139.5}});
      histos.add("MixingQA/hMECollisionBins", ";bin;Entries", kTH1F, {{140, -0.5, 139.5}});
      histos.add("MixingQA/hMEpvz1", ";pvz;Entries", kTH1F, {{30, -15, 15}});
      histos.add("MixingQA/hMEpvz2", ";pvz;Entries", kTH1F, {{30, -15, 15}});

      // Event QA
      histos.add("EventQA/hMixingQA", "mixing QA", kTH1F, {{2, -0.5, 1.5}});
      histos.add("EventQA/hMult", "Multiplicity", kTH1F, {axesConfigurations.axisMult});
      histos.add("EventQA/hPvz", ";pvz;Entries", kTH1F, {{30, -15, 15}});
      histos.add("EventQA/hMultFT0vsTPC", ";centFT0M;multNTracksPVeta1", kTH2F, {{100, 0, 100}, {300, 0, 300}});
    }
    bool hStrange = false;
    for (int i = 0; i < AssocParticleTypes; i++) {
      if (TESTBIT(doCorrelation, i)) {
        if (masterConfigurations.doFullCorrelationStudy) {
          histos.add(fmt::format("sameEvent/Signal/{}", Particlenames[i]).c_str(), "", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
        }
        if (doDeltaPhiStarCheck && masterConfigurations.doFullCorrelationStudy) {
          histos.add(fmt::format("sameEvent/Signal/{}DeltaPhiStar", Particlenames[i]).c_str(), "", kTH3F, {{100, -0.3, 0.3}, {50, -0.05, 0.05}, {2, -1, 1}}); // -1 oposite charge, 1 same charge
        }
        if (i < IndexPion) {
          if (!efficiencyFlags.applyEffAsFunctionOfMultAndPhi) {
            histos.add(fmt::format("h{}EtaVsPtVsPhi", Particlenames[i]).c_str(), "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi});
            histos.add(fmt::format("h{}EtaVsPtVsPhiBg", Particlenames[i]).c_str(), "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi});
          } else {
            histos.add(fmt::format("h{}EtaVsPtVsPhiVsCent", Particlenames[i]).c_str(), "", kTHnF, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi, axesConfigurations.axisMult});
            histos.add(fmt::format("h{}EtaVsPtVsPhiVsCentBg", Particlenames[i]).c_str(), "", kTHnF, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi, axesConfigurations.axisMult});
          }
          if (masterConfigurations.doLocalDensityStudy && i < AssocV0Types) {
            histos.add(fmt::format("h{}EtaVsPtVsLocalDensity", Particlenames[i]).c_str(), "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisLocalDensity});
            histos.add(fmt::format("h{}EtaVsPtVsLocalDensityBg", Particlenames[i]).c_str(), "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisLocalDensity});
          }
          histos.add(fmt::format("h3d{}Spectrum", Particlenames[i]).c_str(), fmt::format("h3d{}Spectrum", Particlenames[i]).c_str(), kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisMult, axesConfigurations.axisMassNSigma});
          histos.add(fmt::format("h3d{}SpectrumY", Particlenames[i]).c_str(), fmt::format("h3d{}SpectrumY", Particlenames[i]).c_str(), kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisMult, axesConfigurations.axisMassNSigma});
          hStrange = true;
          if (doITSClustersQA) {
            histos.add(fmt::format("hITSClusters{}NegativeDaughterToward", Particlenames[i]).c_str(), "", kTH3F, {axesConfigurations.axisPtAssoc, {8, -0.5, 7.5}, {20, 0, 10}});
            histos.add(fmt::format("hITSClusters{}PositiveDaughterToward", Particlenames[i]).c_str(), "", kTH3F, {axesConfigurations.axisPtAssoc, {8, -0.5, 7.5}, {20, 0, 10}});
            histos.add(fmt::format("hITSClusters{}NegativeDaughterTransverse", Particlenames[i]).c_str(), "", kTH3F, {axesConfigurations.axisPtAssoc, {8, -0.5, 7.5}, {20, 0, 10}});
            histos.add(fmt::format("hITSClusters{}PositiveDaughterTransverse", Particlenames[i]).c_str(), "", kTH3F, {axesConfigurations.axisPtAssoc, {8, -0.5, 7.5}, {20, 0, 10}});
          }
        }
      }
    }
    if (TESTBIT(doCorrelation, 7) && doprocessSameEventHPions) {
      histos.add("hPionEtaVsPtAllSelected", "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
      histos.add("hPionEtaVsPt", "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
      histos.add("hPositivePionEtaVsPt", "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
      histos.add("hNegativePionEtaVsPt", "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
    }
    if (TESTBIT(doCorrelation, 8) && doprocessSameEventHHadrons) {
      histos.add("hAsssocTrackEtaVsPtVsPhi", "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi});
      histos.add("hAssocPrimaryEtaVsPt", "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
      histos.add("hAssocHadronsAllSelectedEtaVsPt", "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
      histos.add("hAssocPtResolution", ";p_{T}^{reconstructed} (GeV/c); p_{T}^{generated} (GeV/c)", kTH2F, {axesConfigurations.axisPtQA, axesConfigurations.axisPtQA});
    }

    if (hStrange && masterConfigurations.doFullCorrelationStudy) {
      histos.addClone("sameEvent/Signal/", "sameEvent/LeftBg/");
      histos.addClone("sameEvent/Signal/", "sameEvent/RightBg/");
    }

    if (masterConfigurations.doFullCorrelationStudy && (doprocessSameEventHV0s || doprocessPairLossK0MC) && masterConfigurations.doCorrelationsHadronV0daughter) {
      if (TESTBIT(doCorrelation, 0)) {
        histos.add("sameEvent/Signal/K0Short_hSameSign", "", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
        histos.add("sameEvent/Signal/K0Short_hOppositeSign", "", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
      }
      if (TESTBIT(doCorrelation, 1)) {
        histos.add("sameEvent/Signal/Lambda_hSameSign", "", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
        histos.add("sameEvent/Signal/Lambda_hOppositeSign", "", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
      }
      if (TESTBIT(doCorrelation, 2)) {
        histos.add("sameEvent/Signal/AntiLambda_hSameSign", "", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
        histos.add("sameEvent/Signal/AntiLambda_hOppositeSign", "", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
      }
    }
    if (TESTBIT(doCorrelation, 0) && (doprocessSameEventHV0s || doprocessPairLossK0MC) && masterConfigurations.doMassSpectrumCheck) {
      histos.add("hK0ShortPtVsMass", "", kTH2F, {axesConfigurations.axisPtQA, axesConfigurations.axisK0ShortMass});
    }
    if (TESTBIT(doCorrelation, 1) && (doprocessSameEventHV0s || doprocessPairLossK0MC) && masterConfigurations.doMassSpectrumCheck) {
      histos.add("hLambdaPtVsMass", "", kTH2F, {axesConfigurations.axisPtQA, axesConfigurations.axisLambdaMass});
    }
    if (TESTBIT(doCorrelation, 2) && (doprocessSameEventHV0s || doprocessPairLossK0MC) && masterConfigurations.doMassSpectrumCheck) {
      histos.add("hAntiLambdaPtVsMass", "", kTH2F, {axesConfigurations.axisPtQA, axesConfigurations.axisLambdaMass});
    }
    LOGF(info, "Init THnFs done");
    // mixed-event correlation functions
    if ((doprocessMixedEventHV0sInBuffer || doprocessMixedEventHCascadesInBuffer || doprocessMixedEventHV0s || doprocessMixedEventHCascades || doprocessMixedEventHPions || doprocessMixedEventHHadrons) && masterConfigurations.doFullCorrelationStudy) {
      histos.addClone("sameEvent/", "mixedEvent/");
    }
    if (masterConfigurations.doFullCorrelationStudy && (doprocessSameEventHV0s || doprocessPairLossK0MC)) {
      if (TESTBIT(doCorrelation, 0)) {
        histos.add("sameEvent/InvariantMass/K0Short/hNearSide", "", kTH3F, {axesConfigurations.axisPtAssoc, axesConfigurations.axisPtTrigger, axesConfigurations.axisK0ShortMass});
        histos.add("sameEvent/InvariantMass/K0Short/hAwaySide", "", kTH3F, {axesConfigurations.axisPtAssoc, axesConfigurations.axisPtTrigger, axesConfigurations.axisK0ShortMass});
        histos.add("sameEvent/InvariantMass/K0Short/hUnderlyingEvent", "", kTH3F, {axesConfigurations.axisPtAssoc, axesConfigurations.axisPtTrigger, axesConfigurations.axisK0ShortMass});
      }
      if (TESTBIT(doCorrelation, 1)) {
        histos.add("sameEvent/InvariantMass/Lambda/hNearSide", "", kTH3F, {axesConfigurations.axisPtAssoc, axesConfigurations.axisPtTrigger, axesConfigurations.axisLambdaMass});
        histos.add("sameEvent/InvariantMass/Lambda/hAwaySide", "", kTH3F, {axesConfigurations.axisPtAssoc, axesConfigurations.axisPtTrigger, axesConfigurations.axisLambdaMass});
        histos.add("sameEvent/InvariantMass/Lambda/hUnderlyingEvent", "", kTH3F, {axesConfigurations.axisPtAssoc, axesConfigurations.axisPtTrigger, axesConfigurations.axisLambdaMass});
      }
      if (TESTBIT(doCorrelation, 2)) {
        histos.add("sameEvent/InvariantMass/AntiLambda/hNearSide", "", kTH3F, {axesConfigurations.axisPtAssoc, axesConfigurations.axisPtTrigger, axesConfigurations.axisLambdaMass});
        histos.add("sameEvent/InvariantMass/AntiLambda/hAwaySide", "", kTH3F, {axesConfigurations.axisPtAssoc, axesConfigurations.axisPtTrigger, axesConfigurations.axisLambdaMass});
        histos.add("sameEvent/InvariantMass/AntiLambda/hUnderlyingEvent", "", kTH3F, {axesConfigurations.axisPtAssoc, axesConfigurations.axisPtTrigger, axesConfigurations.axisLambdaMass});
      }
    }
    if (doprocessSameEventHHadrons && masterConfigurations.doFullCorrelationStudy) {
      histos.add("sameEvent/TriggerParticlesHadron", "TriggersHadron", kTH2F, {axesConfigurations.axisPtQA, axesConfigurations.axisMult});
    }
    if ((doprocessSameEventHV0s || doprocessPairLossK0MC) && masterConfigurations.doFullCorrelationStudy) {
      histos.add("sameEvent/TriggerParticlesV0", "TriggersV0", kTH2F, {axesConfigurations.axisPtQA, axesConfigurations.axisMult});
    }
    if (doprocessSameEventHCascades && masterConfigurations.doFullCorrelationStudy) {
      histos.add("sameEvent/TriggerParticlesCascade", "TriggersCascade", kTH2F, {axesConfigurations.axisPtQA, axesConfigurations.axisMult});
    }
    if (doprocessSameEventHPions && masterConfigurations.doFullCorrelationStudy) {
      histos.add("sameEvent/TriggerParticlesPion", "TriggersPion", kTH2F, {axesConfigurations.axisPtQA, axesConfigurations.axisMult});
    }

    // MC generated plots
    if (doprocessMCGenerated) {
      if (!efficiencyFlags.applyEffAsFunctionOfMultAndPhi) {
        histos.add("Generated/hTrigger", "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
      } else {
        histos.add("Generated/hTrigger", "", kTHnF, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi, axesConfigurations.axisMult});
      }
      histos.add("Generated/hPositiveTrigger", "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
      histos.add("Generated/hNegativeTrigger", "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
      for (int i = 0; i < AssocParticleTypes; i++) {
        if (!efficiencyFlags.applyEffAsFunctionOfMultAndPhi) {
          histos.add(fmt::format("Generated/h{}", Particlenames[i]).c_str(), "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
        } else {
          histos.add(fmt::format("Generated/h{}", Particlenames[i]).c_str(), "", kTHnF, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi, axesConfigurations.axisMult});
        }
        if (i == IndexPion) {
          histos.add(fmt::format("Generated/hPositive{}", Particlenames[i]).c_str(), "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
          histos.add(fmt::format("Generated/hNegative{}", Particlenames[i]).c_str(), "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
        }
      }
      histos.addClone("Generated/", "GeneratedWithPV/");

      // The density axis is the reconstructed one on both sides: these generated histograms are counted
      // from the tracks of the best collision, exactly like their reconstructed counterparts.
      if (masterConfigurations.doLocalDensityStudy && masterConfigurations.doPPAnalysis) {
        histos.add("GeneratedWithPV/hTriggerLocalDensity", "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisLocalDensity});
        for (int i = 0; i < AssocParticleTypesNoHadron; i++) {
          histos.add(fmt::format("GeneratedWithPV/h{}LocalDensity", Particlenames[i]).c_str(), "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisLocalDensity});
        }
      }

      // histograms within |y|<0.5, vs multiplicity
      for (int i = 0; i < AssocParticleTypesNoHadron; i++) {
        histos.add(fmt::format("GeneratedWithPV/h{}_MidYVsMult", Particlenames[i]).c_str(), "", kTH2F, {axesConfigurations.axisPtQA, axesConfigurations.axisMult});
        histos.add(fmt::format("GeneratedWithPV/h{}_MidYVsMult_TwoPVsOrMore", Particlenames[i]).c_str(), "", kTH2F, {axesConfigurations.axisPtQA, axesConfigurations.axisMult});
      }
    }
    if (doprocessClosureTest) {
      if (pairLossK0Configurations.doClosureTestStages) {
        // Naming inside ClosureTest/PairLossK0: each folder is one reconstruction
        // requirement imposed on the same truth h-K0 pair. "Any" means the object
        // must have at least one reconstructed counterpart (a track with a
        // matching MC label for the trigger, a V0 candidate with a matching MC
        // core for the K0) in any reconstructed collision associated with this MC
        // collision, with no quality selection whatsoever.
        //   folder        trigger  K0      reconstruction requirement
        //   Truth         truth    truth   none (PairLossGenPair)
        //   AnyTrack      any      truth   trigger has a track in any collision
        //   AnyTrackK0    truth    any     K0 has a V0 candidate in any collision
        //   AnyTrackBoth  any      any     both requirements at once
        // The "any reconstructed collision" level exists only here: the
        // processPairLossK0MC stage ladder is evaluated in the best collision.
        //   Final         final    final   fully selected, both in one collision
        // "final" means the object has a reconstructed counterpart that survives
        // every selection the reconstructed correlation applies, and for the pair
        // both counterparts must live in the same reconstructed collision. Final/
        // is therefore the direct truth-coordinate counterpart of the ordinary
        // reconstructed correlation: Rec/Final isolates what is left once pair
        // loss is divided out -- duplicate reconstructed objects, fakes and bin
        // migration.
        // Every folder has the same three objects -- sameEvent/K0Short, hTrigger,
        // hK0Short -- and each of them is filled at the level its own folder
        // prescribes, so a folder can be normalised without looking at any other.
        // Consequence: the single-particle spectra repeat across folders in pairs
        // (Truth/hTrigger == AnyTrackK0/hTrigger, AnyTrack/hTrigger ==
        // AnyTrackBoth/hTrigger, Truth/hK0Short == AnyTrack/hK0Short,
        // AnyTrackK0/hK0Short == AnyTrackBoth/hK0Short). That is intended.
        histos.add("ClosureTest/PairLossK0/Truth/sameEvent/K0Short", "truth h-K0 pairs with the processPairLossK0MC selections", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
        histos.add("ClosureTest/PairLossK0/AnyTrack/sameEvent/K0Short", "truth h-K0 pairs whose truth trigger has a reconstructed-track match in any associated collision", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
        histos.add("ClosureTest/PairLossK0/AnyTrackK0/sameEvent/K0Short", "truth h-K0 pairs whose truth K0 has a V0-candidate match in any associated collision", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
        histos.add("ClosureTest/PairLossK0/AnyTrackBoth/sameEvent/K0Short", "truth h-K0 pairs with both the trigger track match and the K0 V0-candidate match in any associated collision", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
        histos.add("ClosureTest/PairLossK0/Truth/hTrigger", "truth triggers with the processPairLossK0MC selections;#it{p}_{T}^{truth} (GeV/#it{c});#eta^{truth};#varphi^{truth}", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi});
        histos.add("ClosureTest/PairLossK0/Truth/hK0Short", "truth K0s with the processPairLossK0MC selections;#it{p}_{T}^{truth} (GeV/#it{c});#eta^{truth};#varphi^{truth}", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi});
        histos.add("ClosureTest/PairLossK0/AnyTrack/hTrigger", "truth triggers with a reconstructed-track match in any associated collision;#it{p}_{T}^{truth} (GeV/#it{c});#eta^{truth};#varphi^{truth}", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi});
        histos.add("ClosureTest/PairLossK0/AnyTrack/hK0Short", "truth K0s; the K0 stays at truth level in this folder;#it{p}_{T}^{truth} (GeV/#it{c});#eta^{truth};#varphi^{truth}", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi});
        histos.add("ClosureTest/PairLossK0/AnyTrackK0/hTrigger", "truth triggers; the trigger stays at truth level in this folder;#it{p}_{T}^{truth} (GeV/#it{c});#eta^{truth};#varphi^{truth}", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi});
        histos.add("ClosureTest/PairLossK0/AnyTrackK0/hK0Short", "truth K0s with a V0-candidate match in any associated collision;#it{p}_{T}^{truth} (GeV/#it{c});#eta^{truth};#varphi^{truth}", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi});
        histos.add("ClosureTest/PairLossK0/AnyTrackBoth/hTrigger", "truth triggers with a reconstructed-track match in any associated collision;#it{p}_{T}^{truth} (GeV/#it{c});#eta^{truth};#varphi^{truth}", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi});
        histos.add("ClosureTest/PairLossK0/AnyTrackBoth/hK0Short", "truth K0s with a V0-candidate match in any associated collision;#it{p}_{T}^{truth} (GeV/#it{c});#eta^{truth};#varphi^{truth}", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi});
        histos.add("ClosureTest/PairLossK0/Final/sameEvent/K0Short", "truth h-K0 pairs whose trigger and K0 both have a fully selected reconstructed counterpart in the same reconstructed collision", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
        histos.add("ClosureTest/PairLossK0/Final/hTrigger", "truth triggers with a fully selected reconstructed counterpart;#it{p}_{T}^{truth} (GeV/#it{c});#eta^{truth};#varphi^{truth}", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi});
        histos.add("ClosureTest/PairLossK0/Final/hK0Short", "truth K0s with a fully selected reconstructed counterpart;#it{p}_{T}^{truth} (GeV/#it{c});#eta^{truth};#varphi^{truth}", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi});
      }
      for (int i = 0; i < AssocParticleTypes; i++) {
        if (TESTBIT(doCorrelation, i)) {
          histos.add(fmt::format("ClosureTest/sameEvent/{}", Particlenames[i]).c_str(), "", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
          if (masterConfigurations.doClosureTestTriggerWithRecoTrackMatch) {
            histos.add(fmt::format("ClosureTest/TriggerWithRecoTrackMatch/sameEvent/{}", Particlenames[i]).c_str(), "truth pairs whose trigger has at least one reconstructed track with a matching MC label", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
          }
          if (masterConfigurations.doCorrelationsHadronV0daughter) {
            histos.add(fmt::format("ClosureTest/sameEvent/{}_SameSignDaughter", Particlenames[i]).c_str(), "", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
            histos.add(fmt::format("ClosureTest/sameEvent/{}_OppSignDaughter", Particlenames[i]).c_str(), "", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
            if (masterConfigurations.doClosureTestTriggerWithRecoTrackMatch) {
              histos.add(fmt::format("ClosureTest/TriggerWithRecoTrackMatch/sameEvent/{}_SameSignDaughter", Particlenames[i]).c_str(), "truth trigger-daughter pairs whose trigger has a reconstructed-track match; same-sign daughter", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
              histos.add(fmt::format("ClosureTest/TriggerWithRecoTrackMatch/sameEvent/{}_OppSignDaughter", Particlenames[i]).c_str(), "truth trigger-daughter pairs whose trigger has a reconstructed-track match; opposite-sign daughter", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
            }
          }
        }
        if (TESTBIT(doCorrelation, i)) {
          histos.add(fmt::format("ClosureTest/h{}", Particlenames[i]).c_str(), "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi});
        }
      }
      histos.add("ClosureTest/hTrigger", "Trigger Tracks", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
      if (masterConfigurations.doClosureTestTriggerWithRecoTrackMatch) {
        histos.add("ClosureTest/TriggerWithRecoTrackMatch/hTrigger", "Truth triggers with at least one reconstructed track with a matching MC label;#it{p}_{T}^{truth} (GeV/c);#eta^{truth};centrality (%)", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
      }
    }
    if (doprocessFeedDown) {
      histos.add("hLambdaXiMinusFeeddownMatrix", "hLambdaXiMinusFeeddownMatrix", kTH2F, {axisPtLambda, axisPtCascade});
      histos.add("hLambdaXiZeroFeeddownMatrix", "hLambdaXiZeroFeeddownMatrix", kTH2F, {axisPtLambda, axisPtCascade});
      histos.add("hLambdaOmegaFeeddownMatrix", "hLambdaOmegaFeeddownMatrix", kTH2F, {axisPtLambda, axisPtCascade});
      histos.add("hAntiLambdaXiPlusFeeddownMatrix", "hAntiLambdaXiPlusFeeddownMatrix", kTH2F, {axisPtLambda, axisPtCascade});
      histos.add("hAntiLambdaXiZeroFeeddownMatrix", "hAntiLambdaXiZeroFeeddownMatrix", kTH2F, {axisPtLambda, axisPtCascade});
      histos.add("hAntiLambdaOmegaFeeddownMatrix", "hAntiLambdaOmegaFeeddownMatrix", kTH2F, {axisPtLambda, axisPtCascade});
      histos.add("hLambdaFromXiMinusEtaVsPtVsPhi", "hLambdaFromXiMinusEtaVsPtVsPhi", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi});
      histos.add("hLambdaFromXiZeroEtaVsPtVsPhi", "hLambdaFromXiZeroEtaVsPtVsPhi", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi});
      histos.add("hAntiLambdaFromXiPlusEtaVsPtVsPhi", "hAntiLambdaFromXiPlusEtaVsPtVsPhi", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi});
      histos.add("hAntiLambdaFromXiZeroEtaVsPtVsPhi", "hAntiLambdaFromXiZeroEtaVsPtVsPhi", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi});
      histos.add("GeneratedWithPV/hLambdaFromXiZero", "", kTH2F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta});
      histos.add("GeneratedWithPV/hLambdaFromXiMinus", "", kTH2F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta});
      histos.add("GeneratedWithPV/hAntiLambdaFromXiZero", "", kTH2F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta});
      histos.add("GeneratedWithPV/hAntiLambdaFromXiPlus", "", kTH2F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta});
    }
    if (doprocessPrediction) {
      mCounter.mPdgDatabase = pdgDB.service;
      mCounter.mSelectPrimaries = doAssocPhysicalPrimary.value;
      histos.add("Prediction/hEventSelection", "hEventSelection", kTH1F, {{3, 0, 3}});
      std::array<TString, 3> eventSelLabel = {"Read", "INELgt0", "|Z|<10"};
      for (int i = 1; i <= histos.get<TH1>(HIST("Prediction/hEventSelection"))->GetNbinsX(); i++) {
        histos.get<TH1>(HIST("Prediction/hEventSelection"))->GetXaxis()->SetBinLabel(i, eventSelLabel[i - 1]);
      }
      if (masterConfigurations.useCentralityinPrediction) {
        if (masterConfigurations.doSeparateFT0Prediction) {
          histos.add("Prediction/hTriggerFT0A", "Trigger Tracks FT0A", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
          histos.add("Prediction/hTriggerFT0C", "Trigger Tracks FT0C", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
        }
        histos.add("Prediction/hTrigger", "Trigger Tracks", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
      } else {
        if (masterConfigurations.doSeparateFT0Prediction) {
          histos.add("Prediction/hTriggerFT0A", "Trigger Tracks FT0A", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMultiplicity});
          histos.add("Prediction/hFT0AvsNchEta08", "Nch in 0.8 vs FT0A multiplicity", kTH2F, {axesConfigurations.axisMultiplicity, axesConfigurations.axisMidrapidityMultiplicity});
          histos.add("Prediction/hFT0AvsNchEta05", "Nch in 0.5 vs FT0A multiplicity", kTH2F, {axesConfigurations.axisMultiplicity, axesConfigurations.axisMidrapidityMultiplicity});
          histos.add("Prediction/hTriggerFT0C", "Trigger Tracks FT0C", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMultiplicity});
          histos.add("Prediction/hFT0CvsNchEta08", "Nch in 0.8 vs FT0C multiplicity", kTH2F, {axesConfigurations.axisMultiplicity, axesConfigurations.axisMidrapidityMultiplicity});
          histos.add("Prediction/hFT0CvsNchEta05", "Nch in 0.5 vs FT0C multiplicity", kTH2F, {axesConfigurations.axisMultiplicity, axesConfigurations.axisMidrapidityMultiplicity});
        }
        histos.add("Prediction/hTrigger", "Trigger Tracks", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMultiplicity});
        histos.add("Prediction/hFT0MvsNchEta08", "Nch in 0.8 vs FT0M multiplicity", kTH2F, {axesConfigurations.axisMultiplicity, axesConfigurations.axisMidrapidityMultiplicity});
        histos.add("Prediction/hFT0MvsNchEta05", "Nch in 0.5 vs FT0M multiplicity", kTH2F, {axesConfigurations.axisMultiplicity, axesConfigurations.axisMidrapidityMultiplicity});
      }
      for (int i = 0; i < AssocParticleTypes; i++) {
        if (TESTBIT(doCorrelation, i)) {
          histos.add(fmt::format("Prediction/h{}", Particlenames[i]).c_str(), "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi});
        }
        if (masterConfigurations.useCentralityinPrediction) {
          if (TESTBIT(doCorrelation, i)) {
            histos.add(fmt::format("Prediction/sameEvent/{}", Particlenames[i]).c_str(), "", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
          }
        } else {
          if (TESTBIT(doCorrelation, i)) {
            histos.add(fmt::format("Prediction/sameEvent/{}", Particlenames[i]).c_str(), "", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultiplicityNDim});
          }
        }
      }
      if (masterConfigurations.doSeparateFT0Prediction) {
        histos.addClone("Prediction/sameEvent/", "Prediction/sameEventFT0A/");
        histos.addClone("Prediction/sameEvent/", "Prediction/sameEventFT0C/");
      }
    }
    // visual inspection of sizes
    histos.print();

    // initialize CCDB *only* if efficiency correction requested
    // skip if not requested, saves a bit of time
    if (efficiencyFlags.applyEfficiencyCorrection || doprocessPairLossK0MC) {
      ccdb->setURL(ccdburl);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      ccdb->setFatalWhenNull(false);
    }
  }

  // this function allows for all event selections to be done in a modular way
  template <typename TCollision>
  bool isCollisionSelected(TCollision const& collision)
  {
    // ________________________________________________
    // Perform basic event selection
    if (!collision.sel8()) {
      return false;
    }
    if (!collision.selection_bit(aod::evsel::kIsTriggerTVX) && masterConfigurations.requireGoodTriggerTVX) {
      // FT0 vertex (acceptable FT0C-FT0A time difference) collisions
      return false;
    }
    if (std::abs(collision.posZ()) > masterConfigurations.zVertexCut) {
      return false;
    }
    if (collision.centFT0M() > axisRanges[5][1] || collision.centFT0M() < axisRanges[5][0]) {
      return false;
    }
    if (!collision.isInelGt0() && masterConfigurations.selectINELgtZERO) {
      return false;
    }
    if (!collision.isInelGt1() && masterConfigurations.selectINELgtONE) {
      return false;
    }
    if (!collision.selection_bit(aod::evsel::kIsGoodITSLayersAll) && masterConfigurations.requireAllGoodITSLayers) {
      return false;
    }
    if (!collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV) && masterConfigurations.requireGoodZvtxFT0vsPV) {
      // removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference
      // use this cut at low multiplicities with caution
      return false;
    }
    if (!collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup) && masterConfigurations.rejectSameBunchPileup) {
      // rejects collisions which are associated with the same "found-by-T0" bunch crossing
      return false;
    }
    if (zorroMask.value != "") {
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initZorro(bc);
      bool zorroSelected = zorro.isSelected(collision.template bc_as<aod::BCsWithTimestamps>().globalBC()); /// Just let Zorro do the accounting
      if (!zorroSelected) {
        return false;
      }
    }
    return true;
  }

  // event selections in Pb-Pb
  template <typename TCollision>
  bool isCollisionSelectedPbPb(TCollision const& collision, bool fillHists)
  {
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 0.5 /* all collisions */);
    }

    // Perform basic event selection
    if (!collision.sel8()) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 1.5 /* collisions  after sel8*/);
    }

    if (!collision.selection_bit(aod::evsel::kIsTriggerTVX) && masterConfigurations.requireGoodTriggerTVX) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 2.5 /* FT0 vertex (acceptable FT0C-FT0A time difference) collisions */);
    }

    if (std::abs(collision.posZ()) > masterConfigurations.zVertexCut) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 3.5 /* collisions  after sel pvz sel*/);
    }

    if (!collision.selection_bit(aod::evsel::kIsGoodITSLayersAll) && masterConfigurations.requireAllGoodITSLayers) {
      // cut time intervals with dead ITS staves
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 4.5 /* collisions  after cut time intervals with dead ITS staves*/);
    }

    if (!collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV) && masterConfigurations.requireGoodZvtxFT0vsPV) {
      // removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference
      // use this cut at low multiplicities with caution
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 5.5 /* removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference*/);
    }

    auto occupancy = collision.trackOccupancyInTimeRange();
    if (occupancy < cfgCutOccupancyLow || occupancy > cfgCutOccupancyHigh) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 6.5 /* Below min occupancy and Above max occupancy*/);
    }

    /*
    if (collision.alias_bit(kTVXinTRD)) {
      // TRD triggered
      return false;
    }
    */

    if (!collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      // reject collisions close to Time Frame borders
      // O2-4623
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 7.5 /* reject collisions close to Time Frame borders*/);
    }

    if (!collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      // reject events affected by the ITS ROF border
      // O2-4309
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 8.5 /* reject events affected by the ITS ROF border*/);
    }

    if (!collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      // rejects collisions which are associated with the same "found-by-T0" bunch crossing
      // https://indico.cern.ch/event/1396220/#1-event-selection-with-its-rof
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 9.5 /* rejects collisions which are associated with the same "found-by-T0" bunch crossing*/);
    }
    return true;
  }

  double getMagneticField(uint64_t timestamp)
  {
    static parameters::GRPMagField* grpo = nullptr;
    if (grpo == nullptr) {
      grpo = ccdb->getForTimeStamp<parameters::GRPMagField>("GLO/Config/GRPMagField", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 0;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
    }

    return 0.1 * (grpo->getNominalL3Field()); // 1 T = 10 kG
  }

  double getPairLossMagneticField(int runNumber, uint64_t timestamp)
  {
    if (mPairLossRunNumber == runNumber) {
      return mPairLossMagneticField;
    }
    auto* grpo = ccdb->getForTimeStamp<parameters::GRPMagField>("GLO/Config/GRPMagField", timestamp);
    if (grpo == nullptr) {
      LOGF(fatal, "GRP object not found for K0 pair-loss diagnostic at timestamp %llu", timestamp);
      return 0.0;
    }
    mPairLossRunNumber = runNumber;
    mPairLossMagneticField = 0.1 * grpo->getNominalL3Field(); // 1 T = 10 kG
    LOGF(info, "K0 pair-loss diagnostic loaded %.1f T for run %d", mPairLossMagneticField, runNumber);
    return mPairLossMagneticField;
  }
  // if this process function is enabled, it will be such that only events with trigger particles within a given
  // trigger pt bin are taken for the entire processing. This allows for the calculation of e.g. efficiencies
  // within an event class that has a trigger (which may differ with respect to other cases, to be checked)

  // for map determining which trigger bins are present and which aren't
  std::vector<uint32_t> triggerPresenceMap;

  void processSelectEventWithTrigger(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults> const& collisions,
                                     aod::TriggerTracks const& triggerTracks, TracksComplete const&, aod::BCsWithTimestamps const&)
  {
    // setup
    triggerPresenceMap.clear();
    triggerPresenceMap.resize(collisions.size(), 0);

    for (auto const& collision : collisions) {
      // ________________________________________________
      // Perform basic event selection
      if (!isCollisionSelected(collision)) {
        continue;
      }

      // do not forget to re-group ...
      auto slicedTriggerTracks = triggerTracks.sliceBy(collisionSliceTracks, collision.globalIndex());

      for (auto const& triggerTrack : slicedTriggerTracks) {
        auto track = triggerTrack.track_as<TracksComplete>();
        if (!isValidTrigger(track, triggerTrack.isLeading())) {
          continue;
        }
        if (trackSelection.checkForITSTPCMissmatchMC && bitcheck(triggerTrack.mcMask(), 13)) {
          continue;
        }
        auto binNumber = histos.get<TH1>(HIST("axes/hPtTriggerAxis"))->FindFixBin(track.pt()) - 1;
        SETBIT(triggerPresenceMap[collision.globalIndex()], binNumber);
      }
    }
  }

  void processSameEventHHadrons(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults>::iterator const& collision,
                                aod::AssocHadrons const& assocHadrons, aod::TriggerTracks const& triggerTracks,
                                TracksComplete const&, aod::BCsWithTimestamps const&)
  {
    BinningTypePP colBinning{{axesConfigurations.axisVtxZ, axesConfigurations.axisMult}, true}; // true is for 'ignore overflows' (true by default). Underflows and overflows will have bin -1.

    // ________________________________________________
    // skip if desired trigger not found
    if (triggerPresenceMap.size() > 0 && !TESTBIT(triggerPresenceMap[collision.globalIndex()], triggerBinToSelect)) {
      return;
    }

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    auto bField = getMagneticField(bc.timestamp());

    if (efficiencyFlags.applyEfficiencyCorrection) {
      initEfficiencyFromCCDB(bc);
    }

    // ________________________________________________
    // Perform basic event selection
    if (!isCollisionSelected(collision)) {
      return;
    }
    // ________________________________________________
    if (!doprocessSameEventHCascades && !doprocessSameEventHV0s && !doprocessSameEventHPions && doMixingQAandEventQA) {
      histos.fill(HIST("MixingQA/hSECollisionBins"), colBinning.getBin({collision.posZ(), collision.centFT0M()}));
      histos.fill(HIST("EventQA/hMult"), collision.centFT0M());
      histos.fill(HIST("EventQA/hPvz"), collision.posZ());
    }

    // Do basic QA
    if (!doprocessSameEventHCascades && !doprocessSameEventHV0s && !doprocessSameEventHPions) {
      for (auto const& triggerTrack : triggerTracks) {
        auto track = triggerTrack.track_as<TracksComplete>();
        if (!isValidTrigger(track, triggerTrack.isLeading())) {
          continue;
        }
        if (trackSelection.checkForITSTPCMissmatchMC && bitcheck(triggerTrack.mcMask(), 13)) {
          continue;
        }
        histos.fill(HIST("hDCAzTriggerHadron"), track.dcaZ(), track.pt());
        histos.fill(HIST("hDCAxyTriggerHadron"), track.dcaXY(), track.pt());
        float efficiency = 1.0f;
        if (efficiencyFlags.applyEfficiencyCorrection) {
          efficiency = hEfficiencyTrigger->Interpolate(track.pt(), track.eta());
        }
        if (efficiency == 0) { // check for zero efficiency, do not apply if the case
          efficiency = 1;
        }
        float weight = efficiencyFlags.applyEfficiencyCorrection ? 1. / efficiency : 1.0f;
        histos.fill(HIST("hTriggerAllSelectedEtaVsPt"), track.pt(), track.eta(), collision.centFT0M());
        histos.fill(HIST("hTriggerPtResolution"), track.pt(), triggerTrack.mcOriginalPt());
        if (masterConfigurations.doTriggPhysicalPrimary && !triggerTrack.mcPhysicalPrimary()) {
          continue;
        }
        histos.fill(HIST("hTriggerPrimaryEtaVsPt"), track.pt(), track.eta(), collision.centFT0M());
        histos.fill(HIST("hTrackEtaVsPtVsPhi"), track.pt(), track.eta(), track.phi(), weight);
      }
    }
    for (auto const& assocTrack : assocHadrons) {
      auto assoc = assocTrack.track_as<TracksComplete>();
      if (!isValidAssocHadron(assoc)) {
        continue;
      }
      if (trackSelection.checkForITSTPCMissmatchMC && bitcheck(assocTrack.mcMask(), 13)) {
        continue;
      }
      float efficiency = 1.0f;
      float purity = 1.0f;
      histos.fill(HIST("hDCAzAssociatedHadron"), assoc.dcaZ(), assoc.pt());
      histos.fill(HIST("hDCAxyAssociatedHadron"), assoc.dcaXY(), assoc.pt());
      if (efficiencyFlags.applyEfficiencyCorrection) {
        efficiency = hEfficiencyHadron->Interpolate(assoc.pt(), assoc.eta());
        if (efficiencyFlags.applyPurityHadron) {
          purity = hPurityHadron->Interpolate(assoc.pt());
        }
      }
      if (efficiency == 0) { // check for zero efficiency, do not apply if the case
        efficiency = 1;
      }
      float weight = efficiencyFlags.applyEfficiencyCorrection ? purity / efficiency : 1.0f;
      histos.fill(HIST("hAssocHadronsAllSelectedEtaVsPt"), assoc.pt(), assoc.eta(), collision.centFT0M(), weight);
      histos.fill(HIST("hAssocPtResolution"), assoc.pt(), assocTrack.mcOriginalPt());
      if (doAssocPhysicalPrimary && !assocTrack.mcPhysicalPrimary()) {
        continue;
      }
      histos.fill(HIST("hAssocPrimaryEtaVsPt"), assoc.pt(), assoc.eta(), collision.centFT0M());
      histos.fill(HIST("hAsssocTrackEtaVsPtVsPhi"), assoc.pt(), assoc.eta(), assoc.phi(), weight);
    }

    // ________________________________________________
    // Do hadron - hadron correlations
    if (masterConfigurations.doFullCorrelationStudy) {
      fillCorrelationsHadron(triggerTracks, assocHadrons, false, collision.posZ(), collision.centFT0M(), bField);
    }
  }

  void runSameEventHV0sCore(aod::AssocV0s const& associatedV0s, aod::TriggerTracks const& triggerTracks, TracksComplete const& tracks,
                            float pvx, float pvy, float pvz, float cent, float multNTracksPVeta1, double bField)
  {
    std::variant<BinningTypePP, BinningTypePbPb> colBinning =
      masterConfigurations.doPPAnalysis
        ? std::variant<BinningTypePP, BinningTypePbPb>{
            BinningTypePP{{axesConfigurations.axisVtxZ, axesConfigurations.axisMult}, true}}
        : std::variant<BinningTypePP, BinningTypePbPb>{BinningTypePbPb{{axesConfigurations.axisVtxZ, axesConfigurations.axisMult}, true}};

    if (!doprocessSameEventHCascades && doMixingQAandEventQA) {
      std::visit([&](auto const& binning) {
        histos.fill(HIST("MixingQA/hSECollisionBins"), binning.getBin({pvz, cent}));
      },
                 colBinning);
      histos.fill(HIST("EventQA/hMult"), cent);
      histos.fill(HIST("EventQA/hPvz"), pvz);
      histos.fill(HIST("EventQA/hMultFT0vsTPC"), cent, multNTracksPVeta1);
    }

    // Do basic QA
    std::array<TH2F*, AssocV0Types> hEfficiencyV0{nullptr, nullptr, nullptr};
    std::array<THnF*, AssocV0Types> hEfficiencyV0MultVsPhi{nullptr, nullptr, nullptr};
    if (efficiencyFlags.applyEffAsFunctionOfMultAndPhi) {
      hEfficiencyV0MultVsPhi[0] = hEfficiencyK0ShortMultVsPhi;
      hEfficiencyV0MultVsPhi[1] = hEfficiencyLambdaMultVsPhi;
      hEfficiencyV0MultVsPhi[2] = hEfficiencyAntiLambdaMultVsPhi;
    } else {
      hEfficiencyV0[0] = hEfficiencyK0Short;
      hEfficiencyV0[1] = hEfficiencyLambda;
      hEfficiencyV0[2] = hEfficiencyAntiLambda;
    }

    for (auto const& v0 : associatedV0s) {
      auto v0Data = v0.v0Core_as<V0DatasWithoutTrackX>();

      //---] track quality check [---
      auto postrack = v0Data.posTrack_as<TracksComplete>();
      auto negtrack = v0Data.negTrack_as<TracksComplete>();
      if (postrack.tpcNClsCrossedRows() < trackSelection.minTPCNCrossedRowsAssociated || negtrack.tpcNClsCrossedRows() < trackSelection.minTPCNCrossedRowsAssociated) {
        continue;
      }
      if (trackSelection.checksRequireTPCChi2 && (postrack.tpcChi2NCl() < trackSelection.minTPCChi2PerClusterAssociated || negtrack.tpcChi2NCl() < trackSelection.minTPCChi2PerClusterAssociated)) {
        continue;
      }
      if (trackSelection.requireClusterInITS && (postrack.itsNCls() < trackSelection.minITSClustersForDaughterTracks || negtrack.itsNCls() < trackSelection.minITSClustersForDaughterTracks)) {
        continue;
      }
      //---] syst cuts [---
      if (masterConfigurations.doPPAnalysis && (v0Data.v0radius() < v0Selection.v0RadiusMin || v0Data.v0radius() > v0Selection.v0RadiusMax ||
                                                std::abs(v0Data.dcapostopv()) < v0Selection.dcapostopv || std::abs(v0Data.dcanegtopv()) < v0Selection.dcanegtopv ||
                                                v0Data.v0cosPA() < v0Selection.v0cospa || v0Data.dcaV0daughters() > v0Selection.dcaV0dau)) {
        continue;
      }
      if (!masterConfigurations.doPPAnalysis && !v0SelectedPbPb(v0Data)) {
        continue;
      }
      uint64_t selMap = v0selectionBitmap(v0Data, pvx, pvy, pvz);
      const int localDensityV0 = masterConfigurations.doLocalDensityStudy ? computeLocalDensity(v0Data.eta(), v0Data.phi(), tracks, {v0Data.posTrackId(), v0Data.negTrackId()}) : 0;

      static_for<0, 2>([&](auto i) {
        constexpr int Index = i.value;
        float efficiency = 1.0f;
        if (efficiencyFlags.applyEffAsFunctionOfMultAndPhi) {
          if (efficiencyFlags.applyEfficiencyCorrection) {
            std::array<double, 4> bin = {v0Data.pt(), v0Data.eta(), v0Data.phi(), cent};
            efficiency = hEfficiencyV0MultVsPhi[Index]->GetBinContent(hEfficiencyV0MultVsPhi[Index]->GetBin(bin.data()));
          }
        } else {
          if (efficiencyFlags.applyEfficiencyCorrection) {
            efficiency = hEfficiencyV0[Index]->Interpolate(v0Data.pt(), v0Data.eta());
          }
        }
        if (efficiency == 0) { // check for zero efficiency, do not apply if the case
          efficiency = 1;
        }
        float weight = efficiencyFlags.applyEfficiencyCorrection ? 1. / efficiency : 1.0f;
        if (v0.compatible(Index, trackSelection.dEdxCompatibility) && (!masterConfigurations.doMCassociation || v0.mcTrue(Index)) && (!doAssocPhysicalPrimary || v0.mcPhysicalPrimary()) && (!efficiencyFlags.applyEfficiencyCorrection || efficiency != 0)) {
          if ((TESTBIT(doCorrelation, Index)) && (masterConfigurations.doPPAnalysis || (TESTBIT(selMap, Index) && TESTBIT(selMap, Index + 3)))) {
            histos.fill(HIST("h3d") + HIST(V0names[Index]) + HIST("Spectrum"), v0Data.pt(), cent, v0.invMassNSigma(Index), weight);
            if (masterConfigurations.doMassSpectrumCheck) {
              histos.fill(HIST("h") + HIST(V0names[Index]) + HIST("PtVsMass"), v0Data.pt(), getV0InvariantMass<Index>(v0Data));
            }
            if (std::abs(v0Data.rapidity(Index)) < ySel) {
              histos.fill(HIST("h3d") + HIST(V0names[Index]) + HIST("SpectrumY"), v0Data.pt(), cent, v0.invMassNSigma(Index), weight);
            }
            if (!efficiencyFlags.applyEffAsFunctionOfMultAndPhi) {
              if ((-massWindowConfigurations.maxBgNSigma < v0.invMassNSigma(Index) && v0.invMassNSigma(Index) < -massWindowConfigurations.minBgNSigma) || (+massWindowConfigurations.minBgNSigma < v0.invMassNSigma(Index) && v0.invMassNSigma(Index) < +massWindowConfigurations.maxBgNSigma)) {
                histos.fill(HIST("h") + HIST(V0names[Index]) + HIST("EtaVsPtVsPhiBg"), v0Data.pt(), v0Data.eta(), v0Data.phi(), weight);
              }
              if (-massWindowConfigurations.maxPeakNSigma < v0.invMassNSigma(Index) && v0.invMassNSigma(Index) < +massWindowConfigurations.maxPeakNSigma) {
                histos.fill(HIST("h") + HIST(V0names[Index]) + HIST("EtaVsPtVsPhi"), v0Data.pt(), v0Data.eta(), v0Data.phi(), weight);
              }
            } else {
              if ((-massWindowConfigurations.maxBgNSigma < v0.invMassNSigma(Index) && v0.invMassNSigma(Index) < -massWindowConfigurations.minBgNSigma) || (+massWindowConfigurations.minBgNSigma < v0.invMassNSigma(Index) && v0.invMassNSigma(Index) < +massWindowConfigurations.maxBgNSigma)) {
                histos.fill(HIST("h") + HIST(V0names[Index]) + HIST("EtaVsPtVsPhiVsCentBg"), v0Data.pt(), v0Data.eta(), v0Data.phi(), cent, weight);
              }
              if (-massWindowConfigurations.maxPeakNSigma < v0.invMassNSigma(Index) && v0.invMassNSigma(Index) < +massWindowConfigurations.maxPeakNSigma) {
                histos.fill(HIST("h") + HIST(V0names[Index]) + HIST("EtaVsPtVsPhiVsCent"), v0Data.pt(), v0Data.eta(), v0Data.phi(), cent, weight);
              }
            }
            if (masterConfigurations.doLocalDensityStudy) {
              if ((-massWindowConfigurations.maxBgNSigma < v0.invMassNSigma(Index) && v0.invMassNSigma(Index) < -massWindowConfigurations.minBgNSigma) || (+massWindowConfigurations.minBgNSigma < v0.invMassNSigma(Index) && v0.invMassNSigma(Index) < +massWindowConfigurations.maxBgNSigma)) {
                histos.fill(HIST("h") + HIST(V0names[Index]) + HIST("EtaVsPtVsLocalDensityBg"), v0Data.pt(), v0Data.eta(), localDensityV0, weight);
              }
              if (-massWindowConfigurations.maxPeakNSigma < v0.invMassNSigma(Index) && v0.invMassNSigma(Index) < +massWindowConfigurations.maxPeakNSigma) {
                histos.fill(HIST("h") + HIST(V0names[Index]) + HIST("EtaVsPtVsLocalDensity"), v0Data.pt(), v0Data.eta(), localDensityV0, weight);
              }
            }
          }
        }
      });
    }
    if (!doprocessSameEventHCascades) {
      for (auto const& triggerTrack : triggerTracks) {
        auto track = triggerTrack.track_as<TracksComplete>();
        if (!isValidTrigger(track, triggerTrack.isLeading())) {
          continue;
        }
        if (trackSelection.checkForITSTPCMissmatchMC && bitcheck(triggerTrack.mcMask(), 13)) {
          continue;
        }
        histos.fill(HIST("hTriggerAllSelectedEtaVsPt"), track.pt(), track.eta(), cent);
        histos.fill(HIST("hTriggerPtResolution"), track.pt(), triggerTrack.mcOriginalPt());
        if (masterConfigurations.doTriggPhysicalPrimary && !triggerTrack.mcPhysicalPrimary()) {
          continue;
        }
        if (!efficiencyFlags.applyEffAsFunctionOfMultAndPhi) {
          histos.fill(HIST("hTriggerPrimaryEtaVsPt"), track.pt(), track.eta(), cent);
        } else {
          histos.fill(HIST("hTriggerPrimaryEtaVsPt"), track.pt(), track.eta(), track.phi(), cent);
        }
        histos.fill(HIST("hTrackEtaVsPtVsPhi"), track.pt(), track.eta(), track.phi());
        if (masterConfigurations.doLocalDensityStudy) {
          histos.fill(HIST("hTrackEtaVsPtVsLocalDensity"), track.pt(), track.eta(), computeLocalDensity(track.eta(), track.phi(), tracks, {track.globalIndex()}));
        }
      }
    }

    // ________________________________________________
    // Do hadron - V0 correlations
    if (masterConfigurations.doFullCorrelationStudy) {
      fillCorrelationsV0(triggerTracks, associatedV0s, false, false, pvx, pvy, pvz, cent, bField);
    }
  }

  template <typename TCollision>
  void runSameEventHV0s(TCollision const& collision, aod::AssocV0s const& associatedV0s, aod::TriggerTracks const& triggerTracks, TracksComplete const& tracks)
  {
    const float cent = masterConfigurations.doPPAnalysis ? collision.centFT0M() : collision.centFT0C();

    // Keep the trigger-presence and reconstructed-event decisions identical
    // for the ordinary Rec process and for its PairLoss replica.
    if (triggerPresenceMap.size() > 0 && !TESTBIT(triggerPresenceMap[collision.globalIndex()], triggerBinToSelect)) {
      return;
    }
    if ((masterConfigurations.doPPAnalysis && !isCollisionSelected(collision)) ||
        (!masterConfigurations.doPPAnalysis && !isCollisionSelectedPbPb(collision, true))) {
      return;
    }

    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    const auto bField = getMagneticField(bc.timestamp());
    if (efficiencyFlags.applyEfficiencyCorrection) {
      initEfficiencyFromCCDB(bc);
    }
    runSameEventHV0sCore(associatedV0s, triggerTracks, tracks,
                         collision.posX(), collision.posY(), collision.posZ(), cent, collision.multNTracksPVeta1(), bField);
  }

  void processSameEventHV0s(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::PVMults>::iterator const& collision,
                            aod::AssocV0s const& associatedV0s, aod::TriggerTracks const& triggerTracks,
                            V0DatasWithoutTrackX const&, TracksComplete const& tracks, aod::BCsWithTimestamps const&)
  {
    runSameEventHV0s(collision, associatedV0s, triggerTracks, tracks);
  }

  void processSameEventHCascades(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::PVMults>::iterator const& collision,
                                 aod::AssocV0s const&, aod::AssocCascades const& associatedCascades, aod::TriggerTracks const& triggerTracks,
                                 V0DatasWithoutTrackX const&, aod::CascDatas const&, TracksComplete const&, aod::BCsWithTimestamps const&)
  {
    std::variant<BinningTypePP, BinningTypePbPb> colBinning =
      masterConfigurations.doPPAnalysis
        ? std::variant<BinningTypePP, BinningTypePbPb>{
            BinningTypePP{{axesConfigurations.axisVtxZ, axesConfigurations.axisMult}, true}}
        : std::variant<BinningTypePP, BinningTypePbPb>{BinningTypePbPb{{axesConfigurations.axisVtxZ, axesConfigurations.axisMult}, true}};

    double cent = masterConfigurations.doPPAnalysis ? collision.centFT0M() : collision.centFT0C();
    // ________________________________________________
    // skip if desired trigger not found
    if (triggerPresenceMap.size() > 0 && !TESTBIT(triggerPresenceMap[collision.globalIndex()], triggerBinToSelect)) {
      return;
    }

    // ________________________________________________
    // Perform basic event selection
    if (((masterConfigurations.doPPAnalysis && !isCollisionSelected(collision))) || (!masterConfigurations.doPPAnalysis && !isCollisionSelectedPbPb(collision, true))) {
      return;
    }
    // ________________________________________________
    if (doMixingQAandEventQA) {
      std::visit([&](auto const& binning) {
        histos.fill(HIST("MixingQA/hSECollisionBins"), binning.getBin({collision.posZ(), cent}));
      },
                 colBinning);
      histos.fill(HIST("EventQA/hMult"), cent);
      histos.fill(HIST("EventQA/hPvz"), collision.posZ());
    }
    // Do basic QA
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    auto bField = getMagneticField(bc.timestamp());
    if (efficiencyFlags.applyEfficiencyCorrection) {
      initEfficiencyFromCCDB(bc);
    }
    std::array<TH2F*, AssocCascadeTypes> hEfficiencyCascade{nullptr, nullptr, nullptr, nullptr};
    std::array<THnF*, AssocCascadeTypes> hEfficiencyCascadeMultVsPhi{nullptr, nullptr, nullptr, nullptr};
    if (efficiencyFlags.applyEfficiencyCorrection) {
      if (efficiencyFlags.applyEffAsFunctionOfMultAndPhi) {
        hEfficiencyCascadeMultVsPhi[0] = hEfficiencyXiMinusMultVsPhi;
        hEfficiencyCascadeMultVsPhi[1] = hEfficiencyXiPlusMultVsPhi;
        hEfficiencyCascadeMultVsPhi[2] = hEfficiencyOmegaMinusMultVsPhi;
        hEfficiencyCascadeMultVsPhi[3] = hEfficiencyOmegaPlusMultVsPhi;
      } else {
        hEfficiencyCascade[0] = hEfficiencyXiMinus;
        hEfficiencyCascade[1] = hEfficiencyXiPlus;
        hEfficiencyCascade[2] = hEfficiencyOmegaMinus;
        hEfficiencyCascade[3] = hEfficiencyOmegaPlus;
      }
    }
    for (auto const& casc : associatedCascades) {
      auto cascData = casc.cascData();

      //---] syst cuts [---
      if (masterConfigurations.doPPAnalysis && (std::abs(cascData.dcapostopv()) < v0Selection.dcapostopv ||
                                                std::abs(cascData.dcanegtopv()) < v0Selection.dcanegtopv ||
                                                std::abs(cascData.dcabachtopv()) < cascadeSelections.dcaBachToPV ||
                                                cascData.dcaV0daughters() > cascadeSelections.cascdcaV0dau ||
                                                cascData.dcacascdaughters() > cascadeSelections.dcaCascDaughters ||
                                                cascData.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < cascadeSelections.cascv0cospa ||
                                                cascData.casccosPA(collision.posX(), collision.posY(), collision.posZ()) < cascadeSelections.cascCospa ||
                                                cascData.cascradius() < cascadeSelections.cascRadius ||
                                                std::abs(cascData.dcav0topv(collision.posX(), collision.posY(), collision.posZ())) < cascadeSelections.cascdcaV0ToPV ||
                                                std::abs(cascData.mLambda() - o2::constants::physics::MassLambda0) > cascadeSelections.cascV0masswindow)) {
        continue;
      }
      if (!masterConfigurations.doPPAnalysis && !cascadeSelectedPbPb(cascData, collision.posX(), collision.posY(), collision.posZ())) {
        continue;
      }
      uint64_t cascselMap = cascadeselectionBitmap(cascData, collision.posX(), collision.posY(), collision.posZ());
      //---] track quality check [---
      auto postrack = cascData.posTrack_as<TracksComplete>();
      auto negtrack = cascData.negTrack_as<TracksComplete>();
      auto bachtrack = cascData.bachelor_as<TracksComplete>();
      if (postrack.tpcNClsCrossedRows() < trackSelection.minTPCNCrossedRowsAssociated || negtrack.tpcNClsCrossedRows() < trackSelection.minTPCNCrossedRowsAssociated || bachtrack.tpcNClsCrossedRows() < trackSelection.minTPCNCrossedRowsAssociated) {
        continue;
      }
      if (trackSelection.checksRequireTPCChi2 && (postrack.tpcChi2NCl() < trackSelection.minTPCChi2PerClusterAssociated || negtrack.tpcChi2NCl() < trackSelection.minTPCChi2PerClusterAssociated || bachtrack.tpcChi2NCl() < trackSelection.minTPCChi2PerClusterAssociated)) {
        continue;
      }
      if (trackSelection.requireClusterInITS && (postrack.itsNCls() < trackSelection.minITSClustersForDaughterTracks || negtrack.itsNCls() < trackSelection.minITSClustersForDaughterTracks)) {
        continue;
      }

      static_for<0, 3>([&](auto i) {
        constexpr int Index = i.value;
        if ((Index == IndexOmegaMinus || Index == IndexOmegaPlus) && casc.compatible(Index, trackSelection.dEdxCompatibility) && std::abs(casc.invMassNSigma(Index - 2)) < massWindowConfigurations.nSigmaNearXiMassCenter) {
          return;
        }
        float efficiency = 1.0f;
        if (efficiencyFlags.applyEfficiencyCorrection) {
          if (efficiencyFlags.applyEffAsFunctionOfMultAndPhi) {
            std::array<double, 4> bin = {cascData.pt(), cascData.eta(), cascData.phi(), cent};
            efficiency = hEfficiencyCascadeMultVsPhi[Index]->GetBinContent(hEfficiencyCascadeMultVsPhi[Index]->GetBin(bin.data()));
          } else {
            efficiency = hEfficiencyCascade[Index]->Interpolate(cascData.pt(), cascData.eta());
          }
        }
        if (efficiency == 0) { // check for zero efficiency, do not apply if the case
          efficiency = 1;
        }
        float weight = efficiencyFlags.applyEfficiencyCorrection ? 1. / efficiency : 1.0f;
        if (casc.compatible(Index, trackSelection.dEdxCompatibility) && (!masterConfigurations.doMCassociation || casc.mcTrue(Index)) && (!doAssocPhysicalPrimary || casc.mcPhysicalPrimary()) && (!efficiencyFlags.applyEfficiencyCorrection || efficiency != 0)) {
          if (TESTBIT(doCorrelation, Index + 3) && (masterConfigurations.doPPAnalysis || (TESTBIT(cascselMap, Index) && TESTBIT(cascselMap, Index + 4) && TESTBIT(cascselMap, Index + 8) && TESTBIT(cascselMap, Index + 12)))) {
            histos.fill(HIST("h3d") + HIST(Cascadenames[Index]) + HIST("Spectrum"), cascData.pt(), cent, casc.invMassNSigma(Index), weight);
            if (std::abs(cascData.rapidity(Index)) < ySel) {
              histos.fill(HIST("h3d") + HIST(Cascadenames[Index]) + HIST("SpectrumY"), cascData.pt(), cent, casc.invMassNSigma(Index), weight);
            }
            if (!efficiencyFlags.applyEffAsFunctionOfMultAndPhi) {
              if (-massWindowConfigurations.maxPeakNSigma < casc.invMassNSigma(Index) && casc.invMassNSigma(Index) < +massWindowConfigurations.maxPeakNSigma) {
                histos.fill(HIST("h") + HIST(Cascadenames[Index]) + HIST("EtaVsPtVsPhi"), cascData.pt(), cascData.eta(), cascData.phi(), weight);
              }
              if ((-massWindowConfigurations.maxBgNSigma < casc.invMassNSigma(Index) && casc.invMassNSigma(Index) < -massWindowConfigurations.minBgNSigma) || (+massWindowConfigurations.minBgNSigma < casc.invMassNSigma(Index) && casc.invMassNSigma(Index) < +massWindowConfigurations.maxBgNSigma)) {
                histos.fill(HIST("h") + HIST(Cascadenames[Index]) + HIST("EtaVsPtVsPhiBg"), cascData.pt(), cascData.eta(), cascData.phi(), weight);
              }
            } else {
              if (-massWindowConfigurations.maxPeakNSigma < casc.invMassNSigma(Index) && casc.invMassNSigma(Index) < +massWindowConfigurations.maxPeakNSigma) {
                histos.fill(HIST("h") + HIST(Cascadenames[Index]) + HIST("EtaVsPtVsPhiVsCent"), cascData.pt(), cascData.eta(), cascData.phi(), cent, weight);
              }
              if ((-massWindowConfigurations.maxBgNSigma < casc.invMassNSigma(Index) && casc.invMassNSigma(Index) < -massWindowConfigurations.minBgNSigma) || (+massWindowConfigurations.minBgNSigma < casc.invMassNSigma(Index) && casc.invMassNSigma(Index) < +massWindowConfigurations.maxBgNSigma)) {
                histos.fill(HIST("h") + HIST(Cascadenames[Index]) + HIST("EtaVsPtVsPhiVsCentBg"), cascData.pt(), cascData.eta(), cascData.phi(), cent, weight);
              }
            }
          }
        }
      });
    }
    for (auto const& triggerTrack : triggerTracks) {
      auto track = triggerTrack.track_as<TracksComplete>();
      if (!isValidTrigger(track, triggerTrack.isLeading())) {
        continue;
      }
      if (trackSelection.checkForITSTPCMissmatchMC && bitcheck(triggerTrack.mcMask(), 13)) {
        continue;
      }
      histos.fill(HIST("hTriggerAllSelectedEtaVsPt"), track.pt(), track.eta(), cent);
      histos.fill(HIST("hTriggerPtResolution"), track.pt(), triggerTrack.mcOriginalPt());
      if (masterConfigurations.doTriggPhysicalPrimary && !triggerTrack.mcPhysicalPrimary()) {
        continue;
      }
      if (!efficiencyFlags.applyEffAsFunctionOfMultAndPhi) {
        histos.fill(HIST("hTriggerPrimaryEtaVsPt"), track.pt(), track.eta(), cent);
      } else {
        histos.fill(HIST("hTriggerPrimaryEtaVsPt"), track.pt(), track.eta(), track.phi(), cent);
      }
      if (track.sign() > 0) {
        histos.fill(HIST("hPositiveTriggerPrimaryEtaVsPt"), track.pt(), track.eta(), cent);
      } else {
        histos.fill(HIST("hNegativeTriggerPrimaryEtaVsPt"), track.pt(), track.eta(), cent);
      }
      histos.fill(HIST("hTrackEtaVsPtVsPhi"), track.pt(), track.eta(), track.phi());
    }

    // ________________________________________________
    // Do hadron - cascade correlations
    if (masterConfigurations.doFullCorrelationStudy) {
      fillCorrelationsCascade(triggerTracks, associatedCascades, false, false, collision.posX(), collision.posY(), collision.posZ(), cent, bField);
    }
  }
  void processSameEventHPions(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults>::iterator const& collision,
                              soa::Join<aod::AssocHadrons, aod::AssocPID> const& associatedPions, soa::Join<aod::TriggerTracks, aod::TriggerTrackExtras> const& triggerTracks,
                              TracksComplete const&, aod::BCsWithTimestamps const&)
  {
    BinningTypePP colBinning{{axesConfigurations.axisVtxZ, axesConfigurations.axisMult}, true};
    // ________________________________________________
    // skip if desired trigger not found
    if (triggerPresenceMap.size() > 0 && !TESTBIT(triggerPresenceMap[collision.globalIndex()], triggerBinToSelect)) {
      return;
    }
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    auto bField = getMagneticField(bc.timestamp());

    if (efficiencyFlags.applyEfficiencyCorrection) {
      initEfficiencyFromCCDB(bc);
    }

    // ________________________________________________
    // Perform basic event selection
    if (!isCollisionSelected(collision)) {
      return;
    }
    // ________________________________________________
    if (!doprocessSameEventHCascades && !doprocessSameEventHV0s && doMixingQAandEventQA) {
      histos.fill(HIST("MixingQA/hSECollisionBins"), colBinning.getBin({collision.posZ(), collision.centFT0M()}));
      histos.fill(HIST("EventQA/hMult"), collision.centFT0M());
      histos.fill(HIST("EventQA/hPvz"), collision.posZ());
    }
    // Do basic QA
    for (auto const& pion : associatedPions) {
      auto pionTrack = pion.track_as<TracksComplete>();
      if (!isValidAssocHadron(pionTrack)) {
        continue;
      }
      if (trackSelection.checkForITSTPCMissmatchMC && bitcheck(pion.mcMask(), 13)) {
        continue;
      }

      histos.fill(HIST("hPionEtaVsPtAllSelected"), pionTrack.pt(), pionTrack.eta(), collision.centFT0M());
      if (doAssocPhysicalPrimary && !pion.mcPhysicalPrimary()) {
        continue;
      }
      if (masterConfigurations.doMCassociation && std::abs(pion.pdgCode()) != PdgCodes[IndexPion]) {
        continue;
      }
      histos.fill(HIST("hPionEtaVsPt"), pionTrack.pt(), pionTrack.eta(), collision.centFT0M());
      if (pionTrack.sign() > 0) {
        histos.fill(HIST("hPositivePionEtaVsPt"), pionTrack.pt(), pionTrack.eta(), collision.centFT0M());
      } else {
        histos.fill(HIST("hNegativePionEtaVsPt"), pionTrack.pt(), pionTrack.eta(), collision.centFT0M());
      }
    }
    if (!doprocessSameEventHCascades && !doprocessSameEventHV0s) {
      for (auto const& triggerTrack : triggerTracks) {
        auto track = triggerTrack.track_as<TracksComplete>();
        if (!isValidTrigger(track, triggerTrack.isLeading())) {
          continue;
        }
        if (trackSelection.checkForITSTPCMissmatchMC && bitcheck(triggerTrack.mcMask(), 13)) {
          continue;
        }
        histos.fill(HIST("hTriggerAllSelectedEtaVsPt"), track.pt(), track.eta(), collision.centFT0M());
        histos.fill(HIST("hTriggerPtResolution"), track.pt(), triggerTrack.mcOriginalPt());
        if (masterConfigurations.doTriggPhysicalPrimary && !triggerTrack.mcPhysicalPrimary()) {
          continue;
        }
        histos.fill(HIST("hTriggerPrimaryEtaVsPt"), track.pt(), track.eta(), collision.centFT0M());
        histos.fill(HIST("hTrackEtaVsPtVsPhi"), track.pt(), track.eta(), track.phi());
      }
    }

    // ________________________________________________
    // Do hadron - Pion correlations
    if (masterConfigurations.doFullCorrelationStudy) {
      fillCorrelationsHadron(triggerTracks, associatedPions, false, collision.posZ(), collision.centFT0M(), bField);
    }
  }

  void processMixedEventHHadrons(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults> const& collisions,
                                 aod::AssocHadrons const& assocHadrons, aod::TriggerTracks const& triggerTracks,
                                 TracksComplete const&, aod::BCsWithTimestamps const&)
  {
    BinningTypePP colBinning{{axesConfigurations.axisVtxZ, axesConfigurations.axisMult}, true};
    for (auto const& [collision1, collision2] : soa::selfCombinations(colBinning, masterConfigurations.mixingParameter, -1, collisions, collisions)) {
      auto bc = collision1.bc_as<aod::BCsWithTimestamps>();
      auto bField = getMagneticField(bc.timestamp());
      // ________________________________________________
      if (efficiencyFlags.applyEfficiencyCorrection) {
        initEfficiencyFromCCDB(bc);
      }
      // ________________________________________________
      // skip if desired trigger not found
      if (triggerPresenceMap.size() > 0 && (!TESTBIT(triggerPresenceMap[collision1.globalIndex()], triggerBinToSelect) || !TESTBIT(triggerPresenceMap[collision2.globalIndex()], triggerBinToSelect))) {
        return;
      }

      // ________________________________________________
      // Perform basic event selection on both collisions
      if (!isCollisionSelected(collision1) || !isCollisionSelected(collision2)) {
        continue;
      }
      if (collision1.centFT0M() > axisRanges[5][1] || collision1.centFT0M() < axisRanges[5][0]) {
        continue;
      }
      if (collision2.centFT0M() > axisRanges[5][1] || collision2.centFT0M() < axisRanges[5][0]) {
        continue;
      }
      if (doMixingQAandEventQA) {
        if (collision1.globalIndex() == collision2.globalIndex()) {
          histos.fill(HIST("MixingQA/hMixingQA"), 0.0f); // same-collision pair counting
        }
        histos.fill(HIST("MixingQA/hMEpvz1"), collision1.posZ());
        histos.fill(HIST("MixingQA/hMEpvz2"), collision2.posZ());
        histos.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), collision1.centFT0M()}));
      }
      // ________________________________________________
      // Do slicing
      auto slicedTriggerTracks = triggerTracks.sliceBy(collisionSliceTracks, collision1.globalIndex());
      auto slicedAssocHadrons = assocHadrons.sliceBy(collisionSliceHadrons, collision2.globalIndex());
      // ________________________________________________
      // Do hadron - hadron correlations
      if (masterConfigurations.doFullCorrelationStudy) {
        fillCorrelationsHadron(slicedTriggerTracks, slicedAssocHadrons, true, collision1.posZ(), collision1.centFT0M(), bField);
      }
    }
  }

  void processMixedEventHV0s(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::PVMults> const& collisions,
                             aod::AssocV0s const& associatedV0s, aod::TriggerTracks const& triggerTracks,
                             V0DatasWithoutTrackX const&, TracksComplete const&, aod::BCsWithTimestamps const&)
  {
    std::variant<BinningTypePP, BinningTypePbPb> colBinning =
      masterConfigurations.doPPAnalysis
        ? std::variant<BinningTypePP, BinningTypePbPb>{
            BinningTypePP{{axesConfigurations.axisVtxZ, axesConfigurations.axisMult}, true}}
        : std::variant<BinningTypePP, BinningTypePbPb>{BinningTypePbPb{{axesConfigurations.axisVtxZ, axesConfigurations.axisMult}, true}};

    std::visit([&](auto const& binning) {
      for (auto const& [collision1, collision2] : soa::selfCombinations(binning, masterConfigurations.mixingParameter, -1, collisions, collisions)) {
        double cent1 = masterConfigurations.doPPAnalysis ? collision1.centFT0M() : collision1.centFT0C();
        double cent2 = masterConfigurations.doPPAnalysis ? collision2.centFT0M() : collision2.centFT0C();
        auto bc = collision1.template bc_as<aod::BCsWithTimestamps>();
        auto bField = getMagneticField(bc.timestamp());
        // ________________________________________________
        if (efficiencyFlags.applyEfficiencyCorrection) {
          initEfficiencyFromCCDB(bc);
        }
        // ________________________________________________
        // skip if desired trigger not found
        if (triggerPresenceMap.size() > 0 && (!TESTBIT(triggerPresenceMap[collision1.globalIndex()], triggerBinToSelect) || !TESTBIT(triggerPresenceMap[collision2.globalIndex()], triggerBinToSelect))) {
          continue;
        }

        // Perform basic event selection on both collisions
        if ((masterConfigurations.doPPAnalysis && (!isCollisionSelected(collision1) || !isCollisionSelected(collision2))) || (!masterConfigurations.doPPAnalysis && (!isCollisionSelectedPbPb(collision1, false) || (!isCollisionSelectedPbPb(collision2, false))))) {
          continue;
        }
        if (cent1 > axisRanges[5][1] || cent1 < axisRanges[5][0]) {
          continue;
        }
        if (cent2 > axisRanges[5][1] || cent2 < axisRanges[5][0]) {
          continue;
        }

        if (!doprocessMixedEventHCascades && doMixingQAandEventQA) {
          if (collision1.globalIndex() == collision2.globalIndex()) {
            histos.fill(HIST("MixingQA/hMixingQA"), 0.0f); // same-collision pair counting
          }
          histos.fill(HIST("MixingQA/hMEpvz1"), collision1.posZ());
          histos.fill(HIST("MixingQA/hMEpvz2"), collision2.posZ());
          histos.fill(HIST("MixingQA/hMECollisionBins"), binning.getBin({collision1.posZ(), cent1}));
        }
        // ________________________________________________
        // Do slicing
        auto slicedTriggerTracks = triggerTracks.sliceBy(collisionSliceTracks, collision1.globalIndex());
        auto slicedAssocV0s = associatedV0s.sliceBy(collisionSliceV0s, collision2.globalIndex());
        // ________________________________________________
        // Do hadron - V0 correlations
        if (masterConfigurations.doFullCorrelationStudy) {
          fillCorrelationsV0(slicedTriggerTracks, slicedAssocV0s, true, false, collision1.posX(), collision1.posY(), collision1.posZ(), cent1, bField);
        }
      }
    },
               colBinning);
  }

  void processMixedEventHCascades(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::PVMults> const& collisions,
                                  aod::AssocV0s const&, aod::AssocCascades const& associatedCascades, aod::TriggerTracks const& triggerTracks,
                                  V0DatasWithoutTrackX const&, aod::CascDatas const&, TracksComplete const&, aod::BCsWithTimestamps const&)
  {
    std::variant<BinningTypePP, BinningTypePbPb> colBinning =
      masterConfigurations.doPPAnalysis
        ? std::variant<BinningTypePP, BinningTypePbPb>{
            BinningTypePP{{axesConfigurations.axisVtxZ, axesConfigurations.axisMult}, true}}
        : std::variant<BinningTypePP, BinningTypePbPb>{BinningTypePbPb{{axesConfigurations.axisVtxZ, axesConfigurations.axisMult}, true}};

    std::visit([&](auto const& binning) {
      for (auto const& [collision1, collision2] : soa::selfCombinations(binning, masterConfigurations.mixingParameter, -1, collisions, collisions)) {
        double cent1 = masterConfigurations.doPPAnalysis ? collision1.centFT0M() : collision1.centFT0C();
        double cent2 = masterConfigurations.doPPAnalysis ? collision2.centFT0M() : collision2.centFT0C();
        // ________________________________________________
        auto bc = collision1.template bc_as<aod::BCsWithTimestamps>();
        auto bField = getMagneticField(bc.timestamp());
        if (efficiencyFlags.applyEfficiencyCorrection) {
          initEfficiencyFromCCDB(bc);
        }
        // ________________________________________________
        // skip if desired trigger not found
        if (triggerPresenceMap.size() > 0 && (!TESTBIT(triggerPresenceMap[collision1.globalIndex()], triggerBinToSelect) || !TESTBIT(triggerPresenceMap[collision2.globalIndex()], triggerBinToSelect))) {
          continue;
        }

        // Perform basic event selection on both collisions
        if ((masterConfigurations.doPPAnalysis && (!isCollisionSelected(collision1) || !isCollisionSelected(collision2))) || (!masterConfigurations.doPPAnalysis && (!isCollisionSelectedPbPb(collision1, false) || (!isCollisionSelectedPbPb(collision2, false))))) {
          continue;
        }
        if (cent1 > axisRanges[5][1] || cent1 < axisRanges[5][0]) {
          continue;
        }
        if (cent2 > axisRanges[5][1] || cent2 < axisRanges[5][0]) {
          continue;
        }
        if (doMixingQAandEventQA) {
          if (collision1.globalIndex() == collision2.globalIndex()) {
            histos.fill(HIST("MixingQA/hMixingQA"), 0.0f); // same-collision pair counting
          }
          histos.fill(HIST("MixingQA/hMEpvz1"), collision1.posZ());
          histos.fill(HIST("MixingQA/hMEpvz2"), collision2.posZ());
          histos.fill(HIST("MixingQA/hMECollisionBins"), binning.getBin({collision1.posZ(), cent1}));
        }
        // ________________________________________________
        // Do slicing
        auto slicedTriggerTracks = triggerTracks.sliceBy(collisionSliceTracks, collision1.globalIndex());
        auto slicedAssocCascades = associatedCascades.sliceBy(collisionSliceCascades, collision2.globalIndex());
        // ________________________________________________
        // Do hadron - cascade correlations
        if (masterConfigurations.doFullCorrelationStudy) {
          fillCorrelationsCascade(slicedTriggerTracks, slicedAssocCascades, true, false, collision1.posX(), collision1.posY(), collision1.posZ(), cent1, bField);
        }
      }
    },
               colBinning);
  }

  void processMixedEventHPions(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults> const& collisions,
                               soa::Join<aod::AssocHadrons, aod::AssocPID> const& assocPions, soa::Join<aod::TriggerTracks, aod::TriggerTrackExtras> const& triggerTracks,
                               TracksComplete const&, aod::BCsWithTimestamps const&)
  {
    BinningTypePP colBinning{{axesConfigurations.axisVtxZ, axesConfigurations.axisMult}, true};
    for (auto const& [collision1, collision2] : soa::selfCombinations(colBinning, masterConfigurations.mixingParameter, -1, collisions, collisions)) {
      auto bc = collision1.bc_as<aod::BCsWithTimestamps>();
      auto bField = getMagneticField(bc.timestamp());
      // ________________________________________________
      if (efficiencyFlags.applyEfficiencyCorrection) {
        initEfficiencyFromCCDB(bc);
      }
      // ________________________________________________
      // skip if desired trigger not found
      if (triggerPresenceMap.size() > 0 && (!TESTBIT(triggerPresenceMap[collision1.globalIndex()], triggerBinToSelect) || !TESTBIT(triggerPresenceMap[collision2.globalIndex()], triggerBinToSelect))) {
        continue;
      }

      // ________________________________________________
      // Perform basic event selection on both collisions
      if (!isCollisionSelected(collision1) || !isCollisionSelected(collision2)) {
        continue;
      }
      if (collision1.centFT0M() > axisRanges[5][1] || collision1.centFT0M() < axisRanges[5][0]) {
        continue;
      }
      if (collision2.centFT0M() > axisRanges[5][1] || collision2.centFT0M() < axisRanges[5][0]) {
        continue;
      }
      if (doMixingQAandEventQA) {
        if (collision1.globalIndex() == collision2.globalIndex()) {
          histos.fill(HIST("MixingQA/hMixingQA"), 0.0f); // same-collision pair counting
        }
        histos.fill(HIST("MixingQA/hMEpvz1"), collision1.posZ());
        histos.fill(HIST("MixingQA/hMEpvz2"), collision2.posZ());
        histos.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), collision1.centFT0M()}));
      }
      // ________________________________________________
      // Do slicing
      auto slicedTriggerTracks = triggerTracks.sliceBy(collisionSliceTracks, collision1.globalIndex());
      auto slicedAssocPions = assocPions.sliceBy(collisionSliceHadrons, collision2.globalIndex());
      // ________________________________________________
      // Do hadron - cascade correlations
      if (masterConfigurations.doFullCorrelationStudy) {
        fillCorrelationsHadron(slicedTriggerTracks, slicedAssocPions, true, collision1.posZ(), collision1.centFT0M(), bField);
      }
    }
  }

  void processMCGenerated(aod::McCollision const& /*mcCollision*/, soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::PVMults>> const& collisions, aod::McParticles const& mcParticles, TracksCompleteMC const& tracks)
  {
    histos.fill(HIST("hClosureTestEventCounter"), 2.5f);

    for (auto const& mcParticle : mcParticles) {
      double geta = mcParticle.eta();
      if (std::abs(geta) > etaSel) {
        continue;
      }
      double gpt = mcParticle.pt();
      if (std::abs(mcParticle.pdgCode()) == PDG_t::kPiPlus || std::abs(mcParticle.pdgCode()) == PDG_t::kKPlus || std::abs(mcParticle.pdgCode()) == PDG_t::kProton || std::abs(mcParticle.pdgCode()) == PDG_t::kElectron || std::abs(mcParticle.pdgCode()) == PDG_t::kMuonMinus) {
        if (!masterConfigurations.doTriggPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
          histos.fill(HIST("hGeneratedQAPtTrigger"), gpt, 0.0f); // step 1: before all selections
        }
      }

      if (!doAssocPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
        if (mcParticle.pdgCode() == PDG_t::kK0Short && masterConfigurations.doCorrelationK0Short) {
          histos.fill(HIST("hGeneratedQAPtAssociatedK0"), gpt, 0.0f); // step 1: before all selections
        }
      }
    }

    for (auto const& mcParticle : mcParticles) {
      if (doAssocPhysicalPrimaryInGen && !mcParticle.isPhysicalPrimary()) {
        continue;
      }
      static_for<0, 7>([&](auto i) {
        constexpr int Index = i.value;
        if (i == IndexPion && mcParticle.pdgCode() > Neutral) {
          histos.fill(HIST("Generated/hPositive") + HIST(Particlenames[Index]), mcParticle.pt(), mcParticle.eta(), 1);
        } else if (i == IndexPion && mcParticle.pdgCode() < Neutral) {
          histos.fill(HIST("Generated/hNegative") + HIST(Particlenames[Index]), mcParticle.pt(), mcParticle.eta(), 1);
        } else if (mcParticle.pdgCode() == PdgCodes[i]) {
          if (efficiencyFlags.applyEffAsFunctionOfMultAndPhi) {
            histos.fill(HIST("Generated/h") + HIST(Particlenames[Index]), mcParticle.pt(), mcParticle.eta(), mcParticle.phi(), 1);
          } else {
            histos.fill(HIST("Generated/h") + HIST(Particlenames[Index]), mcParticle.pt(), mcParticle.eta(), 1);
          }
        }
      });
    }
    if (collisions.size() < 1) {
      return;
    }

    // determine best collision properties
    int biggestNContribs = -1;
    int64_t bestCollisionId = -1;
    float bestCollisionFT0Mpercentile = -1;
    float bestCollisionFT0Cpercentile = -1;
    float bestCollisionVtxZ = 0.0f;
    bool bestCollisionSel8 = false;
    bool bestCollisionINELgtZERO = false;
    bool bestCollisionINELgtONE = false;
    bool bestCollisionNoSameBunchPileup = false;
    bool bestCollisionGoodTriggerTVX = false;
    bool bestCollisionGoodZvtxFT0vsPV = false;
    bool isCollisionSelect = false;
    uint32_t bestCollisionTriggerPresenceMap = 0;

    for (auto const& collision : collisions) {
      if (biggestNContribs < collision.numContrib()) {
        biggestNContribs = collision.numContrib();
        bestCollisionId = collision.globalIndex();
        bestCollisionFT0Mpercentile = collision.centFT0M();
        bestCollisionFT0Cpercentile = collision.centFT0C();
        if (masterConfigurations.applyNewMCSelection) {
          isCollisionSelect = ((masterConfigurations.doPPAnalysis && isCollisionSelected(collision)) || (!masterConfigurations.doPPAnalysis && isCollisionSelectedPbPb(collision, false)));
        } else {
          bestCollisionSel8 = collision.sel8();
          bestCollisionVtxZ = collision.posZ();
          bestCollisionINELgtZERO = collision.isInelGt0();
          bestCollisionINELgtONE = collision.isInelGt1();
          bestCollisionNoSameBunchPileup = collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup);
          bestCollisionGoodTriggerTVX = collision.selection_bit(aod::evsel::kIsTriggerTVX);
          bestCollisionGoodZvtxFT0vsPV = collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV);
        }
        if (triggerPresenceMap.size() > 0) {
          bestCollisionTriggerPresenceMap = triggerPresenceMap[collision.globalIndex()];
        }
      }
    }

    if (collisions.size() > 1) {
      for (auto const& mcParticle : mcParticles) {
        if (doAssocPhysicalPrimaryInGen && !mcParticle.isPhysicalPrimary()) {
          continue;
        }
        if (std::abs(mcParticle.y()) > ySel) {
          continue;
        }
        static_for<0, 7>([&](auto i) {
          constexpr int Index = i.value;

          if (mcParticle.pdgCode() == PdgCodes[i]) {
            histos.fill(HIST("GeneratedWithPV/h") + HIST(Particlenames[Index]) + HIST("_MidYVsMult_TwoPVsOrMore"), mcParticle.pt(), bestCollisionFT0Mpercentile);
          }
        });
      }
    }

    // do selections on best collision
    // WARNING: if 2 PV case large, this will not necessarily be fine!
    //          caution advised!

    // ________________________________________________
    // skip if desired trigger not found
    if (triggerPresenceMap.size() > 0 && !TESTBIT(bestCollisionTriggerPresenceMap, triggerBinToSelect)) {
      return;
    }
    if (masterConfigurations.applyNewMCSelection) {
      if (!isCollisionSelect) {
        return;
      }
    } else {
      if (!bestCollisionSel8) {
        return;
      }
      if (std::abs(bestCollisionVtxZ) > masterConfigurations.zVertexCut) {
        return;
      }
      if (!bestCollisionINELgtZERO) {
        return;
      }
      if (masterConfigurations.selectINELgtONE && !bestCollisionINELgtONE) {
        return;
      }
      if (masterConfigurations.rejectSameBunchPileup && !bestCollisionNoSameBunchPileup) {
        return;
      }
      if (masterConfigurations.requireGoodTriggerTVX && !bestCollisionGoodTriggerTVX) {
        return;
      }
      if (masterConfigurations.requireGoodZvtxFT0vsPV && !bestCollisionGoodZvtxFT0vsPV) {
        return;
      }
    }

    histos.fill(HIST("hClosureTestEventCounter"), 3.5f);

    for (auto const& mcParticle : mcParticles) {
      double geta = mcParticle.eta();
      if (std::abs(geta) > etaSel) {
        continue;
      }
      double gpt = mcParticle.pt();
      if (std::abs(mcParticle.pdgCode()) == PDG_t::kPiPlus || std::abs(mcParticle.pdgCode()) == PDG_t::kKPlus || std::abs(mcParticle.pdgCode()) == PDG_t::kProton || std::abs(mcParticle.pdgCode()) == PDG_t::kElectron || std::abs(mcParticle.pdgCode()) == PDG_t::kMuonMinus) {
        if (!masterConfigurations.doTriggPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
          histos.fill(HIST("hGeneratedQAPtTrigger"), gpt, 1.0f); // step 2: after event selection
        }
      }

      if (!doAssocPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
        if (mcParticle.pdgCode() == PDG_t::kK0Short && masterConfigurations.doCorrelationK0Short) {
          histos.fill(HIST("hGeneratedQAPtAssociatedK0"), gpt, 1.0f); // step 2: before all selections
        }
      }
    }

    // The local density of a generated particle is counted from the tracks of the best collision,
    // the same collision whose centrality already labels these generated histograms.
    const auto bestCollisionTracks = tracks.sliceBy(pairLossTracksPerCollision, bestCollisionId);

    for (auto const& mcParticle : mcParticles) {
      if (doAssocPhysicalPrimaryInGen && !mcParticle.isPhysicalPrimary()) {
        continue;
      }
      double geta = mcParticle.eta();
      double gpt = mcParticle.pt();
      if (std::abs(mcParticle.pdgCode()) == PDG_t::kPiPlus || std::abs(mcParticle.pdgCode()) == PDG_t::kKPlus || std::abs(mcParticle.pdgCode()) == PDG_t::kProton || std::abs(mcParticle.pdgCode()) == PDG_t::kElectron || std::abs(mcParticle.pdgCode()) == PDG_t::kMuonMinus) {
        if (efficiencyFlags.applyEffAsFunctionOfMultAndPhi) {
          histos.fill(HIST("GeneratedWithPV/hTrigger"), gpt, geta, mcParticle.phi(), bestCollisionFT0Cpercentile);
        } else {
          histos.fill(HIST("GeneratedWithPV/hTrigger"), gpt, geta, bestCollisionFT0Mpercentile);
        }
        if (masterConfigurations.doLocalDensityStudy && masterConfigurations.doPPAnalysis) {
          histos.fill(HIST("GeneratedWithPV/hTriggerLocalDensity"), gpt, geta, computeLocalDensityGen(mcParticle, bestCollisionTracks));
        }
        if (mcParticle.pdgCode() > 0) {
          histos.fill(HIST("GeneratedWithPV/hPositiveTrigger"), gpt, geta, bestCollisionFT0Mpercentile);
        } else {
          histos.fill(HIST("GeneratedWithPV/hNegativeTrigger"), gpt, geta, bestCollisionFT0Mpercentile);
        }
      }

      if (mcParticle.pdgCode() == PDG_t::kLambda0 && !doAssocPhysicalPrimaryInGen && !mcParticle.isPhysicalPrimary()) {
        if (std::abs(geta) > etaSel) {
          continue;
        }
        auto lamMothers = mcParticle.mothers_as<aod::McParticles>();
        if (lamMothers.size() == 1) {
          for (const auto& lamParticleMother : lamMothers) {
            if (std::abs(lamParticleMother.eta()) > etaSel) {
              continue;
            }
            if (doprocessFeedDown && lamParticleMother.pdgCode() == PDG_t::kXiMinus) // Xi Minus Mother Matched
            {
              histos.fill(HIST("GeneratedWithPV/hLambdaFromXiMinus"), gpt, geta);
            }
            if (doprocessFeedDown && lamParticleMother.pdgCode() == o2::constants::physics::Pdg::kXi0) // Xi Zero Mother Matched
            {
              histos.fill(HIST("GeneratedWithPV/hLambdaFromXiZero"), gpt, geta);
            }
          }
        }
      }
      if (mcParticle.pdgCode() == PDG_t::kLambda0Bar && !doAssocPhysicalPrimaryInGen && !mcParticle.isPhysicalPrimary()) {
        if (std::abs(geta) > etaSel) {
          continue;
        }
        auto lamMothers = mcParticle.mothers_as<aod::McParticles>();
        if (lamMothers.size() == 1) {
          for (const auto& lamParticleMother : lamMothers) {
            if (std::abs(lamParticleMother.eta()) > etaSel) {
              continue;
            }
            if (doprocessFeedDown && lamParticleMother.pdgCode() == PDG_t::kXiPlusBar) {
              histos.fill(HIST("GeneratedWithPV/hAntiLambdaFromXiPlus"), gpt, geta);
            }
            if (doprocessFeedDown && lamParticleMother.pdgCode() == -o2::constants::physics::Pdg::kXi0) // Xi Zero Mother Matched
            {
              histos.fill(HIST("GeneratedWithPV/hAntiLambdaFromXiZero"), gpt, geta);
            }
          }
        }
      }
      static_for<0, 7>([&](auto i) {
        constexpr int Index = i.value;
        if (i == IndexPion && mcParticle.pdgCode() == PdgCodes[i] && mcParticle.pdgCode() > Neutral) {
          histos.fill(HIST("GeneratedWithPV/hPositive") + HIST(Particlenames[Index]), mcParticle.pt(), mcParticle.eta(), bestCollisionFT0Mpercentile);
        } else if (i == IndexPion && mcParticle.pdgCode() == PdgCodes[i] && mcParticle.pdgCode() < Neutral) {
          histos.fill(HIST("GeneratedWithPV/hNegative") + HIST(Particlenames[Index]), mcParticle.pt(), mcParticle.eta(), bestCollisionFT0Mpercentile);
        }

        if (mcParticle.pdgCode() == PdgCodes[i]) {
          if (efficiencyFlags.applyEffAsFunctionOfMultAndPhi) {
            histos.fill(HIST("GeneratedWithPV/h") + HIST(Particlenames[Index]), gpt, geta, mcParticle.phi(), bestCollisionFT0Cpercentile);
          } else {
            histos.fill(HIST("GeneratedWithPV/h") + HIST(Particlenames[Index]), gpt, geta, bestCollisionFT0Mpercentile);
          }
          if (std::abs(mcParticle.y()) < ySel) {
            histos.fill(HIST("GeneratedWithPV/h") + HIST(Particlenames[Index]) + HIST("_MidYVsMult"), gpt, bestCollisionFT0Mpercentile);
          }
          if (masterConfigurations.doLocalDensityStudy && masterConfigurations.doPPAnalysis) {
            histos.fill(HIST("GeneratedWithPV/h") + HIST(Particlenames[Index]) + HIST("LocalDensity"), gpt, geta, computeLocalDensityGen(mcParticle, bestCollisionTracks));
          }
        }
      });
    }
  }
  void processPairLossK0MC(aod::McCollision const& mcCollision,
                           soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::PVMults>> const& recCollisions,
                           aod::McParticles const& mcParticles,
                           aod::McTrackLabels const& trackLabels,
                           aod::V0MCCores const& v0MCCores,
                           aod::TriggerTracks const& triggerTracks,
                           aod::AssocV0s const& associatedV0s,
                           aod::BCsWithTimestamps const&,
                           V0DatasWithoutTrackX const& v0Candidates,
                           TracksComplete const& tracks)
  {
    // Part 3: a self-contained generator-level study. It deliberately touches no
    // reconstructed quantity in its event selection or in any of its axes: the
    // event is selected on generated INEL>0 and the generated vertex only, the
    // multiplicity is counted from generated particles, and every object is
    // filled with generated coordinates.
    //
    // Reconstruction enters in exactly one place -- whether a generated object
    // has a reconstructed counterpart at all -- and that splits the very same
    // generated sample into Reconstructed/ and NotReconstructed/. Because all
    // three folders share generated coordinates, Gen == Reconstructed +
    // NotReconstructed bin by bin, so NotReconstructed/Gen reads directly as
    // "in which pT, eta, phi and multiplicity region do generated objects fail
    // to be reconstructed".
    //
    // "Reconstructed" is the loosest possible statement, with no quality
    // selection of any kind: for a trigger, some track in some reconstructed
    // collision of this MC collision carries its MC label; for a K0, some V0
    // candidate carries its MC core. Objects belonging to an MC collision that
    // produced no reconstructed collision at all therefore land in
    // NotReconstructed/ too; hEventCounter and hNRecoCollisions are there so
    // that contribution can be separated out afterwards.
    auto runGenLevelStudy = [&]() {
      histos.fill(HIST("PairLossK0/GenStudy/hEventCounter"), 0.0f);

      // Generated-level event selection. No reconstructed variable is used.
      // INEL>1 implies INEL>0, so only the tighter enabled selection has to be evaluated
      if (masterConfigurations.selectINELgtONE) {
        if (!o2::pwglf::isINELgt1mc(mcParticles, pdgDB)) {
          return;
        }
      } else if (masterConfigurations.selectINELgtZERO) {
        if (!o2::pwglf::isINELgt0mc(mcParticles, pdgDB)) {
          return;
        }
      }
      histos.fill(HIST("PairLossK0/GenStudy/hEventCounter"), 1.0f);
      if (std::abs(mcCollision.posZ()) > masterConfigurations.zVertexCut) {
        return;
      }
      histos.fill(HIST("PairLossK0/GenStudy/hEventCounter"), 2.0f);
      if (recCollisions.size() > 0) {
        histos.fill(HIST("PairLossK0/GenStudy/hEventCounter"), 3.0f);
      }
      histos.fill(HIST("PairLossK0/GenStudy/hNRecoCollisions"), recCollisions.size());

      // Multiplicity of this MC collision: generated charged physical primaries
      // within |eta| < 0.8. Primaries always, independent of every analysis
      // configurable, so that the multiplicity axis keeps one fixed meaning.
      //
      // Deliberately a local counter rather than a member: the shared mCounter
      // only gets its PDG database wired up when processPrediction runs, and its
      // mSelectPrimaries follows doAssocPhysicalPrimary. Adding a second counter
      // as a task member is not an option either -- the struct is already at the
      // member limit that Framework/StructToTuple.h can destructure. The object
      // is a bool and a pointer, so building it per MC collision costs nothing.
      o2::pwglf::ParticleCounter<o2::framework::O2DatabasePDG> genStudyCounter;
      genStudyCounter.mPdgDatabase = pdgDB.service;
      genStudyCounter.mSelectPrimaries = true;
      const float generatedNch = genStudyCounter.countEta08(mcParticles);
      histos.fill(HIST("PairLossK0/GenStudy/hNch"), generatedNch);

      // Reconstructed-object bookkeeping. The framework has already grouped
      // recCollisions by this MC collision, so dereferencing a reconstructed
      // collision back to its MC collision needs no extra work here, and the
      // generated event selection above is by construction identical for all of
      // them.
      std::unordered_set<int64_t> reconstructedTrackMcIds;
      std::unordered_set<int64_t> reconstructedV0McIds;
      std::vector<GenStudyPairObject> genStudyTriggers;
      std::vector<GenStudyPairObject> genStudyK0s;
      for (auto const& collision : recCollisions) {
        const auto trackSlice = tracks.sliceBy(pairLossTracksPerCollision, collision.globalIndex());
        for (auto const& track : trackSlice) {
          const auto trackLabel = trackLabels.iteratorAt(track.globalIndex());
          if (trackLabel.has_mcParticle()) {
            reconstructedTrackMcIds.insert(trackLabel.mcParticleId());
          }
        }
        const auto v0Slice = v0Candidates.sliceBy(pairLossV0sPerCollision, collision.globalIndex());
        for (auto const& v0 : v0Slice) {
          const auto v0MC = v0MCCores.iteratorAt(v0.globalIndex());
          if (v0MC.particleIdMC() < 0 || v0MC.pdgCode() != PDG_t::kK0Short) {
            continue;
          }
          reconstructedV0McIds.insert(v0MC.particleIdMC());
        }
      }

      for (auto const& mcParticle : mcParticles) {
        const float genPt = mcParticle.pt();
        const float genEta = mcParticle.eta();
        const float genPhi = mcParticle.phi();
        if (std::abs(genEta) > etaSel) {
          continue;
        }

        if (isPairLossTriggerPdg(mcParticle.pdgCode()) &&
            genPt >= axisRanges[3][0] && genPt <= axisRanges[3][1] &&
            (!masterConfigurations.doTriggPhysicalPrimary || mcParticle.isPhysicalPrimary())) {
          // Same charge requirement the other two parts apply, so that the three
          // trigger definitions stay comparable.
          auto const* pdgParticle = pdgDB->GetParticle(mcParticle.pdgCode());
          const double charge = pdgParticle != nullptr ? pdgParticle->Charge() : 0.0;
          const int sign = charge > 0.0 ? 1 : (charge < 0.0 ? -1 : 0);
          if (!((triggerParticleCharge > 0 && sign < 0) || (triggerParticleCharge < 0 && sign > 0) || sign == 0)) {
            histos.fill(HIST("PairLossK0/GenStudy/Gen/hTrigger"), genPt, genEta, genPhi, generatedNch);
            if (reconstructedTrackMcIds.count(mcParticle.globalIndex()) > 0) {
              histos.fill(HIST("PairLossK0/GenStudy/Reconstructed/hTrigger"), genPt, genEta, genPhi, generatedNch);
            } else {
              histos.fill(HIST("PairLossK0/GenStudy/NotReconstructed/hTrigger"), genPt, genEta, genPhi, generatedNch);
            }
            genStudyTriggers.push_back(GenStudyPairObject{
              .pt = genPt,
              .eta = genEta,
              .phi = genPhi,
              .globalIndex = static_cast<int64_t>(mcParticle.globalIndex()),
              .motherIndex = mcParticle.has_mothers() ? static_cast<int64_t>(mcParticle.mothers_first_as<aod::McParticles>().globalIndex()) : -1,
              .reconstructed = reconstructedTrackMcIds.count(mcParticle.globalIndex()) > 0});
          }
        }

        if (mcParticle.pdgCode() == PDG_t::kK0Short &&
            genPt >= axisRanges[2][0] && genPt <= axisRanges[2][1] &&
            (!doAssocPhysicalPrimary || mcParticle.isPhysicalPrimary())) {
          // Same findability definition the stage ladder uses: decayed to
          // pi+ pi- with both charged daughters inside the tracking acceptance.
          // A K0 that is not findable could never have been reconstructed, so
          // splitting on it is what makes NotReconstructed/ interpretable --
          // without it the folder is dominated by decays whose daughters simply
          // left the acceptance.
          bool hasPositiveDaughter = false;
          bool hasNegativeDaughter = false;
          bool daughtersInAcceptance = true;
          for (auto const& daughter : mcParticle.daughters_as<aod::McParticles>()) {
            if (daughter.pdgCode() != PDG_t::kPiPlus && daughter.pdgCode() != -PDG_t::kPiPlus) {
              continue;
            }
            if (daughter.pdgCode() == PDG_t::kPiPlus) {
              hasPositiveDaughter = true;
            } else {
              hasNegativeDaughter = true;
            }
            if (daughter.pt() < pairLossK0Configurations.daughterPtMin ||
                std::abs(daughter.eta()) > pairLossK0Configurations.daughterEtaMax) {
              daughtersInAcceptance = false;
            }
          }
          const float k0Findable = (hasPositiveDaughter && hasNegativeDaughter && daughtersInAcceptance) ? 1.0f : 0.0f;

          histos.fill(HIST("PairLossK0/GenStudy/Gen/hK0Short"), genPt, genEta, genPhi, generatedNch, k0Findable);
          if (reconstructedV0McIds.count(mcParticle.globalIndex()) > 0) {
            histos.fill(HIST("PairLossK0/GenStudy/Reconstructed/hK0Short"), genPt, genEta, genPhi, generatedNch, k0Findable);
          } else {
            histos.fill(HIST("PairLossK0/GenStudy/NotReconstructed/hK0Short"), genPt, genEta, genPhi, generatedNch, k0Findable);
          }
          genStudyK0s.push_back(GenStudyPairObject{
            .pt = genPt,
            .eta = genEta,
            .phi = genPhi,
            .globalIndex = static_cast<int64_t>(mcParticle.globalIndex()),
            .motherIndex = -1,
            .reconstructed = reconstructedV0McIds.count(mcParticle.globalIndex()) > 0});
        }
      }

      // h-K0 correlations of the objects collected above, in generated coordinates.
      // Same delta-phi / delta-eta convention as every other correlation in this task
      // (trigger minus associated), and the same autocorrelation rejection: a trigger
      // that is a decay product of the K0 it would be paired with is skipped.
      // Every pair goes into Gen/ and into exactly one of the four exclusive classes.
      for (auto const& trigger : genStudyTriggers) {
        for (auto const& k0 : genStudyK0s) {
          if (trigger.globalIndex == k0.globalIndex || trigger.motherIndex == k0.globalIndex) {
            continue;
          }
          const float deltaPhi = computeDeltaPhi(trigger.phi, k0.phi);
          const float deltaEta = trigger.eta - k0.eta;
          histos.fill(HIST("PairLossK0/GenStudy/Gen/hCorrelation"), deltaEta, deltaPhi, trigger.pt, k0.pt, generatedNch);
          if (trigger.reconstructed && k0.reconstructed) {
            histos.fill(HIST("PairLossK0/GenStudy/Reconstructed/hCorrelation"), deltaEta, deltaPhi, trigger.pt, k0.pt, generatedNch);
          } else if (trigger.reconstructed) {
            histos.fill(HIST("PairLossK0/GenStudy/OnlyTriggerReconstructed/hCorrelation"), deltaEta, deltaPhi, trigger.pt, k0.pt, generatedNch);
          } else if (k0.reconstructed) {
            histos.fill(HIST("PairLossK0/GenStudy/OnlyK0Reconstructed/hCorrelation"), deltaEta, deltaPhi, trigger.pt, k0.pt, generatedNch);
          } else {
            histos.fill(HIST("PairLossK0/GenStudy/NotReconstructed/hCorrelation"), deltaEta, deltaPhi, trigger.pt, k0.pt, generatedNch);
          }
        }
      }
    };

    // Part 2. Wrapped in a lambda so that its own early exits leave the other
    // two parts free to run: the three parts are independent, not exclusive.
    auto runRecComparison = [&]() {
      if (recCollisions.size() == 0) {
        return;
      }

      int64_t bestCollisionId = -1;
      int largestNContributors = -1;
      float bestCollisionVtxX = 0.0f;
      float bestCollisionVtxY = 0.0f;
      float bestCollisionVtxZ = 0.0f;
      float bestCollisionMultiplicity = -1.0f;
      bool pairLossEventSelected = false;
      for (auto const& collision : recCollisions) {
        if (collision.numContrib() <= largestNContributors) {
          continue;
        }
        largestNContributors = collision.numContrib();
        bestCollisionId = collision.globalIndex();
        bestCollisionVtxX = collision.posX();
        bestCollisionVtxY = collision.posY();
        bestCollisionVtxZ = collision.posZ();
        bestCollisionMultiplicity = masterConfigurations.doPPAnalysis ? collision.centFT0M() : collision.centFT0C();
        pairLossEventSelected = !pairLossK0Configurations.applyRecoEventSelection ||
                                (masterConfigurations.doPPAnalysis ? isCollisionSelected(collision) : isCollisionSelectedPbPb(collision, false));
      }
      if (bestCollisionId < 0) {
        return;
      }

      // First build inclusive truth-acceptance lists for the early Rec control
      // stages. Exact old PairLoss truth lists are derived separately below so
      // that configured primary and leading-trigger requirements are applied
      // in precisely their original order.
      std::vector<PairLossTruthTrackInfo> comparisonTruthTriggers;
      std::vector<PairLossTruthK0Info> comparisonTruthK0s;
      for (auto const& mcParticle : mcParticles) {
        if (isPairLossTriggerPdg(mcParticle.pdgCode()) && std::abs(mcParticle.eta()) <= etaSel &&
            mcParticle.pt() >= axisRanges[3][0] && mcParticle.pt() <= axisRanges[3][1]) {
          auto const* pdgParticle = pdgDB->GetParticle(mcParticle.pdgCode());
          const double charge = pdgParticle != nullptr ? pdgParticle->Charge() : 0.0;
          const int sign = charge > 0.0 ? 1 : (charge < 0.0 ? -1 : 0);
          if (!((triggerParticleCharge > 0 && sign < 0) || (triggerParticleCharge < 0 && sign > 0) || sign == 0)) {
            comparisonTruthTriggers.push_back(PairLossTruthTrackInfo{
              .globalIndex = static_cast<int64_t>(mcParticle.globalIndex()),
              .pt = mcParticle.pt(),
              .eta = mcParticle.eta(),
              .phi = mcParticle.phi(),
              .sign = sign,
              .physicalPrimary = mcParticle.isPhysicalPrimary()});
          }
        }

        if (mcParticle.pdgCode() != PDG_t::kK0Short || std::abs(mcParticle.eta()) > etaSel ||
            mcParticle.pt() < axisRanges[2][0] || mcParticle.pt() > axisRanges[2][1]) {
          continue;
        }
        PairLossTruthK0Info truthK0{
          .globalIndex = static_cast<int64_t>(mcParticle.globalIndex()),
          .pt = mcParticle.pt(),
          .eta = mcParticle.eta(),
          .phi = mcParticle.phi(),
          .decayRadius = -1.0f,
          .findable = false,
          .physicalPrimary = mcParticle.isPhysicalPrimary(),
          .positiveDaughter = {},
          .negativeDaughter = {}};
        for (auto const& daughter : mcParticle.daughters_as<aod::McParticles>()) {
          PairLossTruthTrackInfo daughterInfo{
            .globalIndex = static_cast<int64_t>(daughter.globalIndex()),
            .pt = daughter.pt(),
            .eta = daughter.eta(),
            .phi = daughter.phi(),
            .sign = daughter.pdgCode() > 0 ? 1 : -1,
            .physicalPrimary = daughter.isPhysicalPrimary()};
          if (daughter.pdgCode() == PDG_t::kPiPlus) {
            truthK0.positiveDaughter = daughterInfo;
          } else if (daughter.pdgCode() == -PDG_t::kPiPlus) {
            truthK0.negativeDaughter = daughterInfo;
          }
        }
        comparisonTruthK0s.push_back(truthK0);
      }

      auto primaryTruthTriggers = comparisonTruthTriggers;
      if (masterConfigurations.doTriggPhysicalPrimary) {
        std::erase_if(primaryTruthTriggers, [](auto const& trigger) { return !trigger.physicalPrimary; });
      }
      auto truthTriggers = primaryTruthTriggers;
      if (useTheLeadingParticleAsTrigger && truthTriggers.size() > 1) {
        const auto leadingTrigger = std::max_element(truthTriggers.begin(), truthTriggers.end(), [](auto const& lhs, auto const& rhs) { return lhs.pt < rhs.pt; });
        const auto leadingTriggerCopy = *leadingTrigger;
        truthTriggers.clear();
        truthTriggers.push_back(leadingTriggerCopy);
      }
      auto truthK0s = comparisonTruthK0s;
      if (doAssocPhysicalPrimary) {
        std::erase_if(truthK0s, [](auto const& k0) { return !k0.physicalPrimary; });
      }

      // Truth: identical pair definition to the original PairLoss truth stage.
      if (pairLossEventSelected) {
        for (auto const& truthTrigger : truthTriggers) {
          for (auto const& truthK0 : truthK0s) {
            if (truthTrigger.globalIndex == truthK0.positiveDaughter.globalIndex || truthTrigger.globalIndex == truthK0.negativeDaughter.globalIndex) {
              continue;
            }
            const float truthDeltaPhi = computeDeltaPhi(truthTrigger.phi, truthK0.phi);
            float truthDeltaEta = truthTrigger.eta - truthK0.eta;
            if (masterConfigurations.doMirroringInDelataEta) {
              truthDeltaEta = std::abs(truthDeltaEta);
            }
            if (truthDeltaPhi < axisRanges[0][0] || truthDeltaPhi > axisRanges[0][1] ||
                truthDeltaEta < axisRanges[1][0] || truthDeltaEta > axisRanges[1][1]) {
              continue;
            }
            histos.fill(HIST("PairLossK0/Comparison/Truth"), truthDeltaPhi, truthDeltaEta, truthK0.pt, truthTrigger.pt, bestCollisionVtxZ, bestCollisionMultiplicity);
          }
        }
      }

      // Gen: reproduce the ClosureTest Gen event and particle/pair definition.
      float genBestCollisionMultiplicity = -1.0f;
      float genBestCollisionVtxZ = 0.0f;
      bool genBestCollisionSel8 = false;
      bool genBestCollisionINELgtZERO = false;
      bool genBestCollisionINELgtONE = false;
      bool genBestCollisionNoSameBunchPileup = false;
      bool genBestCollisionGoodTriggerTVX = false;
      bool genBestCollisionGoodZvtxFT0vsPV = false;
      bool genCollisionSelected = false;
      int genLargestNContributors = -1;
      uint32_t genBestCollisionTriggerPresenceMap = 0;
      for (auto const& recCollision : recCollisions) {
        if (genLargestNContributors >= recCollision.numContrib()) {
          continue;
        }
        genLargestNContributors = recCollision.numContrib();
        genBestCollisionMultiplicity = masterConfigurations.doPPAnalysis ? recCollision.centFT0M() : recCollision.centFT0C();
        if (masterConfigurations.applyNewMCSelection) {
          genCollisionSelected = (masterConfigurations.doPPAnalysis && isCollisionSelected(recCollision)) ||
                                 (!masterConfigurations.doPPAnalysis && isCollisionSelectedPbPb(recCollision, false));
        } else {
          genBestCollisionSel8 = recCollision.sel8();
          genBestCollisionVtxZ = recCollision.posZ();
          genBestCollisionINELgtZERO = recCollision.isInelGt0();
          genBestCollisionINELgtONE = recCollision.isInelGt1();
          genBestCollisionNoSameBunchPileup = recCollision.selection_bit(o2::aod::evsel::kNoSameBunchPileup);
          genBestCollisionGoodTriggerTVX = recCollision.selection_bit(aod::evsel::kIsTriggerTVX);
          genBestCollisionGoodZvtxFT0vsPV = recCollision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV);
        }
        if (triggerPresenceMap.size() > 0) {
          genBestCollisionTriggerPresenceMap = triggerPresenceMap[recCollision.globalIndex()];
        }
      }
      bool genEventSelected = triggerPresenceMap.size() == 0 || TESTBIT(genBestCollisionTriggerPresenceMap, triggerBinToSelect);
      if (masterConfigurations.applyNewMCSelection) {
        genEventSelected = genEventSelected && genCollisionSelected;
      } else if (masterConfigurations.doGenEventSelection) {
        genEventSelected = genEventSelected && genBestCollisionSel8 && std::abs(genBestCollisionVtxZ) <= masterConfigurations.zVertexCut &&
                           genBestCollisionINELgtZERO && (!masterConfigurations.selectINELgtONE || genBestCollisionINELgtONE) &&
                           (!masterConfigurations.rejectSameBunchPileup || genBestCollisionNoSameBunchPileup) &&
                           (!masterConfigurations.requireGoodTriggerTVX || genBestCollisionGoodTriggerTVX) &&
                           (!masterConfigurations.requireGoodZvtxFT0vsPV || genBestCollisionGoodZvtxFT0vsPV) &&
                           genBestCollisionMultiplicity >= axisRanges[5][0] && genBestCollisionMultiplicity <= axisRanges[5][1];
      }

      if (genEventSelected && masterConfigurations.doCorrelationK0Short && TESTBIT(doCorrelation, IndexK0)) {
        std::vector<uint32_t> genTriggerIndices;
        std::vector<uint32_t> genK0Indices;
        int iteratorNumber = -1;
        for (auto const& mcParticle : mcParticles) {
          ++iteratorNumber;
          if (std::abs(mcParticle.eta()) > etaSel) {
            continue;
          }
          if (isPairLossTriggerPdg(mcParticle.pdgCode()) && (!masterConfigurations.doTriggPhysicalPrimary || mcParticle.isPhysicalPrimary())) {
            genTriggerIndices.push_back(iteratorNumber);
          }
          if (mcParticle.pdgCode() == PDG_t::kK0Short && masterConfigurations.doCorrelationK0Short && (!doAssocPhysicalPrimary || mcParticle.isPhysicalPrimary())) {
            genK0Indices.push_back(iteratorNumber);
          }
        }
        for (auto const& triggerIndex : genTriggerIndices) {
          auto triggerParticle = mcParticles.iteratorAt(triggerIndex);
          if (triggerParticle.pt() > axisRanges[3][1] || triggerParticle.pt() < axisRanges[3][0]) {
            continue;
          }
          auto const& triggerMother = triggerParticle.mothers_first_as<aod::McParticles>();
          const auto triggerMotherIndex = triggerMother.globalIndex();
          for (auto const& k0Index : genK0Indices) {
            auto k0Particle = mcParticles.iteratorAt(k0Index);
            if (triggerIndex == k0Index || triggerMotherIndex == k0Particle.globalIndex()) {
              continue;
            }
            const float deltaPhi = computeDeltaPhi(triggerParticle.phi(), k0Particle.phi());
            const float deltaEta = triggerParticle.eta() - k0Particle.eta();
            if (deltaPhi < axisRanges[0][0] || deltaPhi > axisRanges[0][1] || deltaEta < axisRanges[1][0] || deltaEta > axisRanges[1][1] ||
                k0Particle.pt() < axisRanges[2][0] || k0Particle.pt() > axisRanges[2][1]) {
              continue;
            }
            histos.fill(HIST("PairLossK0/Comparison/Gen"), deltaPhi, deltaEta, k0Particle.pt(), triggerParticle.pt(), genBestCollisionVtxZ, genBestCollisionMultiplicity);
          }
        }
      }

      // Old Final control: reproduce the original PairLoss final-pair
      // definition directly.  This is deliberately independent of the Rec
      // control chain below: it uses the old best-collision object maps,
      // reconstructed-pair range and autocorrelation decision, then fills
      // each accepted truth pair once with unit weight and truth coordinates.
      PairLossTrackMap oldFinalTriggers;
      const auto oldFinalTriggerSlice = triggerTracks.sliceBy(collisionSliceTracks, bestCollisionId);
      for (auto const& triggerEntry : oldFinalTriggerSlice) {
        auto track = triggerEntry.track_as<TracksComplete>();
        const auto trackLabel = trackLabels.iteratorAt(track.globalIndex());
        if (!trackLabel.has_mcParticle()) {
          continue;
        }
        if (!isValidTrigger(track, triggerEntry.isLeading())) {
          continue;
        }
        if (trackSelection.checkForITSTPCMissmatchMC && bitcheck(triggerEntry.mcMask(), 13)) {
          continue;
        }
        if (masterConfigurations.doTriggPhysicalPrimary && !triggerEntry.mcPhysicalPrimary()) {
          continue;
        }
        oldFinalTriggers[trackLabel.mcParticleId()].push_back(makePairLossTrackInfo(track));
      }

      PairLossV0Map oldFinalV0s;
      const auto oldFinalV0Slice = associatedV0s.sliceBy(collisionSliceV0s, bestCollisionId);
      for (auto const& assocEntry : oldFinalV0Slice) {
        auto v0 = assocEntry.v0Core_as<V0DatasWithoutTrackX>();
        const auto v0MC = v0MCCores.iteratorAt(v0.globalIndex());
        if (v0MC.particleIdMC() < 0 || v0MC.pdgCode() != PDG_t::kK0Short || !assocEntry.mcTrue(IndexK0)) {
          continue;
        }
        auto positiveTrack = v0.posTrack_as<TracksComplete>();
        auto negativeTrack = v0.negTrack_as<TracksComplete>();
        bool passesFinalSelection = true;
        if (masterConfigurations.doPPAnalysis) {
          passesFinalSelection = v0.v0radius() >= v0Selection.v0RadiusMin && v0.v0radius() <= v0Selection.v0RadiusMax &&
                                 std::abs(v0.dcapostopv()) >= v0Selection.dcapostopv && std::abs(v0.dcanegtopv()) >= v0Selection.dcanegtopv &&
                                 v0.v0cosPA() >= v0Selection.v0cospa && v0.dcaV0daughters() <= v0Selection.dcaV0dau;
        } else {
          const float dcaCut = v0Selection.dcaDaugToPVForK0s == 0.0f ? v0Selection.dcaMesonToPV : v0Selection.dcaDaugToPVForK0s;
          const bool passesLifetime = v0.distovertotmom(bestCollisionVtxX, bestCollisionVtxY, bestCollisionVtxZ) *
                                        o2::constants::physics::MassK0Short <
                                      v0Selection.lifetimecutK0S;
          const bool passesDaughterDcaAndArmenteros = std::abs(v0.dcapostopv()) > dcaCut && std::abs(v0.dcanegtopv()) > dcaCut && v0.qtarm() * v0Selection.armPodCut > std::abs(v0.alpha());
          passesFinalSelection = v0SelectedPbPb(v0) && passesLifetime && passesDaughterDcaAndArmenteros;
        }
        passesFinalSelection = passesFinalSelection && positiveTrack.tpcNClsCrossedRows() >= trackSelection.minTPCNCrossedRowsAssociated &&
                               negativeTrack.tpcNClsCrossedRows() >= trackSelection.minTPCNCrossedRowsAssociated;
        if (trackSelection.checksRequireTPCChi2) {
          passesFinalSelection = passesFinalSelection && positiveTrack.tpcChi2NCl() >= trackSelection.minTPCChi2PerClusterAssociated &&
                                 negativeTrack.tpcChi2NCl() >= trackSelection.minTPCChi2PerClusterAssociated;
        }
        if (trackSelection.requireClusterInITS) {
          passesFinalSelection = passesFinalSelection && positiveTrack.itsNCls() >= trackSelection.minITSClustersForDaughterTracks &&
                                 negativeTrack.itsNCls() >= trackSelection.minITSClustersForDaughterTracks;
        }
        passesFinalSelection = passesFinalSelection && assocEntry.compatible(IndexK0, trackSelection.dEdxCompatibility) &&
                               (!doAssocPhysicalPrimary || assocEntry.mcPhysicalPrimary()) &&
                               assocEntry.invMassNSigma(IndexK0) > -massWindowConfigurations.maxPeakNSigma &&
                               assocEntry.invMassNSigma(IndexK0) < massWindowConfigurations.maxPeakNSigma &&
                               v0.pt() >= axisRanges[2][0] && v0.pt() <= axisRanges[2][1];
        if (!passesFinalSelection) {
          continue;
        }
        oldFinalV0s[v0MC.particleIdMC()].push_back(PairLossV0Info{
          .globalIndex = static_cast<int64_t>(v0.globalIndex()),
          .positiveTrackId = static_cast<int64_t>(positiveTrack.globalIndex()),
          .negativeTrackId = static_cast<int64_t>(negativeTrack.globalIndex()),
          .pt = v0.pt(),
          .eta = v0.eta(),
          .phi = v0.phi(),
          .radius = v0.v0radius(),
          .cosPA = v0.v0cosPA(),
          .dcaDaughters = v0.dcaV0daughters(),
          .massNSigma = assocEntry.invMassNSigma(IndexK0)});
      }

      // Final pairs are recorded here and compared with Rec variant 14 after
      // the exact Rec replica below has run.
      struct PairLossFinalRecord {
        PairLossPairKey key{};
        float deltaPhi = 0.0f;
        float deltaEta = 0.0f;
        float k0Pt = 0.0f;
        float triggerPt = 0.0f;
      };
      std::vector<PairLossFinalRecord> finalPairRecords;

      if (pairLossEventSelected) {
        for (auto const& truthTrigger : truthTriggers) {
          if (oldFinalTriggers.find(truthTrigger.globalIndex) == oldFinalTriggers.end()) {
            continue;
          }
          histos.fill(HIST("PairLossK0/Comparison/TriggerVariants"), 9, truthTrigger.pt, bestCollisionVtxZ, bestCollisionMultiplicity);
        }

        for (auto const& truthTrigger : truthTriggers) {
          const auto triggerMatches = oldFinalTriggers.find(truthTrigger.globalIndex);
          if (triggerMatches == oldFinalTriggers.end()) {
            continue;
          }
          for (auto const& truthK0 : truthK0s) {
            if (truthTrigger.globalIndex == truthK0.positiveDaughter.globalIndex || truthTrigger.globalIndex == truthK0.negativeDaughter.globalIndex) {
              continue;
            }
            float truthDeltaEta = truthTrigger.eta - truthK0.eta;
            if (masterConfigurations.doMirroringInDelataEta) {
              truthDeltaEta = std::abs(truthDeltaEta);
            }
            const float truthDeltaPhi = computeDeltaPhi(truthTrigger.phi, truthK0.phi);
            if (truthDeltaPhi < axisRanges[0][0] || truthDeltaPhi > axisRanges[0][1] ||
                truthDeltaEta < axisRanges[1][0] || truthDeltaEta > axisRanges[1][1]) {
              continue;
            }
            const auto v0Matches = oldFinalV0s.find(truthK0.globalIndex);
            if (v0Matches == oldFinalV0s.end()) {
              continue;
            }
            bool oldFinalPair = false;
            for (auto const& reconstructedTrigger : triggerMatches->second) {
              for (auto const& reconstructedV0 : v0Matches->second) {
                float reconstructedDeltaEta = reconstructedTrigger.eta - reconstructedV0.eta;
                if (masterConfigurations.doMirroringInDelataEta) {
                  reconstructedDeltaEta = std::abs(reconstructedDeltaEta);
                }
                const float reconstructedDeltaPhi = computeDeltaPhi(reconstructedTrigger.phi, reconstructedV0.phi);
                if (reconstructedDeltaPhi < axisRanges[0][0] || reconstructedDeltaPhi > axisRanges[0][1] ||
                    reconstructedDeltaEta < axisRanges[1][0] || reconstructedDeltaEta > axisRanges[1][1]) {
                  continue;
                }
                if (doAutocorrelationRejection && (reconstructedTrigger.globalIndex == reconstructedV0.positiveTrackId || reconstructedTrigger.globalIndex == reconstructedV0.negativeTrackId)) {
                  continue;
                }
                oldFinalPair = true;
                break;
              }
              if (oldFinalPair) {
                break;
              }
            }
            if (oldFinalPair) {
              histos.fill(HIST("PairLossK0/Comparison/Final"), truthDeltaPhi, truthDeltaEta, truthK0.pt, truthTrigger.pt, bestCollisionVtxZ, bestCollisionMultiplicity);
              histos.fill(HIST("PairLossK0/Comparison/RecVariants"), 15, truthDeltaPhi, truthDeltaEta, truthK0.pt, truthTrigger.pt, bestCollisionVtxZ, bestCollisionMultiplicity);
              finalPairRecords.push_back(PairLossFinalRecord{
                PairLossPairKey{truthTrigger.globalIndex, truthK0.globalIndex},
                truthDeltaPhi, truthDeltaEta, truthK0.pt, truthTrigger.pt});
            }
          }
        }
      }

      // Prepare MC-label lookup for the hook in the exact Rec signal branch.
      pairLossComparison.clear();
      for (auto const& collision : recCollisions) {
        for (auto const& track : tracks.sliceBy(pairLossTracksPerCollision, collision.globalIndex())) {
          const auto trackLabel = trackLabels.iteratorAt(track.globalIndex());
          if (trackLabel.has_mcParticle()) {
            pairLossComparison.trackToMc[track.globalIndex()] = trackLabel.mcParticleId();
          }
        }
        for (auto const& v0 : v0Candidates.sliceBy(pairLossV0sPerCollision, collision.globalIndex())) {
          const auto v0MC = v0MCCores.iteratorAt(v0.globalIndex());
          if (v0MC.particleIdMC() >= 0 && v0MC.pdgCode() == PDG_t::kK0Short) {
            pairLossComparison.v0ToMc[v0.globalIndex()] = v0MC.particleIdMC();
          }
        }
      }
      for (auto const& truthTrigger : comparisonTruthTriggers) {
        pairLossComparison.truthTriggers[truthTrigger.globalIndex] = truthTrigger;
      }
      for (auto const& truthK0 : comparisonTruthK0s) {
        pairLossComparison.truthK0s[truthK0.globalIndex] = truthK0;
      }
      for (auto const& truthTrigger : primaryTruthTriggers) {
        pairLossComparison.primaryTruthTriggerIds.insert(truthTrigger.globalIndex);
      }
      for (auto const& truthK0 : truthK0s) {
        pairLossComparison.primaryTruthK0Ids.insert(truthK0.globalIndex);
      }
      for (auto const& truthTrigger : truthTriggers) {
        pairLossComparison.oldTruthTriggerIds.insert(truthTrigger.globalIndex);
      }
      for (auto const& truthK0 : truthK0s) {
        pairLossComparison.oldTruthK0Ids.insert(truthK0.globalIndex);
      }
      for (auto const& record : finalPairRecords) {
        pairLossComparison.finalTruthPairs.insert(record.key);
      }
      pairLossComparison.bestCollisionId = bestCollisionId;
      pairLossComparison.vtxZ = bestCollisionVtxZ;
      pairLossComparison.multiplicity = bestCollisionMultiplicity;
      pairLossComparison.active = true;
      for (auto const& recCollision : recCollisions) {
        const auto recTriggerSlice = triggerTracks.sliceBy(collisionSliceTracks, recCollision.globalIndex());
        const auto recV0Slice = associatedV0s.sliceBy(collisionSliceV0s, recCollision.globalIndex());
        const auto recTrackSlice = tracks.sliceBy(pairLossTracksPerCollision, recCollision.globalIndex());
        runSameEventHV0s(recCollision, recV0Slice, recTriggerSlice, recTrackSlice);
      }

      // Final-minus-Rec counterpart of RecNotInFinal. Both set differences use
      // exactly the same truth-pair key and the same unit-weight truth axes.
      for (auto const& record : finalPairRecords) {
        if (pairLossComparison.countedPairs.contains(record.key)) {
          continue;
        }
        histos.fill(HIST("PairLossK0/Comparison/FinalNotInRec"), record.deltaPhi, record.deltaEta, record.k0Pt, record.triggerPt, pairLossComparison.vtxZ, pairLossComparison.multiplicity);
      }

      pairLossComparison.clear();
    };

    if (pairLossK0Configurations.doGenLevelStudy) {
      runGenLevelStudy();
    }
    if (pairLossK0Configurations.doRecComparison) {
      runRecComparison();
    }
    if (!pairLossK0Configurations.doStageDiagnostics) {
      return;
    }

    // Part 1 follows.
    histos.fill(HIST("PairLossK0/Event/hCounter"), 0.0f);
    histos.fill(HIST("PairLossK0/Event/hNRecoCollisions"), recCollisions.size());

    if (recCollisions.size() == 0) {
      return;
    }
    histos.fill(HIST("PairLossK0/Event/hCounter"), 1.0f);

    int64_t bestCollisionId = -1;
    int largestNContributors = -1;
    for (auto const& collision : recCollisions) {
      if (collision.numContrib() > largestNContributors) {
        largestNContributors = collision.numContrib();
        bestCollisionId = collision.globalIndex();
      }
    }
    if (bestCollisionId < 0) {
      return;
    }

    PairLossTrackMap tracksAnyCollision;
    PairLossV0Map v0sAnyCollision;
    for (auto const& collision : recCollisions) {
      const auto trackSlice = tracks.sliceBy(pairLossTracksPerCollision, collision.globalIndex());
      for (auto const& track : trackSlice) {
        const auto trackLabel = trackLabels.iteratorAt(track.globalIndex());
        if (!trackLabel.has_mcParticle()) {
          continue;
        }
        tracksAnyCollision[trackLabel.mcParticleId()].push_back(makePairLossTrackInfo(track));
      }
      const auto v0Slice = v0Candidates.sliceBy(pairLossV0sPerCollision, collision.globalIndex());
      for (auto const& v0 : v0Slice) {
        const auto v0MC = v0MCCores.iteratorAt(v0.globalIndex());
        if (v0MC.particleIdMC() < 0 || v0MC.pdgCode() != PDG_t::kK0Short) {
          continue;
        }
        v0sAnyCollision[v0MC.particleIdMC()].push_back(PairLossV0Info{
          .globalIndex = static_cast<int64_t>(v0.globalIndex()),
          .positiveTrackId = static_cast<int64_t>(v0.posTrackId()),
          .negativeTrackId = static_cast<int64_t>(v0.negTrackId()),
          .pt = v0.pt(),
          .eta = v0.eta(),
          .phi = v0.phi(),
          .radius = v0.v0radius(),
          .cosPA = v0.v0cosPA(),
          .dcaDaughters = v0.dcaV0daughters(),
          .massNSigma = 0.0f});
      }
    }

    for (auto const& collision : recCollisions) {
      if (static_cast<int64_t>(collision.globalIndex()) != bestCollisionId) {
        continue;
      }

      const bool collisionSelected = !pairLossK0Configurations.applyRecoEventSelection ||
                                     (masterConfigurations.doPPAnalysis ? isCollisionSelected(collision) : isCollisionSelectedPbPb(collision, false));
      if (!collisionSelected) {
        return;
      }
      histos.fill(HIST("PairLossK0/Event/hCounter"), 2.0f);

      const auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      const double magneticField = getPairLossMagneticField(bc.runNumber(), bc.timestamp());
      const int magneticFieldSign = magneticField > 0.0 ? 1 : (magneticField < 0.0 ? -1 : 0);
      const float multiplicity = masterConfigurations.doPPAnalysis ? collision.centFT0M() : collision.centFT0C();

      PairLossTrackMap tracksBestCollision;
      const auto bestTrackSlice = tracks.sliceBy(pairLossTracksPerCollision, bestCollisionId);
      for (auto const& track : bestTrackSlice) {
        const auto trackLabel = trackLabels.iteratorAt(track.globalIndex());
        if (!trackLabel.has_mcParticle()) {
          continue;
        }
        tracksBestCollision[trackLabel.mcParticleId()].push_back(makePairLossTrackInfo(track));
      }

      PairLossV0Map v0sBestCollision;
      const auto bestV0Slice = v0Candidates.sliceBy(pairLossV0sPerCollision, bestCollisionId);
      for (auto const& v0 : bestV0Slice) {
        const auto v0MC = v0MCCores.iteratorAt(v0.globalIndex());
        if (v0MC.particleIdMC() < 0 || v0MC.pdgCode() != PDG_t::kK0Short) {
          continue;
        }
        v0sBestCollision[v0MC.particleIdMC()].push_back(PairLossV0Info{
          .globalIndex = static_cast<int64_t>(v0.globalIndex()),
          .positiveTrackId = static_cast<int64_t>(v0.posTrackId()),
          .negativeTrackId = static_cast<int64_t>(v0.negTrackId()),
          .pt = v0.pt(),
          .eta = v0.eta(),
          .phi = v0.phi(),
          .radius = v0.v0radius(),
          .cosPA = v0.v0cosPA(),
          .dcaDaughters = v0.dcaV0daughters(),
          .massNSigma = 0.0f});
      }

      PairLossTrackMap triggersInTable;
      PairLossTrackMap triggersFinal;
      PairLossV0Map v0sInTable;
      PairLossV0Map v0sFinal;
      for (auto const& finalCollision : recCollisions) {
        const int64_t finalCollisionId = static_cast<int64_t>(finalCollision.globalIndex());
        const bool isBestCollision = finalCollisionId == bestCollisionId;
        if (!isBestCollision) {
          if (pairLossK0Configurations.finalPairUseBestCollisionOnly) {
            continue;
          }
          if (pairLossK0Configurations.applyRecoEventSelection &&
              !(masterConfigurations.doPPAnalysis ? isCollisionSelected(finalCollision) : isCollisionSelectedPbPb(finalCollision, false))) {
            continue;
          }
        }

        const auto triggerSlice = triggerTracks.sliceBy(collisionSliceTracks, finalCollisionId);
        for (auto const& triggerEntry : triggerSlice) {
          auto track = triggerEntry.track_as<TracksComplete>();
          const auto trackLabel = trackLabels.iteratorAt(track.globalIndex());
          if (!trackLabel.has_mcParticle()) {
            continue;
          }
          const int64_t mcParticleId = trackLabel.mcParticleId();
          auto trackInfo = makePairLossTrackInfo(track);
          trackInfo.collisionId = finalCollisionId;
          if (isBestCollision) {
            triggersInTable[mcParticleId].push_back(trackInfo);
          }
          if (!isValidTrigger(track, triggerEntry.isLeading())) {
            continue;
          }
          if (trackSelection.checkForITSTPCMissmatchMC && bitcheck(triggerEntry.mcMask(), 13)) {
            continue;
          }
          if (masterConfigurations.doTriggPhysicalPrimary && !triggerEntry.mcPhysicalPrimary()) {
            continue;
          }
          triggersFinal[mcParticleId].push_back(trackInfo);
        }

        const auto assocV0Slice = associatedV0s.sliceBy(collisionSliceV0s, finalCollisionId);
        for (auto const& assocEntry : assocV0Slice) {
          auto v0 = assocEntry.v0Core_as<V0DatasWithoutTrackX>();
          const auto v0MC = v0MCCores.iteratorAt(v0.globalIndex());
          if (v0MC.particleIdMC() < 0 || v0MC.pdgCode() != PDG_t::kK0Short || !assocEntry.mcTrue(IndexK0)) {
            continue;
          }
          auto positiveTrack = v0.posTrack_as<TracksComplete>();
          auto negativeTrack = v0.negTrack_as<TracksComplete>();
          const PairLossV0Info v0Info{
            .globalIndex = static_cast<int64_t>(v0.globalIndex()),
            .positiveTrackId = static_cast<int64_t>(positiveTrack.globalIndex()),
            .negativeTrackId = static_cast<int64_t>(negativeTrack.globalIndex()),
            .pt = v0.pt(),
            .eta = v0.eta(),
            .phi = v0.phi(),
            .radius = v0.v0radius(),
            .cosPA = v0.v0cosPA(),
            .dcaDaughters = v0.dcaV0daughters(),
            .massNSigma = assocEntry.invMassNSigma(IndexK0),
            .collisionId = finalCollisionId};
          if (isBestCollision) {
            v0sInTable[v0MC.particleIdMC()].push_back(v0Info);
          }

          bool passesFinalSelection = true;
          if (masterConfigurations.doPPAnalysis) {
            passesFinalSelection = v0.v0radius() >= v0Selection.v0RadiusMin && v0.v0radius() <= v0Selection.v0RadiusMax &&
                                   std::abs(v0.dcapostopv()) >= v0Selection.dcapostopv && std::abs(v0.dcanegtopv()) >= v0Selection.dcanegtopv &&
                                   v0.v0cosPA() >= v0Selection.v0cospa && v0.dcaV0daughters() <= v0Selection.dcaV0dau;
          } else {
            const float dcaCut = v0Selection.dcaDaugToPVForK0s == 0.0f ? v0Selection.dcaMesonToPV : v0Selection.dcaDaugToPVForK0s;
            const bool passesLifetime = v0.distovertotmom(finalCollision.posX(), finalCollision.posY(), finalCollision.posZ()) * o2::constants::physics::MassK0Short < v0Selection.lifetimecutK0S;
            const bool passesDaughterDcaAndArmenteros = std::abs(v0.dcapostopv()) > dcaCut && std::abs(v0.dcanegtopv()) > dcaCut && v0.qtarm() * v0Selection.armPodCut > std::abs(v0.alpha());
            passesFinalSelection = v0SelectedPbPb(v0) && passesLifetime && passesDaughterDcaAndArmenteros;
          }
          passesFinalSelection = passesFinalSelection && positiveTrack.tpcNClsCrossedRows() >= trackSelection.minTPCNCrossedRowsAssociated &&
                                 negativeTrack.tpcNClsCrossedRows() >= trackSelection.minTPCNCrossedRowsAssociated;
          if (trackSelection.checksRequireTPCChi2) {
            passesFinalSelection = passesFinalSelection && positiveTrack.tpcChi2NCl() >= trackSelection.minTPCChi2PerClusterAssociated &&
                                   negativeTrack.tpcChi2NCl() >= trackSelection.minTPCChi2PerClusterAssociated;
          }
          if (trackSelection.requireClusterInITS) {
            passesFinalSelection = passesFinalSelection && positiveTrack.itsNCls() >= trackSelection.minITSClustersForDaughterTracks &&
                                   negativeTrack.itsNCls() >= trackSelection.minITSClustersForDaughterTracks;
          }
          passesFinalSelection = passesFinalSelection && assocEntry.compatible(IndexK0, trackSelection.dEdxCompatibility) &&
                                 (!doAssocPhysicalPrimary || assocEntry.mcPhysicalPrimary()) &&
                                 assocEntry.invMassNSigma(IndexK0) > -massWindowConfigurations.maxPeakNSigma &&
                                 assocEntry.invMassNSigma(IndexK0) < massWindowConfigurations.maxPeakNSigma &&
                                 v0.pt() >= axisRanges[2][0] && v0.pt() <= axisRanges[2][1];
          if (passesFinalSelection) {
            v0sFinal[v0MC.particleIdMC()].push_back(v0Info);
          }
        }
      }

      std::vector<PairLossTruthTrackInfo> truthTriggers;
      std::vector<PairLossTruthK0Info> truthK0s;
      for (auto const& mcParticle : mcParticles) {
        if (isPairLossTriggerPdg(mcParticle.pdgCode()) && std::abs(mcParticle.eta()) <= etaSel &&
            mcParticle.pt() >= axisRanges[3][0] && mcParticle.pt() <= axisRanges[3][1] &&
            (!masterConfigurations.doTriggPhysicalPrimary || mcParticle.isPhysicalPrimary())) {
          auto const* pdgParticle = pdgDB->GetParticle(mcParticle.pdgCode());
          const double charge = pdgParticle != nullptr ? pdgParticle->Charge() : 0.0;
          const int sign = charge > 0.0 ? 1 : (charge < 0.0 ? -1 : 0);
          if ((triggerParticleCharge > 0 && sign < 0) || (triggerParticleCharge < 0 && sign > 0) || sign == 0) {
            continue;
          }
          truthTriggers.push_back(PairLossTruthTrackInfo{
            .globalIndex = static_cast<int64_t>(mcParticle.globalIndex()),
            .pt = mcParticle.pt(),
            .eta = mcParticle.eta(),
            .phi = mcParticle.phi(),
            .sign = sign});
        }

        if (mcParticle.pdgCode() != PDG_t::kK0Short || std::abs(mcParticle.eta()) > etaSel ||
            mcParticle.pt() < axisRanges[2][0] || mcParticle.pt() > axisRanges[2][1] ||
            (doAssocPhysicalPrimary && !mcParticle.isPhysicalPrimary())) {
          continue;
        }
        PairLossTruthK0Info truthK0{
          .globalIndex = static_cast<int64_t>(mcParticle.globalIndex()),
          .pt = mcParticle.pt(),
          .eta = mcParticle.eta(),
          .phi = mcParticle.phi(),
          .decayRadius = -1.0f,
          .findable = false,
          .positiveDaughter = {},
          .negativeDaughter = {}};
        const auto daughters = mcParticle.daughters_as<aod::McParticles>();
        for (auto const& daughter : daughters) {
          PairLossTruthTrackInfo daughterInfo{
            .globalIndex = static_cast<int64_t>(daughter.globalIndex()),
            .pt = daughter.pt(),
            .eta = daughter.eta(),
            .phi = daughter.phi(),
            .sign = daughter.pdgCode() > 0 ? 1 : -1};
          if (daughter.pdgCode() == PDG_t::kPiPlus) {
            truthK0.positiveDaughter = daughterInfo;
          } else if (daughter.pdgCode() == -PDG_t::kPiPlus) {
            truthK0.negativeDaughter = daughterInfo;
          } else {
            continue;
          }
          const float daughterDecayRadius = std::hypot(daughter.vx() - mcCollision.posX(), daughter.vy() - mcCollision.posY());
          if (truthK0.decayRadius < 0.0f) {
            truthK0.decayRadius = daughterDecayRadius;
          }
        }
        const bool hasPionDaughters = truthK0.positiveDaughter.globalIndex >= 0 && truthK0.negativeDaughter.globalIndex >= 0;
        truthK0.findable = hasPionDaughters && truthK0.positiveDaughter.pt >= pairLossK0Configurations.daughterPtMin &&
                           truthK0.negativeDaughter.pt >= pairLossK0Configurations.daughterPtMin &&
                           std::abs(truthK0.positiveDaughter.eta) <= pairLossK0Configurations.daughterEtaMax &&
                           std::abs(truthK0.negativeDaughter.eta) <= pairLossK0Configurations.daughterEtaMax;
        truthK0s.push_back(truthK0);
      }

      if (useTheLeadingParticleAsTrigger && truthTriggers.size() > 1) {
        const auto leadingTrigger = std::max_element(truthTriggers.begin(), truthTriggers.end(), [](auto const& lhs, auto const& rhs) { return lhs.pt < rhs.pt; });
        const auto leadingTriggerCopy = *leadingTrigger;
        truthTriggers.clear();
        truthTriggers.push_back(leadingTriggerCopy);
      }
      if (!truthTriggers.empty()) {
        histos.fill(HIST("PairLossK0/Event/hCounter"), 3.0f);
      }
      if (!truthK0s.empty()) {
        histos.fill(HIST("PairLossK0/Event/hCounter"), 4.0f);
      }

      auto contains = [](auto const& map, int64_t key) {
        auto const iterator = map.find(key);
        return iterator != map.end() && !iterator->second.empty();
      };
      auto numberOfMatches = [](auto const& map, int64_t key) -> size_t {
        auto const iterator = map.find(key);
        return iterator == map.end() ? 0 : iterator->second.size();
      };
      auto bestTrack = [](PairLossTrackMap const& map, int64_t key) -> PairLossTrackInfo const* {
        auto const iterator = map.find(key);
        if (iterator == map.end() || iterator->second.empty()) {
          return nullptr;
        }
        return &*std::max_element(iterator->second.begin(), iterator->second.end(), [](auto const& lhs, auto const& rhs) {
          if (lhs.tpcCrossedRows != rhs.tpcCrossedRows) {
            return lhs.tpcCrossedRows < rhs.tpcCrossedRows;
          }
          return lhs.itsClusters < rhs.itsClusters;
        });
      };
      auto bestV0 = [](PairLossV0Map const& map, int64_t key) -> PairLossV0Info const* {
        auto const iterator = map.find(key);
        if (iterator == map.end() || iterator->second.empty()) {
          return nullptr;
        }
        return &*std::max_element(iterator->second.begin(), iterator->second.end(), [](auto const& lhs, auto const& rhs) { return lhs.cosPA < rhs.cosPA; });
      };

      // Separate from the pair loop, which would count a trigger once per truth K0 and drop it entirely in events without one.
      // Each stage holds the trigger-side condition of the same stage of hCounts/hPhysics, so the ratio keeps one trigger population.
      for (auto const& truthTrigger : truthTriggers) {
        const std::array<bool, PairLossK0NStages> triggerStagePassed = {
          true,                                                    // Gen pair
          true,                                                    // Findable K0->pi+pi-: K0 side only
          contains(tracksBestCollision, truthTrigger.globalIndex), // Trigger pure reco (best collision)
          contains(triggersInTable, truthTrigger.globalIndex),     // Trigger in TriggerTracks
          contains(triggersFinal, truthTrigger.globalIndex),       // Trigger final selection
          true,                                                    // V0 pure reco: K0 side only
          true,                                                    // V0 candidate: K0 side only
          true,                                                    // V0 in AssocV0s: K0 side only
          true,                                                    // V0 final selection: K0 side only
          contains(triggersFinal, truthTrigger.globalIndex)};      // Final reconstructed pair
        std::array<int, PairLossK0NStages> triggerStageFills;
        triggerStageFills.fill(1);
        if (!pairLossK0Configurations.fillFinalPairOnce) {
          triggerStageFills[PairLossTriggerPureReco] = static_cast<int>(numberOfMatches(tracksBestCollision, truthTrigger.globalIndex));
          triggerStageFills[PairLossTriggerInTable] = static_cast<int>(numberOfMatches(triggersInTable, truthTrigger.globalIndex));
          triggerStageFills[PairLossTriggerFinal] = static_cast<int>(numberOfMatches(triggersFinal, truthTrigger.globalIndex));
          triggerStageFills[PairLossFinalPair] = static_cast<int>(numberOfMatches(triggersFinal, truthTrigger.globalIndex));
        }
        for (int stage = 0; stage < PairLossK0NStages; ++stage) {
          if (!triggerStagePassed[stage]) {
            continue;
          }
          for (int fill = 0; fill < triggerStageFills[stage]; ++fill) {
            histos.fill(HIST("PairLossK0/Stage/hTriggers"), stage, truthTrigger.pt, truthTrigger.eta, truthTrigger.phi, collision.posZ(), multiplicity);
          }
        }
      }

      bool hasTruthPair = false;
      for (auto const& truthTrigger : truthTriggers) {
        for (auto const& truthK0 : truthK0s) {
          if (truthTrigger.globalIndex == truthK0.positiveDaughter.globalIndex || truthTrigger.globalIndex == truthK0.negativeDaughter.globalIndex) {
            continue;
          }
          const float truthDeltaPhi = computeDeltaPhi(truthTrigger.phi, truthK0.phi);
          float truthDeltaEta = truthTrigger.eta - truthK0.eta;
          if (masterConfigurations.doMirroringInDelataEta) {
            truthDeltaEta = std::abs(truthDeltaEta);
          }
          if (truthDeltaPhi < axisRanges[0][0] || truthDeltaPhi > axisRanges[0][1] || truthDeltaEta < axisRanges[1][0] || truthDeltaEta > axisRanges[1][1]) {
            continue;
          }
          hasTruthPair = true;

          const auto positiveDeltaPhiStar = calculateMinimumDeltaPhiStar(truthTrigger, truthK0.positiveDaughter, magneticField, truthK0.decayRadius);
          const auto negativeDeltaPhiStar = calculateMinimumDeltaPhiStar(truthTrigger, truthK0.negativeDaughter, magneticField, truthK0.decayRadius);
          PairLossDeltaPhiStarInfo closestDeltaPhiStar;
          PairLossTruthTrackInfo const* closestDaughter = nullptr;
          if (positiveDeltaPhiStar.valid && (!negativeDeltaPhiStar.valid || positiveDeltaPhiStar.minAbs <= negativeDeltaPhiStar.minAbs)) {
            closestDeltaPhiStar = positiveDeltaPhiStar;
            closestDaughter = &truthK0.positiveDaughter;
          } else if (negativeDeltaPhiStar.valid) {
            closestDeltaPhiStar = negativeDeltaPhiStar;
            closestDaughter = &truthK0.negativeDaughter;
          }
          const float closestDeltaEta = closestDaughter != nullptr ? truthTrigger.eta - closestDaughter->eta : 999.0f;
          const int closestChargeProduct = closestDaughter != nullptr ? truthTrigger.sign * closestDaughter->sign : 0;

          if (truthK0.findable) {
            histos.fill(HIST("PairLossK0/Geometry/hDeltaPhiStarValidity"), closestDeltaPhiStar.valid, truthK0.pt);
            if (positiveDeltaPhiStar.valid) {
              histos.fill(HIST("PairLossK0/Geometry/hGeneratedDaughters"), positiveDeltaPhiStar.signedAtMin, truthTrigger.eta - truthK0.positiveDaughter.eta, truthK0.pt, truthTrigger.pt, magneticFieldSign, truthTrigger.sign * truthK0.positiveDaughter.sign);
            }
            if (negativeDeltaPhiStar.valid) {
              histos.fill(HIST("PairLossK0/Geometry/hGeneratedDaughters"), negativeDeltaPhiStar.signedAtMin, truthTrigger.eta - truthK0.negativeDaughter.eta, truthK0.pt, truthTrigger.pt, magneticFieldSign, truthTrigger.sign * truthK0.negativeDaughter.sign);
            }
            if (closestDeltaPhiStar.valid) {
              histos.fill(HIST("PairLossK0/Geometry/hRawDeltaPhiVsSignedDeltaPhiStar"), computeDeltaPhi(truthTrigger.phi, closestDaughter->phi), closestDeltaPhiStar.signedAtMin, magneticFieldSign);
              histos.fill(HIST("PairLossK0/Geometry/hGeneratedMinDeltaPhiStarVsDecayRadius"), closestDeltaPhiStar.minAbs, truthK0.decayRadius, truthK0.pt);
              histos.fill(HIST("PairLossK0/Geometry/hGeneratedRadiusAtMinimum"), closestDeltaPhiStar.radiusAtMin, closestDeltaPhiStar.minAbs, truthK0.pt);
            }
          }

          const bool triggerBestCollision = contains(tracksBestCollision, truthTrigger.globalIndex);
          const bool triggerInTable = contains(triggersInTable, truthTrigger.globalIndex);
          const bool triggerFinal = contains(triggersFinal, truthTrigger.globalIndex);
          const bool positiveDaughterBestCollision = contains(tracksBestCollision, truthK0.positiveDaughter.globalIndex);
          const bool negativeDaughterBestCollision = contains(tracksBestCollision, truthK0.negativeDaughter.globalIndex);
          const bool bothDaughtersBestCollision = positiveDaughterBestCollision && negativeDaughterBestCollision;
          const bool v0BestCollision = contains(v0sBestCollision, truthK0.globalIndex);
          const bool v0InTable = contains(v0sInTable, truthK0.globalIndex);
          const bool v0Final = contains(v0sFinal, truthK0.globalIndex);

          histos.fill(HIST("PairLossK0/Matching/hNMatches"), 0.0f, numberOfMatches(tracksAnyCollision, truthTrigger.globalIndex));
          histos.fill(HIST("PairLossK0/Matching/hNMatches"), 1.0f, numberOfMatches(tracksBestCollision, truthTrigger.globalIndex));
          histos.fill(HIST("PairLossK0/Matching/hNMatches"), 2.0f, numberOfMatches(triggersFinal, truthTrigger.globalIndex));
          histos.fill(HIST("PairLossK0/Matching/hNMatches"), 3.0f, numberOfMatches(v0sAnyCollision, truthK0.globalIndex));
          histos.fill(HIST("PairLossK0/Matching/hNMatches"), 4.0f, numberOfMatches(v0sBestCollision, truthK0.globalIndex));
          histos.fill(HIST("PairLossK0/Matching/hNMatches"), 5.0f, numberOfMatches(v0sInTable, truthK0.globalIndex));
          histos.fill(HIST("PairLossK0/Matching/hNMatches"), 6.0f, numberOfMatches(v0sFinal, truthK0.globalIndex));

          bool pairBeforeAutocorrelation = false;
          bool rejectedByAutocorrelation = false;
          std::vector<std::pair<PairLossTrackInfo, PairLossV0Info>> finalPairCombinations;
          if (triggerFinal && v0Final) {
            for (auto const& reconstructedTrigger : triggersFinal.at(truthTrigger.globalIndex)) {
              for (auto const& reconstructedV0 : v0sFinal.at(truthK0.globalIndex)) {
                if (reconstructedTrigger.collisionId != reconstructedV0.collisionId) {
                  continue;
                }
                const float reconstructedDeltaPhi = computeDeltaPhi(reconstructedTrigger.phi, reconstructedV0.phi);
                float reconstructedDeltaEta = reconstructedTrigger.eta - reconstructedV0.eta;
                if (masterConfigurations.doMirroringInDelataEta) {
                  reconstructedDeltaEta = std::abs(reconstructedDeltaEta);
                }
                if (reconstructedDeltaPhi < axisRanges[0][0] || reconstructedDeltaPhi > axisRanges[0][1] || reconstructedDeltaEta < axisRanges[1][0] || reconstructedDeltaEta > axisRanges[1][1]) {
                  continue;
                }
                pairBeforeAutocorrelation = true;
                if (doAutocorrelationRejection && (reconstructedTrigger.globalIndex == reconstructedV0.positiveTrackId || reconstructedTrigger.globalIndex == reconstructedV0.negativeTrackId)) {
                  rejectedByAutocorrelation = true;
                  continue;
                }
                finalPairCombinations.emplace_back(reconstructedTrigger, reconstructedV0);
                if (pairLossK0Configurations.fillFinalPairOnce) {
                  break;
                }
              }
              if (pairLossK0Configurations.fillFinalPairOnce && !finalPairCombinations.empty()) {
                break;
              }
            }
          }
          const bool finalPair = !finalPairCombinations.empty();
          const int finalPairFills = static_cast<int>(finalPairCombinations.size());

          // Order must match PairLossK0Stage / PairLossK0StageNames one to one.
          // positiveDaughterBestCollision and negativeDaughterBestCollision are
          // intentionally absent: the per-daughter breakdown lives in
          // State/hDaughterTrackStateClose, and pairBeforeAutocorrelation in
          // Geometry/hAutocorrelationRejected.
          const std::array<bool, PairLossK0NStages> stagePassed = {
            true,
            truthK0.findable,
            triggerBestCollision,
            triggerInTable,
            triggerFinal,
            bothDaughtersBestCollision,
            v0BestCollision,
            v0InTable,
            v0Final,
            finalPair};
          for (int stage = 0; stage < PairLossK0NStages; ++stage) {
            if (!stagePassed[stage]) {
              continue;
            }
            const int stageFills = (stage == PairLossFinalPair && !pairLossK0Configurations.fillFinalPairOnce) ? finalPairFills : 1;
            for (int fill = 0; fill < stageFills; ++fill) {
              histos.fill(HIST("PairLossK0/Stage/hCounts"), stage);
              histos.fill(HIST("PairLossK0/Stage/hPhysics"), stage, truthDeltaPhi, truthDeltaEta, truthK0.pt, truthTrigger.pt, multiplicity);
              if (truthK0.findable) {
                histos.fill(HIST("PairLossK0/Stage/hCountsFindable"), stage);
                histos.fill(HIST("PairLossK0/Stage/hPhysicsFindable"), stage, truthDeltaPhi, truthDeltaEta, truthK0.pt, truthTrigger.pt, multiplicity);
                if (closestDeltaPhiStar.valid) {
                  histos.fill(HIST("PairLossK0/Stage/hClose"), stage, closestDeltaPhiStar.minAbs, closestDeltaEta, truthK0.pt, truthTrigger.pt, magneticFieldSign, closestChargeProduct);
                }
              }
            }
          }

          // Attribution of the stage3 -> stage4 gate. Only meaningful for triggers that already
          // exist in the best collision.
          //
          // Key convention: the Passed bin is decided by the real triggerInTable, not by the
          // replica. An earlier version ran the replica on bestTrack() and took that as the final
          // answer, which left the Passed count of this plot more than ten percent away from
          // stage 4 (the replica uses this task's triggerTracks* configurables, while the real
          // table is produced by hStrangeCorrelationFilter.cxx from its own configurables; the two
          // are kept aligned by hand through comments and diverge as soon as the json changes only
          // the filter).
          //
          // The present definition is aligned by construction:
          //   count(Passed)            == count(stage 4)
          //   count(Passed) + count(Failed bins) == count(stage 3)
          // The replica is used for attribution only once the real table has already said "not in";
          // the two Unexplained/Rejects bins capture configuration drift, and a non-zero entry in
          // either means this task's triggerTracks* have to be realigned with the filter.
          if (triggerBestCollision) {
            int failureReason = PairLossTriggerTracksPassed;
            if (triggerInTable) {
              // Really passed. Also check whether the replica agrees, to monitor configuration drift.
              auto const* bestCollisionTrigger = bestTrack(tracksBestCollision, truthTrigger.globalIndex);
              if (bestCollisionTrigger != nullptr &&
                  classifyTriggerTracksFailure(*bestCollisionTrigger) != PairLossTriggerTracksPassed) {
                failureReason = PairLossTriggerTracksInTableButCutRejects;
              }
            } else {
              // Really failed. Attribute using the candidate track in the best collision that got
              // "furthest": classifyTriggerTracksFailure returns following the early-return order of
              // isValidTrigger, so a larger return value means the track was rejected later, and the
              // maximum therefore corresponds to the candidate closest to passing. The semantics are
              // "even the best candidate fell at this cut".
              failureReason = PairLossTriggerTracksNotInTableUnexplained;
              auto const iterator = tracksBestCollision.find(truthTrigger.globalIndex);
              if (iterator != tracksBestCollision.end()) {
                int furthest = -1;
                for (auto const& candidate : iterator->second) {
                  const int reason = classifyTriggerTracksFailure(candidate);
                  if (reason == PairLossTriggerTracksPassed) {
                    // The replica says this track should be in the table, but it is not:
                    // configuration drift.
                    furthest = -1;
                    break;
                  }
                  furthest = std::max(furthest, reason);
                }
                if (furthest > 0) {
                  failureReason = furthest;
                }
              }
            }
            histos.fill(HIST("PairLossK0/Stage/hTriggerTracksFailureReason"), failureReason, truthDeltaPhi, truthDeltaEta, truthK0.pt, truthTrigger.pt);
          }

          const int finalObjectState = (triggerFinal ? 2 : 0) + (v0Final ? 1 : 0);
          const int trackLevelState = (triggerBestCollision ? 2 : 0) + (bothDaughtersBestCollision ? 1 : 0);
          const int daughterTrackState = (positiveDaughterBestCollision ? 2 : 0) + (negativeDaughterBestCollision ? 1 : 0);
          histos.fill(HIST("PairLossK0/State/hFinalObjectStatePhysics"), finalObjectState, truthDeltaPhi, truthDeltaEta, truthK0.pt, truthTrigger.pt, magneticFieldSign);
          if (truthK0.findable && closestDeltaPhiStar.valid) {
            histos.fill(HIST("PairLossK0/State/hFinalObjectStateClose"), finalObjectState, closestDeltaPhiStar.minAbs, closestDeltaEta, truthK0.pt, truthTrigger.pt, magneticFieldSign, closestChargeProduct);
            histos.fill(HIST("PairLossK0/State/hTrackLevelStateClose"), trackLevelState, closestDeltaPhiStar.minAbs, closestDeltaEta, truthK0.pt, truthTrigger.pt, magneticFieldSign, closestChargeProduct);
            histos.fill(HIST("PairLossK0/State/hDaughterTrackStateClose"), daughterTrackState, closestDeltaPhiStar.minAbs, closestDeltaEta, truthK0.pt, truthTrigger.pt, magneticFieldSign, closestChargeProduct);

            if (auto const* matchedTrigger = bestTrack(tracksBestCollision, truthTrigger.globalIndex)) {
              histos.fill(HIST("PairLossK0/TrackQA/hTriggerTPCPairWeighted"), closestDeltaPhiStar.minAbs, matchedTrigger->tpcCrossedRows, matchedTrigger->tpcSharedClusters);
              histos.fill(HIST("PairLossK0/TrackQA/hTriggerITSPairWeighted"), closestDeltaPhiStar.minAbs, matchedTrigger->itsClusters, matchedTrigger->hasLayer0);
            }
            if (auto const* matchedPositiveDaughter = bestTrack(tracksBestCollision, truthK0.positiveDaughter.globalIndex)) {
              histos.fill(HIST("PairLossK0/TrackQA/hPositiveDaughterTPCPairWeighted"), closestDeltaPhiStar.minAbs, matchedPositiveDaughter->tpcCrossedRows, matchedPositiveDaughter->tpcSharedClusters);
              histos.fill(HIST("PairLossK0/TrackQA/hPositiveDaughterITSPairWeighted"), closestDeltaPhiStar.minAbs, matchedPositiveDaughter->itsClusters, matchedPositiveDaughter->hasLayer0);
            }
            if (auto const* matchedNegativeDaughter = bestTrack(tracksBestCollision, truthK0.negativeDaughter.globalIndex)) {
              histos.fill(HIST("PairLossK0/TrackQA/hNegativeDaughterTPCPairWeighted"), closestDeltaPhiStar.minAbs, matchedNegativeDaughter->tpcCrossedRows, matchedNegativeDaughter->tpcSharedClusters);
              histos.fill(HIST("PairLossK0/TrackQA/hNegativeDaughterITSPairWeighted"), closestDeltaPhiStar.minAbs, matchedNegativeDaughter->itsClusters, matchedNegativeDaughter->hasLayer0);
            }
            if (auto const* matchedV0 = bestV0(v0sBestCollision, truthK0.globalIndex)) {
              histos.fill(HIST("PairLossK0/V0QA/hCosPAVsMinDeltaPhiStar"), closestDeltaPhiStar.minAbs, matchedV0->cosPA);
              histos.fill(HIST("PairLossK0/V0QA/hDcaDaughtersVsMinDeltaPhiStar"), closestDeltaPhiStar.minAbs, matchedV0->dcaDaughters);
              histos.fill(HIST("PairLossK0/V0QA/hRadiusVsMinDeltaPhiStar"), closestDeltaPhiStar.minAbs, matchedV0->radius);
            }
            if (auto const* matchedAssocV0 = bestV0(v0sInTable, truthK0.globalIndex)) {
              histos.fill(HIST("PairLossK0/V0QA/hMassNSigmaVsMinDeltaPhiStar"), closestDeltaPhiStar.minAbs, matchedAssocV0->massNSigma, truthK0.pt);
            }
            if (pairBeforeAutocorrelation && !finalPair && rejectedByAutocorrelation) {
              histos.fill(HIST("PairLossK0/Geometry/hAutocorrelationRejected"), closestDeltaPhiStar.minAbs, closestDeltaEta, truthK0.pt, truthTrigger.pt, magneticFieldSign);
            }
          }

          for (auto const& [selectedTrigger, selectedV0] : finalPairCombinations) {
            histos.fill(HIST("PairLossK0/Response/hTriggerPt"), truthTrigger.pt, selectedTrigger.pt);
            histos.fill(HIST("PairLossK0/Response/hK0Pt"), truthK0.pt, selectedV0.pt);
            histos.fill(HIST("PairLossK0/Response/hDeltaPhi"), truthDeltaPhi, computeDeltaPhi(selectedTrigger.phi, selectedV0.phi));
            float reconstructedDeltaEta = selectedTrigger.eta - selectedV0.eta;
            if (masterConfigurations.doMirroringInDelataEta) {
              reconstructedDeltaEta = std::abs(reconstructedDeltaEta);
            }
            histos.fill(HIST("PairLossK0/Response/hDeltaEta"), truthDeltaEta, reconstructedDeltaEta);
            if (truthK0.findable) {
              if (positiveDeltaPhiStar.valid) {
                histos.fill(HIST("PairLossK0/Geometry/hFinalDaughters"), positiveDeltaPhiStar.signedAtMin, truthTrigger.eta - truthK0.positiveDaughter.eta, truthK0.pt, truthTrigger.pt, magneticFieldSign, truthTrigger.sign * truthK0.positiveDaughter.sign);
              }
              if (negativeDeltaPhiStar.valid) {
                histos.fill(HIST("PairLossK0/Geometry/hFinalDaughters"), negativeDeltaPhiStar.signedAtMin, truthTrigger.eta - truthK0.negativeDaughter.eta, truthK0.pt, truthTrigger.pt, magneticFieldSign, truthTrigger.sign * truthK0.negativeDaughter.sign);
              }
              if (closestDeltaPhiStar.valid) {
                histos.fill(HIST("PairLossK0/Geometry/hFinalMinDeltaPhiStarVsDecayRadius"), closestDeltaPhiStar.minAbs, truthK0.decayRadius, truthK0.pt);
                histos.fill(HIST("PairLossK0/Geometry/hFinalRadiusAtMinimum"), closestDeltaPhiStar.radiusAtMin, closestDeltaPhiStar.minAbs, truthK0.pt);
              }
            }
          }
        }
      }
      if (hasTruthPair) {
        histos.fill(HIST("PairLossK0/Event/hCounter"), 5.0f);
      }
      return;
    }
  }

  void processClosureTest(aod::McCollision const& /*mcCollision*/,
                          soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::PVMults>> const& recCollisions,
                          aod::McParticles const& mcParticles,
                          aod::V0MCCores const& v0MCCores,
                          V0DatasWithoutTrackX const& v0Candidates,
                          TracksCompleteMC const& tracks,
                          aod::TriggerTracks const& triggerTracks,
                          aod::AssocV0s const& associatedV0s)
  {

    // Reproduce the first processPairLossK0MC stages without changing the
    // selections or output of the pre-existing closure-test analysis below. All
    // histograms remain anchored to the same truth h-K0 pair; the folders differ
    // only in which reconstruction requirement is imposed:
    //   AnyTrack     the truth trigger must have at least one reconstructed track
    //                pointing back to it through its MC label, in any
    //                reconstructed collision associated with this MC collision
    //   AnyTrackK0   the truth K0 must have at least one reconstructed V0
    //                candidate pointing back to it, in any associated collision
    //   AnyTrackBoth both requirements at the same time
    //   Final        the strictest stage: both the trigger and the K0 must have a
    //                *fully selected* reconstructed counterpart -- the very
    //                TriggerTracks and AssocV0s entries, passing the very final
    //                selections, that the reconstructed correlation is built from --
    //                and both in one and the same reconstructed collision
    auto fillPairLossK0TruthAndAnyTrack = [&]() {
      if (!pairLossK0Configurations.doClosureTestStages) {
        return;
      }
      if (recCollisions.size() == 0) {
        return;
      }

      int64_t pairLossBestCollisionId = -1;
      int pairLossLargestNContributors = -1;
      for (auto const& collision : recCollisions) {
        if (collision.numContrib() > pairLossLargestNContributors) {
          pairLossLargestNContributors = collision.numContrib();
          pairLossBestCollisionId = collision.globalIndex();
        }
      }
      if (pairLossBestCollisionId < 0) {
        return;
      }

      std::unordered_set<int64_t> pairLossAnyTrackMcParticleIds;
      std::unordered_set<int64_t> pairLossAnyV0McParticleIds;
      for (auto const& collision : recCollisions) {
        const auto trackSlice = tracks.sliceBy(pairLossTracksPerCollision, collision.globalIndex());
        for (auto const& track : trackSlice) {
          if (track.has_mcParticle()) {
            pairLossAnyTrackMcParticleIds.insert(track.mcParticleId());
          }
        }
        // Same V0-side bookkeeping as tracksAnyCollision/v0sAnyCollision in
        // processPairLossK0MC: any V0 candidate in any associated collision whose
        // MC core is a true K0 short, with no candidate-quality selection at all.
        const auto v0Slice = v0Candidates.sliceBy(pairLossV0sPerCollision, collision.globalIndex());
        for (auto const& v0 : v0Slice) {
          const auto v0MC = v0MCCores.iteratorAt(v0.globalIndex());
          if (v0MC.particleIdMC() < 0 || v0MC.pdgCode() != PDG_t::kK0Short) {
            continue;
          }
          pairLossAnyV0McParticleIds.insert(v0MC.particleIdMC());
        }
      }

      // Final-stage bookkeeping. A truth object counts here only if it has at least
      // one fully selected reconstructed counterpart, where "fully selected" means
      // the exact same conditions the reconstructed correlation imposes in
      // fillCorrelationsV0(): the object is in the TriggerTracks / AssocV0s table
      // to begin with, and it then passes isValidTrigger() resp. the whole V0
      // selection chain including the dE/dx compatibility bit and the peak mass
      // window. The selections are read from the filter tables instead of being
      // recomputed from raw tracks so that this stage cannot drift away from the
      // reconstructed analysis it exists to be compared against.
      //
      // The stage is restricted to the best reconstructed collision and keeps the
      // reconstructed kinematics of every matched object, so that the pair test below
      // can also impose the reconstructed-coordinate angular range and the
      // reconstructed autocorrelation rejection. Both are conditions the reconstructed
      // same-event correlation imposes as well, so this Final stage is by construction
      // the same object as PairLossK0/Comparison/Final: the truth pairs that really do
      // end up in the reconstructed correlation.
      PairLossTrackMap pairLossFinalTriggers;
      PairLossV0Map pairLossFinalV0s;
      for (auto const& collision : recCollisions) {
        if (static_cast<int64_t>(collision.globalIndex()) != pairLossBestCollisionId) {
          continue;
        }

        const auto finalTriggerSlice = triggerTracks.sliceBy(collisionSliceTracks, collision.globalIndex());
        for (auto const& triggerEntry : finalTriggerSlice) {
          auto track = triggerEntry.track_as<TracksCompleteMC>();
          if (!track.has_mcParticle()) {
            continue;
          }
          if (!isValidTrigger(track, triggerEntry.isLeading())) {
            continue;
          }
          if (trackSelection.checkForITSTPCMissmatchMC && bitcheck(triggerEntry.mcMask(), 13)) {
            continue;
          }
          if (masterConfigurations.doTriggPhysicalPrimary && !triggerEntry.mcPhysicalPrimary()) {
            continue;
          }
          pairLossFinalTriggers[track.mcParticleId()].push_back(makePairLossTrackInfo(track));
        }

        const auto finalV0Slice = associatedV0s.sliceBy(collisionSliceV0s, collision.globalIndex());
        for (auto const& assocEntry : finalV0Slice) {
          auto v0 = assocEntry.v0Core_as<V0DatasWithoutTrackX>();
          const auto v0MC = v0MCCores.iteratorAt(v0.globalIndex());
          if (v0MC.particleIdMC() < 0 || v0MC.pdgCode() != PDG_t::kK0Short || !assocEntry.mcTrue(IndexK0)) {
            continue;
          }
          auto positiveTrack = v0.posTrack_as<TracksCompleteMC>();
          auto negativeTrack = v0.negTrack_as<TracksCompleteMC>();
          bool passesFinalSelection = true;
          if (masterConfigurations.doPPAnalysis) {
            passesFinalSelection = v0.v0radius() >= v0Selection.v0RadiusMin && v0.v0radius() <= v0Selection.v0RadiusMax &&
                                   std::abs(v0.dcapostopv()) >= v0Selection.dcapostopv && std::abs(v0.dcanegtopv()) >= v0Selection.dcanegtopv &&
                                   v0.v0cosPA() >= v0Selection.v0cospa && v0.dcaV0daughters() <= v0Selection.dcaV0dau;
          } else {
            const float dcaCut = v0Selection.dcaDaugToPVForK0s == 0.0f ? v0Selection.dcaMesonToPV : v0Selection.dcaDaugToPVForK0s;
            const bool passesLifetime = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) *
                                          o2::constants::physics::MassK0Short <
                                        v0Selection.lifetimecutK0S;
            const bool passesDaughterDcaAndArmenteros = std::abs(v0.dcapostopv()) > dcaCut && std::abs(v0.dcanegtopv()) > dcaCut &&
                                                        v0.qtarm() * v0Selection.armPodCut > std::abs(v0.alpha());
            passesFinalSelection = v0SelectedPbPb(v0) && passesLifetime && passesDaughterDcaAndArmenteros;
          }
          passesFinalSelection = passesFinalSelection && positiveTrack.tpcNClsCrossedRows() >= trackSelection.minTPCNCrossedRowsAssociated &&
                                 negativeTrack.tpcNClsCrossedRows() >= trackSelection.minTPCNCrossedRowsAssociated;
          if (trackSelection.checksRequireTPCChi2) {
            passesFinalSelection = passesFinalSelection && positiveTrack.tpcChi2NCl() >= trackSelection.minTPCChi2PerClusterAssociated &&
                                   negativeTrack.tpcChi2NCl() >= trackSelection.minTPCChi2PerClusterAssociated;
          }
          if (trackSelection.requireClusterInITS) {
            passesFinalSelection = passesFinalSelection && positiveTrack.itsNCls() >= trackSelection.minITSClustersForDaughterTracks &&
                                   negativeTrack.itsNCls() >= trackSelection.minITSClustersForDaughterTracks;
          }
          passesFinalSelection = passesFinalSelection && assocEntry.compatible(IndexK0, trackSelection.dEdxCompatibility) &&
                                 (!doAssocPhysicalPrimary || assocEntry.mcPhysicalPrimary()) &&
                                 assocEntry.invMassNSigma(IndexK0) > -massWindowConfigurations.maxPeakNSigma &&
                                 assocEntry.invMassNSigma(IndexK0) < massWindowConfigurations.maxPeakNSigma &&
                                 v0.pt() >= axisRanges[2][0] && v0.pt() <= axisRanges[2][1];
          if (!passesFinalSelection) {
            continue;
          }
          pairLossFinalV0s[v0MC.particleIdMC()].push_back(PairLossV0Info{
            .globalIndex = static_cast<int64_t>(v0.globalIndex()),
            .positiveTrackId = static_cast<int64_t>(positiveTrack.globalIndex()),
            .negativeTrackId = static_cast<int64_t>(negativeTrack.globalIndex()),
            .pt = v0.pt(),
            .eta = v0.eta(),
            .phi = v0.phi(),
            .radius = v0.v0radius(),
            .cosPA = v0.v0cosPA(),
            .dcaDaughters = v0.dcaV0daughters(),
            .massNSigma = assocEntry.invMassNSigma(IndexK0)});
        }
      }

      // Object-level membership, used only for the single-particle spectra: the object
      // has a fully selected reconstructed counterpart in the best collision. The pair
      // histogram uses pairLossHasFinalPair() instead, which is stricter: it also
      // requires the reconstructed pair itself to fall in the reconstructed angular
      // range and to survive the reconstructed autocorrelation rejection.
      auto pairLossHasFinalTrigger = [&](int64_t mcId) {
        return pairLossFinalTriggers.find(mcId) != pairLossFinalTriggers.end();
      };
      auto pairLossHasFinalK0 = [&](int64_t mcId) {
        return pairLossFinalV0s.find(mcId) != pairLossFinalV0s.end();
      };
      auto pairLossHasFinalPair = [&](int64_t triggerMcId, int64_t k0McId) {
        const auto triggerMatches = pairLossFinalTriggers.find(triggerMcId);
        if (triggerMatches == pairLossFinalTriggers.end()) {
          return false;
        }
        const auto v0Matches = pairLossFinalV0s.find(k0McId);
        if (v0Matches == pairLossFinalV0s.end()) {
          return false;
        }
        for (auto const& reconstructedTrigger : triggerMatches->second) {
          for (auto const& reconstructedV0 : v0Matches->second) {
            float reconstructedDeltaEta = reconstructedTrigger.eta - reconstructedV0.eta;
            if (masterConfigurations.doMirroringInDelataEta) {
              reconstructedDeltaEta = std::abs(reconstructedDeltaEta);
            }
            const float reconstructedDeltaPhi = computeDeltaPhi(reconstructedTrigger.phi, reconstructedV0.phi);
            if (reconstructedDeltaPhi < axisRanges[0][0] || reconstructedDeltaPhi > axisRanges[0][1] ||
                reconstructedDeltaEta < axisRanges[1][0] || reconstructedDeltaEta > axisRanges[1][1]) {
              continue;
            }
            if (doAutocorrelationRejection && (reconstructedTrigger.globalIndex == reconstructedV0.positiveTrackId || reconstructedTrigger.globalIndex == reconstructedV0.negativeTrackId)) {
              continue;
            }
            return true;
          }
        }
        return false;
      };

      for (auto const& collision : recCollisions) {
        if (static_cast<int64_t>(collision.globalIndex()) != pairLossBestCollisionId) {
          continue;
        }

        const bool collisionSelected = !pairLossK0Configurations.applyRecoEventSelection ||
                                       (masterConfigurations.doPPAnalysis ? isCollisionSelected(collision) : isCollisionSelectedPbPb(collision, false));
        if (!collisionSelected) {
          return;
        }

        const float pairLossBestCollisionVtxZ = collision.posZ();
        const float pairLossBestCollisionMultiplicity = masterConfigurations.doPPAnalysis ? collision.centFT0M() : collision.centFT0C();
        std::vector<PairLossTruthTrackInfo> pairLossTruthTriggers;
        std::vector<PairLossTruthK0Info> pairLossTruthK0s;

        for (auto const& mcParticle : mcParticles) {
          if (isPairLossTriggerPdg(mcParticle.pdgCode()) && std::abs(mcParticle.eta()) <= etaSel &&
              mcParticle.pt() >= axisRanges[3][0] && mcParticle.pt() <= axisRanges[3][1] &&
              (!masterConfigurations.doTriggPhysicalPrimary || mcParticle.isPhysicalPrimary())) {
            auto const* pdgParticle = pdgDB->GetParticle(mcParticle.pdgCode());
            const double charge = pdgParticle != nullptr ? pdgParticle->Charge() : 0.0;
            const int sign = charge > 0.0 ? 1 : (charge < 0.0 ? -1 : 0);
            if (!((triggerParticleCharge > 0 && sign < 0) || (triggerParticleCharge < 0 && sign > 0) || sign == 0)) {
              pairLossTruthTriggers.push_back(PairLossTruthTrackInfo{
                .globalIndex = static_cast<int64_t>(mcParticle.globalIndex()),
                .pt = mcParticle.pt(),
                .eta = mcParticle.eta(),
                .phi = mcParticle.phi(),
                .sign = sign});
            }
          }

          if (mcParticle.pdgCode() != PDG_t::kK0Short || std::abs(mcParticle.eta()) > etaSel ||
              mcParticle.pt() < axisRanges[2][0] || mcParticle.pt() > axisRanges[2][1] ||
              (doAssocPhysicalPrimary && !mcParticle.isPhysicalPrimary())) {
            continue;
          }

          PairLossTruthK0Info truthK0{
            .globalIndex = static_cast<int64_t>(mcParticle.globalIndex()),
            .pt = mcParticle.pt(),
            .eta = mcParticle.eta(),
            .phi = mcParticle.phi(),
            .decayRadius = -1.0f,
            .findable = false,
            .positiveDaughter = {},
            .negativeDaughter = {}};
          for (auto const& daughter : mcParticle.daughters_as<aod::McParticles>()) {
            PairLossTruthTrackInfo daughterInfo{
              .globalIndex = static_cast<int64_t>(daughter.globalIndex()),
              .pt = daughter.pt(),
              .eta = daughter.eta(),
              .phi = daughter.phi(),
              .sign = daughter.pdgCode() > 0 ? 1 : -1};
            if (daughter.pdgCode() == PDG_t::kPiPlus) {
              truthK0.positiveDaughter = daughterInfo;
            } else if (daughter.pdgCode() == -PDG_t::kPiPlus) {
              truthK0.negativeDaughter = daughterInfo;
            }
          }
          pairLossTruthK0s.push_back(truthK0);
        }

        if (useTheLeadingParticleAsTrigger && pairLossTruthTriggers.size() > 1) {
          const auto leadingTrigger = std::max_element(pairLossTruthTriggers.begin(), pairLossTruthTriggers.end(), [](auto const& lhs, auto const& rhs) { return lhs.pt < rhs.pt; });
          const auto leadingTriggerCopy = *leadingTrigger;
          pairLossTruthTriggers.clear();
          pairLossTruthTriggers.push_back(leadingTriggerCopy);
        }

        // One entry per object per folder, at the level that folder prescribes:
        // the trigger is at truth level in Truth/ and AnyTrackK0/, at any level in
        // AnyTrack/ and AnyTrackBoth/, at fully-selected level in Final/; the K0 is
        // at truth level in Truth/ and AnyTrack/, at any level in AnyTrackK0/ and
        // AnyTrackBoth/, at fully-selected level in Final/.
        for (auto const& truthTrigger : pairLossTruthTriggers) {
          histos.fill(HIST("ClosureTest/PairLossK0/Truth/hTrigger"), truthTrigger.pt, truthTrigger.eta, truthTrigger.phi);
          histos.fill(HIST("ClosureTest/PairLossK0/AnyTrackK0/hTrigger"), truthTrigger.pt, truthTrigger.eta, truthTrigger.phi);
          if (pairLossAnyTrackMcParticleIds.find(truthTrigger.globalIndex) != pairLossAnyTrackMcParticleIds.end()) {
            histos.fill(HIST("ClosureTest/PairLossK0/AnyTrack/hTrigger"), truthTrigger.pt, truthTrigger.eta, truthTrigger.phi);
            histos.fill(HIST("ClosureTest/PairLossK0/AnyTrackBoth/hTrigger"), truthTrigger.pt, truthTrigger.eta, truthTrigger.phi);
          }
          if (pairLossHasFinalTrigger(truthTrigger.globalIndex)) {
            histos.fill(HIST("ClosureTest/PairLossK0/Final/hTrigger"), truthTrigger.pt, truthTrigger.eta, truthTrigger.phi);
          }
        }
        for (auto const& truthK0 : pairLossTruthK0s) {
          histos.fill(HIST("ClosureTest/PairLossK0/Truth/hK0Short"), truthK0.pt, truthK0.eta, truthK0.phi);
          histos.fill(HIST("ClosureTest/PairLossK0/AnyTrack/hK0Short"), truthK0.pt, truthK0.eta, truthK0.phi);
          if (pairLossAnyV0McParticleIds.find(truthK0.globalIndex) != pairLossAnyV0McParticleIds.end()) {
            histos.fill(HIST("ClosureTest/PairLossK0/AnyTrackK0/hK0Short"), truthK0.pt, truthK0.eta, truthK0.phi);
            histos.fill(HIST("ClosureTest/PairLossK0/AnyTrackBoth/hK0Short"), truthK0.pt, truthK0.eta, truthK0.phi);
          }
          if (pairLossHasFinalK0(truthK0.globalIndex)) {
            histos.fill(HIST("ClosureTest/PairLossK0/Final/hK0Short"), truthK0.pt, truthK0.eta, truthK0.phi);
          }
        }

        for (auto const& truthTrigger : pairLossTruthTriggers) {
          const bool triggerHasAnyTrack = pairLossAnyTrackMcParticleIds.find(truthTrigger.globalIndex) != pairLossAnyTrackMcParticleIds.end();
          for (auto const& truthK0 : pairLossTruthK0s) {
            if (truthTrigger.globalIndex == truthK0.positiveDaughter.globalIndex || truthTrigger.globalIndex == truthK0.negativeDaughter.globalIndex) {
              continue;
            }
            const float truthDeltaPhi = computeDeltaPhi(truthTrigger.phi, truthK0.phi);
            float truthDeltaEta = truthTrigger.eta - truthK0.eta;
            if (masterConfigurations.doMirroringInDelataEta) {
              truthDeltaEta = std::abs(truthDeltaEta);
            }
            if (truthDeltaPhi < axisRanges[0][0] || truthDeltaPhi > axisRanges[0][1] ||
                truthDeltaEta < axisRanges[1][0] || truthDeltaEta > axisRanges[1][1]) {
              continue;
            }

            const bool k0HasAnyV0 = pairLossAnyV0McParticleIds.find(truthK0.globalIndex) != pairLossAnyV0McParticleIds.end();

            histos.fill(HIST("ClosureTest/PairLossK0/Truth/sameEvent/K0Short"), truthDeltaPhi, truthDeltaEta, truthK0.pt, truthTrigger.pt, pairLossBestCollisionVtxZ, pairLossBestCollisionMultiplicity);
            if (triggerHasAnyTrack) {
              histos.fill(HIST("ClosureTest/PairLossK0/AnyTrack/sameEvent/K0Short"), truthDeltaPhi, truthDeltaEta, truthK0.pt, truthTrigger.pt, pairLossBestCollisionVtxZ, pairLossBestCollisionMultiplicity);
            }
            if (k0HasAnyV0) {
              histos.fill(HIST("ClosureTest/PairLossK0/AnyTrackK0/sameEvent/K0Short"), truthDeltaPhi, truthDeltaEta, truthK0.pt, truthTrigger.pt, pairLossBestCollisionVtxZ, pairLossBestCollisionMultiplicity);
            }
            if (triggerHasAnyTrack && k0HasAnyV0) {
              histos.fill(HIST("ClosureTest/PairLossK0/AnyTrackBoth/sameEvent/K0Short"), truthDeltaPhi, truthDeltaEta, truthK0.pt, truthTrigger.pt, pairLossBestCollisionVtxZ, pairLossBestCollisionMultiplicity);
            }
            if (pairLossHasFinalPair(truthTrigger.globalIndex, truthK0.globalIndex)) {
              histos.fill(HIST("ClosureTest/PairLossK0/Final/sameEvent/K0Short"), truthDeltaPhi, truthDeltaEta, truthK0.pt, truthTrigger.pt, pairLossBestCollisionVtxZ, pairLossBestCollisionMultiplicity);
            }
          }
        }
        return;
      }
    };
    fillPairLossK0TruthAndAnyTrack();

    std::vector<uint32_t> triggerIndices;
    std::vector<std::vector<uint32_t>> associatedIndices;
    std::vector<uint32_t> assocHadronIndices;
    std::vector<uint32_t> piIndices;
    std::vector<uint32_t> k0ShortIndices;
    std::vector<uint32_t> lambdaIndices;
    std::vector<uint32_t> antiLambdaIndices;
    std::vector<uint32_t> xiMinusIndices;
    std::vector<uint32_t> xiPlusIndices;
    std::vector<uint32_t> omegaMinusIndices;
    std::vector<uint32_t> omegaPlusIndices;

    // A truth trigger belongs to this set when at least one reconstructed track
    // in any reconstructed collision associated with the current MC collision
    // points back to it through its MC label. No best-collision, track-quality,
    // or final-trigger selection is imposed here: this is the "any track" stage.
    std::unordered_set<int64_t> mcParticleIdsWithRecoTrackMatch;
    if (masterConfigurations.doClosureTestTriggerWithRecoTrackMatch) {
      for (auto const& recCollision : recCollisions) {
        const auto trackSlice = tracks.sliceBy(pairLossTracksPerCollision, recCollision.globalIndex());
        for (auto const& track : trackSlice) {
          if (track.has_mcParticle()) {
            mcParticleIdsWithRecoTrackMatch.insert(track.mcParticleId());
          }
        }
      }
    }

    for (auto const& mcParticle : mcParticles) {
      double geta = mcParticle.eta();
      if (std::abs(geta) > etaSel) {
        continue;
      }
      double gpt = mcParticle.pt();
      if (std::abs(mcParticle.pdgCode()) == PDG_t::kPiPlus || std::abs(mcParticle.pdgCode()) == PDG_t::kKPlus || std::abs(mcParticle.pdgCode()) == PDG_t::kProton || std::abs(mcParticle.pdgCode()) == PDG_t::kElectron || std::abs(mcParticle.pdgCode()) == PDG_t::kMuonMinus) {
        if (!masterConfigurations.doTriggPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
          histos.fill(HIST("hClosureQAPtTrigger"), gpt, 0.0f); // step 1: no event selection whatsoever
        }
      }

      if (!doAssocPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
        if (mcParticle.pdgCode() == PDG_t::kK0Short && masterConfigurations.doCorrelationK0Short) {
          histos.fill(HIST("hClosureQAPtAssociatedK0"), gpt, 0.0f); // step 1: no event selection whatsoever
        }
      }
    }

    histos.fill(HIST("hClosureTestEventCounter"), 0.5f);

    float bestCollisionCentpercentile = -1;
    float bestCollisionVtxZ = 0.0f;
    bool bestCollisionSel8 = false;
    bool bestCollisionINELgtZERO = false;
    bool bestCollisionINELgtONE = false;
    bool bestCollisionNoSameBunchPileup = false;
    bool bestCollisionGoodTriggerTVX = false;
    bool bestCollisionGoodZvtxFT0vsPV = false;
    bool isCollisionSelect = false;
    int biggestNContribs = -1;
    uint32_t bestCollisionTriggerPresenceMap = 0;

    for (auto const& recCollision : recCollisions) {
      if (biggestNContribs < recCollision.numContrib()) {
        biggestNContribs = recCollision.numContrib();
        bestCollisionCentpercentile = masterConfigurations.doPPAnalysis ? recCollision.centFT0M() : recCollision.centFT0C();
        if (masterConfigurations.applyNewMCSelection) {
          isCollisionSelect = ((masterConfigurations.doPPAnalysis && isCollisionSelected(recCollision)) || (!masterConfigurations.doPPAnalysis && isCollisionSelectedPbPb(recCollision, false)));
        } else {
          bestCollisionSel8 = recCollision.sel8();
          bestCollisionVtxZ = recCollision.posZ();
          bestCollisionINELgtZERO = recCollision.isInelGt0();
          bestCollisionINELgtONE = recCollision.isInelGt1();
          bestCollisionNoSameBunchPileup = recCollision.selection_bit(o2::aod::evsel::kNoSameBunchPileup);
          bestCollisionGoodTriggerTVX = recCollision.selection_bit(aod::evsel::kIsTriggerTVX);
          bestCollisionGoodZvtxFT0vsPV = recCollision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV);
        }
        if (triggerPresenceMap.size() > 0) {
          bestCollisionTriggerPresenceMap = triggerPresenceMap[recCollision.globalIndex()];
        }
      }
    }
    // ________________________________________________
    // skip if desired trigger not found
    if (triggerPresenceMap.size() > 0 && !TESTBIT(bestCollisionTriggerPresenceMap, triggerBinToSelect)) {
      return;
    }

    if (masterConfigurations.applyNewMCSelection) {
      if (!isCollisionSelect) {
        return;
      }
    } else {
      if (masterConfigurations.doGenEventSelection) {
        if (!bestCollisionSel8) {
          return;
        }
        if (std::abs(bestCollisionVtxZ) > masterConfigurations.zVertexCut) {
          return;
        }
        if (!bestCollisionINELgtZERO) {
          return;
        }
        if (masterConfigurations.selectINELgtONE && !bestCollisionINELgtONE) {
          return;
        }
        if (masterConfigurations.rejectSameBunchPileup && !bestCollisionNoSameBunchPileup) {
          return;
        }
        if (masterConfigurations.requireGoodTriggerTVX && !bestCollisionGoodTriggerTVX) {
          return;
        }
        if (masterConfigurations.requireGoodZvtxFT0vsPV && !bestCollisionGoodZvtxFT0vsPV) {
          return;
        }
        if (bestCollisionCentpercentile > axisRanges[5][1] || bestCollisionCentpercentile < axisRanges[5][0]) {
          return;
        }
      }
    }
    histos.fill(HIST("hClosureTestEventCounter"), 1.5f);

    for (auto const& mcParticle : mcParticles) {
      double geta = mcParticle.eta();
      if (std::abs(geta) > etaSel) {
        continue;
      }
      double gpt = mcParticle.pt();
      if (std::abs(mcParticle.pdgCode()) == PDG_t::kPiPlus || std::abs(mcParticle.pdgCode()) == PDG_t::kKPlus || std::abs(mcParticle.pdgCode()) == PDG_t::kProton || std::abs(mcParticle.pdgCode()) == PDG_t::kElectron || std::abs(mcParticle.pdgCode()) == PDG_t::kMuonMinus) {
        if (!masterConfigurations.doTriggPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
          histos.fill(HIST("hClosureQAPtTrigger"), gpt, 1.0f); // step 2: after event selection
        }
      }

      if (!doAssocPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
        if (mcParticle.pdgCode() == PDG_t::kK0Short && masterConfigurations.doCorrelationK0Short) {
          histos.fill(HIST("hClosureQAPtAssociatedK0"), gpt, 1.0f); // step 2: after event selection
        }
      }
    }

    int iteratorNum = -1;
    for (auto const& mcParticle : mcParticles) {
      iteratorNum = iteratorNum + 1;
      double geta = mcParticle.eta();
      double gpt = mcParticle.pt();
      double gphi = mcParticle.phi();
      if (std::abs(geta) > etaSel) {
        continue;
      }
      if (std::abs(mcParticle.pdgCode()) == PDG_t::kPiPlus || std::abs(mcParticle.pdgCode()) == PDG_t::kKPlus || std::abs(mcParticle.pdgCode()) == PDG_t::kProton || std::abs(mcParticle.pdgCode()) == PDG_t::kElectron || std::abs(mcParticle.pdgCode()) == PDG_t::kMuonMinus) {
        if (!masterConfigurations.doTriggPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
          triggerIndices.emplace_back(iteratorNum);
          histos.fill(HIST("ClosureTest/hTrigger"), gpt, geta, bestCollisionCentpercentile);
          if (masterConfigurations.doClosureTestTriggerWithRecoTrackMatch && mcParticleIdsWithRecoTrackMatch.find(mcParticle.globalIndex()) != mcParticleIdsWithRecoTrackMatch.end()) {
            histos.fill(HIST("ClosureTest/TriggerWithRecoTrackMatch/hTrigger"), gpt, geta, bestCollisionCentpercentile);
          }
        }
        if (masterConfigurations.doCorrelationHadron) {
          if (!doAssocPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
            assocHadronIndices.emplace_back(iteratorNum);
            histos.fill(HIST("ClosureTest/hHadron"), gpt, geta, gphi);
          }
        }
      }
      if (!doAssocPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
        if (std::abs(mcParticle.pdgCode()) == PDG_t::kPiPlus && masterConfigurations.doCorrelationPion) {
          piIndices.emplace_back(iteratorNum);
          histos.fill(HIST("ClosureTest/hPion"), gpt, geta, gphi);
        }
        if (mcParticle.pdgCode() == PDG_t::kK0Short && masterConfigurations.doCorrelationK0Short) {
          k0ShortIndices.emplace_back(iteratorNum);
          histos.fill(HIST("ClosureTest/hK0Short"), gpt, geta, gphi);
        }
        if (mcParticle.pdgCode() == PDG_t::kLambda0 && masterConfigurations.doCorrelationLambda) {
          lambdaIndices.emplace_back(iteratorNum);
          histos.fill(HIST("ClosureTest/hLambda"), gpt, geta, gphi);
        }
        if (mcParticle.pdgCode() == PDG_t::kLambda0Bar && masterConfigurations.doCorrelationAntiLambda) {
          antiLambdaIndices.emplace_back(iteratorNum);
          histos.fill(HIST("ClosureTest/hAntiLambda"), gpt, geta, gphi);
        }
        if (mcParticle.pdgCode() == PDG_t::kXiMinus && masterConfigurations.doCorrelationXiMinus) {
          xiMinusIndices.emplace_back(iteratorNum);
          histos.fill(HIST("ClosureTest/hXiMinus"), gpt, geta, gphi);
        }
        if (mcParticle.pdgCode() == PDG_t::kXiPlusBar && masterConfigurations.doCorrelationXiPlus) {
          xiPlusIndices.emplace_back(iteratorNum);
          histos.fill(HIST("ClosureTest/hXiPlus"), gpt, geta, gphi);
        }
        if (mcParticle.pdgCode() == PDG_t::kOmegaMinus && masterConfigurations.doCorrelationOmegaMinus) {
          omegaMinusIndices.emplace_back(iteratorNum);
          histos.fill(HIST("ClosureTest/hOmegaMinus"), gpt, geta, gphi);
        }
        if (mcParticle.pdgCode() == PDG_t::kOmegaPlusBar && masterConfigurations.doCorrelationOmegaPlus) {
          omegaPlusIndices.emplace_back(iteratorNum);
          histos.fill(HIST("ClosureTest/hOmegaPlus"), gpt, geta, gphi);
        }
      }
    }

    associatedIndices.emplace_back(k0ShortIndices);
    associatedIndices.emplace_back(lambdaIndices);
    associatedIndices.emplace_back(antiLambdaIndices);
    associatedIndices.emplace_back(xiMinusIndices);
    associatedIndices.emplace_back(xiPlusIndices);
    associatedIndices.emplace_back(omegaMinusIndices);
    associatedIndices.emplace_back(omegaPlusIndices);
    associatedIndices.emplace_back(piIndices);
    associatedIndices.emplace_back(assocHadronIndices);

    for (std::size_t iTrigger = 0; iTrigger < triggerIndices.size(); iTrigger++) {
      auto triggerParticle = mcParticles.iteratorAt(triggerIndices[iTrigger]);
      // check range of trigger particle
      if (triggerParticle.pt() > axisRanges[3][1] || triggerParticle.pt() < axisRanges[3][0]) {
        continue;
      }
      double getatrigger = triggerParticle.eta();
      double gphitrigger = triggerParticle.phi();
      double pttrigger = triggerParticle.pt();
      const bool triggerHasRecoTrackMatch = masterConfigurations.doClosureTestTriggerWithRecoTrackMatch && mcParticleIdsWithRecoTrackMatch.find(triggerParticle.globalIndex()) != mcParticleIdsWithRecoTrackMatch.end();
      auto const* triggerPdg = pdgDB->GetParticle(triggerParticle.pdgCode());
      const double triggerCharge = triggerPdg ? triggerPdg->Charge() : 0.;
      auto const& mother = triggerParticle.mothers_first_as<aod::McParticles>();
      auto globalIndex = mother.globalIndex();
      static_for<0, 8>([&](auto i) { // associated loop
        constexpr int Index = i.value;
        for (std::size_t iassoc = 0; iassoc < associatedIndices[Index].size(); iassoc++) {
          auto assocParticle = mcParticles.iteratorAt(associatedIndices[Index][iassoc]);
          if (triggerIndices[iTrigger] != associatedIndices[Index][iassoc] && globalIndex != assocParticle.globalIndex()) { // avoid self
            double getaassoc = assocParticle.eta();
            double gphiassoc = assocParticle.phi();
            double ptassoc = assocParticle.pt();
            double deltaphi = computeDeltaPhi(gphitrigger, gphiassoc);
            double deltaeta = getatrigger - getaassoc;

            // skip if basic ranges not met
            if (deltaphi < axisRanges[0][0] || deltaphi > axisRanges[0][1]) {
              continue;
            }
            if (deltaeta < axisRanges[1][0] || deltaeta > axisRanges[1][1]) {
              continue;
            }
            if (ptassoc < axisRanges[2][0] || ptassoc > axisRanges[2][1]) {
              continue;
            }
            if (TESTBIT(doCorrelation, i)) {
              histos.fill(HIST("ClosureTest/sameEvent/") + HIST(Particlenames[Index]), computeDeltaPhi(gphitrigger, gphiassoc), deltaeta, ptassoc, pttrigger, bestCollisionVtxZ, bestCollisionCentpercentile);
              if (triggerHasRecoTrackMatch) {
                histos.fill(HIST("ClosureTest/TriggerWithRecoTrackMatch/sameEvent/") + HIST(Particlenames[Index]), computeDeltaPhi(gphitrigger, gphiassoc), deltaeta, ptassoc, pttrigger, bestCollisionVtxZ, bestCollisionCentpercentile);
              }
            }
            if (i < 3 && TESTBIT(doCorrelation, i) && masterConfigurations.doCorrelationsHadronV0daughter) {
              auto const assocParticleDaughters = assocParticle.daughters_as<aod::McParticles>();
              for (const auto& assocParticleDaughter : assocParticleDaughters) {
                auto const* daughterPdg = pdgDB->GetParticle(assocParticleDaughter.pdgCode());
                const double daughterCharge = daughterPdg ? daughterPdg->Charge() : 0.;
                // Neutral daughters, including pi0, do not correspond to reconstructed V0 daughter tracks.
                if (triggerCharge == 0. || daughterCharge == 0.) {
                  continue;
                }
                if (triggerCharge * daughterCharge > 0.) {
                  histos.fill(HIST("ClosureTest/sameEvent/") + HIST(Particlenames[Index]) + HIST("_SameSignDaughter"), computeDeltaPhi(gphitrigger, assocParticleDaughter.phi()), getatrigger - assocParticleDaughter.eta(), assocParticleDaughter.pt(), pttrigger, bestCollisionVtxZ, bestCollisionCentpercentile);
                  if (triggerHasRecoTrackMatch) {
                    histos.fill(HIST("ClosureTest/TriggerWithRecoTrackMatch/sameEvent/") + HIST(Particlenames[Index]) + HIST("_SameSignDaughter"), computeDeltaPhi(gphitrigger, assocParticleDaughter.phi()), getatrigger - assocParticleDaughter.eta(), assocParticleDaughter.pt(), pttrigger, bestCollisionVtxZ, bestCollisionCentpercentile);
                  }
                } else {
                  histos.fill(HIST("ClosureTest/sameEvent/") + HIST(Particlenames[Index]) + HIST("_OppSignDaughter"), computeDeltaPhi(gphitrigger, assocParticleDaughter.phi()), getatrigger - assocParticleDaughter.eta(), assocParticleDaughter.pt(), pttrigger, bestCollisionVtxZ, bestCollisionCentpercentile);
                  if (triggerHasRecoTrackMatch) {
                    histos.fill(HIST("ClosureTest/TriggerWithRecoTrackMatch/sameEvent/") + HIST(Particlenames[Index]) + HIST("_OppSignDaughter"), computeDeltaPhi(gphitrigger, assocParticleDaughter.phi()), getatrigger - assocParticleDaughter.eta(), assocParticleDaughter.pt(), pttrigger, bestCollisionVtxZ, bestCollisionCentpercentile);
                  }
                }
              }
            }
          }
        }
      });
    }
  }

  void processFeedDown(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults>::iterator const& collision, aod::AssocV0s const& associatedV0s, aod::McParticles const&, V0DatasWithoutTrackXMC const&, TracksComplete const&, aod::BCsWithTimestamps const&)
  {

    // ________________________________________________
    // Perform basic event selection
    if (!isCollisionSelected(collision)) {
      return;
    }

    for (auto const& v0 : associatedV0s) {
      auto v0Data = v0.v0Core_as<V0DatasWithoutTrackXMC>();

      //---] track quality check [---
      auto postrack = v0Data.posTrack_as<TracksComplete>();
      auto negtrack = v0Data.negTrack_as<TracksComplete>();
      if (postrack.tpcNClsCrossedRows() < trackSelection.minTPCNCrossedRowsAssociated || negtrack.tpcNClsCrossedRows() < trackSelection.minTPCNCrossedRowsAssociated) {
        continue;
      }

      //---] syst cuts [---
      if (v0Data.v0radius() < v0Selection.v0RadiusMin || v0Data.v0radius() > v0Selection.v0RadiusMax ||
          std::abs(v0Data.dcapostopv()) < v0Selection.dcapostopv || std::abs(v0Data.dcanegtopv()) < v0Selection.dcanegtopv ||
          v0Data.v0cosPA() < v0Selection.v0cospa || v0Data.dcaV0daughters() > v0Selection.dcaV0dau) {
        continue;
      }

      if (v0Data.has_mcParticle()) {
        auto v0mcParticle = v0Data.mcParticle_as<aod::McParticles>();
        int mcParticlePdg = v0mcParticle.pdgCode();
        if (mcParticlePdg == PDG_t::kLambda0 && !v0mcParticle.isPhysicalPrimary()) {
          auto v0mothers = v0mcParticle.mothers_as<aod::McParticles>();
          if (v0mothers.size() == 1) {
            for (const auto& v0mcParticleMother : v0mothers) {
              // auto& v0mcParticleMother = v0mothers.front();
              if (std::abs(v0mcParticleMother.eta()) > etaSel) {
                continue;
              }
              if (v0mcParticleMother.pdgCode() == PDG_t::kXiMinus) // Xi Minus Mother Matched
              {
                histos.fill(HIST("hLambdaXiMinusFeeddownMatrix"), v0mcParticle.pt(), v0mcParticleMother.pt());
                histos.fill(HIST("hLambdaFromXiMinusEtaVsPtVsPhi"), v0mcParticle.pt(), v0mcParticle.eta(), v0mcParticle.phi());
              }
              if (v0mcParticleMother.pdgCode() == o2::constants::physics::Pdg::kXi0) // Xi Zero Mother Matched
              {
                histos.fill(HIST("hLambdaXiZeroFeeddownMatrix"), v0mcParticle.pt(), v0mcParticleMother.pt());
                histos.fill(HIST("hLambdaFromXiZeroEtaVsPtVsPhi"), v0mcParticle.pt(), v0mcParticle.eta(), v0mcParticle.phi());
              }
              if (v0mcParticleMother.pdgCode() == PDG_t::kOmegaMinus) // Omega Mother Matched
              {
                histos.fill(HIST("hLambdaOmegaFeeddownMatrix"), v0mcParticle.pt(), v0mcParticleMother.pt());
              }
            }
          }
        }
        if (mcParticlePdg == PDG_t::kLambda0Bar && !v0mcParticle.isPhysicalPrimary()) {
          auto v0mothers = v0mcParticle.mothers_as<aod::McParticles>();
          if (v0mothers.size() == 1) {
            for (const auto& v0mcParticleMother : v0mothers) {
              if (std::abs(v0mcParticleMother.eta()) > etaSel) {
                continue;
              }
              if (v0mcParticleMother.pdgCode() == PDG_t::kXiPlusBar) // Xi Plus Mother Matched
              {
                histos.fill(HIST("hAntiLambdaXiPlusFeeddownMatrix"), v0mcParticle.pt(), v0mcParticleMother.pt());
                histos.fill(HIST("hAntiLambdaFromXiPlusEtaVsPtVsPhi"), v0mcParticle.pt(), v0mcParticle.eta(), v0mcParticle.phi());
              }
              if (v0mcParticleMother.pdgCode() == -o2::constants::physics::Pdg::kXi0) // Anti Xi Zero Mother Matched
              {
                histos.fill(HIST("hAntiLambdaXiZeroFeeddownMatrix"), v0mcParticle.pt(), v0mcParticleMother.pt());
                histos.fill(HIST("hAntiLambdaFromXiZeroEtaVsPtVsPhi"), v0mcParticle.pt(), v0mcParticle.eta(), v0mcParticle.phi());
              }
              if (v0mcParticleMother.pdgCode() == PDG_t::kOmegaPlusBar) // Omega Mother Matched
              {
                histos.fill(HIST("hAntiLambdaOmegaFeeddownMatrix"), v0mcParticle.pt(), v0mcParticleMother.pt());
              }
            }
          }
        }
      }
    }
  }
  void processMixedEventHV0sInBuffer(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::PVMults>::iterator const& collision,
                                     aod::AssocV0s const& associatedV0s, aod::TriggerTracks const& triggerTracks,
                                     V0DatasWithoutTrackX const&, TracksComplete const&, aod::BCsWithTimestamps const&)
  {

    double cent = masterConfigurations.doPPAnalysis ? collision.centFT0M() : collision.centFT0C();

    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    auto bField = getMagneticField(bc.timestamp());
    // ________________________________________________
    if (efficiencyFlags.applyEfficiencyCorrection) {
      initEfficiencyFromCCDB(bc);
    }
    // ________________________________________________
    // skip if desired trigger not found
    if (triggerPresenceMap.size() > 0 && !TESTBIT(triggerPresenceMap[collision.globalIndex()], triggerBinToSelect)) {
      return;
    }

    // Perform basic event selection on both collisions
    if (((masterConfigurations.doPPAnalysis && !isCollisionSelected(collision))) || (!masterConfigurations.doPPAnalysis && !isCollisionSelectedPbPb(collision, false))) {
      return;
    }
    if (cent > axisRanges[5][1] || cent < axisRanges[5][0]) {
      return;
    }

    // ________________________________________________
    if (masterConfigurations.doFullCorrelationStudy) {
      fillCorrelationsV0(triggerTracks, associatedV0s, true, true, collision.posX(), collision.posY(), collision.posZ(), cent, bField);
    }
  }
  void processMixedEventHCascadesInBuffer(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::PVMults>::iterator const& collision,
                                          aod::AssocV0s const&, aod::AssocCascades const& associatedCascades, aod::TriggerTracks const& triggerTracks,
                                          V0DatasWithoutTrackX const&, aod::CascDatas const&, TracksComplete const&, aod::BCsWithTimestamps const&)
  {
    double cent = masterConfigurations.doPPAnalysis ? collision.centFT0M() : collision.centFT0C();
    // ________________________________________________
    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    auto bField = getMagneticField(bc.timestamp());
    if (efficiencyFlags.applyEfficiencyCorrection) {
      initEfficiencyFromCCDB(bc);
    }
    // ________________________________________________
    // skip if desired trigger not found
    if (triggerPresenceMap.size() > 0 && !TESTBIT(triggerPresenceMap[collision.globalIndex()], triggerBinToSelect)) {
      return;
    }
    // Perform basic event selection on both collisions
    if ((masterConfigurations.doPPAnalysis && !isCollisionSelected(collision)) || (!masterConfigurations.doPPAnalysis && (!isCollisionSelectedPbPb(collision, false)))) {
      return;
    }
    if (cent > axisRanges[5][1] || cent < axisRanges[5][0]) {
      return;
    }
    // ________________________________________________
    // Do hadron - cascade correlations
    if (masterConfigurations.doFullCorrelationStudy) {
      fillCorrelationsCascade(triggerTracks, associatedCascades, true, true, collision.posX(), collision.posY(), collision.posZ(), cent, bField);
    }
  }
  void processPrediction(soa::Join<aod::McCollisions, aod::McCentFT0Ms, aod::McCentFT0Cs, aod::McCentFT0As>::iterator const& mcCollision, aod::McParticles const& mcParticles)
  {
    std::vector<uint32_t> triggerIndices;
    std::vector<std::vector<uint32_t>> associatedIndices;
    std::vector<uint32_t> assocHadronIndices;
    std::vector<uint32_t> piIndices;
    std::vector<uint32_t> k0ShortIndices;
    std::vector<uint32_t> lambdaIndices;
    std::vector<uint32_t> antiLambdaIndices;
    std::vector<uint32_t> xiMinusIndices;
    std::vector<uint32_t> xiPlusIndices;
    std::vector<uint32_t> omegaMinusIndices;
    std::vector<uint32_t> omegaPlusIndices;
    float centMultFT0M = -1;
    float centMultFT0A = -1;
    float centMultFT0C = -1;
    float multFT0M = -1;
    float multFT0A = -1;
    float multFT0C = -1;
    float multEta08 = -1;
    float multEta05 = -1;
    histos.fill(HIST("Prediction/hEventSelection"), 0.5);
    // INEL>1 implies INEL>0, so only the tighter enabled selection has to be evaluated
    if (masterConfigurations.selectINELgtONE) {
      if (!o2::pwglf::isINELgt1mc(mcParticles, pdgDB)) {
        return;
      }
    } else if (masterConfigurations.selectINELgtZERO) {
      if (!o2::pwglf::isINELgt0mc(mcParticles, pdgDB)) {
        return;
      }
    }
    histos.fill(HIST("Prediction/hEventSelection"), 1.5);
    if (std::abs(mcCollision.posZ()) > masterConfigurations.zVertexCut) {
      return;
    }
    histos.fill(HIST("Prediction/hEventSelection"), 2.5);
    if (masterConfigurations.useCentralityinPrediction) {
      centMultFT0M = mcCollision.centFT0M();
      centMultFT0A = mcCollision.centFT0A();
      centMultFT0C = mcCollision.centFT0C();
    } else {
      multFT0M = mCounter.countFT0A(mcParticles) + mCounter.countFT0C(mcParticles);
      multFT0A = mCounter.countFT0A(mcParticles);
      multFT0C = mCounter.countFT0C(mcParticles);
      multEta08 = mCounter.countEta08(mcParticles);
      multEta05 = mCounter.countEta05(mcParticles);
    }
    if (!masterConfigurations.useCentralityinPrediction) {
      if (masterConfigurations.doSeparateFT0Prediction) {
        histos.fill(HIST("Prediction/hFT0AvsNchEta08"), multFT0A, multEta08);
        histos.fill(HIST("Prediction/hFT0AvsNchEta05"), multFT0A, multEta05);
        histos.fill(HIST("Prediction/hFT0CvsNchEta08"), multFT0C, multEta08);
        histos.fill(HIST("Prediction/hFT0CvsNchEta05"), multFT0C, multEta05);
      }
      histos.fill(HIST("Prediction/hFT0MvsNchEta08"), multFT0M, multEta08);
      histos.fill(HIST("Prediction/hFT0MvsNchEta05"), multFT0M, multEta05);
    }
    int iteratorNum = -1;
    for (auto const& mcParticle : mcParticles) {
      iteratorNum = iteratorNum + 1;
      double geta = mcParticle.eta();
      double gpt = mcParticle.pt();
      double gphi = mcParticle.phi();
      if (std::abs(geta) > etaSel) {
        continue;
      }
      if (std::abs(mcParticle.pdgCode()) == PDG_t::kPiPlus || std::abs(mcParticle.pdgCode()) == PDG_t::kKPlus || std::abs(mcParticle.pdgCode()) == PDG_t::kProton || std::abs(mcParticle.pdgCode()) == PDG_t::kElectron || std::abs(mcParticle.pdgCode()) == PDG_t::kMuonMinus) {
        if (!masterConfigurations.doTriggPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
          triggerIndices.emplace_back(iteratorNum);
          if (masterConfigurations.useCentralityinPrediction) {
            if (masterConfigurations.doSeparateFT0Prediction) {
              histos.fill(HIST("Prediction/hTriggerFT0A"), gpt, geta, centMultFT0A);
              histos.fill(HIST("Prediction/hTriggerFT0C"), gpt, geta, centMultFT0C);
            }
            histos.fill(HIST("Prediction/hTrigger"), gpt, geta, centMultFT0M);
          } else {
            if (masterConfigurations.doSeparateFT0Prediction) {
              histos.fill(HIST("Prediction/hTriggerFT0A"), gpt, geta, multFT0A);
              histos.fill(HIST("Prediction/hTriggerFT0C"), gpt, geta, multFT0C);
            }
            histos.fill(HIST("Prediction/hTrigger"), gpt, geta, multFT0M);
          }
        }
        if (masterConfigurations.doCorrelationHadron) {
          if (!doAssocPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
            assocHadronIndices.emplace_back(iteratorNum);
            histos.fill(HIST("Prediction/hHadron"), gpt, geta, gphi);
          }
        }
      }
      if (!doAssocPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
        if (std::abs(mcParticle.pdgCode()) == PDG_t::kPiPlus && masterConfigurations.doCorrelationPion) {
          piIndices.emplace_back(iteratorNum);
          histos.fill(HIST("Prediction/hPion"), gpt, geta, gphi);
        }
        if (mcParticle.pdgCode() == PDG_t::kK0Short && masterConfigurations.doCorrelationK0Short) {
          k0ShortIndices.emplace_back(iteratorNum);
          histos.fill(HIST("Prediction/hK0Short"), gpt, geta, gphi);
        }
        if (mcParticle.pdgCode() == PDG_t::kLambda0 && masterConfigurations.doCorrelationLambda) {
          lambdaIndices.emplace_back(iteratorNum);
          histos.fill(HIST("Prediction/hLambda"), gpt, geta, gphi);
        }
        if (mcParticle.pdgCode() == PDG_t::kLambda0Bar && masterConfigurations.doCorrelationAntiLambda) {
          antiLambdaIndices.emplace_back(iteratorNum);
          histos.fill(HIST("Prediction/hAntiLambda"), gpt, geta, gphi);
        }
        if (mcParticle.pdgCode() == PDG_t::kXiMinus && masterConfigurations.doCorrelationXiMinus) {
          xiMinusIndices.emplace_back(iteratorNum);
          histos.fill(HIST("Prediction/hXiMinus"), gpt, geta, gphi);
        }
        if (mcParticle.pdgCode() == PDG_t::kXiPlusBar && masterConfigurations.doCorrelationXiPlus) {
          xiPlusIndices.emplace_back(iteratorNum);
          histos.fill(HIST("Prediction/hXiPlus"), gpt, geta, gphi);
        }
        if (mcParticle.pdgCode() == PDG_t::kOmegaMinus && masterConfigurations.doCorrelationOmegaMinus) {
          omegaMinusIndices.emplace_back(iteratorNum);
          histos.fill(HIST("Prediction/hOmegaMinus"), gpt, geta, gphi);
        }
        if (mcParticle.pdgCode() == PDG_t::kOmegaPlusBar && masterConfigurations.doCorrelationOmegaPlus) {
          omegaPlusIndices.emplace_back(iteratorNum);
          histos.fill(HIST("Prediction/hOmegaPlus"), gpt, geta, gphi);
        }
      }
    }

    associatedIndices.emplace_back(k0ShortIndices);
    associatedIndices.emplace_back(lambdaIndices);
    associatedIndices.emplace_back(antiLambdaIndices);
    associatedIndices.emplace_back(xiMinusIndices);
    associatedIndices.emplace_back(xiPlusIndices);
    associatedIndices.emplace_back(omegaMinusIndices);
    associatedIndices.emplace_back(omegaPlusIndices);
    associatedIndices.emplace_back(piIndices);
    associatedIndices.emplace_back(assocHadronIndices);
    for (std::size_t iTrigger = 0; iTrigger < triggerIndices.size(); iTrigger++) {
      auto triggerParticle = mcParticles.iteratorAt(triggerIndices[iTrigger]);
      // check range of trigger particle
      if (triggerParticle.pt() > axisRanges[3][1] || triggerParticle.pt() < axisRanges[3][0]) {
        continue;
      }
      double getatrigger = triggerParticle.eta();
      double gphitrigger = triggerParticle.phi();
      double pttrigger = triggerParticle.pt();
      auto const& mother = triggerParticle.mothers_first_as<aod::McParticles>();
      auto globalIndex = mother.globalIndex();
      static_for<0, 8>([&](auto i) { // associated loop
        constexpr int Index = i.value;
        for (std::size_t iassoc = 0; iassoc < associatedIndices[Index].size(); iassoc++) {
          auto assocParticle = mcParticles.iteratorAt(associatedIndices[Index][iassoc]);
          if (triggerIndices[iTrigger] != associatedIndices[Index][iassoc] && globalIndex != assocParticle.globalIndex()) { // avoid self
            double getaassoc = assocParticle.eta();
            double gphiassoc = assocParticle.phi();
            double ptassoc = assocParticle.pt();
            double deltaphi = computeDeltaPhi(gphitrigger, gphiassoc);
            double deltaeta = getatrigger - getaassoc;

            // skip if basic ranges not met
            if (deltaphi < axisRanges[0][0] || deltaphi > axisRanges[0][1]) {
              continue;
            }
            if (deltaeta < axisRanges[1][0] || deltaeta > axisRanges[1][1]) {
              continue;
            }
            if (ptassoc < axisRanges[2][0] || ptassoc > axisRanges[2][1]) {
              continue;
            }
            if (TESTBIT(doCorrelation, i)) {
              if (masterConfigurations.useCentralityinPrediction) {
                histos.fill(HIST("Prediction/sameEvent/") + HIST(Particlenames[Index]), computeDeltaPhi(gphitrigger, gphiassoc), deltaeta, ptassoc, pttrigger, mcCollision.posZ(), centMultFT0M);
                if (masterConfigurations.doSeparateFT0Prediction) {
                  histos.fill(HIST("Prediction/sameEventFT0A/") + HIST(Particlenames[Index]), computeDeltaPhi(gphitrigger, gphiassoc), deltaeta, ptassoc, pttrigger, mcCollision.posZ(), centMultFT0A);
                  histos.fill(HIST("Prediction/sameEventFT0C/") + HIST(Particlenames[Index]), computeDeltaPhi(gphitrigger, gphiassoc), deltaeta, ptassoc, pttrigger, mcCollision.posZ(), centMultFT0C);
                }
              } else {
                histos.fill(HIST("Prediction/sameEvent/") + HIST(Particlenames[Index]), computeDeltaPhi(gphitrigger, gphiassoc), deltaeta, ptassoc, pttrigger, mcCollision.posZ(), multFT0M);
                if (masterConfigurations.doSeparateFT0Prediction) {
                  histos.fill(HIST("Prediction/sameEventFT0A/") + HIST(Particlenames[Index]), computeDeltaPhi(gphitrigger, gphiassoc), deltaeta, ptassoc, pttrigger, mcCollision.posZ(), multFT0A);
                  histos.fill(HIST("Prediction/sameEventFT0C/") + HIST(Particlenames[Index]), computeDeltaPhi(gphitrigger, gphiassoc), deltaeta, ptassoc, pttrigger, mcCollision.posZ(), multFT0C);
                }
              }
            }
          }
        }
      });
    }
  }
  PROCESS_SWITCH(HStrangeCorrelation, processSelectEventWithTrigger, "Select events with trigger only", true);
  PROCESS_SWITCH(HStrangeCorrelation, processSameEventHV0s, "Process same events, h-V0s", true);
  PROCESS_SWITCH(HStrangeCorrelation, processSameEventHCascades, "Process same events, h-Cascades", true);
  PROCESS_SWITCH(HStrangeCorrelation, processSameEventHPions, "Process same events, h-Pion", true);
  PROCESS_SWITCH(HStrangeCorrelation, processSameEventHHadrons, "Process same events, h-h", true);

  PROCESS_SWITCH(HStrangeCorrelation, processMixedEventHV0s, "Process mixed events, h-V0s", true);
  PROCESS_SWITCH(HStrangeCorrelation, processMixedEventHCascades, "Process mixed events, h-Cascades", true);
  PROCESS_SWITCH(HStrangeCorrelation, processMixedEventHPions, "Process mixed events, h-Pion", true);
  PROCESS_SWITCH(HStrangeCorrelation, processMixedEventHHadrons, "Process mixed events, h-h", true);
  PROCESS_SWITCH(HStrangeCorrelation, processMixedEventHV0sInBuffer, "Process mixed events in buffer, h-h", true);
  PROCESS_SWITCH(HStrangeCorrelation, processMixedEventHCascadesInBuffer, "Process mixed events in buffer, h-h", true);

  PROCESS_SWITCH(HStrangeCorrelation, processMCGenerated, "Process MC generated", false);
  PROCESS_SWITCH(HStrangeCorrelation, processClosureTest, "Process Closure Test", false);
  PROCESS_SWITCH(HStrangeCorrelation, processPairLossK0MC, "Diagnose truth-matched h-K0 pair loss", false);
  PROCESS_SWITCH(HStrangeCorrelation, processFeedDown, "process Feed Down", false);
  PROCESS_SWITCH(HStrangeCorrelation, processPrediction, "process model prediction", false);
};

thread_local HStrangeCorrelation::PairLossComparisonContext HStrangeCorrelation::pairLossComparison{};

WorkflowSpec defineDataProcessing(ConfigContext const& context)
{
  return WorkflowSpec{
    adaptAnalysisTask<HStrangeCorrelation>(context)};
}
