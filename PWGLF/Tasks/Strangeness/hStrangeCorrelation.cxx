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

#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/StaticFor.h"
#include "Framework/runDataProcessing.h"
#include <Framework/Configurable.h>

#include <TPDGCode.h>

#include <string>
#include <vector>

using namespace o2;
using namespace o2::constants::math;
using namespace o2::framework;
using namespace o2::framework::expressions;

using TracksComplete = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA>;
using V0DatasWithoutTrackX = soa::Join<aod::V0Indices, aod::V0Cores>;
using V0DatasWithoutTrackXMC = soa::Join<aod::V0Indices, aod::V0Cores, aod::McV0Labels>;

struct HStrangeCorrelation {
  // for efficiency corrections if requested
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  Service<o2::framework::O2DatabasePDG> pdgDB;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // event filtering
  Configurable<std::string> zorroMask{"zorroMask", "", "zorro trigger class to select on (empty: none)"};

  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  // master analysis switches
  Configurable<bool> doPPAnalysis{"doPPAnalysis", true, "if in pp, set to true"};
  Configurable<bool> doFullCorrelationStudy{"doFullCorrelationStudy", true, "if true, do full correlation study by creating all THnSparse histograms for the correlation function"};
  Configurable<bool> doCorrelationHadron{"doCorrelationHadron", false, "do Hadron correlation"};
  Configurable<bool> doCorrelationK0Short{"doCorrelationK0Short", true, "do K0Short correlation"};
  Configurable<bool> doCorrelationLambda{"doCorrelationLambda", false, "do Lambda correlation"};
  Configurable<bool> doCorrelationAntiLambda{"doCorrelationAntiLambda", false, "do AntiLambda correlation"};
  Configurable<bool> doCorrelationXiMinus{"doCorrelationXiMinus", false, "do XiMinus correlation"};
  Configurable<bool> doCorrelationXiPlus{"doCorrelationXiPlus", false, "do XiPlus correlation"};
  Configurable<bool> doCorrelationOmegaMinus{"doCorrelationOmegaMinus", false, "do OmegaMinus correlation"};
  Configurable<bool> doCorrelationOmegaPlus{"doCorrelationOmegaPlus", false, "do OmegaPlus correlation"};
  Configurable<bool> doCorrelationPion{"doCorrelationPion", false, "do Pion correlation"};
  Configurable<bool> doGenEventSelection{"doGenEventSelection", true, "use event selections when performing closure test for the gen events"};
  Configurable<bool> selectINELgtZERO{"selectINELgtZERO", true, "select INEL>0 events"};
  Configurable<float> zVertexCut{"zVertexCut", 10, "Cut on PV position"};
  Configurable<bool> requireAllGoodITSLayers{"requireAllGoodITSLayers", false, " require that in the event all ITS are good"};
  Configurable<bool> skipUnderOverflowInTHn{"skipUnderOverflowInTHn", false, "skip under/overflow in THns"};
  Configurable<int> mixingParameter{"mixingParameter", 10, "how many events are mixed"};
  Configurable<bool> doMCassociation{"doMCassociation", false, "fill everything only for MC associated"};
  Configurable<bool> doTriggPhysicalPrimary{"doTriggPhysicalPrimary", false, "require physical primary for trigger particles"};
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

  // used for event selections in Pb-Pb
  Configurable<int> cfgCutOccupancyHigh{"cfgCutOccupancyHigh", 3000, "High cut on TPC occupancy"};
  Configurable<int> cfgCutOccupancyLow{"cfgCutOccupancyLow", 0, "Low cut on TPC occupancy"};

  // Axes - configurable for smaller sizes
  struct : ConfigurableGroup {
    ConfigurableAxis axisMult{"axisMult", {VARIABLE_WIDTH, 0.0f, 0.01f, 1.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 70.0f, 100.0f}, "Mixing bins - multiplicity"};
    ConfigurableAxis axisVtxZ{"axisVtxZ", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
    ConfigurableAxis axisPhi{"axisPhi", {72, 0, TwoPI}, "#phi"};
    ConfigurableAxis axisEta{"axisEta", {80, -0.8, +0.8}, "#eta"};
    ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {72, -PIHalf, PIHalf * 3}, "delta #varphi axis for histograms"};
    ConfigurableAxis axisDeltaEta{"axisDeltaEta", {50, -1.6, 1.6}, "delta eta axis for histograms"};
    ConfigurableAxis axisPtAssoc{"axisPtAssoc", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt associated axis for histograms"};
    ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.0, 1.0, 2.0, 3.0, 100}, "pt associated axis for histograms"};
    ConfigurableAxis axisPtQA{"axisPtQA", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for QA histograms"};
    ConfigurableAxis axisMultCount{"axisMultCount", {VARIABLE_WIDTH, 0, 200, 400, 600, 800, 1000, 1400, 1800, 2300, 2800, 3300, 4000, 5000, 6000}, "Mixing bins - multiplicity"};
    ConfigurableAxis axisMassNSigma{"axisMassNSigma", {40, -2, 2}, "Axis for mass Nsigma"};
  } axesConfigurations;

  // for topo var QA
  struct : ConfigurableGroup {
    Configurable<float> maxPeakNSigma{"maxPeakNSigma", 5, "Peak region edge definition (in sigma)"};
    Configurable<float> minBgNSigma{"minBgNSigma", 5, "Bg region edge closest to peak (in sigma)"};
    Configurable<float> maxBgNSigma{"maxBgNSigma", 10, "Bg region edge furthest to peak (in sigma)"};
    Configurable<float> nSigmaNearXiMassCenter{"nSigmaNearXiMassCenter", 1, "for Oemga analysis only, to check if candidate mass is around Xi"};
  } massWindowConfigurations; // allows for gap between peak and bg in case someone wants to

  // Implementation of on-the-spot efficiency correction
  struct : ConfigurableGroup {
    Configurable<bool> applyEfficiencyCorrection{"applyEfficiencyCorrection", false, "apply efficiency correction"};
    Configurable<bool> applyEfficiencyForTrigger{"applyEfficiencyForTrigger", false, "apply efficiency correction for the trigger particle"};
    Configurable<bool> applyEfficiencyPropagation{"applyEfficiencyPropagation", false, "propagate also the efficiency uncertainty"};
    Configurable<bool> applyPurityHadron{"applyPurityHadron", false, "apply the purity correction for associated hadrons"};
    Configurable<bool> applyPurityTrigger{"applyPurityTrigger", false, "apply the purity correction for trigger particle"};
    Configurable<bool> applyEffAsFunctionOfMult{"applyEffAsFunctionOfMult", false, "apply efficiency as a function of multiplicity as well"};
  } efficiencyFlags;
  Configurable<std::string> ccdburl{"ccdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository to use"};
  Configurable<std::string> efficiencyCCDBPath{"efficiencyCCDBPath", "GLO/Config/GeometryAligned", "Path of the efficiency corrections"};

  // Configurables for doing subwagon systematics
  // Group all settings necessary for systematics in a specific ConfigurableGroup
  struct : ConfigurableGroup {
    std::string prefix = "systematics";
    // --- Track quality variations (single track, both trigger and assoc daughters)
    Configurable<int> minTPCNCrossedRowsTrigger{"minTPCNCrossedRowsTrigger", 70, "Minimum TPC crossed rows (trigger)"};
    Configurable<int> minTPCNCrossedRowsAssociated{"minTPCNCrossedRowsAssociated", 70, "Minimum TPC crossed rows (associated)"};
    Configurable<bool> triggerRequireITS{"triggerRequireITS", true, "require ITS signal in trigger tracks"};
    Configurable<bool> assocRequireITS{"assocRequireITS", true, "require ITS signal in associated primary tracks"};
    Configurable<int> triggerMaxTPCSharedClusters{"triggerMaxTPCSharedClusters", 200, "maximum number of shared TPC clusters (inclusive)"};
    Configurable<int> assocMaxTPCSharedClusters{"assocMaxTPCSharedClusters", 200, "maximum number of shared TPC clusters (inclusive) for assoc primary tracks"};
    Configurable<bool> triggerRequireL0{"triggerRequireL0", false, "require ITS L0 cluster for trigger"};
    Configurable<bool> assocRequireL0{"assocRequireL0", true, "require ITS L0 cluster for assoc primary track"};
    // Track quality in PbPb
    Configurable<float> tpcPidNsigmaCut{"tpcPidNsigmaCut", 5, "tpcPidNsigmaCut"};

    // --- Trigger: DCA variation from basic formula: |DCAxy| <  0.004f + (0.013f / pt)
    Configurable<float> dcaXYconstant{"dcaXYconstant", 0.004, "[0] in |DCAxy| < [0]+[1]/pT"};
    Configurable<float> dcaXYpTdep{"dcaXYpTdep", 0.013, "[1] in |DCAxy| < [0]+[1]/pT"};

    Configurable<float> dcaXYconstantAssoc{"dcaXYconstantAssoc", 0.004, "[0] in |DCAxy| < [0]+[1]/pT"};
    Configurable<float> dcaXYpTdepAssoc{"dcaXYpTdepAssoc", 0.013, "[1] in |DCAxy| < [0]+[1]/pT"};

    // --- Associated: topological variable variation (OK to vary all-at-once, at least for first study)
    Configurable<double> v0cospa{"v0cospa", 0.97, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0)
    Configurable<float> dcaV0dau{"dcaV0dau", 1.0, "DCA V0 Daughters"};
    Configurable<float> dcanegtopv{"dcanegtopv", 0.06, "DCA Neg To PV"};
    Configurable<float> dcapostopv{"dcapostopv", 0.06, "DCA Pos To PV"};
    Configurable<float> v0RadiusMin{"v0RadiusMin", 0.5, "v0radius"};
    Configurable<float> v0RadiusMax{"v0RadiusMax", 200, "v0radius"};
    // more V0 selections in PbPb
    Configurable<float> lifetimecutK0S{"lifetimecutK0S", 20, "lifetimecutK0S"};
    Configurable<float> lifetimecutLambda{"lifetimecutLambda", 30, "lifetimecutLambda"};
    Configurable<float> dcanegtopvK0S{"dcanegtopvK0S", 0.1, "DCA Neg To PV"};
    Configurable<float> dcapostopvK0S{"dcapostopvK0S", 0.1, "DCA Pos To PV"};
    Configurable<float> dcanegtopvLambda{"dcanegtopvLambda", 0.05, "DCA Neg To PV"};
    Configurable<float> dcapostopvLambda{"dcapostopvLambda", 0.2, "DCA Pos To PV"};
    Configurable<float> dcanegtopvAntiLambda{"dcanegtopvAntiLambda", 0.2, "DCA Neg To PV"};
    Configurable<float> dcapostopvAntiLambda{"dcapostopvAntiLambda", 0.05, "DCA Pos To PV"};
    // original equation: lArmPt*2>TMath::Abs(lArmAlpha) only for K0S
    Configurable<float> armPodCut{"armPodCut", 5.0f, "pT * (cut) > |alpha|, AP cut. Negative: no cut"};

    // cascade selections
    Configurable<double> cascCospa{"cascCospa", 0.95, "cascCospa"};
    Configurable<float> cascDcacascdau{"cascDcacascdau", 1.0, "cascDcacascdau"};
    Configurable<float> cascDcabachtopv{"cascDcabachtopv", 0.1, "cascDcabachtopv"};
    Configurable<float> cascRadius{"cascRadius", 0.5, "cascRadius"};
    Configurable<float> cascV0masswindow{"cascV0masswindow", 0.01, "cascV0masswindow"};
    Configurable<float> cascMindcav0topv{"cascMindcav0topv", 0.01, "cascMindcav0topv"};
    // more cascade selections in PbPb
    Configurable<float> bachBaryonCosPA{"bachBaryonCosPA", 0.9999, "Bachelor baryon CosPA"};
    Configurable<float> bachBaryonDCAxyToPV{"bachBaryonDCAxyToPV", 0.08, "DCA bachelor baryon to PV"};
    Configurable<float> dcaBaryonToPV{"dcaBaryonToPV", 0.05, "DCA of baryon doughter track To PV"};
    Configurable<float> dcaMesonToPV{"dcaMesonToPV", 0.1, "DCA of meson doughter track To PV"};
    Configurable<float> dcaBachToPV{"dcaBachToPV", 0.07, "DCA Bach To PV"};
    Configurable<float> cascdcaV0dau{"cascdcaV0dau", 0.5, "DCA V0 Daughters"};
    Configurable<float> dcaCacsDauPar0{"dcaCacsDauPar0", 0.8, " par for pt dep DCA cascade daughter cut, p_T < 1 GeV/c"};
    Configurable<float> dcaCacsDauPar1{"dcaCacsDauPar1", 0.5, " par for pt dep DCA cascade daughter cut, 1< p_T < 4 GeV/c"};
    Configurable<float> dcaCacsDauPar2{"dcaCacsDauPar2", 0.2, " par for pt dep DCA cascade daughter cut, p_T > 4 GeV/c"};
    Configurable<float> cascdcaV0ToPV{"cascdcaV0ToPV", 0.06, "DCA V0 To PV"};
    Configurable<double> cascv0cospa{"cascv0cospa", 0.98, "V0 CosPA"};
    Configurable<float> cascv0RadiusMin{"cascv0RadiusMin", 2.5, "v0radius"};
    Configurable<float> proplifetime{"proplifetime", 3, "ctau/<ctau>"};
    Configurable<float> lambdaMassWin{"lambdaMassWin", 0.005, "V0 Mass window limit"};
    Configurable<float> rejcomp{"rejcomp", 0.008, "Competing Cascade rejection"};
    Configurable<float> rapCut{"rapCut", 0.8, "Rapidity acceptance"};

    // dE/dx for associated daughters
    Configurable<int> dEdxCompatibility{"dEdxCompatibility", 1, "0: loose, 1: normal, 2: tight. Defined in HStrangeCorrelationFilter"};

    // on the fly correction instead of mixingParameter
    Configurable<bool> doOnTheFlyFlattening{"doOnTheFlyFlattening", 0, "enable an on-the-fly correction instead of using mixing"};

    // (N.B.: sources that can be investigated in post are not listed!)
  } systCuts;

  // objects to use for efficiency corrections
  TH2F* hEfficiencyTrigger;
  TH3F* hEfficiencyTriggerMult;
  TH2F* hEfficiencyPion;
  TH2F* hEfficiencyK0Short;
  TH2F* hEfficiencyLambda;
  TH2F* hEfficiencyAntiLambda;
  TH2F* hEfficiencyXiMinus;
  TH2F* hEfficiencyXiPlus;
  TH2F* hEfficiencyOmegaMinus;
  TH2F* hEfficiencyOmegaPlus;
  TH2F* hEfficiencyHadron;
  TH3F* hEfficiencyHadronMult;
  TH1F* hPurityHadron;
  TH2F* hPurityHadronMult;
  // objects to propagate the efficiency uncertainty
  TH2F* hEfficiencyUncertaintyTrigger;
  TH3F* hEfficiencyUncertaintyTriggerMult;
  TH2F* hEfficiencyUncertaintyPion;
  TH2F* hEfficiencyUncertaintyK0Short;
  TH2F* hEfficiencyUncertaintyLambda;
  TH2F* hEfficiencyUncertaintyAntiLambda;
  TH2F* hEfficiencyUncertaintyXiMinus;
  TH2F* hEfficiencyUncertaintyXiPlus;
  TH2F* hEfficiencyUncertaintyOmegaMinus;
  TH2F* hEfficiencyUncertaintyOmegaPlus;
  TH2F* hEfficiencyUncertaintyHadron;
  TH3F* hEfficiencyUncertaintyHadronMult;
  TH1F* hPurityUncertaintyHadron;
  TH2F* hPurityUncertaintyHadronMult;

  using BinningTypePP = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;
  using BinningTypePbPb = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;

  // collision slicing for mixed events
  Preslice<aod::TriggerTracks> collisionSliceTracks = aod::triggerTracks::collisionId;
  Preslice<aod::AssocV0s> collisionSliceV0s = aod::assocV0s::collisionId;
  Preslice<aod::AssocCascades> collisionSliceCascades = aod::assocCascades::collisionId;
  // Preslice<aod::AssocPions> collisionSlicePions = aod::assocHadrons::collisionId;
  Preslice<aod::AssocHadrons> collisionSliceHadrons = aod::assocHadrons::collisionId;
  Preslice<aod::McParticles> perCollision = aod::mcparticle::mcCollisionId;

  static constexpr std::string_view kV0names[] = {"K0Short", "Lambda", "AntiLambda"};
  static constexpr std::string_view kCascadenames[] = {"XiMinus", "XiPlus", "OmegaMinus", "OmegaPlus"};
  static constexpr std::string_view kParticlenames[] = {"K0Short", "Lambda", "AntiLambda", "XiMinus", "XiPlus", "OmegaMinus", "OmegaPlus", "Pion", "Hadron"};
  static constexpr int kPdgCodes[] = {310, 3122, -3122, 3312, -3312, 3334, -3334, 211};

  uint16_t doCorrelation;
  int mRunNumber;
  int mRunNumberZorro;

  std::vector<std::vector<float>> axisRanges;

  const float ctauxiPDG = 4.91;     // from PDG
  const float ctauomegaPDG = 2.461; // from PDG

  /// Function to aid in calculating delta-phi
  /// \param phi1 first phi value
  /// \param phi2 second phi value
  double computeDeltaPhi(double phi1, double phi2)
  {
    double deltaPhi = phi1 - phi2;
    double shiftedDeltaPhi = RecoDecay::constrainAngle(deltaPhi, -PIHalf);
    return shiftedDeltaPhi;
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

    TList* listEfficiencies = ccdb->getForTimeStamp<TList>(efficiencyCCDBPath, timeStamp);

    if (!listEfficiencies) {
      LOG(fatal) << "Problem getting TList object with efficiencies!";
    }

    hEfficiencyTrigger = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyTrigger"));
    hEfficiencyTriggerMult = static_cast<TH3F*>(listEfficiencies->FindObject("hEfficiencyTriggerMult"));
    hEfficiencyK0Short = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyK0Short"));
    hEfficiencyLambda = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyLambda"));
    hEfficiencyAntiLambda = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyAntiLambda"));
    hEfficiencyXiMinus = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyXiMinus"));
    hEfficiencyXiPlus = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyXiPlus"));
    hEfficiencyOmegaMinus = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyOmegaMinus"));
    hEfficiencyOmegaPlus = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyOmegaPlus"));
    hEfficiencyHadron = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyHadron"));
    hEfficiencyHadronMult = static_cast<TH3F*>(listEfficiencies->FindObject("hEfficiencyHadronMult"));
    hEfficiencyPion = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyPion"));
    hPurityHadron = static_cast<TH1F*>(listEfficiencies->FindObject("hPurityHadron"));
    hPurityHadronMult = static_cast<TH2F*>(listEfficiencies->FindObject("hPurityHadronMult"));
    hEfficiencyUncertaintyTrigger = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyUncertaintyTrigger"));
    hEfficiencyUncertaintyTriggerMult = static_cast<TH3F*>(listEfficiencies->FindObject("hEfficiencyUncertaintyTriggerMult"));
    hEfficiencyUncertaintyK0Short = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyUncertaintyK0Short"));
    hEfficiencyUncertaintyLambda = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyUncertaintyLambda"));
    hEfficiencyUncertaintyAntiLambda = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyUncertaintyAntiLambda"));
    hEfficiencyUncertaintyXiMinus = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyUncertaintyXiMinus"));
    hEfficiencyUncertaintyXiPlus = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyUncertaintyXiPlus"));
    hEfficiencyUncertaintyOmegaMinus = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyUncertaintyOmegaMinus"));
    hEfficiencyUncertaintyOmegaPlus = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyUncertaintyOmegaPlus"));
    hEfficiencyUncertaintyPion = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyUncertaintyPion"));
    hEfficiencyUncertaintyHadron = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyUncertaintyHadron"));
    hEfficiencyUncertaintyHadronMult = static_cast<TH3F*>(listEfficiencies->FindObject("hEfficiencyUncertaintyHadronMult"));
    hPurityUncertaintyHadron = static_cast<TH1F*>(listEfficiencies->FindObject("hPurityUncertaintyHadron"));
    hPurityUncertaintyHadronMult = static_cast<TH2F*>(listEfficiencies->FindObject("hPurityUncertaintyHadronMult"));
    if (efficiencyFlags.applyEfficiencyPropagation && !hEfficiencyUncertaintyTrigger)
      LOG(fatal) << "Problem getting hEfficiencyUncertaintyTrigger!";
    LOG(info) << "Efficiencies now loaded for " << mRunNumber;
  }

  template <typename TV0>
  uint64_t V0selectionBitmap(TV0 v0, float pvx, float pvy, float pvz)
  // precalculate this information so that a check is one mask operation, not many
  {
    uint64_t bitMap = 0;
    // proper lifetime , DCA daughter to prim.vtx
    if (doCorrelationK0Short) {
      // proper lifetime
      if (v0.distovertotmom(pvx, pvy, pvz) * o2::constants::physics::MassK0Short < systCuts.lifetimecutK0S)
        SETBIT(bitMap, 0);
      // DCA daughter to prim.vtx and armenteros
      if (std::abs(v0.dcapostopv()) > systCuts.dcapostopvK0S && std::abs(v0.dcanegtopv()) > systCuts.dcanegtopvK0S && v0.qtarm() * systCuts.armPodCut > std::abs(v0.alpha()))
        SETBIT(bitMap, 3);
    }
    if (doCorrelationLambda) {
      // proper lifetime
      if (v0.distovertotmom(pvx, pvy, pvz) * o2::constants::physics::MassLambda0 < systCuts.lifetimecutLambda)
        SETBIT(bitMap, 1);
      // DCA daughter to prim.vtx
      if (std::abs(v0.dcapostopv()) > systCuts.dcapostopvLambda && std::abs(v0.dcanegtopv()) > systCuts.dcanegtopvLambda)
        SETBIT(bitMap, 4);
    }
    if (doCorrelationAntiLambda) {
      // proper lifetime
      if (v0.distovertotmom(pvx, pvy, pvz) * o2::constants::physics::MassLambda0 < systCuts.lifetimecutLambda)
        SETBIT(bitMap, 2);
      // DCA daughter to prim.vtx
      if (std::abs(v0.dcapostopv()) > systCuts.dcapostopvAntiLambda && std::abs(v0.dcanegtopv()) > systCuts.dcanegtopvAntiLambda)
        SETBIT(bitMap, 5);
    }
    return bitMap;
  }

  template <typename TCascade>
  uint64_t CascadeselectionBitmap(TCascade casc, float pvx, float pvy, float pvz)
  {
    uint64_t bitMap = 0;
    float cascpos = std::hypot(casc.x() - pvx, casc.y() - pvy, casc.z() - pvz);
    float cascptotmom = std::hypot(casc.px(), casc.py(), casc.pz());
    float ctauXi = o2::constants::physics::MassXiMinus * cascpos / ((cascptotmom + 1e-13) * ctauxiPDG);
    float ctauOmega = o2::constants::physics::MassOmegaMinus * cascpos / ((cascptotmom + 1e-13) * ctauomegaPDG);
    // TPC PID and DCA daughter to prim.vtx and comopeting casc.rej and life time
    if (doCorrelationXiMinus) {
      // DCA daughter to prim.vtx
      if (std::abs(casc.dcabachtopv()) > systCuts.dcaBachToPV && std::abs(casc.dcapostopv()) > systCuts.dcaBaryonToPV &&
          std::abs(casc.dcanegtopv()) > systCuts.dcaMesonToPV)
        SETBIT(bitMap, 0);
      // comopeting casc.rej
      if (std::abs(casc.mOmega() - o2::constants::physics::MassOmegaMinus) > systCuts.rejcomp)
        SETBIT(bitMap, 4);
      if (ctauXi < systCuts.proplifetime)
        SETBIT(bitMap, 8);
      // y cut
      if (std::abs(casc.yXi()) < systCuts.rapCut)
        SETBIT(bitMap, 12);
    }
    if (doCorrelationXiPlus) {
      // DCA daughter to prim.vtx
      if (std::abs(casc.dcabachtopv()) > systCuts.dcaBachToPV && std::abs(casc.dcapostopv()) > systCuts.dcaMesonToPV &&
          std::abs(casc.dcanegtopv()) > systCuts.dcaBaryonToPV)
        SETBIT(bitMap, 1);
      // comopeting casc.rej
      if (std::abs(casc.mOmega() - o2::constants::physics::MassOmegaMinus) > systCuts.rejcomp)
        SETBIT(bitMap, 5);
      // life time
      if (ctauXi < systCuts.proplifetime)
        SETBIT(bitMap, 9);
      // y cut
      if (std::abs(casc.yXi()) > systCuts.rapCut)
        SETBIT(bitMap, 13);
    }
    if (doCorrelationOmegaMinus) {
      // DCA daughter to prim.vtx
      if (std::abs(casc.dcabachtopv()) > systCuts.dcaBachToPV && std::abs(casc.dcapostopv()) > systCuts.dcaBaryonToPV &&
          std::abs(casc.dcanegtopv()) > systCuts.dcaMesonToPV)
        SETBIT(bitMap, 2);
      // comopeting casc.rej
      if (std::abs(casc.mXi() - o2::constants::physics::MassXiMinus) > systCuts.rejcomp)
        SETBIT(bitMap, 6);
      // life time
      if (ctauOmega < systCuts.proplifetime)
        SETBIT(bitMap, 10);
      // y cut
      if (std::abs(casc.yOmega()) < systCuts.rapCut)
        SETBIT(bitMap, 14);
    }
    if (doCorrelationOmegaPlus) {
      // DCA daughter to prim.vtx
      if (std::abs(casc.dcabachtopv()) > systCuts.dcaBachToPV && std::abs(casc.dcapostopv()) > systCuts.dcaMesonToPV &&
          std::abs(casc.dcanegtopv()) > systCuts.dcaBaryonToPV)
        SETBIT(bitMap, 3);
      // comopeting casc.rej
      if (std::abs(casc.mXi() - o2::constants::physics::MassXiMinus) > systCuts.rejcomp)
        SETBIT(bitMap, 7);
      // life time
      if (ctauOmega < systCuts.proplifetime)
        SETBIT(bitMap, 11);
      // y cut
      if (std::abs(casc.yOmega()) > systCuts.rapCut)
        SETBIT(bitMap, 15);
    }
    return bitMap;
  }

  template <class TTrack>
  bool isValidTrigger(TTrack track)
  {
    if (track.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRowsTrigger) {
      return false; // crossed rows
    }
    if (!track.hasITS() && systCuts.triggerRequireITS) {
      return false; // skip, doesn't have ITS signal (skips lots of TPC-only!)
    }
    if (track.tpcNClsShared() > systCuts.triggerMaxTPCSharedClusters) {
      return false; // skip, has shared clusters
    }
    if (!(TESTBIT(track.itsClusterMap(), 0)) && systCuts.triggerRequireL0) {
      return false; // skip, doesn't have cluster in ITS L0
    }
    // systematic variations: trigger DCAxy
    if (std::abs(track.dcaXY()) > systCuts.dcaXYconstant + systCuts.dcaXYpTdep * std::abs(track.signed1Pt())) {
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
    return true;
  }
  template <class TTrack>
  bool isValidAssocHadron(TTrack track)
  {
    if (track.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRowsAssociated) {
      return false; // crossed rows
    }
    if (!track.hasITS() && systCuts.assocRequireITS) {
      return false; // skip, doesn't have ITS signal (skips lots of TPC-only!)
    }
    if (track.tpcNClsShared() > systCuts.assocMaxTPCSharedClusters) {
      return false; // skip, has shared clusters
    }
    if (!(TESTBIT(track.itsClusterMap(), 0)) && systCuts.assocRequireL0) {
      return false; // skip, doesn't have cluster in ITS L0
    }
    // systematic variations: trigger DCAxy
    if (std::abs(track.dcaXY()) > systCuts.dcaXYconstantAssoc + systCuts.dcaXYpTdepAssoc * std::abs(track.signed1Pt())) {
      return false;
    }
    if (track.pt() > axisRanges[2][1] || track.pt() < axisRanges[2][0]) {
      return false;
    }
    return true;
  }
  // V0selection in PbPb
  template <typename TV0>
  bool V0SelectedPbPb(TV0 v0)
  {
    // v0radius
    if (v0.v0radius() < systCuts.v0RadiusMin)
      return false;
    if (v0.v0radius() > systCuts.v0RadiusMax)
      return false;
    // v0cosPA
    if (v0.v0cosPA() < systCuts.v0cospa)
      return false;
    // dcaV0daughters
    if (v0.dcaV0daughters() > systCuts.dcaV0dau)
      return false;
    return true;
  }

  // cascadeselection in PbPb
  template <typename TCascade>
  bool CascadeSelectedPbPb(TCascade casc, float pvx, float pvy, float pvz)
  {
    // bachBaryonCosPA
    if (casc.bachBaryonCosPA() < systCuts.bachBaryonCosPA)
      return false;
    // bachBaryonDCAxyToPV
    if (std::abs(casc.bachBaryonDCAxyToPV()) > systCuts.bachBaryonDCAxyToPV)
      return false;
    // casccosPA
    if (casc.casccosPA(pvx, pvy, pvz) < systCuts.cascCospa)
      return false;
    // dcacascdaughters
    float ptDepCut = systCuts.dcaCacsDauPar0;
    if (casc.pt() > 1 && casc.pt() < 4)
      ptDepCut = systCuts.dcaCacsDauPar1;
    else if (casc.pt() > 4)
      ptDepCut = systCuts.dcaCacsDauPar2;
    if (casc.dcacascdaughters() > ptDepCut)
      return false;
    // dcaV0daughters
    if (casc.dcaV0daughters() > systCuts.dcaV0dau)
      return false;
    // dcav0topv
    if (std::abs(casc.dcav0topv(pvx, pvy, pvz)) < systCuts.cascdcaV0ToPV)
      return false;
    // cascradius
    if (casc.cascradius() < systCuts.cascRadius)
      return false;
    // v0radius
    if (casc.v0radius() < systCuts.cascv0RadiusMin)
      return false;
    // v0cosPA
    if (casc.v0cosPA(casc.x(), casc.y(), casc.z()) < systCuts.cascv0cospa)
      return false;
    // lambdaMassWin
    if (std::abs(casc.mLambda() - o2::constants::physics::MassLambda0) > systCuts.lambdaMassWin)
      return false;
    return true;
  }
  double calculateAverageDeltaPhiStar(double* trigg, double* assoc, double B)
  {
    double dPhiStar = 0;
    double dPhiStarMean = 0;

    double dPhi = assoc[0] - trigg[0];
    double phaseProton = (-0.3 * B * assoc[2]) / (2 * assoc[1]);
    double phaseTrack = (-0.3 * B * trigg[2]) / (2 * trigg[1]);

    for (double r = 0.8; r <= 2.5; r += 0.05) {
      dPhiStar = dPhi + std::asin(phaseProton * r) - std::asin(phaseTrack * r);
      dPhiStarMean += (dPhiStar / 34);
    }

    return dPhiStarMean;
  }
  void fillTriggerHistogram(std::shared_ptr<TH2> hist, double pt, double mult, float eff, float effUncert, float purity, float purityErr)
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
  void fillCorrelationHistogram(std::shared_ptr<THn> hist, double binFillThn[], float etaWeight, float efficiency, float totalEffUncert, float purity, float totalPurityUncert)
  {
    float previousContent, previousError2, currentContent, currentError2;
    int bin = hist->GetBin(binFillThn);
    previousContent = hist->GetBinContent(bin);
    previousError2 = hist->GetBinError2(bin);
    currentContent = previousContent + etaWeight * purity / (efficiency);
    currentError2 = previousError2 + std::pow(etaWeight * purity / (efficiency), 2) + std::pow(etaWeight * totalPurityUncert / (efficiency), 2) + std::pow(totalEffUncert * purity * etaWeight, 2) / std::pow(efficiency, 4);
    hist->SetBinContent(bin, currentContent);
    hist->SetBinError2(bin, currentError2);
  }
  void fillCorrelationsV0(aod::TriggerTracks const& triggers, aod::AssocV0s const& assocs, bool mixing, float pvx, float pvy, float pvz, float mult, double bField)
  {
    for (auto const& triggerTrack : triggers) {
      if (doTriggPhysicalPrimary && !triggerTrack.mcPhysicalPrimary())
        continue;
      auto trigg = triggerTrack.track_as<TracksComplete>();
      if (!isValidTrigger(trigg))
        continue;
      float efficiencyTrigg = 1.0f;
      float efficiencyTriggError = 0.0f;
      float purityTrigg = 1.0f;
      float purityTriggErr = 0.0;
      if (efficiencyFlags.applyEfficiencyForTrigger) {
        efficiencyTrigg = hEfficiencyTrigger->Interpolate(trigg.pt(), trigg.eta());
        if (efficiencyFlags.applyPurityTrigger)
          purityTrigg = hPurityHadron->Interpolate(trigg.pt());
        if (efficiencyFlags.applyEfficiencyPropagation) {
          efficiencyTriggError = hEfficiencyUncertaintyTrigger->Interpolate(trigg.pt(), trigg.eta());
          if (efficiencyFlags.applyPurityTrigger)
            purityTriggErr = hPurityHadron->Interpolate(trigg.pt());
        }
        if (efficiencyTrigg == 0) { // check for zero efficiency, do not apply if the case
          efficiencyTrigg = 1;
          efficiencyTriggError = 0;
        }
      }
      if (!mixing) {

        fillTriggerHistogram(histos.get<TH2>(HIST("sameEvent/TriggerParticlesV0")), trigg.pt(), mult, efficiencyTrigg, efficiencyTriggError, purityTrigg, purityTriggErr);
      }

      double triggSign = trigg.sign();
      double triggForDeltaPhiStar[] = {trigg.phi(), trigg.pt(), triggSign};

      for (auto const& assocCandidate : assocs) {
        auto assoc = assocCandidate.v0Core_as<V0DatasWithoutTrackX>();

        //---] syst cuts [---
        if ((doPPAnalysis && (assoc.v0radius() < systCuts.v0RadiusMin || assoc.v0radius() > systCuts.v0RadiusMax ||
                              std::abs(assoc.dcapostopv()) < systCuts.dcapostopv || std::abs(assoc.dcanegtopv()) < systCuts.dcanegtopv ||
                              assoc.v0cosPA() < systCuts.v0cospa || assoc.dcaV0daughters() > systCuts.dcaV0dau)))
          continue;

        if (!doPPAnalysis && !V0SelectedPbPb(assoc))
          continue;

        uint64_t selMap = V0selectionBitmap(assoc, pvx, pvy, pvz);

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
        if (postrack.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRowsAssociated || negtrack.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRowsAssociated)
          continue;

        float deltaphi = computeDeltaPhi(trigg.phi(), assoc.phi());
        float deltaeta = trigg.eta() - assoc.eta();
        float ptassoc = assoc.pt();
        float pttrigger = trigg.pt();

        // skip if basic ranges not met
        if (deltaphi < axisRanges[0][0] || deltaphi > axisRanges[0][1])
          continue;
        if (deltaeta < axisRanges[1][0] || deltaeta > axisRanges[1][1])
          continue;
        if (ptassoc < axisRanges[2][0] || ptassoc > axisRanges[2][1])
          continue;

        TH2F* hEfficiencyV0[3];
        hEfficiencyV0[0] = hEfficiencyK0Short;
        hEfficiencyV0[1] = hEfficiencyLambda;
        hEfficiencyV0[2] = hEfficiencyAntiLambda;

        TH2F* hEfficiencyUncertaintyV0[3];
        hEfficiencyUncertaintyV0[0] = hEfficiencyUncertaintyK0Short;
        hEfficiencyUncertaintyV0[1] = hEfficiencyUncertaintyLambda;
        hEfficiencyUncertaintyV0[2] = hEfficiencyUncertaintyAntiLambda;

        float etaWeight = 1;
        if (systCuts.doOnTheFlyFlattening) {
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
        if (assocCandidate.compatible(2, systCuts.dEdxCompatibility)) {
          phiProton = negtrack.phi();
          etaProton = negtrack.eta();
          ptProton = negtrack.pt();
          signProton = negtrack.sign();
        }
        double assocForDeltaPhiStar[] = {phiProton, ptProton, signProton};
        double assocForDeltaPhiStarPion[] = {phiPion, ptPion, -1};

        static_for<0, 2>([&](auto i) {
          constexpr int Index = i.value;
          float efficiency = 1.0f;
          float totalEffUncert = 0.0;
          float efficiencyError = 0.0f;
          if (efficiencyFlags.applyEfficiencyCorrection) {
            efficiency = hEfficiencyV0[Index]->Interpolate(ptassoc, assoc.eta());
            if (efficiencyFlags.applyEfficiencyPropagation)
              efficiencyError = hEfficiencyUncertaintyV0[Index]->Interpolate(ptassoc, assoc.eta());
          }
          if (efficiency == 0) { // check for zero efficiency, do not apply if the case
            efficiency = 1;
            efficiencyError = 0;
          }
          if (efficiencyFlags.applyEfficiencyPropagation) {
            totalEffUncert = std::sqrt(std::pow(efficiencyTrigg * efficiencyError, 2) + std::pow(efficiencyTriggError * efficiency, 2));
          }
          double binFillThn[6] = {deltaphi, deltaeta, ptassoc, pttrigger, pvz, mult};
          if (TESTBIT(doCorrelation, Index) && (!efficiencyFlags.applyEfficiencyCorrection || efficiency != 0) && (doPPAnalysis || (TESTBIT(selMap, Index) && TESTBIT(selMap, Index + 3)))) {
            if (assocCandidate.compatible(Index, systCuts.dEdxCompatibility) && (!doMCassociation || assocCandidate.mcTrue(Index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && !mixing && -massWindowConfigurations.maxBgNSigma < assocCandidate.invMassNSigma(Index) && assocCandidate.invMassNSigma(Index) < -massWindowConfigurations.minBgNSigma) {
              fillCorrelationHistogram(histos.get<THn>(HIST("sameEvent/LeftBg/") + HIST(kV0names[Index])), binFillThn, etaWeight, efficiency * efficiencyTrigg, totalEffUncert, purityTrigg, purityTriggErr);
              if (doDeltaPhiStarCheck) {
                double deltaPhiStar = calculateAverageDeltaPhiStar(triggForDeltaPhiStar, assocForDeltaPhiStar, bField);
                double deltaPhiStarPion = calculateAverageDeltaPhiStar(triggForDeltaPhiStar, assocForDeltaPhiStarPion, bField);
                if ((Index == 0 && triggSign > 0) || (Index == 1 && triggSign > 0) || (Index == 2 && triggSign < 0)) {
                  histos.fill(HIST("sameEvent/LeftBg/") + HIST(kV0names[Index]) + HIST("DeltaPhiStar"), deltaPhiStar, trigg.eta() - etaProton, 0.5);
                  if (Index == 0) {
                    histos.fill(HIST("sameEvent/LeftBg/") + HIST(kV0names[Index]) + HIST("DeltaPhiStar"), deltaPhiStarPion, trigg.eta() - etaPion, -0.5);
                  }
                } else {
                  histos.fill(HIST("sameEvent/LeftBg/") + HIST(kV0names[Index]) + HIST("DeltaPhiStar"), deltaPhiStar, trigg.eta() - etaProton, -0.5);
                  if (Index == 0) {
                    histos.fill(HIST("sameEvent/LeftBg/") + HIST(kV0names[Index]) + HIST("DeltaPhiStar"), deltaPhiStarPion, trigg.eta() - etaPion, 0.5);
                  }
                }
              }
            }
            if (assocCandidate.compatible(Index, systCuts.dEdxCompatibility) && (!doMCassociation || assocCandidate.mcTrue(Index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && !mixing && -massWindowConfigurations.maxPeakNSigma < assocCandidate.invMassNSigma(Index) && assocCandidate.invMassNSigma(Index) < +massWindowConfigurations.maxPeakNSigma) {
              fillCorrelationHistogram(histos.get<THn>(HIST("sameEvent/Signal/") + HIST(kV0names[Index])), binFillThn, etaWeight, efficiency * efficiencyTrigg, totalEffUncert, purityTrigg, purityTriggErr);
              if (std::abs(deltaphi) < 0.8 && doITSClustersQA) {
                histos.fill(HIST("hITSClusters") + HIST(kV0names[Index]) + HIST("NegativeDaughterToward"), ptassoc, negtrack.itsNCls(), assoc.v0radius());
                histos.fill(HIST("hITSClusters") + HIST(kV0names[Index]) + HIST("PositiveDaughterToward"), ptassoc, postrack.itsNCls(), assoc.v0radius());
              }
              if (std::abs(deltaphi) > 1 && std::abs(deltaphi) < 2 && doITSClustersQA) {
                histos.fill(HIST("hITSClusters") + HIST(kV0names[Index]) + HIST("NegativeDaughterTransverse"), ptassoc, negtrack.itsNCls(), assoc.v0radius());
                histos.fill(HIST("hITSClusters") + HIST(kV0names[Index]) + HIST("PositiveDaughterTransverse"), ptassoc, postrack.itsNCls(), assoc.v0radius());
              }
              if (doDeltaPhiStarCheck) {
                double deltaPhiStar = calculateAverageDeltaPhiStar(triggForDeltaPhiStar, assocForDeltaPhiStar, bField);
                double deltaPhiStarPion = calculateAverageDeltaPhiStar(triggForDeltaPhiStar, assocForDeltaPhiStarPion, bField);
                if ((Index == 0 && triggSign > 0) || (Index == 1 && triggSign > 0) || (Index == 2 && triggSign < 0)) {
                  histos.fill(HIST("sameEvent/Signal/") + HIST(kV0names[Index]) + HIST("DeltaPhiStar"), deltaPhiStar, trigg.eta() - etaProton, 0.5);
                  if (Index == 0) {
                    histos.fill(HIST("sameEvent/Signal/") + HIST(kV0names[Index]) + HIST("DeltaPhiStar"), deltaPhiStarPion, trigg.eta() - etaPion, -0.5);
                  }
                } else {
                  histos.fill(HIST("sameEvent/Signal/") + HIST(kV0names[Index]) + HIST("DeltaPhiStar"), deltaPhiStar, trigg.eta() - etaProton, -0.5);
                  if (Index == 0) {
                    histos.fill(HIST("sameEvent/Signal/") + HIST(kV0names[Index]) + HIST("DeltaPhiStar"), deltaPhiStarPion, trigg.eta() - etaPion, 0.5);
                  }
                }
              }
            }
            if (assocCandidate.compatible(Index, systCuts.dEdxCompatibility) && (!doMCassociation || assocCandidate.mcTrue(Index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && !mixing && +massWindowConfigurations.minBgNSigma < assocCandidate.invMassNSigma(Index) && assocCandidate.invMassNSigma(Index) < +massWindowConfigurations.maxBgNSigma) {
              fillCorrelationHistogram(histos.get<THn>(HIST("sameEvent/RightBg/") + HIST(kV0names[Index])), binFillThn, etaWeight, efficiency * efficiencyTrigg, totalEffUncert, purityTrigg, purityTriggErr);
              if (doDeltaPhiStarCheck) {
                double deltaPhiStar = calculateAverageDeltaPhiStar(triggForDeltaPhiStar, assocForDeltaPhiStar, bField);
                double deltaPhiStarPion = calculateAverageDeltaPhiStar(triggForDeltaPhiStar, assocForDeltaPhiStarPion, bField);
                if ((Index == 0 && triggSign > 0) || (Index == 1 && triggSign > 0) || (Index == 2 && triggSign < 0)) {
                  histos.fill(HIST("sameEvent/RightBg/") + HIST(kV0names[Index]) + HIST("DeltaPhiStar"), deltaPhiStar, trigg.eta() - etaProton, 0.5);
                  if (Index == 0) {
                    histos.fill(HIST("sameEvent/RightBg/") + HIST(kV0names[Index]) + HIST("DeltaPhiStar"), deltaPhiStarPion, trigg.eta() - etaPion, -0.5);
                  }
                } else {
                  histos.fill(HIST("sameEvent/RightBg/") + HIST(kV0names[Index]) + HIST("DeltaPhiStar"), deltaPhiStar, trigg.eta() - etaProton, -0.5);
                  if (Index == 0) {
                    histos.fill(HIST("sameEvent/RightBg/") + HIST(kV0names[Index]) + HIST("DeltaPhiStar"), deltaPhiStarPion, trigg.eta() - etaPion, 0.5);
                  }
                }
              }
            }
            if (assocCandidate.compatible(Index, systCuts.dEdxCompatibility) && (!doMCassociation || assocCandidate.mcTrue(Index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && mixing && -massWindowConfigurations.maxBgNSigma < assocCandidate.invMassNSigma(Index) && assocCandidate.invMassNSigma(Index) < -massWindowConfigurations.minBgNSigma)
              fillCorrelationHistogram(histos.get<THn>(HIST("mixedEvent/LeftBg/") + HIST(kV0names[Index])), binFillThn, 1, efficiency * efficiencyTrigg, totalEffUncert, purityTrigg, purityTriggErr);
            if (assocCandidate.compatible(Index, systCuts.dEdxCompatibility) && (!doMCassociation || assocCandidate.mcTrue(Index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && mixing && -massWindowConfigurations.maxPeakNSigma < assocCandidate.invMassNSigma(Index) && assocCandidate.invMassNSigma(Index) < +massWindowConfigurations.maxPeakNSigma)
              fillCorrelationHistogram(histos.get<THn>(HIST("mixedEvent/Signal/") + HIST(kV0names[Index])), binFillThn, 1, efficiency * efficiencyTrigg, totalEffUncert, purityTrigg, purityTriggErr);
            if (assocCandidate.compatible(Index, systCuts.dEdxCompatibility) && (!doMCassociation || assocCandidate.mcTrue(Index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && mixing && +massWindowConfigurations.minBgNSigma < assocCandidate.invMassNSigma(Index) && assocCandidate.invMassNSigma(Index) < +massWindowConfigurations.maxBgNSigma)
              fillCorrelationHistogram(histos.get<THn>(HIST("mixedEvent/RightBg/") + HIST(kV0names[Index])), binFillThn, 1, efficiency * efficiencyTrigg, totalEffUncert, purityTrigg, purityTriggErr);
          }
        });
      }
    }
  }

  void fillCorrelationsCascade(aod::TriggerTracks const& triggers, aod::AssocCascades const& assocs, bool mixing, float pvx, float pvy, float pvz, float mult, double bField)
  {
    for (auto const& triggerTrack : triggers) {
      if (doTriggPhysicalPrimary && !triggerTrack.mcPhysicalPrimary())
        continue;
      auto trigg = triggerTrack.track_as<TracksComplete>();
      if (!isValidTrigger(trigg))
        continue;

      float efficiencyTrigg = 1.0f;
      float efficiencyTriggError = 0.0f;
      float purityTrigg = 1.0f;
      float purityTriggErr = 0.0f;
      if (efficiencyFlags.applyEfficiencyForTrigger) {
        if (efficiencyFlags.applyEffAsFunctionOfMult) {
          efficiencyTrigg = hEfficiencyTriggerMult->Interpolate(trigg.pt(), trigg.eta(), mult);
        } else {
          efficiencyTrigg = hEfficiencyTrigger->Interpolate(trigg.pt(), trigg.eta());
        }
        if (efficiencyFlags.applyPurityTrigger) {
          if (efficiencyFlags.applyEffAsFunctionOfMult)
            purityTrigg = hPurityHadronMult->Interpolate(trigg.pt(), mult);
          else
            purityTrigg = hPurityHadron->Interpolate(trigg.pt());
        }
        if (efficiencyFlags.applyEfficiencyPropagation) {
          if (efficiencyFlags.applyEffAsFunctionOfMult)
            efficiencyTriggError = hEfficiencyUncertaintyTriggerMult->Interpolate(trigg.pt(), trigg.eta(), mult);
          else
            efficiencyTriggError = hEfficiencyUncertaintyTrigger->Interpolate(trigg.pt(), trigg.eta());
          if (efficiencyFlags.applyPurityTrigger) {
            if (efficiencyFlags.applyEffAsFunctionOfMult)
              purityTriggErr = hPurityUncertaintyHadronMult->Interpolate(trigg.pt(), mult);
            else
              purityTriggErr = hPurityUncertaintyHadron->Interpolate(trigg.pt());
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
      double triggForDeltaPhiStar[] = {trigg.phi(), trigg.pt(), triggSign};

      for (auto const& assocCandidate : assocs) {
        auto assoc = assocCandidate.cascData();

        //---] syst cuts [---
        if (doPPAnalysis && (std::abs(assoc.dcapostopv()) < systCuts.dcapostopv ||
                             std::abs(assoc.dcanegtopv()) < systCuts.dcanegtopv ||
                             std::abs(assoc.dcabachtopv()) < systCuts.cascDcabachtopv ||
                             assoc.dcaV0daughters() > systCuts.dcaV0dau ||
                             assoc.dcacascdaughters() > systCuts.cascDcacascdau ||
                             assoc.v0cosPA(pvx, pvy, pvz) < systCuts.v0cospa ||
                             assoc.casccosPA(pvx, pvy, pvz) < systCuts.cascCospa ||
                             assoc.cascradius() < systCuts.cascRadius ||
                             std::abs(assoc.dcav0topv(pvx, pvy, pvz)) < systCuts.cascMindcav0topv ||
                             std::abs(assoc.mLambda() - o2::constants::physics::MassLambda0) > systCuts.cascV0masswindow))
          continue;
        if (!doPPAnalysis && !CascadeSelectedPbPb(assoc, pvx, pvy, pvz))
          continue;
        uint64_t CascselMap = CascadeselectionBitmap(assoc, pvx, pvy, pvz);
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
        double assocForDeltaPhiStar[] = {phiProton, ptProton, signProton};
        //---] track quality check [---
        if (postrack.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRowsAssociated || negtrack.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRowsAssociated || bachtrack.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRowsAssociated)
          continue;

        float deltaphi = computeDeltaPhi(trigg.phi(), assoc.phi());
        float deltaeta = trigg.eta() - assoc.eta();
        float ptassoc = assoc.pt();
        float pttrigger = trigg.pt();

        // skip if basic ranges not met
        if (deltaphi < axisRanges[0][0] || deltaphi > axisRanges[0][1])
          continue;
        if (deltaeta < axisRanges[1][0] || deltaeta > axisRanges[1][1])
          continue;
        if (ptassoc < axisRanges[2][0] || ptassoc > axisRanges[2][1])
          continue;

        TH2F* hEfficiencyCascade[4];
        hEfficiencyCascade[0] = hEfficiencyXiMinus;
        hEfficiencyCascade[1] = hEfficiencyXiPlus;
        hEfficiencyCascade[2] = hEfficiencyOmegaMinus;
        hEfficiencyCascade[3] = hEfficiencyOmegaPlus;

        TH2F* hEfficiencyUncertaintyCascade[4];
        hEfficiencyUncertaintyCascade[0] = hEfficiencyUncertaintyXiMinus;
        hEfficiencyUncertaintyCascade[1] = hEfficiencyUncertaintyXiPlus;
        hEfficiencyUncertaintyCascade[2] = hEfficiencyUncertaintyOmegaMinus;
        hEfficiencyUncertaintyCascade[3] = hEfficiencyUncertaintyOmegaPlus;

        float etaWeight = 1;
        if (systCuts.doOnTheFlyFlattening) {
          float preWeight = 1 - std::abs(deltaeta) / 1.6;
          etaWeight = preWeight != 0 ? 1.0f / preWeight : 1.0f;
        }

        static_for<0, 3>([&](auto i) {
          constexpr int Index = i.value;
          float efficiency = 1.0f;
          float totalEffUncert = 0.0;
          float efficiencyError = 0.0f;
          if (efficiencyFlags.applyEfficiencyCorrection) {
            efficiency = hEfficiencyCascade[Index]->Interpolate(ptassoc, assoc.eta());
            if (efficiencyFlags.applyEfficiencyPropagation)
              efficiencyError = hEfficiencyUncertaintyCascade[Index]->Interpolate(ptassoc, assoc.eta());
          }
          if (efficiency == 0) { // check for zero efficiency, do not apply if the case
            efficiency = 1;
            efficiencyError = 0;
          }
          if (efficiencyFlags.applyEfficiencyPropagation) {
            totalEffUncert = std::sqrt(std::pow(efficiencyTrigg * efficiencyError, 2) + std::pow(efficiencyTriggError * efficiency, 2));
          }
          double binFillThn[6] = {deltaphi, deltaeta, ptassoc, pttrigger, pvz, mult};
          if (TESTBIT(doCorrelation, Index + 3) && (!efficiencyFlags.applyEfficiencyCorrection || efficiency != 0) && (doPPAnalysis || (TESTBIT(CascselMap, Index) && TESTBIT(CascselMap, Index + 4) && TESTBIT(CascselMap, Index + 8) && TESTBIT(CascselMap, Index + 12)))) {
            if (assocCandidate.compatible(Index, systCuts.dEdxCompatibility) && (!doMCassociation || assocCandidate.mcTrue(Index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && !mixing && -massWindowConfigurations.maxBgNSigma < assocCandidate.invMassNSigma(Index) && assocCandidate.invMassNSigma(Index) < -massWindowConfigurations.minBgNSigma) {
              fillCorrelationHistogram(histos.get<THn>(HIST("sameEvent/LeftBg/") + HIST(kCascadenames[Index])), binFillThn, etaWeight, efficiency * efficiencyTrigg, totalEffUncert, purityTrigg, purityTriggErr);
              if (doDeltaPhiStarCheck) {
                double deltaPhiStar = calculateAverageDeltaPhiStar(triggForDeltaPhiStar, assocForDeltaPhiStar, bField);
                if ((Index == 0 && triggSign > 0) || (Index == 1 && triggSign < 0) || (Index == 2 && triggSign > 0) || (Index == 3 && triggSign < 0))
                  histos.fill(HIST("sameEvent/LeftBg/") + HIST(kCascadenames[Index]) + HIST("DeltaPhiStar"), deltaPhiStar, trigg.eta() - etaProton, 0.5);
                else
                  histos.fill(HIST("sameEvent/LeftBg/") + HIST(kCascadenames[Index]) + HIST("DeltaPhiStar"), deltaPhiStar, trigg.eta() - etaProton, -0.5);
              }
            }
            if (assocCandidate.compatible(Index, systCuts.dEdxCompatibility) && (!doMCassociation || assocCandidate.mcTrue(Index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && !mixing && -massWindowConfigurations.maxPeakNSigma < assocCandidate.invMassNSigma(Index) && assocCandidate.invMassNSigma(Index) < +massWindowConfigurations.maxPeakNSigma) {
              fillCorrelationHistogram(histos.get<THn>(HIST("sameEvent/Signal/") + HIST(kCascadenames[Index])), binFillThn, etaWeight, efficiency * efficiencyTrigg, totalEffUncert, purityTrigg, purityTriggErr);
              if (doDeltaPhiStarCheck) {
                double deltaPhiStar = calculateAverageDeltaPhiStar(triggForDeltaPhiStar, assocForDeltaPhiStar, bField);
                if ((Index == 0 && triggSign > 0) || (Index == 1 && triggSign < 0) || (Index == 2 && triggSign > 0) || (Index == 3 && triggSign < 0))
                  histos.fill(HIST("sameEvent/Signal/") + HIST(kCascadenames[Index]) + HIST("DeltaPhiStar"), deltaPhiStar, trigg.eta() - etaProton, 0.5);
                else
                  histos.fill(HIST("sameEvent/Signal/") + HIST(kCascadenames[Index]) + HIST("DeltaPhiStar"), deltaPhiStar, trigg.eta() - etaProton, -0.5);
              }
            }
            if (assocCandidate.compatible(Index, systCuts.dEdxCompatibility) && (!doMCassociation || assocCandidate.mcTrue(Index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && !mixing && +massWindowConfigurations.minBgNSigma < assocCandidate.invMassNSigma(Index) && assocCandidate.invMassNSigma(Index) < +massWindowConfigurations.maxBgNSigma) {
              fillCorrelationHistogram(histos.get<THn>(HIST("sameEvent/RightBg/") + HIST(kCascadenames[Index])), binFillThn, etaWeight, efficiency * efficiencyTrigg, totalEffUncert, purityTrigg, purityTriggErr);
              if (doDeltaPhiStarCheck) {
                double deltaPhiStar = calculateAverageDeltaPhiStar(triggForDeltaPhiStar, assocForDeltaPhiStar, bField);
                if ((Index == 0 && triggSign > 0) || (Index == 1 && triggSign < 0) || (Index == 2 && triggSign > 0) || (Index == 3 && triggSign < 0))
                  histos.fill(HIST("sameEvent/RightBg/") + HIST(kCascadenames[Index]) + HIST("DeltaPhiStar"), deltaPhiStar, trigg.eta() - etaProton, 0.5);
                else
                  histos.fill(HIST("sameEvent/RightBg/") + HIST(kCascadenames[Index]) + HIST("DeltaPhiStar"), deltaPhiStar, trigg.eta() - etaProton, -0.5);
              }
            }
            if (assocCandidate.compatible(Index, systCuts.dEdxCompatibility) && (!doMCassociation || assocCandidate.mcTrue(Index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && mixing && -massWindowConfigurations.maxBgNSigma < assocCandidate.invMassNSigma(Index) && assocCandidate.invMassNSigma(Index) < -massWindowConfigurations.minBgNSigma)
              fillCorrelationHistogram(histos.get<THn>(HIST("mixedEvent/LeftBg/") + HIST(kCascadenames[Index])), binFillThn, 1, efficiency * efficiencyTrigg, totalEffUncert, purityTrigg, purityTriggErr);
            if (assocCandidate.compatible(Index, systCuts.dEdxCompatibility) && (!doMCassociation || assocCandidate.mcTrue(Index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && mixing && -massWindowConfigurations.maxPeakNSigma < assocCandidate.invMassNSigma(Index) && assocCandidate.invMassNSigma(Index) < +massWindowConfigurations.maxPeakNSigma)
              fillCorrelationHistogram(histos.get<THn>(HIST("mixedEvent/Signal/") + HIST(kCascadenames[Index])), binFillThn, 1, efficiency * efficiencyTrigg, totalEffUncert, purityTrigg, purityTriggErr);
            if (assocCandidate.compatible(Index, systCuts.dEdxCompatibility) && (!doMCassociation || assocCandidate.mcTrue(Index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && mixing && +massWindowConfigurations.minBgNSigma < assocCandidate.invMassNSigma(Index) && assocCandidate.invMassNSigma(Index) < +massWindowConfigurations.maxBgNSigma)
              fillCorrelationHistogram(histos.get<THn>(HIST("mixedEvent/RightBg/") + HIST(kCascadenames[Index])), binFillThn, 1, efficiency * efficiencyTrigg, totalEffUncert, purityTrigg, purityTriggErr);
          }
        });
      }
    }
  }
  template <typename TTriggers, typename THadrons>
  void fillCorrelationsHadron(TTriggers const& triggers, THadrons const& assocs, bool mixing, float pvz, float mult, double bField)
  {

    for (auto const& triggerTrack : triggers) {
      if (doTriggPhysicalPrimary && !triggerTrack.mcPhysicalPrimary())
        continue;
      auto trigg = triggerTrack.template track_as<TracksComplete>();
      if (!isValidTrigger(trigg))
        continue;

      float efficiencyTrigger = 1.0f;
      float efficiencyTriggerError = 0.0f;
      float purityTrigger = 1.0f;
      float purityTriggerError = 0.0f;
      if (efficiencyFlags.applyEfficiencyForTrigger) {
        if (efficiencyFlags.applyEffAsFunctionOfMult)
          efficiencyTrigger = hEfficiencyTriggerMult->Interpolate(trigg.pt(), trigg.eta(), mult);
        else
          efficiencyTrigger = hEfficiencyTrigger->Interpolate(trigg.pt(), trigg.eta());
        if (efficiencyFlags.applyPurityTrigger) {
          if (efficiencyFlags.applyEffAsFunctionOfMult)
            purityTrigger = hPurityHadronMult->Interpolate(trigg.pt(), mult);
          else
            purityTrigger = hPurityHadron->Interpolate(trigg.pt());
        }
        if (efficiencyFlags.applyEfficiencyPropagation) {
          if (efficiencyFlags.applyEffAsFunctionOfMult)
            efficiencyTriggerError = hEfficiencyUncertaintyTriggerMult->Interpolate(trigg.pt(), trigg.eta(), mult);
          else
            efficiencyTriggerError = hEfficiencyUncertaintyTrigger->Interpolate(trigg.pt(), trigg.eta());
          if (efficiencyFlags.applyPurityTrigger) {
            if (efficiencyFlags.applyEffAsFunctionOfMult)
              purityTriggerError = hPurityUncertaintyHadronMult->Interpolate(trigg.pt(), mult);
            else
              purityTriggerError = hPurityUncertaintyHadron->Interpolate(trigg.pt());
          }
        }
        if (efficiencyTrigger == 0) { // check for zero efficiency, do not apply if the case
          efficiencyTrigger = 1;
          efficiencyTriggerError = 0;
        }
      }
      if (!mixing) {
        if constexpr (requires { triggerTrack.extra(); })
          fillTriggerHistogram(histos.get<TH2>(HIST("sameEvent/TriggerParticlesPion")), trigg.pt(), mult, efficiencyTrigger, efficiencyTriggerError, purityTrigger, purityTriggerError);
        else
          fillTriggerHistogram(histos.get<TH2>(HIST("sameEvent/TriggerParticlesHadron")), trigg.pt(), mult, efficiencyTrigger, efficiencyTriggerError, purityTrigger, purityTriggerError);
      }
      double triggSign = trigg.sign();
      double triggForDeltaPhiStar[] = {trigg.phi(), trigg.pt(), triggSign};
      for (auto const& assocTrack : assocs) {
        auto assoc = assocTrack.template track_as<TracksComplete>();

        //---] removing autocorrelations [---
        if (doAutocorrelationRejection) {
          if (trigg.globalIndex() == assoc.globalIndex()) {
            if constexpr (requires { assocTrack.nSigmaTPCPi(); })
              histos.fill(HIST("hNumberOfRejectedPairsPion"), 0.5);
            else
              histos.fill(HIST("hNumberOfRejectedPairsHadron"), 0.5);
            continue;
          }
        }
        //---] track quality check [---
        if (!isValidAssocHadron(assoc))
          continue;
        if (doAssocPhysicalPrimary && !assocTrack.mcPhysicalPrimary()) {
          continue;
        }
        float deltaphi = computeDeltaPhi(trigg.phi(), assoc.phi());
        float deltaeta = trigg.eta() - assoc.eta();
        float ptassoc = assoc.pt();
        float pttrigger = trigg.pt();

        double assocSign = assoc.sign();
        double assocForDeltaPhiStar[] = {assoc.phi(), assoc.pt(), assocSign};

        float etaWeight = 1.;
        if (systCuts.doOnTheFlyFlattening) {
          float preWeight = 1 - std::abs(deltaeta) / 1.6;
          etaWeight = preWeight != 0 ? 1.0f / preWeight : 1.0f;
        }

        // skip if basic ranges not met
        if (deltaphi < axisRanges[0][0] || deltaphi > axisRanges[0][1])
          continue;
        if (deltaeta < axisRanges[1][0] || deltaeta > axisRanges[1][1])
          continue;
        if (ptassoc < axisRanges[2][0] || ptassoc > axisRanges[2][1])
          continue;

        float efficiency = 1;
        float purity = 1.0f;
        float purityUncertainty = 0.0f;
        float totalEffUncert = 0.0;
        float efficiencyUncertainty = 0.0f;
        float totalPurityUncert = 0.0;
        if (efficiencyFlags.applyEfficiencyCorrection) {
          if constexpr (requires { assocTrack.nSigmaTPCPi(); }) {
            efficiency = hEfficiencyPion->Interpolate(ptassoc, assoc.eta());
            if (efficiencyFlags.applyEfficiencyPropagation)
              efficiencyUncertainty = hEfficiencyUncertaintyPion->Interpolate(ptassoc, assoc.eta());
          } else {
            if (efficiencyFlags.applyEffAsFunctionOfMult)
              efficiency = hEfficiencyHadronMult->Interpolate(ptassoc, assoc.eta(), mult);
            else
              efficiency = hEfficiencyHadron->Interpolate(ptassoc, assoc.eta());
            if (efficiencyFlags.applyPurityHadron) {
              if (efficiencyFlags.applyEffAsFunctionOfMult)
                purity = hPurityHadronMult->Interpolate(ptassoc, mult);
              else
                purity = hPurityHadron->Interpolate(ptassoc);
            }
            if (efficiencyFlags.applyEfficiencyPropagation) {
              if (efficiencyFlags.applyEffAsFunctionOfMult)
                efficiencyUncertainty = hEfficiencyUncertaintyHadronMult->Interpolate(ptassoc, assoc.eta(), mult);
              else
                efficiencyUncertainty = hEfficiencyUncertaintyHadron->Interpolate(ptassoc, assoc.eta());
              if (efficiencyFlags.applyPurityHadron) {
                if (efficiencyFlags.applyEffAsFunctionOfMult)
                  purityUncertainty = hPurityUncertaintyHadronMult->Interpolate(ptassoc, mult);
                else
                  purityUncertainty = hPurityUncertaintyHadron->Interpolate(ptassoc);
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
        double binFillThn[6] = {deltaphi, deltaeta, ptassoc, pttrigger, pvz, mult};
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
    mRunNumber = 0;
    mRunNumberZorro = 0;
    hEfficiencyPion = 0x0;
    hEfficiencyK0Short = 0x0;
    hEfficiencyLambda = 0x0;
    hEfficiencyAntiLambda = 0x0;
    hEfficiencyXiMinus = 0x0;
    hEfficiencyXiPlus = 0x0;
    hEfficiencyOmegaMinus = 0x0;
    hEfficiencyOmegaPlus = 0x0;
    hEfficiencyUncertaintyTrigger = 0x0;
    hEfficiencyUncertaintyXiMinus = 0x0;
    hEfficiencyUncertaintyXiPlus = 0x0;
    hEfficiencyUncertaintyOmegaMinus = 0x0;
    hEfficiencyUncertaintyOmegaPlus = 0x0;
    hEfficiencyUncertaintyPion = 0x0;
    hEfficiencyUncertaintyK0Short = 0x0;
    hEfficiencyUncertaintyLambda = 0x0;
    hEfficiencyUncertaintyAntiLambda = 0x0;

    hEfficiencyHadron = 0x0;
    hPurityHadron = 0x0;
    hPurityUncertaintyHadron = 0x0;
    hEfficiencyUncertaintyHadron = 0x0;

    // set bitmap for convenience
    doCorrelation = 0;
    if (doCorrelationK0Short)
      SETBIT(doCorrelation, 0);
    if (doCorrelationLambda)
      SETBIT(doCorrelation, 1);
    if (doCorrelationAntiLambda)
      SETBIT(doCorrelation, 2);
    if (doCorrelationXiMinus)
      SETBIT(doCorrelation, 3);
    if (doCorrelationXiPlus)
      SETBIT(doCorrelation, 4);
    if (doCorrelationOmegaMinus)
      SETBIT(doCorrelation, 5);
    if (doCorrelationOmegaPlus)
      SETBIT(doCorrelation, 6);
    if (doCorrelationPion)
      SETBIT(doCorrelation, 7);
    if (doCorrelationHadron)
      SETBIT(doCorrelation, 8);

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

    // store the original axes in specific TH1Cs for completeness
    histos.add("axes/hDeltaPhiAxis", "", kTH1C, {preAxisDeltaPhi});
    histos.add("axes/hDeltaEtaAxis", "", kTH1C, {preAxisDeltaEta});
    histos.add("axes/hPtAssocAxis", "", kTH1C, {preAxisPtAssoc});
    histos.add("axes/hPtTriggerAxis", "", kTH1C, {preAxisPtTrigger});
    histos.add("axes/hVertexZAxis", "", kTH1C, {preAxisVtxZ});
    histos.add("axes/hMultAxis", "", kTH1C, {preAxisMult});

    std::vector<double> edgesDeltaPhiOrig = preAxisDeltaPhi.binEdges;
    std::vector<double> edgesDeltaEtaOrig = preAxisDeltaEta.binEdges;
    std::vector<double> edgesPtAssocOrig = preAxisPtAssoc.binEdges;
    std::vector<double> edgesPtTriggerOrig = preAxisPtTrigger.binEdges;
    std::vector<double> edgesVtxZOrig = preAxisVtxZ.binEdges;
    std::vector<double> edgesMultOrig = preAxisMult.binEdges;

    std::vector<float> rangesDeltaPhi = {static_cast<float>(edgesDeltaPhiOrig[0]), static_cast<float>(edgesDeltaPhiOrig[edgesDeltaPhiOrig.size() - 1])};
    std::vector<float> rangesDeltaEta = {static_cast<float>(edgesDeltaEtaOrig[0]), static_cast<float>(edgesDeltaEtaOrig[edgesDeltaEtaOrig.size() - 1])};
    std::vector<float> rangesPtAssoc = {static_cast<float>(edgesPtAssocOrig[0]), static_cast<float>(edgesPtAssocOrig[edgesPtAssocOrig.size() - 1])};
    std::vector<float> rangesPtTrigger = {static_cast<float>(edgesPtTriggerOrig[0]), static_cast<float>(edgesPtTriggerOrig[edgesPtTriggerOrig.size() - 1])};
    std::vector<float> rangesVtxZ = {static_cast<float>(edgesVtxZOrig[0]), static_cast<float>(edgesVtxZOrig[edgesVtxZOrig.size() - 1])};
    std::vector<float> rangesMult = {static_cast<float>(edgesMultOrig[0]), static_cast<float>(edgesMultOrig[edgesMultOrig.size() - 1])};

    axisRanges.emplace_back(rangesDeltaPhi);
    axisRanges.emplace_back(rangesDeltaEta);
    axisRanges.emplace_back(rangesPtAssoc);
    axisRanges.emplace_back(rangesPtTrigger);
    axisRanges.emplace_back(rangesVtxZ);
    axisRanges.emplace_back(rangesMult);

    std::vector<double> edgesDeltaPhi;
    std::vector<double> edgesDeltaEta;
    std::vector<double> edgesPtAssoc;
    std::vector<double> edgesPtTrigger;
    std::vector<double> edgesVtxZ;
    std::vector<double> edgesMult;

    // v--- skipUnderOverflowInTHn ---v
    //
    // if enabled, this will change the axes such that they will solely cover the interval from
    // edge[1] to edge[n-1]; this will mean that the bin 1 and bin N will be stored in
    // under / overflow bins and will have to be manually unpacked. Do not forget to do the manual
    // unpacking a posteriori!
    //
    // this feature is meant to save memory conveniently.
    // it should actually be implemented centrally in ROOT but ok, this will do it for now.

    int offset = skipUnderOverflowInTHn ? 1 : 0;
    // ===] delta-phi [===
    if (!preAxisDeltaPhi.nBins.has_value()) {
      // variable binning, use bins provided
      for (int i = offset; i < static_cast<int>(edgesDeltaPhiOrig.size()) - offset; i++)
        edgesDeltaPhi.emplace_back(edgesDeltaPhiOrig[i]);
    } else {
      // fixed binning, generate the bin edges on-the-spot
      double min = edgesDeltaPhiOrig[0];
      double delta = (edgesDeltaPhiOrig[1] - edgesDeltaPhiOrig[0]) / preAxisDeltaPhi.nBins.value();
      for (int i = offset; i < preAxisDeltaPhi.nBins.value() + 1 - offset; i++)
        edgesDeltaPhi.emplace_back(min + static_cast<double>(i) * delta);
    }
    // ===] delta-eta [===
    if (!preAxisDeltaEta.nBins.has_value()) {
      // variable binning, use bins provided
      for (int i = offset; i < static_cast<int>(edgesDeltaEtaOrig.size()) - offset; i++)
        edgesDeltaEta.emplace_back(edgesDeltaEtaOrig[i]);
    } else {
      // fixed binning, generate the bin edges on-the-spot
      double min = edgesDeltaEtaOrig[0];
      double delta = (edgesDeltaEtaOrig[1] - edgesDeltaEtaOrig[0]) / preAxisDeltaEta.nBins.value();
      for (int i = offset; i < preAxisDeltaEta.nBins.value() + 1 - offset; i++)
        edgesDeltaEta.emplace_back(min + static_cast<double>(i) * delta);
    }
    // ===] pt assoc [===
    if (!preAxisPtAssoc.nBins.has_value()) {
      // variable binning, use bins provided
      for (int i = offset; i < static_cast<int>(edgesPtAssocOrig.size()) - offset; i++)
        edgesPtAssoc.emplace_back(edgesPtAssocOrig[i]);
    } else {
      // fixed binning, generate the bin edges on-the-spot
      double min = edgesPtAssocOrig[0];
      double delta = (edgesPtAssocOrig[1] - edgesPtAssocOrig[0]) / preAxisPtAssoc.nBins.value();
      for (int i = offset; i < preAxisPtAssoc.nBins.value() + 1 - offset; i++)
        edgesPtAssoc.emplace_back(min + static_cast<double>(i) * delta);
    }
    // ===] pt trigger [===
    if (!preAxisPtTrigger.nBins.has_value()) {
      // variable binning, use bins provided
      for (int i = offset; i < static_cast<int>(edgesPtTriggerOrig.size()) - offset; i++)
        edgesPtTrigger.emplace_back(edgesPtTriggerOrig[i]);
    } else {
      // fixed binning, generate the bin edges on-the-spot
      double min = edgesPtTriggerOrig[0];
      double delta = (edgesPtTriggerOrig[1] - edgesPtTriggerOrig[0]) / preAxisPtTrigger.nBins.value();
      for (int i = offset; i < preAxisPtTrigger.nBins.value() + 1 - offset; i++)
        edgesPtTrigger.emplace_back(min + static_cast<double>(i) * delta);
    }
    // ===] vtx Z [===
    if (!preAxisVtxZ.nBins.has_value()) {
      // variable binning, use bins provided
      for (int i = offset; i < static_cast<int>(edgesVtxZOrig.size()) - offset; i++)
        edgesVtxZ.emplace_back(edgesVtxZOrig[i]);
    } else {
      // fixed binning, generate the bin edges on-the-spot
      double min = edgesVtxZOrig[0];
      double delta = (edgesVtxZOrig[1] - edgesVtxZOrig[0]) / preAxisVtxZ.nBins.value();
      for (int i = offset; i < preAxisVtxZ.nBins.value() + 1 - offset; i++)
        edgesVtxZ.emplace_back(min + static_cast<double>(i) * delta);
    }
    // ===] mult percentile [===
    if (!preAxisMult.nBins.has_value()) {
      // variable binning, use bins provided
      for (int i = offset; i < static_cast<int>(edgesMultOrig.size()) - offset; i++)
        edgesMult.emplace_back(edgesMultOrig[i]);
    } else {
      // fixed binning, generate the bin edges on-the-spot
      double min = edgesMultOrig[0];
      double delta = (edgesMultOrig[1] - edgesMultOrig[0]) / preAxisMult.nBins.value();
      for (int i = offset; i < preAxisMult.nBins.value() + 1 - offset; i++)
        edgesMult.emplace_back(min + static_cast<double>(i) * delta);
    }

    LOGF(info, "Initialized THnF axis delta-phi with %i bins.", edgesDeltaPhi.size() - 1);
    LOGF(info, "Initialized THnF axis delta-eta with %i bins.", edgesDeltaEta.size() - 1);
    LOGF(info, "Initialized THnF axis pTassoc with %i bins.", edgesPtAssoc.size() - 1);
    LOGF(info, "Initialized THnF axis pTtrigger with %i bins.", edgesPtTrigger.size() - 1);
    LOGF(info, "Initialized THnF axis vertex-Z with %i bins.", edgesVtxZ.size() - 1);
    LOGF(info, "Initialized THnF axis multiplicity with %i bins.", edgesMult.size() - 1);

    const AxisSpec axisDeltaPhiNDim{edgesDeltaPhi, "#Delta#varphi"};
    const AxisSpec axisDeltaEtaNDim{edgesDeltaEta, "#Delta#eta"};
    const AxisSpec axisPtAssocNDim{edgesPtAssoc, "#it{p}_{T}^{assoc} (GeV/c)"};
    const AxisSpec axisPtTriggerNDim{edgesPtTrigger, "#it{p}_{T}^{trigger} (GeV/c)"};
    const AxisSpec axisVtxZNDim{edgesVtxZ, "vertex Z (cm)"};
    const AxisSpec axisMultNDim{edgesMult, "mult percentile"};
    if (!doPPAnalysis) {
      // event selections in Pb-Pb
      histos.add("hEventSelection", "hEventSelection", kTH1F, {{10, 0, 10}});
      TString eventSelLabel[] = {"all", "sel8", "kIsTriggerTVX", "PV_{z}", "kIsGoodITSLayersAll", "kIsGoodZvtxFT0vsPV", "OccupCut", "kNoITSROFrameBorder", "kNoSameBunchPileup ", " kNoCollInTimeRangeStandard"};
      for (int i = 1; i <= histos.get<TH1>(HIST("hEventSelection"))->GetNbinsX(); i++) {
        histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(i, eventSelLabel[i - 1]);
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
    if (doprocessSameEventHV0s || doprocessSameEventHCascades || doprocessSameEventHPions || doprocessSameEventHHadrons) {
      histos.add("hTriggerAllSelectedEtaVsPt", "hTriggerAllSelectedEtaVsPt", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
      histos.add("hPositiveTriggerPrimaryEtaVsPt", "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
      histos.add("hNegativeTriggerPrimaryEtaVsPt", "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
      // QA and THn Histograms
      histos.add("hTriggerPtResolution", ";p_{T}^{reconstructed} (GeV/c); p_{T}^{generated} (GeV/c)", kTH2F, {axesConfigurations.axisPtQA, axesConfigurations.axisPtQA});
      histos.add("hTriggerPrimaryEtaVsPt", "hTriggerPrimaryEtaVsPt", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
      histos.add("hTrackEtaVsPtVsPhi", "hTrackEtaVsPtVsPhi", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi});
      histos.add("hAssocTrackEtaVsPtVsPhi", "hAssocTrackEtaVsPtVsPhi", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi});
      // histos.add("hTrackAttempt", "Attempt", kTH3F, {axisPtQA, axisEta, axisPhi});
    }
    if (doprocessSameEventHPions || doprocessSameEventHHadrons || doprocessMixedEventHPions || doprocessMixedEventHHadrons) {
      histos.add("hNumberOfRejectedPairsHadron", "hNumberOfRejectedPairsHadron", kTH1F, {{1, 0, 1}});
      histos.add("hNumberOfRejectedPairsPion", "hNumberOfRejectedPairsPion", kTH1F, {{1, 0, 1}});
    }
    if (doprocessSameEventHV0s || doprocessMixedEventHV0s) {
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
    for (int i = 0; i < 9; i++) {
      if (TESTBIT(doCorrelation, i)) {
        if (doFullCorrelationStudy)
          histos.add(fmt::format("sameEvent/Signal/{}", kParticlenames[i]).c_str(), "", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
        if (doDeltaPhiStarCheck && doFullCorrelationStudy) {
          histos.add(fmt::format("sameEvent/Signal/{}DeltaPhiStar", kParticlenames[i]).c_str(), "", kTH3F, {{100, -0.3, 0.3}, {50, -0.05, 0.05}, {2, -1, 1}}); // -1 oposite charge, 1 same charge
        }
        if (i < 7) {
          histos.add(fmt::format("h{}EtaVsPtVsPhi", kParticlenames[i]).c_str(), "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi});
          histos.add(fmt::format("h3d{}Spectrum", kParticlenames[i]).c_str(), fmt::format("h3d{}Spectrum", kParticlenames[i]).c_str(), kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisMult, axesConfigurations.axisMassNSigma});
          histos.add(fmt::format("h3d{}SpectrumY", kParticlenames[i]).c_str(), fmt::format("h3d{}SpectrumY", kParticlenames[i]).c_str(), kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisMult, axesConfigurations.axisMassNSigma});
          hStrange = true;
          histos.add(fmt::format("h{}EtaVsPtVsPhiBg", kParticlenames[i]).c_str(), "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi});
          if (doITSClustersQA) {
            histos.add(fmt::format("hITSClusters{}NegativeDaughterToward", kParticlenames[i]).c_str(), "", kTH3F, {axesConfigurations.axisPtAssoc, {8, -0.5, 7.5}, {20, 0, 10}});
            histos.add(fmt::format("hITSClusters{}PositiveDaughterToward", kParticlenames[i]).c_str(), "", kTH3F, {axesConfigurations.axisPtAssoc, {8, -0.5, 7.5}, {20, 0, 10}});
            histos.add(fmt::format("hITSClusters{}NegativeDaughterTransverse", kParticlenames[i]).c_str(), "", kTH3F, {axesConfigurations.axisPtAssoc, {8, -0.5, 7.5}, {20, 0, 10}});
            histos.add(fmt::format("hITSClusters{}PositiveDaughterTransverse", kParticlenames[i]).c_str(), "", kTH3F, {axesConfigurations.axisPtAssoc, {8, -0.5, 7.5}, {20, 0, 10}});
          }
        }
      }
    }
    if (TESTBIT(doCorrelation, 7)) {
      histos.add("hPionEtaVsPtAllSelected", "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
      histos.add("hPionEtaVsPt", "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
      histos.add("hPositivePionEtaVsPt", "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
      histos.add("hNegativePionEtaVsPt", "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
    }
    if (TESTBIT(doCorrelation, 8)) {
      histos.add("hAsssocTrackEtaVsPtVsPhi", "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi});
      histos.add("hAssocPrimaryEtaVsPt", "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
      histos.add("hAssocHadronsAllSelectedEtaVsPt", "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
      histos.add("hAssocPtResolution", ";p_{T}^{reconstructed} (GeV/c); p_{T}^{generated} (GeV/c)", kTH2F, {axesConfigurations.axisPtQA, axesConfigurations.axisPtQA});
    }

    if (hStrange && doFullCorrelationStudy) {
      histos.addClone("sameEvent/Signal/", "sameEvent/LeftBg/");
      histos.addClone("sameEvent/Signal/", "sameEvent/RightBg/");
    }

    LOGF(info, "Init THnFs done");
    // mixed-event correlation functions
    if ((doprocessMixedEventHV0s || doprocessMixedEventHCascades || doprocessMixedEventHPions || doprocessMixedEventHHadrons) && doFullCorrelationStudy) {
      histos.addClone("sameEvent/", "mixedEvent/");
    }
    if (doprocessSameEventHHadrons && doFullCorrelationStudy)
      histos.add("sameEvent/TriggerParticlesHadron", "TriggersHadron", kTH2F, {axesConfigurations.axisPtQA, axesConfigurations.axisMult});
    if (doprocessSameEventHV0s && doFullCorrelationStudy)
      histos.add("sameEvent/TriggerParticlesV0", "TriggersV0", kTH2F, {axesConfigurations.axisPtQA, axesConfigurations.axisMult});
    if (doprocessSameEventHCascades && doFullCorrelationStudy)
      histos.add("sameEvent/TriggerParticlesCascade", "TriggersCascade", kTH2F, {axesConfigurations.axisPtQA, axesConfigurations.axisMult});
    if (doprocessSameEventHPions && doFullCorrelationStudy)
      histos.add("sameEvent/TriggerParticlesPion", "TriggersPion", kTH2F, {axesConfigurations.axisPtQA, axesConfigurations.axisMult});

    // MC generated plots
    if (doprocessMCGenerated) {
      histos.add("Generated/hTrigger", "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
      histos.add("Generated/hPositiveTrigger", "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
      histos.add("Generated/hNegativeTrigger", "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
      for (int i = 0; i < 9; i++) {
        histos.add(fmt::format("Generated/h{}", kParticlenames[i]).c_str(), "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
        if (i == 7) {
          histos.add(fmt::format("Generated/hPositive{}", kParticlenames[i]).c_str(), "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
          histos.add(fmt::format("Generated/hNegative{}", kParticlenames[i]).c_str(), "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
        }
      }
      histos.addClone("Generated/", "GeneratedWithPV/");

      // histograms within |y|<0.5, vs multiplicity
      for (int i = 0; i < 8; i++) {
        histos.add(fmt::format("GeneratedWithPV/h{}_MidYVsMult", kParticlenames[i]).c_str(), "", kTH2F, {axesConfigurations.axisPtQA, axesConfigurations.axisMult});
        histos.add(fmt::format("GeneratedWithPV/h{}_MidYVsMult_TwoPVsOrMore", kParticlenames[i]).c_str(), "", kTH2F, {axesConfigurations.axisPtQA, axesConfigurations.axisMult});
      }
    }
    if (doprocessClosureTest) {
      for (int i = 0; i < 9; i++) {
        if (TESTBIT(doCorrelation, i))
          histos.add(fmt::format("ClosureTest/sameEvent/{}", kParticlenames[i]).c_str(), "", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisPtTriggerNDim, axisVtxZNDim, axisMultNDim});
        if (TESTBIT(doCorrelation, i))
          histos.add(fmt::format("ClosureTest/h{}", kParticlenames[i]).c_str(), "", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisPhi});
      }
      histos.add("ClosureTest/hTrigger", "Trigger Tracks", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisEta, axesConfigurations.axisMult});
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

    // visual inspection of sizes
    histos.print();

    // initialize CCDB *only* if efficiency correction requested
    // skip if not requested, saves a bit of time
    if (efficiencyFlags.applyEfficiencyCorrection) {
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
    if (std::abs(collision.posZ()) > zVertexCut) {
      return false;
    }
    if (collision.centFT0M() > axisRanges[5][1] || collision.centFT0M() < axisRanges[5][0]) {
      return false;
    }
    if (!collision.isInelGt0() && selectINELgtZERO) {
      return false;
    }
    if (!collision.selection_bit(aod::evsel::kIsGoodITSLayersAll) && requireAllGoodITSLayers) {
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
  bool isCollisionSelectedPbPb(TCollision collision, bool fillHists)
  {
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 0.5 /* all collisions */);

    // Perform basic event selection
    if (!collision.sel8()) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 1.5 /* collisions  after sel8*/);

    if (!collision.selection_bit(aod::evsel::kIsTriggerTVX)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 2.5 /* FT0 vertex (acceptable FT0C-FT0A time difference) collisions */);

    if (std::abs(collision.posZ()) > zVertexCut) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 3.5 /* collisions  after sel pvz sel*/);

    if (!collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      // cut time intervals with dead ITS staves
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 4.5 /* collisions  after cut time intervals with dead ITS staves*/);

    if (!collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      // removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference
      // use this cut at low multiplicities with caution
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 5.5 /* removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference*/);

    auto occupancy = collision.trackOccupancyInTimeRange();
    if (occupancy < cfgCutOccupancyLow || occupancy > cfgCutOccupancyHigh)
      return false;
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 6.5 /* Below min occupancy and Above max occupancy*/);

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
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 7.5 /* reject collisions close to Time Frame borders*/);

    if (!collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      // reject events affected by the ITS ROF border
      // O2-4309
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 8.5 /* reject events affected by the ITS ROF border*/);

    if (!collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      // rejects collisions which are associated with the same "found-by-T0" bunch crossing
      // https://indico.cern.ch/event/1396220/#1-event-selection-with-its-rof
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 9.5 /* rejects collisions which are associated with the same "found-by-T0" bunch crossing*/);
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
        if (!isValidTrigger(track)) {
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
        if (!isValidTrigger(track))
          continue;
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
        if (doTriggPhysicalPrimary && !triggerTrack.mcPhysicalPrimary())
          continue;
        histos.fill(HIST("hTriggerPrimaryEtaVsPt"), track.pt(), track.eta(), collision.centFT0M());
        histos.fill(HIST("hTrackEtaVsPtVsPhi"), track.pt(), track.eta(), track.phi(), weight);
      }
    }
    for (auto const& assocTrack : assocHadrons) {
      auto assoc = assocTrack.track_as<TracksComplete>();
      if (!isValidAssocHadron(assoc))
        continue;
      float efficiency = 1.0f;
      float purity = 1.0f;
      if (efficiencyFlags.applyEfficiencyCorrection) {
        efficiency = hEfficiencyHadron->Interpolate(assoc.pt(), assoc.eta());
        if (efficiencyFlags.applyPurityHadron)
          purity = hPurityHadron->Interpolate(assoc.pt());
      }
      if (efficiency == 0) { // check for zero efficiency, do not apply if the case
        efficiency = 1;
      }
      float weight = efficiencyFlags.applyEfficiencyCorrection ? purity / efficiency : 1.0f;
      histos.fill(HIST("hAssocHadronsAllSelectedEtaVsPt"), assoc.pt(), assoc.eta(), collision.centFT0M(), weight);
      histos.fill(HIST("hAssocPtResolution"), assoc.pt(), assocTrack.mcOriginalPt());
      if (doAssocPhysicalPrimary && !assocTrack.mcPhysicalPrimary())
        continue;
      histos.fill(HIST("hAssocPrimaryEtaVsPt"), assoc.pt(), assoc.eta(), collision.centFT0M());
      histos.fill(HIST("hAsssocTrackEtaVsPtVsPhi"), assoc.pt(), assoc.eta(), assoc.phi(), weight);
    }

    // ________________________________________________
    // Do hadron - hadron correlations
    if (doFullCorrelationStudy)
      fillCorrelationsHadron(triggerTracks, assocHadrons, false, collision.posZ(), collision.centFT0M(), bField);
  }

  void processSameEventHV0s(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::PVMults>::iterator const& collision,
                            aod::AssocV0s const& associatedV0s, aod::TriggerTracks const& triggerTracks,
                            V0DatasWithoutTrackX const&, TracksComplete const&, aod::BCsWithTimestamps const&)
  {
    std::variant<BinningTypePP, BinningTypePbPb> colBinning =
      doPPAnalysis
        ? std::variant<BinningTypePP, BinningTypePbPb>{
            BinningTypePP{{axesConfigurations.axisVtxZ, axesConfigurations.axisMult}, true}}
        : std::variant<BinningTypePP, BinningTypePbPb>{BinningTypePbPb{{axesConfigurations.axisVtxZ, axesConfigurations.axisMult}, true}};

    double cent = doPPAnalysis ? collision.centFT0M() : collision.centFT0C();
    // ________________________________________________
    // skip if desired trigger not found
    if (triggerPresenceMap.size() > 0 && !TESTBIT(triggerPresenceMap[collision.globalIndex()], triggerBinToSelect)) {
      return;
    }

    // ________________________________________________
    // Perform basic event selection
    if (((doPPAnalysis && !isCollisionSelected(collision))) || (!doPPAnalysis && !isCollisionSelectedPbPb(collision, true))) {
      return;
    }
    // ________________________________________________
    if (!doprocessSameEventHCascades && doMixingQAandEventQA) {
      std::visit([&](auto const& binning) {
        histos.fill(HIST("MixingQA/hSECollisionBins"), binning.getBin({collision.posZ(), cent}));
      },
                 colBinning);
      histos.fill(HIST("EventQA/hMult"), cent);
      histos.fill(HIST("EventQA/hPvz"), collision.posZ());
      histos.fill(HIST("EventQA/hMultFT0vsTPC"), cent, collision.multNTracksPVeta1());
    }
    // Do basic QA
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    auto bField = getMagneticField(bc.timestamp());
    if (efficiencyFlags.applyEfficiencyCorrection) {
      initEfficiencyFromCCDB(bc);
    }
    TH2F* hEfficiencyV0[3];
    hEfficiencyV0[0] = hEfficiencyK0Short;
    hEfficiencyV0[1] = hEfficiencyLambda;
    hEfficiencyV0[2] = hEfficiencyAntiLambda;

    for (auto const& v0 : associatedV0s) {
      auto v0Data = v0.v0Core_as<V0DatasWithoutTrackX>();

      //---] track quality check [---
      auto postrack = v0Data.posTrack_as<TracksComplete>();
      auto negtrack = v0Data.negTrack_as<TracksComplete>();
      if (postrack.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRowsAssociated || negtrack.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRowsAssociated)
        continue;

      //---] syst cuts [---
      if (doPPAnalysis && (v0Data.v0radius() < systCuts.v0RadiusMin || v0Data.v0radius() > systCuts.v0RadiusMax ||
                           std::abs(v0Data.dcapostopv()) < systCuts.dcapostopv || std::abs(v0Data.dcanegtopv()) < systCuts.dcanegtopv ||
                           v0Data.v0cosPA() < systCuts.v0cospa || v0Data.dcaV0daughters() > systCuts.dcaV0dau))
        continue;
      if (!doPPAnalysis && !V0SelectedPbPb(v0Data))
        continue;
      uint64_t selMap = V0selectionBitmap(v0Data, collision.posX(), collision.posY(), collision.posZ());

      static_for<0, 2>([&](auto i) {
        constexpr int Index = i.value;
        float efficiency = 1.0f;
        if (efficiencyFlags.applyEfficiencyCorrection) {
          efficiency = hEfficiencyV0[Index]->Interpolate(v0Data.pt(), v0Data.eta());
        }
        if (efficiency == 0) { // check for zero efficiency, do not apply if the case
          efficiency = 1;
        }
        float weight = efficiencyFlags.applyEfficiencyCorrection ? 1. / efficiency : 1.0f;
        if (v0.compatible(Index, systCuts.dEdxCompatibility) && (!doMCassociation || v0.mcTrue(Index)) && (!doAssocPhysicalPrimary || v0.mcPhysicalPrimary()) && (!efficiencyFlags.applyEfficiencyCorrection || efficiency != 0)) {
          if ((TESTBIT(doCorrelation, Index)) && (doPPAnalysis || (TESTBIT(selMap, Index) && TESTBIT(selMap, Index + 3)))) {
            histos.fill(HIST("h3d") + HIST(kV0names[Index]) + HIST("Spectrum"), v0Data.pt(), cent, v0.invMassNSigma(Index), weight);
            if (std::abs(v0Data.rapidity(Index)) < ySel) {
              histos.fill(HIST("h3d") + HIST(kV0names[Index]) + HIST("SpectrumY"), v0Data.pt(), cent, v0.invMassNSigma(Index), weight);
            }
            if ((-massWindowConfigurations.maxBgNSigma < v0.invMassNSigma(Index) && v0.invMassNSigma(Index) < -massWindowConfigurations.minBgNSigma) || (+massWindowConfigurations.minBgNSigma < v0.invMassNSigma(Index) && v0.invMassNSigma(Index) < +massWindowConfigurations.maxBgNSigma)) {
              histos.fill(HIST("h") + HIST(kV0names[Index]) + HIST("EtaVsPtVsPhiBg"), v0Data.pt(), v0Data.eta(), v0Data.phi(), weight);
            }
            if (-massWindowConfigurations.maxPeakNSigma < v0.invMassNSigma(Index) && v0.invMassNSigma(Index) < +massWindowConfigurations.maxPeakNSigma) {
              histos.fill(HIST("h") + HIST(kV0names[Index]) + HIST("EtaVsPtVsPhi"), v0Data.pt(), v0Data.eta(), v0Data.phi(), weight);
            }
          }
        }
      });
    }
    if (!doprocessSameEventHCascades) {
      for (auto const& triggerTrack : triggerTracks) {
        auto track = triggerTrack.track_as<TracksComplete>();
        if (!isValidTrigger(track))
          continue;
        histos.fill(HIST("hTriggerAllSelectedEtaVsPt"), track.pt(), track.eta(), cent);
        histos.fill(HIST("hTriggerPtResolution"), track.pt(), triggerTrack.mcOriginalPt());
        if (doTriggPhysicalPrimary && !triggerTrack.mcPhysicalPrimary())
          continue;
        histos.fill(HIST("hTriggerPrimaryEtaVsPt"), track.pt(), track.eta(), cent);
        histos.fill(HIST("hTrackEtaVsPtVsPhi"), track.pt(), track.eta(), track.phi());
      }
    }

    // ________________________________________________
    // Do hadron - V0 correlations
    if (doFullCorrelationStudy)
      fillCorrelationsV0(triggerTracks, associatedV0s, false, collision.posX(), collision.posY(), collision.posZ(), cent, bField);
  }

  void processSameEventHCascades(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::PVMults>::iterator const& collision,
                                 aod::AssocV0s const&, aod::AssocCascades const& associatedCascades, aod::TriggerTracks const& triggerTracks,
                                 V0DatasWithoutTrackX const&, aod::CascDatas const&, TracksComplete const&, aod::BCsWithTimestamps const&)
  {
    std::variant<BinningTypePP, BinningTypePbPb> colBinning =
      doPPAnalysis
        ? std::variant<BinningTypePP, BinningTypePbPb>{
            BinningTypePP{{axesConfigurations.axisVtxZ, axesConfigurations.axisMult}, true}}
        : std::variant<BinningTypePP, BinningTypePbPb>{BinningTypePbPb{{axesConfigurations.axisVtxZ, axesConfigurations.axisMult}, true}};

    double cent = doPPAnalysis ? collision.centFT0M() : collision.centFT0C();
    // ________________________________________________
    // skip if desired trigger not found
    if (triggerPresenceMap.size() > 0 && !TESTBIT(triggerPresenceMap[collision.globalIndex()], triggerBinToSelect)) {
      return;
    }

    // ________________________________________________
    // Perform basic event selection
    if (((doPPAnalysis && !isCollisionSelected(collision))) || (!doPPAnalysis && !isCollisionSelectedPbPb(collision, true))) {
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
    TH2F* hEfficiencyCascade[4];
    hEfficiencyCascade[0] = hEfficiencyXiMinus;
    hEfficiencyCascade[1] = hEfficiencyXiPlus;
    hEfficiencyCascade[2] = hEfficiencyOmegaMinus;
    hEfficiencyCascade[3] = hEfficiencyOmegaPlus;

    for (auto const& casc : associatedCascades) {
      auto cascData = casc.cascData();

      //---] syst cuts [---
      if (doPPAnalysis && (std::abs(cascData.dcapostopv()) < systCuts.dcapostopv ||
                           std::abs(cascData.dcanegtopv()) < systCuts.dcanegtopv ||
                           std::abs(cascData.dcabachtopv()) < systCuts.cascDcabachtopv ||
                           cascData.dcaV0daughters() > systCuts.dcaV0dau ||
                           cascData.dcacascdaughters() > systCuts.cascDcacascdau ||
                           cascData.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < systCuts.v0cospa ||
                           cascData.casccosPA(collision.posX(), collision.posY(), collision.posZ()) < systCuts.cascCospa ||
                           cascData.cascradius() < systCuts.cascRadius ||
                           std::abs(cascData.dcav0topv(collision.posX(), collision.posY(), collision.posZ())) < systCuts.cascMindcav0topv ||
                           std::abs(cascData.mLambda() - o2::constants::physics::MassLambda0) > systCuts.cascV0masswindow))
        continue;
      if (!doPPAnalysis && !CascadeSelectedPbPb(cascData, collision.posX(), collision.posY(), collision.posZ()))
        continue;
      uint64_t CascselMap = CascadeselectionBitmap(cascData, collision.posX(), collision.posY(), collision.posZ());
      //---] track quality check [---
      auto postrack = cascData.posTrack_as<TracksComplete>();
      auto negtrack = cascData.negTrack_as<TracksComplete>();
      auto bachtrack = cascData.bachelor_as<TracksComplete>();
      if (postrack.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRowsAssociated || negtrack.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRowsAssociated || bachtrack.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRowsAssociated)
        continue;

      static_for<0, 3>([&](auto i) {
        constexpr int Index = i.value;
        if ((Index == 2 || Index == 3) && casc.compatible(Index, systCuts.dEdxCompatibility) && std::abs(casc.invMassNSigma(Index - 2)) < massWindowConfigurations.nSigmaNearXiMassCenter) {
          return;
        }
        float efficiency = 1.0f;
        if (efficiencyFlags.applyEfficiencyCorrection) {
          efficiency = hEfficiencyCascade[Index]->Interpolate(cascData.pt(), cascData.eta());
        }
        if (efficiency == 0) { // check for zero efficiency, do not apply if the case
          efficiency = 1;
        }
        float weight = efficiencyFlags.applyEfficiencyCorrection ? 1. / efficiency : 1.0f;
        if (casc.compatible(Index, systCuts.dEdxCompatibility) && (!doMCassociation || casc.mcTrue(Index)) && (!doAssocPhysicalPrimary || casc.mcPhysicalPrimary()) && (!efficiencyFlags.applyEfficiencyCorrection || efficiency != 0)) {
          if (TESTBIT(doCorrelation, Index + 3) && (doPPAnalysis || (TESTBIT(CascselMap, Index) && TESTBIT(CascselMap, Index + 4) && TESTBIT(CascselMap, Index + 8) && TESTBIT(CascselMap, Index + 12)))) {
            histos.fill(HIST("h3d") + HIST(kCascadenames[Index]) + HIST("Spectrum"), cascData.pt(), cent, casc.invMassNSigma(Index), weight);
            if (std::abs(cascData.rapidity(Index)) < ySel) {
              histos.fill(HIST("h3d") + HIST(kCascadenames[Index]) + HIST("SpectrumY"), cascData.pt(), cent, casc.invMassNSigma(Index), weight);
            }
            if (-massWindowConfigurations.maxPeakNSigma < casc.invMassNSigma(Index) && casc.invMassNSigma(Index) < +massWindowConfigurations.maxPeakNSigma) {
              histos.fill(HIST("h") + HIST(kCascadenames[Index]) + HIST("EtaVsPtVsPhi"), cascData.pt(), cascData.eta(), cascData.phi(), weight);
            }
            if ((-massWindowConfigurations.maxBgNSigma < casc.invMassNSigma(Index) && casc.invMassNSigma(Index) < -massWindowConfigurations.minBgNSigma) || (+massWindowConfigurations.minBgNSigma < casc.invMassNSigma(Index) && casc.invMassNSigma(Index) < +massWindowConfigurations.maxBgNSigma)) {
              histos.fill(HIST("h") + HIST(kCascadenames[Index]) + HIST("EtaVsPtVsPhiBg"), cascData.pt(), cascData.eta(), cascData.phi(), weight);
            }
          }
        }
      });
    }
    for (auto const& triggerTrack : triggerTracks) {
      auto track = triggerTrack.track_as<TracksComplete>();
      if (!isValidTrigger(track))
        continue;
      histos.fill(HIST("hTriggerAllSelectedEtaVsPt"), track.pt(), track.eta(), cent);
      histos.fill(HIST("hTriggerPtResolution"), track.pt(), triggerTrack.mcOriginalPt());
      if (doTriggPhysicalPrimary && !triggerTrack.mcPhysicalPrimary())
        continue;
      histos.fill(HIST("hTriggerPrimaryEtaVsPt"), track.pt(), track.eta(), cent);
      if (track.sign() > 0)
        histos.fill(HIST("hPositiveTriggerPrimaryEtaVsPt"), track.pt(), track.eta(), cent);
      else
        histos.fill(HIST("hNegativeTriggerPrimaryEtaVsPt"), track.pt(), track.eta(), cent);
      histos.fill(HIST("hTrackEtaVsPtVsPhi"), track.pt(), track.eta(), track.phi());
    }

    // ________________________________________________
    // Do hadron - cascade correlations
    if (doFullCorrelationStudy)
      fillCorrelationsCascade(triggerTracks, associatedCascades, false, collision.posX(), collision.posY(), collision.posZ(), cent, bField);
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
      if (!isValidAssocHadron(pionTrack))
        continue;

      histos.fill(HIST("hPionEtaVsPtAllSelected"), pionTrack.pt(), pionTrack.eta(), collision.centFT0M());
      if (doAssocPhysicalPrimary && !pion.mcPhysicalPrimary())
        continue;
      if (doMCassociation && std::abs(pion.pdgCode()) != 211)
        continue;
      histos.fill(HIST("hPionEtaVsPt"), pionTrack.pt(), pionTrack.eta(), collision.centFT0M());
      if (pionTrack.sign() > 0)
        histos.fill(HIST("hPositivePionEtaVsPt"), pionTrack.pt(), pionTrack.eta(), collision.centFT0M());
      else
        histos.fill(HIST("hNegativePionEtaVsPt"), pionTrack.pt(), pionTrack.eta(), collision.centFT0M());
    }
    if (!doprocessSameEventHCascades && !doprocessSameEventHV0s) {
      for (auto const& triggerTrack : triggerTracks) {
        auto track = triggerTrack.track_as<TracksComplete>();
        if (!isValidTrigger(track))
          continue;
        histos.fill(HIST("hTriggerAllSelectedEtaVsPt"), track.pt(), track.eta(), collision.centFT0M());
        histos.fill(HIST("hTriggerPtResolution"), track.pt(), triggerTrack.mcOriginalPt());
        if (doTriggPhysicalPrimary && !triggerTrack.mcPhysicalPrimary())
          continue;
        histos.fill(HIST("hTriggerPrimaryEtaVsPt"), track.pt(), track.eta(), collision.centFT0M());
        histos.fill(HIST("hTrackEtaVsPtVsPhi"), track.pt(), track.eta(), track.phi());
      }
    }

    // ________________________________________________
    // Do hadron - Pion correlations
    if (doFullCorrelationStudy)
      fillCorrelationsHadron(triggerTracks, associatedPions, false, collision.posZ(), collision.centFT0M(), bField);
  }

  void processMixedEventHHadrons(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults> const& collisions,
                                 aod::AssocHadrons const& assocHadrons, aod::TriggerTracks const& triggerTracks,
                                 TracksComplete const&, aod::BCsWithTimestamps const&)
  {
    BinningTypePP colBinning{{axesConfigurations.axisVtxZ, axesConfigurations.axisMult}, true};
    for (auto const& [collision1, collision2] : soa::selfCombinations(colBinning, mixingParameter, -1, collisions, collisions)) {
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
      if (collision1.centFT0M() > axisRanges[5][1] || collision1.centFT0M() < axisRanges[5][0])
        continue;
      if (collision2.centFT0M() > axisRanges[5][1] || collision2.centFT0M() < axisRanges[5][0])
        continue;
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
      if (doFullCorrelationStudy)
        fillCorrelationsHadron(slicedTriggerTracks, slicedAssocHadrons, true, collision1.posZ(), collision1.centFT0M(), bField);
    }
  }

  void processMixedEventHV0s(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::PVMults> const& collisions,
                             aod::AssocV0s const& associatedV0s, aod::TriggerTracks const& triggerTracks,
                             V0DatasWithoutTrackX const&, TracksComplete const&, aod::BCsWithTimestamps const&)
  {
    std::variant<BinningTypePP, BinningTypePbPb> colBinning =
      doPPAnalysis
        ? std::variant<BinningTypePP, BinningTypePbPb>{
            BinningTypePP{{axesConfigurations.axisVtxZ, axesConfigurations.axisMult}, true}}
        : std::variant<BinningTypePP, BinningTypePbPb>{BinningTypePbPb{{axesConfigurations.axisVtxZ, axesConfigurations.axisMult}, true}};

    std::visit([&](auto const& binning) {
      for (auto const& [collision1, collision2] : soa::selfCombinations(binning, mixingParameter, -1, collisions, collisions)) {
        double cent1 = doPPAnalysis ? collision1.centFT0M() : collision1.centFT0C();
        double cent2 = doPPAnalysis ? collision2.centFT0M() : collision2.centFT0C();
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
        if ((doPPAnalysis && (!isCollisionSelected(collision1) || !isCollisionSelected(collision2))) || (!doPPAnalysis && (!isCollisionSelectedPbPb(collision1, true) || (!isCollisionSelectedPbPb(collision2, true))))) {
          continue;
        }
        if (cent1 > axisRanges[5][1] || cent1 < axisRanges[5][0])
          continue;
        if (cent2 > axisRanges[5][1] || cent2 < axisRanges[5][0])
          continue;

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
        if (doFullCorrelationStudy)
          fillCorrelationsV0(slicedTriggerTracks, slicedAssocV0s, true, collision1.posX(), collision1.posY(), collision1.posZ(), cent1, bField);
      }
    },
               colBinning);
  }

  void processMixedEventHCascades(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::PVMults> const& collisions,
                                  aod::AssocV0s const&, aod::AssocCascades const& associatedCascades, aod::TriggerTracks const& triggerTracks,
                                  V0DatasWithoutTrackX const&, aod::CascDatas const&, TracksComplete const&, aod::BCsWithTimestamps const&)
  {
    std::variant<BinningTypePP, BinningTypePbPb> colBinning =
      doPPAnalysis
        ? std::variant<BinningTypePP, BinningTypePbPb>{
            BinningTypePP{{axesConfigurations.axisVtxZ, axesConfigurations.axisMult}, true}}
        : std::variant<BinningTypePP, BinningTypePbPb>{BinningTypePbPb{{axesConfigurations.axisVtxZ, axesConfigurations.axisMult}, true}};

    std::visit([&](auto const& binning) {
      for (auto const& [collision1, collision2] : soa::selfCombinations(binning, mixingParameter, -1, collisions, collisions)) {
        double cent1 = doPPAnalysis ? collision1.centFT0M() : collision1.centFT0C();
        double cent2 = doPPAnalysis ? collision2.centFT0M() : collision2.centFT0C();
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
        if ((doPPAnalysis && (!isCollisionSelected(collision1) || !isCollisionSelected(collision2))) || (!doPPAnalysis && (!isCollisionSelectedPbPb(collision1, true) || (!isCollisionSelectedPbPb(collision2, true))))) {
          continue;
        }
        if (cent1 > axisRanges[5][1] || cent1 < axisRanges[5][0])
          continue;
        if (cent2 > axisRanges[5][1] || cent2 < axisRanges[5][0])
          continue;
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
        if (doFullCorrelationStudy)
          fillCorrelationsCascade(slicedTriggerTracks, slicedAssocCascades, true, collision1.posX(), collision1.posY(), collision1.posZ(), cent1, bField);
      }
    },
               colBinning);
  }

  void processMixedEventHPions(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults> const& collisions,
                               soa::Join<aod::AssocHadrons, aod::AssocPID> const& assocPions, soa::Join<aod::TriggerTracks, aod::TriggerTrackExtras> const& triggerTracks,
                               TracksComplete const&, aod::BCsWithTimestamps const&)
  {
    BinningTypePP colBinning{{axesConfigurations.axisVtxZ, axesConfigurations.axisMult}, true};
    for (auto const& [collision1, collision2] : soa::selfCombinations(colBinning, mixingParameter, -1, collisions, collisions)) {
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
      if (collision1.centFT0M() > axisRanges[5][1] || collision1.centFT0M() < axisRanges[5][0])
        continue;
      if (collision2.centFT0M() > axisRanges[5][1] || collision2.centFT0M() < axisRanges[5][0])
        continue;
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
      if (doFullCorrelationStudy)
        fillCorrelationsHadron(slicedTriggerTracks, slicedAssocPions, true, collision1.posZ(), collision1.centFT0M(), bField);
    }
  }

  void processMCGenerated(aod::McCollision const& /*mcCollision*/, soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults>> const& collisions, aod::McParticles const& mcParticles)
  {
    histos.fill(HIST("hClosureTestEventCounter"), 2.5f);

    for (auto const& mcParticle : mcParticles) {
      double geta = mcParticle.eta();
      if (std::abs(geta) > etaSel) {
        continue;
      }
      double gpt = mcParticle.pt();
      if (std::abs(mcParticle.pdgCode()) == PDG_t::kPiPlus || std::abs(mcParticle.pdgCode()) == PDG_t::kKPlus || std::abs(mcParticle.pdgCode()) == PDG_t::kProton || std::abs(mcParticle.pdgCode()) == PDG_t::kElectron || std::abs(mcParticle.pdgCode()) == PDG_t::kMuonMinus) {
        if (!doTriggPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
          histos.fill(HIST("hGeneratedQAPtTrigger"), gpt, 0.0f); // step 1: before all selections
        }
      }

      if (!doAssocPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
        if (mcParticle.pdgCode() == PDG_t::kK0Short && doCorrelationK0Short) {
          histos.fill(HIST("hGeneratedQAPtAssociatedK0"), gpt, 0.0f); // step 1: before all selections
        }
      }
    }

    for (auto const& mcParticle : mcParticles) {
      if (doAssocPhysicalPrimaryInGen && !mcParticle.isPhysicalPrimary())
        continue;
      static_for<0, 7>([&](auto i) {
        constexpr int Index = i.value;
        if (i == 0 || i == 7) {
          if (std::abs(mcParticle.pdgCode()) == kPdgCodes[i]) {
            histos.fill(HIST("Generated/h") + HIST(kParticlenames[Index]), mcParticle.pt(), mcParticle.eta(), 1);
            if (i == 7 && mcParticle.pdgCode() > 0)
              histos.fill(HIST("Generated/hPositive") + HIST(kParticlenames[Index]), mcParticle.pt(), mcParticle.eta(), 1);
            else if (i == 7 && mcParticle.pdgCode() < 0)
              histos.fill(HIST("Generated/hNegative") + HIST(kParticlenames[Index]), mcParticle.pt(), mcParticle.eta(), 1);
          }
        } else {
          if (mcParticle.pdgCode() == kPdgCodes[i])
            histos.fill(HIST("Generated/h") + HIST(kParticlenames[Index]), mcParticle.pt(), mcParticle.eta(), 1);
        }
      });
    }
    if (collisions.size() < 1)
      return;

    // determine best collision properties
    int biggestNContribs = -1;
    int bestCollisionFT0Mpercentile = -1;
    float bestCollisionVtxZ = 0.0f;
    bool bestCollisionSel8 = false;
    bool bestCollisionINELgtZERO = false;
    uint32_t bestCollisionTriggerPresenceMap = 0;

    for (auto const& collision : collisions) {
      if (biggestNContribs < collision.numContrib()) {
        biggestNContribs = collision.numContrib();
        bestCollisionFT0Mpercentile = collision.centFT0M();
        bestCollisionSel8 = collision.sel8();
        bestCollisionVtxZ = collision.posZ();
        bestCollisionINELgtZERO = collision.isInelGt0();
        if (triggerPresenceMap.size() > 0)
          bestCollisionTriggerPresenceMap = triggerPresenceMap[collision.globalIndex()];
      }
    }

    if (collisions.size() > 1) {
      for (auto const& mcParticle : mcParticles) {
        if (doAssocPhysicalPrimaryInGen && !mcParticle.isPhysicalPrimary())
          continue;
        if (std::abs(mcParticle.y()) > ySel)
          continue;
        static_for<0, 7>([&](auto i) {
          constexpr int Index = i.value;
          if (i == 0 || i == 7) {
            if (std::abs(mcParticle.pdgCode()) == kPdgCodes[i]) {
              histos.fill(HIST("GeneratedWithPV/h") + HIST(kParticlenames[Index]) + HIST("_MidYVsMult_TwoPVsOrMore"), mcParticle.pt(), bestCollisionFT0Mpercentile);
            }
          } else {
            if (mcParticle.pdgCode() == kPdgCodes[i])
              histos.fill(HIST("GeneratedWithPV/h") + HIST(kParticlenames[Index]) + HIST("_MidYVsMult_TwoPVsOrMore"), mcParticle.pt(), bestCollisionFT0Mpercentile);
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
    if (!bestCollisionSel8)
      return;
    if (std::abs(bestCollisionVtxZ) > 10.0f)
      return;
    if (!bestCollisionINELgtZERO)
      return;

    histos.fill(HIST("hClosureTestEventCounter"), 3.5f);

    for (auto const& mcParticle : mcParticles) {
      double geta = mcParticle.eta();
      if (std::abs(geta) > etaSel) {
        continue;
      }
      double gpt = mcParticle.pt();
      if (std::abs(mcParticle.pdgCode()) == PDG_t::kPiPlus || std::abs(mcParticle.pdgCode()) == PDG_t::kKPlus || std::abs(mcParticle.pdgCode()) == PDG_t::kProton || std::abs(mcParticle.pdgCode()) == PDG_t::kElectron || std::abs(mcParticle.pdgCode()) == PDG_t::kMuonMinus) {
        if (!doTriggPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
          histos.fill(HIST("hGeneratedQAPtTrigger"), gpt, 1.0f); // step 2: after event selection
        }
      }

      if (!doAssocPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
        if (mcParticle.pdgCode() == PDG_t::kK0Short && doCorrelationK0Short) {
          histos.fill(HIST("hGeneratedQAPtAssociatedK0"), gpt, 1.0f); // step 2: before all selections
        }
      }
    }

    for (auto const& mcParticle : mcParticles) {
      if (doAssocPhysicalPrimaryInGen && !mcParticle.isPhysicalPrimary()) {
        continue;
      }
      double geta = mcParticle.eta();
      double gpt = mcParticle.pt();
      if (std::abs(mcParticle.pdgCode()) == PDG_t::kPiPlus || std::abs(mcParticle.pdgCode()) == PDG_t::kKPlus || std::abs(mcParticle.pdgCode()) == PDG_t::kProton || std::abs(mcParticle.pdgCode()) == PDG_t::kElectron || std::abs(mcParticle.pdgCode()) == PDG_t::kMuonMinus) {
        histos.fill(HIST("GeneratedWithPV/hTrigger"), gpt, geta, bestCollisionFT0Mpercentile);
        if (mcParticle.pdgCode() > 0)
          histos.fill(HIST("GeneratedWithPV/hPositiveTrigger"), gpt, geta, bestCollisionFT0Mpercentile);
        else
          histos.fill(HIST("GeneratedWithPV/hNegativeTrigger"), gpt, geta, bestCollisionFT0Mpercentile);
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
            if (lamParticleMother.pdgCode() == PDG_t::kXiMinus) // Xi Minus Mother Matched
            {
              histos.fill(HIST("GeneratedWithPV/hLambdaFromXiMinus"), gpt, geta);
            }
            if (lamParticleMother.pdgCode() == o2::constants::physics::Pdg::kXi0) // Xi Zero Mother Matched
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
            if (lamParticleMother.pdgCode() == PDG_t::kXiPlusBar) {
              histos.fill(HIST("GeneratedWithPV/hAntiLambdaFromXiPlus"), gpt, geta);
            }
            if (lamParticleMother.pdgCode() == -o2::constants::physics::Pdg::kXi0) // Xi Zero Mother Matched
            {
              histos.fill(HIST("GeneratedWithPV/hAntiLambdaFromXiZero"), gpt, geta);
            }
          }
        }
      }
      static_for<0, 7>([&](auto i) {
        constexpr int Index = i.value;
        if (i == 0 || i == 7) {
          if (std::abs(mcParticle.pdgCode()) == kPdgCodes[i]) {
            histos.fill(HIST("GeneratedWithPV/h") + HIST(kParticlenames[Index]), gpt, geta, bestCollisionFT0Mpercentile);
            if (std::abs(mcParticle.y()) < ySel)
              histos.fill(HIST("GeneratedWithPV/h") + HIST(kParticlenames[Index]) + HIST("_MidYVsMult"), gpt, bestCollisionFT0Mpercentile);
            if (i == 7 && mcParticle.pdgCode() > 0)
              histos.fill(HIST("GeneratedWithPV/hPositive") + HIST(kParticlenames[Index]), mcParticle.pt(), mcParticle.eta(), bestCollisionFT0Mpercentile);
            else if (i == 7 && mcParticle.pdgCode() < 0)
              histos.fill(HIST("GeneratedWithPV/hNegative") + HIST(kParticlenames[Index]), mcParticle.pt(), mcParticle.eta(), bestCollisionFT0Mpercentile);
          }

        } else {
          if (mcParticle.pdgCode() == kPdgCodes[i]) {
            histos.fill(HIST("GeneratedWithPV/h") + HIST(kParticlenames[Index]), gpt, geta, bestCollisionFT0Mpercentile);
            if (std::abs(mcParticle.y()) < ySel)
              histos.fill(HIST("GeneratedWithPV/h") + HIST(kParticlenames[Index]) + HIST("_MidYVsMult"), gpt, bestCollisionFT0Mpercentile);
          }
        }
      });
    }
  }
  void processClosureTest(aod::McCollision const& /*mcCollision*/, soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults>> const& recCollisions, aod::McParticles const& mcParticles)
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

    for (auto const& mcParticle : mcParticles) {
      double geta = mcParticle.eta();
      if (std::abs(geta) > etaSel) {
        continue;
      }
      double gpt = mcParticle.pt();
      if (std::abs(mcParticle.pdgCode()) == PDG_t::kPiPlus || std::abs(mcParticle.pdgCode()) == PDG_t::kKPlus || std::abs(mcParticle.pdgCode()) == PDG_t::kProton || std::abs(mcParticle.pdgCode()) == PDG_t::kElectron || std::abs(mcParticle.pdgCode()) == PDG_t::kMuonMinus) {
        if (!doTriggPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
          histos.fill(HIST("hClosureQAPtTrigger"), gpt, 0.0f); // step 1: no event selection whatsoever
        }
      }

      if (!doAssocPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
        if (mcParticle.pdgCode() == PDG_t::kK0Short && doCorrelationK0Short) {
          histos.fill(HIST("hClosureQAPtAssociatedK0"), gpt, 0.0f); // step 1: no event selection whatsoever
        }
      }
    }

    histos.fill(HIST("hClosureTestEventCounter"), 0.5f);

    int bestCollisionFT0Mpercentile = -1;
    float bestCollisionVtxZ = 0.0f;
    bool bestCollisionSel8 = false;
    bool bestCollisionINELgtZERO = false;
    int biggestNContribs = -1;
    uint32_t bestCollisionTriggerPresenceMap = 0;

    for (auto const& recCollision : recCollisions) {
      if (biggestNContribs < recCollision.numContrib()) {
        biggestNContribs = recCollision.numContrib();
        bestCollisionFT0Mpercentile = recCollision.centFT0M();
        bestCollisionSel8 = recCollision.sel8();
        bestCollisionVtxZ = recCollision.posZ();
        bestCollisionINELgtZERO = recCollision.isInelGt0();
        if (triggerPresenceMap.size() > 0)
          bestCollisionTriggerPresenceMap = triggerPresenceMap[recCollision.globalIndex()];
      }
    }
    // ________________________________________________
    // skip if desired trigger not found
    if (triggerPresenceMap.size() > 0 && !TESTBIT(bestCollisionTriggerPresenceMap, triggerBinToSelect)) {
      return;
    }

    if (doGenEventSelection) {
      if (!bestCollisionSel8)
        return;
      if (std::abs(bestCollisionVtxZ) > zVertexCut)
        return;
      if (!bestCollisionINELgtZERO)
        return;
      if (bestCollisionFT0Mpercentile > axisRanges[5][1] || bestCollisionFT0Mpercentile < axisRanges[5][0]) {
        return;
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
        if (!doTriggPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
          histos.fill(HIST("hClosureQAPtTrigger"), gpt, 1.0f); // step 2: after event selection
        }
      }

      if (!doAssocPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
        if (mcParticle.pdgCode() == PDG_t::kK0Short && doCorrelationK0Short) {
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
        if (!doTriggPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
          triggerIndices.emplace_back(iteratorNum);
          histos.fill(HIST("ClosureTest/hTrigger"), gpt, geta, bestCollisionFT0Mpercentile);
        }
        if (doCorrelationHadron) {
          if (!doAssocPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
            assocHadronIndices.emplace_back(iteratorNum);
            histos.fill(HIST("ClosureTest/hHadron"), gpt, geta, gphi);
          }
        }
      }
      if (!doAssocPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
        if (std::abs(mcParticle.pdgCode()) == PDG_t::kPiPlus && doCorrelationPion) {
          piIndices.emplace_back(iteratorNum);
          histos.fill(HIST("ClosureTest/hPion"), gpt, geta, gphi);
        }
        if (mcParticle.pdgCode() == PDG_t::kK0Short && doCorrelationK0Short) {
          k0ShortIndices.emplace_back(iteratorNum);
          histos.fill(HIST("ClosureTest/hK0Short"), gpt, geta, gphi);
        }
        if (mcParticle.pdgCode() == PDG_t::kLambda0 && doCorrelationLambda) {
          lambdaIndices.emplace_back(iteratorNum);
          histos.fill(HIST("ClosureTest/hLambda"), gpt, geta, gphi);
        }
        if (mcParticle.pdgCode() == PDG_t::kLambda0Bar && doCorrelationAntiLambda) {
          antiLambdaIndices.emplace_back(iteratorNum);
          histos.fill(HIST("ClosureTest/hAntiLambda"), gpt, geta, gphi);
        }
        if (mcParticle.pdgCode() == PDG_t::kXiMinus && doCorrelationXiMinus) {
          xiMinusIndices.emplace_back(iteratorNum);
          histos.fill(HIST("ClosureTest/hXiMinus"), gpt, geta, gphi);
        }
        if (mcParticle.pdgCode() == PDG_t::kXiPlusBar && doCorrelationXiPlus) {
          xiPlusIndices.emplace_back(iteratorNum);
          histos.fill(HIST("ClosureTest/hXiPlus"), gpt, geta, gphi);
        }
        if (mcParticle.pdgCode() == PDG_t::kOmegaMinus && doCorrelationOmegaMinus) {
          omegaMinusIndices.emplace_back(iteratorNum);
          histos.fill(HIST("ClosureTest/hOmegaMinus"), gpt, geta, gphi);
        }
        if (mcParticle.pdgCode() == PDG_t::kOmegaPlusBar && doCorrelationOmegaPlus) {
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
            if (deltaphi < axisRanges[0][0] || deltaphi > axisRanges[0][1])
              continue;
            if (deltaeta < axisRanges[1][0] || deltaeta > axisRanges[1][1])
              continue;
            if (ptassoc < axisRanges[2][0] || ptassoc > axisRanges[2][1])
              continue;
            if (TESTBIT(doCorrelation, i))
              histos.fill(HIST("ClosureTest/sameEvent/") + HIST(kParticlenames[Index]), computeDeltaPhi(gphitrigger, gphiassoc), deltaeta, ptassoc, pttrigger, bestCollisionVtxZ, bestCollisionFT0Mpercentile);
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
      if (postrack.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRowsAssociated || negtrack.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRowsAssociated)
        continue;

      //---] syst cuts [---
      if (v0Data.v0radius() < systCuts.v0RadiusMin || v0Data.v0radius() > systCuts.v0RadiusMax ||
          std::abs(v0Data.dcapostopv()) < systCuts.dcapostopv || std::abs(v0Data.dcanegtopv()) < systCuts.dcanegtopv ||
          v0Data.v0cosPA() < systCuts.v0cospa || v0Data.dcaV0daughters() > systCuts.dcaV0dau)
        continue;

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
  PROCESS_SWITCH(HStrangeCorrelation, processSelectEventWithTrigger, "Select events with trigger only", true);
  PROCESS_SWITCH(HStrangeCorrelation, processSameEventHV0s, "Process same events, h-V0s", true);
  PROCESS_SWITCH(HStrangeCorrelation, processSameEventHCascades, "Process same events, h-Cascades", true);
  PROCESS_SWITCH(HStrangeCorrelation, processSameEventHPions, "Process same events, h-Pion", true);
  PROCESS_SWITCH(HStrangeCorrelation, processSameEventHHadrons, "Process same events, h-h", true);

  PROCESS_SWITCH(HStrangeCorrelation, processMixedEventHV0s, "Process mixed events, h-V0s", true);
  PROCESS_SWITCH(HStrangeCorrelation, processMixedEventHCascades, "Process mixed events, h-Cascades", true);
  PROCESS_SWITCH(HStrangeCorrelation, processMixedEventHPions, "Process mixed events, h-Pion", true);
  PROCESS_SWITCH(HStrangeCorrelation, processMixedEventHHadrons, "Process mixed events, h-h", true);

  PROCESS_SWITCH(HStrangeCorrelation, processMCGenerated, "Process MC generated", false);
  PROCESS_SWITCH(HStrangeCorrelation, processClosureTest, "Process Closure Test", false);
  PROCESS_SWITCH(HStrangeCorrelation, processFeedDown, "process Feed Down", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HStrangeCorrelation>(cfgc)};
}
