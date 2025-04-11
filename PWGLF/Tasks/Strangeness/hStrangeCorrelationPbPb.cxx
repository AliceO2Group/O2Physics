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
/// \file hStrangeCorrelationPbPb.cxx
/// \brief This task serves to do hadron-(strange hadron) correlation studies.
///  The yield will be calculated using the two-particle correlation method.
///  Trigger particle : Hadrons
///  Associated Particles : V0s or Cascades
///  this task requires the hStrangeCorrelationFilter to have been run before.
///
/// \author YaZhen Lin (yalin@cern.ch)

#include <string>
#include <vector>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFHStrangeCorrelationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Framework/StaticFor.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsTPC/TrackTPC.h"
#include "EventFiltering/Zorro.h"
#include "EventFiltering/ZorroSummary.h"
#include <TPDGCode.h>

using namespace o2;
using namespace o2::constants::math;
using namespace o2::framework;
using namespace o2::framework::expressions;

// using TracksComplete = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA>;
using V0DatasWithoutTrackX = soa::Join<aod::V0Indices, aod::V0Cores>;
using TracksComplete = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::TracksDCA>;
// using DauTracks = soa::Join<aod::DauTrackExtras, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr>;
using TrackTPC = o2::tpc::TrackTPC;

struct HStrangeCorrelation {
  // for efficiency corrections if requested
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  Service<o2::framework::O2DatabasePDG> pdgDB;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // event filtering
  Configurable<string> zorroMask{"zorroMask", "", "zorro trigger class to select on (empty: none)"};

  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  Configurable<bool> doCorrelationK0Short{"doCorrelationK0Short", true, "do K0Short correlation"};
  Configurable<bool> doCorrelationLambda{"doCorrelationLambda", false, "do Lambda correlation"};
  Configurable<bool> doCorrelationAntiLambda{"doCorrelationAntiLambda", false, "do AntiLambda correlation"};
  Configurable<bool> doCorrelationXiMinus{"doCorrelationXiMinus", false, "do XiMinus correlation"};
  Configurable<bool> doCorrelationXiPlus{"doCorrelationXiPlus", false, "do XiPlus correlation"};
  Configurable<bool> doCorrelationOmegaMinus{"doCorrelationOmegaMinus", false, "do OmegaMinus correlation"};
  Configurable<bool> doCorrelationOmegaPlus{"doCorrelationOmegaPlus", false, "do OmegaPlus correlation"};
  Configurable<bool> doGenEventSelection{"doGenEventSelection", true, "use event selections when performing closure test for the gen events"};

  Configurable<bool> skipUnderOverflowInTHn{"skipUnderOverflowInTHn", false, "skip under/overflow in THns"};
  Configurable<int> mixingParameter{"mixingParameter", 10, "how many events are mixed"};
  Configurable<bool> doMCassociation{"doMCassociation", false, "fill everything only for MC associated"};
  Configurable<bool> doTriggPhysicalPrimary{"doTriggPhysicalPrimary", false, "require physical primary for trigger particles"};
  Configurable<bool> doAssocPhysicalPrimary{"doAssocPhysicalPrimary", false, "require physical primary for associated particles"};
  Configurable<bool> doLambdaPrimary{"doLambdaPrimary", false, "do primary selection for lambda"};
  Configurable<bool> doAutocorrelationRejection{"doAutocorrelationRejection", true, "reject pairs where trigger Id is the same as daughter particle Id"};

  Configurable<int> triggerBinToSelect{"triggerBinToSelect", 0, "trigger bin to select on if processSelectEventWithTrigger enabled"};
  Configurable<int> triggerParticleCharge{"triggerParticleCharge", 0, "For checks, if 0 all charged tracks, if -1 only neg., if 1 only positive"};

  // used for event selection
  Configurable<bool> doEventSelected{"doEventSelected", true, "do event selection"};
  Configurable<bool> doCentSelectedhigh{"doCentSelectedhigh", true, "do centrality selection 0/%-5/%"};
  Configurable<bool> doCentSelectedmed{"doCentSelectedmed", false, "do centrality selection 5/%-10/%"};
  Configurable<bool> doCentSelectedlow{"doCentSelectedlow", false, "do centrality selection 60/%-80/%"};
  Configurable<float> zVertexCut{"zVertexCut", 10, "Cut on PV position"};
  Configurable<int> cfgCutOccupancyHigh{"cfgCutOccupancyHigh", 3000, "High cut on TPC occupancy"};
  Configurable<int> cfgCutOccupancyLow{"cfgCutOccupancyLow", 0, "Low cut on TPC occupancy"};
  Configurable<int> centhighMin{"centhighMin", 0, "centrality selection 0/%-5/% Minimal accepted centrality"};
  Configurable<int> centhighMax{"centhighMax", 5, "centrality selection 0/%-5/% Maximal accepted centrality"};
  Configurable<int> centmedMin{"centmedMin", 5, "centrality selection 5/%-10/% Minimal accepted centrality"};
  Configurable<int> centmedMax{"centmedMax", 10, "centrality selection 5/%-10/% Maximal accepted centrality"};
  Configurable<int> centlowMin{"centlowMin", 60, "centrality selection 60/%-80/% Minimal accepted centrality"};
  Configurable<int> centlowMax{"centlowMax", 80, "centrality selection 60/%-80/% Maximal accepted centrality"};

  // Axes - configurable for smaller sizes
  ConfigurableAxis axisMult{"axisMult", {VARIABLE_WIDTH, 0.0f, 0.01f, 1.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 70.0f, 100.0f}, "Mixing bins - multiplicity"};
  ConfigurableAxis axisVtxZ{"axisVtxZ", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis axisPhi{"axisPhi", {72, 0, TwoPI}, "#phi"};
  ConfigurableAxis axisEta{"axisEta", {40, -0.8, +0.8}, "#eta"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {72, -PIHalf, PIHalf * 3}, "delta #varphi axis for histograms"};
  ConfigurableAxis axisDeltaEta{"axisDeltaEta", {80, -1.6, 1.6}, "delta eta axis for histograms"};
  ConfigurableAxis axisPtAssoc{"axisPtAssoc", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt associated axis for histograms"};
  ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 10.0, 12.0, 15.0, 20.0, 30.0, 40.0, 70.0, 100}, "pt associated axis for histograms"};
  ConfigurableAxis axisPtQA{"axisPtQA", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for QA histograms"};
  ConfigurableAxis axisMultCount{"axisMultCount", {VARIABLE_WIDTH, 0, 200, 400, 600, 800, 1000, 1400, 1800, 2300, 2800, 3300, 4000, 5000, 6000}, "Mixing bins - multiplicity"};
  ConfigurableAxis axisMassNSigma{"axisMassNSigma", {40, -2, 2}, "Axis for mass Nsigma"};
  ConfigurableAxis axisAlpha{"axisAlpha", {40, -1, 1}, "Axis for armpodCut alpha"};
  ConfigurableAxis axisptarm{"axisptarm", {40, 0, 0.3}, "Axis for armpodCut ptarm"};

  // for topo var QA
  struct : ConfigurableGroup {
    Configurable<float> maxPeakNSigma{"maxPeakNSigma", 5, "Peak region edge definition (in sigma)"};
    Configurable<float> minBgNSigma{"minBgNSigma", 5, "Bg region edge closest to peak (in sigma)"};
    Configurable<float> maxBgNSigma{"maxBgNSigma", 10, "Bg region edge furthest to peak (in sigma)"};
  } massWindowConfigurations; // allows for gap between peak and bg in case someone wants to

  // Implementation of on-the-spot efficiency correction
  Configurable<bool> applyEfficiencyCorrection{"applyEfficiencyCorrection", false, "apply efficiency correction"};
  Configurable<bool> applyEfficiencyForTrigger{"applyEfficiencyForTrigger", false, "apply efficiency correction for the trigger particle"};
  Configurable<std::string> ccdburl{"ccdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository to use"};
  Configurable<std::string> efficiencyCCDBPath{"efficiencyCCDBPath", "GLO/Config/GeometryAligned", "Path of the efficiency corrections"};

  // Configurables for doing subwagon systematics
  // Group all settings necessary for systematics in a specific ConfigurableGroup
  struct : ConfigurableGroup {
    std::string prefix = "systematics";

    // Additional selection on the AP plot (exclusive for K0Short)
    // original equation: lArmPt*2>TMath::Abs(lArmAlpha)
    Configurable<float> armPodCut{"armPodCut", 5.0f, "pT * (cut) > |alpha|, AP cut. Negative: no cut"};

    // --- Track quality variations (single track, both trigger and assoc daughters)
    Configurable<int> minTPCNCrossedRowsTrigger{"minTPCNCrossedRowsTrigger", 70, "Minimum TPC crossed rows (trigger)"}; // casc 80
    Configurable<int> minTPCNCrossedRowsAssociated{"minTPCNCrossedRowsAssociated", 70, "Minimum TPC crossed rows (associated)"};
    Configurable<bool> triggerRequireITS{"triggerRequireITS", true, "require ITS signal in trigger tracks"};
    Configurable<int> triggerMaxTPCSharedClusters{"triggerMaxTPCSharedClusters", 200, "maximum number of shared TPC clusters (inclusive)"};
    Configurable<bool> triggerRequireL0{"triggerRequireL0", false, "require ITS L0 cluster for trigger"};
    Configurable<float> tpcPidNsigmaCut{"tpcPidNsigmaCut", 5, "tpcPidNsigmaCut"};

    // --- Trigger: DCA variation from basic formula: |DCAxy| <  0.004f + (0.013f / pt)
    Configurable<float> dcaXYconstant{"dcaXYconstant", 0.004, "[0] in |DCAxy| < [0]+[1]/pT"};
    Configurable<float> dcaXYpTdep{"dcaXYpTdep", 0.013, "[1] in |DCAxy| < [0]+[1]/pT"};

    // --- Associated: topological variable variation (OK to vary all-at-once, at least for first study)
    Configurable<double> v0cospa{"v0cospa", 0.995, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0)
    Configurable<float> dcaV0dau{"dcaV0dau", 1.0, "DCA V0 Daughters"};
    Configurable<float> v0RadiusMax{"v0RadiusMax", 1E5, "v0radius"};
    Configurable<float> lifetimecutK0S{"lifetimecutK0S", 20, "lifetimecutK0S"};
    Configurable<float> lifetimecutLambda{"lifetimecutLambda", 30, "lifetimecutLambda"};
    Configurable<float> dcanegtopvK0S{"dcanegtopvK0S", 0.1, "DCA Neg To PV"};
    Configurable<float> dcapostopvK0S{"dcapostopvK0S", 0.1, "DCA Pos To PV"};
    Configurable<float> dcanegtopvLambda{"dcanegtopvLambda", 0.05, "DCA Neg To PV"};
    Configurable<float> dcapostopvLambda{"dcapostopvLambda", 0.2, "DCA Pos To PV"};
    Configurable<float> dcanegtopvAntiLambda{"dcanegtopvAntiLambda", 0.2, "DCA Neg To PV"};
    Configurable<float> dcapostopvAntiLambda{"dcapostopvAntiLambda", 0.05, "DCA Pos To PV"};
    Configurable<float> v0RadiusMin{"v0RadiusMin", 1.2, "v0radius"};

    // cascade selection
    Configurable<float> bachBaryonCosPA{"bachBaryonCosPA", 0.9999, "Bachelor baryon CosPA"};
    Configurable<float> bachBaryonDCAxyToPV{"bachBaryonDCAxyToPV", 0.08, "DCA bachelor baryon to PV"};
    Configurable<double> cascCospa{"cascCospa", 0.98, "cascCospa"};

    Configurable<float> dcapostopvXiminus{"dcapostopvXiminus", 0.05, "DCA of negative doughter track To PV"};
    Configurable<float> dcanegtopvXiminus{"dcanegtopvXiminus", 0.1, "DCA of positive doughter track To PV"};
    Configurable<float> dcabachToPVXiminus{"dcabachToPVXiminus", 0.07, "DCA Bach To PV"};

    Configurable<float> dcapostopvXiplus{"dcapostopvXiplus", 0.1, "DCA of negative doughter track To PV"};
    Configurable<float> dcanegtopvXiplus{"dcanegtopvXiplus", 0.05, "DCA of positive doughter track To PV"};
    Configurable<float> dcabachToPVXiplus{"dcabachToPVXiplus", 0.07, "DCA Bach To PV"};

    Configurable<float> dcapostopvOmegaMinus{"dcapostopvOmegaMinus", 0.05, "DCA of negative doughter track To PV"};
    Configurable<float> dcanegtopvOmegaMinus{"dcanegtopvOmegaMinus", 0.1, "DCA of positive doughter track To PV"};
    Configurable<float> dcabachToPVOmegaMinus{"dcabachToPVOmegaMinus", 0.07, "DCA Bach To PV"};

    Configurable<float> dcapostopvOmegaPlus{"dcapostopvOmegaPlus", 0.1, "DCA of negative doughter track To PV"};
    Configurable<float> dcanegtopvOmegaPlus{"dcanegtopvOmegaPlus", 0.05, "DCA of positive doughter track To PV"};
    Configurable<float> dcabachToPVOmegaPlus{"dcabachToPVOmegaPlus", 0.07, "DCA Bach To PV"};

    Configurable<float> cascdcaV0dau{"cascdcaV0dau", 0.5, "DCA V0 Daughters"};
    Configurable<float> cascDcacascdau{"cascDcacascdau", 0.8, "cascDcacascdau"};
    Configurable<float> dcaCacsDauPar0{"dcaCacsDauPar0", 0.8, " par for pt dep DCA cascade daughter cut, p_T < 1 GeV/c"};
    Configurable<float> dcaCacsDauPar1{"dcaCacsDauPar1", 0.5, " par for pt dep DCA cascade daughter cut, 1< p_T < 4 GeV/c"};
    Configurable<float> dcaCacsDauPar2{"dcaCacsDauPar2", 0.2, " par for pt dep DCA cascade daughter cut, p_T > 4 GeV/c"};
    Configurable<float> cascdcaV0ToPV{"cascdcaV0ToPV", 0.06, "DCA V0 To PV"};
    Configurable<float> cascRadius{"cascRadius", 1.3, "cascRadius"};
    Configurable<double> cascv0cospa{"cascv0cospa", 0.98, "V0 CosPA"};
    Configurable<float> cascv0RadiusMin{"cascv0RadiusMin", 2.5, "v0radius"};
    Configurable<float> proplifetime{"proplifetime", 3, "ctau/<ctau>"};
    Configurable<float> lambdaMassWin{"lambdaMassWin", 0.005, "V0 Mass window limit"};
    Configurable<float> rejcomp{"rejcomp", 0.008, "Competing Cascade rejection"};
    Configurable<float> rapCut{"rapCut", 0.8, "Rapidity acceptance"};

    // dE/dx for associated daughters
    Configurable<int> dEdxCompatibility{"dEdxCompatibility", 1, "0: loose, 1: normal, 2: tight. Defined in HStrangeCorrelationFilter"};

    // (N.B.: sources that can be investigated in post are not listed!)
  } systCuts;

  // objects to use for efficiency corrections
  TH2F* hEfficiencyTrigger;
  TH2F* hEfficiencyPion;
  TH2F* hEfficiencyK0Short;
  TH2F* hEfficiencyLambda;
  TH2F* hEfficiencyAntiLambda;
  TH2F* hEfficiencyXiMinus;
  TH2F* hEfficiencyXiPlus;
  TH2F* hEfficiencyOmegaMinus;
  TH2F* hEfficiencyOmegaPlus;
  TH2F* hEfficiencyHadron;

  using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
  BinningType colBinning{{axisVtxZ, axisMult}, true}; // true is for 'ignore overflows' (true by default). Underflows and overflows will have bin -1.

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
  uint64_t bitMap;
  int mRunNumber;
  int mRunNumberZorro;

  const float ctauxiPDG = 4.91;     // from PDG
  const float ctauomegaPDG = 2.461; // from PDG

  std::vector<std::vector<float>> axisRanges;

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
    hEfficiencyK0Short = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyK0Short"));
    hEfficiencyLambda = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyLambda"));
    hEfficiencyAntiLambda = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyAntiLambda"));
    hEfficiencyXiMinus = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyXiMinus"));
    hEfficiencyXiPlus = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyXiPlus"));
    hEfficiencyOmegaMinus = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyOmegaMinus"));
    hEfficiencyOmegaPlus = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyOmegaPlus"));
    hEfficiencyHadron = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyHadron"));
    hEfficiencyPion = static_cast<TH2F*>(listEfficiencies->FindObject("hEfficiencyPion"));
    LOG(info) << "Efficiencies now loaded for " << mRunNumber;
  }

  template <typename TV0>
  uint64_t V0selectionBitmap(TV0 v0, float pvx, float pvy, float pvz, bool fillHists)
  // precalculate this information so that a check is one mask operation, not many
  {
    uint64_t bitMap = 0;
    auto postrack = v0.template posTrack_as<TracksComplete>();
    auto negtrack = v0.template negTrack_as<TracksComplete>();
    // proper lifetime and TPC PID and DCA daughter to prim.vtx
    if (doCorrelationK0Short) {
      // proper lifetime
      if (fillHists)
        histos.fill(HIST("V0selection/lifetimeK0Short"), v0.distovertotmom(pvx, pvy, pvz) * o2::constants::physics::MassK0Short);
      if (v0.distovertotmom(pvx, pvy, pvz) * o2::constants::physics::MassK0Short < systCuts.lifetimecutK0S)
        SETBIT(bitMap, 0);
      // TPC PID
      if (fillHists) {
        histos.fill(HIST("V0selection/tpcNsigmaPosK0Short"), postrack.tpcNSigmaPi());
        histos.fill(HIST("V0selection/tpcNsigmaNegK0Short"), negtrack.tpcNSigmaPi());
      }
      if (std::fabs(postrack.tpcNSigmaPi()) < systCuts.tpcPidNsigmaCut && std::fabs(negtrack.tpcNSigmaPi()) < systCuts.tpcPidNsigmaCut)
        SETBIT(bitMap, 3);
      // DCA daughter to prim.vtx
      if (fillHists) {
        histos.fill(HIST("V0selection/dcapostopvK0Short"), v0.dcapostopv());
        histos.fill(HIST("V0selection/dcanegtopvK0Short"), v0.dcanegtopv());
      }
      if (std::abs(v0.dcapostopv()) > systCuts.dcapostopvK0S && std::abs(v0.dcanegtopv()) > systCuts.dcanegtopvK0S)
        SETBIT(bitMap, 6);
      // armenteros
      if (fillHists)
        histos.fill(HIST("V0selection/hArmenteros-PodolanskiK0Short"), v0.alpha(), v0.qtarm());
      if (v0.qtarm() * systCuts.armPodCut > std::abs(v0.alpha()))
        SETBIT(bitMap, 9);
    }
    if (doCorrelationLambda) {
      // proper lifetime
      if (fillHists)
        histos.fill(HIST("V0selection/lifetimeLambda"), v0.distovertotmom(pvx, pvy, pvz) * o2::constants::physics::MassLambda0);
      if (v0.distovertotmom(pvx, pvy, pvz) * o2::constants::physics::MassLambda0 < systCuts.lifetimecutLambda)
        SETBIT(bitMap, 1);
      // TPC PID
      if (fillHists) {
        histos.fill(HIST("V0selection/tpcNsigmaPosLambda"), postrack.tpcNSigmaPr());
        histos.fill(HIST("V0selection/tpcNsigmaNegLambda"), negtrack.tpcNSigmaPi());
      }
      if (std::fabs(postrack.tpcNSigmaPr()) < systCuts.tpcPidNsigmaCut && std::fabs(negtrack.tpcNSigmaPi()) < systCuts.tpcPidNsigmaCut)
        SETBIT(bitMap, 4);
      // DCA daughter to prim.vtx and armenteros
      if (fillHists) {
        histos.fill(HIST("V0selection/dcapostopvLambda"), v0.dcapostopv());
        histos.fill(HIST("V0selection/dcanegtopvLambda"), v0.dcanegtopv());
      }
      if (std::abs(v0.dcapostopv()) > systCuts.dcapostopvLambda && std::abs(v0.dcanegtopv()) > systCuts.dcanegtopvLambda)
        SETBIT(bitMap, 7);
      // armenteros
      if (fillHists)
        histos.fill(HIST("V0selection/hArmenteros-PodolanskiLambda"), v0.alpha(), v0.qtarm());
      if (v0.qtarm() * systCuts.armPodCut < std::abs(v0.alpha())) // lambda
        SETBIT(bitMap, 10);
    }
    if (doCorrelationAntiLambda) {
      // proper lifetime
      if (fillHists)
        histos.fill(HIST("V0selection/lifetimeAntiLambda"), v0.distovertotmom(pvx, pvy, pvz) * o2::constants::physics::MassLambda0);
      if (v0.distovertotmom(pvx, pvy, pvz) * o2::constants::physics::MassLambda0 < systCuts.lifetimecutLambda)
        SETBIT(bitMap, 2);
      // TPC PID
      if (fillHists) {
        histos.fill(HIST("V0selection/tpcNsigmaPosAntiLambda"), postrack.tpcNSigmaPi());
        histos.fill(HIST("V0selection/tpcNsigmaNegAntiLambda"), negtrack.tpcNSigmaPr());
      }
      if (std::fabs(postrack.tpcNSigmaPi()) < systCuts.tpcPidNsigmaCut && std::fabs(negtrack.tpcNSigmaPr()) < systCuts.tpcPidNsigmaCut)
        SETBIT(bitMap, 5);
      // DCA daughter to prim.vtx
      if (fillHists) {
        histos.fill(HIST("V0selection/dcapostopvAntiLambda"), v0.dcapostopv());
        histos.fill(HIST("V0selection/dcanegtopvAntiLambda"), v0.dcanegtopv());
      }
      if (std::abs(v0.dcapostopv()) > systCuts.dcapostopvAntiLambda && std::abs(v0.dcanegtopv()) > systCuts.dcanegtopvAntiLambda)
        SETBIT(bitMap, 8);
      // armenteros
      if (fillHists)
        histos.fill(HIST("V0selection/hArmenteros-PodolanskiAntiLambda"), v0.alpha(), v0.qtarm());
      if (v0.qtarm() * systCuts.armPodCut < std::abs(v0.alpha())) // lambda
        SETBIT(bitMap, 11);
    }
    return bitMap;
  }

  template <typename TCascade>
  uint64_t CascadeselectionBitmap(TCascade casc, float pvx, float pvy, float pvz, bool fillHists)
  {
    uint64_t bitMap = 0;
    float cascpos = std::hypot(casc.x() - pvx, casc.y() - pvy, casc.z() - pvz);
    float cascptotmom = std::hypot(casc.px(), casc.py(), casc.pz());
    float ctauXi = o2::constants::physics::MassXiMinus * cascpos / ((cascptotmom + 1e-13) * ctauxiPDG);
    float ctauOmega = o2::constants::physics::MassOmegaMinus * cascpos / ((cascptotmom + 1e-13) * ctauomegaPDG);
    auto postrack = casc.template posTrack_as<TracksComplete>();
    auto negtrack = casc.template negTrack_as<TracksComplete>();
    auto bachtrack = casc.template bachelor_as<TracksComplete>();

    // TPC PID and DCA daughter to prim.vtx and comopeting casc.rej and life time
    if (doCorrelationXiMinus) {
      // TPC PID
      if (fillHists) {
        histos.fill(HIST("Cascadeselection/tpcNSigmaPosXiMinus"), postrack.tpcNSigmaPr());
        histos.fill(HIST("Cascadeselection/tpcNSigmaNegXiMinus"), negtrack.tpcNSigmaPi());
        histos.fill(HIST("Cascadeselection/tpcNSigmaBachXiMinus"), bachtrack.tpcNSigmaPi());
      }
      if (std::fabs(postrack.tpcNSigmaPr()) < systCuts.tpcPidNsigmaCut && std::fabs(negtrack.tpcNSigmaPi()) < systCuts.tpcPidNsigmaCut &&
          std::fabs(bachtrack.tpcNSigmaPi()) < systCuts.tpcPidNsigmaCut)
        SETBIT(bitMap, 0);
      // DCA daughter to prim.vtx
      if (fillHists) {
        histos.fill(HIST("Cascadeselection/dcapostopvXiMinus"), casc.dcapostopv());
        histos.fill(HIST("Cascadeselection/dcanegtopvXiMinus"), casc.dcanegtopv());
        histos.fill(HIST("Cascadeselection/dcabachtopvXiMinus"), casc.dcabachtopv());
      }
      if (std::abs(casc.dcabachtopv()) > systCuts.dcabachToPVXiminus && std::abs(casc.dcapostopv()) > systCuts.dcapostopvXiminus &&
          std::abs(casc.dcanegtopv()) > systCuts.dcanegtopvXiminus)
        SETBIT(bitMap, 4);
      // comopeting casc.rej
      if (fillHists)
        histos.fill(HIST("Cascadeselection/rejcompXiMinus"), casc.mOmega() - o2::constants::physics::MassOmegaMinus);
      if (std::abs(casc.mOmega() - o2::constants::physics::MassOmegaMinus) > systCuts.rejcomp)
        SETBIT(bitMap, 8);
      // life time
      if (fillHists)
        histos.fill(HIST("Cascadeselection/proplifetimeXiMinus"), ctauXi);
      if (ctauXi < systCuts.proplifetime)
        SETBIT(bitMap, 12);
      // y cut
      if (fillHists)
        histos.fill(HIST("Cascadeselection/ycutXiMinus"), casc.yXi());
      if (std::abs(casc.yXi()) < systCuts.rapCut)
        SETBIT(bitMap, 16);
    }
    if (doCorrelationXiPlus) {
      // TPC PID
      if (fillHists) {
        histos.fill(HIST("Cascadeselection/tpcNSigmaPosXiPlus"), postrack.tpcNSigmaPi());
        histos.fill(HIST("Cascadeselection/tpcNSigmaNegXiPlus"), negtrack.tpcNSigmaPr());
        histos.fill(HIST("Cascadeselection/tpcNSigmaBachXiPlus"), bachtrack.tpcNSigmaPi());
      }
      if (std::fabs(postrack.tpcNSigmaPi()) < systCuts.tpcPidNsigmaCut && std::fabs(negtrack.tpcNSigmaPr()) < systCuts.tpcPidNsigmaCut &&
          std::fabs(bachtrack.tpcNSigmaPi()) < systCuts.tpcPidNsigmaCut)
        SETBIT(bitMap, 1);
      // DCA daughter to prim.vtx
      if (fillHists) {
        histos.fill(HIST("Cascadeselection/dcapostopvXiPlus"), casc.dcapostopv());
        histos.fill(HIST("Cascadeselection/dcanegtopvXiPlus"), casc.dcanegtopv());
        histos.fill(HIST("Cascadeselection/dcabachtopvXiPlus"), casc.dcabachtopv());
      }
      if (std::abs(casc.dcabachtopv()) > systCuts.dcabachToPVXiplus && std::abs(casc.dcapostopv()) > systCuts.dcapostopvXiplus &&
          std::abs(casc.dcanegtopv()) > systCuts.dcanegtopvXiplus)
        SETBIT(bitMap, 5);
      // comopeting casc.rej
      if (fillHists)
        histos.fill(HIST("Cascadeselection/rejcompXiPlus"), casc.mOmega() - o2::constants::physics::MassOmegaMinus);
      if (std::abs(casc.mOmega() - o2::constants::physics::MassOmegaMinus) > systCuts.rejcomp)
        SETBIT(bitMap, 9);
      // life time
      if (fillHists)
        histos.fill(HIST("Cascadeselection/proplifetimeXiPlus"), ctauXi);
      if (ctauXi < systCuts.proplifetime)
        SETBIT(bitMap, 13);
      // y cut
      if (fillHists)
        histos.fill(HIST("Cascadeselection/ycutXiPlus"), casc.yXi());
      if (std::abs(casc.yXi()) > systCuts.rapCut)
        SETBIT(bitMap, 17);
    }
    if (doCorrelationOmegaMinus) {
      // TPC PID
      if (fillHists) {
        histos.fill(HIST("Cascadeselection/tpcNSigmaPosOmegaMinus"), postrack.tpcNSigmaPr());
        histos.fill(HIST("Cascadeselection/tpcNSigmaNegOmegaMinus"), negtrack.tpcNSigmaPi());
        histos.fill(HIST("Cascadeselection/tpcNSigmaBachOmegaMinus"), bachtrack.tpcNSigmaKa());
      }
      if (std::fabs(postrack.tpcNSigmaPr()) < systCuts.tpcPidNsigmaCut && std::fabs(negtrack.tpcNSigmaPi()) < systCuts.tpcPidNsigmaCut &&
          std::fabs(bachtrack.tpcNSigmaKa()) < systCuts.tpcPidNsigmaCut)
        SETBIT(bitMap, 2);
      // DCA daughter to prim.vtx
      if (fillHists) {
        histos.fill(HIST("Cascadeselection/dcapostopvOmegaMinus"), casc.dcapostopv());
        histos.fill(HIST("Cascadeselection/dcanegtopvOmegaMinus"), casc.dcanegtopv());
        histos.fill(HIST("Cascadeselection/dcabachtopvOmegaMinus"), casc.dcabachtopv());
      }
      if (std::abs(casc.dcabachtopv()) > systCuts.dcabachToPVOmegaMinus && std::abs(casc.dcapostopv()) > systCuts.dcapostopvOmegaMinus &&
          std::abs(casc.dcanegtopv()) > systCuts.dcanegtopvOmegaMinus)
        SETBIT(bitMap, 6);
      // comopeting casc.rej
      if (fillHists)
        histos.fill(HIST("Cascadeselection/rejcompOmegaMinus"), casc.mXi() - o2::constants::physics::MassXiMinus);
      if (std::abs(casc.mXi() - o2::constants::physics::MassXiMinus) > systCuts.rejcomp)
        SETBIT(bitMap, 10);
      // life time
      if (fillHists)
        histos.fill(HIST("Cascadeselection/proplifetimeOmegaMinus"), ctauOmega);
      if (ctauOmega < systCuts.proplifetime)
        SETBIT(bitMap, 14);
      // y cut
      if (fillHists)
        histos.fill(HIST("Cascadeselection/ycutOmegaMinus"), casc.yOmega());
      if (std::abs(casc.yOmega()) < systCuts.rapCut)
        SETBIT(bitMap, 18);
    }
    if (doCorrelationOmegaPlus) {
      // TPC PID
      if (fillHists) {
        histos.fill(HIST("Cascadeselection/tpcNSigmaPosOmegaPlus"), postrack.tpcNSigmaPi());
        histos.fill(HIST("Cascadeselection/tpcNSigmaNegOmegaPlus"), negtrack.tpcNSigmaPr());
        histos.fill(HIST("Cascadeselection/tpcNSigmaBachOmegaPlus"), bachtrack.tpcNSigmaKa());
      }
      if (std::fabs(postrack.tpcNSigmaPi()) < systCuts.tpcPidNsigmaCut && std::fabs(negtrack.tpcNSigmaPr()) < systCuts.tpcPidNsigmaCut &&
          std::fabs(bachtrack.tpcNSigmaKa()) < systCuts.tpcPidNsigmaCut)
        SETBIT(bitMap, 3);
      // DCA daughter to prim.vtx
      if (fillHists) {
        histos.fill(HIST("Cascadeselection/dcapostopvOmegaPlus"), casc.dcapostopv());
        histos.fill(HIST("Cascadeselection/dcanegtopvOmegaPlus"), casc.dcanegtopv());
        histos.fill(HIST("Cascadeselection/dcabachtopvOmegaPlus"), casc.dcabachtopv());
      }
      if (std::abs(casc.dcabachtopv()) > systCuts.dcabachToPVOmegaPlus && std::abs(casc.dcapostopv()) > systCuts.dcapostopvOmegaPlus &&
          std::abs(casc.dcanegtopv()) > systCuts.dcanegtopvOmegaPlus)
        SETBIT(bitMap, 7);
      // comopeting casc.rej
      if (fillHists)
        histos.fill(HIST("Cascadeselection/rejcompOmegaPlus"), casc.mXi() - o2::constants::physics::MassXiMinus);
      if (std::abs(casc.mXi() - o2::constants::physics::MassXiMinus) > systCuts.rejcomp)
        SETBIT(bitMap, 11);
      // life time
      if (fillHists)
        histos.fill(HIST("Cascadeselection/proplifetimeOmegaPlus"), ctauOmega);
      if (ctauOmega < systCuts.proplifetime)
        SETBIT(bitMap, 15);
      // y cut
      if (fillHists)
        histos.fill(HIST("Cascadeselection/ycutOmegaPlus"), casc.yOmega());
      if (std::abs(casc.yOmega()) > systCuts.rapCut)
        SETBIT(bitMap, 19);
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
    if (track.pt() < axisRanges[3][0]) {
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

  // V0selection
  template <typename TV0>
  bool V0Selected(TV0 v0, bool fillHists)
  {
    auto postrack = v0.template posTrack_as<TracksComplete>();
    auto negtrack = v0.template negTrack_as<TracksComplete>();
    // v0radius
    if (fillHists)
      histos.fill(HIST("V0selection/hv0radius"), v0.v0radius());
    if (v0.v0radius() < systCuts.v0RadiusMin) {
      return false;
    }
    if (v0.v0radius() > systCuts.v0RadiusMax) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("V0selection/hv0radius_after"), v0.v0radius());
    // v0cosPA
    if (fillHists)
      histos.fill(HIST("V0selection/hv0cosPA"), v0.v0cosPA());
    if (v0.v0cosPA() < systCuts.v0cospa) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("V0selection/hv0cosPA_after"), v0.v0cosPA());
    // dcaV0daughters
    if (fillHists)
      histos.fill(HIST("V0selection/hdcaV0daughters"), v0.dcaV0daughters());
    if (v0.dcaV0daughters() > systCuts.dcaV0dau) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("V0selection/hdcaV0daughters_after"), v0.dcaV0daughters());
    // tpcNClsCrossedRows
    if (fillHists) {
      histos.fill(HIST("V0selection/htpcNClsCrossedRowspos"), postrack.tpcNClsCrossedRows());
      histos.fill(HIST("V0selection/htpcNClsCrossedRowsneg"), negtrack.tpcNClsCrossedRows());
    }
    if (postrack.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRowsAssociated || negtrack.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRowsAssociated)
      return false;
    if (fillHists) {
      histos.fill(HIST("V0selection/htpcNClsCrossedRowspos_after"), postrack.tpcNClsCrossedRows());
      histos.fill(HIST("V0selection/htpcNClsCrossedRowsneg_after"), negtrack.tpcNClsCrossedRows());
    }
    return true;
  }

  // cascadeselection
  template <typename TCascade>
  bool CascadeSelected(TCascade casc, float pvx, float pvy, float pvz, bool fillHists)
  {
    auto postrack = casc.template posTrack_as<TracksComplete>();
    auto negtrack = casc.template negTrack_as<TracksComplete>();
    auto bachtrack = casc.template bachelor_as<TracksComplete>();
    // bachBaryonCosPA
    if (fillHists)
      histos.fill(HIST("Cascadeselection/hbachBaryonCosPA"), casc.bachBaryonCosPA());
    if (casc.bachBaryonCosPA() < systCuts.bachBaryonCosPA) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("Cascadeselection/hbachBaryonCosPA_after"), casc.bachBaryonCosPA());
    // bachBaryonDCAxyToPV
    if (fillHists)
      histos.fill(HIST("Cascadeselection/hbachBaryonDCAxyToPV"), casc.bachBaryonDCAxyToPV());
    if (std::abs(casc.bachBaryonDCAxyToPV()) > systCuts.bachBaryonDCAxyToPV) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("Cascadeselection/hbachBaryonDCAxyToPV_after"), casc.bachBaryonDCAxyToPV());
    // casccosPA
    if (fillHists)
      histos.fill(HIST("Cascadeselection/hcasccosPA"), casc.casccosPA(pvx, pvy, pvz));
    if (casc.casccosPA(pvx, pvy, pvz) < systCuts.cascCospa) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("Cascadeselection/hcasccosPA_after"), casc.casccosPA(pvx, pvy, pvz));
    // dcacascdaughters
    if (fillHists)
      histos.fill(HIST("Cascadeselection/hdcacascdaughters"), casc.dcacascdaughters());
    float ptDepCut = systCuts.dcaCacsDauPar0;
    if (casc.pt() > 1 && casc.pt() < 4)
      ptDepCut = systCuts.dcaCacsDauPar1;
    else if (casc.pt() > 4)
      ptDepCut = systCuts.dcaCacsDauPar2;
    if (casc.dcacascdaughters() > ptDepCut)
      return false;
    if (fillHists)
      histos.fill(HIST("Cascadeselection/hdcacascdaughters_after"), casc.dcacascdaughters());
    // dcaV0daughters
    if (fillHists)
      histos.fill(HIST("Cascadeselection/hdcaV0daughters"), casc.dcaV0daughters());
    if (casc.dcaV0daughters() > systCuts.dcaV0dau)
      return false;
    if (fillHists)
      histos.fill(HIST("Cascadeselection/hdcaV0daughters_after"), casc.dcaV0daughters());
    // dcav0topv
    if (fillHists)
      histos.fill(HIST("Cascadeselection/hdcav0topv"), casc.dcav0topv(pvx, pvy, pvz));
    if (std::abs(casc.dcav0topv(pvx, pvy, pvz)) < systCuts.cascdcaV0ToPV)
      return false;
    if (fillHists)
      histos.fill(HIST("Cascadeselection/hdcav0topv_after"), casc.dcav0topv(pvx, pvy, pvz));
    // cascradius
    if (fillHists)
      histos.fill(HIST("Cascadeselection/hcascradius"), casc.cascradius());
    if (casc.cascradius() < systCuts.cascRadius)
      return false;
    if (fillHists)
      histos.fill(HIST("Cascadeselection/hcascradius_after"), casc.cascradius());
    // v0radius
    if (fillHists)
      histos.fill(HIST("Cascadeselection/hv0radius"), casc.v0radius());
    if (casc.v0radius() < systCuts.cascv0RadiusMin)
      return false;
    if (fillHists)
      histos.fill(HIST("Cascadeselection/hv0radius_after"), casc.v0radius());
    // v0cosPA
    if (fillHists)
      histos.fill(HIST("Cascadeselection/hv0cosPA"), casc.v0cosPA(casc.x(), casc.y(), casc.z()));
    if (casc.v0cosPA(casc.x(), casc.y(), casc.z()) < systCuts.cascv0cospa)
      return false;
    if (fillHists)
      histos.fill(HIST("Cascadeselection/hv0cosPA_after"), casc.v0cosPA(casc.x(), casc.y(), casc.z()));
    // lambdaMassWin
    if (fillHists)
      histos.fill(HIST("Cascadeselection/hmLambdadiff"), casc.mLambda() - o2::constants::physics::MassLambda0);
    if (std::abs(casc.mLambda() - o2::constants::physics::MassLambda0) > systCuts.lambdaMassWin)
      return false;
    if (fillHists)
      histos.fill(HIST("Cascadeselection/hmLambdadiff_after"), casc.mLambda() - o2::constants::physics::MassLambda0);
    //---] track quality check [---
    if (fillHists) {
      histos.fill(HIST("Cascadeselection/htpcNClsCrossedRowspos"), postrack.tpcNClsCrossedRows());
      histos.fill(HIST("Cascadeselection/htpcNClsCrossedRowsneg"), negtrack.tpcNClsCrossedRows());
      histos.fill(HIST("Cascadeselection/htpcNClsCrossedRowsbach"), bachtrack.tpcNClsCrossedRows());
    }
    if (postrack.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRowsAssociated || negtrack.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRowsAssociated || bachtrack.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRowsAssociated)
      return false;
    if (fillHists) {
      histos.fill(HIST("Cascadeselection/htpcNClsCrossedRowspos_after"), postrack.tpcNClsCrossedRows());
      histos.fill(HIST("Cascadeselection/htpcNClsCrossedRowsneg_after"), negtrack.tpcNClsCrossedRows());
      histos.fill(HIST("Cascadeselection/htpcNClsCrossedRowsbach_after"), bachtrack.tpcNClsCrossedRows());
    }
    return true;
  }

  void fillCorrelationsV0(aod::TriggerTracks const& triggers, aod::AssocV0s const& assocs, bool mixing, float pvx, float pvy, float pvz, float mult)
  {
    for (auto const& triggerTrack : triggers) {
      if (doTriggPhysicalPrimary && !triggerTrack.mcPhysicalPrimary())
        continue;
      auto trigg = triggerTrack.track_as<TracksComplete>();
      if (!isValidTrigger(trigg))
        continue;

      if (!mixing) {
        float efficiency = 1.0f;
        if (applyEfficiencyForTrigger) {
          efficiency = hEfficiencyTrigger->Interpolate(trigg.pt(), trigg.eta());
        }
        float weight = (applyEfficiencyForTrigger) ? 1. / efficiency : 1.0f;
        histos.fill(HIST("sameEvent/TriggerParticlesV0"), trigg.pt(), mult, weight);
      }

      for (auto const& assocCandidate : assocs) {
        auto assoc = assocCandidate.v0Core_as<V0DatasWithoutTrackX>();
        //---] syst cuts [---
        if (!V0Selected(assoc, false))
          continue;
        uint64_t selMap = V0selectionBitmap(assoc, pvx, pvy, pvz, false);

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
        static_for<0, 2>([&](auto i) {
          constexpr int Index = i.value;
          float efficiency = 1.0f;
          if (applyEfficiencyCorrection) {
            efficiency = hEfficiencyV0[Index]->Interpolate(ptassoc, assoc.eta());
          }
          if (applyEfficiencyForTrigger) {
            efficiency = efficiency * hEfficiencyTrigger->Interpolate(pttrigger, trigg.eta());
          }

          float weight = (applyEfficiencyCorrection || applyEfficiencyForTrigger) ? 1. / efficiency : 1.0f;
          if (TESTBIT(doCorrelation, Index) && TESTBIT(selMap, Index) && TESTBIT(selMap, Index + 3) && TESTBIT(selMap, Index + 6) && TESTBIT(selMap, Index + 9) && (!applyEfficiencyCorrection || efficiency != 0)) {
            if (assocCandidate.compatible(Index, systCuts.dEdxCompatibility) && (!doMCassociation || assocCandidate.mcTrue(Index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && !mixing && -massWindowConfigurations.maxBgNSigma < assocCandidate.invMassNSigma(Index) && assocCandidate.invMassNSigma(Index) < -massWindowConfigurations.minBgNSigma) {
              if (doCentSelectedhigh && mult < centhighMax && mult > centhighMin)
                histos.fill(HIST("sameEvent/LeftBg/CentHigh/") + HIST(kV0names[Index]), deltaphi, deltaeta, ptassoc, pvz, mult, weight);
              if (doCentSelectedmed && mult < centmedMax && mult > centmedMin)
                histos.fill(HIST("sameEvent/LeftBg/CentMed/") + HIST(kV0names[Index]), deltaphi, deltaeta, ptassoc, pvz, mult, weight);
              if (doCentSelectedlow && mult < centlowMax && mult > centlowMin)
                histos.fill(HIST("sameEvent/LeftBg/CentLow/") + HIST(kV0names[Index]), deltaphi, deltaeta, ptassoc, pvz, mult, weight);
            }
            if (assocCandidate.compatible(Index, systCuts.dEdxCompatibility) && (!doMCassociation || assocCandidate.mcTrue(Index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && !mixing && -massWindowConfigurations.maxPeakNSigma < assocCandidate.invMassNSigma(Index) && assocCandidate.invMassNSigma(Index) < +massWindowConfigurations.maxPeakNSigma) {
              if (doCentSelectedhigh && mult < centhighMax && mult > centhighMin)
                histos.fill(HIST("sameEvent/Signal/CentHigh/") + HIST(kV0names[Index]), deltaphi, deltaeta, ptassoc, pvz, mult, weight);
              if (doCentSelectedmed && mult < centmedMax && mult > centmedMin)
                histos.fill(HIST("sameEvent/Signal/CentMed/") + HIST(kV0names[Index]), deltaphi, deltaeta, ptassoc, pvz, mult, weight);
              if (doCentSelectedlow && mult < centlowMax && mult > centlowMin)
                histos.fill(HIST("sameEvent/Signal/CentLow/") + HIST(kV0names[Index]), deltaphi, deltaeta, ptassoc, pvz, mult, weight);
            }
            if (assocCandidate.compatible(Index, systCuts.dEdxCompatibility) && (!doMCassociation || assocCandidate.mcTrue(Index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && !mixing && +massWindowConfigurations.minBgNSigma < assocCandidate.invMassNSigma(Index) && assocCandidate.invMassNSigma(Index) < +massWindowConfigurations.maxBgNSigma) {
              if (doCentSelectedhigh && mult < centhighMax && mult > centhighMin)
                histos.fill(HIST("sameEvent/RightBg/CentHigh/") + HIST(kV0names[Index]), deltaphi, deltaeta, ptassoc, pvz, mult, weight);
              if (doCentSelectedmed && mult < centmedMax && mult > centmedMin)
                histos.fill(HIST("sameEvent/RightBg/CentMed/") + HIST(kV0names[Index]), deltaphi, deltaeta, ptassoc, pvz, mult, weight);
              if (doCentSelectedlow && mult < centlowMax && mult > centlowMin)
                histos.fill(HIST("sameEvent/RightBg/CentLow/") + HIST(kV0names[Index]), deltaphi, deltaeta, ptassoc, pvz, mult, weight);
            }

            if (assocCandidate.compatible(Index, systCuts.dEdxCompatibility) && (!doMCassociation || assocCandidate.mcTrue(Index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && mixing && -massWindowConfigurations.maxBgNSigma < assocCandidate.invMassNSigma(Index) && assocCandidate.invMassNSigma(Index) < -massWindowConfigurations.minBgNSigma) {
              if (doCentSelectedhigh && mult < centhighMax && mult > centhighMin)
                histos.fill(HIST("mixedEvent/LeftBg/CentHigh/") + HIST(kV0names[Index]), deltaphi, deltaeta, ptassoc, pvz, mult, weight);
              if (doCentSelectedmed && mult < centmedMax && mult > centmedMin)
                histos.fill(HIST("mixedEvent/LeftBg/CentMed/") + HIST(kV0names[Index]), deltaphi, deltaeta, ptassoc, pvz, mult, weight);
              if (doCentSelectedlow && mult < centlowMax && mult > centlowMin)
                histos.fill(HIST("mixedEvent/LeftBg/CentLow/") + HIST(kV0names[Index]), deltaphi, deltaeta, ptassoc, pvz, mult, weight);
            }
            if (assocCandidate.compatible(Index, systCuts.dEdxCompatibility) && (!doMCassociation || assocCandidate.mcTrue(Index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && mixing && -massWindowConfigurations.maxPeakNSigma < assocCandidate.invMassNSigma(Index) && assocCandidate.invMassNSigma(Index) < +massWindowConfigurations.maxPeakNSigma) {
              if (doCentSelectedhigh && mult < centhighMax && mult > centhighMin)
                histos.fill(HIST("mixedEvent/Signal/CentHigh/") + HIST(kV0names[Index]), deltaphi, deltaeta, ptassoc, pvz, mult, weight);
              if (doCentSelectedmed && mult < centmedMax && mult > centmedMin)
                histos.fill(HIST("mixedEvent/Signal/CentMed/") + HIST(kV0names[Index]), deltaphi, deltaeta, ptassoc, pvz, mult, weight);
              if (doCentSelectedlow && mult < centlowMax && mult > centlowMin)
                histos.fill(HIST("mixedEvent/Signal/CentLow/") + HIST(kV0names[Index]), deltaphi, deltaeta, ptassoc, pvz, mult, weight);
            }
            if (assocCandidate.compatible(Index, systCuts.dEdxCompatibility) && (!doMCassociation || assocCandidate.mcTrue(Index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && mixing && +massWindowConfigurations.minBgNSigma < assocCandidate.invMassNSigma(Index) && assocCandidate.invMassNSigma(Index) < +massWindowConfigurations.maxBgNSigma) {
              if (doCentSelectedhigh && mult < centhighMax && mult > centhighMin)
                histos.fill(HIST("mixedEvent/RightBg/CentHigh/") + HIST(kV0names[Index]), deltaphi, deltaeta, ptassoc, pvz, mult, weight);
              if (doCentSelectedmed && mult < centmedMax && mult > centmedMin)
                histos.fill(HIST("mixedEvent/RightBg/CentMed/") + HIST(kV0names[Index]), deltaphi, deltaeta, ptassoc, pvz, mult, weight);
              if (doCentSelectedlow && mult < centlowMax && mult > centlowMin)
                histos.fill(HIST("mixedEvent/RightBg/CentLow/") + HIST(kV0names[Index]), deltaphi, deltaeta, ptassoc, pvz, mult, weight);
            }
          }
        });
      }
    }
  }

  void fillCorrelationsCascade(aod::TriggerTracks const& triggers, aod::AssocCascades const& assocs, bool mixing, float pvx, float pvy, float pvz, float mult)
  {
    for (auto const& triggerTrack : triggers) {
      if (doTriggPhysicalPrimary && !triggerTrack.mcPhysicalPrimary())
        continue;
      auto trigg = triggerTrack.track_as<TracksComplete>();
      if (!isValidTrigger(trigg))
        continue;

      if (!mixing) {
        float efficiency = 1.0f;
        if (applyEfficiencyForTrigger) {
          efficiency = hEfficiencyTrigger->Interpolate(trigg.pt(), trigg.eta());
        }
        float weight = (applyEfficiencyForTrigger) ? 1. / efficiency : 1.0f;

        histos.fill(HIST("sameEvent/TriggerParticlesCascade"), trigg.pt(), mult, weight);
      }
      for (auto const& assocCandidate : assocs) {
        auto assoc = assocCandidate.cascData();

        //---] syst cuts [---
        if (!CascadeSelected(assoc, pvx, pvy, pvz, false))
          continue;
        uint64_t CascselMap = CascadeselectionBitmap(assoc, pvx, pvy, pvz, false);

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

        static_for<0, 3>([&](auto i) {
          constexpr int Index = i.value;
          float efficiency = 1.0f;
          if (applyEfficiencyCorrection) {
            efficiency = hEfficiencyCascade[Index]->Interpolate(ptassoc, assoc.eta());
          }
          if (applyEfficiencyForTrigger) {
            efficiency = efficiency * hEfficiencyTrigger->Interpolate(pttrigger, trigg.eta());
          }
          float weight = (applyEfficiencyCorrection || applyEfficiencyForTrigger) ? 1. / efficiency : 1.0f;
          if (TESTBIT(doCorrelation, Index + 3) && TESTBIT(CascselMap, Index) && TESTBIT(CascselMap, Index + 4) && TESTBIT(CascselMap, Index + 8) && TESTBIT(CascselMap, Index + 12) && TESTBIT(CascselMap, Index + 16) && (!applyEfficiencyCorrection || efficiency != 0)) {
            if (assocCandidate.compatible(Index, systCuts.dEdxCompatibility) && (!doMCassociation || assocCandidate.mcTrue(Index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && !mixing && -massWindowConfigurations.maxBgNSigma < assocCandidate.invMassNSigma(Index) && assocCandidate.invMassNSigma(Index) < -massWindowConfigurations.minBgNSigma) {
              if (doCentSelectedhigh && mult < centhighMax && mult > centhighMin)
                histos.fill(HIST("sameEvent/LeftBg/CentHigh/") + HIST(kCascadenames[Index]), deltaphi, deltaeta, ptassoc, pvz, mult, weight);
              if (doCentSelectedmed && mult < centmedMax && mult > centmedMin)
                histos.fill(HIST("sameEvent/LeftBg/CentMed/") + HIST(kCascadenames[Index]), deltaphi, deltaeta, ptassoc, pvz, mult, weight);
              if (doCentSelectedlow && mult < centlowMax && mult > centlowMin)
                histos.fill(HIST("sameEvent/LeftBg/CentLow/") + HIST(kCascadenames[Index]), deltaphi, deltaeta, ptassoc, pvz, mult, weight);
            }
            if (assocCandidate.compatible(Index, systCuts.dEdxCompatibility) && (!doMCassociation || assocCandidate.mcTrue(Index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && !mixing && -massWindowConfigurations.maxPeakNSigma < assocCandidate.invMassNSigma(Index) && assocCandidate.invMassNSigma(Index) < +massWindowConfigurations.maxPeakNSigma) {
              if (doCentSelectedhigh && mult < centhighMax && mult > centhighMin)
                histos.fill(HIST("sameEvent/Signal/CentHigh/") + HIST(kCascadenames[Index]), deltaphi, deltaeta, ptassoc, pvz, mult, weight);
              if (doCentSelectedmed && mult < centmedMax && mult > centmedMin)
                histos.fill(HIST("sameEvent/Signal/CentMed/") + HIST(kCascadenames[Index]), deltaphi, deltaeta, ptassoc, pvz, mult, weight);
              if (doCentSelectedlow && mult < centlowMax && mult > centlowMin)
                histos.fill(HIST("sameEvent/Signal/CentLow/") + HIST(kCascadenames[Index]), deltaphi, deltaeta, ptassoc, pvz, mult, weight);
            }
            if (assocCandidate.compatible(Index, systCuts.dEdxCompatibility) && (!doMCassociation || assocCandidate.mcTrue(Index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && !mixing && +massWindowConfigurations.minBgNSigma < assocCandidate.invMassNSigma(Index) && assocCandidate.invMassNSigma(Index) < +massWindowConfigurations.maxBgNSigma) {
              if (doCentSelectedhigh && mult < centhighMax && mult > centhighMin)
                histos.fill(HIST("sameEvent/RightBg/CentHigh/") + HIST(kCascadenames[Index]), deltaphi, deltaeta, ptassoc, pvz, mult, weight);
              if (doCentSelectedmed && mult < centmedMax && mult > centmedMin)
                histos.fill(HIST("sameEvent/RightBg/CentMed/") + HIST(kCascadenames[Index]), deltaphi, deltaeta, ptassoc, pvz, mult, weight);
              if (doCentSelectedlow && mult < centlowMax && mult > centlowMin)
                histos.fill(HIST("sameEvent/RightBg/CentLow/") + HIST(kCascadenames[Index]), deltaphi, deltaeta, ptassoc, pvz, mult, weight);
            }

            if (assocCandidate.compatible(Index, systCuts.dEdxCompatibility) && (!doMCassociation || assocCandidate.mcTrue(Index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && mixing && -massWindowConfigurations.maxBgNSigma < assocCandidate.invMassNSigma(Index) && assocCandidate.invMassNSigma(Index) < -massWindowConfigurations.minBgNSigma) {
              if (doCentSelectedhigh && mult < centhighMax && mult > centhighMin)
                histos.fill(HIST("mixedEvent/LeftBg/CentHigh/") + HIST(kCascadenames[Index]), deltaphi, deltaeta, ptassoc, pvz, mult, weight);
              if (doCentSelectedmed && mult < centmedMax && mult > centmedMin)
                histos.fill(HIST("mixedEvent/LeftBg/CentMed/") + HIST(kCascadenames[Index]), deltaphi, deltaeta, ptassoc, pvz, mult, weight);
              if (doCentSelectedlow && mult < centlowMax && mult > centlowMin)
                histos.fill(HIST("mixedEvent/LeftBg/CentLow/") + HIST(kCascadenames[Index]), deltaphi, deltaeta, ptassoc, pvz, mult, weight);
            }
            if (assocCandidate.compatible(Index, systCuts.dEdxCompatibility) && (!doMCassociation || assocCandidate.mcTrue(Index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && mixing && -massWindowConfigurations.maxPeakNSigma < assocCandidate.invMassNSigma(Index) && assocCandidate.invMassNSigma(Index) < +massWindowConfigurations.maxPeakNSigma) {
              if (doCentSelectedhigh && mult < centhighMax && mult > centhighMin)
                histos.fill(HIST("mixedEvent/Signal/CentHigh/") + HIST(kCascadenames[Index]), deltaphi, deltaeta, ptassoc, pvz, mult, weight);
              if (doCentSelectedmed && mult < centmedMax && mult > centmedMin)
                histos.fill(HIST("mixedEvent/Signal/CentMed/") + HIST(kCascadenames[Index]), deltaphi, deltaeta, ptassoc, pvz, mult, weight);
              if (doCentSelectedlow && mult < centlowMax && mult > centlowMin)
                histos.fill(HIST("mixedEvent/Signal/CentLow/") + HIST(kCascadenames[Index]), deltaphi, deltaeta, ptassoc, pvz, mult, weight);
            }
            if (assocCandidate.compatible(Index, systCuts.dEdxCompatibility) && (!doMCassociation || assocCandidate.mcTrue(Index)) && (!doAssocPhysicalPrimary || assocCandidate.mcPhysicalPrimary()) && mixing && +massWindowConfigurations.minBgNSigma < assocCandidate.invMassNSigma(Index) && assocCandidate.invMassNSigma(Index) < +massWindowConfigurations.maxBgNSigma) {
              if (doCentSelectedhigh && mult < centhighMax && mult > centhighMin)
                histos.fill(HIST("mixedEvent/RightBg/CentHigh/") + HIST(kCascadenames[Index]), deltaphi, deltaeta, ptassoc, pvz, mult, weight);
              if (doCentSelectedmed && mult < centmedMax && mult > centmedMin)
                histos.fill(HIST("mixedEvent/RightBg/CentMed/") + HIST(kCascadenames[Index]), deltaphi, deltaeta, ptassoc, pvz, mult, weight);
              if (doCentSelectedlow && mult < centlowMax && mult > centlowMin)
                histos.fill(HIST("mixedEvent/RightBg/CentLow/") + HIST(kCascadenames[Index]), deltaphi, deltaeta, ptassoc, pvz, mult, weight);
            }
          }
        });
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

    hEfficiencyHadron = 0x0;

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

    // Store axis ranges to prevent spurious filling
    // axis status:
    // --- Delta-phi is safe -> math forbids insanity
    // --- Delta-eta depends on pre-filter -> check
    // --- pT assoc depends on binning -> check
    // --- vertex Z is safe -> skipped at evsel level
    // --- multiplicity -> check

    // grab axis edge from ConfigurableAxes
    const AxisSpec preAxisDeltaPhi{axisDeltaPhi, "#Delta#varphi"};
    const AxisSpec preAxisDeltaEta{axisDeltaEta, "#Delta#eta"};
    const AxisSpec preAxisPtAssoc{axisPtAssoc, "#it{p}_{T}^{assoc} (GeV/c)"};
    const AxisSpec preAxisPtTrigger{axisPtTrigger, "#it{p}_{T}^{trigger} (GeV/c)"};
    const AxisSpec preAxisVtxZ{axisVtxZ, "vertex Z (cm)"};
    const AxisSpec preAxisMult{axisMult, "mult percentile"};

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

    // event selection
    histos.add("hEventSelection", "hEventSelection", kTH1F, {{9, 0, 9}});
    TString eventSelLabel[] = {"all", "sel8", "PV_{z}", "kIsGoodITSLayersAll", "kIsGoodZvtxFT0vsPV", "OccupCut", "kNoITSROFrameBorder", "kNoSameBunchPileup ", " kNoCollInTimeRangeStandard"};
    for (int i = 1; i <= histos.get<TH1>(HIST("hEventSelection"))->GetNbinsX(); i++) {
      histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(i, eventSelLabel[i - 1]);
    }

    // V0selection cut
    histos.add("V0selection/hv0radius", ";v0radius", kTH1F, {{100, 0, 20}});
    histos.add("V0selection/hv0radius_after", ";v0radius", kTH1F, {{100, 0, 20}});
    histos.add("V0selection/hv0cosPA", ";v0cosPA", kTH1F, {{100, 0.99, 1.0}});
    histos.add("V0selection/hv0cosPA_after", ";v0cosPA", kTH1F, {{100, 0.99, 1.0}});
    histos.add("V0selection/hdcaV0daughters", ";dcaV0daughters", kTH1F, {{100, 0., 2.0}});
    histos.add("V0selection/hdcaV0daughters_after", ";dcaV0daughters", kTH1F, {{100, 0., 2.0}});
    histos.add("V0selection/htpcNClsCrossedRowspos", ";tpcNClsCrossedRowspos", kTH1F, {{200, 0, 200}});
    histos.add("V0selection/htpcNClsCrossedRowspos_after", ";tpcNClsCrossedRowspos", kTH1F, {{200, 0, 200}});
    histos.add("V0selection/htpcNClsCrossedRowsneg", ";tpcNClsCrossedRowsneg", kTH1F, {{200, 0, 200}});
    histos.add("V0selection/htpcNClsCrossedRowsneg_after", ";tpcNClsCrossedRowsneg", kTH1F, {{200, 0, 200}});
    for (int i = 0; i < 3; i++) {
      if (TESTBIT(doCorrelation, i)) {
        histos.add(fmt::format("V0selection/lifetime{}", kParticlenames[i]).c_str(), ";lifetime", kTH1F, {{50, 0, 50}});
        histos.add(fmt::format("V0selection/lifetime{}_after", kParticlenames[i]).c_str(), ";lifetime", kTH1F, {{50, 0, 50}});
        histos.add(fmt::format("V0selection/tpcNsigmaPos{}", kParticlenames[i]).c_str(), ";tpcNsigmaPos", kTH1F, {{100, -10, 10}});
        histos.add(fmt::format("V0selection/tpcNsigmaPos{}_after", kParticlenames[i]).c_str(), ";tpcNsigmaPos", kTH1F, {{100, -10, 10}});
        histos.add(fmt::format("V0selection/tpcNsigmaNeg{}", kParticlenames[i]).c_str(), ";tpcNsigmaNeg", kTH1F, {{100, -10, 10}});
        histos.add(fmt::format("V0selection/tpcNsigmaNeg{}_after", kParticlenames[i]).c_str(), ";tpcNsigmaNeg", kTH1F, {{100, -10, 10}});
        histos.add(fmt::format("V0selection/dcapostopv{}", kParticlenames[i]).c_str(), ";dcapostopv", kTH1F, {{100, -2, 2}});
        histos.add(fmt::format("V0selection/dcapostopv{}_after", kParticlenames[i]).c_str(), ";dcapostopv", kTH1F, {{100, -2, 2}});
        histos.add(fmt::format("V0selection/dcanegtopv{}", kParticlenames[i]).c_str(), ";dcanegtopv", kTH1F, {{100, -2, 2}});
        histos.add(fmt::format("V0selection/dcanegtopv{}_after", kParticlenames[i]).c_str(), ";dcanegtopv", kTH1F, {{100, -2, 2}});
        histos.add(fmt::format("V0selection/hArmenteros-Podolanski{}", kParticlenames[i]).c_str(), ";#alpha;p_{T}", kTH2F, {axisAlpha, axisptarm});
        histos.add(fmt::format("V0selection/hArmenteros-Podolanski{}_after", kParticlenames[i]).c_str(), ";#alpha;p_{T}", kTH2F, {axisAlpha, axisptarm});
      }
    }

    // cascadeselection cut
    histos.add("Cascadeselection/hbachBaryonCosPA", ";bachBaryonCosPA", kTH1F, {{1000, 0.98, 1.02}});
    histos.add("Cascadeselection/hbachBaryonCosPA_after", ";bachBaryonCosPA", kTH1F, {{100, 0.98, 1.02}});
    histos.add("Cascadeselection/hbachBaryonDCAxyToPV", ";bachBaryonDCAxyToPV", kTH1F, {{100, -0.1, 0.1}});
    histos.add("Cascadeselection/hbachBaryonDCAxyToPV_after", ";bachBaryonDCAxyToPV", kTH1F, {{100, -0.1, 0.1}});
    histos.add("Cascadeselection/hcasccosPA", ";casccosPA", kTH1F, {{200, 0.95, 1.0}});
    histos.add("Cascadeselection/hcasccosPA_after", ";casccosPA", kTH1F, {{200, 0.95, 1.0}});
    histos.add("Cascadeselection/hdcacascdaughters", ";dcacascdaughters", kTH1F, {{100, 0., 2.0}});
    histos.add("Cascadeselection/hdcacascdaughters_after", ";dcacascdaughters", kTH1F, {{100, 0., 2.0}});
    histos.add("Cascadeselection/hdcaV0daughters", ";dcaV0daughters", kTH1F, {{100, 0., 2.0}});
    histos.add("Cascadeselection/hdcaV0daughters_after", ";dcaV0daughters", kTH1F, {{100, 0., 2.0}});
    histos.add("Cascadeselection/hdcav0topv", ";dcav0topv", kTH1F, {{100, 0., 2.0}});
    histos.add("Cascadeselection/hdcav0topv_after", ";dcav0topv", kTH1F, {{100, 0., 2.0}});
    histos.add("Cascadeselection/hcascradius", ";cascradius", kTH1F, {{100, 0, 20}});
    histos.add("Cascadeselection/hcascradius_after", ";cascradius", kTH1F, {{100, 0, 20}});
    histos.add("Cascadeselection/hv0radius", ";v0radius", kTH1F, {{100, 0, 20}});
    histos.add("Cascadeselection/hv0radius_after", ";v0radius", kTH1F, {{100, 0, 20}});
    histos.add("Cascadeselection/hv0cosPA", ";v0cosPA", kTH1F, {{100, 0.98, 1.0}});
    histos.add("Cascadeselection/hv0cosPA_after", ";v0cosPA", kTH1F, {{100, 0.98, 1.0}});
    histos.add("Cascadeselection/hmLambdadiff", ";mLambdadiff", kTH1F, {{40, -0.2, 0.2}});
    histos.add("Cascadeselection/hmLambdadiff_after", ";mLambdadiff", kTH1F, {{40, -0.2, 0.2}});
    histos.add("Cascadeselection/htpcNClsCrossedRowspos", ";tpcNClsCrossedRowspos", kTH1F, {{200, 0, 200}});
    histos.add("Cascadeselection/htpcNClsCrossedRowspos_after", ";tpcNClsCrossedRowspos", kTH1F, {{200, 0, 200}});
    histos.add("Cascadeselection/htpcNClsCrossedRowsneg", ";tpcNClsCrossedRowsneg", kTH1F, {{200, 0, 200}});
    histos.add("Cascadeselection/htpcNClsCrossedRowsneg_after", ";htpcNClsCrossedRowsneg", kTH1F, {{200, 0, 200}});
    histos.add("Cascadeselection/htpcNClsCrossedRowsbach", ";tpcNClsCrossedRowsbach", kTH1F, {{200, 0, 200}});
    histos.add("Cascadeselection/htpcNClsCrossedRowsbach_after", ";tpcNClsCrossedRowsbach", kTH1F, {{200, 0, 200}});
    for (int i = 3; i < 7; i++) {
      if (TESTBIT(doCorrelation, i)) {
        histos.add(fmt::format("Cascadeselection/tpcNSigmaPos{}", kParticlenames[i]).c_str(), ";tpcNSigmaPos", kTH1F, {{100, -10, 10}});
        histos.add(fmt::format("Cascadeselection/tpcNSigmaPos{}_after", kParticlenames[i]).c_str(), ";tpcNSigmaPos", kTH1F, {{100, -10, 10}});
        histos.add(fmt::format("Cascadeselection/tpcNSigmaNeg{}", kParticlenames[i]).c_str(), ";tpcNSigmaNeg", kTH1F, {{100, -10, 10}});
        histos.add(fmt::format("Cascadeselection/tpcNSigmaNeg{}_after", kParticlenames[i]).c_str(), ";tpcNSigmaNeg", kTH1F, {{100, -10, 10}});
        histos.add(fmt::format("Cascadeselection/tpcNSigmaBach{}", kParticlenames[i]).c_str(), ";tpcNSigmaBach", kTH1F, {{100, -10, 10}});
        histos.add(fmt::format("Cascadeselection/tpcNSigmaBach{}_after", kParticlenames[i]).c_str(), ";tpcNSigmaBach", kTH1F, {{100, -10, 10}});
        histos.add(fmt::format("Cascadeselection/dcapostopv{}", kParticlenames[i]).c_str(), ";dcapostopv", kTH1F, {{100, -4, 4}});
        histos.add(fmt::format("Cascadeselection/dcapostopv{}_after", kParticlenames[i]).c_str(), ";dcapostopv", kTH1F, {{100, -4, 4}});
        histos.add(fmt::format("Cascadeselection/dcanegtopv{}", kParticlenames[i]).c_str(), ";dcanegtopv", kTH1F, {{100, -4, 4}});
        histos.add(fmt::format("Cascadeselection/dcanegtopv{}_after", kParticlenames[i]).c_str(), ";dcanegtopv", kTH1F, {{100, -4, 4}});
        histos.add(fmt::format("Cascadeselection/dcabachtopv{}", kParticlenames[i]).c_str(), ";dcabachtopv", kTH1F, {{100, -4, 4}});
        histos.add(fmt::format("Cascadeselection/dcabachtopv{}_after", kParticlenames[i]).c_str(), ";dcabachtopv", kTH1F, {{100, -4, 4}});
        histos.add(fmt::format("Cascadeselection/rejcomp{}", kParticlenames[i]).c_str(), ";rejcomp", kTH1F, {{20, -0.05, 0.05}});
        histos.add(fmt::format("Cascadeselection/rejcomp{}_after", kParticlenames[i]).c_str(), ";rejcomp", kTH1F, {{20, -0.05, 0.05}});
        histos.add(fmt::format("Cascadeselection/proplifetime{}", kParticlenames[i]).c_str(), ";proplifetime", kTH1F, {{20, 0, 8}});
        histos.add(fmt::format("Cascadeselection/proplifetime{}_after", kParticlenames[i]).c_str(), ";proplifetime", kTH1F, {{20, 0, 8}});
        histos.add(fmt::format("Cascadeselection/ycut{}", kParticlenames[i]).c_str(), ";ycut", kTH1F, {{20, -1.0, 1.0}});
        histos.add(fmt::format("Cascadeselection/ycut{}_after", kParticlenames[i]).c_str(), ";ycutXi", kTH1F, {{20, -1.0, 1.0}});
      }
    }

    // multNTracksPV vs centrality
    histos.add("hCentralityVsmultNTracksPV", "hCentralityVsmultNTracksPV", kTH2F, {axisMult, axisMultCount});
    histos.add("hCentralityVsmultNTracks", "hCentralityVsmultNTracks", kTH2F, {axisMult, axisMultCount});

    // Some QA plots
    histos.add("hGeneratedQAPtTrigger", "hGeneratedQAPtTrigger", kTH2F, {axisPtQA, {5, -0.5f, 4.5f}});
    histos.add("hGeneratedQAPtAssociatedK0", "hGeneratedQAPtAssociatedK0", kTH2F, {axisPtQA, {5, -0.5f, 4.5f}});
    histos.add("hClosureQAPtTrigger", "hClosureQAPtTrigger", kTH2F, {axisPtQA, {5, -0.5f, 4.5f}});
    histos.add("hClosureQAPtAssociatedK0", "hClosureQAPtAssociatedK0", kTH2F, {axisPtQA, {5, -0.5f, 4.5f}});

    histos.add("hTriggerAllSelectedEtaVsPt", "hTriggerAllSelectedEtaVsPt", kTH3F, {axisPtQA, axisEta, axisMult});

    histos.add("hClosureTestEventCounter", "hClosureTestEventCounter", kTH1F, {{10, 0, 10}});

    histos.add("hNumberOfRejectedPairsHadron", "hNumberOfRejectedPairsHadron", kTH1F, {{1, 0, 1}});
    histos.add("hNumberOfRejectedPairsV0", "hNumberOfRejectedPairsV0", kTH1F, {{1, 0, 1}});
    histos.add("hNumberOfRejectedPairsCascades", "hNumberOfRejectedPairsCascades", kTH1F, {{1, 0, 1}});
    histos.add("hNumberOfRejectedPairsPion", "hNumberOfRejectedPairsPion", kTH1F, {{1, 0, 1}});

    // mixing QA
    histos.add("MixingQA/hSECollisionBins", ";bin;Entries", kTH1F, {{140, -0.5, 139.5}});
    histos.add("MixingQA/hMECollisionBins", ";bin;Entries", kTH1F, {{140, -0.5, 139.5}});
    histos.add("MixingQA/hMEpvz1", ";pvz;Entries", kTH1F, {{30, -15, 15}});
    histos.add("MixingQA/hMEpvz2", ";pvz;Entries", kTH1F, {{30, -15, 15}});

    // Event QA
    histos.add("EventQA/hMixingQA", "mixing QA", kTH1F, {{2, -0.5, 1.5}});
    histos.add("EventQA/hMult", "Multiplicity", kTH1F, {axisMult});
    histos.add("EventQA/hPvz", ";pvz;Entries", kTH1F, {{30, -15, 15}});
    histos.add("EventQA/hPvz_aftersel", ";pvzaftersel;Entries", kTH1F, {{30, -15, 15}});
    histos.add("EventQA/hMultFT0vsTPC", ";centFT0C;multNTracksPVeta1", kTH2F, {{100, 0, 100}, {300, 0, 300}});
    histos.add("EventQA/TPCPID", "EventQA/TPCPID", kTH2F, {{20, 0, 2}, {300, 0, 300}});
    // after cut
    histos.add("EventQA/hMult_cut", "Multiplicity", kTH1F, {axisMult});

    // QA and THn Histograms
    histos.add("hTriggerPtResolution", ";p_{T}^{reconstructed} (GeV/c); p_{T}^{generated} (GeV/c)", kTH2F, {axisPtQA, axisPtQA});
    histos.add("hTriggerPrimaryEtaVsPt", "hTriggerPrimaryEtaVsPt", kTH3F, {axisPtQA, axisEta, axisPhi});
    histos.add("hTrackEtaVsPtVsPhi", "hTrackEtaVsPtVsPhi", kTH3F, {axisPtQA, axisEta, axisPhi});
    histos.add("hTrackAttempt", "Attempt", kTH3F, {axisPtQA, axisEta, axisPhi});

    bool hStrange = false;
    for (int i = 0; i < 9; i++) {
      if (TESTBIT(doCorrelation, i)) {
        if (doCentSelectedhigh) {
          histos.add(fmt::format("h{}EtaVsPtVsPhi/CentHigh", kParticlenames[i]).c_str(), "", kTH3F, {axisPtQA, axisEta, axisPhi});
          histos.add(fmt::format("h3d{}Spectrum/CentHigh", kParticlenames[i]).c_str(), fmt::format("h3d{}Spectrum/CentHigh", kParticlenames[i]).c_str(), kTH3F, {axisPtQA, axisMult, axisMassNSigma});
          histos.add(fmt::format("h3d{}SpectrumY/CentHigh", kParticlenames[i]).c_str(), fmt::format("h3d{}SpectrumY/CentHigh", kParticlenames[i]).c_str(), kTH3F, {axisPtQA, axisMult, axisMassNSigma});
          histos.add(fmt::format("sameEvent/Signal/CentHigh/{}", kParticlenames[i]).c_str(), "", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisVtxZNDim, axisMultNDim});
        }
        if (doCentSelectedmed) {
          histos.add(fmt::format("h{}EtaVsPtVsPhi/CentMed", kParticlenames[i]).c_str(), "", kTH3F, {axisPtQA, axisEta, axisPhi});
          histos.add(fmt::format("h3d{}Spectrum/CentMed", kParticlenames[i]).c_str(), fmt::format("h3d{}Spectrum/CentMed", kParticlenames[i]).c_str(), kTH3F, {axisPtQA, axisMult, axisMassNSigma});
          histos.add(fmt::format("h3d{}SpectrumY/CentMed", kParticlenames[i]).c_str(), fmt::format("h3d{}SpectrumY/CentMed", kParticlenames[i]).c_str(), kTH3F, {axisPtQA, axisMult, axisMassNSigma});
          histos.add(fmt::format("sameEvent/Signal/CentMed/{}", kParticlenames[i]).c_str(), "", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisVtxZNDim, axisMultNDim});
        }
        if (doCentSelectedlow) {
          histos.add(fmt::format("h{}EtaVsPtVsPhi/CentLow", kParticlenames[i]).c_str(), "", kTH3F, {axisPtQA, axisEta, axisPhi});
          histos.add(fmt::format("h3d{}Spectrum/CentLow", kParticlenames[i]).c_str(), fmt::format("h3d{}Spectrum/CentLow", kParticlenames[i]).c_str(), kTH3F, {axisPtQA, axisMult, axisMassNSigma});
          histos.add(fmt::format("h3d{}SpectrumY/CentLow", kParticlenames[i]).c_str(), fmt::format("h3d{}SpectrumY/CentLow", kParticlenames[i]).c_str(), kTH3F, {axisPtQA, axisMult, axisMassNSigma});
          histos.add(fmt::format("sameEvent/Signal/CentLow/{}", kParticlenames[i]).c_str(), "", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisVtxZNDim, axisMultNDim});
        }

        if (i < 7) {
          hStrange = true;
          if (doCentSelectedhigh)
            histos.add(fmt::format("h{}EtaVsPtVsPhiBg/CentHigh", kParticlenames[i]).c_str(), "", kTH3F, {axisPtQA, axisEta, axisPhi});
          if (doCentSelectedmed)
            histos.add(fmt::format("h{}EtaVsPtVsPhiBg/CentMed", kParticlenames[i]).c_str(), "", kTH3F, {axisPtQA, axisEta, axisPhi});
          if (doCentSelectedlow)
            histos.add(fmt::format("h{}EtaVsPtVsPhiBg/CentLow", kParticlenames[i]).c_str(), "", kTH3F, {axisPtQA, axisEta, axisPhi});
        }
      }
    }
    if (hStrange) {
      histos.addClone("sameEvent/Signal/", "sameEvent/LeftBg/");
      histos.addClone("sameEvent/Signal/", "sameEvent/RightBg/");
    }

    LOGF(info, "Init THnFs done");
    // mixed-event correlation functions
    if (doprocessMixedEventHV0s || doprocessMixedEventHCascades) {
      histos.addClone("sameEvent/", "mixedEvent/");
    }
    if (doprocessSameEventHV0s)
      histos.add("sameEvent/TriggerParticlesV0", "TriggersV0", kTH2F, {axisPtQA, axisMult});
    if (doprocessSameEventHCascades)
      histos.add("sameEvent/TriggerParticlesCascade", "TriggersCascade", kTH2F, {axisPtQA, axisMult});

    // MC generated plots
    if (doprocessMCGenerated) {
      histos.add("Generated/hTrigger", "", kTH2F, {axisPtQA, axisEta});
      for (int i = 0; i < 8; i++) {
        histos.add(fmt::format("Generated/h{}", kParticlenames[i]).c_str(), "", kTH2F, {axisPtQA, axisEta});
      }

      histos.addClone("Generated/", "GeneratedWithPV/");

      // histograms within |y|<0.5, vs multiplicity
      for (int i = 0; i < 8; i++) {
        histos.add(fmt::format("GeneratedWithPV/h{}_MidYVsMult", kParticlenames[i]).c_str(), "", kTH2F, {axisPtQA, axisMult});
        histos.add(fmt::format("GeneratedWithPV/h{}_MidYVsMult_TwoPVsOrMore", kParticlenames[i]).c_str(), "", kTH2F, {axisPtQA, axisMult});
      }
    }
    if (doprocessClosureTest) {
      for (int i = 0; i < 8; i++) {
        if (TESTBIT(doCorrelation, i))
          histos.add(fmt::format("ClosureTest/sameEvent/{}", kParticlenames[i]).c_str(), "", kTHnF, {axisDeltaPhiNDim, axisDeltaEtaNDim, axisPtAssocNDim, axisVtxZNDim, axisMultNDim});
        if (TESTBIT(doCorrelation, i))
          histos.add(fmt::format("ClosureTest/h{}", kParticlenames[i]).c_str(), "", kTH3F, {axisPtQA, axisEta, axisPhi});
      }
      histos.add("ClosureTest/hTrigger", "Trigger Tracks", kTH3F, {axisPtQA, axisEta, axisPhi});
    }
    // visual inspection of sizes
    histos.print();
    // initialize CCDB *only* if efficiency correction requested
    // skip if not requested, saves a bit of time
    if (applyEfficiencyCorrection) {
      ccdb->setURL(ccdburl);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      ccdb->setFatalWhenNull(false);
    }
  }

  // event selection
  template <typename TCollision>
  bool eventSelected(TCollision collision, bool fillHists)
  {
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 0.5 /* all collisions */);

    // Perform basic event selection
    if (!collision.sel8()) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 1.5 /* collisions  after sel8*/);

    float vtxz = -999;
    if (collision.numContrib() > 1) {
      vtxz = collision.posZ();
      float zRes = std::sqrt(collision.covZZ());
      if (zRes > 0.25 && collision.numContrib() < 20)
        vtxz = -999;
    }
    if (std::fabs(vtxz) > zVertexCut)
      return false;
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 2.5 /* collisions  after sel pvz sel*/);

    if (!collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      // cut time intervals with dead ITS staves
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 3.5 /* collisions  after cut time intervals with dead ITS staves*/);

    if (!collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      // removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference
      // use this cut at low multiplicities with caution
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 4.5 /* removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference*/);

    auto occupancy = collision.trackOccupancyInTimeRange();
    if (occupancy < cfgCutOccupancyLow || occupancy > cfgCutOccupancyHigh)
      return false;
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 5.5 /* Below min occupancy and Above max occupancy*/);

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
      histos.fill(HIST("hEventSelection"), 6.5 /* reject collisions close to Time Frame borders*/);

    if (!collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      // reject events affected by the ITS ROF border
      // O2-4309
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 7.5 /* reject events affected by the ITS ROF border*/);

    if (!collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      // rejects collisions which are associated with the same "found-by-T0" bunch crossing
      // https://indico.cern.ch/event/1396220/#1-event-selection-with-its-rof
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 8.5 /* rejects collisions which are associated with the same "found-by-T0" bunch crossing*/);
    return true;
  }

  void processSameEventHV0s(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::PVMults>::iterator const& collision,
                            aod::AssocV0s const& associatedV0s, aod::TriggerTracks const& triggerTracks,
                            V0DatasWithoutTrackX const&, aod::V0sLinked const&, TracksComplete const&, aod::BCsWithTimestamps const&)
  {
    const auto cent = collision.centFT0C();

    // ________________________________________________
    // Perform event selection
    if (!eventSelected(collision, true) || !doEventSelected) {
      return;
    }
    // ________________________________________________
    if (!doprocessSameEventHCascades) {
      histos.fill(HIST("MixingQA/hSECollisionBins"), colBinning.getBin({collision.posZ(), cent}));
      histos.fill(HIST("EventQA/hMult"), cent);
      histos.fill(HIST("EventQA/hPvz"), collision.posZ());
      histos.fill(HIST("EventQA/hMultFT0vsTPC"), cent, collision.multNTracksPVeta1());
    }
    // Do basic QA

    if (applyEfficiencyCorrection) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
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
      if (!V0Selected(v0Data, true))
        continue;
      uint64_t selMap = V0selectionBitmap(v0Data, collision.posX(), collision.posY(), collision.posZ(), true);

      static_for<0, 2>([&](auto i) {
        constexpr int Index = i.value;
        float efficiency = 1.0f;
        if (applyEfficiencyCorrection) {
          efficiency = hEfficiencyV0[Index]->Interpolate(v0Data.pt(), v0Data.eta());
        }
        float weight = applyEfficiencyCorrection ? 1. / efficiency : 1.0f;
        if (v0.compatible(Index, systCuts.dEdxCompatibility) && (!doMCassociation || v0.mcTrue(Index)) && (!doAssocPhysicalPrimary || v0.mcPhysicalPrimary()) && (!applyEfficiencyCorrection || efficiency != 0)) {
          if (TESTBIT(doCorrelation, Index) && TESTBIT(selMap, Index) && TESTBIT(selMap, Index + 3) && TESTBIT(selMap, Index + 6) && TESTBIT(selMap, Index + 9)) {
            if (doCentSelectedhigh && cent < centhighMax && cent > centhighMin)
              histos.fill(HIST("h3d") + HIST(kV0names[Index]) + HIST("Spectrum/CentHigh"), v0Data.pt(), cent, v0.invMassNSigma(Index), weight);
            if (doCentSelectedmed && cent < centmedMax && cent > centmedMin)
              histos.fill(HIST("h3d") + HIST(kV0names[Index]) + HIST("Spectrum/CentMed"), v0Data.pt(), cent, v0.invMassNSigma(Index), weight);
            if (doCentSelectedlow && cent < centlowMax && cent > centlowMin)
              histos.fill(HIST("h3d") + HIST(kV0names[Index]) + HIST("Spectrum/CentLow"), v0Data.pt(), cent, v0.invMassNSigma(Index), weight);
            if (std::abs(v0Data.rapidity(Index)) < 0.5) {
              if (doCentSelectedhigh && cent < centhighMax && cent > centhighMin)
                histos.fill(HIST("h3d") + HIST(kV0names[Index]) + HIST("SpectrumY/CentHigh"), v0Data.pt(), cent, v0.invMassNSigma(Index), weight);
              if (doCentSelectedmed && cent < centmedMax && cent > centmedMin)
                histos.fill(HIST("h3d") + HIST(kV0names[Index]) + HIST("SpectrumY/CentMed"), v0Data.pt(), cent, v0.invMassNSigma(Index), weight);
              if (doCentSelectedlow && cent < centlowMax && cent > centlowMin)
                histos.fill(HIST("h3d") + HIST(kV0names[Index]) + HIST("SpectrumY/CentLow"), v0Data.pt(), cent, v0.invMassNSigma(Index), weight);
            }
            if ((-massWindowConfigurations.maxBgNSigma < v0.invMassNSigma(Index) && v0.invMassNSigma(Index) < -massWindowConfigurations.minBgNSigma) || (+massWindowConfigurations.minBgNSigma < v0.invMassNSigma(Index) && v0.invMassNSigma(Index) < +massWindowConfigurations.maxBgNSigma)) {
              if (doCentSelectedhigh && cent < centhighMax && cent > centhighMin)
                histos.fill(HIST("h") + HIST(kV0names[Index]) + HIST("EtaVsPtVsPhiBg/CentHigh"), v0Data.pt(), v0Data.eta(), v0Data.phi(), weight);
              if (doCentSelectedmed && cent < centmedMax && cent > centmedMin)
                histos.fill(HIST("h") + HIST(kV0names[Index]) + HIST("EtaVsPtVsPhiBg/CentMed"), v0Data.pt(), v0Data.eta(), v0Data.phi(), weight);
              if (doCentSelectedlow && cent < centlowMax && cent > centlowMin)
                histos.fill(HIST("h") + HIST(kV0names[Index]) + HIST("EtaVsPtVsPhiBg/CentLow"), v0Data.pt(), v0Data.eta(), v0Data.phi(), weight);
            }
            if (-massWindowConfigurations.maxPeakNSigma < v0.invMassNSigma(Index) && v0.invMassNSigma(Index) < +massWindowConfigurations.maxPeakNSigma) {
              if (doCentSelectedhigh && cent < centhighMax && cent > centhighMin)
                histos.fill(HIST("h") + HIST(kV0names[Index]) + HIST("EtaVsPtVsPhi/CentHigh"), v0Data.pt(), v0Data.eta(), v0Data.phi(), weight);
              if (doCentSelectedmed && cent < centmedMax && cent > centmedMin)
                histos.fill(HIST("h") + HIST(kV0names[Index]) + HIST("EtaVsPtVsPhi/CentMed"), v0Data.pt(), v0Data.eta(), v0Data.phi(), weight);
              if (doCentSelectedlow && cent < centlowMax && cent > centlowMin)
                histos.fill(HIST("h") + HIST(kV0names[Index]) + HIST("EtaVsPtVsPhi/CentLow"), v0Data.pt(), v0Data.eta(), v0Data.phi(), weight);
            }
            // proper lifetime
            histos.fill(HIST("V0selection/lifetime") + HIST(kV0names[Index]) + HIST("_after"), v0Data.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short);
            // TPC PID
            histos.fill(HIST("V0selection/tpcNsigmaPos") + HIST(kV0names[Index]) + HIST("_after"), postrack.tpcNSigmaPr());
            histos.fill(HIST("V0selection/tpcNsigmaNeg") + HIST(kV0names[Index]) + HIST("_after"), negtrack.tpcNSigmaPi());
            // DCA daughter to prim.vtx
            histos.fill(HIST("V0selection/dcapostopv") + HIST(kV0names[Index]) + HIST("_after"), v0Data.dcapostopv());
            histos.fill(HIST("V0selection/dcanegtopv") + HIST(kV0names[Index]) + HIST("_after"), v0Data.dcanegtopv());
            // Armenteros-Podolanski
            histos.fill(HIST("V0selection/hArmenteros-Podolanski") + HIST(kV0names[Index]) + HIST("_after"), v0Data.alpha(), v0Data.qtarm());
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
    fillCorrelationsV0(triggerTracks, associatedV0s, false, collision.posX(), collision.posY(), collision.posZ(), cent);
  }

  void processSameEventHCascades(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::PVMults>::iterator const& collision,
                                 aod::AssocV0s const&, aod::AssocCascades const& associatedCascades, aod::TriggerTracks const& triggerTracks,
                                 V0DatasWithoutTrackX const&, aod::V0sLinked const&, aod::CascDatas const&, TracksComplete const&, aod::BCsWithTimestamps const&)
  {
    const auto cent = collision.centFT0C();
    // ________________________________________________
    // skip if desired trigger not found
    /*
    if (triggerPresenceMap.size() > 0 && !TESTBIT(triggerPresenceMap[collision.globalIndex()], triggerBinToSelect)) {
      return;
    }
    */

    // ________________________________________________
    // Perform event selection
    if (!eventSelected(collision, true) || !doEventSelected) {
      return;
    }

    // ________________________________________________
    histos.fill(HIST("MixingQA/hSECollisionBins"), colBinning.getBin({collision.posZ(), cent}));
    histos.fill(HIST("EventQA/hMult"), cent);
    histos.fill(HIST("EventQA/hPvz"), collision.posZ());
    // Do basic QA
    if (applyEfficiencyCorrection) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initEfficiencyFromCCDB(bc);
    }
    TH2F* hEfficiencyCascade[4];
    hEfficiencyCascade[0] = hEfficiencyXiMinus;
    hEfficiencyCascade[1] = hEfficiencyXiPlus;
    hEfficiencyCascade[2] = hEfficiencyOmegaMinus;
    hEfficiencyCascade[3] = hEfficiencyOmegaPlus;

    for (auto const& casc : associatedCascades) {
      auto assoc = casc.cascData();
      auto postrack = assoc.posTrack_as<TracksComplete>();
      auto negtrack = assoc.negTrack_as<TracksComplete>();
      auto bachtrack = assoc.bachelor_as<TracksComplete>();
      float cascpos = std::hypot(assoc.x() - collision.posX(), assoc.y() - collision.posY(), assoc.z() - collision.posZ());
      float cascptotmom = std::hypot(assoc.px(), assoc.py(), assoc.pz());
      float ctauXi = o2::constants::physics::MassXiMinus * cascpos / ((cascptotmom + 1e-13) * ctauxiPDG);
      float ctauOmega = o2::constants::physics::MassOmegaMinus * cascpos / ((cascptotmom + 1e-13) * ctauomegaPDG);
      //---] syst cuts [---
      if (!CascadeSelected(assoc, collision.posX(), collision.posY(), collision.posZ(), true))
        continue;
      uint64_t CascselMap = CascadeselectionBitmap(assoc, collision.posX(), collision.posY(), collision.posZ(), true);

      static_for<0, 3>([&](auto i) {
        constexpr int Index = i.value;
        float efficiency = 1.0f;
        if (applyEfficiencyCorrection) {
          efficiency = hEfficiencyCascade[Index]->Interpolate(assoc.pt(), assoc.eta());
        }
        float weight = applyEfficiencyCorrection ? 1. / efficiency : 1.0f;
        if (casc.compatible(Index, systCuts.dEdxCompatibility) && (!doMCassociation || casc.mcTrue(Index)) && (!doAssocPhysicalPrimary || casc.mcPhysicalPrimary()) && (!applyEfficiencyCorrection || efficiency != 0)) {
          if (TESTBIT(doCorrelation, Index + 3) && TESTBIT(CascselMap, Index) && TESTBIT(CascselMap, Index + 4) && TESTBIT(CascselMap, Index + 8) && TESTBIT(CascselMap, Index + 12) && TESTBIT(CascselMap, Index + 16)) {
            if (doCentSelectedhigh && cent < centhighMax && cent > centhighMin)
              histos.fill(HIST("h3d") + HIST(kCascadenames[Index]) + HIST("Spectrum/CentHigh"), assoc.pt(), cent, casc.invMassNSigma(Index), weight);
            if (doCentSelectedmed && cent < centmedMax && cent > centmedMin)
              histos.fill(HIST("h3d") + HIST(kCascadenames[Index]) + HIST("Spectrum/CentMed"), assoc.pt(), cent, casc.invMassNSigma(Index), weight);
            if (doCentSelectedlow && cent < centlowMax && cent > centlowMin)
              histos.fill(HIST("h3d") + HIST(kCascadenames[Index]) + HIST("Spectrum/CentLow"), assoc.pt(), cent, casc.invMassNSigma(Index), weight);
            if (std::abs(assoc.rapidity(Index)) < 0.5) {
              if (doCentSelectedhigh && cent < centhighMax && cent > centhighMin)
                histos.fill(HIST("h3d") + HIST(kCascadenames[Index]) + HIST("SpectrumY/CentHigh"), assoc.pt(), cent, casc.invMassNSigma(Index), weight);
              if (doCentSelectedmed && cent < centmedMax && cent > centmedMin)
                histos.fill(HIST("h3d") + HIST(kCascadenames[Index]) + HIST("SpectrumY/CentMed"), assoc.pt(), cent, casc.invMassNSigma(Index), weight);
              if (doCentSelectedlow && cent < centlowMax && cent > centlowMin)
                histos.fill(HIST("h3d") + HIST(kCascadenames[Index]) + HIST("SpectrumY/CentLow"), assoc.pt(), cent, casc.invMassNSigma(Index), weight);
            }
            if ((-massWindowConfigurations.maxBgNSigma < casc.invMassNSigma(Index) && casc.invMassNSigma(Index) < -massWindowConfigurations.minBgNSigma) || (+massWindowConfigurations.minBgNSigma < casc.invMassNSigma(Index) && casc.invMassNSigma(Index) < +massWindowConfigurations.maxBgNSigma)) {
              if (doCentSelectedhigh && cent < centhighMax && cent > centhighMin)
                histos.fill(HIST("h") + HIST(kCascadenames[Index]) + HIST("EtaVsPtVsPhiBg/CentHigh"), assoc.pt(), assoc.eta(), assoc.phi(), weight);
              if (doCentSelectedmed && cent < centmedMax && cent > centmedMin)
                histos.fill(HIST("h") + HIST(kCascadenames[Index]) + HIST("EtaVsPtVsPhiBg/CentMed"), assoc.pt(), assoc.eta(), assoc.phi(), weight);
              if (doCentSelectedlow && cent < centlowMax && cent > centlowMin)
                histos.fill(HIST("h") + HIST(kCascadenames[Index]) + HIST("EtaVsPtVsPhiBg/CentLow"), assoc.pt(), assoc.eta(), assoc.phi(), weight);
            }
            if (-massWindowConfigurations.maxPeakNSigma < casc.invMassNSigma(Index) && casc.invMassNSigma(Index) < +massWindowConfigurations.maxPeakNSigma) {
              if (doCentSelectedhigh && cent < centhighMax && cent > centhighMin)
                histos.fill(HIST("h") + HIST(kCascadenames[Index]) + HIST("EtaVsPtVsPhi/CentHigh"), assoc.pt(), assoc.eta(), assoc.phi(), weight);
              if (doCentSelectedmed && cent < centmedMax && cent > centmedMin)
                histos.fill(HIST("h") + HIST(kCascadenames[Index]) + HIST("EtaVsPtVsPhi/CentMed"), assoc.pt(), assoc.eta(), assoc.phi(), weight);
              if (doCentSelectedlow && cent < centlowMax && cent > centlowMin)
                histos.fill(HIST("h") + HIST(kCascadenames[Index]) + HIST("EtaVsPtVsPhi/CentLow"), assoc.pt(), assoc.eta(), assoc.phi(), weight);
            }
            // TPC PID  HIST("Cascadeselection/tpcNSigmaBach") + HIST(kCascadenames[Index]) + HIST("_after")
            histos.fill(HIST("Cascadeselection/tpcNSigmaPos") + HIST(kCascadenames[Index]) + HIST("_after"), postrack.tpcNSigmaPr());
            histos.fill(HIST("Cascadeselection/tpcNSigmaNeg") + HIST(kCascadenames[Index]) + HIST("_after"), negtrack.tpcNSigmaPi());
            histos.fill(HIST("Cascadeselection/tpcNSigmaBach") + HIST(kCascadenames[Index]) + HIST("_after"), bachtrack.tpcNSigmaPi());
            // DCA daughter to prim.vtx   HIST("Cascadeselection/ycut") + HIST(kCascadenames[Index]) + HIST("_after")
            histos.fill(HIST("Cascadeselection/dcapostopv") + HIST(kCascadenames[Index]) + HIST("_after"), assoc.dcapostopv());
            histos.fill(HIST("Cascadeselection/dcanegtopv") + HIST(kCascadenames[Index]) + HIST("_after"), assoc.dcanegtopv());
            histos.fill(HIST("Cascadeselection/dcabachtopv") + HIST(kCascadenames[Index]) + HIST("_after"), assoc.dcabachtopv());
            // y cut
            histos.fill(HIST("Cascadeselection/ycut") + HIST(kCascadenames[Index]) + HIST("_after"), assoc.yXi());
            // comopeting casc.rej and life time
            if (doCorrelationXiMinus) {
              histos.fill(HIST("Cascadeselection/rejcompXiMinus_after"), assoc.mOmega() - o2::constants::physics::MassOmegaMinus);
              histos.fill(HIST("Cascadeselection/proplifetimeXiMinus_after"), ctauXi);
            }
            if (doCorrelationXiPlus) {
              histos.fill(HIST("Cascadeselection/rejcompXiPlus_after"), assoc.mOmega() - o2::constants::physics::MassOmegaMinus);
              histos.fill(HIST("Cascadeselection/proplifetimeXiPlus_after"), ctauXi);
            }
            if (doCorrelationOmegaMinus) {
              histos.fill(HIST("Cascadeselection/rejcompOmegaMinus_after"), assoc.mXi() - o2::constants::physics::MassXiMinus);
              histos.fill(HIST("Cascadeselection/proplifetimeOmegaMinus_after"), ctauOmega);
            }
            if (doCorrelationOmegaPlus) {
              histos.fill(HIST("Cascadeselection/rejcompOmegaPlus_after"), assoc.mXi() - o2::constants::physics::MassXiMinus);
              histos.fill(HIST("Cascadeselection/proplifetimeOmegaPlus_after"), ctauOmega);
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
      histos.fill(HIST("hTrackEtaVsPtVsPhi"), track.pt(), track.eta(), track.phi());
    }

    // ________________________________________________
    // Do hadron - cascade correlations
    fillCorrelationsCascade(triggerTracks, associatedCascades, false, collision.posX(), collision.posY(), collision.posZ(), collision.centFT0C());
  }

  void processMixedEventHV0s(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::PVMults> const& collisions,
                             aod::AssocV0s const& associatedV0s, aod::TriggerTracks const& triggerTracks,
                             V0DatasWithoutTrackX const&, aod::V0sLinked const&, TracksComplete const&, aod::BCsWithTimestamps const&)
  {
    for (auto const& [collision1, collision2] : soa::selfCombinations(colBinning, mixingParameter, -1, collisions, collisions)) {
      // ________________________________________________
      if (applyEfficiencyCorrection) {
        auto bc = collision1.bc_as<aod::BCsWithTimestamps>();
        initEfficiencyFromCCDB(bc);
      }
      /*
      // ________________________________________________
      // skip if desired trigger not found
      if (triggerPresenceMap.size() > 0) {
        continue;
      }
      */
      // Perform event selection on both collisions
      if (!eventSelected(collision1, false) || !eventSelected(collision2, false) || !doEventSelected) {
        continue;
      }

      if (!doprocessMixedEventHCascades) {
        if (collision1.globalIndex() == collision2.globalIndex()) {
          histos.fill(HIST("MixingQA/hMixingQA"), 0.0f); // same-collision pair counting
        }
        histos.fill(HIST("MixingQA/hMEpvz1"), collision1.posZ());
        histos.fill(HIST("MixingQA/hMEpvz2"), collision2.posZ());
        histos.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), collision1.centFT0C()}));
      }
      // ________________________________________________
      // Do hadron - V0 correlations
      fillCorrelationsV0(triggerTracks, associatedV0s, true, collision1.posX(), collision1.posY(), collision1.posZ(), collision1.centFT0C());
    }
  }
  void processMixedEventHCascades(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::PVMults> const& collisions,
                                  aod::AssocV0s const&, aod::AssocCascades const& associatedCascades, aod::TriggerTracks const& triggerTracks,
                                  V0DatasWithoutTrackX const&, aod::V0sLinked const&, aod::CascDatas const&, TracksComplete const&, aod::BCsWithTimestamps const&)
  {
    for (auto const& [collision1, collision2] : soa::selfCombinations(colBinning, mixingParameter, -1, collisions, collisions)) {
      // ________________________________________________
      if (applyEfficiencyCorrection) {
        auto bc = collision1.bc_as<aod::BCsWithTimestamps>();
        initEfficiencyFromCCDB(bc);
      }
      /*
      // ________________________________________________
      // skip if desired trigger not found
      if (triggerPresenceMap.size() > 0 && (!TESTBIT(triggerPresenceMap[collision1.globalIndex()], triggerBinToSelect) || !TESTBIT(triggerPresenceMap[collision2.globalIndex()], triggerBinToSelect))) {
        continue;
      }
      */

      // Perform event selection on both collisions
      if (!eventSelected(collision1, false) || !eventSelected(collision2, false) || !doEventSelected) {
        continue;
      }

      if (collision1.globalIndex() == collision2.globalIndex()) {
        histos.fill(HIST("MixingQA/hMixingQA"), 0.0f); // same-collision pair counting
      }

      histos.fill(HIST("MixingQA/hMEpvz1"), collision1.posZ());
      histos.fill(HIST("MixingQA/hMEpvz2"), collision2.posZ());
      histos.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), collision1.centFT0C()}));
      // ________________________________________________
      // Do slicing
      auto slicedTriggerTracks = triggerTracks.sliceBy(collisionSliceTracks, collision1.globalIndex());
      auto slicedAssocCascades = associatedCascades.sliceBy(collisionSliceCascades, collision2.globalIndex());
      // ________________________________________________
      // Do hadron - cascade correlations
      fillCorrelationsCascade(slicedTriggerTracks, slicedAssocCascades, true, collision1.posX(), collision1.posY(), collision1.posZ(), collision1.centFT0C());
    }
  }

  void processMCGenerated(aod::McCollision const& /*mcCollision*/, soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::PVMults>> const& collisions, aod::McParticles const& mcParticles)
  {
    histos.fill(HIST("hClosureTestEventCounter"), 2.5f);

    for (auto const& mcParticle : mcParticles) {
      double geta = mcParticle.eta();
      if (std::abs(geta) > 0.8f) {
        continue;
      }
      double gpt = mcParticle.pt();
      if (std::abs(mcParticle.pdgCode()) == 211 || std::abs(mcParticle.pdgCode()) == 321 || std::abs(mcParticle.pdgCode()) == 2212 || std::abs(mcParticle.pdgCode()) == 11 || std::abs(mcParticle.pdgCode()) == 13) {
        if (!doTriggPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
          histos.fill(HIST("hGeneratedQAPtTrigger"), gpt, 0.0f); // step 1: before all selections
        }
      }

      if (!doAssocPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
        if (std::abs(mcParticle.pdgCode()) == 310 && doCorrelationK0Short) {
          histos.fill(HIST("hGeneratedQAPtAssociatedK0"), gpt, 0.0f); // step 1: before all selections
        }
      }
    }

    for (auto const& mcParticle : mcParticles) {
      if (!mcParticle.isPhysicalPrimary())
        continue;
      static_for<0, 7>([&](auto i) {
        constexpr int Index = i.value;
        if (i == 0 || i == 7) {
          if (std::abs(mcParticle.pdgCode()) == kPdgCodes[i])
            histos.fill(HIST("Generated/h") + HIST(kParticlenames[Index]), mcParticle.pt(), mcParticle.eta());
        } else {
          if (mcParticle.pdgCode() == kPdgCodes[i])
            histos.fill(HIST("Generated/h") + HIST(kParticlenames[Index]), mcParticle.pt(), mcParticle.eta());
        }
      });
    }
    if (collisions.size() < 1)
      return;

    // determine best collision properties
    int biggestNContribs = -1;
    int bestCollisionFT0Cpercentile = -1;
    float bestCollisionVtxZ = 0.0f;
    bool bestCollisionSel8 = false;
    // bool bestCollisionINELgtZERO = false;
    // uint32_t bestCollisionTriggerPresenceMap = 0;

    for (auto const& collision : collisions) {
      if (biggestNContribs < collision.numContrib()) {
        biggestNContribs = collision.numContrib();
        bestCollisionFT0Cpercentile = collision.centFT0C();
        bestCollisionSel8 = collision.sel8();
        bestCollisionVtxZ = collision.posZ();
        // bestCollisionINELgtZERO = collision.isInelGt0();
        /*
        if (triggerPresenceMap.size() > 0)
          bestCollisionTriggerPresenceMap = triggerPresenceMap[collision.globalIndex()]; */
      }
    }

    if (collisions.size() > 1) {
      for (auto const& mcParticle : mcParticles) {
        if (!mcParticle.isPhysicalPrimary())
          continue;
        if (std::abs(mcParticle.y()) > 0.5)
          continue;
        static_for<0, 7>([&](auto i) {
          constexpr int Index = i.value;
          if (i == 0 || i == 7) {
            if (std::abs(mcParticle.pdgCode()) == kPdgCodes[i])
              histos.fill(HIST("GeneratedWithPV/h") + HIST(kParticlenames[Index]) + HIST("_MidYVsMult_TwoPVsOrMore"), mcParticle.pt(), bestCollisionFT0Cpercentile);
          } else {
            if (mcParticle.pdgCode() == kPdgCodes[i])
              histos.fill(HIST("GeneratedWithPV/h") + HIST(kParticlenames[Index]) + HIST("_MidYVsMult_TwoPVsOrMore"), mcParticle.pt(), bestCollisionFT0Cpercentile);
          }
        });
      }
    }

    // do selections on best collision
    if (!bestCollisionSel8)
      return;
    if (std::abs(bestCollisionVtxZ) > 10.0f)
      return;

    histos.fill(HIST("hClosureTestEventCounter"), 3.5f);

    for (auto const& mcParticle : mcParticles) {
      double geta = mcParticle.eta();
      if (std::abs(geta) > 0.8f) {
        continue;
      }
      double gpt = mcParticle.pt();
      if (std::abs(mcParticle.pdgCode()) == 211 || std::abs(mcParticle.pdgCode()) == 321 || std::abs(mcParticle.pdgCode()) == 2212 || std::abs(mcParticle.pdgCode()) == 11 || std::abs(mcParticle.pdgCode()) == 13) {
        if (!doTriggPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
          histos.fill(HIST("hGeneratedQAPtTrigger"), gpt, 1.0f); // step 2: after event selection
        }
      }

      if (!doAssocPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
        if (std::abs(mcParticle.pdgCode()) == 310 && doCorrelationK0Short) {
          histos.fill(HIST("hGeneratedQAPtAssociatedK0"), gpt, 1.0f); // step 2: before all selections
        }
      }
    }

    for (auto const& mcParticle : mcParticles) {
      if (!mcParticle.isPhysicalPrimary()) {
        continue;
      }
      double geta = mcParticle.eta();
      double gpt = mcParticle.pt();
      if (std::abs(mcParticle.pdgCode()) == 211 || std::abs(mcParticle.pdgCode()) == 321 || std::abs(mcParticle.pdgCode()) == 2212 || std::abs(mcParticle.pdgCode()) == 11 || std::abs(mcParticle.pdgCode()) == 13)
        histos.fill(HIST("GeneratedWithPV/hTrigger"), gpt, geta);
      static_for<0, 7>([&](auto i) {
        constexpr int Index = i.value;
        if (i == 0 || i == 7) {
          if (std::abs(mcParticle.pdgCode()) == kPdgCodes[i]) {
            histos.fill(HIST("GeneratedWithPV/h") + HIST(kParticlenames[Index]), gpt, geta);
            if (std::abs(mcParticle.y()) < 0.5)
              histos.fill(HIST("GeneratedWithPV/h") + HIST(kParticlenames[Index]) + HIST("_MidYVsMult"), gpt, bestCollisionFT0Cpercentile);
          }

        } else {
          if (mcParticle.pdgCode() == kPdgCodes[i]) {
            histos.fill(HIST("GeneratedWithPV/h") + HIST(kParticlenames[Index]), gpt, geta);
            if (std::abs(mcParticle.y()) < 0.5)
              histos.fill(HIST("GeneratedWithPV/h") + HIST(kParticlenames[Index]) + HIST("_MidYVsMult"), gpt, bestCollisionFT0Cpercentile);
          }
        }
      });
    }
  }
  void processClosureTest(aod::McCollision const& /*mcCollision*/, soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::PVMults>> const& recCollisions, aod::McParticles const& mcParticles)
  {
    std::vector<uint32_t> triggerIndices;
    std::vector<std::vector<uint32_t>> associatedIndices;
    std::vector<uint32_t> k0ShortIndices;
    std::vector<uint32_t> lambdaIndices;
    std::vector<uint32_t> antiLambdaIndices;
    std::vector<uint32_t> xiMinusIndices;
    std::vector<uint32_t> xiPlusIndices;
    std::vector<uint32_t> omegaMinusIndices;
    std::vector<uint32_t> omegaPlusIndices;

    for (auto const& mcParticle : mcParticles) {
      double geta = mcParticle.eta();
      if (std::abs(geta) > 0.8f) {
        continue;
      }
      double gpt = mcParticle.pt();
      if (std::abs(mcParticle.pdgCode()) == 211 || std::abs(mcParticle.pdgCode()) == 321 || std::abs(mcParticle.pdgCode()) == 2212 || std::abs(mcParticle.pdgCode()) == 11 || std::abs(mcParticle.pdgCode()) == 13) {
        if (!doTriggPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
          histos.fill(HIST("hClosureQAPtTrigger"), gpt, 0.0f); // step 1: no event selection whatsoever
        }
      }

      if (!doAssocPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
        if (mcParticle.pdgCode() == 310 && doCorrelationK0Short) {
          histos.fill(HIST("hClosureQAPtAssociatedK0"), gpt, 0.0f); // step 1: no event selection whatsoever
        }
      }
    }

    histos.fill(HIST("hClosureTestEventCounter"), 0.5f);

    int bestCollisionFT0Cpercentile = -1;
    float bestCollisionVtxZ = 0.0f;
    bool bestCollisionSel8 = false;
    int biggestNContribs = -1;

    for (auto const& recCollision : recCollisions) {
      if (biggestNContribs < recCollision.numContrib()) {
        biggestNContribs = recCollision.numContrib();
        bestCollisionFT0Cpercentile = recCollision.centFT0C();
        bestCollisionSel8 = recCollision.sel8();
        bestCollisionVtxZ = recCollision.posZ();
      }
    }
    // ________________________________________________

    if (doGenEventSelection) {
      if (!bestCollisionSel8)
        return;
      if (std::abs(bestCollisionVtxZ) > zVertexCut)
        return;
      if (bestCollisionFT0Cpercentile > axisRanges[5][1] || bestCollisionFT0Cpercentile < axisRanges[5][0])
        return;
    }

    histos.fill(HIST("hClosureTestEventCounter"), 1.5f);

    for (auto const& mcParticle : mcParticles) {
      double geta = mcParticle.eta();
      if (std::abs(geta) > 0.8f) {
        continue;
      }
      double gpt = mcParticle.pt();
      if (std::abs(mcParticle.pdgCode()) == 211 || std::abs(mcParticle.pdgCode()) == 321 || std::abs(mcParticle.pdgCode()) == 2212 || std::abs(mcParticle.pdgCode()) == 11 || std::abs(mcParticle.pdgCode()) == 13) {
        if (!doTriggPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
          histos.fill(HIST("hClosureQAPtTrigger"), gpt, 1.0f); // step 2: after event selection
        }
      }

      if (!doAssocPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
        if (mcParticle.pdgCode() == 310 && doCorrelationK0Short) {
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
      if (std::abs(geta) > 0.8f) {
        continue;
      }
      if (std::abs(mcParticle.pdgCode()) == 211 || std::abs(mcParticle.pdgCode()) == 321 || std::abs(mcParticle.pdgCode()) == 2212 || std::abs(mcParticle.pdgCode()) == 11 || std::abs(mcParticle.pdgCode()) == 13) {
        if (!doTriggPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
          triggerIndices.emplace_back(iteratorNum);
          histos.fill(HIST("ClosureTest/hTrigger"), gpt, geta, gphi);
        }
      }
      if (!doAssocPhysicalPrimary || mcParticle.isPhysicalPrimary()) {
        if (mcParticle.pdgCode() == 310 && doCorrelationK0Short) {
          k0ShortIndices.emplace_back(iteratorNum);
          histos.fill(HIST("ClosureTest/hK0Short"), gpt, geta, gphi);
        }
        if (mcParticle.pdgCode() == 3122 && doCorrelationLambda) {
          lambdaIndices.emplace_back(iteratorNum);
          histos.fill(HIST("ClosureTest/hLambda"), gpt, geta, gphi);
        }
        if (mcParticle.pdgCode() == -3122 && doCorrelationAntiLambda) {
          antiLambdaIndices.emplace_back(iteratorNum);
          histos.fill(HIST("ClosureTest/hAntiLambda"), gpt, geta, gphi);
        }
        if (mcParticle.pdgCode() == 3312 && doCorrelationXiMinus) {
          xiMinusIndices.emplace_back(iteratorNum);
          histos.fill(HIST("ClosureTest/hXiMinus"), gpt, geta, gphi);
        }
        if (mcParticle.pdgCode() == -3312 && doCorrelationXiPlus) {
          xiPlusIndices.emplace_back(iteratorNum);
          histos.fill(HIST("ClosureTest/hXiPlus"), gpt, geta, gphi);
        }
        if (mcParticle.pdgCode() == 3334 && doCorrelationOmegaMinus) {
          omegaMinusIndices.emplace_back(iteratorNum);
          histos.fill(HIST("ClosureTest/hOmegaMinus"), gpt, geta, gphi);
        }
        if (mcParticle.pdgCode() == -3334 && doCorrelationOmegaPlus) {
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

    for (std::size_t iTrigger = 0; iTrigger < triggerIndices.size(); iTrigger++) {
      auto triggerParticle = mcParticles.iteratorAt(triggerIndices[iTrigger]);
      // check range of trigger particle
      if (triggerParticle.pt() < axisRanges[3][0]) {
        continue;
      }
      double getatrigger = triggerParticle.eta();
      double gphitrigger = triggerParticle.phi();
      // double pttrigger = triggerParticle.pt();
      auto const& mother = triggerParticle.mothers_first_as<aod::McParticles>();
      auto globalIndex = mother.globalIndex();
      static_for<0, 6>([&](auto i) { // associated loop
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
              histos.fill(HIST("ClosureTest/sameEvent/") + HIST(kParticlenames[Index]), computeDeltaPhi(gphitrigger, gphiassoc), deltaeta, ptassoc, bestCollisionVtxZ, bestCollisionFT0Cpercentile);
          }
        }
      });
    }
  }

  PROCESS_SWITCH(HStrangeCorrelation, processSameEventHV0s, "Process same events, h-V0s", true);
  PROCESS_SWITCH(HStrangeCorrelation, processSameEventHCascades, "Process same events, h-Cascades", false);
  PROCESS_SWITCH(HStrangeCorrelation, processMixedEventHV0s, "Process mixed events, h-V0s", true);
  PROCESS_SWITCH(HStrangeCorrelation, processMixedEventHCascades, "Process mixed events, h-Cascades", false);
  PROCESS_SWITCH(HStrangeCorrelation, processMCGenerated, "Process MC generated", false);
  PROCESS_SWITCH(HStrangeCorrelation, processClosureTest, "Process Closure Test", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HStrangeCorrelation>(cfgc)};
}
