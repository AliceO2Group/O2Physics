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

#include "PWGLF/DataModel/DrvCollisions.h"
#include "PWGLF/DataModel/LFHStrangeCorrelationTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "TF1.h"
#include <TPDGCode.h>

#include <string>
#include <unordered_set>

using namespace o2;
using namespace o2::constants::math;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define BIT_SET(var, nbit) ((var) |= (1 << (nbit)))
#define BIT_CHECK(var, nbit) ((var) & (1 << (nbit)))

struct HStrangeDerivedData {

  // master analysis switches
  Configurable<bool> doPPAnalysis{"doPPAnalysis", true, "if in pp, set to true"};
  Configurable<bool> useParameterization{"useParameterization", true, "ture for parameterization method, false for hist method"};
  // Operational
  Configurable<bool> fillTableOnlyWithCompatible{"fillTableOnlyWithCompatible", true, "pre-apply dE/dx, broad mass window in table filling"};
  Configurable<float> strangedEdxNSigmaLoose{"strangedEdxNSigmaLoose", 5, "Nsigmas for strange decay daughters"};
  Configurable<float> strangedEdxNSigma{"strangedEdxNSigma", 4, "Nsigmas for strange decay daughters"};
  Configurable<float> strangedEdxNSigmaTight{"strangedEdxNSigmaTight", 3, "Nsigmas for strange decay daughters"};
  Configurable<float> nSigmaNearXiMassCenter{"nSigmaNearXiMassCenter", 0, "for Oemga analysis only, to check if candidate mass is around Xi"};

  // used for event selections in Pb-Pb
  Configurable<int> cfgCutOccupancyHigh{"cfgCutOccupancyHigh", 3000, "High cut on TPC occupancy"};
  Configurable<int> cfgCutOccupancyLow{"cfgCutOccupancyLow", 0, "Low cut on TPC occupancy"};

  struct : ConfigurableGroup {
    // event filtering
    Configurable<float> zVertexCut{"zVertexCut", 10, "Cut on PV position"};
    Configurable<bool> selectINELgtZERO{"selectINELgtZERO", true, "select INEL>0 events"};
    Configurable<bool> requireAllGoodITSLayers{"requireAllGoodITSLayers", false, " require that in the event all ITS are good"};
  } eventSelections;

  struct : ConfigurableGroup {
    // Trigger particle selections in phase space
    Configurable<float> triggerEtaMin{"triggerEtaMin", -0.8, "triggeretamin"};
    Configurable<float> triggerEtaMax{"triggerEtaMax", 0.8, "triggeretamax"};
    Configurable<float> triggerPtCutMin{"triggerPtCutMin", 8, "triggerptmin"};
    Configurable<float> triggerPtCutMax{"triggerPtCutMax", 20, "triggerptmax"};

    // Track quality
    Configurable<int> minTPCNCrossedRows{"minTPCNCrossedRows", 70, "Minimum TPC crossed rows"};
    Configurable<bool> triggerRequireITS{"triggerRequireITS", true, "require ITS signal in trigger tracks"};
    Configurable<bool> assocRequireITS{"assocRequireITS", true, "require ITS signal in assoc tracks"};
    Configurable<int> triggerMaxTPCSharedClusters{"triggerMaxTPCSharedClusters", 200, "maximum number of shared TPC clusters (inclusive)"};
    Configurable<bool> triggerRequireL0{"triggerRequireL0", false, "require ITS L0 cluster for trigger"};

    // Associated particle selections in phase space
    Configurable<float> assocEtaMin{"assocEtaMin", -0.8, "triggeretamin"};
    Configurable<float> assocEtaMax{"assocEtaMax", 0.8, "triggeretamax"};
    Configurable<float> assocPtCutMin{"assocPtCutMin", 0.2, "assocptmin"};
    Configurable<float> assocPtCutMax{"assocPtCutMax", 10, "assocptmax"};

    // Associated pion identification
    Configurable<float> pionMinBayesProb{"pionMinBayesProb", 0.95, "minimal Bayesian probability for pion ID"};
    Configurable<float> assocPionNSigmaTPCFOF{"assocPionNSigmaTPCFOF", 3, "minimal n sigma in TOF and TPC for Pion ID"};
    Configurable<float> rejectSigma{"rejectSigma", 1, "n sigma for rejecting pion candidates"};

    // V0 selections
    Configurable<float> v0Cospa{"v0Cospa", 0.97, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0)
    Configurable<float> dcaV0dau{"dcaV0dau", 1.0, "DCA V0 Daughters"};
    Configurable<float> dcaNegtopv{"dcaNegtopv", 0.06, "DCA Neg To PV"};
    Configurable<float> dcaPostopv{"dcaPostopv", 0.06, "DCA Pos To PV"};
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

    // specific selections
    Configurable<float> lambdaCospa{"lambdaCospa", 0.995, "CosPA for lambda"}; // allows for tighter selection for Lambda

    // primary particle DCAxy selections
    // formula: |DCAxy| <  0.004f + (0.013f / pt)
    Configurable<float> dcaXYconstant{"dcaXYconstant", 0.004, "[0] in |DCAxy| < [0]+[1]/pT"};
    Configurable<float> dcaXYpTdep{"dcaXYpTdep", 0.013, "[1] in |DCAxy| < [0]+[1]/pT"};

    // cascade selections
    Configurable<float> cascCospa{"cascCospa", 0.95, "cascCospa"};
    Configurable<float> cascRadius{"cascRadius", 0.5, "cascRadius"};
    Configurable<float> dcaCascdau{"dcaCascdau", 1.0, "dcaCascdau"};
    Configurable<float> dcaBachtopv{"dcaBachtopv", 0.1, "dcaBachtopv"};
    Configurable<float> cascV0masswindow{"cascV0masswindow", 0.01, "cascV0masswindow"};
    Configurable<float> cascMindcav0topv{"cascMindcav0topv", 0.01, "cascMindcav0topv"};
  } systCuts;
  struct : ConfigurableGroup {
    // cascade selections in PbPb
    Configurable<float> cascDcacascdau{"cascDcacascdau", 1.0, "cascDcacascdau"};
    Configurable<float> cascDcabachtopv{"cascDcabachtopv", 0.1, "cascDcabachtopv"};
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
    Configurable<float> cascv0cospa{"cascv0cospa", 0.98, "V0 CosPA"};
    Configurable<float> cascv0RadiusMin{"cascv0RadiusMin", 2.5, "v0radius"};
    Configurable<float> proplifetime{"proplifetime", 3, "ctau/<ctau>"};
    Configurable<float> lambdaMassWin{"lambdaMassWin", 0.005, "V0 Mass window limit"};
    Configurable<float> rejcomp{"rejcomp", 0.008, "Competing Cascade rejection"};
    Configurable<float> rapCut{"rapCut", 0.8, "Rapidity acceptance"};
  } MorePbPbsystCuts;
  // must include windows for background and peak
  Configurable<float> maxMassNSigma{"maxMassNSigma", 12.0f, "max mass region to be considered for further analysis"};

  // For extracting strangeness mass QA plots
  struct : ConfigurableGroup {
    ConfigurableAxis axisPtQA{"axisPtQA", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for QA histograms"};
    ConfigurableAxis axisK0ShortMass{"axisK0ShortMass", {200, 0.400f, 0.600f}, "Inv. Mass (GeV/c^{2})"};
    ConfigurableAxis axisLambdaMass{"axisLambdaMass", {200, 1.01f, 1.21f}, "Inv. Mass (GeV/c^{2})"};
    ConfigurableAxis axisXiMass{"axisXiMass", {200, 1.22f, 1.42f}, "Inv. Mass (GeV/c^{2})"};
    ConfigurableAxis axisOmegaMass{"axisOmegaMass", {200, 1.57f, 1.77f}, "Inv. Mass (GeV/c^{2})"};
    ConfigurableAxis axisMult{"axisMult", {VARIABLE_WIDTH, 0.0f, 0.01f, 1.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 70.0f, 100.0f}, "Centrality percentile bins"};
  } axesConfigurations;

  // QA
  Configurable<bool> doTrueSelectionInMass{"doTrueSelectionInMass", false, "Fill mass histograms only with true primary Particles for MC"};
  // Do declarative selections for DCAs, if possible
  Filter preFilterTracks = nabs(aod::track::dcaXY) < systCuts.dcaXYconstant + systCuts.dcaXYpTdep * nabs(aod::track::signed1Pt);
  Filter preFilterV0 = nabs(aod::v0data::dcapostopv) > systCuts.dcaPostopv&&
                                                         nabs(aod::v0data::dcanegtopv) > systCuts.dcaNegtopv&& aod::v0data::dcaV0daughters < systCuts.dcaV0dau;
  Filter preFilterCascade =
    nabs(aod::cascdata::dcapostopv) > systCuts.dcaPostopv&& nabs(aod::cascdata::dcanegtopv) > systCuts.dcaNegtopv&& nabs(aod::cascdata::dcabachtopv) > systCuts.dcaBachtopv&& aod::cascdata::dcaV0daughters < systCuts.dcaV0dau&& aod::cascdata::dcacascdaughters < systCuts.dcaCascdau;

  // using V0LinkedTagged = soa::Join<aod::V0sLinked, aod::V0Tags>;
  // using CascadesLinkedTagged = soa::Join<aod::CascadesLinked, aod::CascTags>;
  using FullTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA>;
  using FullTracksMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels>;
  using DauTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::TracksDCA>;
  using DauTracksMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::TracksDCA, aod::McTrackLabels>;
  // using IDTracks= soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidBayesPi, aod::pidBayesKa, aod::pidBayesPr, aod::TOFSignal>; // prepared for Bayesian PID
  using IDTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::pidTPCFullEl, aod::pidTOFFullEl, aod::TOFSignal, aod::TracksDCA>;
  using IDTracksMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::pidTPCFullEl, aod::pidTOFFullEl, aod::TOFSignal, aod::TracksDCA, aod::McTrackLabels>;
  using V0DatasWithoutTrackX = soa::Join<aod::V0Indices, aod::V0Cores>;
  using V0DatasWithoutTrackXMC = soa::Join<aod::V0Indices, aod::V0Cores, aod::V0MCCores>;
  using CascDatasMC = soa::Join<aod::CascDatas, aod::CascMCCores>;

  Produces<aod::TriggerTracks> triggerTrack;
  Produces<aod::TriggerTrackExtras> triggerTrackExtra;
  Produces<aod::AssocV0s> assocV0;
  Produces<aod::AssocCascades> assocCascades;
  Produces<aod::AssocHadrons> assocHadrons;
  Produces<aod::AssocPID> assocPID;

  Produces<aod::V0CoresBaseWithDua> v0cores;
  Produces<aod::DrvTracks> drvtracks;
  Produces<aod::BCandTime> bcAndTime;
  Produces<aod::DrvCollisions> newCollisions;

  struct : ConfigurableGroup {
    // invariant mass parametrizations
    Configurable<std::vector<float>> massParsK0Mean{"massParsK0Mean", {0.495967, 0.000095, 0.001120, 0.800000}, "pars in [0]+[1]*x+[2]*std::exp(-[3]*x)"};
    Configurable<std::vector<float>> massParsK0Width{"massParsK0Width", {0.002324, 0.000600, 0.005076, 1.644687}, "pars in [0]+[1]*x+[2]*std::exp(-[3]*x)"};

    Configurable<std::vector<float>> massParsLambdaMean{"massParsLambdaMean", {1.115554, 0.000002, -0.000311, 1.303969}, "pars in [0]+[1]*x+[2]*std::exp(-[3]*x)"};
    Configurable<std::vector<float>> massParsLambdaWidth{"massParsLambdaWidth", {0.001066, 0.000168, 0.001893, 1.407199}, "pars in [0]+[1]*x+[2]*std::exp(-[3]*x)"};

    Configurable<std::vector<float>> massParsCascadeMean{"massParsCascadeMean", {1.322150, -0.000087, -0.000761, 0.316391}, "pars in [0]+[1]*x+[2]*std::exp(-[3]*x)"};
    Configurable<std::vector<float>> massParsCascadeWidth{"massParsCascadeWidth", {0.001269, 0.000249, 0.002790, 1.128544}, "pars in [0]+[1]*x+[2]*std::exp(-[3]*x)"};

    Configurable<std::vector<float>> massParsOmegaMean{"massParsOmegaMean", {1.671908, -0.000000, 0.001027, 2.263832}, "pars in [0]+[1]*x+[2]*std::exp(-[3]*x)"};
    Configurable<std::vector<float>> massParsOmegaWidth{"massParsOmegaWidth", {0.001223, 0.000300, 0.040718, 2.826750}, "pars in [0]+[1]*x+[2]*std::exp(-[3]*x)"};
  } parameters;
  TF1* fK0Mean = new TF1("fK0Mean", "[0]+[1]*x+[2]*std::exp(-[3]*x)");
  TF1* fK0Width = new TF1("fK0Width", "[0]+[1]*x+[2]*std::exp(-[3]*x)");
  TF1* fLambdaMean = new TF1("fLambdaMean", "[0]+[1]*x+[2]*std::exp(-[3]*x)");
  TF1* fLambdaWidth = new TF1("fLambdaWidth", "[0]+[1]*x+[2]*std::exp(-[3]*x)");
  TF1* fXiMean = new TF1("fXiMean", "[0]+[1]*x+[2]*std::exp(-[3]*x)");
  TF1* fXiWidth = new TF1("fXiWidth", "[0]+[1]*x+[2]*std::exp(-[3]*x)");
  TF1* fOmegaMean = new TF1("fomegaMean", "[0]+[1]*x+[2]*std::exp(-[3]*x)");
  TF1* fOmegaWidth = new TF1("fomegaWidth", "[0]+[1]*x+[2]*std::exp(-[3]*x)");
  void init(InitContext const&)
  {
    if (useParameterization) {
      fK0Mean->SetParameters(parameters.massParsK0Mean->at(0), parameters.massParsK0Mean->at(1), parameters.massParsK0Mean->at(2), parameters.massParsK0Mean->at(3));
      fK0Width->SetParameters(parameters.massParsK0Width->at(0), parameters.massParsK0Width->at(1), parameters.massParsK0Width->at(2), parameters.massParsK0Width->at(3));
      fLambdaMean->SetParameters(parameters.massParsLambdaMean->at(0), parameters.massParsLambdaMean->at(1), parameters.massParsLambdaMean->at(2), parameters.massParsLambdaMean->at(3));
      fLambdaWidth->SetParameters(parameters.massParsLambdaWidth->at(0), parameters.massParsLambdaWidth->at(1), parameters.massParsLambdaWidth->at(2), parameters.massParsLambdaWidth->at(3));
      fXiMean->SetParameters(parameters.massParsCascadeMean->at(0), parameters.massParsCascadeMean->at(1), parameters.massParsCascadeMean->at(2), parameters.massParsCascadeMean->at(3));
      fXiWidth->SetParameters(parameters.massParsCascadeWidth->at(0), parameters.massParsCascadeWidth->at(1), parameters.massParsCascadeWidth->at(2), parameters.massParsCascadeWidth->at(3));
      fOmegaMean->SetParameters(parameters.massParsOmegaMean->at(0), parameters.massParsOmegaMean->at(1), parameters.massParsOmegaMean->at(2), parameters.massParsOmegaMean->at(3));
      fOmegaWidth->SetParameters(parameters.massParsOmegaWidth->at(0), parameters.massParsOmegaWidth->at(1), parameters.massParsOmegaWidth->at(2), parameters.massParsOmegaWidth->at(3));
    }
  }
  template <class TTrack>
  bool isValidTrigger(TTrack track)
  {
    if (track.eta() > systCuts.triggerEtaMax || track.eta() < systCuts.triggerEtaMin) {
      return false;
    }
    if (track.pt() > systCuts.triggerPtCutMax || track.pt() < systCuts.triggerPtCutMin) {
      return false;
    }
    if (track.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRows) {
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
    return true;
  }
  template <class TTracks>
  bool hasValidTrigger(TTracks tracks)
  {
    for (auto const& track : tracks) {
      if (!isValidTrigger(track))
        continue;
      return true;
    }
    return false;
  }
  template <typename TCollision>
  bool isCollisionSelectedPbPb(TCollision collision)
  {
    if (!collision.selection_bit(aod::evsel::kIsTriggerTVX)) /* FT0 vertex (acceptable FT0C-FT0A time difference) collisions */
      return false;
    if (!collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) // cut time intervals with dead ITS staves
      return false;
    if (!collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) // removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference
      return false;
    auto occupancy = collision.trackOccupancyInTimeRange();
    if (occupancy < cfgCutOccupancyLow || occupancy > cfgCutOccupancyHigh) /* Below min occupancy and Above max occupancy*/
      return false;
    if (!collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) // reject collisions close to Time Frame borders
      return false;
    if (!collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) // reject events affected by the ITS ROF border
      return false;
    if (!collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) // rejects collisions which are associated with the same "found-by-T0" bunch crossing
      return false;
    return true;
  }
  template <typename TCollision, typename TV0, typename DauTracks = DauTracks>
  bool isValidV0(TCollision collision, TV0 const& v0)
  {
    if (v0.v0radius() < systCuts.v0RadiusMin || v0.v0radius() > systCuts.v0RadiusMax || v0.eta() > systCuts.assocEtaMax || v0.eta() < systCuts.assocEtaMin || v0.v0cosPA() < systCuts.v0Cospa) {
      return false;
    }
    if (v0.pt() > systCuts.assocPtCutMax || v0.pt() < systCuts.assocPtCutMin) {
      return false;
    }
    // check dE/dx compatibility
    int compatibleK0Short = 0;
    int compatibleLambda = 0;
    int compatibleAntiLambda = 0;

    auto posdau = v0.template posTrack_as<DauTracks>();
    auto negdau = v0.template negTrack_as<DauTracks>();

    if (negdau.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRows)
      return false;
    if (posdau.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRows)
      return false;
    if (v0.dcaV0daughters() > systCuts.dcaV0dau)
      return false;
    if (std::abs(posdau.tpcNSigmaPi()) < strangedEdxNSigmaLoose && std::abs(negdau.tpcNSigmaPi()) < strangedEdxNSigmaLoose) {
      if (doPPAnalysis || (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short < systCuts.lifetimecutK0S &&
                           std::abs(v0.dcapostopv()) > systCuts.dcapostopvK0S && std::abs(v0.dcanegtopv()) > systCuts.dcanegtopvK0S &&
                           v0.qtarm() * systCuts.armPodCut > std::abs(v0.alpha()))) {
        BIT_SET(compatibleK0Short, 0);
      }
    }

    if (std::abs(posdau.tpcNSigmaPr()) < strangedEdxNSigmaLoose && std::abs(negdau.tpcNSigmaPi()) < strangedEdxNSigmaLoose) {
      if (v0.v0cosPA() > systCuts.lambdaCospa) {
        if (doPPAnalysis || (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short < systCuts.lifetimecutLambda &&
                             std::abs(v0.dcapostopv()) > systCuts.dcapostopvLambda && std::abs(v0.dcanegtopv()) > systCuts.dcanegtopvLambda)) {
          BIT_SET(compatibleLambda, 0);
        }
      }
    }

    if (std::abs(posdau.tpcNSigmaPi()) < strangedEdxNSigmaLoose && std::abs(negdau.tpcNSigmaPr()) < strangedEdxNSigmaLoose) {
      if (v0.v0cosPA() > systCuts.lambdaCospa) {
        if (doPPAnalysis || (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short < systCuts.lifetimecutLambda &&
                             std::abs(v0.dcapostopv()) > systCuts.dcapostopvAntiLambda && std::abs(v0.dcanegtopv()) > systCuts.dcanegtopvAntiLambda)) {
          BIT_SET(compatibleAntiLambda, 0);
        }
      }
    }
    float massNSigmaK0Short = -20.0f;
    float massNSigmaLambda = -20.0f;
    float massNSigmaAntiLambda = -20.0f;

    massNSigmaK0Short = (v0.mK0Short() - fK0Mean->Eval(v0.pt())) / (fK0Width->Eval(v0.pt()) + 1e-6);
    massNSigmaLambda = (v0.mLambda() - fLambdaMean->Eval(v0.pt())) / (fLambdaWidth->Eval(v0.pt()) + 1e-6);
    massNSigmaAntiLambda = (v0.mAntiLambda() - fLambdaMean->Eval(v0.pt())) / (fLambdaWidth->Eval(v0.pt()) + 1e-6);

    if ((compatibleK0Short > 0 && std::abs(massNSigmaK0Short) < maxMassNSigma) ||
        (compatibleLambda > 0 && std::abs(massNSigmaLambda) < maxMassNSigma) ||
        (compatibleAntiLambda > 0 && std::abs(massNSigmaAntiLambda) < maxMassNSigma))
      return true;

    return false;
  }
  template <typename TCollision, typename TV0s, typename DauTracks = DauTracks>
  bool hasValidV0(TCollision collision, TV0s const& V0s)
  {
    for (auto const& v0 : V0s) {
      if (isValidV0(collision, v0)) {
        return true;
      }
    }
    return false;
  }
  void processProduceDerivedData(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::PVMults>::iterator const& collision, DauTracks const& tracks, soa::Filtered<V0DatasWithoutTrackX> const& v0s, aod::BCsWithTimestamps const&)
  {
    if (!isCollisionSelectedPbPb(collision))
      return;
    if (!hasValidTrigger(tracks))
      return;
    if (!hasValidV0(collision, v0s))
      return;

    newCollisions(collision.bcId(),
                  collision.posX(), collision.posY(), collision.posZ(),
                  collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ(),
                  collision.centFT0C(), collision.centFT0M(),
                  collision.multNTracksPVeta1(),
                  collision.isInelGt0(), collision.isInelGt1(),
                  collision.flags(),
                  collision.chi2(),
                  collision.numContrib(),
                  collision.collisionTime(),
                  collision.collisionTimeRes());
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    bcAndTime(
      bc.timestamp(),
      bc.runNumber(),
      bc.globalBC());

    std::unordered_set<int64_t> selectedTrackIndexSet;

    for (auto const& track : tracks) {
      if (!isValidTrigger(track))
        continue;
      selectedTrackIndexSet.insert(track.globalIndex());
      drvtracks(
        track.collisionId(),
        track.pt(),
        track.eta(),
        track.phi(),
        track.sign(),
        track.signed1Pt(),
        track.px(),
        track.py(),
        track.pz(),
        track.tpcNClsCrossedRows(),
        track.hasITS(),
        track.tpcNClsShared(),
        track.itsClusterMap(),
        track.dcaXY(),
        track.tpcNSigmaPi(),
        track.tpcNSigmaPr());
    }
    for (auto const& v0 : v0s) {
      if (!isValidV0(collision, v0))
        continue;
      auto posdau = v0.posTrack_as<DauTracks>();
      auto negdau = v0.negTrack_as<DauTracks>();

      if (selectedTrackIndexSet.count(posdau.globalIndex()) > 0 || selectedTrackIndexSet.count(negdau.globalIndex()) > 0) {
        continue;
      }

      v0cores(v0.x(), v0.y(), v0.z(),
              v0.pxpos(), v0.pypos(), v0.pzpos(),
              v0.pxneg(), v0.pyneg(), v0.pzneg(),
              v0.dcaV0daughters(),
              v0.dcapostopv(),
              v0.dcanegtopv(),
              v0.v0cosPA(),
              v0.dcav0topv(),
              v0.v0Type(),
              v0.posTrackId(),
              v0.negTrackId(),
              v0.collisionId());
      // drvtracks(
      //   posdau.collisionId(),
      //   posdau.pt(),
      //   posdau.eta(),
      //   posdau.phi(),
      //   posdau.sign(),
      //   posdau.signed1Pt(),
      //   posdau.px(),
      //   posdau.py(),
      //   posdau.pz(),
      //   posdau.tpcNClsCrossedRows(),
      //   posdau.hasITS(),
      //   posdau.tpcNClsShared(),
      //   posdau.itsClusterMap(),
      //   posdau.dcaXY(),
      //   posdau.tpcNSigmaPi(),
      //   posdau.tpcNSigmaPr());
      // drvtracks(
      //   negdau.collisionId(),
      //   negdau.pt(),
      //   negdau.eta(),
      //   negdau.phi(),
      //   negdau.sign(),
      //   negdau.signed1Pt(),
      //   negdau.px(),
      //   negdau.py(),
      //   negdau.pz(),
      //   negdau.tpcNClsCrossedRows(),
      //   negdau.hasITS(),
      //   negdau.tpcNClsShared(),
      //   negdau.itsClusterMap(),
      //   negdau.dcaXY(),
      //   negdau.tpcNSigmaPi(),
      //   negdau.tpcNSigmaPr());
    }
  }
  PROCESS_SWITCH(HStrangeDerivedData, processProduceDerivedData, "Produce DerivedData", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HStrangeDerivedData>(cfgc)};
}
