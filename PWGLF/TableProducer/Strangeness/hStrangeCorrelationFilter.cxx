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
/// \brief This task pre-filters tracks, V0s and cascades to do h-strangeness
///        correlations with an analysis task.
///
/// \file hStrangeCorrelationFilter.cxx
/// \author Kai Cui (kaicui@mails.ccnu.edu.cn)
/// \author Lucia Anna Tarasovicova (lucia.anna.husova@cern.ch)
/// \author David Dobrigkeit Chinellato (david.dobrigkeit.chinellato@cern.ch)
/// \author Zhongbao Yin (Zhong-Bao.Yin@cern.ch)

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

#include "CCDB/BasicCCDBManager.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "TF1.h"
#include <TPDGCode.h>

#include <string>

using namespace o2;
using namespace o2::constants::math;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define BIT_SET(var, nbit) ((var) |= (1 << (nbit)))
#define BIT_CHECK(var, nbit) ((var) & (1 << (nbit)))

struct HStrangeCorrelationFilter {
  const float ctauxiPDG = 4.91;     // from PDG
  const float ctauomegaPDG = 2.461; // from PDG

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // master analysis switches
  Configurable<bool> doPPAnalysis{"doPPAnalysis", true, "if in pp, set to true"};
  Configurable<bool> useParameterization{"useParameterization", true, "ture for parameterization method, false for hist method"};
  // Operational
  Configurable<bool> fillTableOnlyWithCompatible{"fillTableOnlyWithCompatible", true, "pre-apply dE/dx, broad mass window in table filling"};
  Configurable<float> strangedEdxNSigmaLoose{"strangedEdxNSigmaLoose", 5, "Nsigmas for strange decay daughters"};
  Configurable<float> strangedEdxNSigma{"strangedEdxNSigma", 4, "Nsigmas for strange decay daughters"};
  Configurable<float> strangedEdxNSigmaTight{"strangedEdxNSigmaTight", 3, "Nsigmas for strange decay daughters"};
  Configurable<std::string> zorroMask{"zorroMask", "", "zorro trigger class to select on (empty: none)"};
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
    Configurable<float> triggerPtCutMin{"triggerPtCutMin", 3, "triggerptmin"};
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
  Configurable<std::string> ccdburl{"ccdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository to use"};
  Configurable<std::string> parameterCCDBPath{"parameterCCDBPath", "Users/k/kcui/LHC25b4a/parameter", "Path of the mean and sigma"};

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
  TH1F* hK0ShortMean;
  TH1F* hK0ShortWidth;
  TH1F* hLambdaMean;
  TH1F* hLambdaWidth;
  TH1F* hXiMean;
  TH1F* hXiWidth;
  TH1F* hOmegaMean;
  TH1F* hOmegaWidth;
  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};
  int mRunNumber;

  void init(InitContext const&)
  {
    zorroSummary.setObject(zorro.getZorroSummary());
    mRunNumber = -1;
    if (useParameterization) {
      fK0Mean->SetParameters(parameters.massParsK0Mean->at(0), parameters.massParsK0Mean->at(1), parameters.massParsK0Mean->at(2), parameters.massParsK0Mean->at(3));
      fK0Width->SetParameters(parameters.massParsK0Width->at(0), parameters.massParsK0Width->at(1), parameters.massParsK0Width->at(2), parameters.massParsK0Width->at(3));
      fLambdaMean->SetParameters(parameters.massParsLambdaMean->at(0), parameters.massParsLambdaMean->at(1), parameters.massParsLambdaMean->at(2), parameters.massParsLambdaMean->at(3));
      fLambdaWidth->SetParameters(parameters.massParsLambdaWidth->at(0), parameters.massParsLambdaWidth->at(1), parameters.massParsLambdaWidth->at(2), parameters.massParsLambdaWidth->at(3));
      fXiMean->SetParameters(parameters.massParsCascadeMean->at(0), parameters.massParsCascadeMean->at(1), parameters.massParsCascadeMean->at(2), parameters.massParsCascadeMean->at(3));
      fXiWidth->SetParameters(parameters.massParsCascadeWidth->at(0), parameters.massParsCascadeWidth->at(1), parameters.massParsCascadeWidth->at(2), parameters.massParsCascadeWidth->at(3));
      fOmegaMean->SetParameters(parameters.massParsOmegaMean->at(0), parameters.massParsOmegaMean->at(1), parameters.massParsOmegaMean->at(2), parameters.massParsOmegaMean->at(3));
      fOmegaWidth->SetParameters(parameters.massParsOmegaWidth->at(0), parameters.massParsOmegaWidth->at(1), parameters.massParsOmegaWidth->at(2), parameters.massParsOmegaWidth->at(3));
    } else {
      hK0ShortMean = 0x0;
      hK0ShortWidth = 0x0;
      hLambdaMean = 0x0;
      hLambdaWidth = 0x0;
      hXiMean = 0x0;
      hXiWidth = 0x0;
      hOmegaMean = 0x0;
      hOmegaWidth = 0x0;
    }
    if (doprocessV0s || doprocessV0sMC) {
      histos.add("h3dMassK0Short", "h3dMassK0Short", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisK0ShortMass, axesConfigurations.axisMult});
      histos.add("h3dMassLambda", "h3dMassLambda", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisLambdaMass, axesConfigurations.axisMult});
      histos.add("h3dMassAntiLambda", "h3dMassAntiLambda", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisLambdaMass, axesConfigurations.axisMult});
    }
    if (doprocessCascades || doprocessCascadesMC) {
      histos.add("h3dMassXiMinus", "h3dMassXiMinus", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisXiMass, axesConfigurations.axisMult});
      histos.add("h3dMassXiPlus", "h3dMassXiPlus", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisXiMass, axesConfigurations.axisMult});
      histos.add("h3dMassOmegaMinus", "h3dMassOmegaMinus", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisOmegaMass, axesConfigurations.axisMult});
      histos.add("h3dMassOmegaPlus", "h3dMassOmegaPlus", kTH3F, {axesConfigurations.axisPtQA, axesConfigurations.axisOmegaMass, axesConfigurations.axisMult});
    }
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), zorroMask.value);
    zorro.populateHistRegistry(histos, bc.runNumber());

    mRunNumber = bc.runNumber();
  }

  void initParametersFromCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    mRunNumber = bc.runNumber();
    LOG(info) << "Loading mean and sigma from CCDB for run " << mRunNumber << " now...";
    auto timeStamp = bc.timestamp();

    TList* listParameters = ccdb->getForTimeStamp<TList>(parameterCCDBPath, timeStamp);

    if (!listParameters) {
      LOG(fatal) << "Problem getting TList object with parameters!";
    }
    if (doprocessV0s || doprocessV0sMC) {
      hK0ShortMean = static_cast<TH1F*>(listParameters->FindObject("hK0ShortMean"));
      hK0ShortWidth = static_cast<TH1F*>(listParameters->FindObject("hK0ShortWidth"));
      hLambdaMean = static_cast<TH1F*>(listParameters->FindObject("hLambdaMean"));
      hLambdaWidth = static_cast<TH1F*>(listParameters->FindObject("hLambdaWidth"));
    }
    if (doprocessCascades || doprocessCascadesMC) {
      hXiMean = static_cast<TH1F*>(listParameters->FindObject("hXiMean"));
      hXiWidth = static_cast<TH1F*>(listParameters->FindObject("hXiWidth"));
      hOmegaMean = static_cast<TH1F*>(listParameters->FindObject("hOmegaMean"));
      hOmegaWidth = static_cast<TH1F*>(listParameters->FindObject("hOmegaWidth"));
    }
    LOG(info) << "parameters now loaded for " << mRunNumber;
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
    if (std::abs(collision.posZ()) > eventSelections.zVertexCut) {
      return false;
    }
    if (collision.centFT0M() > 100 || collision.centFT0M() < 0) {
      return false;
    }
    if (!collision.isInelGt0() && eventSelections.selectINELgtZERO) {
      return false;
    }
    if (!collision.selection_bit(aod::evsel::kIsGoodITSLayersAll) && eventSelections.requireAllGoodITSLayers) {
      return false;
    }
    if (zorroMask.value != "") {
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      bool zorroSelected = zorro.isSelected(collision.template bc_as<aod::BCsWithTimestamps>().globalBC()); /// Just let Zorro do the accounting
      if (!zorroSelected) {
        return false;
      }
    }
    return true;
  }

  // more event selections in Pb-Pb
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

  // reco-level trigger quality checks (N.B.: DCA is filtered, not selected)
  template <class TTrack>
  bool isValidTrigger(TTrack track)
  {
    if (track.eta() > systCuts.triggerEtaMax || track.eta() < systCuts.triggerEtaMin) {
      return false;
    }
    // if (track.sign()= 1 ) {continue;}
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
    if (!(BIT_CHECK(track.itsClusterMap(), 0)) && systCuts.triggerRequireL0) {
      return false; // skip, doesn't have cluster in ITS L0
    }
    return true;
  }
  template <class TTrack>
  bool isValidAssocTrack(TTrack assoc)
  {
    if (assoc.eta() > systCuts.assocEtaMax || assoc.eta() < systCuts.assocEtaMin) {
      return false;
    }
    if (assoc.pt() > systCuts.assocPtCutMax || assoc.pt() < systCuts.assocPtCutMin) {
      return false;
    }
    if (assoc.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRows) {
      return false; // crossed rows
    }
    if (!assoc.hasITS() && systCuts.assocRequireITS) {
      return false; // skip, doesn't have ITS signal (skips lots of TPC-only!)
    }

    // do this only if information is available
    float nSigmaTPCTOF[8] = {-10, -10, -10, -10, -10, -10, -10, -10};
    if constexpr (requires { assoc.tofSignal(); } && !requires { assoc.mcParticle(); }) {
      if (assoc.tofSignal() > 0) {
        if (std::sqrt(assoc.tofNSigmaPi() * assoc.tofNSigmaPi() + assoc.tpcNSigmaPi() * assoc.tpcNSigmaPi()) > systCuts.assocPionNSigmaTPCFOF)
          return false;
        if (assoc.tofNSigmaPr() < systCuts.rejectSigma)
          return false;
        if (assoc.tpcNSigmaPr() < systCuts.rejectSigma)
          return false;
        if (assoc.tofNSigmaKa() < systCuts.rejectSigma)
          return false;
        if (assoc.tpcNSigmaKa() < systCuts.rejectSigma)
          return false;
        nSigmaTPCTOF[4] = assoc.tofNSigmaPi();
        nSigmaTPCTOF[5] = assoc.tofNSigmaKa();
        nSigmaTPCTOF[6] = assoc.tofNSigmaPr();
        nSigmaTPCTOF[7] = assoc.tofNSigmaEl();
      } else {
        if (assoc.tpcNSigmaPi() > systCuts.assocPionNSigmaTPCFOF)
          return false;
        if (assoc.tpcNSigmaPr() < systCuts.rejectSigma)
          return false;
        if (assoc.tpcNSigmaKa() < systCuts.rejectSigma)
          return false;
      }
      nSigmaTPCTOF[0] = assoc.tpcNSigmaPi();
      nSigmaTPCTOF[1] = assoc.tpcNSigmaKa();
      nSigmaTPCTOF[2] = assoc.tpcNSigmaPr();
      nSigmaTPCTOF[3] = assoc.tpcNSigmaEl();
    }

    bool physicalPrimary = false;
    float origPt = -1;
    float pdgCode = -9999;
    if constexpr (requires { assoc.mcParticle(); }) {
      if (assoc.has_mcParticle()) {
        auto mcParticle = assoc.mcParticle();
        physicalPrimary = mcParticle.isPhysicalPrimary();
        origPt = mcParticle.pt();
        pdgCode = mcParticle.pdgCode();
      }
    }

    assocHadrons(
      assoc.collisionId(),
      physicalPrimary,
      assoc.globalIndex(),
      origPt,
      pdgCode);
    assocPID(
      nSigmaTPCTOF[0],
      nSigmaTPCTOF[1],
      nSigmaTPCTOF[2],
      nSigmaTPCTOF[3],
      nSigmaTPCTOF[4],
      nSigmaTPCTOF[5],
      nSigmaTPCTOF[6],
      nSigmaTPCTOF[7]);
    return true;
  }

  // cascadeselection in PbPb
  template <typename TCascade>
  bool CascadeSelectedPbPb(TCascade casc, float pvx, float pvy, float pvz)
  {
    // bachBaryonCosPA
    if (casc.bachBaryonCosPA() < MorePbPbsystCuts.bachBaryonCosPA)
      return false;
    // bachBaryonDCAxyToPV
    if (std::abs(casc.bachBaryonDCAxyToPV()) > MorePbPbsystCuts.bachBaryonDCAxyToPV)
      return false;
    // casccosPA
    if (casc.casccosPA(pvx, pvy, pvz) < systCuts.cascCospa)
      return false;
    // dcacascdaughters
    float ptDepCut = MorePbPbsystCuts.dcaCacsDauPar0;
    if (casc.pt() > 1 && casc.pt() < 4)
      ptDepCut = MorePbPbsystCuts.dcaCacsDauPar1;
    else if (casc.pt() > 4)
      ptDepCut = MorePbPbsystCuts.dcaCacsDauPar2;
    if (casc.dcacascdaughters() > ptDepCut)
      return false;
    // dcaV0daughters
    if (casc.dcaV0daughters() > systCuts.dcaV0dau)
      return false;
    // dcav0topv
    if (std::abs(casc.dcav0topv(pvx, pvy, pvz)) < MorePbPbsystCuts.cascdcaV0ToPV)
      return false;
    // cascradius
    if (casc.cascradius() < systCuts.cascRadius)
      return false;
    // v0radius
    if (casc.v0radius() < MorePbPbsystCuts.cascv0RadiusMin)
      return false;
    // v0cosPA
    if (casc.v0cosPA(casc.x(), casc.y(), casc.z()) < MorePbPbsystCuts.cascv0cospa)
      return false;
    // lambdaMassWin
    if (std::abs(casc.mLambda() - o2::constants::physics::MassLambda0) > MorePbPbsystCuts.lambdaMassWin)
      return false;
    return true;
  }

  // for real data processing
  void processTriggers(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::PVMults>::iterator const& collision, soa::Filtered<FullTracks> const& tracks, aod::BCsWithTimestamps const&)
  {
    if (((doPPAnalysis && !isCollisionSelected(collision))) || (!doPPAnalysis && !isCollisionSelectedPbPb(collision))) {
      return;
    }

    /// _________________________________________________
    /// Step 1: Populate table with trigger tracks
    for (auto const& track : tracks) {
      if (!isValidTrigger(track))
        continue;
      triggerTrack(
        track.collisionId(),
        false, // if you decide to check real data for primaries, you'll have a hard time
        track.globalIndex(),
        0);
      triggerTrackExtra(1);
    }
  }

  // for MC processing
  void processTriggersMC(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::PVMults>::iterator const& collision, soa::Filtered<FullTracksMC> const& tracks, aod::McParticles const&, aod::BCsWithTimestamps const&)
  {
    if (((doPPAnalysis && !isCollisionSelected(collision))) || (!doPPAnalysis && !isCollisionSelectedPbPb(collision))) {
      return;
    }

    /// _________________________________________________
    /// Step 1: Populate table with trigger tracks
    for (auto const& track : tracks) {
      if (!isValidTrigger(track))
        continue;
      bool physicalPrimary = false;
      float origPt = -1;
      if (track.has_mcParticle()) {
        auto mcParticle = track.mcParticle();
        physicalPrimary = mcParticle.isPhysicalPrimary();
        origPt = mcParticle.pt();
      }
      triggerTrack(
        track.collisionId(),
        physicalPrimary,
        track.globalIndex(),
        origPt);
      triggerTrackExtra(1);
    }
  }

  void processAssocPions(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<IDTracks> const& tracks, aod::BCsWithTimestamps const&)
  {
    // Load parameters for sideband subtraction
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    // Perform basic event selection
    if (!collision.sel8()) {
      return;
    }
    // No need to correlate stuff that's in far collisions
    if (std::abs(collision.posZ()) > 10.0) {
      return;
    }
    if (zorroMask.value != "") {
      initCCDB(bc);
      bool zorroSelected = zorro.isSelected(collision.bc_as<aod::BCsWithTimestamps>().globalBC()); /// Just let Zorro do the accounting
      if (!zorroSelected) {
        return;
      }
    }

    /// _________________________________________________
    /// Step 1: Populate table with trigger tracks
    for (auto const& track : tracks) {
      if (!isValidAssocTrack(track))
        continue;
    }
  }

  void processAssocPionsMC(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<IDTracksMC> const& tracks, aod::McParticles const&, aod::BCsWithTimestamps const&)
  {
    // Load parameters for sideband subtraction
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    // Perform basic event selection
    if (!collision.sel8()) {
      return;
    }
    // No need to correlate stuff that's in far collisions
    if (std::abs(collision.posZ()) > 10.0) {
      return;
    }
    if (zorroMask.value != "") {
      initCCDB(bc);
      bool zorroSelected = zorro.isSelected(collision.bc_as<aod::BCsWithTimestamps>().globalBC()); /// Just let Zorro do the accounting
      if (!zorroSelected) {
        return;
      }
    }

    /// _________________________________________________
    /// Step 1: Populate table with trigger tracks
    for (auto const& track : tracks) {
      if (!isValidAssocTrack(track))
        continue;
    }
  }

  void processAssocHadrons(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<FullTracks> const& tracks, aod::BCsWithTimestamps const&)
  {
    // Load parameters for sideband subtraction
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    // Perform basic event selection
    if (!collision.sel8()) {
      return;
    }
    // No need to correlate stuff that's in far collisions
    if (std::abs(collision.posZ()) > 10.0) {
      return;
    }
    if (zorroMask.value != "") {
      initCCDB(bc);
      bool zorroSelected = zorro.isSelected(collision.bc_as<aod::BCsWithTimestamps>().globalBC()); /// Just let Zorro do the accounting
      if (!zorroSelected) {
        return;
      }
    }

    /// _________________________________________________
    /// Step 1: Populate table with trigger tracks
    for (auto const& track : tracks) {
      if (!isValidAssocTrack(track))
        continue;
    }
  }
  void processAssocHadronsMC(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<FullTracksMC> const& tracks, aod::McParticles const&, aod::BCsWithTimestamps const&)
  {
    // Load parameters for sideband subtraction
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    // Perform basic event selection
    if (!collision.sel8()) {
      return;
    }
    // No need to correlate stuff that's in far collisions
    if (std::abs(collision.posZ()) > 10.0) {
      return;
    }
    if (zorroMask.value != "") {
      initCCDB(bc);
      bool zorroSelected = zorro.isSelected(collision.bc_as<aod::BCsWithTimestamps>().globalBC()); /// Just let Zorro do the accounting
      if (!zorroSelected) {
        return;
      }
    }

    /// _________________________________________________
    /// Step 1: Populate table with trigger tracks
    for (auto const& track : tracks) {
      if (!isValidAssocTrack(track))
        continue;
    }
  }

  void processV0s(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::PVMults>::iterator const& collision, DauTracks const&, soa::Filtered<V0DatasWithoutTrackX> const& V0s, aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    double cent = doPPAnalysis ? collision.centFT0M() : collision.centFT0C();
    if (((doPPAnalysis && !isCollisionSelected(collision))) || (!doPPAnalysis && !isCollisionSelectedPbPb(collision))) {
      return;
    }
    /// _________________________________________________
    /// Populate table with associated V0s
    for (auto const& v0 : V0s) {
      if (v0.v0radius() < systCuts.v0RadiusMin || v0.v0radius() > systCuts.v0RadiusMax || v0.eta() > systCuts.assocEtaMax || v0.eta() < systCuts.assocEtaMin || v0.v0cosPA() < systCuts.v0Cospa) {
        continue;
      }
      if (v0.pt() > systCuts.assocPtCutMax || v0.pt() < systCuts.assocPtCutMin) {
        continue;
      }
      // check dE/dx compatibility
      int compatibleK0Short = 0;
      int compatibleLambda = 0;
      int compatibleAntiLambda = 0;

      auto posdau = v0.posTrack_as<DauTracks>();
      auto negdau = v0.negTrack_as<DauTracks>();

      if (negdau.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRows)
        continue;
      if (posdau.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRows)
        continue;

      if (std::abs(posdau.tpcNSigmaPi()) < strangedEdxNSigmaLoose && std::abs(negdau.tpcNSigmaPi()) < strangedEdxNSigmaLoose) {
        if (doPPAnalysis || (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short < systCuts.lifetimecutK0S &&
                             std::abs(v0.dcapostopv()) > systCuts.dcapostopvK0S && std::abs(v0.dcanegtopv()) > systCuts.dcanegtopvK0S &&
                             v0.qtarm() * systCuts.armPodCut > std::abs(v0.alpha()))) {
          BIT_SET(compatibleK0Short, 0);
        }
      }
      if (std::abs(posdau.tpcNSigmaPi()) < strangedEdxNSigma && std::abs(negdau.tpcNSigmaPi()) < strangedEdxNSigma) {
        if (doPPAnalysis || (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short < systCuts.lifetimecutK0S &&
                             std::abs(v0.dcapostopv()) > systCuts.dcapostopvK0S && std::abs(v0.dcanegtopv()) > systCuts.dcanegtopvK0S &&
                             v0.qtarm() * systCuts.armPodCut > std::abs(v0.alpha()))) {
          BIT_SET(compatibleK0Short, 1);
        }
      }
      if (std::abs(posdau.tpcNSigmaPi()) < strangedEdxNSigmaTight && std::abs(negdau.tpcNSigmaPi()) < strangedEdxNSigmaTight) {
        if (doPPAnalysis || (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short < systCuts.lifetimecutK0S &&
                             std::abs(v0.dcapostopv()) > systCuts.dcapostopvK0S && std::abs(v0.dcanegtopv()) > systCuts.dcanegtopvK0S &&
                             v0.qtarm() * systCuts.armPodCut > std::abs(v0.alpha()))) {
          BIT_SET(compatibleK0Short, 2);
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
      if (std::abs(posdau.tpcNSigmaPr()) < strangedEdxNSigma && std::abs(negdau.tpcNSigmaPi()) < strangedEdxNSigma) {
        if (v0.v0cosPA() > systCuts.lambdaCospa) {
          if (doPPAnalysis || (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short < systCuts.lifetimecutLambda &&
                               std::abs(v0.dcapostopv()) > systCuts.dcapostopvLambda && std::abs(v0.dcanegtopv()) > systCuts.dcanegtopvLambda)) {
            BIT_SET(compatibleLambda, 1);
          }
        }
      }
      if (std::abs(posdau.tpcNSigmaPr()) < strangedEdxNSigmaTight && std::abs(negdau.tpcNSigmaPi()) < strangedEdxNSigmaTight) {
        if (v0.v0cosPA() > systCuts.lambdaCospa) {
          if (doPPAnalysis || (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short < systCuts.lifetimecutLambda &&
                               std::abs(v0.dcapostopv()) > systCuts.dcapostopvLambda && std::abs(v0.dcanegtopv()) > systCuts.dcanegtopvLambda)) {
            BIT_SET(compatibleLambda, 2);
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
      if (std::abs(posdau.tpcNSigmaPi()) < strangedEdxNSigma && std::abs(negdau.tpcNSigmaPr()) < strangedEdxNSigma) {
        if (v0.v0cosPA() > systCuts.lambdaCospa) {
          if (doPPAnalysis || (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short < systCuts.lifetimecutLambda &&
                               std::abs(v0.dcapostopv()) > systCuts.dcapostopvAntiLambda && std::abs(v0.dcanegtopv()) > systCuts.dcanegtopvAntiLambda)) {
            BIT_SET(compatibleAntiLambda, 1);
          }
        }
      }
      if (std::abs(posdau.tpcNSigmaPi()) < strangedEdxNSigmaTight && std::abs(negdau.tpcNSigmaPr()) < strangedEdxNSigmaTight) {
        if (v0.v0cosPA() > systCuts.lambdaCospa) {
          if (doPPAnalysis || (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short < systCuts.lifetimecutLambda &&
                               std::abs(v0.dcapostopv()) > systCuts.dcapostopvAntiLambda && std::abs(v0.dcanegtopv()) > systCuts.dcanegtopvAntiLambda)) {
            BIT_SET(compatibleAntiLambda, 2);
          }
        }
      }
      float massNSigmaK0Short = -20.0f;
      float massNSigmaLambda = -20.0f;
      float massNSigmaAntiLambda = -20.0f;
      if (useParameterization) {
        massNSigmaK0Short = (v0.mK0Short() - fK0Mean->Eval(v0.pt())) / (fK0Width->Eval(v0.pt()) + 1e-6);
        massNSigmaLambda = (v0.mLambda() - fLambdaMean->Eval(v0.pt())) / (fLambdaWidth->Eval(v0.pt()) + 1e-6);
        massNSigmaAntiLambda = (v0.mAntiLambda() - fLambdaMean->Eval(v0.pt())) / (fLambdaWidth->Eval(v0.pt()) + 1e-6);
      } else {
        // Load parameters for sideband subtraction
        initParametersFromCCDB(bc);
        // simplified handling: calculate NSigma in mass here
        if (v0.pt() < 0.2f || v0.pt() > 14.5f) {
          massNSigmaK0Short = (v0.mK0Short() - hK0ShortMean->GetBinContent(hK0ShortMean->FindBin(v0.pt()))) / (hK0ShortWidth->GetBinContent(hK0ShortWidth->FindBin(v0.pt())) + 1e-6);
          massNSigmaLambda = (v0.mLambda() - hLambdaMean->GetBinContent(hLambdaMean->FindBin(v0.pt()))) / (hLambdaWidth->GetBinContent(hLambdaMean->FindBin(v0.pt())) + 1e-6);
          massNSigmaAntiLambda = (v0.mAntiLambda() - hLambdaMean->GetBinContent(hLambdaMean->FindBin(v0.pt()))) / (hLambdaWidth->GetBinContent(hLambdaMean->FindBin(v0.pt())) + 1e-6);
        } else {
          massNSigmaK0Short = (v0.mK0Short() - hK0ShortMean->Interpolate(v0.pt())) / (hK0ShortWidth->Interpolate(v0.pt()) + 1e-6);
          massNSigmaLambda = (v0.mLambda() - hLambdaMean->Interpolate(v0.pt())) / (hLambdaWidth->Interpolate(v0.pt()) + 1e-6);
          massNSigmaAntiLambda = (v0.mAntiLambda() - hLambdaMean->Interpolate(v0.pt())) / (hLambdaWidth->Interpolate(v0.pt()) + 1e-6);
        }
      }
      if (compatibleK0Short)
        histos.fill(HIST("h3dMassK0Short"), v0.pt(), v0.mK0Short(), cent);
      if (compatibleLambda)
        histos.fill(HIST("h3dMassLambda"), v0.pt(), v0.mLambda(), cent);
      if (compatibleAntiLambda)
        histos.fill(HIST("h3dMassAntiLambda"), v0.pt(), v0.mAntiLambda(), cent);

      if (!fillTableOnlyWithCompatible ||
          ( // start major condition check
            (compatibleK0Short > 0 && std::abs(massNSigmaK0Short) < maxMassNSigma) ||
            (compatibleLambda > 0 && std::abs(massNSigmaLambda) < maxMassNSigma) ||
            (compatibleAntiLambda > 0 && std::abs(massNSigmaAntiLambda) < maxMassNSigma)) // end major condition check
      ) {
        assocV0(v0.collisionId(), v0.globalIndex(),
                compatibleK0Short, compatibleLambda, compatibleAntiLambda,
                false, false, false, false,
                massNSigmaK0Short, massNSigmaLambda, massNSigmaAntiLambda);
      }
    }
  }

  void processV0sMC(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::PVMults>::iterator const& collision, DauTracksMC const&, soa::Filtered<V0DatasWithoutTrackXMC> const& V0s, aod::McParticles const&, aod::BCsWithTimestamps const&)
  {
    double cent = doPPAnalysis ? collision.centFT0M() : collision.centFT0C();
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    if (((doPPAnalysis && !isCollisionSelected(collision))) || (!doPPAnalysis && !isCollisionSelectedPbPb(collision))) {
      return;
    }
    /// _________________________________________________
    /// Populate table with associated V0s

    for (auto const& v0 : V0s) {
      if (v0.v0radius() < systCuts.v0RadiusMin || v0.v0radius() > systCuts.v0RadiusMax || v0.eta() > systCuts.assocEtaMax || v0.eta() < systCuts.assocEtaMin || v0.v0cosPA() < systCuts.v0Cospa) {
        continue;
      }
      if (v0.pt() > systCuts.assocPtCutMax || v0.pt() < systCuts.assocPtCutMin) {
        continue;
      }
      // check dE/dx compatibility
      int compatibleK0Short = 0;
      int compatibleLambda = 0;
      int compatibleAntiLambda = 0;

      auto posdau = v0.posTrack_as<DauTracksMC>();
      auto negdau = v0.negTrack_as<DauTracksMC>();

      if (negdau.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRows)
        continue;
      if (posdau.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRows)
        continue;

      if (std::abs(posdau.tpcNSigmaPi()) < strangedEdxNSigmaLoose && std::abs(negdau.tpcNSigmaPi()) < strangedEdxNSigmaLoose) {
        if (doPPAnalysis || (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short < systCuts.lifetimecutK0S &&
                             std::abs(v0.dcapostopv()) > systCuts.dcapostopvK0S && std::abs(v0.dcanegtopv()) > systCuts.dcanegtopvK0S &&
                             v0.qtarm() * systCuts.armPodCut > std::abs(v0.alpha()))) {
          BIT_SET(compatibleK0Short, 0);
        }
      }
      if (std::abs(posdau.tpcNSigmaPi()) < strangedEdxNSigma && std::abs(negdau.tpcNSigmaPi()) < strangedEdxNSigma) {
        if (doPPAnalysis || (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short < systCuts.lifetimecutK0S &&
                             std::abs(v0.dcapostopv()) > systCuts.dcapostopvK0S && std::abs(v0.dcanegtopv()) > systCuts.dcanegtopvK0S &&
                             v0.qtarm() * systCuts.armPodCut > std::abs(v0.alpha()))) {
          BIT_SET(compatibleK0Short, 1);
        }
      }
      if (std::abs(posdau.tpcNSigmaPi()) < strangedEdxNSigmaTight && std::abs(negdau.tpcNSigmaPi()) < strangedEdxNSigmaTight) {
        if (doPPAnalysis || (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short < systCuts.lifetimecutK0S &&
                             std::abs(v0.dcapostopv()) > systCuts.dcapostopvK0S && std::abs(v0.dcanegtopv()) > systCuts.dcanegtopvK0S &&
                             v0.qtarm() * systCuts.armPodCut > std::abs(v0.alpha()))) {
          BIT_SET(compatibleK0Short, 2);
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
      if (std::abs(posdau.tpcNSigmaPr()) < strangedEdxNSigma && std::abs(negdau.tpcNSigmaPi()) < strangedEdxNSigma) {
        if (v0.v0cosPA() > systCuts.lambdaCospa) {
          if (doPPAnalysis || (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short < systCuts.lifetimecutLambda &&
                               std::abs(v0.dcapostopv()) > systCuts.dcapostopvLambda && std::abs(v0.dcanegtopv()) > systCuts.dcanegtopvLambda)) {
            BIT_SET(compatibleLambda, 1);
          }
        }
      }
      if (std::abs(posdau.tpcNSigmaPr()) < strangedEdxNSigmaTight && std::abs(negdau.tpcNSigmaPi()) < strangedEdxNSigmaTight) {
        if (v0.v0cosPA() > systCuts.lambdaCospa) {
          if (doPPAnalysis || (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short < systCuts.lifetimecutLambda &&
                               std::abs(v0.dcapostopv()) > systCuts.dcapostopvLambda && std::abs(v0.dcanegtopv()) > systCuts.dcanegtopvLambda)) {
            BIT_SET(compatibleLambda, 2);
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
      if (std::abs(posdau.tpcNSigmaPi()) < strangedEdxNSigma && std::abs(negdau.tpcNSigmaPr()) < strangedEdxNSigma) {
        if (v0.v0cosPA() > systCuts.lambdaCospa) {
          if (doPPAnalysis || (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short < systCuts.lifetimecutLambda &&
                               std::abs(v0.dcapostopv()) > systCuts.dcapostopvAntiLambda && std::abs(v0.dcanegtopv()) > systCuts.dcanegtopvAntiLambda)) {
            BIT_SET(compatibleAntiLambda, 1);
          }
        }
      }
      if (std::abs(posdau.tpcNSigmaPi()) < strangedEdxNSigmaTight && std::abs(negdau.tpcNSigmaPr()) < strangedEdxNSigmaTight) {
        if (v0.v0cosPA() > systCuts.lambdaCospa) {
          if (doPPAnalysis || (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short < systCuts.lifetimecutLambda &&
                               std::abs(v0.dcapostopv()) > systCuts.dcapostopvAntiLambda && std::abs(v0.dcanegtopv()) > systCuts.dcanegtopvAntiLambda)) {
            BIT_SET(compatibleAntiLambda, 2);
          }
        }
      }

      float massNSigmaK0Short = -20.0f;
      float massNSigmaLambda = -20.0f;
      float massNSigmaAntiLambda = -20.0f;
      if (useParameterization) {
        massNSigmaK0Short = (v0.mK0Short() - fK0Mean->Eval(v0.pt())) / (fK0Width->Eval(v0.pt()) + 1e-6);
        massNSigmaLambda = (v0.mLambda() - fLambdaMean->Eval(v0.pt())) / (fLambdaWidth->Eval(v0.pt()) + 1e-6);
        massNSigmaAntiLambda = (v0.mAntiLambda() - fLambdaMean->Eval(v0.pt())) / (fLambdaWidth->Eval(v0.pt()) + 1e-6);
      } else {
        // Load parameters for sideband subtraction
        initParametersFromCCDB(bc);
        // simplified handling: calculate NSigma in mass here
        if (v0.pt() < 0.2f || v0.pt() > 14.5f) {
          massNSigmaK0Short = (v0.mK0Short() - hK0ShortMean->GetBinContent(hK0ShortMean->FindBin(v0.pt()))) / (hK0ShortWidth->GetBinContent(hK0ShortWidth->FindBin(v0.pt())) + 1e-6);
          massNSigmaLambda = (v0.mLambda() - hLambdaMean->GetBinContent(hLambdaMean->FindBin(v0.pt()))) / (hLambdaWidth->GetBinContent(hLambdaMean->FindBin(v0.pt())) + 1e-6);
          massNSigmaAntiLambda = (v0.mAntiLambda() - hLambdaMean->GetBinContent(hLambdaMean->FindBin(v0.pt()))) / (hLambdaWidth->GetBinContent(hLambdaMean->FindBin(v0.pt())) + 1e-6);
        } else {
          massNSigmaK0Short = (v0.mK0Short() - hK0ShortMean->Interpolate(v0.pt())) / (hK0ShortWidth->Interpolate(v0.pt()) + 1e-6);
          massNSigmaLambda = (v0.mLambda() - hLambdaMean->Interpolate(v0.pt())) / (hLambdaWidth->Interpolate(v0.pt()) + 1e-6);
          massNSigmaAntiLambda = (v0.mAntiLambda() - hLambdaMean->Interpolate(v0.pt())) / (hLambdaWidth->Interpolate(v0.pt()) + 1e-6);
        }
      }
      bool v0PhysicalPrimary = false;
      bool trueK0Short = false;
      bool trueLambda = false;
      bool trueAntiLambda = false;
      v0PhysicalPrimary = v0.isPhysicalPrimary();
      if (v0.pdgCode() == 310)
        trueK0Short = true;
      if (v0.pdgCode() == 3122)
        trueLambda = true;
      if (v0.pdgCode() == -3122)
        trueAntiLambda = true;
      if (compatibleK0Short && (!doTrueSelectionInMass || (trueK0Short && v0PhysicalPrimary)))
        histos.fill(HIST("h3dMassK0Short"), v0.pt(), v0.mK0Short(), cent);
      if (compatibleLambda && (!doTrueSelectionInMass || (trueLambda && v0PhysicalPrimary)))
        histos.fill(HIST("h3dMassLambda"), v0.pt(), v0.mLambda(), cent);
      if (compatibleAntiLambda && (!doTrueSelectionInMass || (trueAntiLambda && v0PhysicalPrimary)))
        histos.fill(HIST("h3dMassAntiLambda"), v0.pt(), v0.mAntiLambda(), cent);

      if (!fillTableOnlyWithCompatible ||
          ( // start major condition check
            (compatibleK0Short > 0 && std::abs(massNSigmaK0Short) < maxMassNSigma) ||
            (compatibleLambda > 0 && std::abs(massNSigmaLambda) < maxMassNSigma) ||
            (compatibleAntiLambda > 0 && std::abs(massNSigmaAntiLambda) < maxMassNSigma)) // end major condition check
      ) {
        assocV0(v0.collisionId(), v0.globalIndex(),
                compatibleK0Short, compatibleLambda, compatibleAntiLambda,
                trueK0Short, trueLambda, trueAntiLambda, v0PhysicalPrimary,
                massNSigmaK0Short, massNSigmaLambda, massNSigmaAntiLambda);
      }
    }
  }

  void processCascades(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::PVMults>::iterator const& collision, DauTracks const&, soa::Filtered<V0DatasWithoutTrackX> const& /*V0s*/, soa::Filtered<aod::CascDatas> const& Cascades, aod::BCsWithTimestamps const&)
  {
    double cent = doPPAnalysis ? collision.centFT0M() : collision.centFT0C();
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    if (((doPPAnalysis && !isCollisionSelected(collision))) || (!doPPAnalysis && !isCollisionSelectedPbPb(collision))) {
      return;
    }
    /// _________________________________________________
    /// Step 3: Populate table with associated Cascades
    for (auto const& casc : Cascades) {
      if (casc.eta() > systCuts.assocEtaMax || casc.eta() < systCuts.assocEtaMin) {
        continue;
      }
      if (casc.pt() > systCuts.assocPtCutMax || casc.pt() < systCuts.assocPtCutMin) {
        continue;
      }
      if (doPPAnalysis && (casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < systCuts.v0Cospa ||
                           casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) < systCuts.cascCospa ||
                           casc.cascradius() < systCuts.cascRadius ||
                           std::abs(casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ())) < systCuts.cascMindcav0topv ||
                           std::abs(casc.mLambda() - o2::constants::physics::MassLambda0) > systCuts.cascV0masswindow))
        continue;
      auto bachTrackCast = casc.bachelor_as<DauTracks>();
      auto posTrackCast = casc.posTrack_as<DauTracks>();
      auto negTrackCast = casc.negTrack_as<DauTracks>();

      // minimum TPC crossed rows
      if (bachTrackCast.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRows)
        continue;
      if (posTrackCast.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRows)
        continue;
      if (negTrackCast.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRows)
        continue;
      if (!doPPAnalysis && !CascadeSelectedPbPb(casc, collision.posX(), collision.posY(), collision.posZ()))
        continue;
      // check dE/dx compatibility
      int compatibleXiMinus = 0;
      int compatibleXiPlus = 0;
      int compatibleOmegaMinus = 0;
      int compatibleOmegaPlus = 0;
      float cascpos = std::hypot(casc.x() - collision.posX(), casc.y() - collision.posY(), casc.z() - collision.posZ());
      float cascptotmom = std::hypot(casc.px(), casc.py(), casc.pz());
      float ctauXi = o2::constants::physics::MassXiMinus * cascpos / ((cascptotmom + 1e-13) * ctauxiPDG);
      float ctauOmega = o2::constants::physics::MassOmegaMinus * cascpos / ((cascptotmom + 1e-13) * ctauomegaPDG);

      if (std::abs(posTrackCast.tpcNSigmaPr()) < strangedEdxNSigmaLoose && std::abs(negTrackCast.tpcNSigmaPi()) < strangedEdxNSigmaLoose && std::abs(bachTrackCast.tpcNSigmaPi()) < strangedEdxNSigmaLoose && casc.sign() < 0) {
        if (doPPAnalysis || (std::abs(casc.dcabachtopv()) > MorePbPbsystCuts.dcaBachToPV && std::abs(casc.dcapostopv()) > MorePbPbsystCuts.dcaBaryonToPV &&
                             std::abs(casc.dcanegtopv()) > MorePbPbsystCuts.dcaMesonToPV && std::abs(casc.mOmega() - o2::constants::physics::MassOmegaMinus) > MorePbPbsystCuts.rejcomp &&
                             ctauXi < MorePbPbsystCuts.proplifetime && std::abs(casc.yXi()) < MorePbPbsystCuts.rapCut)) {
          BIT_SET(compatibleXiMinus, 0);
        }
      }
      if (std::abs(posTrackCast.tpcNSigmaPr()) < strangedEdxNSigma && std::abs(negTrackCast.tpcNSigmaPi()) < strangedEdxNSigma && std::abs(bachTrackCast.tpcNSigmaPi()) < strangedEdxNSigma && casc.sign() < 0) {
        if (doPPAnalysis || (std::abs(casc.dcabachtopv()) > MorePbPbsystCuts.dcaBachToPV && std::abs(casc.dcapostopv()) > MorePbPbsystCuts.dcaBaryonToPV &&
                             std::abs(casc.dcanegtopv()) > MorePbPbsystCuts.dcaMesonToPV && std::abs(casc.mOmega() - o2::constants::physics::MassOmegaMinus) > MorePbPbsystCuts.rejcomp &&
                             ctauXi < MorePbPbsystCuts.proplifetime && std::abs(casc.yXi()) < MorePbPbsystCuts.rapCut)) {
          BIT_SET(compatibleXiMinus, 1);
        }
      }
      if (std::abs(posTrackCast.tpcNSigmaPr()) < strangedEdxNSigmaTight && std::abs(negTrackCast.tpcNSigmaPi()) < strangedEdxNSigmaTight && std::abs(bachTrackCast.tpcNSigmaPi()) < strangedEdxNSigmaTight && casc.sign() < 0) {
        if (doPPAnalysis || (std::abs(casc.dcabachtopv()) > MorePbPbsystCuts.dcaBachToPV && std::abs(casc.dcapostopv()) > MorePbPbsystCuts.dcaBaryonToPV &&
                             std::abs(casc.dcanegtopv()) > MorePbPbsystCuts.dcaMesonToPV && std::abs(casc.mOmega() - o2::constants::physics::MassOmegaMinus) > MorePbPbsystCuts.rejcomp &&
                             ctauXi < MorePbPbsystCuts.proplifetime && std::abs(casc.yXi()) < MorePbPbsystCuts.rapCut)) {
          BIT_SET(compatibleXiMinus, 2);
        }
      }

      if (std::abs(posTrackCast.tpcNSigmaPi()) < strangedEdxNSigmaLoose && std::abs(negTrackCast.tpcNSigmaPr()) < strangedEdxNSigmaLoose && std::abs(bachTrackCast.tpcNSigmaPi()) < strangedEdxNSigmaLoose && casc.sign() > 0) {
        if (doPPAnalysis || (std::abs(casc.dcabachtopv()) > MorePbPbsystCuts.dcaBachToPV && std::abs(casc.dcapostopv()) > MorePbPbsystCuts.dcaBaryonToPV &&
                             std::abs(casc.dcanegtopv()) > MorePbPbsystCuts.dcaMesonToPV && std::abs(casc.mOmega() - o2::constants::physics::MassOmegaMinus) > MorePbPbsystCuts.rejcomp &&
                             ctauXi < MorePbPbsystCuts.proplifetime && std::abs(casc.yXi()) < MorePbPbsystCuts.rapCut)) {
          BIT_SET(compatibleXiPlus, 0);
        }
      }
      if (std::abs(posTrackCast.tpcNSigmaPi()) < strangedEdxNSigma && std::abs(negTrackCast.tpcNSigmaPr()) < strangedEdxNSigma && std::abs(bachTrackCast.tpcNSigmaPi()) < strangedEdxNSigma && casc.sign() > 0) {
        if (doPPAnalysis || (std::abs(casc.dcabachtopv()) > MorePbPbsystCuts.dcaBachToPV && std::abs(casc.dcapostopv()) > MorePbPbsystCuts.dcaBaryonToPV &&
                             std::abs(casc.dcanegtopv()) > MorePbPbsystCuts.dcaMesonToPV && std::abs(casc.mOmega() - o2::constants::physics::MassOmegaMinus) > MorePbPbsystCuts.rejcomp &&
                             ctauXi < MorePbPbsystCuts.proplifetime && std::abs(casc.yXi()) < MorePbPbsystCuts.rapCut)) {
          BIT_SET(compatibleXiPlus, 1);
        }
      }
      if (std::abs(posTrackCast.tpcNSigmaPi()) < strangedEdxNSigmaTight && std::abs(negTrackCast.tpcNSigmaPr()) < strangedEdxNSigmaTight && std::abs(bachTrackCast.tpcNSigmaPi()) < strangedEdxNSigmaTight && casc.sign() > 0) {
        if (doPPAnalysis || (std::abs(casc.dcabachtopv()) > MorePbPbsystCuts.dcaBachToPV && std::abs(casc.dcapostopv()) > MorePbPbsystCuts.dcaBaryonToPV &&
                             std::abs(casc.dcanegtopv()) > MorePbPbsystCuts.dcaMesonToPV && std::abs(casc.mOmega() - o2::constants::physics::MassOmegaMinus) > MorePbPbsystCuts.rejcomp &&
                             ctauXi < MorePbPbsystCuts.proplifetime && std::abs(casc.yXi()) < MorePbPbsystCuts.rapCut)) {
          BIT_SET(compatibleXiPlus, 2);
        }
      }

      if (std::abs(posTrackCast.tpcNSigmaPr()) < strangedEdxNSigmaLoose && std::abs(negTrackCast.tpcNSigmaPi()) < strangedEdxNSigmaLoose && std::abs(bachTrackCast.tpcNSigmaKa()) < strangedEdxNSigmaLoose && casc.sign() < 0) {
        if (doPPAnalysis || (std::abs(casc.dcabachtopv()) > MorePbPbsystCuts.dcaBachToPV && std::abs(casc.dcapostopv()) > MorePbPbsystCuts.dcaBaryonToPV &&
                             std::abs(casc.dcanegtopv()) > MorePbPbsystCuts.dcaMesonToPV && std::abs(casc.mXi() - o2::constants::physics::MassXiMinus) > MorePbPbsystCuts.rejcomp &&
                             ctauOmega < MorePbPbsystCuts.proplifetime && std::abs(casc.yOmega()) < MorePbPbsystCuts.rapCut)) {
          BIT_SET(compatibleOmegaMinus, 0);
        }
      }
      if (std::abs(posTrackCast.tpcNSigmaPr()) < strangedEdxNSigma && std::abs(negTrackCast.tpcNSigmaPi()) < strangedEdxNSigma && std::abs(bachTrackCast.tpcNSigmaKa()) < strangedEdxNSigma && casc.sign() < 0) {
        if (doPPAnalysis || (std::abs(casc.dcabachtopv()) > MorePbPbsystCuts.dcaBachToPV && std::abs(casc.dcapostopv()) > MorePbPbsystCuts.dcaBaryonToPV &&
                             std::abs(casc.dcanegtopv()) > MorePbPbsystCuts.dcaMesonToPV && std::abs(casc.mXi() - o2::constants::physics::MassXiMinus) > MorePbPbsystCuts.rejcomp &&
                             ctauOmega < MorePbPbsystCuts.proplifetime && std::abs(casc.yOmega()) < MorePbPbsystCuts.rapCut)) {
          BIT_SET(compatibleOmegaMinus, 1);
        }
      }
      if (std::abs(posTrackCast.tpcNSigmaPr()) < strangedEdxNSigmaTight && std::abs(negTrackCast.tpcNSigmaPi()) < strangedEdxNSigmaTight && std::abs(bachTrackCast.tpcNSigmaKa()) < strangedEdxNSigmaTight && casc.sign() < 0) {
        if (doPPAnalysis || (std::abs(casc.dcabachtopv()) > MorePbPbsystCuts.dcaBachToPV && std::abs(casc.dcapostopv()) > MorePbPbsystCuts.dcaBaryonToPV &&
                             std::abs(casc.dcanegtopv()) > MorePbPbsystCuts.dcaMesonToPV && std::abs(casc.mXi() - o2::constants::physics::MassXiMinus) > MorePbPbsystCuts.rejcomp &&
                             ctauOmega < MorePbPbsystCuts.proplifetime && std::abs(casc.yOmega()) < MorePbPbsystCuts.rapCut)) {
          BIT_SET(compatibleOmegaMinus, 2);
        }
      }

      if (std::abs(posTrackCast.tpcNSigmaPi()) < strangedEdxNSigmaLoose && std::abs(negTrackCast.tpcNSigmaPr()) < strangedEdxNSigmaLoose && std::abs(bachTrackCast.tpcNSigmaKa()) < strangedEdxNSigmaLoose && casc.sign() > 0) {
        if (doPPAnalysis || (std::abs(casc.dcabachtopv()) > MorePbPbsystCuts.dcaBachToPV && std::abs(casc.dcapostopv()) > MorePbPbsystCuts.dcaMesonToPV &&
                             std::abs(casc.dcanegtopv()) > MorePbPbsystCuts.dcaBaryonToPV && std::abs(casc.mXi() - o2::constants::physics::MassXiMinus) > MorePbPbsystCuts.rejcomp &&
                             ctauOmega < MorePbPbsystCuts.proplifetime && std::abs(casc.yOmega()) < MorePbPbsystCuts.rapCut)) {
          BIT_SET(compatibleOmegaPlus, 0);
        }
      }
      if (std::abs(posTrackCast.tpcNSigmaPi()) < strangedEdxNSigma && std::abs(negTrackCast.tpcNSigmaPr()) < strangedEdxNSigma && std::abs(bachTrackCast.tpcNSigmaKa()) < strangedEdxNSigma && casc.sign() > 0) {
        if (doPPAnalysis || (std::abs(casc.dcabachtopv()) > MorePbPbsystCuts.dcaBachToPV && std::abs(casc.dcapostopv()) > MorePbPbsystCuts.dcaMesonToPV &&
                             std::abs(casc.dcanegtopv()) > MorePbPbsystCuts.dcaBaryonToPV && std::abs(casc.mXi() - o2::constants::physics::MassXiMinus) > MorePbPbsystCuts.rejcomp &&
                             ctauOmega < MorePbPbsystCuts.proplifetime && std::abs(casc.yOmega()) < MorePbPbsystCuts.rapCut)) {
          BIT_SET(compatibleOmegaPlus, 1);
        }
      }
      if (std::abs(posTrackCast.tpcNSigmaPi()) < strangedEdxNSigmaTight && std::abs(negTrackCast.tpcNSigmaPr()) < strangedEdxNSigmaTight && std::abs(bachTrackCast.tpcNSigmaKa()) < strangedEdxNSigmaTight && casc.sign() > 0) {
        if (doPPAnalysis || (std::abs(casc.dcabachtopv()) > MorePbPbsystCuts.dcaBachToPV && std::abs(casc.dcapostopv()) > MorePbPbsystCuts.dcaMesonToPV &&
                             std::abs(casc.dcanegtopv()) > MorePbPbsystCuts.dcaBaryonToPV && std::abs(casc.mXi() - o2::constants::physics::MassXiMinus) > MorePbPbsystCuts.rejcomp &&
                             ctauOmega < MorePbPbsystCuts.proplifetime && std::abs(casc.yOmega()) < MorePbPbsystCuts.rapCut)) {
          BIT_SET(compatibleOmegaPlus, 2);
        }
      }
      float massNSigmaXi = -20.0f;
      float massNSigmaOmega = -20.0f;
      if (useParameterization) {
        massNSigmaXi = (casc.mXi() - fXiMean->Eval(casc.pt())) / (fXiWidth->Eval(casc.pt()) + 1e-6);
        massNSigmaOmega = (casc.mOmega() - fOmegaMean->Eval(casc.pt())) / (fOmegaWidth->Eval(casc.pt()) + 1e-6);
      } else {
        // Load parameters for sideband subtraction
        initParametersFromCCDB(bc);
        if (casc.pt() < 0.2f || casc.pt() > 14.5f) {
          massNSigmaXi = (casc.mXi() - hXiMean->GetBinContent(hXiMean->FindBin(casc.pt()))) / (hXiWidth->GetBinContent(hXiWidth->FindBin(casc.pt())) + 1e-6);
          massNSigmaOmega = (casc.mOmega() - hOmegaMean->GetBinContent(hOmegaMean->FindBin(casc.pt()))) / (hOmegaWidth->GetBinContent(hOmegaWidth->FindBin(casc.pt())) + 1e-6);
        } else {
          massNSigmaXi = (casc.mXi() - hXiMean->Interpolate(casc.pt())) / (hXiWidth->Interpolate(casc.pt()) + 1e-6);
          massNSigmaOmega = (casc.mOmega() - hOmegaMean->Interpolate(casc.pt())) / (hOmegaWidth->Interpolate(casc.pt()) + 1e-6);
        }
      }
      if (compatibleXiMinus)
        histos.fill(HIST("h3dMassXiMinus"), casc.pt(), casc.mXi(), cent);
      if (compatibleXiPlus)
        histos.fill(HIST("h3dMassXiPlus"), casc.pt(), casc.mXi(), cent);
      if (compatibleOmegaMinus && std::abs(massNSigmaXi) > nSigmaNearXiMassCenter)
        histos.fill(HIST("h3dMassOmegaMinus"), casc.pt(), casc.mOmega(), cent);
      if (compatibleOmegaPlus && std::abs(massNSigmaXi) > nSigmaNearXiMassCenter)
        histos.fill(HIST("h3dMassOmegaPlus"), casc.pt(), casc.mOmega(), cent);

      if (!fillTableOnlyWithCompatible ||
          ( // start major condition check
            ((compatibleXiMinus > 0 || compatibleXiPlus > 0) && std::abs(massNSigmaXi) < maxMassNSigma) ||
            ((compatibleOmegaMinus > 0 || compatibleOmegaPlus > 0) && std::abs(massNSigmaOmega) < maxMassNSigma && std::abs(massNSigmaXi) > nSigmaNearXiMassCenter)) // end major condition check
      ) {
        assocCascades(casc.collisionId(), casc.globalIndex(),
                      compatibleXiMinus, compatibleXiPlus, compatibleOmegaMinus, compatibleOmegaPlus,
                      false, false, false, false, false,
                      massNSigmaXi, massNSigmaOmega);
      }
    }
  }

  void processCascadesMC(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::PVMults>::iterator const& collision, DauTracks const&, soa::Filtered<V0DatasWithoutTrackXMC> const& /*V0s*/, soa::Filtered<CascDatasMC> const& Cascades, aod::McParticles const&, aod::BCsWithTimestamps const&)
  {
    double cent = doPPAnalysis ? collision.centFT0M() : collision.centFT0C();
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    // Perform basic event selection
    if (((doPPAnalysis && !isCollisionSelected(collision))) || (!doPPAnalysis && !isCollisionSelectedPbPb(collision))) {
      return;
    }
    /// _________________________________________________
    /// Step 3: Populate table with associated Cascades
    for (auto const& casc : Cascades) {
      if (casc.eta() > systCuts.assocEtaMax || casc.eta() < systCuts.assocEtaMin) {
        continue;
      }
      if (casc.pt() > systCuts.assocPtCutMax || casc.pt() < systCuts.assocPtCutMin) {
        continue;
      }
      if (doPPAnalysis && (casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < systCuts.v0Cospa ||
                           casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) < systCuts.cascCospa ||
                           casc.cascradius() < systCuts.cascRadius ||
                           std::abs(casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ())) < systCuts.cascMindcav0topv ||
                           std::abs(casc.mLambda() - o2::constants::physics::MassLambda0) > systCuts.cascV0masswindow))
        continue;

      auto bachTrackCast = casc.bachelor_as<DauTracks>();
      auto posTrackCast = casc.posTrack_as<DauTracks>();
      auto negTrackCast = casc.negTrack_as<DauTracks>();

      // minimum TPC crossed rows
      if (bachTrackCast.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRows)
        continue;
      if (posTrackCast.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRows)
        continue;
      if (negTrackCast.tpcNClsCrossedRows() < systCuts.minTPCNCrossedRows)
        continue;
      if (!doPPAnalysis && !CascadeSelectedPbPb(casc, collision.posX(), collision.posY(), collision.posZ()))
        continue;

      // check dE/dx compatibility
      int compatibleXiMinus = 0;
      int compatibleXiPlus = 0;
      int compatibleOmegaMinus = 0;
      int compatibleOmegaPlus = 0;
      float cascpos = std::hypot(casc.x() - collision.posX(), casc.y() - collision.posY(), casc.z() - collision.posZ());
      float cascptotmom = std::hypot(casc.px(), casc.py(), casc.pz());
      float ctauXi = o2::constants::physics::MassXiMinus * cascpos / ((cascptotmom + 1e-13) * ctauxiPDG);
      float ctauOmega = o2::constants::physics::MassOmegaMinus * cascpos / ((cascptotmom + 1e-13) * ctauomegaPDG);

      if (std::abs(posTrackCast.tpcNSigmaPr()) < strangedEdxNSigmaLoose && std::abs(negTrackCast.tpcNSigmaPi()) < strangedEdxNSigmaLoose && std::abs(bachTrackCast.tpcNSigmaPi()) < strangedEdxNSigmaLoose && casc.sign() < 0) {
        if (doPPAnalysis || (std::abs(casc.dcabachtopv()) > MorePbPbsystCuts.dcaBachToPV && std::abs(casc.dcapostopv()) > MorePbPbsystCuts.dcaBaryonToPV &&
                             std::abs(casc.dcanegtopv()) > MorePbPbsystCuts.dcaMesonToPV && std::abs(casc.mOmega() - o2::constants::physics::MassOmegaMinus) > MorePbPbsystCuts.rejcomp &&
                             ctauXi < MorePbPbsystCuts.proplifetime && std::abs(casc.yXi()) < MorePbPbsystCuts.rapCut)) {
          BIT_SET(compatibleXiMinus, 0);
        }
      }
      if (std::abs(posTrackCast.tpcNSigmaPr()) < strangedEdxNSigma && std::abs(negTrackCast.tpcNSigmaPi()) < strangedEdxNSigma && std::abs(bachTrackCast.tpcNSigmaPi()) < strangedEdxNSigma && casc.sign() < 0) {
        if (doPPAnalysis || (std::abs(casc.dcabachtopv()) > MorePbPbsystCuts.dcaBachToPV && std::abs(casc.dcapostopv()) > MorePbPbsystCuts.dcaBaryonToPV &&
                             std::abs(casc.dcanegtopv()) > MorePbPbsystCuts.dcaMesonToPV && std::abs(casc.mOmega() - o2::constants::physics::MassOmegaMinus) > MorePbPbsystCuts.rejcomp &&
                             ctauXi < MorePbPbsystCuts.proplifetime && std::abs(casc.yXi()) < MorePbPbsystCuts.rapCut)) {
          BIT_SET(compatibleXiMinus, 1);
        }
      }
      if (std::abs(posTrackCast.tpcNSigmaPr()) < strangedEdxNSigmaTight && std::abs(negTrackCast.tpcNSigmaPi()) < strangedEdxNSigmaTight && std::abs(bachTrackCast.tpcNSigmaPi()) < strangedEdxNSigmaTight && casc.sign() < 0) {
        if (doPPAnalysis || (std::abs(casc.dcabachtopv()) > MorePbPbsystCuts.dcaBachToPV && std::abs(casc.dcapostopv()) > MorePbPbsystCuts.dcaBaryonToPV &&
                             std::abs(casc.dcanegtopv()) > MorePbPbsystCuts.dcaMesonToPV && std::abs(casc.mOmega() - o2::constants::physics::MassOmegaMinus) > MorePbPbsystCuts.rejcomp &&
                             ctauXi < MorePbPbsystCuts.proplifetime && std::abs(casc.yXi()) < MorePbPbsystCuts.rapCut)) {
          BIT_SET(compatibleXiMinus, 2);
        }
      }

      if (std::abs(posTrackCast.tpcNSigmaPi()) < strangedEdxNSigmaLoose && std::abs(negTrackCast.tpcNSigmaPr()) < strangedEdxNSigmaLoose && std::abs(bachTrackCast.tpcNSigmaPi()) < strangedEdxNSigmaLoose && casc.sign() > 0) {
        if (doPPAnalysis || (std::abs(casc.dcabachtopv()) > MorePbPbsystCuts.dcaBachToPV && std::abs(casc.dcapostopv()) > MorePbPbsystCuts.dcaBaryonToPV &&
                             std::abs(casc.dcanegtopv()) > MorePbPbsystCuts.dcaMesonToPV && std::abs(casc.mOmega() - o2::constants::physics::MassOmegaMinus) > MorePbPbsystCuts.rejcomp &&
                             ctauXi < MorePbPbsystCuts.proplifetime && std::abs(casc.yXi()) < MorePbPbsystCuts.rapCut)) {
          BIT_SET(compatibleXiPlus, 0);
        }
      }
      if (std::abs(posTrackCast.tpcNSigmaPi()) < strangedEdxNSigma && std::abs(negTrackCast.tpcNSigmaPr()) < strangedEdxNSigma && std::abs(bachTrackCast.tpcNSigmaPi()) < strangedEdxNSigma && casc.sign() > 0) {
        if (doPPAnalysis || (std::abs(casc.dcabachtopv()) > MorePbPbsystCuts.dcaBachToPV && std::abs(casc.dcapostopv()) > MorePbPbsystCuts.dcaBaryonToPV &&
                             std::abs(casc.dcanegtopv()) > MorePbPbsystCuts.dcaMesonToPV && std::abs(casc.mOmega() - o2::constants::physics::MassOmegaMinus) > MorePbPbsystCuts.rejcomp &&
                             ctauXi < MorePbPbsystCuts.proplifetime && std::abs(casc.yXi()) < MorePbPbsystCuts.rapCut)) {
          BIT_SET(compatibleXiPlus, 1);
        }
      }
      if (std::abs(posTrackCast.tpcNSigmaPi()) < strangedEdxNSigmaTight && std::abs(negTrackCast.tpcNSigmaPr()) < strangedEdxNSigmaTight && std::abs(bachTrackCast.tpcNSigmaPi()) < strangedEdxNSigmaTight && casc.sign() > 0) {
        if (doPPAnalysis || (std::abs(casc.dcabachtopv()) > MorePbPbsystCuts.dcaBachToPV && std::abs(casc.dcapostopv()) > MorePbPbsystCuts.dcaBaryonToPV &&
                             std::abs(casc.dcanegtopv()) > MorePbPbsystCuts.dcaMesonToPV && std::abs(casc.mOmega() - o2::constants::physics::MassOmegaMinus) > MorePbPbsystCuts.rejcomp &&
                             ctauXi < MorePbPbsystCuts.proplifetime && std::abs(casc.yXi()) < MorePbPbsystCuts.rapCut)) {
          BIT_SET(compatibleXiPlus, 2);
        }
      }

      if (std::abs(posTrackCast.tpcNSigmaPr()) < strangedEdxNSigmaLoose && std::abs(negTrackCast.tpcNSigmaPi()) < strangedEdxNSigmaLoose && std::abs(bachTrackCast.tpcNSigmaKa()) < strangedEdxNSigmaLoose && casc.sign() < 0) {
        if (doPPAnalysis || (std::abs(casc.dcabachtopv()) > MorePbPbsystCuts.dcaBachToPV && std::abs(casc.dcapostopv()) > MorePbPbsystCuts.dcaBaryonToPV &&
                             std::abs(casc.dcanegtopv()) > MorePbPbsystCuts.dcaMesonToPV && std::abs(casc.mXi() - o2::constants::physics::MassXiMinus) > MorePbPbsystCuts.rejcomp &&
                             ctauOmega < MorePbPbsystCuts.proplifetime && std::abs(casc.yOmega()) < MorePbPbsystCuts.rapCut)) {
          BIT_SET(compatibleOmegaMinus, 0);
        }
      }
      if (std::abs(posTrackCast.tpcNSigmaPr()) < strangedEdxNSigma && std::abs(negTrackCast.tpcNSigmaPi()) < strangedEdxNSigma && std::abs(bachTrackCast.tpcNSigmaKa()) < strangedEdxNSigma && casc.sign() < 0) {
        if (doPPAnalysis || (std::abs(casc.dcabachtopv()) > MorePbPbsystCuts.dcaBachToPV && std::abs(casc.dcapostopv()) > MorePbPbsystCuts.dcaBaryonToPV &&
                             std::abs(casc.dcanegtopv()) > MorePbPbsystCuts.dcaMesonToPV && std::abs(casc.mXi() - o2::constants::physics::MassXiMinus) > MorePbPbsystCuts.rejcomp &&
                             ctauOmega < MorePbPbsystCuts.proplifetime && std::abs(casc.yOmega()) < MorePbPbsystCuts.rapCut)) {
          BIT_SET(compatibleOmegaMinus, 1);
        }
      }
      if (std::abs(posTrackCast.tpcNSigmaPr()) < strangedEdxNSigmaTight && std::abs(negTrackCast.tpcNSigmaPi()) < strangedEdxNSigmaTight && std::abs(bachTrackCast.tpcNSigmaKa()) < strangedEdxNSigmaTight && casc.sign() < 0) {
        if (doPPAnalysis || (std::abs(casc.dcabachtopv()) > MorePbPbsystCuts.dcaBachToPV && std::abs(casc.dcapostopv()) > MorePbPbsystCuts.dcaBaryonToPV &&
                             std::abs(casc.dcanegtopv()) > MorePbPbsystCuts.dcaMesonToPV && std::abs(casc.mXi() - o2::constants::physics::MassXiMinus) > MorePbPbsystCuts.rejcomp &&
                             ctauOmega < MorePbPbsystCuts.proplifetime && std::abs(casc.yOmega()) < MorePbPbsystCuts.rapCut)) {
          BIT_SET(compatibleOmegaMinus, 2);
        }
      }

      if (std::abs(posTrackCast.tpcNSigmaPi()) < strangedEdxNSigmaLoose && std::abs(negTrackCast.tpcNSigmaPr()) < strangedEdxNSigmaLoose && std::abs(bachTrackCast.tpcNSigmaKa()) < strangedEdxNSigmaLoose && casc.sign() > 0) {
        if (doPPAnalysis || (std::abs(casc.dcabachtopv()) > MorePbPbsystCuts.dcaBachToPV && std::abs(casc.dcapostopv()) > MorePbPbsystCuts.dcaMesonToPV &&
                             std::abs(casc.dcanegtopv()) > MorePbPbsystCuts.dcaBaryonToPV && std::abs(casc.mXi() - o2::constants::physics::MassXiMinus) > MorePbPbsystCuts.rejcomp &&
                             ctauOmega < MorePbPbsystCuts.proplifetime && std::abs(casc.yOmega()) < MorePbPbsystCuts.rapCut)) {
          BIT_SET(compatibleOmegaPlus, 0);
        }
      }
      if (std::abs(posTrackCast.tpcNSigmaPi()) < strangedEdxNSigma && std::abs(negTrackCast.tpcNSigmaPr()) < strangedEdxNSigma && std::abs(bachTrackCast.tpcNSigmaKa()) < strangedEdxNSigma && casc.sign() > 0) {
        if (doPPAnalysis || (std::abs(casc.dcabachtopv()) > MorePbPbsystCuts.dcaBachToPV && std::abs(casc.dcapostopv()) > MorePbPbsystCuts.dcaMesonToPV &&
                             std::abs(casc.dcanegtopv()) > MorePbPbsystCuts.dcaBaryonToPV && std::abs(casc.mXi() - o2::constants::physics::MassXiMinus) > MorePbPbsystCuts.rejcomp &&
                             ctauOmega < MorePbPbsystCuts.proplifetime && std::abs(casc.yOmega()) < MorePbPbsystCuts.rapCut)) {
          BIT_SET(compatibleOmegaPlus, 1);
        }
      }
      if (std::abs(posTrackCast.tpcNSigmaPi()) < strangedEdxNSigmaTight && std::abs(negTrackCast.tpcNSigmaPr()) < strangedEdxNSigmaTight && std::abs(bachTrackCast.tpcNSigmaKa()) < strangedEdxNSigmaTight && casc.sign() > 0) {
        if (doPPAnalysis || (std::abs(casc.dcabachtopv()) > MorePbPbsystCuts.dcaBachToPV && std::abs(casc.dcapostopv()) > MorePbPbsystCuts.dcaMesonToPV &&
                             std::abs(casc.dcanegtopv()) > MorePbPbsystCuts.dcaBaryonToPV && std::abs(casc.mXi() - o2::constants::physics::MassXiMinus) > MorePbPbsystCuts.rejcomp &&
                             ctauOmega < MorePbPbsystCuts.proplifetime && std::abs(casc.yOmega()) < MorePbPbsystCuts.rapCut)) {
          BIT_SET(compatibleOmegaPlus, 2);
        }
      }
      float massNSigmaXi = 20.0f;
      float massNSigmaOmega = 20.0f;
      if (useParameterization) {
        massNSigmaXi = (casc.mXi() - fXiMean->Eval(casc.pt())) / (fXiWidth->Eval(casc.pt()) + 1e-6);
        massNSigmaOmega = (casc.mOmega() - fOmegaMean->Eval(casc.pt())) / (fOmegaWidth->Eval(casc.pt()) + 1e-6);
      } else {
        // Load parameters for sideband subtraction
        initParametersFromCCDB(bc);
        if (casc.pt() < 0.2f || casc.pt() > 14.5f) {
          massNSigmaXi = (casc.mXi() - hXiMean->GetBinContent(hXiMean->FindBin(casc.pt()))) / (hXiWidth->GetBinContent(hXiWidth->FindBin(casc.pt())) + 1e-6);
          massNSigmaOmega = (casc.mOmega() - hOmegaMean->GetBinContent(hOmegaMean->FindBin(casc.pt()))) / (hOmegaWidth->GetBinContent(hOmegaWidth->FindBin(casc.pt())) + 1e-6);
        } else {
          massNSigmaXi = (casc.mXi() - hXiMean->Interpolate(casc.pt())) / (hXiWidth->Interpolate(casc.pt()) + 1e-6);
          massNSigmaOmega = (casc.mOmega() - hOmegaMean->Interpolate(casc.pt())) / (hOmegaWidth->Interpolate(casc.pt()) + 1e-6);
        }
      }

      bool cascPhysicalPrimary = false;
      bool trueXiMinus = false;
      bool trueXiPlus = false;
      bool trueOmegaMinus = false;
      bool trueOmegaPlus = false;
      cascPhysicalPrimary = casc.isPhysicalPrimary();
      if (casc.pdgCode() == 3312)
        trueXiMinus = true;
      if (casc.pdgCode() == -3312)
        trueXiPlus = true;
      if (casc.pdgCode() == 3334)
        trueOmegaMinus = true;
      if (casc.pdgCode() == -3334)
        trueOmegaPlus = true;
      if (compatibleXiMinus && (!doTrueSelectionInMass || (trueXiMinus && cascPhysicalPrimary)))
        histos.fill(HIST("h3dMassXiMinus"), casc.pt(), casc.mXi(), cent);
      if (compatibleXiPlus && (!doTrueSelectionInMass || (trueXiPlus && cascPhysicalPrimary)))
        histos.fill(HIST("h3dMassXiPlus"), casc.pt(), casc.mXi(), cent);
      if (compatibleOmegaMinus && (!doTrueSelectionInMass || (trueOmegaMinus && cascPhysicalPrimary)) && std::abs(massNSigmaXi) > nSigmaNearXiMassCenter)
        histos.fill(HIST("h3dMassOmegaMinus"), casc.pt(), casc.mOmega(), cent);
      if (compatibleOmegaPlus && (!doTrueSelectionInMass || (trueOmegaPlus && cascPhysicalPrimary)) && std::abs(massNSigmaXi) > nSigmaNearXiMassCenter)
        histos.fill(HIST("h3dMassOmegaPlus"), casc.pt(), casc.mOmega(), cent);

      if (!fillTableOnlyWithCompatible ||
          ( // start major condition check
            ((compatibleXiMinus > 0 || compatibleXiPlus > 0) && std::abs(massNSigmaXi) < maxMassNSigma) ||
            ((compatibleOmegaMinus > 0 || compatibleOmegaPlus > 0) && std::abs(massNSigmaOmega) < maxMassNSigma && std::abs(massNSigmaXi) > nSigmaNearXiMassCenter)) // end major condition check
      ) {
        assocCascades(casc.collisionId(), casc.globalIndex(),
                      compatibleXiMinus, compatibleXiPlus, compatibleOmegaMinus, compatibleOmegaPlus,
                      trueXiMinus, trueXiPlus,
                      trueOmegaMinus, trueOmegaPlus,
                      cascPhysicalPrimary,
                      massNSigmaXi, massNSigmaOmega);
      }
    }
  }
  PROCESS_SWITCH(HStrangeCorrelationFilter, processTriggers, "Produce trigger tables", true);
  PROCESS_SWITCH(HStrangeCorrelationFilter, processTriggersMC, "Produce trigger tables for MC", false);
  PROCESS_SWITCH(HStrangeCorrelationFilter, processV0s, "Produce associated V0 tables", true);
  PROCESS_SWITCH(HStrangeCorrelationFilter, processV0sMC, "Produce associated V0 tables for MC", false);
  PROCESS_SWITCH(HStrangeCorrelationFilter, processAssocPions, "Produce associated Pion tables", false);
  PROCESS_SWITCH(HStrangeCorrelationFilter, processAssocPionsMC, "Produce associated Pion tables for MC", false);
  PROCESS_SWITCH(HStrangeCorrelationFilter, processCascades, "Produce associated cascade tables", true);
  PROCESS_SWITCH(HStrangeCorrelationFilter, processCascadesMC, "Produce associated cascade tables for MC", false);
  PROCESS_SWITCH(HStrangeCorrelationFilter, processAssocHadrons, "Produce associated Hadron tables", true);
  PROCESS_SWITCH(HStrangeCorrelationFilter, processAssocHadronsMC, "Produce associated Hadron tables for MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HStrangeCorrelationFilter>(cfgc)};
}
