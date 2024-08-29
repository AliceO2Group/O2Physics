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
/// \author Kai Cui (kaicui@mails.ccnu.edu.cn)
/// \author Lucia Anna Tarasovicova (lucia.anna.husova@cern.ch)
/// \author David Dobrigkeit Chinellato (david.dobrigkeit.chinellato@cern.ch)
/// \author Zhongbao Yin (Zhong-Bao.Yin@cern.ch)

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFHStrangeCorrelationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "TF1.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define bitset(var, nbit) ((var) |= (1 << (nbit)))
#define bitcheck(var, nbit) ((var) & (1 << (nbit)))

struct hstrangecorrelationfilter {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Operational
  Configurable<bool> fillTableOnlyWithCompatible{"fillTableOnlyWithCompatible", true, "pre-apply dE/dx, broad mass window in table filling"};
  Configurable<float> strangedEdxNSigmaLoose{"strangedEdxNSigmaLoose", 5, "Nsigmas for strange decay daughters"};
  Configurable<float> strangedEdxNSigma{"strangedEdxNSigma", 4, "Nsigmas for strange decay daughters"};
  Configurable<float> strangedEdxNSigmaTight{"strangedEdxNSigmaTight", 3, "Nsigmas for strange decay daughters"};

  // Trigger particle selections in phase space
  Configurable<float> triggerEtaMin{"triggerEtaCutMin", -0.8, "triggeretamin"};
  Configurable<float> triggerEtaMax{"triggerEtaCutMax", 0.8, "triggeretamax"};
  Configurable<float> triggerPtCutMin{"triggerPtCutMin", 3, "triggerptmin"};
  Configurable<float> triggerPtCutMax{"triggerPtCutMax", 20, "triggerptmax"};

  // Track quality
  Configurable<int> minTPCNCrossedRows{"minTPCNCrossedRows", 70, "Minimum TPC crossed rows"};
  Configurable<bool> triggerRequireITS{"triggerRequireITS", true, "require ITS signal in trigger tracks"};
  Configurable<int> triggerMaxTPCSharedClusters{"triggerMaxTPCSharedClusters", 200, "maximum number of shared TPC clusters (inclusive)"};
  Configurable<bool> triggerRequireL0{"triggerRequireL0", false, "require ITS L0 cluster for trigger"};

  // Associated particle selections in phase space
  Configurable<float> assocEtaMin{"assocEtaCutMin", -0.8, "triggeretamin"};
  Configurable<float> assocEtaMax{"assocEtaCutMax", 0.8, "triggeretamax"};
  Configurable<float> assocPtCutMin{"assocPtCutMin", 0.2, "assocptmin"};
  Configurable<float> assocPtCutMax{"assocPtCutMax", 10, "assocptmax"};

  // Associated pion identification
  Configurable<float> pionMinBayesProb{"pionMinBayesProb", 0.95, "minimal Bayesian probability for pion ID"};
  Configurable<float> assocPionNSigmaTPCFOF{"assocPionNSigmaTPCFOF", 3, "minimal n sigma in TOF and TPC for Pion ID"};
  Configurable<float> rejectSigma{"rejectSigma", 1, "n sigma for rejecting pion candidates"};

  // V0 selections
  Configurable<double> v0Cospa{"v0cospa", 0.97, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<float> dcaV0dau{"dcav0dau", 1.0, "DCA V0 Daughters"};
  Configurable<float> dcaNegtopv{"dcanegtopv", 0.06, "DCA Neg To PV"};
  Configurable<float> dcaPostopv{"dcapostopv", 0.06, "DCA Pos To PV"};
  Configurable<float> v0RadiusMin{"v0radiusmin", 0.5, "v0radius"};
  Configurable<float> v0RadiusMax{"v0radiusmax", 200, "v0radius"};

  // specific selections
  Configurable<double> lambdaCospa{"lambdaCospa", 0.995, "CosPA for lambda"}; // allows for tighter selection for Lambda

  // primary particle DCAxy selections
  // formula: |DCAxy| <  0.004f + (0.013f / pt)
  Configurable<float> dcaXYconstant{"dcaXYconstant", 0.004, "[0] in |DCAxy| < [0]+[1]/pT"};
  Configurable<float> dcaXYpTdep{"dcaXYpTdep", 0.013, "[1] in |DCAxy| < [0]+[1]/pT"};

  // cascade selections
  Configurable<double> cascadesetting_cospa{"cascadesetting_cospa", 0.95, "cascadesetting_cospa"};
  Configurable<float> cascadesetting_dcacascdau{"cascadesetting_dcacascdau", 1.0, "cascadesetting_dcacascdau"};
  Configurable<float> cascadesetting_dcabachtopv{"cascadesetting_dcabachtopv", 0.1, "cascadesetting_dcabachtopv"};
  Configurable<float> cascadesetting_cascradius{"cascadesetting_cascradius", 0.5, "cascadesetting_cascradius"};
  Configurable<float> cascadesetting_v0masswindow{"cascadesetting_v0masswindow", 0.01, "cascadesetting_v0masswindow"};
  Configurable<float> cascadesetting_mindcav0topv{"cascadesetting_mindcav0topv", 0.01, "cascadesetting_mindcav0topv"};

  // invariant mass parametrizations
  Configurable<std::vector<float>> massParsK0Mean{"massParsK0Mean", {0.495, 0.000250, 0.0, 0.0}, "pars in [0]+[1]*x+[2]*TMath::Exp(-[3]*x)"};
  Configurable<std::vector<float>> massParsK0Width{"massParsK0Width", {0.00354, 0.000609, 0.0, 0.0}, "pars in [0]+[1]*x+[2]*TMath::Exp(-[3]*x)"};

  Configurable<std::vector<float>> massParsLambdaMean{"massParsLambdaMean", {1.114, 0.000314, 0.140, 11.9}, "pars in [0]+[1]*x+[2]*TMath::Exp(-[3]*x)"};
  Configurable<std::vector<float>> massParsLambdaWidth{"massParsLambdaWidth", {0.00127, 0.000172, 0.00261, 2.02}, "pars in [0]+[1]*x+[2]*TMath::Exp(-[3]*x)"};

  Configurable<std::vector<float>> massParsCascadeMean{"massParsCascadeMean", {1.32, 0.000278, 0.0, 0.0}, "pars in [0]+[1]*x+[2]*TMath::Exp(-[3]*x)"};
  Configurable<std::vector<float>> massParsCascadeWidth{"massParsCascadeWidth", {0.00189, 0.000227, 0.00370, 1.635}, "pars in [0]+[1]*x+[2]*TMath::Exp(-[3]*x)"};

  Configurable<std::vector<float>> massParsOmegaMean{"massParsOmegaMean", {1.67, 0.000298, 0.0, 0.0}, "pars in [0]+[1]*x+[2]*TMath::Exp(-[3]*x)"};
  Configurable<std::vector<float>> massParsOmegaWidth{"massParsOmegaWidth", {0.00189, 0.000325, 0.00606, 1.77}, "pars in [0]+[1]*x+[2]*TMath::Exp(-[3]*x)"};

  // must include windows for background and peak
  Configurable<float> maxMassNSigma{"maxMassNSigma", 12.0f, "max mass region to be considered for further analysis"};

  // For extracting strangeness mass QA plots
  ConfigurableAxis axisPtQA{"axisPtQA", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for QA histograms"};
  ConfigurableAxis axisK0ShortMass{"axisK0ShortMass", {200, 0.400f, 0.600f}, "Inv. Mass (GeV/c^{2})"};
  ConfigurableAxis axisLambdaMass{"axisLambdaMass", {200, 1.01f, 1.21f}, "Inv. Mass (GeV/c^{2})"};
  ConfigurableAxis axisXiMass{"axisXiMass", {200, 1.22f, 1.42f}, "Inv. Mass (GeV/c^{2})"};
  ConfigurableAxis axisOmegaMass{"axisOmegaMass", {200, 1.57f, 1.77f}, "Inv. Mass (GeV/c^{2})"};
  ConfigurableAxis axisMult{"axisMult", {VARIABLE_WIDTH, 0.0f, 0.01f, 1.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 70.0f, 100.0f}, "Centrality percentile bins"};

  // QA
  Configurable<bool> doTrueSelectionInMass{"doTrueSelectionInMass", false, "Fill mass histograms only with true primary Particles for MC"};
  // Do declarative selections for DCAs, if possible
  Filter preFilterTracks = nabs(aod::track::dcaXY) < dcaXYconstant + dcaXYpTdep * nabs(aod::track::signed1Pt);
  Filter preFilterV0 = nabs(aod::v0data::dcapostopv) > dcaPostopv&&
                                                         nabs(aod::v0data::dcanegtopv) > dcaNegtopv&& aod::v0data::dcaV0daughters < dcaV0dau;
  Filter preFilterCascade =
    nabs(aod::cascdata::dcapostopv) > dcaPostopv&& nabs(aod::cascdata::dcanegtopv) > dcaNegtopv&& nabs(aod::cascdata::dcabachtopv) > cascadesetting_dcabachtopv&& aod::cascdata::dcaV0daughters < dcaV0dau&& aod::cascdata::dcacascdaughters < cascadesetting_dcacascdau;

  // histogram defined with HistogramRegistry
  HistogramRegistry registry{
    "registry",
    {}};
  using V0LinkedTagged = soa::Join<aod::V0sLinked, aod::V0Tags>;
  using CascadesLinkedTagged = soa::Join<aod::CascadesLinked, aod::CascTags>;
  using DauTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::TracksDCA>;
  using DauTracksMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::TracksDCA, aod::McTrackLabels>;
  // using IDTracks= soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidBayesPi, aod::pidBayesKa, aod::pidBayesPr, aod::TOFSignal>; // prepared for Bayesian PID
  using IDTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::TOFSignal, aod::TracksDCA>;
  using V0DatasWithoutTrackX = soa::Join<aod::V0Indices, aod::V0Cores>;

  Produces<aod::TriggerTracks> triggerTrack;
  Produces<aod::AssocPions> assocPion;
  Produces<aod::AssocV0s> assocV0;
  Produces<aod::AssocCascades> assocCascades;

  TF1* fK0Mean = new TF1("fK0Mean", "[0]+[1]*x+[2]*TMath::Exp(-[3]*x)");
  TF1* fK0Width = new TF1("fK0Width", "[0]+[1]*x+[2]*TMath::Exp(-[3]*x)");
  TF1* fLambdaMean = new TF1("fLambdaMean", "[0]+[1]*x+[2]*TMath::Exp(-[3]*x)");
  TF1* fLambdaWidth = new TF1("fLambdaWidth", "[0]+[1]*x+[2]*TMath::Exp(-[3]*x)");
  TF1* fXiMean = new TF1("fXiMean", "[0]+[1]*x+[2]*TMath::Exp(-[3]*x)");
  TF1* fXiWidth = new TF1("fXiWidth", "[0]+[1]*x+[2]*TMath::Exp(-[3]*x)");
  TF1* fOmegaMean = new TF1("fomegaMean", "[0]+[1]*x+[2]*TMath::Exp(-[3]*x)");
  TF1* fOmegaWidth = new TF1("fomegaWidth", "[0]+[1]*x+[2]*TMath::Exp(-[3]*x)");

  void init(InitContext const&)
  {
    fK0Mean->SetParameters(massParsK0Mean->at(0), massParsK0Mean->at(1), massParsK0Mean->at(2), massParsK0Mean->at(3));
    fK0Width->SetParameters(massParsK0Width->at(0), massParsK0Width->at(1), massParsK0Width->at(2), massParsK0Width->at(3));
    fLambdaMean->SetParameters(massParsLambdaMean->at(0), massParsLambdaMean->at(1), massParsLambdaMean->at(2), massParsLambdaMean->at(3));
    fLambdaWidth->SetParameters(massParsLambdaWidth->at(0), massParsLambdaWidth->at(1), massParsLambdaWidth->at(2), massParsLambdaWidth->at(3));
    fXiMean->SetParameters(massParsCascadeMean->at(0), massParsCascadeMean->at(1), massParsCascadeMean->at(2), massParsCascadeMean->at(3));
    fXiWidth->SetParameters(massParsCascadeWidth->at(0), massParsCascadeWidth->at(1), massParsCascadeWidth->at(2), massParsCascadeWidth->at(3));
    fOmegaMean->SetParameters(massParsOmegaMean->at(0), massParsOmegaMean->at(1), massParsOmegaMean->at(2), massParsOmegaMean->at(3));
    fOmegaWidth->SetParameters(massParsOmegaWidth->at(0), massParsOmegaWidth->at(1), massParsOmegaWidth->at(2), massParsOmegaWidth->at(3));

    histos.add("h3dMassK0Short", "h3dMassK0Short", kTH3F, {axisPtQA, axisK0ShortMass, axisMult});
    histos.add("h3dMassLambda", "h3dMassLambda", kTH3F, {axisPtQA, axisLambdaMass, axisMult});
    histos.add("h3dMassAntiLambda", "h3dMassAntiLambda", kTH3F, {axisPtQA, axisLambdaMass, axisMult});
    histos.add("h3dMassXiMinus", "h3dMassXiMinus", kTH3F, {axisPtQA, axisXiMass, axisMult});
    histos.add("h3dMassXiPlus", "h3dMassXiPlus", kTH3F, {axisPtQA, axisXiMass, axisMult});
    histos.add("h3dMassOmegaMinus", "h3dMassOmegaMinus", kTH3F, {axisPtQA, axisOmegaMass, axisMult});
    histos.add("h3dMassOmegaPlus", "h3dMassOmegaPlus", kTH3F, {axisPtQA, axisOmegaMass, axisMult});
  }

  // reco-level trigger quality checks (N.B.: DCA is filtered, not selected)
  template <class TTrack>
  bool isValidTrigger(TTrack track)
  {
    if (track.eta() > triggerEtaMax || track.eta() < triggerEtaMin) {
      return false;
    }
    // if (track.sign()= 1 ) {continue;}
    if (track.pt() > triggerPtCutMax || track.pt() < triggerPtCutMin) {
      return false;
    }
    if (track.tpcNClsCrossedRows() < minTPCNCrossedRows) {
      return false; // crossed rows
    }
    if (!track.hasITS() && triggerRequireITS) {
      return false; // skip, doesn't have ITS signal (skips lots of TPC-only!)
    }
    if (track.tpcNClsShared() > triggerMaxTPCSharedClusters) {
      return false; // skip, has shared clusters
    }
    if (!(bitcheck(track.itsClusterMap(), 0)) && triggerRequireL0) {
      return false; // skip, doesn't have cluster in ITS L0
    }
    return true;
  }

  // for real data processing
  void processTriggers(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<DauTracks> const& tracks)
  {
    // Perform basic event selection
    if (!collision.sel8()) {
      return;
    }
    // No need to correlate stuff that's in far collisions
    if (TMath::Abs(collision.posZ()) > 10.0) {
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
        track.globalIndex());
    }
  }

  // for MC processing
  void processTriggersMC(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<DauTracksMC> const& tracks, aod::McParticles const&)
  {
    // Perform basic event selection
    if (!collision.sel8()) {
      return;
    }
    // No need to correlate stuff that's in far collisions
    if (TMath::Abs(collision.posZ()) > 10.0) {
      return;
    }

    /// _________________________________________________
    /// Step 1: Populate table with trigger tracks
    for (auto const& track : tracks) {
      if (!isValidTrigger(track))
        continue;
      bool physicalPrimary = false;
      if (track.has_mcParticle()) {
        auto mcParticle = track.mcParticle();
        physicalPrimary = mcParticle.isPhysicalPrimary();
      }
      triggerTrack(
        track.collisionId(),
        physicalPrimary,
        track.globalIndex());
    }
  }

  void processAssocPions(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<IDTracks> const& tracks)
  {
    // Perform basic event selection
    if (!collision.sel8()) {
      return;
    }
    // No need to correlate stuff that's in far collisions
    if (TMath::Abs(collision.posZ()) > 10.0) {
      return;
    }

    /// _________________________________________________
    /// Step 1: Populate table with trigger tracks
    for (auto const& track : tracks) {
      if (track.eta() > assocEtaMax || track.eta() < assocEtaMin) {
        continue;
      }
      // if (track.sign()= 1 ) {continue;}
      if (track.pt() > assocPtCutMax || track.pt() < assocPtCutMin) {
        continue;
      }
      if (track.tpcNClsCrossedRows() < minTPCNCrossedRows) {
        continue; // crossed rows
      }
      if (!track.hasITS() && triggerRequireITS) {
        continue; // skip, doesn't have ITS signal (skips lots of TPC-only!)
      }
      // prepared for Bayesian PID
      //  if (!track.bayesPi() > pionMinBayesProb) {
      //    continue;
      //  }
      //  if (track.bayesPi() < track.bayesPr() || track.bayesPi() < track.bayesKa()){
      //    continue;
      //  }
      //  if (track.tpcNSigmaPi() < assocPionNSigmaTPCFOF){
      //    continue;
      //  }
      //  if (track.tofSignal() > 0 && track.tofNSigmaPi() < assocPionNSigmaTPCFOF){
      //    continue;
      //  }
      if (track.tofSignal() > 0) {
        if (std::sqrt(track.tofNSigmaPi() * track.tofNSigmaPi() + track.tpcNSigmaPi() * track.tpcNSigmaPi()) > assocPionNSigmaTPCFOF)
          continue;
        if (track.tofNSigmaPr() < rejectSigma)
          continue;
        if (track.tpcNSigmaPr() < rejectSigma)
          continue;
        if (track.tofNSigmaKa() < rejectSigma)
          continue;
        if (track.tpcNSigmaKa() < rejectSigma)
          continue;
      } else {
        if (track.tpcNSigmaPi() > assocPionNSigmaTPCFOF)
          continue;
        if (track.tpcNSigmaPr() < rejectSigma)
          continue;
        if (track.tpcNSigmaKa() < rejectSigma)
          continue;
      }

      assocPion(
        track.collisionId(),
        track.globalIndex());
    }
  }

  void processV0s(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>::iterator const& collision, DauTracks const&, soa::Filtered<V0DatasWithoutTrackX> const& V0s, V0LinkedTagged const&)
  {
    // Perform basic event selection
    if (!collision.sel8()) {
      return;
    }
    // No need to correlate stuff that's in far collisions
    if (TMath::Abs(collision.posZ()) > 10.0) {
      return;
    }
    /// _________________________________________________
    /// Populate table with associated V0s
    for (auto const& v0 : V0s) {
      if (v0.v0radius() < v0RadiusMin || v0.v0radius() > v0RadiusMax || v0.eta() > assocEtaMax || v0.eta() < assocEtaMin || v0.v0cosPA() < v0Cospa) {
        continue;
      }
      // check dE/dx compatibility
      int compatibleK0Short = 0;
      int compatibleLambda = 0;
      int compatibleAntiLambda = 0;

      auto posdau = v0.posTrack_as<DauTracks>();
      auto negdau = v0.negTrack_as<DauTracks>();
      auto origV0entry = v0.v0_as<V0LinkedTagged>(); // retrieve tags

      if (negdau.tpcNClsCrossedRows() < minTPCNCrossedRows)
        continue;
      if (posdau.tpcNClsCrossedRows() < minTPCNCrossedRows)
        continue;

      if (TMath::Abs(posdau.tpcNSigmaPi()) < strangedEdxNSigmaLoose && TMath::Abs(negdau.tpcNSigmaPi()) < strangedEdxNSigmaLoose)
        bitset(compatibleK0Short, 0);
      if (TMath::Abs(posdau.tpcNSigmaPi()) < strangedEdxNSigma && TMath::Abs(negdau.tpcNSigmaPi()) < strangedEdxNSigma)
        bitset(compatibleK0Short, 1);
      if (TMath::Abs(posdau.tpcNSigmaPi()) < strangedEdxNSigmaTight && TMath::Abs(negdau.tpcNSigmaPi()) < strangedEdxNSigmaTight)
        bitset(compatibleK0Short, 2);

      if (TMath::Abs(posdau.tpcNSigmaPr()) < strangedEdxNSigmaLoose && TMath::Abs(negdau.tpcNSigmaPi()) < strangedEdxNSigmaLoose)
        if (v0.v0cosPA() > lambdaCospa)
          bitset(compatibleLambda, 0);
      if (TMath::Abs(posdau.tpcNSigmaPr()) < strangedEdxNSigma && TMath::Abs(negdau.tpcNSigmaPi()) < strangedEdxNSigma)
        if (v0.v0cosPA() > lambdaCospa)
          bitset(compatibleLambda, 1);
      if (TMath::Abs(posdau.tpcNSigmaPr()) < strangedEdxNSigmaTight && TMath::Abs(negdau.tpcNSigmaPi()) < strangedEdxNSigmaTight)
        if (v0.v0cosPA() > lambdaCospa)
          bitset(compatibleLambda, 2);

      if (TMath::Abs(posdau.tpcNSigmaPi()) < strangedEdxNSigmaLoose && TMath::Abs(negdau.tpcNSigmaPr()) < strangedEdxNSigmaLoose)
        if (v0.v0cosPA() > lambdaCospa)
          bitset(compatibleAntiLambda, 0);
      if (TMath::Abs(posdau.tpcNSigmaPi()) < strangedEdxNSigma && TMath::Abs(negdau.tpcNSigmaPr()) < strangedEdxNSigma)
        if (v0.v0cosPA() > lambdaCospa)
          bitset(compatibleAntiLambda, 1);
      if (TMath::Abs(posdau.tpcNSigmaPi()) < strangedEdxNSigmaTight && TMath::Abs(negdau.tpcNSigmaPr()) < strangedEdxNSigmaTight)
        if (v0.v0cosPA() > lambdaCospa)
          bitset(compatibleAntiLambda, 2);

      // simplified handling: calculate NSigma in mass here
      float massNSigmaK0Short = (v0.mK0Short() - fK0Mean->Eval(v0.pt())) / (fK0Width->Eval(v0.pt()) + 1e-6);
      float massNSigmaLambda = (v0.mLambda() - fLambdaMean->Eval(v0.pt())) / (fLambdaWidth->Eval(v0.pt()) + 1e-6);
      float massNSigmaAntiLambda = (v0.mAntiLambda() - fLambdaMean->Eval(v0.pt())) / (fLambdaWidth->Eval(v0.pt()) + 1e-6);

      if (compatibleK0Short && (!doTrueSelectionInMass || (origV0entry.isTrueK0Short() && origV0entry.isPhysicalPrimary())))
        histos.fill(HIST("h3dMassK0Short"), v0.pt(), v0.mK0Short(), collision.centFT0M());
      if (compatibleLambda && (!doTrueSelectionInMass || (origV0entry.isTrueLambda() && origV0entry.isPhysicalPrimary())))
        histos.fill(HIST("h3dMassLambda"), v0.pt(), v0.mLambda(), collision.centFT0M());
      if (compatibleAntiLambda && (!doTrueSelectionInMass || (origV0entry.isTrueAntiLambda() && origV0entry.isPhysicalPrimary())))
        histos.fill(HIST("h3dMassAntiLambda"), v0.pt(), v0.mAntiLambda(), collision.centFT0M());

      if (!fillTableOnlyWithCompatible ||
          ( // start major condition check
            (compatibleK0Short > 0 && std::abs(massNSigmaK0Short) < maxMassNSigma) ||
            (compatibleLambda > 0 && std::abs(massNSigmaLambda) < maxMassNSigma) ||
            (compatibleAntiLambda > 0 && std::abs(massNSigmaAntiLambda) < maxMassNSigma)) // end major condition check
      ) {
        assocV0(v0.collisionId(), v0.globalIndex(),
                compatibleK0Short, compatibleLambda, compatibleAntiLambda,
                origV0entry.isTrueK0Short(), origV0entry.isTrueLambda(), origV0entry.isTrueAntiLambda(), origV0entry.isPhysicalPrimary(),
                massNSigmaK0Short, massNSigmaLambda, massNSigmaAntiLambda);
      }
    }
  }
  void processCascades(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>::iterator const& collision, DauTracks const&, soa::Filtered<V0DatasWithoutTrackX> const& /*V0s*/, soa::Filtered<aod::CascDatas> const& Cascades, aod::V0sLinked const&, CascadesLinkedTagged const&)
  {
    // Perform basic event selection
    if (!collision.sel8()) {
      return;
    }
    // No need to correlate stuff that's in far collisions
    if (TMath::Abs(collision.posZ()) > 10.0) {
      return;
    }
    /// _________________________________________________
    /// Step 3: Populate table with associated Cascades
    for (auto const& casc : Cascades) {
      auto bachTrackCast = casc.bachelor_as<DauTracks>();
      auto posTrackCast = casc.posTrack_as<DauTracks>();
      auto negTrackCast = casc.negTrack_as<DauTracks>();
      auto origCascadeEntry = casc.cascade_as<CascadesLinkedTagged>();

      // minimum TPC crossed rows
      if (bachTrackCast.tpcNClsCrossedRows() < minTPCNCrossedRows)
        continue;
      if (posTrackCast.tpcNClsCrossedRows() < minTPCNCrossedRows)
        continue;
      if (negTrackCast.tpcNClsCrossedRows() < minTPCNCrossedRows)
        continue;

      // check dE/dx compatibility
      int compatibleXiMinus = 0;
      int compatibleXiPlus = 0;
      int compatibleOmegaMinus = 0;
      int compatibleOmegaPlus = 0;

      if (TMath::Abs(posTrackCast.tpcNSigmaPr()) < strangedEdxNSigmaLoose && TMath::Abs(negTrackCast.tpcNSigmaPi()) < strangedEdxNSigmaLoose && TMath::Abs(bachTrackCast.tpcNSigmaPi()) < strangedEdxNSigmaLoose && casc.sign() < 0)
        bitset(compatibleXiMinus, 0);
      if (TMath::Abs(posTrackCast.tpcNSigmaPr()) < strangedEdxNSigma && TMath::Abs(negTrackCast.tpcNSigmaPi()) < strangedEdxNSigma && TMath::Abs(bachTrackCast.tpcNSigmaPi()) < strangedEdxNSigma && casc.sign() < 0)
        bitset(compatibleXiMinus, 1);
      if (TMath::Abs(posTrackCast.tpcNSigmaPr()) < strangedEdxNSigmaTight && TMath::Abs(negTrackCast.tpcNSigmaPi()) < strangedEdxNSigmaTight && TMath::Abs(bachTrackCast.tpcNSigmaPi()) < strangedEdxNSigmaTight && casc.sign() < 0)
        bitset(compatibleXiMinus, 2);

      if (TMath::Abs(posTrackCast.tpcNSigmaPi()) < strangedEdxNSigmaLoose && TMath::Abs(negTrackCast.tpcNSigmaPr()) < strangedEdxNSigmaLoose && TMath::Abs(bachTrackCast.tpcNSigmaPi()) < strangedEdxNSigmaLoose && casc.sign() > 0)
        bitset(compatibleXiPlus, 0);
      if (TMath::Abs(posTrackCast.tpcNSigmaPi()) < strangedEdxNSigma && TMath::Abs(negTrackCast.tpcNSigmaPr()) < strangedEdxNSigma && TMath::Abs(bachTrackCast.tpcNSigmaPi()) < strangedEdxNSigma && casc.sign() > 0)
        bitset(compatibleXiPlus, 1);
      if (TMath::Abs(posTrackCast.tpcNSigmaPi()) < strangedEdxNSigmaTight && TMath::Abs(negTrackCast.tpcNSigmaPr()) < strangedEdxNSigmaTight && TMath::Abs(bachTrackCast.tpcNSigmaPi()) < strangedEdxNSigmaTight && casc.sign() > 0)
        bitset(compatibleXiPlus, 2);

      if (TMath::Abs(posTrackCast.tpcNSigmaPr()) < strangedEdxNSigmaLoose && TMath::Abs(negTrackCast.tpcNSigmaPi()) < strangedEdxNSigmaLoose && TMath::Abs(bachTrackCast.tpcNSigmaKa()) < strangedEdxNSigmaLoose && casc.sign() < 0)
        bitset(compatibleOmegaMinus, 0);
      if (TMath::Abs(posTrackCast.tpcNSigmaPr()) < strangedEdxNSigma && TMath::Abs(negTrackCast.tpcNSigmaPi()) < strangedEdxNSigma && TMath::Abs(bachTrackCast.tpcNSigmaKa()) < strangedEdxNSigma && casc.sign() < 0)
        bitset(compatibleOmegaMinus, 1);
      if (TMath::Abs(posTrackCast.tpcNSigmaPr()) < strangedEdxNSigmaTight && TMath::Abs(negTrackCast.tpcNSigmaPi()) < strangedEdxNSigmaTight && TMath::Abs(bachTrackCast.tpcNSigmaKa()) < strangedEdxNSigmaTight && casc.sign() < 0)
        bitset(compatibleOmegaMinus, 2);

      if (TMath::Abs(posTrackCast.tpcNSigmaPi()) < strangedEdxNSigmaLoose && TMath::Abs(negTrackCast.tpcNSigmaPr()) < strangedEdxNSigmaLoose && TMath::Abs(bachTrackCast.tpcNSigmaKa()) < strangedEdxNSigmaLoose && casc.sign() > 0)
        bitset(compatibleOmegaPlus, 0);
      if (TMath::Abs(posTrackCast.tpcNSigmaPi()) < strangedEdxNSigma && TMath::Abs(negTrackCast.tpcNSigmaPr()) < strangedEdxNSigma && TMath::Abs(bachTrackCast.tpcNSigmaKa()) < strangedEdxNSigma && casc.sign() > 0)
        bitset(compatibleOmegaPlus, 1);
      if (TMath::Abs(posTrackCast.tpcNSigmaPi()) < strangedEdxNSigmaTight && TMath::Abs(negTrackCast.tpcNSigmaPr()) < strangedEdxNSigmaTight && TMath::Abs(bachTrackCast.tpcNSigmaKa()) < strangedEdxNSigmaTight && casc.sign() > 0)
        bitset(compatibleOmegaPlus, 2);

      float massNSigmaXi = (casc.mXi() - fXiMean->Eval(casc.pt())) / (fXiWidth->Eval(casc.pt()) + 1e-6);
      float massNSigmaOmega = (casc.mOmega() - fOmegaMean->Eval(casc.pt())) / (fOmegaWidth->Eval(casc.pt()) + 1e-6);

      if (compatibleXiMinus && (!doTrueSelectionInMass || (origCascadeEntry.isTrueXiMinus() && origCascadeEntry.isPhysicalPrimary())))
        histos.fill(HIST("h3dMassXiMinus"), casc.pt(), casc.mXi(), collision.centFT0M());
      if (compatibleXiPlus && (!doTrueSelectionInMass || (origCascadeEntry.isTrueXiPlus() && origCascadeEntry.isPhysicalPrimary())))
        histos.fill(HIST("h3dMassXiPlus"), casc.pt(), casc.mXi(), collision.centFT0M());
      if (compatibleOmegaMinus && (!doTrueSelectionInMass || (origCascadeEntry.isTrueOmegaMinus() && origCascadeEntry.isPhysicalPrimary())))
        histos.fill(HIST("h3dMassOmegaMinus"), casc.pt(), casc.mOmega(), collision.centFT0M());
      if (compatibleOmegaPlus && (!doTrueSelectionInMass || (origCascadeEntry.isTrueOmegaPlus() && origCascadeEntry.isPhysicalPrimary())))
        histos.fill(HIST("h3dMassOmegaPlus"), casc.pt(), casc.mOmega(), collision.centFT0M());

      if (!fillTableOnlyWithCompatible ||
          ( // start major condition check
            ((compatibleXiMinus > 0 || compatibleXiPlus > 0) && std::abs(massNSigmaXi) < maxMassNSigma) ||
            ((compatibleOmegaMinus > 0 || compatibleOmegaPlus > 0) && std::abs(massNSigmaOmega) < maxMassNSigma)) // end major condition check
      ) {
        assocCascades(casc.collisionId(), casc.globalIndex(),
                      compatibleXiMinus, compatibleXiPlus, compatibleOmegaMinus, compatibleOmegaPlus,
                      origCascadeEntry.isTrueXiMinus(), origCascadeEntry.isTrueXiPlus(),
                      origCascadeEntry.isTrueOmegaMinus(), origCascadeEntry.isTrueOmegaPlus(),
                      origCascadeEntry.isPhysicalPrimary(),
                      massNSigmaXi, massNSigmaOmega);
      }
    }
  }

  PROCESS_SWITCH(hstrangecorrelationfilter, processTriggers, "Produce trigger tables", true);
  PROCESS_SWITCH(hstrangecorrelationfilter, processTriggersMC, "Produce trigger tables for MC", false);
  PROCESS_SWITCH(hstrangecorrelationfilter, processV0s, "Produce associated V0 tables", true);
  PROCESS_SWITCH(hstrangecorrelationfilter, processAssocPions, "Produce associated Pion tables", true);
  PROCESS_SWITCH(hstrangecorrelationfilter, processCascades, "Produce associated cascade tables", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<hstrangecorrelationfilter>(cfgc)};
}
