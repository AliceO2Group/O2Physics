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
#include "TF1.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct hstrangecorrelationfilter {
  // Trigger particle selections in phase space
  Configurable<float> triggerEtaMin{"triggerEtaCutMin", -0.8, "triggeretamin"};
  Configurable<float> triggerEtaMax{"triggerEtaCutMax", 0.8, "triggeretamax"};
  Configurable<float> triggerPtCutMin{"triggerPtCutMin", 3, "triggerptmin"};
  Configurable<float> triggerPtCutMax{"triggerPtCutMax", 20, "triggerptmax"};

  // Track quality
  Configurable<int> minTPCNCrossedRows{"minTPCNCrossedRows", 50, "Minimum TPC crossed rows"};
  Configurable<bool> triggerRequireITS{"triggerRequireITS", true, "require ITS signal in trigger tracks"};

  // Associated particle selections in phase space
  Configurable<float> assocEtaMin{"assocEtaCutMin", -0.8, "triggeretamin"};
  Configurable<float> assocEtaMax{"assocEtaCutMax", 0.8, "triggeretamax"};
  Configurable<float> assocPtCutMin{"assocPtCutMin", 0.2, "assocptmin"};
  Configurable<float> assocPtCutMax{"assocPtCutMax", 10, "assocptmax"};

  // Associated pion identification
  Configurable<float> pionMinBayesProb{"pionMinBayesProb", 0.95, "minimal Bayesian probablity for pion ID"};
  Configurable<float> assocPionNSigmaTPCFOF{"assocPionNSigmaTPCFOF", 3, "minimal n sigma in TOF and TPC for Pion ID"};
  Configurable<float> rejectSigma{"rejectSigma", 1, "n sigma for rejecting pion candidates"};

  // V0 selections
  Configurable<double> v0Cospa{"v0cospa", 0.97, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<float> dcaV0dau{"dcav0dau", 1.0, "DCA V0 Daughters"};
  Configurable<float> dcaNegtopv{"dcanegtopv", 0.06, "DCA Neg To PV"};
  Configurable<float> dcaPostopv{"dcapostopv", 0.06, "DCA Pos To PV"};
  Configurable<float> v0RadiusMin{"v0radiusmin", 0.5, "v0radius"};
  Configurable<float> v0RadiusMax{"v0radiusmax", 200, "v0radius"};

  // primary particle DCAxy selections
  // formula: |DCAxy| <  0.004f + (0.013f / pt)
  Configurable<std::vector<float>> dcaXYpars{"dcaXYpars", {0.004, 0.013}, "pars in |DCAxy| < [0]+[1]/pT"};
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


  // definitions of peak and background
  Configurable<float> peakNsigma{"peakNsigma", 3.0f, "peak region is +/- this many sigmas away"};
  Configurable<float> backgroundNsigma{"backgroundNsigma", 6.0f, "bg region is +/- this many sigmas away (minus peak)"};

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
  using DauTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::TracksDCA>;
  // using IDTracks= soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidBayesPi, aod::pidBayesKa, aod::pidBayesPr, aod::TOFSignal>; // prepared for Bayesian PID
  using IDTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::TOFSignal, aod::TracksDCA>;

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
  }

  void processTriggers(soa::Join<aod::Collisions, aod::EvSels, aod::Mults>::iterator const& collision, soa::Filtered<DauTracks> const& tracks)
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
      if (track.eta() > triggerEtaMax || track.eta() < triggerEtaMin) {
        continue;
      }
      // if (track.sign()= 1 ) {continue;}
      if (track.pt() > triggerPtCutMax || track.pt() < triggerPtCutMin) {
        continue;
      }
      if (track.tpcNClsCrossedRows() < minTPCNCrossedRows) {
        continue; // crossed rows
      }
      if (!track.hasITS() && triggerRequireITS) {
        continue; // skip, doesn't have ITS signal (skips lots of TPC-only!)
      }
      triggerTrack(
        track.collisionId(),
        track.globalIndex());
    }
  }
  void processAssocPions(soa::Join<aod::Collisions, aod::EvSels, aod::Mults>::iterator const& collision, soa::Filtered<IDTracks> const& tracks)
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

  void processV0s(soa::Join<aod::Collisions, aod::EvSels, aod::Mults>::iterator const& collision, DauTracks const&, soa::Filtered<aod::V0Datas> const& V0s, aod::V0sLinked const&)
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
      if (v0.v0radius() < v0RadiusMin || v0.v0radius() > v0RadiusMax || v0.eta() > assocEtaMax || v0.eta() < assocEtaMin || v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < v0Cospa) {
        continue;
      }
      // check dE/dx compatibility
      bool compatibleK0Short = false;
      bool compatibleLambda = false;
      bool compatibleAntiLambda = false;

      auto posdau = v0.posTrack_as<DauTracks>();
      auto negdau = v0.negTrack_as<DauTracks>();

      if (negdau.tpcNClsCrossedRows() < minTPCNCrossedRows)
        continue;
      if (posdau.tpcNClsCrossedRows() < minTPCNCrossedRows)
        continue;

      if (TMath::Abs(posdau.tpcNSigmaPi()) < 5 && TMath::Abs(negdau.tpcNSigmaPi()) < 5) {
        compatibleK0Short = true;
      }
      if (TMath::Abs(posdau.tpcNSigmaPr()) < 5 && TMath::Abs(negdau.tpcNSigmaPi()) < 5) {
        compatibleLambda = true;
      }
      if (TMath::Abs(posdau.tpcNSigmaPi()) < 5 && TMath::Abs(negdau.tpcNSigmaPr()) < 5) {
        compatibleAntiLambda = true;
      }
      // check whether V0s are in the regin
      int massRegK0Short = -1;
      if (TMath::Abs(v0.mK0Short() - fK0Mean->Eval(v0.pt()) < peakNsigma * fK0Width->Eval(v0.pt()))) {
        massRegK0Short = 2;
      }
      if (v0.mK0Short() > fK0Mean->Eval(v0.pt()) + peakNsigma * fK0Width->Eval(v0.pt()) && v0.mK0Short() < fK0Mean->Eval(v0.pt()) + backgroundNsigma * fK0Width->Eval(v0.pt())) {
        massRegK0Short = 3;
      }
      if (v0.mK0Short() < fK0Mean->Eval(v0.pt()) - peakNsigma * fK0Width->Eval(v0.pt()) && v0.mK0Short() > fK0Mean->Eval(v0.pt()) - backgroundNsigma * fK0Width->Eval(v0.pt())) {
        massRegK0Short = 1;
      }
      if (v0.mK0Short() > fK0Mean->Eval(v0.pt()) + backgroundNsigma * fK0Width->Eval(v0.pt())) {
        massRegK0Short = 4;
      }
      if (v0.mK0Short() < fK0Mean->Eval(v0.pt()) - backgroundNsigma * fK0Width->Eval(v0.pt())) {
        massRegK0Short = 0;
      }

      int massRegLambda = -1;
      if (TMath::Abs(v0.mLambda() - fLambdaMean->Eval(v0.pt()) < peakNsigma * fLambdaWidth->Eval(v0.pt()))) {
        massRegLambda = 2;
      }
      if (v0.mLambda() > fLambdaMean->Eval(v0.pt()) + peakNsigma * fLambdaWidth->Eval(v0.pt()) && v0.mLambda() < fLambdaMean->Eval(v0.pt()) + backgroundNsigma * fLambdaWidth->Eval(v0.pt())) {
        massRegLambda = 3;
      }
      if (v0.mLambda() < fLambdaMean->Eval(v0.pt()) - peakNsigma * fLambdaWidth->Eval(v0.pt()) && v0.mLambda() > fLambdaMean->Eval(v0.pt()) - backgroundNsigma * fLambdaWidth->Eval(v0.pt())) {
        massRegLambda = 1;
      }
      if (v0.mLambda() > fLambdaMean->Eval(v0.pt()) + backgroundNsigma * fLambdaWidth->Eval(v0.pt())) {
        massRegLambda = 4;
      }
      if (v0.mLambda() < fLambdaMean->Eval(v0.pt()) - backgroundNsigma * fLambdaWidth->Eval(v0.pt())) {
        massRegLambda = 0;
      }

      int massRegAntiLambda = -1;
      if (TMath::Abs(v0.mAntiLambda() - fLambdaMean->Eval(v0.pt()) < peakNsigma * fLambdaWidth->Eval(v0.pt()))) {
        massRegAntiLambda = 2;
      }
      if (v0.mAntiLambda() > fLambdaMean->Eval(v0.pt()) + peakNsigma * fLambdaWidth->Eval(v0.pt()) && v0.mAntiLambda() < fLambdaMean->Eval(v0.pt()) + backgroundNsigma * fLambdaWidth->Eval(v0.pt())) {
        massRegAntiLambda = 3;
      }
      if (v0.mAntiLambda() < fLambdaMean->Eval(v0.pt()) - peakNsigma * fLambdaWidth->Eval(v0.pt()) && v0.mAntiLambda() > fLambdaMean->Eval(v0.pt()) - backgroundNsigma * fLambdaWidth->Eval(v0.pt())) {
        massRegAntiLambda = 1;
      }
      if (v0.mAntiLambda() > fLambdaMean->Eval(v0.pt()) + backgroundNsigma * fLambdaWidth->Eval(v0.pt())) {
        massRegAntiLambda = 4;
      }
      if (v0.mAntiLambda() < fLambdaMean->Eval(v0.pt()) - backgroundNsigma * fLambdaWidth->Eval(v0.pt())) {
        massRegAntiLambda = 0;
      }
      assocV0(v0.collisionId(), v0.globalIndex(), compatibleK0Short, compatibleLambda, compatibleAntiLambda, massRegK0Short, massRegLambda, massRegAntiLambda);
    }
  }
  void processCascades(soa::Join<aod::Collisions, aod::EvSels, aod::Mults>::iterator const& collision, DauTracks const&, soa::Filtered<aod::V0Datas> const& V0s, soa::Filtered<aod::CascDatas> const& Cascades, aod::V0sLinked const&)
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
      auto v0 = casc.v0_as<o2::aod::V0sLinked>();
      if (!(v0.has_v0Data())) {
        return; // skip those cascades for which V0 doesn't exist
      }
      auto v0data = v0.v0Data(); // de-reference index to correct v0data in case it exists
      auto bachTrackCast = casc.bachelor_as<DauTracks>();
      auto posTrackCast = v0data.posTrack_as<DauTracks>();
      auto negTrackCast = v0data.negTrack_as<DauTracks>();

      // minimum TPC crossed rows
      if (bachTrackCast.tpcNClsCrossedRows() < minTPCNCrossedRows)
        continue;
      if (posTrackCast.tpcNClsCrossedRows() < minTPCNCrossedRows)
        continue;
      if (negTrackCast.tpcNClsCrossedRows() < minTPCNCrossedRows)
        continue;

      // check dE/dx compatibility
      bool compatibleXiMinus = false;
      bool compatibleXiPlus = false;
      bool compatibleOmegaMinus = false;
      bool compatibleOmegaPlus = false;

      if (TMath::Abs(posTrackCast.tpcNSigmaPr()) < 5 && TMath::Abs(negTrackCast.tpcNSigmaPi()) < 5 && TMath::Abs(bachTrackCast.tpcNSigmaPi()) < 5) {
        compatibleXiMinus = true;
      }
      if (TMath::Abs(posTrackCast.tpcNSigmaPi()) < 5 && TMath::Abs(negTrackCast.tpcNSigmaPr()) < 5 && TMath::Abs(bachTrackCast.tpcNSigmaPi()) < 5) {
        compatibleXiPlus = true;
      }
      if (TMath::Abs(posTrackCast.tpcNSigmaPr()) < 5 && TMath::Abs(negTrackCast.tpcNSigmaPi()) < 5 && TMath::Abs(bachTrackCast.tpcNSigmaKa()) < 5) {
        compatibleOmegaMinus = true;
      }
      if (TMath::Abs(posTrackCast.tpcNSigmaPi()) < 5 && TMath::Abs(negTrackCast.tpcNSigmaPr()) < 5 && TMath::Abs(bachTrackCast.tpcNSigmaKa()) < 5) {
        compatibleOmegaPlus = true;
      }

      int massRegXi = -1;
      if (TMath::Abs(casc.mXi() - fXiMean->Eval(casc.pt()) < peakNsigma * fXiWidth->Eval(casc.pt()))) {
        massRegXi = 2;
      }
      if (casc.mXi() > fXiMean->Eval(casc.pt()) + peakNsigma * fXiWidth->Eval(casc.pt()) && casc.mXi() < fXiMean->Eval(casc.pt()) + backgroundNsigma * fXiWidth->Eval(casc.pt())) {
        massRegXi = 3;
      }
      if (casc.mXi() < fXiMean->Eval(casc.pt()) - peakNsigma * fXiWidth->Eval(casc.pt()) && casc.mXi() > fXiMean->Eval(casc.pt()) - backgroundNsigma * fXiWidth->Eval(casc.pt())) {
        massRegXi = 1;
      }
      if (casc.mXi() > fXiMean->Eval(casc.pt()) + backgroundNsigma * fXiWidth->Eval(casc.pt())) {
        massRegXi = 4;
      }
      if (casc.mXi() < fXiMean->Eval(casc.pt()) - backgroundNsigma * fXiWidth->Eval(casc.pt())) {
        massRegXi = 0;
      }

      int massRegOmega = -1;
      if (TMath::Abs(casc.mOmega() - fOmegaMean->Eval(casc.pt()) < peakNsigma * fOmegaWidth->Eval(casc.pt()))) {
        massRegOmega = 2;
      }
      if (casc.mOmega() > fOmegaMean->Eval(casc.pt()) + peakNsigma * fOmegaWidth->Eval(casc.pt()) && casc.mOmega() < fOmegaMean->Eval(casc.pt()) + backgroundNsigma * fOmegaWidth->Eval(casc.pt())) {
        massRegOmega = 3;
      }
      if (casc.mOmega() < fOmegaMean->Eval(casc.pt()) - peakNsigma * fOmegaWidth->Eval(casc.pt()) && casc.mOmega() > fOmegaMean->Eval(casc.pt()) - backgroundNsigma * fOmegaWidth->Eval(casc.pt())) {
        massRegOmega = 1;
      }
      if (casc.mOmega() > fOmegaMean->Eval(casc.pt()) + backgroundNsigma * fOmegaWidth->Eval(casc.pt())) {
        massRegOmega = 4;
      }
      if (casc.mOmega() < fOmegaMean->Eval(casc.pt()) - backgroundNsigma * fOmegaWidth->Eval(casc.pt())) {
        massRegOmega = 0;
      }
      assocCascades(casc.collisionId(), casc.globalIndex(), compatibleXiMinus, compatibleXiPlus, compatibleOmegaMinus, compatibleOmegaPlus, massRegXi, massRegOmega);
    }
  }

  PROCESS_SWITCH(hstrangecorrelationfilter, processTriggers, "Produce trigger tables", true);
  PROCESS_SWITCH(hstrangecorrelationfilter, processV0s, "Produce associated V0 tables", true);
  PROCESS_SWITCH(hstrangecorrelationfilter, processAssocPions, "Produce associated Pion tables", true);
  PROCESS_SWITCH(hstrangecorrelationfilter, processCascades, "Produce associated cascade tables", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<hstrangecorrelationfilter>(cfgc)};
}
