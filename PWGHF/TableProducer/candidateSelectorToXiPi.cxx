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

/// \file candidateSelectorOmegac.cxx
/// \brief Xic and Omegac → Xi Pi selection task
/// \author Federica Zanone <federica.zanone@cern.ch>, Heidelberg University & GSI

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectorPID.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::analysis::pdg;

enum pidInfoStored {
  kPiFromLam = 0,
  kPrFromLam,
  kPiFromCasc,
  kPiFromCharm
};

/// Struct for applying Omegac selection cuts
struct HfCandidateSelectorToXiPi {
  Produces<aod::HfSelToXiPi> hfSelToXiPi;

  // LF analysis selections
  // zPV -> can be already set in HFeventselection -> 10 cm
  // sel8 -> can be already set in HFeventselection -> true
  Configurable<double> radiusCascMin{"radiusCascMin", 0.5, "Min cascade radius"};
  Configurable<double> radiusV0Min{"radiusV0Min", 2, "Min V0 radius"};
  Configurable<double> cosPAV0Min{"cosPAV0Min", 0.95, "Min valueCosPA V0"};
  Configurable<double> cosPACascMin{"cosPACascMin", 0.95, "Min value CosPA cascade"};
  Configurable<double> dcaCascDauMax{"dcaCascDauMax", 5.0, "Max DCA cascade daughters"};
  Configurable<double> dcaV0DauMax{"dcaV0DauMax", 5.0, "Max DCA V0 daughters"};

  Configurable<double> dcaOmegacDauMax{"dcaOmegacDauMax", 5.0, "Max DCA omegac daughters"};

  // limit charm baryon invariant mass spectrum
  Configurable<double> invMassOmegacMin{"invMassOmegacMin", 2.0, "Lower limit invariant mass spectrum charm baryon"}; // 2.4 Omegac0 only
  Configurable<double> invMassOmegacMax{"invMassOmegacMax", 3.1, "Upper limit invariant mass spectrum charm baryon"};

  // kinematic selections
  Configurable<double> etaTrackCharmBachMax{"etaTrackCharmBachMax", 0.8, "Max absolute value of eta for charm baryon bachelor"};
  Configurable<double> etaTrackLFDauMax{"etaTrackLFDauMax", 1.0, "Max absolute value of eta for V0 and cascade daughters"};
  Configurable<double> ptPiFromCascMin{"ptPiFromCascMin", 0.15, "Min pT pi <- casc"};
  Configurable<double> ptPiFromOmeMin{"ptPiFromOmeMin", 0.2, "Min pT pi <- omegac"};

  Configurable<double> impactParameterXYPriPiMin{"impactParameterXYPriPiMin", 0., "Min dcaxy primary pi track to PV"};
  Configurable<double> impactParameterXYPriPiMax{"impactParameterXYPriPiMax", 10., "Max dcaxy primary pi track to PV"};
  Configurable<double> impactParameterZPriPiMin{"impactParameterZPriPiMin", 0., "Min dcaz primary pi track to PV"};
  Configurable<double> impactParameterZPriPiMax{"impactParameterZPriPiMax", 10., "Max dcaz primary pi track to PV"};

  Configurable<double> impactParameterXYCascMin{"impactParameterXYCascMin", 0., "Min dcaxy cascade track to PV"};
  Configurable<double> impactParameterXYCascMax{"impactParameterXYCascMax", 10., "Max dcaxy cascade track to PV"};
  Configurable<double> impactParameterZCascMin{"impactParameterZCascMin", 0., "Min dcaz cascade track to PV"};
  Configurable<double> impactParameterZCascMax{"impactParameterZCascMax", 10., "Max dcaz cascade track to PV"};

  Configurable<double> ptCandMin{"ptCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> ptCandMax{"ptCandMax", 50., "Upper bound of candidate pT"};

  // PID options
  Configurable<bool> usePidTpcOnly{"usePidTpcOnly", true, "Perform PID using only TPC"};
  Configurable<bool> usePidTpcTofCombined{"usePidTpcTofCombined", false, "Perform PID using TPC & TOF"};

  // PID - TPC selections
  Configurable<double> ptPiPidTpcMin{"ptPiPidTpcMin", 0.15, "Lower bound of track pT for TPC PID for pion selection"};
  Configurable<double> ptPiPidTpcMax{"ptPiPidTpcMax", 5., "Upper bound of track pT for TPC PID for pion selection"};
  Configurable<double> nSigmaTpcPiMax{"nSigmaTpcPiMax", 3., "Nsigma cut on TPC only for pion selection"};
  Configurable<double> nSigmaTpcCombinedPiMax{"nSigmaTpcCombinedPiMax", 5., "Nsigma cut on TPC combined with TOF for pion selection"};

  Configurable<double> ptPrPidTpcMin{"ptPrPidTpcMin", 0.15, "Lower bound of track pT for TPC PID for proton selection"};
  Configurable<double> ptPrPidTpcMax{"ptPrPidTpcMax", 5., "Upper bound of track pT for TPC PID for proton selection"};
  Configurable<double> nSigmaTpcPrMax{"nSigmaTpcPrMax", 3., "Nsigma cut on TPC only for proton selection"};
  Configurable<double> nSigmaTpcCombinedPrMax{"nSigmaTpcCombinedPrMax", 5., "Nsigma cut on TPC combined with TOF for proton selection"};

  // PID - TOF selections
  Configurable<double> ptPiPidTofMin{"ptPiPidTofMin", 0.15, "Lower bound of track pT for TOF PID for pion selection"};
  Configurable<double> ptPiPidTofMax{"ptPiPidTofMax", 5., "Upper bound of track pT for TOF PID for pion selection"};
  Configurable<double> nSigmaTofPiMax{"nSigmaTofPiMax", 3., "Nsigma cut on TOF only for pion selection"};
  Configurable<double> nSigmaTofCombinedPiMax{"nSigmaTofCombinedPiMax", 5., "Nsigma cut on TOF combined with TPC for pion selection"};

  Configurable<double> ptPrPidTofMin{"ptPrPidTofMin", 0.15, "Lower bound of track pT for TOF PID for proton selection"};
  Configurable<double> ptPrPidTofMax{"ptPrPidTofMax", 5., "Upper bound of track pT for TOF PID for proton selection"};
  Configurable<double> nSigmaTofPrMax{"nSigmaTofPrMax", 3., "Nsigma cut on TOF only for proton selection"};
  Configurable<double> nSigmaTofCombinedPrMax{"nSigmaTofCombinedPrMax", 5., "Nsigma cut on TOF combined with TPC for proton selection"};

  // invariant mass cuts
  Configurable<double> sigmaInvMassLambda{"sigmaInvMassLambda", 0.0025, "Invariant mass cut for lambda (sigma)"};
  Configurable<double> sigmaInvMassCascade{"sigmaInvMassCascade", 0.0025, "Invariant mass cut for cascade (sigma)"};
  Configurable<int> nSigmaInvMassCut{"nSigmaInvMassCut", 4, "Number of sigma for invariant mass cut"};

  // detector clusters selections
  Configurable<int> nClustersTpcMin{"nClustersTpcMin", 70, "Minimum number of TPC clusters requirement"};
  Configurable<int> nTpcCrossedRowsMin{"nTpcCrossedRowsMin", 70, "Minimum number of TPC crossed rows requirement"};
  Configurable<double> tpcCrossedRowsOverFindableClustersRatioMin{"tpcCrossedRowsOverFindableClustersRatioMin", 0.8, "Minimum ratio TPC crossed rows over findable clusters requirement"};
  Configurable<int> nClustersItsMin{"nClustersItsMin", 3, "Minimum number of ITS clusters requirement for pi <- Omegac"};
  Configurable<int> nClustersItsInnBarrMin{"nClustersItsInnBarrMin", 1, "Minimum number of ITS clusters in inner barrel requirement for pi <- Omegac"};

  TrackSelectorPi selectorPion;
  TrackSelectorPr selectorProton;

  using TracksSel = soa::Join<aod::TracksWDcaExtra, aod::TracksPidPi, aod::TracksPidPr>;

  HistogramRegistry registry{"registry"}; // for QA of selections

  OutputObj<TH1F> hInvMassOmegac{TH1F("hInvMassOmegac", "Omegac invariant mass;inv mass;entries", 500, 2.2, 3.1)};
  OutputObj<TH1F> hNEventsSaved{TH1F("hNEventsSaved", "Events with a charmed baryon candidate;Events source;N. events", 3, 0, 3)};

  void init(InitContext const&)
  {
    selectorPion.setRangePtTpc(ptPiPidTpcMin, ptPiPidTpcMax);
    selectorPion.setRangeNSigmaTpc(-nSigmaTpcPiMax, nSigmaTpcPiMax);
    selectorPion.setRangeNSigmaTpcCondTof(-nSigmaTpcCombinedPiMax, nSigmaTpcCombinedPiMax);
    selectorPion.setRangePtTof(ptPiPidTofMin, ptPiPidTofMax);
    selectorPion.setRangeNSigmaTof(-nSigmaTofPiMax, nSigmaTofPiMax);
    selectorPion.setRangeNSigmaTofCondTpc(-nSigmaTofCombinedPiMax, nSigmaTofCombinedPiMax);

    selectorProton.setRangePtTpc(ptPrPidTpcMin, ptPrPidTpcMax);
    selectorProton.setRangeNSigmaTpc(-nSigmaTpcPrMax, nSigmaTpcPrMax);
    selectorProton.setRangeNSigmaTpcCondTof(-nSigmaTpcCombinedPrMax, nSigmaTpcCombinedPrMax);
    selectorProton.setRangePtTof(ptPrPidTofMin, ptPrPidTofMax);
    selectorProton.setRangeNSigmaTof(-nSigmaTofPrMax, nSigmaTofPrMax);
    selectorProton.setRangeNSigmaTofCondTpc(-nSigmaTofCombinedPrMax, nSigmaTofCombinedPrMax);

    registry.add("hSelPID", "hSelPID;status;entries", {HistType::kTH1F, {{12, 0., 12.}}});
    registry.add("hTest", "Test status consecutive;status;entries", {HistType::kTH1F, {{12, 0., 12.}}});

    // for QA of the selections (bin 0 -> candidates that did not pass the selection, bin 1 -> candidates that passed the selection)
    registry.add("hSelSignDec", "hSelSignDec;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelEtaPosV0Dau", "hSelEtaPosV0Dau;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelEtaNegV0Dau", "hSelEtaNegV0Dau;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelEtaPiFromCasc", "hSelEtaPiFromCasc;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelEtaPiFromOme", "hSelEtaPiFromOme;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelRadCasc", "hSelRadCasc;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelRadV0", "hSelRadV0;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelCosPACasc", "hSelCosPACasc;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelCosPAV0", "hSelCosPAV0;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelDCACascDau", "hSelDCACascDau;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelDCAV0Dau", "hSelDCAV0Dau;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelDCAOmeDau", "hSelDCAOmeDau;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelDCAXYPrimPi", "hSelDCAXYPrimPi;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelDCAZPrimPi", "hSelDCAZPrimPi;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelDCAXYCasc", "hSelDCAXYCasc;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelDCAZCasc", "hSelDCAZCasc;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelPtPiFromCasc", "hSelPtPiFromCasc;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelPtPiFromOme", "hSelPtPiFromOme;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelNClsTPCPiFromOme", "hSelNClsTPCPiFromOme;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelNClsTPCPiFromLam", "hSelNClsTPCPiFromLam;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelNClsTPCPrFromLam", "hSelNClsTPCPrFromLam;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelNClsTPCPiFromCasc", "hSelNClsTPCPiFromCasc;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelNCrossRowsTPCPiFromOme", "hSelNCrossRowsTPCPiFromOme;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelNCrossRowsTPCPiFromLam", "hSelNCrossRowsTPCPiFromLam;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelNCrossRowsTPCPrFromLam", "hSelNCrossRowsTPCPrFromLam;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelNCrossRowsTPCPiFromCasc", "hSelNCrossRowsTPCPiFromCasc;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelCrossRowsOverFindClsTPCAllTracks", "hSelCrossRowsOverFindClsTPCAllTracks;status;entries", {HistType::kTH1F, {{10, 0., 10.}}});
    registry.add("hSelNClsITSPiFromOme", "hSelNClsITSPiFromOme;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelNClsITSInnerPiFromOme", "hSelNClsITSInnerPiFromOme;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelMassLam", "hSelMassLam;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelMassCasc", "hSelMassCasc;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelMassOme", "hSelMassOme;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
  }

  void process(aod::HfCandToXiPi const& candidates,
               TracksSel const&)
  {
    double massLambdaFromPDG = o2::analysis::pdg::MassLambda0;
    double massXiFromPDG = o2::analysis::pdg::MassXiMinus;

    int collId = -999;

    // looping over omegac candidates
    for (const auto& candidate : candidates) {

      bool resultSelections = true; // True if the candidate passes all the selections, False otherwise

      auto trackV0PosDau = candidate.posTrack_as<TracksSel>();    // positive V0 daughter
      auto trackV0NegDau = candidate.negTrack_as<TracksSel>();    // negative V0 daughter
      auto trackPiFromCasc = candidate.bachelor_as<TracksSel>();  // pion <- cascade
      auto trackPiFromOmeg = candidate.primaryPi_as<TracksSel>(); // pion <- omegac

      auto trackPiFromLam = trackV0NegDau;
      auto trackPrFromLam = trackV0PosDau;

      int8_t signDecay = candidate.signDecay(); // sign of pi <- cascade

      if (signDecay > 0) {
        trackPiFromLam = trackV0PosDau;
        trackPrFromLam = trackV0NegDau;
        registry.fill(HIST("hSelSignDec"), 1);
      } else if (signDecay == 0) {
        resultSelections = false;
        registry.fill(HIST("hSelSignDec"), 0);
      } else {
        registry.fill(HIST("hSelSignDec"), 1);
      }

      // eta selection
      double etaV0PosDau = candidate.etaV0PosDau();
      double etaV0NegDau = candidate.etaV0NegDau();
      double etaPiFromCasc = candidate.etaPiFromCasc();
      double etaPiFromOme = candidate.etaPiFromOme();
      if (std::abs(etaV0PosDau) > etaTrackLFDauMax) {
        resultSelections = false;
        registry.fill(HIST("hSelEtaPosV0Dau"), 0);
      } else {
        registry.fill(HIST("hSelEtaPosV0Dau"), 1);
      }
      if (std::abs(etaV0NegDau) > etaTrackLFDauMax) {
        resultSelections = false;
        registry.fill(HIST("hSelEtaNegV0Dau"), 0);
      } else {
        registry.fill(HIST("hSelEtaNegV0Dau"), 1);
      }
      if (std::abs(etaPiFromCasc) > etaTrackLFDauMax) {
        resultSelections = false;
        registry.fill(HIST("hSelEtaPiFromCasc"), 0);
      } else {
        registry.fill(HIST("hSelEtaPiFromCasc"), 1);
      }
      if (std::abs(etaPiFromOme) > etaTrackCharmBachMax) {
        resultSelections = false;
        registry.fill(HIST("hSelEtaPiFromOme"), 0);
      } else {
        registry.fill(HIST("hSelEtaPiFromOme"), 1);
      }

      // minimum radius cut (LFcut)
      if (RecoDecay::sqrtSumOfSquares(candidate.xDecayVtxCascade(), candidate.yDecayVtxCascade()) < radiusCascMin) {
        resultSelections = false;
        registry.fill(HIST("hSelRadCasc"), 0);
      } else {
        registry.fill(HIST("hSelRadCasc"), 1);
      }
      if (RecoDecay::sqrtSumOfSquares(candidate.xDecayVtxV0(), candidate.yDecayVtxV0()) < radiusV0Min) {
        resultSelections = false;
        registry.fill(HIST("hSelRadV0"), 0);
      } else {
        registry.fill(HIST("hSelRadV0"), 1);
      }

      // cosPA (LFcut)
      if (candidate.cosPACasc() < cosPACascMin) {
        resultSelections = false;
        registry.fill(HIST("hSelCosPACasc"), 0);
      } else {
        registry.fill(HIST("hSelCosPACasc"), 1);
      }
      if (candidate.cosPAV0() < cosPAV0Min) {
        resultSelections = false;
        registry.fill(HIST("hSelCosPAV0"), 0);
      } else {
        registry.fill(HIST("hSelCosPAV0"), 1);
      }

      // cascade and v0 daughters dca cut (LF cut)
      if (candidate.dcaCascDau() > dcaCascDauMax) {
        resultSelections = false;
        registry.fill(HIST("hSelDCACascDau"), 0);
      } else {
        registry.fill(HIST("hSelDCACascDau"), 1);
      }

      if (candidate.dcaV0Dau() > dcaV0DauMax) {
        resultSelections = false;
        registry.fill(HIST("hSelDCAV0Dau"), 0);
      } else {
        registry.fill(HIST("hSelDCAV0Dau"), 1);
      }

      // dca omegac daughters cut
      if (candidate.dcaOmegacDau() > dcaOmegacDauMax) {
        resultSelections = false;
        registry.fill(HIST("hSelDCAOmeDau"), 0);
      } else {
        registry.fill(HIST("hSelDCAOmeDau"), 1);
      }

      // cut on primary pion dcaXY and dcaZ
      if ((std::abs(candidate.impactParPrimaryPiXY()) < impactParameterXYPriPiMin) || (std::abs(candidate.impactParPrimaryPiXY()) > impactParameterXYPriPiMax)) {
        resultSelections = false;
        registry.fill(HIST("hSelDCAXYPrimPi"), 0);
      } else {
        registry.fill(HIST("hSelDCAXYPrimPi"), 1);
      }
      if ((std::abs(candidate.impactParPrimaryPiZ()) < impactParameterZPriPiMin) || (std::abs(candidate.impactParPrimaryPiZ()) > impactParameterZPriPiMax)) {
        resultSelections = false;
        registry.fill(HIST("hSelDCAZPrimPi"), 0);
      } else {
        registry.fill(HIST("hSelDCAZPrimPi"), 1);
      }

      // cut on cascade dcaXY and dcaZ
      if ((std::abs(candidate.impactParCascXY()) < impactParameterXYCascMin) || (std::abs(candidate.impactParCascXY()) > impactParameterXYCascMax)) {
        resultSelections = false;
        registry.fill(HIST("hSelDCAXYCasc"), 0);
      } else {
        registry.fill(HIST("hSelDCAXYCasc"), 1);
      }
      if ((std::abs(candidate.impactParCascZ()) < impactParameterZCascMin) || (std::abs(candidate.impactParCascZ()) > impactParameterZCascMax)) {
        resultSelections = false;
        registry.fill(HIST("hSelDCAZCasc"), 0);
      } else {
        registry.fill(HIST("hSelDCAZCasc"), 1);
      }

      // pT selections
      double ptPiFromCasc = RecoDecay::sqrtSumOfSquares(candidate.pxPiFromCasc(), candidate.pyPiFromCasc());
      double ptPiFromOme = RecoDecay::sqrtSumOfSquares(candidate.pxPrimaryPi(), candidate.pyPrimaryPi());
      if (std::abs(ptPiFromCasc) < ptPiFromCascMin) {
        resultSelections = false;
        registry.fill(HIST("hSelPtPiFromCasc"), 0);
      } else {
        registry.fill(HIST("hSelPtPiFromCasc"), 1);
      }
      if (std::abs(ptPiFromOme) < ptPiFromOmeMin) {
        resultSelections = false;
        registry.fill(HIST("hSelPtPiFromOme"), 0);
      } else {
        registry.fill(HIST("hSelPtPiFromOme"), 1);
      }

      //  TPC clusters selections
      if (trackPiFromOmeg.tpcNClsFound() < nClustersTpcMin) {
        resultSelections = false;
        registry.fill(HIST("hSelNClsTPCPiFromOme"), 0);
      } else {
        registry.fill(HIST("hSelNClsTPCPiFromOme"), 1);
      }
      if (trackPiFromLam.tpcNClsFound() < nClustersTpcMin) {
        resultSelections = false;
        registry.fill(HIST("hSelNClsTPCPiFromLam"), 0);
      } else {
        registry.fill(HIST("hSelNClsTPCPiFromLam"), 1);
      }
      if (trackPrFromLam.tpcNClsFound() < nClustersTpcMin) {
        resultSelections = false;
        registry.fill(HIST("hSelNClsTPCPrFromLam"), 0);
      } else {
        registry.fill(HIST("hSelNClsTPCPrFromLam"), 1);
      }
      if (trackPiFromCasc.tpcNClsFound() < nClustersTpcMin) {
        resultSelections = false;
        registry.fill(HIST("hSelNClsTPCPiFromCasc"), 0);
      } else {
        registry.fill(HIST("hSelNClsTPCPiFromCasc"), 1);
      }

      // TPC crossed rows selection
      if (trackPiFromOmeg.tpcNClsCrossedRows() < nTpcCrossedRowsMin) {
        resultSelections = false;
        registry.fill(HIST("hSelNCrossRowsTPCPiFromOme"), 0);
      } else {
        registry.fill(HIST("hSelNCrossRowsTPCPiFromOme"), 1);
      }
      if (trackPiFromLam.tpcNClsCrossedRows() < nTpcCrossedRowsMin) {
        resultSelections = false;
        registry.fill(HIST("hSelNCrossRowsTPCPiFromLam"), 0);
      } else {
        registry.fill(HIST("hSelNCrossRowsTPCPiFromLam"), 1);
      }
      if (trackPrFromLam.tpcNClsCrossedRows() < nTpcCrossedRowsMin) {
        resultSelections = false;
        registry.fill(HIST("hSelNCrossRowsTPCPrFromLam"), 0);
      } else {
        registry.fill(HIST("hSelNCrossRowsTPCPrFromLam"), 1);
      }
      if (trackPiFromCasc.tpcNClsCrossedRows() < nTpcCrossedRowsMin) {
        resultSelections = false;
        registry.fill(HIST("hSelNCrossRowsTPCPiFromCasc"), 0);
      } else {
        registry.fill(HIST("hSelNCrossRowsTPCPiFromCasc"), 1);
      }

      // further TPC selection
      if (trackPiFromOmeg.tpcCrossedRowsOverFindableCls() < tpcCrossedRowsOverFindableClustersRatioMin) {
        resultSelections = false;
        registry.fill(HIST("hSelCrossRowsOverFindClsTPCAllTracks"), 0);
      } else {
        registry.fill(HIST("hSelCrossRowsOverFindClsTPCAllTracks"), 1);
      }
      if (trackPiFromCasc.tpcCrossedRowsOverFindableCls() < tpcCrossedRowsOverFindableClustersRatioMin) {
        resultSelections = false;
        registry.fill(HIST("hSelCrossRowsOverFindClsTPCAllTracks"), 2);
      } else {
        registry.fill(HIST("hSelCrossRowsOverFindClsTPCAllTracks"), 3);
      }
      if (trackPiFromLam.tpcCrossedRowsOverFindableCls() < tpcCrossedRowsOverFindableClustersRatioMin) {
        resultSelections = false;
        registry.fill(HIST("hSelCrossRowsOverFindClsTPCAllTracks"), 4);
      } else {
        registry.fill(HIST("hSelCrossRowsOverFindClsTPCAllTracks"), 5);
      }
      if (trackPrFromLam.tpcCrossedRowsOverFindableCls() < tpcCrossedRowsOverFindableClustersRatioMin) {
        resultSelections = false;
        registry.fill(HIST("hSelCrossRowsOverFindClsTPCAllTracks"), 6);
      } else {
        registry.fill(HIST("hSelCrossRowsOverFindClsTPCAllTracks"), 7);
      }

      //  ITS clusters selection
      if (trackPiFromOmeg.itsNCls() < nClustersItsMin) {
        resultSelections = false;
        registry.fill(HIST("hSelNClsITSPiFromOme"), 0);
      } else {
        registry.fill(HIST("hSelNClsITSPiFromOme"), 1);
      }
      if (trackPiFromOmeg.itsNClsInnerBarrel() < nClustersItsInnBarrMin) {
        resultSelections = false;
        registry.fill(HIST("hSelNClsITSInnerPiFromOme"), 0);
      } else {
        registry.fill(HIST("hSelNClsITSInnerPiFromOme"), 1);
      }

      // track-level PID selection

      // for TrackSelectorPID
      int statusPidPrFromLam = -999;
      int statusPidPiFromLam = -999;
      int statusPidPiFromCasc = -999;
      int statusPidPiFromCharm = -999;

      bool statusPidLambda = false;
      bool statusPidCascade = false;
      bool statusPidCharm = false;

      int infoTpcStored = 0;
      int infoTofStored = 0;

      if (usePidTpcOnly == usePidTpcTofCombined) {
        LOGF(fatal, "Check the PID configurables, usePidTpcOnly and usePidTpcTofCombined can't have the same value");
      }

      if (trackPiFromLam.hasTPC()) {
        SETBIT(infoTpcStored, kPiFromLam);
      }
      if (trackPrFromLam.hasTPC()) {
        SETBIT(infoTpcStored, kPrFromLam);
      }
      if (trackPiFromCasc.hasTPC()) {
        SETBIT(infoTpcStored, kPiFromCasc);
      }
      if (trackPiFromOmeg.hasTPC()) {
        SETBIT(infoTpcStored, kPiFromCharm);
      }
      if (trackPiFromLam.hasTOF()) {
        SETBIT(infoTofStored, kPiFromLam);
      }
      if (trackPrFromLam.hasTOF()) {
        SETBIT(infoTofStored, kPrFromLam);
      }
      if (trackPiFromCasc.hasTOF()) {
        SETBIT(infoTofStored, kPiFromCasc);
      }
      if (trackPiFromOmeg.hasTOF()) {
        SETBIT(infoTofStored, kPiFromCharm);
      }

      if (usePidTpcOnly) {
        statusPidPrFromLam = selectorProton.statusTpc(trackPrFromLam);
        statusPidPiFromLam = selectorPion.statusTpc(trackPiFromLam);
        statusPidPiFromCasc = selectorPion.statusTpc(trackPiFromCasc);
        statusPidPiFromCharm = selectorPion.statusTpc(trackPiFromOmeg);
      } else if (usePidTpcTofCombined) {
        statusPidPrFromLam = selectorProton.statusTpcAndTof(trackPrFromLam);
        statusPidPiFromLam = selectorPion.statusTpcAndTof(trackPiFromLam);
        statusPidPiFromCasc = selectorPion.statusTpcAndTof(trackPiFromCasc);
        statusPidPiFromCharm = selectorPion.statusTpcAndTof(trackPiFromOmeg);
      }

      if (statusPidPrFromLam == TrackSelectorPID::Accepted && statusPidPiFromLam == TrackSelectorPID::Accepted) {
        statusPidLambda = true;
        if (resultSelections) {
          registry.fill(HIST("hTest"), 0.5);
        }
      }

      if (statusPidPrFromLam == TrackSelectorPID::Accepted && statusPidPiFromLam == TrackSelectorPID::Accepted && statusPidPiFromCasc == TrackSelectorPID::Accepted) {
        statusPidCascade = true;
        if (resultSelections) {
          registry.fill(HIST("hTest"), 1.5);
        }
      }

      if (statusPidPrFromLam == TrackSelectorPID::Accepted && statusPidPiFromLam == TrackSelectorPID::Accepted && statusPidPiFromCasc == TrackSelectorPID::Accepted && statusPidPiFromCharm == TrackSelectorPID::Accepted) {
        statusPidCharm = true;
        if (resultSelections) {
          registry.fill(HIST("hTest"), 2.5);
        }
      }

      // invariant mass cuts
      bool statusInvMassLambda = false;
      bool statusInvMassCascade = false;
      bool statusInvMassOmegac = false;

      double invMassLambda = candidate.invMassLambda();
      double invMassCascade = candidate.invMassCascade();
      double invMassOmegac = candidate.invMassOmegac();

      if (std::abs(invMassLambda - massLambdaFromPDG) < (nSigmaInvMassCut * sigmaInvMassLambda)) {
        statusInvMassLambda = true;
        registry.fill(HIST("hSelMassLam"), 1);
        if (statusPidLambda && statusPidCascade && statusPidCharm && resultSelections) {
          registry.fill(HIST("hTest"), 3.5);
        }
      } else {
        registry.fill(HIST("hSelMassLam"), 0);
      }

      if (std::abs(invMassCascade - massXiFromPDG) < (nSigmaInvMassCut * sigmaInvMassCascade)) {
        statusInvMassCascade = true;
        registry.fill(HIST("hSelMassCasc"), 1);
        if (statusPidLambda && statusPidCascade && statusPidCharm && statusInvMassLambda && resultSelections) {
          registry.fill(HIST("hTest"), 4.5);
        }
      } else {
        registry.fill(HIST("hSelMassCasc"), 0);
      }

      if ((invMassOmegac >= invMassOmegacMin) && (invMassOmegac <= invMassOmegacMax)) {
        statusInvMassOmegac = true;
        registry.fill(HIST("hSelMassOme"), 1);
        if (statusPidLambda && statusPidCascade && statusPidCharm && statusInvMassLambda && statusInvMassCascade && resultSelections) {
          registry.fill(HIST("hTest"), 5.5);
        }
      } else {
        registry.fill(HIST("hSelMassOme"), 0);
      }

      hfSelToXiPi(statusPidLambda, statusPidCascade, statusPidCharm, statusInvMassLambda, statusInvMassCascade, statusInvMassOmegac, resultSelections, infoTpcStored, infoTofStored,
                  trackPiFromOmeg.tpcNSigmaPi(), trackPiFromCasc.tpcNSigmaPi(), trackPiFromLam.tpcNSigmaPi(), trackPrFromLam.tpcNSigmaPr(),
                  trackPiFromOmeg.tofNSigmaPi(), trackPiFromCasc.tofNSigmaPi(), trackPiFromLam.tofNSigmaPi(), trackPrFromLam.tofNSigmaPr());

      if (resultSelections) {
        if (!statusPidLambda) {
          registry.fill(HIST("hSelPID"), 0.5);
        }
        if (statusPidLambda) {
          registry.fill(HIST("hSelPID"), 1.5);
        }
        if (!statusPidCascade) {
          registry.fill(HIST("hSelPID"), 2.5);
        }
        if (statusPidCascade) {
          registry.fill(HIST("hSelPID"), 3.5);
        }
        if (!statusPidCharm) {
          registry.fill(HIST("hSelPID"), 4.5);
        }
        if (statusPidCharm) {
          registry.fill(HIST("hSelPID"), 5.5);
        }
        if (!statusInvMassLambda) {
          registry.fill(HIST("hSelPID"), 6.5);
        }
        if (statusInvMassLambda) {
          registry.fill(HIST("hSelPID"), 7.5);
        }
        if (!statusInvMassCascade) {
          registry.fill(HIST("hSelPID"), 8.5);
        }
        if (statusInvMassCascade) {
          registry.fill(HIST("hSelPID"), 9.5);
        }
        if (!statusInvMassOmegac) {
          registry.fill(HIST("hSelPID"), 10.5);
        }
        if (statusInvMassOmegac) {
          registry.fill(HIST("hSelPID"), 11.5);
        }
      }

      if (statusPidLambda && statusPidCascade && statusPidCharm && statusInvMassLambda && statusInvMassCascade && statusInvMassOmegac && resultSelections) {
        hInvMassOmegac->Fill(invMassOmegac);

        if (candidate.collisionId() != collId) {
          hNEventsSaved->Fill(0.5);
          collId = candidate.collisionId();
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateSelectorToXiPi>(cfgc)};
}
