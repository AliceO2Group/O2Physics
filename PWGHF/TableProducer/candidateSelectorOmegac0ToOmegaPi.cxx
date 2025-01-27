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

/// \file candidateSelectorToOmegaPi.cxx
/// \brief Omegac0 → Omega Pi selection task
/// \author Federica Zanone <federica.zanone@cern.ch>, Heidelberg University
/// \author Ruiqi Yin <ruiqi.yin@cern.ch>, Fudan University

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectorPID.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsAnalysis.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::analysis;

enum pidInfoStored {
  kPiFromLam = 0,
  kPrFromLam,
  kKaFromCasc,
  kPiFromCharm
};

/// Struct for applying Omegac0 -> Omega pi selection cuts
struct HfCandidateSelectorToOmegaPi {
  Produces<aod::HfSelToOmegaPi> hfSelToOmegaPi;

  // LF analysis selections
  Configurable<double> radiusCascMin{"radiusCascMin", 0.5, "Min cascade radius"};
  Configurable<double> radiusV0Min{"radiusV0Min", 1.1, "Min V0 radius"};
  Configurable<double> cosPAV0Min{"cosPAV0Min", 0.97, "Min valueCosPA V0"};
  Configurable<double> cosPACascMin{"cosPACascMin", 0.97, "Min value CosPA cascade"};
  Configurable<double> dcaCascDauMax{"dcaCascDauMax", 1.0, "Max DCA cascade daughters"};
  Configurable<double> dcaV0DauMax{"dcaV0DauMax", 1.0, "Max DCA V0 daughters"};
  Configurable<float> dcaBachToPvMin{"dcaBachToPvMin", 0.04, "DCA Bach To PV"};
  Configurable<float> dcaNegToPvMin{"dcaNegToPvMin", 0.06, "DCA Neg To PV"};
  Configurable<float> dcaPosToPvMin{"dcaPosToPvMin", 0.06, "DCA Pos To PV"};
  Configurable<float> v0MassWindow{"v0MassWindow", 0.01, "V0 mass window"};
  Configurable<float> cascadeMassWindow{"cascadeMassWindow", 0.01, "Cascade mass window"};
  Configurable<bool> applyTrkSelLf{"applyTrkSelLf", true, "Apply track selection for LF daughters"};

  // limit charm baryon invariant mass spectrum
  Configurable<double> invMassCharmBaryonMin{"invMassCharmBaryonMin", 2.3, "Lower limit invariant mass spectrum charm baryon"}; // 2.4 Omegac0 only
  Configurable<double> invMassCharmBaryonMax{"invMassCharmBaryonMax", 3.1, "Upper limit invariant mass spectrum charm baryon"};

  // kinematic selections
  Configurable<double> etaTrackCharmBachMax{"etaTrackCharmBachMax", 0.8, "Max absolute value of eta for charm baryon bachelor"};
  Configurable<double> etaTrackLFDauMax{"etaTrackLFDauMax", 1.0, "Max absolute value of eta for V0 and cascade daughters"};
  Configurable<double> ptKaFromCascMin{"ptKaFromCascMin", 0.15, "Min pT kaon <- casc"};
  Configurable<double> ptPiFromCharmBaryonMin{"ptPiFromCharmBaryonMin", 0.2, "Min pT pi <- charm baryon"};

  Configurable<double> impactParameterXYPiFromCharmBaryonMin{"impactParameterXYPiFromCharmBaryonMin", 0., "Min dcaxy pi from charm baryon track to PV"};
  Configurable<double> impactParameterXYPiFromCharmBaryonMax{"impactParameterXYPiFromCharmBaryonMax", 10., "Max dcaxy pi from charm baryon track to PV"};
  Configurable<double> impactParameterZPiFromCharmBaryonMin{"impactParameterZPiFromCharmBaryonMin", 0., "Min dcaz pi from charm baryon track to PV"};
  Configurable<double> impactParameterZPiFromCharmBaryonMax{"impactParameterZPiFromCharmBaryonMax", 10., "Max dcaz pi from charm baryon track to PV"};

  Configurable<double> impactParameterXYCascMin{"impactParameterXYCascMin", 0., "Min dcaxy cascade track to PV"};
  Configurable<double> impactParameterXYCascMax{"impactParameterXYCascMax", 10., "Max dcaxy cascade track to PV"};
  Configurable<double> impactParameterZCascMin{"impactParameterZCascMin", 0., "Min dcaz cascade track to PV"};
  Configurable<double> impactParameterZCascMax{"impactParameterZCascMax", 10., "Max dcaz cascade track to PV"};

  Configurable<double> ptCandMin{"ptCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> ptCandMax{"ptCandMax", 50., "Upper bound of candidate pT"};

  Configurable<double> dcaCharmBaryonDauMax{"dcaCharmBaryonDauMax", 2.0, "Max DCA charm baryon daughters"};

  // PID options
  Configurable<bool> usePidTpcOnly{"usePidTpcOnly", false, "Perform PID using only TPC"};
  Configurable<bool> usePidTpcTofCombined{"usePidTpcTofCombined", true, "Perform PID using TPC & TOF"};

  // PID - TPC selections
  Configurable<double> ptPiPidTpcMin{"ptPiPidTpcMin", -1, "Lower bound of track pT for TPC PID for pion selection"};
  Configurable<double> ptPiPidTpcMax{"ptPiPidTpcMax", 9999.9, "Upper bound of track pT for TPC PID for pion selection"};
  Configurable<double> nSigmaTpcPiMax{"nSigmaTpcPiMax", 3., "Nsigma cut on TPC only for pion selection"};
  Configurable<double> nSigmaTpcCombinedPiMax{"nSigmaTpcCombinedPiMax", 0., "Nsigma cut on TPC combined with TOF for pion selection"};

  Configurable<double> ptPrPidTpcMin{"ptPrPidTpcMin", -1, "Lower bound of track pT for TPC PID for proton selection"};
  Configurable<double> ptPrPidTpcMax{"ptPrPidTpcMax", 9999.9, "Upper bound of track pT for TPC PID for proton selection"};
  Configurable<double> nSigmaTpcPrMax{"nSigmaTpcPrMax", 3., "Nsigma cut on TPC only for proton selection"};
  Configurable<double> nSigmaTpcCombinedPrMax{"nSigmaTpcCombinedPrMax", 0., "Nsigma cut on TPC combined with TOF for proton selection"};

  Configurable<double> ptKaPidTpcMin{"ptKaPidTpcMin", -1, "Lower bound of track pT for TPC PID for kaon selection"};
  Configurable<double> ptKaPidTpcMax{"ptKaPidTpcMax", 9999.9, "Upper bound of track pT for TPC PID for kaon selection"};
  Configurable<double> nSigmaTpcKaMax{"nSigmaTpcKaMax", 3., "Nsigma cut on TPC only for kaon selection"};
  Configurable<double> nSigmaTpcCombinedKaMax{"nSigmaTpcCombinedKaMax", 0., "Nsigma cut on TPC combined with TOF for kaon selection"};

  // PID - TOF selections
  Configurable<double> ptPiPidTofMin{"ptPiPidTofMin", -1, "Lower bound of track pT for TOF PID for pion selection"};
  Configurable<double> ptPiPidTofMax{"ptPiPidTofMax", 9999.9, "Upper bound of track pT for TOF PID for pion selection"};
  Configurable<double> nSigmaTofPiMax{"nSigmaTofPiMax", 3., "Nsigma cut on TOF only for pion selection"};
  Configurable<double> nSigmaTofCombinedPiMax{"nSigmaTofCombinedPiMax", 0., "Nsigma cut on TOF combined with TPC for pion selection"};

  Configurable<double> ptPrPidTofMin{"ptPrPidTofMin", -1, "Lower bound of track pT for TOF PID for proton selection"};
  Configurable<double> ptPrPidTofMax{"ptPrPidTofMax", 9999.9, "Upper bound of track pT for TOF PID for proton selection"};
  Configurable<double> nSigmaTofPrMax{"nSigmaTofPrMax", 3., "Nsigma cut on TOF only for proton selection"};
  Configurable<double> nSigmaTofCombinedPrMax{"nSigmaTofCombinedPrMax", 0., "Nsigma cut on TOF combined with TPC for proton selection"};

  Configurable<double> ptKaPidTofMin{"ptKaPidTofMin", -1, "Lower bound of track pT for TOF PID for kaon selection"};
  Configurable<double> ptKaPidTofMax{"ptKaPidTofMax", 9999.9, "Upper bound of track pT for TOF PID for kaon selection"};
  Configurable<double> nSigmaTofKaMax{"nSigmaTofKaMax", 3., "Nsigma cut on TOF only for kaon selection"};
  Configurable<double> nSigmaTofCombinedKaMax{"nSigmaTofCombinedKaMax", 0., "Nsigma cut on TOF combined with TOF for kaon selection"};

  // detector clusters selections
  Configurable<int> nClustersTpcMin{"nClustersTpcMin", 70, "Minimum number of TPC clusters requirement"};
  Configurable<int> nTpcCrossedRowsMin{"nTpcCrossedRowsMin", 70, "Minimum number of TPC crossed rows requirement"};
  Configurable<double> tpcCrossedRowsOverFindableClustersRatioMin{"tpcCrossedRowsOverFindableClustersRatioMin", 0.8, "Minimum ratio TPC crossed rows over findable clusters requirement"};
  Configurable<float> tpcChi2PerClusterMax{"tpcChi2PerClusterMax", 4, "Maximum value of chi2 fit over TPC clusters"};
  Configurable<int> nClustersItsMin{"nClustersItsMin", 3, "Minimum number of ITS clusters requirement for pi <- charm baryon"};
  Configurable<int> nClustersItsInnBarrMin{"nClustersItsInnBarrMin", 1, "Minimum number of ITS clusters in inner barrel requirement for pi <- charm baryon"};
  Configurable<float> itsChi2PerClusterMax{"itsChi2PerClusterMax", 36, "Maximum value of chi2 fit over ITS clusters for pi <- charm baryon"};

  // KF selection
  Configurable<bool> applyKFpreselections{"applyKFpreselections", false, "Apply KFParticle related rejection"};
  Configurable<bool> applyCompetingCascRejection{"applyCompetingCascRejection", false, "Apply competing Xi(for Omegac0) rejection"};
  Configurable<float> cascadeRejMassWindow{"cascadeRejMassWindow", 0.01, "competing Xi(for Omegac0) rejection mass window"};
  Configurable<float> v0LdlMin{"v0LdlMin", 3., "Minimum value of l/dl of V0"}; // l/dl and Chi2 are to be determined
  Configurable<float> cascLdlMin{"cascLdlMin", 1., "Minimum value of l/dl of casc"};
  Configurable<float> omegacLdlMax{"omegacLdlMax", 5., "Maximum value of l/dl of Omegac"};
  Configurable<float> cTauOmegacMax{"cTauOmegacMax", 0.4, "lifetime τ of Omegac"};
  Configurable<float> v0Chi2OverNdfMax{"v0Chi2OverNdfMax", 100., "Maximum chi2Geo/NDF of V0"};
  Configurable<float> cascChi2OverNdfMax{"cascChi2OverNdfMax", 100., "Maximum chi2Geo/NDF of casc"};
  Configurable<float> omegacChi2OverNdfMax{"omegacChi2OverNdfMax", 100., "Maximum chi2Geo/NDF of Omegac"};
  Configurable<float> chi2TopoV0ToCascMax{"chi2TopoV0ToCascMax", 100., "Maximum chi2Topo/NDF of V0ToCas"};
  Configurable<float> chi2TopoOmegacToPvMax{"chi2TopoOmegacToPvMax", 100., "Maximum chi2Topo/NDF of OmegacToPv"};
  Configurable<float> chi2TopoCascToOmegacMax{"chi2TopoCascToOmegacMax", 100., "Maximum chi2Topo/NDF of CascToOmegac"};
  Configurable<float> chi2TopoCascToPvMax{"chi2TopoCascToPvMax", 100., "Maximum chi2Topo/NDF of CascToPv"};
  Configurable<float> decayLenXYOmegacMax{"decayLenXYOmegacMax", 1.5, "Maximum decay lengthXY of Omegac"};
  Configurable<float> decayLenXYCascMin{"decayLenXYCascMin", 1., "Minimum decay lengthXY of Cascade"};
  Configurable<float> decayLenXYLambdaMin{"decayLenXYLambdaMin", 0., "Minimum decay lengthXY of V0"};
  Configurable<float> cosPaCascToOmegacMin{"cosPaCascToOmegacMin", 0.995, "Minimum cosPA of cascade<-Omegac"};
  Configurable<float> cosPaV0ToCascMin{"cosPaV0ToCascMin", 0.99, "Minimum cosPA of V0<-cascade"};

  TrackSelectorPi selectorPion;
  TrackSelectorPr selectorProton;
  TrackSelectorKa selectorKaon;

  using TracksSel = soa::Join<aod::TracksWDcaExtra, aod::TracksPidPi, aod::TracksPidPr, aod::TracksPidKa>;
  using TracksSelLf = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksPidPi, aod::TracksPidPr, aod::TracksPidKa>;

  HistogramRegistry registry{"registry"}; // for QA of selections

  OutputObj<TH1D> hInvMassCharmBaryon{TH1D("hInvMassCharmBaryon", "Charm baryon invariant mass;inv mass;entries", 500, 2.3, 3.1)};

  void init(InitContext const&)
  {
    std::array<bool, 2> processesSelector = {doprocessOmegac0SelectorWithDCAFitter, doprocessOmegac0SelectorWithKFParticle};
    const int nProcessesSelector = std::accumulate(processesSelector.begin(), processesSelector.end(), 0);
    if (nProcessesSelector != 1) {
      LOGP(fatal, "At most one process function for selector can be enabled at a time.");
    }

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

    selectorKaon.setRangePtTpc(ptKaPidTpcMin, ptKaPidTpcMax);
    selectorKaon.setRangeNSigmaTpc(-nSigmaTpcKaMax, nSigmaTpcKaMax);
    selectorKaon.setRangeNSigmaTpcCondTof(-nSigmaTpcCombinedKaMax, nSigmaTpcCombinedKaMax);
    selectorKaon.setRangePtTof(ptKaPidTofMin, ptKaPidTofMax);
    selectorKaon.setRangeNSigmaTof(-nSigmaTofKaMax, nSigmaTofKaMax);
    selectorKaon.setRangeNSigmaTofCondTpc(-nSigmaTofCombinedKaMax, nSigmaTofCombinedKaMax);

    const AxisSpec axisSel{2, -0.5, 1.5, "status"};

    registry.add("hSelPID", "hSelPID;status;entries", {HistType::kTH1D, {{12, 0., 12.}}});
    registry.add("hStatusCheck", "Check consecutive selections status;status;entries", {HistType::kTH1D, {{12, 0., 12.}}});

    // for QA of the selections (bin 0 -> candidates that did not pass the selection, bin 1 -> candidates that passed the selection)
    registry.add("hSelSignDec", "hSelSignDec;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelEtaPosV0Dau", "hSelEtaPosV0Dau;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelEtaNegV0Dau", "hSelEtaNegV0Dau;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelEtaKaFromCasc", "hSelEtaKaFromCasc;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelEtaPiFromCharm", "hSelEtaPiFromCharm;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelRadCasc", "hSelRadCasc;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelRadV0", "hSelRadV0;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelCosPACasc", "hSelCosPACasc;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelCosPAV0", "hSelCosPAV0;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelDCACascDau", "hSelDCACascDau;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelDCAV0Dau", "hSelDCAV0Dau;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelDCACharmDau", "hSelDCACharmDau;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelDCAXYPrimPi", "hSelDCAXYPrimPi;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelDCAZPrimPi", "hSelDCAZPrimPi;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelDCAXYCasc", "hSelDCAXYCasc;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelDCAZCasc", "hSelDCAZCasc;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelPtKaFromCasc", "hSelPtKaFromCasc;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelPtPiFromCharm", "hSelPtPiFromCharm;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelTPCQualityPiFromCharm", "hSelTPCQualityPiFromCharm;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelTPCQualityPiFromLam", "hSelTPCQualityPiFromLam;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelTPCQualityPrFromLam", "hSelTPCQualityPrFromLam;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelTPCQualityKaFromCasc", "hSelTPCQualityKaFromCasc;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelITSQualityPiFromCharm", "hSelITSQualityPiFromCharm;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelMassLam", "hSelMassLam;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelMassCasc", "hSelMassCasc;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelMassCharmBaryon", "hSelMassCharmBaryon;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelDcaXYToPvV0Daughters", "hSelDcaXYToPvV0Daughters;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelDcaXYToPvKaFromCasc", "hSelDcaXYToPvKaFromCasc;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelCompetingCasc", "hSelCompetingCasc;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelKFstatus", "hSelKFstatus;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelV0_Casc_Omegacldl", "hSelV0_Casc_Omegacldl;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelctauOmegac", "hSelctauOmegac;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelChi2GeooverNDFV0_Casc_Omegac", "hSelChi2GeooverNDFV0_Casc_Omegac;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelChi2TopooverNDFV0_Casc_Omegac", "hSelChi2TopooverNDFV0_Casc_Omegac;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSeldecayLenXYOmegac_Casc_V0", "hSeldecayLenXYOmegac_Casc_V0;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelcosPaCascToOmegac_V0ToCasc", "hSelcosPaCascToOmegac_V0ToCasc;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hInvMassXiMinus_rej_cut", "hInvMassXiMinus_rej_cut", kTH1D, {{1000, 1.25f, 1.65f}});
  }

  template <int ConstructMethod, typename Candidates>
  void runOmegac0Selector(const Candidates& candidates,
                          TracksSel const& tracks,
                          TracksSelLf const& lfTracks)
  {
    // looping over charm baryon candidates
    for (const auto& candidate : candidates) {

      bool resultSelections = true; // True if the candidate passes all the selections, False otherwise

      auto trackV0PosDauId = candidate.posTrackId();                   // positive V0 daughter
      auto trackV0NegDauId = candidate.negTrackId();                   // negative V0 daughter
      auto trackKaFromCascId = candidate.bachelorId();                 // kaon <- cascade
      auto trackPiFromCharmId = candidate.bachelorFromCharmBaryonId(); // pion <- charm baryon
      auto trackV0PosDau = lfTracks.rawIteratorAt(trackV0PosDauId);
      auto trackV0NegDau = lfTracks.rawIteratorAt(trackV0NegDauId);
      auto trackKaFromCasc = lfTracks.rawIteratorAt(trackKaFromCascId);
      auto trackPiFromCharm = tracks.rawIteratorAt(trackPiFromCharmId);

      auto trackPiFromLam = trackV0NegDau;
      auto trackPrFromLam = trackV0PosDau;

      int8_t signDecay = candidate.signDecay(); // sign of pi <- cascade

      if (signDecay > 0) {
        trackPiFromLam = trackV0PosDau;
        trackPrFromLam = trackV0NegDau;
        registry.fill(HIST("hSelSignDec"), 1); // anti-particle decay
      } else if (signDecay < 0) {
        registry.fill(HIST("hSelSignDec"), 0); // particle decay
      }

      // eta selection
      double etaV0PosDau = candidate.etaV0PosDau();
      double etaV0NegDau = candidate.etaV0NegDau();
      double etaKaFromCasc = candidate.etaBachFromCasc();
      double etaPiFromCharmBaryon = candidate.etaBachFromCharmBaryon();
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
      if (std::abs(etaKaFromCasc) > etaTrackLFDauMax) {
        resultSelections = false;
        registry.fill(HIST("hSelEtaKaFromCasc"), 0);
      } else {
        registry.fill(HIST("hSelEtaKaFromCasc"), 1);
      }
      if (std::abs(etaPiFromCharmBaryon) > etaTrackCharmBachMax) {
        resultSelections = false;
        registry.fill(HIST("hSelEtaPiFromCharm"), 0);
      } else {
        registry.fill(HIST("hSelEtaPiFromCharm"), 1);
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

      // dca charm baryon daughters cut
      if (candidate.dcaCharmBaryonDau() > dcaCharmBaryonDauMax) {
        resultSelections = false;
        registry.fill(HIST("hSelDCACharmDau"), 0);
      } else {
        registry.fill(HIST("hSelDCACharmDau"), 1);
      }

      // dcaXY v0 daughters to PV cut
      if (std::abs(candidate.dcaXYToPvV0Dau0()) < dcaPosToPvMin || std::abs(candidate.dcaXYToPvV0Dau1()) < dcaNegToPvMin) {
        resultSelections = false;
        registry.fill(HIST("hSelDcaXYToPvV0Daughters"), 0);
      } else {
        registry.fill(HIST("hSelDcaXYToPvV0Daughters"), 1);
      }

      // dcaXY ka <-- cascade to PV cut
      if (std::abs(candidate.dcaXYToPvCascDau()) < dcaBachToPvMin) {
        resultSelections = false;
        registry.fill(HIST("hSelDcaXYToPvKaFromCasc"), 0);
      } else {
        registry.fill(HIST("hSelDcaXYToPvKaFromCasc"), 1);
      }

      // cut on charm bachelor pion dcaXY and dcaZ
      if ((std::abs(candidate.impactParBachFromCharmBaryonXY()) < impactParameterXYPiFromCharmBaryonMin) || (std::abs(candidate.impactParBachFromCharmBaryonXY()) > impactParameterXYPiFromCharmBaryonMax)) {
        resultSelections = false;
        registry.fill(HIST("hSelDCAXYPrimPi"), 0);
      } else {
        registry.fill(HIST("hSelDCAXYPrimPi"), 1);
      }
      if ((std::abs(candidate.impactParBachFromCharmBaryonZ()) < impactParameterZPiFromCharmBaryonMin) || (std::abs(candidate.impactParBachFromCharmBaryonZ()) > impactParameterZPiFromCharmBaryonMax)) {
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
      double ptKaFromCasc = RecoDecay::sqrtSumOfSquares(candidate.pxBachFromCasc(), candidate.pyBachFromCasc());
      double ptPiFromCharmBaryon = RecoDecay::sqrtSumOfSquares(candidate.pxBachFromCharmBaryon(), candidate.pyBachFromCharmBaryon());
      if (std::abs(ptKaFromCasc) < ptKaFromCascMin) {
        resultSelections = false;
        registry.fill(HIST("hSelPtKaFromCasc"), 0);
      } else {
        registry.fill(HIST("hSelPtKaFromCasc"), 1);
      }
      if (std::abs(ptPiFromCharmBaryon) < ptPiFromCharmBaryonMin) {
        resultSelections = false;
        registry.fill(HIST("hSelPtPiFromCharm"), 0);
      } else {
        registry.fill(HIST("hSelPtPiFromCharm"), 1);
      }

      if constexpr (ConstructMethod == hf_cand_casc_lf::ConstructMethod::KfParticle) {
        // KFParticle Preselections(kfsel)
        if (applyKFpreselections) {

          bool inputKF = false;
          if (resultSelections) {
            inputKF = true;
            registry.fill(HIST("hSelKFstatus"), 0);
          }

          //  Competing Ξ rejection(KF)  Try to reject cases in which the candidate has a an inv. mass compatibler to Xi (bachelor pion) instead of Omega (bachelor kaon)
          if (applyCompetingCascRejection) {
            if (std::abs(candidate.cascRejectInvmass() - o2::constants::physics::MassXiMinus) < cascadeRejMassWindow) {
              resultSelections = false;
              registry.fill(HIST("hSelCompetingCasc"), 0);
            } else {
              registry.fill(HIST("hSelCompetingCasc"), 1);
              registry.fill(HIST("hInvMassXiMinus_rej_cut"), candidate.cascRejectInvmass());
            }
          }

          // v0&Casc&Omegac ldl selection
          if ((candidate.v0ldl() < v0LdlMin) || (candidate.cascldl() < cascLdlMin) || (candidate.omegacldl() > omegacLdlMax)) {
            resultSelections = false;
            registry.fill(HIST("hSelV0_Casc_Omegacldl"), 0);
          } else {
            registry.fill(HIST("hSelV0_Casc_Omegacldl"), 1);
          }

          // Omegac ctau selsection
          if (candidate.ctauOmegac() > cTauOmegacMax) {
            resultSelections = false;
            registry.fill(HIST("hSelctauOmegac"), 0);
          } else {
            registry.fill(HIST("hSelctauOmegac"), 1);
          }

          // Chi2Geo/NDF V0&Casc&Omegac selection
          if ((candidate.v0Chi2OverNdf() > v0Chi2OverNdfMax) || (candidate.cascChi2OverNdf() > cascChi2OverNdfMax) || (candidate.omegacChi2OverNdf() > omegacChi2OverNdfMax)) {
            resultSelections = false;
            registry.fill(HIST("hSelChi2GeooverNDFV0_Casc_Omegac"), 0);
          } else {
            registry.fill(HIST("hSelChi2GeooverNDFV0_Casc_Omegac"), 1);
          }

          // Chi2Topo/NDF (chi2TopoV0ToCasc chi2TopoOmegacToPv chi2TopoCascToOmegac chi2TopoCascToPv) selection  (???????????/NDF of which particle????????)
          if ((candidate.chi2TopoV0ToCasc() > chi2TopoV0ToCascMax) || (candidate.chi2TopoOmegacToPv() > chi2TopoOmegacToPvMax) || (candidate.chi2TopoCascToOmegac() > chi2TopoCascToOmegacMax) || (candidate.chi2TopoCascToPv() > chi2TopoCascToPvMax)) {
            resultSelections = false;
            registry.fill(HIST("hSelChi2TopooverNDFV0_Casc_Omegac"), 0);
          } else {
            registry.fill(HIST("hSelChi2TopooverNDFV0_Casc_Omegac"), 1);
          }

          // DecaylengthXY of Omegac&Casc&V0 selection
          if ((std::abs(candidate.decayLenXYOmegac()) > decayLenXYOmegacMax) || (std::abs(candidate.decayLenXYCasc()) < decayLenXYCascMin) || (std::abs(candidate.decayLenXYLambda()) < decayLenXYLambdaMin)) {
            resultSelections = false;
            registry.fill(HIST("hSeldecayLenXYOmegac_Casc_V0"), 0);
          } else {
            registry.fill(HIST("hSeldecayLenXYOmegac_Casc_V0"), 1);
          }

          // KFPA cut cosPaCascToOmegac cosPaV0ToCasc
          if ((candidate.cosPaCascToOmegac() < cosPaCascToOmegacMin) || (candidate.cosPaV0ToCasc() < cosPaV0ToCascMin)) {
            resultSelections = false;
            registry.fill(HIST("hSelcosPaCascToOmegac_V0ToCasc"), 0);
          } else {
            registry.fill(HIST("hSelcosPaCascToOmegac_V0ToCasc"), 1);
          }

          if (resultSelections && inputKF) {
            registry.fill(HIST("hSelKFstatus"), 1);
          }
        }
      }

      //  TPC clusters selections
      if (applyTrkSelLf) {
        if (!isSelectedTrackTpcQuality(trackPiFromLam, nClustersTpcMin, nTpcCrossedRowsMin, tpcCrossedRowsOverFindableClustersRatioMin, tpcChi2PerClusterMax)) {
          resultSelections = false;
          registry.fill(HIST("hSelTPCQualityPiFromLam"), 0);
        } else {
          registry.fill(HIST("hSelTPCQualityPiFromLam"), 1);
        }
        if (!isSelectedTrackTpcQuality(trackPrFromLam, nClustersTpcMin, nTpcCrossedRowsMin, tpcCrossedRowsOverFindableClustersRatioMin, tpcChi2PerClusterMax)) {
          resultSelections = false;
          registry.fill(HIST("hSelTPCQualityPrFromLam"), 0);
        } else {
          registry.fill(HIST("hSelTPCQualityPrFromLam"), 1);
        }
        if (!isSelectedTrackTpcQuality(trackKaFromCasc, nClustersTpcMin, nTpcCrossedRowsMin, tpcCrossedRowsOverFindableClustersRatioMin, tpcChi2PerClusterMax)) {
          resultSelections = false;
          registry.fill(HIST("hSelTPCQualityKaFromCasc"), 0);
        } else {
          registry.fill(HIST("hSelTPCQualityKaFromCasc"), 1);
        }
      }
      if (!isSelectedTrackTpcQuality(trackPiFromCharm, nClustersTpcMin, nTpcCrossedRowsMin, tpcCrossedRowsOverFindableClustersRatioMin, tpcChi2PerClusterMax)) {
        resultSelections = false;
        registry.fill(HIST("hSelTPCQualityPiFromCharm"), 0);
      } else {
        registry.fill(HIST("hSelTPCQualityPiFromCharm"), 1);
      }

      //  ITS clusters selection
      if (!isSelectedTrackItsQuality(trackPiFromCharm, nClustersItsMin, itsChi2PerClusterMax) || trackPiFromCharm.itsNClsInnerBarrel() < nClustersItsInnBarrMin) {
        resultSelections = false;
        registry.fill(HIST("hSelITSQualityPiFromCharm"), 0);
      } else {
        registry.fill(HIST("hSelITSQualityPiFromCharm"), 1);
      }

      // track-level PID selection

      // for TrackSelectorPID
      int statusPidPrFromLam = -999;
      int statusPidPiFromLam = -999;
      int statusPidKaFromCasc = -999;
      int statusPidPiFromCharmBaryon = -999;

      bool statusPidLambda = false;
      bool statusPidCascade = false;
      bool statusPidCharmBaryon = false;

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
      if (trackKaFromCasc.hasTPC()) {
        SETBIT(infoTpcStored, kKaFromCasc);
      }
      if (trackPiFromCharm.hasTPC()) {
        SETBIT(infoTpcStored, kPiFromCharm);
      }
      if (trackPiFromLam.hasTOF()) {
        SETBIT(infoTofStored, kPiFromLam);
      }
      if (trackPrFromLam.hasTOF()) {
        SETBIT(infoTofStored, kPrFromLam);
      }
      if (trackKaFromCasc.hasTOF()) {
        SETBIT(infoTofStored, kKaFromCasc);
      }
      if (trackPiFromCharm.hasTOF()) {
        SETBIT(infoTofStored, kPiFromCharm);
      }

      if (usePidTpcOnly) {
        statusPidPrFromLam = selectorProton.statusTpc(trackPrFromLam);
        statusPidPiFromLam = selectorPion.statusTpc(trackPiFromLam);
        statusPidKaFromCasc = selectorKaon.statusTpc(trackKaFromCasc);
        statusPidPiFromCharmBaryon = selectorPion.statusTpc(trackPiFromCharm);
      } else if (usePidTpcTofCombined) {
        statusPidPrFromLam = selectorProton.statusTpcOrTof(trackPrFromLam);
        statusPidPiFromLam = selectorPion.statusTpcOrTof(trackPiFromLam);
        statusPidKaFromCasc = selectorKaon.statusTpcOrTof(trackKaFromCasc);
        statusPidPiFromCharmBaryon = selectorPion.statusTpcOrTof(trackPiFromCharm);
      }

      if (statusPidPrFromLam == TrackSelectorPID::Accepted && statusPidPiFromLam == TrackSelectorPID::Accepted) {
        statusPidLambda = true;
        if (resultSelections) {
          registry.fill(HIST("hStatusCheck"), 0.5);
        }
      }

      if (statusPidPrFromLam == TrackSelectorPID::Accepted && statusPidPiFromLam == TrackSelectorPID::Accepted && statusPidKaFromCasc == TrackSelectorPID::Accepted) {
        statusPidCascade = true;
        if (resultSelections) {
          registry.fill(HIST("hStatusCheck"), 1.5);
        }
      }

      if (statusPidPrFromLam == TrackSelectorPID::Accepted && statusPidPiFromLam == TrackSelectorPID::Accepted && statusPidKaFromCasc == TrackSelectorPID::Accepted && statusPidPiFromCharmBaryon == TrackSelectorPID::Accepted) {
        statusPidCharmBaryon = true;
        if (resultSelections) {
          registry.fill(HIST("hStatusCheck"), 2.5);
        }
      }

      // invariant mass cuts
      bool statusInvMassLambda = false;
      bool statusInvMassCascade = false;
      bool statusInvMassCharmBaryon = false;

      double invMassLambda = candidate.invMassLambda();
      double invMassCascade = candidate.invMassCascade();
      double invMassCharmBaryon = candidate.invMassCharmBaryon();

      if (std::abs(invMassLambda - o2::constants::physics::MassLambda0) < v0MassWindow) {
        statusInvMassLambda = true;
        registry.fill(HIST("hSelMassLam"), 1);
        if (statusPidLambda && statusPidCascade && statusPidCharmBaryon && resultSelections) {
          registry.fill(HIST("hStatusCheck"), 3.5);
        }
      } else {
        registry.fill(HIST("hSelMassLam"), 0);
      }

      if (std::abs(invMassCascade - o2::constants::physics::MassOmegaMinus) < cascadeMassWindow) {
        statusInvMassCascade = true;
        registry.fill(HIST("hSelMassCasc"), 1);
        if (statusPidLambda && statusPidCascade && statusPidCharmBaryon && statusInvMassLambda && resultSelections) {
          registry.fill(HIST("hStatusCheck"), 4.5);
        }
      } else {
        registry.fill(HIST("hSelMassCasc"), 0);
      }

      if ((invMassCharmBaryon >= invMassCharmBaryonMin) && (invMassCharmBaryon <= invMassCharmBaryonMax)) {
        statusInvMassCharmBaryon = true;
        registry.fill(HIST("hSelMassCharmBaryon"), 1);
        if (statusPidLambda && statusPidCascade && statusPidCharmBaryon && statusInvMassLambda && statusInvMassCascade && resultSelections) {
          registry.fill(HIST("hStatusCheck"), 5.5);
        }
      } else {
        registry.fill(HIST("hSelMassCharmBaryon"), 0);
      }

      hfSelToOmegaPi(statusPidLambda, statusPidCascade, statusPidCharmBaryon, statusInvMassLambda, statusInvMassCascade, statusInvMassCharmBaryon, resultSelections, infoTpcStored, infoTofStored,
                     trackPiFromCharm.tpcNSigmaPi(), trackKaFromCasc.tpcNSigmaKa(), trackPiFromLam.tpcNSigmaPi(), trackPrFromLam.tpcNSigmaPr(),
                     trackPiFromCharm.tofNSigmaPi(), trackKaFromCasc.tofNSigmaKa(), trackPiFromLam.tofNSigmaPi(), trackPrFromLam.tofNSigmaPr());

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
        if (!statusPidCharmBaryon) {
          registry.fill(HIST("hSelPID"), 4.5);
        }
        if (statusPidCharmBaryon) {
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
        if (!statusInvMassCharmBaryon) {
          registry.fill(HIST("hSelPID"), 10.5);
        }
        if (statusInvMassCharmBaryon) {
          registry.fill(HIST("hSelPID"), 11.5);
        }
      }

      if (statusPidLambda && statusPidCascade && statusPidCharmBaryon && statusInvMassLambda && statusInvMassCascade && statusInvMassCharmBaryon && resultSelections) {
        hInvMassCharmBaryon->Fill(invMassCharmBaryon);
      }
    }
  } // end process

  void processOmegac0SelectorWithKFParticle(soa::Join<aod::HfCandToOmegaPi, aod::HfOmegacKf> const& candidates,
                                            TracksSel const& tracks,
                                            TracksSelLf const& lfTracks)
  {
    runOmegac0Selector<hf_cand_casc_lf::ConstructMethod::KfParticle>(candidates, tracks, lfTracks);
  }
  PROCESS_SWITCH(HfCandidateSelectorToOmegaPi, processOmegac0SelectorWithKFParticle, "Run Omegac0 to Omega pi selector with both DCA and KFParticle related selection.", false);

  void processOmegac0SelectorWithDCAFitter(aod::HfCandToOmegaPi const& candidates,
                                           TracksSel const& tracks,
                                           TracksSelLf const& lfTracks)
  {
    runOmegac0Selector<hf_cand_casc_lf::ConstructMethod::DcaFitter>(candidates, tracks, lfTracks);
  }
  PROCESS_SWITCH(HfCandidateSelectorToOmegaPi, processOmegac0SelectorWithDCAFitter, "Run Omegac0 to Omega pi selector with only DCA related selection.", true);

}; // end struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateSelectorToOmegaPi>(cfgc)};
}
