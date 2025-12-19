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

/// \file candidateSelectorXic0ToXiPiKf.cxx
/// \brief Xic0 → Xi Pi selection task
/// \author Ran Tu <ran.tu@cern.ch>, Fudan University
/// \author Tao Fang <tao.fang@cern.ch>, Central China Normal University

#include "PWGHF/Core/HfMlResponseXic0ToXiPiKf.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsAnalysis.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelectorPID.h"

#include <CCDB/CcdbApi.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>

#include <Rtypes.h>

#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::analysis;

enum PidInfoStored {
  PiFromLam = 0,
  PrFromLam,
  PiFromCasc,
  PiFromCharm
};

/// Struct for applying Xic0 -> Xi pi selection cuts
struct HfCandidateSelectorXic0ToXiPiKf {
  Produces<aod::HfSelToXiPiKf> hfSelToXiPi;
  Produces<aod::HfMlToXiPi> hfMlToXiPi;

  // kinematic selections
  Configurable<double> etaTrackCharmBachMax{"etaTrackCharmBachMax", 0.8, "Max absolute value of eta for charm baryon bachelor"};
  Configurable<double> etaTrackLFDauMax{"etaTrackLFDauMax", 0.8, "Max absolute value of eta for V0 and cascade daughters"};
  Configurable<double> ptPiFromCascMin{"ptPiFromCascMin", 0.15, "Min pT pion <- casc"};

  // minimum radius cut (LFcut)
  Configurable<double> radiusCascMin{"radiusCascMin", 0.5, "Min cascade radius"};
  Configurable<double> radiusV0Min{"radiusV0Min", 1.1, "Min V0 radius"};

  Configurable<float> v0MassWindow{"v0MassWindow", 0.01, "V0 mass window"};
  Configurable<float> cascadeMassWindow{"cascadeMassWindow", 0.01, "Cascade mass window"};
  Configurable<bool> applyTrkSelLf{"applyTrkSelLf", true, "Apply track selection for LF daughters"};

  // limit charm baryon invariant mass spectrum
  Configurable<double> invMassCharmBaryonMin{"invMassCharmBaryonMin", 2.0, "Lower limit invariant mass spectrum charm baryon"}; // 2.4 Omegac0 only
  Configurable<double> invMassCharmBaryonMax{"invMassCharmBaryonMax", 3.1, "Upper limit invariant mass spectrum charm baryon"};

  // PID options
  Configurable<bool> usePidTpcOnly{"usePidTpcOnly", false, "Perform PID using only TPC"};
  Configurable<bool> usePidTpcTofCombined{"usePidTpcTofCombined", true, "Perform PID using TPC & TOF"};

  // PID - TPC selections

  Configurable<double> ptPrPidTpcMin{"ptPrPidTpcMin", -1, "Lower bound of track pT for TPC PID for proton selection"};
  Configurable<double> ptPrPidTpcMax{"ptPrPidTpcMax", 9999.9, "Upper bound of track pT for TPC PID for proton selection"};
  Configurable<double> nSigmaTpcPrMax{"nSigmaTpcPrMax", 3., "Nsigma cut on TPC only for proton selection"};
  Configurable<double> nSigmaTpcCombinedPrMax{"nSigmaTpcCombinedPrMax", 0., "Nsigma cut on TPC combined with TOF for proton selection"};

  Configurable<double> ptPiPidTpcMin{"ptPiPidTpcMin", -1, "Lower bound of track pT for TPC PID for pion selection"};
  Configurable<double> ptPiPidTpcMax{"ptPiPidTpcMax", 9999.9, "Upper bound of track pT for TPC PID for pion selection"};
  Configurable<double> nSigmaTpcPiMax{"nSigmaTpcPiMax", 3., "Nsigma cut on TPC only for pion selection"};
  Configurable<double> nSigmaTpcCombinedPiMax{"nSigmaTpcCombinedPiMax", 0., "Nsigma cut on TPC combined with TOF for pion selection"};

  // PID - TOF selections

  Configurable<double> ptPrPidTofMin{"ptPrPidTofMin", -1, "Lower bound of track pT for TOF PID for proton selection"};
  Configurable<double> ptPrPidTofMax{"ptPrPidTofMax", 9999.9, "Upper bound of track pT for TOF PID for proton selection"};
  Configurable<double> nSigmaTofPrMax{"nSigmaTofPrMax", 5., "Nsigma cut on TOF only for proton selection"};
  Configurable<double> nSigmaTofCombinedPrMax{"nSigmaTofCombinedPrMax", 0., "Nsigma cut on TOF combined with TPC for proton selection"};

  Configurable<double> ptPiPidTofMin{"ptPiPidTofMin", -1, "Lower bound of track pT for TOF PID for pion selection"};
  Configurable<double> ptPiPidTofMax{"ptPiPidTofMax", 9999.9, "Upper bound of track pT for TOF PID for pion selection"};
  Configurable<double> nSigmaTofPiMax{"nSigmaTofPiMax", 5., "Nsigma cut on TOF only for pion selection"};
  Configurable<double> nSigmaTofCombinedPiMax{"nSigmaTofCombinedPiMax", 0., "Nsigma cut on TOF combined with TOF for pion selection"};

  // detector clusters selections
  Configurable<int> nClustersTpcMin{"nClustersTpcMin", 70, "Minimum number of TPC clusters requirement"};
  Configurable<int> nTpcCrossedRowsMin{"nTpcCrossedRowsMin", 70, "Minimum number of TPC crossed rows requirement"};
  Configurable<double> tpcCrossedRowsOverFindableClustersRatioMin{"tpcCrossedRowsOverFindableClustersRatioMin", 0.8, "Minimum ratio TPC crossed rows over findable clusters requirement"};
  Configurable<float> tpcChi2PerClusterMax{"tpcChi2PerClusterMax", 4, "Maximum value of chi2 fit over TPC clusters"};
  Configurable<int> nClustersItsMin{"nClustersItsMin", 3, "Minimum number of ITS clusters requirement for pi <- charm baryon"};
  Configurable<int> nClustersItsInnBarrMin{"nClustersItsInnBarrMin", 1, "Minimum number of ITS clusters in inner barrel requirement for pi <- charm baryon"};
  Configurable<float> itsChi2PerClusterMax{"itsChi2PerClusterMax", 36, "Maximum value of chi2 fit over ITS clusters for pi <- charm baryon"};

  // topological cuts
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_xic_to_xi_pi::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_xic_to_xi_pi::Cuts[0], hf_cuts_xic_to_xi_pi::NBinsPt, hf_cuts_xic_to_xi_pi::NCutVars, hf_cuts_xic_to_xi_pi::labelsPt, hf_cuts_xic_to_xi_pi::labelsCutVar}, "Xic0 candidate selection per pT bin"};

  // ML inference
  Configurable<bool> applyMl{"applyMl", true, "Flag to apply ML selections"};
  Configurable<std::vector<double>> binsPtMl{"binsPtMl", std::vector<double>{hf_cuts_ml::vecBinsPt}, "pT bin limits for ML application"};
  Configurable<std::vector<int>> cutDirMl{"cutDirMl", std::vector<int>{hf_cuts_ml::vecCutDir}, "Whether to reject score values greater or smaller than the threshold"};
  Configurable<LabeledArray<double>> cutsMl{"cutsMl", {hf_cuts_ml::Cuts[0], hf_cuts_ml::NBinsPt, hf_cuts_ml::NCutScores, hf_cuts_ml::labelsPt, hf_cuts_ml::labelsCutScore}, "ML selections per pT bin"};
  Configurable<int> nClassesMl{"nClassesMl", static_cast<int>(hf_cuts_ml::NCutScores), "Number of classes in ML model"};
  Configurable<std::vector<std::string>> namesInputFeatures{"namesInputFeatures", std::vector<std::string>{"feature1", "feature2"}, "Names of ML model input features"};

  // CCDB configuration
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::vector<std::string>> modelPathsCCDB{"modelPathsCCDB", std::vector<std::string>{"EventFiltering/PWGHF/BDTXic0ToXipiKf"}, "Paths of models on CCDB"};
  Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"ModelHandler_onnx_Xic0ToXipiKf.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};

  o2::analysis::HfMlResponseXic0ToXiPiKf<float> hfMlResponse;
  std::vector<float> outputMlXic0ToXiPiKf;
  o2::ccdb::CcdbApi ccdbApi;

  TrackSelectorPr selectorProton;
  TrackSelectorPi selectorPion;

  using TracksSel = soa::Join<aod::TracksWDcaExtra, aod::TracksPidPi, aod::TracksPidPr, aod::TracksPidPi>;
  using TracksSelLf = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksPidPi, aod::TracksPidPr, aod::TracksPidPi>;

  HistogramRegistry registry{"registry"}; // for QA of selections

  OutputObj<TH1D> hInvMassCharmBaryon{TH1D("hInvMassCharmBaryon", "Charm baryon invariant mass;inv mass;entries", 500, 2.3, 3.1)};

  void init(InitContext const&)
  {
    selectorProton.setRangePtTpc(ptPrPidTpcMin, ptPrPidTpcMax);
    selectorProton.setRangeNSigmaTpc(-nSigmaTpcPrMax, nSigmaTpcPrMax);
    selectorProton.setRangeNSigmaTpcCondTof(-nSigmaTpcCombinedPrMax, nSigmaTpcCombinedPrMax);
    selectorProton.setRangePtTof(ptPrPidTofMin, ptPrPidTofMax);
    selectorProton.setRangeNSigmaTof(-nSigmaTofPrMax, nSigmaTofPrMax);
    selectorProton.setRangeNSigmaTofCondTpc(-nSigmaTofCombinedPrMax, nSigmaTofCombinedPrMax);

    selectorPion.setRangePtTpc(ptPiPidTpcMin, ptPiPidTpcMax);
    selectorPion.setRangeNSigmaTpc(-nSigmaTpcPiMax, nSigmaTpcPiMax);
    selectorPion.setRangeNSigmaTpcCondTof(-nSigmaTpcCombinedPiMax, nSigmaTpcCombinedPiMax);
    selectorPion.setRangePtTof(ptPiPidTofMin, ptPiPidTofMax);
    selectorPion.setRangeNSigmaTof(-nSigmaTofPiMax, nSigmaTofPiMax);
    selectorPion.setRangeNSigmaTofCondTpc(-nSigmaTofCombinedPiMax, nSigmaTofCombinedPiMax);

    const AxisSpec axisSel{2, -0.5, 1.5, "status"};
    registry.add("hStatusCheck", "Check consecutive selections status;status;entries", {HistType::kTH1D, {{3, 0., 3.}}});
    // sign of candidates
    registry.add("hSelSignDec", "hSelSignDec;status;entries", {HistType::kTH1D, {axisSel}});
    // basic selections (bin 0 -> candidates that did not pass the selection, bin 1 -> candidates that passed the selection)
    registry.add("hSelPtPiFromCasc", "hSelPtPiFromCasc;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelRadCasc", "hSelRadCasc;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelRadV0", "hSelRadV0;status;entries", {HistType::kTH1D, {axisSel}});

    // TPC and ITS selections (bin 0 -> candidates that did not pass the selection, bin 1 -> candidates that passed the selection)
    registry.add("hSelTPCQualityPiFromCharm", "hSelTPCQualityPiFromCharm;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelTPCQualityPiFromLam", "hSelTPCQualityPiFromLam;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelTPCQualityPrFromLam", "hSelTPCQualityPrFromLam;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelTPCQualityPiFromCasc", "hSelTPCQualityPiFromCasc;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelITSQualityPiFromCharm", "hSelITSQualityPiFromCharm;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelPID", "hSelPID;status;entries", {HistType::kTH1D, {{12, 0., 12.}}});
    registry.add("hSelMassLam", "hSelMassLam;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelMassCasc", "hSelMassCasc;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelMassCharmBaryon", "hSelMassCharmBaryon;status;entries", {HistType::kTH1D, {axisSel}});
    registry.add("hSelMlXic0", "hSelMlXic0;status;entries", {HistType::kTH1D, {axisSel}});

    if (applyMl) {
      hfMlResponse.configure(binsPtMl, cutsMl, cutDirMl, nClassesMl);
      if (loadModelsFromCCDB) {
        ccdbApi.init(ccdbUrl);
        hfMlResponse.setModelPathsCCDB(onnxFileNames, ccdbApi, modelPathsCCDB, timestampCCDB);
      } else {
        hfMlResponse.setModelPathsLocal(onnxFileNames);
      }
      hfMlResponse.cacheInputFeaturesIndices(namesInputFeatures);
      hfMlResponse.init();
    }
  }

  void process(aod::HfCandToXiPiKf const& candidates,
               TracksSel const& tracks,
               TracksSelLf const& lfTracks)
  {

    // looping over charm baryon candidates
    for (const auto& candidate : candidates) {

      outputMlXic0ToXiPiKf.clear();

      auto ptCand = RecoDecay::pt(candidate.pxCharmBaryon(), candidate.pyCharmBaryon());

      bool resultSelections = true; // True if the candidate passes all the selections, False otherwise

      auto trackV0PosDauId = candidate.posTrackId();                   // positive V0 daughter
      auto trackV0NegDauId = candidate.negTrackId();                   // negative V0 daughter
      auto trackPiFromCascId = candidate.bachelorId();                 // pion <- cascade
      auto trackPiFromCharmId = candidate.bachelorFromCharmBaryonId(); // pion <- charm baryon
      auto trackV0PosDau = lfTracks.rawIteratorAt(trackV0PosDauId);
      auto trackV0NegDau = lfTracks.rawIteratorAt(trackV0NegDauId);
      auto trackPiFromCasc = lfTracks.rawIteratorAt(trackPiFromCascId);
      auto trackPiFromCharm = tracks.rawIteratorAt(trackPiFromCharmId);

      auto trackPiFromLam = trackV0NegDau;
      auto trackPrFromLam = trackV0PosDau;

      int8_t const signDecay = candidate.signDecay(); // sign of pi <- cascade

      if (signDecay > 0) {
        trackPiFromLam = trackV0PosDau;
        trackPrFromLam = trackV0NegDau;
        registry.fill(HIST("hSelSignDec"), 1); // anti-particle decay
      } else if (signDecay < 0) {
        registry.fill(HIST("hSelSignDec"), 0); // particle decay
      }

      // eta selection
      double const etaV0PosDau = candidate.etaV0PosDau();
      double const etaV0NegDau = candidate.etaV0NegDau();
      double const etaPiFromCasc = candidate.etaBachFromCasc();
      double const etaPiFromCharmBaryon = candidate.etaBachFromCharmBaryon();
      if (std::abs(etaV0PosDau) > etaTrackLFDauMax || std::abs(etaV0NegDau) > etaTrackLFDauMax || std::abs(etaPiFromCasc) > etaTrackLFDauMax || std::abs(etaPiFromCharmBaryon) > etaTrackCharmBachMax) {
        resultSelections = false;
      }
      double const ptPiFromCasc = RecoDecay::sqrtSumOfSquares(candidate.pxBachFromCasc(), candidate.pyBachFromCasc());
      if (std::abs(ptPiFromCasc) < ptPiFromCascMin) {
        resultSelections = false;
        registry.fill(HIST("hSelPtPiFromCasc"), 0);
      } else {
        registry.fill(HIST("hSelPtPiFromCasc"), 1);
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
        if (!isSelectedTrackTpcQuality(trackPiFromCasc, nClustersTpcMin, nTpcCrossedRowsMin, tpcCrossedRowsOverFindableClustersRatioMin, tpcChi2PerClusterMax)) {
          resultSelections = false;
          registry.fill(HIST("hSelTPCQualityPiFromCasc"), 0);
        } else {
          registry.fill(HIST("hSelTPCQualityPiFromCasc"), 1);
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
      int statusPidPiFromCasc = -999;
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
        SETBIT(infoTpcStored, PiFromLam);
      }
      if (trackPrFromLam.hasTPC()) {
        SETBIT(infoTpcStored, PiFromLam);
      }
      if (trackPiFromCasc.hasTPC()) {
        SETBIT(infoTpcStored, PiFromCasc);
      }
      if (trackPiFromCharm.hasTPC()) {
        SETBIT(infoTpcStored, PiFromCharm);
      }
      if (trackPiFromLam.hasTOF()) {
        SETBIT(infoTofStored, PiFromLam);
      }
      if (trackPrFromLam.hasTOF()) {
        SETBIT(infoTofStored, PiFromLam);
      }
      if (trackPiFromCasc.hasTOF()) {
        SETBIT(infoTofStored, PiFromCasc);
      }
      if (trackPiFromCharm.hasTOF()) {
        SETBIT(infoTofStored, PiFromCharm);
      }

      if (usePidTpcOnly) {
        statusPidPrFromLam = selectorProton.statusTpc(trackPrFromLam);
        statusPidPiFromLam = selectorPion.statusTpc(trackPiFromLam);
        statusPidPiFromCasc = selectorPion.statusTpc(trackPiFromCasc);
        statusPidPiFromCharmBaryon = selectorPion.statusTpc(trackPiFromCharm);
      } else if (usePidTpcTofCombined) {
        statusPidPrFromLam = selectorProton.statusTpcOrTof(trackPrFromLam);
        statusPidPiFromLam = selectorPion.statusTpcOrTof(trackPiFromLam);
        statusPidPiFromCasc = selectorPion.statusTpcOrTof(trackPiFromCasc);
        statusPidPiFromCharmBaryon = selectorPion.statusTpcOrTof(trackPiFromCharm);
      }

      if (statusPidPrFromLam == TrackSelectorPID::Accepted && statusPidPiFromLam == TrackSelectorPID::Accepted) {
        statusPidLambda = true;
        if (resultSelections) {
          registry.fill(HIST("hStatusCheck"), 0.5);
        }
      }

      if (statusPidPrFromLam == TrackSelectorPID::Accepted && statusPidPiFromLam == TrackSelectorPID::Accepted && statusPidPiFromCasc == TrackSelectorPID::Accepted) {
        statusPidCascade = true;
        if (resultSelections) {
          registry.fill(HIST("hStatusCheck"), 1.5);
        }
      }

      if (statusPidPrFromLam == TrackSelectorPID::Accepted && statusPidPiFromLam == TrackSelectorPID::Accepted && statusPidPiFromCasc == TrackSelectorPID::Accepted && statusPidPiFromCharmBaryon == TrackSelectorPID::Accepted) {
        statusPidCharmBaryon = true;
      } else {
        registry.fill(HIST("hStatusCheck"), 2.5);
        resultSelections = false;
      }

      // invariant mass cuts
      bool statusInvMassLambda = false;
      bool statusInvMassCascade = false;
      bool statusInvMassCharmBaryon = false;

      double const invMassLambda = candidate.invMassLambda();
      double const invMassCascade = candidate.invMassCascade();
      double const invMassCharmBaryon = candidate.invMassCharmBaryon();

      if (resultSelections) {
        resultSelections = selectionTopolKf(candidate);
        if (std::abs(invMassLambda - o2::constants::physics::MassLambda0) < v0MassWindow) {
          statusInvMassLambda = true;
          registry.fill(HIST("hSelMassLam"), 1);
        } else {
          resultSelections = false;
          registry.fill(HIST("hSelMassLam"), 0);
        }

        if (std::abs(invMassCascade - o2::constants::physics::MassXiMinus) < cascadeMassWindow) {
          statusInvMassCascade = true;
          registry.fill(HIST("hSelMassCasc"), 1);
        } else {
          resultSelections = false;
          registry.fill(HIST("hSelMassCasc"), 0);
        }

        if ((invMassCharmBaryon >= invMassCharmBaryonMin) && (invMassCharmBaryon <= invMassCharmBaryonMax)) {
          statusInvMassCharmBaryon = true;
          registry.fill(HIST("hSelMassCharmBaryon"), 1);
        } else {
          resultSelections = false;
          registry.fill(HIST("hSelMassCharmBaryon"), 0);
        }
      }
      // ML selections
      if (applyMl) {
        bool isSelectedMlXic0 = false;
        std::vector<float> inputFeaturesXic0 = hfMlResponse.getInputFeatures(candidate, trackPiFromLam, trackPiFromCasc, trackPiFromCharm);
        if (!resultSelections) {
          hfMlToXiPi(outputMlXic0ToXiPiKf);
        } else {
          isSelectedMlXic0 = hfMlResponse.isSelectedMl(inputFeaturesXic0, ptCand, outputMlXic0ToXiPiKf);
          registry.fill(HIST("hSelMlXic0"), isSelectedMlXic0);
          hfMlToXiPi(outputMlXic0ToXiPiKf);
        }
      }

      hfSelToXiPi(resultSelections,
                  trackPiFromCharm.tpcNSigmaPi(), trackPiFromCasc.tpcNSigmaPi(), trackPiFromLam.tpcNSigmaPi(), trackPrFromLam.tpcNSigmaPr(),
                  trackPiFromCharm.tofNSigmaPi(), trackPiFromCasc.tofNSigmaPi(), trackPiFromLam.tofNSigmaPi(), trackPrFromLam.tofNSigmaPr());

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

  /// \param candidate is candidate
  /// \return true if candidate passes all cuts
  template <typename T>
  bool selectionTopolKf(const T& candidate)
  {
    auto candpT = RecoDecay::pt(candidate.pxCharmBaryon(), candidate.pyCharmBaryon());

    int const pTBin = findBin(binsPt, candpT);
    if (pTBin == -1) {
      return false;
    }

    double const ptPiFromCharmBaryon = RecoDecay::sqrtSumOfSquares(candidate.pxBachFromCharmBaryon(), candidate.pyBachFromCharmBaryon());
    if (ptPiFromCharmBaryon <= cuts->get(pTBin, "ptPiFromCharmBaryon")) {
      return false;
    }

    // cosPA (LFcut)
    if (candidate.cosPACasc() < cuts->get(pTBin, "cosPACasc")) {
      return false;
    }
    if (candidate.cosPAV0() < cuts->get(pTBin, "cosPAV0")) {
      return false;
    }
    if (candidate.cosPaCascToXic() < cuts->get(pTBin, "cosPaCascToXic")) {
      return false;
    }
    if (candidate.cosPaV0ToCasc() < cuts->get(pTBin, "cosPaV0ToCasc")) {
      return false;
    }

    // dca cut
    if (candidate.dcaCharmBaryonDau() > cuts->get(pTBin, "dcaCharmBaryonDau")) {
      return false;
    }
    if (candidate.dcaCascDau() > cuts->get(pTBin, "dcaCascDau")) {
      return false;
    }

    if (candidate.dcaV0Dau() > cuts->get(pTBin, "dcaV0Dau")) {
      return false;
    }

    // dcaXY pion <-- cascade to PV cut
    if (std::abs(candidate.dcaXYToPvCascDau()) < cuts->get(pTBin, "dcaXYToPvCascDau")) {
      return false;
    }

    // dcaXY v0 daughters to PV cut
    if (std::abs(candidate.dcaXYToPvV0Dau0()) < cuts->get(pTBin, "dcaXYToPvV0Dau0") || std::abs(candidate.dcaXYToPvV0Dau1()) < cuts->get(pTBin, "dcaXYToPvV0Dau0")) {
      return false;
    }

    // dacXY pion <-- Xic0 to PV cut
    if (std::abs(candidate.kfDcaXYPiFromXic()) > cuts->get(pTBin, "kfDcaXYPiFromXic")) {
      return false;
    }

    // dacXY cascade to PV cut
    if (std::abs(candidate.kfDcaXYCascToPv()) > cuts->get(pTBin, "kfDcaXYCascToPv")) {
      return false;
    }

    // Chi2Geo
    if (candidate.chi2GeoXic() < 0 || candidate.chi2GeoXic() > cuts->get(pTBin, "chi2GeoXic")) {
      return false;
    }
    if (candidate.chi2GeoCasc() < 0 || candidate.chi2GeoCasc() > cuts->get(pTBin, "chi2GeoCasc")) {
      return false;
    }
    if (candidate.chi2GeoV0() < 0 || candidate.chi2GeoV0() > cuts->get(pTBin, "chi2GeoV0")) {
      return false;
    }

    // Chi2Topo
    if (candidate.chi2TopoXicToPv() < 0 || candidate.chi2TopoXicToPv() > cuts->get(pTBin, "chi2TopoXicToPv")) {
      return false;
    }
    if (candidate.chi2TopoPiFromXicToPv() < 0 || candidate.chi2TopoPiFromXicToPv() > cuts->get(pTBin, "chi2TopoPiFromXicToPv")) {
      return false;
    }
    if (candidate.chi2TopoCascToPv() < 0 || candidate.chi2TopoCascToPv() > cuts->get(pTBin, "chi2TopoCascToPv")) {
      return false;
    }
    if (candidate.chi2TopoV0ToPv() > 0 && candidate.chi2TopoV0ToPv() < cuts->get(pTBin, "chi2TopoV0ToPv")) {
      return false;
    }
    if (candidate.chi2TopoV0ToCasc() < 0 || candidate.chi2TopoV0ToCasc() > cuts->get(pTBin, "chi2TopoV0ToCasc")) {
      return false;
    }
    if (candidate.chi2TopoCascToXic() < 0 || candidate.chi2TopoCascToXic() > cuts->get(pTBin, "chi2TopoCascToXic")) {
      return false;
    }

    // ldl
    if (candidate.cascldl() < cuts->get(pTBin, "cascldl")) {
      return false;
    }
    if (candidate.v0ldl() < cuts->get(pTBin, "v0ldl")) {
      return false;
    }
    // decay length
    if (std::abs(candidate.decayLenXYXic()) > cuts->get(pTBin, "decayLenXYXic")) {
      return false;
    }
    if (std::abs(candidate.decayLenXYCasc()) < cuts->get(pTBin, "decayLenXYCasc")) {
      return false;
    }
    if (std::abs(candidate.decayLenXYLambda()) < cuts->get(pTBin, "decayLenXYLambda")) {
      return false;
    }
    // ctau
    if (std::abs(candidate.cTauXic()) > cuts->get(pTBin, "cTauXic")) {
      return false;
    }
    return true;
  }

}; // end struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateSelectorXic0ToXiPiKf>(cfgc)};
}
