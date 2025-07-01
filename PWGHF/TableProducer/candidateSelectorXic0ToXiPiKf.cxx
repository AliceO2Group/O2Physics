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
  Produces<aod::HfMlToXiPiKf> hfMlToXiPi;

  // kinematic selections
  Configurable<double> etaTrackCharmBachMax{"etaTrackCharmBachMax", 0.8, "Max absolute value of eta for charm baryon bachelor"};
  Configurable<double> etaTrackLFDauMax{"etaTrackLFDauMax", 0.8, "Max absolute value of eta for V0 and cascade daughters"};
  Configurable<double> ptPiFromCascMin{"ptPiFromCascMin", 0.15, "Min pT pion <- casc"};

  // minimum radius cut (LFcut)
  Configurable<double> radiusCascMin{"radiusCascMin", 0.5, "Min cascade radius"};
  Configurable<double> radiusV0Min{"radiusV0Min", 1.1, "Min V0 radius"};

  // cosPA
  Configurable<double> cosPACascToPvMin{"cosPACascToPvMin", 0.99, "Min value CosPA cascade to primary vertex"};
  Configurable<double> cosPAV0ToPvMin{"cosPAV0ToPvMin", 0.97, "Min value CosPA V0 to ptimary vertex"};
  Configurable<double> cosPACascToXicMin{"cosPACascToXicMin", 0.99, "Min value CosPA cascade to charm baryon"};
  Configurable<double> cosPAV0ToCascMin{"cosPAV0ToCascMin", 0.99, "Min value CosPA V0 to cascade"};

  // DCA
  Configurable<double> dcaCharmBaryonDauMax{"dcaCharmBaryonDauMax", 0.1, "Max DCA charm baryon daughters"};
  Configurable<double> dcaCascDauMax{"dcaCascDauMax", 0.2, "Max DCA cascade daughters"};
  Configurable<double> dcaV0DauMax{"dcaV0DauMax", 1.0, "Max DCA V0 daughters"};
  Configurable<float> dcaBachToPvMin{"dcaBachToPvMin", 0.04, "DCA Bach To PV"};
  Configurable<float> dcaNegToPvMin{"dcaNegToPvMin", 0.06, "DCA Neg To PV"};
  Configurable<float> dcaPosToPvMin{"dcaPosToPvMin", 0.06, "DCA Pos To PV"};
  Configurable<float> kfDcaXYPiFromXicMax{"kfDcaXYPiFromXicMax", 0.05, "DCA Pion From Xic To PV"};
  Configurable<float> kfDcaXYCascToPvMax{"kfDcaXYCascToPvMax", 0.3, "DCA cascade To PV"};
  // Chi2Geo
  Configurable<float> chi2GeoXicMax{"chi2GeoXicMax", 70, "Geometrical Chi2 of Xic0"};
  Configurable<float> chi2GeoCascMax{"chi2GeoCascMax", 60, "Geometrical Chi2 of cascade"};
  Configurable<float> chi2GeoV0Max{"chi2GeoV0Max", 100, "Geometrical Chi2 of V0"};
  // Chi2Topo
  Configurable<float> chi2TopoXicToPvMax{"chi2TopoXicToPvMax", 120, "Topological Chi2 of Xic0 to primary vertex"};
  Configurable<float> chi2TopoPiFromXicToPvMax{"chi2TopoPiFromXicToPvMax", 250, "Topological Chi2 of pion <-- Xic0 to primary vertex"};
  Configurable<float> chi2TopoCascToPvMax{"chi2TopoCascToPvMax", 250, "Topological Chi2 of cascade to primary vertex"};
  Configurable<float> chi2TopoV0ToPvMin{"chi2TopoV0ToPvMin", 0.4, "Topological Chi2 of V0 to primary vertex"};
  Configurable<float> chi2TopoV0ToCascMax{"chi2TopoV0ToCascMax", 100, "Topological Chi2 of V0 to cascade"};
  Configurable<float> chi2TopoCascToXicMax{"chi2TopoCascToXicMax", 300, "Topological Chi2 of cascade to Xic0"};
  // ldl
  Configurable<float> cascldlMin{"cascldlMin", 0.0, "cascade ldl"};
  Configurable<float> v0ldlMin{"v0ldlMin", 0.0, "V0 ldl"};
  // decay length
  Configurable<float> decayLenXYXicMax{"decayLenXYXicMax", 1.5, "Xic0 decay length"};
  Configurable<float> decayLenXYCascMin{"decayLenXYCascMin", 0.0, "cascade decay length"};
  Configurable<float> decayLenXYLambdaMin{"decayLenXYLambdaMin", 0.0, "cascade decay length"};

  // ctau
  Configurable<float> ctXicMax{"ctXicMax", 0.4, "Xic0 ctau"};

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
  std::vector<float> outputMlXic0ToXiPiKf = {};
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
    registry.add("hCandSel", "hCandSel;selection type;entries", {HistType::kTH1D, {{28, 0., 28.}}});
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
      double etaPiFromCasc = candidate.etaBachFromCasc();
      double etaPiFromCharmBaryon = candidate.etaBachFromCharmBaryon();
      if (std::abs(etaV0PosDau) > etaTrackLFDauMax || std::abs(etaV0NegDau) > etaTrackLFDauMax || std::abs(etaPiFromCasc) > etaTrackLFDauMax || std::abs(etaPiFromCharmBaryon) > etaTrackCharmBachMax) {
        resultSelections = false;
      }
      double ptPiFromCasc = RecoDecay::sqrtSumOfSquares(candidate.pxBachFromCasc(), candidate.pyBachFromCasc());
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

      double invMassLambda = candidate.invMassLambda();
      double invMassCascade = candidate.invMassCascade();
      double invMassCharmBaryon = candidate.invMassCharmBaryon();

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

    int pTBin = findBin(binsPt, candpT);
    if (pTBin == -1) {
      return false;
    }

    double ptPiFromCharmBaryon = RecoDecay::sqrtSumOfSquares(candidate.pxBachFromCharmBaryon(), candidate.pyBachFromCharmBaryon());
    if (ptPiFromCharmBaryon <= cuts->get(pTBin, "pT pi from Xic")) {
      return false;
    }

    // fill all candidates before preselection
    registry.fill(HIST("hCandSel"), 0.5);

    // cosPA (LFcut)
    if (candidate.cosPACasc() < cosPACascToPvMin) {
      return false;
    } else {
      registry.fill(HIST("hCandSel"), 1.5);
    }
    if (candidate.cosPAV0() < cosPAV0ToPvMin) {
      return false;
    } else {
      registry.fill(HIST("hCandSel"), 2.5);
    }
    if (candidate.cosPaCascToXic() < cosPACascToXicMin) {
      return false;
    } else {
      registry.fill(HIST("hCandSel"), 3.5);
    }
    if (candidate.cosPaV0ToCasc() < cosPAV0ToCascMin) {
      return false;
    } else {
      registry.fill(HIST("hCandSel"), 4.5);
    }

    // dca cut
    if (candidate.dcaCharmBaryonDau() > dcaCharmBaryonDauMax) {
      return false;
    } else {
      registry.fill(HIST("hCandSel"), 5.5);
    }
    if (candidate.dcaCascDau() > dcaCascDauMax) {
      return false;
    } else {
      registry.fill(HIST("hCandSel"), 6.5);
    }

    if (candidate.dcaV0Dau() > dcaV0DauMax) {
      return false;
    } else {
      registry.fill(HIST("hCandSel"), 7.5);
    }

    // dcaXY pion <-- cascade to PV cut
    if (std::abs(candidate.dcaXYToPvCascDau()) < dcaBachToPvMin) {
      return false;
    } else {
      registry.fill(HIST("hCandSel"), 8.5);
    }

    // dcaXY v0 daughters to PV cut
    if (std::abs(candidate.dcaXYToPvV0Dau0()) < dcaPosToPvMin || std::abs(candidate.dcaXYToPvV0Dau1()) < dcaNegToPvMin) {
      return false;
    } else {
      registry.fill(HIST("hCandSel"), 9.5);
    }

    // dacXY pion <-- Xic0 to PV cut
    if (std::abs(candidate.kfDcaXYPiFromXic()) > kfDcaXYPiFromXicMax) {
      return false;
    } else {
      registry.fill(HIST("hCandSel"), 10.5);
    }

    // dacXY cascade to PV cut
    if (std::abs(candidate.kfDcaXYCascToPv()) > kfDcaXYCascToPvMax) {
      return false;
    } else {
      registry.fill(HIST("hCandSel"), 11.5);
    }

    // Chi2Geo
    if (candidate.chi2GeoXic() < 0 || candidate.chi2GeoXic() > chi2GeoXicMax) {
      return false;
    } else {
      registry.fill(HIST("hCandSel"), 12.5);
    }
    if (candidate.chi2GeoCasc() < 0 || candidate.chi2GeoCasc() > chi2GeoCascMax) {
      return false;
    } else {
      registry.fill(HIST("hCandSel"), 13.5);
    }
    if (candidate.chi2GeoV0() < 0 || candidate.chi2GeoV0() > chi2GeoV0Max) {
      return false;
    } else {
      registry.fill(HIST("hCandSel"), 14.5);
    }

    // Chi2Topo
    if (candidate.chi2TopoXicToPv() < 0 || candidate.chi2TopoXicToPv() > chi2TopoXicToPvMax) {
      return false;
    } else {
      registry.fill(HIST("hCandSel"), 15.5);
    }
    if (candidate.chi2TopoPiFromXicToPv() < 0 || candidate.chi2TopoPiFromXicToPv() > chi2TopoPiFromXicToPvMax) {
      return false;
    } else {
      registry.fill(HIST("hCandSel"), 16.5);
    }
    if (candidate.chi2TopoCascToPv() < 0 || candidate.chi2TopoCascToPv() > chi2TopoCascToPvMax) {
      return false;
    } else {
      registry.fill(HIST("hCandSel"), 17.5);
    }
    if (candidate.chi2TopoV0ToPv() > 0 && candidate.chi2TopoV0ToPv() < chi2TopoV0ToPvMin) {
      return false;
    } else {
      registry.fill(HIST("hCandSel"), 18.5);
    }
    if (candidate.chi2TopoV0ToCasc() < 0 || candidate.chi2TopoV0ToCasc() > chi2TopoV0ToCascMax) {
      return false;
    } else {
      registry.fill(HIST("hCandSel"), 19.5);
    }
    if (candidate.chi2TopoCascToXic() < 0 || candidate.chi2TopoCascToXic() > chi2TopoCascToXicMax) {
      return false;
    } else {
      registry.fill(HIST("hCandSel"), 20.5);
    }

    // ldl
    if (candidate.cascldl() < cascldlMin) {
      return false;
    } else {
      registry.fill(HIST("hCandSel"), 21.5);
    }
    if (candidate.v0ldl() < v0ldlMin) {
      return false;
    } else {
      registry.fill(HIST("hCandSel"), 22.5);
    }
    // decay length
    if (std::abs(candidate.decayLenXYXic()) > decayLenXYXicMax) {
      return false;
    } else {
      registry.fill(HIST("hCandSel"), 23.5);
    }
    if (std::abs(candidate.decayLenXYCasc()) < decayLenXYCascMin) {
      return false;
    } else {
      registry.fill(HIST("hCandSel"), 24.5);
    }
    if (std::abs(candidate.decayLenXYLambda()) < decayLenXYLambdaMin) {
      return false;
    } else {
      registry.fill(HIST("hCandSel"), 25.5);
    }
    // ctau
    if (std::abs(candidate.cTauXic()) > ctXicMax) {
      return false;
    } else {
      registry.fill(HIST("hCandSel"), 26.5);
    }
    return true;
  }

}; // end struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateSelectorXic0ToXiPiKf>(cfgc)};
}
