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

/// \file
/// \brief  Omegac0 → Omega Pi selection task
/// \author

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectorPID.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/HfMlResponse.h"
#include "PWGHF/Core/HfMlResponseOmegaCToOmegaPi.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::analysis;
int toCounter = 0;
enum pidInfoStored {
  kPiFromLam = 0,
  kPrFromLam,
  kKaFromCasc,
  kPiFromCharm
};

/// Struct for applying Omegac0 selection cuts
struct HfCandidateSelectorOmegacToOmegaPiWithKfp {
  Produces<aod::HfSelOmegacToOmegaPi> HfSelOmegacToOmegaPi;
  Produces<aod::HfMlSelOmegacToOmegaPi> hfMlToOmegaPiCandidate;

  // LF analysis selections
  Configurable<double> radiusCascMin{"radiusCascMin", 0.6, "Min cascade radius"};
  Configurable<double> radiusV0Min{"radiusV0Min", 1.2, "Min V0 radius"};
  Configurable<double> cosPAV0Min{"cosPAV0Min", 0.97, "Min valueCosPA V0"};
  Configurable<double> cosPACascMin{"cosPACascMin", 0.97, "Min value CosPA cascade"};
  Configurable<double> dcaCascDauMax{"dcaCascDauMax", 1.0, "Max DCA cascade daughters"};
  Configurable<double> dcaV0DauMax{"dcaV0DauMax", 1.0, "Max DCA V0 daughters"};
  Configurable<float> dcaBachToPvMin{"dcaBachToPvMin", 0.04, "DCA Bach To PV"};
  Configurable<float> dcaNegToPvMin{"dcaNegToPvMin", 0.06, "DCA Neg To PV"};
  Configurable<float> dcaPosToPvMin{"dcaPosToPvMin", 0.06, "DCA Pos To PV"};
  Configurable<float> v0MassWindow{"v0MassWindow", 0.01, "V0 mass window"};
  Configurable<float> cascadeMassWindow{"cascadeMassWindow", 0.01, "Cascade mass window"};

  // limit charm baryon invariant mass spectrum
  Configurable<double> invMassCharmBaryonMin{"invMassCharmBaryonMin", 2.0, "Lower limit invariant mass spectrum charm baryon"}; // 2.4 Omegac0 only
  Configurable<double> invMassCharmBaryonMax{"invMassCharmBaryonMax", 3.1, "Upper limit invariant mass spectrum charm baryon"};

  // kinematic selections
  Configurable<double> etaTrackCharmBachMax{"etaTrackCharmBachMax", 0.8, "Max absolute value of eta for charm baryon bachelor"};
  Configurable<double> etaTrackLFDauMax{"etaTrackLFDauMax", 1.0, "Max absolute value of eta for V0 and cascade daughters"};
  Configurable<double> ptKaFromCascMin{"ptKaFromCascMin", 0.15, "Min pT Kaon <- casc"};
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
  Configurable<double> nSigmaTpcKaMax{"nSigmaTpcKarMax", 3., "Nsigma cut on TPC only for pkaon selection"};
  Configurable<double> nSigmaTpcCombinedKaMax{"nSigmaTpcCombinedKarMax", 0., "Nsigma cut on TPC combined with TOF for kaon selection"};
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
  Configurable<double> nSigmaTofCombinedKaMax{"nSigmaTofCombinedKaMax", 0., "Nsigma cut on TOF combined with TPC for kaon selection"};

  // detector clusters selections
  Configurable<int> nClustersTpcMin{"nClustersTpcMin", 70, "Minimum number of TPC clusters requirement"};
  Configurable<int> nTpcCrossedRowsMin{"nTpcCrossedRowsMin", 70, "Minimum number of TPC crossed rows requirement"};
  Configurable<double> tpcCrossedRowsOverFindableClustersRatioMin{"tpcCrossedRowsOverFindableClustersRatioMin", 0.8, "Minimum ratio TPC crossed rows over findable clusters requirement"};
  Configurable<int> nClustersItsMin{"nClustersItsMin", 3, "Minimum number of ITS clusters requirement for pi <- charm baryon"};
  Configurable<int> nClustersItsInnBarrMin{"nClustersItsInnBarrMin", 1, "Minimum number of ITS clusters in inner barrel requirement for pi <- charm baryon"};

  // ML inference
  Configurable<bool> applyMl{"applyMl", false, "Flag to apply ML selections"};
  Configurable<std::vector<double>> binsPtMl{"binsPtMl", std::vector<double>{hf_cuts_ml::vecBinsPt}, "pT bin limits for ML application"};
  Configurable<std::vector<int>> cutDirMl{"cutDirMl", std::vector<int>{hf_cuts_ml::vecCutDir}, "Whether to reject score values greater or smaller than the threshold"};
  Configurable<LabeledArray<double>> cutsMl{"cutsMl", {hf_cuts_ml::cuts[0], hf_cuts_ml::nBinsPt, hf_cuts_ml::nCutScores, hf_cuts_ml::labelsPt, hf_cuts_ml::labelsCutScore}, "ML selections per pT bin"};
  Configurable<int8_t> nClassesMl{"nClassesMl", (int8_t)hf_cuts_ml::nCutScores, "Number of classes in ML model"};
  Configurable<bool> enableDebugMl{"enableDebugMl", false, "Flag to enable histograms to monitor BDT application"};
  Configurable<std::vector<std::string>> namesInputFeatures{"namesInputFeatures", std::vector<std::string>{"feature1", "feature2"}, "Names of ML model input features"};
  // CCDB configuration
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::vector<std::string>> modelPathsCCDB{"modelPathsCCDB", std::vector<std::string>{"EventFiltering/PWGHF/BDTD0"}, "Paths of models on CCDB"};
  Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"ModelHandler_onnx_Omegac0ToOmegaPi.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};

  o2::analysis::HfMlResponseOmegaCToOmegaPi<float> hfMlResponse;
  std::vector<float> outputMlOmegaC = {};
  std::vector<float> outputMlOmegaCbar = {};
  o2::ccdb::CcdbApi ccdbApi;

  HfHelper hfHelper;

  TrackSelectorPi selectorPion;
  TrackSelectorPr selectorProton;
  TrackSelectorKa selectorKaon;

  using TracksSel = soa::Join<aod::TracksWDcaExtra, aod::TracksPidPi, aod::TracksPidPr, aod::TracksPidKa>;
  // using TracksSel = soa::Join<aod::TracksWDcaExtra, aod::TracksPidPi, aod::TracksPidPr, aod::TracksPidKa, aod::PidTpcTofFullPi, aod::PidTpcTofFullPr, aod::PidTpcTofFullKa>;
  HistogramRegistry registry{"registry"}; // for QA of selections

  OutputObj<TH1F> hInvMassCharmBaryon{TH1F("hInvMassCharmBaryon", "Charm baryon invariant mass;inv mass;entries", 500, 2.2, 3.1)};
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

    selectorKaon.setRangePtTpc(ptKaPidTpcMin, ptKaPidTpcMax);
    selectorKaon.setRangeNSigmaTpc(-nSigmaTpcKaMax, nSigmaTpcKaMax);
    selectorKaon.setRangeNSigmaTpcCondTof(-nSigmaTpcCombinedKaMax, nSigmaTpcCombinedKaMax);
    selectorKaon.setRangePtTof(ptKaPidTofMin, ptKaPidTofMax);
    selectorKaon.setRangeNSigmaTof(-nSigmaTofKaMax, nSigmaTofKaMax);
    selectorKaon.setRangeNSigmaTofCondTpc(-nSigmaTofCombinedKaMax, nSigmaTofCombinedKaMax);

    registry.add("hSelPID", "hSelPID;status;entries", {HistType::kTH1F, {{12, 0., 12.}}});
    registry.add("hStatusCheck", "Check consecutive selections status;status;entries", {HistType::kTH1F, {{14, -2., 12.}}});

    // for QA of the selections (bin 0 -> candidates that did not pass the selection, bin 1 -> candidates that passed the selection)
    registry.add("hSelSignDec", "hSelSignDec;status;entries", {HistType::kTH1F, {{3, 0., 3.}}});
    registry.add("hSelEtaPosV0Dau", "hSelEtaPosV0Dau;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelEtaNegV0Dau", "hSelEtaNegV0Dau;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelEtaKaFromCasc", "hSelEtaKaFromCasc;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelEtaPiFromCharm", "hSelEtaPiFromCharm;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelRadCasc", "hSelRadCasc;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelRadV0", "hSelRadV0;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelCosPACasc", "hSelCosPACasc;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelCosPAV0", "hSelCosPAV0;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelDCACascDau", "hSelDCACascDau;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelDCAV0Dau", "hSelDCAV0Dau;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelDCACharmDau", "hSelDCACharmDau;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelDCAXYPrimPi", "hSelDCAXYPrimPi;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelDCAZPrimPi", "hSelDCAZPrimPi;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelDCAXYCasc", "hSelDCAXYCasc;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelDCAZCasc", "hSelDCAZCasc;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelPtKaFromCasc", "hSelPtKaFromCasc;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelPtPiFromCharm", "hSelPtPiFromCharm;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelNClsTPCPiFromCharm", "hSelNClsTPCPiFromCharm;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelNClsTPCPiFromLam", "hSelNClsTPCPiFromLam;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelNClsTPCPrFromLam", "hSelNClsTPCPrFromLam;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelNClsTPCKaFromCasc", "hSelNClsTPCKaFromCasc;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelNCrossRowsTPCPiFromCharm", "hSelNCrossRowsTPCPiFromCharm;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelNCrossRowsTPCPiFromLam", "hSelNCrossRowsTPCPiFromLam;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelNCrossRowsTPCPrFromLam", "hSelNCrossRowsTPCPrFromLam;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelNCrossRowsTPCKaFromCasc", "hSelNCrossRowsTPCKaFromCasc;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelCrossRowsOverFindClsTPCAllTracks", "hSelCrossRowsOverFindClsTPCAllTracks;status;entries", {HistType::kTH1F, {{10, 0., 10.}}});
    registry.add("hSelNClsITSPiFromCharm", "hSelNClsITSPiFromCharm;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelNClsITSInnerPiFromCharm", "hSelNClsITSInnerPiFromCharm;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelMassLam", "hSelMassLam;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelMassCasc", "hSelMassCasc;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelMassCharmBaryon", "hSelMassCharmBaryon;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelDcaXYToPvV0Daughters", "hSelDcaXYToPvV0Daughters;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hSelDcaXYToPvKaFromCasc", "hSelDcaXYToPvKaFromCasc;status;entries", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hMassCascBachMass", "hMassCascBachMass", {HistType::kTH1F, {{3000, 0.0f, 1.0f, "Inv. Mass (GeV/c^{2})"}}});
    registry.add("hBDTScoreTest1", "hBDTScoreTest1", {HistType::kTH1F, {{100, 0.0f, 1.0f, "score"}}});
    registry.add("hBDTScoreTest2", "hBDTScoreTest2", {HistType::kTH1F, {{100, 0.0f, 1.0f, " score"}}});
    registry.add("hBDTScoreTest3", "hBDTScoreTest3", {HistType::kTH1F, {{100, 0.0f, 1.0f, " score"}}});

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

  void process(aod::HfCandOmegaC const& candidates,
               TracksSel const&)
  {

    double massLambdaFromPDG = o2::constants::physics::MassLambda0;
    double massOmegaFromPDG = o2::constants::physics::MassOmegaMinus;

    // looping over charm baryon candidates
    for (const auto& candidate : candidates) {

      int statusOmegaC = 0;
      int statusOmegaCBar = 0;
      bool resultSelections = true; // True if the candidate passes all the selections, False otherwise

      auto trackV0PosDau = candidate.posTrack_as<TracksSel>();             // positive V0 daughter
      auto trackV0NegDau = candidate.negTrack_as<TracksSel>();             // negative V0 daughter
      auto trackKaFromCasc = candidate.bachelor_as<TracksSel>();           // kaon <- cascade
      auto trackPiFromCharm = candidate.piFromCharmBaryon_as<TracksSel>(); // pion <- charm baryon

      auto trackPiFromLam = trackV0NegDau;
      auto trackPrFromLam = trackV0PosDau;

      auto ptCand = candidate.ptCharmBaryon();

      int8_t signDecay = candidate.signDecay(); // sign of Kaon <- cascade

      outputMlOmegaC.clear();    // Ml score
      outputMlOmegaCbar.clear(); // Ml score

      if (signDecay > 0) {
        trackPiFromLam = trackV0PosDau;
        trackPrFromLam = trackV0NegDau;
        registry.fill(HIST("hSelSignDec"), 1); // anti-particle decay
      } else if (signDecay < 0) {
        registry.fill(HIST("hSelSignDec"), 0); // particle decay
      }
      { // FIXME: could be done better...
      }
      // eta selection
      double etaV0PosDau = candidate.etaV0PosDau();
      double etaV0NegDau = candidate.etaV0NegDau();
      double etaKaFromCasc = candidate.etaKaFromCasc();
      double etaPiFromCharmBaryon = candidate.etaPiFromCharmBaryon();
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
      if (candidate.dcaXYToPvV0Dau0() < dcaPosToPvMin || candidate.dcaXYToPvV0Dau1() < dcaNegToPvMin) {
        resultSelections = false;
        registry.fill(HIST("hSelDcaXYToPvV0Daughters"), 0);
      } else {
        registry.fill(HIST("hSelDcaXYToPvV0Daughters"), 1);
      }

      // dcaXY Kaon <-- cascade to PV cut
      if (candidate.dcaXYToPvCascDau() < dcaBachToPvMin) {
        resultSelections = false;
        registry.fill(HIST("hSelDcaXYToPvKaFromCasc"), 0);
      } else {
        registry.fill(HIST("hSelDcaXYToPvKaFromCasc"), 1);
      }

      // cut on charm bachelor pion dcaXY and dcaZ
      if ((std::abs(candidate.impactParPiFromCharmBaryonXY()) < impactParameterXYPiFromCharmBaryonMin) || (std::abs(candidate.impactParPiFromCharmBaryonXY()) > impactParameterXYPiFromCharmBaryonMax)) {
        resultSelections = false;
        registry.fill(HIST("hSelDCAXYPrimPi"), 0);
      } else {
        registry.fill(HIST("hSelDCAXYPrimPi"), 1);
      }
      if ((std::abs(candidate.impactParPiFromCharmBaryonZ()) < impactParameterZPiFromCharmBaryonMin) || (std::abs(candidate.impactParPiFromCharmBaryonZ()) > impactParameterZPiFromCharmBaryonMax)) {
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
      double ptKaFromCasc = RecoDecay::sqrtSumOfSquares(candidate.pxKaFromCasc(), candidate.pyKaFromCasc());
      double ptPiFromCharmBaryon = RecoDecay::sqrtSumOfSquares(candidate.pxPiFromCharmBaryon(), candidate.pyPiFromCharmBaryon());
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

      //  TPC clusters selections
      if (trackPiFromCharm.tpcNClsFound() < nClustersTpcMin) {
        resultSelections = false;
        registry.fill(HIST("hSelNClsTPCPiFromCharm"), 0);
      } else {
        registry.fill(HIST("hSelNClsTPCPiFromCharm"), 1);
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
      if (trackKaFromCasc.tpcNClsFound() < nClustersTpcMin) {
        resultSelections = false;
        registry.fill(HIST("hSelNClsTPCKaFromCasc"), 0);
      } else {
        registry.fill(HIST("hSelNClsTPCKaFromCasc"), 1);
      }

      // TPC crossed rows selection
      if (trackPiFromCharm.tpcNClsCrossedRows() < nTpcCrossedRowsMin) {
        resultSelections = false;
        registry.fill(HIST("hSelNCrossRowsTPCPiFromCharm"), 0);
      } else {
        registry.fill(HIST("hSelNCrossRowsTPCPiFromCharm"), 1);
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
      if (trackKaFromCasc.tpcNClsCrossedRows() < nTpcCrossedRowsMin) {
        resultSelections = false;
        registry.fill(HIST("hSelNCrossRowsTPCKaFromCasc"), 0);
      } else {
        registry.fill(HIST("hSelNCrossRowsTPCKaFromCasc"), 1);
      }

      // further TPC selection
      if (trackPiFromCharm.tpcCrossedRowsOverFindableCls() < tpcCrossedRowsOverFindableClustersRatioMin) {
        resultSelections = false;
        registry.fill(HIST("hSelCrossRowsOverFindClsTPCAllTracks"), 0);
      } else {
        registry.fill(HIST("hSelCrossRowsOverFindClsTPCAllTracks"), 1);
      }
      if (trackKaFromCasc.tpcCrossedRowsOverFindableCls() < tpcCrossedRowsOverFindableClustersRatioMin) {
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
      if (trackPiFromCharm.itsNCls() < nClustersItsMin) {
        resultSelections = false;
        registry.fill(HIST("hSelNClsITSPiFromCharm"), 0);
      } else {
        registry.fill(HIST("hSelNClsITSPiFromCharm"), 1);
      }
      if (trackPiFromCharm.itsNClsInnerBarrel() < nClustersItsInnBarrMin) {
        resultSelections = false;
        registry.fill(HIST("hSelNClsITSInnerPiFromCharm"), 0);
      } else {
        registry.fill(HIST("hSelNClsITSInnerPiFromCharm"), 1);
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
      if (resultSelections) {
        registry.fill(HIST("hStatusCheck"), -0.5);
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
      if (std::abs(invMassLambda - massLambdaFromPDG) < v0MassWindow) {
        statusInvMassLambda = true;
        registry.fill(HIST("hSelMassLam"), 1);
        if (statusPidLambda && statusPidCascade && statusPidCharmBaryon && resultSelections) {
          registry.fill(HIST("hStatusCheck"), 3.5);
        }
      } else {
        registry.fill(HIST("hSelMassLam"), 0);
      }

      if (std::abs(invMassCascade - massOmegaFromPDG) < cascadeMassWindow) {
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
      //############################# ML selections #############################
      if (applyMl) {

        bool isSelectedMlOmegac = false;
        bool isSelectedMlOmegacBar = false;

        if (signDecay > 0) {
          statusOmegaC = 1;
        }

        if (signDecay < 0) {
          statusOmegaCBar = 1;
        }

        if (ptCand <= 1.f || ptCand >= 6.0f) {
          hfMlToOmegaPiCandidate(outputMlOmegaC, outputMlOmegaCbar);
        } else {
          if (statusOmegaC > 0) {
            std::vector<float> inputFeaturesOmegaC = hfMlResponse.getInputFeatures(candidate,
                                                                                   trackV0PosDau, trackV0NegDau, trackKaFromCasc, trackPiFromCharm, o2::constants::physics::kOmegaC0);
            isSelectedMlOmegac = hfMlResponse.isSelectedMl(inputFeaturesOmegaC, ptCand, outputMlOmegaC);
            registry.fill(HIST("hBDTScoreTest1"), outputMlOmegaC[0]);
            registry.fill(HIST("hBDTScoreTest2"), outputMlOmegaC[1]);
            registry.fill(HIST("hBDTScoreTest3"), outputMlOmegaC[2]);
          }
          if (statusOmegaCBar > 0) {
            std::vector<float> inputFeaturesOmegaCBar = hfMlResponse.getInputFeatures(candidate,
                                                                                      trackV0PosDau, trackV0NegDau, trackKaFromCasc, trackPiFromCharm, o2::constants::physics::kOmegaC0);
            isSelectedMlOmegacBar = hfMlResponse.isSelectedMl(inputFeaturesOmegaCBar, ptCand, outputMlOmegaCbar);
          }
          hfMlToOmegaPiCandidate(outputMlOmegaC, outputMlOmegaCbar);
        }
        if (!isSelectedMlOmegac) {
          statusOmegaC = 0;
        }
        if (!isSelectedMlOmegacBar) {
          statusOmegaCBar = 0;
        }
      }
      //#####################End#########################################
      HfSelOmegacToOmegaPi(statusPidLambda, statusPidCascade, statusPidCharmBaryon, statusInvMassLambda,
                           statusInvMassCascade, statusInvMassCharmBaryon, resultSelections, infoTpcStored, infoTofStored,
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
};  // end struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateSelectorOmegacToOmegaPiWithKfp>(cfgc)};
}
