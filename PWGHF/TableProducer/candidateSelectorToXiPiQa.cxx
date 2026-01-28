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

/// \file candidateSelectorToXiPiQa.cxx
/// \brief Selection of Xic0 and Xicp candidates
///
/// \author Jinhyun Park <jinhyun.park@cern.ch>, Pusan National University
/// \author Krista Smith <krista.lizbeth.smith@cern.ch>, Pusan National University

#include "PWGHF/Core/HfMlResponseXic0ToXiPi.h"
#include "PWGHF/Core/HfMlResponseXic0ToXiPiKf.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsAnalysis.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelectorPID.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
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

/// Struct for applying Omegac0/Xic0 selection cuts
struct HfCandidateSelectorToXiPiQa {

  // DCAFitter
  Produces<aod::HfSelToXiPi> hfSelToXiPi;
  // KFParticle
  Produces<aod::HfSelToXiPiKf> hfSelToXiPiKf;
  // ML selection - Filled with both DCAFitter and KFParticle
  Produces<aod::HfMlToXiPi> hfMlToXiPi;

  // cuts from SelectorCuts.h - pT dependent cuts
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_xic_to_xi_pi::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_xic_to_xi_pi::Cuts[0], hf_cuts_xic_to_xi_pi::NBinsPt, hf_cuts_xic_to_xi_pi::NCutVars, hf_cuts_xic_to_xi_pi::labelsPt, hf_cuts_xic_to_xi_pi::labelsCutVar}, "Xic0 candidate selection per Pt Bin"};
  // ML inference
  Configurable<bool> applyMl{"applyMl", false, "Flag to apply ML selections"};
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
  // LF analysis selections
  Configurable<double> etaTrackCharmBachMax{"etaTrackCharmBachMax", 0.8, "Max absolute value of eta for charm baryon bachelor"};
  Configurable<double> etaTrackLFDauMax{"etaTrackLFDauMax", 0.8, "Max absolute value of eta for V0 and cascade daughters"};
  Configurable<double> ptPiFromCascMin{"ptPiFromCascMin", 0.15, "Min pT pi <--Casc"};
  Configurable<double> ptPiFromCharmBaryonMin{"ptPiFromCharmBaryonMin", 0.2, "Min pT pi <--Casc"};
  Configurable<double> radiusCascMin{"radiusCascMin", 0.5, "Min Cascade radius"};
  Configurable<double> radiusV0Min{"radiusV0Min", 1.1, "Min V0 radius"};
  Configurable<double> impactParXYPiFromCharmBaryonMin{"impactParXYPiFromCharmBaryonMin", 0., "Min dcaxy pi from charm baryon track to pV"};
  Configurable<double> impactParXYPiFromCharmBaryonMax{"impactParXYPiFromCharmBaryonMax", 10., "Max dcaxy pi from charm baryon track to pV"};
  Configurable<double> impactParXYCascMin{"impactParXYCascMin", 0., "Min dcaxy casc track to pV"};
  Configurable<double> impactParXYCascMax{"impactParXYCascMax", 10., "Max dcaxy casc track to pV"};
  Configurable<double> impactParZPiFromCharmBaryonMin{"impactParZPiFromCharmBaryonMin", 0., "Min dcaz pi from charm baryon track to pV"};
  Configurable<double> impactParZPiFromCharmBaryonMax{"impactParZPiFromCharmBaryonMax", 10., "Max dcaz pi from charm baryon track to pV"};
  Configurable<double> impactParZCascMin{"impactParZCascMin", 0., "Min dcaz casc track to pV"};
  Configurable<double> impactParZCascMax{"impactParZCascMax", 10., "Max dcaz casc track to pV"};
  Configurable<bool> applyTrkSelLf{"applyTrkSelLf", true, "Apply track selection for LF daughters"};
  // Mass window
  Configurable<double> v0MassWindow{"v0MassWindow", 0.01, "v0 mass window"};
  Configurable<double> cascMassWindow{"cascMassWindow", 0.01, "cascade mass window"};
  Configurable<double> invMassCharmBaryonMin{"invMassCharmBaryonMin", 2.0, "Lower limit of invariant mass spectrum charm baryon"};
  Configurable<double> invMassCharmBaryonMax{"invMassCharmBaryonMax", 3.1, "Lower limit of invariant mass spectrum charm baryon"};
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
  // PID - TOF selections
  Configurable<double> ptPiPidTofMin{"ptPiPidTofMin", -1, "Lower bound of track pT for TOF PID for pion selection"};
  Configurable<double> ptPiPidTofMax{"ptPiPidTofMax", 9999.9, "Upper bound of track pT for TOF PID for pion selection"};
  Configurable<double> nSigmaTofPiMax{"nSigmaTofPiMax", 3., "Nsigma cut on TOF only for pion selection"};
  Configurable<double> nSigmaTofCombinedPiMax{"nSigmaTofCombinedPiMax", 0., "Nsigma cut on TOF combined with TPC for pion selection"};
  Configurable<double> ptPrPidTofMin{"ptPrPidTofMin", -1, "Lower bound of track pT for TOF PID for proton selection"};
  Configurable<double> ptPrPidTofMax{"ptPrPidTofMax", 9999.9, "Upper bound of track pT for TOF PID for proton selection"};
  Configurable<double> nSigmaTofPrMax{"nSigmaTofPrMax", 3., "Nsigma cut on TOF only for proton selection"};
  Configurable<double> nSigmaTofCombinedPrMax{"nSigmaTofCombinedPrMax", 0., "Nsigma cut on TOF combined with TPC for proton selection"};
  // detector track quality selections
  Configurable<int> nClustersTpcMin{"nClustersTpcMin", 70, "Minimum number of TPC clusters requirement"};
  Configurable<int> nTpcCrossedRowsMin{"nTpcCrossedRowsMin", 70, "Minimum number of TPC crossed rows requirement"};
  Configurable<double> tpcCrossedRowsOverFindableClustersRatioMin{"tpcCrossedRowsOverFindableClustersRatioMin", 0.8, "Minimum ratio TPC crossed rows over findable clusters requirement"};
  Configurable<float> tpcChi2PerClusterMax{"tpcChi2PerClusterMax", 4, "Maximum value of chi2 fit over TPC clusters"};
  Configurable<int> nClustersItsMin{"nClustersItsMin", 3, "Minimum number of ITS clusters requirement for pi <- charm baryon"};
  Configurable<int> nClustersItsInnBarrMin{"nClustersItsInnBarrMin", 1, "Minimum number of ITS clusters in inner barrel requirement for pi <- charm baryon"};
  Configurable<float> itsChi2PerClusterMax{"itsChi2PerClusterMax", 36, "Maximum value of chi2 fit over ITS clusters for pi <- charm baryon"};

  o2::analysis::HfMlResponseXic0ToXiPi<float> hfMlResponseDca;
  o2::analysis::HfMlResponseXic0ToXiPiKf<float> hfMlResponseKf;
  std::vector<float> outputMlXic0ToXiPi = {};
  o2::ccdb::CcdbApi ccdbApi;

  TrackSelectorPi selectorPion;
  TrackSelectorPr selectorProton;

  using TracksSel = soa::Join<aod::TracksWDcaExtra, aod::TracksPidPi, aod::TracksPidPr>;

  HistogramRegistry registry{"registry"}; // for QA of selections

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

    const AxisSpec axisSel{2, -0.5, 1.5, "status"};

    registry.add("hSelPID", "hSelPID;status;entries", {HistType::kTH1F, {{12, 0., 12.}}});
    registry.add("hStatusCheck", "Check consecutive selections status;status;entries", {HistType::kTH1F, {{12, 0., 12.}}});
    registry.add("hSelSignDec", "hSelSignDec;status;entries", {HistType::kTH1F, {axisSel}});

    // for QA of the selections (bin 0 -> candidates that did not pass the selection, bin 1 -> candidates that passed the selection)
    registry.add("hSelEtaPosV0Dau", "hSelEtaPosV0Dau;status;entries", {HistType::kTH1F, {axisSel}});
    registry.add("hSelEtaNegV0Dau", "hSelEtaNegV0Dau;status;entries", {HistType::kTH1F, {axisSel}});
    registry.add("hSelEtaPiFromCasc", "hSelEtaPiFromCasc;status;entries", {HistType::kTH1F, {axisSel}});
    registry.add("hSelEtaPiFromCharm", "hSelEtaPiFromCharm;status;entries", {HistType::kTH1F, {axisSel}});
    registry.add("hSelRadCasc", "hSelRadCasc;status;entries", {HistType::kTH1F, {axisSel}});
    registry.add("hSelRadV0", "hSelRadV0;status;entries", {HistType::kTH1F, {axisSel}});
    registry.add("hSelCosPACasc", "hSelCosPACasc;status;entries", {HistType::kTH1F, {axisSel}});
    registry.add("hSelCosPAV0", "hSelCosPAV0;status;entries", {HistType::kTH1F, {axisSel}});
    registry.add("hSelDCACascDau", "hSelDCACascDau;status;entries", {HistType::kTH1F, {axisSel}});
    registry.add("hSelDCAV0Dau", "hSelDCAV0Dau;status;entries", {HistType::kTH1F, {axisSel}});
    registry.add("hSelDCACharmDau", "hSelDCACharmDau;status;entries", {HistType::kTH1F, {axisSel}});
    registry.add("hSelDCAXYPrimPi", "hSelDCAXYPrimPi;status;entries", {HistType::kTH1F, {axisSel}});
    registry.add("hSelDCAZPrimPi", "hSelDCAZPrimPi;status;entries", {HistType::kTH1F, {axisSel}});
    registry.add("hSelDCAXYCasc", "hSelDCAXYCasc;status;entries", {HistType::kTH1F, {axisSel}});
    registry.add("hSelDCAZCasc", "hSelDCAZCasc;status;entries", {HistType::kTH1F, {axisSel}});
    registry.add("hSelPtPiFromCasc", "hSelPtPiFromCasc;status;entries", {HistType::kTH1F, {axisSel}});
    registry.add("hSelPtPiFromCharm", "hSelPtPiFromCharm;status;entries", {HistType::kTH1F, {axisSel}});
    registry.add("hSelTPCQualityPiFromCharm", "hSelTPCQualityPiFromCharm;status;entries", {HistType::kTH1F, {axisSel}});
    registry.add("hSelTPCQualityPiFromLam", "hSelTPCQualityPiFromLam;status;entries", {HistType::kTH1F, {axisSel}});
    registry.add("hSelTPCQualityPrFromLam", "hSelTPCQualityPrFromLam;status;entries", {HistType::kTH1F, {axisSel}});
    registry.add("hSelTPCQualityPiFromCasc", "hSelTPCQualityPiFromCasc;status;entries", {HistType::kTH1F, {axisSel}});
    registry.add("hSelITSQualityPiFromCharm", "hSelITSQualityPiFromCharm;status;entries", {HistType::kTH1F, {axisSel}});
    registry.add("hSelMassLam", "hSelMassLam;status;entries", {HistType::kTH1F, {axisSel}});
    registry.add("hSelMassCasc", "hSelMassCasc;status;entries", {HistType::kTH1F, {axisSel}});
    registry.add("hSelMassCharmBaryon", "hSelMassCharmBaryon;status;entries", {HistType::kTH1F, {axisSel}});
    registry.add("hSelDcaXYToPvV0Daughters", "hSelDcaXYToPvV0Daughters;status;entries", {HistType::kTH1F, {axisSel}});
    registry.add("hSelDcaXYToPvPiFromCasc", "hSelDcaXYToPvPiFromCasc;status;entries", {HistType::kTH1F, {axisSel}});

    // temporary check histogram for KFParticle
    const AxisSpec axisKf{20, -0.5, 19.5, "Pass status"};
    registry.add("hStatusCheckKf", "hStatusCheck;kf status;entries", {HistType::kTH1F, {axisKf}});
    registry.add("hStatusCheckKf_res", "hStatusCheck;kf status;entries", {HistType::kTH1F, {axisSel}});

    // invarinat mass histograms
    registry.add("hInvMassCharmBaryon", "Charm baryon invariant mass; int mass; entries", {HistType::kTH1F, {{1500, 1.5, 4.5}}});
    registry.add("hInvMassCharmBaryonBkg", "Charm baryon invariant mass, rejected; int mass; entries", {HistType::kTH1F, {{1500, 1.5, 4.5}}});

    // HfMlResponse initialization
    if (applyMl) {
      if (doprocessSelectionDCAFitter) {
        hfMlResponseDca.configure(binsPtMl, cutsMl, cutDirMl, nClassesMl);
        if (loadModelsFromCCDB) {
          ccdbApi.init(ccdbUrl);
          hfMlResponseDca.setModelPathsCCDB(onnxFileNames, ccdbApi, modelPathsCCDB, timestampCCDB);
        } else {
          hfMlResponseDca.setModelPathsLocal(onnxFileNames);
        }
        hfMlResponseDca.cacheInputFeaturesIndices(namesInputFeatures);
        hfMlResponseDca.init();
      } else {
        hfMlResponseKf.configure(binsPtMl, cutsMl, cutDirMl, nClassesMl);
        if (loadModelsFromCCDB) {
          ccdbApi.init(ccdbUrl);
          hfMlResponseKf.setModelPathsCCDB(onnxFileNames, ccdbApi, modelPathsCCDB, timestampCCDB);
        } else {
          hfMlResponseKf.setModelPathsLocal(onnxFileNames);
        }
        hfMlResponseKf.cacheInputFeaturesIndices(namesInputFeatures);
        hfMlResponseKf.init();
      }
    }
  }

  // KF Specific cuts
  template <typename T>
  bool ExtraSelectionKf(const T& candidate)
  {

    bool selectionResult = true;

    auto candPt = RecoDecay::pt(candidate.pxCharmBaryon(), candidate.pyCharmBaryon());
    const int pTBin = findBin(binsPt, candPt);

    if (pTBin == -1) {
      selectionResult = false;
    } else {
      registry.fill(HIST("hStatusCheckKf"), 0.0);
    }

    // cosPa (KF Specific)
    if (candidate.cosPaCascToXic() < cuts->get(pTBin, "cosPaCascToXic")) {
      selectionResult = false;
    } else {
      registry.fill(HIST("hStatusCheckKf"), 1.0);
    }
    if (candidate.cosPaV0ToCasc() < cuts->get(pTBin, "cosPaV0ToCasc")) {
      selectionResult = false;
    } else {
      registry.fill(HIST("hStatusCheckKf"), 2.0);
    }

    // dcaXY pion0 <-- Xic0 to PV cut
    if (std::abs(candidate.kfDcaXYPiFromXic()) > cuts->get(pTBin, "kfDcaXYPiFromXic")) {
      selectionResult = false;
    } else {
      registry.fill(HIST("hStatusCheckKf"), 3.0);
    }

    // dcaXY cascade tp PV cut
    if (std::abs(candidate.kfDcaXYCascToPv()) > cuts->get(pTBin, "kfDcaXYCascToPv")) {
      selectionResult = false;
    } else {
      registry.fill(HIST("hStatusCheckKf"), 4.0);
    }

    // Chi2Ge0
    if (candidate.chi2GeoXic() < 0 || candidate.chi2GeoXic() > cuts->get(pTBin, "chi2GeoXic")) {
      selectionResult = false;
    } else {
      registry.fill(HIST("hStatusCheckKf"), 5.0);
    }
    if (candidate.chi2GeoCasc() < 0 || candidate.chi2GeoCasc() > cuts->get(pTBin, "chi2GeoCasc")) {
      selectionResult = false;
    } else {
      registry.fill(HIST("hStatusCheckKf"), 6.0);
    }
    if (candidate.chi2GeoV0() < 0 || candidate.chi2GeoV0() > cuts->get(pTBin, "chi2GeoV0")) {
      selectionResult = false;
    } else {
      registry.fill(HIST("hStatusCheckKf"), 7.0);
    }

    // Chi2 Topo
    if (candidate.chi2TopoXicToPv() < 0 || candidate.chi2TopoXicToPv() > cuts->get(pTBin, "chi2TopoXicToPv")) {
      selectionResult = false;
    } else {
      registry.fill(HIST("hStatusCheckKf"), 8.0);
    }
    if (candidate.chi2TopoPiFromXicToPv() < 0 || candidate.chi2TopoPiFromXicToPv() > cuts->get(pTBin, "chi2TopoPiFromXicToPv")) {
      selectionResult = false;
    } else {
      registry.fill(HIST("hStatusCheckKf"), 9.0);
    }
    if (candidate.chi2TopoCascToPv() < 0 || candidate.chi2TopoCascToPv() > cuts->get(pTBin, "chi2TopoCascToPv")) {
      selectionResult = false;
    } else {
      registry.fill(HIST("hStatusCheckKf"), 10.0);
    }
    if (candidate.chi2TopoV0ToPv() < 0 || candidate.chi2TopoV0ToPv() < cuts->get(pTBin, "chi2TopoV0ToPv")) {
      selectionResult = false;
    } else {
      registry.fill(HIST("hStatusCheckKf"), 11.0);
    }
    if (candidate.chi2TopoV0ToCasc() < 0 || candidate.chi2TopoV0ToCasc() > cuts->get(pTBin, "chi2TopoV0ToCasc")) {
      selectionResult = false;
    } else {
      registry.fill(HIST("hStatusCheckKf"), 12.0);
    }
    if (candidate.chi2TopoCascToXic() < 0 || candidate.chi2TopoCascToXic() > cuts->get(pTBin, "chi2TopoCascToXic")) {
      selectionResult = false;
    } else {
      registry.fill(HIST("hStatusCheckKf"), 13.0);
    }

    // ldl
    if (candidate.cascldl() < cuts->get(pTBin, "cascldl")) {
      selectionResult = false;
    } else {
      registry.fill(HIST("hStatusCheckKf"), 14.0);
    }
    if (candidate.v0ldl() < cuts->get(pTBin, "v0ldl")) {
      selectionResult = false;
    } else {
      registry.fill(HIST("hStatusCheckKf"), 15.0);
    }

    // decay length
    if (std::abs(candidate.decayLenXYXic()) > cuts->get(pTBin, "decayLenXYXic")) {
      selectionResult = false;
    } else {
      registry.fill(HIST("hStatusCheckKf"), 16.0);
    }
    if (std::abs(candidate.decayLenXYCasc()) < cuts->get(pTBin, "decayLenXYCasc")) {
      selectionResult = false;
    } else {
      registry.fill(HIST("hStatusCheckKf"), 17.0);
    }
    if (std::abs(candidate.decayLenXYLambda()) < cuts->get(pTBin, "decayLenXYLambda")) {
      selectionResult = false;
    } else {
      registry.fill(HIST("hStatusCheckKf"), 18.0);
    }

    // ctau
    if (std::abs(candidate.cTauXic()) > cuts->get(pTBin, "cTauXic")) {
      selectionResult = false;
    } else {
      registry.fill(HIST("hStatusCheckKf"), 19.0);
    }

    // Fill the kf output result
    if (selectionResult) {
      registry.fill(HIST("hStatusCheckKf_res"), 1);
    } else {
      registry.fill(HIST("hStatusCheckKf_res"), 0);
    }

    return selectionResult;
  }

  // DCA specific cuts
  template <typename T>
  bool ExtraSelectionDca(const T& candidate)
  {

    bool selectionResult = true;

    auto candPt = RecoDecay::pt(candidate.pxCharmBaryon(), candidate.pyCharmBaryon());
    const int pTBin = findBin(binsPt, candPt);

    if (pTBin == -1) {
      selectionResult = false;
    }

    // impact parmeter
    if (std::abs(candidate.impactParBachFromCharmBaryonXY()) < impactParXYPiFromCharmBaryonMin || std::abs(candidate.impactParBachFromCharmBaryonXY()) > impactParXYPiFromCharmBaryonMax) {
      selectionResult = false;
    }

    if (std::abs(candidate.impactParBachFromCharmBaryonZ()) < impactParZPiFromCharmBaryonMin || std::abs(candidate.impactParBachFromCharmBaryonZ()) > impactParZPiFromCharmBaryonMax) {
      selectionResult = false;
    }

    if (std::abs(candidate.impactParCascXY()) < impactParXYCascMin || std::abs(candidate.impactParCascXY()) > impactParXYCascMax) {
      selectionResult = false;
    }

    if (std::abs(candidate.impactParCascZ()) < impactParZCascMin || std::abs(candidate.impactParCascZ()) > impactParZCascMax) {
      selectionResult = false;
    }

    return selectionResult;
  }

  template <bool dokf, typename TCandTable>
  void runSelection(TCandTable const& candidates, TracksSel const&)
  {
    double massLambdaFromPDG = o2::constants::physics::MassLambda0;
    double massXiFromPDG = o2::constants::physics::MassXiMinus;

    // looping over charm baryon candidates
    for (const auto& candidate : candidates) {

      bool resultSelections = true; // True if the candidate passes all the selections, False otherwise
      outputMlXic0ToXiPi.clear();

      auto trackV0PosDau = candidate.template posTrack_as<TracksSel>();
      auto trackV0NegDau = candidate.template negTrack_as<TracksSel>();
      auto trackPiFromCasc = candidate.template bachelor_as<TracksSel>();
      auto trackPiFromCharm = candidate.template bachelorFromCharmBaryon_as<TracksSel>();

      int8_t signDecay = candidate.signDecay(); // sign of pi <- cascade

      auto trackPiFromLam = trackV0NegDau;
      auto trackPrFromLam = trackV0PosDau;

      if (signDecay > 0) {
        trackPiFromLam = trackV0PosDau;
        trackPrFromLam = trackV0NegDau;
        registry.fill(HIST("hSelSignDec"), 1); // anti-particle decay
      } else {
        registry.fill(HIST("hSelSignDec"), 0); // particle decay
      }

      // pT selection
      auto ptCandXic0 = RecoDecay::pt(candidate.pxCharmBaryon(), candidate.pyCharmBaryon());
      int pTBin = findBin(binsPt, ptCandXic0);
      if (pTBin == -1) {
        resultSelections = false;
      }

      // eta selection
      double etaV0PosDau = candidate.etaV0PosDau();
      double etaV0NegDau = candidate.etaV0NegDau();
      double etaPiFromCasc = candidate.etaBachFromCasc();
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
      if (std::abs(etaPiFromCasc) > etaTrackLFDauMax) {
        resultSelections = false;
        registry.fill(HIST("hSelEtaPiFromCasc"), 0);
      } else {
        registry.fill(HIST("hSelEtaPiFromCasc"), 1);
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
      if (candidate.cosPACasc() < cuts->get(pTBin, "cosPACasc")) {
        resultSelections = false;
        registry.fill(HIST("hSelCosPACasc"), 0);
      } else {
        registry.fill(HIST("hSelCosPACasc"), 1);
      }
      if (candidate.cosPAV0() < cuts->get(pTBin, "cosPAV0")) {
        resultSelections = false;
        registry.fill(HIST("hSelCosPAV0"), 0);
      } else {
        registry.fill(HIST("hSelCosPAV0"), 1);
      }

      // cascade and v0 daughters dca cut (LF cut)
      if (candidate.dcaCascDau() > cuts->get(pTBin, "dcaCascDau")) {
        resultSelections = false;
        registry.fill(HIST("hSelDCACascDau"), 0);
      } else {
        registry.fill(HIST("hSelDCACascDau"), 1);
      }

      if (candidate.dcaV0Dau() > cuts->get(pTBin, "dcaV0Dau")) {
        resultSelections = false;
        registry.fill(HIST("hSelDCAV0Dau"), 0);
      } else {
        registry.fill(HIST("hSelDCAV0Dau"), 1);
      }

      // dca charm baryon daughters cut
      if (candidate.dcaCharmBaryonDau() > cuts->get(pTBin, "dcaCharmBaryonDau")) {
        resultSelections = false;
        registry.fill(HIST("hSelDCACharmDau"), 0);
      } else {
        registry.fill(HIST("hSelDCACharmDau"), 1);
      }

      // dcaXY v0 daughters to PV cut
      if (std::abs(candidate.dcaXYToPvV0Dau0()) < cuts->get(pTBin, "dcaXYToPvV0Dau0") || std::abs(candidate.dcaXYToPvV0Dau1()) < cuts->get(pTBin, "dcaXYToPvV0Dau1")) {
        resultSelections = false;
        registry.fill(HIST("hSelDcaXYToPvV0Daughters"), 0);
      } else {
        registry.fill(HIST("hSelDcaXYToPvV0Daughters"), 1);
      }

      // dcaXY pi <-- cascade to PV cut
      if (std::abs(candidate.dcaXYToPvCascDau()) < cuts->get(pTBin, "dcaXYToPvCascDau")) {
        resultSelections = false;
        registry.fill(HIST("hSelDcaXYToPvPiFromCasc"), 0);
      } else {
        registry.fill(HIST("hSelDcaXYToPvPiFromCasc"), 1);
      }

      // pT selections
      double ptPiFromCasc = RecoDecay::sqrtSumOfSquares(candidate.pxBachFromCasc(), candidate.pyBachFromCasc());
      double ptPiFromCharmBaryon = RecoDecay::sqrtSumOfSquares(candidate.pxBachFromCharmBaryon(), candidate.pyBachFromCharmBaryon());

      if (std::abs(ptPiFromCasc) < ptPiFromCascMin) {
        resultSelections = false;
        registry.fill(HIST("hSelPtPiFromCasc"), 0);
      } else {
        registry.fill(HIST("hSelPtPiFromCasc"), 1);
      }
      if (std::abs(ptPiFromCharmBaryon) <= cuts->get(pTBin, "ptPiFromCharmBaryon")) {
        resultSelections = false;
        registry.fill(HIST("hSelPtPiFromCharm"), 0);
      } else {
        registry.fill(HIST("hSelPtPiFromCharm"), 1);
      }

      // DCAFitter && KFParticle specific selection
      if constexpr (dokf == false) {
        resultSelections = ExtraSelectionDca(candidate);
      } else {
        resultSelections = ExtraSelectionKf(candidate);
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
        SETBIT(infoTpcStored, PrFromLam);
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
        SETBIT(infoTofStored, PrFromLam);
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

      if (std::abs(invMassLambda - massLambdaFromPDG) < v0MassWindow) {
        statusInvMassLambda = true;
        registry.fill(HIST("hSelMassLam"), 1);
        if (statusPidLambda && statusPidCascade && statusPidCharmBaryon && resultSelections) {
          registry.fill(HIST("hStatusCheck"), 3.5);
        }
      } else {
        resultSelections = false;
        registry.fill(HIST("hSelMassLam"), 0);
      }

      if (std::abs(invMassCascade - massXiFromPDG) < cascMassWindow) {
        statusInvMassCascade = true;
        registry.fill(HIST("hSelMassCasc"), 1);
        if (statusPidLambda && statusPidCascade && statusPidCharmBaryon && statusInvMassLambda && resultSelections) {
          registry.fill(HIST("hStatusCheck"), 4.5);
        }
      } else {
        resultSelections = false;
        registry.fill(HIST("hSelMassCasc"), 0);
      }

      if ((invMassCharmBaryon >= invMassCharmBaryonMin) && (invMassCharmBaryon <= invMassCharmBaryonMax)) {
        statusInvMassCharmBaryon = true;
        registry.fill(HIST("hSelMassCharmBaryon"), 1);
        if (statusPidLambda && statusPidCascade && statusPidCharmBaryon && statusInvMassLambda && statusInvMassCascade && resultSelections) {
          registry.fill(HIST("hStatusCheck"), 5.5);
        }
      } else {
        resultSelections = false;
        registry.fill(HIST("hSelMassCharmBaryon"), 0);
      }

      // ML BDT selection
      if (applyMl) {
        bool isSelectedMlXic0 = false;
        std::vector<float> inputFeaturesXic0 = {};

        if constexpr (dokf == false) {
          inputFeaturesXic0 = hfMlResponseDca.getInputFeatures(candidate, trackPiFromLam, trackPiFromCasc, trackPiFromCharm);
          isSelectedMlXic0 = hfMlResponseDca.isSelectedMl(inputFeaturesXic0, ptCandXic0, outputMlXic0ToXiPi);
        } else {
          inputFeaturesXic0 = hfMlResponseKf.getInputFeatures(candidate, trackPiFromLam, trackPiFromCasc, trackPiFromCharm);
          isSelectedMlXic0 = hfMlResponseKf.isSelectedMl(inputFeaturesXic0, ptCandXic0, outputMlXic0ToXiPi);
        }

        if (!isSelectedMlXic0) {
          continue;
        }

        hfMlToXiPi(outputMlXic0ToXiPi);
      }

      // Fill in selection result
      if constexpr (dokf == false) {
        hfSelToXiPi(statusPidLambda, statusPidCascade, statusPidCharmBaryon, statusInvMassLambda, statusInvMassCascade, statusInvMassCharmBaryon, resultSelections, infoTpcStored, infoTofStored,
                    trackPiFromCharm.tpcNSigmaPi(), trackPiFromCasc.tpcNSigmaPi(), trackPiFromLam.tpcNSigmaPi(), trackPrFromLam.tpcNSigmaPr(),
                    trackPiFromCharm.tofNSigmaPi(), trackPiFromCasc.tofNSigmaPi(), trackPiFromLam.tofNSigmaPi(), trackPrFromLam.tofNSigmaPr());
      } else {
        // Set result of selection to false if one of the statusPID or statusInvMass is false
        // This is required because selection table for KF does not store any information about statusPID or statusInvMass while DCAFitter does.
        if (!(statusPidLambda && statusPidCascade && statusPidCharmBaryon && statusInvMassCharmBaryon && statusInvMassCascade && statusInvMassLambda)) {
          resultSelections = false;
        }

        hfSelToXiPiKf(resultSelections,
                      // statusPidCharmBaryon, statusPidCascade, statusPidLambda, statusInvMassCharmBaryon, statusInvMassCascade, statusInvMassLambda, infoTpcStored, infoTofStored,
                      trackPiFromCharm.tpcNSigmaPi(), trackPiFromCasc.tpcNSigmaPi(), trackPiFromLam.tpcNSigmaPi(), trackPrFromLam.tpcNSigmaPr(),
                      trackPiFromCharm.tofNSigmaPi(), trackPiFromCasc.tofNSigmaPi(), trackPiFromLam.tofNSigmaPi(), trackPrFromLam.tofNSigmaPr());
      }

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

      // Fill in invariant mass histogram
      if (statusPidLambda && statusPidCascade && statusPidCharmBaryon && statusInvMassLambda && statusInvMassCascade && statusInvMassCharmBaryon && resultSelections) {
        registry.fill(HIST("hInvMassCharmBaryon"), invMassCharmBaryon);
      } else {
        registry.fill(HIST("hInvMassCharmBaryonBkg"), invMassCharmBaryon);
      }

    } // end of candidate loop
  } // end run fuction

  ///////////////////////////////////
  ///    Process with DCAFitter    //
  ///////////////////////////////////
  void processSelectionDCAFitter(aod::HfCandToXiPi const& candidates, TracksSel const& tracks)
  {
    runSelection<false>(candidates, tracks);
  }
  PROCESS_SWITCH(HfCandidateSelectorToXiPiQa, processSelectionDCAFitter, "Xic0 candidate selection with DCAFitter output", true);

  ////////////////////////////////////
  ///    Process with KFParticle    //
  ////////////////////////////////////
  void processSelectionKFParticle(aod::HfCandToXiPiKf const& candidates, TracksSel const& tracks)
  {
    runSelection<true>(candidates, tracks);
  }
  PROCESS_SWITCH(HfCandidateSelectorToXiPiQa, processSelectionKFParticle, "Xic0 candidate selection with KFParticle output", false);

}; // end struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateSelectorToXiPiQa>(cfgc)};
}
