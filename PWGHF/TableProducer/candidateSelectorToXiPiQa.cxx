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

enum {
  doDcaFitter = 0,
  doKfParticle
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
  Configurable<double> etaTrackLFDauMax{"etaTrackLFDauMax", 1.0, "Max absolute value of eta for V0 and cascade daughters"};
  Configurable<double> ptPiFromCascMin{"ptPiFromCascMin", 0.15, "Min pT pi <--Casc"};
  Configurable<double> ptPiFromCharmBaryonMin{"ptPiFromCharmBaryonMin", 0.2, "Min pT pi <--Casc"};
  Configurable<double> radiusCascMin{"radiusCascMin", 0.6, "Min Cascade radius"};
  Configurable<double> radiusV0Min{"radiusV0Min", 1.2, "Min V0 radius"};
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
  using TracksSelLf = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksPidPi, aod::TracksPidPr>;

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
    const AxisSpec axisSelOnLfDca{14, -0.5, 13.5, "status"};
    const AxisSpec axisSelOnLfKf{23, -0.5, 22.5, "status"};
    const AxisSpec axisSelOnHfDca{6, -0.5, 5.5, "status"};
    const AxisSpec axisSelOnHfKf{11, -0.5, 10.5, "status"};

    registry.add("hSelSignDec", "hSelSignDec;status;entries", {HistType::kTH1F, {axisSel}});
    registry.add("hSelStatusCluster", "hSelStatusCluster:# of events Passed;;", {HistType::kTH1F, {{6, -0.5, 5.5}}});
    registry.get<TH1>(HIST("hSelStatusCluster"))->GetXaxis()->SetBinLabel(1, "All");
    registry.get<TH1>(HIST("hSelStatusCluster"))->GetXaxis()->SetBinLabel(2, "TpcCluster PiFromV0");
    registry.get<TH1>(HIST("hSelStatusCluster"))->GetXaxis()->SetBinLabel(3, "TpcCluster PrFromV0");
    registry.get<TH1>(HIST("hSelStatusCluster"))->GetXaxis()->SetBinLabel(4, "TpcCluster PiFromCasc");
    registry.get<TH1>(HIST("hSelStatusCluster"))->GetXaxis()->SetBinLabel(5, "TpcCluster PiFromCharm");
    registry.get<TH1>(HIST("hSelStatusCluster"))->GetXaxis()->SetBinLabel(6, "ItsCluster PiFromCharm");

    registry.add("hSelStatusPID", "hSelStatusPID;# of events Passed;;", {HistType::kTH1F, {{4, -0.5, 3.5}}});
    registry.get<TH1>(HIST("hSelStatusPID"))->GetXaxis()->SetBinLabel(1, "All");
    registry.get<TH1>(HIST("hSelStatusPID"))->GetXaxis()->SetBinLabel(2, "Lambda");
    registry.get<TH1>(HIST("hSelStatusPID"))->GetXaxis()->SetBinLabel(3, "Cascade");
    registry.get<TH1>(HIST("hSelStatusPID"))->GetXaxis()->SetBinLabel(4, "CharmBaryon");

    // For QA of LF & HF selection
    if (doprocessSelectionDCAFitter) {
      registry.add("hSelStatusLf", "hSelStatusLf;# of candidate passed;", {HistType::kTH1F, {axisSelOnLfDca}});
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(1, "All");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(2, "etaV0Dau");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(3, "radiusV0");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(4, "radiusCasc");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(5, "cosPAV0");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(6, "cosPACasc");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(7, "dcaV0Dau");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(8, "dcaCascDau");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(9, "dcaXYToPvV0Dau0");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(10, "dcaXYToPvV0Dau1");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(11, "dcaXYToPvCascDau");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(12, "ptPiFromCasc");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(13, "impactParCascXY");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(14, "impactParCascZ");

      registry.add("hSelStatusHf", "hSelStatusHf;# of candidate passed;", {HistType::kTH1F, {axisSelOnHfDca}});
      registry.get<TH1>(HIST("hSelStatusHf"))->GetXaxis()->SetBinLabel(1, "All");
      registry.get<TH1>(HIST("hSelStatusHf"))->GetXaxis()->SetBinLabel(2, "etaTrackCharmBach");
      registry.get<TH1>(HIST("hSelStatusHf"))->GetXaxis()->SetBinLabel(3, "dcaCharmBaryonDau");
      registry.get<TH1>(HIST("hSelStatusHf"))->GetXaxis()->SetBinLabel(4, "ptPiFromCharmBaryon");
      registry.get<TH1>(HIST("hSelStatusHf"))->GetXaxis()->SetBinLabel(5, "impactParBachFromCharmXY");
      registry.get<TH1>(HIST("hSelStatusHf"))->GetXaxis()->SetBinLabel(6, "impactParBachFromCharmZ");
    }

    if (doprocessSelectionKFParticle) {
      registry.add("hSelStatusLf", "hSelStatusLf;# of candidate passed;", {HistType::kTH1F, {axisSelOnLfKf}});
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(1, "All");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(2, "etaV0Dau");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(3, "radiusV0");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(4, "radiusCasc");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(5, "cosPAV0");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(6, "cosPACasc");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(7, "dcaV0Dau");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(8, "dcaCascDau");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(9, "dcaXYToPvV0Dau0");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(10, "dcaXYToPvV0Dau1");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(11, "dcaXYToPvCascDau");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(12, "ptPiFromCasc");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(13, "cosPAV0ToCasc");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(14, "kfDcaXYCascToPv");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(15, "chi2GeoV0");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(16, "chi2GeoCasc");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(17, "chi2TopoV0ToPv");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(18, "chi2TopoCascToPv");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(19, "chi2TopoV0ToCasc");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(20, "v0ldl");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(21, "cascldl");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(22, "decayLenXYLambda");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(23, "decayLenXYCasc");

      registry.add("hSelStatusHf", "hSelStatusHf;# of candidate passed;", {HistType::kTH1F, {axisSelOnHfKf}});
      registry.get<TH1>(HIST("hSelStatusHf"))->GetXaxis()->SetBinLabel(1, "All");
      registry.get<TH1>(HIST("hSelStatusHf"))->GetXaxis()->SetBinLabel(2, "etaTrackCharmBach");
      registry.get<TH1>(HIST("hSelStatusHf"))->GetXaxis()->SetBinLabel(3, "dcaCharmBaryonDau");
      registry.get<TH1>(HIST("hSelStatusHf"))->GetXaxis()->SetBinLabel(4, "ptPiFromCharmBaryon");
      registry.get<TH1>(HIST("hSelStatusHf"))->GetXaxis()->SetBinLabel(5, "cosPaCascToXic");
      registry.get<TH1>(HIST("hSelStatusHf"))->GetXaxis()->SetBinLabel(6, "kfDcaXYPiFromXic");
      registry.get<TH1>(HIST("hSelStatusHf"))->GetXaxis()->SetBinLabel(7, "chi2GeoXic");
      registry.get<TH1>(HIST("hSelStatusHf"))->GetXaxis()->SetBinLabel(8, "chi2TopoPiFromXicToPv");
      registry.get<TH1>(HIST("hSelStatusHf"))->GetXaxis()->SetBinLabel(9, "chi2TopoCascToXic");
      registry.get<TH1>(HIST("hSelStatusHf"))->GetXaxis()->SetBinLabel(10, "decayLenXYXic");
      registry.get<TH1>(HIST("hSelStatusHf"))->GetXaxis()->SetBinLabel(11, "cTauXic");
    }

    // invarinat mass histograms
    registry.add("hInvMassCharmBaryonWoPidInvMassCut", "Charm baryon invariant mass; int mass; entries", {HistType::kTH1F, {{1500, 1.5, 4.5}}});
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

  // LF cuts - Cuts on LF tracks reco
  // Selection on LF related informations
  // returns true if all cuts are passed
  template <int svReco, typename T>
  bool SelectOnLF(const T& candidate, const int& inputPtBin)
  {

    registry.fill(HIST("hSelStatusLf"), 0.0);

    // Eta selection of V0, Cascade daughters
    double etaV0PosDau = candidate.etaV0PosDau();
    double etaV0NegDau = candidate.etaV0NegDau();
    double etaPiFromCasc = candidate.etaBachFromCasc();

    if (std::abs(etaV0PosDau) > etaTrackLFDauMax || std::abs(etaV0NegDau) > etaTrackLFDauMax || std::abs(etaPiFromCasc) > etaTrackLFDauMax) {
      return false;
    }
    registry.fill(HIST("hSelStatusLf"), 1.0);

    // Minimum radius cut
    double radiusV0 = RecoDecay::sqrtSumOfSquares(candidate.xDecayVtxV0(), candidate.yDecayVtxV0());
    double radiusCasc = RecoDecay::sqrtSumOfSquares(candidate.xDecayVtxCascade(), candidate.yDecayVtxCascade());

    if (radiusV0 < radiusV0Min) {
      return false;
    }
    registry.fill(HIST("hSelStatusLf"), 2.0);
    if (radiusCasc < radiusCascMin) {
      return false;
    }
    registry.fill(HIST("hSelStatusLf"), 3.0);

    // Cosine of pointing angle
    if (candidate.cosPAV0() < cuts->get(inputPtBin, "cosPAV0")) {
      return false;
    }
    registry.fill(HIST("hSelStatusLf"), 4.0);
    if (candidate.cosPACasc() < cuts->get(inputPtBin, "cosPACasc")) {
      return false;
    }
    registry.fill(HIST("hSelStatusLf"), 5.0);

    // Distance of Closest Approach(DCA)
    if (candidate.dcaV0Dau() > cuts->get(inputPtBin, "dcaV0Dau")) {
      return false;
    }
    registry.fill(HIST("hSelStatusLf"), 6.0);

    if (candidate.dcaCascDau() > cuts->get(inputPtBin, "dcaCascDau")) {
      return false;
    }
    registry.fill(HIST("hSelStatusLf"), 7.0);

    if (std::abs(candidate.dcaXYToPvV0Dau0()) < cuts->get(inputPtBin, "dcaXYToPvV0Dau0")) {
      return false;
    }
    registry.fill(HIST("hSelStatusLf"), 8.0);

    if (std::abs(candidate.dcaXYToPvV0Dau1()) < cuts->get(inputPtBin, "dcaXYToPvV0Dau1")) {
      return false;
    }
    registry.fill(HIST("hSelStatusLf"), 9.0);

    if (std::abs(candidate.dcaXYToPvCascDau()) < cuts->get(inputPtBin, "dcaXYToPvCascDau")) {
      return false;
    }
    registry.fill(HIST("hSelStatusLf"), 10.0);

    // pT: Bachelor
    double ptPiFromCasc = RecoDecay::sqrtSumOfSquares(candidate.pxBachFromCasc(), candidate.pyBachFromCasc());
    if (std::abs(ptPiFromCasc) < ptPiFromCascMin) {
      return false;
    }
    registry.fill(HIST("hSelStatusLf"), 11.0);

    // Extra cuts for KFParticle
    if constexpr (svReco == doKfParticle) {
      // Cosine of Pointing angle
      if (candidate.cosPaV0ToCasc() < cuts->get(inputPtBin, "cosPaV0ToCasc")) {
        return false;
      }
      registry.fill(HIST("hSelStatusLf"), 12.0);

      // DCA
      if (candidate.kfDcaXYCascToPv() > cuts->get(inputPtBin, "kfDcaXYCascToPv")) {
        return false;
      }
      registry.fill(HIST("hSelStatusLf"), 13.0);

      // Chi2
      if (candidate.chi2GeoV0() < 0 || candidate.chi2GeoV0() > cuts->get(inputPtBin, "chi2GeoV0")) {
        return false;
      }
      registry.fill(HIST("hSelStatusLf"), 14.0);
      if (candidate.chi2GeoCasc() < 0 || candidate.chi2GeoCasc() > cuts->get(inputPtBin, "chi2GeoCasc")) {
        return false;
      }
      registry.fill(HIST("hSelStatusLf"), 15.0);
      if (candidate.chi2TopoV0ToPv() > 0 && candidate.chi2TopoV0ToPv() < cuts->get(inputPtBin, "chi2TopoV0ToPv")) {
        return false;
      }
      registry.fill(HIST("hSelStatusLf"), 16.0);
      if (candidate.chi2TopoCascToPv() < 0 || candidate.chi2TopoCascToPv() > cuts->get(inputPtBin, "chi2TopoCascToPv")) {
        return false;
      }
      registry.fill(HIST("hSelStatusLf"), 17.0);
      if (candidate.chi2TopoV0ToCasc() < 0 || candidate.chi2TopoV0ToCasc() > cuts->get(inputPtBin, "chi2TopoV0ToCasc")) {
        return false;
      }
      registry.fill(HIST("hSelStatusLf"), 18.0);

      // ldl
      if (candidate.v0ldl() < cuts->get(inputPtBin, "v0ldl")) {
        return false;
      }
      registry.fill(HIST("hSelStatusLf"), 19.0);
      if (candidate.cascldl() < cuts->get(inputPtBin, "cascldl")) {
        return false;
      }
      registry.fill(HIST("hSelStatusLf"), 20.0);

      // Decay length
      if (std::abs(candidate.decayLenXYLambda()) < cuts->get(inputPtBin, "decayLenXYLambda")) {
        return false;
      }
      registry.fill(HIST("hSelStatusLf"), 21.0);
      if (std::abs(candidate.decayLenXYCasc()) < cuts->get(inputPtBin, "decayLenXYCasc")) {
        return false;
      }
      registry.fill(HIST("hSelStatusLf"), 22.0);

    } else {
      // Impact parameter(DCA?)
      if (std::abs(candidate.impactParCascXY()) < impactParXYCascMin || std::abs(candidate.impactParCascXY()) > impactParXYCascMax) {
        return false;
      }
      registry.fill(HIST("hSelStatusLf"), 12.0);
      if (std::abs(candidate.impactParCascZ()) < impactParZCascMin || std::abs(candidate.impactParCascZ()) > impactParZCascMax) {
        return false;
      }
      registry.fill(HIST("hSelStatusLf"), 13.0);
    }

    // If passes all cuts, return true
    return true;
  }

  // HF cuts - Cuts on Charm baryon reco
  // Apply cuts with charm baryon & charm bachelor related informations
  // returns true if all cuts are passed
  template <int svReco, typename T>
  bool SelectOnHF(const T& candidate, const int& inputPtBin)
  {

    registry.fill(HIST("hSelStatusHf"), 0.0);

    // eta selection on charm bayron bachelor
    if (std::abs(candidate.etaBachFromCharmBaryon()) > etaTrackCharmBachMax) {
      return false;
    }
    registry.fill(HIST("hSelStatusHf"), 1.0);
    // Distance of Closest Approach(DCA)
    if (candidate.dcaCharmBaryonDau() > cuts->get(inputPtBin, "dcaCharmBaryonDau")) {
      return false;
    }
    registry.fill(HIST("hSelStatusHf"), 2.0);

    // pT: Charm Bachelor
    double ptPiFromCharmBaryon = RecoDecay::pt(candidate.pxBachFromCharmBaryon(), candidate.pyBachFromCharmBaryon());
    if (std::abs(ptPiFromCharmBaryon) < ptPiFromCharmBaryonMin) {
      return false;
    }
    registry.fill(HIST("hSelStatusHf"), 3.0);

    // specific selections with KFParticle output
    if constexpr (svReco == doKfParticle) {
      // Cosine of pointing angle
      if (candidate.cosPaCascToXic() < cuts->get(inputPtBin, "cosPaCascToXic")) {
        return false;
      }
      registry.fill(HIST("hSelStatusHf"), 4.0);

      // DCA
      if (std::abs(candidate.kfDcaXYPiFromXic()) > cuts->get(inputPtBin, "kfDcaXYPiFromXic")) {
        return false;
      }
      registry.fill(HIST("hSelStatusHf"), 5.0);

      // Chi2
      if (candidate.chi2GeoXic() < 0 || candidate.chi2GeoXic() > cuts->get(inputPtBin, "chi2GeoXic")) {
        return false;
      }
      registry.fill(HIST("hSelStatusHf"), 6.0);
      if (candidate.chi2TopoPiFromXicToPv() < 0 || candidate.chi2TopoPiFromXicToPv() > cuts->get(inputPtBin, "chi2TopoXicToPv")) {
        return false;
      }
      registry.fill(HIST("hSelStatusHf"), 7.0);
      if (candidate.chi2TopoCascToXic() < 0 || candidate.chi2TopoCascToXic() > cuts->get(inputPtBin, "chi2TopoCascToXic")) {
        return false;
      }
      registry.fill(HIST("hSelStatusHf"), 8.0);

      // Decay Length
      if (std::abs(candidate.decayLenXYXic()) > cuts->get(inputPtBin, "decayLenXYXic")) {
        return false;
      }
      registry.fill(HIST("hSelStatusHf"), 9.0);

      // ctau
      if (std::abs(candidate.cTauXic()) > cuts->get(inputPtBin, "cTauXic")) {
        return false;
      }
      registry.fill(HIST("hSelStatusHf"), 10.0);
    } else {
      // Impact parameter(DCA?)
      if ((std::abs(candidate.impactParBachFromCharmBaryonXY()) < impactParXYPiFromCharmBaryonMin) || (std::abs(candidate.impactParBachFromCharmBaryonXY()) > impactParXYPiFromCharmBaryonMax)) {
        return false;
      }
      registry.fill(HIST("hSelStatusHf"), 4.0);
      if ((std::abs(candidate.impactParBachFromCharmBaryonZ()) < impactParZPiFromCharmBaryonMin) || (std::abs(candidate.impactParBachFromCharmBaryonZ()) > impactParZPiFromCharmBaryonMax)) {
        return false;
      }
      registry.fill(HIST("hSelStatusHf"), 5.0);
    }

    // If passes all cuts, return true
    return true;
  }

  template <int svReco, typename TCandTable>
  void runSelection(TCandTable const& candidates,
                    TracksSel const& tracks,
                    TracksSelLf const& lfTracks)
  {
    // looping over charm baryon candidates
    for (const auto& candidate : candidates) {

      bool resultSelections = true; // True if the candidate passes all the selections, False otherwise
      outputMlXic0ToXiPi.clear();

      auto trackV0PosDau = lfTracks.rawIteratorAt(candidate.posTrackId());
      auto trackV0NegDau = lfTracks.rawIteratorAt(candidate.negTrackId());
      auto trackPiFromCasc = lfTracks.rawIteratorAt(candidate.bachelorId());
      auto trackPiFromCharm = tracks.rawIteratorAt(candidate.bachelorFromCharmBaryonId());

      auto trackPiFromLam = trackV0NegDau;
      auto trackPrFromLam = trackV0PosDau;

      int8_t signDecay = candidate.signDecay(); // sign of pi <- cascade

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

      // Topological selection
      const bool selectionResOnLF = SelectOnLF<svReco>(candidate, pTBin);
      const bool selectionResOnHF = SelectOnHF<svReco>(candidate, pTBin);
      if (!selectionResOnLF || !selectionResOnHF) {
        resultSelections = false;
      }

      //  TPC clusters selections
      if (resultSelections) {
        registry.fill(HIST("hSelStatusCluster"), 0.0);
      }
      if (applyTrkSelLf) {
        if (!isSelectedTrackTpcQuality(trackPiFromLam, nClustersTpcMin, nTpcCrossedRowsMin, tpcCrossedRowsOverFindableClustersRatioMin, tpcChi2PerClusterMax)) {
          resultSelections = false;
        } else {
          if (resultSelections) {
            registry.fill(HIST("hSelStatusCluster"), 1.0);
          }
        }

        if (!isSelectedTrackTpcQuality(trackPrFromLam, nClustersTpcMin, nTpcCrossedRowsMin, tpcCrossedRowsOverFindableClustersRatioMin, tpcChi2PerClusterMax)) {
          resultSelections = false;
        } else {
          if (resultSelections) {
            registry.fill(HIST("hSelStatusCluster"), 2.0);
          }
        }

        if (!isSelectedTrackTpcQuality(trackPiFromCasc, nClustersTpcMin, nTpcCrossedRowsMin, tpcCrossedRowsOverFindableClustersRatioMin, tpcChi2PerClusterMax)) {
          resultSelections = false;
        } else {
          if (resultSelections) {
            registry.fill(HIST("hSelStatusCluster"), 3.0);
          }
        }
      }

      if (!isSelectedTrackTpcQuality(trackPiFromCharm, nClustersTpcMin, nTpcCrossedRowsMin, tpcCrossedRowsOverFindableClustersRatioMin, tpcChi2PerClusterMax)) {
        resultSelections = false;
      } else {
        if (resultSelections) {
          registry.fill(HIST("hSelStatusCluster"), 4.0);
        }
      }

      //  ITS clusters selection
      if (!isSelectedTrackItsQuality(trackPiFromCharm, nClustersItsMin, itsChi2PerClusterMax) || trackPiFromCharm.itsNClsInnerBarrel() < nClustersItsInnBarrMin) {
        resultSelections = false;
      } else {
        if (resultSelections) {
          registry.fill(HIST("hSelStatusCluster"), 5.0);
        }
      }

      // Track level PID selection
      if (resultSelections) {
        registry.fill(HIST("hSelStatusPID"), 0.0);
      }
      int statusPidPrFromLam = -999;
      int statusPidPiFromLam = -999;
      int statusPidPiFromCasc = -999;
      int statusPidPiFromCharmBaryon = -999;

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

      bool statusPidLambda = (statusPidPrFromLam == TrackSelectorPID::Accepted) && (statusPidPiFromLam == TrackSelectorPID::Accepted);
      if (statusPidLambda && resultSelections) {
        registry.fill(HIST("hSelStatusPID"), 1.0);
      }
      bool statusPidCascade = (statusPidLambda && statusPidPiFromCasc == TrackSelectorPID::Accepted);
      if (statusPidCascade && resultSelections) {
        registry.fill(HIST("hSelStatusPID"), 2.0);
      }
      bool statusPidCharmBaryon = (statusPidCascade && statusPidPiFromCharmBaryon == TrackSelectorPID::Accepted);
      if (statusPidCharmBaryon && resultSelections) {
        registry.fill(HIST("hSelStatusPID"), 3.0);
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
      }
      if (std::abs(invMassCascade - o2::constants::physics::MassXiMinus) < cascMassWindow) {
        statusInvMassCascade = true;
      }
      if ((invMassCharmBaryon >= invMassCharmBaryonMin) && (invMassCharmBaryon <= invMassCharmBaryonMax)) {
        statusInvMassCharmBaryon = true;
      }

      // ML BDT selection
      if (applyMl) {
        bool isSelectedMlXic0 = false;
        std::vector<float> inputFeaturesXic0 = {};
        if constexpr (svReco == doDcaFitter) {
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
      if constexpr (svReco == doDcaFitter) {
        hfSelToXiPi(statusPidLambda, statusPidCascade, statusPidCharmBaryon, statusInvMassLambda, statusInvMassCascade, statusInvMassCharmBaryon, resultSelections, infoTpcStored, infoTofStored,
                    trackPiFromCharm.tpcNSigmaPi(), trackPiFromCasc.tpcNSigmaPi(), trackPiFromLam.tpcNSigmaPi(), trackPrFromLam.tpcNSigmaPr(),
                    trackPiFromCharm.tofNSigmaPi(), trackPiFromCasc.tofNSigmaPi(), trackPiFromLam.tofNSigmaPi(), trackPrFromLam.tofNSigmaPr());
      } else {
        if (!statusPidCharmBaryon || !statusInvMassCharmBaryon) {
          resultSelections = false;
        }
        hfSelToXiPiKf(resultSelections,
                      trackPiFromCharm.tpcNSigmaPi(), trackPiFromCasc.tpcNSigmaPi(), trackPiFromLam.tpcNSigmaPi(), trackPrFromLam.tpcNSigmaPr(),
                      trackPiFromCharm.tofNSigmaPi(), trackPiFromCasc.tofNSigmaPi(), trackPiFromLam.tofNSigmaPi(), trackPrFromLam.tofNSigmaPr());
      }

      // Fill in invariant mass histogram
      if (resultSelections) {
        registry.fill(HIST("hInvMassCharmBaryonWoPidInvMassCut"), invMassCharmBaryon);
      }
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
  void processSelectionDCAFitter(aod::HfCandToXiPi const& candidates, TracksSel const& tracks, TracksSelLf const& lfTracks)
  {
    runSelection<doDcaFitter>(candidates, tracks, lfTracks);
  }
  PROCESS_SWITCH(HfCandidateSelectorToXiPiQa, processSelectionDCAFitter, "Xic0 candidate selection with DCAFitter output", true);

  ////////////////////////////////////
  ///    Process with KFParticle    //
  ////////////////////////////////////
  void processSelectionKFParticle(aod::HfCandToXiPiKf const& candidates, TracksSel const& tracks, TracksSelLf const& lfTracks)
  {
    runSelection<doKfParticle>(candidates, tracks, lfTracks);
  }
  PROCESS_SWITCH(HfCandidateSelectorToXiPiQa, processSelectionKFParticle, "Xic0 candidate selection with KFParticle output", false);

}; // end struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateSelectorToXiPiQa>(cfgc)};
}
