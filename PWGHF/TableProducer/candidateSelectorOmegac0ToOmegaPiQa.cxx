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

/// \file candidateSelectorOmegac0ToOmegaPiQa.cxx
/// \brief Omegac0 → Omega Pi selection task

/// \author Fabio Catalano <fabio.catalano@cern.ch>, University of Houston
/// \author Maria Fernanda Torres Cabrera <maria.fernanda.torres.cabrera@cern.ch>, University of Houston

#include "PWGHF/Core/HfMlResponseOmegacToOmegaPiQa.h"
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

#include <array>
#include <cstdint>
#include <cstdlib>
#include <numeric>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::analysis;

enum PidInfoStored {
  PiFromLam = 0,
  PrFromLam,
  KaFromCasc,
  PiFromCharm
};

enum {
  doDcaFitter = 0,
  doKfParticle
};

/// Struct for applying Omegac0 -> Omega pi selection cuts
struct HfCandidateSelectorToOmegaPiQa {
  // DCAFitter and KFParticle
  Produces<aod::HfSelToOmegaPi> hfSelToOmegaPi;
  // ML selection - filled for both DCAFitter and KFParticle
  Produces<aod::HfMlSelOmegacToOmegaPi> hfMlSelToOmegaPi;

  // cuts from SelectorCuts.h  - pT dependent cuts
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_omegac_to_omega_pi::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_omegac_to_omega_pi::Cuts[0], hf_cuts_omegac_to_omega_pi::NBinsPt, hf_cuts_omegac_to_omega_pi::NCutVars, hf_cuts_omegac_to_omega_pi::labelsPt, hf_cuts_omegac_to_omega_pi::labelsCutVar}, "OmegaC0 candidate selection per pT bin"};

  // ML inference
  Configurable<bool> applyMl{"applyMl", false, "Flag to apply ML selections"};
  Configurable<std::vector<double>> binsPtMl{"binsPtMl", std::vector<double>{hf_cuts_ml::vecBinsPt}, "pT bin limits for ML application"};
  Configurable<std::vector<int>> cutDirMl{"cutDirMl", std::vector<int>{hf_cuts_ml::vecCutDir}, "Whether to reject score values greater or smaller than the threshold"};
  Configurable<LabeledArray<double>> cutsMl{"cutsMl", {hf_cuts_ml::Cuts[0], hf_cuts_ml::NBinsPt, hf_cuts_ml::NCutScores, hf_cuts_ml::labelsPt, hf_cuts_ml::labelsCutScore}, "ML selections per pT bin"};
  Configurable<int> nClassesMl{"nClassesMl", static_cast<int>(hf_cuts_ml::NCutScores), "Number of classes in ML model"};
  Configurable<std::vector<std::string>> namesInputFeatures{"namesInputFeatures", std::vector<std::string>{"feature1", "feature2"}, "Names of ML model input features"};

  // CCDB configuration
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::vector<std::string>> modelPathsCCDB{"modelPathsCCDB", std::vector<std::string>{"EventFiltering/PWGHF/BDTOmegac0"}, "Paths of models on CCDB"};
  Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"ModelHandler_onnx_Omegac0ToOmegaPi.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};

  // LF analysis selections
  Configurable<double> radiusCascMin{"radiusCascMin", 0.5, "Min cascade radius"};
  Configurable<double> radiusV0Min{"radiusV0Min", 1.1, "Min V0 radius"};
  Configurable<double> cosPAV0Min{"cosPAV0Min", 0.97, "Min value CosPA V0"};
  Configurable<double> cosPACascMin{"cosPACascMin", 0.97, "Min value CosPA cascade"};
  Configurable<double> dcaV0DauMax{"dcaV0DauMax", 1.0, "Max DCA V0 daughters"};
  Configurable<double> dcaCascDauMax{"dcaCascDauMax", 1.0, "Max DCA cascade daughters"};
  Configurable<float> dcaNegToPvMin{"dcaNegToPvMin", 0.06, "DCA Neg To PV"};
  Configurable<float> dcaPosToPvMin{"dcaPosToPvMin", 0.06, "DCA Pos To PV"};
  Configurable<float> dcaBachToPvMin{"dcaBachToPvMin", 0.04, "DCA Bach To PV"};
  Configurable<bool> applyTrkSelLf{"applyTrkSelLf", true, "Apply track selection for LF daughters"};

  // Mass window
  Configurable<double> v0MassWindow{"v0MassWindow", 0.01, "V0 mass window"};
  Configurable<double> cascadeMassWindow{"cascadeMassWindow", 0.01, "Cascade mass window"};
  Configurable<double> invMassCharmBaryonMin{"invMassCharmBaryonMin", 2.3, "Lower limit invariant mass spectrum charm baryon"};
  Configurable<double> invMassCharmBaryonMax{"invMassCharmBaryonMax", 3.1, "Upper limit invariant mass spectrum charm baryon"};

  // kinematic selections
  Configurable<double> etaTrackCharmBachMax{"etaTrackCharmBachMax", 0.8, "Max absolute value of eta for charm baryon bachelor"};
  Configurable<double> etaTrackLFDauMax{"etaTrackLFDauMax", 1.0, "Max absolute value of eta for V0 and cascade daughters"};
  Configurable<double> ptKaFromCascMin{"ptKaFromCascMin", 0.15, "Min pT kaon <- casc"};

  Configurable<double> impactParameterXYPiFromCharmBaryonMin{"impactParameterXYPiFromCharmBaryonMin", 0., "Min dcaxy pi from charm baryon track to PV"};
  Configurable<double> impactParameterXYPiFromCharmBaryonMax{"impactParameterXYPiFromCharmBaryonMax", 10., "Max dcaxy pi from charm baryon track to PV"};
  Configurable<double> impactParameterXYCascMin{"impactParameterXYCascMin", 0., "Min dcaxy cascade track to PV"};
  Configurable<double> impactParameterXYCascMax{"impactParameterXYCascMax", 10., "Max dcaxy cascade track to PV"};
  Configurable<double> impactParameterZPiFromCharmBaryonMin{"impactParameterZPiFromCharmBaryonMin", 0., "Min dcaz pi from charm baryon track to PV"};
  Configurable<double> impactParameterZPiFromCharmBaryonMax{"impactParameterZPiFromCharmBaryonMax", 10., "Max dcaz pi from charm baryon track to PV"};
  Configurable<double> impactParameterZCascMin{"impactParameterZCascMin", 0., "Min dcaz cascade track to PV"};
  Configurable<double> impactParameterZCascMax{"impactParameterZCascMax", 10., "Max dcaz cascade track to PV"};

  Configurable<double> ptCandMin{"ptCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> ptCandMax{"ptCandMax", 50., "Upper bound of candidate pT"};

  Configurable<double> dcaCharmBaryonDauMax{"dcaCharmBaryonDauMax", 2.0, "Max DCA charm baryon daughters"};

  // PID options
  Configurable<bool> usePidTpcOnly{"usePidTpcOnly", true, "Perform PID using only TPC"};
  Configurable<bool> usePidTpcTofCombined{"usePidTpcTofCombined", false, "Perform PID using TPC or TOF"};

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
  Configurable<double> nSigmaTofCombinedKaMax{"nSigmaTofCombinedKaMax", 0., "Nsigma cut on TOF combined with TPC for kaon selection"};

  // detector clusters selections
  Configurable<int> nClustersTpcMin{"nClustersTpcMin", 70, "Minimum number of TPC clusters requirement"};
  Configurable<int> nTpcCrossedRowsMin{"nTpcCrossedRowsMin", 70, "Minimum number of TPC crossed rows requirement"};
  Configurable<double> tpcCrossedRowsOverFindableClustersRatioMin{"tpcCrossedRowsOverFindableClustersRatioMin", 0.8, "Minimum ratio TPC crossed rows over findable clusters requirement"};
  Configurable<float> tpcChi2PerClusterMax{"tpcChi2PerClusterMax", 4, "Maximum value of chi2 fit over TPC clusters"};
  Configurable<int> nClustersItsMin{"nClustersItsMin", 3, "Minimum number of ITS clusters requirement for pi <- charm baryon"};
  Configurable<int> nClustersItsInnBarrMin{"nClustersItsInnBarrMin", 1, "Minimum number of ITS clusters in inner barrel requirement for pi <- charm baryon"};
  Configurable<float> itsChi2PerClusterMax{"itsChi2PerClusterMax", 36, "Maximum value of chi2 fit over ITS clusters for pi <- charm baryon"};

  o2::analysis::HfMlResponseOmegacToOmegaPi<float, aod::hf_cand_casc_lf::ConstructMethod::DcaFitter> hfMlResponseDca;
  o2::analysis::HfMlResponseOmegacToOmegaPi<float, aod::hf_cand_casc_lf::ConstructMethod::KfParticle> hfMlResponseKf;
  std::vector<float> outputMlOmegac;
  o2::ccdb::CcdbApi ccdbApi;

  TrackSelectorPi selectorPion;
  TrackSelectorPr selectorProton;
  TrackSelectorKa selectorKaon;

  using TracksSel = soa::Join<aod::TracksWDcaExtra, aod::TracksPidPi, aod::TracksPidPr, aod::TracksPidKa>;
  using TracksSelLf = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksPidPi, aod::TracksPidPr, aod::TracksPidKa>;

  HistogramRegistry registry{"registry"}; // for QA of selections

  struct : ConfigurableGroup {
    //// KF selection
    std::string prefix = "kfSel";
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
  } kfConfigurableGroup;

  void init(InitContext const&)
  {
    std::array<bool, 2> processesSelector = {doprocessOmegac0SelectorWithDCAFitter, doprocessOmegac0SelectorWithKFParticle};
    const int nProcessesSelector = std::accumulate(processesSelector.begin(), processesSelector.end(), 0);
    if (nProcessesSelector != 1) {
      LOGP(fatal, "Exactly one process function for selector can be enabled at a time.");
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
    const AxisSpec axisSelOnLfDca{16, -0.5, 15.5, "status"};
    const AxisSpec axisSelOnLfKf{23, -0.5, 22.5, "status"};
    const AxisSpec axisSelOnHfDca{6, -0.5, 5.5, "status"};
    const AxisSpec axisSelOnHfKf{12, -0.5, 11.5, "status"};

    // for QA of the selections (bin 0 -> candidates that did not pass the selection, bin 1 -> candidates that passed the selection)
    registry.add("hSelSignDec", "hSelSignDec;status;entries", {HistType::kTH1F, {axisSel}});
    registry.add("hSelStatusCluster", "hSelStatusCluster:# of events Passed;;", {HistType::kTH1F, {{6, -0.5, 5.5}}});
    registry.get<TH1>(HIST("hSelStatusCluster"))->GetXaxis()->SetBinLabel(1, "All");
    registry.get<TH1>(HIST("hSelStatusCluster"))->GetXaxis()->SetBinLabel(2, "TpcCluster PiFromV0");
    registry.get<TH1>(HIST("hSelStatusCluster"))->GetXaxis()->SetBinLabel(3, "TpcCluster PrFromV0");
    registry.get<TH1>(HIST("hSelStatusCluster"))->GetXaxis()->SetBinLabel(4, "TpcCluster KaFromCasc");
    registry.get<TH1>(HIST("hSelStatusCluster"))->GetXaxis()->SetBinLabel(5, "TpcCluster PiFromCharm");
    registry.get<TH1>(HIST("hSelStatusCluster"))->GetXaxis()->SetBinLabel(6, "ItsCluster PiFromCharm");

    registry.add("hSelStatusPID", "hSelStatusPID;# of events Passed;;", {HistType::kTH1F, {{4, -0.5, 3.5}}});
    registry.get<TH1>(HIST("hSelStatusPID"))->GetXaxis()->SetBinLabel(1, "All");
    registry.get<TH1>(HIST("hSelStatusPID"))->GetXaxis()->SetBinLabel(2, "Lambda");
    registry.get<TH1>(HIST("hSelStatusPID"))->GetXaxis()->SetBinLabel(3, "Cascade");
    registry.get<TH1>(HIST("hSelStatusPID"))->GetXaxis()->SetBinLabel(4, "CharmBaryon");

    // For QA of LF & HF selection
    if (doprocessOmegac0SelectorWithDCAFitter) {
      registry.add("hSelStatusLf", "hSelStatusLf;# of candidate passed;", {HistType::kTH1F, {axisSelOnLfDca}});
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(1, "All");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(2, "etaV0PosDau");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(3, "etaV0NegDau");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(4, "etaKaFromCasc");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(5, "radiusV0");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(6, "radiusCasc");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(7, "cosPAV0");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(8, "cosPACasc");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(9, "dcaV0Dau");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(10, "dcaCascDau");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(11, "dcaXYToPvV0Dau0");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(12, "dcaXYToPvV0Dau1");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(13, "dcaXYToPvCascDau");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(14, "ptKaFromCasc");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(15, "impactParCascXY");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(16, "impactParCascZ");

      registry.add("hSelStatusHf", "hSelStatusHf;# of candidate passed;", {HistType::kTH1F, {axisSelOnHfDca}});
      registry.get<TH1>(HIST("hSelStatusHf"))->GetXaxis()->SetBinLabel(1, "All");
      registry.get<TH1>(HIST("hSelStatusHf"))->GetXaxis()->SetBinLabel(2, "etaTrackCharmBach");
      registry.get<TH1>(HIST("hSelStatusHf"))->GetXaxis()->SetBinLabel(3, "dcaCharmBaryonDau");
      registry.get<TH1>(HIST("hSelStatusHf"))->GetXaxis()->SetBinLabel(4, "ptPiFromCharmBaryon");
      registry.get<TH1>(HIST("hSelStatusHf"))->GetXaxis()->SetBinLabel(5, "impactParBachFromCharmXY");
      registry.get<TH1>(HIST("hSelStatusHf"))->GetXaxis()->SetBinLabel(6, "impactParBachFromCharmZ");
    }

    if (doprocessOmegac0SelectorWithKFParticle) {
      registry.add("hSelStatusLf", "hSelStatusLf;# of candidate passed;", {HistType::kTH1F, {axisSelOnLfKf}});
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(1, "All");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(2, "etaV0PosDau");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(3, "etaV0NegDau");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(4, "etaKaFromCasc");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(5, "radiusV0");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(6, "radiusCasc");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(7, "cosPAV0");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(8, "cosPACasc");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(9, "dcaV0Dau");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(10, "dcaCascDau");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(11, "dcaXYToPvV0Dau0");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(12, "dcaXYToPvV0Dau1");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(13, "dcaXYToPvCascDau");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(14, "ptKaFromCasc");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(15, "cosPaV0ToCasc");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(16, "v0Chi2OverNdf");
      registry.get<TH1>(HIST("hSelStatusLf"))->GetXaxis()->SetBinLabel(17, "cascChi2OverNdf");
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
      registry.get<TH1>(HIST("hSelStatusHf"))->GetXaxis()->SetBinLabel(5, "kfptOmegac");
      registry.get<TH1>(HIST("hSelStatusHf"))->GetXaxis()->SetBinLabel(6, "cosPaCascToOmegac");
      registry.get<TH1>(HIST("hSelStatusHf"))->GetXaxis()->SetBinLabel(7, "omegacChi2OverNdf");
      registry.get<TH1>(HIST("hSelStatusHf"))->GetXaxis()->SetBinLabel(8, "chi2TopoOmegacToPv");
      registry.get<TH1>(HIST("hSelStatusHf"))->GetXaxis()->SetBinLabel(9, "chi2TopoCascToOmegac");
      registry.get<TH1>(HIST("hSelStatusHf"))->GetXaxis()->SetBinLabel(10, "decayLenXYOmegac");
      registry.get<TH1>(HIST("hSelStatusHf"))->GetXaxis()->SetBinLabel(11, "cTauOmegac");
      registry.get<TH1>(HIST("hSelStatusHf"))->GetXaxis()->SetBinLabel(12, "omegacldl");
    }

    registry.add("hInvMassCharmBaryon", "Charm baryon invariant mass; inv. mass; entries", {HistType::kTH1F, {{500, 2.3, 3.1}}});
    registry.add("hPtCharmBaryon", "Charm baryon transverse momentum; p_{T} (GeV/#it{c}); entries", {HistType::kTH1F, {{8000, 0., 80.}}});

    if (doprocessOmegac0SelectorWithKFParticle) {
      registry.add("hSelCompetingCasc", "hSelCompetingCasc;status;entries", {HistType::kTH1F, {axisSel}});
      registry.add("hInvMassXiMinus_rej_cut", "hInvMassXiMinus_rej_cut;m_{#Lambda#pi} under Xi hypothesis (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{1000, 1.25f, 1.65f}}});
    }

    // HfMlResponse initialization
    if (applyMl) {
      if (doprocessOmegac0SelectorWithKFParticle) {
        registry.add("hBDTScoreKF", "hBDTScoreKF", {HistType::kTH1D, {{100, 0.0f, 1.0f, "score"}}});
        hfMlResponseKf.configure(binsPtMl, cutsMl, cutDirMl, nClassesMl);
        if (loadModelsFromCCDB) {
          ccdbApi.init(ccdbUrl);
          hfMlResponseKf.setModelPathsCCDB(onnxFileNames, ccdbApi, modelPathsCCDB, timestampCCDB);
        } else {
          hfMlResponseKf.setModelPathsLocal(onnxFileNames);
        }
        hfMlResponseKf.cacheInputFeaturesIndices(namesInputFeatures);
        hfMlResponseKf.init();
      } else if (doprocessOmegac0SelectorWithDCAFitter) {
        registry.add("hBDTScoreDCA", "hBDTScoreDCA", {HistType::kTH1D, {{100, 0.0f, 1.0f, "score"}}});
        hfMlResponseDca.configure(binsPtMl, cutsMl, cutDirMl, nClassesMl);
        if (loadModelsFromCCDB) {
          ccdbApi.init(ccdbUrl);
          hfMlResponseDca.setModelPathsCCDB(onnxFileNames, ccdbApi, modelPathsCCDB, timestampCCDB);
        } else {
          hfMlResponseDca.setModelPathsLocal(onnxFileNames);
        }
        hfMlResponseDca.cacheInputFeaturesIndices(namesInputFeatures);
        hfMlResponseDca.init();
      }
    }
  }

  // LF cuts - Cuts on LF tracks reco
  // Selection on LF related informations
  // returns true if all cuts are passed
  template <int svReco, typename T>
  bool selectOnLf(const T& candidate)
  {

    registry.fill(HIST("hSelStatusLf"), 0.0);

    // Eta selection of V0, Cascade daughters
    double etaV0PosDau = candidate.etaV0PosDau();
    double etaV0NegDau = candidate.etaV0NegDau();
    double etaKaFromCasc = candidate.etaBachFromCasc();

    if (std::abs(etaV0PosDau) > etaTrackLFDauMax) {
      return false;
    }
    registry.fill(HIST("hSelStatusLf"), 1.0);

    if (std::abs(etaV0NegDau) > etaTrackLFDauMax) {
      return false;
    }
    registry.fill(HIST("hSelStatusLf"), 2.0);

    if (std::abs(etaKaFromCasc) > etaTrackLFDauMax) {
      return false;
    }
    registry.fill(HIST("hSelStatusLf"), 3.0);

    // Minimum radius cut
    double radiusV0 = RecoDecay::sqrtSumOfSquares(candidate.xDecayVtxV0(), candidate.yDecayVtxV0());
    double radiusCasc = RecoDecay::sqrtSumOfSquares(candidate.xDecayVtxCascade(), candidate.yDecayVtxCascade());

    if (radiusV0 < radiusV0Min) {
      return false;
    }
    registry.fill(HIST("hSelStatusLf"), 4.0);
    if (radiusCasc < radiusCascMin) {
      return false;
    }
    registry.fill(HIST("hSelStatusLf"), 5.0);

    // Cosine of pointing angle
    if (candidate.cosPAV0() < cosPAV0Min) {
      return false;
    }
    registry.fill(HIST("hSelStatusLf"), 6.0);
    if (candidate.cosPACasc() < cosPACascMin) {
      return false;
    }
    registry.fill(HIST("hSelStatusLf"), 7.0);

    // Distance of Closest Approach(DCA)
    if (candidate.dcaV0Dau() > dcaV0DauMax) {
      return false;
    }
    registry.fill(HIST("hSelStatusLf"), 8.0);

    if (candidate.dcaCascDau() > dcaCascDauMax) {
      return false;
    }
    registry.fill(HIST("hSelStatusLf"), 9.0);

    if (std::abs(candidate.dcaXYToPvV0Dau0()) < dcaPosToPvMin) {
      return false;
    }
    registry.fill(HIST("hSelStatusLf"), 10.0);

    if (std::abs(candidate.dcaXYToPvV0Dau1()) < dcaNegToPvMin) {
      return false;
    }
    registry.fill(HIST("hSelStatusLf"), 11.0);

    if (std::abs(candidate.dcaXYToPvCascDau()) < dcaBachToPvMin) {
      return false;
    }
    registry.fill(HIST("hSelStatusLf"), 12.0);

    // pT: Bachelor
    double ptKaFromCasc = RecoDecay::sqrtSumOfSquares(candidate.pxBachFromCasc(), candidate.pyBachFromCasc());
    if (std::abs(ptKaFromCasc) < ptKaFromCascMin) {
      return false;
    }
    registry.fill(HIST("hSelStatusLf"), 13.0);

    // Extra cuts for KFParticle
    if constexpr (svReco == doKfParticle) {
      if (kfConfigurableGroup.applyKFpreselections) {
        // Cosine of Pointing angle
        if (candidate.cosPaV0ToCasc() < kfConfigurableGroup.cosPaV0ToCascMin) {
          return false;
        }
        registry.fill(HIST("hSelStatusLf"), 14.0);

        // Chi2
        if (candidate.v0Chi2OverNdf() < 0 || candidate.v0Chi2OverNdf() > kfConfigurableGroup.v0Chi2OverNdfMax) {
          return false;
        }
        registry.fill(HIST("hSelStatusLf"), 15.0);
        if (candidate.cascChi2OverNdf() < 0 || candidate.cascChi2OverNdf() > kfConfigurableGroup.cascChi2OverNdfMax) {
          return false;
        }
        registry.fill(HIST("hSelStatusLf"), 16.0);
        if (candidate.chi2TopoCascToPv() < 0 || candidate.chi2TopoCascToPv() > kfConfigurableGroup.chi2TopoCascToPvMax) {
          return false;
        }
        registry.fill(HIST("hSelStatusLf"), 17.0);
        if (candidate.chi2TopoV0ToCasc() < 0 || candidate.chi2TopoV0ToCasc() > kfConfigurableGroup.chi2TopoV0ToCascMax) {
          return false;
        }
        registry.fill(HIST("hSelStatusLf"), 18.0);

        // ldl
        if (candidate.v0ldl() < kfConfigurableGroup.v0LdlMin) {
          return false;
        }
        registry.fill(HIST("hSelStatusLf"), 19.0);
        if (candidate.cascldl() < kfConfigurableGroup.cascLdlMin) {
          return false;
        }
        registry.fill(HIST("hSelStatusLf"), 20.0);

        // Decay length
        if (std::abs(candidate.decayLenXYLambda()) < kfConfigurableGroup.decayLenXYLambdaMin) {
          return false;
        }
        registry.fill(HIST("hSelStatusLf"), 21.0);
        if (std::abs(candidate.decayLenXYCasc()) < kfConfigurableGroup.decayLenXYCascMin) {
          return false;
        }
        registry.fill(HIST("hSelStatusLf"), 22.0);

        // Competing Xi rejection
        if (kfConfigurableGroup.applyCompetingCascRejection) {
          const auto invMassXiHypothesis = candidate.cascRejectInvmass();

          if (std::abs(invMassXiHypothesis - o2::constants::physics::MassXiMinus) < kfConfigurableGroup.cascadeRejMassWindow) {
            registry.fill(HIST("hSelCompetingCasc"), 0.0);
            return false;
          }

          registry.fill(HIST("hSelCompetingCasc"), 1.0);
          registry.fill(HIST("hInvMassXiMinus_rej_cut"), invMassXiHypothesis);
        }
      }
    } else {
      // Impact parameter
      if (std::abs(candidate.impactParCascXY()) < impactParameterXYCascMin || std::abs(candidate.impactParCascXY()) > impactParameterXYCascMax) {
        return false;
      }
      registry.fill(HIST("hSelStatusLf"), 14.0);
      if (std::abs(candidate.impactParCascZ()) < impactParameterZCascMin || std::abs(candidate.impactParCascZ()) > impactParameterZCascMax) {
        return false;
      }
      registry.fill(HIST("hSelStatusLf"), 15.0);
    }

    // If passes all cuts, return true
    return true;
  }

  // HF cuts - Cuts on Charm baryon reco
  // Apply cuts with charm baryon & charm bachelor related informations
  // returns true if all cuts are passed
  template <int svReco, typename T>
  bool selectOnHf(const T& candidate, const int& inputPtBin)
  {
    registry.fill(HIST("hSelStatusHf"), 0.0);

    // eta selection on charm bayron bachelor
    if (std::abs(candidate.etaBachFromCharmBaryon()) > etaTrackCharmBachMax) {
      return false;
    }
    registry.fill(HIST("hSelStatusHf"), 1.0);
    // Distance of Closest Approach(DCA)
    if (candidate.dcaCharmBaryonDau() > dcaCharmBaryonDauMax) {
      return false;
    }
    registry.fill(HIST("hSelStatusHf"), 2.0);

    // pT: Charm Bachelor
    double ptPiFromCharmBaryon = RecoDecay::sqrtSumOfSquares(candidate.pxBachFromCharmBaryon(), candidate.pyBachFromCharmBaryon());
    if (inputPtBin < 0 || ptPiFromCharmBaryon < cuts->get(inputPtBin, "pT pi from Omegac")) {
      return false;
    }
    registry.fill(HIST("hSelStatusHf"), 3.0);

    // specific selections with KFParticle output
    if constexpr (svReco == doKfParticle) {
      if (kfConfigurableGroup.applyKFpreselections) {
        // Omegac Pt selection
        if (std::abs(candidate.kfptOmegac()) <= ptCandMin || std::abs(candidate.kfptOmegac()) >= ptCandMax) {
          return false;
        }
        registry.fill(HIST("hSelStatusHf"), 4.0);

        // Cosine of pointing angle
        if (candidate.cosPaCascToOmegac() < kfConfigurableGroup.cosPaCascToOmegacMin) {
          return false;
        }
        registry.fill(HIST("hSelStatusHf"), 5.0);

        // Chi2
        if (candidate.omegacChi2OverNdf() < 0 || candidate.omegacChi2OverNdf() > kfConfigurableGroup.omegacChi2OverNdfMax) {
          return false;
        }
        registry.fill(HIST("hSelStatusHf"), 6.0);
        if (candidate.chi2TopoOmegacToPv() < 0 || candidate.chi2TopoOmegacToPv() > kfConfigurableGroup.chi2TopoOmegacToPvMax) {
          return false;
        }
        registry.fill(HIST("hSelStatusHf"), 7.0);
        if (candidate.chi2TopoCascToOmegac() < 0 || candidate.chi2TopoCascToOmegac() > kfConfigurableGroup.chi2TopoCascToOmegacMax) {
          return false;
        }
        registry.fill(HIST("hSelStatusHf"), 8.0);

        // Decay Length
        if (std::abs(candidate.decayLenXYOmegac()) > kfConfigurableGroup.decayLenXYOmegacMax) {
          return false;
        }
        registry.fill(HIST("hSelStatusHf"), 9.0);

        // ctau
        if (std::abs(candidate.cTauOmegac()) > kfConfigurableGroup.cTauOmegacMax) {
          return false;
        }
        registry.fill(HIST("hSelStatusHf"), 10.0);

        // Omegac l/dl
        if (candidate.omegacldl() > kfConfigurableGroup.omegacLdlMax) {
          return false;
        }
        registry.fill(HIST("hSelStatusHf"), 11.0);
      }
    } else {
      // Impact parameter
      if ((std::abs(candidate.impactParBachFromCharmBaryonXY()) < impactParameterXYPiFromCharmBaryonMin) || (std::abs(candidate.impactParBachFromCharmBaryonXY()) > impactParameterXYPiFromCharmBaryonMax)) {
        return false;
      }
      registry.fill(HIST("hSelStatusHf"), 4.0);
      if ((std::abs(candidate.impactParBachFromCharmBaryonZ()) < impactParameterZPiFromCharmBaryonMin) || (std::abs(candidate.impactParBachFromCharmBaryonZ()) > impactParameterZPiFromCharmBaryonMax)) {
        return false;
      }
      registry.fill(HIST("hSelStatusHf"), 5.0);
    }

    // If passes all cuts, return true
    return true;
  }

  template <int svReco, typename TCandTable>
  void runOmegac0Selector(TCandTable const& candidates,
                          TracksSel const& tracks,
                          TracksSelLf const& lfTracks)
  {
    // looping over charm baryon candidates
    for (const auto& candidate : candidates) {

      bool resultSelections = true; // True if the candidate passes all the selections, False otherwise
      outputMlOmegac.clear();

      auto trackV0PosDau = lfTracks.rawIteratorAt(candidate.posTrackId());
      auto trackV0NegDau = lfTracks.rawIteratorAt(candidate.negTrackId());
      auto trackKaFromCasc = lfTracks.rawIteratorAt(candidate.bachelorId());
      auto trackPiFromCharm = tracks.rawIteratorAt(candidate.bachelorFromCharmBaryonId());

      auto trackPiFromLam = trackV0NegDau;
      auto trackPrFromLam = trackV0PosDau;

      int8_t const signDecay = candidate.signDecay(); // sign of pi <- cascade

      if (signDecay > 0) {
        trackPiFromLam = trackV0PosDau;
        trackPrFromLam = trackV0NegDau;
        registry.fill(HIST("hSelSignDec"), 1); // anti-particle decay
      } else {
        registry.fill(HIST("hSelSignDec"), 0); // particle decay
      }

      // pT selection
      auto ptCandOmegac = RecoDecay::pt(candidate.pxCharmBaryon(), candidate.pyCharmBaryon());

      if (ptCandOmegac < ptCandMin || ptCandOmegac > ptCandMax) {
        resultSelections = false;
      }

      int pTBin = findBin(binsPt, ptCandOmegac);
      if (pTBin == -1) {
        resultSelections = false;
      }

      // Topological selection
      const bool selectionResOnLF = selectOnLf<svReco>(candidate);
      const bool selectionResOnHF = selectOnHf<svReco>(candidate, pTBin);
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

        if (!isSelectedTrackTpcQuality(trackKaFromCasc, nClustersTpcMin, nTpcCrossedRowsMin, tpcCrossedRowsOverFindableClustersRatioMin, tpcChi2PerClusterMax)) {
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
      int statusPidKaFromCasc = -999;
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
      if (trackKaFromCasc.hasTPC()) {
        SETBIT(infoTpcStored, KaFromCasc);
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
      if (trackKaFromCasc.hasTOF()) {
        SETBIT(infoTofStored, KaFromCasc);
      }
      if (trackPiFromCharm.hasTOF()) {
        SETBIT(infoTofStored, PiFromCharm);
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

      bool statusPidLambda = (statusPidPrFromLam == TrackSelectorPID::Accepted) && (statusPidPiFromLam == TrackSelectorPID::Accepted);
      if (statusPidLambda && resultSelections) {
        registry.fill(HIST("hSelStatusPID"), 1.0);
      }
      bool statusPidCascade = (statusPidLambda && statusPidKaFromCasc == TrackSelectorPID::Accepted);
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

      double const invMassLambda = candidate.invMassLambda();
      double const invMassCascade = candidate.invMassCascade();
      double const invMassCharmBaryon = candidate.invMassCharmBaryon();

      if (std::abs(invMassLambda - o2::constants::physics::MassLambda0) < v0MassWindow) {
        statusInvMassLambda = true;
      }
      if (std::abs(invMassCascade - o2::constants::physics::MassOmegaMinus) < cascadeMassWindow) {
        statusInvMassCascade = true;
      }
      if ((invMassCharmBaryon >= invMassCharmBaryonMin) && (invMassCharmBaryon <= invMassCharmBaryonMax)) {
        statusInvMassCharmBaryon = true;
      }

      // Fill in selection result
      if (!statusPidLambda || !statusPidCascade || !statusPidCharmBaryon ||
          !statusInvMassLambda || !statusInvMassCascade || !statusInvMassCharmBaryon) {
        resultSelections = false;
      }

      // Check candidate pT range for ML inference
      if (applyMl && findBin(binsPtMl, ptCandOmegac) == -1) {
        resultSelections = false;
      }

      // ML BDT selection
      if (applyMl && resultSelections) {
        bool isSelectedMlOmegac = false;
        std::vector<float> inputFeaturesOmegaC = {};

        if constexpr (svReco == doKfParticle) {
          inputFeaturesOmegaC = hfMlResponseKf.getInputFeatures(candidate, trackPiFromLam, trackKaFromCasc, trackPiFromCharm);
          isSelectedMlOmegac = hfMlResponseKf.isSelectedMl(inputFeaturesOmegaC, ptCandOmegac, outputMlOmegac);
          if (isSelectedMlOmegac) {
            registry.fill(HIST("hBDTScoreKF"), outputMlOmegac[0]);
          } else {
            resultSelections = false;
          }
        } else if constexpr (svReco == doDcaFitter) {
          inputFeaturesOmegaC = hfMlResponseDca.getInputFeatures(candidate, trackPiFromLam, trackKaFromCasc, trackPiFromCharm);
          isSelectedMlOmegac = hfMlResponseDca.isSelectedMl(inputFeaturesOmegaC, ptCandOmegac, outputMlOmegac);
          if (isSelectedMlOmegac) {
            registry.fill(HIST("hBDTScoreDCA"), outputMlOmegac[0]);
          } else {
            resultSelections = false;
          }
        }
      }

      if (applyMl) {
        hfMlSelToOmegaPi(outputMlOmegac);
      }

      hfSelToOmegaPi(statusPidLambda, statusPidCascade, statusPidCharmBaryon, statusInvMassLambda, statusInvMassCascade, statusInvMassCharmBaryon, resultSelections, infoTpcStored, infoTofStored,
                     trackPiFromCharm.tpcNSigmaPi(), trackKaFromCasc.tpcNSigmaKa(), trackPiFromLam.tpcNSigmaPi(), trackPrFromLam.tpcNSigmaPr(),
                     trackPiFromCharm.tofNSigmaPi(), trackKaFromCasc.tofNSigmaKa(), trackPiFromLam.tofNSigmaPi(), trackPrFromLam.tofNSigmaPr());

      // Fill in invariant mass histogram
      if (statusPidLambda && statusPidCascade && statusPidCharmBaryon && statusInvMassLambda && statusInvMassCascade && statusInvMassCharmBaryon && resultSelections) {
        registry.fill(HIST("hInvMassCharmBaryon"), invMassCharmBaryon);
        if constexpr (svReco == doKfParticle) {
          registry.fill(HIST("hPtCharmBaryon"), candidate.kfptOmegac());
        } else {
          registry.fill(HIST("hPtCharmBaryon"), ptCandOmegac);
        }
      }

    } // end of candidate loop
  } // end run function

  ///////////////////////////////////
  ///    Process with DCAFitter    //
  ///////////////////////////////////
  void processOmegac0SelectorWithDCAFitter(aod::HfCandToOmegaPi const& candidates,
                                           TracksSel const& tracks,
                                           TracksSelLf const& lfTracks)
  {
    runOmegac0Selector<doDcaFitter>(candidates, tracks, lfTracks);
  }
  PROCESS_SWITCH(HfCandidateSelectorToOmegaPiQa, processOmegac0SelectorWithDCAFitter, "Omegac0 candidate selection with DCAFitter output", true);

  ////////////////////////////////////
  ///    Process with KFParticle    //
  ////////////////////////////////////
  void processOmegac0SelectorWithKFParticle(soa::Join<aod::HfCandToOmegaPi, aod::HfOmegacKf> const& candidates,
                                            TracksSel const& tracks,
                                            TracksSelLf const& lfTracks)
  {
    runOmegac0Selector<doKfParticle>(candidates, tracks, lfTracks);
  }
  PROCESS_SWITCH(HfCandidateSelectorToOmegaPiQa, processOmegac0SelectorWithKFParticle, "Omegac0 candidate selection with KFParticle output", false);

}; // end struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateSelectorToOmegaPiQa>(cfgc)};
}
