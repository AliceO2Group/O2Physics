// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file mlBasedTrackSelector.cxx
/// \brief D± → π± K∓ π± and D_s± → K± K∓ π± skimming with track-level ML
/// TEST ONLY, PLEASE DO NOT USE FOR ANALYSIS
/// \author Fabrizio Chinu <fabrizio.chinu@cern.ch>, Universita and INFN Torino

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/HfMlResponseHfTracks.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"
#include "PWGHF/Utils/utilsAnalysis.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"
#include "PWGHF/Utils/utilsEvSelHf.h"

#include "Common/CCDB/TriggerAliases.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Tools/ML/MlResponse.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DCAFitter/DCAFitterN.h>
#include <DetectorsBase/MatLayerCylSet.h>
#include <DetectorsBase/Propagator.h>
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
#include <ReconstructionDataFormats/DCA.h>
#include <ReconstructionDataFormats/Track.h>
#include <ReconstructionDataFormats/Vertex.h>

#include <TH1.h>

#include <Rtypes.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::analysis;
using namespace o2::hf_evsel;
using namespace o2::aod;
using namespace o2::hf_centrality;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

/// Bit positions of aod::HfSelTrack::isSelProng (only Cand3Prong is ever set here).
enum CandidateType {
  Cand2Prong = 0,
  Cand3Prong,
  CandV0bachelor,
  CandDstar,
  CandCascadeBachelor,
  NCandidateTypes
};

/// The two 3-prong channels handled by this task.
enum Channel3Prong {
  ChannelDplusToPiKPi = 0,
  ChannelDsToKKPi,
  NChannels3Prong
};

constexpr std::array<int, NChannels3Prong> DecayTypeOfChannel{
  hf_cand_3prong::DecayType::DplusToPiKPi,
  hf_cand_3prong::DecayType::DsToKKPi};

/// Bit positions of aod::HfSelTrack::isIdentifiedPid filled by this task: the decision of the
/// track-level model, one bit per (channel, role). A track can be pion-like for one channel and
/// kaon-like for the other, so the four bits are independent.
enum TrackMlRole {
  RolePionDplus = 0,
  RoleKaonDplus,
  RolePionDs,
  RoleKaonDs,
  NTrackMlRoles
};

constexpr std::array<int, NChannels3Prong> RolePionOfChannel{RolePionDplus, RolePionDs};
constexpr std::array<int, NChannels3Prong> RoleKaonOfChannel{RoleKaonDplus, RoleKaonDs};

constexpr int NProngs = 3;    // prongs of the candidates built here
constexpr int NMassHypos = 2; // mass hypotheses per channel, i.e. the two orderings of the same-sign pair

/// Whether each prong is a kaon, per channel and mass hypothesis.
constexpr std::array<std::array<std::array<bool, NProngs>, NMassHypos>, NChannels3Prong> IsKaonProng{{
  {{{{false, true, false}}, {{false, true, false}}}},
  {{{{true, true, false}}, {{false, true, true}}}},
}};

constexpr double PtMaxModel = 1.e10;

constexpr uint32_t MaskSameSignAny = (1u << RolePionDplus) | (1u << RolePionDs) | (1u << RoleKaonDs);
constexpr uint32_t MaskOppSignAny = (1u << RoleKaonDplus) | (1u << RoleKaonDs);

/// Where the combinatorics spends its time, accumulated in seconds into hTiming.
enum TimingStep {
  TimeLoopTotal = 0,  ///< everything inside the triple loop
  TimeFit2Prong,      ///< the 2-track vertex used to skip pairs early
  TimeFit3Prong,      ///< the 3-track vertex of the surviving triplets
  TimeCacheFill,      ///< building the per-collision prong cache, propagation included
  TimeCachePropagate, ///< only the re-propagation inside that build
  NTimingSteps
};

/// How often each stage of the triple loop runs, accumulated into hLoopCounters.
enum LoopCounter {
  CountPairsSeen = 0,      ///< (same-sign, opposite-sign) pairs reached
  CountPairsCandidateOk,   ///< ... of those, surviving the ML pair test
  CountFit2Prong,          ///< ... of those, handed to the 2-prong fitter
  CountPairsRejected2P,    ///< ... of those, dropped by the 2-track vertex
  CountTripletsEnumerated, ///< triplets built from the surviving pairs
  CountTripletsWritten,    ///< ... of those, written to the skim
  CountTracksCached,       ///< associations put into the prong cache
  CountTracksPropagated,   ///< ... of those, needing a re-propagation to this collision
  NLoopCounters
};

/// Where the ML track selection spends its time, accumulated in seconds into its own hTiming.
enum TrackTimingStep {
  TrackTimeTotal = 0,  ///< the whole per-collision loop: slicing, features, models, table filling
  TrackTimeFeatures,   ///< building the input features, re-propagation included
  TrackTimeModelDplus, ///< evaluating the D+ model: input vector, ONNX call, thresholds
  TrackTimeModelDs,    ///< evaluating the Ds model
  NTrackTimingSteps
};

inline std::chrono::steady_clock::time_point tickIf(const bool enabled)
{
  return enabled ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point{};
}
inline double elapsedSeconds(std::chrono::steady_clock::time_point const& start)
{
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

/// Event selection
struct HfTrackSelectorTagSelCollisions {
  Produces<aod::HfSelCollision> rowSelectedCollision;

  Configurable<bool> fillHistograms{"fillHistograms", true, "fill histograms"};
  Configurable<std::string> triggerClassName{"triggerClassName", "kINT7", "Run 2 trigger class, only for Run 2 converted data"};
  HfEventSelection hfEvSel;                   // event selection and monitoring
  Service<o2::ccdb::BasicCCDBManager> ccdb{}; // needed for evSelection

  // QA histos
  HistogramRegistry registry{"registry"};
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  void init(InitContext const&)
  {
    const std::array<bool, 7> doProcess = {doprocessTrigAndCentFT0ASel, doprocessTrigAndCentFT0CSel, doprocessTrigAndCentFT0MSel, doprocessTrigAndCentFV0ASel, doprocessTrigSel, doprocessNoTrigSel, doprocessUpcSel};
    if (std::count(doProcess.begin(), doProcess.end(), true) != 1) {
      LOGP(fatal, "One and only one process function for collision selection can be enabled at a time!");
    }

    // set numerical value of the Run 2 trigger class
    auto* const triggerAlias = std::find(aliasLabels, aliasLabels + kNaliases, triggerClassName.value.data());
    if (triggerAlias != aliasLabels + kNaliases) {
      hfEvSel.triggerClass.value = std::distance(aliasLabels, triggerAlias);
    }

    hfEvSel.init(registry, &zorroSummary); // collision monitoring
    if (fillHistograms) {
      if (doprocessTrigAndCentFT0ASel || doprocessTrigAndCentFT0CSel || doprocessTrigAndCentFT0MSel || doprocessTrigAndCentFV0ASel) {
        const AxisSpec axisCentrality{200, 0., 100., "centrality percentile"};
        registry.add("hCentralitySelected", "Centrality percentile of selected events in the centrality interval; centrality percentile;entries", {HistType::kTH1D, {axisCentrality}});
        registry.add("hCentralityRejected", "Centrality percentile of selected events outside the centrality interval; centrality percentile;entries", {HistType::kTH1D, {axisCentrality}});
      }
    }
  }

  /// Collision selection
  /// \param collision  collision table with
  template <bool ApplyTrigSel, bool ApplyUpcSel, o2::hf_centrality::CentralityEstimator CentEstimator, typename Col, typename BCsType>
  void selectCollision(const Col& collision,
                       const BCsType& bcs)
  {
    float centrality{-1.f};
    o2::hf_evsel::HfCollisionRejectionMask rejectionMask{};

    if constexpr (ApplyUpcSel) {
      rejectionMask = hfEvSel.getHfCollisionRejectionMaskWithUpc<ApplyTrigSel, CentEstimator>(
        collision, centrality, ccdb, registry, bcs);
    } else {
      rejectionMask = hfEvSel.getHfCollisionRejectionMask<ApplyTrigSel, CentEstimator, BCsType>(
        collision, centrality, ccdb, registry);
    }

    if (fillHistograms) {
      hfEvSel.fillHistograms(collision, rejectionMask, centrality);
      // additional centrality histos
      if constexpr (CentEstimator != o2::hf_centrality::None) {
        if (rejectionMask == 0) {
          registry.fill(HIST("hCentralitySelected"), centrality);
        } else if (rejectionMask == BIT(EventRejection::Centrality)) { // rejected by centrality only
          registry.fill(HIST("hCentralityRejected"), centrality);
        }
      }
    }

    // fill table row
    rowSelectedCollision(rejectionMask);
  }

  /// Event selection with trigger and FT0A centrality selection
  void processTrigAndCentFT0ASel(soa::Join<aod::Collisions,
                                           aod::EvSels, aod::PVMults, aod::CentFT0As>::iterator const& collision,
                                 aod::BcFullInfos const& bcs)
  {
    selectCollision<true, false, CentralityEstimator::FT0A>(collision, bcs);
  }
  PROCESS_SWITCH(HfTrackSelectorTagSelCollisions, processTrigAndCentFT0ASel, "Use trigger and centrality selection with FT0A", false);

  /// Event selection with trigger and FT0C centrality selection
  void processTrigAndCentFT0CSel(soa::Join<aod::Collisions,
                                           aod::EvSels, aod::PVMults, aod::CentFT0Cs>::iterator const& collision,
                                 aod::BcFullInfos const& bcs)
  {
    selectCollision<true, false, CentralityEstimator::FT0C>(collision, bcs);
  }
  PROCESS_SWITCH(HfTrackSelectorTagSelCollisions, processTrigAndCentFT0CSel, "Use trigger and centrality selection with FT0C", false);

  /// Event selection with trigger and FT0M centrality selection
  void processTrigAndCentFT0MSel(soa::Join<aod::Collisions,
                                           aod::EvSels, aod::PVMults, aod::CentFT0Ms>::iterator const& collision,
                                 aod::BcFullInfos const& bcs)
  {
    selectCollision<true, false, CentralityEstimator::FT0M>(collision, bcs);
  }
  PROCESS_SWITCH(HfTrackSelectorTagSelCollisions, processTrigAndCentFT0MSel, "Use trigger and centrality selection with FT0M", false);

  /// Event selection with trigger and FV0A centrality selection
  void processTrigAndCentFV0ASel(soa::Join<aod::Collisions,
                                           aod::EvSels, aod::PVMults, aod::CentFV0As>::iterator const& collision,
                                 aod::BcFullInfos const& bcs)
  {
    selectCollision<true, false, CentralityEstimator::FV0A>(collision, bcs);
  }
  PROCESS_SWITCH(HfTrackSelectorTagSelCollisions, processTrigAndCentFV0ASel, "Use trigger and centrality selection with FV0A", false);

  /// Event selection with trigger selection
  void processTrigSel(soa::Join<aod::Collisions,
                                aod::EvSels, aod::PVMults>::iterator const& collision,
                      aod::BcFullInfos const& bcs)
  {
    selectCollision<true, false, CentralityEstimator::None>(collision, bcs);
  }
  PROCESS_SWITCH(HfTrackSelectorTagSelCollisions, processTrigSel, "Use trigger selection", false);

  /// Event selection without trigger selection
  void processNoTrigSel(soa::Join<aod::Collisions, aod::PVMults>::iterator const& collision,
                        aod::BcFullInfos const& bcs)
  {
    selectCollision<false, false, CentralityEstimator::None>(collision, bcs);
  }
  PROCESS_SWITCH(HfTrackSelectorTagSelCollisions, processNoTrigSel, "Do not use trigger selection", true);

  /// Event selection with UPC
  void processUpcSel(soa::Join<aod::Collisions, aod::EvSels, aod::PVMults>::iterator const& collision,
                     aod::BcFullInfos const& bcs,
                     aod::FT0s const& /*ft0s*/,
                     aod::FV0As const& /*fv0as*/,
                     aod::FDDs const& /*fdds*/,
                     aod::Zdcs const& /*zdcs*/)
  {
    selectCollision<true, true, CentralityEstimator::None>(collision, bcs);
  }
  PROCESS_SWITCH(HfTrackSelectorTagSelCollisions, processUpcSel, "Use UPC event selection", false);
};

/// Track selection, ML based
struct HfTrackSelectorTagSelTracks {
  Produces<aod::HfSelTrack> rowSelectedTrack;

  struct : ConfigurableGroup {
    Configurable<bool> testAcknowledgement{"testAcknowledgement", false, "test acknowledgement"};
    Configurable<bool> fillHistograms{"fillHistograms", true, "fill histograms"};
    Configurable<float> ptMinTrack{"ptMinTrack", 0.3f, "min. track pT entering the charm combinatorics"};
    Configurable<bool> enableTiming{"enableTiming", false, "fill hTiming with the CPU of the feature building and of each model evaluation (adds two clock reads per call)"};
    // D+ model
    Configurable<bool> applyMlDplus{"applyMlDplus", true, "evaluate the D+ track model"};
    Configurable<std::string> onnxFileNameDplus{"onnxFileNameDplus", "ModelHandler_DplusTracks.onnx", "ONNX file name of the D+ track model"};
    Configurable<std::vector<std::string>> inputFeaturesDplus{"inputFeaturesDplus", std::vector<std::string>{"pt", "eta", "dcaXY", "dcaZ", "sigmaDcaXY", "sigmaDcaZ", "normDcaXY", "normDcaZ", "signed1Pt", "tgl", "sign", "isPvContributor", "itsNCls", "itsNClsInnerBarrel", "itsChi2NCl", "tpcNClsFound", "tpcCrossedRowsOverFindableCls", "tpcChi2NCl", "tpcFractionSharedCls", "tpcNSigmaPi", "tpcNSigmaKa"}, "input features of the D+ model, in the order it expects"};
    Configurable<double> thresholdScorePionDplus{"thresholdScorePionDplus", 0., "min. D+ pion-class score"};
    Configurable<double> thresholdScoreKaonDplus{"thresholdScoreKaonDplus", 0., "min. D+ kaon-class score"};
    // Ds model
    Configurable<bool> applyMlDs{"applyMlDs", true, "evaluate the Ds track model"};
    Configurable<std::string> onnxFileNameDs{"onnxFileNameDs", "ModelHandler_DsTracks.onnx", "ONNX file name of the Ds track model"};
    Configurable<std::vector<std::string>> inputFeaturesDs{"inputFeaturesDs", std::vector<std::string>{"pt", "eta", "dcaXY", "dcaZ", "sigmaDcaXY", "sigmaDcaZ", "normDcaXY", "normDcaZ", "signed1Pt", "tgl", "sign", "isPvContributor", "itsNCls", "itsNClsInnerBarrel", "itsChi2NCl", "tpcNClsFound", "tpcCrossedRowsOverFindableCls", "tpcChi2NCl", "tpcFractionSharedCls", "tpcNSigmaPi", "tpcNSigmaKa"}, "input features of the Ds model, in the order it expects"};
    Configurable<double> thresholdScorePionDs{"thresholdScorePionDs", 0., "min. Ds pion-class score"};
    Configurable<double> thresholdScoreKaonDs{"thresholdScoreKaonDs", 0., "min. Ds kaon-class score"};
    // ONNX runtime
    Configurable<bool> loadModelsFromCcdb{"loadModelsFromCcdb", false, "load the ONNX models from CCDB instead of a local path"};
    Configurable<std::string> mlModelPathCcdb{"mlModelPathCcdb", "path/to/ml/models", "CCDB path of the ML models"};
    Configurable<int64_t> timestampCcdbForMlModels{"timestampCcdbForMlModels", -1, "timestamp of the ONNX files to be queried in CCDB"};
    Configurable<bool> enableOnnxOptimizations{"enableOnnxOptimizations", true, "enable the ONNX graph optimisations"};
    Configurable<int> onnxThreads{"onnxThreads", 1, "number of threads used by the ONNX runtime (0 = let onnxruntime decide)"};
    // CCDB
    Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
    Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
    Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};
  } config;

  /// One track-level model: the response, its score thresholds and whether it runs.
  struct HfTrackModel {
    HfMlResponseHfTracks<float> response{};
    double thresholdPion{0.};
    double thresholdKaon{0.};
    bool enabled{false};
  };
  std::array<HfTrackModel, NChannels3Prong> models;

  Service<o2::ccdb::BasicCCDBManager> ccdb{};
  o2::ccdb::CcdbApi ccdbApi;
  o2::base::MatLayerCylSet* lut{};
  o2::base::Propagator::MatCorrType noMatCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  int runNumber{};

  // running variables for the re-propagation of tracks to a non-default collision
  o2::dataformats::DCA dcaInfoCov;
  o2::dataformats::VertexBase vtx;

  using TracksWithSelAndPid = soa::Join<aod::TracksWCovDcaExtra, aod::TracksDCACov, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTPCFullKa>;
  using CollisionsWithSel = soa::Join<aod::Collisions, aod::HfSelCollision>;
  using CollisionsWithCentFT0CAndSel = soa::Join<aod::Collisions, aod::CentFT0Cs, aod::HfSelCollision>;

  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;

  // per-data-frame accumulators for hTiming, filled into it at the end of each process call
  std::array<double, NTrackTimingSteps> timing{};
  bool doTiming{false}; ///< cached config.enableTiming: read once, used on the hot path

  HistogramRegistry registry{"registry"};

  /// Add this data frame's accumulated seconds to hTiming and reset the accumulators.
  void flushTiming()
  {
    if (doTiming) {
      // Fill(bin, weight) so the bins accumulate seconds over the whole run
      for (int iStep = 0; iStep < NTrackTimingSteps; iStep++) {
        registry.fill(HIST("hTiming"), iStep, timing[iStep]);
      }
    }
    timing.fill(0.);
  }

  /// Configure one model from its configurables.
  /// \param model is the model to configure
  /// \param onnxFile is the ONNX file name
  /// \param features are the input feature names, in the order the model expects
  /// \param thrPion, thrKaon are the score thresholds
  /// \param name is used in the log message
  void configureModel(HfTrackModel& model,
                      std::string const& onnxFile,
                      std::vector<std::string> const& features,
                      const double thrPion,
                      const double thrKaon,
                      const char* name)
  {
    model.thresholdPion = thrPion;
    model.thresholdKaon = thrKaon;

    const std::vector<double> binsPtSingle{0., PtMaxModel};
    const std::vector<std::string> onnxFiles{onnxFile};
    const auto nOut = static_cast<std::size_t>(HfTrackMlClass::NClasses);
    const std::vector<int> cutDir(nOut, o2::cuts_ml::CutDirection::CutNot);
    const std::vector<double> dummyCutValues(nOut, 0.);
    const o2::framework::LabeledArray<double> dummyCuts{dummyCutValues.data(), 1u, static_cast<uint32_t>(nOut), {}, {}};

    model.response.configure(binsPtSingle, dummyCuts, cutDir, static_cast<uint8_t>(nOut));
    if (config.loadModelsFromCcdb) {
      ccdbApi.init(config.ccdbUrl);
      model.response.setModelPathsCCDB(onnxFiles, ccdbApi, std::vector<std::string>{config.mlModelPathCcdb.value}, config.timestampCcdbForMlModels);
    } else {
      model.response.setModelPathsLocal(onnxFiles);
    }
    model.response.cacheInputFeaturesIndices(features);
    model.response.init(config.enableOnnxOptimizations, config.onnxThreads);
    model.enabled = true;
    LOGP(info, "{}: configured from {} with {} input features and {} output classes (pion score > {}, kaon score > {})",
         name, onnxFile, features.size(), nOut, thrPion, thrKaon);
  }

  void init(InitContext const&)
  {
    if (!config.testAcknowledgement) {
      LOGF(fatal, "ml-based-track-selector is for test only, please do not use for analysis. If you are aware of what you are doing, set testAcknowledgement to true in the configuration.");
    }

    const std::array<bool, 2> doProcess = {doprocessTracks, doprocessTracksWithCentFT0C};
    if (std::count(doProcess.begin(), doProcess.end(), true) != 1) {
      LOGP(fatal, "One and only one process function of HfTrackSelectorTagSelTracks can be enabled at a time!");
    }

    if (config.applyMlDplus) {
      configureModel(models[ChannelDplusToPiKPi], config.onnxFileNameDplus, config.inputFeaturesDplus,
                     config.thresholdScorePionDplus, config.thresholdScoreKaonDplus, "D+ track model");
    }
    if (config.applyMlDs) {
      configureModel(models[ChannelDsToKKPi], config.onnxFileNameDs, config.inputFeaturesDs,
                     config.thresholdScorePionDs, config.thresholdScoreKaonDs, "Ds track model");
    }
    if (!models[ChannelDplusToPiKPi].enabled && !models[ChannelDsToKKPi].enabled) {
      LOGP(fatal, "At least one of the D+ and Ds track models must be enabled!");
    }

    ccdb->setURL(config.ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(config.ccdbPathLut));
    runNumber = 0;

    if (config.fillHistograms) {
      const AxisSpec axisPtProng{360, 0., 36., "#it{p}_{T}^{track} (GeV/#it{c})"};
      const AxisSpec axisScore{100, 0., 1., "ML score"};
      const AxisSpec axisEta{100, -1., 1., "#it{#eta}"};

      registry.add("hPtNoCuts", "all track associations;#it{p}_{T}^{track} (GeV/#it{c});entries", {HistType::kTH1D, {axisPtProng}});
      registry.add("hPtQuality", "track associations passing the pT floor and the quality cut, scored;#it{p}_{T}^{track} (GeV/#it{c});entries", {HistType::kTH1D, {axisPtProng}});
      registry.add("hPtQualityRejColl", "track associations passing the pT floor and the quality cut, not scored: collision fails the HF event selection;#it{p}_{T}^{track} (GeV/#it{c});entries", {HistType::kTH1D, {axisPtProng}});
      registry.add("hPtSelected", "track associations selected by at least one model;#it{p}_{T}^{track} (GeV/#it{c});entries", {HistType::kTH1D, {axisPtProng}});
      registry.add("hEtaSelected", "track associations selected by at least one model;#it{#eta};entries", {HistType::kTH1D, {axisEta}});
      registry.add("hScorePionDplus", "D^{#plus} pion-class score;#it{p}_{T}^{track} (GeV/#it{c});score;entries", {HistType::kTH2D, {axisPtProng, axisScore}});
      registry.add("hScoreKaonDplus", "D^{#plus} kaon-class score;#it{p}_{T}^{track} (GeV/#it{c});score;entries", {HistType::kTH2D, {axisPtProng, axisScore}});
      registry.add("hScorePionDs", "D_{s}^{#plus} pion-class score;#it{p}_{T}^{track} (GeV/#it{c});score;entries", {HistType::kTH2D, {axisPtProng, axisScore}});
      registry.add("hScoreKaonDs", "D_{s}^{#plus} kaon-class score;#it{p}_{T}^{track} (GeV/#it{c});score;entries", {HistType::kTH2D, {axisPtProng, axisScore}});
      // seconds, accumulated over the run via Fill(bin, weight); per-call cost = bin / hPtQuality entries
      auto hTiming = registry.add<TH1>("hTiming", "CPU of the ML track selection;;seconds", {HistType::kTH1D, {{NTrackTimingSteps, -0.5f, static_cast<float>(NTrackTimingSteps) - 0.5f}}});
      hTiming->GetXaxis()->SetBinLabel(TrackTimeTotal + 1, "track loop total");
      hTiming->GetXaxis()->SetBinLabel(TrackTimeFeatures + 1, "features");
      hTiming->GetXaxis()->SetBinLabel(TrackTimeModelDplus + 1, "D+ model");
      hTiming->GetXaxis()->SetBinLabel(TrackTimeModelDs + 1, "Ds model");
    }
    doTiming = config.enableTiming && config.fillHistograms;
  }

  /// Fill the container of ML input features for one (track, collision) association.
  /// The DCA and its resolution are recomputed by propagating to the collision under study
  /// whenever that is not the track's own collision.
  /// \param collision is the collision the track is associated to
  /// \param track is the track
  /// \param centrality is the centrality percentile of the collision
  /// \param features is the container to be filled
  template <typename TCollision, typename TTrack>
  void fillFeatures(TCollision const& collision, TTrack const& track, float centrality, HfTrackMlFeatures& features)
  {
    float trackPt = track.pt();
    float trackEta = track.eta();
    float dcaXY = track.dcaXY();
    float dcaZ = track.dcaZ();
    float sigmaDcaXY = std::sqrt(std::abs(track.sigmaDcaXY2()));
    float sigmaDcaZ = std::sqrt(std::abs(track.sigmaDcaZ2()));

    if (track.collisionId() != collision.globalIndex()) {
      // not the default collision of this track: re-propagate to get DCA and its covariance
      auto trackParCov = getTrackParCov(track);
      dcaInfoCov.set(999.f, 999.f, 999.f, 999.f, 999.f);
      vtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
      vtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
      if (o2::base::Propagator::Instance()->propagateToDCABxByBz(vtx, trackParCov, 2.f, noMatCorr, &dcaInfoCov)) {
        trackPt = trackParCov.getPt();
        trackEta = trackParCov.getEta();
        dcaXY = dcaInfoCov.getY();
        dcaZ = dcaInfoCov.getZ();
        sigmaDcaXY = std::sqrt(std::abs(dcaInfoCov.getSigmaY2()));
        sigmaDcaZ = std::sqrt(std::abs(dcaInfoCov.getSigmaZ2()));
      }
    }

    constexpr float MinSigmaDca = 1.e-6f; // guard against a null resolution before dividing

    features.pt = trackPt;
    features.eta = trackEta;
    features.dcaXY = dcaXY;
    features.dcaZ = dcaZ;
    features.sigmaDcaXY = sigmaDcaXY;
    features.sigmaDcaZ = sigmaDcaZ;
    features.normDcaXY = dcaXY / std::max(sigmaDcaXY, MinSigmaDca);
    features.normDcaZ = dcaZ / std::max(sigmaDcaZ, MinSigmaDca);
    features.signed1Pt = track.signed1Pt();
    features.tgl = track.tgl();
    features.sign = static_cast<float>(track.sign());
    features.isPvContributor = static_cast<float>(track.isPVContributor());
    features.itsNCls = static_cast<float>(track.itsNCls());
    features.itsNClsInnerBarrel = static_cast<float>(track.itsNClsInnerBarrel());
    features.itsChi2NCl = track.itsChi2NCl();
    features.tpcNClsFound = static_cast<float>(track.tpcNClsFound());
    features.tpcCrossedRowsOverFindableCls = track.tpcCrossedRowsOverFindableCls();
    features.tpcChi2NCl = track.tpcChi2NCl();
    features.tpcFractionSharedCls = track.tpcFractionSharedCls();
    features.tpcNSigmaPi = track.tpcNSigmaPi();
    features.tpcNSigmaKa = track.tpcNSigmaKa();
    features.centrality = centrality;
  }

  /// Evaluate one model and turn its output into the two role decisions.
  /// \param model is the model to evaluate
  /// \param features are the track features
  /// \param scores is filled with (background, pion-role, kaon-role); -1 if not evaluated
  /// \param isSelPion, isSelKaon are filled with the threshold decisions
  /// \return true if the model was actually evaluated
  bool evaluateModel(HfTrackModel& model, HfTrackMlFeatures const& features,
                     std::array<float, 3>& scores, bool& isSelPion, bool& isSelKaon)
  {
    scores = {-1.f, -1.f, -1.f};
    isSelPion = false;
    isSelKaon = false;
    if (!model.enabled) {
      return false;
    }
    auto inputFeatures = model.response.getInputFeatures(features);
    // one model for the whole pT range, so the bin index is always 0
    const auto output = model.response.getModelOutput(inputFeatures, 0);
    if (output.size() != scores.size()) {
      LOGP(fatal, "the model emitted {} scores, expected {}: only 3-class (bkg/pi/K) models are supported", output.size(), scores.size());
    }
    std::copy(output.begin(), output.end(), scores.begin());

    isSelPion = scores[1] > model.thresholdPion;
    isSelKaon = scores[2] > model.thresholdKaon;
    return true;
  }

  /// Selection tag for tracks
  /// \param collision is the collision iterator
  /// \param trackIndicesCollision are the track indices associated to this collision
  /// \param centrality is the centrality percentile of the collision
  /// \param isCollisionSelected is whether the collision passes the HF event selection
  template <typename TTracks, typename TCollision, typename GroupedTrackIndices>
  void runTagSelTracks(TCollision const& collision,
                       TTracks const& /*tracks*/,
                       GroupedTrackIndices const& trackIndicesCollision,
                       float centrality,
                       const bool isCollisionSelected)
  {
    const bool scoreCollision = isCollisionSelected;

    for (const auto& trackId : trackIndicesCollision) {
      const auto track = trackId.template track_as<TTracks>();

      // One row per track association, always: aod::HfSelTrack is joined with aod::TrackAssoc
      // downstream, so the two tables have to stay index-aligned. Rejected associations are
      // written with an empty selection mask, not skipped.
      uint32_t statusProng{0};
      uint32_t isIdentifiedPid{0};

      if (config.fillHistograms) {
        registry.fill(HIST("hPtNoCuts"), track.pt());
      }

      const bool passesQuality = track.pt() >= config.ptMinTrack && track.isGlobalTrackWoDCA();
      if (passesQuality && !scoreCollision && config.fillHistograms) {
        registry.fill(HIST("hPtQualityRejColl"), track.pt());
      }

      if (passesQuality && scoreCollision) {
        if (config.fillHistograms) {
          registry.fill(HIST("hPtQuality"), track.pt());
        }

        HfTrackMlFeatures features{};
        auto tStart = tickIf(doTiming);
        fillFeatures(collision, track, centrality, features);
        if (doTiming) {
          timing[TrackTimeFeatures] += elapsedSeconds(tStart);
        }

        std::array<float, 3> scores{};
        bool isSelPion{false};
        bool isSelKaon{false};

        tStart = tickIf(doTiming);
        const bool evaluatedDplus = evaluateModel(models[ChannelDplusToPiKPi], features, scores, isSelPion, isSelKaon);
        if (doTiming) {
          timing[TrackTimeModelDplus] += elapsedSeconds(tStart);
        }
        if (evaluatedDplus) {
          if (isSelPion) {
            SETBIT(isIdentifiedPid, RolePionDplus);
          }
          if (isSelKaon) {
            SETBIT(isIdentifiedPid, RoleKaonDplus);
          }
          if (config.fillHistograms) {
            registry.fill(HIST("hScorePionDplus"), features.pt, scores[1]);
            registry.fill(HIST("hScoreKaonDplus"), features.pt, scores[2]);
          }
        }

        tStart = tickIf(doTiming);
        const bool evaluatedDs = evaluateModel(models[ChannelDsToKKPi], features, scores, isSelPion, isSelKaon);
        if (doTiming) {
          timing[TrackTimeModelDs] += elapsedSeconds(tStart);
        }
        if (evaluatedDs) {
          if (isSelPion) {
            SETBIT(isIdentifiedPid, RolePionDs);
          }
          if (isSelKaon) {
            SETBIT(isIdentifiedPid, RoleKaonDs);
          }
          if (config.fillHistograms) {
            registry.fill(HIST("hScorePionDs"), features.pt, scores[1]);
            registry.fill(HIST("hScoreKaonDs"), features.pt, scores[2]);
          }
        }

        // the track enters the combinatorics if it can play any role in any of the two channels
        if (isIdentifiedPid != 0u) {
          SETBIT(statusProng, CandidateType::Cand3Prong);
          if (config.fillHistograms) {
            registry.fill(HIST("hPtSelected"), features.pt);
            registry.fill(HIST("hEtaSelected"), features.eta);
          }
        }
      }

      const bool isPositive = track.sign() > 0;
      rowSelectedTrack(statusProng, isIdentifiedPid, isPositive);
    }
  }

  /// Set the magnetic field and the material LUT for the current collision, needed by the
  /// re-propagation of tracks associated to a collision that is not their own.
  template <typename TCollision>
  void initMagneticField(TCollision const& collision)
  {
    const auto bc = collision.template bc_as<o2::aod::BCsWithTimestamps>();
    initCCDB(bc, runNumber, ccdb, config.ccdbPathGrpMag, lut, false);
  }

  void processTracks(CollisionsWithSel const& collisions,
                     aod::TrackAssoc const& trackIndices,
                     TracksWithSelAndPid const& tracks,
                     aod::BCsWithTimestamps const&)
  {
    rowSelectedTrack.reserve(trackIndices.size());
    const auto tStart = tickIf(doTiming);
    for (const auto& collision : collisions) {
      initMagneticField(collision);
      const auto groupedTrackIndices = trackIndices.sliceBy(trackIndicesPerCollision, collision.globalIndex());
      runTagSelTracks<TracksWithSelAndPid>(collision, tracks, groupedTrackIndices, -1.f, collision.whyRejectColl() == 0);
    }
    if (doTiming) {
      timing[TrackTimeTotal] += elapsedSeconds(tStart);
    }
    flushTiming();
  }
  PROCESS_SWITCH(HfTrackSelectorTagSelTracks, processTracks, "Process tracks without centrality", true);

  void processTracksWithCentFT0C(CollisionsWithCentFT0CAndSel const& collisions,
                                 aod::TrackAssoc const& trackIndices,
                                 TracksWithSelAndPid const& tracks,
                                 aod::BCsWithTimestamps const&)
  {
    rowSelectedTrack.reserve(trackIndices.size());
    const auto tStart = tickIf(doTiming);
    for (const auto& collision : collisions) {
      initMagneticField(collision);
      const auto groupedTrackIndices = trackIndices.sliceBy(trackIndicesPerCollision, collision.globalIndex());
      runTagSelTracks<TracksWithSelAndPid>(collision, tracks, groupedTrackIndices, collision.centFT0C(), collision.whyRejectColl() == 0);
    }
    if (doTiming) {
      timing[TrackTimeTotal] += elapsedSeconds(tStart);
    }
    flushTiming();
  }
  PROCESS_SWITCH(HfTrackSelectorTagSelTracks, processTracksWithCentFT0C, "Process tracks with FT0C centrality as ML input feature", false);
};

/// Pre-selection of 3-prong secondary vertices
struct HfMlBasedTrackSelector {
  Produces<aod::Hf3Prongs> rowTrackIndexProng3;

  struct : ConfigurableGroup {
    Configurable<bool> fillHistograms{"fillHistograms", true, "fill histograms"};
    // preselection
    Configurable<double> ptTolerance{"ptTolerance", 0.1, "pT tolerance in GeV/c for applying preselections before vertex reconstruction"};
    // preselection of 3-prongs using the decay length computed only with the first two tracks
    Configurable<double> minTwoTrackDecayLengthFor3Prongs{"minTwoTrackDecayLengthFor3Prongs", 0., "Minimum decay length computed with 2 tracks for 3-prongs to speedup combinatorial"};
    Configurable<double> maxTwoTrackChi2PcaFor3Prongs{"maxTwoTrackChi2PcaFor3Prongs", 1.e10, "Maximum chi2 pca computed with 2 tracks for 3-prongs to speedup combinatorial"};
    Configurable<bool> enableTiming{"enableTiming", false, "fill hTiming/hLoopCounters with the CPU and rates of the triple loop (adds two clock reads per fit)"};
    // vertexing
    Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
    Configurable<bool> useAbsDCA{"useAbsDCA", false, "Minimise abs. distance rather than chi2"};
    Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
    Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
    Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
    Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
    Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations if chi2/chi2old > this"};
    // CCDB
    Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
    Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
    Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};

    // D+ cuts
    Configurable<std::vector<double>> binsPtDplusToPiKPi{"binsPtDplusToPiKPi", std::vector<double>{hf_cuts_presel_3prong::vecBinsPt}, "pT bin limits for D+->piKpi pT-dependent cuts"};
    Configurable<LabeledArray<double>> cutsDplusToPiKPi{"cutsDplusToPiKPi", {hf_cuts_presel_3prong::Cuts[0], hf_cuts_presel_3prong::NBinsPt, hf_cuts_presel_3prong::NCutVars, hf_cuts_presel_3prong::labelsPt, hf_cuts_presel_3prong::labelsCutVar}, "D+->piKpi selections per pT bin"};
    // Ds+ cuts
    Configurable<std::vector<double>> binsPtDsToKKPi{"binsPtDsToKKPi", std::vector<double>{hf_cuts_presel_ds::vecBinsPt}, "pT bin limits for Ds+->KKPi pT-dependent cuts"};
    Configurable<LabeledArray<double>> cutsDsToKKPi{"cutsDsToKKPi", {hf_cuts_presel_ds::Cuts[0], hf_cuts_presel_ds::NBinsPt, hf_cuts_presel_ds::NCutVars, hf_cuts_presel_ds::labelsPt, hf_cuts_presel_ds::labelsCutVar}, "Ds+->KKPi selections per pT bin"};
  } config;

  SliceCache cache;
  o2::vertexing::DCAFitterN<2> df2; // 2-prong vertex fitter, only used for the 3-prong speed-up
  o2::vertexing::DCAFitterN<3> df3; // 3-prong vertex fitter
  Service<o2::ccdb::BasicCCDBManager> ccdb{};
  o2::base::MatLayerCylSet* lut{};
  o2::base::Propagator::MatCorrType noMatCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  int runNumber{};

  // masses of the two mass hypotheses of each channel, in the (same-sign, opposite-sign,
  // same-sign) prong ordering used by the combinatorics
  std::array<std::array<std::array<double, NProngs>, NMassHypos>, NChannels3Prong> arrMass3Prong{};
  // cuts and pT binning, one entry per channel
  std::array<LabeledArray<double>, NChannels3Prong> cut3Prong{};
  std::array<std::vector<double>, NChannels3Prong> binsPt3Prong{};

  /// One track of the collision under study, already propagated to that collision's PV.
  struct HfProngCandidate {
    o2::track::TrackParCov parCov;
    std::array<float, 3> pVec{};
    uint32_t mask{};       ///< aod::HfSelTrack::isIdentifiedPid
    int64_t globalIndex{}; ///< index in the track table, written to the skim
  };

  /// Cache of propagated tracks for one collision, split by charge.
  struct HfProngCache {
    std::vector<HfProngCandidate> prongs;
    std::vector<int> sameSign; ///< can sit in a same-sign slot: pi(D+), pi(Ds) or K(Ds)
    std::vector<int> oppSign;  ///< can sit in the opposite-sign slot: K(D+) or K(Ds)

    void clear()
    {
      prongs.clear();
      sameSign.clear();
      oppSign.clear();
    }
  };
  HfProngCache cachePos;
  HfProngCache cacheNeg;

  // per-collision accumulators for hTiming / hLoopCounters, reset at each collision
  std::array<double, NTimingSteps> timing{};
  std::array<int64_t, NLoopCounters> counters{};
  bool doTiming{false}; ///< cached config.enableTiming: read once, used on the hot path

  using SelectedCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::HfSelCollision>>;
  using FilteredTrackAssocSel = soa::Filtered<soa::Join<aod::TrackAssoc, aod::HfSelTrack>>;

  // filter collisions
  Filter filterSelectCollisions = (aod::hf_sel_collision::whyRejectColl == static_cast<o2::hf_evsel::HfCollisionRejectionMask>(0));
  // filter track indices
  Filter filterSelectTrackIds = ((aod::hf_sel_track::isSelProng & static_cast<uint32_t>(BIT(CandidateType::Cand3Prong))) != 0u);

  Preslice<FilteredTrackAssocSel> trackIndicesPerCollision = aod::track_association::collisionId;

  Partition<FilteredTrackAssocSel> positiveFor3Prongs = aod::hf_sel_track::isPositive == true && ((aod::hf_sel_track::isSelProng & static_cast<uint32_t>(BIT(CandidateType::Cand3Prong))) != 0u);
  Partition<FilteredTrackAssocSel> negativeFor3Prongs = aod::hf_sel_track::isPositive == false && ((aod::hf_sel_track::isSelProng & static_cast<uint32_t>(BIT(CandidateType::Cand3Prong))) != 0u);

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    if (!doprocess3Prongs) {
      return;
    }
    doTiming = config.enableTiming && config.fillHistograms;

    arrMass3Prong[ChannelDplusToPiKPi] = std::array{std::array{MassPiPlus, MassKPlus, MassPiPlus},
                                                    std::array{MassPiPlus, MassKPlus, MassPiPlus}};

    arrMass3Prong[ChannelDsToKKPi] = std::array{std::array{MassKPlus, MassKPlus, MassPiPlus},
                                                std::array{MassPiPlus, MassKPlus, MassKPlus}};

    // cuts retrieved by json, in the order of Channel3Prong
    cut3Prong = {config.cutsDplusToPiKPi, config.cutsDsToKKPi};
    binsPt3Prong = {config.binsPtDplusToPiKPi, config.binsPtDsToKKPi};

    df2.setPropagateToPCA(config.propagateToPCA);
    df2.setMaxR(config.maxR);
    df2.setMaxDZIni(config.maxDZIni);
    df2.setMinParamChange(config.minParamChange);
    df2.setMinRelChi2Change(config.minRelChi2Change);
    df2.setUseAbsDCA(config.useAbsDCA);
    df2.setWeightedFinalPCA(config.useWeightedFinalPCA);

    df3.setPropagateToPCA(config.propagateToPCA);
    df3.setMaxR(config.maxR);
    df3.setMaxDZIni(config.maxDZIni);
    df3.setMinParamChange(config.minParamChange);
    df3.setMinRelChi2Change(config.minRelChi2Change);
    df3.setUseAbsDCA(config.useAbsDCA);
    df3.setWeightedFinalPCA(config.useWeightedFinalPCA);

    ccdb->setURL(config.ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(config.ccdbPathLut));
    runNumber = 0;

    if (config.fillHistograms) {
      const AxisSpec axisNumTracks{500, -0.5f, 499.5f, "Number of tracks"};
      const AxisSpec axisNumCands{1000, -0.5f, 999.5f, "Number of candidates"};
      registry.add("hNTracks", "Number of selected tracks;# of selected tracks;entries", {HistType::kTH1D, {axisNumTracks}});
      // seconds and raw counts, accumulated over the run via Fill(bin, weight)
      auto hTiming = registry.add<TH1>("hTiming", "CPU inside the triple loop;;seconds", {HistType::kTH1D, {{NTimingSteps, -0.5f, static_cast<float>(NTimingSteps) - 0.5f}}});
      hTiming->GetXaxis()->SetBinLabel(TimeLoopTotal + 1, "triple loop total");
      hTiming->GetXaxis()->SetBinLabel(TimeFit2Prong + 1, "2-prong vertex fit");
      hTiming->GetXaxis()->SetBinLabel(TimeFit3Prong + 1, "3-prong vertex fit");
      hTiming->GetXaxis()->SetBinLabel(TimeCacheFill + 1, "prong cache fill");
      hTiming->GetXaxis()->SetBinLabel(TimeCachePropagate + 1, "cache re-propagation");
      auto hLoopCounters = registry.add<TH1>("hLoopCounters", "Triple loop stages;;entries", {HistType::kTH1D, {{NLoopCounters, -0.5f, static_cast<float>(NLoopCounters) - 0.5f}}});
      hLoopCounters->GetXaxis()->SetBinLabel(CountPairsSeen + 1, "pairs seen");
      hLoopCounters->GetXaxis()->SetBinLabel(CountPairsCandidateOk + 1, "pairs passing ML roles");
      hLoopCounters->GetXaxis()->SetBinLabel(CountFit2Prong + 1, "2-prong fits");
      hLoopCounters->GetXaxis()->SetBinLabel(CountPairsRejected2P + 1, "pairs rejected by 2-prong");
      hLoopCounters->GetXaxis()->SetBinLabel(CountTripletsEnumerated + 1, "triplets enumerated");
      hLoopCounters->GetXaxis()->SetBinLabel(CountTripletsWritten + 1, "triplets written");
      hLoopCounters->GetXaxis()->SetBinLabel(CountTracksCached + 1, "tracks cached");
      hLoopCounters->GetXaxis()->SetBinLabel(CountTracksPropagated + 1, "tracks re-propagated");
      registry.add("hVtx3ProngX", "3-prong candidates;#it{x}_{sec. vtx.} (cm);entries", {HistType::kTH1D, {{1000, -2., 2.}}});
      registry.add("hVtx3ProngY", "3-prong candidates;#it{y}_{sec. vtx.} (cm);entries", {HistType::kTH1D, {{1000, -2., 2.}}});
      registry.add("hVtx3ProngZ", "3-prong candidates;#it{z}_{sec. vtx.} (cm);entries", {HistType::kTH1D, {{1000, -20., 20.}}});
      registry.add("hNCand3Prong", "3-prong candidates preselected;# of candidates;entries", {HistType::kTH1D, {axisNumCands}});
      registry.add("hNCand3ProngVsNTracks", "3-prong candidates preselected;# of selected tracks;# of candidates;entries", {HistType::kTH2D, {axisNumTracks, axisNumCands}});
      registry.add("hMassDPlusToPiKPi", "D^{#plus} candidates;inv. mass (#pi K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{500, 0., 5.}}});
      registry.add("hMassDsToKKPi", "D_{s}^{#plus} candidates;inv. mass (K K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{500, 0., 5.}}});
    }
  }

  /// Check whether the two prongs can satisfy any mass hypothesis of either channel, based on their ML role masks.
  /// \param maskOpp is the isIdentifiedPid mask of the opposite-sign prong
  /// \param maskSame is the isIdentifiedPid mask of the first same-sign prong
  static bool isCandidateAllowed(const uint32_t maskOpp, const uint32_t maskSame)
  {
    const bool okDplus = TESTBIT(maskOpp, RoleKaonDplus) && TESTBIT(maskSame, RolePionDplus);
    const bool okDs = TESTBIT(maskOpp, RoleKaonDs) && (TESTBIT(maskSame, RoleKaonDs) || TESTBIT(maskSame, RolePionDs));
    return okDplus || okDs;
  }

  /// Check whether the three ML masks can satisfy any mass hypothesis of either channel
  /// \param maskOpp is the mask of the opposite-sign prong (prong 1)
  /// \param mask1, mask2 are the masks of the same-sign prongs (prongs 0 and 2)
  static bool tripletRoleOk(const uint32_t maskOpp, const uint32_t mask1, const uint32_t mask2)
  {
    const std::array<uint32_t, NProngs> masks{mask1, maskOpp, mask2};
    for (int iChannel = 0; iChannel < NChannels3Prong; iChannel++) {
      for (int iHypo = 0; iHypo < NMassHypos; iHypo++) {
        bool isHypoAlive = true;
        for (int iProng = 0; iProng < NProngs; iProng++) {
          const int role = IsKaonProng[iChannel][iHypo][iProng] ? RoleKaonOfChannel[iChannel] : RolePionOfChannel[iChannel];
          if (!TESTBIT(masks[iProng], role)) {
            isHypoAlive = false;
            break;
          }
        }
        if (isHypoAlive) {
          return true;
        }
      }
    }
    return false;
  }

  /// Require the ML role of each prong to match the mass hypothesis under test.
  /// The opposite-sign prong is the kaon in both channels; the same-sign pair is (pi, pi) for the
  /// D+ and (K, pi) or (pi, K) for the Ds, which is exactly what the two mass hypotheses encode.
  /// \param isIdentifiedPid0,1,2 are the aod::HfSelTrack::isIdentifiedPid masks of the prongs, in
  ///        the (same-sign, opposite-sign, same-sign) ordering
  /// \param whichHypo information of the mass hypotheses that are still alive
  /// \param isSelected is a bitmap with selection outcome
  void applyMlRoleSelection(const uint32_t isIdentifiedPid0, const uint32_t isIdentifiedPid1, const uint32_t isIdentifiedPid2,
                            std::array<int, NChannels3Prong>& whichHypo, auto& isSelected)
  {
    const std::array<uint32_t, NProngs> isIdentifiedPid{isIdentifiedPid0, isIdentifiedPid1, isIdentifiedPid2};

    for (int iChannel = 0; iChannel < NChannels3Prong; iChannel++) {
      for (int iHypo = 0; iHypo < NMassHypos; iHypo++) {
        for (int iProng = 0; iProng < NProngs; iProng++) {
          const int role = IsKaonProng[iChannel][iHypo][iProng] ? RoleKaonOfChannel[iChannel] : RolePionOfChannel[iChannel];
          if (!TESTBIT(isIdentifiedPid[iProng], role)) {
            CLRBIT(whichHypo[iChannel], iHypo);
            break;
          }
        }
      }
      if (whichHypo[iChannel] == 0) {
        CLRBIT(isSelected, iChannel);
      }
    }
  }

  /// Method to perform selections on difference from nominal mass for phi decay
  /// \param binPt pt bin for the cuts
  /// \param pVecTrack0 is the momentum array of the first daughter track
  /// \param pVecTrack1 is the momentum array of the second daughter track
  /// \param pVecTrack2 is the momentum array of the third daughter track
  /// \param whichHypo information of the mass hypoteses that were selected
  /// \param isSelected is a bitmap with selection outcome
  template <typename T1>
  void applyPreselectionPhiDecay(const int binPt, T1 const& pVecTrack0, T1 const& pVecTrack1, T1 const& pVecTrack2,
                                 std::array<int, NChannels3Prong>& whichHypo, auto& isSelected)
  {
    const double deltaMassMax = cut3Prong[ChannelDsToKKPi].get(binPt, 5u); // 5u == "deltaMassKK"
    if (TESTBIT(whichHypo[ChannelDsToKKPi], 0)) {
      const double mass2PhiKKPi = RecoDecay::m2(std::array{pVecTrack0, pVecTrack1}, std::array{arrMass3Prong[ChannelDsToKKPi][0][0], arrMass3Prong[ChannelDsToKKPi][0][1]});
      if (mass2PhiKKPi > (MassPhi + deltaMassMax) * (MassPhi + deltaMassMax) || (deltaMassMax < MassPhi && mass2PhiKKPi < (MassPhi - deltaMassMax) * (MassPhi - deltaMassMax))) {
        CLRBIT(whichHypo[ChannelDsToKKPi], 0);
      }
    }
    if (TESTBIT(whichHypo[ChannelDsToKKPi], 1)) {
      const double mass2PhiPiKK = RecoDecay::m2(std::array{pVecTrack1, pVecTrack2}, std::array{arrMass3Prong[ChannelDsToKKPi][1][1], arrMass3Prong[ChannelDsToKKPi][1][2]});
      if (mass2PhiPiKK > (MassPhi + deltaMassMax) * (MassPhi + deltaMassMax) || (deltaMassMax < MassPhi && mass2PhiPiKK < (MassPhi - deltaMassMax) * (MassPhi - deltaMassMax))) {
        CLRBIT(whichHypo[ChannelDsToKKPi], 1);
      }
    }
    if (whichHypo[ChannelDsToKKPi] == 0) {
      CLRBIT(isSelected, ChannelDsToKKPi);
    }
  }

  /// Method to perform selections for 3-prong candidates before vertex reconstruction
  /// \param pVecTrack0 is the momentum array of the first daughter track
  /// \param pVecTrack1 is the momentum array of the second daughter track
  /// \param pVecTrack2 is the momentum array of the third daughter track
  /// \param whichHypo information of the mass hypotheses that were selected
  /// \param isSelected is a bitmap with selection outcome
  template <typename T2>
  void applyPreselection3Prong(T2 const& pVecTrack0, T2 const& pVecTrack1, T2 const& pVecTrack2,
                               std::array<int, NChannels3Prong>& whichHypo, auto& isSelected)
  {
    const auto pt = RecoDecay::pt(pVecTrack0, pVecTrack1, pVecTrack2) + config.ptTolerance; // add tolerance because of no reco decay vertex

    for (int iChannel = 0; iChannel < NChannels3Prong; iChannel++) {
      if (!TESTBIT(isSelected, iChannel)) {
        whichHypo[iChannel] = 0;
        continue;
      }

      // pT
      const auto binPt = findBin(&binsPt3Prong[iChannel], pt);
      // return immediately if it is outside the defined pT bins
      if (binPt == -1) {
        CLRBIT(isSelected, iChannel);
        whichHypo[iChannel] = 0;
        continue;
      }

      // invariant mass
      const double minMass = cut3Prong[iChannel].get(binPt, 0u);
      const double maxMass = cut3Prong[iChannel].get(binPt, 1u);
      if (minMass >= 0. && maxMass > 0.) {
        std::array<double, 2> massHypos = {0., 0.};
        const std::array arrMom{pVecTrack0, pVecTrack1, pVecTrack2};
        const double min2 = minMass * minMass;
        const double max2 = maxMass * maxMass;
        massHypos[0] = RecoDecay::m2(arrMom, arrMass3Prong[iChannel][0]);
        massHypos[1] = (iChannel != ChannelDplusToPiKPi) ? RecoDecay::m2(arrMom, arrMass3Prong[iChannel][1]) : massHypos[0];
        if (massHypos[0] < min2 || massHypos[0] >= max2) {
          CLRBIT(whichHypo[iChannel], 0);
        }
        if (massHypos[1] < min2 || massHypos[1] >= max2) {
          CLRBIT(whichHypo[iChannel], 1);
        }
        if (whichHypo[iChannel] == 0) {
          CLRBIT(isSelected, iChannel);
          continue;
        }
      }

      // prong pT
      const auto ptProngMin = cut3Prong[iChannel].get(binPt, 4u); // 4u == "ptProngMin"
      const auto pt2ProngMin = ptProngMin * ptProngMin;
      if (RecoDecay::pt2(pVecTrack0) < pt2ProngMin || RecoDecay::pt2(pVecTrack1) < pt2ProngMin || RecoDecay::pt2(pVecTrack2) < pt2ProngMin) {
        CLRBIT(isSelected, iChannel);
        continue;
      }

      if (iChannel == ChannelDsToKKPi) {
        applyPreselectionPhiDecay(binPt, pVecTrack0, pVecTrack1, pVecTrack2, whichHypo, isSelected);
      }
    }
  }

  /// Method to perform selections for a 2-track vertex used to speed up the 3-prong combinatorics
  /// \param secVtx is the secondary vertex
  /// \param primVtx is the primary vertex
  /// \param dcaFitter is the DCAFitter used for the 2-track vertex
  /// \returns true if the candidate is selected
  template <typename T1, typename T2, typename T3>
  bool isTwoTrackVertexSelectedFor3Prongs(const T1& secVtx, const T2& primVtx, const T3& dcaFitter)
  {
    if (dcaFitter.getChi2AtPCACandidate() > config.maxTwoTrackChi2PcaFor3Prongs) {
      return false;
    }
    const auto decLen = RecoDecay::distance(primVtx, secVtx);
    return static_cast<bool>(decLen >= config.minTwoTrackDecayLengthFor3Prongs);
  }

  /// Method to perform selections for 3-prong candidates after vertex reconstruction
  /// \param pVecCand is the array for the candidate momentum after reconstruction of secondary vertex
  /// \param secVtx is the secondary vertex
  /// \param primVtx is the primary vertex
  /// \param isSelected is a bitmap with selection outcome
  template <typename T1, typename T2, typename T3>
  void applySelection3Prong(const T1& pVecCand, const T2& secVtx, const T3& primVtx, auto& isSelected)
  {
    if (isSelected == 0) {
      return;
    }

    const auto pt = RecoDecay::pt(pVecCand);
    for (int iChannel = 0; iChannel < NChannels3Prong; iChannel++) {
      if (!TESTBIT(isSelected, iChannel)) {
        continue;
      }

      // pT
      const auto binPt = findBin(&binsPt3Prong[iChannel], pt);
      if (binPt == -1) { // cut if it is outside the defined pT bins
        CLRBIT(isSelected, iChannel);
        continue;
      }

      // cos of pointing angle
      const auto cpa = RecoDecay::cpa(primVtx, secVtx, pVecCand);
      if (cpa < cut3Prong[iChannel].get(binPt, 2u)) { // 2u == "cosp"
        CLRBIT(isSelected, iChannel);
        continue;
      }

      // decay length
      const auto decayLength = RecoDecay::distance(primVtx, secVtx);
      if (decayLength < cut3Prong[iChannel].get(binPt, 3u)) { // 3u == "decL"
        CLRBIT(isSelected, iChannel);
      }
    }
  }

  /// Translate the compact per-channel selection bitmap into the hfflag written to aod::Hf3Prongs,
  /// which uses the hf_cand_3prong::DecayType bit positions.
  uint8_t hfFlagOfSelection(const uint32_t isSelected) const
  {
    uint8_t hfFlag{0};
    for (int iChannel = 0; iChannel < NChannels3Prong; iChannel++) {
      if (TESTBIT(isSelected, iChannel)) {
        SETBIT(hfFlag, DecayTypeOfChannel[iChannel]);
      }
    }
    return hfFlag;
  }

  /// Reconstruct and preselect one 3-prong candidate, and fill its table row if it survives.
  /// \param collision is the collision the candidate belongs to
  /// \param trackParVar0,1,2 are the track parametrisations, in the (same-sign, opposite-sign,
  ///        same-sign) prong ordering
  /// \param pVecTrack0,1,2 are the corresponding momenta at the primary vertex
  /// \param isIdentifiedPid0,1,2 are the corresponding ML role masks
  /// \param globalIndex0,1,2 are the corresponding track global indices
  template <typename TTrackParVar>
  void processTriplet(SelectedCollisions::iterator const& collision,
                      TTrackParVar const& trackParVar0, TTrackParVar const& trackParVar1, TTrackParVar const& trackParVar2,
                      std::array<float, 3> const& pVecTrack0, std::array<float, 3> const& pVecTrack1, std::array<float, 3> const& pVecTrack2,
                      const uint32_t isIdentifiedPid0, const uint32_t isIdentifiedPid1, const uint32_t isIdentifiedPid2,
                      const int64_t globalIndex0, const int64_t globalIndex1, const int64_t globalIndex2)
  {
    uint32_t isSelected3ProngCand = BIT(NChannels3Prong) - 1;
    std::array<int, NChannels3Prong> whichHypo3Prong{};
    whichHypo3Prong.fill(BIT(NMassHypos) - 1); // all mass hypotheses alive

    applyMlRoleSelection(isIdentifiedPid0, isIdentifiedPid1, isIdentifiedPid2, whichHypo3Prong, isSelected3ProngCand);
    if (isSelected3ProngCand == 0) {
      return;
    }

    applyPreselection3Prong(pVecTrack0, pVecTrack1, pVecTrack2, whichHypo3Prong, isSelected3ProngCand);
    if (isSelected3ProngCand == 0) {
      return;
    }

    // reconstruct the 3-prong secondary vertex
    const auto tStart3Prong = tickIf(doTiming);
    int nVtxFrom3ProngFitter = 0;
    try {
      nVtxFrom3ProngFitter = df3.process(trackParVar0, trackParVar1, trackParVar2);
    } catch (...) {
      nVtxFrom3ProngFitter = 0;
    }
    if (doTiming) {
      timing[TimeFit3Prong] += elapsedSeconds(tStart3Prong);
    }
    if (nVtxFrom3ProngFitter == 0) {
      return;
    }

    // get secondary vertex
    const auto& secondaryVertex3 = df3.getPCACandidate();
    // get track momenta
    std::array<float, 3> pvec0{};
    std::array<float, 3> pvec1{};
    std::array<float, 3> pvec2{};
    df3.getTrack(0).getPxPyPzGlo(pvec0);
    df3.getTrack(1).getPxPyPzGlo(pvec1);
    df3.getTrack(2).getPxPyPzGlo(pvec2);
    const auto pVecCandProng3 = RecoDecay::pVec(pvec0, pvec1, pvec2);

    // 3-prong selections after secondary vertex
    const std::array pvCoord{collision.posX(), collision.posY(), collision.posZ()};
    applySelection3Prong(pVecCandProng3, secondaryVertex3, pvCoord, isSelected3ProngCand);
    if (isSelected3ProngCand == 0) {
      return;
    }

    // fill table row
    rowTrackIndexProng3(collision.globalIndex(), globalIndex0, globalIndex1, globalIndex2, hfFlagOfSelection(isSelected3ProngCand));

    // fill histograms
    if (config.fillHistograms) {
      registry.fill(HIST("hVtx3ProngX"), secondaryVertex3[0]);
      registry.fill(HIST("hVtx3ProngY"), secondaryVertex3[1]);
      registry.fill(HIST("hVtx3ProngZ"), secondaryVertex3[2]);
      const std::array arr3Mom{pvec0, pvec1, pvec2};
      for (int iChannel = 0; iChannel < NChannels3Prong; iChannel++) {
        if (!TESTBIT(isSelected3ProngCand, iChannel)) {
          continue;
        }
        if (TESTBIT(whichHypo3Prong[iChannel], 0)) {
          const auto mass3Prong = RecoDecay::m(arr3Mom, arrMass3Prong[iChannel][0]);
          if (iChannel == ChannelDplusToPiKPi) {
            registry.fill(HIST("hMassDPlusToPiKPi"), mass3Prong);
          } else {
            registry.fill(HIST("hMassDsToKKPi"), mass3Prong);
          }
        }
        // the two D+ hypotheses are identical, only the Ds has a second one worth filling
        if (iChannel == ChannelDsToKKPi && TESTBIT(whichHypo3Prong[iChannel], 1)) {
          registry.fill(HIST("hMassDsToKKPi"), RecoDecay::m(arr3Mom, arrMass3Prong[iChannel][1]));
        }
      }
    }
  }

  /// Build the per-collision cache of prongs propagated to the collision's primary vertex.
  /// \param collision is the collision under study
  /// \param trackIndices are the track associations of this collision
  /// \param cacheThis is the prong cache to be filled
  template <typename TTrackIndices>
  void fillProngCache(SelectedCollisions::iterator const& collision,
                      TTrackIndices const& trackIndices, HfProngCache& cacheThis)
  {
    const auto tStartCache = tickIf(doTiming);
    const auto thisCollId = collision.globalIndex();
    for (const auto& trackIndex : trackIndices) {
      const auto track = trackIndex.template track_as<aod::TracksWCovDca>();
      HfProngCandidate prong{getTrackParCov(track), track.pVector(), trackIndex.isIdentifiedPid(), track.globalIndex()};
      ++counters[CountTracksCached];
      if (thisCollId != track.collisionId()) { // this is not the "default" collision for this track, we have to re-propagate it
        const auto tStartPropagate = tickIf(doTiming);
        ++counters[CountTracksPropagated];
        std::array dcaInfo{track.dcaXY(), track.dcaZ()};
        o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, prong.parCov, 2.f, noMatCorr, &dcaInfo);
        getPxPyPz(prong.parCov, prong.pVec);
        if (doTiming) {
          timing[TimeCachePropagate] += elapsedSeconds(tStartPropagate);
        }
      }
      const int index = static_cast<int>(cacheThis.prongs.size());
      cacheThis.prongs.push_back(prong);
      if ((prong.mask & MaskSameSignAny) != 0u) {
        cacheThis.sameSign.push_back(index);
      }
      if ((prong.mask & MaskOppSignAny) != 0u) {
        cacheThis.oppSign.push_back(index);
      }
    }
    if (doTiming) {
      timing[TimeCacheFill] += elapsedSeconds(tStartCache);
    }
  }

  /// Run the combinatorial loop over the two same-sign prongs and the opposite-sign prong.
  /// \param collision is the collision under study
  /// \param sameSignCache holds the prongs of the charge the candidate carries; its sameSign
  ///        list supplies the two same-sign slots
  /// \param oppSignCache holds the opposite charge; its oppSign list supplies the kaon slot
  void runCombinatorics(SelectedCollisions::iterator const& collision,
                        HfProngCache const& sameSignCache, HfProngCache const& oppSignCache)
  {
    const std::vector<int>& sameSignList = sameSignCache.sameSign;
    const std::vector<int>& oppSignList = oppSignCache.oppSign;
    constexpr std::size_t NSameSignProngs = NProngs - 1; // two of the three prongs carry the candidate charge
    if (sameSignList.size() < NSameSignProngs || oppSignList.empty()) {
      return;
    }
    const auto tStartLoop = tickIf(doTiming);
    const std::array pvCoord2Prong{collision.posX(), collision.posY(), collision.posZ()};
    const bool useTwoTrackVertex = config.minTwoTrackDecayLengthFor3Prongs > 0.f || config.maxTwoTrackChi2PcaFor3Prongs < 1.e9f; // o2-linter: disable="magic-number" (default maxTwoTrackChi2PcaFor3Prongs is 1.e10)

    for (const int iOpp : oppSignList) { // o2-linter: disable=const-ref-in-for-loop (int elements)
      const auto& prongOpp = oppSignCache.prongs[iOpp];

      for (std::size_t i1 = 0; i1 + 1 < sameSignList.size(); ++i1) {
        const auto& prong1 = sameSignCache.prongs[sameSignList[i1]];
        ++counters[CountPairsSeen];
        if (!isCandidateAllowed(prongOpp.mask, prong1.mask)) {
          continue;
        }
        ++counters[CountPairsCandidateOk];

        // optional 2-track vertex to skip the pair early
        if (useTwoTrackVertex) {
          const auto tStart2Prong = tickIf(doTiming);
          ++counters[CountFit2Prong];
          int nVtxFrom2ProngFitter = 0;
          try {
            nVtxFrom2ProngFitter = df2.process(prong1.parCov, prongOpp.parCov);
          } catch (...) {
          }
          const bool pairOk = (nVtxFrom2ProngFitter != 0) &&
                              isTwoTrackVertexSelectedFor3Prongs(df2.getPCACandidate(), pvCoord2Prong, df2);
          if (doTiming) {
            timing[TimeFit2Prong] += elapsedSeconds(tStart2Prong);
          }
          if (!pairOk) {
            ++counters[CountPairsRejected2P];
            continue;
          }
        }

        for (std::size_t i2 = i1 + 1; i2 < sameSignList.size(); ++i2) {
          const auto& prong2 = sameSignCache.prongs[sameSignList[i2]];
          if (!tripletRoleOk(prongOpp.mask, prong1.mask, prong2.mask)) {
            continue;
          }
          ++counters[CountTripletsEnumerated];

          processTriplet(collision, prong1.parCov, prongOpp.parCov, prong2.parCov,
                         prong1.pVec, prongOpp.pVec, prong2.pVec,
                         prong1.mask, prongOpp.mask, prong2.mask,
                         prong1.globalIndex, prongOpp.globalIndex, prong2.globalIndex);
        }
      }
    }
    if (doTiming) {
      timing[TimeLoopTotal] += elapsedSeconds(tStartLoop);
    }
  }

  void process3Prongs(SelectedCollisions const& collisions,
                      aod::BCsWithTimestamps const&,
                      FilteredTrackAssocSel const&,
                      aod::TracksWCovDca const& /*tracks*/)
  {
    for (const auto& collision : collisions) {
      // set the magnetic field from CCDB
      const auto bc = collision.bc_as<o2::aod::BCsWithTimestamps>();
      initCCDB(bc, runNumber, ccdb, config.ccdbPathGrpMag, lut, false);
      df2.setBz(o2::base::Propagator::Instance()->getNominalBz());
      df3.setBz(o2::base::Propagator::Instance()->getNominalBz());

      // used to calculate number of candidates per event
      auto nCand3 = rowTrackIndexProng3.lastIndex();

      timing.fill(0.);
      counters.fill(0);

      const auto thisCollId = collision.globalIndex();

      const auto allPos = positiveFor3Prongs->sliceByCached(aod::track::collisionId, thisCollId, cache);
      const auto allNeg = negativeFor3Prongs->sliceByCached(aod::track::collisionId, thisCollId, cache);
      const int nTracks = allPos.size() + allNeg.size();

      // We fill the prong cache for each collision, so that the prongs are propagated to the PV only once.
      cachePos.clear();
      cacheNeg.clear();
      fillProngCache(collision, allPos, cachePos);
      fillProngCache(collision, allNeg, cacheNeg);

      // D+ -> pi+ K- pi+ / Ds -> K+ K- pi+ and their charge conjugates: the same-sign pair carries
      // the charge of the candidate, the opposite-sign prong is the kaon in both channels
      runCombinatorics(collision, cachePos, cacheNeg);
      runCombinatorics(collision, cacheNeg, cachePos);

      nCand3 = rowTrackIndexProng3.lastIndex() - nCand3; // number of 3-prong candidates in this collision
      counters[CountTripletsWritten] = nCand3;

      if (config.fillHistograms) {
        // Fill(bin, weight) so the bins accumulate seconds and counts over the whole run
        for (int iStep = 0; iStep < NTimingSteps; iStep++) {
          registry.fill(HIST("hTiming"), iStep, timing[iStep]);
        }
        for (int iStep = 0; iStep < NLoopCounters; iStep++) {
          registry.fill(HIST("hLoopCounters"), iStep, static_cast<double>(counters[iStep]));
        }
        registry.fill(HIST("hNTracks"), nTracks);
        registry.fill(HIST("hNCand3Prong"), nCand3);
        registry.fill(HIST("hNCand3ProngVsNTracks"), nTracks, nCand3);
      }
    }
  }
  PROCESS_SWITCH(HfMlBasedTrackSelector, process3Prongs, "Process 3-prong skim", true);

  void processNo3Prongs(SelectedCollisions const&)
  {
    // dummy
  }
  PROCESS_SWITCH(HfMlBasedTrackSelector, processNo3Prongs, "Do not process 3-prongs", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<HfTrackSelectorTagSelCollisions>(cfgc));
  workflow.push_back(adaptAnalysisTask<HfTrackSelectorTagSelTracks>(cfgc));
  workflow.push_back(adaptAnalysisTask<HfMlBasedTrackSelector>(cfgc));
  return workflow;
}
