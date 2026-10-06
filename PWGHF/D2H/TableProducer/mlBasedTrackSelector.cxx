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

/// Bit positions of aod::HfSelTrack::isSelProng (only Cand2Prong, Cand3Prong and CandDstar are ever set here).
enum CandidateType {
  Cand2Prong = 0,
  Cand3Prong,
  CandV0bachelor,
  CandDstar,
  CandCascadeBachelor,
  NCandidateTypes
};

/// The 2-prong channels handled by this task.
enum Channels2Prong {
  ChannelD0ToPiK = 0,
  NChannels2Prong
};

/// The 3-prong channels handled by this task.
enum Channels3Prong {
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
  RolePionD0 = 0,
  RoleKaonD0,
  RolePionDplus,
  RoleKaonDplus,
  RolePionDs,
  RoleKaonDs,
  RoleSoftPiDstar,
  NTrackMlRoles
};

constexpr std::array<int, NChannels3Prong> RolePionOfChannel{RolePionDplus, RolePionDs};
constexpr std::array<int, NChannels3Prong> RoleKaonOfChannel{RoleKaonDplus, RoleKaonDs};

constexpr int N2Prongs = 2;   // Number of prongs for 2-prong candidates
constexpr int N3Prongs = 3;   // Number of prongs for 3-prong candidates
constexpr int NMassHypos = 2; // mass hypotheses per channel, i.e. the two orderings of the same-sign pair

/// Role of each track in the candidate, per channel and mass hypothesis.
/// The ordering of the channels should strictly follow Channels2Prong
constexpr TrackMlRole trackRoles2Prongs[NChannels2Prong][NMassHypos][N2Prongs] = {
  {
    {RolePionD0, RoleKaonD0},  // D0 → π± K∓
    {RoleKaonD0, RolePionD0}   // D0 → π± K∓
  }
};

/// Role of each track in the candidate, per channel and mass hypothesis.
/// The ordering of the channels should strictly follow Channels3Prong
constexpr TrackMlRole trackRoles3Prongs[NChannels3Prong][NMassHypos][N3Prongs] = {
  {
    {RolePionDplus, RoleKaonDplus, RolePionDplus},  // D± → π± K∓ π±
    {RolePionDplus, RoleKaonDplus, RolePionDplus}   // D± → π± K∓ π± (swapped, no effect for D+)
  },
  {
    {RoleKaonDs, RoleKaonDs, RolePionDs}, // D_s± → K± K∓ π±
    {RolePionDs, RoleKaonDs, RoleKaonDs}  // D_s± → K± K∓ π± (swapped)
  }
};

constexpr double PtMaxModel = 1.e10;

constexpr uint32_t MaskSameSignAny = (1u << RolePionD0) | (1u << RolePionDplus) | (1u << RolePionDs) | (1u << RoleKaonDs) | (1u << RoleSoftPiDstar);
constexpr uint32_t MaskOppSignAny = (1u << RoleKaonD0) | (1u << RoleKaonDplus) | (1u << RoleKaonDs);

/// Where the combinatorics spends its time, accumulated in seconds into hTiming.
enum TimingStep {
  TimeLoopTotal = 0,  ///< everything inside the triple loop
  TimeFit2Prong,      ///< the 2-track vertex used to skip pairs early
  TimeProcessD0,      ///< D0 candidate processing time
  TimeFit3Prong,      ///< the 3-track vertex of the surviving triplets
  TimeProcessDstar,   ///< D* candidate processing time
  TimeCacheFill,      ///< building the per-collision prong cache, propagation included
  TimeCachePropagate, ///< only the re-propagation inside that build
  NTimingSteps
};

/// How often each stage of the triple loop runs, accumulated into hLoopCounters.
enum LoopCounter {
  CountPairsSeen = 0,      ///< (same-sign, opposite-sign) pairs reached
  CountPairsCandidateOk,   ///< ... of those, surviving the ML pair test
  CountFit2Prong,          ///< ... of those, handed to the 2-prong fitter
  CountPairsWritten,       ///< ... of those, accepted as 2-prong candidates
  CountPairsRejected2P,    ///< ... of those, dropped by the 2-track vertex
  CountTripletsEnumerated, ///< triplets built from the surviving pairs
  CountTripletsWritten,    ///< ... of those, written to the skim
  CountDstarsEnumerated,   ///< triplets built from the surviving pairs
  CountDstarsWritten,      ///< ... of those, written to the skim
  CountTracksCached,       ///< associations put into the prong cache
  CountTracksPropagated,   ///< ... of those, needing a re-propagation to this collision
  NLoopCounters
};

/// Where the ML track selection spends its time, accumulated in seconds into its own hTiming.
enum TrackTimingStep {
  TrackTimeTotal = 0,  ///< the whole per-collision loop: slicing, features, models, table filling
  TrackTimeFeatures,   ///< building the input features, re-propagation included
  TrackTimeModelD0,    ///< evaluating the D0 model: input vector, ONNX call, thresholds
  TrackTimeModelDplus, ///< evaluating the D+ model: input vector, ONNX call, thresholds
  TrackTimeModelDs,    ///< evaluating the Ds model
  TrackTimeModelDstar, ///< evaluating the D* model
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
    Configurable<float> etaMaxTrack{"etaMaxTrack", 0.8f, "max. track eta entering the charm combinatorics"};
    Configurable<float> ptMinSoftPi{"ptMinSoftPi", 0.1f, "min. soft pion track pT entering the charm combinatorics"};
    Configurable<bool> enableTiming{"enableTiming", false, "fill hTiming with the CPU of the feature building and of each model evaluation (adds two clock reads per call)"};
    // D0 model
    Configurable<bool> applyMlD0{"applyMlD0", true, "evaluate the D0 track model"};
    Configurable<std::string> onnxFileNameD0{"onnxFileNameD0", "ModelHandler_D0Tracks.onnx", "ONNX file name of the D0 track model"};
    Configurable<std::vector<std::string>> inputFeaturesD0{"inputFeaturesD0", std::vector<std::string>{"pt", "eta", "dcaXY", "dcaZ", "sigmaDcaXY", "sigmaDcaZ", "normDcaXY", "normDcaZ", "signed1Pt", "tgl", "sign", "isPvContributor", "itsNCls", "itsNClsInnerBarrel", "itsChi2NCl", "tpcNClsFound", "tpcCrossedRowsOverFindableCls", "tpcChi2NCl", "tpcFractionSharedCls", "tpcNSigmaPi", "tpcNSigmaKa"}, "input features of the D0 model, in the order it expects"};
    Configurable<double> thresholdScorePionD0{"thresholdScorePionD0", 0., "min. D0 pion-class score"};
    Configurable<double> thresholdScoreKaonD0{"thresholdScoreKaonD0", 0., "min. D0 kaon-class score"};
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
    // D* model
    Configurable<bool> applyMlDstar{"applyMlDstar", true, "evaluate the D* track model"};
    Configurable<std::string> onnxFileNameDstar{"onnxFileNameDstar", "ModelHandler_DstarTracks.onnx", "ONNX file name of the D* track model"};
    Configurable<std::vector<std::string>> inputFeaturesDstar{"inputFeaturesDstar", std::vector<std::string>{"pt", "eta", "dcaXY", "dcaZ", "sigmaDcaXY", "sigmaDcaZ", "normDcaXY", "normDcaZ", "signed1Pt", "tgl", "sign", "isPvContributor", "itsNCls", "itsNClsInnerBarrel", "itsChi2NCl", "tpcNClsFound", "tpcCrossedRowsOverFindableCls", "tpcChi2NCl", "tpcFractionSharedCls", "tpcNSigmaPi", "tpcNSigmaKa"}, "input features of the D* model, in the order it expects"};
    Configurable<double> thresholdScorePionDstar{"thresholdScorePionDstar", 0., "min. D* pion-class score"};
    // ONNX runtime
    Configurable<bool> loadModelsFromCcdb{"loadModelsFromCcdb", false, "load the ONNX models from CCDB instead of a local path"};
    Configurable<std::string> mlModelPathCcdbD0{"mlModelPathCcdbD0", "path/to/ml/models/D0", "CCDB path of the D0 track model"};
    Configurable<std::string> mlModelPathCcdbDplus{"mlModelPathCcdbDplus", "path/to/ml/models/Dplus", "CCDB path of the D+ track model"};
    Configurable<std::string> mlModelPathCcdbDs{"mlModelPathCcdbDs", "path/to/ml/models/Ds", "CCDB path of the Ds track model"};
    Configurable<std::string> mlModelPathCcdbDstar{"mlModelPathCcdbDstar", "path/to/ml/models/Dstar", "CCDB path of the D* track model"};
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
  std::array<HfTrackModel, NChannels2Prong> models2Prong;
  std::array<HfTrackModel, NChannels3Prong> models3Prong;
  HfTrackModel modelSoftPiDstar;

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
  /// \param ccdbPath is the CCDB path of the model, used if the models are loaded from CCDB
  /// \param features are the input feature names, in the order the model expects
  /// \param thrPion, thrKaon are the score thresholds
  /// \param name is used in the log message
  void configureModel(HfTrackModel& model,
                      std::string const& onnxFile,
                      std::string const& ccdbPath,
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
      model.response.setModelPathsCCDB(onnxFiles, ccdbApi, std::vector<std::string>{ccdbPath}, config.timestampCcdbForMlModels);
    } else {
      model.response.setModelPathsLocal(onnxFiles);
    }
    model.response.cacheInputFeaturesIndices(features);
    model.response.init(config.enableOnnxOptimizations, config.onnxThreads);
    model.enabled = true;
    LOGP(info, "{}: configured from {} with {} input features and {} output classes (pion score > {}, kaon score > {})",
         name, config.loadModelsFromCcdb ? "CCDB " + ccdbPath : onnxFile, features.size(), nOut, thrPion, thrKaon);
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

    if (config.applyMlD0) {
      configureModel(models2Prong[ChannelD0ToPiK], config.onnxFileNameD0, config.mlModelPathCcdbD0, config.inputFeaturesD0,
                     config.thresholdScorePionD0, config.thresholdScoreKaonD0, "D0 track model");
    }
    if (config.applyMlDplus) {
      configureModel(models3Prong[ChannelDplusToPiKPi], config.onnxFileNameDplus, config.mlModelPathCcdbDplus, config.inputFeaturesDplus,
                     config.thresholdScorePionDplus, config.thresholdScoreKaonDplus, "D+ track model");
    }
    if (config.applyMlDs) {
      configureModel(models3Prong[ChannelDsToKKPi], config.onnxFileNameDs, config.mlModelPathCcdbDs, config.inputFeaturesDs,
                     config.thresholdScorePionDs, config.thresholdScoreKaonDs, "Ds track model");
    }
    if (config.applyMlDstar) {
      double thresholdScoreKaonDstar = 0; // D* model has no kaon score, but the configureModel() function expects a value
      configureModel(modelSoftPiDstar, config.onnxFileNameDstar, config.mlModelPathCcdbDstar, config.inputFeaturesDstar,
                     config.thresholdScorePionDstar, thresholdScoreKaonDstar, "D* track model");
    }
    if (!models2Prong[ChannelD0ToPiK].enabled &&
        !models3Prong[ChannelDplusToPiKPi].enabled &&
        !models3Prong[ChannelDsToKKPi].enabled &&
        !modelSoftPiDstar.enabled) {
      LOGP(fatal, "At least one of the D0, D+ and Ds track models must be enabled!");
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
      registry.add("hPtSoftPiQuality", "track associations passing the pT floor and the quality cut, scored;#it{p}_{T}^{track} (GeV/#it{c});entries", {HistType::kTH1D, {axisPtProng}});
      registry.add("hPtQualityRejColl", "track associations passing the pT floor and the quality cut, not scored: collision fails the HF event selection;#it{p}_{T}^{track} (GeV/#it{c});entries", {HistType::kTH1D, {axisPtProng}});
      registry.add("hPtQualitySoftPiRejColl", "track associations passing the pT floor and the quality cut, not scored: collision fails the HF event selection;#it{p}_{T}^{track} (GeV/#it{c});entries", {HistType::kTH1D, {axisPtProng}});
      registry.add("hPtSelected2Prong", "track associations selected by at least one model;#it{p}_{T}^{track} (GeV/#it{c});entries", {HistType::kTH1D, {axisPtProng}});
      registry.add("hEtaSelected2Prong", "track associations selected by at least one model;#it{#eta};entries", {HistType::kTH1D, {axisEta}});
      registry.add("hPtSelected3Prong", "track associations selected by at least one model;#it{p}_{T}^{track} (GeV/#it{c});entries", {HistType::kTH1D, {axisPtProng}});
      registry.add("hEtaSelected3Prong", "track associations selected by at least one model;#it{#eta};entries", {HistType::kTH1D, {axisEta}});
      registry.add("hPtSelectedSoftPi", "track associations selected as D* soft pions;#it{p}_{T}^{track} (GeV/#it{c});entries", {HistType::kTH1D, {axisPtProng}});
      registry.add("hEtaSelectedSoftPi", "track associations selected as D* soft pions;#it{#eta};entries", {HistType::kTH1D, {axisEta}});
      registry.add("hScorePionD0", "D^{0} pion-class score;#it{p}_{T}^{track} (GeV/#it{c});score;entries", {HistType::kTH2D, {axisPtProng, axisScore}});
      registry.add("hScoreKaonD0", "D^{0} kaon-class score;#it{p}_{T}^{track} (GeV/#it{c});score;entries", {HistType::kTH2D, {axisPtProng, axisScore}});
      registry.add("hScorePionDplus", "D^{#plus} pion-class score;#it{p}_{T}^{track} (GeV/#it{c});score;entries", {HistType::kTH2D, {axisPtProng, axisScore}});
      registry.add("hScoreKaonDplus", "D^{#plus} kaon-class score;#it{p}_{T}^{track} (GeV/#it{c});score;entries", {HistType::kTH2D, {axisPtProng, axisScore}});
      registry.add("hScorePionDs", "D_{s}^{#plus} pion-class score;#it{p}_{T}^{track} (GeV/#it{c});score;entries", {HistType::kTH2D, {axisPtProng, axisScore}});
      registry.add("hScoreKaonDs", "D_{s}^{#plus} kaon-class score;#it{p}_{T}^{track} (GeV/#it{c});score;entries", {HistType::kTH2D, {axisPtProng, axisScore}});
      registry.add("hScoreSoftPiDstar", "D^{*} soft pion-class score;#it{p}_{T}^{track} (GeV/#it{c});score;entries", {HistType::kTH2D, {axisPtProng, axisScore}});
      // seconds, accumulated over the run via Fill(bin, weight); per-call cost = bin / hPtQuality entries
      auto hTiming = registry.add<TH1>("hTiming", "CPU of the ML track selection;;seconds", {HistType::kTH1D, {{NTrackTimingSteps, -0.5f, static_cast<float>(NTrackTimingSteps) - 0.5f}}});
      hTiming->GetXaxis()->SetBinLabel(TrackTimeTotal + 1, "track loop total");
      hTiming->GetXaxis()->SetBinLabel(TrackTimeFeatures + 1, "features");
      hTiming->GetXaxis()->SetBinLabel(TrackTimeModelD0 + 1, "D0 model");
      hTiming->GetXaxis()->SetBinLabel(TrackTimeModelDplus + 1, "D+ model");
      hTiming->GetXaxis()->SetBinLabel(TrackTimeModelDs + 1, "Ds model");
      hTiming->GetXaxis()->SetBinLabel(TrackTimeModelDstar + 1, "D* model");
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
    features.tpcNClsCrossedRows = static_cast<float>(track.tpcNClsCrossedRows());
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
      LOGP(fatal, "the model emitted {} scores, expected {}: only 3-class (bkg/pi/K) models are supported for D0, D+, Ds, only 2-class models are supported for D* candidates", output.size(), scores.size());
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
      const bool isPositive = track.sign() > 0;
      if (config.fillHistograms) {
        registry.fill(HIST("hPtNoCuts"), track.pt());
      }

      // One row per track association, always: aod::HfSelTrack is joined with aod::TrackAssoc
      // downstream, so the two tables have to stay index-aligned. Rejected associations are
      // written with an empty selection mask, not skipped.
      uint32_t statusProng{0};
      uint32_t isIdentifiedPid{0};

      const bool passesQuality = track.pt() >= config.ptMinTrack && std::abs(track.eta()) <= config.etaMaxTrack && track.isGlobalTrackWoDCA();
      const bool passesQualitySoftPi = track.pt() >= config.ptMinSoftPi && std::abs(track.eta()) <= config.etaMaxTrack && track.isGlobalTrackWoDCA();
      if (!scoreCollision) {
        if (config.fillHistograms && passesQuality) {
          registry.fill(HIST("hPtQualityRejColl"), track.pt());
        }
        if (config.fillHistograms && passesQualitySoftPi) {
          registry.fill(HIST("hPtQualitySoftPiRejColl"), track.pt());
        }
        rowSelectedTrack(statusProng, isIdentifiedPid, isPositive);
        continue;
      }

      HfTrackMlFeatures features{};
      auto tStart = tickIf(doTiming);
      fillFeatures(collision, track, centrality, features);
      if (doTiming) {
        timing[TrackTimeFeatures] += elapsedSeconds(tStart);
      }

      std::array<float, 3> scores{};
      bool isSelPion{false}, isSelKaon{false};
      bool from2Prong{false}, from3Prong{false}, fromDstar{false};

      if (passesQualitySoftPi && modelSoftPiDstar.enabled) {
        if (config.fillHistograms) {
          registry.fill(HIST("hPtSoftPiQuality"), track.pt());
        }

        tStart = tickIf(doTiming);
        const bool evaluatedDstar = evaluateModel(modelSoftPiDstar, features, scores, isSelPion, isSelKaon);
        if (doTiming) {
          timing[TrackTimeModelDstar] += elapsedSeconds(tStart);
        }
        if (evaluatedDstar) {
          if (isSelPion) {
            SETBIT(isIdentifiedPid, RoleSoftPiDstar);
            fromDstar = true;
          }
          if (config.fillHistograms) {
            registry.fill(HIST("hScoreSoftPiDstar"), features.pt, scores[1]);
          }
        }
      }

      if (passesQuality) {
        if (config.fillHistograms) {
          registry.fill(HIST("hPtQuality"), track.pt());
        }

        tStart = tickIf(doTiming);
        const bool evaluatedD0 = evaluateModel(models2Prong[ChannelD0ToPiK], features, scores, isSelPion, isSelKaon);
        if (doTiming) {
          timing[TrackTimeModelD0] += elapsedSeconds(tStart);
        }
        if (evaluatedD0) {
          if (isSelPion) {
            SETBIT(isIdentifiedPid, RolePionD0);
            from2Prong = true;
          }
          if (isSelKaon) {
            SETBIT(isIdentifiedPid, RoleKaonD0);
            from2Prong = true;
          }
          if (config.fillHistograms) {
            registry.fill(HIST("hScorePionD0"), features.pt, scores[1]);
            registry.fill(HIST("hScoreKaonD0"), features.pt, scores[2]);
          }
        }

        tStart = tickIf(doTiming);
        const bool evaluatedDplus = evaluateModel(models3Prong[ChannelDplusToPiKPi], features, scores, isSelPion, isSelKaon);
        if (doTiming) {
          timing[TrackTimeModelDplus] += elapsedSeconds(tStart);
        }
        if (evaluatedDplus) {
          if (isSelPion) {
            SETBIT(isIdentifiedPid, RolePionDplus);
            from3Prong = true;
          }
          if (isSelKaon) {
            SETBIT(isIdentifiedPid, RoleKaonDplus);
            from3Prong = true;
          }
          if (config.fillHistograms) {
            registry.fill(HIST("hScorePionDplus"), features.pt, scores[1]);
            registry.fill(HIST("hScoreKaonDplus"), features.pt, scores[2]);
          }
        }

        tStart = tickIf(doTiming);
        const bool evaluatedDs = evaluateModel(models3Prong[ChannelDsToKKPi], features, scores, isSelPion, isSelKaon);
        if (doTiming) {
          timing[TrackTimeModelDs] += elapsedSeconds(tStart);
        }
        if (evaluatedDs) {
          if (isSelPion) {
            SETBIT(isIdentifiedPid, RolePionDs);
            from3Prong = true;
          }
          if (isSelKaon) {
            SETBIT(isIdentifiedPid, RoleKaonDs);
            from3Prong = true;
          }
          if (config.fillHistograms) {
            registry.fill(HIST("hScorePionDs"), features.pt, scores[1]);
            registry.fill(HIST("hScoreKaonDs"), features.pt, scores[2]);
          }
        }
      }

      // the track enters the combinatorics if it can play any role in any of the channels
      if (from2Prong) {
        SETBIT(statusProng, CandidateType::Cand2Prong);
        if (config.fillHistograms) {
          registry.fill(HIST("hPtSelected2Prong"), features.pt);
          registry.fill(HIST("hEtaSelected2Prong"), features.eta);
        }
      }
      if (from3Prong) {
        SETBIT(statusProng, CandidateType::Cand3Prong);
        if (config.fillHistograms) {
          registry.fill(HIST("hPtSelected3Prong"), features.pt);
          registry.fill(HIST("hEtaSelected3Prong"), features.eta);
        }
      }
      if (fromDstar) {
        SETBIT(statusProng, CandidateType::CandDstar);
        if (config.fillHistograms) {
          registry.fill(HIST("hPtSelectedSoftPi"), features.pt);
          registry.fill(HIST("hEtaSelectedSoftPi"), features.eta);
        }
      }
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
  Produces<aod::Hf2Prongs> rowTrackIndexProng2;
  Produces<aod::Hf3Prongs> rowTrackIndexProng3;
  Produces<aod::HfDstars> rowTrackIndexDstar;

  struct : ConfigurableGroup {
    Configurable<bool> fillHistograms{"fillHistograms", true, "fill histograms"};
    Configurable<bool> do2Prongs{"do2Prongs", true, "store 2-prong candidates"};
    Configurable<bool> doDstar{"doDstar", true, "store D* candidates"};
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

    // D0 cuts
    Configurable<std::vector<double>> binsPtD0ToPiK{"binsPtD0ToPiK", std::vector<double>{hf_cuts_presel_2prong::vecBinsPt}, "pT bin limits for D0->piK pT-dependent cuts"};
    Configurable<LabeledArray<double>> cutsD0ToPiK{"cutsD0ToPiK", {hf_cuts_presel_2prong::Cuts[0], hf_cuts_presel_2prong::NBinsPt, hf_cuts_presel_2prong::NCutVars, hf_cuts_presel_2prong::labelsPt, hf_cuts_presel_2prong::labelsCutVar}, "D0->piK selections per pT bin"};
    // D+ cuts
    Configurable<std::vector<double>> binsPtDplusToPiKPi{"binsPtDplusToPiKPi", std::vector<double>{hf_cuts_presel_3prong::vecBinsPt}, "pT bin limits for D+->piKpi pT-dependent cuts"};
    Configurable<LabeledArray<double>> cutsDplusToPiKPi{"cutsDplusToPiKPi", {hf_cuts_presel_3prong::Cuts[0], hf_cuts_presel_3prong::NBinsPt, hf_cuts_presel_3prong::NCutVars, hf_cuts_presel_3prong::labelsPt, hf_cuts_presel_3prong::labelsCutVar}, "D+->piKpi selections per pT bin"};
    // Ds+ cuts
    Configurable<std::vector<double>> binsPtDsToKKPi{"binsPtDsToKKPi", std::vector<double>{hf_cuts_presel_ds::vecBinsPt}, "pT bin limits for Ds+->KKPi pT-dependent cuts"};
    Configurable<LabeledArray<double>> cutsDsToKKPi{"cutsDsToKKPi", {hf_cuts_presel_ds::Cuts[0], hf_cuts_presel_ds::NBinsPt, hf_cuts_presel_ds::NCutVars, hf_cuts_presel_ds::labelsPt, hf_cuts_presel_ds::labelsCutVar}, "Ds+->KKPi selections per pT bin"};
    // D*+ cuts
    Configurable<std::vector<double>> binsPtDstarToD0Pi{"binsPtDstarToD0Pi", std::vector<double>{hf_cuts_presel_dstar::vecBinsPt}, "pT bin limits for D*+->D0pi pT-dependent cuts"};
    Configurable<LabeledArray<double>> cutsDstarToD0Pi{"cutsDstarToD0Pi", {hf_cuts_presel_dstar::Cuts[0], hf_cuts_presel_dstar::NBinsPt, hf_cuts_presel_dstar::NCutVars, hf_cuts_presel_dstar::labelsPt, hf_cuts_presel_dstar::labelsCutVar}, "D*+->D0pi selections per pT bin"};
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
  std::array<std::array<std::array<double, N3Prongs>, NMassHypos>, NChannels3Prong> arrMass3Prong{};
  std::array<std::array<std::array<double, N2Prongs>, NMassHypos>, NChannels2Prong> arrMass2Prong{};
  // cuts and pT binning, one entry per channel
  std::array<LabeledArray<double>, NChannels2Prong> cut2Prong{};
  std::array<std::vector<double>, NChannels2Prong> binsPt2Prong{};
  LabeledArray<double> cutDstar{};
  std::vector<double> binsPtDstar{};
  std::array<LabeledArray<double>, NChannels3Prong> cut3Prong{};
  std::array<std::vector<double>, NChannels3Prong> binsPt3Prong{};

  /// One track of the collision under study, already propagated to that collision's PV.
  struct HfProngCandidate {
    o2::track::TrackParCov parCov;
    std::array<float, 3> pVec{};
    std::array<float, 2> dcaInfo{};   ///< DCA (xy, z) to the primary vertex
    uint32_t mask{};       ///< aod::HfSelTrack::isIdentifiedPid
    int64_t globalIndex{}; ///< index in the track table, written to the skim
  };

  /// Cache of propagated tracks for one collision, split by charge.
  struct HfProngCache {
    std::vector<HfProngCandidate> prongs;
    std::vector<int> sameSign; ///< can sit in a same-sign slot: pi(D0), pi(D+), pi(Ds) or K(Ds), softpi(D*)
    std::vector<int> oppSign;  ///< can sit in the opposite-sign slot: K(D0) K(D+) or K(Ds)

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
  Filter filterSelectTrackIds = ( (aod::hf_sel_track::isSelProng & static_cast<uint32_t>(BIT(CandidateType::Cand2Prong))) != 0u ||
                                  (aod::hf_sel_track::isSelProng & static_cast<uint32_t>(BIT(CandidateType::Cand3Prong))) != 0u ||
                                  (aod::hf_sel_track::isSelProng & static_cast<uint32_t>(BIT(CandidateType::CandDstar))) != 0u );

  Preslice<FilteredTrackAssocSel> trackIndicesPerCollision = aod::track_association::collisionId;

  Partition<FilteredTrackAssocSel> positiveHfTracks = (aod::hf_sel_track::isPositive == true);
  Partition<FilteredTrackAssocSel> negativeHfTracks = (aod::hf_sel_track::isPositive == false);

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    if (!doprocessCandidates) {
      return;
    }
    doTiming = config.enableTiming && config.fillHistograms;

    arrMass2Prong[ChannelD0ToPiK] = std::array{std::array{MassPiPlus, MassKPlus},
                                               std::array{MassKPlus, MassPiPlus}};

    arrMass3Prong[ChannelDplusToPiKPi] = std::array{std::array{MassPiPlus, MassKPlus, MassPiPlus},
                                                    std::array{MassPiPlus, MassKPlus, MassPiPlus}};

    arrMass3Prong[ChannelDsToKKPi] = std::array{std::array{MassKPlus, MassKPlus, MassPiPlus},
                                                std::array{MassPiPlus, MassKPlus, MassKPlus}};

    // cuts retrieved by json, in the order of Channel3Prong
    cut2Prong = {config.cutsD0ToPiK};
    binsPt2Prong = {config.binsPtD0ToPiK};
    cutDstar = {config.cutsDstarToD0Pi};
    binsPtDstar = {config.binsPtDstarToD0Pi};

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
      hTiming->GetXaxis()->SetBinLabel(TimeProcessD0 + 1, "D0 processing");
      hTiming->GetXaxis()->SetBinLabel(TimeFit3Prong + 1, "3-prong vertex fit");
      hTiming->GetXaxis()->SetBinLabel(TimeProcessDstar + 1, "D* processing");
      hTiming->GetXaxis()->SetBinLabel(TimeCacheFill + 1, "prong cache fill");
      hTiming->GetXaxis()->SetBinLabel(TimeCachePropagate + 1, "cache re-propagation");
      auto hLoopCounters = registry.add<TH1>("hLoopCounters", "Triple loop stages;;entries", {HistType::kTH1D, {{NLoopCounters, -0.5f, static_cast<float>(NLoopCounters) - 0.5f}}});
      hLoopCounters->GetXaxis()->SetBinLabel(CountPairsSeen + 1, "pairs seen");
      hLoopCounters->GetXaxis()->SetBinLabel(CountPairsCandidateOk + 1, "pairs passing ML roles");
      hLoopCounters->GetXaxis()->SetBinLabel(CountFit2Prong + 1, "2-prong fits");
      hLoopCounters->GetXaxis()->SetBinLabel(CountPairsWritten + 1, "pairs written");
      hLoopCounters->GetXaxis()->SetBinLabel(CountPairsRejected2P + 1, "pairs rejected by 2-prong");
      hLoopCounters->GetXaxis()->SetBinLabel(CountTripletsEnumerated + 1, "triplets enumerated");
      hLoopCounters->GetXaxis()->SetBinLabel(CountTripletsWritten + 1, "triplets written");
      hLoopCounters->GetXaxis()->SetBinLabel(CountDstarsEnumerated + 1, "D* enumerated");
      hLoopCounters->GetXaxis()->SetBinLabel(CountDstarsWritten + 1, "D* written");
      hLoopCounters->GetXaxis()->SetBinLabel(CountTracksCached + 1, "tracks cached");
      hLoopCounters->GetXaxis()->SetBinLabel(CountTracksPropagated + 1, "tracks re-propagated");
      registry.add("hVtx2ProngX", "2-prong candidates;#it{x}_{sec. vtx.} (cm);entries", {HistType::kTH1D, {{1000, -2., 2.}}});
      registry.add("hVtx2ProngY", "2-prong candidates;#it{y}_{sec. vtx.} (cm);entries", {HistType::kTH1D, {{1000, -2., 2.}}});
      registry.add("hVtx2ProngZ", "2-prong candidates;#it{z}_{sec. vtx.} (cm);entries", {HistType::kTH1D, {{1000, -20., 20.}}});
      registry.add("hNCand2Prong", "2-prong candidates preselected;# of candidates;entries", {HistType::kTH1D, {axisNumCands}});
      registry.add("hNCand2ProngVsNTracks", "2-prong candidates preselected;# of selected tracks;# of candidates;entries", {HistType::kTH2D, {axisNumTracks, axisNumCands}});
      registry.add("hMassD0ToPiK", "D^{0} candidates;inv. mass (#pi K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{500, 0., 5.}}});
      registry.add("hVtx3ProngX", "3-prong candidates;#it{x}_{sec. vtx.} (cm);entries", {HistType::kTH1D, {{1000, -2., 2.}}});
      registry.add("hVtx3ProngY", "3-prong candidates;#it{y}_{sec. vtx.} (cm);entries", {HistType::kTH1D, {{1000, -2., 2.}}});
      registry.add("hVtx3ProngZ", "3-prong candidates;#it{z}_{sec. vtx.} (cm);entries", {HistType::kTH1D, {{1000, -20., 20.}}});
      registry.add("hNCand3Prong", "3-prong candidates preselected;# of candidates;entries", {HistType::kTH1D, {axisNumCands}});
      registry.add("hNCand3ProngVsNTracks", "3-prong candidates preselected;# of selected tracks;# of candidates;entries", {HistType::kTH2D, {axisNumTracks, axisNumCands}});
      registry.add("hNCandDstar", "D* candidates preselected;# of candidates;entries", {HistType::kTH1D, {axisNumCands}});
      registry.add("hMassDPlusToPiKPi", "D^{#plus} candidates;inv. mass (#pi K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{500, 0., 5.}}});
      registry.add("hMassDsToKKPi", "D_{s}^{#plus} candidates;inv. mass (K K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{500, 0., 5.}}});
      registry.add("hMassDstarToD0Pi", "D*^{#plus} candidates;inv. mass (D0 #pi) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {{130, 0.135, 0.2}}});
    }
  }

  /// Check whether the two prongs can satisfy any mass hypothesis of either channel, based on their ML role masks.
  /// \param maskOpp is the isIdentifiedPid mask of the opposite-sign prong
  /// \param maskSame is the isIdentifiedPid mask of the first same-sign prong
  static void checkCandidateAllowed(const uint32_t maskOpp, const uint32_t maskSame, bool& canBe2Prong, bool& canBe3Prong)
  {
    canBe2Prong = TESTBIT(maskOpp, RoleKaonD0) && TESTBIT(maskSame, RolePionD0);
    // For D* we check that we can have a D0 and later apply the selection for the soft pion

    const bool okDplus = TESTBIT(maskOpp, RoleKaonDplus) && TESTBIT(maskSame, RolePionDplus);
    const bool okDs = TESTBIT(maskOpp, RoleKaonDs) && (TESTBIT(maskSame, RoleKaonDs) || TESTBIT(maskSame, RolePionDs));
    canBe3Prong = okDplus || okDs;
  }

  /// Check whether the three ML masks can satisfy any mass hypothesis of either channel
  /// \param maskOpp is the mask of the opposite-sign prong (prong 1)
  /// \param mask1, mask2 are the masks of the same-sign prongs (prongs 0 and 2)
  static bool pairRoleOk(const uint32_t maskOpp, const uint32_t maskSame)
  {
    const std::array<uint32_t, N2Prongs> masks{maskSame, maskOpp};
    for (int iChannel = 0; iChannel < NChannels2Prong; iChannel++) {
      for (int iHypo = 0; iHypo < NMassHypos; iHypo++) {
        bool isHypoAlive = true;
        for (int iProng = 0; iProng < N2Prongs; iProng++) {
          const int role = trackRoles2Prongs[iChannel][iHypo][iProng];
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

  /// Check whether the three ML masks can satisfy any mass hypothesis of either channel
  /// \param maskOpp is the mask of the opposite-sign prong (prong 1)
  /// \param mask1, mask2 are the masks of the same-sign prongs (prongs 0 and 2)
  static bool tripletRoleOk(const uint32_t maskOpp, const uint32_t mask1, const uint32_t mask2)
  {
    const std::array<uint32_t, N3Prongs> masks{mask1, maskOpp, mask2};
    for (int iChannel = 0; iChannel < NChannels3Prong; iChannel++) {
      for (int iHypo = 0; iHypo < NMassHypos; iHypo++) {
        bool isHypoAlive = true;
        for (int iProng = 0; iProng < N3Prongs; iProng++) {
          const int role = trackRoles3Prongs[iChannel][iHypo][iProng];
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
  /// \param isIdentifiedPid is the aod::HfSelTrack::isIdentifiedPid mask of the prong, in
  ///        the (same-sign, opposite-sign, same-sign) ordering
  /// \param whichHypo information of the mass hypotheses that are still alive
  /// \param isSelected is a bitmap with selection outcome
  template <std::size_t NProngs, std::size_t NCandChannels>
  void applyMlRoleSelection(const std::array<uint32_t, NProngs>& isIdentifiedPid,
                            std::array<int, NCandChannels>& whichHypo,
                            auto& isSelected)
  {
    for (size_t iChannel = 0; iChannel < NCandChannels; iChannel++) {
      for (size_t iHypo = 0; iHypo < NMassHypos; iHypo++) {
        for (size_t iProng = 0; iProng < NProngs; iProng++) {
          int role{-1};
          if constexpr (NProngs == 2) {
            role = trackRoles2Prongs[iChannel][iHypo][iProng];
          } else if constexpr (NProngs == 3) {
            role = trackRoles3Prongs[iChannel][iHypo][iProng];
          }
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

    // Precompute, no dependence on channel or mass hypothesis
    const auto pt = RecoDecay::pt(pVecCand);
    const auto cpa = RecoDecay::cpa(primVtx, secVtx, pVecCand);
    const auto decayLength = RecoDecay::distance(primVtx, secVtx);

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
      if (cpa < cut3Prong[iChannel].get(binPt, 2u)) { // 2u == "cosp"
        CLRBIT(isSelected, iChannel);
        continue;
      }

      // decay length
      if (decayLength < cut3Prong[iChannel].get(binPt, 3u)) { // 3u == "decL"
        CLRBIT(isSelected, iChannel);
      }
    }
  }

  /// Translate the compact per-channel selection bitmap into the hfflag written to aod::Hf3Prongs,
  /// which uses the hf_cand_3prong::DecayType bit positions.
  uint8_t hfFlagOfSelection(const uint32_t isSelected, bool is3Prong) const
  {
    uint8_t hfFlag{0};
    size_t nChannels = is3Prong ? static_cast<size_t>(NChannels3Prong) : static_cast<size_t>(NChannels2Prong);
    for (size_t iChannel = 0; iChannel < nChannels; iChannel++) {
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

    std::array<uint32_t, N3Prongs> isIdentifiedPid{isIdentifiedPid0, isIdentifiedPid1, isIdentifiedPid2};
    applyMlRoleSelection(isIdentifiedPid, whichHypo3Prong, isSelected3ProngCand);
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
    rowTrackIndexProng3(collision.globalIndex(), globalIndex0, globalIndex1, globalIndex2, hfFlagOfSelection(isSelected3ProngCand, true));

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

  /// Method to perform selections for 2-prong candidates before vertex reconstruction
  /// \param pVecTrack0 is the momentum array of the first daughter track
  /// \param pVecTrack1 is the momentum array of the second daughter track
  /// \param dcaTrack0 is the dcaXY of the first daughter track
  /// \param dcaTrack1 is the dcaXY of the second daughter track
  /// \param cutStatus is a 2D array with outcome of each selection (filled only in debug mode)
  /// \param whichHypo information of the mass hypoteses that were selected
  /// \param isSelected is a bitmap with selection outcome
  /// \param pt2Prong is the pt of the 2-prong candidate
  // void applySelection2Prong(const T1& pVecCand, const T2& secVtx, const T3& primVtx, T4& cutStatus, auto& isSelected)
  template <typename T1, typename T2, typename T3, typename T4, typename T5>
  void applySelection2Prong(const T1& pVecCand, T2 const& prong0, T2 const& prong1, const T3& secVtx, const T4& primVtx, T5& whichHypo, auto& isSelected)
  {

    const auto pt2Prong = RecoDecay::pt(pVecCand[0], pVecCand[1]);
    const auto impParProduct = prong0.dcaInfo[0] * prong1.dcaInfo[0];
    const auto cpa = RecoDecay::cpa(primVtx, secVtx, pVecCand);

    for (int iDecay2P = 0; iDecay2P < NChannels2Prong; iDecay2P++) {

      // return immediately if it is outside the defined pT bins
      const auto binPt = findBin(&binsPt2Prong[iDecay2P], pt2Prong);
      if (binPt == -1) {
        CLRBIT(isSelected, iDecay2P);
        continue;
      }

      // invariant mass
      double massHypos[2] = {0., 0.};

      if (TESTBIT(isSelected, iDecay2P)) {
        const double minMass = cut2Prong[iDecay2P].get(binPt, 0u);
        const double maxMass = cut2Prong[iDecay2P].get(binPt, 1u);
        if (minMass >= 0. && maxMass > 0.) {
          const std::array<std::array<float, 3>, N2Prongs> arr2Mom{prong0.pVec, prong1.pVec};
          massHypos[0] = RecoDecay::m2(arr2Mom, arrMass2Prong[iDecay2P][0]);
          massHypos[1] = RecoDecay::m2(arr2Mom, arrMass2Prong[iDecay2P][1]);
          const double min2 = minMass * minMass;
          const double max2 = maxMass * maxMass;
          if (massHypos[0] < min2 || massHypos[0] >= max2) {
            CLRBIT(whichHypo[iDecay2P], 0);
          }
          if (massHypos[1] < min2 || massHypos[1] >= max2) {
            CLRBIT(whichHypo[iDecay2P], 1);
          }
          if (whichHypo[iDecay2P] == 0) {
            CLRBIT(isSelected, iDecay2P);
          }
        }
      }

      // imp. par. product cut
      if (TESTBIT(isSelected, iDecay2P)) {
        if (impParProduct > cut2Prong[iDecay2P].get(binPt, 3u)) {
          CLRBIT(isSelected, iDecay2P);
        }
      }

      // cos of pointing angle
      if (TESTBIT(isSelected, iDecay2P)) {
        if (cpa < cut2Prong[iDecay2P].get(binPt, 2u)) { // 2u == "cospIndex[iDecay2P]"
          CLRBIT(isSelected, iDecay2P);
        }
      }
    }
  }

  /// Reconstruct and preselect one 2-prong candidate, and fill its table row if it survives.
  /// \param collision is the collision the candidate belongs to
  /// \param prong0,1 are the track informations
  template <typename TProngInfo>
  void processPair(SelectedCollisions::iterator const& collision,
                   TProngInfo const& prong0, TProngInfo const& prong1)
  {
    uint32_t isSelected2ProngCand = BIT(NChannels2Prong) - 1;
    std::array<int, NChannels2Prong> whichHypo2Prong{};
    whichHypo2Prong.fill(BIT(NMassHypos) - 1); // 2 bits on, all mass hypotheses alive

    std::array<uint32_t, N2Prongs> isIdentifiedPid{prong0.mask, prong1.mask};
    applyMlRoleSelection(isIdentifiedPid, whichHypo2Prong, isSelected2ProngCand);
    if (isSelected2ProngCand == 0) {
      return;
    }

    // get secondary vertex
    const auto& secondaryVertex2 = df2.getPCACandidate();
    // get track momenta
    std::array<float, 3> pvec0{}, pvec1{};
    df2.getTrack(0).getPxPyPzGlo(pvec0);
    df2.getTrack(1).getPxPyPzGlo(pvec1);
    const auto pVecCandProng2 = RecoDecay::pVec(pvec0, pvec1);

    // 2-prong selections
    const std::array pvCoord{collision.posX(), collision.posY(), collision.posZ()};
    applySelection2Prong(pVecCandProng2, prong0, prong1, secondaryVertex2, pvCoord, whichHypo2Prong, isSelected2ProngCand);

    // Only D0 for now
    if (!TESTBIT(isSelected2ProngCand, ChannelD0ToPiK)) {
      return;
    }

    // fill table row
    rowTrackIndexProng2(collision.globalIndex(), prong0.globalIndex, prong1.globalIndex, hfFlagOfSelection(isSelected2ProngCand, false));
    ++counters[CountPairsWritten];

    // fill histograms
    if (config.fillHistograms) {
      registry.fill(HIST("hVtx2ProngX"), secondaryVertex2[0]);
      registry.fill(HIST("hVtx2ProngY"), secondaryVertex2[1]);
      registry.fill(HIST("hVtx2ProngZ"), secondaryVertex2[2]);
      const std::array arr2Mom{pvec0, pvec1};
      registry.fill(HIST("hMassD0ToPiK"), RecoDecay::m(arr2Mom, arrMass2Prong[ChannelD0ToPiK][0]));
      registry.fill(HIST("hMassD0ToPiK"), RecoDecay::m(arr2Mom, arrMass2Prong[ChannelD0ToPiK][1]));
    }
  }

  /// Method to perform selections for D* candidates before vertex reconstruction
  /// \param collision is the collision the candidate belongs to
  /// \param lastFilledD0 is the index of the last filled D0 candidate in the table
  /// \param softPiGlobalIndex is the global index of the soft pion track
  /// \param pVecTrack0 is the momentum array of the first daughter track (same charge)
  /// \param pVecTrack1 is the momentum array of the second daughter track (opposite charge)
  /// \param pVecTrack2 is the momentum array of the third daughter track (same charge)
  void processDstar(SelectedCollisions::iterator const& collision, int lastFilledD0, int64_t softPiGlobalIndex,
                    std::array<float, 3> const& pVecTrack0, std::array<float, 3> const& pVecTrack1, std::array<float, 3> const& pVecTrack2)
  {
    const std::array arrMom{pVecTrack0, pVecTrack1, pVecTrack2};
    const std::array arrMomD0{pVecTrack0, pVecTrack1};
    const auto pt = RecoDecay::pt(pVecTrack0, pVecTrack1, pVecTrack2) + config.ptTolerance; // add tolerance because of no reco decay vertex

    // pT
    const auto binPt = findBin(config.binsPtDstarToD0Pi, pt);
    // return immediately if it is outside the defined pT bins
    if (binPt == -1) {
      return;
    }

    // D0 mass
    const double deltaMassD0 = config.cutsDstarToD0Pi->get(binPt, 1u); // 1u == deltaMassD0Index
    const double invMassD0 = RecoDecay::m(arrMomD0, std::array{MassPiPlus, MassKPlus});
    if (std::abs(invMassD0 - MassD0) > deltaMassD0) {
      return;
    }

    // D*+ mass
    const double maxDeltaMass = config.cutsDstarToD0Pi->get(binPt, 0u); // 0u == deltaMassIndex
    const double invMassDstar = RecoDecay::m(arrMom, std::array{MassPiPlus, MassKPlus, MassPiPlus});
    const double deltaMass = invMassDstar - invMassD0;
    if (deltaMass > maxDeltaMass) {
      return;
    }

    // D* candidate is accepted

    // fill table row and fill histograms
    rowTrackIndexDstar(collision.globalIndex(), softPiGlobalIndex, lastFilledD0);
    if (config.fillHistograms) {
      registry.fill(HIST("hMassDstarToD0Pi"), deltaMass);
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
      HfProngCandidate prong{getTrackParCov(track), track.pVector(), {track.dcaXY(), track.dcaZ()}, trackIndex.isIdentifiedPid(), track.globalIndex()};
      ++counters[CountTracksCached];
      if (thisCollId != track.collisionId()) { // this is not the "default" collision for this track, we have to re-propagate it
        const auto tStartPropagate = tickIf(doTiming);
        ++counters[CountTracksPropagated];
        o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, prong.parCov, 2.f, noMatCorr, &prong.dcaInfo);
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
    constexpr std::size_t NSameSign3Prongs = N3Prongs - 1; // two of the three prongs carry the candidate charge
    if ((!config.do2Prongs && !config.doDstar && sameSignList.size() < NSameSign3Prongs) || oppSignList.empty()) {
      return;
    }
    const auto tStartLoop = tickIf(doTiming);
    const std::array pvCoord2Prong{collision.posX(), collision.posY(), collision.posZ()};
    const bool useTwoTrackVertex = config.minTwoTrackDecayLengthFor3Prongs > 0.f || config.maxTwoTrackChi2PcaFor3Prongs < 1.e9f; // o2-linter: disable="magic-number" (default maxTwoTrackChi2PcaFor3Prongs is 1.e10)

    for (const int iOpp : oppSignList) { // o2-linter: disable=const-ref-in-for-loop (int elements)
      const auto& prongOpp = oppSignCache.prongs[iOpp];

      for (std::size_t i1 = 0; i1 < sameSignList.size(); ++i1) {
        const auto& prong1 = sameSignCache.prongs[sameSignList[i1]];
        ++counters[CountPairsSeen];
        bool canBe2Prong{false}, canBe3Prong{false};
        bool isTwoProngVtxGoodFor3Prongs{true};
        checkCandidateAllowed(prongOpp.mask, prong1.mask, canBe2Prong, canBe3Prong);
        if (!canBe2Prong && !canBe3Prong) {
          continue;
        }
        if (!canBe3Prong && !config.do2Prongs && !config.doDstar) {
          continue;
        }
        ++counters[CountPairsCandidateOk];

        int lastFilledD0 = rowTrackIndexProng2.lastIndex();
        // optional 2-track vertex to skip the pair early
        if (useTwoTrackVertex || config.do2Prongs || config.doDstar) {
          const auto tStart2Prong = tickIf(doTiming);
          ++counters[CountFit2Prong];
          int nVtxFrom2ProngFitter = 0;
          try {
            nVtxFrom2ProngFitter = df2.process(prong1.parCov, prongOpp.parCov);
          } catch (...) {
          }
          if (doTiming) {
            timing[TimeFit2Prong] += elapsedSeconds(tStart2Prong);
          }
          const bool hasVertices = (nVtxFrom2ProngFitter != 0);
          if (!hasVertices) {
            ++counters[CountPairsRejected2P];
            continue;
          }
          const auto tStartD0 = tickIf(doTiming);
          if (config.do2Prongs || config.doDstar) {
            processPair(collision, prong1, prongOpp);
          }
          if (doTiming) {
            timing[TimeProcessD0] += elapsedSeconds(tStartD0);
          }
          isTwoProngVtxGoodFor3Prongs = isTwoTrackVertexSelectedFor3Prongs(df2.getPCACandidate(), pvCoord2Prong, df2);
        }

        // D* -> D0 pi+ candidates
        const auto tStartDstar = tickIf(doTiming);
        if (config.doDstar &&
            lastFilledD0 != rowTrackIndexProng2.lastIndex()) { // we have a new D0 candidate, so we can try to form a D* candidate

          lastFilledD0 = rowTrackIndexProng2.lastIndex();
          // Loop over all tracks to search for soft pions
          for (std::size_t iSoftPi = 0; iSoftPi < sameSignList.size(); ++iSoftPi) {
            const auto& prongSoftPi = sameSignCache.prongs[sameSignList[iSoftPi]];

            if (iSoftPi == i1) {
              continue; // skip the same track
            }
            if (!TESTBIT(prongSoftPi.mask, RoleSoftPiDstar)) {
              continue; // skip if the second same-sign prong is not a soft pion candidate
            }
            ++counters[CountDstarsEnumerated];

            // Update the index of the last filled D0 candidate
            processDstar(collision, lastFilledD0, prongSoftPi.globalIndex,
                         prong1.pVec, prongOpp.pVec, prongSoftPi.pVec);
          }
        }
        if (doTiming) {
          timing[TimeProcessDstar] += elapsedSeconds(tStartDstar);
        }

        if (!canBe3Prong) { continue; }
        if (!isTwoProngVtxGoodFor3Prongs) {
          ++counters[CountPairsRejected2P];
          continue;
        }
        // 3-prong candidates
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

  void processCandidates(SelectedCollisions const& collisions,
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
      auto nCand2 = rowTrackIndexProng2.lastIndex();
      auto nCand3 = rowTrackIndexProng3.lastIndex();
      auto nCandDstars = rowTrackIndexDstar.lastIndex();

      timing.fill(0.);
      counters.fill(0);

      const auto thisCollId = collision.globalIndex();

      const auto allPos = positiveHfTracks->sliceByCached(aod::track::collisionId, thisCollId, cache);
      const auto allNeg = negativeHfTracks->sliceByCached(aod::track::collisionId, thisCollId, cache);
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

      nCand2 = rowTrackIndexProng2.lastIndex() - nCand2; // number of 2-prong candidates in this collision
      counters[CountPairsWritten] = nCand2;
      nCand3 = rowTrackIndexProng3.lastIndex() - nCand3; // number of 3-prong candidates in this collision
      counters[CountTripletsWritten] = nCand3;
      nCandDstars = rowTrackIndexDstar.lastIndex() - nCandDstars; // number of D* candidates in this collision
      counters[CountDstarsWritten] = nCandDstars;

      if (config.fillHistograms) {
        // Fill(bin, weight) so the bins accumulate seconds and counts over the whole run
        for (int iStep = 0; iStep < NTimingSteps; iStep++) {
          registry.fill(HIST("hTiming"), iStep, timing[iStep]);
        }
        for (int iStep = 0; iStep < NLoopCounters; iStep++) {
          registry.fill(HIST("hLoopCounters"), iStep, static_cast<double>(counters[iStep]));
        }
        registry.fill(HIST("hNTracks"), nTracks);
        registry.fill(HIST("hNCand2Prong"), nCand2);
        registry.fill(HIST("hNCand3Prong"), nCand3);
        registry.fill(HIST("hNCandDstar"), nCandDstars);
        registry.fill(HIST("hNCand3ProngVsNTracks"), nTracks, nCand3);
      }
    }
    LOG(info) << "Processed " << collisions.size() << " collisions for 3-prong candidates";
  }
  PROCESS_SWITCH(HfMlBasedTrackSelector, processCandidates, "Process 3-prong skim", true);

  void processDummy(SelectedCollisions const&)
  {
    // dummy
  }
  PROCESS_SWITCH(HfMlBasedTrackSelector, processDummy, "Do not process 3-prongs", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<HfTrackSelectorTagSelCollisions>(cfgc));
  workflow.push_back(adaptAnalysisTask<HfTrackSelectorTagSelTracks>(cfgc));
  workflow.push_back(adaptAnalysisTask<HfMlBasedTrackSelector>(cfgc));
  return workflow;
}
