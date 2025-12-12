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

/// \file tagAndProbeDmesons.cxx
/// \brief Task for tracking efficiency studies with tag-and-probe using 3-prong D-meson decays
///
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

#include "PWGHF/Utils/utilsAnalysis.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Tools/ML/MlResponse.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DCAFitter/DCAFitterN.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include <TPDGCode.h>

#include <algorithm>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace tagandprobe
{
enum TagChannels : uint8_t {
  DplusToKPiPi = 0,
  DsOrDplusToKKPi,
  DstarPlusToDzeroPi,
  DstarMinusToDzeroBarPi,
  DstarToDzeroToKK,
  NTagChannels
};

enum TrackTypes : uint8_t {
  GlobalWoDca = 0,
  GlobalWoDcaWoItsIb,
  GlobalWoDcaWoIts,
  GlobalWoDcaWoTpc,
  NTrackTypes
};

enum SignalFlags : uint8_t {
  Bkg = 0,
  Prompt,
  NonPrompt,
  Resonant,
  BkgFromNoHf
};

static constexpr int nBinsPt = 7;
static constexpr int nCutVars = 9;
constexpr double binsPt[nBinsPt + 1] = {0., 1., 2., 4., 6., 10., 20., 1000.};
auto vecBinsPt = std::vector<double>{binsPt, binsPt + nBinsPt + 1};

// default values for the cuts
constexpr double cuts[nBinsPt][nCutVars] = {{1.815f, 1.915f, 0.01f, 0.01f, 2.f, 2.f, 0.f, 0.90f, 0.90f},
                                            {1.815f, 1.915f, 0.01f, 0.01f, 2.f, 2.f, 0.f, 0.90f, 0.90f},
                                            {1.815f, 1.915f, 0.02f, 0.02f, 2.f, 2.f, 0.f, 0.90f, 0.90f},
                                            {1.815f, 1.915f, 0.02f, 0.02f, 2.f, 2.f, 0.f, 0.90f, 0.90f},
                                            {1.815f, 1.915f, 0.04f, 0.04f, 2.f, 2.f, 0.f, 0.95f, 0.95f},
                                            {1.815f, 1.915f, 0.04f, 0.04f, 2.f, 2.f, 0.f, 0.95f, 0.95f},
                                            {1.815f, 1.915f, 0.06f, 0.06f, 2.f, 2.f, 0.f, 0.95f, 0.95f}};

constexpr double mlCuts[nBinsPt][3] = {{0.f, 0.f, 0.f},
                                       {0.f, 0.f, 0.f},
                                       {0.f, 0.f, 0.f},
                                       {0.f, 0.f, 0.f},
                                       {0.f, 0.f, 0.f},
                                       {0.f, 0.f, 0.f},
                                       {0.f, 0.f, 0.f}};

static const std::vector<std::string> labelsEmpty{};
static const std::vector<std::string> labelsCutVar = {"min mass", "max mass", "min decayLength", "min decayLengthXY", "min normDecayLength", "min normDecayLengthXY", "max impParProd", "min cosPointing", "min cosPointingXY"};
static const std::vector<std::string> labelsMlScores = {"max bkg score", "min prompt score", "min non-prompt score"};

DECLARE_SOA_INDEX_COLUMN(Collision, collision);                   //! Collision index
DECLARE_SOA_INDEX_COLUMN_FULL(Track0, track0, int, Tracks, "_0"); //! Index to first track
DECLARE_SOA_INDEX_COLUMN_FULL(Track1, track1, int, Tracks, "_1"); //! Index to second track
// Topological variables
DECLARE_SOA_COLUMN(TagPt, tagPt, float);                                     //! Tag's pT
DECLARE_SOA_COLUMN(TagInvMass, tagInvMass, float);                           //! Tag's invMass
DECLARE_SOA_COLUMN(DecayLength, decayLength, float);                         //! Decay length of the tag (cm)
DECLARE_SOA_COLUMN(DecayLengthXY, decayLengthXY, float);                     //! Transverse decay length of the tag (cm)
DECLARE_SOA_COLUMN(DecayLengthNormalised, decayLengthNormalised, float);     //! Normalised decay length of the tag
DECLARE_SOA_COLUMN(DecayLengthXYNormalised, decayLengthXYNormalised, float); //! Normalised transverse decay length of the tag
DECLARE_SOA_COLUMN(TrackDcaXY0, trackDcaXY0, float);                         //! DCAxy of the first tag track
DECLARE_SOA_COLUMN(TrackDcaXY1, trackDcaXY1, float);                         //! DCAxy of the second tag track
DECLARE_SOA_COLUMN(ProductTrackDcaXY, productTrackDcaXY, float);             //! Product of DCAxy of the two tag tracks
DECLARE_SOA_COLUMN(Cpa, cpa, float);                                         //! Cosine pointing angle of the tag
DECLARE_SOA_COLUMN(CpaXY, cpaXY, float);                                     //! Cosine of the pointing angle in XY of the tag
DECLARE_SOA_COLUMN(Radius, radius, float);                                   //! Radius of the tag
DECLARE_SOA_COLUMN(DecChannel, decChannel, uint8_t);                         //! Flag the selected decay channel
// MC info
DECLARE_SOA_COLUMN(IsSignal, isSignal, uint8_t);                     //! Flag for a signal
DECLARE_SOA_INDEX_COLUMN_FULL(Mother, mother, int, McParticles, ""); //! Index to MC particle mother of the tag tracks
// ML scores
DECLARE_SOA_COLUMN(MlScores, mlScores, std::vector<float>); //! ML scores (bkg, prompt, non-prompt)
} // namespace tagandprobe

DECLARE_SOA_TABLE(PiPiFromDpTags, "AOD", "PIPIFROMDPTAG", //! Table for same sign 2-pion vertices used as tags
                  soa::Index<>,
                  aod::tagandprobe::CollisionId,
                  aod::tagandprobe::Track0Id,
                  aod::tagandprobe::Track1Id,
                  aod::tagandprobe::MlScores,
                  aod::tagandprobe::Radius);
DECLARE_SOA_TABLE(PiPiFromDpMcTags, "AOD", "PIPIFROMDPMCTAG", //! Table with MC truth for same sign 2-pion vertices used as tags
                  aod::tagandprobe::IsSignal,
                  aod::tagandprobe::MotherId);
DECLARE_SOA_TABLE(KaKaFromDspTags, "AOD", "KAKAFROMDSPTAG", //! Table for opposite sign 2-kaon vertices used as tags
                  soa::Index<>,
                  aod::tagandprobe::CollisionId,
                  aod::tagandprobe::Track0Id,
                  aod::tagandprobe::Track1Id,
                  aod::tagandprobe::MlScores,
                  aod::tagandprobe::Radius,
                  soa::Marker<1>);
DECLARE_SOA_TABLE(KaKaFromDsMcTags, "AOD", "KAKAFROMDSMCTAG", //! Table with MC truth for opposite sign 2-kaon vertices used as tags
                  aod::tagandprobe::IsSignal,
                  aod::tagandprobe::MotherId,
                  soa::Marker<1>);
DECLARE_SOA_TABLE(PiKaFromDzTags, "AOD", "PIKAFROMDZTAG", //! Table for opposite sign pion(+)-kaon(-) vertices used as tags
                  soa::Index<>,
                  aod::tagandprobe::CollisionId,
                  aod::tagandprobe::Track0Id,
                  aod::tagandprobe::Track1Id,
                  aod::tagandprobe::MlScores,
                  aod::tagandprobe::Radius,
                  soa::Marker<2>);
DECLARE_SOA_TABLE(PiKaFromDzMcTags, "AOD", "PIKAFROMDZMCTAG", //! Table with MC truth for opposite sign pion(+)-kaon(-) vertices used as tags
                  aod::tagandprobe::IsSignal,
                  aod::tagandprobe::MotherId,
                  soa::Marker<2>);
DECLARE_SOA_TABLE(KaPiFromDzTags, "AOD", "KAPIFROMDZTAG", //! Table for opposite sign kaon(+)-pion(-) vertices used as tags
                  soa::Index<>,
                  aod::tagandprobe::CollisionId,
                  aod::tagandprobe::Track0Id,
                  aod::tagandprobe::Track1Id,
                  aod::tagandprobe::MlScores,
                  aod::tagandprobe::Radius,
                  soa::Marker<3>);
DECLARE_SOA_TABLE(KaPiFromDzMcTags, "AOD", "KAPIFROMDZMCTAG", //! Table with MC truth for opposite sign kaon(+)-pion(-) vertices used as tags
                  aod::tagandprobe::IsSignal,
                  aod::tagandprobe::MotherId,
                  soa::Marker<3>);
DECLARE_SOA_TABLE(TagTopoVariables, "AOD", "TAGTOPOVARIABLE", //! Table for the Tags' Topological variables
                  aod::tagandprobe::TagPt,
                  aod::tagandprobe::TagInvMass,
                  aod::tagandprobe::DecayLength,
                  aod::tagandprobe::DecayLengthXY,
                  aod::tagandprobe::DecayLengthNormalised,
                  aod::tagandprobe::DecayLengthXYNormalised,
                  aod::tagandprobe::TrackDcaXY0,
                  aod::tagandprobe::TrackDcaXY1,
                  aod::tagandprobe::ProductTrackDcaXY,
                  aod::tagandprobe::Cpa,
                  aod::tagandprobe::CpaXY,
                  aod::tagandprobe::IsSignal,
                  aod::tagandprobe::DecChannel);
} // namespace o2::aod

/// Reconstruction of 2-prong displaced vertices (very good quality and purity)
/// 1) K∓K± for φ from Ds± or D± → φπ± decays
/// 2) π±π± for D± → K∓π±π± decays
/// 3) K∓π± for D0 from D±* → D0π± decays
struct TagTwoProngDisplacedVertices {

  Produces<aod::PiPiFromDpTags> tagPiPiTable;
  Produces<aod::PiPiFromDpMcTags> tagPiPiMcTable;
  Produces<aod::KaKaFromDspTags> tagKaKaTable;
  Produces<aod::KaKaFromDsMcTags> tagKaKaMcTable;
  Produces<aod::KaPiFromDzTags> tagKaPiTable;
  Produces<aod::KaPiFromDzMcTags> tagKaPiMcTable;
  Produces<aod::PiKaFromDzTags> tagPiKaTable;
  Produces<aod::PiKaFromDzMcTags> tagPiKaMcTable;
  Produces<aod::TagTopoVariables> tagVarsTable;

  SliceCache cache;
  Configurable<int> fillTopoVarsTable{"fillTopoVarsTable", 0, "flag to fill tag table with topological variables (0 -> disabled, 1 -> signal only, 2 -> bkg only, 3 -> bkg from no HF only, 4 -> all)"};
  Configurable<float> downsamplingForTopoVarTable{"downsamplingForTopoVarTable", 1.1, "fraction of tag candidates to downscale in filling table with topological variables"};
  Configurable<float> ptTagMaxForDownsampling{"ptTagMaxForDownsampling", 5., "maximum pT for downscaling of tag candidates in filling table with topological variables"};
  Configurable<bool> applyTofPid{"applyTofPid", true, "flag to enable TOF PID selection"};
  Configurable<bool> studyDzeroReflections{"studyDzeroReflections", false, "flag to study Dzero reflections"};
  Configurable<float> trackNumSigmaTof{"trackNumSigmaTof", 3.f, "number of sigma for TOF PID compatibility"};
  Configurable<float> trackNumSigmaTpc{"trackNumSigmaTpc", 3.f, "number of sigma for TOF PID compatibility"};
  Configurable<float> trackDcaXyMin{"trackDcaXyMin", 0.002f, "minimum DCAxy for tracks with pT < 2 GeV/c"};
  Configurable<float> trackPtMin{"trackPtMin", 0.4f, "minimum track pT"};

  Configurable<std::vector<double>> binsPtPiPiFromDplus{"binsPtPiPiFromDplus", std::vector<double>{aod::tagandprobe::vecBinsPt}, "pT bin limits for pipi pairs from D+ decays"};
  Configurable<std::vector<double>> binsPtKaKaFromDsOrDplus{"binsPtKaKaFromDsOrDplus", std::vector<double>{aod::tagandprobe::vecBinsPt}, "pT bin limits for KK pairs from Ds or D+ decays"};
  Configurable<std::vector<double>> binsPtDzeroFromDstar{"binsPtDzeroFromDstar", std::vector<double>{aod::tagandprobe::vecBinsPt}, "pT bin limits for Kpi pairs from D0 <- D*+ decays"};
  Configurable<std::vector<double>> binsPtDzeroKaKaFromDstar{"binsPtDzeroKaKaFromDstar", std::vector<double>{aod::tagandprobe::vecBinsPt}, "pT bin limits for KK pairs from D0 <- D*+ decays"};

  Configurable<LabeledArray<double>> cutsPiPiFromDplus{"cutsPiPiFromDplus", {aod::tagandprobe::cuts[0], aod::tagandprobe::nBinsPt, aod::tagandprobe::nCutVars, aod::tagandprobe::labelsEmpty, aod::tagandprobe::labelsCutVar}, "Selections for pipi pairs from D+ decays"};
  Configurable<LabeledArray<double>> cutsKaKaFromDsOrDplus{"cutsKaKaFromDsOrDplus", {aod::tagandprobe::cuts[0], aod::tagandprobe::nBinsPt, aod::tagandprobe::nCutVars, aod::tagandprobe::labelsEmpty, aod::tagandprobe::labelsCutVar}, "Selections for KK pairs from Ds or D+ decays"};
  Configurable<LabeledArray<double>> cutsDzeroFromDstar{"cutsDzeroFromDstar", {aod::tagandprobe::cuts[0], aod::tagandprobe::nBinsPt, aod::tagandprobe::nCutVars, aod::tagandprobe::labelsEmpty, aod::tagandprobe::labelsCutVar}, "Selections for Kpi pairs from D0 <- D*+ decays"};
  Configurable<LabeledArray<double>> cutsDzeroKaKaFromDstar{"cutsDzeroKaKaFromDstar", {aod::tagandprobe::cuts[0], aod::tagandprobe::nBinsPt, aod::tagandprobe::nCutVars, aod::tagandprobe::labelsEmpty, aod::tagandprobe::labelsCutVar}, "Selections for KK pairs from D0 <- D*+ decays"};

  // ML models for triggers
  struct : ConfigurableGroup {
    std::string prefix = "ML";
    Configurable<bool> applyMlPiPiFromDplus{"applyMlPiPiFromDplus", false, "Flag to enable ML application for pipi pairs from D+ decays"};
    Configurable<bool> applyMlKaKaFromDsOrDplus{"applyMlKaKaFromDsOrDplus", false, "Flag to enable ML application for KK pairs from Ds or D+ decays"};
    Configurable<bool> applyMlDzeroFromDstar{"applyMlDzeroFromDstar", false, "Flag to enable ML application for Kpi pairs from D0 <- D*+ decays"};
    Configurable<bool> applyMlDzeroKaKaFromDstar{"applyMlDzeroKaKaFromDstar", false, "Flag to enable ML application for KK pairs from D0 <- D*+ decays"};
    Configurable<int64_t> timestampCcdbForMlModels{"timestampCcdbForMlModels", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
    Configurable<bool> loadMlModelsFromCcdb{"loadMlModelsFromCcdb", true, "Flag to enable or disable the loading of ML models from CCDB"};
    // ML models (one per pT bin)
    Configurable<std::vector<std::string>> modelPathsCcdbPiPiFromDplus{"modelPathsCcdbPiPiFromDplus", std::vector<std::string>{"/Users/f/fgrosa/TagAndProbe/DplusPt2to3"}, "Paths of models on CCDB for pipi pairs from D+ decays"};
    Configurable<std::vector<std::string>> modelPathsCcdbKaKaFromDsOrDplus{"modelPathsCcdbKaKaFromDsOrDplus", std::vector<std::string>{"/Users/f/fgrosa/TagAndProbe/DsPt2to3"}, "Paths of models on CCDB for KK pairs from Ds or D+ decays"};
    Configurable<std::vector<std::string>> modelPathsCcdbDzeroFromDstar{"modelPathsCcdbDzeroFromDstar", std::vector<std::string>{"/Users/f/fgrosa/TagAndProbe/DzeroPt2to3"}, "Paths of models on CCDB for Kpi pairs from D0 <- D*+ decays"};
    Configurable<std::vector<std::string>> modelPathsCcdbDzeroKaKaFromDstar{"modelPathsCcdbDzeroKaKaFromDstar", std::vector<std::string>{"/Users/f/fgrosa/TagAndProbe/DzeroToKKPt2to3"}, "Paths of models on CCDB for KK pairs from D0 <- D*+ decays"};
    Configurable<std::vector<std::string>> onnxFileNamesPiPiFromDplus{"onnxFileNamesPiPiFromDplus", std::vector<std::string>{"ModelHandler_onnx_DplusToKPiPi.onnx"}, "ONNX file names for pipi pairs from D+ decays"};
    Configurable<std::vector<std::string>> onnxFileNamesKaKaFromDsOrDplus{"onnxFileNamesKaKaFromDsOrDplus", std::vector<std::string>{"ModelHandler_onnx_DsToKKPi.onnx"}, "ONNX file names for KK pairs from Ds or D+ decays"};
    Configurable<std::vector<std::string>> onnxFileNamesDzeroFromDstar{"onnxFileNamesDzeroFromDstar", std::vector<std::string>{"ModelHandler_onnx_DzeroToKPi.onnx"}, "ONNX file names for Kpi pairs from D0 <- D*+ decays"};
    Configurable<std::vector<std::string>> onnxFileNamesDzeroKaKaFromDstar{"onnxFileNamesDzeroKaKaFromDstar", std::vector<std::string>{"ModelHandler_onnx_DzeroToKK.onnx"}, "ONNX file names for KK pairs from D0 <- D*+ decays"};
    // ML cuts
    Configurable<int> numMlClasses{"numMlClasses", 3, "Number of classes for the ML models"};
    Configurable<std::vector<int>> cutDirMl{"cutDirMl", std::vector<int>{o2::cuts_ml::CutDirection::CutGreater, o2::cuts_ml::CutDirection::CutSmaller, o2::cuts_ml::CutDirection::CutSmaller}, "Whether to reject score values greater or smaller than the threshold"};
    Configurable<LabeledArray<double>> mlCutsPiPiFromDplus{"mlCutsPiPiFromDplus", {aod::tagandprobe::mlCuts[0], aod::tagandprobe::nBinsPt, 3, aod::tagandprobe::labelsEmpty, aod::tagandprobe::labelsMlScores}, "ML Selections for pipi pairs from D+ decays"};
    Configurable<LabeledArray<double>> mlCutsKaKaFromDsOrDplus{"mlCutsKaKaFromDsOrDplus", {aod::tagandprobe::mlCuts[0], aod::tagandprobe::nBinsPt, 3, aod::tagandprobe::labelsEmpty, aod::tagandprobe::labelsMlScores}, "ML Selections for KK pairs from Ds or D+ decays"};
    Configurable<LabeledArray<double>> mlCutsDzeroFromDstar{"mlCutsDzeroFromDstar", {aod::tagandprobe::mlCuts[0], aod::tagandprobe::nBinsPt, 3, aod::tagandprobe::labelsEmpty, aod::tagandprobe::labelsMlScores}, "ML Selections for Kpi pairs from D0 <- D*+ decays"};
    Configurable<LabeledArray<double>> mlCutsDzeroKaKaFromDstar{"mlCutsDzeroKaKaFromDstar", {aod::tagandprobe::mlCuts[0], aod::tagandprobe::nBinsPt, 3, aod::tagandprobe::labelsEmpty, aod::tagandprobe::labelsMlScores}, "ML Selections for KK pairs from D0 <- D*+ decays"};
  } mlConfig;

  using TracksWithSelAndDca = soa::Join<aod::Tracks, aod::TracksCov, aod::TracksDCA, aod::TracksExtra, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa>;
  using TracksWithSelAndDcaMc = soa::Join<TracksWithSelAndDca, aod::McTrackLabels>;
  using CollisionsWithEvSel = soa::Join<aod::Collisions, aod::EvSels>;

  Filter evSelFilter = aod::evsel::sel8 == true;                                                                                                              // simple event selection
  Filter collisionFilter = nabs(aod::collision::posZ) < 10.f;                                                                                                 // simple event selection
  Filter trackFilter = requireGlobalTrackWoDCAInFilter() && aod::track::pt > trackPtMin && (nabs(aod::track::dcaXY) > trackDcaXyMin || aod::track::pt > 2.f); // for the tag, we only consider global tracks with large dcaXY (low pT only)
  using TracksWithSelAndDcaFiltered = soa::Filtered<TracksWithSelAndDca>;
  using TracksWithSelAndDcaMcFiltered = soa::Filtered<TracksWithSelAndDcaMc>;
  using CollisionsFiltered = soa::Filtered<CollisionsWithEvSel>;

  // in the partition we only apply TPC PID
  Preslice<TracksWithSelAndDcaFiltered> perCollision = aod::track::collisionId;
  Partition<TracksWithSelAndDcaFiltered> positivePions = aod::track::signed1Pt > 0.f && nabs(aod::pidtpc::tpcNSigmaPi) < trackNumSigmaTpc;
  Partition<TracksWithSelAndDcaFiltered> negativePions = aod::track::signed1Pt < 0.f && nabs(aod::pidtpc::tpcNSigmaPi) < trackNumSigmaTpc;
  Partition<TracksWithSelAndDcaFiltered> positiveKaons = aod::track::signed1Pt > 0.f && nabs(aod::pidtpc::tpcNSigmaKa) < trackNumSigmaTpc;
  Partition<TracksWithSelAndDcaFiltered> negativeKaons = aod::track::signed1Pt < 0.f && nabs(aod::pidtpc::tpcNSigmaKa) < trackNumSigmaTpc;
  Partition<TracksWithSelAndDcaMcFiltered> positivePionsMc = aod::track::signed1Pt > 0.f && nabs(aod::pidtpc::tpcNSigmaPi) < trackNumSigmaTpc;
  Partition<TracksWithSelAndDcaMcFiltered> negativePionsMc = aod::track::signed1Pt < 0.f && nabs(aod::pidtpc::tpcNSigmaPi) < trackNumSigmaTpc;
  Partition<TracksWithSelAndDcaMcFiltered> positiveKaonsMc = aod::track::signed1Pt > 0.f && nabs(aod::pidtpc::tpcNSigmaKa) < trackNumSigmaTpc;
  Partition<TracksWithSelAndDcaMcFiltered> negativeKaonsMc = aod::track::signed1Pt < 0.f && nabs(aod::pidtpc::tpcNSigmaKa) < trackNumSigmaTpc;

  std::array<o2::analysis::MlResponse<float>, aod::tagandprobe::TagChannels::NTagChannels> mlResponse{};
  std::array<bool, aod::tagandprobe::TagChannels::NTagChannels> applyMl{};
  ccdb::CcdbApi ccdbApi;
  Service<ccdb::BasicCCDBManager> ccdb;
  vertexing::DCAFitterN<2> vertexer;
  int runNumber{0};

  std::array<std::array<double, 2>, aod::tagandprobe::TagChannels::NTagChannels> masses = {std::array{constants::physics::MassPionCharged, constants::physics::MassPionCharged},
                                                                                           std::array{constants::physics::MassKaonCharged, constants::physics::MassKaonCharged},
                                                                                           std::array{constants::physics::MassPionCharged, constants::physics::MassKaonCharged},
                                                                                           std::array{constants::physics::MassKaonCharged, constants::physics::MassPionCharged},
                                                                                           std::array{constants::physics::MassKaonCharged, constants::physics::MassKaonCharged}};

  std::array<LabeledArray<double>, aod::tagandprobe::TagChannels::NTagChannels> topologicalCuts{};
  std::array<std::vector<double>, aod::tagandprobe::TagChannels::NTagChannels> ptBinsForTopologicalCuts{};
  std::vector<std::shared_ptr<TH2>> hBkgMlScore{};
  std::vector<std::shared_ptr<TH2>> hPromptMlScore{};
  std::vector<std::shared_ptr<TH2>> hNonPromptMlScore{};
  std::vector<std::shared_ptr<TH2>> hDataMlScore{};

  HistogramRegistry registry{"registry"};

  void init(InitContext&)
  {
    if ((doprocessPiPiFromDplus && doprocessPiPiFromDplusMc) || (doprocessKaKaFromDsOrDplus && doprocessKaKaFromDsOrDplusMc) || (doprocessKaPiFromDstar && doprocessKaPiFromDstarMc)) {
      LOGP(fatal, "The process functions for the same channel with and without MC truth cannot be enabled at the same time! Please check your configuration");
    }

    std::string ccdbUrl = "http://alice-ccdb.cern.ch";
    ccdb->setURL(ccdbUrl.data());
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    ccdbApi.init(ccdbUrl.data());

    vertexer.setPropagateToPCA(true);
    vertexer.setMaxR(200.f);
    vertexer.setMaxDZIni(4.f);
    vertexer.setMinParamChange(1.e-3);
    vertexer.setMinRelChi2Change(0.9f);
    vertexer.setUseAbsDCA(false);

    topologicalCuts = {cutsPiPiFromDplus, cutsKaKaFromDsOrDplus, cutsDzeroFromDstar, cutsDzeroFromDstar, cutsDzeroKaKaFromDstar};
    ptBinsForTopologicalCuts = {binsPtPiPiFromDplus, binsPtKaKaFromDsOrDplus, binsPtDzeroFromDstar, binsPtDzeroFromDstar, binsPtDzeroKaKaFromDstar};

    const AxisSpec axisPt{250, 0.f, 50.f};
    const AxisSpec axisPtDzeroRefl{{0.f, 0.5f, 0.75f, 1.0f, 1.25f, 1.5f, 2.0f, 2.5f, 3.0f, 3.5f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 10.0f, 12.0f, 14.0f, 16.0f, 18.0f, 20.0f, 24.0f, 30.0f, 40.0f}};
    const AxisSpec axisMassPiPi{250, 0.f, 2.5f};
    const AxisSpec axisReflFlag{3, 0.5f, 3.5f};
    const AxisSpec axisMassKaKa{200, constants::physics::MassPhi - 0.05f, constants::physics::MassPhi + 0.05f};
    const AxisSpec axisMassKaPi{400, constants::physics::MassD0 - 0.2f, constants::physics::MassD0 + 0.2f};
    const AxisSpec axisMlScore{1000, 0.f, 1.f};

    if (doprocessPiPiFromDplus || doprocessPiPiFromDplusMc) {
      registry.add<TH2>("hMassPiPiVsPt", ";#it{p}_{T}(#pi#pi) (GeV/#it{c}); #it{M}(#pi#pi) (GeV/#it{c}^{2})", HistType::kTH2D, {axisPt, axisMassPiPi});
    }
    if (doprocessKaKaFromDsOrDplus || doprocessKaKaFromDsOrDplusMc) {
      registry.add<TH2>("hMassKaKaVsPt", ";#it{p}_{T}(KK) (GeV/#it{c}); #it{M}(KK) (GeV/#it{c}^{2})", HistType::kTH2D, {axisPt, axisMassKaKa});
    }
    if (doprocessKaPiFromDstar || doprocessKaPiFromDstarMc) {
      if (!studyDzeroReflections) {
        registry.add<TH2>("hMassKaPiVsPt", ";#it{p}_{T}(K#pi) (GeV/#it{c}); #it{M}(K#pi) (GeV/#it{c}^{2})", HistType::kTH2D, {axisPt, axisMassKaPi});
      } else {
        registry.add<THnSparse>("hMassKaPiVsPt", ";#it{p}_{T}(K#pi) (GeV/#it{c}); #it{M}(K#pi) (GeV/#it{c}^{2}); #it{M}(#piK) (GeV/#it{c}^{2}); ReflFag", HistType::kTHnSparseF, {axisPtDzeroRefl, axisMassKaPi, axisMassKaPi, axisReflFlag});
      }
    }
    if (doprocessKaKaFromDzero) {
      registry.add<TH2>("hMassDzeroKaKaVsPt", ";#it{p}_{T}(K#pi) (GeV/#it{c}); #it{M}(K#pi) (GeV/#it{c}^{2})", HistType::kTH2D, {axisPt, axisMassKaPi});
    }

    if (mlConfig.applyMlPiPiFromDplus || mlConfig.applyMlDzeroFromDstar || mlConfig.applyMlKaKaFromDsOrDplus || mlConfig.applyMlDzeroKaKaFromDstar) {
      if (doprocessPiPiFromDplusMc || doprocessKaKaFromDsOrDplusMc || doprocessKaPiFromDstarMc) {
        for (int iScore{0}; iScore < mlConfig.numMlClasses; ++iScore) {
          hBkgMlScore.push_back(registry.add<TH2>(Form("hBkgMlScore%d", iScore), Form(";#it{p}_{T}(tag) (GeV/#it{c});ML score %d; counts", iScore), HistType::kTH2D, {axisPt, axisMlScore}));
          hPromptMlScore.push_back(registry.add<TH2>(Form("hPromptMlScore%d", iScore), Form(";#it{p}_{T}(tag) (GeV/#it{c});ML score %d; counts", iScore), HistType::kTH2D, {axisPt, axisMlScore}));
          hNonPromptMlScore.push_back(registry.add<TH2>(Form("hNonPromptMlScore%d", iScore), Form(";#it{p}_{T}(tag) (GeV/#it{c});ML score %d; counts", iScore), HistType::kTH2D, {axisPt, axisMlScore}));
        }
      } else {
        for (int iScore{0}; iScore < mlConfig.numMlClasses; ++iScore) {
          hDataMlScore.push_back(registry.add<TH2>(Form("hMlScore%d", iScore), Form(";#it{p}_{T}(tag) (GeV/#it{c});ML score %d; counts", iScore), HistType::kTH2D, {axisPt, axisMlScore}));
        }
      }
    }

    const std::array<LabeledArray<double>, aod::tagandprobe::TagChannels::NTagChannels> mlCuts = {mlConfig.mlCutsPiPiFromDplus, mlConfig.mlCutsKaKaFromDsOrDplus, mlConfig.mlCutsDzeroFromDstar, mlConfig.mlCutsDzeroFromDstar, mlConfig.mlCutsDzeroKaKaFromDstar};
    const std::array<std::vector<std::string>, aod::tagandprobe::TagChannels::NTagChannels> onnxFileNames = {mlConfig.onnxFileNamesPiPiFromDplus, mlConfig.onnxFileNamesKaKaFromDsOrDplus, mlConfig.onnxFileNamesDzeroFromDstar, mlConfig.onnxFileNamesDzeroFromDstar, mlConfig.onnxFileNamesDzeroKaKaFromDstar};
    const std::array<std::vector<std::string>, aod::tagandprobe::TagChannels::NTagChannels> modelPathsCcdb = {mlConfig.modelPathsCcdbPiPiFromDplus, mlConfig.modelPathsCcdbKaKaFromDsOrDplus, mlConfig.modelPathsCcdbDzeroFromDstar, mlConfig.modelPathsCcdbDzeroFromDstar, mlConfig.modelPathsCcdbDzeroKaKaFromDstar};
    applyMl = {mlConfig.applyMlPiPiFromDplus, mlConfig.applyMlKaKaFromDsOrDplus, mlConfig.applyMlDzeroFromDstar, mlConfig.applyMlDzeroFromDstar, mlConfig.applyMlDzeroKaKaFromDstar};
    for (auto iChannel{0u}; iChannel < aod::tagandprobe::TagChannels::NTagChannels; ++iChannel) {
      if (applyMl[iChannel]) {
        mlResponse[iChannel].configure(ptBinsForTopologicalCuts[iChannel], mlCuts[iChannel], mlConfig.cutDirMl, mlConfig.numMlClasses);
        if (mlConfig.loadMlModelsFromCcdb) {
          mlResponse[iChannel].setModelPathsCCDB(onnxFileNames[iChannel], ccdbApi, modelPathsCcdb[iChannel], mlConfig.timestampCcdbForMlModels);
        } else {
          mlResponse[iChannel].setModelPathsLocal(onnxFileNames[iChannel]);
        }
        mlResponse[iChannel].init();
      }
    }
  }

  /// Fill a vector with the Mothers pdg codes
  /// \param pdgDecayMothers vector pdg codes of possible mothers
  /// \param pdgResonances vector pdg codes of possible resonanced in the decays
  /// \param channel decay channel
  void pdgMothersDecayChannel(std::vector<int>& pdgDecayMothers, std::vector<int>& pdgResonances, const uint8_t channel)
  {
    pdgDecayMothers.clear();
    pdgResonances.clear();
    if (channel == aod::tagandprobe::TagChannels::DplusToKPiPi) {
      pdgDecayMothers.push_back(constants::physics::Pdg::kDPlus);
      pdgResonances.push_back(313);    // K*(892)0
      pdgResonances.push_back(10313);  // K1(1270)0
      pdgResonances.push_back(100313); // K*(1410)0
      pdgResonances.push_back(10311);  // K*0(1430)0
      pdgResonances.push_back(100311); // K*(1460)0
      pdgResonances.push_back(20313);  // K1(1400)0
      pdgResonances.push_back(30313);  // K*(1680)0
    } else if (channel == aod::tagandprobe::TagChannels::DsOrDplusToKKPi) {
      pdgDecayMothers.push_back(constants::physics::Pdg::kDPlus);
      pdgDecayMothers.push_back(constants::physics::Pdg::kDS);
    } else if (channel == aod::tagandprobe::TagChannels::DstarPlusToDzeroPi || channel == aod::tagandprobe::TagChannels::DstarMinusToDzeroBarPi || channel == aod::tagandprobe::TagChannels::DstarToDzeroToKK) {
      pdgDecayMothers.push_back(constants::physics::Pdg::kDStar);
    }
  }

  /// Check if the given tag tracks belong to a D meson
  /// \param firstTrack candidate
  /// \param SecondTrack candidate
  /// \param mcParticles McParticles table
  /// \param channel decay channel
  /// \param pdgDecayMothers vector pdg codes of possible mothers
  /// \param pdgResonances vector pdg codes of possible resonanced in the decays
  /// \param motherIdx particle mother index
  /// \return a flag that contains the information of MC truth (see aod::tagandprobe::SignalFlags)
  template <typename PParticles, typename TTrack>
  uint8_t getTagOrigin(TTrack const& firsTrack,
                       TTrack const& secondTrack,
                       PParticles const& mcParticles,
                       const uint8_t channel,
                       std::vector<int>& pdgDecayMothers,
                       std::vector<int>& pdgResonances,
                       int& motherIdx)
  {
    int pdgTagMother{0};
    int pdgProbeParticle{-1};
    uint8_t signalFlag = 0;

    if (channel == aod::tagandprobe::TagChannels::DplusToKPiPi) {
      pdgTagMother = constants::physics::Pdg::kDPlus;
      pdgProbeParticle = 321; // Ka
    } else if (channel == aod::tagandprobe::TagChannels::DsOrDplusToKKPi) {
      pdgTagMother = constants::physics::Pdg::kPhi;
      pdgProbeParticle = 211; // Pi
    } else {
      pdgTagMother = constants::physics::Pdg::kD0;
      pdgProbeParticle = 211; // Pi
    }

    if (!firsTrack.has_mcParticle() || !secondTrack.has_mcParticle()) {
      SETBIT(signalFlag, aod::tagandprobe::SignalFlags::Bkg);
      SETBIT(signalFlag, aod::tagandprobe::SignalFlags::BkgFromNoHf);
      return signalFlag;
    } else {
      auto firstMcTrack = firsTrack.template mcParticle_as<PParticles>();
      auto secondMcTrack = secondTrack.template mcParticle_as<PParticles>();
      auto firstTrackMotherId = RecoDecay::getMother(mcParticles, firstMcTrack, pdgTagMother, true);
      auto secondTrackMotherId = RecoDecay::getMother(mcParticles, secondMcTrack, pdgTagMother, true);

      bool isTaggedAsSignal{false}, isResonant{false};
      if ((firstTrackMotherId == secondTrackMotherId) && (firstTrackMotherId != -1)) {
        auto particleMother = mcParticles.rawIteratorAt(firstTrackMotherId);

        /// π±π± for D± → K∓π±π± decays
        if (channel == aod::tagandprobe::TagChannels::DplusToKPiPi) {
          auto particleMother = mcParticles.rawIteratorAt(firstTrackMotherId);
          auto daughters = particleMother.template daughters_as<PParticles>();

          // Check if the probe is within the mother's particle daughters
          if (daughters.size() == 3) { // non-resonant decay
            for (auto& daughter : daughters) {
              if (std::abs(daughter.pdgCode()) == pdgProbeParticle) {
                isTaggedAsSignal = true;
                motherIdx = firstTrackMotherId;
                break;
              }
            }
          } else if (daughters.size() == 2) { // resonant decay
            for (auto& daughter : daughters) {
              auto absPdg = std::abs(daughter.pdgCode());
              if (std::find(pdgResonances.begin(), pdgResonances.end(), absPdg) != pdgResonances.end()) {
                isTaggedAsSignal = true;
                isResonant = true;
                motherIdx = firstTrackMotherId;
                break;
              }
            }
          }
        } else {
          ///  K∓K± for φ from Ds± or D± → φπ± decays
          ///  K∓π± for D0 from D±* → D0π± decays
          for (auto pdgGrandMother : pdgDecayMothers) {
            auto grandMotherId = RecoDecay::getMother(mcParticles, particleMother, pdgGrandMother, true);
            if (grandMotherId != -1) {
              auto particleGrandMother = mcParticles.rawIteratorAt(grandMotherId);
              auto daughters = particleGrandMother.template daughters_as<PParticles>();
              // Check if the probe is within the GrandMother's particle daughters
              if (daughters.size() == 2) { // exclude undesired decays, such as Ds± → φπ±π±π∓
                for (auto& daughter : daughters) {
                  if (std::abs(daughter.pdgCode()) == pdgProbeParticle) {
                    isTaggedAsSignal = true;
                    motherIdx = grandMotherId;
                    break;
                  }
                }
              }
            }
          }
        }
      }

      // check if it is non-prompt from beauty
      if (isTaggedAsSignal) {
        if (RecoDecay::getCharmHadronOrigin(mcParticles, mcParticles.rawIteratorAt(firstTrackMotherId)) == RecoDecay::OriginType::NonPrompt) {
          SETBIT(signalFlag, aod::tagandprobe::SignalFlags::NonPrompt);
        } else {
          SETBIT(signalFlag, aod::tagandprobe::SignalFlags::Prompt);
        }
        if (isResonant) {
          SETBIT(signalFlag, aod::tagandprobe::SignalFlags::Resonant);
        }
        return signalFlag;
      }

      // if not signal, it must be background
      SETBIT(signalFlag, aod::tagandprobe::SignalFlags::Bkg);

      auto originFirstTrack = RecoDecay::getCharmHadronOrigin(mcParticles, firstMcTrack, true);
      auto originSecondTrack = RecoDecay::getCharmHadronOrigin(mcParticles, secondMcTrack, true);
      if (originFirstTrack == RecoDecay::OriginType::None && originSecondTrack == RecoDecay::OriginType::None) {
        SETBIT(signalFlag, aod::tagandprobe::SignalFlags::BkgFromNoHf);
      }

      return signalFlag;
    }
  }

  template <typename Pvec>
  bool isSelectedInvariantMass(const Pvec& pVecTrackFirst,
                               const Pvec& pVecTrackSecond,
                               const uint8_t channel,
                               float& invMass2,
                               const int& ptBin)
  {
    auto arrMomentum = std::array{pVecTrackFirst, pVecTrackSecond};

    auto invMassMin = topologicalCuts[channel].get(ptBin, 0u);
    auto invMassMax = topologicalCuts[channel].get(ptBin, 1u);
    invMass2 = RecoDecay::m2(arrMomentum, masses[channel]);
    if (invMass2 > invMassMax * invMassMax || invMass2 < invMassMin * invMassMin) {
      return false;
    }
    return true;
  }

  template <typename TTrack>
  bool isSelectedPidTof(const TTrack& track,
                        const uint8_t channel)
  {
    if (!track.hasTOF()) { // TOF not forced anyway
      return true;
    }

    switch (channel) {
      case aod::tagandprobe::TagChannels::DplusToKPiPi: {
        if (std::abs(track.tofNSigmaPi()) < trackNumSigmaTof) {
          return true;
        }
        break;
      }
      case aod::tagandprobe::TagChannels::DsOrDplusToKKPi: {
        if (std::abs(track.tofNSigmaKa()) < trackNumSigmaTof) {
          return true;
        }
        break;
      }
      case aod::tagandprobe::TagChannels::DstarPlusToDzeroPi: {
        if ((track.signed1Pt() > 0 && std::abs(track.tofNSigmaPi()) < trackNumSigmaTof) || (track.signed1Pt() < 0 && std::abs(track.tofNSigmaKa()) < trackNumSigmaTof)) {
          return true;
        }
        break;
      }
      case aod::tagandprobe::TagChannels::DstarMinusToDzeroBarPi: {
        if ((track.signed1Pt() < 0 && std::abs(track.tofNSigmaPi()) < trackNumSigmaTof) || (track.signed1Pt() > 0 && std::abs(track.tofNSigmaKa()) < trackNumSigmaTof)) {
          return true;
        }
        break;
      }
      case aod::tagandprobe::TagChannels::DstarToDzeroToKK: {
        if (std::abs(track.tofNSigmaKa()) < trackNumSigmaTof) {
          return true;
        }
        break;
      }
    }
    return false;
  }

  template <typename PV, typename SV, typename CovMatSV, typename PVec>
  bool isSelectedTopology(const PV& primVtx,
                          const SV& secVtx,
                          const CovMatSV& covMatrixSecVtx,
                          const PVec& pVec,
                          std::array<float, 2>& trackDcaXy,
                          const uint8_t channel,
                          const int& ptBin,
                          std::vector<float>& topoVars)
  {
    topoVars.clear();
    std::array<float, 3> pvCoord = {primVtx.getX(), primVtx.getY(), primVtx.getZ()};
    auto decLen = RecoDecay::distance(pvCoord, secVtx);
    if (decLen < topologicalCuts[channel].get(ptBin, 2u)) {
      return false;
    }
    topoVars.push_back(decLen);

    auto covMatrixPV = primVtx.getCov();

    auto decLenXy = RecoDecay::distanceXY(pvCoord, secVtx);
    if (decLenXy < topologicalCuts[channel].get(ptBin, 3u)) {
      return false;
    }
    topoVars.push_back(decLenXy);

    float phi, theta;
    getPointDirection(pvCoord, secVtx, phi, theta);
    auto errorDecLen = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixSecVtx, phi, theta));
    if (decLen / errorDecLen < topologicalCuts[channel].get(ptBin, 4u)) {
      return false;
    }
    topoVars.push_back(decLen / errorDecLen);

    auto errorDecLenXy = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.f) + getRotatedCovMatrixXX(covMatrixSecVtx, phi, 0.f));
    if (decLenXy / errorDecLenXy < topologicalCuts[channel].get(ptBin, 5u)) {
      return false;
    }
    topoVars.push_back(decLenXy / errorDecLenXy);

    if (trackDcaXy[0] * trackDcaXy[1] > topologicalCuts[channel].get(ptBin, 6u)) {
      return false;
    }
    topoVars.push_back(trackDcaXy[0]);
    topoVars.push_back(trackDcaXy[1]);
    topoVars.push_back(trackDcaXy[0] * trackDcaXy[1]);

    auto cpa = RecoDecay::cpa(pvCoord, secVtx, pVec);
    if (cpa < topologicalCuts[channel].get(ptBin, 7u)) {
      return false;
    }
    topoVars.push_back(cpa);

    auto cpaXy = RecoDecay::cpaXY(pvCoord, secVtx, pVec);
    if (cpaXy < topologicalCuts[channel].get(ptBin, 8u)) {
      return false;
    }
    topoVars.push_back(cpaXy);

    return true;
  }

  template <bool doMc, typename CCollision, typename TTracks, typename PParticles>
  void computeCombinatorialSameCharge(CCollision const& collision,
                                      TTracks const& tracks, // pool of tracks
                                      const uint8_t channel,
                                      float& /*bz*/,
                                      std::vector<int>& pdgDecayMothers,
                                      std::vector<int>& pdgResonances,
                                      PParticles const& mcParticles)
  {
    for (auto trackFirst = tracks.begin(); trackFirst != tracks.end(); ++trackFirst) {

      if (applyTofPid && !isSelectedPidTof(trackFirst, channel)) {
        continue;
      }

      for (auto trackSecond = trackFirst + 1; trackSecond != tracks.end(); ++trackSecond) {

        if (applyTofPid && !isSelectedPidTof(trackSecond, channel)) {
          continue;
        }

        float invMass2{0.f};
        std::array<float, 3> pVecTrackFirst{trackFirst.px(), trackFirst.py(), trackFirst.pz()};
        std::array<float, 3> pVecTrackSecond{trackSecond.px(), trackSecond.py(), trackSecond.pz()};

        auto pVec = RecoDecay::pVec(pVecTrackFirst, pVecTrackSecond);
        auto ptTag = RecoDecay::pt(pVec);
        auto ptBin = o2::analysis::findBin(&ptBinsForTopologicalCuts[channel], ptTag);
        if (ptBin == -1) {
          continue;
        }

        if (!isSelectedInvariantMass(pVecTrackFirst, pVecTrackSecond, channel, invMass2, ptBin)) {
          continue;
        }

        auto trackParCovFirst = getTrackParCov(trackFirst);
        auto trackParCovSecond = getTrackParCov(trackSecond);

        int nVertices{0};
        try {
          nVertices = vertexer.process(trackParCovFirst, trackParCovSecond);
        } catch (...) {
          LOG(error) << "Exception caught in DCA fitter process call!";
          continue;
        }
        if (nVertices != 1) {
          continue;
        }

        auto primVtx = getPrimaryVertex(collision);
        const auto& secVtx = vertexer.getPCACandidate();
        const auto& covMatrixPCA = vertexer.calcPCACovMatrixFlat();
        std::array<float, 2> trackDcaXy{trackFirst.dcaXY(), trackSecond.dcaXY()};
        std::vector<float> topoVars{};
        if (!isSelectedTopology(primVtx, secVtx, covMatrixPCA, pVec, trackDcaXy, channel, ptBin, topoVars)) {
          continue;
        }

        uint8_t isSignal{0u};
        int motherIdx{-1};
        if constexpr (doMc) {
          isSignal = getTagOrigin(trackFirst, trackSecond, mcParticles, channel, pdgDecayMothers, pdgResonances, motherIdx);
        }

        std::vector<float> mlScoresTag{};
        if (applyMl[channel]) {
          bool isMlSelected = mlResponse[channel].isSelectedMl(topoVars, ptTag, mlScoresTag);
          // we fill control histograms
          if constexpr (doMc) {
            if (TESTBIT(isSignal, aod::tagandprobe::SignalFlags::Bkg)) {
              for (int iScore{0}; iScore < mlConfig.numMlClasses; ++iScore) {
                hBkgMlScore.at(iScore)->Fill(ptTag, mlScoresTag.at(iScore));
              }
            } else if (TESTBIT(isSignal, aod::tagandprobe::SignalFlags::Prompt)) {
              for (int iScore{0}; iScore < mlConfig.numMlClasses; ++iScore) {
                hPromptMlScore.at(iScore)->Fill(ptTag, mlScoresTag.at(iScore));
              }
            } else if (TESTBIT(isSignal, aod::tagandprobe::SignalFlags::NonPrompt)) {
              for (int iScore{0}; iScore < mlConfig.numMlClasses; ++iScore) {
                hNonPromptMlScore.at(iScore)->Fill(ptTag, mlScoresTag.at(iScore));
              }
            }
          } else {
            for (int iScore{0}; iScore < mlConfig.numMlClasses; ++iScore) {
              hDataMlScore.at(iScore)->Fill(ptTag, mlScoresTag.at(iScore));
            }
          }
          if (!isMlSelected) { // for the time being all the topological variables used for all channels (decLen, decLenXy, normDecLen, normDecLenXy, cosp, cospXy, dcaXyTrack0, dcaXyTrack1, dcaProd)
            continue;
          }
        }

        float invMass{std::sqrt(invMass2)};
        registry.fill(HIST("hMassPiPiVsPt"), ptTag, invMass); // only channel with same sign tracks for the moment

        if (fillTopoVarsTable) {
          bool fillTable{true};
          if (fillTopoVarsTable == 1 && !(TESTBIT(isSignal, aod::tagandprobe::SignalFlags::Prompt) || TESTBIT(isSignal, aod::tagandprobe::SignalFlags::NonPrompt))) { // only signal
            fillTable = false;
          } else if (fillTopoVarsTable == 2 && !TESTBIT(isSignal, aod::tagandprobe::SignalFlags::Bkg)) { // only background
            fillTable = false;
          } else if (fillTopoVarsTable == 3 && !TESTBIT(isSignal, aod::tagandprobe::SignalFlags::BkgFromNoHf)) { // only background excluding tracks from other HF decays
            fillTable = false;
          }
          float pseudoRndm = trackFirst.pt() * 1000. - static_cast<int64_t>(trackFirst.pt() * 1000);
          if (ptTag < ptTagMaxForDownsampling && pseudoRndm >= downsamplingForTopoVarTable) {
            fillTable = false;
          }
          if (fillTable) {
            tagVarsTable(ptTag, invMass, topoVars[0], topoVars[1], topoVars[2], topoVars[3], trackDcaXy[0], trackDcaXy[1], topoVars[6], topoVars[7], topoVars[8], isSignal, channel);
          }
        } else {
          float radius = std::hypot(secVtx[0], secVtx[1]);
          tagPiPiTable(trackFirst.collisionId(), trackFirst.globalIndex(), trackSecond.globalIndex(), mlScoresTag, radius);
          if constexpr (doMc) {
            tagPiPiMcTable(isSignal, motherIdx);
          }
        }
      }
    }
  }

  template <bool doMc, typename CCollision, typename TTracks, typename PParticles>
  void computeCombinatorialOppositeCharge(CCollision const& collision,
                                          TTracks const& tracksPos,
                                          TTracks const& tracksNeg,
                                          const uint8_t channel,
                                          float& /*bz*/,
                                          std::vector<int>& pdgDecayMothers,
                                          std::vector<int>& pdgResonances,
                                          PParticles const& mcParticles)
  {
    for (const auto& trackPos : tracksPos) {

      if (applyTofPid && !isSelectedPidTof(trackPos, channel)) {
        continue;
      }

      for (const auto& trackNeg : tracksNeg) {

        if (applyTofPid && !isSelectedPidTof(trackNeg, channel)) {
          continue;
        }

        float invMass2{0.f};
        std::array<float, 3> pVecTrackPos{trackPos.px(), trackPos.py(), trackPos.pz()};
        std::array<float, 3> pVecTrackNeg{trackNeg.px(), trackNeg.py(), trackNeg.pz()};

        auto pVec = RecoDecay::pVec(pVecTrackPos, pVecTrackNeg);
        auto ptTag = RecoDecay::pt(pVec);
        auto ptBin = o2::analysis::findBin(&ptBinsForTopologicalCuts[channel], ptTag);
        if (ptBin == -1) {
          continue;
        }

        if (!isSelectedInvariantMass(pVecTrackPos, pVecTrackNeg, channel, invMass2, ptBin)) {
          continue;
        }

        auto trackParCovPos = getTrackParCov(trackPos);
        auto trackParCovNeg = getTrackParCov(trackNeg);

        int nVertices{0};
        try {
          nVertices = vertexer.process(trackParCovPos, trackParCovNeg);
        } catch (...) {
          LOG(error) << "Exception caught in DCA fitter process call!";
          continue;
        }
        if (nVertices != 1) {
          continue;
        }

        auto primVtx = getPrimaryVertex(collision);
        const auto& secVtx = vertexer.getPCACandidate();
        const auto& covMatrixPCA = vertexer.calcPCACovMatrixFlat();
        std::array<float, 2> trackDcaXy{trackPos.dcaXY(), trackNeg.dcaXY()};
        std::vector<float> topoVars{};
        if (!isSelectedTopology(primVtx, secVtx, covMatrixPCA, pVec, trackDcaXy, channel, ptBin, topoVars)) {
          continue;
        }

        uint8_t isSignal{0u};
        int motherIdx{-1};
        if constexpr (doMc) {
          isSignal = getTagOrigin(trackPos, trackNeg, mcParticles, channel, pdgDecayMothers, pdgResonances, motherIdx);
        }

        std::vector<float> mlScoresTag{};
        if (applyMl[channel]) {
          bool isMlSelected = mlResponse[channel].isSelectedMl(topoVars, ptTag, mlScoresTag);
          // we fill control histograms
          if constexpr (doMc) {
            if (TESTBIT(isSignal, aod::tagandprobe::SignalFlags::Bkg)) {
              for (int iScore{0}; iScore < mlConfig.numMlClasses; ++iScore) {
                hBkgMlScore.at(iScore)->Fill(ptTag, mlScoresTag.at(iScore));
              }
            } else if (TESTBIT(isSignal, aod::tagandprobe::SignalFlags::Prompt)) {
              for (int iScore{0}; iScore < mlConfig.numMlClasses; ++iScore) {
                hPromptMlScore.at(iScore)->Fill(ptTag, mlScoresTag.at(iScore));
              }
            } else if (TESTBIT(isSignal, aod::tagandprobe::SignalFlags::NonPrompt)) {
              for (int iScore{0}; iScore < mlConfig.numMlClasses; ++iScore) {
                hNonPromptMlScore.at(iScore)->Fill(ptTag, mlScoresTag.at(iScore));
              }
            }
          } else {
            for (int iScore{0}; iScore < mlConfig.numMlClasses; ++iScore) {
              hDataMlScore.at(iScore)->Fill(ptTag, mlScoresTag.at(iScore));
            }
          }
          if (!isMlSelected) { // for the time being all the topological variables used for all channels (decLen, decLenXy, normDecLen, normDecLenXy, cosp, cospXy, dcaXyTrack0, dcaXyTrack1, dcaProd)
            continue;
          }
        }

        float invMass{std::sqrt(invMass2)};
        if (channel == aod::tagandprobe::TagChannels::DsOrDplusToKKPi) {
          registry.fill(HIST("hMassKaKaVsPt"), ptTag, invMass);
        } else if (channel == aod::tagandprobe::TagChannels::DstarPlusToDzeroPi) {
          if (!studyDzeroReflections) {
            registry.fill(HIST("hMassKaPiVsPt"), ptTag, invMass);
          } else {
            float invMassRefl{0.f};
            int isDzero = 1;
            if (std::abs(trackPos.tpcNSigmaKa()) < trackNumSigmaTpc && (std::abs(trackNeg.tpcNSigmaPi()) < trackNumSigmaTpc)) {
              isDzero = 3;
              if (applyTofPid) {
                if (!isSelectedPidTof(trackNeg, aod::tagandprobe::TagChannels::DstarMinusToDzeroBarPi) || !isSelectedPidTof(trackPos, aod::tagandprobe::TagChannels::DstarMinusToDzeroBarPi))
                  isDzero = 1;
              }
            }
            if (isDzero == 3) {
              auto arrMomentum = std::array{pVecTrackNeg, pVecTrackPos};
              invMassRefl = std::sqrt(RecoDecay::m2(arrMomentum, masses[channel]));
            }
            registry.fill(HIST("hMassKaPiVsPt"), ptTag, invMass, invMassRefl, isDzero);
          }
        } else if (channel == aod::tagandprobe::TagChannels::DstarMinusToDzeroBarPi) {
          if (!studyDzeroReflections) {
            registry.fill(HIST("hMassKaPiVsPt"), ptTag, invMass);
          } else {
            float invMassRefl{0.f};
            int isDzero = 2;
            if (std::abs(trackNeg.tpcNSigmaKa()) < trackNumSigmaTpc && (std::abs(trackPos.tpcNSigmaPi()) < trackNumSigmaTpc)) {
              isDzero = 3;
              if (applyTofPid) {
                if (!isSelectedPidTof(trackNeg, aod::tagandprobe::TagChannels::DstarPlusToDzeroPi) || !isSelectedPidTof(trackPos, aod::tagandprobe::TagChannels::DstarPlusToDzeroPi))
                  isDzero = 2;
              }
            }
            if (isDzero == 3) {
              auto arrMomentum = std::array{pVecTrackNeg, pVecTrackPos};
              invMassRefl = std::sqrt(RecoDecay::m2(arrMomentum, masses[channel]));
            }
            registry.fill(HIST("hMassKaPiVsPt"), ptTag, invMass, invMassRefl, isDzero);
          }
        } else if (channel == aod::tagandprobe::TagChannels::DstarToDzeroToKK) {
          registry.fill(HIST("hMassDzeroKaKaVsPt"), ptTag, invMass);
        }

        if (fillTopoVarsTable) {
          bool fillTable{true};
          if (fillTopoVarsTable == 1 && !(TESTBIT(isSignal, aod::tagandprobe::SignalFlags::Prompt) || TESTBIT(isSignal, aod::tagandprobe::SignalFlags::NonPrompt))) { // only signal
            fillTable = false;
          } else if (fillTopoVarsTable == 2 && !TESTBIT(isSignal, aod::tagandprobe::SignalFlags::Bkg)) { // only background
            fillTable = false;
          } else if (fillTopoVarsTable == 3 && !TESTBIT(isSignal, aod::tagandprobe::SignalFlags::BkgFromNoHf)) { // only background excluding tracks from other HF decays
            fillTable = false;
          }
          float pseudoRndm = trackPos.pt() * 1000. - static_cast<int64_t>(trackPos.pt() * 1000);
          if (ptTag < ptTagMaxForDownsampling && pseudoRndm >= downsamplingForTopoVarTable) {
            fillTable = false;
          }
          if (fillTable) {
            tagVarsTable(ptTag, invMass, topoVars[0], topoVars[1], topoVars[2], topoVars[3], trackDcaXy[0], trackDcaXy[1], topoVars[6], topoVars[7], topoVars[8], isSignal, channel);
          }
        } else {
          if (channel == aod::tagandprobe::TagChannels::DsOrDplusToKKPi) {
            float radius = std::hypot(secVtx[0], secVtx[1]);
            tagKaKaTable(trackPos.collisionId(), trackPos.globalIndex(), trackNeg.globalIndex(), mlScoresTag, radius);
            if constexpr (doMc) {
              tagKaKaMcTable(isSignal, motherIdx);
            }
          } else if (channel == aod::tagandprobe::TagChannels::DstarPlusToDzeroPi) {
            float radius = std::hypot(secVtx[0], secVtx[1]);
            tagPiKaTable(trackPos.collisionId(), trackPos.globalIndex(), trackNeg.globalIndex(), mlScoresTag, radius);
            if constexpr (doMc) {
              tagPiKaMcTable(isSignal, motherIdx);
            }
          } else if (channel == aod::tagandprobe::TagChannels::DstarMinusToDzeroBarPi) {
            float radius = std::hypot(secVtx[0], secVtx[1]);
            tagKaPiTable(trackPos.collisionId(), trackPos.globalIndex(), trackNeg.globalIndex(), mlScoresTag, radius);
            if constexpr (doMc) {
              tagKaPiMcTable(isSignal, motherIdx);
            }
          }
        }
      }
    }
  }

  void processPiPiFromDplusMc(CollisionsFiltered::iterator const& collision,
                              TracksWithSelAndDcaMcFiltered const&,
                              aod::BCsWithTimestamps const&,
                              aod::McParticles const& mcParticles)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    float bz{0};
    if (runNumber != bc.runNumber()) {
      parameters::GRPMagField* grpo = ccdb->getForTimeStamp<parameters::GRPMagField>("GLO/Config/GRPMagField", bc.timestamp());
      if (grpo != nullptr) {
        base::Propagator::initFieldFromGRP(grpo);
        bz = base::Propagator::Instance()->getNominalBz();
      } else {
        LOGF(fatal, "GRP object is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
      }
      runNumber = bc.runNumber();
    }

    std::vector<int> pdgDecayMothers{};
    std::vector<int> pdgResonances{};
    pdgMothersDecayChannel(pdgDecayMothers, pdgResonances, aod::tagandprobe::TagChannels::DplusToKPiPi);

    auto groupPositive = positivePionsMc->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    computeCombinatorialSameCharge<true>(collision, groupPositive, aod::tagandprobe::TagChannels::DplusToKPiPi, bz, pdgDecayMothers, pdgResonances, mcParticles);

    auto groupNegative = negativePionsMc->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    computeCombinatorialSameCharge<true>(collision, groupNegative, aod::tagandprobe::TagChannels::DplusToKPiPi, bz, pdgDecayMothers, pdgResonances, mcParticles);
  }
  PROCESS_SWITCH(TagTwoProngDisplacedVertices, processPiPiFromDplusMc, "Process pipi combinatorial to tag pion pairs from D+ decays Mc", false);

  void processPiPiFromDplus(CollisionsFiltered::iterator const& collision,
                            TracksWithSelAndDcaFiltered const& tracks,
                            aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    float bz{0};
    if (runNumber != bc.runNumber()) {
      parameters::GRPMagField* grpo = ccdb->getForTimeStamp<parameters::GRPMagField>("GLO/Config/GRPMagField", bc.timestamp());
      if (grpo != nullptr) {
        base::Propagator::initFieldFromGRP(grpo);
        bz = base::Propagator::Instance()->getNominalBz();
      } else {
        LOGF(fatal, "GRP object is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
      }
      runNumber = bc.runNumber();
    }

    std::vector<int> pdgDecayMothers{};
    std::vector<int> pdgResonances{};
    pdgMothersDecayChannel(pdgDecayMothers, pdgResonances, aod::tagandprobe::TagChannels::DplusToKPiPi);

    auto groupPositive = positivePions->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    computeCombinatorialSameCharge<false>(collision, groupPositive, aod::tagandprobe::TagChannels::DplusToKPiPi, bz, pdgDecayMothers, pdgResonances, tracks);

    auto groupNegative = negativePions->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    computeCombinatorialSameCharge<false>(collision, groupNegative, aod::tagandprobe::TagChannels::DplusToKPiPi, bz, pdgDecayMothers, pdgResonances, tracks);
  }
  PROCESS_SWITCH(TagTwoProngDisplacedVertices, processPiPiFromDplus, "Process pipi combinatorial to tag pion pairs from D+ decays", false);

  void processKaKaFromDsOrDplusMc(CollisionsFiltered::iterator const& collision,
                                  TracksWithSelAndDcaMcFiltered const&,
                                  aod::BCsWithTimestamps const&,
                                  aod::McParticles const& mcParticles)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    float bz{0};
    if (runNumber != bc.runNumber()) {
      parameters::GRPMagField* grpo = ccdb->getForTimeStamp<parameters::GRPMagField>("GLO/Config/GRPMagField", bc.timestamp());
      if (grpo != nullptr) {
        base::Propagator::initFieldFromGRP(grpo);
        bz = base::Propagator::Instance()->getNominalBz();
      } else {
        LOGF(fatal, "GRP object is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
      }
      runNumber = bc.runNumber();
    }

    std::vector<int> pdgDecayMothers{};
    std::vector<int> pdgResonances{};
    pdgMothersDecayChannel(pdgDecayMothers, pdgResonances, aod::tagandprobe::TagChannels::DsOrDplusToKKPi);

    auto groupPositive = positiveKaonsMc->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto groupNegative = negativeKaonsMc->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    computeCombinatorialOppositeCharge<true>(collision, groupPositive, groupNegative, aod::tagandprobe::TagChannels::DsOrDplusToKKPi, bz, pdgDecayMothers, pdgResonances, mcParticles);
  }
  PROCESS_SWITCH(TagTwoProngDisplacedVertices, processKaKaFromDsOrDplusMc, "Process KK combinatorial to tag kaon pairs from Ds+/D+ decays Mc", false);

  void processKaKaFromDsOrDplus(CollisionsFiltered::iterator const& collision,
                                TracksWithSelAndDcaFiltered const& tracks,
                                aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    float bz{0};
    if (runNumber != bc.runNumber()) {
      parameters::GRPMagField* grpo = ccdb->getForTimeStamp<parameters::GRPMagField>("GLO/Config/GRPMagField", bc.timestamp());
      if (grpo != nullptr) {
        base::Propagator::initFieldFromGRP(grpo);
        bz = base::Propagator::Instance()->getNominalBz();
      } else {
        LOGF(fatal, "GRP object is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
      }
      runNumber = bc.runNumber();
    }

    std::vector<int> pdgDecayMothers{};
    std::vector<int> pdgResonances{};
    pdgMothersDecayChannel(pdgDecayMothers, pdgResonances, aod::tagandprobe::TagChannels::DsOrDplusToKKPi);

    auto groupPositive = positiveKaons->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto groupNegative = negativeKaons->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    computeCombinatorialOppositeCharge<false>(collision, groupPositive, groupNegative, aod::tagandprobe::TagChannels::DsOrDplusToKKPi, bz, pdgDecayMothers, pdgResonances, tracks);
  }
  PROCESS_SWITCH(TagTwoProngDisplacedVertices, processKaKaFromDsOrDplus, "Process KK combinatorial to tag kaon pairs from Ds+/D+ decays", false);

  void processKaKaFromDzero(CollisionsFiltered::iterator const& collision,
                            TracksWithSelAndDcaFiltered const& tracks,
                            aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    float bz{0};
    if (runNumber != bc.runNumber()) {
      parameters::GRPMagField* grpo = ccdb->getForTimeStamp<parameters::GRPMagField>("GLO/Config/GRPMagField", bc.timestamp());
      if (grpo != nullptr) {
        base::Propagator::initFieldFromGRP(grpo);
        bz = base::Propagator::Instance()->getNominalBz();
      } else {
        LOGF(fatal, "GRP object is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
      }
      runNumber = bc.runNumber();
    }

    std::vector<int> pdgDecayMothers{};
    std::vector<int> pdgResonances{};
    pdgMothersDecayChannel(pdgDecayMothers, pdgResonances, aod::tagandprobe::TagChannels::DstarToDzeroToKK);

    auto groupPositive = positiveKaons->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto groupNegative = negativeKaons->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    computeCombinatorialOppositeCharge<false>(collision, groupPositive, groupNegative, aod::tagandprobe::TagChannels::DstarToDzeroToKK, bz, pdgDecayMothers, pdgResonances, tracks);
  }
  PROCESS_SWITCH(TagTwoProngDisplacedVertices, processKaKaFromDzero, "Process KK combinatorial to tag kaon pairs from Dzero decays", false);

  void processKaPiFromDstar(CollisionsFiltered::iterator const& collision,
                            TracksWithSelAndDcaFiltered const& tracks,
                            aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    float bz{0};
    if (runNumber != bc.runNumber()) {
      parameters::GRPMagField* grpo = ccdb->getForTimeStamp<parameters::GRPMagField>("GLO/Config/GRPMagField", bc.timestamp());
      if (grpo != nullptr) {
        base::Propagator::initFieldFromGRP(grpo);
        bz = base::Propagator::Instance()->getNominalBz();
      } else {
        LOGF(fatal, "GRP object is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
      }
      runNumber = bc.runNumber();
    }

    std::vector<int> pdgDecayMothers{};
    std::vector<int> pdgResonances{};
    pdgMothersDecayChannel(pdgDecayMothers, pdgResonances, aod::tagandprobe::TagChannels::DstarPlusToDzeroPi);

    auto groupPionPositive = positivePions->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto groupPionNegative = negativePions->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto groupKaonPositive = positiveKaons->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto groupKaonNegative = negativeKaons->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    computeCombinatorialOppositeCharge<false>(collision, groupPionPositive, groupKaonNegative, aod::tagandprobe::TagChannels::DstarPlusToDzeroPi, bz, pdgDecayMothers, pdgResonances, tracks);
    computeCombinatorialOppositeCharge<false>(collision, groupKaonPositive, groupPionNegative, aod::tagandprobe::TagChannels::DstarMinusToDzeroBarPi, bz, pdgDecayMothers, pdgResonances, tracks);
  }
  PROCESS_SWITCH(TagTwoProngDisplacedVertices, processKaPiFromDstar, "Process Kpi combinatorial to tag D0 from D*+ decays", false);

  void processKaPiFromDstarMc(CollisionsFiltered::iterator const& collision,
                              TracksWithSelAndDcaMcFiltered const&,
                              aod::BCsWithTimestamps const&,
                              aod::McParticles const& mcParticles)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    float bz{0};
    if (runNumber != bc.runNumber()) {
      parameters::GRPMagField* grpo = ccdb->getForTimeStamp<parameters::GRPMagField>("GLO/Config/GRPMagField", bc.timestamp());
      if (grpo != nullptr) {
        base::Propagator::initFieldFromGRP(grpo);
        bz = base::Propagator::Instance()->getNominalBz();
      } else {
        LOGF(fatal, "GRP object is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
      }
      runNumber = bc.runNumber();
    }

    std::vector<int> pdgDecayMothers{};
    std::vector<int> pdgResonances{};
    pdgMothersDecayChannel(pdgDecayMothers, pdgResonances, aod::tagandprobe::TagChannels::DstarPlusToDzeroPi);

    auto groupPionPositive = positivePionsMc->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto groupPionNegative = negativePionsMc->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto groupKaonPositive = positiveKaonsMc->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto groupKaonNegative = negativeKaonsMc->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    computeCombinatorialOppositeCharge<true>(collision, groupPionPositive, groupKaonNegative, aod::tagandprobe::TagChannels::DstarPlusToDzeroPi, bz, pdgDecayMothers, pdgResonances, mcParticles);
    computeCombinatorialOppositeCharge<true>(collision, groupKaonPositive, groupPionNegative, aod::tagandprobe::TagChannels::DstarMinusToDzeroBarPi, bz, pdgDecayMothers, pdgResonances, mcParticles);
  }
  PROCESS_SWITCH(TagTwoProngDisplacedVertices, processKaPiFromDstarMc, "Process Kpi combinatorial to tag D0 from D*+ decays", false);
};

/// Probe third track reconstruction efficiency with different selections
struct ProbeThirdTrack {

  enum EvSelITSLayers : uint8_t {
    None = 0,
    GoodITSLayer3,
    GoodITSLayer0123,
    GoodITSLayersAll,
    NEvSelITSLayers
  };

  // ML models for triggers
  struct : ConfigurableGroup {
    std::string prefix = "ML";
    Configurable<bool> applyMlPiPiFromDplus{"applyMlPiPiFromDplus", false, "Flag to enable ML application for pipi pairs from D+ decays"};
    Configurable<bool> applyMlKaKaFromDsOrDplus{"applyMlKaKaFromDsOrDplus", false, "Flag to enable ML application for KK pairs from Ds or D+ decays"};
    Configurable<bool> applyMlDzeroFromDstar{"applyMlDzeroFromDstar", false, "Flag to enable ML application for Kpi pairs from D0 <- D*+ decays"};
    // pt bins
    Configurable<std::vector<double>> binsPtPiPiFromDplus{"binsPtPiPiFromDplus", std::vector<double>{aod::tagandprobe::vecBinsPt}, "pT bin limits for pipi pairs from D+ decays"};
    Configurable<std::vector<double>> binsPtKaKaFromDsOrDplus{"binsPtKaKaFromDsOrDplus", std::vector<double>{aod::tagandprobe::vecBinsPt}, "pT bin limits for KK pairs from Ds or D+ decays"};
    Configurable<std::vector<double>> binsPtDzeroFromDstar{"binsPtDzeroFromDstar", std::vector<double>{aod::tagandprobe::vecBinsPt}, "pT bin limits for Kpi pairs from D0 <- D*+ decays"};
    // ML cuts
    Configurable<LabeledArray<double>> mlCutsPiPiFromDplus{"mlCutsPiPiFromDplus", {aod::tagandprobe::mlCuts[0], aod::tagandprobe::nBinsPt, 3, aod::tagandprobe::labelsEmpty, aod::tagandprobe::labelsMlScores}, "ML Selections for pipi pairs from D+ decays"};
    Configurable<LabeledArray<double>> mlCutsKaKaFromDsOrDplus{"mlCutsKaKaFromDsOrDplus", {aod::tagandprobe::mlCuts[0], aod::tagandprobe::nBinsPt, 3, aod::tagandprobe::labelsEmpty, aod::tagandprobe::labelsMlScores}, "ML Selections for KK pairs from Ds or D+ decays"};
    Configurable<LabeledArray<double>> mlCutsDzeroFromDstar{"mlCutsDzeroFromDstar", {aod::tagandprobe::mlCuts[0], aod::tagandprobe::nBinsPt, 3, aod::tagandprobe::labelsEmpty, aod::tagandprobe::labelsMlScores}, "ML Selections for Kpi pairs from D0 <- D*+ decays"};
  } mlConfig;
  Configurable<float> ptCandMin{"ptCandMin", 0.f, "Minimum candidate pt for THnSparse filling"};
  Configurable<bool> fillTpcOnlyCase{"fillTpcOnlyCase", true, "Fill output for TPC only case (not needed for thinned data or Pb-Pb)"};
  Configurable<uint8_t> requireCollisionsGoodITS{"requireCollisionsGoodITS", 0u, "Event selection for ITS full acceptance (0: none, 1: layer 3, 2: layers 0-1-2-3, 3: all layers)"};

  ConfigurableAxis axisPtProbe{"axisPtProbe", {VARIABLE_WIDTH, 0.05f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.2f, 1.5f, 2.0f, 2.5f, 3.0f, 3.5f, 4.0f, 4.5f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f, 10.f, 12.f, 15.f, 20.f, 25.f, 30.f}, "Axis for pt Probe"};
  ConfigurableAxis axisPtTag{"axisPtTag", {VARIABLE_WIDTH, 0.05f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.2f, 1.5f, 2.0f, 2.5f, 3.0f, 3.5f, 4.0f, 4.5f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f, 10.f, 12.f, 15.f, 20.f, 25.f, 30.f}, "Axis for pt Tag"};
  ConfigurableAxis axisPtD{"axisPtD", {VARIABLE_WIDTH, 0.f, 0.5f, 1.f, 1.5f, 2.0f, 2.5f, 3.0f, 3.5f, 4.0f, 4.5f, 5.0f, 5.5f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 8.5f, 9.0f, 9.5f, 10.f, 11.f, 12.f, 14.f, 16.f, 20.f, 24.f, 36.f, 50.f}, "Axis for pt D"};
  ConfigurableAxis axisYD{"axisYD", {20, -1.f, 1.f}, "Axis for YD"};
  ConfigurableAxis axisEtaProbe{"axisEtaProbe", {20, -1.f, 1.f}, "Axis for Eta Probe"};
  ConfigurableAxis axisNumCrossRowTpc{"axisNumCrossRowTpc", {51, 49.5f, 100.5f}, "Axis for Number of CrossRowTpc"};
  ConfigurableAxis axisTpcChi2PerClus{"axisTpcChi2PerClus", {8, 2.f, 10.f}, "Axis for TpcChi2 Per Cluster"};
  ConfigurableAxis axisNumCluIts{"axisNumCluIts", {5, 2.5f, 7.5f}, "Axis for Number of Cluster ITS"};
  ConfigurableAxis axisPtMinTagdaught{"axisPtMinTagdaught", {10, 0.f, 1.f}, "Axis for Pt Min of Tag daughter"};
  ConfigurableAxis axisAbsEtaMaxTagdaught{"axisAbsEtaMaxTagdaught", {10, 0.f, 1.f}, "Axis for AbsEtaMax for Tag daughter"};
  ConfigurableAxis axisRadiusTag{"axisRadiusTag", {100, 0.f, 0.5f}, "Axis for Tag Radius (cm)"};

  Filter tagMcFilter = aod::tagandprobe::isSignal > static_cast<uint8_t>(0);

  using TracksWithDca = soa::Join<aod::Tracks, aod::TracksDCA, aod::TracksExtra>;
  using TracksWithDcaMc = soa::Join<TracksWithDca, aod::McTrackLabels>;
  using CollisionsWithEvSel = soa::Join<aod::Collisions, aod::EvSels>;
  using FilteredPiPiFromDpMcTags = soa::Filtered<soa::Join<aod::PiPiFromDpTags, aod::PiPiFromDpMcTags>>;
  using FilteredKaKaFromDspMcTags = soa::Filtered<soa::Join<aod::KaKaFromDspTags, aod::KaKaFromDsMcTags>>;
  using FilteredPiKaFromDzMcTags = soa::Filtered<soa::Join<aod::PiKaFromDzTags, aod::PiKaFromDzMcTags>>;
  using FilteredKaPiFromDzMcTags = soa::Filtered<soa::Join<aod::KaPiFromDzTags, aod::KaPiFromDzMcTags>>;

  Preslice<aod::PiPiFromDpTags> tagsPiPiPerCollision = aod::tagandprobe::collisionId;
  Preslice<aod::KaKaFromDspTags> tagsKaKaPerCollision = aod::tagandprobe::collisionId;
  Preslice<aod::PiKaFromDzTags> tagsPiKaPerCollision = aod::tagandprobe::collisionId;
  Preslice<aod::KaPiFromDzTags> tagsKaPiPerCollision = aod::tagandprobe::collisionId;
  Preslice<FilteredPiPiFromDpMcTags> tagsPiPiMcPerCollision = aod::tagandprobe::collisionId;
  Preslice<FilteredKaKaFromDspMcTags> tagsKaKaMcPerCollision = aod::tagandprobe::collisionId;
  Preslice<FilteredPiKaFromDzMcTags> tagsPiKaMcPerCollision = aod::tagandprobe::collisionId;
  Preslice<FilteredKaPiFromDzMcTags> tagsKaPiMcPerCollision = aod::tagandprobe::collisionId;
  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;

  std::array<std::array<double, 3>, aod::tagandprobe::TagChannels::NTagChannels> masses = {std::array{constants::physics::MassPionCharged, constants::physics::MassPionCharged, constants::physics::MassKaonCharged},
                                                                                           std::array{constants::physics::MassKaonCharged, constants::physics::MassKaonCharged, constants::physics::MassPionCharged},
                                                                                           std::array{constants::physics::MassPionCharged, constants::physics::MassKaonCharged, constants::physics::MassPionCharged},
                                                                                           std::array{constants::physics::MassKaonCharged, constants::physics::MassPionCharged, constants::physics::MassPionCharged},
                                                                                           std::array{constants::physics::MassKaonCharged, constants::physics::MassKaonCharged, constants::physics::MassPionCharged}};

  std::array<TrackSelection, aod::tagandprobe::TrackTypes::NTrackTypes> trackSelector{}; // define the track selectors
  std::array<bool, aod::tagandprobe::TagChannels::NTagChannels> applyMl{};
  std::array<float, aod::tagandprobe::TagChannels::NTagChannels> minInvMass{};
  std::array<float, aod::tagandprobe::TagChannels::NTagChannels> maxInvMass{};

  std::array<std::array<std::shared_ptr<THnSparse>, aod::tagandprobe::TrackTypes::NTrackTypes>, aod::tagandprobe::TagChannels::NTagChannels> histos{};
  std::array<std::shared_ptr<THnSparse>, aod::tagandprobe::TagChannels::NTagChannels> histosGen{};
  HistogramRegistry registry{"registry"};

  void init(InitContext&)
  {
    if ((doprocessCombinatorialDplusToKaPiPi && doprocessCombinatorialDplusToKaPiPiMc) || (doprocessCombinatorialDsToPhiPi && doprocessCombinatorialDsToPhiPiMc) || (doprocessCombinatorialDstarToDzeroPi && doprocessCombinatorialDstarToDzeroPiMc)) {
      LOGP(fatal, "The process functions for the same channel with and without MC truth cannot be enabled at the same time! Please check your configuration");
    }

    // ITS-TPC tracks (global tracks)
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDca].SetTrackType(o2::aod::track::TrackTypeEnum::Track);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDca].SetPtRange(0.05f, 1e10f);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDca].SetEtaRange(-1.f, 1.f);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDca].SetRequireITSRefit(true);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDca].SetRequireTPCRefit(true);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDca].SetMinNCrossedRowsTPC(50);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDca].SetMaxChi2PerClusterTPC(10.f);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDca].SetRequireHitsInITSLayers(1, {0, 1, 2}); // one hit in any IB layer
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDca].SetMaxChi2PerClusterITS(36.f);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDca].SetMaxDcaZ(2.f);

    // TPC tracks (global tracks without ITS IB requirement)
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoItsIb].SetTrackType(o2::aod::track::TrackTypeEnum::Track);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoItsIb].SetPtRange(0.05f, 1e10f);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoItsIb].SetEtaRange(-1.f, 1.f);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoItsIb].SetRequireTPCRefit(true);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoItsIb].SetMinNCrossedRowsTPC(50);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoItsIb].SetMaxChi2PerClusterTPC(10.f);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoItsIb].SetRequireHitsInITSLayers(3, {0, 1, 2, 3, 4, 5, 6}); // at least three hits in whatever layer
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoItsIb].SetMaxChi2PerClusterITS(36.f);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoItsIb].SetMaxDcaZ(2.f);

    // TPC tracks (global tracks without ITS requirements)
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoIts].SetTrackType(o2::aod::track::TrackTypeEnum::Track);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoIts].SetPtRange(0.05f, 1e10f);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoIts].SetEtaRange(-1.f, 1.f);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoIts].SetRequireTPCRefit(true);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoIts].SetMinNCrossedRowsTPC(50);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoIts].SetMaxChi2PerClusterTPC(10.f);

    // ITS tracks (global tracks without TPC requirements)
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoTpc].SetTrackType(o2::aod::track::TrackTypeEnum::Track);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoTpc].SetPtRange(0.05f, 1e10f);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoTpc].SetEtaRange(-1.f, 1.f);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoTpc].SetRequireITSRefit(true);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoTpc].SetRequireHitsInITSLayers(1, {0, 1, 2}); // one hit in any SPD layer
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoTpc].SetMaxChi2PerClusterITS(36.f);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoTpc].SetMaxDcaZ(2.f);

    std::array<AxisSpec, aod::tagandprobe::TagChannels::NTagChannels> axisMass = {AxisSpec{225, 1.65f, 2.10f}, AxisSpec{225, 1.65f, 2.10f}, AxisSpec{350, 0.135f, 0.17f}, AxisSpec{350, 0.135f, 0.17f}, AxisSpec{350, 0.135f, 0.17f}};
    std::array<AxisSpec, aod::tagandprobe::TagChannels::NTagChannels> axisMassTag = {AxisSpec{125, 0.f, 2.5f}, AxisSpec{100, constants::physics::MassPhi - 0.05f, constants::physics::MassPhi + 0.05f}, AxisSpec{200, constants::physics::MassD0 - 0.2f, constants::physics::MassD0 + 0.2f}, AxisSpec{200, constants::physics::MassD0 - 0.2f, constants::physics::MassD0 + 0.2f}, AxisSpec{200, constants::physics::MassD0 - 0.2f, constants::physics::MassD0 + 0.2f}};

    std::string trackTypes[aod::tagandprobe::TrackTypes::NTrackTypes] = {"ItsTpc", "ItsTpcNoIb", "Tpc", "Its"};
    std::string tagChannels[aod::tagandprobe::TagChannels::NTagChannels] = {"DplusToKPiPi", "DsOrDplusToKKPi", "DstarPlusToDzeroPi", "DstarMinusToDzeroBarPi", "DstarChargedToDzeroToKK"};

    for (int iChannel{0}; iChannel < aod::tagandprobe::TagChannels::NTagChannels; ++iChannel) {
      for (int iTrackType{0}; iTrackType < aod::tagandprobe::TrackTypes::NTrackTypes; ++iTrackType) {
        if (!fillTpcOnlyCase && iTrackType == aod::tagandprobe::TrackTypes::GlobalWoDcaWoIts) {
          continue;
        }
        histos[iChannel][iTrackType] = registry.add<THnSparse>(Form("h%sVsPtProbeTag_%s", tagChannels[iChannel].data(), trackTypes[iTrackType].data()),
                                                               "; #it{p}_{T}(D) (GeV/#it{c}); #it{p}_{T}(tag) (GeV/#it{c}); #it{p}_{T}(probe) (GeV/#it{c}); #it{p}_{T}^{TPC in}(probe) (GeV/#it{c}); #it{M}(D) (GeV/#it{c}^{2}); #it{M}(tag) (GeV/#it{c}^{2}); #it{#eta}(probe); #it{N}_{cross rows}^{TPC}(probe); #chi^{2}/#it{N}_{clusters}^{TPC}(probe); #it{N}_{clusters}^{ITS}(probe); radius(tag) (cm)",
                                                               HistType::kTHnSparseF, {axisPtD, axisPtTag, axisPtProbe, axisPtProbe, axisMass[iChannel], axisMassTag[iChannel], axisEtaProbe, axisNumCrossRowTpc, axisTpcChi2PerClus, axisNumCluIts, axisRadiusTag});
        auto invMassBins = axisMass[iChannel].binEdges;
        minInvMass[iChannel] = invMassBins.front();
        maxInvMass[iChannel] = invMassBins.back();
      }
    }
    for (int iChannel{0}; iChannel < aod::tagandprobe::TagChannels::NTagChannels; ++iChannel) {
      histosGen[iChannel] = registry.add<THnSparse>(Form("hGen%s", tagChannels[iChannel].data()), ";#it{p}_{T}(D_{parent}) (GeV/#it{c}); #it{y}(D_{parent});#it{p}_{T}(tag) (GeV/#it{c}); #it{y}(tag);#it{p}_{T}(probe) (GeV/#it{c}); #it{#eta}(probe);#it{p}_{T}^{min}(tag daughters);#it{#eta}_{max}(tag daughters); radius(tag) (cm)", HistType::kTHnSparseF, {axisPtD, axisYD, axisPtTag, axisYD, axisPtProbe, axisEtaProbe, axisPtMinTagdaught, axisAbsEtaMaxTagdaught, axisRadiusTag});
    }
    applyMl = {mlConfig.applyMlPiPiFromDplus, mlConfig.applyMlKaKaFromDsOrDplus, mlConfig.applyMlDzeroFromDstar};
  }

  template <typename TTrack>
  void computeInvariantMass(TTrack const& trackFirst, TTrack const& trackSecond, TTrack const& trackThird, const uint8_t channel, float& ptTag, float& invMassTag, float& ptD, float& invMass)
  {
    std::array<float, 3> pVecTrackFirst{trackFirst.px(), trackFirst.py(), trackFirst.pz()};
    std::array<float, 3> pVecTrackSecond{trackSecond.px(), trackSecond.py(), trackSecond.pz()};
    std::array<float, 3> pVecTrackThird{trackThird.px(), trackThird.py(), trackThird.pz()};
    auto arrMomentum = std::array{pVecTrackFirst, pVecTrackSecond, pVecTrackThird};
    auto arrMomentumTag = std::array{pVecTrackFirst, pVecTrackSecond};
    ptTag = RecoDecay::pt(RecoDecay::pVec(pVecTrackFirst, pVecTrackSecond));
    ptD = RecoDecay::pt(RecoDecay::pVec(pVecTrackFirst, pVecTrackSecond, pVecTrackThird));
    invMass = RecoDecay::m(arrMomentum, masses[channel]);
    auto massesTagDau = std::array{masses[channel][0], masses[channel][1]};
    invMassTag = RecoDecay::m(arrMomentumTag, massesTagDau);

    if (channel == aod::tagandprobe::TagChannels::DstarPlusToDzeroPi || channel == aod::tagandprobe::TagChannels::DstarMinusToDzeroBarPi) {
      invMass -= invMassTag;
    }
  }

  template <uint8_t channel, bool doMc, typename TTrackIndices, typename TTrack, typename TTracks, typename PParticles>
  void loopOverThirdTrack(TTrackIndices const& groupedTrackThirdIndices, TTracks const& /*tracks*/, TTrack const& trackFirst, TTrack const& trackSecond, PParticles const mcParticles, const int motherIdxTag, const float radius)
  {
    for (const auto& trackIndex : groupedTrackThirdIndices) {
      auto trackThird = trackIndex.template track_as<TTracks>();
      if constexpr (doMc) {
        if (!trackThird.has_mcParticle()) {
          continue;
        }
        int motherIdxProbe{-1};
        int pdgMotherFirst{0}, pdgMotherSecond{0};
        if constexpr (channel == aod::tagandprobe::TagChannels::DplusToKPiPi) {
          pdgMotherFirst = constants::physics::Pdg::kDPlus;
        } else if constexpr (channel == aod::tagandprobe::TagChannels::DsOrDplusToKKPi) {
          pdgMotherFirst = constants::physics::Pdg::kDS;
          pdgMotherSecond = constants::physics::Pdg::kDPlus;
        } else if constexpr (channel == aod::tagandprobe::TagChannels::DstarPlusToDzeroPi || channel == aod::tagandprobe::TagChannels::DstarMinusToDzeroBarPi) {
          pdgMotherFirst = constants::physics::Pdg::kDStar;
        }
        auto particleProbe = trackThird.template mcParticle_as<PParticles>();
        motherIdxProbe = RecoDecay::getMother(mcParticles, particleProbe, pdgMotherFirst, true);
        if constexpr (channel == aod::tagandprobe::TagChannels::DsOrDplusToKKPi) {
          if (motherIdxProbe < 0) {
            motherIdxProbe = RecoDecay::getMother(mcParticles, particleProbe, pdgMotherSecond, true);
          }
        }
        if (motherIdxProbe < 0 || motherIdxTag != motherIdxProbe) {
          continue;
        }
      }

      if (trackThird.globalIndex() == trackFirst.globalIndex() || trackThird.globalIndex() == trackSecond.globalIndex()) {
        continue;
      }
      if constexpr (channel == aod::tagandprobe::TagChannels::DplusToKPiPi) { // must be opposite sign
        if (trackThird.signed1Pt() * trackFirst.signed1Pt() > 0.) {
          continue;
        }
      } else if constexpr (channel == aod::tagandprobe::TagChannels::DstarPlusToDzeroPi) { // must be positive
        if (trackThird.signed1Pt() < 0.) {
          continue;
        }
      } else if constexpr (channel == aod::tagandprobe::TagChannels::DstarMinusToDzeroBarPi) { // must be negative
        if (trackThird.signed1Pt() > 0.) {
          continue;
        }
      }
      auto ptTrackThird = trackThird.pt();
      auto ptTpcInnerTrackThird = trackThird.tpcInnerParam() / std::sqrt(1.f + trackThird.tgl() * trackThird.tgl());
      auto etaTrackThird = trackThird.eta();
      auto numTpcCrossRowTrackThird = trackThird.tpcNClsCrossedRows();
      auto numTpcChi2NumCluTrackThird = trackThird.tpcChi2NCl();
      auto numItsCluTrackThird = trackThird.itsNCls();
      float invMass{-1.f}, invMassTag{-1.f}, ptTag{-1.f}, ptD{-1.f};
      computeInvariantMass(trackFirst, trackSecond, trackThird, channel, ptTag, invMassTag, ptD, invMass);
      if (invMass < minInvMass[channel] || invMass > maxInvMass[channel]) {
        continue;
      }
      if (ptD < ptCandMin) {
        /// candidate pt lower than the minimum allowed value, let's skip it
        continue;
      }
      for (int iTrackType{0}; iTrackType < aod::tagandprobe::TrackTypes::NTrackTypes; ++iTrackType) {
        if (!fillTpcOnlyCase && iTrackType == aod::tagandprobe::TrackTypes::GlobalWoDcaWoIts) {
          continue;
        }
        if (trackSelector[iTrackType].IsSelected(trackThird)) {
          histos[channel][iTrackType]->Fill(ptD, ptTag, ptTrackThird, ptTpcInnerTrackThird, invMass, invMassTag, etaTrackThird, numTpcCrossRowTrackThird, numTpcChi2NumCluTrackThird, numItsCluTrackThird, radius);
        }
      }
    }
  }

  bool isCollisionSelected(CollisionsWithEvSel::iterator const& collision)
  {
    switch (requireCollisionsGoodITS.value) {
      case None:
        return true;
        break;
      case GoodITSLayer3:
        return collision.selection_bit(o2::aod::evsel::kIsGoodITSLayer3);
        break;
      case GoodITSLayer0123:
        return collision.selection_bit(o2::aod::evsel::kIsGoodITSLayer0123);
        break;
      case GoodITSLayersAll:
        return collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll);
        break;
      default:
        LOGP(fatal, "Event selection flag for ITS acceptance {} not properly set.", requireCollisionsGoodITS.value);
        break;
    }
  }

  template <uint8_t channel, bool doMc, typename TTags, typename TTrackIndices, typename TTracks, typename PParticles>
  void runCombinatorialThirdTrack(TTags const& groupedTags,
                                  TTrackIndices const& groupedTrackIndices,
                                  TTracks const& tracks,
                                  PParticles const& mcParticles)
  {
    for (const auto& tag : groupedTags) {
      auto trackFirst = tag.template track0_as<TTracks>();
      auto trackSecond = tag.template track1_as<TTracks>();
      if (applyMl[channel] && tag.mlScores().size() == 3) {
        std::array<float, 3> pVecTrackFirst{trackFirst.px(), trackFirst.py(), trackFirst.pz()};
        std::array<float, 3> pVecTrackSecond{trackSecond.px(), trackSecond.py(), trackSecond.pz()};
        auto ptTag = RecoDecay::pt(RecoDecay::pVec(pVecTrackFirst, pVecTrackSecond));
        auto ptBin = o2::analysis::findBin(&mlConfig.binsPtPiPiFromDplus.value, ptTag);
        if (tag.mlScores()[0] > mlConfig.mlCutsPiPiFromDplus->get(ptBin, 0u) || tag.mlScores()[1] < mlConfig.mlCutsPiPiFromDplus->get(ptBin, 1u) || tag.mlScores()[2] < mlConfig.mlCutsPiPiFromDplus->get(ptBin, 2u)) {
          continue;
        }
      }
      int motherIdxTag{-1};
      if constexpr (doMc) {
        motherIdxTag = tag.motherId();
      }
      loopOverThirdTrack<channel, doMc>(groupedTrackIndices, tracks, trackFirst, trackSecond, mcParticles, motherIdxTag, tag.radius());
    }
  }

  void processCombinatorialDplusToKaPiPi(CollisionsWithEvSel const& collisions,
                                         aod::PiPiFromDpTags const& tagsPiPi,
                                         aod::TrackAssoc const& trackIndices,
                                         TracksWithDca const& tracks)
  {
    for (const auto& collision : collisions) {
      if (!isCollisionSelected(collision)) {
        continue;
      }
      auto thisCollId = collision.globalIndex();
      auto groupedTrackIndices = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      // D+ -> pi+pi+K- and c.c.
      auto groupedTagsPiPi = tagsPiPi.sliceBy(tagsPiPiPerCollision, thisCollId);
      runCombinatorialThirdTrack<aod::tagandprobe::TagChannels::DplusToKPiPi, false>(groupedTagsPiPi, groupedTrackIndices, tracks, tracks);
    }
  }
  PROCESS_SWITCH(ProbeThirdTrack, processCombinatorialDplusToKaPiPi, "Process combinatorial of tagged 2-pion vertices with additional track", true);

  void processCombinatorialDplusToKaPiPiMc(CollisionsWithEvSel const& collisions,
                                           FilteredPiPiFromDpMcTags const& tagsPiPi,
                                           aod::TrackAssoc const& trackIndices,
                                           TracksWithDcaMc const& tracks,
                                           aod::McParticles const& mcParticles)
  {
    for (const auto& collision : collisions) {
      if (!isCollisionSelected(collision)) {
        continue;
      }
      auto thisCollId = collision.globalIndex();
      auto groupedTrackIndices = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      // D+ -> pi+pi+K- and c.c.
      auto groupedTagsPiPi = tagsPiPi.sliceBy(tagsPiPiMcPerCollision, thisCollId);
      runCombinatorialThirdTrack<aod::tagandprobe::TagChannels::DplusToKPiPi, true>(groupedTagsPiPi, groupedTrackIndices, tracks, mcParticles);
    }
  }
  PROCESS_SWITCH(ProbeThirdTrack, processCombinatorialDplusToKaPiPiMc, "Process combinatorial of tagged 2-pion vertices with additional track using MC truth", false);

  void processCombinatorialDsToPhiPi(CollisionsWithEvSel const& collisions,
                                     aod::KaKaFromDspTags const& tagsKaKa,
                                     aod::TrackAssoc const& trackIndices,
                                     TracksWithDca const& tracks)
  {
    for (const auto& collision : collisions) {
      if (!isCollisionSelected(collision)) {
        continue;
      }
      auto thisCollId = collision.globalIndex();
      // Ds+/D+ -> phi(->K+K-)pi+ and c.c.
      auto groupedTrackIndices = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      auto groupedTagsKaKa = tagsKaKa.sliceBy(tagsKaKaPerCollision, thisCollId);
      runCombinatorialThirdTrack<aod::tagandprobe::TagChannels::DsOrDplusToKKPi, false>(groupedTagsKaKa, groupedTrackIndices, tracks, tracks);
    }
  }
  PROCESS_SWITCH(ProbeThirdTrack, processCombinatorialDsToPhiPi, "Process combinatorial of tagged 2-kaon (phi) vertices with additional track", true);

  void processCombinatorialDsToPhiPiMc(CollisionsWithEvSel const& collisions,
                                       FilteredKaKaFromDspMcTags const& tagsKaKa,
                                       aod::TrackAssoc const& trackIndices,
                                       TracksWithDcaMc const& tracks,
                                       aod::McParticles const& mcParticles)
  {
    for (const auto& collision : collisions) {
      if (!isCollisionSelected(collision)) {
        continue;
      }
      auto thisCollId = collision.globalIndex();
      // Ds+/D+ -> phi(->K+K-)pi+ and c.c.
      auto groupedTrackIndices = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      auto groupedTagsKaKa = tagsKaKa.sliceBy(tagsKaKaMcPerCollision, thisCollId);
      runCombinatorialThirdTrack<aod::tagandprobe::TagChannels::DsOrDplusToKKPi, true>(groupedTagsKaKa, groupedTrackIndices, tracks, mcParticles);
    }
  }
  PROCESS_SWITCH(ProbeThirdTrack, processCombinatorialDsToPhiPiMc, "Process combinatorial of tagged 2-kaon (phi) vertices with additional track using MC truth", false);

  void processCombinatorialDstarToDzeroPi(CollisionsWithEvSel const& collisions,
                                          aod::PiKaFromDzTags const& tagsPiKa,
                                          aod::KaPiFromDzTags const& tagsKaPi,
                                          aod::TrackAssoc const& trackIndices,
                                          TracksWithDca const& tracks)
  {
    for (const auto& collision : collisions) {
      if (!isCollisionSelected(collision)) {
        continue;
      }
      auto thisCollId = collision.globalIndex();
      auto groupedTrackIndices = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      // D*+ -> D0(->pi+K-)pi+
      auto groupedTagsPiKa = tagsPiKa.sliceBy(tagsPiKaPerCollision, thisCollId);
      runCombinatorialThirdTrack<aod::tagandprobe::TagChannels::DstarPlusToDzeroPi, false>(groupedTagsPiKa, groupedTrackIndices, tracks, tracks);
      // D*- -> D0bar(->K+pi-)pi-
      auto groupedTagsKaPi = tagsKaPi.sliceBy(tagsKaPiPerCollision, thisCollId);
      runCombinatorialThirdTrack<aod::tagandprobe::TagChannels::DstarMinusToDzeroBarPi, false>(groupedTagsKaPi, groupedTrackIndices, tracks, tracks);
    }
  }
  PROCESS_SWITCH(ProbeThirdTrack, processCombinatorialDstarToDzeroPi, "Process combinatorial of tagged pion-kaon (D0) vertices with additional track", true);

  void processCombinatorialDstarToDzeroPiMc(CollisionsWithEvSel const& collisions,
                                            FilteredPiKaFromDzMcTags const& tagsPiKa,
                                            FilteredKaPiFromDzMcTags const& tagsKaPi,
                                            aod::TrackAssoc const& trackIndices,
                                            TracksWithDcaMc const& tracks,
                                            aod::McParticles const& mcParticles)
  {
    for (const auto& collision : collisions) {
      if (!isCollisionSelected(collision)) {
        continue;
      }
      auto thisCollId = collision.globalIndex();
      auto groupedTrackIndices = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      // D*+ -> D0(->pi+K-)pi+
      auto groupedTagsPiKa = tagsPiKa.sliceBy(tagsPiKaMcPerCollision, thisCollId);
      runCombinatorialThirdTrack<aod::tagandprobe::TagChannels::DstarPlusToDzeroPi, true>(groupedTagsPiKa, groupedTrackIndices, tracks, mcParticles);
      // D*- -> D0bar(->K+pi-)pi-
      auto groupedTagsKaPi = tagsKaPi.sliceBy(tagsKaPiMcPerCollision, thisCollId);
      runCombinatorialThirdTrack<aod::tagandprobe::TagChannels::DstarMinusToDzeroBarPi, true>(groupedTagsKaPi, groupedTrackIndices, tracks, mcParticles);
    }
  }
  PROCESS_SWITCH(ProbeThirdTrack, processCombinatorialDstarToDzeroPiMc, "Process combinatorial of tagged pion-kaon (D0) vertices with additional track using MC truth", false);

  void processDummy(aod::Collisions const&) {}
  PROCESS_SWITCH(ProbeThirdTrack, processDummy, "Dummy process function that does nothing", false);

  void processGeneratedDstarToDzeroPi(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles)
  {

    if (std::abs(mcCollision.posZ()) > 10.)
      return;
    std::array<int, 3> arrDstar = {kPiPlus, kKMinus, kPiPlus};
    int8_t* sign = nullptr;
    std::vector<int> listIndexDaughters;
    float ptDzero = -1, yDzero = -999, ptSoftPion = -1, etaSoftPion = -999, ptminTagDaughers = 9999., etamaxTagDaugthers = 0.;
    int indexProbe;
    for (auto const& mcPart : mcParticles) {
      //  LOGP(info, "particle id: {}", mcPart.pdgCode());
      if (RecoDecay::isMatchedMCGen<true, true, 3>(mcParticles, mcPart, constants::physics::Pdg::kDStar, arrDstar, true, sign, 2, &listIndexDaughters)) {
        // LOGP(info, "Selected particle id: {}", mcPart.pdgCode());
        ptDzero = -1;
        yDzero = -999;
        ptSoftPion = -1;
        etaSoftPion = -999;
        ptminTagDaughers = 9999.;
        etamaxTagDaugthers = 0.;
        indexProbe = -1;
        for (auto const& iDaughtIndex : mcPart.daughtersIds()) {
          //  Printf("mcpart.daugthersIds, index: %d",idaughtindex);
          auto mcPartDstarDaught = mcParticles.rawIteratorAt(iDaughtIndex - mcParticles.offset());
          if (std::abs(mcPartDstarDaught.pdgCode()) == constants::physics::Pdg::kD0) {
            ptDzero = mcPartDstarDaught.pt();
            yDzero = mcPartDstarDaught.y();
          } else if (std::abs(mcPartDstarDaught.pdgCode()) == kPiPlus) {
            ptSoftPion = mcPartDstarDaught.pt();
            etaSoftPion = mcPartDstarDaught.eta();
            indexProbe = iDaughtIndex;
          }
        }
        for (auto const& idx : listIndexDaughters) {
          // LOGP(info, "listIndexDaughters, index: {}", idx);
          if (idx == indexProbe) {
            continue;
          }
          auto mcPartDaught = mcParticles.rawIteratorAt(idx - mcParticles.offset());
          ptminTagDaughers = std::min(mcPartDaught.pt(), ptminTagDaughers);
          etamaxTagDaugthers = std::max(std::abs(mcPartDaught.eta()), etamaxTagDaugthers);
        }
        // registry.fill(HIST(Form("hGen%s",tagChannels[aod::tagandprobe::DstarPlusToDzeroPi].data())),
        if (mcPart.pdgCode() > 0) {
          histosGen[aod::tagandprobe::DstarPlusToDzeroPi]->Fill(mcPart.pt(), mcPart.y(), ptDzero, yDzero, ptSoftPion, etaSoftPion, ptminTagDaughers, etamaxTagDaugthers);
        } else {
          histosGen[aod::tagandprobe::DstarMinusToDzeroBarPi]->Fill(mcPart.pt(), mcPart.y(), ptDzero, yDzero, ptSoftPion, etaSoftPion, ptminTagDaughers, etamaxTagDaugthers);
        }
      }
    }
  }
  PROCESS_SWITCH(ProbeThirdTrack, processGeneratedDstarToDzeroPi, "Count generated particles", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<TagTwoProngDisplacedVertices>(cfgc));
  workflow.push_back(adaptAnalysisTask<ProbeThirdTrack>(cfgc));
  return workflow;
}
