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
///
/// \file cascadeflow.cxx
///
/// \brief Task to create derived data for cascade flow analyses
/// \author Chiara De Martin (chiara.de.martin@cern.ch)
/// \author Maximiliano Puccio (maximiliano.puccio@cern.ch)

#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/cascqaanalysis.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Tools/ML/MlResponse.h"

#include "CCDB/BasicCCDBManager.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

#include "Math/Vector3D.h"
#include "TRandom3.h"

#include <memory>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::analysis;
using namespace o2::framework::expressions;
using std::array;

using DauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;
using CollEventPlane = soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraFT0CQVs, aod::StraFT0CQVsEv, aod::StraTPCQVs, aod::StraStamps>::iterator;
using CollEventPlaneCentralFW = soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraFT0CQVs, aod::StraFT0AQVs, aod::StraFT0MQVs, aod::StraFV0AQVs, aod::StraTPCQVs, aod::StraStamps>::iterator;
using CollEventPlaneCentralFWOnlyFT0C = soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraFT0CQVs, aod::StraTPCQVs, aod::StraStamps>::iterator;
using CollEventAndSpecPlane = soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraFT0CQVs, aod::StraFT0CQVsEv, aod::StraTPCQVs, aod::StraZDCSP, aod::StraStamps>::iterator;
using CollEventAndSpecPlaneCentralFW = soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraFT0CQVs, aod::StraTPCQVs, aod::StraZDCSP, aod::StraStamps>::iterator;
using MCCollisionsStra = soa::Join<aod::StraMCCollisions, aod::StraMCCollMults>;
using V0Candidates = soa::Join<aod::V0CollRefs, aod::V0Cores, aod::V0Extras>;
using V0MCCandidates = soa::Join<aod::V0CollRefs, aod::V0Cores, aod::V0Extras, aod::V0CoreMCLabels>;
using CascCandidates = soa::Join<aod::CascCollRefs, aod::CascCores, aod::CascExtras, aod::CascBBs>;
using CascMCCandidates = soa::Join<aod::CascCollRefs, aod::CascCores, aod::CascExtras, aod::CascBBs, aod::CascCoreMCLabels>;

const int nParticles = 2; // Xi, Omega
const int nCharges = 2;   // Lambda, AntiLambda
const int nParameters = 4;

namespace cascadev2
{
enum species { Xi = 0,
               Omega = 1 };
constexpr double massSigmaParameters[nParameters][nParticles]{
  {4.9736e-3, 0.006815},
  {-2.39594, -2.257},
  {1.8064e-3, 0.00138},
  {1.03468e-1, 0.1898}};
static const std::vector<std::string> massSigmaParameterNames{"p0", "p1", "p2", "p3"};
static const std::vector<std::string> speciesNames{"Xi", "Omega"};
const double AlphaXi[2] = {-0.390, 0.371};     // decay parameter of XiMinus and XiPlus
const double AlphaOmega[2] = {0.0154, -0.018}; // decay parameter of OmegaMinus and OmegaPlus
const double AlphaLambda[2] = {0.747, -0.757}; // decay parameter of Lambda and AntiLambda

std::shared_ptr<TH2> hMassBeforeSelVsPt[nParticles];
std::shared_ptr<TH2> hMassAfterSelVsPt[nParticles];
std::shared_ptr<TH1> hSignalScoreBeforeSel[nParticles];
std::shared_ptr<TH1> hBkgScoreBeforeSel[nParticles];
std::shared_ptr<TH1> hSignalScoreAfterSel[nParticles];
std::shared_ptr<TH1> hBkgScoreAfterSel[nParticles];
std::shared_ptr<THn> hSparseV2C[nParticles];
} // namespace cascadev2

namespace lambdav2
{
enum species { Lambda = 0,
               AntiLambda = 1 };
static const std::vector<std::string> speciesNames{"Lambda", "AntiLambda"};
const double AlphaLambda[2] = {0.747, -0.757}; // decay parameter of Lambda and AntiLambda

std::shared_ptr<TH2> hMassBeforeSelVsPt[nCharges];
std::shared_ptr<TH2> hMassAfterSelVsPt[nCharges];
std::shared_ptr<TH2> hMassAfterSelVsPtTrue[nCharges];
} // namespace lambdav2

namespace cascade_flow_cuts_ml
{
// direction of the cut
enum CutDirection {
  CutGreater = 0, // require score < cut value
  CutSmaller,     // require score > cut value
  CutNot          // do not cut on score
};

static constexpr int nBinsPt = 8;
static constexpr int nCutScores = 2;
// default values for the pT bin edges, offset by 1 from the bin numbers in cuts array
constexpr double binsPt[nBinsPt + 1] = {
  0.6,
  1.,
  2.,
  3.,
  4.,
  5.,
  6.,
  8.,
  10.};
auto vecBinsPt = std::vector<double>{binsPt, binsPt + nBinsPt + 1};

// default values for the ML model paths, one model per pT bin
static const std::vector<std::string> modelPaths = {""};

// default values for the cut directions
constexpr int cutDir[nCutScores] = {CutSmaller, CutNot}; // CutSmaller selects values > fixed value to signal BDT score, CutNot does not apply any selection to the background BDT score
auto vecCutDir = std::vector<int>{cutDir, cutDir + nCutScores};

// default values for the cuts
constexpr double cuts[nBinsPt][nCutScores] = { // background, signal
  {0., 0.9},
  {0., 0.9},
  {0., 0.9},
  {0., 0.9},
  {0., 0.9},
  {0., 0.9},
  {0., 0.9},
  {0., 0.9}};

// row labels
static const std::vector<std::string> labelsPt = {
  "pT bin 0",
  "pT bin 1",
  "pT bin 2",
  "pT bin 3",
  "pT bin 4",
  "pT bin 5",
  "pT bin 6",
  "pT bin 7"};

// column labels
static const std::vector<std::string> labelsCutScore = {"Background score", "Signal score"};
} // namespace cascade_flow_cuts_ml

struct cascadeFlow {

  Configurable<bool> isQVecT0C{"isQVecT0C", 1, ""};
  Configurable<bool> isQVecT0A{"isQVecT0A", 0, ""};
  Configurable<bool> isQVecT0M{"isQVecT0M", 0, ""};
  Configurable<bool> isQVecV0A{"isQVecV0A", 0, ""};
  Configurable<bool> isCollisionCentrality{"isCollisionCentrality", 0, ""}; // 0: FT0C, 1: FT0M (implemented only for Lambda analysis in OO)

  // Output filling criteria
  struct : ConfigurableGroup {
    Configurable<bool> isFillTree{"isFillTree", 1, ""};
    Configurable<bool> isFillTHNXi{"isFillTHNXi", 1, ""};
    Configurable<bool> isFillTHNXi_PzVsPsi{"isFillTHNXi_PzVsPsi", 1, ""};
    Configurable<bool> isFillTHNOmega{"isFillTHNOmega", 1, ""};
    Configurable<bool> isFillTHNOmega_PzVsPsi{"isFillTHNOmega_PzVsPsi", 1, ""};
    Configurable<bool> isFillTHNLambda{"isFillTHNLambda", 1, ""};
    Configurable<bool> isFillTHNLambda_PzVsPsi{"isFillTHNLambda_PzVsPsi", 1, ""};
    Configurable<bool> isFillTHN_V2{"isFillTHN_V2", 1, ""};
    Configurable<bool> isFillTHN_Pz{"isFillTHN_Pz", 1, ""};
    Configurable<bool> isFillTHN_PzFromLambda{"isFillTHN_PzFromLambda", 1, ""};
    Configurable<bool> isFillTHN_Acc{"isFillTHN_Acc", 1, ""};
    Configurable<bool> isFillTHN_AccFromLambdaVsCasc{"isFillTHN_AccFromLambdaVsCasc", 1, ""};
    Configurable<bool> isFillTHN_AccFromLambdaVsLambda{"isFillTHN_AccFromLambdaVsLambda", 1, ""};
  } fillingConfigs;

  // axes
  ConfigurableAxis axisQVs{"axisQVs", {500, -10.f, 10.f}, "axisQVs"};
  ConfigurableAxis axisQVsNorm{"axisQVsNorm", {200, -1.f, 1.f}, "axisQVsNorm"};

  // Configurable for shift correction
  struct : ConfigurableGroup {
    Configurable<bool> cfgShiftCorr{"cfgShiftCorr", 0, ""};
    Configurable<std::string> cfgShiftPathFT0C{"cfgShiftPathFT0C", "Users/c/chdemart/OOpass2Shift/ShiftFT0C", "Path for Shift"};
    Configurable<std::string> cfgShiftPathFV0A{"cfgShiftPathFV0A", "Users/c/chdemart/OOpass2Shift/ShiftFV0A", "Path for Shift"};
    Configurable<std::string> cfgShiftPathFT0A{"cfgShiftPathFT0A", "Users/c/chdemart/OOpass2Shift/ShiftFT0A", "Path for Shift"};
    Configurable<std::string> cfgShiftPathTPCL{"cfgShiftPathTPCL", "Users/c/chdemart/OOpass2Shift/ShiftTPCL", "Path for Shift"};
    Configurable<std::string> cfgShiftPathTPCR{"cfgShiftPathTPCR", "Users/c/chdemart/OOpass2Shift/ShiftTPCR", "Path for Shift"};
  } ShiftConfigs;
  //  Configurable<float> cfgHarmonic{"cfgHarmonic", 2, "Harmonic for event plane calculation"};

  // THN axes
  struct : ConfigurableGroup {
    ConfigurableAxis thnConfigAxisFT0C{"thnConfigAxisFT0C", {8, 0, 80}, "FT0C centrality (%)"};
    ConfigurableAxis thnConfigAxisEta{"thnConfigAxisEta", {8, -0.8, 0.8}, "pseudorapidity"};
    ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {VARIABLE_WIDTH, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2, 2.25, 2.5, 2.75, 3, 3.5, 4, 5, 6, 8, 10}, "#it{p}_{T} (GeV/#it{c})"};
    ConfigurableAxis thnConfigAxisPtLambda{"thnConfigAxisPtLambda", {VARIABLE_WIDTH, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2, 2.25, 2.5, 2.75, 3, 3.5, 4, 5, 6, 8, 10, 20}, "#it{p}_{T} (GeV/#it{c})"};
    ConfigurableAxis thnConfigAxisCharge{"thnConfigAxisCharge", {2, 0, 2}, ""};
    ConfigurableAxis thnConfigAxisPsiDiff{"thnConfigAxisPsiDiff", {100, 0, o2::constants::math::TwoPI}, ""};
    ConfigurableAxis thnConfigAxisMassXi{"thnConfigAxisMassXi", {45, 1.300, 1.345}, ""};
    ConfigurableAxis thnConfigAxisMassOmega{"thnConfigAxisMassOmega", {45, 1.655, 1.690}, ""};
    ConfigurableAxis thnConfigAxisMassLambda{"thnConfigAxisMassLambda", {60, 1.1, 1.13}, ""};
    ConfigurableAxis thnConfigAxisBDTScore{"thnConfigAxisBDTScore", {15, 0.4, 1}, ""};
    ConfigurableAxis thnConfigAxisV2{"thnConfigAxiV2", {100, -1., 1.}, ""};
    ConfigurableAxis thnConfigAxisPzs2Xi{"thnConfigAxiPzs2Xi", {200, -2.8, 2.8}, ""};
    ConfigurableAxis thnConfigAxisPzs2Omega{"thnConfigAxiPzs2Omega", {200, -70, 70}, ""};
    ConfigurableAxis thnConfigAxisPzs2Lambda{"thnConfigAxiPzs2Lambda", {200, -2, 2}, ""};
    ConfigurableAxis thnConfigAxisCos2Theta{"thnConfigAxiCos2Theta", {100, 0, 1}, ""};
    ConfigurableAxis thnConfigAxisCos2ThetaL{"thnConfigAxiCos2ThetaL", {100, 0, 1}, ""};
    ConfigurableAxis thnConfigAxisCosThetaXiAlpha{"thnConfigAxisCosThetaXiAlpha", {200, -2.8, 2.8}, ""};
    ConfigurableAxis thnConfigAxisCosThetaOmegaAlpha{"thnConfigAxisCosThetaOmegaAlpha", {200, -70, 70}, ""};
    ConfigurableAxis thnConfigAxisCosThetaProtonAlpha{"thnConfigAxisCosThetaProtonAlpha", {200, -2, 2}, ""};
  } thnAxisConfigs;

  // Event selection criteria
  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range (cm)"};
  Configurable<bool> sel8{"sel8", 1, "Apply sel8 event selection"};
  Configurable<bool> isNoSameBunchPileupCut{"isNoSameBunchPileupCut", 1, "Same found-by-T0 bunch crossing rejection"};
  Configurable<bool> isGoodZvtxFT0vsPVCut{"isGoodZvtxFT0vsPVCut", 1, "z of PV by tracks and z of PV from FT0 A-C time difference cut"};
  Configurable<bool> isGoodEventEP{"isGoodEventEP", 1, "Event is used to calibrate event plane"};
  Configurable<bool> isTrackOccupancySel{"isTrackOccupancySel", 0, "isTrackOccupancySel"};
  Configurable<int> MinOccupancy{"MinOccupancy", 0, "MinOccupancy"};
  Configurable<int> MaxOccupancy{"MaxOccupancy", 500, "MaxOccupancy"};
  Configurable<bool> isFT0OccupancySel{"isFT0OccupancySel", 0, "isFT0OccupancySel"};
  Configurable<int> MinOccupancyFT0{"MinOccupancyFT0", 0, "MinOccupancyFT0"};
  Configurable<int> MaxOccupancyFT0{"MaxOccupancyFT0", 5000, "MaxOccupancyFT0"};
  Configurable<bool> isNoCollInStandardTimeRange{"isNoCollInStandardTimeRange", 1, "To remove collisions in +-10 micros time range"};
  Configurable<bool> isNoCollInNarrowTimeRange{"isNoCollInNarrowTimeRange", 0, "To remove collisions in +-2 micros time range"};
  Configurable<bool> isNoCollInRofStandard{"isNoCollInRofStandard", 0, "To remove collisions in the same ITS ROF and with a multiplicity above a certain threshold"};
  Configurable<bool> isNoTVXinTRD{"isNoTVXinTRD", 0, "To remove collisions with trigger in TRD"};

  struct : ConfigurableGroup {
    Configurable<float> MinPt{"MinPt", 0.6, "Min pt of cascade"};
    Configurable<float> MaxPt{"MaxPt", 10, "Max pt of cascade"};
    Configurable<float> MinPtLambda{"MinPtLambda", 0.4, "Min pt of daughter lambda"};
    Configurable<float> MaxPtLambda{"MaxPtLambda", 10, "Max pt of daughter lambda"};
    Configurable<float> etaCasc{"etaCasc", 0.8, "etaCasc"};
    Configurable<float> etaLambdaMax{"etaLambdaMax", 0.8, "etaLambdaMax"};
    Configurable<float> MinLambdaMass{"MinLambdaMass", 1.1, ""};
    Configurable<float> MaxLambdaMass{"MaxLambdaMass", 1.13, ""};
    Configurable<float> MinXiMass{"MinXiMass", 1.300, ""};
    Configurable<float> MaxXiMass{"MaxXiMass", 1.345, ""};
    Configurable<float> MinOmegaMass{"MinOmegaMass", 1.655, ""};
    Configurable<float> MaxOmegaMass{"MaxOmegaMass", 1.690, ""};
  } CandidateConfigs;

  struct : ConfigurableGroup {
    Configurable<float> MinPtV0{"MinPtV0", 0.2, "Min pt of v0"};
    Configurable<float> MaxPtV0{"MaxPtV0", 10, "Max pt of v0"};
    Configurable<float> MinMassLambda{"MinMassLambda", 1.105, ""};
    Configurable<float> MaxMassLambda{"MaxMassLambda", 1.125, ""};
    Configurable<float> MinMassLambdaInTree{"MinMassLambdaInTree", 1.09, ""};
    Configurable<float> MaxMassLambdaInTree{"MaxMassLambdaInTree", 1.14, ""};
    Configurable<float> etaV0{"etaV0", 0.8, "etaV0"};
    Configurable<float> v0cospa{"v0cospa", 0.97, "min V0 CosPA"};
    Configurable<float> dcav0dau{"dcav0dau", 1.0, "max DCA V0 Daughters (cm)"};
    Configurable<float> dcanegtopv{"dcanegtopv", .05, "min DCA Neg To PV (cm)"};
    Configurable<float> dcapostopv{"dcapostopv", .05, "min DCA Pos To PV (cm)"};
    Configurable<float> v0radius{"v0radius", 1.2, "minimum V0 radius (cm)"};
    Configurable<float> v0radiusMax{"v0radiusMax", 1E5, "maximum V0 radius (cm)"};
    Configurable<float> rapidityLambda{"rapidityLambda", 0.5, "rapidityLambda"};
    Configurable<float> etaLambda{"etaLambda", 0.8, "etaLambda"};
    Configurable<float> dauTrackV0Eta{"dauTrackV0Eta", 0.8, "dauTrackV0Eta"};
  } V0Configs;

  Configurable<double> sideBandStart{"sideBandStart", 5, "Start of the sideband region in number of sigmas"};
  Configurable<double> sideBandEnd{"sideBandEnd", 7, "End of the sideband region in number of sigmas"};
  Configurable<double> downsample{"downsample", 1., "Downsample training output tree"};
  Configurable<bool> doNTPCSigmaCut{"doNTPCSigmaCut", 1, "doNtpcSigmaCut"};
  Configurable<float> nsigmatpcPr{"nsigmatpcPr", 5, "nsigmatpcPr"};
  Configurable<float> nsigmatpcPi{"nsigmatpcPi", 5, "nsigmatpcPi"};
  Configurable<float> mintpccrrows{"mintpccrrows", 70, "mintpccrrows"};

  Configurable<bool> isStoreTrueCascOnly{"isStoreTrueCascOnly", 1, ""};
  Configurable<float> etaCascMCGen{"etaCascMCGen", 0.8, "etaCascMCGen"};
  Configurable<float> yCascMCGen{"yCascMCGen", 0.5, "yCascMCGen"};

  struct : ConfigurableGroup {
    Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::vector<std::string>> modelPathsCCDBXi{"modelPathsCCDBXi", std::vector<std::string>{"Users/c/chdemart/CascadesFlow"}, "Paths of models on CCDB"};
    Configurable<std::vector<std::string>> modelPathsCCDBOmega{"modelPathsCCDBOmega", std::vector<std::string>{"Users/c/chdemart/CascadesFlow"}, "Paths of models on CCDB"};
    Configurable<std::vector<std::string>> onnxFileNamesXi{"onnxFileNamesXi", std::vector<std::string>{"model_onnx.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
    Configurable<std::vector<std::string>> onnxFileNamesOmega{"onnxFileNamesOmega", std::vector<std::string>{"model_onnx.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
    Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
    Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", true, "Flag to enable or disable the loading of models from CCDB"};
    Configurable<std::string> acceptancePathsCCDBXi{"acceptancePathsCCDBXi", "Users/c/chdemart/AcceptanceXi", "Paths of Xi acceptance on CCDB"};
    Configurable<std::string> acceptancePathsCCDBOmega{"acceptancePathsCCDBOmega", "Users/c/chdemart/AcceptanceOmega", "Paths of Omega acceptance on CCDB"};
    Configurable<std::string> acceptancePathsCCDBLambda{"acceptancePathsCCDBLambda", "Users/c/chdemart/AcceptanceLambda", "Paths of Lambda acceptance on CCDB"};
    Configurable<std::string> acceptancePathsCCDBPrimaryLambda{"acceptancePathsCCDBPrimaryLambda", "Users/c/chdemart/AcceptanceLambda", "Paths of PrimaryLambda acceptance on CCDB"};
    Configurable<std::string> acceptanceHistoNameCasc{"acceptanceHistoNameCasc", "histoCos2ThetaNoFit2D", "Histo name of acceptance on CCDB"};
    Configurable<std::string> acceptanceHistoNameLambda{"acceptanceHistoNameLambda", "histoCos2ThetaLambdaFromCNoFit2D", "Histo name of acceptance on CCDB"};
    Configurable<std::string> acceptanceHistoNamePrimaryLambda{"acceptanceHistoNamePrimaryLambda", "histoCos2ThetaLambdaFromCNoFit2D", "Histo name of acceptance on CCDB"};
    Configurable<std::string> resoPaths{"resoPath", "Users/c/chdemart/Resolution/", "Paths of resolution"};
    Configurable<std::string> resoHistoName{"resoHistoName", "hResoPerCentBinsV0A", "Histo name of resolution"};
    Configurable<std::string> centWeightPaths{"centWeightPath", "Users/c/chdemart/CentralityWeight/", "Paths of centrality weight"};
    Configurable<std::string> centWeightHistoName{"centWeightHistoName", "hCentWeight", "Histo name of centrality weight"};
  } ccdbConfigs;

  // ML inference
  Configurable<bool> isApplyML{"isApplyML", 1, "Flag to apply ML selections"};
  Configurable<std::vector<double>> binsPtMl{"binsPtMl", std::vector<double>{cascade_flow_cuts_ml::vecBinsPt}, "pT bin limits for ML application"};
  Configurable<std::vector<int>> cutDirMl{"cutDirMl", std::vector<int>{cascade_flow_cuts_ml::vecCutDir}, "Whether to reject score values greater or smaller than the threshold"};
  Configurable<LabeledArray<double>> cutsMl{"cutsMl", {cascade_flow_cuts_ml::cuts[0], cascade_flow_cuts_ml::nBinsPt, cascade_flow_cuts_ml::nCutScores, cascade_flow_cuts_ml::labelsPt, cascade_flow_cuts_ml::labelsCutScore}, "ML selections per pT bin"};
  Configurable<int> nClassesMl{"nClassesMl", static_cast<int>(cascade_flow_cuts_ml::nCutScores), "Number of classes in ML model"};

  // acceptance crrection
  Configurable<bool> applyAcceptanceCorrection{"applyAcceptanceCorrection", false, "apply acceptance correction"};
  Configurable<bool> applyResoCorrection{"applyResoCorrection", false, "apply resolution correction"};
  Configurable<bool> applyCentWeightCorrection{"applyCentWeightCorrection", false, "apply centrality weight correction"};

  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // Add objects needed for ML inference
  o2::analysis::MlResponse<float> mlResponseXi;
  o2::analysis::MlResponse<float> mlResponseOmega;

  template <typename TCollision>
  bool AcceptEvent(TCollision const& collision, bool isFillHisto)
  {
    if (isFillHisto) {
      histos.fill(HIST("hNEvents"), 0.5);
      histos.fill(HIST("hEventNchCorrelationBefCuts"), collision.multNTracksPVeta1(), collision.multNTracksGlobal());
      histos.fill(HIST("hEventPVcontributorsVsCentralityBefCuts"), collision.centFT0C(), collision.multNTracksPVeta1());
      histos.fill(HIST("hEventGlobalTracksVsCentralityBefCuts"), collision.centFT0C(), collision.multNTracksGlobal());
    }

    // Event selection if required
    if (sel8 && !collision.sel8()) {
      return false;
    }
    if (isFillHisto)
      histos.fill(HIST("hNEvents"), 1.5);

    // Z vertex selection
    if (std::abs(collision.posZ()) > cutzvertex) {
      return false;
    }
    if (isFillHisto)
      histos.fill(HIST("hNEvents"), 2.5);

    // kNoSameBunchPileup selection
    if (isNoSameBunchPileupCut && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if (isFillHisto)
      histos.fill(HIST("hNEvents"), 3.5);

    // kIsGoodZvtxFT0vsPV selection
    if (isGoodZvtxFT0vsPVCut && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if (isFillHisto)
      histos.fill(HIST("hNEvents"), 4.5);

    // occupancy cut
    int occupancy = collision.trackOccupancyInTimeRange();
    if (isTrackOccupancySel && (occupancy < MinOccupancy || occupancy > MaxOccupancy)) {
      return false;
    }
    // occupancy cut based on FT0C
    int occupancyFT0 = collision.ft0cOccupancyInTimeRange();
    if (isFT0OccupancySel && (occupancyFT0 < MinOccupancyFT0 || occupancyFT0 > MaxOccupancyFT0)) {
      return false;
    }

    if (isFillHisto)
      histos.fill(HIST("hNEvents"), 5.5);

    // time-based event selection
    if (isNoCollInStandardTimeRange && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return false;
    }
    if (isNoCollInNarrowTimeRange && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict)) {
      return false;
    }

    if (isFillHisto)
      histos.fill(HIST("hNEvents"), 6.5);

    // In-ROF event selection
    if (isNoCollInRofStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      return false;
    }

    if (isFillHisto)
      histos.fill(HIST("hNEvents"), 7.5);

    // TVX in TRD
    //  if (isNoTVXinTRD && collision.alias_bit(kTVXinTRD)){
    //   return false;
    //  }

    if (isFillHisto)
      histos.fill(HIST("hNEvents"), 8.5);

    if (isFillHisto) {
      histos.fill(HIST("hEventNchCorrelation"), collision.multNTracksPVeta1(), collision.multNTracksGlobal());
      histos.fill(HIST("hEventPVcontributorsVsCentrality"), collision.centFT0C(), collision.multNTracksPVeta1());
      histos.fill(HIST("hEventGlobalTracksVsCentrality"), collision.centFT0C(), collision.multNTracksGlobal());
    }

    return true;
  }

  template <typename TCascade, typename TDaughter>
  bool IsCascAccepted(TCascade casc, TDaughter negExtra, TDaughter posExtra, TDaughter bachExtra, int& counter) // loose cuts on topological selections of cascades
  {
    // TPC cuts as those implemented for the training of the signal
    if (doNTPCSigmaCut) {
      if (casc.sign() < 0) {
        if (std::abs(posExtra.tpcNSigmaPr()) > nsigmatpcPr || std::abs(negExtra.tpcNSigmaPi()) > nsigmatpcPi)
          return false;
      } else if (casc.sign() > 0) {
        if (std::abs(posExtra.tpcNSigmaPi()) > nsigmatpcPi || std::abs(negExtra.tpcNSigmaPr()) > nsigmatpcPr)
          return false;
      }
    }
    counter++;

    if (posExtra.tpcCrossedRows() < mintpccrrows || negExtra.tpcCrossedRows() < mintpccrrows || bachExtra.tpcCrossedRows() < mintpccrrows)
      return false;

    counter++;
    return true;
  }

  template <typename TDaughter>
  bool isLambdaAccepted(TDaughter negExtra, TDaughter posExtra, int& counter) // loose cuts on topological selections of v0s
  {
    // TPC cuts as those implemented for the training of the signal
    if (doNTPCSigmaCut) {
      if (std::abs(posExtra.tpcNSigmaPr()) > nsigmatpcPr || std::abs(negExtra.tpcNSigmaPi()) > nsigmatpcPi)
        return false;
    }
    counter++;

    if (posExtra.tpcCrossedRows() < mintpccrrows || negExtra.tpcCrossedRows() < mintpccrrows)
      return false;

    counter++;

    // eta daughters)
    //     if (abs(posExtra.eta()) > V0Configs.dauTrackV0Eta || abs(negExtra.y()) > V0Configs.dauTrackV0Eta) return false;

    return true;
  }
  template <typename TDaughter>
  bool isAntiLambdaAccepted(TDaughter negExtra, TDaughter posExtra, int& counter) // loose cuts on topological selections of v0s
  {
    // TPC cuts as those implemented for the training of the signal
    if (doNTPCSigmaCut) {
      if (std::abs(negExtra.tpcNSigmaPr()) > nsigmatpcPr || std::abs(posExtra.tpcNSigmaPi()) > nsigmatpcPi)
        return false;
    }
    counter++;

    if (posExtra.tpcCrossedRows() < mintpccrrows || negExtra.tpcCrossedRows() < mintpccrrows)
      return false;

    counter++;
    return true;
  }

  template <typename TV0>
  bool isV0TopoAccepted(TV0 v0)
  {
    // topological selections
    if (v0.v0radius() < V0Configs.v0radius)
      return false;
    if (v0.v0radius() > V0Configs.v0radiusMax)
      return false;
    if (std::abs(v0.dcapostopv()) < V0Configs.dcapostopv)
      return false;
    if (std::abs(v0.dcanegtopv()) < V0Configs.dcanegtopv)
      return false;
    if (v0.v0cosPA() < V0Configs.v0cospa)
      return false;
    if (v0.dcaV0daughters() > V0Configs.dcav0dau)
      return false;
    // rapidity selection
    if (std::abs(v0.yLambda()) > V0Configs.rapidityLambda)
      return false;
    if (std::abs(v0.eta()) > V0Configs.etaLambda)
      return false;

    return true;
  }

  double GetPhiInRange(double phi)
  {
    while (phi < 0) {
      phi += o2::constants::math::PI;
    }
    while (phi > o2::constants::math::PI) {
      phi -= o2::constants::math::PI;
    }
    return phi;
  }

  int currentRunNumber = -999;
  int lastRunNumber = -999;
  TProfile3D* shiftprofile;
  TProfile3D* shiftprofileFT0C;
  TProfile3D* shiftprofileFV0A;
  TProfile3D* shiftprofileFT0A;
  TProfile3D* shiftprofileTPCL;
  TProfile3D* shiftprofileTPCR;
  std::string fullCCDBShiftCorrPath;
  std::string fullCCDBShiftCorrPathFT0C;
  std::string fullCCDBShiftCorrPathFV0A;
  std::string fullCCDBShiftCorrPathFT0A;
  std::string fullCCDBShiftCorrPathTPCL;
  std::string fullCCDBShiftCorrPathTPCR;

  template <typename TCollision>
  double ApplyShiftCorrection(TCollision coll, double psiT0C, TProfile3D* shiftprofile)
  {
    auto deltapsiFT0C = 0.0;
    int nmode = 2;

    for (int ishift = 1; ishift <= 10; ishift++) {
      auto coeffshiftxFT0C = shiftprofile->GetBinContent(shiftprofile->FindBin(coll.centFT0C(), 0.5, ishift - 0.5));
      auto coeffshiftyFT0C = shiftprofile->GetBinContent(shiftprofile->FindBin(coll.centFT0C(), 1.5, ishift - 0.5));

      deltapsiFT0C += ((1 / (1.0 * ishift)) * (-coeffshiftxFT0C * std::cos(ishift * static_cast<float>(nmode) * psiT0C) + coeffshiftyFT0C * TMath::Sin(ishift * static_cast<float>(nmode) * psiT0C)));
    }
    return psiT0C + deltapsiFT0C;
  }

  template <typename TCollision>
  double ComputeEPResolutionwShifts(TCollision coll, double psiT0C, double psiV0A, double psiT0A, double psiTPCA, double psiTPCC, TProfile3D* shiftprofileA, TProfile3D* shiftprofileB, TProfile3D* shiftprofileC, TProfile3D* shiftprofileD, TProfile3D* shiftprofileE)
  {
    int nmode = 2;
    auto deltapsiFT0C = 0.0;
    auto deltapsiFV0A = 0.0;
    auto deltapsiFT0A = 0.0;
    auto deltapsiTPCA = 0.0;
    auto deltapsiTPCC = 0.0;
    for (int ishift = 1; ishift <= 10; ishift++) {
      auto coeffshiftxFT0C = shiftprofileA->GetBinContent(shiftprofileA->FindBin(coll.centFT0C(), 0.5, ishift - 0.5));
      auto coeffshiftyFT0C = shiftprofileA->GetBinContent(shiftprofileA->FindBin(coll.centFT0C(), 1.5, ishift - 0.5));
      auto coeffshiftxTPCA = shiftprofileB->GetBinContent(shiftprofileB->FindBin(coll.centFT0C(), 0.5, ishift - 0.5));
      auto coeffshiftyTPCA = shiftprofileB->GetBinContent(shiftprofileB->FindBin(coll.centFT0C(), 1.5, ishift - 0.5));
      auto coeffshiftxTPCC = shiftprofileC->GetBinContent(shiftprofileC->FindBin(coll.centFT0C(), 0.5, ishift - 0.5));
      auto coeffshiftyTPCC = shiftprofileC->GetBinContent(shiftprofileC->FindBin(coll.centFT0C(), 1.5, ishift - 0.5));
      auto coeffshiftxFV0A = shiftprofileD->GetBinContent(shiftprofileD->FindBin(coll.centFT0C(), 0.5, ishift - 0.5));
      auto coeffshiftyFV0A = shiftprofileD->GetBinContent(shiftprofileD->FindBin(coll.centFT0C(), 1.5, ishift - 0.5));
      auto coeffshiftxFT0A = shiftprofileE->GetBinContent(shiftprofileE->FindBin(coll.centFT0C(), 0.5, ishift - 0.5));
      auto coeffshiftyFT0A = shiftprofileE->GetBinContent(shiftprofileE->FindBin(coll.centFT0C(), 1.5, ishift - 0.5));
      deltapsiFT0C += ((1 / (1.0 * ishift)) * (-coeffshiftxFT0C * std::cos(ishift * static_cast<float>(nmode) * psiT0C) + coeffshiftyFT0C * TMath::Sin(ishift * static_cast<float>(nmode) * psiT0C)));
      deltapsiFV0A += ((1 / (1.0 * ishift)) * (-coeffshiftxFV0A * std::cos(ishift * static_cast<float>(nmode) * psiV0A) + coeffshiftyFV0A * TMath::Sin(ishift * static_cast<float>(nmode) * psiV0A)));
      deltapsiFT0A += ((1 / (1.0 * ishift)) * (-coeffshiftxFT0A * std::cos(ishift * static_cast<float>(nmode) * psiT0A) + coeffshiftyFT0A * TMath::Sin(ishift * static_cast<float>(nmode) * psiT0A)));
      deltapsiTPCA += ((1 / (1.0 * ishift)) * (-coeffshiftxTPCA * std::cos(ishift * static_cast<float>(nmode) * psiTPCA) + coeffshiftyTPCA * TMath::Sin(ishift * static_cast<float>(nmode) * psiTPCA)));
      deltapsiTPCC += ((1 / (1.0 * ishift)) * (-coeffshiftxTPCC * std::cos(ishift * static_cast<float>(nmode) * psiTPCC) + coeffshiftyTPCC * TMath::Sin(ishift * static_cast<float>(nmode) * psiTPCC)));
    }
    histos.fill(HIST("Psi_EP_FT0C_shifted"), coll.centFT0C(), psiT0C + deltapsiFT0C);
    histos.fill(HIST("Psi_EP_FV0A_shifted"), coll.centFT0C(), psiV0A + deltapsiFV0A);
    histos.fill(HIST("Psi_EP_FT0A_shifted"), coll.centFT0C(), psiT0A + deltapsiFT0A);
    histos.fill(HIST("Psi_EP_TPCA_shifted"), coll.centFT0C(), psiTPCA + deltapsiTPCA);
    histos.fill(HIST("Psi_EP_TPCC_shifted"), coll.centFT0C(), psiTPCC + deltapsiTPCC);
    resolution.fill(HIST("QVectorsT0CTPCA_Shifted"), std::cos(static_cast<float>(nmode) * (psiT0C + deltapsiFT0C - psiTPCA - deltapsiTPCA)), coll.centFT0C());
    resolution.fill(HIST("QVectorsT0CTPCC_Shifted"), std::cos(static_cast<float>(nmode) * (psiT0C + deltapsiFT0C - psiTPCC - deltapsiTPCC)), coll.centFT0C());
    resolution.fill(HIST("QVectorsT0CV0A_Shifted"), std::cos(static_cast<float>(nmode) * (psiT0C + deltapsiFT0C - psiV0A - deltapsiFV0A)), coll.centFT0C());
    resolution.fill(HIST("QVectorsT0CT0A_Shifted"), std::cos(static_cast<float>(nmode) * (psiT0C + deltapsiFT0C - psiT0A - deltapsiFT0A)), coll.centFT0C());
    resolution.fill(HIST("QVectorsV0ATPCC_Shifted"), std::cos(static_cast<float>(nmode) * (psiV0A + deltapsiFV0A - psiTPCC - deltapsiTPCC)), coll.centFT0C());
    resolution.fill(HIST("QVectorsV0ATPCA_Shifted"), std::cos(static_cast<float>(nmode) * (psiV0A + deltapsiFV0A - psiTPCA - deltapsiTPCA)), coll.centFT0C());
    resolution.fill(HIST("QVectorsT0ATPCC_Shifted"), std::cos(static_cast<float>(nmode) * (psiT0A + deltapsiFT0A - psiTPCC - deltapsiTPCC)), coll.centFT0C());
    resolution.fill(HIST("QVectorsT0ATPCA_Shifted"), std::cos(static_cast<float>(nmode) * (psiT0A + deltapsiFT0A - psiTPCA - deltapsiTPCA)), coll.centFT0C());
    resolution.fill(HIST("QVectorsTPCAC_Shifted"), std::cos(static_cast<float>(nmode) * (psiTPCA + deltapsiTPCA - psiTPCC - deltapsiTPCC)), coll.centFT0C());
    return true;
  }

  // objects to use for acceptance correction
  TH2F* hAcceptanceXi;
  TH2F* hAcceptanceOmega;
  TH2F* hAcceptanceLambda;
  TH2F* hAcceptancePrimaryLambda;

  // objects to use for resolution correction
  TH1F* hReso;

  // objects to use for centrality weight
  TH1F* hCentWeight;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry histosMCGen{"histosMCGen", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry resolution{"resolution", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};

  // Tables to produce
  Produces<aod::CascTraining> trainingSample;
  Produces<aod::CascAnalysis> analysisSample;
  Produces<aod::LambdaAnalysis> analysisLambdaSample;
  Configurable<LabeledArray<double>> parSigmaMass{
    "parSigmaMass",
    {cascadev2::massSigmaParameters[0], nParameters, nParticles,
     cascadev2::massSigmaParameterNames, cascadev2::speciesNames},
    "Mass resolution parameters: [0]*exp([1]*x)+[2]*exp([3]*x)"};

  float getNsigmaMass(const cascadev2::species s, const float pt, const float nsigma = 6)
  {
    const auto sigma = parSigmaMass->get(0u, s) * std::exp(parSigmaMass->get(1, s) * pt) + parSigmaMass->get(2, s) * std::exp(parSigmaMass->get(3, s) * pt);
    return nsigma * sigma;
  }

  template <class collision_t, class cascade_t>
  void fillTrainingTable(collision_t coll, cascade_t casc, int pdgCode)
  {
    trainingSample(coll.centFT0C(),
                   casc.sign(),
                   casc.pt(),
                   casc.eta(),
                   casc.mXi(),
                   casc.mOmega(),
                   casc.mLambda(),
                   casc.cascradius(),
                   casc.v0radius(),
                   casc.casccosPA(coll.posX(), coll.posY(), coll.posZ()),
                   casc.v0cosPA(coll.posX(), coll.posY(), coll.posZ()),
                   casc.dcapostopv(),
                   casc.dcanegtopv(),
                   casc.dcabachtopv(),
                   casc.dcacascdaughters(),
                   casc.dcaV0daughters(),
                   casc.dcav0topv(coll.posX(), coll.posY(), coll.posZ()),
                   casc.bachBaryonCosPA(),
                   casc.bachBaryonDCAxyToPV(),
                   pdgCode);
  }

  template <class collision_t, class cascade_t>
  void fillAnalysedTable(collision_t coll, bool hasEventPlane, bool hasSpectatorPlane, cascade_t casc, float v2CSP, float v2CEP, float v1SP_ZDCA, float v1SP_ZDCC, float PsiT0C, float BDTresponseXi, float BDTresponseOmega, int pdgCode)
  {
    double masses[nParticles]{o2::constants::physics::MassXiMinus, o2::constants::physics::MassOmegaMinus};
    ROOT::Math::PxPyPzMVector cascadeVector[nParticles], lambdaVector, protonVector;
    float cosThetaStarLambda[nParticles], cosThetaStarProton;
    lambdaVector.SetCoordinates(casc.pxlambda(), casc.pylambda(), casc.pzlambda(), o2::constants::physics::MassLambda);
    ROOT::Math::Boost lambdaBoost{lambdaVector.BoostToCM()};
    if (casc.sign() > 0) {
      protonVector.SetCoordinates(casc.pxneg(), casc.pyneg(), casc.pzneg(), o2::constants::physics::MassProton);
    } else {
      protonVector.SetCoordinates(casc.pxpos(), casc.pypos(), casc.pzpos(), o2::constants::physics::MassProton);
    }
    auto boostedProton{lambdaBoost(protonVector)};
    cosThetaStarProton = boostedProton.Pz() / boostedProton.P();
    for (int i{0}; i < nParticles; ++i) {
      cascadeVector[i].SetCoordinates(casc.px(), casc.py(), casc.pz(), masses[i]);
      ROOT::Math::Boost cascadeBoost{cascadeVector[i].BoostToCM()};
      auto boostedLambda{cascadeBoost(lambdaVector)};
      cosThetaStarLambda[i] = boostedLambda.Pz() / boostedLambda.P();
    }

    // time-based event selection
    bool isNoCollInTimeRangeStd = 0;
    if (coll.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard))
      isNoCollInTimeRangeStd = 1;

    // IN-ROF pile-up rejection
    bool isNoCollInRofStd = 0;
    if (coll.selection_bit(o2::aod::evsel::kNoCollInRofStandard))
      isNoCollInRofStd = 1;

    // TVX in TRD
    // bool isTVXinTRD = 0;
    //    if (coll.alias_bit(kTVXinTRD)) isTVXinTRD = 1;

    analysisSample(coll.centFT0C(),
                   isNoCollInTimeRangeStd,
                   isNoCollInRofStd,
                   hasEventPlane,
                   hasSpectatorPlane,
                   casc.sign(),
                   casc.pt(),
                   casc.eta(),
                   casc.phi(),
                   casc.mLambda(),
                   casc.mXi(),
                   casc.mOmega(),
                   v2CSP,
                   v2CEP,
                   v1SP_ZDCA,
                   v1SP_ZDCC,
                   PsiT0C,
                   BDTresponseXi,
                   BDTresponseOmega,
                   cosThetaStarLambda[0],
                   cosThetaStarLambda[1],
                   cosThetaStarProton,
                   pdgCode);
  }

  template <class collision_t, class v0_t>
  void fillAnalysedLambdaTable(collision_t coll, bool hasEventPlane, bool hasSpectatorPlane, int chargeIndex, v0_t v0, float v2CEP, float psiT0C, double pzs2Lambda, double cos2ThetaLambda, double cosThetaLambda)
  {
    double invMassLambda = 0;
    if (chargeIndex == 0)
      invMassLambda = v0.mLambda();
    else if (chargeIndex == 1)
      invMassLambda = v0.mAntiLambda();
    else
      invMassLambda = v0.mLambda();
    analysisLambdaSample(coll.centFT0C(),
                         hasEventPlane,
                         hasSpectatorPlane,
                         chargeIndex,
                         v0.pt(),
                         v0.phi(),
                         v0.eta(),
                         invMassLambda,
                         v0.v0radius(),
                         v0.dcapostopv(),
                         v0.dcanegtopv(),
                         v0.v0cosPA(),
                         v0.dcaV0daughters(),
                         v2CEP,
                         psiT0C,
                         pzs2Lambda,
                         cos2ThetaLambda,
                         cosThetaLambda);
  }

  void initAcceptanceFromCCDB()
  {
    LOG(info) << "Loading acceptance from CCDB ";
    TList* listAcceptanceXi = ccdb->get<TList>(ccdbConfigs.acceptancePathsCCDBXi);
    if (!listAcceptanceXi)
      LOG(fatal) << "Problem getting TList object with acceptance for Xi!";
    TList* listAcceptanceOmega = ccdb->get<TList>(ccdbConfigs.acceptancePathsCCDBOmega);
    if (!listAcceptanceOmega)
      LOG(fatal) << "Problem getting TList object with acceptance for Omega!";
    TList* listAcceptanceLambda = ccdb->get<TList>(ccdbConfigs.acceptancePathsCCDBLambda);
    if (!listAcceptanceLambda)
      LOG(fatal) << "Problem getting TList object with acceptance for Lambda!";
    TList* listAcceptancePrimaryLambda = ccdb->get<TList>(ccdbConfigs.acceptancePathsCCDBPrimaryLambda);
    if (!listAcceptancePrimaryLambda)
      LOG(fatal) << "Problem getting TList object with acceptance for Primary Lambda!";

    hAcceptanceXi = static_cast<TH2F*>(listAcceptanceXi->FindObject(Form("%s", ccdbConfigs.acceptanceHistoNameCasc->data())));
    if (!hAcceptanceXi) {
      LOG(fatal) << "The histogram for Xi is not there!";
    }
    hAcceptanceXi->SetName("hAcceptanceXi");
    hAcceptanceOmega = static_cast<TH2F*>(listAcceptanceOmega->FindObject(Form("%s", ccdbConfigs.acceptanceHistoNameCasc->data())));
    if (!hAcceptanceOmega) {
      LOG(fatal) << "The histogram for omega is not there!";
    }
    hAcceptanceOmega->SetName("hAcceptanceOmega");
    hAcceptanceLambda = static_cast<TH2F*>(listAcceptanceLambda->FindObject(Form("%s", ccdbConfigs.acceptanceHistoNameLambda->data())));
    if (!hAcceptanceLambda) {
      LOG(fatal) << "The histogram for Lambda is not there!";
    }
    hAcceptanceLambda->SetName("hAcceptanceLambda");
    hAcceptancePrimaryLambda = static_cast<TH2F*>(listAcceptancePrimaryLambda->FindObject(Form("%s", ccdbConfigs.acceptanceHistoNamePrimaryLambda->data())));
    if (!hAcceptancePrimaryLambda) {
      LOG(fatal) << "The histogram for Primary Lambda is not there!";
    }
    hAcceptancePrimaryLambda->SetName("hAcceptancePrimaryLambda");
    LOG(info) << "Acceptance now loaded";
  }
  void initResoFromCCDB()
  {
    LOG(info) << "Loading resolution from CCDB ";
    TList* listReso = ccdb->get<TList>(ccdbConfigs.resoPaths);
    if (!listReso)
      LOG(fatal) << "Problem getting TList object with resolution!";

    hReso = static_cast<TH1F*>(listReso->FindObject(Form("%s", ccdbConfigs.resoHistoName->data())));
    if (!hReso) {
      LOG(fatal) << "The histogram for resolution is not there";
    }
    hReso->SetName("hReso");
    LOG(info) << "Resolution now loaded";
  }

  void initCentWeightFromCCDB()
  {
    LOG(info) << "Loading resolution from CCDB ";
    TList* listCentWeight = ccdb->get<TList>(ccdbConfigs.centWeightPaths);
    if (!listCentWeight)
      LOG(fatal) << "Problem getting TList object with resolution!";

    hCentWeight = static_cast<TH1F*>(listCentWeight->FindObject(Form("%s", ccdbConfigs.centWeightHistoName->data())));
    if (!hCentWeight) {
      LOG(fatal) << "The histogram for resolution is not there";
    }
    hCentWeight->SetName("hCentWeight");
    LOG(info) << "Centrality weight now loaded";
  }

  void init(InitContext const&)
  {

    float minMass[2]{1.28, 1.6};
    float maxMass[2]{1.36, 1.73};
    float minMassLambda[2]{1.09, 1.09};
    float maxMassLambda[2]{1.14, 1.14};
    const AxisSpec shiftAxis = {10, 0, 10, "shift"};
    const AxisSpec basisAxis = {2, 0, 2, "basis"};
    const AxisSpec massCascAxis[2]{{static_cast<int>((maxMass[0] - minMass[0]) / 0.001f), minMass[0], maxMass[0], "#Xi candidate mass (GeV/c^{2})"},
                                   {static_cast<int>((maxMass[1] - minMass[1]) / 0.001f), minMass[1], maxMass[1], "#Omega candidate mass (GeV/c^{2})"}};
    const AxisSpec massLambdaAxis[2]{{static_cast<int>((maxMassLambda[0] - minMassLambda[0]) / 0.001f), minMassLambda[0], maxMassLambda[0], "#Lambda candidate mass (GeV/c^{2})"},
                                     {static_cast<int>((maxMassLambda[1] - minMassLambda[1]) / 0.001f), minMassLambda[1], maxMassLambda[1], "#bar{#Lambda} candidate mass (GeV/c^{2})"}};
    const AxisSpec ptAxisCasc{static_cast<int>((CandidateConfigs.MaxPt - CandidateConfigs.MinPt) / 0.2), CandidateConfigs.MinPt, CandidateConfigs.MaxPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec ptAxisLambda{static_cast<int>((V0Configs.MaxPtV0 - V0Configs.MinPtV0) / 0.2), V0Configs.MinPtV0, V0Configs.MaxPtV0, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec v2Axis{200, -1., 1., "#it{v}_{2}"};
    const AxisSpec CentAxis{20, 0., 100., "FT0C centrality percentile"};
    const AxisSpec CentAxisPerCent{100, 0., 100., "Percent FT0C centrality percentile"};
    TString hNEventsLabels[10] = {"All", "sel8", "z vrtx", "kNoSameBunchPileup", "kIsGoodZvtxFT0vsPV", "trackOccupancyInTimeRange", "kNoCollInTimeRange", "kNoCollInROF", "kTVXinTRD", "kIsGoodEventEP"};
    TString hNEventsLabelsMC[6] = {"All", "z vtx", ">=1RecoColl", "1Reco", "2Reco", "EvSelected"};
    TString hNCascLabelsMC[8] = {"All Xi", "all Omega", "Xi: has MC coll", "Om: has MC coll", "Xi: isPrimary", "Om: is Primary", "Xi: |eta|<0.8", "Om: |eta| < 0.8"};

    resolution.add("QVectorsT0CTPCA", "QVectorsT0CTPCA", HistType::kTH2F, {axisQVs, CentAxisPerCent});
    resolution.add("QVectorsT0CTPCC", "QVectorsT0CTPCC", HistType::kTH2F, {axisQVs, CentAxisPerCent});
    resolution.add("QVectorsTPCAC", "QVectorsTPCAC", HistType::kTH2F, {axisQVs, CentAxisPerCent});
    resolution.add("QVectorsT0CV0A", "QVectorsT0CV0A", HistType::kTH2F, {axisQVs, CentAxisPerCent});
    resolution.add("QVectorsV0ATPCC", "QVectorsV0ATPCC", HistType::kTH2F, {axisQVs, CentAxisPerCent});
    resolution.add("QVectorsV0ATPCA", "QVectorsV0ATPCA", HistType::kTH2F, {axisQVs, CentAxisPerCent});
    resolution.add("QVectorsT0CT0A", "QVectorsT0CT0A", HistType::kTH2F, {axisQVs, CentAxisPerCent});
    resolution.add("QVectorsT0ATPCC", "QVectorsT0ATPCC", HistType::kTH2F, {axisQVs, CentAxisPerCent});
    resolution.add("QVectorsT0ATPCA", "QVectorsT0ATPCA", HistType::kTH2F, {axisQVs, CentAxisPerCent});
    resolution.add("EP_T0CTPCA", "EP_T0CTPCA", HistType::kTH2F, {axisQVs, CentAxisPerCent});
    resolution.add("EP_T0CTPCC", "EP_T0CTPCC", HistType::kTH2F, {axisQVs, CentAxisPerCent});
    resolution.add("EP_TPCAC", "EP_TPCAC", HistType::kTH2F, {axisQVs, CentAxisPerCent});
    resolution.add("EP_T0CV0A", "EP_T0CV0A", HistType::kTH2F, {axisQVs, CentAxisPerCent});
    resolution.add("EP_V0ATPCC", "EP_V0ATPCC", HistType::kTH2F, {axisQVs, CentAxisPerCent});
    resolution.add("EP_V0ATPCA", "EP_V0ATPCA", HistType::kTH2F, {axisQVs, CentAxisPerCent});
    resolution.add("EP_T0CT0A", "EP_T0CT0A", HistType::kTH2F, {axisQVs, CentAxisPerCent});
    resolution.add("EP_T0ATPCC", "EP_T0ATPCC", HistType::kTH2F, {axisQVs, CentAxisPerCent});
    resolution.add("EP_T0ATPCA", "EP_T0ATPCA", HistType::kTH2F, {axisQVs, CentAxisPerCent});
    resolution.add("QVectorsNormT0CTPCA", "QVectorsNormT0CTPCA", HistType::kTH2F, {axisQVsNorm, CentAxisPerCent});
    resolution.add("QVectorsNormT0CTPCC", "QVectorsNormT0CTPCC", HistType::kTH2F, {axisQVsNorm, CentAxisPerCent});
    resolution.add("QVectorsNormTPCAC", "QVectorsNormTPCCB", HistType::kTH2F, {axisQVsNorm, CentAxisPerCent});
    resolution.add("QVectorsNormT0CV0A", "QVectorsNormT0CV0A", HistType::kTH2F, {axisQVs, CentAxisPerCent});
    resolution.add("QVectorsNormV0ATPCC", "QVectorsNormV0ATPCC", HistType::kTH2F, {axisQVs, CentAxisPerCent});
    resolution.add("QVectorsNormV0ATPCA", "QVectorsNormV0ATPCA", HistType::kTH2F, {axisQVs, CentAxisPerCent});
    resolution.add("QVectorsNormT0CT0A", "QVectorsNormT0CT0A", HistType::kTH2F, {axisQVs, CentAxisPerCent});
    resolution.add("QVectorsNormT0ATPCC", "QVectorsNormT0ATPCC", HistType::kTH2F, {axisQVs, CentAxisPerCent});
    resolution.add("QVectorsNormT0ATPCA", "QVectorsNormT0ATPCA", HistType::kTH2F, {axisQVs, CentAxisPerCent});
    resolution.add("QVectorsSpecPlane", "QVectorsSpecPlane", HistType::kTH2F, {axisQVsNorm, CentAxisPerCent});
    resolution.add("QVectorsT0CTPCA_Shifted", "QVectorsT0CTPCA_Shifted", HistType::kTH2F, {axisQVs, CentAxisPerCent});
    resolution.add("QVectorsT0CTPCC_Shifted", "QVectorsT0CTPCC_Shifted", HistType::kTH2F, {axisQVs, CentAxisPerCent});
    resolution.add("QVectorsTPCAC_Shifted", "QVectorsTPCAC_Shifted", HistType::kTH2F, {axisQVs, CentAxisPerCent});
    resolution.add("QVectorsT0CV0A_Shifted", "QVectorsT0CV0A_Shifted", HistType::kTH2F, {axisQVs, CentAxisPerCent});
    resolution.add("QVectorsV0ATPCC_Shifted", "QVectorsV0ATPCC_Shifted", HistType::kTH2F, {axisQVs, CentAxisPerCent});
    resolution.add("QVectorsV0ATPCA_Shifted", "QVectorsV0ATPCA_Shifted", HistType::kTH2F, {axisQVs, CentAxisPerCent});
    resolution.add("QVectorsT0CT0A_Shifted", "QVectorsT0CT0A_Shifted", HistType::kTH2F, {axisQVs, CentAxisPerCent});
    resolution.add("QVectorsT0ATPCC_Shifted", "QVectorsT0ATPCC_Shifted", HistType::kTH2F, {axisQVs, CentAxisPerCent});
    resolution.add("QVectorsT0ATPCA_Shifted", "QVectorsT0ATPCA_Shifted", HistType::kTH2F, {axisQVs, CentAxisPerCent});

    histos.add("ShiftFT0C", "ShiftFT0C", kTProfile3D, {CentAxis, basisAxis, shiftAxis});
    histos.add("ShiftFV0A", "ShiftFV0A", kTProfile3D, {CentAxis, basisAxis, shiftAxis});
    histos.add("ShiftFT0A", "ShiftFT0A", kTProfile3D, {CentAxis, basisAxis, shiftAxis});
    histos.add("ShiftTPCL", "ShiftTPCL", kTProfile3D, {CentAxis, basisAxis, shiftAxis});
    histos.add("ShiftTPCR", "ShiftTPCR", kTProfile3D, {CentAxis, basisAxis, shiftAxis});

    histos.add("hNEvents", "hNEvents", {HistType::kTH1D, {{10, 0.f, 10.f}}});
    for (Int_t n = 1; n <= histos.get<TH1>(HIST("hNEvents"))->GetNbinsX(); n++) {
      histos.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(n, hNEventsLabels[n - 1]);
    }
    histos.add("hEventVertexZ", "hEventVertexZ", kTH1F, {{120, -12., 12.}});
    histos.add("hEventCentrality", "hEventCentrality", kTH1F, {{101, 0, 101}});
    histos.add("hEventCentralityBefEvSel", "hEventCentralityBefEvSel", kTH1F, {{101, 0, 101}});
    histos.add("hEventCentralityBefEPSel", "hEventCentralityBefEPSel", kTH1F, {{101, 0, 101}});
    histos.add("hEventCentralityT0M", "hEventCentralityT0M", kTH1F, {{101, 0, 101}});
    histos.add("hEventCentralityBefEvSelT0M", "hEventCentralityBefEvSelT0M", kTH1F, {{101, 0, 101}});
    histos.add("hEventCentralityBefEPSelT0M", "hEventCentralityBefEPSelT0M", kTH1F, {{101, 0, 101}});
    histos.add("hPsiT0C", "hPsiT0C", HistType::kTH1D, {{100, -o2::constants::math::PI, o2::constants::math::PI}});
    histos.add("hPsiT0CvsCentFT0C", "hPsiT0CvsCentFT0C", HistType::kTH2D, {CentAxis, {100, -o2::constants::math::PI, o2::constants::math::PI}});
    histos.add("Psi_EP_FT0C_shifted", "Psi_EP_FT0C_shifted", HistType::kTH2D, {CentAxis, {100, -o2::constants::math::PI, o2::constants::math::PI}});
    histos.add("Psi_EP_FV0A_shifted", "Psi_EP_FT0C_shifted", HistType::kTH2D, {CentAxis, {100, -o2::constants::math::PI, o2::constants::math::PI}});
    histos.add("Psi_EP_FT0A_shifted", "Psi_EP_FT0C_shifted", HistType::kTH2D, {CentAxis, {100, -o2::constants::math::PI, o2::constants::math::PI}});
    histos.add("Psi_EP_TPCA_shifted", "Psi_EP_FT0C_shifted", HistType::kTH2D, {CentAxis, {100, -o2::constants::math::PI, o2::constants::math::PI}});
    histos.add("Psi_EP_TPCC_shifted", "Psi_EP_FT0C_shifted", HistType::kTH2D, {CentAxis, {100, -o2::constants::math::PI, o2::constants::math::PI}});
    histos.add("hPsiZDCA_vs_ZDCC", "hPsiZDCA_vs_ZDCC", HistType::kTH2D, {{100, -o2::constants::math::PI, o2::constants::math::PI}, {100, -o2::constants::math::PI, o2::constants::math::PI}});
    histos.add("hEventNchCorrelation", "hEventNchCorrelation", kTH2F, {{5000, 0, 5000}, {5000, 0, 2500}});
    histos.add("hEventPVcontributorsVsCentrality", "hEventPVcontributorsVsCentrality", kTH2F, {{100, 0, 100}, {5000, 0, 5000}});
    histos.add("hEventGlobalTracksVsCentrality", "hEventGlobalTracksVsCentrality", kTH2F, {{100, 0, 100}, {2500, 0, 2500}});
    histos.add("hEventNchCorrelationBefCuts", "hEventNchCorrelationBefCuts", kTH2F, {{5000, 0, 5000}, {2500, 0, 2500}});
    histos.add("hEventPVcontributorsVsCentralityBefCuts", "hEventPVcontributorsVsCentralityBefCuts", kTH2F, {{100, 0, 100}, {5000, 0, 5000}});
    histos.add("hEventGlobalTracksVsCentralityBefCuts", "hEventGlobalTracksVsCentralityBefCuts", kTH2F, {{100, 0, 100}, {2500, 0, 2500}});
    histos.add("hEventNchCorrelationAfterEP", "hEventNchCorrelationAfterEP", kTH2F, {{5000, 0, 5000}, {2500, 0, 2500}});
    histos.add("hEventPVcontributorsVsCentralityAfterEP", "hEventPVcontributorsVsCentralityAfterEP", kTH2F, {{100, 0, 100}, {5000, 0, 5000}});
    histos.add("hEventGlobalTracksVsCentralityAfterEP", "hEventGlobalTracksVsCentralityAfterEP", kTH2F, {{100, 0, 100}, {2500, 0, 2500}});
    histos.add("hMultNTracksITSTPCVsCentrality", "hMultNTracksITSTPCVsCentrality", kTH2F, {{100, 0, 100}, {1000, 0, 5000}});

    histos.add("hCandidate", "hCandidate", HistType::kTH1F, {{22, -0.5, 21.5}});
    histos.add("hLambdaCandidate", "hLambdaCandidate", HistType::kTH1F, {{5, -0.5, 4.5}});
    histos.add("hCascadeSignal", "hCascadeSignal", HistType::kTH1F, {{6, -0.5, 5.5}});
    histos.add("hCascade", "hCascade", HistType::kTH1F, {{6, -0.5, 5.5}});
    histos.add("hCascadeDauSel", "hCascadeDauSel", HistType::kTH1F, {{2, -0.5, 1.5}});
    histos.add("hLambdaDauSel", "hLambdaDauSel", HistType::kTH1F, {{3, -0.5, 2.5}});
    histos.add("hALambdaDauSel", "hALambdaDauSel", HistType::kTH1F, {{3, -0.5, 2.5}});
    histos.add("hXiPtvsCent", "hXiPtvsCent", HistType::kTH2F, {{100, 0, 100}, {400, 0, 20}});
    histos.add("hXiPtvsCentEta08", "hXiPtvsCentEta08", HistType::kTH2F, {{100, 0, 100}, {400, 0, 20}});
    histos.add("hXiPtvsCentY05", "hXiPtvsCentY05", HistType::kTH2F, {{100, 0, 100}, {400, 0, 20}});
    histos.add("hOmegaPtvsCent", "hOmegaPtvsCent", HistType::kTH2F, {{100, 0, 100}, {400, 0, 20}});
    histos.add("hOmegaPtvsCentEta08", "hOmegaPtvsCentEta08", HistType::kTH2F, {{100, 0, 100}, {400, 0, 20}});
    histos.add("hOmegaPtvsCentY05", "hOmegaPtvsCentY05", HistType::kTH2F, {{100, 0, 100}, {400, 0, 20}});
    histos.add("hCascadePhi", "hCascadePhi", HistType::kTH1F, {{100, 0, o2::constants::math::TwoPI}});
    histos.add("hcascminuspsiT0C", "hcascminuspsiT0C", HistType::kTH1F, {{100, 0, o2::constants::math::PI}});
    histos.add("hLambdaPhi", "hLambdaPhi", HistType::kTH1F, {{100, 0, o2::constants::math::TwoPI}});
    histos.add("hlambdaminuspsiT0C", "hlambdaminuspsiT0C", HistType::kTH1F, {{100, 0, o2::constants::math::PI}});
    histos.add("hv2CEPvsFT0C", "hv2CEPvsFT0C", HistType::kTH2F, {CentAxis, {100, -1, 1}});
    histos.add("hv2CEPvsv2CSP", "hv2CEPvsV2CSP", HistType::kTH2F, {{100, -1, 1}, {100, -1, 1}});
    histos.add("hv1EPvsv1SP", "hV1EPvsV1SP", HistType::kTH2F, {{100, -1, 1}, {100, -1, 1}});
    histos.add("hv1SP_ZDCA_vs_ZDCC", "hv1SP_ZDCA_vs_ZDCC", HistType::kTH2F, {{100, -1, 1}, {100, -1, 1}});

    const AxisSpec thnAxisFT0C{thnAxisConfigs.thnConfigAxisFT0C, "FT0C (%)"};
    const AxisSpec thnAxisEta{thnAxisConfigs.thnConfigAxisEta, "#eta"};
    const AxisSpec thnAxisPt{thnAxisConfigs.thnConfigAxisPt, "p_{T}"};
    const AxisSpec thnAxisPtLambda{thnAxisConfigs.thnConfigAxisPtLambda, "p_{T, #Lambda}"};
    const AxisSpec thnAxisCharge{thnAxisConfigs.thnConfigAxisCharge, "Charge"};
    const AxisSpec thnAxisPsiDiff{thnAxisConfigs.thnConfigAxisPsiDiff, "2(phi-Psi)"};
    const AxisSpec thnAxisMassXi{thnAxisConfigs.thnConfigAxisMassXi, "inv. mass (#Lambda #pi) (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisMassOmega{thnAxisConfigs.thnConfigAxisMassOmega, "inv. mass (#Lambda K) (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisMassLambda{thnAxisConfigs.thnConfigAxisMassLambda, "inv. mass (p #pi) (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisBDTScore{thnAxisConfigs.thnConfigAxisBDTScore, "BDT score"};
    const AxisSpec thnAxisV2{thnAxisConfigs.thnConfigAxisV2, "v_{2}"};
    const AxisSpec thnAxisPzs2Xi{thnAxisConfigs.thnConfigAxisPzs2Xi, "Pzs2Xi"};
    const AxisSpec thnAxisPzs2Omega{thnAxisConfigs.thnConfigAxisPzs2Omega, "Pzs2Omega"};
    const AxisSpec thnAxisPzs2Lambda{thnAxisConfigs.thnConfigAxisPzs2Lambda, "Pzs2Lambda"};
    const AxisSpec thnAxisCos2Theta{thnAxisConfigs.thnConfigAxisCos2Theta, "Cos2Theta"};
    const AxisSpec thnAxisCos2ThetaL{thnAxisConfigs.thnConfigAxisCos2ThetaL, "Cos2ThetaL"};
    const AxisSpec thnAxisCosThetaXiAlpha{thnAxisConfigs.thnConfigAxisCosThetaXiAlpha, "CosThetaXiWithAlpha"};
    const AxisSpec thnAxisCosThetaOmegaAlpha{thnAxisConfigs.thnConfigAxisCosThetaOmegaAlpha, "CosThetaOmegaWithAlpha"};
    const AxisSpec thnAxisCosThetaProtonAlpha{thnAxisConfigs.thnConfigAxisCosThetaProtonAlpha, "CosThetaProtonWithAlpha"};

    histos.add("hCentvsPtvsPrimaryFracLambda", "hCentvsPtvsPrimaryFracLambda", HistType::kTH3F, {{100, 0, 100}, thnAxisPtLambda, {4, -0.5, 3.5}});
    histos.add("hCentvsPrimaryFracLambda", "hCentvsPrimaryFracLambda", HistType::kTH2F, {{100, 0, 100}, {4, -0.5, 3.5}});
    histos.add("massXi_ProtonAcc", "massXi", HistType::kTH1F, {thnAxisMassXi});
    histos.add("massOmega_ProtonAcc", "massOmega", HistType::kTH1F, {thnAxisMassOmega});

    if (fillingConfigs.isFillTHNXi) {
      if (fillingConfigs.isFillTHN_V2)
        histos.add("hXiV2", "THn for v2 of Xi", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisPt, thnAxisMassXi, thnAxisBDTScore, thnAxisV2});
      if (fillingConfigs.isFillTHN_Pz)
        histos.add("hXiPzs2", "THn for Pzs2 of Xi", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisPt, thnAxisMassXi, thnAxisBDTScore, thnAxisPzs2Xi});
      if (fillingConfigs.isFillTHN_PzFromLambda)
        histos.add("hXiPzs2FromLambda", "THn for Pzs2 of Xi", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisPt, thnAxisMassXi, thnAxisBDTScore, thnAxisPzs2Lambda});
      if (fillingConfigs.isFillTHN_Acc)
        histos.add("hXiCos2Theta", "THn for Cos2Theta of Xi", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisEta, thnAxisPt, thnAxisMassXi, thnAxisBDTScore, thnAxisCos2Theta});
      if (fillingConfigs.isFillTHN_AccFromLambdaVsCasc)
        histos.add("hXiCos2ThetaFromLambda", "THn for Cos2Theta of Lambda vs Xi mass and pt", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisEta, thnAxisPt, thnAxisMassXi, thnAxisBDTScore, thnAxisCos2Theta});
      if (fillingConfigs.isFillTHN_AccFromLambdaVsLambda)
        histos.add("hXiCos2ThetaFromLambdaL", "THn for Cos2Theta of Lambda vs Lambda mass and pt", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisEta, thnAxisPtLambda, thnAxisMassLambda, thnAxisBDTScore, thnAxisCos2ThetaL});
    }
    if (fillingConfigs.isFillTHNXi_PzVsPsi) {
      if (fillingConfigs.isFillTHN_Pz)
        histos.add("hXiPzVsPsi", "THn for cosTheta of Xi", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisPt, thnAxisMassXi, thnAxisBDTScore, thnAxisCosThetaXiAlpha, thnAxisPsiDiff});
      if (fillingConfigs.isFillTHN_PzFromLambda)
        histos.add("hXiPzVsPsiFromLambda", "THn for cosTheta of Xi", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisPt, thnAxisMassXi, thnAxisBDTScore, thnAxisCosThetaProtonAlpha, thnAxisPsiDiff});
      if (fillingConfigs.isFillTHN_Acc)
        histos.add("hXiCos2ThetaVsPsi", "THn for cos2Theta of Xi", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisEta, thnAxisPt, thnAxisMassXi, thnAxisBDTScore, thnAxisCos2Theta, thnAxisPsiDiff});
      if (fillingConfigs.isFillTHN_AccFromLambdaVsCasc)
        histos.add("hXiCos2ThetaVsPsiFromLambda", "THn for cos2Theta of Lambda vs Xi mass and pt", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisEta, thnAxisPt, thnAxisMassXi, thnAxisBDTScore, thnAxisCos2Theta, thnAxisPsiDiff});
      if (fillingConfigs.isFillTHN_AccFromLambdaVsLambda)
        histos.add("hXiCos2ThetaVsPsiFromLambdaL", "THn for cos2Theta of Lambda vs Lambda mass and pt", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisEta, thnAxisPtLambda, thnAxisMassLambda, thnAxisBDTScore, thnAxisCos2ThetaL, thnAxisPsiDiff});
    }
    if (fillingConfigs.isFillTHNOmega) {
      if (fillingConfigs.isFillTHN_V2)
        histos.add("hOmegaV2", "THn for v2 of Omega", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisPt, thnAxisMassOmega, thnAxisBDTScore, thnAxisV2});
      if (fillingConfigs.isFillTHN_Pz)
        histos.add("hOmegaPzs2", "THn for Pzs2 of Omega", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisPt, thnAxisMassOmega, thnAxisBDTScore, thnAxisPzs2Omega});
      if (fillingConfigs.isFillTHN_PzFromLambda)
        histos.add("hOmegaPzs2FromLambda", "THn for Pzs2 of Omega", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisPt, thnAxisMassOmega, thnAxisBDTScore, thnAxisPzs2Lambda});
      if (fillingConfigs.isFillTHN_Acc)
        histos.add("hOmegaCos2Theta", "THn for Cos2Theta of Omega", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisEta, thnAxisPt, thnAxisMassOmega, thnAxisBDTScore, thnAxisCos2Theta});
      if (fillingConfigs.isFillTHN_AccFromLambdaVsCasc)
        histos.add("hOmegaCos2ThetaFromLambda", "THn for Cos2Theta of Lambda vs Omega mass and pt", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisEta, thnAxisPt, thnAxisMassOmega, thnAxisBDTScore, thnAxisCos2Theta});
      if (fillingConfigs.isFillTHN_AccFromLambdaVsLambda)
        histos.add("hOmegaCos2ThetaFromLambdaL", "THn for Cos2Theta of Lambda vs Lambda mass and pt", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisEta, thnAxisPtLambda, thnAxisMassLambda, thnAxisBDTScore, thnAxisCos2ThetaL});
    }
    if (fillingConfigs.isFillTHNOmega_PzVsPsi) {
      if (fillingConfigs.isFillTHN_Pz)
        histos.add("hOmegaPzVsPsi", "THn for cosTheta of Omega", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisPt, thnAxisMassOmega, thnAxisBDTScore, thnAxisCosThetaOmegaAlpha, thnAxisPsiDiff});
      if (fillingConfigs.isFillTHN_PzFromLambda)
        histos.add("hOmegaPzVsPsiFromLambda", "THn for cosTheta of Omega", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisPt, thnAxisMassOmega, thnAxisBDTScore, thnAxisCosThetaProtonAlpha, thnAxisPsiDiff});
      if (fillingConfigs.isFillTHN_Acc)
        histos.add("hOmegaCos2ThetaVsPsi", "THn for cos2Theta of Omega", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisEta, thnAxisPt, thnAxisMassOmega, thnAxisBDTScore, thnAxisCos2Theta, thnAxisPsiDiff});
      if (fillingConfigs.isFillTHN_AccFromLambdaVsCasc)
        histos.add("hOmegaCos2ThetaVsPsiFromLambda", "THn for cos2Theta of Lambda vs Omega mass and pt", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisEta, thnAxisPt, thnAxisMassOmega, thnAxisBDTScore, thnAxisCos2Theta, thnAxisPsiDiff});
      if (fillingConfigs.isFillTHN_AccFromLambdaVsLambda)
        histos.add("hOmegaCos2ThetaVsPsiFromLambdaL", "THn for cos2Theta of Lambda vs Lambda mass and pt", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisEta, thnAxisPtLambda, thnAxisMassLambda, thnAxisBDTScore, thnAxisCos2ThetaL, thnAxisPsiDiff});
    }
    if (fillingConfigs.isFillTHNLambda) {
      if (fillingConfigs.isFillTHN_V2)
        histos.add("hLambdaV2", "THn for v2 of Lambda", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisPtLambda, thnAxisMassLambda, thnAxisV2});
      if (fillingConfigs.isFillTHN_Pz)
        histos.add("hLambdaPzs2", "THn for Pzs2 of Lambda", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisPtLambda, thnAxisMassLambda, thnAxisPzs2Lambda});
      if (fillingConfigs.isFillTHN_Acc)
        histos.add("hLambdaCos2Theta", "THn for Cos2Theta of Lambda", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisEta, thnAxisPtLambda, thnAxisMassLambda, thnAxisCos2Theta});
    }
    if (fillingConfigs.isFillTHNLambda_PzVsPsi) {
      if (fillingConfigs.isFillTHN_Pz)
        histos.add("hLambdaPzVsPsi", "THn for cosTheta of Lambda", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisPtLambda, thnAxisMassLambda, thnAxisCosThetaProtonAlpha, thnAxisPsiDiff});
      if (fillingConfigs.isFillTHN_Acc)
        histos.add("hLambdaCos2ThetaVsPsi", "THn for cos2Theta of Lambda", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisEta, thnAxisPtLambda, thnAxisMassLambda, thnAxisCos2Theta, thnAxisPsiDiff});
    }

    histosMCGen.add("h2DGenXiEta08", "h2DGenXiEta08", HistType::kTH2F, {{100, 0, 100}, {400, 0, 20}});
    histosMCGen.add("h2DGenOmegaEta08", "h2DGenOmegaEta08", HistType::kTH2F, {{100, 0, 100}, {400, 0, 20}});
    histosMCGen.add("h2DGenXiY05", "h2DGenXiY05", HistType::kTH2F, {{100, 0, 100}, {400, 0, 20}});
    histosMCGen.add("h2DGenOmegaY05", "h2DGenOmegaY05", HistType::kTH2F, {{100, 0, 100}, {400, 0, 20}});
    histosMCGen.add("hGenXiY", "hGenXiY", HistType::kTH1F, {{100, -1, 1}});
    histosMCGen.add("hGenOmegaY", "hGenOmegaY", HistType::kTH1F, {{100, -1, 1}});
    histosMCGen.add("hZvertexGen", "hZvertexGen", HistType::kTH1F, {{100, -20, 20}});
    histosMCGen.add("hNEventsMC", "hNEventsMC", {HistType::kTH1F, {{6, 0.f, 6.f}}});
    for (Int_t n = 1; n <= histosMCGen.get<TH1>(HIST("hNEventsMC"))->GetNbinsX(); n++) {
      histosMCGen.get<TH1>(HIST("hNEventsMC"))->GetXaxis()->SetBinLabel(n, hNEventsLabelsMC[n - 1]);
    }
    histosMCGen.add("hNCascGen", "hNCascGen", {HistType::kTH1F, {{8, 0.f, 8.f}}});
    for (Int_t n = 1; n <= histosMCGen.get<TH1>(HIST("hNCascGen"))->GetNbinsX(); n++) {
      histosMCGen.get<TH1>(HIST("hNCascGen"))->GetXaxis()->SetBinLabel(n, hNCascLabelsMC[n - 1]);
    }

    for (int iS{0}; iS < nParticles; ++iS) {
      cascadev2::hMassBeforeSelVsPt[iS] = histos.add<TH2>(Form("hMassBeforeSelVsPt%s", cascadev2::speciesNames[iS].data()), "hMassBeforeSelVsPt", HistType::kTH2F, {massCascAxis[iS], ptAxisCasc});
      cascadev2::hMassAfterSelVsPt[iS] = histos.add<TH2>(Form("hMassAfterSelVsPt%s", cascadev2::speciesNames[iS].data()), "hMassAfterSelVsPt", HistType::kTH2F, {massCascAxis[iS], ptAxisCasc});
      cascadev2::hSignalScoreBeforeSel[iS] = histos.add<TH1>(Form("hSignalScoreBeforeSel%s", cascadev2::speciesNames[iS].data()), "Signal score before selection;BDT first score;entries", HistType::kTH1F, {{100, 0., 1.}});
      cascadev2::hBkgScoreBeforeSel[iS] = histos.add<TH1>(Form("hBkgScoreBeforeSel%s", cascadev2::speciesNames[iS].data()), "Bkg score before selection;BDT first score;entries", HistType::kTH1F, {{100, 0., 1.}});
      cascadev2::hSignalScoreAfterSel[iS] = histos.add<TH1>(Form("hSignalScoreAfterSel%s", cascadev2::speciesNames[iS].data()), "Signal score after selection;BDT first score;entries", HistType::kTH1F, {{100, 0., 1.}});
      cascadev2::hBkgScoreAfterSel[iS] = histos.add<TH1>(Form("hBkgScoreAfterSel%s", cascadev2::speciesNames[iS].data()), "Bkg score after selection;BDT first score;entries", HistType::kTH1F, {{100, 0., 1.}});
      cascadev2::hSparseV2C[iS] = histos.add<THn>(Form("hSparseV2C%s", cascadev2::speciesNames[iS].data()), "hSparseV2C", HistType::kTHnF, {massCascAxis[iS], ptAxisCasc, v2Axis, CentAxis});
    }
    for (int iS{0}; iS < nCharges; ++iS) {
      lambdav2::hMassBeforeSelVsPt[iS] = histos.add<TH2>(Form("hMassBeforeSelVsPt%s", lambdav2::speciesNames[iS].data()), "hMassBeforeSelVsPt", HistType::kTH2F, {massLambdaAxis[iS], ptAxisLambda});
      lambdav2::hMassAfterSelVsPt[iS] = histos.add<TH2>(Form("hMassAfterSelVsPt%s", lambdav2::speciesNames[iS].data()), "hMassAfterSelVsPt", HistType::kTH2F, {massLambdaAxis[iS], ptAxisLambda});
      lambdav2::hMassAfterSelVsPtTrue[iS] = histos.add<TH2>(Form("hMassAfterSelVsPtTrue%s", lambdav2::speciesNames[iS].data()), "hMassAfterSelVsPtTrue", HistType::kTH2F, {massLambdaAxis[iS], ptAxisLambda});
    }

    if (isApplyML) {
      // Configure and initialise the ML class
      mlResponseXi.configure(binsPtMl, cutsMl, cutDirMl, nClassesMl);
      mlResponseOmega.configure(binsPtMl, cutsMl, cutDirMl, nClassesMl);
      // Bonus: retrieve the model from CCDB (needed for ML application on the GRID)
      if (ccdbConfigs.loadModelsFromCCDB) {
        ccdbApi.init(ccdbConfigs.ccdbUrl);
        mlResponseXi.setModelPathsCCDB(ccdbConfigs.onnxFileNamesXi, ccdbApi, ccdbConfigs.modelPathsCCDBXi, ccdbConfigs.timestampCCDB);
        mlResponseOmega.setModelPathsCCDB(ccdbConfigs.onnxFileNamesOmega, ccdbApi, ccdbConfigs.modelPathsCCDBOmega, ccdbConfigs.timestampCCDB); // TODO: use different model for Xi and Omega
      } else {
        mlResponseXi.setModelPathsLocal(ccdbConfigs.onnxFileNamesXi);
        mlResponseOmega.setModelPathsLocal(ccdbConfigs.onnxFileNamesOmega);
      }
      mlResponseXi.init();
      mlResponseOmega.init();
    }
    if (applyAcceptanceCorrection) {
      ccdb->setURL(ccdbConfigs.ccdbUrl);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      ccdb->setFatalWhenNull(false);
      initAcceptanceFromCCDB();
    }
    if (applyResoCorrection) {
      ccdb->setURL(ccdbConfigs.ccdbUrl);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      ccdb->setFatalWhenNull(false);
      initResoFromCCDB();
    }
    if (applyCentWeightCorrection) {
      ccdb->setURL(ccdbConfigs.ccdbUrl);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      ccdb->setFatalWhenNull(false);
      initCentWeightFromCCDB();
    }
  }

  void processTrainingBackground(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels>::iterator const& coll, soa::Join<aod::CascCollRefs, aod::CascCores, aod::CascExtras, aod::CascBBs> const& Cascades, DauTracks const&)
  {

    int counter = 0;

    if (!AcceptEvent(coll, 1)) {
      return;
    }

    histos.fill(HIST("hEventCentrality"), coll.centFT0C());
    histos.fill(HIST("hEventVertexZ"), coll.posZ());

    for (auto const& casc : Cascades) {
      if (gRandom->Uniform() > downsample) {
        continue;
      }

      float sigmaRangeXi[2]{getNsigmaMass(cascadev2::Xi, casc.pt(), sideBandStart), getNsigmaMass(cascadev2::Xi, casc.pt(), sideBandEnd)};
      float sigmaRangeOmega[2]{getNsigmaMass(cascadev2::Omega, casc.pt(), sideBandStart), getNsigmaMass(cascadev2::Omega, casc.pt(), sideBandEnd)};

      if ((std::abs(casc.mXi() - constants::physics::MassXiMinus) < sigmaRangeXi[0] ||
           std::abs(casc.mXi() - constants::physics::MassXiMinus) > sigmaRangeXi[1]) &&
          (std::abs(casc.mOmega() - constants::physics::MassOmegaMinus) < sigmaRangeOmega[0] ||
           std::abs(casc.mOmega() - constants::physics::MassOmegaMinus) > sigmaRangeOmega[1])) {
        continue;
      }

      /// Add some minimal cuts for single track variables (min number of TPC clusters)
      auto negExtra = casc.negTrackExtra_as<DauTracks>();
      auto posExtra = casc.posTrackExtra_as<DauTracks>();
      auto bachExtra = casc.bachTrackExtra_as<DauTracks>();

      if (doNTPCSigmaCut) {
        if (casc.sign() < 0) {
          if (std::abs(posExtra.tpcNSigmaPr()) > nsigmatpcPr || std::abs(negExtra.tpcNSigmaPi()) > nsigmatpcPi)
            continue;
          histos.fill(HIST("hCandidate"), ++counter);
        } else if (casc.sign() > 0) {
          if (std::abs(posExtra.tpcNSigmaPi()) > nsigmatpcPi || std::abs(negExtra.tpcNSigmaPr()) > nsigmatpcPr)
            continue;
          histos.fill(HIST("hCandidate"), ++counter);
        }
      } else {
        ++counter;
      }
      if (posExtra.tpcCrossedRows() < mintpccrrows || negExtra.tpcCrossedRows() < mintpccrrows || bachExtra.tpcCrossedRows() < mintpccrrows)
        continue;
      histos.fill(HIST("hCandidate"), ++counter);

      fillTrainingTable(coll, casc, 0);
    }
  }

  void processTrainingSignal(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels>::iterator const& coll, CascMCCandidates const& Cascades, DauTracks const&, soa::Join<aod::CascMCCores, aod::CascMCCollRefs> const&)
  {

    if (!AcceptEvent(coll, 1)) {
      return;
    }
    histos.fill(HIST("hEventCentrality"), coll.centFT0C());
    histos.fill(HIST("hEventVertexZ"), coll.posZ());

    for (auto const& casc : Cascades) {
      if (!casc.has_cascMCCore())
        continue;

      auto cascMC = casc.cascMCCore_as<soa::Join<aod::CascMCCores, aod::CascMCCollRefs>>();
      int pdgCode{cascMC.pdgCode()};
      if (!(std::abs(pdgCode) == PDG_t::kXiMinus && std::abs(cascMC.pdgCodeV0()) == PDG_t::kLambda0 && std::abs(cascMC.pdgCodeBachelor()) == PDG_t::kPiPlus) // Xi
          && !(std::abs(pdgCode) == PDG_t::kOmegaMinus && std::abs(cascMC.pdgCodeV0()) == PDG_t::kLambda0 && std::abs(cascMC.pdgCodeBachelor()) == kKPlus))  // Omega
        continue;

      auto negExtra = casc.negTrackExtra_as<DauTracks>();
      auto posExtra = casc.posTrackExtra_as<DauTracks>();
      auto bachExtra = casc.bachTrackExtra_as<DauTracks>();

      int counter = 0;
      bool isCascCandidate = 0;
      isCascCandidate = IsCascAccepted(casc, negExtra, posExtra, bachExtra, counter);
      histos.fill(HIST("hCascadeSignal"), counter);

      // PDG cascades
      if (isCascCandidate)
        fillTrainingTable(coll, casc, pdgCode); // I only store cascades that passed PID and track quality selections
    }
  }

  void processAnalyseData(CollEventAndSpecPlane const& coll, CascCandidates const& Cascades, DauTracks const&)
  {

    if (!AcceptEvent(coll, 1)) {
      return;
    }

    // select only events used for the calibration of the event plane
    if (isGoodEventEP && !coll.triggereventep()) {
      return;
    }

    // event has event plane
    bool hasEventPlane = 0;
    if (coll.triggereventep())
      hasEventPlane = 1;

    // event has spectator plane
    bool hasSpectatorPlane = 0;
    if (coll.triggereventsp())
      hasSpectatorPlane = 1;

    histos.fill(HIST("hNEvents"), 9.5);
    histos.fill(HIST("hEventNchCorrelationAfterEP"), coll.multNTracksPVeta1(), coll.multNTracksGlobal());
    histos.fill(HIST("hEventPVcontributorsVsCentralityAfterEP"), coll.centFT0C(), coll.multNTracksPVeta1());
    histos.fill(HIST("hEventGlobalTracksVsCentralityAfterEP"), coll.centFT0C(), coll.multNTracksGlobal());

    histos.fill(HIST("hEventCentrality"), coll.centFT0C());
    histos.fill(HIST("hEventVertexZ"), coll.posZ());

    ROOT::Math::XYZVector eventplaneVecT0C{coll.qvecFT0CRe(), coll.qvecFT0CIm(), 0};
    ROOT::Math::XYZVector eventplaneVecTPCA{coll.qvecBPosRe(), coll.qvecBPosIm(), 0};
    ROOT::Math::XYZVector eventplaneVecTPCC{coll.qvecBNegRe(), coll.qvecBNegIm(), 0};
    ROOT::Math::XYZVector spectatorplaneVecZDCA{std::cos(coll.psiZDCA()), std::sin(coll.psiZDCA()), 0}; // eta positive = projectile
    ROOT::Math::XYZVector spectatorplaneVecZDCC{std::cos(coll.psiZDCC()), std::sin(coll.psiZDCC()), 0}; // eta negative = target

    const float psiT0C = std::atan2(coll.qvecFT0CIm(), coll.qvecFT0CRe()) * 0.5f;
    const float psiTPCA = std::atan2(coll.qvecBPosIm(), coll.qvecBPosRe()) * 0.5f;
    const float psiTPCC = std::atan2(coll.qvecBNegIm(), coll.qvecBNegRe()) * 0.5f;
    float psiT0CCorr = psiT0C;

    for (int ishift = 1; ishift <= 10; ishift++) {
      histos.fill(HIST("ShiftFT0C"), coll.centFT0C(), 0.5, ishift - 0.5, std::sin(ishift * 2 * psiT0C));
      histos.fill(HIST("ShiftFT0C"), coll.centFT0C(), 1.5, ishift - 0.5, std::cos(ishift * 2 * psiT0C));

      histos.fill(HIST("ShiftTPCL"), coll.centFT0C(), 0.5, ishift - 0.5, std::sin(ishift * 2 * psiTPCA));
      histos.fill(HIST("ShiftTPCL"), coll.centFT0C(), 1.5, ishift - 0.5, std::cos(ishift * 2 * psiTPCA));

      histos.fill(HIST("ShiftTPCR"), coll.centFT0C(), 0.5, ishift - 0.5, std::sin(ishift * 2 * psiTPCC));
      histos.fill(HIST("ShiftTPCR"), coll.centFT0C(), 1.5, ishift - 0.5, std::cos(ishift * 2 * psiTPCC));
    }

    if (ShiftConfigs.cfgShiftCorr) {
      currentRunNumber = coll.runNumber();
      if (currentRunNumber != lastRunNumber) {
        fullCCDBShiftCorrPathFT0C = ShiftConfigs.cfgShiftPathFT0C;
        fullCCDBShiftCorrPathTPCL = ShiftConfigs.cfgShiftPathTPCL;
        fullCCDBShiftCorrPathTPCR = ShiftConfigs.cfgShiftPathTPCR;
        shiftprofileFT0C = ccdb->getForTimeStamp<TProfile3D>(fullCCDBShiftCorrPathFT0C, coll.timestamp());
        shiftprofileTPCL = ccdb->getForTimeStamp<TProfile3D>(fullCCDBShiftCorrPathTPCL, coll.timestamp());
        shiftprofileTPCR = ccdb->getForTimeStamp<TProfile3D>(fullCCDBShiftCorrPathTPCR, coll.timestamp());
        lastRunNumber = currentRunNumber;
      }
    }

    if (ShiftConfigs.cfgShiftCorr) {
      psiT0CCorr = ApplyShiftCorrection(coll, psiT0C, shiftprofileFT0C);
      ComputeEPResolutionwShifts(coll, psiT0C, psiT0C, psiT0C, psiTPCA, psiTPCC, shiftprofileFT0C, shiftprofileTPCL, shiftprofileTPCR, shiftprofileFT0C, shiftprofileFT0C);
    }

    histos.fill(HIST("hPsiT0C"), psiT0CCorr);
    histos.fill(HIST("hPsiZDCA_vs_ZDCC"), coll.psiZDCC(), coll.psiZDCA());
    histos.fill(HIST("hPsiT0CvsCentFT0C"), coll.centFT0C(), psiT0CCorr);

    resolution.fill(HIST("QVectorsT0CTPCA"), eventplaneVecT0C.Dot(eventplaneVecTPCA), coll.centFT0C());
    resolution.fill(HIST("QVectorsT0CTPCC"), eventplaneVecT0C.Dot(eventplaneVecTPCC), coll.centFT0C());
    resolution.fill(HIST("QVectorsTPCAC"), eventplaneVecTPCA.Dot(eventplaneVecTPCC), coll.centFT0C());
    resolution.fill(HIST("EP_T0CTPCA"), std::cos(2 * (psiT0C - psiTPCA)), coll.centFT0C());
    resolution.fill(HIST("EP_T0CTPCC"), std::cos(2 * (psiT0C - psiTPCC)), coll.centFT0C());
    resolution.fill(HIST("EP_TPCAC"), std::cos(2 * (psiTPCA - psiTPCC)), coll.centFT0C());
    resolution.fill(HIST("QVectorsNormT0CTPCA"), eventplaneVecT0C.Dot(eventplaneVecTPCA) / (coll.qTPCR() * coll.sumAmplFT0C()), coll.centFT0C());
    resolution.fill(HIST("QVectorsNormT0CTPCC"), eventplaneVecT0C.Dot(eventplaneVecTPCC) / (coll.qTPCL() * coll.sumAmplFT0C()), coll.centFT0C());
    resolution.fill(HIST("QVectorsNormTPCAC"), eventplaneVecTPCA.Dot(eventplaneVecTPCC) / (coll.qTPCR() * coll.qTPCL()), coll.centFT0C());
    resolution.fill(HIST("QVectorsSpecPlane"), spectatorplaneVecZDCC.Dot(spectatorplaneVecZDCA), coll.centFT0C());

    std::vector<float> bdtScore[nParticles];
    for (auto const& casc : Cascades) {

      /// Add some minimal cuts for single track variables (min number of TPC clusters)
      auto negExtra = casc.negTrackExtra_as<DauTracks>();
      auto posExtra = casc.posTrackExtra_as<DauTracks>();
      auto bachExtra = casc.bachTrackExtra_as<DauTracks>();

      int counter = 0;
      bool isCascCandidate = 0;
      isCascCandidate = IsCascAccepted(casc, negExtra, posExtra, bachExtra, counter);
      histos.fill(HIST("hCascade"), counter);
      histos.fill(HIST("hCascadeDauSel"), (int)isCascCandidate);
      if (!isCascCandidate)
        continue;

      // ML selections
      bool isSelectedCasc[2]{false, false};

      std::vector<float> inputFeaturesCasc{casc.cascradius(),
                                           casc.v0radius(),
                                           casc.casccosPA(coll.posX(), coll.posY(), coll.posZ()),
                                           casc.v0cosPA(coll.posX(), coll.posY(), coll.posZ()),
                                           casc.dcapostopv(),
                                           casc.dcanegtopv(),
                                           casc.dcabachtopv(),
                                           casc.dcacascdaughters(),
                                           casc.dcaV0daughters(),
                                           casc.dcav0topv(coll.posX(), coll.posY(), coll.posZ()),
                                           casc.bachBaryonCosPA(),
                                           casc.bachBaryonDCAxyToPV()};

      float massCasc[2]{casc.mXi(), casc.mOmega()};

      // pt cut
      if (casc.pt() < CandidateConfigs.MinPt || casc.pt() > CandidateConfigs.MaxPt) {
        continue;
      }

      cascadev2::hMassBeforeSelVsPt[0]->Fill(massCasc[0], casc.pt());
      cascadev2::hMassBeforeSelVsPt[1]->Fill(massCasc[1], casc.pt());

      if (isApplyML) {
        // Retrieve model output and selection outcome
        isSelectedCasc[0] = mlResponseXi.isSelectedMl(inputFeaturesCasc, casc.pt(), bdtScore[0]);
        isSelectedCasc[1] = mlResponseOmega.isSelectedMl(inputFeaturesCasc, casc.pt(), bdtScore[1]);

        for (int iS{0}; iS < nParticles; ++iS) {
          // Fill BDT score histograms before selection
          cascadev2::hSignalScoreBeforeSel[iS]->Fill(bdtScore[0][1]);
          cascadev2::hBkgScoreBeforeSel[iS]->Fill(bdtScore[1][0]);

          // Fill histograms for selected candidates
          if (isSelectedCasc[iS]) {
            cascadev2::hSignalScoreAfterSel[iS]->Fill(bdtScore[0][1]);
            cascadev2::hBkgScoreAfterSel[iS]->Fill(bdtScore[1][0]);
            cascadev2::hMassAfterSelVsPt[iS]->Fill(massCasc[iS], casc.pt());
          }
        }
      } else {
        isSelectedCasc[0] = true;
        isSelectedCasc[1] = true;
      }

      ROOT::Math::XYZVector cascQvec{std::cos(2 * casc.phi()), std::sin(2 * casc.phi()), 0};
      auto v2CSP = cascQvec.Dot(eventplaneVecT0C); // not normalised by amplitude
      auto cascminuspsiT0C = GetPhiInRange(casc.phi() - psiT0CCorr);
      auto v2CEP = std::cos(2.0 * cascminuspsiT0C);
      ROOT::Math::XYZVector cascUvec{std::cos(casc.phi()), std::sin(casc.phi()), 0};
      auto v1SP_ZDCA = cascUvec.Dot(spectatorplaneVecZDCA);
      auto v1SP_ZDCC = cascUvec.Dot(spectatorplaneVecZDCC);
      auto v1EP_ZDCA = std::cos(casc.phi() - coll.psiZDCA());
      auto v1EP_ZDCC = std::cos(casc.phi() - coll.psiZDCC());
      float v1SP = 0.5 * (v1SP_ZDCA - v1SP_ZDCC);
      float v1EP = 0.5 * (v1EP_ZDCA - v1EP_ZDCC); // same as v1SP

      // polarization variables
      double masses[2]{o2::constants::physics::MassXiMinus, o2::constants::physics::MassOmegaMinus};
      ROOT::Math::PxPyPzMVector cascadeVector[2], lambdaVector, protonVector;
      float cosThetaStarLambda[2], cosThetaStarProton;
      lambdaVector.SetCoordinates(casc.pxlambda(), casc.pylambda(), casc.pzlambda(), o2::constants::physics::MassLambda);
      ROOT::Math::Boost lambdaBoost{lambdaVector.BoostToCM()};
      if (casc.sign() > 0) {
        protonVector.SetCoordinates(casc.pxneg(), casc.pyneg(), casc.pzneg(), o2::constants::physics::MassProton);
      } else {
        protonVector.SetCoordinates(casc.pxpos(), casc.pypos(), casc.pzpos(), o2::constants::physics::MassProton);
      }
      auto boostedProton{lambdaBoost(protonVector)};
      cosThetaStarProton = boostedProton.Pz() / boostedProton.P();
      for (int i{0}; i < nParticles; ++i) {
        cascadeVector[i].SetCoordinates(casc.px(), casc.py(), casc.pz(), masses[i]);
        ROOT::Math::Boost cascadeBoost{cascadeVector[i].BoostToCM()};
        auto boostedLambda{cascadeBoost(lambdaVector)};
        cosThetaStarLambda[i] = boostedLambda.Pz() / boostedLambda.P();
      }

      double ptLambda = std::sqrt(std::pow(casc.pxlambda(), 2) + std::pow(casc.pylambda(), 2));
      auto etaLambda = RecoDecay::eta(std::array{casc.pxlambda(), casc.pylambda(), casc.pzlambda()});

      // acceptance values if requested
      double meanCos2ThetaLambdaFromXi = 1;
      double meanCos2ThetaLambdaFromOmega = 1;
      double meanCos2ThetaProtonFromLambda = 1;
      if (applyAcceptanceCorrection) {
        if (ptLambda < CandidateConfigs.MinPtLambda || ptLambda > CandidateConfigs.MaxPtLambda) {
          continue;
        }
        if (std::abs(casc.eta()) > CandidateConfigs.etaCasc)
          continue;
        if (std::abs(etaLambda) > CandidateConfigs.etaLambdaMax)
          continue;
        int bin2DXi = hAcceptanceXi->FindBin(casc.pt(), casc.eta());
        int bin2DOmega = hAcceptanceOmega->FindBin(casc.pt(), casc.eta());
        int bin2DLambda = hAcceptanceLambda->FindBin(ptLambda, etaLambda);
        meanCos2ThetaLambdaFromXi = hAcceptanceXi->GetBinContent(bin2DXi);
        meanCos2ThetaLambdaFromOmega = hAcceptanceOmega->GetBinContent(bin2DOmega);
        meanCos2ThetaProtonFromLambda = hAcceptanceLambda->GetBinContent(bin2DLambda);
      }

      int chargeIndex = 0;
      if (casc.sign() > 0)
        chargeIndex = 1;
      double pzs2Xi = cosThetaStarLambda[0] * std::sin(2 * (casc.phi() - psiT0CCorr)) / cascadev2::AlphaXi[chargeIndex] / meanCos2ThetaLambdaFromXi;
      double pzs2Omega = cosThetaStarLambda[1] * std::sin(2 * (casc.phi() - psiT0CCorr)) / cascadev2::AlphaOmega[chargeIndex] / meanCos2ThetaLambdaFromOmega;
      double cos2ThetaXi = cosThetaStarLambda[0] * cosThetaStarLambda[0];
      double cos2ThetaOmega = cosThetaStarLambda[1] * cosThetaStarLambda[1];
      double pzs2LambdaFromCasc = cosThetaStarProton * std::sin(2 * (casc.phi() - psiT0CCorr)) / cascadev2::AlphaLambda[chargeIndex] / meanCos2ThetaProtonFromLambda;
      double cos2ThetaLambda = cosThetaStarProton * cosThetaStarProton;

      double cosThetaXiWithAlpha = cosThetaStarLambda[0] / cascadev2::AlphaXi[chargeIndex];
      double cosThetaOmegaWithAlpha = cosThetaStarLambda[1] / cascadev2::AlphaOmega[chargeIndex];
      double cosThetaProtonWithAlpha = cosThetaStarProton / cascadev2::AlphaLambda[chargeIndex];

      histos.fill(HIST("hv2CEPvsFT0C"), coll.centFT0C(), v2CEP);
      histos.fill(HIST("hv2CEPvsv2CSP"), v2CSP, v2CEP);
      histos.fill(HIST("hv1EPvsv1SP"), v1SP, v1EP);
      histos.fill(HIST("hv1SP_ZDCA_vs_ZDCC"), v1SP_ZDCC, v1SP_ZDCA);
      histos.fill(HIST("hCascadePhi"), casc.phi());
      histos.fill(HIST("hcascminuspsiT0C"), cascminuspsiT0C);

      double values[4]{casc.mXi(), casc.pt(), v2CSP, coll.centFT0C()};
      if (isSelectedCasc[0]) {
        cascadev2::hSparseV2C[0]->Fill(values);
      }
      if (isSelectedCasc[1]) {
        values[0] = casc.mOmega();
        cascadev2::hSparseV2C[0]->Fill(values);
      }

      float BDTresponse[nParticles]{0.f, 0.f};
      if (isApplyML) {
        BDTresponse[0] = bdtScore[0][1];
        BDTresponse[1] = bdtScore[1][1];
      }

      if (std::abs(casc.eta()) < CandidateConfigs.etaCasc) {
        if (fillingConfigs.isFillTHNXi) {
          if (fillingConfigs.isFillTHN_V2)
            histos.get<THn>(HIST("hXiV2"))->Fill(coll.centFT0C(), chargeIndex, casc.pt(), casc.mXi(), BDTresponse[0], v2CEP);
          if (fillingConfigs.isFillTHN_Pz)
            histos.get<THn>(HIST("hXiPzs2"))->Fill(coll.centFT0C(), chargeIndex, casc.pt(), casc.mXi(), BDTresponse[0], pzs2Xi);
          if (casc.mLambda() > CandidateConfigs.MinLambdaMass && casc.mLambda() < CandidateConfigs.MaxLambdaMass) {
            if (fillingConfigs.isFillTHN_PzFromLambda)
              histos.get<THn>(HIST("hXiPzs2FromLambda"))->Fill(coll.centFT0C(), chargeIndex, casc.pt(), casc.mXi(), BDTresponse[0], pzs2LambdaFromCasc);
          }
          if (fillingConfigs.isFillTHN_Acc)
            histos.get<THn>(HIST("hXiCos2Theta"))->Fill(coll.centFT0C(), chargeIndex, casc.eta(), casc.pt(), casc.mXi(), BDTresponse[0], cos2ThetaXi);
          if (fillingConfigs.isFillTHN_AccFromLambdaVsCasc)
            histos.get<THn>(HIST("hXiCos2ThetaFromLambda"))->Fill(coll.centFT0C(), chargeIndex, casc.eta(), casc.pt(), casc.mXi(), BDTresponse[0], cos2ThetaLambda);
          if (casc.mXi() > CandidateConfigs.MinXiMass && casc.mXi() < CandidateConfigs.MaxXiMass) {
            if (fillingConfigs.isFillTHN_AccFromLambdaVsLambda)
              histos.get<THn>(HIST("hXiCos2ThetaFromLambdaL"))->Fill(coll.centFT0C(), chargeIndex, etaLambda, ptLambda, casc.mLambda(), BDTresponse[0], cos2ThetaLambda);
            histos.get<TH1>(HIST("massXi_ProtonAcc"))->Fill(casc.mXi());
          }
        }
        if (fillingConfigs.isFillTHNXi_PzVsPsi) {
          if (fillingConfigs.isFillTHN_Pz)
            histos.get<THn>(HIST("hXiPzVsPsi"))->Fill(coll.centFT0C(), chargeIndex, casc.pt(), casc.mXi(), BDTresponse[0], cosThetaXiWithAlpha, 2 * cascminuspsiT0C);
          if (fillingConfigs.isFillTHN_PzFromLambda)
            histos.get<THn>(HIST("hXiPzVsPsiFromLambda"))->Fill(coll.centFT0C(), chargeIndex, casc.pt(), casc.mXi(), BDTresponse[0], cosThetaProtonWithAlpha, 2 * cascminuspsiT0C);
          if (casc.mLambda() > CandidateConfigs.MinLambdaMass && casc.mLambda() < CandidateConfigs.MaxLambdaMass) {
            if (fillingConfigs.isFillTHN_Acc)
              histos.get<THn>(HIST("hXiCos2ThetaVsPsi"))->Fill(coll.centFT0C(), chargeIndex, casc.eta(), casc.pt(), casc.mXi(), BDTresponse[0], cos2ThetaXi, 2 * cascminuspsiT0C);
          }
          if (fillingConfigs.isFillTHN_AccFromLambdaVsCasc)
            histos.get<THn>(HIST("hXiCos2ThetaVsPsiFromLambda"))->Fill(coll.centFT0C(), chargeIndex, casc.eta(), casc.pt(), casc.mXi(), BDTresponse[0], cos2ThetaLambda, 2 * cascminuspsiT0C);
          if (casc.mXi() > CandidateConfigs.MinXiMass && casc.mXi() < CandidateConfigs.MaxXiMass) {
            if (fillingConfigs.isFillTHN_AccFromLambdaVsLambda)
              histos.get<THn>(HIST("hXiCos2ThetaVsPsiFromLambdaL"))->Fill(coll.centFT0C(), chargeIndex, etaLambda, ptLambda, casc.mLambda(), BDTresponse[0], cos2ThetaLambda, 2 * cascminuspsiT0C);
            histos.get<TH1>(HIST("massXi_ProtonAcc"))->Fill(casc.mXi());
          }
        }
        if (fillingConfigs.isFillTHNOmega) {
          if (fillingConfigs.isFillTHN_V2)
            histos.get<THn>(HIST("hOmegaV2"))->Fill(coll.centFT0C(), chargeIndex, casc.pt(), casc.mOmega(), BDTresponse[1], v2CEP);
          if (fillingConfigs.isFillTHN_Pz)
            histos.get<THn>(HIST("hOmegaPzs2"))->Fill(coll.centFT0C(), chargeIndex, casc.pt(), casc.mOmega(), BDTresponse[1], pzs2Omega);
          if (casc.mLambda() > CandidateConfigs.MinLambdaMass && casc.mLambda() < CandidateConfigs.MaxLambdaMass) {
            if (fillingConfigs.isFillTHN_PzFromLambda)
              histos.get<THn>(HIST("hOmegaPzs2FromLambda"))->Fill(coll.centFT0C(), chargeIndex, casc.pt(), casc.mOmega(), BDTresponse[1], pzs2LambdaFromCasc);
          }
          if (fillingConfigs.isFillTHN_Acc)
            histos.get<THn>(HIST("hOmegaCos2Theta"))->Fill(coll.centFT0C(), chargeIndex, casc.eta(), casc.pt(), casc.mOmega(), BDTresponse[1], cos2ThetaOmega);
          if (fillingConfigs.isFillTHN_AccFromLambdaVsCasc)
            histos.get<THn>(HIST("hOmegaCos2ThetaFromLambda"))->Fill(coll.centFT0C(), chargeIndex, casc.eta(), casc.pt(), casc.mOmega(), BDTresponse[1], cos2ThetaLambda);
          if (casc.mOmega() > CandidateConfigs.MinOmegaMass && casc.mOmega() < CandidateConfigs.MaxOmegaMass) {
            if (fillingConfigs.isFillTHN_AccFromLambdaVsLambda)
              histos.get<THn>(HIST("hOmegaCos2ThetaFromLambdaL"))->Fill(coll.centFT0C(), chargeIndex, etaLambda, ptLambda, casc.mLambda(), BDTresponse[1], cos2ThetaLambda);
            histos.get<TH1>(HIST("massOmega_ProtonAcc"))->Fill(casc.mOmega());
          }
        }
        if (fillingConfigs.isFillTHNOmega_PzVsPsi) {
          if (fillingConfigs.isFillTHN_Pz)
            histos.get<THn>(HIST("hOmegaPzVsPsi"))->Fill(coll.centFT0C(), chargeIndex, casc.pt(), casc.mOmega(), BDTresponse[1], cosThetaOmegaWithAlpha, 2 * cascminuspsiT0C);
          if (fillingConfigs.isFillTHN_PzFromLambda)
            histos.get<THn>(HIST("hOmegaPzVsPsiFromLambda"))->Fill(coll.centFT0C(), chargeIndex, casc.pt(), casc.mOmega(), BDTresponse[1], cosThetaProtonWithAlpha, 2 * cascminuspsiT0C);
          if (casc.mLambda() > CandidateConfigs.MinLambdaMass && casc.mLambda() < CandidateConfigs.MaxLambdaMass) {
            if (fillingConfigs.isFillTHN_Acc)
              histos.get<THn>(HIST("hOmegaCos2ThetaVsPsi"))->Fill(coll.centFT0C(), chargeIndex, casc.eta(), casc.pt(), casc.mOmega(), BDTresponse[1], cos2ThetaOmega, 2 * cascminuspsiT0C);
          }
          if (fillingConfigs.isFillTHN_AccFromLambdaVsCasc)
            histos.get<THn>(HIST("hOmegaCos2ThetaVsPsiFromLambda"))->Fill(coll.centFT0C(), chargeIndex, casc.eta(), casc.pt(), casc.mOmega(), BDTresponse[1], cos2ThetaLambda, 2 * cascminuspsiT0C);
          if (casc.mOmega() > CandidateConfigs.MinOmegaMass && casc.mOmega() < CandidateConfigs.MaxOmegaMass) {
            if (fillingConfigs.isFillTHN_AccFromLambdaVsLambda)
              histos.get<THn>(HIST("hOmegaCos2ThetaVsPsiFromLambdaL"))->Fill(coll.centFT0C(), chargeIndex, etaLambda, ptLambda, casc.mLambda(), BDTresponse[1], cos2ThetaLambda, 2 * cascminuspsiT0C);
            histos.get<TH1>(HIST("massOmega_ProtonAcc"))->Fill(casc.mOmega());
          }
        }
      }

      if (isSelectedCasc[0] || isSelectedCasc[1]) {
        if (fillingConfigs.isFillTree)
          fillAnalysedTable(coll, hasEventPlane, hasSpectatorPlane, casc, v2CSP, v2CEP, v1SP_ZDCA, v1SP_ZDCC, psiT0CCorr, BDTresponse[0], BDTresponse[1], 0);
      }
    }
  }

  void processAnalyseDataEP2CentralFW(CollEventPlaneCentralFWOnlyFT0C const& coll, CascCandidates const& Cascades, DauTracks const&)
  {

    if (!AcceptEvent(coll, 1)) {
      return;
    }

    // select only events used for the calibration of the event plane
    if (isGoodEventEP) {
      if (std::abs(coll.qvecFT0CRe()) > 990 || std::abs(coll.qvecFT0CIm()) > 990 || std::abs(coll.qvecBNegRe()) > 990 || std::abs(coll.qvecBNegIm()) > 990 || std::abs(coll.qvecBPosRe()) > 990 || std::abs(coll.qvecBPosIm()) > 990) {
        return;
      }
    }

    // event has FT0C event plane
    bool hasEventPlane = 0;
    if (std::abs(coll.qvecFT0CRe()) < 990 && std::abs(coll.qvecFT0CIm()) < 990)
      hasEventPlane = 1;

    histos.fill(HIST("hNEvents"), 9.5);
    histos.fill(HIST("hEventNchCorrelationAfterEP"), coll.multNTracksPVeta1(), coll.multNTracksGlobal());
    histos.fill(HIST("hEventPVcontributorsVsCentralityAfterEP"), coll.centFT0C(), coll.multNTracksPVeta1());
    histos.fill(HIST("hEventGlobalTracksVsCentralityAfterEP"), coll.centFT0C(), coll.multNTracksGlobal());

    histos.fill(HIST("hEventCentrality"), coll.centFT0C());
    histos.fill(HIST("hEventVertexZ"), coll.posZ());

    ROOT::Math::XYZVector eventplaneVecT0C{coll.qvecFT0CRe(), coll.qvecFT0CIm(), 0};
    ROOT::Math::XYZVector eventplaneVecTPCA{coll.qvecBPosRe(), coll.qvecBPosIm(), 0};
    ROOT::Math::XYZVector eventplaneVecTPCC{coll.qvecBNegRe(), coll.qvecBNegIm(), 0};

    const float psiT0C = std::atan2(coll.qvecFT0CIm(), coll.qvecFT0CRe()) * 0.5f;
    const float psiTPCA = std::atan2(coll.qvecBPosIm(), coll.qvecBPosRe()) * 0.5f;
    const float psiTPCC = std::atan2(coll.qvecBNegIm(), coll.qvecBNegRe()) * 0.5f;
    float psiT0CCorr = psiT0C;
    for (int ishift = 1; ishift <= 10; ishift++) {
      histos.fill(HIST("ShiftFT0C"), coll.centFT0C(), 0.5, ishift - 0.5, std::sin(ishift * 2 * psiT0C));
      histos.fill(HIST("ShiftFT0C"), coll.centFT0C(), 1.5, ishift - 0.5, std::cos(ishift * 2 * psiT0C));

      histos.fill(HIST("ShiftTPCL"), coll.centFT0C(), 0.5, ishift - 0.5, std::sin(ishift * 2 * psiTPCA));
      histos.fill(HIST("ShiftTPCL"), coll.centFT0C(), 1.5, ishift - 0.5, std::cos(ishift * 2 * psiTPCA));

      histos.fill(HIST("ShiftTPCR"), coll.centFT0C(), 0.5, ishift - 0.5, std::sin(ishift * 2 * psiTPCC));
      histos.fill(HIST("ShiftTPCR"), coll.centFT0C(), 1.5, ishift - 0.5, std::cos(ishift * 2 * psiTPCC));
    }

    if (ShiftConfigs.cfgShiftCorr) {
      currentRunNumber = coll.runNumber();
      if (currentRunNumber != lastRunNumber) {
        fullCCDBShiftCorrPathFT0C = ShiftConfigs.cfgShiftPathFT0C;
        fullCCDBShiftCorrPathTPCL = ShiftConfigs.cfgShiftPathTPCL;
        fullCCDBShiftCorrPathTPCR = ShiftConfigs.cfgShiftPathTPCR;
        shiftprofileFT0C = ccdb->getForTimeStamp<TProfile3D>(fullCCDBShiftCorrPathFT0C, coll.timestamp());
        shiftprofileTPCL = ccdb->getForTimeStamp<TProfile3D>(fullCCDBShiftCorrPathTPCL, coll.timestamp());
        shiftprofileTPCR = ccdb->getForTimeStamp<TProfile3D>(fullCCDBShiftCorrPathTPCR, coll.timestamp());
        lastRunNumber = currentRunNumber;
      }
    }

    if (ShiftConfigs.cfgShiftCorr) {
      psiT0CCorr = ApplyShiftCorrection(coll, psiT0C, shiftprofileFT0C);
      ComputeEPResolutionwShifts(coll, psiT0C, psiT0C, psiT0C, psiTPCA, psiTPCC, shiftprofileFT0C, shiftprofileTPCL, shiftprofileTPCR, shiftprofileFT0C, shiftprofileFT0C);
    }

    histos.fill(HIST("hPsiT0C"), psiT0CCorr);
    histos.fill(HIST("hPsiT0CvsCentFT0C"), coll.centFT0C(), psiT0CCorr);

    resolution.fill(HIST("QVectorsT0CTPCA"), eventplaneVecT0C.Dot(eventplaneVecTPCA), coll.centFT0C());
    resolution.fill(HIST("QVectorsT0CTPCC"), eventplaneVecT0C.Dot(eventplaneVecTPCC), coll.centFT0C());
    resolution.fill(HIST("QVectorsTPCAC"), eventplaneVecTPCA.Dot(eventplaneVecTPCC), coll.centFT0C());
    resolution.fill(HIST("EP_T0CTPCA"), std::cos(2 * (psiT0C - psiTPCA)), coll.centFT0C());
    resolution.fill(HIST("EP_T0CTPCC"), std::cos(2 * (psiT0C - psiTPCC)), coll.centFT0C());
    resolution.fill(HIST("EP_TPCAC"), std::cos(2 * (psiTPCA - psiTPCC)), coll.centFT0C());
    resolution.fill(HIST("QVectorsNormT0CTPCA"), eventplaneVecT0C.Dot(eventplaneVecTPCA) / (coll.qTPCR() * coll.sumAmplFT0C()), coll.centFT0C());
    resolution.fill(HIST("QVectorsNormT0CTPCC"), eventplaneVecT0C.Dot(eventplaneVecTPCC) / (coll.qTPCL() * coll.sumAmplFT0C()), coll.centFT0C());
    resolution.fill(HIST("QVectorsNormTPCAC"), eventplaneVecTPCA.Dot(eventplaneVecTPCC) / (coll.qTPCR() * coll.qTPCL()), coll.centFT0C());

    std::vector<float> bdtScore[nParticles];
    for (auto const& casc : Cascades) {

      /// Add some minimal cuts for single track variables (min number of TPC clusters)
      auto negExtra = casc.negTrackExtra_as<DauTracks>();
      auto posExtra = casc.posTrackExtra_as<DauTracks>();
      auto bachExtra = casc.bachTrackExtra_as<DauTracks>();

      int counter = 0;
      bool isCascCandidate = 0;
      isCascCandidate = IsCascAccepted(casc, negExtra, posExtra, bachExtra, counter);
      histos.fill(HIST("hCascade"), counter);
      histos.fill(HIST("hCascadeDauSel"), (int)isCascCandidate);
      if (!isCascCandidate)
        continue;

      // ML selections
      bool isSelectedCasc[nParticles]{false, false};

      std::vector<float> inputFeaturesCasc{casc.cascradius(),
                                           casc.v0radius(),
                                           casc.casccosPA(coll.posX(), coll.posY(), coll.posZ()),
                                           casc.v0cosPA(coll.posX(), coll.posY(), coll.posZ()),
                                           casc.dcapostopv(),
                                           casc.dcanegtopv(),
                                           casc.dcabachtopv(),
                                           casc.dcacascdaughters(),
                                           casc.dcaV0daughters(),
                                           casc.dcav0topv(coll.posX(), coll.posY(), coll.posZ()),
                                           casc.bachBaryonCosPA(),
                                           casc.bachBaryonDCAxyToPV()};

      float massCasc[nParticles]{casc.mXi(), casc.mOmega()};

      // pt cut
      if (casc.pt() < CandidateConfigs.MinPt || casc.pt() > CandidateConfigs.MaxPt) {
        continue;
      }

      cascadev2::hMassBeforeSelVsPt[0]->Fill(massCasc[0], casc.pt());
      cascadev2::hMassBeforeSelVsPt[1]->Fill(massCasc[1], casc.pt());

      if (isApplyML) {
        // Retrieve model output and selection outcome
        isSelectedCasc[0] = mlResponseXi.isSelectedMl(inputFeaturesCasc, casc.pt(), bdtScore[0]);
        isSelectedCasc[1] = mlResponseOmega.isSelectedMl(inputFeaturesCasc, casc.pt(), bdtScore[1]);

        for (int iS{0}; iS < nParticles; ++iS) {
          // Fill BDT score histograms before selection
          cascadev2::hSignalScoreBeforeSel[iS]->Fill(bdtScore[0][1]);
          cascadev2::hBkgScoreBeforeSel[iS]->Fill(bdtScore[1][0]);

          // Fill histograms for selected candidates
          if (isSelectedCasc[iS]) {
            cascadev2::hSignalScoreAfterSel[iS]->Fill(bdtScore[0][1]);
            cascadev2::hBkgScoreAfterSel[iS]->Fill(bdtScore[1][0]);
            cascadev2::hMassAfterSelVsPt[iS]->Fill(massCasc[iS], casc.pt());
          }
        }
      } else {
        isSelectedCasc[0] = true;
        isSelectedCasc[1] = true;
      }

      ROOT::Math::XYZVector cascQvec{std::cos(2 * casc.phi()), std::sin(2 * casc.phi()), 0};
      auto v2CSP = cascQvec.Dot(eventplaneVecT0C); // not normalised by amplitude
      auto cascminuspsiT0C = GetPhiInRange(casc.phi() - psiT0CCorr);
      auto v2CEP = std::cos(2.0 * cascminuspsiT0C);
      ROOT::Math::XYZVector cascUvec{std::cos(casc.phi()), std::sin(casc.phi()), 0};

      // polarization variables
      double masses[nParticles]{o2::constants::physics::MassXiMinus, o2::constants::physics::MassOmegaMinus};
      ROOT::Math::PxPyPzMVector cascadeVector[nParticles], lambdaVector, protonVector;
      float cosThetaStarLambda[nParticles], cosThetaStarProton;
      lambdaVector.SetCoordinates(casc.pxlambda(), casc.pylambda(), casc.pzlambda(), o2::constants::physics::MassLambda);
      ROOT::Math::Boost lambdaBoost{lambdaVector.BoostToCM()};
      if (casc.sign() > 0) {
        protonVector.SetCoordinates(casc.pxneg(), casc.pyneg(), casc.pzneg(), o2::constants::physics::MassProton);
      } else {
        protonVector.SetCoordinates(casc.pxpos(), casc.pypos(), casc.pzpos(), o2::constants::physics::MassProton);
      }
      auto boostedProton{lambdaBoost(protonVector)};
      cosThetaStarProton = boostedProton.Pz() / boostedProton.P();
      for (int i{0}; i < nParticles; ++i) {
        cascadeVector[i].SetCoordinates(casc.px(), casc.py(), casc.pz(), masses[i]);
        ROOT::Math::Boost cascadeBoost{cascadeVector[i].BoostToCM()};
        auto boostedLambda{cascadeBoost(lambdaVector)};
        cosThetaStarLambda[i] = boostedLambda.Pz() / boostedLambda.P();
      }

      double ptLambda = std::sqrt(std::pow(casc.pxlambda(), 2) + std::pow(casc.pylambda(), 2));
      auto etaLambda = RecoDecay::eta(std::array{casc.pxlambda(), casc.pylambda(), casc.pzlambda()});

      // acceptance values if requested
      double meanCos2ThetaLambdaFromXi = 1;
      double meanCos2ThetaLambdaFromOmega = 1;
      double meanCos2ThetaProtonFromLambda = 1;
      if (applyAcceptanceCorrection) {
        if (ptLambda < CandidateConfigs.MinPtLambda || ptLambda > CandidateConfigs.MaxPtLambda) {
          continue;
        }
        if (std::abs(casc.eta()) > CandidateConfigs.etaCasc)
          continue;
        if (std::abs(etaLambda) > CandidateConfigs.etaLambdaMax)
          continue;
        int bin2DXi = hAcceptanceXi->FindBin(casc.pt(), casc.eta());
        int bin2DOmega = hAcceptanceOmega->FindBin(casc.pt(), casc.eta());
        int bin2DLambda = hAcceptanceLambda->FindBin(ptLambda, etaLambda);
        meanCos2ThetaLambdaFromXi = hAcceptanceXi->GetBinContent(bin2DXi);
        meanCos2ThetaLambdaFromOmega = hAcceptanceOmega->GetBinContent(bin2DOmega);
        meanCos2ThetaProtonFromLambda = hAcceptanceLambda->GetBinContent(bin2DLambda);
      }

      int chargeIndex = 0;
      if (casc.sign() > 0)
        chargeIndex = 1;
      double pzs2Xi = cosThetaStarLambda[0] * std::sin(2 * (casc.phi() - psiT0CCorr)) / cascadev2::AlphaXi[chargeIndex] / meanCos2ThetaLambdaFromXi;
      double pzs2Omega = cosThetaStarLambda[1] * std::sin(2 * (casc.phi() - psiT0CCorr)) / cascadev2::AlphaOmega[chargeIndex] / meanCos2ThetaLambdaFromOmega;
      double cos2ThetaXi = cosThetaStarLambda[0] * cosThetaStarLambda[0];
      double cos2ThetaOmega = cosThetaStarLambda[1] * cosThetaStarLambda[1];
      double pzs2LambdaFromCasc = cosThetaStarProton * std::sin(2 * (casc.phi() - psiT0CCorr)) / cascadev2::AlphaLambda[chargeIndex] / meanCos2ThetaProtonFromLambda;
      double cos2ThetaLambda = cosThetaStarProton * cosThetaStarProton;

      double cosThetaXiWithAlpha = cosThetaStarLambda[0] / cascadev2::AlphaXi[chargeIndex];
      double cosThetaOmegaWithAlpha = cosThetaStarLambda[1] / cascadev2::AlphaOmega[chargeIndex];
      double cosThetaProtonWithAlpha = cosThetaStarProton / cascadev2::AlphaLambda[chargeIndex];

      histos.fill(HIST("hv2CEPvsFT0C"), coll.centFT0C(), v2CEP);
      histos.fill(HIST("hv2CEPvsv2CSP"), v2CSP, v2CEP);
      histos.fill(HIST("hCascadePhi"), casc.phi());
      histos.fill(HIST("hcascminuspsiT0C"), cascminuspsiT0C);

      double values[4]{casc.mXi(), casc.pt(), v2CSP, coll.centFT0C()};
      if (isSelectedCasc[0]) {
        cascadev2::hSparseV2C[0]->Fill(values);
      }
      if (isSelectedCasc[1]) {
        values[0] = casc.mOmega();
        cascadev2::hSparseV2C[0]->Fill(values);
      }

      float BDTresponse[nParticles]{0.f, 0.f};
      if (isApplyML) {
        BDTresponse[0] = bdtScore[0][1];
        BDTresponse[1] = bdtScore[1][1];
      }

      if (std::abs(casc.eta()) < CandidateConfigs.etaCasc) {
        if (fillingConfigs.isFillTHNXi) {
          if (fillingConfigs.isFillTHN_V2)
            histos.get<THn>(HIST("hXiV2"))->Fill(coll.centFT0C(), chargeIndex, casc.pt(), casc.mXi(), BDTresponse[0], v2CEP);
          if (fillingConfigs.isFillTHN_Pz)
            histos.get<THn>(HIST("hXiPzs2"))->Fill(coll.centFT0C(), chargeIndex, casc.pt(), casc.mXi(), BDTresponse[0], pzs2Xi);
          if (casc.mLambda() > CandidateConfigs.MinLambdaMass && casc.mLambda() < CandidateConfigs.MaxLambdaMass) {
            if (fillingConfigs.isFillTHN_PzFromLambda)
              histos.get<THn>(HIST("hXiPzs2FromLambda"))->Fill(coll.centFT0C(), chargeIndex, casc.pt(), casc.mXi(), BDTresponse[0], pzs2LambdaFromCasc);
          }
          if (fillingConfigs.isFillTHN_Acc)
            histos.get<THn>(HIST("hXiCos2Theta"))->Fill(coll.centFT0C(), chargeIndex, casc.eta(), casc.pt(), casc.mXi(), BDTresponse[0], cos2ThetaXi);
          if (fillingConfigs.isFillTHN_AccFromLambdaVsCasc)
            histos.get<THn>(HIST("hXiCos2ThetaFromLambda"))->Fill(coll.centFT0C(), chargeIndex, casc.eta(), casc.pt(), casc.mXi(), BDTresponse[0], cos2ThetaLambda);
          if (casc.mXi() > CandidateConfigs.MinXiMass && casc.mXi() < CandidateConfigs.MaxXiMass) {
            if (fillingConfigs.isFillTHN_AccFromLambdaVsLambda)
              histos.get<THn>(HIST("hXiCos2ThetaFromLambdaL"))->Fill(coll.centFT0C(), chargeIndex, etaLambda, ptLambda, casc.mLambda(), BDTresponse[0], cos2ThetaLambda);
            histos.get<TH1>(HIST("massXi_ProtonAcc"))->Fill(casc.mXi());
          }
        }
        if (fillingConfigs.isFillTHNXi_PzVsPsi) {
          if (fillingConfigs.isFillTHN_Pz)
            histos.get<THn>(HIST("hXiPzVsPsi"))->Fill(coll.centFT0C(), chargeIndex, casc.pt(), casc.mXi(), BDTresponse[0], cosThetaXiWithAlpha, 2 * cascminuspsiT0C);
          if (fillingConfigs.isFillTHN_PzFromLambda)
            histos.get<THn>(HIST("hXiPzVsPsiFromLambda"))->Fill(coll.centFT0C(), chargeIndex, casc.pt(), casc.mXi(), BDTresponse[0], cosThetaProtonWithAlpha, 2 * cascminuspsiT0C);
          if (casc.mLambda() > CandidateConfigs.MinLambdaMass && casc.mLambda() < CandidateConfigs.MaxLambdaMass) {
            if (fillingConfigs.isFillTHN_Acc)
              histos.get<THn>(HIST("hXiCos2ThetaVsPsi"))->Fill(coll.centFT0C(), chargeIndex, casc.eta(), casc.pt(), casc.mXi(), BDTresponse[0], cos2ThetaXi, 2 * cascminuspsiT0C);
          }
          if (fillingConfigs.isFillTHN_AccFromLambdaVsCasc)
            histos.get<THn>(HIST("hXiCos2ThetaVsPsiFromLambda"))->Fill(coll.centFT0C(), chargeIndex, casc.eta(), casc.pt(), casc.mXi(), BDTresponse[0], cos2ThetaLambda, 2 * cascminuspsiT0C);
          if (casc.mXi() > CandidateConfigs.MinXiMass && casc.mXi() < CandidateConfigs.MaxXiMass) {
            if (fillingConfigs.isFillTHN_AccFromLambdaVsLambda)
              histos.get<THn>(HIST("hXiCos2ThetaVsPsiFromLambdaL"))->Fill(coll.centFT0C(), chargeIndex, etaLambda, ptLambda, casc.mLambda(), BDTresponse[0], cos2ThetaLambda, 2 * cascminuspsiT0C);
            histos.get<TH1>(HIST("massXi_ProtonAcc"))->Fill(casc.mXi());
          }
        }
        if (fillingConfigs.isFillTHNOmega) {
          if (fillingConfigs.isFillTHN_V2)
            histos.get<THn>(HIST("hOmegaV2"))->Fill(coll.centFT0C(), chargeIndex, casc.pt(), casc.mOmega(), BDTresponse[1], v2CEP);
          if (fillingConfigs.isFillTHN_Pz)
            histos.get<THn>(HIST("hOmegaPzs2"))->Fill(coll.centFT0C(), chargeIndex, casc.pt(), casc.mOmega(), BDTresponse[1], pzs2Omega);
          if (casc.mLambda() > CandidateConfigs.MinLambdaMass && casc.mLambda() < CandidateConfigs.MaxLambdaMass) {
            if (fillingConfigs.isFillTHN_PzFromLambda)
              histos.get<THn>(HIST("hOmegaPzs2FromLambda"))->Fill(coll.centFT0C(), chargeIndex, casc.pt(), casc.mOmega(), BDTresponse[1], pzs2LambdaFromCasc);
          }
          if (fillingConfigs.isFillTHN_Acc)
            histos.get<THn>(HIST("hOmegaCos2Theta"))->Fill(coll.centFT0C(), chargeIndex, casc.eta(), casc.pt(), casc.mOmega(), BDTresponse[1], cos2ThetaOmega);
          if (fillingConfigs.isFillTHN_AccFromLambdaVsCasc)
            histos.get<THn>(HIST("hOmegaCos2ThetaFromLambda"))->Fill(coll.centFT0C(), chargeIndex, casc.eta(), casc.pt(), casc.mOmega(), BDTresponse[1], cos2ThetaLambda);
          if (casc.mOmega() > CandidateConfigs.MinOmegaMass && casc.mOmega() < CandidateConfigs.MaxOmegaMass) {
            if (fillingConfigs.isFillTHN_AccFromLambdaVsLambda)
              histos.get<THn>(HIST("hOmegaCos2ThetaFromLambdaL"))->Fill(coll.centFT0C(), chargeIndex, etaLambda, ptLambda, casc.mLambda(), BDTresponse[1], cos2ThetaLambda);
            histos.get<TH1>(HIST("massOmega_ProtonAcc"))->Fill(casc.mOmega());
          }
        }
        if (fillingConfigs.isFillTHNOmega_PzVsPsi) {
          if (fillingConfigs.isFillTHN_Pz)
            histos.get<THn>(HIST("hOmegaPzVsPsi"))->Fill(coll.centFT0C(), chargeIndex, casc.pt(), casc.mOmega(), BDTresponse[1], cosThetaOmegaWithAlpha, 2 * cascminuspsiT0C);
          if (fillingConfigs.isFillTHN_PzFromLambda)
            histos.get<THn>(HIST("hOmegaPzVsPsiFromLambda"))->Fill(coll.centFT0C(), chargeIndex, casc.pt(), casc.mOmega(), BDTresponse[1], cosThetaProtonWithAlpha, 2 * cascminuspsiT0C);
          if (casc.mLambda() > CandidateConfigs.MinLambdaMass && casc.mLambda() < CandidateConfigs.MaxLambdaMass) {
            if (fillingConfigs.isFillTHN_Acc)
              histos.get<THn>(HIST("hOmegaCos2ThetaVsPsi"))->Fill(coll.centFT0C(), chargeIndex, casc.eta(), casc.pt(), casc.mOmega(), BDTresponse[1], cos2ThetaOmega, 2 * cascminuspsiT0C);
          }
          if (fillingConfigs.isFillTHN_AccFromLambdaVsCasc)
            histos.get<THn>(HIST("hOmegaCos2ThetaVsPsiFromLambda"))->Fill(coll.centFT0C(), chargeIndex, casc.eta(), casc.pt(), casc.mOmega(), BDTresponse[1], cos2ThetaLambda, 2 * cascminuspsiT0C);
          if (casc.mOmega() > CandidateConfigs.MinOmegaMass && casc.mOmega() < CandidateConfigs.MaxOmegaMass) {
            if (fillingConfigs.isFillTHN_AccFromLambdaVsLambda)
              histos.get<THn>(HIST("hOmegaCos2ThetaVsPsiFromLambdaL"))->Fill(coll.centFT0C(), chargeIndex, etaLambda, ptLambda, casc.mLambda(), BDTresponse[1], cos2ThetaLambda, 2 * cascminuspsiT0C);
            histos.get<TH1>(HIST("massOmega_ProtonAcc"))->Fill(casc.mOmega());
          }
        }
      }

      if (isSelectedCasc[0] || isSelectedCasc[1]) {
        if (fillingConfigs.isFillTree)
          fillAnalysedTable(coll, hasEventPlane, 0, casc, v2CSP, v2CEP, 0, 0, psiT0CCorr, BDTresponse[0], BDTresponse[1], 0);
      }
    }
  }

  void processAnalyseLambdaEP2CentralFW(CollEventPlaneCentralFW const& coll, V0Candidates const& V0s, DauTracks const&)
  {

    histos.fill(HIST("hEventCentralityBefEvSel"), coll.centFT0C());
    histos.fill(HIST("hEventCentralityBefEvSelT0M"), coll.centFT0M());

    Float_t collisionCentrality = 0;
    if (isCollisionCentrality == 0) { // T0C
      collisionCentrality = coll.centFT0C();
    } else if (isCollisionCentrality == 1) { // T0M
      collisionCentrality = coll.centFT0M();
    }

    if (!AcceptEvent(coll, 1)) {
      return;
    }

    double qvecRe = 0;
    double qvecIm = 0;

    if (isQVecT0C) {
      qvecRe = coll.qvecFT0CRe();
      qvecIm = coll.qvecFT0CIm();
    } else if (isQVecT0A) {
      qvecRe = coll.qvecFT0ARe();
      qvecIm = coll.qvecFT0AIm();
    } else if (isQVecT0M) {
      qvecRe = coll.qvecFT0MRe();
      qvecIm = coll.qvecFT0MIm();
    } else if (isQVecV0A) {
      qvecRe = coll.qvecFV0ARe();
      qvecIm = coll.qvecFV0AIm();
    }

    double qvecReV0A = coll.qvecFV0ARe();
    double qvecImV0A = coll.qvecFV0AIm();
    double qvecReT0A = coll.qvecFT0ARe();
    double qvecImT0A = coll.qvecFT0AIm();

    histos.fill(HIST("hEventCentralityBefEPSel"), collisionCentrality);
    histos.fill(HIST("hEventCentralityBefEPSelT0M"), coll.centFT0M());
    // select only events used for the calibration of the event plane
    if (isGoodEventEP) {
      if (std::abs(qvecRe) > 990 || std::abs(qvecIm) > 990 || std::abs(coll.qvecBNegRe()) > 990 || std::abs(coll.qvecBNegIm()) > 990 || std::abs(coll.qvecBPosRe()) > 990 || std::abs(coll.qvecBPosIm()) > 990) {
        return;
      }
    }

    bool hasSpectatorPlane = 0;
    bool hasEventPlane = 1;

    histos.fill(HIST("hNEvents"), 9.5);
    histos.fill(HIST("hEventNchCorrelationAfterEP"), coll.multNTracksPVeta1(), coll.multNTracksGlobal());
    histos.fill(HIST("hEventPVcontributorsVsCentralityAfterEP"), collisionCentrality, coll.multNTracksPVeta1());
    histos.fill(HIST("hEventGlobalTracksVsCentralityAfterEP"), collisionCentrality, coll.multNTracksGlobal());

    histos.fill(HIST("hEventCentrality"), collisionCentrality);
    histos.fill(HIST("hEventCentralityT0M"), coll.centFT0M());
    histos.fill(HIST("hEventVertexZ"), coll.posZ());

    ROOT::Math::XYZVector eventplaneVecT0C{qvecRe, qvecIm, 0};
    ROOT::Math::XYZVector eventplaneVecV0A{qvecReV0A, qvecImV0A, 0};
    ROOT::Math::XYZVector eventplaneVecT0A{qvecReT0A, qvecImT0A, 0};
    ROOT::Math::XYZVector eventplaneVecTPCA{coll.qvecBPosRe(), coll.qvecBPosIm(), 0};
    ROOT::Math::XYZVector eventplaneVecTPCC{coll.qvecBNegRe(), coll.qvecBNegIm(), 0};

    const float psiT0C = std::atan2(qvecIm, qvecRe) * 0.5f;
    const float psiV0A = std::atan2(qvecImV0A, qvecReV0A) * 0.5f;
    const float psiT0A = std::atan2(qvecImT0A, qvecReT0A) * 0.5f;
    const float psiTPCA = std::atan2(coll.qvecBPosIm(), coll.qvecBPosRe()) * 0.5f;
    const float psiTPCC = std::atan2(coll.qvecBNegIm(), coll.qvecBNegRe()) * 0.5f;
    float psiT0CCorr = psiT0C;
    for (int ishift = 1; ishift <= 10; ishift++) {
      histos.fill(HIST("ShiftFT0C"), collisionCentrality, 0.5, ishift - 0.5, std::sin(ishift * 2 * psiT0C));
      histos.fill(HIST("ShiftFT0C"), collisionCentrality, 1.5, ishift - 0.5, std::cos(ishift * 2 * psiT0C));

      histos.fill(HIST("ShiftFV0A"), collisionCentrality, 0.5, ishift - 0.5, std::sin(ishift * 2 * psiV0A));
      histos.fill(HIST("ShiftFV0A"), collisionCentrality, 1.5, ishift - 0.5, std::cos(ishift * 2 * psiV0A));

      histos.fill(HIST("ShiftFT0A"), collisionCentrality, 0.5, ishift - 0.5, std::sin(ishift * 2 * psiT0A));
      histos.fill(HIST("ShiftFT0A"), collisionCentrality, 1.5, ishift - 0.5, std::cos(ishift * 2 * psiT0A));

      histos.fill(HIST("ShiftTPCL"), collisionCentrality, 0.5, ishift - 0.5, std::sin(ishift * 2 * psiTPCA));
      histos.fill(HIST("ShiftTPCL"), collisionCentrality, 1.5, ishift - 0.5, std::cos(ishift * 2 * psiTPCA));

      histos.fill(HIST("ShiftTPCR"), collisionCentrality, 0.5, ishift - 0.5, std::sin(ishift * 2 * psiTPCC));
      histos.fill(HIST("ShiftTPCR"), collisionCentrality, 1.5, ishift - 0.5, std::cos(ishift * 2 * psiTPCC));
    }

    if (ShiftConfigs.cfgShiftCorr) {
      currentRunNumber = coll.runNumber();
      if (currentRunNumber != lastRunNumber) {
        fullCCDBShiftCorrPathFT0C = ShiftConfigs.cfgShiftPathFT0C;
        fullCCDBShiftCorrPathTPCL = ShiftConfigs.cfgShiftPathTPCL;
        fullCCDBShiftCorrPathTPCR = ShiftConfigs.cfgShiftPathTPCR;
        fullCCDBShiftCorrPathFV0A = ShiftConfigs.cfgShiftPathFV0A;
        fullCCDBShiftCorrPathFT0A = ShiftConfigs.cfgShiftPathFT0A;
        shiftprofileFT0C = ccdb->getForTimeStamp<TProfile3D>(fullCCDBShiftCorrPathFT0C, coll.timestamp());
        shiftprofileTPCL = ccdb->getForTimeStamp<TProfile3D>(fullCCDBShiftCorrPathTPCL, coll.timestamp());
        shiftprofileTPCR = ccdb->getForTimeStamp<TProfile3D>(fullCCDBShiftCorrPathTPCR, coll.timestamp());
        shiftprofileFV0A = ccdb->getForTimeStamp<TProfile3D>(fullCCDBShiftCorrPathFV0A, coll.timestamp());
        shiftprofileFT0A = ccdb->getForTimeStamp<TProfile3D>(fullCCDBShiftCorrPathFT0A, coll.timestamp());
        lastRunNumber = currentRunNumber;
      }
    }

    if (ShiftConfigs.cfgShiftCorr) {
      psiT0CCorr = ApplyShiftCorrection(coll, psiT0C, shiftprofileFT0C);
      ComputeEPResolutionwShifts(coll, psiT0C, psiV0A, psiT0A, psiTPCA, psiTPCC, shiftprofileFT0C, shiftprofileTPCL, shiftprofileTPCR, shiftprofileFV0A, shiftprofileFT0A);
    }

    histos.fill(HIST("hPsiT0C"), psiT0CCorr);
    histos.fill(HIST("hPsiT0CvsCentFT0C"), collisionCentrality, psiT0CCorr);

    resolution.fill(HIST("QVectorsT0CTPCA"), eventplaneVecT0C.Dot(eventplaneVecTPCA), collisionCentrality);
    resolution.fill(HIST("QVectorsT0CTPCC"), eventplaneVecT0C.Dot(eventplaneVecTPCC), collisionCentrality);
    resolution.fill(HIST("QVectorsTPCAC"), eventplaneVecTPCA.Dot(eventplaneVecTPCC), collisionCentrality);
    resolution.fill(HIST("QVectorsT0CV0A"), eventplaneVecT0C.Dot(eventplaneVecV0A), collisionCentrality);
    resolution.fill(HIST("QVectorsV0ATPCC"), eventplaneVecV0A.Dot(eventplaneVecTPCC), collisionCentrality);
    resolution.fill(HIST("QVectorsV0ATPCA"), eventplaneVecV0A.Dot(eventplaneVecTPCA), collisionCentrality);
    resolution.fill(HIST("QVectorsT0CT0A"), eventplaneVecT0C.Dot(eventplaneVecT0A), collisionCentrality);
    resolution.fill(HIST("QVectorsT0ATPCC"), eventplaneVecT0A.Dot(eventplaneVecTPCC), collisionCentrality);
    resolution.fill(HIST("QVectorsT0ATPCA"), eventplaneVecT0A.Dot(eventplaneVecTPCA), collisionCentrality);

    resolution.fill(HIST("EP_T0CTPCA"), std::cos(2 * (psiT0C - psiTPCA)), coll.centFT0C());
    resolution.fill(HIST("EP_T0CTPCC"), std::cos(2 * (psiT0C - psiTPCC)), coll.centFT0C());
    resolution.fill(HIST("EP_TPCAC"), std::cos(2 * (psiTPCA - psiTPCC)), coll.centFT0C());
    resolution.fill(HIST("EP_T0CV0A"), std::cos(2 * (psiT0C - psiV0A)), coll.centFT0C());
    resolution.fill(HIST("EP_V0ATPCC"), std::cos(2 * (psiV0A - psiTPCC)), coll.centFT0C());
    resolution.fill(HIST("EP_V0ATPCA"), std::cos(2 * (psiV0A - psiTPCA)), coll.centFT0C());
    resolution.fill(HIST("EP_T0CT0A"), std::cos(2 * (psiT0C - psiT0A)), coll.centFT0C());
    resolution.fill(HIST("EP_T0ATPCC"), std::cos(2 * (psiT0A - psiTPCC)), coll.centFT0C());
    resolution.fill(HIST("EP_T0ATPCA"), std::cos(2 * (psiT0A - psiTPCA)), coll.centFT0C());

    resolution.fill(HIST("QVectorsNormT0CTPCA"), eventplaneVecT0C.Dot(eventplaneVecTPCA) / (coll.qTPCR() * coll.sumAmplFT0C()), collisionCentrality);
    resolution.fill(HIST("QVectorsNormT0CTPCC"), eventplaneVecT0C.Dot(eventplaneVecTPCC) / (coll.qTPCL() * coll.sumAmplFT0C()), collisionCentrality);
    resolution.fill(HIST("QVectorsNormTPCAC"), eventplaneVecTPCA.Dot(eventplaneVecTPCC) / (coll.qTPCR() * coll.qTPCL()), collisionCentrality);
    resolution.fill(HIST("QVectorsNormT0CV0A"), eventplaneVecT0C.Dot(eventplaneVecV0A) / (coll.sumAmplFT0C() * coll.sumAmplFV0A()), collisionCentrality);
    resolution.fill(HIST("QVectorsNormV0ATPCC"), eventplaneVecV0A.Dot(eventplaneVecTPCC) / (coll.qTPCL() * coll.sumAmplFV0A()), collisionCentrality);
    resolution.fill(HIST("QVectorsNormV0ATPCA"), eventplaneVecV0A.Dot(eventplaneVecTPCA) / (coll.qTPCR() * coll.sumAmplFV0A()), collisionCentrality);
    resolution.fill(HIST("QVectorsNormT0CT0A"), eventplaneVecT0C.Dot(eventplaneVecT0A) / (coll.sumAmplFT0C() * coll.sumAmplFT0A()), collisionCentrality);
    resolution.fill(HIST("QVectorsNormT0ATPCC"), eventplaneVecT0A.Dot(eventplaneVecTPCC) / (coll.qTPCL() * coll.sumAmplFT0A()), collisionCentrality);
    resolution.fill(HIST("QVectorsNormT0ATPCA"), eventplaneVecT0A.Dot(eventplaneVecTPCA) / (coll.qTPCR() * coll.sumAmplFT0A()), collisionCentrality);

    double EPresolution = 1;
    if (applyResoCorrection) {
      int centBin = hReso->FindBin(collisionCentrality);
      EPresolution = hReso->GetBinContent(centBin);
    }
    double centWeight = 1;
    if (applyCentWeightCorrection) {
      int centBin = hCentWeight->FindBin(collisionCentrality);
      centWeight = hCentWeight->GetBinContent(centBin);
    }

    std::vector<float> bdtScore[nParticles];
    for (auto const& v0 : V0s) {

      /// Add some minimal cuts for single track variables (min number of TPC clusters)
      auto negExtra = v0.negTrackExtra_as<DauTracks>();
      auto posExtra = v0.posTrackExtra_as<DauTracks>();

      int counterLambda = 0;
      int counterALambda = 0;
      bool isLambdaCandidate = 0;
      bool isALambdaCandidate = 0;
      if (isLambdaAccepted(negExtra, posExtra, counterLambda))
        isLambdaCandidate = 1;
      if (isAntiLambdaAccepted(negExtra, posExtra, counterALambda))
        isALambdaCandidate = 1;
      histos.fill(HIST("hLambdaDauSel"), counterLambda);
      histos.fill(HIST("hALambdaDauSel"), counterALambda);

      // pt cut
      if (v0.pt() < V0Configs.MinPtV0 || v0.pt() > V0Configs.MaxPtV0) {
        continue;
      }

      float massV0[nCharges]{v0.mLambda(), v0.mAntiLambda()};
      lambdav2::hMassBeforeSelVsPt[0]->Fill(massV0[0], v0.pt());
      lambdav2::hMassBeforeSelVsPt[1]->Fill(massV0[1], v0.pt());

      bool isSelectedV0[2]{false, false};
      if (isV0TopoAccepted(v0) && isLambdaCandidate)
        isSelectedV0[0] = true;
      if (isV0TopoAccepted(v0) && isALambdaCandidate)
        isSelectedV0[1] = true;

      int chargeIndex = -1;
      if (isSelectedV0[0] && !isSelectedV0[1]) { // Lambdas
        histos.fill(HIST("hLambdaCandidate"), 0);
        chargeIndex = 0;
      }
      if (isSelectedV0[1] && !isSelectedV0[0]) { // AntiLambdas
        histos.fill(HIST("hLambdaCandidate"), 1);
        chargeIndex = 1;
      }
      if (isSelectedV0[0] && isSelectedV0[1]) {
        histos.fill(HIST("hLambdaCandidate"), 2);
        if (v0.mLambda() > V0Configs.MinMassLambda && v0.mLambda() < V0Configs.MaxMassLambda && v0.mAntiLambda() > V0Configs.MinMassLambda && v0.mAntiLambda() < V0Configs.MaxMassLambda) {
          histos.fill(HIST("hLambdaCandidate"), 3);
          continue; // in case of ambiguity between Lambda and AntiLambda, I skip the particle; checked to be zero in range 1.105 - 1.125
        }
        if (v0.mLambda() > V0Configs.MinMassLambda && v0.mLambda() < V0Configs.MaxMassLambda)
          chargeIndex = 0;
        else if (v0.mAntiLambda() > V0Configs.MinMassLambda && v0.mAntiLambda() < V0Configs.MaxMassLambda)
          chargeIndex = 1;
        else {
          chargeIndex = 2; // these are bkg candidates
          histos.fill(HIST("hLambdaCandidate"), 4);
        }
      }
      if (!isSelectedV0[0] && !isSelectedV0[1])
        continue;

      ROOT::Math::XYZVector lambdaQvec{std::cos(2 * v0.phi()), std::sin(2 * v0.phi()), 0};
      auto v2CSP = lambdaQvec.Dot(eventplaneVecT0C); // not normalised by amplitude
      auto lambdaminuspsiT0C = GetPhiInRange(v0.phi() - psiT0CCorr);
      auto v2CEP = std::cos(2.0 * lambdaminuspsiT0C);
      ROOT::Math::XYZVector lambdaUvec{std::cos(v0.phi()), std::sin(v0.phi()), 0};

      // polarization variables
      double massLambda = o2::constants::physics::MassLambda;
      float cosThetaStarProton[nCharges];
      ROOT::Math::PxPyPzMVector lambdaVector, protonVector[nCharges];
      lambdaVector.SetCoordinates(v0.px(), v0.py(), v0.pz(), massLambda);
      ROOT::Math::Boost lambdaBoost{lambdaVector.BoostToCM()};
      for (int i{0}; i < nCharges; ++i) {
        if (i == 0)
          protonVector[i].SetCoordinates(v0.pxpos(), v0.pypos(), v0.pzpos(), o2::constants::physics::MassProton);
        else
          protonVector[i].SetCoordinates(v0.pxneg(), v0.pyneg(), v0.pzneg(), o2::constants::physics::MassProton);
        auto boostedProton{lambdaBoost(protonVector[i])};
        cosThetaStarProton[i] = boostedProton.Pz() / boostedProton.P();
      }

      // acceptance values if requested
      double meanCos2ThetaProtonFromLambda = 1;
      if (applyAcceptanceCorrection) {
        int bin2DLambda = hAcceptancePrimaryLambda->FindBin(v0.pt(), v0.eta());
        meanCos2ThetaProtonFromLambda = hAcceptancePrimaryLambda->GetBinContent(bin2DLambda);
      }

      double pzs2Lambda = 0;
      double cos2ThetaLambda = 0;
      double cosThetaLambda = 0;
      if (chargeIndex == 0) {
        pzs2Lambda = cosThetaStarProton[0] * std::sin(2 * (v0.phi() - psiT0CCorr)) / lambdav2::AlphaLambda[0] / meanCos2ThetaProtonFromLambda / EPresolution;
        cos2ThetaLambda = cosThetaStarProton[0] * cosThetaStarProton[0];
        cosThetaLambda = cosThetaStarProton[0] / cascadev2::AlphaLambda[0] / meanCos2ThetaProtonFromLambda / EPresolution;
      } else if (chargeIndex == 1) {
        pzs2Lambda = cosThetaStarProton[1] * std::sin(2 * (v0.phi() - psiT0CCorr)) / lambdav2::AlphaLambda[1] / meanCos2ThetaProtonFromLambda / EPresolution;
        cos2ThetaLambda = cosThetaStarProton[1] * cosThetaStarProton[1];
        cosThetaLambda = cosThetaStarProton[1] / cascadev2::AlphaLambda[1] / meanCos2ThetaProtonFromLambda / EPresolution;
      } else { // I treat these bkg candidates as Lambdas for the purpose of calculating Pz
        pzs2Lambda = cosThetaStarProton[0] * std::sin(2 * (v0.phi() - psiT0CCorr)) / lambdav2::AlphaLambda[0] / meanCos2ThetaProtonFromLambda / EPresolution;
        cos2ThetaLambda = cosThetaStarProton[0] * cosThetaStarProton[0];
        cosThetaLambda = cosThetaStarProton[0] / cascadev2::AlphaLambda[0] / meanCos2ThetaProtonFromLambda / EPresolution;
      }

      histos.fill(HIST("hv2CEPvsFT0C"), collisionCentrality, v2CEP);
      histos.fill(HIST("hv2CEPvsv2CSP"), v2CSP, v2CEP);
      histos.fill(HIST("hLambdaPhi"), v0.phi());
      histos.fill(HIST("hlambdaminuspsiT0C"), lambdaminuspsiT0C);

      if (fillingConfigs.isFillTHNLambda) {
        if (fillingConfigs.isFillTHN_V2)
          histos.get<THn>(HIST("hLambdaV2"))->Fill(collisionCentrality, chargeIndex, v0.pt(), v0.mLambda(), v2CEP);
        if (fillingConfigs.isFillTHN_Pz) {
          //          histos.get<THn>(HIST("hLambdaPzs2"))->Fill(collisionCentrality, chargeIndex, v0.pt(), v0.mLambda(), pzs2Lambda);
          histos.get<THn>(HIST("hLambdaPzs2"))->Fill(collisionCentrality, chargeIndex, v0.pt(), v0.mLambda(), pzs2Lambda, centWeight);
        }
        if (fillingConfigs.isFillTHN_Acc)
          histos.get<THn>(HIST("hLambdaCos2Theta"))->Fill(collisionCentrality, chargeIndex, v0.eta(), v0.pt(), v0.mLambda(), cos2ThetaLambda);
      }
      if (fillingConfigs.isFillTHNLambda_PzVsPsi) {
        if (fillingConfigs.isFillTHN_Pz)
          histos.get<THn>(HIST("hLambdaPzVsPsi"))->Fill(collisionCentrality, chargeIndex, v0.pt(), v0.mLambda(), cosThetaLambda, 2 * lambdaminuspsiT0C, centWeight);
        if (fillingConfigs.isFillTHN_Acc)
          histos.get<THn>(HIST("hLambdaCos2ThetaVsPsi"))->Fill(collisionCentrality, chargeIndex, v0.eta(), v0.pt(), v0.mLambda(), cos2ThetaLambda, 2 * lambdaminuspsiT0C);
      }

      double invMassLambda = 0;
      if (chargeIndex == 0)
        invMassLambda = v0.mLambda();
      else if (chargeIndex == 1)
        invMassLambda = v0.mAntiLambda();
      else
        invMassLambda = v0.mLambda();

      // mass selection
      if (invMassLambda < V0Configs.MinMassLambdaInTree || invMassLambda > V0Configs.MaxMassLambdaInTree)
        continue;

      if (fillingConfigs.isFillTree)
        fillAnalysedLambdaTable(coll, hasEventPlane, hasSpectatorPlane, chargeIndex, v0, v2CEP, psiT0CCorr, pzs2Lambda, cos2ThetaLambda, cosThetaLambda);
    }
  }

  void processAnalyseDataEPCentralFW(CollEventAndSpecPlaneCentralFW const& coll, CascCandidates const& Cascades, DauTracks const&)
  {

    if (!AcceptEvent(coll, 1)) {
      return;
    }

    // select only events used for the calibration of the event plane
    if (isGoodEventEP) {
      if (std::abs(coll.qvecFT0CRe()) > 990 || std::abs(coll.qvecFT0CIm()) > 990 || std::abs(coll.qvecBNegRe()) > 990 || std::abs(coll.qvecBNegIm()) > 990 || std::abs(coll.qvecBPosRe()) > 990 || std::abs(coll.qvecBPosIm()) > 990) {
        return;
      }
    }

    // event has FT0C event plane
    bool hasEventPlane = 0;
    if (std::abs(coll.qvecFT0CRe()) < 990 && std::abs(coll.qvecFT0CIm()) < 990)
      hasEventPlane = 1;

    // event has spectator plane
    bool hasSpectatorPlane = 0;
    if (coll.triggereventsp())
      hasSpectatorPlane = 1;

    histos.fill(HIST("hNEvents"), 9.5);
    histos.fill(HIST("hEventNchCorrelationAfterEP"), coll.multNTracksPVeta1(), coll.multNTracksGlobal());
    histos.fill(HIST("hEventPVcontributorsVsCentralityAfterEP"), coll.centFT0C(), coll.multNTracksPVeta1());
    histos.fill(HIST("hEventGlobalTracksVsCentralityAfterEP"), coll.centFT0C(), coll.multNTracksGlobal());

    histos.fill(HIST("hEventCentrality"), coll.centFT0C());
    histos.fill(HIST("hEventVertexZ"), coll.posZ());

    ROOT::Math::XYZVector eventplaneVecT0C{coll.qvecFT0CRe(), coll.qvecFT0CIm(), 0};
    ROOT::Math::XYZVector eventplaneVecTPCA{coll.qvecBPosRe(), coll.qvecBPosIm(), 0};
    ROOT::Math::XYZVector eventplaneVecTPCC{coll.qvecBNegRe(), coll.qvecBNegIm(), 0};
    ROOT::Math::XYZVector spectatorplaneVecZDCA{std::cos(coll.psiZDCA()), std::sin(coll.psiZDCA()), 0}; // eta positive = projectile
    ROOT::Math::XYZVector spectatorplaneVecZDCC{std::cos(coll.psiZDCC()), std::sin(coll.psiZDCC()), 0}; // eta negative = target

    float NormQvT0C = std::sqrt(eventplaneVecT0C.Dot(eventplaneVecT0C));
    float NormQvTPCA = std::sqrt(eventplaneVecTPCA.Dot(eventplaneVecTPCA));
    float NormQvTPCC = std::sqrt(eventplaneVecTPCC.Dot(eventplaneVecTPCC));

    const float psiT0C = std::atan2(coll.qvecFT0CIm(), coll.qvecFT0CRe()) * 0.5f;
    const float psiTPCA = std::atan2(coll.qvecBPosIm(), coll.qvecBPosRe()) * 0.5f;
    const float psiTPCC = std::atan2(coll.qvecBNegIm(), coll.qvecBNegRe()) * 0.5f;
    float psiT0CCorr = psiT0C;
    for (int ishift = 1; ishift <= 10; ishift++) {
      histos.fill(HIST("ShiftFT0C"), coll.centFT0C(), 0.5, ishift - 0.5, std::sin(ishift * 2 * psiT0C));
      histos.fill(HIST("ShiftFT0C"), coll.centFT0C(), 1.5, ishift - 0.5, std::cos(ishift * 2 * psiT0C));

      histos.fill(HIST("ShiftTPCL"), coll.centFT0C(), 0.5, ishift - 0.5, std::sin(ishift * 2 * psiTPCA));
      histos.fill(HIST("ShiftTPCL"), coll.centFT0C(), 1.5, ishift - 0.5, std::cos(ishift * 2 * psiTPCA));

      histos.fill(HIST("ShiftTPCR"), coll.centFT0C(), 0.5, ishift - 0.5, std::sin(ishift * 2 * psiTPCC));
      histos.fill(HIST("ShiftTPCR"), coll.centFT0C(), 1.5, ishift - 0.5, std::cos(ishift * 2 * psiTPCC));
    }
    if (ShiftConfigs.cfgShiftCorr) {
      currentRunNumber = coll.runNumber();
      if (currentRunNumber != lastRunNumber) {
        fullCCDBShiftCorrPathFT0C = ShiftConfigs.cfgShiftPathFT0C;
        fullCCDBShiftCorrPathTPCL = ShiftConfigs.cfgShiftPathTPCL;
        fullCCDBShiftCorrPathTPCR = ShiftConfigs.cfgShiftPathTPCR;
        shiftprofileFT0C = ccdb->getForTimeStamp<TProfile3D>(fullCCDBShiftCorrPathFT0C, coll.timestamp());
        shiftprofileTPCL = ccdb->getForTimeStamp<TProfile3D>(fullCCDBShiftCorrPathTPCL, coll.timestamp());
        shiftprofileTPCR = ccdb->getForTimeStamp<TProfile3D>(fullCCDBShiftCorrPathTPCR, coll.timestamp());
        lastRunNumber = currentRunNumber;
      }
    }

    if (ShiftConfigs.cfgShiftCorr) {
      psiT0CCorr = ApplyShiftCorrection(coll, psiT0C, shiftprofileFT0C);
      ComputeEPResolutionwShifts(coll, psiT0C, psiT0C, psiT0C, psiTPCA, psiTPCC, shiftprofileFT0C, shiftprofileTPCL, shiftprofileTPCR, shiftprofileFT0C, shiftprofileFT0C);
    }

    histos.fill(HIST("hpsiT0C"), psiT0CCorr);
    histos.fill(HIST("hpsiT0CvsCentFT0C"), coll.centFT0C(), psiT0CCorr);

    resolution.fill(HIST("QVectorsT0CTPCA"), eventplaneVecT0C.Dot(eventplaneVecTPCA), coll.centFT0C());
    resolution.fill(HIST("QVectorsT0CTPCC"), eventplaneVecT0C.Dot(eventplaneVecTPCC), coll.centFT0C());
    resolution.fill(HIST("QVectorsTPCAC"), eventplaneVecTPCA.Dot(eventplaneVecTPCC), coll.centFT0C());

    resolution.fill(HIST("EP_T0CTPCA"), std::cos(2 * (psiT0C - psiTPCA)), coll.centFT0C());
    resolution.fill(HIST("EP_T0CTPCC"), std::cos(2 * (psiT0C - psiTPCC)), coll.centFT0C());
    resolution.fill(HIST("EP_TPCAC"), std::cos(2 * (psiTPCA - psiTPCC)), coll.centFT0C());

    resolution.fill(HIST("QVectorsNormT0CTPCA"), eventplaneVecT0C.Dot(eventplaneVecTPCA) / (NormQvT0C * NormQvTPCA), coll.centFT0C());
    resolution.fill(HIST("QVectorsNormT0CTPCC"), eventplaneVecT0C.Dot(eventplaneVecTPCC) / (NormQvT0C * NormQvTPCC), coll.centFT0C());
    resolution.fill(HIST("QVectorsNormTPCAC"), eventplaneVecTPCA.Dot(eventplaneVecTPCC) / (NormQvTPCA * NormQvTPCC), coll.centFT0C());
    resolution.fill(HIST("QVectorsSpecPlane"), spectatorplaneVecZDCC.Dot(spectatorplaneVecZDCA), coll.centFT0C());

    std::vector<float> bdtScore[nParticles];
    for (auto const& casc : Cascades) {

      /// Add some minimal cuts for single track variables (min number of TPC clusters)
      auto negExtra = casc.negTrackExtra_as<DauTracks>();
      auto posExtra = casc.posTrackExtra_as<DauTracks>();
      auto bachExtra = casc.bachTrackExtra_as<DauTracks>();

      int counter = 0;
      bool isCascCandidate = 0;
      isCascCandidate = IsCascAccepted(casc, negExtra, posExtra, bachExtra, counter);
      histos.fill(HIST("hCascade"), counter);
      histos.fill(HIST("hCascadeDauSel"), (int)isCascCandidate);
      if (!isCascCandidate)
        continue;

      // ML selections
      bool isSelectedCasc[nParticles]{false, false};

      std::vector<float> inputFeaturesCasc{casc.cascradius(),
                                           casc.v0radius(),
                                           casc.casccosPA(coll.posX(), coll.posY(), coll.posZ()),
                                           casc.v0cosPA(coll.posX(), coll.posY(), coll.posZ()),
                                           casc.dcapostopv(),
                                           casc.dcanegtopv(),
                                           casc.dcabachtopv(),
                                           casc.dcacascdaughters(),
                                           casc.dcaV0daughters(),
                                           casc.dcav0topv(coll.posX(), coll.posY(), coll.posZ()),
                                           casc.bachBaryonCosPA(),
                                           casc.bachBaryonDCAxyToPV()};

      float massCasc[nParticles]{casc.mXi(), casc.mOmega()};

      // inv mass loose cut
      if (casc.pt() < CandidateConfigs.MinPt || casc.pt() > CandidateConfigs.MaxPt) {
        continue;
      }

      cascadev2::hMassBeforeSelVsPt[0]->Fill(massCasc[0], casc.pt());
      cascadev2::hMassBeforeSelVsPt[1]->Fill(massCasc[1], casc.pt());

      if (isApplyML) {
        // Retrieve model output and selection outcome
        isSelectedCasc[0] = mlResponseXi.isSelectedMl(inputFeaturesCasc, casc.pt(), bdtScore[0]);
        isSelectedCasc[1] = mlResponseOmega.isSelectedMl(inputFeaturesCasc, casc.pt(), bdtScore[1]);

        for (int iS{0}; iS < nParticles; ++iS) {
          // Fill BDT score histograms before selection
          cascadev2::hSignalScoreBeforeSel[iS]->Fill(bdtScore[0][1]);
          cascadev2::hBkgScoreBeforeSel[iS]->Fill(bdtScore[1][0]);

          // Fill histograms for selected candidates
          if (isSelectedCasc[iS]) {
            cascadev2::hSignalScoreAfterSel[iS]->Fill(bdtScore[0][1]);
            cascadev2::hBkgScoreAfterSel[iS]->Fill(bdtScore[1][0]);
            cascadev2::hMassAfterSelVsPt[iS]->Fill(massCasc[iS], casc.pt());
          }
        }
      } else {
        isSelectedCasc[0] = true;
        isSelectedCasc[1] = true;
      }

      ROOT::Math::XYZVector cascQvec{std::cos(2 * casc.phi()), std::sin(2 * casc.phi()), 0};
      auto v2CSP = cascQvec.Dot(eventplaneVecT0C);
      auto cascminuspsiT0C = GetPhiInRange(casc.phi() - psiT0CCorr);
      auto v2CEP = std::cos(2.0 * cascminuspsiT0C);
      ROOT::Math::XYZVector cascUvec{std::cos(casc.phi()), std::sin(casc.phi()), 0};
      auto v1SP_ZDCA = cascUvec.Dot(spectatorplaneVecZDCA);
      auto v1SP_ZDCC = cascUvec.Dot(spectatorplaneVecZDCC);
      auto v1EP_ZDCA = std::cos(casc.phi() - coll.psiZDCA());
      auto v1EP_ZDCC = std::cos(casc.phi() - coll.psiZDCC());
      float v1SP = 0.5 * (v1SP_ZDCA - v1SP_ZDCC);
      float v1EP = 0.5 * (v1EP_ZDCA - v1EP_ZDCC); // same as v1SP

      histos.fill(HIST("hv2CEPvsFT0C"), coll.centFT0C(), v2CEP);
      histos.fill(HIST("hv2CEPvsv2CSP"), v2CSP, v2CEP);
      histos.fill(HIST("hv1EPvsv1SP"), v1SP, v1EP);
      histos.fill(HIST("hv1SP_ZDCA_vs_ZDCC"), v1SP_ZDCC, v1SP_ZDCA);
      histos.fill(HIST("hCascadePhi"), casc.phi());
      histos.fill(HIST("hcascminuspsiT0C"), cascminuspsiT0C);
      double values[4]{casc.mXi(), casc.pt(), v2CSP, coll.centFT0C()};
      if (isSelectedCasc[0]) {
        cascadev2::hSparseV2C[0]->Fill(values);
      }
      if (isSelectedCasc[1]) {
        values[0] = casc.mOmega();
        cascadev2::hSparseV2C[0]->Fill(values);
      }

      float BDTresponse[nParticles]{0.f, 0.f};
      if (isApplyML) {
        BDTresponse[0] = bdtScore[0][1];
        BDTresponse[1] = bdtScore[1][1];
      }
      if (isSelectedCasc[0] || isSelectedCasc[1])
        if (fillingConfigs.isFillTree)
          fillAnalysedTable(coll, hasEventPlane, hasSpectatorPlane, casc, v2CSP, v2CEP, v1SP_ZDCA, v1SP_ZDCC, psiT0CCorr, BDTresponse[0], BDTresponse[1], 0);
    }
  }

  void processAnalyseMC(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels>::iterator const& coll, CascMCCandidates const& Cascades, DauTracks const&, soa::Join<aod::CascMCCores, aod::CascMCCollRefs> const&)
  {

    if (!AcceptEvent(coll, 1)) {
      return;
    }

    bool hasEventPlane = 0;     // no info at the moment
    bool hasSpectatorPlane = 0; // no info at the moment

    histos.fill(HIST("hNEvents"), 9.5);
    histos.fill(HIST("hEventNchCorrelationAfterEP"), coll.multNTracksPVeta1(), coll.multNTracksGlobal());
    histos.fill(HIST("hEventPVcontributorsVsCentralityAfterEP"), coll.centFT0C(), coll.multNTracksPVeta1());
    histos.fill(HIST("hEventGlobalTracksVsCentralityAfterEP"), coll.centFT0C(), coll.multNTracksGlobal());
    histos.fill(HIST("hEventCentrality"), coll.centFT0C());
    histos.fill(HIST("hEventVertexZ"), coll.posZ());
    histos.fill(HIST("hMultNTracksITSTPCVsCentrality"), coll.centFT0C(), coll.multNTracksITSTPC());

    std::vector<float> bdtScore[nParticles];
    for (auto const& casc : Cascades) {

      if (!casc.has_cascMCCore())
        continue;

      auto cascMC = casc.cascMCCore_as<soa::Join<aod::CascMCCores, aod::CascMCCollRefs>>();

      int pdgCode{cascMC.pdgCode()};
      if (!(std::abs(pdgCode) == PDG_t::kXiMinus && std::abs(cascMC.pdgCodeV0()) == PDG_t::kLambda0 && std::abs(cascMC.pdgCodeBachelor()) == PDG_t::kPiPlus)       // Xi
          && !(std::abs(pdgCode) == PDG_t::kOmegaMinus && std::abs(cascMC.pdgCodeV0()) == PDG_t::kLambda0 && std::abs(cascMC.pdgCodeBachelor()) == PDG_t::kKPlus)) // Omega
      {
        pdgCode = 0;
      }

      // rapidity definition
      float XiY = RecoDecay::y(std::array{casc.px(), casc.py(), casc.pz()}, constants::physics::MassXiMinus);
      float OmegaY = RecoDecay::y(std::array{casc.px(), casc.py(), casc.pz()}, constants::physics::MassOmegaMinus);
      // true reco cascades before applying any selection
      if (std::abs(pdgCode) == PDG_t::kXiMinus && std::abs(cascMC.pdgCodeV0()) == PDG_t::kLambda0 && std::abs(cascMC.pdgCodeBachelor()) == PDG_t::kPiPlus) {
        histos.fill(HIST("hXiPtvsCent"), coll.centFT0C(), casc.pt());
        if (std::abs(casc.eta()) < 0.8)
          histos.fill(HIST("hXiPtvsCentEta08"), coll.centFT0C(), casc.pt());
        if (std::abs(XiY) < 0.5)
          histos.fill(HIST("hXiPtvsCentY05"), coll.centFT0C(), casc.pt());
      } else if (std::abs(pdgCode) == PDG_t::kOmegaMinus && std::abs(cascMC.pdgCodeV0()) == PDG_t::kLambda0 && std::abs(cascMC.pdgCodeBachelor()) == PDG_t::kKPlus) {
        histos.fill(HIST("hOmegaPtvsCent"), coll.centFT0C(), casc.pt());
        if (std::abs(casc.eta()) < 0.8)
          histos.fill(HIST("hOmegaPtvsCentEta08"), coll.centFT0C(), casc.pt());
        if (std::abs(OmegaY) < 0.5)
          histos.fill(HIST("hOmegaPtvsCentY05"), coll.centFT0C(), casc.pt());
      }

      /// Add some minimal cuts for single track variables (min number of TPC clusters)
      auto negExtra = casc.negTrackExtra_as<DauTracks>();
      auto posExtra = casc.posTrackExtra_as<DauTracks>();
      auto bachExtra = casc.bachTrackExtra_as<DauTracks>();

      int counter = 0;
      bool isCascCandidate = 0;
      isCascCandidate = IsCascAccepted(casc, negExtra, posExtra, bachExtra, counter);
      histos.fill(HIST("hCascade"), counter);
      histos.fill(HIST("hCascadeDauSel"), (int)isCascCandidate);
      if (!isCascCandidate)
        continue;

      // ML selections
      bool isSelectedCasc[nParticles]{false, false};

      std::vector<float> inputFeaturesCasc{casc.cascradius(),
                                           casc.v0radius(),
                                           casc.casccosPA(coll.posX(), coll.posY(), coll.posZ()),
                                           casc.v0cosPA(coll.posX(), coll.posY(), coll.posZ()),
                                           casc.dcapostopv(),
                                           casc.dcanegtopv(),
                                           casc.dcabachtopv(),
                                           casc.dcacascdaughters(),
                                           casc.dcaV0daughters(),
                                           casc.dcav0topv(coll.posX(), coll.posY(), coll.posZ()),
                                           casc.bachBaryonCosPA(),
                                           casc.bachBaryonDCAxyToPV()};

      float massCasc[nParticles]{casc.mXi(), casc.mOmega()};

      if (casc.pt() < CandidateConfigs.MinPt || casc.pt() > CandidateConfigs.MaxPt) {
        continue;
      }

      cascadev2::hMassBeforeSelVsPt[0]->Fill(massCasc[0], casc.pt());
      cascadev2::hMassBeforeSelVsPt[1]->Fill(massCasc[1], casc.pt());

      if (isApplyML) {
        // Retrieve model output and selection outcome
        isSelectedCasc[0] = mlResponseXi.isSelectedMl(inputFeaturesCasc, casc.pt(), bdtScore[0]);
        isSelectedCasc[1] = mlResponseOmega.isSelectedMl(inputFeaturesCasc, casc.pt(), bdtScore[1]);

        for (int iS{0}; iS < nParticles; ++iS) {
          // Fill BDT score histograms before selection
          cascadev2::hSignalScoreBeforeSel[iS]->Fill(bdtScore[0][1]);
          cascadev2::hBkgScoreBeforeSel[iS]->Fill(bdtScore[1][0]);

          // Fill histograms for selected candidates
          if (isSelectedCasc[iS]) {
            cascadev2::hSignalScoreAfterSel[iS]->Fill(bdtScore[0][1]);
            cascadev2::hBkgScoreAfterSel[iS]->Fill(bdtScore[1][0]);
            cascadev2::hMassAfterSelVsPt[iS]->Fill(massCasc[iS], casc.pt());
          }
        }
      } else {
        isSelectedCasc[0] = true;
        isSelectedCasc[1] = true;
      }

      histos.fill(HIST("hCascadePhi"), casc.phi());

      float BDTresponse[nParticles]{0.f, 0.f};
      const float psiT0C = 0; // not defined in MC for now
      auto v2CSP = 0;         // not defined in MC for now
      auto v2CEP = 0;         // not defined in MC for now
      auto v1SP_ZDCA = 0;     // not defined in MC for now
      auto v1SP_ZDCC = 0;     // not defined in MC for now

      if (isApplyML) {
        BDTresponse[0] = bdtScore[0][1];
        BDTresponse[1] = bdtScore[1][1];
      }
      if (isStoreTrueCascOnly) {
        if (pdgCode == 0)
          continue;
      }
      if (isSelectedCasc[0] || isSelectedCasc[1])
        fillAnalysedTable(coll, hasEventPlane, hasSpectatorPlane, casc, v2CSP, v2CEP, v1SP_ZDCA, v1SP_ZDCC, psiT0C, BDTresponse[0], BDTresponse[1], pdgCode);
    }
  }

  void processMCGen(MCCollisionsStra::iterator const& mcCollision, const soa::SmallGroups<soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraCollLabels>>& collisions, const soa::SmallGroups<soa::Join<aod::CascMCCores, aod::CascMCCollRefs>>& cascMC)
  {

    histosMCGen.fill(HIST("hZvertexGen"), mcCollision.posZ());
    histosMCGen.fill(HIST("hNEventsMC"), 0.5);
    // Generated with accepted z vertex
    if (std::abs(mcCollision.posZ()) > cutzvertex) {
      return;
    }
    histosMCGen.fill(HIST("hNEventsMC"), 1.5);
    // Check if there is at least one of the reconstructed collisions associated to this MC collision
    if (collisions.size() < 1)
      return;
    histosMCGen.fill(HIST("hNEventsMC"), 2.5);
    if (collisions.size() == 1)
      histosMCGen.fill(HIST("hNEventsMC"), 3.5);
    else if (collisions.size() == 2)
      histosMCGen.fill(HIST("hNEventsMC"), 4.5);

    int biggestNContribs = -1;
    float centrality = 100.5f;
    int nCollisions = 0;
    for (auto const& coll : collisions) {
      if (!AcceptEvent(coll, 0)) {
        continue;
      }
      if (biggestNContribs < coll.multPVTotalContributors()) {
        biggestNContribs = coll.multPVTotalContributors();
        centrality = coll.centFT0C();
      }
      nCollisions++;
    }
    if (nCollisions < 1) {
      return;
    }

    histosMCGen.fill(HIST("hNEventsMC"), 5.5);

    for (auto const& cascmc : cascMC) {
      if (std::abs(cascmc.pdgCode()) == PDG_t::kXiMinus)
        histosMCGen.fill(HIST("hNCascGen"), 0.5);
      else if (std::abs(cascmc.pdgCode()) == PDG_t::kOmegaMinus)
        histosMCGen.fill(HIST("hNCascGen"), 1.5);
      if (!cascmc.has_straMCCollision())
        continue;
      if (std::abs(cascmc.pdgCode()) == PDG_t::kXiMinus)
        histosMCGen.fill(HIST("hNCascGen"), 2.5);
      else if (std::abs(cascmc.pdgCode()) == PDG_t::kOmegaMinus)
        histosMCGen.fill(HIST("hNCascGen"), 3.5);
      if (!cascmc.isPhysicalPrimary())
        continue;
      if (std::abs(cascmc.pdgCode()) == PDG_t::kXiMinus)
        histosMCGen.fill(HIST("hNCascGen"), 4.5);
      else if (std::abs(cascmc.pdgCode()) == PDG_t::kOmegaMinus)
        histosMCGen.fill(HIST("hNCascGen"), 5.5);

      float ptmc = RecoDecay::sqrtSumOfSquares(cascmc.pxMC(), cascmc.pyMC());

      float theta = std::atan(ptmc / cascmc.pzMC()); //-pi/2 < theta < pi/2

      float theta1 = 0;

      // if pz is positive (i.e. positive rapidity): 0 < theta < pi/2
      if (theta > 0)
        theta1 = theta; // 0 < theta1/2 < pi/4 --> 0 < tan (theta1/2) < 1 --> positive eta
      // if pz is negative (i.e. negative rapidity): -pi/2 < theta < 0 --> we need 0 < theta1/2 < pi/2 for the ln to be defined
      else
        theta1 = o2::constants::math::PI + theta; // pi/2 < theta1 < pi --> pi/4 < theta1/2 <  pi/2 --> 1 < tan (theta1/2) --> negative eta

      float cascMCeta = -std::log(std::tan(theta1 / 2));
      float cascMCy = 0;
      if (std::abs(cascmc.pdgCode()) == PDG_t::kXiMinus) {
        cascMCy = RecoDecay::y(std::array{cascmc.pxMC(), cascmc.pyMC(), cascmc.pzMC()}, constants::physics::MassXiMinus);
        if (std::abs(cascMCeta) < etaCascMCGen) {
          histosMCGen.fill(HIST("h2DGenXiEta08"), centrality, ptmc);
          histosMCGen.fill(HIST("hNCascGen"), 6.5);
        }
        if (std::abs(cascMCy) < yCascMCGen)
          histosMCGen.fill(HIST("h2DGenXiY05"), centrality, ptmc);
        histosMCGen.fill(HIST("hGenXiY"), cascMCy);
      } else if (std::abs(cascmc.pdgCode()) == PDG_t::kOmegaMinus) {
        cascMCy = RecoDecay::y(std::array{cascmc.pxMC(), cascmc.pyMC(), cascmc.pzMC()}, constants::physics::MassOmegaMinus);
        if (std::abs(cascMCeta) < etaCascMCGen) {
          histosMCGen.fill(HIST("h2DGenOmegaEta08"), centrality, ptmc);
          histosMCGen.fill(HIST("hNCascGen"), 7.5);
        }
        if (std::abs(cascMCy) < yCascMCGen)
          histosMCGen.fill(HIST("h2DGenOmegaY05"), centrality, ptmc);
        histosMCGen.fill(HIST("hGenOmegaY"), cascMCy);
      }
    }
  }

  void processMCPrimaryLambdaFraction(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels>::iterator const& coll, V0MCCandidates const& V0s, DauTracks const&, soa::Join<aod::V0MCCores, aod::V0MCCollRefs> const&)
  {

    Float_t collisionCentrality = 0;
    if (isCollisionCentrality == 0) { // T0C
      collisionCentrality = coll.centFT0C();
    } else if (isCollisionCentrality == 1) { // T0M
      collisionCentrality = coll.centFT0M();
    }

    histos.fill(HIST("hEventCentralityBefEvSel"), collisionCentrality);

    if (!AcceptEvent(coll, 1)) {
      return;
    }

    // no EP requirements as in MC we do not have EP info
    histos.fill(HIST("hNEvents"), 9.5);
    histos.fill(HIST("hEventCentrality"), collisionCentrality);
    histos.fill(HIST("hEventVertexZ"), coll.posZ());

    for (auto const& v0 : V0s) {

      /// Add some minimal cuts for single track variables (min number of TPC clusters)
      auto negExtra = v0.negTrackExtra_as<DauTracks>();
      auto posExtra = v0.posTrackExtra_as<DauTracks>();

      int counterLambda = 0;
      int counterALambda = 0;
      bool isLambdaCandidate = 0;
      bool isALambdaCandidate = 0;

      // check if v0 has MC info
      if (!v0.has_v0MCCore())
        continue;

      //-- selections ------------------------------------------------------------
      if (isLambdaAccepted(negExtra, posExtra, counterLambda))
        isLambdaCandidate = 1;
      if (isAntiLambdaAccepted(negExtra, posExtra, counterALambda))
        isALambdaCandidate = 1;
      histos.fill(HIST("hLambdaDauSel"), counterLambda);
      histos.fill(HIST("hALambdaDauSel"), counterALambda);

      if (v0.pt() < V0Configs.MinPtV0 || v0.pt() > V0Configs.MaxPtV0) {
        continue;
      }

      float massV0[nCharges]{v0.mLambda(), v0.mAntiLambda()};
      lambdav2::hMassBeforeSelVsPt[0]->Fill(massV0[0], v0.pt());
      lambdav2::hMassBeforeSelVsPt[1]->Fill(massV0[1], v0.pt());

      bool isSelectedV0[2]{false, false};
      if (isV0TopoAccepted(v0) && isLambdaCandidate)
        isSelectedV0[0] = true;
      if (isV0TopoAccepted(v0) && isALambdaCandidate)
        isSelectedV0[1] = true;

      int chargeIndex = -1;
      if (isSelectedV0[0] && !isSelectedV0[1]) { // Lambdas
        histos.fill(HIST("hLambdaCandidate"), 0);
        chargeIndex = 0;
      }
      if (isSelectedV0[1] && !isSelectedV0[0]) { // AntiLambdas
        histos.fill(HIST("hLambdaCandidate"), 1);
        chargeIndex = 1;
      }
      if (isSelectedV0[0] && isSelectedV0[1]) {
        histos.fill(HIST("hLambdaCandidate"), 2);
        if (v0.mLambda() > V0Configs.MinMassLambda && v0.mLambda() < V0Configs.MaxMassLambda && v0.mAntiLambda() > V0Configs.MinMassLambda && v0.mAntiLambda() < V0Configs.MaxMassLambda) {
          histos.fill(HIST("hLambdaCandidate"), 3);
          continue; // in case of ambiguity between Lambda and AntiLambda, I skip the particle; checked to be zero in range 1.105 - 1.125
        }
        if (v0.mLambda() > V0Configs.MinMassLambda && v0.mLambda() < V0Configs.MaxMassLambda)
          chargeIndex = 0;
        else if (v0.mAntiLambda() > V0Configs.MinMassLambda && v0.mAntiLambda() < V0Configs.MaxMassLambda)
          chargeIndex = 1;
        else {
          chargeIndex = 2; // these are bkg candidates
          histos.fill(HIST("hLambdaCandidate"), 4);
        }
      }
      if (!isSelectedV0[0] && !isSelectedV0[1])
        continue;

      lambdav2::hMassAfterSelVsPt[0]->Fill(massV0[0], v0.pt());
      lambdav2::hMassAfterSelVsPt[1]->Fill(massV0[1], v0.pt());
      //--------------------------------------------------------------

      auto v0MC = v0.v0MCCore_as<soa::Join<aod::V0MCCores, aod::V0MCCollRefs>>();
      int pdgCode{v0MC.pdgCode()};
      // select true lambdas
      bool isTrueLambda = 0;
      bool isTrueALambda = 0;
      if (pdgCode == PDG_t::kLambda0 && v0MC.pdgCodePositive() == PDG_t::kProton && v0MC.pdgCodeNegative() == PDG_t::kPiMinus)
        isTrueLambda = 1;
      else if (pdgCode == PDG_t::kLambda0Bar && v0MC.pdgCodePositive() == PDG_t::kPiPlus && v0MC.pdgCodeNegative() == PDG_t::kProtonBar)
        isTrueALambda = 1;
      if (!isTrueLambda && !isTrueALambda)
        continue;
      if (isTrueLambda)
        lambdav2::hMassAfterSelVsPtTrue[0]->Fill(massV0[0], v0.pt());
      if (isTrueALambda)
        lambdav2::hMassAfterSelVsPtTrue[1]->Fill(massV0[1], v0.pt());

      bool isPrimary = v0MC.isPhysicalPrimary();

      // histo for primary fraction
      if (isTrueLambda) {
        if (isPrimary) {
          histos.fill(HIST("hCentvsPtvsPrimaryFracLambda"), collisionCentrality, v0.pt(), 0);
          histos.fill(HIST("hCentvsPrimaryFracLambda"), collisionCentrality, 0);
        } else {
          histos.fill(HIST("hCentvsPtvsPrimaryFracLambda"), collisionCentrality, v0.pt(), 1);
          histos.fill(HIST("hCentvsPrimaryFracLambda"), collisionCentrality, 1);
        }
      } else if (isTrueALambda) {
        if (isPrimary) {
          histos.fill(HIST("hCentvsPtvsPrimaryFracLambda"), collisionCentrality, v0.pt(), 2);
          histos.fill(HIST("hCentvsPrimaryFracLambda"), collisionCentrality, 2);
        } else {
          histos.fill(HIST("hCentvsPtvsPrimaryFracLambda"), collisionCentrality, v0.pt(), 3);
          histos.fill(HIST("hCentvsPrimaryFracLambda"), collisionCentrality, 3);
        }
      }
    }
  }

  PROCESS_SWITCH(cascadeFlow, processTrainingBackground, "Process to create the training dataset for the background", true);
  PROCESS_SWITCH(cascadeFlow, processTrainingSignal, "Process to create the training dataset for the signal", false);
  PROCESS_SWITCH(cascadeFlow, processAnalyseData, "Process to apply ML model to the data", false);
  PROCESS_SWITCH(cascadeFlow, processAnalyseDataEP2CentralFW, "Process to apply ML model to the data - second order event plane calibration from central framework", false);
  PROCESS_SWITCH(cascadeFlow, processAnalyseDataEPCentralFW, "Process to apply ML model to the data - event plane calibration from central framework", false);
  PROCESS_SWITCH(cascadeFlow, processAnalyseLambdaEP2CentralFW, "Process to measure flow and polarization of Lambda - event plane calibration from central framework", false);
  PROCESS_SWITCH(cascadeFlow, processAnalyseMC, "Process to apply ML model to the MC", false);
  PROCESS_SWITCH(cascadeFlow, processMCGen, "Process to store MC generated particles", false);
  PROCESS_SWITCH(cascadeFlow, processMCPrimaryLambdaFraction, "Process to compute primary lambda fraction", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<cascadeFlow>(cfgc, TaskName{"lf-cascade-flow"})};
}
