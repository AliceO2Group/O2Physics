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
/// \brief Task to create derived data for cascade flow analyses
/// \authors: Chiara De Martin (chiara.de.martin@cern.ch), Maximiliano Puccio (maximiliano.puccio@cern.ch)

#include <vector>
#include <string>
#include <memory>
#include "Math/Vector3D.h"
#include "TRandom3.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/cascqaanalysis.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "Tools/ML/MlResponse.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::analysis;
using namespace o2::framework::expressions;
using std::array;

using DauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;
using CollEventPlane = soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraFT0CQVs, aod::StraFT0CQVsEv, aod::StraTPCQVs>::iterator;
using CollEventAndSpecPlane = soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraFT0CQVs, aod::StraFT0CQVsEv, aod::StraTPCQVs, aod::StraZDCSP>::iterator;
using CollEventPlaneCentralFW = soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraFT0CQVs, aod::StraTPCQVs, aod::StraZDCSP>::iterator;
using MCCollisionsStra = soa::Join<aod::StraMCCollisions, aod::StraMCCollMults>;
using CascCandidates = soa::Join<aod::CascCollRefs, aod::CascCores, aod::CascExtras, aod::CascBBs>;
using CascMCCandidates = soa::Join<aod::CascCollRefs, aod::CascCores, aod::CascExtras, aod::CascBBs, aod::CascCoreMCLabels>;

namespace cascadev2
{
enum species { Xi = 0,
               Omega = 1 };
constexpr double massSigmaParameters[4][2]{
  {4.9736e-3, 0.006815},
  {-2.39594, -2.257},
  {1.8064e-3, 0.00138},
  {1.03468e-1, 0.1898}};
static const std::vector<std::string> massSigmaParameterNames{"p0", "p1", "p2", "p3"};
static const std::vector<std::string> speciesNames{"Xi", "Omega"};
const double AlphaXi[2] = {-0.390, 0.371};     // decay parameter of XiMinus and XiPlus
const double AlphaOmega[2] = {0.0154, -0.018}; // decay parameter of OmegaMinus and OmegaPlus
const double AlphaLambda[2] = {0.747, -0.757}; // decay parameter of Lambda and AntiLambda

std::shared_ptr<TH2> hMassBeforeSelVsPt[2];
std::shared_ptr<TH2> hMassAfterSelVsPt[2];
std::shared_ptr<TH1> hSignalScoreBeforeSel[2];
std::shared_ptr<TH1> hBkgScoreBeforeSel[2];
std::shared_ptr<TH1> hSignalScoreAfterSel[2];
std::shared_ptr<TH1> hBkgScoreAfterSel[2];
std::shared_ptr<THn> hSparseV2C[2];
} // namespace cascadev2

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

  // Output filling criteria
  Configurable<bool> isFillTree{"isFillTree", 1, ""};
  Configurable<bool> isFillTHNXi{"isFillTHNXi", 1, ""};
  Configurable<bool> isFillTHNXi_PzVsPsi{"isFillTHNXi_PzVsPsi", 1, ""};
  Configurable<bool> isFillTHNOmega{"isFillTHNOmega", 1, ""};
  Configurable<bool> isFillTHNOmega_PzVsPsi{"isFillTHNOmega_PzVsPsi", 1, ""};

  // axes
  ConfigurableAxis axisQVs{"axisQVs", {500, -10.f, 10.f}, "axisQVs"};
  ConfigurableAxis axisQVsNorm{"axisQVsNorm", {200, -1.f, 1.f}, "axisQVsNorm"};

  // THN axes
  ConfigurableAxis thnConfigAxisFT0C{"thnConfigAxisFT0C", {8, 0, 80}, "FT0C centrality (%)"};
  ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {VARIABLE_WIDTH, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2, 2.25, 2.5, 2.75, 3, 3.5, 4, 5, 6, 8, 10}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis thnConfigAxisCharge{"thnConfigAxisCharge", {2, 0, 2}, ""};
  ConfigurableAxis thnConfigAxisPsiDiff{"thnConfigAxisPsiDiff", {100, 0, 2 * TMath::Pi()}, ""};
  ConfigurableAxis thnConfigAxisMassXi{"thnConfigAxisMassXi", {45, 1.300, 1.345}, ""};
  ConfigurableAxis thnConfigAxisMassOmega{"thnConfigAxisMassOmega", {45, 1.655, 1.690}, ""};
  ConfigurableAxis thnConfigAxisMassLambda{"thnConfigAxisMassLambda", {60, 1.1, 1.13}, ""};
  ConfigurableAxis thnConfigAxisBDTScore{"thnConfigAxisBDTScore", {15, 0.4, 1}, ""};
  ConfigurableAxis thnConfigAxisV2{"thnConfigAxiV2", {100, -1., 1.}, ""};
  ConfigurableAxis thnConfigAxisPzs2Xi{"thnConfigAxiPzs2Xi", {200, -2.8, 2.8}, ""};
  ConfigurableAxis thnConfigAxisPzs2Omega{"thnConfigAxiPzs2Omega", {200, -70, 70}, ""};
  ConfigurableAxis thnConfigAxisPzs2Lambda{"thnConfigAxiPzs2Lambda", {200, -2, 2}, ""};
  ConfigurableAxis thnConfigAxisCos2Theta{"thnConfigAxiCos2Theta", {100, 0, 1}, ""};
  ConfigurableAxis thnConfigAxisCosThetaXiAlpha{"thnConfigAxisCosThetaXiAlpha", {200, -2.8, 2.8}, ""};
  ConfigurableAxis thnConfigAxisCosThetaOmegaAlpha{"thnConfigAxisCosThetaOmegaAlpha", {200, -70, 70}, ""};
  ConfigurableAxis thnConfigAxisCosThetaProtonAlpha{"thnConfigAxisCosThetaProtonAlpha", {200, -2, 2}, ""};

  // Event selection criteria
  Configurable<bool> isStoreTrueCascOnly{"isStoreTrueCascOnly", 1, ""};
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

  Configurable<float> MinPt{"MinPt", 0.6, "Min pt of cascade"};
  Configurable<float> MaxPt{"MaxPt", 10, "Max pt of cascade"};
  Configurable<double> sideBandStart{"sideBandStart", 5, "Start of the sideband region in number of sigmas"};
  Configurable<double> sideBandEnd{"sideBandEnd", 7, "End of the sideband region in number of sigmas"};
  Configurable<double> downsample{"downsample", 1., "Downsample training output tree"};
  Configurable<bool> doNTPCSigmaCut{"doNTPCSigmaCut", 1, "doNtpcSigmaCut"};
  Configurable<float> nsigmatpcPr{"nsigmatpcPr", 5, "nsigmatpcPr"};
  Configurable<float> nsigmatpcPi{"nsigmatpcPi", 5, "nsigmatpcPi"};
  Configurable<float> mintpccrrows{"mintpccrrows", 70, "mintpccrrows"};
  Configurable<float> etaCascMCGen{"etaCascMCGen", 0.8, "etaCascMCGen"};
  Configurable<float> etaCasc{"etaCasc", 0.8, "etaCasc"};
  Configurable<float> yCascMCGen{"yCascMCGen", 0.5, "yCascMCGen"};

  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::vector<std::string>> modelPathsCCDBXi{"modelPathsCCDBXi", std::vector<std::string>{"Users/c/chdemart/CascadesFlow"}, "Paths of models on CCDB"};
  Configurable<std::vector<std::string>> modelPathsCCDBOmega{"modelPathsCCDBOmega", std::vector<std::string>{"Users/c/chdemart/CascadesFlow"}, "Paths of models on CCDB"};
  Configurable<std::vector<std::string>> onnxFileNamesXi{"onnxFileNamesXi", std::vector<std::string>{"model_onnx.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
  Configurable<std::vector<std::string>> onnxFileNamesOmega{"onnxFileNamesOmega", std::vector<std::string>{"model_onnx.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", true, "Flag to enable or disable the loading of models from CCDB"};

  // ML inference
  Configurable<bool> isApplyML{"isApplyML", 1, "Flag to apply ML selections"};
  Configurable<std::vector<double>> binsPtMl{"binsPtMl", std::vector<double>{cascade_flow_cuts_ml::vecBinsPt}, "pT bin limits for ML application"};
  Configurable<std::vector<int>> cutDirMl{"cutDirMl", std::vector<int>{cascade_flow_cuts_ml::vecCutDir}, "Whether to reject score values greater or smaller than the threshold"};
  Configurable<LabeledArray<double>> cutsMl{"cutsMl", {cascade_flow_cuts_ml::cuts[0], cascade_flow_cuts_ml::nBinsPt, cascade_flow_cuts_ml::nCutScores, cascade_flow_cuts_ml::labelsPt, cascade_flow_cuts_ml::labelsCutScore}, "ML selections per pT bin"};
  Configurable<int> nClassesMl{"nClassesMl", static_cast<int>(cascade_flow_cuts_ml::nCutScores), "Number of classes in ML model"};

  o2::ccdb::CcdbApi ccdbApi;

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
    if (TMath::Abs(collision.posZ()) > cutzvertex) {
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

  double GetPhiInRange(double phi)
  {
    while (phi < 0) {
      phi += TMath::Pi();
    }
    while (phi > TMath::Pi()) {
      phi -= TMath::Pi();
    }
    return phi;
  }

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry histosMCGen{"histosMCGen", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry resolution{"resolution", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};

  // Tables to produce
  Produces<aod::CascTraining> trainingSample;
  Produces<aod::CascAnalysis> analysisSample;
  Configurable<LabeledArray<double>> parSigmaMass{
    "parSigmaMass",
    {cascadev2::massSigmaParameters[0], 4, 2,
     cascadev2::massSigmaParameterNames, cascadev2::speciesNames},
    "Mass resolution parameters: [0]*exp([1]*x)+[2]*exp([3]*x)"};

  float getNsigmaMass(const cascadev2::species s, const float pt, const float nsigma = 6)
  {
    const auto sigma = parSigmaMass->get(0u, s) * exp(parSigmaMass->get(1, s) * pt) + parSigmaMass->get(2, s) * exp(parSigmaMass->get(3, s) * pt);
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
    for (int i{0}; i < 2; ++i) {
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

  void init(InitContext const&)
  {

    float minMass[2]{1.28, 1.6};
    float maxMass[2]{1.36, 1.73};
    const AxisSpec massCascAxis[2]{{static_cast<int>((maxMass[0] - minMass[0]) / 0.001f), minMass[0], maxMass[0], "#Xi candidate mass (GeV/c^{2})"},
                                   {static_cast<int>((maxMass[1] - minMass[1]) / 0.001f), minMass[1], maxMass[1], "#Omega candidate mass (GeV/c^{2})"}};
    const AxisSpec ptAxis{static_cast<int>((MaxPt - MinPt) / 0.2), MinPt, MaxPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec v2Axis{200, -1., 1., "#it{v}_{2}"};
    const AxisSpec CentAxis{18, 0., 90., "FT0C centrality percentile"};
    TString hNEventsLabels[10] = {"All", "sel8", "z vrtx", "kNoSameBunchPileup", "kIsGoodZvtxFT0vsPV", "trackOccupancyInTimeRange", "kNoCollInTimeRange", "kNoCollInROF", "kTVXinTRD", "kIsGoodEventEP"};
    TString hNEventsLabelsMC[6] = {"All", "z vtx", ">=1RecoColl", "1Reco", "2Reco", "EvSelected"};
    TString hNCascLabelsMC[8] = {"All Xi", "all Omega", "Xi: has MC coll", "Om: has MC coll", "Xi: isPrimary", "Om: is Primary", "Xi: |eta|<0.8", "Om: |eta| < 0.8"};

    resolution.add("QVectorsT0CTPCA", "QVectorsT0CTPCA", HistType::kTH2F, {axisQVs, CentAxis});
    resolution.add("QVectorsT0CTPCC", "QVectorsT0CTPCC", HistType::kTH2F, {axisQVs, CentAxis});
    resolution.add("QVectorsTPCAC", "QVectorsTPCAC", HistType::kTH2F, {axisQVs, CentAxis});
    resolution.add("QVectorsNormT0CTPCA", "QVectorsNormT0CTPCA", HistType::kTH2F, {axisQVsNorm, CentAxis});
    resolution.add("QVectorsNormT0CTPCC", "QVectorsNormT0CTPCC", HistType::kTH2F, {axisQVsNorm, CentAxis});
    resolution.add("QVectorsNormTPCAC", "QVectorsNormTPCCB", HistType::kTH2F, {axisQVsNorm, CentAxis});
    resolution.add("QVectorsSpecPlane", "QVectorsSpecPlane", HistType::kTH2F, {axisQVsNorm, CentAxis});

    histos.add("hNEvents", "hNEvents", {HistType::kTH1F, {{10, 0.f, 10.f}}});
    for (Int_t n = 1; n <= histos.get<TH1>(HIST("hNEvents"))->GetNbinsX(); n++) {
      histos.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(n, hNEventsLabels[n - 1]);
    }
    histos.add("hEventVertexZ", "hEventVertexZ", kTH1F, {{120, -12., 12.}});
    histos.add("hEventCentrality", "hEventCentrality", kTH1F, {{101, 0, 101}});
    histos.add("hPsiT0C", "hPsiT0C", HistType::kTH1D, {{100, -TMath::Pi(), TMath::Pi()}});
    histos.add("hPsiT0CvsCentFT0C", "hPsiT0CvsCentFT0C", HistType::kTH2D, {CentAxis, {100, -TMath::Pi(), TMath::Pi()}});
    histos.add("hPsiZDCA_vs_ZDCC", "hPsiZDCA_vs_ZDCC", HistType::kTH2D, {{100, -TMath::Pi(), TMath::Pi()}, {100, -TMath::Pi(), TMath::Pi()}});
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
    histos.add("hCascadeSignal", "hCascadeSignal", HistType::kTH1F, {{6, -0.5, 5.5}});
    histos.add("hCascade", "hCascade", HistType::kTH1F, {{6, -0.5, 5.5}});
    histos.add("hXiPtvsCent", "hXiPtvsCent", HistType::kTH2F, {{100, 0, 100}, {400, 0, 20}});
    histos.add("hXiPtvsCentEta08", "hXiPtvsCentEta08", HistType::kTH2F, {{100, 0, 100}, {400, 0, 20}});
    histos.add("hXiPtvsCentY05", "hXiPtvsCentY05", HistType::kTH2F, {{100, 0, 100}, {400, 0, 20}});
    histos.add("hOmegaPtvsCent", "hOmegaPtvsCent", HistType::kTH2F, {{100, 0, 100}, {400, 0, 20}});
    histos.add("hOmegaPtvsCentEta08", "hOmegaPtvsCentEta08", HistType::kTH2F, {{100, 0, 100}, {400, 0, 20}});
    histos.add("hOmegaPtvsCentY05", "hOmegaPtvsCentY05", HistType::kTH2F, {{100, 0, 100}, {400, 0, 20}});
    histos.add("hCascadePhi", "hCascadePhi", HistType::kTH1F, {{100, 0, 2 * TMath::Pi()}});
    histos.add("hcascminuspsiT0C", "hcascminuspsiT0C", HistType::kTH1F, {{100, 0, TMath::Pi()}});
    histos.add("hv2CEPvsFT0C", "hv2CEPvsFT0C", HistType::kTH2F, {CentAxis, {100, -1, 1}});
    histos.add("hv2CEPvsv2CSP", "hv2CEPvsV2CSP", HistType::kTH2F, {{100, -1, 1}, {100, -1, 1}});
    histos.add("hv1EPvsv1SP", "hV1EPvsV1SP", HistType::kTH2F, {{100, -1, 1}, {100, -1, 1}});
    histos.add("hv1SP_ZDCA_vs_ZDCC", "hv1SP_ZDCA_vs_ZDCC", HistType::kTH2F, {{100, -1, 1}, {100, -1, 1}});

    const AxisSpec thnAxisFT0C{thnConfigAxisFT0C, "FT0C (%)"};
    const AxisSpec thnAxisPt{thnConfigAxisPt, "p_{T}"};
    const AxisSpec thnAxisCharge{thnConfigAxisCharge, "Charge"};
    const AxisSpec thnAxisPsiDiff{thnConfigAxisPsiDiff, "2(phi-Psi)"};
    const AxisSpec thnAxisMassXi{thnConfigAxisMassXi, "inv. mass (#Lambda #pi) (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisMassOmega{thnConfigAxisMassOmega, "inv. mass (#Lambda K) (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisMassLambda{thnConfigAxisMassLambda, "inv. mass (p #pi) (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisBDTScore{thnConfigAxisBDTScore, "BDT score"};
    const AxisSpec thnAxisV2{thnConfigAxisV2, "v_{2}"};
    const AxisSpec thnAxisPzs2Xi{thnConfigAxisPzs2Xi, "Pzs2Xi"};
    const AxisSpec thnAxisPzs2Omega{thnConfigAxisPzs2Omega, "Pzs2Omega"};
    const AxisSpec thnAxisPzs2Lambda{thnConfigAxisPzs2Lambda, "Pzs2Lambda"};
    const AxisSpec thnAxisCos2Theta{thnConfigAxisCos2Theta, "Cos2Theta"};
    const AxisSpec thnAxisCosThetaXiAlpha{thnConfigAxisCosThetaXiAlpha, "CosThetaXiWithAlpha"};
    const AxisSpec thnAxisCosThetaOmegaAlpha{thnConfigAxisCosThetaOmegaAlpha, "CosThetaOmegaWithAlpha"};
    const AxisSpec thnAxisCosThetaProtonAlpha{thnConfigAxisCosThetaProtonAlpha, "CosThetaProtonWithAlpha"};

    if (isFillTHNXi) {
      histos.add("hXiV2", "THn for v2 of Xi", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisPt, thnAxisMassXi, thnAxisBDTScore, thnAxisV2});
      histos.add("hXiPzs2", "THn for Pzs2 of Xi", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisPt, thnAxisMassXi, thnAxisBDTScore, thnAxisPzs2Xi});
      histos.add("hXiPzs2FromLambda", "THn for Pzs2 of Xi", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisPt, thnAxisMassXi, thnAxisBDTScore, thnAxisPzs2Lambda});
      histos.add("hXiCos2Theta", "THn for Cos2Theta of Xi", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisPt, thnAxisMassXi, thnAxisBDTScore, thnAxisCos2Theta});
      histos.add("hXiCos2ThetaFromLambda", "THn for Cos2Theta of Xi", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisPt, thnAxisMassXi, thnAxisBDTScore, thnAxisCos2Theta});
    }
    if (isFillTHNXi_PzVsPsi) {
      histos.add("hXiPzVsPsi", "THn for cosTheta of Xi", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisPt, thnAxisMassXi, thnAxisBDTScore, thnAxisCosThetaXiAlpha, thnAxisPsiDiff});
      histos.add("hXiPzVsPsiFromLambda", "THn for cosTheta of Xi", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisPt, thnAxisMassXi, thnAxisBDTScore, thnAxisCosThetaProtonAlpha, thnAxisPsiDiff});
      histos.add("hXiCos2ThetaVsPsi", "THn for cos2Theta of Xi", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisPt, thnAxisMassXi, thnAxisBDTScore, thnAxisCos2Theta, thnAxisPsiDiff});
      histos.add("hXiCos2ThetaVsPsiFromLambda", "THn for cos2Theta of Xi", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisPt, thnAxisMassXi, thnAxisBDTScore, thnAxisCos2Theta, thnAxisPsiDiff});
    }
    if (isFillTHNOmega) {
      histos.add("hOmegaV2", "THn for v2 of Omega", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisPt, thnAxisMassOmega, thnAxisBDTScore, thnAxisV2});
      histos.add("hOmegaPzs2", "THn for Pzs2 of Omega", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisPt, thnAxisMassOmega, thnAxisBDTScore, thnAxisPzs2Omega});
      histos.add("hOmegaPzs2FromLambda", "THn for Pzs2 of Omega", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisPt, thnAxisMassOmega, thnAxisBDTScore, thnAxisPzs2Lambda});
      histos.add("hOmegaCos2Theta", "THn for Cos2Theta of Omega", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisPt, thnAxisMassOmega, thnAxisBDTScore, thnAxisCos2Theta});
      histos.add("hOmegaCos2ThetaFromLambda", "THn for Cos2Theta of Omega", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisPt, thnAxisMassOmega, thnAxisBDTScore, thnAxisCos2Theta});
    }
    if (isFillTHNOmega_PzVsPsi) {
      histos.add("hOmegaPzVsPsi", "THn for cosTheta of Omega", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisPt, thnAxisMassOmega, thnAxisBDTScore, thnAxisCosThetaOmegaAlpha, thnAxisPsiDiff});
      histos.add("hOmegaPzVsPsiFromLambda", "THn for cosTheta of Omega", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisPt, thnAxisMassOmega, thnAxisBDTScore, thnAxisCosThetaProtonAlpha, thnAxisPsiDiff});
      histos.add("hOmegaCos2ThetaVsPsi", "THn for cos2Theta of Omega", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisPt, thnAxisMassOmega, thnAxisBDTScore, thnAxisCos2Theta, thnAxisPsiDiff});
      histos.add("hOmegaCos2ThetaVsPsiFromLambda", "THn for cos2Theta of Omega", HistType::kTHnF, {thnAxisFT0C, thnAxisCharge, thnAxisPt, thnAxisMassOmega, thnAxisBDTScore, thnAxisCos2Theta, thnAxisPsiDiff});
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

    for (int iS{0}; iS < 2; ++iS) {
      cascadev2::hMassBeforeSelVsPt[iS] = histos.add<TH2>(Form("hMassBeforeSelVsPt%s", cascadev2::speciesNames[iS].data()), "hMassBeforeSelVsPt", HistType::kTH2F, {massCascAxis[iS], ptAxis});
      cascadev2::hMassAfterSelVsPt[iS] = histos.add<TH2>(Form("hMassAfterSelVsPt%s", cascadev2::speciesNames[iS].data()), "hMassAfterSelVsPt", HistType::kTH2F, {massCascAxis[iS], ptAxis});
      cascadev2::hSignalScoreBeforeSel[iS] = histos.add<TH1>(Form("hSignalScoreBeforeSel%s", cascadev2::speciesNames[iS].data()), "Signal score before selection;BDT first score;entries", HistType::kTH1F, {{100, 0., 1.}});
      cascadev2::hBkgScoreBeforeSel[iS] = histos.add<TH1>(Form("hBkgScoreBeforeSel%s", cascadev2::speciesNames[iS].data()), "Bkg score before selection;BDT first score;entries", HistType::kTH1F, {{100, 0., 1.}});
      cascadev2::hSignalScoreAfterSel[iS] = histos.add<TH1>(Form("hSignalScoreAfterSel%s", cascadev2::speciesNames[iS].data()), "Signal score after selection;BDT first score;entries", HistType::kTH1F, {{100, 0., 1.}});
      cascadev2::hBkgScoreAfterSel[iS] = histos.add<TH1>(Form("hBkgScoreAfterSel%s", cascadev2::speciesNames[iS].data()), "Bkg score after selection;BDT first score;entries", HistType::kTH1F, {{100, 0., 1.}});
      cascadev2::hSparseV2C[iS] = histos.add<THn>(Form("hSparseV2C%s", cascadev2::speciesNames[iS].data()), "hSparseV2C", HistType::kTHnF, {massCascAxis[iS], ptAxis, v2Axis, CentAxis});
    }
    if (isApplyML) {
      // Configure and initialise the ML class
      mlResponseXi.configure(binsPtMl, cutsMl, cutDirMl, nClassesMl);
      mlResponseOmega.configure(binsPtMl, cutsMl, cutDirMl, nClassesMl);
      // Bonus: retrieve the model from CCDB (needed for ML application on the GRID)
      if (loadModelsFromCCDB) {
        ccdbApi.init(ccdbUrl);
        mlResponseXi.setModelPathsCCDB(onnxFileNamesXi, ccdbApi, modelPathsCCDBXi, timestampCCDB);
        mlResponseOmega.setModelPathsCCDB(onnxFileNamesOmega, ccdbApi, modelPathsCCDBOmega, timestampCCDB); // TODO: use different model for Xi and Omega
      } else {
        mlResponseXi.setModelPathsLocal(onnxFileNamesXi);
        mlResponseOmega.setModelPathsLocal(onnxFileNamesOmega);
      }
      mlResponseXi.init();
      mlResponseOmega.init();
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

    for (auto& casc : Cascades) {
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

    for (auto& casc : Cascades) {
      if (!casc.has_cascMCCore())
        continue;

      auto cascMC = casc.cascMCCore_as<soa::Join<aod::CascMCCores, aod::CascMCCollRefs>>();
      int pdgCode{cascMC.pdgCode()};
      if (!(std::abs(pdgCode) == 3312 && std::abs(cascMC.pdgCodeV0()) == 3122 && std::abs(cascMC.pdgCodeBachelor()) == 211)     // Xi
          && !(std::abs(pdgCode) == 3334 && std::abs(cascMC.pdgCodeV0()) == 3122 && std::abs(cascMC.pdgCodeBachelor()) == 321)) // Omega
        continue;

      auto negExtra = casc.negTrackExtra_as<DauTracks>();
      auto posExtra = casc.posTrackExtra_as<DauTracks>();
      auto bachExtra = casc.bachTrackExtra_as<DauTracks>();

      int counter = 0;
      IsCascAccepted(casc, negExtra, posExtra, bachExtra, counter);
      histos.fill(HIST("hCascadeSignal"), counter);

      // PDG cascades
      fillTrainingTable(coll, casc, pdgCode);
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

    const float PsiT0C = std::atan2(coll.qvecFT0CIm(), coll.qvecFT0CRe()) * 0.5f;
    histos.fill(HIST("hPsiT0C"), PsiT0C);
    histos.fill(HIST("hPsiZDCA_vs_ZDCC"), coll.psiZDCC(), coll.psiZDCA());
    histos.fill(HIST("hPsiT0CvsCentFT0C"), coll.centFT0C(), PsiT0C);

    resolution.fill(HIST("QVectorsT0CTPCA"), eventplaneVecT0C.Dot(eventplaneVecTPCA), coll.centFT0C());
    resolution.fill(HIST("QVectorsT0CTPCC"), eventplaneVecT0C.Dot(eventplaneVecTPCC), coll.centFT0C());
    resolution.fill(HIST("QVectorsTPCAC"), eventplaneVecTPCA.Dot(eventplaneVecTPCC), coll.centFT0C());
    resolution.fill(HIST("QVectorsNormT0CTPCA"), eventplaneVecT0C.Dot(eventplaneVecTPCA) / (coll.qTPCR() * coll.sumAmplFT0C()), coll.centFT0C());
    resolution.fill(HIST("QVectorsNormT0CTPCC"), eventplaneVecT0C.Dot(eventplaneVecTPCC) / (coll.qTPCL() * coll.sumAmplFT0C()), coll.centFT0C());
    resolution.fill(HIST("QVectorsNormTPCAC"), eventplaneVecTPCA.Dot(eventplaneVecTPCC) / (coll.qTPCR() * coll.qTPCL()), coll.centFT0C());
    resolution.fill(HIST("QVectorsSpecPlane"), spectatorplaneVecZDCC.Dot(spectatorplaneVecZDCA), coll.centFT0C());

    std::vector<float> bdtScore[2];
    for (auto& casc : Cascades) {

      /// Add some minimal cuts for single track variables (min number of TPC clusters)
      auto negExtra = casc.negTrackExtra_as<DauTracks>();
      auto posExtra = casc.posTrackExtra_as<DauTracks>();
      auto bachExtra = casc.bachTrackExtra_as<DauTracks>();

      int counter = 0;
      IsCascAccepted(casc, negExtra, posExtra, bachExtra, counter);
      histos.fill(HIST("hCascade"), counter);

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

      // inv mass loose cut
      if (casc.pt() < MinPt || casc.pt() > MaxPt) {
        continue;
      }

      cascadev2::hMassBeforeSelVsPt[0]->Fill(massCasc[0], casc.pt());
      cascadev2::hMassBeforeSelVsPt[1]->Fill(massCasc[1], casc.pt());

      if (isApplyML) {
        // Retrieve model output and selection outcome
        isSelectedCasc[0] = mlResponseXi.isSelectedMl(inputFeaturesCasc, casc.pt(), bdtScore[0]);
        isSelectedCasc[1] = mlResponseOmega.isSelectedMl(inputFeaturesCasc, casc.pt(), bdtScore[1]);

        for (int iS{0}; iS < 2; ++iS) {
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
      auto cascminuspsiT0C = GetPhiInRange(casc.phi() - PsiT0C);
      auto v2CEP = TMath::Cos(2.0 * cascminuspsiT0C);
      ROOT::Math::XYZVector cascUvec{std::cos(casc.phi()), std::sin(casc.phi()), 0};
      auto v1SP_ZDCA = cascUvec.Dot(spectatorplaneVecZDCA);
      auto v1SP_ZDCC = cascUvec.Dot(spectatorplaneVecZDCC);
      auto v1EP_ZDCA = TMath::Cos(casc.phi() - coll.psiZDCA());
      auto v1EP_ZDCC = TMath::Cos(casc.phi() - coll.psiZDCC());
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
      for (int i{0}; i < 2; ++i) {
        cascadeVector[i].SetCoordinates(casc.px(), casc.py(), casc.pz(), masses[i]);
        ROOT::Math::Boost cascadeBoost{cascadeVector[i].BoostToCM()};
        auto boostedLambda{cascadeBoost(lambdaVector)};
        cosThetaStarLambda[i] = boostedLambda.Pz() / boostedLambda.P();
      }
      int ChargeIndex = 0;
      if (casc.sign() > 0)
        ChargeIndex = 1;
      double Pzs2Xi = cosThetaStarLambda[0] * std::sin(2 * (casc.phi() - PsiT0C)) / cascadev2::AlphaXi[ChargeIndex];
      double Pzs2Omega = cosThetaStarLambda[1] * std::sin(2 * (casc.phi() - PsiT0C)) / cascadev2::AlphaOmega[ChargeIndex];
      double Cos2ThetaXi = cosThetaStarLambda[0] * cosThetaStarLambda[0];
      double Cos2ThetaOmega = cosThetaStarLambda[1] * cosThetaStarLambda[1];
      double Pzs2LambdaFromCasc = cosThetaStarProton * std::sin(2 * (casc.phi() - PsiT0C)) / cascadev2::AlphaLambda[ChargeIndex];
      double Cos2ThetaLambda = cosThetaStarProton * cosThetaStarProton;

      double CosThetaXiWithAlpha = cosThetaStarLambda[0] / cascadev2::AlphaXi[ChargeIndex];
      double CosThetaOmegaWithAlpha = cosThetaStarLambda[1] / cascadev2::AlphaOmega[ChargeIndex];
      double CosThetaProtonWithAlpha = cosThetaStarProton / cascadev2::AlphaLambda[ChargeIndex];

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

      float BDTresponse[2]{0.f, 0.f};
      if (isApplyML) {
        BDTresponse[0] = bdtScore[0][1];
        BDTresponse[1] = bdtScore[1][1];
      }

      if (std::abs(casc.eta()) < etaCasc) {
        if (isFillTHNXi) {
          histos.get<THn>(HIST("hXiV2"))->Fill(coll.centFT0C(), ChargeIndex, casc.pt(), casc.mXi(), BDTresponse[0], v2CEP);
          histos.get<THn>(HIST("hXiPzs2"))->Fill(coll.centFT0C(), ChargeIndex, casc.pt(), casc.mXi(), BDTresponse[0], Pzs2Xi);
          histos.get<THn>(HIST("hXiPzs2FromLambda"))->Fill(coll.centFT0C(), ChargeIndex, casc.pt(), casc.mXi(), BDTresponse[0], Pzs2LambdaFromCasc);
          histos.get<THn>(HIST("hXiCos2Theta"))->Fill(coll.centFT0C(), ChargeIndex, casc.pt(), casc.mXi(), BDTresponse[0], Cos2ThetaXi);
          histos.get<THn>(HIST("hXiCos2ThetaFromLambda"))->Fill(coll.centFT0C(), ChargeIndex, casc.pt(), casc.mXi(), BDTresponse[0], Cos2ThetaLambda);
        }
        if (isFillTHNXi_PzVsPsi) {
          histos.get<THn>(HIST("hXiPzVsPsi"))->Fill(coll.centFT0C(), ChargeIndex, casc.pt(), casc.mXi(), BDTresponse[0], CosThetaXiWithAlpha, 2 * cascminuspsiT0C);
          histos.get<THn>(HIST("hXiPzVsPsiFromLambda"))->Fill(coll.centFT0C(), ChargeIndex, casc.pt(), casc.mXi(), BDTresponse[0], CosThetaProtonWithAlpha, 2 * cascminuspsiT0C);
          histos.get<THn>(HIST("hXiCos2ThetaVsPsi"))->Fill(coll.centFT0C(), ChargeIndex, casc.pt(), casc.mXi(), BDTresponse[0], Cos2ThetaXi, 2 * cascminuspsiT0C);
          histos.get<THn>(HIST("hXiCos2ThetaVsPsiFromLambda"))->Fill(coll.centFT0C(), ChargeIndex, casc.pt(), casc.mXi(), BDTresponse[0], Cos2ThetaLambda, 2 * cascminuspsiT0C);
        }
        if (isFillTHNOmega) {
          histos.get<THn>(HIST("hOmegaV2"))->Fill(coll.centFT0C(), ChargeIndex, casc.pt(), casc.mOmega(), BDTresponse[1], v2CEP);
          histos.get<THn>(HIST("hOmegaPzs2"))->Fill(coll.centFT0C(), ChargeIndex, casc.pt(), casc.mOmega(), BDTresponse[1], Pzs2Omega);
          histos.get<THn>(HIST("hOmegaPzs2FromLambda"))->Fill(coll.centFT0C(), ChargeIndex, casc.pt(), casc.mOmega(), BDTresponse[1], Pzs2LambdaFromCasc);
          histos.get<THn>(HIST("hOmegaCos2Theta"))->Fill(coll.centFT0C(), ChargeIndex, casc.pt(), casc.mOmega(), BDTresponse[1], Cos2ThetaOmega);
          histos.get<THn>(HIST("hOmegaCos2ThetaFromLambda"))->Fill(coll.centFT0C(), ChargeIndex, casc.pt(), casc.mOmega(), BDTresponse[1], Cos2ThetaLambda);
        }
        if (isFillTHNOmega_PzVsPsi) {
          histos.get<THn>(HIST("hOmegaPzVsPsi"))->Fill(coll.centFT0C(), ChargeIndex, casc.pt(), casc.mOmega(), BDTresponse[0], CosThetaOmegaWithAlpha, 2 * cascminuspsiT0C);
          histos.get<THn>(HIST("hOmegaPzVsPsiFromLambda"))->Fill(coll.centFT0C(), ChargeIndex, casc.pt(), casc.mOmega(), BDTresponse[0], CosThetaProtonWithAlpha, 2 * cascminuspsiT0C);
          histos.get<THn>(HIST("hOmegaCos2ThetaVsPsi"))->Fill(coll.centFT0C(), ChargeIndex, casc.pt(), casc.mOmega(), BDTresponse[0], Cos2ThetaOmega, 2 * cascminuspsiT0C);
          histos.get<THn>(HIST("hOmegaCos2ThetaVsPsiFromLambda"))->Fill(coll.centFT0C(), ChargeIndex, casc.pt(), casc.mOmega(), BDTresponse[0], Cos2ThetaLambda, 2 * cascminuspsiT0C);
        }
      }

      if (isSelectedCasc[0] || isSelectedCasc[1]) {
        if (isFillTree)
          fillAnalysedTable(coll, hasEventPlane, hasSpectatorPlane, casc, v2CSP, v2CEP, v1SP_ZDCA, v1SP_ZDCC, PsiT0C, BDTresponse[0], BDTresponse[1], 0);
      }
    }
  }

  void processAnalyseDataEPCentralFW(CollEventPlaneCentralFW const& coll, CascCandidates const& Cascades, DauTracks const&)
  {

    if (!AcceptEvent(coll, 1)) {
      return;
    }

    // select only events used for the calibration of the event plane
    if (isGoodEventEP) {
      if (abs(coll.qvecFT0CRe()) > 990 || abs(coll.qvecFT0CIm()) > 990 || abs(coll.qvecBNegRe()) > 990 || abs(coll.qvecBNegIm()) > 990 || abs(coll.qvecBPosRe()) > 990 || abs(coll.qvecBPosIm()) > 990) {
        return;
      }
    }

    // event has event plane
    bool hasEventPlane = 0;
    if (abs(coll.qvecFT0CRe()) > 990 || abs(coll.qvecFT0CIm()) > 990 || abs(coll.qvecBNegRe()) > 990 || abs(coll.qvecBNegIm()) > 990 || abs(coll.qvecBPosRe()) > 990 || abs(coll.qvecBPosIm()) > 990)
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

    float NormQvT0C = sqrt(eventplaneVecT0C.Dot(eventplaneVecT0C));
    float NormQvTPCA = sqrt(eventplaneVecTPCA.Dot(eventplaneVecTPCA));
    float NormQvTPCC = sqrt(eventplaneVecTPCC.Dot(eventplaneVecTPCC));

    const float PsiT0C = std::atan2(coll.qvecFT0CIm(), coll.qvecFT0CRe()) * 0.5f;
    histos.fill(HIST("hPsiT0C"), PsiT0C);
    histos.fill(HIST("hPsiT0CvsCentFT0C"), coll.centFT0C(), PsiT0C);

    resolution.fill(HIST("QVectorsT0CTPCA"), eventplaneVecT0C.Dot(eventplaneVecTPCA), coll.centFT0C());
    resolution.fill(HIST("QVectorsT0CTPCC"), eventplaneVecT0C.Dot(eventplaneVecTPCC), coll.centFT0C());
    resolution.fill(HIST("QVectorsTPCAC"), eventplaneVecTPCA.Dot(eventplaneVecTPCC), coll.centFT0C());
    resolution.fill(HIST("QVectorsNormT0CTPCA"), eventplaneVecT0C.Dot(eventplaneVecTPCA) / (NormQvT0C * NormQvTPCA), coll.centFT0C());
    resolution.fill(HIST("QVectorsNormT0CTPCC"), eventplaneVecT0C.Dot(eventplaneVecTPCC) / (NormQvT0C * NormQvTPCC), coll.centFT0C());
    resolution.fill(HIST("QVectorsNormTPCAC"), eventplaneVecTPCA.Dot(eventplaneVecTPCC) / (NormQvTPCA * NormQvTPCC), coll.centFT0C());
    resolution.fill(HIST("QVectorsSpecPlane"), spectatorplaneVecZDCC.Dot(spectatorplaneVecZDCA), coll.centFT0C());

    std::vector<float> bdtScore[2];
    for (auto& casc : Cascades) {

      /// Add some minimal cuts for single track variables (min number of TPC clusters)
      auto negExtra = casc.negTrackExtra_as<DauTracks>();
      auto posExtra = casc.posTrackExtra_as<DauTracks>();
      auto bachExtra = casc.bachTrackExtra_as<DauTracks>();

      int counter = 0;
      IsCascAccepted(casc, negExtra, posExtra, bachExtra, counter);
      histos.fill(HIST("hCascade"), counter);

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

      // inv mass loose cut
      if (casc.pt() < MinPt || casc.pt() > MaxPt) {
        continue;
      }

      cascadev2::hMassBeforeSelVsPt[0]->Fill(massCasc[0], casc.pt());
      cascadev2::hMassBeforeSelVsPt[1]->Fill(massCasc[1], casc.pt());

      if (isApplyML) {
        // Retrieve model output and selection outcome
        isSelectedCasc[0] = mlResponseXi.isSelectedMl(inputFeaturesCasc, casc.pt(), bdtScore[0]);
        isSelectedCasc[1] = mlResponseOmega.isSelectedMl(inputFeaturesCasc, casc.pt(), bdtScore[1]);

        for (int iS{0}; iS < 2; ++iS) {
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
      auto cascminuspsiT0C = GetPhiInRange(casc.phi() - PsiT0C);
      auto v2CEP = TMath::Cos(2.0 * cascminuspsiT0C);
      ROOT::Math::XYZVector cascUvec{std::cos(casc.phi()), std::sin(casc.phi()), 0};
      auto v1SP_ZDCA = cascUvec.Dot(spectatorplaneVecZDCA);
      auto v1SP_ZDCC = cascUvec.Dot(spectatorplaneVecZDCC);
      auto v1EP_ZDCA = TMath::Cos(casc.phi() - coll.psiZDCA());
      auto v1EP_ZDCC = TMath::Cos(casc.phi() - coll.psiZDCC());
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

      float BDTresponse[2]{0.f, 0.f};
      if (isApplyML) {
        BDTresponse[0] = bdtScore[0][1];
        BDTresponse[1] = bdtScore[1][1];
      }
      if (isSelectedCasc[0] || isSelectedCasc[1])
        if (isFillTree)
          fillAnalysedTable(coll, hasEventPlane, hasSpectatorPlane, casc, v2CSP, v2CEP, v1SP_ZDCA, v1SP_ZDCC, PsiT0C, BDTresponse[0], BDTresponse[1], 0);
    }
  }

  void processAnalyseMC(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels>::iterator const& coll, CascMCCandidates const& Cascades, DauTracks const&, soa::Join<aod::CascMCCores, aod::CascMCCollRefs> const&)
  {

    if (!AcceptEvent(coll, 1)) {
      return;
    }

    bool hasEventPlane = 0;     // no info at the moment
    bool hasSpectatorPlane = 0; // no infor at the moment

    histos.fill(HIST("hNEvents"), 9.5);
    histos.fill(HIST("hEventNchCorrelationAfterEP"), coll.multNTracksPVeta1(), coll.multNTracksGlobal());
    histos.fill(HIST("hEventPVcontributorsVsCentralityAfterEP"), coll.centFT0C(), coll.multNTracksPVeta1());
    histos.fill(HIST("hEventGlobalTracksVsCentralityAfterEP"), coll.centFT0C(), coll.multNTracksGlobal());
    histos.fill(HIST("hEventCentrality"), coll.centFT0C());
    histos.fill(HIST("hEventVertexZ"), coll.posZ());
    histos.fill(HIST("hMultNTracksITSTPCVsCentrality"), coll.centFT0C(), coll.multNTracksITSTPC());

    std::vector<float> bdtScore[2];
    for (auto& casc : Cascades) {

      if (!casc.has_cascMCCore())
        continue;

      auto cascMC = casc.cascMCCore_as<soa::Join<aod::CascMCCores, aod::CascMCCollRefs>>();

      int pdgCode{cascMC.pdgCode()};
      if (!(std::abs(pdgCode) == 3312 && std::abs(cascMC.pdgCodeV0()) == 3122 && std::abs(cascMC.pdgCodeBachelor()) == 211)     // Xi
          && !(std::abs(pdgCode) == 3334 && std::abs(cascMC.pdgCodeV0()) == 3122 && std::abs(cascMC.pdgCodeBachelor()) == 321)) // Omega
      {
        pdgCode = 0;
      }

      // rapidity definition
      float XiY = RecoDecay::y(std::array{casc.px(), casc.py(), casc.pz()}, constants::physics::MassXiMinus);
      float OmegaY = RecoDecay::y(std::array{casc.px(), casc.py(), casc.pz()}, constants::physics::MassOmegaMinus);
      // true reco cascades before applying any selection
      if (std::abs(pdgCode) == 3312 && std::abs(cascMC.pdgCodeV0()) == 3122 && std::abs(cascMC.pdgCodeBachelor()) == 211) {
        histos.fill(HIST("hXiPtvsCent"), coll.centFT0C(), casc.pt());
        if (std::abs(casc.eta()) < 0.8)
          histos.fill(HIST("hXiPtvsCentEta08"), coll.centFT0C(), casc.pt());
        if (std::abs(XiY) < 0.5)
          histos.fill(HIST("hXiPtvsCentY05"), coll.centFT0C(), casc.pt());
      } else if (std::abs(pdgCode) == 3334 && std::abs(cascMC.pdgCodeV0()) == 3122 && std::abs(cascMC.pdgCodeBachelor()) == 321) {
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
      IsCascAccepted(casc, negExtra, posExtra, bachExtra, counter);
      histos.fill(HIST("hCascade"), counter);

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

      if (casc.pt() < MinPt || casc.pt() > MaxPt) {
        continue;
      }

      cascadev2::hMassBeforeSelVsPt[0]->Fill(massCasc[0], casc.pt());
      cascadev2::hMassBeforeSelVsPt[1]->Fill(massCasc[1], casc.pt());

      if (isApplyML) {
        // Retrieve model output and selection outcome
        isSelectedCasc[0] = mlResponseXi.isSelectedMl(inputFeaturesCasc, casc.pt(), bdtScore[0]);
        isSelectedCasc[1] = mlResponseOmega.isSelectedMl(inputFeaturesCasc, casc.pt(), bdtScore[1]);

        for (int iS{0}; iS < 2; ++iS) {
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

      float BDTresponse[2]{0.f, 0.f};
      const float PsiT0C = 0; // not defined in MC for now
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
        fillAnalysedTable(coll, hasEventPlane, hasSpectatorPlane, casc, v2CSP, v2CEP, v1SP_ZDCA, v1SP_ZDCC, PsiT0C, BDTresponse[0], BDTresponse[1], pdgCode);
    }
  }

  void processMCGen(MCCollisionsStra::iterator const& mcCollision, const soa::SmallGroups<soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraCollLabels>>& collisions, const soa::SmallGroups<soa::Join<aod::CascMCCores, aod::CascMCCollRefs>>& cascMC)
  {

    histosMCGen.fill(HIST("hZvertexGen"), mcCollision.posZ());
    histosMCGen.fill(HIST("hNEventsMC"), 0.5);
    // Generated with accepted z vertex
    if (TMath::Abs(mcCollision.posZ()) > cutzvertex) {
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
      if (TMath::Abs(cascmc.pdgCode()) == 3312)
        histosMCGen.fill(HIST("hNCascGen"), 0.5);
      else if (TMath::Abs(cascmc.pdgCode()) == 3334)
        histosMCGen.fill(HIST("hNCascGen"), 1.5);
      if (!cascmc.has_straMCCollision())
        continue;
      if (TMath::Abs(cascmc.pdgCode()) == 3312)
        histosMCGen.fill(HIST("hNCascGen"), 2.5);
      else if (TMath::Abs(cascmc.pdgCode()) == 3334)
        histosMCGen.fill(HIST("hNCascGen"), 3.5);
      if (!cascmc.isPhysicalPrimary())
        continue;
      if (TMath::Abs(cascmc.pdgCode()) == 3312)
        histosMCGen.fill(HIST("hNCascGen"), 4.5);
      else if (TMath::Abs(cascmc.pdgCode()) == 3334)
        histosMCGen.fill(HIST("hNCascGen"), 5.5);

      float ptmc = RecoDecay::sqrtSumOfSquares(cascmc.pxMC(), cascmc.pyMC());

      float theta = std::atan(ptmc / cascmc.pzMC()); //-pi/2 < theta < pi/2

      float theta1 = 0;

      // if pz is positive (i.e. positive rapidity): 0 < theta < pi/2
      if (theta > 0)
        theta1 = theta; // 0 < theta1/2 < pi/4 --> 0 < tan (theta1/2) < 1 --> positive eta
      // if pz is negative (i.e. negative rapidity): -pi/2 < theta < 0 --> we need 0 < theta1/2 < pi/2 for the ln to be defined
      else
        theta1 = TMath::Pi() + theta; // pi/2 < theta1 < pi --> pi/4 < theta1/2 <  pi/2 --> 1 < tan (theta1/2) --> negative eta

      float cascMCeta = -log(std::tan(theta1 / 2));
      float cascMCy = 0;

      if (TMath::Abs(cascmc.pdgCode()) == 3312) {
        cascMCy = RecoDecay::y(std::array{cascmc.pxMC(), cascmc.pyMC(), cascmc.pzMC()}, constants::physics::MassXiMinus);
        if (TMath::Abs(cascMCeta) < etaCascMCGen) {
          histosMCGen.fill(HIST("h2DGenXiEta08"), centrality, ptmc);
          histosMCGen.fill(HIST("hNCascGen"), 6.5);
        }
        if (TMath::Abs(cascMCy) < yCascMCGen)
          histosMCGen.fill(HIST("h2DGenXiY05"), centrality, ptmc);
        histosMCGen.fill(HIST("hGenXiY"), cascMCy);
      } else if (TMath::Abs(cascmc.pdgCode() == 3334)) {
        cascMCy = RecoDecay::y(std::array{cascmc.pxMC(), cascmc.pyMC(), cascmc.pzMC()}, constants::physics::MassOmegaMinus);
        if (TMath::Abs(cascMCeta) < etaCascMCGen) {
          histosMCGen.fill(HIST("h2DGenOmegaEta08"), centrality, ptmc);
          histosMCGen.fill(HIST("hNCascGen"), 7.5);
        }
        if (TMath::Abs(cascMCy) < yCascMCGen)
          histosMCGen.fill(HIST("h2DGenOmegaY05"), centrality, ptmc);
        histosMCGen.fill(HIST("hGenOmegaY"), cascMCy);
      }
    }
  }

  PROCESS_SWITCH(cascadeFlow, processTrainingBackground, "Process to create the training dataset for the background", true);
  PROCESS_SWITCH(cascadeFlow, processTrainingSignal, "Process to create the training dataset for the signal", false);
  PROCESS_SWITCH(cascadeFlow, processAnalyseData, "Process to apply ML model to the data", false);
  PROCESS_SWITCH(cascadeFlow, processAnalyseDataEPCentralFW, "Process to apply ML model to the data - event plane calibration from central framework", false);
  PROCESS_SWITCH(cascadeFlow, processAnalyseMC, "Process to apply ML model to the MC", false);
  PROCESS_SWITCH(cascadeFlow, processMCGen, "Process to store MC generated particles", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<cascadeFlow>(cfgc, TaskName{"lf-cascade-flow"})};
}
