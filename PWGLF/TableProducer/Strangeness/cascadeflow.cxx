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
using CollEventPlaneCentralFW = soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraFT0CQVs, aod::StraTPCQVs>::iterator;
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

std::shared_ptr<TH2> hMassBeforeSelVsPt[2];
std::shared_ptr<TH2> hMassAfterSelVsPt[2];
std::shared_ptr<TH1> hSignalScoreBeforeSel[2];
std::shared_ptr<TH1> hBkgScoreBeforeSel[2];
std::shared_ptr<TH1> hSignalScoreAfterSel[2];
std::shared_ptr<TH1> hBkgScoreAfterSel[2];
std::shared_ptr<THnSparse> hSparseV2C[2];
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

  // axes
  ConfigurableAxis axisQVs{"axisQVs", {500, -10.f, 10.f}, "axisQVs"};
  ConfigurableAxis axisQVsNorm{"axisQVsNorm", {200, -1.f, 1.f}, "axisQVsNorm"};

  // Event selection criteria
  Configurable<bool> isStoreTrueCascOnly{"isStoreTrueCascOnly", 1, ""};
  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range (cm)"};
  Configurable<bool> sel8{"sel8", 1, "Apply sel8 event selection"};
  Configurable<bool> isNoSameBunchPileupCut{"isNoSameBunchPileupCut", 1, "Same found-by-T0 bunch crossing rejection"};
  Configurable<bool> isGoodZvtxFT0vsPVCut{"isGoodZvtxFT0vsPVCut", 1, "z of PV by tracks and z of PV from FT0 A-C time difference cut"};
  Configurable<bool> isGoodEventEP{"isGoodEventEP", 1, "Event is used to calibrate event plane"};
  Configurable<int> MinOccupancy{"MinOccupancy", 0, "MinOccupancy"};
  Configurable<int> MaxOccupancy{"MaxOccupancy", 500, "MaxOccupancy"};
  Configurable<bool> isCollInStandardTimeRange{"isCollInStandardTimeRange", 1, "To remove collisions in +-10 micros time range"};
  Configurable<bool> isCollInNarrowTimeRange{"isCollInNarrowTimeRange", 0, "To remove collisions in +-4 micros time range"};

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
    if (occupancy < MinOccupancy || occupancy > MaxOccupancy) {
      return false;
    }
    if (isFillHisto)
      histos.fill(HIST("hNEvents"), 5.5);

    if (isCollInStandardTimeRange && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return false;
    }

    if (isCollInNarrowTimeRange && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow)) {
      return false;
    }
    if (isFillHisto)
      histos.fill(HIST("hNEvents"), 6.5);

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
  void fillAnalysedTable(collision_t coll, cascade_t casc, float v2CSP, float v2CEP, float PsiT0C, float BDTresponseXi, float BDTresponseOmega, int pdgCode)
  {
    double masses[2]{o2::constants::physics::MassXiMinus, o2::constants::physics::MassOmegaMinus};
    ROOT::Math::PxPyPzMVector cascadeVector[2], lambdaVector, protonVector;
    float cosThetaStarLambda[2], cosThetaStarProton;
    lambdaVector.SetCoordinates(casc.pxlambda(), casc.pylambda(), casc.pzlambda(), o2::constants::physics::MassLambda);
    ROOT::Math::Boost lambdaBoost{lambdaVector.BoostToCM()};
    if (casc.sign() < 0) {
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

    analysisSample(coll.centFT0C(),
                   casc.sign(),
                   casc.pt(),
                   casc.eta(),
                   casc.phi(),
                   casc.mXi(),
                   casc.mOmega(),
                   v2CSP,
                   v2CEP,
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
    TString hNEventsLabels[8] = {"All", "sel8", "z vrtx", "kNoSameBunchPileup", "kIsGoodZvtxFT0vsPV", "trackOccupancyInTimeRange", "kNoCollInTimeRange", "kIsGoodEventEP"};
    TString hNEventsLabelsMC[6] = {"All", "z vtx", ">=1RecoColl", "1Reco", "2Reco", "EvSelected"};
    TString hNCascLabelsMC[8] = {"All Xi", "all Omega", "Xi: has MC coll", "Om: has MC coll", "Xi: isPrimary", "Om: is Primary", "Xi: |eta|<0.8", "Om: |eta| < 0.8"};

    resolution.add("QVectorsT0CTPCA", "QVectorsT0CTPCA", HistType::kTH2F, {axisQVs, CentAxis});
    resolution.add("QVectorsT0CTPCC", "QVectorsT0CTPCC", HistType::kTH2F, {axisQVs, CentAxis});
    resolution.add("QVectorsTPCAC", "QVectorsTPCAC", HistType::kTH2F, {axisQVs, CentAxis});
    resolution.add("QVectorsNormT0CTPCA", "QVectorsNormT0CTPCA", HistType::kTH2F, {axisQVsNorm, CentAxis});
    resolution.add("QVectorsNormT0CTPCC", "QVectorsNormT0CTPCC", HistType::kTH2F, {axisQVsNorm, CentAxis});
    resolution.add("QVectorsNormTPCAC", "QVectorsNormTPCCB", HistType::kTH2F, {axisQVsNorm, CentAxis});

    histos.add("hNEvents", "hNEvents", {HistType::kTH1F, {{8, 0.f, 8.f}}});
    for (Int_t n = 1; n <= histos.get<TH1>(HIST("hNEvents"))->GetNbinsX(); n++) {
      histos.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(n, hNEventsLabels[n - 1]);
    }
    histos.add("hEventVertexZ", "hEventVertexZ", kTH1F, {{120, -12., 12.}});
    histos.add("hEventCentrality", "hEventCentrality", kTH1F, {{101, 0, 101}});
    histos.add("hPsiT0C", "hPsiT0C", HistType::kTH1D, {{100, -TMath::Pi(), TMath::Pi()}});
    histos.add("hPsiT0CvsCentFT0C", "hPsiT0CvsCentFT0C", HistType::kTH2D, {CentAxis, {100, -TMath::Pi(), TMath::Pi()}});
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
    histos.add("hXiPtvsCent", "hXiPtvsCent", HistType::kTH2F, {{100, 0, 100}, {200, 0, 20}});
    histos.add("hOmegaPtvsCent", "hOmegaPtvsCent", HistType::kTH2F, {{100, 0, 100}, {200, 0, 20}});
    histos.add("hCascadePhi", "hCascadePhi", HistType::kTH1F, {{100, 0, 2 * TMath::Pi()}});
    histos.add("hcascminuspsiT0C", "hcascminuspsiT0C", HistType::kTH1F, {{100, 0, TMath::Pi()}});
    histos.add("hv2CEPvsFT0C", "hv2CEPvsFT0C", HistType::kTH2F, {CentAxis, {100, -1, 1}});
    histos.add("hv2CEPvsv2CSP", "hv2CEPvsV2CSP", HistType::kTH2F, {{100, -1, 1}, {100, -1, 1}});

    histosMCGen.add("h2DGenXiEta08", "h2DGenXiEta08", HistType::kTH2F, {{100, 0, 100}, {200, 0, 20}});
    histosMCGen.add("h2DGenOmegaEta08", "h2DGenOmegaEta08", HistType::kTH2F, {{100, 0, 100}, {200, 0, 20}});
    histosMCGen.add("h2DGenXiY05", "h2DGenXiY05", HistType::kTH2F, {{100, 0, 100}, {200, 0, 20}});
    histosMCGen.add("h2DGenOmegaY05", "h2DGenOmegaY05", HistType::kTH2F, {{100, 0, 100}, {200, 0, 20}});
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
      cascadev2::hSparseV2C[iS] = histos.add<THnSparse>(Form("hSparseV2C%s", cascadev2::speciesNames[iS].data()), "hSparseV2C", HistType::kTHnSparseF, {massCascAxis[iS], ptAxis, v2Axis, CentAxis});
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

  void processAnalyseData(CollEventPlane const& coll, CascCandidates const& Cascades, DauTracks const&)
  {

    if (!AcceptEvent(coll, 1)) {
      return;
    }

    // select only events used for the calibration of the event plane
    if (isGoodEventEP && !coll.triggereventep()) {
      return;
    }
    histos.fill(HIST("hNEvents"), 7.5);
    histos.fill(HIST("hEventNchCorrelationAfterEP"), coll.multNTracksPVeta1(), coll.multNTracksGlobal());
    histos.fill(HIST("hEventPVcontributorsVsCentralityAfterEP"), coll.centFT0C(), coll.multNTracksPVeta1());
    histos.fill(HIST("hEventGlobalTracksVsCentralityAfterEP"), coll.centFT0C(), coll.multNTracksGlobal());

    histos.fill(HIST("hEventCentrality"), coll.centFT0C());
    histos.fill(HIST("hEventVertexZ"), coll.posZ());

    ROOT::Math::XYZVector eventplaneVecT0C{coll.qvecFT0CRe(), coll.qvecFT0CIm(), 0};
    ROOT::Math::XYZVector eventplaneVecTPCA{coll.qvecBPosRe(), coll.qvecBPosIm(), 0};
    ROOT::Math::XYZVector eventplaneVecTPCC{coll.qvecBNegRe(), coll.qvecBNegIm(), 0};

    const float PsiT0C = std::atan2(coll.qvecFT0CIm(), coll.qvecFT0CRe()) * 0.5f;
    histos.fill(HIST("hPsiT0C"), PsiT0C);
    histos.fill(HIST("hPsiT0CvsCentFT0C"), coll.centFT0C(), PsiT0C);

    resolution.fill(HIST("QVectorsT0CTPCA"), eventplaneVecT0C.Dot(eventplaneVecTPCA), coll.centFT0C());
    resolution.fill(HIST("QVectorsT0CTPCC"), eventplaneVecT0C.Dot(eventplaneVecTPCC), coll.centFT0C());
    resolution.fill(HIST("QVectorsTPCAC"), eventplaneVecTPCA.Dot(eventplaneVecTPCC), coll.centFT0C());
    resolution.fill(HIST("QVectorsNormT0CTPCA"), eventplaneVecT0C.Dot(eventplaneVecTPCA) / (coll.qTPCR() * coll.sumAmplFT0C()), coll.centFT0C());
    resolution.fill(HIST("QVectorsNormT0CTPCC"), eventplaneVecT0C.Dot(eventplaneVecTPCC) / (coll.qTPCL() * coll.sumAmplFT0C()), coll.centFT0C());
    resolution.fill(HIST("QVectorsNormTPCAC"), eventplaneVecTPCA.Dot(eventplaneVecTPCC) / (coll.qTPCR() * coll.qTPCL()), coll.centFT0C());

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

      float BDTresponse[2]{0.f, 0.f};
      if (isApplyML) {
        BDTresponse[0] = bdtScore[0][1];
        BDTresponse[1] = bdtScore[1][1];
      }
      if (isSelectedCasc[0] || isSelectedCasc[1])
        fillAnalysedTable(coll, casc, v2CSP, v2CEP, PsiT0C, BDTresponse[0], BDTresponse[1], 0);
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

    histos.fill(HIST("hNEvents"), 7.5);
    histos.fill(HIST("hEventNchCorrelationAfterEP"), coll.multNTracksPVeta1(), coll.multNTracksGlobal());
    histos.fill(HIST("hEventPVcontributorsVsCentralityAfterEP"), coll.centFT0C(), coll.multNTracksPVeta1());
    histos.fill(HIST("hEventGlobalTracksVsCentralityAfterEP"), coll.centFT0C(), coll.multNTracksGlobal());

    histos.fill(HIST("hEventCentrality"), coll.centFT0C());
    histos.fill(HIST("hEventVertexZ"), coll.posZ());

    ROOT::Math::XYZVector eventplaneVecT0C{coll.qvecFT0CRe(), coll.qvecFT0CIm(), 0};
    ROOT::Math::XYZVector eventplaneVecTPCA{coll.qvecBPosRe(), coll.qvecBPosIm(), 0};
    ROOT::Math::XYZVector eventplaneVecTPCC{coll.qvecBNegRe(), coll.qvecBNegIm(), 0};
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

      float BDTresponse[2]{0.f, 0.f};
      if (isApplyML) {
        BDTresponse[0] = bdtScore[0][1];
        BDTresponse[1] = bdtScore[1][1];
      }
      if (isSelectedCasc[0] || isSelectedCasc[1])
        fillAnalysedTable(coll, casc, v2CSP, v2CEP, PsiT0C, BDTresponse[0], BDTresponse[1], 0);
    }
  }

  void processAnalyseMC(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels>::iterator const& coll, CascMCCandidates const& Cascades, DauTracks const&, soa::Join<aod::CascMCCores, aod::CascMCCollRefs> const&)
  {

    if (!AcceptEvent(coll, 1)) {
      return;
    }

    histos.fill(HIST("hNEvents"), 7.5);
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

      // true reco cascades before applying any selection
      if (std::abs(pdgCode) == 3312 && std::abs(cascMC.pdgCodeV0()) == 3122 && std::abs(cascMC.pdgCodeBachelor()) == 211) {
        histos.fill(HIST("hXiPtvsCent"), coll.centFT0C(), casc.pt());
      } else if (std::abs(pdgCode) == 3334 && std::abs(cascMC.pdgCodeV0()) == 3122 && std::abs(cascMC.pdgCodeBachelor()) == 321) {
        histos.fill(HIST("hOmegaPtvsCent"), coll.centFT0C(), casc.pt());
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

      if (isApplyML) {
        BDTresponse[0] = bdtScore[0][1];
        BDTresponse[1] = bdtScore[1][1];
      }
      if (isStoreTrueCascOnly) {
        if (pdgCode == 0)
          continue;
      }
      if (isSelectedCasc[0] || isSelectedCasc[1])
        fillAnalysedTable(coll, casc, v2CSP, v2CEP, PsiT0C, BDTresponse[0], BDTresponse[1], pdgCode);
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
    int bestCollisionIndex = -1;
    float centrality = 100.5f;
    int nCollisions = 0;
    for (auto const& coll : collisions) {
      if (!AcceptEvent(coll, 0)) {
        continue;
      }
      if (biggestNContribs < coll.multPVTotalContributors()) {
        biggestNContribs = coll.multPVTotalContributors();
        bestCollisionIndex = coll.globalIndex();
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
