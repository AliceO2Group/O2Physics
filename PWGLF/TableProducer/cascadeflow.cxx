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
using CollEventPlane = soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraFT0AQVs, aod::StraFT0CQVs, aod::StraFT0MQVs>::iterator;

ROOT::Math::XYZVector eventplaneVecT0A, eventplaneVecT0C, cascphiVec;

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
constexpr int cutDir[nCutScores] = {CutGreater, CutNot};
auto vecCutDir = std::vector<int>{cutDir, cutDir + nCutScores};

// default values for the cuts
constexpr double cuts[nBinsPt][nCutScores] = {
  {0.5, 0.5},
  {0.5, 0.5},
  {0.5, 0.5},
  {0.5, 0.5},
  {0.5, 0.5},
  {0.5, 0.5},
  {0.5, 0.5},
  {0.5, 0.5}};

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
static const std::vector<std::string> labelsCutScore = {"score class 1", "score class 2"};
static const std::vector<std::string> labelsDmesCutScore = {"ML score signal", "ML score bkg"};
} // namespace cascade_flow_cuts_ml

struct cascadeFlow {

  Configurable<bool> isOmega{"isOmega", 0, "Xi or Omega"};
  Configurable<float> MinPt{"MinPt", 0.6, "Min pt of cascade"};
  Configurable<float> MaxPt{"MaxPt", 10, "Max pt of cascade"};
  Configurable<bool> isApplyML{"isApplyML", 1, ""};
  Configurable<double> sideBandStart{"sideBandStart", 5, "Start of the sideband region in number of sigmas"};
  Configurable<double> sideBandEnd{"sideBandEnd", 7, "End of the sideband region in number of sigmas"};
  Configurable<double> downsample{"downsample", 1., "Downsample training output tree"};
  Configurable<bool> doNTPCSigmaCut{"doNTPCSigmaCut", 1, "doNtpcSigmaCut"};
  Configurable<float> nsigmatpcPr{"nsigmatpcPr", 5, "nsigmatpcPr"};
  Configurable<float> nsigmatpcPi{"nsigmatpcPi", 5, "nsigmatpcPi"};
  Configurable<float> mintpccrrows{"mintpccrrows", 3, "mintpccrrows"};

  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::vector<std::string>> modelPathsCCDB{"modelPathsCCDB", std::vector<std::string>{"Users/c/chdemart/CascadesFlow"}, "Paths of models on CCDB"};
  Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"model_onnx.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};

  // ML inference
  Configurable<bool> applyMl{"applyMl", true, "Flag to apply ML selections"};
  Configurable<std::vector<double>> binsPtMl{"binsPtMl", std::vector<double>{cascade_flow_cuts_ml::vecBinsPt}, "pT bin limits for ML application"};
  Configurable<std::vector<int>> cutDirMl{"cutDirMl", std::vector<int>{cascade_flow_cuts_ml::vecCutDir}, "Whether to reject score values greater or smaller than the threshold"};
  Configurable<LabeledArray<double>> cutsMl{"cutsMl", {cascade_flow_cuts_ml::cuts[0], cascade_flow_cuts_ml::nBinsPt, cascade_flow_cuts_ml::nCutScores, cascade_flow_cuts_ml::labelsPt, cascade_flow_cuts_ml::labelsCutScore}, "ML selections per pT bin"};
  Configurable<int8_t> nClassesMl{"nClassesMl", (int8_t)cascade_flow_cuts_ml::nCutScores, "Number of classes in ML model"};
  //  Configurable<std::vector<std::string>> namesInputFeatures{"namesInputFeatures", std::vector<std::string>{"feature1", "feature2"}, "Names of ML model input features"};

  o2::ccdb::CcdbApi ccdbApi;

  // Add objects needed for ML inference
  o2::analysis::MlResponse<float> mlResponse;
  std::vector<float> outputMl = {};

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
    double result = phi;
    while (result < 0) {
      result = result + 2. * TMath::Pi() / 2;
    }
    while (result > 2. * TMath::Pi() / 2) {
      result = result - 2. * TMath::Pi() / 2;
    }
    return result;
  }

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

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
    trainingSample(coll.centFT0M(),
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
  //  void fillAnalisedTable(collision_t coll, cascade_t casc, float BDTresponse)
  void fillAnalysedTable(collision_t coll, cascade_t casc, float v2A, float v2C, float BDTresponse)
  {
    analysisSample(coll.centFT0M(),
                   casc.sign(),
                   casc.pt(),
                   casc.eta(),
                   casc.mXi(),
                   casc.mOmega(),
                   v2A,
                   v2C,
                   BDTresponse);
  }

  void init(InitContext const&)
  {

    float MinMass = 1.28;
    if (isOmega)
      MinMass = 1.6;
    float MaxMass = 1.36;
    if (isOmega)
      MaxMass = 1.73;
    ConfigurableAxis vertexZ{"vertexZ", {20, -10, 10}, "vertex axis for histograms"};
    ConfigurableAxis massCascAxis{"massCascAxis", {100, MinMass, MaxMass}, ""};
    ConfigurableAxis ptAxis{"ptAxis", {100, 0, 10}, ""};
    histos.add("hEventVertexZ", "hEventVertexZ", kTH1F, {vertexZ});
    histos.add("hEventCentrality", "hEventCentrality", kTH1F, {{101, 0, 101}});
    histos.add("hCandidate", "hCandidate", HistType::kTH1D, {{22, -0.5, 21.5}});
    histos.add("hCascadeSignal", "hCascadeSignal", HistType::kTH1D, {{6, -0.5, 5.5}});
    histos.add("hCascade", "hCascade", HistType::kTH1D, {{6, -0.5, 5.5}});
    histos.add("hCascadePhi", "hCascadePhi", HistType::kTH1D, {{100, 0, 2 * TMath::Pi()}});
    histos.add("hcascminuspsiT0A", "hcascminuspsiT0A", HistType::kTH1D, {{100, 0, 2 * TMath::Pi()}});
    histos.add("hcascminuspsiT0C", "hcascminuspsiT0C", HistType::kTH1D, {{100, 0, 2 * TMath::Pi()}});
    histos.add("hPsiT0A", "hPsiT0A", HistType::kTH1D, {{100, 0, 2 * TMath::Pi()}});
    histos.add("hPsiT0C", "hPsiT0C", HistType::kTH1D, {{100, 0, 2 * TMath::Pi()}});
    histos.add("hFT0ARe", "hFT0ARe", HistType::kTH1D, {{100, -1, 1}});
    histos.add("hFT0AIm", "hFT0AIm", HistType::kTH1D, {{100, -1, 1}});
    histos.add("hFT0A", "hFT0A", HistType::kTH2D, {{100, -1, 1}, {100, -1, 1}});
    histos.add("hv2Norm", "hv2Norm", HistType::kTH1D, {{100, 0, 1}});
    histos.add("hMassBeforeSel", "hMassBeforeSel", HistType::kTH1F, {massCascAxis});
    histos.add("hMassAfterSel", "hMassAfterSel", HistType::kTH1F, {massCascAxis});
    histos.add("hMassBeforeSelVsPt", "hMassBeforeSelVsPt", HistType::kTH2F, {{massCascAxis}, {ptAxis}});
    histos.add("hMassAfterSelVsPt", "hMassAfterSelVsPt", HistType::kTH2F, {{massCascAxis}, {ptAxis}});
    histos.add("hSignalScoreBeforeSel", "Signal score before selection;BDT first score;entries", HistType::kTH1F, {{100, 0., 1.}});
    histos.add("hBkgScoreBeforeSel", "Bkg score before selection;BDT first score;entries", HistType::kTH1F, {{100, 0., 1.}});
    histos.add("hSignalScoreAfterSel", "Signal score after selection;BDT first score;entries", HistType::kTH1F, {{100, 0., 1.}});
    histos.add("hBkgScoreAfterSel", "Bkg score after selection;BDT first score;entries", HistType::kTH1F, {{100, 0., 1.}});

    if (isApplyML) {
      // Configure and initialise the ML class
      mlResponse.configure(binsPtMl, cutsMl, cutDirMl, nClassesMl);
      // Bonus: retrieve the model from CCDB (needed for ML application on the GRID)
      if (loadModelsFromCCDB) {
        ccdbApi.init(ccdbUrl);
        mlResponse.setModelPathsCCDB(onnxFileNames, ccdbApi, modelPathsCCDB, timestampCCDB);
      } else {
        mlResponse.setModelPathsLocal(onnxFileNames);
      }
      mlResponse.init();
    }
  }

  void processTrainingBackground(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels>::iterator const& coll, soa::Join<aod::CascCollRefs, aod::CascCores, aod::CascExtras, aod::CascBBs> const& Cascades, DauTracks const&)
  {

    int counter = 0;

    if (!coll.sel8() || std::abs(coll.posZ()) > 10.) {
      return;
    }

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
      } else
        ++counter;

      if (posExtra.tpcCrossedRows() < mintpccrrows || negExtra.tpcCrossedRows() < mintpccrrows || bachExtra.tpcCrossedRows() < mintpccrrows)
        continue;
      histos.fill(HIST("hCandidate"), ++counter);

      fillTrainingTable(coll, casc, 0);
    }
  }

  void processTrainingSignal(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels>::iterator const& coll, soa::Join<aod::CascCollRefs, aod::CascCores, aod::CascMCCores, aod::CascExtras, aod::CascBBs> const& Cascades, DauTracks const&)
  {

    if (!coll.sel8() || std::abs(coll.posZ()) > 10.) {
      return;
    }

    for (auto& casc : Cascades) {
      int pdgCode{casc.pdgCode()};
      if (!(std::abs(pdgCode) == 3312 && std::abs(casc.pdgCodeV0()) == 3122 && std::abs(casc.pdgCodeBachelor()) == 211)     // Xi
          && !(std::abs(pdgCode) == 3334 && std::abs(casc.pdgCodeV0()) == 3122 && std::abs(casc.pdgCodeBachelor()) == 321)) // Omega
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

  void processAnalyseData(CollEventPlane const& coll, soa::Join<aod::CascCollRefs, aod::CascCores, aod::CascExtras, aod::CascBBs> const& Cascades, DauTracks const&)
  {

    // std::cout << "onnxFileNames " << onnxFileNames.value[0] << std::endl;

    if (!coll.sel8() || std::abs(coll.posZ()) > 10.) {
      return;
    }

    histos.fill(HIST("hEventCentrality"), coll.centFT0M());
    histos.fill(HIST("hEventVertexZ"), coll.posZ());

    eventplaneVecT0A = ROOT::Math::XYZVector(coll.qvecFT0ARe(), coll.qvecFT0AIm(), 0);
    eventplaneVecT0C = ROOT::Math::XYZVector(coll.qvecFT0CRe(), coll.qvecFT0CIm(), 0);
    float eventplaneVecT0ANorm = sqrt(eventplaneVecT0A.Dot(eventplaneVecT0A));
    float eventplaneVecT0CNorm = sqrt(eventplaneVecT0C.Dot(eventplaneVecT0C));
    //    float v2Norm = eventplaneVecT0A.Dot(eventplaneVecT0C)/coll.sumAmplFT0A()/coll.sumAmplFT0C();
    float v2Norm = eventplaneVecT0A.Dot(eventplaneVecT0C) / eventplaneVecT0ANorm / eventplaneVecT0CNorm;
    // std::cout << "eventplaneVecT0ANorm: " << eventplaneVecT0ANorm <<  " = " << sqrt(pow(coll.qvecFT0ARe(), 2) + pow(coll.qvecFT0AIm(), 2)) << std::endl; //these two are equal
    // std::cout << "stored norm value: " << coll.sumAmplFT0A() << std::endl;
    histos.fill(HIST("hv2Norm"), v2Norm);

    float PsiT0A = TMath::ACos(coll.qvecFT0ARe()) / 2;
    float PsiT0C = TMath::ACos(coll.qvecFT0CRe()) / 2;
    if (coll.qvecFT0AIm() < 0)
      PsiT0A = -PsiT0A + TMath::Pi() / 2; // to get dstribution between 0 and Pi
    if (coll.qvecFT0CIm() < 0)
      PsiT0C = -PsiT0C + TMath::Pi() / 2; // to get dstribution between 0 and Pi

    histos.fill(HIST("hFT0ARe"), coll.qvecFT0ARe());
    histos.fill(HIST("hFT0AIm"), coll.qvecFT0AIm());
    histos.fill(HIST("hFT0A"), coll.qvecFT0ARe(), coll.qvecFT0AIm());
    histos.fill(HIST("hPsiT0A"), PsiT0A);
    histos.fill(HIST("hPsiT0C"), PsiT0C);

    for (auto& casc : Cascades) {

      /// Add some minimal cuts for single track variables (min number of TPC clusters)
      auto negExtra = casc.negTrackExtra_as<DauTracks>();
      auto posExtra = casc.posTrackExtra_as<DauTracks>();
      auto bachExtra = casc.bachTrackExtra_as<DauTracks>();

      int counter = 0;
      IsCascAccepted(casc, negExtra, posExtra, bachExtra, counter);
      histos.fill(HIST("hCascade"), counter);

      // ML selections
      bool isSelectedCasc = false;

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

      float massCasc = casc.mXi();
      if (isOmega)
        massCasc = casc.mOmega();

      // inv mass loose cut
      if (casc.pt() < MinPt)
        continue;
      if (casc.pt() > MaxPt)
        continue;
      if (isOmega) {
        if (massCasc < 1.6 || massCasc > 1.73)
          continue;
      } else {
        if (massCasc < 1.28 || massCasc > 1.36)
          continue;
      }

      histos.fill(HIST("hMassBeforeSel"), massCasc);
      histos.fill(HIST("hMassBeforeSelVsPt"), massCasc, casc.pt());

      if (isApplyML) {
        // Retrieve model output and selection outcome
        isSelectedCasc = mlResponse.isSelectedMl(inputFeaturesCasc, casc.pt(), outputMl);
        // Fill BDT score histograms before selection
        histos.fill(HIST("hSignalScoreBeforeSel"), outputMl[0]);
        histos.fill(HIST("hBkgScoreBeforeSel"), outputMl[1]);

        // Fill histograms for selected candidates
        if (isSelectedCasc) {
          histos.fill(HIST("hSignalScoreAfterSel"), outputMl[0]);
          histos.fill(HIST("hBkgScoreAfterSel"), outputMl[1]);
        }

        outputMl.clear(); // not necessary in this case but for good measure
      } else {
        isSelectedCasc = true;
      }

      if (isSelectedCasc) {
        histos.fill(HIST("hMassAfterSel"), massCasc);
        histos.fill(HIST("hMassAfterSelVsPt"), massCasc, casc.pt());
      }

      cascphiVec = ROOT::Math::XYZVector(cos(2 * casc.phi()), sin(2 * casc.phi()), 0);
      auto v2A = cascphiVec.Dot(eventplaneVecT0A) / sqrt(eventplaneVecT0A.Dot(eventplaneVecT0A));
      auto v2C = cascphiVec.Dot(eventplaneVecT0C) / sqrt(eventplaneVecT0C.Dot(eventplaneVecT0C));
      auto cascminuspsiT0A = GetPhiInRange(casc.phi() - PsiT0A);
      auto cascminuspsiT0C = GetPhiInRange(casc.phi() - PsiT0C);
      // float v2A_diff =  TMath::Cos(2.0 * cascminuspsiT0A);
      // float v2C_diff =  TMath::Cos(2.0 * cascminuspsiT0C);
      // std::cout << "v2A: " << v2A << " v2C " << v2C << std::endl;
      // std::cout << "v2A_diff: " << v2A_diff << " v2C_diff " << v2C_diff << std::endl;
      histos.fill(HIST("hCascadePhi"), casc.phi());
      histos.fill(HIST("hcascminuspsiT0A"), cascminuspsiT0A);
      histos.fill(HIST("hcascminuspsiT0C"), cascminuspsiT0C);

      float BDTresponse = 0;
      if (isApplyML)
        BDTresponse = outputMl[0];
      if (isSelectedCasc)
        fillAnalysedTable(coll, casc, v2A, v2C, BDTresponse);
    }
  }

  PROCESS_SWITCH(cascadeFlow, processTrainingBackground, "Process to create the training dataset for the background", true);
  PROCESS_SWITCH(cascadeFlow, processTrainingSignal, "Process to create the training dataset for the signal", false);
  PROCESS_SWITCH(cascadeFlow, processAnalyseData, "Process to apply ML model to the data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<cascadeFlow>(cfgc, TaskName{"lf-cascade-flow"})};
}
