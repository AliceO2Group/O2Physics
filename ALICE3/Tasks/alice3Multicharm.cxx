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

/// \file alice3Multicharm.cxx
/// \brief consumer task for alice 3 multicharm studies
/// \author Jesper Karlsson Gumprecht <jesper.gumprecht@cern.ch>

//  *+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
//   Decay finder task for ALICE 3
//  *+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
//
//    Uses specific ALICE 3 PID and performance for studying
//    HF decays. Work in progress: use at your own risk!

#include "ALICE3/DataModel/OTFMulticharm.h"
#include "ALICE3/ML/MulticharmMlResponse.h"
#include "Tools/ML/MlResponse.h"

#include <CCDB/CcdbApi.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct Alice3Multicharm {
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  o2::analysis::MulticharmMlResponse<float> mlResponse;
  o2::ccdb::CcdbApi ccdbApi;
  static constexpr float ToMicrons = 1e+4;

  struct : ConfigurableGroup {
    std::string prefix = "ml";
    Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> localPath{"localPath", "MCharm_BDTModel.onnx", "(std::string) Path to the local .onnx file."};
    Configurable<int64_t> timestampCCDB{"timestampCCDB", 1695750420200, "timestamp of the ONNX file for ML model used to query in CCDB. Please use 1695750420200"};
    Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};
    Configurable<bool> enableBDT{"enableBDT", false, "Enables bdt model"};
    Configurable<std::vector<std::string>> modelPathsCCDB{"modelPathsCCDB", std::vector<std::string>{"Users/j/jekarlss/MLModels"}, "Paths of models on CCDB"};
    Configurable<std::vector<double>> ptBinEdges{"ptBinEdges", {0, 2, 4, 6, 8, 10, 15}, "Bin edges for pT dependant BDT"};
    Configurable<LabeledArray<double>> scoreCuts{"scoreCuts", {multi_charm_ml::Cuts[0].data(), multi_charm_ml::NBinsPt, multi_charm_ml::NClasses, multi_charm_ml::labelsPt, multi_charm_ml::labelsCutScore}, "ML selections per pT bin"};
    Configurable<std::vector<int>> cutDir{"cutDir", std::vector<int>{o2::cuts_ml::CutDirection::CutNot}, "Whether to reject score values greater or smaller than the threshold"};
    Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"MCharm_BDTModel.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
    Configurable<std::vector<std::string>> namesInputFeatures{"namesInputFeatures", std::vector<std::string>(multi_charm_ml::namesInputFeatures), "Names of ML model input features"};
  } ml;

  ConfigurableAxis axisEta{"axisEta", {80, -4.0f, +4.0f}, "#eta"};
  ConfigurableAxis axisXicMass{"axisXicMass", {200, 2.368f, 2.568f}, "Xic Inv Mass (GeV/c^{2})"};
  ConfigurableAxis axisXiccMass{"axisXiccMass", {200, 3.521f, 3.721f}, "Xicc Inv Mass (GeV/c^{2})"};
  ConfigurableAxis axisDCA{"axisDCA", {400, 0, 400}, "DCA (#mum)"};
  ConfigurableAxis axisRadiusLarge{"axisRadiusLarge", {1000, 0, 20}, "Decay radius (cm)"};
  ConfigurableAxis axisRadius{"axisRadius", {10000, 0, 10000}, "Decay radius (#mum)"};
  ConfigurableAxis axisNSigma{"axisNSigma", {21, -10, 10}, "nsigma"};
  ConfigurableAxis axisDecayLength{"axisDecayLength", {2000, 0, 2000}, "Decay lenght (#mum)"};
  ConfigurableAxis axisDcaDaughters{"axisDcaDaughters", {200, 0, 100}, "DCA (mum)"};
  ConfigurableAxis axisBDTScore{"axisBDTScore", {100, 0, 1}, "BDT Score"};
  ConfigurableAxis axisBDTScoreFine{"axisBDTScoreFine", {1000, 0, 1}, "BDT Score for 1D histogram"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for QA histograms"};

  struct : ConfigurableGroup {
    std::string prefix = "selVals";
    Configurable<float> xiMinConstDCAxy{"xiMinConstDCAxy", -1, "[0] in |DCAxy| > [0]+[1]/pT"};
    Configurable<float> xiMinPtDepDCAxy{"xiMinPtDepDCAxy", 0, "[1] in |DCAxy| > [0]+[1]/pT"};
    Configurable<float> xiMinConstDCAz{"xiMinConstDCAz", -1, "[0] in |DCAz| > [0]+[1]/pT"};
    Configurable<float> xiMinPtDepDCAz{"xiMinPtDepDCAz", 0, "[1] in |DCAz| > [0]+[1]/pT"};
    Configurable<float> xiMinRadius{"xiMinRadius", -1, "Minimum R2D for Xi decay (cm)"};

    Configurable<float> picMinConstDCAxy{"picMinConstDCAxy", -1, "[0] in |DCAz| > [0]+[1]/pT"};
    Configurable<float> picMinPtDepDCAxy{"picMinPtDepDCAxy", 0, "[1] in |DCAz| > [0]+[1]/pT"};
    Configurable<float> picMinConstDCAz{"picMinConstDCAz", -1, "[0] in |DCAxy| > [0]+[1]/pT"};
    Configurable<float> picMinPtDepDCAz{"picMinPtDepDCAz", 0, "[1] in |DCAxy| > [0]+[1]/pT"};
    Configurable<float> picMinPt{"picMinPt", -1, "Minimum pT for Xic pions"};

    Configurable<float> piccMinConstDCAxy{"piccMinConstDCAxy", -1, "[0] in |DCAxy| > [0]+[1]/pT"};
    Configurable<float> piccMinPtDepDCAxy{"piccMinPtDepDCAxy", 0, "[1] in |DCAxy| > [0]+[1]/pT"};
    Configurable<float> piccMinConstDCAz{"piccMinConstDCAz", -1, "[0] in |DCAz| > [0]+[1]/pT"};
    Configurable<float> piccMinPtDepDCAz{"piccMinPtDepDCAz", 0, "[1] in |DCAz| > [0]+[1]/pT"};
    Configurable<float> piccMinPt{"piccMinPt", -1, "Minimum pT for Xicc pions"};

    Configurable<float> xicMinConstDCAxy{"xicMinConstDCAxy", -1, "[0] in |DCAz| > [0]+[1]/pT"};
    Configurable<float> xicMinConstDCAz{"xicMinConstDCAz", -1, "[0] in |DCAxy| > [0]+[1]/pT"};
    Configurable<float> xicMinPtDepDCAxy{"xicMinPtDepDCAxy", -1, "[0] in |DCAz| > [0]+[1]/pT"};
    Configurable<float> xicMinPtDepDCAz{"xicMinPtDepDCAz", -1, "[0] in |DCAxy| > [0]+[1]/pT"};
    Configurable<float> xicMaxDauDCA{"xicMaxDauDCA", 1e+4, "DCA between Xic daughters (cm)"};
    Configurable<float> xicMinRadius{"xicMinRadius", -1, "Minimum R2D for Xic decay (cm)"};
    Configurable<float> xicMinDecayDistanceFromPV{"xicMinDecayDistanceFromPV", -1, "Minimum distance for Xic decay from PV (cm)"};
    Configurable<float> xicMinProperLength{"xicMinProperLength", -1, "Minimum proper length for Xic decay (cm)"};
    Configurable<float> xicMaxProperLength{"xicMaxProperLength", 1e+4, "Minimum proper length for Xic decay (cm)"};

    Configurable<float> xiccMaxDCAxy{"xiccMaxDCAxy", 1e+4, "Maximum DCAxy"};
    Configurable<float> xiccMaxDCAz{"xiccMaxDCAz", 1e+4, "Maximum DCAz"};
    Configurable<float> xiccMaxDauDCA{"xiccMaxDauDCA", 1e+4, "DCA between Xicc daughters (cm)"};
    Configurable<float> xiccMinRadius{"xiccMinRadius", -1, "Minimum R2D for Xicc decay (cm)"};
    Configurable<float> xiccMinProperLength{"xiccMinProperLength", -1, "Minimum proper length for Xicc decay (cm)"};
    Configurable<float> xiccMaxProperLength{"xiccMaxProperLength", 1e+4, "Minimum proper length for Xicc decay (cm)"};
  } selVals;

  struct : ConfigurableGroup {
    std::string prefix = "selFlags";
    Configurable<bool> applyXiMinDCAxy{"applyXiMinDCAxy", false, "Apply |DCAxy| > [0]+[1]/pT"};
    Configurable<bool> applyXiMinDCAz{"applyXiMinDCAz", false, "Apply |DCAz| > [0]+[1]/pT"};
    Configurable<bool> applyXiMinRadius{"applyXiMinRadius", false, "Apply min radius"};

    Configurable<bool> applyPicMinDCAxy{"applyPicMinDCAxy", false, "Apply |DCAxy| > [0]+[1]/pT"};
    Configurable<bool> applyPicMinDCAz{"applyPicMinDCAz", false, "Apply |DCAz| > [0]+[1]/pT"};
    Configurable<bool> applyPiccMinDCAxy{"applyPiccMinDCAxy", false, "Apply |DCAxy| > [0]+[1]/pT"};
    Configurable<bool> applyPiccMinDCAz{"applyPiccMinDCAz", false, "Apply |DCAz| > [0]+[1]/pT"};

    Configurable<bool> applyXicMinDCAxy{"applyXicMinDCAxy", false, "Apply |DCAxy| > [0]+[1]/pT"};
    Configurable<bool> applyXicMinDCAz{"applyXicMinDCAz", false, "Apply |DCAz| > [0]+[1]/pT"};
    Configurable<bool> applyXicMinRadius{"applyXicMinRadius", false, "Apply min radius"};
    Configurable<bool> applyXicMaxDauDCA{"applyXicMaxDauDCA", false, "Apply max dau dca"};
    Configurable<bool> applyXicMinDistanceFromPV{"applyXicMinDistanceFromPV", false, "Apply min distance from PV (3D)"};
    Configurable<bool> applyXicMinProperLength{"applyXicMinProperLength", false, "Apply min proper length"};
    Configurable<bool> applyXicMaxProperLength{"applyXicMaxProperLength", false, "Apply max proper length"};

    Configurable<bool> applyXiccMinDCAxy{"applyXiccMinDCAxy", false, "Apply |DCAxy| > [0]+[1]/pT"};
    Configurable<bool> applyXiccMinDCAz{"applyXiccMinDCAz", false, "Apply |DCAz| > [0]+[1]/pT"};
    Configurable<bool> applyXiccMinRadius{"applyXiccMinRadius", false, "Apply min radius"};
    Configurable<bool> applyXiccMaxDauDCA{"applyXiccMaxDauDCA", false, "Apply max dau dca"};
    Configurable<bool> applyXiccMinProperLength{"applyXiccMinProperLength", false, "Apply min proper length"};
    Configurable<bool> applyXiccMaxProperLength{"applyXiccMaxProperLength", false, "Apply max proper length"};
  } selFlags;

  Configurable<int> otfConfig{"otfConfig", 0, "OTF configuration flag"};
  Filter configFilter = (aod::otfmulticharm::lutConfigId == otfConfig);

  void init(InitContext&)
  {
    histos.add("BeforeSel/hDCAXicDaughters", "hDCAXicDaughters; DCA between Xic daughters (#mum)", kTH1D, {axisDcaDaughters});
    histos.add("BeforeSel/hDCAXiccDaughters", "hDCAXiccDaughters; DCA between Xicc daughters (#mum)", kTH1D, {axisDcaDaughters});
    histos.add("BeforeSel/hDCAxyXi", "hDCAxyXi; Xi DCAxy to PV (#mum)", kTH1D, {axisDCA});
    histos.add("BeforeSel/hDCAzXi", "hDCAzXi; Xi DCAz to PV (#mum)", kTH1D, {axisDCA});
    histos.add("BeforeSel/hDCAxyXic", "hDCAxyXic; Xic DCAxy to PV (#mum)", kTH1D, {axisDCA});
    histos.add("BeforeSel/hDCAzXic", "hDCAzXic; Xic DCAz to PV (#mum)", kTH1D, {axisDCA});
    histos.add("BeforeSel/hDCAxyXicc", "hDCAxyXicc; Xicc DCAxy to PV (#mum)", kTH1D, {axisDCA});
    histos.add("BeforeSel/hDCAzXicc", "hDCAzXicc; Xicc DCAz to PV (#mum)", kTH1D, {axisDCA});
    histos.add("BeforeSel/hDecayRadiusXic", "hDecayRadiusXic; Distance (#mum)", kTH1D, {axisRadius});
    histos.add("BeforeSel/hDecayRadiusXicc", "hDecayRadiusXicc; Distance (#mum)", kTH1D, {axisRadius});
    histos.add("BeforeSel/hDecayDistanceFromPVXic", "hDecayDistanceFromPVXic; Distance (#mum)", kTH1D, {axisDecayLength});
    histos.add("BeforeSel/hProperLengthXic", "hProperLengthXic; Distance (#mum)", kTH1D, {axisDecayLength});
    histos.add("BeforeSel/hProperLengthXicc", "hProperLengthXicc; Distance (#mum)", kTH1D, {axisDecayLength});
    histos.add("BeforeSel/hPi1cDCAxy", "hPi1cDCAxy; Pi1c DCAxy (#mum)", kTH1D, {axisDCA});
    histos.add("BeforeSel/hPi1cDCAz", "hPi1cDCAz; Pi1c DCAz (#mum)", kTH1D, {axisDCA});
    histos.add("BeforeSel/hPi2cDCAxy", "hPi2cDCAxy; Pi2c DCAxy (#mum)", kTH1D, {axisDCA});
    histos.add("BeforeSel/hPi2cDCAz", "hPi2cDCAz; Pi2c DCAz (#mum)", kTH1D, {axisDCA});
    histos.add("BeforeSel/hPiccDCAxy", "hPiccDCAxy; Picc DCAxy (#mum)", kTH1D, {axisDCA});
    histos.add("BeforeSel/hPiccDCAz", "hPiccDCAz; Picc DCAz (#mum)", kTH1D, {axisDCA});
    histos.add("BeforeSel/hPi1cPt", "hPi1cPt; Pi1c pT (Gev/#it(c))", kTH1D, {axisPt});
    histos.add("BeforeSel/hPi2cPt", "hPi2cPt; Pi2c pT (Gev/#it(c))", kTH1D, {axisPt});
    histos.add("BeforeSel/hPiccPt", "hPiccPt; Picc pT (Gev/#it(c))", kTH1D, {axisPt});
    histos.addClone("BeforeSel/", "AfterSel/");

    auto hMCharmBuilding = histos.add<TH1>("hMCharmBuilding", "hMCharmBuilding", kTH1D, {{22, -0.5, 21.5}});
    hMCharmBuilding->GetXaxis()->SetBinLabel(1, "nTotalCandidates");
    hMCharmBuilding->GetXaxis()->SetBinLabel(2, "xicMaxDauDCA");
    hMCharmBuilding->GetXaxis()->SetBinLabel(3, "xiccMaxDauDCA");
    hMCharmBuilding->GetXaxis()->SetBinLabel(4, "xiMinConstDCAxy");
    hMCharmBuilding->GetXaxis()->SetBinLabel(5, "xiMinConstDCAz");
    hMCharmBuilding->GetXaxis()->SetBinLabel(6, "pi1cMinDCAxy");
    hMCharmBuilding->GetXaxis()->SetBinLabel(7, "pi1cMinDCAz");
    hMCharmBuilding->GetXaxis()->SetBinLabel(8, "pi2cMinDCAxy");
    hMCharmBuilding->GetXaxis()->SetBinLabel(9, "pi2cMinDCAz");
    hMCharmBuilding->GetXaxis()->SetBinLabel(10, "piccMinDCAxy");
    hMCharmBuilding->GetXaxis()->SetBinLabel(11, "piccMinDCAz");
    hMCharmBuilding->GetXaxis()->SetBinLabel(12, "xicMinDCAxy");
    hMCharmBuilding->GetXaxis()->SetBinLabel(13, "xicMinDCAz");
    hMCharmBuilding->GetXaxis()->SetBinLabel(14, "xiccMaxDCAxy");
    hMCharmBuilding->GetXaxis()->SetBinLabel(15, "xiccMaxDCAz");
    hMCharmBuilding->GetXaxis()->SetBinLabel(16, "xicMinRadius");
    hMCharmBuilding->GetXaxis()->SetBinLabel(17, "xiccMinRadius");
    hMCharmBuilding->GetXaxis()->SetBinLabel(18, "xicMinProperLength");
    hMCharmBuilding->GetXaxis()->SetBinLabel(19, "xicMaxProperLength");
    hMCharmBuilding->GetXaxis()->SetBinLabel(20, "xiccMinProperLength");
    hMCharmBuilding->GetXaxis()->SetBinLabel(21, "xiccMaxProperLength");
    hMCharmBuilding->GetXaxis()->SetBinLabel(22, "xicMinDecayDistanceFromPV");

    histos.add("hXiccMass", "hXiccMass", kTH1D, {axisXiccMass});
    histos.add("hXicMass", "hXicMass", kTH1D, {axisXicMass});
    histos.add("hXiccPt", "hXiccPt", kTH1D, {axisPt});
    histos.add("hXicPt", "hXicPt", kTH1D, {axisPt});
    histos.add("h3dXicc", "h3dXicc; Xicc pT (GeV/#it(c)); Xicc #eta; Xicc mass (GeV/#it(c)^{2})", kTH3D, {axisPt, axisEta, axisXiccMass});
    histos.add("hConfigId", "hConfigId", kTH1D, {{11, -0.5, 10.5}});

    if (ml.enableBDT) {
      mlResponse.configure(ml.ptBinEdges, ml.scoreCuts, ml.cutDir, multi_charm_ml::NClasses);
      if (ml.loadModelsFromCCDB) {
        ccdbApi.init(ml.ccdbUrl);
        mlResponse.setModelPathsCCDB(ml.onnxFileNames, ccdbApi, ml.modelPathsCCDB, ml.timestampCCDB);
      } else {
        mlResponse.setModelPathsLocal(ml.onnxFileNames);
      }

      mlResponse.cacheInputFeaturesIndices(ml.namesInputFeatures);
      mlResponse.init();

      histos.add("BDT/BeforeSel/hBDTScoreSignalFine", "hBDTScoreSignalFine", kTH1D, {axisBDTScoreFine});
      histos.add("BDT/BeforeSel/hBDTScoreSignal", "hBDTScoreSignal", kTH1D, {axisBDTScore});
      histos.add("BDT/BeforeSel/hBDTScoreVsXiccMassSignal", "hBDTScoreVsXiccMassSignal", kTH2D, {axisXiccMass, axisBDTScore});
      histos.add("BDT/BeforeSel/hBDTScoreVsXiccPtSignal", "hBDTScoreVsXiccPtSignal", kTH2D, {axisPt, axisBDTScore});
      histos.add("BDT/BeforeSel/h3dBDTScoreSignal", "h3dBDTScoreSignal", kTH3D, {axisPt, axisXiccMass, axisBDTScore});
      histos.add("BDT/BeforeSel/hBDTScoreBackgroundFine", "hBDTScoreBackgroundFine", kTH1D, {axisBDTScoreFine});
      histos.add("BDT/BeforeSel/hBDTScoreBackground", "hBDTScoreBackground", kTH1D, {axisBDTScore});
      histos.add("BDT/BeforeSel/hBDTScoreVsXiccMassBackground", "hBDTScoreVsXiccMassBackground", kTH2D, {axisXiccMass, axisBDTScore});
      histos.add("BDT/BeforeSel/hBDTScoreVsXiccPtBackground", "hBDTScoreVsXiccPtBackground", kTH2D, {axisPt, axisBDTScore});
      histos.add("BDT/BeforeSel/h3dBDTScoreBackground", "h3dBDTScoreBackground", kTH3D, {axisPt, axisXiccMass, axisBDTScore});
      histos.add("BDT/BeforeSel/hDCAXicDaughters", "hDCAXicDaughters", kTH2D, {{axisBDTScore, axisDcaDaughters}});
      histos.add("BDT/BeforeSel/hDCAXiccDaughters", "hDCAXiccDaughters", kTH2D, {{axisBDTScore, axisDcaDaughters}});
      histos.add("BDT/BeforeSel/hDCAxyXi", "hDCAxyXi", kTH2D, {{axisBDTScore, axisDCA}});
      histos.add("BDT/BeforeSel/hDCAzXi", "hDCAzXi", kTH2D, {{axisBDTScore, axisDCA}});
      histos.add("BDT/BeforeSel/hDCAxyXic", "hDCAxyXic", kTH2D, {{axisBDTScore, axisDCA}});
      histos.add("BDT/BeforeSel/hDCAzXic", "hDCAzXic", kTH2D, {{axisBDTScore, axisDCA}});
      histos.add("BDT/BeforeSel/hDCAxyXicc", "hDCAxyXicc", kTH2D, {{axisBDTScore, axisDCA}});
      histos.add("BDT/BeforeSel/hDCAzXicc", "hDCAzXicc", kTH2D, {{axisBDTScore, axisDCA}});
      histos.add("BDT/BeforeSel/hDecayRadiusXic", "hDecayRadiusXic", kTH2D, {{axisBDTScore, axisRadius}});
      histos.add("BDT/BeforeSel/hDecayRadiusXicc", "hDecayRadiusXicc", kTH2D, {{axisBDTScore, axisRadius}});
      histos.add("BDT/BeforeSel/hDecayDistanceFromPVXic", "hDecayDistanceFromPVXic", kTH2D, {{axisBDTScore, axisDecayLength}});
      histos.add("BDT/BeforeSel/hProperLengthXic", "hProperLengthXic", kTH2D, {{axisBDTScore, axisDecayLength}});
      histos.add("BDT/BeforeSel/hProperLengthXicc", "hProperLengthXicc", kTH2D, {{axisBDTScore, axisDecayLength}});
      histos.add("BDT/BeforeSel/hPi1cDCAxy", "hPi1cDCAxy", kTH2D, {{axisBDTScore, axisDCA}});
      histos.add("BDT/BeforeSel/hPi1cDCAz", "hPi1cDCAz", kTH2D, {{axisBDTScore, axisDCA}});
      histos.add("BDT/BeforeSel/hPi2cDCAxy", "hPi2cDCAxy", kTH2D, {{axisBDTScore, axisDCA}});
      histos.add("BDT/BeforeSel/hPi2cDCAz", "hPi2cDCAz", kTH2D, {{axisBDTScore, axisDCA}});
      histos.add("BDT/BeforeSel/hPiccDCAxy", "hPiccDCAxy", kTH2D, {{axisBDTScore, axisDCA}});
      histos.add("BDT/BeforeSel/hPiccDCAz", "hPiccDCAz", kTH2D, {{axisBDTScore, axisDCA}});
      histos.add("BDT/BeforeSel/hPi1cPt", "hPi1cPt", kTH2D, {{axisBDTScore, axisPt}});
      histos.add("BDT/BeforeSel/hPi2cPt", "hPi2cPt", kTH2D, {{axisBDTScore, axisPt}});
      histos.add("BDT/BeforeSel/hPiccPt", "hPiccPt", kTH2D, {{axisBDTScore, axisPt}});
      histos.addClone("BDT/BeforeSel/", "BDT/AfterSel/");
    }

    histos.print();
  }

  void processXicc(soa::Filtered<aod::MCharmCores> const& xiccCands)
  {
    for (const auto& xiccCand : xiccCands) {
      float bdtPredictedBackground = -999.f;
      float bdtPredictedSignal = -999.f;

      histos.fill(HIST("hConfigId"), xiccCand.lutConfigId());
      histos.fill(HIST("BeforeSel/hDCAXicDaughters"), xiccCand.xicDauDCA() * ToMicrons);
      histos.fill(HIST("BeforeSel/hDCAXiccDaughters"), xiccCand.xiccDauDCA() * ToMicrons);
      histos.fill(HIST("BeforeSel/hDCAxyXi"), std::fabs(xiccCand.xiDCAxy() * ToMicrons));
      histos.fill(HIST("BeforeSel/hDCAzXi"), std::fabs(xiccCand.xiDCAz() * ToMicrons));
      histos.fill(HIST("BeforeSel/hDCAxyXic"), std::fabs(xiccCand.xicDCAxy() * ToMicrons));
      histos.fill(HIST("BeforeSel/hDCAzXic"), std::fabs(xiccCand.xicDCAz() * ToMicrons));
      histos.fill(HIST("BeforeSel/hDCAxyXicc"), std::fabs(xiccCand.xiccDCAxy() * ToMicrons));
      histos.fill(HIST("BeforeSel/hDCAzXicc"), std::fabs(xiccCand.xiccDCAz() * ToMicrons));
      histos.fill(HIST("BeforeSel/hDecayRadiusXic"), xiccCand.xicDecayRadius2D() * ToMicrons);
      histos.fill(HIST("BeforeSel/hDecayRadiusXicc"), xiccCand.xiccDecayRadius2D() * ToMicrons);
      histos.fill(HIST("BeforeSel/hDecayDistanceFromPVXic"), xiccCand.xicDistanceFromPV() * ToMicrons);
      histos.fill(HIST("BeforeSel/hProperLengthXic"), xiccCand.xicProperLength() * ToMicrons);
      histos.fill(HIST("BeforeSel/hProperLengthXicc"), xiccCand.xiccProperLength() * ToMicrons);
      histos.fill(HIST("BeforeSel/hPi1cDCAxy"), xiccCand.pi1cDCAxy() * ToMicrons);
      histos.fill(HIST("BeforeSel/hPi1cDCAz"), xiccCand.pi1cDCAz() * ToMicrons);
      histos.fill(HIST("BeforeSel/hPi2cDCAxy"), xiccCand.pi2cDCAxy() * ToMicrons);
      histos.fill(HIST("BeforeSel/hPi2cDCAz"), xiccCand.pi2cDCAz() * ToMicrons);
      histos.fill(HIST("BeforeSel/hPiccDCAxy"), xiccCand.piccDCAxy() * ToMicrons);
      histos.fill(HIST("BeforeSel/hPiccDCAz"), xiccCand.piccDCAz() * ToMicrons);
      histos.fill(HIST("BeforeSel/hPi1cPt"), xiccCand.pi1cPt());
      histos.fill(HIST("BeforeSel/hPi2cPt"), xiccCand.pi2cPt());
      histos.fill(HIST("BeforeSel/hPiccPt"), xiccCand.piccPt());

      if (ml.enableBDT) {
        std::vector<float> outputScore;
        std::vector<float> inputFeatures{
          xiccCand.xicDauDCA(),
          xiccCand.xiccDauDCA(),
          xiccCand.xiDCAxy(),
          xiccCand.xicDCAxy(),
          xiccCand.xiccDCAxy(),
          xiccCand.xiDCAz(),
          xiccCand.xicDCAz(),
          xiccCand.xiccDCAz(),
          xiccCand.pi1cDCAxy(),
          xiccCand.pi2cDCAxy(),
          xiccCand.piccDCAxy(),
          xiccCand.pi1cDCAz(),
          xiccCand.pi2cDCAz(),
          xiccCand.piccDCAz(),
          xiccCand.xicDecayRadius2D(),
          xiccCand.xiccDecayRadius2D(),
          xiccCand.xicProperLength(),
          xiccCand.xicDistanceFromPV(),
          xiccCand.xiccProperLength()};

        mlResponse.isSelectedMl(inputFeatures, xiccCand.xiccPt(), outputScore);
        bdtPredictedBackground = static_cast<float>(outputScore[0]);
        bdtPredictedSignal = static_cast<float>(outputScore[1]);

        histos.fill(HIST("BDT/BeforeSel/hBDTScoreSignal"), bdtPredictedSignal);
        histos.fill(HIST("BDT/BeforeSel/hBDTScoreSignalFine"), bdtPredictedSignal);
        histos.fill(HIST("BDT/BeforeSel/hBDTScoreVsXiccMassSignal"), xiccCand.xiccMass(), bdtPredictedSignal);
        histos.fill(HIST("BDT/BeforeSel/hBDTScoreVsXiccPtSignal"), xiccCand.xiccPt(), bdtPredictedSignal);
        histos.fill(HIST("BDT/BeforeSel/h3dBDTScoreSignal"), xiccCand.xiccPt(), xiccCand.xiccMass(), bdtPredictedSignal);
        histos.fill(HIST("BDT/BeforeSel/hBDTScoreBackground"), bdtPredictedBackground);
        histos.fill(HIST("BDT/BeforeSel/hBDTScoreBackgroundFine"), bdtPredictedBackground);
        histos.fill(HIST("BDT/BeforeSel/hBDTScoreVsXiccMassBackground"), xiccCand.xiccMass(), bdtPredictedBackground);
        histos.fill(HIST("BDT/BeforeSel/hBDTScoreVsXiccPtBackground"), xiccCand.xiccPt(), bdtPredictedBackground);
        histos.fill(HIST("BDT/BeforeSel/h3dBDTScoreBackground"), xiccCand.xiccPt(), xiccCand.xiccMass(), bdtPredictedBackground);
        histos.fill(HIST("BDT/BeforeSel/hDCAXicDaughters"), bdtPredictedSignal, xiccCand.xicDauDCA() * ToMicrons);
        histos.fill(HIST("BDT/BeforeSel/hDCAXiccDaughters"), bdtPredictedSignal, xiccCand.xiccDauDCA() * ToMicrons);
        histos.fill(HIST("BDT/BeforeSel/hDCAxyXi"), bdtPredictedSignal, std::fabs(xiccCand.xiDCAxy() * ToMicrons));
        histos.fill(HIST("BDT/BeforeSel/hDCAzXi"), bdtPredictedSignal, std::fabs(xiccCand.xiDCAz() * ToMicrons));
        histos.fill(HIST("BDT/BeforeSel/hDCAxyXic"), bdtPredictedSignal, std::fabs(xiccCand.xicDCAxy() * ToMicrons));
        histos.fill(HIST("BDT/BeforeSel/hDCAzXic"), bdtPredictedSignal, std::fabs(xiccCand.xicDCAz() * ToMicrons));
        histos.fill(HIST("BDT/BeforeSel/hDCAxyXicc"), bdtPredictedSignal, std::fabs(xiccCand.xiccDCAxy() * ToMicrons));
        histos.fill(HIST("BDT/BeforeSel/hDCAzXicc"), bdtPredictedSignal, std::fabs(xiccCand.xiccDCAz() * ToMicrons));
        histos.fill(HIST("BDT/BeforeSel/hDecayRadiusXic"), bdtPredictedSignal, xiccCand.xicDecayRadius2D() * ToMicrons);
        histos.fill(HIST("BDT/BeforeSel/hDecayRadiusXicc"), bdtPredictedSignal, xiccCand.xiccDecayRadius2D() * ToMicrons);
        histos.fill(HIST("BDT/BeforeSel/hDecayDistanceFromPVXic"), bdtPredictedSignal, xiccCand.xicDistanceFromPV() * ToMicrons);
        histos.fill(HIST("BDT/BeforeSel/hProperLengthXic"), bdtPredictedSignal, xiccCand.xicProperLength() * ToMicrons);
        histos.fill(HIST("BDT/BeforeSel/hProperLengthXicc"), bdtPredictedSignal, xiccCand.xiccProperLength() * ToMicrons);
        histos.fill(HIST("BDT/BeforeSel/hPi1cDCAxy"), bdtPredictedSignal, xiccCand.pi1cDCAxy() * ToMicrons);
        histos.fill(HIST("BDT/BeforeSel/hPi1cDCAz"), bdtPredictedSignal, xiccCand.pi1cDCAz() * ToMicrons);
        histos.fill(HIST("BDT/BeforeSel/hPi2cDCAxy"), bdtPredictedSignal, xiccCand.pi2cDCAxy() * ToMicrons);
        histos.fill(HIST("BDT/BeforeSel/hPi2cDCAz"), bdtPredictedSignal, xiccCand.pi2cDCAz() * ToMicrons);
        histos.fill(HIST("BDT/BeforeSel/hPiccDCAxy"), bdtPredictedSignal, xiccCand.piccDCAxy() * ToMicrons);
        histos.fill(HIST("BDT/BeforeSel/hPiccDCAz"), bdtPredictedSignal, xiccCand.piccDCAz() * ToMicrons);
        histos.fill(HIST("BDT/BeforeSel/hPi1cPt"), bdtPredictedSignal, xiccCand.pi1cPt());
        histos.fill(HIST("BDT/BeforeSel/hPi2cPt"), bdtPredictedSignal, xiccCand.pi2cPt());
        histos.fill(HIST("BDT/BeforeSel/hPiccPt"), bdtPredictedSignal, xiccCand.piccPt());
      }

      histos.fill(HIST("hMCharmBuilding"), 0);
      if (selFlags.applyXicMaxDauDCA && xiccCand.xicDauDCA() > selVals.xicMaxDauDCA) {
        continue;
      }

      histos.fill(HIST("hMCharmBuilding"), 1);
      if (selFlags.applyXiccMaxDauDCA && xiccCand.xiccDauDCA() > selVals.xiccMaxDauDCA) {
        continue;
      }

      histos.fill(HIST("hMCharmBuilding"), 2);
      if (selFlags.applyXiMinDCAxy && std::fabs(xiccCand.xiDCAxy()) < selVals.xiMinConstDCAxy + (selVals.xiMinPtDepDCAxy / xiccCand.xiPt())) {
        continue;
      }

      histos.fill(HIST("hMCharmBuilding"), 3);
      if (selFlags.applyXiMinDCAz && std::fabs(xiccCand.xiDCAz()) < selVals.xiMinConstDCAz + (selVals.xiMinPtDepDCAz / xiccCand.xiPt())) {
        continue;
      }

      histos.fill(HIST("hMCharmBuilding"), 4);
      if (selFlags.applyPicMinDCAxy && std::fabs(xiccCand.pi1cDCAxy()) < selVals.picMinConstDCAxy + (selVals.picMinPtDepDCAxy / xiccCand.pi1cPt())) {
        continue;
      }

      histos.fill(HIST("hMCharmBuilding"), 5);
      if (selFlags.applyPicMinDCAz && std::fabs(xiccCand.pi1cDCAz()) < selVals.picMinConstDCAz + (selVals.picMinPtDepDCAz / xiccCand.pi1cPt())) {
        continue;
      }

      histos.fill(HIST("hMCharmBuilding"), 6);
      if (selFlags.applyPicMinDCAxy && std::fabs(xiccCand.pi2cDCAxy()) < selVals.picMinConstDCAxy + (selVals.picMinPtDepDCAxy / xiccCand.pi2cPt())) {
        continue;
      }

      histos.fill(HIST("hMCharmBuilding"), 7);
      if (selFlags.applyPicMinDCAz && std::fabs(xiccCand.pi2cDCAz()) < selVals.picMinConstDCAz + (selVals.picMinPtDepDCAz / xiccCand.pi2cPt())) {
        continue;
      }

      histos.fill(HIST("hMCharmBuilding"), 8);
      if (selFlags.applyPiccMinDCAxy && std::fabs(xiccCand.piccDCAxy()) < selVals.piccMinConstDCAxy + (selVals.piccMinPtDepDCAxy / xiccCand.piccPt())) {
        continue;
      }

      histos.fill(HIST("hMCharmBuilding"), 9);
      if (selFlags.applyPiccMinDCAz && std::fabs(xiccCand.piccDCAz()) < selVals.piccMinConstDCAz + (selVals.piccMinPtDepDCAz / xiccCand.piccPt())) {
        continue;
      }

      histos.fill(HIST("hMCharmBuilding"), 10);
      if (selFlags.applyXicMinDCAxy && std::fabs(xiccCand.xicDCAxy()) < selVals.xicMinConstDCAxy + (selVals.xicMinPtDepDCAxy / xiccCand.xicPt())) {
        continue;
      }

      histos.fill(HIST("hMCharmBuilding"), 11);
      if (selFlags.applyXicMinDCAz && std::fabs(xiccCand.xicDCAz()) < selVals.xicMinConstDCAz + (selVals.xicMinPtDepDCAz / xiccCand.xicPt())) {
        continue;
      }

      histos.fill(HIST("hMCharmBuilding"), 12);
      if (selFlags.applyXiccMinDCAxy && std::fabs(xiccCand.xiccDCAxy()) > selVals.xiccMaxDCAxy) {
        continue;
      }

      histos.fill(HIST("hMCharmBuilding"), 13);
      if (selFlags.applyXiccMinDCAz && std::fabs(xiccCand.xiccDCAz()) > selVals.xiccMaxDCAz) {
        continue;
      }

      histos.fill(HIST("hMCharmBuilding"), 14);
      if (selFlags.applyXicMinRadius && xiccCand.xicDecayRadius2D() < selVals.xicMinRadius) {
        continue;
      }

      histos.fill(HIST("hMCharmBuilding"), 15);
      if (selFlags.applyXiccMinRadius && xiccCand.xiccDecayRadius2D() < selVals.xiccMinRadius) {
        continue;
      }

      histos.fill(HIST("hMCharmBuilding"), 16);
      if (selFlags.applyXicMinProperLength && xiccCand.xicProperLength() < selVals.xicMinProperLength) {
        continue;
      }

      histos.fill(HIST("hMCharmBuilding"), 17);
      if (selFlags.applyXicMaxProperLength && xiccCand.xicProperLength() > selVals.xicMaxProperLength) {
        continue;
      }

      histos.fill(HIST("hMCharmBuilding"), 18);
      if (selFlags.applyXiccMinProperLength && xiccCand.xiccProperLength() < selVals.xiccMinProperLength) {
        continue;
      }

      histos.fill(HIST("hMCharmBuilding"), 19);
      if (selFlags.applyXiccMaxProperLength && xiccCand.xiccProperLength() > selVals.xiccMaxProperLength) {
        continue;
      }

      histos.fill(HIST("hMCharmBuilding"), 20);
      if (selFlags.applyXicMinDistanceFromPV && xiccCand.xicDistanceFromPV() < selVals.xicMinDecayDistanceFromPV) {
        continue;
      }

      histos.fill(HIST("hMCharmBuilding"), 21);
      histos.fill(HIST("AfterSel/hDCAXicDaughters"), xiccCand.xicDauDCA() * ToMicrons);
      histos.fill(HIST("AfterSel/hDCAXiccDaughters"), xiccCand.xiccDauDCA() * ToMicrons);
      histos.fill(HIST("AfterSel/hDCAxyXi"), std::fabs(xiccCand.xiDCAxy() * ToMicrons));
      histos.fill(HIST("AfterSel/hDCAzXi"), std::fabs(xiccCand.xiDCAz() * ToMicrons));
      histos.fill(HIST("AfterSel/hDCAxyXic"), std::fabs(xiccCand.xicDCAxy() * ToMicrons));
      histos.fill(HIST("AfterSel/hDCAzXic"), std::fabs(xiccCand.xicDCAz() * ToMicrons));
      histos.fill(HIST("AfterSel/hDCAxyXicc"), std::fabs(xiccCand.xiccDCAxy() * ToMicrons));
      histos.fill(HIST("AfterSel/hDCAzXicc"), std::fabs(xiccCand.xiccDCAz() * ToMicrons));
      histos.fill(HIST("AfterSel/hDecayRadiusXic"), xiccCand.xicDecayRadius2D() * ToMicrons);
      histos.fill(HIST("AfterSel/hDecayRadiusXicc"), xiccCand.xiccDecayRadius2D() * ToMicrons);
      histos.fill(HIST("AfterSel/hDecayDistanceFromPVXic"), xiccCand.xicDistanceFromPV() * ToMicrons);
      histos.fill(HIST("AfterSel/hProperLengthXic"), xiccCand.xicProperLength() * ToMicrons);
      histos.fill(HIST("AfterSel/hProperLengthXicc"), xiccCand.xiccProperLength() * ToMicrons);
      histos.fill(HIST("AfterSel/hPi1cDCAxy"), xiccCand.pi1cDCAxy() * ToMicrons);
      histos.fill(HIST("AfterSel/hPi1cDCAz"), xiccCand.pi1cDCAz() * ToMicrons);
      histos.fill(HIST("AfterSel/hPi2cDCAxy"), xiccCand.pi2cDCAxy() * ToMicrons);
      histos.fill(HIST("AfterSel/hPi2cDCAz"), xiccCand.pi2cDCAz() * ToMicrons);
      histos.fill(HIST("AfterSel/hPiccDCAxy"), xiccCand.piccDCAxy() * ToMicrons);
      histos.fill(HIST("AfterSel/hPiccDCAz"), xiccCand.piccDCAz() * ToMicrons);
      histos.fill(HIST("AfterSel/hPi1cPt"), xiccCand.pi1cPt());
      histos.fill(HIST("AfterSel/hPi2cPt"), xiccCand.pi2cPt());
      histos.fill(HIST("AfterSel/hPiccPt"), xiccCand.piccPt());
      histos.fill(HIST("hXiccMass"), xiccCand.xiccMass());
      histos.fill(HIST("hXicMass"), xiccCand.xicMass());
      histos.fill(HIST("hXiccPt"), xiccCand.xiccPt());
      histos.fill(HIST("hXicPt"), xiccCand.xicPt());
      histos.fill(HIST("h3dXicc"), xiccCand.xiccPt(), xiccCand.xiccEta(), xiccCand.xiccMass());
      if (ml.enableBDT) {
        histos.fill(HIST("BDT/AfterSel/hBDTScoreSignal"), bdtPredictedSignal);
        histos.fill(HIST("BDT/AfterSel/hBDTScoreSignalFine"), bdtPredictedSignal);
        histos.fill(HIST("BDT/AfterSel/hBDTScoreVsXiccMassSignal"), xiccCand.xiccMass(), bdtPredictedSignal);
        histos.fill(HIST("BDT/AfterSel/hBDTScoreVsXiccPtSignal"), xiccCand.xiccPt(), bdtPredictedSignal);
        histos.fill(HIST("BDT/AfterSel/h3dBDTScoreSignal"), xiccCand.xiccPt(), xiccCand.xiccMass(), bdtPredictedSignal);
        histos.fill(HIST("BDT/AfterSel/hBDTScoreBackground"), bdtPredictedBackground);
        histos.fill(HIST("BDT/AfterSel/hBDTScoreBackgroundFine"), bdtPredictedBackground);
        histos.fill(HIST("BDT/AfterSel/hBDTScoreVsXiccMassBackground"), xiccCand.xiccMass(), bdtPredictedBackground);
        histos.fill(HIST("BDT/AfterSel/hBDTScoreVsXiccPtBackground"), xiccCand.xiccPt(), bdtPredictedBackground);
        histos.fill(HIST("BDT/AfterSel/h3dBDTScoreBackground"), xiccCand.xiccPt(), xiccCand.xiccMass(), bdtPredictedBackground);
        histos.fill(HIST("BDT/AfterSel/hDCAXicDaughters"), bdtPredictedSignal, xiccCand.xicDauDCA() * ToMicrons);
        histos.fill(HIST("BDT/AfterSel/hDCAXiccDaughters"), bdtPredictedSignal, xiccCand.xiccDauDCA() * ToMicrons);
        histos.fill(HIST("BDT/AfterSel/hDCAxyXi"), bdtPredictedSignal, std::fabs(xiccCand.xiDCAxy() * ToMicrons));
        histos.fill(HIST("BDT/AfterSel/hDCAzXi"), bdtPredictedSignal, std::fabs(xiccCand.xiDCAz() * ToMicrons));
        histos.fill(HIST("BDT/AfterSel/hDCAxyXic"), bdtPredictedSignal, std::fabs(xiccCand.xicDCAxy() * ToMicrons));
        histos.fill(HIST("BDT/AfterSel/hDCAzXic"), bdtPredictedSignal, std::fabs(xiccCand.xicDCAz() * ToMicrons));
        histos.fill(HIST("BDT/AfterSel/hDCAxyXicc"), bdtPredictedSignal, std::fabs(xiccCand.xiccDCAxy() * ToMicrons));
        histos.fill(HIST("BDT/AfterSel/hDCAzXicc"), bdtPredictedSignal, std::fabs(xiccCand.xiccDCAz() * ToMicrons));
        histos.fill(HIST("BDT/AfterSel/hDecayRadiusXic"), bdtPredictedSignal, xiccCand.xicDecayRadius2D() * ToMicrons);
        histos.fill(HIST("BDT/AfterSel/hDecayRadiusXicc"), bdtPredictedSignal, xiccCand.xiccDecayRadius2D() * ToMicrons);
        histos.fill(HIST("BDT/AfterSel/hDecayDistanceFromPVXic"), bdtPredictedSignal, xiccCand.xicDistanceFromPV() * ToMicrons);
        histos.fill(HIST("BDT/AfterSel/hProperLengthXic"), bdtPredictedSignal, xiccCand.xicProperLength() * ToMicrons);
        histos.fill(HIST("BDT/AfterSel/hProperLengthXicc"), bdtPredictedSignal, xiccCand.xiccProperLength() * ToMicrons);
        histos.fill(HIST("BDT/AfterSel/hPi1cDCAxy"), bdtPredictedSignal, xiccCand.pi1cDCAxy() * ToMicrons);
        histos.fill(HIST("BDT/AfterSel/hPi1cDCAz"), bdtPredictedSignal, xiccCand.pi1cDCAz() * ToMicrons);
        histos.fill(HIST("BDT/AfterSel/hPi2cDCAxy"), bdtPredictedSignal, xiccCand.pi2cDCAxy() * ToMicrons);
        histos.fill(HIST("BDT/AfterSel/hPi2cDCAz"), bdtPredictedSignal, xiccCand.pi2cDCAz() * ToMicrons);
        histos.fill(HIST("BDT/AfterSel/hPiccDCAxy"), bdtPredictedSignal, xiccCand.piccDCAxy() * ToMicrons);
        histos.fill(HIST("BDT/AfterSel/hPiccDCAz"), bdtPredictedSignal, xiccCand.piccDCAz() * ToMicrons);
        histos.fill(HIST("BDT/AfterSel/hPi1cPt"), bdtPredictedSignal, xiccCand.pi1cPt());
        histos.fill(HIST("BDT/AfterSel/hPi2cPt"), bdtPredictedSignal, xiccCand.pi2cPt());
        histos.fill(HIST("BDT/AfterSel/hPiccPt"), bdtPredictedSignal, xiccCand.piccPt());
      }
    }
  }

  PROCESS_SWITCH(Alice3Multicharm, processXicc, "find Xicc baryons", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<Alice3Multicharm>(cfgc)};
}
