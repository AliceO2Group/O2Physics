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

#include "ALICE3/DataModel/A3DecayFinderTables.h"
#include "ALICE3/DataModel/OTFMulticharm.h"
#include "ALICE3/DataModel/OTFRICH.h"
#include "ALICE3/DataModel/OTFStrangeness.h"
#include "ALICE3/DataModel/OTFTOF.h"
#include "ALICE3/DataModel/tracksAlice3.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Tools/ML/MlResponse.h"
#include "Tools/ML/model.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DCAFitter/DCAFitterN.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsVertexing/PVertexer.h"
#include "DetectorsVertexing/PVertexerHelpers.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <iterator>
#include <map>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::ml;
using namespace o2::framework;
using namespace o2::framework::expressions;

using MultiCharmTracksPID = soa::Join<aod::MCharmCores, aod::MCharmIndices>;
using MultiCharmTracksFull = soa::Join<aod::MCharmCores, aod::MCharmIndices, aod::MCharmExtra>;

struct Alice3Multicharm {
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  std::map<int, int> pdgToBin;
  o2::ml::OnnxModel bdtMCharm;
  std::map<std::string, std::string> metadata;
  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  struct : ConfigurableGroup {
    std::string prefix = "bdt"; // JSON group name
    Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> localPath{"localPath", "MCharm_BDTModel.onnx", "(std::string) Path to the local .onnx file."};
    Configurable<std::string> pathCCDB{"pathCCDB", "Users/j/jekarlss/MLModels", "Path on CCDB"};
    Configurable<int64_t> timestampCCDB{"timestampCCDB", 1695750420200, "timestamp of the ONNX file for ML model used to query in CCDB. Please use 1695750420200"};
    Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};
    Configurable<bool> enableOptimizations{"enableOptimizations", false, "Enables the ONNX extended model-optimization: sessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_EXTENDED)"};
    Configurable<bool> enableML{"enableML", false, "Enables bdt model"};
  } bdt;

  ConfigurableAxis axisEta{"axisEta", {80, -4.0f, +4.0f}, "#eta"};
  ConfigurableAxis axisXicMass{"axisXicMass", {200, 2.368f, 2.568f}, "Xic Inv Mass (GeV/c^{2})"};
  ConfigurableAxis axisXiccMass{"axisXiccMass", {200, 3.521f, 3.721f}, "Xicc Inv Mass (GeV/c^{2})"};
  ConfigurableAxis axisDCA{"axisDCA", {400, 0, 400}, "DCA (#mum)"};
  ConfigurableAxis axisRadiusLarge{"axisRadiusLarge", {1000, 0, 20}, "Decay radius (cm)"};
  ConfigurableAxis axisRadius{"axisRadius", {10000, 0, 10000}, "Decay radius (#mum)"};
  ConfigurableAxis axisTofTrackDelta{"axisTofTrackDelta", {200, 0, 1000}, "TOF track time"};
  ConfigurableAxis axisNSigma{"axisNSigma", {21, -10, 10}, "nsigma"};
  ConfigurableAxis axisDecayLength{"axisDecayLength", {2000, 0, 2000}, "Decay lenght (#mum)"};
  ConfigurableAxis axisDcaDaughters{"axisDcaDaughters", {200, 0, 100}, "DCA (mum)"};
  ConfigurableAxis axisBDTScore{"axisBDTScore", {100, 0, 1}, "BDT Score"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for QA histograms"};

  Configurable<float> xiMinDCAxy{"xiMinDCAxy", -1, "[0] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> xiMinDCAz{"xiMinDCAz", -1, "[0] in |DCAz| > [0]+[1]/pT"};
  Configurable<float> xiMinRadius{"xiMinRadius", -1, "Minimum R2D for Xic decay (cm)"};

  Configurable<float> picMinDCAxy{"picMinDCAxy", -1, "[0] in |DCAz| > [0]+[1]/pT"};
  Configurable<float> picMinDCAz{"picMinDCAz", -1, "[0] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> picMinPt{"picMinPt", -1, "Minimum pT for Xic pions"};

  Configurable<float> piccMinDCAxy{"piccMinDCAxy", -1, "[0] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> piccMinDCAz{"piccMinDCAz", -1, "[0] in |DCAz| > [0]+[1]/pT"};
  Configurable<float> piccMinPt{"piccMinPt", -1, "Minimum pT for Xicc pions"};

  Configurable<float> xicMaxDauDCA{"xicMaxDauDCA", 1e+4, "DCA between Xic daughters (cm)"};
  Configurable<float> xicMinDCAxy{"xicMinDCAxy", -1, "[0] in |DCAz| > [0]+[1]/pT"};
  Configurable<float> xicMinDCAz{"xicMinDCAz", -1, "[0] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> xiccMaxDCAxy{"xiccMaxDCAxy", 1e+4, "Maximum DCAxy"};
  Configurable<float> xiccMaxDCAz{"xiccMaxDCAz", 1e+4, "Maximum DCAz"};
  Configurable<float> xicMinRadius{"xicMinRadius", -1, "Minimum R2D for Xic decay (cm)"};
  Configurable<float> xicMinDecayDistanceFromPV{"xicMinDecayDistanceFromPV", -1, "Minimum distance for Xic decay from PV (cm)"};
  Configurable<float> xicMinProperLength{"xicMinProperLength", -1, "Minimum proper length for Xic decay (cm)"};
  Configurable<float> xicMaxProperLength{"xicMaxProperLength", 1e+4, "Minimum proper length for Xic decay (cm)"};

  Configurable<float> xiccMaxDauDCA{"xiccMaxDauDCA", 1e+4, "DCA between Xicc daughters (cm)"};
  Configurable<float> xiccMinRadius{"xiccMinRadius", -1, "Minimum R2D for Xicc decay (cm)"};
  Configurable<float> xiccMinProperLength{"xiccMinProperLength", -1, "Minimum proper length for Xicc decay (cm)"};
  Configurable<float> xiccMaxProperLength{"xiccMaxProperLength", 1e+4, "Minimum proper length for Xicc decay (cm)"};
  Configurable<int> otfConfig{"otfConfig", 0, "OTF configuration flag"};

  Filter configFilter = (aod::otfmulticharm::lutConfigId == otfConfig);

  void init(InitContext&)
  {
    histos.add("SelectionQA/hDCAXicDaughters", "hDCAXicDaughters; DCA between Xic daughters (#mum)", kTH1D, {axisDcaDaughters});
    histos.add("SelectionQA/hDCAXiccDaughters", "hDCAXiccDaughters; DCA between Xicc daughters (#mum)", kTH1D, {axisDcaDaughters});
    histos.add("SelectionQA/hDCAxyXi", "hDCAxyXi; Xi DCAxy to PV (#mum)", kTH1D, {axisDCA});
    histos.add("SelectionQA/hDCAzXi", "hDCAzXi; Xi DCAz to PV (#mum)", kTH1D, {axisDCA});
    histos.add("SelectionQA/hDCAxyXic", "hDCAxyXic; Xic DCAxy to PV (#mum)", kTH1D, {axisDCA});
    histos.add("SelectionQA/hDCAzXic", "hDCAzXic; Xic DCAz to PV (#mum)", kTH1D, {axisDCA});
    histos.add("SelectionQA/hDCAxyXicc", "hDCAxyXicc; Xicc DCAxy to PV (#mum)", kTH1D, {axisDCA});
    histos.add("SelectionQA/hDCAzXicc", "hDCAzXicc; Xicc DCAz to PV (#mum)", kTH1D, {axisDCA});
    histos.add("SelectionQA/hDecayRadiusXic", "hDecayRadiusXic; Distance (#mum)", kTH1D, {axisRadius});
    histos.add("SelectionQA/hDecayRadiusXicc", "hDecayRadiusXicc; Distance (#mum)", kTH1D, {axisRadius});
    histos.add("SelectionQA/hDecayDistanceFromPVXic", "hDecayDistanceFromPVXic; Distance (#mum)", kTH1D, {axisDecayLength});
    histos.add("SelectionQA/hProperLengthXic", "hProperLengthXic; Distance (#mum)", kTH1D, {axisDecayLength});
    histos.add("SelectionQA/hProperLengthXicc", "hProperLengthXicc; Distance (#mum)", kTH1D, {axisDecayLength});
    histos.add("SelectionQA/hPi1cDCAxy", "hPi1cDCAxy; Pi1c DCAxy (#mum)", kTH1D, {axisDCA});
    histos.add("SelectionQA/hPi1cDCAz", "hPi1cDCAz; Pi1c DCAz (#mum)", kTH1D, {axisDCA});
    histos.add("SelectionQA/hPi2cDCAxy", "hPi2cDCAxy; Pi2c DCAxy (#mum)", kTH1D, {axisDCA});
    histos.add("SelectionQA/hPi2cDCAz", "hPi2cDCAz; Pi2c DCAz (#mum)", kTH1D, {axisDCA});
    histos.add("SelectionQA/hPiccDCAxy", "hPiccDCAxy; Picc DCAxy (#mum)", kTH1D, {axisDCA});
    histos.add("SelectionQA/hPiccDCAz", "hPiccDCAz; Picc DCAz (#mum)", kTH1D, {axisDCA});
    histos.add("SelectionQA/hPi1cPt", "hPi1cPt; Pi1c pT (Gev/#it(c))", kTH1D, {axisPt});
    histos.add("SelectionQA/hPi2cPt", "hPi2cPt; Pi2c pT (Gev/#it(c))", kTH1D, {axisPt});
    histos.add("SelectionQA/hPiccPt", "hPiccPt; Picc pT (Gev/#it(c))", kTH1D, {axisPt});

    auto hMCharmBuilding = histos.add<TH1>("hMCharmBuilding", "hMCharmBuilding", kTH1D, {{22, -0.5, 21.5}});
    hMCharmBuilding->GetXaxis()->SetBinLabel(1, "nTotalCandidates");
    hMCharmBuilding->GetXaxis()->SetBinLabel(2, "xicMaxDauDCA");
    hMCharmBuilding->GetXaxis()->SetBinLabel(3, "xiccMaxDauDCA");
    hMCharmBuilding->GetXaxis()->SetBinLabel(4, "xiMinDCAxy");
    hMCharmBuilding->GetXaxis()->SetBinLabel(5, "xiMinDCAz");
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

    if (doprocessXiccExtra) {
      histos.add("XiccProngs/h3dPos", "h3dPos; Xicc pT (GeV/#it(c)); Pos pT (GeV/#it(c)); Pos #eta", kTH3D, {axisPt, axisPt, axisEta});
      histos.add("XiccProngs/h3dNeg", "h3dNeg; Xicc pT (GeV/#it(c)); Neg pT (GeV/#it(c)); Neg #eta", kTH3D, {axisPt, axisPt, axisEta});
      histos.add("XiccProngs/h3dBach", "h3dBach; Xicc pT (GeV/#it(c)); Bach pT (GeV/#it(c)); Bach #eta", kTH3D, {axisPt, axisPt, axisEta});
      histos.add("XiccProngs/h3dPi1c", "h3dPi1c; Xicc pT (GeV/#it(c)); Pi1c pT (GeV/#it(c)); Pi1c #eta", kTH3D, {axisPt, axisPt, axisEta});
      histos.add("XiccProngs/h3dPi2c", "h3dPi2c; Xicc pT (GeV/#it(c)); Pi2c pT (GeV/#it(c)); Pi2c #eta", kTH3D, {axisPt, axisPt, axisEta});
      histos.add("XiccProngs/h3dPicc", "h3dPicc; Xicc pT (GeV/#it(c)); Picc pT (GeV/#it(c)); Picc #eta", kTH3D, {axisPt, axisPt, axisEta});
    }

    histos.add("hXiccMass", "hXiccMass", kTH1D, {axisXiccMass});
    histos.add("hXicMass", "hXicMass", kTH1D, {axisXicMass});
    histos.add("hXiccPt", "hXiccPt", kTH1D, {axisPt});
    histos.add("hXicPt", "hXicPt", kTH1D, {axisPt});
    histos.add("h3dXicc", "h3dXicc; Xicc pT (GeV/#it(c)); Xicc #eta; Xicc mass (GeV/#it(c)^{2})", kTH3D, {axisPt, axisEta, axisXiccMass});
    histos.add("hConfigId", "hConfigId", kTH1D, {{11, -0.5, 10.5}});

    if (bdt.enableML) {
      ccdb->setURL(bdt.ccdbUrl.value);
      if (bdt.loadModelsFromCCDB) {
        ccdbApi.init(bdt.ccdbUrl);
        LOG(info) << "Fetching model for timestamp: " << bdt.timestampCCDB.value;
        bool retrieveSuccessMCharm = ccdbApi.retrieveBlob(bdt.pathCCDB.value, ".", metadata, bdt.timestampCCDB.value, false, bdt.localPath.value);

        if (retrieveSuccessMCharm) {
          bdtMCharm.initModel(bdt.localPath.value, bdt.enableOptimizations.value);
        } else {
          LOG(fatal) << "Error encountered while fetching/loading the MCharm model from CCDB! Maybe the model doesn't exist yet for this runnumber/timestamp?";
        }
      } else {
        bdtMCharm.initModel(bdt.localPath.value, bdt.enableOptimizations.value);
      }

      histos.add("hBDTScore", "hBDTScore", kTH1D, {axisBDTScore});
      histos.add("hBDTScoreVsXiccMass", "hBDTScoreVsXiccMass", kTH2D, {axisXiccMass, axisBDTScore});
      histos.add("hBDTScoreVsXiccPt", "hBDTScoreVsXiccPt", kTH2D, {axisPt, axisBDTScore});
      histos.add("h3dBDTScore", "h3dBDTScore", kTH3D, {axisPt, axisXiccMass, axisBDTScore});
      histos.add("hDCAXicDaughters", "hDCAXicDaughters", kTH2D, {{axisBDTScore, axisDcaDaughters}});
      histos.add("hDCAXiccDaughters", "hDCAXiccDaughters", kTH2D, {{axisBDTScore, axisDcaDaughters}});
      histos.add("hDCAxyXi", "hDCAxyXi", kTH2D, {{axisBDTScore, axisDCA}});
      histos.add("hDCAzXi", "hDCAzXi", kTH2D, {{axisBDTScore, axisDCA}});
      histos.add("hDCAxyXic", "hDCAxyXic", kTH2D, {{axisBDTScore, axisDCA}});
      histos.add("hDCAzXic", "hDCAzXic", kTH2D, {{axisBDTScore, axisDCA}});
      histos.add("hDCAxyXicc", "hDCAxyXicc", kTH2D, {{axisBDTScore, axisDCA}});
      histos.add("hDCAzXicc", "hDCAzXicc", kTH2D, {{axisBDTScore, axisDCA}});
      histos.add("hDecayRadiusXic", "hDecayRadiusXic", kTH2D, {{axisBDTScore, axisRadius}});
      histos.add("hDecayRadiusXicc", "hDecayRadiusXicc", kTH2D, {{axisBDTScore, axisRadius}});
      histos.add("hDecayDistanceFromPVXic", "hDecayDistanceFromPVXic", kTH2D, {{axisBDTScore, axisDecayLength}});
      histos.add("hProperLengthXic", "hProperLengthXic", kTH2D, {{axisBDTScore, axisDecayLength}});
      histos.add("hProperLengthXicc", "hProperLengthXicc", kTH2D, {{axisBDTScore, axisDecayLength}});
      histos.add("hPi1cDCAxy", "hPi1cDCAxy", kTH2D, {{axisBDTScore, axisDCA}});
      histos.add("hPi1cDCAz", "hPi1cDCAz", kTH2D, {{axisBDTScore, axisDCA}});
      histos.add("hPi2cDCAxy", "hPi2cDCAxy", kTH2D, {{axisBDTScore, axisDCA}});
      histos.add("hPi2cDCAz", "hPi2cDCAz", kTH2D, {{axisBDTScore, axisDCA}});
      histos.add("hPiccDCAxy", "hPiccDCAxy", kTH2D, {{axisBDTScore, axisDCA}});
      histos.add("hPiccDCAz", "hPiccDCAz", kTH2D, {{axisBDTScore, axisDCA}});
      histos.add("hPi1cPt", "hPi1cPt", kTH2D, {{axisBDTScore, axisPt}});
      histos.add("hPi2cPt", "hPi2cPt", kTH2D, {{axisBDTScore, axisPt}});
      histos.add("hPiccPt", "hPiccPt", kTH2D, {{axisBDTScore, axisPt}});
      histos.add("hXiccMass", "hXiccMass", kTH2D, {{axisBDTScore, axisXiccMass}});
      histos.add("hXicMass", "hXicMass", kTH2D, {{axisBDTScore, axisXicMass}});
      histos.add("hXiccPt", "hXiccPt", kTH2D, {{axisBDTScore, axisPt}});
      histos.add("hXicPt", "hXicPt", kTH2D, {{axisBDTScore, axisPt}});
    }
  }

  template <typename TMCharmCands>
  void genericProcessXicc(TMCharmCands const& xiccCands)
  {
    for (const auto& xiccCand : xiccCands) {
      int icfg = xiccCand.lutConfigId();
      histos.fill(HIST("hConfigId"), icfg);
      if (bdt.enableML) {
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

        float* probabilityMCharm = bdtMCharm.evalModel(inputFeatures);
        float bdtScore = probabilityMCharm[1];

        histos.fill(HIST("hBDTScore"), bdtScore);
        histos.fill(HIST("hBDTScoreVsXiccMass"), xiccCand.xiccMass(), bdtScore);
        histos.fill(HIST("hBDTScoreVsXiccPt"), xiccCand.xiccPt(), bdtScore);
        histos.fill(HIST("h3dBDTScore"), xiccCand.xiccPt(), xiccCand.xiccMass(), bdtScore);
        histos.fill(HIST("hDCAXicDaughters"), xiccCand.xiccPt(), bdtScore);
        histos.fill(HIST("hDCAXiccDaughters"), xiccCand.xiccPt(), bdtScore);
        histos.fill(HIST("hDCAxyXi"), xiccCand.xiccPt(), bdtScore);
        histos.fill(HIST("hDCAzXi"), xiccCand.xiccPt(), bdtScore);
        histos.fill(HIST("hDCAxyXic"), xiccCand.xiccPt(), bdtScore);
        histos.fill(HIST("hDCAzXic"), xiccCand.xiccPt(), bdtScore);
        histos.fill(HIST("hDCAxyXicc"), xiccCand.xiccPt(), bdtScore);
        histos.fill(HIST("hDCAzXicc"), xiccCand.xiccPt(), bdtScore);
        histos.fill(HIST("hDecayRadiusXic"), xiccCand.xiccPt(), bdtScore);
        histos.fill(HIST("hDecayRadiusXicc"), xiccCand.xiccPt(), bdtScore);
        histos.fill(HIST("hDecayDistanceFromPVXic"), xiccCand.xiccPt(), bdtScore);
        histos.fill(HIST("hProperLengthXic"), xiccCand.xiccPt(), bdtScore);
        histos.fill(HIST("hProperLengthXicc"), xiccCand.xiccPt(), bdtScore);
        histos.fill(HIST("hPi1cDCAxy"), xiccCand.xiccPt(), bdtScore);
        histos.fill(HIST("hPi1cDCAz"), xiccCand.xiccPt(), bdtScore);
        histos.fill(HIST("hPi2cDCAxy"), xiccCand.xiccPt(), bdtScore);
        histos.fill(HIST("hPi2cDCAz"), xiccCand.xiccPt(), bdtScore);
        histos.fill(HIST("hPiccDCAxy"), xiccCand.xiccPt(), bdtScore);
        histos.fill(HIST("hPiccDCAz"), xiccCand.xiccPt(), bdtScore);
        histos.fill(HIST("hPi1cPt"), xiccCand.xiccPt(), bdtScore);
        histos.fill(HIST("hPi2cPt"), xiccCand.xiccPt(), bdtScore);
        histos.fill(HIST("hPiccPt"), xiccCand.xiccPt(), bdtScore);
        histos.fill(HIST("hXiccMass"), xiccCand.xiccPt(), bdtScore);
        histos.fill(HIST("hXicMass"), xiccCand.xiccPt(), bdtScore);
        histos.fill(HIST("hXicPt"), xiccCand.xiccPt(), bdtScore);
        histos.fill(HIST("hXiccPt"), xiccCand.xiccPt(), bdtScore);
        histos.fill(HIST("h3dXicc"), xiccCand.xiccPt(), bdtScore);
      }

      histos.fill(HIST("hMCharmBuilding"), 0);
      if (xiccCand.xicDauDCA() > xicMaxDauDCA) {
        continue;
      } else {
        histos.fill(HIST("hMCharmBuilding"), 1);
      }

      if (xiccCand.xiccDauDCA() > xiccMaxDauDCA) {
        continue;
      } else {
        histos.fill(HIST("hMCharmBuilding"), 2);
      }

      if (std::fabs(xiccCand.xiDCAxy()) < xiMinDCAxy) {
        continue;
      } else {
        histos.fill(HIST("hMCharmBuilding"), 3);
      }

      if (std::fabs(xiccCand.xiDCAz()) < xiMinDCAz) {
        continue;
      } else {
        histos.fill(HIST("hMCharmBuilding"), 4);
      }

      if (std::fabs(xiccCand.pi1cDCAxy()) < picMinDCAxy) {
        continue;
      } else {
        histos.fill(HIST("hMCharmBuilding"), 5);
      }

      if (std::fabs(xiccCand.pi1cDCAz()) < picMinDCAz) {
        continue;
      } else {
        histos.fill(HIST("hMCharmBuilding"), 6);
      }

      if (std::fabs(xiccCand.pi2cDCAxy()) < picMinDCAxy) {
        continue;
      } else {
        histos.fill(HIST("hMCharmBuilding"), 7);
      }

      if (std::fabs(xiccCand.pi2cDCAz()) < picMinDCAz) {
        continue;
      } else {
        histos.fill(HIST("hMCharmBuilding"), 8);
      }

      if (std::fabs(xiccCand.piccDCAxy()) < piccMinDCAxy) {
        continue;
      } else {
        histos.fill(HIST("hMCharmBuilding"), 9);
      }

      if (std::fabs(xiccCand.piccDCAz()) < piccMinDCAz) {
        continue;
      } else {
        histos.fill(HIST("hMCharmBuilding"), 10);
      }

      if (std::fabs(xiccCand.xicDCAxy()) < xicMinDCAxy) {
        continue;
      } else {
        histos.fill(HIST("hMCharmBuilding"), 11);
      }

      if (std::fabs(xiccCand.xicDCAz()) < xicMinDCAz) {
        continue;
      } else {
        histos.fill(HIST("hMCharmBuilding"), 12);
      }

      if (std::fabs(xiccCand.xiccDCAxy()) > xiccMaxDCAxy) {
        continue;
      } else {
        histos.fill(HIST("hMCharmBuilding"), 13);
      }

      if (std::fabs(xiccCand.xiccDCAz()) > xiccMaxDCAz) {
        continue;
      } else {
        histos.fill(HIST("hMCharmBuilding"), 14);
      }

      if (xiccCand.xicDecayRadius2D() < xicMinRadius) {
        continue;
      } else {
        histos.fill(HIST("hMCharmBuilding"), 15);
      }

      if (xiccCand.xiccDecayRadius2D() < xiccMinRadius) {
        continue;
      } else {
        histos.fill(HIST("hMCharmBuilding"), 16);
      }

      if (xiccCand.xicProperLength() < xicMinProperLength) {
        continue;
      } else {
        histos.fill(HIST("hMCharmBuilding"), 17);
      }

      if (xiccCand.xicProperLength() > xicMaxProperLength) {
        continue;
      } else {
        histos.fill(HIST("hMCharmBuilding"), 18);
      }

      if (xiccCand.xiccProperLength() < xiccMinProperLength) {
        continue;
      } else {
        histos.fill(HIST("hMCharmBuilding"), 19);
      }

      if (xiccCand.xiccProperLength() > xiccMaxProperLength) {
        continue;
      } else {
        histos.fill(HIST("hMCharmBuilding"), 20);
      }

      if (xiccCand.xicDistanceFromPV() < xicMinDecayDistanceFromPV) {
        continue;
      } else {
        histos.fill(HIST("hMCharmBuilding"), 21);
      }

      histos.fill(HIST("SelectionQA/hDCAXicDaughters"), xiccCand.xicDauDCA() * 1e+4);
      histos.fill(HIST("SelectionQA/hDCAXiccDaughters"), xiccCand.xiccDauDCA() * 1e+4);
      histos.fill(HIST("SelectionQA/hDCAxyXi"), std::fabs(xiccCand.xiDCAxy() * 1e+4));
      histos.fill(HIST("SelectionQA/hDCAzXi"), std::fabs(xiccCand.xiDCAz() * 1e+4));
      histos.fill(HIST("SelectionQA/hDCAxyXic"), std::fabs(xiccCand.xicDCAxy() * 1e+4));
      histos.fill(HIST("SelectionQA/hDCAzXic"), std::fabs(xiccCand.xicDCAz() * 1e+4));
      histos.fill(HIST("SelectionQA/hDCAxyXicc"), std::fabs(xiccCand.xiccDCAxy() * 1e+4));
      histos.fill(HIST("SelectionQA/hDCAzXicc"), std::fabs(xiccCand.xiccDCAz() * 1e+4));
      histos.fill(HIST("SelectionQA/hDecayRadiusXic"), xiccCand.xicDecayRadius2D() * 1e+4);
      histos.fill(HIST("SelectionQA/hDecayRadiusXicc"), xiccCand.xiccDecayRadius2D() * 1e+4);
      histos.fill(HIST("SelectionQA/hDecayDistanceFromPVXic"), xiccCand.xicDistanceFromPV() * 1e+4);
      histos.fill(HIST("SelectionQA/hProperLengthXic"), xiccCand.xicProperLength() * 1e+4);
      histos.fill(HIST("SelectionQA/hProperLengthXicc"), xiccCand.xiccProperLength() * 1e+4);
      histos.fill(HIST("SelectionQA/hPi1cDCAxy"), xiccCand.pi1cDCAxy() * 1e+4);
      histos.fill(HIST("SelectionQA/hPi1cDCAz"), xiccCand.pi1cDCAz() * 1e+4);
      histos.fill(HIST("SelectionQA/hPi2cDCAxy"), xiccCand.pi2cDCAxy() * 1e+4);
      histos.fill(HIST("SelectionQA/hPi2cDCAz"), xiccCand.pi2cDCAz() * 1e+4);
      histos.fill(HIST("SelectionQA/hPiccDCAxy"), xiccCand.piccDCAxy() * 1e+4);
      histos.fill(HIST("SelectionQA/hPiccDCAz"), xiccCand.piccDCAz() * 1e+4);
      histos.fill(HIST("SelectionQA/hPi1cPt"), xiccCand.pi1cPt());
      histos.fill(HIST("SelectionQA/hPi2cPt"), xiccCand.pi2cPt());
      histos.fill(HIST("SelectionQA/hPiccPt"), xiccCand.piccPt());

      if constexpr (requires { xiccCand.negPt(); }) { // if extra table
        histos.fill(HIST("XiccProngs/h3dNeg"), xiccCand.xiccPt(), xiccCand.negPt(), xiccCand.negEta());
        histos.fill(HIST("XiccProngs/h3dPos"), xiccCand.xiccPt(), xiccCand.posPt(), xiccCand.posEta());
        histos.fill(HIST("XiccProngs/h3dBach"), xiccCand.xiccPt(), xiccCand.bachPt(), xiccCand.bachEta());
        histos.fill(HIST("XiccProngs/h3dPi1c"), xiccCand.xiccPt(), xiccCand.pi1cPt(), xiccCand.pi1cEta());
        histos.fill(HIST("XiccProngs/h3dPi2c"), xiccCand.xiccPt(), xiccCand.pi2cPt(), xiccCand.pi2cEta());
        histos.fill(HIST("XiccProngs/h3dPicc"), xiccCand.xiccPt(), xiccCand.piccPt(), xiccCand.piccEta());
      }

      histos.fill(HIST("hXiccMass"), xiccCand.xiccMass());
      histos.fill(HIST("hXicMass"), xiccCand.xicMass());
      histos.fill(HIST("hXiccPt"), xiccCand.xiccPt());
      histos.fill(HIST("hXicPt"), xiccCand.xicPt());
      histos.fill(HIST("h3dXicc"), xiccCand.xiccPt(), xiccCand.xiccEta(), xiccCand.xiccMass());
    }
  }

  void processXicc(soa::Filtered<aod::MCharmCores> const& multiCharmTracks)
  {
    genericProcessXicc(multiCharmTracks);
  }

  void processXiccExtra(soa::Filtered<MultiCharmTracksFull> const& multiCharmTracks)
  {
    genericProcessXicc(multiCharmTracks);
  }

  PROCESS_SWITCH(Alice3Multicharm, processXicc, "find Xicc baryons", true);
  PROCESS_SWITCH(Alice3Multicharm, processXiccExtra, "find Xicc baryons with all QA", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<Alice3Multicharm>(cfgc)};
}
