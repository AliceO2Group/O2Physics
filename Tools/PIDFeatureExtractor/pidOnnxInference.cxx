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

/// \file pidOnnxInference.cxx
/// \brief Run the FSE PID ONNX model (loaded from CCDB, or a local file for
///        testing) over the ROOT TTree produced by pidFeatureExtractor.cxx,
///        and write per-track class probabilities to a new ROOT file
///        (and/or CSV).
///
///        Built as a plain DataProcessorSpec/AlgorithmSpec, NOT
///        adaptAnalysisTask<> - this is a standalone batch job with no AOD
///        table subscription, and adaptAnalysisTask's process() reflection
///        specifically requires an AOD table/iterator argument (confirmed
///        by a real compiler error: a raw ProcessingContext& process()
///        signature is rejected outright). Using the lower-level DPL
///        primitives directly avoids that machinery entirely, rather than
///        trying to satisfy a reflection mechanism built for a different
///        kind of task. Configurable<> member auto-binding is an
///        AnalysisTask-specific convenience, so options here are declared
///        explicitly and read via ic.options().get<T>(...) instead.
///
///        Uses o2::analysis::MlResponse (Tools/ML/MlResponse.h) -
///        O2Physics's generic ONNX/CCDB inference wrapper - for model
///        loading and execution.
///
/// \author Robert Forynski

#include "Tools/ML/MlResponse.h"

#include <Framework/ConfigParamSpec.h>
#include <Framework/ControlService.h>
#include <Framework/DataProcessorSpec.h>
#include <Framework/runDataProcessing.h>

#include <TFile.h>
#include <TTree.h>

#include <cstdint>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;

namespace
{
constexpr int kNumClasses = 4; // pi, ka, pr, el - fixed order throughout, matches the paper's model

/// itsClusterSizes packs 7 ITS layers into 4 bits each; a derived cluster
/// count is a far more sensible model input than the raw packed value.
int getItsNClusters(uint32_t v)
{
  int n = 0;
  for (int layer = 0; layer < 7; layer++) {
    if ((v >> (layer * 4)) & 0xF) {
      n++;
    }
  }
  return n;
}

int argmax4(std::vector<float> const& v)
{
  int best = 0;
  for (int i = 1; i < kNumClasses; i++) {
    if (v[i] > v[best]) {
      best = i;
    }
  }
  return best;
}

/// All the real work: load the model, read the whole input tree, run
/// inference row by row, write predictions. Called once from init().
///
/// Feature order fed to the model - THIS MUST MATCH YOUR TRAINING SCRIPT'S
/// COLUMN ORDER EXACTLY. Reasonable default (every reconstructed feature
/// except vz/centFT0C/sign/trackType and the Bayesian columns, which are a
/// comparison baseline, not a model input), followed by a 7-length group
/// mask (TPC/TOF/TRD/ITS/EMCal/HMPID/centrality; ITS and centrality are
/// assumed always-present). NOT verified against your actual training code.
void runInference(std::string const& inputRootFile, std::string const& inputTreeName,
                   std::string const& outputPath, bool exportCsv,
                   o2::analysis::MlResponse<float>& mlResponse)
{
  std::unique_ptr<TFile> inFile(TFile::Open(inputRootFile.c_str(), "READ"));
  if (!inFile || inFile->IsZombie()) {
    LOG(fatal) << "Could not open input file " << inputRootFile;
    return;
  }
  auto* tree = dynamic_cast<TTree*>(inFile->Get(inputTreeName.c_str()));
  if (!tree) {
    LOG(fatal) << "Tree " << inputTreeName << " not found in " << inputRootFile;
    return;
  }

  // Bind the input branches actually used below - must match
  // pidFeatureExtractor.cxx's branch names exactly.
  float p = 0, pt = 0, px = 0, py = 0, pz = 0, eta = 0, phi = 0;
  float dcaXY = 0, dcaZ = 0;
  bool hasTPC = false;
  float tpcSignal = 0, tpcNSigmaPi = 0, tpcNSigmaKa = 0, tpcNSigmaPr = 0, tpcNSigmaEl = 0;
  int tpcNClsFound = 0;
  float tpcChi2NCl = 0;
  bool hasTOF = false;
  float tofMass = 0, beta = 0, tofNSigmaPi = 0, tofNSigmaKa = 0, tofNSigmaPr = 0, tofNSigmaEl = 0;
  bool hasTRD = false;
  float trdSignal = 0, trdChi2 = 0;
  int trdPattern = 0;
  int itsClusterSizes = 0;
  float itsChi2NCl = 0;
  bool hasEMCal = false;
  float trackEtaEmcal = 0, trackPhiEmcal = 0;
  bool hasHMPID = false;
  float hmpidSignal = 0, hmpidQMip = 0;
  int hmpidNPhotons = 0, hmpidClusSize = 0;
  float hmpidMom = 0;

  tree->SetBranchAddress("p", &p);
  tree->SetBranchAddress("pt", &pt);
  tree->SetBranchAddress("px", &px);
  tree->SetBranchAddress("py", &py);
  tree->SetBranchAddress("pz", &pz);
  tree->SetBranchAddress("eta", &eta);
  tree->SetBranchAddress("phi", &phi);
  tree->SetBranchAddress("dcaXY", &dcaXY);
  tree->SetBranchAddress("dcaZ", &dcaZ);
  tree->SetBranchAddress("hasTPC", &hasTPC);
  tree->SetBranchAddress("tpcSignal", &tpcSignal);
  tree->SetBranchAddress("tpcNSigmaPi", &tpcNSigmaPi);
  tree->SetBranchAddress("tpcNSigmaKa", &tpcNSigmaKa);
  tree->SetBranchAddress("tpcNSigmaPr", &tpcNSigmaPr);
  tree->SetBranchAddress("tpcNSigmaEl", &tpcNSigmaEl);
  tree->SetBranchAddress("tpcNClsFound", &tpcNClsFound);
  tree->SetBranchAddress("tpcChi2NCl", &tpcChi2NCl);
  tree->SetBranchAddress("hasTOF", &hasTOF);
  tree->SetBranchAddress("tofMass", &tofMass);
  tree->SetBranchAddress("beta", &beta);
  tree->SetBranchAddress("tofNSigmaPi", &tofNSigmaPi);
  tree->SetBranchAddress("tofNSigmaKa", &tofNSigmaKa);
  tree->SetBranchAddress("tofNSigmaPr", &tofNSigmaPr);
  tree->SetBranchAddress("tofNSigmaEl", &tofNSigmaEl);
  tree->SetBranchAddress("hasTRD", &hasTRD);
  tree->SetBranchAddress("trdSignal", &trdSignal);
  tree->SetBranchAddress("trdChi2", &trdChi2);
  tree->SetBranchAddress("trdPattern", &trdPattern);
  tree->SetBranchAddress("itsClusterSizes", &itsClusterSizes);
  tree->SetBranchAddress("itsChi2NCl", &itsChi2NCl);
  tree->SetBranchAddress("hasEMCal", &hasEMCal);
  tree->SetBranchAddress("trackEtaEmcal", &trackEtaEmcal);
  tree->SetBranchAddress("trackPhiEmcal", &trackPhiEmcal);
  tree->SetBranchAddress("hasHMPID", &hasHMPID);
  tree->SetBranchAddress("hmpidSignal", &hmpidSignal);
  tree->SetBranchAddress("hmpidQMip", &hmpidQMip);
  tree->SetBranchAddress("hmpidNPhotons", &hmpidNPhotons);
  tree->SetBranchAddress("hmpidClusSize", &hmpidClusSize);
  tree->SetBranchAddress("hmpidMom", &hmpidMom);

  std::unique_ptr<TFile> outFile(TFile::Open((outputPath + ".root").c_str(), "RECREATE"));
  TTree outTree("pid_predictions", "PID ML predictions");
  float mlProbPi = 0, mlProbKa = 0, mlProbPr = 0, mlProbEl = 0;
  int mlPredictedClass = 0;
  outTree.Branch("mlProbPi", &mlProbPi);
  outTree.Branch("mlProbKa", &mlProbKa);
  outTree.Branch("mlProbPr", &mlProbPr);
  outTree.Branch("mlProbEl", &mlProbEl);
  outTree.Branch("mlPredictedClass", &mlPredictedClass);

  std::ofstream csv;
  if (exportCsv) {
    csv.open(outputPath + ".csv");
    csv << "mlProbPi,mlProbKa,mlProbPr,mlProbEl,mlPredictedClass\n";
  }

  std::vector<float> x;
  std::vector<float> mlOutput;
  Long64_t nEntries = tree->GetEntries();
  for (Long64_t i = 0; i < nEntries; i++) {
    tree->GetEntry(i);

    x.clear();
    x.reserve(38 + 7);
    x.push_back(p);
    x.push_back(pt);
    x.push_back(px);
    x.push_back(py);
    x.push_back(pz);
    x.push_back(eta);
    x.push_back(phi);
    x.push_back(dcaXY);
    x.push_back(dcaZ);
    x.push_back(static_cast<float>(hasTPC));
    x.push_back(tpcSignal);
    x.push_back(tpcNSigmaPi);
    x.push_back(tpcNSigmaKa);
    x.push_back(tpcNSigmaPr);
    x.push_back(tpcNSigmaEl);
    x.push_back(static_cast<float>(tpcNClsFound));
    x.push_back(tpcChi2NCl);
    x.push_back(static_cast<float>(hasTOF));
    x.push_back(tofMass);
    x.push_back(beta);
    x.push_back(tofNSigmaPi);
    x.push_back(tofNSigmaKa);
    x.push_back(tofNSigmaPr);
    x.push_back(tofNSigmaEl);
    x.push_back(static_cast<float>(hasTRD));
    x.push_back(trdSignal);
    x.push_back(trdChi2);
    x.push_back(static_cast<float>(trdPattern));
    x.push_back(static_cast<float>(getItsNClusters(static_cast<uint32_t>(itsClusterSizes))));
    x.push_back(itsChi2NCl);
    x.push_back(static_cast<float>(hasEMCal));
    x.push_back(trackEtaEmcal);
    x.push_back(trackPhiEmcal);
    x.push_back(static_cast<float>(hasHMPID));
    x.push_back(hmpidSignal);
    x.push_back(hmpidQMip);
    x.push_back(static_cast<float>(hmpidNPhotons));
    x.push_back(static_cast<float>(hmpidClusSize));
    x.push_back(hmpidMom);
    // 7-length group mask
    x.push_back(static_cast<float>(hasTPC));
    x.push_back(static_cast<float>(hasTOF));
    x.push_back(static_cast<float>(hasTRD));
    x.push_back(1.f); // ITS
    x.push_back(static_cast<float>(hasEMCal));
    x.push_back(static_cast<float>(hasHMPID));
    x.push_back(1.f); // centrality

    mlResponse.isSelectedMl(x, pt, mlOutput); // return value (selection) unused; mlOutput carries the 4 raw scores
    mlProbPi = mlOutput[0];
    mlProbKa = mlOutput[1];
    mlProbPr = mlOutput[2];
    mlProbEl = mlOutput[3];
    mlPredictedClass = argmax4(mlOutput);
    outTree.Fill();

    if (exportCsv) {
      csv << mlProbPi << ',' << mlProbKa << ',' << mlProbPr << ',' << mlProbEl << ',' << mlPredictedClass << '\n';
    }
  }

  outFile->cd();
  outTree.Write();
  outFile->Close();
  if (exportCsv) {
    csv.close();
  }

  LOG(info) << "PidOnnxInference: wrote " << nEntries << " predictions to " << outputPath << ".root";
}
} // namespace

WorkflowSpec defineDataProcessing(ConfigContext const&)
{
  DataProcessorSpec spec{
    "pid-onnx-inference",
    Inputs{},
    Outputs{},
    AlgorithmSpec{[](InitContext& ic) {
      auto inputRootFile = ic.options().get<std::string>("inputRootFile");
      auto inputTreeName = ic.options().get<std::string>("inputTreeName");
      auto outputPath = ic.options().get<std::string>("outputPath");
      auto exportCsv = ic.options().get<bool>("exportCsv");
      auto loadModelFromCcdb = ic.options().get<bool>("loadModelFromCcdb");
      auto ccdbUrl = ic.options().get<std::string>("ccdbUrl");
      auto modelPathsCcdb = ic.options().get<std::vector<std::string>>("modelPathsCcdb");
      auto timestampCcdb = ic.options().get<int64_t>("timestampCcdb");
      auto onnxFileNames = ic.options().get<std::vector<std::string>>("onnxFileNames");
      auto binsPtMl = ic.options().get<std::vector<double>>("binsPtMl");
      auto nClassesMl = static_cast<int8_t>(ic.options().get<int>("nClassesMl"));

      // Unused thresholds (CutNot everywhere) - this task always reports
      // all four probabilities rather than applying a selection cut, so
      // cutsMl/cutDirMl don't need to be user-configurable; hardcoded here
      // rather than exposed as options.
      static constexpr double kDefaultCutsMl[1][kNumClasses] = {{0., 0., 0., 0.}};
      LabeledArray<double> cutsMl{kDefaultCutsMl[0], 1, kNumClasses, {"pT bin 0"}, {"prob pi", "prob ka", "prob pr", "prob el"}};
      std::vector<int> cutDirMl{cuts_ml::CutNot, cuts_ml::CutNot, cuts_ml::CutNot, cuts_ml::CutNot};

      auto mlResponse = std::make_shared<o2::analysis::MlResponse<float>>();
      mlResponse->configure(binsPtMl, cutsMl, cutDirMl, nClassesMl);
      if (loadModelFromCcdb) {
        auto ccdbApi = std::make_shared<o2::ccdb::CcdbApi>();
        ccdbApi->init(ccdbUrl);
        mlResponse->setModelPathsCCDB(onnxFileNames, *ccdbApi, modelPathsCcdb, timestampCcdb);
      } else {
        mlResponse->setModelPathsLocal(onnxFileNames);
      }
      mlResponse->init();

      runInference(inputRootFile, inputTreeName, outputPath, exportCsv, *mlResponse);

      return [](ProcessingContext& pc) {
        // One-shot batch job: all work already happened in init(). Signal
        // completion immediately rather than waiting for input that will
        // never arrive.
        pc.services().get<ControlService>().endOfStream();
        pc.services().get<ControlService>().readyToQuit(QuitRequest::Me);
      };
    }},
    Options{
      {"inputRootFile", VariantType::String, "pid_features_data.root", {"ROOT file produced by pidFeatureExtractor.cxx"}},
      {"inputTreeName", VariantType::String, "pid_features", {"Name of the TTree inside inputRootFile"}},
      {"outputPath", VariantType::String, "pid_predictions", {"Output file base name (no extension)"}},
      {"exportCsv", VariantType::Bool, false, {"Also write predictions to CSV alongside the ROOT output"}},
      {"loadModelFromCcdb", VariantType::Bool, true, {"Load the ONNX model from CCDB (else from onnxFileNames as a local path)"}},
      {"ccdbUrl", VariantType::String, "http://alice-ccdb.cern.ch", {"CCDB URL"}},
      {"modelPathsCcdb", VariantType::ArrayString, std::vector<std::string>{"Users/YOURNAME/PidFeatureExtractor/model"}, {"CCDB path to the model"}},
      {"timestampCcdb", VariantType::Int64, static_cast<int64_t>(-1), {"CCDB query timestamp for the model, -1 = latest"}},
      {"onnxFileNames", VariantType::ArrayString, std::vector<std::string>{"pid_feature_model.onnx"}, {"Local ONNX file path(s), used when loadModelFromCcdb is false"}},
      {"binsPtMl", VariantType::ArrayDouble, std::vector<double>{-1., 9999.}, {"pT bin edges for MlResponse (single bin = model isn't pT-binned)"}},
      {"nClassesMl", VariantType::Int, kNumClasses, {"Number of model output classes"}},
    }};

  return WorkflowSpec{spec};
}
