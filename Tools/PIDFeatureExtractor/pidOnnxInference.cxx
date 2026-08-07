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
///        Deliberately NOT an AOD-table-subscribing DPL task: it reads the
///        input file directly via plain TFile/TTree in init(), same as
///        pidFeatureExtractor.cxx now writes its output - no
///        DECLARE_SOA_TABLE, no Produces<>, avoiding the framework issue
///        that broke the earlier table-based version of this pair of
///        tasks tonight. Uses o2::analysis::MlResponse
///        (Tools/ML/MlResponse.h) - O2Physics's generic ONNX/CCDB
///        inference wrapper - for model loading and execution.
///
/// \author Robert Forynski

#include "Tools/ML/MlResponse.h"

#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/ControlService.h>
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
} // namespace

/// PidOnnxInference: applies the FSE ONNX model to the features written by
/// pidFeatureExtractor.cxx and writes out per-track class probabilities.
///
/// This is a one-shot batch task, not a per-timeframe DPL processor: all
/// real work happens in init() (load model, read the whole input tree,
/// write predictions), since the input is a finished file from a previous
/// run rather than live AO2D data. A trivial empty process() is kept so
/// the task still registers normally with adaptAnalysisTask, matching
/// every other task in this project - NOT VERIFIED that
/// ControlService::endOfStream()/readyToQuit() is the correct way to make
/// the workflow terminate cleanly after init() finishes; if the workflow
/// hangs instead of exiting, this is the first thing to revisit.
///
/// - cutDirMl defaults to cuts_ml::CutNot (value 2) for every class, so
///   MlResponse never rejects a track - this task always reports all four
///   probabilities rather than making a pass/fail decision.
/// - A single pT bin is used by default (model isn't pT-binned); lower
///   edge is -1, not 0, so no track's pt() can land exactly on the
///   boundary (MlResponse::findBin() treats that as out-of-range).
/// - buildModelInput()'s feature order is a REASONABLE DEFAULT, not
///   verified against the actual training code - see the comment there.
struct PidOnnxInference {
  // --- input: the ROOT file/tree pidFeatureExtractor.cxx wrote -------------
  Configurable<std::string> inputRootFile{"inputRootFile", "pid_features_data.root", "ROOT file produced by pidFeatureExtractor.cxx"};
  Configurable<std::string> inputTreeName{"inputTreeName", "pid_features", "Name of the TTree inside inputRootFile"};

  // --- output ------------------------------------------------------------------
  Configurable<std::string> outputPath{"outputPath", "pid_predictions", "Output file base name (no extension)"};
  Configurable<bool> exportCsv{"exportCsv", false, "Also write predictions to CSV alongside the ROOT output"};

  // --- model location ----------------------------------------------------------
  Configurable<bool> loadModelFromCcdb{"loadModelFromCcdb", true, "Load the ONNX model from CCDB (else from onnxFileNames as a local path)"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "CCDB URL"};
  Configurable<std::vector<std::string>> modelPathsCcdb{"modelPathsCcdb", std::vector<std::string>{"Users/YOURNAME/PidFeatureExtractor/model"}, "CCDB path to the model"};
  Configurable<int64_t> timestampCcdb{"timestampCcdb", -1, "CCDB query timestamp for the model, -1 = latest"};
  Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"pid_feature_model.onnx"}, "Local ONNX file path(s), used when loadModelFromCcdb is false"};

  // --- MlResponse plumbing: a single pT bin, no selection cut applied --------
  Configurable<std::vector<double>> binsPtMl{"binsPtMl", std::vector<double>{-1., 9999.}, "pT bin edges for MlResponse (single bin = model isn't pT-binned)"};

  static constexpr double kDefaultCutsMl[1][kNumClasses] = {{0., 0., 0., 0.}};
  Configurable<LabeledArray<double>> cutsMl{"cutsMl", {kDefaultCutsMl[0], 1, kNumClasses, {"pT bin 0"}, {"prob pi", "prob ka", "prob pr", "prob el"}}, "Unused thresholds (CutNot everywhere) - required by MlResponse's interface"};
  Configurable<int8_t> nClassesMl{"nClassesMl", static_cast<int8_t>(kNumClasses), "Number of model output classes"};

  o2::ccdb::CcdbApi ccdbApi;
  o2::analysis::MlResponse<float> mlResponse;

  void init(InitContext& ic)
  {
    mlResponse.configure(binsPtMl, cutsMl, std::vector<int>{cuts_ml::CutNot, cuts_ml::CutNot, cuts_ml::CutNot, cuts_ml::CutNot}, nClassesMl);
    if (loadModelFromCcdb.value) {
      ccdbApi.init(ccdbUrl.value);
      mlResponse.setModelPathsCCDB(onnxFileNames, ccdbApi, modelPathsCcdb.value, timestampCcdb.value);
    } else {
      mlResponse.setModelPathsLocal(onnxFileNames);
    }
    mlResponse.init();

    runInferenceOverFile();

    // One-shot batch job: nothing to subscribe to, so tell DPL we're done
    // rather than waiting for input that will never arrive.
    ic.services().get<ControlService>().endOfStream();
    ic.services().get<ControlService>().readyToQuit(QuitRequest::Me);
  }

  /// Trivial, intentionally empty - exists only so this task registers
  /// like every other task in the project. All real work is in init().
  void process(ProcessingContext&) {}

  void runInferenceOverFile()
  {
    std::unique_ptr<TFile> inFile(TFile::Open(inputRootFile.value.c_str(), "READ"));
    if (!inFile || inFile->IsZombie()) {
      LOG(fatal) << "Could not open input file " << inputRootFile.value;
      return;
    }
    auto* tree = dynamic_cast<TTree*>(inFile->Get(inputTreeName.value.c_str()));
    if (!tree) {
      LOG(fatal) << "Tree " << inputTreeName.value << " not found in " << inputRootFile.value;
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

    std::unique_ptr<TFile> outFile(TFile::Open((outputPath.value + ".root").c_str(), "RECREATE"));
    TTree outTree("pid_predictions", "PID ML predictions");
    float mlProbPi = 0, mlProbKa = 0, mlProbPr = 0, mlProbEl = 0;
    int mlPredictedClass = 0;
    outTree.Branch("mlProbPi", &mlProbPi);
    outTree.Branch("mlProbKa", &mlProbKa);
    outTree.Branch("mlProbPr", &mlProbPr);
    outTree.Branch("mlProbEl", &mlProbEl);
    outTree.Branch("mlPredictedClass", &mlPredictedClass);

    std::ofstream csv;
    if (exportCsv.value) {
      csv.open(outputPath.value + ".csv");
      csv << "mlProbPi,mlProbKa,mlProbPr,mlProbEl,mlPredictedClass\n";
    }

    // --------------------------------------------------------------------------
    // Feature order fed to the model - THIS MUST MATCH YOUR TRAINING SCRIPT'S
    // COLUMN ORDER EXACTLY. Reasonable default (every reconstructed feature
    // except vz/centFT0C/sign/trackType and the Bayesian columns, which are a
    // comparison baseline, not a model input), followed by a 7-length group
    // mask (TPC/TOF/TRD/ITS/EMCal/HMPID/centrality; ITS and centrality are
    // assumed always-present). NOT verified against your actual training code.
    // --------------------------------------------------------------------------
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

      if (exportCsv.value) {
        csv << mlProbPi << ',' << mlProbKa << ',' << mlProbPr << ',' << mlProbEl << ',' << mlPredictedClass << '\n';
      }
    }

    outFile->cd();
    outTree.Write();
    outFile->Close();
    if (exportCsv.value) {
      csv.close();
    }

    LOG(info) << "PidOnnxInference: wrote " << nEntries << " predictions to " << outputPath.value << ".root";
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PidOnnxInference>(cfgc)};
}
