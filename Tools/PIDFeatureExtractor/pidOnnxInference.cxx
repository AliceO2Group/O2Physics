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
///        A normal Configurable<>-based adaptAnalysisTask, as required by
///        this repository's conventions (workflow topology must be
///        expressed via process function switches / Configurable<>, not
///        hand-built DataProcessorSpec Options{}). All real work still
///        happens once, in init() - it reads the extractor's output file
///        directly via plain TFile/TTree, not an AOD table. process() is
///        intentionally a no-op with a minimal, always-valid AOD argument
///        (aod::Collisions), present only so this registers as a normal
///        analysis task; it does nothing per-collision. Running this task
///        therefore still requires a valid AO2D input to satisfy the
///        pipeline, even though its content is unused.
///
///        Uses o2::analysis::MlResponse (Tools/ML/MlResponse.h) -
///        O2Physics's generic ONNX/CCDB inference wrapper - for model
///        loading and execution.
///
/// \author Robert Forynski

#include "Tools/ML/MlResponse.h"

#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/runDataProcessing.h>

#include <TFile.h>
#include <TTree.h>

#include <cstdint>
#include <fstream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;

namespace
{
constexpr int kNumClasses = 4; // pi, ka, pr, el - fixed order throughout, matches the paper's model
constexpr float kNaN = std::numeric_limits<float>::quiet_NaN();
constexpr int kNumItsLayers = 7;
constexpr int kBitsPerItsLayer = 4;
constexpr uint32_t kItsLayerMask = 0xF;

/// itsClusterSizes packs 7 ITS layers into 4 bits each; a derived cluster
/// count is a far more sensible model input than the raw packed value.
int getItsNClusters(uint32_t v)
{
  int n = 0;
  for (int layer = 0; layer < kNumItsLayers; layer++) {
    if ((v >> (layer * kBitsPerItsLayer)) & kItsLayerMask) {
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

/// Per-group enable/disable, independent of what the input tree's hasXXX
/// flags say. Default is everything enabled (true) - the normal case,
/// using each track's real detector coverage as-is. Turning a group off
/// forces its features to the same "absent" sentinel used when the
/// detector genuinely didn't fire, and clears its mask bit - useful for
/// testing how the model behaves with a detector deliberately excluded,
/// independent of the data itself.
struct GroupToggles {
  bool useTPC = true;
  bool useTOF = true;
  bool useTRD = true;
  bool useITS = true;
  bool useEMCal = true;
  bool useHMPID = true;
  bool useCentrality = true;
};

/// All the real work: load the model, read the whole input tree, run
/// inference row by row, write predictions. Called once from init().
///
/// Feature order fed to the model - THIS MUST MATCH YOUR TRAINING SCRIPT'S
/// COLUMN ORDER EXACTLY. Reasonable default (every reconstructed feature
/// except vz/centFT0C/sign/trackType and the Bayesian columns, which are a
/// comparison baseline, not a model input), followed by a 7-length group
/// mask (TPC/TOF/TRD/ITS/EMCal/HMPID/centrality). Each group can be
/// disabled via GroupToggles regardless of what the data says - see there
/// for details. ITS and centrality have no hasXXX flag in the input tree
/// (assumed always-present in the data itself), so their toggle is the
/// only way to exclude them.
void runInference(std::string const& inputRootFile, std::string const& inputTreeName,
                  std::string const& outputPath, bool exportCsv,
                  o2::analysis::MlResponse<float>& mlResponse, GroupToggles const& groups)
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

    // Effective presence = what the data says AND the group is enabled.
    bool effTPC = hasTPC && groups.useTPC;
    bool effTOF = hasTOF && groups.useTOF;
    bool effTRD = hasTRD && groups.useTRD;
    bool effEMCal = hasEMCal && groups.useEMCal;
    bool effHMPID = hasHMPID && groups.useHMPID;

    x.clear();
    x.reserve(39 + 7);
    x.push_back(p);
    x.push_back(pt);
    x.push_back(px);
    x.push_back(py);
    x.push_back(pz);
    x.push_back(eta);
    x.push_back(phi);
    x.push_back(dcaXY);
    x.push_back(dcaZ);
    x.push_back(static_cast<float>(effTPC));
    x.push_back(effTPC ? tpcSignal : kNaN);
    x.push_back(effTPC ? tpcNSigmaPi : kNaN);
    x.push_back(effTPC ? tpcNSigmaKa : kNaN);
    x.push_back(effTPC ? tpcNSigmaPr : kNaN);
    x.push_back(effTPC ? tpcNSigmaEl : kNaN);
    x.push_back(effTPC ? static_cast<float>(tpcNClsFound) : 0.f);
    x.push_back(effTPC ? tpcChi2NCl : kNaN);
    x.push_back(static_cast<float>(effTOF));
    x.push_back(effTOF ? tofMass : kNaN);
    x.push_back(effTOF ? beta : kNaN);
    x.push_back(effTOF ? tofNSigmaPi : kNaN);
    x.push_back(effTOF ? tofNSigmaKa : kNaN);
    x.push_back(effTOF ? tofNSigmaPr : kNaN);
    x.push_back(effTOF ? tofNSigmaEl : kNaN);
    x.push_back(static_cast<float>(effTRD));
    x.push_back(effTRD ? trdSignal : kNaN);
    x.push_back(effTRD ? trdChi2 : kNaN);
    x.push_back(effTRD ? static_cast<float>(trdPattern) : 0.f);
    x.push_back(groups.useITS ? static_cast<float>(getItsNClusters(static_cast<uint32_t>(itsClusterSizes))) : 0.f);
    x.push_back(groups.useITS ? itsChi2NCl : kNaN);
    x.push_back(static_cast<float>(effEMCal));
    x.push_back(effEMCal ? trackEtaEmcal : kNaN);
    x.push_back(effEMCal ? trackPhiEmcal : kNaN);
    x.push_back(static_cast<float>(effHMPID));
    x.push_back(effHMPID ? hmpidSignal : kNaN);
    x.push_back(effHMPID ? hmpidQMip : kNaN);
    x.push_back(effHMPID ? static_cast<float>(hmpidNPhotons) : 0.f);
    x.push_back(effHMPID ? static_cast<float>(hmpidClusSize) : 0.f);
    x.push_back(effHMPID ? hmpidMom : kNaN);
    // 7-length group mask
    x.push_back(static_cast<float>(effTPC));
    x.push_back(static_cast<float>(effTOF));
    x.push_back(static_cast<float>(effTRD));
    x.push_back(static_cast<float>(groups.useITS));
    x.push_back(static_cast<float>(effEMCal));
    x.push_back(static_cast<float>(effHMPID));
    x.push_back(static_cast<float>(groups.useCentrality));

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

/// PidOnnxInference: applies the FSE ONNX model to the features written by
/// pidFeatureExtractor.cxx and writes out per-track class probabilities.
///
/// All real work happens once, in init() - see runInference() above.
/// process() is intentionally empty; it exists only so this task has a
/// valid AOD-subscribing signature, as required by this repository's
/// conventions.
struct PidOnnxInference {
  Configurable<std::string> inputRootFile{"inputRootFile", "pid_features_data.root", "ROOT file produced by pidFeatureExtractor.cxx"};
  Configurable<std::string> inputTreeName{"inputTreeName", "pid_features", "Name of the TTree inside inputRootFile"};
  Configurable<std::string> outputPath{"outputPath", "pid_predictions", "Output file base name (no extension)"};
  Configurable<bool> exportCsv{"exportCsv", false, "Also write predictions to CSV alongside the ROOT output"};

  Configurable<bool> loadModelFromCcdb{"loadModelFromCcdb", true, "Load the ONNX model from CCDB (else from onnxFileNames as a local path)"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "CCDB URL"};
  Configurable<std::vector<std::string>> modelPathsCcdb{"modelPathsCcdb", std::vector<std::string>{"Users/YOURNAME/PidFeatureExtractor/model"}, "CCDB path to the model"};
  Configurable<int64_t> timestampCcdb{"timestampCcdb", -1, "CCDB query timestamp for the model, -1 = latest"};
  Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"pid_feature_model.onnx"}, "Local ONNX file path(s), used when loadModelFromCcdb is false"};

  Configurable<std::vector<double>> binsPtMl{"binsPtMl", std::vector<double>{-1., 9999.}, "pT bin edges for MlResponse (single bin = model isn't pT-binned)"};
  Configurable<int> nClassesMl{"nClassesMl", static_cast<int>(kNumClasses), "Number of model output classes"};

  Configurable<bool> useTPC{"useTPC", true, "Include TPC in inference. Default true (all detectors present); set false to force TPC excluded regardless of the data"};
  Configurable<bool> useTOF{"useTOF", true, "Include TOF in inference"};
  Configurable<bool> useTRD{"useTRD", true, "Include TRD in inference"};
  Configurable<bool> useITS{"useITS", true, "Include ITS in inference"};
  Configurable<bool> useEMCal{"useEMCal", true, "Include EMCal in inference"};
  Configurable<bool> useHMPID{"useHMPID", true, "Include HMPID in inference"};
  Configurable<bool> useCentrality{"useCentrality", true, "Include centrality in inference"};

  o2::ccdb::CcdbApi ccdbApi;
  o2::analysis::MlResponse<float> mlResponse;

  void init(InitContext&)
  {
    GroupToggles groups;
    groups.useTPC = useTPC.value;
    groups.useTOF = useTOF.value;
    groups.useTRD = useTRD.value;
    groups.useITS = useITS.value;
    groups.useEMCal = useEMCal.value;
    groups.useHMPID = useHMPID.value;
    groups.useCentrality = useCentrality.value;

    // Unused thresholds (CutNot everywhere) - this task always reports
    // all four probabilities rather than applying a selection cut, so
    // cutsMl/cutDirMl don't need to be user-configurable.
    static constexpr double kDefaultCutsMl[1][kNumClasses] = {{0., 0., 0., 0.}};
    LabeledArray<double> cutsMl{kDefaultCutsMl[0], 1, kNumClasses, {"pT bin 0"}, {"prob pi", "prob ka", "prob pr", "prob el"}};
    std::vector<int> cutDirMl{cuts_ml::CutNot, cuts_ml::CutNot, cuts_ml::CutNot, cuts_ml::CutNot};

    mlResponse.configure(binsPtMl.value, cutsMl, cutDirMl, static_cast<int8_t>(nClassesMl.value));
    if (loadModelFromCcdb.value) {
      ccdbApi.init(ccdbUrl.value);
      mlResponse.setModelPathsCCDB(onnxFileNames.value, ccdbApi, modelPathsCcdb.value, timestampCcdb.value);
    } else {
      mlResponse.setModelPathsLocal(onnxFileNames.value);
    }
    mlResponse.init();

    runInference(inputRootFile.value, inputTreeName.value, outputPath.value, exportCsv.value, mlResponse, groups);
  }

  /// Intentionally empty - all real work happens once in init(). Present
  /// only so this task has a valid, AOD-subscribing process() signature.
  void process(aod::Collisions const&) {}
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PidOnnxInference>(cfgc)};
}
