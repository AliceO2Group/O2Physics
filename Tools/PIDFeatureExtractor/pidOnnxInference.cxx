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
///        testing) over the tables produced by pidFeatureExtractor.cxx, and
///        write out per-track class probabilities.
///
///        Uses Tools/ML/MlResponse - O2Physics's generic ONNX/CCDB inference
///        wrapper (not PID-specific, used across several PWG groups) -
///        rather than hand-rolled ONNX/CCDB loading, so this task stays
///        small. It depends only on that shared ML infrastructure and on
///        pidFeatureExtractor.h in this same folder - no other analysis
///        task.
///
/// \author Robert Forynski

#include "pidFeatureExtractor.h"
//
#include "Tools/ML/MlResponse.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/runDataProcessing.h>

#include <cstdint>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;

// ============================================================================
// OUTPUT TABLE
// ----------------------------------------------------------------------------
// Declared directly here since nothing else needs to depend on it.
// ============================================================================
namespace o2::aod
{
namespace pidpred
{
DECLARE_SOA_COLUMN(MlProbPi, mlProbPi, float);   //!
DECLARE_SOA_COLUMN(MlProbKa, mlProbKa, float);   //!
DECLARE_SOA_COLUMN(MlProbPr, mlProbPr, float);   //!
DECLARE_SOA_COLUMN(MlProbEl, mlProbEl, float);   //!
DECLARE_SOA_COLUMN(MlPredictedClass, mlPredictedClass, int); //! argmax of the 4 probs above: 0=pi,1=ka,2=pr,3=el
} // namespace pidpred

DECLARE_SOA_TABLE(PidMlPredictions, "AOD", "PIDMLPRED", //!
                  o2::soa::Index<>,
                  pidpred::MlProbPi, pidpred::MlProbKa, pidpred::MlProbPr, pidpred::MlProbEl,
                  pidpred::MlPredictedClass);
} // namespace o2::aod

namespace
{
constexpr int kNumClasses = 4; // pi, ka, pr, el - fixed order throughout, matches the paper's model

// ----------------------------------------------------------------------------
// Feature order fed to the model.
// ----------------------------------------------------------------------------
// THIS MUST MATCH YOUR TRAINING SCRIPT'S COLUMN ORDER EXACTLY - a silent
// mismatch here is the single most likely way this task produces wrong
// predictions without any error. The list below is a reasonable default
// (every reconstructed feature in PidFeaturesData/PidFeaturesMc except
// vz/centFT0C/sign/trackType and the Bayesian columns, which are a
// comparison baseline, not a model input) - it is NOT verified against your
// actual training code. Reorder, add, or drop entries to match exactly
// before trusting the output.
//
// If your ONNX export takes features and mask as two separate input
// tensors rather than one concatenated vector, split this function
// accordingly.

/// itsClusterSizes packs 7 ITS layers into 4 bits each; a derived cluster
/// count is a far more sensible model input than the raw packed value.
/// Small and duplicated here rather than shared with pidFeatureExtractor.cxx
/// so this task stays self-contained.
template <typename TRow>
int getItsNClusters(TRow const& row)
{
  auto v = static_cast<uint32_t>(row.itsClusterSizes());
  int n = 0;
  for (int layer = 0; layer < 7; layer++) {
    if ((v >> (layer * 4)) & 0xF) {
      n++;
    }
  }
  return n;
}

template <typename TRow>
std::vector<float> buildModelInput(TRow const& row)
{
  std::vector<float> x;
  x.reserve(38 + 7);

  // Kinematics
  x.push_back(row.p());
  x.push_back(row.pt());
  x.push_back(row.px());
  x.push_back(row.py());
  x.push_back(row.pz());
  x.push_back(row.eta());
  x.push_back(row.phi());
  // Impact parameters
  x.push_back(row.dcaXY());
  x.push_back(row.dcaZ());
  // TPC
  x.push_back(static_cast<float>(row.hasTPC()));
  x.push_back(row.tpcSignal());
  x.push_back(row.tpcNSigmaPi());
  x.push_back(row.tpcNSigmaKa());
  x.push_back(row.tpcNSigmaPr());
  x.push_back(row.tpcNSigmaEl());
  x.push_back(static_cast<float>(row.tpcNClsFound()));
  x.push_back(row.tpcChi2NCl());
  // TOF
  x.push_back(static_cast<float>(row.hasTOF()));
  x.push_back(row.tofMass());
  x.push_back(row.beta());
  x.push_back(row.tofNSigmaPi());
  x.push_back(row.tofNSigmaKa());
  x.push_back(row.tofNSigmaPr());
  x.push_back(row.tofNSigmaEl());
  // TRD
  x.push_back(static_cast<float>(row.hasTRD()));
  x.push_back(row.trdSignal());
  x.push_back(row.trdChi2());
  x.push_back(static_cast<float>(row.trdPattern()));
  // ITS
  x.push_back(static_cast<float>(getItsNClusters(row)));
  x.push_back(row.itsChi2NCl());
  // EMCal
  x.push_back(static_cast<float>(row.hasEMCal()));
  x.push_back(row.trackEtaEmcal());
  x.push_back(row.trackPhiEmcal());
  // HMPID
  x.push_back(static_cast<float>(row.hasHMPID()));
  x.push_back(row.hmpidSignal());
  x.push_back(row.hmpidQMip());
  x.push_back(static_cast<float>(row.hmpidNPhotons()));
  x.push_back(static_cast<float>(row.hmpidClusSize()));
  x.push_back(row.hmpidMom());

  // 7-length group mask: TPC, TOF, TRD, ITS, EMCal, HMPID, centrality.
  // ITS is assumed always present (global track requires it); centrality is
  // assumed always present. Both are real assumptions, not derived facts -
  // adjust if your training data ever has either group absent.
  x.push_back(static_cast<float>(row.hasTPC()));
  x.push_back(static_cast<float>(row.hasTOF()));
  x.push_back(static_cast<float>(row.hasTRD()));
  x.push_back(1.f); // ITS
  x.push_back(static_cast<float>(row.hasEMCal()));
  x.push_back(static_cast<float>(row.hasHMPID()));
  x.push_back(1.f); // centrality

  return x;
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

/// PidOnnxInference: applies the FSE ONNX model to features produced by
/// pidFeatureExtractor.cxx and writes out per-track class probabilities.
///
/// Model loading (local file or CCDB) and ONNX execution are entirely
/// handled by o2::analysis::MlResponse - this task only builds the input
/// feature vector and reads back the output. A single global pT bin is used
/// by default since the paper's model isn't pT-binned; MlResponse's usual
/// per-class cut mechanism is disabled (cutDirMl = CutNot for every class,
/// confirmed value 2 in the real o2::cuts_ml enum) - this task always
/// reports all four probabilities rather than applying a pass/fail cut.
///
/// Two things MlResponse enforces that are worth knowing before debugging a
/// failure here: getModelOutput() calls LOG(fatal) if the input vector's
/// length doesn't match the ONNX model's declared input node count (unless
/// that node is a dynamic axis), and separately if pt lands outside
/// binsPtMl's range entirely (see the comment on binsPtMl below).
struct PidOnnxInference {
  Produces<aod::PidMlPredictions> pidMlPredictions;

  // --- model location ----------------------------------------------------------
  Configurable<bool> loadModelFromCcdb{"loadModelFromCcdb", true, "Load the ONNX model from CCDB (else from onnxFileNames as a local path)"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "CCDB URL"};
  Configurable<std::vector<std::string>> modelPathsCcdb{"modelPathsCcdb", std::vector<std::string>{"Users/YOURNAME/PidFeatureExtractor/model"}, "CCDB path to the model"};
  Configurable<int64_t> timestampCcdb{"timestampCcdb", -1, "CCDB query timestamp for the model, -1 = latest"};
  Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"pid_feature_model.onnx"}, "Local ONNX file path(s), used when loadModelFromCcdb is false"};

  // --- MlResponse plumbing: a single pT bin, no selection cut applied --------
  // Lower edge is -1 (not 0): MlResponse::findBin() rejects value < front()
  // as out-of-range (fatal in getModelOutput), and track.pt() could in
  // principle be exactly 0 - keeping the edge below any physical pT avoids
  // that boundary case entirely.
  Configurable<std::vector<double>> binsPtMl{"binsPtMl", std::vector<double>{-1., 9999.}, "pT bin edges for MlResponse (single bin = model isn't pT-binned)"};
  Configurable<std::vector<int>> cutDirMl{"cutDirMl", std::vector<int>{cuts_ml::CutNot, cuts_ml::CutNot, cuts_ml::CutNot, cuts_ml::CutNot}, "Per-class cut direction; CutNot = always accept, this task doesn't select"};

  static constexpr double kDefaultCutsMl[1][kNumClasses] = {{0., 0., 0., 0.}};
  Configurable<LabeledArray<double>> cutsMl{"cutsMl", {kDefaultCutsMl[0], 1, kNumClasses, {"pT bin 0"}, {"prob pi", "prob ka", "prob pr", "prob el"}}, "Unused thresholds (CutNot everywhere) - required by MlResponse's interface"};
  Configurable<int8_t> nClassesMl{"nClassesMl", static_cast<int8_t>(kNumClasses), "Number of model output classes"};

  o2::ccdb::CcdbApi ccdbApi;
  o2::analysis::MlResponse<float> mlResponse;
  std::vector<float> mlOutput;

  void init(InitContext const&)
  {
    mlResponse.configure(binsPtMl, cutsMl, cutDirMl, nClassesMl);
    if (loadModelFromCcdb) {
      ccdbApi.init(ccdbUrl.value);
      mlResponse.setModelPathsCCDB(onnxFileNames, ccdbApi, modelPathsCcdb.value, timestampCcdb.value);
    } else {
      mlResponse.setModelPathsLocal(onnxFileNames);
    }
    mlResponse.init();
  }

  template <typename TTable>
  void runInference(TTable const& rows)
  {
    for (auto const& row : rows) {
      auto x = buildModelInput(row);
      mlResponse.isSelectedMl(x, row.pt(), mlOutput); // return value (selection) unused; mlOutput carries the 4 raw scores
      pidMlPredictions(mlOutput[0], mlOutput[1], mlOutput[2], mlOutput[3], argmax4(mlOutput));
      mlOutput.clear();
    }
  }

  void processData(aod::PidFeaturesData const& rows)
  {
    runInference(rows);
  }
  PROCESS_SWITCH(PidOnnxInference, processData, "Run inference on PidFeaturesData", true);

  void processMc(aod::PidFeaturesMc const& rows)
  {
    runInference(rows);
  }
  PROCESS_SWITCH(PidOnnxInference, processMc, "Run inference on PidFeaturesMc", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PidOnnxInference>(cfgc)};
}
