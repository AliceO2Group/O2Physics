// Copyright 2019-2023 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file   applyMlSelection.cxx
/// \brief  Showcase usage of trained ML model exported to ONNX for candidate selection in analysis
///
/// \author Fabio Catalano <fabio.catalano@cern.ch>, CERN

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "Tools/ML/MlResponse.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;

static constexpr double defaultCutsMl[1][3] = {{0.5, 0.5, 0.5}};

struct applyMlSelection {
  // Analysis configuration
  Configurable<double> ptCandMin{"ptCandMin", 1., "Lower bound of candidate pT"};
  Configurable<double> ptCandMax{"ptCandMax", 36., "Upper bound of candidate pT"};
  Configurable<double> yCandRecoMax{"yCandRecoMax", 0.8, "max. cand. rapidity"};
  // ML inference
  Configurable<std::vector<double>> binsPtMl{"binsPtMl", std::vector<double>{1., 36.}, "pT bin limits for ML application"};
  Configurable<std::vector<int>> cutDirMl{"cutDirMl", std::vector<int>{cuts_ml::CutSmaller, cuts_ml::CutNot, cuts_ml::CutNot}, "Whether to reject score values greater or smaller than the threshold"};
  Configurable<LabeledArray<double>> cutsMl{"cutsMl", {defaultCutsMl[0], 1, 3, {"pT bin 0"}, {"score prompt", "score non-prompt", "score bkg"}}, "ML selections per pT bin"};
  Configurable<int8_t> nClassesMl{"nClassesMl", (int8_t)3, "Number of classes in ML model"};
  // Model file names
  Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"model_onnx.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
  // Bonus: CCDB configuration (needed for ML application on the GRID)
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> modelPathsCCDB{"modelPathsCCDB", "Users/f/fcatalan/O2AT3/MlInference", "Path on CCDB"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};

  Filter filterDsFlag = (o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(BIT(aod::hf_cand_3prong::DecayType::DsToKKPi))) != static_cast<uint8_t>(0);

  HfHelper hfHelper;
  o2::ccdb::CcdbApi ccdbApi;
  int nCandidates = 0;

  // Add objects needed for ML inference
  std::vector<float> outputMl = {};

  // Add histograms for other BDT scores and for distributions after selections
  HistogramRegistry registry{
    "registry",
    {{"hMassBeforeSel", "Ds candidates before selection;inv. mass (KK#pi) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{100, 1.77, 2.17}}}},
     {"hPromptScoreBeforeSel", "Prompt score before selection;BDT first score;entries", {HistType::kTH1F, {{100, 0., 1.}}}}}};

  void init(InitContext const&)
  {
    // Add histograms vs pT (only for selected candidates)
    auto vbins = (std::vector<double>)binsPtMl;
    registry.add("hMassAfterSelVsPt", "Ds candidates;inv. mass (KK#pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{100, 1.77, 2.17}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPromptScoreAfterSelVsPt", "Prompt score after selection;BDT first score;entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    // Configure and initialise the ML class

    // Bonus: retrieve the model from CCDB (needed for ML application on the GRID)
  }

  void process(soa::Filtered<aod::HfCand3Prong> const& candidates)
  {
    // Looping over Ds candidates
    for (const auto& candidate : candidates) {

      nCandidates++;
      if (nCandidates % 1000 == 0) {
        LOG(info) << "Candidates processed: " << nCandidates;
      }

      auto candpT = candidate.pt();

      // Check that the candidate pT is within the analysis range
      if (candpT < ptCandMin || candpT > ptCandMax) {
        continue;
      }
      // Check that the candidate rapidity is within the analysis range
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yDs(candidate)) > yCandRecoMax) {
        continue;
      }

      // Fill mass histogram before ML selections
      registry.fill(HIST("hMassBeforeSel"), hfHelper.invMassDsToPiKK(candidate));
      registry.fill(HIST("hMassBeforeSel"), hfHelper.invMassDsToKKPi(candidate));

      // Perform ML selections for one mass hypothesis (Ds -> PhiPi -> PiKK)
      std::vector<float> inputFeaturesPiKK{candidate.cpa(),
                                           candidate.cpaXY(),
                                           candidate.decayLength(),
                                           candidate.decayLengthXY(),
                                           static_cast<float>(hfHelper.deltaMassPhiDsToPiKK(candidate)),
                                           candidate.impactParameterXY(),
                                           static_cast<float>(hfHelper.cos3PiKDsToPiKK(candidate)),
                                           candidate.maxNormalisedDeltaIP()};

      // Retrieve model output and selection outcome

      // Fill BDT score histograms before selection
      registry.fill(HIST("hPromptScoreBeforeSel"), outputMl[0]);

      // Fill histograms for selected candidates
      bool isSelectedMlPiKK = true;
      if (isSelectedMlPiKK) {
        registry.fill(HIST("hMassAfterSelVsPt"), hfHelper.invMassDsToPiKK(candidate), candidate.pt());
        registry.fill(HIST("hPromptScoreAfterSelVsPt"), outputMl[0], candidate.pt());
      }

      outputMl.clear(); // not necessary in this case but for good measure

      // Perform ML selections for other mass hypothesis (Ds -> PhiPi -> KKPi)
      std::vector<float> inputFeaturesKKPi{candidate.cpa(),
                                           candidate.cpaXY(),
                                           candidate.decayLength(),
                                           candidate.decayLengthXY(),
                                           static_cast<float>(hfHelper.deltaMassPhiDsToKKPi(candidate)),
                                           candidate.impactParameterXY(),
                                           static_cast<float>(hfHelper.cos3PiKDsToKKPi(candidate)),
                                           candidate.maxNormalisedDeltaIP()};

      // Retrieve model output and selection outcome

      // Fill BDT score histograms before selection

      // Fill histograms for selected candidates

      outputMl.clear(); // not necessary in this case but for good measure
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<applyMlSelection>(cfgc)};
}
