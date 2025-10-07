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
/// \file alice3PidEvaluation.cxx
///
/// \brief This task computes purity and efficiency from the OTF PID tables for multiple detectors.
///        Analyzes individual detectors (Tracker, TOF Inner, TOF Outer, RICH),
///        as well as combined detector performance using quadrature combination
///        of nSigma values.
///
/// \author Henrik Fribert TUM
/// \since  August 14, 2025
///

#include "ALICE3/DataModel/OTFPIDTrk.h"
#include "ALICE3/DataModel/OTFRICH.h"
#include "ALICE3/DataModel/OTFTOF.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonUtils/NameConf.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TVector3.h"

#include <array>
#include <memory>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;

struct Alice3PidEvaluation {

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  static constexpr float kInvalidNSigmaValue = 999.0f;

  Configurable<float> maxNSigmaForIdentification{"maxNSigmaForIdentification", 3.0f, "Maximum |nSigma| allowed for particle identification (closest-hypothesis rule)"};
  Configurable<int> numLogBins{"numLogBins", 200, "Number of logarithmic momentum bins"};
  Configurable<bool> useClosestHypothesisRule{"useClosestHypothesisRule", true, "Use closest-hypothesis rule: assign track to hypothesis with smallest |nSigma|"};
  Configurable<bool> useMinimalIdentification{"useMinimalIdentification", false, "Require that only one hypothesis is within the cutoff"};
  Configurable<bool> includeTrackerInCombined{"includeTrackerInCombined", true, "Include Tracker in combined analysis"};
  Configurable<bool> includeTofInnerInCombined{"includeTofInnerInCombined", true, "Include TOF Inner in combined analysis"};
  Configurable<bool> includeTofOuterInCombined{"includeTofOuterInCombined", true, "Include TOF Outer in combined analysis"};
  Configurable<bool> includeRichInCombined{"includeRichInCombined", false, "Include RICH in combined analysis"};

  std::vector<double> mLogBins;

  enum PidHypothesis { kElectron,
                       kMuon,
                       kPion,
                       kKaon,
                       kProton,
                       kDeuteron,
                       kTriton,
                       kHelium3,
                       kAlpha,
                       kCount };
  static constexpr std::array<int, PidHypothesis::kCount> kHypothesisPdg = {11, 13, 211, 321, 2212, 1000010020, 1000010030, 1000020030, 1000020040};
  static constexpr std::array<const char*, PidHypothesis::kCount> kHypothesisNames = {"Electron", "Muon", "Pion", "Kaon", "Proton", "Deuteron", "Triton", "Helium3", "Alpha"};

  struct DetectorHistograms {
    std::array<std::shared_ptr<TH1>, PidHypothesis::kCount> hTotalTrue;
    std::array<std::shared_ptr<TProfile>, PidHypothesis::kCount> hEfficiency;
    std::array<std::shared_ptr<TProfile>, PidHypothesis::kCount> hPurityAsHypothesis;
  };

  std::shared_ptr<TH2> hDetectorParticipation2D;

  DetectorHistograms trackerHists;
  DetectorHistograms tofInnerHists;
  DetectorHistograms tofOuterHists;
  DetectorHistograms richHists;
  DetectorHistograms combinedHists;
  DetectorHistograms combinedNoTrackerHists;
  void init(o2::framework::InitContext&)
  {
    LOG(info) << "Initializing multi-detector PID evaluation using closest-hypothesis rule";
    LOG(info) << "Maximum |nSigma| for identification: " << maxNSigmaForIdentification.value;
    LOG(info) << "Closest-hypothesis rule: " << (useClosestHypothesisRule.value ? "ENABLED" : "DISABLED");
    LOG(info) << "Require unique identification: " << (useMinimalIdentification.value ? "YES" : "NO");
    LOG(info) << "Combined analysis includes: "
              << (includeTrackerInCombined.value ? "Tracker " : "")
              << (includeTofInnerInCombined.value ? "TOF_Inner " : "")
              << (includeTofOuterInCombined.value ? "TOF_Outer " : "")
              << (includeRichInCombined.value ? "RICH " : "");

    mLogBins.clear();
    double pMin = 0.05;
    double pMax = 10;
    double logMin = std::log10(pMin);
    double logMax = std::log10(pMax);
    double dLog = (logMax - logMin) / numLogBins.value;
    for (int i = 0; i <= numLogBins.value; ++i) {
      mLogBins.push_back(std::pow(10, logMin + i * dLog));
    }
    const AxisSpec axisMomentum{mLogBins, "#it{p} (GeV/#it{c})"};

    auto createDetectorHistograms = [&](DetectorHistograms& detHists, const std::string& detectorName) {
      for (int trueIdx = 0; trueIdx < PidHypothesis::kCount; ++trueIdx) {
        const auto& trueName = kHypothesisNames[trueIdx];
        detHists.hTotalTrue[trueIdx] = histos.add<TH1>(Form("%s/hTotalTrue%s", detectorName.c_str(), trueName),
                                                       Form("%s: Total True %s; #it{p} (GeV/#it{c})", detectorName.c_str(), trueName),
                                                       kTH1F, {axisMomentum});
        detHists.hEfficiency[trueIdx] = histos.add<TProfile>(Form("%s/hEfficiency%s", detectorName.c_str(), trueName),
                                                             Form("%s: PID Efficiency for %s; #it{p} (GeV/#it{c}); Efficiency", detectorName.c_str(), trueName),
                                                             kTProfile, {axisMomentum});
      }

      for (int hypIdx = 0; hypIdx < PidHypothesis::kCount; ++hypIdx) {
        const auto& hypName = kHypothesisNames[hypIdx];
        detHists.hPurityAsHypothesis[hypIdx] = histos.add<TProfile>(Form("%s/hPurityAs%s", detectorName.c_str(), hypName),
                                                                    Form("%s: Purity when selecting as %s; #it{p} (GeV/#it{c}); Purity", detectorName.c_str(), hypName),
                                                                    kTProfile, {axisMomentum});
      }
    };

    createDetectorHistograms(trackerHists, "Tracker");
    createDetectorHistograms(tofInnerHists, "TOF_Inner");
    createDetectorHistograms(tofOuterHists, "TOF_Outer");
    createDetectorHistograms(richHists, "RICH");
    createDetectorHistograms(combinedHists, "Combined");
    createDetectorHistograms(combinedNoTrackerHists, "Combined_NoTracker");

    const AxisSpec axisDetectorCount{5, -0.5, 4.5, "Number of detectors"};
    hDetectorParticipation2D = histos.add<TH2>("Combined/hDetectorParticipation2D",
                                               "Detector participation vs momentum; #it{p} (GeV/#it{c}); Number of detectors",
                                               kTH2F, {axisMomentum, axisDetectorCount});
  }

  void process(soa::Join<aod::Tracks, aod::TracksCov, aod::McTrackLabels, aod::UpgradeTrkPids, aod::UpgradeTrkPidSignals, aod::UpgradeTofs, aod::UpgradeRichs> const& tracks,
               aod::McParticles const& /*mcParticles*/)
  {

    auto isValidNSigma = [](float nSigma) -> bool {
      return (nSigma < kInvalidNSigmaValue && nSigma > -kInvalidNSigmaValue);
    };

    auto computeCombinedNSigma = [&](const std::vector<std::array<float, PidHypothesis::kCount>>& detectorNSigmas, float p) -> std::array<float, PidHypothesis::kCount> {
      std::array<float, PidHypothesis::kCount> combinedNSigma;
      int totalValidDetectors = 0;
      for (const auto& detNSigma : detectorNSigmas) {
        bool detectorHasValidMeasurement = false;
        for (int hypIdx = 0; hypIdx < PidHypothesis::kCount; ++hypIdx) {
          if (isValidNSigma(detNSigma[hypIdx])) {
            detectorHasValidMeasurement = true;
            break;
          }
        }
        if (detectorHasValidMeasurement) {
          totalValidDetectors++;
        }
      }

      for (int hypIdx = 0; hypIdx < PidHypothesis::kCount; ++hypIdx) {
        float sumSquares = 0.0f;
        int validDetectors = 0;

        for (const auto& detNSigma : detectorNSigmas) {
          if (isValidNSigma(detNSigma[hypIdx])) {
            sumSquares += detNSigma[hypIdx] * detNSigma[hypIdx];
            validDetectors++;
          }
        }

        if (validDetectors > 0) {
          combinedNSigma[hypIdx] = std::sqrt(sumSquares);
        } else {
          combinedNSigma[hypIdx] = kInvalidNSigmaValue;
        }
      }
      if (totalValidDetectors > 0) {
        hDetectorParticipation2D->Fill(p, totalValidDetectors);
      }

      return combinedNSigma;
    };

    auto analyzeDetector = [&](DetectorHistograms& detHists, const std::array<float, PidHypothesis::kCount>& nSigmaValues,
                               int trueParticleIndex, float p) {
      detHists.hTotalTrue[trueParticleIndex]->Fill(p);
      bool hasValidNSigma = false;
      for (int i = 0; i < PidHypothesis::kCount; ++i) {
        if (isValidNSigma(nSigmaValues[i])) {
          hasValidNSigma = true;
          break;
        }
      }
      if (!hasValidNSigma) {
        return;
      }

      bool correctlyIdentified = false;
      int selectedHypothesis = -1;

      if (useClosestHypothesisRule.value) {
        float minAbsNSigma = kInvalidNSigmaValue;
        int bestHypothesis = -1;
        int validHypothesesCount = 0;

        for (int hypIdx = 0; hypIdx < PidHypothesis::kCount; ++hypIdx) {
          if (isValidNSigma(nSigmaValues[hypIdx])) {
            float absNSigma = std::fabs(nSigmaValues[hypIdx]);
            if (absNSigma < minAbsNSigma) {
              minAbsNSigma = absNSigma;
              bestHypothesis = hypIdx;
            }
            if (absNSigma < maxNSigmaForIdentification.value) {
              validHypothesesCount++;
            }
          }
        }

        if (bestHypothesis >= 0 && minAbsNSigma < maxNSigmaForIdentification.value) {
          if (useMinimalIdentification.value && validHypothesesCount > 1) {
            selectedHypothesis = -1;
          } else {
            selectedHypothesis = bestHypothesis;
          }
        }
        correctlyIdentified = (selectedHypothesis == trueParticleIndex);
      } else {
        correctlyIdentified = (std::fabs(nSigmaValues[trueParticleIndex]) < maxNSigmaForIdentification.value);
        for (int hypIdx = 0; hypIdx < PidHypothesis::kCount; ++hypIdx) {
          if (std::fabs(nSigmaValues[hypIdx]) < maxNSigmaForIdentification.value) {
            bool isCorrect = (hypIdx == trueParticleIndex);
            detHists.hPurityAsHypothesis[hypIdx]->Fill(p, isCorrect ? 1.0 : 0.0);
          }
        }
      }

      detHists.hEfficiency[trueParticleIndex]->Fill(p, correctlyIdentified ? 1.0 : 0.0);

      if (useClosestHypothesisRule.value && selectedHypothesis >= 0) {
        bool isCorrect = (selectedHypothesis == trueParticleIndex);
        detHists.hPurityAsHypothesis[selectedHypothesis]->Fill(p, isCorrect ? 1.0 : 0.0);
      }
    };

    for (const auto& track : tracks) {
      if (!track.has_mcParticle()) {
        continue;
      }

      const auto& mcParticle = track.mcParticle();
      const float p = mcParticle.p();
      const int truePdg = std::abs(mcParticle.pdgCode());

      int trueParticleIndex = -1;
      for (int i = 0; i < PidHypothesis::kCount; ++i) {
        if (kHypothesisPdg[i] == truePdg) {
          trueParticleIndex = i;
          break;
        }
      }
      if (trueParticleIndex == -1) {
        continue;
      }

      std::array<float, PidHypothesis::kCount> trackerNSigma;
      trackerNSigma[kElectron] = track.nSigmaTrkEl();
      trackerNSigma[kMuon] = track.nSigmaTrkMu();
      trackerNSigma[kPion] = track.nSigmaTrkPi();
      trackerNSigma[kKaon] = track.nSigmaTrkKa();
      trackerNSigma[kProton] = track.nSigmaTrkPr();
      trackerNSigma[kDeuteron] = track.nSigmaTrkDe();
      trackerNSigma[kTriton] = track.nSigmaTrkTr();
      trackerNSigma[kHelium3] = track.nSigmaTrkHe();
      trackerNSigma[kAlpha] = track.nSigmaTrkAl();

      analyzeDetector(trackerHists, trackerNSigma, trueParticleIndex, p);

      std::array<float, PidHypothesis::kCount> tofInnerNSigma;
      tofInnerNSigma[kElectron] = track.nSigmaElectronInnerTOF();
      tofInnerNSigma[kMuon] = track.nSigmaMuonInnerTOF();
      tofInnerNSigma[kPion] = track.nSigmaPionInnerTOF();
      tofInnerNSigma[kKaon] = track.nSigmaKaonInnerTOF();
      tofInnerNSigma[kProton] = track.nSigmaProtonInnerTOF();
      tofInnerNSigma[kDeuteron] = track.nSigmaDeuteronInnerTOF();
      tofInnerNSigma[kTriton] = track.nSigmaTritonInnerTOF();
      tofInnerNSigma[kHelium3] = track.nSigmaHelium3InnerTOF();
      tofInnerNSigma[kAlpha] = track.nSigmaAlphaInnerTOF();

      analyzeDetector(tofInnerHists, tofInnerNSigma, trueParticleIndex, p);

      std::array<float, PidHypothesis::kCount> tofOuterNSigma;
      tofOuterNSigma[kElectron] = track.nSigmaElectronOuterTOF();
      tofOuterNSigma[kMuon] = track.nSigmaMuonOuterTOF();
      tofOuterNSigma[kPion] = track.nSigmaPionOuterTOF();
      tofOuterNSigma[kKaon] = track.nSigmaKaonOuterTOF();
      tofOuterNSigma[kProton] = track.nSigmaProtonOuterTOF();
      tofOuterNSigma[kDeuteron] = track.nSigmaDeuteronOuterTOF();
      tofOuterNSigma[kTriton] = track.nSigmaTritonOuterTOF();
      tofOuterNSigma[kHelium3] = track.nSigmaHelium3OuterTOF();
      tofOuterNSigma[kAlpha] = track.nSigmaAlphaOuterTOF();

      analyzeDetector(tofOuterHists, tofOuterNSigma, trueParticleIndex, p);

      std::array<float, PidHypothesis::kCount> richNSigma;
      richNSigma[kElectron] = track.nSigmaElectronRich();
      richNSigma[kMuon] = track.nSigmaMuonRich();
      richNSigma[kPion] = track.nSigmaPionRich();
      richNSigma[kKaon] = track.nSigmaKaonRich();
      richNSigma[kProton] = track.nSigmaProtonRich();
      richNSigma[kDeuteron] = track.nSigmaDeuteronRich();
      richNSigma[kTriton] = track.nSigmaTritonRich();
      richNSigma[kHelium3] = track.nSigmaHelium3Rich();
      richNSigma[kAlpha] = track.nSigmaAlphaRich();

      analyzeDetector(richHists, richNSigma, trueParticleIndex, p);

      std::vector<std::array<float, PidHypothesis::kCount>> allDetectorNSigmas;

      if (includeTrackerInCombined.value) {
        allDetectorNSigmas.push_back(trackerNSigma);
      }
      if (includeTofInnerInCombined.value) {
        allDetectorNSigmas.push_back(tofInnerNSigma);
      }
      if (includeTofOuterInCombined.value) {
        allDetectorNSigmas.push_back(tofOuterNSigma);
      }
      if (includeRichInCombined.value) {
        allDetectorNSigmas.push_back(richNSigma);
      }
      if (!allDetectorNSigmas.empty()) {
        std::array<float, PidHypothesis::kCount> combinedNSigma = computeCombinedNSigma(allDetectorNSigmas, p);
        analyzeDetector(combinedHists, combinedNSigma, trueParticleIndex, p);
      }

      std::vector<std::array<float, PidHypothesis::kCount>> noTrackerDetectorNSigmas;

      if (includeTofInnerInCombined.value) {
        noTrackerDetectorNSigmas.push_back(tofInnerNSigma);
      }
      if (includeTofOuterInCombined.value) {
        noTrackerDetectorNSigmas.push_back(tofOuterNSigma);
      }
      if (includeRichInCombined.value) {
        noTrackerDetectorNSigmas.push_back(richNSigma);
      }
      if (!noTrackerDetectorNSigmas.empty()) {
        std::array<float, PidHypothesis::kCount> combinedNoTrackerNSigma = computeCombinedNSigma(noTrackerDetectorNSigmas, p);
        analyzeDetector(combinedNoTrackerHists, combinedNoTrackerNSigma, trueParticleIndex, p);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<Alice3PidEvaluation>(cfgc)};
}
