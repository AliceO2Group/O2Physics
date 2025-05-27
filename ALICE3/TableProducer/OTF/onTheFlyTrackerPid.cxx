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
/// \file onTheFlyTrackerPid.cxx
///
/// \brief This task produces the PID information that can be obtained from the tracker layers (i.e. cluster size and ToT).
///        It currently contemplates 5 particle types: electrons, muons, pions, kaons and protons.
///
/// \author Berkin Ulukutlu TUM
/// \author Henrik Fribert TUM
/// \author Nicolò Jacazio Università del Piemonte Orientale
/// \since  May 22, 2025
///

#include <utility>
#include <map>
#include <string>
#include <algorithm>
#include <vector>

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "ALICE3/Core/TrackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "CommonUtils/NameConf.h"
#include "CCDB/CcdbApi.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "CommonConstants/GeomConstants.h"
#include "CommonConstants/PhysicsConstants.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TH2F.h"
#include "TVector3.h"
#include "TString.h"
#include "ALICE3/DataModel/OTFRICH.h"
#include "DetectorsVertexing/HelixHelper.h"
#include "TableHelper.h"
#include "ALICE3/Core/DelphesO2TrackSmearer.h"
#include "ALICE3/DataModel/OTFPIDTrk.h"

using namespace o2;
using namespace o2::framework;

struct OnTheFlyTrackerPid {
  Produces<aod::UpgradeTrkPidSignals> tableUpgradeTrkPidSignals;
  Produces<aod::UpgradeTrkPids> tableUpgradeTrkPids;

  // necessary for particle charges
  Service<o2::framework::O2DatabasePDG> pdg;

  static constexpr int kMaxBarrelLayers = 8;
  static constexpr int kMaxForwardLayers = 9;

  struct : ConfigurableGroup {
    Configurable<std::string> efficiencyFormula{"efficiencyFormula", "1.0/(1.0+exp(-(x-0.01)/0.2))", "ROOT TF1 formula for efficiency"};
    Configurable<std::string> landauFormula{"landauFormula", "TMath::Landau(x, 1, 1, true)", "ROOT TF1 formula for Landau distribution (e.g. ToT response)"};
    Configurable<int> averageMethod{"averageMethod", 0, "Method to average the ToT and cluster size. 0: truncated mean"};
  } simConfig;

  TF1* mEfficiency = nullptr;
  static constexpr int kEtaBins = 50;
  static constexpr float kEtaMin = -2.5;
  static constexpr float kEtaMax = 2.5;
  static constexpr int kPtBins = 200;
  static constexpr float kPtMin = 0.0;
  static constexpr float kPtMax = 20.0;

  std::array<std::array<TF1*, kPtBins>, kEtaBins> mElossPi;

  void init(o2::framework::InitContext&)
  {

    for (int i = 0; i < kEtaBins; i++) {
      for (int j = 0; j < kPtBins; j++) {
        mElossPi[i][j] = new TF1(Form("mElossPi_%d_%d", i, j), simConfig.landauFormula.value.c_str(), 0, 20);
      }
    }
    mEfficiency = new TF1("mEfficiency", simConfig.efficiencyFormula.value.c_str(), 0, 20);
  }

  void process(soa::Join<aod::Collisions, aod::McCollisionLabels>::iterator const&,
               soa::Join<aod::Tracks, aod::TracksCov, aod::McTrackLabels> const& tracks,
               aod::McParticles const&,
               aod::McCollisions const&)
  {
    std::array<float, kMaxBarrelLayers> timeOverThresholdBarrel;
    std::array<float, kMaxBarrelLayers> clusterSizeBarrel;
    // std::array<float, kMaxForwardLayers> timeOverThresholdForward;
    // std::array<float, kMaxForwardLayers> clusterSizeForward;

    auto noSignalTrack = [&]() {
      tableUpgradeTrkPidSignals(0.f, 0.f);          // no PID information
      tableUpgradeTrkPids(0.f, 0.f, 0.f, 0.f, 0.f); // no PID information
    };

    for (const auto& track : tracks) {
      if (!track.has_mcParticle()) {
        noSignalTrack();
        continue;
      }
      const auto& mcParticle = track.mcParticle();
      const auto& pdgInfo = pdg->GetParticle(mcParticle.pdgCode());
      if (!pdgInfo) {
        LOG(warning) << "PDG code " << mcParticle.pdgCode() << " not found in the database";
        noSignalTrack();
        continue;
      }
      const float pt = mcParticle.pt();
      const float eta = mcParticle.eta();

      const int binnedPt = static_cast<int>((pt - kPtMin) / kPtBins);
      const int binnedEta = static_cast<int>((eta - kEtaMin) / kEtaBins);
      if (binnedPt < 0 || binnedPt >= kPtBins || binnedEta < 0 || binnedEta >= kEtaBins) {
        noSignalTrack();
        continue;
      }
      for (int i = 0; i < kMaxBarrelLayers; i++) {
        timeOverThresholdBarrel[i] = -1;
        clusterSizeBarrel[i] = -1;

        // Check if layer is efficient
        if (mEfficiency->Eval(pt) > gRandom->Uniform(0, 1)) {
          timeOverThresholdBarrel[i] = mElossPi[binnedEta][binnedPt]->GetRandom(); // Simulate ToT
          clusterSizeBarrel[i] = mElossPi[binnedEta][binnedPt]->GetRandom();       // Simulate cluster size
        }
      }

      // Now we do the average
      switch (simConfig.averageMethod) {
        case 0: { // truncated mean
          float meanToT = 0;
          float meanClusterSize = 0;
          // Order them by ToT
          std::sort(timeOverThresholdBarrel.begin(), timeOverThresholdBarrel.end());
          std::sort(clusterSizeBarrel.begin(), clusterSizeBarrel.end());
          static constexpr int kTruncatedMean = 5;
          // Take the mean of the first 5 values
          for (int i = 0; i < kTruncatedMean; i++) {
            meanToT += timeOverThresholdBarrel[i];
            meanClusterSize += clusterSizeBarrel[i];
          }
          meanToT /= kTruncatedMean;
          meanClusterSize /= kTruncatedMean;
          // Fill the table
          tableUpgradeTrkPidSignals(meanToT, meanClusterSize);
        } break;

        default:
          LOG(fatal) << "Unknown average method " << simConfig.averageMethod;
          break;
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<OnTheFlyTrackerPid>(cfgc)}; }
