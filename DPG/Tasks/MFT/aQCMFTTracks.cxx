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
/// \file   aQCMFTTracks.cxx
/// \brief  This task runs over AO2Ds and fills some basic objects
///         needed in the asynchronous QC checks on MFT tracks.
/// \author David Grund
/// \since

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/LHCConstants.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/DataTypes.h"
#include "Framework/TimingInfo.h"
#include "Framework/runDataProcessing.h"

#include <TH1F.h>
#include <TH2F.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::aod;

struct CheckMFT {
  HistogramRegistry registry{"registry",
                             {// 2d histograms
                              {"mMFTTrackEtaPhi_5_MinClusters", "Track #eta , #phi (NCls >= 5); #eta; #phi", {HistType::kTH2F, {{50, -4, -2}, {100, -3.2, 3.2}}}},
                              {"mMFTTrackXY_5_MinClusters", "Track Position (NCls >= 5); x; y", {HistType::kTH2F, {{320, -16, 16}, {320, -16, 16}}}},
                              {"mMFTTrackEtaPhi_7_MinClusters", "Track #eta , #phi (NCls >= 7); #eta; #phi", {HistType::kTH2F, {{50, -4, -2}, {100, -3.2, 3.2}}}},
                              {"mMFTTrackXY_7_MinClusters", "Track Position (NCls >= 7); x; y", {HistType::kTH2F, {{320, -16, 16}, {320, -16, 16}}}},
                              {"mMFTTrackEtaPhi_8_MinClusters", "Track #eta , #phi (NCls >= 8); #eta; #phi", {HistType::kTH2F, {{50, -4, -2}, {100, -3.2, 3.2}}}},
                              {"mMFTTrackXY_8_MinClusters", "Track Position (NCls >= 8); x; y", {HistType::kTH2F, {{320, -16, 16}, {320, -16, 16}}}},
                              // 1d histograms
                              {"mMFTTrackEta", "Track #eta; #eta; # entries", {HistType::kTH1F, {{50, -4, -2}}}},
                              {"mMFTTrackNumberOfClusters", "Number Of Clusters Per Track; # clusters; # entries", {HistType::kTH1F, {{10, 0.5, 10.5}}}},
                              {"mMFTTrackPhi", "Track #phi; #phi; # entries", {HistType::kTH1F, {{100, -3.2, 3.2}}}},
                              {"mMFTTrackTanl", "Track tan #lambda; tan #lambda; # entries", {HistType::kTH1F, {{100, -25, 0}}}},
                              {"mMFTTrackInvQPt", "Track q/p_{T}; q/p_{T} [1/GeV]; # entries", {HistType::kTH1F, {{250, -10, 10}}}}}};
  Configurable<bool> avClsPlots{"avClsPlots", false, "Enable average cluster plots"};

  void init(o2::framework::InitContext&)
  {
    if (avClsPlots) {
      registry.add("mMFTTrackAvgClusters", "Average number of clusters per track; p;# clusters; # entries", {HistType::kTH2F, {{100, 0, 100}, {100, 0, 100}}});
      registry.add("mMFTTrackAvgClustersTru", "Average number of clusters per track; p;# clusters; # entries", {HistType::kTH2F, {{100, 0, 100}, {100, 0, 100}}});
      if (doprocessMC) {
        registry.add("mMFTTrackAvgClustersHe", "Average number of clusters per track; p;# clusters; # entries", {HistType::kTH2F, {{100, 0, 100}, {100, 0, 100}}});
        registry.add("mMFTTrackAvgClustersTruHe", "Average number of clusters per track; p;# clusters; # entries", {HistType::kTH2F, {{100, 0, 100}, {100, 0, 100}}});
      }
    }
  }
  void process(aod::MFTTracks const& mfttracks)
  {
    for (auto& track : mfttracks) {
      // 2d histograms
      float x = track.x();
      float y = track.y();
      float eta = track.eta();
      float phi = track.phi();
      float nCls = track.nClusters();
      if (nCls >= 5) {
        registry.fill(HIST("mMFTTrackXY_5_MinClusters"), x, y);
        registry.fill(HIST("mMFTTrackEtaPhi_5_MinClusters"), eta, phi);
        if (nCls >= 7) {
          registry.fill(HIST("mMFTTrackXY_7_MinClusters"), x, y);
          registry.fill(HIST("mMFTTrackEtaPhi_7_MinClusters"), eta, phi);
          if (nCls >= 8) {
            registry.fill(HIST("mMFTTrackXY_8_MinClusters"), x, y);
            registry.fill(HIST("mMFTTrackEtaPhi_8_MinClusters"), eta, phi);
          }
        }
      }
      if (avClsPlots) {
        std::array<float, 10> clsSize;
        for (unsigned int layer = 0; layer < 10; layer++) {
          clsSize[layer] = (track.mftClusterSizesAndTrackFlags() >> (layer * 6)) & 0x3f;
          // LOG(info) << "Layer " << layer << ": " << clsSize[layer];
        }
        float avgCls = 0;
        for (unsigned int layer = 0; layer < 10; layer++) {
          avgCls += clsSize[layer];
        }
        avgCls /= track.nClusters();

        std::sort(clsSize.begin(), clsSize.end());
        float truncatedAvgCls = 0;
        int ncls = 0;
        for (unsigned int layer = 0; layer < 10; layer++) {
          if (clsSize[layer] > 0) {
            truncatedAvgCls += clsSize[layer];
            ncls++;
            if (ncls >= 3) {
              break; // we take the average of the first 5 non-zero clusters
            }
          }
        }
        truncatedAvgCls /= ncls;

        registry.fill(HIST("mMFTTrackAvgClusters"), track.p(), avgCls);
        registry.fill(HIST("mMFTTrackAvgClustersTru"), track.p(), truncatedAvgCls);
      }
      // 1d histograms
      registry.fill(HIST("mMFTTrackEta"), eta);
      registry.fill(HIST("mMFTTrackNumberOfClusters"), nCls);
      registry.fill(HIST("mMFTTrackPhi"), phi);
      registry.fill(HIST("mMFTTrackTanl"), track.tgl());
      registry.fill(HIST("mMFTTrackInvQPt"), track.signed1Pt());
    }
  }

  void processMC(soa::Join<aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks,
                 aod::McParticles const&)
  {
    for (auto& track : mfttracks) {
      if (avClsPlots) {
        std::array<float, 10> clsSize;
        for (unsigned int layer = 0; layer < 10; layer++) {
          clsSize[layer] = (track.mftClusterSizesAndTrackFlags() >> (layer * 6)) & 0x3f;
          // LOG(info) << "Layer " << layer << ": " << clsSize[layer];
        }
        float avgCls = 0;
        for (unsigned int layer = 0; layer < 10; layer++) {
          avgCls += clsSize[layer];
        }
        avgCls /= track.nClusters();

        std::sort(clsSize.begin(), clsSize.end());
        float truncatedAvgCls = 0;
        int ncls = 0;
        for (unsigned int layer = 0; layer < 10; layer++) {
          if (clsSize[layer] > 0) {
            truncatedAvgCls += clsSize[layer];
            ncls++;
            if (ncls >= 3) {
              break; // we take the average of the first 5 non-zero clusters
            }
          }
        }
        truncatedAvgCls /= ncls;

        if (track.has_mcParticle()) {
          const auto& mcParticle = track.mcParticle();
          if (std::abs(mcParticle.pdgCode()) == 1000020040) { // He4
            registry.fill(HIST("mMFTTrackAvgClustersHe"), track.p(), avgCls);
            registry.fill(HIST("mMFTTrackAvgClustersTruHe"), track.p(), truncatedAvgCls);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(CheckMFT, processMC, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CheckMFT>(cfgc),
  };
}
