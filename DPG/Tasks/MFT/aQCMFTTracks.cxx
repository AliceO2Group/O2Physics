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

#include "CCDB/BasicCCDBManager.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"

#include "Framework/DataTypes.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "CommonConstants/LHCConstants.h"
#include "Framework/TimingInfo.h"
#include "DataFormatsITSMFT/ROFRecord.h"

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
      // 1d histograms
      registry.fill(HIST("mMFTTrackEta"), eta);
      registry.fill(HIST("mMFTTrackNumberOfClusters"), nCls);
      registry.fill(HIST("mMFTTrackPhi"), phi);
      registry.fill(HIST("mMFTTrackTanl"), track.tgl());
      registry.fill(HIST("mMFTTrackInvQPt"), track.signed1Pt());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CheckMFT>(cfgc),
  };
}
