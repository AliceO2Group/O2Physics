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
/// \file   qaTrackSelection.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \since  04/11/2022
/// \brief  Task to check how many tracks pass the cuts
///

// O2 includes
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/TrackSelectionTables.h"

using namespace o2::framework;

struct QaTrackCuts {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    const AxisSpec axisSelections{30, 0.5, 30.5f, "Selection"};
    // histos.add("events", "events", kTH1F, {axisSelections});
    auto h = histos.add<TH1>("tracks", "Tracks selection", kTH1F, {axisSelections});
    h->GetXaxis()->SetBinLabel(1, "Tracks read");
    h->GetXaxis()->SetBinLabel(2, "isGlobalTrackSDD");
    h->GetXaxis()->SetBinLabel(3, "passedTrackType");
    h->GetXaxis()->SetBinLabel(4, "passedPtRange");
    h->GetXaxis()->SetBinLabel(5, "passedEtaRange");
    h->GetXaxis()->SetBinLabel(6, "passedTPCNCls");
    h->GetXaxis()->SetBinLabel(7, "passedTPCCrossedRows");
    h->GetXaxis()->SetBinLabel(8, "passedTPCCrossedRowsOverNCls");
    h->GetXaxis()->SetBinLabel(9, "passedTPCChi2NDF");
    h->GetXaxis()->SetBinLabel(10, "passedTPCRefit");
    h->GetXaxis()->SetBinLabel(11, "passedITSNCls");
    h->GetXaxis()->SetBinLabel(12, "passedITSChi2NDF");
    h->GetXaxis()->SetBinLabel(13, "passedITSRefit");
    h->GetXaxis()->SetBinLabel(14, "passedITSHits");
    h->GetXaxis()->SetBinLabel(15, "passedGoldenChi2");
    h->GetXaxis()->SetBinLabel(16, "passedDCAxy");
    h->GetXaxis()->SetBinLabel(17, "passedDCAz");
    h->GetXaxis()->SetBinLabel(18, "isQualityTrack");
    h->GetXaxis()->SetBinLabel(19, "isPrimaryTrack");
    h->GetXaxis()->SetBinLabel(20, "isInAcceptanceTrack");
    h->GetXaxis()->SetBinLabel(21, "isGlobalTrack");
    h->GetXaxis()->SetBinLabel(22, "isGlobalTrackWoPtEta");
    h->GetXaxis()->SetBinLabel(23, "isGlobalTrackWoDCA");

    h = histos.add<TH1>("isQualityTrack", "Tracks selection for isQualityTrack", kTH1F, {axisSelections});
    h->GetXaxis()->SetBinLabel(1, "Tracks read");
    h->GetXaxis()->SetBinLabel(2, "passedTrackType");
    h->GetXaxis()->SetBinLabel(3, "passedTPCNCls");
    h->GetXaxis()->SetBinLabel(4, "passedTPCCrossedRows");
    h->GetXaxis()->SetBinLabel(5, "passedTPCCrossedRowsOverNCls");
    h->GetXaxis()->SetBinLabel(6, "passedTPCChi2NDF");
    h->GetXaxis()->SetBinLabel(7, "passedTPCRefit");
    h->GetXaxis()->SetBinLabel(8, "passedITSNCls");
    h->GetXaxis()->SetBinLabel(9, "passedITSChi2NDF");
    h->GetXaxis()->SetBinLabel(10, "passedITSRefit");
    h->GetXaxis()->SetBinLabel(11, "passedITSHits");
  }

  void process(const o2::soa::Join<o2::aod::Tracks, o2::aod::TrackSelection>& tracks,
               o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels> const&)
  {
    for (const auto& track : tracks) {
      histos.fill(HIST("tracks"), 1.f);
      if (track.isGlobalTrackSDD()) {
        histos.fill(HIST("tracks"), 2.f);
      }
      if (track.passedTrackType()) {
        histos.fill(HIST("tracks"), 3.f);
      }
      if (track.passedPtRange()) {
        histos.fill(HIST("tracks"), 4.f);
      }
      if (track.passedEtaRange()) {
        histos.fill(HIST("tracks"), 5.f);
      }
      if (track.passedTPCNCls()) {
        histos.fill(HIST("tracks"), 6.f);
      }
      if (track.passedTPCCrossedRows()) {
        histos.fill(HIST("tracks"), 7.f);
      }
      if (track.passedTPCCrossedRowsOverNCls()) {
        histos.fill(HIST("tracks"), 8.f);
      }
      if (track.passedTPCChi2NDF()) {
        histos.fill(HIST("tracks"), 9.f);
      }
      if (track.passedTPCRefit()) {
        histos.fill(HIST("tracks"), 10.f);
      }
      if (track.passedITSNCls()) {
        histos.fill(HIST("tracks"), 11.f);
      }
      if (track.passedITSChi2NDF()) {
        histos.fill(HIST("tracks"), 12.f);
      }
      if (track.passedITSRefit()) {
        histos.fill(HIST("tracks"), 13.f);
      }
      if (track.passedITSHits()) {
        histos.fill(HIST("tracks"), 14.f);
      }
      if (track.passedGoldenChi2()) {
        histos.fill(HIST("tracks"), 15.f);
      }
      if (track.passedDCAxy()) {
        histos.fill(HIST("tracks"), 16.f);
      }
      if (track.passedDCAz()) {
        histos.fill(HIST("tracks"), 17.f);
      }
      if (track.isQualityTrack()) {
        histos.fill(HIST("tracks"), 18.f);
      }
      if (track.isPrimaryTrack()) {
        histos.fill(HIST("tracks"), 19.f);
      }
      if (track.isInAcceptanceTrack()) {
        histos.fill(HIST("tracks"), 20.f);
      }
      if (track.isGlobalTrack()) {
        histos.fill(HIST("tracks"), 21.f);
      }
      if (track.isGlobalTrackWoPtEta()) {
        histos.fill(HIST("tracks"), 22.f);
      }
      if (track.isGlobalTrackWoDCA()) {
        histos.fill(HIST("tracks"), 23.f);
      }
      // isQualityTrack
      histos.fill(HIST("isQualityTrack"), 1.f);
      if (track.passedTrackType()) {
        histos.fill(HIST("isQualityTrack"), 2.f);
        if (track.passedTPCNCls()) {
          histos.fill(HIST("isQualityTrack"), 3.f);
          if (track.passedTPCCrossedRows()) {
            histos.fill(HIST("isQualityTrack"), 4.f);
            if (track.passedTPCCrossedRowsOverNCls()) {
              histos.fill(HIST("isQualityTrack"), 5.f);
              if (track.passedTPCChi2NDF()) {
                histos.fill(HIST("isQualityTrack"), 6.f);
                if (track.passedTPCRefit()) {
                  histos.fill(HIST("isQualityTrack"), 7.f);
                  if (track.passedITSNCls()) {
                    histos.fill(HIST("isQualityTrack"), 8.f);
                    if (track.passedITSChi2NDF()) {
                      histos.fill(HIST("isQualityTrack"), 9.f);
                      if (track.passedITSRefit()) {
                        histos.fill(HIST("isQualityTrack"), 10.f);
                        if (track.passedITSHits()) {
                          histos.fill(HIST("isQualityTrack"), 11.f);
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<QaTrackCuts>(cfgc)};
}
