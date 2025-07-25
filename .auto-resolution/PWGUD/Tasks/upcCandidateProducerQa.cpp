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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/EventSelection.h"
#include "CommonConstants/LHCConstants.h"
#include "PWGUD/Core/UPCCutparHolder.h"
#include "PWGUD/Core/UPCHelpers.h"
#include "PWGUD/DataModel/UDTables.h"

using namespace o2::framework;
using namespace o2::framework::expressions;

struct UpcCandProducerQa {
  HistogramRegistry histRegistry{"HistRegistry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    const AxisSpec axisPt{500, 0., 5., ""};
    const AxisSpec axisEta{600, -6., 6., ""};
    const AxisSpec axisPhi{628, 0., 6.28, ""};
    const AxisSpec axisColNContrib{1001, -1., 1000., ""};

    histRegistry.add("TracksQA/Barrel/PtVsColNContrib/TPC", "", kTH2F, {axisPt, axisColNContrib});
    histRegistry.add("TracksQA/Barrel/PtVsColNContrib/TPC_TOF", "", kTH2F, {axisPt, axisColNContrib});
    histRegistry.add("TracksQA/Barrel/PtVsColNContrib/TPC_ITS", "", kTH2F, {axisPt, axisColNContrib});
    histRegistry.add("TracksQA/Barrel/PtVsColNContrib/TPC_TRD", "", kTH2F, {axisPt, axisColNContrib});
    histRegistry.add("TracksQA/Barrel/PtVsColNContrib/TPC_TRD_TOF", "", kTH2F, {axisPt, axisColNContrib});
    histRegistry.add("TracksQA/Barrel/PtVsColNContrib/TPC_ITS_TOF", "", kTH2F, {axisPt, axisColNContrib});
    histRegistry.add("TracksQA/Barrel/PtVsColNContrib/TPC_ITS_TRD", "", kTH2F, {axisPt, axisColNContrib});
    histRegistry.add("TracksQA/Barrel/PtVsColNContrib/TPC_ITS_TRD_TOF", "", kTH2F, {axisPt, axisColNContrib});

    histRegistry.add("TracksQA/Barrel/Pt/TPC", "", kTH1F, {axisPt});
    histRegistry.add("TracksQA/Barrel/Pt/TPC_TOF", "", kTH1F, {axisPt});
    histRegistry.add("TracksQA/Barrel/Pt/TPC_ITS", "", kTH1F, {axisPt});
    histRegistry.add("TracksQA/Barrel/Pt/TPC_TRD", "", kTH1F, {axisPt});
    histRegistry.add("TracksQA/Barrel/Pt/TPC_TRD_TOF", "", kTH1F, {axisPt});
    histRegistry.add("TracksQA/Barrel/Pt/TPC_ITS_TOF", "", kTH1F, {axisPt});
    histRegistry.add("TracksQA/Barrel/Pt/TPC_ITS_TRD", "", kTH1F, {axisPt});
    histRegistry.add("TracksQA/Barrel/Pt/TPC_ITS_TRD_TOF", "", kTH1F, {axisPt});

    histRegistry.add("TracksQA/Barrel/Pt_Gen/TPC", "", kTH1F, {axisPt});
    histRegistry.add("TracksQA/Barrel/Pt_Gen/TPC_TOF", "", kTH1F, {axisPt});
    histRegistry.add("TracksQA/Barrel/Pt_Gen/TPC_ITS", "", kTH1F, {axisPt});
    histRegistry.add("TracksQA/Barrel/Pt_Gen/TPC_TRD", "", kTH1F, {axisPt});
    histRegistry.add("TracksQA/Barrel/Pt_Gen/TPC_TRD_TOF", "", kTH1F, {axisPt});
    histRegistry.add("TracksQA/Barrel/Pt_Gen/TPC_ITS_TOF", "", kTH1F, {axisPt});
    histRegistry.add("TracksQA/Barrel/Pt_Gen/TPC_ITS_TRD", "", kTH1F, {axisPt});
    histRegistry.add("TracksQA/Barrel/Pt_Gen/TPC_ITS_TRD_TOF", "", kTH1F, {axisPt});

    histRegistry.add("TracksQA/Barrel/Eta/TPC", "", kTH1F, {axisEta});
    histRegistry.add("TracksQA/Barrel/Eta/TPC_TOF", "", kTH1F, {axisEta});
    histRegistry.add("TracksQA/Barrel/Eta/TPC_ITS", "", kTH1F, {axisEta});
    histRegistry.add("TracksQA/Barrel/Eta/TPC_TRD", "", kTH1F, {axisEta});
    histRegistry.add("TracksQA/Barrel/Eta/TPC_TRD_TOF", "", kTH1F, {axisEta});
    histRegistry.add("TracksQA/Barrel/Eta/TPC_ITS_TOF", "", kTH1F, {axisEta});
    histRegistry.add("TracksQA/Barrel/Eta/TPC_ITS_TRD", "", kTH1F, {axisEta});
    histRegistry.add("TracksQA/Barrel/Eta/TPC_ITS_TRD_TOF", "", kTH1F, {axisEta});

    histRegistry.add("TracksQA/Barrel/Phi/TPC", "", kTH1F, {axisPhi});
    histRegistry.add("TracksQA/Barrel/Phi/TPC_TOF", "", kTH1F, {axisPhi});
    histRegistry.add("TracksQA/Barrel/Phi/TPC_ITS", "", kTH1F, {axisPhi});
    histRegistry.add("TracksQA/Barrel/Phi/TPC_TRD", "", kTH1F, {axisPhi});
    histRegistry.add("TracksQA/Barrel/Phi/TPC_TRD_TOF", "", kTH1F, {axisPhi});
    histRegistry.add("TracksQA/Barrel/Phi/TPC_ITS_TOF", "", kTH1F, {axisPhi});
    histRegistry.add("TracksQA/Barrel/Phi/TPC_ITS_TRD", "", kTH1F, {axisPhi});
    histRegistry.add("TracksQA/Barrel/Phi/TPC_ITS_TRD_TOF", "", kTH1F, {axisPhi});

    histRegistry.add("TracksQA/Muon/Pt/MCH_MID", "", kTH1F, {axisPt});
    histRegistry.add("TracksQA/Muon/Pt_Gen/MCH_MID", "", kTH1F, {axisPt});

    histRegistry.add("TracksQA/MC/PtEl", "", kTH1F, {axisPt});
    histRegistry.add("TracksQA/MC/PtMu", "", kTH1F, {axisPt});
  }

  using BarrelTracks = o2::soa::Join<o2::aod::Tracks, o2::aod::TracksExtra>;
  using FwdTracks = o2::aod::FwdTracks;

  void updateBarrelTrackQA(const BarrelTracks::iterator& track, int32_t colNContrib)
  {
    // check basic cuts
    // only tracks in TOF acceptance
    if (track.pt() < 0.3 || std::abs(track.eta()) > 0.8 || track.tpcNClsCrossedRows() < 70)
      return;

    int8_t mask = 0;
    if (track.hasTPC())
      SETBIT(mask, 0);
    if (track.hasITS())
      SETBIT(mask, 1);
    if (track.hasTOF())
      SETBIT(mask, 2);
    if (track.hasTRD())
      SETBIT(mask, 3);

    float pt = track.pt();
    float eta = track.eta();
    float phi = track.phi();

    if (mask == 1) {
      histRegistry.fill(HIST("TracksQA/Barrel/PtVsColNContrib/TPC"), pt, colNContrib);
      histRegistry.fill(HIST("TracksQA/Barrel/Pt/TPC"), pt);
      histRegistry.fill(HIST("TracksQA/Barrel/Eta/TPC"), eta);
      histRegistry.fill(HIST("TracksQA/Barrel/Phi/TPC"), phi);
    }
    if (mask == 3) {
      histRegistry.fill(HIST("TracksQA/Barrel/PtVsColNContrib/TPC_ITS"), pt, colNContrib);
      histRegistry.fill(HIST("TracksQA/Barrel/Pt/TPC_ITS"), pt);
      histRegistry.fill(HIST("TracksQA/Barrel/Eta/TPC_ITS"), eta);
      histRegistry.fill(HIST("TracksQA/Barrel/Phi/TPC_ITS"), phi);
    }
    if (mask == 5) {
      histRegistry.fill(HIST("TracksQA/Barrel/PtVsColNContrib/TPC_TOF"), pt, colNContrib);
      histRegistry.fill(HIST("TracksQA/Barrel/Pt/TPC_TOF"), pt);
      histRegistry.fill(HIST("TracksQA/Barrel/Eta/TPC_TOF"), eta);
      histRegistry.fill(HIST("TracksQA/Barrel/Phi/TPC_TOF"), phi);
    }
    if (mask == 7) {
      histRegistry.fill(HIST("TracksQA/Barrel/PtVsColNContrib/TPC_ITS_TOF"), pt, colNContrib);
      histRegistry.fill(HIST("TracksQA/Barrel/Pt/TPC_ITS_TOF"), pt);
      histRegistry.fill(HIST("TracksQA/Barrel/Eta/TPC_ITS_TOF"), eta);
      histRegistry.fill(HIST("TracksQA/Barrel/Phi/TPC_ITS_TOF"), phi);
    }
    if (mask == 9) {
      histRegistry.fill(HIST("TracksQA/Barrel/PtVsColNContrib/TPC_TRD"), pt, colNContrib);
      histRegistry.fill(HIST("TracksQA/Barrel/Pt/TPC_TRD"), pt);
      histRegistry.fill(HIST("TracksQA/Barrel/Eta/TPC_TRD"), eta);
      histRegistry.fill(HIST("TracksQA/Barrel/Phi/TPC_TRD"), phi);
    }
    if (mask == 11) {
      histRegistry.fill(HIST("TracksQA/Barrel/PtVsColNContrib/TPC_ITS_TRD"), pt, colNContrib);
      histRegistry.fill(HIST("TracksQA/Barrel/Pt/TPC_ITS_TRD"), pt);
      histRegistry.fill(HIST("TracksQA/Barrel/Eta/TPC_ITS_TRD"), eta);
      histRegistry.fill(HIST("TracksQA/Barrel/Phi/TPC_ITS_TRD"), phi);
    }
    if (mask == 13) {
      histRegistry.fill(HIST("TracksQA/Barrel/PtVsColNContrib/TPC_TRD_TOF"), pt, colNContrib);
      histRegistry.fill(HIST("TracksQA/Barrel/Pt/TPC_TRD_TOF"), pt);
      histRegistry.fill(HIST("TracksQA/Barrel/Eta/TPC_TRD_TOF"), eta);
      histRegistry.fill(HIST("TracksQA/Barrel/Phi/TPC_TRD_TOF"), phi);
    }
    if (mask == 15) {
      histRegistry.fill(HIST("TracksQA/Barrel/PtVsColNContrib/TPC_ITS_TRD_TOF"), pt, colNContrib);
      histRegistry.fill(HIST("TracksQA/Barrel/Pt/TPC_ITS_TRD_TOF"), pt);
      histRegistry.fill(HIST("TracksQA/Barrel/Eta/TPC_ITS_TRD_TOF"), eta);
      histRegistry.fill(HIST("TracksQA/Barrel/Phi/TPC_ITS_TRD_TOF"), phi);
    }
  }

  void updateFwdTrackQA(const FwdTracks::iterator& fwdtrack)
  {
    if (fwdtrack.eta() < -4. || fwdtrack.eta() > -2.5)
      return;
    histRegistry.fill(HIST("TracksQA/Muon/Pt/MCH_MID"), fwdtrack.pt());
  }

  void process(o2::aod::Collisions const&,
               BarrelTracks const& tracks, o2::aod::AmbiguousTracks const& ambTracks,
               FwdTracks const& fwdtracks)
  {

    std::unordered_set<int64_t> ambTrIds;
    for (const auto& ambTrk : ambTracks) {
      auto trkId = ambTrk.trackId();
      ambTrIds.insert(trkId);
    }

    for (const auto& track : tracks) {
      int32_t nContrib = -1;
      if (ambTrIds.find(track.globalIndex()) == ambTrIds.end()) {
        const auto& col = track.collision();
        nContrib = col.numContrib();
      }
      updateBarrelTrackQA(track, nContrib);
    }

    ambTrIds.clear();

    for (const auto& fwdtrack : fwdtracks) {
      if (fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)
        continue;
      updateFwdTrackQA(fwdtrack);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UpcCandProducerQa>(cfgc)};
}
