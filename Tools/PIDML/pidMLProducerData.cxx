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

/// \file pidMLProducerData.cxx
/// \brief Produce PID ML skimmed data from data files.
///
/// \author Maja Kabus <mkabus@cern.ch>
///
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Tools/PIDML/pidML.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"

struct PidMlProducerData {
  Produces<aod::PidTracksDataMl> pidTracksTableML;
  Produces<aod::PidTracksData> pidTracksTable;

  Filter trackFilter = requireGlobalTrackInFilter();
  using BigTracksML = soa::Filtered<soa::Join<aod::FullTracks, aod::TracksDCA, aod::pidTOFbeta, aod::TrackSelection, aod::TOFSignal>>;
  using MyCollisionML = aod::Collisions::iterator;
  using BigTracks = soa::Filtered<soa::Join<aod::FullTracks, aod::TracksDCA, aod::pidTOFbeta, aod::pidTPCFullEl, aod::pidTOFFullEl, aod::pidTPCFullMu, aod::pidTOFFullMu, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::TrackSelection, aod::TOFSignal>>;
  using MyCollision = soa::Join<aod::Collisions, aod::CentRun2V0Ms, aod::Mults>::iterator;

  static constexpr float kEps = 1e-10f;
  static constexpr float kMissingBeta = -999.0f;
  static constexpr float kMissingTOFSignal = -999.0f;

  HistogramRegistry registry{
    "registry",
    {
      {"minus/hTPCSigvsP", "TPC signal vs #it{p};#it{p} (GeV/#it{c});TPC signal", {HistType::kTH2F, {{500, 0., 10.}, {1000, 0., 600.}}}},
      {"minus/hTOFBetavsP", "TOF beta vs #it{p};#it{p} (GeV/#it{c});TOF #it{beta}", {HistType::kTH2F, {{500, 0., 10.}, {6000, -3., 3.}}}},
      {"minus/hTOFSigvsP", "TOF signal vs #it{p};#it{p} (GeV/#it{c});TOF signal", {HistType::kTH2F, {{500, 0., 10.}, {10000, -5000., 80000.}}}},
      {"minus/filtered/hTOFSigvsP", "TOF signal (filtered) vs #it{p};#it{p} (GeV/#it{c});TOF signal", {HistType::kTH2F, {{500, 0., 10.}, {10000, -5000., 80000.}}}},
      {"minus/hTRDPattvsP", "TRD pattern vs #it{p};#it{p} (GeV/#it{c});TRD pattern", {HistType::kTH2F, {{500, 0., 10.}, {110, -10., 100.}}}},
      {"minus/hTRDSigvsP", "TRD signal vs #it{p};#it{p} (GeV/#it{c});TRD signal", {HistType::kTH2F, {{500, 0., 10.}, {2500, -2., 100.}}}},
      {"minus/hP", "#it{p};#it{p} (GeV/#it{c})", {HistType::kTH1F, {{500, 0., 6.}}}},
      {"minus/hPt", "#it{p}_{t};#it{p}_{t} (GeV/#it{c})", {HistType::kTH1F, {{500, 0., 6.}}}},
      {"minus/hPx", "#it{p}_{x};#it{p}_{x} (GeV/#it{c})", {HistType::kTH1F, {{1000, -6., 6.}}}},
      {"minus/hPy", "#it{p}_{y};#it{p}_{y} (GeV/#it{c})", {HistType::kTH1F, {{1000, -6., 6.}}}},
      {"minus/hPz", "#it{p}_{z};#it{p}_{z} (GeV/#it{c})", {HistType::kTH1F, {{1000, -6., 6.}}}},
      {"minus/hX", "#it{x};#it{x}", {HistType::kTH1F, {{1000, -2., 2.}}}},
      {"minus/hY", "#it{y};#it{y}", {HistType::kTH1F, {{1000, -2., 2.}}}},
      {"minus/hZ", "#it{z};#it{z}", {HistType::kTH1F, {{1000, -10., 10.}}}},
      {"minus/hAlpha", "#{alpha};#{alpha}", {HistType::kTH1F, {{1000, -5., 5.}}}},
      {"minus/hTrackType", "Track Type;Track Type", {HistType::kTH1F, {{300, 0., 300.}}}},
      {"minus/hTPCNClsShared", "hTPCNClsShared;hTPCNClsShared", {HistType::kTH1F, {{100, 0., 100.}}}},
      {"minus/hDcaXY", "#it{DcaXY};#it{DcaXY}", {HistType::kTH1F, {{1000, -1., 1.}}}},
      {"minus/hDcaZ", "#it{DcaZ};#it{DcaZ}", {HistType::kTH1F, {{1000, -1., 1.}}}},
      {"plus/hTPCSigvsP", "TPC signal vs #it{p};#it{p} (GeV/#it{c});TPC signal", {HistType::kTH2F, {{500, 0., 10.}, {1000, 0., 600.}}}},
      {"plus/hTOFBetavsP", "TOF beta vs #it{p};#it{p} (GeV/#it{c});TOF #it{beta}", {HistType::kTH2F, {{500, 0., 10.}, {6000, -3., 3.}}}},
      {"plus/hTOFSigvsP", "TOF signal vs #it{p};#it{p} (GeV/#it{c});TOF signal", {HistType::kTH2F, {{500, 0., 10.}, {10000, -5000., 80000.}}}},
      {"plus/filtered/hTOFSigvsP", "TOF signal (filtered) vs #it{p};#it{p} (GeV/#it{c});TOF signal", {HistType::kTH2F, {{500, 0., 10.}, {10000, -5000., 80000.}}}},
      {"plus/hTRDPattvsP", "TRD pattern vs #it{p};#it{p} (GeV/#it{c});TRD pattern", {HistType::kTH2F, {{500, 0., 10.}, {110, -10., 100.}}}},
      {"plus/hTRDSigvsP", "TRD signal vs #it{p};#it{p} (GeV/#it{c});TRD signal", {HistType::kTH2F, {{500, 0., 10.}, {2500, -2., 100.}}}},
      {"plus/hP", "#it{p};#it{p} (GeV/#it{c})", {HistType::kTH1F, {{500, 0., 6.}}}},
      {"plus/hPt", "#it{p}_{t};#it{p}_{t} (GeV/#it{c})", {HistType::kTH1F, {{500, 0., 6.}}}},
      {"plus/hPx", "#it{p}_{x};#it{p}_{x} (GeV/#it{c})", {HistType::kTH1F, {{1000, -6., 6.}}}},
      {"plus/hPy", "#it{p}_{y};#it{p}_{y} (GeV/#it{c})", {HistType::kTH1F, {{1000, -6., 6.}}}},
      {"plus/hPz", "#it{p}_{z};#it{p}_{z} (GeV/#it{c})", {HistType::kTH1F, {{1000, -6., 6.}}}},
      {"plus/hX", "#it{x};#it{x}", {HistType::kTH1F, {{1000, -2., 2.}}}},
      {"plus/hY", "#it{y};#it{y}", {HistType::kTH1F, {{1000, -2., 2.}}}},
      {"plus/hZ", "#it{z};#it{z}", {HistType::kTH1F, {{1000, -10., 10.}}}},
      {"plus/hAlpha", "#{alpha};#{alpha}", {HistType::kTH1F, {{1000, -5., 5.}}}},
      {"plus/hTrackType", "Track Type;Track Type", {HistType::kTH1F, {{300, 0., 300.}}}},
      {"plus/hTPCNClsShared", "hTPCNClsShared;hTPCNClsShared", {HistType::kTH1F, {{100, 0., 100.}}}},
      {"plus/hDcaXY", "#it{DcaXY};#it{DcaXY}", {HistType::kTH1F, {{1000, -1., 1.}}}},
      {"plus/hDcaZ", "#it{DcaZ};#it{DcaZ}", {HistType::kTH1F, {{1000, -1., 1.}}}},
    }};

  template <typename T>
  void fillHistos(const T& track)
  {
    if (track.sign() < 0) {
      registry.fill(HIST("minus/hTPCSigvsP"), track.p(), track.tpcSignal());
      registry.fill(HIST("minus/hTOFBetavsP"), track.p(), track.beta());
      registry.fill(HIST("minus/hTOFSigvsP"), track.p(), track.tofSignal());
      if (TMath::Abs(track.beta() - kMissingBeta) >= kEps) {
        registry.fill(HIST("minus/filtered/hTOFSigvsP"), track.p(), track.tofSignal());
      } else {
        registry.fill(HIST("minus/filtered/hTOFSigvsP"), track.p(), kMissingTOFSignal);
      }
      registry.fill(HIST("minus/hTRDPattvsP"), track.p(), track.trdPattern());
      registry.fill(HIST("minus/hTRDSigvsP"), track.p(), track.trdSignal());
      registry.fill(HIST("minus/hP"), track.p());
      registry.fill(HIST("minus/hPt"), track.pt());
      registry.fill(HIST("minus/hPx"), track.px());
      registry.fill(HIST("minus/hPy"), track.py());
      registry.fill(HIST("minus/hPz"), track.pz());
      registry.fill(HIST("minus/hX"), track.x());
      registry.fill(HIST("minus/hY"), track.y());
      registry.fill(HIST("minus/hZ"), track.z());
      registry.fill(HIST("minus/hAlpha"), track.alpha());
      registry.fill(HIST("minus/hTrackType"), track.trackType());
      registry.fill(HIST("minus/hTPCNClsShared"), track.tpcNClsShared());
      registry.fill(HIST("minus/hDcaXY"), track.dcaXY());
      registry.fill(HIST("minus/hDcaZ"), track.dcaZ());
    } else {
      registry.fill(HIST("plus/hTPCSigvsP"), track.p(), track.tpcSignal());
      registry.fill(HIST("plus/hTOFBetavsP"), track.p(), track.beta());
      registry.fill(HIST("plus/hTOFSigvsP"), track.p(), track.tofSignal());
      if (TMath::Abs(track.beta() - kMissingBeta) >= kEps) {
        registry.fill(HIST("plus/filtered/hTOFSigvsP"), track.p(), track.tofSignal());
      } else {
        registry.fill(HIST("plus/filtered/hTOFSigvsP"), track.p(), kMissingTOFSignal);
      }
      registry.fill(HIST("plus/hTRDPattvsP"), track.p(), track.trdPattern());
      registry.fill(HIST("plus/hTRDSigvsP"), track.p(), track.trdSignal());
      registry.fill(HIST("plus/hP"), track.p());
      registry.fill(HIST("plus/hPt"), track.pt());
      registry.fill(HIST("plus/hPx"), track.px());
      registry.fill(HIST("plus/hPy"), track.py());
      registry.fill(HIST("plus/hPz"), track.pz());
      registry.fill(HIST("plus/hX"), track.x());
      registry.fill(HIST("plus/hY"), track.y());
      registry.fill(HIST("plus/hZ"), track.z());
      registry.fill(HIST("plus/hAlpha"), track.alpha());
      registry.fill(HIST("plus/hTrackType"), track.trackType());
      registry.fill(HIST("plus/hTPCNClsShared"), track.tpcNClsShared());
      registry.fill(HIST("plus/hDcaXY"), track.dcaXY());
      registry.fill(HIST("plus/hDcaZ"), track.dcaZ());
    }
  }

  void processML(MyCollisionML const& /*collision*/, BigTracksML const& tracks)
  {
    for (const auto& track : tracks) {
      pidTracksTableML(track.tpcSignal(), track.trdSignal(), track.trdPattern(),
                       track.tofSignal(), track.beta(),
                       track.p(), track.pt(), track.px(), track.py(), track.pz(),
                       track.sign(),
                       track.x(), track.y(), track.z(),
                       track.alpha(),
                       track.trackType(),
                       track.tpcNClsShared(),
                       track.dcaXY(), track.dcaZ());

      fillHistos(track);
    }
  }
  PROCESS_SWITCH(PidMlProducerData, processML, "Produce only ML real data", true);

  void processAll(MyCollision const& collision, BigTracks const& tracks)
  {
    for (const auto& track : tracks) {
      pidTracksTable(collision.centRun2V0M(),
                     collision.multFV0A(), collision.multFV0C(), collision.multFV0M(),
                     collision.multFT0A(), collision.multFT0C(), collision.multFT0M(),
                     collision.multZNA(), collision.multZNC(),
                     collision.multTracklets(), collision.multTPC(),
                     track.tpcSignal(), track.trdSignal(), track.trdPattern(),
                     track.trackEtaEmcal(), track.trackPhiEmcal(),
                     track.tofSignal(), track.beta(),
                     track.p(), track.pt(), track.px(), track.py(), track.pz(),
                     track.sign(),
                     track.x(), track.y(), track.z(),
                     track.alpha(),
                     track.trackType(),
                     track.tpcNClsShared(),
                     track.dcaXY(), track.dcaZ(),
                     track.tpcNSigmaEl(), track.tpcExpSigmaEl(), track.tpcExpSignalDiffEl(),
                     track.tofNSigmaEl(), track.tofExpSigmaEl(), track.tofExpSignalDiffEl(),
                     track.tpcNSigmaMu(), track.tpcExpSigmaMu(), track.tpcExpSignalDiffMu(),
                     track.tofNSigmaMu(), track.tofExpSigmaMu(), track.tofExpSignalDiffMu(),
                     track.tpcNSigmaPi(), track.tpcExpSigmaPi(), track.tpcExpSignalDiffPi(),
                     track.tofNSigmaPi(), track.tofExpSigmaPi(), track.tofExpSignalDiffPi(),
                     track.tpcNSigmaKa(), track.tpcExpSigmaKa(), track.tpcExpSignalDiffKa(),
                     track.tofNSigmaKa(), track.tofExpSigmaKa(), track.tofExpSignalDiffKa(),
                     track.tpcNSigmaPr(), track.tpcExpSigmaPr(), track.tpcExpSignalDiffPr(),
                     track.tofNSigmaPr(), track.tofExpSigmaPr(), track.tofExpSignalDiffPr());

      fillHistos(track);
    }
  }
  PROCESS_SWITCH(PidMlProducerData, processAll, "Produce all real data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PidMlProducerData>(cfgc)};
}
