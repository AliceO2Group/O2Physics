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

/// \file AntiNucleiTask.cxx
/// \brief A task to analyse Anti-nuclei
/// \author Arkaprabha Saha <arkaprabha.saha@cern.ch>

#include "Common/Core/PID/TPCPIDResponse.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "DataFormatsTPC/BetheBlochAleph.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include <TParameter.h>

#include <cmath>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using CollisionWithEvSel = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>;
using TotalTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTOFDe>;

namespace
{
static const std::vector<std::string> particleName{"d"};
static const double kBetheBlochDefault[6]{-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32};
static const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};
static const float maxEtaCut = 0.8f;
static const int minTpcCrossedRowsCut = 70;
static const float maxVertexZCut = 10.f;
} // namespace

struct AntiNucleiTask {
  // Histogram registry: for holding histograms
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Configurable track cuts
  Configurable<float> trackNclusTPCcut{"trackNclusTPCcut", 70.0f, "min number of TPC clusters"};
  Configurable<float> trackNclusITScut{"trackNclusITScut", 4.0f, "min number of ITS clusters"};
  Configurable<float> maxChi2TPC{"maxChi2TPC", 4.0f, "max chi2 per cluster TPC"};
  Configurable<float> minChi2TPC{"minChi2TPC", 0.0f, "min chi2 per cluster TPC"};
  Configurable<float> chi2ITS{"chi2ITS", 36.0f, "max chi2 per cluster ITS"};
  Configurable<float> trackDCAz{"trackDCAz", 0.1f, "maxDCAz"};
  Configurable<float> trackDCAxy{"trackDCAxy", 0.1f, "maxDCAxy"};
  Configurable<float> tpcNSigmaCut{"tpcNSigmaCut", 3.0f, "tpcNSigmaCut"};
  Configurable<LabeledArray<double>> cfgBetheBlochParams{"cfgBetheBlochParams", {kBetheBlochDefault, 1, 6, particleName, betheBlochParNames}, "TPC Bethe-Bloch parameterisation for deuteron"};

  void init(InitContext const&)
  { // Defining the Histogram Axes
    ConfigurableAxis etaAxis{"etaAxis", {16, -0.8, +0.8}, "#eta"};
    ConfigurableAxis phiAxis{"phiAxis", {70, 0.f, 7.f}, "#phi"};
    ConfigurableAxis zVtxAxis{"zVtxAxis", {100, -20.f, 20.f}, "Primary Vertex z (cm)"};
    ConfigurableAxis nSigmaAxis{"nSigmaAxis", {50, -5.f, 5.f}, "N_{#sigma}"};
    ConfigurableAxis ptAxis{"ptAxis", {200, -10.0f, 10.0f}, "p_{T} (GeV/c)"};
    ConfigurableAxis centAxis{"centAxis", {100, 0, 100.0f}, "Centrality"};
    ConfigurableAxis momAxis{"momAxis", {5.e2, 0.f, 5.f}, "momentum axis binning"};
    ConfigurableAxis tpcAxis{"tpcAxis", {4.e2, 0.f, 4.e3f}, "tpc signal axis binning"};

    // Creating histograms
    histos.add("RawzVtx", "RawzVtx", kTH1F, {{zVtxAxis, "Primary Vertex z (cm)"}});
    histos.add("zVtx", "zVtx", kTH1F, {{zVtxAxis, "Primary Vertex z (cm)"}});
    histos.add("RawEta", "RawEta", kTH1F, {{etaAxis, "#eta"}});
    histos.add("Eta", "Eta", kTH1F, {{etaAxis, "#eta"}});
    histos.add("RawPhi", "RawPhi", kTH1F, {{phiAxis, "#phi (rad)"}});
    histos.add("Phi", "Phi", kTH1F, {{phiAxis, "#phi (rad)"}});
    histos.add("RawPt", "RawPt", kTH1F, {{ptAxis, "#it{p}_{T} (GeV/#it{c})"}});
    histos.add("Pt", "Pt", kTH1F, {{ptAxis, "#it{p}_{T} (GeV/#it{c})"}});
    histos.add("TpcSignal", "TpcSignal", kTH2F, {{momAxis, "#it{p}_{TPC} (GeV/#it{c})"}, {tpcAxis, "d#it{E}/d#it{x}_{TPC}"}});
    histos.add("RawtpcNSigma", "RawtpcNSigma", kTH3F, {{centAxis, "Centrality"}, {ptAxis, "#it{p}_{T} (GeV/#it{c})"}, {nSigmaAxis, "N_{#sigma}"}});
    histos.add("tpcNSigma", "tpcNSigma", kTH3F, {{centAxis, "Centrality"}, {ptAxis, "#it{p}_{T} (GeV/#it{c})"}, {nSigmaAxis, "N_{#sigma}"}});
    histos.add("RawtofNSigma", "RawtofNSigma", kTH3F, {{centAxis, "Centrality"}, {ptAxis, "#it{p}_{T} (GeV/#it{c})"}, {nSigmaAxis, "N_{#sigma}"}});
    histos.add("tofNSigma", "tofNSigma", kTH3F, {{centAxis, "Centrality"}, {ptAxis, "#it{p}_{T} (GeV/#it{c})"}, {nSigmaAxis, "N_{#sigma}"}});
  }

  // Function to apply track cuts
  template <typename T>
  bool isGoodTrack(const T& track)
  {
    if (track.eta() > maxEtaCut)
      return false;
    if (track.tpcNClsFound() < trackNclusTPCcut)
      return false;
    if (track.tpcNClsCrossedRows() < minTpcCrossedRowsCut)
      return false;
    if (track.itsNCls() < trackNclusITScut)
      return false;
    if (track.tpcChi2NCl() > maxChi2TPC)
      return false;
    if (track.tpcChi2NCl() < minChi2TPC)
      return false;
    if (track.itsChi2NCl() > chi2ITS)
      return false;
    if (std::abs(track.dcaXY()) > trackDCAxy)
      return false;
    if (std::abs(track.dcaZ()) > trackDCAz)
      return false;

    return true;
  }

  // The process function
  void process(CollisionWithEvSel::iterator const& collision, TotalTracks const& tracks)
  {
    // Event Selection
    if (std::abs(collision.posZ()) > maxVertexZCut) {
      return;
    }

    // Filling the z-vertex histogram before the event selection cuts.
    histos.fill(HIST("RawzVtx"), collision.posZ());

    // Applying the built-in O2 event selection (sel8).
    if (!collision.sel8()) {
      return;
    }

    // Filling the z-vertex histogram after the event selection cuts.
    histos.fill(HIST("zVtx"), collision.posZ());

    // Track Selection
    for (const auto& track : tracks) {

      double expBethe{tpc::BetheBlochAleph(static_cast<double>(track.tpcInnerParam()), cfgBetheBlochParams->get("p0"), cfgBetheBlochParams->get("p1"), cfgBetheBlochParams->get("p2"), cfgBetheBlochParams->get("p3"), cfgBetheBlochParams->get("p4"))};
      double expSigma{expBethe * cfgBetheBlochParams->get("resolution")};
      float tpcNSigma = static_cast<float>((track.tpcSignal() - expBethe) / expSigma);

      float pt = track.sign() > 0 ? 2 * track.pt() : -2 * track.pt();
      // Filling histograms with track data before applying any cuts.
      histos.fill(HIST("RawEta"), track.eta());
      histos.fill(HIST("RawPhi"), track.phi());
      histos.fill(HIST("RawPt"), pt);
      histos.fill(HIST("RawtpcNSigma"), collision.centFT0C(), pt, tpcNSigma);
      histos.fill(HIST("RawtofNSigma"), collision.centFT0C(), pt, track.tofNSigmaDe());

      // If the track is good, fill the "after cuts" histograms.
      if (isGoodTrack(track)) {
        histos.fill(HIST("Eta"), track.eta());
        histos.fill(HIST("Phi"), track.phi());
        histos.fill(HIST("Pt"), pt);
        histos.fill(HIST("tpcNSigma"), collision.centFT0C(), pt, tpcNSigma);
        histos.fill(HIST("TpcSignal"), track.tpcInnerParam(), track.tpcSignal());

        if (std::abs(tpcNSigma) < tpcNSigmaCut) {
          histos.fill(HIST("tofNSigma"), collision.centFT0C(), pt, track.tofNSigmaDe());
        }
      }
    }
  }

  PROCESS_SWITCH(AntiNucleiTask, process, "process", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<AntiNucleiTask>(cfgc)};
}
