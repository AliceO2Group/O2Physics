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
/// \brief This task creates basic histograms for Pb-Pb multiplicity analysis.
///
/// \author hhesouno

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/ASoAHelpers.h"

#include "ReconstructionDataFormats/GlobalTrackID.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct multiplicityPbPb {

  // Configurable<float> estimatorEta{"estimatorEta", 1.0, "eta range for INEL>0 sample definition"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> nBinsEta{"nBinsEta", 100, "N bins in eta histo"};
  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histo"};
  Configurable<int> nBinsDCA{"nBinsDCA", 1000, "N bins in DCA histo"};
  Configurable<int> nBinsMult{"nBinsMult", 100, "N bins in Multiplicity histo"};
  Configurable<int> nBinsPhi{"nBinsPhi", 100, "N bins in phi histo"};
  Configurable<int> nBinsZvtx{"nBinsZvtx", 100, "N bins in Zvtx histo"};

  Filter trackDCA = nabs(aod::track::dcaXY) < 0.2f; // makes a big difference in etaHistogram

  using myCompleteTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA>;
  // using myCompleteTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,  aod::McTrackLabels>;
  using myFilteredTracks = soa::Filtered<myCompleteTracks>;

  // Preslice<aod::Tracks> perCollision = aod::track::collisionId;

  void init(InitContext const&)
  {
    // define axes you want to use
    const AxisSpec axisEta{nBinsEta, -2, 2, "#eta"};
    const AxisSpec axisPt{nBinsPt, 0, 5, "p_T"};
    const AxisSpec axisDeltaPt{100, -1.0, +1.0, "#Delta(p_{T})"};

    const AxisSpec axisCounter{1, 0, +1, ""};

    const AxisSpec axisDCAxy{nBinsDCA, -0.22, 0.22, "DCA_{xy} (cm)"};
    const AxisSpec axisDCAz{nBinsDCA, -0.22, 0.22, "DCA_{z} (cm)"};

    const AxisSpec axisNtrk{nBinsMult, 0, 300, "N_{trk}"};

    const AxisSpec axisPhi{nBinsPhi, -0.4, 6.8, "#phi"};
    const AxisSpec axisZvtx{nBinsZvtx, -30, 30, "Z_{vtx} (cm)"};

    histos.add("etaHistogram", "; ", kTH1F, {axisEta});
    histos.add("ptHistogram", "; ", kTH1F, {axisPt});
    //
    histos.add("eventCounter", "eventCounter", kTH1F, {axisCounter});

    histos.add("DCAxy", "; DCA_{xy} (cm)", kTH1F, {axisDCAxy});
    histos.add("DCAz", "; DCA_{z} (cm)", kTH1F, {axisDCAz});
    // do not know how:
    histos.add("Multiplicity", "; tracks; events", kTH1F, {axisNtrk});

    histos.add("Anton/PhiTracks", "; #phi; tracks", kTH1F, {axisPhi});
    histos.add("Anton/ZvtxEvents", "; Z_{vtx} (cm); events", kTH1F, {axisZvtx});

    histos.add("Anton/EtaZvtxTracks", "; #eta; Z_{vtx} (cm); tracks", kTH2F, {axisEta, axisZvtx});
    histos.add("Anton/NtrkZvtxEvents", "; N_{trk}; Z_{vtx} (cm); events", kTH2F, {axisNtrk, axisZvtx});
    histos.add("Anton/PhiEtaTracks", "; #phi; #eta; tracks", kTH2F, {axisPhi, axisEta});
  }

  // void process(aod::Collision const& collision, soa::Filtered<myCompleteTracks> const& tracks, aod::McParticles const&)
  void process(aod::Collision const& collision, soa::Filtered<myCompleteTracks> const& tracks)
  {
    int trackCounter = 0;

    // auto groupedTracks = tracks.sliceBy(perCollision, collision.globalIndex());

    histos.fill(HIST("eventCounter"), 0.5);

    // histos.fill(HIST("Multiplicity"), groupedTracks.size());

    histos.fill(HIST("Anton/ZvtxEvents"), collision.posZ());

    for (auto& track : tracks) {
      if (track.tpcNClsCrossedRows() < 70)
        continue; // badly tracked

      ++trackCounter;

      histos.fill(HIST("etaHistogram"), track.eta());
      histos.fill(HIST("ptHistogram"), track.pt());

      histos.fill(HIST("DCAxy"), track.dcaXY());
      histos.fill(HIST("DCAz"), track.dcaZ());

      histos.fill(HIST("Anton/PhiTracks"), track.phi());

      histos.fill(HIST("Anton/EtaZvtxTracks"), track.eta(), collision.posZ());
      histos.fill(HIST("Anton/PhiEtaTracks"), track.phi(), track.eta());
    }

    histos.fill(HIST("Multiplicity"), trackCounter);

    histos.fill(HIST("Anton/NtrkZvtxEvents"), trackCounter, collision.posZ());
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<multiplicityPbPb>(cfgc)};
}
