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
/// \brief This task creates basic histograms for Pb-Pb multiplicity analysis. It produces fake real and testing simulated data.
///
/// \author hhesouno

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"

#include "ReconstructionDataFormats/GlobalTrackID.h"

#include "Common/DataModel/EventSelection.h"

#include "Framework/O2DatabasePDGPlugin.h"

#include <TRandom.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct mcdata {

  Service<o2::framework::O2DatabasePDG> pdg;

  // Configurable<float> estimatorEta{"estimatorEta", 1.0, "eta range for INEL>0 sample definition"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> nBinsEta{"nBinsEta", 100, "N bins in eta histo"};
  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histo"};
  Configurable<int> nBinsDCA{"nBinsDCA", 1000, "N bins in DCA histo"};
  Configurable<int> nBinsMult{"nBinsMult", 100, "N bins in Multiplicity histo"};
  Configurable<int> nBinsPhi{"nBinsPhi", 100, "N bins in phi histo"};
  Configurable<int> nBinsZvtx{"nBinsZvtx", 100, "N bins in Zvtx histo"};

  Filter trackDCA = nabs(aod::track::dcaXY) < 0.2f; // makes a big difference in etaHistogram

  using myCompleteTracks = soa::Join<aod::Tracks, aod::TracksDCA>;
  using myFilteredTracks = soa::Filtered<myCompleteTracks>;

  void init(InitContext const&)
  {
    // define axes you want to use
    const AxisSpec axisEta{nBinsEta, -2, 2, "#eta"};
    const AxisSpec axisCounter{1, 0, +1, ""};
    const AxisSpec axisNtrk{nBinsMult, 0, 300, "N_{trk}"};
    const AxisSpec axisPhi{nBinsPhi, -0.4, 6.8, "#phi"};
    const AxisSpec axisZvtx{nBinsZvtx, -30, 30, "Z_{vtx} (cm)"};

    histos.add("etaHistogram", "; ", kTH1D, {axisEta});
    histos.add("MCGENetaHistogram", "; ", kTH1D, {axisEta});
    histos.add("FakeDataetaHistogram", "; ", kTH1D, {axisEta});
    histos.add("MCGENFakeDataetaHistogram", "; ", kTH1D, {axisEta});
    histos.add("eventCounter", "eventCounter", kTH1D, {axisCounter});
    histos.add("FakeDataeventCounter", "eventCounter", kTH1D, {axisCounter});
    histos.add("MCGENeventCounter", "eventCounter", kTH1D, {axisCounter});
    histos.add("MCGENFakeDataeventCounter", "eventCounter", kTH1D, {axisCounter});
    histos.add("Multiplicity", "; tracks; events", kTH1D, {axisNtrk});
    histos.add("FakeDataMultiplicity", "; tracks; events", kTH1D, {axisNtrk});
    histos.add("MCGENMultiplicity", "; tracks; events", kTH1D, {axisNtrk});
    histos.add("PhiTracks", "; #phi; tracks", kTH1D, {axisPhi});
    histos.add("FakeDataPhiTracks", "; #phi; tracks", kTH1D, {axisPhi});

    histos.add("ZvtxEvents", "; Z_{vtx} (cm); events", kTH1D, {axisZvtx});
    histos.add("FakeDataZvtxEvents", "; Z_{vtx} (cm); events", kTH1D, {axisZvtx});
    histos.add("MCGENZvtxEvents", "; Z_{vtx} (cm); events", kTH1D, {axisZvtx});

    histos.add("EtaZvtxTracks", "; #eta; Z_{vtx} (cm); tracks", kTH2D, {axisEta, axisZvtx});
    histos.add("FakeDataEtaZvtxTracks", "; #eta; Z_{vtx} (cm); tracks", kTH2D, {axisEta, axisZvtx});
    histos.add("MCGENEtaZvtxTracks", "; #eta; Z_{vtx} (cm); tracks", kTH2D, {axisEta, axisZvtx});
  }

  // void process(aod::Collision const& collision, soa::Filtered<myCompleteTracks> const& tracks, aod::McParticles const&)
  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<myCompleteTracks> const& tracks)
  {
    if (!collision.sel8()) {
      return;
    }

    int random = gRandom->Integer(2);
    if (random == 0) {
      int trackCounter = 0;

      histos.fill(HIST("eventCounter"), 0.5);

      histos.fill(HIST("ZvtxEvents"), collision.posZ());

      for (auto& track : tracks) {
        ++trackCounter;

        histos.fill(HIST("etaHistogram"), track.eta());
        histos.fill(HIST("PhiTracks"), track.phi());
        histos.fill(HIST("EtaZvtxTracks"), track.eta(), collision.posZ());
      }

      histos.fill(HIST("Multiplicity"), trackCounter);
    }

    if (random == 1) {
      int faketrackCounter = 0;

      histos.fill(HIST("FakeDataeventCounter"), 0.5);

      histos.fill(HIST("FakeDataZvtxEvents"), collision.posZ());

      for (auto& track : tracks) {
        ++faketrackCounter;

        histos.fill(HIST("FakeDataetaHistogram"), track.eta());
        histos.fill(HIST("FakeDataPhiTracks"), track.phi());
        histos.fill(HIST("FakeDataEtaZvtxTracks"), track.eta(), collision.posZ());
      }

      histos.fill(HIST("FakeDataMultiplicity"), faketrackCounter);
    }
  }

  void processMCGEN(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles)
  {
    int random = gRandom->Integer(2);
    if (random == 0) {
      int MCparticleCounter = 0;

      histos.fill(HIST("MCGENeventCounter"), 0.5);
      histos.fill(HIST("MCGENZvtxEvents"), mcCollision.posZ());

      for (auto& mcParticle : mcParticles) {
        ++MCparticleCounter;
        if (mcParticle.isPhysicalPrimary()) {
          auto pdgparticle = pdg->GetParticle(mcParticle.pdgCode());
          if (pdgparticle != nullptr) {
            if (std::abs(pdgparticle->Charge()) < 3)
              continue;
          }
          histos.fill(HIST("MCGENetaHistogram"), mcParticle.eta());
          histos.fill(HIST("MCGENEtaZvtxTracks"), mcParticle.eta(), mcCollision.posZ());
        }
      }
      histos.fill(HIST("MCGENMultiplicity"), MCparticleCounter);
    }
    if (random == 1) {
      histos.fill(HIST("MCGENFakeDataeventCounter"), 0.5);

      for (auto& mcParticle : mcParticles) {
        if (mcParticle.isPhysicalPrimary()) {
          auto pdgparticle = pdg->GetParticle(mcParticle.pdgCode());
          if (pdgparticle != nullptr) {
            if (std::abs(pdgparticle->Charge()) < 3)
              continue;
          }
          histos.fill(HIST("MCGENFakeDataetaHistogram"), mcParticle.eta());
        }
      }
    }
  }
  PROCESS_SWITCH(mcdata, processMCGEN, "process for GEN MC data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<mcdata>(cfgc)};
}
