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
#include "Framework/AnalysisDataModel.h"

#include "ReconstructionDataFormats/GlobalTrackID.h"

#include "Common/DataModel/EventSelection.h"

#include "Framework/O2DatabasePDGPlugin.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct multiplicityPbPb {

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
    const AxisSpec axisPt{nBinsPt, 0, 5, "p_T"};

    const AxisSpec axisCounter{1, 0, +1, ""};

    const AxisSpec axisDCAxy{nBinsDCA, -0.22, 0.22, "DCA_{xy} (cm)"};
    const AxisSpec axisDCAz{nBinsDCA, -0.22, 0.22, "DCA_{z} (cm)"};

    const AxisSpec axisNtrk{nBinsMult, 0, 300, "N_{trk}"};

    const AxisSpec axisPhi{nBinsPhi, -0.4, 6.8, "#phi"};
    const AxisSpec axisZvtx{nBinsZvtx, -30, 30, "Z_{vtx} (cm)"};

    histos.add("etaHistogram", "; ", kTH1D, {axisEta});
    histos.add("MCGENetaHistogram", "; ", kTH1D, {axisEta});
    histos.add("ptHistogram", "; ", kTH1D, {axisPt});
    histos.add("MCGENptHistogram", "; ", kTH1D, {axisPt});

    histos.add("eventCounter", "eventCounter", kTH1D, {axisCounter});
    histos.add("MCGENeventCounter", "eventCounter", kTH1D, {axisCounter});

    histos.add("DCAxy", "; DCA_{xy} (cm)", kTH1D, {axisDCAxy});
    histos.add("DCAz", "; DCA_{z} (cm)", kTH1D, {axisDCAz});

    histos.add("Multiplicity", "; tracks; events", kTH1D, {axisNtrk});
    histos.add("MCGENMultiplicity", "; tracks; events", kTH1D, {axisNtrk});

    histos.add("Multiplicity_01", "; tracks; events", kTH1D, {axisNtrk});         // |eta|<0.1
    histos.add("Multiplicity_01_10", "; tracks; events", kTH1D, {axisNtrk});      // |eta|<0.1, |Zvtx|<10 cm
    histos.add("MCGENMultiplicity_01", "; tracks; events", kTH1D, {axisNtrk});    // |eta|<0.1
    histos.add("MCGENMultiplicity_01_10", "; tracks; events", kTH1D, {axisNtrk}); // |eta|<0.1, |Zvtx|<10 cm

    histos.add("PhiTracks", "; #phi; tracks", kTH1D, {axisPhi});

    histos.add("ZvtxEvents", "; Z_{vtx} (cm); events", kTH1D, {axisZvtx});
    histos.add("ZvtxEvents_01", "; Z_{vtx} (cm); events", kTH1D, {axisZvtx}); // |eta|<0.1
    histos.add("MCGENZvtxEvents", "; Z_{vtx} (cm); events", kTH1D, {axisZvtx});
    histos.add("MCGENZvtxEvents_01", "; Z_{vtx} (cm); events", kTH1D, {axisZvtx}); // |eta|<0.1

    histos.add("EtaZvtxTracks", "; #eta; Z_{vtx} (cm); tracks", kTH2D, {axisEta, axisZvtx});
    histos.add("NtrkZvtxEvents", "; N_{trk}; Z_{vtx} (cm); events", kTH2D, {axisNtrk, axisZvtx});

    histos.add("MCGENEtaZvtxTracks", "; #eta; Z_{vtx} (cm); tracks", kTH2D, {axisEta, axisZvtx});
    histos.add("MCGENNtrkZvtxEvents", "; N_{trk}; Z_{vtx} (cm); events", kTH2D, {axisNtrk, axisZvtx});

    histos.add("PhiEtaTracks", "; #phi; #eta; tracks", kTH2D, {axisPhi, axisEta});
    histos.add("MCGENPhiEtaTracks", "; #phi; #eta; tracks", kTH2D, {axisPhi, axisEta});
  }

  // void process(aod::Collision const& collision, soa::Filtered<myCompleteTracks> const& tracks, aod::McParticles const&)
  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<myCompleteTracks> const& tracks)
  {
    if (!collision.sel8()) {
      return;
    }

    int trackCounter = 0;
    int trackCounter_01 = 0;
    int trackCounter_01_10 = 0;

    histos.fill(HIST("eventCounter"), 0.5);

    histos.fill(HIST("ZvtxEvents"), collision.posZ());

    for (auto& track : tracks) {
      ++trackCounter;

      histos.fill(HIST("etaHistogram"), track.eta());
      histos.fill(HIST("ptHistogram"), track.pt());

      histos.fill(HIST("DCAxy"), track.dcaXY());
      histos.fill(HIST("DCAz"), track.dcaZ());

      histos.fill(HIST("PhiTracks"), track.phi());

      histos.fill(HIST("EtaZvtxTracks"), track.eta(), collision.posZ());
      histos.fill(HIST("PhiEtaTracks"), track.phi(), track.eta());

      if (std::abs(track.eta()) < 0.1)
        ++trackCounter_01;
      if (std::abs(track.eta()) < 0.1 && std::abs(collision.posZ()) < 10)
        ++trackCounter_01_10;
    }

    histos.fill(HIST("Multiplicity"), trackCounter);
    histos.fill(HIST("Multiplicity_01"), trackCounter_01); // to be deleted
    histos.fill(HIST("Multiplicity_01_10"), trackCounter_01_10);

    histos.fill(HIST("NtrkZvtxEvents"), trackCounter, collision.posZ());
  }

  void processMCGEN(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles)
  {
    int MCparticleCounter = 0;
    int MCparticleCounter_01 = 0;
    int MCparticleCounter_01_10 = 0;

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
        histos.fill(HIST("MCGENptHistogram"), mcParticle.pt());

        histos.fill(HIST("MCGENEtaZvtxTracks"), mcParticle.eta(), mcCollision.posZ());
        histos.fill(HIST("MCGENPhiEtaTracks"), mcParticle.phi(), mcParticle.eta());

        if (std::abs(mcParticle.eta()) < 0.1)
          ++MCparticleCounter_01;
        if (std::abs(mcParticle.eta()) < 0.1 && std::abs(mcCollision.posZ()) < 10)
          ++MCparticleCounter_01_10;
      }
    }
    histos.fill(HIST("MCGENMultiplicity"), MCparticleCounter);
    histos.fill(HIST("MCGENMultiplicity_01"), MCparticleCounter_01);
    histos.fill(HIST("MCGENMultiplicity_01_10"), MCparticleCounter_01_10);

    histos.fill(HIST("MCGENNtrkZvtxEvents"), MCparticleCounter, mcCollision.posZ());
  }
  PROCESS_SWITCH(multiplicityPbPb, processMCGEN, "process for GEN MC data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<multiplicityPbPb>(cfgc)};
}
