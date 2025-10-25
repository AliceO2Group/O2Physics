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

/// \file treeCreatorElectronDCA.cxx
/// \brief Basic electron DCA analysis task
///
/// \author Martin Voelkl <martin.andreas.volkl@cern.ch>, University of Birmingham

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/AliasTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TPDGCode.h>

#include <cstdlib>

using namespace o2;
using namespace o2::aod;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace hf_ele_mc_red
{
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(SourcePdg, sourcePdg, int);
DECLARE_SOA_COLUMN(DcaXY, dcaXY, float);
DECLARE_SOA_COLUMN(ProductionRadius, productionRadius, float);

} // namespace hf_ele_mc_red
DECLARE_SOA_TABLE(HFeleMCRedTable, "AOD", "HFELERED",
                  hf_ele_mc_red::Eta, hf_ele_mc_red::Phi, hf_ele_mc_red::Pt, hf_ele_mc_red::SourcePdg, hf_ele_mc_red::DcaXY, hf_ele_mc_red::ProductionRadius);
} // namespace o2::aod

/// Electron DCA analysis task
struct HfTreeCreatorElectronDCA {
  Produces<o2::aod::HFeleMCRedTable> hfEleTable;

  Configurable<float> etaRange{"etaRange", 0.5, "pseudorapidity range"};
  Configurable<float> pTMin{"pTMin", 0.5, "min pT"};

  HfHelper hfHelper;
  Service<o2::framework::O2DatabasePDG> pdg;

  using TracksWExt = soa::Join<o2::aod::Tracks, o2::aod::TracksExtra, o2::aod::TracksDCA, aod::TrackSelection, o2::aod::TrackSelectionExtension, aod::TracksPidPi, aod::TracksPidKa>;
  using TracksWExtMc = soa::Join<o2::aod::Tracks, o2::aod::TracksExtra, o2::aod::TracksDCA, aod::TrackSelection, o2::aod::TrackSelectionExtension, aod::TracksPidPi, aod::TracksPidKa, McTrackLabels>;

  HistogramRegistry registry{
    "registry",
    {{"hZVertex", "z Vertex;z_{vtx};counts", {HistType::kTH1F, {{100, -20., 20.}}}},
     {"hpTTracks", "pT of tracks; p_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 10.}}}},
     {"hpTElectrons", "pT of electrons; p_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 10.}}}}}};

  void init(InitContext&)
  {
  }

  void processData(aod::Collisions::iterator const& collision)
  {
    registry.get<TH1>(HIST("hZVertex"))->Fill(collision.posZ());
  }
  PROCESS_SWITCH(HfTreeCreatorElectronDCA, processData, "Process Data", true);

  void processMc(aod::Collisions::iterator const& collision,
                 TracksWExtMc const& tracks,
                 aod::McParticles const&)
  {
    registry.get<TH1>(HIST("hZVertex"))->Fill(collision.posZ());
    int pdgCode = 0, absPDGCode = 0, sourcePDG = 0;
    for (const auto& track : tracks) {
      if (!track.trackCutFlagFb3()) {
        continue;
      }
      registry.get<TH1>(HIST("hpTTracks"))->Fill(track.pt());
      if (track.pt() < pTMin) {
        continue;
      }
      if (std::abs(track.eta()) > etaRange) {
        continue;
      }
      if (track.mcParticleId() < 1) {
        continue;
      }
      auto mcTrack = track.mcParticle();
      if (std::abs(mcTrack.pdgCode()) == kElectron) {
        bool isConversion = false;
        bool isBeauty = false;
        bool isCharm = false;
        double productionRadius = RecoDecay::sqrtSumOfSquares(mcTrack.vx(), mcTrack.vy());
        registry.get<TH1>(HIST("hpTElectrons"))->Fill(track.pt());
        auto motherTracks = mcTrack.mothers_as<aod::McParticles>();
        int numberOfMothers = motherTracks.size();
        // Categorise the electron sources
        int const firstMotherPDG = motherTracks[0].pdgCode();
        if (firstMotherPDG == kGamma) {
          isConversion = true;
        }
        while (numberOfMothers == 1) // loop through all generations
        {
          pdgCode = motherTracks[0].pdgCode();
          absPDGCode = std::abs(pdgCode);
          if ((absPDGCode / 100) == 4 || (absPDGCode / 1000) == 4) {
            isCharm = true;
            sourcePDG = pdgCode;
          }
          if ((absPDGCode / 100) == 5 || (absPDGCode / 1000) == 5) {
            isBeauty = true;
            sourcePDG = pdgCode; // already in order, since beauty would decay to charm
          }
          auto firstMother = motherTracks[0];
          motherTracks = firstMother.mothers_as<aod::McParticles>();
          numberOfMothers = motherTracks.size();
        }
        if (!isBeauty && !isCharm) {
          if (isConversion) {
            sourcePDG = kGamma;
          } else {
            sourcePDG = firstMotherPDG;
          }
        }
        hfEleTable(track.eta(), track.phi(), track.pt(), sourcePDG, track.dcaXY(), productionRadius);
      }
    }
  }
  PROCESS_SWITCH(HfTreeCreatorElectronDCA, processMc, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTreeCreatorElectronDCA>(cfgc)};
}
