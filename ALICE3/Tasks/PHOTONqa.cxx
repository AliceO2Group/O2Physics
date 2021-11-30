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
/// \file   PHOTONqa.cxx
/// \author  Ana Marin
/// \since  12/10/2021
/// \brief  Task to use the ALICE3 PHOTON table
///

// O2 includes
#include "Framework/AnalysisTask.h"
#include "ALICE3/DataModel/PHOTON.h"
#include "Common/Core/MC.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/PID/PIDResponse.h"
#include "ReconstructionDataFormats/PID.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::track;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace indices
{
DECLARE_SOA_INDEX_COLUMN(Track, track);
DECLARE_SOA_INDEX_COLUMN(Photon, photon);
DECLARE_SOA_INDEX_COLUMN(McParticle, mcParticle);
} // namespace indices
DECLARE_SOA_INDEX_TABLE_USER(PhotonTracksIndex, Tracks, "PhotonTrk", indices::TrackId, indices::PhotonId);
DECLARE_SOA_INDEX_TABLE_USER(PhotonMcPartIndex, McParticles, "PhotonPart", indices::McParticleId, indices::PhotonId);
} // namespace o2::aod

struct photonIndexBuilder { // Builder of the PHOTON-track index linkage
  Builds<o2::aod::PhotonTracksIndex> ind;
  Builds<o2::aod::PhotonMcPartIndex> indPart;

  void init(o2::framework::InitContext&)
  {
  }
};

struct photonQaMc {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::QAObject};

  void init(o2::framework::InitContext&)
  {
    histos.add("photonp", "; #it{p}_{PHOTON} (GeV/#it{c});Entries", HistType::kTH1F, {{100, 0, 100}});
    histos.add("photonpVsp", ";#it{p}_{PHOTON} (GeV/#it{c});#it{p} (GeV/#it{c});Entries", HistType::kTH2F, {{100, 0, 100}, {100, 0, 100}});
    histos.add("photonpxVspx", ";#it{p}_{PHOTON x} (GeV/#it{c});#it{p} (GeV/#it{c});Entries", HistType::kTH2F, {{100, 0, 100}, {100, 0, 100}});
    histos.add("photonpyVspy", ";#it{p}_{PHOTON y} (GeV/#it{c});#it{p} (GeV/#it{c});Entries", HistType::kTH2F, {{100, 0, 100}, {100, 0, 100}});
    histos.add("photonpzVspz", ";#it{p}_{PHOTON z} (GeV/#it{c});#it{p} (GeV/#it{c});Entries", HistType::kTH2F, {{100, 0, 100}, {100, 0, 100}});
    histos.add("PDGs", "Particle PDGs;PDG Code", kTH1D, {{100, 0.f, 100.f}});
  }

  using Trks = soa::Join<aod::Tracks, aod::PhotonTracksIndex, aod::TracksExtra>;
  void process(const soa::Join<aod::McParticles, aod::PhotonMcPartIndex>& mcParticles,
               const Trks& tracks,
               const aod::McTrackLabels& labels,
               const aod::Photons& photons,
               const aod::Collisions& colls)
  {
    for (auto& particle : mcParticles) {
      if (!particle.has_photon())
        continue;
      const float photonp = RecoDecay::sqrtSumOfSquares(particle.photon().px(), particle.photon().py(), particle.photon().pz());
      histos.fill(HIST("photonpVsp"), photonp, particle.p());
      histos.fill(HIST("photonpxVspx"), particle.photon().px(), particle.px());
      histos.fill(HIST("photonpyVspy"), particle.photon().py(), particle.py());
      histos.fill(HIST("photonpzVspz"), particle.photon().pz(), particle.pz());
      histos.get<TH1>(HIST("PDGs"))->Fill(Form("%i", particle.pdgCode()), 1.f);
    }
    for (auto& photon : photons) {
      const float photonp = RecoDecay::sqrtSumOfSquares(photon.px(), photon.py(), photon.pz());
      LOGF(debug, "Photon momentum %f", photonp);
      histos.fill(HIST("photonp"), photonp);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  auto workflow = WorkflowSpec{adaptAnalysisTask<photonIndexBuilder>(cfg)};
  workflow.push_back(adaptAnalysisTask<photonQaMc>(cfg));
  return workflow;
}
