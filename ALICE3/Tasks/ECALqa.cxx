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
/// \file   ECALqa.cxx
/// \author Nicolo' Jacazio
/// \since  14/09/2021
/// \brief  Task to use the ALICE3 ECAL table
///

// O2 includes
#include "ALICE3/DataModel/ECAL.h"

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/PID.h"

using namespace o2;
using namespace o2::track;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{

namespace indices
{
DECLARE_SOA_INDEX_COLUMN(Track, track);
DECLARE_SOA_INDEX_COLUMN(ECAL, ecal);

DECLARE_SOA_INDEX_COLUMN(McParticle, mcparticle);

} // namespace indices

DECLARE_SOA_INDEX_TABLE_USER(ECALTracksIndex, Tracks, "ECALTRK", indices::TrackId, indices::ECALId);
DECLARE_SOA_INDEX_TABLE_USER(ECALMcPartIndex, McParticles, "ECALPART", indices::McParticleId, indices::ECALId);
} // namespace o2::aod

struct ecalIndexBuilder { // Builder of the ECAL-track index linkage
  Builds<o2::aod::ECALTracksIndex> ind;
  Builds<o2::aod::ECALMcPartIndex> indPart;
  void init(o2::framework::InitContext&)
  {
  }
};

struct ecalQaMc {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::QAObject};

  void init(o2::framework::InitContext&)
  {
    histos.add("energy", ";Energy (GeV/#it{c});Entries", HistType::kTH1F, {{100, 0, 100}});
    histos.add("energyVsp", ";Energy (GeV/#it{c});#it{p} (GeV/#it{c});Entries", HistType::kTH2F, {{100, 0, 100}, {100, 0, 100}});
    histos.add("ecalpVsp", ";#it{p}_{ECAL} (GeV/#it{c});#it{p} (GeV/#it{c});Entries", HistType::kTH2F, {{100, 0, 100}, {100, 0, 100}});
    histos.add("ecalpxVspx", ";#it{p}_{ECAL x} (GeV/#it{c});#it{p} (GeV/#it{c});Entries", HistType::kTH2F, {{100, 0, 100}, {100, 0, 100}});
    histos.add("ecalpyVspy", ";#it{p}_{ECAL y} (GeV/#it{c});#it{p} (GeV/#it{c});Entries", HistType::kTH2F, {{100, 0, 100}, {100, 0, 100}});
    histos.add("ecalpzVspz", ";#it{p}_{ECAL z} (GeV/#it{c});#it{p} (GeV/#it{c});Entries", HistType::kTH2F, {{100, 0, 100}, {100, 0, 100}});
    histos.add("energyVsecalp", ";Energy (GeV/#it{c});#it{p}_{ECAL} (GeV/#it{c});Entries", HistType::kTH2F, {{100, 0, 100}, {100, 0, 100}});
    histos.add("energyVsecalpx", ";Energy (GeV/#it{c});#it{p}_{ECAL x} (GeV/#it{c});Entries", HistType::kTH2F, {{100, 0, 100}, {100, 0, 100}});
    histos.add("energyVsecalpy", ";Energy (GeV/#it{c});#it{p}_{ECAL y} (GeV/#it{c});Entries", HistType::kTH2F, {{100, 0, 100}, {100, 0, 100}});
    histos.add("energyVsecalpz", ";Energy (GeV/#it{c});#it{p}_{ECAL z} (GeV/#it{c});Entries", HistType::kTH2F, {{100, 0, 100}, {100, 0, 100}});
    histos.add("PDGs", "Particle PDGs;PDG Code", kTH1D, {{100, 0.f, 100.f}});
  }

  using Trks = soa::Join<aod::Tracks, aod::ECALTracksIndex, aod::TracksExtra>;
  void process(const soa::Join<aod::McParticles, aod::ECALMcPartIndex>& mcParticles,
               const Trks&,
               const aod::McTrackLabels&,
               const aod::ECALs&,
               const aod::Collisions&)
  {
    for (auto& particle : mcParticles) {
      if (!particle.has_ecal())
        continue;
      histos.fill(HIST("energy"), particle.ecal().e());
      const float ecalp = std::sqrt(particle.ecal().px() * particle.ecal().px() + particle.ecal().py() * particle.ecal().py() + particle.ecal().pz() * particle.ecal().pz());
      histos.fill(HIST("energyVsp"), particle.ecal().e(), particle.p());
      histos.fill(HIST("ecalpVsp"), ecalp, particle.p());
      histos.fill(HIST("ecalpxVspx"), particle.ecal().px(), particle.px());
      histos.fill(HIST("ecalpyVspy"), particle.ecal().py(), particle.py());
      histos.fill(HIST("ecalpzVspz"), particle.ecal().pz(), particle.pz());
      histos.fill(HIST("energyVsecalp"), particle.ecal().e(), ecalp);
      histos.fill(HIST("energyVsecalpx"), particle.ecal().e(), particle.ecal().px());
      histos.fill(HIST("energyVsecalpy"), particle.ecal().e(), particle.ecal().py());
      histos.fill(HIST("energyVsecalpz"), particle.ecal().e(), particle.ecal().pz());
      histos.get<TH1>(HIST("PDGs"))->Fill(Form("%i", particle.pdgCode()), 1.f);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  auto workflow = WorkflowSpec{adaptAnalysisTask<ecalIndexBuilder>(cfg)};
  workflow.push_back(adaptAnalysisTask<ecalQaMc>(cfg));
  return workflow;
}
