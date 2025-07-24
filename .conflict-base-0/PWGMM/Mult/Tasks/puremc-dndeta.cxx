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

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include <TDatabasePDG.h>
#include <TPDGCode.h>

#include "Axes.h"
#include "Histograms.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace pwgmm::mult;
using namespace pwgmm::mult::histograms;

auto static constexpr mincharge = 3.f;

struct PureMcMultiplicityCounter {
  SliceCache cache;
  Preslice<aod::Tracks> perCol = aod::track::collisionId;
  Preslice<aod::McParticles> perMcCol = aod::mcparticle::mcCollisionId;
  Service<o2::framework::O2DatabasePDG> pdg;
  Configurable<float> etaRange{"eta-range", 1.0f, "Eta range to consider for multiplicity definition"};
  ConfigurableAxis multBinning{"multiplicity-binning", {301, -0.5, 300.5}, ""};

  Configurable<bool> useProcessId{"use-process-id", true, "Use generator-provided process ID"};
  Configurable<bool> useScale{"use-event-scale", true, "Use generator-provided event scale"};
  Configurable<bool> useNMPI{"use-n-mpi", true, "Use generator-provided N MPIs"};

  HistogramRegistry registry{
    "registry",
    {
      {MCVertex.data(), "; X (cm) ; Y (cm) ; Z (cm) ; events", {HistType::kTHnSparseF, {ZAxis, ZAxis, ZAxis}}},
      {MCParticles.data(), " ; p_{T} (GeV/c); #eta ; Y ; Z_{vtx} (cm) ; species", {HistType::kTHnSparseF, {PtAxis_wide, EtaAxis, RapidityAxis, ZAxis, SpeciesAxis}}},
    } //
  };

  std::shared_ptr<THnSparse> mcCor = nullptr;
  std::shared_ptr<THnSparse> mcRes = nullptr;
  std::vector<double> fillVector;

  void init(InitContext const&)
  {
    auto histo = registry.get<THnSparse>(HIST(MCParticles.data()));
    auto* axis = histo->GetAxis(4);
    for (auto i = 0U; i < speciesIds.size(); ++i) {
      axis->SetBinLabel(i + 1, species[i].data());
    }
    axis->SetBinLabel(5, "Other");

    AxisSpec MultAxis = {multBinning};
    std::vector<AxisSpec> correlateAxes{{MultAxis, ZAxis}};
    std::vector<AxisSpec> responseAxes{{MultAxis, MultAxis, ZAxis}};
    std::string correlateLabels = " ; N_{p} ; Z_{vtx} (cm)";
    std::string responseLabels = " ; N_{p} ; N_{trk} ; Z_{vtx} (cm)";
    if (useProcessId) {
      correlateAxes.emplace_back(ProcAxis);
      correlateLabels += " ; process ID";
      responseAxes.emplace_back(ProcAxis);
      responseLabels += " ; process ID";
    }
    if (useScale) {
      correlateAxes.emplace_back(ScaleAxis);
      correlateLabels += " ; scale";
      responseAxes.emplace_back(ScaleAxis);
      responseLabels += " ; scale";
    }
    if (useNMPI) {
      correlateAxes.emplace_back(MPIAxis);
      correlateLabels += " ; N_{MPI}";
      responseAxes.emplace_back(MPIAxis);
      responseLabels += " ; N_{MPI}";
    }

    auto h = registry.add(MCCorrelates.data(), correlateLabels.data(), {HistType::kTHnSparseF, correlateAxes});
    mcCor = std::get<1>(h);

    if (doprocessReco) {
      registry.add({RecoCorrelates.data(), " ; N_{trk}; Z_{vtx} (cm); events", {HistType::kTHnSparseF, {MultAxis, ZAxis}}});
      registry.add({RecoVertex.data(), "; X (cm) ; Y (cm) ; Z (cm) ; events", {HistType::kTHnSparseF, {ZAxis, ZAxis, ZAxis}}});
      registry.add({RecoTracks.data(), " ; p_{T} (GeV/c); #eta ; DCA_{XY} (cm) ; DCA_{Z} (cm) ; Z_{vtx} (cm)", {HistType::kTHnSparseF, {PtAxis_wide, EtaAxis, DCAAxis, DCAAxis, ZAxis}}});
    }
    if (doprocessResponse) {
      auto hl = registry.add({MCResponse.data(), responseLabels.data(), {HistType::kTHnSparseF, responseAxes}});
      mcRes = std::get<1>(hl);

      registry.add({MCEfficiency.data(), "", {HistType::kTH1F, {{static_cast<int>(EvEffBins::kSelectedPVgt0), 0.5, static_cast<float>(EvEffBins::kSelectedPVgt0) + 0.5}}}});
      auto eff = registry.get<TH1>(HIST(MCEfficiency.data()));
      auto* x = eff->GetXaxis();
      x->SetBinLabel(static_cast<int>(EvEffBins::kGen), EvEffBinLabels[static_cast<int>(EvEffBins::kGen)].data());
      x->SetBinLabel(static_cast<int>(EvEffBins::kGengt0), EvEffBinLabels[static_cast<int>(EvEffBins::kGengt0)].data());
      x->SetBinLabel(static_cast<int>(EvEffBins::kRec), EvEffBinLabels[static_cast<int>(EvEffBins::kRec)].data());
      x->SetBinLabel(static_cast<int>(EvEffBins::kSelected), EvEffBinLabels[static_cast<int>(EvEffBins::kSelected)].data());
      x->SetBinLabel(static_cast<int>(EvEffBins::kSelectedgt0), EvEffBinLabels[static_cast<int>(EvEffBins::kSelectedgt0)].data());
      x->SetBinLabel(static_cast<int>(EvEffBins::kSelectedPVgt0), EvEffBinLabels[static_cast<int>(EvEffBins::kSelectedPVgt0)].data());
    }
    if (doprocessEfficiency) {
      registry.add({MCParticlesEfficiencyN.data(), " ; p_{T} (GeV/c); #eta ; Y ; Z_{vtx} (cm) ; species", {HistType::kTHnSparseF, {PtAxis_wide, EtaAxis, RapidityAxis, ZAxis, SpeciesAxis}}});
      histo = registry.get<THnSparse>(HIST(MCParticlesEfficiencyN.data()));
      axis = histo->GetAxis(4);
      for (auto i = 0U; i < speciesIds.size(); ++i) {
        axis->SetBinLabel(i + 1, species[i].data());
      }
      axis->SetBinLabel(5, "Other");
      registry.add({MCParticlesEfficiencyD.data(), " ; p_{T} (GeV/c); #eta ; Y ; Z_{vtx} (cm) ; species", {HistType::kTHnSparseF, {PtAxis_wide, EtaAxis, RapidityAxis, ZAxis, SpeciesAxis}}});
      histo = registry.get<THnSparse>(HIST(MCParticlesEfficiencyD.data()));
      axis = histo->GetAxis(4);
      for (auto i = 0U; i < speciesIds.size(); ++i) {
        axis->SetBinLabel(i + 1, species[i].data());
      }
      axis->SetBinLabel(5, "Other");
    }
  }

  template <typename C>
  void fillMcCorrelates(C const& collision, int Np)
  {
    fillVector.clear();
    fillVector.push_back(Np);
    fillVector.push_back(collision.posZ());
    if (useProcessId) {
      fillVector.push_back(collision.processId());
    }
    if (useScale) {
      fillVector.push_back(collision.ptHard());
    }
    if (useNMPI) {
      fillVector.push_back(collision.nMPI());
    }
    mcCor->Fill(fillVector.data());
  }

  void process(soa::Join<aod::McCollisions, aod::HepMCXSections, aod::HepMCPdfInfos>::iterator const& collision, aod::McParticles const& particles)
  {
    registry.fill(HIST(MCVertex.data()), collision.posX(), collision.posY(), collision.posZ());

    auto Np = 0;
    for (auto const& particle : particles) {
      if (!particle.isPhysicalPrimary()) {
        continue;
      }
      if (!particle.producedByGenerator()) {
        continue;
      }
      auto charge = 0.;
      auto* p = pdg->GetParticle(particle.pdgCode());
      if (p != nullptr) {
        charge = p->Charge();
      }
      if (std::abs(charge) < mincharge) {
        continue;
      }

      auto pos = std::distance(speciesIds.begin(), std::find(speciesIds.begin(), speciesIds.end(), particle.pdgCode())) + 1;
      registry.fill(HIST(MCParticles.data()), particle.pt(), particle.eta(), particle.y(), collision.posZ(), pos);

      if (std::abs(particle.eta()) >= etaRange) {
        continue;
      }
      ++Np;
    }

    fillMcCorrelates(collision, Np);
  }

  void processReco(aod::Collision const& collision, soa::Join<aod::Tracks, aod::TracksDCA> const& tracks)
  {
    registry.fill(HIST(RecoVertex.data()), collision.posX(), collision.posY(), collision.posZ());
    auto Ntrk = 0;
    for (auto& track : tracks) {
      registry.fill(HIST(RecoTracks.data()), track.pt(), track.eta(), track.dcaXY(), track.dcaZ(), collision.posZ());
      if (std::abs(track.eta()) < etaRange) {
        ++Ntrk;
      }
    }
    registry.fill(HIST(RecoCorrelates.data()), Ntrk, collision.posZ());
  }

  PROCESS_SWITCH(PureMcMultiplicityCounter, processReco, "Process smeared tracks", false);

  template <typename C>
  void fillMcResponse(C const& collision, int Np, int Ntrk)
  {
    fillVector.clear();
    fillVector.push_back(Np);
    fillVector.push_back(Ntrk);
    fillVector.push_back(collision.posZ());
    if (useProcessId) {
      fillVector.push_back(collision.processId());
    }
    if (useScale) {
      fillVector.push_back(collision.ptHard());
    }
    if (useNMPI) {
      fillVector.push_back(collision.nMPI());
    }
    mcRes->Fill(fillVector.data());
  }

  void processResponse(soa::Join<aod::McCollisions, aod::HepMCXSections, aod::HepMCPdfInfos>::iterator const& mccollision, soa::SmallGroups<soa::Join<aod::Collisions, aod::McCollisionLabels>> const& collisions, aod::McParticles const& particles, soa::Join<aod::Tracks, aod::TracksDCA, aod::McTrackLabels> const& tracks)
  {
    registry.fill(HIST(MCEfficiency.data()), 0);
    auto Np = 0;
    for (auto const& particle : particles) {
      if (!particle.isPhysicalPrimary()) {
        continue;
      }
      if (!particle.producedByGenerator()) {
        continue;
      }
      auto charge = 0.;
      auto* p = pdg->GetParticle(particle.pdgCode());
      if (p != nullptr) {
        charge = p->Charge();
      }
      if (std::abs(charge) < mincharge) {
        continue;
      }
      if (std::abs(particle.eta()) >= etaRange) {
        continue;
      }
      ++Np;
    }

    for (auto& collision : collisions) {
      auto Ntrk = 0;
      registry.fill(HIST(MCEfficiency.data()), 2);
      auto tracksample = tracks.sliceBy(perCol, collision.globalIndex());
      for (auto& track : tracksample) {
        if (std::abs(track.eta()) < etaRange) {
          ++Ntrk;
        }
      }
      if (Ntrk > 0) {
        registry.fill(HIST(MCEfficiency.data()), 3);
      }
      fillMcResponse(mccollision, Np, Ntrk);
    }
  }

  PROCESS_SWITCH(PureMcMultiplicityCounter, processResponse, "Process response", false);

  void processEfficiency(soa::Join<aod::Collisions, aod::McCollisionLabels>::iterator const& collision, aod::McCollisions const&, aod::McParticles const& particles, soa::Join<aod::Tracks, aod::McTrackLabels> const& tracks)
  {
    if (!collision.has_mcCollision()) {
      return;
    }
    auto psample = particles.sliceBy(perMcCol, collision.mcCollisionId());

    for (auto& particle : psample) {
      if (!particle.isPhysicalPrimary()) {
        continue;
      }
      if (!particle.producedByGenerator()) {
        continue;
      }
      auto charge = 0.;
      auto p = pdg->GetParticle(particle.pdgCode());
      if (p != nullptr) {
        charge = p->Charge();
      }
      if (std::abs(charge) < 3.) {
        continue;
      }
      if (std::abs(particle.eta()) >= etaRange) {
        continue;
      }
      auto pos = std::distance(speciesIds.begin(), std::find(speciesIds.begin(), speciesIds.end(), particle.pdgCode())) + 1;
      registry.fill(HIST(MCParticlesEfficiencyD.data()), particle.pt(), particle.eta(), particle.y(), collision.posZ(), pos);
    }

    for (auto& track : tracks) {
      if (std::abs(track.eta()) >= etaRange) {
        continue;
      }
      if (track.has_mcParticle()) {
        auto particle = track.mcParticle();
        if (particle.mcCollisionId() != collision.mcCollisionId()) {
          continue;
        }
        if (particle.isPhysicalPrimary()) {
          auto pos = std::distance(speciesIds.begin(), std::find(speciesIds.begin(), speciesIds.end(), particle.pdgCode())) + 1;
          registry.fill(HIST(MCParticlesEfficiencyN.data()), particle.pt(), particle.eta(), particle.y(), collision.posZ(), pos);
        }
      }
    }
  }

  PROCESS_SWITCH(PureMcMultiplicityCounter, processEfficiency, "Process efficiency", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return {adaptAnalysisTask<PureMcMultiplicityCounter>(cfgc)};
}
