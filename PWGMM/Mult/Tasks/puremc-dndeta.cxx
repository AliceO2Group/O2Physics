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
#include "Framework/StaticFor.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include <TDatabasePDG.h>
#include <TPDGCode.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

AxisSpec ZAxis = {301, -30.1, 30.1};
AxisSpec DeltaZAxis = {61, -6.1, 6.1};
AxisSpec DCAAxis = {601, -3.01, 3.01};
AxisSpec EtaAxis = {62, -6.2, 6.2};
AxisSpec RapidityAxis = {102, -10.2, 10.2};
AxisSpec PhiAxis = {629, 0, 2 * M_PI};
AxisSpec PtAxis = {2401, -0.005, 24.005};
AxisSpec PtAxis_wide = {1041, -0.05, 104.05};
AxisSpec PtAxisEff = {{0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6,
                       1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0}};
AxisSpec ScaleAxis = {121, -0.5, 120.5};
AxisSpec MPIAxis = {51, -0.5, 50.5};
AxisSpec ProcAxis = {21, 89.5, 110.5};

auto static constexpr mincharge = 3.f;

static constexpr std::string_view species[] = {"pi", "p", "e", "K"};
static constexpr std::array<int, 4> speciesIds{kPiPlus, kProton, kElectron, kKPlus};

struct PureMcMultiplicityCounter {
  SliceCache cache;
  Preslice<aod::Tracks> perCol = aod::track::collisionId;
  Preslice<aod::McParticles> perMcCol = aod::mcparticle::mcCollisionId;
  Service<o2::framework::O2DatabasePDG> pdg;
  Configurable<float> etaRange{"eta-range", 1.0f, "Eta range to consider"};
  ConfigurableAxis multBinning{"multiplicity-binning", {301, -0.5, 300.5}, ""};

  HistogramRegistry registry{
    "registry",
    {
      {"MCEvents/Vertex/X", "; X (cm); events", {HistType::kTH1F, {ZAxis}}},
      {"MCEvents/Vertex/Y", "; Y (cm); events", {HistType::kTH1F, {ZAxis}}},
      {"MCEvents/Vertex/Z", "; Z (cm); events", {HistType::kTH1F, {ZAxis}}},
      {"Particles/Primaries/Pt", " ;p_{T} (GeV/c)", {HistType::kTH1F, {PtAxis_wide}}},
      {"Particles/Primaries/Eta", " ; #eta", {HistType::kTH1F, {EtaAxis}}},
      {"Particles/Primaries/Y", " ; y", {HistType::kTH1F, {RapidityAxis}}},
      {"Particles/Primaries/EtaZvtx", " ; #eta; Z_{vtx} (cm); particles", {HistType::kTH2F, {EtaAxis, ZAxis}}},
      {"MCEvents/Properties/ProcessID", "; process ID", {HistType::kTH1F, {ProcAxis}}},
      {"MCEvents/Properties/ScalePerProcessID", " ; scale (GeV); process ID", {HistType::kTH2F, {ScaleAxis, ProcAxis}}},
      {"MCEvents/Properties/NMPIsPerProcessID", " ; N_{MPI}; process ID", {HistType::kTH2F, {MPIAxis, ProcAxis}}},
      {"MCEvents/Properties/ScaleVsNMPIsPerProcessID", " ;scale (GeV); N_{MPI}; process ID", {HistType::kTHnSparseF, {MPIAxis, ScaleAxis, ProcAxis}}},
    } //
  };

  void init(InitContext const&)
  {
    for (auto i = 0u; i < speciesIds.size(); ++i) {
      registry.add({(std::string("Particles/Primaries/") + std::string(species[i]) + "/Pt").c_str(), " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxis_wide}}});
      registry.add({(std::string("Particles/Primaries/") + std::string(species[i]) + "/Eta").c_str(), " ; #eta", {HistType::kTH1F, {EtaAxis}}});
      registry.add({(std::string("Particles/Primaries/") + std::string(species[i]) + "/Y").c_str(), " ; y", {HistType::kTH1F, {RapidityAxis}}});
    }

    AxisSpec MultAxis = {multBinning};
    registry.add({"MCEvents/NtrkZvtx", " ; N_{trk}; Z_{vtx} (cm); events", {HistType::kTH2F, {MultAxis, ZAxis}}});
    registry.add({"MCEvents/Properties/MultiplicityPerProcessID", " ; N_{p}; process ID; events", {HistType::kTH2F, {MultAxis, ProcAxis}}});
    registry.add({"MCEvents/Properties/MultiplicityVsScalePerProcessID", " ; scale; N_{p}; process ID", {HistType::kTHnSparseF, {ScaleAxis, MultAxis, ProcAxis}}});
    registry.add({"MCEvents/Properties/MultiplicityVsMPIsPerProcessID", " ; N_{MPI}; N_{p}; process ID", {HistType::kTHnSparseF, {MPIAxis, MultAxis, ProcAxis}}});

    if (doprocessReco) {
      registry.add({"Events/Properties/Multiplicity", " ; N_{p}; events", {HistType::kTH1F, {MultAxis}}});
      registry.add({"Events/Vertex/X", "; X (cm); events", {HistType::kTH1F, {ZAxis}}});
      registry.add({"Events/Vertex/Y", "; Y (cm); events", {HistType::kTH1F, {ZAxis}}});
      registry.add({"Events/Vertex/Z", "; Z (cm); events", {HistType::kTH1F, {ZAxis}}});
      registry.add({"Tracks/Pt", " ;p_{T} (GeV/c)", {HistType::kTH1F, {PtAxis_wide}}});
      registry.add({"Tracks/Eta", " ; #eta", {HistType::kTH1F, {EtaAxis}}});
      registry.add({"Tracks/DCAXY", "; DCA_{XY} (cm)", {HistType::kTH1F, {DCAAxis}}});
      registry.add({"Tracks/DCAZ", "; DCA_{Z} (cm)", {HistType::kTH1F, {DCAAxis}}});
      registry.add({"Tracks/EtaZvtx", " ; #eta; Z_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}});
      registry.add({"Events/NtrkZvtx", " ; N_{trk}; Z_{vtx} (cm); events", {HistType::kTH2F, {MultAxis, ZAxis}}});
    }
    if (doprocessResponse) {
      registry.add({"MCEvents/VertexCorrelation", " ; Z_{vtx}^{gen} (cm); Z_{vtx}^{rec} (cm)", {HistType::kTH2F, {ZAxis, ZAxis}}});
      registry.add({"MCEvents/Efficiency", "", {HistType::kTH1F, {{6, -0.5, 5.5}}}});
      registry.add({"MCEvents/Response", " ; N_{gen}; N_{rec}", {HistType::kTH2F, {MultAxis, MultAxis}}});
      auto eff = registry.get<TH1>(HIST("MCEvents/Efficiency"));
      auto* x = eff->GetXaxis();
      x->SetBinLabel(1, "Generated");
      x->SetBinLabel(2, "Generate INEL>0");
      x->SetBinLabel(3, "Reconstructed");
      x->SetBinLabel(4, "Reconstructed INEL>0");
      x->SetBinLabel(5, "Selected");
      x->SetBinLabel(6, "Selected INEL>0");

      registry.add({"MCEvents/EfficiencyMultiplicityN", " ; N_{gen}", {HistType::kTH1F, {MultAxis}}});
      registry.add({"MCEvents/EfficiencyMultiplicityNExtra", " ; N_{gen}", {HistType::kTH1F, {MultAxis}}});
      registry.add({"MCEvents/EfficiencyMultiplicityD", " ; N_{gen}", {HistType::kTH1F, {MultAxis}}});
    }
    if (doprocessEfficiency) {
      registry.add({"Particles/Primaries/EfficiencyN", " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxisEff}}});
      registry.add({"Particles/Primaries/EfficiencyD", " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxisEff}}});
      registry.add({"Particles/Secondaries/EfficiencyN", " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxisEff}}});
      registry.add({"Particles/Primaries/PtCorrelation", " ; p_{T}^{gen} (GeV/c); p_{T}^{rec} (GeV/c)", {HistType::kTH2F, {PtAxis_wide, PtAxis_wide}}});
      registry.add({"Particles/Primaries/EtaCorrelation", " ; #eta^{gen}; #eta^{rec}", {HistType::kTH2F, {EtaAxis, EtaAxis}}});
      registry.add({"Particles/Secondaries/PtCorrelation", " ; p_{T}^{gen} (GeV/c); p_{T}^{rec} (GeV/c)", {HistType::kTH2F, {PtAxis_wide, PtAxis_wide}}});
      registry.add({"Particles/Secondaries/EtaCorrelation", " ; #eta^{gen}; #eta^{rec}", {HistType::kTH2F, {EtaAxis, EtaAxis}}});

      for (auto i = 0u; i < speciesIds.size(); ++i) {
        registry.add({(std::string("Particles/Primaries/") + std::string(species[i]) + "/EfficiencyN").c_str(), " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxis_wide}}});
        registry.add({(std::string("Particles/Primaries/") + std::string(species[i]) + "/EfficiencyD").c_str(), " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxis_wide}}});
      }
    }
  }

  void process(soa::Join<aod::McCollisions, aod::HepMCXSections, aod::HepMCPdfInfos>::iterator const& collision, aod::McParticles const& particles)
  {
    registry.fill(HIST("MCEvents/Vertex/X"), collision.posX());
    registry.fill(HIST("MCEvents/Vertex/Y"), collision.posY());
    registry.fill(HIST("MCEvents/Vertex/Z"), collision.posZ());

    registry.fill(HIST("MCEvents/Properties/ProcessID"), collision.processId());
    registry.fill(HIST("MCEvents/Properties/ScalePerProcessID"), collision.ptHard(), collision.processId());
    registry.fill(HIST("MCEvents/Properties/NMPIsPerProcessID"), collision.nMPI(), collision.processId());
    registry.fill(HIST("MCEvents/Properties/ScaleVsNMPIsPerProcessID"), collision.ptHard(), collision.nMPI(), collision.processId());

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
      registry.fill(HIST("Particles/Primaries/Eta"), particle.eta());
      registry.fill(HIST("Particles/Primaries/Y"), particle.y());
      registry.fill(HIST("Particles/Primaries/EtaZvtx"), particle.eta(), collision.posZ());

      static_for<0, 3>(
        [&](auto idx) {
          constexpr int i = idx.value;
          if (particle.pdgCode() == speciesIds[i]) {
            registry.fill(HIST("Particles/Primaries/") + HIST(species[i]) + HIST("/Eta"), particle.eta());
            registry.fill(HIST("Particles/Primaries/") + HIST(species[i]) + HIST("/Y"), particle.y());
          }
        });

      if (std::abs(particle.eta()) >= etaRange) {
        continue;
      }
      ++Np;
      registry.fill(HIST("Particles/Primaries/Pt"), particle.pt());
      static_for<0, 3>(
        [&](auto idx) {
          constexpr int i = idx.value;
          if (particle.pdgCode() == speciesIds[i]) {
            registry.fill(HIST("Particles/Primaries/") + HIST(species[i]) + HIST("/Pt"), particle.pt());
          }
        });
    }
    registry.fill(HIST("MCEvents/Properties/MultiplicityVsMPIsPerProcessID"), collision.nMPI(), Np, collision.processId());
    registry.fill(HIST("MCEvents/Properties/MultiplicityVsScalePerProcessID"), collision.ptHard(), Np, collision.processId());
    registry.fill(HIST("MCEvents/Properties/MultiplicityPerProcessID"), Np, collision.processId());
    registry.fill(HIST("MCEvents/NtrkZvtx"), Np, collision.posZ());
  }

  void processReco(aod::Collision const& collision, soa::Join<aod::Tracks, aod::TracksDCA> const& tracks)
  {
    registry.fill(HIST("Events/Vertex/X"), collision.posX());
    registry.fill(HIST("Events/Vertex/Y"), collision.posY());
    registry.fill(HIST("Events/Vertex/Z"), collision.posZ());
    auto Ntrk = 0;
    for (auto& track : tracks) {
      registry.fill(HIST("Tracks/DCAXY"), track.dcaXY());
      registry.fill(HIST("Tracks/DCAZ"), track.dcaZ());
      registry.fill(HIST("Tracks/Pt"), track.pt());
      registry.fill(HIST("Tracks/Eta"), track.eta());
      registry.fill(HIST("Tracks/EtaZvtx"), track.eta(), collision.posZ());
      if (std::abs(track.eta()) < etaRange) {
        ++Ntrk;
      }
    }
    registry.fill(HIST("Events/Properties/Multiplicity"), Ntrk);
    registry.fill(HIST("Events/NtrkZvtx"), Ntrk, collision.posZ());
  }

  PROCESS_SWITCH(PureMcMultiplicityCounter, processReco, "Process smeared tracks", false);

  void processResponse(aod::McCollision const& mccollision, soa::SmallGroups<soa::Join<aod::Collisions, aod::McCollisionLabels>> const& collisions, aod::McParticles const& particles, soa::Join<aod::Tracks, aod::TracksDCA, aod::McTrackLabels> const& tracks)
  {
    registry.fill(HIST("MCEvents/Efficiency"), 0);
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
    registry.fill(HIST("MCEvents/EfficiencyMultiplicityD"), Np);
    if (Np > 0) {
      registry.fill(HIST("MCEvents/Efficiency"), 1);
    }
    if (collisions.size() > 0) {
      registry.fill(HIST("MCEvents/EfficiencyMultiplicityN"), Np);
    }
    if (collisions.size() > 1) {
      for (auto i = 1; i < collisions.size(); ++i) {
        registry.fill(HIST("MCEvents/EfficiencyMultiplicityNExtra"), Np);
      }
    }
    for (auto& collision : collisions) {
      auto Ntrk = 0;
      registry.fill(HIST("MCEvents/Efficiency"), 2);
      auto tracksample = tracks.sliceBy(perCol, collision.globalIndex());
      for (auto& track : tracksample) {
        if (std::abs(track.eta()) < etaRange) {
          ++Ntrk;
        }
      }
      if (Ntrk > 0) {
        registry.fill(HIST("MCEvents/Efficiency"), 3);
      }
      registry.fill(HIST("MCEvents/Response"), Np, Ntrk);
      registry.fill(HIST("MCEvents/VertexCorrelation"), mccollision.posZ(), collision.posZ());
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
      registry.fill(HIST("Particles/Primaries/EfficiencyD"), particle.pt());
      static_for<0, 3>(
        [&](auto idx) {
          constexpr int i = idx.value;
          if (particle.pdgCode() == speciesIds[i]) {
            registry.fill(HIST("Particles/Primaries/") + HIST(species[i]) + HIST("/EfficiencyD"), particle.pt());
          }
        });
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
          registry.fill(HIST("Particles/Primaries/EfficiencyN"), particle.pt());
          registry.fill(HIST("Particles/Primaries/PtCorrelation"), particle.pt(), track.pt());
          registry.fill(HIST("Particles/Primaries/EtaCorrelation"), particle.eta(), track.eta());

          static_for<0, 3>(
            [&](auto idx) {
              constexpr int i = idx.value;
              if (particle.pdgCode() == speciesIds[i]) {
                registry.fill(HIST("Particles/Primaries/") + HIST(species[i]) + HIST("/EfficiencyN"), particle.pt());
              }
            });
        } else {
          registry.fill(HIST("Particles/Secondaries/EfficiencyN"), particle.pt());
          registry.fill(HIST("Particles/Secondaries/PtCorrelation"), particle.pt(), track.pt());
          registry.fill(HIST("Particles/Secondaries/EtaCorrelation"), particle.eta(), track.eta());
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
