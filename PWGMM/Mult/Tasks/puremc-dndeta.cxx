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
#include <TDatabasePDG.h>
#include <TPDGCode.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

AxisSpec ZAxis = {301, -30.1, 30.1};
AxisSpec DeltaZAxis = {61, -6.1, 6.1};
AxisSpec DCAAxis = {601, -3.01, 3.01};
AxisSpec EtaAxis = {22, -2.2, 2.2};
AxisSpec RapidityAxis = {102, -10.2, 10.2};
AxisSpec PhiAxis = {629, 0, 2 * M_PI};
AxisSpec PtAxis = {2401, -0.005, 24.005};
AxisSpec PtAxis_wide = {1041, -0.05, 104.05};
AxisSpec PtAxisEff = {{0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6,
                       1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0}};

auto static constexpr mincharge = 3.f;

static constexpr std::string_view species[] = {"pi", "p", "e", "K"};
static constexpr std::array<int, 4> speciesIds{kPiPlus, kProton, kElectron, kKPlus};

struct PureMcMultiplicityCounter {
  SliceCache cache;
  Service<o2::framework::O2DatabasePDG> pdg;
  Configurable<float> etaRange{"eta-range", 1.0f, "Eta range to consider"};
  ConfigurableAxis multBinning{"multBinning", {301, -0.5, 300.5}, ""};

  HistogramRegistry registry{
    "registry",
    {
      {"Events/Vertex/X", "; X (cm); events", {HistType::kTH1F, {ZAxis}}},
      {"Events/Vertex/Y", "; Y (cm); events", {HistType::kTH1F, {ZAxis}}},
      {"Events/Vertex/Z", "; Z (cm); events", {HistType::kTH1F, {ZAxis}}},
      {"Particles/Primaries/Pt", " ;p_{T} (GeV/c)", {HistType::kTH1F, {PtAxis_wide}}},
      {"Particles/Primaries/Eta", " ; #eta", {HistType::kTH1F, {EtaAxis}}},
      {"Particles/Primaries/Y", " ; y", {HistType::kTH1F, {RapidityAxis}}},
    } //
  };

  void init(InitContext const&)
  {
    AxisSpec MultAxis = {multBinning, "N_{p}"};
    registry.add({"Events/Properties/Multiplicity", " ; N_{p}; events", {HistType::kTH1F, {MultAxis}}});
    for (auto i = 0u; i < speciesIds.size(); ++i) {
      registry.add({(std::string("Particles/Primaries/") + std::string(species[i]) + "/Pt").c_str(), " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxis_wide}}});
      registry.add({(std::string("Particles/Primaries/") + std::string(species[i]) + "/Eta").c_str(), " ; #eta", {HistType::kTH1F, {EtaAxis}}});
      registry.add({(std::string("Particles/Primaries/") + std::string(species[i]) + "/Y").c_str(), " ; y", {HistType::kTH1F, {RapidityAxis}}});
    }
  }

  void process(aod::McCollision const& collision, aod::McParticles const& particles)
  {
    registry.fill(HIST("Events/Vertex/X"), collision.posX());
    registry.fill(HIST("Events/Vertex/Y"), collision.posY());
    registry.fill(HIST("Events/Vertex/Z"), collision.posZ());

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
    registry.fill(HIST("Events/Properties/Multiplicity"), Np);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return {adaptAnalysisTask<PureMcMultiplicityCounter>(cfgc)};
}
