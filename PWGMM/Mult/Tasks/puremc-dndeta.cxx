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

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

AxisSpec ZAxis = {301, -30.1, 30.1};
AxisSpec DeltaZAxis = {61, -6.1, 6.1};
AxisSpec DCAAxis = {601, -3.01, 3.01};
AxisSpec EtaAxis = {22, -2.2, 2.2};
AxisSpec PhiAxis = {629, 0, 2 * M_PI};
AxisSpec PtAxis = {2401, -0.005, 24.005};
AxisSpec PtAxis_wide = {1041, -0.05, 104.05};

auto static constexpr mincharge = 3.f;

struct PureMcMultiplicityCounter {
  SliceCache cache;
  Service<o2::framework::O2DatabasePDG> pdg;
  HistogramRegistry histos{
    "histos",
    {
      {"Events/Vertex/X", "; X (cm); events", {HistType::kTH1F, {ZAxis}}},
      {"Events/Vertex/Y", "; Y (cm); events", {HistType::kTH1F, {ZAxis}}},
      {"Events/Vertex/Z", "; Z (cm); events", {HistType::kTH1F, {ZAxis}}},
      {"Particles/Primaries/Pt", " ;p_{T} (GeV/c)", {HistType::kTH1F, {PtAxis_wide}}},
    }};
  Configurable<float> etaRange{"eta-range", 1.0f, "Eta range to consider"};

  void Init(InitContext const&) {}
  void process(aod::McCollision const& collision, aod::McParticles const& particles)
  {
    histos.fill(HIST("Events/Vertex/X"), collision.posX());
    histos.fill(HIST("Events/Vertex/Y"), collision.posY());
    histos.fill(HIST("Events/Vertex/Z"), collision.posZ());

    for (auto const& particle : particles) {
      if (!particle.isPhysicalPrimary()) {
        continue;
      }
      if (std::abs(particle.eta()) >= etaRange) {
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

      histos.fill(HIST("Particles/Primaries/Pt"), particle.pt());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return {adaptAnalysisTask<PureMcMultiplicityCounter>(cfgc)};
}
