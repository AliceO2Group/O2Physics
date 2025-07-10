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
//
/// \file testMcStdTabsRl.cxx
/// \brief task to test the Monte Carlo UD production generatorIDs on hyperloop
///
/// \author Roman Lavicka <roman.lavicka@cern.ch>, Austrian Academy of Sciences & SMI
/// \since  12.02.2025
//

// C++ headers
#include <algorithm>
#include <set>
#include <utility>
#include <vector>

// O2 headers
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

// O2Physics headers
#include "PWGUD/Core/UPCTauCentralBarrelHelperRL.h"

// ROOT headers
#include "Math/Vector4D.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

struct TestMcStdTabsRl {

  // Global varialbes
  Service<o2::framework::O2DatabasePDG> pdg;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // declare configurables
  Configurable<bool> useMcgenidGetGeneratorID{"useMcgenidGetGeneratorID", true, {"Use o2::mcgenid::getGeneratorId instead of o2::mccollision::getGeneratorId; default it true."}};

  struct : ConfigurableGroup {
    ConfigurableAxis zzAxisNtracks{"zzAxisNtracks", {100, -0.5, 99.5}, "Number of tracks in collision"};
    ConfigurableAxis zzAxisNparticles{"zzAxisNparticles", {100, -0.5, 99.5}, "Number of particles in collision"};
    ConfigurableAxis zzAxisNprocesses{"zzAxisNprocesses", {1000, -0.5, 999.5}, "Number of processes"};
    ConfigurableAxis zzAxisInvMassWide{"zzAxisInvMassWide", {1000, 0., 10.}, "Invariant mass (GeV/c^{2}), wider range"};
    ConfigurableAxis zzAxisPt{"zzAxisPt", {400, 0., 2.}, "Transversal momentum (GeV/c)"};
    ConfigurableAxis zzAxisRap{"zzAxisRap", {50, -1.2, 1.2}, "Rapidity (a.u.)"};
  } confAxis;

  // init
  void init(InitContext&)
  {
    histos.add("Events/Truth/hGenIDvsCountCollisions", ";Process ID", HistType::kTH2D, {confAxis.zzAxisNprocesses, {1, 0.5, 1.5}});
    histos.add("Events/Truth/hGenIDvsPDGcodesAll", ";Process ID ;PDG codes of all particles (-)", HistType::kTH2D, {confAxis.zzAxisNprocesses, {2001, -1000, 1000}});
    histos.add("Events/Truth/hGenIDvsPDGcodesNoMother", ";Process ID ;PDG codes of particles without mother (-)", HistType::kTH2D, {confAxis.zzAxisNprocesses, {2001, -1000, 1000}});
    histos.add("Events/Truth/hGenIDvsPDGcodesDaughters", ";Process ID ;PDG codes of daughters of particles without mother (-)", HistType::kTH2D, {confAxis.zzAxisNprocesses, {2001, -1000, 1000}});
    histos.add("Events/Truth/hGenIDvsNparticles", ";Process ID ;Number of particles in a collision (-)", HistType::kTH2D, {confAxis.zzAxisNprocesses, confAxis.zzAxisNparticles});
    histos.add("Events/Truth/hGenIDvsNdaughters", ";Process ID ;Number of daughters of no-mother particle in a collision (-)", HistType::kTH2D, {confAxis.zzAxisNprocesses, confAxis.zzAxisNparticles});
    histos.add("Events/Truth/hGenIDvsMotherMass", ";Process ID ;Mother invariant mass (GeV/c^{2})", HistType::kTH2D, {confAxis.zzAxisNprocesses, confAxis.zzAxisInvMassWide});
    histos.add("Events/Truth/hGenIDvsMotherPt", ";Process ID ;Mother p_{T} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisNprocesses, confAxis.zzAxisPt});
    histos.add("Events/Truth/hGenIDvsMotherRap", ";Process ID ;Mother rapidity (-)", HistType::kTH2D, {confAxis.zzAxisNprocesses, confAxis.zzAxisRap});

  } // end init

  void processMCgen(aod::McCollision const& collision, aod::McParticles const& particles)
  {

    const auto genID = useMcgenidGetGeneratorID ? o2::mcgenid::getGeneratorId(collision.getGeneratorId()) : collision.getGeneratorId();

    histos.get<TH2>(HIST("Events/Truth/hGenIDvsCountCollisions"))->Fill(genID, 1);
    histos.get<TH2>(HIST("Events/Truth/hGenIDvsNparticles"))->Fill(genID, particles.size());

    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> mother;
    for (const auto& particle : particles) {
      histos.get<TH2>(HIST("Events/Truth/hGenIDvsPDGcodesAll"))->Fill(genID, particle.pdgCode());
      //        if (!particle.isPhysicalPrimary()) continue;
      if (particle.has_mothers())
        continue;
      mother.SetPxPyPzE(particle.px(), particle.py(), particle.pz(), energy(pdg->Mass(particle.pdgCode()), particle.px(), particle.py(), particle.pz()));
      histos.get<TH2>(HIST("Events/Truth/hGenIDvsPDGcodesNoMother"))->Fill(genID, particle.pdgCode());
      histos.get<TH2>(HIST("Events/Truth/hGenIDvsMotherMass"))->Fill(genID, mother.M());
      histos.get<TH2>(HIST("Events/Truth/hGenIDvsMotherPt"))->Fill(genID, particle.pt());
      histos.get<TH2>(HIST("Events/Truth/hGenIDvsMotherRap"))->Fill(genID, particle.y());
      const auto& daughters = particle.daughters_as<aod::McParticles>();
      histos.get<TH2>(HIST("Events/Truth/hGenIDvsNdaughters"))->Fill(genID, daughters.size());
      for (const auto& daughter : daughters) {
        histos.get<TH2>(HIST("Events/Truth/hGenIDvsPDGcodesDaughters"))->Fill(genID, daughter.pdgCode());
      }
    }

  } // end processMCgenDG

  PROCESS_SWITCH(TestMcStdTabsRl, processMCgen, "Iterate Monte Carlo UD tables with truth data.", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TestMcStdTabsRl>(cfgc)};
}
