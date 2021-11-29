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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

#include "TDatabasePDG.h"
#include "Common/Core/MC.h"

using namespace o2;
using namespace o2::framework;

using Particles = aod::McParticles;

//First approach to analysing an AO2D.root file
//written thanks to https://aliceo2group.github.io/analysis-framework/docs/tutorials/analysistask.html

struct analyseMFTTracks {
  int icoll = 0;

  HistogramRegistry registry{
    "registry",
    {
      {"TracksPhiEta_in_coll", "; #varphi; #eta; tracks", {HistType::kTH2F, {{600, -M_PI, M_PI}, {35, -4.5, -1.}}}},
      {"TracksPhiEta", "; #varphi; #eta; tracks", {HistType::kTH2F, {{600, -M_PI, M_PI}, {35, -4.5, -1.}}}},                                                           //
      {"TracksPhiEtaGen", "; #varphi; #eta; tracks", {HistType::kTH2F, {{600, 0, 2 * M_PI}, {35, -4.5, -1.}}}},                                                        //
      {"TracksEtaZvtx", "; #eta; Z_{vtx}; tracks", {HistType::kTH2F, {{35, -4.5, -1.}, {201, -20.1, 20.1}}}},                                                          //
      {"NtrkZvtx", "; N_{trk}; Z_{vtx}; events", {HistType::kTH2F, {{301, -0.5, 300.5}, {201, -20.1, 20.1}}}},                                                         //
      {"NtrkEta", "; N_{trk}; #eta; events", {HistType::kTH2F, {{301, -0.5, 300.5}, {35, -4.5, -1.}}}},                                                                //
      {"Multiplicity", "alibi_nightlies/O2DPG_pp_minbias_testbeam.sh/15-11-2021-18:00 - tf13; collisionID; N_{trk}^{MFT}", {HistType::kTH1F, {{101, 0., 100.}}}, true} //
    }                                                                                                                                                                  //
  };

  void process(o2::aod::Collision const& collision, o2::aod::MFTTracks const& tracks)
  {

    auto z = collision.posZ();
    registry.fill(HIST("NtrkZvtx"), tracks.size(), z);

    for (auto& track : tracks) {
      registry.fill(HIST("TracksPhiEta_in_coll"), track.phi(), track.eta());
      registry.fill(HIST("TracksEtaZvtx"), track.eta(), z);
      registry.fill(HIST("Multiplicity"), icoll);
      registry.fill(HIST("NtrkEta"), tracks.size(), track.eta());
    }
    icoll++;
  }
  //end of process

  //aod::McCollisions
};
//end of MyTask

struct analyseGenTracks {

  HistogramRegistry registryGen{
    "registryGen",
    {
      {"TracksPhiEtaGen", "; #varphi; #eta; tracks", {HistType::kTH2F, {{600, -M_PI, M_PI}, {35, -4.5, -1.}}}},   //
      {"TracksEtaZvtxGen", "; #eta; Z_{vtx}; tracks", {HistType::kTH2F, {{35, -4.5, -1.}, {201, -20.1, 20.1}}}},  //
      {"NtrkZvtxGen", "; N_{trk}; Z_{vtx}; events", {HistType::kTH2F, {{301, -0.5, 300.5}, {201, -20.1, 20.1}}}}, //
      {"NtrkEtaGen", "; N_{trk}; #eta; events", {HistType::kTH2F, {{301, -0.5, 300.5}, {35, -4.5, -1.}}}},        //
    }                                                                                                             //
  };

  expressions::Filter posZFilterMC = (aod::mccollision::posZ < 15) && (aod::mccollision::posZ > -15);

  void process(soa::Filtered<aod::McCollisions>::iterator const& mcCollision, Particles const& particles)
  {
    int nChargedPrimaryParticles = 0;
    auto z = mcCollision.posZ();

    for (auto& particle : particles) {
      auto p = TDatabasePDG::Instance()->GetParticle(particle.pdgCode());
      int charge = 0;
      if (p == nullptr) {
        // unknown particles will be skipped
        if (particle.pdgCode() > 1000000000) {
          //          auto x = (std::trunc(particle.pdgCode() / 10000) - 100000);
          //          charge = x - std::trunc(x / 1000) * 1000;
          LOGF(debug, "[{}] Nucleus with PDG code {}", particle.globalIndex(), particle.pdgCode() /*, charge*/); // (charge %d)
        } else {
          LOGF(debug, "[{}] Unknown particle with PDG code {}", particle.globalIndex(), particle.pdgCode());
        }
      } else {
        charge = p->Charge();
      }
      if (charge != 0 && particle.isPhysicalPrimary()) {
        registryGen.fill(HIST("TracksEtaZvtxGen"), particle.eta(), z);
        registryGen.fill(HIST("TracksPhiEtaGen"), particle.phi(), particle.eta());

        registryGen.fill(HIST("TracksPhiEtaGen"), particle.phi(), particle.eta());
        nChargedPrimaryParticles++;
      }
    }

    //registryGen.fill(HIST("NtrkZvtxGen"), nChargedPrimaryParticles, mcCollision.posZ());
  }
};
//end of the gen task

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<analyseMFTTracks>(cfgc),
    //adaptAnalysisTask<analyseGenTracks>(cfgc),
  };
}
