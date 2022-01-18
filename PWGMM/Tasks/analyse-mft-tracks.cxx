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

using namespace o2;
using namespace o2::framework;

using Particles = aod::McParticles;

//the task analyseMFTTracks loops over MFT tracks and generated particles and fills basic histograms

struct analyseMFTTracks {
  int icoll = 0;
  Service<TDatabasePDG> pdg;
  HistogramRegistry registry{
    "registry",
    {
      {"TracksPhiEta_in_coll", "; #varphi; #eta; ", {HistType::kTH2F, {{600, -M_PI, M_PI}, {35, -4.5, -1.}}}}, //
      {"TracksPhiEta", "; #varphi; #eta; ", {HistType::kTH2F, {{600, -M_PI, M_PI}, {35, -4.5, -1.}}}}, //

      {"TracksPhiEta_5", "#eta_{kin} versus #varphi_{kin} of MFT tracks with ncls #\geq 5; #varphi; #eta; ", {HistType::kTH2F, {{600, -M_PI, M_PI}, {35, -4.5, -1.}}}}, //
      {"TracksPhiEta_6", "#eta_{kin} versus #varphi_{kin} of MFT tracks with ncls #\geq 6; #varphi; #eta; ", {HistType::kTH2F, {{600, -M_PI, M_PI}, {35, -4.5, -1.}}}}, //
      {"TracksPhiEta_7", "#eta_{kin} versus #varphi_{kin} of MFT tracks with ncls #\geq 7; #varphi; #eta; ", {HistType::kTH2F, {{600, -M_PI, M_PI}, {35, -4.5, -1.}}}}, //
      {"TracksPhiEta_8", "#eta_{kin} versus #varphi_{kin} of MFT tracks with ncls #\geq 8; #varphi; #eta; ", {HistType::kTH2F, {{600, -M_PI, M_PI}, {35, -4.5, -1.}}}}, //
      {"TracksPhiEta_9", "#eta_{kin} versus #varphi_{kin} of MFT tracks with ncls #\geq 9; #varphi; #eta; ", {HistType::kTH2F, {{600, -M_PI, M_PI}, {35, -4.5, -1.}}}}, //

      {"TracksEtaZvtx", "; #eta; Z_{vtx}; ", {HistType::kTH2F, {{35, -4.5, -1.}, {201, -20.1, 20.1}}}},        //
      {"NtrkZvtx", "; N_{trk}; Z_{vtx}; ", {HistType::kTH2F, {{301, -0.5, 300.5}, {201, -20.1, 20.1}}}},       //
      {"NtrkEta", "#eta; N_{trk}; ", {HistType::kTH1F, {{35, -4.5, -1.}}}},                                    //
      {"TracksXLastClsYLastCls", "; x_{LastCls}; y_{LastCls}; ", {HistType::kTH2F, {{400, -20, 20.}, {400, -20, 20.}}}},        //
      {"TracksXLastClsYLastCls_firstFace", "; x_{LastCls}; y_{LastCls}; ", {HistType::kTH2F, {{400, -20, 20.}, {400, -20, 20.}}}},        //
      {"TracksXLastClsYLastCls_secondFace", "; x_{LastCls}; y_{LastCls}; ", {HistType::kTH2F, {{400, -20, 20.}, {400, -20, 20.}}}},
      {"TracksPhiEta_LastCls", "; #varphi; #eta; ", {HistType::kTH2F, {{600, -M_PI, M_PI}, {35, -4.5, -1.}}}}, //
      {"TracksPhiEta_LastClsLastDisk", "; #varphi; #eta; ", {HistType::kTH2F, {{600, -M_PI, M_PI}, {35, -4.5, -1.}}}}, //
      {"TracksPhiEta_LastClsLastDisk_firstFace", "; #varphi; #eta; ", {HistType::kTH2F, {{600, -M_PI, M_PI}, {35, -4.5, -1.}}}}, //
      {"TracksPhiEta_LastClsLastDisk_secondFace", "; #varphi; #eta; ", {HistType::kTH2F, {{600, -M_PI, M_PI}, {35, -4.5, -1.}}}}, //



      {"TracksPhiEtaGen", "; #varphi; #eta; ", {HistType::kTH2F, {{600, 0, 2 * M_PI}, {35, -4.5, -1.}}}},      //
      {"TracksEtaZvtxGen", "; #eta; Z_{vtx}; ", {HistType::kTH2F, {{35, -4.5, -1.}, {201, -20.1, 20.1}}}},     //
      {"NtrkZLastCls", "z_{LastCls}; count; ", {HistType::kTH1F, {{26, -79, -66}}}},
      //{"NtrkZvtxGen", "; N_{trk}; Z_{vtx}; events", {HistType::kTH2F, {{301, -0.5, 300.5}, {201, -20.1, 20.1}}}},            //
      {"NtrkEtaGen", "#eta; N_{trk}; ", {HistType::kTH1F, {{35, -4.5, -1.}}}}, //
    }                                                                                //
  };

  void processRec(o2::aod::Collision const& collision, o2::aod::MFTTracks const& tracks)
  {

    auto z = collision.posZ();
    registry.fill(HIST("NtrkZvtx"), tracks.size(), z);
    TVector3 v;

    for (auto& track : tracks) {
      double zLastCls = track.z();
      v.SetXYZ(track.x(),track.y(),zLastCls);
      double etaLastCls = v.Eta();
      double phiLastCls = v.Phi();
      registry.fill(HIST("TracksPhiEta_in_coll"), track.phi(), track.eta());
      registry.fill(HIST("TracksEtaZvtx"), track.eta(), z);
      registry.fill(HIST("NtrkEta"), track.eta());
      if (zLastCls < -70)//disk 4
      {
        registry.fill(HIST("TracksXLastClsYLastCls"), track.x(), track.y());
        registry.fill(HIST("TracksPhiEta_LastClsLastDisk"), phiLastCls, etaLastCls);

        if (zLastCls > -77)//layer 8
        {
          registry.fill(HIST("TracksXLastClsYLastCls_firstFace"), track.x(), track.y());
          registry.fill(HIST("TracksPhiEta_LastClsLastDisk_firstFace"), phiLastCls, etaLastCls);
        }

      }
      if (zLastCls < -77)//layer 9
      {
        registry.fill(HIST("TracksXLastClsYLastCls_secondFace"), track.x(), track.y());
        registry.fill(HIST("TracksPhiEta_LastClsLastDisk_secondFace"), phiLastCls, etaLastCls);
      }
      registry.fill(HIST("TracksPhiEta_LastCls"), phiLastCls, etaLastCls);
      registry.fill(HIST("NtrkZLastCls"), track.z());
    }
    icoll++;
  }
  //end of processRec
  PROCESS_SWITCH(analyseMFTTracks, processRec, "Process rec level", true);



  void processTracks(o2::aod::MFTTracks const& tracks)
  {

    for (auto& track : tracks)
    {
      registry.fill(HIST("TracksPhiEta"), track.phi(), track.eta());

      if (track.nClusters()>=5)
      {
        registry.fill(HIST("TracksPhiEta_5"), track.phi(), track.eta());
      }
      if (track.nClusters()>=6)
      {
        registry.fill(HIST("TracksPhiEta_6"), track.phi(), track.eta());
      }
      if (track.nClusters()>=7)
      {
        registry.fill(HIST("TracksPhiEta_7"), track.phi(), track.eta());
      }
      if (track.nClusters()>=8)
      {
        registry.fill(HIST("TracksPhiEta_8"), track.phi(), track.eta());
      }
      if (track.nClusters()>=9)
      {
        registry.fill(HIST("TracksPhiEta_9"), track.phi(), track.eta());
      }
    }//end of tracks loop
  }
  //end of processTracks
  PROCESS_SWITCH(analyseMFTTracks, processTracks, "Process rec level for all tracks", true);





  void processGen(aod::McCollisions::iterator const& mcCollision, Particles const& particles)
  {

    int nChargedPrimaryParticles = 0;
    auto z = mcCollision.posZ();

    for (auto& particle : particles) {
      auto p = pdg->GetParticle(particle.pdgCode());
      //auto p = TDatabasePDG::Instance()->GetParticle(particle.pdgCode());
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
        registry.fill(HIST("TracksEtaZvtxGen"), particle.eta(), z);
        registry.fill(HIST("TracksPhiEtaGen"), particle.phi(), particle.eta());

        registry.fill(HIST("TracksPhiEtaGen"), particle.phi(), particle.eta());
        registry.fill(HIST("NtrkEtaGen"), particle.eta());
        nChargedPrimaryParticles++;
      }
    }

    //registry.fill(HIST("NtrkZvtxGen"), nChargedPrimaryParticles, mcCollision.posZ());
  }

  PROCESS_SWITCH(analyseMFTTracks, processGen, "Process gen level", false);
};
//end of the task analyseMFTTracks

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<analyseMFTTracks>(cfgc),
  };
}
