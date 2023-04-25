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
// ========================
//
// This code loops over v0 photons and makes pairs for neutral mesons analyses.
//    Please write to: daiki.sekihata@cern.ch

#include "TString.h"
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include "Math/LorentzRotation.h"
#include "Math/Rotation3D.h"
#include "Math/AxisAngle.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/RecoDecay.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/PCMUtilities.h"
#include "PWGEM/PhotonMeson/Utils/MCUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

using MyV0Photons = soa::Join<aod::V0Photons, aod::V0RecalculationAndKF>;
using MyV0Photon = MyV0Photons::iterator;

struct Pi0EtaToGammaGammaMC {
  enum PairType {
    kPCMPCM = 0,
    kPHOSPHOS = 1,
    kEMCEMC = 2,
    kPCMPHOS = 3,
    kPCMEMC = 4,
    kPHOSEMC = 5,
  };

  using MyMCV0Legs = soa::Join<aod::V0Legs, aod::EMMCParticleLabels>;

  HistogramRegistry registry{"Pi0EtaToGammaGammaMC"};

  void init(InitContext& context)
  {
    addhistograms();
  }

  static constexpr std::string_view pairnames[6] = {"PCMPCM", "PHOSPHOS", "EMCEMC", "PCMPHOS", "PCMEMC", "PHOSEMC"};
  static constexpr std::string_view parnames[3] = {"Gamma", "Pi0", "Eta"};
  void addhistograms()
  {
    for (int i = 0; i < 6; i++) {
      registry.add(Form("%s/hCollisionCounter", pairnames[i].data()), "Collision counter", HistType::kTH1F, {{5, 0.5f, 5.5f}});
      for (int j = 0; j < 3; j++) {
        registry.add(Form("%s/h2MggPt_%s", pairnames[i].data(), parnames[j].data()), "M_{#gamma#gamma} vs. p_{T};m_{#gamma#gamma} (GeV/c^{2});p_{T,#gamma#gamma} (GeV/c)", HistType::kTH2F, {{400, 0, 0.8}, {400, 0.0f, 40}}, true);
      }
    }

    registry.add("Generated/hCollisionCounter", "Collision counter", HistType::kTH1F, {{5, 0.5f, 5.5f}});
    registry.add("Generated/hRecCollision", "reconstructed collision per MC collision;N_{collision}^{rec.} per MC collision;counts", HistType::kTH1F, {{101, -0.5f, 100.5f}});
    for (int j = 0; j < 3; j++) {
      registry.add(Form("Generated/hPt_%s", parnames[j].data()), "generated p_{T};p_{T} (GeV/c)", HistType::kTH1F, {{400, 0, 40}}, true);
      registry.add(Form("Generated/hY_%s", parnames[j].data()), "generated rapidity; rapidity y", HistType::kTH1F, {{40, -2, +2}}, true);
      registry.add(Form("Generated/hPhi_%s", parnames[j].data()), "generated azimuthal angle;#varphi (rad.)", HistType::kTH1F, {{90, 0, TMath::TwoPi()}}, true);
    }
  }

  Preslice<MyV0Photons> perCollision_pcm = aod::v0photon::collisionId;

  template <PairType pairtype, typename TEvents, typename TPhotons1, typename TPhotons2, typename TPreslice1, typename TPreslice2, typename TV0Legs, typename TMCParticles>
  void TruePairing(TEvents const& collisions, TPhotons1 const& photons1, TPhotons2 const& photons2, TPreslice1 const& perCollision1, TPreslice2 const& perCollision2, TV0Legs const& v0legs, TMCParticles const& mcparticles)
  {
    constexpr int itmp = pairtype;
    for (auto& collision : collisions) {

      if ((pairtype == kPHOSPHOS || pairtype == kPCMPHOS) && !collision.isPHOSCPVreadout()) {
        continue;
      }
      if ((pairtype == kEMCEMC || pairtype == kPCMEMC) && !collision.isEMCreadout()) {
        continue;
      }

      registry.fill(HIST(pairnames[itmp]) + HIST("/hCollisionCounter"), 1.0); // all
      if (!collision.sel8()) {
        continue;
      }
      registry.fill(HIST(pairnames[itmp]) + HIST("/hCollisionCounter"), 2.0); // FT0VX i.e. FT0and

      if (collision.numContrib() < 0.5) {
        continue;
      }
      registry.fill(HIST(pairnames[itmp]) + HIST("/hCollisionCounter"), 3.0); // Ncontrib > 0

      if (abs(collision.posZ()) > 10.0) {
        continue;
      }
      registry.fill(HIST(pairnames[itmp]) + HIST("/hCollisionCounter"), 4.0); //|Zvtx| < 10 cm

      auto photons1_coll = photons1.sliceBy(perCollision1, collision.collisionId());
      auto photons2_coll = photons2.sliceBy(perCollision2, collision.collisionId());

      int pi0id = -1;
      int etaid = -1;
      if (pairtype == PairType::kPCMPCM || pairtype == PairType::kPHOSPHOS || pairtype == PairType::kEMCEMC) {
        for (auto& [g1, g2] : combinations(CombinationsStrictlyUpperIndexPolicy(photons1_coll, photons2_coll))) {

          if constexpr (pairtype == PairType::kPCMPCM) { // check 2 legs
            auto pos1 = g1.template posTrack_as<MyMCV0Legs>();
            auto ele1 = g1.template negTrack_as<MyMCV0Legs>();
            auto pos2 = g2.template posTrack_as<MyMCV0Legs>();
            auto ele2 = g2.template negTrack_as<MyMCV0Legs>();

            auto pos1mc = pos1.template emmcparticle_as<aod::EMMCParticles>();
            auto ele1mc = ele1.template emmcparticle_as<aod::EMMCParticles>();
            auto pos2mc = pos2.template emmcparticle_as<aod::EMMCParticles>();
            auto ele2mc = ele2.template emmcparticle_as<aod::EMMCParticles>();
            // LOGF(info,"pos1mc.globalIndex() = %d , ele1mc.globalIndex() = %d , pos2mc.globalIndex() = %d , ele2mc.globalIndex() = %d", pos1mc.globalIndex(), ele1mc.globalIndex(), pos2mc.globalIndex(), ele2mc.globalIndex());

            int photonid1 = FindCommonMotherFrom2Prongs(pos1mc, ele1mc, -11, 11, 22, mcparticles);
            int photonid2 = FindCommonMotherFrom2Prongs(pos2mc, ele2mc, -11, 11, 22, mcparticles);

            if (photonid1 < 0) { // check swap, true electron is reconstructed as positron and vice versa.
              photonid1 = FindCommonMotherFrom2Prongs(pos1mc, ele1mc, 11, -11, 22, mcparticles);
            }
            if (photonid2 < 0) { // check swap, true electron is reconstructed as positron and vice versa.
              photonid2 = FindCommonMotherFrom2Prongs(pos2mc, ele2mc, 11, -11, 22, mcparticles);
            }

            // LOGF(info,"photonid1 = %d , photonid2 = %d", photonid1, photonid2);
            if (photonid1 < 0 || photonid2 < 0) {
              continue;
            }

            // if(photonid1 == photonid2) {
            //   continue;
            // }
            auto g1mc = mcparticles.iteratorAt(photonid1);
            auto g2mc = mcparticles.iteratorAt(photonid2);
            pi0id = FindCommonMotherFrom2Prongs(g1mc, g2mc, 22, 22, 111, mcparticles);
            etaid = FindCommonMotherFrom2Prongs(g1mc, g2mc, 22, 22, 221, mcparticles);
          }

          ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
          ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

          if (pi0id > 0) {
            registry.fill(HIST(pairnames[itmp]) + HIST("/h2MggPt_Pi0"), v12.M(), v12.Pt());
          }
          if (etaid > 0) {
            registry.fill(HIST(pairnames[itmp]) + HIST("/h2MggPt_Eta"), v12.M(), v12.Pt());
          }
        } // end of combination

      } else {
        for (auto& [g1, g2] : combinations(CombinationsFullIndexPolicy(photons1_coll, photons2_coll))) {

          ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
          ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          registry.fill(HIST(pairnames[itmp]) + HIST("/h2MggPt_Same"), v12.M(), v12.Pt());

        } // end of combination
      }

    } // end of collision loop
  }

  // using MyCollisions = soa::Join<aod::EMReducedEvents, aod::EMReducedMCEventLabels>;

  void processPCMPCM(soa::Join<aod::EMReducedEvents, aod::EMReducedMCEventLabels> const& collisions, MyV0Photons const& v0photons, MyMCV0Legs const& v0legs, aod::EMMCParticles const& mcparticles)
  {
    TruePairing<PairType::kPCMPCM>(collisions, v0photons, v0photons, perCollision_pcm, perCollision_pcm, v0legs, mcparticles);
  }

  Configurable<float> maxYgen{"maxYgen", 0.9, "maximum rapidity for generated particles"};

  // Preslice<aod::EMMCParticles> perMcCollision = aod::emmcparticle::emreducedmceventId;
  // Preslice<soa::Join<aod::EMReducedEvents, aod::EMReducedMCEventLabels>> rec_perMcCollision = aod::emmceventlabel::emreducedmceventId;
  void processGen(soa::Join<aod::EMReducedEvents, aod::EMReducedMCEventLabels> const& collisions, aod::EMReducedMCEvents const& mccollisions, aod::EMMCParticles const& mcparticles)
  {
    // loop over mc stack and fill histograms for pure MC truth signals
    // all MC tracks which belong to the MC event corresponding to the current reconstructed event

    for (auto& mccollision : mccollisions) {
      int nrec = 0;
      // auto collision_per_mccoll = collisions.sliceBy(rec_perMcCollision, mccollision.globalIndex());
      // int nrec_per_mc = collision_per_mccoll.size();
      for (auto& collision : collisions) {
        if (mccollision.globalIndex() == collision.emreducedmceventId()) {
          nrec++;
        }
      }
      registry.fill(HIST("Generated/hRecCollision"), nrec);
    }

    for (auto& collision : collisions) {
      registry.fill(HIST("Generated/hCollisionCounter"), 1.0); // all
      if (!collision.sel8()) {
        continue;
      }
      registry.fill(HIST("Generated/hCollisionCounter"), 2.0); // FT0VX i.e. FT0and

      if (collision.numContrib() < 0.5) {
        continue;
      }
      registry.fill(HIST("Generated/hCollisionCounter"), 3.0); // Ncontrib > 0

      if (abs(collision.posZ()) > 10.0) {
        continue;
      }
      registry.fill(HIST("Generated/hCollisionCounter"), 4.0); //|Zvtx| < 10 cm
      auto mccollision = collision.emreducedmcevent();

      // auto mctracks_coll = mcparticles.sliceBy(perMcCollision, mccollision.globalIndex());
      // for (auto& mctrack : mctracks_coll) {
      for (auto& mctrack : mcparticles) {

        if (mctrack.emreducedmceventId() != mccollision.globalIndex()) {
          continue;
        }

        if (abs(mctrack.y()) > maxYgen) {
          continue;
        }
        int pdg = mctrack.pdgCode();

        if (abs(pdg) == 22 && IsPhysicalPrimary(mctrack.emreducedmcevent(), mctrack, mcparticles)) {
          registry.fill(HIST("Generated/hPt_Gamma"), mctrack.pt());
          registry.fill(HIST("Generated/hY_Gamma"), mctrack.y());
          registry.fill(HIST("Generated/hPhi_Gamma"), mctrack.phi());
        }
        if (abs(pdg) == 111 && IsPhysicalPrimary(mctrack.emreducedmcevent(), mctrack, mcparticles)) {
          registry.fill(HIST("Generated/hPt_Pi0"), mctrack.pt());
          registry.fill(HIST("Generated/hY_Pi0"), mctrack.y());
          registry.fill(HIST("Generated/hPhi_Pi0"), mctrack.phi());
        }
        if (abs(pdg) == 221 && IsPhysicalPrimary(mctrack.emreducedmcevent(), mctrack, mcparticles)) {
          registry.fill(HIST("Generated/hPt_Eta"), mctrack.pt());
          registry.fill(HIST("Generated/hY_Eta"), mctrack.y());
          registry.fill(HIST("Generated/hPhi_Eta"), mctrack.phi());
        }

      } // end of mc track loop
    }   // end of collision loop
  }

  void processDummy(aod::EMReducedEvents::iterator const& collision) {}

  PROCESS_SWITCH(Pi0EtaToGammaGammaMC, processPCMPCM, "true pairing PCM-PCM", false);
  PROCESS_SWITCH(Pi0EtaToGammaGammaMC, processGen, "process generated information", false);
  PROCESS_SWITCH(Pi0EtaToGammaGammaMC, processDummy, "Dummy function", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<Pi0EtaToGammaGammaMC>(cfgc, TaskName{"pi0eta-to-gammagamma-mc"})};
}
