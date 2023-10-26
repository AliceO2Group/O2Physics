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
// This code runs loop over v0 photons for PCM QC.
//    Please write to: daiki.sekihata@cern.ch

#include <array>
#include "TString.h"
#include "THashList.h"
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
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/Core/CutsLibrary.h"
#include "PWGEM/PhotonMeson/Core/HistogramsLibrary.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using std::array;

using MyV0Photons = soa::Join<aod::V0PhotonsKF, aod::V0Recalculation, aod::V0KFEMReducedEventIds>;
using MyV0Photon = MyV0Photons::iterator;

struct PCMQCMC {
  Configurable<std::string> fConfigPCMCuts{"cfgPCMCuts", "analysis,qc,nocut", "Comma separated list of v0 photon cuts"};
  Configurable<float> maxY{"maxY", 0.9, "maximum rapidity for generated particles"};

  std::vector<V0PhotonCut> fPCMCuts;

  OutputObj<THashList> fOutputEvent{"Event"};
  OutputObj<THashList> fOutputV0Leg{"V0Leg"};
  OutputObj<THashList> fOutputV0{"V0"};
  OutputObj<THashList> fOutputGen{"Generated"};
  THashList* fMainList = new THashList();

  void addhistograms()
  {
    fMainList->SetOwner(true);
    fMainList->SetName("fMainList");

    // create sub lists first.
    o2::aod::emphotonhistograms::AddHistClass(fMainList, "Event");
    THashList* list_ev = reinterpret_cast<THashList*>(fMainList->FindObject("Event"));
    o2::aod::emphotonhistograms::DefineHistograms(list_ev, "Event");

    o2::aod::emphotonhistograms::AddHistClass(fMainList, "V0Leg");
    THashList* list_v0leg = reinterpret_cast<THashList*>(fMainList->FindObject("V0Leg"));

    o2::aod::emphotonhistograms::AddHistClass(fMainList, "V0");
    THashList* list_v0 = reinterpret_cast<THashList*>(fMainList->FindObject("V0"));

    o2::aod::emphotonhistograms::AddHistClass(fMainList, "Generated");
    THashList* list_gen = reinterpret_cast<THashList*>(fMainList->FindObject("Generated"));
    o2::aod::emphotonhistograms::DefineHistograms(list_gen, "Generated", "Photon");

    for (const auto& cut : fPCMCuts) {
      const char* cutname = cut.GetName();
      o2::aod::emphotonhistograms::AddHistClass(list_v0leg, cutname);
      o2::aod::emphotonhistograms::AddHistClass(list_v0, cutname);
    }

    // for single tracks
    for (auto& cut : fPCMCuts) {
      std::string_view cutname = cut.GetName();
      THashList* list = reinterpret_cast<THashList*>(fMainList->FindObject("V0Leg")->FindObject(cutname.data()));
      o2::aod::emphotonhistograms::DefineHistograms(list, "V0Leg");
    }

    // for V0s
    for (auto& cut : fPCMCuts) {
      std::string_view cutname = cut.GetName();
      THashList* list = reinterpret_cast<THashList*>(fMainList->FindObject("V0")->FindObject(cutname.data()));
      o2::aod::emphotonhistograms::DefineHistograms(list, "V0", "mc");
    }
  }

  void DefineCuts()
  {
    TString cutNamesStr = fConfigPCMCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        const char* cutname = objArray->At(icut)->GetName();
        LOGF(info, "add cut : %s", cutname);
        fPCMCuts.push_back(*pcmcuts::GetCut(cutname));
      }
    }
    LOGF(info, "Number of PCM cuts = %d", fPCMCuts.size());
  }

  void init(InitContext& context)
  {
    DefineCuts();
    addhistograms(); // please call this after DefinCuts();

    fOutputEvent.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("Event")));
    fOutputV0Leg.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("V0Leg")));
    fOutputV0.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("V0")));
    fOutputGen.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("Generated")));
  }

  Preslice<MyV0Photons> perCollision = aod::v0photon::emreducedeventId;
  using MyMCV0Legs = soa::Join<aod::V0Legs, aod::V0LegMCLabels>;
  void processQCMC(soa::Join<aod::EMReducedEvents, aod::EMReducedMCEventLabels> const& collisions, MyV0Photons const& v0photons, MyMCV0Legs const& v0legs, aod::EMMCParticles const& mcparticles, aod::EMReducedMCEvents const&)
  {
    THashList* list_ev = static_cast<THashList*>(fMainList->FindObject("Event"));
    THashList* list_v0 = static_cast<THashList*>(fMainList->FindObject("V0"));
    THashList* list_v0leg = static_cast<THashList*>(fMainList->FindObject("V0Leg"));

    for (auto& collision : collisions) {
      reinterpret_cast<TH1F*>(fMainList->FindObject("Event")->FindObject("hZvtx_before"))->Fill(collision.posZ());
      reinterpret_cast<TH1F*>(fMainList->FindObject("Event")->FindObject("hCollisionCounter"))->Fill(1.0);
      if (!collision.sel8()) {
        continue;
      }
      reinterpret_cast<TH1F*>(fMainList->FindObject("Event")->FindObject("hCollisionCounter"))->Fill(2.0);

      if (collision.numContrib() < 0.5) {
        continue;
      }
      reinterpret_cast<TH1F*>(fMainList->FindObject("Event")->FindObject("hCollisionCounter"))->Fill(3.0);

      if (abs(collision.posZ()) > 10.0) {
        continue;
      }
      reinterpret_cast<TH1F*>(fMainList->FindObject("Event")->FindObject("hCollisionCounter"))->Fill(4.0);
      reinterpret_cast<TH1F*>(fMainList->FindObject("Event")->FindObject("hZvtx_after"))->Fill(collision.posZ());
      o2::aod::emphotonhistograms::FillHistClass<EMHistType::kEvent>(list_ev, "", collision);
      auto V0Photons_coll = v0photons.sliceBy(perCollision, collision.collisionId());

      for (const auto& cut : fPCMCuts) {
        THashList* list_v0_cut = static_cast<THashList*>(list_v0->FindObject(cut.GetName()));
        THashList* list_v0leg_cut = static_cast<THashList*>(list_v0leg->FindObject(cut.GetName()));
        int nv0 = 0;
        for (auto& v0 : V0Photons_coll) {
          auto pos = v0.posTrack_as<MyMCV0Legs>();
          auto ele = v0.negTrack_as<MyMCV0Legs>();
          auto posmc = pos.template emmcparticle_as<aod::EMMCParticles>();
          auto elemc = ele.template emmcparticle_as<aod::EMMCParticles>();

          int photonid = FindCommonMotherFrom2Prongs(posmc, elemc, -11, 11, 22, mcparticles);
          if (photonid < 0) { // check swap, true electron is reconstructed as positron and vice versa.
            photonid = FindCommonMotherFrom2Prongs(posmc, elemc, 11, -11, 22, mcparticles);
          }

          if (photonid < 0) {
            continue;
          }
          auto mcphoton = mcparticles.iteratorAt(photonid);

          if (cut.IsSelected<MyMCV0Legs>(v0)) {
            o2::aod::emphotonhistograms::FillHistClass<EMHistType::kV0>(list_v0_cut, "", v0);
            if (IsPhysicalPrimary(mcphoton.emreducedmcevent(), mcphoton, mcparticles)) {
              reinterpret_cast<TH1F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hPt_Photon_Primary"))->Fill(v0.pt());
              reinterpret_cast<TH2F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hEtaPhi_Photon_Primary"))->Fill(v0.phi(), v0.eta());

              reinterpret_cast<TH2F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hPtGen_DeltaPtOverPtGen"))->Fill(mcphoton.pt(), (v0.pt() - mcphoton.pt()) / mcphoton.pt());
              reinterpret_cast<TH2F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hPtGen_DeltaEta"))->Fill(mcphoton.pt(), v0.eta() - mcphoton.eta());
              reinterpret_cast<TH2F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hPtGen_DeltaPhi"))->Fill(mcphoton.pt(), v0.phi() - mcphoton.phi());

            } else if (IsFromWD(mcphoton.emreducedmcevent(), mcphoton, mcparticles)) {
              reinterpret_cast<TH1F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hPt_Photon_FromWD"))->Fill(v0.pt());
              reinterpret_cast<TH2F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hEtaPhi_Photon_FromWD"))->Fill(v0.phi(), v0.eta());
            } else {
              // int mother_pdg = 0;
              // if(mcphoton.has_mothers()){
              //   auto mp = mcphoton.template mothers_first_as<aod::EMMCParticles>();
              //   mother_pdg = mp.pdgCode();
              // }
              // LOGF(info, "mcphoton.vx() = %f, mcphoton.vy() = %f, mcphoton.vz() = %f, mother_pdg = %d", mcphoton.vx(), mcphoton.vy(), mcphoton.vz(), mother_pdg);
              float rxy = sqrt(mcphoton.vx() * mcphoton.vx() + mcphoton.vy() * mcphoton.vy());
              reinterpret_cast<TH1F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hRZ_Photon_hs"))->Fill(mcphoton.vz(), rxy);
              reinterpret_cast<TH1F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hPt_Photon_hs"))->Fill(v0.pt());
              reinterpret_cast<TH1F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hEta_Photon_hs"))->Fill(v0.eta());
              reinterpret_cast<TH1F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hPhi_Photon_hs"))->Fill(v0.phi());
            }

            reinterpret_cast<TH1F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hConvPoint_diffX"))->Fill(elemc.vx(), v0.vx() - elemc.vx());
            reinterpret_cast<TH1F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hConvPoint_diffY"))->Fill(elemc.vy(), v0.vy() - elemc.vy());
            reinterpret_cast<TH1F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hConvPoint_diffZ"))->Fill(elemc.vz(), v0.vz() - elemc.vz());
            reinterpret_cast<TH1F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hConvPoint_diffX_recalc"))->Fill(elemc.vx(), v0.recalculatedVtxX() - elemc.vx());
            reinterpret_cast<TH1F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hConvPoint_diffY_recalc"))->Fill(elemc.vy(), v0.recalculatedVtxY() - elemc.vy());
            reinterpret_cast<TH1F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hConvPoint_diffZ_recalc"))->Fill(elemc.vz(), v0.recalculatedVtxZ() - elemc.vz());

            nv0++;
            reinterpret_cast<TH1F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hCorrTgl"))->Fill(ele.tgl(), pos.tgl());
            reinterpret_cast<TH1F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hCorrZ"))->Fill(ele.z(), pos.z());
            for (auto& leg : {pos, ele}) {
              o2::aod::emphotonhistograms::FillHistClass<EMHistType::kV0Leg>(list_v0leg_cut, "", leg);
            }
          }
        } // end of v0 loop
        reinterpret_cast<TH1F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hNgamma"))->Fill(nv0);
      } // end of cut loop
    }   // end of collision loop
  }     // end of process

  PresliceUnsorted<aod::EMMCParticles> perMcCollision = aod::emmcparticle::emreducedmceventId;
  void processGen(soa::Join<aod::EMReducedEvents, aod::EMReducedMCEventLabels> const& collisions, aod::EMReducedMCEvents const&, aod::EMMCParticles const& mcparticles)
  {
    // loop over mc stack and fill histograms for pure MC truth signals
    // all MC tracks which belong to the MC event corresponding to the current reconstructed event

    for (auto& collision : collisions) {
      auto mccollision = collision.emreducedmcevent();
      // LOGF(info, "mccollision.globalIndex() = %d", mccollision.globalIndex());

      reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hCollisionCounter"))->Fill(1.0);
      reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hZvtx_before"))->Fill(mccollision.posZ());
      if (!collision.sel8()) {
        continue;
      }
      reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hCollisionCounter"))->Fill(2.0);

      if (collision.numContrib() < 0.5) {
        continue;
      }
      reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hCollisionCounter"))->Fill(3.0);

      if (abs(collision.posZ()) > 10.0) {
        continue;
      }
      reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hCollisionCounter"))->Fill(4.0);
      reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hZvtx_after"))->Fill(mccollision.posZ());

      auto mctracks_coll = mcparticles.sliceBy(perMcCollision, mccollision.globalIndex());
      for (auto& mctrack : mctracks_coll) {
        if (abs(mctrack.y()) > maxY) {
          continue;
        }

        if (abs(mctrack.pdgCode()) == 22 && IsPhysicalPrimary(mctrack.emreducedmcevent(), mctrack, mcparticles)) {
          reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hPt_Photon"))->Fill(mctrack.pt());
          reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hY_Photon"))->Fill(mctrack.y());
          reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hPhi_Photon"))->Fill(mctrack.phi());
        }
      }
    }
  }

  void processDummy(aod::EMReducedEvents::iterator const& collision)
  {
    // do nothing
  }

  PROCESS_SWITCH(PCMQCMC, processQCMC, "run PCM QC in MC", false);
  PROCESS_SWITCH(PCMQCMC, processGen, "run generated information", false);
  PROCESS_SWITCH(PCMQCMC, processDummy, "Dummy function", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<PCMQCMC>(cfgc, TaskName{"pcm-qc-mc"})};
}
