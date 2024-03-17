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
using namespace o2::aod::pwgem::mcutil;
using namespace o2::aod::pwgem::photon;
using std::array;

using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent, aod::EMMCEventLabels>;
using MyCollision = MyCollisions::iterator;

using MyV0Photons = soa::Join<aod::V0PhotonsKF, aod::V0KFEMEventIds>;
using MyV0Photon = MyV0Photons::iterator;

struct PCMQCMC {
  using MyMCV0Legs = soa::Join<aod::V0Legs, aod::V0LegMCLabels>;

  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};

  Configurable<std::string> fConfigPCMCuts{"cfgPCMCuts", "qc,qc_ITSTPC,qc_ITSonly,qc_TPConly,wwire_ib,nocut", "Comma separated list of v0 photon cuts"};
  Configurable<float> maxY{"maxY", 0.9f, "maximum rapidity for generated particles"};
  Configurable<float> maxRgen{"maxRgen", 100.f, "maximum radius for generated particles"};
  Configurable<float> margin_z_mc{"margin_z_mc", 7.0, "margin for z cut in cm for MC"};

  std::vector<V0PhotonCut> fPCMCuts;

  Configurable<std::string> fConfigEMEventCut{"cfgEMEventCut", "minbias", "em event cut"}; // only 1 event cut per wagon
  EMEventCut fEMEventCut;
  static constexpr std::string_view event_types[2] = {"before", "after"};

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
    o2::aod::pwgem::photon::histogram::AddHistClass(fMainList, "Event");
    THashList* list_ev = reinterpret_cast<THashList*>(fMainList->FindObject("Event"));
    for (const auto& evtype : event_types) {
      THashList* list_ev_type = reinterpret_cast<THashList*>(o2::aod::pwgem::photon::histogram::AddHistClass(list_ev, evtype.data()));
      o2::aod::pwgem::photon::histogram::DefineHistograms(list_ev_type, "Event", evtype.data());
    }

    o2::aod::pwgem::photon::histogram::AddHistClass(fMainList, "V0Leg");
    THashList* list_v0leg = reinterpret_cast<THashList*>(fMainList->FindObject("V0Leg"));

    o2::aod::pwgem::photon::histogram::AddHistClass(fMainList, "V0");
    THashList* list_v0 = reinterpret_cast<THashList*>(fMainList->FindObject("V0"));

    o2::aod::pwgem::photon::histogram::AddHistClass(fMainList, "Generated");
    THashList* list_gen = reinterpret_cast<THashList*>(fMainList->FindObject("Generated"));
    o2::aod::pwgem::photon::histogram::DefineHistograms(list_gen, "Generated", "Photon");

    for (const auto& cut : fPCMCuts) {
      const char* cutname = cut.GetName();
      o2::aod::pwgem::photon::histogram::AddHistClass(list_v0leg, cutname);
      o2::aod::pwgem::photon::histogram::AddHistClass(list_v0, cutname);
    }

    // for single tracks
    for (auto& cut : fPCMCuts) {
      std::string_view cutname = cut.GetName();
      THashList* list = reinterpret_cast<THashList*>(fMainList->FindObject("V0Leg")->FindObject(cutname.data()));
      o2::aod::pwgem::photon::histogram::DefineHistograms(list, "V0Leg", "mc");
    }

    // for V0s
    for (auto& cut : fPCMCuts) {
      std::string_view cutname = cut.GetName();
      THashList* list = reinterpret_cast<THashList*>(fMainList->FindObject("V0")->FindObject(cutname.data()));
      o2::aod::pwgem::photon::histogram::DefineHistograms(list, "V0", "mc");
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
    TString ev_cut_name = fConfigEMEventCut.value;
    fEMEventCut = *eventcuts::GetCut(ev_cut_name.Data());

    fOutputEvent.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("Event")));
    fOutputV0Leg.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("V0Leg")));
    fOutputV0.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("V0")));
    fOutputGen.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("Generated")));
  }

  Partition<MyCollisions> grouped_collisions = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax); // this goes to same event.

  Preslice<MyV0Photons> perCollision = aod::v0photonkf::emeventId;
  void processQCMC(MyCollisions const& collisions, MyV0Photons const& v0photons, MyMCV0Legs const& v0legs, aod::EMMCParticles const& mcparticles, aod::EMMCEvents const&)
  {
    THashList* list_ev_before = static_cast<THashList*>(fMainList->FindObject("Event")->FindObject(event_types[0].data()));
    THashList* list_ev_after = static_cast<THashList*>(fMainList->FindObject("Event")->FindObject(event_types[1].data()));
    THashList* list_v0 = static_cast<THashList*>(fMainList->FindObject("V0"));
    THashList* list_v0leg = static_cast<THashList*>(fMainList->FindObject("V0Leg"));

    for (auto& collision : grouped_collisions) {
      const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      o2::aod::pwgem::photon::histogram::FillHistClass<EMHistType::kEvent>(list_ev_before, "", collision);
      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      o2::aod::pwgem::photon::histogram::FillHistClass<EMHistType::kEvent>(list_ev_after, "", collision);
      reinterpret_cast<TH1F*>(list_ev_before->FindObject("hCollisionCounter"))->Fill("accepted", 1.f);
      reinterpret_cast<TH1F*>(list_ev_after->FindObject("hCollisionCounter"))->Fill("accepted", 1.f);

      auto V0Photons_coll = v0photons.sliceBy(perCollision, collision.globalIndex());
      for (const auto& cut : fPCMCuts) {
        THashList* list_v0_cut = static_cast<THashList*>(list_v0->FindObject(cut.GetName()));
        THashList* list_v0leg_cut = static_cast<THashList*>(list_v0leg->FindObject(cut.GetName()));
        int nv0 = 0;
        for (auto& v0 : V0Photons_coll) {
          auto pos = v0.posTrack_as<MyMCV0Legs>();
          auto ele = v0.negTrack_as<MyMCV0Legs>();
          auto posmc = pos.template emmcparticle_as<aod::EMMCParticles>();
          auto elemc = ele.template emmcparticle_as<aod::EMMCParticles>();

          if (cut.IsSelected<MyMCV0Legs>(v0)) {
            reinterpret_cast<TH1F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hPt_Photon_Candidate"))->Fill(v0.pt());
            reinterpret_cast<TH2F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hEtaPhi_Photon_Candidate"))->Fill(v0.phi(), v0.eta());

            int photonid = FindCommonMotherFrom2Prongs(posmc, elemc, -11, 11, 22, mcparticles);
            if (photonid < 0) {
              continue;
            }
            auto mcphoton = mcparticles.iteratorAt(photonid);
            float rxy_rec = sqrt(v0.vx() * v0.vx() + v0.vy() * v0.vy());
            float rxy_mc = sqrt(posmc.vx() * posmc.vx() + posmc.vy() * posmc.vy());
            float eta_cp = std::atanh(v0.vz() / sqrt(pow(v0.vx(), 2) + pow(v0.vy(), 2) + pow(v0.vz(), 2)));

            o2::aod::pwgem::photon::histogram::FillHistClass<EMHistType::kV0>(list_v0_cut, "", v0);
            if (IsPhysicalPrimary(mcphoton.emmcevent(), mcphoton, mcparticles)) {
              reinterpret_cast<TH1F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hPt_Photon_Primary"))->Fill(v0.pt());
              reinterpret_cast<TH2F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hEtaPhi_Photon_Primary"))->Fill(v0.phi(), v0.eta());
              reinterpret_cast<TH2F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hXY_Photon_Primary"))->Fill(v0.vx(), v0.vy());
              reinterpret_cast<TH2F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hXY_Photon_Primary_MC"))->Fill(posmc.vx(), posmc.vy());
              reinterpret_cast<TH2F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hRZ_Photon_Primary"))->Fill(v0.vz(), rxy_rec);
              reinterpret_cast<TH2F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hRZ_Photon_Primary_MC"))->Fill(posmc.vz(), rxy_mc);

              reinterpret_cast<TH2F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hPtGen_DeltaPtOverPtGen"))->Fill(mcphoton.pt(), (v0.pt() - mcphoton.pt()) / mcphoton.pt());
              reinterpret_cast<TH2F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hPtGen_DeltaEta"))->Fill(mcphoton.pt(), v0.eta() - mcphoton.eta());
              reinterpret_cast<TH2F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hPtGen_DeltaPhi"))->Fill(mcphoton.pt(), v0.phi() - mcphoton.phi());

              reinterpret_cast<TH2F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hEtaRec_DeltaPtOverPtGen"))->Fill(eta_cp, (v0.pt() - mcphoton.pt()) / mcphoton.pt());
              reinterpret_cast<TH2F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hEtaRec_DeltaEta"))->Fill(eta_cp, v0.eta() - mcphoton.eta());
              reinterpret_cast<TH2F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hEtaRec_DeltaPhi"))->Fill(eta_cp, v0.phi() - mcphoton.phi());

            } else if (IsFromWD(mcphoton.emmcevent(), mcphoton, mcparticles)) {
              reinterpret_cast<TH1F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hPt_Photon_FromWD"))->Fill(v0.pt());
              reinterpret_cast<TH2F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hEtaPhi_Photon_FromWD"))->Fill(v0.phi(), v0.eta());
            } else {
              // int mother_pdg = 0;
              // if(mcphoton.has_mothers()){
              //   auto mp = mcphoton.template mothers_first_as<aod::EMMCParticles>();
              //   mother_pdg = mp.pdgCode();
              // }
              // LOGF(info, "mcphoton.vx() = %f, mcphoton.vy() = %f, mcphoton.vz() = %f, mother_pdg = %d", mcphoton.vx(), mcphoton.vy(), mcphoton.vz(), mother_pdg);
              float rxy_photon_hs = sqrt(mcphoton.vx() * mcphoton.vx() + mcphoton.vy() * mcphoton.vy());
              reinterpret_cast<TH2F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hRZ_Photon_hs"))->Fill(mcphoton.vz(), rxy_photon_hs);
              reinterpret_cast<TH1F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hPt_Photon_hs"))->Fill(v0.pt());
              reinterpret_cast<TH1F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hEta_Photon_hs"))->Fill(v0.eta());
              reinterpret_cast<TH1F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hPhi_Photon_hs"))->Fill(v0.phi());
            }

            reinterpret_cast<TH1F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hConvPoint_diffX"))->Fill(elemc.vx(), v0.vx() - elemc.vx());
            reinterpret_cast<TH1F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hConvPoint_diffY"))->Fill(elemc.vy(), v0.vy() - elemc.vy());
            reinterpret_cast<TH1F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hConvPoint_diffZ"))->Fill(elemc.vz(), v0.vz() - elemc.vz());
            nv0++;

            for (auto& leg : {pos, ele}) {
              o2::aod::pwgem::photon::histogram::FillHistClass<EMHistType::kV0Leg>(list_v0leg_cut, "", leg);
              auto mcleg = leg.template emmcparticle_as<aod::EMMCParticles>();
              reinterpret_cast<TH2F*>(fMainList->FindObject("V0Leg")->FindObject(cut.GetName())->FindObject("hPtGen_DeltaPtOverPtGen"))->Fill(mcleg.pt(), (leg.pt() - mcleg.pt()) / mcleg.pt());
              reinterpret_cast<TH2F*>(fMainList->FindObject("V0Leg")->FindObject(cut.GetName())->FindObject("hPtGen_DeltaEta"))->Fill(mcleg.pt(), leg.eta() - mcleg.eta());
              reinterpret_cast<TH2F*>(fMainList->FindObject("V0Leg")->FindObject(cut.GetName())->FindObject("hPtGen_DeltaPhi"))->Fill(mcleg.pt(), leg.phi() - mcleg.phi());
              reinterpret_cast<TH2F*>(fMainList->FindObject("V0Leg")->FindObject(cut.GetName())->FindObject("hEtaRec_DeltaPtOverPtGen"))->Fill(eta_cp, (leg.pt() - mcleg.pt()) / mcleg.pt());
              reinterpret_cast<TH2F*>(fMainList->FindObject("V0Leg")->FindObject(cut.GetName())->FindObject("hEtaRec_DeltaEta"))->Fill(eta_cp, leg.eta() - mcleg.eta());
              reinterpret_cast<TH2F*>(fMainList->FindObject("V0Leg")->FindObject(cut.GetName())->FindObject("hEtaRec_DeltaPhi"))->Fill(eta_cp, leg.phi() - mcleg.phi());
            }
          }
        } // end of v0 loop
        reinterpret_cast<TH1F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hNgamma"))->Fill(nv0);
      } // end of cut loop
    }   // end of collision loop
  }     // end of process

  PresliceUnsorted<aod::EMMCParticles> perMcCollision = aod::emmcparticle::emmceventId;
  void processGen(MyCollisions const& collisions, aod::EMMCEvents const&, aod::EMMCParticles const& mcparticles)
  {
    // loop over mc stack and fill histograms for pure MC truth signals
    // all MC tracks which belong to the MC event corresponding to the current reconstructed event

    for (auto& collision : grouped_collisions) {
      const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }
      auto mccollision = collision.emmcevent();
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

      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }

      reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hCollisionCounter"))->Fill(4.0);
      reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hZvtx_after"))->Fill(mccollision.posZ());

      auto mctracks_coll = mcparticles.sliceBy(perMcCollision, mccollision.globalIndex());
      for (auto& mctrack : mctracks_coll) {
        if (abs(mctrack.y()) > maxY) {
          continue;
        }

        if (abs(mctrack.pdgCode()) == 22 && IsPhysicalPrimary(mctrack.emmcevent(), mctrack, mcparticles)) {
          reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hPt_Photon"))->Fill(mctrack.pt());
          reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hY_Photon"))->Fill(mctrack.y());
          reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hPhi_Photon"))->Fill(mctrack.phi());

          bool is_ele_fromPC = false;
          bool is_pos_fromPC = false;
          auto daughtersIds = mctrack.daughtersIds(); // always size = 2. first and last index. one should run loop from the first index to the last index.
          for (auto& daughterId : daughtersIds) {
            if (daughterId < 0) {
              continue;
            }
            auto daughter = mcparticles.iteratorAt(daughterId); // always electron and positron
            float rxy_gen_e = sqrt(pow(daughter.vx(), 2) + pow(daughter.vy(), 2));
            if (rxy_gen_e > maxRgen || rxy_gen_e < abs(daughter.vz()) * TMath::Tan(2 * TMath::ATan(TMath::Exp(-maxY))) - margin_z_mc) {
              continue;
            }

            if (daughter.pdgCode() == 11) { // electron from photon conversion
              is_ele_fromPC = true;
            } else if (daughter.pdgCode() == -11) { // positron from photon conversion
              is_pos_fromPC = true;
            }
          }                                     // end of daughter loop
          if (is_ele_fromPC && is_pos_fromPC) { // ele and pos from photon conversion
            reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hPt_ConvertedPhoton"))->Fill(mctrack.pt());
            reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hY_ConvertedPhoton"))->Fill(mctrack.y());
            reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hPhi_ConvertedPhoton"))->Fill(mctrack.phi());

            auto daughter = mcparticles.iteratorAt(daughtersIds[0]); // choose ele or pos.
            float rxy_gen_e = sqrt(pow(daughter.vx(), 2) + pow(daughter.vy(), 2));
            reinterpret_cast<TH2F*>(fMainList->FindObject("Generated")->FindObject("hPhotonRZ"))->Fill(daughter.vz(), rxy_gen_e);
            reinterpret_cast<TH2F*>(fMainList->FindObject("Generated")->FindObject("hPhotonRxy"))->Fill(daughter.vx(), daughter.vy());
            reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hPhotonPhivsRxy"))->Fill(daughter.phi(), rxy_gen_e);
          }
        }
      } // end of mctrack loop per collision
    }   // end of collision loop
  }

  void processDummy(MyCollisions const& collisions)
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
