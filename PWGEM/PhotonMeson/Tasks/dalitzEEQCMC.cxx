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
// This code runs loop over dalitz ee table for dalitz QC.
//    Please write to: daiki.sekihata@cern.ch

#include <array>
#include "TString.h"
#include "THashList.h"
#include "TDirectory.h"
#include "Math/Vector4D.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"

#include "Common/Core/RecoDecay.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Core/DalitzEECut.h"
#include "PWGEM/PhotonMeson/Core/CutsLibrary.h"
#include "PWGEM/PhotonMeson/Core/HistogramsLibrary.h"
#include "PWGEM/PhotonMeson/Utils/MCUtilities.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using std::array;

using MyCollisions = soa::Join<aod::EMReducedEvents, aod::EMReducedEventsMult, aod::EMReducedEventsCent, aod::EMReducedEventsNee, aod::EMReducedMCEventLabels>;
using MyCollision = MyCollisions::iterator;

using MyDalitzEEs = soa::Join<aod::DalitzEEs, aod::DalitzEEEMReducedEventIds>;
using MyDalitzEE = MyDalitzEEs::iterator;

using MyTracks = soa::Join<aod::EMPrimaryElectrons, aod::EMPrimaryElectronEMReducedEventIds, aod::EMPrimaryElectronsPrefilterBit>;
using MyMCTracks = soa::Join<MyTracks, aod::EMPrimaryElectronMCLabels>;

struct DalitzEEQCMC {
  Configurable<std::string> fConfigDalitzEECuts{"cfgDalitzEECuts", "mee_all_tpchadrejortofreq_lowB,nocut", "Comma separated list of dalitz ee cuts"};
  std::vector<DalitzEECut> fDalitzEECuts;

  OutputObj<THashList> fOutputEvent{"Event"};
  OutputObj<THashList> fOutputTrack{"Track"};
  OutputObj<THashList> fOutputDalitzEE{"DalitzEE"};
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

    o2::aod::emphotonhistograms::AddHistClass(fMainList, "Track");
    THashList* list_track = reinterpret_cast<THashList*>(fMainList->FindObject("Track"));

    o2::aod::emphotonhistograms::AddHistClass(fMainList, "DalitzEE");
    THashList* list_dalitzee = reinterpret_cast<THashList*>(fMainList->FindObject("DalitzEE"));

    o2::aod::emphotonhistograms::AddHistClass(fMainList, "Generated");
    THashList* list_gen = reinterpret_cast<THashList*>(fMainList->FindObject("Generated"));
    o2::aod::emphotonhistograms::DefineHistograms(list_gen, "Generated", "dielectron");

    for (const auto& cut : fDalitzEECuts) {
      const char* cutname = cut.GetName();
      o2::aod::emphotonhistograms::AddHistClass(list_track, cutname);
      o2::aod::emphotonhistograms::AddHistClass(list_dalitzee, cutname);
    }

    // for single tracks
    for (auto& cut : fDalitzEECuts) {
      std::string_view cutname = cut.GetName();
      THashList* list = reinterpret_cast<THashList*>(fMainList->FindObject("Track")->FindObject(cutname.data()));
      o2::aod::emphotonhistograms::DefineHistograms(list, "Track");
    }

    // for DalitzEEs
    for (auto& cut : fDalitzEECuts) {
      std::string_view cutname = cut.GetName();
      THashList* list = reinterpret_cast<THashList*>(fMainList->FindObject("DalitzEE")->FindObject(cutname.data()));
      o2::aod::emphotonhistograms::DefineHistograms(list, "DalitzEE", "mc");
    }
  }

  void DefineCuts()
  {
    TString cutNamesStr = fConfigDalitzEECuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        const char* cutname = objArray->At(icut)->GetName();
        LOGF(info, "add cut : %s", cutname);
        fDalitzEECuts.push_back(*dalitzeecuts::GetCut(cutname));
      }
    }
    LOGF(info, "Number of Dalitz cuts = %d", fDalitzEECuts.size());
  }

  void init(InitContext& context)
  {
    DefineCuts();
    addhistograms(); // please call this after DefinCuts();

    fOutputEvent.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("Event")));
    fOutputTrack.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("Track")));
    fOutputDalitzEE.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("DalitzEE")));
    fOutputGen.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("Generated")));
  }

  template <typename TTrack, typename TMCTracks>
  int FindLF(TTrack const& posmc, TTrack const& elemc, TMCTracks const& mcparticles)
  {
    int arr[] = {
      FindCommonMotherFrom2Prongs(posmc, elemc, -11, 11, 111, mcparticles), FindCommonMotherFrom2Prongs(posmc, elemc, -11, 11, 221, mcparticles), FindCommonMotherFrom2Prongs(posmc, elemc, -11, 11, 331, mcparticles), FindCommonMotherFrom2Prongs(posmc, elemc, -11, 11, 113, mcparticles), FindCommonMotherFrom2Prongs(posmc, elemc, -11, 11, 223, mcparticles), FindCommonMotherFrom2Prongs(posmc, elemc, -11, 11, 333, mcparticles)};
    int size = sizeof(arr) / sizeof(*arr);
    int max = *std::max_element(arr, arr + size);
    return max;
  }

  Partition<MyDalitzEEs> uls_pairs = o2::aod::dalitzee::sign == 0;

  SliceCache cache;
  Preslice<MyDalitzEEs> perCollision = aod::dalitzee::emreducedeventId;

  std::vector<uint64_t> used_trackIds;

  void processQCMC(MyCollisions const& collisions, MyDalitzEEs const& dileptons, MyMCTracks const& tracks, aod::EMMCParticles const& mcparticles, aod::EMReducedMCEvents const&)
  {
    THashList* list_ev = static_cast<THashList*>(fMainList->FindObject("Event"));
    THashList* list_dalitzee = static_cast<THashList*>(fMainList->FindObject("DalitzEE"));
    THashList* list_track = static_cast<THashList*>(fMainList->FindObject("Track"));
    double values[4] = {0, 0, 0, 0};

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

      auto uls_pairs_per_coll = uls_pairs->sliceByCached(o2::aod::dalitzee::emreducedeventId, collision.globalIndex(), cache);

      for (const auto& cut : fDalitzEECuts) {
        THashList* list_dalitzee_cut = static_cast<THashList*>(list_dalitzee->FindObject(cut.GetName()));
        THashList* list_track_cut = static_cast<THashList*>(list_track->FindObject(cut.GetName()));
        used_trackIds.reserve(uls_pairs_per_coll.size() * 2);

        int nuls = 0;
        for (auto& uls_pair : uls_pairs_per_coll) {
          auto pos = uls_pair.template posTrack_as<MyMCTracks>();
          auto ele = uls_pair.template negTrack_as<MyMCTracks>();

          if (!cut.IsSelected<MyMCTracks>(uls_pair)) {
            continue;
          }

          auto posmc = pos.template emmcparticle_as<aod::EMMCParticles>();
          auto elemc = ele.template emmcparticle_as<aod::EMMCParticles>();
          int mother_id = FindLF(posmc, elemc, mcparticles);
          int photonid = FindCommonMotherFrom2Prongs(posmc, elemc, -11, 11, 22, mcparticles);
          if (mother_id < 0 && photonid < 0) {
            continue;
          }
          if (mother_id > 0) {
            auto mcmother = mcparticles.iteratorAt(mother_id);
            if (IsPhysicalPrimary(mcmother.emreducedmcevent(), mcmother, mcparticles)) {
              values[0] = uls_pair.mass();
              values[1] = uls_pair.pt();
              values[2] = uls_pair.dcaXY();
              values[3] = uls_pair.phiv();
              reinterpret_cast<THnSparseF*>(list_dalitzee_cut->FindObject("hs_dilepton_uls_same"))->Fill(values);

              if (mcmother.pdgCode() == 111) {
                reinterpret_cast<TH2F*>(list_dalitzee_cut->FindObject("hMvsPhiV_Pi0"))->Fill(uls_pair.phiv(), uls_pair.mass());
              }

              nuls++;
              for (auto& track : {pos, ele}) {
                if (std::find(used_trackIds.begin(), used_trackIds.end(), track.globalIndex()) == used_trackIds.end()) {
                  o2::aod::emphotonhistograms::FillHistClass<EMHistType::kTrack>(list_track_cut, "", track);
                  used_trackIds.emplace_back(track.globalIndex());
                }
              }
            } // end of LF
          } else if (photonid > 0) {
            auto mcphoton = mcparticles.iteratorAt(photonid);
            if (IsPhysicalPrimary(mcphoton.emreducedmcevent(), mcphoton, mcparticles) && IsEleFromPC(elemc, mcparticles) && IsEleFromPC(posmc, mcparticles)) {
              reinterpret_cast<TH2F*>(list_dalitzee_cut->FindObject("hMvsPhiV_Photon"))->Fill(uls_pair.phiv(), uls_pair.mass());
            }
          }
        } // end of uls pair loop
        reinterpret_cast<TH1F*>(list_dalitzee_cut->FindObject("hNpair_uls"))->Fill(nuls);

        used_trackIds.clear();
        used_trackIds.shrink_to_fit();
      } // end of cut loop
    }   // end of collision loop
  }     // end of process
  PROCESS_SWITCH(DalitzEEQCMC, processQCMC, "run Dalitz QC", true);

  Configurable<float> min_mcPt{"min_mcPt", 0.05, "min. MC pT"};
  Configurable<float> max_mcPt{"max_mcPt", 1e+10, "max. MC pT"};
  Configurable<float> max_mcEta{"max_mcEta", 0.9, "max. MC eta"};
  Partition<aod::EMMCParticles> posTracks = o2::aod::mcparticle::pdgCode == -11; // e+
  Partition<aod::EMMCParticles> negTracks = o2::aod::mcparticle::pdgCode == +11; // e-
  PresliceUnsorted<aod::EMMCParticles> perMcCollision = aod::emmcparticle::emreducedmceventId;
  void processGen(MyCollisions const& collisions, aod::EMReducedMCEvents const&, aod::EMMCParticles const& mcparticles)
  {
    // loop over mc stack and fill histograms for pure MC truth signals
    // all MC tracks which belong to the MC event corresponding to the current reconstructed event
    for (auto& collision : collisions) {
      auto mccollision = collision.emreducedmcevent();
      // LOGF(info, "mccollision.globalIndex() = %d", mccollision.globalIndex());
      // auto mctracks_coll = mcparticles.sliceBy(perMcCollision, mccollision.globalIndex());

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
      auto posTracks_per_coll = posTracks->sliceByCachedUnsorted(o2::aod::emmcparticle::emreducedmceventId, mccollision.globalIndex(), cache);
      auto negTracks_per_coll = negTracks->sliceByCachedUnsorted(o2::aod::emmcparticle::emreducedmceventId, mccollision.globalIndex(), cache);
      for (auto& [t1, t2] : combinations(CombinationsFullIndexPolicy(posTracks_per_coll, negTracks_per_coll))) {
        // LOGF(info, "pdg1 = %d, pdg2 = %d", t1.pdgCode(), t2.pdgCode());

        if (t1.pt() < min_mcPt || max_mcPt < t1.pt() || abs(t1.eta()) > max_mcEta) {
          continue;
        }
        if (t2.pt() < min_mcPt || max_mcPt < t2.pt() || abs(t2.eta()) > max_mcEta) {
          continue;
        }

        int mother_id = FindLF(t1, t2, mcparticles);
        if (mother_id < 0) {
          continue;
        }
        auto mcmother = mcparticles.iteratorAt(mother_id);
        if (IsPhysicalPrimary(mcmother.emreducedmcevent(), mcmother, mcparticles)) {
          ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassElectron);
          ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassElectron);
          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          reinterpret_cast<TH2F*>(fMainList->FindObject("Generated")->FindObject("hMvsPt"))->Fill(v12.M(), v12.Pt());
        } // end of LF
      }   // end of true ULS pair loop
    }     // end of collision loop
  }
  PROCESS_SWITCH(DalitzEEQCMC, processGen, "run genrated info", true);

  void processDummy(MyCollisions const& collisions) {}
  PROCESS_SWITCH(DalitzEEQCMC, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DalitzEEQCMC>(cfgc, TaskName{"dalitz-ee-qc-mc"})};
}
