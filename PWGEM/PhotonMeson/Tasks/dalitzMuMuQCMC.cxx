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

using MyCollisions = soa::Join<aod::EMReducedEvents, aod::EMReducedEventsMult, aod::EMReducedEventsCent, aod::EMReducedEventsNmumu, aod::EMReducedMCEventLabels>;
using MyCollision = MyCollisions::iterator;

using MyDalitzMuMus = soa::Join<aod::DalitzMuMus, aod::DalitzMuMuEMReducedEventIds>;
using MyDalitzMuMu = MyDalitzMuMus::iterator;

using MyTracks = soa::Join<aod::EMPrimaryMuons, aod::EMPrimaryMuonEMReducedEventIds, aod::EMPrimaryMuonsPrefilterBit>;
using MyMCTracks = soa::Join<MyTracks, aod::EMPrimaryMuonMCLabels>;

struct DalitzMuMuQCMC {
  Configurable<std::string> fConfigDalitzMuMuCuts{"cfgDalitzMuMuCuts", "nocut", "Comma separated list of dalitz mumu cuts"};
  std::vector<DalitzEECut> fDalitzMuMuCuts;

  OutputObj<THashList> fOutputEvent{"Event"};
  OutputObj<THashList> fOutputTrack{"Track"};
  OutputObj<THashList> fOutputDalitzMuMu{"DalitzMuMu"};
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

    o2::aod::emphotonhistograms::AddHistClass(fMainList, "DalitzMuMu");
    THashList* list_dalitzmumu = reinterpret_cast<THashList*>(fMainList->FindObject("DalitzMuMu"));

    o2::aod::emphotonhistograms::AddHistClass(fMainList, "Generated");
    THashList* list_gen = reinterpret_cast<THashList*>(fMainList->FindObject("Generated"));
    o2::aod::emphotonhistograms::DefineHistograms(list_gen, "Generated", "dimuon");

    for (const auto& cut : fDalitzMuMuCuts) {
      const char* cutname = cut.GetName();
      o2::aod::emphotonhistograms::AddHistClass(list_track, cutname);
      o2::aod::emphotonhistograms::AddHistClass(list_dalitzmumu, cutname);
    }

    // for single tracks
    for (auto& cut : fDalitzMuMuCuts) {
      std::string_view cutname = cut.GetName();
      THashList* list = reinterpret_cast<THashList*>(fMainList->FindObject("Track")->FindObject(cutname.data()));
      o2::aod::emphotonhistograms::DefineHistograms(list, "Track", "Mu");
    }

    // for DalitzMuMus
    for (auto& cut : fDalitzMuMuCuts) {
      std::string_view cutname = cut.GetName();
      THashList* list = reinterpret_cast<THashList*>(fMainList->FindObject("DalitzMuMu")->FindObject(cutname.data()));
      o2::aod::emphotonhistograms::DefineHistograms(list, "DalitzMuMu", "mc");
    }
  }

  void DefineCuts()
  {
    TString cutNamesStr = fConfigDalitzMuMuCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        const char* cutname = objArray->At(icut)->GetName();
        LOGF(info, "add cut : %s", cutname);
        fDalitzMuMuCuts.push_back(*dalitzeecuts::GetCut(cutname));
      }
    }
    LOGF(info, "Number of Dalitz cuts = %d", fDalitzMuMuCuts.size());
  }

  void init(InitContext& context)
  {
    DefineCuts();
    addhistograms(); // please call this after DefineCuts();

    fOutputEvent.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("Event")));
    fOutputTrack.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("Track")));
    fOutputDalitzMuMu.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("DalitzMuMu")));
    fOutputGen.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("Generated")));
  }

  template <typename TTrack, typename TMCTracks>
  int FindLF(TTrack const& posmc, TTrack const& elemc, TMCTracks const& mcparticles)
  {
    int arr[] = {
      FindCommonMotherFrom2Prongs(posmc, elemc, -13, 13, 111, mcparticles), FindCommonMotherFrom2Prongs(posmc, elemc, -13, 13, 221, mcparticles), FindCommonMotherFrom2Prongs(posmc, elemc, -13, 13, 331, mcparticles), FindCommonMotherFrom2Prongs(posmc, elemc, -13, 13, 113, mcparticles), FindCommonMotherFrom2Prongs(posmc, elemc, -13, 13, 223, mcparticles), FindCommonMotherFrom2Prongs(posmc, elemc, -13, 13, 333, mcparticles)};
    int size = sizeof(arr) / sizeof(*arr);
    int max = *std::max_element(arr, arr + size);
    return max;
  }

  Partition<MyDalitzMuMus> uls_pairs = o2::aod::dalitzmumu::sign == 0;

  SliceCache cache;
  Preslice<MyDalitzMuMus> perCollision = aod::dalitzmumu::emreducedeventId;

  std::vector<uint64_t> used_trackIds;

  void processQCMC(MyCollisions const& collisions, MyDalitzMuMus const& dileptons, MyMCTracks const& tracks, aod::EMMCParticles const& mcparticles, aod::EMReducedMCEvents const&)
  {
    THashList* list_ev = static_cast<THashList*>(fMainList->FindObject("Event"));
    THashList* list_dalitzmumu = static_cast<THashList*>(fMainList->FindObject("DalitzMuMu"));
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

      auto uls_pairs_per_coll = uls_pairs->sliceByCached(o2::aod::dalitzmumu::emreducedeventId, collision.globalIndex(), cache);

      for (const auto& cut : fDalitzMuMuCuts) {
        THashList* list_dalitzmumu_cut = static_cast<THashList*>(list_dalitzmumu->FindObject(cut.GetName()));
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

          // LOGF(info, "posmc.pdgCode() = %d, , posmc.has_mothers()() = %d | elemc.pdgCode() = %d, elemc.has_mothers() = %d", posmc.pdgCode(), posmc.has_mothers(), elemc.pdgCode(), elemc.has_mothers());

          int mother_id = FindLF(posmc, elemc, mcparticles);
          if (mother_id < 0) {
            continue;
          }
          if (mother_id > 0) {
            auto mcmother = mcparticles.iteratorAt(mother_id);
            if (IsPhysicalPrimary(mcmother.emreducedmcevent(), mcmother, mcparticles)) {
              values[0] = uls_pair.mass();
              values[1] = uls_pair.pt();
              values[2] = uls_pair.dcaXY();
              values[3] = uls_pair.phiv();
              reinterpret_cast<THnSparseF*>(list_dalitzmumu_cut->FindObject("hs_dilepton_uls_same"))->Fill(values);

              nuls++;
              for (auto& track : {pos, ele}) {
                if (std::find(used_trackIds.begin(), used_trackIds.end(), track.globalIndex()) == used_trackIds.end()) {
                  o2::aod::emphotonhistograms::FillHistClass<EMHistType::kTrack>(list_track_cut, "", track);
                  used_trackIds.emplace_back(track.globalIndex());
                }
              }
            } // end of LF
          }
        } // end of uls pair loop
        reinterpret_cast<TH1F*>(list_dalitzmumu_cut->FindObject("hNpair_uls"))->Fill(nuls);

        used_trackIds.clear();
        used_trackIds.shrink_to_fit();
      } // end of cut loop
    }   // end of collision loop
  }     // end of process
  PROCESS_SWITCH(DalitzMuMuQCMC, processQCMC, "run Dalitz QC", true);

  Configurable<float> min_mcPt{"min_mcPt", 0.05, "min. MC pT"};
  Configurable<float> max_mcPt{"max_mcPt", 0.6, "max. MC pT"};
  Configurable<float> max_mcEta{"max_mcEta", 0.9, "max. MC eta"};
  Partition<aod::EMMCParticles> posTracks = o2::aod::mcparticle::pdgCode == -13; // mu+
  Partition<aod::EMMCParticles> negTracks = o2::aod::mcparticle::pdgCode == +13; // mu-
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
          ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassMuon);
          ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassMuon);
          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          reinterpret_cast<TH2F*>(fMainList->FindObject("Generated")->FindObject("hMvsPt"))->Fill(v12.M(), v12.Pt());
        } // end of LF
      }   // end of true ULS pair loop
    }     // end of collision loop
  }
  PROCESS_SWITCH(DalitzMuMuQCMC, processGen, "run genrated info", true);

  void processDummy(MyCollisions const& collisions) {}
  PROCESS_SWITCH(DalitzMuMuQCMC, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DalitzMuMuQCMC>(cfgc, TaskName{"dalitz-mumu-qc-mc"})};
}
