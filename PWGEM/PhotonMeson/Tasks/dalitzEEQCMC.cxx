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
#include "Framework/AnalysisDataModel.h"
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

using MyDalitzEEs = soa::Join<aod::DalitzEEs, aod::DalitzEEEMReducedEventIds>;
using MyDalitzEE = MyDalitzEEs::iterator;

struct DalitzEEQCMC {
  Configurable<std::string> fConfigDalitzEECuts{"cfgDalitzEECuts", "mee_all_tpchadrejortofreq_lowB,nocut", "Comma separated list of dalitz ee cuts"};
  std::vector<DalitzEECut> fDalitzEECuts;

  OutputObj<THashList> fOutputEvent{"Event"};
  OutputObj<THashList> fOutputTrack{"Track"};
  OutputObj<THashList> fOutputDalitzEE{"DalitzEE"};
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
      o2::aod::emphotonhistograms::DefineHistograms(list, "DalitzEE");
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
  using MyMCTracks = soa::Join<aod::EMPrimaryTracks, aod::EMPrimaryTrackMCLabels>;
  void processQCMC(aod::EMReducedEvents const& collisions, MyDalitzEEs const& dileptons, MyMCTracks const& tracks, aod::EMMCParticles const& mcparticles, aod::EMReducedMCEvents const&)
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

        int nuls = 0;
        for (auto& uls_pair : uls_pairs_per_coll) {
          auto pos = uls_pair.posTrack_as<MyMCTracks>();
          auto ele = uls_pair.negTrack_as<MyMCTracks>();
          auto posmc = pos.template emmcparticle_as<aod::EMMCParticles>();
          auto elemc = ele.template emmcparticle_as<aod::EMMCParticles>();

          int mother_id = FindLF(posmc, elemc, mcparticles);
          if (mother_id < 0) {
            continue;
          }
          auto mcmother = mcparticles.iteratorAt(mother_id);

          if (!IsPhysicalPrimary(mcmother.emreducedmcevent(), mcmother, mcparticles)) {
            continue;
          }
          if (cut.IsSelected<MyMCTracks>(uls_pair)) {
            values[0] = uls_pair.mee();
            values[1] = uls_pair.pt();
            values[2] = uls_pair.dcaeeXY();
            values[3] = uls_pair.phiv();
            reinterpret_cast<THnSparseF*>(list_dalitzee_cut->FindObject("hs_dilepton_uls"))->Fill(values);
            nuls++;
            for (auto& leg : {pos, ele}) {
              o2::aod::emphotonhistograms::FillHistClass<EMHistType::kTrack>(list_track_cut, "", leg);
            }
          }
        } // end of uls pair loop
        reinterpret_cast<TH1F*>(list_dalitzee_cut->FindObject("hNpair_uls"))->Fill(nuls);
      } // end of cut loop
    }   // end of collision loop
  }     // end of process
  PROCESS_SWITCH(DalitzEEQCMC, processQCMC, "run Dalitz QC", true);

  void processDummy(aod::EMReducedEvents::iterator const& collision) {}
  PROCESS_SWITCH(DalitzEEQCMC, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DalitzEEQCMC>(cfgc, TaskName{"dalitz-ee-qc-mc"})};
}
