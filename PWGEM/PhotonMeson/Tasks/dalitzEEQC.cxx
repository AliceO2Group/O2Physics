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

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using std::array;

using MyDalitzEEs = soa::Join<aod::DalitzEEs, aod::DalitzEEEMReducedEventIds>;
using MyDalitzEE = MyDalitzEEs::iterator;

struct DalitzEEQC {
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

  Partition<MyDalitzEEs> uls_pairs = o2::aod::dalitzee::sign == 0;
  Partition<MyDalitzEEs> lspp_pairs = o2::aod::dalitzee::sign == +1;
  Partition<MyDalitzEEs> lsmm_pairs = o2::aod::dalitzee::sign == -1;

  SliceCache cache;
  Preslice<MyDalitzEEs> perCollision = aod::dalitzee::emreducedeventId;
  void processQC(aod::EMReducedEvents const& collisions, MyDalitzEEs const& dileptons, aod::EMPrimaryTracks const& tracks)
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

      // auto dileptons_coll = dileptons.sliceBy(perCollision, collision.collisionId());

      auto uls_pairs_per_coll = uls_pairs->sliceByCached(o2::aod::dalitzee::emreducedeventId, collision.globalIndex(), cache);
      auto lspp_pairs_per_coll = lspp_pairs->sliceByCached(o2::aod::dalitzee::emreducedeventId, collision.globalIndex(), cache);
      auto lsmm_pairs_per_coll = lsmm_pairs->sliceByCached(o2::aod::dalitzee::emreducedeventId, collision.globalIndex(), cache);

      for (const auto& cut : fDalitzEECuts) {
        THashList* list_dalitzee_cut = static_cast<THashList*>(list_dalitzee->FindObject(cut.GetName()));
        THashList* list_track_cut = static_cast<THashList*>(list_track->FindObject(cut.GetName()));

        int nuls = 0, nlspp = 0, nlsmm = 0;
        for (auto& uls_pair : uls_pairs_per_coll) {
          auto pos = uls_pair.posTrack_as<aod::EMPrimaryTracks>();
          auto ele = uls_pair.negTrack_as<aod::EMPrimaryTracks>();
          if (cut.IsSelected<aod::EMPrimaryTracks>(uls_pair)) {
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

        for (auto& lspp_pair : lspp_pairs_per_coll) {
          if (cut.IsSelected<aod::EMPrimaryTracks>(lspp_pair)) {
            values[0] = lspp_pair.mee();
            values[1] = lspp_pair.pt();
            values[2] = lspp_pair.dcaeeXY();
            values[3] = lspp_pair.phiv();
            reinterpret_cast<THnSparseF*>(list_dalitzee_cut->FindObject("hs_dilepton_lspp"))->Fill(values);
            nlspp++;
          }
        } // end of lspp pair loop
        reinterpret_cast<TH1F*>(list_dalitzee_cut->FindObject("hNpair_lspp"))->Fill(nlspp);

        for (auto& lsmm_pair : lsmm_pairs_per_coll) {
          if (cut.IsSelected<aod::EMPrimaryTracks>(lsmm_pair)) {
            values[0] = lsmm_pair.mee();
            values[1] = lsmm_pair.pt();
            values[2] = lsmm_pair.dcaeeXY();
            values[3] = lsmm_pair.phiv();
            reinterpret_cast<THnSparseF*>(list_dalitzee_cut->FindObject("hs_dilepton_lsmm"))->Fill(values);
            nlsmm++;
          }
        } // end of lsmm pair loop
        reinterpret_cast<TH1F*>(list_dalitzee_cut->FindObject("hNpair_lsmm"))->Fill(nlsmm);

      } // end of cut loop
    }   // end of collision loop
  }     // end of process
  PROCESS_SWITCH(DalitzEEQC, processQC, "run Dalitz QC", true);

  void processDummy(aod::EMReducedEvents::iterator const& collision) {}
  PROCESS_SWITCH(DalitzEEQC, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DalitzEEQC>(cfgc, TaskName{"dalitz-ee-qc"})};
}
