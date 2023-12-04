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
#include "PWGEM/PhotonMeson/Utils/PCMUtilities.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using std::array;

using MyCollisions = soa::Join<aod::EMReducedEvents, aod::EMReducedEventsMult, aod::EMReducedEventsCent, aod::EMReducedEventsBz, aod::EMReducedEventsNee>;
using MyCollision = MyCollisions::iterator;

using MyDalitzEEs = soa::Join<aod::DalitzEEs, aod::DalitzEEEMReducedEventIds>;
using MyDalitzEE = MyDalitzEEs::iterator;

using MyTracks = soa::Join<aod::EMPrimaryElectrons, aod::EMPrimaryElectronEMReducedEventIds, aod::EMPrimaryElectronsPrefilterBit>;
using MyTrack = MyTracks::iterator;

struct DalitzEEQC {
  Configurable<std::string> fConfigDalitzEECuts{"cfgDalitzEECuts", "mee_all_tpchadrejortofreq_lowB,nocut", "Comma separated list of dalitz ee cuts"};
  std::vector<DalitzEECut> fDalitzEECuts;

  OutputObj<THashList> fOutputEvent{"Event"};
  OutputObj<THashList> fOutputTrack{"Track"};
  OutputObj<THashList> fOutputDalitzEE{"DalitzEE"};
  THashList* fMainList = new THashList();
  bool doMix = false;

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
      if (doMix) {
        o2::aod::emphotonhistograms::DefineHistograms(list, "DalitzEE", "mix");
      } else {
        o2::aod::emphotonhistograms::DefineHistograms(list, "DalitzEE", "");
      }
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
    if (context.mOptions.get<bool>("processEventMixing")) {
      doMix = true;
    }

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

  std::vector<uint64_t> used_trackIds;

  void processQC(MyCollisions const& collisions, MyDalitzEEs const& dileptons, MyTracks const& tracks)
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
      auto lspp_pairs_per_coll = lspp_pairs->sliceByCached(o2::aod::dalitzee::emreducedeventId, collision.globalIndex(), cache);
      auto lsmm_pairs_per_coll = lsmm_pairs->sliceByCached(o2::aod::dalitzee::emreducedeventId, collision.globalIndex(), cache);

      for (const auto& cut : fDalitzEECuts) {
        THashList* list_dalitzee_cut = static_cast<THashList*>(list_dalitzee->FindObject(cut.GetName()));
        THashList* list_track_cut = static_cast<THashList*>(list_track->FindObject(cut.GetName()));
        used_trackIds.reserve(uls_pairs_per_coll.size() * 2);

        int nuls = 0, nlspp = 0, nlsmm = 0;
        for (auto& uls_pair : uls_pairs_per_coll) {
          auto pos = uls_pair.template posTrack_as<MyTracks>();
          auto ele = uls_pair.template negTrack_as<MyTracks>();
          if (cut.IsSelected<MyTracks>(uls_pair)) {
            values[0] = uls_pair.mass();
            values[1] = uls_pair.pt();
            values[2] = uls_pair.dcaXY();
            values[3] = uls_pair.phiv();
            reinterpret_cast<THnSparseF*>(list_dalitzee_cut->FindObject("hs_dilepton_uls_same"))->Fill(values);
            nuls++;
            for (auto& track : {pos, ele}) {
              if (std::find(used_trackIds.begin(), used_trackIds.end(), track.globalIndex()) == used_trackIds.end()) {
                o2::aod::emphotonhistograms::FillHistClass<EMHistType::kTrack>(list_track_cut, "", track);
                used_trackIds.emplace_back(track.globalIndex());
              }
            }
          }
        } // end of uls pair loop
        reinterpret_cast<TH1F*>(list_dalitzee_cut->FindObject("hNpair_uls"))->Fill(nuls);

        for (auto& lspp_pair : lspp_pairs_per_coll) {
          if (cut.IsSelected<MyTracks>(lspp_pair)) {
            values[0] = lspp_pair.mass();
            values[1] = lspp_pair.pt();
            values[2] = lspp_pair.dcaXY();
            values[3] = lspp_pair.phiv();
            reinterpret_cast<THnSparseF*>(list_dalitzee_cut->FindObject("hs_dilepton_lspp_same"))->Fill(values);
            nlspp++;
          }
        } // end of lspp pair loop
        reinterpret_cast<TH1F*>(list_dalitzee_cut->FindObject("hNpair_lspp"))->Fill(nlspp);

        for (auto& lsmm_pair : lsmm_pairs_per_coll) {
          if (cut.IsSelected<MyTracks>(lsmm_pair)) {
            values[0] = lsmm_pair.mass();
            values[1] = lsmm_pair.pt();
            values[2] = lsmm_pair.dcaXY();
            values[3] = lsmm_pair.phiv();
            reinterpret_cast<THnSparseF*>(list_dalitzee_cut->FindObject("hs_dilepton_lsmm_same"))->Fill(values);
            nlsmm++;
          }
        } // end of lsmm pair loop
        reinterpret_cast<TH1F*>(list_dalitzee_cut->FindObject("hNpair_lsmm"))->Fill(nlsmm);

        used_trackIds.clear();
        used_trackIds.shrink_to_fit();
      } // end of cut loop
    }   // end of collision loop
  }     // end of process
  PROCESS_SWITCH(DalitzEEQC, processQC, "run Dalitz QC", true);

  Configurable<int> ndepth{"ndepth", 10, "depth for event mixing"};
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis ConfCentBins{"ConfCentBins", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 7.0f, 80.0f, 90.0f, 100.0f, 999.f}, "Mixing bins - centrality"};
  using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;
  BinningType colBinning{{ConfVtxBins, ConfCentBins}, true};
  Filter collisionFilter_common = nabs(o2::aod::collision::posZ) < 10.f && o2::aod::collision::numContrib > (uint16_t)0 && o2::aod::evsel::sel8 == true;
  Filter collisionFilter_subsys = (o2::aod::emreducedevent::neeuls >= 1) || (o2::aod::emreducedevent::neelspp >= 1) || (o2::aod::emreducedevent::neelsmm >= 1);
  using MyFilteredCollisions = soa::Filtered<MyCollisions>; // this goes to mixed event.

  // e+, e- enter to event mixing, only if any pair exists. If you want to do mixed event, please store LS for ee
  void processEventMixing(MyFilteredCollisions const& collisions, MyDalitzEEs const& dileptons, MyTracks const& tracks)
  {
    THashList* list_dalitzee = static_cast<THashList*>(fMainList->FindObject("DalitzEE"));
    double values[4] = {0, 0, 0, 0};
    ROOT::Math::PtEtaPhiMVector v1, v2, v12;
    float phiv = 0;

    for (auto& [collision1, collision2] : soa::selfCombinations(colBinning, ndepth, -1, collisions, collisions)) { // internally, CombinationsStrictlyUpperIndexPolicy(collisions, collisions) is called.
      auto dileptons_coll1 = dileptons.sliceBy(perCollision, collision1.globalIndex());
      auto dileptons_coll2 = dileptons.sliceBy(perCollision, collision2.globalIndex());
      // LOGF(info, "collision1.globalIndex() = %d, collision1: posZ = %f, sel8 = %d | collision2.globalIndex() = %d, collision2: posZ = %f, sel8 = %d",collision1.globalIndex(), collision1.posZ(), collision1.sel8(), collision2.globalIndex(), collision2.posZ(), collision2.sel8());

      for (auto& cut : fDalitzEECuts) {
        THashList* list_dalitzee_cut = static_cast<THashList*>(list_dalitzee->FindObject(cut.GetName()));
        for (auto& [dl1, dl2] : combinations(soa::CombinationsFullIndexPolicy(dileptons_coll1, dileptons_coll2))) {
          if (!cut.IsSelected<MyTracks>(dl1) || !cut.IsSelected<MyTracks>(dl2)) {
            continue;
          }

          auto pos1 = dl1.template posTrack_as<MyTracks>();
          auto ele1 = dl1.template negTrack_as<MyTracks>();
          auto pos2 = dl2.template posTrack_as<MyTracks>();
          auto ele2 = dl2.template negTrack_as<MyTracks>();

          // float dcaxy1 = t1.dcaXY() / sqrt(t1.cYY());
          // float dcaxy2 = t2.dcaXY() / sqrt(t2.cYY());
          // float dcamumuxy = sqrt((pow(dcaxy1, 2) + pow(dcaxy2, 2)) / 2.);
          // float dcaz1 = t1.dcaZ() / sqrt(t1.cZZ());
          // float dcaz2 = t2.dcaZ() / sqrt(t2.cZZ());
          // float dcamumuz = sqrt((pow(dcaz1, 2) + pow(dcaz2, 2)) / 2.);

          // mix uls1
          // ROOT::Math::PtEtaPhiMVector v1(pos1.pt(), pos1.eta(), pos1.phi(), o2::constants::physics::MassElectron);
          // ROOT::Math::PtEtaPhiMVector v2(ele2.pt(), ele2.eta(), ele2.phi(), o2::constants::physics::MassElectron);
          // ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          v1 = ROOT::Math::PtEtaPhiMVector(pos1.pt(), pos1.eta(), pos1.phi(), o2::constants::physics::MassElectron);
          v2 = ROOT::Math::PtEtaPhiMVector(ele2.pt(), ele2.eta(), ele2.phi(), o2::constants::physics::MassElectron);
          v12 = v1 + v2;
          phiv = getPhivPair(pos1.px(), pos1.py(), pos1.pz(), ele2.px(), ele2.py(), ele2.pz(), pos1.sign(), ele2.sign(), collision1.bz());
          values[0] = v12.M();
          values[1] = v12.Pt();
          values[2] = sqrt((pow(pos1.dcaXY() / sqrt(pos1.cYY()), 2) + pow(ele2.dcaXY() / sqrt(ele2.cYY()), 2)) / 2.); // pair DCAxy
          values[3] = phiv;
          if (cut.IsSelectedTrack(pos1) && cut.IsSelectedTrack(ele2) && cut.IsSelectedPair(v12.M(), phiv)) {
            reinterpret_cast<THnSparseF*>(list_dalitzee_cut->FindObject("hs_dilepton_uls_mix"))->Fill(values);
          }

          // mix uls2
          v1 = ROOT::Math::PtEtaPhiMVector(pos2.pt(), pos2.eta(), pos2.phi(), o2::constants::physics::MassElectron);
          v2 = ROOT::Math::PtEtaPhiMVector(ele1.pt(), ele1.eta(), ele1.phi(), o2::constants::physics::MassElectron);
          v12 = v1 + v2;
          phiv = getPhivPair(pos2.px(), pos2.py(), pos2.pz(), ele1.px(), ele1.py(), ele1.pz(), pos2.sign(), ele1.sign(), collision1.bz());
          values[0] = v12.M();
          values[1] = v12.Pt();
          values[2] = sqrt((pow(pos2.dcaXY() / sqrt(pos2.cYY()), 2) + pow(ele1.dcaXY() / sqrt(ele1.cYY()), 2)) / 2.); // pair DCAxy
          values[3] = phiv;
          if (cut.IsSelectedTrack(pos2) && cut.IsSelectedTrack(ele1) && cut.IsSelectedPair(v12.M(), phiv)) {
            reinterpret_cast<THnSparseF*>(list_dalitzee_cut->FindObject("hs_dilepton_uls_mix"))->Fill(values);
          }

          // mix lspp
          v1 = ROOT::Math::PtEtaPhiMVector(pos1.pt(), pos1.eta(), pos1.phi(), o2::constants::physics::MassElectron);
          v2 = ROOT::Math::PtEtaPhiMVector(pos2.pt(), pos2.eta(), pos2.phi(), o2::constants::physics::MassElectron);
          v12 = v1 + v2;
          phiv = getPhivPair(pos1.px(), pos1.py(), pos1.pz(), pos2.px(), pos2.py(), pos2.pz(), pos1.sign(), pos2.sign(), collision1.bz());
          values[0] = v12.M();
          values[1] = v12.Pt();
          values[2] = sqrt((pow(pos1.dcaXY() / sqrt(pos1.cYY()), 2) + pow(pos2.dcaXY() / sqrt(pos2.cYY()), 2)) / 2.); // pair DCAxy
          values[3] = phiv;
          if (cut.IsSelectedTrack(pos1) && cut.IsSelectedTrack(pos2) && cut.IsSelectedPair(v12.M(), phiv)) {
            reinterpret_cast<THnSparseF*>(list_dalitzee_cut->FindObject("hs_dilepton_lspp_mix"))->Fill(values);
          }

          // mix lsmm
          v1 = ROOT::Math::PtEtaPhiMVector(ele1.pt(), ele1.eta(), ele1.phi(), o2::constants::physics::MassElectron);
          v2 = ROOT::Math::PtEtaPhiMVector(ele2.pt(), ele2.eta(), ele2.phi(), o2::constants::physics::MassElectron);
          v12 = v1 + v2;
          phiv = getPhivPair(ele1.px(), ele1.py(), ele1.pz(), ele2.px(), ele2.py(), ele2.pz(), ele1.sign(), ele2.sign(), collision1.bz());
          values[0] = v12.M();
          values[1] = v12.Pt();
          values[2] = sqrt((pow(ele1.dcaXY() / sqrt(ele1.cYY()), 2) + pow(ele2.dcaXY() / sqrt(ele2.cYY()), 2)) / 2.); // pair DCAxy
          values[3] = phiv;
          if (cut.IsSelectedTrack(ele1) && cut.IsSelectedTrack(ele2) && cut.IsSelectedPair(v12.M(), phiv)) {
            reinterpret_cast<THnSparseF*>(list_dalitzee_cut->FindObject("hs_dilepton_lsmm_mix"))->Fill(values);
          }

        } // end of different dileptn combinations
      }   // end of cut loop
    }     // end of different collision combinations
  }
  PROCESS_SWITCH(DalitzEEQC, processEventMixing, "run Dalitz EE QC event mixing", true);

  void processDummy(MyCollisions const& collisions) {}
  PROCESS_SWITCH(DalitzEEQC, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DalitzEEQC>(cfgc, TaskName{"dalitz-ee-qc"})};
}
