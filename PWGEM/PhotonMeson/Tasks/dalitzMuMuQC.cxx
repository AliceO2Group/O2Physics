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

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::photon;
using std::array;

using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent, aod::EMEventsNmumu>;
using MyCollision = MyCollisions::iterator;

using MyDalitzMuMus = soa::Join<aod::DalitzMuMus, aod::DalitzMuMuEMEventIds>;
using MyDalitzMuMu = MyDalitzMuMus::iterator;

using MyTracks = soa::Join<aod::EMPrimaryMuons, aod::EMPrimaryMuonEMEventIds, aod::EMPrimaryMuonsPrefilterBit>;
using MyTrack = MyTracks::iterator;

struct DalitzMuMuQC {
  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};

  Configurable<std::string> fConfigDalitzMuMuCuts{"cfgDalitzMuMuCuts", "nocut", "Comma separated list of dalitz mumu cuts"};
  std::vector<DalitzEECut> fDalitzMuMuCuts;

  Configurable<std::string> fConfigEMEventCut{"cfgEMEventCut", "minbias", "em event cut"}; // only 1 event cut per wagon
  EMEventCut fEMEventCut;
  static constexpr std::string_view event_types[2] = {"before", "after"};

  OutputObj<THashList> fOutputEvent{"Event"};
  OutputObj<THashList> fOutputTrack{"Track"};
  OutputObj<THashList> fOutputDalitzMuMu{"DalitzMuMu"};
  THashList* fMainList = new THashList();
  bool doMix = false;

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

    o2::aod::pwgem::photon::histogram::AddHistClass(fMainList, "Track");
    THashList* list_track = reinterpret_cast<THashList*>(fMainList->FindObject("Track"));

    o2::aod::pwgem::photon::histogram::AddHistClass(fMainList, "DalitzMuMu");
    THashList* list_dalitzmumu = reinterpret_cast<THashList*>(fMainList->FindObject("DalitzMuMu"));

    for (const auto& cut : fDalitzMuMuCuts) {
      const char* cutname = cut.GetName();
      o2::aod::pwgem::photon::histogram::AddHistClass(list_track, cutname);
      o2::aod::pwgem::photon::histogram::AddHistClass(list_dalitzmumu, cutname);
    }

    // for single tracks
    for (auto& cut : fDalitzMuMuCuts) {
      std::string_view cutname = cut.GetName();
      THashList* list = reinterpret_cast<THashList*>(fMainList->FindObject("Track")->FindObject(cutname.data()));
      o2::aod::pwgem::photon::histogram::DefineHistograms(list, "Track", "Mu");
    }

    // for DalitzMuMus
    for (auto& cut : fDalitzMuMuCuts) {
      std::string_view cutname = cut.GetName();
      THashList* list = reinterpret_cast<THashList*>(fMainList->FindObject("DalitzMuMu")->FindObject(cutname.data()));

      if (doMix) {
        o2::aod::pwgem::photon::histogram::DefineHistograms(list, "DalitzMuMu", "mix");
      } else {
        o2::aod::pwgem::photon::histogram::DefineHistograms(list, "DalitzMuMu", "");
      }
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
    if (context.mOptions.get<bool>("processEventMixing")) {
      doMix = true;
    }

    DefineCuts();
    addhistograms(); // please call this after DefinCuts();
    TString ev_cut_name = fConfigEMEventCut.value;
    fEMEventCut = *eventcuts::GetCut(ev_cut_name.Data());

    fOutputEvent.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("Event")));
    fOutputTrack.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("Track")));
    fOutputDalitzMuMu.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("DalitzMuMu")));
  }

  Partition<MyDalitzMuMus> uls_pairs = o2::aod::dalitzmumu::sign == 0;
  Partition<MyDalitzMuMus> lspp_pairs = o2::aod::dalitzmumu::sign == +1;
  Partition<MyDalitzMuMus> lsmm_pairs = o2::aod::dalitzmumu::sign == -1;

  SliceCache cache;
  Preslice<MyDalitzMuMus> perCollision = aod::dalitzmumu::emeventId;
  Preslice<MyTracks> perCollision_track = aod::emprimarymuon::emeventId;
  Partition<MyCollisions> grouped_collisions = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax); // this goes to same event.

  std::vector<uint64_t> used_trackIds;

  void processQC(MyCollisions const&, MyDalitzMuMus const&, MyTracks const&)
  {
    THashList* list_ev_before = static_cast<THashList*>(fMainList->FindObject("Event")->FindObject(event_types[0].data()));
    THashList* list_ev_after = static_cast<THashList*>(fMainList->FindObject("Event")->FindObject(event_types[1].data()));
    THashList* list_dalitzmumu = static_cast<THashList*>(fMainList->FindObject("DalitzMuMu"));
    THashList* list_track = static_cast<THashList*>(fMainList->FindObject("Track"));
    double values[3] = {0, 0, 0};
    float dca_pos_3d = 999.f, dca_ele_3d = 999.f, dca_ee_3d = 999.f;

    for (auto& collision : grouped_collisions) {
      float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
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

      auto uls_pairs_per_coll = uls_pairs->sliceByCached(o2::aod::dalitzmumu::emeventId, collision.globalIndex(), cache);
      auto lspp_pairs_per_coll = lspp_pairs->sliceByCached(o2::aod::dalitzmumu::emeventId, collision.globalIndex(), cache);
      auto lsmm_pairs_per_coll = lsmm_pairs->sliceByCached(o2::aod::dalitzmumu::emeventId, collision.globalIndex(), cache);

      for (const auto& cut : fDalitzMuMuCuts) {
        THashList* list_dalitzmumu_cut = static_cast<THashList*>(list_dalitzmumu->FindObject(cut.GetName()));
        THashList* list_track_cut = static_cast<THashList*>(list_track->FindObject(cut.GetName()));
        used_trackIds.reserve(uls_pairs_per_coll.size() * 2);

        int nuls = 0, nlspp = 0, nlsmm = 0;
        for (auto& uls_pair : uls_pairs_per_coll) {
          auto pos = uls_pair.template posTrack_as<MyTracks>();
          auto ele = uls_pair.template negTrack_as<MyTracks>();

          if (cut.IsSelected<MyTracks>(uls_pair)) {
            dca_pos_3d = pos.dca3DinSigma();
            dca_ele_3d = ele.dca3DinSigma();
            dca_ee_3d = std::sqrt((dca_pos_3d * dca_pos_3d + dca_ele_3d * dca_ele_3d) / 2.);

            values[0] = uls_pair.mass();
            values[1] = uls_pair.pt();
            values[2] = dca_ee_3d;
            reinterpret_cast<THnSparseF*>(list_dalitzmumu_cut->FindObject("hs_dilepton_uls_same"))->Fill(values);
            nuls++;
            for (auto& track : {pos, ele}) {
              if (std::find(used_trackIds.begin(), used_trackIds.end(), track.globalIndex()) == used_trackIds.end()) {
                o2::aod::pwgem::photon::histogram::FillHistClass<EMHistType::kTrack>(list_track_cut, "", track);
                used_trackIds.emplace_back(track.globalIndex());
              }
            }
          }
        } // end of uls pair loop
        reinterpret_cast<TH1F*>(list_dalitzmumu_cut->FindObject("hNpair_uls"))->Fill(nuls);

        for (auto& lspp_pair : lspp_pairs_per_coll) {
          auto pos = lspp_pair.template posTrack_as<MyTracks>();
          auto ele = lspp_pair.template negTrack_as<MyTracks>();
          if (cut.IsSelected<MyTracks>(lspp_pair)) {
            dca_pos_3d = pos.dca3DinSigma();
            dca_ele_3d = ele.dca3DinSigma();
            dca_ee_3d = std::sqrt((dca_pos_3d * dca_pos_3d + dca_ele_3d * dca_ele_3d) / 2.);

            values[0] = lspp_pair.mass();
            values[1] = lspp_pair.pt();
            values[2] = dca_ee_3d;
            reinterpret_cast<THnSparseF*>(list_dalitzmumu_cut->FindObject("hs_dilepton_lspp_same"))->Fill(values);
            nlspp++;
          }
        } // end of lspp pair loop
        reinterpret_cast<TH1F*>(list_dalitzmumu_cut->FindObject("hNpair_lspp"))->Fill(nlspp);

        for (auto& lsmm_pair : lsmm_pairs_per_coll) {
          auto pos = lsmm_pair.template posTrack_as<MyTracks>();
          auto ele = lsmm_pair.template negTrack_as<MyTracks>();
          if (cut.IsSelected<MyTracks>(lsmm_pair)) {
            dca_pos_3d = pos.dca3DinSigma();
            dca_ele_3d = ele.dca3DinSigma();
            dca_ee_3d = std::sqrt((dca_pos_3d * dca_pos_3d + dca_ele_3d * dca_ele_3d) / 2.);

            values[0] = lsmm_pair.mass();
            values[1] = lsmm_pair.pt();
            values[2] = dca_ee_3d;
            reinterpret_cast<THnSparseF*>(list_dalitzmumu_cut->FindObject("hs_dilepton_lsmm_same"))->Fill(values);
            nlsmm++;
          }
        } // end of lsmm pair loop
        reinterpret_cast<TH1F*>(list_dalitzmumu_cut->FindObject("hNpair_lsmm"))->Fill(nlsmm);

        used_trackIds.clear();
        used_trackIds.shrink_to_fit();
      } // end of cut loop
    }   // end of collision loop
  }     // end of process
  PROCESS_SWITCH(DalitzMuMuQC, processQC, "run Dalitz QC", true);

  Configurable<int> ndepth{"ndepth", 10, "depth for event mixing"};
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis ConfCentBins{"ConfCentBins", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 999.f}, "Mixing bins - centrality"};
  using BinningType_M = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;
  using BinningType_A = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0A>;
  using BinningType_C = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
  BinningType_M colBinning_M{{ConfVtxBins, ConfCentBins}, true};
  BinningType_A colBinning_A{{ConfVtxBins, ConfCentBins}, true};
  BinningType_C colBinning_C{{ConfVtxBins, ConfCentBins}, true};

  Filter collisionFilter_common = nabs(o2::aod::collision::posZ) < 10.f && o2::aod::collision::numContrib > (uint16_t)0 && o2::aod::evsel::sel8 == true;
  Filter collisionFilter_centrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  Filter collisionFilter_subsys = (o2::aod::emevent::nmumuuls >= 1) || (o2::aod::emevent::nmumulspp >= 1) || (o2::aod::emevent::nmumulsmm >= 1);
  using MyFilteredCollisions = soa::Filtered<MyCollisions>; // this goes to mixed event.

  // mu+, mu- enter to event mixing, only if any pair exists. If you want to do mixed event, please store LS for mumu
  template <typename TEvents, typename TTracks, typename TMixedBinning>
  void MixedEventPairing(TEvents const& collisions, TTracks const& tracks, TMixedBinning const& colBinning)
  {
    THashList* list_dalitzmumu = static_cast<THashList*>(fMainList->FindObject("DalitzMuMu"));
    double values[3] = {0, 0, 0};
    ROOT::Math::PtEtaPhiMVector v1, v2, v12;
    float phiv = 0;
    float dca_pos_3d = 999.f, dca_ele_3d = 999.f, dca_ee_3d = 999.f;

    for (auto& [collision1, collision2] : soa::selfCombinations(colBinning, ndepth, -1, collisions, collisions)) { // internally, CombinationsStrictlyUpperIndexPolicy(collisions, collisions) is called.

      const float centralities1[3] = {collision1.centFT0M(), collision1.centFT0A(), collision1.centFT0C()};
      const float centralities2[3] = {collision2.centFT0M(), collision2.centFT0A(), collision2.centFT0C()};

      if (centralities1[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities1[cfgCentEstimator]) {
        continue;
      }
      if (centralities2[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities2[cfgCentEstimator]) {
        continue;
      }
      if (!fEMEventCut.IsSelected(collision1) || !fEMEventCut.IsSelected(collision2)) {
        continue;
      }

      auto tracks_coll1 = tracks.sliceBy(perCollision_track, collision1.globalIndex());
      auto tracks_coll2 = tracks.sliceBy(perCollision_track, collision2.globalIndex());

      for (auto& cut : fDalitzMuMuCuts) {
        THashList* list_dalitzmumu_cut = static_cast<THashList*>(list_dalitzmumu->FindObject(cut.GetName()));
        for (auto& [t1, t2] : combinations(soa::CombinationsFullIndexPolicy(tracks_coll1, tracks_coll2))) {
          if (t1.trackId() == t2.trackId()) { // this is protection against pairing identical 2 tracks. This happens, when TTCA is used. TTCA can assign a track to several possible collisions.
            continue;
          }

          v1 = ROOT::Math::PtEtaPhiMVector(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassMuon);
          v2 = ROOT::Math::PtEtaPhiMVector(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassMuon);
          v12 = v1 + v2;

          dca_pos_3d = t1.dca3DinSigma();
          dca_ele_3d = t2.dca3DinSigma();
          dca_ee_3d = std::sqrt((dca_pos_3d * dca_pos_3d + dca_ele_3d * dca_ele_3d) / 2.);

          values[0] = v12.M();
          values[1] = v12.Pt();
          values[2] = dca_ee_3d;

          if (cut.IsSelectedTrack(t1) && cut.IsSelectedTrack(t2) && cut.IsSelectedPair(v12.M(), dca_ee_3d, phiv)) {
            if (t1.sign() * t2.sign() < 0) {
              reinterpret_cast<THnSparseF*>(list_dalitzmumu_cut->FindObject("hs_dilepton_uls_mix"))->Fill(values);
            } else if (t1.sign() > 0 && t2.sign() > 0) {
              reinterpret_cast<THnSparseF*>(list_dalitzmumu_cut->FindObject("hs_dilepton_lspp_mix"))->Fill(values);
            } else if (t1.sign() < 0 && t2.sign() < 0) {
              reinterpret_cast<THnSparseF*>(list_dalitzmumu_cut->FindObject("hs_dilepton_lsmm_mix"))->Fill(values);
            } else {
              LOGF(info, "This should not happen.");
            }
          }
        } // end of different dileptn combinations
      }   // end of cut loop
    }     // end of different collision combinations
  }

  void processEventMixing(MyFilteredCollisions const& collisions, MyTracks const& tracks)
  {
    if (cfgCentEstimator == 0) {
      MixedEventPairing(collisions, tracks, colBinning_M);
    } else if (cfgCentEstimator == 1) {
      MixedEventPairing(collisions, tracks, colBinning_A);
    } else if (cfgCentEstimator == 2) {
      MixedEventPairing(collisions, tracks, colBinning_C);
    }
  }
  PROCESS_SWITCH(DalitzMuMuQC, processEventMixing, "run Dalitz MuMu QC event mixing", true);

  void processDummy(MyCollisions const&) {}
  PROCESS_SWITCH(DalitzMuMuQC, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DalitzMuMuQC>(cfgc, TaskName{"dalitz-mumu-qc"})};
}
