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
using namespace o2::aod::pwgem::photon;
using std::array;

using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent, aod::EMEventsBz, aod::EMEventsNee>;
using MyCollision = MyCollisions::iterator;

using MyDalitzEEs = soa::Join<aod::DalitzEEs, aod::DalitzEEEMEventIds>;
using MyDalitzEE = MyDalitzEEs::iterator;

using MyTracks = soa::Join<aod::EMPrimaryElectrons, aod::EMPrimaryElectronEMEventIds, aod::EMPrimaryElectronsPrefilterBit>;
using MyTrack = MyTracks::iterator;

struct DalitzEEQC {

  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};

  Configurable<std::string> fConfigDalitzEECuts{"cfgDalitzEECuts", "mee_all_tpchadrejortofreq_lowB,nocut", "Comma separated list of dalitz ee cuts"};
  std::vector<DalitzEECut> fDalitzEECuts;

  Configurable<bool> cfgDoDCAstudy{"cfgDoDCAstudy", false, "flag to fill histograms for DCA"};

  Configurable<std::string> fConfigEMEventCut{"cfgEMEventCut", "minbias", "em event cut"}; // only 1 event cut per wagon
  EMEventCut fEMEventCut;
  static constexpr std::string_view event_types[2] = {"before", "after"};

  static constexpr std::string_view cent_det_names[3] = {"FT0M", "FT0A", "FT0C"};
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
    THashList* list_ev = reinterpret_cast<THashList*>(o2::aod::pwgem::photon::histogram::AddHistClass(fMainList, "Event"));
    for (const auto& evtype : event_types) {
      THashList* list_ev_type = reinterpret_cast<THashList*>(o2::aod::pwgem::photon::histogram::AddHistClass(list_ev, evtype.data()));
      o2::aod::pwgem::photon::histogram::DefineHistograms(list_ev_type, "Event", evtype.data());
    }

    THashList* list_track = reinterpret_cast<THashList*>(o2::aod::pwgem::photon::histogram::AddHistClass(fMainList, "Track"));
    THashList* list_dalitzee = reinterpret_cast<THashList*>(o2::aod::pwgem::photon::histogram::AddHistClass(fMainList, "DalitzEE"));

    for (const auto& cut : fDalitzEECuts) {
      const char* cutname = cut.GetName();
      THashList* list_track_cut = reinterpret_cast<THashList*>(o2::aod::pwgem::photon::histogram::AddHistClass(list_track, cutname));
      o2::aod::pwgem::photon::histogram::DefineHistograms(list_track_cut, "Track");

      THashList* list_dalitzee_cut = reinterpret_cast<THashList*>(o2::aod::pwgem::photon::histogram::AddHistClass(list_dalitzee, cutname));
      std::string histo_sub_group = "";
      if (doMix) {
        histo_sub_group += "mix";
      }
      if (cfgDoDCAstudy) {
        histo_sub_group += "dca";
      }
      o2::aod::pwgem::photon::histogram::DefineHistograms(list_dalitzee_cut, "DalitzEE", histo_sub_group.data());
    } // end of cut loop
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
    TString ev_cut_name = fConfigEMEventCut.value;
    fEMEventCut = *eventcuts::GetCut(ev_cut_name.Data());

    fOutputEvent.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("Event")));
    fOutputTrack.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("Track")));
    fOutputDalitzEE.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("DalitzEE")));
  }

  Partition<MyDalitzEEs> uls_pairs = o2::aod::dalitzee::sign == 0;
  Partition<MyDalitzEEs> lspp_pairs = o2::aod::dalitzee::sign == +1;
  Partition<MyDalitzEEs> lsmm_pairs = o2::aod::dalitzee::sign == -1;

  SliceCache cache;
  Preslice<MyDalitzEEs> perCollision = aod::dalitzee::emeventId;
  Preslice<MyTracks> perCollision_track = aod::emprimaryelectron::emeventId;
  Partition<MyCollisions> grouped_collisions = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax); // this goes to same event.

  void processQC(MyCollisions const&, MyDalitzEEs const&, MyTracks const&)
  {
    THashList* list_ev_before = static_cast<THashList*>(fMainList->FindObject("Event")->FindObject(event_types[0].data()));
    THashList* list_ev_after = static_cast<THashList*>(fMainList->FindObject("Event")->FindObject(event_types[1].data()));
    THashList* list_dalitzee = static_cast<THashList*>(fMainList->FindObject("DalitzEE"));
    THashList* list_track = static_cast<THashList*>(fMainList->FindObject("Track"));
    double values[3] = {0, 0, 0};
    double values_single[4] = {0, 0, 0, 0};
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

      auto uls_pairs_per_coll = uls_pairs->sliceByCached(o2::aod::dalitzee::emeventId, collision.globalIndex(), cache);
      auto lspp_pairs_per_coll = lspp_pairs->sliceByCached(o2::aod::dalitzee::emeventId, collision.globalIndex(), cache);
      auto lsmm_pairs_per_coll = lsmm_pairs->sliceByCached(o2::aod::dalitzee::emeventId, collision.globalIndex(), cache);

      for (const auto& cut : fDalitzEECuts) {
        THashList* list_dalitzee_cut = static_cast<THashList*>(list_dalitzee->FindObject(cut.GetName()));
        THashList* list_track_cut = static_cast<THashList*>(list_track->FindObject(cut.GetName()));

        std::vector<uint64_t> used_trackIds;
        used_trackIds.reserve(uls_pairs_per_coll.size() * 2);

        int nuls = 0, nlspp = 0, nlsmm = 0;
        for (auto& uls_pair : uls_pairs_per_coll) {
          auto pos = uls_pair.template posTrack_as<MyTracks>();
          auto ele = uls_pair.template negTrack_as<MyTracks>();
          if (cut.IsSelected<MyTracks>(uls_pair)) {
            dca_pos_3d = pos.dca3DinSigma();
            dca_ele_3d = ele.dca3DinSigma();
            dca_ee_3d = std::sqrt((dca_pos_3d * dca_pos_3d + dca_ele_3d * dca_ele_3d) / 2.);

            if (cfgDoDCAstudy) {
              values_single[0] = uls_pair.mass();
              values_single[1] = dca_pos_3d;
              values_single[2] = dca_ele_3d;
              values_single[3] = dca_ee_3d;
              reinterpret_cast<THnSparseF*>(list_dalitzee_cut->FindObject("hs_dilepton_uls_dca_same"))->Fill(values_single);
            }

            values[0] = uls_pair.mass();
            values[1] = uls_pair.pt();
            values[2] = dca_ee_3d;
            reinterpret_cast<THnSparseF*>(list_dalitzee_cut->FindObject("hs_dilepton_uls_same"))->Fill(values);
            reinterpret_cast<TH2F*>(list_dalitzee_cut->FindObject("hMvsPhiV_uls_same"))->Fill(uls_pair.phiv(), uls_pair.mass());
            nuls++;
            for (auto& track : {pos, ele}) {
              if (std::find(used_trackIds.begin(), used_trackIds.end(), track.globalIndex()) == used_trackIds.end()) {
                o2::aod::pwgem::photon::histogram::FillHistClass<EMHistType::kTrack>(list_track_cut, "", track);
                used_trackIds.emplace_back(track.globalIndex());
              }
            }
          }
        } // end of uls pair loop
        reinterpret_cast<TH1F*>(list_dalitzee_cut->FindObject("hNpair_uls"))->Fill(nuls);

        for (auto& lspp_pair : lspp_pairs_per_coll) {
          auto pos = lspp_pair.template posTrack_as<MyTracks>();
          auto ele = lspp_pair.template negTrack_as<MyTracks>();
          if (cut.IsSelected<MyTracks>(lspp_pair)) {
            dca_pos_3d = pos.dca3DinSigma();
            dca_ele_3d = ele.dca3DinSigma();
            dca_ee_3d = std::sqrt((dca_pos_3d * dca_pos_3d + dca_ele_3d * dca_ele_3d) / 2.);

            if (cfgDoDCAstudy) {
              values_single[0] = lspp_pair.mass();
              values_single[1] = dca_pos_3d;
              values_single[2] = dca_ele_3d;
              values_single[3] = dca_ee_3d;
              reinterpret_cast<THnSparseF*>(list_dalitzee_cut->FindObject("hs_dilepton_lspp_dca_same"))->Fill(values_single);
            }
            values[0] = lspp_pair.mass();
            values[1] = lspp_pair.pt();
            values[2] = dca_ee_3d;
            reinterpret_cast<THnSparseF*>(list_dalitzee_cut->FindObject("hs_dilepton_lspp_same"))->Fill(values);
            reinterpret_cast<TH2F*>(list_dalitzee_cut->FindObject("hMvsPhiV_lspp_same"))->Fill(lspp_pair.phiv(), lspp_pair.mass());
            nlspp++;
          }
        } // end of lspp pair loop
        reinterpret_cast<TH1F*>(list_dalitzee_cut->FindObject("hNpair_lspp"))->Fill(nlspp);

        for (auto& lsmm_pair : lsmm_pairs_per_coll) {
          auto pos = lsmm_pair.template posTrack_as<MyTracks>();
          auto ele = lsmm_pair.template negTrack_as<MyTracks>();
          if (cut.IsSelected<MyTracks>(lsmm_pair)) {
            dca_pos_3d = pos.dca3DinSigma();
            dca_ele_3d = ele.dca3DinSigma();
            dca_ee_3d = std::sqrt((dca_pos_3d * dca_pos_3d + dca_ele_3d * dca_ele_3d) / 2.);

            if (cfgDoDCAstudy) {
              values_single[0] = lsmm_pair.mass();
              values_single[1] = dca_pos_3d;
              values_single[2] = dca_ele_3d;
              values_single[3] = dca_ee_3d;
              reinterpret_cast<THnSparseF*>(list_dalitzee_cut->FindObject("hs_dilepton_lsmm_dca_same"))->Fill(values_single);
            }

            values[0] = lsmm_pair.mass();
            values[1] = lsmm_pair.pt();
            values[2] = dca_ee_3d;
            reinterpret_cast<THnSparseF*>(list_dalitzee_cut->FindObject("hs_dilepton_lsmm_same"))->Fill(values);
            reinterpret_cast<TH2F*>(list_dalitzee_cut->FindObject("hMvsPhiV_lsmm_same"))->Fill(lsmm_pair.phiv(), lsmm_pair.mass());
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
  ConfigurableAxis ConfCentBins{"ConfCentBins", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 999.f}, "Mixing bins - centrality"};
  using BinningType_M = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;
  using BinningType_A = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0A>;
  using BinningType_C = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
  BinningType_M colBinning_M{{ConfVtxBins, ConfCentBins}, true};
  BinningType_A colBinning_A{{ConfVtxBins, ConfCentBins}, true};
  BinningType_C colBinning_C{{ConfVtxBins, ConfCentBins}, true};

  Filter collisionFilter_common = nabs(o2::aod::collision::posZ) < 10.f && o2::aod::collision::numContrib > (uint16_t)0 && o2::aod::evsel::sel8 == true;
  Filter collisionFilter_centrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  Filter collisionFilter_subsys = (o2::aod::emevent::neeuls >= 1) || (o2::aod::emevent::neelspp >= 1) || (o2::aod::emevent::neelsmm >= 1);
  using MyFilteredCollisions = soa::Filtered<MyCollisions>; // this goes to mixed event.

  // e+, e- enter to event mixing, only if any pair exists. If you want to do mixed event, please store LS for ee
  template <typename TEvents, typename TTracks, typename TMixedBinning>
  void MixedEventPairing(TEvents const& collisions, TTracks const& tracks, TMixedBinning const& colBinning)
  {
    THashList* list_dalitzee = static_cast<THashList*>(fMainList->FindObject("DalitzEE"));
    double values[3] = {0, 0, 0};
    ROOT::Math::PtEtaPhiMVector v1, v2, v12;
    float phiv = 0;
    double values_single[4] = {0, 0, 0, 0};
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
      // LOGF(info, "collision1.globalIndex() = %d, collision1: posZ = %f, sel8 = %d, centFT0C = %f, ndl1 = %d | collision2.globalIndex() = %d, collision2: posZ = %f, sel8 = %d, centFT0C = %f, ndl2 = %d",collision1.globalIndex(), collision1.posZ(), collision1.sel8(), collision1.centFT0C(), tracks_coll1.size(), collision2.globalIndex(), collision2.posZ(), collision2.sel8(), collision2.centFT0C(), tracks_coll2.size());

      for (auto& cut : fDalitzEECuts) {
        THashList* list_dalitzee_cut = static_cast<THashList*>(list_dalitzee->FindObject(cut.GetName()));
        for (auto& [t1, t2] : combinations(soa::CombinationsFullIndexPolicy(tracks_coll1, tracks_coll2))) {
          v1 = ROOT::Math::PtEtaPhiMVector(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassElectron);
          v2 = ROOT::Math::PtEtaPhiMVector(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassElectron);
          v12 = v1 + v2;
          phiv = getPhivPair(t1.px(), t1.py(), t1.pz(), t2.px(), t2.py(), t2.pz(), t1.sign(), t2.sign(), collision1.bz());
          // if (!std::isfinite(phiv)) {
          //   LOGF(info, "t1.px() = %f, t1.py() = %f, t1.pz() = %f, t2.px() = %f, t2.py() = %f, t2.pz() = %f, phiv = %f", t1.px(), t1.py(), t1.pz(), t2.px(), t2.py(), t2.pz(), phiv);
          // }

          if (t1.trackId() == t2.trackId()) { // this is protection against pairing identical 2 tracks. This happens, when TTCA is used. TTCA can assign a track to several possible collisions.
            continue;
          }

          dca_pos_3d = t1.dca3DinSigma();
          dca_ele_3d = t2.dca3DinSigma();
          dca_ee_3d = std::sqrt((dca_pos_3d * dca_pos_3d + dca_ele_3d * dca_ele_3d) / 2.);

          values_single[0] = v12.M();
          values_single[1] = dca_pos_3d;
          values_single[2] = dca_ele_3d;
          values_single[3] = dca_ee_3d;

          values[0] = v12.M();
          values[1] = v12.Pt();
          values[2] = dca_ee_3d;

          if (cut.IsSelectedTrack(t1) && cut.IsSelectedTrack(t2) && cut.IsSelectedPair(v12.M(), dca_ee_3d, phiv)) {
            if (t1.sign() * t2.sign() < 0) {
              reinterpret_cast<THnSparseF*>(list_dalitzee_cut->FindObject("hs_dilepton_uls_mix"))->Fill(values);
              reinterpret_cast<TH2F*>(list_dalitzee_cut->FindObject("hMvsPhiV_uls_mix"))->Fill(phiv, v12.M());
            } else if (t1.sign() > 0 && t2.sign() > 0) {
              reinterpret_cast<THnSparseF*>(list_dalitzee_cut->FindObject("hs_dilepton_lspp_mix"))->Fill(values);
              reinterpret_cast<TH2F*>(list_dalitzee_cut->FindObject("hMvsPhiV_lspp_mix"))->Fill(phiv, v12.M());
            } else if (t1.sign() < 0 && t2.sign() < 0) {
              reinterpret_cast<THnSparseF*>(list_dalitzee_cut->FindObject("hs_dilepton_lsmm_mix"))->Fill(values);
              reinterpret_cast<TH2F*>(list_dalitzee_cut->FindObject("hMvsPhiV_lsmm_mix"))->Fill(phiv, v12.M());
            } else {
              LOGF(info, "This should not happen.");
            }

            if (cfgDoDCAstudy) {
              if (t1.sign() * t2.sign() < 0) {
                reinterpret_cast<THnSparseF*>(list_dalitzee_cut->FindObject("hs_dilepton_uls_dca_mix"))->Fill(values_single);
              } else if (t1.sign() > 0 && t2.sign() > 0) {
                reinterpret_cast<THnSparseF*>(list_dalitzee_cut->FindObject("hs_dilepton_lspp_dca_mix"))->Fill(values_single);
              } else if (t1.sign() < 0 && t2.sign() < 0) {
                reinterpret_cast<THnSparseF*>(list_dalitzee_cut->FindObject("hs_dilepton_lsmm_dca_mix"))->Fill(values_single);
              } else {
                LOGF(info, "This should not happen.");
              }
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
  PROCESS_SWITCH(DalitzEEQC, processEventMixing, "run Dalitz EE QC event mixing", true);

  void processDummy(MyCollisions const&) {}
  PROCESS_SWITCH(DalitzEEQC, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DalitzEEQC>(cfgc, TaskName{"dalitz-ee-qc"})};
}
