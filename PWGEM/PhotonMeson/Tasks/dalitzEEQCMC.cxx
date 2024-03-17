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
using namespace o2::aod::pwgem::mcutil;
using namespace o2::aod::pwgem::photon;
using std::array;

using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent, aod::EMEventsNee, aod::EMMCEventLabels>;
using MyCollision = MyCollisions::iterator;

using MyDalitzEEs = soa::Join<aod::DalitzEEs, aod::DalitzEEEMEventIds>;
using MyDalitzEE = MyDalitzEEs::iterator;

using MyTracks = soa::Join<aod::EMPrimaryElectrons, aod::EMPrimaryElectronEMEventIds, aod::EMPrimaryElectronsPrefilterBit>;
using MyMCTracks = soa::Join<MyTracks, aod::EMPrimaryElectronMCLabels>;

struct DalitzEEQCMC {
  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};

  Configurable<std::string> fConfigDalitzEECuts{"cfgDalitzEECuts", "mee_all_tpchadrejortofreq_lowB,nocut", "Comma separated list of dalitz ee cuts"};
  std::vector<DalitzEECut> fDalitzEECuts;

  Configurable<bool> cfgDoDCAstudy{"cfgDoDCAstudy", false, "flag to fill histograms for DCA"};

  Configurable<std::string> fConfigEMEventCut{"cfgEMEventCut", "minbias", "em event cut"}; // only 1 event cut per wagon
  EMEventCut fEMEventCut;
  static constexpr std::string_view event_types[2] = {"before", "after"};

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
    o2::aod::pwgem::photon::histogram::AddHistClass(fMainList, "Event");
    THashList* list_ev = reinterpret_cast<THashList*>(fMainList->FindObject("Event"));
    for (const auto& evtype : event_types) {
      THashList* list_ev_type = reinterpret_cast<THashList*>(o2::aod::pwgem::photon::histogram::AddHistClass(list_ev, evtype.data()));
      o2::aod::pwgem::photon::histogram::DefineHistograms(list_ev_type, "Event", evtype.data());
    }

    o2::aod::pwgem::photon::histogram::AddHistClass(fMainList, "Track");
    THashList* list_track = reinterpret_cast<THashList*>(fMainList->FindObject("Track"));

    o2::aod::pwgem::photon::histogram::AddHistClass(fMainList, "DalitzEE");
    THashList* list_dalitzee = reinterpret_cast<THashList*>(fMainList->FindObject("DalitzEE"));

    o2::aod::pwgem::photon::histogram::AddHistClass(fMainList, "Generated");
    THashList* list_gen = reinterpret_cast<THashList*>(fMainList->FindObject("Generated"));
    o2::aod::pwgem::photon::histogram::DefineHistograms(list_gen, "Generated", "dielectron");

    for (const auto& cut : fDalitzEECuts) {
      const char* cutname = cut.GetName();
      o2::aod::pwgem::photon::histogram::AddHistClass(list_track, cutname);
      o2::aod::pwgem::photon::histogram::AddHistClass(list_dalitzee, cutname);
    }

    // for single tracks
    for (auto& cut : fDalitzEECuts) {
      std::string_view cutname = cut.GetName();
      THashList* list = reinterpret_cast<THashList*>(fMainList->FindObject("Track")->FindObject(cutname.data()));
      o2::aod::pwgem::photon::histogram::DefineHistograms(list, "Track", "mc");
    }

    // for DalitzEEs
    for (auto& cut : fDalitzEECuts) {
      std::string_view cutname = cut.GetName();
      THashList* list = reinterpret_cast<THashList*>(fMainList->FindObject("DalitzEE")->FindObject(cutname.data()));
      // o2::aod::pwgem::photon::histogram::DefineHistograms(list, "DalitzEE", "mc");
      std::string histo_sub_group = "mc,";
      if (cfgDoDCAstudy) {
        histo_sub_group += "dca";
      }
      o2::aod::pwgem::photon::histogram::DefineHistograms(list, "DalitzEE", histo_sub_group.data());
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
    TString ev_cut_name = fConfigEMEventCut.value;
    fEMEventCut = *eventcuts::GetCut(ev_cut_name.Data());

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
  Preslice<MyDalitzEEs> perCollision = aod::dalitzee::emeventId;
  Partition<MyCollisions> grouped_collisions = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax); // this goes to same event.

  std::vector<uint64_t> used_trackIds;

  void processQCMC(MyCollisions const& collisions, MyDalitzEEs const& dileptons, MyMCTracks const& tracks, aod::EMMCParticles const& mcparticles, aod::EMMCEvents const&)
  {
    THashList* list_ev_before = static_cast<THashList*>(fMainList->FindObject("Event")->FindObject(event_types[0].data()));
    THashList* list_ev_after = static_cast<THashList*>(fMainList->FindObject("Event")->FindObject(event_types[1].data()));
    THashList* list_dalitzee = static_cast<THashList*>(fMainList->FindObject("DalitzEE"));
    THashList* list_track = static_cast<THashList*>(fMainList->FindObject("Track"));
    double values[4] = {0, 0, 0, 0};
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
            if (IsPhysicalPrimary(mcmother.emmcevent(), mcmother, mcparticles)) {
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
              values[3] = uls_pair.phiv();
              reinterpret_cast<THnSparseF*>(list_dalitzee_cut->FindObject("hs_dilepton_uls_same"))->Fill(values);

              if (mcmother.pdgCode() == 111) {
                reinterpret_cast<TH2F*>(list_dalitzee_cut->FindObject("hMvsPhiV_Pi0"))->Fill(uls_pair.phiv(), uls_pair.mass());
                reinterpret_cast<TH2F*>(list_dalitzee_cut->FindObject("hMvsOPA_Pi0"))->Fill(uls_pair.opangle(), uls_pair.mass());
              } else if (mcmother.pdgCode() == 221) {
                reinterpret_cast<TH2F*>(list_dalitzee_cut->FindObject("hMvsPhiV_Eta"))->Fill(uls_pair.phiv(), uls_pair.mass());
                reinterpret_cast<TH2F*>(list_dalitzee_cut->FindObject("hMvsOPA_Eta"))->Fill(uls_pair.opangle(), uls_pair.mass());
              }

              nuls++;
              for (auto& track : {pos, ele}) {
                if (std::find(used_trackIds.begin(), used_trackIds.end(), track.globalIndex()) == used_trackIds.end()) {
                  o2::aod::pwgem::photon::histogram::FillHistClass<EMHistType::kTrack>(list_track_cut, "", track);
                  used_trackIds.emplace_back(track.globalIndex());
                  auto mctrack = track.template emmcparticle_as<aod::EMMCParticles>();
                  reinterpret_cast<TH2F*>(fMainList->FindObject("Track")->FindObject(cut.GetName())->FindObject("hPtGen_DeltaPtOverPtGen"))->Fill(mctrack.pt(), (track.pt() - mctrack.pt()) / mctrack.pt());
                  reinterpret_cast<TH2F*>(fMainList->FindObject("Track")->FindObject(cut.GetName())->FindObject("hPtGen_DeltaEta"))->Fill(mctrack.pt(), track.eta() - mctrack.eta());
                  reinterpret_cast<TH2F*>(fMainList->FindObject("Track")->FindObject(cut.GetName())->FindObject("hPtGen_DeltaPhi"))->Fill(mctrack.pt(), track.phi() - mctrack.phi());
                }
              }
            } // end of LF
          } else if (photonid > 0) {
            auto mcphoton = mcparticles.iteratorAt(photonid);
            if (IsPhysicalPrimary(mcphoton.emmcevent(), mcphoton, mcparticles) && IsEleFromPC(elemc, mcparticles) && IsEleFromPC(posmc, mcparticles)) {
              reinterpret_cast<TH2F*>(list_dalitzee_cut->FindObject("hMvsPhiV_Photon"))->Fill(uls_pair.phiv(), uls_pair.mass());
              reinterpret_cast<TH2F*>(list_dalitzee_cut->FindObject("hMvsOPA_Photon"))->Fill(uls_pair.opangle(), uls_pair.mass());
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
  PresliceUnsorted<aod::EMMCParticles> perMcCollision = aod::emmcparticle::emmceventId;
  void processGen(MyCollisions const& collisions, aod::EMMCEvents const&, aod::EMMCParticles const& mcparticles)
  {
    // loop over mc stack and fill histograms for pure MC truth signals
    // all MC tracks which belong to the MC event corresponding to the current reconstructed event
    for (auto& collision : grouped_collisions) {
      float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }
      auto mccollision = collision.emmcevent();
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

      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }

      reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hCollisionCounter"))->Fill(4.0);
      reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hZvtx_after"))->Fill(mccollision.posZ());
      auto posTracks_per_coll = posTracks->sliceByCachedUnsorted(o2::aod::emmcparticle::emmceventId, mccollision.globalIndex(), cache);
      auto negTracks_per_coll = negTracks->sliceByCachedUnsorted(o2::aod::emmcparticle::emmceventId, mccollision.globalIndex(), cache);
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
        if (IsPhysicalPrimary(mcmother.emmcevent(), mcmother, mcparticles)) {
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
