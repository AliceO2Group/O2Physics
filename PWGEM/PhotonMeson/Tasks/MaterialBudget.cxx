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
// This code loops over v0 photons for studying material budget.
//    Please write to: daiki.sekihata@cern.ch

#include <cstring>
#include <iterator>

#include "TString.h"
#include "Math/Vector4D.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/RecoDecay.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/PairUtilities.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/Core/PairCut.h"
#include "PWGEM/PhotonMeson/Core/CutsLibrary.h"
#include "PWGEM/PhotonMeson/Core/HistogramsLibrary.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::photonpair;

using MyCollisions = soa::Join<aod::EMReducedEvents, aod::EMReducedEventsMult, aod::EMReducedEventsCent, aod::EMReducedEventsNgPCM>;
using MyCollision = MyCollisions::iterator;

using MyV0Photons = soa::Join<aod::V0PhotonsKF, aod::V0Recalculation, aod::V0KFEMReducedEventIds>;
using MyV0Photon = MyV0Photons::iterator;

struct MaterialBudget {

  Configurable<float> CentMin{"CentMin", -1, "min. centrality"};
  Configurable<float> CentMax{"CentMax", 999, "max. centrality"};
  Configurable<std::string> CentEstimator{"CentEstimator", "FT0M", "centrality estimator"};

  Configurable<std::string> fConfigTagCuts{"cfgTagCuts", "qc", "Comma separated list of V0 photon cuts for tag"};
  Configurable<std::string> fConfigProbeCuts{"cfgProbeCuts", "qc,wwire_ib", "Comma separated list of V0 photon cuts for probe"};
  Configurable<std::string> fConfigPairCuts{"cfgPairCuts", "nocut", "Comma separated list of pair cuts"};
  Configurable<bool> fDoMixing{"DoMixing", false, "do event mixing"};

  OutputObj<THashList> fOutputEvent{"Event"};
  OutputObj<THashList> fOutputV0{"V0"};
  OutputObj<THashList> fOutputPair{"Pair"}; // 2-photon pair
  THashList* fMainList = new THashList();

  std::vector<V0PhotonCut> fTagCuts;
  std::vector<V0PhotonCut> fProbeCuts;
  std::vector<PairCut> fPairCuts;

  std::vector<std::string> fPairNames;
  void init(InitContext& context)
  {
    if (context.mOptions.get<bool>("processMB")) {
      fPairNames.push_back("PCMPCM");
    }

    DefineTagCuts();
    DefineProbeCuts();
    DefinePairCuts();
    addhistograms();

    fOutputEvent.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("Event")));
    fOutputV0.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("V0")));
    fOutputPair.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("Pair")));
  }

  template <typename TCuts1, typename TCuts2, typename TCuts3>
  void add_pair_histograms(THashList* list_pair, const std::string pairname, TCuts1 const& tagcuts, TCuts2 const& probecuts, TCuts3 const& cuts3)
  {
    for (auto& cut1 : tagcuts) {
      for (auto& cut2 : probecuts) {
        std::string cutname1 = cut1.GetName();
        std::string cutname2 = cut2.GetName();

        THashList* list_pair_subsys = reinterpret_cast<THashList*>(list_pair->FindObject(pairname.data()));
        std::string photon_cut_name = cutname1 + "_" + cutname2;
        o2::aod::emphotonhistograms::AddHistClass(list_pair_subsys, photon_cut_name.data());
        THashList* list_pair_subsys_photoncut = reinterpret_cast<THashList*>(list_pair_subsys->FindObject(photon_cut_name.data()));

        for (auto& cut3 : cuts3) {
          std::string pair_cut_name = cut3.GetName();
          o2::aod::emphotonhistograms::AddHistClass(list_pair_subsys_photoncut, pair_cut_name.data());
          THashList* list_pair_subsys_paircut = reinterpret_cast<THashList*>(list_pair_subsys_photoncut->FindObject(pair_cut_name.data()));
          o2::aod::emphotonhistograms::DefineHistograms(list_pair_subsys_paircut, "material_budget_study", "Pair");
        } // end of cut3 loop pair cut
      }   // end of cut2 loop
    }     // end of cut1 loop
  }

  static constexpr std::string_view pairnames[9] = {"PCMPCM", "PHOSPHOS", "EMCEMC", "PCMPHOS", "PCMEMC", "PCMDalitzEE", "PCMDalitzMuMu", "PHOSEMC", "DalitzEEDalitzEE"};
  void addhistograms()
  {
    fMainList->SetOwner(true);
    fMainList->SetName("fMainList");

    // create sub lists first.
    o2::aod::emphotonhistograms::AddHistClass(fMainList, "Event");
    THashList* list_ev = reinterpret_cast<THashList*>(fMainList->FindObject("Event"));

    o2::aod::emphotonhistograms::AddHistClass(fMainList, "Pair");
    THashList* list_pair = reinterpret_cast<THashList*>(fMainList->FindObject("Pair"));

    o2::aod::emphotonhistograms::AddHistClass(fMainList, "V0");
    THashList* list_v0 = reinterpret_cast<THashList*>(fMainList->FindObject("V0"));

    // for V0s
    for (const auto& cut : fProbeCuts) {
      const char* cutname = cut.GetName();
      THashList* list_v0_cut = o2::aod::emphotonhistograms::AddHistClass(list_v0, cutname);
      o2::aod::emphotonhistograms::DefineHistograms(list_v0_cut, "material_budget_study", "V0");
    }

    for (auto& pairname : fPairNames) {
      LOGF(info, "Enabled pairs = %s", pairname.data());

      o2::aod::emphotonhistograms::AddHistClass(list_ev, pairname.data());
      THashList* list_ev_pair = reinterpret_cast<THashList*>(list_ev->FindObject(pairname.data()));
      o2::aod::emphotonhistograms::DefineHistograms(list_ev_pair, "Event");

      o2::aod::emphotonhistograms::AddHistClass(list_pair, pairname.data());

      if (pairname == "PCMPCM") {
        add_pair_histograms(list_pair, pairname, fTagCuts, fProbeCuts, fPairCuts);
      }

    } // end of pair name loop
  }

  void DefineTagCuts()
  {
    TString cutNamesStr = fConfigTagCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        const char* cutname = objArray->At(icut)->GetName();
        LOGF(info, "add cut : %s", cutname);
        fTagCuts.push_back(*pcmcuts::GetCut(cutname));
      }
    }
    LOGF(info, "Number of Tag cuts = %d", fTagCuts.size());
  }

  void DefineProbeCuts()
  {
    TString cutNamesStr = fConfigProbeCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        const char* cutname = objArray->At(icut)->GetName();
        LOGF(info, "add cut : %s", cutname);
        fProbeCuts.push_back(*pcmcuts::GetCut(cutname));
      }
    }
    LOGF(info, "Number of Probe cuts = %d", fProbeCuts.size());
  }

  void DefinePairCuts()
  {
    TString cutNamesStr = fConfigPairCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        const char* cutname = objArray->At(icut)->GetName();
        LOGF(info, "add cut : %s", cutname);
        fPairCuts.push_back(*paircuts::GetCut(cutname));
      }
    }
    LOGF(info, "Number of Pair cuts = %d", fPairCuts.size());
  }

  Preslice<MyV0Photons> perCollision_pcm = aod::v0photonkf::emreducedeventId;

  template <PairType pairtype, typename TG1, typename TG2, typename TCut1, typename TCut2>
  bool IsSelectedPair(TG1 const& g1, TG2 const& g2, TCut1 const& cut1, TCut2 const& cut2)
  {
    return o2::aod::photonpair::IsSelectedPair<aod::V0Legs, aod::V0Legs>(g1, g2, cut1, cut2);
  }

  template <PairType pairtype, typename TEvents, typename TPhotons, typename TPreslice, typename TCuts, typename TLegs>
  void fillsinglephoton(TEvents const& collisions, TPhotons const& photons, TPreslice const& perCollision, TCuts const& cuts, TLegs const& legs)
  {
    THashList* list_ev_pair = static_cast<THashList*>(fMainList->FindObject("Event")->FindObject(pairnames[pairtype].data()));
    THashList* list_v0 = static_cast<THashList*>(fMainList->FindObject("V0"));
    double value[4] = {0.f};
    for (auto& collision : collisions) {
      reinterpret_cast<TH1F*>(fMainList->FindObject("Event")->FindObject(pairnames[pairtype].data())->FindObject("hZvtx_before"))->Fill(collision.posZ());
      reinterpret_cast<TH1F*>(fMainList->FindObject("Event")->FindObject(pairnames[pairtype].data())->FindObject("hCollisionCounter"))->Fill(1.0); // all
      if (!collision.sel8()) {
        continue;
      }
      reinterpret_cast<TH1F*>(fMainList->FindObject("Event")->FindObject(pairnames[pairtype].data())->FindObject("hCollisionCounter"))->Fill(2.0); // FT0VX i.e. FT0and
      if (collision.numContrib() < 0.5) {
        continue;
      }
      reinterpret_cast<TH1F*>(fMainList->FindObject("Event")->FindObject(pairnames[pairtype].data())->FindObject("hCollisionCounter"))->Fill(3.0); // Ncontrib > 0
      if (abs(collision.posZ()) > 10.0) {
        continue;
      }
      reinterpret_cast<TH1F*>(fMainList->FindObject("Event")->FindObject(pairnames[pairtype].data())->FindObject("hZvtx_after"))->Fill(collision.posZ());
      reinterpret_cast<TH1F*>(fMainList->FindObject("Event")->FindObject(pairnames[pairtype].data())->FindObject("hCollisionCounter"))->Fill(4.0); // |Zvtx| < 10 cm
      o2::aod::emphotonhistograms::FillHistClass<EMHistType::kEvent>(list_ev_pair, "", collision);

      auto photons_coll = photons.sliceBy(perCollision, collision.globalIndex());
      for (auto& cut : cuts) {
        for (auto& photon : photons_coll) {
          if (!cut.template IsSelected<TLegs>(photon)) {
            continue;
          }

          float phi_cp = atan2(photon.vy(), photon.vx());
          float eta_cp = std::atanh(photon.vz() / sqrt(pow(photon.vx(), 2) + pow(photon.vy(), 2) + pow(photon.vz(), 2)));
          value[0] = photon.pt();
          value[1] = photon.v0radius();
          value[2] = phi_cp > 0 ? phi_cp : phi_cp + TMath::TwoPi();
          value[3] = eta_cp;
          reinterpret_cast<THnSparseF*>(list_v0->FindObject(cut.GetName())->FindObject("hs_conv_point"))->Fill(value);

        } // end of photon loop
      }   // end of cut loop

    } // end of collision loop
  }

  template <PairType pairtype, typename TEvents, typename TPhotons1, typename TPhotons2, typename TPreslice1, typename TPreslice2, typename TCuts1, typename TCuts2, typename TPairCuts, typename TLegs>
  void SameEventPairing(TEvents const& collisions, TPhotons1 const& photons1, TPhotons2 const& photons2, TPreslice1 const& perCollision1, TPreslice2 const& perCollision2, TCuts1 const& tagcuts, TCuts2 const& probecuts, TPairCuts const& paircuts, TLegs const& legs)
  {
    THashList* list_pair_ss = static_cast<THashList*>(fMainList->FindObject("Pair")->FindObject(pairnames[pairtype].data()));

    for (auto& collision : collisions) {
      auto photons1_coll = photons1.sliceBy(perCollision1, collision.globalIndex());
      auto photons2_coll = photons2.sliceBy(perCollision2, collision.globalIndex());

      double value[6] = {0.f};
      float phi_cp2 = 0.f, eta_cp2 = 0.f;
      for (auto& tagcut : tagcuts) {
        for (auto& probecut : probecuts) {
          for (auto& g1 : photons1_coll) {
            for (auto& g2 : photons2_coll) {
              if (g1.globalIndex() == g2.globalIndex()) {
                continue;
              }
              if (!IsSelectedPair<pairtype>(g1, g2, tagcut, probecut)) {
                continue;
              }
              for (auto& paircut : paircuts) {
                if (!paircut.IsSelected(g1, g2)) {
                  continue;
                }

                ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.); // tag
                ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.); // probe
                ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
                phi_cp2 = atan2(g2.vy(), g2.vx());
                eta_cp2 = std::atanh(g2.vz() / sqrt(pow(g2.vx(), 2) + pow(g2.vy(), 2) + pow(g2.vz(), 2)));
                value[0] = v12.M();
                value[1] = g1.pt();
                value[2] = g2.pt();
                value[3] = g2.v0radius();
                value[4] = phi_cp2 > 0.f ? phi_cp2 : phi_cp2 + TMath::TwoPi();
                value[5] = eta_cp2;
                reinterpret_cast<THnSparseF*>(list_pair_ss->FindObject(Form("%s_%s", tagcut.GetName(), probecut.GetName()))->FindObject(paircut.GetName())->FindObject("hs_conv_point_same"))->Fill(value);
              } // end of pair cut loop
            }   // end of g2 loop
          }     // end of g1 loop
        }       // end of probecut loop
      }         // end of tagcut loop
    }           // end of collision loop
  }

  Configurable<int> ndepth{"ndepth", 10, "depth for event mixing"};
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis ConfMultBins{"ConfMultBins", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 200.0f, 1e+10f}, "Mixing bins - multiplicity"};
  using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultNTracksPV>;
  BinningType colBinning{{ConfVtxBins, ConfMultBins}, true};

  template <PairType pairtype, typename TEvents, typename TPhotons1, typename TPhotons2, typename TPreslice1, typename TPreslice2, typename TCuts1, typename TCuts2, typename TPairCuts, typename TLegs>
  void MixedEventPairing(TEvents const& collisions, TPhotons1 const& photons1, TPhotons2 const& photons2, TPreslice1 const& perCollision1, TPreslice2 const& perCollision2, TCuts1 const& tagcuts, TCuts2 const& probecuts, TPairCuts const& paircuts, TLegs const& legs)
  {
    THashList* list_pair_ss = static_cast<THashList*>(fMainList->FindObject("Pair")->FindObject(pairnames[pairtype].data()));
    // LOGF(info, "Number of collisions after filtering: %d", collisions.size());
    for (auto& [collision1, collision2] : soa::selfCombinations(colBinning, ndepth, -1, collisions, collisions)) { // internally, CombinationsStrictlyUpperIndexPolicy(collisions, collisions) is called.

      // LOGF(info, "Mixed event globalIndex: (%d, %d) , ngpcm: (%d, %d), ngphos: (%d, %d), ngemc: (%d, %d)",
      //     collision1.globalIndex(), collision2.globalIndex(), collision1.ngpcm(), collision2.ngpcm(), collision1.ngphos(), collision2.ngphos(), collision1.ngemc(), collision2.ngemc());

      auto photons_coll1 = photons1.sliceBy(perCollision1, collision1.globalIndex());
      auto photons_coll2 = photons2.sliceBy(perCollision2, collision2.globalIndex());
      // LOGF(info, "collision1: posZ = %f, numContrib = %d , sel8 = %d | collision2: posZ = %f, numContrib = %d , sel8 = %d",
      //     collision1.posZ(), collision1.numContrib(), collision1.sel8(), collision2.posZ(), collision2.numContrib(), collision2.sel8());

      double value[6] = {0.f};
      float phi_cp2 = 0.f, eta_cp2 = 0.f;
      for (auto& tagcut : tagcuts) {
        for (auto& probecut : probecuts) {
          for (auto& g1 : photons_coll1) {
            for (auto& g2 : photons_coll2) {
              if (!IsSelectedPair<pairtype>(g1, g2, tagcut, probecut)) {
                continue;
              }
              // LOGF(info, "Mixed event photon pair: (%d, %d) from events (%d, %d), photon event: (%d, %d)", g1.index(), g2.index(), collision1.index(), collision2.index(), g1.globalIndex(), g2.globalIndex());

              for (auto& paircut : paircuts) {
                if (!paircut.IsSelected(g1, g2)) {
                  continue;
                }

                ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.); // tag
                ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.); // probe
                ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
                phi_cp2 = atan2(g2.vy(), g2.vx());
                eta_cp2 = std::atanh(g2.vz() / sqrt(pow(g2.vx(), 2) + pow(g2.vy(), 2) + pow(g2.vz(), 2)));
                value[0] = v12.M();
                value[1] = g1.pt();
                value[2] = g2.pt();
                value[3] = g2.v0radius();
                value[4] = phi_cp2 > 0.f ? phi_cp2 : phi_cp2 + TMath::TwoPi();
                value[5] = eta_cp2;
                reinterpret_cast<THnSparseF*>(list_pair_ss->FindObject(Form("%s_%s", tagcut.GetName(), probecut.GetName()))->FindObject(paircut.GetName())->FindObject("hs_conv_point_mix"))->Fill(value);

              } // end of pair cut loop
            }   // end of g2 loop
          }     // end of g1 loop
        }       // end of probecut loop
      }         // end of tagcut loop
    }           // end of different collision combinations
  }

  Partition<MyCollisions> grouped_collisions = CentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < CentMax; // this goes to same event.
  Filter collisionFilter_common = nabs(o2::aod::collision::posZ) < 10.f && o2::aod::collision::numContrib > (uint16_t)0 && o2::aod::evsel::sel8 == true && CentMin < o2::aod::cent::centFT0M&& o2::aod::cent::centFT0M < CentMax;
  Filter collisionFilter_subsys = o2::aod::emreducedevent::ngpcm >= 2;
  using MyFilteredCollisions = soa::Filtered<MyCollisions>; // this goes to mixed event.

  void processMB(MyCollisions const& collisions, MyFilteredCollisions const& filtered_collisions, MyV0Photons const& v0photons, aod::V0Legs const& legs)
  {
    fillsinglephoton<PairType::kPCMPCM>(grouped_collisions, v0photons, perCollision_pcm, fProbeCuts, legs);
    SameEventPairing<PairType::kPCMPCM>(filtered_collisions, v0photons, v0photons, perCollision_pcm, perCollision_pcm, fTagCuts, fProbeCuts, fPairCuts, legs);
    if (fDoMixing) {
      MixedEventPairing<PairType::kPCMPCM>(filtered_collisions, v0photons, v0photons, perCollision_pcm, perCollision_pcm, fTagCuts, fProbeCuts, fPairCuts, legs);
    }
  }

  void processDummy(MyCollisions::iterator const& collision) {}

  PROCESS_SWITCH(MaterialBudget, processMB, "process material budget", false);
  PROCESS_SWITCH(MaterialBudget, processDummy, "Dummy function", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MaterialBudget>(cfgc, TaskName{"material-budget"})};
}
