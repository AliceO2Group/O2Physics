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
#include "PWGEM/PhotonMeson/Utils/MCUtilities.h"
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

using MyMCV0Legs = soa::Join<aod::V0Legs, aod::V0LegMCLabels>;
using MyMCV0Leg = MyMCV0Legs::iterator;

struct MaterialBudgetMC {

  Configurable<float> CentMin{"CentMin", -1, "min. centrality"};
  Configurable<float> CentMax{"CentMax", 999, "max. centrality"};
  Configurable<std::string> CentEstimator{"CentEstimator", "FT0M", "centrality estimator"};

  Configurable<float> maxY{"maxY", 0.9, "maximum rapidity for generated particles"};
  Configurable<std::string> fConfigTagCuts{"cfgTagCuts", "qc", "Comma separated list of V0 photon cuts for tag"};
  Configurable<std::string> fConfigProbeCuts{"cfgProbeCuts", "qc,wwire_ib", "Comma separated list of V0 photon cuts for probe"};
  Configurable<std::string> fConfigPairCuts{"cfgPairCuts", "nocut", "Comma separated list of pair cuts"};

  OutputObj<THashList> fOutputEvent{"Event"};
  OutputObj<THashList> fOutputV0{"V0"};
  OutputObj<THashList> fOutputPair{"Pair"}; // 2-photon pair
  OutputObj<THashList> fOutputGen{"Generated"};
  THashList* fMainList = new THashList();

  std::vector<V0PhotonCut> fTagCuts;
  std::vector<V0PhotonCut> fProbeCuts;
  std::vector<PairCut> fPairCuts;

  std::vector<std::string> fPairNames;
  void init(InitContext& context)
  {
    if (context.mOptions.get<bool>("processMBMC")) {
      fPairNames.push_back("PCMPCM");
    }

    DefineTagCuts();
    DefineProbeCuts();
    DefinePairCuts();
    addhistograms();

    fOutputEvent.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("Event")));
    fOutputV0.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("V0")));
    fOutputPair.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("Pair")));
    fOutputGen.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("Generated")));
  }

  template <typename TCuts1, typename TCuts2, typename TCuts3>
  void add_pair_histograms(THashList* list_pair, const std::string pairname, TCuts1 const& tagcuts, TCuts2 const& probecuts, TCuts3 const& cuts3)
  {
    for (auto& tagcut : tagcuts) {
      for (auto& probecut : probecuts) {
        std::string cutname1 = tagcut.GetName();
        std::string cutname2 = probecut.GetName();

        // if (cutname1 == cutname2) {
        //   continue;
        // }

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
      }   // end of probecut loop
    }     // end of tagcut loop
  }

  static constexpr std::string_view pairnames[8] = {"PCMPCM", "PHOSPHOS", "EMCEMC", "PCMPHOS", "PCMEMC", "PCMDalitzEE", "PCMDalitzMuMu", "PHOSEMC"};
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

    o2::aod::emphotonhistograms::AddHistClass(fMainList, "Generated");
    THashList* list_gen = reinterpret_cast<THashList*>(fMainList->FindObject("Generated"));
    o2::aod::emphotonhistograms::DefineHistograms(list_gen, "Generated", "Photon");
    o2::aod::emphotonhistograms::DefineHistograms(list_gen, "Generated", "ConversionStudy");
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
    LOGF(info, "Number of Tag PCM cuts = %d", fTagCuts.size());
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
    LOGF(info, "Number of Probe PCM cuts = %d", fProbeCuts.size());
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
  bool IsSelectedPair(TG1 const& g1, TG2 const& g2, TCut1 const& tagcut, TCut2 const& probecut)
  {
    return o2::aod::photonpair::IsSelectedPair<MyMCV0Legs, MyMCV0Legs>(g1, g2, tagcut, probecut);
  }

  template <PairType pairtype, typename TEvents, typename TPhotons, typename TPreslice, typename TCuts, typename TLegs, typename TMCParticles, typename TMCEvents>
  void fillsinglephoton(TEvents const& collisions, TPhotons const& photons, TPreslice const& perCollision, TCuts const& cuts, TLegs const& legs, TMCParticles const& mcparticles, TMCEvents const&)
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

          auto pos = photon.template posTrack_as<MyMCV0Legs>();
          auto ele = photon.template negTrack_as<MyMCV0Legs>();
          auto posmc = pos.template emmcparticle_as<aod::EMMCParticles>();
          auto elemc = ele.template emmcparticle_as<aod::EMMCParticles>();
          int photonid = FindCommonMotherFrom2Prongs(posmc, elemc, -11, 11, 22, mcparticles);

          bool is_photon_physical_primary = false;
          if (photonid > 0) {
            auto mcphoton = mcparticles.iteratorAt(photonid);
            is_photon_physical_primary = IsPhysicalPrimary(mcphoton.emreducedmcevent(), mcphoton, mcparticles);
          }
          if (!is_photon_physical_primary) {
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

  template <PairType pairtype, typename TEvents, typename TPhotons1, typename TPhotons2, typename TPreslice1, typename TPreslice2, typename TCuts1, typename TCuts2, typename TPairCuts, typename TLegs, typename TMCParticles, typename TMCEvents>
  void TruePairing(TEvents const& collisions, TPhotons1 const& photons1, TPhotons2 const& photons2, TPreslice1 const& perCollision1, TPreslice2 const& perCollision2, TCuts1 const& tagcuts, TCuts2 const& probecuts, TPairCuts const& paircuts, TLegs const& legs, TMCParticles const& mcparticles, TMCEvents const&)
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

                auto pos1 = g1.template posTrack_as<MyMCV0Legs>();
                auto ele1 = g1.template negTrack_as<MyMCV0Legs>();
                auto pos2 = g2.template posTrack_as<MyMCV0Legs>();
                auto ele2 = g2.template negTrack_as<MyMCV0Legs>();

                auto pos1mc = pos1.template emmcparticle_as<aod::EMMCParticles>();
                auto ele1mc = ele1.template emmcparticle_as<aod::EMMCParticles>();
                auto pos2mc = pos2.template emmcparticle_as<aod::EMMCParticles>();
                auto ele2mc = ele2.template emmcparticle_as<aod::EMMCParticles>();
                // LOGF(info,"pos1mc.globalIndex() = %d , ele1mc.globalIndex() = %d , pos2mc.globalIndex() = %d , ele2mc.globalIndex() = %d", pos1mc.globalIndex(), ele1mc.globalIndex(), pos2mc.globalIndex(), ele2mc.globalIndex());

                int photonid1 = FindCommonMotherFrom2Prongs(pos1mc, ele1mc, -11, 11, 22, mcparticles);
                if (photonid1 < 0) {
                  continue;
                }
                auto g1mc = mcparticles.iteratorAt(photonid1);

                int photonid2 = FindCommonMotherFrom2Prongs(pos2mc, ele2mc, -11, 11, 22, mcparticles);
                if (photonid2 < 0) {
                  continue;
                }
                auto g2mc = mcparticles.iteratorAt(photonid2);

                int pi0id = FindCommonMotherFrom2Prongs(g1mc, g2mc, 22, 22, 111, mcparticles);
                if (pi0id < 0) {
                  continue;
                }
                auto pi0mc = mcparticles.iteratorAt(pi0id);
                if (!IsPhysicalPrimary(pi0mc.emreducedmcevent(), pi0mc, mcparticles)) {
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

  Partition<MyCollisions> grouped_collisions = CentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < CentMax; // this goes to same event.
  Filter collisionFilter_common = nabs(o2::aod::collision::posZ) < 10.f && o2::aod::collision::numContrib > (uint16_t)0 && o2::aod::evsel::sel8 == true && CentMin < o2::aod::cent::centFT0M&& o2::aod::cent::centFT0M < CentMax;
  Filter collisionFilter_subsys = o2::aod::emreducedevent::ngpcm >= 2;
  using MyFilteredCollisions = soa::Filtered<MyCollisions>; // this goes to mixed event.

  void processMBMC(MyCollisions const& collisions, MyFilteredCollisions const& filtered_collisions, MyV0Photons const& v0photons, MyMCV0Legs const& legs, aod::EMMCParticles const& mcparticles, aod::EMReducedMCEvents const& mccollisions)
  {
    fillsinglephoton<PairType::kPCMPCM>(grouped_collisions, v0photons, perCollision_pcm, fProbeCuts, legs, mcparticles, mccollisions);
    TruePairing<PairType::kPCMPCM>(filtered_collisions, v0photons, v0photons, perCollision_pcm, perCollision_pcm, fTagCuts, fProbeCuts, fPairCuts, legs, mcparticles, mccollisions);
  }

  PresliceUnsorted<aod::EMMCParticles> perMcCollision = aod::emmcparticle::emreducedmceventId;
  void processGen(soa::Join<MyCollisions, aod::EMReducedMCEventLabels> const& collisions, aod::EMReducedMCEvents const&, aod::EMMCParticles const& mcparticles)
  {
    // loop over mc stack and fill histograms for pure MC truth signals
    // all MC tracks which belong to the MC event corresponding to the current reconstructed event

    for (auto& collision : collisions) {
      auto mccollision = collision.emreducedmcevent();
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
      reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hCollisionCounter"))->Fill(4.0);
      reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hZvtx_after"))->Fill(mccollision.posZ());

      auto mctracks_coll = mcparticles.sliceBy(perMcCollision, mccollision.globalIndex());
      for (auto& mctrack : mctracks_coll) {
        if (abs(mctrack.y()) > maxY) {
          continue;
        }
        if (abs(mctrack.pdgCode()) == 22 && IsPhysicalPrimary(mctrack.emreducedmcevent(), mctrack, mcparticles)) {
          reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hPt_Photon"))->Fill(mctrack.pt());
          reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hY_Photon"))->Fill(mctrack.y());
          reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hPhi_Photon"))->Fill(mctrack.phi());
        }

        int photonid = IsEleFromPC(mctrack, mcparticles);
        if (photonid > 0) {
          auto mcphoton = mcparticles.iteratorAt(photonid);
          if (!IsPhysicalPrimary(mcphoton.emreducedmcevent(), mcphoton, mcparticles)) {
            continue;
          }
          float rxy = sqrt(pow(mctrack.vx(), 2) + pow(mctrack.vy(), 2));
          reinterpret_cast<TH2F*>(fMainList->FindObject("Generated")->FindObject("hPhotonRZ"))->Fill(mctrack.vz(), rxy);
          reinterpret_cast<TH2F*>(fMainList->FindObject("Generated")->FindObject("hPhotonRxy"))->Fill(mctrack.vx(), mctrack.vy());
          reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hPhotonPhivsRxy"))->Fill(mctrack.phi(), rxy);
        }
      }
    }
  }

  void processDummy(MyCollisions::iterator const& collision) {}

  PROCESS_SWITCH(MaterialBudgetMC, processMBMC, "process material budget", false);
  PROCESS_SWITCH(MaterialBudgetMC, processGen, "process generated information", false);
  PROCESS_SWITCH(MaterialBudgetMC, processDummy, "Dummy function", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MaterialBudgetMC>(cfgc, TaskName{"material-budget-mc"})};
}
