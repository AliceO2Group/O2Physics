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

using MyV0Photons = soa::Join<aod::V0PhotonsKF, aod::V0Recalculation, aod::V0KFEMReducedEventIds>;
using MyV0Photon = MyV0Photons::iterator;

struct MaterialBudgetMC {
  using MyMCV0Legs = soa::Join<aod::V0Legs, aod::EMMCParticleLabels>;

  Configurable<float> CentMin{"CentMin", -1, "min. centrality"};
  Configurable<float> CentMax{"CentMax", 999, "max. centrality"};
  Configurable<std::string> CentEstimator{"CentEstimator", "FT0M", "centrality estimator"};

  Configurable<float> maxY{"maxY", 999.f, "maximum rapidity for generated particles"};
  Configurable<std::string> fConfigTagPCMCuts{"cfgTagPCMCuts", "qc_ITSTPC", "Comma separated list of V0 photon cuts for tag"};
  Configurable<std::string> fConfigProbePCMCuts{"cfgProbePCMCuts", "qc_ITSTPC,qc_TPConly,qc_ITSonly", "Comma separated list of V0 photon cuts for probe"};
  Configurable<std::string> fConfigPairCuts{"cfgPairCuts", "nocut", "Comma separated list of pair cuts"};

  OutputObj<THashList> fOutputEvent{"Event"};
  OutputObj<THashList> fOutputV0{"V0"};
  OutputObj<THashList> fOutputPair{"Pair"}; // 2-photon pair
  OutputObj<THashList> fOutputGen{"Generated"};
  THashList* fMainList = new THashList();

  std::vector<V0PhotonCut> fTagPCMCuts;
  std::vector<V0PhotonCut> fProbePCMCuts;
  std::vector<PairCut> fPairCuts;

  std::vector<std::string> fPairNames;
  void init(InitContext& context)
  {
    if (context.mOptions.get<bool>("processMBMC")) {
      fPairNames.push_back("PCMPCM");
    }

    DefineTagPCMCuts();
    DefineProbePCMCuts();
    DefinePairCuts();
    addhistograms();

    fOutputEvent.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("Event")));
    fOutputV0.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("V0")));
    fOutputPair.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("Pair")));
    fOutputGen.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("Generated")));
  }

  template <typename TCuts1, typename TCuts2, typename TCuts3>
  void add_pair_histograms(THashList* list_pair, const std::string pairname, TCuts1 const& cuts1, TCuts2 const& cuts2, TCuts3 const& cuts3)
  {
    for (auto& cut1 : cuts1) {
      for (auto& cut2 : cuts2) {
        std::string cutname1 = cut1.GetName();
        std::string cutname2 = cut2.GetName();

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
      }   // end of cut2 loop
    }     // end of cut1 loop
  }

  static constexpr std::string_view pairnames[6] = {"PCMPCM", "PHOSPHOS", "EMCEMC", "PCMPHOS", "PCMEMC", "PHOSEMC"};
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
    for (const auto& cut : fProbePCMCuts) {
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
        add_pair_histograms(list_pair, pairname, fTagPCMCuts, fProbePCMCuts, fPairCuts);
      }

    } // end of pair name loop

    o2::aod::emphotonhistograms::AddHistClass(fMainList, "Generated");
    THashList* list_gen = reinterpret_cast<THashList*>(fMainList->FindObject("Generated"));
    o2::aod::emphotonhistograms::DefineHistograms(list_gen, "Generated", "Photon");
    o2::aod::emphotonhistograms::DefineHistograms(list_gen, "Generated", "ConversionStudy");
  }

  void DefineTagPCMCuts()
  {
    TString cutNamesStr = fConfigTagPCMCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        const char* cutname = objArray->At(icut)->GetName();
        LOGF(info, "add cut : %s", cutname);
        fTagPCMCuts.push_back(*pcmcuts::GetCut(cutname));
      }
    }
    LOGF(info, "Number of Tag PCM cuts = %d", fTagPCMCuts.size());
  }

  void DefineProbePCMCuts()
  {
    TString cutNamesStr = fConfigProbePCMCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        const char* cutname = objArray->At(icut)->GetName();
        LOGF(info, "add cut : %s", cutname);
        fProbePCMCuts.push_back(*pcmcuts::GetCut(cutname));
      }
    }
    LOGF(info, "Number of Probe PCM cuts = %d", fProbePCMCuts.size());
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
    return o2::aod::photonpair::IsSelectedPair<MyMCV0Legs, MyMCV0Legs>(g1, g2, cut1, cut2);
  }

  template <typename TEvents, typename TPhotons, typename TPreslice, typename TCuts, typename TLegs, typename TMCParticles, typename TMCEvents>
  void fillsinglephoton(TEvents const& collisions, TPhotons const& photons, TPreslice const& perCollision, TCuts const& cuts, TLegs const& legs, TMCParticles const& mcparticles, TMCEvents const&)
  {
    THashList* list_v0 = static_cast<THashList*>(fMainList->FindObject("V0"));
    double value[4] = {0.f};
    for (auto& collision : collisions) {
      if (!collision.sel8()) {
        continue;
      }
      if (collision.numContrib() < 0.5) {
        continue;
      }
      if (abs(collision.posZ()) > 10.0) {
        continue;
      }

      auto photons_coll = photons.sliceBy(perCollision, collision.collisionId());
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
          if (photonid < 0) { // check swap, true electron is reconstructed as positron and vice versa.
            photonid = FindCommonMotherFrom2Prongs(posmc, elemc, 11, -11, 22, mcparticles);
          }

          bool is_photon_physical_primary = false;
          if (photonid > 0) {
            auto mcphoton = mcparticles.iteratorAt(photonid);
            is_photon_physical_primary = IsPhysicalPrimary(mcphoton.emreducedmcevent(), mcphoton, mcparticles);
          }
          if (!is_photon_physical_primary) {
            continue;
          }

          float phi_cp = atan2(photon.recalculatedVtxY(), photon.recalculatedVtxX());
          float eta_cp = std::atanh(photon.recalculatedVtxZ() / sqrt(pow(photon.recalculatedVtxX(), 2) + pow(photon.recalculatedVtxY(), 2) + pow(photon.recalculatedVtxZ(), 2)));
          value[0] = photon.pt();
          value[1] = photon.recalculatedVtxR();
          value[2] = phi_cp > 0 ? phi_cp : phi_cp + TMath::TwoPi();
          value[3] = eta_cp;
          reinterpret_cast<THnSparseF*>(list_v0->FindObject(cut.GetName())->FindObject("hs_conv_point"))->Fill(value);

        } // end of photon loop
      }   // end of cut loop

    } // end of collision loop
  }

  template <PairType pairtype, typename TEvents, typename TPhotons1, typename TPhotons2, typename TPreslice1, typename TPreslice2, typename TCuts1, typename TCuts2, typename TPairCuts, typename TLegs, typename TMCParticles, typename TMCEvents>
  void SameEventPairing(TEvents const& collisions, TPhotons1 const& photons1, TPhotons2 const& photons2, TPreslice1 const& perCollision1, TPreslice2 const& perCollision2, TCuts1 const& cuts1, TCuts2 const& cuts2, TPairCuts const& paircuts, TLegs const& legs, TMCParticles const& mcparticles, TMCEvents const&)
  {
    THashList* list_ev_pair = static_cast<THashList*>(fMainList->FindObject("Event")->FindObject(pairnames[pairtype].data()));
    THashList* list_pair_ss = static_cast<THashList*>(fMainList->FindObject("Pair")->FindObject(pairnames[pairtype].data()));

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

      auto photons1_coll = photons1.sliceBy(perCollision1, collision.collisionId());
      auto photons2_coll = photons2.sliceBy(perCollision2, collision.collisionId());

      double value[9] = {0.f};
      float phi_cp1 = 0.f, eta_cp1 = 0.f;
      float phi_cp2 = 0.f, eta_cp2 = 0.f;
      for (auto& cut1 : cuts1) {
        for (auto& cut2 : cuts2) {
          // if (std::string(cut1.GetName()) == std::string(cut2.GetName())) {
          //   continue;
          // }
          for (auto& g1 : photons1_coll) {
            for (auto& g2 : photons2_coll) {
              if (g1.globalIndex() == g2.globalIndex()) {
                continue;
              }
              if (!IsSelectedPair<pairtype>(g1, g2, cut1, cut2)) {
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
                int photonid2 = FindCommonMotherFrom2Prongs(pos2mc, ele2mc, -11, 11, 22, mcparticles);

                if (photonid1 < 0) { // check swap, true electron is reconstructed as positron and vice versa.
                  photonid1 = FindCommonMotherFrom2Prongs(pos1mc, ele1mc, 11, -11, 22, mcparticles);
                }
                if (photonid2 < 0) { // check swap, true electron is reconstructed as positron and vice versa.
                  photonid2 = FindCommonMotherFrom2Prongs(pos2mc, ele2mc, 11, -11, 22, mcparticles);
                }

                // LOGF(info,"photonid1 = %d , photonid2 = %d", photonid1, photonid2);
                if (photonid1 < 0 || photonid2 < 0) {
                  continue;
                }

                auto g1mc = mcparticles.iteratorAt(photonid1);
                auto g2mc = mcparticles.iteratorAt(photonid2);
                int pi0id = FindCommonMotherFrom2Prongs(g1mc, g2mc, 22, 22, 111, mcparticles);

                bool is_pi0_physical_primary = false;
                if (pi0id > 0) {
                  auto pi0mc = mcparticles.iteratorAt(pi0id);
                  is_pi0_physical_primary = IsPhysicalPrimary(pi0mc.emreducedmcevent(), pi0mc, mcparticles);
                }
                if (!is_pi0_physical_primary) {
                  continue;
                }

                ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.); // tag
                ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.); // probe
                ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
                phi_cp1 = atan2(g1.recalculatedVtxY(), g1.recalculatedVtxX());
                eta_cp1 = std::atanh(g1.recalculatedVtxZ() / sqrt(pow(g1.recalculatedVtxX(), 2) + pow(g1.recalculatedVtxY(), 2) + pow(g1.recalculatedVtxZ(), 2)));
                phi_cp2 = atan2(g2.recalculatedVtxY(), g2.recalculatedVtxX());
                eta_cp2 = std::atanh(g2.recalculatedVtxZ() / sqrt(pow(g2.recalculatedVtxX(), 2) + pow(g2.recalculatedVtxY(), 2) + pow(g2.recalculatedVtxZ(), 2)));
                value[0] = v12.M();
                value[1] = g1.pt();
                value[2] = g1.recalculatedVtxR();
                value[3] = phi_cp1 > 0.f ? phi_cp1 : phi_cp1 + TMath::TwoPi();
                value[4] = eta_cp1;
                value[5] = g2.pt();
                value[6] = g2.recalculatedVtxR();
                value[7] = phi_cp2 > 0.f ? phi_cp2 : phi_cp2 + TMath::TwoPi();
                value[8] = eta_cp2;

                reinterpret_cast<THnSparseF*>(list_pair_ss->FindObject(Form("%s_%s", cut1.GetName(), cut2.GetName()))->FindObject(paircut.GetName())->FindObject("hs_conv_point_same"))->Fill(value);

              } // end of pair cut loop
            }   // end of g2 loop
          }     // end of g1 loop
        }       // end of cut2 loop
      }         // end of cut1 loop
    }           // end of collision loop
  }

  Partition<aod::EMReducedEvents> grouped_collisions = CentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < CentMax; // this goes to same event.

  void processMBMC(aod::EMReducedEvents const& collisions, MyV0Photons const& v0photons, MyMCV0Legs const& legs, aod::EMMCParticles const& mcparticles, aod::EMReducedMCEvents const& mccollisions)
  {
    fillsinglephoton(grouped_collisions, v0photons, perCollision_pcm, fProbePCMCuts, legs, mcparticles, mccollisions);
    SameEventPairing<PairType::kPCMPCM>(grouped_collisions, v0photons, v0photons, perCollision_pcm, perCollision_pcm, fTagPCMCuts, fProbePCMCuts, fPairCuts, legs, mcparticles, mccollisions);
  }

  PresliceUnsorted<aod::EMMCParticles> perMcCollision = aod::emmcparticle::emreducedmceventId;
  void processGen(soa::Join<aod::EMReducedEvents, aod::EMReducedMCEventLabels> const& collisions, aod::EMReducedMCEvents const&, aod::EMMCParticles const& mcparticles)
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

  void processDummy(aod::EMReducedEvents::iterator const& collision)
  {
    // do nothing
  }

  PROCESS_SWITCH(MaterialBudgetMC, processMBMC, "process material budget", false);
  PROCESS_SWITCH(MaterialBudgetMC, processGen, "process generated information", false);
  PROCESS_SWITCH(MaterialBudgetMC, processDummy, "Dummy function", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MaterialBudgetMC>(cfgc, TaskName{"material-budget-mc"})};
}
