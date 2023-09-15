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
// This code loops over photons and makes pairs for neutral mesons analyses.
//    Please write to: daiki.sekihata@cern.ch

#include "TString.h"
#include "Math/Vector4D.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/RecoDecay.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/PairUtilities.h"
#include "PWGEM/PhotonMeson/Utils/MCUtilities.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/Core/PHOSPhotonCut.h"
#include "PWGEM/PhotonMeson/Core/EMCPhotonCut.h"
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

struct Pi0EtaToGammaGammaMC {
  using MyMCV0Legs = soa::Join<aod::V0Legs, aod::EMMCParticleLabels>;

  // HistogramRegistry registry{"Pi0EtaToGammaGammaMC"};
  Configurable<float> maxY{"maxY", 0.9, "maximum rapidity for generated particles"};
  Configurable<std::string> fConfigPCMCuts{"cfgPCMCuts", "analysis,qc,nocut", "Comma separated list of V0 photon cuts"};
  Configurable<std::string> fConfigPairCuts{"cfgPairCuts", "nocut,asym08", "Comma separated list of pair cuts"};

  OutputObj<THashList> fOutputEvent{"Event"};
  OutputObj<THashList> fOutputPair{"Pair"}; // 2-photon pair
  OutputObj<THashList> fOutputGen{"Generated"};
  THashList* fMainList = new THashList();

  std::vector<V0PhotonCut> fPCMCuts;
  std::vector<PHOSPhotonCut> fPHOSCuts;
  std::vector<EMCPhotonCut> fEMCCuts;
  std::vector<PairCut> fPairCuts;

  std::vector<std::string> fPairNames;
  void init(InitContext& context)
  {
    if (context.mOptions.get<bool>("processPCMPCM")) {
      fPairNames.push_back("PCMPCM");
    }
    if (context.mOptions.get<bool>("processPHOSPHOS")) {
      fPairNames.push_back("PHOSPHOS");
    }
    if (context.mOptions.get<bool>("processEMCEMC")) {
      fPairNames.push_back("EMCEMC");
    }
    if (context.mOptions.get<bool>("processPCMPHOS")) {
      fPairNames.push_back("PCMPHOS");
    }
    if (context.mOptions.get<bool>("processPCMEMC")) {
      fPairNames.push_back("PCMEMC");
    }
    if (context.mOptions.get<bool>("processPHOSEMC")) {
      fPairNames.push_back("PHOSEMC");
    }

    DefinePCMCuts();
    DefinePairCuts();
    addhistograms();

    fOutputEvent.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("Event")));
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

        if ((pairname == "PCMPCM" || pairname == "PHOSPHOS" || pairname == "EMCEMC") && (cutname1 != cutname2))
          continue;

        THashList* list_pair_subsys = reinterpret_cast<THashList*>(list_pair->FindObject(pairname.data()));
        std::string photon_cut_name = cutname1 + "_" + cutname2;
        o2::aod::emphotonhistograms::AddHistClass(list_pair_subsys, photon_cut_name.data());
        THashList* list_pair_subsys_photoncut = reinterpret_cast<THashList*>(list_pair_subsys->FindObject(photon_cut_name.data()));

        for (auto& cut3 : cuts3) {
          std::string pair_cut_name = cut3.GetName();
          o2::aod::emphotonhistograms::AddHistClass(list_pair_subsys_photoncut, pair_cut_name.data());
          THashList* list_pair_subsys_paircut = reinterpret_cast<THashList*>(list_pair_subsys_photoncut->FindObject(pair_cut_name.data()));
          o2::aod::emphotonhistograms::DefineHistograms(list_pair_subsys_paircut, "gammagamma_mass_pt_mc");
        } // end of cut3 loop
      }   // end of cut2 loop
    }     // end of cut1 loop
  }

  static constexpr std::string_view pairnames[6] = {"PCMPCM", "PHOSPHOS", "EMCEMC", "PCMPHOS", "PCMEMC", "PHOSEMC"};
  static constexpr std::string_view parnames[3] = {"Gamma", "Pi0", "Eta"};
  void addhistograms()
  {
    fMainList->SetOwner(true);
    fMainList->SetName("fMainList");

    // create sub lists first.
    o2::aod::emphotonhistograms::AddHistClass(fMainList, "Event");
    THashList* list_ev = reinterpret_cast<THashList*>(fMainList->FindObject("Event"));

    o2::aod::emphotonhistograms::AddHistClass(fMainList, "Pair");
    THashList* list_pair = reinterpret_cast<THashList*>(fMainList->FindObject("Pair"));

    o2::aod::emphotonhistograms::AddHistClass(fMainList, "Generated");
    THashList* list_gen = reinterpret_cast<THashList*>(fMainList->FindObject("Generated"));
    o2::aod::emphotonhistograms::DefineHistograms(list_gen, "Generated", "Pi0Eta");

    for (auto& pairname : fPairNames) {
      LOGF(info, "Enabled pairs = %s", pairname.data());

      o2::aod::emphotonhistograms::AddHistClass(list_ev, pairname.data());
      THashList* list_ev_pair = reinterpret_cast<THashList*>(list_ev->FindObject(pairname.data()));
      o2::aod::emphotonhistograms::DefineHistograms(list_ev_pair, "Event");

      o2::aod::emphotonhistograms::AddHistClass(list_pair, pairname.data());

      if (pairname == "PCMPCM") {
        add_pair_histograms(list_pair, pairname, fPCMCuts, fPCMCuts, fPairCuts);
      }
      if (pairname == "PHOSPHOS") {
        add_pair_histograms(list_pair, pairname, fPHOSCuts, fPHOSCuts, fPairCuts);
      }
      if (pairname == "EMCEMC") {
        add_pair_histograms(list_pair, pairname, fEMCCuts, fEMCCuts, fPairCuts);
      }
      if (pairname == "PCMPHOS") {
        add_pair_histograms(list_pair, pairname, fPCMCuts, fPHOSCuts, fPairCuts);
      }
      if (pairname == "PCMEMC") {
        add_pair_histograms(list_pair, pairname, fPCMCuts, fEMCCuts, fPairCuts);
      }
      if (pairname == "PHOSEMC") {
        add_pair_histograms(list_pair, pairname, fPHOSCuts, fEMCCuts, fPairCuts);
      }

    } // end of pair name loop
  }

  void DefinePCMCuts()
  {
    TString cutNamesStr = fConfigPCMCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        const char* cutname = objArray->At(icut)->GetName();
        LOGF(info, "add cut : %s", cutname);
        fPCMCuts.push_back(*pcmcuts::GetCut(cutname));
      }
    }
    LOGF(info, "Number of PCM cuts = %d", fPCMCuts.size());
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

  template <PairType pairtype, typename TG1, typename TG2, typename TCut1, typename TCut2>
  bool IsSelectedPair(TG1 const& g1, TG2 const& g2, TCut1 const& cut1, TCut2 const& cut2)
  {
    bool is_selected_pair = false;
    if constexpr (pairtype == PairType::kPCMPCM) {
      is_selected_pair = o2::aod::photonpair::IsSelectedPair<MyMCV0Legs, MyMCV0Legs>(g1, g2, cut1, cut2);
    } else if constexpr (pairtype == PairType::kPHOSPHOS) {
      is_selected_pair = o2::aod::photonpair::IsSelectedPair<int, int>(g1, g2, cut1, cut2); // dummy, because track matching is not ready.
    } else if constexpr (pairtype == PairType::kEMCEMC) {
      is_selected_pair = o2::aod::photonpair::IsSelectedPair<aod::SkimEMCMTs, aod::SkimEMCMTs>(g1, g2, cut1, cut2);
    } else if constexpr (pairtype == PairType::kPCMPHOS) {
      is_selected_pair = o2::aod::photonpair::IsSelectedPair<MyMCV0Legs, int>(g1, g2, cut1, cut2);
    } else if constexpr (pairtype == PairType::kPCMEMC) {
      is_selected_pair = o2::aod::photonpair::IsSelectedPair<MyMCV0Legs, aod::SkimEMCMTs>(g1, g2, cut1, cut2);
    } else if constexpr (pairtype == PairType::kPHOSEMC) {
      is_selected_pair = o2::aod::photonpair::IsSelectedPair<int, aod::SkimEMCMTs>(g1, g2, cut1, cut2);
    } else {
      is_selected_pair = true;
    }
    return is_selected_pair;
  }

  Preslice<MyV0Photons> perCollision_pcm = aod::v0photon::emreducedeventId;

  template <PairType pairtype, typename TEvents, typename TPhotons1, typename TPhotons2, typename TPreslice1, typename TPreslice2, typename TCuts1, typename TCuts2, typename TPairCuts, typename TV0Legs, typename TMCParticles, typename TMCEvents>
  void TruePairing(TEvents const& collisions, TPhotons1 const& photons1, TPhotons2 const& photons2, TPreslice1 const& perCollision1, TPreslice2 const& perCollision2, TCuts1 const& cuts1, TCuts2 const& cuts2, TPairCuts const& paircuts, TV0Legs const& v0legs, TMCParticles const& mcparticles, TMCEvents const& mcevents)
  {
    THashList* list_ev_pair = static_cast<THashList*>(fMainList->FindObject("Event")->FindObject(pairnames[pairtype].data()));
    THashList* list_pair_ss = static_cast<THashList*>(fMainList->FindObject("Pair")->FindObject(pairnames[pairtype].data()));

    for (auto& collision : collisions) {

      if ((pairtype == kPHOSPHOS || pairtype == kPCMPHOS) && !collision.isPHOSCPVreadout()) {
        continue;
      }
      if ((pairtype == kEMCEMC || pairtype == kPCMEMC) && !collision.isEMCreadout()) {
        continue;
      }

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

      int pi0id = -1;
      int etaid = -1;
      if (pairtype == PairType::kPCMPCM || pairtype == PairType::kPHOSPHOS || pairtype == PairType::kEMCEMC) {
        for (auto& cut : cuts1) {
          for (auto& paircut : paircuts) {
            for (auto& [g1, g2] : combinations(CombinationsStrictlyUpperIndexPolicy(photons1_coll, photons2_coll))) {
              if (!IsSelectedPair<pairtype>(g1, g2, cut, cut)) {
                continue;
              }
              if (!paircut.IsSelected(g1, g2)) {
                continue;
              }

              if constexpr (pairtype == PairType::kPCMPCM) { // check 2 legs
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
                pi0id = FindCommonMotherFrom2Prongs(g1mc, g2mc, 22, 22, 111, mcparticles);
                etaid = FindCommonMotherFrom2Prongs(g1mc, g2mc, 22, 22, 221, mcparticles);
              }

              ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
              ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
              ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
              if (abs(v12.Rapidity()) > maxY) {
                continue;
              }

              if (pi0id > 0) {
                auto pi0mc = mcparticles.iteratorAt(pi0id);
                if (IsPhysicalPrimary(pi0mc.emreducedmcevent(), pi0mc, mcparticles)) {
                  reinterpret_cast<TH2F*>(list_pair_ss->FindObject(Form("%s_%s", cut.GetName(), cut.GetName()))->FindObject(paircut.GetName())->FindObject("hMggPt_Pi0_Primary"))->Fill(v12.M(), v12.Pt());
                } else if (IsFromWD(pi0mc.emreducedmcevent(), pi0mc, mcparticles)) {
                  reinterpret_cast<TH2F*>(list_pair_ss->FindObject(Form("%s_%s", cut.GetName(), cut.GetName()))->FindObject(paircut.GetName())->FindObject("hMggPt_Pi0_FromWD"))->Fill(v12.M(), v12.Pt());
                }
              }
              if (etaid > 0) {
                reinterpret_cast<TH2F*>(list_pair_ss->FindObject(Form("%s_%s", cut.GetName(), cut.GetName()))->FindObject(paircut.GetName())->FindObject("hMggPt_Eta_Primary"))->Fill(v12.M(), v12.Pt());
              }
            } // end of combination
          }   // end of paircutloop
        }     // end of cut loop

      } else { // different subsystem pairs
        for (auto& cut1 : cuts1) {
          for (auto& cut2 : cuts2) {
            for (auto& paircut : paircuts) {
              for (auto& [g1, g2] : combinations(CombinationsFullIndexPolicy(photons1_coll, photons2_coll))) {
                if (!IsSelectedPair<pairtype>(g1, g2, cut1, cut2)) {
                  continue;
                }
                if (!paircut.IsSelected(g1, g2)) {
                  continue;
                }

                if constexpr (pairtype == PairType::kPCMPHOS || pairtype == PairType::kPCMEMC) {
                  auto pos = g1.template posTrack_as<MyMCV0Legs>();
                  auto ele = g1.template negTrack_as<MyMCV0Legs>();
                  if constexpr (pairtype == PairType::kPCMPHOS) {
                    if (o2::aod::photonpair::DoesV0LegMatchWithCluster(pos, g2, 0.02, 0.4, 0.2) || o2::aod::photonpair::DoesV0LegMatchWithCluster(ele, g2, 0.02, 0.4, 0.2)) {
                      continue;
                    }
                  } else if constexpr (pairtype == PairType::kPCMEMC) {
                    if (o2::aod::photonpair::DoesV0LegMatchWithCluster(pos, g2, 0.02, 0.4, 0.5) || o2::aod::photonpair::DoesV0LegMatchWithCluster(ele, g2, 0.02, 0.4, 0.5)) {
                      continue;
                    }
                  }
                }

                ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
                ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
                ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
                if (abs(v12.Rapidity()) > maxY) {
                  continue;
                }
                reinterpret_cast<TH2F*>(list_pair_ss->FindObject(Form("%s_%s", cut1.GetName(), cut2.GetName()))->FindObject(paircut.GetName())->FindObject("hMggPt_Pi0_Primary"))->Fill(v12.M(), v12.Pt());
                reinterpret_cast<TH2F*>(list_pair_ss->FindObject(Form("%s_%s", cut1.GetName(), cut2.GetName()))->FindObject(paircut.GetName())->FindObject("hMggPt_Eta_Primary"))->Fill(v12.M(), v12.Pt());

              } // end of combination
            }   // end of paircutloop
          }     // end of cut2 loop
        }       // end of cut1 loop
      }

    } // end of collision loop
  }

  using MyCollisions = soa::Join<aod::EMReducedEvents, aod::EMReducedMCEventLabels>;

  void processPCMPCM(MyCollisions const& collisions, MyV0Photons const& v0photons, MyMCV0Legs const& v0legs, aod::EMMCParticles const& mcparticles, aod::EMReducedMCEvents const& mccollisions)
  {
    TruePairing<PairType::kPCMPCM>(collisions, v0photons, v0photons, perCollision_pcm, perCollision_pcm, fPCMCuts, fPCMCuts, fPairCuts, v0legs, mcparticles, mccollisions);
  }
  void processPHOSPHOS(MyCollisions const& collisions) {}
  void processEMCEMC(MyCollisions const& collisions) {}
  void processPCMPHOS(MyCollisions const& collisions) {}
  void processPCMEMC(MyCollisions const& collisions) {}
  void processPHOSEMC(MyCollisions const& collisions) {}

  PresliceUnsorted<aod::EMMCParticles> perMcCollision = aod::emmcparticle::emreducedmceventId;
  PresliceUnsorted<MyCollisions> rec_perMcCollision = aod::emmceventlabel::emreducedmceventId;
  void processGen(MyCollisions const& collisions, aod::EMReducedMCEvents const& mccollisions, aod::EMMCParticles const& mcparticles)
  {
    // loop over mc stack and fill histograms for pure MC truth signals
    // all MC tracks which belong to the MC event corresponding to the current reconstructed event
    for (auto& mccollision : mccollisions) {
      auto collision_per_mccoll = collisions.sliceBy(rec_perMcCollision, mccollision.globalIndex());
      int nrec_per_mc = collision_per_mccoll.size();
      reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hNrecPerMCCollision"))->Fill(nrec_per_mc); // all
    }

    for (auto& collision : collisions) {

      auto mccollision = collision.emreducedmcevent();
      reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hZvtx_before"))->Fill(mccollision.posZ());
      reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hCollisionCounter"))->Fill(1.0); // all
      if (!collision.sel8()) {
        continue;
      }
      reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hCollisionCounter"))->Fill(2.0); // FT0VX i.e. FT0and

      if (collision.numContrib() < 0.5) {
        continue;
      }
      reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hCollisionCounter"))->Fill(3.0); // Ncontrib > 0

      if (abs(collision.posZ()) > 10.0) {
        continue;
      }
      reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hZvtx_after"))->Fill(mccollision.posZ());
      reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hCollisionCounter"))->Fill(4.0); // |Zvtx| < 10 cm

      auto mctracks_coll = mcparticles.sliceBy(perMcCollision, mccollision.globalIndex());
      for (auto& mctrack : mctracks_coll) {
        if (abs(mctrack.y()) > maxY) {
          continue;
        }
        int pdg = mctrack.pdgCode();

        if (abs(pdg) == 111 && IsPhysicalPrimary(mctrack.emreducedmcevent(), mctrack, mcparticles)) {
          reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hPt_Pi0"))->Fill(mctrack.pt());
          reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hY_Pi0"))->Fill(mctrack.y());
          reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hPhi_Pi0"))->Fill(mctrack.phi());
        }
        if (abs(pdg) == 221 && IsPhysicalPrimary(mctrack.emreducedmcevent(), mctrack, mcparticles)) {
          reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hPt_Eta"))->Fill(mctrack.pt());
          reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hY_Eta"))->Fill(mctrack.y());
          reinterpret_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hPhi_Eta"))->Fill(mctrack.phi());
        }

      } // end of mc track loop
    }   // end of collision loop
  }

  void processDummy(aod::EMReducedEvents::iterator const& collision) {}

  PROCESS_SWITCH(Pi0EtaToGammaGammaMC, processPCMPCM, "true pairing PCM-PCM", false);
  PROCESS_SWITCH(Pi0EtaToGammaGammaMC, processPHOSPHOS, "true pairing PHOS-PHOS", false);
  PROCESS_SWITCH(Pi0EtaToGammaGammaMC, processEMCEMC, "true pairing EMC-EMC", false);
  PROCESS_SWITCH(Pi0EtaToGammaGammaMC, processPCMPHOS, "true pairing PCM-PHOS", false);
  PROCESS_SWITCH(Pi0EtaToGammaGammaMC, processPCMEMC, "true pairing PCM-EMC", false);
  PROCESS_SWITCH(Pi0EtaToGammaGammaMC, processPHOSEMC, "true pairing PHOS-EMC", false);
  PROCESS_SWITCH(Pi0EtaToGammaGammaMC, processGen, "process generated information", false);
  PROCESS_SWITCH(Pi0EtaToGammaGammaMC, processDummy, "Dummy function", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<Pi0EtaToGammaGammaMC>(cfgc, TaskName{"pi0eta-to-gammagamma-mc"})};
}
