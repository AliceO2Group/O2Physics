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
// This code runs loop over photons with PCM and PHOS for direct photon analysis.
//    Please write to: daiki.sekihata@cern.ch

#include <cstring>
#include <iterator>

#include "TString.h"
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include "Math/LorentzRotation.h"
#include "Math/Rotation3D.h"
#include "Math/AxisAngle.h"
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
#include "PWGEM/PhotonMeson/Utils/PairUtilities.h"
#include "PWGEM/PhotonMeson/Utils/MCUtilities.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
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

using MyV0Photons = soa::Join<aod::V0PhotonsKF, aod::V0Recalculation>;
using MyV0Photon = MyV0Photons::iterator;

struct TaggingPi0MC {
  using MyMCV0Legs = soa::Join<aod::V0Legs, aod::EMMCParticleLabels>;

  Configurable<float> maxY{"maxY", 0.9, "maximum rapidity for reconstructed particles"};
  Configurable<std::string> fConfigPCMCuts{"cfgPCMCuts", "analysis1690", "Comma separated list of V0 photon cuts"};
  Configurable<std::string> fConfigPCMibwCuts{"cfgPCMibwCuts", "wwire_ib", "Comma separated list of V0 photon cuts"};
  Configurable<std::string> fConfigPHOSCuts{"cfgPHOSCuts", "test02,test03", "Comma separated list of PHOS photon cuts"};
  Configurable<std::string> fConfigEMCCuts{"fConfigEMCCuts", "standard", "Comma separated list of EMCal photon cuts"};
  Configurable<std::string> fConfigPairCuts{"cfgPairCuts", "nocut", "Comma separated list of pair cuts"};

  // Configurable for EMCal cuts
  Configurable<float> EMC_minTime{"EMC_minTime", -20., "Minimum cluster time for EMCal time cut"};
  Configurable<float> EMC_maxTime{"EMC_maxTime", +25., "Maximum cluster time for EMCal time cut"};
  Configurable<float> EMC_minM02{"EMC_minM02", 0.1, "Minimum M02 for EMCal M02 cut"};
  Configurable<float> EMC_maxM02{"EMC_maxM02", 0.7, "Maximum M02 for EMCal M02 cut"};
  Configurable<float> EMC_minE{"EMC_minE", 0.7, "Minimum cluster energy for EMCal energy cut"};
  Configurable<int> EMC_minNCell{"EMC_minNCell", 1, "Minimum number of cells per cluster for EMCal NCell cut"};
  Configurable<std::vector<float>> EMC_TM_Eta{"EMC_TM_Eta", {0.01f, 4.07f, -2.5f}, "|eta| <= [0]+(pT+[1])^[2] for EMCal track matching"};
  Configurable<std::vector<float>> EMC_TM_Phi{"EMC_TM_Phi", {0.015f, 3.65f, -2.f}, "|phi| <= [0]+(pT+[1])^[2] for EMCal track matching"};
  Configurable<float> EMC_Eoverp{"EMC_Eoverp", 1.75, "Minimum cluster energy over track momentum for EMCal track matching"};
  Configurable<bool> EMC_UseExoticCut{"EMC_UseExoticCut", true, "FLag to use the EMCal exotic cluster cut"};

  OutputObj<THashList> fOutputEvent{"Event"};
  OutputObj<THashList> fOutputPair{"Pair"}; // 2-photon pair
  OutputObj<THashList> fOutputPCM{"PCM"};   // v0-photon
  THashList* fMainList = new THashList();

  std::vector<V0PhotonCut> fPCMCuts;
  std::vector<V0PhotonCut> fPCMibwCuts;
  std::vector<PHOSPhotonCut> fPHOSCuts;
  std::vector<EMCPhotonCut> fEMCCuts;
  std::vector<PairCut> fPairCuts;

  std::vector<std::string> fPairNames;
  void init(InitContext& context)
  {
    if (context.mOptions.get<bool>("processPCMPCMibw")) {
      fPairNames.push_back("PCMPCMibw");
    }
    if (context.mOptions.get<bool>("processPCMPHOS")) {
      fPairNames.push_back("PCMPHOS");
    }
    if (context.mOptions.get<bool>("processPCMEMC")) {
      fPairNames.push_back("PCMEMC");
    }

    DefinePCMCuts();
    DefinePCMibwCuts();
    DefinePHOSCuts();
    DefineEMCCuts();
    DefinePairCuts();
    addhistograms();

    fOutputEvent.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("Event")));
    fOutputPair.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("Pair")));
    fOutputPCM.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("PCM")));
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
          o2::aod::emphotonhistograms::DefineHistograms(list_pair_subsys_paircut, "tagging_pi0_mc", "pair");
        } // end of pair cut loop
      }   // end of cut2 loop
    }     // end of cut1 loop
  }

  static constexpr std::string_view pairnames[7] = {"PCMPCM", "PHOSPHOS", "EMCEMC", "PCMPHOS", "PCMEMC", "PHOSEMC", "PCMPCMibw"};
  void addhistograms()
  {
    fMainList->SetOwner(true);
    fMainList->SetName("fMainList");

    // create sub lists first.
    o2::aod::emphotonhistograms::AddHistClass(fMainList, "Event");
    THashList* list_ev = reinterpret_cast<THashList*>(fMainList->FindObject("Event"));

    o2::aod::emphotonhistograms::AddHistClass(fMainList, "PCM");
    THashList* list_pcm = reinterpret_cast<THashList*>(fMainList->FindObject("PCM"));
    for (auto& cut : fPCMCuts) {
      THashList* list_pcm_cut = o2::aod::emphotonhistograms::AddHistClass(list_pcm, cut.GetName());
      o2::aod::emphotonhistograms::DefineHistograms(list_pcm_cut, "tagging_pi0_mc", "pcm");
    }

    o2::aod::emphotonhistograms::AddHistClass(fMainList, "Pair");
    THashList* list_pair = reinterpret_cast<THashList*>(fMainList->FindObject("Pair"));

    for (auto& pairname : fPairNames) {
      LOGF(info, "Enabled pairs = %s", pairname.data());

      o2::aod::emphotonhistograms::AddHistClass(list_ev, pairname.data());
      THashList* list_ev_pair = reinterpret_cast<THashList*>(list_ev->FindObject(pairname.data()));
      o2::aod::emphotonhistograms::DefineHistograms(list_ev_pair, "Event");

      o2::aod::emphotonhistograms::AddHistClass(list_pair, pairname.data());

      if (pairname == "PCMPCMibw") {
        add_pair_histograms(list_pair, pairname, fPCMCuts, fPCMibwCuts, fPairCuts);
      }
      if (pairname == "PCMPHOS") {
        add_pair_histograms(list_pair, pairname, fPCMCuts, fPHOSCuts, fPairCuts);
      }
      if (pairname == "PCMEMC") {
        add_pair_histograms(list_pair, pairname, fPCMCuts, fEMCCuts, fPairCuts);
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
  void DefinePCMibwCuts()
  {
    TString cutNamesStr = fConfigPCMibwCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        const char* cutname = objArray->At(icut)->GetName();
        LOGF(info, "add cut : %s", cutname);
        fPCMibwCuts.push_back(*pcmcuts::GetCut(cutname));
      }
    }
    LOGF(info, "Number of PCMibw cuts = %d", fPCMibwCuts.size());
  }
  void DefinePHOSCuts()
  {
    TString cutNamesStr = fConfigPHOSCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        const char* cutname = objArray->At(icut)->GetName();
        LOGF(info, "add cut : %s", cutname);
        fPHOSCuts.push_back(*phoscuts::GetCut(cutname));
      }
    }
    LOGF(info, "Number of PHOS cuts = %d", fPHOSCuts.size());
  }

  void DefineEMCCuts()
  {
    const float a = EMC_TM_Eta->at(0);
    const float b = EMC_TM_Eta->at(1);
    const float c = EMC_TM_Eta->at(2);

    const float d = EMC_TM_Phi->at(0);
    const float e = EMC_TM_Phi->at(1);
    const float f = EMC_TM_Phi->at(2);
    LOGF(info, "EMCal track matching parameters : a = %f, b = %f, c = %f, d = %f, e = %f, f = %f", a, b, c, d, e, f);

    TString cutNamesStr = fConfigEMCCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        const char* cutname = objArray->At(icut)->GetName();
        LOGF(info, "add cut : %s", cutname);
        if (std::strcmp(cutname, "custom") == 0) {
          EMCPhotonCut* custom_cut = new EMCPhotonCut(cutname, cutname);
          custom_cut->SetMinE(EMC_minE);
          custom_cut->SetMinNCell(EMC_minNCell);
          custom_cut->SetM02Range(EMC_minM02, EMC_maxM02);
          custom_cut->SetTimeRange(EMC_minTime, EMC_maxTime);

          custom_cut->SetTrackMatchingEta([&a, &b, &c](float pT) {
            return a + pow(pT + b, c);
          });
          custom_cut->SetTrackMatchingPhi([&d, &e, &f](float pT) {
            return d + pow(pT + e, f);
          });

          custom_cut->SetMinEoverP(EMC_Eoverp);
          custom_cut->SetUseExoticCut(EMC_UseExoticCut);
          fEMCCuts.push_back(*custom_cut);
        } else {
          fEMCCuts.push_back(*aod::emccuts::GetCut(cutname));
        }
      }
    }
    LOGF(info, "Number of EMCal cuts = %d", fEMCCuts.size());
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
    if constexpr (pairtype == PairType::kPCMPCMibw) {
      is_selected_pair = o2::aod::photonpair::IsSelectedPair<MyMCV0Legs, MyMCV0Legs>(g1, g2, cut1, cut2);
    } else if constexpr (pairtype == PairType::kPCMPHOS) {
      is_selected_pair = o2::aod::photonpair::IsSelectedPair<MyMCV0Legs, int>(g1, g2, cut1, cut2);
    } else if constexpr (pairtype == PairType::kPCMEMC) {
      is_selected_pair = o2::aod::photonpair::IsSelectedPair<MyMCV0Legs, aod::SkimEMCMTs>(g1, g2, cut1, cut2);
    } else {
      is_selected_pair = true;
    }
    return is_selected_pair;
  }

  template <PairType pairtype, typename TEvents, typename TPhotons1, typename TPhotons2, typename TPreslice1, typename TPreslice2, typename TCuts1, typename TCuts2, typename TPairCuts, typename TLegs, typename TMCParticles, typename TMCEvents>
  void TruePairing(TEvents const& collisions, TPhotons1 const& photons1, TPhotons2 const& photons2, TPreslice1 const& perCollision1, TPreslice2 const& perCollision2, TCuts1 const& cuts1, TCuts2 const& cuts2, TPairCuts const& paircuts, TLegs const& legs, TMCParticles const& mcparticles, TMCEvents const& mcevents)
  {
    THashList* list_ev_pair = static_cast<THashList*>(fMainList->FindObject("Event")->FindObject(pairnames[pairtype].data()));
    THashList* list_pair_ss = static_cast<THashList*>(fMainList->FindObject("Pair")->FindObject(pairnames[pairtype].data()));
    THashList* list_pcm = static_cast<THashList*>(fMainList->FindObject("PCM"));

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

      for (auto& cut1 : cuts1) {
        for (auto& g1 : photons1_coll) {

          if (!cut1.template IsSelected<MyMCV0Legs>(g1)) {
            continue;
          }
          if (abs(g1.eta()) > maxY) { // photon is massless particle. rapidity = pseudo-rapidity
            continue;
          }

          auto pos1 = g1.template posTrack_as<MyMCV0Legs>();
          auto ele1 = g1.template negTrack_as<MyMCV0Legs>();
          auto pos1mc = pos1.template emmcparticle_as<aod::EMMCParticles>();
          auto ele1mc = ele1.template emmcparticle_as<aod::EMMCParticles>();

          int photonid1 = FindCommonMotherFrom2Prongs(pos1mc, ele1mc, -11, 11, 22, mcparticles);
          if (photonid1 < 0) { // check swap, true electron is reconstructed as positron and vice versa.
            photonid1 = FindCommonMotherFrom2Prongs(pos1mc, ele1mc, 11, -11, 22, mcparticles);
          }

          if (photonid1 > 0) {
            auto mcphoton1 = mcparticles.iteratorAt(photonid1);
            int pi0id1 = IsXFromY(mcphoton1, mcparticles, 22, 111);
            if (pi0id1 > 0) {
              auto mcpi01 = mcparticles.iteratorAt(pi0id1);
              if (IsPhysicalPrimary(mcpi01.emreducedmcevent(), mcpi01, mcparticles)) {
                reinterpret_cast<TH1F*>(list_pcm->FindObject(Form("%s", cut1.GetName()))->FindObject("hPt_v0photon_Pi0_Primary"))->Fill(g1.pt());
              } else if (IsFromWD(mcpi01.emreducedmcevent(), mcpi01, mcparticles)) {
                reinterpret_cast<TH1F*>(list_pcm->FindObject(Form("%s", cut1.GetName()))->FindObject("hPt_v0photon_Pi0_FromWD"))->Fill(g1.pt());
              } else {
                reinterpret_cast<TH1F*>(list_pcm->FindObject(Form("%s", cut1.GetName()))->FindObject("hPt_v0photon_Pi0_hs"))->Fill(g1.pt());
              }
            }
          }

        } // end of pcm photon loop
      }   // end of cut loop

      for (auto& cut1 : cuts1) {
        for (auto& cut2 : cuts2) {
          for (auto& paircut : paircuts) {
            for (auto& [g1, g2] : combinations(CombinationsFullIndexPolicy(photons1_coll, photons2_coll))) {
              if (!IsSelectedPair<pairtype>(g1, g2, cut1, cut2)) {
                continue;
              }

              auto pos1 = g1.template posTrack_as<MyMCV0Legs>();
              auto ele1 = g1.template negTrack_as<MyMCV0Legs>();

              if constexpr (pairtype == PairType::kPCMPHOS || pairtype == PairType::kPCMEMC) {
                if constexpr (pairtype == PairType::kPCMPHOS) {
                  if (o2::aod::photonpair::DoesV0LegMatchWithCluster(pos1, g2, 0.02, 0.4, 0.2) || o2::aod::photonpair::DoesV0LegMatchWithCluster(ele1, g2, 0.02, 0.4, 0.2)) {
                    continue;
                  }
                } else if constexpr (pairtype == PairType::kPCMEMC) {
                  if (o2::aod::photonpair::DoesV0LegMatchWithCluster(pos1, g2, 0.02, 0.4, 0.5) || o2::aod::photonpair::DoesV0LegMatchWithCluster(ele1, g2, 0.02, 0.4, 0.5)) {
                    continue;
                  }
                }
              }

              int pi0id = -1;
              if constexpr (pairtype == PairType::kPCMPCMibw) {

                auto pos1mc = pos1.template emmcparticle_as<aod::EMMCParticles>();
                auto ele1mc = ele1.template emmcparticle_as<aod::EMMCParticles>();

                int photonid1 = FindCommonMotherFrom2Prongs(pos1mc, ele1mc, -11, 11, 22, mcparticles);
                if (photonid1 < 0) { // check swap, true electron is reconstructed as positron and vice versa.
                  photonid1 = FindCommonMotherFrom2Prongs(pos1mc, ele1mc, 11, -11, 22, mcparticles);
                }

                auto pos2 = g2.template posTrack_as<MyMCV0Legs>();
                auto ele2 = g2.template negTrack_as<MyMCV0Legs>();

                auto pos2mc = pos2.template emmcparticle_as<aod::EMMCParticles>();
                auto ele2mc = ele2.template emmcparticle_as<aod::EMMCParticles>();

                int photonid2 = FindCommonMotherFrom2Prongs(pos2mc, ele2mc, -11, 11, 22, mcparticles);
                if (photonid2 < 0) { // check swap, true electron is reconstructed as positron and vice versa.
                  photonid2 = FindCommonMotherFrom2Prongs(pos2mc, ele2mc, 11, -11, 22, mcparticles);
                }

                if (photonid1 < 0 || photonid2 < 0) {
                  continue;
                }

                auto g1mc = mcparticles.iteratorAt(photonid1);
                auto g2mc = mcparticles.iteratorAt(photonid2);
                pi0id = FindCommonMotherFrom2Prongs(g1mc, g2mc, 22, 22, 111, mcparticles);
              }

              ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.); // pcm
              ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.); // phos or pcmibw
              ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
              // if (abs(v12.Rapidity()) > maxY) {
              //   continue;
              // }

              if (pi0id > 0) {
                auto mcpi0 = mcparticles.iteratorAt(pi0id);
                if (IsPhysicalPrimary(mcpi0.emreducedmcevent(), mcpi0, mcparticles)) {
                  reinterpret_cast<TH2F*>(list_pair_ss->FindObject(Form("%s_%s", cut1.GetName(), cut2.GetName()))->FindObject(paircut.GetName())->FindObject("hMggPt_Pi0_Primary"))->Fill(v12.M(), v1.Pt());
                } else if (IsFromWD(mcpi0.emreducedmcevent(), mcpi0, mcparticles)) {
                  reinterpret_cast<TH2F*>(list_pair_ss->FindObject(Form("%s_%s", cut1.GetName(), cut2.GetName()))->FindObject(paircut.GetName())->FindObject("hMggPt_Pi0_FromWD"))->Fill(v12.M(), v1.Pt());
                } else {
                  reinterpret_cast<TH2F*>(list_pair_ss->FindObject(Form("%s_%s", cut1.GetName(), cut2.GetName()))->FindObject(paircut.GetName())->FindObject("hMggPt_Pi0_hs"))->Fill(v12.M(), v1.Pt());
                }
              }
            } // end of combination
          }   // end of pair cut loop
        }     // end of cut2 loop
      }       // end of cut1 loop
    }         // end of collision loop
  }

  Preslice<MyV0Photons> perCollision_pcm = aod::v0photon::collisionId;
  Preslice<aod::PHOSClusters> perCollision_phos = aod::skimmedcluster::collisionId;
  Preslice<aod::SkimEMCClusters> perCollision_emc = aod::skimmedcluster::collisionId;

  void processPCMPCMibw(aod::EMReducedEvents const& collisions, MyV0Photons const& v0photons, MyMCV0Legs const& legs, aod::EMMCParticles const& mcparticles, aod::EMReducedMCEvents const& mccollisions)
  {
    TruePairing<PairType::kPCMPCMibw>(collisions, v0photons, v0photons, perCollision_pcm, perCollision_pcm, fPCMCuts, fPCMibwCuts, fPairCuts, legs, mcparticles, mccollisions);
  }

  void processPCMPHOS(aod::EMReducedEvents const& collisions, MyV0Photons const& v0photons, aod::PHOSClusters const& phosclusters, MyMCV0Legs const& legs, aod::EMMCParticles const& mcparticles, aod::EMReducedMCEvents const& mccollisions)
  {
    TruePairing<PairType::kPCMPHOS>(collisions, v0photons, phosclusters, perCollision_pcm, perCollision_phos, fPCMCuts, fPHOSCuts, fPairCuts, legs, mcparticles, mccollisions);
  }

  void processPCMEMC(aod::EMReducedEvents const& collisions, MyV0Photons const& v0photons, aod::SkimEMCClusters const& emcclusters, MyMCV0Legs const& legs, aod::EMMCParticles const& mcparticles, aod::EMReducedMCEvents const& mccollisions)
  {
    TruePairing<PairType::kPCMEMC>(collisions, v0photons, emcclusters, perCollision_pcm, perCollision_emc, fPCMCuts, fEMCCuts, fPairCuts, legs, mcparticles, mccollisions);
  }

  void processDummy(aod::EMReducedEvents::iterator const& collision) {}

  PROCESS_SWITCH(TaggingPi0MC, processPCMPCMibw, "pairing PCM-PCMibw", false);
  PROCESS_SWITCH(TaggingPi0MC, processPCMPHOS, "pairing PCM-PHOS", false);
  PROCESS_SWITCH(TaggingPi0MC, processPCMEMC, "pairing PCM-EMCal", false);
  PROCESS_SWITCH(TaggingPi0MC, processDummy, "Dummy function", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TaggingPi0MC>(cfgc, TaskName{"tagging-pi0-mc"})};
}
