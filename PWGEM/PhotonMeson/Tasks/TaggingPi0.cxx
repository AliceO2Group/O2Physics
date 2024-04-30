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
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Common/Core/RecoDecay.h"
#include "PWGEM/PhotonMeson/Utils/PairUtilities.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/Core/DalitzEECut.h"
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
using namespace o2::aod::pwgem::photon;

using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent, aod::EMEventsNgPCM, aod::EMEventsNgPHOS, aod::EMEventsNgEMC, aod::EMEventsNee>;
using MyCollision = MyCollisions::iterator;

using MyV0Photons = soa::Join<aod::V0PhotonsKF, aod::V0KFEMEventIds>;
using MyV0Photon = MyV0Photons::iterator;

using MyDalitzEEs = soa::Join<aod::DalitzEEs, aod::DalitzEEEMEventIds>;
using MyDalitzEE = MyDalitzEEs::iterator;

using MyPrimaryElectrons = soa::Join<aod::EMPrimaryElectrons, aod::EMPrimaryElectronEMEventIds, aod::EMPrimaryElectronsPrefilterBit>;
using MyPrimaryElectron = MyPrimaryElectrons::iterator;

struct TaggingPi0 {

  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999, "max. centrality"};

  // Configurable<float> maxY{"maxY", 0.9, "maximum rapidity for reconstructed particles"};
  Configurable<std::string> fConfigPCMCuts{"cfgPCMCuts", "qc", "Comma separated list of V0 photon cuts"};
  Configurable<std::string> fConfigDalitzEECuts{"cfgDalitzEECuts", "mee_0_120_tpchadrejortofreq,mee_0_120_tpchadrejortofreq_lowB", "Comma separated list of Dalitz ee cuts"};
  Configurable<std::string> fConfigPHOSCuts{"cfgPHOSCuts", "test02,test03", "Comma separated list of PHOS photon cuts"};
  Configurable<std::string> fConfigEMCCuts{"cfgEMCCuts", "standard", "Comma separated list of EMCal photon cuts"};
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

  Configurable<std::string> fConfigEMEventCut{"cfgEMEventCut", "minbias", "em event cut"}; // only 1 event cut per wagon
  EMEventCut fEMEventCut;
  static constexpr std::string_view event_types[2] = {"before", "after"};

  OutputObj<THashList> fOutputEvent{"Event"};
  OutputObj<THashList> fOutputPair{"Pair"}; // 2-photon pair
  THashList* fMainList = new THashList();

  std::vector<V0PhotonCut> fPCMCuts;
  std::vector<DalitzEECut> fDalitzEECuts;
  std::vector<PHOSPhotonCut> fPHOSCuts;
  std::vector<EMCPhotonCut> fEMCCuts;
  std::vector<PairCut> fPairCuts;

  std::vector<std::string> fPairNames;
  void init(InitContext& context)
  {
    if (context.mOptions.get<bool>("processPCMDalitzEE")) {
      fPairNames.push_back("PCMDalitzEE");
    }
    if (context.mOptions.get<bool>("processPCMPHOS")) {
      fPairNames.push_back("PCMPHOS");
    }
    if (context.mOptions.get<bool>("processPCMEMC")) {
      fPairNames.push_back("PCMEMC");
    }

    DefinePCMCuts();
    DefineDalitzEECuts();
    DefinePHOSCuts();
    DefineEMCCuts();
    DefinePairCuts();
    addhistograms();
    TString ev_cut_name = fConfigEMEventCut.value;
    fEMEventCut = *eventcuts::GetCut(ev_cut_name.Data());

    fOutputEvent.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("Event")));
    fOutputPair.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("Pair")));
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
        o2::aod::pwgem::photon::histogram::AddHistClass(list_pair_subsys, photon_cut_name.data());
        THashList* list_pair_subsys_photoncut = reinterpret_cast<THashList*>(list_pair_subsys->FindObject(photon_cut_name.data()));

        for (auto& cut3 : cuts3) {
          std::string pair_cut_name = cut3.GetName();
          o2::aod::pwgem::photon::histogram::AddHistClass(list_pair_subsys_photoncut, pair_cut_name.data());
          THashList* list_pair_subsys_paircut = reinterpret_cast<THashList*>(list_pair_subsys_photoncut->FindObject(pair_cut_name.data()));
          o2::aod::pwgem::photon::histogram::DefineHistograms(list_pair_subsys_paircut, "tagging_pi0");
        } // end of cut3 loop
      }   // end of cut2 loop
    }     // end of cut1 loop
  }

  static constexpr std::string_view pairnames[9] = {"PCMPCM", "PHOSPHOS", "EMCEMC", "PCMPHOS", "PCMEMC", "PCMDalitzEE", "PCMDalitzMuMu", "PHOSEMC", "DalitzEEDalitzEE"};
  void addhistograms()
  {
    fMainList->SetOwner(true);
    fMainList->SetName("fMainList");

    // create sub lists first.
    o2::aod::pwgem::photon::histogram::AddHistClass(fMainList, "Event");
    THashList* list_ev = reinterpret_cast<THashList*>(fMainList->FindObject("Event"));

    o2::aod::pwgem::photon::histogram::AddHistClass(fMainList, "Pair");
    THashList* list_pair = reinterpret_cast<THashList*>(fMainList->FindObject("Pair"));

    for (auto& pairname : fPairNames) {
      LOGF(info, "Enabled pairs = %s", pairname.data());

      THashList* list_ev_pair = reinterpret_cast<THashList*>(o2::aod::pwgem::photon::histogram::AddHistClass(list_ev, pairname.data()));
      for (const auto& evtype : event_types) {
        THashList* list_ev_type = reinterpret_cast<THashList*>(o2::aod::pwgem::photon::histogram::AddHistClass(list_ev_pair, evtype.data()));
        o2::aod::pwgem::photon::histogram::DefineHistograms(list_ev_type, "Event", evtype.data());
      }

      o2::aod::pwgem::photon::histogram::AddHistClass(list_pair, pairname.data());

      if (pairname == "PCMPHOS") {
        add_pair_histograms(list_pair, pairname, fPCMCuts, fPHOSCuts, fPairCuts);
      }
      if (pairname == "PCMEMC") {
        add_pair_histograms(list_pair, pairname, fPCMCuts, fEMCCuts, fPairCuts);
      }
      if (pairname == "PCMDalitzEE") {
        add_pair_histograms(list_pair, pairname, fPCMCuts, fDalitzEECuts, fPairCuts);
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

  void DefineDalitzEECuts()
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
    LOGF(info, "Number of DalitzEE cuts = %d", fDalitzEECuts.size());
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
          fEMCCuts.push_back(*emccuts::GetCut(cutname));
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
    if constexpr (pairtype == PairType::kPCMPHOS) {
      is_selected_pair = o2::aod::photonpair::IsSelectedPair<aod::V0Legs, int>(g1, g2, cut1, cut2);
    } else if constexpr (pairtype == PairType::kPCMEMC) {
      is_selected_pair = o2::aod::photonpair::IsSelectedPair<aod::V0Legs, aod::SkimEMCMTs>(g1, g2, cut1, cut2);
    } else if constexpr (pairtype == PairType::kPCMDalitzEE) {
      is_selected_pair = o2::aod::photonpair::IsSelectedPair<aod::V0Legs, MyPrimaryElectrons>(g1, g2, cut1, cut2);
    } else {
      is_selected_pair = true;
    }
    return is_selected_pair;
  }

  template <PairType pairtype, typename TEvents, typename TPhotons1, typename TPhotons2, typename TPreslice1, typename TPreslice2, typename TCuts1, typename TCuts2, typename TPairCuts, typename TLegs, typename TEMPrimaryElectrons>
  void SameEventPairing(TEvents const& collisions, TPhotons1 const& photons1, TPhotons2 const& photons2, TPreslice1 const& perCollision1, TPreslice2 const& perCollision2, TCuts1 const& cuts1, TCuts2 const& cuts2, TPairCuts const& paircuts, TLegs const& /*legs*/, TEMPrimaryElectrons const& /*emprimaryelectrons*/)
  {
    THashList* list_ev_pair_before = static_cast<THashList*>(fMainList->FindObject("Event")->FindObject(pairnames[pairtype].data())->FindObject(event_types[0].data()));
    THashList* list_ev_pair_after = static_cast<THashList*>(fMainList->FindObject("Event")->FindObject(pairnames[pairtype].data())->FindObject(event_types[1].data()));
    THashList* list_pair_ss = static_cast<THashList*>(fMainList->FindObject("Pair")->FindObject(pairnames[pairtype].data()));

    for (auto& collision : collisions) {

      if ((pairtype == kPHOSPHOS || pairtype == kPCMPHOS) && !collision.alias_bit(triggerAliases::kTVXinPHOS)) {
        continue;
      }
      if ((pairtype == kEMCEMC || pairtype == kPCMEMC) && !collision.alias_bit(triggerAliases::kTVXinEMC)) {
        continue;
      }

      const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      o2::aod::pwgem::photon::histogram::FillHistClass<EMHistType::kEvent>(list_ev_pair_before, "", collision);
      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      o2::aod::pwgem::photon::histogram::FillHistClass<EMHistType::kEvent>(list_ev_pair_after, "", collision);
      reinterpret_cast<TH1F*>(list_ev_pair_before->FindObject("hCollisionCounter"))->Fill("accepted", 1.f);
      reinterpret_cast<TH1F*>(list_ev_pair_after->FindObject("hCollisionCounter"))->Fill("accepted", 1.f);

      auto photons1_coll = photons1.sliceBy(perCollision1, collision.globalIndex());
      auto photons2_coll = photons2.sliceBy(perCollision2, collision.globalIndex());

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
                auto pos = g1.template posTrack_as<aod::V0Legs>();
                auto ele = g1.template negTrack_as<aod::V0Legs>();
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
              ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.); // pcm
              ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.); // phos or emc or dalitzee
              if constexpr (pairtype == PairType::kPCMDalitzEE) {
                v2.SetM(g2.mass());
              }
              ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
              // if (abs(v12.Rapidity()) > maxY) {
              //   continue;
              // }
              reinterpret_cast<TH2F*>(list_pair_ss->FindObject(Form("%s_%s", cut1.GetName(), cut2.GetName()))->FindObject(paircut.GetName())->FindObject("hMggPt_Same"))->Fill(v12.M(), v1.Pt());
            } // end of combination
          }   // end of pair cut loop
        }     // end of cut2 loop
      }       // end of cut1 loop
    }         // end of collision loop
  }

  Configurable<int> ndepth{"ndepth", 10, "depth for event mixing"};
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis ConfCentBins{"ConfCentBins", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 999.f}, "Mixing bins - centrality"};
  using BinningType_M = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;
  using BinningType_A = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0A>;
  using BinningType_C = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
  BinningType_M colBinning_M{{ConfVtxBins, ConfCentBins}, true};
  BinningType_A colBinning_A{{ConfVtxBins, ConfCentBins}, true};
  BinningType_C colBinning_C{{ConfVtxBins, ConfCentBins}, true};

  template <PairType pairtype, typename TEvents, typename TPhotons1, typename TPhotons2, typename TPreslice1, typename TPreslice2, typename TCuts1, typename TCuts2, typename TPairCuts, typename TLegs, typename TEMPrimaryElectrons, typename TMixedBinning>
  void MixedEventPairing(TEvents const& collisions, TPhotons1 const& photons1, TPhotons2 const& photons2, TPreslice1 const& perCollision1, TPreslice2 const& perCollision2, TCuts1 const& cuts1, TCuts2 const& cuts2, TPairCuts const& paircuts, TLegs const& /*legs*/, TEMPrimaryElectrons const& /*emprimaryelectrons*/, TMixedBinning const& colBinning)
  {
    THashList* list_pair_ss = static_cast<THashList*>(fMainList->FindObject("Pair")->FindObject(pairnames[pairtype].data()));

    // LOGF(info, "Number of collisions after filtering: %d", collisions.size());
    for (auto& [collision1, collision2] : soa::selfCombinations(colBinning, ndepth, -1, collisions, collisions)) { // internally, CombinationsStrictlyUpperIndexPolicy(collisions, collisions) is called.

      // LOGF(info, "Mixed event globalIndex: (%d, %d) , counter = %d, ngpcm: (%d, %d), ngphos: (%d, %d), ngemc: (%d, %d)",
      //     collision1.globalIndex(), collision2.globalIndex(), nev, collision1.ngpcm(), collision2.ngpcm(), collision1.ngphos(), collision2.ngphos(), collision1.ngemc(), collision2.ngemc());

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

      auto photons_coll1 = photons1.sliceBy(perCollision1, collision1.globalIndex());
      auto photons_coll2 = photons2.sliceBy(perCollision2, collision2.globalIndex());
      // LOGF(info, "collision1: posZ = %f, numContrib = %d , sel8 = %d | collision2: posZ = %f, numContrib = %d , sel8 = %d",
      //     collision1.posZ(), collision1.numContrib(), collision1.sel8(), collision2.posZ(), collision2.numContrib(), collision2.sel8());

      for (auto& cut1 : cuts1) {
        for (auto& cut2 : cuts2) {
          for (auto& paircut : paircuts) {
            for (auto& [g1, g2] : combinations(soa::CombinationsFullIndexPolicy(photons_coll1, photons_coll2))) {
              // LOGF(info, "Mixed event photon pair: (%d, %d) from events (%d, %d), photon event: (%d, %d)", g1.index(), g2.index(), collision1.index(), collision2.index(), g1.globalIndex(), g2.globalIndex());

              if ((pairtype == PairType::kPCMPCM || pairtype == PairType::kPHOSPHOS || pairtype == PairType::kEMCEMC) && (TString(cut1.GetName()) != TString(cut2.GetName()))) {
                continue;
              }
              if (!IsSelectedPair<pairtype>(g1, g2, cut1, cut2)) {
                continue;
              }
              if (!paircut.IsSelected(g1, g2)) {
                continue;
              }

              ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.); // pcm
              ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.); // phos or emc or dalitzee
              if constexpr (pairtype == PairType::kPCMDalitzEE) {
                v2.SetM(g2.mass());
                auto pos_sv = g1.template posTrack_as<aod::V0Legs>();
                auto ele_sv = g1.template negTrack_as<aod::V0Legs>();
                auto pos_pv = g2.template posTrack_as<MyPrimaryElectrons>();
                auto ele_pv = g2.template negTrack_as<MyPrimaryElectrons>();
                if (pos_sv.trackId() == pos_pv.trackId() || ele_sv.trackId() == ele_pv.trackId()) {
                  continue;
                }
              }
              ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
              // if (abs(v12.Rapidity()) > maxY) {
              //   continue;
              // }
              reinterpret_cast<TH2F*>(list_pair_ss->FindObject(Form("%s_%s", cut1.GetName(), cut2.GetName()))->FindObject(paircut.GetName())->FindObject("hMggPt_Mixed"))->Fill(v12.M(), v1.Pt());

            } // end of different photon combinations
          }   // end of pair cut loop
        }     // end of cut2 loop
      }       // end of cut1 loop
    }         // end of different collision combinations
  }

  Partition<MyCollisions> grouped_collisions = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax); // this goes to same event.
  Filter collisionFilter_common = nabs(o2::aod::collision::posZ) < 10.f && o2::aod::collision::numContrib > (uint16_t)0 && o2::aod::evsel::sel8 == true;
  Filter collisionFilter_centrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  Filter collisionFilter_subsys = (o2::aod::emevent::ngpcm >= 1) || (o2::aod::emevent::ngphos >= 1) || (o2::aod::emevent::ngemc >= 1) || (o2::aod::emevent::neeuls >= 1);
  using MyFilteredCollisions = soa::Filtered<MyCollisions>;

  Filter DalitzEEFilter = o2::aod::dalitzee::sign == 0; // analyze only uls
  using MyFilteredDalitzEEs = soa::Filtered<MyDalitzEEs>;

  Preslice<MyV0Photons> perCollision_pcm = aod::v0photonkf::emeventId;
  Preslice<MyDalitzEEs> perCollision_dalitz = aod::dalitzee::emeventId;
  Preslice<aod::PHOSClusters> perCollision_phos = aod::skimmedcluster::collisionId;
  Preslice<aod::SkimEMCClusters> perCollision_emc = aod::skimmedcluster::collisionId;

  void processPCMDalitzEE(MyCollisions const&, MyFilteredCollisions const& filtered_collisions, MyV0Photons const& v0photons, aod::V0Legs const& legs, MyFilteredDalitzEEs const& dielectrons, MyPrimaryElectrons const& emprimaryelectrons)
  {
    SameEventPairing<PairType::kPCMDalitzEE>(grouped_collisions, v0photons, dielectrons, perCollision_pcm, perCollision_dalitz, fPCMCuts, fDalitzEECuts, fPairCuts, legs, emprimaryelectrons);
    if (cfgCentEstimator == 0) {
      MixedEventPairing<PairType::kPCMDalitzEE>(filtered_collisions, v0photons, dielectrons, perCollision_pcm, perCollision_dalitz, fPCMCuts, fDalitzEECuts, fPairCuts, legs, emprimaryelectrons, colBinning_M);
    } else if (cfgCentEstimator == 1) {
      MixedEventPairing<PairType::kPCMDalitzEE>(filtered_collisions, v0photons, dielectrons, perCollision_pcm, perCollision_dalitz, fPCMCuts, fDalitzEECuts, fPairCuts, legs, emprimaryelectrons, colBinning_A);
    } else if (cfgCentEstimator == 2) {
      MixedEventPairing<PairType::kPCMDalitzEE>(filtered_collisions, v0photons, dielectrons, perCollision_pcm, perCollision_dalitz, fPCMCuts, fDalitzEECuts, fPairCuts, legs, emprimaryelectrons, colBinning_C);
    }
  }

  void processPCMPHOS(MyCollisions const&, MyFilteredCollisions const& filtered_collisions, MyV0Photons const& v0photons, aod::PHOSClusters const& phosclusters, aod::V0Legs const& legs)
  {
    SameEventPairing<PairType::kPCMPHOS>(grouped_collisions, v0photons, phosclusters, perCollision_pcm, perCollision_phos, fPCMCuts, fPHOSCuts, fPairCuts, legs, nullptr);
    MixedEventPairing<PairType::kPCMPHOS>(filtered_collisions, v0photons, phosclusters, perCollision_pcm, perCollision_phos, fPCMCuts, fPHOSCuts, fPairCuts, legs, nullptr, colBinning_C);
  }

  void processPCMEMC(MyCollisions const&, MyFilteredCollisions const& filtered_collisions, MyV0Photons const& v0photons, aod::SkimEMCClusters const& emcclusters, aod::V0Legs const& legs)
  {
    SameEventPairing<PairType::kPCMEMC>(grouped_collisions, v0photons, emcclusters, perCollision_pcm, perCollision_emc, fPCMCuts, fEMCCuts, fPairCuts, legs, nullptr);
    MixedEventPairing<PairType::kPCMEMC>(filtered_collisions, v0photons, emcclusters, perCollision_pcm, perCollision_emc, fPCMCuts, fEMCCuts, fPairCuts, legs, nullptr, colBinning_C);
  }

  void processDummy(MyCollisions const&) {}

  PROCESS_SWITCH(TaggingPi0, processPCMDalitzEE, "pairing PCM-Dalitz", false);
  PROCESS_SWITCH(TaggingPi0, processPCMPHOS, "pairing PCM-PHOS", false);
  PROCESS_SWITCH(TaggingPi0, processPCMEMC, "pairing PCM-EMCal", false);
  PROCESS_SWITCH(TaggingPi0, processDummy, "Dummy function", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TaggingPi0>(cfgc, TaskName{"tagging-pi0"})};
}
