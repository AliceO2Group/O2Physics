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

#include "Common/Core/RecoDecay.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/PairUtilities.h"
#include "PWGEM/PhotonMeson/Utils/MCUtilities.h"
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
using namespace o2::aod::pwgem::mcutil;
using namespace o2::aod::pwgem::photon;

using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent, aod::EMEventsQvec, aod::EMMCEventLabels>;
using MyCollision = MyCollisions::iterator;

using MyV0Photons = soa::Join<aod::V0PhotonsKF, aod::V0KFEMEventIds>;
using MyV0Photon = MyV0Photons::iterator;

using MyDalitzEEs = soa::Join<aod::DalitzEEs, aod::DalitzEEEMEventIds>;
using MyDalitzEE = MyDalitzEEs::iterator;

using MyDalitzMuMus = soa::Join<aod::DalitzMuMus, aod::DalitzMuMuEMEventIds>;
using MyDalitzMuMu = MyDalitzMuMus::iterator;

using MyEMCClusters = soa::Join<aod::SkimEMCClusters, aod::EMEMCClusterMCLabels, aod::EMCEMEventIds>;
using MyEMCCluster = MyEMCClusters::iterator;

struct Pi0EtaToGammaGammaMC {
  using MyMCV0Legs = soa::Join<aod::V0Legs, aod::V0LegMCLabels>;
  using MyMCElectrons = soa::Join<aod::EMPrimaryElectrons, aod::EMPrimaryElectronMCLabels, aod::EMPrimaryElectronsPrefilterBit>;
  using MyMCMuons = soa::Join<aod::EMPrimaryMuons, aod::EMPrimaryMuonMCLabels, aod::EMPrimaryMuonsPrefilterBit>;

  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999, "max. centrality"};

  Configurable<float> maxY_track{"maxY_track", 0.9, "maximum rapidity for generated particles"};                 // for PCM and dielectron
  Configurable<float> minPhi_track{"minPhi_track", 0, "minimum azimuthal angle for generated particles"};        // for PCM and dielectron
  Configurable<float> maxPhi_track{"maxPhi_track", 2 * M_PI, "maximum azimuthal angle for generated particles"}; // for PCM and dielectron
  Configurable<float> maxRgen{"maxRgen", 90.f, "maximum radius for generated particles"};
  Configurable<float> margin_z_mc{"margin_z_mc", 7.0, "margin for z cut in cm for MC"};

  Configurable<float> maxY_phos{"maxY_phos", 0.9, "maximum rapidity for generated particles"};                     // for EMC
  Configurable<float> minPhi_phos{"minPhi_phos", 3 / 2 * M_PI, "minimum azimuthal angle for generated particles"}; // for PCM and dielectron
  Configurable<float> maxPhi_phos{"maxPhi_phos", 2 * M_PI, "maximum azimuthal angle for generated particles"};     // for PCM and dielectron

  Configurable<float> maxY_emc{"maxY_emc", 0.9, "maximum rapidity for generated particles"};                  // for EMC
  Configurable<float> minPhi_emc{"minPhi_emc", M_PI / 2., "minimum azimuthal angle for generated particles"}; // for PCM and dielectron
  Configurable<float> maxPhi_emc{"maxPhi_emc", M_PI, "maximum azimuthal angle for generated particles"};      // for PCM and dielectron

  Configurable<std::string> fConfigPCMCuts{"cfgPCMCuts", "qc,nocut", "Comma separated list of V0 photon cuts"};
  Configurable<std::string> fConfigDalitzEECuts{"cfgDalitzEECuts", "mee_0_120_tpchadrejortofreq_lowB,mee_120_500_tpchadrejortofreq_lowB,mee_0_500_tpchadrejortofreq_lowB", "Comma separated list of Dalitz ee cuts"};
  Configurable<std::string> fConfigDalitzMuMuCuts{"cfgDalitzMuMuCuts", "mmumu_0_500_tpctof_lowB", "Comma separated list of Dalitz mumu cuts"};
  Configurable<std::string> fConfigEMCCuts{"cfgEMCCuts", "custom,standard,nocut", "Comma separated list of EMCal photon cuts"};

  // Configurable for EMCal cuts
  Configurable<bool> requireCaloReadout{"requireCaloReadout", true, "Require calorimeters readout when analyzing EMCal/PHOS"};
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

  Configurable<std::string> fConfigPairCuts{"cfgPairCuts", "nocut,asym08", "Comma separated list of pair cuts"};

  Configurable<std::string> fConfigEMEventCut{"cfgEMEventCut", "minbias", "em event cut"}; // only 1 event cut per wagon
  EMEventCut fEMEventCut;
  static constexpr std::string_view event_types[2] = {"before", "after"};

  OutputObj<THashList> fOutputEvent{"Event"};
  OutputObj<THashList> fOutputPair{"Pair"}; // 2-photon pair
  OutputObj<THashList> fOutputGen{"Generated"};
  THashList* fMainList = new THashList();

  std::vector<V0PhotonCut> fPCMCuts;
  std::vector<DalitzEECut> fDalitzEECuts;
  std::vector<DalitzEECut> fDalitzMuMuCuts;
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
    if (context.mOptions.get<bool>("processPCMDalitzEE")) {
      fPairNames.push_back("PCMDalitzEE");
    }
    if (context.mOptions.get<bool>("processPCMDalitzMuMu")) {
      fPairNames.push_back("PCMDalitzMuMu");
    }
    if (context.mOptions.get<bool>("processPHOSEMC")) {
      fPairNames.push_back("PHOSEMC");
    }

    DefinePCMCuts();
    DefineDalitzEECuts();
    DefineDalitzMuMuCuts();
    DefineEMCCuts();
    DefinePairCuts();
    addhistograms();
    TString ev_cut_name = fConfigEMEventCut.value;
    fEMEventCut = *eventcuts::GetCut(ev_cut_name.Data());

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
        o2::aod::pwgem::photon::histogram::AddHistClass(list_pair_subsys, photon_cut_name.data());
        THashList* list_pair_subsys_photoncut = reinterpret_cast<THashList*>(list_pair_subsys->FindObject(photon_cut_name.data()));

        for (auto& cut3 : cuts3) {
          std::string pair_cut_name = cut3.GetName();
          o2::aod::pwgem::photon::histogram::AddHistClass(list_pair_subsys_photoncut, pair_cut_name.data());
          THashList* list_pair_subsys_paircut = reinterpret_cast<THashList*>(list_pair_subsys_photoncut->FindObject(pair_cut_name.data()));
          o2::aod::pwgem::photon::histogram::DefineHistograms(list_pair_subsys_paircut, "gammagamma_mass_pt_mc");
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

    o2::aod::pwgem::photon::histogram::AddHistClass(fMainList, "Generated");
    THashList* list_gen = reinterpret_cast<THashList*>(fMainList->FindObject("Generated"));

    for (auto& pairname : fPairNames) {
      LOGF(info, "Enabled pairs = %s", pairname.data());

      // for events
      THashList* list_ev_pair = reinterpret_cast<THashList*>(o2::aod::pwgem::photon::histogram::AddHistClass(list_ev, pairname.data()));
      for (const auto& evtype : event_types) {
        THashList* list_ev_type = reinterpret_cast<THashList*>(o2::aod::pwgem::photon::histogram::AddHistClass(list_ev_pair, evtype.data()));
        o2::aod::pwgem::photon::histogram::DefineHistograms(list_ev_type, "Event", evtype.data());
      }

      // for generated particles
      THashList* list_gen_pair = reinterpret_cast<THashList*>(o2::aod::pwgem::photon::histogram::AddHistClass(list_gen, pairname.data()));
      o2::aod::pwgem::photon::histogram::DefineHistograms(list_gen_pair, "Generated", "Pi0Eta");

      // for truely reconstructed particles
      o2::aod::pwgem::photon::histogram::AddHistClass(list_pair, pairname.data());
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
      if (pairname == "PCMDalitzEE") {
        add_pair_histograms(list_pair, pairname, fPCMCuts, fDalitzEECuts, fPairCuts);
      }
      if (pairname == "PCMDalitzMuMu") {
        add_pair_histograms(list_pair, pairname, fPCMCuts, fDalitzMuMuCuts, fPairCuts);
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

  void DefineDalitzMuMuCuts()
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
    LOGF(info, "Number of DalitzMuMu cuts = %d", fDalitzMuMuCuts.size());
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
    } else if constexpr (pairtype == PairType::kPCMDalitzEE) {
      is_selected_pair = o2::aod::photonpair::IsSelectedPair<MyMCV0Legs, MyMCElectrons>(g1, g2, cut1, cut2);
    } else if constexpr (pairtype == PairType::kPCMDalitzMuMu) {
      is_selected_pair = o2::aod::photonpair::IsSelectedPair<MyMCV0Legs, MyMCMuons>(g1, g2, cut1, cut2);
    } else if constexpr (pairtype == PairType::kPHOSEMC) {
      is_selected_pair = o2::aod::photonpair::IsSelectedPair<int, aod::SkimEMCMTs>(g1, g2, cut1, cut2);
    } else {
      is_selected_pair = true;
    }
    return is_selected_pair;
  }

  template <PairType pairtype, typename TMCParticle, typename TMCParticles>
  bool IsInAcceptance(TMCParticle const& mcparticle, TMCParticles const& mcparticles)
  {
    if (pairtype == PairType::kPCMPCM) {
      return o2::aod::pwgem::mcutil::IsInAcceptance(mcparticle, mcparticles, std::vector<int>{22, 22}, -maxY_track, +maxY_track, minPhi_track, maxPhi_track);
    } else if (pairtype == PairType::kPCMDalitzEE) {
      return o2::aod::pwgem::mcutil::IsInAcceptance(mcparticle, mcparticles, std::vector<int>{-11, 11, 22}, -maxY_track, +maxY_track, minPhi_track, maxPhi_track);
    } else if (pairtype == PairType::kPCMDalitzMuMu) {
      return o2::aod::pwgem::mcutil::IsInAcceptance(mcparticle, mcparticles, std::vector<int>{-13, 13, 22}, -maxY_track, +maxY_track, minPhi_track, maxPhi_track);
    } else if (pairtype == PairType::kPHOSPHOS) {
      return o2::aod::pwgem::mcutil::IsInAcceptance(mcparticle, mcparticles, std::vector<int>{-13, 13, 22}, -maxY_phos, +maxY_phos, minPhi_phos, maxPhi_phos);
    } else if (pairtype == PairType::kEMCEMC) {
      return o2::aod::pwgem::mcutil::IsInAcceptance(mcparticle, mcparticles, std::vector<int>{22, 22}, -maxY_emc, +maxY_emc, minPhi_emc, maxPhi_emc);
    }
    return true;
  }

  Preslice<MyV0Photons> perCollision_pcm = aod::v0photonkf::emeventId;
  Preslice<MyDalitzEEs> perCollision_dalitzee = aod::dalitzee::emeventId;
  Preslice<MyDalitzMuMus> perCollision_dalitzmumu = aod::dalitzmumu::emeventId;
  Preslice<MyEMCClusters> perCollision_emc = aod::emccluster::emeventId;

  template <PairType pairtype, typename TEvents, typename TPhotons1, typename TPhotons2, typename TPreslice1, typename TPreslice2, typename TCuts1, typename TCuts2, typename TPairCuts, typename TV0Legs, typename TEMPrimaryElectrons, typename TEMPrimaryMuons, typename TMCEvents, typename TMCParticles>
  void TruePairing(TEvents const& collisions, TPhotons1 const& photons1, TPhotons2 const& photons2, TPreslice1 const& perCollision1, TPreslice2 const& perCollision2, TCuts1 const& cuts1, TCuts2 const& cuts2, TPairCuts const& paircuts, TV0Legs const& /*v0legs*/, TEMPrimaryElectrons const& /*emprimaryelectrons*/, TEMPrimaryMuons const& /*emprimarymuons*/, TMCEvents const& /*mcevents*/, TMCParticles const& mcparticles)
  {
    THashList* list_ev_pair_before = static_cast<THashList*>(fMainList->FindObject("Event")->FindObject(pairnames[pairtype].data())->FindObject(event_types[0].data()));
    THashList* list_ev_pair_after = static_cast<THashList*>(fMainList->FindObject("Event")->FindObject(pairnames[pairtype].data())->FindObject(event_types[1].data()));
    THashList* list_pair_ss = static_cast<THashList*>(fMainList->FindObject("Pair")->FindObject(pairnames[pairtype].data()));

    for (auto& collision : collisions) {

      if ((pairtype == kPHOSPHOS || pairtype == kPCMPHOS) && !collision.alias_bit(triggerAliases::kTVXinPHOS)) {
        continue;
      }
      if ((pairtype == kEMCEMC || pairtype == kPCMEMC) && ((!collision.alias_bit(triggerAliases::kTVXinEMC) && requireCaloReadout) || collision.ncollsPerBC() != 1)) {
        continue;
      }

      float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
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

      int pi0id = -1;
      int etaid = -1;
      if constexpr (pairtype == PairType::kPCMPCM || pairtype == PairType::kPHOSPHOS || pairtype == PairType::kEMCEMC) {
        for (auto& cut : cuts1) {
          for (auto& paircut : paircuts) {
            for (auto& [g1, g2] : combinations(CombinationsStrictlyUpperIndexPolicy(photons1_coll, photons2_coll))) {
              if (!IsSelectedPair<pairtype>(g1, g2, cut, cut)) {
                continue;
              }
              if (!paircut.IsSelected(g1, g2)) {
                continue;
              }

              int photonid1 = -1;
              int photonid2 = -1;
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

                photonid1 = FindCommonMotherFrom2Prongs(pos1mc, ele1mc, -11, 11, 22, mcparticles);
                photonid2 = FindCommonMotherFrom2Prongs(pos2mc, ele2mc, -11, 11, 22, mcparticles);
              } else if constexpr (pairtype == PairType::kEMCEMC) {
                auto cluster1mcparticle = mcparticles.iteratorAt(g1.emmcparticleId());
                auto cluster2mcparticle = mcparticles.iteratorAt(g2.emmcparticleId());

                photonid1 = FindMotherInChain(cluster1mcparticle, mcparticles, std::vector<int>{111, 221});
                photonid2 = FindMotherInChain(cluster2mcparticle, mcparticles, std::vector<int>{111, 221});
              }

              if (photonid1 < 0 || photonid2 < 0) {
                continue;
              }
              auto g1mc = mcparticles.iteratorAt(photonid1);
              auto g2mc = mcparticles.iteratorAt(photonid2);

              pi0id = FindCommonMotherFrom2Prongs(g1mc, g2mc, 22, 22, 111, mcparticles);
              etaid = FindCommonMotherFrom2Prongs(g1mc, g2mc, 22, 22, 221, mcparticles);

              if (pi0id < 0 && etaid < 0) {
                continue;
              }
              ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
              ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
              ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
              if (abs(v12.Rapidity()) > maxY_track) {
                continue;
              }

              if (pi0id > 0) {
                auto pi0mc = mcparticles.iteratorAt(pi0id);
                if (pi0mc.isPhysicalPrimary() || pi0mc.producedByGenerator()) {
                  if constexpr (pairtype == PairType::kPCMPCM) {
                    if (!IsConversionPointInAcceptance(g1mc, maxRgen, maxY_track, margin_z_mc, mcparticles) || !IsConversionPointInAcceptance(g2mc, maxRgen, maxY_track, margin_z_mc, mcparticles)) {
                      continue;
                    }
                  }
                  reinterpret_cast<TH2F*>(list_pair_ss->FindObject(Form("%s_%s", cut.GetName(), cut.GetName()))->FindObject(paircut.GetName())->FindObject("hMggPt_Pi0_Primary"))->Fill(v12.M(), v12.Pt());
                } else if (IsFromWD(pi0mc.emmcevent(), pi0mc, mcparticles)) {
                  reinterpret_cast<TH2F*>(list_pair_ss->FindObject(Form("%s_%s", cut.GetName(), cut.GetName()))->FindObject(paircut.GetName())->FindObject("hMggPt_Pi0_FromWD"))->Fill(v12.M(), v12.Pt());
                }
              } else if (etaid > 0) {
                auto etamc = mcparticles.iteratorAt(etaid);
                if (etamc.isPhysicalPrimary() || etamc.producedByGenerator()) {
                  if constexpr (pairtype == PairType::kPCMPCM) {
                    if (!IsConversionPointInAcceptance(g1mc, maxRgen, maxY_track, margin_z_mc, mcparticles) || !IsConversionPointInAcceptance(g2mc, maxRgen, maxY_track, margin_z_mc, mcparticles)) {
                      continue;
                    }
                  }
                  reinterpret_cast<TH2F*>(list_pair_ss->FindObject(Form("%s_%s", cut.GetName(), cut.GetName()))->FindObject(paircut.GetName())->FindObject("hMggPt_Eta_Primary"))->Fill(v12.M(), v12.Pt());
                }
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

                int photonid1 = -1;
                int photonid2 = -1;

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
                } else if constexpr (pairtype == PairType::kPCMDalitzEE) { // check 4 legs
                  auto pos1 = g1.template posTrack_as<MyMCV0Legs>();
                  auto ele1 = g1.template negTrack_as<MyMCV0Legs>();
                  auto pos2 = g2.template posTrack_as<MyMCElectrons>();
                  auto ele2 = g2.template negTrack_as<MyMCElectrons>();
                  if (pos1.trackId() == pos2.trackId() || ele1.trackId() == ele2.trackId()) {
                    continue;
                  }

                  auto pos1mc = pos1.template emmcparticle_as<aod::EMMCParticles>();
                  auto ele1mc = ele1.template emmcparticle_as<aod::EMMCParticles>();
                  auto pos2mc = pos2.template emmcparticle_as<aod::EMMCParticles>();
                  auto ele2mc = ele2.template emmcparticle_as<aod::EMMCParticles>();
                  // LOGF(info,"pos1mc.globalIndex() = %d , ele1mc.globalIndex() = %d , pos2mc.globalIndex() = %d , ele2mc.globalIndex() = %d", pos1mc.globalIndex(), ele1mc.globalIndex(), pos2mc.globalIndex(), ele2mc.globalIndex());

                  photonid1 = FindCommonMotherFrom2Prongs(pos1mc, ele1mc, -11, 11, 22, mcparticles); // real photon
                  if (photonid1 < 0) {
                    continue;
                  }
                  auto g1mc = mcparticles.iteratorAt(photonid1);

                  if (cut2.IsPhotonConversionSelected()) {                                             // v0photon + photon conversion on ITSib stored in dielectron table. pi0 -> gamma gamma
                    photonid2 = FindCommonMotherFrom2Prongs(pos2mc, ele2mc, -11, 11, 22, mcparticles); // photon conversion stored in dielectron table
                    if (photonid2 < 0) {
                      continue;
                    }
                    auto g2mc = mcparticles.iteratorAt(photonid2);
                    pi0id = FindCommonMotherFrom2Prongs(g1mc, g2mc, 22, 22, 111, mcparticles);
                    etaid = FindCommonMotherFrom2Prongs(g1mc, g2mc, 22, 22, 221, mcparticles);
                  } else { // pi0/eta -> ee gamma, dalitz decay
                    pi0id = FindCommonMotherFrom3Prongs(g1mc, pos2mc, ele2mc, 22, -11, 11, 111, mcparticles);
                    etaid = FindCommonMotherFrom3Prongs(g1mc, pos2mc, ele2mc, 22, -11, 11, 221, mcparticles);
                  }
                } else if constexpr (pairtype == PairType::kPCMDalitzMuMu) { // check 4 legs
                  auto pos1 = g1.template posTrack_as<MyMCV0Legs>();
                  auto ele1 = g1.template negTrack_as<MyMCV0Legs>();
                  auto pos2 = g2.template posTrack_as<MyMCMuons>();
                  auto ele2 = g2.template negTrack_as<MyMCMuons>();
                  if (pos1.trackId() == pos2.trackId() || ele1.trackId() == ele2.trackId()) {
                    continue;
                  }

                  auto pos1mc = pos1.template emmcparticle_as<aod::EMMCParticles>();
                  auto ele1mc = ele1.template emmcparticle_as<aod::EMMCParticles>();
                  auto pos2mc = pos2.template emmcparticle_as<aod::EMMCParticles>();
                  auto ele2mc = ele2.template emmcparticle_as<aod::EMMCParticles>();
                  // LOGF(info,"pos1mc.globalIndex() = %d , ele1mc.globalIndex() = %d , pos2mc.globalIndex() = %d , ele2mc.globalIndex() = %d", pos1mc.globalIndex(), ele1mc.globalIndex(), pos2mc.globalIndex(), ele2mc.globalIndex());

                  photonid1 = FindCommonMotherFrom2Prongs(pos1mc, ele1mc, -11, 11, 22, mcparticles); // real photon
                  if (photonid1 < 0) {
                    continue;
                  }
                  auto g1mc = mcparticles.iteratorAt(photonid1);
                  pi0id = -1;
                  etaid = FindCommonMotherFrom3Prongs(g1mc, pos2mc, ele2mc, 22, -13, 13, 221, mcparticles);
                }

                if (pi0id < 0 && etaid < 0) {
                  continue;
                }

                ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
                ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
                if constexpr (pairtype == PairType::kPCMDalitzEE || pairtype == PairType::kPCMDalitzMuMu) {
                  v2.SetM(g2.mass());
                }
                ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
                if (abs(v12.Rapidity()) > maxY_track) {
                  continue;
                }
                if (pi0id > 0) {
                  auto pi0mc = mcparticles.iteratorAt(pi0id);
                  if (pi0mc.isPhysicalPrimary() || pi0mc.producedByGenerator()) {
                    if constexpr (pairtype == PairType::kPCMDalitzEE || pairtype == kPCMDalitzMuMu) {
                      auto g1mc = mcparticles.iteratorAt(photonid1);
                      if (!IsConversionPointInAcceptance(g1mc, maxRgen, maxY_track, margin_z_mc, mcparticles)) {
                        continue;
                      }
                    }
                    reinterpret_cast<TH2F*>(list_pair_ss->FindObject(Form("%s_%s", cut1.GetName(), cut2.GetName()))->FindObject(paircut.GetName())->FindObject("hMggPt_Pi0_Primary"))->Fill(v12.M(), v12.Pt());
                  } else if (IsFromWD(pi0mc.emmcevent(), pi0mc, mcparticles)) {
                    reinterpret_cast<TH2F*>(list_pair_ss->FindObject(Form("%s_%s", cut1.GetName(), cut2.GetName()))->FindObject(paircut.GetName())->FindObject("hMggPt_Pi0_FromWD"))->Fill(v12.M(), v12.Pt());
                  }
                } else if (etaid > 0) {
                  auto etamc = mcparticles.iteratorAt(etaid);
                  if (etamc.isPhysicalPrimary() || etamc.producedByGenerator()) {
                    if constexpr (pairtype == PairType::kPCMDalitzEE || pairtype == kPCMDalitzMuMu) {
                      auto g1mc = mcparticles.iteratorAt(photonid1);
                      if (!IsConversionPointInAcceptance(g1mc, maxRgen, maxY_track, margin_z_mc, mcparticles)) {
                        continue;
                      }
                    }
                    reinterpret_cast<TH2F*>(list_pair_ss->FindObject(Form("%s_%s", cut1.GetName(), cut2.GetName()))->FindObject(paircut.GetName())->FindObject("hMggPt_Eta_Primary"))->Fill(v12.M(), v12.Pt());
                  }
                }

              } // end of combination
            }   // end of paircutloop
          }     // end of cut2 loop
        }       // end of cut1 loop
      }

    } // end of collision loop
  }

  PresliceUnsorted<aod::EMMCParticles> perMcCollision = aod::emmcparticle::emmceventId;
  PresliceUnsorted<MyCollisions> rec_perMcCollision = aod::emmceventlabel::emmceventId;

  template <PairType pairtype, typename TCollisions, typename TMCCollisions, typename TMCParticles>
  void runGenInfo(TCollisions const& collisions, TMCCollisions const& mccollisions, TMCParticles const& mcparticles)
  {
    THashList* list_gen_pair = static_cast<THashList*>(fMainList->FindObject("Generated")->FindObject(pairnames[pairtype].data()));
    // loop over mc stack and fill histograms for pure MC truth signals
    // all MC tracks which belong to the MC event corresponding to the current reconstructed event
    for (auto& mccollision : mccollisions) {
      auto collision_per_mccoll = collisions.sliceBy(rec_perMcCollision, mccollision.globalIndex());
      int nrec_per_mc = collision_per_mccoll.size();
      reinterpret_cast<TH1F*>(list_gen_pair->FindObject("hNrecPerMCCollision"))->Fill(nrec_per_mc); // all
    }

    for (auto& collision : collisions) {
      if ((pairtype == kPHOSPHOS || pairtype == kPCMPHOS) && !collision.alias_bit(triggerAliases::kTVXinPHOS)) {
        continue; // I don't know why this is necessary in simulation.
      }
      if ((pairtype == kEMCEMC || pairtype == kPCMEMC) && ((!collision.alias_bit(triggerAliases::kTVXinEMC) && requireCaloReadout) || collision.ncollsPerBC() != 1)) {
        continue; // I don't know why this is necessary in simulation.
      }

      float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }

      auto mccollision = collision.emmcevent();
      auto mctracks_coll = mcparticles.sliceBy(perMcCollision, mccollision.globalIndex());
      for (auto& mctrack : mctracks_coll) {
        if (abs(mctrack.y()) > maxY_track) {
          continue;
        }
        int pdg = mctrack.pdgCode();

        if (abs(pdg) == 111 && (mctrack.isPhysicalPrimary() || mctrack.producedByGenerator())) {
          reinterpret_cast<TH1F*>(list_gen_pair->FindObject("hPt_Pi0"))->Fill(mctrack.pt());
          reinterpret_cast<TH1F*>(list_gen_pair->FindObject("hY_Pi0"))->Fill(mctrack.y());
          reinterpret_cast<TH1F*>(list_gen_pair->FindObject("hPhi_Pi0"))->Fill(mctrack.phi());
          if (IsInAcceptance<pairtype>(mctrack, mcparticles)) {
            reinterpret_cast<TH1F*>(list_gen_pair->FindObject("hPt_Pi0_Acc"))->Fill(mctrack.pt());
            reinterpret_cast<TH1F*>(list_gen_pair->FindObject("hY_Pi0_Acc"))->Fill(mctrack.y());
            reinterpret_cast<TH1F*>(list_gen_pair->FindObject("hPhi_Pi0_Acc"))->Fill(mctrack.phi());
          }
        } else if (abs(pdg) == 221 && (mctrack.isPhysicalPrimary() || mctrack.producedByGenerator())) {
          reinterpret_cast<TH1F*>(list_gen_pair->FindObject("hPt_Eta"))->Fill(mctrack.pt());
          reinterpret_cast<TH1F*>(list_gen_pair->FindObject("hY_Eta"))->Fill(mctrack.y());
          reinterpret_cast<TH1F*>(list_gen_pair->FindObject("hPhi_Eta"))->Fill(mctrack.phi());
          if (IsInAcceptance<pairtype>(mctrack, mcparticles)) {
            reinterpret_cast<TH1F*>(list_gen_pair->FindObject("hPt_Eta_Acc"))->Fill(mctrack.pt());
            reinterpret_cast<TH1F*>(list_gen_pair->FindObject("hY_Eta_Acc"))->Fill(mctrack.y());
            reinterpret_cast<TH1F*>(list_gen_pair->FindObject("hPhi_Eta_Acc"))->Fill(mctrack.phi());
          }
        }
      } // end of mc track loop
    }   // end of collision loop
  }

  Partition<MyCollisions> grouped_collisions = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax); // this goes to same event.

  void processPCMPCM(MyCollisions const&, MyV0Photons const& v0photons, MyMCV0Legs const& v0legs, aod::EMMCEvents const& mccollisions, aod::EMMCParticles const& mcparticles)
  {
    TruePairing<PairType::kPCMPCM>(grouped_collisions, v0photons, v0photons, perCollision_pcm, perCollision_pcm, fPCMCuts, fPCMCuts, fPairCuts, v0legs, nullptr, nullptr, mccollisions, mcparticles);
    runGenInfo<PairType::kPCMPCM>(grouped_collisions, mccollisions, mcparticles);
  }
  void processPHOSPHOS(MyCollisions const&) {}
  void processEMCEMC(MyCollisions const&, MyEMCClusters const& emcclusters, aod::EMMCEvents const& mccollisions, aod::EMMCParticles const& mcparticles)
  {
    TruePairing<PairType::kEMCEMC>(grouped_collisions, emcclusters, emcclusters, perCollision_emc, perCollision_emc, fEMCCuts, fEMCCuts, fPairCuts, nullptr, nullptr, nullptr, mccollisions, mcparticles);
    runGenInfo<PairType::kEMCEMC>(grouped_collisions, mccollisions, mcparticles);
  }
  void processPCMPHOS(MyCollisions const&) {}
  void processPCMEMC(MyCollisions const&) {}

  void processPCMDalitzEE(MyCollisions const&, MyV0Photons const& v0photons, MyMCV0Legs const& v0legs, MyDalitzEEs const& dileptons, MyMCElectrons const& emprimaryelectrons, aod::EMMCEvents const& mccollisions, aod::EMMCParticles const& mcparticles)
  {
    TruePairing<PairType::kPCMDalitzEE>(grouped_collisions, v0photons, dileptons, perCollision_pcm, perCollision_dalitzee, fPCMCuts, fDalitzEECuts, fPairCuts, v0legs, emprimaryelectrons, nullptr, mccollisions, mcparticles);
    runGenInfo<PairType::kPCMDalitzEE>(grouped_collisions, mccollisions, mcparticles);
  }

  void processPCMDalitzMuMu(MyCollisions const&, MyV0Photons const& v0photons, MyMCV0Legs const& v0legs, MyDalitzMuMus const& dileptons, MyMCMuons const& emprimarymuons, aod::EMMCEvents const& mccollisions, aod::EMMCParticles const& mcparticles)
  {
    TruePairing<PairType::kPCMDalitzMuMu>(grouped_collisions, v0photons, dileptons, perCollision_pcm, perCollision_dalitzmumu, fPCMCuts, fDalitzMuMuCuts, fPairCuts, v0legs, nullptr, emprimarymuons, mccollisions, mcparticles);
    runGenInfo<PairType::kPCMDalitzMuMu>(grouped_collisions, mccollisions, mcparticles);
  }

  void processPHOSEMC(MyCollisions const&) {}

  void processDummy(MyCollisions const&) {}

  PROCESS_SWITCH(Pi0EtaToGammaGammaMC, processPCMPCM, "true pairing PCM-PCM", false);
  PROCESS_SWITCH(Pi0EtaToGammaGammaMC, processPHOSPHOS, "true pairing PHOS-PHOS", false);
  PROCESS_SWITCH(Pi0EtaToGammaGammaMC, processEMCEMC, "true pairing EMC-EMC", false);
  PROCESS_SWITCH(Pi0EtaToGammaGammaMC, processPCMPHOS, "true pairing PCM-PHOS", false);
  PROCESS_SWITCH(Pi0EtaToGammaGammaMC, processPCMEMC, "true pairing PCM-EMC", false);
  PROCESS_SWITCH(Pi0EtaToGammaGammaMC, processPCMDalitzEE, "true pairing PCM-DalitzEE", false);
  PROCESS_SWITCH(Pi0EtaToGammaGammaMC, processPCMDalitzMuMu, "true pairing PCM-DalitzMuMu", false);
  PROCESS_SWITCH(Pi0EtaToGammaGammaMC, processPHOSEMC, "true pairing PHOS-EMC", false);
  PROCESS_SWITCH(Pi0EtaToGammaGammaMC, processDummy, "Dummy function", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<Pi0EtaToGammaGammaMC>(cfgc, TaskName{"pi0eta-to-gammagamma-mc"})};
}
