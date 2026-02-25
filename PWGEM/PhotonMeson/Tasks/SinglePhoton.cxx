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
// This code loops over photon candidate and fill histograms
//    Please write to: daiki.sekihata@cern.ch

#include "EMPhotonEventCut.h"

#include "PWGEM/PhotonMeson/Core/CutsLibrary.h"
#include "PWGEM/PhotonMeson/Core/EMCPhotonCut.h"
#include "PWGEM/PhotonMeson/Core/HistogramsLibrary.h"
#include "PWGEM/PhotonMeson/Core/PHOSPhotonCut.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"

#include "Common/CCDB/TriggerAliases.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TH2.h>
#include <THashList.h>
#include <TString.h>

#include <array>
#include <cmath>
#include <cstring>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::photon;

using MyCollisions = soa::Join<aod::EMEvents_004, aod::EMEventsAlias, aod::EMEventsMult_000, aod::EMEventsCent_000, aod::EMEventsQvec_001>;
using MyCollision = MyCollisions::iterator;

using MyV0Photons = soa::Join<aod::V0PhotonsKF, aod::V0KFEMEventIds>;
using MyV0Photon = MyV0Photons::iterator;

struct SinglePhoton {
  enum EMDetType {
    kPCM = 0,
    kPHOS = 1,
    kEMC = 2,
  };

  Configurable<bool> cfgDoFlow{"cfgDoFlow", false, "flag to analyze vn"};

  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", -1, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999, "max. centrality"};

  Configurable<float> maxY{"maxY", 0.9, "maximum rapidity for reconstructed particles"};
  Configurable<std::string> fConfigPCMCuts{"cfgPCMCuts", "qc,qc_ITSTPC,qc_ITSonly,qc_TPConly,wwire_ib", "Comma separated list of V0 photon cuts"};
  Configurable<std::string> fConfigPHOSCuts{"cfgPHOSCuts", "test02,test03", "Comma separated list of PHOS photon cuts"};
  Configurable<std::string> fConfigEMCCuts{"fConfigEMCCuts", "custom,standard,nocut", "Comma separated list of EMCal photon cuts"};

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
  EMPhotonEventCut fEMEventCut;
  static constexpr std::string_view event_types[2] = {"before", "after"};

  OutputObj<THashList> fOutputEvent{"Event"};
  OutputObj<THashList> fOutputPhoton{"Photon"}; // single photon
  THashList* fMainList = new THashList();

  std::vector<V0PhotonCut> fPCMCuts;
  std::vector<PHOSPhotonCut> fPHOSCuts;
  std::vector<EMCPhotonCut> fEMCCuts;

  std::vector<std::string> fDetNames;
  void init(InitContext& context)
  {
    if (context.mOptions.get<bool>("processPCM")) {
      fDetNames.push_back("PCM");
    }
    // if (context.mOptions.get<bool>("processPHOS")) {
    //   fDetNames.push_back("PHOS");
    // }
    // if (context.mOptions.get<bool>("processEMC")) {
    //   fDetNames.push_back("EMC");
    // }

    DefinePCMCuts();
    DefinePHOSCuts();
    DefineEMCCuts();
    addhistograms();
    TString ev_cut_name = fConfigEMEventCut.value;
    fEMEventCut = *eventcuts::GetCut(ev_cut_name.Data());

    fOutputEvent.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("Event")));
    fOutputPhoton.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("Photon")));
  }

  template <typename TCuts1>
  void add_histograms(THashList* list_photon, const std::string detname, TCuts1 const& cuts1)
  {
    for (auto& cut1 : cuts1) {
      std::string cutname1 = cut1.GetName();

      THashList* list_photon_subsys = reinterpret_cast<THashList*>(list_photon->FindObject(detname.data()));
      o2::aod::pwgem::photon::histogram::AddHistClass(list_photon_subsys, cutname1.data());
      THashList* list_photon_subsys_cut = reinterpret_cast<THashList*>(list_photon_subsys->FindObject(cutname1.data()));
      if (cfgDoFlow) {
        o2::aod::pwgem::photon::histogram::DefineHistograms(list_photon_subsys_cut, "singlephoton", "qvector");
      } else {
        o2::aod::pwgem::photon::histogram::DefineHistograms(list_photon_subsys_cut, "singlephoton", "");
      }

    } // end of cut1 loop
  }

  static constexpr std::string_view detnames[3] = {"PCM", "PHOS", "EMC"};
  void addhistograms()
  {
    fMainList->SetOwner(true);
    fMainList->SetName("fMainList");

    // create sub lists first.
    o2::aod::pwgem::photon::histogram::AddHistClass(fMainList, "Event");
    THashList* list_ev = reinterpret_cast<THashList*>(fMainList->FindObject("Event"));

    o2::aod::pwgem::photon::histogram::AddHistClass(fMainList, "Photon");
    THashList* list_photon = reinterpret_cast<THashList*>(fMainList->FindObject("Photon"));

    for (auto& detname : fDetNames) {
      LOGF(info, "Enabled detector = %s", detname.data());

      o2::aod::pwgem::photon::histogram::AddHistClass(list_ev, detname.data());
      THashList* list_ev_det = reinterpret_cast<THashList*>(list_ev->FindObject(detname.data()));
      for (const auto& evtype : event_types) {
        THashList* list_ev_det_type = reinterpret_cast<THashList*>(o2::aod::pwgem::photon::histogram::AddHistClass(list_ev_det, evtype.data()));
        if (cfgDoFlow) {
          o2::aod::pwgem::photon::histogram::DefineHistograms(list_ev_det_type, "Event", "qvector");
        } else {
          o2::aod::pwgem::photon::histogram::DefineHistograms(list_ev_det_type, "Event", "");
        }
      }

      o2::aod::pwgem::photon::histogram::AddHistClass(list_photon, detname.data());

      if (detname == "PCM") {
        add_histograms(list_photon, detname, fPCMCuts);
      }
      if (detname == "PHOS") {
        add_histograms(list_photon, detname, fPHOSCuts);
      }
      if (detname == "EMC") {
        add_histograms(list_photon, detname, fEMCCuts);
      }

    } // end of detector name loop
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

          custom_cut->SetTrackMatchingEtaParams(EMC_TM_Eta->at(0), EMC_TM_Eta->at(1), EMC_TM_Eta->at(2));
          custom_cut->SetTrackMatchingPhiParams(EMC_TM_Phi->at(0), EMC_TM_Phi->at(1), EMC_TM_Phi->at(2));

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

  Preslice<MyV0Photons> perCollision = aod::v0photonkf::emphotoneventId;
  // Preslice<aod::PHOSClusters> perCollision_phos = aod::skimmedcluster::collisionId;
  // Preslice<aod::SkimEMCClusters> perCollision_emc = aod::skimmedcluster::collisionId;

  template <EMDetType photontype, typename TG1, typename TCut1>
  bool IsSelected(TG1 const& g1, TCut1 const& cut1)
  {
    bool is_selected = false;
    if constexpr (photontype == EMDetType::kPCM) {
      is_selected = cut1.template IsSelected<decltype(g1), aod::V0Legs>(g1);
    } else if constexpr (photontype == EMDetType::kPHOS) {
      is_selected = cut1.template IsSelected<int>(g1); // dummy, because track matching is not ready.
      //} else if constexpr (photontype == EMDetType::kEMC) {
      //  is_selected = cut1.template IsSelected<aod::SkimEMCMTs>(g1);
    } else {
      is_selected = true;
    }
    return is_selected;
  }

  template <EMDetType photontype, typename TEvents, typename TPhotons1, typename TPreslice1, typename TCuts1, typename TV0Legs>
  void FillPhoton(TEvents const& collisions, TPhotons1 const& photons1, TPreslice1 const& perCollision1, TCuts1 const& cuts1, TV0Legs const&)
  {
    THashList* list_ev_before = static_cast<THashList*>(fMainList->FindObject("Event")->FindObject(detnames[photontype].data())->FindObject(event_types[0].data()));
    THashList* list_ev_after = static_cast<THashList*>(fMainList->FindObject("Event")->FindObject(detnames[photontype].data())->FindObject(event_types[1].data()));
    THashList* list_photon_det = static_cast<THashList*>(fMainList->FindObject("Photon")->FindObject(detnames[photontype].data()));

    for (auto& collision : collisions) {
      if (photontype == EMDetType::kPHOS && !collision.alias_bit(triggerAliases::kTVXinPHOS)) {
        continue;
      }
      if (photontype == EMDetType::kEMC && !collision.alias_bit(triggerAliases::kTVXinEMC)) {
        continue;
      }

      float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      if (cfgDoFlow) {
        o2::aod::pwgem::photon::histogram::FillHistClass<EMHistType::kEvent_Cent_Qvec>(list_ev_before, "", collision);
      } else {
        o2::aod::pwgem::photon::histogram::FillHistClass<EMHistType::kEvent>(list_ev_before, "", collision);
      }

      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }

      if (cfgDoFlow) {
        o2::aod::pwgem::photon::histogram::FillHistClass<EMHistType::kEvent_Cent_Qvec>(list_ev_after, "", collision);
      } else {
        o2::aod::pwgem::photon::histogram::FillHistClass<EMHistType::kEvent>(list_ev_after, "", collision);
      }
      reinterpret_cast<TH1F*>(list_ev_before->FindObject("hCollisionCounter"))->Fill("accepted", 1.f);
      reinterpret_cast<TH1F*>(list_ev_after->FindObject("hCollisionCounter"))->Fill("accepted", 1.f);
      std::array<float, 2> q2ft0m = {collision.q2xft0m(), collision.q2yft0m()};
      std::array<float, 2> q2ft0a = {collision.q2xft0a(), collision.q2yft0a()};
      std::array<float, 2> q2ft0c = {collision.q2xft0c(), collision.q2yft0c()};
      // std::array<float, 2> q2fv0a = {collision.q2xfv0a(), collision.q2yfv0a()};

      auto photons1_coll = photons1.sliceBy(perCollision1, collision.globalIndex());
      for (auto& cut : cuts1) {
        THashList* list_photon_det_cut = static_cast<THashList*>(list_photon_det->FindObject(cut.GetName()));
        for (auto& photon : photons1_coll) {
          if (!IsSelected<photontype>(photon, cut)) {
            continue;
          }
          if (abs(photon.eta()) > maxY) {
            continue;
          }
          if (cfgDoFlow) {
            std::array<float, 2> u2_photon = {std::cos(2 * photon.phi()), std::sin(2 * photon.phi())};
            reinterpret_cast<TH2F*>(list_photon_det_cut->FindObject("hPt_SPQ2FT0M"))->Fill(photon.pt(), RecoDecay::dotProd(u2_photon, q2ft0m));
            reinterpret_cast<TH2F*>(list_photon_det_cut->FindObject("hPt_SPQ2FT0A"))->Fill(photon.pt(), RecoDecay::dotProd(u2_photon, q2ft0a));
            reinterpret_cast<TH2F*>(list_photon_det_cut->FindObject("hPt_SPQ2FT0C"))->Fill(photon.pt(), RecoDecay::dotProd(u2_photon, q2ft0c));
            // reinterpret_cast<TH2F*>(list_photon_det_cut->FindObject("hPt_SPQ2FV0A"))->Fill(photon.pt(), RecoDecay::dotProd(u2_photon, q2fv0a));
          } else {
            reinterpret_cast<TH1F*>(list_photon_det_cut->FindObject("hPt"))->Fill(photon.pt());
          }

          reinterpret_cast<TH1F*>(list_photon_det_cut->FindObject("hY"))->Fill(photon.eta());
          reinterpret_cast<TH1F*>(list_photon_det_cut->FindObject("hPhi"))->Fill(photon.phi());
        } // end of photon loop
      } // end of cut loop
    } // end of collision loop
  }

  Partition<MyCollisions> grouped_collisions = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax); // this goes to same event.

  void processPCM(MyCollisions const&, MyV0Photons const& v0photons, aod::V0Legs const& legs)
  {
    FillPhoton<EMDetType::kPCM>(grouped_collisions, v0photons, perCollision, fPCMCuts, legs);
  }

  // void processPHOS(MyCollisions const& collisions, aod::PHOSClusters const& phosclusters)
  // {
  //   FillPhoton<EMDetType::kPHOS>(grouped_collisions, phosclusters, perCollision_phos, fPHOSCuts, nullptr);
  // }

  //  void processEMC(MyCollisions const& collisions, aod::SkimEMCClusters const& emcclusters)
  //  {
  //    FillPhoton<EMDetType::kEMC>(grouped_collisions, emcclusters, perCollision_emc, fEMCCuts, nullptr);
  //  }

  void processDummy(MyCollisions::iterator const&) {}

  PROCESS_SWITCH(SinglePhoton, processPCM, "single photon with PCM", false);
  // PROCESS_SWITCH(SinglePhoton, processPHOS, "single photon with PHOS", false);
  //  PROCESS_SWITCH(SinglePhoton, processEMC, "single photon with EMC", false);
  PROCESS_SWITCH(SinglePhoton, processDummy, "Dummy function", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SinglePhoton>(cfgc, TaskName{"single-photon"})};
}
