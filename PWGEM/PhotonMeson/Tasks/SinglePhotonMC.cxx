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

/// \file SinglePhotonMC.cxx
/// \brief HeadThis code loops over photon candidate and fill histograms
/// \author D. Sekihata, daiki.sekihata@cern.ch

#include "EMPhotonEventCut.h"

#include "PWGEM/Dilepton/Utils/MCUtilities.h"
#include "PWGEM/PhotonMeson/Core/CutsLibrary.h"
#include "PWGEM/PhotonMeson/Core/HistogramsLibrary.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/DataModel/EventTables.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/MCUtilities.h"

#include "Common/CCDB/TriggerAliases.h"
#include "Common/DataModel/Centrality.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <THashList.h>
#include <TString.h>

#include <ranges>
#include <string>
#include <string_view>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::photonmeson::utils::mcutil;
using namespace o2::aod::pwgem::dilepton::utils::mcutil;
using namespace o2::aod::pwgem::photon;

using MyCollisions = soa::Join<aod::PMEvents, aod::EMEventsAlias, aod::EMEventsMult_000, aod::EMEventsCent_000, aod::EMEventsQvec_001, aod::EMMCEventLabels>;
using MyCollision = MyCollisions::iterator;

using MyV0Photons = soa::Join<aod::V0PhotonsKF, aod::V0KFEMEventIds>;
using MyV0Photon = MyV0Photons::iterator;

using MyMCV0Legs = soa::Join<aod::V0Legs, aod::V0LegMCLabels>;

struct SinglePhotonMC {
  enum EMDetType {
    kPCM = 0,
    kPHOS = 1,
    kEMC = 2,
  };

  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", -1, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999, "max. centrality"};

  Configurable<float> maxY{"maxY", 0.9, "maximum rapidity for reconstructed particles"};
  Configurable<float> maxRgen{"maxRgen", 90.f, "maximum radius for generated particles"};
  Configurable<float> margin_z_mc{"margin_z_mc", 7.0, "margin for z cut in cm for MC"};

  Configurable<std::string> fConfigPCMCuts{"cfgPCMCuts", "analysis,qc,nocut", "Comma separated list of V0 photon cuts"};

  Configurable<std::string> fConfigEMEventCut{"cfgEMEventCut", "minbias", "em event cut"}; // only 1 event cut per wagon
  EMPhotonEventCut fEMEventCut;
  static constexpr std::array<std::string_view, 2> event_types = {"before", "after"};

  OutputObj<THashList> fOutputEvent{"Event"};
  OutputObj<THashList> fOutputPhoton{"Photon"}; // single photon
  OutputObj<THashList> fOutputGen{"Generated"};
  THashList* fMainList = new THashList();

  std::vector<V0PhotonCut> fPCMCuts;
  // std::vector<PHOSPhotonCut> fPHOSCuts;
  // std::vector<EMCPhotonCut> fEMCCuts;

  std::vector<std::string> fDetNames;
  void init(InitContext& context)
  {
    if (context.mOptions.get<bool>("processPCM")) {
      fDetNames.push_back("PCM");
    }

    DefinePCMCuts();
    addhistograms();
    TString ev_cut_name = fConfigEMEventCut.value;
    fEMEventCut = *eventcuts::GetCut(ev_cut_name.Data());

    fOutputEvent.setObject(dynamic_cast<THashList*>(fMainList->FindObject("Event")));
    fOutputPhoton.setObject(dynamic_cast<THashList*>(fMainList->FindObject("Photon")));
    fOutputGen.setObject(dynamic_cast<THashList*>(fMainList->FindObject("Generated")));
  }

  template <typename TCuts1>
  void add_photon_histograms(THashList* list_photon, std::string const& detname, TCuts1 const& cuts1)
  {
    for (auto& cut1 : cuts1) {
      std::string cutname1 = cut1.GetName();

      auto* list_photon_subsys = dynamic_cast<THashList*>(list_photon->FindObject(detname.data()));
      o2::aod::pwgem::photon::histogram::AddHistClass(list_photon_subsys, cutname1.data());
      auto* list_photon_subsys_cut = dynamic_cast<THashList*>(list_photon_subsys->FindObject(cutname1.data()));
      o2::aod::pwgem::photon::histogram::DefineHistograms(list_photon_subsys_cut, "singlephoton", "mc");
    } // end of cut1 loop
  }

  static constexpr std::array<std::string_view, 3> detnames = {"PCM", "PHOS", "EMC"};
  void addhistograms()
  {
    fMainList->SetOwner(true);
    fMainList->SetName("fMainList");

    // create sub lists first.
    o2::aod::pwgem::photon::histogram::AddHistClass(fMainList, "Event");
    auto* list_ev = dynamic_cast<THashList*>(fMainList->FindObject("Event"));

    o2::aod::pwgem::photon::histogram::AddHistClass(fMainList, "Photon");
    auto* list_photon = dynamic_cast<THashList*>(fMainList->FindObject("Photon"));

    o2::aod::pwgem::photon::histogram::AddHistClass(fMainList, "Generated");
    auto* list_gen = dynamic_cast<THashList*>(fMainList->FindObject("Generated"));
    o2::aod::pwgem::photon::histogram::DefineHistograms(list_gen, "Generated", "Photon");

    for (const auto& detname : fDetNames) {
      LOGF(info, "Enabled detector = %s", detname.data());

      o2::aod::pwgem::photon::histogram::AddHistClass(list_ev, detname.data());
      auto* list_ev_det = dynamic_cast<THashList*>(list_ev->FindObject(detname.data()));
      for (const auto& evtype : event_types) {
        THashList* list_ev_det_type = o2::aod::pwgem::photon::histogram::AddHistClass(list_ev_det, evtype.data());
        o2::aod::pwgem::photon::histogram::DefineHistograms(list_ev_det_type, "Event", evtype.data());
      }

      o2::aod::pwgem::photon::histogram::AddHistClass(list_photon, detname.data());

      if (detname == "PCM") {
        add_photon_histograms(list_photon, detname, fPCMCuts);
      }
    } // end of detector name loop
  }

  void DefinePCMCuts()
  {
    if (fConfigPCMCuts.value.empty()) {
      return;
    }

    std::string_view namesView(fConfigPCMCuts.value);

    for (auto name : namesView | std::views::split(',')) {
      std::string cutString(name.begin(), name.end());
      const char* cutname = cutString.c_str();
      LOGF(info, "add tag cut : %s", cutname);
      fPCMCuts.push_back(*pcmcuts::GetCut(cutname));
    }
    LOGF(info, "Number of PCM cuts = %d", fPCMCuts.size());
  }

  Preslice<MyV0Photons> perCollision = aod::v0photonkf::pmeventId;

  template <EMDetType photontype, typename TG1, typename TCut1>
  bool IsSelected(TG1 const& g1, TCut1 const& cut1)
  {
    bool is_selected = false;
    if constexpr (photontype == EMDetType::kPCM) {
      is_selected = cut1.template IsSelected<decltype(g1), MyMCV0Legs>(g1);
    } else if constexpr (photontype == EMDetType::kPHOS) {
      is_selected = cut1.template IsSelected<int>(g1); // dummy, because track matching is not ready.
      //} else if constexpr (photontype == EMDetType::kEMC) {
      //  is_selected = cut1.template IsSelected<aod::SkimEMCMTs>(g1);
    } else {
      is_selected = true;
    }
    return is_selected;
  }

  template <EMDetType photontype, typename TEvents, typename TPhotons1, typename TPreslice1, typename TCuts1, typename TV0Legs, typename TMCParticles, typename TMCEvents>
  void FillTruePhoton(TEvents const& collisions, TPhotons1 const& photons1, TPreslice1 const& perCollision1, TCuts1 const& cuts1, TV0Legs const&, TMCParticles const& mcparticles, TMCEvents const&)
  {
    auto* list_ev_before = dynamic_cast<THashList*>(fMainList->FindObject("Event")->FindObject(detnames[photontype].data())->FindObject(event_types[0].data()));
    auto* list_ev_after = dynamic_cast<THashList*>(fMainList->FindObject("Event")->FindObject(detnames[photontype].data())->FindObject(event_types[1].data()));
    auto* list_photon_det = dynamic_cast<THashList*>(fMainList->FindObject("Photon")->FindObject(detnames[photontype].data()));

    for (auto& collision : collisions) {
      if (photontype == EMDetType::kPHOS && !collision.alias_bit(triggerAliases::kTVXinPHOS)) {
        continue;
      }
      if (photontype == EMDetType::kEMC && !collision.alias_bit(triggerAliases::kTVXinEMC)) {
        continue;
      }

      std::array<float, 3> centralities{collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      o2::aod::pwgem::photon::histogram::FillHistClass<EMHistType::kEvent>(list_ev_before, "", collision);
      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      o2::aod::pwgem::photon::histogram::FillHistClass<EMHistType::kEvent>(list_ev_after, "", collision);
      dynamic_cast<TH1F*>(list_ev_before->FindObject("hCollisionCounter"))->Fill("accepted", 1.f);
      dynamic_cast<TH1F*>(list_ev_after->FindObject("hCollisionCounter"))->Fill("accepted", 1.f);

      auto photons1_coll = photons1.sliceBy(perCollision1, collision.globalIndex());
      for (auto& cut : cuts1) {
        auto* list_photon_det_cut = dynamic_cast<THashList*>(list_photon_det->FindObject(cut.getName().c_str()));
        for (auto& photon : photons1_coll) {
          if (!IsSelected<photontype>(photon, cut)) {
            continue;
          }
          if (abs(photon.eta()) > maxY) {
            continue;
          }
          dynamic_cast<TH1F*>(list_photon_det_cut->FindObject("hPt"))->Fill(photon.pt());
          dynamic_cast<TH1F*>(list_photon_det_cut->FindObject("hY"))->Fill(photon.eta());
          dynamic_cast<TH1F*>(list_photon_det_cut->FindObject("hPhi"))->Fill(photon.phi());

          int photonid = -1;
          if constexpr (photontype == EMDetType::kPCM) {
            auto pos = photon.template posTrack_as<MyMCV0Legs>();
            auto ele = photon.template negTrack_as<MyMCV0Legs>();
            auto posmc = pos.template emmcparticle_as<aod::EMMCParticles>();
            auto elemc = ele.template emmcparticle_as<aod::EMMCParticles>();
            photonid = FindCommonMotherFrom2Prongs(posmc, elemc, -11, 11, 22, mcparticles);
            if (photonid < 0) { // check swap, true electron is reconstructed as positron and vice versa.
              photonid = FindCommonMotherFrom2Prongs(posmc, elemc, 11, -11, 22, mcparticles);
            }
          }
          if (photonid < 0) {
            continue;
          }

          auto mcphoton = mcparticles.iteratorAt(photonid);

          if (mcphoton.isPhysicalPrimary() || mcphoton.producedByGenerator()) {
            if constexpr (photontype == EMDetType::kPCM) {
              if (!IsConversionPointInAcceptance(mcphoton, maxRgen, maxY, margin_z_mc, mcparticles)) {
                continue;
              }
            }
            dynamic_cast<TH1F*>(list_photon_det_cut->FindObject("hPt_Photon_Primary"))->Fill(photon.pt());
            dynamic_cast<TH1F*>(list_photon_det_cut->FindObject("hY_Photon_Primary"))->Fill(photon.eta());
            dynamic_cast<TH1F*>(list_photon_det_cut->FindObject("hPhi_Photon_Primary"))->Fill(photon.phi());
          } else if (IsFromWD(mcphoton.emmcevent(), mcphoton, mcparticles) > 0) {
            dynamic_cast<TH1F*>(list_photon_det_cut->FindObject("hPt_Photon_FromWD"))->Fill(photon.pt());
            dynamic_cast<TH1F*>(list_photon_det_cut->FindObject("hY_Photon_FromWD"))->Fill(photon.eta());
            dynamic_cast<TH1F*>(list_photon_det_cut->FindObject("hPhi_Photon_FromWD"))->Fill(photon.phi());
          } else {
            dynamic_cast<TH1F*>(list_photon_det_cut->FindObject("hPt_Photon_hs"))->Fill(photon.pt());
            dynamic_cast<TH1F*>(list_photon_det_cut->FindObject("hY_Photon_hs"))->Fill(photon.eta());
            dynamic_cast<TH1F*>(list_photon_det_cut->FindObject("hPhi_Photon_hs"))->Fill(photon.phi());
          }

        } // end of photon loop
      } // end of cut loop
    } // end of collision loop
  }

  Partition<MyCollisions> grouped_collisions = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax); // this goes to same event.

  void processPCM(MyCollisions const&, MyV0Photons const& v0photons, MyMCV0Legs const& legs, aod::EMMCParticles const& mcparticles, aod::EMMCEvents const& mccollisions)
  {
    FillTruePhoton<EMDetType::kPCM>(grouped_collisions, v0photons, perCollision, fPCMCuts, legs, mcparticles, mccollisions);
  }

  // void processPHOS(MyCollisions const& collisions, aod::PHOSClusters const& phosclusters, aod::EMMCParticles const& mcparticles, aod::EMMCEvents const& mccollisions)
  // {
  //   FillTruePhoton<EMDetType::kPHOS>(grouped_collisions, phosclusters, perCollision_phos, fPHOSCuts, nullptr, mcparticles, mccollisions);
  // }

  // void processEMC(MyCollisions const& collisions, aod::SkimEMCClusters const& emcclusters, aod::EMMCParticles const& mcparticles, aod::EMMCEvents const& mccollisions)
  // {
  //   FillTruePhoton<EMDetType::kEMC>(grouped_collisions, emcclusters, perCollision_emc, fEMCCuts, nullptr, mcparticles, mccollisions);
  // }

  PresliceUnsorted<aod::EMMCParticles> perMcCollision = aod::emmcparticle::emmceventId;
  void processGen(MyCollisions const&, aod::EMMCEvents const&, aod::EMMCParticles const& mcparticles)
  {
    // loop over mc stack and fill histograms for pure MC truth signals
    // all MC tracks which belong to the MC event corresponding to the current reconstructed event
    for (auto& collision : grouped_collisions) {
      std::array<float, 3> centralities{collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }
      auto mccollision = collision.emmcevent();

      dynamic_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hCollisionCounter"))->Fill(1.0);
      dynamic_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hZvtx_before"))->Fill(mccollision.posZ());
      if (!collision.sel8()) {
        continue;
      }
      dynamic_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hCollisionCounter"))->Fill(2.0);

      if (collision.numContrib() < 0.5) {
        continue;
      }
      dynamic_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hCollisionCounter"))->Fill(3.0);

      if (abs(collision.posZ()) > 10.0) {
        continue;
      }
      dynamic_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hCollisionCounter"))->Fill(4.0);
      dynamic_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hZvtx_after"))->Fill(mccollision.posZ());

      auto mctracks_coll = mcparticles.sliceBy(perMcCollision, collision.emmceventId());
      for (auto& mctrack : mctracks_coll) {
        if (abs(mctrack.y()) > maxY) {
          continue;
        }

        if (abs(mctrack.pdgCode()) == 22 && (mctrack.isPhysicalPrimary() || mctrack.producedByGenerator())) {
          dynamic_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hPt_Photon"))->Fill(mctrack.pt());
          dynamic_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hY_Photon"))->Fill(mctrack.y());
          dynamic_cast<TH1F*>(fMainList->FindObject("Generated")->FindObject("hPhi_Photon"))->Fill(mctrack.phi());
        }
      }
    }
  }
  void processDummy(MyCollisions const&) {}

  PROCESS_SWITCH(SinglePhotonMC, processPCM, "single photon with PCM", false);
  // PROCESS_SWITCH(SinglePhotonMC, processPHOS, "single photon with PHOS", false);
  // PROCESS_SWITCH(SinglePhotonMC, processEMC, "single photon with EMC", false);
  PROCESS_SWITCH(SinglePhotonMC, processGen, "analyze MC truth information", false);
  PROCESS_SWITCH(SinglePhotonMC, processDummy, "Dummy function", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& context)
{
  return WorkflowSpec{
    adaptAnalysisTask<SinglePhotonMC>(context, TaskName{"single-photon-mc"})};
}
