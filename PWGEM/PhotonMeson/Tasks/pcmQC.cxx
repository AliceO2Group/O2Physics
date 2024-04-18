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
// This code runs loop over v0 photons for PCM QC.
//    Please write to: daiki.sekihata@cern.ch

#include <array>
#include "TString.h"
#include "THashList.h"
#include "TDirectory.h"
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
#include "PWGEM/PhotonMeson/Utils/PCMUtilities.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/Core/CutsLibrary.h"
#include "PWGEM/PhotonMeson/Core/HistogramsLibrary.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::photon;
using std::array;

using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent>;
using MyCollision = MyCollisions::iterator;

using MyV0Photons = soa::Join<aod::V0PhotonsKF, aod::V0KFEMEventIds>;
using MyV0Photon = MyV0Photons::iterator;

struct PCMQC {
  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};

  Configurable<std::string> fConfigPCMCuts{"cfgPCMCuts", "qc,qc_ITSTPC,qc_ITSonly,qc_TPConly,wwire_ib,nocut", "Comma separated list of v0 photon cuts"};
  std::vector<V0PhotonCut> fPCMCuts;

  Configurable<std::string> fConfigEMEventCut{"cfgEMEventCut", "minbias", "em event cut"}; // only 1 event cut per wagon
  EMEventCut fEMEventCut;
  static constexpr std::string_view event_types[2] = {"before", "after"};

  OutputObj<THashList> fOutputEvent{"Event"};
  OutputObj<THashList> fOutputV0Leg{"V0Leg"};
  OutputObj<THashList> fOutputV0{"V0"};
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

    o2::aod::pwgem::photon::histogram::AddHistClass(fMainList, "V0Leg");
    THashList* list_v0leg = reinterpret_cast<THashList*>(fMainList->FindObject("V0Leg"));

    o2::aod::pwgem::photon::histogram::AddHistClass(fMainList, "V0");
    THashList* list_v0 = reinterpret_cast<THashList*>(fMainList->FindObject("V0"));

    for (const auto& cut : fPCMCuts) {
      const char* cutname = cut.GetName();
      o2::aod::pwgem::photon::histogram::AddHistClass(list_v0leg, cutname);
      o2::aod::pwgem::photon::histogram::AddHistClass(list_v0, cutname);
    }

    // for single tracks
    for (auto& cut : fPCMCuts) {
      std::string_view cutname = cut.GetName();
      THashList* list = reinterpret_cast<THashList*>(fMainList->FindObject("V0Leg")->FindObject(cutname.data()));
      o2::aod::pwgem::photon::histogram::DefineHistograms(list, "V0Leg");
    }

    // for V0s
    for (auto& cut : fPCMCuts) {
      std::string_view cutname = cut.GetName();
      THashList* list = reinterpret_cast<THashList*>(fMainList->FindObject("V0")->FindObject(cutname.data()));
      o2::aod::pwgem::photon::histogram::DefineHistograms(list, "V0");
    }
  }

  void DefineCuts()
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

  void init(InitContext&)
  {
    DefineCuts();
    addhistograms(); // please call this after DefinCuts();
    TString ev_cut_name = fConfigEMEventCut.value;
    fEMEventCut = *eventcuts::GetCut(ev_cut_name.Data());
    fOutputEvent.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("Event")));
    fOutputV0Leg.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("V0Leg")));
    fOutputV0.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("V0")));
  }

  Preslice<MyV0Photons> perCollision = aod::v0photonkf::emeventId;
  Partition<MyCollisions> grouped_collisions = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax); // this goes to same event.

  void processQC(MyCollisions const&, MyV0Photons const& v0photons, aod::V0Legs const&)
  {
    THashList* list_ev_before = static_cast<THashList*>(fMainList->FindObject("Event")->FindObject(event_types[0].data()));
    THashList* list_ev_after = static_cast<THashList*>(fMainList->FindObject("Event")->FindObject(event_types[1].data()));
    THashList* list_v0 = static_cast<THashList*>(fMainList->FindObject("V0"));
    THashList* list_v0leg = static_cast<THashList*>(fMainList->FindObject("V0Leg"));

    for (auto& collision : grouped_collisions) {
      const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
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

      auto V0Photons_coll = v0photons.sliceBy(perCollision, collision.globalIndex());
      for (const auto& cut : fPCMCuts) {
        THashList* list_v0_cut = static_cast<THashList*>(list_v0->FindObject(cut.GetName()));
        THashList* list_v0leg_cut = static_cast<THashList*>(list_v0leg->FindObject(cut.GetName()));

        int nv0 = 0;
        for (auto& v0 : V0Photons_coll) {
          auto pos = v0.posTrack_as<aod::V0Legs>();
          auto ele = v0.negTrack_as<aod::V0Legs>();
          if (cut.IsSelected<aod::V0Legs>(v0)) {
            o2::aod::pwgem::photon::histogram::FillHistClass<EMHistType::kV0>(list_v0_cut, "", v0);
            nv0++;

            for (auto& leg : {pos, ele}) {
              o2::aod::pwgem::photon::histogram::FillHistClass<EMHistType::kV0Leg>(list_v0leg_cut, "", leg);
            }
          }
        } // end of v0 loop
        reinterpret_cast<TH1F*>(fMainList->FindObject("V0")->FindObject(cut.GetName())->FindObject("hNgamma"))->Fill(nv0);
      } // end of cut loop
    }   // end of collision loop
  }     // end of process

  void processDummy(MyCollisions const&) {}

  PROCESS_SWITCH(PCMQC, processQC, "run PCM QC", true);
  PROCESS_SWITCH(PCMQC, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<PCMQC>(cfgc, TaskName{"pcm-qc"})};
}
