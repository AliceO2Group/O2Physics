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
// This code runs loop over PHOS clusters for PHOS QC.
//    Please write to: daiki.sekihata@cern.ch

#include "PWGEM/PhotonMeson/Core/CutsLibrary.h"
#include "PWGEM/PhotonMeson/Core/HistogramsLibrary.h"
#include "PWGEM/PhotonMeson/Core/PHOSPhotonCut.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"

#include "Common/CCDB/TriggerAliases.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>

#include <THashList.h>
#include <TString.h>

#include <cstdlib>
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

using MyCollisions = soa::Join<aod::EMEvents_004, aod::EMEventsAlias, aod::EMEventsMult_000, aod::EMEventsCent_000>;
using MyCollision = MyCollisions::iterator;

struct phosQC {

  Configurable<std::string> fConfigPHOSCuts{"cfgPHOSCuts", "test02,test03,test05", "Comma separated list of phos photon cuts"};
  std::vector<PHOSPhotonCut> fPHOSCuts;

  OutputObj<THashList> fOutputEvent{"Event"};
  OutputObj<THashList> fOutputCluster{"Cluster"};
  THashList* fMainList = new THashList();

  // static constexpr std::string_view ambtracktypes[2] = {"NonAmb", "Amb"};
  void addhistograms()
  {
    fMainList->SetOwner(true);
    fMainList->SetName("fMainList");

    // create sub lists first.
    o2::aod::pwgem::photon::histogram::AddHistClass(fMainList, "Event");
    THashList* list_ev = reinterpret_cast<THashList*>(fMainList->FindObject("Event"));
    o2::aod::pwgem::photon::histogram::DefineHistograms(list_ev, "Event");

    o2::aod::pwgem::photon::histogram::AddHistClass(fMainList, "Cluster");
    THashList* list_cluster = reinterpret_cast<THashList*>(fMainList->FindObject("Cluster"));

    for (const auto& cut : fPHOSCuts) {
      const char* cutname = cut.GetName();
      o2::aod::pwgem::photon::histogram::AddHistClass(list_cluster, cutname);
    }

    // for Clusters
    for (auto& cut : fPHOSCuts) {
      std::string_view cutname = cut.GetName();
      THashList* list = reinterpret_cast<THashList*>(fMainList->FindObject("Cluster")->FindObject(cutname.data()));
      o2::aod::pwgem::photon::histogram::DefineHistograms(list, "Cluster", "PHOS");
    }
  }

  void DefineCuts()
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

  void init(InitContext&)
  {
    DefineCuts();
    addhistograms(); // please call this after DefinCuts();

    fOutputEvent.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("Event")));
    fOutputCluster.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("Cluster")));
  }

  Preslice<aod::PHOSClusters> perCollision = aod::skimmedcluster::collisionId;

  void processQC(MyCollisions const& collisions, aod::PHOSClusters const& clusters)
  {
    THashList* list_ev = static_cast<THashList*>(fMainList->FindObject("Event"));
    THashList* list_cluster = static_cast<THashList*>(fMainList->FindObject("Cluster"));

    for (auto& collision : collisions) {

      if (!collision.alias_bit(triggerAliases::kTVXinPHOS)) {
        continue;
      }

      reinterpret_cast<TH1F*>(fMainList->FindObject("Event")->FindObject("hCollisionCounter"))->Fill(1.0);
      if (!collision.sel8()) {
        continue;
      }
      reinterpret_cast<TH1F*>(fMainList->FindObject("Event")->FindObject("hCollisionCounter"))->Fill(2.0);

      if (collision.numContrib() < 0.5) {
        continue;
      }
      reinterpret_cast<TH1F*>(fMainList->FindObject("Event")->FindObject("hCollisionCounter"))->Fill(3.0);

      if (abs(collision.posZ()) > 10.0) {
        continue;
      }
      reinterpret_cast<TH1F*>(fMainList->FindObject("Event")->FindObject("hCollisionCounter"))->Fill(4.0);
      reinterpret_cast<TH1F*>(fMainList->FindObject("Event")->FindObject("hZvtx"))->Fill(collision.posZ());
      o2::aod::pwgem::photon::histogram::FillHistClass<EMHistType::kEvent>(list_ev, "", collision);

      auto clusters_per_coll = clusters.sliceBy(perCollision, collision.collisionId());
      for (const auto& cut : fPHOSCuts) {
        THashList* list_cluster_cut = static_cast<THashList*>(list_cluster->FindObject(cut.GetName()));
        int ng = 0;
        for (auto& cluster : clusters_per_coll) {

          if (cut.IsSelected(cluster)) {
            o2::aod::pwgem::photon::histogram::FillHistClass<EMHistType::kPHOSCluster>(list_cluster_cut, "", cluster);
            ng++;
          }
        } // end of v0 loop
        reinterpret_cast<TH1F*>(fMainList->FindObject("Cluster")->FindObject(cut.GetName())->FindObject("hNgamma"))->Fill(ng);
      } // end of cut loop
    } // end of collision loop
  } // end of process

  void processDummy(MyCollisions const&) {}

  PROCESS_SWITCH(phosQC, processQC, "run PHOS QC", false);
  PROCESS_SWITCH(phosQC, processDummy, "Dummy function", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<phosQC>(cfgc, TaskName{"phos-qc"})};
}
