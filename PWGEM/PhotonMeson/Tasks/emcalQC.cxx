// Copyright 2019-2024 CERN and copyright holders of ALICE O2.
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
// This code runs loop over EMCal clusters for EMCal QC.
//    Please write to: nicolas.strangmann@cern.ch

#include <array>
#include "TString.h"
#include "THashList.h"
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
#include "PWGEM/PhotonMeson/Core/EMCPhotonCut.h"
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

struct emcalQC {

  Configurable<bool> requireCaloReadout{"requireCaloReadout", true, "Require calorimeters readout when analyzing EMCal/PHOS"};
  Configurable<std::string> fConfigEMCCuts{"cfgEMCCuts", "custom,standard,nocut", "Comma separated list of EMCal photon cuts"};
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
  Configurable<bool> fDo2DQC{"Do2DQC", true, "Flag to output 2D QC histograms displaying the energy dependence on the second axis"};

  std::vector<EMCPhotonCut> fEMCCuts;

  OutputObj<THashList> fOutputEvent{"Event"};
  OutputObj<THashList> fOutputCluster{"Cluster"};
  THashList* fMainList = new THashList();

  void addhistograms()
  {
    fMainList->SetOwner(true);
    fMainList->SetName("fMainList");

    // create sub lists first.
    o2::aod::pwgem::photon::histogram::AddHistClass(fMainList, "Event");
    THashList* list_ev = reinterpret_cast<THashList*>(fMainList->FindObject("Event"));
    o2::aod::pwgem::photon::histogram::DefineHistograms(list_ev, "Event");

    list_ev->Add(new TH1F("hNEvents", "hNEvents", 6, 0.5f, 6.5f));
    reinterpret_cast<TH2F*>(list_ev->FindObject("hNEvents"))->GetXaxis()->SetBinLabel(1, "all cols");
    reinterpret_cast<TH2F*>(list_ev->FindObject("hNEvents"))->GetXaxis()->SetBinLabel(2, "sel8");
    reinterpret_cast<TH2F*>(list_ev->FindObject("hNEvents"))->GetXaxis()->SetBinLabel(3, "emc readout");
    reinterpret_cast<TH2F*>(list_ev->FindObject("hNEvents"))->GetXaxis()->SetBinLabel(4, "1+ Contributor");
    reinterpret_cast<TH2F*>(list_ev->FindObject("hNEvents"))->GetXaxis()->SetBinLabel(5, "z<10cm");
    reinterpret_cast<TH2F*>(list_ev->FindObject("hNEvents"))->GetXaxis()->SetBinLabel(6, "unique col");

    o2::aod::pwgem::photon::histogram::AddHistClass(fMainList, "Cluster");
    THashList* list_cluster = reinterpret_cast<THashList*>(fMainList->FindObject("Cluster"));

    for (const auto& cut : fEMCCuts) {
      const char* cutname = cut.GetName();
      o2::aod::pwgem::photon::histogram::AddHistClass(list_cluster, cutname);
    }

    // for Clusters
    for (auto& cut : fEMCCuts) {
      std::string_view cutname = cut.GetName();
      THashList* list = reinterpret_cast<THashList*>(fMainList->FindObject("Cluster")->FindObject(cutname.data()));
      o2::aod::pwgem::photon::histogram::DefineHistograms(list, "Cluster", fDo2DQC ? "2D_EMC" : "EMC");
    }
  }

  void DefineCuts()
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
    LOGF(info, "Number of EMC cuts = %d", fEMCCuts.size());
  }

  void init(InitContext&)
  {
    DefineCuts();
    addhistograms(); // please call this after DefinCuts();

    fOutputEvent.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("Event")));
    fOutputCluster.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("Cluster")));
  }

  Preslice<aod::SkimEMCClusters> perCollision = aod::skimmedcluster::collisionId;

  void processQC(MyCollisions const& collisions, aod::SkimEMCClusters const& clusters)
  {
    THashList* list_ev = static_cast<THashList*>(fMainList->FindObject("Event"));
    THashList* list_cluster = static_cast<THashList*>(fMainList->FindObject("Cluster"));

    for (auto& collision : collisions) {
      reinterpret_cast<TH1F*>(fMainList->FindObject("Event")->FindObject("hNEvents"))->Fill(1.0);
      if (!collision.sel8()) {
        continue;
      }
      reinterpret_cast<TH1F*>(fMainList->FindObject("Event")->FindObject("hNEvents"))->Fill(2.0);
      if (collision.alias_bit(triggerAliases::kTVXinEMC) == 0 && requireCaloReadout) {
        continue;
      }
      reinterpret_cast<TH1F*>(fMainList->FindObject("Event")->FindObject("hNEvents"))->Fill(3.0);

      if (collision.numContrib() < 0.5) {
        continue;
      }
      reinterpret_cast<TH1F*>(fMainList->FindObject("Event")->FindObject("hNEvents"))->Fill(4.0);

      if (abs(collision.posZ()) > 10.0) {
        continue;
      }
      reinterpret_cast<TH1F*>(fMainList->FindObject("Event")->FindObject("hNEvents"))->Fill(5.0);
      reinterpret_cast<TH1F*>(fMainList->FindObject("Event")->FindObject("hZvtx"))->Fill(collision.posZ());

      if (collision.ncollsPerBC() != 1) { // Check that the collision is unique (the only one in the bc)
        continue;
      }
      reinterpret_cast<TH1F*>(fMainList->FindObject("Event")->FindObject("hNEvents"))->Fill(6.0);

      o2::aod::pwgem::photon::histogram::FillHistClass<EMHistType::kEvent>(list_ev, "", collision); // Fill event histograms

      auto clusters_per_coll = clusters.sliceBy(perCollision, collision.collisionId());
      for (const auto& cut : fEMCCuts) {
        THashList* list_cluster_cut = static_cast<THashList*>(list_cluster->FindObject(cut.GetName()));
        int ng = 0;
        for (auto& cluster : clusters_per_coll) {
          // Fill the cluster properties before applying any cuts

          // Apply cuts one by one and fill in hClusterQualityCuts histogram
          reinterpret_cast<TH2F*>(fMainList->FindObject("Cluster")->FindObject(cut.GetName())->FindObject("hClusterQualityCuts"))->Fill(0., cluster.e());
          auto track = nullptr;

          // Define two boleans to see, whether the cluster "survives" the EMC cluster cuts to later check, whether the cuts in this task align with the ones in EMCPhotonCut.h:
          bool survivesIsSelectedEMCalCuts = true;                    // Survives "manual" cuts listed in this task
          bool survivesIsSelectedCuts = cut.IsSelected<int>(cluster); // Survives the cutlist defines in EMCPhotonCut.h, which is also used in the Pi0Eta task

          for (int icut = 0; icut < static_cast<int>(EMCPhotonCut::EMCPhotonCuts::kNCuts); icut++) { // Loop through different cut observables
            EMCPhotonCut::EMCPhotonCuts specificcut = static_cast<EMCPhotonCut::EMCPhotonCuts>(icut);
            if (!cut.IsSelectedEMCal(specificcut, cluster, track)) { // Check whether cluster passes this cluster requirement, if not, fill why in the next row
              reinterpret_cast<TH1F*>(fMainList->FindObject("Cluster")->FindObject(cut.GetName())->FindObject("hClusterQualityCuts"))->Fill(icut + 1, cluster.e());
              survivesIsSelectedEMCalCuts = false;
            }
          }

          if (survivesIsSelectedCuts != survivesIsSelectedEMCalCuts) {
            LOGF(info, "Cummulative application of IsSelectedEMCal cuts does not equal the IsSelected result");
          }

          if (survivesIsSelectedCuts) {
            o2::aod::pwgem::photon::histogram::FillHistClass<EMHistType::kEMCCluster>(list_cluster_cut, fDo2DQC ? "2D" : "1D", cluster);
            reinterpret_cast<TH2F*>(fMainList->FindObject("Cluster")->FindObject(cut.GetName())->FindObject("hClusterQualityCuts"))->Fill(7., cluster.e());
            ng++;
          }
        }
        reinterpret_cast<TH1F*>(fMainList->FindObject("Cluster")->FindObject(cut.GetName())->FindObject("hNgamma"))->Fill(ng);
      } // end of cut loop
    }   // end of collision loop
  }     // end of process

  void processDummy(MyCollisions const&) {}

  PROCESS_SWITCH(emcalQC, processQC, "run EMCal QC", false);
  PROCESS_SWITCH(emcalQC, processDummy, "Dummy function", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<emcalQC>(cfgc, TaskName{"emcal-qc"})};
}
