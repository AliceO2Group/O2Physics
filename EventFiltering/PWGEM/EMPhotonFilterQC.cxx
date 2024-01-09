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

// \brief Quick QC task for CEFP for EM photon
// \author daiki.sekihata@cern.ch

#include "Math/Vector4D.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/CaloClusters.h"
#include "DataFormatsPHOS/TriggerRecord.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/Core/CutsLibrary.h"
#include "EventFiltering/filterTables.h"
#include "Framework/HistogramRegistry.h"

using namespace o2;
using namespace o2::soa;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::PhotonFilters>;
using MyCollision = MyCollisions::iterator;

struct EMPhotonFilterQC {
  HistogramRegistry registry{"registry"};

  void init(o2::framework::InitContext&)
  {
    addhistograms();
  }

  void addhistograms()
  {
    auto hEventCounter = registry.add<TH1>("hEventCounter", "hEventCounter", kTH1F, {{20, 0.5f, 20.5f}});
    hEventCounter->GetXaxis()->SetBinLabel(1, "all");
    hEventCounter->GetXaxis()->SetBinLabel(2, "sel8 && |Z_{vtx}| < 10 cm");
    hEventCounter->GetXaxis()->SetBinLabel(3, "PCM High p_{T} photon");
    hEventCounter->GetXaxis()->SetBinLabel(4, "PCM High p_{T} photon && sel8 && |Z_{vtx}| < 10 cm");
    hEventCounter->GetXaxis()->SetBinLabel(5, "PCM MB calibration");
    hEventCounter->GetXaxis()->SetBinLabel(6, "PCM MB calibration && sel8 && |Z_{vtx}| < 10 cm");
    hEventCounter->GetXaxis()->SetBinLabel(7, "PCM #eta #rightarrow ee#gamma");
    hEventCounter->GetXaxis()->SetBinLabel(8, "PCM #eta #rightarrow ee#gamma && sel8 && |Z_{vtx}| < 10 cm");
    hEventCounter->GetXaxis()->SetBinLabel(9, "PCM #eta #rightarrow #gamma#gamma");
    hEventCounter->GetXaxis()->SetBinLabel(10, "PCM #eta #rightarrow #gamma#gamma && sel8 && |Z_{vtx}| < 10 cm");
    hEventCounter->GetXaxis()->SetBinLabel(11, "PCM and ee");
    hEventCounter->GetXaxis()->SetBinLabel(12, "PCM and ee && sel8 && |Z_{vtx}| < 10 cm");

    registry.add<TH1>("PCM_HighPt/hPt", "pT of PCM photon;p_{T,#gamma} (GeV/c)", kTH1F, {{200, 0, 20}});
    registry.add<TH2>("PCM_HighPt/hEtaPhi", "#eta vs. #varphi of PCM photon", kTH2F, {{72, 0, TMath::TwoPi()}, {40, -2, +2}});
    registry.add<TH2>("PCM_MBCalib/hMeePt", "mass ee#gamma;m_{ee} (GeV/c^{2});p_{T,ee} (GeV/c)", kTH2F, {{50, 0.f, 0.5f}, {100, 0, 10}});
    registry.add<TH2>("PCM_MBCalib/hMeegPt", "mass ee#gamma;m_{ee#gamma} (GeV/c^{2});p_{T,#gamma} (GeV/c)", kTH2F, {{200, 0.f, 0.4f}, {100, 0, 10}});
    registry.add<TH2>("PCM_EtaDalitz/hMeePt", "mass ee;m_{ee} (GeV/c^{2});p_{T,ee} (GeV/c)", kTH2F, {{400, 0.f, 4.f}, {100, 0, 10}});
    registry.add<TH2>("PCM_EtaDalitz/hMeegPt", "mass ee#gamma;m_{ee#gamma} (GeV/c^{2});p_{T,ee#gamma} (GeV/c)", kTH2F, {{400, 0.f, 0.8f}, {100, 0, 10}});
    registry.add<TH2>("PCM_EtaGG/hMggPt", "mass ee#gamma;m_{#gamma#gamma} (GeV/c^{2});p_{T,#gamma#gamma} (GeV/c)", kTH2F, {{400, 0.f, 0.8f}, {100, 0, 10}});
    registry.add<TH2>("PCM_EE/hMeePt", "mass ee;m_{ee} (GeV/c^{2});p_{T,ee} (GeV/c)", kTH2F, {{400, 0.f, 4.f}, {100, 0, 10}});
    registry.add<TH2>("PCM_EE/hMeegPt", "mass ee#gamma;m_{ee#gamma} (GeV/c^{2});p_{T,ee#gamma} (GeV/c)", kTH2F, {{400, 0.f, 0.8f}, {100, 0, 10}});
  }

  Preslice<aod::V0PhotonsKF> perCollision_pcm = aod::v0photonkf::collisionId;
  // Preslice<aod::DalitzEEs> perCollision_ee = aod::dalitzee::collisionId;
  // void processPCM(MyCollisions const& collisions, aod::V0PhotonsKF const& v0photons, aod::DalitzEEs const& dielectrons)
  void processPCM(MyCollisions const& collisions, aod::V0PhotonsKF const& v0photons)
  {
    for (auto& collision : collisions) {
      registry.fill(HIST("hEventCounter"), 1);
      if (collision.sel8() && abs(collision.posZ()) < 10.f) {
        registry.fill(HIST("hEventCounter"), 2);
      }
      auto v0photons_coll = v0photons.sliceBy(perCollision_pcm, collision.globalIndex());
      // auto dielectrons_coll = dielectrons.sliceBy(perCollision_ee, collision.globalIndex());

      if (collision.hasPCMHighPtPhoton()) {
        registry.fill(HIST("hEventCounter"), 3);
        if (collision.sel8() && abs(collision.posZ()) < 10.f) {
          registry.fill(HIST("hEventCounter"), 4);
        }
        for (auto& v0photon : v0photons_coll) {
          registry.fill(HIST("PCM_HighPt/hPt"), v0photon.pt());
          registry.fill(HIST("PCM_HighPt/hEtaPhi"), v0photon.phi(), v0photon.eta());
        } // end of v0 photon loop
      }

      // if (collision.hasPCMMatCalib()) {
      //   registry.fill(HIST("hEventCounter"), 5);
      //   if (collision.sel8() && abs(collision.posZ()) < 10.f) {
      //     registry.fill(HIST("hEventCounter"), 6);
      //   }
      //   for (auto& [g1, g2] : combinations(CombinationsFullIndexPolicy(v0photons_coll, dielectrons_coll))) {
      //     registry.fill(HIST("PCM_MBCalib/hMeePt"), g2.mass(), g2.pt());
      //     ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
      //     ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), g2.mass());
      //     ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
      //     registry.fill(HIST("PCM_MBCalib/hMeegPt"), v12.M(), g1.pt());
      //   } // end of dielectron-photon pair loop
      // }

      // if (collision.hasPCMEtaDalitz()) {
      //   registry.fill(HIST("hEventCounter"), 7);
      //   if (collision.sel8() && abs(collision.posZ()) < 10.f) {
      //     registry.fill(HIST("hEventCounter"), 8);
      //   }
      //   for (auto& dielectron : dielectrons_coll) {
      //     registry.fill(HIST("PCM_EtaDalitz/hMeePt"), dielectron.mass(), dielectron.pt());
      //   }
      //   for (auto& [g1, g2] : combinations(CombinationsFullIndexPolicy(v0photons_coll, dielectrons_coll))) {
      //     ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
      //     ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), g2.mass());
      //     ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
      //     registry.fill(HIST("PCM_EtaDalitz/hMeegPt"), v12.M(), v12.Pt());
      //   } // end of dielectron-photon pair loop
      // }

      // if (collision.hasPCMEtaGG()) {
      //   registry.fill(HIST("hEventCounter"), 9);
      //   if (collision.sel8() && abs(collision.posZ()) < 10.f) {
      //     registry.fill(HIST("hEventCounter"), 10);
      //   }
      //   for (auto& [g1, g2] : combinations(CombinationsStrictlyUpperIndexPolicy(v0photons_coll, v0photons_coll))) {
      //     ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
      //     ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
      //     ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
      //     registry.fill(HIST("PCM_EtaGG/hMggPt"), v12.M(), v12.Pt());
      //   } // end of dielectron-photon pair loop
      // }

      // if (collision.hasPCMandEE()) {
      //   registry.fill(HIST("hEventCounter"), 11);
      //   if (collision.sel8() && abs(collision.posZ()) < 10.f) {
      //     registry.fill(HIST("hEventCounter"), 12);
      //   }

      //   for (auto& dielectron : dielectrons_coll) {
      //     registry.fill(HIST("PCM_EE/hMeePt"), dielectron.mass(), dielectron.pt());
      //   }
      //   for (auto& [g1, g2] : combinations(CombinationsFullIndexPolicy(v0photons_coll, dielectrons_coll))) {
      //     ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
      //     ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), g2.mass());
      //     ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
      //     registry.fill(HIST("PCM_EE/hMeegPt"), v12.M(), v12.Pt());
      //   } // end of dielectron-photon pair loop
      // }

    } // end of collision loop
  }

  void processPHOS(MyCollisions const& collisions) {}
  void processEMC(MyCollisions const& collisions) {}

  PROCESS_SWITCH(EMPhotonFilterQC, processPCM, "Process PCM software trigger QC", true);
  PROCESS_SWITCH(EMPhotonFilterQC, processPHOS, "Process PHOS software trigger QC", false);
  PROCESS_SWITCH(EMPhotonFilterQC, processEMC, "Process EMC software trigger QC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{
    adaptAnalysisTask<EMPhotonFilterQC>(cfg, TaskName{"em-photon-filter-qc"})};
}
