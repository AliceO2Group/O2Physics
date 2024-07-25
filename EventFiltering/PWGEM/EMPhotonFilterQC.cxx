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
    auto hEventCounter = registry.add<TH1>("hEventCounter", "hEventCounter", kTH1F, {{10, 0.5f, 10.5f}});
    hEventCounter->GetXaxis()->SetBinLabel(1, "all");
    hEventCounter->GetXaxis()->SetBinLabel(2, "No TFB");
    hEventCounter->GetXaxis()->SetBinLabel(3, "sel8 && |Z_{vtx}| < 10 cm && No TFB");
    hEventCounter->GetXaxis()->SetBinLabel(4, "PCM High p_{T} photon");
    hEventCounter->GetXaxis()->SetBinLabel(5, "PCM High p_{T} photon && sel8 && |Z_{vtx}| < 10 cm && No TFB");
    hEventCounter->GetXaxis()->SetBinLabel(6, "PCM and dielectron");
    hEventCounter->GetXaxis()->SetBinLabel(7, "PCM and dielectron && sel8 && |Z_{vtx}| < 10 cm && No TFB");

    registry.add<TH1>("PCM_HighPt/hPt", "pT of PCM photon;p_{T,#gamma} (GeV/c)", kTH1F, {{200, 0, 20}});
    registry.add<TH2>("PCM_HighPt/hEtaPhi", "#eta vs. #varphi of PCM photon", kTH2F, {{72, 0, TMath::TwoPi()}, {40, -2, +2}});
    registry.add<TH2>("PCM_HighPt/hMeegvsPtg", "meeg vs. pTg;m_{ee#gamma} (GeV/c^{2});p_{T,#gamma} (GeV/c)", kTH2F, {{200, 0, 0.8}, {200, 0, 20}});
    registry.add<TH2>("PCM_EE/hMeePt", "mass ee;m_{ee} (GeV/c^{2});p_{T,ee} (GeV/c)", kTH2F, {{400, 0.f, 4.f}, {100, 0, 10}});
    registry.add<TH2>("PCM_EE/hMeePhiV", "mass ee;#varphi_{V} (rad.);m_{ee} (GeV/c^{2})", kTH2F, {{180, 0.f, TMath::Pi()}, {100, 0, 0.1}});
    registry.add<TH2>("PCM_EE/hTPCdEdx", "n sigma TPC el vs. pin;p_{in} (GeV/c);n #sigma_{e}^{TPC}", kTH2F, {{1000, 0.f, 10.f}, {100, -5, 5}});
    registry.add<TH2>("PCM_EE/hMeegPt", "mass ee#gamma;m_{ee#gamma} (GeV/c^{2});p_{T,ee#gamma} (GeV/c)", kTH2F, {{400, 0.f, 0.8f}, {100, 0, 10}});
  }

  // Preslice<aod::V0PhotonsKF> perCollision_pcm = aod::v0photonkf::collisionId;
  Preslice<aod::DalitzEEs> perCollision_ee = aod::dalitzee::collisionId;
  Preslice<aod::EMSwtInfosPCM> perCollision_swt_pcm = aod::pwgem::photon::swtinfo::collisionId;
  Preslice<aod::EMSwtInfosPair> perCollision_swt_pair = aod::pwgem::photon::swtinfo::collisionId;
  std::vector<uint64_t> stored_dielectronIds;
  std::vector<uint64_t> stored_trackIds;

  void processPCM(MyCollisions const& collisions, aod::V0PhotonsKF const& /*v0photons*/, aod::V0Legs const&, aod::DalitzEEs const& dielectrons, aod::EMPrimaryElectrons const&, aod::EMSwtInfosPCM const& swt_pcm, aod::EMSwtInfosPair const& swt_pair)
  {
    stored_dielectronIds.reserve(swt_pair.size());
    stored_trackIds.reserve(swt_pair.size() * 2);

    for (auto& collision : collisions) {
      registry.fill(HIST("hEventCounter"), 1);
      if (collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
        registry.fill(HIST("hEventCounter"), 2);
      }
      if (collision.sel8() && abs(collision.posZ()) < 10.f && collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
        registry.fill(HIST("hEventCounter"), 3);
      }
      // auto v0photons_coll = v0photons.sliceBy(perCollision_pcm, collision.globalIndex());
      auto dielectrons_coll = dielectrons.sliceBy(perCollision_ee, collision.globalIndex());
      auto triggers_pcm_coll = swt_pcm.sliceBy(perCollision_swt_pcm, collision.globalIndex());
      auto triggers_eeg_coll = swt_pair.sliceBy(perCollision_swt_pair, collision.globalIndex());

      if (collision.hasPCMHighPtPhoton()) {
        registry.fill(HIST("hEventCounter"), 4);
        if (collision.sel8() && abs(collision.posZ()) < 10.f && collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
          registry.fill(HIST("hEventCounter"), 5);
        }

        for (auto& trigger_pcm_id : triggers_pcm_coll) {
          auto v0photon = trigger_pcm_id.template triggerV0PhotonHighPt_as<aod::V0PhotonsKF>();
          registry.fill(HIST("PCM_HighPt/hPt"), v0photon.pt());
          registry.fill(HIST("PCM_HighPt/hEtaPhi"), v0photon.phi(), v0photon.eta());
        } // end of v0 photon loop

        for (auto& [trigger_pcm_id, g2] : combinations(CombinationsFullIndexPolicy(triggers_pcm_coll, dielectrons_coll))) {
          auto g1 = trigger_pcm_id.template triggerV0PhotonHighPt_as<aod::V0PhotonsKF>();
          ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
          ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), g2.mass());
          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          registry.fill(HIST("PCM_HighPt/hMeegvsPtg"), v12.M(), v1.Pt());
        } // end of eeg pair loop for high pT
      }

      if (collision.hasPCMandEE()) {
        registry.fill(HIST("hEventCounter"), 6);
        if (collision.sel8() && abs(collision.posZ()) < 10.f && collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
          registry.fill(HIST("hEventCounter"), 7);
        }

        for (auto& trigger_eeg_id : triggers_eeg_coll) {
          auto v0photon = trigger_eeg_id.template triggerV0PhotonPair_as<aod::V0PhotonsKF>();
          auto dielectron = trigger_eeg_id.template triggerDielectronPair_as<aod::DalitzEEs>();
          auto pos_pv = dielectron.template posTrack_as<aod::EMPrimaryElectrons>();
          auto ele_pv = dielectron.template negTrack_as<aod::EMPrimaryElectrons>();

          if (std::find(stored_trackIds.begin(), stored_trackIds.end(), pos_pv.globalIndex()) == stored_trackIds.end()) {
            registry.fill(HIST("PCM_EE/hTPCdEdx"), pos_pv.tpcInnerParam(), pos_pv.tpcNSigmaEl());
            stored_trackIds.emplace_back(pos_pv.globalIndex());
          }
          if (std::find(stored_trackIds.begin(), stored_trackIds.end(), ele_pv.globalIndex()) == stored_trackIds.end()) {
            registry.fill(HIST("PCM_EE/hTPCdEdx"), ele_pv.tpcInnerParam(), ele_pv.tpcNSigmaEl());
            stored_trackIds.emplace_back(ele_pv.globalIndex());
          }
          if (std::find(stored_dielectronIds.begin(), stored_dielectronIds.end(), dielectron.globalIndex()) == stored_dielectronIds.end()) {
            registry.fill(HIST("PCM_EE/hMeePt"), dielectron.mass(), dielectron.pt());
            registry.fill(HIST("PCM_EE/hMeePhiV"), dielectron.phiv(), dielectron.mass());
            stored_dielectronIds.emplace_back(dielectron.globalIndex());
          }
          ROOT::Math::PtEtaPhiMVector v1(v0photon.pt(), v0photon.eta(), v0photon.phi(), 0.);
          ROOT::Math::PtEtaPhiMVector v2(dielectron.pt(), dielectron.eta(), dielectron.phi(), dielectron.mass());
          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          registry.fill(HIST("PCM_EE/hMeegPt"), v12.M(), v12.Pt());
        } // end of eeg pair loop
      }
    } // end of collision loop
    stored_trackIds.clear();
    stored_trackIds.shrink_to_fit();
    stored_dielectronIds.clear();
    stored_dielectronIds.shrink_to_fit();
  }

  void processPHOS(MyCollisions const&) {}
  void processEMC(MyCollisions const&) {}
  void processDummy(MyCollisions const&) {}

  PROCESS_SWITCH(EMPhotonFilterQC, processPCM, "Process PCM software trigger QC", false);
  PROCESS_SWITCH(EMPhotonFilterQC, processPHOS, "Process PHOS software trigger QC", false);
  PROCESS_SWITCH(EMPhotonFilterQC, processEMC, "Process EMC software trigger QC", false);
  PROCESS_SWITCH(EMPhotonFilterQC, processDummy, "Process dummy", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{
    adaptAnalysisTask<EMPhotonFilterQC>(cfg, TaskName{"em-photon-filter-qc"})};
}
