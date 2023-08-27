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
// This code produces reduced events for photon analyses.
//    Please write to: daiki.sekihata@cern.ch

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct createEMReducedEvent {
  Produces<o2::aod::EMReducedEvents> event;

  enum SubSystem {
    kPCM = 0x1,
    kPHOS = 0x2,
    kEMC = 0x4,
    kUndef = -1,
  };

  // Configurable for filter/cuts
  Configurable<int> minN_PCM{"minN_PCM", 0, "Minimum number of V0s for PCM. Events are saved if either minimum number condition is met"};
  Configurable<int> minN_PHOS{"minN_PHOS", 0, "Minimum number of clusters for PHOS. Events are saved if either minimum number condition is met"};
  Configurable<int> minN_EMC{"minN_EMC", 0, "Minimum number of clusters for EMCal. Events are saved if either minimum number condition is met"};

  HistogramRegistry registry{"registry"};

  Preslice<aod::V0PhotonsKF> perCollision_pcm = aod::v0photon::collisionId;
  Preslice<aod::PHOSClusters> perCollision_phos = aod::skimmedcluster::collisionId;
  Preslice<aod::SkimEMCClusters> perCollision_emc = aod::skimmedcluster::collisionId;

  void init(o2::framework::InitContext&)
  {
    auto hEventCounter = registry.add<TH1>("hEventCounter", "hEventCounter", kTH1I, {{6, 0.5f, 6.5f}});
    hEventCounter->GetXaxis()->SetBinLabel(1, "all");
    hEventCounter->GetXaxis()->SetBinLabel(2, "has > minN_PCM");
    hEventCounter->GetXaxis()->SetBinLabel(3, "has > minN_PHOS");
    hEventCounter->GetXaxis()->SetBinLabel(4, "has > minN_EMC");
    hEventCounter->GetXaxis()->SetBinLabel(5, "has > minN_any");
    hEventCounter->GetXaxis()->SetBinLabel(6, "sel8");
    registry.add<TH1>("hNGammas_PCM", ";#it{N}_{#gamma,PCM};#it{count}", kTH1I, {{21, -0.5, 20.5}});
    registry.add<TH1>("hNGammas_PHOS", ";#it{N}_{#gamma,PHOS};#it{count}", kTH1I, {{21, -0.5, 20.5}});
    registry.add<TH1>("hNGammas_EMC", ";#it{N}_{#gamma,EMC};#it{count}", kTH1I, {{21, -0.5, 20.5}});
  }

  using MyCollisions = soa::Join<aod::Collisions, aod::Mults, aod::EvSels, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentNTPVs>;

  template <uint8_t system, typename TPCMs, typename TPHOSs, typename TEMCs>
  void process(MyCollisions const& collisions, aod::BCs const&, TPCMs const& v0photons, TPHOSs const& phosclusters, TEMCs const& emcclusters)
  {
    for (auto& collision : collisions) {
      registry.fill(HIST("hEventCounter"), 1);

      // auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      bool is_phoscpv_readout = collision.alias_bit(kTVXinPHOS);
      bool is_emc_readout = collision.alias_bit(kTVXinEMC);

      int ng_pcm = 0;
      int ng_phos = 0;
      int ng_emc = 0;
      bool minN_any = false;

      if constexpr (static_cast<bool>(system & kPCM)) {
        auto v0photons_coll = v0photons.sliceBy(perCollision_pcm, collision.globalIndex());
        ng_pcm = v0photons_coll.size();
        registry.fill(HIST("hNGammas_PCM"), ng_pcm);
      }
      if constexpr (static_cast<bool>(system & kPHOS)) {
        auto phos_coll = phosclusters.sliceBy(perCollision_phos, collision.globalIndex());
        ng_phos = phos_coll.size();
        registry.fill(HIST("hNGammas_PHOS"), ng_phos);
      }
      if constexpr (static_cast<bool>(system & kEMC)) {
        auto emc_coll = emcclusters.sliceBy(perCollision_emc, collision.globalIndex());
        ng_emc = emc_coll.size();
        registry.fill(HIST("hNGammas_EMC"), ng_emc);
      }
      if (ng_pcm >= minN_PCM) {
        minN_any = true;
        registry.fill(HIST("hEventCounter"), 2);
      }
      if (ng_phos >= minN_PHOS) {
        minN_any = true;
        registry.fill(HIST("hEventCounter"), 3);
      }
      if (ng_emc >= minN_EMC) {
        minN_any = true;
        registry.fill(HIST("hEventCounter"), 4);
      }
      if (minN_any) {
        registry.fill(HIST("hEventCounter"), 5);
      }
      if (collision.sel8()) {
        registry.fill(HIST("hEventCounter"), 6);
      }

      // store event selection decisions
      uint64_t tag = collision.selection_raw();
      event(collision.globalIndex(), tag, collision.bc().runNumber(), collision.bc().triggerMask(), collision.sel8(),
            is_phoscpv_readout, is_emc_readout,
            collision.posX(), collision.posY(), collision.posZ(),
            collision.numContrib(), collision.collisionTime(), collision.collisionTimeRes(),
            collision.multTPC(), collision.multFV0A(), collision.multFV0C(), collision.multFT0A(), collision.multFT0C(),
            collision.multFDDA(), collision.multFDDC(), collision.multZNA(), collision.multZNC(), collision.multTracklets(), collision.multNTracksPV(), collision.multNTracksPVeta1(),
            collision.centFT0M(), collision.centFT0A(), collision.centFT0C(), collision.centNTPV(),
            ng_pcm, ng_phos, ng_emc); // ng is needed for event mixing to filter events that contain at least 1 photon.

    } // end of collision loop

  } // end of process

  void process_PCM(MyCollisions const& collisions, aod::BCs const& bcs, aod::V0PhotonsKF const& v0photons)
  {
    process<kPCM>(collisions, bcs, v0photons, nullptr, nullptr);
  }
  void process_PHOS(MyCollisions const& collisions, aod::BCs const& bcs, aod::PHOSClusters const& phosclusters)
  {
    process<kPHOS>(collisions, bcs, nullptr, phosclusters, nullptr);
  }
  void process_EMC(MyCollisions const& collisions, aod::BCs const& bcs, aod::SkimEMCClusters const& emcclusters)
  {
    process<kEMC>(collisions, bcs, nullptr, nullptr, emcclusters);
  }
  void process_PCM_PHOS(MyCollisions const& collisions, aod::BCs const& bcs, aod::V0PhotonsKF const& v0photons, aod::PHOSClusters const& phosclusters)
  {
    const uint8_t sysflag = kPCM | kPHOS;
    process<sysflag>(collisions, bcs, v0photons, phosclusters, nullptr);
  }
  void process_PCM_EMC(MyCollisions const& collisions, aod::BCs const& bcs, aod::V0PhotonsKF const& v0photons, aod::SkimEMCClusters const& emcclusters)
  {
    const uint8_t sysflag = kPCM | kEMC;
    process<sysflag>(collisions, bcs, v0photons, nullptr, emcclusters);
  }
  void process_PHOS_EMC(MyCollisions const& collisions, aod::BCs const& bcs, aod::PHOSClusters const& phosclusters, aod::SkimEMCClusters const& emcclusters)
  {
    const uint8_t sysflag = kPHOS | kEMC;
    process<sysflag>(collisions, bcs, nullptr, phosclusters, emcclusters);
  }
  void process_PCM_PHOS_EMC(MyCollisions const& collisions, aod::BCs const& bcs, aod::V0PhotonsKF const& v0photons, aod::PHOSClusters const& phosclusters, aod::SkimEMCClusters const& emcclusters)
  {
    const uint8_t sysflag = kPCM | kPHOS | kEMC;
    process<sysflag>(collisions, bcs, v0photons, phosclusters, emcclusters);
  }

  void processDummy(MyCollisions const& collisions) {}

  PROCESS_SWITCH(createEMReducedEvent, process_PCM, "create em event table for PCM", false);
  PROCESS_SWITCH(createEMReducedEvent, process_PHOS, "create em event table for PHOS", false);
  PROCESS_SWITCH(createEMReducedEvent, process_EMC, "create em event table for EMCal", false);
  PROCESS_SWITCH(createEMReducedEvent, process_PCM_PHOS, "create em event table for PCM, PHOS", false);
  PROCESS_SWITCH(createEMReducedEvent, process_PCM_EMC, "create em event table for PCM, EMCal", false);
  PROCESS_SWITCH(createEMReducedEvent, process_PHOS_EMC, "create em event table for PHOS, EMCal", false);
  PROCESS_SWITCH(createEMReducedEvent, process_PCM_PHOS_EMC, "create em event table for PCM, PHOS, EMCal", false);
  PROCESS_SWITCH(createEMReducedEvent, processDummy, "processDummy", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<createEMReducedEvent>(cfgc, TaskName{"create-emreduced-event"})};
}
