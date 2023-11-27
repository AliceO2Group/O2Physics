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

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>;
using MyCollisions_Cent = soa::Join<MyCollisions, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentNTPVs>; // centrality table has dependency on multiplicity table.

using MyCollisionsMC = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::McCollisionLabels>;
using MyCollisionsMC_Cent = soa::Join<MyCollisionsMC, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentNTPVs>; // centrality table has dependency on multiplicity table.

struct CreateEMEvent {
  Produces<o2::aod::EMReducedEvents> event;
  Produces<o2::aod::EMReducedEventsMult> event_mult;
  Produces<o2::aod::EMReducedEventsCent> event_cent;

  enum class EMEventType : int {
    kEvent = 0,
    kEvent_Cent = 1,
  };

  HistogramRegistry registry{"registry"};
  void init(o2::framework::InitContext&)
  {
    auto hEventCounter = registry.add<TH1>("hEventCounter", "hEventCounter", kTH1I, {{7, 0.5f, 7.5f}});
    hEventCounter->GetXaxis()->SetBinLabel(1, "all");
    hEventCounter->GetXaxis()->SetBinLabel(2, "sel8");
  }

  //! Please don't skip any event!
  template <bool isMC, EMEventType eventype, typename TEvents>
  void skimEvent(TEvents const& collisions, aod::BCs const&)
  {
    for (auto& collision : collisions) {
      if constexpr (isMC) {
        if (!collision.has_mcCollision()) {
          continue;
        }
      }

      registry.fill(HIST("hEventCounter"), 1);

      // auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      bool is_phoscpv_readout = collision.alias_bit(kTVXinPHOS);
      bool is_emc_readout = collision.alias_bit(kTVXinEMC);

      if (collision.sel8()) {
        registry.fill(HIST("hEventCounter"), 2);
      }

      uint64_t tag = collision.selection_raw();
      event(collision.globalIndex(), tag, collision.bc().runNumber(), collision.bc().triggerMask(), collision.sel8(),
            is_phoscpv_readout, is_emc_readout,
            collision.posX(), collision.posY(), collision.posZ(),
            collision.numContrib(), collision.collisionTime(), collision.collisionTimeRes());

      event_mult(collision.multTPC(),
                 collision.multFV0A(), collision.multFV0C(), collision.multFT0A(), collision.multFT0C(), collision.multFDDA(), collision.multFDDC(),
                 collision.multZNA(), collision.multZNC(),
                 collision.multTracklets(), collision.multNTracksPV(), collision.multNTracksPVeta1());

      if constexpr (eventype == EMEventType::kEvent) {
        event_cent(105.f, 105.f, 105.f, 105.f);
      } else if constexpr (eventype == EMEventType::kEvent_Cent) {
        event_cent(collision.centFT0M(), collision.centFT0A(), collision.centFT0C(), collision.centNTPV());
      } else {
        event_cent(105.f, 105.f, 105.f, 105.f);
      }
    } // end of collision loop
  }   // end of skimEvent

  void processEvent(MyCollisions const& collisions, aod::BCs const& bcs)
  {
    skimEvent<false, EMEventType::kEvent>(collisions, bcs);
  }
  PROCESS_SWITCH(CreateEMEvent, processEvent, "process event info", false);

  void processEventMC(MyCollisionsMC const& collisions, aod::BCs const& bcs)
  {
    skimEvent<true, EMEventType::kEvent>(collisions, bcs);
  }
  PROCESS_SWITCH(CreateEMEvent, processEventMC, "process event info", false);

  void processEvent_Cent(MyCollisions_Cent const& collisions, aod::BCs const& bcs)
  {
    skimEvent<false, EMEventType::kEvent_Cent>(collisions, bcs);
  }
  PROCESS_SWITCH(CreateEMEvent, processEvent_Cent, "process event info", false);

  void processEventMC_Cent(MyCollisionsMC_Cent const& collisions, aod::BCs const& bcs)
  {
    skimEvent<true, EMEventType::kEvent_Cent>(collisions, bcs);
  }
  PROCESS_SWITCH(CreateEMEvent, processEventMC_Cent, "process event info", false);

  void processDummy(aod::Collisions const& collisions) {}
  PROCESS_SWITCH(CreateEMEvent, processDummy, "processDummy", true);
};
struct AssociatePhotonToEMEvent {
  Produces<o2::aod::V0KFEMReducedEventIds> v0kfeventid;
  Produces<o2::aod::EMPrimaryElectronEMReducedEventIds> prmeleventid;
  Produces<o2::aod::EMPrimaryMuonEMReducedEventIds> prmmueventid;
  Produces<o2::aod::PHOSEMReducedEventIds> phoseventid;
  Produces<o2::aod::EMCEMReducedEventIds> emceventid;

  Produces<o2::aod::EMReducedEventsNgPCM> event_ng_pcm;
  Produces<o2::aod::EMReducedEventsNgPHOS> event_ng_phos;
  Produces<o2::aod::EMReducedEventsNgEMC> event_ng_emc;

  Preslice<aod::V0PhotonsKF> perCollision_pcm = aod::v0photonkf::collisionId;
  Preslice<aod::EMPrimaryElectrons> perCollision_el = aod::emprimaryelectron::collisionId;
  Preslice<aod::EMPrimaryMuons> perCollision_mu = aod::emprimarymuon::collisionId;
  Preslice<aod::PHOSClusters> perCollision_phos = aod::skimmedcluster::collisionId;
  Preslice<aod::SkimEMCClusters> perCollision_emc = aod::skimmedcluster::collisionId;

  bool doPCM = false, doDalitzEE = false, doDalitzMuMu = false, doPHOS = false, doEMC = false;
  void init(o2::framework::InitContext& context)
  {
    if (context.mOptions.get<bool>("processPCM")) {
      doPCM = true;
    }
    if (context.mOptions.get<bool>("processDalitzEE")) {
      doDalitzEE = true;
    }
    if (context.mOptions.get<bool>("processDalitzMuMu")) {
      doDalitzMuMu = true;
    }
    if (context.mOptions.get<bool>("processPHOS")) {
      doPHOS = true;
    }
    if (context.mOptions.get<bool>("processEMC")) {
      doEMC = true;
    }
  }

  template <typename TCollisions, typename TEventNg>
  void zero_padding(TCollisions const& collisions, TEventNg& event_ng)
  {
    int nc = collisions.size();
    for (int ic = 0; ic < nc; ic++) {
      event_ng(0);
    } // end of collision loop
  }

  template <typename TCollisions, typename TPhotons, typename TEventIds, typename TEventNg, typename TPreslice>
  void fillEventId_Ng(TCollisions const& collisions, TPhotons const& photons, TEventIds& eventIds, TEventNg& event_ng, TPreslice const& perCollision)
  {
    for (auto& collision : collisions) {
      auto photons_coll = photons.sliceBy(perCollision, collision.collisionId());
      int ng = photons_coll.size();
      for (int ig = 0; ig < ng; ig++) {
        eventIds(collision.globalIndex());
      } // end of photon loop
      event_ng(ng);
    } // end of collision loop
  }

  template <typename TCollisions, typename TPhotons, typename TEventIds, typename TPreslice>
  void fillEventId(TCollisions const& collisions, TPhotons const& photons, TEventIds& eventIds, TPreslice const& perCollision)
  {
    for (auto& collision : collisions) {
      auto photons_coll = photons.sliceBy(perCollision, collision.collisionId());
      int ng = photons_coll.size();
      for (int ig = 0; ig < ng; ig++) {
        eventIds(collision.globalIndex());
      } // end of photon loop
    }   // end of collision loop
  }

  // This struct is for both data and MC.
  // Note that reconstructed collisions without mc collisions are already rejected in CreateEMEvent in MC.

  void processPCM(aod::EMReducedEvents const& collisions, aod::V0PhotonsKF const& photons)
  {
    fillEventId_Ng(collisions, photons, v0kfeventid, event_ng_pcm, perCollision_pcm);
  }

  void processDalitzEE(aod::EMReducedEvents const& collisions, aod::EMPrimaryElectrons const& tracks)
  {
    fillEventId(collisions, tracks, prmeleventid, perCollision_el);
  }

  void processDalitzMuMu(aod::EMReducedEvents const& collisions, aod::EMPrimaryMuons const& tracks)
  {
    fillEventId(collisions, tracks, prmmueventid, perCollision_mu);
  }

  void processPHOS(aod::EMReducedEvents const& collisions, aod::PHOSClusters const& photons)
  {
    fillEventId_Ng(collisions, photons, phoseventid, event_ng_phos, perCollision_phos);
  }

  void processEMC(aod::EMReducedEvents const& collisions, aod::SkimEMCClusters const& photons)
  {
    fillEventId_Ng(collisions, photons, emceventid, event_ng_emc, perCollision_emc);
  }

  void processZeroPadding(aod::EMReducedEvents const& collisions)
  {
    if (!doPCM) {
      zero_padding(collisions, event_ng_pcm);
    }
    if (!doPHOS) {
      zero_padding(collisions, event_ng_phos);
    }
    if (!doEMC) {
      zero_padding(collisions, event_ng_emc);
    }
  }

  PROCESS_SWITCH(AssociatePhotonToEMEvent, processPCM, "process pcm-event indexing", false);
  PROCESS_SWITCH(AssociatePhotonToEMEvent, processDalitzEE, "process dalitzee-event indexing", false);
  PROCESS_SWITCH(AssociatePhotonToEMEvent, processDalitzMuMu, "process dalitzmumu-event indexing", false);
  PROCESS_SWITCH(AssociatePhotonToEMEvent, processPHOS, "process phos-event indexing", false);
  PROCESS_SWITCH(AssociatePhotonToEMEvent, processEMC, "process emc-event indexing", false);
  PROCESS_SWITCH(AssociatePhotonToEMEvent, processZeroPadding, "process zero padding.", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CreateEMEvent>(cfgc, TaskName{"create-emevent"}),
    adaptAnalysisTask<AssociatePhotonToEMEvent>(cfgc, TaskName{"associate-photon-to-emevent"}),
  };
}
