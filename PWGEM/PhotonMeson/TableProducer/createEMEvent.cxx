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

#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"

#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

using MyBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>;
using MyCollisions_Cent = soa::Join<MyCollisions, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentNTPVs>; // centrality table has dependency on multiplicity table.
using MyCollisions_Cent_Qvec = soa::Join<MyCollisions_Cent, aod::QvectorFT0Cs, aod::QvectorFT0As, aod::QvectorFT0Ms, aod::QvectorFV0As, aod::QvectorBPoss, aod::QvectorBNegs>;

using MyCollisionsMC = soa::Join<MyCollisions, aod::McCollisionLabels>;
using MyCollisionsMC_Cent = soa::Join<MyCollisionsMC, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentNTPVs>; // centrality table has dependency on multiplicity table.
using MyCollisionsMC_Cent_Qvec = soa::Join<MyCollisionsMC_Cent, aod::QvectorFT0Cs, aod::QvectorFT0As, aod::QvectorFT0Ms, aod::QvectorFV0As, aod::QvectorBPoss, aod::QvectorBNegs>;

struct CreateEMEvent {
  Produces<o2::aod::EMEvents> event;
  Produces<o2::aod::EMEventsCov> eventcov;
  Produces<o2::aod::EMEventsMult> event_mult;
  Produces<o2::aod::EMEventsCent> event_cent;
  Produces<o2::aod::EMEventsQvec> event_qvec;

  enum class EMEventType : int {
    kEvent = 0,
    kEvent_Cent = 1,
    kEvent_Cent_Qvec = 2,
  };

  // CCDB options
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};
  Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};

  HistogramRegistry registry{"registry"};
  void init(o2::framework::InitContext&)
  {
    auto hEventCounter = registry.add<TH1>("hEventCounter", "hEventCounter", kTH1I, {{7, 0.5f, 7.5f}});
    hEventCounter->GetXaxis()->SetBinLabel(1, "all");
    hEventCounter->GetXaxis()->SetBinLabel(2, "sel8");
  }

  int mRunNumber;
  float d_bz;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  template <typename TBC>
  void initCCDB(TBC const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    // In case override, don't proceed, please - no CCDB access required
    if (d_bz_input > -990) {
      d_bz = d_bz_input;
      o2::parameters::GRPMagField grpmag;
      if (fabs(d_bz) > 1e-5) {
        grpmag.setL3Current(30000.f / (d_bz / 5.0f));
      }
      mRunNumber = bc.runNumber();
      return;
    }

    auto run3grp_timestamp = bc.timestamp();
    o2::parameters::GRPObject* grpo = 0x0;
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (!skipGRPOquery)
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    if (grpo) {
      // Fetch magnetic field from ccdb for current collision
      d_bz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      // Fetch magnetic field from ccdb for current collision
      d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    }
    mRunNumber = bc.runNumber();
  }

  PresliceUnsorted<MyCollisions> preslice_collisions_per_bc = o2::aod::evsel::foundBCId;
  std::unordered_map<uint64_t, int> map_ncolls_per_bc;

  //! Please don't skip any event!
  template <bool isMC, EMEventType eventype, typename TCollisions, typename TBCs>
  void skimEvent(TCollisions const& collisions, TBCs const& bcs)
  {
    // first count the number of collisions per bc
    for (auto& bc : bcs) {
      auto collisions_per_bc = collisions.sliceBy(preslice_collisions_per_bc, bc.globalIndex());
      map_ncolls_per_bc[bc.globalIndex()] = collisions_per_bc.size();
      // LOGF(info, "bc-loop | bc.globalIndex() = %d , collisions_per_bc.size() = %d", bc.globalIndex(), collisions_per_bc.size());
    }

    for (auto& collision : collisions) {
      if constexpr (isMC) {
        if (!collision.has_mcCollision()) {
          continue;
        }
      }

      auto bc = collision.template foundBC_as<TBCs>();
      initCCDB(bc);

      // LOGF(info, "collision-loop | bc.globalIndex() = %d, ncolls_per_bc = %d", bc.globalIndex(), map_ncolls_per_bc[bc.globalIndex()]);
      registry.fill(HIST("hEventCounter"), 1);

      if (collision.sel8()) {
        registry.fill(HIST("hEventCounter"), 2);
      }

      // uint64_t tag = collision.selection_raw();
      event(collision.globalIndex(), bc.runNumber(), bc.globalBC(), collision.sel8(), collision.alias_raw(), collision.selection_raw(), bc.timestamp(), map_ncolls_per_bc[bc.globalIndex()],
            collision.posX(), collision.posY(), collision.posZ(),
            collision.numContrib(), collision.trackOccupancyInTimeRange());

      eventcov(collision.covXX(), collision.covXY(), collision.covXZ(), collision.covYY(), collision.covYZ(), collision.covZZ(), collision.chi2());

      event_mult(collision.multFT0A(), collision.multFT0C(), collision.multTPC(), collision.multNTracksPV(), collision.multNTracksPVeta1(), collision.multNTracksPVetaHalf());

      if constexpr (eventype == EMEventType::kEvent) {
        event_cent(105.f, 105.f, 105.f, 105.f);
        event_qvec(999.f, 999.f, 999.f, 999.f, 999.f, 999.f,
                   // 999.f, 999.f,
                   999.f, 999.f, 999.f, 999.f,
                   999.f, 999.f, 999.f, 999.f, 999.f, 999.f,
                   // 999.f, 999.f,
                   999.f, 999.f, 999.f, 999.f);
      } else if constexpr (eventype == EMEventType::kEvent_Cent) {
        event_cent(collision.centFT0M(), collision.centFT0A(), collision.centFT0C(), collision.centNTPV());
        event_qvec(999.f, 999.f, 999.f, 999.f, 999.f, 999.f,
                   // 999.f, 999.f,
                   999.f, 999.f, 999.f, 999.f,
                   999.f, 999.f, 999.f, 999.f, 999.f, 999.f,
                   // 999.f, 999.f,
                   999.f, 999.f, 999.f, 999.f);
      } else if constexpr (eventype == EMEventType::kEvent_Cent_Qvec) {
        event_cent(collision.centFT0M(), collision.centFT0A(), collision.centFT0C(), collision.centNTPV());
        event_qvec(collision.qvecFT0MRe(), collision.qvecFT0MIm(), collision.qvecFT0ARe(), collision.qvecFT0AIm(), collision.qvecFT0CRe(), collision.qvecFT0CIm(),
                   // collision.qvecFV0ARe(), collision.qvecFV0AIm(),
                   collision.qvecBPosRe(), collision.qvecBPosIm(), collision.qvecBNegRe(), collision.qvecBNegIm(),
                   999.f, 999.f, 999.f, 999.f, 999.f, 999.f,
                   // 999.f, 999.f,
                   999.f, 999.f, 999.f, 999.f); // as of 20240416, only 2nd harmonics is supported.
      } else {
        event_cent(105.f, 105.f, 105.f, 105.f);
        event_qvec(999.f, 999.f, 999.f, 999.f, 999.f, 999.f,
                   // 999.f, 999.f,
                   999.f, 999.f, 999.f, 999.f,
                   999.f, 999.f, 999.f, 999.f, 999.f, 999.f,
                   // 999.f, 999.f,
                   999.f, 999.f, 999.f, 999.f);
      }
    } // end of collision loop
    map_ncolls_per_bc.clear();
  } // end of skimEvent

  void processEvent(MyCollisions const& collisions, MyBCs const& bcs)
  {
    skimEvent<false, EMEventType::kEvent>(collisions, bcs);
  }
  PROCESS_SWITCH(CreateEMEvent, processEvent, "process event info", false);

  void processEventMC(MyCollisionsMC const& collisions, MyBCs const& bcs)
  {
    skimEvent<true, EMEventType::kEvent>(collisions, bcs);
  }
  PROCESS_SWITCH(CreateEMEvent, processEventMC, "process event info", false);

  void processEvent_Cent(MyCollisions_Cent const& collisions, MyBCs const& bcs)
  {
    skimEvent<false, EMEventType::kEvent_Cent>(collisions, bcs);
  }
  PROCESS_SWITCH(CreateEMEvent, processEvent_Cent, "process event info", false);

  void processEvent_Cent_Qvec(MyCollisions_Cent_Qvec const& collisions, MyBCs const& bcs)
  {
    skimEvent<false, EMEventType::kEvent_Cent_Qvec>(collisions, bcs);
  }
  PROCESS_SWITCH(CreateEMEvent, processEvent_Cent_Qvec, "process event info", false);

  void processEventMC_Cent(MyCollisionsMC_Cent const& collisions, MyBCs const& bcs)
  {
    skimEvent<true, EMEventType::kEvent_Cent>(collisions, bcs);
  }
  PROCESS_SWITCH(CreateEMEvent, processEventMC_Cent, "process event info", false);

  void processEventMC_Cent_Qvec(MyCollisionsMC_Cent_Qvec const& collisions, MyBCs const& bcs)
  {
    skimEvent<true, EMEventType::kEvent_Cent_Qvec>(collisions, bcs);
  }
  PROCESS_SWITCH(CreateEMEvent, processEventMC_Cent_Qvec, "process event info", false);

  void processDummy(aod::Collisions const&) {}
  PROCESS_SWITCH(CreateEMEvent, processDummy, "processDummy", true);
};
struct AssociatePhotonToEMEvent {
  Produces<o2::aod::V0KFEMEventIds> v0kfeventid;
  Produces<o2::aod::EMPrimaryElectronEMEventIds> prmeleventid;
  Produces<o2::aod::EMPrimaryMuonEMEventIds> prmmueventid;
  Produces<o2::aod::PHOSEMEventIds> phoseventid;
  Produces<o2::aod::EMCEMEventIds> emceventid;

  Produces<o2::aod::EMEventsNgPCM> event_ng_pcm;
  Produces<o2::aod::EMEventsNgPHOS> event_ng_phos;
  Produces<o2::aod::EMEventsNgEMC> event_ng_emc;

  Preslice<aod::V0PhotonsKF> perCollision_pcm = aod::v0photonkf::collisionId;
  PresliceUnsorted<aod::EMPrimaryElectrons> perCollision_el = aod::emprimaryelectron::collisionId;
  PresliceUnsorted<aod::EMPrimaryMuons> perCollision_mu = aod::emprimarymuon::collisionId;
  Preslice<aod::PHOSClusters> perCollision_phos = aod::skimmedcluster::collisionId;
  Preslice<aod::SkimEMCClusters> perCollision_emc = aod::skimmedcluster::collisionId;

  bool doPCM = false, doElectron = false, doFwdMuon = false, doPHOS = false, doEMC = false;
  void init(o2::framework::InitContext& context)
  {
    if (context.mOptions.get<bool>("processPCM")) {
      doPCM = true;
    }
    if (context.mOptions.get<bool>("processElectron")) {
      doElectron = true;
    }
    if (context.mOptions.get<bool>("processFwdMuon")) {
      doFwdMuon = true;
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

  void processPCM(aod::EMEvents const& collisions, aod::V0PhotonsKF const& photons)
  {
    fillEventId_Ng(collisions, photons, v0kfeventid, event_ng_pcm, perCollision_pcm);
  }

  void processElectron(aod::EMEvents const& collisions, aod::EMPrimaryElectrons const& tracks)
  {
    fillEventId(collisions, tracks, prmeleventid, perCollision_el);
  }

  void processFwdMuon(aod::EMEvents const& collisions, aod::EMPrimaryMuons const& tracks)
  {
    fillEventId(collisions, tracks, prmmueventid, perCollision_mu);
  }

  void processPHOS(aod::EMEvents const& collisions, aod::PHOSClusters const& photons)
  {
    fillEventId_Ng(collisions, photons, phoseventid, event_ng_phos, perCollision_phos);
  }

  void processEMC(aod::EMEvents const& collisions, aod::SkimEMCClusters const& photons)
  {
    fillEventId_Ng(collisions, photons, emceventid, event_ng_emc, perCollision_emc);
  }

  void processZeroPadding(aod::EMEvents const& collisions)
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
  PROCESS_SWITCH(AssociatePhotonToEMEvent, processElectron, "process dalitzee-event indexing", false);
  PROCESS_SWITCH(AssociatePhotonToEMEvent, processFwdMuon, "process forward muon indexing", false);
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
