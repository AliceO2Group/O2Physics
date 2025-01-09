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

// #include <string>
// #include <unordered_map>
#include <chrono>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"

#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "Common/CCDB/TriggerAliases.h"

#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

using MyBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;
using MyQvectors = soa::Join<aod::QvectorFT0CVecs, aod::QvectorFT0AVecs, aod::QvectorFT0MVecs, aod::QvectorBPosVecs, aod::QvectorBNegVecs, aod::QvectorBTotVecs>;

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>;
using MyCollisions_Cent = soa::Join<MyCollisions, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>; // centrality table has dependency on multiplicity table.
using MyCollisions_Cent_Qvec = soa::Join<MyCollisions_Cent, MyQvectors>;

using MyCollisionsMC = soa::Join<MyCollisions, aod::McCollisionLabels>;
using MyCollisionsMC_Cent = soa::Join<MyCollisionsMC, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>; // centrality table has dependency on multiplicity table.
using MyCollisionsMC_Cent_Qvec = soa::Join<MyCollisionsMC_Cent, MyQvectors>;

struct CreateEMEvent {
  Produces<o2::aod::EMEvents> event;
  Produces<o2::aod::EMEventsCov> eventcov;
  Produces<o2::aod::EMEventsMult> event_mult;
  Produces<o2::aod::EMEventsCent> event_cent;
  Produces<o2::aod::EMEventsQvec> event_qvec;
  Produces<o2::aod::EMEventsWeight> event_weights;

  enum class EMEventType : int {
    kEvent = 0,
    kEvent_Cent = 1,
    kEvent_Cent_Qvec = 2,
    kEvent_JJ = 3,
  };

  // CCDB options
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};
  Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};
  Configurable<bool> applyEveSel_at_skimming{"applyEveSel_at_skimming", false, "flag to apply minimal event selection at the skimming level"};
  Configurable<bool> needEMCTrigger{"needEMCTrigger", false, "flag to only save events which have kTVXinEMC trigger bit. To reduce PbPb derived data size"};
  Configurable<bool> needPHSTrigger{"needPHSTrigger", false, "flag to only save events which have kTVXinPHOS trigger bit. To reduce PbPb derived data size"};

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

  // PresliceUnsorted<MyCollisions> preslice_collisions_per_bc = o2::aod::evsel::foundBCId;
  // std::unordered_map<uint64_t, int> map_ncolls_per_bc;

  template <bool isMC, EMEventType eventype, typename TCollisions, typename TBCs>
  void skimEvent(TCollisions const& collisions, TBCs const& bcs)
  {
    auto start = std::chrono::high_resolution_clock::now();
    // first count the number of collisions per bc
    // for (auto& bc : bcs) {
    //   auto collisions_per_bc = collisions.sliceBy(preslice_collisions_per_bc, bc.globalIndex());
    //   map_ncolls_per_bc[bc.globalIndex()] = collisions_per_bc.size();
    //   // LOGF(info, "bc-loop | bc.globalIndex() = %d , collisions_per_bc.size() = %d", bc.globalIndex(), collisions_per_bc.size());
    // }

    for (auto& collision : collisions) {
      if constexpr (isMC) {
        if (!collision.has_mcCollision()) {
          continue;
        }
      }

      auto bc = collision.template foundBC_as<TBCs>();
      initCCDB(bc);

      if (applyEveSel_at_skimming && (!collision.selection_bit(o2::aod::evsel::kIsTriggerTVX) || !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder))) {
        continue;
      }
      if (needEMCTrigger && !collision.alias_bit(kTVXinEMC)) {
        continue;
      }
      if (needPHSTrigger && !collision.alias_bit(kTVXinPHOS)) {
        continue;
      }

      float qDefault = 999.f; // default value for q vectors if not obtained

      // LOGF(info, "collision-loop | bc.globalIndex() = %d, ncolls_per_bc = %d", bc.globalIndex(), map_ncolls_per_bc[bc.globalIndex()]);
      registry.fill(HIST("hEventCounter"), 1);

      if (collision.sel8()) {
        registry.fill(HIST("hEventCounter"), 2);
      }

      // uint64_t tag = collision.selection_raw();
      event(collision.globalIndex(), bc.runNumber(), bc.globalBC(), collision.alias_raw(), collision.selection_raw(), bc.timestamp(),
            collision.posX(), collision.posY(), collision.posZ(),
            collision.numContrib(), collision.trackOccupancyInTimeRange(), collision.ft0cOccupancyInTimeRange());

      eventcov(collision.covXX(), collision.covXY(), collision.covXZ(), collision.covYY(), collision.covYZ(), collision.covZZ(), collision.chi2());

      event_mult(collision.multFT0A(), collision.multFT0C(), collision.multNTracksPV(), collision.multNTracksPVeta1(), collision.multNTracksPVetaHalf());

      if constexpr (eventype != EMEventType::kEvent_JJ) {
        event_weights(1.f);
      }

      if constexpr (eventype == EMEventType::kEvent) {
        event_cent(105.f, 105.f, 105.f);
        event_qvec(qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault,
                   qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault);
      } else if constexpr (eventype == EMEventType::kEvent_Cent) {
        event_cent(collision.centFT0M(), collision.centFT0A(), collision.centFT0C());
        event_qvec(qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault,
                   qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault);
      } else if constexpr (eventype == EMEventType::kEvent_Cent_Qvec) {
        event_cent(collision.centFT0M(), collision.centFT0A(), collision.centFT0C());
        const size_t qvecSize = collision.qvecFT0CReVec().size();
        if (qvecSize >= 2) { // harmonics 2,3
          event_qvec(collision.qvecFT0MReVec()[0], collision.qvecFT0MImVec()[0], collision.qvecFT0AReVec()[0], collision.qvecFT0AImVec()[0], collision.qvecFT0CReVec()[0], collision.qvecFT0CImVec()[0], collision.qvecBPosReVec()[0], collision.qvecBPosImVec()[0], collision.qvecBNegReVec()[0], collision.qvecBNegImVec()[0], collision.qvecBTotReVec()[0], collision.qvecBTotImVec()[0],
                    collision.qvecFT0MReVec()[1], collision.qvecFT0MImVec()[1], collision.qvecFT0AReVec()[1], collision.qvecFT0AImVec()[1], collision.qvecFT0CReVec()[1], collision.qvecFT0CImVec()[1], collision.qvecBPosReVec()[1], collision.qvecBPosImVec()[1], collision.qvecBNegReVec()[1], collision.qvecBNegImVec()[1], collision.qvecBTotReVec()[1], collision.qvecBTotImVec()[1]);
        } else if (qvecSize >= 1) { // harmonics 2
          event_qvec(collision.qvecFT0MReVec()[0], collision.qvecFT0MImVec()[0], collision.qvecFT0AReVec()[0], collision.qvecFT0AImVec()[0], collision.qvecFT0CReVec()[0], collision.qvecFT0CImVec()[0], collision.qvecBPosReVec()[0], collision.qvecBPosImVec()[0], collision.qvecBNegReVec()[0], collision.qvecBNegImVec()[0], collision.qvecBTotReVec()[0], collision.qvecBTotImVec()[0],
                   qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault);
        }
      } else {
        event_cent(105.f, 105.f, 105.f);
        event_qvec(qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault,
                   qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault);
      }
    } // end of collision loop
    // map_ncolls_per_bc.clear();
    auto end = std::chrono::high_resolution_clock::now();
    LOG(info) << "Time elapsed: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << " ns";
  } // end of skimEvent

  void fillEventWeights(MyCollisionsMC const& collisions, aod::McCollisions const&, MyBCs const&)
  {
    for (auto& collision : collisions) {
      if (!collision.has_mcCollision()) {
        continue;
      }

      auto bc = collision.template foundBC_as<MyBCs>();
      initCCDB(bc);

      if (applyEveSel_at_skimming && (!collision.selection_bit(o2::aod::evsel::kIsTriggerTVX) || !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder))) {
        continue;
      }
      auto mcCollision = collision.mcCollision();
      event_weights(mcCollision.weight());
    }
  }

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

  void processEventJJMC(MyCollisionsMC const& collisions, aod::McCollisions const& mcCollisions, MyBCs const& bcs)
  {
    skimEvent<true, EMEventType::kEvent_JJ>(collisions, bcs);
    fillEventWeights(collisions, mcCollisions, bcs);
  }
  PROCESS_SWITCH(CreateEMEvent, processEventJJMC, "process event info", false);

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

  Preslice<aod::V0PhotonsKF> perCollision_pcm = aod::v0photonkf::collisionId;
  PresliceUnsorted<aod::EMPrimaryElectronsFromDalitz> perCollision_el = aod::emprimaryelectron::collisionId;
  PresliceUnsorted<aod::EMPrimaryMuons> perCollision_mu = aod::emprimarymuon::collisionId;
  Preslice<aod::PHOSClusters> perCollision_phos = aod::skimmedcluster::collisionId;
  Preslice<aod::SkimEMCClusters> perCollision_emc = aod::skimmedcluster::collisionId;

  void init(o2::framework::InitContext&) {}

  template <typename TCollisions, typename TPhotons, typename TEventIds, typename TPreslice>
  void fillEventId(TCollisions const& collisions, TPhotons const& photons, TEventIds& eventIds, TPreslice const& perCollision)
  {
    for (auto& collision : collisions) {
      auto photons_coll = photons.sliceBy(perCollision, collision.collisionId());
      int ng = photons_coll.size();
      for (int ig = 0; ig < ng; ig++) {
        eventIds(collision.globalIndex());
      } // end of photon loop
    } // end of collision loop
  }

  // This struct is for both data and MC.
  // Note that reconstructed collisions without mc collisions are already rejected in CreateEMEvent in MC.

  void processPCM(aod::EMEvents const& collisions, aod::V0PhotonsKF const& photons)
  {
    fillEventId(collisions, photons, v0kfeventid, perCollision_pcm);
  }

  void processElectronFromDalitz(aod::EMEvents const& collisions, aod::EMPrimaryElectronsFromDalitz const& tracks)
  {
    fillEventId(collisions, tracks, prmeleventid, perCollision_el);
  }

  void processPHOS(aod::EMEvents const& collisions, aod::PHOSClusters const& photons)
  {
    fillEventId(collisions, photons, phoseventid, perCollision_phos);
  }

  void processEMC(aod::EMEvents const& collisions, aod::SkimEMCClusters const& photons)
  {
    fillEventId(collisions, photons, emceventid, perCollision_emc);
  }

  void processDummy(aod::EMEvents const&) {}

  PROCESS_SWITCH(AssociatePhotonToEMEvent, processPCM, "process pcm-event indexing", false);
  PROCESS_SWITCH(AssociatePhotonToEMEvent, processElectronFromDalitz, "process dalitzee-event indexing", false);
  PROCESS_SWITCH(AssociatePhotonToEMEvent, processPHOS, "process phos-event indexing", false);
  PROCESS_SWITCH(AssociatePhotonToEMEvent, processEMC, "process emc-event indexing", false);
  PROCESS_SWITCH(AssociatePhotonToEMEvent, processDummy, "process dummy", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CreateEMEvent>(cfgc, TaskName{"create-emevent-photon"}),
    adaptAnalysisTask<AssociatePhotonToEMEvent>(cfgc, TaskName{"associate-photon-to-emevent"}),
  };
}
