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
/// \file createEMEventDilepton.cxx
/// \brief This code produces reduced events for dilepton analyses.
/// \author Daiki Sekihata, daiki.sekihata@cern.ch

#include "PWGEM/Dilepton/DataModel/dileptonTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Qvectors.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>

#include <algorithm>
#include <cstdint>
#include <iterator>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

using MyBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;
// using MyMults = soa::Join<aod::Mults, /*aod::MultsGlobal,*/ aod::FT0MultZeqs, aod::PVMultZeqs/*, aod::GlobalMultZeqs*/>;
using MyCents = soa::Join<aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentNTPVs, aod::CentNGlobals>; // centrality table has dependency on multiplicity table.
using MyQvectors = soa::Join<aod::QvectorFT0CVecs, aod::QvectorFT0AVecs, aod::QvectorFT0MVecs, aod::QvectorFV0AVecs, aod::QvectorBPosVecs, aod::QvectorBNegVecs, aod::QvectorBTotVecs>;

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::EMEvSels, aod::EMEoIs, aod::Mults>;
using MyCollisions_Cent = soa::Join<MyCollisions, MyCents>;
using MyCollisions_Cent_Qvec = soa::Join<MyCollisions, MyCents, MyQvectors>;

using MyCollisionsWithSWT = soa::Join<MyCollisions, aod::EMSWTriggerBitsTMP>;
using MyCollisionsWithSWT_Cent = soa::Join<MyCollisionsWithSWT, MyCents>;
using MyCollisionsWithSWT_Cent_Qvec = soa::Join<MyCollisionsWithSWT, MyCents, MyQvectors>;

using MyCollisionsMC = soa::Join<MyCollisions, aod::McCollisionLabels>;
using MyCollisionsMC_Cent = soa::Join<MyCollisionsMC, MyCents>;
using MyCollisionsMC_Cent_Qvec = soa::Join<MyCollisionsMC, MyCents, MyQvectors>;

struct CreateEMEventDilepton {
  Produces<o2::aod::EMBCs_001> embc;
  Produces<o2::aod::EMEvents> event;
  Produces<o2::aod::EMEventsXY> eventXY;
  // Produces<o2::aod::EMEventsCov> eventcov;
  Produces<o2::aod::EMEventsMult> event_mult;
  Produces<o2::aod::EMEventsCent> event_cent;
  // Produces<o2::aod::EMEventsQvec> event_qvec;
  Produces<o2::aod::EMEventsQvec2> event_qvec2;
  Produces<o2::aod::EMEventsQvec3> event_qvec3;
  Produces<o2::aod::EMEventNormInfos> event_norm_info;

  enum class EMEventType : int {
    kEvent = 0,
    kEvent_Cent = 1,
    kEvent_Cent_Qvec = 2,
  };

  // // CCDB options
  // Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  HistogramRegistry registry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  void init(o2::framework::InitContext&)
  {
    // ccdb->setURL(ccdburl);
    // ccdb->setCaching(true);
    // ccdb->setLocalObjectValidityChecking();
    // ccdb->setFatalWhenNull(false);

    auto hEventCounter = registry.add<TH1>("hEventCounter", "hEventCounter", kTH1D, {{7, 0.5f, 7.5f}});
    hEventCounter->GetXaxis()->SetBinLabel(1, "all");
    hEventCounter->GetXaxis()->SetBinLabel(2, "sel8");

    registry.add("hCentFT0M", "hCentFT0M;centrality FT0M (%);Number of collisions", kTH1F, {{110, 0, 110}}, false);
    registry.add("hCentFT0A", "hCentFT0A;centrality FT0A (%);Number of collisions", kTH1F, {{110, 0, 110}}, false);
    registry.add("hCentFT0C", "hCentFT0C;centrality FT0C (%);Number of collisions", kTH1F, {{110, 0, 110}}, false);
    registry.add("hCentNTPV", "hCentNTPV;centrality NTPV (%);Number of collisions", kTH1F, {{110, 0, 110}}, false);
    registry.add("hCentNGlobal", "hCentNGlobal;centrality NGlobal (%);Number of collisions", kTH1F, {{110, 0, 110}}, false);
  }

  ~CreateEMEventDilepton() {}

  int mRunNumber{0};
  // Service<o2::ccdb::BasicCCDBManager> ccdb;

  template <bool isMC, bool isTriggerAnalysis, EMEventType eventtype, typename TCollisions, typename TBCs>
  void skimEvent(TCollisions const& collisions, TBCs const& bcs)
  {
    for (const auto& bc : bcs) {
      if (bc.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
        embc(o2::aod::emevsel::reduceSelectionBit(bc), bc.rct_raw()); // TVX is fired.
        // embc(bc.selection_raw(), bc.rct_raw()); // TVX is fired.
      }
    } // end of bc loop

    for (const auto& collision : collisions) {
      if constexpr (isMC) {
        if (!collision.has_mcCollision()) {
          continue;
        }
      }
      registry.fill(HIST("hEventCounter"), 1);

      // auto bc = collision.template foundBC_as<TBCs>();
      auto bc = collision.template bc_as<TBCs>(); // use this for Zorro

      if (collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
        int8_t posZint8 = static_cast<int8_t>(collision.posZ() * 2.f);
        if (posZint8 == 0) {
          if (collision.posZ() < 0.f) {
            posZint8 = -1;
          } else {
            posZint8 = +1;
          }
        }
        if constexpr (eventtype == EMEventType::kEvent) {
          event_norm_info(o2::aod::emevsel::reduceSelectionBit(collision), collision.rct_raw(), posZint8, static_cast<uint8_t>(105.f + 110.f), static_cast<uint8_t>(105.f + 110.f), static_cast<uint8_t>(105.f + 110.f) /*, static_cast<uint8_t>(105.f + 110.f)*/);
        } else if constexpr (eventtype == EMEventType::kEvent_Cent || eventtype == EMEventType::kEvent_Cent_Qvec) {
          uint8_t centFT0Muint8 = collision.centFT0M() < 1.f ? static_cast<uint8_t>(collision.centFT0M() * 100.f) : static_cast<uint8_t>(collision.centFT0M() + 110.f);
          uint8_t centFT0Cuint8 = collision.centFT0C() < 1.f ? static_cast<uint8_t>(collision.centFT0C() * 100.f) : static_cast<uint8_t>(collision.centFT0C() + 110.f);
          uint8_t centNTPVuint8 = collision.centNTPV() < 1.f ? static_cast<uint8_t>(collision.centNTPV() * 100.f) : static_cast<uint8_t>(collision.centNTPV() + 110.f);
          // uint8_t centNGlobaluint8 = collision.centNGlobal() < 1.f ? static_cast<uint8_t>(collision.centNGlobal() * 100.f) : static_cast<uint8_t>(collision.centNGlobal() + 110.f);
          event_norm_info(o2::aod::emevsel::reduceSelectionBit(collision), collision.rct_raw(), posZint8, centFT0Muint8, centFT0Cuint8, centNTPVuint8 /*, centNGlobaluint8*/);
        } else {
          event_norm_info(o2::aod::emevsel::reduceSelectionBit(collision), collision.rct_raw(), posZint8, static_cast<uint8_t>(105.f + 110.f), static_cast<uint8_t>(105.f + 110.f), static_cast<uint8_t>(105.f + 110.f) /*, static_cast<uint8_t>(105.f + 110.f)*/);
        }
      }

      if (!collision.isSelected()) { // minimal cut for MB
        continue;
      }

      if (!collision.isEoI()) { // events with at least 1 lepton for data reduction.
        continue;
      }

      if constexpr (isTriggerAnalysis) {
        if (collision.swtaliastmp_raw() == 0) {
          continue;
        }
      }

      registry.fill(HIST("hEventCounter"), 2);

      event(collision.globalIndex(), bc.runNumber(), bc.globalBC(), o2::aod::emevsel::reduceSelectionBit(collision), collision.rct_raw(), bc.timestamp(),
            collision.posZ(),
            collision.numContrib(), collision.trackOccupancyInTimeRange(), collision.ft0cOccupancyInTimeRange());

      eventXY(collision.posX(), collision.posY());

      // eventcov(collision.covXX(), collision.covXY(), collision.covXZ(), collision.covYY(), collision.covYZ(), collision.covZZ(), collision.chi2());

      event_mult(collision.multFT0A(), collision.multFT0C(), collision.multNTracksPV() /*, collision.multNTracksGlobal()*/);

      if constexpr (eventtype == EMEventType::kEvent) {
        event_cent(105.f, 105.f, 105.f, 105.f);
        event_qvec2(999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f);
        event_qvec3(999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f);
        // event_qvec(
        //   999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f,
        //   999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f);
      } else if constexpr (eventtype == EMEventType::kEvent_Cent) {
        registry.fill(HIST("hCentFT0M"), collision.centFT0M());
        registry.fill(HIST("hCentFT0A"), collision.centFT0A());
        registry.fill(HIST("hCentFT0C"), collision.centFT0C());
        registry.fill(HIST("hCentNTPV"), collision.centNTPV());
        registry.fill(HIST("hCentNGlobal"), collision.centFT0M());

        event_cent(collision.centFT0M(), collision.centFT0A(), collision.centFT0C(), collision.centNTPV() /*, collision.centNGlobal()*/);
        event_qvec2(999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f);
        event_qvec3(999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f);
        // event_qvec(
        //   999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f,
        //   999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f);
      } else if constexpr (eventtype == EMEventType::kEvent_Cent_Qvec) {
        registry.fill(HIST("hCentFT0M"), collision.centFT0M());
        registry.fill(HIST("hCentFT0A"), collision.centFT0A());
        registry.fill(HIST("hCentFT0C"), collision.centFT0C());
        registry.fill(HIST("hCentNTPV"), collision.centNTPV());
        registry.fill(HIST("hCentNGlobal"), collision.centFT0M());

        event_cent(collision.centFT0M(), collision.centFT0A(), collision.centFT0C(), collision.centNTPV() /*, collision.centNGlobal()*/);
        float q2xft0m = 999.f, q2yft0m = 999.f, q2xft0a = 999.f, q2yft0a = 999.f, q2xft0c = 999.f, q2yft0c = 999.f, q2xfv0a = 999.f, q2yfv0a = 999.f, q2xbpos = 999.f, q2ybpos = 999.f, q2xbneg = 999.f, q2ybneg = 999.f, q2xbtot = 999.f, q2ybtot = 999.f;
        float q3xft0m = 999.f, q3yft0m = 999.f, q3xft0a = 999.f, q3yft0a = 999.f, q3xft0c = 999.f, q3yft0c = 999.f, q3xfv0a = 999.f, q3yfv0a = 999.f, q3xbpos = 999.f, q3ybpos = 999.f, q3xbneg = 999.f, q3ybneg = 999.f, q3xbtot = 999.f, q3ybtot = 999.f;

        if (collision.qvecFT0CReVec().size() >= 2) { // harmonics 2,3
          q2xft0m = collision.qvecFT0MReVec()[0], q2xft0a = collision.qvecFT0AReVec()[0], q2xft0c = collision.qvecFT0CReVec()[0], q2xfv0a = collision.qvecFV0AReVec()[0], q2xbpos = collision.qvecBPosReVec()[0], q2xbneg = collision.qvecBNegReVec()[0], q2xbtot = collision.qvecBTotReVec()[0];
          q2yft0m = collision.qvecFT0MImVec()[0], q2yft0a = collision.qvecFT0AImVec()[0], q2yft0c = collision.qvecFT0CImVec()[0], q2yfv0a = collision.qvecFV0AImVec()[0], q2ybpos = collision.qvecBPosImVec()[0], q2ybneg = collision.qvecBNegImVec()[0], q2ybtot = collision.qvecBTotImVec()[0];
          q3xft0m = collision.qvecFT0MReVec()[1], q3xft0a = collision.qvecFT0AReVec()[1], q3xft0c = collision.qvecFT0CReVec()[1], q3xfv0a = collision.qvecFV0AReVec()[1], q3xbpos = collision.qvecBPosReVec()[1], q3xbneg = collision.qvecBNegReVec()[1], q3xbtot = collision.qvecBTotReVec()[1];
          q3yft0m = collision.qvecFT0MImVec()[1], q3yft0a = collision.qvecFT0AImVec()[1], q3yft0c = collision.qvecFT0CImVec()[1], q3yfv0a = collision.qvecFV0AImVec()[1], q3ybpos = collision.qvecBPosImVec()[1], q3ybneg = collision.qvecBNegImVec()[1], q3ybtot = collision.qvecBTotImVec()[1];
        } else if (collision.qvecFT0CReVec().size() >= 1) { // harmonics 2
          q2xft0m = collision.qvecFT0MReVec()[0], q2xft0a = collision.qvecFT0AReVec()[0], q2xft0c = collision.qvecFT0CReVec()[0], q2xfv0a = collision.qvecFV0AReVec()[0], q2xbpos = collision.qvecBPosReVec()[0], q2xbneg = collision.qvecBNegReVec()[0], q2xbtot = collision.qvecBTotReVec()[0];
          q2yft0m = collision.qvecFT0MImVec()[0], q2yft0a = collision.qvecFT0AImVec()[0], q2yft0c = collision.qvecFT0CImVec()[0], q2yfv0a = collision.qvecFV0AImVec()[0], q2ybpos = collision.qvecBPosImVec()[0], q2ybneg = collision.qvecBNegImVec()[0], q2ybtot = collision.qvecBTotImVec()[0];
        }
        // event_qvec(q2xft0m, q2yft0m, q2xft0a, q2yft0a, q2xft0c, q2yft0c, q2xfv0a, q2yfv0a, q2xbpos, q2ybpos, q2xbneg, q2ybneg, q2xbtot, q2ybtot,
        //   q3xft0m, q3yft0m, q3xft0a, q3yft0a, q3xft0c, q3yft0c, q3xfv0a, q3yfv0a, q3xbpos, q3ybpos, q3xbneg, q3ybneg, q3xbtot, q3ybtot);
        event_qvec2(q2xft0m, q2yft0m, q2xft0a, q2yft0a, q2xft0c, q2yft0c, q2xfv0a, q2yfv0a, q2xbpos, q2ybpos, q2xbneg, q2ybneg, q2xbtot, q2ybtot);
        event_qvec3(q3xft0m, q3yft0m, q3xft0a, q3yft0a, q3xft0c, q3yft0c, q3xfv0a, q3yfv0a, q3xbpos, q3ybpos, q3xbneg, q3ybneg, q3xbtot, q3ybtot);
      } else {
        event_cent(105.f, 105.f, 105.f, 105.f);
        // event_qvec( 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f,
        //   999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f);
        event_qvec2(999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f);
        event_qvec3(999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f);
      }
    } // end of collision loop
  } // end of skimEvent

  //---------- for data ----------

  void processEvent(MyCollisions const& collisions, MyBCs const& bcs)
  {
    skimEvent<false, false, EMEventType::kEvent>(collisions, bcs);
  }
  PROCESS_SWITCH(CreateEMEventDilepton, processEvent, "process event info", false);

  void processEvent_Cent(MyCollisions_Cent const& collisions, MyBCs const& bcs)
  {
    skimEvent<false, false, EMEventType::kEvent_Cent>(collisions, bcs);
  }
  PROCESS_SWITCH(CreateEMEventDilepton, processEvent_Cent, "process event info", false);

  void processEvent_Cent_Qvec(MyCollisions_Cent_Qvec const& collisions, MyBCs const& bcs)
  {
    skimEvent<false, false, EMEventType::kEvent_Cent_Qvec>(collisions, bcs);
  }
  PROCESS_SWITCH(CreateEMEventDilepton, processEvent_Cent_Qvec, "process event info", false);

  //---------- for data with swt ----------

  void processEvent_SWT(MyCollisionsWithSWT const& collisions, MyBCs const& bcs)
  {
    skimEvent<false, true, EMEventType::kEvent>(collisions, bcs);
  }
  PROCESS_SWITCH(CreateEMEventDilepton, processEvent_SWT, "process event info", false);

  void processEvent_SWT_Cent(MyCollisionsWithSWT_Cent const& collisions, MyBCs const& bcs)
  {
    skimEvent<false, true, EMEventType::kEvent_Cent>(collisions, bcs);
  }
  PROCESS_SWITCH(CreateEMEventDilepton, processEvent_SWT_Cent, "process event info", false);

  void processEvent_SWT_Cent_Qvec(MyCollisionsWithSWT_Cent_Qvec const& collisions, MyBCs const& bcs)
  {
    skimEvent<false, true, EMEventType::kEvent_Cent_Qvec>(collisions, bcs);
  }
  PROCESS_SWITCH(CreateEMEventDilepton, processEvent_SWT_Cent_Qvec, "process event info", false);

  //---------- for MC ----------

  void processEventMC(MyCollisionsMC const& collisions, MyBCs const& bcs)
  {
    skimEvent<true, false, EMEventType::kEvent>(collisions, bcs);
  }
  PROCESS_SWITCH(CreateEMEventDilepton, processEventMC, "process event info", false);

  void processEventMC_Cent(MyCollisionsMC_Cent const& collisions, MyBCs const& bcs)
  {
    skimEvent<true, false, EMEventType::kEvent_Cent>(collisions, bcs);
  }
  PROCESS_SWITCH(CreateEMEventDilepton, processEventMC_Cent, "process event info", false);

  void processEventMC_Cent_Qvec(MyCollisionsMC_Cent_Qvec const& collisions, MyBCs const& bcs)
  {
    skimEvent<true, false, EMEventType::kEvent_Cent_Qvec>(collisions, bcs);
  }
  PROCESS_SWITCH(CreateEMEventDilepton, processEventMC_Cent_Qvec, "process event info", false);

  void processDummy(aod::Collisions const&) {}
  PROCESS_SWITCH(CreateEMEventDilepton, processDummy, "processDummy", true);
};
struct AssociateDileptonToEMEvent {
  Produces<o2::aod::EMPrimaryElectronEMEventIds> prmeleventid;
  Produces<o2::aod::EMPrimaryMuonEMEventIds> prmmueventid;
  Produces<o2::aod::EMPrimaryTrackEMEventIds> prmtrackeventid;

  PresliceUnsorted<aod::EMPrimaryElectrons> perCollision_el = aod::emprimaryelectron::collisionId;
  PresliceUnsorted<aod::EMPrimaryMuons> perCollision_mu = aod::emprimarymuon::collisionId;
  Preslice<aod::EMPrimaryTracks> perCollision_track = aod::emprimarytrack::collisionId;

  void init(o2::framework::InitContext&) {}

  template <typename TCollisions, typename TLeptons, typename TEventIds, typename TPreslice>
  void fillEventId(TCollisions const& collisions, TLeptons const& leptons, TEventIds& eventIds, TPreslice const& perCollision)
  {
    for (const auto& collision : collisions) {
      auto leptons_coll = leptons.sliceBy(perCollision, collision.collisionId());
      int nl = leptons_coll.size();
      // LOGF(info, "collision.collisionId() = %d , nl = %d", collision.collisionId(), nl);
      for (int il = 0; il < nl; il++) {
        eventIds(collision.globalIndex());
      } // end of photon loop
    } // end of collision loop
  }

  // This struct is for both data and MC.
  // Note that reconstructed collisions without mc collisions are already rejected in CreateEMEventDilepton in MC.

  void processElectron(aod::EMEvents const& collisions, aod::EMPrimaryElectrons const& tracks)
  {
    fillEventId(collisions, tracks, prmeleventid, perCollision_el);
  }

  void processFwdMuon(aod::EMEvents const& collisions, aod::EMPrimaryMuons const& tracks)
  {
    fillEventId(collisions, tracks, prmmueventid, perCollision_mu);
  }

  void processChargedTrack(aod::EMEvents const& collisions, aod::EMPrimaryTracks const& tracks)
  {
    fillEventId(collisions, tracks, prmtrackeventid, perCollision_track);
  }

  void processDummy(aod::EMEvents const&) {}

  PROCESS_SWITCH(AssociateDileptonToEMEvent, processElectron, "process indexing for electrons", false);
  PROCESS_SWITCH(AssociateDileptonToEMEvent, processFwdMuon, "process indexing for forward muons", false);
  PROCESS_SWITCH(AssociateDileptonToEMEvent, processChargedTrack, "process indexing for charged tracks", false);
  PROCESS_SWITCH(AssociateDileptonToEMEvent, processDummy, "process dummy", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CreateEMEventDilepton>(cfgc, TaskName{"create-emevent-dilepton"}),
    adaptAnalysisTask<AssociateDileptonToEMEvent>(cfgc, TaskName{"associate-dilepton-to-emevent"}),
  };
}
