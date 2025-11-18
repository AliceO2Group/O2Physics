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

/// \file createEMEventPhoton.cxx
/// \brief This code produces reduced events for photon analyses.
/// \author Daiki Sekihata, daiki.sekihata@cern.ch

#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
//
#include "PWGJE/DataModel/Jet.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Qvectors.h"

#include <CCDB/BasicCCDBManager.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DataFormatsParameters/GRPObject.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

using MyBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;
using MyQvectors = soa::Join<aod::QvectorFT0CVecs, aod::QvectorFT0AVecs, aod::QvectorFT0MVecs, aod::QvectorFV0AVecs, aod::QvectorBPosVecs, aod::QvectorBNegVecs, aod::QvectorBTotVecs>;

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::EMEvSels, aod::EMEoIs, aod::Mults>;
using MyCollisionsCent = soa::Join<MyCollisions, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>; // centrality table has dependency on multiplicity table.
using MyCollisionsCentQvec = soa::Join<MyCollisionsCent, MyQvectors>;

using MyCollisionsWithSWT = soa::Join<MyCollisions, aod::EMSWTriggerBitsTMP>;
using MyCollisionsWithSWT_Cent = soa::Join<MyCollisionsWithSWT, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>; // centrality table has dependency on multiplicity table.
using MyCollisionsWithSWT_Cent_Qvec = soa::Join<MyCollisionsWithSWT_Cent, MyQvectors>;

using MyCollisionsMC = soa::Join<MyCollisions, aod::McCollisionLabels>;
using MyCollisionsMCCent = soa::Join<MyCollisionsMC, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>; // centrality table has dependency on multiplicity table.
using MyCollisionsMCCentQvec = soa::Join<MyCollisionsMCCent, MyQvectors>;

struct CreateEMEventPhoton {
  Produces<o2::aod::EMBCs> embc;
  Produces<o2::aod::EMEvents> event;
  // Produces<o2::aod::EMEventsCov> eventCov;
  Produces<o2::aod::EMEventsMult> eventMult;
  Produces<o2::aod::EMEventsCent> eventCent;
  Produces<o2::aod::EMEventsQvec> eventQvec;
  Produces<o2::aod::EMSWTriggerBits> emswtbit;
  Produces<o2::aod::EMEventNormInfos> event_norm_info;
  Produces<o2::aod::EMEventsWeight> eventWeights;

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
  Configurable<double> dBzInput{"d_bz", -999, "bz field, -999 is automatic"};
  Configurable<bool> needEMCTrigger{"needEMCTrigger", false, "flag to only save events which have kTVXinEMC trigger bit. To reduce PbPb derived data size"};
  Configurable<bool> needPHSTrigger{"needPHSTrigger", false, "flag to only save events which have kTVXinPHOS trigger bit. To reduce PbPb derived data size"};
  Configurable<bool> enableJJHistograms{"enableJJHistograms", false, "flag to fill JJ QA histograms for outlier rejection"};
  Configurable<float> maxpTJetOverpTHard{"maxpTJetOverpTHard", 2., "set weight to 0 for JJ events with larger pTJet/pTHard"};

  HistogramRegistry registry{"registry"};
  void init(o2::framework::InitContext&)
  {
    auto hEventCounter = registry.add<TH1>("hEventCounter", "hEventCounter", kTH1I, {{7, 0.5f, 7.5f}});
    hEventCounter->GetXaxis()->SetBinLabel(1, "all");
    hEventCounter->GetXaxis()->SetBinLabel(2, "sel8");

    if (enableJJHistograms) {
      auto hJJ_pTHardVsJetpT = registry.add<TH2>("hJJ_pTHardVsJetpT", "hJJ_pTHardVsJetpT;#bf{#it{p}_{T}^{hard}};#bf{#it{p}_{T}^{leading jet}}", kTH2F, {{500, 0, 1000}, {500, 0, 1000}});
    }
  }

  int mRunNumber;
  float dBz;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  template <typename TBC>
  void initCCDB(TBC const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    // In case override, don't proceed, please - no CCDB access required
    if (dBzInput > -990) {
      dBz = dBzInput;
      o2::parameters::GRPMagField grpmag;
      if (std::fabs(dBz) > 1e-5) {
        grpmag.setL3Current(30000.f / (dBz / 5.0f));
      }
      mRunNumber = bc.runNumber();
      return;
    }

    auto run3GRPTimestamp = bc.timestamp();
    o2::parameters::GRPObject* grpo = 0x0;
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (!skipGRPOquery)
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3GRPTimestamp);
    if (grpo) {
      // Fetch magnetic field from ccdb for current collision
      dBz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3GRPTimestamp << " with magnetic field of " << dBz << " kZG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3GRPTimestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3GRPTimestamp;
      }
      // Fetch magnetic field from ccdb for current collision
      dBz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for timestamp " << run3GRPTimestamp << " with magnetic field of " << dBz << " kZG";
    }
    mRunNumber = bc.runNumber();
  }

  template <bool isMC, bool isTriggerAnalysis, EMEventType eventtype, typename TCollisions, typename TBCs>
  void skimEvent(TCollisions const& collisions, TBCs const& bcs)
  {
    for (const auto& bc : bcs) {
      if (bc.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
        embc(bc.alias_raw(), bc.selection_raw(), bc.rct_raw()); // TVX is fired.
      }
    } // end of bc loop

    for (const auto& collision : collisions) {
      if constexpr (isMC) {
        if (!collision.has_mcCollision()) {
          continue;
        }
      }

      auto bc = collision.template foundBC_as<TBCs>();
      initCCDB(bc);

      if (needEMCTrigger && !collision.alias_bit(kTVXinEMC)) {
        continue;
      }
      if (needPHSTrigger && !collision.alias_bit(kTVXinPHOS)) {
        continue;
      }

      if (collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
        if constexpr (eventtype == EMEventType::kEvent) {
          event_norm_info(collision.alias_raw(), collision.selection_raw(), collision.rct_raw(), static_cast<int16_t>(10.f * collision.posZ()), 105.f);
        } else if constexpr (eventtype == EMEventType::kEvent_Cent || eventtype == EMEventType::kEvent_Cent_Qvec) {
          event_norm_info(collision.alias_raw(), collision.selection_raw(), collision.rct_raw(), static_cast<int16_t>(10.f * collision.posZ()), collision.centFT0C());
        } else {
          event_norm_info(collision.alias_raw(), collision.selection_raw(), collision.rct_raw(), static_cast<int16_t>(10.f * collision.posZ()), 105.f);
        }
      }

      if (!collision.isSelected()) {
        continue;
      }

      if (!collision.isEoI()) { // events with at least 1 photon for data reduction.
        continue;
      }

      if constexpr (isTriggerAnalysis) {
        if (collision.swtaliastmp_raw() == 0) {
          continue;
        } else {
          emswtbit(collision.swtaliastmp_raw());
        }
      }

      const float qDefault = 999.f; // default value for q vectors if not obtained

      registry.fill(HIST("hEventCounter"), 1);

      if (collision.sel8()) {
        registry.fill(HIST("hEventCounter"), 2);
      }

      event(collision.globalIndex(), bc.runNumber(), bc.globalBC(), collision.alias_raw(), collision.selection_raw(), collision.rct_raw(), bc.timestamp(),
            collision.posZ(),
            collision.numContrib(), collision.trackOccupancyInTimeRange(), collision.ft0cOccupancyInTimeRange());

      // eventCov(collision.covXX(), collision.covXY(), collision.covXZ(), collision.covYY(), collision.covYZ(), collision.covZZ(), collision.chi2());

      eventMult(collision.multFT0A(), collision.multFT0C(), collision.multNTracksPV(), collision.multNTracksPVeta1(), collision.multNTracksPVetaHalf());

      if constexpr (eventtype != EMEventType::kEvent_JJ) {
        eventWeights(1.f);
      }

      if constexpr (eventtype == EMEventType::kEvent) {
        eventCent(105.f, 105.f, 105.f);
        eventQvec(qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault,
                  qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault);
      } else if constexpr (eventtype == EMEventType::kEvent_Cent) {
        eventCent(collision.centFT0M(), collision.centFT0A(), collision.centFT0C());
        eventQvec(qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault,
                  qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault);
      } else if constexpr (eventtype == EMEventType::kEvent_Cent_Qvec) {
        eventCent(collision.centFT0M(), collision.centFT0A(), collision.centFT0C());
        const size_t qvecSize = collision.qvecFT0CReVec().size();
        if (qvecSize >= 2) { // harmonics 2,3
          eventQvec(collision.qvecFT0MReVec()[0], collision.qvecFT0MImVec()[0], collision.qvecFT0AReVec()[0], collision.qvecFT0AImVec()[0], collision.qvecFT0CReVec()[0], collision.qvecFT0CImVec()[0], collision.qvecFV0AReVec()[0], collision.qvecFV0AImVec()[0], collision.qvecBPosReVec()[0], collision.qvecBPosImVec()[0], collision.qvecBNegReVec()[0], collision.qvecBNegImVec()[0], collision.qvecBTotReVec()[0], collision.qvecBTotImVec()[0],
                    collision.qvecFT0MReVec()[1], collision.qvecFT0MImVec()[1], collision.qvecFT0AReVec()[1], collision.qvecFT0AImVec()[1], collision.qvecFT0CReVec()[1], collision.qvecFT0CImVec()[1], collision.qvecFV0AReVec()[1], collision.qvecFV0AImVec()[1], collision.qvecBPosReVec()[1], collision.qvecBPosImVec()[1], collision.qvecBNegReVec()[1], collision.qvecBNegImVec()[1], collision.qvecBTotReVec()[1], collision.qvecBTotImVec()[1]);
        } else if (qvecSize >= 1) { // harmonics 2
          eventQvec(collision.qvecFT0MReVec()[0], collision.qvecFT0MImVec()[0], collision.qvecFT0AReVec()[0], collision.qvecFT0AImVec()[0], collision.qvecFT0CReVec()[0], collision.qvecFT0CImVec()[0], collision.qvecFV0AReVec()[0], collision.qvecFV0AImVec()[0], collision.qvecBPosReVec()[0], collision.qvecBPosImVec()[0], collision.qvecBNegReVec()[0], collision.qvecBNegImVec()[0], collision.qvecBTotReVec()[0], collision.qvecBTotImVec()[0],
                    qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault);
        }
      } else {
        eventCent(105.f, 105.f, 105.f);
        eventQvec(qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault,
                  qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault, qDefault);
      }
    } // end of collision loop

  } // end of skimEvent

  Preslice<aod::FullMCParticleLevelJets> perCollision_jet = aod::jet::mcCollisionId;

  using MyJJCollisions = soa::Join<aod::McCollisions, aod::HepMCXSections>;

  void fillEventWeights(MyCollisionsMC const& collisions, MyJJCollisions const&, MyBCs const&, aod::FullMCParticleLevelJets const& jets)
  {
    for (const auto& collision : collisions) {
      if (!collision.has_mcCollision()) {
        continue;
      }
      if (!collision.isSelected()) {
        continue;
      }

      auto bc = collision.template foundBC_as<MyBCs>();
      initCCDB(bc);

      auto mcCollision = collision.mcCollision_as<MyJJCollisions>();

      // Outlier rejection: Set weight to 0 for events with large pTJet/pTHard
      // ----------------------------------------------------------------------
      auto jetsInThisCollision = jets.sliceBy(perCollision_jet, mcCollision.globalIndex());
      float collisionWeight = mcCollision.weight();
      for (const auto& jet : jetsInThisCollision) {
        if (jet.pt() > maxpTJetOverpTHard * mcCollision.ptHard())
          collisionWeight = 0.f;
        registry.fill(HIST("hJJ_pTHardVsJetpT"), mcCollision.ptHard(), jet.pt());
      }

      eventWeights(collisionWeight);
    }
  }

  // for data
  void processEvent(MyCollisions const& collisions, MyBCs const& bcs)
  {
    skimEvent<false, false, EMEventType::kEvent>(collisions, bcs);
  }
  PROCESS_SWITCH(CreateEMEventPhoton, processEvent, "process event info", false);

  void processEvent_Cent(MyCollisionsCent const& collisions, MyBCs const& bcs)
  {
    skimEvent<false, false, EMEventType::kEvent_Cent>(collisions, bcs);
  }
  PROCESS_SWITCH(CreateEMEventPhoton, processEvent_Cent, "process event info", false);

  void processEvent_Cent_Qvec(MyCollisionsCentQvec const& collisions, MyBCs const& bcs)
  {
    skimEvent<false, false, EMEventType::kEvent_Cent_Qvec>(collisions, bcs);
  }
  PROCESS_SWITCH(CreateEMEventPhoton, processEvent_Cent_Qvec, "process event info", false);

  void processEvent_SWT(MyCollisionsWithSWT const& collisions, MyBCs const& bcs)
  {
    skimEvent<false, true, EMEventType::kEvent>(collisions, bcs);
  }
  PROCESS_SWITCH(CreateEMEventPhoton, processEvent_SWT, "process event info", false);

  void processEvent_SWT_Cent(MyCollisionsWithSWT_Cent const& collisions, MyBCs const& bcs)
  {
    skimEvent<false, true, EMEventType::kEvent_Cent>(collisions, bcs);
  }
  PROCESS_SWITCH(CreateEMEventPhoton, processEvent_SWT_Cent, "process event info", false);

  void processEvent_SWT_Cent_Qvec(MyCollisionsWithSWT_Cent_Qvec const& collisions, MyBCs const& bcs)
  {
    skimEvent<false, true, EMEventType::kEvent_Cent_Qvec>(collisions, bcs);
  }
  PROCESS_SWITCH(CreateEMEventPhoton, processEvent_SWT_Cent_Qvec, "process event info", false);

  // for MC
  void processEventMC(MyCollisionsMC const& collisions, MyBCs const& bcs)
  {
    skimEvent<true, false, EMEventType::kEvent>(collisions, bcs);
  }
  PROCESS_SWITCH(CreateEMEventPhoton, processEventMC, "process event info", false);

  void processEventJJMC(MyCollisionsMC const& collisions, MyJJCollisions const& mcCollisions, MyBCs const& bcs, aod::FullMCParticleLevelJets const& jets)
  {
    skimEvent<true, false, EMEventType::kEvent_JJ>(collisions, bcs);
    fillEventWeights(collisions, mcCollisions, bcs, jets);
  }
  PROCESS_SWITCH(CreateEMEventPhoton, processEventJJMC, "process event info", false);

  void processEventMC_Cent(MyCollisionsMCCent const& collisions, MyBCs const& bcs)
  {
    skimEvent<true, false, EMEventType::kEvent_Cent>(collisions, bcs);
  }
  PROCESS_SWITCH(CreateEMEventPhoton, processEventMC_Cent, "process event info", false);

  void processEventMC_Cent_Qvec(MyCollisionsMCCentQvec const& collisions, MyBCs const& bcs)
  {
    skimEvent<true, false, EMEventType::kEvent_Cent_Qvec>(collisions, bcs);
  }
  PROCESS_SWITCH(CreateEMEventPhoton, processEventMC_Cent_Qvec, "process event info", false);

  void processDummy(aod::Collisions const&) {}
  PROCESS_SWITCH(CreateEMEventPhoton, processDummy, "processDummy", true);
};
struct AssociatePhotonToEMEvent {
  Produces<o2::aod::V0KFEMEventIds> v0kfeventid;
  Produces<o2::aod::EMPrimaryElectronEMEventIds> prmeleventid;
  Produces<o2::aod::PHOSEMEventIds> phoseventid;
  Produces<o2::aod::EMCEMEventIds> emceventid;
  Produces<o2::aod::EMPrimaryTrackEMEventIds> prmtrackeventid;

  Preslice<aod::V0PhotonsKF> perCollisionPCM = aod::v0photonkf::collisionId;
  PresliceUnsorted<aod::EMPrimaryElectronsFromDalitz> perCollisionEl = aod::emprimaryelectron::collisionId;
  Preslice<aod::PHOSClusters> perCollisionPHOS = aod::skimmedcluster::collisionId;
  Preslice<aod::SkimEMCClusters> perCollisionEMC = aod::skimmedcluster::collisionId;
  Preslice<aod::EMPrimaryTracks> perCollision_track = aod::emprimarytrack::collisionId;

  void init(o2::framework::InitContext&) {}

  template <typename TCollisions, typename TPhotons, typename TEventIds, typename TPreslice>
  void fillEventId(TCollisions const& collisions, TPhotons const& photons, TEventIds& eventIds, TPreslice const& perCollision)
  {
    for (const auto& collision : collisions) {
      auto photonsColl = photons.sliceBy(perCollision, collision.collisionId());
      int ng = photonsColl.size();
      for (int ig = 0; ig < ng; ig++) {
        eventIds(collision.globalIndex());
      } // end of photon loop
    } // end of collision loop
  }

  // This struct is for both data and MC.
  // Note that reconstructed collisions without mc collisions are already rejected in CreateEMEventPhoton in MC.

  void processPCM(aod::EMEvents const& collisions, aod::V0PhotonsKF const& photons)
  {
    fillEventId(collisions, photons, v0kfeventid, perCollisionPCM);
  }

  void processElectronFromDalitz(aod::EMEvents const& collisions, aod::EMPrimaryElectronsFromDalitz const& tracks)
  {
    fillEventId(collisions, tracks, prmeleventid, perCollisionEl);
  }

  void processPHOS(aod::EMEvents const& collisions, aod::PHOSClusters const& photons)
  {
    fillEventId(collisions, photons, phoseventid, perCollisionPHOS);
  }

  void processEMC(aod::EMEvents const& collisions, aod::SkimEMCClusters const& photons)
  {
    fillEventId(collisions, photons, emceventid, perCollisionEMC);
  }

  void processChargedTrack(aod::EMEvents const& collisions, aod::EMPrimaryTracks const& tracks)
  {
    fillEventId(collisions, tracks, prmtrackeventid, perCollision_track);
  }

  void processDummy(aod::EMEvents const&) {}

  PROCESS_SWITCH(AssociatePhotonToEMEvent, processPCM, "process pcm-event indexing", false);
  PROCESS_SWITCH(AssociatePhotonToEMEvent, processElectronFromDalitz, "process dalitzee-event indexing", false);
  PROCESS_SWITCH(AssociatePhotonToEMEvent, processPHOS, "process phos-event indexing", false);
  PROCESS_SWITCH(AssociatePhotonToEMEvent, processEMC, "process emc-event indexing", false);
  PROCESS_SWITCH(AssociatePhotonToEMEvent, processChargedTrack, "process indexing for charged tracks", false);
  PROCESS_SWITCH(AssociatePhotonToEMEvent, processDummy, "process dummy", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CreateEMEventPhoton>(cfgc, TaskName{"create-emevent-photon"}),
    adaptAnalysisTask<AssociatePhotonToEMEvent>(cfgc, TaskName{"associate-photon-to-emevent"}),
  };
}
