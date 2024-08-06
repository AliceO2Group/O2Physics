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
#include "EventFiltering/Zorro.h"
#include "Common/Core/TableHelper.h"

#include "PWGEM/Dilepton/DataModel/dileptonTables.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

using MyBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;
using MyQvectors = soa::Join<aod::QvectorFT0CVecs, aod::QvectorFT0AVecs, aod::QvectorFT0MVecs, aod::QvectorBPosVecs, aod::QvectorBNegVecs, aod::QvectorBTotVecs>;

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>;
using MyCollisions_Cent = soa::Join<MyCollisions, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentNTPVs>; // centrality table has dependency on multiplicity table.
using MyCollisions_Cent_Qvec = soa::Join<MyCollisions_Cent, MyQvectors>;

using MyCollisionsMC = soa::Join<MyCollisions, aod::McCollisionLabels>;
using MyCollisionsMC_Cent = soa::Join<MyCollisionsMC, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentNTPVs>; // centrality table has dependency on multiplicity table.
using MyCollisionsMC_Cent_Qvec = soa::Join<MyCollisionsMC_Cent, MyQvectors>;

struct CreateEMEventDilepton {
  Produces<o2::aod::EMEvents> event;
  // Produces<o2::aod::EMEventsCov> eventcov;
  Produces<o2::aod::EMEventsMult> event_mult;
  Produces<o2::aod::EMEventsCent> event_cent;
  Produces<o2::aod::EMEventsQvec> event_qvec;
  Produces<o2::aod::EMSWTriggerInfos> emswtbit;

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
  Configurable<bool> applyEveSel_at_skimming{"applyEveSel_at_skimming", false, "flag to apply minimal event selection at the skimming level"};
  Configurable<bool> enable_swt{"enable_swt", false, "flag to process skimmed data (swt triggered)"};
  Configurable<std::string> cfg_swt_names{"cfg_swt_names", "fHighTrackMult,fHighFt0Mult", "comma-separated software trigger names"}; // !trigger names have to be pre-registered in dileptonTable.h for bit operation!

  HistogramRegistry registry{"registry"};
  void init(o2::framework::InitContext&)
  {
    auto hEventCounter = registry.add<TH1>("hEventCounter", "hEventCounter", kTH1I, {{7, 0.5f, 7.5f}});
    hEventCounter->GetXaxis()->SetBinLabel(1, "all");
    hEventCounter->GetXaxis()->SetBinLabel(2, "sel8");

    registry.add("hNInspectedTVX", "N inspected TVX;run number;N_{TVX}", kTProfile, {{80000, 520000.5, 600000.5}}, true);
  }

  ~CreateEMEventDilepton()
  {
    swt_names.clear();
    swt_names.shrink_to_fit();
  }

  Zorro zorro;
  std::vector<int> mTOIidx;
  std::vector<std::string> swt_names;
  uint64_t mNinspectedTVX{0};

  int mRunNumber;
  float d_bz;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  template <typename TBC>
  void initCCDB(TBC const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    if (enable_swt) {
      LOGF(info, "enable software triggers : %s", cfg_swt_names.value.data());
      mTOIidx = zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), cfg_swt_names.value);
      std::stringstream tokenizer(cfg_swt_names.value);
      std::string token;
      while (std::getline(tokenizer, token, ',')) {
        swt_names.emplace_back(token);
      }
      for (auto& idx : mTOIidx) {
        LOGF(info, "Trigger of Interest : index = %d", idx);
      }
      mNinspectedTVX = zorro.getInspectedTVX()->GetBinContent(1);
      LOGF(info, "total inspected TVX events = %d in run number %d", mNinspectedTVX, bc.runNumber());
      registry.fill(HIST("hNInspectedTVX"), bc.runNumber(), mNinspectedTVX);
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

  Preslice<aod::V0PhotonsKF> perCollision_pcm = aod::v0photonkf::collisionId;
  PresliceUnsorted<aod::EMPrimaryElectrons> perCollision_el = aod::emprimaryelectron::collisionId;
  PresliceUnsorted<aod::EMPrimaryMuons> perCollision_mu = aod::emprimarymuon::collisionId;

  template <bool isMC, EMEventType eventype, typename TCollisions, typename TBCs>
  void skimEvent(TCollisions const& collisions, TBCs const&)
  {
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

      if (enable_swt) {
        if (zorro.isSelected(bc.globalBC())) {     // triggered event
          auto swt_bitset = zorro.getLastResult(); // this has to be called after zorro::isSelected, or simply call zorro.fetch
          // LOGF(info, "swt_bitset.to_string().c_str() = %s", swt_bitset.to_string().c_str());
          uint16_t trigger_bitmap = 0;
          for (size_t idx = 0; idx < mTOIidx.size(); idx++) {
            if (swt_bitset.test(mTOIidx[idx])) {
              auto swtname = swt_names[idx];
              trigger_bitmap |= BIT(o2::aod::pwgem::dilepton::swt::aliasLabels.at(swtname));
              // LOGF(info, "swtname = %s is fired. swt index in original swt table = %d, swt index for EM table = %d", swtname.data(), mTOIidx[idx], o2::aod::pwgem::dilepton::swt::aliasLabels.at(swtname));
            }
          }
          emswtbit(trigger_bitmap, mNinspectedTVX);
        } else { // rejected
          continue;
        }
      }

      // LOGF(info, "collision.multNTracksPV() = %d, collision.multFT0A() = %f, collision.multFT0C() = %f", collision.multNTracksPV(), collision.multFT0A(), collision.multFT0C());

      registry.fill(HIST("hEventCounter"), 1);

      event(collision.globalIndex(), bc.runNumber(), bc.globalBC(), collision.alias_raw(), collision.selection_raw(), bc.timestamp(),
            collision.posX(), collision.posY(), collision.posZ(),
            collision.numContrib(), collision.trackOccupancyInTimeRange());

      // eventcov(collision.covXX(), collision.covXY(), collision.covXZ(), collision.covYY(), collision.covYZ(), collision.covZZ(), collision.chi2());

      event_mult(collision.multFT0A(), collision.multFT0C(), collision.multTPC(), collision.multNTracksPV(), collision.multNTracksPVeta1(), collision.multNTracksPVetaHalf());

      if constexpr (eventype == EMEventType::kEvent) {
        event_cent(105.f, 105.f, 105.f, 105.f);
        event_qvec(
          999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f,
          999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f,
          999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f);
      } else if constexpr (eventype == EMEventType::kEvent_Cent) {
        event_cent(collision.centFT0M(), collision.centFT0A(), collision.centFT0C(), collision.centNTPV());
        event_qvec(
          999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f,
          999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f,
          999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f);
      } else if constexpr (eventype == EMEventType::kEvent_Cent_Qvec) {
        event_cent(collision.centFT0M(), collision.centFT0A(), collision.centFT0C(), collision.centNTPV());
        float q2xft0m = 999.f, q2yft0m = 999.f, q2xft0a = 999.f, q2yft0a = 999.f, q2xft0c = 999.f, q2yft0c = 999.f, q2xbpos = 999.f, q2ybpos = 999.f, q2xbneg = 999.f, q2ybneg = 999.f, q2xbtot = 999.f, q2ybtot = 999.f;
        float q3xft0m = 999.f, q3yft0m = 999.f, q3xft0a = 999.f, q3yft0a = 999.f, q3xft0c = 999.f, q3yft0c = 999.f, q3xbpos = 999.f, q3ybpos = 999.f, q3xbneg = 999.f, q3ybneg = 999.f, q3xbtot = 999.f, q3ybtot = 999.f;
        float q4xft0m = 999.f, q4yft0m = 999.f, q4xft0a = 999.f, q4yft0a = 999.f, q4xft0c = 999.f, q4yft0c = 999.f, q4xbpos = 999.f, q4ybpos = 999.f, q4xbneg = 999.f, q4ybneg = 999.f, q4xbtot = 999.f, q4ybtot = 999.f;

        if (collision.qvecFT0CReVec().size() >= 3) { // harmonics 2,3,4
          q2xft0m = collision.qvecFT0MReVec()[0], q2xft0a = collision.qvecFT0AReVec()[0], q2xft0c = collision.qvecFT0CReVec()[0], q2xbpos = collision.qvecBPosReVec()[0], q2xbneg = collision.qvecBNegReVec()[0], q2xbtot = collision.qvecBTotReVec()[0];
          q2yft0m = collision.qvecFT0MImVec()[0], q2yft0a = collision.qvecFT0AImVec()[0], q2yft0c = collision.qvecFT0CImVec()[0], q2ybpos = collision.qvecBPosImVec()[0], q2ybneg = collision.qvecBNegImVec()[0], q2ybtot = collision.qvecBTotImVec()[0];
          q3xft0m = collision.qvecFT0MReVec()[1], q3xft0a = collision.qvecFT0AReVec()[1], q3xft0c = collision.qvecFT0CReVec()[1], q3xbpos = collision.qvecBPosReVec()[1], q3xbneg = collision.qvecBNegReVec()[1], q3xbtot = collision.qvecBTotReVec()[1];
          q3yft0m = collision.qvecFT0MImVec()[1], q3yft0a = collision.qvecFT0AImVec()[1], q3yft0c = collision.qvecFT0CImVec()[1], q3ybpos = collision.qvecBPosImVec()[1], q3ybneg = collision.qvecBNegImVec()[1], q3ybtot = collision.qvecBTotImVec()[1];
          q4xft0m = collision.qvecFT0MReVec()[2], q4xft0a = collision.qvecFT0AReVec()[2], q4xft0c = collision.qvecFT0CReVec()[2], q4xbpos = collision.qvecBPosReVec()[2], q4xbneg = collision.qvecBNegReVec()[2], q4xbtot = collision.qvecBTotReVec()[2];
          q4yft0m = collision.qvecFT0MImVec()[2], q4yft0a = collision.qvecFT0AImVec()[2], q4yft0c = collision.qvecFT0CImVec()[2], q4ybpos = collision.qvecBPosImVec()[2], q4ybneg = collision.qvecBNegImVec()[2], q4ybtot = collision.qvecBTotImVec()[2];
        } else if (collision.qvecFT0CReVec().size() >= 2) { // harmonics 2,3
          q2xft0m = collision.qvecFT0MReVec()[0], q2xft0a = collision.qvecFT0AReVec()[0], q2xft0c = collision.qvecFT0CReVec()[0], q2xbpos = collision.qvecBPosReVec()[0], q2xbneg = collision.qvecBNegReVec()[0], q2xbtot = collision.qvecBTotReVec()[0];
          q2yft0m = collision.qvecFT0MImVec()[0], q2yft0a = collision.qvecFT0AImVec()[0], q2yft0c = collision.qvecFT0CImVec()[0], q2ybpos = collision.qvecBPosImVec()[0], q2ybneg = collision.qvecBNegImVec()[0], q2ybtot = collision.qvecBTotImVec()[0];
          q3xft0m = collision.qvecFT0MReVec()[1], q3xft0a = collision.qvecFT0AReVec()[1], q3xft0c = collision.qvecFT0CReVec()[1], q3xbpos = collision.qvecBPosReVec()[1], q3xbneg = collision.qvecBNegReVec()[1], q3xbtot = collision.qvecBTotReVec()[1];
          q3yft0m = collision.qvecFT0MImVec()[1], q3yft0a = collision.qvecFT0AImVec()[1], q3yft0c = collision.qvecFT0CImVec()[1], q3ybpos = collision.qvecBPosImVec()[1], q3ybneg = collision.qvecBNegImVec()[1], q3ybtot = collision.qvecBTotImVec()[1];
        } else if (collision.qvecFT0CReVec().size() >= 1) { // harmonics 2
          q2xft0m = collision.qvecFT0MReVec()[0], q2xft0a = collision.qvecFT0AReVec()[0], q2xft0c = collision.qvecFT0CReVec()[0], q2xbpos = collision.qvecBPosReVec()[0], q2xbneg = collision.qvecBNegReVec()[0], q2xbtot = collision.qvecBTotReVec()[0];
          q2yft0m = collision.qvecFT0MImVec()[0], q2yft0a = collision.qvecFT0AImVec()[0], q2yft0c = collision.qvecFT0CImVec()[0], q2ybpos = collision.qvecBPosImVec()[0], q2ybneg = collision.qvecBNegImVec()[0], q2ybtot = collision.qvecBTotImVec()[0];
        }
        event_qvec(
          q2xft0m, q2yft0m, q2xft0a, q2yft0a, q2xft0c, q2yft0c, q2xbpos, q2ybpos, q2xbneg, q2ybneg, q2xbtot, q2ybtot,
          q3xft0m, q3yft0m, q3xft0a, q3yft0a, q3xft0c, q3yft0c, q3xbpos, q3ybpos, q3xbneg, q3ybneg, q3xbtot, q3ybtot,
          q4xft0m, q4yft0m, q4xft0a, q4yft0a, q4xft0c, q4yft0c, q4xbpos, q4ybpos, q4xbneg, q4ybneg, q4xbtot, q4ybtot);
      } else {
        event_cent(105.f, 105.f, 105.f, 105.f);
        event_qvec(
          999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f,
          999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f,
          999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f);
      }
    } // end of collision loop
  }   // end of skimEvent

  void processEvent(MyCollisions const& collisions, MyBCs const& bcs)
  {
    skimEvent<false, EMEventType::kEvent>(collisions, bcs);
  }
  PROCESS_SWITCH(CreateEMEventDilepton, processEvent, "process event info", false);

  void processEvent_Cent(MyCollisions_Cent const& collisions, MyBCs const& bcs)
  {
    skimEvent<false, EMEventType::kEvent_Cent>(collisions, bcs);
  }
  PROCESS_SWITCH(CreateEMEventDilepton, processEvent_Cent, "process event info", false);

  void processEvent_Cent_Qvec(MyCollisions_Cent_Qvec const& collisions, MyBCs const& bcs)
  {
    skimEvent<false, EMEventType::kEvent_Cent_Qvec>(collisions, bcs);
  }
  PROCESS_SWITCH(CreateEMEventDilepton, processEvent_Cent_Qvec, "process event info", false);

  void processEventMC(MyCollisionsMC const& collisions, MyBCs const& bcs)
  {
    skimEvent<true, EMEventType::kEvent>(collisions, bcs);
  }
  PROCESS_SWITCH(CreateEMEventDilepton, processEventMC, "process event info", false);

  void processEventMC_Cent(MyCollisionsMC_Cent const& collisions, MyBCs const& bcs)
  {
    skimEvent<true, EMEventType::kEvent_Cent>(collisions, bcs);
  }
  PROCESS_SWITCH(CreateEMEventDilepton, processEventMC_Cent, "process event info", false);

  void processEventMC_Cent_Qvec(MyCollisionsMC_Cent_Qvec const& collisions, MyBCs const& bcs)
  {
    skimEvent<true, EMEventType::kEvent_Cent_Qvec>(collisions, bcs);
  }
  PROCESS_SWITCH(CreateEMEventDilepton, processEventMC_Cent_Qvec, "process event info", false);

  void processDummy(aod::Collisions const&) {}
  PROCESS_SWITCH(CreateEMEventDilepton, processDummy, "processDummy", true);
};
struct AssociateDileptonToEMEvent {
  Produces<o2::aod::V0KFEMEventIds> v0kfeventid;
  Produces<o2::aod::EMPrimaryElectronEMEventIds> prmeleventid;
  Produces<o2::aod::EMPrimaryMuonEMEventIds> prmmueventid;

  Preslice<aod::V0PhotonsKF> perCollision_pcm = aod::v0photonkf::collisionId;
  PresliceUnsorted<aod::EMPrimaryElectrons> perCollision_el = aod::emprimaryelectron::collisionId;
  PresliceUnsorted<aod::EMPrimaryMuons> perCollision_mu = aod::emprimarymuon::collisionId;

  void init(o2::framework::InitContext&) {}

  template <typename TCollisions, typename TLeptons, typename TEventIds, typename TPreslice>
  void fillEventId(TCollisions const& collisions, TLeptons const& leptons, TEventIds& eventIds, TPreslice const& perCollision)
  {
    for (auto& collision : collisions) {
      auto leptons_coll = leptons.sliceBy(perCollision, collision.collisionId());
      int nl = leptons_coll.size();
      // LOGF(info, "collision.collisionId() = %d , nl = %d", collision.collisionId(), nl);
      for (int il = 0; il < nl; il++) {
        eventIds(collision.globalIndex());
      } // end of photon loop
    }   // end of collision loop
  }

  // This struct is for both data and MC.
  // Note that reconstructed collisions without mc collisions are already rejected in CreateEMEventDilepton in MC.

  void processPCM(aod::EMEvents const& collisions, aod::V0PhotonsKF const& photons)
  {
    fillEventId(collisions, photons, v0kfeventid, perCollision_pcm);
  }

  void processElectron(aod::EMEvents const& collisions, aod::EMPrimaryElectrons const& tracks)
  {
    fillEventId(collisions, tracks, prmeleventid, perCollision_el);
  }

  void processFwdMuon(aod::EMEvents const& collisions, aod::EMPrimaryMuons const& tracks)
  {
    fillEventId(collisions, tracks, prmmueventid, perCollision_mu);
  }

  void processDummy(aod::EMEvents const&) {}

  PROCESS_SWITCH(AssociateDileptonToEMEvent, processPCM, "process pcm-event indexing", false);
  PROCESS_SWITCH(AssociateDileptonToEMEvent, processElectron, "process dalitzee-event indexing", false);
  PROCESS_SWITCH(AssociateDileptonToEMEvent, processFwdMuon, "process forward muon indexing", false);
  PROCESS_SWITCH(AssociateDileptonToEMEvent, processDummy, "process dummy", true);
};
struct EMEventPropertyTask {
  SliceCache cache;
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;
  using Run3Tracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA>;
  Zorro zorro;

  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<bool> applyEveSel_at_skimming{"applyEveSel_at_skimming", false, "flag to apply minimal event selection at the skimming level"};
  Configurable<bool> enable_swt{"enable_swt", false, "flag to process skimmed data (swt triggered)"};
  Configurable<std::string> cfg_swt_names{"cfg_swt_names", "", "comma-separated software trigger names"};
  Configurable<bool> fillQAHistogram{"fillQAHistogram", false, "flag to fill QA histograms"};

  Produces<aod::EMEventsProperty> evprop;
  struct : ConfigurableGroup {
    std::string prefix = "spherocity_cutgroup";
    Configurable<bool> require_isPVContributor{"require_isPVContributor", false, "require tracks to be PV contributors"};
    Configurable<int> min_ntrack{"min_ntrack", 3, "min. number of tracks"};
    Configurable<float> min_pt{"min_pt", 0.15, "min. pT of track in GeV/c"};
    Configurable<float> min_eta{"min_eta", -0.8, "min. eta of track"};
    Configurable<float> max_eta{"max_eta", +0.8, "max. eta of track"};
    Configurable<float> max_dcaxy{"max_dcaxy", 2.4, "max. DCAxy of track in cm"};
    Configurable<float> max_dcaz{"max_dcaz", 3.2, "max. DCAz of track in cm"};
    Configurable<int> min_ncluster_tpc{"min_ncluster_tpc", 50, "min. number of TPC clusters"};
    Configurable<float> max_chi2tpc{"max_chi2tpc", 4.0, "max. chi2/ncls TPC"};
  } spherocity_cuts;

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  void init(InitContext& initContext)
  {
    getTaskOptionValue(initContext, "create-emevent-dilepton", "applyEveSel_at_skimming", applyEveSel_at_skimming.value, true); // for EM users.
    getTaskOptionValue(initContext, "create-emevent-dilepton", "enable_swt", enable_swt.value, true);                           // for EM users.
    getTaskOptionValue(initContext, "create-emevent-dilepton", "cfg_swt_names", cfg_swt_names.value, true);                     // for EM users.

    if (fillQAHistogram) {
      fRegistry.add("Spherocity/hPt", "pT;p_{T} (GeV/c)", kTH1F, {{200, 0.0f, 10}}, false);
      fRegistry.add("Spherocity/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{180, 0, 2 * M_PI}, {40, -1.0f, 1.0f}}, false);
      fRegistry.add("Spherocity/hSpherocity_ptweighted", "spherocity;Number of used tracks;spherocity", kTH2F, {{101, -0.5, 100.5}, {100, 0.0f, 1}}, false);
      fRegistry.add("Spherocity/hSpherocity_ptunweighted", "spherocity;Number of used tracks;spherocity", kTH2F, {{101, -0.5, 100.5}, {100, 0.0f, 1}}, false);
    }
  }

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber = 0;
  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    if (enable_swt) {
      zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), cfg_swt_names.value);
    }
    mRunNumber = bc.runNumber();
  }

  template <typename TTracks>
  int getSpherocity(TTracks const& tracks, float& spherocity_ptweighted, float& spherocity_ptunweighted)
  {
    // Reference for spherocity : https://arxiv.org/pdf/1905.07208, https://arxiv.org/abs/2310.20406, https://arxiv.org/abs/1205.3963
    int used_ntrack_spherocity = 0;
    float sum_pt = 0.f, sum_ntrack = 0.f;

    for (auto const& track : tracks) {
      if (spherocity_cuts.require_isPVContributor && !track.isPVContributor()) {
        continue;
      }
      if (track.tpcNClsFound() < spherocity_cuts.min_ncluster_tpc) {
        continue;
      }

      sum_pt += track.pt();
      sum_ntrack += 1.f;

      if (fillQAHistogram) {
        fRegistry.fill(HIST("Spherocity/hPt"), track.pt());
        fRegistry.fill(HIST("Spherocity/hEtaPhi"), track.phi(), track.eta());
      }
      used_ntrack_spherocity++;
    } // end of track loop per collision

    float tempSph = 1.f, tempSph_pt1 = 1.f;
    for (int i = 0; i < 360 / 0.1; i++) {
      float nx = std::cos(M_PI / 180.f * i * 0.1);
      float ny = std::sin(M_PI / 180.f * i * 0.1);
      float sum_crossprod = 0.f, sum_crossprod_pt1 = 0.f;
      for (auto const& track : tracks) {
        if (spherocity_cuts.require_isPVContributor && !track.isPVContributor()) {
          continue;
        }
        if (track.tpcNClsFound() < spherocity_cuts.min_ncluster_tpc) {
          continue;
        }
        float px = track.px();
        float py = track.py();
        sum_crossprod += abs(px * ny - py * nx);
        sum_crossprod_pt1 += abs(std::cos(track.phi()) * ny - std::sin(track.phi()) * nx);
      }
      float sph = std::pow(sum_crossprod / sum_pt, 2);
      float sph_pt1 = std::pow(sum_crossprod_pt1 / sum_ntrack, 2);
      if (sph < tempSph) {
        tempSph = sph;
      }
      if (sph_pt1 < tempSph_pt1) {
        tempSph_pt1 = sph_pt1;
      }
    } // end of track loop per collision
    spherocity_ptweighted = std::pow(M_PI_2, 2) * tempSph;
    spherocity_ptunweighted = std::pow(M_PI_2, 2) * tempSph_pt1;
    if (used_ntrack_spherocity < spherocity_cuts.min_ntrack) {
      spherocity_ptweighted = -1.f;
      spherocity_ptunweighted = -1.f;
    }
    return used_ntrack_spherocity;
  }

  Partition<Run3Tracks> tracks_for_spherocity = spherocity_cuts.min_pt < aod::track::pt && spherocity_cuts.min_eta < o2::aod::track::eta && o2::aod::track::eta < spherocity_cuts.max_eta && nabs(o2::aod::track::dcaXY) < spherocity_cuts.max_dcaxy && nabs(o2::aod::track::dcaZ) < spherocity_cuts.max_dcaz && ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::ITS) == true && ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::TPC) == true && o2::aod::track::tpcChi2NCl < spherocity_cuts.max_chi2tpc; // ITS-TPC matched tracks

  void processProp(soa::Join<aod::Collisions, aod::EvSels> const& collisions, Run3Tracks const&, aod::BCsWithTimestamps const&)
  {
    for (auto& collision : collisions) {

      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (applyEveSel_at_skimming && (!collision.selection_bit(o2::aod::evsel::kIsTriggerTVX) || !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder))) {
        continue;
      }
      if (enable_swt && !zorro.isSelected(bc.globalBC())) {
        continue;
      }

      auto tracks_for_spherocity_per_collision = tracks_for_spherocity->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);
      float spherocity_ptweighted = -1.f, spherocity_ptunweighted = -1.f;
      int ntrack = getSpherocity(tracks_for_spherocity_per_collision, spherocity_ptweighted, spherocity_ptunweighted);
      if (fillQAHistogram) {
        fRegistry.fill(HIST("Spherocity/hSpherocity_ptweighted"), ntrack, spherocity_ptweighted);
        fRegistry.fill(HIST("Spherocity/hSpherocity_ptunweighted"), ntrack, spherocity_ptunweighted);
      }
      evprop(spherocity_ptweighted, spherocity_ptunweighted, ntrack);
    } // end of collision loop
  }
  PROCESS_SWITCH(EMEventPropertyTask, processProp, "process event property", true);

  void processDummy(soa::Join<aod::Collisions, aod::EvSels> const& collisions, aod::BCsWithTimestamps const&)
  {
    for (auto& collision : collisions) {
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (applyEveSel_at_skimming && (!collision.selection_bit(o2::aod::evsel::kIsTriggerTVX) || !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder))) {
        continue;
      }
      if (enable_swt && !zorro.isSelected(bc.globalBC())) {
        continue;
      }
      evprop(-1.f, -1.f, 0);
    } // end of collision loop
  }
  PROCESS_SWITCH(EMEventPropertyTask, processDummy, "process dummy", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CreateEMEventDilepton>(cfgc, TaskName{"create-emevent-dilepton"}),
    adaptAnalysisTask<AssociateDileptonToEMEvent>(cfgc, TaskName{"associate-dilepton-to-emevent"}),
    adaptAnalysisTask<EMEventPropertyTask>(cfgc, TaskName{"emevent-property"}),
  };
}
