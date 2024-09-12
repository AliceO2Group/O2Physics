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
#include "Common/Core/TableHelper.h"

#include "PWGEM/Dilepton/DataModel/dileptonTables.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

using MyBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;
using MyQvectors = soa::Join<aod::QvectorFT0CVecs, aod::QvectorFT0AVecs, aod::QvectorFT0MVecs, aod::QvectorBPosVecs, aod::QvectorBNegVecs, aod::QvectorBTotVecs>;

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::EMEvSels, aod::Mults>;
using MyCollisions_Cent = soa::Join<MyCollisions, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentNTPVs>; // centrality table has dependency on multiplicity table.
using MyCollisions_Cent_Qvec = soa::Join<MyCollisions_Cent, MyQvectors>;

using MyCollisionsWithSWT = soa::Join<MyCollisions, aod::EMSWTriggerInfosTMP>;
using MyCollisionsWithSWT_Cent = soa::Join<MyCollisionsWithSWT, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentNTPVs>; // centrality table has dependency on multiplicity table.
using MyCollisionsWithSWT_Cent_Qvec = soa::Join<MyCollisionsWithSWT_Cent, MyQvectors>;

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

  HistogramRegistry registry{"registry"};
  void init(o2::framework::InitContext&)
  {
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    auto hEventCounter = registry.add<TH1>("hEventCounter", "hEventCounter", kTH1I, {{7, 0.5f, 7.5f}});
    hEventCounter->GetXaxis()->SetBinLabel(1, "all");
    hEventCounter->GetXaxis()->SetBinLabel(2, "sel8");

    registry.add("hNInspectedTVX", "N inspected TVX;run number;N_{TVX}", kTProfile, {{80000, 520000.5, 600000.5}}, true);
  }

  ~CreateEMEventDilepton() {}

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

  Preslice<aod::V0PhotonsKF> perCollision_pcm = aod::v0photonkf::collisionId;
  PresliceUnsorted<aod::EMPrimaryElectrons> perCollision_el = aod::emprimaryelectron::collisionId;
  PresliceUnsorted<aod::EMPrimaryMuons> perCollision_mu = aod::emprimarymuon::collisionId;

  template <bool isMC, bool isTriggerAnalysis, EMEventType eventype, typename TCollisions, typename TBCs>
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

      if (!collision.isSelected()) {
        continue;
      }

      if constexpr (isTriggerAnalysis) {
        if (collision.swtaliastmp_raw() == 0) {
          continue;
        } else {
          emswtbit(collision.swtaliastmp_raw(), collision.nInspectedTVX());
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
          999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f);
      } else if constexpr (eventype == EMEventType::kEvent_Cent) {
        event_cent(collision.centFT0M(), collision.centFT0A(), collision.centFT0C(), collision.centNTPV());
        event_qvec(
          999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f,
          999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f);
      } else if constexpr (eventype == EMEventType::kEvent_Cent_Qvec) {
        event_cent(collision.centFT0M(), collision.centFT0A(), collision.centFT0C(), collision.centNTPV());
        float q2xft0m = 999.f, q2yft0m = 999.f, q2xft0a = 999.f, q2yft0a = 999.f, q2xft0c = 999.f, q2yft0c = 999.f, q2xbpos = 999.f, q2ybpos = 999.f, q2xbneg = 999.f, q2ybneg = 999.f, q2xbtot = 999.f, q2ybtot = 999.f;
        float q3xft0m = 999.f, q3yft0m = 999.f, q3xft0a = 999.f, q3yft0a = 999.f, q3xft0c = 999.f, q3yft0c = 999.f, q3xbpos = 999.f, q3ybpos = 999.f, q3xbneg = 999.f, q3ybneg = 999.f, q3xbtot = 999.f, q3ybtot = 999.f;

        if (collision.qvecFT0CReVec().size() >= 2) { // harmonics 2,3
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
          q3xft0m, q3yft0m, q3xft0a, q3yft0a, q3xft0c, q3yft0c, q3xbpos, q3ybpos, q3xbneg, q3ybneg, q3xbtot, q3ybtot);
      } else {
        event_cent(105.f, 105.f, 105.f, 105.f);
        event_qvec(
          999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f,
          999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f, 999.f);
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

  //---------- for data with swt----------

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
    } // end of collision loop
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
  void init(InitContext&)
  {
    if (fillQAHistogram) {
      fRegistry.add("Spherocity/hPt", "pT;p_{T} (GeV/c)", kTH1F, {{200, 0.0f, 10}}, false);
      fRegistry.add("Spherocity/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{180, 0, 2 * M_PI}, {40, -1.0f, 1.0f}}, false);
      fRegistry.add("Spherocity/hSpherocity_ptweighted", "spherocity;Number of used tracks;spherocity", kTH2F, {{101, -0.5, 100.5}, {100, 0.0f, 1}}, false);
      fRegistry.add("Spherocity/hSpherocity_ptunweighted", "spherocity;Number of used tracks;spherocity", kTH2F, {{101, -0.5, 100.5}, {100, 0.0f, 1}}, false);
    }
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

  void processProp(aod::EMEvents const& collisions, Run3Tracks const&)
  {
    for (auto& collision : collisions) {
      auto tracks_for_spherocity_per_collision = tracks_for_spherocity->sliceByCached(o2::aod::track::collisionId, collision.collisionId(), cache);
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

  void processDummy(aod::EMEvents const& collisions)
  {
    for (int i = 0; i < collisions.size(); i++) {
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
