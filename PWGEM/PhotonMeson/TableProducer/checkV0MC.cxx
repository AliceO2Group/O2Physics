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

/// \brief check MC true V0s
/// \author daiki.sekihata@cern.ch

#include <TMath.h>
#include <TVector2.h>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/Core/trackUtilities.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/PCMUtilities.h"
#include "PWGEM/PhotonMeson/Utils/MCUtilities.h"

using namespace o2;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;

using MyTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksDCA, aod::TracksCovIU>;
using MyTracksMC = soa::Join<MyTracks, aod::McTrackLabels>;

struct checkV0MC {

  Configurable<float> maxX{"maxX", 83.1, "maximum X (starting point of track X)"};
  Configurable<float> maxY{"maxY", 20.0, "maximum X (starting point of track X)"};
  Configurable<float> minpt{"minpt", 0.01, "min pt for track"};
  Configurable<float> maxeta{"maxeta", 0.9, "eta acceptance"};
  Configurable<float> dcamin{"dcamin", 0.1, "dcamin"};
  Configurable<float> dcamax{"dcamax", 1e+10, "dcamax"};

  HistogramRegistry fRegistry{
    "fRegistry",
    {
      {"hDiffCollId", "difference in collision id between ele/pos", {HistType::kTH1F, {{101, -50.5f, +50.5f}}}},
      {"hOrphanTrack", "Number of orphan track from true photon conversions", {HistType::kTH1F, {{1, 0.5f, +1.5f}}}},
      {"hV0PhotonCorrX", "V0 photon track iu X correlation;track iu X for ele;track iu X for pos", {HistType::kTH2F, {{100, 0, +100}, {100, 0, +100}}}},
      {"hV0PhotonCorrY", "V0 photon track iu Y correlation;track iu Y for ele;track iu Y diff", {HistType::kTH2F, {{200, -100, +100}, {400, -20, +20}}}},
      {"hV0PhotonCorrZ", "V0 photon track iu Z correlation;track iu Z for ele;track iu Z diff", {HistType::kTH2F, {{200, -100, +100}, {400, -20, +20}}}},
      {"hV0PhotonCorrYvsX", "V0 photon track iu Y correlation;track iu Y for ele;track iu Y diff", {HistType::kTH2F, {{100, 0, +100}, {400, -20, +20}}}},
      {"hV0PhotonCorrZvsX", "V0 photon track iu Z correlation;track iu Z for ele;track iu Z diff", {HistType::kTH2F, {{100, 0, +100}, {400, -20, +20}}}},
    },
  };

  template <typename TTrack>
  bool isSelected(TTrack const& track)
  {
    if (track.pt() < minpt || abs(track.eta()) > maxeta) {
      return false;
    }
    if (abs(track.dcaXY()) < dcamin || dcamax < abs(track.dcaXY())) {
      return false;
    }
    if (!track.hasITS() && !track.hasTPC()) {
      return false;
    }

    if (track.hasTPC()) {
      if (track.tpcSignal() < 40.f || 110.f < track.tpcSignal()) {
        return false;
      }
    }
    return true;
  }

  Filter trackFilter = o2::aod::track::x < maxX && nabs(o2::aod::track::y) < maxY && o2::aod::track::pt > minpt&& nabs(o2::aod::track::eta) < maxeta&& dcamin < nabs(o2::aod::track::dcaXY) && nabs(o2::aod::track::dcaXY) < dcamax;
  using MyFilteredTracksMC = soa::Filtered<MyTracksMC>;
  Partition<MyFilteredTracksMC> posTracks = o2::aod::track::signed1Pt > 0.f;
  Partition<MyFilteredTracksMC> negTracks = o2::aod::track::signed1Pt < 0.f;
  Partition<MyFilteredTracksMC> orphan_posTracks = o2::aod::track::signed1Pt > 0.f && o2::aod::track::collisionId < 0;
  Partition<MyFilteredTracksMC> orphan_negTracks = o2::aod::track::signed1Pt < 0.f && o2::aod::track::collisionId < 0;

  int ndf = 0;

  using MyCollisions = soa::Join<aod::McCollisionLabels, aod::Collisions>;
  void processMC(MyCollisions const& collisions, aod::McCollisions const&, aod::BCsWithTimestamps const& bcs, MyFilteredTracksMC const& tracks, aod::McParticles const& mcTracks)
  {
    LOGF(info, "n pos tracks = %d , n neg tracks = %d", posTracks.size(), negTracks.size());
    ndf++;
    if (ndf > 3) { // 3 DFs are enough.
      return;
    }

    for (auto& [ele, pos] : combinations(CombinationsFullIndexPolicy(negTracks, posTracks))) {
      if (!ele.has_mcParticle() || !pos.has_mcParticle()) {
        continue;
      }

      if (!isSelected(ele) || !isSelected(pos)) {
        continue;
      }

      auto elemc = ele.template mcParticle_as<aod::McParticles>();
      auto posmc = pos.template mcParticle_as<aod::McParticles>();

      int photonid = FindCommonMotherFrom2Prongs(posmc, elemc, -11, 11, 22, mcTracks);
      if (photonid < 0) { // check swap, true electron is reconstructed as positron and vice versa.
        photonid = FindCommonMotherFrom2Prongs(posmc, elemc, 11, -11, 22, mcTracks);
      }
      if (photonid < 0) {
        continue;
      }
      // At this point, there are only true photon conversions.
      auto photonmc = mcTracks.iteratorAt(photonid);
      if (!IsPhysicalPrimary(photonmc.mcCollision(), photonmc, mcTracks)) {
        continue;
      }

      bool has_coll_ele = ele.has_collision();
      bool has_coll_pos = pos.has_collision();

      if (!has_coll_ele) {
        fRegistry.fill(HIST("hOrphanTrack"), 1);
      }
      if (!has_coll_pos) {
        fRegistry.fill(HIST("hOrphanTrack"), 1);
      }

      if (!has_coll_ele || !has_coll_pos) {
        continue;
      }

      auto collision_ele = ele.collision();
      auto collision_pos = pos.collision();
      int diff = collision_ele.globalIndex() - collision_pos.globalIndex();
      fRegistry.fill(HIST("hDiffCollId"), diff);

      fRegistry.fill(HIST("hV0PhotonCorrX"), ele.x(), pos.x());
      fRegistry.fill(HIST("hV0PhotonCorrY"), ele.y(), pos.y() - ele.y());
      fRegistry.fill(HIST("hV0PhotonCorrZ"), ele.z(), pos.z() - ele.z());
      fRegistry.fill(HIST("hV0PhotonCorrYvsX"), ele.x(), pos.y() - ele.y());
      fRegistry.fill(HIST("hV0PhotonCorrZvsX"), ele.x(), pos.z() - ele.z());

    } // end of pairing loop
  }
  PROCESS_SWITCH(checkV0MC, processMC, "process mc truth info", true);

  void processDummy(soa::Join<aod::McCollisionLabels, aod::Collisions> const& collisions) {}
  PROCESS_SWITCH(checkV0MC, processDummy, "process dummy", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<checkV0MC>(cfgc, TaskName{"check-v0-mc"})};
}
