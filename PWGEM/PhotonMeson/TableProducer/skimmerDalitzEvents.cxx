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

/// \file skimmerDalitzEvents.cxx
/// \brief write tables for photons and electrons for dalitz decay
/// \author josuha.konig@cern.ch

#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"

#include <Framework/runDataProcessing.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::common::core;

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::PMEvSels>;

struct skimmerDalitzEvents {

  o2::framework::SliceCache cache;

  PresliceUnsorted<aod::EMPrimaryElectronsFromDalitz> perCollisionEl = aod::emprimaryelectron::collisionId;
  Preslice<aod::V0PhotonsKF> perCol_pcm = o2::aod::v0photonkf::collisionId;

  Produces<aod::V0PhotonsKF> v0photonskf;
  Produces<aod::V0Legs> v0legs;
  Produces<aod::V0LegsXYZ> v0legsXYZ;
  Produces<aod::V0LegsDeDxMC> v0legsDeDxMC;
  Produces<aod::V0PhotonsPhiVPsi> v0photonsphivpsi;

  Produces<aod::EMPrimaryElectronsFromDalitz> emprimaryelectrons;
  Produces<aod::EMPrimaryElectronsDeDxMC> emprimaryelectronsDeDxMC;

  // ---------- for data ----------
  template <bool isMC>
  void process(MyCollisions const& collisions, aod::EMPrimaryElectronsFromDalitzTmp const& emPrimaryElecTmp, aod::V0PhotonsKFTmp const& v0photonskfTmp, aod::V0LegsTmp const& v0LegsTmp, aod::V0LegsXYZTmp const& v0legsXYZTmp, aod::V0PhotonsPhiVPsiTmp const& v0PhiVPsiTmp, aod::EMPrimaryElectronsDeDxMCTmp const* electronsDeDxMCTmp, aod::V0LegsDeDxMCTmp const* v0LegsDeDxMCTmp)
  {
    for (const auto& collision : collisions) {
      // auto el_per_coll = emPrimaryElecTmp.sliceBy(perCollisionEl, collision.globalIndex());
      auto tracks = emPrimaryElecTmp.sliceByCached(o2::aod::emprimaryelectronda::pmeventId, collision.globalIndex(), cache);

      if (tracks.size() < 2) {
        continue;
      }

      auto photonskf = v0photonskfTmp.sliceByCached(o2::aod::emprimaryelectronda::pmeventId, collision.globalIndex(), cache);

      if (photonskf.size() == 0) {
        continue;
      }

      // Write electron tracks
      for (const auto& track : tracks) {
        emprimaryelectrons(track.collisionId(), track.trackId(), track.sign(),
                           track.pt(), track.eta(), track.phi(),
                           track.dcaXY(), track.dcaZ(), track.cYY(), track.cZY(), track.cZZ(),
                           track.tpcNClsFindable(), track.tpcNClsFindableMinusFound(), track.tpcNClsFindableMinusCrossedRows(), track.tpcNClsShared(),
                           track.tpcChi2NCl(), track.tpcInnerParam(),
                           track.tpcSignal(), track.tpcNSigmaEl(), track.tpcNSigmaPi(),
                           track.beta(), track.tofNSigmaEl(),
                           track.itsClusterSizes(), track.itsChi2NCl(), track.tofChi2(), track.detectorMap());
      }

      for (const auto& v0 : photonskf) {
        v0photonskf(v0.collisionId(), v0.v0Id(), v0.posTrack(), v0.negTrack(),
                    v0.vx(), v0.vy(), v0.vz(),
                    v0.px(), v0.py(), v0.pz(),
                    v0.mGamma(),
                    v0.dcaXYtopv(), v0.dcaZtopv(),
                    v0.cospa(), v0.cospaXY(), v0.cospaRZ(), v0.pca(),
                    v0.alpha(), v0.qtarm(),
                    v0.chiSquareNDF());
      }

      // Now load all other tables

      auto v0LegsTmpPerColl = v0LegsTmp.sliceByCached(o2::aod::emprimaryelectronda::pmeventId, collision.globalIndex(), cache);
      for (const auto& leg : v0LegsTmp) {
        v0legs(
          leg.collisionId(), leg.trackId(), leg.sign(),
          leg.px(), leg.py(), leg.pz(),
          leg.dcaXY(), leg.dcaZ(),
          leg.tpcNClsFindable(), leg.tpcNClsFindableMinusFound(), leg.tpcNClsFindableMinusCrossedRows(), leg.tpcNClsShared(),
          leg.tpcChi2NCl(), leg.tpcInnerParam(),
          leg.tpcSignal(), leg.tpcNSigmaEl(), leg.tpcNSigmaPi(),
          leg.itsClusterSizes(), leg.itsChi2NCl(), leg.detectorMap());
      }

      auto v0legsXYZTmpPerColl = v0legsXYZTmp.sliceByCached(o2::aod::emprimaryelectronda::pmeventId, collision.globalIndex(), cache);
      for (const auto& leg : v0legsXYZTmpPerColl) {
        v0legsXYZ(
          leg.x(), leg.y(), leg.z());
      }

      auto v0legsPhiVPsiTmpPerColl = v0PhiVPsiTmp.sliceByCached(o2::aod::emprimaryelectronda::pmeventId, collision.globalIndex(), cache);
      for (const auto& leg : v0legsPhiVPsiTmpPerColl) {
        v0photonsphivpsi(
          leg.phiv(), leg.psipair());
      }

      if constexpr (isMC) {
        auto electronsDeDxMCTmpPerColl = (*electronsDeDxMCTmp).sliceByCached(o2::aod::emprimaryelectronda::pmeventId, collision.globalIndex(), cache);
        for (const auto& track : electronsDeDxMCTmpPerColl) {
          emprimaryelectronsDeDxMC(
            track.mcTunedTPCSignal());
        }

        auto v0legsDeDxMCTmpPerColl = (*v0LegsDeDxMCTmp).sliceByCached(o2::aod::emprimaryelectronda::pmeventId, collision.globalIndex(), cache);
        for (const auto& leg : v0legsDeDxMCTmpPerColl) {
          v0legsDeDxMC(
            leg.mcTunedTPCSignal());
        }
      }
    }
  }

  void processRec(MyCollisions const& collisions, aod::EMPrimaryElectronsFromDalitzTmp const& emPrimaryElecTmp, aod::V0PhotonsKFTmp const& v0photonskfTmp, aod::V0LegsTmp const& v0LegsTmp, aod::V0LegsXYZTmp const& v0legsXYZTmp, aod::V0PhotonsPhiVPsiTmp const& v0PhiVPsiTmp)
  {
    process<false>(collisions, emPrimaryElecTmp, v0photonskfTmp, v0LegsTmp, v0legsXYZTmp, v0PhiVPsiTmp, nullptr, nullptr);
  }

  void processMC(MyCollisions const& collisions, aod::EMPrimaryElectronsFromDalitzTmp const& emPrimaryElecTmp, aod::V0PhotonsKFTmp const& v0photonskfTmp, aod::V0LegsTmp const& v0LegsTmp, aod::V0LegsXYZTmp const& v0legsXYZTmp, aod::V0PhotonsPhiVPsiTmp const& v0PhiVPsiTmp, aod::EMPrimaryElectronsDeDxMCTmp const& electronsDeDxMCTmp, aod::V0LegsDeDxMCTmp const& v0LegsDeDxMCTmp)
  {
    process<true>(collisions, emPrimaryElecTmp, v0photonskfTmp, v0LegsTmp, v0legsXYZTmp, v0PhiVPsiTmp, &electronsDeDxMCTmp, &v0LegsDeDxMCTmp);
  }

  PROCESS_SWITCH(skimmerDalitzEvents, processRec, "process reconstructed info only", false);  // data
  PROCESS_SWITCH(skimmerDalitzEvents, processMC, "process reconstructed and MC info", false); // MC
};

WorkflowSpec defineDataProcessing(ConfigContext const& context)
{
  return WorkflowSpec{
    adaptAnalysisTask<skimmerDalitzEvents>(context, TaskName{"skimmer-dalitz-events"})};
}
