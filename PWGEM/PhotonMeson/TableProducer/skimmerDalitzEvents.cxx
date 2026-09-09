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

#include <Framework/ASoA.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::common::core;

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::PMEvSels>;
using MyElectrons = o2::aod::EMPrimaryElectronsFromDalitzTmp;
using MyV0Photons = o2::soa::Join<o2::aod::V0PhotonsKFTmp, aod::V0PhotonsPhiVPsiTmp>;
using MyV0PhotonsLegs = o2::soa::Join<o2::aod::V0LegsTmp, aod::V0LegsXYZTmp>;

using MyElectronsMC = o2::soa::Join<o2::aod::EMPrimaryElectronsFromDalitzTmp, aod::EMPrimaryElectronsDeDxMCTmp>;
using MyV0PhotonsLegsMC = o2::soa::Join<o2::aod::V0LegsTmp, aod::V0LegsXYZTmp, aod::V0LegsDeDxMCTmp>;

struct skimmerDalitzEvents {

  Produces<aod::V0PhotonsKF> v0photonskf;
  Produces<aod::V0Legs> v0legs;
  Produces<aod::V0LegsXYZ> v0legsXYZ;
  Produces<aod::V0LegsDeDxMC> v0legsDeDxMC;
  Produces<aod::V0PhotonsPhiVPsi> v0photonsphivpsi;

  Produces<aod::EMPrimaryElectronsFromDalitz> emprimaryelectrons;
  Produces<aod::EMPrimaryElectronsDeDxMC> emprimaryelectronsDeDxMC;

  Configurable<unsigned int> minNElecCand{"minNElecCand", 2, "minimum number of electrons/positron candidates in one event"};
  Configurable<unsigned int> minNGammaCand{"minNGammaCand", 1, "minimum number of V0 photon candidates in one event"};

  // ---------- for data ----------
  template <bool isMC, typename TElectrons, typename TV0Photons, typename TV0Legs>
  void process(MyCollisions const& collisions, TElectrons const& emPrimaryElecTmp, TV0Photons const& v0photonskfTmp, TV0Legs const& v0LegsTmp)
  {
    PresliceUnsorted<TElectrons> perCollisionEl = aod::emprimaryelectron::collisionId;
    Preslice<TV0Photons> perCol_pcm = o2::aod::v0photonkf::collisionId;
    Preslice<TV0Legs> perCol_legs = o2::aod::v0leg::collisionId;
    for (const auto& collision : collisions) {
      auto tracks = emPrimaryElecTmp.sliceBy(perCollisionEl, collision.globalIndex()); // o2::aod::track::collisionId // o2::aod::emprimaryelectronda::pmeventId

      if (tracks.size() < minNElecCand) {
        continue;
      }

      auto photonskf = v0photonskfTmp.sliceBy(perCol_pcm, collision.globalIndex());

      if (photonskf.size() < minNGammaCand) {
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

        if constexpr (isMC) {
          emprimaryelectronsDeDxMC(
            track.mcTunedTPCSignal());
        }
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

        v0photonsphivpsi(
          v0.phiv(), v0.psipair());
      }

      auto v0LegsTmpPerColl = v0LegsTmp.sliceBy(perCol_legs, collision.globalIndex());
      for (const auto& leg : v0LegsTmp) {
        v0legs(
          leg.collisionId(), leg.trackId(), leg.sign(),
          leg.px(), leg.py(), leg.pz(),
          leg.dcaXY(), leg.dcaZ(),
          leg.tpcNClsFindable(), leg.tpcNClsFindableMinusFound(), leg.tpcNClsFindableMinusCrossedRows(), leg.tpcNClsShared(),
          leg.tpcChi2NCl(), leg.tpcInnerParam(),
          leg.tpcSignal(), leg.tpcNSigmaEl(), leg.tpcNSigmaPi(),
          leg.itsClusterSizes(), leg.itsChi2NCl(), leg.detectorMap());

        v0legsXYZ(
          leg.x(), leg.y(), leg.z());

        if constexpr (isMC) {
          v0legsDeDxMC(
            leg.mcTunedTPCSignal());
        }
      }
    }
  }

  void processRec(MyCollisions const& collisions, MyElectrons const& emPrimaryElecTmp, MyV0Photons const& v0photonskfTmp, MyV0PhotonsLegs const& v0LegsTmp)
  {
    process<false, MyElectrons, MyV0Photons, MyV0PhotonsLegs>(collisions, emPrimaryElecTmp, v0photonskfTmp, v0LegsTmp);
  }

  void processMC(MyCollisions const& collisions, MyElectronsMC const& emPrimaryElecTmp, MyV0Photons const& v0photonskfTmp, MyV0PhotonsLegsMC const& v0LegsTmp)
  {
    process<true, MyElectronsMC, MyV0Photons, MyV0PhotonsLegsMC>(collisions, emPrimaryElecTmp, v0photonskfTmp, v0LegsTmp);
  }

  PROCESS_SWITCH(skimmerDalitzEvents, processRec, "process reconstructed info only", false);  // data
  PROCESS_SWITCH(skimmerDalitzEvents, processMC, "process reconstructed and MC info", false); // MC
};

WorkflowSpec defineDataProcessing(ConfigContext const& context)
{
  return WorkflowSpec{
    adaptAnalysisTask<skimmerDalitzEvents>(context, TaskName{"skimmer-dalitz-events"})};
}
