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

/// \file photonChargedTriggerProducer.cxx
/// \author Julius Kinner
/// \brief photon-jet angular correlation table producer
///
/// Table producer for photon-jet angular correlation analysis (see photonChargedTriggerCorrelation.cxx)

#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/PCMUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/PhotonChargedTriggerCorrelation.h"

#include "Common/Core/TableHelper.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

#include "Math/Vector4D.h"
#include "TMath.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <deque>
#include <memory>
#include <random>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;

// correlation derived data ===================================================================================================================================================================

struct PhotonChargedTriggerProducer {
  // reco
  Produces<aod::CollisionsExtraCorr> collisionExtraCorrTable;
  Produces<aod::Triggers> triggerTable;
  Produces<aod::Hadrons> hadronTable;
  Produces<aod::Pipms> pipmTable;
  Produces<aod::PhotonPCMs> photonPCMTable;
  Produces<aod::PhotonPCMPairs> photonPCMPairTable;
  // mc
  Produces<aod::McCollisionsExtraCorr> mcCollisionExtraCorrTable;
  Produces<aod::TriggerParticles> triggerParticleTable;

  Configurable<double> zPvMax{"zPvMax", 7, "maximum absZ primary-vertex cut"};
  Configurable<int> occupancyMin{"occupancyMin", 0, "minimum occupancy cut"};
  Configurable<int> occupancyMax{"occupancyMax", 2000, "maximum occupancy cut"};
  Configurable<double> etaMax{"etaMax", 0.8, "maximum absEta cut"};

  Configurable<std::string> eventSelections{"eventSelections", "sel8", "JE framework - event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "JE framework - track selections"};
  Configurable<std::string> triggerMasks{"triggerMasks", "", "JE framework - skimmed data trigger masks (relevent for correlation: fTrackLowPt,fTrackHighPt)"};

  Configurable<double> piPIDLowPt{"piPIDLowPt", 0.5, "max pt value for pipm PID without tof"};
  Configurable<double> piPIDHighPt{"piPIDHighPt", 2.5, "min pt value for pipm PID without tof in relativistic rise of Bethe-Bloch"};
  Configurable<std::vector<double>> nSigmaPiTpcLowPt{"nSigmaPiTpcLowPt", {-2, 2}, "minimum-maximum nSigma for pipm in tpc at low pt"};
  Configurable<std::vector<double>> nSigmaPiTpcMidPt{"nSigmaPiTpcMidPt", {-1, 1}, "minimum-maximum nSigma for pipm in tpc at mid pt"};
  Configurable<std::vector<double>> nSigmaPiTof{"nSigmaPiTof", {-1, 2}, "minimum-maximum nSigma for pipm in tof"};
  Configurable<std::vector<double>> nSigmaPiRelRise{"nSigmaPiRelRise", {0, 2}, "minimum-maximum nSigma pipm tpc at high pt"};

  Configurable<float> ptTrigMin{"ptTrigMin", 5, "minimum pT of triggers"};

  // derivatives of configurables

  std::vector<int> eventSelectionBits;
  int trackSelection = -1;
  std::vector<int> triggerMaskBits;

  // for mc
  Service<framework::O2DatabasePDG> pdg;

  // partitions++
  SliceCache cache;

  Preslice<aod::JetTracks> perColTracks = aod::jtrack::collisionId;
  Preslice<aod::JetParticles> perColMcParticles = aod::jmcparticle::mcCollisionId;

  Preslice<aod::V0PhotonsKF> perColV0Photons = aod::v0photonkf::collisionId;

  // functions ================================================================================================================================================================================

  // selections ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // event selection
  template <typename T_collision>
  bool checkEventSelection(T_collision const& collision)
  {
    if (!jetderiveddatautilities::selectTrigger(collision, triggerMaskBits))
      return false;
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits))
      return false;
    if (std::abs(collision.posZ()) > zPvMax)
      return false;
    if (collision.trackOccupancyInTimeRange() < occupancyMin || collision.trackOccupancyInTimeRange() > occupancyMax)
      return false;
    return true;
  }

  // checks global track cuts
  template <typename T_track>
  bool checkGlobalTrackEta(T_track const& track)
  {
    if (!jetderiveddatautilities::selectTrack(track, trackSelection))
      return false;
    if (!jetderiveddatautilities::applyTrackKinematics(track, 0.1, 1000, -1 * etaMax, etaMax))
      return false;
    return true;
  }

  // checks pipm selection (just PID (no additional track cuts))
  template <typename T_track>
  bool checkPipmTPCTOF(T_track const& track)
  {
    // too low for tof
    if (track.pt() < piPIDLowPt) {
      if (track.tpcNSigmaPi() > nSigmaPiTpcLowPt.value[0] && track.tpcNSigmaPi() < nSigmaPiTpcLowPt.value[1]) {
        return true;
      }
      return false;
    }
    // Bethe-Bloch overlap (-> tpc + tof)
    if (track.pt() < piPIDHighPt) {
      if (track.hasTOF()) { // has to stay inside pt-if due to return-layout of function
        if (track.tpcNSigmaPi() > nSigmaPiTpcMidPt.value[0] && track.tpcNSigmaPi() < nSigmaPiTpcMidPt.value[1] &&
            track.tofNSigmaPi() > nSigmaPiTof.value[0] && track.tofNSigmaPi() < nSigmaPiTof.value[1]) {
          return true;
        }
      }
      return false;
    }
    // Bethe-Bloch rel rise (too high for tof)
    if (track.tpcNSigmaPi() > nSigmaPiRelRise.value[0] && track.tpcNSigmaPi() < nSigmaPiRelRise.value[1]) {
      return true;
    }
    return false;
  }

  // checks pipm selection (just PID (no additional track cuts))
  template <typename T_track>
  bool checkPipmTPC(T_track const& track)
  {
    // Bethe-Bloch rel rise
    if (track.pt() > piPIDHighPt) {
      if (track.tpcNSigmaPi() > nSigmaPiRelRise.value[0] && track.tpcNSigmaPi() < nSigmaPiRelRise.value[1]) {
        return true;
      }
    }
    return false;
  }

  // analysis /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void init(InitContext const&)
  {
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
    triggerMaskBits = jetderiveddatautilities::initialiseTriggerMaskBits(triggerMasks);
  }

  void processRecoCollisionTrigger(aod::JetCollision const& collision, aod::JetTracks const& tracks)
  {
    // event selection
    const bool isSelectedEvent = checkEventSelection(collision);
    // trigger event check
    bool isTriggerEvent = false;
    // number global tracks
    int nGlobalTracks = 0;

    // count global tracks (for independence of multiplicity task (uses only JE derieved data))
    for (auto const& track : tracks) {
      // track selection
      if (!checkGlobalTrackEta(track))
        continue;

      nGlobalTracks++;

      if (!isSelectedEvent)
        continue;
      if (track.pt() < ptTrigMin)
        continue;

      isTriggerEvent = true;

      // trigger info
      triggerTable(track.collisionId(), track.globalIndex(), track.pt(), track.phi(), track.eta());
    }

    // collision info
    collisionExtraCorrTable(isSelectedEvent, isTriggerEvent, nGlobalTracks);
  }
  PROCESS_SWITCH(PhotonChargedTriggerProducer, processRecoCollisionTrigger, "process correlation collision_extra and trigger table (reconstructed)", false);

  void processRecoPipmTPCTOF(aod::JetCollision const& collision,
                             soa::Join<aod::JetTracks, aod::JTrackPIs> const& tracks, soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCPi, aod::pidTOFPi> const&)
  {
    // event selection
    if (!checkEventSelection(collision))
      return;

    // hadron/pipm
    for (auto const& track : tracks) {
      // track selection
      if (!checkGlobalTrackEta(track))
        continue;

      // hadron
      hadronTable(track.collisionId(), track.globalIndex(), track.pt(), track.phi(), track.eta());

      // pipm selection
      auto const& trackPID = track.track_as<soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCPi, aod::pidTOFPi>>();
      if (!checkPipmTPCTOF(trackPID))
        continue;

      // pipm
      pipmTable(track.collisionId(), track.globalIndex(), track.pt(), track.phi(), track.eta());
    }
  }
  PROCESS_SWITCH(PhotonChargedTriggerProducer, processRecoPipmTPCTOF, "process pipm (TPC-TOF) table (reconstructed)", false);

  void processRecoPipmTPC(aod::JetCollision const& collision,
                          soa::Join<aod::JetTracks, aod::JTrackPIs> const& tracks, soa::Join<aod::Tracks, aod::pidTPCPi> const&)
  {
    // event selection
    if (!checkEventSelection(collision))
      return;

    // hadron/pipm
    for (auto const& track : tracks) {
      // track selection
      if (!checkGlobalTrackEta(track))
        continue;

      // hadron
      hadronTable(track.collisionId(), track.globalIndex(), track.pt(), track.phi(), track.eta());

      // pipm selection
      auto const& trackPID = track.track_as<soa::Join<aod::Tracks, aod::pidTPCPi>>();
      if (!checkPipmTPC(trackPID))
        continue;

      // pipm
      pipmTable(track.collisionId(), track.globalIndex(), track.pt(), track.phi(), track.eta());
    }
  }
  PROCESS_SWITCH(PhotonChargedTriggerProducer, processRecoPipmTPC, "process pipm (TPC) table (reconstructed)", false);

  void processRecoPhotonPCM(soa::Join<aod::JetCollisions, aod::JCollisionPIs>::iterator const& collision,
                            aod::V0PhotonsKF const& v0Photons, aod::V0Legs const&)
  {
    // event selection
    if (!checkEventSelection(collision))
      return;

    // photonsPCM (for some reason collsionId not an index column (?))
    auto const v0PhotonsThisEvent = v0Photons.sliceBy(perColV0Photons, collision.collisionId());

    // photonPCM
    for (auto const& v0Photon : v0PhotonsThisEvent) {
      // photon selection
      if (std::abs(v0Photon.eta()) > etaMax)
        continue;

      // photon PCM
      photonPCMTable(v0Photon.collisionId(), v0Photon.globalIndex(),
                     v0Photon.posTrack().trackId(), v0Photon.negTrack().trackId(), v0Photon.pt(), v0Photon.phi(), v0Photon.eta());
    }

    // photonPCm pairs
    for (auto const& [v0Photon1, v0Photon2] : soa::combinations(soa::CombinationsStrictlyUpperIndexPolicy(v0PhotonsThisEvent, v0PhotonsThisEvent))) {
      // get kinematics
      ROOT::Math::PtEtaPhiMVector const p4V0PCM1(v0Photon1.pt(), v0Photon1.eta(), v0Photon1.phi(), 0.);
      ROOT::Math::PtEtaPhiMVector const p4V0PCM2(v0Photon2.pt(), v0Photon2.eta(), v0Photon2.phi(), 0.);
      ROOT::Math::PtEtaPhiMVector const p4V0PCMPair = p4V0PCM1 + p4V0PCM2;

      // pi0 selection
      if (std::abs(p4V0PCMPair.Eta()) > etaMax)
        continue;

      // save info
      photonPCMPairTable(v0Photon1.collisionId(), v0Photon1.globalIndex(), v0Photon2.globalIndex(),
                         v0Photon1.posTrack().trackId(), v0Photon1.negTrack().trackId(), v0Photon2.posTrack().trackId(), v0Photon2.negTrack().trackId(),
                         p4V0PCMPair.Pt(), RecoDecay::constrainAngle(p4V0PCMPair.Phi(), 0), p4V0PCMPair.Eta(), p4V0PCMPair.M());
    }
  }
  PROCESS_SWITCH(PhotonChargedTriggerProducer, processRecoPhotonPCM, "process photonPCM table (reconstructed)", false);

  void processMcCorrTables(aod::JetMcCollision const&, aod::JetParticles const& mcParticles)
  {
    // trigger event check
    bool isTriggerEvent = false;
    // number charged particles in eta range
    int nCharged = 0;

    // particle loop
    for (auto const& mcParticle : mcParticles) {
      // track selection
      auto const pdgParticle = pdg->GetParticle(mcParticle.pdgCode());
      if (!pdgParticle || pdgParticle->Charge() == 0)
        continue;
      if (!mcParticle.isPhysicalPrimary())
        continue;
      if (std::abs(mcParticle.eta()) > etaMax)
        continue;

      nCharged++;

      // trigger selection
      if (mcParticle.pt() < ptTrigMin)
        continue;

      isTriggerEvent = true;

      // trigger info
      triggerParticleTable(mcParticle.mcCollisionId(), mcParticle.globalIndex(), mcParticle.pt(), mcParticle.phi(), mcParticle.eta());
    }

    // collision info
    mcCollisionExtraCorrTable(isTriggerEvent, nCharged);
  }
  PROCESS_SWITCH(PhotonChargedTriggerProducer, processMcCorrTables, "process table production (mc)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& configContext)
{
  return WorkflowSpec{adaptAnalysisTask<PhotonChargedTriggerProducer>(configContext)};
}
