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

/// \file photonChargedTriggerCorrelation.cxx
/// \author Julius Kinner
/// \brief photon-jet correlation analysis
///
/// Analysis for angular correlations between jets and photons via two-particle correlations with charged high-pt triggers
/// Associated hadrons (tracks), pipm, photons (PCM), pi0 (PCM)
/// Also contains checks and monte-carlo (efficiency, purity, mc-true correlation,...)
/// End goal of studying correlations between direct photons and jets

#include <cmath>
#include <deque>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>
#include <memory>
#include <random>

#include "TMath.h"
#include "Math/Vector4D.h"

#include "CCDB/BasicCCDBManager.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/TableHelper.h"

#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/PCMUtilities.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/PhotonChargedTriggerCorrelation.h"

const double absEtaMax = 0.8;
#define DPHI_SCALE constants::math::TwoPI - constants::math::PIHalf
#define DETA_SCALE 4 * absEtaMax - 2 * absEtaMax

using namespace o2;
using namespace o2::framework;

using CorrCollisions = soa::Join<aod::JetCollisions, aod::CollisionsExtraCorr, aod::MultsGlobal>;
using CorrCollision = CorrCollisions::iterator;
using CorrMcDCollisions = soa::Join<aod::JetCollisionsMCD, aod::CollisionsExtraCorr, aod::MultsGlobal>;
using CorrMcDCollision = CorrMcDCollisions::iterator;
using CorrMcCollisions = soa::Join<aod::JetMcCollisions, aod::McCollisionsExtraCorr, aod::MultMCExtras>;
using CorrMcCollision = CorrMcCollisions::iterator;

using BinningZPvMult = ColumnBinningPolicy<aod::jcollision::PosZ, aod::mult::MultNTracksGlobal>;

// correlation derived data ===================================================================================================================================================================

struct CorrelationTableProducer {
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
  Configurable<double> etaMax{"etaMax", 1 * absEtaMax, "maximum absEta cut"};

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
  Partition<aod::JetTracks> partitionTriggerTracks = aod::jtrack::pt > ptTrigMin;
  Partition<aod::JetParticles> partitionTriggerParticles = aod::jmcparticle::pt > ptTrigMin;

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

  void processRecoCollisionTrigger(aod::JetCollision const& collision, aod::JetTracks const&)
  {
    // event selection
    const bool isSelectedEvent = checkEventSelection(collision);
    // trigger event check
    bool isTriggerEvent = false;

    if (isSelectedEvent) {
      // group collision
      auto const triggers = partitionTriggerTracks->sliceByCached(aod::jtrack::collisionId, collision.globalIndex(), cache);

      // trigger loop
      for (auto const& trigger : triggers) {
        // track selection
        if (!checkGlobalTrackEta(trigger))
          continue;

        // detect trigger event
        isTriggerEvent = true;

        // trigger info
        triggerTable(trigger.collisionId(), trigger.globalIndex(), trigger.pt(), trigger.phi(), trigger.eta());
      }
    }

    // collision info
    collisionExtraCorrTable(isSelectedEvent, isTriggerEvent);
  }
  PROCESS_SWITCH(CorrelationTableProducer, processRecoCollisionTrigger, "process correlation collision_extra and trigger table (reconstructed)", false);

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
  PROCESS_SWITCH(CorrelationTableProducer, processRecoPipmTPCTOF, "process pipm (TPC-TOF) table (reconstructed)", false);

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
  PROCESS_SWITCH(CorrelationTableProducer, processRecoPipmTPC, "process pipm (TPC) table (reconstructed)", false);

  void processRecoPhotonPCM(soa::Join<aod::JetCollisions, aod::JCollisionPIs>::iterator const& collision, aod::Collisions const&,
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
                         p4V0PCMPair.Pt(), p4V0PCMPair.Phi() + constants::math::PI, p4V0PCMPair.Eta(), p4V0PCMPair.M());
    }
  }
  PROCESS_SWITCH(CorrelationTableProducer, processRecoPhotonPCM, "process photonPCM table (reconstructed)", false);

  void processMcCorrTables(aod::JetMcCollision const& mcCollision, aod::JetParticles const&)
  {
    // group collision
    auto const triggers = partitionTriggerParticles->sliceByCached(aod::jmcparticle::mcCollisionId, mcCollision.globalIndex(), cache);
    // trigger event check
    bool isTriggerEvent = false;

    // trigger loop
    for (auto const& trigger : triggers) {
      // track selection
      auto const pdgParticle = pdg->GetParticle(trigger.pdgCode());
      if (!pdgParticle || pdgParticle->Charge() == 0)
        continue;
      if (!trigger.isPhysicalPrimary())
        continue;
      if (std::abs(trigger.eta()) > etaMax)
        continue;

      // detect trigger event
      isTriggerEvent = true;

      // trigger info
      triggerParticleTable(mcCollision.globalIndex(), trigger.globalIndex(), trigger.pt(), trigger.phi(), trigger.eta());
    }

    // collision info
    mcCollisionExtraCorrTable(isTriggerEvent);
  }
  PROCESS_SWITCH(CorrelationTableProducer, processMcCorrTables, "process table production (mc)", false);
};

// correlation analysis =======================================================================================================================================================================

struct PhotonChargedTriggerCorrelation {
  // configurables

  // general (kenobi)
  Configurable<std::string> pathCcdbEff{"pathCcdbEff", "Users/j/jkinner/efficiency/set_in_config", "base path to the ccdb efficiencies"};
  Configurable<std::string> urlCcdb{"urlCcdb", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<int64_t> noLaterThanCcdb{"noLaterThanCcdb",
                                        std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(),
                                        "latest acceptable timestamp of creation for the object"};

  // analysis
  Configurable<bool> doEffCorrectionTrigger{"doEffCorrectionTrigger", false, "whether to do on-the-fly mixing correction for triggers"};
  Configurable<bool> doEffCorrectionHadron{"doEffCorrectionHadron", false, "whether to do on-the-fly mixing correction for hadrons"};
  Configurable<bool> doEffCorrectionPipm{"doEffCorrectionPipm", false, "whether to do on-the-fly mixing correction for pipm"};
  Configurable<bool> doEffCorrectionPhotonPCM{"doEffCorrectionPhotonPCM", false, "whether to do on-the-fly mixing correction for photonPCM"};

  Configurable<bool> doTrigEvMixing{"doTrigEvMixing", false, "whether to use trigger events for trigger mixing"};
  Configurable<bool> doTrigEvEff{"doTrigEvEff", false, "whether to use trigger events for efficiency histograms"};
  Configurable<int> nTriggerSavedForMixing{"nTriggerSavedForMixing", 2048, "number of triggers that are saved for mixing with other events"};
  Configurable<int> nTriggerMixingHadron{"nTriggerMixingHadron", 64, "number of triggers that are used for hadron mixing"};
  Configurable<int> nTriggerMixingPipm{"nTriggerMixingPipm", 64, "number of triggers that are used for pipm mixing"};
  Configurable<int> nTriggerMixingPhotonPCM{"nTriggerMixingPhotonPCM", 64, "number of triggers that are saved for photonPCM mixing"};
  Configurable<int> nTriggerMixingPi0PCM{"nTriggerMixingPi0PCM", 64, "number of triggers that are saved for pi0PCM mixing"};
  Configurable<int> nNeighboursMixingPi0PCMPair{"nNeighboursMixingPi0PCMPair", 64, "number neighbours used for for pi0PCM pair mixing"};
  Configurable<std::vector<double>> pi0PCMMassRange{"pi0PCMMassRange", {0.10, 0.15}, "photon-pair mass integration range for pi0PCM"};
  Configurable<std::vector<double>> pi0PCMSideMassRange{"pi0PCMSideMassRange", {0.16, 0.24}, "photon-pair mass integration range outside outside pi0PCM region"};

  Configurable<bool> requireSingleCollisionPurity{"requireSingleCollisionPurity", true, "whether particle from single chosen MC-col associated to reco-col (else just type/kin match)"};

  // for histograms
  Configurable<int> nBinsZPv{"nBinsZPv", 100, "number zPv bins in histos for QA"};
  Configurable<int> nBinsZPvSmol{"nBinsZPvSmol", 28, "number zPv bins but smaller"};
  Configurable<int> nBinsMult{"nBinsMult", 200, "number multiplicity bins in histos for QA"};
  Configurable<int> nBinsMultSmol{"nBinsMultSmol", 20, "number multiplicity bins but smaller"};
  Configurable<int> nBinsOccupancy{"nBinsOccupancy", 2000, "number occupancy bins in histos for QA"};

  Configurable<int> nBinsPhi{"nBinsPhi", 72, "number phi bins"};
  Configurable<int> nBinsEta{"nBinsEta", 40, "number eta bins"};
  Configurable<int> nBinsMgg{"nBinsMgg", 160, "number mass-photon-pair bins"};

  Configurable<std::vector<double>> binsPtTrig{"binsPtTrig", {5, 10, 25, 50}, "correlation ptTrig bins"};
  Configurable<std::vector<double>> binsPtAssoc{"binsPtAssoc",
                                                {0.2, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10, 12.5, 15, 17.5, 20, 30, 40},
                                                "correlation ptAssoc bins"};
  Configurable<std::vector<double>> binsDPhi{"binsDPhi",
                                             {0.00 * DPHI_SCALE,
                                              0.04 * DPHI_SCALE, 0.08 * DPHI_SCALE, 0.11 * DPHI_SCALE, 0.14 * DPHI_SCALE,
                                              0.16 * DPHI_SCALE, 0.18 * DPHI_SCALE, 0.20 * DPHI_SCALE, 0.22 * DPHI_SCALE,
                                              0.23 * DPHI_SCALE, 0.24 * DPHI_SCALE, 0.25 * DPHI_SCALE, 0.26 * DPHI_SCALE, 0.27 * DPHI_SCALE, 0.28 * DPHI_SCALE,
                                              0.30 * DPHI_SCALE, 0.32 * DPHI_SCALE, 0.34 * DPHI_SCALE, 0.36 * DPHI_SCALE,
                                              0.39 * DPHI_SCALE, 0.42 * DPHI_SCALE, 0.46 * DPHI_SCALE, 0.50 * DPHI_SCALE,
                                              0.54 * DPHI_SCALE, 0.58 * DPHI_SCALE, 0.61 * DPHI_SCALE, 0.64 * DPHI_SCALE,
                                              0.66 * DPHI_SCALE, 0.68 * DPHI_SCALE, 0.70 * DPHI_SCALE, 0.72 * DPHI_SCALE,
                                              0.74 * DPHI_SCALE, 0.76 * DPHI_SCALE, 0.78 * DPHI_SCALE,
                                              0.80 * DPHI_SCALE, 0.82 * DPHI_SCALE, 0.84 * DPHI_SCALE, 0.86 * DPHI_SCALE,
                                              0.89 * DPHI_SCALE, 0.92 * DPHI_SCALE, 0.96 * DPHI_SCALE, 1.00 * DPHI_SCALE},
                                             "correlation bins DeltaPhi"};
  Configurable<std::vector<double>> binsDEta{"binsDEta",
                                             {0 / 32. * DETA_SCALE,
                                              1 / 32. * DETA_SCALE, 2 / 32. * DETA_SCALE, 3 / 32. * DETA_SCALE, 4 / 32. * DETA_SCALE,
                                              5 / 32. * DETA_SCALE, 6 / 32. * DETA_SCALE, 7 / 32. * DETA_SCALE, 8 / 32. * DETA_SCALE,
                                              9 / 32. * DETA_SCALE, 10 / 32. * DETA_SCALE, 11 / 32. * DETA_SCALE, 12 / 32. * DETA_SCALE, 13 / 32. * DETA_SCALE, 14 / 32. * DETA_SCALE,
                                              59 / 128. * DETA_SCALE, 62 / 128. * DETA_SCALE, 64 / 128. * DETA_SCALE, 66 / 128. * DETA_SCALE, 69 / 128. * DETA_SCALE, 18 / 32. * DETA_SCALE,
                                              19 / 32. * DETA_SCALE, 20 / 32. * DETA_SCALE, 21 / 32. * DETA_SCALE, 22 / 32. * DETA_SCALE, 23 / 32. * DETA_SCALE, 24 / 32. * DETA_SCALE,
                                              25 / 32. * DETA_SCALE, 26 / 32. * DETA_SCALE, 27 / 32. * DETA_SCALE, 28 / 32. * DETA_SCALE,
                                              29 / 32. * DETA_SCALE, 30 / 32. * DETA_SCALE, 31 / 32. * DETA_SCALE, 32 / 32. * DETA_SCALE},
                                             "correlation bins DeltaEta"};
  Configurable<std::vector<double>> binsZPv{"binsZPv",
                                            {-7, -5, -3, -1, 1, 3, 5, 7},
                                            "zPv mixing bins"};
  Configurable<std::vector<double>> binsMult{"binsMult",
                                             {-0.5, 9.5, 14.5, 19.5, 25.5, 32},
                                             "multiplicity mixing bins"};

  // configurables from other tasks

  double etaMax;

  // objects to hold histograms
  HistogramRegistry histos{"histogramRegistry", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  // ccdb calls
  Service<ccdb::BasicCCDBManager> ccdb;
  // for mc
  Service<framework::O2DatabasePDG> pdg;

  // random number generation
  static constexpr unsigned int SeedRandomEngine = 12345;
  std::mt19937 randomEngine{SeedRandomEngine};

  // partitions
  SliceCache cache;

  // prepare for per collision slicing
  Preslice<aod::JetTracks> perColTracks = aod::jtrack::collisionId;
  Preslice<aod::Triggers> perColTriggers = aod::corr_particle::jetCollisionId;
  Preslice<aod::Hadrons> perColHadrons = aod::corr_particle::jetCollisionId;
  Preslice<aod::Pipms> perColPipms = aod::corr_particle::jetCollisionId;
  Preslice<aod::PhotonPCMs> perColPhotonPCMs = aod::corr_particle::jetCollisionId;
  Preslice<aod::JetParticles> perColMcParticles = aod::jmcparticle::mcCollisionId;
  Preslice<aod::TriggerParticles> perColTriggerParticles = aod::corr_particle::jetMcCollisionId;

  // combinations binning
  // cumbersome, but still better than having extra configurable or figuring out how to init binningZPvMult later while declaring it here
  std::function<std::vector<double>(std::vector<double> const&, double)> prependValueToVector =
    [](std::vector<double> const& vec, double const value) {
      std::vector<double> resultVec = {value};
      resultVec.insert(resultVec.end(), vec.begin(), vec.end());
      return resultVec;
    };
  BinningZPvMult binningZPvMult{{prependValueToVector(binsZPv.value, VARIABLE_WIDTH), prependValueToVector(binsMult.value, VARIABLE_WIDTH)}, true};

  // declare analysis variables

  // efficiency histograms
  TH1D* h1PtInvEffTrigger;
  TH1D* h1PtInvEffHadron;
  TH1D* h1PtInvEffPipm;
  TH1D* h1PtInvEffPhotonPCM;

  // mixing trigger memory
  int nTriggersThisDataFrame;
  // organised as zPv- and mult-bin matrix of deques to save trigger info beyond single dataframe
  // extra bin for mult overflow
  // with ajusted zVtx (see triggerBinValuesZPv in init) and mult overflow -> all events accounted for
  // (possibly replace by some advanced derived data method and O2 event mixing in future?)
  std::vector<double> triggerBinValuesZPv;
  std::vector<double> triggerBinValuesMult;
  std::vector<std::vector<std::deque<float>>> savedTriggersZPvMultPt;
  std::vector<std::vector<std::deque<float>>> savedTriggersZPvMultPhi;
  std::vector<std::vector<std::deque<float>>> savedTriggersZPvMultEta;

  // functions ================================================================================================================================================================================

  // general (kenobi) /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // get histograms from ccdb
  // save efficiencies from ccdb in histogram registry
  void initCcdbHistograms()
  {
    // trigger
    h1PtInvEffTrigger = nullptr;
    if (doEffCorrectionTrigger) {
      h1PtInvEffTrigger = ccdb->getForTimeStamp<TH1D>(pathCcdbEff.value + "/trigger", noLaterThanCcdb.value);

      const double* effBinsTrigger = h1PtInvEffTrigger->GetXaxis()->GetXbins()->GetArray();
      const AxisSpec axisPtEffTrigger{std::vector<double>(effBinsTrigger, effBinsTrigger + h1PtInvEffTrigger->GetNbinsX() + 1), "#it{p}_{T}"};
      histos.add("usedEff/h1_pt_invEff_trigger_ccdb", "h1_pt_invEff_trigger_ccdb", kTH1D, {axisPtEffTrigger}, true);
      for (int iBin = 1; iBin <= h1PtInvEffTrigger->GetNbinsX(); iBin++) {
        histos.get<TH1>(HIST("usedEff/h1_pt_invEff_trigger_ccdb"))->SetBinContent(iBin, h1PtInvEffTrigger->GetBinContent(iBin));
        histos.get<TH1>(HIST("usedEff/h1_pt_invEff_trigger_ccdb"))->SetBinError(iBin, h1PtInvEffTrigger->GetBinError(iBin));
      }
    }
    // hadron
    h1PtInvEffHadron = nullptr;
    if (doEffCorrectionHadron) {
      h1PtInvEffHadron = ccdb->getForTimeStamp<TH1D>(pathCcdbEff.value + "/hadron", noLaterThanCcdb.value);

      const double* effBinsHadron = h1PtInvEffHadron->GetXaxis()->GetXbins()->GetArray();
      const AxisSpec axisPtEffHadron{std::vector<double>(effBinsHadron, effBinsHadron + h1PtInvEffHadron->GetNbinsX() + 1), "#it{p}_{T}"};
      histos.add("usedEff/h1_pt_invEff_hadron_ccdb", "h1_pt_invEff_hadron_ccdb", kTH1D, {axisPtEffHadron}, true);
      for (int iBin = 1; iBin <= h1PtInvEffHadron->GetNbinsX(); iBin++) {
        histos.get<TH1>(HIST("usedEff/h1_pt_invEff_hadron_ccdb"))->SetBinContent(iBin, h1PtInvEffHadron->GetBinContent(iBin));
        histos.get<TH1>(HIST("usedEff/h1_pt_invEff_hadron_ccdb"))->SetBinError(iBin, h1PtInvEffHadron->GetBinError(iBin));
      }
    }
    // pipm
    h1PtInvEffPipm = nullptr;
    if (doEffCorrectionPipm) {
      h1PtInvEffPipm = ccdb->getForTimeStamp<TH1D>(pathCcdbEff.value + "/pipm", noLaterThanCcdb.value);

      const double* effBinsPipm = h1PtInvEffPipm->GetXaxis()->GetXbins()->GetArray();
      const AxisSpec axisPtEffPipm{std::vector<double>(effBinsPipm, effBinsPipm + h1PtInvEffPipm->GetNbinsX() + 1), "#it{p}_{T}"};
      histos.add("usedEff/h1_pt_invEff_pipm_ccdb", "h1_pt_invEff_pipm_ccdb", kTH1D, {axisPtEffPipm}, true);
      for (int iBin = 1; iBin <= h1PtInvEffPipm->GetNbinsX(); iBin++) {
        histos.get<TH1>(HIST("usedEff/h1_pt_invEff_pipm_ccdb"))->SetBinContent(iBin, h1PtInvEffPipm->GetBinContent(iBin));
        histos.get<TH1>(HIST("usedEff/h1_pt_invEff_pipm_ccdb"))->SetBinError(iBin, h1PtInvEffPipm->GetBinError(iBin));
      }
    }
    // photonPCM
    h1PtInvEffPhotonPCM = nullptr;
    if (doEffCorrectionPhotonPCM) {
      h1PtInvEffPhotonPCM = ccdb->getForTimeStamp<TH1D>(pathCcdbEff.value + "/photonPCM", noLaterThanCcdb.value);

      const double* effBinsPhotonPCM = h1PtInvEffPhotonPCM->GetXaxis()->GetXbins()->GetArray();
      const AxisSpec axisPtEffPhotonPCM{std::vector<double>(effBinsPhotonPCM, effBinsPhotonPCM + h1PtInvEffPhotonPCM->GetNbinsX() + 1), "#it{p}_{T}"};
      histos.add("usedEff/h1_pt_invEff_photonPCM_ccdb", "h1_pt_invEff_photonPCM_ccdb", kTH1D, {axisPtEffPhotonPCM}, true);
      for (int iBin = 1; iBin <= h1PtInvEffPhotonPCM->GetNbinsX(); iBin++) {
        histos.get<TH1>(HIST("usedEff/h1_pt_invEff_photonPCM_ccdb"))->SetBinContent(iBin, h1PtInvEffPhotonPCM->GetBinContent(iBin));
        histos.get<TH1>(HIST("usedEff/h1_pt_invEff_photonPCM_ccdb"))->SetBinError(iBin, h1PtInvEffPhotonPCM->GetBinError(iBin));
      }
    }
  }

  // create histograms
  void initHistograms()
  {
    // define axes
    const AxisSpec axisN{1, 0., 1., "#it{N}_{something}"};
    const AxisSpec axisCategories{16, 0., 16., "categories"};

    const AxisSpec axisZPv{nBinsZPv, -10, 10, "#it{z}_{pv}"};
    const AxisSpec axisZPvSmol{nBinsZPvSmol, -7, 7, "#it{z}_{pv}"};
    const AxisSpec axisMult{nBinsMult + 1, -0.5, nBinsMult + 0.5, "multiplicity"};
    const AxisSpec axisMultSmol{nBinsMultSmol + 1, -0.5, nBinsMultSmol + 0.5, "multiplicity"};
    const AxisSpec axisOccupancy{nBinsOccupancy + 1, -0.5, nBinsOccupancy + 0.5, "occupancy"};

    const AxisSpec axisPhi{nBinsPhi, 0, constants::math::TwoPI, "#it{#varphi}"};
    const AxisSpec axisEta{nBinsEta, -etaMax, etaMax, "#it{#eta}"};
    const AxisSpec axisMgg{nBinsMgg, 0, 0.8, "#it{m}_{#gamma#gamma}"};

    const AxisSpec axisPtTrig{binsPtTrig, "#it{p}_{T}^{trig}"};
    const AxisSpec axisPtAssoc{binsPtAssoc, "#it{p}_{T}^{assoc}"};
    const AxisSpec axisDPhi{binsDPhi, "#Delta#it{#varphi}"};
    const AxisSpec axisDEta{binsDEta, "#Delta#it{#eta}"};
    const AxisSpec axisZPvBinning{binsZPv, "#it{z}_{pv} correlation binning"};
    const AxisSpec axisMultBinning{binsMult, "multiplicity correlation binning"};

    // reco info
    histos.add("reco/info/h1_nEvents", "h1_nEvents", kTH1D, {axisCategories});
    histos.get<TH1>(HIST("reco/info/h1_nEvents"))->GetXaxis()->SetBinLabel(1, "#it{N}_{ev}^{sel}");
    histos.get<TH1>(HIST("reco/info/h1_nEvents"))->GetXaxis()->SetBinLabel(2, "#it{N}_{ev}");
    histos.get<TH1>(HIST("reco/info/h1_nEvents"))->GetXaxis()->SetBinLabel(3, "#it{N}_{ev}^{trig}");

    histos.add("reco/info/h2_zPvMult", "h2_zPvMult", kTHnSparseD, {axisZPv, axisMult}, true);
    histos.add("reco/info/h1_occupancy", "h1_occupancy", kTH1D, {axisOccupancy}, true);

    // reco (correlation) analysis
    histos.add("reco/info/h2_zPvMult_trigEv", "h2_zPvMult_trigEv", kTHnSparseD, {axisZPv, axisMult}, true);
    histos.add("reco/info/h1_occupancy_trigEv", "h1_occupancy_trigEv", kTH1D, {axisOccupancy}, true);
    histos.add("reco/corr/h3_ptPhiEta_trig", "h3_ptPhiEta_trig", kTHnSparseD, {axisPtAssoc, axisPhi, axisEta}, true);

    // hadron
    histos.add("reco/plain/h3_ptPhiEta_hadron", "h3_ptPhiEta_hadron", kTHnSparseD, {axisPtAssoc, axisPhi, axisEta}, true);
    histos.add("reco/corr/h3_ptPhiEta_assoc_hadron", "h3_ptPhiEta_assoc_hadron", kTHnSparseD, {axisPtAssoc, axisPhi, axisEta}, true);
    histos.add("reco/corr/h6_corr_hadron", "h6_corr_hadron",
               kTHnSparseF, {axisDPhi, axisDEta, axisPtTrig, axisPtAssoc, axisZPvBinning, axisMultBinning}, true);
    histos.add("reco/corr/h6_mix_hadron", "h6_mix_hadron",
               kTHnSparseF, {axisDPhi, axisDEta, axisPtTrig, axisPtAssoc, axisZPvBinning, axisMultBinning}, true);
    // pipm
    histos.add("reco/plain/h3_ptPhiEta_pipm", "h3_ptPhiEta_pipm", kTHnSparseD, {axisPtAssoc, axisPhi, axisEta}, true);
    histos.add("reco/corr/h3_ptPhiEta_assoc_pipm", "h3_ptPhiEta_assoc_pipm", kTHnSparseD, {axisPtAssoc, axisPhi, axisEta}, true);
    histos.add("reco/corr/h6_corr_pipm", "h6_corr_pipm",
               kTHnSparseF, {axisDPhi, axisDEta, axisPtTrig, axisPtAssoc, axisZPvBinning, axisMultBinning}, true);
    histos.add("reco/corr/h6_mix_pipm", "h6_mix_pipm",
               kTHnSparseF, {axisDPhi, axisDEta, axisPtTrig, axisPtAssoc, axisZPvBinning, axisMultBinning}, true);
    // photonPCM
    histos.add("reco/plain/h3_ptPhiEta_photonPCM", "h3_ptPhiEta_photonPCM", kTHnSparseD, {axisPtAssoc, axisPhi, axisEta}, true);
    histos.add("reco/corr/h3_ptPhiEta_assoc_photonPCM", "h3_ptPhiEta_assoc_photonPCM", kTHnSparseD, {axisPtAssoc, axisPhi, axisEta}, true);
    histos.add("reco/corr/h6_corr_photonPCM", "h6_corr_photonPCM",
               kTHnSparseF, {axisDPhi, axisDEta, axisPtTrig, axisPtAssoc, axisZPvBinning, axisMultBinning}, true);
    histos.add("reco/corr/h6_mix_photonPCM", "h6_mix_photonPCM",
               kTHnSparseF, {axisDPhi, axisDEta, axisPtTrig, axisPtAssoc, axisZPvBinning, axisMultBinning}, true);
    // photonPCM pairs
    histos.add("reco/plain/h4_ptMggZPvMult_photonPCMPair", "h4_ptMggZPvMult_photonPCMPair", kTHnSparseD, {axisPtAssoc, axisMgg, axisZPvBinning, axisMultBinning}, true);
    histos.add("reco/plain/h3_ptPhiEta_pi0PCMPeak", "h3_ptPhiEta_pi0PCMPeak", kTHnSparseD, {axisPtAssoc, axisPhi, axisEta}, true);
    histos.add("reco/corr/h4_ptMggZPvMult_assoc_photonPCMPair", "h4_ptMggZPvMult_assoc_photonPCMPair", kTHnSparseD, {axisPtAssoc, axisMgg, axisZPvBinning, axisMultBinning}, true);
    histos.add("reco/corr/h3_ptPhiEta_assoc_pi0PCMPeak", "h3_ptPhiEta_assoc_pi0PCMPeak", kTHnSparseD, {axisPtAssoc, axisPhi, axisEta}, true);
    // peak (mgg)
    histos.add("reco/corr/h6_corr_pi0PCMPeak", "h6_corr_pi0PCMPeak",
               kTHnSparseF, {axisDPhi, axisDEta, axisPtTrig, axisPtAssoc, axisZPvBinning, axisMultBinning}, true);
    histos.add("reco/corr/h6_mix_pi0PCMPeak", "h6_mix_pi0PCMPeak",
               kTHnSparseF, {axisDPhi, axisDEta, axisPtTrig, axisPtAssoc, axisZPvBinning, axisMultBinning}, true);
    // side (mgg)
    histos.add("reco/corr/h6_corr_pi0PCMSide", "h6_corr_pi0PCMSide",
               kTHnSparseF, {axisDPhi, axisDEta, axisPtTrig, axisPtAssoc, axisZPvBinning, axisMultBinning}, true);
    histos.add("reco/corr/h6_mix_pi0PCMSide", "h6_mix_pi0PCMSide",
               kTHnSparseF, {axisDPhi, axisDEta, axisPtTrig, axisPtAssoc, axisZPvBinning, axisMultBinning}, true);
    // event mixing for photon pairs
    histos.add("reco/plain/h2_zPvMult_photonPCMPair_evMix", "h2_zPvMult_photonPCMPair_evMix", kTHnSparseD, {axisZPv, axisMult}, true);
    histos.add("reco/plain/h4_ptMggZPvMult_photonPCMPair_evMix", "h4_ptMggZPvMult_photonPCMPair_evMix", kTHnSparseD, {axisPtAssoc, axisMgg, axisZPvBinning, axisMultBinning}, true);
    histos.add("reco/plain/h3_ptPhiEta_pi0PCMPeak_evMix", "h3_ptPhiEta_pi0PCMPeak_evMix", kTHnSparseD, {axisPtAssoc, axisPhi, axisEta}, true);

    // mc info
    histos.add("mc/info/h1_nEvents_mcTrue", "h1_nEvents_mcTrue", kTH1D, {axisN});
    histos.add("mc/info/h1_nTriggerEvents_mcTrue", "h1_nTriggerEvents_mcTrue", kTH1D, {axisN});

    histos.add("mc/info/h1_zPv_mcTrue", "h1_zPv_mcTrue", kTH1D, {axisZPv}, true);
    histos.add("mc/info/h1_mult_mcTrue", "h1_mult_mcTrue", kTH1D, {axisMult}, true);

    // reco and true collision correlations
    for (auto const& collision_type : {"true", "true_reco"}) {
      histos.add(std::format("mc/{}/corr/h3_ptPhiEta_trig", collision_type).data(), "h3_ptPhiEta_trig", kTHnSparseD, {axisPtAssoc, axisPhi, axisEta}, true);
      // hadron
      histos.add(std::format("mc/{}/corr/h3_ptPhiEta_assoc_hadron", collision_type).data(), "h3_ptPhiEta_assoc_hadron",
                 kTHnSparseD, {axisPtAssoc, axisPhi, axisEta}, true);
      histos.add(std::format("mc/{}/corr/h4_corr_hadron", collision_type).data(), "h4_corr_hadron",
                 kTHnSparseD, {axisDPhi, axisDEta, axisPtTrig, axisPtAssoc}, true);
      // pipm
      histos.add(std::format("mc/{}/corr/h3_ptPhiEta_assoc_pipm", collision_type).data(), "h3_ptPhiEta_assoc_pipm",
                 kTHnSparseD, {axisPtAssoc, axisPhi, axisEta}, true);
      histos.add(std::format("mc/{}/corr/h4_corr_pipm", collision_type).data(), "h4_corr_pipm",
                 kTHnSparseD, {axisDPhi, axisDEta, axisPtTrig, axisPtAssoc}, true);
      // photon
      histos.add(std::format("mc/{}/corr/h3_ptPhiEta_assoc_photon", collision_type).data(), "h3_ptPhiEta_assoc_photon",
                 kTHnSparseD, {axisPtAssoc, axisPhi, axisEta}, true);
      histos.add(std::format("mc/{}/corr/h4_corr_photon", collision_type).data(), "h4_corr_photon",
                 kTHnSparseD, {axisDPhi, axisDEta, axisPtTrig, axisPtAssoc}, true);
      // pi0
      histos.add(std::format("mc/{}/corr/h3_ptPhiEta_assoc_pi0", collision_type).data(), "h3_ptPhiEta_assoc_pi0",
                 kTHnSparseD, {axisPtAssoc, axisPhi, axisEta}, true);
      histos.add(std::format("mc/{}/corr/h4_corr_pi0", collision_type).data(), "h4_corr_pi0",
                 kTHnSparseD, {axisDPhi, axisDEta, axisPtTrig, axisPtAssoc}, true);
    }

    // mc efficiency/purity
    std::function<void(std::string)> add_effHists =
      [&](std::string name_id) {
        histos.add(std::format("mc/eff/h3_ptPhiEta_{}", name_id).data(), "h3_ptPhiEta_mcReco_hadron",
                   kTHnSparseD, {axisPtAssoc, axisPhi, axisEta}, true);
        histos.add(std::format("mc/eff/h3_ptZPvMult_{}", name_id).data(), "h3_ptZPvMult_mcReco_hadron",
                   kTHnSparseD, {axisPtAssoc, axisZPvSmol, axisMultSmol}, true);
      };
    // mc tracks
    add_effHists("mcReco_hadron");
    add_effHists("mcReco_hasCorrectMc_hadron");
    add_effHists("mcTrue_hadron");
    add_effHists("mcTrue_recoCol_hadron");
    // mc pipm PID
    add_effHists("mcReco_pipm");
    add_effHists("mcReco_hasCorrectMc_pipm");
    add_effHists("mcTrue_pipm");
    add_effHists("mcTrue_recoCol_pipm");
    // mc photonPCM
    add_effHists("mcReco_photonPCM");
    add_effHists("mcReco_hasCorrectMc_photonPCM");
    add_effHists("mcTrue_photon");
    add_effHists("mcTrue_recoCol_photon");
    // mc pi0PCM
    add_effHists("mcReco_pi0PCM");
    add_effHists("mcReco_hasCorrectMc_pi0PCM");
    add_effHists("mcTrue_pi0");
    add_effHists("mcTrue_recoCol_pi0");

    // test of the test while testing another test. featuring a test
    histos.add("test/h2_mult_comp", "h2_mult_comp", kTH2D, {axisMult, axisMult}, true);
    histos.add("test/h2_tracks_zPvMultDep", "h2_tracks_zPvMultDep", kTH2D, {axisZPv, axisMult}, true);
    histos.add("test/h2_globalTracks_zPvMultDep", "h2_globalTracks_zPvMultDep", kTH2D, {axisZPv, axisMult}, true);
  }

  // selections ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // checks if mcParticle is charged
  template <typename T_mcParticle>
  bool checkChargedMc(T_mcParticle const& mcParticle)
  {
    auto const pdgParticle = pdg->GetParticle(mcParticle.pdgCode());
    if (!pdgParticle || pdgParticle->Charge() == 0)
      return false;
    return true;
  }
  // checks if mcParticle should be detected (physicalPrimary, |eta|)
  template <typename T_mcParticle>
  bool checkPrimaryEtaMc(T_mcParticle const& mcParticle)
  {
    if (!mcParticle.isPhysicalPrimary())
      return false;
    if (std::abs(mcParticle.eta()) > etaMax)
      return false;
    return true;
  }
  // checks if mcParticle should be detected as primary track (physicalPrimary, charge, |eta|)
  template <typename T_mcParticle>
  bool checkPrimaryTrackMc(T_mcParticle const& mcParticle)
  {
    if (!checkPrimaryEtaMc(mcParticle))
      return false;
    if (!checkChargedMc(mcParticle))
      return false;
    return true;
  }
  // checks if mcParticle should be detected as 'primary' pi0->gg (|eta| not checked)
  template <typename T_mcParticle>
  bool checkPi0ToGG(T_mcParticle const& mcParticle)
  {
    if (mcParticle.pdgCode() != PDG_t::kPi0)
      return false;
    // identify primary pi0 (account for 0 daughters for some reason)
    if (mcParticle.template daughters_as<aod::JetParticles>().size() == 0)
      return false;
    for (auto const& pi0_daughter : mcParticle.template daughters_as<aod::JetParticles>()) {
      if (!pi0_daughter.isPhysicalPrimary())
        return false;
    }
    // select pi0 -> gg
    constexpr int NDaughtersPi0ToGG = 2;
    if (mcParticle.template daughters_as<aod::JetParticles>().size() != NDaughtersPi0ToGG)
      return false;
    return true;
  }

  // analysis helpers /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  template <typename T_h1>
  double getH1ValueAt(T_h1 const* const h1, double const value)
  {
    return h1->GetBinContent(h1->FindFixBin(value));
  }
  // efficiency helpers
  // define enum class for particle type
  enum class ParticleType { Trigger,
                            Hadron,
                            Pipm,
                            PhotonPCM };
  // efficiency function
  template <ParticleType T>
  double getInvEff(double const value)
  {
    if constexpr (T == ParticleType::Trigger) {
      return doEffCorrectionTrigger ? getH1ValueAt(h1PtInvEffTrigger, value) : 1;
    } else if constexpr (T == ParticleType::Hadron) {
      return doEffCorrectionHadron ? getH1ValueAt(h1PtInvEffHadron, value) : 1;
    } else if constexpr (T == ParticleType::Pipm) {
      return doEffCorrectionPipm ? getH1ValueAt(h1PtInvEffPipm, value) : 1;
    } else if constexpr (T == ParticleType::PhotonPCM) {
      return doEffCorrectionPhotonPCM ? getH1ValueAt(h1PtInvEffPhotonPCM, value) : 1;
    } else {
      return 1;
    }
  }

  // performs 'phi1 - phi2' and pushes it into the interval [-pi/2, 3pi/2]
  inline double getDeltaPhi(double const phi1, double const phi2)
  {
    return RecoDecay::constrainAngle(phi1 - phi2, -1 * constants::math::PIHalf);
  }

  // finds bin that value belongs to (assumes ordered bins) (starts at 0; includes underflow (return -1) and overlflow (return bins.size() - 1))
  // should be faster than some std binary search due to small number of bins (zPv, mult)
  int findIntervalBin(double value, const std::vector<double>& bins)
  {
    const int n = bins.size() - 1;
    if (value < bins[0])
      return -1; // underflow
    for (int i_bin = 0; i_bin < n; i_bin++)
      if (value < bins[i_bin + 1])
        return i_bin;
    return n; // overflow
  }

  // checks that two values belong to the same category (assumes ordered bins)
  // returns -1 for negative result (also for under/overflow values) and bin number (starting at 0) otherwise
  int checkSameBin(double const value1, double const value2, std::vector<double> const& bins)
  {
    // reject underflow
    if (value1 < bins[0])
      return -1;
    // loop over bins
    const int n = bins.size() - 1;
    for (int i_bin = 0; i_bin < n; i_bin++) {
      if (value1 < bins[i_bin + 1]) {
        if (value2 < bins[i_bin + 1] && value2 >= bins[i_bin]) {
          return i_bin;
        }
        return -1;
      }
    }
    // reject overflow
    return -1;
  }

  // analysis /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // generalised correlation functions
  // per collision

  // plain info
  template <typename T_collision, typename T_associatedThisEvent,
            typename T_funcPlain>
  void corrProcessPlain(T_collision const& collision, T_associatedThisEvent const& associatedThisEvent,
                        T_funcPlain&& funcPlain)
  {
    // normal spectra (per event - not per trigger)
    for (auto const& associated : associatedThisEvent) {
      funcPlain(collision, associated);
    }
  }

  // correlation
  template <typename T_collision, typename T_triggersThisEvent, typename T_associatedThisEvent,
            typename T_funcCorrelation>
  void corrProcessCorrelation(T_collision const& collision, T_triggersThisEvent const& triggersThisEvent, T_associatedThisEvent const& associatedThisEvent,
                              T_funcCorrelation&& funcCorrelation)
  {
    // correlation combinations
    for (auto const& [trigger, associated] : soa::combinations(soa::CombinationsFullIndexPolicy(triggersThisEvent, associatedThisEvent))) {
      funcCorrelation(collision, trigger, associated);
    }
  }

  // mixing
  template <typename T_collision, typename T_associatedThisEvent,
            typename T_funcMixing>
  void corrProcessMixing(T_collision const& collision, T_associatedThisEvent const& associatedThisEvent,
                         T_funcMixing&& funcMixing,
                         size_t const nTriggerMixing)
  {
    // skip if event does not contain valid trigger
    if (doTrigEvMixing && !collision.trigEv())
      return;

    // mixing loops (more efficient than O2 mixing (for now))
    // prepare zPv-mult binned saved triggers
    const int iBinCorrZPv = findIntervalBin(collision.posZ(), triggerBinValuesZPv);
    const int iBinCorrMult = findIntervalBin(collision.multNTracksGlobal(), triggerBinValuesMult);
    auto const& savedTriggersPt = savedTriggersZPvMultPt[iBinCorrZPv][iBinCorrMult];
    auto const& savedTriggersPhi = savedTriggersZPvMultPhi[iBinCorrZPv][iBinCorrMult];
    auto const& savedTriggersEta = savedTriggersZPvMultEta[iBinCorrZPv][iBinCorrMult];
    // number of triggers
    const int mixUpToTriggerN = std::min(savedTriggersPt.size(), nTriggerMixing + nTriggersThisDataFrame);
    const float perTriggerWeight = 1. / (mixUpToTriggerN - nTriggersThisDataFrame); // mixUpToTriggerN <= nTriggersThisDataFrame not problematic since no loop then
    // mixing loops
    for (int i_mixingTrigger = nTriggersThisDataFrame; i_mixingTrigger < mixUpToTriggerN; i_mixingTrigger++) {
      for (auto const& associated : associatedThisEvent) {
        funcMixing(collision, savedTriggersPt[i_mixingTrigger], savedTriggersPhi[i_mixingTrigger], savedTriggersEta[i_mixingTrigger], associated, perTriggerWeight);
      }
    }
  }

  void init(InitContext& initContext)
  {
    // analysis info
    ccdb->setURL(urlCcdb.value);
    // enabling object caching (otherwise each call goes to CCDB server)
    ccdb->setCaching(true);
    // ccdb->setLocalObjectValidityChecking();
    // not later than now, will be replaced by the value of train creation (avoids replacing objects while a train is running)
    ccdb->setCreatedNotAfter(noLaterThanCcdb.value);

    // init analysis variables

    // get variabels from other tasks
    getTaskOptionValue(initContext, "correlation-table-producer", "etaMax", etaMax, false);

    // mixing trigger memory
    triggerBinValuesZPv = binsZPv;
    triggerBinValuesMult = binsMult;
    // prevent rounding errors in bin finding (multiplicity accounted for by it going to 0 and already considering overflow separately)
    triggerBinValuesZPv.front() *= 1.0001;
    triggerBinValuesZPv.back() *= 1.0001;
    // init correct size of zPv-mult matrix
    savedTriggersZPvMultPt.resize(binsZPv.value.size() - 1);
    savedTriggersZPvMultPhi.resize(binsZPv.value.size() - 1);
    savedTriggersZPvMultEta.resize(binsZPv.value.size() - 1);
    for (size_t i_zPv = 0; i_zPv < binsZPv.value.size() - 1; i_zPv++) {
      savedTriggersZPvMultPt[i_zPv].resize(binsMult.value.size());
      savedTriggersZPvMultPhi[i_zPv].resize(binsMult.value.size());
      savedTriggersZPvMultEta[i_zPv].resize(binsMult.value.size());
    }

    // histograms from ccdb
    initCcdbHistograms();

    // create analysis histograms
    initHistograms();
  }

  // reconstructed ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void processInfo(CorrCollision const& collision)
  {
    // trigger events
    if (collision.trigEv()) {
      histos.fill(HIST("reco/info/h1_nEvents"), 2.5);
    }

    // event selection
    histos.fill(HIST("reco/info/h1_nEvents"), 1.5);
    if (!collision.selEv())
      return;
    histos.fill(HIST("reco/info/h1_nEvents"), 0.5);

    // QA
    histos.fill(HIST("reco/info/h2_zPvMult"), collision.posZ(), collision.multNTracksGlobal());
    histos.fill(HIST("reco/info/h1_occupancy"), collision.trackOccupancyInTimeRange());
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processInfo, "process general info on collisions and tracks for analysis and qa", false);

  void processCorrFirst(CorrCollisions const& collisions, aod::Triggers const& triggers)
  {
    // do at beginning of each data frame (before other correlation process functions)
    // (PROCESS_SWITCH of this process has to be declared first)

    // set trigger counter
    nTriggersThisDataFrame = triggers.size();

    for (auto const& collision : collisions) {
      // event selection
      if (!collision.selEv())
        continue;

      // group collision
      auto const triggersThisEvent = triggers.sliceBy(perColTriggers, collision.globalIndex());

      // trigger loop
      for (auto const& trigger : triggersThisEvent) {
        // trigger info
        histos.fill(HIST("reco/corr/h3_ptPhiEta_trig"), trigger.pt(), trigger.phi(), trigger.eta(),
                    getInvEff<ParticleType::Trigger>(trigger.pt()));

        // save triggers for mixing
        const int iBinCorrZPv = findIntervalBin(collision.posZ(), triggerBinValuesZPv);
        const int iBinCorrMult = findIntervalBin(collision.multNTracksGlobal(), triggerBinValuesMult);
        // special cases (floating point precision errors, mult overflow) should be taken care of by triggerBinValuesZPv and triggerBinValuesMult
        savedTriggersZPvMultPt[iBinCorrZPv][iBinCorrMult].push_front(trigger.pt());
        savedTriggersZPvMultPhi[iBinCorrZPv][iBinCorrMult].push_front(trigger.phi());
        savedTriggersZPvMultEta[iBinCorrZPv][iBinCorrMult].push_front(trigger.eta());
        if (static_cast<int>(savedTriggersZPvMultPt[iBinCorrZPv][iBinCorrMult].size()) > nTriggerSavedForMixing) {
          savedTriggersZPvMultPt[iBinCorrZPv][iBinCorrMult].pop_back();
          savedTriggersZPvMultPhi[iBinCorrZPv][iBinCorrMult].pop_back();
          savedTriggersZPvMultEta[iBinCorrZPv][iBinCorrMult].pop_back();
        }
      }

      // trigger event info
      if (collision.trigEv()) {
        histos.fill(HIST("reco/info/h2_zPvMult_trigEv"), collision.posZ(), collision.multNTracksGlobal());
        histos.fill(HIST("reco/info/h1_occupancy_trigEv"), collision.trackOccupancyInTimeRange());
      }
    }
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processCorrFirst, "process to gather info before correlation processes", false);

  void processCorrHadron(CorrCollision const& collision, aod::Triggers const& triggers, aod::Hadrons const& hadrons)
  {
    // event selection
    if (!collision.selEv())
      return;

    auto const funcPlain = [this]([[maybe_unused]] auto const& collision, auto const& associated) {
      histos.fill(HIST("reco/plain/h3_ptPhiEta_hadron"),
                  associated.pt(), associated.phi(), associated.eta(),
                  getInvEff<ParticleType::Hadron>(associated.pt()));
    };
    corrProcessPlain(collision, hadrons, funcPlain);

    auto const funcCorrelation = [this](auto const& collision, auto const& trigger, auto const& associated) {
      // exclude self correlation
      if (trigger.jetTrackId() == associated.jetTrackId())
        return;

      histos.fill(HIST("reco/corr/h3_ptPhiEta_assoc_hadron"),
                  associated.pt(), associated.phi(), associated.eta(),
                  getInvEff<ParticleType::Trigger>(trigger.pt()) * getInvEff<ParticleType::Hadron>(associated.pt()));
      histos.fill(HIST("reco/corr/h6_corr_hadron"),
                  getDeltaPhi(trigger.phi(), associated.phi()),
                  trigger.eta() - associated.eta(),
                  trigger.pt(), associated.pt(), collision.posZ(), collision.multNTracksGlobal(),
                  getInvEff<ParticleType::Trigger>(trigger.pt()) * getInvEff<ParticleType::Hadron>(associated.pt()));
    };
    corrProcessCorrelation(collision, triggers, hadrons, funcCorrelation);

    auto const funcMixing = [this](auto const& collision,
                                   float const mixingTriggerPt, float const mixingTriggerPhi, float const mixingTriggerEta, auto const& associated, auto const perTriggerWeight) {
      histos.fill(HIST("reco/corr/h6_mix_hadron"),
                  getDeltaPhi(mixingTriggerPhi, associated.phi()),
                  mixingTriggerEta - associated.eta(),
                  mixingTriggerPt, associated.pt(), collision.posZ(), collision.multNTracksGlobal(),
                  perTriggerWeight * getInvEff<ParticleType::Trigger>(mixingTriggerPt) * getInvEff<ParticleType::Hadron>(associated.pt()));
    };
    corrProcessMixing(collision, hadrons, funcMixing, nTriggerMixingHadron);
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processCorrHadron, "process standard correlation for associated hardons", false);

  void processCorrPipm(CorrCollision const& collision, aod::Triggers const& triggers, aod::Pipms const& pipms)
  {
    // event selection
    if (!collision.selEv())
      return;

    auto const funcPlain = [this]([[maybe_unused]] auto const& collision, auto const& associated) {
      histos.fill(HIST("reco/plain/h3_ptPhiEta_pipm"),
                  associated.pt(), associated.phi(), associated.eta(),
                  getInvEff<ParticleType::Pipm>(associated.pt()));
    };
    corrProcessPlain(collision, pipms, funcPlain);

    auto const funcCorrelation = [this](auto const& collision, auto const& trigger, auto const& associated) {
      // exclude self correlation
      if (trigger.jetTrackId() == associated.jetTrackId())
        return;

      histos.fill(HIST("reco/corr/h3_ptPhiEta_assoc_pipm"),
                  associated.pt(), associated.phi(), associated.eta(),
                  getInvEff<ParticleType::Trigger>(trigger.pt()) * getInvEff<ParticleType::Pipm>(associated.pt()));
      histos.fill(HIST("reco/corr/h6_corr_pipm"),
                  getDeltaPhi(trigger.phi(), associated.phi()),
                  trigger.eta() - associated.eta(),
                  trigger.pt(), associated.pt(), collision.posZ(), collision.multNTracksGlobal(),
                  getInvEff<ParticleType::Trigger>(trigger.pt()) * getInvEff<ParticleType::Pipm>(associated.pt()));
    };
    corrProcessCorrelation(collision, triggers, pipms, funcCorrelation);

    auto const funcMixing = [this](auto const& collision,
                                   float const mixingTriggerPt, float const mixingTriggerPhi, float const mixingTriggerEta, auto const& associated, auto const perTriggerWeight) {
      histos.fill(HIST("reco/corr/h6_mix_pipm"),
                  getDeltaPhi(mixingTriggerPhi, associated.phi()),
                  mixingTriggerEta - associated.eta(),
                  mixingTriggerPt, associated.pt(), collision.posZ(), collision.multNTracksGlobal(),
                  perTriggerWeight * getInvEff<ParticleType::Trigger>(mixingTriggerPt) * getInvEff<ParticleType::Pipm>(associated.pt()));
    };
    corrProcessMixing(collision, pipms, funcMixing, nTriggerMixingPipm);
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processCorrPipm, "process standard correlation for associated pipm", false);

  void processCorrPhotonPCM(CorrCollision const& collision, aod::Triggers const& triggers, aod::PhotonPCMs const& photonPCMs)
  {
    // event selection
    if (!collision.selEv())
      return;

    auto const funcPlain = [this]([[maybe_unused]] auto const& collision, auto const& associated) {
      histos.fill(HIST("reco/plain/h3_ptPhiEta_photonPCM"),
                  associated.pt(), associated.phi(), associated.eta(),
                  getInvEff<ParticleType::PhotonPCM>(associated.pt()));
    };
    corrProcessPlain(collision, photonPCMs, funcPlain);

    auto const funcCorrelation = [this](auto const& collision, auto const& trigger, auto const& associated) {
      // exclude self correlation
      if (trigger.jetTrackId() == associated.posTrackId() || trigger.jetTrackId() == associated.negTrackId())
        return;

      histos.fill(HIST("reco/corr/h3_ptPhiEta_assoc_photonPCM"),
                  associated.pt(), associated.phi(), associated.eta(),
                  getInvEff<ParticleType::Trigger>(trigger.pt()) * getInvEff<ParticleType::PhotonPCM>(associated.pt()));
      histos.fill(HIST("reco/corr/h6_corr_photonPCM"),
                  getDeltaPhi(trigger.phi(), associated.phi()),
                  trigger.eta() - associated.eta(),
                  trigger.pt(), associated.pt(), collision.posZ(), collision.multNTracksGlobal(),
                  getInvEff<ParticleType::Trigger>(trigger.pt()) * getInvEff<ParticleType::PhotonPCM>(associated.pt()));
    };
    corrProcessCorrelation(collision, triggers, photonPCMs, funcCorrelation);

    auto const funcMixing = [this](auto const& collision,
                                   float const mixingTriggerPt, float const mixingTriggerPhi, float const mixingTriggerEta, auto const& associated, auto const perTriggerWeight) {
      histos.fill(HIST("reco/corr/h6_mix_photonPCM"),
                  getDeltaPhi(mixingTriggerPhi, associated.phi()),
                  mixingTriggerEta - associated.eta(),
                  mixingTriggerPt, associated.pt(), collision.posZ(), collision.multNTracksGlobal(),
                  perTriggerWeight * getInvEff<ParticleType::Trigger>(mixingTriggerPt) * getInvEff<ParticleType::PhotonPCM>(associated.pt()));
    };
    corrProcessMixing(collision, photonPCMs, funcMixing, nTriggerMixingPhotonPCM);
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processCorrPhotonPCM, "process standard correlation for associated photonPCM", false);

  void processCorrPi0PCM(CorrCollision const& collision, aod::Triggers const& triggers, aod::PhotonPCMPairs const& photonPCMPairs)
  {
    // event selection
    if (!collision.selEv())
      return;

    auto const funcPlain = [this](auto const& collision, auto const& associated) {
      histos.fill(HIST("reco/plain/h4_ptMggZPvMult_photonPCMPair"), associated.pt(), associated.mgg(), collision.posZ(), collision.multNTracksGlobal());
      // pi0 mass range
      if (associated.mgg() > pi0PCMMassRange.value[0] && associated.mgg() < pi0PCMMassRange.value[1]) {
        histos.fill(HIST("reco/plain/h3_ptPhiEta_pi0PCMPeak"), associated.pt(), associated.phi(), associated.eta());
      }
    };
    corrProcessPlain(collision, photonPCMPairs, funcPlain);

    auto const funcCorrelation = [this](auto const& collision, auto const& trigger, auto const& associated) {
      // exclude self correlation
      if (trigger.jetTrackId() == associated.posTrack1Id() || trigger.jetTrackId() == associated.negTrack1Id() ||
          trigger.jetTrackId() == associated.negTrack2Id() || trigger.jetTrackId() == associated.posTrack2Id())
        return;

      histos.fill(HIST("reco/corr/h4_ptMggZPvMult_assoc_photonPCMPair"),
                  associated.pt(), associated.mgg(), collision.posZ(), collision.multNTracksGlobal(),
                  getInvEff<ParticleType::Trigger>(trigger.pt()));

      if (associated.mgg() > pi0PCMMassRange.value[0] && associated.mgg() < pi0PCMMassRange.value[1]) {
        // pi0 mass range
        histos.fill(HIST("reco/corr/h3_ptPhiEta_assoc_pi0PCMPeak"),
                    associated.pt(), associated.phi(), associated.eta(),
                    getInvEff<ParticleType::Trigger>(trigger.pt()));
        histos.fill(HIST("reco/corr/h6_corr_pi0PCMPeak"),
                    getDeltaPhi(trigger.phi(), associated.phi()),
                    trigger.eta() - associated.eta(),
                    trigger.pt(), associated.pt(), collision.posZ(), collision.multNTracksGlobal(),
                    getInvEff<ParticleType::Trigger>(trigger.pt()));
      } else if (associated.mgg() > pi0PCMSideMassRange.value[0] && associated.mgg() < pi0PCMSideMassRange.value[1]) {
        // pi0 mass side range
        histos.fill(HIST("reco/corr/h6_corr_pi0PCMSide"),
                    getDeltaPhi(trigger.phi(), associated.phi()),
                    trigger.eta() - associated.eta(),
                    trigger.pt(), associated.pt(), collision.posZ(), collision.multNTracksGlobal(),
                    getInvEff<ParticleType::Trigger>(trigger.pt()));
      }
    };
    corrProcessCorrelation(collision, triggers, photonPCMPairs, funcCorrelation);

    auto const funcMixing = [this](auto const& collision,
                                   float const mixingTriggerPt, float const mixingTriggerPhi, float const mixingTriggerEta, auto const& associated, auto const perTriggerWeight) {
      if (associated.mgg() > pi0PCMMassRange.value[0] && associated.mgg() < pi0PCMMassRange.value[1]) {
        // pi0 mass range
        histos.fill(HIST("reco/corr/h6_mix_pi0PCMPeak"),
                    getDeltaPhi(mixingTriggerPhi, associated.phi()),
                    mixingTriggerEta - associated.eta(),
                    mixingTriggerPt, associated.pt(), collision.posZ(), collision.multNTracksGlobal(),
                    perTriggerWeight * getInvEff<ParticleType::Trigger>(mixingTriggerPt));
      } else if (associated.mgg() > pi0PCMSideMassRange.value[0] && associated.mgg() < pi0PCMSideMassRange.value[1]) {
        // pi0 mass side range
        histos.fill(HIST("reco/corr/h6_mix_pi0PCMSide"),
                    getDeltaPhi(mixingTriggerPhi, associated.phi()),
                    mixingTriggerEta - associated.eta(),
                    mixingTriggerPt, associated.pt(), collision.posZ(), collision.multNTracksGlobal(),
                    perTriggerWeight * getInvEff<ParticleType::Trigger>(mixingTriggerPt));
      }
    };
    corrProcessMixing(collision, photonPCMPairs, funcMixing, nTriggerMixingPi0PCM);
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processCorrPi0PCM, "process standard correlation for associated pi0PCM", false);

  void processCorrPi0PCMMix(CorrCollisions const& collisions, aod::PhotonPCMs const& photonPCMs)
  {
    auto photonPCMsTuple = std::make_tuple(photonPCMs);
    SameKindPair<CorrCollisions, aod::PhotonPCMs, BinningZPvMult> pairs{binningZPvMult, nNeighboursMixingPi0PCMPair, -1, collisions, photonPCMsTuple, &cache};

    // mixed events
    for (auto pair = pairs.begin(); pair != pairs.end(); pair++) {
      auto const& [collision1, photonPCMs1, collision2, photonPCMs2] = *pair;

      // // check that current und mixing-trigger event are from the same zPv/mult bins
      // if (checkSameBin(collision1.posZ(), collision2.posZ(), binsZPv) == -1) {
      //   std::printf("ERROR: zPv bins do not match\n"); continue;
      // }
      // if (checkSameBin(collision1.multNTracksGlobal(), collision2.multNTracksGlobal(), binsMult) == -1) {
      //   std::printf("ERROR: multiplicity bins do not match\n"); continue;
      // }

      // event selection
      if (!collision1.selEv())
        continue;
      if (!collision2.selEv())
        continue;

      // event info
      histos.fill(HIST("reco/plain/h2_zPvMult_photonPCMPair_evMix"), collision1.posZ(), collision1.multNTracksGlobal());

      // mixing loop
      for (auto const& [photonPCM1, photonPCM2] : soa::combinations(soa::CombinationsFullIndexPolicy(photonPCMs1, photonPCMs2))) {
        ROOT::Math::PtEtaPhiMVector const p4photonPCM1(photonPCM1.pt(), photonPCM1.eta(), photonPCM1.phi(), 0.);
        ROOT::Math::PtEtaPhiMVector const p4photonPCM2(photonPCM2.pt(), photonPCM2.eta(), photonPCM2.phi(), 0.);
        ROOT::Math::PtEtaPhiMVector const p4photonPCMPair = p4photonPCM1 + p4photonPCM2;

        // plain
        histos.fill(HIST("reco/plain/h4_ptMggZPvMult_photonPCMPair_evMix"), p4photonPCMPair.pt(), p4photonPCMPair.M(), collision1.posZ(), collision1.multNTracksGlobal());
        // pi0 mass range
        if (p4photonPCMPair.M() > pi0PCMMassRange.value[0] && p4photonPCMPair.M() < pi0PCMMassRange.value[1]) {
          histos.fill(HIST("reco/plain/h3_ptPhiEta_pi0PCMPeak_evMix"), p4photonPCMPair.pt(), p4photonPCMPair.phi() + constants::math::PI, p4photonPCMPair.eta());
        }
      }
    }
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processCorrPi0PCMMix, "process gamma-gamma mixing for photonPCM", false);

  // mc ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void processMcInfo(CorrMcCollision const& mcCollision)
  {
    // event counter
    histos.fill(HIST("mc/info/h1_nEvents_mcTrue"), 0.5);
    // trigger events
    if (mcCollision.trigEv()) {
      histos.fill(HIST("mc/info/h1_nTriggerEvents_mcTrue"), 0.5);
    }

    // QA
    histos.fill(HIST("mc/info/h1_zPv_mcTrue"), mcCollision.posZ());
    histos.fill(HIST("mc/info/h1_mult_mcTrue"), mcCollision.multMCNParticlesEta08());
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processMcInfo, "process general info on mc collisions and tracks for analysis and qa", false);

  void processMcTrueCorr(CorrMcCollision const&, aod::TriggerParticles const& triggerParticles, aod::JetParticles const& mcParticles)
  {
    // trigger pairing loop
    for (auto const& trigger : triggerParticles) {
      // trigger info
      histos.fill(HIST("mc/true/corr/h3_ptPhiEta_trig"), trigger.pt(), trigger.phi(), trigger.eta());

      // hadrons (tracks) and pipm
      for (auto const& associated : mcParticles) {
        // exclude self correlation
        if (trigger.jetMcParticleId() == associated.globalIndex())
          continue;

        // standard particles (marked physical primary)
        if (checkPrimaryEtaMc(associated)) {
          // charged primary ('hadron') selection
          if (checkChargedMc(associated)) {
            histos.fill(HIST("mc/true/corr/h3_ptPhiEta_assoc_hadron"), associated.pt(), associated.phi(), associated.eta());
            histos.fill(HIST("mc/true/corr/h4_corr_hadron"),
                        getDeltaPhi(trigger.phi(), associated.phi()),
                        trigger.eta() - associated.eta(),
                        trigger.pt(), associated.pt());
          }

          // pipm selection
          if (std::abs(associated.pdgCode()) == PDG_t::kPiPlus) {
            histos.fill(HIST("mc/true/corr/h3_ptPhiEta_assoc_pipm"), associated.pt(), associated.phi(), associated.eta());
            histos.fill(HIST("mc/true/corr/h4_corr_pipm"),
                        getDeltaPhi(trigger.phi(), associated.phi()),
                        trigger.eta() - associated.eta(),
                        trigger.pt(), associated.pt());
          }

          // photon selection
          if (associated.pdgCode() == PDG_t::kGamma) {
            histos.fill(HIST("mc/true/corr/h3_ptPhiEta_assoc_photon"), associated.pt(), associated.phi(), associated.eta());
            histos.fill(HIST("mc/true/corr/h4_corr_photon"),
                        getDeltaPhi(trigger.phi(), associated.phi()),
                        trigger.eta() - associated.eta(),
                        trigger.pt(), associated.pt());
          }
        }

        // decaying particles (not marked physical primary)
        if ((std::abs(associated.eta()) < etaMax)) {
          // pi0 selection
          if (checkPi0ToGG(associated)) {
            histos.fill(HIST("mc/true/corr/h3_ptPhiEta_assoc_pi0"), associated.pt(), associated.phi(), associated.eta());
            histos.fill(HIST("mc/true/corr/h4_corr_pi0"),
                        getDeltaPhi(trigger.phi(), associated.phi()),
                        trigger.eta() - associated.eta(),
                        trigger.pt(), associated.pt());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processMcTrueCorr, "process mc-true (all collisions) correlation for multiple associated particles", false);

  void processMcTrueRecoColCorr(CorrMcDCollision const& collision, aod::JetMcCollisions const&, aod::TriggerParticles const& triggerParticles, aod::JetParticles const& mcParticles)
  {
    // event selection
    if (!collision.selEv())
      return;

    // group collision
    auto const triggerParticlesThisEvent = triggerParticles.sliceBy(perColTriggerParticles, collision.mcCollisionId());
    auto const mcParticlesThisEvent = mcParticles.sliceBy(perColMcParticles, collision.mcCollisionId());

    // trigger pairing loop
    for (auto const& trigger : triggerParticlesThisEvent) {
      // trigger info
      histos.fill(HIST("mc/true_reco/corr/h3_ptPhiEta_trig"), trigger.pt(), trigger.phi(), trigger.eta());

      // hadrons (tracks) and pipm
      for (auto const& associated : mcParticlesThisEvent) {
        // exclude self correlation
        if (trigger.jetMcParticleId() == associated.globalIndex())
          continue;

        // standard particles (marked physical primary)
        if (checkPrimaryEtaMc(associated)) {
          // charged primary ('hadron') selection
          if (checkChargedMc(associated)) {
            histos.fill(HIST("mc/true_reco/corr/h3_ptPhiEta_assoc_hadron"), associated.pt(), associated.phi(), associated.eta());
            histos.fill(HIST("mc/true_reco/corr/h4_corr_hadron"),
                        getDeltaPhi(trigger.phi(), associated.phi()),
                        trigger.eta() - associated.eta(),
                        trigger.pt(), associated.pt());
          }

          // pipm selection
          if (std::abs(associated.pdgCode()) == PDG_t::kPiPlus) {
            histos.fill(HIST("mc/true_reco/corr/h3_ptPhiEta_assoc_pipm"), associated.pt(), associated.phi(), associated.eta());
            histos.fill(HIST("mc/true_reco/corr/h4_corr_pipm"),
                        getDeltaPhi(trigger.phi(), associated.phi()),
                        trigger.eta() - associated.eta(),
                        trigger.pt(), associated.pt());
          }

          // photon selection
          if (associated.pdgCode() == PDG_t::kGamma) {
            histos.fill(HIST("mc/true_reco/corr/h3_ptPhiEta_assoc_photon"), associated.pt(), associated.phi(), associated.eta());
            histos.fill(HIST("mc/true_reco/corr/h4_corr_photon"),
                        getDeltaPhi(trigger.phi(), associated.phi()),
                        trigger.eta() - associated.eta(),
                        trigger.pt(), associated.pt());
          }
        }

        // decaying particles (not marked physical primary)
        if ((std::abs(associated.eta()) < etaMax)) {
          // pi0 selection
          if (checkPi0ToGG(associated)) {
            histos.fill(HIST("mc/true_reco/corr/h3_ptPhiEta_assoc_pi0"), associated.pt(), associated.phi(), associated.eta());
            histos.fill(HIST("mc/true_reco/corr/h4_corr_pi0"),
                        getDeltaPhi(trigger.phi(), associated.phi()),
                        trigger.eta() - associated.eta(),
                        trigger.pt(), associated.pt());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processMcTrueRecoColCorr, "process mc-true (reco collisions) correlation for multiple associated particles", false);

  void processMcTrueEff(CorrMcCollision const& mcCollision, aod::JetParticles const& mcParticles)
  {
    // event selection
    if (doTrigEvEff && !mcCollision.trigEv())
      return;

    for (auto const& mcParticle : mcParticles) {
      // standard particles (marked physical primary)
      if (checkPrimaryEtaMc(mcParticle)) {
        // hadrons
        if (checkChargedMc(mcParticle)) {
          histos.fill(HIST("mc/eff/h3_ptPhiEta_mcTrue_hadron"), mcParticle.pt(), mcParticle.phi(), mcParticle.eta());
          histos.fill(HIST("mc/eff/h3_ptZPvMult_mcTrue_hadron"), mcParticle.pt(), mcCollision.posZ(), mcCollision.multMCNParticlesEta08());
        }
        // pipm
        if (std::abs(mcParticle.pdgCode()) == PDG_t::kPiPlus) {
          histos.fill(HIST("mc/eff/h3_ptPhiEta_mcTrue_pipm"), mcParticle.pt(), mcParticle.phi(), mcParticle.eta());
          histos.fill(HIST("mc/eff/h3_ptZPvMult_mcTrue_pipm"), mcParticle.pt(), mcCollision.posZ(), mcCollision.multMCNParticlesEta08());
        }
        // photons
        if (mcParticle.pdgCode() == PDG_t::kGamma) {
          histos.fill(HIST("mc/eff/h3_ptPhiEta_mcTrue_photon"), mcParticle.pt(), mcParticle.phi(), mcParticle.eta());
          histos.fill(HIST("mc/eff/h3_ptZPvMult_mcTrue_photon"), mcParticle.pt(), mcCollision.posZ(), mcCollision.multMCNParticlesEta08());
        }
      }

      // decaying particles (not marked physical primary)
      if ((std::abs(mcParticle.eta()) < etaMax)) {
        // pi0
        if (checkPi0ToGG(mcParticle)) {
          histos.fill(HIST("mc/eff/h3_ptPhiEta_mcTrue_pi0"), mcParticle.pt(), mcParticle.phi(), mcParticle.eta());
          histos.fill(HIST("mc/eff/h3_ptZPvMult_mcTrue_pi0"), mcParticle.pt(), mcCollision.posZ(), mcCollision.multMCNParticlesEta08());
        }
      }
    }
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processMcTrueEff, "process MC-true data (all collisions) to calculate efficiencies", false);

  void processMcRecoColEff(CorrMcDCollision const& collision, aod::JetMcCollisions const&, aod::JetTracksMCD const& tracks,
                           aod::Triggers const& triggers, aod::Hadrons const& hadrons, aod::Pipms const& pipms,
                           aod::PhotonPCMs const& photonPCMs, aod::PhotonPCMPairs const& photonPCMPairs,
                           aod::JetParticles const& mcParticles)
  {
    int excludeTriggerTrackId = -1;
    int excludeTriggerParticleId = -1;

    // event selection
    if (!collision.selEv())
      return;
    if (doTrigEvEff && !collision.trigEv())
      return;

    auto const mcParticlesThisEvent = mcParticles.sliceBy(perColMcParticles, collision.mcCollisionId());

    // random trigger
    if (doTrigEvEff) {
      std::uniform_int_distribution<int> intDistribution(0, static_cast<int>(triggers.size()) - 1);
      auto const& excludeTrigger = triggers.rawIteratorAt(intDistribution(randomEngine));
      if (excludeTrigger.jetTrack_as<aod::JetTracksMCD>().has_mcParticle()) {
        excludeTriggerParticleId = excludeTrigger.jetTrack_as<aod::JetTracksMCD>().mcParticleId();
        excludeTriggerTrackId = excludeTrigger.jetTrack_as<aod::JetTracksMCD>().globalIndex();
      }
    }

    // hadrons
    for (auto const& hadron : hadrons) {
      if (doTrigEvEff && hadron.jetTrackId() == excludeTriggerTrackId)
        continue;
      histos.fill(HIST("mc/eff/h3_ptPhiEta_mcReco_hadron"), hadron.pt(), hadron.phi(), hadron.eta());
      histos.fill(HIST("mc/eff/h3_ptZPvMult_mcReco_hadron"), hadron.pt(), collision.posZ(), collision.multNTracksGlobal());
      // purity
      if (!hadron.jetTrack_as<aod::JetTracksMCD>().has_mcParticle())
        continue;
      auto const hadronParticle = hadron.jetTrack_as<aod::JetTracksMCD>().mcParticle();
      if (!checkPrimaryTrackMc(hadronParticle))
        continue;
      if (requireSingleCollisionPurity && hadronParticle.mcCollisionId() != collision.mcCollisionId())
        continue;

      histos.fill(HIST("mc/eff/h3_ptPhiEta_mcReco_hasCorrectMc_hadron"), hadron.pt(), hadron.phi(), hadron.eta());
      histos.fill(HIST("mc/eff/h3_ptZPvMult_mcReco_hasCorrectMc_hadron"), hadron.pt(), collision.posZ(), collision.multNTracksGlobal());
    }

    // pipm
    for (auto const& pipm : pipms) {
      if (doTrigEvEff && pipm.jetTrackId() == excludeTriggerTrackId)
        continue;
      histos.fill(HIST("mc/eff/h3_ptPhiEta_mcReco_pipm"), pipm.pt(), pipm.phi(), pipm.eta());
      histos.fill(HIST("mc/eff/h3_ptZPvMult_mcReco_pipm"), pipm.pt(), collision.posZ(), collision.multNTracksGlobal());
      // purity
      if (!pipm.jetTrack_as<aod::JetTracksMCD>().has_mcParticle())
        continue;
      auto const pipmParticle = pipm.jetTrack_as<aod::JetTracksMCD>().mcParticle();
      if (std::abs(pipmParticle.pdgCode()) != PDG_t::kPiPlus || !checkPrimaryEtaMc(pipmParticle))
        continue;
      if (requireSingleCollisionPurity && pipmParticle.mcCollisionId() != collision.mcCollisionId())
        continue;

      histos.fill(HIST("mc/eff/h3_ptPhiEta_mcReco_hasCorrectMc_pipm"), pipm.pt(), pipm.phi(), pipm.eta());
      histos.fill(HIST("mc/eff/h3_ptZPvMult_mcReco_hasCorrectMc_pipm"), pipm.pt(), collision.posZ(), collision.multNTracksGlobal());
    }

    // photon mc checks

    auto const isConversionPhoton = [&](auto const& posTrack, auto const& negTrack) {
      // check same mother
      auto const& posMothers = posTrack.mcParticle().template mothers_as<aod::JetParticles>();
      auto const& negMothers = negTrack.mcParticle().template mothers_as<aod::JetParticles>();
      if (posMothers.size() != 1 || negMothers.size() != 1)
        return false;
      if (posMothers.begin()->globalIndex() != negMothers.begin()->globalIndex())
        return false;
      // check photon
      if (posMothers.begin()->pdgCode() != PDG_t::kGamma)
        return false;

      return true;
    };
    auto const isGGFromPi0 = [&](auto const& posTrack1, auto const& negTrack1, auto const& posTrack2, auto const& negTrack2) {
      if (!isConversionPhoton(posTrack1, negTrack1) || !isConversionPhoton(posTrack2, negTrack2))
        return false;
      // check same mother
      auto const& mothers1 = (*(posTrack1.mcParticle().template mothers_as<aod::JetParticles>().begin())).template mothers_as<aod::JetParticles>();
      auto const& mothers2 = (*(posTrack2.mcParticle().template mothers_as<aod::JetParticles>().begin())).template mothers_as<aod::JetParticles>();
      constexpr int NMothersPhotonFromPi0 = 2; // for some reason two mothers (same particle) for pi0 decays (contradicts PYTHIA documentation, but whatever)
      if (mothers1.size() != NMothersPhotonFromPi0 || mothers2.size() != NMothersPhotonFromPi0)
        return false;
      if (mothers1.begin()->globalIndex() != mothers2.begin()->globalIndex())
        return false;
      // check pi0
      if (mothers1.begin()->pdgCode() != PDG_t::kPi0)
        return false;

      return true;
    };

    // photonPCM
    for (auto const& photonPCM : photonPCMs) {
      histos.fill(HIST("mc/eff/h3_ptPhiEta_mcReco_photonPCM"), photonPCM.pt(), photonPCM.phi(), photonPCM.eta());
      histos.fill(HIST("mc/eff/h3_ptZPvMult_mcReco_photonPCM"), photonPCM.pt(), collision.posZ(), collision.multNTracksGlobal());

      // purity
      // (V0Legs does not have the tracks reference as index column (just int)??)
      auto const& posTrack = tracks.rawIteratorAt(photonPCM.posTrackId() - tracks.offset());
      auto const& negTrack = tracks.rawIteratorAt(photonPCM.negTrackId() - tracks.offset());
      if (!posTrack.has_mcParticle() || !negTrack.has_mcParticle())
        continue;
      if (!isConversionPhoton(posTrack, negTrack) || !checkPrimaryEtaMc(*(posTrack.mcParticle().mothers_as<aod::JetParticles>().begin())))
        continue;
      if (requireSingleCollisionPurity && posTrack.mcParticle().mcCollisionId() != collision.mcCollisionId())
        continue;

      histos.fill(HIST("mc/eff/h3_ptPhiEta_mcReco_hasCorrectMc_photonPCM"), photonPCM.pt(), photonPCM.phi(), photonPCM.eta());
      histos.fill(HIST("mc/eff/h3_ptZPvMult_mcReco_hasCorrectMc_photonPCM"), photonPCM.pt(), collision.posZ(), collision.multNTracksGlobal());
    }

    // pi0PCM
    for (auto const& photonPCMPair : photonPCMPairs) {
      if (photonPCMPair.mgg() < pi0PCMMassRange.value[0] || photonPCMPair.mgg() > pi0PCMMassRange.value[1])
        continue;

      histos.fill(HIST("mc/eff/h3_ptPhiEta_mcReco_pi0PCM"), photonPCMPair.pt(), photonPCMPair.phi(), photonPCMPair.eta());
      histos.fill(HIST("mc/eff/h3_ptZPvMult_mcReco_pi0PCM"), photonPCMPair.pt(), collision.posZ(), collision.multNTracksGlobal());

      // purity
      auto const& posTrack1 = tracks.rawIteratorAt(photonPCMPair.posTrack1Id() - tracks.offset());
      auto const& negTrack1 = tracks.rawIteratorAt(photonPCMPair.negTrack1Id() - tracks.offset());
      auto const& posTrack2 = tracks.rawIteratorAt(photonPCMPair.posTrack2Id() - tracks.offset());
      auto const& negTrack2 = tracks.rawIteratorAt(photonPCMPair.negTrack2Id() - tracks.offset());
      if (!posTrack1.has_mcParticle() || !negTrack1.has_mcParticle() || !posTrack2.has_mcParticle() || !negTrack2.has_mcParticle())
        continue;
      if (!isGGFromPi0(posTrack1, negTrack1, posTrack2, negTrack2) ||
          std::abs((*(posTrack1.mcParticle().mothers_as<aod::JetParticles>().begin())).mothers_as<aod::JetParticles>().begin()->eta()) > etaMax)
        continue;
      if (requireSingleCollisionPurity &&
          (posTrack1.mcParticle().mcCollisionId() != collision.mcCollisionId() || posTrack2.mcParticle().mcCollisionId() != collision.mcCollisionId()))
        continue;

      histos.fill(HIST("mc/eff/h3_ptPhiEta_mcReco_hasCorrectMc_pi0PCM"), photonPCMPair.pt(), photonPCMPair.phi(), photonPCMPair.eta());
      histos.fill(HIST("mc/eff/h3_ptZPvMult_mcReco_hasCorrectMc_pi0PCM"), photonPCMPair.pt(), collision.posZ(), collision.multNTracksGlobal());
    }

    // mcParticle loop
    for (auto const& mcParticle : mcParticlesThisEvent) {
      // standard particles (marked physical primary)
      if (checkPrimaryEtaMc(mcParticle)) {
        // hadrons
        if (checkChargedMc(mcParticle) && (!doTrigEvEff || mcParticle.globalIndex() != excludeTriggerParticleId)) {
          histos.fill(HIST("mc/eff/h3_ptPhiEta_mcTrue_recoCol_hadron"), mcParticle.pt(), mcParticle.phi(), mcParticle.eta());
          histos.fill(HIST("mc/eff/h3_ptZPvMult_mcTrue_recoCol_hadron"), mcParticle.pt(), collision.mcCollision().posZ(), collision.multNTracksGlobal());
        }
        // pipm
        if (std::abs(mcParticle.pdgCode()) == PDG_t::kPiPlus && (!doTrigEvEff || mcParticle.globalIndex() != excludeTriggerParticleId)) {
          histos.fill(HIST("mc/eff/h3_ptPhiEta_mcTrue_recoCol_pipm"), mcParticle.pt(), mcParticle.phi(), mcParticle.eta());
          histos.fill(HIST("mc/eff/h3_ptZPvMult_mcTrue_recoCol_pipm"), mcParticle.pt(), collision.mcCollision().posZ(), collision.multNTracksGlobal());
        }
        // photons
        if (mcParticle.pdgCode() == PDG_t::kGamma) {
          histos.fill(HIST("mc/eff/h3_ptPhiEta_mcTrue_recoCol_photon"), mcParticle.pt(), mcParticle.phi(), mcParticle.eta());
          histos.fill(HIST("mc/eff/h3_ptZPvMult_mcTrue_recoCol_photon"), mcParticle.pt(), collision.mcCollision().posZ(), collision.multNTracksGlobal());
        }
      }

      // decaying particles (not marked physical primary)
      if ((std::abs(mcParticle.eta()) < etaMax)) {
        // pi0
        if (checkPi0ToGG(mcParticle)) {
          histos.fill(HIST("mc/eff/h3_ptPhiEta_mcTrue_recoCol_pi0"), mcParticle.pt(), mcParticle.phi(), mcParticle.eta());
          histos.fill(HIST("mc/eff/h3_ptZPvMult_mcTrue_recoCol_pi0"), mcParticle.pt(), collision.mcCollision().posZ(), collision.multNTracksGlobal());
        }
      }
    }
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processMcRecoColEff, "process MC data to calculate efficiencies and purities", false);

  // test /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void processTest(CorrCollision const& collision,
                   soa::Join<aod::JetTracks, aod::JTrackPIs> const& tracks, soa::Join<aod::Tracks, aod::TracksExtra> const&,
                   aod::Hadrons const& hadrons)
  {
    // event selection
    if (!collision.selEv())
      return;

    histos.fill(HIST("test/h2_mult_comp"), collision.multNTracksGlobal(), hadrons.size());

    for (auto const& track : tracks) {
      auto const fullTrack = track.track_as<soa::Join<aod::Tracks, aod::TracksExtra>>();

      constexpr float Mincrossedrows = 40;
      constexpr float Maxchi2tpc = 5.0;
      constexpr float Maxchi2its = 6.0;
      constexpr float MaxR = 83.1;
      constexpr float MinPtTrackiu = 0.1;

      if (!fullTrack.hasITS() && !fullTrack.hasTPC())
        continue;
      if (fullTrack.x() * fullTrack.x() + fullTrack.y() * fullTrack.y() > MaxR * MaxR || fullTrack.pt() < MinPtTrackiu)
        continue;
      if (fullTrack.hasTPC()) {
        if (fullTrack.tpcNClsCrossedRows() < Mincrossedrows || fullTrack.tpcChi2NCl() > Maxchi2tpc)
          continue;
      }
      if (fullTrack.hasITS()) {
        if (fullTrack.itsChi2NCl() > Maxchi2its)
          continue;
      }

      histos.fill(HIST("test/h2_tracks_zPvMultDep"), collision.posZ(), collision.multNTracksGlobal());
    }

    histos.fill(HIST("test/h2_globalTracks_zPvMultDep"), collision.posZ(), collision.multNTracksGlobal(), hadrons.size());
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processTest, "process just to test things", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& configContext)
{
  return WorkflowSpec{
    adaptAnalysisTask<CorrelationTableProducer>(configContext),
    adaptAnalysisTask<PhotonChargedTriggerCorrelation>(configContext)};
}
