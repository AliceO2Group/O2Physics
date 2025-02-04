// Copyright 2019-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file femtoUniverseProducerReducedTask.cxx
/// \brief Tasks that produces the track tables used for the pairing (tracks only)
/// \author Luca Barioglio, TU München, andreas.mathis@ph.tum.de
/// \author Georgios Mantzaridis, TU München, georgios.mantzaridis@tum.de
/// \author Anton Riedel, TU München, anton.riedel@tum.de
/// \author Zuzanna Chochulska, WUT Warsaw & CTU Prague, zchochul@cern.ch

#include <CCDB/BasicCCDBManager.h>
#include <vector>

#include "TMath.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseCollisionSelection.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseTrackSelection.h"
#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "Math/Vector4D.h"
#include "PWGCF/FemtoUniverse/Core/femtoUtils.h"

using namespace o2;
using namespace o2::analysis::femto_universe;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{

using FemtoFullCollision = soa::Join<aod::Collisions,
                                     aod::EvSels,
                                     aod::Mults>::iterator;
using FemtoFullCollisionMC = soa::Join<aod::Collisions,
                                       aod::EvSels,
                                       aod::Mults,
                                       aod::McCollisionLabels>::iterator;

using FemtoFullTracks = soa::Join<aod::FullTracks, aod::TracksDCA,
                                  aod::pidTPCEl, aod::pidTPCMu, aod::pidTPCPi,
                                  aod::pidTPCKa, aod::pidTPCPr, aod::pidTPCDe,
                                  aod::pidTOFEl, aod::pidTOFMu, aod::pidTOFPi,
                                  aod::pidTOFKa, aod::pidTOFPr, aod::pidTOFDe>;
} // namespace o2::aod

struct femtoUniverseProducerReducedTask {

  Produces<aod::FdCollisions> outputCollision;
  Produces<aod::FDParticles> outputParts;
  Produces<aod::FdMCParticles> outputPartsMC;
  Produces<aod::FDExtParticles> outputDebugParts;
  Produces<aod::FDMCLabels> outputPartsMCLabels;
  Produces<aod::FDExtMCParticles> outputDebugPartsMC;

  Configurable<bool> ConfDebugOutput{"ConfDebugOutput", true, "Debug output"};
  Configurable<bool> ConfIsTrigger{"ConfIsTrigger", false, "Store all collisions"}; // Choose if filtering or skimming version is run
  Configurable<bool> ConfIsRun3{"ConfIsRun3", false, "Running on Run3 or Run2"};    // Choose if running on converted data or Run3
  Configurable<bool> ConfIsMC{"ConfIsMC", false, "Running on MC; implemented only for Run3"};

  // Event cuts
  FemtoUniverseCollisionSelection colCuts;
  Configurable<bool> ConfUseTPCmult{"ConfUseTPCmult", false, "Use multiplicity based on the number of tracks with TPC information"};
  Configurable<float> ConfEvtZvtx{"ConfEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};
  Configurable<bool> ConfEvtTriggerCheck{"ConfEvtTriggerCheck", true, "Evt sel: check for trigger"};
  Configurable<int> ConfEvtTriggerSel{"ConfEvtTriggerSel", kINT7, "Evt sel: trigger"};
  Configurable<bool> ConfEvtOfflineCheck{"ConfEvtOfflineCheck", false, "Evt sel: check for offline selection"};
  Configurable<float> ConfCentFT0Min{"ConfCentFT0Min", 0.f, "Min CentFT0 value for centrality selection"};
  Configurable<float> ConfCentFT0Max{"ConfCentFT0Max", 200.f, "Max CentFT0 value for centrality selection"};

  Configurable<bool> ConfRejectNotPropagatedTracks{"ConfRejectNotPropagatedTracks", false, "True: reject not propagated tracks"};
  Configurable<bool> ConfRejectITSHitandTOFMissing{"ConfRejectITSHitandTOFMissing", false, "True: reject if neither ITS hit nor TOF timing satisfied"};

  Configurable<bool> ConfForceGRP{"ConfForceGRP", false, "Set true if the magnetic field configuration is not available in the usual CCDB directory (e.g. for Run 2 converted data or unanchorad Monte Carlo)"};
  Configurable<int> ConfPDGCodeTrack{"ConfPDGCodeTrack", 2212, "PDG code of the selected track for Monte Carlo truth"};
  // Track cuts
  FemtoUniverseTrackSelection trackCuts;
  Configurable<std::vector<float>> ConfTrkCharge{FemtoUniverseTrackSelection::getSelectionName(femto_universe_track_selection::kSign, "ConfTrk"), std::vector<float>{-1, 1}, FemtoUniverseTrackSelection::getSelectionHelper(femto_universe_track_selection::kSign, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkPtmin{FemtoUniverseTrackSelection::getSelectionName(femto_universe_track_selection::kpTMin, "ConfTrk"), std::vector<float>{0.4f, 0.6f, 0.5f}, FemtoUniverseTrackSelection::getSelectionHelper(femto_universe_track_selection::kpTMin, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkPtmax{FemtoUniverseTrackSelection::getSelectionName(femto_universe_track_selection::kpTMax, "ConfTrk"), std::vector<float>{5.4f, 5.6f, 5.5f}, FemtoUniverseTrackSelection::getSelectionHelper(femto_universe_track_selection::kpTMax, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkEta{FemtoUniverseTrackSelection::getSelectionName(femto_universe_track_selection::kEtaMax, "ConfTrk"), std::vector<float>{0.8f, 0.7f, 0.9f}, FemtoUniverseTrackSelection::getSelectionHelper(femto_universe_track_selection::kEtaMax, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkTPCnclsMin{FemtoUniverseTrackSelection::getSelectionName(femto_universe_track_selection::kTPCnClsMin, "ConfTrk"), std::vector<float>{80.f, 70.f, 60.f}, FemtoUniverseTrackSelection::getSelectionHelper(femto_universe_track_selection::kTPCnClsMin, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkTPCfCls{FemtoUniverseTrackSelection::getSelectionName(femto_universe_track_selection::kTPCfClsMin, "ConfTrk"), std::vector<float>{0.7f, 0.83f, 0.9f}, FemtoUniverseTrackSelection::getSelectionHelper(femto_universe_track_selection::kTPCfClsMin, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkTPCcRowsMin{FemtoUniverseTrackSelection::getSelectionName(femto_universe_track_selection::kTPCcRowsMin, "ConfTrk"), std::vector<float>{70.f, 60.f, 80.f}, FemtoUniverseTrackSelection::getSelectionHelper(femto_universe_track_selection::kTPCcRowsMin, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkTPCsCls{FemtoUniverseTrackSelection::getSelectionName(femto_universe_track_selection::kTPCsClsMax, "ConfTrk"), std::vector<float>{0.1f, 160.f}, FemtoUniverseTrackSelection::getSelectionHelper(femto_universe_track_selection::kTPCsClsMax, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkTPCfracsCls{FemtoUniverseTrackSelection::getSelectionName(femto_universe_track_selection::kTPCfracsClsMax, "ConfTrk"), std::vector<float>{0.1f, 160.f}, FemtoUniverseTrackSelection::getSelectionHelper(femto_universe_track_selection::kTPCfracsClsMax, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkITSnclsMin{FemtoUniverseTrackSelection::getSelectionName(femto_universe_track_selection::kITSnClsMin, "ConfTrk"), std::vector<float>{-1.f, 2.f, 4.f}, FemtoUniverseTrackSelection::getSelectionHelper(femto_universe_track_selection::kITSnClsMin, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkITSnclsIbMin{FemtoUniverseTrackSelection::getSelectionName(femto_universe_track_selection::kITSnClsIbMin, "ConfTrk"), std::vector<float>{-1.f, 1.f}, FemtoUniverseTrackSelection::getSelectionHelper(femto_universe_track_selection::kITSnClsIbMin, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkDCAxyMax{FemtoUniverseTrackSelection::getSelectionName(femto_universe_track_selection::kDCAxyMax, "ConfTrk"), std::vector<float>{0.1f, 0.5f}, FemtoUniverseTrackSelection::getSelectionHelper(femto_universe_track_selection::kDCAxyMax, "Track selection: ")}; /// here we need an open cut to do the DCA fits later on!
  Configurable<std::vector<float>> ConfTrkDCAzMax{FemtoUniverseTrackSelection::getSelectionName(femto_universe_track_selection::kDCAzMax, "ConfTrk"), std::vector<float>{0.2f, 0.5f}, FemtoUniverseTrackSelection::getSelectionHelper(femto_universe_track_selection::kDCAzMax, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkPIDnSigmaMax{FemtoUniverseTrackSelection::getSelectionName(femto_universe_track_selection::kPIDnSigmaMax, "Conf"), std::vector<float>{3.5f, 3.f, 2.5f}, FemtoUniverseTrackSelection::getSelectionHelper(femto_universe_track_selection::kPIDnSigmaMax, "Track selection: ")};
  // off set the center of the nsigma distribution to deal with bad TPC/TOF calibration
  Configurable<float> ConfPIDnSigmaOffsetTPC{"ConfPIDnSigmaOffsetTPC", 0., "Offset for TPC nSigma because of bad calibration"};
  Configurable<float> ConfPIDnSigmaOffsetTOF{"ConfPIDnSigmaOffsetTOF", 0., "Offset for TOF nSigma because of bad calibration"};
  Configurable<std::vector<int>> ConfPIDspecies{"ConfPIDspecies", std::vector<int>{o2::track::PID::Pion, o2::track::PID::Kaon, o2::track::PID::Proton, o2::track::PID::Deuteron}, "Trk sel: Particles species for PID"};

  HistogramRegistry qaRegistry{"QAHistos", {}, OutputObjHandlingPolicy::QAObject};

  int mRunNumber;
  float mMagField;
  Service<o2::ccdb::BasicCCDBManager> ccdb; /// Accessing the CCDB

  void init(InitContext&)
  {
    colCuts.setCuts(ConfEvtZvtx, ConfEvtTriggerCheck, ConfEvtTriggerSel, ConfEvtOfflineCheck, ConfIsRun3, ConfCentFT0Min, ConfCentFT0Max);
    colCuts.init(&qaRegistry);

    trackCuts.setSelection(ConfTrkCharge, femto_universe_track_selection::kSign, femto_universe_selection::kEqual);
    trackCuts.setSelection(ConfTrkPtmin, femto_universe_track_selection::kpTMin, femto_universe_selection::kLowerLimit);
    trackCuts.setSelection(ConfTrkPtmax, femto_universe_track_selection::kpTMax, femto_universe_selection::kUpperLimit);
    trackCuts.setSelection(ConfTrkEta, femto_universe_track_selection::kEtaMax, femto_universe_selection::kAbsUpperLimit);
    trackCuts.setSelection(ConfTrkTPCnclsMin, femto_universe_track_selection::kTPCnClsMin, femto_universe_selection::kLowerLimit);
    trackCuts.setSelection(ConfTrkTPCfCls, femto_universe_track_selection::kTPCfClsMin, femto_universe_selection::kLowerLimit);
    trackCuts.setSelection(ConfTrkTPCcRowsMin, femto_universe_track_selection::kTPCcRowsMin, femto_universe_selection::kLowerLimit);
    trackCuts.setSelection(ConfTrkTPCsCls, femto_universe_track_selection::kTPCsClsMax, femto_universe_selection::kUpperLimit);
    trackCuts.setSelection(ConfTrkTPCfracsCls, femto_universe_track_selection::kTPCfracsClsMax, femto_universe_selection::kUpperLimit);
    trackCuts.setSelection(ConfTrkITSnclsMin, femto_universe_track_selection::kITSnClsMin, femto_universe_selection::kLowerLimit);
    trackCuts.setSelection(ConfTrkITSnclsIbMin, femto_universe_track_selection::kITSnClsIbMin, femto_universe_selection::kLowerLimit);
    trackCuts.setSelection(ConfTrkDCAxyMax, femto_universe_track_selection::kDCAxyMax, femto_universe_selection::kAbsUpperLimit);
    trackCuts.setSelection(ConfTrkDCAzMax, femto_universe_track_selection::kDCAzMax, femto_universe_selection::kAbsUpperLimit);
    trackCuts.setSelection(ConfTrkPIDnSigmaMax, femto_universe_track_selection::kPIDnSigmaMax, femto_universe_selection::kAbsUpperLimit);
    trackCuts.setPIDSpecies(ConfPIDspecies);
    trackCuts.setnSigmaPIDOffset(ConfPIDnSigmaOffsetTPC, ConfPIDnSigmaOffsetTOF);
    trackCuts.init<aod::femtouniverseparticle::ParticleType::kTrack,
                   aod::femtouniverseparticle::TrackType::kNoChild,
                   aod::femtouniverseparticle::CutContainerType>(&qaRegistry);
    mRunNumber = 0;
    mMagField = 0.0;
    /// Initializing CCDB
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);
  }

  /// Function to retrieve the nominal magnetic field in kG (0.1T) and convert it directly to T
  void getMagneticFieldTesla(aod::BCsWithTimestamps::iterator bc)
  {
    // TODO done only once (and not per run). Will be replaced by CCDBConfigurable
    // get magnetic field for run
    if (mRunNumber == bc.runNumber())
      return;
    auto timestamp = bc.timestamp();
    float output = -999;

    if (ConfIsRun3 && !ConfForceGRP) {
      static o2::parameters::GRPMagField* grpo = nullptr;
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with L3 ", timestamp, grpo->getL3Current());
      // taken from GRP onject definition of getNominalL3Field; update later to something smarter (mNominalL3Field = std::lround(5.f * mL3Current / 30000.f);)
      auto NominalL3Field = std::lround(5.f * grpo->getL3Current() / 30000.f);
      output = 0.1 * (NominalL3Field);

    } else {

      static o2::parameters::GRPObject* grpo = nullptr;
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>("GLO/GRP/GRP", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
      output = 0.1 * (grpo->getNominalL3Field());
    }
    mMagField = output;
    mRunNumber = bc.runNumber();
  }

  template <typename CandidateType>
  void fillMCParticle(CandidateType const& candidate, o2::aod::femtouniverseparticle::ParticleType fdparttype)
  {
    if (candidate.has_mcParticle()) {
      // get corresponding MC particle and its info
      auto particleMC = candidate.mcParticle();
      auto pdgCode = particleMC.pdgCode();

      int particleOrigin = 99;
      auto motherparticleMC = particleMC.template mothers_as<aod::McParticles>().front();

      if (abs(pdgCode) == abs(ConfPDGCodeTrack.value)) {

        if (particleMC.isPhysicalPrimary()) {
          particleOrigin = aod::femtouniverse_mc_particle::ParticleOriginMCTruth::kPrimary;
        } else if (motherparticleMC.producedByGenerator()) {
          particleOrigin = checkDaughterType(fdparttype, motherparticleMC.pdgCode());
        } else {
          particleOrigin = aod::femtouniverse_mc_particle::ParticleOriginMCTruth::kMaterial;
        }

      } else {

        particleOrigin = aod::femtouniverse_mc_particle::ParticleOriginMCTruth::kFake;
      }

      outputPartsMC(particleOrigin, pdgCode, particleMC.pt(), particleMC.eta(), particleMC.phi());
      outputPartsMCLabels(outputPartsMC.lastIndex());
    } else {
      outputPartsMCLabels(-1);
    }
  }

  template <bool isMC, typename CollisionType, typename TrackType>
  void fillCollisionsAndTracks(CollisionType const& col, TrackType const& tracks) /// \todo with FilteredFullV0s
  {
    // get magnetic field for run

    const auto vtxZ = col.posZ();
    const auto spher = colCuts.computeSphericity(col, tracks);

    int mult = 0;
    int multNtr = 0;
    if (ConfIsRun3) {
      mult = col.multFV0M();
      multNtr = col.multNTracksPV();
    } else {
      mult = 0.5 * (col.multFV0M()); /// For benchmarking on Run 2, V0M in FemtoUniverseRun2 is defined V0M/2
      multNtr = col.multTracklets();
    }
    if (ConfUseTPCmult) {
      multNtr = col.multTPC();
    }

    /// First thing to do is to check whether the basic event selection criteria are fulfilled
    // If the basic selection is NOT fulfilled:
    // in case of skimming run - don't store such collisions
    // in case of trigger run - store such collisions but don't store any particle candidates for such collisions
    if (!colCuts.isSelected(col)) {
      if (ConfIsTrigger) {
        outputCollision(vtxZ, mult, multNtr, colCuts.computeSphericity(col, tracks), mMagField);
      }
      return;
    }

    colCuts.fillQA(col);
    // now the table is filled
    outputCollision(vtxZ, mult, multNtr, spher, mMagField);

    // these IDs are necessary to keep track of the children
    // since this producer only produces the tables for tracks, there are no children
    std::vector<int> childIDs = {0, 0};
    for (auto& track : tracks) {
      /// if the most open selection criteria are not fulfilled there is no point looking further at the track
      if (!trackCuts.isSelectedMinimal(track)) {
        continue;
      }
      trackCuts.fillQA<aod::femtouniverseparticle::ParticleType::kTrack, aod::femtouniverseparticle::TrackType::kNoChild>(track);
      // an array of two bit-wise containers of the systematic variations is obtained
      // one container for the track quality cuts and one for the PID cuts
      auto cutContainer = trackCuts.getCutContainer<aod::femtouniverseparticle::CutContainerType>(track);

      // now the table is filled
      outputParts(outputCollision.lastIndex(),
                  track.pt(),
                  track.eta(),
                  track.phi(),
                  aod::femtouniverseparticle::ParticleType::kTrack,
                  cutContainer.at(femto_universe_track_selection::TrackContainerPosition::kCuts),
                  cutContainer.at(femto_universe_track_selection::TrackContainerPosition::kPID),
                  track.dcaXY(), childIDs, 0, 0);

      if constexpr (isMC) {
        fillMCParticle(track, o2::aod::femtouniverseparticle::ParticleType::kTrack);
      }

      if (ConfDebugOutput) {
        outputDebugParts(track.sign(),
                         (uint8_t)track.tpcNClsFound(),
                         track.tpcNClsFindable(),
                         (uint8_t)track.tpcNClsCrossedRows(),
                         track.tpcNClsShared(),
                         track.tpcFractionSharedCls(),
                         track.tpcInnerParam(),
                         track.itsNCls(),
                         track.itsNClsInnerBarrel(),
                         track.dcaXY(),
                         track.dcaZ(),
                         track.tpcSignal(),
                         track.tpcNSigmaStoreEl(),
                         track.tpcNSigmaStorePi(),
                         track.tpcNSigmaStoreKa(),
                         track.tpcNSigmaStorePr(),
                         track.tpcNSigmaStoreDe(),
                         track.tofNSigmaStoreEl(),
                         track.tofNSigmaStorePi(),
                         track.tofNSigmaStoreKa(),
                         track.tofNSigmaStorePr(),
                         track.tofNSigmaStoreDe(),
                         -999.,
                         -999.,
                         -999.,
                         -999.,
                         -999.,
                         -999.);
      }
    }
  }

  void processData(aod::FemtoFullCollision const& col, aod::BCsWithTimestamps const&, aod::FemtoFullTracks const& tracks)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    fillCollisionsAndTracks<false>(col, tracks);
  }
  PROCESS_SWITCH(femtoUniverseProducerReducedTask, processData, "Provide experimental data", true);

  void processMC(aod::FemtoFullCollisionMC const& col,
                 aod::BCsWithTimestamps const&,
                 soa::Join<aod::FemtoFullTracks, aod::McTrackLabels> const& tracks,
                 aod::McCollisions const&, aod::McParticles const&)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    fillCollisionsAndTracks<true>(col, tracks);
  }
  PROCESS_SWITCH(femtoUniverseProducerReducedTask, processMC, "Provide MC data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<femtoUniverseProducerReducedTask>(cfgc)};
  return workflow;
}
