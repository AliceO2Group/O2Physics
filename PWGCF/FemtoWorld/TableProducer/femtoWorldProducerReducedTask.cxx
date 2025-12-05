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

/// \file femtoWorldProducerReducedTask.cxx
/// \brief Tasks that produces the track tables used for the pairing (tracks only)
/// \author Luca Barioglio, TU München, andreas.mathis@ph.tum.de

#include "PWGCF/FemtoWorld/Core/FemtoWorldCollisionSelection.h"
#include "PWGCF/FemtoWorld/Core/FemtoWorldTrackSelection.h"
#include "PWGCF/FemtoWorld/DataModel/FemtoWorldDerived.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"
#include <CCDB/BasicCCDBManager.h>

#include "Math/Vector4D.h"
#include "TMath.h"

using namespace o2;
using namespace o2::analysis::femtoWorld;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{

using FemtoFullCollision = soa::Join<aod::Collisions,
                                     aod::EvSels,
                                     aod::Mults>::iterator;
using FemtoFullTracks = soa::Join<aod::FullTracks,
                                  aod::TracksDCA, aod::TOFSignal,
                                  aod::pidTPCEl, aod::pidTPCMu, aod::pidTPCPi,
                                  aod::pidTPCKa, aod::pidTPCPr, aod::pidTPCDe,
                                  aod::pidTOFEl, aod::pidTOFMu, aod::pidTOFPi,
                                  aod::pidTOFKa, aod::pidTOFPr, aod::pidTOFDe>;
} // namespace o2::aod

struct femtoWorldProducerReducedTask {

  Produces<aod::FemtoWorldCollisions> outputCollision;
  Produces<aod::FemtoWorldParticles> outputTracks;
  Produces<aod::FemtoWorldDebugParticles> outputDebugTracks;

  Configurable<bool> ConfDebugOutput{"ConfDebugOutput", true, "Debug output"};

  // Choose if filtering or skimming version is run
  Configurable<bool> ConfIsTrigger{"ConfIsTrigger", false, "Store all collisions"};

  // Choose if running on converted data or pilot beam
  Configurable<bool> ConfIsRun3{"ConfIsRun3", false, "Running on Pilot beam"};

  /// Event cuts
  FemtoWorldCollisionSelection colCuts;
  Configurable<float> ConfEvtZvtx{"ConfEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};
  Configurable<bool> ConfEvtTriggerCheck{"ConfEvtTriggerCheck", true, "Evt sel: check for trigger"};
  Configurable<int> ConfEvtTriggerSel{"ConfEvtTriggerSel", kINT7, "Evt sel: trigger"};
  Configurable<bool> ConfEvtOfflineCheck{"ConfEvtOfflineCheck", false, "Evt sel: check for offline selection"};

  FemtoWorldTrackSelection trackCuts;
  Configurable<std::vector<float>> ConfTrkCharge{FemtoWorldTrackSelection::getSelectionName(femtoWorldTrackSelection::kSign, "ConfTrk"), std::vector<float>{-1, 1}, FemtoWorldTrackSelection::getSelectionHelper(femtoWorldTrackSelection::kSign, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkPtmin{FemtoWorldTrackSelection::getSelectionName(femtoWorldTrackSelection::kpTMin, "ConfTrk"), std::vector<float>{0.4f, 0.6f, 0.5f}, FemtoWorldTrackSelection::getSelectionHelper(femtoWorldTrackSelection::kpTMin, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkEta{FemtoWorldTrackSelection::getSelectionName(femtoWorldTrackSelection::kEtaMax, "ConfTrk"), std::vector<float>{0.8f, 0.7f, 0.9f}, FemtoWorldTrackSelection::getSelectionHelper(femtoWorldTrackSelection::kEtaMax, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkTPCnclsMin{FemtoWorldTrackSelection::getSelectionName(femtoWorldTrackSelection::kTPCnClsMin, "ConfTrk"), std::vector<float>{80.f, 70.f, 60.f}, FemtoWorldTrackSelection::getSelectionHelper(femtoWorldTrackSelection::kTPCnClsMin, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkTPCfCls{FemtoWorldTrackSelection::getSelectionName(femtoWorldTrackSelection::kTPCfClsMin, "ConfTrk"), std::vector<float>{0.7f, 0.83f, 0.9f}, FemtoWorldTrackSelection::getSelectionHelper(femtoWorldTrackSelection::kTPCfClsMin, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkTPCcRowsMin{FemtoWorldTrackSelection::getSelectionName(femtoWorldTrackSelection::kTPCcRowsMin, "ConfTrk"), std::vector<float>{70.f, 60.f, 80.f}, FemtoWorldTrackSelection::getSelectionHelper(femtoWorldTrackSelection::kTPCcRowsMin, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkTPCsCls{FemtoWorldTrackSelection::getSelectionName(femtoWorldTrackSelection::kTPCsClsMax, "ConfTrk"), std::vector<float>{0.1f, 160.f}, FemtoWorldTrackSelection::getSelectionHelper(femtoWorldTrackSelection::kTPCsClsMax, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkITSnclsMin{FemtoWorldTrackSelection::getSelectionName(femtoWorldTrackSelection::kITSnClsMin, "ConfTrk"), std::vector<float>{-1.f, 2.f, 4.f}, FemtoWorldTrackSelection::getSelectionHelper(femtoWorldTrackSelection::kITSnClsMin, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkITSnclsIbMin{FemtoWorldTrackSelection::getSelectionName(femtoWorldTrackSelection::kITSnClsIbMin, "ConfTrk"), std::vector<float>{-1.f, 1.f}, FemtoWorldTrackSelection::getSelectionHelper(femtoWorldTrackSelection::kITSnClsIbMin, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkDCAxyMax{FemtoWorldTrackSelection::getSelectionName(femtoWorldTrackSelection::kDCAxyMax, "ConfTrk"), std::vector<float>{0.1f, 3.5f}, FemtoWorldTrackSelection::getSelectionHelper(femtoWorldTrackSelection::kDCAxyMax, "Track selection: ")}; /// here we need an open cut to do the DCA fits later on!
  Configurable<std::vector<float>> ConfTrkDCAzMax{FemtoWorldTrackSelection::getSelectionName(femtoWorldTrackSelection::kDCAzMax, "ConfTrk"), std::vector<float>{0.2f, 3.5f}, FemtoWorldTrackSelection::getSelectionHelper(femtoWorldTrackSelection::kDCAzMax, "Track selection: ")};
  /// \todo Reintegrate PID to the general selection container
  Configurable<std::vector<float>> ConfTrkPIDnSigmaMax{FemtoWorldTrackSelection::getSelectionName(femtoWorldTrackSelection::kPIDnSigmaMax, "ConfTrk"), std::vector<float>{3.5f, 3.f, 2.5f}, FemtoWorldTrackSelection::getSelectionHelper(femtoWorldTrackSelection::kPIDnSigmaMax, "Track selection: ")};
  Configurable<std::vector<int>> ConfTrkTPIDspecies{"ConfTrkTPIDspecies", std::vector<int>{o2::track::PID::Pion, o2::track::PID::Kaon, o2::track::PID::Proton, o2::track::PID::Deuteron}, "Trk sel: Particles species for PID"};

  HistogramRegistry qaRegistry{"QAHistos", {}, OutputObjHandlingPolicy::QAObject};

  int mRunNumber;
  float mMagField;
  Service<o2::ccdb::BasicCCDBManager> ccdb; /// Accessing the CCDB

  void init(InitContext&)
  {
    colCuts.setCuts(ConfEvtZvtx, ConfEvtTriggerCheck, ConfEvtTriggerSel, ConfEvtOfflineCheck, ConfIsRun3);
    colCuts.init(&qaRegistry);

    trackCuts.setSelection(ConfTrkCharge, femtoWorldTrackSelection::kSign, femtoWorldSelection::kEqual);
    trackCuts.setSelection(ConfTrkPtmin, femtoWorldTrackSelection::kpTMin, femtoWorldSelection::kLowerLimit);
    trackCuts.setSelection(ConfTrkEta, femtoWorldTrackSelection::kEtaMax, femtoWorldSelection::kAbsUpperLimit);
    trackCuts.setSelection(ConfTrkTPCnclsMin, femtoWorldTrackSelection::kTPCnClsMin, femtoWorldSelection::kLowerLimit);
    if (!ConfIsRun3) {
      trackCuts.setSelection(ConfTrkTPCfCls, femtoWorldTrackSelection::kTPCfClsMin, femtoWorldSelection::kLowerLimit);
    }
    trackCuts.setSelection(ConfTrkTPCcRowsMin, femtoWorldTrackSelection::kTPCcRowsMin, femtoWorldSelection::kLowerLimit);
    trackCuts.setSelection(ConfTrkTPCsCls, femtoWorldTrackSelection::kTPCsClsMax, femtoWorldSelection::kUpperLimit);
    if (ConfIsRun3) {
      trackCuts.setSelection(ConfTrkITSnclsMin, femtoWorldTrackSelection::kITSnClsMin, femtoWorldSelection::kLowerLimit);
      trackCuts.setSelection(ConfTrkITSnclsIbMin, femtoWorldTrackSelection::kITSnClsIbMin, femtoWorldSelection::kLowerLimit);
    }
    trackCuts.setSelection(ConfTrkDCAxyMax, femtoWorldTrackSelection::kDCAxyMax, femtoWorldSelection::kAbsUpperLimit);
    trackCuts.setSelection(ConfTrkDCAzMax, femtoWorldTrackSelection::kDCAzMax, femtoWorldSelection::kAbsUpperLimit);
    trackCuts.setSelection(ConfTrkPIDnSigmaMax, femtoWorldTrackSelection::kPIDnSigmaMax, femtoWorldSelection::kAbsUpperLimit);
    trackCuts.setPIDSpecies(ConfTrkTPIDspecies);
    trackCuts.init<aod::femtoworldparticle::ParticleType::kTrack, aod::femtoworldparticle::TrackType::kNoChild, aod::femtoworldparticle::cutContainerType>(&qaRegistry);
    mRunNumber = 0;
    mMagField = 0.0;
    /// Initializing CCDB
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    long now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);
  }

  /// Function to retrieve the nominal mgnetic field in kG (0.1T) and convert it directly to T
  float getMagneticFieldTesla(uint64_t timestamp)
  {
    // TODO done only once (and not per run). Will be replaced by CCDBConfigurable
    static o2::parameters::GRPObject* grpo = nullptr;
    if (grpo == nullptr) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>("GLO/GRP/GRP", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 0;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
    }
    float output = 0.1 * (grpo->getNominalL3Field());
    return output;
  }

  void process(aod::FemtoFullCollision const& col, aod::BCsWithTimestamps const&, aod::FemtoFullTracks const& tracks) /// \todo with FilteredFullV0s
  {
    // get magnetic field for run
    auto bc = col.bc_as<aod::BCsWithTimestamps>();
    if (mRunNumber != bc.runNumber()) {
      mMagField = getMagneticFieldTesla(bc.timestamp());
      mRunNumber = bc.runNumber();
    }
    const auto vtxZ = col.posZ();
    const auto spher = colCuts.computeSphericity(col, tracks);
    /// For benchmarking on Run 2, V0M in FemtoWorldRun2 is defined V0M/2
    int mult = 0;
    if (ConfIsRun3) {
      mult = col.multFT0M(); /// Mult based on T0, temporary storing to be fixed and checked
    } else {
      mult = 0.5 * (col.multFV0M());
    }
    /// First thing to do is to check whether the basic event selection criteria are fulfilled
    // If the basic selection is NOT fulfilled:
    // in case of skimming run - don't store such collisions
    // in case of trigger run - store such collisions but don't store any particle candidates for such collisions
    if (!colCuts.isSelected(col)) {
      if (ConfIsTrigger) {
        outputCollision(col.posZ(), mult, colCuts.computeSphericity(col, tracks), mMagField);
      }
      return;
    }

    colCuts.fillQA(col);
    // now the table is filled
    outputCollision(vtxZ, mult, spher, mMagField);

    int childIDs[2] = {0, 0}; // these IDs are necessary to keep track of the children
    for (auto& track : tracks) {
      /// if the most open selection criteria are not fulfilled there is no point looking further at the track
      if (!trackCuts.isSelectedMinimal(track)) {
        continue;
      }
      trackCuts.fillQA<aod::femtoworldparticle::ParticleType::kTrack, aod::femtoworldparticle::TrackType::kNoChild>(track);
      // the bit-wise container of the systematic variations is obtained
      auto cutContainer = trackCuts.getCutContainer<aod::femtoworldparticle::cutContainerType>(track);

      // now the table is filled
      outputTracks(outputCollision.lastIndex(),
                   track.pt(),
                   track.eta(),
                   track.phi(),
                   aod::femtoworldparticle::ParticleType::kTrack,
                   cutContainer.at(femtoWorldTrackSelection::TrackContainerPosition::kCuts),
                   cutContainer.at(femtoWorldTrackSelection::TrackContainerPosition::kPID),
                   track.dcaXY(), childIDs, 0, 0);
      if (ConfDebugOutput) {
        outputDebugTracks(track.sign(),
                          (uint8_t)track.tpcNClsFound(),
                          track.tpcNClsFindable(),
                          (uint8_t)track.tpcNClsCrossedRows(),
                          track.tpcNClsShared(),
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
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<femtoWorldProducerReducedTask>(cfgc)};
  return workflow;
}
