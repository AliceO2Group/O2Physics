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

/// \file femtoDreamProducerTask.cxx
/// \brief Tasks that produces the track tables used for the pairing
/// \author Laura Serksnyte, TU München, laura.serksnyte@tum.de

#include <CCDB/BasicCCDBManager.h>
#include <cstdint>
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/Core/trackUtilities.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGCF/DataModel/FemtoDerived.h"
#include "ReconstructionDataFormats/Track.h"
#include "Math/Vector4D.h"
#include "TMath.h"
#include "FemtoDreamCollisionSelection.h"
#include "FemtoDreamTrackSelection.h"
#include "FemtoDreamV0Selection.h"
#include "FemtoUtils.h"

using namespace o2;
using namespace o2::analysis::femtoDream;
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

using FemtoFullTracks = soa::Join<aod::FullTracks,
                                  aod::TracksDCA, aod::TOFSignal,
                                  aod::pidTPCEl, aod::pidTPCMu, aod::pidTPCPi,
                                  aod::pidTPCKa, aod::pidTPCPr, aod::pidTPCDe,
                                  aod::pidTOFEl, aod::pidTOFMu, aod::pidTOFPi,
                                  aod::pidTOFKa, aod::pidTOFPr, aod::pidTOFDe>;

// using FilteredFullV0s = soa::Filtered<aod::V0Datas>; /// predefined Join table for o2::aod::V0s = soa::Join<o2::aod::TransientV0s, o2::aod::StoredV0s> to be used when we add v0Filter
} // namespace o2::aod

/// \todo fix how to pass array to setSelection, getRow() passing a different type!
// static constexpr float arrayV0Sel[3][3] = {{100.f, 100.f, 100.f}, {0.2f, 0.2f, 0.2f}, {100.f, 100.f, 100.f}};
// unsigned int rows = sizeof(arrayV0Sel) / sizeof(arrayV0Sel[0]);
// unsigned int columns = sizeof(arrayV0Sel[0]) / sizeof(arrayV0Sel[0][0]);

template <typename T>
int getRowDaughters(int daughID, T const& vecID)
{
  int rowInPrimaryTrackTableDaugh = -1;
  for (size_t i = 0; i < vecID.size(); i++) {
    if (vecID.at(i) == daughID) {
      rowInPrimaryTrackTableDaugh = i;
      break;
    }
  }
  return rowInPrimaryTrackTableDaugh;
}

struct femtoDreamProducerTask {

  Produces<aod::FemtoDreamCollisions> outputCollision;
  Produces<aod::FemtoDreamParticles> outputParts;
  Produces<aod::FemtoDreamMCParticles> outputPartsMC;
  Produces<aod::FemtoDreamDebugParticles> outputDebugParts;
  Produces<aod::FemtoDreamMCLabels> outputPartsMCLabels;
  Produces<aod::FemtoDreamDebugMCParticles> outputDebugPartsMC;

  Configurable<bool> ConfDebugOutput{"ConfDebugOutput", true, "Debug output"};

  // Choose if filtering or skimming version is run

  Configurable<bool> ConfIsTrigger{"ConfIsTrigger", false, "Store all collisions"};

  // Choose if running on converted data or Run3  / Pilot
  Configurable<bool> ConfIsRun3{"ConfIsRun3", false, "Running on Run3 or pilot"};
  Configurable<bool> ConfIsMC{"ConfIsMC", false, "Running on MC; implemented only for Run3"};

  /// Event cuts
  FemtoDreamCollisionSelection colCuts;
  Configurable<bool> ConfUseTPCmult{"ConfUseTPCmult", false, "Use multiplicity based on the number of tracks with TPC information"};
  Configurable<float> ConfEvtZvtx{"ConfEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};
  Configurable<bool> ConfEvtTriggerCheck{"ConfEvtTriggerCheck", true, "Evt sel: check for trigger"};
  Configurable<int> ConfEvtTriggerSel{"ConfEvtTriggerSel", kINT7, "Evt sel: trigger"};
  Configurable<bool> ConfEvtOfflineCheck{"ConfEvtOfflineCheck", false, "Evt sel: check for offline selection"};

  Configurable<bool> ConfStoreV0{"ConfStoreV0", true, "True: store V0 table"};
  // just sanity check to make sure in case there are problems in conversion or MC production it does not affect results
  Configurable<bool> ConfRejectNotPropagatedTracks{"ConfRejectNotPropagatedTracks", false, "True: reject not propagated tracks"};
  Configurable<bool> ConfRejectITSHitandTOFMissing{"ConfRejectITSHitandTOFMissing", false, "True: reject if neither ITS hit nor TOF timing satisfied"};

  Configurable<int> ConfPDGCodeTrack{"ConfPDGCodeTrack", 2212, "PDG code of the selected track for Monte Carlo truth"};
  FemtoDreamTrackSelection trackCuts;
  Configurable<std::vector<float>> ConfTrkCharge{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kSign, "ConfTrk"), std::vector<float>{-1, 1}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kSign, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkPtmin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kpTMin, "ConfTrk"), std::vector<float>{0.4f, 0.6f, 0.5f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kpTMin, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkPtmax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kpTMax, "ConfTrk"), std::vector<float>{5.4f, 5.6f, 5.5f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kpTMax, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkEta{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kEtaMax, "ConfTrk"), std::vector<float>{0.8f, 0.7f, 0.9f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kEtaMax, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkTPCnclsMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCnClsMin, "ConfTrk"), std::vector<float>{80.f, 70.f, 60.f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCnClsMin, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkTPCfCls{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCfClsMin, "ConfTrk"), std::vector<float>{0.7f, 0.83f, 0.9f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCfClsMin, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkTPCcRowsMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCcRowsMin, "ConfTrk"), std::vector<float>{70.f, 60.f, 80.f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCcRowsMin, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkTPCsCls{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCsClsMax, "ConfTrk"), std::vector<float>{0.1f, 160.f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCsClsMax, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkITSnclsMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kITSnClsMin, "ConfTrk"), std::vector<float>{-1.f, 2.f, 4.f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kITSnClsMin, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkITSnclsIbMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kITSnClsIbMin, "ConfTrk"), std::vector<float>{-1.f, 1.f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kITSnClsIbMin, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkDCAxyMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kDCAxyMax, "ConfTrk"), std::vector<float>{0.1f, 3.5f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kDCAxyMax, "Track selection: ")}; /// here we need an open cut to do the DCA fits later on!
  Configurable<std::vector<float>> ConfTrkDCAzMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kDCAzMax, "ConfTrk"), std::vector<float>{0.2f, 3.5f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kDCAzMax, "Track selection: ")};
  /// \todo Reintegrate PID to the general selection container
  Configurable<std::vector<float>> ConfTrkPIDnSigmaMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kPIDnSigmaMax, "ConfTrk"), std::vector<float>{3.5f, 3.f, 2.5f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kPIDnSigmaMax, "Track selection: ")};
  Configurable<float> ConfTrkPIDnSigmaOffsetTPC{"ConfTrkPIDnSigmaOffsetTPC", 0., "Offset for TPC nSigma because of bad calibration"};
  Configurable<float> ConfTrkPIDnSigmaOffsetTOF{"ConfTrkPIDnSigmaOffsetTOF", 0., "Offset for TOF nSigma because of bad calibration"};
  Configurable<std::vector<int>> ConfTrkTPIDspecies{"ConfTrkTPIDspecies", std::vector<int>{o2::track::PID::Pion, o2::track::PID::Kaon, o2::track::PID::Proton, o2::track::PID::Deuteron}, "Trk sel: Particles species for PID"};

  FemtoDreamV0Selection v0Cuts;
  TrackSelection* o2PhysicsTrackSelection;
  /// \todo Labeled array (see Track-Track task)

  Configurable<std::vector<float>> ConfV0Sign{FemtoDreamV0Selection::getSelectionName(femtoDreamV0Selection::kV0Sign, "ConfV0"), std::vector<float>{-1, 1}, FemtoDreamV0Selection::getSelectionHelper(femtoDreamV0Selection::kV0Sign, "V0 selection: ")};
  Configurable<std::vector<float>> ConfV0PtMin{FemtoDreamV0Selection::getSelectionName(femtoDreamV0Selection::kpTV0Min, "ConfV0"), std::vector<float>{0.3f, 0.4f, 0.5f}, FemtoDreamV0Selection::getSelectionHelper(femtoDreamV0Selection::kpTV0Min, "V0 selection: ")};
  Configurable<std::vector<float>> ConfDCAV0DaughMax{FemtoDreamV0Selection::getSelectionName(femtoDreamV0Selection::kDCAV0DaughMax, "ConfV0"), std::vector<float>{1.2f, 1.5f}, FemtoDreamV0Selection::getSelectionHelper(femtoDreamV0Selection::kDCAV0DaughMax, "V0 selection: ")};
  Configurable<std::vector<float>> ConfCPAV0Min{FemtoDreamV0Selection::getSelectionName(femtoDreamV0Selection::kCPAV0Min, "ConfV0"), std::vector<float>{0.99f, 0.995f}, FemtoDreamV0Selection::getSelectionHelper(femtoDreamV0Selection::kCPAV0Min, "V0 selection: ")};

  Configurable<std::vector<float>> V0TranRadV0Min{FemtoDreamV0Selection::getSelectionName(femtoDreamV0Selection::kTranRadV0Min, "ConfV0"), std::vector<float>{0.2f}, FemtoDreamV0Selection::getSelectionHelper(femtoDreamV0Selection::kTranRadV0Min, "V0 selection: ")};
  Configurable<std::vector<float>> V0TranRadV0Max{FemtoDreamV0Selection::getSelectionName(femtoDreamV0Selection::kTranRadV0Max, "ConfV0"), std::vector<float>{100.f}, FemtoDreamV0Selection::getSelectionHelper(femtoDreamV0Selection::kTranRadV0Max, "V0 selection: ")};
  Configurable<std::vector<float>> V0DecVtxMax{FemtoDreamV0Selection::getSelectionName(femtoDreamV0Selection::kDecVtxMax, "ConfV0"), std::vector<float>{100.f}, FemtoDreamV0Selection::getSelectionHelper(femtoDreamV0Selection::kDecVtxMax, "V0 selection: ")};

  Configurable<std::vector<float>> ConfV0DaughCharge{"ConfV0DaughCharge", std::vector<float>{-1, 1}, "V0 Daugh sel: Charge"};
  Configurable<std::vector<float>> ConfDaughEta{"ConfDaughEta", std::vector<float>{0.8f}, "V0 Daugh sel: max eta"};
  Configurable<std::vector<float>> ConfV0DaughTPCnclsMin{"ConfV0DaughTPCnclsMin", std::vector<float>{80.f, 70.f, 60.f}, "V0 Daugh sel: Min. nCls TPC"};
  Configurable<std::vector<float>> ConfV0DaughDCAMin{"ConfV0DaughDCAMin", std::vector<float>{0.05f, 0.06f}, "V0 Daugh sel:  Max. DCA Daugh to PV (cm)"};
  Configurable<std::vector<float>> ConfV0DaughPIDnSigmaMax{"ConfV0DaughPIDnSigmaMax", std::vector<float>{5.f, 4.f}, "V0 Daugh sel: Max. PID nSigma TPC"};

  Configurable<std::vector<int>> ConfV0DaughTPIDspecies{"ConfV0DaughTPIDspecies", std::vector<int>{o2::track::PID::Pion, o2::track::PID::Proton}, "V0 Daugh sel: Particles species for PID"};

  Configurable<float> ConfInvMassLowLimit{"ConfInvMassLowLimit", 1.05, "Lower limit of the V0 invariant mass"};
  Configurable<float> ConfInvMassUpLimit{"ConfInvMassUpLimit", 1.30, "Upper limit of the V0 invariant mass"};

  Configurable<bool> ConfRejectKaons{"ConfRejectKaons", false, "Switch to reject kaons"};
  Configurable<float> ConfInvKaonMassLowLimit{"ConfInvKaonMassLowLimit", 0.48, "Lower limit of the V0 invariant mass for Kaon rejection"};
  Configurable<float> ConfInvKaonMassUpLimit{"ConfInvKaonMassUpLimit", 0.515, "Upper limit of the V0 invariant mass for Kaon rejection"};

  /// \todo should we add filter on min value pT/eta of V0 and daughters?
  /*Filter v0Filter = (nabs(aod::v0data::x) < V0DecVtxMax.value) &&
                    (nabs(aod::v0data::y) < V0DecVtxMax.value) &&
                    (nabs(aod::v0data::z) < V0DecVtxMax.value);*/
  // (aod::v0data::v0radius > V0TranRadV0Min.value); to be added, not working for now do not know why

  HistogramRegistry qaRegistry{"QAHistos", {}, OutputObjHandlingPolicy::QAObject};

  int mRunNumber;
  float mMagField;
  Service<o2::ccdb::BasicCCDBManager> ccdb; /// Accessing the CCDB

  void init(InitContext&)
  {
    if (doprocessData == false && doprocessMC == false) {
      LOGF(fatal, "Neither processData nor processMC enabled. Please choose one.");
    }
    if (doprocessData == true && doprocessMC == true) {
      LOGF(fatal, "Cannot enable processData and processMC at the same time. Please choose one.");
    }

    colCuts.setCuts(ConfEvtZvtx, ConfEvtTriggerCheck, ConfEvtTriggerSel, ConfEvtOfflineCheck, ConfIsRun3);
    colCuts.init(&qaRegistry);

    Configurable<int> ConfPDGCodeTrack{"ConfPDGCodeTrack", 2212, "PDG code of the selected track for Monte Carlo truth"};
    trackCuts.setSelection(ConfTrkCharge, femtoDreamTrackSelection::kSign, femtoDreamSelection::kEqual);
    trackCuts.setSelection(ConfTrkPtmin, femtoDreamTrackSelection::kpTMin, femtoDreamSelection::kLowerLimit);
    trackCuts.setSelection(ConfTrkPtmax, femtoDreamTrackSelection::kpTMax, femtoDreamSelection::kUpperLimit);
    trackCuts.setSelection(ConfTrkEta, femtoDreamTrackSelection::kEtaMax, femtoDreamSelection::kAbsUpperLimit);
    trackCuts.setSelection(ConfTrkTPCnclsMin, femtoDreamTrackSelection::kTPCnClsMin, femtoDreamSelection::kLowerLimit);
    trackCuts.setSelection(ConfTrkTPCfCls, femtoDreamTrackSelection::kTPCfClsMin, femtoDreamSelection::kLowerLimit);
    trackCuts.setSelection(ConfTrkTPCcRowsMin, femtoDreamTrackSelection::kTPCcRowsMin, femtoDreamSelection::kLowerLimit);
    trackCuts.setSelection(ConfTrkTPCsCls, femtoDreamTrackSelection::kTPCsClsMax, femtoDreamSelection::kUpperLimit);
    trackCuts.setSelection(ConfTrkITSnclsMin, femtoDreamTrackSelection::kITSnClsMin, femtoDreamSelection::kLowerLimit);
    trackCuts.setSelection(ConfTrkITSnclsIbMin, femtoDreamTrackSelection::kITSnClsIbMin, femtoDreamSelection::kLowerLimit);
    trackCuts.setSelection(ConfTrkDCAxyMax, femtoDreamTrackSelection::kDCAxyMax, femtoDreamSelection::kAbsUpperLimit);
    trackCuts.setSelection(ConfTrkDCAzMax, femtoDreamTrackSelection::kDCAzMax, femtoDreamSelection::kAbsUpperLimit);
    trackCuts.setSelection(ConfTrkPIDnSigmaMax, femtoDreamTrackSelection::kPIDnSigmaMax, femtoDreamSelection::kAbsUpperLimit);
    trackCuts.setPIDSpecies(ConfTrkTPIDspecies);
    trackCuts.setnSigmaPIDOffset(ConfTrkPIDnSigmaOffsetTPC, ConfTrkPIDnSigmaOffsetTOF);
    trackCuts.init<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::TrackType::kNoChild, aod::femtodreamparticle::cutContainerType>(&qaRegistry);

    /// \todo fix how to pass array to setSelection, getRow() passing a different type!
    // v0Cuts.setSelection(ConfV0Selection->getRow(0), femtoDreamV0Selection::kDecVtxMax, femtoDreamSelection::kAbsUpperLimit);
    if (ConfStoreV0) {
      v0Cuts.setSelection(ConfV0Sign, femtoDreamV0Selection::kV0Sign, femtoDreamSelection::kEqual);
      v0Cuts.setSelection(ConfV0PtMin, femtoDreamV0Selection::kpTV0Min, femtoDreamSelection::kLowerLimit);
      v0Cuts.setSelection(ConfDCAV0DaughMax, femtoDreamV0Selection::kDCAV0DaughMax, femtoDreamSelection::kUpperLimit);
      v0Cuts.setSelection(ConfCPAV0Min, femtoDreamV0Selection::kCPAV0Min, femtoDreamSelection::kLowerLimit);

      v0Cuts.setChildCuts(femtoDreamV0Selection::kPosTrack, ConfV0DaughCharge, femtoDreamTrackSelection::kSign, femtoDreamSelection::kEqual);
      v0Cuts.setChildCuts(femtoDreamV0Selection::kPosTrack, ConfDaughEta, femtoDreamTrackSelection::kEtaMax, femtoDreamSelection::kAbsUpperLimit);
      v0Cuts.setChildCuts(femtoDreamV0Selection::kPosTrack, ConfV0DaughTPCnclsMin, femtoDreamTrackSelection::kTPCnClsMin, femtoDreamSelection::kLowerLimit);
      v0Cuts.setChildCuts(femtoDreamV0Selection::kPosTrack, ConfV0DaughDCAMin, femtoDreamTrackSelection::kDCAMin, femtoDreamSelection::kAbsLowerLimit);
      v0Cuts.setChildCuts(femtoDreamV0Selection::kPosTrack, ConfV0DaughPIDnSigmaMax, femtoDreamTrackSelection::kPIDnSigmaMax, femtoDreamSelection::kAbsUpperLimit);
      v0Cuts.setChildCuts(femtoDreamV0Selection::kNegTrack, ConfV0DaughCharge, femtoDreamTrackSelection::kSign, femtoDreamSelection::kEqual);
      v0Cuts.setChildCuts(femtoDreamV0Selection::kNegTrack, ConfDaughEta, femtoDreamTrackSelection::kEtaMax, femtoDreamSelection::kAbsUpperLimit);
      v0Cuts.setChildCuts(femtoDreamV0Selection::kNegTrack, ConfV0DaughTPCnclsMin, femtoDreamTrackSelection::kTPCnClsMin, femtoDreamSelection::kLowerLimit);
      v0Cuts.setChildCuts(femtoDreamV0Selection::kNegTrack, ConfV0DaughDCAMin, femtoDreamTrackSelection::kDCAMin, femtoDreamSelection::kAbsLowerLimit);
      v0Cuts.setChildCuts(femtoDreamV0Selection::kNegTrack, ConfV0DaughPIDnSigmaMax, femtoDreamTrackSelection::kPIDnSigmaMax, femtoDreamSelection::kAbsUpperLimit);
      v0Cuts.setChildPIDSpecies(femtoDreamV0Selection::kPosTrack, ConfV0DaughTPIDspecies);
      v0Cuts.setChildPIDSpecies(femtoDreamV0Selection::kNegTrack, ConfV0DaughTPIDspecies);
      v0Cuts.init<aod::femtodreamparticle::ParticleType::kV0, aod::femtodreamparticle::ParticleType::kV0Child, aod::femtodreamparticle::cutContainerType>(&qaRegistry);
      v0Cuts.setInvMassLimits(ConfInvMassLowLimit, ConfInvMassUpLimit);

      v0Cuts.setChildRejectNotPropagatedTracks(femtoDreamV0Selection::kPosTrack, ConfRejectNotPropagatedTracks);
      v0Cuts.setChildRejectNotPropagatedTracks(femtoDreamV0Selection::kNegTrack, ConfRejectNotPropagatedTracks);

      v0Cuts.setnSigmaPIDOffsetTPC(ConfTrkPIDnSigmaOffsetTPC);
      v0Cuts.setChildnSigmaPIDOffset(femtoDreamV0Selection::kPosTrack, ConfTrkPIDnSigmaOffsetTPC, ConfTrkPIDnSigmaOffsetTOF);
      v0Cuts.setChildnSigmaPIDOffset(femtoDreamV0Selection::kNegTrack, ConfTrkPIDnSigmaOffsetTPC, ConfTrkPIDnSigmaOffsetTOF);

      if (ConfRejectKaons) {
        v0Cuts.setKaonInvMassLimits(ConfInvKaonMassLowLimit, ConfInvKaonMassUpLimit);
      }
      if (ConfRejectITSHitandTOFMissing) {
        o2PhysicsTrackSelection = new TrackSelection(getGlobalTrackSelection());
        o2PhysicsTrackSelection->SetRequireHitsInITSLayers(1, {0, 1, 2, 3});
      }
    }

    mRunNumber = 0;
    mMagField = 0.0;
    /// Initializing CCDB
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);
  }

  /// Function to retrieve the nominal mgnetic field in kG (0.1T) and convert it directly to T
  void getMagneticFieldTesla(aod::BCsWithTimestamps::iterator bc)
  {
    // TODO done only once (and not per run). Will be replaced by CCDBConfigurable
    // get magnetic field for run
    if (mRunNumber == bc.runNumber())
      return;
    auto timestamp = bc.timestamp();
    float output = -999;

    if (ConfIsRun3 && !ConfIsMC) {
      static o2::parameters::GRPMagField* grpo = nullptr;
      if (grpo == nullptr) {
        grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", timestamp);
        if (grpo == nullptr) {
          LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
          return;
        }
        LOGF(info, "Retrieved GRP for timestamp %llu with L3 ", timestamp, grpo->getL3Current());
      }
      // taken from GRP onject definition of getNominalL3Field; update later to something smarter (mNominalL3Field = std::lround(5.f * mL3Current / 30000.f);)
      auto NominalL3Field = std::lround(5.f * grpo->getL3Current() / 30000.f);
      output = 0.1 * (NominalL3Field);
    }

    if (!ConfIsRun3 || (ConfIsRun3 && ConfIsMC)) {
      static o2::parameters::GRPObject* grpo = nullptr;
      if (grpo == nullptr) {
        grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>("GLO/GRP/GRP", timestamp);
        if (grpo == nullptr) {
          LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
          return;
        }
        LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
      }
      output = 0.1 * (grpo->getNominalL3Field());
    }
    mMagField = output;
    mRunNumber = bc.runNumber();
  }

  template <bool isTrackOrV0, typename ParticleType>
  void fillDebugParticle(ParticleType const& particle)
  {
    if constexpr (isTrackOrV0) {
      outputDebugParts(particle.sign(),
                       (uint8_t)particle.tpcNClsFound(),
                       particle.tpcNClsFindable(),
                       (uint8_t)particle.tpcNClsCrossedRows(),
                       particle.tpcNClsShared(),
                       particle.tpcInnerParam(),
                       particle.itsNCls(),
                       particle.itsNClsInnerBarrel(),
                       particle.dcaXY(),
                       particle.dcaZ(),
                       particle.tpcSignal(),
                       particle.tpcNSigmaStoreEl(),
                       particle.tpcNSigmaStorePi(),
                       particle.tpcNSigmaStoreKa(),
                       particle.tpcNSigmaStorePr(),
                       particle.tpcNSigmaStoreDe(),
                       particle.tofNSigmaStoreEl(),
                       particle.tofNSigmaStorePi(),
                       particle.tofNSigmaStoreKa(),
                       particle.tofNSigmaStorePr(),
                       particle.tofNSigmaStoreDe(),
                       -999., -999., -999., -999., -999., -999.);
    } else {
      outputDebugParts(-999., -999., -999., -999., -999., -999., -999.,
                       -999., -999., -999., -999., -999., -999., -999.,
                       -999., -999., -999., -999., -999., -999., -999.,
                       particle.dcaV0daughters(),
                       particle.v0radius(),
                       particle.x(),
                       particle.y(),
                       particle.z(),
                       particle.mK0Short()); // QA for v0
    }
  }

  template <typename ParticleType>
  void fillMCParticle(ParticleType const& particle, o2::aod::femtodreamparticle::ParticleType fdparttype)
  {
    if (particle.has_mcParticle()) {
      // get corresponding MC particle and its info
      auto particleMC = particle.mcParticle();
      auto pdgCode = particleMC.pdgCode();

      int particleOrigin = 99;
      auto motherparticleMC = particleMC.template mothers_as<aod::McParticles>().front();

      if (abs(pdgCode) == abs(ConfPDGCodeTrack.value)) {

        if (particleMC.isPhysicalPrimary()) {
          particleOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kPrimary;
        } else if (motherparticleMC.producedByGenerator()) {
          particleOrigin = checkDaughterType(fdparttype, motherparticleMC.pdgCode());
        } else {
          particleOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kMaterial;
        }

      } else {

        particleOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kFake;
      }
      outputPartsMC(particleOrigin, pdgCode, particleMC.pt(), particleMC.eta(), particleMC.phi());
      outputPartsMCLabels(outputPartsMC.lastIndex());
    } else {
      outputPartsMCLabels(-1);
    }
  }

  template <bool isMC, typename V0Type, typename TrackType, typename CollisionType>
  void fillCollisionsAndTracksAndV0(CollisionType const& col, TrackType const& tracks, V0Type const& fullV0s)
  {

    const auto vtxZ = col.posZ();
    const auto spher = colCuts.computeSphericity(col, tracks);
    int mult = 0;
    int multNtr = 0;
    if (ConfIsRun3) {
      mult = col.multFV0M();
      multNtr = col.multNTracksPV();
    } else {
      mult = 0.5 * (col.multFV0M()); /// For benchmarking on Run 2, V0M in FemtoDreamRun2 is defined V0M/2
      multNtr = col.multTracklets();
    }
    if (ConfUseTPCmult) {
      multNtr = col.multTPC();
    }

    // check whether the basic event selection criteria are fulfilled
    // if the basic selection is NOT fulfilled:
    // in case of skimming run - don't store such collisions
    // in case of trigger run - store such collisions but don't store any particle candidates for such collisions
    if (!colCuts.isSelected(col)) {
      if (ConfIsTrigger) {
        outputCollision(vtxZ, mult, multNtr, spher, mMagField);
      }
      return;
    }

    colCuts.fillQA(col);
    outputCollision(vtxZ, mult, multNtr, spher, mMagField);

    int childIDs[2] = {0, 0};    // these IDs are necessary to keep track of the children
    std::vector<int> tmpIDtrack; // this vector keeps track of the matching of the primary track table row <-> aod::track table global index

    for (auto& track : tracks) {
      /// if the most open selection criteria are not fulfilled there is no point looking further at the track
      if (!trackCuts.isSelectedMinimal(track)) {
        continue;
      }
      trackCuts.fillQA<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::TrackType::kNoChild>(track);
      // the bit-wise container of the systematic variations is obtained
      auto cutContainer = trackCuts.getCutContainer<aod::femtodreamparticle::cutContainerType>(track);

      // now the table is filled
      outputParts(outputCollision.lastIndex(),
                  track.pt(),
                  track.eta(),
                  track.phi(),
                  aod::femtodreamparticle::ParticleType::kTrack,
                  cutContainer.at(femtoDreamTrackSelection::TrackContainerPosition::kCuts),
                  cutContainer.at(femtoDreamTrackSelection::TrackContainerPosition::kPID),
                  track.dcaXY(),
                  childIDs, 0, 0);
      tmpIDtrack.push_back(track.globalIndex());
      if (ConfDebugOutput) {
        fillDebugParticle<true>(track);
      }

      if constexpr (isMC) {
        fillMCParticle(track, o2::aod::femtodreamparticle::ParticleType::kTrack);
      }
    }

    if (ConfStoreV0) {
      for (auto& v0 : fullV0s) {
        auto postrack = v0.template posTrack_as<TrackType>();
        auto negtrack = v0.template negTrack_as<TrackType>(); ///\tocheck funnily enough if we apply the filter the sign of Pos and Neg track is always negative
        // const auto dcaXYpos = postrack.dcaXY();
        // const auto dcaZpos = postrack.dcaZ();
        // const auto dcapos = std::sqrt(pow(dcaXYpos, 2.) + pow(dcaZpos, 2.));
        v0Cuts.fillLambdaQA(col, v0, postrack, negtrack);

        if (!v0Cuts.isSelectedMinimal(col, v0, postrack, negtrack)) {
          continue;
        }

        if (ConfRejectITSHitandTOFMissing) {
          // Uncomment only when TOF timing is solved
          // bool itsHit = o2PhysicsTrackSelection->IsSelected(postrack, TrackSelection::TrackCuts::kITSHits);
          // bool itsHit = o2PhysicsTrackSelection->IsSelected(negtrack, TrackSelection::TrackCuts::kITSHits);
        }

        v0Cuts.fillQA<aod::femtodreamparticle::ParticleType::kV0, aod::femtodreamparticle::ParticleType::kV0Child>(col, v0, postrack, negtrack); ///\todo fill QA also for daughters
        auto cutContainerV0 = v0Cuts.getCutContainer<aod::femtodreamparticle::cutContainerType>(col, v0, postrack, negtrack);

        if ((cutContainerV0.at(femtoDreamV0Selection::V0ContainerPosition::kV0) > 0) && (cutContainerV0.at(femtoDreamV0Selection::V0ContainerPosition::kPosCuts) > 0) && (cutContainerV0.at(femtoDreamV0Selection::V0ContainerPosition::kNegCuts) > 0)) {
          int postrackID = v0.posTrackId();
          int rowInPrimaryTrackTablePos = -1;
          rowInPrimaryTrackTablePos = getRowDaughters(postrackID, tmpIDtrack);
          childIDs[0] = rowInPrimaryTrackTablePos;
          childIDs[1] = 0;
          outputParts(outputCollision.lastIndex(), v0.positivept(), v0.positiveeta(), v0.positivephi(), aod::femtodreamparticle::ParticleType::kV0Child, cutContainerV0.at(femtoDreamV0Selection::V0ContainerPosition::kPosCuts), cutContainerV0.at(femtoDreamV0Selection::V0ContainerPosition::kPosPID), 0., childIDs, 0, 0);
          const int rowOfPosTrack = outputParts.lastIndex();
          if constexpr (isMC) {
            fillMCParticle(postrack, o2::aod::femtodreamparticle::ParticleType::kV0Child);
          }
          int negtrackID = v0.negTrackId();
          int rowInPrimaryTrackTableNeg = -1;
          rowInPrimaryTrackTableNeg = getRowDaughters(negtrackID, tmpIDtrack);
          childIDs[0] = 0;
          childIDs[1] = rowInPrimaryTrackTableNeg;
          outputParts(outputCollision.lastIndex(), v0.negativept(), v0.negativeeta(), v0.negativephi(), aod::femtodreamparticle::ParticleType::kV0Child, cutContainerV0.at(femtoDreamV0Selection::V0ContainerPosition::kNegCuts), cutContainerV0.at(femtoDreamV0Selection::V0ContainerPosition::kNegPID), 0., childIDs, 0, 0);
          const int rowOfNegTrack = outputParts.lastIndex();
          if constexpr (isMC) {
            fillMCParticle(negtrack, o2::aod::femtodreamparticle::ParticleType::kV0Child);
          }
          int indexChildID[2] = {rowOfPosTrack, rowOfNegTrack};
          outputParts(outputCollision.lastIndex(), v0.pt(), v0.eta(), v0.phi(), aod::femtodreamparticle::ParticleType::kV0, cutContainerV0.at(femtoDreamV0Selection::V0ContainerPosition::kV0), 0, v0.v0cosPA(col.posX(), col.posY(), col.posZ()), indexChildID, v0.mLambda(), v0.mAntiLambda());
          if (ConfDebugOutput) {
            fillDebugParticle<true>(postrack); // QA for positive daughter
            fillDebugParticle<true>(negtrack); // QA for negative daughter
            fillDebugParticle<false>(v0);      // QA for v0
          }
          if constexpr (isMC) {
            fillMCParticle(v0, o2::aod::femtodreamparticle::ParticleType::kV0);
          }
        }
      }
    }
  }

  void processData(aod::FemtoFullCollision const& col, aod::BCsWithTimestamps const&, aod::FemtoFullTracks const& tracks,
                   o2::aod::V0Datas const& fullV0s) /// \todo with FilteredFullV0s
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    fillCollisionsAndTracksAndV0<false>(col, tracks, fullV0s);
  }
  PROCESS_SWITCH(femtoDreamProducerTask, processData, "Provide experimental data", true);

  void processMC(aod::FemtoFullCollisionMC const& col,
                 aod::BCsWithTimestamps const&,
                 soa::Join<aod::FemtoFullTracks, aod::McTrackLabels> const& tracks,
                 aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles,
                 soa::Join<o2::aod::V0Datas, aod::McV0Labels> const& fullV0s) /// \todo with FilteredFullV0s
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    fillCollisionsAndTracksAndV0<true>(col, tracks, fullV0s);
  }
  PROCESS_SWITCH(femtoDreamProducerTask, processMC, "Provide MC data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<femtoDreamProducerTask>(cfgc)};
  return workflow;
}
