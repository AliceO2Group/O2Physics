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
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de

#include "FemtoDreamCollisionSelection.h"
#include "FemtoDreamTrackSelection.h"
#include "FemtoDreamV0Selection.h"
#include "PWGCF/DataModel/FemtoDerived.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/PID/PIDResponse.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/StrangenessTables.h"
#include "Math/Vector4D.h"
#include "TMath.h"

using namespace o2;
using namespace o2::analysis::femtoDream;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{

using FemtoFullCollision = soa::Join<aod::Collisions,
                                     aod::EvSels,
                                     aod::Mults>::iterator;
using FemtoFullTracks = soa::Join<aod::FullTracks,
                                  aod::TracksExtended, aod::TOFSignal,
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
  Produces<aod::FemtoDreamDebugParticles> outputDebugParts;

  Configurable<bool> ConfDebugOutput{"ConfDebugOutput", true, "Debug output"};

  // Choose if filtering or skimming version is run

  Configurable<bool> ConfIsTrigger{"ConfIsTrigger", false, "Store all collisions"};

  // Choose if running on converted data or pilot beam
  Configurable<bool> ConfIsRun3{"ConfIsRun3", false, "Running on Pilot beam"};

  /// Event cuts
  FemtoDreamCollisionSelection colCuts;
  Configurable<float> ConfEvtZvtx{"ConfEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};
  Configurable<bool> ConfEvtTriggerCheck{"ConfEvtTriggerCheck", true, "Evt sel: check for trigger"};
  Configurable<int> ConfEvtTriggerSel{"ConfEvtTriggerSel", kINT7, "Evt sel: trigger"};
  Configurable<bool> ConfEvtOfflineCheck{"ConfEvtOfflineCheck", false, "Evt sel: check for offline selection"};

  Configurable<bool> ConfStoreV0{"ConfStoreV0", true, "True: store V0 table"};

  FemtoDreamTrackSelection trackCuts;
  Configurable<std::vector<float>> ConfTrkCharge{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kSign, "ConfTrk"), std::vector<float>{-1, 1}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kSign, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkPtmin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kpTMin, "ConfTrk"), std::vector<float>{0.4f, 0.6f, 0.5f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kpTMin, "Track selection: ")};
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
  Configurable<std::vector<int>> ConfTrkTPIDspecies{"ConfTrkTPIDspecies", std::vector<int>{o2::track::PID::Pion, o2::track::PID::Kaon, o2::track::PID::Proton, o2::track::PID::Deuteron}, "Trk sel: Particles species for PID"};

  FemtoDreamV0Selection v0Cuts;
  /// \todo Labeled array (see Track-Track task)

  Configurable<std::vector<float>> ConfV0Sign{FemtoDreamV0Selection::getSelectionName(femtoDreamV0Selection::kV0Sign, "ConfV0"), std::vector<float>{-1, 1}, FemtoDreamV0Selection::getSelectionHelper(femtoDreamV0Selection::kV0Sign, "V0 selection: ")};
  Configurable<std::vector<float>> ConfV0PtMin{FemtoDreamV0Selection::getSelectionName(femtoDreamV0Selection::kpTV0Min, "ConfV0"), std::vector<float>{0.3f, 0.4f, 0.5f}, FemtoDreamV0Selection::getSelectionHelper(femtoDreamV0Selection::kpTV0Min, "V0 selection: ")};
  Configurable<std::vector<float>> ConfDCAV0DaughMax{FemtoDreamV0Selection::getSelectionName(femtoDreamV0Selection::kDCAV0DaughMax, "ConfV0"), std::vector<float>{1.2f, 1.5f}, FemtoDreamV0Selection::getSelectionHelper(femtoDreamV0Selection::kDCAV0DaughMax, "V0 selection: ")};
  Configurable<std::vector<float>> ConfCPAV0Min{FemtoDreamV0Selection::getSelectionName(femtoDreamV0Selection::kCPAV0Min, "ConfV0"), std::vector<float>{0.99f, 0.995f}, FemtoDreamV0Selection::getSelectionHelper(femtoDreamV0Selection::kCPAV0Min, "V0 selection: ")};

  MutableConfigurable<float> V0TranRadV0Min{FemtoDreamV0Selection::getSelectionName(femtoDreamV0Selection::kTranRadV0Min, "ConfV0"), 0.2f, FemtoDreamV0Selection::getSelectionHelper(femtoDreamV0Selection::kTranRadV0Min, "V0 selection: ")}; // QUESTION: why mutable?
  MutableConfigurable<float> V0TranRadV0Max{FemtoDreamV0Selection::getSelectionName(femtoDreamV0Selection::kTranRadV0Max, "ConfV0"), 100.f, FemtoDreamV0Selection::getSelectionHelper(femtoDreamV0Selection::kTranRadV0Max, "V0 selection: ")};
  MutableConfigurable<float> V0DecVtxMax{FemtoDreamV0Selection::getSelectionName(femtoDreamV0Selection::kDecVtxMax, "ConfV0"), 100.f, FemtoDreamV0Selection::getSelectionHelper(femtoDreamV0Selection::kDecVtxMax, "V0 selection: ")};

  Configurable<std::vector<float>> ConfV0DaughTPCnclsMin{"ConfV0DaughTPCnclsMin", std::vector<float>{80.f, 70.f, 60.f}, "V0 Daugh sel: Min. nCls TPC"};
  Configurable<std::vector<float>> ConfV0DaughDCAMin{"ConfV0DaughDCAMin", std::vector<float>{0.05f, 0.06f}, "V0 Daugh sel:  Max. DCA Daugh to PV (cm)"};
  Configurable<std::vector<float>> ConfV0DaughPIDnSigmaMax{"ConfV0DaughPIDnSigmaMax", std::vector<float>{5.f, 4.f}, "V0 Daugh sel: Max. PID nSigma TPC"};

  Configurable<float> ConfInvMassLowLimit{"ConfInvMassLowLimit", 1.1, "Lower limit of the V0 invariant mass"};
  Configurable<float> ConfInvMassUpLimit{"ConfInvMassUpLimit", 1.13, "Upper limit of the V0 invariant mass"};

  /// \todo should we add filter on min value pT/eta of V0 and daughters?
  /*Filter v0Filter = (nabs(aod::v0data::x) < V0DecVtxMax.value) &&
                    (nabs(aod::v0data::y) < V0DecVtxMax.value) &&
                    (nabs(aod::v0data::z) < V0DecVtxMax.value);*/
  // (aod::v0data::v0radius > V0TranRadV0Min.value); to be added, not working for now do not know why

  HistogramRegistry qaRegistry{"QAHistos", {}, OutputObjHandlingPolicy::QAObject};

  void init(InitContext&)
  {
    colCuts.setCuts(ConfEvtZvtx, ConfEvtTriggerCheck, ConfEvtTriggerSel, ConfEvtOfflineCheck, ConfIsRun3);
    colCuts.init(&qaRegistry);

    trackCuts.setSelection(ConfTrkCharge, femtoDreamTrackSelection::kSign, femtoDreamSelection::kEqual);
    trackCuts.setSelection(ConfTrkPtmin, femtoDreamTrackSelection::kpTMin, femtoDreamSelection::kLowerLimit);
    trackCuts.setSelection(ConfTrkEta, femtoDreamTrackSelection::kEtaMax, femtoDreamSelection::kAbsUpperLimit);
    trackCuts.setSelection(ConfTrkTPCnclsMin, femtoDreamTrackSelection::kTPCnClsMin, femtoDreamSelection::kLowerLimit);
    trackCuts.setSelection(ConfTrkTPCfCls, femtoDreamTrackSelection::kTPCfClsMin, femtoDreamSelection::kLowerLimit);
    trackCuts.setSelection(ConfTrkTPCcRowsMin, femtoDreamTrackSelection::kTPCcRowsMin, femtoDreamSelection::kLowerLimit);
    trackCuts.setSelection(ConfTrkTPCsCls, femtoDreamTrackSelection::kTPCsClsMax, femtoDreamSelection::kUpperLimit);
    if (ConfIsRun3) {
      trackCuts.setSelection(ConfTrkITSnclsMin, femtoDreamTrackSelection::kITSnClsMin, femtoDreamSelection::kLowerLimit);
      trackCuts.setSelection(ConfTrkITSnclsIbMin, femtoDreamTrackSelection::kITSnClsIbMin, femtoDreamSelection::kLowerLimit);
    }
    trackCuts.setSelection(ConfTrkDCAxyMax, femtoDreamTrackSelection::kDCAxyMax, femtoDreamSelection::kAbsUpperLimit);
    trackCuts.setSelection(ConfTrkDCAzMax, femtoDreamTrackSelection::kDCAzMax, femtoDreamSelection::kAbsUpperLimit);
    trackCuts.setSelection(ConfTrkPIDnSigmaMax, femtoDreamTrackSelection::kPIDnSigmaMax, femtoDreamSelection::kAbsUpperLimit);
    trackCuts.setPIDSpecies(ConfTrkTPIDspecies);
    trackCuts.init<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::cutContainerType>(&qaRegistry);

    /// \todo fix how to pass array to setSelection, getRow() passing a different type!
    // v0Cuts.setSelection(ConfV0Selection->getRow(0), femtoDreamV0Selection::kDecVtxMax, femtoDreamSelection::kAbsUpperLimit);
    if (ConfStoreV0) {
      v0Cuts.setSelection(ConfV0Sign, femtoDreamV0Selection::kV0Sign, femtoDreamSelection::kEqual);
      v0Cuts.setSelection(ConfV0PtMin, femtoDreamV0Selection::kpTV0Min, femtoDreamSelection::kLowerLimit);
      v0Cuts.setSelection(ConfDCAV0DaughMax, femtoDreamV0Selection::kDCAV0DaughMax, femtoDreamSelection::kUpperLimit);
      v0Cuts.setSelection(ConfCPAV0Min, femtoDreamV0Selection::kCPAV0Min, femtoDreamSelection::kLowerLimit);

      v0Cuts.setChildCuts(femtoDreamV0Selection::kPosTrack, ConfTrkCharge, femtoDreamTrackSelection::kSign, femtoDreamSelection::kEqual);
      v0Cuts.setChildCuts(femtoDreamV0Selection::kPosTrack, ConfTrkEta, femtoDreamTrackSelection::kEtaMax, femtoDreamSelection::kAbsUpperLimit);
      v0Cuts.setChildCuts(femtoDreamV0Selection::kPosTrack, ConfV0DaughTPCnclsMin, femtoDreamTrackSelection::kTPCnClsMin, femtoDreamSelection::kLowerLimit);
      // v0Cuts.setChildCuts(femtoDreamV0Selection::kPosTrack, ConfV0DaughDCAMin, femtoDreamTrackSelection::kDCAMin, femtoDreamSelection::kAbsLowerLimit);
      v0Cuts.setChildCuts(femtoDreamV0Selection::kPosTrack, ConfV0DaughPIDnSigmaMax, femtoDreamTrackSelection::kPIDnSigmaMax, femtoDreamSelection::kAbsUpperLimit);
      v0Cuts.setChildCuts(femtoDreamV0Selection::kNegTrack, ConfTrkCharge, femtoDreamTrackSelection::kSign, femtoDreamSelection::kEqual);
      v0Cuts.setChildCuts(femtoDreamV0Selection::kNegTrack, ConfTrkEta, femtoDreamTrackSelection::kEtaMax, femtoDreamSelection::kAbsUpperLimit);
      v0Cuts.setChildCuts(femtoDreamV0Selection::kNegTrack, ConfV0DaughTPCnclsMin, femtoDreamTrackSelection::kTPCnClsMin, femtoDreamSelection::kLowerLimit);
      // v0Cuts.setChildCuts(femtoDreamV0Selection::kNegTrack, ConfV0DaughDCAMin, femtoDreamTrackSelection::kDCAMin, femtoDreamSelection::kAbsLowerLimit);
      v0Cuts.setChildCuts(femtoDreamV0Selection::kNegTrack, ConfV0DaughPIDnSigmaMax, femtoDreamTrackSelection::kPIDnSigmaMax, femtoDreamSelection::kAbsUpperLimit);
      v0Cuts.setChildPIDSpecies(femtoDreamV0Selection::kPosTrack, ConfTrkTPIDspecies);
      v0Cuts.setChildPIDSpecies(femtoDreamV0Selection::kNegTrack, ConfTrkTPIDspecies);
      v0Cuts.init<aod::femtodreamparticle::ParticleType::kV0, aod::femtodreamparticle::ParticleType::kV0Child, aod::femtodreamparticle::cutContainerType>(&qaRegistry);
      v0Cuts.setInvMassLimits(ConfInvMassLowLimit, ConfInvMassUpLimit);
    }
  }

  void process(aod::FemtoFullCollision const& col, aod::BCsWithTimestamps const&, aod::FemtoFullTracks const& tracks,
               o2::aod::V0Datas const& fullV0s) /// \todo with FilteredFullV0s
  {
    auto bc = col.bc_as<aod::BCsWithTimestamps>(); /// adding timestamp to access magnetic field later
    /// First thing to do is to check whether the basic event selection criteria are fulfilled
    // If the basic selection is NOT fullfilled:
    // in case of skimming run - don't store such collisions
    // in case of trigger run - store such collisions but don't store any particle candidates for such collisions
    if (!colCuts.isSelected(col)) {
      if (ConfIsTrigger) {
        outputCollision(col.posZ(), col.multFV0M(), colCuts.computeSphericity(col, tracks), bc.timestamp());
      }
      return;
    }

    const auto vtxZ = col.posZ();
    const auto mult = col.multFV0M();
    const auto spher = colCuts.computeSphericity(col, tracks);
    colCuts.fillQA(col);

    // now the table is filled
    if (ConfIsRun3) {
      outputCollision(vtxZ, col.multFT0M(), spher, bc.timestamp());
    } else {
      outputCollision(vtxZ, mult, spher, bc.timestamp());
    }

    int childIDs[2] = {0, 0};    // these IDs are necessary to keep track of the children
    std::vector<int> tmpIDtrack; // this vector keeps track of the matching of the primary track table row <-> aod::track table global index

    for (auto& track : tracks) {
      /// if the most open selection criteria are not fulfilled there is no point looking further at the track
      if (!trackCuts.isSelectedMinimal(track)) {
        continue;
      }
      trackCuts.fillQA<aod::femtodreamparticle::ParticleType::kTrack>(track);
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
        outputDebugParts(track.sign(),
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
                         -999.);
      }
    }

    if (ConfStoreV0) {
      for (auto& v0 : fullV0s) {
        auto postrack = v0.posTrack_as<aod::FemtoFullTracks>();
        auto negtrack = v0.negTrack_as<aod::FemtoFullTracks>(); ///\tocheck funnily enough if we apply the filter the sign of Pos and Neg track is always negative
        // const auto dcaXYpos = postrack.dcaXY();
        // const auto dcaZpos = postrack.dcaZ();
        // const auto dcapos = std::sqrt(pow(dcaXYpos, 2.) + pow(dcaZpos, 2.));
        if (!v0Cuts.isSelectedMinimal(col, v0, postrack, negtrack)) {
          continue;
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
          int negtrackID = v0.negTrackId();
          int rowInPrimaryTrackTableNeg = -1;
          rowInPrimaryTrackTableNeg = getRowDaughters(negtrackID, tmpIDtrack);
          childIDs[0] = 0;
          childIDs[1] = rowInPrimaryTrackTableNeg;
          outputParts(outputCollision.lastIndex(), v0.negativept(), v0.negativeeta(), v0.negativephi(), aod::femtodreamparticle::ParticleType::kV0Child, cutContainerV0.at(femtoDreamV0Selection::V0ContainerPosition::kNegCuts), cutContainerV0.at(femtoDreamV0Selection::V0ContainerPosition::kNegPID), 0., childIDs, 0, 0);
          const int rowOfNegTrack = outputParts.lastIndex();
          int indexChildID[2] = {rowOfPosTrack, rowOfNegTrack};
          outputParts(outputCollision.lastIndex(), v0.pt(), v0.eta(), v0.phi(), aod::femtodreamparticle::ParticleType::kV0, cutContainerV0.at(femtoDreamV0Selection::V0ContainerPosition::kV0), 0, v0.v0cosPA(col.posX(), col.posY(), col.posZ()), indexChildID, v0.mLambda(), v0.mAntiLambda());
          if (ConfDebugOutput) {
            outputDebugParts(postrack.sign(),
                             (uint8_t)postrack.tpcNClsFound(),
                             postrack.tpcNClsFindable(),
                             (uint8_t)postrack.tpcNClsCrossedRows(),
                             postrack.tpcNClsShared(),
                             postrack.tpcInnerParam(),
                             postrack.itsNCls(),
                             postrack.itsNClsInnerBarrel(),
                             postrack.dcaXY(),
                             postrack.dcaZ(),
                             postrack.tpcSignal(),
                             postrack.tpcNSigmaStoreEl(),
                             postrack.tpcNSigmaStorePi(),
                             postrack.tpcNSigmaStoreKa(),
                             postrack.tpcNSigmaStorePr(),
                             postrack.tpcNSigmaStoreDe(),
                             postrack.tofNSigmaStoreEl(),
                             postrack.tofNSigmaStorePi(),
                             postrack.tofNSigmaStoreKa(),
                             postrack.tofNSigmaStorePr(),
                             postrack.tofNSigmaStoreDe(),
                             -999.,
                             -999.,
                             -999.,
                             -999.,
                             -999.); // QA for positive daughter
            outputDebugParts(negtrack.sign(),
                             (uint8_t)negtrack.tpcNClsFound(),
                             negtrack.tpcNClsFindable(),
                             (uint8_t)negtrack.tpcNClsCrossedRows(),
                             negtrack.tpcNClsShared(),
                             negtrack.tpcInnerParam(),
                             negtrack.itsNCls(),
                             negtrack.itsNClsInnerBarrel(),
                             negtrack.dcaXY(),
                             negtrack.dcaZ(),
                             negtrack.tpcSignal(),
                             negtrack.tpcNSigmaStoreEl(),
                             negtrack.tpcNSigmaStorePi(),
                             negtrack.tpcNSigmaStoreKa(),
                             negtrack.tpcNSigmaStorePr(),
                             negtrack.tpcNSigmaStoreDe(),
                             negtrack.tofNSigmaStoreEl(),
                             negtrack.tofNSigmaStorePi(),
                             negtrack.tofNSigmaStoreKa(),
                             negtrack.tofNSigmaStorePr(),
                             negtrack.tofNSigmaStoreDe(),
                             -999.,
                             -999.,
                             -999.,
                             -999.,
                             -999.); // QA for negative daughter
            outputDebugParts(-999.,
                             -999.,
                             -999.,
                             -999.,
                             -999.,
                             -999.,
                             -999.,
                             -999.,
                             -999.,
                             -999.,
                             -999.,
                             -999.,
                             -999.,
                             -999.,
                             -999.,
                             -999.,
                             -999.,
                             -999.,
                             -999.,
                             -999.,
                             -999.,
                             v0.dcaV0daughters(),
                             v0.v0radius(),
                             v0.x(),
                             v0.y(),
                             v0.z()); // QA for V0
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<femtoDreamProducerTask>(cfgc)};
  return workflow;
}
