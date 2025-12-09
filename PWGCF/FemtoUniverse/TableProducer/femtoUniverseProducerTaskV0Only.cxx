// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
/// \file femtoUniverseProducerTaskV0Only.cxx
/// \brief Tasks that produces the track tables used for the pairing
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Zuzanna Chochulska, WUT Warsaw & CTU Prague, zchochul@cern.ch

#include "PWGCF/FemtoUniverse/Core/FemtoUniverseCollisionSelection.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseTrackSelection.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseV0Selection.h"
#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "DataFormatsParameters/GRPMagField.h"
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

#include <vector>

using namespace o2;
using namespace o2::analysis::femto_universe;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{

using FemtoFullCollision =
  soa::Join<aod::Collisions, aod::EvSels, aod::Mults>::iterator;
using FemtoFullTracks =
  soa::Join<aod::FullTracks, aod::TracksDCA, aod::TOFSignal, aod::pidTPCEl,
            aod::pidTPCMu, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr,
            aod::pidTPCDe, aod::pidTOFEl, aod::pidTOFMu, aod::pidTOFPi,
            aod::pidTOFKa, aod::pidTOFPr, aod::pidTOFDe>;
// using FilteredFullV0s = soa::Filtered<aod::V0Datas>; /// predefined Join
// table for o2::aod::V0s = soa::Join<o2::aod::TransientV0s, o2::aod::StoredV0s>
// to be used when we add v0Filter
} // namespace o2::aod

/// \todo fix how to pass array to setSelection, getRow() passing a different
/// type!
// static constexpr float arrayV0Sel[3][3] = {{100.f, 100.f, 100.f}, {0.2f,
// 0.2f, 0.2f}, {100.f, 100.f, 100.f}}; unsigned int rows = sizeof(arrayV0Sel) /
// sizeof(arrayV0Sel[0]); unsigned int columns = sizeof(arrayV0Sel[0]) /
// sizeof(arrayV0Sel[0][0]);

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

struct femtoUniverseProducerTaskV0Only {

  Produces<aod::FdCollisions> outputCollision;
  Produces<aod::FDParticles> outputParts;
  Produces<aod::FDExtParticles> outputDebugParts;

  Configurable<bool> ConfDebugOutput{"ConfDebugOutput", true, "Debug output"};

  // Choose if filtering or skimming version is run

  Configurable<bool> ConfIsTrigger{"ConfIsTrigger", false, "Store all collisions"};

  // Choose if running on converted data or Run3  / Pilot
  Configurable<bool> ConfIsRun3{"ConfIsRun3", false, "Running on Run3 or pilot"};
  Configurable<bool> ConfIsMC{"ConfIsMC", false, "Running on MC; implemented only for Run3"};

  /// Event cuts
  FemtoUniverseCollisionSelection colCuts;
  Configurable<bool> ConfUseTPCmult{"ConfUseTPCmult", false, "Use multiplicity based on the number of tracks with TPC information"};
  Configurable<float> ConfEvtZvtx{"ConfEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};
  Configurable<bool> ConfEvtTriggerCheck{"ConfEvtTriggerCheck", true, "Evt sel: check for trigger"};
  Configurable<int> ConfEvtTriggerSel{"ConfEvtTriggerSel", kINT7, "Evt sel: trigger"};
  Configurable<bool> ConfEvtOfflineCheck{"ConfEvtOfflineCheck", false, "Evt sel: check for offline selection"};
  Configurable<float> ConfCentFT0Min{"ConfCentFT0Min", 0.f, "Min CentFT0 value for centrality selection"};
  Configurable<float> ConfCentFT0Max{"ConfCentFT0Max", 200.f, "Max CentFT0 value for centrality selection"};

  Configurable<bool> ConfStoreV0{"ConfStoreV0", true, "True: store V0 table"};
  // just sanity check to make sure in case there are problems in conversion or
  // MC production it does not affect results
  Configurable<bool> ConfRejectNotPropagatedTracks{"ConfRejectNotPropagatedTracks", false, "True: reject not propagated tracks"};
  FemtoUniverseV0Selection v0Cuts;
  /// \todo Labeled array (see Track-Track task)

  Configurable<std::vector<float>> ConfV0Sign{
    FemtoUniverseV0Selection::getSelectionName(femto_universe_v0_selection::kV0Sign,
                                               "ConfV0"),
    std::vector<float>{-1, 1},
    FemtoUniverseV0Selection::getSelectionHelper(femto_universe_v0_selection::kV0Sign,
                                                 "V0 selection: ")};
  Configurable<std::vector<float>> ConfV0PtMin{
    FemtoUniverseV0Selection::getSelectionName(femto_universe_v0_selection::kV0pTMin,
                                               "ConfV0"),
    std::vector<float>{0.3f},
    FemtoUniverseV0Selection::getSelectionHelper(femto_universe_v0_selection::kV0pTMin,
                                                 "V0 selection: ")};
  Configurable<std::vector<float>> ConfV0PtMax{
    FemtoUniverseV0Selection::getSelectionName(femto_universe_v0_selection::kV0pTMax,
                                               "ConfV0"),
    std::vector<float>{6.f},
    FemtoUniverseV0Selection::getSelectionHelper(femto_universe_v0_selection::kV0pTMax,
                                                 "V0 selection: ")};
  Configurable<std::vector<float>> ConfV0EtaMax{
    FemtoUniverseV0Selection::getSelectionName(femto_universe_v0_selection::kV0etaMax,
                                               "ConfV0"),
    std::vector<float>{6.f},
    FemtoUniverseV0Selection::getSelectionHelper(femto_universe_v0_selection::kV0etaMax,
                                                 "V0 selection: ")};
  Configurable<std::vector<float>> ConfDCAV0DaughMax{
    FemtoUniverseV0Selection::getSelectionName(
      femto_universe_v0_selection::kV0DCADaughMax, "ConfV0"),
    std::vector<float>{1.5f},
    FemtoUniverseV0Selection::getSelectionHelper(
      femto_universe_v0_selection::kV0DCADaughMax, "V0 selection: ")};
  Configurable<std::vector<float>> ConfCPAV0Min{
    FemtoUniverseV0Selection::getSelectionName(femto_universe_v0_selection::kV0CPAMin,
                                               "ConfV0"),
    std::vector<float>{0.99f},
    FemtoUniverseV0Selection::getSelectionHelper(
      femto_universe_v0_selection::kV0CPAMin, "V0 selection: ")};

  Configurable<std::vector<float>> V0TranRadV0Min{
    FemtoUniverseV0Selection::getSelectionName(
      femto_universe_v0_selection::kV0TranRadMin, "ConfV0"),
    std::vector<float>{0.2f},
    FemtoUniverseV0Selection::getSelectionHelper(
      femto_universe_v0_selection::kV0TranRadMin, "V0 selection: ")};
  Configurable<std::vector<float>> V0TranRadV0Max{
    FemtoUniverseV0Selection::getSelectionName(
      femto_universe_v0_selection::kV0TranRadMax, "ConfV0"),
    std::vector<float>{100.f},
    FemtoUniverseV0Selection::getSelectionHelper(
      femto_universe_v0_selection::kV0TranRadMax, "V0 selection: ")};
  Configurable<std::vector<float>> V0DecVtxMax{
    FemtoUniverseV0Selection::getSelectionName(
      femto_universe_v0_selection::kV0DecVtxMax, "ConfV0"),
    std::vector<float>{100.f},
    FemtoUniverseV0Selection::getSelectionHelper(
      femto_universe_v0_selection::kV0DecVtxMax, "V0 selection: ")};

  Configurable<std::vector<float>> ConfV0DaughCharge{
    "ConfV0DaughCharge", std::vector<float>{-1, 1}, "V0 Daugh sel: Charge"};
  Configurable<std::vector<float>> ConfDaughEta{
    "ConfDaughEta", std::vector<float>{0.8f}, "V0 Daugh sel: max eta"};
  Configurable<std::vector<float>> ConfV0DaughTPCnclsMin{
    "ConfV0DaughTPCnclsMin", std::vector<float>{70.f},
    "V0 Daugh sel: Min. nCls TPC"};
  Configurable<std::vector<float>> ConfV0DaughDCAMin{
    "ConfV0DaughDCAMin", std::vector<float>{0.05f},
    "V0 Daugh sel:  Max. DCA Daugh to PV (cm)"};
  Configurable<std::vector<float>> ConfV0DaughPIDnSigmaMax{
    "ConfV0DaughPIDnSigmaMax", std::vector<float>{5.f},
    "V0 Daugh sel: Max. PID nSigma TPC"};

  Configurable<std::vector<int>> ConfV0DaughTPIDspecies{
    "ConfV0DaughTPIDspecies",
    std::vector<int>{o2::track::PID::Pion, o2::track::PID::Kaon,
                     o2::track::PID::Proton},
    "V0 Daugh sel: Particles species for PID"};

  Configurable<float> ConfInvMassLowLimit{
    "ConfInvMassLowLimit", 1.05, "Lower limit of the V0 invariant mass"};
  Configurable<float> ConfInvMassUpLimit{
    "ConfInvMassUpLimit", 1.30, "Upper limit of the V0 invariant mass"};

  Configurable<bool> ConfRejectKaons{"ConfRejectKaons", false,
                                     "Switch to reject kaons"};
  Configurable<float> ConfInvKaonMassLowLimit{
    "ConfInvKaonMassLowLimit", 0.48,
    "Lower limit of the V0 invariant mass for Kaon rejection"};
  Configurable<float> ConfInvKaonMassUpLimit{
    "ConfInvKaonMassUpLimit", 0.515,
    "Upper limit of the V0 invariant mass for Kaon rejection"};

  /// \todo should we add filter on min value pT/eta of V0 and daughters?
  /*Filter v0Filter = (nabs(aod::v0data::x) < V0DecVtxMax.value) &&
                    (nabs(aod::v0data::y) < V0DecVtxMax.value) &&
                    (nabs(aod::v0data::z) < V0DecVtxMax.value);*/
  // (aod::v0data::v0radius > V0TranRadV0Min.value); to be added, not working
  // for now do not know why

  HistogramRegistry qaRegistry{
    "QAHistos",
    {},
    OutputObjHandlingPolicy::QAObject};

  int mRunNumber;
  float mMagField;
  Service<o2::ccdb::BasicCCDBManager> ccdb; /// Accessing the CCDB

  void init(InitContext&)
  {
    colCuts.setCuts(ConfEvtZvtx, ConfEvtTriggerCheck, ConfEvtTriggerSel,
                    ConfEvtOfflineCheck, ConfIsRun3, ConfCentFT0Min, ConfCentFT0Max);
    colCuts.init(&qaRegistry);

    /// \todo fix how to pass array to setSelection, getRow() passing a
    /// different type!
    // v0Cuts.setSelection(ConfV0Selection->getRow(0),
    // femto_universe_v0_selection::kDecVtxMax, femto_universe_selection::kAbsUpperLimit);
    if (ConfStoreV0) {
      v0Cuts.setSelection(ConfV0Sign, femto_universe_v0_selection::kV0Sign,
                          femto_universe_selection::kEqual);
      v0Cuts.setSelection(ConfV0PtMin, femto_universe_v0_selection::kV0pTMin,
                          femto_universe_selection::kLowerLimit);
      v0Cuts.setSelection(ConfV0PtMax, femto_universe_v0_selection::kV0pTMax,
                          femto_universe_selection::kUpperLimit);
      v0Cuts.setSelection(ConfV0EtaMax, femto_universe_v0_selection::kV0etaMax,
                          femto_universe_selection::kUpperLimit);
      v0Cuts.setSelection(ConfDCAV0DaughMax,
                          femto_universe_v0_selection::kV0DCADaughMax,
                          femto_universe_selection::kUpperLimit);
      v0Cuts.setSelection(ConfCPAV0Min, femto_universe_v0_selection::kV0CPAMin,
                          femto_universe_selection::kLowerLimit);

      v0Cuts.setChildCuts(femto_universe_v0_selection::kPosTrack, ConfV0DaughCharge,
                          femto_universe_track_selection::kSign,
                          femto_universe_selection::kEqual);
      v0Cuts.setChildCuts(femto_universe_v0_selection::kPosTrack, ConfDaughEta,
                          femto_universe_track_selection::kEtaMax,
                          femto_universe_selection::kAbsUpperLimit);
      v0Cuts.setChildCuts(femto_universe_v0_selection::kPosTrack,
                          ConfV0DaughTPCnclsMin,
                          femto_universe_track_selection::kTPCnClsMin,
                          femto_universe_selection::kLowerLimit);
      v0Cuts.setChildCuts(femto_universe_v0_selection::kPosTrack, ConfV0DaughDCAMin,
                          femto_universe_track_selection::kDCAMin,
                          femto_universe_selection::kAbsLowerLimit);
      v0Cuts.setChildCuts(femto_universe_v0_selection::kPosTrack,
                          ConfV0DaughPIDnSigmaMax,
                          femto_universe_track_selection::kPIDnSigmaMax,
                          femto_universe_selection::kAbsUpperLimit);
      v0Cuts.setChildCuts(femto_universe_v0_selection::kNegTrack, ConfV0DaughCharge,
                          femto_universe_track_selection::kSign,
                          femto_universe_selection::kEqual);
      v0Cuts.setChildCuts(femto_universe_v0_selection::kNegTrack, ConfDaughEta,
                          femto_universe_track_selection::kEtaMax,
                          femto_universe_selection::kAbsUpperLimit);
      v0Cuts.setChildCuts(femto_universe_v0_selection::kNegTrack,
                          ConfV0DaughTPCnclsMin,
                          femto_universe_track_selection::kTPCnClsMin,
                          femto_universe_selection::kLowerLimit);
      v0Cuts.setChildCuts(femto_universe_v0_selection::kNegTrack, ConfV0DaughDCAMin,
                          femto_universe_track_selection::kDCAMin,
                          femto_universe_selection::kAbsLowerLimit);
      v0Cuts.setChildCuts(femto_universe_v0_selection::kNegTrack,
                          ConfV0DaughPIDnSigmaMax,
                          femto_universe_track_selection::kPIDnSigmaMax,
                          femto_universe_selection::kAbsUpperLimit);
      v0Cuts.setChildPIDSpecies(femto_universe_v0_selection::kPosTrack,
                                ConfV0DaughTPIDspecies);
      v0Cuts.setChildPIDSpecies(femto_universe_v0_selection::kNegTrack,
                                ConfV0DaughTPIDspecies);
      v0Cuts.init<aod::femtouniverseparticle::ParticleType::kV0,
                  aod::femtouniverseparticle::ParticleType::kV0Child,
                  aod::femtouniverseparticle::CutContainerType>(&qaRegistry);
      v0Cuts.setInvMassLimits(ConfInvMassLowLimit, ConfInvMassUpLimit);
      v0Cuts.setChildRejectNotPropagatedTracks(femto_universe_v0_selection::kPosTrack,
                                               ConfRejectNotPropagatedTracks);
      v0Cuts.setChildRejectNotPropagatedTracks(femto_universe_v0_selection::kNegTrack,
                                               ConfRejectNotPropagatedTracks);
      if (ConfRejectKaons) {
        v0Cuts.setKaonInvMassLimits(ConfInvKaonMassLowLimit,
                                    ConfInvKaonMassUpLimit);
      }
    }
    mRunNumber = 0;
    mMagField = 0.0;
    /// Initializing CCDB
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(
                    std::chrono::system_clock::now().time_since_epoch())
                    .count();
    ccdb->setCreatedNotAfter(now);
  }

  /// Function to retrieve the nominal mgnetic field in kG (0.1T) and convert it
  /// directly to T
  float getMagneticFieldTesla(uint64_t timestamp)
  {
    // TODO done only once (and not per run). Will be replaced by
    // CCDBConfigurable
    float output = -999;

    if (ConfIsRun3 && !ConfIsMC) {
      static o2::parameters::GRPMagField* grpo = nullptr;
      if (grpo == nullptr) {
        grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(
          "GLO/Config/GRPMagField", timestamp);
        if (grpo == nullptr) {
          LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
          return 0;
        }
        LOGF(info, "Retrieved GRP for timestamp %llu with L3 ", timestamp,
             grpo->getL3Current());
      }
      // taken from GRP onject definition of getNominalL3Field; update later to
      // something smarter (mNominalL3Field = std::lround(5.f * mL3Current /
      // 30000.f);)
      auto NominalL3Field = std::lround(5.f * grpo->getL3Current() / 30000.f);
      output = 0.1 * (NominalL3Field);
    }

    if (!ConfIsRun3 || (ConfIsRun3 && ConfIsMC)) {
      static o2::parameters::GRPObject* grpo = nullptr;
      if (grpo == nullptr) {
        grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>("GLO/GRP/GRP",
                                                                timestamp);
        if (grpo == nullptr) {
          LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
          return 0;
        }
        LOGF(info,
             "Retrieved GRP for timestamp %llu with magnetic field of %d kG",
             timestamp, grpo->getNominalL3Field());
      }
      output = 0.1 * (grpo->getNominalL3Field());
    }
    return output;
  }

  void process(aod::FemtoFullCollision const& col,
               aod::BCsWithTimestamps const&,
               aod::FemtoFullTracks const& tracks,
               o2::aod::V0Datas const& fullV0s) /// \todo with FilteredFullV0s
  {
    // get magnetic field for run
    auto bc = col.bc_as<aod::BCsWithTimestamps>();
    if (mRunNumber != bc.runNumber()) {
      mMagField = getMagneticFieldTesla(bc.timestamp());
      mRunNumber = bc.runNumber();
    }

    const auto vtxZ = col.posZ();
    const auto spher = colCuts.computeSphericity(col, tracks);

    int mult = 0;
    int multNtr = 0;
    if (ConfIsRun3) {
      mult = col.multFV0M();
      multNtr = col.multNTracksPV();
    } else {
      mult = 0.5 * (col.multFV0M()); /// For benchmarking on Run 2, V0M in
                                     /// FemtoUniverseRun2 is defined V0M/2
      multNtr = col.multTracklets();
    }
    if (ConfUseTPCmult) {
      multNtr = col.multTPC();
    }

    /// First thing to do is to check whether the basic event selection criteria
    /// are fulfilled
    // If the basic selection is NOT fulfilled:
    // in case of skimming run - don't store such collisions
    // in case of trigger run - store such collisions but don't store any
    // particle candidates for such collisions
    if (!colCuts.isSelected(col)) {
      if (ConfIsTrigger) {
        outputCollision(vtxZ, mult, multNtr, spher, mMagField);
      }
      return;
    }

    colCuts.fillQA(col);
    // now the table is filled
    outputCollision(vtxZ, mult, multNtr, spher, mMagField);

    std::vector<int> childIDs = {0, 0}; // these IDs are necessary to keep track of the children
    std::vector<int>
      tmpIDtrack; // this vector keeps track of the matching of the primary
                  // track table row <-> aod::track table global index

    if (ConfStoreV0) {
      for (auto& v0 : fullV0s) {
        auto postrack = v0.posTrack_as<aod::FemtoFullTracks>();
        auto negtrack =
          v0.negTrack_as<aod::FemtoFullTracks>(); ///\tocheck funnily enough
                                                  /// if we apply the filter
                                                  /// the sign of Pos and Neg
                                                  /// track is always negative
        // const auto dcaXYpos = postrack.dcaXY();
        // const auto dcaZpos = postrack.dcaZ();
        // const auto dcapos = std::sqrt(pow(dcaXYpos, 2.) + pow(dcaZpos, 2.));
        v0Cuts.fillLambdaQA(col, v0, postrack, negtrack);

        if (!v0Cuts.isSelectedMinimal(col, v0, postrack, negtrack)) {
          continue;
        }
        v0Cuts.fillQA<aod::femtouniverseparticle::ParticleType::kV0,
                      aod::femtouniverseparticle::ParticleType::kV0Child>(
          col, v0, postrack, negtrack); ///\todo fill QA also for daughters
        auto cutContainerV0 =
          v0Cuts.getCutContainer<aod::femtouniverseparticle::CutContainerType>(
            col, v0, postrack, negtrack);

        if ((cutContainerV0.at(
               femto_universe_v0_selection::V0ContainerPosition::kV0) > 0) &&
            (cutContainerV0.at(
               femto_universe_v0_selection::V0ContainerPosition::kPosCuts) > 0) &&
            (cutContainerV0.at(
               femto_universe_v0_selection::V0ContainerPosition::kNegCuts) > 0)) {
          int postrackID = v0.posTrackId();
          int rowInPrimaryTrackTablePos = -1;
          rowInPrimaryTrackTablePos = getRowDaughters(postrackID, tmpIDtrack);
          childIDs[0] = rowInPrimaryTrackTablePos;
          childIDs[1] = 0;
          outputParts(outputCollision.lastIndex(), v0.positivept(),
                      v0.positiveeta(), v0.positivephi(),
                      aod::femtouniverseparticle::ParticleType::kV0Child,
                      cutContainerV0.at(
                        femto_universe_v0_selection::V0ContainerPosition::kPosCuts),
                      cutContainerV0.at(
                        femto_universe_v0_selection::V0ContainerPosition::kPosPID),
                      0., childIDs, 0, 0);
          const int rowOfPosTrack = outputParts.lastIndex();
          int negtrackID = v0.negTrackId();
          int rowInPrimaryTrackTableNeg = -1;
          rowInPrimaryTrackTableNeg = getRowDaughters(negtrackID, tmpIDtrack);
          childIDs[0] = 0;
          childIDs[1] = rowInPrimaryTrackTableNeg;
          outputParts(outputCollision.lastIndex(), v0.negativept(),
                      v0.negativeeta(), v0.negativephi(),
                      aod::femtouniverseparticle::ParticleType::kV0Child,
                      cutContainerV0.at(
                        femto_universe_v0_selection::V0ContainerPosition::kNegCuts),
                      cutContainerV0.at(
                        femto_universe_v0_selection::V0ContainerPosition::kNegPID),
                      0., childIDs, 0, 0);
          const int rowOfNegTrack = outputParts.lastIndex();
          std::vector<int> indexChildID = {rowOfPosTrack, rowOfNegTrack};
          outputParts(outputCollision.lastIndex(), v0.pt(), v0.eta(), v0.phi(),
                      aod::femtouniverseparticle::ParticleType::kV0,
                      cutContainerV0.at(
                        femto_universe_v0_selection::V0ContainerPosition::kV0),
                      0, v0.v0cosPA(),
                      indexChildID, v0.mLambda(), v0.mAntiLambda());
          if (ConfDebugOutput) {
            outputDebugParts(
              postrack.sign(), (uint8_t)postrack.tpcNClsFound(),
              postrack.tpcNClsFindable(),
              (uint8_t)postrack.tpcNClsCrossedRows(),
              postrack.tpcNClsShared(), postrack.tpcFractionSharedCls(), postrack.tpcInnerParam(),
              postrack.itsNCls(), postrack.itsNClsInnerBarrel(),
              postrack.dcaXY(), postrack.dcaZ(), postrack.tpcSignal(),
              postrack.tpcNSigmaStoreEl(), postrack.tpcNSigmaStorePi(),
              postrack.tpcNSigmaStoreKa(), postrack.tpcNSigmaStorePr(),
              postrack.tpcNSigmaStoreDe(), postrack.tofNSigmaStoreEl(),
              postrack.tofNSigmaStorePi(), postrack.tofNSigmaStoreKa(),
              postrack.tofNSigmaStorePr(), postrack.tofNSigmaStoreDe(), -999.,
              -999., -999., -999., -999.,
              -999.); // QA for positive daughter
            outputDebugParts(
              negtrack.sign(), (uint8_t)negtrack.tpcNClsFound(),
              negtrack.tpcNClsFindable(),
              (uint8_t)negtrack.tpcNClsCrossedRows(),
              negtrack.tpcNClsShared(), negtrack.tpcFractionSharedCls(), negtrack.tpcInnerParam(),
              negtrack.itsNCls(), negtrack.itsNClsInnerBarrel(),
              negtrack.dcaXY(), negtrack.dcaZ(), negtrack.tpcSignal(),
              negtrack.tpcNSigmaStoreEl(), negtrack.tpcNSigmaStorePi(),
              negtrack.tpcNSigmaStoreKa(), negtrack.tpcNSigmaStorePr(),
              negtrack.tpcNSigmaStoreDe(), negtrack.tofNSigmaStoreEl(),
              negtrack.tofNSigmaStorePi(), negtrack.tofNSigmaStoreKa(),
              negtrack.tofNSigmaStorePr(), negtrack.tofNSigmaStoreDe(), -999.,
              -999., -999., -999., -999.,
              -999.); // QA for negative daughter
            outputDebugParts(-999., -999., -999., -999., -999., -999., -999., -999.,
                             -999., -999., -999., -999., -999., -999., -999.,
                             -999., -999., -999., -999., -999., -999., -999.,
                             v0.dcaV0daughters(), v0.v0radius(), v0.x(), v0.y(),
                             v0.z(),
                             v0.mK0Short()); // QA for V0
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<femtoUniverseProducerTaskV0Only>(cfgc)};
  return workflow;
}
