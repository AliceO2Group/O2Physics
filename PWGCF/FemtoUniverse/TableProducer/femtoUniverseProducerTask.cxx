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

/// \file femtoUniverseProducerTask.cxx
/// \brief Tasks that produces the track tables used for the pairing
/// \author Laura Serksnyte, TU München, laura.serksnyte@tum.de
/// \author Zuzanna Chochulska, WUT Warsaw & CTU Prague, zchochul@cern.ch
/// \author Malgorzata Janik, WUT Warsaw, majanik@cern.ch
/// \author Pritam Chakraborty, WUT Warsaw, pritam.chakraborty@cern.ch
/// \author Shirajum Monira, WUT Warsaw, shirajum.monira@cern.ch

#include "PWGCF/FemtoUniverse/Core/FemtoUniverseCascadeSelection.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseCollisionSelection.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniversePhiSelection.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseTrackSelection.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseV0Selection.h"
#include "PWGCF/FemtoUniverse/Core/femtoUtils.h"
#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"
#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"
#include <CCDB/BasicCCDBManager.h>

#include "Math/Vector4D.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include <TPDGCode.h>

#include <algorithm>
#include <experimental/type_traits>
#include <set>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::analysis::femto_universe;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

namespace o2::aod
{

using FemtoFullCollision = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>::iterator;
using FemtoFullCollisionCentPP = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::Mults>::iterator;
using FemtoFullCollisionCentRun2 = soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms, aod::Mults>::iterator;
using FemtoFullCollisionCentRun3 = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::Mults>::iterator;
using FemtoFullCollisionMC = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::McCollisionLabels>::iterator;
using FemtoFullCollisionCentRun3MCs = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::Mults, aod::McCollisionLabels>;
using FemtoFullCollisionCentRun3MC = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::Mults, aod::McCollisionLabels>::iterator;
using FemtoFullTracks = soa::Join<aod::FullTracks, aod::TracksDCA, aod::TOFSignal, aod::TrackSelection,
                                  aod::pidTPCEl, aod::pidTPCMu, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTPCDe,
                                  aod::pidTOFEl, aod::pidTOFMu, aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr, aod::pidTOFDe>;

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

struct FemtoUniverseProducerTask {
  Produces<aod::FdCollisions> outputCollision;
  Produces<aod::FDExtCollisions> outputCollExtra;
  Produces<aod::FDParticles> outputParts;
  Produces<aod::FdMCParticles> outputPartsMC;
  Produces<aod::FDExtParticles> outputDebugParts;
  Produces<aod::FDMCLabels> outputPartsMCLabels;
  Produces<aod::FDExtMCParticles> outputDebugPartsMC;
  Produces<aod::FDCascParticles> outputCascParts;

  Configurable<bool> confIsDebug{"confIsDebug", true, "Enable Debug tables"};
  Configurable<bool> confIsUseCutculator{"confIsUseCutculator", true, "Enable cutculator for track cuts"};
  // Choose if filtering or skimming version is run
  // Configurable<bool> confIsTrigger{"confIsTrigger", false, "Store all collisions"}; //Commented: not used configurable
  // Choose if running on converted data or Run3  / Pilot
  Configurable<bool> confIsRun3{"confIsRun3", true, "Running on Run3 or pilot"};
  // Configurable<bool> confIsMC{"confIsMC", false, "Running on MC; implemented only for Run3"}; //Commented: not used configurable

  Configurable<bool> confIsForceGRP{"confIsForceGRP", false, "Set true if the magnetic field configuration is not available in the usual CCDB directory (e.g. for Run 2 converted data or unanchorad Monte Carlo)"};

  Configurable<bool> confDoSpher{"confDoSpher", false, "Calculate sphericity. If false sphericity will take value of 2."};
  Configurable<bool> confStoreMCmothers{"confStoreMCmothers", false, "MC truth: Fill with not only primary particles and store mothers' PDG in tempFitVar."};
  Configurable<bool> confFillCollExt{"confFillCollExt", false, "Option to fill collision extended table"};

  Configurable<bool> confCollMCTruthOnlyReco{"confCollMCTruthOnlyReco", false, "Fill only MC truth collisions that were reconstructed and selected"};
  Configurable<bool> confFillMCTruthV0Daugh{"confFillMCTruthV0Daugh", true, "Fill MC truth daughters of V0"};

  /// Event filtering (used for v0-cascade analysis)
  Configurable<std::string> zorroMask{"zorroMask", "", "zorro trigger class to select on (empty: none)"};

  /// Event cuts
  FemtoUniverseCollisionSelection colCuts;
  struct : o2::framework::ConfigurableGroup {
    Configurable<bool> confEvtUseTPCmult{"confEvtUseTPCmult", false, "Use multiplicity based on the number of tracks with TPC information"};
    Configurable<float> confEvtZvtx{"confEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};
    Configurable<bool> confEvtTriggerCheck{"confEvtTriggerCheck", true, "Evt sel: check for trigger"};
    Configurable<int> confEvtTriggerSel{"confEvtTriggerSel", kINT7, "Evt sel: trigger"};
    Configurable<bool> confEvtOfflineCheck{"confEvtOfflineCheck", false, "Evt sel: check for offline selection"};
    Configurable<bool> confIsActivateV0{"confIsActivateV0", false, "Activate filling of V0 into femtouniverse tables"};
    Configurable<bool> confActivateSecondaries{"confActivateSecondaries", false, "Fill secondary MC gen particles that were reconstructed"};
    Configurable<bool> confIsActivateCascade{"confIsActivateCascade", false, "Activate filling of Cascade into femtouniverse tables"};
    Configurable<bool> confIsActivatePhi{"confIsActivatePhi", false, "Activate filling of Phi into femtouniverse tables"};
    Configurable<bool> confIsActiveD0{"confIsActiveD0", false, "Activate filling FU tables for D0/D0bar mesons"};
    Configurable<bool> confMCTruthAnalysisWithPID{"confMCTruthAnalysisWithPID", true, "1: take only particles with specified PDG, 0: all particles (for MC Truth)"};
    Configurable<std::vector<int>> confMCTruthPDGCodes{"confMCTruthPDGCodes", std::vector<int>{211, -211, 2212, -2212, 333}, "PDG of particles to be stored"};
    Configurable<float> confCentFT0Min{"confCentFT0Min", 0.f, "Min CentFT0 value for centrality selection"};
    Configurable<float> confCentFT0Max{"confCentFT0Max", 200.f, "Max CentFT0 value for centrality selection"};
    Configurable<bool> confEvIsGoodZvtxFT0vsPV{"confEvIsGoodZvtxFT0vsPV", true, "Require kIsGoodZvtxFT0vsPV selection on Events."};
    Configurable<bool> confEvNoSameBunchPileup{"confEvNoSameBunchPileup", true, "Require kNoSameBunchPileup selection on Events."};
    Configurable<bool> confIsUsePileUp{"confIsUsePileUp", true, "Required for choosing whether to run the pile-up cuts"};
    Configurable<bool> confEvIsVertexITSTPC{"confEvIsVertexITSTPC", true, "Require kIsVertexITSTPC selection on Events"};
    Configurable<bool> confIsGoodITSLayersAll{"confIsGoodITSLayersAll", true, "Require IsGoodITSLayersAll selection on Events."};
    Configurable<bool> confNoITSROFrameBorder{"confNoITSROFrameBorder", true, "Require NoITSROFrameBorder selection on Events."};
    Configurable<bool> confNoTimeFrameBorder{"confNoTimeFrameBorder", true, "Require kNoTimeFrameBorder selection on Events."};
    Configurable<bool> confNoCollInRofStandard{"confNoCollInRofStandard", true, "Require NoCollInRofStandard selection on Events."};
    Configurable<bool> confNoHighMultCollInPrevRof{"confNoHighMultCollInPrevRof", true, "Require NoHighMultCollInPrevRof selection on Events."};
    Configurable<bool> confNoCollInTimeRangeStandard{"confNoCollInTimeRangeStandard", true, "Require NoCollInTimeRangeStandard selection on Events."};
    Configurable<int> confTPCOccupancyMin{"confTPCOccupancyMin", 0, "Minimum value for TPC Occupancy selection"};
    Configurable<int> confTPCOccupancyMax{"confTPCOccupancyMax", 500, "Maximum value for TPC Occupancy selection"};
    Configurable<bool> confIsCent{"confIsCent", true, "Centrality or multiplicity selection"};
  } ConfGeneral;
  Filter customCollCentFilter = (aod::cent::centFT0C > ConfGeneral.confCentFT0Min) &&
                                (aod::cent::centFT0C < ConfGeneral.confCentFT0Max);

  // just sanity check to make sure in case there are problems in conversion or
  // MC production it does not affect results
  Configurable<bool> confTrkRejectNotPropagated{"confTrkRejectNotPropagated", false, "True: reject not propagated tracks"};
  // Configurable<bool> ConfRejectITSHitandTOFMissing{
  //     "ConfRejectITSHitandTOFMissing", false,
  //     "True: reject if neither ITS hit nor TOF timing satisfied"};

  FemtoUniverseTrackSelection trackCuts;
  struct : o2::framework::ConfigurableGroup {
    Configurable<std::vector<float>> confTrkCharge{FemtoUniverseTrackSelection::getSelectionName(femto_universe_track_selection::kSign, "ConfTrk"), std::vector<float>{-1, 1}, FemtoUniverseTrackSelection::getSelectionHelper(femto_universe_track_selection::kSign, "Track selection: ")};
    Configurable<std::vector<float>> confTrkPtmin{FemtoUniverseTrackSelection::getSelectionName(femto_universe_track_selection::kpTMin, "ConfTrk"), std::vector<float>{0.5f, 0.4f, 0.6f}, FemtoUniverseTrackSelection::getSelectionHelper(femto_universe_track_selection::kpTMin, "Track selection: ")};
    Configurable<std::vector<float>> confTrkPtmax{FemtoUniverseTrackSelection::getSelectionName(femto_universe_track_selection::kpTMax, "ConfTrk"), std::vector<float>{5.4f, 5.6f, 5.5f}, FemtoUniverseTrackSelection::getSelectionHelper(femto_universe_track_selection::kpTMax, "Track selection: ")};
    Configurable<std::vector<float>> confTrkEta{FemtoUniverseTrackSelection::getSelectionName(femto_universe_track_selection::kEtaMax, "ConfTrk"), std::vector<float>{0.8f, 0.7f, 0.9f}, FemtoUniverseTrackSelection::getSelectionHelper(femto_universe_track_selection::kEtaMax, "Track selection: ")};
    Configurable<std::vector<float>> confTrkTPCnclsMin{FemtoUniverseTrackSelection::getSelectionName(femto_universe_track_selection::kTPCnClsMin, "ConfTrk"), std::vector<float>{70.f}, FemtoUniverseTrackSelection::getSelectionHelper(femto_universe_track_selection::kTPCnClsMin, "Track selection: ")};
    Configurable<std::vector<float>> confTrkTPCfCls{FemtoUniverseTrackSelection::getSelectionName(femto_universe_track_selection::kTPCfClsMin, "ConfTrk"), std::vector<float>{0.83f}, FemtoUniverseTrackSelection::getSelectionHelper(femto_universe_track_selection::kTPCfClsMin, "Track selection: ")};
    Configurable<std::vector<float>> confTrkTPCcRowsMin{FemtoUniverseTrackSelection::getSelectionName(femto_universe_track_selection::kTPCcRowsMin, "ConfTrk"), std::vector<float>{70.f, 60.f, 80.f}, FemtoUniverseTrackSelection::getSelectionHelper(femto_universe_track_selection::kTPCcRowsMin, "Track selection: ")};
    Configurable<std::vector<float>> confTrkTPCsCls{FemtoUniverseTrackSelection::getSelectionName(femto_universe_track_selection::kTPCsClsMax, "ConfTrk"), std::vector<float>{0.1f, 160.f}, FemtoUniverseTrackSelection::getSelectionHelper(femto_universe_track_selection::kTPCsClsMax, "Track selection: ")};
    Configurable<std::vector<float>> confTrkTPCfracsCls{FemtoUniverseTrackSelection::getSelectionName(femto_universe_track_selection::kTPCfracsClsMax, "ConfTrk"), std::vector<float>{0.1f, 160.f}, FemtoUniverseTrackSelection::getSelectionHelper(femto_universe_track_selection::kTPCfracsClsMax, "Track selection: ")};
    Configurable<std::vector<float>> confTrkITSnclsMin{FemtoUniverseTrackSelection::getSelectionName(femto_universe_track_selection::kITSnClsMin, "ConfTrk"), std::vector<float>{-1.f, 2.f, 4.f}, FemtoUniverseTrackSelection::getSelectionHelper(femto_universe_track_selection::kITSnClsMin, "Track selection: ")};
    Configurable<std::vector<float>> confTrkITSnclsIbMin{FemtoUniverseTrackSelection::getSelectionName(femto_universe_track_selection::kITSnClsIbMin, "ConfTrk"), std::vector<float>{-1.f, 1.f}, FemtoUniverseTrackSelection::getSelectionHelper(femto_universe_track_selection::kITSnClsIbMin, "Track selection: ")};
    Configurable<std::vector<float>> confTrkDCAxyMax{FemtoUniverseTrackSelection::getSelectionName(femto_universe_track_selection::kDCAxyMax, "ConfTrk"), std::vector<float>{0.1f, 3.5f}, FemtoUniverseTrackSelection::getSelectionHelper(femto_universe_track_selection::kDCAxyMax, "Track selection: ")};
    Configurable<std::vector<float>> confTrkDCAzMax{FemtoUniverseTrackSelection::getSelectionName(femto_universe_track_selection::kDCAzMax, "ConfTrk"), std::vector<float>{0.2f}, FemtoUniverseTrackSelection::getSelectionHelper(femto_universe_track_selection::kDCAzMax, "Track selection: ")}; /// \todo Reintegrate PID to the general selection container
    Configurable<std::vector<float>> confTrkPIDnSigmaMax{FemtoUniverseTrackSelection::getSelectionName(femto_universe_track_selection::kPIDnSigmaMax, "ConfTrk"), std::vector<float>{3.5f, 3.f, 2.5f}, FemtoUniverseTrackSelection::getSelectionHelper(femto_universe_track_selection::kPIDnSigmaMax, "Track selection: ")};
    Configurable<std::vector<int>> confTrkPIDspecies{"confTrkPIDspecies", std::vector<int>{o2::track::PID::Pion, o2::track::PID::Kaon, o2::track::PID::Proton, o2::track::PID::Deuteron}, "Trk sel: Particles species for PID (Pion=2, Kaon=3, Proton=4, Deuteron=5)"};
    Configurable<bool> confIsOnlyMCTrack{"confIsOnlyMCTrack", false, "Enable filling of only MC Tracks"};
    // Numbers from ~/alice/O2/DataFormats/Reconstruction/include/ReconstructionDataFormats/PID.h //static constexpr ID Pion = 2; static constexpr ID Kaon = 3; static constexpr ID Proton = 4; static constexpr ID Deuteron = 5;
  } ConfTrkSelection;

  Configurable<float> confTrkPIDnSigmaOffsetTPC{"confTrkPIDnSigmaOffsetTPC", 0., "Offset for TPC nSigma because of bad calibration"};
  Configurable<float> confTrkPIDnSigmaOffsetTOF{"confTrkPIDnSigmaOffsetTOF", 0., "Offset for TOF nSigma because of bad calibration"};
  Configurable<float> confTOFpTmin{"confTOFpTmin", 500, "TOF pT min"};

  // TrackSelection *o2PhysicsTrackSelection;
  /// \todo Labeled array (see Track-Track task)

  // V0
  FemtoUniverseV0Selection v0Cuts;
  struct : o2::framework::ConfigurableGroup {
    Configurable<std::vector<float>> confV0Sign{FemtoUniverseV0Selection::getSelectionName(femto_universe_v0_selection::kV0Sign, "ConfV0"), std::vector<float>{-1, 1}, FemtoUniverseV0Selection::getSelectionHelper(femto_universe_v0_selection::kV0Sign, "V0 selection: ")};
    Configurable<std::vector<float>> confV0PtMin{FemtoUniverseV0Selection::getSelectionName(femto_universe_v0_selection::kV0pTMin, "ConfV0"), std::vector<float>{0.3f, 0.4f, 0.5f}, FemtoUniverseV0Selection::getSelectionHelper(femto_universe_v0_selection::kV0pTMin, "V0 selection: ")};
    Configurable<std::vector<float>> confV0PtMax{FemtoUniverseV0Selection::getSelectionName(femto_universe_v0_selection::kV0pTMax, "ConfV0"), std::vector<float>{3.3f, 3.4f, 3.5f}, FemtoUniverseV0Selection::getSelectionHelper(femto_universe_v0_selection::kV0pTMax, "V0 selection: ")};
    Configurable<std::vector<float>> confV0EtaMax{FemtoUniverseV0Selection::getSelectionName(femto_universe_v0_selection::kV0etaMax, "ConfV0"), std::vector<float>{0.8f, 0.7f, 0.9f}, FemtoUniverseV0Selection::getSelectionHelper(femto_universe_v0_selection::kV0etaMax, "V0 selection: ")};
    Configurable<std::vector<float>> confV0DCADaughMax{FemtoUniverseV0Selection::getSelectionName(femto_universe_v0_selection::kV0DCADaughMax, "ConfV0"), std::vector<float>{1.2f, 1.5f}, FemtoUniverseV0Selection::getSelectionHelper(femto_universe_v0_selection::kV0DCADaughMax, "V0 selection: ")};
    Configurable<std::vector<float>> confV0CPAMin{FemtoUniverseV0Selection::getSelectionName(femto_universe_v0_selection::kV0CPAMin, "ConfV0"), std::vector<float>{0.99f, 0.995f}, FemtoUniverseV0Selection::getSelectionHelper(femto_universe_v0_selection::kV0CPAMin, "V0 selection: ")};
    Configurable<std::vector<float>> confV0TranRadMin{FemtoUniverseV0Selection::getSelectionName(femto_universe_v0_selection::kV0TranRadMin, "ConfV0"), std::vector<float>{0.2f}, FemtoUniverseV0Selection::getSelectionHelper(femto_universe_v0_selection::kV0TranRadMin, "V0 selection: ")};
    Configurable<std::vector<float>> confV0TranRadMax{FemtoUniverseV0Selection::getSelectionName(femto_universe_v0_selection::kV0TranRadMax, "ConfV0"), std::vector<float>{100.f}, FemtoUniverseV0Selection::getSelectionHelper(femto_universe_v0_selection::kV0TranRadMax, "V0 selection: ")};
    Configurable<std::vector<float>> confV0DecVtxMax{FemtoUniverseV0Selection::getSelectionName(femto_universe_v0_selection::kV0DecVtxMax, "ConfV0"), std::vector<float>{100.f}, FemtoUniverseV0Selection::getSelectionHelper(femto_universe_v0_selection::kV0DecVtxMax, "V0 selection: ")};

    Configurable<std::vector<float>> confChildCharge{"confChildCharge", std::vector<float>{-1, 1}, "V0 Child sel: Charge"};
    Configurable<std::vector<float>> confChildEtaMax{"confChildEtaMax", std::vector<float>{0.8f}, "V0 Child sel: max eta"};
    Configurable<std::vector<float>> confChildTPCnClsMin{"confChildTPCnClsMin", std::vector<float>{80.f, 70.f, 60.f}, "V0 Child sel: Min. nCls TPC"};
    Configurable<std::vector<float>> confChildDCAMin{"confChildDCAMin", std::vector<float>{0.05f, 0.06f}, "V0 Child sel:  Max. DCA Daugh to PV (cm)"};
    Configurable<std::vector<float>> confChildPIDnSigmaMax{"confChildPIDnSigmaMax", std::vector<float>{5.f, 4.f}, "V0 Child sel: Max. PID nSigma TPC"};
    Configurable<std::vector<int>> confChildPIDspecies{"confChildPIDspecies", std::vector<int>{o2::track::PID::Pion, o2::track::PID::Proton}, "V0 Child sel: Particles species for PID"};

    Configurable<float> confV0InvMassLowLimit{"confV0InvMassLowLimit", 1.05, "Lower limit of the V0 invariant mass"};
    Configurable<float> confV0InvMassUpLimit{"confV0InvMassUpLimit", 1.30, "Upper limit of the V0 invariant mass"};

    Configurable<bool> confV0RejectKaons{"confV0RejectKaons", false, "Switch to reject kaons"};
    Configurable<float> confV0InvKaonMassLowLimit{"confV0InvKaonMassLowLimit", 0.48, "Lower limit of the V0 invariant mass for Kaon rejection"};
    Configurable<float> confV0InvKaonMassUpLimit{"confV0InvKaonMassUpLimit", 0.515, "Upper limit of the V0 invariant mass for Kaon rejection"};

    Configurable<std::vector<int>> confV0PDGMCTruth{"confV0PDGMCTruth", std::vector<int>{2212, -211, 3122}, "PDG codes of V0 daughters and mother, the order must be as follows -- positive daughter, negative daughter, mother"};
  } ConfV0Selection;

  struct : o2::framework::ConfigurableGroup {
    Configurable<float> confPtLowFilterCut{"confPtLowFilterCut", 0.14, "Lower limit for Pt for the global track"};   // pT low
    Configurable<float> confPtHighFilterCut{"confPtHighFilterCut", 5.0, "Higher limit for Pt for the global track"}; // pT high
    Configurable<float> confEtaFilterCut{"confEtaFilterCut", 0.8, "Eta cut for the global track"};                   // eta
    Configurable<float> confDcaXYCustom1FilterCut{"confDcaXYCustom1FilterCut", 0.0105, "Value for [1] custom DCAxy cut -> |DCAxy| < [1] + [2]/pT"};
    Configurable<float> confDcaXYCustom2FilterCut{"confDcaXYCustom2FilterCut", 0.035, "Value for [2] custom DCAxy cut -> |DCAxy| < [1] + [2]/pT"};
    Configurable<float> confDcaZCustom1FilterCut{"confDcaZCustom1FilterCut", 0.02, "Value for [1] custom cut on DCAz -> |DCAz| < [1] + [2]/pT"};
    Configurable<float> confDcaZCustom2FilterCut{"confDcaZCustom2FilterCut", 0.0, "Value for [2] custom cut on DCAz -> |DCAz| < [1] + [2]/pT"};
    Configurable<bool> confIsApplyTrkCutMCTruth{"confIsApplyTrkCutMCTruth", false, "Apply eta, pT selection cut on MCTruth tracks "};
    Configurable<bool> confIsOnlyPrimary{"confIsOnlyPrimary", false, "Select only primaries"};
  } ConfFilterCuts;

  Filter globalCutFilter = requireGlobalTrackWoDCAInFilter();
  Filter customTrackFilter = (aod::track::pt > ConfFilterCuts.confPtLowFilterCut) &&
                             (aod::track::pt < ConfFilterCuts.confPtHighFilterCut) &&
                             (nabs(aod::track::eta) < ConfFilterCuts.confEtaFilterCut) &&
                             (nabs(aod::track::dcaZ) < (ConfFilterCuts.confDcaZCustom1FilterCut + ConfFilterCuts.confDcaZCustom2FilterCut / aod::track::pt)) &&
                             (nabs(aod::track::dcaXY) < (ConfFilterCuts.confDcaXYCustom1FilterCut + ConfFilterCuts.confDcaXYCustom2FilterCut / aod::track::pt));

  // CASCADE
  FemtoUniverseCascadeSelection cascadeCuts;
  struct : o2::framework::ConfigurableGroup {
    Configurable<std::vector<float>> confCascSign{FemtoUniverseCascadeSelection::getSelectionName(femto_universe_cascade_selection::kCascadeSign, "ConfCasc"), std::vector<float>{-1, 1}, FemtoUniverseCascadeSelection::getSelectionHelper(femto_universe_cascade_selection::kCascadeSign, "Cascade selection: ")};
    Configurable<std::vector<float>> confCascPtMin{FemtoUniverseCascadeSelection::getSelectionName(femto_universe_cascade_selection::kCascadepTMin, "ConfCasc"), std::vector<float>{0.3f, 0.4f, 0.5f}, FemtoUniverseCascadeSelection::getSelectionHelper(femto_universe_cascade_selection::kCascadepTMin, "Cascade selection: ")};
    Configurable<std::vector<float>> confCascPtMax{FemtoUniverseCascadeSelection::getSelectionName(femto_universe_cascade_selection::kCascadepTMax, "ConfCasc"), std::vector<float>{3.3f, 3.4f, 3.5f}, FemtoUniverseCascadeSelection::getSelectionHelper(femto_universe_cascade_selection::kCascadepTMax, "Cascade selection: ")};
    Configurable<std::vector<float>> confCascEtaMax{FemtoUniverseCascadeSelection::getSelectionName(femto_universe_cascade_selection::kCascadeetaMax, "ConfCasc"), std::vector<float>{0.8f, 0.7f, 0.9f}, FemtoUniverseCascadeSelection::getSelectionHelper(femto_universe_cascade_selection::kCascadeetaMax, "Cascade selection: ")};
    Configurable<std::vector<float>> confCascV0DCADaughMax{FemtoUniverseCascadeSelection::getSelectionName(femto_universe_cascade_selection::kCascadeV0DCADaughMax, "ConfCasc"), std::vector<float>{1.f, 1.2f, 1.5f}, FemtoUniverseCascadeSelection::getSelectionHelper(femto_universe_cascade_selection::kCascadeV0DCADaughMax, "Cascade selection: ")};
    Configurable<std::vector<float>> confCascV0CPAMin{FemtoUniverseCascadeSelection::getSelectionName(femto_universe_cascade_selection::kCascadeV0CPAMin, "ConfCasc"), std::vector<float>{0.99f, 0.95f}, FemtoUniverseCascadeSelection::getSelectionHelper(femto_universe_cascade_selection::kCascadeV0CPAMin, "Cascade selection: ")};
    Configurable<std::vector<float>> confCascV0TranRadMin{FemtoUniverseCascadeSelection::getSelectionName(femto_universe_cascade_selection::kCascadeV0TranRadMin, "ConfCasc"), std::vector<float>{0.2f}, FemtoUniverseCascadeSelection::getSelectionHelper(femto_universe_cascade_selection::kCascadeV0TranRadMin, "Cascade selection: ")};
    Configurable<std::vector<float>> confCascV0TranRadMax{FemtoUniverseCascadeSelection::getSelectionName(femto_universe_cascade_selection::kCascadeV0TranRadMax, "ConfCasc"), std::vector<float>{100.f}, FemtoUniverseCascadeSelection::getSelectionHelper(femto_universe_cascade_selection::kCascadeV0TranRadMax, "Cascade selection: ")};
    Configurable<std::vector<float>> confCascV0DecVtxMax{FemtoUniverseCascadeSelection::getSelectionName(femto_universe_cascade_selection::kCascadeV0DecVtxMax, "ConfCasc"), std::vector<float>{100.f}, FemtoUniverseCascadeSelection::getSelectionHelper(femto_universe_cascade_selection::kCascadeV0DecVtxMax, "Cascade selection: ")};
    Configurable<std::vector<float>> confCascDCADaughMax{FemtoUniverseCascadeSelection::getSelectionName(femto_universe_cascade_selection::kCascadeDCADaughMax, "ConfCasc"), std::vector<float>{1.f, 1.2f, 1.5f}, FemtoUniverseCascadeSelection::getSelectionHelper(femto_universe_cascade_selection::kCascadeDCADaughMax, "Cascade selection: ")};
    Configurable<std::vector<float>> confCascCPAMin{FemtoUniverseCascadeSelection::getSelectionName(femto_universe_cascade_selection::kCascadeCPAMin, "ConfCasc"), std::vector<float>{0.99f, 0.95f}, FemtoUniverseCascadeSelection::getSelectionHelper(femto_universe_cascade_selection::kCascadeCPAMin, "Cascade selection: ")};
    Configurable<std::vector<float>> confCascTranRadMin{FemtoUniverseCascadeSelection::getSelectionName(femto_universe_cascade_selection::kCascadeTranRadMin, "ConfCasc"), std::vector<float>{0.2f, 0.5f}, FemtoUniverseCascadeSelection::getSelectionHelper(femto_universe_cascade_selection::kCascadeTranRadMin, "Cascade selection: ")};
    Configurable<std::vector<float>> confCascTranRadMax{FemtoUniverseCascadeSelection::getSelectionName(femto_universe_cascade_selection::kCascadeTranRadMax, "ConfCasc"), std::vector<float>{100.f}, FemtoUniverseCascadeSelection::getSelectionHelper(femto_universe_cascade_selection::kCascadeTranRadMax, "Cascade selection: ")};
    Configurable<std::vector<float>> confCascDecVtxMax{FemtoUniverseCascadeSelection::getSelectionName(femto_universe_cascade_selection::kCascadeDecVtxMax, "ConfCasc"), std::vector<float>{100.f}, FemtoUniverseCascadeSelection::getSelectionHelper(femto_universe_cascade_selection::kCascadeDecVtxMax, "Cascade selection: ")};

    Configurable<std::vector<float>> confCascDCAPosToPV{FemtoUniverseCascadeSelection::getSelectionName(femto_universe_cascade_selection::kCascadeDCAPosToPV, "ConfCasc"), std::vector<float>{0.1f}, FemtoUniverseCascadeSelection::getSelectionHelper(femto_universe_cascade_selection::kCascadeDCAPosToPV, "Cascade selection: ")};
    Configurable<std::vector<float>> confCascDCANegToPV{FemtoUniverseCascadeSelection::getSelectionName(femto_universe_cascade_selection::kCascadeDCANegToPV, "ConfCasc"), std::vector<float>{0.1f}, FemtoUniverseCascadeSelection::getSelectionHelper(femto_universe_cascade_selection::kCascadeDCANegToPV, "Cascade selection: ")};
    Configurable<std::vector<float>> confCascDCABachToPV{FemtoUniverseCascadeSelection::getSelectionName(femto_universe_cascade_selection::kCascadeDCABachToPV, "ConfCasc"), std::vector<float>{0.1f}, FemtoUniverseCascadeSelection::getSelectionHelper(femto_universe_cascade_selection::kCascadeDCABachToPV, "Cascade selection: ")};
    Configurable<std::vector<float>> confCascDCAV0ToPV{FemtoUniverseCascadeSelection::getSelectionName(femto_universe_cascade_selection::kCascadeDCAV0ToPV, "ConfCasc"), std::vector<float>{0.01f}, FemtoUniverseCascadeSelection::getSelectionHelper(femto_universe_cascade_selection::kCascadeDCAV0ToPV, "Cascade selection: ")};
    Configurable<std::vector<float>> confCascV0MassLowLimit{FemtoUniverseCascadeSelection::getSelectionName(femto_universe_cascade_selection::kCascadeV0MassMin, "ConfCasc"), std::vector<float>{1.05f}, FemtoUniverseCascadeSelection::getSelectionHelper(femto_universe_cascade_selection::kCascadeV0MassMin, "Cascade selection: ")};
    Configurable<std::vector<float>> confCascV0MassUpLimit{FemtoUniverseCascadeSelection::getSelectionName(femto_universe_cascade_selection::kCascadeV0MassMax, "ConfCasc"), std::vector<float>{1.30f}, FemtoUniverseCascadeSelection::getSelectionHelper(femto_universe_cascade_selection::kCascadeV0MassMax, "Cascade selection: ")};

    Configurable<std::vector<float>> confCascChildCharge{"confCascChildCharge", std::vector<float>{-1, 1}, "Cascade Child sel: Charge"};
    Configurable<std::vector<float>> confCascChildEtaMax{"confCascChildEtaMax", std::vector<float>{0.8f}, "Cascade Child sel: max eta"};
    Configurable<std::vector<float>> confCascChildTPCnClsMin{"confCascChildTPCnClsMin", std::vector<float>{80.f, 70.f, 60.f}, "Cascade Child sel: Min. nCls TPC"};
    // Configurable<std::vector<float>> confCascChildDCAMin{"confCascChildDCAMin", std::vector<float>{0.05f, 0.06f}, "Cascade Child sel:  Max. DCA Daugh to PV (cm)"}; //Commented: not used variable
    Configurable<std::vector<float>> confCascChildPIDnSigmaMax{"confCascChildPIDnSigmaMax", std::vector<float>{3.f, 4.f}, "Cascade Child sel: Max. PID nSigma TPC"};
    Configurable<std::vector<int>> confCascChildPIDspecies{"confCascChildPIDspecies", std::vector<int>{o2::track::PID::Pion, o2::track::PID::Proton}, "Cascade Child sel: particle species for PID"};

    Configurable<float> confXiInvMassLowLimit{"confXiInvMassLowLimit", 1.25, "Lower limit of the Xi invariant mass"};
    Configurable<float> confXiInvMassUpLimit{"confXiInvMassUpLimit", 1.40, "Upper limit of the Xi invariant mass"};
    Configurable<float> confOmegaInvMassLowLimit{"confOmegaInvMassLowLimit", 1.60, "Lower limit of the Omega invariant mass"};
    Configurable<float> confOmegaInvMassUpLimit{"confOmegaInvMassUpLimit", 1.80, "Upper limit of the Omega invariant mass"};
  } ConfCascadeSelection;

  // PHI
  FemtoUniversePhiSelection phiCuts;
  struct : o2::framework::ConfigurableGroup {
    /// Phi meson
    Configurable<float> confPhiPtLowLimit{"confPhiPtLowLimit", 0.8, "Lower limit of the Phi pT."};
    Configurable<float> confPhiPtHighLimit{"confPhiPtHighLimit", 4.0, "Higher limit of the Phi pT."};
    Configurable<float> confPhiEtaHighLimit{"confPhiEtaHighLimit", 0.8, "Maximum eta value of the Phi"};
    Configurable<float> confPhiInvMassLowLimit{"confPhiInvMassLowLimit", 1.011, "Lower limit of the Phi invariant mass"};
    Configurable<float> confPhiInvMassUpLimit{"confPhiInvMassUpLimit", 1.027, "Upper limit of the Phi invariant mass"};
    // Phi meson daughters
    Configurable<float> confPhiKaonRejectPionNsigma{"confPhiKaonRejectPionNsigma", 3.0, "Reject if particle could be a Pion combined nsigma value."};
    Configurable<float> confPhiKaonRejectProtonNsigma{"confPhiKaonRejectProtonNsigma", 3.0, "Reject if particle could be a Proton combined nsigma value."};
    // Kaons
    Configurable<bool> confPhiDoLFPID4Kaons{"confPhiDoLFPID4Kaons", true, "Switch on do PID for Kaons as in LF"};
    Configurable<float> confNSigmaTPCKaonLF{"confNSigmaTPCKaonLF", 3.0, "TPC Kaon Sigma as in LF"};
    Configurable<float> confNSigmaCombKaonLF{"confNSigmaCombKaonLF", 3.0, "TPC and TOF Kaon Sigma (combined) as in LF"};
    Configurable<float> confMomKaonLF{"confMomKaonLF", 0.5, "Momentum threshold for kaon identification as in LF"};
    Configurable<float> confMomKaonRejected{"confMomKaonRejected", 0.5, "Momentum threshold for rejected kaon"};
    Configurable<float> confMomKaon03{"confMomKaon03", 0.3, "Momentum threshold for kaon identification pT = 0.3 GeV/c"};
    Configurable<float> confMomKaon045{"confMomKaon045", 0.45, "Momentum threshold for kaon identification pT = 0.45 GeV/c"};
    Configurable<float> confMomKaon055{"confMomKaon055", 0.55, "Momentum threshold for kaon identification pT = 0.55 GeV/c"};
    Configurable<float> confMomKaon15{"confMomKaon15", 1.5, "Momentum threshold for kaon identification pT = 1.5 GeV/c"};
    Configurable<float> confPhiKaonNsigmaTPCfrom00to03{"confPhiKaonNsigmaTPCfrom00to03", 3.0, "Reject if Kaons in 0.0-0.3 are have TPC n sigma above this value."};
    Configurable<float> confPhiKaonNsigmaTPCfrom03to045{"confPhiKaonNsigmaTPCfrom03to045", 2.0, "Reject if Kaons in 0.3-0.45 are have TPC n sigma above this value."};
    Configurable<float> confPhiKaonNsigmaTPCfrom045to055{"confPhiKaonNsigmaTPCfrom045to055", 1.0, "Reject if Kaons in 0.45-0.55 are have TPC n sigma above this value."};
    Configurable<float> confPhiKaonNsigmaTPCfrom055to15{"confPhiKaonNsigmaTPCfrom055to15", 3.0, "Reject if Kaons in 0.55-1.5 are have TPC n sigma above this value."};
    Configurable<float> confPhiKaonNsigmaTOFfrom055to15{"confPhiKaonNsigmaTOFfrom055to15", 3.0, "Reject if Kaons in 0.55-1.5 are have TOF n sigma above this value."};
    Configurable<float> confPhiKaonNsigmaTPCfrom15{"confPhiKaonNsigmaTPCfrom15", 3.0, "Reject if Kaons above 1.5 are have TPC n sigma above this value."};
    Configurable<float> confPhiKaonNsigmaTOFfrom15{"confPhiKaonNsigmaTOFfrom15", 3.0, "Reject if Kaons above 1.5 are have TOF n sigma above this value."};
  } ConfPhiSelection;

  // PDG codes for fillMCParticle function
  Configurable<int> confPDGCodePartOne{"confPDGCodePartOne", 321, "Particle 1 - PDG code"};
  Configurable<int> confPDGCodePartTwo{"confPDGCodePartTwo", 321, "Particle 2 - PDG code"};

  // D0/D0bar mesons
  struct : o2::framework::ConfigurableGroup {
    Configurable<float> trackD0CandEtaMax{"trackD0CandEtaMax", 0.8, "max. track/D0 cand. pseudorapidity"};
    Configurable<double> yD0CandGenMax{"yD0CandGenMax", 0.5, "max. gen. D0 cand. rapidity"};
    Configurable<double> yD0CandMax{"yD0CandMax", 0.8, "max. D0 cand. rapidity"};
    Configurable<float> trackD0pTGenMin{"trackD0pTGenMin", 0.0, "MC Truth, min. pT for tracks and D0/D0bar cand."};
    Configurable<float> trackD0pTGenMax{"trackD0pTGenMax", 24.0, "MC Truth, max. pT for tracks and D0/D0bar cand."};
    Configurable<bool> useYCutD0Cand{"useYCutD0Cand", true, "True - apply cut on y of D0 cand./false - apply cut on eta"};
    Configurable<bool> storeD0D0barDoubleMassHypo{"storeD0D0barDoubleMassHypo", false, "Store D0/D0bar cand. which pass selection criteria for both, D0 and D0bar"};
    Configurable<std::vector<int>> classMlD0D0bar{"classMlD0D0bar", {0, 1, 2}, "Indexes of ML scores to be stored. Three indexes max."};
  } ConfD0Selection;

  // PID bitmask configurables
  struct : o2::framework::ConfigurableGroup {
    Configurable<float> confMinMomTOF{"confMinMomTOF", 0.75, "momentum threshold for particle identification using TOF"};
    Configurable<float> confNsigmaTPCParticleChild{"confNsigmaTPCParticleChild", 3.0, "TPC Sigma for particle (daugh & bach) momentum < Confmom"};
    Configurable<float> confNsigmaTOFParticleChild{"confNsigmaTOFParticleChild", 3.0, "TOF Sigma for particle (daugh & bach) momentum > Confmom"};
    Configurable<float> confNsigmaTPCParticle{"confNsigmaTPCParticle", 3.0, "TPC Sigma for particle (track) momentum < Confmom"};
    Configurable<float> confNsigmaCombinedParticle{"confNsigmaCombinedParticle", 3.0, "TPC and TOF Sigma (combined) for particle (track) momentum > Confmom"};
  } ConfPIDBitmask;

  HfHelper hfHelper;
  bool isKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK)
  {

    if (mom < ConfPhiSelection.confMomKaon03) { // 0.0-0.3
      if (std::abs(nsigmaTPCK) < ConfPhiSelection.confPhiKaonNsigmaTPCfrom00to03) {
        return true;
      } else {
        return false;
      }
    } else if (mom < ConfPhiSelection.confMomKaon045) { // 0.30 - 0.45
      if (std::abs(nsigmaTPCK) < ConfPhiSelection.confPhiKaonNsigmaTPCfrom03to045) {
        return true;
      } else {
        return false;
      }
    } else if (mom < ConfPhiSelection.confMomKaon055) { // 0.45-0.55
      if (std::abs(nsigmaTPCK) < ConfPhiSelection.confPhiKaonNsigmaTPCfrom045to055) {
        return true;
      } else {
        return false;
      }
    } else if (mom < ConfPhiSelection.confMomKaon15) { // 0.55-1.5 (now we use TPC and TOF)
      if ((std::abs(nsigmaTOFK) < ConfPhiSelection.confPhiKaonNsigmaTOFfrom055to15) && (std::abs(nsigmaTPCK) < ConfPhiSelection.confPhiKaonNsigmaTPCfrom055to15)) {
        {
          return true;
        }
      } else {
        return false;
      }
    } else if (mom > ConfPhiSelection.confMomKaon15) { // 1.5 -
      if ((std::abs(nsigmaTOFK) < ConfPhiSelection.confPhiKaonNsigmaTOFfrom15) && (std::abs(nsigmaTPCK) < ConfPhiSelection.confPhiKaonNsigmaTPCfrom15)) {
        return true;
      } else {
        return false;
      }
    } else {
      return false;
    }
  }

  bool isKaonNSigmaLF(float mom, float nsigmaTPCK, float nsigmaTOFK, bool hasTOF)
  {
    if (mom < ConfPhiSelection.confMomKaonLF) {
      if (std::abs(nsigmaTPCK) < ConfPhiSelection.confNSigmaTPCKaonLF) {
        return true;
      } else {
        return false;
      }
    } else if (mom >= ConfPhiSelection.confMomKaonLF) { // 0.5-1.5 (now we use TPC and TOF)
      if (!hasTOF) {
        return false;
      } else {
        if (std::sqrt(nsigmaTPCK * nsigmaTPCK + nsigmaTOFK * nsigmaTOFK) < ConfPhiSelection.confNSigmaCombKaonLF) {
          return true;
        } else {
          return false;
        }
      }
    } else {
      return false;
    }
  }

  bool isKaonRejected(float mom, float nsigmaTPCPr, float nsigmaTOFPr, float nsigmaTPCPi, float nsigmaTOFPi)
  {
    if (mom < ConfPhiSelection.confMomKaonRejected) {
      if (std::abs(nsigmaTPCPi) < ConfPhiSelection.confPhiKaonRejectPionNsigma.value) {
        return true;
      } else if (std::abs(nsigmaTPCPr) < ConfPhiSelection.confPhiKaonRejectProtonNsigma.value) {
        return true;
      }
    }
    if (mom > ConfPhiSelection.confMomKaonRejected) {
      if (std::hypot(nsigmaTOFPi, nsigmaTPCPi) < ConfPhiSelection.confPhiKaonRejectPionNsigma.value) {
        return true;
      } else if (std::hypot(nsigmaTOFPr, nsigmaTPCPr) < ConfPhiSelection.confPhiKaonRejectProtonNsigma.value) {
        return true;
      } else {
        return false;
      }
    } else {
      return false;
    }
  }

  bool isNSigmaTPC(float nsigmaTPCParticle)
  {
    return (std::abs(nsigmaTPCParticle) < ConfPIDBitmask.confNsigmaTPCParticleChild);
  }

  bool isNSigmaTOF(float mom, float nsigmaTOFParticle, bool hasTOF)
  {
    // Cut only on daughter and bachelor tracks, that have TOF signal
    if (mom > ConfPIDBitmask.confMinMomTOF && hasTOF) {
      return (std::abs(nsigmaTOFParticle) < ConfPIDBitmask.confNsigmaTOFParticleChild);
    } else {
      return true;
    }
  }

  bool isNSigmaCombined(float mom, float nsigmaTPCParticle, float nsigmaTOFParticle)
  {
    if (mom <= ConfPIDBitmask.confMinMomTOF) {
      return (std::abs(nsigmaTPCParticle) < ConfPIDBitmask.confNsigmaTPCParticle);
    } else {
      return (std::hypot(nsigmaTOFParticle, nsigmaTPCParticle) < ConfPIDBitmask.confNsigmaCombinedParticle);
    }
  }

  template <typename TrackType>
  aod::femtouniverseparticle::CutContainerType PIDBitmask(const TrackType& track)
  {
    static const o2::track::PID pids[] = {o2::track::PID::Proton, o2::track::PID::Pion, o2::track::PID::Kaon};
    aod::femtouniverseparticle::CutContainerType mask = 0u;
    for (UInt_t i = 0; i < 3; ++i) {
      if (isNSigmaTPC(trackCuts.getNsigmaTPC(track, pids[i])))
        mask |= (1u << i);
      if (isNSigmaTOF(track.p(), trackCuts.getNsigmaTOF(track, pids[i]), track.hasTOF()))
        mask |= (8u << i);
      if (isNSigmaCombined(track.p(), trackCuts.getNsigmaTPC(track, pids[i]), trackCuts.getNsigmaTOF(track, pids[i])))
        mask |= (64u << i);
    }
    if (track.hasTOF())
      mask |= (512u);
    return mask;
  }

  template <class T>
  using hasStrangeTOFinV0 = decltype(std::declval<T&>().tofNSigmaLaPr());

  /// bitmask to save strangeness TOF for V0 analysis
  template <typename V0Type>
  aod::femtouniverseparticle::CutContainerType PIDStrangeTOFBitmaskV0(const V0Type& v0)
  {
    aod::femtouniverseparticle::CutContainerType mask = 0u;
    if constexpr (std::experimental::is_detected<hasStrangeTOFinV0, V0Type>::value) {
      if (v0.tofNSigmaLaPr() < ConfPIDBitmask.confNsigmaTOFParticleChild)
        mask |= (1u);
      if (v0.tofNSigmaLaPi() < ConfPIDBitmask.confNsigmaTOFParticleChild)
        mask |= (2u);
      if (v0.tofNSigmaALaPr() < ConfPIDBitmask.confNsigmaTOFParticleChild)
        mask |= (4u);
      if (v0.tofNSigmaALaPi() < ConfPIDBitmask.confNsigmaTOFParticleChild)
        mask |= (8u);
      if (v0.tofNSigmaK0PiPlus() < ConfPIDBitmask.confNsigmaTOFParticleChild)
        mask |= (16u);
      if (v0.tofNSigmaK0PiMinus() < ConfPIDBitmask.confNsigmaTOFParticleChild)
        mask |= (32u);
    }
    return mask;
  }

  template <class T>
  using hasStrangeTOFinCasc = decltype(std::declval<T&>().tofNSigmaXiLaPi());

  /// bitmask to save strangeness TOF for cascade analysis
  template <typename CascType>
  aod::femtouniverseparticle::CutContainerType PIDStrangeTOFBitmaskCasc(const CascType& casc)
  {
    aod::femtouniverseparticle::CutContainerType mask = 0u;
    if constexpr (std::experimental::is_detected<hasStrangeTOFinCasc, CascType>::value) {
      if (casc.tofNSigmaXiLaPi() < ConfPIDBitmask.confNsigmaTOFParticleChild)
        mask |= (1u);
      if (casc.tofNSigmaXiLaPr() < ConfPIDBitmask.confNsigmaTOFParticleChild)
        mask |= (2u);
      if (casc.tofNSigmaXiPi() < ConfPIDBitmask.confNsigmaTOFParticleChild)
        mask |= (4u);
      if (casc.tofNSigmaOmLaPi() < ConfPIDBitmask.confNsigmaTOFParticleChild)
        mask |= (8u);
      if (casc.tofNSigmaOmLaPr() < ConfPIDBitmask.confNsigmaTOFParticleChild)
        mask |= (16u);
      if (casc.tofNSigmaOmKa() < ConfPIDBitmask.confNsigmaTOFParticleChild)
        mask |= (32u);
    }
    return mask;
  }

  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};
  int mRunNumberZorro = 0;

  HistogramRegistry qaRegistry{"QAHistos", {}, OutputObjHandlingPolicy::QAObject};
  HistogramRegistry cascadeQaRegistry{"CascadeQAHistos", {}, OutputObjHandlingPolicy::QAObject};

  void initZorro(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumberZorro == bc.runNumber())
      return;

    zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), zorroMask.value);
    zorro.populateHistRegistry(qaRegistry, bc.runNumber());
    mRunNumberZorro = bc.runNumber();
  }

  /// \todo should we add filter on min value pT/eta of V0 and daughters?
  /*Filter v0Filter = (nabs(aod::v0data::x) < V0DecVtxMax.value) &&
                    (nabs(aod::v0data::y) < V0DecVtxMax.value) &&
                    (nabs(aod::v0data::z) < V0DecVtxMax.value);*/
  // (aod::v0data::v0radius > V0TranRadV0Min.value); to be added, not working
  // for now do not know why

  int mRunNumber = 0;
  float mMagField;
  Service<o2::ccdb::BasicCCDBManager> ccdb; /// Accessing the CCDB
  ctpRateFetcher mRateFetcher;              // inspired by zdcSP.cxx in PWGLF

  void init(InitContext&)
  {
    if ((doprocessFullData || doprocessTrackPhiData || doprocessTrackData || doprocessTrackV0 || doprocessTrackCascadeData || doprocessTrackV0Cascade || doprocessTrackD0mesonData || doprocessTrackD0DataML || doprocessTrackCentRun2Data || doprocessTrackV0CentRun2Data || doprocessTrackCentRun3Data || doprocessV0CentRun3Data || doprocessCascadeCentRun3Data || doprocessTrackDataCentPP) == false && (doprocessFullMC || doprocessTrackMC || doprocessTrackMCTruth || doprocessTrackMCGen || doprocessTruthAndFullMCV0 || doprocessTrackD0MC || doprocessTruthAndFullMCCasc || doprocessFullMCCent || doprocessTrackCentRun3DataMC || doprocessTruthAndFullMCCentRun3 || doprocessTruthAndFullMCCentRun3V0) == false) {
      LOGF(fatal, "Neither processFullData nor processFullMC enabled. Please choose one.");
    }
    if ((doprocessFullData || doprocessTrackPhiData || doprocessTrackData || doprocessTrackV0 || doprocessTrackCascadeData || doprocessTrackV0Cascade || doprocessTrackD0mesonData || doprocessTrackD0DataML || doprocessTrackCentRun2Data || doprocessTrackV0CentRun2Data || doprocessTrackCentRun3Data || doprocessV0CentRun3Data || doprocessCascadeCentRun3Data || doprocessTrackDataCentPP) == true && (doprocessFullMC || doprocessTrackMC || doprocessTrackMCTruth || doprocessTrackMCGen || doprocessTruthAndFullMCV0 || doprocessTrackD0MC || doprocessTruthAndFullMCCasc || doprocessFullMCCent || doprocessTrackCentRun3DataMC || doprocessTruthAndFullMCCentRun3 || doprocessTruthAndFullMCCentRun3V0) == true) {
      LOGF(fatal,
           "Cannot enable process Data and process MC at the same time. "
           "Please choose one.");
    }

    zorroSummary.setObject(zorro.getZorroSummary());

    colCuts.setCuts(ConfGeneral.confEvtZvtx, ConfGeneral.confEvtTriggerCheck, ConfGeneral.confEvtTriggerSel, ConfGeneral.confEvtOfflineCheck, confIsRun3, ConfGeneral.confCentFT0Min, ConfGeneral.confCentFT0Max);
    colCuts.init(&qaRegistry);

    trackCuts.setSelection(ConfTrkSelection.confTrkCharge, femto_universe_track_selection::kSign, femto_universe_selection::kEqual);
    trackCuts.setSelection(ConfTrkSelection.confTrkPtmin, femto_universe_track_selection::kpTMin, femto_universe_selection::kLowerLimit);
    trackCuts.setSelection(ConfTrkSelection.confTrkPtmax, femto_universe_track_selection::kpTMax, femto_universe_selection::kUpperLimit);
    trackCuts.setSelection(ConfTrkSelection.confTrkEta, femto_universe_track_selection::kEtaMax, femto_universe_selection::kAbsUpperLimit);
    trackCuts.setSelection(ConfTrkSelection.confTrkTPCnclsMin, femto_universe_track_selection::kTPCnClsMin, femto_universe_selection::kLowerLimit);
    trackCuts.setSelection(ConfTrkSelection.confTrkTPCfCls, femto_universe_track_selection::kTPCfClsMin, femto_universe_selection::kLowerLimit);
    trackCuts.setSelection(ConfTrkSelection.confTrkTPCcRowsMin, femto_universe_track_selection::kTPCcRowsMin, femto_universe_selection::kLowerLimit);
    trackCuts.setSelection(ConfTrkSelection.confTrkTPCsCls, femto_universe_track_selection::kTPCsClsMax, femto_universe_selection::kUpperLimit);
    trackCuts.setSelection(ConfTrkSelection.confTrkTPCfracsCls, femto_universe_track_selection::kTPCfracsClsMax, femto_universe_selection::kUpperLimit);
    trackCuts.setSelection(ConfTrkSelection.confTrkITSnclsMin, femto_universe_track_selection::kITSnClsMin, femto_universe_selection::kLowerLimit);
    trackCuts.setSelection(ConfTrkSelection.confTrkITSnclsIbMin, femto_universe_track_selection::kITSnClsIbMin, femto_universe_selection::kLowerLimit);
    trackCuts.setSelection(ConfTrkSelection.confTrkDCAxyMax, femto_universe_track_selection::kDCAxyMax, femto_universe_selection::kAbsUpperLimit);
    trackCuts.setSelection(ConfTrkSelection.confTrkDCAzMax, femto_universe_track_selection::kDCAzMax, femto_universe_selection::kAbsUpperLimit);
    trackCuts.setSelection(ConfTrkSelection.confTrkPIDnSigmaMax, femto_universe_track_selection::kPIDnSigmaMax, femto_universe_selection::kAbsUpperLimit);
    trackCuts.setPIDSpecies(ConfTrkSelection.confTrkPIDspecies);
    trackCuts.setnSigmaPIDOffset(confTrkPIDnSigmaOffsetTPC, confTrkPIDnSigmaOffsetTOF);
    trackCuts.init<aod::femtouniverseparticle::ParticleType::kTrack, aod::femtouniverseparticle::TrackType::kNoChild, aod::femtouniverseparticle::CutContainerType>(&qaRegistry);

    /// \todo fix how to pass array to setSelection, getRow() passing a
    /// different type!
    // v0Cuts.setSelection(ConfV0Selection->getRow(0),
    // femto_universe_v0_selection::kDecVtxMax, femto_universe_selection::kAbsUpperLimit);
    if (ConfGeneral.confIsActivateV0) {
      // initializing for V0
      v0Cuts.setSelection(ConfV0Selection.confV0Sign, femto_universe_v0_selection::kV0Sign, femto_universe_selection::kEqual);
      v0Cuts.setSelection(ConfV0Selection.confV0PtMin, femto_universe_v0_selection::kV0pTMin, femto_universe_selection::kLowerLimit);
      v0Cuts.setSelection(ConfV0Selection.confV0PtMax, femto_universe_v0_selection::kV0pTMax, femto_universe_selection::kUpperLimit);
      v0Cuts.setSelection(ConfV0Selection.confV0EtaMax, femto_universe_v0_selection::kV0etaMax, femto_universe_selection::kAbsUpperLimit);
      v0Cuts.setSelection(ConfV0Selection.confV0DCADaughMax, femto_universe_v0_selection::kV0DCADaughMax, femto_universe_selection::kUpperLimit);
      v0Cuts.setSelection(ConfV0Selection.confV0CPAMin, femto_universe_v0_selection::kV0CPAMin, femto_universe_selection::kLowerLimit);
      v0Cuts.setSelection(ConfV0Selection.confV0TranRadMin, femto_universe_v0_selection::kV0TranRadMin, femto_universe_selection::kLowerLimit);
      v0Cuts.setSelection(ConfV0Selection.confV0TranRadMax, femto_universe_v0_selection::kV0TranRadMax, femto_universe_selection::kUpperLimit);
      v0Cuts.setSelection(ConfV0Selection.confV0DecVtxMax, femto_universe_v0_selection::kV0DecVtxMax, femto_universe_selection::kUpperLimit);
      v0Cuts.setChildCuts(femto_universe_v0_selection::kPosTrack, ConfV0Selection.confChildCharge, femto_universe_track_selection::kSign, femto_universe_selection::kEqual);
      v0Cuts.setChildCuts(femto_universe_v0_selection::kPosTrack, ConfV0Selection.confChildEtaMax, femto_universe_track_selection::kEtaMax, femto_universe_selection::kAbsUpperLimit);
      v0Cuts.setChildCuts(femto_universe_v0_selection::kPosTrack, ConfV0Selection.confChildTPCnClsMin, femto_universe_track_selection::kTPCnClsMin, femto_universe_selection::kLowerLimit);
      v0Cuts.setChildCuts(femto_universe_v0_selection::kPosTrack, ConfV0Selection.confChildDCAMin, femto_universe_track_selection::kDCAMin, femto_universe_selection::kAbsLowerLimit);
      v0Cuts.setChildCuts(femto_universe_v0_selection::kPosTrack, ConfV0Selection.confChildPIDnSigmaMax, femto_universe_track_selection::kPIDnSigmaMax, femto_universe_selection::kAbsUpperLimit);
      v0Cuts.setChildCuts(femto_universe_v0_selection::kNegTrack, ConfV0Selection.confChildCharge, femto_universe_track_selection::kSign, femto_universe_selection::kEqual);
      v0Cuts.setChildCuts(femto_universe_v0_selection::kNegTrack, ConfV0Selection.confChildEtaMax, femto_universe_track_selection::kEtaMax, femto_universe_selection::kAbsUpperLimit);
      v0Cuts.setChildCuts(femto_universe_v0_selection::kNegTrack, ConfV0Selection.confChildTPCnClsMin, femto_universe_track_selection::kTPCnClsMin, femto_universe_selection::kLowerLimit);
      v0Cuts.setChildCuts(femto_universe_v0_selection::kNegTrack, ConfV0Selection.confChildDCAMin, femto_universe_track_selection::kDCAMin, femto_universe_selection::kAbsLowerLimit);
      v0Cuts.setChildCuts(femto_universe_v0_selection::kNegTrack, ConfV0Selection.confChildPIDnSigmaMax, femto_universe_track_selection::kPIDnSigmaMax, femto_universe_selection::kAbsUpperLimit);
      v0Cuts.setChildPIDSpecies(femto_universe_v0_selection::kPosTrack, ConfV0Selection.confChildPIDspecies);
      v0Cuts.setChildPIDSpecies(femto_universe_v0_selection::kNegTrack, ConfV0Selection.confChildPIDspecies);
      v0Cuts.init<aod::femtouniverseparticle::ParticleType::kV0, aod::femtouniverseparticle::ParticleType::kV0Child, aod::femtouniverseparticle::CutContainerType>(&qaRegistry);
      v0Cuts.setInvMassLimits(ConfV0Selection.confV0InvMassLowLimit, ConfV0Selection.confV0InvMassUpLimit);

      v0Cuts.setChildRejectNotPropagatedTracks(femto_universe_v0_selection::kPosTrack, confTrkRejectNotPropagated);
      v0Cuts.setChildRejectNotPropagatedTracks(femto_universe_v0_selection::kNegTrack, confTrkRejectNotPropagated);

      v0Cuts.setnSigmaPIDOffsetTPC(confTrkPIDnSigmaOffsetTPC);
      v0Cuts.setChildnSigmaPIDOffset(femto_universe_v0_selection::kPosTrack, confTrkPIDnSigmaOffsetTPC, confTrkPIDnSigmaOffsetTOF);
      v0Cuts.setChildnSigmaPIDOffset(femto_universe_v0_selection::kNegTrack, confTrkPIDnSigmaOffsetTPC, confTrkPIDnSigmaOffsetTOF);

      if (ConfV0Selection.confV0RejectKaons) {
        v0Cuts.setKaonInvMassLimits(ConfV0Selection.confV0InvKaonMassLowLimit, ConfV0Selection.confV0InvKaonMassUpLimit);
      }
      // if (ConfRejectITSHitandTOFMissing) {
      //   o2PhysicsTrackSelection = new
      //   TrackSelection(getGlobalTrackSelection());
      //   o2PhysicsTrackSelection->SetRequireHitsInITSLayers(1, {0, 1, 2, 3});
      // }
    }

    if (ConfGeneral.confIsActivateCascade) {
      // initializing for cascades
      cascadeCuts.setSelection(ConfCascadeSelection.confCascSign, femto_universe_cascade_selection::kCascadeSign, femto_universe_selection::kEqual);
      cascadeCuts.setSelection(ConfCascadeSelection.confCascPtMin, femto_universe_cascade_selection::kCascadepTMin, femto_universe_selection::kLowerLimit);
      cascadeCuts.setSelection(ConfCascadeSelection.confCascPtMax, femto_universe_cascade_selection::kCascadepTMax, femto_universe_selection::kUpperLimit);
      cascadeCuts.setSelection(ConfCascadeSelection.confCascEtaMax, femto_universe_cascade_selection::kCascadeetaMax, femto_universe_selection::kAbsUpperLimit);
      // v0 child cuts
      cascadeCuts.setSelection(ConfCascadeSelection.confCascV0DCADaughMax, femto_universe_cascade_selection::kCascadeV0DCADaughMax, femto_universe_selection::kUpperLimit);
      cascadeCuts.setSelection(ConfCascadeSelection.confCascV0CPAMin, femto_universe_cascade_selection::kCascadeV0CPAMin, femto_universe_selection::kLowerLimit);
      cascadeCuts.setSelection(ConfCascadeSelection.confCascV0TranRadMin, femto_universe_cascade_selection::kCascadeV0TranRadMin, femto_universe_selection::kLowerLimit);
      cascadeCuts.setSelection(ConfCascadeSelection.confCascV0TranRadMax, femto_universe_cascade_selection::kCascadeV0TranRadMax, femto_universe_selection::kUpperLimit);
      cascadeCuts.setSelection(ConfCascadeSelection.confCascV0DecVtxMax, femto_universe_cascade_selection::kCascadeV0DecVtxMax, femto_universe_selection::kUpperLimit);
      // cascade cuts
      cascadeCuts.setSelection(ConfCascadeSelection.confCascDCADaughMax, femto_universe_cascade_selection::kCascadeDCADaughMax, femto_universe_selection::kUpperLimit);
      cascadeCuts.setSelection(ConfCascadeSelection.confCascCPAMin, femto_universe_cascade_selection::kCascadeCPAMin, femto_universe_selection::kLowerLimit);
      cascadeCuts.setSelection(ConfCascadeSelection.confCascTranRadMin, femto_universe_cascade_selection::kCascadeTranRadMin, femto_universe_selection::kLowerLimit);
      cascadeCuts.setSelection(ConfCascadeSelection.confCascTranRadMax, femto_universe_cascade_selection::kCascadeTranRadMax, femto_universe_selection::kUpperLimit);
      cascadeCuts.setSelection(ConfCascadeSelection.confCascDecVtxMax, femto_universe_cascade_selection::kCascadeDecVtxMax, femto_universe_selection::kUpperLimit);
      cascadeCuts.setSelection(ConfCascadeSelection.confCascDCAPosToPV, femto_universe_cascade_selection::kCascadeDCAPosToPV, femto_universe_selection::kLowerLimit);
      cascadeCuts.setSelection(ConfCascadeSelection.confCascDCANegToPV, femto_universe_cascade_selection::kCascadeDCANegToPV, femto_universe_selection::kLowerLimit);
      cascadeCuts.setSelection(ConfCascadeSelection.confCascDCABachToPV, femto_universe_cascade_selection::kCascadeDCABachToPV, femto_universe_selection::kLowerLimit);
      cascadeCuts.setSelection(ConfCascadeSelection.confCascDCAV0ToPV, femto_universe_cascade_selection::kCascadeDCAV0ToPV, femto_universe_selection::kLowerLimit);
      cascadeCuts.setSelection(ConfCascadeSelection.confCascV0MassLowLimit, femto_universe_cascade_selection::kCascadeV0MassMin, femto_universe_selection::kLowerLimit);
      cascadeCuts.setSelection(ConfCascadeSelection.confCascV0MassUpLimit, femto_universe_cascade_selection::kCascadeV0MassMax, femto_universe_selection::kUpperLimit);
      // children cuts
      cascadeCuts.setChildCuts(femto_universe_cascade_selection::kPosTrack, ConfCascadeSelection.confCascChildCharge, femto_universe_track_selection::kSign, femto_universe_selection::kEqual);
      cascadeCuts.setChildCuts(femto_universe_cascade_selection::kPosTrack, ConfCascadeSelection.confCascChildEtaMax, femto_universe_track_selection::kEtaMax, femto_universe_selection::kAbsUpperLimit);
      cascadeCuts.setChildCuts(femto_universe_cascade_selection::kPosTrack, ConfCascadeSelection.confCascChildTPCnClsMin, femto_universe_track_selection::kTPCnClsMin, femto_universe_selection::kLowerLimit);
      cascadeCuts.setChildCuts(femto_universe_cascade_selection::kPosTrack, ConfCascadeSelection.confCascChildPIDnSigmaMax, femto_universe_track_selection::kPIDnSigmaMax, femto_universe_selection::kAbsUpperLimit);
      cascadeCuts.setChildCuts(femto_universe_cascade_selection::kNegTrack, ConfCascadeSelection.confCascChildCharge, femto_universe_track_selection::kSign, femto_universe_selection::kEqual);
      cascadeCuts.setChildCuts(femto_universe_cascade_selection::kNegTrack, ConfCascadeSelection.confCascChildEtaMax, femto_universe_track_selection::kEtaMax, femto_universe_selection::kAbsUpperLimit);
      cascadeCuts.setChildCuts(femto_universe_cascade_selection::kNegTrack, ConfCascadeSelection.confCascChildTPCnClsMin, femto_universe_track_selection::kTPCnClsMin, femto_universe_selection::kLowerLimit);
      cascadeCuts.setChildCuts(femto_universe_cascade_selection::kNegTrack, ConfCascadeSelection.confCascChildPIDnSigmaMax, femto_universe_track_selection::kPIDnSigmaMax, femto_universe_selection::kAbsUpperLimit);
      cascadeCuts.setChildCuts(femto_universe_cascade_selection::kBachTrack, ConfCascadeSelection.confCascChildCharge, femto_universe_track_selection::kSign, femto_universe_selection::kEqual);
      cascadeCuts.setChildCuts(femto_universe_cascade_selection::kBachTrack, ConfCascadeSelection.confCascChildEtaMax, femto_universe_track_selection::kEtaMax, femto_universe_selection::kAbsUpperLimit);
      cascadeCuts.setChildCuts(femto_universe_cascade_selection::kBachTrack, ConfCascadeSelection.confCascChildTPCnClsMin, femto_universe_track_selection::kTPCnClsMin, femto_universe_selection::kLowerLimit);
      cascadeCuts.setChildCuts(femto_universe_cascade_selection::kBachTrack, ConfCascadeSelection.confCascChildPIDnSigmaMax, femto_universe_track_selection::kPIDnSigmaMax, femto_universe_selection::kAbsUpperLimit);

      // TODO
      cascadeCuts.setChildPIDSpecies(femto_universe_cascade_selection::kPosTrack, ConfV0Selection.confChildPIDspecies);
      cascadeCuts.setChildPIDSpecies(femto_universe_cascade_selection::kNegTrack, ConfV0Selection.confChildPIDspecies);
      cascadeCuts.setChildPIDSpecies(femto_universe_cascade_selection::kBachTrack, ConfCascadeSelection.confCascChildPIDspecies);

      // check if works correctly for bachelor track
      cascadeCuts.init<aod::femtouniverseparticle::ParticleType::kCascade, aod::femtouniverseparticle::ParticleType::kV0Child, aod::femtouniverseparticle::ParticleType::kCascadeBachelor, aod::femtouniverseparticle::CutContainerType>(&cascadeQaRegistry);
      // invmass cuts
      cascadeCuts.setInvMassLimits(ConfCascadeSelection.confXiInvMassLowLimit, ConfCascadeSelection.confOmegaInvMassLowLimit, ConfCascadeSelection.confXiInvMassUpLimit, ConfCascadeSelection.confOmegaInvMassUpLimit);
    }

    if (ConfGeneral.confIsActivatePhi) {
      // initializing for Phi meson
      phiCuts.init<aod::femtouniverseparticle::ParticleType::kPhi, aod::femtouniverseparticle::ParticleType::kPhiChild, aod::femtouniverseparticle::CutContainerType>(&qaRegistry);
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

  /// Function to retrieve the nominal magnetic field in kG (0.1T) and convert it directly to T
  void getMagneticFieldTesla(aod::BCsWithTimestamps::iterator bc)
  {
    // TODO done only once (and not per run). Will be replaced by CCDBConfigurable
    // get magnetic field for run
    if (mRunNumber == bc.runNumber())
      return;
    auto timestamp = bc.timestamp();
    float output = -999;

    if (confIsRun3 && !confIsForceGRP) {
      static o2::parameters::GRPMagField* grpo = nullptr;
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with L3 ", timestamp, grpo->getL3Current());
      // taken from GRP onject definition of getnominalL3Field; update later to something smarter (mnominalL3Field = std::lround(5.f * mL3Current / 30000.f);)
      auto nominalL3Field = std::lround(5.f * grpo->getL3Current() / 30000.f);
      output = 0.1 * (nominalL3Field);

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

  template <bool isTrackOrV0, bool isPhiOrD0, bool isCasc, typename ParticleType>
  void fillDebugParticle(ParticleType const& particle)
  {
    if constexpr (isTrackOrV0) {
      if constexpr (std::experimental::is_detected<hasStrangeTOFinV0, ParticleType>::value) {
        outputDebugParts(particle.sign(), (uint8_t)particle.tpcNClsFound(),
                         particle.tpcNClsFindable(),
                         (uint8_t)particle.tpcNClsCrossedRows(),
                         particle.tpcNClsShared(), particle.tpcFractionSharedCls(), particle.tpcInnerParam(),
                         particle.itsNCls(), particle.itsNClsInnerBarrel(),
                         particle.dcaXY(), particle.dcaZ(), particle.tpcSignal(),
                         particle.tpcNSigmaStoreEl(), particle.tpcNSigmaStorePi(),
                         particle.tpcNSigmaStoreKa(), particle.tpcNSigmaStorePr(),
                         particle.tpcNSigmaStoreDe(), particle.tofNSigmaStoreEl(),
                         particle.tofNSigmaStorePi(), particle.tofNSigmaStoreKa(),
                         particle.tofNSigmaStorePr(), particle.tofNSigmaStoreDe(),
                         particle.tofNSigmaLaPr(), particle.tofNSigmaLaPi(), particle.tofNSigmaALaPr(),
                         particle.tofNSigmaALaPi(), particle.tofNSigmaK0PiPlus(), particle.tofNSigmaK0PiMinus());
      } else {
        outputDebugParts(particle.sign(), (uint8_t)particle.tpcNClsFound(),
                         particle.tpcNClsFindable(),
                         (uint8_t)particle.tpcNClsCrossedRows(),
                         particle.tpcNClsShared(), particle.tpcFractionSharedCls(), particle.tpcInnerParam(),
                         particle.itsNCls(), particle.itsNClsInnerBarrel(),
                         particle.dcaXY(), particle.dcaZ(), particle.tpcSignal(),
                         particle.tpcNSigmaStoreEl(), particle.tpcNSigmaStorePi(),
                         particle.tpcNSigmaStoreKa(), particle.tpcNSigmaStorePr(),
                         particle.tpcNSigmaStoreDe(), particle.tofNSigmaStoreEl(),
                         particle.tofNSigmaStorePi(), particle.tofNSigmaStoreKa(),
                         particle.tofNSigmaStorePr(), particle.tofNSigmaStoreDe(),
                         -999., -999., -999., -999., -999., -999.);
      }
    } else if constexpr (isPhiOrD0) {
      outputDebugParts(-999., -999., -999., -999., -999., -999., -999., -999., -999.,
                       -999., -999., -999., -999., -999., -999., -999., -999.,
                       -999., -999., -999., -999., -999.,
                       -999., -999.,
                       -999., -999., -999.,
                       -999.); // QA for phi or D0/D0bar children
    } else if constexpr (isCasc) {
      if constexpr (std::experimental::is_detected<hasStrangeTOFinCasc, ParticleType>::value) {
        outputDebugParts(-999., -999., -999., -999., -999., -999., -999., -999., -999.,
                         -999., -999., -999., -999., -999., -999., -999., particle.tofNSigmaXiLaPi(),
                         particle.tofNSigmaXiLaPr(), particle.tofNSigmaXiPi(), particle.tofNSigmaOmLaPi(),
                         particle.tofNSigmaOmLaPr(), particle.tofNSigmaOmKa(),
                         particle.dcacascdaughters(), particle.cascradius(),
                         particle.x(), particle.y(), particle.z(), -999.);
      } else {
        outputDebugParts(-999., -999., -999., -999., -999., -999., -999., -999., -999.,
                         -999., -999., -999., -999., -999., -999., -999., -999.,
                         -999., -999., -999., -999., -999.,
                         particle.dcacascdaughters(), particle.cascradius(),
                         particle.x(), particle.y(), particle.z(), -999.);
      }
    } else {
      // LOGF(info, "isTrack0orV0: %d, isPhi: %d", isTrackOrV0, isPhiOrD0);
      outputDebugParts(-999., -999., -999., -999., -999., -999., -999., -999., -999.,
                       particle.dcav0topv(), -999., -999., -999., -999., -999., -999., -999.,
                       -999., -999., -999., -999., -999.,
                       particle.dcaV0daughters(), particle.v0radius(),
                       particle.x(), particle.y(), particle.z(),
                       particle.mK0Short()); // QA for v0
    }
  }

  template <bool isD0ML, bool isD0barML, typename ParticleType>
  void fillDebugD0D0barML(ParticleType const& particle)
  {
    if constexpr (isD0ML) {
      outputDebugParts(-999., -999., -999., -999., -999., -999., -999., -999., -999.,
                       -999., -999., -999., -999., -999., -999., -999., -999.,
                       -999., -999., -999., -999., -999.,
                       -999.,
                       hfHelper.yD0(particle), // getter transRadius
                       particle.mlProbD0()[0], // getter decayVtxX
                       particle.mlProbD0()[1], // getter decayVtxY
                       particle.mlProbD0()[2], // getter decayVtxZ
                       -999.);                 // Additional info for D0/D0bar
    } else if constexpr (isD0barML) {
      outputDebugParts(-999., -999., -999., -999., -999., -999., -999., -999., -999.,
                       -999., -999., -999., -999., -999., -999., -999., -999.,
                       -999., -999., -999., -999., -999.,
                       -999.,
                       hfHelper.yD0(particle),    // getter transRadius
                       particle.mlProbD0bar()[0], // getter decayVtxX
                       particle.mlProbD0bar()[1], // getter decayVtxY
                       particle.mlProbD0bar()[2], // getter decayVtxZ
                       -999.);                    // Additional info for D0/D0bar
    } else {
      outputDebugParts(-999., -999., -999., -999., -999., -999., -999., -999., -999.,
                       -999., -999., -999., -999., -999., -999., -999., -999.,
                       -999., -999., -999., -999., -999.,
                       -999., -999., -999., -999., -999., -999.);
    }
  }

  template <bool isD0ML, bool isD0barML, typename ParticleType>
  void fillDebugD0D0barMcMl(ParticleType const& particle)
  {
    int8_t originMcReco = 2; // 0 - prompt, 1 - non-prompt, 2 - default/else
    if (particle.originMcRec() == RecoDecay::OriginType::Prompt) {
      originMcReco = 0;
    } else if (particle.originMcRec() == RecoDecay::OriginType::NonPrompt) {
      originMcReco = 1;
    } else {
      originMcReco = 2;
    }
    if constexpr (isD0ML) {
      outputDebugParts(particle.flagMcMatchRec(), // getter sign
                       originMcReco, -999., -999., -999., -999., -999., -999., -999.,
                       -999., -999., -999., -999., -999., -999., -999., -999.,
                       -999., -999., -999., -999., -999.,
                       -999.,
                       hfHelper.yD0(particle), // getter transRadius
                       particle.mlProbD0()[0], // getter decayVtxX
                       particle.mlProbD0()[1], // getter decayVtxY
                       particle.mlProbD0()[2], // getter decayVtxZ
                       -999.);                 // Additional info for D0/D0bar
    } else if constexpr (isD0barML) {
      outputDebugParts(particle.flagMcMatchRec(), -999., -999., -999., -999., -999., -999., -999., -999.,
                       originMcReco, -999., -999., -999., -999., -999., -999., -999.,
                       -999., -999., -999., -999., -999.,
                       -999.,
                       hfHelper.yD0(particle),    // getter transRadius
                       particle.mlProbD0bar()[0], // getter decayVtxX
                       particle.mlProbD0bar()[1], // getter decayVtxY
                       particle.mlProbD0bar()[2], // getter decayVtxZ
                       -999.);                    // Additional info for D0/D0bar
    } else {
      outputDebugParts(-999., -999., -999., -999., -999., -999., -999., -999., -999.,
                       -999., -999., -999., -999., -999., -999., -999., -999.,
                       -999., -999., -999., -999., -999.,
                       -999., -999., -999., -999., -999., -999.);
    }
  }

  template <typename ParticleType>
  int32_t getMotherPDG(ParticleType particle)
  {
    auto motherparticlesMC = particle.template mothers_as<aod::McParticles>();
    if (!motherparticlesMC.empty()) {
      auto motherparticleMC = motherparticlesMC.front();
      return particle.isPhysicalPrimary() ? 0 : motherparticleMC.pdgCode();
    } else {
      return 9999;
    }
  }

  template <typename ParticleType>
  void fillDebugParticleMC(ParticleType const& particle)
  {
    outputDebugPartsMC(getMotherPDG(particle));
  }

  template <typename ParticleType>
  void fillMCParticle(ParticleType const& particle, o2::aod::femtouniverseparticle::ParticleType fdparttype)
  {
    if (particle.has_mcParticle()) {
      // get corresponding MC particle and its info
      auto particleMC = particle.mcParticle();
      auto pdgCode = particleMC.pdgCode();
      int particleOrigin = 99;
      auto motherparticlesMC = particleMC.template mothers_as<aod::McParticles>();

      if (std::abs(pdgCode) == std::abs(confPDGCodePartOne.value) || std::abs(pdgCode) == std::abs(confPDGCodePartTwo.value)) {
        if (particleMC.isPhysicalPrimary()) {
          particleOrigin = aod::femtouniverse_mc_particle::ParticleOriginMCTruth::kPrimary;
        } else if (!motherparticlesMC.empty()) {
          auto motherparticleMC = motherparticlesMC.front();
          if (motherparticleMC.producedByGenerator()) {
            particleOrigin = checkDaughterType(fdparttype, motherparticleMC.pdgCode());
          } else {
            particleOrigin = aod::femtouniverse_mc_particle::ParticleOriginMCTruth::kMaterial;
          }
        }
      } else {
        particleOrigin = aod::femtouniverse_mc_particle::ParticleOriginMCTruth::kFake;
      }

      outputPartsMC(particleOrigin, pdgCode, particleMC.pt(), particleMC.eta(), particleMC.phi());
      fillDebugParticleMC(particleMC);
      outputPartsMCLabels(outputPartsMC.lastIndex());
    } else {
      outputPartsMCLabels(-1);
      outputDebugPartsMC(9999);
    }
  }

  template <typename MCParticleType>
  void fillMCTruthParticle(MCParticleType const& mcparticle, o2::aod::femtouniverseparticle::ParticleType fdparttype)
  {
    auto pdgCode = mcparticle.pdgCode();
    int particleOrigin = 99;

    // Determine particle origin
    if (std::abs(pdgCode) == std::abs(confPDGCodePartOne.value) || std::abs(pdgCode) == std::abs(confPDGCodePartTwo.value)) {
      if (mcparticle.isPhysicalPrimary()) {
        particleOrigin = aod::femtouniverse_mc_particle::ParticleOriginMCTruth::kPrimary;
      } else {
        auto mothers = mcparticle.template mothers_as<aod::McParticles>();
        if (!mothers.empty()) {
          auto mother = mothers.front();
          if (mother.producedByGenerator()) {
            particleOrigin = checkDaughterType(fdparttype, mother.pdgCode());
          } else {
            particleOrigin = aod::femtouniverse_mc_particle::ParticleOriginMCTruth::kMaterial;
          }
        }
      }
    } else {
      particleOrigin = aod::femtouniverse_mc_particle::ParticleOriginMCTruth::kFake;
    }

    // Fill MC companion tables
    outputPartsMC(particleOrigin, pdgCode, mcparticle.pt(), mcparticle.eta(), mcparticle.phi());
    fillDebugParticleMC(mcparticle);
    outputPartsMCLabels(outputPartsMC.lastIndex());

    // Artificial fill of a debug table -- solves conflicts between FDParticles and FDExtParticles tables
    outputDebugParts(-111., -111., -111., -111., -111., -111., -111., -111., -111.,
                     -111., -111., -111., -111., -111., -111., -111., -111.,
                     -111., -111., -111., -111., -111.,
                     -111., -111., -111., -111., -111., -111.);
  }

  template <typename ParticleType>
  void fillMCParticlePhi(ParticleType const& kaon1, ParticleType const& kaon2)
  {
    if (kaon1.has_mcParticle() && kaon2.has_mcParticle()) {
      // get corresponding MC particle and its info
      auto kaon1MC = kaon1.mcParticle();
      auto kaon2MC = kaon2.mcParticle();
      auto pdgCode1 = kaon1MC.pdgCode();
      auto pdgCode2 = kaon2MC.pdgCode();

      int phiOrigin = 99;
      auto motherskaon1MC = kaon1MC.template mothers_as<aod::McParticles>();
      auto motherskaon2MC = kaon2MC.template mothers_as<aod::McParticles>();

      if (std::abs(pdgCode1) == std::abs(321) || std::abs(pdgCode2) == std::abs(-321)) {
        if ((kaon1MC.isPhysicalPrimary() && kaon2MC.isPhysicalPrimary()) && (!motherskaon1MC.empty() && !motherskaon2MC.empty())) {
          for (const auto& particleMotherOfNeg : motherskaon1MC) {
            for (const auto& particleMotherOfPos : motherskaon2MC) {
              if (particleMotherOfNeg == particleMotherOfPos && particleMotherOfNeg.pdgCode() == Pdg::kPhi) {
                phiOrigin = aod::femtouniverse_mc_particle::ParticleOriginMCTruth::kPrimary;
              } else {
                phiOrigin = aod::femtouniverse_mc_particle::ParticleOriginMCTruth::kFake;
              }
            }
          }
        } else {
          phiOrigin = aod::femtouniverse_mc_particle::ParticleOriginMCTruth::kFake;
        }
      } else {
        phiOrigin = aod::femtouniverse_mc_particle::ParticleOriginMCTruth::kFake;
      }

      TLorentzVector part1Vec;
      TLorentzVector part2Vec;

      const auto mMassOne = o2::constants::physics::MassKPlus;  // FIXME: Get from the PDG service of the common header
      const auto mMassTwo = o2::constants::physics::MassKMinus; // FIXME: Get from the PDG service of the common header

      part1Vec.SetPtEtaPhiM(kaon1MC.pt(), kaon1MC.eta(), kaon1MC.phi(), mMassOne);
      part2Vec.SetPtEtaPhiM(kaon2MC.pt(), kaon2MC.eta(), kaon2MC.phi(), mMassTwo);

      TLorentzVector sumVec(part1Vec);
      sumVec += part2Vec;

      float phiEta = sumVec.Eta();
      float phiPt = sumVec.Pt();
      float phiPhi = sumVec.Phi();

      outputPartsMC(phiOrigin, 333, phiPt, phiEta, phiPhi);
      outputPartsMCLabels(outputPartsMC.lastIndex());
    } else {
      outputPartsMCLabels(-1);
    }
  }

  template <bool isMC, typename CollisionType, typename TrackType>
  bool fillCollisions(CollisionType const& col, TrackType const& tracks)
  {
    const auto vtxZ = col.posZ();
    float mult = 0;
    int multNtr = 0;
    if (confIsRun3) {
      mult = col.multFV0M();
      multNtr = col.multNTracksPV();
    } else {
      mult = 0.5 * (col.multFV0M()); /// For benchmarking on Run 2, V0M in
                                     /// FemtoUniverseRun2 is defined V0M/2
      multNtr = col.multTracklets();
    }
    if (ConfGeneral.confEvtUseTPCmult) {
      multNtr = col.multTPC();
    }

    // check whether the basic event selection criteria are fulfilled
    // if the basic selection is NOT fulfilled:
    // in case of skimming run - don't store such collisions
    // in case of trigger run - store such collisions but don't store any
    // particle candidates for such collisions
    if (!colCuts.isSelected(col)) {
      return false;
    }

    if (zorroMask.value != "") {
      auto bc = col.template bc_as<aod::BCsWithTimestamps>();
      initZorro(bc);
      if (!zorro.isSelected(col.template bc_as<aod::BCsWithTimestamps>().globalBC()))
        return false;
    }

    if (!ConfGeneral.confIsUsePileUp) {
      outputCollision(vtxZ, mult, multNtr, confDoSpher ? colCuts.computeSphericity(col, tracks) : 2, mMagField);
      colCuts.fillQA(col);
      return true;
    } else if ((!ConfGeneral.confEvNoSameBunchPileup || col.selection_bit(aod::evsel::kNoSameBunchPileup)) && (!ConfGeneral.confEvIsGoodZvtxFT0vsPV || col.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) && (!ConfGeneral.confEvIsVertexITSTPC || col.selection_bit(aod::evsel::kIsVertexITSTPC))) {
      outputCollision(vtxZ, mult, multNtr, confDoSpher ? colCuts.computeSphericity(col, tracks) : 2, mMagField);
      colCuts.fillQA(col);
      return true;
    } else {
      return false;
    }
  }

  template <bool isMC, typename CollisionType, typename TrackType>
  bool fillCollisionsCentPP(CollisionType const& col, TrackType const& tracks)
  {
    const auto vtxZ = col.posZ();
    float mult = 0;
    int multNtr = 0;
    if (confIsRun3) {
      mult = col.centFT0M();
      multNtr = col.multNTracksPV();
    } else {
      mult = 0.5 * (col.multFV0M()); /// For benchmarking on Run 2, V0M in
                                     /// FemtoUniverseRun2 is defined V0M/2
      multNtr = col.multTracklets();
    }
    if (ConfGeneral.confEvtUseTPCmult) {
      multNtr = col.multTPC();
    }

    // check whether the basic event selection criteria are fulfilled
    // if the basic selection is NOT fulfilled:
    // in case of skimming run - don't store such collisions
    // in case of trigger run - store such collisions but don't store any
    // particle candidates for such collisions
    if (!colCuts.isSelected(col)) {
      return false;
    }
    if (!ConfGeneral.confIsUsePileUp) {
      outputCollision(vtxZ, mult, multNtr, confDoSpher ? colCuts.computeSphericity(col, tracks) : 2, mMagField);
      colCuts.fillQA(col);
      return true;
    } else if ((!ConfGeneral.confEvNoSameBunchPileup || col.selection_bit(aod::evsel::kNoSameBunchPileup)) && (!ConfGeneral.confEvIsGoodZvtxFT0vsPV || col.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) && (!ConfGeneral.confEvIsVertexITSTPC || col.selection_bit(aod::evsel::kIsVertexITSTPC))) {
      outputCollision(vtxZ, mult, multNtr, confDoSpher ? colCuts.computeSphericity(col, tracks) : 2, mMagField);
      colCuts.fillQA(col);
      return true;
    } else {
      return false;
    }
  }

  template <typename CollisionType, typename TrackType>
  void fillMCTruthCollisions(CollisionType const& col, TrackType const& tracks)
  {
    for (const auto& c : col) {
      const auto vtxZ = c.posZ();
      float mult = confIsRun3 ? c.multFV0M() : 0.5 * (c.multFV0M());
      int multNtr = confIsRun3 ? c.multNTracksPV() : c.multTracklets();

      if (std::abs(vtxZ) > ConfGeneral.confEvtZvtx) {
        continue;
      }

      if (confDoSpher) {
        outputCollision(vtxZ, mult, multNtr, colCuts.computeSphericity(col, tracks), mMagField);
      } else {
        outputCollision(vtxZ, mult, multNtr, 2, mMagField);
      }
    }
  }

  template <bool isMC, typename CollisionType>
  bool fillCollisionsCentRun2(CollisionType const& col)
  {
    const auto vtxZ = col.posZ();
    const auto cent = col.centRun2V0M();
    const auto multNtr = col.multNTracksPV();

    // check whether the basic event selection criteria are fulfilled
    // if the basic selection is NOT fulfilled:
    // in case of skimming run - don't store such collisions
    // in case of trigger run - store such collisions but don't store any
    // particle candidates for such collisions
    if (!colCuts.isSelected(col)) {
      return false;
    } else {
      outputCollision(vtxZ, cent, multNtr, 2, mMagField);
      return true;
    }
  }

  template <bool isMC, typename CollisionType>
  bool fillCollisionsCentRun3(CollisionType const& col)
  {
    const auto vtxZ = col.posZ();
    const auto multNtr = col.multNTracksPV();
    const auto cent = col.centFT0C();
    const auto occupancy = col.trackOccupancyInTimeRange();

    // check whether the basic event selection criteria are fulfilled
    // if the basic selection is NOT fulfilled:
    // in case of skimming run - don't store such collisions
    // in case of trigger run - store such collisions but don't store any
    // particle candidates for such collisions

    if (!colCuts.isSelectedRun3(col) || (occupancy < ConfGeneral.confTPCOccupancyMin || occupancy > ConfGeneral.confTPCOccupancyMax)) {
      return false;
    } else {
      if ((!ConfGeneral.confEvNoSameBunchPileup || col.selection_bit(aod::evsel::kNoSameBunchPileup)) &&
          (!ConfGeneral.confEvIsGoodZvtxFT0vsPV || col.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) &&
          (!ConfGeneral.confIsGoodITSLayersAll || col.selection_bit(aod::evsel::kIsGoodITSLayersAll)) &&
          (!ConfGeneral.confNoCollInRofStandard || col.selection_bit(aod::evsel::kNoCollInRofStandard)) &&
          (!ConfGeneral.confNoHighMultCollInPrevRof || col.selection_bit(aod::evsel::kNoHighMultCollInPrevRof)) &&
          (!ConfGeneral.confEvIsVertexITSTPC || col.selection_bit(aod::evsel::kIsVertexITSTPC)) &&
          (!ConfGeneral.confNoCollInTimeRangeStandard || col.selection_bit(aod::evsel::kNoCollInTimeRangeStandard)) &&
          (!ConfGeneral.confNoITSROFrameBorder || col.selection_bit(aod::evsel::kNoITSROFrameBorder)) &&
          (!ConfGeneral.confNoTimeFrameBorder || col.selection_bit(aod::evsel::kNoTimeFrameBorder))) {
        outputCollision(vtxZ, cent, multNtr, 2, mMagField);
        return true;
      } else {
        return false;
      }
    }
  }

  template <typename CollisionType>
  bool fillMCTruthCollisionsCentRun3(CollisionType const& col)
  {
    const auto vtxZ = col.posZ();

    if (std::abs(vtxZ) > ConfGeneral.confEvtZvtx) {
      return false;
    } else {
      outputCollision(vtxZ, 0, 0, 2, mMagField);
      return true;
    }
  }

  template <bool isMC, typename CollisionType>
  void fillCollisionsCentRun3ColExtra(CollisionType const& col, double irrate)
  {
    const auto occupancy = col.trackOccupancyInTimeRange();
    outputCollExtra(irrate, occupancy);
  }

  template <bool isMC, typename TrackType>
  void fillTracks(TrackType const& tracks)
  {
    std::vector<int> childIDs = {0, 0}; // these IDs are necessary to keep track of the children
    std::vector<int> tmpIDtrack;        // this vector keeps track of the matching of the primary track table row <-> aod::track table global index

    for (const auto& track : tracks) {
      /// if the most open selection criteria are not fulfilled there is no
      /// point looking further at the track
      if (!trackCuts.isSelectedMinimal(track)) {
        continue;
      }

      if (track.pt() > confTOFpTmin) {
        if (!track.hasTOF()) {
          continue;
        }
      }

      trackCuts.fillQA<aod::femtouniverseparticle::ParticleType::kTrack,
                       aod::femtouniverseparticle::TrackType::kNoChild>(track);
      // the bit-wise container of the systematic variations is obtained
      auto cutContainer = trackCuts.getCutContainer<aod::femtouniverseparticle::CutContainerType>(track);

      // now the table is filled
      outputParts(outputCollision.lastIndex(), track.pt(), track.eta(),
                  track.phi(), aod::femtouniverseparticle::ParticleType::kTrack,
                  cutContainer.at(
                    femto_universe_track_selection::TrackContainerPosition::kCuts),
                  confIsUseCutculator ? cutContainer.at(
                                          femto_universe_track_selection::TrackContainerPosition::kPID)
                                      : PIDBitmask(track),
                  track.dcaXY(), childIDs,
                  track.hasTOF(), // hasTOF getter is mLambda()
                  track.sign());  // sign getter is mAntiLambda()

      tmpIDtrack.push_back(track.globalIndex());
      if (confIsDebug) {
        fillDebugParticle<true, false, false>(track);
      }

      if constexpr (isMC) {
        fillMCParticle(track, o2::aod::femtouniverseparticle::ParticleType::kTrack);
      }
    }
  }

  template <bool isMC, typename CollisionType, typename V0Type, typename TrackType>
  void fillV0(CollisionType const& col, V0Type const& fullV0s, TrackType const&)
  {
    std::vector<int> childIDs = {0, 0}; // these IDs are necessary to keep track of the children
    std::vector<int> tmpIDtrack;        // this vector keeps track of the matching of the primary track table row <-> aod::track table global index
    for (const auto& v0 : fullV0s) {
      auto postrack = v0.template posTrack_as<TrackType>();
      auto negtrack = v0.template negTrack_as<TrackType>();
      ///\tocheck funnily enough if we apply the filter the
      /// sign of Pos and Neg track is always negative
      // const auto dcaXYpos = postrack.dcaXY();
      // const auto dcaZpos = postrack.dcaZ();
      // const auto dcapos = std::sqrt(pow(dcaXYpos, 2.) + pow(dcaZpos, 2.));
      v0Cuts.fillLambdaQA(col, v0, postrack, negtrack);

      if (!v0Cuts.isSelectedMinimal(col, v0, postrack, negtrack)) {
        continue;
      }

      // if (ConfRejectITSHitandTOFMissing) {
      // Uncomment only when TOF timing is solved
      // bool itsHit = o2PhysicsTrackSelection->IsSelected(postrack,
      // TrackSelection::TrackCuts::kITSHits); bool itsHit =
      // o2PhysicsTrackSelection->IsSelected(negtrack,
      // TrackSelection::TrackCuts::kITSHits);
      // }

      v0Cuts.fillQA<aod::femtouniverseparticle::ParticleType::kV0, aod::femtouniverseparticle::ParticleType::kV0Child>(col, v0, postrack, negtrack); ///\todo fill QA also for daughters
      auto cutContainerV0 = v0Cuts.getCutContainer<aod::femtouniverseparticle::CutContainerType>(col, v0, postrack, negtrack);

      int postrackID = v0.posTrackId();
      int rowInPrimaryTrackTablePos = -1;
      rowInPrimaryTrackTablePos = getRowDaughters(postrackID, tmpIDtrack);
      childIDs[0] = rowInPrimaryTrackTablePos;
      childIDs[1] = 0;
      outputParts(outputCollision.lastIndex(), v0.positivept(),
                  v0.positiveeta(), v0.positivephi(),
                  aod::femtouniverseparticle::ParticleType::kV0Child,
                  cutContainerV0.at(femto_universe_v0_selection::V0ContainerPosition::kPosCuts),
                  confIsUseCutculator ? cutContainerV0.at(femto_universe_v0_selection::V0ContainerPosition::kPosPID) : PIDBitmask(postrack),
                  0.,
                  childIDs,
                  0,
                  postrack.sign());
      const int rowOfPosTrack = outputParts.lastIndex();
      if constexpr (isMC) {
        fillMCParticle(postrack, o2::aod::femtouniverseparticle::ParticleType::kV0Child);
      }
      int negtrackID = v0.negTrackId();
      int rowInPrimaryTrackTableNeg = -1;
      rowInPrimaryTrackTableNeg = getRowDaughters(negtrackID, tmpIDtrack);
      childIDs[0] = 0;
      childIDs[1] = rowInPrimaryTrackTableNeg;
      outputParts(outputCollision.lastIndex(),
                  v0.negativept(),
                  v0.negativeeta(),
                  v0.negativephi(),
                  aod::femtouniverseparticle::ParticleType::kV0Child,
                  cutContainerV0.at(femto_universe_v0_selection::V0ContainerPosition::kNegCuts),
                  confIsUseCutculator ? cutContainerV0.at(femto_universe_v0_selection::V0ContainerPosition::kNegPID) : PIDBitmask(negtrack),
                  0.,
                  childIDs,
                  0,
                  negtrack.sign());
      const int rowOfNegTrack = outputParts.lastIndex();
      if constexpr (isMC) {
        fillMCParticle(negtrack, o2::aod::femtouniverseparticle::ParticleType::kV0Child);
      }
      std::vector<int> indexChildID = {rowOfPosTrack, rowOfNegTrack};
      outputParts(outputCollision.lastIndex(),
                  v0.pt(),
                  v0.eta(),
                  v0.phi(),
                  aod::femtouniverseparticle::ParticleType::kV0,
                  cutContainerV0.at(femto_universe_v0_selection::V0ContainerPosition::kV0),
                  PIDStrangeTOFBitmaskV0(v0),
                  v0.v0cosPA(),
                  indexChildID,
                  v0.mLambda(),
                  v0.mAntiLambda());
      if (confIsDebug) {
        fillDebugParticle<true, false, false>(postrack); // QA for positive daughter
        fillDebugParticle<true, false, false>(negtrack); // QA for negative daughter
        fillDebugParticle<false, false, false>(v0);      // QA for v0
      }
      if constexpr (isMC) {
        fillMCParticle(v0, o2::aod::femtouniverseparticle::ParticleType::kV0);
      }
    }
  }

  template <typename MCParticlesType>
  void fillV0MCTruth(MCParticlesType const& mcParticles)
  {
    for (const auto& mc : mcParticles) { // Loop over all MC Truth particles
      if (mc.pdgCode() != ConfV0Selection.confV0PDGMCTruth.value[2])
        continue; // Artificially single out V0s

      auto daughters = mc.template daughters_as<aod::McParticles>(); // Access daughters (no differentiation of signs, it needs to be checked separately)

      bool foundPos = false, foundNeg = false;

      std::vector<int> childIDs = {0, 0};
      int rowPos = 0;
      int rowNeg = 0;

      for (auto const& d : daughters) {                                 // Loop over daughters
        if (d.pdgCode() == ConfV0Selection.confV0PDGMCTruth.value[0]) { // Check for a positive child
          foundPos = true;

          outputParts(outputCollision.lastIndex(),
                      d.pt(),
                      d.eta(),
                      d.phi(),
                      aod::femtouniverseparticle::ParticleType::kMCTruthTrack,
                      0,
                      0,
                      d.pdgCode(),
                      childIDs, // {0, 0}
                      0,
                      0);
          rowPos = outputParts.lastIndex();
          fillMCTruthParticle(d, aod::femtouniverseparticle::ParticleType::kMCTruthTrack);
        } else if (d.pdgCode() == ConfV0Selection.confV0PDGMCTruth.value[1]) { // Check for a negative child
          foundNeg = true;

          outputParts(outputCollision.lastIndex(),
                      d.pt(),
                      d.eta(),
                      d.phi(),
                      aod::femtouniverseparticle::ParticleType::kMCTruthTrack,
                      0,
                      0,
                      d.pdgCode(),
                      childIDs, // {0, 0}
                      0,
                      0);
          rowNeg = outputParts.lastIndex();
          fillMCTruthParticle(d, aod::femtouniverseparticle::ParticleType::kMCTruthTrack);
        }
      }

      if (!foundPos || !foundNeg)
        continue;

      childIDs[0] = rowPos;
      childIDs[1] = rowNeg;
      outputParts(outputCollision.lastIndex(),
                  mc.pt(),
                  mc.eta(),
                  mc.phi(),
                  aod::femtouniverseparticle::ParticleType::kMCTruthTrack,
                  0,
                  0,
                  mc.pdgCode(),
                  childIDs,
                  0,
                  0);
      fillMCTruthParticle(mc, aod::femtouniverseparticle::ParticleType::kMCTruthTrack);
    }
  }

  template <typename MCParticlesType>
  void fillTracksMCTruth(MCParticlesType const& mcParticles)
  {
    for (const auto& mc : mcParticles) { // Loop over all MC Truth particles
      if (ConfFilterCuts.confIsApplyTrkCutMCTruth) {
        if (std::abs(mc.eta()) > ConfFilterCuts.confEtaFilterCut || mc.pt() < ConfFilterCuts.confPtLowFilterCut || mc.pt() > ConfFilterCuts.confPtHighFilterCut) {
          continue;
        }
      }

      if (ConfFilterCuts.confIsOnlyPrimary) {
        if (!mc.isPhysicalPrimary()) {
          continue;
        }
      }

      std::vector<int> childIDs = {0, 0};
      outputParts(outputCollision.lastIndex(),
                  mc.pt(),
                  mc.eta(),
                  mc.phi(),
                  aod::femtouniverseparticle::ParticleType::kMCTruthTrack,
                  0,
                  0,
                  mc.pdgCode(),
                  childIDs,
                  0,
                  0);
      fillMCTruthParticle(mc, aod::femtouniverseparticle::ParticleType::kMCTruthTrack);
    }
  }

  template <bool isMC, typename CollisionType, typename CascadeType, typename TrackType>
  void fillCascade(CollisionType const& col, CascadeType const& fullCascades, TrackType const&)
  {
    std::vector<int> childIDs = {0, 0, 0}; // child1, child2, bachelor; these IDs are necessary to keep track of the children
    std::vector<int> tmpIDtrack;           // this vector keeps track of the matching of the primary track table row <-> aod::track table global index

    for (const auto& casc : fullCascades) {
      const auto& posTrackCasc = casc.template posTrack_as<TrackType>();
      const auto& negTrackCasc = casc.template negTrack_as<TrackType>();
      const auto& bachTrackCasc = casc.template bachelor_as<TrackType>();

      if (!cascadeCuts.isSelectedMinimal(col, casc, posTrackCasc, negTrackCasc, bachTrackCasc)) {
        continue;
      }

      cascadeCuts.fillCascadeQA(col, casc, posTrackCasc, negTrackCasc);
      cascadeCuts.fillQA(col, casc, posTrackCasc, negTrackCasc, bachTrackCasc); // fill QA for daughters

      int postrackID = casc.posTrackId();
      int rowInPrimaryTrackTablePos = -1;
      rowInPrimaryTrackTablePos = getRowDaughters(postrackID, tmpIDtrack);
      childIDs[0] = rowInPrimaryTrackTablePos; // pos
      childIDs[1] = 0;                         // neg
      childIDs[2] = 0;                         // bachelor
      float hasTOF = posTrackCasc.hasTOF() ? 1 : 0;
      outputParts(outputCollision.lastIndex(),
                  casc.positivept(),
                  casc.positiveeta(),
                  casc.positivephi(),
                  aod::femtouniverseparticle::ParticleType::kV0Child,
                  0,
                  PIDBitmask(posTrackCasc),
                  hasTOF,
                  childIDs,
                  0,
                  posTrackCasc.sign());
      const int rowOfPosTrack = outputParts.lastIndex();
      if constexpr (isMC) {
        fillMCParticle(posTrackCasc, o2::aod::femtouniverseparticle::ParticleType::kV0Child);
      }
      int negtrackID = casc.negTrackId();
      int rowInPrimaryTrackTableNeg = -1;
      rowInPrimaryTrackTableNeg = getRowDaughters(negtrackID, tmpIDtrack);
      childIDs[0] = 0;                         // pos
      childIDs[1] = rowInPrimaryTrackTableNeg; // neg
      childIDs[2] = 0;                         // bachelor
      hasTOF = negTrackCasc.hasTOF() ? 1 : 0;
      outputParts(outputCollision.lastIndex(),
                  casc.negativept(),
                  casc.negativeeta(),
                  casc.negativephi(),
                  aod::femtouniverseparticle::ParticleType::kV0Child,
                  0,
                  PIDBitmask(negTrackCasc),
                  hasTOF,
                  childIDs,
                  0,
                  negTrackCasc.sign());
      const int rowOfNegTrack = outputParts.lastIndex();
      if constexpr (isMC) {
        fillMCParticle(negTrackCasc, o2::aod::femtouniverseparticle::ParticleType::kV0Child);
      }
      //  bachelor
      int bachtrackID = casc.bachelorId();
      int rowInPrimaryTrackTableBach = -1;
      rowInPrimaryTrackTableBach = getRowDaughters(bachtrackID, tmpIDtrack);
      childIDs[0] = 0;                          // pos
      childIDs[1] = 0;                          // neg
      childIDs[2] = rowInPrimaryTrackTableBach; // bachelor
      hasTOF = bachTrackCasc.hasTOF() ? 1 : 0;
      outputParts(outputCollision.lastIndex(),
                  casc.bachelorpt(),
                  casc.bacheloreta(),
                  casc.bachelorphi(),
                  aod::femtouniverseparticle::ParticleType::kCascadeBachelor,
                  0,
                  PIDBitmask(bachTrackCasc),
                  hasTOF,
                  childIDs,
                  0,
                  bachTrackCasc.sign());
      const int rowOfBachTrack = outputParts.lastIndex();
      if constexpr (isMC) {
        fillMCParticle(bachTrackCasc, o2::aod::femtouniverseparticle::ParticleType::kCascadeBachelor);
      }
      // cascade
      std::vector<int> indexCascChildID = {rowOfPosTrack, rowOfNegTrack, rowOfBachTrack};
      outputParts(outputCollision.lastIndex(),
                  casc.pt(),
                  casc.eta(),
                  casc.phi(),
                  aod::femtouniverseparticle::ParticleType::kCascade,
                  0,
                  PIDStrangeTOFBitmaskCasc(casc),
                  0,
                  indexCascChildID,
                  casc.mXi(),
                  casc.mOmega());
      outputCascParts(outputCollision.lastIndex(),
                      outputParts.lastIndex(),
                      casc.dcaV0daughters(),
                      casc.v0cosPA(col.posX(), col.posY(), col.posZ()),
                      casc.v0radius(),
                      casc.casccosPA(col.posX(), col.posY(), col.posZ()),
                      casc.dcacascdaughters(),
                      casc.cascradius(),
                      casc.dcapostopv(),
                      casc.dcanegtopv(),
                      casc.dcabachtopv(),
                      casc.dcav0topv(col.posX(), col.posY(), col.posZ()));
      if (confIsDebug) {
        fillDebugParticle<true, false, false>(posTrackCasc);  // QA for positive daughter
        fillDebugParticle<true, false, false>(negTrackCasc);  // QA for negative daughter
        fillDebugParticle<true, false, false>(bachTrackCasc); // QA for negative daughter
        fillDebugParticle<false, false, true>(casc);          // QA for cascade
      }
      if constexpr (isMC) {
        fillMCParticle(casc, o2::aod::femtouniverseparticle::ParticleType::kCascade);
      }
    }
  }

  template <bool isMC, typename HfCandidate, typename TrackType, typename CollisionType>
  void fillD0D0barData(CollisionType const&, TrackType const&, HfCandidate const& hfCands)
  {
    std::vector<int> childIDs = {0, 0}; // these IDs are necessary to keep track of the children
    std::vector<int> tmpIDtrack;        // this vector keeps track of the matching of the primary track table row <-> aod::track table global index
    double invMassD0 = 0.0;
    double invMassD0bar = 0.0;
    bool isD0D0bar = false;
    uint8_t daughFlag = 0; // flag = 0 (daugh of D0 or D0bar), 1 (daug of D0), -1 (daugh of D0bar)

    for (const auto& hfCand : hfCands) {

      if (!(hfCand.hfflag() & 1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        continue;
      }

      if (ConfD0Selection.useYCutD0Cand && std::abs(hfHelper.yD0(hfCand)) > ConfD0Selection.yD0CandMax) {
        continue;
      }

      if (!(ConfD0Selection.useYCutD0Cand) && std::abs(hfCand.eta()) > ConfD0Selection.trackD0CandEtaMax) {
        continue;
      }

      // int postrackID = hfCand.prong0().globalIndex();
      int postrackID = hfCand.prong0Id(); // Index to first prong
      int rowInPrimaryTrackTablePos = -1;
      rowInPrimaryTrackTablePos = getRowDaughters(postrackID, tmpIDtrack);
      childIDs[0] = rowInPrimaryTrackTablePos;
      childIDs[1] = 0;
      auto postrack = hfCand.template prong0_as<TrackType>();
      auto negtrack = hfCand.template prong1_as<TrackType>();

      if (hfCand.isSelD0() == 1 && hfCand.isSelD0bar() == 0) {
        invMassD0 = hfHelper.invMassD0ToPiK(hfCand);
        invMassD0bar = -hfHelper.invMassD0barToKPi(hfCand);
        isD0D0bar = true;
        daughFlag = 1;
      } else if (hfCand.isSelD0() == 0 && hfCand.isSelD0bar() == 1) {
        invMassD0 = -hfHelper.invMassD0ToPiK(hfCand);
        invMassD0bar = hfHelper.invMassD0barToKPi(hfCand);
        isD0D0bar = true;
        daughFlag = -1;
      } else if (hfCand.isSelD0() == 1 && hfCand.isSelD0bar() == 1) {
        invMassD0 = hfHelper.invMassD0ToPiK(hfCand);
        invMassD0bar = hfHelper.invMassD0barToKPi(hfCand);
        if (ConfD0Selection.storeD0D0barDoubleMassHypo) {
          isD0D0bar = true;
          daughFlag = 0;
        } else {
          isD0D0bar = false;
          daughFlag = 0;
        }
      } else {
        invMassD0 = 0.0;
        invMassD0bar = 0.0;
        isD0D0bar = false;
      }

      if (isD0D0bar) {
        outputParts(outputCollision.lastIndex(),
                    hfCand.ptProng0(),
                    RecoDecay::eta(std::array{hfCand.pxProng0(), hfCand.pyProng0(), hfCand.pzProng0()}), // eta
                    RecoDecay::phi(hfCand.pxProng0(), hfCand.pyProng0()),                                // phi
                    aod::femtouniverseparticle::ParticleType::kD0Child,
                    -999, // cutContainerV0.at(femto_universe_v0_selection::V0ContainerPosition::kPosCuts),
                    -999, // cutContainerV0.at(femto_universe_v0_selection::V0ContainerPosition::kPosPID),
                    -999,
                    childIDs,
                    postrack.sign(), // D0 mass -> positive daughter of D0/D0bar
                    daughFlag);      // D0bar mass -> sign that the daugh is from D0 or D0 decay
        const int rowOfPosTrack = outputParts.lastIndex();
        if constexpr (isMC) {
          fillMCParticle(postrack, o2::aod::femtouniverseparticle::ParticleType::kD0Child);
        }
        // int negtrackID = hfCand.prong1().globalIndex();
        int negtrackID = hfCand.prong1Id();
        int rowInPrimaryTrackTableNeg = -1;
        rowInPrimaryTrackTableNeg = getRowDaughters(negtrackID, tmpIDtrack);
        childIDs[0] = 0;
        childIDs[1] = rowInPrimaryTrackTableNeg;

        outputParts(outputCollision.lastIndex(),
                    hfCand.ptProng1(),
                    RecoDecay::eta(std::array{hfCand.pxProng1(), hfCand.pyProng1(), hfCand.pzProng1()}), // eta
                    RecoDecay::phi(hfCand.pxProng1(), hfCand.pyProng1()),                                // phi
                    aod::femtouniverseparticle::ParticleType::kD0Child,
                    -999, // cutContainerV0.at(femto_universe_v0_selection::V0ContainerPosition::kNegCuts),
                    -999, // cutContainerV0.at(femto_universe_v0_selection::V0ContainerPosition::kNegPID),
                    -999,
                    childIDs,
                    negtrack.sign(), // negative daughter of D0/D0bar
                    daughFlag);      // sign that the daugh is from D0 or D0 decay
        const int rowOfNegTrack = outputParts.lastIndex();
        if constexpr (isMC) {
          fillMCParticle(negtrack, o2::aod::femtouniverseparticle::ParticleType::kD0Child);
        }
        std::vector<int> indexChildID = {rowOfPosTrack, rowOfNegTrack};

        outputParts(outputCollision.lastIndex(),
                    hfCand.pt(),
                    hfCand.eta(),
                    hfCand.phi(),
                    aod::femtouniverseparticle::ParticleType::kD0,
                    -999, // cut, CutContainerType
                    -999, // PID, CutContainerType
                    -999.,
                    indexChildID,
                    invMassD0,     // D0 mass (mLambda)
                    invMassD0bar); // D0bar mass (mAntiLambda)

        if (confIsDebug) {
          fillDebugParticle<false, true, false>(postrack); // QA for positive daughter
          fillDebugParticle<false, true, false>(negtrack); // QA for negative daughter
          fillDebugParticle<false, true, false>(hfCand);   // QA for D0/D0bar
        }
        if constexpr (isMC) {
          fillMCParticle(hfCand, o2::aod::femtouniverseparticle::ParticleType::kD0);
        }
      }
    }
  }

  template <bool isMC, typename HfCandidate, typename TrackType, typename CollisionType>
  void fillD0D0barDataMl(CollisionType const&, TrackType const&, HfCandidate const& hfCands)
  {
    std::vector<int> childIDs = {0, 0}; // these IDs are necessary to keep track of the children
    std::vector<int> tmpIDtrack;        // this vector keeps track of the matching of the primary track table row <-> aod::track table global index
    double invMassD0 = 0.0;
    double invMassD0bar = 0.0;
    bool isD0D0bar = false;
    int8_t daughFlag = 0; // flag = 0 (daugh of D0 or D0bar), 1 (daug of D0), -1 (daugh of D0bar)

    for (const auto& hfCand : hfCands) {

      if (!(hfCand.hfflag() & 1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        continue;
      }

      if (ConfD0Selection.useYCutD0Cand && std::abs(hfHelper.yD0(hfCand)) > ConfD0Selection.yD0CandMax) {
        continue;
      }

      if (!(ConfD0Selection.useYCutD0Cand) && std::abs(hfCand.eta()) > ConfD0Selection.trackD0CandEtaMax) {
        continue;
      }

      int postrackID = hfCand.prong0Id(); // Index to first prong
      int rowInPrimaryTrackTablePos = -1;
      rowInPrimaryTrackTablePos = getRowDaughters(postrackID, tmpIDtrack);
      childIDs[0] = rowInPrimaryTrackTablePos;
      childIDs[1] = 0;
      auto postrack = hfCand.template prong0_as<TrackType>();
      auto negtrack = hfCand.template prong1_as<TrackType>();

      if (hfCand.isSelD0() == 1 && hfCand.isSelD0bar() == 0) {
        invMassD0 = hfHelper.invMassD0ToPiK(hfCand);
        invMassD0bar = -hfHelper.invMassD0barToKPi(hfCand);
        isD0D0bar = true;
        daughFlag = 1;
      } else if (hfCand.isSelD0() == 0 && hfCand.isSelD0bar() == 1) {
        invMassD0 = -hfHelper.invMassD0ToPiK(hfCand);
        invMassD0bar = hfHelper.invMassD0barToKPi(hfCand);
        isD0D0bar = true;
        daughFlag = -1;
      } else if (hfCand.isSelD0() == 1 && hfCand.isSelD0bar() == 1) {
        invMassD0 = hfHelper.invMassD0ToPiK(hfCand);
        invMassD0bar = hfHelper.invMassD0barToKPi(hfCand);
        if (ConfD0Selection.storeD0D0barDoubleMassHypo) {
          isD0D0bar = true;
          daughFlag = 0;
        } else {
          isD0D0bar = false;
          daughFlag = 0;
        }
      } else {
        invMassD0 = 0.0;
        invMassD0bar = 0.0;
        isD0D0bar = false;
      }

      if (isD0D0bar) {
        outputParts(outputCollision.lastIndex(),
                    hfCand.ptProng0(),
                    RecoDecay::eta(std::array{hfCand.pxProng0(), hfCand.pyProng0(), hfCand.pzProng0()}), // eta
                    RecoDecay::phi(hfCand.pxProng0(), hfCand.pyProng0()),                                // phi
                    aod::femtouniverseparticle::ParticleType::kD0Child,
                    -999, // cutContainerV0.at(femto_universe_v0_selection::V0ContainerPosition::kPosCuts),
                    -999, // cutContainerV0.at(femto_universe_v0_selection::V0ContainerPosition::kPosPID),
                    -999,
                    childIDs,
                    postrack.sign(), // D0 mass -> positive daughter of D0/D0bar
                    daughFlag);      // D0bar mass -> sign that the daugh is from D0 or D0 decay
        const int rowOfPosTrack = outputParts.lastIndex();
        if constexpr (isMC) {
          fillMCParticle(postrack, o2::aod::femtouniverseparticle::ParticleType::kD0Child);
        }
        // int negtrackID = hfCand.prong1().globalIndex();
        int negtrackID = hfCand.prong1Id();
        int rowInPrimaryTrackTableNeg = -1;
        rowInPrimaryTrackTableNeg = getRowDaughters(negtrackID, tmpIDtrack);
        childIDs[0] = 0;
        childIDs[1] = rowInPrimaryTrackTableNeg;

        outputParts(outputCollision.lastIndex(),
                    hfCand.ptProng1(),
                    RecoDecay::eta(std::array{hfCand.pxProng1(), hfCand.pyProng1(), hfCand.pzProng1()}), // eta
                    RecoDecay::phi(hfCand.pxProng1(), hfCand.pyProng1()),                                // phi
                    aod::femtouniverseparticle::ParticleType::kD0Child,
                    -999, // cutContainerV0.at(femto_universe_v0_selection::V0ContainerPosition::kNegCuts),
                    -999, // cutContainerV0.at(femto_universe_v0_selection::V0ContainerPosition::kNegPID),
                    -999,
                    childIDs,
                    negtrack.sign(), // negative daughter of D0/D0bar
                    daughFlag);      // sign that the daugh is from D0 or D0 decay
        const int rowOfNegTrack = outputParts.lastIndex();
        if constexpr (isMC) {
          fillMCParticle(negtrack, o2::aod::femtouniverseparticle::ParticleType::kD0Child);
        }
        std::vector<int> indexChildID = {rowOfPosTrack, rowOfNegTrack};

        outputParts(outputCollision.lastIndex(),
                    hfCand.pt(),
                    hfCand.eta(),
                    hfCand.phi(),
                    aod::femtouniverseparticle::ParticleType::kD0,
                    -999,  // cut, CutContainerType
                    -999,  // PID, CutContainerType
                    -999., // tempFitVar
                    indexChildID,
                    invMassD0,     // D0 mass (mLambda)
                    invMassD0bar); // D0bar mass (mAntiLambda)

        if (confIsDebug) {
          fillDebugParticle<false, true, false>(postrack); // QA for positive daughter
          fillDebugParticle<false, true, false>(negtrack); // QA for negative daughter
          if (hfCand.isSelD0() == 1) {
            fillDebugD0D0barML<true, false>(hfCand); // QA for D0/D0bar
          } else if (hfCand.isSelD0bar() == 1) {
            fillDebugD0D0barML<false, true>(hfCand);
          } else {
            fillDebugD0D0barML<false, false>(hfCand);
          }
        }
        if constexpr (isMC) {
          fillMCParticleD0(hfCand);
        }
      }
    }
  }

  template <bool isMC, typename McPart, typename HfCandidate, typename TrackType, typename CollisionType>
  void fillD0D0barMcMl(CollisionType const&, TrackType const&, HfCandidate const& hfCands, McPart const& mcParticles)
  {
    std::vector<int> childIDs = {0, 0}; // these IDs are necessary to keep track of the children
    std::vector<int> tmpIDtrack;        // this vector keeps track of the matching of the primary track table row <-> aod::track table global index
    double invMassD0 = 0.0;
    double invMassD0bar = 0.0;
    bool isD0D0bar = false;
    int indexMcRec = -1;
    int8_t sign = 0;
    int8_t daughFlag = 0; // flag = 0 (daugh of D0 or D0bar), 1 (daug of D0), -1 (daugh of D0bar)

    for (const auto& hfCand : hfCands) {

      if (!(hfCand.hfflag() & 1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        continue;
      }

      if (ConfD0Selection.useYCutD0Cand && std::abs(hfHelper.yD0(hfCand)) > ConfD0Selection.yD0CandMax) {
        continue;
      }

      if (!(ConfD0Selection.useYCutD0Cand) && std::abs(hfCand.eta()) > ConfD0Selection.trackD0CandEtaMax) {
        continue;
      }

      // Check whether the D0 candidate has the corresponding MC particle
      auto postrack = hfCand.template prong0_as<TrackType>();
      auto negtrack = hfCand.template prong1_as<TrackType>();
      auto arrayDaughters = std::array{postrack, negtrack};
      indexMcRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, Pdg::kD0, std::array{+kPiPlus, -kKPlus}, true, &sign);

      if (!(indexMcRec > -1)) {
        continue;
      }

      if (std::abs(hfCand.flagMcMatchRec()) == o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) {
        int postrackID = hfCand.prong0Id(); // Index to first prong
        int rowInPrimaryTrackTablePos = -1;
        rowInPrimaryTrackTablePos = getRowDaughters(postrackID, tmpIDtrack);
        childIDs[0] = rowInPrimaryTrackTablePos;
        childIDs[1] = 0;

        if (hfCand.isSelD0() == 1 && hfCand.isSelD0bar() == 0) {
          invMassD0 = hfHelper.invMassD0ToPiK(hfCand);
          invMassD0bar = -hfHelper.invMassD0barToKPi(hfCand);
          isD0D0bar = true;
          daughFlag = 1;
        } else if (hfCand.isSelD0() == 0 && hfCand.isSelD0bar() == 1) {
          invMassD0 = -hfHelper.invMassD0ToPiK(hfCand);
          invMassD0bar = hfHelper.invMassD0barToKPi(hfCand);
          isD0D0bar = true;
          daughFlag = -1;
        } else if (hfCand.isSelD0() == 1 && hfCand.isSelD0bar() == 1) {
          invMassD0 = hfHelper.invMassD0ToPiK(hfCand);
          invMassD0bar = hfHelper.invMassD0barToKPi(hfCand);
          if (ConfD0Selection.storeD0D0barDoubleMassHypo) {
            isD0D0bar = true;
            daughFlag = 0;
          } else {
            isD0D0bar = false;
            daughFlag = 0;
          }
        } else {
          invMassD0 = 0.0;
          invMassD0bar = 0.0;
          isD0D0bar = false;
        }

        if (isD0D0bar) {
          outputParts(outputCollision.lastIndex(),
                      hfCand.ptProng0(),
                      RecoDecay::eta(std::array{hfCand.pxProng0(), hfCand.pyProng0(), hfCand.pzProng0()}), // eta
                      RecoDecay::phi(hfCand.pxProng0(), hfCand.pyProng0()),                                // phi
                      aod::femtouniverseparticle::ParticleType::kD0Child,
                      -999, // cutContainerV0.at(femto_universe_v0_selection::V0ContainerPosition::kPosCuts),
                      -999, // cutContainerV0.at(femto_universe_v0_selection::V0ContainerPosition::kPosPID),
                      -999,
                      childIDs,
                      postrack.sign(), // D0 mass -> positive daughter of D0/D0bar
                      daughFlag);      // D0bar mass -> sign that the daugh is from D0 or D0 decay
          const int rowOfPosTrack = outputParts.lastIndex();
          if constexpr (isMC) {
            fillMCParticle(postrack, o2::aod::femtouniverseparticle::ParticleType::kD0Child);
          }
          // int negtrackID = hfCand.prong1().globalIndex();
          int negtrackID = hfCand.prong1Id();
          int rowInPrimaryTrackTableNeg = -1;
          rowInPrimaryTrackTableNeg = getRowDaughters(negtrackID, tmpIDtrack);
          childIDs[0] = 0;
          childIDs[1] = rowInPrimaryTrackTableNeg;

          outputParts(outputCollision.lastIndex(),
                      hfCand.ptProng1(),
                      RecoDecay::eta(std::array{hfCand.pxProng1(), hfCand.pyProng1(), hfCand.pzProng1()}), // eta
                      RecoDecay::phi(hfCand.pxProng1(), hfCand.pyProng1()),                                // phi
                      aod::femtouniverseparticle::ParticleType::kD0Child,
                      -999, // cutContainerV0.at(femto_universe_v0_selection::V0ContainerPosition::kNegCuts),
                      -999, // cutContainerV0.at(femto_universe_v0_selection::V0ContainerPosition::kNegPID),
                      -999,
                      childIDs,
                      negtrack.sign(), // negative daughter of D0/D0bar
                      daughFlag);      // sign that the daugh is from D0 or D0 decay
          const int rowOfNegTrack = outputParts.lastIndex();
          if constexpr (isMC) {
            fillMCParticle(negtrack, o2::aod::femtouniverseparticle::ParticleType::kD0Child);
          }
          std::vector<int> indexChildID = {rowOfPosTrack, rowOfNegTrack};

          outputParts(outputCollision.lastIndex(),
                      hfCand.pt(),
                      hfCand.eta(),
                      hfCand.phi(),
                      aod::femtouniverseparticle::ParticleType::kD0,
                      -999,  // cut, CutContainerType
                      -999,  // PID, CutContainerType
                      -999., // tempFitVar
                      indexChildID,
                      invMassD0,     // D0 mass (mLambda)
                      invMassD0bar); // D0bar mass (mAntiLambda)

          if (confIsDebug) {
            fillDebugParticle<false, true, false>(postrack); // QA for positive daughter
            fillDebugParticle<false, true, false>(negtrack); // QA for negative daughter
            if (hfCand.isSelD0() == 1) {
              fillDebugD0D0barMcMl<true, false>(hfCand); // QA for D0/D0bar
            } else if (hfCand.isSelD0bar() == 1) {
              fillDebugD0D0barMcMl<false, true>(hfCand);
            } else {
              fillDebugD0D0barMcMl<false, false>(hfCand);
            }
          }
          if constexpr (isMC) {
            auto particleMother = mcParticles.rawIteratorAt(indexMcRec);                        // gen. level pT
            auto yGen = RecoDecay::y(particleMother.pVector(), o2::constants::physics::MassD0); // gen. level y
            outputPartsMC(0, particleMother.pdgCode(), particleMother.pt(), yGen, particleMother.phi());
            outputPartsMCLabels(outputPartsMC.lastIndex());
          }
        }
      }
    }
  }

  template <bool isMC, typename TrackType, typename CollisionType>
  void fillPhi(CollisionType const& col, TrackType const& tracks)
  {
    std::vector<int> childIDs = {0, 0}; // these IDs are necessary to keep track of the children
    std::vector<int> tmpIDtrack;        // this vector keeps track of the matching of the primary track table row <-> aod::track table global index
    // lorentz vectors and filling the tables
    for (const auto& [p1, p2] : combinations(soa::CombinationsFullIndexPolicy(tracks, tracks))) {
      if (!trackCuts.isSelectedMinimal(p1) || !trackCuts.isSelectedMinimal(p1)) {
        continue;
      }
      if ((!(p1.sign() == 1)) || (!(p2.sign() == -1))) {
        continue;
      }
      // implementing PID cuts for phi children
      if (ConfPhiSelection.confPhiDoLFPID4Kaons) {
        if ((p1.isGlobalTrackWoDCA() == false) || (p2.isGlobalTrackWoDCA() == false) || (p1.isPVContributor() == false) || (p2.isPVContributor() == false)) {
          continue;
        }
        if (!(isKaonNSigmaLF(p1.pt(), trackCuts.getNsigmaTPC(p1, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p1, o2::track::PID::Kaon), p1.hasTOF()))) {
          continue;
        }
        if (!(isKaonNSigmaLF(p2.pt(), trackCuts.getNsigmaTPC(p2, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p2, o2::track::PID::Kaon), p2.hasTOF()))) {
          continue;
        }
      } else {
        if (!(isKaonNSigma(p1.pt(), trackCuts.getNsigmaTPC(p1, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p1, o2::track::PID::Kaon)))) {
          continue;
        }
        if (!(isKaonNSigma(p2.pt(), trackCuts.getNsigmaTPC(p2, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p2, o2::track::PID::Kaon)))) {
          continue;
        }
      }
      if (isKaonRejected(p1.p(), trackCuts.getNsigmaTPC(p1, o2::track::PID::Proton), trackCuts.getNsigmaTOF(p1, o2::track::PID::Proton), trackCuts.getNsigmaTPC(p1, o2::track::PID::Pion), trackCuts.getNsigmaTOF(p1, o2::track::PID::Pion))) {
        continue;
      }
      if (isKaonRejected(p2.p(), trackCuts.getNsigmaTPC(p2, o2::track::PID::Proton), trackCuts.getNsigmaTOF(p2, o2::track::PID::Proton), trackCuts.getNsigmaTPC(p2, o2::track::PID::Pion), trackCuts.getNsigmaTOF(p2, o2::track::PID::Pion))) {
        continue;
      }

      TLorentzVector part1Vec;
      TLorentzVector part2Vec;

      const auto mMassOne = o2::constants::physics::MassKPlus;  // FIXME: Get from the PDG service of the common header
      const auto mMassTwo = o2::constants::physics::MassKMinus; // FIXME: Get from the PDG service of the common header

      part1Vec.SetPtEtaPhiM(p1.pt(), p1.eta(), p1.phi(), mMassOne);
      part2Vec.SetPtEtaPhiM(p2.pt(), p2.eta(), p2.phi(), mMassTwo);

      TLorentzVector sumVec(part1Vec);
      sumVec += part2Vec;

      float phiEta = sumVec.Eta();
      if (std::abs(phiEta) > ConfPhiSelection.confPhiEtaHighLimit.value) {
        continue;
      }

      float phiPt = sumVec.Pt();
      if ((phiPt < ConfPhiSelection.confPhiPtLowLimit.value) || (phiPt > ConfPhiSelection.confPhiPtHighLimit.value)) {
        continue;
      }

      float phiPhi = RecoDecay::constrainAngle(sumVec.Phi(), 0);
      float phiM = sumVec.M();

      if ((phiM < ConfPhiSelection.confPhiInvMassLowLimit) || (phiM > ConfPhiSelection.confPhiInvMassUpLimit))
        continue;

      phiCuts.fillQA<aod::femtouniverseparticle::ParticleType::kPhi, aod::femtouniverseparticle::ParticleType::kPhiChild>(col, p1, p1, p2, 321, -321); ///\todo fill QA also for daughters

      int postrackID = p1.globalIndex();
      int rowInPrimaryTrackTablePos = -1; // does it do anything?
      rowInPrimaryTrackTablePos = getRowDaughters(postrackID, tmpIDtrack);
      childIDs[0] = rowInPrimaryTrackTablePos;
      childIDs[1] = 0;

      outputParts(outputCollision.lastIndex(), p1.pt(),
                  p1.eta(), p1.phi(),
                  aod::femtouniverseparticle::ParticleType::kPhiChild,
                  -999,       // cutContainer
                  -999,       // cutContainer
                  p1.dcaXY(), // tempFitVar
                  childIDs,
                  0,
                  p1.sign());
      const int rowOfPosTrack = outputParts.lastIndex();
      if constexpr (isMC) {
        fillMCParticle(p1, o2::aod::femtouniverseparticle::ParticleType::kPhiChild);
      }
      int negtrackID = p2.globalIndex();
      int rowInPrimaryTrackTableNeg = -1;
      rowInPrimaryTrackTableNeg = getRowDaughters(negtrackID, tmpIDtrack);
      childIDs[0] = 0;
      childIDs[1] = rowInPrimaryTrackTableNeg;
      outputParts(outputCollision.lastIndex(),
                  p2.pt(),
                  p2.eta(),
                  p2.phi(),
                  aod::femtouniverseparticle::ParticleType::kPhiChild,
                  -999,       // cutContainer
                  -999,       // cutContainer
                  p2.dcaXY(), // tempFitVar
                  childIDs,
                  0,
                  p2.sign()); // sign, workaround for now
      const int rowOfNegTrack = outputParts.lastIndex();
      if constexpr (isMC) {
        fillMCParticle(p2, o2::aod::femtouniverseparticle::ParticleType::kPhiChild);
      }
      std::vector<int> indexChildID = {rowOfPosTrack, rowOfNegTrack};

      outputParts(outputCollision.lastIndex(),
                  phiPt,
                  phiEta,
                  phiPhi,
                  aod::femtouniverseparticle::ParticleType::kPhi,
                  -999, // cutContainer
                  0,
                  phiM, // tempFitVar
                  indexChildID,
                  phiM, // phi.mLambda(), //for now it will have a mLambda getter, maybe we will change it in the future so it's more logical
                  0);   // phi.mAntiLambda() <- sign

      if (confIsDebug) {
        fillDebugParticle<true, false, false>(p1); // QA for positive daughter
        fillDebugParticle<true, false, false>(p2); // QA for negative daughter
        fillDebugParticle<false, true, false>(p1); // QA for phi
      }
      if constexpr (isMC) {
        fillMCParticlePhi(p1, p2);
      }
    }
  }

  template <typename TrackType, bool transientLabels = false, bool resolveDaughs = false>
  void fillParticles(TrackType const& tracks, std::optional<std::reference_wrapper<const std::set<int>>> recoMcIds = std::nullopt)
  {
    std::vector<int> childIDs = {0, 0}; // these IDs are necessary to keep track of the children
    std::vector<int> tmpIDtrack;

    for (const auto& particle : tracks) {
      /// if the most open selection criteria are not fulfilled there is no
      /// point looking further at the track

      if (particle.eta() < -ConfFilterCuts.confEtaFilterCut || particle.eta() > ConfFilterCuts.confEtaFilterCut)
        continue;
      if (particle.pt() < ConfFilterCuts.confPtLowFilterCut || particle.pt() > ConfFilterCuts.confPtHighFilterCut)
        continue;

      uint32_t pdgCode = static_cast<uint32_t>(particle.pdgCode());

      if (ConfGeneral.confMCTruthAnalysisWithPID) {
        bool pass = false;
        std::vector<int> tmpPDGCodes = ConfGeneral.confMCTruthPDGCodes; // necessary due to some features of the Configurable
        for (auto const& pdg : tmpPDGCodes) {
          if (static_cast<int>(pdg) == static_cast<int>(pdgCode)) {
            if (pdgCode == Pdg::kPhi) { // && (recoMcIds && recoMcIds->get().contains(particle.globalIndex()))) { // ATTENTION: all Phi mesons are NOT primary particles
              pass = true;
            } else if (pdgCode == Pdg::kD0) {
              pass = true;
            } else {
              if (confStoreMCmothers || particle.isPhysicalPrimary() || (ConfGeneral.confActivateSecondaries && recoMcIds && recoMcIds->get().contains(particle.globalIndex())))
                pass = true;
            }
          }
        }
        if (!pass)
          continue;
      }

      // we cannot use isSelectedMinimal since it takes Ncls
      // if (!trackCuts.isSelectedMinimal(track)) {
      //   continue;
      // }

      // trackCuts.fillQA<aod::femtouniverseparticle::ParticleType::kTrack,
      //                  aod::femtouniverseparticle::TrackType::kNoChild>(track);
      //  the bit-wise container of the systematic variations is obtained
      // auto cutContainer = trackCuts.getCutContainer<aod::femtouniverseparticle::CutContainerType>(track);
      // instead of the bitmask, the PDG of the particle is stored as uint32_t

      int32_t variablePDG = confStoreMCmothers ? getMotherPDG(particle) : particle.pdgCode();

      // now the table is filled
      if constexpr (resolveDaughs) {
        tmpIDtrack.push_back(particle.globalIndex());
        continue;
      }

      if (ConfGeneral.confIsActivateCascade)
        childIDs.push_back(0);
      outputParts(outputCollision.lastIndex(),
                  particle.pt(),
                  particle.eta(),
                  particle.phi(),
                  aod::femtouniverseparticle::ParticleType::kMCTruthTrack,
                  0,
                  pdgCode,
                  variablePDG,
                  childIDs,
                  0,
                  0);

      if (confIsDebug) {
        fillDebugParticle<false, true, false>(particle);
      }

      // Workaround to keep the FDParticles and MC label tables
      // aligned, so that they can be joined in the task.
      if constexpr (transientLabels) {
        outputPartsMCLabels(-1);
        outputDebugPartsMC(9999);
      }
    }
    if constexpr (resolveDaughs) {
      childIDs[0] = 0;
      childIDs[1] = 0;
      auto minDaughs = 2ul;
      if (ConfGeneral.confIsActivateCascade) {
        childIDs.push_back(0);
        minDaughs = 3ul;
      }
      for (std::size_t i = 0; i < tmpIDtrack.size(); i++) {
        const auto& particle = tracks.iteratorAt(tmpIDtrack[i] - tracks.begin().globalIndex());
        for (int daughIndex = 0, n = std::min(minDaughs, particle.daughtersIds().size()); daughIndex < n; daughIndex++) {
          // loop to find the corresponding index of the daughters
          for (std::size_t j = 0; j < tmpIDtrack.size(); j++) {
            if (tmpIDtrack[j] == particle.daughtersIds()[daughIndex]) {
              childIDs[daughIndex] = i - j;
              break;
            }
          }
        }

        int32_t variablePDG = confStoreMCmothers ? getMotherPDG(particle) : particle.pdgCode();

        outputParts(outputCollision.lastIndex(),
                    particle.pt(),
                    particle.eta(),
                    particle.phi(),
                    aod::femtouniverseparticle::ParticleType::kMCTruthTrack,
                    0,
                    static_cast<uint32_t>(particle.pdgCode()),
                    variablePDG,
                    childIDs,
                    0,
                    0);

        if (confIsDebug) {
          fillDebugParticle<false, true, false>(particle);
        }

        // Workaround to keep the FDParticles and MC label tables
        // aligned, so that they can be joined in the task.
        if constexpr (transientLabels) {
          outputPartsMCLabels(-1);
          outputDebugPartsMC(9999);
        }
      }
    }
  }

  template <typename TrackType, bool transientLabels = false>
  void fillMCTruthParticlesD0(TrackType const& mcParts)
  {
    std::vector<int> childIDs = {0, 0}; // these IDs are necessary to keep track of the children
    std::vector<int> tmpIDtrack;

    for (const auto& particle : mcParts) {

      if (particle.pt() < ConfD0Selection.trackD0pTGenMin || particle.pt() > ConfD0Selection.trackD0pTGenMax)
        continue;

      int pdgCode = particle.pdgCode();

      if (ConfGeneral.confMCTruthAnalysisWithPID) {
        bool pass = false;
        std::vector<int> tmpPDGCodes = ConfGeneral.confMCTruthPDGCodes; // necessary due to some features of the Configurable
        for (auto const& pdg : tmpPDGCodes) {
          if (static_cast<int>(pdg) == static_cast<int>(pdgCode)) {
            if (pdgCode == Pdg::kPhi) {
              pass = true;
            } else if (std::abs(pdgCode) == Pdg::kD0) {
              pass = true;
            } else {
              if (particle.isPhysicalPrimary())
                pass = true;
            }
          }
        }
        if (!pass)
          continue;
      }

      /// check if we end-up with the correct final state using MC info
      int8_t sign = 0;
      int8_t origin = -99;
      int8_t mcGenFlag = -99;
      if (std::abs(particle.pdgCode()) == Pdg::kD0 && !RecoDecay::isMatchedMCGen(mcParts, particle, Pdg::kD0, std::array{+kPiPlus, -kKPlus}, true, &sign)) {
        /// check if we have D0(bar) → π± K∓
        continue;
      }

      if (std::abs(particle.pdgCode()) == Pdg::kD0) {
        if (std::abs(particle.y()) > ConfD0Selection.yD0CandGenMax) {
          continue;
        } else {
          origin = RecoDecay::getCharmHadronOrigin(mcParts, particle);
          mcGenFlag = particle.flagMcMatchGen();
        }
      } else {
        if (std::abs(particle.eta()) > ConfD0Selection.trackD0CandEtaMax) {
          continue;
        }
      }

      outputParts(outputCollision.lastIndex(),
                  particle.pt(),
                  particle.eta(),
                  particle.phi(),
                  aod::femtouniverseparticle::ParticleType::kMCTruthTrack,
                  -999.,
                  pdgCode,
                  pdgCode, // getter tempFitVar
                  childIDs,
                  mcGenFlag,
                  origin); // D0(bar) origin

      if (confIsDebug) {
        fillDebugParticle<false, true, false>(particle);
      }

      // Workaround to keep the FDParticles and MC label tables
      // aligned, so that they can be joined in the task.
      if constexpr (transientLabels) {
        outputPartsMCLabels(-1);
        outputDebugPartsMC(-999);
      }
    }
  }

  template <bool isMC, typename V0Type, typename TrackType,
            typename CollisionType>
  void fillCollisionsAndTracksAndV0AndPhi(CollisionType const& col, TrackType const& tracks, V0Type const& fullV0s)
  {
    const auto colcheck = fillCollisions<isMC>(col, tracks);
    if (colcheck) {
      fillTracks<isMC>(tracks);
      if (ConfGeneral.confIsActivateV0) {
        fillV0<isMC>(col, fullV0s, tracks);
      }
      if (ConfGeneral.confIsActivatePhi) {
        fillPhi<isMC>(col, tracks);
      }
    }
  }

  void processFullData(aod::FemtoFullCollision const& col,
                       aod::BCsWithTimestamps const&,
                       aod::FemtoFullTracks const& tracks,
                       o2::aod::V0Datas const& fullV0s)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    fillCollisionsAndTracksAndV0AndPhi<false>(col, tracks, fullV0s);
  }
  PROCESS_SWITCH(FemtoUniverseProducerTask, processFullData, "Provide experimental data", false);

  void processTrackV0(aod::FemtoFullCollision const& col,
                      aod::BCsWithTimestamps const&,
                      soa::Filtered<aod::FemtoFullTracks> const& tracks,
                      o2::aod::V0Datas const& fullV0s)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    fillCollisionsAndTracksAndV0AndPhi<false>(col, tracks, fullV0s);
  }
  PROCESS_SWITCH(FemtoUniverseProducerTask, processTrackV0, "Provide experimental data for track v0", false);

  void processTrackCascadeData(aod::FemtoFullCollision const& col,
                               aod::BCsWithTimestamps const&,
                               soa::Filtered<aod::FemtoFullTracks> const& tracks,
                               o2::aod::CascDatas const& fullCascades)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    const auto colcheck = fillCollisions<false>(col, tracks);
    if (colcheck) {
      fillTracks<false>(tracks);
      fillCascade<false>(col, fullCascades, tracks);
    }
  }
  PROCESS_SWITCH(FemtoUniverseProducerTask, processTrackCascadeData, "Provide experimental data for track cascades", false);

  void processTrackV0Cascade(aod::FemtoFullCollision const& col,
                             aod::BCsWithTimestamps const&,
                             soa::Filtered<aod::FemtoFullTracks> const& tracks,
                             soa::Join<o2::aod::V0Datas, o2::aod::V0TOFNSigmas> const& fullV0s,
                             soa::Join<o2::aod::CascDatas, o2::aod::CascTOFNSigmas> const& fullCascades)
  {
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    const auto colcheck = fillCollisions<false>(col, tracks);
    if (colcheck) {
      fillTracks<false>(tracks);
      if (ConfGeneral.confIsActivateV0) {
        fillV0<false>(col, fullV0s, tracks);
      }
      if (ConfGeneral.confIsActivateCascade) {
        fillCascade<false>(col, fullCascades, tracks);
      }
    }
  }
  PROCESS_SWITCH(FemtoUniverseProducerTask, processTrackV0Cascade, "Provide experimental data for track, v0 and cascades", false);

  /*void processTrackV0CentRun3(aod::FemtoFullCollisionCentRun3 const& col,
                              aod::BCsWithTimestamps const&,
                              soa::Filtered<aod::FemtoFullTracks> const& tracks,
                              o2::aod::V0Datas const& fullV0s)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    fillCollisionsCentRun3<false>(col, tracks, fullV0s);
    fillTracks<false>(tracks);
    fillV0<false>(col, fullV0s, tracks);
  }
  PROCESS_SWITCH(FemtoUniverseProducerTask, processTrackV0CentRun3, "Provide experimental data for track v0", false);*/

  void processFullMC(aod::FemtoFullCollisionMC const& col,
                     aod::BCsWithTimestamps const&,
                     soa::Join<aod::FemtoFullTracks, aod::McTrackLabels> const& tracks,
                     aod::McCollisions const&,
                     aod::McParticles const&,
                     soa::Join<o2::aod::V0Datas, aod::McV0Labels> const& fullV0s)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    fillCollisionsAndTracksAndV0AndPhi<true>(col, tracks, fullV0s);
  }
  PROCESS_SWITCH(FemtoUniverseProducerTask, processFullMC, "Provide MC data (tracks, V0, Phi)", false);

  void processTrackMC(aod::FemtoFullCollisionMC const& col,
                      aod::BCsWithTimestamps const&,
                      soa::Join<aod::FemtoFullTracks, aod::McTrackLabels> const& tracks,
                      aod::McCollisions const&,
                      aod::McParticles const&)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    const auto colcheck = fillCollisions<true>(col, tracks);
    if (colcheck) {
      fillTracks<true>(tracks);
    }
  }
  PROCESS_SWITCH(FemtoUniverseProducerTask, processTrackMC, "Provide MC data for track analysis", false);

  void processTrackPhiMC(aod::FemtoFullCollisionMC const& col,
                         aod::BCsWithTimestamps const&,
                         soa::Join<aod::FemtoFullTracks, aod::McTrackLabels> const& tracks,
                         aod::McCollisions const&,
                         aod::McParticles const&)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    const auto colcheck = fillCollisions<true>(col, tracks);
    if (colcheck) {
      fillTracks<true>(tracks);
      fillPhi<true>(col, tracks);
    }
  }
  PROCESS_SWITCH(FemtoUniverseProducerTask, processTrackPhiMC, "Provide MC data for track Phi analysis", false);

  void processTrackData(aod::FemtoFullCollision const& col,
                        aod::BCsWithTimestamps const&,
                        aod::FemtoFullTracks const& tracks)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    const double ir = 1.0; // fetch IR
    // fill the tables
    const auto colcheck = fillCollisions<false>(col, tracks);
    if (colcheck) {
      if (confFillCollExt) {
        fillCollisionsCentRun3ColExtra<false>(col, ir);
      }
      fillTracks<false>(tracks);
    }
  }
  PROCESS_SWITCH(FemtoUniverseProducerTask, processTrackData,
                 "Provide experimental data for track track", true);

  void processTrackDataCentPP(aod::FemtoFullCollisionCentPP const& col,
                              aod::BCsWithTimestamps const&,
                              aod::FemtoFullTracks const& tracks)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    const double ir = 1.0; // fetch IR
    // fill the tables
    const auto colcheck = fillCollisionsCentPP<false>(col, tracks);
    if (colcheck) {
      if (confFillCollExt) {
        fillCollisionsCentRun3ColExtra<false>(col, ir);
      }
      fillTracks<false>(tracks);
    }
  }
  PROCESS_SWITCH(FemtoUniverseProducerTask, processTrackDataCentPP,
                 "Provide experimental data for track track", false);

  // using FilteredFemtoFullTracks = soa::Filtered<FemtoFullTracks>;
  void processTrackPhiData(aod::FemtoFullCollision const& col,
                           aod::BCsWithTimestamps const&,
                           soa::Filtered<aod::FemtoFullTracks> const& tracks)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    const auto colcheck = fillCollisions<false>(col, tracks);
    if (colcheck) {
      fillTracks<false>(tracks);
      fillPhi<false>(col, tracks);
    }
  }
  PROCESS_SWITCH(FemtoUniverseProducerTask, processTrackPhiData,
                 "Provide experimental data for track phi", false);

  void processTrackD0mesonData(aod::FemtoFullCollision const& col,
                               aod::BCsWithTimestamps const&,
                               soa::Filtered<aod::FemtoFullTracks> const& tracks,
                               soa::Join<aod::HfCand2Prong, aod::HfSelD0> const& candidates)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    const auto colcheck = fillCollisions<false>(col, tracks);
    if (colcheck) {
      fillTracks<false>(tracks);
      fillD0D0barData<false>(col, tracks, candidates);
    }
  }
  PROCESS_SWITCH(FemtoUniverseProducerTask, processTrackD0mesonData,
                 "Provide experimental data for track D0 meson", false);

  void processTrackD0DataML(aod::FemtoFullCollision const& col,
                            aod::BCsWithTimestamps const&,
                            soa::Filtered<aod::FemtoFullTracks> const& tracks,
                            soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfMlD0> const& candidates)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    const auto colcheck = fillCollisions<false>(col, tracks);
    if (colcheck) {
      fillTracks<false>(tracks);
      fillD0D0barDataMl<false>(col, tracks, candidates);
    }
  }
  PROCESS_SWITCH(FemtoUniverseProducerTask, processTrackD0DataML,
                 "Provide experimental data for track D0 meson using ML selection for D0s", false);

  void processTrackMCTruth(aod::McCollision const&,
                           soa::SmallGroups<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::McCollisionLabels>> const& collisions,
                           aod::McParticles const& mcParticles,
                           aod::BCsWithTimestamps const&)
  {
    // magnetic field for run not needed for mc truth
    // fill the tables
    fillMCTruthCollisions(collisions, mcParticles);
    fillParticles(mcParticles);
  }
  PROCESS_SWITCH(FemtoUniverseProducerTask, processTrackMCTruth, "Provide MC data for MC truth track analysis", false);

  void processTrackMCGen(aod::McCollision const& col,
                         aod::McParticles const& mcParticles)
  {
    outputCollision(col.posZ(), 0, 0, 2, 0);
    fillParticles(mcParticles);
  }
  PROCESS_SWITCH(FemtoUniverseProducerTask, processTrackMCGen, "Provide MC Generated for model comparisons", false);

  template <class T>
  using HasBachelor = decltype(std::declval<T&>().bachelorphi());

  Preslice<aod::McParticles> perMCCollision = aod::mcparticle::mcCollisionId;
  PresliceUnsorted<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::McCollisionLabels>> recoCollsPerMCColl = aod::mcparticle::mcCollisionId;
  PresliceUnsorted<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::McCollisionLabels>> recoCollsPerMCCollCentPbPb = aod::mcparticle::mcCollisionId;
  Preslice<soa::Join<aod::FemtoFullTracks, aod::McTrackLabels>> perCollisionTracks = aod::track::collisionId;

  template <class StrangePartType>

  void processTruthAndFullMC(
    aod::McCollisions const& mccols,
    aod::McParticles const& mcParticles,
    soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::McCollisionLabels> const& collisions,
    soa::Filtered<soa::Join<aod::FemtoFullTracks, aod::McTrackLabels>> const& tracks,
    StrangePartType const& strangeParts,
    aod::BCsWithTimestamps const&,
    Preslice<StrangePartType>& ps)
  {
    // recos
    static std::set<int> recoMcIds;
    static std::set<int> mcColIds;
    recoMcIds.clear();
    mcColIds.clear();

    for (const auto& col : collisions) {
      auto groupedTracks = tracks.sliceBy(perCollisionTracks, col.globalIndex());
      auto groupedStrageParts = strangeParts.sliceBy(ps, col.globalIndex());
      getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
      if constexpr (std::experimental::is_detected<HasBachelor, typename StrangePartType::iterator>::value) {
        const auto colcheck = fillCollisions<true>(col, groupedTracks);
        if (colcheck) {
          mcColIds.insert(col.mcCollisionId());
          fillTracks<true>(groupedTracks);
          fillCascade<true>(col, groupedStrageParts, groupedTracks);
        }
      } else {
        mcColIds.insert(col.mcCollisionId());
        fillCollisionsAndTracksAndV0AndPhi<true>(col, groupedTracks, groupedStrageParts);
      }
      for (const auto& track : groupedTracks) {
        if (trackCuts.isSelectedMinimal(track))
          recoMcIds.insert(track.mcParticleId());
      }
    }

    // truth
    for (const auto& mccol : mccols) {
      auto groupedMCParticles = mcParticles.sliceBy(perMCCollision, mccol.globalIndex());
      if (confCollMCTruthOnlyReco && !mcColIds.contains(mccol.globalIndex())) {
        continue;
      }
      auto groupedCollisions = collisions.sliceBy(recoCollsPerMCColl, mccol.globalIndex());
      fillMCTruthCollisions(groupedCollisions, groupedMCParticles);                           // fills the reco collisions for mc collision
      fillParticles<decltype(groupedMCParticles), true, true>(groupedMCParticles, recoMcIds); // fills mc particles
    }
  }

  void processTruthAndFullMCCentRun3(aod::McCollisions const& mccols,
                                     aod::McParticles const& mcParticles,
                                     aod::FemtoFullCollisionCentRun3MCs const& collisions,
                                     soa::Filtered<soa::Join<aod::FemtoFullTracks, aod::McTrackLabels>> const& tracks,
                                     aod::BCsWithTimestamps const&)
  {
    // recos
    std::set<int> recoMcIds;
    std::set<int> mcCollisions;
    for (const auto& col : collisions) {
      auto groupedTracks = tracks.sliceBy(perCollisionTracks, col.globalIndex());
      auto bc = col.bc_as<aod::BCsWithTimestamps>();
      getMagneticFieldTesla(bc);
      const auto ir = mRateFetcher.fetch(ccdb.service, bc.timestamp(), mRunNumber, "ZNC hadronic") * 1.e-3; // fetch IR
      bool colcheck = false;
      // fill the tables

      if (ConfGeneral.confIsCent) {
        colcheck = fillCollisionsCentRun3<true>(col);
      } else {
        colcheck = fillCollisions<true>(col, groupedTracks);
      }

      if (colcheck) {
        mcCollisions.insert(col.mcCollisionId());
        fillCollisionsCentRun3ColExtra<true>(col, ir);
        fillTracks<true>(groupedTracks);
      }
      for (const auto& track : groupedTracks) {
        if (trackCuts.isSelectedMinimal(track))
          recoMcIds.insert(track.mcParticleId());
      }
    }

    // truth
    for (const auto& mccol : mccols) {
      if (confCollMCTruthOnlyReco && !mcCollisions.contains(mccol.globalIndex())) {
        continue;
      }
      auto groupedCollisions = collisions.sliceBy(recoCollsPerMCCollCentPbPb, mccol.globalIndex());
      for (const auto& col : groupedCollisions) {
        const auto colcheck = fillMCTruthCollisionsCentRun3(col); // fills the reco collisions for mc collision
        if (colcheck) {
          auto groupedMCParticles = mcParticles.sliceBy(perMCCollision, mccol.globalIndex());
          outputCollExtra(1.0, 1.0);
          if (!ConfTrkSelection.confIsOnlyMCTrack) {
            fillParticles<decltype(groupedMCParticles), true, true>(groupedMCParticles, recoMcIds); // fills mc particles
          } else {
            fillTracksMCTruth(groupedMCParticles);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(FemtoUniverseProducerTask, processTruthAndFullMCCentRun3, "Provide both MC truth and reco for tracks in Pb-Pb", false);

  Preslice<soa::Join<o2::aod::V0Datas, aod::V0TOFNSigmas, aod::McV0Labels>> perCollisionV0s = aod::track::collisionId;
  void processTruthAndFullMCV0(
    aod::McCollisions const& mccols,
    aod::McParticles const& mcParticles,
    soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::McCollisionLabels> const& collisions,
    soa::Filtered<soa::Join<aod::FemtoFullTracks, aod::McTrackLabels>> const& tracks,
    soa::Join<o2::aod::V0Datas, aod::V0TOFNSigmas, aod::McV0Labels> const& fullV0s,
    aod::BCsWithTimestamps const& bcs)
  {
    processTruthAndFullMC(mccols, mcParticles, collisions, tracks, fullV0s, bcs, perCollisionV0s);
  }
  PROCESS_SWITCH(FemtoUniverseProducerTask, processTruthAndFullMCV0, "Provide both MC truth and reco for tracks and V0s", false);

  Preslice<soa::Join<o2::aod::CascDatas, aod::CascTOFNSigmas, aod::McCascLabels>> perCollisionCascs = aod::track::collisionId;
  void processTruthAndFullMCCasc(
    aod::McCollisions const& mccols,
    aod::McParticles const& mcParticles,
    soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::McCollisionLabels> const& collisions,
    soa::Filtered<soa::Join<aod::FemtoFullTracks, aod::McTrackLabels>> const& tracks,
    soa::Join<o2::aod::CascDatas, aod::CascTOFNSigmas, aod::McCascLabels> const& fullCascades,
    aod::BCsWithTimestamps const& bcs)
  {
    processTruthAndFullMC(mccols, mcParticles, collisions, tracks, fullCascades, bcs, perCollisionCascs);
  }
  PROCESS_SWITCH(FemtoUniverseProducerTask, processTruthAndFullMCCasc, "Provide both MC truth and reco for tracks and Cascades", false);

  void processTruthAndFullMCCentRun3Casc(
    aod::McCollisions const& mccols,
    aod::McParticles const& mcParticles,
    aod::FemtoFullCollisionCentRun3MCs const& collisions,
    soa::Filtered<soa::Join<aod::FemtoFullTracks, aod::McTrackLabels>> const& tracks,
    soa::Join<o2::aod::CascDatas, aod::McCascLabels> const& fullCascades,
    aod::BCsWithTimestamps const&)
  {

    // recos
    std::set<int> recoMcIds;
    for (const auto& col : collisions) {
      auto groupedTracks = tracks.sliceBy(perCollisionTracks, col.globalIndex());
      auto groupedCascParts = fullCascades.sliceBy(perCollisionCascs, col.globalIndex());
      getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
      const auto colcheck = fillCollisionsCentRun3<true>(col);
      if (colcheck) {
        fillTracks<true>(groupedTracks);
        fillCascade<true>(col, groupedCascParts, groupedTracks);
      }
      for (const auto& track : groupedTracks) {
        if (trackCuts.isSelectedMinimal(track))
          recoMcIds.insert(track.mcParticleId());
      }
    }

    // truth
    for (const auto& mccol : mccols) {
      auto groupedCollisions = collisions.sliceBy(recoCollsPerMCCollCentPbPb, mccol.globalIndex());
      for (const auto& col : groupedCollisions) {
        const auto colcheck = fillMCTruthCollisionsCentRun3(col); // fills the reco collisions for mc collision
        if (colcheck) {
          auto groupedMCParticles = mcParticles.sliceBy(perMCCollision, mccol.globalIndex());
          outputCollExtra(1.0, 1.0);
          fillParticles<decltype(groupedMCParticles), true, true>(groupedMCParticles, recoMcIds); // fills mc particles
        }
      }
    }
  }
  PROCESS_SWITCH(FemtoUniverseProducerTask, processTruthAndFullMCCentRun3Casc, "Provide both MC truth and reco for tracks and cascades with centrality", false);

  Preslice<soa::Join<aod::McParticles, aod::HfCand2ProngMcGen>> mcPartPerMcColl = aod::mcparticle::mcCollisionId;
  Preslice<soa::Join<aod::HfCand2Prong, aod::HfCand2ProngMcRec, aod::HfSelD0, aod::HfMlD0>> perCollisionD0s = aod::track::collisionId;
  void processTrackD0MC(aod::McCollisions const& mccols,
                        soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::McCollisionLabels> const& collisions,
                        soa::Filtered<soa::Join<aod::FemtoFullTracks, aod::McTrackLabels>> const& tracks,
                        soa::Join<aod::McParticles, aod::HfCand2ProngMcGen> const& hfMcGenCands,
                        soa::Join<aod::HfCand2Prong, aod::HfCand2ProngMcRec, aod::HfSelD0, aod::HfMlD0> const& hfMcRecoCands,
                        aod::BCsWithTimestamps const&,
                        aod::McParticles const& mcParts)
  {
    // MC Reco
    for (const auto& col : collisions) {
      auto groupedTracks = tracks.sliceBy(perCollisionTracks, col.globalIndex());
      auto groupedD0s = hfMcRecoCands.sliceBy(perCollisionD0s, col.globalIndex());
      // get magnetic field for run
      getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
      // fill the tables
      const auto colcheck = fillCollisions<true>(col, tracks);
      if (colcheck) {
        fillTracks<true>(groupedTracks);
        fillD0D0barMcMl<true>(col, groupedTracks, groupedD0s, mcParts);
      }
    }
    // MC Truth
    for (const auto& mccol : mccols) {
      auto groupedMCParticles = hfMcGenCands.sliceBy(mcPartPerMcColl, mccol.globalIndex());
      auto groupedCollisions = collisions.sliceBy(recoCollsPerMCColl, mccol.globalIndex());
      fillMCTruthCollisions(groupedCollisions, groupedMCParticles);                   // fills the reco collisions for mc collision
      fillMCTruthParticlesD0<decltype(groupedMCParticles), true>(groupedMCParticles); // fills mc particles
    }
  }
  PROCESS_SWITCH(FemtoUniverseProducerTask, processTrackD0MC, "Provide MC data for track D0 analysis", false);

  void processFullMCCent(aod::FemtoFullCollisionCentRun3 const& col,
                         aod::BCsWithTimestamps const&,
                         soa::Join<aod::FemtoFullTracks, aod::McTrackLabels> const& tracks,
                         aod::McCollisions const&,
                         aod::McParticles const&)
  {
    // get magnetic field for run
    auto bc = col.bc_as<aod::BCsWithTimestamps>();
    getMagneticFieldTesla(bc);
    const double ir = mRateFetcher.fetch(ccdb.service, bc.timestamp(), mRunNumber, "ZNC hadronic") * 1.e-3; // fetch IR

    // fill the tables
    const auto colcheck = fillCollisionsCentRun3<true>(col);
    if (colcheck) {
      fillCollisionsCentRun3ColExtra<true>(col, ir);
      fillTracks<true>(tracks);
    }
  }
  PROCESS_SWITCH(FemtoUniverseProducerTask, processFullMCCent, "Provide MC data with centrality bins", false);

  void processTrackCentRun2Data(aod::FemtoFullCollisionCentRun2 const& col,
                                aod::BCsWithTimestamps const&,
                                soa::Filtered<aod::FemtoFullTracks> const& tracks)
  {
    // get magnetic field for run
    auto bc = col.bc_as<aod::BCsWithTimestamps>();
    getMagneticFieldTesla(bc);
    const double ir = 0.0; // fetch IR

    // fill the tables
    const auto colcheck = fillCollisionsCentRun2<false>(col);
    if (colcheck) {
      fillCollisionsCentRun3ColExtra<false>(col, ir);
      fillTracks<false>(tracks);
    }
  }
  PROCESS_SWITCH(FemtoUniverseProducerTask, processTrackCentRun2Data, "Provide experimental data for Run 2 with centrality for track track", false);

  void processTrackV0CentRun2Data(aod::FemtoFullCollisionCentRun2 const& col,
                                  aod::BCsWithTimestamps const&,
                                  soa::Filtered<aod::FemtoFullTracks> const& tracks,
                                  aod::V0Datas const& fullV0s)
  {
    // get magnetic field for run
    auto bc = col.bc_as<aod::BCsWithTimestamps>();
    getMagneticFieldTesla(bc);
    const double ir = 0.0; // fetch IR

    // fill the tables
    const auto colcheck = fillCollisionsCentRun2<false>(col);
    if (colcheck) {
      fillCollisionsCentRun3ColExtra<false>(col, ir);
      fillTracks<false>(tracks);
      fillV0<false>(col, fullV0s, tracks);
    }
  }
  PROCESS_SWITCH(FemtoUniverseProducerTask, processTrackV0CentRun2Data, "Provide experimental data for Run 2 with centrality for track V0", false);

  void processTrackCentRun3Data(aod::FemtoFullCollisionCentRun3 const& col,
                                aod::BCsWithTimestamps const&,
                                soa::Filtered<aod::FemtoFullTracks> const& tracks)
  {
    // get magnetic field for run
    auto bc = col.bc_as<aod::BCsWithTimestamps>();
    getMagneticFieldTesla(bc);
    const auto ir = mRateFetcher.fetch(ccdb.service, bc.timestamp(), mRunNumber, "ZNC hadronic") * 1.e-3; // fetch IR

    // fill the tables
    const auto colcheck = fillCollisionsCentRun3<false>(col);
    if (colcheck) {
      fillCollisionsCentRun3ColExtra<false>(col, ir);
      fillTracks<false>(tracks);
    }
  }
  PROCESS_SWITCH(FemtoUniverseProducerTask, processTrackCentRun3Data, "Provide experimental data for Run 3 with centrality for track track", false);

  void processTrackCentRun3DataMC(aod::FemtoFullCollisionCentRun3MC const& col,
                                  aod::BCsWithTimestamps const&,
                                  soa::Join<aod::FemtoFullTracks, aod::McTrackLabels> const& tracks,
                                  aod::McCollisions const&,
                                  aod::McParticles const&)
  {
    // get magnetic field for run
    auto bc = col.bc_as<aod::BCsWithTimestamps>();
    getMagneticFieldTesla(bc);
    const auto ir = mRateFetcher.fetch(ccdb.service, bc.timestamp(), mRunNumber, "ZNC hadronic") * 1.e-3; // fetch IR

    // fill the tables
    const auto colcheck = fillCollisionsCentRun3<true>(col);
    if (colcheck) {
      fillCollisionsCentRun3ColExtra<true>(col, ir);
      fillTracks<true>(tracks);
    }
  }
  PROCESS_SWITCH(FemtoUniverseProducerTask, processTrackCentRun3DataMC, "Provide MC data for track analysis", false);

  void processV0CentRun3Data(aod::FemtoFullCollisionCentRun3 const& col,
                             aod::BCsWithTimestamps const&,
                             soa::Filtered<aod::FemtoFullTracks> const& tracks,
                             aod::V0Datas const& fullV0s)
  {
    // get magnetic field for run
    auto bc = col.bc_as<aod::BCsWithTimestamps>();
    getMagneticFieldTesla(bc);
    const auto ir = mRateFetcher.fetch(ccdb.service, bc.timestamp(), mRunNumber, "ZNC hadronic") * 1.e-3; // fetch IR

    // fill the tables
    const auto colcheck = fillCollisionsCentRun3<false>(col);
    if (colcheck) {
      fillCollisionsCentRun3ColExtra<false>(col, ir);
      fillTracks<false>(tracks);
      fillV0<false>(col, fullV0s, tracks);
    }
  }
  PROCESS_SWITCH(FemtoUniverseProducerTask, processV0CentRun3Data, "Provide experimental data for Run 3 with centrality for track track", false);

  void processTruthAndFullMCCentRun3V0(
    aod::McCollisions const& mccols,
    aod::McParticles const& mcParticles,
    aod::FemtoFullCollisionCentRun3MCs const& collisions,
    soa::Filtered<soa::Join<aod::FemtoFullTracks, aod::McTrackLabels>> const& tracks,
    soa::Join<o2::aod::V0Datas, aod::McV0Labels> const& fullV0s,
    aod::BCsWithTimestamps const&)
  {

    // MCReco
    std::set<int> recoMcIds;
    std::set<int> mcCollisions;
    for (const auto& col : collisions) {                                          // loop over collisions
      auto groupedTracks = tracks.sliceBy(perCollisionTracks, col.globalIndex()); // slicing for tracks
      auto groupedV0Parts = fullV0s.sliceBy(perCollisionV0s, col.globalIndex());  // slicing for V0
      getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
      const auto colcheck = fillCollisionsCentRun3<false>(col);
      if (colcheck) {
        mcCollisions.insert(col.mcCollisionId());
        fillTracks<true>(groupedTracks);
        fillV0<true>(col, groupedV0Parts, groupedTracks);
      }
      for (const auto& track : groupedTracks) {
        if (trackCuts.isSelectedMinimal(track))
          recoMcIds.insert(track.mcParticleId());
      }
    }

    // MCTruth
    for (const auto& mccol : mccols) {
      if (confCollMCTruthOnlyReco && !mcCollisions.contains(mccol.globalIndex())) {
        continue;
      }
      auto groupedCollisions = collisions.sliceBy(recoCollsPerMCCollCentPbPb, mccol.globalIndex()); // slicing for MC collisions
      auto groupedMCParticles = mcParticles.sliceBy(perMCCollision, mccol.globalIndex());           // slicing for MC particles
      for (const auto& col : groupedCollisions) {
        const auto colcheck = fillMCTruthCollisionsCentRun3(col);
        if (colcheck) {
          outputCollExtra(1.0, 1.0);
          if (confFillMCTruthV0Daugh) {
            fillV0MCTruth(groupedMCParticles); // fills MC V0s and its daughters
          } else {
            fillParticles<decltype(groupedMCParticles), true, true>(groupedMCParticles, recoMcIds); // fills mc particles
          }
        }
      }
    }
  }
  PROCESS_SWITCH(FemtoUniverseProducerTask, processTruthAndFullMCCentRun3V0, "Provide both MC truth and reco for tracks and V0s with centrality", false);

  void processCascadeCentRun3Data(aod::FemtoFullCollisionCentRun3 const& col,
                                  aod::BCsWithTimestamps const&,
                                  soa::Filtered<aod::FemtoFullTracks> const& tracks,
                                  aod::CascDatas const& fullCascades)
  {
    // get magnetic field for run
    auto bc = col.bc_as<aod::BCsWithTimestamps>();
    getMagneticFieldTesla(bc);
    const auto ir = mRateFetcher.fetch(ccdb.service, bc.timestamp(), mRunNumber, "ZNC hadronic") * 1.e-3; // fetch IR

    // fill the tables
    const auto colcheck = fillCollisionsCentRun3<false>(col);
    if (colcheck) {
      fillCollisionsCentRun3ColExtra<false>(col, ir);
      fillTracks<false>(tracks);
      fillCascade<false>(col, fullCascades, tracks);
    }
  }
  PROCESS_SWITCH(FemtoUniverseProducerTask, processCascadeCentRun3Data, "Provide experimental data for Run 3 with centrality for track track", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<FemtoUniverseProducerTask>(cfgc)};
  return workflow;
}
