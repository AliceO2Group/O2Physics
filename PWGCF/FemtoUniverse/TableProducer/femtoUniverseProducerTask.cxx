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

/// \file femtoUniverseProducerTask.cxx
/// \brief Tasks that produces the track tables used for the pairing
/// \author Laura Serksnyte, TU München, laura.serksnyte@tum.de
/// \author Zuzanna Chochulska, WUT Warsaw & CTU Prague, zchochul@cern.ch
/// \author Malgorzata Janik, WUT Warsaw, majanik@cern.ch
/// \author Pritam Chakraborty, WUT Warsaw, pritam.chakraborty@cern.ch

#include <CCDB/BasicCCDBManager.h>
#include <TPDGCode.h>
#include <vector>
#include <algorithm>
#include <set>

#include "CommonConstants/PhysicsConstants.h"
#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseCollisionSelection.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseTrackSelection.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseV0Selection.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseCascadeSelection.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniversePhiSelection.h"
#include "PWGCF/FemtoUniverse/Core/femtoUtils.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Core/HfHelper.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Math/Vector4D.h"
#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "ReconstructionDataFormats/Track.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "Framework/O2DatabasePDGPlugin.h"

using namespace o2;
using namespace o2::analysis::femto_universe;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

namespace o2::aod
{

using FemtoFullCollision =
  soa::Join<aod::Collisions, aod::EvSels, aod::Mults>::iterator;
using FemtoFullCollisionCentPP =
  soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::Mults>::iterator;
using FemtoFullCollisionCentRun2 =
  soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms, aod::Mults>::iterator;
using FemtoFullCollisionCentRun3 =
  soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::Mults>::iterator;
using FemtoFullCollisionMC = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::McCollisionLabels>::iterator;
using FemtoFullCollisionCentRun3MC =
  soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::Mults, aod::McCollisionLabels>::iterator;
using FemtoFullTracks =
  soa::Join<aod::FullTracks, aod::TracksDCA, aod::TOFSignal, aod::pidTPCEl, aod::TrackSelection,
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
  // Choose if filtering or skimming version is run
  // Configurable<bool> confIsTrigger{"confIsTrigger", false, "Store all collisions"}; //Commented: not used configurable
  // Choose if running on converted data or Run3  / Pilot
  Configurable<bool> confIsRun3{"confIsRun3", true, "Running on Run3 or pilot"};
  // Configurable<bool> confIsMC{"confIsMC", false, "Running on MC; implemented only for Run3"}; //Commented: not used configurable

  Configurable<bool> confIsForceGRP{"confIsForceGRP", false, "Set true if the magnetic field configuration is not available in the usual CCDB directory (e.g. for Run 2 converted data or unanchorad Monte Carlo)"};

  Configurable<bool> confDoSpher{"confDoSpher", false, "Calculate sphericity. If false sphericity will take value of 2."};
  Configurable<bool> confFillCollExt{"confFillCollExt", false, "Option to fill collision extended table"};

  /// Event cuts
  FemtoUniverseCollisionSelection colCuts;
  Configurable<bool> confEvtUseTPCmult{"confEvtUseTPCmult", false, "Use multiplicity based on the number of tracks with TPC information"};
  Configurable<float> confEvtZvtx{"confEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};
  Configurable<bool> confEvtTriggerCheck{"confEvtTriggerCheck", true, "Evt sel: check for trigger"};
  Configurable<int> confEvtTriggerSel{"confEvtTriggerSel", kINT7, "Evt sel: trigger"};
  Configurable<bool> confEvtOfflineCheck{"confEvtOfflineCheck", false, "Evt sel: check for offline selection"};
  Configurable<bool> confIsActivateV0{"confIsActivateV0", false, "Activate filling of V0 into femtouniverse tables"};
  Configurable<bool> confActivateSecondaries{"confActivateSecondaries", false, "Fill secondary MC gen particles that were reconstructed"};
  Configurable<bool> confIsActivateCascade{"confIsActivateCascade", false, "Activate filling of Cascade into femtouniverse tables"};
  Configurable<bool> confIsSelectCascOmega{"confIsSelectCascOmega", false, "Select Omegas for cascade analysis"};
  Configurable<bool> confIsActivatePhi{"confIsActivatePhi", false, "Activate filling of Phi into femtouniverse tables"};
  Configurable<bool> confMCTruthAnalysisWithPID{"confMCTruthAnalysisWithPID", true, "1: take only particles with specified PDG, 0: all particles (for MC Truth)"};
  Configurable<std::vector<int>> confMCTruthPDGCodes{"confMCTruthPDGCodes", std::vector<int>{211, -211, 2212, -2212, 333}, "PDG of particles to be stored"};
  Configurable<float> confCentFT0Min{"confCentFT0Min", 0.f, "Min CentFT0 value for centrality selection"};
  Configurable<float> confCentFT0Max{"confCentFT0Max", 200.f, "Max CentFT0 value for centrality selection"};
  Configurable<bool> confEvIsGoodZvtxFT0vsPV{"confEvIsGoodZvtxFT0vsPV", true, "Require kIsGoodZvtxFT0vsPV selection on Events."};
  Configurable<bool> confEvNoSameBunchPileup{"confEvNoSameBunchPileup", true, "Require kNoSameBunchPileup selection on Events."};
  Configurable<bool> confIsUsePileUp{"confIsUsePileUp", true, "Required for choosing whether to run the pile-up cuts"};
  Configurable<bool> confEvIsVertexITSTPC{"confEvIsVertexITSTPC", true, "Require kIsVertexITSTPC selection on Events"};
  Configurable<int> confTPCOccupancyMin{"confTPCOccupancyMin", 0, "Minimum value for TPC Occupancy selection"};
  Configurable<int> confTPCOccupancyMax{"confTPCOccupancyMax", 500, "Maximum value for TPC Occupancy selection"};

  Filter customCollCentFilter = (aod::cent::centFT0C > confCentFT0Min) &&
                                (aod::cent::centFT0C < confCentFT0Max);

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
    // Configurable<bool> confIsFillV0s{"confIsFillV0s", false, "Choice to fill V0s"}; //Commented: not used configurable
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
  } ConfV0Selection;

  struct : o2::framework::ConfigurableGroup {
    Configurable<float> confPtLowFilterCut{"confPtLowFilterCut", 0.14, "Lower limit for Pt for the global track"};   // pT low
    Configurable<float> confPtHighFilterCut{"confPtHighFilterCut", 5.0, "Higher limit for Pt for the global track"}; // pT high
    Configurable<float> confEtaFilterCut{"confEtaFilterCut", 0.8, "Eta cut for the global track"};                   // eta
    Configurable<bool> confDxaXYCustom0Cut{"confDxaXYCustom0Cut", false, "Enable Custom Dcaxy < [0] cut."};
    Configurable<float> confDcaXYFilterCut{"confDcaXYFilterCut", 2.4, "Value for DCA_XY for the global track"}; // max dca to vertex XY
    Configurable<float> confDcaZFilterCut{"confDcaZFilterCut", 3.2, "Value for DCA_Z for the global track"};    // max dca to vertex Z
    Configurable<bool> confDcaXYCustom1Cut{"confDcaXYCustom1Cut", true, "Enable Custom |DCAxy| < [1] + [2]/pt cut."};
    Configurable<float> confDcaXYCustom11FilterCut{"confDcaXYCustom11FilterCut", 0.004, "Value for [1] custom DCAxy cut -> |DCAxy| < [1] + [2]/pT"};
    Configurable<float> confDcaXYCustom12FilterCut{"confDcaXYCustom12FilterCut", 0.013, "Value for [2] custom DCAxy cut -> |DCAxy| < [1] + [2]/pT"};
  } ConfFilterCuts;

  Filter globalCutFilter = requireGlobalTrackInFilter();
  Filter customTrackFilter = (aod::track::pt > ConfFilterCuts.confPtLowFilterCut) &&
                             (aod::track::pt < ConfFilterCuts.confPtHighFilterCut) &&
                             (nabs(aod::track::eta) < ConfFilterCuts.confEtaFilterCut) &&
                             (!ConfFilterCuts.confDxaXYCustom0Cut || (aod::track::dcaXY < ConfFilterCuts.confDcaXYFilterCut)) && // true if configurable set to false or if configurable is true and it passes the selection
                             (aod::track::dcaZ < ConfFilterCuts.confDcaZFilterCut) &&
                             (!ConfFilterCuts.confDcaXYCustom1Cut || (nabs(aod::track::dcaXY) < ConfFilterCuts.confDcaXYCustom11FilterCut + ConfFilterCuts.confDcaXYCustom12FilterCut / aod::track::pt)); // same logic here

  // CASCADE
  FemtoUniverseCascadeSelection cascadeCuts;
  struct : o2::framework::ConfigurableGroup {
    // Configurable<bool> confIsFillCascades{"confIsFillCascades", false, "Choice to fill cascades"}; //Commented: not used configurable
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

    Configurable<float> confCascInvMassLowLimit{"confCascInvMassLowLimit", 1.25, "Lower limit of the cascade invariant mass"};
    Configurable<float> confCascInvMassUpLimit{"confCascInvMassUpLimit", 1.40, "Upper limit of the cascade invariant mass"};

    Configurable<bool> confCascRejectCompetingMass{"confCascRejectCompetingMass", false, "Switch on to reject Omegas (for Xi) or Xis (for Omegas)"};
    Configurable<float> confCascInvCompetingMassLowLimit{"confCascInvCompetingMassLowLimit", 1.66, "Lower limit of the cascade invariant mass for competing mass rejection"};
    Configurable<float> confCascInvCompetingMassUpLimit{"confCascInvCompetingMassUpLimit", 1.68, "Upper limit of the cascade invariant mass for competing mass rejection"};
  } ConfCascadeSelection;

  // PHI
  FemtoUniversePhiSelection phiCuts;
  struct : o2::framework::ConfigurableGroup {
    /// Phi meson
    Configurable<float> confInvMassLowLimitPhi{"confInvMassLowLimitPhi", 1.011, "Lower limit of the Phi invariant mass"}; // change that to do invariant mass cut
    Configurable<float> confInvMassUpLimitPhi{"confInvMassUpLimitPhi", 1.027, "Upper limit of the Phi invariant mass"};
    Configurable<float> confPtLowLimitPhi{"confPtLowLimitPhi", 0.8, "Lower limit of the Phi pT."};
    Configurable<float> confPtHighLimitPhi{"confPtHighLimitPhi", 4.0, "Higher limit of the Phi pT."};
    // Phi meson daughters
    Configurable<bool> confLooseTPCNSigma{"confLooseTPCNSigma", false, "Use loose TPC N sigmas for Kaon PID."};
    Configurable<float> confLooseTPCNSigmaValue{"confLooseTPCNSigmaValue", 10, "Value for the loose TPC N Sigma for Kaon PID."};
    Configurable<bool> confLooseTOFNSigma{"confLooseTOFNSigma", false, "Use loose TPC N sigmas for Kaon PID."};
    Configurable<float> confNsigmaRejectPion{"confNsigmaRejectPion", 3.0, "Reject if particle could be a Pion combined nsigma value."};
    Configurable<float> confNsigmaRejectProton{"confNsigmaRejectProton", 3.0, "Reject if particle could be a Proton combined nsigma value."};
    Configurable<float> confLooseTOFNSigmaValue{"confLooseTOFNSigmaValue", 10, "Value for the loose TOF N Sigma for Kaon PID."};
  } ConfPhiSelection;

  // PDG codes for fillMCParticle function
  Configurable<int> confPDGCodePartOne{"confPDGCodePartOne", 321, "Particle 1 - PDG code"};
  Configurable<int> confPDGCodePartTwo{"confPDGCodePartTwo", 321, "Particle 2 - PDG code"};

  // D0/D0bar mesons
  struct : o2::framework::ConfigurableGroup {
    Configurable<float> confD0D0barCandMaxY{"confD0D0barCandMaxY", -1., "max. cand. rapidity"};
    Configurable<float> confD0D0barCandEtaCut{"confD0D0barCandEtaCut", 0.8, "max. cand. pseudorapidity"};
    Configurable<bool> storeD0D0barDoubleMassHypo{"storeD0D0barDoubleMassHypo", false, "Store D0/D0bar cand. which pass selection criteria for both, D0 and D0bar"};
    Configurable<bool> applyMLSelD0D0bar{"applyMLSelD0D0bar", false, "Use ML D0/D0bar selection"};
    Configurable<std::vector<int>> classMlD0D0bar{"classMlD0D0bar", {0, 1, 2}, "Indexes of ML scores to be stored. Three indexes max."};
  } ConfD0Selection;

  HfHelper hfHelper;
  bool isKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK)
  {

    if (mom < 0.3) { // 0.0-0.3
      if (std::abs(nsigmaTPCK) < 3.0) {
        return true;
      } else {
        return false;
      }
    } else if (mom < 0.45) { // 0.30 - 0.45
      if (std::abs(nsigmaTPCK) < 2.0) {
        return true;
      } else {
        return false;
      }
    } else if (mom < 0.55) { // 0.45-0.55
      if (std::abs(nsigmaTPCK) < 1.0) {
        return true;
      } else {
        return false;
      }
    } else if (mom < 1.5) { // 0.55-1.5 (now we use TPC and TOF)
      if ((std::abs(nsigmaTOFK) < 3.0) && (std::abs(nsigmaTPCK) < 3.0)) {
        {
          return true;
        }
      } else {
        return false;
      }
    } else if (mom > 1.5) { // 1.5 -
      if ((std::abs(nsigmaTOFK) < 2.0) && (std::abs(nsigmaTPCK) < 3.0)) {
        return true;
      } else {
        return false;
      }
    } else {
      return false;
    }
  }

  bool isKaonNSigmaTPCLoose(float mom, float nsigmaTPCK, float nsigmaTOFK)
  {

    if (mom < 0.3) { // 0.0-0.3
      if (std::abs(nsigmaTPCK) < ConfPhiSelection.confLooseTPCNSigmaValue.value) {
        return true;
      } else {
        return false;
      }
    } else if (mom < 0.45) { // 0.30 - 0.45
      if (std::abs(nsigmaTPCK) < ConfPhiSelection.confLooseTPCNSigmaValue.value) {
        return true;
      } else {
        return false;
      }
    } else if (mom < 0.55) { // 0.45-0.55
      if (std::abs(nsigmaTPCK) < ConfPhiSelection.confLooseTPCNSigmaValue.value) {
        return true;
      } else {
        return false;
      }
    } else if (mom < 1.5) { // 0.55-1.5 (now we use TPC and TOF)
      if ((std::abs(nsigmaTOFK) < 3.0) && (std::abs(nsigmaTPCK) < ConfPhiSelection.confLooseTPCNSigmaValue.value)) {
        return true;
      } else {
        return false;
      }
    } else if (mom > 1.5) { // 1.5 -
      if ((std::abs(nsigmaTOFK) < 2.0) && (std::abs(nsigmaTPCK) < ConfPhiSelection.confLooseTPCNSigmaValue.value)) {
        return true;
      } else {
        return false;
      }
    } else {
      return false;
    }
  }

  bool isKaonNSigmaTOFLoose(float mom, float nsigmaTPCK, float nsigmaTOFK)
  {
    if (mom < 0.3) { // 0.0-0.3
      if (std::abs(nsigmaTPCK) < 3.0) {
        return true;
      } else {
        return false;
      }
    } else if (mom < 0.45) { // 0.30 - 0.45
      if (std::abs(nsigmaTPCK) < 2.0) {
        return true;
      } else {
        return false;
      }
    } else if (mom < 0.55) { // 0.45-0.55
      if (std::abs(nsigmaTPCK) < 1.0) {
        return true;
      } else {
        return false;
      }
    } else if (mom < 1.5) { // 0.55-1.5 (now we use TPC and TOF)
      if ((std::abs(nsigmaTOFK) < ConfPhiSelection.confLooseTOFNSigmaValue.value) && (std::abs(nsigmaTPCK) < 3.0)) {
        {
          return true;
        }
      } else {
        return false;
      }
    } else if (mom > 1.5) { // 1.5 -
      if ((std::abs(nsigmaTOFK) < ConfPhiSelection.confLooseTOFNSigmaValue.value) && (std::abs(nsigmaTPCK) < 3.0)) {
        return true;
      } else {
        return false;
      }
    } else {
      return false;
    }
  }

  bool isKaonRejected(float mom, float nsigmaTPCPr, float nsigmaTOFPr, float nsigmaTPCPi, float nsigmaTOFPi)
  {
    if (mom < 0.5) {
      if (std::abs(nsigmaTPCPi) < ConfPhiSelection.confNsigmaRejectPion.value) {
        return true;
      } else if (std::abs(nsigmaTPCPr) < ConfPhiSelection.confNsigmaRejectProton.value) {
        return true;
      }
    }
    if (mom > 0.5) {
      if (std::hypot(nsigmaTOFPi, nsigmaTPCPi) < ConfPhiSelection.confNsigmaRejectPion.value) {
        return true;
      } else if (std::hypot(nsigmaTOFPr, nsigmaTPCPr) < ConfPhiSelection.confNsigmaRejectProton.value) {
        return true;
      } else {
        return false;
      }
    } else {
      return false;
    }
  }

  /// \todo should we add filter on min value pT/eta of V0 and daughters?
  /*Filter v0Filter = (nabs(aod::v0data::x) < V0DecVtxMax.value) &&
                    (nabs(aod::v0data::y) < V0DecVtxMax.value) &&
                    (nabs(aod::v0data::z) < V0DecVtxMax.value);*/
  // (aod::v0data::v0radius > V0TranRadV0Min.value); to be added, not working
  // for now do not know why

  HistogramRegistry qaRegistry{"QAHistos", {}, OutputObjHandlingPolicy::QAObject};
  HistogramRegistry cascadeQaRegistry{"CascadeQAHistos", {}, OutputObjHandlingPolicy::QAObject};

  int mRunNumber = 0;
  float mMagField;
  Service<o2::ccdb::BasicCCDBManager> ccdb; /// Accessing the CCDB
  ctpRateFetcher mRateFetcher;              // inspired by zdcSP.cxx in PWGLF

  void init(InitContext&)
  {
    if ((doprocessFullData || doprocessTrackPhiData || doprocessTrackData || doprocessTrackV0 || doprocessTrackCascadeData || doprocessTrackD0mesonData || doprocessTrackD0DataML || doprocessTrackCentRun2Data || doprocessTrackV0CentRun2Data || doprocessTrackCentRun3Data || doprocessV0CentRun3Data || doprocessCascadeCentRun3Data || doprocessTrackDataCentPP) == false && (doprocessFullMC || doprocessTrackMC || doprocessTrackMCTruth || doprocessTrackMCGen || doprocessTruthAndFullMC || doprocessFullMCCent || doprocessTrackCentRun3DataMC) == false) {
      LOGF(fatal, "Neither processFullData nor processFullMC enabled. Please choose one.");
    }
    if ((doprocessFullData || doprocessTrackPhiData || doprocessTrackData || doprocessTrackV0 || doprocessTrackCascadeData || doprocessTrackD0mesonData || doprocessTrackD0DataML || doprocessTrackCentRun2Data || doprocessTrackV0CentRun2Data || doprocessTrackCentRun3Data || doprocessV0CentRun3Data || doprocessCascadeCentRun3Data || doprocessTrackDataCentPP) == true && (doprocessFullMC || doprocessTrackMC || doprocessTrackMCTruth || doprocessTrackMCGen || doprocessTruthAndFullMC || doprocessFullMCCent || doprocessTrackCentRun3DataMC) == true) {
      LOGF(fatal,
           "Cannot enable process Data and process MC at the same time. "
           "Please choose one.");
    }

    colCuts.setCuts(confEvtZvtx, confEvtTriggerCheck, confEvtTriggerSel, confEvtOfflineCheck, confIsRun3, confCentFT0Min, confCentFT0Max);
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
    if (confIsActivateV0) {
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

    if (confIsActivateCascade) {
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
      cascadeCuts.init<aod::femtouniverseparticle::ParticleType::kCascade, aod::femtouniverseparticle::ParticleType::kV0Child, aod::femtouniverseparticle::ParticleType::kCascadeBachelor, aod::femtouniverseparticle::CutContainerType>(&cascadeQaRegistry, confIsSelectCascOmega);
      // invmass cuts
      cascadeCuts.setInvMassLimits(ConfCascadeSelection.confCascInvMassLowLimit, ConfCascadeSelection.confCascInvMassUpLimit);

      if (ConfCascadeSelection.confCascRejectCompetingMass) {
        cascadeCuts.setCompetingInvMassLimits(ConfCascadeSelection.confCascInvCompetingMassLowLimit, ConfCascadeSelection.confCascInvCompetingMassUpLimit);
      }
    }

    if (confIsActivatePhi) {
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

  template <bool isTrackOrV0, bool isPhiOrD0, bool isXi, typename ParticleType>
  void fillDebugParticle(ParticleType const& particle)
  {
    if constexpr (isTrackOrV0) {
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
    } else if constexpr (isPhiOrD0) {
      outputDebugParts(-999., -999., -999., -999., -999., -999., -999., -999., -999.,
                       -999., -999., -999., -999., -999., -999., -999., -999.,
                       -999., -999., -999., -999., -999.,
                       -999., -999.,
                       -999., -999., -999.,
                       -999.); // QA for phi or D0/D0bar
    } else if constexpr (isXi) {
      outputDebugParts(-999., -999., -999., -999., -999., -999., -999., -999., -999.,
                       -999., -999., -999., -999., -999., -999., -999., -999.,
                       -999., -999., -999., -999., -999.,
                       particle.dcacascdaughters(), particle.cascradius(),
                       particle.x(), particle.y(), particle.z(),
                       particle.mOmega()); // QA for Xi Cascades (later do the same for Omegas)
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
          if (motherparticleMC.producedByGenerator())
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
              if (particleMotherOfNeg == particleMotherOfPos && particleMotherOfNeg.pdgCode() == 333) {
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
    if (confEvtUseTPCmult) {
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
    if (!confIsUsePileUp) {
      outputCollision(vtxZ, mult, multNtr, confDoSpher ? colCuts.computeSphericity(col, tracks) : 2, mMagField);
      colCuts.fillQA(col);
      return true;
    } else if ((!confEvNoSameBunchPileup || col.selection_bit(aod::evsel::kNoSameBunchPileup)) && (!confEvIsGoodZvtxFT0vsPV || col.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) && (!confEvIsVertexITSTPC || col.selection_bit(aod::evsel::kIsVertexITSTPC))) {
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
    if (confEvtUseTPCmult) {
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
    if (!confIsUsePileUp) {
      outputCollision(vtxZ, mult, multNtr, confDoSpher ? colCuts.computeSphericity(col, tracks) : 2, mMagField);
      colCuts.fillQA(col);
      return true;
    } else if ((!confEvNoSameBunchPileup || col.selection_bit(aod::evsel::kNoSameBunchPileup)) && (!confEvIsGoodZvtxFT0vsPV || col.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) && (!confEvIsVertexITSTPC || col.selection_bit(aod::evsel::kIsVertexITSTPC))) {
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
      float mult = 0;
      int multNtr = 0;

      if (std::abs(vtxZ) > confEvtZvtx) {
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

    if (!colCuts.isSelectedRun3(col) || (occupancy < confTPCOccupancyMin || occupancy > confTPCOccupancyMax)) {
      return false;
    } else {
      if (col.selection_bit(aod::evsel::kNoSameBunchPileup) &&
          col.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) &&
          col.selection_bit(aod::evsel::kIsGoodITSLayersAll) &&
          col.selection_bit(aod::evsel::kNoCollInRofStandard) &&
          col.selection_bit(aod::evsel::kNoHighMultCollInPrevRof) &&
          col.selection_bit(aod::evsel::kNoCollInTimeRangeStandard)) {
        outputCollision(vtxZ, cent, multNtr, 2, mMagField);
        return true;
      } else {
        return false;
      }
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
      if (!confIsActivateCascade) {
        outputParts(outputCollision.lastIndex(), track.pt(), track.eta(),
                    track.phi(), aod::femtouniverseparticle::ParticleType::kTrack,
                    cutContainer.at(
                      femto_universe_track_selection::TrackContainerPosition::kCuts),
                    cutContainer.at(
                      femto_universe_track_selection::TrackContainerPosition::kPID),
                    track.dcaXY(), childIDs, 0,
                    track.sign()); // sign getter is mAntiLambda()
      } else {
        outputCascParts(outputCollision.lastIndex(), track.pt(), track.eta(),
                        track.phi(), aod::femtouniverseparticle::ParticleType::kTrack,
                        cutContainer.at(
                          femto_universe_track_selection::TrackContainerPosition::kCuts),
                        cutContainer.at(
                          femto_universe_track_selection::TrackContainerPosition::kPID),
                        track.dcaXY(), childIDs, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
      }
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
                  cutContainerV0.at(femto_universe_v0_selection::V0ContainerPosition::kPosPID),
                  0.,
                  childIDs,
                  0,
                  0);
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
                  cutContainerV0.at(femto_universe_v0_selection::V0ContainerPosition::kNegPID),
                  0.,
                  childIDs,
                  0,
                  0);
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
                  0,
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
      outputCascParts(outputCollision.lastIndex(),
                      casc.positivept(),
                      casc.positiveeta(),
                      casc.positivephi(),
                      aod::femtouniverseparticle::ParticleType::kV0Child,
                      0, // cutContainerV0.at(femto_universe_v0_selection::V0ContainerPosition::kPosCuts),
                      0, // cutContainerV0.at(femto_universe_v0_selection::V0ContainerPosition::kPosPID),
                      0.,
                      childIDs,
                      0,
                      0,
                      0,
                      0,
                      0,
                      0,
                      0,
                      0,
                      0,
                      0,
                      0,
                      0);
      const int rowOfPosTrack = outputCascParts.lastIndex();
      // if constexpr (isMC) {
      //   fillMCParticle(postrack, o2::aod::femtouniverseparticle::ParticleType::kV0Child);
      // }
      int negtrackID = casc.negTrackId();
      int rowInPrimaryTrackTableNeg = -1;
      rowInPrimaryTrackTableNeg = getRowDaughters(negtrackID, tmpIDtrack);
      childIDs[0] = 0;                         // pos
      childIDs[1] = rowInPrimaryTrackTableNeg; // neg
      childIDs[2] = 0;                         // bachelor
      outputCascParts(outputCollision.lastIndex(),
                      casc.negativept(),
                      casc.negativeeta(),
                      casc.negativephi(),
                      aod::femtouniverseparticle::ParticleType::kV0Child,
                      0, // cutContainerV0.at(femto_universe_v0_selection::V0ContainerPosition::kNegCuts),
                      0, // cutContainerV0.at(femto_universe_v0_selection::V0ContainerPosition::kNegPID),
                      0.,
                      childIDs,
                      0,
                      0,
                      0,
                      0,
                      0,
                      0,
                      0,
                      0,
                      0,
                      0,
                      0,
                      0);
      const int rowOfNegTrack = outputCascParts.lastIndex();
      // if constexpr (isMC) {
      //   fillMCParticle(negtrack, o2::aod::femtouniverseparticle::ParticleType::kV0Child);
      // }
      //  bachelor
      int bachtrackID = casc.bachelorId();
      int rowInPrimaryTrackTableBach = -1;
      rowInPrimaryTrackTableBach = getRowDaughters(bachtrackID, tmpIDtrack);
      childIDs[0] = 0;                          // pos
      childIDs[1] = 0;                          // neg
      childIDs[2] = rowInPrimaryTrackTableBach; // bachelor
      outputCascParts(outputCollision.lastIndex(),
                      casc.bachelorpt(),
                      casc.bacheloreta(),
                      casc.bachelorphi(),
                      aod::femtouniverseparticle::ParticleType::kCascadeBachelor,
                      0, // cutContainerV0.at(femto_universe_v0_selection::V0ContainerPosition::kNegCuts),
                      0, // cutContainerV0.at(femto_universe_v0_selection::V0ContainerPosition::kNegPID),
                      0.,
                      childIDs,
                      0,
                      0,
                      0,
                      0,
                      0,
                      0,
                      0,
                      0,
                      0,
                      0,
                      0,
                      0);
      const int rowOfBachTrack = outputCascParts.lastIndex();
      // cascade
      std::vector<int> indexCascChildID = {rowOfPosTrack, rowOfNegTrack, rowOfBachTrack};
      outputCascParts(outputCollision.lastIndex(),
                      casc.pt(),
                      casc.eta(),
                      casc.phi(),
                      aod::femtouniverseparticle::ParticleType::kCascade,
                      0, // cutContainerV0.at(femto_universe_v0_selection::V0ContainerPosition::kV0),
                      0,
                      0,
                      indexCascChildID,
                      confIsSelectCascOmega ? casc.mOmega() : casc.mXi(),
                      confIsSelectCascOmega ? casc.mOmega() : casc.mXi(),
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
    }
  }

  template <bool isMC, typename HfCandidate, typename TrackType, typename CollisionType>
  void fillD0mesons(CollisionType const&, TrackType const&, HfCandidate const& hfCands)
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

      if (ConfD0Selection.confD0D0barCandMaxY >= 0. && std::abs(hfHelper.yD0(hfCand)) > ConfD0Selection.confD0D0barCandMaxY) {
        continue;
      }

      if (std::abs(hfCand.eta()) > ConfD0Selection.confD0D0barCandEtaCut) {
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
  void fillD0D0barUsingML(CollisionType const&, TrackType const&, HfCandidate const& hfCands)
  {
    std::vector<int> childIDs = {0, 0}; // these IDs are necessary to keep track of the children
    std::vector<int> tmpIDtrack;        // this vector keeps track of the matching of the primary track table row <-> aod::track table global index
    double invMassD0 = 0.0;
    double invMassD0bar = 0.0;
    bool isD0D0bar = false;
    std::vector<float> outputMlD0D0bar = {-1., -1., -1.}; // this vector keeps the probabilities from the ML model for D0/D0bar
    uint8_t daughFlag = 0;                                // flag = 0 (daugh of D0 or D0bar), 1 (daug of D0), -1 (daugh of D0bar)

    for (const auto& hfCand : hfCands) {

      if (!(hfCand.hfflag() & 1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        continue;
      }

      if (ConfD0Selection.confD0D0barCandMaxY >= 0. && std::abs(hfHelper.yD0(hfCand)) > ConfD0Selection.confD0D0barCandMaxY) {
        continue;
      }

      if (std::abs(hfCand.eta()) > ConfD0Selection.confD0D0barCandEtaCut) {
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
        if (ConfD0Selection.applyMLSelD0D0bar) {
          for (unsigned int iClass = 0; iClass < ConfD0Selection.classMlD0D0bar->size(); iClass++) {
            outputMlD0D0bar[iClass] = hfCand.mlProbD0()[ConfD0Selection.classMlD0D0bar->at(iClass)];
          }
        }
        isD0D0bar = true;
        daughFlag = 1;
      } else if (hfCand.isSelD0() == 0 && hfCand.isSelD0bar() == 1) {
        invMassD0 = -hfHelper.invMassD0ToPiK(hfCand);
        invMassD0bar = hfHelper.invMassD0barToKPi(hfCand);
        if (ConfD0Selection.applyMLSelD0D0bar) {
          for (unsigned int iClass = 0; iClass < ConfD0Selection.classMlD0D0bar->size(); iClass++) {
            outputMlD0D0bar[iClass] = hfCand.mlProbD0bar()[ConfD0Selection.classMlD0D0bar->at(iClass)];
          }
        }
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
                    -999,               // cut, CutContainerType
                    -999,               // PID, CutContainerType
                    outputMlD0D0bar[0], // saving only the probaility for store class 1 - background
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
      // implementing PID cuts for phi children
      if (ConfPhiSelection.confLooseTPCNSigma.value) {
        if (!(isKaonNSigmaTPCLoose(p1.pt(), trackCuts.getNsigmaTPC(p1, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p1, o2::track::PID::Kaon)))) {
          continue;
        }
        if (!(isKaonNSigmaTPCLoose(p2.pt(), trackCuts.getNsigmaTPC(p2, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p2, o2::track::PID::Kaon)))) {
          continue;
        }
      }
      if (ConfPhiSelection.confLooseTOFNSigma.value) {
        if (!(isKaonNSigmaTOFLoose(p1.pt(), trackCuts.getNsigmaTPC(p1, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p1, o2::track::PID::Kaon)))) {
          continue;
        }
        if (!(isKaonNSigmaTOFLoose(p2.pt(), trackCuts.getNsigmaTPC(p2, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p2, o2::track::PID::Kaon)))) {
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

      if ((!(p1.sign() == 1)) || (!(p2.sign() == -1))) {
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
      if (std::abs(phiEta) > 0.8) {
        continue;
      }

      float phiPt = sumVec.Pt();
      if ((phiPt < ConfPhiSelection.confPtLowLimitPhi.value) || (phiPt > ConfPhiSelection.confPtHighLimitPhi.value)) {
        continue;
      }

      float phiPhi = RecoDecay::constrainAngle(sumVec.Phi(), 0);
      float phiM = sumVec.M();

      if (((phiM < ConfPhiSelection.confInvMassLowLimitPhi.value) || (phiM > ConfPhiSelection.confInvMassUpLimitPhi.value))) {
        continue;
      }

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

      if (confMCTruthAnalysisWithPID) {
        bool pass = false;
        std::vector<int> tmpPDGCodes = confMCTruthPDGCodes; // necessary due to some features of the Configurable
        for (auto const& pdg : tmpPDGCodes) {
          if (static_cast<int>(pdg) == static_cast<int>(pdgCode)) {
            if (pdgCode == 333) { // && (recoMcIds && recoMcIds->get().contains(particle.globalIndex()))) { // ATTENTION: all Phi mesons are NOT primary particles
              pass = true;
            } else {
              if (particle.isPhysicalPrimary() || (confActivateSecondaries && recoMcIds && recoMcIds->get().contains(particle.globalIndex())))
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

      // now the table is filled
      if constexpr (resolveDaughs) {
        tmpIDtrack.push_back(particle.globalIndex());
        continue;
      }
      if (!confIsActivateCascade) {
        outputParts(outputCollision.lastIndex(),
                    particle.pt(),
                    particle.eta(),
                    particle.phi(),
                    aod::femtouniverseparticle::ParticleType::kMCTruthTrack,
                    0,
                    pdgCode,
                    pdgCode,
                    childIDs,
                    0,
                    0);
      } else {
        outputCascParts(outputCollision.lastIndex(),
                        particle.pt(),
                        particle.eta(),
                        particle.phi(),
                        aod::femtouniverseparticle::ParticleType::kMCTruthTrack,
                        0,
                        pdgCode,
                        pdgCode,
                        childIDs, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
      }
      if (confIsDebug) {
        fillDebugParticle<false, true, false>(particle);
      }

      // Workaround to keep the FDParticles and MC label tables
      // aligned, so that they can be joined in the task.
      if constexpr (transientLabels) {
        outputPartsMCLabels(-1);
      }
    }
    if constexpr (resolveDaughs) {
      childIDs[0] = 0;
      childIDs[1] = 0;
      for (std::size_t i = 0; i < tmpIDtrack.size(); i++) {
        const auto& particle = tracks.iteratorAt(tmpIDtrack[i] - tracks.begin().globalIndex());
        for (int daughIndex = 0, n = std::min(2ul, particle.daughtersIds().size()); daughIndex < n; daughIndex++) {
          // loop to find the corresponding index of the daughters
          for (std::size_t j = 0; j < tmpIDtrack.size(); j++) {
            if (tmpIDtrack[j] == particle.daughtersIds()[daughIndex]) {
              childIDs[daughIndex] = i - j;
              break;
            }
          }
        }
        outputParts(outputCollision.lastIndex(),
                    particle.pt(),
                    particle.eta(),
                    particle.phi(),
                    aod::femtouniverseparticle::ParticleType::kMCTruthTrack,
                    0,
                    static_cast<uint32_t>(particle.pdgCode()),
                    particle.pdgCode(),
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
        }
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
      if (confIsActivateV0) {
        fillV0<isMC>(col, fullV0s, tracks);
      }
      if (confIsActivatePhi) {
        fillPhi<isMC>(col, tracks);
      }
    }
    // if (confIsActivateCascade) {
    //   fillCascade<false>(col, fullCascades, tracks);
    // }
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

  void processTrackD0MC(aod::FemtoFullCollisionMC const& col,
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
      // fillD0mesons<true>(col, tracks, candidates);
    }
  }
  PROCESS_SWITCH(FemtoUniverseProducerTask, processTrackD0MC, "Provide MC data for track D0 analysis", false);

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
      fillD0mesons<false>(col, tracks, candidates);
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
      fillD0D0barUsingML<false>(col, tracks, candidates);
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

  Preslice<aod::McParticles> perMCCollision = aod::mcparticle::mcCollisionId;
  PresliceUnsorted<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::McCollisionLabels>> recoCollsPerMCColl = aod::mcparticle::mcCollisionId;
  Preslice<soa::Join<aod::FemtoFullTracks, aod::McTrackLabels>> perCollisionTracks = aod::track::collisionId;
  Preslice<soa::Join<o2::aod::V0Datas, aod::McV0Labels>> perCollisionV0s = aod::track::collisionId;
  void processTruthAndFullMC(
    aod::McCollisions const& mccols,
    aod::McParticles const& mcParticles,
    soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::McCollisionLabels> const& collisions,
    soa::Filtered<soa::Join<aod::FemtoFullTracks, aod::McTrackLabels>> const& tracks,
    soa::Join<o2::aod::V0Datas, aod::McV0Labels> const& fullV0s,
    aod::BCsWithTimestamps const&)
  {
    // recos
    std::set<int> recoMcIds;
    for (const auto& col : collisions) {
      auto groupedTracks = tracks.sliceBy(perCollisionTracks, col.globalIndex());
      auto groupedV0s = fullV0s.sliceBy(perCollisionV0s, col.globalIndex());
      getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
      fillCollisionsAndTracksAndV0AndPhi<true>(col, groupedTracks, groupedV0s);
      for (const auto& track : groupedTracks) {
        if (trackCuts.isSelectedMinimal(track))
          recoMcIds.insert(track.mcParticleId());
      }
    }

    // truth
    for (const auto& mccol : mccols) {
      auto groupedMCParticles = mcParticles.sliceBy(perMCCollision, mccol.globalIndex());
      auto groupedCollisions = collisions.sliceBy(recoCollsPerMCColl, mccol.globalIndex());
      fillMCTruthCollisions(groupedCollisions, groupedMCParticles);                           // fills the reco collisions for mc collision
      fillParticles<decltype(groupedMCParticles), true, true>(groupedMCParticles, recoMcIds); // fills mc particles
    }
  }
  PROCESS_SWITCH(FemtoUniverseProducerTask, processTruthAndFullMC, "Provide both MC truth and reco for tracks and V0s", false);

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
