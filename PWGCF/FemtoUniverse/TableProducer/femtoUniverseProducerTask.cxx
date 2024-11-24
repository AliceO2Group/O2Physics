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
#include <TDatabasePDG.h> // FIXME
#include <vector>
#include <algorithm>
#include <set>

#include "Common/Core/trackUtilities.h"
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
#include "PWGCF/FemtoUniverse/Core/FemtoUtils.h"
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
using namespace o2::analysis::femtoUniverse;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{

using FemtoFullCollision =
  soa::Join<aod::Collisions, aod::EvSels, aod::Mults>::iterator;
using FemtoFullCollisionCentRun2 =
  soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms, aod::Mults>::iterator;
using FemtoFullCollisionCentRun3 =
  soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::Mults>::iterator;
using FemtoFullCollisionMC = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::McCollisionLabels>::iterator;

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

struct femtoUniverseProducerTask {
  Produces<aod::FDCollisions> outputCollision;
  Produces<aod::FDParticles> outputParts;
  Produces<aod::FDMCParticles> outputPartsMC;
  Produces<aod::FDExtParticles> outputDebugParts;
  Produces<aod::FDMCLabels> outputPartsMCLabels;
  Produces<aod::FDExtMCParticles> outputDebugPartsMC;
  Produces<aod::FDCascParticles> outputCascParts;

  Configurable<bool> ConfIsDebug{"ConfIsDebug", true, "Enable Debug tables"};
  // Choose if filtering or skimming version is run
  Configurable<bool> ConfIsTrigger{"ConfIsTrigger", false, "Store all collisions"};
  // Choose if running on converted data or Run3  / Pilot
  Configurable<bool> ConfIsRun3{"ConfIsRun3", true, "Running on Run3 or pilot"};
  Configurable<bool> ConfIsMC{"ConfIsMC", false, "Running on MC; implemented only for Run3"};

  Configurable<bool> ConfIsForceGRP{"ConfIsForceGRP", false, "Set true if the magnetic field configuration is not available in the usual CCDB directory (e.g. for Run 2 converted data or unanchorad Monte Carlo)"};

  Configurable<bool> ConfDoSpher{"ConfDoSpher", false, "Calculate sphericity. If false sphericity will take value of 2."};

  /// Event cuts
  FemtoUniverseCollisionSelection colCuts;
  Configurable<bool> ConfEvtUseTPCmult{"ConfEvtUseTPCmult", false, "Use multiplicity based on the number of tracks with TPC information"};
  Configurable<float> ConfEvtZvtx{"ConfEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};
  Configurable<bool> ConfEvtTriggerCheck{"ConfEvtTriggerCheck", true, "Evt sel: check for trigger"};
  Configurable<int> ConfEvtTriggerSel{"ConfEvtTriggerSel", kINT7, "Evt sel: trigger"};
  Configurable<bool> ConfEvtOfflineCheck{"ConfEvtOfflineCheck", false, "Evt sel: check for offline selection"};
  Configurable<bool> ConfIsActivateV0{"ConfIsActivateV0", false, "Activate filling of V0 into femtouniverse tables"};
  Configurable<bool> ConfActivateSecondaries{"ConfActivateSecondaries", false, "Fill secondary MC gen particles that were reconstructed"};
  Configurable<bool> ConfIsActivateCascade{"ConfIsActivateCascade", false, "Activate filling of Cascade into femtouniverse tables"};
  Configurable<bool> ConfIsSelectCascOmega{"ConfIsSelectCascOmega", false, "Select Omegas for cascade analysis"};
  Configurable<bool> ConfIsActivatePhi{"ConfIsActivatePhi", false, "Activate filling of Phi into femtouniverse tables"};
  Configurable<bool> ConfMCTruthAnalysisWithPID{"ConfMCTruthAnalysisWithPID", true, "1: take only particles with specified PDG, 0: all particles (for MC Truth)"};
  Configurable<std::vector<int>> ConfMCTruthPDGCodes{"ConfMCTruthPDGCodes", std::vector<int>{211, -211, 2212, -2212, 333}, "PDG of particles to be stored"};
  Configurable<float> ConfCentFT0Min{"ConfCentFT0Min", 0.f, "Min CentFT0 value for centrality selection"};
  Configurable<float> ConfCentFT0Max{"ConfCentFT0Max", 200.f, "Max CentFT0 value for centrality selection"};
  Configurable<bool> ConfEvIsGoodZvtxFT0vsPV{"ConfEvIsGoodZvtxFT0vsPV", true, "Require kIsGoodZvtxFT0vsPV selection on Events."};
  Configurable<bool> ConfEvNoSameBunchPileup{"ConfEvNoSameBunchPileup", true, "Require kNoSameBunchPileup selection on Events."};
  Configurable<bool> ConfIsUsePileUp{"ConfIsUsePileUp", true, "Required for choosing whether to run the pile-up cuts"};
  Configurable<bool> ConfEvIsVertexITSTPC{"ConfEvIsVertexITSTPC", true, "Require kIsVertexITSTPC selection on Events"};
  Configurable<int> ConfTPCOccupancyMin{"ConfTPCOccupancyMin", 0, "Minimum value for TPC Occupancy selection"};
  Configurable<int> ConfTPCOccupancyMax{"ConfTPCOccupancyMax", 500, "Maximum value for TPC Occupancy selection"};

  Filter CustomCollCentFilter = (aod::cent::centFT0C > ConfCentFT0Min) &&
                                (aod::cent::centFT0C < ConfCentFT0Max);

  // just sanity check to make sure in case there are problems in conversion or
  // MC production it does not affect results
  Configurable<bool> ConfTrkRejectNotPropagated{"ConfTrkRejectNotPropagated", false, "True: reject not propagated tracks"};
  // Configurable<bool> ConfRejectITSHitandTOFMissing{
  //     "ConfRejectITSHitandTOFMissing", false,
  //     "True: reject if neither ITS hit nor TOF timing satisfied"};

  FemtoUniverseTrackSelection trackCuts;
  Configurable<std::vector<float>> ConfTrkCharge{FemtoUniverseTrackSelection::getSelectionName(femtoUniverseTrackSelection::kSign, "ConfTrk"), std::vector<float>{-1, 1}, FemtoUniverseTrackSelection::getSelectionHelper(femtoUniverseTrackSelection::kSign, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkPtmin{FemtoUniverseTrackSelection::getSelectionName(femtoUniverseTrackSelection::kpTMin, "ConfTrk"), std::vector<float>{0.5f, 0.4f, 0.6f}, FemtoUniverseTrackSelection::getSelectionHelper(femtoUniverseTrackSelection::kpTMin, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkPtmax{FemtoUniverseTrackSelection::getSelectionName(femtoUniverseTrackSelection::kpTMax, "ConfTrk"), std::vector<float>{5.4f, 5.6f, 5.5f}, FemtoUniverseTrackSelection::getSelectionHelper(femtoUniverseTrackSelection::kpTMax, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkEta{FemtoUniverseTrackSelection::getSelectionName(femtoUniverseTrackSelection::kEtaMax, "ConfTrk"), std::vector<float>{0.8f, 0.7f, 0.9f}, FemtoUniverseTrackSelection::getSelectionHelper(femtoUniverseTrackSelection::kEtaMax, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkTPCnclsMin{FemtoUniverseTrackSelection::getSelectionName(femtoUniverseTrackSelection::kTPCnClsMin, "ConfTrk"), std::vector<float>{70.f}, FemtoUniverseTrackSelection::getSelectionHelper(femtoUniverseTrackSelection::kTPCnClsMin, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkTPCfCls{FemtoUniverseTrackSelection::getSelectionName(femtoUniverseTrackSelection::kTPCfClsMin, "ConfTrk"), std::vector<float>{0.83f}, FemtoUniverseTrackSelection::getSelectionHelper(femtoUniverseTrackSelection::kTPCfClsMin, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkTPCcRowsMin{FemtoUniverseTrackSelection::getSelectionName(femtoUniverseTrackSelection::kTPCcRowsMin, "ConfTrk"), std::vector<float>{70.f, 60.f, 80.f}, FemtoUniverseTrackSelection::getSelectionHelper(femtoUniverseTrackSelection::kTPCcRowsMin, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkTPCsCls{FemtoUniverseTrackSelection::getSelectionName(femtoUniverseTrackSelection::kTPCsClsMax, "ConfTrk"), std::vector<float>{0.1f, 160.f}, FemtoUniverseTrackSelection::getSelectionHelper(femtoUniverseTrackSelection::kTPCsClsMax, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkITSnclsMin{FemtoUniverseTrackSelection::getSelectionName(femtoUniverseTrackSelection::kITSnClsMin, "ConfTrk"), std::vector<float>{-1.f, 2.f, 4.f}, FemtoUniverseTrackSelection::getSelectionHelper(femtoUniverseTrackSelection::kITSnClsMin, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkITSnclsIbMin{FemtoUniverseTrackSelection::getSelectionName(femtoUniverseTrackSelection::kITSnClsIbMin, "ConfTrk"), std::vector<float>{-1.f, 1.f}, FemtoUniverseTrackSelection::getSelectionHelper(femtoUniverseTrackSelection::kITSnClsIbMin, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkDCAxyMax{FemtoUniverseTrackSelection::getSelectionName(femtoUniverseTrackSelection::kDCAxyMax, "ConfTrk"), std::vector<float>{0.1f, 3.5f}, FemtoUniverseTrackSelection::getSelectionHelper(femtoUniverseTrackSelection::kDCAxyMax, "Track selection: ")};
  Configurable<std::vector<float>> ConfTrkDCAzMax{FemtoUniverseTrackSelection::getSelectionName(femtoUniverseTrackSelection::kDCAzMax, "ConfTrk"), std::vector<float>{0.2f}, FemtoUniverseTrackSelection::getSelectionHelper(femtoUniverseTrackSelection::kDCAzMax, "Track selection: ")}; /// \todo Reintegrate PID to the general selection container
  Configurable<std::vector<float>> ConfTrkPIDnSigmaMax{FemtoUniverseTrackSelection::getSelectionName(femtoUniverseTrackSelection::kPIDnSigmaMax, "ConfTrk"), std::vector<float>{3.5f, 3.f, 2.5f}, FemtoUniverseTrackSelection::getSelectionHelper(femtoUniverseTrackSelection::kPIDnSigmaMax, "Track selection: ")};
  Configurable<float> ConfTrkPIDnSigmaOffsetTPC{"ConfTrkPIDnSigmaOffsetTPC", 0., "Offset for TPC nSigma because of bad calibration"};
  Configurable<float> ConfTrkPIDnSigmaOffsetTOF{"ConfTrkPIDnSigmaOffsetTOF", 0., "Offset for TOF nSigma because of bad calibration"};
  Configurable<std::vector<int>> ConfTrkPIDspecies{"ConfTrkPIDspecies", std::vector<int>{o2::track::PID::Pion, o2::track::PID::Kaon, o2::track::PID::Proton, o2::track::PID::Deuteron}, "Trk sel: Particles species for PID (Pion=2, Kaon=3, Proton=4, Deuteron=5)"};
  // Numbers from ~/alice/O2/DataFormats/Reconstruction/include/ReconstructionDataFormats/PID.h //static constexpr ID Pion = 2; static constexpr ID Kaon = 3; static constexpr ID Proton = 4; static constexpr ID Deuteron = 5;
  Configurable<float> ConfTOFpTmin{"ConfTOFpTmin", 500, "TOF pT min"};

  // TrackSelection *o2PhysicsTrackSelection;
  /// \todo Labeled array (see Track-Track task)

  // V0
  FemtoUniverseV0Selection v0Cuts;
  struct : o2::framework::ConfigurableGroup {
    Configurable<std::vector<float>> ConfV0Sign{FemtoUniverseV0Selection::getSelectionName(femtoUniverseV0Selection::kV0Sign, "ConfV0"), std::vector<float>{-1, 1}, FemtoUniverseV0Selection::getSelectionHelper(femtoUniverseV0Selection::kV0Sign, "V0 selection: ")};
    Configurable<std::vector<float>> ConfV0PtMin{FemtoUniverseV0Selection::getSelectionName(femtoUniverseV0Selection::kV0pTMin, "ConfV0"), std::vector<float>{0.3f, 0.4f, 0.5f}, FemtoUniverseV0Selection::getSelectionHelper(femtoUniverseV0Selection::kV0pTMin, "V0 selection: ")};
    Configurable<std::vector<float>> ConfV0PtMax{FemtoUniverseV0Selection::getSelectionName(femtoUniverseV0Selection::kV0pTMax, "ConfV0"), std::vector<float>{3.3f, 3.4f, 3.5f}, FemtoUniverseV0Selection::getSelectionHelper(femtoUniverseV0Selection::kV0pTMax, "V0 selection: ")};
    Configurable<std::vector<float>> ConfV0EtaMax{FemtoUniverseV0Selection::getSelectionName(femtoUniverseV0Selection::kV0etaMax, "ConfV0"), std::vector<float>{0.8f, 0.7f, 0.9f}, FemtoUniverseV0Selection::getSelectionHelper(femtoUniverseV0Selection::kV0etaMax, "V0 selection: ")};
    Configurable<std::vector<float>> ConfV0DCADaughMax{FemtoUniverseV0Selection::getSelectionName(femtoUniverseV0Selection::kV0DCADaughMax, "ConfV0"), std::vector<float>{1.2f, 1.5f}, FemtoUniverseV0Selection::getSelectionHelper(femtoUniverseV0Selection::kV0DCADaughMax, "V0 selection: ")};
    Configurable<std::vector<float>> ConfV0CPAMin{FemtoUniverseV0Selection::getSelectionName(femtoUniverseV0Selection::kV0CPAMin, "ConfV0"), std::vector<float>{0.99f, 0.995f}, FemtoUniverseV0Selection::getSelectionHelper(femtoUniverseV0Selection::kV0CPAMin, "V0 selection: ")};
    Configurable<std::vector<float>> ConfV0TranRadMin{FemtoUniverseV0Selection::getSelectionName(femtoUniverseV0Selection::kV0TranRadMin, "ConfV0"), std::vector<float>{0.2f}, FemtoUniverseV0Selection::getSelectionHelper(femtoUniverseV0Selection::kV0TranRadMin, "V0 selection: ")};
    Configurable<std::vector<float>> ConfV0TranRadMax{FemtoUniverseV0Selection::getSelectionName(femtoUniverseV0Selection::kV0TranRadMax, "ConfV0"), std::vector<float>{100.f}, FemtoUniverseV0Selection::getSelectionHelper(femtoUniverseV0Selection::kV0TranRadMax, "V0 selection: ")};
    Configurable<std::vector<float>> ConfV0DecVtxMax{FemtoUniverseV0Selection::getSelectionName(femtoUniverseV0Selection::kV0DecVtxMax, "ConfV0"), std::vector<float>{100.f}, FemtoUniverseV0Selection::getSelectionHelper(femtoUniverseV0Selection::kV0DecVtxMax, "V0 selection: ")};

    Configurable<std::vector<float>> ConfChildCharge{"ConfChildSign", std::vector<float>{-1, 1}, "V0 Child sel: Charge"};
    Configurable<std::vector<float>> ConfChildEtaMax{"ConfChildEtaMax", std::vector<float>{0.8f}, "V0 Child sel: max eta"};
    Configurable<std::vector<float>> ConfChildTPCnClsMin{"ConfChildTPCnClsMin", std::vector<float>{80.f, 70.f, 60.f}, "V0 Child sel: Min. nCls TPC"};
    Configurable<std::vector<float>> ConfChildDCAMin{"ConfChildDCAMin", std::vector<float>{0.05f, 0.06f}, "V0 Child sel:  Max. DCA Daugh to PV (cm)"};
    Configurable<std::vector<float>> ConfChildPIDnSigmaMax{"ConfChildPIDnSigmaMax", std::vector<float>{5.f, 4.f}, "V0 Child sel: Max. PID nSigma TPC"};
    Configurable<std::vector<int>> ConfChildPIDspecies{"ConfChildPIDspecies", std::vector<int>{o2::track::PID::Pion, o2::track::PID::Proton}, "V0 Child sel: Particles species for PID"};

    Configurable<float> ConfV0InvMassLowLimit{"ConfV0InvV0MassLowLimit", 1.05, "Lower limit of the V0 invariant mass"};
    Configurable<float> ConfV0InvMassUpLimit{"ConfV0InvV0MassUpLimit", 1.30, "Upper limit of the V0 invariant mass"};

    Configurable<bool> ConfV0RejectKaons{"ConfV0RejectKaons", false, "Switch to reject kaons"};
    Configurable<float> ConfV0InvKaonMassLowLimit{"ConfV0InvKaonMassLowLimit", 0.48, "Lower limit of the V0 invariant mass for Kaon rejection"};
    Configurable<float> ConfV0InvKaonMassUpLimit{"ConfV0InvKaonMassUpLimit", 0.515, "Upper limit of the V0 invariant mass for Kaon rejection"};
  } ConfV0Selection;

  struct : o2::framework::ConfigurableGroup {
    Configurable<float> ConfPtLowFilterCut{"ConfPtLowFilterCut", 0.14, "Lower limit for Pt for the global track"};   // pT low
    Configurable<float> ConfPtHighFilterCut{"ConfPtHighFilterCut", 5.0, "Higher limit for Pt for the global track"}; // pT high
    Configurable<float> ConfEtaFilterCut{"ConfEtaFilterCut", 0.8, "Eta cut for the global track"};                   // eta
    Configurable<bool> ConfDxaXYCustom0Cut{"ConfDxaXYCustom0Cut", false, "Enable Custom Dcaxy < [0] cut."};
    Configurable<float> ConfDcaXYFilterCut{"ConfDcaXYFilterCut", 2.4, "Value for DCA_XY for the global track"}; // max dca to vertex XY
    Configurable<float> ConfDcaZFilterCut{"ConfDcaZFilterCut", 3.2, "Value for DCA_Z for the global track"};    // max dca to vertex Z
    Configurable<bool> ConfDcaXYCustom1Cut{"ConfDcaXYCustom1Cut", true, "Enable Custom |DCAxy| < [1] + [2]/pt cut."};
    Configurable<float> ConfDcaXYCustom11FilterCut{"ConfDcaXY1FilterCut", 0.004, "Value for [1] custom DCAxy cut -> |DCAxy| < [1] + [2]/pT"};
    Configurable<float> ConfDcaXYCustom12FilterCut{"ConfDcaXY2FilterCut", 0.013, "Value for [2] custom DCAxy cut -> |DCAxy| < [1] + [2]/pT"};
  } ConfFilterCuts;

  Filter GlobalCutFilter = requireGlobalTrackInFilter();
  Filter CustomTrackFilter = (aod::track::pt > ConfFilterCuts.ConfPtLowFilterCut) &&
                             (aod::track::pt < ConfFilterCuts.ConfPtHighFilterCut) &&
                             (nabs(aod::track::eta) < ConfFilterCuts.ConfEtaFilterCut) &&
                             (!ConfFilterCuts.ConfDxaXYCustom0Cut || (aod::track::dcaXY < ConfFilterCuts.ConfDcaXYFilterCut)) && // true if configurable set to false or if configurable is true and it passes the selection
                             (aod::track::dcaZ < ConfFilterCuts.ConfDcaZFilterCut) &&
                             (!ConfFilterCuts.ConfDcaXYCustom1Cut || (nabs(aod::track::dcaXY) < ConfFilterCuts.ConfDcaXYCustom11FilterCut + ConfFilterCuts.ConfDcaXYCustom12FilterCut / aod::track::pt)); // same logic here

  // CASCADE
  FemtoUniverseCascadeSelection cascadeCuts;
  struct : o2::framework::ConfigurableGroup {
    Configurable<std::vector<float>> ConfCascSign{FemtoUniverseCascadeSelection::getSelectionName(femtoUniverseCascadeSelection::kCascadeSign, "ConfCasc"), std::vector<float>{-1, 1}, FemtoUniverseCascadeSelection::getSelectionHelper(femtoUniverseCascadeSelection::kCascadeSign, "Cascade selection: ")};
    Configurable<std::vector<float>> ConfCascPtMin{FemtoUniverseCascadeSelection::getSelectionName(femtoUniverseCascadeSelection::kCascadepTMin, "ConfCasc"), std::vector<float>{0.3f, 0.4f, 0.5f}, FemtoUniverseCascadeSelection::getSelectionHelper(femtoUniverseCascadeSelection::kCascadepTMin, "Cascade selection: ")};
    Configurable<std::vector<float>> ConfCascPtMax{FemtoUniverseCascadeSelection::getSelectionName(femtoUniverseCascadeSelection::kCascadepTMax, "ConfCasc"), std::vector<float>{3.3f, 3.4f, 3.5f}, FemtoUniverseCascadeSelection::getSelectionHelper(femtoUniverseCascadeSelection::kCascadepTMax, "Cascade selection: ")};
    Configurable<std::vector<float>> ConfCascEtaMax{FemtoUniverseCascadeSelection::getSelectionName(femtoUniverseCascadeSelection::kCascadeetaMax, "ConfCasc"), std::vector<float>{0.8f, 0.7f, 0.9f}, FemtoUniverseCascadeSelection::getSelectionHelper(femtoUniverseCascadeSelection::kCascadeetaMax, "Cascade selection: ")};
    Configurable<std::vector<float>> ConfCascV0DCADaughMax{FemtoUniverseCascadeSelection::getSelectionName(femtoUniverseCascadeSelection::kCascadeV0DCADaughMax, "ConfCasc"), std::vector<float>{1.f, 1.2f, 1.5f}, FemtoUniverseCascadeSelection::getSelectionHelper(femtoUniverseCascadeSelection::kCascadeV0DCADaughMax, "Cascade selection: ")};
    Configurable<std::vector<float>> ConfCascV0CPAMin{FemtoUniverseCascadeSelection::getSelectionName(femtoUniverseCascadeSelection::kCascadeV0CPAMin, "ConfCasc"), std::vector<float>{0.99f, 0.95f}, FemtoUniverseCascadeSelection::getSelectionHelper(femtoUniverseCascadeSelection::kCascadeV0CPAMin, "Cascade selection: ")};
    Configurable<std::vector<float>> ConfCascV0TranRadMin{FemtoUniverseCascadeSelection::getSelectionName(femtoUniverseCascadeSelection::kCascadeV0TranRadMin, "ConfCasc"), std::vector<float>{0.2f}, FemtoUniverseCascadeSelection::getSelectionHelper(femtoUniverseCascadeSelection::kCascadeV0TranRadMin, "Cascade selection: ")};
    Configurable<std::vector<float>> ConfCascV0TranRadMax{FemtoUniverseCascadeSelection::getSelectionName(femtoUniverseCascadeSelection::kCascadeV0TranRadMax, "ConfCasc"), std::vector<float>{100.f}, FemtoUniverseCascadeSelection::getSelectionHelper(femtoUniverseCascadeSelection::kCascadeV0TranRadMax, "Cascade selection: ")};
    Configurable<std::vector<float>> ConfCascV0DecVtxMax{FemtoUniverseCascadeSelection::getSelectionName(femtoUniverseCascadeSelection::kCascadeV0DecVtxMax, "ConfCasc"), std::vector<float>{100.f}, FemtoUniverseCascadeSelection::getSelectionHelper(femtoUniverseCascadeSelection::kCascadeV0DecVtxMax, "Cascade selection: ")};
    Configurable<std::vector<float>> ConfCascDCADaughMax{FemtoUniverseCascadeSelection::getSelectionName(femtoUniverseCascadeSelection::kCascadeDCADaughMax, "ConfCasc"), std::vector<float>{1.f, 1.2f, 1.5f}, FemtoUniverseCascadeSelection::getSelectionHelper(femtoUniverseCascadeSelection::kCascadeDCADaughMax, "Cascade selection: ")};
    Configurable<std::vector<float>> ConfCascCPAMin{FemtoUniverseCascadeSelection::getSelectionName(femtoUniverseCascadeSelection::kCascadeCPAMin, "ConfCasc"), std::vector<float>{0.99f, 0.95f}, FemtoUniverseCascadeSelection::getSelectionHelper(femtoUniverseCascadeSelection::kCascadeCPAMin, "Cascade selection: ")};
    Configurable<std::vector<float>> ConfCascTranRadMin{FemtoUniverseCascadeSelection::getSelectionName(femtoUniverseCascadeSelection::kCascadeTranRadMin, "ConfCasc"), std::vector<float>{0.2f, 0.5f}, FemtoUniverseCascadeSelection::getSelectionHelper(femtoUniverseCascadeSelection::kCascadeTranRadMin, "Cascade selection: ")};
    Configurable<std::vector<float>> ConfCascTranRadMax{FemtoUniverseCascadeSelection::getSelectionName(femtoUniverseCascadeSelection::kCascadeTranRadMax, "ConfCasc"), std::vector<float>{100.f}, FemtoUniverseCascadeSelection::getSelectionHelper(femtoUniverseCascadeSelection::kCascadeTranRadMax, "Cascade selection: ")};
    Configurable<std::vector<float>> ConfCascDecVtxMax{FemtoUniverseCascadeSelection::getSelectionName(femtoUniverseCascadeSelection::kCascadeDecVtxMax, "ConfCasc"), std::vector<float>{100.f}, FemtoUniverseCascadeSelection::getSelectionHelper(femtoUniverseCascadeSelection::kCascadeDecVtxMax, "Cascade selection: ")};

    Configurable<std::vector<float>> ConfCascDCAPosToPV{FemtoUniverseCascadeSelection::getSelectionName(femtoUniverseCascadeSelection::kCascadeDCAPosToPV, "ConfCasc"), std::vector<float>{0.1f}, FemtoUniverseCascadeSelection::getSelectionHelper(femtoUniverseCascadeSelection::kCascadeDCAPosToPV, "Cascade selection: ")};
    Configurable<std::vector<float>> ConfCascDCANegToPV{FemtoUniverseCascadeSelection::getSelectionName(femtoUniverseCascadeSelection::kCascadeDCANegToPV, "ConfCasc"), std::vector<float>{0.1f}, FemtoUniverseCascadeSelection::getSelectionHelper(femtoUniverseCascadeSelection::kCascadeDCANegToPV, "Cascade selection: ")};
    Configurable<std::vector<float>> ConfCascDCABachToPV{FemtoUniverseCascadeSelection::getSelectionName(femtoUniverseCascadeSelection::kCascadeDCABachToPV, "ConfCasc"), std::vector<float>{0.1f}, FemtoUniverseCascadeSelection::getSelectionHelper(femtoUniverseCascadeSelection::kCascadeDCABachToPV, "Cascade selection: ")};
    Configurable<std::vector<float>> ConfCascDCAV0ToPV{FemtoUniverseCascadeSelection::getSelectionName(femtoUniverseCascadeSelection::kCascadeDCAV0ToPV, "ConfCasc"), std::vector<float>{0.01f}, FemtoUniverseCascadeSelection::getSelectionHelper(femtoUniverseCascadeSelection::kCascadeDCAV0ToPV, "Cascade selection: ")};
    Configurable<std::vector<float>> ConfCascV0MassLowLimit{FemtoUniverseCascadeSelection::getSelectionName(femtoUniverseCascadeSelection::kCascadeV0MassMin, "ConfCasc"), std::vector<float>{1.05f}, FemtoUniverseCascadeSelection::getSelectionHelper(femtoUniverseCascadeSelection::kCascadeV0MassMin, "Cascade selection: ")};
    Configurable<std::vector<float>> ConfCascV0MassUpLimit{FemtoUniverseCascadeSelection::getSelectionName(femtoUniverseCascadeSelection::kCascadeV0MassMax, "ConfCasc"), std::vector<float>{1.30f}, FemtoUniverseCascadeSelection::getSelectionHelper(femtoUniverseCascadeSelection::kCascadeV0MassMax, "Cascade selection: ")};

    Configurable<std::vector<float>> ConfCascChildCharge{"ConfCascChildSign", std::vector<float>{-1, 1}, "Cascade Child sel: Charge"};
    Configurable<std::vector<float>> ConfCascChildEtaMax{"ConfCascChildEtaMax", std::vector<float>{0.8f}, "Cascade Child sel: max eta"};
    Configurable<std::vector<float>> ConfCascChildTPCnClsMin{"ConfCascChildTPCnClsMin", std::vector<float>{80.f, 70.f, 60.f}, "Cascade Child sel: Min. nCls TPC"};
    // Configurable<std::vector<float>> ConfCascChildDCAMin{"ConfCascChildDCAMin", std::vector<float>{0.05f, 0.06f}, "Cascade Child sel:  Max. DCA Daugh to PV (cm)"};
    Configurable<std::vector<float>> ConfCascChildPIDnSigmaMax{"ConfCascChildPIDnSigmaMax", std::vector<float>{3.f, 4.f}, "Cascade Child sel: Max. PID nSigma TPC"};
    Configurable<std::vector<int>> ConfCascChildPIDspecies{"ConfCascChildPIDspecies", std::vector<int>{o2::track::PID::Pion, o2::track::PID::Proton}, "Cascade Child sel: particle species for PID"};

    Configurable<float> ConfCascInvMassLowLimit{"ConfCascInvMassLowLimit", 1.25, "Lower limit of the cascade invariant mass"};
    Configurable<float> ConfCascInvMassUpLimit{"ConfCascInvMassUpLimit", 1.40, "Upper limit of the cascade invariant mass"};

    Configurable<bool> ConfCascRejectCompetingMass{"ConfCascRejectCompetingMass", false, "Switch on to reject Omegas (for Xi) or Xis (for Omegas)"};
    Configurable<float> ConfCascInvCompetingMassLowLimit{"ConfCascInvCompetingMassLowLimit", 1.66, "Lower limit of the cascade invariant mass for competing mass rejection"};
    Configurable<float> ConfCascInvCompetingMassUpLimit{"ConfCascInvCompetingMassUpLimit", 1.68, "Upper limit of the cascade invariant mass for competing mass rejection"};
  } ConfCascadeSelection;

  // PHI
  FemtoUniversePhiSelection phiCuts;
  struct : o2::framework::ConfigurableGroup {
    Configurable<bool> ConfLooseTPCNSigma{"ConfLooseTPCNSigma", false, "Use loose TPC N sigmas for Kaon PID."};
    Configurable<float> ConfLooseTPCNSigmaValue{"ConfLooseTPCNSigmaValue", 10, "Value for the loose TPC N Sigma for Kaon PID."};
    Configurable<bool> ConfLooseTOFNSigma{"ConfLooseTOFNSigma", false, "Use loose TPC N sigmas for Kaon PID."};
    Configurable<float> ConfLooseTOFNSigmaValue{"ConfLooseTOFNSigmaValue", 10, "Value for the loose TOF N Sigma for Kaon PID."};
    Configurable<float> ConfInvMassLowLimitPhi{"ConfInvMassLowLimitPhi", 1.011, "Lower limit of the Phi invariant mass"}; // change that to do invariant mass cut
    Configurable<float> ConfInvMassUpLimitPhi{"ConfInvMassUpLimitPhi", 1.027, "Upper limit of the Phi invariant mass"};
    Configurable<float> ConfPtLowLimitPhi{"ConfPtLowLimitPhi", 0.8, "Lower limit of the Phi pT."};
    Configurable<float> ConfPtHighLimitPhi{"ConfPtHighLimitPhi", 4.0, "Higher limit of the Phi pT."};
    Configurable<float> ConfNsigmaRejectPion{"ConfNsigmaRejectPion", 3.0, "Reject if particle could be a Pion combined nsigma value."};
    Configurable<float> ConfNsigmaRejectProton{"ConfNsigmaRejectProton", 3.0, "Reject if particle could be a Proton combined nsigma value."};
  } ConfPhiSelection;

  Configurable<int> ConfPDGCodePartOne{"ConfPDGCodePartOne", 321, "Particle 1 - PDG code"};
  Configurable<int> ConfPDGCodePartTwo{"ConfPDGCodePartTwo", 321, "Particle 2 - PDG code"};

  // D0/D0bar mesons
  struct : o2::framework::ConfigurableGroup {
    Configurable<float> ConfD0D0barCandMaxY{"ConfD0D0barCandMaxY", -1., "max. cand. rapidity"};
    Configurable<float> ConfD0D0barCandEtaCut{"ConfD0D0barCandEtaCut", 0.8, "max. cand. pseudorapidity"};
  } ConfD0Selection;

  HfHelper hfHelper;

  bool IsKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK)
  {

    if (mom < 0.3) { // 0.0-0.3
      if (TMath::Abs(nsigmaTPCK) < 3.0) {
        return true;
      } else {
        return false;
      }
    } else if (mom < 0.45) { // 0.30 - 0.45
      if (TMath::Abs(nsigmaTPCK) < 2.0) {
        return true;
      } else {
        return false;
      }
    } else if (mom < 0.55) { // 0.45-0.55
      if (TMath::Abs(nsigmaTPCK) < 1.0) {
        return true;
      } else {
        return false;
      }
    } else if (mom < 1.5) { // 0.55-1.5 (now we use TPC and TOF)
      if ((TMath::Abs(nsigmaTOFK) < 3.0) && (TMath::Abs(nsigmaTPCK) < 3.0)) {
        {
          return true;
        }
      } else {
        return false;
      }
    } else if (mom > 1.5) { // 1.5 -
      if ((TMath::Abs(nsigmaTOFK) < 2.0) && (TMath::Abs(nsigmaTPCK) < 3.0)) {
        return true;
      } else {
        return false;
      }
    } else {
      return false;
    }
  }

  bool IsKaonNSigmaTPCLoose(float mom, float nsigmaTPCK, float nsigmaTOFK)
  {

    if (mom < 0.3) { // 0.0-0.3
      if (TMath::Abs(nsigmaTPCK) < ConfPhiSelection.ConfLooseTPCNSigmaValue.value) {
        return true;
      } else {
        return false;
      }
    } else if (mom < 0.45) { // 0.30 - 0.45
      if (TMath::Abs(nsigmaTPCK) < ConfPhiSelection.ConfLooseTPCNSigmaValue.value) {
        return true;
      } else {
        return false;
      }
    } else if (mom < 0.55) { // 0.45-0.55
      if (TMath::Abs(nsigmaTPCK) < ConfPhiSelection.ConfLooseTPCNSigmaValue.value) {
        return true;
      } else {
        return false;
      }
    } else if (mom < 1.5) { // 0.55-1.5 (now we use TPC and TOF)
      if ((TMath::Abs(nsigmaTOFK) < 3.0) && (TMath::Abs(nsigmaTPCK) < ConfPhiSelection.ConfLooseTPCNSigmaValue.value)) {
        return true;
      } else {
        return false;
      }
    } else if (mom > 1.5) { // 1.5 -
      if ((TMath::Abs(nsigmaTOFK) < 2.0) && (TMath::Abs(nsigmaTPCK) < ConfPhiSelection.ConfLooseTPCNSigmaValue.value)) {
        return true;
      } else {
        return false;
      }
    } else {
      return false;
    }
  }

  bool IsKaonNSigmaTOFLoose(float mom, float nsigmaTPCK, float nsigmaTOFK)
  {
    if (mom < 0.3) { // 0.0-0.3
      if (TMath::Abs(nsigmaTPCK) < 3.0) {
        return true;
      } else {
        return false;
      }
    } else if (mom < 0.45) { // 0.30 - 0.45
      if (TMath::Abs(nsigmaTPCK) < 2.0) {
        return true;
      } else {
        return false;
      }
    } else if (mom < 0.55) { // 0.45-0.55
      if (TMath::Abs(nsigmaTPCK) < 1.0) {
        return true;
      } else {
        return false;
      }
    } else if (mom < 1.5) { // 0.55-1.5 (now we use TPC and TOF)
      if ((TMath::Abs(nsigmaTOFK) < ConfPhiSelection.ConfLooseTOFNSigmaValue.value) && (TMath::Abs(nsigmaTPCK) < 3.0)) {
        {
          return true;
        }
      } else {
        return false;
      }
    } else if (mom > 1.5) { // 1.5 -
      if ((TMath::Abs(nsigmaTOFK) < ConfPhiSelection.ConfLooseTOFNSigmaValue.value) && (TMath::Abs(nsigmaTPCK) < 3.0)) {
        return true;
      } else {
        return false;
      }
    } else {
      return false;
    }
  }

  bool IsKaonRejected(float mom, float nsigmaTPCPr, float nsigmaTOFPr, float nsigmaTPCPi, float nsigmaTOFPi)
  {
    if (mom < 0.5) {
      if (TMath::Abs(nsigmaTPCPi) < ConfPhiSelection.ConfNsigmaRejectPion.value) {
        return true;
      } else if (TMath::Abs(nsigmaTPCPr) < ConfPhiSelection.ConfNsigmaRejectProton.value) {
        return true;
      }
    }
    if (mom > 0.5) {
      if (TMath::Hypot(nsigmaTOFPi, nsigmaTPCPi) < ConfPhiSelection.ConfNsigmaRejectPion.value) {
        return true;
      } else if (TMath::Hypot(nsigmaTOFPr, nsigmaTPCPr) < ConfPhiSelection.ConfNsigmaRejectProton.value) {
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

  int mRunNumber;
  float mMagField;
  Service<o2::ccdb::BasicCCDBManager> ccdb; /// Accessing the CCDB

  void init(InitContext&)
  {
    if ((doprocessFullData || doprocessTrackPhiData || doprocessTrackData || doprocessTrackV0 || doprocessTrackCascadeData || doprocessTrackD0mesonData || doprocessTrackCentRun2Data || doprocessTrackCentRun3Data || doprocessTrackV0CentRun3) == false && (doprocessFullMC || doprocessTrackMC || doprocessTrackMCTruth || doprocessTrackMCGen || doprocessTruthAndFullMC || doprocessFullMCCent) == false) {
      LOGF(fatal, "Neither processFullData nor processFullMC enabled. Please choose one.");
    }
    if ((doprocessFullData || doprocessTrackPhiData || doprocessTrackData || doprocessTrackV0 || doprocessTrackCascadeData || doprocessTrackD0mesonData || doprocessTrackCentRun2Data || doprocessTrackCentRun3Data || doprocessTrackV0CentRun3) == true && (doprocessFullMC || doprocessTrackMC || doprocessTrackMCTruth || doprocessTrackMCGen || doprocessTruthAndFullMC || doprocessFullMCCent) == true) {
      LOGF(fatal,
           "Cannot enable process Data and process MC at the same time. "
           "Please choose one.");
    }

    colCuts.setCuts(ConfEvtZvtx, ConfEvtTriggerCheck, ConfEvtTriggerSel, ConfEvtOfflineCheck, ConfIsRun3, ConfCentFT0Min, ConfCentFT0Max);
    colCuts.init(&qaRegistry);

    trackCuts.setSelection(ConfTrkCharge, femtoUniverseTrackSelection::kSign, femtoUniverseSelection::kEqual);
    trackCuts.setSelection(ConfTrkPtmin, femtoUniverseTrackSelection::kpTMin, femtoUniverseSelection::kLowerLimit);
    trackCuts.setSelection(ConfTrkPtmax, femtoUniverseTrackSelection::kpTMax, femtoUniverseSelection::kUpperLimit);
    trackCuts.setSelection(ConfTrkEta, femtoUniverseTrackSelection::kEtaMax, femtoUniverseSelection::kAbsUpperLimit);
    trackCuts.setSelection(ConfTrkTPCnclsMin, femtoUniverseTrackSelection::kTPCnClsMin, femtoUniverseSelection::kLowerLimit);
    trackCuts.setSelection(ConfTrkTPCfCls, femtoUniverseTrackSelection::kTPCfClsMin, femtoUniverseSelection::kLowerLimit);
    trackCuts.setSelection(ConfTrkTPCcRowsMin, femtoUniverseTrackSelection::kTPCcRowsMin, femtoUniverseSelection::kLowerLimit);
    trackCuts.setSelection(ConfTrkTPCsCls, femtoUniverseTrackSelection::kTPCsClsMax, femtoUniverseSelection::kUpperLimit);
    trackCuts.setSelection(ConfTrkITSnclsMin, femtoUniverseTrackSelection::kITSnClsMin, femtoUniverseSelection::kLowerLimit);
    trackCuts.setSelection(ConfTrkITSnclsIbMin, femtoUniverseTrackSelection::kITSnClsIbMin, femtoUniverseSelection::kLowerLimit);
    trackCuts.setSelection(ConfTrkDCAxyMax, femtoUniverseTrackSelection::kDCAxyMax, femtoUniverseSelection::kAbsUpperLimit);
    trackCuts.setSelection(ConfTrkDCAzMax, femtoUniverseTrackSelection::kDCAzMax, femtoUniverseSelection::kAbsUpperLimit);
    trackCuts.setSelection(ConfTrkPIDnSigmaMax, femtoUniverseTrackSelection::kPIDnSigmaMax, femtoUniverseSelection::kAbsUpperLimit);
    trackCuts.setPIDSpecies(ConfTrkPIDspecies);
    trackCuts.setnSigmaPIDOffset(ConfTrkPIDnSigmaOffsetTPC, ConfTrkPIDnSigmaOffsetTOF);
    trackCuts.init<aod::femtouniverseparticle::ParticleType::kTrack, aod::femtouniverseparticle::TrackType::kNoChild, aod::femtouniverseparticle::cutContainerType>(&qaRegistry);

    /// \todo fix how to pass array to setSelection, getRow() passing a
    /// different type!
    // v0Cuts.setSelection(ConfV0Selection->getRow(0),
    // femtoUniverseV0Selection::kDecVtxMax, femtoUniverseSelection::kAbsUpperLimit);
    if (ConfIsActivateV0) {
      // initializing for V0
      v0Cuts.setSelection(ConfV0Selection.ConfV0Sign, femtoUniverseV0Selection::kV0Sign, femtoUniverseSelection::kEqual);
      v0Cuts.setSelection(ConfV0Selection.ConfV0PtMin, femtoUniverseV0Selection::kV0pTMin, femtoUniverseSelection::kLowerLimit);
      v0Cuts.setSelection(ConfV0Selection.ConfV0PtMax, femtoUniverseV0Selection::kV0pTMax, femtoUniverseSelection::kUpperLimit);
      v0Cuts.setSelection(ConfV0Selection.ConfV0EtaMax, femtoUniverseV0Selection::kV0etaMax, femtoUniverseSelection::kAbsUpperLimit);
      v0Cuts.setSelection(ConfV0Selection.ConfV0DCADaughMax, femtoUniverseV0Selection::kV0DCADaughMax, femtoUniverseSelection::kUpperLimit);
      v0Cuts.setSelection(ConfV0Selection.ConfV0CPAMin, femtoUniverseV0Selection::kV0CPAMin, femtoUniverseSelection::kLowerLimit);
      v0Cuts.setSelection(ConfV0Selection.ConfV0TranRadMin, femtoUniverseV0Selection::kV0TranRadMin, femtoUniverseSelection::kLowerLimit);
      v0Cuts.setSelection(ConfV0Selection.ConfV0TranRadMax, femtoUniverseV0Selection::kV0TranRadMax, femtoUniverseSelection::kUpperLimit);
      v0Cuts.setSelection(ConfV0Selection.ConfV0DecVtxMax, femtoUniverseV0Selection::kV0DecVtxMax, femtoUniverseSelection::kUpperLimit);
      v0Cuts.setChildCuts(femtoUniverseV0Selection::kPosTrack, ConfV0Selection.ConfChildCharge, femtoUniverseTrackSelection::kSign, femtoUniverseSelection::kEqual);
      v0Cuts.setChildCuts(femtoUniverseV0Selection::kPosTrack, ConfV0Selection.ConfChildEtaMax, femtoUniverseTrackSelection::kEtaMax, femtoUniverseSelection::kAbsUpperLimit);
      v0Cuts.setChildCuts(femtoUniverseV0Selection::kPosTrack, ConfV0Selection.ConfChildTPCnClsMin, femtoUniverseTrackSelection::kTPCnClsMin, femtoUniverseSelection::kLowerLimit);
      v0Cuts.setChildCuts(femtoUniverseV0Selection::kPosTrack, ConfV0Selection.ConfChildDCAMin, femtoUniverseTrackSelection::kDCAMin, femtoUniverseSelection::kAbsLowerLimit);
      v0Cuts.setChildCuts(femtoUniverseV0Selection::kPosTrack, ConfV0Selection.ConfChildPIDnSigmaMax, femtoUniverseTrackSelection::kPIDnSigmaMax, femtoUniverseSelection::kAbsUpperLimit);
      v0Cuts.setChildCuts(femtoUniverseV0Selection::kNegTrack, ConfV0Selection.ConfChildCharge, femtoUniverseTrackSelection::kSign, femtoUniverseSelection::kEqual);
      v0Cuts.setChildCuts(femtoUniverseV0Selection::kNegTrack, ConfV0Selection.ConfChildEtaMax, femtoUniverseTrackSelection::kEtaMax, femtoUniverseSelection::kAbsUpperLimit);
      v0Cuts.setChildCuts(femtoUniverseV0Selection::kNegTrack, ConfV0Selection.ConfChildTPCnClsMin, femtoUniverseTrackSelection::kTPCnClsMin, femtoUniverseSelection::kLowerLimit);
      v0Cuts.setChildCuts(femtoUniverseV0Selection::kNegTrack, ConfV0Selection.ConfChildDCAMin, femtoUniverseTrackSelection::kDCAMin, femtoUniverseSelection::kAbsLowerLimit);
      v0Cuts.setChildCuts(femtoUniverseV0Selection::kNegTrack, ConfV0Selection.ConfChildPIDnSigmaMax, femtoUniverseTrackSelection::kPIDnSigmaMax, femtoUniverseSelection::kAbsUpperLimit);
      v0Cuts.setChildPIDSpecies(femtoUniverseV0Selection::kPosTrack, ConfV0Selection.ConfChildPIDspecies);
      v0Cuts.setChildPIDSpecies(femtoUniverseV0Selection::kNegTrack, ConfV0Selection.ConfChildPIDspecies);
      v0Cuts.init<aod::femtouniverseparticle::ParticleType::kV0, aod::femtouniverseparticle::ParticleType::kV0Child, aod::femtouniverseparticle::cutContainerType>(&qaRegistry);
      v0Cuts.setInvMassLimits(ConfV0Selection.ConfV0InvMassLowLimit, ConfV0Selection.ConfV0InvMassUpLimit);

      v0Cuts.setChildRejectNotPropagatedTracks(femtoUniverseV0Selection::kPosTrack, ConfTrkRejectNotPropagated);
      v0Cuts.setChildRejectNotPropagatedTracks(femtoUniverseV0Selection::kNegTrack, ConfTrkRejectNotPropagated);

      v0Cuts.setnSigmaPIDOffsetTPC(ConfTrkPIDnSigmaOffsetTPC);
      v0Cuts.setChildnSigmaPIDOffset(femtoUniverseV0Selection::kPosTrack, ConfTrkPIDnSigmaOffsetTPC, ConfTrkPIDnSigmaOffsetTOF);
      v0Cuts.setChildnSigmaPIDOffset(femtoUniverseV0Selection::kNegTrack, ConfTrkPIDnSigmaOffsetTPC, ConfTrkPIDnSigmaOffsetTOF);

      if (ConfV0Selection.ConfV0RejectKaons) {
        v0Cuts.setKaonInvMassLimits(ConfV0Selection.ConfV0InvKaonMassLowLimit, ConfV0Selection.ConfV0InvKaonMassUpLimit);
      }
      // if (ConfRejectITSHitandTOFMissing) {
      //   o2PhysicsTrackSelection = new
      //   TrackSelection(getGlobalTrackSelection());
      //   o2PhysicsTrackSelection->SetRequireHitsInITSLayers(1, {0, 1, 2, 3});
      // }
    }

    if (ConfIsActivateCascade) {
      // initializing for cascades
      cascadeCuts.setSelection(ConfCascadeSelection.ConfCascSign, femtoUniverseCascadeSelection::kCascadeSign, femtoUniverseSelection::kEqual);
      cascadeCuts.setSelection(ConfCascadeSelection.ConfCascPtMin, femtoUniverseCascadeSelection::kCascadepTMin, femtoUniverseSelection::kLowerLimit);
      cascadeCuts.setSelection(ConfCascadeSelection.ConfCascPtMax, femtoUniverseCascadeSelection::kCascadepTMax, femtoUniverseSelection::kUpperLimit);
      cascadeCuts.setSelection(ConfCascadeSelection.ConfCascEtaMax, femtoUniverseCascadeSelection::kCascadeetaMax, femtoUniverseSelection::kAbsUpperLimit);
      // v0 child cuts
      cascadeCuts.setSelection(ConfCascadeSelection.ConfCascV0DCADaughMax, femtoUniverseCascadeSelection::kCascadeV0DCADaughMax, femtoUniverseSelection::kUpperLimit);
      cascadeCuts.setSelection(ConfCascadeSelection.ConfCascV0CPAMin, femtoUniverseCascadeSelection::kCascadeV0CPAMin, femtoUniverseSelection::kLowerLimit);
      cascadeCuts.setSelection(ConfCascadeSelection.ConfCascV0TranRadMin, femtoUniverseCascadeSelection::kCascadeV0TranRadMin, femtoUniverseSelection::kLowerLimit);
      cascadeCuts.setSelection(ConfCascadeSelection.ConfCascV0TranRadMax, femtoUniverseCascadeSelection::kCascadeV0TranRadMax, femtoUniverseSelection::kUpperLimit);
      cascadeCuts.setSelection(ConfCascadeSelection.ConfCascV0DecVtxMax, femtoUniverseCascadeSelection::kCascadeV0DecVtxMax, femtoUniverseSelection::kUpperLimit);
      // cascade cuts
      cascadeCuts.setSelection(ConfCascadeSelection.ConfCascDCADaughMax, femtoUniverseCascadeSelection::kCascadeDCADaughMax, femtoUniverseSelection::kUpperLimit);
      cascadeCuts.setSelection(ConfCascadeSelection.ConfCascCPAMin, femtoUniverseCascadeSelection::kCascadeCPAMin, femtoUniverseSelection::kLowerLimit);
      cascadeCuts.setSelection(ConfCascadeSelection.ConfCascTranRadMin, femtoUniverseCascadeSelection::kCascadeTranRadMin, femtoUniverseSelection::kLowerLimit);
      cascadeCuts.setSelection(ConfCascadeSelection.ConfCascTranRadMax, femtoUniverseCascadeSelection::kCascadeTranRadMax, femtoUniverseSelection::kUpperLimit);
      cascadeCuts.setSelection(ConfCascadeSelection.ConfCascDecVtxMax, femtoUniverseCascadeSelection::kCascadeDecVtxMax, femtoUniverseSelection::kUpperLimit);
      cascadeCuts.setSelection(ConfCascadeSelection.ConfCascDCAPosToPV, femtoUniverseCascadeSelection::kCascadeDCAPosToPV, femtoUniverseSelection::kLowerLimit);
      cascadeCuts.setSelection(ConfCascadeSelection.ConfCascDCANegToPV, femtoUniverseCascadeSelection::kCascadeDCANegToPV, femtoUniverseSelection::kLowerLimit);
      cascadeCuts.setSelection(ConfCascadeSelection.ConfCascDCABachToPV, femtoUniverseCascadeSelection::kCascadeDCABachToPV, femtoUniverseSelection::kLowerLimit);
      cascadeCuts.setSelection(ConfCascadeSelection.ConfCascDCAV0ToPV, femtoUniverseCascadeSelection::kCascadeDCAV0ToPV, femtoUniverseSelection::kLowerLimit);
      cascadeCuts.setSelection(ConfCascadeSelection.ConfCascV0MassLowLimit, femtoUniverseCascadeSelection::kCascadeV0MassMin, femtoUniverseSelection::kLowerLimit);
      cascadeCuts.setSelection(ConfCascadeSelection.ConfCascV0MassUpLimit, femtoUniverseCascadeSelection::kCascadeV0MassMax, femtoUniverseSelection::kUpperLimit);
      // children cuts
      cascadeCuts.setChildCuts(femtoUniverseCascadeSelection::kPosTrack, ConfCascadeSelection.ConfCascChildCharge, femtoUniverseTrackSelection::kSign, femtoUniverseSelection::kEqual);
      cascadeCuts.setChildCuts(femtoUniverseCascadeSelection::kPosTrack, ConfCascadeSelection.ConfCascChildEtaMax, femtoUniverseTrackSelection::kEtaMax, femtoUniverseSelection::kAbsUpperLimit);
      cascadeCuts.setChildCuts(femtoUniverseCascadeSelection::kPosTrack, ConfCascadeSelection.ConfCascChildTPCnClsMin, femtoUniverseTrackSelection::kTPCnClsMin, femtoUniverseSelection::kLowerLimit);
      cascadeCuts.setChildCuts(femtoUniverseCascadeSelection::kPosTrack, ConfCascadeSelection.ConfCascChildPIDnSigmaMax, femtoUniverseTrackSelection::kPIDnSigmaMax, femtoUniverseSelection::kAbsUpperLimit);
      cascadeCuts.setChildCuts(femtoUniverseCascadeSelection::kNegTrack, ConfCascadeSelection.ConfCascChildCharge, femtoUniverseTrackSelection::kSign, femtoUniverseSelection::kEqual);
      cascadeCuts.setChildCuts(femtoUniverseCascadeSelection::kNegTrack, ConfCascadeSelection.ConfCascChildEtaMax, femtoUniverseTrackSelection::kEtaMax, femtoUniverseSelection::kAbsUpperLimit);
      cascadeCuts.setChildCuts(femtoUniverseCascadeSelection::kNegTrack, ConfCascadeSelection.ConfCascChildTPCnClsMin, femtoUniverseTrackSelection::kTPCnClsMin, femtoUniverseSelection::kLowerLimit);
      cascadeCuts.setChildCuts(femtoUniverseCascadeSelection::kNegTrack, ConfCascadeSelection.ConfCascChildPIDnSigmaMax, femtoUniverseTrackSelection::kPIDnSigmaMax, femtoUniverseSelection::kAbsUpperLimit);
      cascadeCuts.setChildCuts(femtoUniverseCascadeSelection::kBachTrack, ConfCascadeSelection.ConfCascChildCharge, femtoUniverseTrackSelection::kSign, femtoUniverseSelection::kEqual);
      cascadeCuts.setChildCuts(femtoUniverseCascadeSelection::kBachTrack, ConfCascadeSelection.ConfCascChildEtaMax, femtoUniverseTrackSelection::kEtaMax, femtoUniverseSelection::kAbsUpperLimit);
      cascadeCuts.setChildCuts(femtoUniverseCascadeSelection::kBachTrack, ConfCascadeSelection.ConfCascChildTPCnClsMin, femtoUniverseTrackSelection::kTPCnClsMin, femtoUniverseSelection::kLowerLimit);
      cascadeCuts.setChildCuts(femtoUniverseCascadeSelection::kBachTrack, ConfCascadeSelection.ConfCascChildPIDnSigmaMax, femtoUniverseTrackSelection::kPIDnSigmaMax, femtoUniverseSelection::kAbsUpperLimit);

      // TODO
      cascadeCuts.setChildPIDSpecies(femtoUniverseCascadeSelection::kPosTrack, ConfV0Selection.ConfChildPIDspecies);
      cascadeCuts.setChildPIDSpecies(femtoUniverseCascadeSelection::kNegTrack, ConfV0Selection.ConfChildPIDspecies);
      cascadeCuts.setChildPIDSpecies(femtoUniverseCascadeSelection::kBachTrack, ConfCascadeSelection.ConfCascChildPIDspecies);

      // check if works correctly for bachelor track
      cascadeCuts.init<aod::femtouniverseparticle::ParticleType::kCascade, aod::femtouniverseparticle::ParticleType::kV0Child, aod::femtouniverseparticle::ParticleType::kCascadeBachelor, aod::femtouniverseparticle::cutContainerType>(&cascadeQaRegistry, ConfIsSelectCascOmega);
      // invmass cuts
      cascadeCuts.setInvMassLimits(ConfCascadeSelection.ConfCascInvMassLowLimit, ConfCascadeSelection.ConfCascInvMassUpLimit);

      if (ConfCascadeSelection.ConfCascRejectCompetingMass) {
        cascadeCuts.setCompetingInvMassLimits(ConfCascadeSelection.ConfCascInvCompetingMassLowLimit, ConfCascadeSelection.ConfCascInvCompetingMassUpLimit);
      }
    }

    if (ConfIsActivatePhi) {
      // initializing for Phi meson
      phiCuts.init<aod::femtouniverseparticle::ParticleType::kPhi, aod::femtouniverseparticle::ParticleType::kPhiChild, aod::femtouniverseparticle::cutContainerType>(&qaRegistry);
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

    if (ConfIsRun3 && !ConfIsForceGRP) {
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

  template <bool isTrackOrV0, bool isPhiOrD0, bool isXi, typename ParticleType>
  void fillDebugParticle(ParticleType const& particle)
  {
    if constexpr (isTrackOrV0) {
      outputDebugParts(particle.sign(), (uint8_t)particle.tpcNClsFound(),
                       particle.tpcNClsFindable(),
                       (uint8_t)particle.tpcNClsCrossedRows(),
                       particle.tpcNClsShared(), particle.tpcInnerParam(),
                       particle.itsNCls(), particle.itsNClsInnerBarrel(),
                       particle.dcaXY(), particle.dcaZ(), particle.tpcSignal(),
                       particle.tpcNSigmaStoreEl(), particle.tpcNSigmaStorePi(),
                       particle.tpcNSigmaStoreKa(), particle.tpcNSigmaStorePr(),
                       particle.tpcNSigmaStoreDe(), particle.tofNSigmaStoreEl(),
                       particle.tofNSigmaStorePi(), particle.tofNSigmaStoreKa(),
                       particle.tofNSigmaStorePr(), particle.tofNSigmaStoreDe(),
                       -999., -999., -999., -999., -999., -999.);
    } else if constexpr (isPhiOrD0) {
      outputDebugParts(-999., -999., -999., -999., -999., -999., -999., -999.,
                       -999., -999., -999., -999., -999., -999., -999., -999.,
                       -999., -999., -999., -999., -999.,
                       -999., -999.,
                       -999., -999., -999.,
                       -999.); // QA for phi or D0/D0bar
    } else if constexpr (isXi) {
      outputDebugParts(-999., -999., -999., -999., -999., -999., -999., -999.,
                       -999., -999., -999., -999., -999., -999., -999., -999.,
                       -999., -999., -999., -999., -999.,
                       particle.dcacascdaughters(), particle.cascradius(),
                       particle.x(), particle.y(), particle.z(),
                       particle.mOmega()); // QA for Xi Cascades (later do the same for Omegas)
    } else {
      // LOGF(info, "isTrack0orV0: %d, isPhi: %d", isTrackOrV0, isPhiOrD0);
      outputDebugParts(-999., -999., -999., -999., -999., -999., -999., -999.,
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

      if (abs(pdgCode) == abs(ConfPDGCodePartOne.value) || abs(pdgCode) == abs(ConfPDGCodePartTwo.value)) {
        if (particleMC.isPhysicalPrimary()) {
          particleOrigin = aod::femtouniverseMCparticle::ParticleOriginMCTruth::kPrimary;
        } else if (!motherparticlesMC.empty()) {
          auto motherparticleMC = motherparticlesMC.front();
          if (motherparticleMC.producedByGenerator())
            particleOrigin = checkDaughterType(fdparttype, motherparticleMC.pdgCode());
        } else {
          particleOrigin = aod::femtouniverseMCparticle::ParticleOriginMCTruth::kMaterial;
        }
      } else {
        particleOrigin = aod::femtouniverseMCparticle::ParticleOriginMCTruth::kFake;
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

      if (abs(pdgCode1) == abs(321) || abs(pdgCode2) == abs(-321)) {
        if ((kaon1MC.isPhysicalPrimary() && kaon2MC.isPhysicalPrimary()) && (!motherskaon1MC.empty() && !motherskaon2MC.empty())) {
          for (auto& particleMotherOfNeg : motherskaon1MC) {
            for (auto& particleMotherOfPos : motherskaon2MC) {
              if (particleMotherOfNeg == particleMotherOfPos && particleMotherOfNeg.pdgCode() == 333) {
                phiOrigin = aod::femtouniverseMCparticle::ParticleOriginMCTruth::kPrimary;
              } else {
                phiOrigin = aod::femtouniverseMCparticle::ParticleOriginMCTruth::kFake;
              }
            }
          }
        } else {
          phiOrigin = aod::femtouniverseMCparticle::ParticleOriginMCTruth::kFake;
        }
      } else {
        phiOrigin = aod::femtouniverseMCparticle::ParticleOriginMCTruth::kFake;
      }

      TLorentzVector part1Vec;
      TLorentzVector part2Vec;

      float mMassOne = TDatabasePDG::Instance()->GetParticle(321)->Mass();  // FIXME: Get from the PDG service of the common header
      float mMassTwo = TDatabasePDG::Instance()->GetParticle(-321)->Mass(); // FIXME: Get from the PDG service of the common header

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
  void fillCollisions(CollisionType const& col, TrackType const& tracks)
  {
    const auto vtxZ = col.posZ();
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
    if (ConfEvtUseTPCmult) {
      multNtr = col.multTPC();
    }

    // check whether the basic event selection criteria are fulfilled
    // if the basic selection is NOT fulfilled:
    // in case of skimming run - don't store such collisions
    // in case of trigger run - store such collisions but don't store any
    // particle candidates for such collisions
    if (!colCuts.isSelected(col)) {
      return;
    }

    if (!ConfIsUsePileUp) {
      if (ConfDoSpher) {
        outputCollision(vtxZ, mult, multNtr, colCuts.computeSphericity(col, tracks), mMagField);
      } else {
        outputCollision(vtxZ, mult, multNtr, 2, mMagField);
      }
    } else {
      if (ConfDoSpher && (!ConfEvNoSameBunchPileup || col.selection_bit(aod::evsel::kNoSameBunchPileup)) && (!ConfEvIsGoodZvtxFT0vsPV || col.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) && (!ConfEvIsVertexITSTPC || col.selection_bit(aod::evsel::kIsVertexITSTPC))) {
        outputCollision(vtxZ, mult, multNtr, colCuts.computeSphericity(col, tracks), mMagField);
      } else {
        outputCollision(vtxZ, mult, multNtr, 2, mMagField);
      }
    }
    colCuts.fillQA(col);
  }

  template <typename CollisionType, typename TrackType>
  void fillMCTruthCollisions(CollisionType const& col, TrackType const& tracks)
  {
    for (auto& c : col) {
      const auto vtxZ = c.posZ();
      int mult = 0;
      int multNtr = 0;

      if (std::abs(vtxZ) > ConfEvtZvtx) {
        continue;
      }

      if (ConfDoSpher) {
        outputCollision(vtxZ, mult, multNtr, colCuts.computeSphericity(col, tracks), mMagField);
      } else {
        outputCollision(vtxZ, mult, multNtr, 2, mMagField);
      }
    }
  }

  template <bool isMC, typename CollisionType, typename TrackType>
  void fillCollisionsCentRun2(CollisionType const& col, TrackType const& tracks)
  {
    const auto vtxZ = col.posZ();
    int cent = 0;
    int multNtr = 0;
    if (!ConfIsRun3) {
      cent = col.centRun2V0M();
      multNtr = col.multNTracksPV();
    }

    // check whether the basic event selection criteria are fulfilled
    // if the basic selection is NOT fulfilled:
    // in case of skimming run - don't store such collisions
    // in case of trigger run - store such collisions but don't store any
    // particle candidates for such collisions
    if (!colCuts.isSelected(col)) {
      return;
    }

    // colCuts.fillQA(col); //for now, TODO: create a configurable so in the FemroUniverseCollisionSelection.h there is an option to plot QA just for the posZ
    if (ConfDoSpher) {
      outputCollision(vtxZ, cent, multNtr, colCuts.computeSphericity(col, tracks), mMagField);
    } else {
      outputCollision(vtxZ, cent, multNtr, 2, mMagField);
    }
  }

  template <bool isMC, typename CollisionType, typename TrackType>
  void fillCollisionsCentRun3(CollisionType const& col, TrackType const& tracks)
  {
    const auto vtxZ = col.posZ();
    int cent = 0;
    int multNtr = 0;
    if (ConfIsRun3) {
      multNtr = col.multNTracksPV();
      cent = col.centFT0C();
    }

    int occupancy = col.trackOccupancyInTimeRange();
    // check whether the basic event selection criteria are fulfilled
    // if the basic selection is NOT fulfilled:
    // in case of skimming run - don't store such collisions
    // in case of trigger run - store such collisions but don't store any
    // particle candidates for such collisions
    if (!colCuts.isSelectedRun3(col)) {
      return;
    }

    // colCuts.fillQA(col); //for now, TODO: create a configurable so in the FemroUniverseCollisionSelection.h there is an option to plot QA just for the posZ
    if ((col.selection_bit(aod::evsel::kNoSameBunchPileup)) && (col.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) && (occupancy > ConfTPCOccupancyMin && occupancy <= ConfTPCOccupancyMax)) {
      if (ConfDoSpher) {
        outputCollision(vtxZ, cent, multNtr, colCuts.computeSphericity(col, tracks), mMagField);
      } else {
        outputCollision(vtxZ, cent, multNtr, 2, mMagField);
      }
    }
  }

  template <bool isMC, typename TrackType>
  void fillTracks(TrackType const& tracks)
  {
    std::vector<int> childIDs = {0, 0}; // these IDs are necessary to keep track of the children
    std::vector<int> tmpIDtrack;        // this vector keeps track of the matching of the primary track table row <-> aod::track table global index

    for (auto& track : tracks) {
      /// if the most open selection criteria are not fulfilled there is no
      /// point looking further at the track
      if (!trackCuts.isSelectedMinimal(track)) {
        continue;
      }

      if (track.pt() > ConfTOFpTmin) {
        if (!track.hasTOF()) {
          continue;
        }
      }

      trackCuts.fillQA<aod::femtouniverseparticle::ParticleType::kTrack,
                       aod::femtouniverseparticle::TrackType::kNoChild>(track);
      // the bit-wise container of the systematic variations is obtained
      auto cutContainer = trackCuts.getCutContainer<aod::femtouniverseparticle::cutContainerType>(track);

      // now the table is filled
      if (!ConfIsActivateCascade) {
        outputParts(outputCollision.lastIndex(), track.pt(), track.eta(),
                    track.phi(), aod::femtouniverseparticle::ParticleType::kTrack,
                    cutContainer.at(
                      femtoUniverseTrackSelection::TrackContainerPosition::kCuts),
                    cutContainer.at(
                      femtoUniverseTrackSelection::TrackContainerPosition::kPID),
                    track.dcaXY(), childIDs, 0, 0);
      } else {
        outputCascParts(outputCollision.lastIndex(), track.pt(), track.eta(),
                        track.phi(), aod::femtouniverseparticle::ParticleType::kTrack,
                        cutContainer.at(
                          femtoUniverseTrackSelection::TrackContainerPosition::kCuts),
                        cutContainer.at(
                          femtoUniverseTrackSelection::TrackContainerPosition::kPID),
                        track.dcaXY(), childIDs, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
      }
      tmpIDtrack.push_back(track.globalIndex());
      if (ConfIsDebug) {
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
    for (auto& v0 : fullV0s) {
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
      auto cutContainerV0 = v0Cuts.getCutContainer<aod::femtouniverseparticle::cutContainerType>(col, v0, postrack, negtrack);

      int postrackID = v0.posTrackId();
      int rowInPrimaryTrackTablePos = -1;
      rowInPrimaryTrackTablePos = getRowDaughters(postrackID, tmpIDtrack);
      childIDs[0] = rowInPrimaryTrackTablePos;
      childIDs[1] = 0;
      outputParts(outputCollision.lastIndex(), v0.positivept(),
                  v0.positiveeta(), v0.positivephi(),
                  aod::femtouniverseparticle::ParticleType::kV0Child,
                  cutContainerV0.at(femtoUniverseV0Selection::V0ContainerPosition::kPosCuts),
                  cutContainerV0.at(femtoUniverseV0Selection::V0ContainerPosition::kPosPID),
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
                  cutContainerV0.at(femtoUniverseV0Selection::V0ContainerPosition::kNegCuts),
                  cutContainerV0.at(femtoUniverseV0Selection::V0ContainerPosition::kNegPID),
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
                  cutContainerV0.at(femtoUniverseV0Selection::V0ContainerPosition::kV0),
                  0,
                  v0.v0cosPA(),
                  indexChildID,
                  v0.mLambda(),
                  v0.mAntiLambda());
      if (ConfIsDebug) {
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
                      0, // cutContainerV0.at(femtoUniverseV0Selection::V0ContainerPosition::kPosCuts),
                      0, // cutContainerV0.at(femtoUniverseV0Selection::V0ContainerPosition::kPosPID),
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
                      0, // cutContainerV0.at(femtoUniverseV0Selection::V0ContainerPosition::kNegCuts),
                      0, // cutContainerV0.at(femtoUniverseV0Selection::V0ContainerPosition::kNegPID),
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
                      0, // cutContainerV0.at(femtoUniverseV0Selection::V0ContainerPosition::kNegCuts),
                      0, // cutContainerV0.at(femtoUniverseV0Selection::V0ContainerPosition::kNegPID),
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
                      0, // cutContainerV0.at(femtoUniverseV0Selection::V0ContainerPosition::kV0),
                      0,
                      0,
                      indexCascChildID,
                      ConfIsSelectCascOmega ? casc.mOmega() : casc.mXi(),
                      ConfIsSelectCascOmega ? casc.mOmega() : casc.mXi(),
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
      if (ConfIsDebug) {
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

    for (auto const& hfCand : hfCands) {

      if (!(hfCand.hfflag() & 1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        continue;
      }

      if (ConfD0Selection.ConfD0D0barCandMaxY >= 0. && std::abs(hfHelper.yD0(hfCand)) > ConfD0Selection.ConfD0D0barCandMaxY) {
        continue;
      }

      if (std::abs(hfCand.eta()) > ConfD0Selection.ConfD0D0barCandEtaCut) {
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
        isD0D0bar = true;
        daughFlag = 0;
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
                    -999, // cutContainerV0.at(femtoUniverseV0Selection::V0ContainerPosition::kPosCuts),
                    -999, // cutContainerV0.at(femtoUniverseV0Selection::V0ContainerPosition::kPosPID),
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
                    -999, // cutContainerV0.at(femtoUniverseV0Selection::V0ContainerPosition::kNegCuts),
                    -999, // cutContainerV0.at(femtoUniverseV0Selection::V0ContainerPosition::kNegPID),
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
                    -999, // cut, cutContainerType
                    -999, // PID, cutContainerType
                    -999,
                    indexChildID,
                    invMassD0,     // D0 mass (mLambda)
                    invMassD0bar); // D0bar mass (mAntiLambda)

        if (ConfIsDebug) {
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
    for (auto& [p1, p2] : combinations(soa::CombinationsFullIndexPolicy(tracks, tracks))) {
      if (!trackCuts.isSelectedMinimal(p1) || !trackCuts.isSelectedMinimal(p1)) {
        continue;
      }
      // implementing PID cuts for phi children
      if (ConfPhiSelection.ConfLooseTPCNSigma.value) {
        if (!(IsKaonNSigmaTPCLoose(p1.pt(), trackCuts.getNsigmaTPC(p1, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p1, o2::track::PID::Kaon)))) {
          continue;
        }
        if (!(IsKaonNSigmaTPCLoose(p2.pt(), trackCuts.getNsigmaTPC(p2, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p2, o2::track::PID::Kaon)))) {
          continue;
        }
      }
      if (ConfPhiSelection.ConfLooseTOFNSigma.value) {
        if (!(IsKaonNSigmaTOFLoose(p1.pt(), trackCuts.getNsigmaTPC(p1, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p1, o2::track::PID::Kaon)))) {
          continue;
        }
        if (!(IsKaonNSigmaTOFLoose(p2.pt(), trackCuts.getNsigmaTPC(p2, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p2, o2::track::PID::Kaon)))) {
          continue;
        }
      } else {
        if (!(IsKaonNSigma(p1.pt(), trackCuts.getNsigmaTPC(p1, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p1, o2::track::PID::Kaon)))) {
          continue;
        }
        if (!(IsKaonNSigma(p2.pt(), trackCuts.getNsigmaTPC(p2, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p2, o2::track::PID::Kaon)))) {
          continue;
        }
      }
      if (IsKaonRejected(p1.p(), trackCuts.getNsigmaTPC(p1, o2::track::PID::Proton), trackCuts.getNsigmaTOF(p1, o2::track::PID::Proton), trackCuts.getNsigmaTPC(p1, o2::track::PID::Pion), trackCuts.getNsigmaTOF(p1, o2::track::PID::Pion))) {
        continue;
      }
      if (IsKaonRejected(p2.p(), trackCuts.getNsigmaTPC(p2, o2::track::PID::Proton), trackCuts.getNsigmaTOF(p2, o2::track::PID::Proton), trackCuts.getNsigmaTPC(p2, o2::track::PID::Pion), trackCuts.getNsigmaTOF(p2, o2::track::PID::Pion))) {
        continue;
      }

      if ((!(p1.sign() == 1)) || (!(p2.sign() == -1))) {
        continue;
      }

      TLorentzVector part1Vec;
      TLorentzVector part2Vec;

      float mMassOne = TDatabasePDG::Instance()->GetParticle(321)->Mass();  // FIXME: Get from the PDG service of the common header
      float mMassTwo = TDatabasePDG::Instance()->GetParticle(-321)->Mass(); // FIXME: Get from the PDG service of the common header

      part1Vec.SetPtEtaPhiM(p1.pt(), p1.eta(), p1.phi(), mMassOne);
      part2Vec.SetPtEtaPhiM(p2.pt(), p2.eta(), p2.phi(), mMassTwo);

      TLorentzVector sumVec(part1Vec);
      sumVec += part2Vec;

      float phiEta = sumVec.Eta();
      if (TMath::Abs(phiEta) > 0.8) {
        continue;
      }

      float phiPt = sumVec.Pt();
      if ((phiPt < ConfPhiSelection.ConfPtLowLimitPhi.value) || (phiPt > ConfPhiSelection.ConfPtHighLimitPhi.value)) {
        continue;
      }

      float phiPhi = sumVec.Phi();
      if (sumVec.Phi() < 0) {
        phiPhi = sumVec.Phi() + 2 * o2::constants::math::PI;
      } else if (sumVec.Phi() >= 0) {
        phiPhi = sumVec.Phi();
      }
      float phiM = sumVec.M();

      if (((phiM < ConfPhiSelection.ConfInvMassLowLimitPhi.value) || (phiM > ConfPhiSelection.ConfInvMassUpLimitPhi.value))) {
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
                  -999, // cutContainerV0.at(femtoUniverseV0Selection::V0ContainerPosition::kPosCuts),
                  -999, // cutContainerV0.at(femtoUniverseV0Selection::V0ContainerPosition::kPosPID),
                  p1.dcaXY(),
                  childIDs,
                  0,
                  1); // sign, workaround for now
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
                  -999, // cutContainerV0.at(femtoUniverseV0Selection::V0ContainerPosition::kNegCuts),
                  -999, // cutContainerV0.at(femtoUniverseV0Selection::V0ContainerPosition::kNegPID),
                  p2.dcaXY(),
                  childIDs,
                  0,
                  -1); // sign, workaround for now
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
                  -999, // cutContainerV0.at(femtoUniverseV0Selection::V0ContainerPosition::kV0),
                  0,
                  phiM, // v0.v0cosPA(),
                  indexChildID,
                  phiM,  // phi.mLambda(), //for now it will have a mLambda getter, maybe we will change it in the future so it's more logical
                  -999); // v0.mAntiLambda()

      if (ConfIsDebug) {
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

    for (auto& particle : tracks) {
      /// if the most open selection criteria are not fulfilled there is no
      /// point looking further at the track

      if (particle.eta() < -ConfFilterCuts.ConfEtaFilterCut || particle.eta() > ConfFilterCuts.ConfEtaFilterCut)
        continue;
      if (particle.pt() < ConfFilterCuts.ConfPtLowFilterCut || particle.pt() > ConfFilterCuts.ConfPtHighFilterCut)
        continue;

      uint32_t pdgCode = static_cast<uint32_t>(particle.pdgCode());

      if (ConfMCTruthAnalysisWithPID) {
        bool pass = false;
        std::vector<int> tmpPDGCodes = ConfMCTruthPDGCodes; // necessary due to some features of the Configurable
        for (uint32_t pdg : tmpPDGCodes) {
          if (static_cast<int>(pdg) == static_cast<int>(pdgCode)) {
            if (pdgCode == 333) { // && (recoMcIds && recoMcIds->get().contains(particle.globalIndex()))) { // ATTENTION: all Phi mesons are NOT primary particles
              pass = true;
            } else {
              if (particle.isPhysicalPrimary() || (ConfActivateSecondaries && recoMcIds && recoMcIds->get().contains(particle.globalIndex())))
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
      // auto cutContainer = trackCuts.getCutContainer<aod::femtouniverseparticle::cutContainerType>(track);
      // instead of the bitmask, the PDG of the particle is stored as uint32_t

      // now the table is filled
      if constexpr (resolveDaughs) {
        tmpIDtrack.push_back(particle.globalIndex());
        continue;
      }
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
      if (ConfIsDebug) {
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
        if (ConfIsDebug) {
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
    fillCollisions<isMC>(col, tracks);
    fillTracks<isMC>(tracks);
    if (ConfIsActivateV0) {
      fillV0<isMC>(col, fullV0s, tracks);
    }
    if (ConfIsActivatePhi) {
      fillPhi<isMC>(col, tracks);
    }
    // if (ConfIsActivateCascade) {
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
  PROCESS_SWITCH(femtoUniverseProducerTask, processFullData, "Provide experimental data", false);

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
  PROCESS_SWITCH(femtoUniverseProducerTask, processTrackV0, "Provide experimental data for track v0", false);

  void processTrackCascadeData(aod::FemtoFullCollision const& col,
                               aod::BCsWithTimestamps const&,
                               soa::Filtered<aod::FemtoFullTracks> const& tracks,
                               o2::aod::CascDatas const& fullCascades)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    fillCollisions<false>(col, tracks);
    fillTracks<false>(tracks);
    fillCascade<false>(col, fullCascades, tracks);
  }
  PROCESS_SWITCH(femtoUniverseProducerTask, processTrackCascadeData, "Provide experimental data for track cascades", false);

  void processTrackV0CentRun3(aod::FemtoFullCollisionCentRun3 const& col,
                              aod::BCsWithTimestamps const&,
                              soa::Filtered<aod::FemtoFullTracks> const& tracks,
                              o2::aod::V0Datas const& fullV0s)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    fillCollisionsCentRun3<false>(col, tracks);
    fillTracks<false>(tracks);
    fillV0<false>(col, fullV0s, tracks);
  }
  PROCESS_SWITCH(femtoUniverseProducerTask, processTrackV0CentRun3, "Provide experimental data for track v0", false);

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
  PROCESS_SWITCH(femtoUniverseProducerTask, processFullMC, "Provide MC data (tracks, V0, Phi)", false);

  void processTrackMC(aod::FemtoFullCollisionMC const& col,
                      aod::BCsWithTimestamps const&,
                      soa::Join<aod::FemtoFullTracks, aod::McTrackLabels> const& tracks,
                      aod::McCollisions const&,
                      aod::McParticles const&)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    fillCollisions<true>(col, tracks);
    fillTracks<true>(tracks);
  }
  PROCESS_SWITCH(femtoUniverseProducerTask, processTrackMC, "Provide MC data for track analysis", false);

  void processTrackPhiMC(aod::FemtoFullCollisionMC const& col,
                         aod::BCsWithTimestamps const&,
                         soa::Join<aod::FemtoFullTracks, aod::McTrackLabels> const& tracks,
                         aod::McCollisions const&,
                         aod::McParticles const&)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    fillCollisions<true>(col, tracks);
    fillTracks<true>(tracks);
    fillPhi<true>(col, tracks);
  }
  PROCESS_SWITCH(femtoUniverseProducerTask, processTrackPhiMC, "Provide MC data for track Phi analysis", false);

  void processTrackD0MC(aod::FemtoFullCollisionMC const& col,
                        aod::BCsWithTimestamps const&,
                        soa::Join<aod::FemtoFullTracks, aod::McTrackLabels> const& tracks,
                        aod::McCollisions const&,
                        aod::McParticles const&)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    fillCollisions<true>(col, tracks);
    fillTracks<true>(tracks);
    // fillD0mesons<true>(col, tracks, candidates);
  }
  PROCESS_SWITCH(femtoUniverseProducerTask, processTrackD0MC, "Provide MC data for track D0 analysis", false);

  void processTrackData(aod::FemtoFullCollision const& col,
                        aod::BCsWithTimestamps const&,
                        aod::FemtoFullTracks const& tracks)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    fillCollisions<false>(col, tracks);
    fillTracks<false>(tracks);
  }
  PROCESS_SWITCH(femtoUniverseProducerTask, processTrackData,
                 "Provide experimental data for track track", true);

  // using FilteredFemtoFullTracks = soa::Filtered<FemtoFullTracks>;
  void processTrackPhiData(aod::FemtoFullCollision const& col,
                           aod::BCsWithTimestamps const&,
                           soa::Filtered<aod::FemtoFullTracks> const& tracks)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    fillCollisions<false>(col, tracks);
    fillTracks<false>(tracks);
    fillPhi<false>(col, tracks);
  }
  PROCESS_SWITCH(femtoUniverseProducerTask, processTrackPhiData,
                 "Provide experimental data for track phi", false);

  void processTrackD0mesonData(aod::FemtoFullCollision const& col,
                               aod::BCsWithTimestamps const&,
                               soa::Filtered<aod::FemtoFullTracks> const& tracks,
                               soa::Join<aod::HfCand2Prong, aod::HfSelD0> const& candidates)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    fillCollisions<false>(col, tracks);
    fillTracks<false>(tracks);
    fillD0mesons<false>(col, tracks, candidates);
  }
  PROCESS_SWITCH(femtoUniverseProducerTask, processTrackD0mesonData,
                 "Provide experimental data for track D0 meson", false);

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
  PROCESS_SWITCH(femtoUniverseProducerTask, processTrackMCTruth, "Provide MC data for MC truth track analysis", false);

  void processTrackMCGen(aod::McCollision const& col,
                         aod::McParticles const& mcParticles)
  {
    outputCollision(col.posZ(), 0, 0, 2, 0);
    fillParticles(mcParticles);
  }
  PROCESS_SWITCH(femtoUniverseProducerTask, processTrackMCGen, "Provide MC Generated for model comparisons", false);

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
    for (auto& col : collisions) {
      auto groupedTracks = tracks.sliceBy(perCollisionTracks, col.globalIndex());
      auto groupedV0s = fullV0s.sliceBy(perCollisionV0s, col.globalIndex());
      getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
      fillCollisionsAndTracksAndV0AndPhi<true>(col, groupedTracks, groupedV0s);
      for (auto& track : groupedTracks) {
        if (trackCuts.isSelectedMinimal(track))
          recoMcIds.insert(track.mcParticleId());
      }
    }

    // truth
    for (auto& mccol : mccols) {
      auto groupedMCParticles = mcParticles.sliceBy(perMCCollision, mccol.globalIndex());
      auto groupedCollisions = collisions.sliceBy(recoCollsPerMCColl, mccol.globalIndex());
      fillMCTruthCollisions(groupedCollisions, groupedMCParticles);                           // fills the reco collisions for mc collision
      fillParticles<decltype(groupedMCParticles), true, true>(groupedMCParticles, recoMcIds); // fills mc particles
    }
  }
  PROCESS_SWITCH(femtoUniverseProducerTask, processTruthAndFullMC, "Provide both MC truth and reco for tracks and V0s", false);

  void processFullMCCent(aod::FemtoFullCollisionCentRun3 const& col,
                         aod::BCsWithTimestamps const&,
                         soa::Join<aod::FemtoFullTracks, aod::McTrackLabels> const& tracks,
                         aod::McCollisions const&,
                         aod::McParticles const&,
                         soa::Join<o2::aod::V0Datas, aod::McV0Labels> const& fullV0s)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    fillCollisionsCentRun3<true>(col, tracks);
    fillTracks<true>(tracks);
    fillV0<true>(col, fullV0s, tracks);
  }
  PROCESS_SWITCH(femtoUniverseProducerTask, processFullMCCent, "Provide MC data with centrality bins", false);

  void processTrackCentRun2Data(aod::FemtoFullCollisionCentRun2 const& col,
                                aod::BCsWithTimestamps const&,
                                aod::FemtoFullTracks const& tracks)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    fillCollisionsCentRun2<false>(col, tracks);
    fillTracks<false>(tracks);
  }
  PROCESS_SWITCH(femtoUniverseProducerTask, processTrackCentRun2Data, "Provide experimental data for Run 2 with centrality for track track", false);

  void processTrackCentRun3Data(aod::FemtoFullCollisionCentRun3 const& col,
                                aod::BCsWithTimestamps const&,
                                soa::Filtered<aod::FemtoFullTracks> const& tracks)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    fillCollisionsCentRun3<false>(col, tracks);
    fillTracks<false>(tracks);
  }
  PROCESS_SWITCH(femtoUniverseProducerTask, processTrackCentRun3Data, "Provide experimental data for Run 3 with centrality for track track", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<femtoUniverseProducerTask>(cfgc)};
  return workflow;
}
