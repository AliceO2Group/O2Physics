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
/// \author Zuzanna Chochulska, WUT Warsaw, zuzanna.chochulska.stud@pw.edu.pl
/// \author Malgorzata Janik, WUT Warsaw, majanik@cern.ch

#include <CCDB/BasicCCDBManager.h>
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseCollisionSelection.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseTrackSelection.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseV0Selection.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniversePhiSelection.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUtils.h"
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

using namespace o2;
using namespace o2::analysis::femtoUniverse;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{

using FemtoFullCollision =
  soa::Join<aod::Collisions, aod::EvSels, aod::Mults>::iterator;
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

  Configurable<bool> ConfIsDebug{"ConfIsDebug", true, "Enable Debug tables"};
  // Choose if filtering or skimming version is run
  Configurable<bool> ConfIsTrigger{"ConfIsTrigger", false, "Store all collisions"};
  // Choose if running on converted data or Run3  / Pilot
  Configurable<bool> ConfIsRun3{"ConfIsRun3", false, "Running on Run3 or pilot"};
  Configurable<bool> ConfIsMC{"ConfIsMC", false, "Running on MC; implemented only for Run3"};

  Configurable<bool> ConfIsForceGRP{"ConfIsForceGRP", false, "Set true if the magnetic field configuration is not available in the usual CCDB directory (e.g. for Run 2 converted data or unanchorad Monte Carlo)"};

  /// Event cuts
  FemtoUniverseCollisionSelection colCuts;
  Configurable<bool> ConfEvtUseTPCmult{"ConfEvtUseTPCmult", false, "Use multiplicity based on the number of tracks with TPC information"};
  Configurable<float> ConfEvtZvtx{"ConfEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};
  Configurable<bool> ConfEvtTriggerCheck{"ConfEvtTriggerCheck", true, "Evt sel: check for trigger"};
  Configurable<int> ConfEvtTriggerSel{"ConfEvtTriggerSel", kINT7, "Evt sel: trigger"};
  Configurable<bool> ConfEvtOfflineCheck{"ConfEvtOfflineCheck", false, "Evt sel: check for offline selection"};
  Configurable<bool> ConfIsActivateV0{"ConfIsActivateV0", true, "Activate filling of V0 into femtouniverse tables"};
  Configurable<bool> ConfIsActivatePhi{"ConfIsActivatePhi", true, "Activate filling of Phi into femtouniverse tables"};
  Configurable<bool> ConfMCTruthAnalysisWithPID{"ConfMCTruthAnalysisWithPID", true, "1: take only particles with specified PDG, 0: all particles (for MC Truth)"};
  Configurable<std::vector<int>> ConfMCTruthPDGCodes{"ConfMCTruthPDGCodes", std::vector<int>{211, -211, 2212, -2212, 333}, "PDG of particles to be stored"};

  // just sanity check to make sure in case there are problems in conversion or
  // MC production it does not affect results
  Configurable<bool> ConfTrkRejectNotPropagated{"ConfTrkRejectNotPropagated", false, "True: reject not propagated tracks"};
  // Configurable<bool> ConfRejectITSHitandTOFMissing{
  //     "ConfRejectITSHitandTOFMissing", false,
  //     "True: reject if neither ITS hit nor TOF timing satisfied"};

  Configurable<int> ConfTrkPDGCode{"ConfTrkPDGCode", 2212, "PDG code of the selected track for Monte Carlo truth"}; // only for checking the particle origin (particles for which PDG does not match are marked as "Fake")
  FemtoUniverseTrackSelection trackCuts;

  struct : o2::framework::ConfigurableGroup {
    Configurable<float> ConfPtLowFilterCut{"ConfPtLowFilterCut", 0.14, "Lower limit for Pt for the global track"};   // pT low
    Configurable<float> ConfPtHighFilterCut{"ConfPtHighFilterCut", 5.0, "Higher limit for Pt for the global track"}; // pT high
    Configurable<float> ConfEtaFilterCut{"ConfEtaFilterCut", 0.8, "Eta cut for the global track"};                   // eta
    Configurable<float> ConfDcaXYFilterCut{"ConfDcaXYFilterCut", 2.4, "Value for DCA_XY for the global track"};      // max dca to vertex XY
    Configurable<float> ConfDcaZFilterCut{"ConfDcaZFilterCut", 3.2, "Value for DCA_Z for the global track"};         // max dca to vertex Z
  } ConfFilterCuts;

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
  Configurable<std::vector<int>> ConfTrkPIDspecies{"ConfTrkPIDspecies", std::vector<int>{o2::track::PID::Pion, o2::track::PID::Kaon, o2::track::PID::Proton, o2::track::PID::Deuteron}, "Trk sel: Particles species for PID"};
  // Numbers from ~/alice/O2/DataFormats/Reconstruction/include/ReconstructionDataFormats/PID.h //static constexpr ID Pion = 2; static constexpr ID Kaon = 3; static constexpr ID Proton = 4; static constexpr ID Deuteron = 5;

  // TrackSelection *o2PhysicsTrackSelection;
  /// \todo Labeled array (see Track-Track task)

  // V0
  FemtoUniverseV0Selection v0Cuts;
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

  Filter GlobalCutFilter = requireGlobalTrackInFilter();

  Filter CustomTrackFilter = (aod::track::pt > ConfFilterCuts.ConfPtLowFilterCut) &&
                             (aod::track::pt < ConfFilterCuts.ConfPtHighFilterCut) &&
                             (nabs(aod::track::eta) < ConfFilterCuts.ConfEtaFilterCut) &&
                             (aod::track::dcaXY < ConfFilterCuts.ConfDcaXYFilterCut) &&
                             (aod::track::dcaZ < ConfFilterCuts.ConfDcaZFilterCut);

  // PHI
  FemtoUniversePhiSelection phiCuts;
  struct : o2::framework::ConfigurableGroup {
    Configurable<std::vector<float>> ConfPhiSign{FemtoUniversePhiSelection::getSelectionName(femtoUniversePhiSelection::kPhiSign, "ConfPhi"), std::vector<float>{-1, 1}, FemtoUniversePhiSelection::getSelectionHelper(femtoUniversePhiSelection::kPhiSign, "Phi selection: ")};
    Configurable<std::vector<float>> ConfPhiPtMin{FemtoUniversePhiSelection::getSelectionName(femtoUniversePhiSelection::kPhipTMin, "ConfPhi"), std::vector<float>{0.3f, 0.4f, 0.5f}, FemtoUniversePhiSelection::getSelectionHelper(femtoUniversePhiSelection::kPhipTMin, "Phi selection: ")};
    Configurable<std::vector<float>> ConfPhiPtMax{FemtoUniversePhiSelection::getSelectionName(femtoUniversePhiSelection::kPhipTMax, "ConfPhi"), std::vector<float>{3.3f, 3.4f, 3.5f}, FemtoUniversePhiSelection::getSelectionHelper(femtoUniversePhiSelection::kPhipTMax, "Phi selection: ")};
    Configurable<std::vector<float>> ConfPhiEtaMax{FemtoUniversePhiSelection::getSelectionName(femtoUniversePhiSelection::kPhietaMax, "ConfPhi"), std::vector<float>{0.8f, 0.7f, 0.9f}, FemtoUniversePhiSelection::getSelectionHelper(femtoUniversePhiSelection::kPhietaMax, "Phi selection: ")};
    Configurable<std::vector<float>> ConfPhiDCADaughMax{FemtoUniversePhiSelection::getSelectionName(femtoUniversePhiSelection::kPhiDCADaughMax, "ConfPhi"), std::vector<float>{1.2f, 1.5f}, FemtoUniversePhiSelection::getSelectionHelper(femtoUniversePhiSelection::kPhiDCADaughMax, "Phi selection: ")};
    Configurable<std::vector<float>> ConfPhiCPAMin{FemtoUniversePhiSelection::getSelectionName(femtoUniversePhiSelection::kPhiCPAMin, "ConfPhi"), std::vector<float>{0.99f, 0.995f}, FemtoUniversePhiSelection::getSelectionHelper(femtoUniversePhiSelection::kPhiCPAMin, "Phi selection: ")};
    Configurable<std::vector<float>> ConfPhiTranRadMin{FemtoUniversePhiSelection::getSelectionName(femtoUniversePhiSelection::kPhiTranRadMin, "ConfPhi"), std::vector<float>{0.2f}, FemtoUniversePhiSelection::getSelectionHelper(femtoUniversePhiSelection::kPhiTranRadMin, "Phi selection: ")};
    Configurable<std::vector<float>> ConfPhiTranRadMax{FemtoUniversePhiSelection::getSelectionName(femtoUniversePhiSelection::kPhiTranRadMax, "ConfPhi"), std::vector<float>{100.f}, FemtoUniversePhiSelection::getSelectionHelper(femtoUniversePhiSelection::kPhiTranRadMax, "Phi selection: ")};
    Configurable<std::vector<float>> ConfPhiDecVtxMax{FemtoUniversePhiSelection::getSelectionName(femtoUniversePhiSelection::kPhiDecVtxMax, "ConfPhi"), std::vector<float>{100.f}, FemtoUniversePhiSelection::getSelectionHelper(femtoUniversePhiSelection::kPhiDecVtxMax, "Phi selection: ")};
  } ConfPhiSelection;

  struct : o2::framework::ConfigurableGroup {
    Configurable<std::vector<float>> ConfPhiChildCharge{"ConfPhiChildSign", std::vector<float>{-1, 1}, "Phi Child sel: Charge"};
    Configurable<std::vector<float>> ConfPhiChildEtaMax{"ConfPhiChildEtaMax", std::vector<float>{0.8f}, "Phi Child sel: max eta"};
    Configurable<std::vector<float>> ConfPhiChildTPCnClsMin{"ConfPhiChildTPCnClsMin", std::vector<float>{80.f, 70.f, 60.f}, "Phi Child sel: Min. nCls TPC"};
    Configurable<std::vector<float>> ConfPhiChildDCAMin{"ConfPhiChildDCAMin", std::vector<float>{0.05f, 0.06f}, "Phi Child sel:  Max. DCA Daugh to PV (cm)"};
    Configurable<std::vector<float>> ConfPhiChildPIDnSigmaMax{"ConfPhiChildPIDnSigmaMax", std::vector<float>{5.f, 4.f}, "Phi Child sel: Max. PID nSigma TPC"};
    Configurable<std::vector<int>> ConfPhiChildPIDspecies{"ConfPhiChildPIDspecies", std::vector<int>{o2::track::PID::Kaon, o2::track::PID::Kaon}, "Phi Child sel: Particles species for PID"};
  } ConfPhiChildSelection;

  struct : o2::framework::ConfigurableGroup {
    Configurable<float> ConfNsigmaCombinedKaon{"ConfNsigmaCombinedKaon", 3.0, "TPC and TOF Kaon Sigma (combined) for momentum > 0.4"};
    Configurable<float> ConfNsigmaTPCKaon{"ConfNsigmaTPCKaon", 3.0, "TPC Kaon Sigma for momentum < 0.4"};
    Configurable<bool> ConfNsigmaTPCTOFKaon{"ConfNsigmaTPCTOFKaon", true, "Use TPC and TOF for PID of Kaons"};
    Configurable<float> ConfInvMassLowLimitPhi{"ConfInvMassLowLimitPhi", 1.011, "Lower limit of the Phi invariant mass"}; // change that to do invariant mass cut
    Configurable<float> ConfInvMassUpLimitPhi{"ConfInvMassUpLimitPhi", 1.027, "Upper limit of the Phi invariant mass"};
  } ConfPhiCommon;
  // PHI child one
  struct : o2::framework::ConfigurableGroup {
    Configurable<int> ConfPDGCodePartOne{"ConfPDGCodePartOne", 321, "Particle 1 - PDG code"};
  } ConfPhiChildOne;
  // PHI child two
  struct : o2::framework::ConfigurableGroup {
    Configurable<int> ConfPDGCodePartTwo{"ConfPDGCodePartTwo", 321, "Particle 2 - PDG code"};
  } ConfPhiChildTwo;

  bool IsKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK)
  {
    //|nsigma_TPC| < 5 for p < 0.4 GeV/c
    //|nsigma_combined| < 5 for p > 0.4

    // using configurables:
    // ConfNsigmaTPCTOFKaon -> are we doing TPC TOF PID for Kaons? (boolean)
    // ConfNsigmaTPCKaon -> TPC Kaon Sigma for momentum < 0.4
    // ConfNsigmaCombinedKaon -> TPC and TOF Kaon Sigma (combined) for momentum > 0.4

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

  /// \todo should we add filter on min value pT/eta of V0 and daughters?
  /*Filter v0Filter = (nabs(aod::v0data::x) < V0DecVtxMax.value) &&
                    (nabs(aod::v0data::y) < V0DecVtxMax.value) &&
                    (nabs(aod::v0data::z) < V0DecVtxMax.value);*/
  // (aod::v0data::v0radius > V0TranRadV0Min.value); to be added, not working
  // for now do not know why

  HistogramRegistry qaRegistry{"QAHistos", {}, OutputObjHandlingPolicy::QAObject};

  int mRunNumber;
  float mMagField;
  Service<o2::ccdb::BasicCCDBManager> ccdb; /// Accessing the CCDB

  void init(InitContext&)
  {
    if ((doprocessFullData || doprocessTrackPhiData || doprocessTrackData) == false && (doprocessFullMC || doprocessTrackMC) == false) {
      LOGF(fatal, "Neither processFullData nor processFullMC enabled. Please choose one.");
    }
    if ((doprocessFullData || doprocessTrackPhiData || doprocessTrackData) == true && (doprocessFullMC || doprocessTrackMC) == true) {
      LOGF(fatal,
           "Cannot enable process Data and process MC at the same time. "
           "Please choose one.");
    }

    colCuts.setCuts(ConfEvtZvtx, ConfEvtTriggerCheck, ConfEvtTriggerSel, ConfEvtOfflineCheck, ConfIsRun3);
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
      v0Cuts.setSelection(ConfV0Sign, femtoUniverseV0Selection::kV0Sign, femtoUniverseSelection::kEqual);
      v0Cuts.setSelection(ConfV0PtMin, femtoUniverseV0Selection::kV0pTMin, femtoUniverseSelection::kLowerLimit);
      v0Cuts.setSelection(ConfV0PtMax, femtoUniverseV0Selection::kV0pTMax, femtoUniverseSelection::kUpperLimit);
      v0Cuts.setSelection(ConfV0EtaMax, femtoUniverseV0Selection::kV0etaMax, femtoUniverseSelection::kAbsUpperLimit);
      v0Cuts.setSelection(ConfV0DCADaughMax, femtoUniverseV0Selection::kV0DCADaughMax, femtoUniverseSelection::kUpperLimit);
      v0Cuts.setSelection(ConfV0CPAMin, femtoUniverseV0Selection::kV0CPAMin, femtoUniverseSelection::kLowerLimit);
      v0Cuts.setSelection(ConfV0TranRadMin, femtoUniverseV0Selection::kV0TranRadMin, femtoUniverseSelection::kLowerLimit);
      v0Cuts.setSelection(ConfV0TranRadMax, femtoUniverseV0Selection::kV0TranRadMax, femtoUniverseSelection::kUpperLimit);
      v0Cuts.setSelection(ConfV0DecVtxMax, femtoUniverseV0Selection::kV0DecVtxMax, femtoUniverseSelection::kUpperLimit);
      v0Cuts.setChildCuts(femtoUniverseV0Selection::kPosTrack, ConfChildCharge, femtoUniverseTrackSelection::kSign, femtoUniverseSelection::kEqual);
      v0Cuts.setChildCuts(femtoUniverseV0Selection::kPosTrack, ConfChildEtaMax, femtoUniverseTrackSelection::kEtaMax, femtoUniverseSelection::kAbsUpperLimit);
      v0Cuts.setChildCuts(femtoUniverseV0Selection::kPosTrack, ConfChildTPCnClsMin, femtoUniverseTrackSelection::kTPCnClsMin, femtoUniverseSelection::kLowerLimit);
      v0Cuts.setChildCuts(femtoUniverseV0Selection::kPosTrack, ConfChildDCAMin, femtoUniverseTrackSelection::kDCAMin, femtoUniverseSelection::kAbsLowerLimit);
      v0Cuts.setChildCuts(femtoUniverseV0Selection::kPosTrack, ConfChildPIDnSigmaMax, femtoUniverseTrackSelection::kPIDnSigmaMax, femtoUniverseSelection::kAbsUpperLimit);
      v0Cuts.setChildCuts(femtoUniverseV0Selection::kNegTrack, ConfChildCharge, femtoUniverseTrackSelection::kSign, femtoUniverseSelection::kEqual);
      v0Cuts.setChildCuts(femtoUniverseV0Selection::kNegTrack, ConfChildEtaMax, femtoUniverseTrackSelection::kEtaMax, femtoUniverseSelection::kAbsUpperLimit);
      v0Cuts.setChildCuts(femtoUniverseV0Selection::kNegTrack, ConfChildTPCnClsMin, femtoUniverseTrackSelection::kTPCnClsMin, femtoUniverseSelection::kLowerLimit);
      v0Cuts.setChildCuts(femtoUniverseV0Selection::kNegTrack, ConfChildDCAMin, femtoUniverseTrackSelection::kDCAMin, femtoUniverseSelection::kAbsLowerLimit);
      v0Cuts.setChildCuts(femtoUniverseV0Selection::kNegTrack, ConfChildPIDnSigmaMax, femtoUniverseTrackSelection::kPIDnSigmaMax, femtoUniverseSelection::kAbsUpperLimit);
      v0Cuts.setChildPIDSpecies(femtoUniverseV0Selection::kPosTrack, ConfChildPIDspecies);
      v0Cuts.setChildPIDSpecies(femtoUniverseV0Selection::kNegTrack, ConfChildPIDspecies);
      v0Cuts.init<aod::femtouniverseparticle::ParticleType::kV0, aod::femtouniverseparticle::ParticleType::kV0Child, aod::femtouniverseparticle::cutContainerType>(&qaRegistry);
      v0Cuts.setInvMassLimits(ConfV0InvMassLowLimit, ConfV0InvMassUpLimit);

      v0Cuts.setChildRejectNotPropagatedTracks(femtoUniverseV0Selection::kPosTrack, ConfTrkRejectNotPropagated);
      v0Cuts.setChildRejectNotPropagatedTracks(femtoUniverseV0Selection::kNegTrack, ConfTrkRejectNotPropagated);

      v0Cuts.setnSigmaPIDOffsetTPC(ConfTrkPIDnSigmaOffsetTPC);
      v0Cuts.setChildnSigmaPIDOffset(femtoUniverseV0Selection::kPosTrack, ConfTrkPIDnSigmaOffsetTPC, ConfTrkPIDnSigmaOffsetTOF);
      v0Cuts.setChildnSigmaPIDOffset(femtoUniverseV0Selection::kNegTrack, ConfTrkPIDnSigmaOffsetTPC, ConfTrkPIDnSigmaOffsetTOF);

      if (ConfV0RejectKaons) {
        v0Cuts.setKaonInvMassLimits(ConfV0InvKaonMassLowLimit, ConfV0InvKaonMassUpLimit);
      }
      // if (ConfRejectITSHitandTOFMissing) {
      //   o2PhysicsTrackSelection = new
      //   TrackSelection(getGlobalTrackSelection());
      //   o2PhysicsTrackSelection->SetRequireHitsInITSLayers(1, {0, 1, 2, 3});
      // }
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

  template <bool isTrackOrV0, bool isPhi, typename ParticleType>
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
    } else if constexpr (isPhi) {
      outputDebugParts(-999., -999., -999., -999., -999., -999., -999., -999.,
                       -999., -999., -999., -999., -999., -999., -999., -999.,
                       -999., -999., -999., -999., -999.,
                       -999., -999.,
                       -999., -999., -999.,
                       -999.); // QA for phi
    } else {
      LOGF(info, "isTrack0orV0: %d, isPhi: %d", isTrackOrV0, isPhi);
      outputDebugParts(-999., -999., -999., -999., -999., -999., -999., -999.,
                       -999., -999., -999., -999., -999., -999., -999., -999.,
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
      auto motherparticleMC = particleMC.template mothers_as<aod::McParticles>().front();

      if (abs(pdgCode) == abs(ConfTrkPDGCode.value)) {
        if (particleMC.isPhysicalPrimary()) {
          particleOrigin = aod::femtouniverseMCparticle::ParticleOriginMCTruth::kPrimary;
        } else if (motherparticleMC.producedByGenerator()) {
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
  template <bool isMC, typename CollisionType, typename TrackType>
  void fillCollisions(CollisionType const& col, TrackType const& tracks)
  {
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
    if (ConfEvtUseTPCmult) {
      multNtr = col.multTPC();
    }

    // check whether the basic event selection criteria are fulfilled
    // if the basic selection is NOT fulfilled:
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
    outputCollision(vtxZ, mult, multNtr, spher, mMagField);
  }

  template <typename CollisionType, typename TrackType>
  void fillMCTruthCollisions(CollisionType const& col, TrackType const& tracks)
  {
    const auto vtxZ = col.posZ();
    const auto spher = 0; // colCuts.computeSphericity(col, tracks);
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

    // CHECK WHAT CUTS SHOULD BE USED FOR MC TRUTH
    //  if (!colCuts.isSelected(col)) {
    //    if (ConfIsTrigger) {
    //      outputCollision(vtxZ, mult, multNtr, spher, mMagField);
    //    }
    //    return;
    //  }

    colCuts.fillQA(col);
    outputCollision(vtxZ, mult, multNtr, spher, mMagField);
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

      trackCuts.fillQA<aod::femtouniverseparticle::ParticleType::kTrack,
                       aod::femtouniverseparticle::TrackType::kNoChild>(track);
      // the bit-wise container of the systematic variations is obtained
      auto cutContainer = trackCuts.getCutContainer<aod::femtouniverseparticle::cutContainerType>(track);

      // now the table is filled
      outputParts(outputCollision.lastIndex(), track.pt(), track.eta(),
                  track.phi(), aod::femtouniverseparticle::ParticleType::kTrack,
                  cutContainer.at(
                    femtoUniverseTrackSelection::TrackContainerPosition::kCuts),
                  cutContainer.at(
                    femtoUniverseTrackSelection::TrackContainerPosition::kPID),
                  track.dcaXY(), childIDs, 0, 0);
      tmpIDtrack.push_back(track.globalIndex());
      if (ConfIsDebug) {
        fillDebugParticle<true, false>(track);
      }

      if constexpr (isMC) {
        fillMCParticle(track, o2::aod::femtouniverseparticle::ParticleType::kTrack);
      }
    }
  }

  template <bool isMC, typename CollisionType, typename V0Type, typename TrackType>
  void fillV0(CollisionType const& col, V0Type const& fullV0s, TrackType const& tracks)
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
                  v0.v0cosPA(col.posX(), col.posY(), col.posZ()),
                  indexChildID,
                  v0.mLambda(),
                  v0.mAntiLambda());
      if (ConfIsDebug) {
        fillDebugParticle<true, false>(postrack); // QA for positive daughter
        fillDebugParticle<true, false>(negtrack); // QA for negative daughter
        fillDebugParticle<false, false>(v0);      // QA for v0
      }
      if constexpr (isMC) {
        fillMCParticle(v0, o2::aod::femtouniverseparticle::ParticleType::kV0);
      }
    }
  }

  template <bool isMC, typename TrackType, typename CollisionType>
  void fillPhi(CollisionType const& col, TrackType const& tracks)
  {
    std::vector<int> childIDs = {0, 0}; // these IDs are necessary to keep track of the children
    std::vector<int> tmpIDtrack;        // this vector keeps track of the matching of the primary track table row <-> aod::track table global index
    // lorentz vectors and filling the tables
    for (auto& [p1, p2] : combinations(soa::CombinationsStrictlyUpperIndexPolicy(tracks, tracks))) {
      // implementing PID cuts for phi children
      if (!(IsKaonNSigma(p1.pt(), trackCuts.getNsigmaTPC(p1, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p1, o2::track::PID::Kaon)))) {
        continue;
      }
      if (!(IsKaonNSigma(p2.pt(), trackCuts.getNsigmaTPC(p2, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p2, o2::track::PID::Kaon)))) {
        continue;
      } else if ((!(p1.sign() == 1)) || (!(p2.sign() == -1))) {
        continue;
      }

      TLorentzVector part1Vec;
      TLorentzVector part2Vec;

      float mMassOne = TDatabasePDG::Instance()->GetParticle(ConfPhiChildOne.ConfPDGCodePartOne)->Mass();
      float mMassTwo = TDatabasePDG::Instance()->GetParticle(ConfPhiChildTwo.ConfPDGCodePartTwo)->Mass();

      part1Vec.SetPtEtaPhiM(p1.pt(), p1.eta(), p1.phi(), mMassOne);
      part2Vec.SetPtEtaPhiM(p2.pt(), p2.eta(), p2.phi(), mMassTwo);

      TLorentzVector sumVec(part1Vec);
      sumVec += part2Vec;

      float phiEta = sumVec.Eta();
      if (TMath::Abs(phiEta) > 0.8) {
        continue;
      }

      float phiPt = sumVec.Pt();
      if ((phiPt < 0.14) || (phiPt > 10.0)) {
        continue;
      }

      float phiPhi = sumVec.Phi();
      if (sumVec.Phi() < 0) {
        phiPhi = sumVec.Phi() + 2 * o2::constants::math::PI;
      } else if (sumVec.Phi() >= 0) {
        phiPhi = sumVec.Phi();
      }
      float phiM = sumVec.M();

      if (((phiM < ConfPhiCommon.ConfInvMassLowLimitPhi) || (phiM > ConfPhiCommon.ConfInvMassUpLimitPhi))) {
        continue;
      }

      phiCuts.fillQA<aod::femtouniverseparticle::ParticleType::kPhi, aod::femtouniverseparticle::ParticleType::kPhiChild>(col, p1, p1, p2, ConfPhiChildOne.ConfPDGCodePartOne, ConfPhiChildTwo.ConfPDGCodePartTwo); ///\todo fill QA also for daughters

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
                  0);
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
                  0);
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
                  phiM, // v0.v0cosPA(col.posX(), col.posY(), col.posZ()),
                  indexChildID,
                  phiM,  // phi.mLambda(), //for now it will have a mLambda getter, maybe we will change it in the future so it's more logical
                  -999); // v0.mAntiLambda()

      if (ConfIsDebug) {
        fillDebugParticle<true, false>(p1); // QA for positive daughter
        fillDebugParticle<true, false>(p2); // QA for negative daughter
        fillDebugParticle<false, true>(p1); // QA for phi
      }
      // if constexpr (isMC) {
      //   fillMCParticle(v0, o2::aod::femtouniverseparticle::ParticleType::kV0);
      // }
    }
  }

  template <typename TrackType>
  void fillParticles(TrackType const& tracks)
  {
    std::vector<int> childIDs = {0, 0}; // these IDs are necessary to keep track of the children

    for (auto& particle : tracks) {
      /// if the most open selection criteria are not fulfilled there is no
      /// point looking further at the track
      if (!particle.isPhysicalPrimary())
        continue;

      uint32_t pdgCode = (uint32_t)particle.pdgCode();

      if (ConfMCTruthAnalysisWithPID) {
        bool pass = false;
        std::vector<int> tmpPDGCodes = ConfMCTruthPDGCodes; // necessary due to some features of the Configurable
        for (uint32_t pdg : tmpPDGCodes) {
          if (static_cast<int>(pdg) == static_cast<int>(pdgCode)) {
            pass = true;
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
  }

  void
    processFullData(aod::FemtoFullCollision const& col,
                    aod::BCsWithTimestamps const&,
                    aod::FemtoFullTracks const& tracks,
                    o2::aod::V0Datas const& fullV0s) /// \todo with FilteredFullV0s
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    fillCollisionsAndTracksAndV0AndPhi<false>(col, tracks, fullV0s);
  }
  PROCESS_SWITCH(femtoUniverseProducerTask, processFullData,
                 "Provide experimental data", false);

  void
    processFullMC(aod::FemtoFullCollisionMC const& col,
                  aod::BCsWithTimestamps const&,
                  soa::Join<aod::FemtoFullTracks, aod::McTrackLabels> const& tracks,
                  aod::McCollisions const& mcCollisions,
                  aod::McParticles const& mcParticles,
                  soa::Join<o2::aod::V0Datas, aod::McV0Labels> const& fullV0s) /// \todo with FilteredFullV0s
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    fillCollisionsAndTracksAndV0AndPhi<true>(col, tracks, fullV0s);
  }
  PROCESS_SWITCH(femtoUniverseProducerTask, processFullMC, "Provide MC data (tracks, V0, Phi)", false);

  void
    processTrackMC(aod::FemtoFullCollisionMC const& col,
                   aod::BCsWithTimestamps const&,
                   soa::Join<aod::FemtoFullTracks, aod::McTrackLabels> const& tracks,
                   aod::McCollisions const& mcCollisions,
                   aod::McParticles const&)
  {
    // get magnetic field for run
    getMagneticFieldTesla(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    fillCollisions<true>(col, tracks);
    fillTracks<true>(tracks);
  }
  PROCESS_SWITCH(femtoUniverseProducerTask, processTrackMC, "Provide MC data for track analysis", false);

  void
    processTrackData(aod::FemtoFullCollision const& col,
                     aod::BCsWithTimestamps const&,
                     aod::FemtoFullTracks const& tracks) /// \todo with FilteredFullV0s
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

  void
    processTrackMCTruth(aod::FemtoFullCollisionMC const& col,
                        aod::BCsWithTimestamps const&,
                        aod::McCollisions const& mcCollisions,
                        aod::McParticles const& mcParticles)
  {
    // magnetic field for run not needed for mc truth
    // fill the tables
    fillMCTruthCollisions(col, mcParticles);
    fillParticles(mcParticles);
  }
  PROCESS_SWITCH(femtoUniverseProducerTask, processTrackMCTruth, "Provide MC data for MC truth track analysis", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<femtoUniverseProducerTask>(cfgc)};
  return workflow;
}
