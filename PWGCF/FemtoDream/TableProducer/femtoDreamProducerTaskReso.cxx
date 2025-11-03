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

/// \file femtoDreamProducerTaskReso.cxx
/// \brief Tasks that produces the track tables used for the pairing
/// \author Laura Serksnyte, TU München, laura.serksnyte@tum.de

#include "PWGCF/DataModel/FemtoDerived.h"
#include "PWGCF/FemtoDream/Core/femtoDreamCascadeSelection.h"
#include "PWGCF/FemtoDream/Core/femtoDreamCollisionSelection.h"
#include "PWGCF/FemtoDream/Core/femtoDreamResoSelection.h"
#include "PWGCF/FemtoDream/Core/femtoDreamTrackSelection.h"
#include "PWGCF/FemtoDream/Core/femtoDreamUtils.h"
#include "PWGCF/FemtoDream/Core/femtoDreamV0SelectionK0Short.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "EventFiltering/Zorro.h"

#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Expressions.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"
#include <CCDB/BasicCCDBManager.h>

#include "Math/Vector4D.h"
#include "TMath.h"

#include <fairlogger/Logger.h>

#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::rctsel;
using namespace o2::analysis::femtoDream;

namespace o2::aod
{
using FemtoFullCollision = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms>::iterator;
using FemtoFullCollisionNoCent = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>::iterator;
using FemtoFullCollisionCentPbPb = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs>::iterator;
using FemtoFullCollisionMC = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms, aod::McCollisionLabels>::iterator;
using FemtoFullCollisionNoCentMC = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::McCollisionLabels>::iterator;
using FemtoFullCollisionMCCentPbPb = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs, aod::McCollisionLabels>::iterator;
using FemtoFullMCgenCollisions = soa::Join<aod::McCollisions, MultsExtraMC>;
using FemtoFullMCgenCollision = FemtoFullMCgenCollisions::iterator;

using FemtoFullTracks =
  soa::Join<aod::FullTracks, aod::TracksDCA,
            aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa,
            aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTPCFullTr, aod::pidTPCFullHe,
            aod::pidTOFFullEl, aod::pidTOFFullPi, aod::pidTOFFullKa,
            aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFFullTr, aod::pidTOFFullHe>;

} // namespace o2::aod

namespace software_triggers
{
static const int nTriggers = 6;
static const std::vector<std::string> triggerNames{"fPPP", "fPPL", "fPLL", "fLLL", "fPD", "fLD"};
static const float triggerSwitches[1][nTriggers]{
  {0, 0, 0, 0, 0, 0}};
} // namespace software_triggers

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

struct FemtoDreamProducerTaskReso {

  SliceCache cache;                                                        // o2::framework, included in ASoAHelpers.h
  Preslice<aod::FemtoFullTracks> perCol = aod::track::collisionId;         // o2::framework included in ASoAHelpers.h
  Partition<aod::FemtoFullTracks> daughter1 = aod::track::signed1Pt > 0.f; // o2::framework included in AnalysisHelper.h
  Partition<aod::FemtoFullTracks> daughter2 = aod::track::signed1Pt < 0.f; // o2::framework included in AnalysisHelper.h

  Zorro zorro;

  Produces<aod::FDCollisions> outputCollision;
  Produces<aod::FDMCCollisions> outputMCCollision;
  Produces<aod::FDMCCollLabels> outputCollsMCLabels;
  Produces<aod::FDParticles> outputParts;
  Produces<aod::FDMCParticles> outputPartsMC;
  Produces<aod::FDExtParticles> outputDebugParts;
  Produces<aod::FDMCLabels> outputPartsMCLabels;
  Produces<aod::FDExtMCParticles> outputDebugPartsMC;
  Produces<aod::FDExtMCLabels> outputPartsExtMCLabels;

  Configurable<bool> confIsDebug{"confIsDebug", true, "Enable Debug tables"};
  Configurable<bool> confUseItsPid{"confUseItsPid", false, "Enable Debug tables"};
  Configurable<bool> confIsRun3{"confIsRun3", false, "Running on Run3 or pilot"}; // true?
  Configurable<bool> confIsForceGRP{"confIsForceGRP", false, "Set true if the magnetic field configuration is not available in the usual CCDB directory (e.g. for Run 2 converted data or unanchorad Monte Carlo)"};
  /// Event cuts
  FemtoDreamCollisionSelection colCuts;
  // Event cuts - Triggers
  Configurable<bool> confEnableTriggerSelection{"confEnableTriggerSelection", false, "Should the trigger selection be enabled for collisions?"};
  Configurable<LabeledArray<float>> confTriggerSwitches{"confTriggerSwitches", {software_triggers::triggerSwitches[0], 1, software_triggers::nTriggers, std::vector<std::string>{"Switch"}, software_triggers::triggerNames}, "Turn on which trigger should be checked for recorded events to pass selection"};
  Configurable<std::string> confBaseCCDBPathForTriggers{"confBaseCCDBPathForTriggers", "Users/m/mpuccio/EventFiltering/OTS/Chunked/", "Provide ccdb path for trigger table; default - trigger coordination"};

  // Event cuts - usual selection criteria
  Configurable<float> confEvtZvtx{"confEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};
  Configurable<bool> confEvtTriggerCheck{"confEvtTriggerCheck", true, "Evt sel: check for trigger"};
  Configurable<int> confEvtTriggerSel{"confEvtTriggerSel", kINT7, "Evt sel: trigger"};
  Configurable<bool> confEvtOfflineCheck{"confEvtOfflineCheck", false, "Evt sel: check for offline selection"};
  Configurable<bool> confEvtAddOfflineCheck{"confEvtAddOfflineCheck", false, "Evt sel: additional checks for offline selection (not part of sel8 yet)"};
  Configurable<bool> confIsActivateV0{"confIsActivateV0", true, "Activate filling of V0 into femtodream tables"};
  Configurable<bool> confIsActivateReso{"confIsActivateReso", true, "Activate filling of sl Resonances into femtodream tables"};
  Configurable<bool> confIsActivatePhi{"confIsActivatePhi", true, "Activates cuts on Phi's and fills tables"};
  Configurable<bool> confIsActivateCascade{"confIsActivateCascade", false, "Activate filling of Cascades into femtodream tables"};
  Configurable<float> confEvtMinSphericity{"confEvtMinSphericity", 0.0f, "Evt sel: Min. sphericity of event"};
  Configurable<float> confEvtSphericityPtmin{"confEvtSphericityPtmin", 0.0f, "Evt sel: Min. Pt for sphericity calculation"};

  Configurable<bool> confTrkRejectNotPropagated{"confTrkRejectNotPropagated", false, "True: reject not propagated tracks"};
  // Configurable<bool> confRejectITSHitandTOFMissing{ "confRejectITSHitandTOFMissing", false, "True: reject if neither ITS hit nor TOF timing satisfied"};
  Configurable<int> confTrkPDGCode{"confTrkPDGCode", 2212, "PDG code of the selected track for Monte Carlo truth"};
  FemtoDreamTrackSelection trackCuts;
  struct : ConfigurableGroup {
    std::string prefix = std::string("Track");
    Configurable<std::vector<float>> confTrkCharge{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kSign, "confTrk"), std::vector<float>{-1, 1}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kSign, "Track selection: ")};
    Configurable<std::vector<float>> confTrkPtmin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kpTMin, "confTrk"), std::vector<float>{0.1f, 0.15f, 0.2f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kpTMin, "Track selection: ")};
    Configurable<std::vector<float>> confTrkPtmax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kpTMax, "confTrk"), std::vector<float>{4.4f, 4.6f, 4.5f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kpTMax, "Track selection: ")};
    Configurable<std::vector<float>> confTrkEta{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kEtaMax, "confTrk"), std::vector<float>{0.8f, 0.85f, 0.9f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kEtaMax, "Track selection: ")};
    Configurable<std::vector<float>> confTrkTPCnclsMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCnClsMin, "confTrk"), std::vector<float>{80.f, 90.f, 100.f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCnClsMin, "Track selection: ")};
    Configurable<std::vector<float>> confTrkTPCfCls{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCfClsMin, "confTrk"), std::vector<float>{0.7f, 0.83f, 0.9f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCfClsMin, "Track selection: ")};
    Configurable<std::vector<float>> confTrkTPCcRowsMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCcRowsMin, "confTrk"), std::vector<float>{70.f, 60.f, 80.f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCcRowsMin, "Track selection: ")};
    Configurable<std::vector<float>> confTrkTPCsCls{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCsClsMax, "confTrk"), std::vector<float>{0.1f, 160.f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCsClsMax, "Track selection: ")};
    Configurable<std::vector<float>> confTrkITSnclsMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kITSnClsMin, "confTrk"), std::vector<float>{-1.f, 2.f, 4.f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kITSnClsMin, "Track selection: ")};
    Configurable<std::vector<float>> confTrkITSnclsIbMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kITSnClsIbMin, "confTrk"), std::vector<float>{-1.f, 1.f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kITSnClsIbMin, "Track selection: ")};
    Configurable<std::vector<float>> confTrkDCAxyMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kDCAxyMax, "confTrk"), std::vector<float>{0.2f, 0.5f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kDCAxyMax, "Track selection: ")};
    Configurable<std::vector<float>> confTrkDCAzMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kDCAzMax, "confTrk"), std::vector<float>{0.2f, 0.5f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kDCAzMax, "Track selection: ")};
    Configurable<std::vector<float>> confTrkPIDnSigmaMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kPIDnSigmaMax, "confTrk"), std::vector<float>{3.5f, 3.f, 2.5f}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kPIDnSigmaMax, "Track selection: ")};
    Configurable<float> confTrkPIDnSigmaOffsetTPC{"confTrkPIDnSigmaOffsetTPC", 0., "Offset for TPC nSigma because of bad calibration"}; // set to zero for run3 or so
    Configurable<float> confTrkPIDnSigmaOffsetTOF{"confTrkPIDnSigmaOffsetTOF", 0., "Offset for TOF nSigma because of bad calibration"};
    Configurable<std::vector<int>> confTrkPIDspecies{"confTrkPIDspecies", std::vector<int>{o2::track::PID::Pion, o2::track::PID::Kaon, o2::track::PID::Proton, o2::track::PID::Deuteron}, "Trk sel: Particles species for PID"};
    // missing DCA Configurable?? because implemented in TrackSelection.h
  } Track;

  FemtoDreamV0Selection v0Cuts;
  Configurable<std::vector<float>> confV0Sign{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0Sign, "confV0"), std::vector<float>{-1, 1}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0Sign, "V0 selection: ")};
  Configurable<std::vector<float>> confV0PtMin{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0pTMin, "confV0"), std::vector<float>{0.3f, 0.4f, 0.5f}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0pTMin, "V0 selection: ")};
  Configurable<std::vector<float>> confV0PtMax{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0pTMax, "confV0"), std::vector<float>{3.3f, 3.4f, 3.5f}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0pTMax, "V0 selection: ")};
  Configurable<std::vector<float>> confV0EtaMax{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0etaMax, "confV0"), std::vector<float>{0.8f, 0.7f, 0.9f}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0etaMax, "V0 selection: ")};
  Configurable<std::vector<float>> confV0DCADaughMax{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0DCADaughMax, "confV0"), std::vector<float>{1.2f, 1.5f}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0DCADaughMax, "V0 selection: ")};
  Configurable<std::vector<float>> confV0CPAMin{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0CPAMin, "confV0"), std::vector<float>{0.99f, 0.995f}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0CPAMin, "V0 selection: ")};
  Configurable<std::vector<float>> confV0TranRadMin{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0TranRadMin, "confV0"), std::vector<float>{0.2f}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0TranRadMin, "V0 selection: ")};
  Configurable<std::vector<float>> confV0TranRadMax{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0TranRadMax, "confV0"), std::vector<float>{100.f}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0TranRadMax, "V0 selection: ")};
  Configurable<std::vector<float>> confV0DecVtxMax{FemtoDreamV0Selection::getSelectionName(femto_dream_v0_selection::kV0DecVtxMax, "confV0"), std::vector<float>{100.f}, FemtoDreamV0Selection::getSelectionHelper(femto_dream_v0_selection::kV0DecVtxMax, "V0 selection: ")};

  Configurable<float> confV0InvMassLowLimit{"confV0InvMassLowLimit", 1.05, "Lower limit of the V0 invariant mass"};
  Configurable<float> confV0InvMassUpLimit{"confV0InvMassUpLimit", 1.30, "Upper limit of the V0 invariant mass"};
  Configurable<bool> confV0RejectKaons{"confV0RejectKaons", false, "Switch to reject kaons"};
  Configurable<bool> confV0RejectLambdas{"confV0RejectLambdas", false, "Switch to reject lambdas (if mother is kaon)"};
  Configurable<float> confV0InvKaonMassLowLimit{"confV0InvKaonMassLowLimit", 0.48, "Lower limit of the V0 invariant mass for Kaon rejection"};
  Configurable<float> confV0InvKaonMassUpLimit{"confV0InvKaonMassUpLimit", 0.515, "Upper limit of the V0 invariant mass for Kaon rejection"};
  Configurable<bool> confV0MotherIsLambda{"confV0MotherIsLambda", true, "True: Lambda, False: K0Short"};

  Configurable<std::vector<float>> confChildSign{"confChildSign", std::vector<float>{-1, 1}, "V0 Child sel: Charge"};
  Configurable<std::vector<float>> confChildEtaMax{"confChildEtaMax", std::vector<float>{0.8f}, "V0 Child sel: max eta"};
  Configurable<std::vector<float>> confChildTPCnClsMin{"confChildTPCnClsMin", std::vector<float>{80.f, 70.f, 60.f}, "V0 Child sel: Min. nCls TPC"};
  Configurable<std::vector<float>> confChildDCAMin{"confChildDCAMin", std::vector<float>{0.05f, 0.06f}, "V0 Child sel:  Max. DCA Daugh to PV (cm)"};
  Configurable<std::vector<float>> confChildPIDnSigmaMax{"confChildPIDnSigmaMax", std::vector<float>{5.f, 4.f}, "V0 Child sel: Max. PID nSigma TPC"};
  Configurable<std::vector<int>> confChildPIDspecies{"confChildPIDspecies", std::vector<int>{o2::track::PID::Pion, o2::track::PID::Proton}, "V0 Child sel: Particles species for PID"};

  FemtoDreamCascadeSelection cascadeCuts;
  struct : o2::framework::ConfigurableGroup {
    Configurable<float> confCascInvMassLowLimit{"confCascInvMassLowLimit", 1.2, "Lower limit of the Cascade invariant mass"};
    Configurable<float> confCascInvMassUpLimit{"confCascInvMassUpLimit", 1.5, "Upper limit of the Cascade invariant mass"};
    Configurable<bool> confCascIsSelectedOmega{"confCascIsSelectedOmega", false, "Select Omegas instead of Xis (invariant mass)"};
    // Cascade
    Configurable<std::vector<float>> confCascadeSign{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadeSign, "confCascade"), std::vector<float>{-1, 1}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadeSign, "Cascade selection: ")};
    Configurable<std::vector<float>> confCascadePtMin{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadePtMin, "confCascade"), std::vector<float>{0.3f, 0.4f}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadePtMin, "Cascade selection: ")};
    Configurable<std::vector<float>> confCascadePtMax{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadePtMax, "confCascade"), std::vector<float>{5.5f, 6.0f}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadePtMax, "Cascade selection: ")};
    Configurable<std::vector<float>> confCascadeEtaMax{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadeEtaMax, "confCascade"), std::vector<float>{0.8f}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadeEtaMax, "Cascade selection: ")};
    Configurable<std::vector<float>> confCascadeDCADaughMax{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadeDCADaughMax, "confCascade"), std::vector<float>{1.f, 1.2f}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadeDCADaughMax, "Cascade selection: ")};
    Configurable<std::vector<float>> confCascadeCPAMin{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadeCPAMin, "confCascade"), std::vector<float>{0.99f, 0.95f}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadeCPAMin, "Cascade selection: ")};
    Configurable<std::vector<float>> confCascadeTranRadMin{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadeTranRadMin, "confCascade"), std::vector<float>{0.2f, 0.5f}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadeTranRadMin, "Cascade selection: ")};
    Configurable<std::vector<float>> confCascadeTranRadMax{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadeTranRadMax, "confCascade"), std::vector<float>{100.f}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadeTranRadMax, "Cascade selection: ")};
    Configurable<std::vector<float>> confCascadeDecVtxMax{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadeDecVtxMax, "confCascade"), std::vector<float>{100.f}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadeDecVtxMax, "Cascade selection: ")};

    // Cascade v0 daughters
    Configurable<std::vector<float>> confCascadeV0DCADaughMax{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadeV0DCADaughMax, "confCascade"), std::vector<float>{1.2f, 1.5f}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadeV0DCADaughMax, "CascV0 selection: ")};
    Configurable<std::vector<float>> confCascadeV0CPAMin{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadeV0CPAMin, "confCascade"), std::vector<float>{0.99f, 0.995f}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadeV0CPAMin, "CascV0 selection: ")};
    Configurable<std::vector<float>> confCascadeV0TranRadMin{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadeV0TranRadMin, "confCascade"), std::vector<float>{0.2f}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadeV0TranRadMin, "CascV0 selection: ")};
    Configurable<std::vector<float>> confCascadeV0TranRadMax{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadeV0TranRadMax, "confCascade"), std::vector<float>{100.f}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadeV0TranRadMax, "CascV0 selection: ")};
    Configurable<std::vector<float>> confCascadeV0DCAtoPVMin{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadeV0DCAtoPVMin, "confCascade"), std::vector<float>{100.f}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadeV0DCAtoPVMin, "CascV0 selection: ")};
    Configurable<std::vector<float>> confCascadeV0DCAtoPVMax{FemtoDreamCascadeSelection::getSelectionName(femtoDreamCascadeSelection::kCascadeV0DCAtoPVMax, "confCascade"), std::vector<float>{100.f}, FemtoDreamCascadeSelection::getSelectionHelper(femtoDreamCascadeSelection::kCascadeV0DCAtoPVMax, "CascV0 selection: ")};
    Configurable<float> confCascV0InvMassLowLimit{"confCascV0InvMassLowLimit", 1.011461, "Lower limit of the Cascade invariant mass"};
    Configurable<float> confCascV0InvMassUpLimit{"confCascV0InvMassUpLimit", 1.027461, "Upper limit of the Cascade invariant mass"};
    // Cascade Daughter Tracks
    Configurable<std::vector<float>> confCascV0ChildSign{"confCascV0ChildSign", std::vector<float>{-1, 1}, "CascV0 Child sel: Charge"};
    Configurable<std::vector<float>> confCascV0ChildPtMin{"confCascV0ChildPtMin", std::vector<float>{0.8f}, "CascV0 Child sel: min pt"};
    Configurable<std::vector<float>> confCascV0ChildEtaMax{"confCascV0ChildEtaMax", std::vector<float>{0.8f}, "CascV0 Child sel: max eta"};
    Configurable<std::vector<float>> confCascV0ChildTPCnClsMin{"confCascV0ChildTPCnClsMin", std::vector<float>{80.f, 70.f, 60.f}, "CascV0 Child sel: Min. nCls TPC"};
    Configurable<std::vector<float>> confCascV0ChildDCAMin{"confCascV0ChildDCAMin", std::vector<float>{0.05f, 0.06f}, "CascV0 Child sel:  Max. DCA Daugh to PV (cm)"};
    Configurable<std::vector<float>> confCascV0ChildPIDnSigmaMax{"confCascV0ChildPIDnSigmaMax", std::vector<float>{5.f, 4.f}, "CascV0 Child sel: Max. PID nSigma TPC"};
    Configurable<std::vector<int>> confCascV0ChildPIDspecies{"confCascV0ChildPIDspecies", std::vector<int>{o2::track::PID::Pion, o2::track::PID::Proton}, "CascV0 Child sel: Particles species for PID"};
    // Cascade Bachelor Track
    Configurable<std::vector<float>> confCascBachelorSign{"confCascBachelorSign", std::vector<float>{-1, 1}, "Cascade Bachelor sel: Charge"};
    Configurable<std::vector<float>> confCascBachelorPtMin{"confCascBachelorPtMin", std::vector<float>{0.8f}, "Cascade Bachelor sel: min pt"};
    Configurable<std::vector<float>> confCascBachelorEtaMax{"confCascBachelorEtaMax", std::vector<float>{0.8f}, "Cascade Bachelor sel: max eta"};
    Configurable<std::vector<float>> confCascBachelorTPCnClsMin{"confCascBachelorTPCnClsMin", std::vector<float>{80.f, 70.f, 60.f}, "Cascade Bachelor sel: Min. nCls TPC"};
    Configurable<std::vector<float>> confCascBachelorDCAMin{"confCascBachelorDCAMin", std::vector<float>{0.05f, 0.06f}, "Cascade Bachelor sel:  Max. DCA Daugh to PV (cm)"};
    Configurable<std::vector<float>> confCascBachelorPIDnSigmaMax{"confCascBachelorPIDnSigmaMax", std::vector<float>{5.f, 4.f}, "Cascade Bachelor sel: Max. PID nSigma TPC"};
    Configurable<std::vector<int>> confCascBachelorPIDspecies{"confCascBachelorPIDspecies", std::vector<int>{o2::track::PID::Pion}, "Cascade Bachelor sel: Particles species for PID"};

    Configurable<bool> confCascRejectCompetingMass{"confCascRejectCompetingMass", false, "Switch on to reject Omegas (for Xi) or Xis (for Omegas)"};
    Configurable<float> confCascInvCompetingMassLowLimit{"confCascInvCompetingMassLowLimit", 1.66, "Lower limit of the cascade invariant mass for competing mass rejection"};
    Configurable<float> confCascInvCompetingMassUpLimit{"confCascInvCompetingMassUpLimit", 1.68, "Upper limit of the cascade invariant mass for competing mass rejection"};

  } confCascSel;

  // Resonances
  FemtoDreamResoSelection resoCuts;
  struct : ConfigurableGroup {
    std::string prefix = std::string("Resonance");

    Configurable<float> confResoInvMass{"confResoInvMass", o2::constants::physics::MassPhi, "Literature value of the Reso invariant mass"};
    Configurable<float> confResoInvMassLowLimit{"confResoInvMassLowLimit", 0.9, "Lower limit of the Reso invariant mass"}; // 1.011461
    Configurable<float> confResoInvMassUpLimit{"confResoInvMassUpLimit", 1.15, "Upper limit of the Reso invariant mass"};  // 1.027461
    Configurable<bool> confUseMassDiff{"confUseMassDiff", false, "(De)activates the usage of invMassDiff when checking combinations"};
    Configurable<bool> confDoLikeSign{"confDoLikeSign", false, "(De)activates likeSign histo filling"};
    Configurable<bool> confUseMassDiffLikeSign{"confUseMassDiffLikeSign", false, "(De)activates the usage of invMassDiff when checking combinations in LikeSign"};
    Configurable<int> confMassQAPart2PID{"confMassQAPart2PID", o2::track::PID::Pion, "Daughter PID for massQAPart2 histograms"};
    Configurable<std::vector<float>> confResoSign{"confResoSign", std::vector<float>{-1., 1.}, "Reso Sign selection"};

    Configurable<std::vector<float>> confDaughterCharge{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kSign, "confDaughter"), std::vector<float>{-1, 1}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kSign, "Reso selection: ")};
    Configurable<std::vector<float>> confDaughterPtMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kpTMin, "confDaughter"), std::vector<float>{0.1, 0.15, 0.2}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kpTMin, "Reso selection: ")};
    Configurable<std::vector<float>> confDaughterPtMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kpTMax, "confDaughter"), std::vector<float>{5.0, 4.0}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kpTMax, "Reso selection: ")};
    Configurable<std::vector<float>> confDaughterEtaMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kEtaMax, "confDaughter"), std::vector<float>{0.8, 0.85, 0.9}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kEtaMax, "Reso selection: ")};
    Configurable<std::vector<float>> confDaughterTPCnClsMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCnClsMin, "confDaughter"), std::vector<float>{75, 85, 100}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCnClsMin, "Reso selection: ")};
    Configurable<std::vector<float>> confDaughterTPCfClsMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCfClsMin, "confDaughter"), std::vector<float>{0.7, 0.8, 0.9}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCfClsMin, "Reso selection: ")};
    Configurable<std::vector<float>> confDaughterTPCcRowsMin{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kTPCcRowsMin, "confDaughter"), std::vector<float>{75, 85, 100}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kTPCcRowsMin, "Reso selection: ")};
    Configurable<std::vector<float>> confDaughterDCAxyMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kDCAxyMax, "confDaughter"), std::vector<float>{0.2, 0.15, 0.1}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kDCAxyMax, "Reso selection: ")};
    Configurable<std::vector<float>> confDaughterDCAzMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kDCAzMax, "confDaughter"), std::vector<float>{0.2, 0.15, 0.1}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kDCAzMax, "Reso selection: ")};
    Configurable<std::vector<float>> confDaughterPIDnSigmaMax{FemtoDreamTrackSelection::getSelectionName(femtoDreamTrackSelection::kPIDnSigmaMax, "confDaughter"), std::vector<float>{3.0, 2.5, 2.0}, FemtoDreamTrackSelection::getSelectionHelper(femtoDreamTrackSelection::kPIDnSigmaMax, "Reso selection: ")};
    // Configurable<float> confResoMassUp{"confResoMassUp", 0.52, "Upper limit for the mass selection of the daughters"};
    // Configurable<float> confResoMassLow{"confResoMassLow", 0.48, "Lower limit for the mass selection of the daughters"};
    Configurable<std::vector<int>> confDaughterPIDspecies{"confDaughterPIDspecies", std::vector<int>{o2::track::PID::Kaon, o2::track::PID::Kaon}, "Reso Daughter sel: Particles species for PID"};
    Configurable<std::vector<float>> confDaughterPTPCThr{"confDaughterPTPCThr", std::vector<float>{0.5, 0.5}, "p_T (GeV/c) & pid dependent Thresholds for case distinction between TPC/TPCTOF"};
  } Resonance;

  /// \todo should we add filter on min value pT/eta of V0 and daughters?
  /*Filter v0Filter = (nabs(aod::v0data::x) < V0DecVtxMax.value) &&
                    (nabs(aod::v0data::y) < V0DecVtxMax.value) &&
                    (nabs(aod::v0data::z) < V0DecVtxMax.value);*/
  // (aod::v0data::v0radius > V0TranRadV0Min.value); to be added, not working
  // for now do not know why

  /// General options
  struct : o2::framework::ConfigurableGroup {
    Configurable<float> confTrkMinChi2PerClusterTPC{"confTrkMinChi2PerClusterTPC", 0.f, "Lower limit for chi2 of TPC; currently for testing only"};
    Configurable<float> confTrkMaxChi2PerClusterTPC{"confTrkMaxChi2PerClusterTPC", 1000.f, "Upper limit for chi2 of TPC; currently for testing only"};
    Configurable<float> confTrkMaxChi2PerClusterITS{"confTrkMaxChi2PerClusterITS", 1000.0f, "Minimal track selection: max allowed chi2 per ITS cluster"}; // 36.0 is default
    Configurable<bool> confTrkTPCRefit{"confTrkTPCRefit", false, "True: require TPC refit"};
    Configurable<bool> confTrkITSRefit{"confTrkITSRefit", false, "True: require ITS refit"};

  } OptionTrackSpecialSelections;

  struct : o2::framework::ConfigurableGroup {
    Configurable<bool> requireRCTFlagChecker{"requireRCTFlagChecker", true, "Check event quality in run condition table"};
    Configurable<std::string> cfgEvtRCTFlagCheckerLabel{"cfgEvtRCTFlagCheckerLabel", "CBT_hadronPID", "Evt sel: RCT flag checker label"};
    Configurable<bool> cfgEvtRCTFlagCheckerLimitAcceptAsBad{"cfgEvtRCTFlagCheckerLimitAcceptAsBad", true, "Evt sel: RCT flag checker treat Limited Acceptance As Bad"};
  } rctCut;

  HistogramRegistry qaRegistry{"QAHistos", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry trackRegistry{"Tracks", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry v0Registry{"V0", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry cascadeRegistry{"Cascade", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resoRegistry{"Reso", {}, OutputObjHandlingPolicy::AnalysisObject};

  int mRunNumber;
  float mMagField;
  std::string zorroTriggerNames = "";
  Service<o2::ccdb::BasicCCDBManager> ccdb; /// Accessing the CCDB
  RCTFlagsChecker rctChecker;

  void init(InitContext&)
  {
    if (doprocessData == false && doprocessData_noCentrality == false && doprocessDataCentPbPb == false && doprocessMC == false && doprocessMCnoCentrality == false && doprocessMCCentPbPb == false) {
      LOGF(fatal, "Neither processData nor processMC enabled. Please choose one.");
    }
    if ((doprocessData == true && doprocessMC == true) || (doprocessData == true && doprocessMCnoCentrality == true) || (doprocessMC == true && doprocessMCnoCentrality == true) || (doprocessData_noCentrality == true && doprocessData == true) || (doprocessData_noCentrality == true && doprocessMC == true) || (doprocessData_noCentrality == true && doprocessMCnoCentrality == true) || (doprocessDataCentPbPb == true && doprocessData == true) || (doprocessDataCentPbPb == true && doprocessData_noCentrality == true) || (doprocessDataCentPbPb == true && doprocessMC == true) || (doprocessDataCentPbPb == true && doprocessMCnoCentrality == true) || (doprocessDataCentPbPb == true && doprocessMCCentPbPb == true)) {
      LOGF(fatal,
           "Cannot enable more than one process switch at the same time. "
           "Please choose one.");
    }

    int cutBits = 8 * sizeof(o2::aod::femtodreamparticle::cutContainerType);
    trackRegistry.add("AnalysisQA/CutCounter", "; Bit; Counter", kTH1F, {{cutBits + 1, -0.5, cutBits + 0.5}});
    trackRegistry.add("AnalysisQA/Chi2ITSTPCperCluster", "; ITS_Chi2; TPC_Chi2", kTH2F, {{100, 0, 50}, {100, 0, 20}});
    trackRegistry.add("AnalysisQA/RefitITSTPC", "; ITS_Refit; TPC_Refit", kTH2F, {{2, 0, 2}, {2, 0, 2}});
    trackRegistry.add("AnalysisQA/getGenStatusCode", "; Bit; Entries", kTH1F, {{200, 0, 200}});
    trackRegistry.add("AnalysisQA/getProcess", "; Bit; Entries", kTH1F, {{200, 0, 200}});
    trackRegistry.add("AnalysisQA/Mother", "; Bit; Entries", kTH1F, {{4000, -4000, 4000}});
    trackRegistry.add("AnalysisQA/Particle", "; Bit; Entries", kTH1F, {{4000, -4000, 4000}});
    v0Registry.add("AnalysisQA/CutCounter", "; Bit; Counter", kTH1F, {{cutBits + 1, -0.5, cutBits + 0.5}});

    resoRegistry.add("AnalysisQA/Reso/InvMass", "Invariant mass V0s;M_{KK};Entries", HistType::kTH1F, {{7000, 0.65, 1.5}});
    resoRegistry.add("AnalysisQA/Reso/InvMassAnti", "Invariant mass V0s;M_{KK};Entries", HistType::kTH1F, {{7000, 0.65, 1.5}});
    resoRegistry.add("AnalysisQA/Reso/InvMass_phi_selected", "Invariant mass V0s;M_{KK};Entries", HistType::kTH1F, {{7000, 0.65, 1.5}});
    resoRegistry.add("AnalysisQA/Reso/InvMassAnti_phi_selected", "Invariant mass V0s;M_{KK};Entries", HistType::kTH1F, {{7000, 0.65, 1.5}});

    // Daughter-histos for normal resonance
    // Histos for PosDaughter
    resoRegistry.add("AnalysisQA/Reso/ResoQA/PosDaughter/Pt", "Transverse momentum of all tracks;p_{T} (GeV/c);Entries", HistType::kTH1F, {{1000, 0, 10}});
    resoRegistry.add("AnalysisQA/Reso/ResoQA/PosDaughter/Eta", "Pseudorapidity of all  tracks;#eta;Entries", HistType::kTH1F, {{1000, -2, 2}});
    resoRegistry.add("AnalysisQA/Reso/ResoQA/PosDaughter/Phi", "Azimuthal angle of all  tracks;#phi;Entries", HistType::kTH1F, {{720, 0, o2::constants::math::TwoPI}});
    resoRegistry.add("AnalysisQA/Reso/ResoQA/PosDaughter/DcaXY", "dcaXY of all  tracks;d_{XY} (cm);Entries", HistType::kTH1F, {{1000, 0, 1}}); // check if cm is correct here
    resoRegistry.add("AnalysisQA/Reso/ResoQA/PosDaughter/DcaZ", "dcaZ of all  tracks;d_{Z} (cm);Entries", HistType::kTH1F, {{1000, 0, 1}});    // check if cm is correct here
    // Histos for NegDaughter
    resoRegistry.add("AnalysisQA/Reso/ResoQA/NegDaughter/Pt", "Transverse momentum of all tracks;p_{T} (GeV/c);Entries", HistType::kTH1F, {{1000, 0, 10}});
    resoRegistry.add("AnalysisQA/Reso/ResoQA/NegDaughter/Eta", "Pseudorapidity of all  tracks;#eta;Entries", HistType::kTH1F, {{1000, -2, 2}});
    resoRegistry.add("AnalysisQA/Reso/ResoQA/NegDaughter/Phi", "Azimuthal angle of all  tracks;#phi;Entries", HistType::kTH1F, {{720, 0, o2::constants::math::TwoPI}});
    resoRegistry.add("AnalysisQA/Reso/ResoQA/NegDaughter/DcaXY", "dcaXY of all  tracks;d_{XY} (cm);Entries", HistType::kTH1F, {{1000, 0, 1}}); // check if cm is correct here
    resoRegistry.add("AnalysisQA/Reso/ResoQA/NegDaughter/DcaZ", "dcaZ of all  tracks;d_{Z} (cm);Entries", HistType::kTH1F, {{1000, 0, 1}});    // check if cm is correct here

    // Daughter-histos for anti resonance
    // Histos for PosDaughter
    resoRegistry.add("AnalysisQA/Reso/AntiResoQA/PosDaughter/Pt", "Transverse momentum of all tracks;p_{T} (GeV/c);Entries", HistType::kTH1F, {{1000, 0, 10}});
    resoRegistry.add("AnalysisQA/Reso/AntiResoQA/PosDaughter/Eta", "Pseudorapidity of all  tracks;#eta;Entries", HistType::kTH1F, {{1000, -2, 2}});
    resoRegistry.add("AnalysisQA/Reso/AntiResoQA/PosDaughter/Phi", "Azimuthal angle of all  tracks;#phi;Entries", HistType::kTH1F, {{720, 0, o2::constants::math::TwoPI}});
    resoRegistry.add("AnalysisQA/Reso/AntiResoQA/PosDaughter/DcaXY", "dcaXY of all  tracks;d_{XY} (cm);Entries", HistType::kTH1F, {{1000, 0, 1}}); // check if cm is correct here
    resoRegistry.add("AnalysisQA/Reso/AntiResoQA/PosDaughter/DcaZ", "dcaZ of all  tracks;d_{Z} (cm);Entries", HistType::kTH1F, {{1000, 0, 1}});    // check if cm is correct here
    // Histos for NegDaughter
    resoRegistry.add("AnalysisQA/Reso/AntiResoQA/NegDaughter/Pt", "Transverse momentum of all tracks;p_{T} (GeV/c);Entries", HistType::kTH1F, {{1000, 0, 10}});
    resoRegistry.add("AnalysisQA/Reso/AntiResoQA/NegDaughter/Eta", "Pseudorapidity of all  tracks;#eta;Entries", HistType::kTH1F, {{1000, -2, 2}});
    resoRegistry.add("AnalysisQA/Reso/AntiResoQA/NegDaughter/Phi", "Azimuthal angle of all  tracks;#phi;Entries", HistType::kTH1F, {{720, 0, o2::constants::math::TwoPI}});
    resoRegistry.add("AnalysisQA/Reso/AntiResoQA/NegDaughter/DcaXY", "dcaXY of all  tracks;d_{XY} (cm);Entries", HistType::kTH1F, {{1000, 0, 1}}); // check if cm is correct here
    resoRegistry.add("AnalysisQA/Reso/AntiResoQA/NegDaughter/DcaZ", "dcaZ of all  tracks;d_{Z} (cm);Entries", HistType::kTH1F, {{1000, 0, 1}});    // check if cm is correct here

    // SelectionQA
    resoRegistry.add("AnalysisQA/Reso/ResoQA/InvMassSwitched", "Invariant mass V0s;M_{KK};Entries", HistType::kTH1F, {{7000, 0.65, 1.5}}); // for opposite mass hypothesis (Reso is anti)
    resoRegistry.add("AnalysisQA/Reso/ResoQA/InvMassBothPID1", "Invariant mass V0s;M_{KK};Entries", HistType::kTH1F, {{7000, 0.65, 1.5}}); // both particles are of type confDaughterPIDspecies[0]
    resoRegistry.add("AnalysisQA/Reso/ResoQA/InvMassBothPID2", "Invariant mass V0s;M_{KK};Entries", HistType::kTH1F, {{7000, 0.65, 1.5}}); // both particles are of type confDaughterPIDspecies[1]

    // AntiSelectionQA
    resoRegistry.add("AnalysisQA/Reso/AntiResoQA/InvMassSwitched", "Invariant mass V0s;M_{KK};Entries", HistType::kTH1F, {{7000, 0.65, 1.5}}); // for opposite mass hypothesis (Reso is normal)
    resoRegistry.add("AnalysisQA/Reso/AntiResoQA/InvMassBothPID1", "Invariant mass V0s;M_{KK};Entries", HistType::kTH1F, {{7000, 0.65, 1.5}}); // both particles are of type confDaughterPIDspecies[0]
    resoRegistry.add("AnalysisQA/Reso/AntiResoQA/InvMassBothPID2", "Invariant mass V0s;M_{KK};Entries", HistType::kTH1F, {{7000, 0.65, 1.5}}); // both particles are of type confDaughterPIDspecies[1]

    // likeSign
    if (Resonance.confDoLikeSign.value) {
      resoRegistry.add("AnalysisQA/ResoLikeSign/ResoQA/InvMass", "Invariant mass V0s;M_{KK};Entries", HistType::kTH1F, {{7000, 0.65, 1.5}});
      resoRegistry.add("AnalysisQA/ResoLikeSign/ResoQA/InvMassSwitched", "Invariant mass V0s;M_{KK};Entries", HistType::kTH1F, {{7000, 0.65, 1.5}}); // for opposite mass hypothesis (Reso is anti)
      resoRegistry.add("AnalysisQA/ResoLikeSign/ResoQA/InvMassBothPID1", "Invariant mass V0s;M_{KK};Entries", HistType::kTH1F, {{7000, 0.65, 1.5}}); // both particles are of type confDaughterPIDspecies[0]
      resoRegistry.add("AnalysisQA/ResoLikeSign/ResoQA/InvMassBothPID2", "Invariant mass V0s;M_{KK};Entries", HistType::kTH1F, {{7000, 0.65, 1.5}}); // both particles are of type confDaughterPIDspecies[1]

      resoRegistry.add("AnalysisQA/ResoLikeSign/AntiResoQA/InvMass", "Invariant mass V0s;M_{KK};Entries", HistType::kTH1F, {{7000, 0.65, 1.5}});
      resoRegistry.add("AnalysisQA/ResoLikeSign/AntiResoQA/InvMassSwitched", "Invariant mass V0s;M_{KK};Entries", HistType::kTH1F, {{7000, 0.65, 1.5}}); // for opposite mass hypothesis (Reso is anti)
      resoRegistry.add("AnalysisQA/ResoLikeSign/AntiResoQA/InvMassBothPID1", "Invariant mass V0s;M_{KK};Entries", HistType::kTH1F, {{7000, 0.65, 1.5}}); // both particles are of type confDaughterPIDspecies[0]
      resoRegistry.add("AnalysisQA/ResoLikeSign/AntiResoQA/InvMassBothPID2", "Invariant mass V0s;M_{KK};Entries", HistType::kTH1F, {{7000, 0.65, 1.5}}); // both particles are of type confDaughterPIDspecies[1]
    }

    resoRegistry.add("AnalysisQA/Reso/Pt_posdaughter_selected", "Transverse momentum of all processed tracks;p_{T} (GeV/c);Entries", HistType::kTH1F, {{1000, 0, 10}});
    resoRegistry.add("AnalysisQA/Reso/Eta_posdaughter_selected", "Pseudorapidity of all processed tracks;#eta;Entries", HistType::kTH1F, {{1000, -2, 2}});
    resoRegistry.add("AnalysisQA/Reso/Phi_posdaughter_selected", "Azimuthal angle of all processed tracks;#phi;Entries", HistType::kTH1F, {{720, 0, o2::constants::math::TwoPI}});
    resoRegistry.add("AnalysisQA/Reso/DCAxy_posdaughter_selected", "dcaXY of all processed tracks;d_{XY} (cm);Entries", HistType::kTH1F, {{1000, 0, 1}}); // check if cm is correct here
    resoRegistry.add("AnalysisQA/Reso/DCAz_posdaughter_selected", "dcaZ of all processed tracks;d_{Z} (cm);Entries", HistType::kTH1F, {{1000, 0, 1}});    // check if cm is correct here
    resoRegistry.add("AnalysisQA/Reso/Eta_negdaughter_selected", "Pseudorapidity of all processed tracks;#eta;Entries", HistType::kTH1F, {{1000, -2, 2}});
    resoRegistry.add("AnalysisQA/Reso/Phi_negdaughter_selected", "Azimuthal angle of all processed tracks;#phi;Entries", HistType::kTH1F, {{720, 0, o2::constants::math::TwoPI}});
    resoRegistry.add("AnalysisQA/Reso/Pt_negdaughter_selected", "Transverse momentum of all processed tracks;p_{T} (GeV/c);Entries", HistType::kTH1F, {{1000, 0, 10}});
    resoRegistry.add("AnalysisQA/Reso/DCAxy_negdaughter_selected", "dcaXY of all processed tracks;d_{XY} (cm);Entries", HistType::kTH1F, {{1000, 0, 1}}); // check if cm is correct here
    resoRegistry.add("AnalysisQA/Reso/DCAz_negdaughter_selected", "dcaZ of all processed tracks;d_{Z} (cm);Entries", HistType::kTH1F, {{1000, 0, 1}});    // check if cm is correct here

    if (confEnableTriggerSelection) {
      for (const auto& triggerName : software_triggers::triggerNames) {
        if (confTriggerSwitches->get("Switch", triggerName.c_str())) {
          zorroTriggerNames += triggerName + ",";
        }
      }
      zorroTriggerNames.pop_back();
    }

    rctChecker.init(rctCut.cfgEvtRCTFlagCheckerLabel, false, rctCut.cfgEvtRCTFlagCheckerLimitAcceptAsBad);

    colCuts.setCuts(confEvtZvtx.value, confEvtTriggerCheck.value, confEvtTriggerSel.value, confEvtOfflineCheck.value, confEvtAddOfflineCheck.value, confIsRun3.value, confEvtMinSphericity.value, confEvtSphericityPtmin.value);
    colCuts.init(&qaRegistry);

    trackCuts.setSelection(Track.confTrkCharge, femtoDreamTrackSelection::kSign, femtoDreamSelection::kEqual);
    trackCuts.setSelection(Track.confTrkPtmin, femtoDreamTrackSelection::kpTMin, femtoDreamSelection::kLowerLimit);
    trackCuts.setSelection(Track.confTrkPtmax, femtoDreamTrackSelection::kpTMax, femtoDreamSelection::kUpperLimit);
    trackCuts.setSelection(Track.confTrkEta, femtoDreamTrackSelection::kEtaMax, femtoDreamSelection::kAbsUpperLimit);
    trackCuts.setSelection(Track.confTrkTPCnclsMin, femtoDreamTrackSelection::kTPCnClsMin, femtoDreamSelection::kLowerLimit);
    trackCuts.setSelection(Track.confTrkTPCfCls, femtoDreamTrackSelection::kTPCfClsMin, femtoDreamSelection::kLowerLimit);
    trackCuts.setSelection(Track.confTrkTPCcRowsMin, femtoDreamTrackSelection::kTPCcRowsMin, femtoDreamSelection::kLowerLimit);
    trackCuts.setSelection(Track.confTrkTPCsCls, femtoDreamTrackSelection::kTPCsClsMax, femtoDreamSelection::kUpperLimit);
    trackCuts.setSelection(Track.confTrkITSnclsMin, femtoDreamTrackSelection::kITSnClsMin, femtoDreamSelection::kLowerLimit);
    trackCuts.setSelection(Track.confTrkITSnclsIbMin, femtoDreamTrackSelection::kITSnClsIbMin, femtoDreamSelection::kLowerLimit);
    trackCuts.setSelection(Track.confTrkDCAxyMax, femtoDreamTrackSelection::kDCAxyMax, femtoDreamSelection::kAbsUpperLimit);
    trackCuts.setSelection(Track.confTrkDCAzMax, femtoDreamTrackSelection::kDCAzMax, femtoDreamSelection::kAbsUpperLimit);
    trackCuts.setSelection(Track.confTrkPIDnSigmaMax, femtoDreamTrackSelection::kPIDnSigmaMax, femtoDreamSelection::kAbsUpperLimit);
    trackCuts.setPIDSpecies(Track.confTrkPIDspecies);
    trackCuts.setnSigmaPIDOffset(Track.confTrkPIDnSigmaOffsetTPC, Track.confTrkPIDnSigmaOffsetTOF);
    trackCuts.init<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::TrackType::kNoChild, aod::femtodreamparticle::cutContainerType>(&qaRegistry, &trackRegistry);

    /// \todo fix how to pass array to setSelection, getRow() passing a
    /// different type!
    // v0Cuts.setSelection(confV0Selection->getRow(0),
    // femto_dream_v0_selection::kDecVtxMax, femtoDreamSelection::kAbsUpperLimit);
    if (confIsActivateV0) {
      v0Cuts.setSelection(confV0Sign, femto_dream_v0_selection::kV0Sign, femtoDreamSelection::kEqual);
      v0Cuts.setSelection(confV0PtMin, femto_dream_v0_selection::kV0pTMin, femtoDreamSelection::kLowerLimit);
      v0Cuts.setSelection(confV0PtMax, femto_dream_v0_selection::kV0pTMax, femtoDreamSelection::kUpperLimit);
      v0Cuts.setSelection(confV0EtaMax, femto_dream_v0_selection::kV0etaMax, femtoDreamSelection::kAbsUpperLimit);
      v0Cuts.setSelection(confV0DCADaughMax, femto_dream_v0_selection::kV0DCADaughMax, femtoDreamSelection::kUpperLimit);
      v0Cuts.setSelection(confV0CPAMin, femto_dream_v0_selection::kV0CPAMin, femtoDreamSelection::kLowerLimit);
      v0Cuts.setSelection(confV0TranRadMin, femto_dream_v0_selection::kV0TranRadMin, femtoDreamSelection::kLowerLimit);
      v0Cuts.setSelection(confV0TranRadMax, femto_dream_v0_selection::kV0TranRadMax, femtoDreamSelection::kUpperLimit);
      v0Cuts.setSelection(confV0DecVtxMax, femto_dream_v0_selection::kV0DecVtxMax, femtoDreamSelection::kUpperLimit);
      v0Cuts.setChildCuts(femto_dream_v0_selection::kPosTrack, confChildSign, femtoDreamTrackSelection::kSign, femtoDreamSelection::kEqual);
      v0Cuts.setChildCuts(femto_dream_v0_selection::kPosTrack, confChildEtaMax, femtoDreamTrackSelection::kEtaMax, femtoDreamSelection::kAbsUpperLimit);
      v0Cuts.setChildCuts(femto_dream_v0_selection::kPosTrack, confChildTPCnClsMin, femtoDreamTrackSelection::kTPCnClsMin, femtoDreamSelection::kLowerLimit);
      v0Cuts.setChildCuts(femto_dream_v0_selection::kPosTrack, confChildDCAMin, femtoDreamTrackSelection::kDCAMin, femtoDreamSelection::kAbsLowerLimit);
      v0Cuts.setChildCuts(femto_dream_v0_selection::kPosTrack, confChildPIDnSigmaMax, femtoDreamTrackSelection::kPIDnSigmaMax, femtoDreamSelection::kAbsUpperLimit);

      v0Cuts.setChildCuts(femto_dream_v0_selection::kNegTrack, confChildSign, femtoDreamTrackSelection::kSign, femtoDreamSelection::kEqual);
      v0Cuts.setChildCuts(femto_dream_v0_selection::kNegTrack, confChildEtaMax, femtoDreamTrackSelection::kEtaMax, femtoDreamSelection::kAbsUpperLimit);
      v0Cuts.setChildCuts(femto_dream_v0_selection::kNegTrack, confChildTPCnClsMin, femtoDreamTrackSelection::kTPCnClsMin, femtoDreamSelection::kLowerLimit);
      v0Cuts.setChildCuts(femto_dream_v0_selection::kNegTrack, confChildDCAMin, femtoDreamTrackSelection::kDCAMin, femtoDreamSelection::kAbsLowerLimit);
      v0Cuts.setChildCuts(femto_dream_v0_selection::kNegTrack, confChildPIDnSigmaMax, femtoDreamTrackSelection::kPIDnSigmaMax, femtoDreamSelection::kAbsUpperLimit);
      v0Cuts.setChildPIDSpecies(femto_dream_v0_selection::kPosTrack, confChildPIDspecies);
      v0Cuts.setChildPIDSpecies(femto_dream_v0_selection::kNegTrack, confChildPIDspecies);
      v0Cuts.init<aod::femtodreamparticle::ParticleType::kV0, aod::femtodreamparticle::ParticleType::kV0Child, aod::femtodreamparticle::cutContainerType>(&qaRegistry, &v0Registry);
      v0Cuts.setInvMassLimits(confV0InvMassLowLimit, confV0InvMassUpLimit);
      v0Cuts.setIsMother(confV0MotherIsLambda);

      v0Cuts.setChildRejectNotPropagatedTracks(femto_dream_v0_selection::kPosTrack, confTrkRejectNotPropagated);
      v0Cuts.setChildRejectNotPropagatedTracks(femto_dream_v0_selection::kNegTrack, confTrkRejectNotPropagated);

      v0Cuts.setnSigmaPIDOffsetTPC(Track.confTrkPIDnSigmaOffsetTPC);
      v0Cuts.setChildnSigmaPIDOffset(femto_dream_v0_selection::kPosTrack, Track.confTrkPIDnSigmaOffsetTPC, Track.confTrkPIDnSigmaOffsetTOF);
      v0Cuts.setChildnSigmaPIDOffset(femto_dream_v0_selection::kNegTrack, Track.confTrkPIDnSigmaOffsetTPC, Track.confTrkPIDnSigmaOffsetTOF);

      if (confV0RejectKaons) {
        v0Cuts.setKaonInvMassLimits(confV0InvKaonMassLowLimit, confV0InvKaonMassUpLimit);
      }
      v0Cuts.setRejectLambda(confV0RejectLambdas);
    }

    if (confIsActivateCascade) {
      // Cascades
      cascadeCuts.setSelection(confCascSel.confCascadeSign, femtoDreamCascadeSelection::kCascadeSign, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadeSign));
      cascadeCuts.setSelection(confCascSel.confCascadePtMin, femtoDreamCascadeSelection::kCascadePtMin, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadePtMin));
      cascadeCuts.setSelection(confCascSel.confCascadePtMax, femtoDreamCascadeSelection::kCascadePtMax, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadePtMax));
      cascadeCuts.setSelection(confCascSel.confCascadeEtaMax, femtoDreamCascadeSelection::kCascadeEtaMax, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadeEtaMax));
      cascadeCuts.setSelection(confCascSel.confCascadeDCADaughMax, femtoDreamCascadeSelection::kCascadeDCADaughMax, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadeDCADaughMax));
      cascadeCuts.setSelection(confCascSel.confCascadeCPAMin, femtoDreamCascadeSelection::kCascadeCPAMin, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadeCPAMin));
      cascadeCuts.setSelection(confCascSel.confCascadeTranRadMin, femtoDreamCascadeSelection::kCascadeTranRadMin, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadeTranRadMin));
      cascadeCuts.setSelection(confCascSel.confCascadeTranRadMax, femtoDreamCascadeSelection::kCascadeTranRadMax, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadeTranRadMax));
      cascadeCuts.setSelection(confCascSel.confCascadeDecVtxMax, femtoDreamCascadeSelection::kCascadeDecVtxMax, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadeDecVtxMax));
      // Cascade v0
      cascadeCuts.setSelection(confCascSel.confCascadeV0DCADaughMax, femtoDreamCascadeSelection::kCascadeV0DCADaughMax, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadeV0DCADaughMax));
      cascadeCuts.setSelection(confCascSel.confCascadeV0CPAMin, femtoDreamCascadeSelection::kCascadeV0CPAMin, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadeV0CPAMin));
      cascadeCuts.setSelection(confCascSel.confCascadeV0TranRadMin, femtoDreamCascadeSelection::kCascadeV0TranRadMin, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadeV0TranRadMin));
      cascadeCuts.setSelection(confCascSel.confCascadeV0TranRadMax, femtoDreamCascadeSelection::kCascadeV0TranRadMax, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadeV0TranRadMax));
      cascadeCuts.setSelection(confCascSel.confCascadeV0DCAtoPVMin, femtoDreamCascadeSelection::kCascadeV0DCAtoPVMin, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadeV0DCAtoPVMin));
      cascadeCuts.setSelection(confCascSel.confCascadeV0DCAtoPVMax, femtoDreamCascadeSelection::kCascadeV0DCAtoPVMax, FemtoDreamCascadeSelection::getSelectionType(femtoDreamCascadeSelection::kCascadeV0DCAtoPVMax));

      // Cascade Daughter Tracks
      cascadeCuts.setChildCuts(femtoDreamCascadeSelection::kPosTrack, confCascSel.confCascV0ChildSign, femtoDreamTrackSelection::kSign, femtoDreamSelection::kEqual);
      cascadeCuts.setChildCuts(femtoDreamCascadeSelection::kPosTrack, confCascSel.confCascV0ChildPtMin, femtoDreamTrackSelection::kpTMin, femtoDreamSelection::kLowerLimit);
      cascadeCuts.setChildCuts(femtoDreamCascadeSelection::kPosTrack, confCascSel.confCascV0ChildEtaMax, femtoDreamTrackSelection::kEtaMax, femtoDreamSelection::kAbsUpperLimit);
      cascadeCuts.setChildCuts(femtoDreamCascadeSelection::kPosTrack, confCascSel.confCascV0ChildTPCnClsMin, femtoDreamTrackSelection::kTPCnClsMin, femtoDreamSelection::kLowerLimit);
      cascadeCuts.setChildCuts(femtoDreamCascadeSelection::kPosTrack, confCascSel.confCascV0ChildDCAMin, femtoDreamTrackSelection::kDCAMin, femtoDreamSelection::kAbsLowerLimit);
      cascadeCuts.setChildCuts(femtoDreamCascadeSelection::kPosTrack, confCascSel.confCascV0ChildPIDnSigmaMax, femtoDreamTrackSelection::kPIDnSigmaMax, femtoDreamSelection::kAbsUpperLimit);
      cascadeCuts.setChildPIDSpecies(femtoDreamCascadeSelection::kPosTrack, confCascSel.confCascV0ChildPIDspecies);

      cascadeCuts.setChildCuts(femtoDreamCascadeSelection::kNegTrack, confCascSel.confCascV0ChildSign, femtoDreamTrackSelection::kSign, femtoDreamSelection::kEqual);
      cascadeCuts.setChildCuts(femtoDreamCascadeSelection::kNegTrack, confCascSel.confCascV0ChildPtMin, femtoDreamTrackSelection::kpTMin, femtoDreamSelection::kLowerLimit);
      cascadeCuts.setChildCuts(femtoDreamCascadeSelection::kNegTrack, confCascSel.confCascV0ChildEtaMax, femtoDreamTrackSelection::kEtaMax, femtoDreamSelection::kAbsUpperLimit);
      cascadeCuts.setChildCuts(femtoDreamCascadeSelection::kNegTrack, confCascSel.confCascV0ChildTPCnClsMin, femtoDreamTrackSelection::kTPCnClsMin, femtoDreamSelection::kLowerLimit);
      cascadeCuts.setChildCuts(femtoDreamCascadeSelection::kNegTrack, confCascSel.confCascV0ChildDCAMin, femtoDreamTrackSelection::kDCAMin, femtoDreamSelection::kAbsLowerLimit);
      cascadeCuts.setChildCuts(femtoDreamCascadeSelection::kNegTrack, confCascSel.confCascV0ChildPIDnSigmaMax, femtoDreamTrackSelection::kPIDnSigmaMax, femtoDreamSelection::kAbsUpperLimit);
      cascadeCuts.setChildPIDSpecies(femtoDreamCascadeSelection::kNegTrack, confCascSel.confCascV0ChildPIDspecies);

      // Cascade Bachelor Track
      cascadeCuts.setChildCuts(femtoDreamCascadeSelection::kBachTrack, confCascSel.confCascBachelorSign, femtoDreamTrackSelection::kSign, femtoDreamSelection::kEqual);
      cascadeCuts.setChildCuts(femtoDreamCascadeSelection::kBachTrack, confCascSel.confCascBachelorPtMin, femtoDreamTrackSelection::kpTMin, femtoDreamSelection::kLowerLimit);
      cascadeCuts.setChildCuts(femtoDreamCascadeSelection::kBachTrack, confCascSel.confCascBachelorEtaMax, femtoDreamTrackSelection::kEtaMax, femtoDreamSelection::kAbsUpperLimit);
      cascadeCuts.setChildCuts(femtoDreamCascadeSelection::kBachTrack, confCascSel.confCascBachelorTPCnClsMin, femtoDreamTrackSelection::kTPCnClsMin, femtoDreamSelection::kLowerLimit);
      cascadeCuts.setChildCuts(femtoDreamCascadeSelection::kBachTrack, confCascSel.confCascBachelorDCAMin, femtoDreamTrackSelection::kDCAMin, femtoDreamSelection::kAbsLowerLimit);
      cascadeCuts.setChildCuts(femtoDreamCascadeSelection::kBachTrack, confCascSel.confCascBachelorPIDnSigmaMax, femtoDreamTrackSelection::kPIDnSigmaMax, femtoDreamSelection::kAbsUpperLimit);
      cascadeCuts.setChildPIDSpecies(femtoDreamCascadeSelection::kBachTrack, confCascSel.confCascBachelorPIDspecies);

      cascadeCuts.init<aod::femtodreamparticle::ParticleType::kCascade, aod::femtodreamparticle::ParticleType::kCascadeV0Child, aod::femtodreamparticle::ParticleType::kCascadeBachelor, aod::femtodreamparticle::cutContainerType>(&qaRegistry, &cascadeRegistry, confCascSel.confCascIsSelectedOmega);
      cascadeCuts.setInvMassLimits(confCascSel.confCascInvMassLowLimit, confCascSel.confCascInvMassUpLimit);
      cascadeCuts.setV0InvMassLimits(confCascSel.confCascV0InvMassLowLimit, confCascSel.confCascV0InvMassUpLimit);
      if (confCascSel.confCascRejectCompetingMass) {
        cascadeCuts.setCompetingInvMassLimits(confCascSel.confCascInvCompetingMassLowLimit, confCascSel.confCascInvCompetingMassUpLimit);
      }
    }

    if (confIsActivateReso.value) {
      resoCuts.setDaughterCuts(femto_dream_reso_selection::kPosdaugh, Resonance.confDaughterCharge, femtoDreamTrackSelection::kSign, femtoDreamSelection::kEqual);
      resoCuts.setDaughterCuts(femto_dream_reso_selection::kPosdaugh, Resonance.confDaughterPtMin, femtoDreamTrackSelection::kpTMin, femtoDreamSelection::kLowerLimit);
      resoCuts.setDaughterCuts(femto_dream_reso_selection::kPosdaugh, Resonance.confDaughterPtMax, femtoDreamTrackSelection::kpTMax, femtoDreamSelection::kUpperLimit);
      resoCuts.setDaughterCuts(femto_dream_reso_selection::kPosdaugh, Resonance.confDaughterEtaMax, femtoDreamTrackSelection::kEtaMax, femtoDreamSelection::kAbsUpperLimit);
      resoCuts.setDaughterCuts(femto_dream_reso_selection::kPosdaugh, Resonance.confDaughterTPCnClsMin, femtoDreamTrackSelection::kTPCnClsMin, femtoDreamSelection::kLowerLimit);
      resoCuts.setDaughterCuts(femto_dream_reso_selection::kPosdaugh, Resonance.confDaughterTPCfClsMin, femtoDreamTrackSelection::kTPCfClsMin, femtoDreamSelection::kLowerLimit);
      resoCuts.setDaughterCuts(femto_dream_reso_selection::kPosdaugh, Resonance.confDaughterTPCcRowsMin, femtoDreamTrackSelection::kTPCcRowsMin, femtoDreamSelection::kLowerLimit);
      resoCuts.setDaughterCuts(femto_dream_reso_selection::kPosdaugh, Resonance.confDaughterDCAxyMax, femtoDreamTrackSelection::kDCAxyMax, femtoDreamSelection::kAbsUpperLimit);
      resoCuts.setDaughterCuts(femto_dream_reso_selection::kPosdaugh, Resonance.confDaughterDCAzMax, femtoDreamTrackSelection::kDCAzMax, femtoDreamSelection::kAbsUpperLimit);
      resoCuts.setDaughterCuts(femto_dream_reso_selection::kPosdaugh, Resonance.confDaughterPIDnSigmaMax, femtoDreamTrackSelection::kPIDnSigmaMax, femtoDreamSelection::kAbsUpperLimit);

      resoCuts.setDaughterCuts(femto_dream_reso_selection::kNegdaugh, Resonance.confDaughterCharge, femtoDreamTrackSelection::kSign, femtoDreamSelection::kEqual);
      resoCuts.setDaughterCuts(femto_dream_reso_selection::kNegdaugh, Resonance.confDaughterPtMin, femtoDreamTrackSelection::kpTMin, femtoDreamSelection::kLowerLimit);
      resoCuts.setDaughterCuts(femto_dream_reso_selection::kNegdaugh, Resonance.confDaughterPtMax, femtoDreamTrackSelection::kpTMax, femtoDreamSelection::kUpperLimit);
      resoCuts.setDaughterCuts(femto_dream_reso_selection::kNegdaugh, Resonance.confDaughterEtaMax, femtoDreamTrackSelection::kEtaMax, femtoDreamSelection::kAbsUpperLimit);
      resoCuts.setDaughterCuts(femto_dream_reso_selection::kNegdaugh, Resonance.confDaughterTPCnClsMin, femtoDreamTrackSelection::kTPCnClsMin, femtoDreamSelection::kLowerLimit);
      resoCuts.setDaughterCuts(femto_dream_reso_selection::kNegdaugh, Resonance.confDaughterTPCfClsMin, femtoDreamTrackSelection::kTPCfClsMin, femtoDreamSelection::kLowerLimit);
      resoCuts.setDaughterCuts(femto_dream_reso_selection::kNegdaugh, Resonance.confDaughterTPCcRowsMin, femtoDreamTrackSelection::kTPCcRowsMin, femtoDreamSelection::kLowerLimit);
      resoCuts.setDaughterCuts(femto_dream_reso_selection::kNegdaugh, Resonance.confDaughterDCAxyMax, femtoDreamTrackSelection::kDCAxyMax, femtoDreamSelection::kAbsUpperLimit);
      resoCuts.setDaughterCuts(femto_dream_reso_selection::kNegdaugh, Resonance.confDaughterDCAzMax, femtoDreamTrackSelection::kDCAzMax, femtoDreamSelection::kAbsUpperLimit);
      resoCuts.setDaughterCuts(femto_dream_reso_selection::kNegdaugh, Resonance.confDaughterPIDnSigmaMax, femtoDreamTrackSelection::kPIDnSigmaMax, femtoDreamSelection::kAbsUpperLimit);

      resoCuts.init<aod::femtodreamparticle::ParticleType::kReso,
                    aod::femtodreamparticle::ParticleType::kResoChild>(&qaRegistry, &v0Registry);

      resoCuts.assign(Resonance.confDaughterPTPCThr); // assigns Configurable value to class member
      resoCuts.setDaughterPIDSpecies(femto_dream_reso_selection::kPosdaugh, Resonance.confDaughterPIDspecies);
      resoCuts.setDaughterPIDSpecies(femto_dream_reso_selection::kNegdaugh, Resonance.confDaughterPIDspecies);
      resoCuts.setDaughternSigmaPIDOffset(femto_dream_reso_selection::kPosdaugh, 0.f, 0.f);
      resoCuts.setDaughternSigmaPIDOffset(femto_dream_reso_selection::kNegdaugh, 0.f, 0.f);

      resoCuts.setSelection(Resonance.confResoSign, femto_dream_reso_selection::kResoSign, femtoDreamSelection::kEqual);

      // resoCuts.init<>();
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
  void initCcdbMagTrig(aod::BCsWithTimestamps::iterator bc)
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

    // Init for zorro to get trigger flags
    if (confEnableTriggerSelection) {
      zorro.setCCDBpath(confBaseCCDBPathForTriggers);
      zorro.initCCDB(ccdb.service, mRunNumber, timestamp, zorroTriggerNames);
    }
  }

  template <bool isTrackOrV0, bool hasItsPid, typename ParticleType, bool isReso = false>
  void fillDebugParticle(ParticleType const& particle)
  {
    if constexpr (isTrackOrV0) {
      if constexpr (hasItsPid) {
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
                         particle.tpcNSigmaEl(),
                         particle.tpcNSigmaPi(),
                         particle.tpcNSigmaKa(),
                         particle.tpcNSigmaPr(),
                         particle.tpcNSigmaDe(),
                         particle.tpcNSigmaTr(),
                         particle.tpcNSigmaHe(),
                         particle.tofNSigmaEl(),
                         particle.tofNSigmaPi(),
                         particle.tofNSigmaKa(),
                         particle.tofNSigmaPr(),
                         particle.tofNSigmaDe(),
                         particle.tofNSigmaTr(),
                         particle.tofNSigmaHe(),
                         o2::analysis::femtoDream::itsSignal(particle),
                         particle.itsNSigmaEl(),
                         particle.itsNSigmaPi(),
                         particle.itsNSigmaKa(),
                         particle.itsNSigmaPr(),
                         particle.itsNSigmaDe(),
                         particle.itsNSigmaTr(),
                         particle.itsNSigmaHe(),
                         -999., -999., -999., -999., -999., -999.,
                         -999., -999., -999., -999., -999., -999., -999.);
      } else {
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
                         particle.tpcNSigmaEl(),
                         particle.tpcNSigmaPi(),
                         particle.tpcNSigmaKa(),
                         particle.tpcNSigmaPr(),
                         particle.tpcNSigmaDe(),
                         particle.tpcNSigmaTr(),
                         particle.tpcNSigmaHe(),
                         particle.tofNSigmaEl(),
                         particle.tofNSigmaPi(),
                         particle.tofNSigmaKa(),
                         particle.tofNSigmaPr(),
                         particle.tofNSigmaDe(),
                         particle.tofNSigmaTr(),
                         particle.tofNSigmaHe(),
                         -999., -999., -999., -999., -999., -999., -999., -999.,
                         -999., -999., -999., -999., -999., -999.,
                         -999., -999., -999., -999., -999., -999., -999.);
      }
    } else if constexpr (isReso) {
      outputDebugParts(-999., -999., -999., -999., -999., -999., -999., -999., -999., -999., // for the moment
                       -999., -999., -999., -999., -999., -999., -999., -999., -999., -999.,
                       -999., -999., -999., -999., -999., -999., -999., -999., -999., -999.,
                       -999., -999., -999., -999., -999., -999., -999., -999., -999., -999.,
                       -999., -999., -999., -999., -999., -999.);
    } else {
      outputDebugParts(-999.,                                                         // sign
                       -999., -999., -999., -999., -999., -999., -999., -999., -999., // track properties (DCA, NCls, crossed rows, etc.)
                       -999., -999., -999., -999., -999., -999., -999., -999.,        // TPC PID (TPC signal + particle hypothesis)
                       -999., -999., -999., -999., -999., -999., -999.,               // TOF PID
                       -999., -999., -999., -999., -999., -999., -999., -999.,        // ITS PID
                       particle.dcaV0daughters(),
                       particle.v0radius(),
                       particle.x(),
                       particle.y(),
                       particle.z(),
                       particle.mK0Short(),
                       -999., -999., -999., -999., -999., -999., -999.); // Cascade properties
    }
  }

  template <typename ParticleType, typename CollisionType>
  void fillDebugCascade(ParticleType const& cascade, CollisionType const& col)
  {
    outputDebugParts(cascade.sign(),                                                // sign
                     -999., -999., -999., -999., -999., -999., -999., -999., -999., // track properties (DCA, NCls, crossed rows, etc.)
                     -999., -999., -999., -999., -999., -999., -999., -999.,        // TPC PID (TPC signal + particle hypothesis)
                     -999., -999., -999., -999., -999., -999., -999.,               // TOF PID
                     -999., -999., -999., -999., -999., -999., -999., -999.,        // ITS PID
                     cascade.dcaV0daughters(),
                     cascade.v0radius(),
                     -999., // DecVtxV0 x
                     -999., // DecVtxV0 y
                     -999., // DecVtxV0 z
                     -999., // mKaon
                     cascade.dcav0topv(col.posX(), col.posY(), col.posZ()),
                     cascade.dcacascdaughters(),
                     cascade.cascradius(),
                     cascade.x(),
                     cascade.y(),
                     cascade.z(),
                     cascade.mOmega()); // QA for Reso
  }

  template <typename CollisionType, typename ParticleType>
  void fillMCParticle(CollisionType const& col, ParticleType const& particle, o2::aod::femtodreamparticle::ParticleType fdparttype)
  {
    if (particle.has_mcParticle()) {

      constexpr int ProcessDirectMother = 4;
      constexpr int ProcessInelasticHadronic = 23;
      constexpr int GenStatusTransport = -1;

      // get corresponding MC particle and its info
      auto particleMC = particle.mcParticle();
      auto pdgCode = particleMC.pdgCode();
      trackRegistry.fill(HIST("AnalysisQA/Particle"), pdgCode);
      int particleOrigin = 99;
      int pdgCodeMother = -1;
      // get list of mothers, but it could be empty (for example in case of injected light nuclei)
      auto motherparticlesMC = particleMC.template mothers_as<aod::McParticles>();
      // check pdg code
      trackRegistry.fill(HIST("AnalysisQA/getGenStatusCode"), particleMC.getGenStatusCode());
      trackRegistry.fill(HIST("AnalysisQA/getProcess"), particleMC.getProcess());
      // if this fails, the particle is a fake
      if (std::abs(pdgCode) == std::abs(confTrkPDGCode.value)) {
        // check first if particle is from pile up
        // check if the collision associated with the particle is the same as the analyzed collision by checking their Ids
        if ((col.has_mcCollision() && (particleMC.mcCollisionId() != col.mcCollisionId())) || !col.has_mcCollision()) {
          particleOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kWrongCollision;
          // check if particle is primary
        } else if (particleMC.isPhysicalPrimary()) {
          particleOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kPrimary;
          // check if particle is secondary
          // particle is from a decay -> getProcess() == 4
          // particle is generated during transport -> getGenStatusCode() == -1
          // list of mothers is not empty
        } else if (particleMC.getProcess() == ProcessDirectMother && particleMC.getGenStatusCode() == GenStatusTransport && !motherparticlesMC.empty()) {
          // get direct mother
          auto motherparticleMC = motherparticlesMC.front();
          pdgCodeMother = motherparticleMC.pdgCode();
          trackRegistry.fill(HIST("AnalysisQA/Mother"), pdgCodeMother);
          particleOrigin = checkDaughterType(fdparttype, motherparticleMC.pdgCode());
          // check if particle is material
          // particle is from inelastic hadronic interaction -> getProcess() == 23
          // particle is generated during transport -> getGenStatusCode() == -1
        } else if (particleMC.getProcess() == ProcessInelasticHadronic && particleMC.getGenStatusCode() == GenStatusTransport) {
          particleOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kMaterial;
          // cross check to see if we missed a case
        } else {
          particleOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kElse;
        }
        // if pdg code is wrong, particle is fake
      } else {
        particleOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kFake;
      }

      outputPartsMC(particleOrigin, pdgCode, particleMC.pt(), particleMC.eta(), particleMC.phi());
      outputPartsMCLabels(outputPartsMC.lastIndex());
      if (confIsDebug) {
        outputPartsExtMCLabels(outputPartsMC.lastIndex());
        outputDebugPartsMC(pdgCodeMother);
      }
    } else {
      outputPartsMCLabels(-1);
      if (confIsDebug) {
        outputPartsExtMCLabels(-1);
      }
    }
  }

  template <typename CollisionType>
  void fillMCCollision(CollisionType const& col)
  {
    if (col.has_mcCollision()) {
      auto genMCcol = col.template mcCollision_as<aod::FemtoFullMCgenCollisions>();
      outputMCCollision(genMCcol.multMCNParticlesEta08());
      outputCollsMCLabels(outputMCCollision.lastIndex());
    } else {
      outputCollsMCLabels(-1);
    }
  }

  template <typename P>
  void fillLikeSign(P const& sliceDaughters)
  {
    for (const auto& track1 : sliceDaughters) {
      if (!resoCuts.daughterSelectionPos(track1) || !resoCuts.isSelectedMinimalPIDPos(track1, Resonance.confDaughterPIDspecies.value)) {
        continue;
      }
      for (const auto& track2 : sliceDaughters) {
        if (!resoCuts.daughterSelectionPos(track2) || !resoCuts.isSelectedMinimalPIDPos(track2, Resonance.confDaughterPIDspecies.value)) {
          continue;
        }

        /// This only works for the case where the mass of opposite charged particles are the same (for example K+/K- have same mass)
        float massPart1 = o2::track::PID::getMass(Resonance.confDaughterPIDspecies.value[0]);
        float massPart2 = o2::track::PID::getMass(Resonance.confDaughterPIDspecies.value[1]);

        /// Resonance
        ROOT::Math::PtEtaPhiMVector tempD1(track1.pt(), track1.eta(), track1.phi(), massPart1);
        ROOT::Math::PtEtaPhiMVector tempD2(track2.pt(), track2.eta(), track2.phi(), massPart2);
        ROOT::Math::PtEtaPhiMVector tempReso = tempD1 + tempD2;
        /// Anti-resonance
        ROOT::Math::PtEtaPhiMVector tempDA1(track1.pt(), track1.eta(), track1.phi(), massPart2);
        ROOT::Math::PtEtaPhiMVector tempDA2(track2.pt(), track2.eta(), track2.phi(), massPart1);
        ROOT::Math::PtEtaPhiMVector tempAntiReso = tempDA1 + tempDA2;

        /// For InvMassQA
        ROOT::Math::PtEtaPhiMVector tempDaughter1MassQA1(track1.pt(), track1.eta(), track1.phi(), massPart1);
        ROOT::Math::PtEtaPhiMVector tempDaughter2MassQA1(track2.pt(), track2.eta(), track2.phi(), massPart1);
        ROOT::Math::PtEtaPhiMVector tempMassQA1 = tempDaughter1MassQA1 + tempDaughter2MassQA1;

        float massQAPart2 = massPart2;
        if (Resonance.confDaughterPIDspecies.value[0] == Resonance.confDaughterPIDspecies.value[1]) {
          massQAPart2 = o2::track::PID::getMass(Resonance.confMassQAPart2PID.value);
        }

        ROOT::Math::PtEtaPhiMVector tempDaughter1MassQA2(track1.pt(), track1.eta(), track1.phi(), massQAPart2);
        ROOT::Math::PtEtaPhiMVector tempDaughter2MassQA2(track2.pt(), track2.eta(), track2.phi(), massQAPart2);
        ROOT::Math::PtEtaPhiMVector tempMassQA2 = tempDaughter1MassQA2 + tempDaughter2MassQA2;

        float massDiff = std::abs(Resonance.confResoInvMass.value - tempReso.M());
        float massDiffAnti = std::abs(Resonance.confResoInvMass.value - tempAntiReso.M());

        bool resoIsNotAnti = true; /// bool for differentianting between particle/antiparticle
        if ((Resonance.confDaughterPIDspecies->size() > 1) && (Resonance.confDaughterPIDspecies.value[0] != Resonance.confDaughterPIDspecies.value[1])) {
          auto [isNormal, WrongCombination] = resoCuts.checkCombination(track1, track2, Resonance.confDaughterPIDspecies.value, massDiff, massDiffAnti, Resonance.confUseMassDiffLikeSign.value);
          if (WrongCombination) {
            continue;
          }
          resoIsNotAnti = isNormal;
        }
        /// Resos, where both daughters have the same PID are defaulted to sign 1. and resoIsNotAnti = true

        if (resoIsNotAnti) {
          resoRegistry.fill(HIST("AnalysisQA/ResoLikeSign/ResoQA/InvMass"), tempReso.M());
          resoRegistry.fill(HIST("AnalysisQA/ResoLikeSign/ResoQA/InvMassSwitched"), tempAntiReso.M());
          resoRegistry.fill(HIST("AnalysisQA/ResoLikeSign/ResoQA/InvMassBothPID1"), tempMassQA1.M());
          resoRegistry.fill(HIST("AnalysisQA/ResoLikeSign/ResoQA/InvMassBothPID2"), tempMassQA2.M());
        } else {
          resoRegistry.fill(HIST("AnalysisQA/ResoLikeSign/AntiResoQA/InvMass"), tempReso.M());
          resoRegistry.fill(HIST("AnalysisQA/ResoLikeSign/AntiResoQA/InvMassSwitched"), tempAntiReso.M());
          resoRegistry.fill(HIST("AnalysisQA/ResoLikeSign/AntiResoQA/InvMassBothPID1"), tempMassQA1.M());
          resoRegistry.fill(HIST("AnalysisQA/ResoLikeSign/AntiResoQA/InvMassBothPID2"), tempMassQA2.M());
        }
      } // for (const &auto track2 : sliceDaughters)
    } // for (const &auto track1 : sliceDaughters)
  }

  template <bool isMC, bool hasItsPid, bool useCentrality, bool analysePbPb, typename CascadeType, typename V0Type, typename TrackType, typename TrackTypeWithItsPid, typename CollisionType>
  void fillCollisionsAndTracksAndV0AndCascade(CollisionType const& col, TrackType const& tracks, TrackTypeWithItsPid const& tracksWithItsPid, V0Type const& fullV0s, CascadeType const& fullCascades)
  {
    // If triggering is enabled, select only events which were triggered wit our triggers
    if (confEnableTriggerSelection) {
      bool zorroSelected = zorro.isSelected(col.template bc_as<aod::BCsWithTimestamps>().globalBC()); /// check if event was selected by triggers of interest
      if (!zorroSelected) {
        return;
      }
    }

    const auto vtxZ = col.posZ();
    const auto spher = colCuts.computeSphericity(col, tracks);
    float mult = 0;
    int multNtr = 0;
    if (confIsRun3) {
      if constexpr (useCentrality) {
        if constexpr (analysePbPb) {
          mult = col.centFT0C();
        } else {
          mult = col.centFT0M();
        }
      } else {
        mult = 0;
      }
      multNtr = col.multNTracksPV();
    } else {
      mult = 1; // multiplicity percentile is know in Run 2
      multNtr = col.multTracklets();
    }

    colCuts.fillQA(col, mult);

    // check whether the basic event selection criteria are fulfilled
    // that included checking if there is at least on usable track or V0
    if (!colCuts.isSelectedCollision(col)) {
      return;
    }

    if (confIsActivateCascade.value) {
      if (colCuts.isEmptyCollision(col, tracks, trackCuts) && colCuts.isCollisionWithoutTrkCasc(col, fullCascades, cascadeCuts, tracks)) {
        return;
      }
    }

    if (confIsActivateV0.value) {
      if (colCuts.isEmptyCollision(col, tracks, trackCuts) && colCuts.isEmptyCollision(col, fullV0s, v0Cuts, tracks)) {
        return;
      }
    } else {
      if (colCuts.isEmptyCollision(col, tracks, trackCuts)) {
        return;
      }
    }

    if (rctCut.requireRCTFlagChecker && !rctChecker(col)) {
      return;
    }

    outputCollision(vtxZ, mult, multNtr, spher, mMagField);
    if constexpr (isMC) {
      fillMCCollision(col);
    }

    std::vector<int> childIDs = {0, 0};           // these IDs are necessary to keep track of the children
    std::vector<int> cascadechildIDs = {0, 0, 0}; // these IDs are necessary to keep track of the children
    std::vector<int> tmpIDtrack;                  // this vector keeps track of the matching of the primary track table row <-> aod::track table global index

    for (const auto& track : tracksWithItsPid) {

      /// if the most open selection criteria are not fulfilled there is no
      /// point looking further at the track
      trackCuts.fillQA<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::TrackType::kNoChild, 0>(track);

      if (track.tpcChi2NCl() < OptionTrackSpecialSelections.confTrkMinChi2PerClusterTPC || track.tpcChi2NCl() > OptionTrackSpecialSelections.confTrkMaxChi2PerClusterTPC) {
        continue;
      }
      if (track.itsChi2NCl() > OptionTrackSpecialSelections.confTrkMaxChi2PerClusterITS) {
        continue;
      }
      if ((OptionTrackSpecialSelections.confTrkTPCRefit && !track.hasTPC()) || (OptionTrackSpecialSelections.confTrkITSRefit && !track.hasITS())) {
        continue;
      }

      if (!trackCuts.isSelectedMinimal(track)) {
        continue;
      }

      trackRegistry.fill(HIST("AnalysisQA/Chi2ITSTPCperCluster"), track.itsChi2NCl(), track.tpcChi2NCl());
      trackRegistry.fill(HIST("AnalysisQA/RefitITSTPC"), track.hasITS(), track.hasTPC());

      trackCuts.fillQA<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::TrackType::kNoChild, 1>(track);
      // the bit-wise container of the systematic variations is obtained
      std::array<o2::aod::femtodreamparticle::cutContainerType, 2> cutContainer;
      cutContainer = trackCuts.getCutContainer<hasItsPid, aod::femtodreamparticle::cutContainerType>(track, track.pt(), track.eta(), sqrtf(powf(track.dcaXY(), 2.f) + powf(track.dcaZ(), 2.f)));

      // now the table is filled
      outputParts(outputCollision.lastIndex(),
                  track.pt(),
                  track.eta(),
                  track.phi(),
                  aod::femtodreamparticle::ParticleType::kTrack,
                  cutContainer.at(femtoDreamTrackSelection::TrackContainerPosition::kCuts),
                  cutContainer.at(femtoDreamTrackSelection::TrackContainerPosition::kPID),
                  track.dcaXY(), childIDs, 0, 0);
      tmpIDtrack.push_back(track.globalIndex());
      if (confIsDebug.value) {
        fillDebugParticle<true, hasItsPid>(track);
      }

      if constexpr (isMC) {
        fillMCParticle(col, track, o2::aod::femtodreamparticle::ParticleType::kTrack);
      }
    }

    if (confIsActivateV0.value) {
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

        // if (confRejectITSHitandTOFMissing) {
        // Uncomment only when TOF timing is solved
        // bool itsHit = o2PhysicsTrackSelection->IsSelected(postrack,
        // TrackSelection::TrackCuts::kITSHits); bool itsHit =
        // o2PhysicsTrackSelection->IsSelected(negtrack,
        // TrackSelection::TrackCuts::kITSHits);
        // }

        v0Cuts.fillQA<aod::femtodreamparticle::ParticleType::kV0, aod::femtodreamparticle::ParticleType::kV0Child>(col, v0, postrack, negtrack); ///\todo fill QA also for daughters
        auto cutContainerV0 = v0Cuts.getCutContainer<aod::femtodreamparticle::cutContainerType>(col, v0, postrack, negtrack);

        int postrackID = v0.posTrackId();
        int rowInPrimaryTrackTablePos = -1;
        rowInPrimaryTrackTablePos = getRowDaughters(postrackID, tmpIDtrack);
        childIDs[0] = rowInPrimaryTrackTablePos;
        childIDs[1] = 0;
        outputParts(outputCollision.lastIndex(),
                    v0.positivept(), v0.positiveeta(), v0.positivephi(),
                    aod::femtodreamparticle::ParticleType::kV0Child,
                    cutContainerV0.at(femto_dream_v0_selection::V0ContainerPosition::kPosCuts),
                    cutContainerV0.at(femto_dream_v0_selection::V0ContainerPosition::kPosPID),
                    postrack.dcaXY(),
                    childIDs,
                    0,
                    0);
        const int rowOfPosTrack = outputParts.lastIndex();
        if constexpr (isMC) {
          fillMCParticle(col, postrack, o2::aod::femtodreamparticle::ParticleType::kV0Child);
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
                    aod::femtodreamparticle::ParticleType::kV0Child,
                    cutContainerV0.at(femto_dream_v0_selection::V0ContainerPosition::kNegCuts),
                    cutContainerV0.at(femto_dream_v0_selection::V0ContainerPosition::kNegPID),
                    negtrack.dcaXY(),
                    childIDs,
                    0,
                    0);
        const int rowOfNegTrack = outputParts.lastIndex();
        if constexpr (isMC) {
          fillMCParticle(col, negtrack, o2::aod::femtodreamparticle::ParticleType::kV0Child);
        }
        std::vector<int> indexChildID = {rowOfPosTrack, rowOfNegTrack};
        float fillMass = 0.;
        float fillMassAnti = 0.;
        if (confV0MotherIsLambda) {
          fillMass = v0.mLambda();
          fillMassAnti = v0.mAntiLambda();
        } else {
          fillMass = v0.mK0Short();
          fillMassAnti = fillMass;
        }
        outputParts(outputCollision.lastIndex(),
                    v0.pt(),
                    v0.eta(),
                    v0.phi(),
                    aod::femtodreamparticle::ParticleType::kV0,
                    cutContainerV0.at(femto_dream_v0_selection::V0ContainerPosition::kV0),
                    0,
                    v0.v0cosPA(),
                    indexChildID,
                    fillMass,
                    fillMassAnti);
        if (confIsDebug.value) {
          fillDebugParticle<true, false>(postrack); // QA for positive daughter
          fillDebugParticle<true, false>(negtrack); // QA for negative daughter
          fillDebugParticle<false, false>(v0);      // QA for v0
        }
        if constexpr (isMC) {
          fillMCParticle(col, v0, o2::aod::femtodreamparticle::ParticleType::kV0);
        }
      }
    }

    if (confIsActivateCascade.value) {
      for (const auto& casc : fullCascades) {
        // get the daughter tracks
        const auto& posTrackCasc = casc.template posTrack_as<TrackType>();
        const auto& negTrackCasc = casc.template negTrack_as<TrackType>();
        const auto& bachTrackCasc = casc.template bachelor_as<TrackType>();

        cascadeCuts.fillQA<0, aod::femtodreamparticle::ParticleType::kCascade>(col, casc, posTrackCasc, negTrackCasc, bachTrackCasc);
        if (!cascadeCuts.isSelectedMinimal(col, casc, posTrackCasc, negTrackCasc, bachTrackCasc)) {
          continue;
        }
        cascadeCuts.fillQA<1, aod::femtodreamparticle::ParticleType::kCascade>(col, casc, posTrackCasc, negTrackCasc, bachTrackCasc);

        // auto cutContainerCasc = cascadeCuts.getCutContainer<aod::femtodreamparticle::cutContainerType>(col, casc, v0daugh, posTrackCasc, negTrackCasc, bachTrackCasc);
        auto cutContainerCasc = cascadeCuts.getCutContainer<aod::femtodreamparticle::cutContainerType>(col, casc, posTrackCasc, negTrackCasc, bachTrackCasc);

        // Fill positive child
        int poscasctrackID = casc.posTrackId();
        int rowInPrimaryTrackTablePosCasc = -1;
        rowInPrimaryTrackTablePosCasc = getRowDaughters(poscasctrackID, tmpIDtrack);
        cascadechildIDs[0] = rowInPrimaryTrackTablePosCasc;
        cascadechildIDs[1] = 0;
        cascadechildIDs[2] = 0;
        outputParts(outputCollision.lastIndex(),
                    posTrackCasc.pt(),
                    posTrackCasc.eta(),
                    posTrackCasc.phi(),
                    aod::femtodreamparticle::ParticleType::kCascadeV0Child,
                    cutContainerCasc.at(femtoDreamCascadeSelection::CascadeContainerPosition::kPosCuts),
                    cutContainerCasc.at(femtoDreamCascadeSelection::CascadeContainerPosition::kPosPID),
                    posTrackCasc.dcaXY(),
                    cascadechildIDs,
                    0,
                    0);
        const int rowOfPosCascadeTrack = outputParts.lastIndex();
        // TODO: include here MC filling
        //------

        // Fill negative child
        int negcasctrackID = casc.negTrackId();
        int rowInPrimaryTrackTableNegCasc = -1;
        rowInPrimaryTrackTableNegCasc = getRowDaughters(negcasctrackID, tmpIDtrack);
        cascadechildIDs[0] = 0;
        cascadechildIDs[1] = rowInPrimaryTrackTableNegCasc;
        cascadechildIDs[2] = 0;
        outputParts(outputCollision.lastIndex(),
                    negTrackCasc.pt(),
                    negTrackCasc.eta(),
                    negTrackCasc.phi(),
                    aod::femtodreamparticle::ParticleType::kCascadeV0Child,
                    cutContainerCasc.at(femtoDreamCascadeSelection::CascadeContainerPosition::kNegCuts),
                    cutContainerCasc.at(femtoDreamCascadeSelection::CascadeContainerPosition::kNegPID),
                    negTrackCasc.dcaXY(),
                    cascadechildIDs,
                    0,
                    0);
        const int rowOfNegCascadeTrack = outputParts.lastIndex();
        // TODO: include here MC filling
        //------

        // Fill bachelor child
        int bachelorcasctrackID = casc.bachelorId();
        int rowInPrimaryTrackTableBachelorCasc = -1;
        rowInPrimaryTrackTableBachelorCasc = getRowDaughters(bachelorcasctrackID, tmpIDtrack);
        cascadechildIDs[0] = 0;
        cascadechildIDs[1] = 0;
        cascadechildIDs[2] = rowInPrimaryTrackTableBachelorCasc;
        outputParts(outputCollision.lastIndex(),
                    bachTrackCasc.pt(),
                    bachTrackCasc.eta(),
                    bachTrackCasc.phi(),
                    aod::femtodreamparticle::ParticleType::kCascadeBachelor,
                    cutContainerCasc.at(femtoDreamCascadeSelection::CascadeContainerPosition::kBachCuts),
                    cutContainerCasc.at(femtoDreamCascadeSelection::CascadeContainerPosition::kBachPID),
                    bachTrackCasc.dcaXY(),
                    cascadechildIDs,
                    0,
                    0);
        const int rowOfBachelorCascadeTrack = outputParts.lastIndex();
        // TODO: include here MC filling
        //------

        // Fill cascades
        float invMassCasc = confCascSel.confCascIsSelectedOmega ? casc.mOmega() : casc.mXi();
        std::vector<int> indexCascadeChildID = {rowOfPosCascadeTrack, rowOfNegCascadeTrack, rowOfBachelorCascadeTrack};
        outputParts(outputCollision.lastIndex(),
                    casc.pt(),
                    casc.eta(),
                    casc.phi(),
                    aod::femtodreamparticle::ParticleType::kCascade,
                    cutContainerCasc.at(femtoDreamCascadeSelection::CascadeContainerPosition::kCascade),
                    0,
                    casc.casccosPA(col.posX(), col.posY(), col.posZ()),
                    indexCascadeChildID,
                    invMassCasc,
                    casc.mLambda());
        // TODO: include here MC filling
        //------

        if (confIsDebug.value) {
          fillDebugParticle<true, false>(posTrackCasc);  // QA for positive daughter
          fillDebugParticle<true, false>(negTrackCasc);  // QA for negative daughter
          fillDebugParticle<true, false>(bachTrackCasc); // QA for negative daughter
          fillDebugCascade(casc, col);                   // QA for Cascade
        }
      }
    }

    if (confIsActivatePhi.value) {

      resoCuts.updateThreshold();

      auto slicePosdaugh = daughter1.sliceByCached(aod::track::collisionId, col.globalIndex(), cache); // o2::framework defined in AnalysisHelper.h
      auto sliceNegdaugh = daughter2.sliceByCached(aod::track::collisionId, col.globalIndex(), cache);

      for (const auto& track1 : slicePosdaugh) {
        if (!resoCuts.daughterSelectionPos(track1) || !resoCuts.isSelectedMinimalPIDPos(track1, Resonance.confDaughterPIDspecies.value)) {
          continue;
        }

        for (const auto& track2 : sliceNegdaugh) {
          if (!resoCuts.daughterSelectionNeg(track2) || !resoCuts.isSelectedMinimalPIDNeg(track2, Resonance.confDaughterPIDspecies.value)) {
            continue; /// loosest cuts for track2
          }

          /// This only works for the case where the mass of opposite charged particles are the same (for example K+/K- have same mass)
          float massPart1 = o2::track::PID::getMass(Resonance.confDaughterPIDspecies.value[0]);
          float massPart2 = o2::track::PID::getMass(Resonance.confDaughterPIDspecies.value[1]);

          /// Resonance
          ROOT::Math::PtEtaPhiMVector tempD1(track1.pt(), track1.eta(), track1.phi(), massPart1);
          ROOT::Math::PtEtaPhiMVector tempD2(track2.pt(), track2.eta(), track2.phi(), massPart2);
          ROOT::Math::PtEtaPhiMVector tempReso = tempD1 + tempD2;
          /// Anti-resonance
          ROOT::Math::PtEtaPhiMVector tempDA1(track1.pt(), track1.eta(), track1.phi(), massPart2);
          ROOT::Math::PtEtaPhiMVector tempDA2(track2.pt(), track2.eta(), track2.phi(), massPart1);
          ROOT::Math::PtEtaPhiMVector tempAntiReso = tempDA1 + tempDA2;

          /// For InvMassQA
          ROOT::Math::PtEtaPhiMVector tempDaughter1MassQA1(track1.pt(), track1.eta(), track1.phi(), massPart1);
          ROOT::Math::PtEtaPhiMVector tempDaughter2MassQA1(track2.pt(), track2.eta(), track2.phi(), massPart1);
          ROOT::Math::PtEtaPhiMVector tempMassQA1 = tempDaughter1MassQA1 + tempDaughter2MassQA1;

          float massQAPart2 = massPart2;
          if (Resonance.confDaughterPIDspecies.value[0] == Resonance.confDaughterPIDspecies.value[1]) {
            massQAPart2 = o2::track::PID::getMass(Resonance.confMassQAPart2PID.value);
          }

          ROOT::Math::PtEtaPhiMVector tempDaughter1MassQA2(track1.pt(), track1.eta(), track1.phi(), massQAPart2);
          ROOT::Math::PtEtaPhiMVector tempDaughter2MassQA2(track2.pt(), track2.eta(), track2.phi(), massQAPart2);
          ROOT::Math::PtEtaPhiMVector tempMassQA2 = tempDaughter1MassQA2 + tempDaughter2MassQA2;

          float massDiff = std::abs(Resonance.confResoInvMass.value - tempReso.M());
          float massDiffAnti = std::abs(Resonance.confResoInvMass.value - tempAntiReso.M());

          bool resoIsNotAnti = true; /// bool for differentianting between particle/antiparticle
          float resoSign = 1.;
          if ((Resonance.confDaughterPIDspecies->size() > 1) && (Resonance.confDaughterPIDspecies.value[0] != Resonance.confDaughterPIDspecies.value[1])) {
            auto [isNormal, WrongCombination] = resoCuts.checkCombination(track1, track2, Resonance.confDaughterPIDspecies.value, massDiff, massDiffAnti, Resonance.confUseMassDiff.value);
            if (WrongCombination) {
              continue;
            }
            if (!isNormal) {
              resoSign = -1.;
            }
            resoIsNotAnti = isNormal;
          }
          /// Resos, where both daughters have the same PID are defaulted to sign 1. and resoIsNotAnti = true

          if (resoIsNotAnti) {
            resoRegistry.fill(HIST("AnalysisQA/Reso/InvMass"), tempReso.M());
            // Fill DaughterQA histos
            resoRegistry.fill(HIST("AnalysisQA/Reso/ResoQA/PosDaughter/Pt"), track1.pt());
            resoRegistry.fill(HIST("AnalysisQA/Reso/ResoQA/NegDaughter/Pt"), track2.pt());
            resoRegistry.fill(HIST("AnalysisQA/Reso/ResoQA/PosDaughter/Eta"), track1.eta());
            resoRegistry.fill(HIST("AnalysisQA/Reso/ResoQA/NegDaughter/Eta"), track2.eta());
            resoRegistry.fill(HIST("AnalysisQA/Reso/ResoQA/PosDaughter/DcaXY"), track1.dcaXY());
            resoRegistry.fill(HIST("AnalysisQA/Reso/ResoQA/NegDaughter/DcaXY"), track2.dcaXY());
            resoRegistry.fill(HIST("AnalysisQA/Reso/ResoQA/PosDaughter/DcaZ"), track1.dcaZ());
            resoRegistry.fill(HIST("AnalysisQA/Reso/ResoQA/NegDaughter/DcaZ"), track2.dcaZ());
            resoRegistry.fill(HIST("AnalysisQA/Reso/ResoQA/PosDaughter/Phi"), track1.phi());
            resoRegistry.fill(HIST("AnalysisQA/Reso/ResoQA/NegDaughter/Phi"), track2.phi());
            // Fill InvMassQA histos
            resoRegistry.fill(HIST("AnalysisQA/Reso/ResoQA/InvMassSwitched"), tempAntiReso.M());
            resoRegistry.fill(HIST("AnalysisQA/Reso/ResoQA/InvMassBothPID1"), tempMassQA1.M());
            resoRegistry.fill(HIST("AnalysisQA/Reso/ResoQA/InvMassBothPID2"), tempMassQA2.M());

            /// MassCut
            if (!(tempReso.M() > Resonance.confResoInvMassLowLimit.value && tempReso.M() < Resonance.confResoInvMassUpLimit.value))
              continue;
            resoRegistry.fill(HIST("AnalysisQA/Reso/InvMass_phi_selected"), tempReso.M());
          } else {
            resoRegistry.fill(HIST("AnalysisQA/Reso/InvMassAnti"), tempAntiReso.M());
            // Fill DaughterQA histos
            resoRegistry.fill(HIST("AnalysisQA/Reso/AntiResoQA/PosDaughter/Pt"), track1.pt());
            resoRegistry.fill(HIST("AnalysisQA/Reso/AntiResoQA/NegDaughter/Pt"), track2.pt());
            resoRegistry.fill(HIST("AnalysisQA/Reso/AntiResoQA/PosDaughter/Eta"), track1.eta());
            resoRegistry.fill(HIST("AnalysisQA/Reso/AntiResoQA/NegDaughter/Eta"), track2.eta());
            resoRegistry.fill(HIST("AnalysisQA/Reso/AntiResoQA/PosDaughter/DcaXY"), track1.dcaXY());
            resoRegistry.fill(HIST("AnalysisQA/Reso/AntiResoQA/NegDaughter/DcaXY"), track2.dcaXY());
            resoRegistry.fill(HIST("AnalysisQA/Reso/AntiResoQA/PosDaughter/DcaZ"), track1.dcaZ());
            resoRegistry.fill(HIST("AnalysisQA/Reso/AntiResoQA/NegDaughter/DcaZ"), track2.dcaZ());
            resoRegistry.fill(HIST("AnalysisQA/Reso/AntiResoQA/PosDaughter/Phi"), track1.phi());
            resoRegistry.fill(HIST("AnalysisQA/Reso/AntiResoQA/NegDaughter/Phi"), track2.phi());
            // Fill InvMassQA histos
            resoRegistry.fill(HIST("AnalysisQA/Reso/AntiResoQA/InvMassSwitched"), tempReso.M());
            resoRegistry.fill(HIST("AnalysisQA/Reso/AntiResoQA/InvMassBothPID1"), tempMassQA1.M());
            resoRegistry.fill(HIST("AnalysisQA/Reso/AntiResoQA/InvMassBothPID2"), tempMassQA2.M());

            /// MassCut
            if (!(tempAntiReso.M() > Resonance.confResoInvMassLowLimit.value && tempAntiReso.M() < Resonance.confResoInvMassUpLimit.value))
              continue;
            resoRegistry.fill(HIST("AnalysisQA/Reso/InvMassAnti_phi_selected"), tempAntiReso.M());
          }

          resoCuts.fillQA<aod::femtodreamparticle::ParticleType::kResoChild,
                          aod::femtodreamparticle::TrackType::kPosChild,
                          aod::femtodreamparticle::TrackType::kNegChild>(track1, track2);

          resoRegistry.fill(HIST("AnalysisQA/Reso/Pt_posdaughter_selected"), track1.pt());
          resoRegistry.fill(HIST("AnalysisQA/Reso/Pt_negdaughter_selected"), track2.pt());
          resoRegistry.fill(HIST("AnalysisQA/Reso/Eta_posdaughter_selected"), track1.eta());
          resoRegistry.fill(HIST("AnalysisQA/Reso/Eta_negdaughter_selected"), track2.eta());
          resoRegistry.fill(HIST("AnalysisQA/Reso/DCAxy_posdaughter_selected"), track1.dcaXY());
          resoRegistry.fill(HIST("AnalysisQA/Reso/DCAxy_negdaughter_selected"), track2.dcaXY());
          resoRegistry.fill(HIST("AnalysisQA/Reso/DCAz_posdaughter_selected"), track1.dcaZ());
          resoRegistry.fill(HIST("AnalysisQA/Reso/DCAz_negdaughter_selected"), track2.dcaZ());
          resoRegistry.fill(HIST("AnalysisQA/Reso/Phi_posdaughter_selected"), track1.phi());
          resoRegistry.fill(HIST("AnalysisQA/Reso/Phi_negdaughter_selected"), track2.phi());

          auto type = resoCuts.getType(track1, track2, resoIsNotAnti); //   kResoPosdaughTPC_NegdaughTPC
                                                                       //   kResoPosdaughTPC_NegdaughTOF
                                                                       //   kResoPosdaughTPC_NegdaughTPC
                                                                       //   kResoPosdaughTOF_NegdaughTOF as possible output

          auto bitmask = resoCuts.getCutContainer<aod::femtodreamparticle::cutContainerType>(track1, track2, resoSign);

          /// Get Variables for Output
          auto outputReso = tempReso;
          auto outputDaugh1 = tempD1;
          auto outputDaugh2 = tempD2;
          if (!resoIsNotAnti) {
            outputReso = tempAntiReso;
            outputDaugh1 = tempDA1;
            outputDaugh2 = tempDA2;
          }

          // fill FDParticles
          int postrkId = track1.globalIndex();
          int rowOfPosTrack = -1;
          rowOfPosTrack = getRowDaughters(postrkId, tmpIDtrack);

          childIDs[0] = rowOfPosTrack; // should give me the row
          childIDs[1] = 0;
          outputParts(outputCollision.lastIndex(),
                      track1.pt(),
                      track1.eta(),
                      track1.phi(),
                      aod::femtodreamparticle::ParticleType::kResoChild,
                      bitmask[1],
                      bitmask[2],
                      track1.dcaXY(),
                      childIDs,
                      outputDaugh1.M(),
                      outputDaugh2.M()); // fill tempFitVar with dcaXY?
          const int rowPosTrk = outputParts.lastIndex();

          int negtrkId = track2.globalIndex();
          int rowOfNegTrack = -1;
          rowOfNegTrack = getRowDaughters(negtrkId, tmpIDtrack);

          childIDs[0] = 0;
          childIDs[1] = rowOfNegTrack;
          outputParts(outputCollision.lastIndex(),
                      track2.pt(),
                      track2.eta(),
                      track2.phi(),
                      aod::femtodreamparticle::ParticleType::kResoChild,
                      bitmask[3],
                      bitmask[4],
                      track2.dcaXY(),
                      childIDs,
                      outputDaugh2.M(),
                      outputDaugh1.M()); // maybe CPA instead of dcaXY()? as tempFitVar?
          const int rowNegTrk = outputParts.lastIndex();

          // Reso
          std::vector<int> indexChildIds = {rowPosTrk, rowNegTrk};
          outputParts(outputCollision.lastIndex(),
                      outputReso.pt(),
                      outputReso.eta(),
                      outputReso.phi(),
                      type,
                      bitmask[0], // PIDBit of neg_daugh merged with sign cutBit
                      bitmask[2], // PIDBit of pos_daugh
                      -999.f,
                      indexChildIds,
                      tempReso.M(),
                      tempAntiReso.M()); // no TempFitVar !!
          // needed?
          if (confIsDebug.value) {
            fillDebugParticle<true, false>(track1); // QA for positive daughter
            fillDebugParticle<true, false>(track2); // QA for negative daughter
            fillDebugParticle<false, false, ROOT::Math::PtEtaPhiMVector, true>(outputReso);
          }
        } // for (const auto& track2 : sliceNegdaugh)
      } // for (const auto& track1 : slicePosdaugh)

      if (Resonance.confDoLikeSign.value) {
        fillLikeSign(slicePosdaugh);
        fillLikeSign(sliceNegdaugh);
      }

    } // if (confIsActivatePhi.value)
  } // void fillCollisionsAndTracksAndV0(...)

  void
    processData(aod::FemtoFullCollision const& col,
                aod::BCsWithTimestamps const&,
                aod::FemtoFullTracks const& tracks,
                o2::aod::V0Datas const& fullV0s,
                o2::aod::CascDatas const& fullCascades)
  {
    // get magnetic field for run
    initCcdbMagTrig(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    auto tracksWithItsPid = soa::Attach<aod::FemtoFullTracks, aod::pidits::ITSNSigmaEl, aod::pidits::ITSNSigmaPi, aod::pidits::ITSNSigmaKa,
                                        aod::pidits::ITSNSigmaPr, aod::pidits::ITSNSigmaDe, aod::pidits::ITSNSigmaTr, aod::pidits::ITSNSigmaHe>(tracks);

    if (confUseItsPid.value) {
      fillCollisionsAndTracksAndV0AndCascade<false, true, true, false>(col, tracks, tracksWithItsPid, fullV0s, fullCascades);
    } else {
      fillCollisionsAndTracksAndV0AndCascade<false, false, true, false>(col, tracks, tracks, fullV0s, fullCascades);
    }
  }
  PROCESS_SWITCH(FemtoDreamProducerTaskReso, processData,
                 "Provide experimental data", true);

  void
    processData_noCentrality(aod::FemtoFullCollisionNoCent const& col,
                             aod::BCsWithTimestamps const&,
                             aod::FemtoFullTracks const& tracks,
                             o2::aod::V0Datas const& fullV0s,
                             o2::aod::CascDatas const& fullCascades)
  {
    // get magnetic field for run
    initCcdbMagTrig(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    auto tracksWithItsPid = soa::Attach<aod::FemtoFullTracks, aod::pidits::ITSNSigmaEl, aod::pidits::ITSNSigmaPi, aod::pidits::ITSNSigmaKa,
                                        aod::pidits::ITSNSigmaPr, aod::pidits::ITSNSigmaDe, aod::pidits::ITSNSigmaTr, aod::pidits::ITSNSigmaHe>(tracks);

    if (confUseItsPid.value) {
      fillCollisionsAndTracksAndV0AndCascade<false, true, false, false>(col, tracks, tracksWithItsPid, fullV0s, fullCascades);
    } else {
      fillCollisionsAndTracksAndV0AndCascade<false, false, false, false>(col, tracks, tracks, fullV0s, fullCascades);
    }
  }
  PROCESS_SWITCH(FemtoDreamProducerTaskReso, processData_noCentrality,
                 "Provide experimental data without centrality information", false);

  void processDataCentPbPb(aod::FemtoFullCollisionCentPbPb const& col,
                           aod::BCsWithTimestamps const&,
                           aod::FemtoFullTracks const& tracks,
                           o2::aod::V0Datas const& fullV0s,
                           o2::aod::CascDatas const& fullCascades)
  {
    // get magnetic field for run
    initCcdbMagTrig(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    auto tracksWithItsPid = soa::Attach<aod::FemtoFullTracks, aod::pidits::ITSNSigmaEl, aod::pidits::ITSNSigmaPi, aod::pidits::ITSNSigmaKa,
                                        aod::pidits::ITSNSigmaPr, aod::pidits::ITSNSigmaDe, aod::pidits::ITSNSigmaTr, aod::pidits::ITSNSigmaHe>(tracks);

    if (confUseItsPid.value) {
      fillCollisionsAndTracksAndV0AndCascade<false, true, true, true>(col, tracks, tracksWithItsPid, fullV0s, fullCascades);
    } else {
      fillCollisionsAndTracksAndV0AndCascade<false, false, true, true>(col, tracks, tracks, fullV0s, fullCascades);
    }
  }
  PROCESS_SWITCH(FemtoDreamProducerTaskReso, processDataCentPbPb,
                 "Provide experimental data with centrality information for PbPb collisions", false);

  void processMC(aod::FemtoFullCollisionMC const& col,
                 aod::BCsWithTimestamps const&,
                 soa::Join<aod::FemtoFullTracks, aod::McTrackLabels> const& tracks,
                 aod::FemtoFullMCgenCollisions const&,
                 aod::McParticles const&,
                 soa::Join<o2::aod::V0Datas, aod::McV0Labels> const& fullV0s, /// \todo with FilteredFullV0s
                 soa::Join<o2::aod::CascDatas, aod::McCascLabels> const& fullCascades)
  {
    // get magnetic field for run
    initCcdbMagTrig(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    fillCollisionsAndTracksAndV0AndCascade<false, false, true, false>(col, tracks, tracks, fullV0s, fullCascades);
  }
  PROCESS_SWITCH(FemtoDreamProducerTaskReso, processMC, "Provide MC data", false);

  void processMCnoCentrality(aod::FemtoFullCollisionNoCentMC const& col,
                             aod::BCsWithTimestamps const&,
                             soa::Join<aod::FemtoFullTracks, aod::McTrackLabels> const& tracks,
                             aod::FemtoFullMCgenCollisions const&,
                             aod::McParticles const&,
                             soa::Join<o2::aod::V0Datas, aod::McV0Labels> const& fullV0s, /// \todo with FilteredFullV0s
                             soa::Join<o2::aod::CascDatas, aod::McCascLabels> const& fullCascades)
  {
    // get magnetic field for run
    initCcdbMagTrig(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    fillCollisionsAndTracksAndV0AndCascade<true, false, false, false>(col, tracks, tracks, fullV0s, fullCascades);
  }
  PROCESS_SWITCH(FemtoDreamProducerTaskReso, processMCnoCentrality, "Provide MC data without requiring a centrality calibration", false);

  void processMCCentPbPb(aod::FemtoFullCollisionMCCentPbPb const& col,
                         aod::BCsWithTimestamps const&,
                         soa::Join<aod::FemtoFullTracks, aod::McTrackLabels> const& tracks,
                         aod::FemtoFullMCgenCollisions const&,
                         aod::McParticles const&,
                         soa::Join<o2::aod::V0Datas, aod::McV0Labels> const& fullV0s, /// \todo with FilteredFullV0s
                         soa::Join<o2::aod::CascDatas, aod::McCascLabels> const& fullCascades)
  {
    // get magnetic field for run
    initCcdbMagTrig(col.bc_as<aod::BCsWithTimestamps>());
    // fill the tables
    fillCollisionsAndTracksAndV0AndCascade<true, false, true, true>(col, tracks, tracks, fullV0s, fullCascades);
  }
  PROCESS_SWITCH(FemtoDreamProducerTaskReso, processMCCentPbPb, "Provide MC data with centrality information for PbPb collisions", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<FemtoDreamProducerTaskReso>(cfgc)};
  return workflow;
}
